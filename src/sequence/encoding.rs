//! Base encoding and decoding.
//!
//! Reads are stored with four bases per byte exactly like the original
//! implementation: A/C/G/T map to 0/1/2/3 and every other character (notably
//! `N`) is stored as `A` while its position is recorded in a side table. The
//! side table is what allows the exact original sequence to be rebuilt, so the
//! encoding has to match byte for byte.

/// `seq_nt6_table`: ASCII to 2-bit code, 5 for anything that is not a base.
pub const SEQUENCE_NT6_TABLE: [u8; 256] = build_nt6_table();

/// `seq_nt4_table` from `htab.cpp`: like [`SEQUENCE_NT6_TABLE`] but with 4 as
/// the "anything else" code, which the bit-vector alignment code relies on
/// because it indexes a five-entry table with the result.
pub const SEQUENCE_NT4_TABLE: [u8; 256] = build_nt4_table();

/// `s_H`: the four bases plus `N`.
pub const BASE_CHARACTERS: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];

/// `rc_Table`: complement of every `seq_nt6_table` value.
pub const COMPLEMENT_CHARACTERS: [u8; 6] = [b'T', b'G', b'C', b'A', b'N', b'N'];

const fn build_nt4_table() -> [u8; 256] {
    let mut table = [4u8; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    // The reference table also maps 'u' to 3.
    table[b'u' as usize] = 3;
    table
}

const fn build_nt6_table() -> [u8; 256] {
    let mut table = [5u8; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
}

/// Two lookup tables used to expand a packed byte back into four bases.
pub struct DecodingTables {
    /// `bit_t_seq_table`: forward expansion of one packed byte.
    pub forward: [[u8; 4]; 256],
    /// `bit_t_seq_table_rc`: reverse complement expansion of one packed byte.
    pub reverse_complement: [[u8; 4]; 256],
}

impl DecodingTables {
    fn new() -> Self {
        let mut forward = [[0u8; 4]; 256];
        let mut reverse_complement = [[0u8; 4]; 256];

        for byte_value in 0..256usize {
            for position in 0..4 {
                let shift = 6 - 2 * position;
                forward[byte_value][position] = BASE_CHARACTERS[(byte_value >> shift) & 3];
            }
            reverse_complement[byte_value][0] = complement_character(forward[byte_value][3]);
            reverse_complement[byte_value][1] = complement_character(forward[byte_value][2]);
            reverse_complement[byte_value][2] = complement_character(forward[byte_value][1]);
            reverse_complement[byte_value][3] = complement_character(forward[byte_value][0]);
        }

        Self {
            forward,
            reverse_complement,
        }
    }
}

/// Complement of an ASCII base, the `RC_CHAR` macro.
pub fn complement_character(base: u8) -> u8 {
    COMPLEMENT_CHARACTERS[SEQUENCE_NT6_TABLE[base as usize] as usize]
}

/// Process wide decoding tables, initialised once like `init_aux_table`.
pub fn decoding_tables() -> &'static DecodingTables {
    use std::sync::OnceLock;
    static TABLES: OnceLock<DecodingTables> = OnceLock::new();
    TABLES.get_or_init(DecodingTables::new)
}

/// A read stored in the packed form used by `All_reads`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PackedSequence {
    /// Four bases per byte, most significant pair first.
    pub packed_bases: Vec<u8>,
    /// Number of bases, which is not `packed_bases.len() * 4`.
    pub length: u64,
    /// Positions that did not encode a plain A/C/G/T, in ascending order.
    pub ambiguous_positions: Vec<u64>,
}

impl PackedSequence {
    /// `ha_compress_base`: packs `source` and records every non A/C/G/T
    /// position.
    pub fn from_ascii(source: &[u8]) -> Self {
        let mut packed_bases = Vec::with_capacity(source.len() / 4 + 1);
        let mut ambiguous_positions = Vec::new();

        let mut position = 0usize;
        while position + 4 <= source.len() {
            let mut packed_byte = 0u8;
            for step in 0..4 {
                let code = encode_base(
                    source[position + step],
                    position + step,
                    &mut ambiguous_positions,
                );
                packed_byte |= code << (6 - 2 * step);
            }
            packed_bases.push(packed_byte);
            position += 4;
        }

        if position < source.len() {
            let mut packed_byte = 0u8;
            let mut shift = 6;
            while position < source.len() {
                let code = encode_base(source[position], position, &mut ambiguous_positions);
                packed_byte |= code << shift;
                shift -= 2;
                position += 1;
            }
            packed_bases.push(packed_byte);
        }

        Self {
            packed_bases,
            length: source.len() as u64,
            ambiguous_positions,
        }
    }

    /// Restores the original sequence.
    pub fn unpack(&self) -> Vec<u8> {
        let tables = decoding_tables();
        let mut sequence = Vec::with_capacity(self.length as usize + 4);
        let mut remaining = self.length as usize;

        for packed_byte in &self.packed_bases {
            let chunk = &tables.forward[*packed_byte as usize];
            let take = remaining.min(4);
            sequence.extend_from_slice(&chunk[..take]);
            remaining -= take;
            if remaining == 0 {
                break;
            }
        }

        for position in &self.ambiguous_positions {
            sequence[*position as usize] = b'N';
        }

        sequence
    }

    /// Restores the reverse complement of the original sequence.
    pub fn unpack_reverse_complement(&self) -> Vec<u8> {
        let mut sequence = self.unpack();
        reverse_complement(&mut sequence);
        sequence
    }

    /// Extracts `[start_position, start_position + length)` on either strand,
    /// the equivalent of `recover_UC_Read_sub_region`.
    ///
    /// The reference reads the packed buffer without checking the end, so a
    /// window that runs past the read yields whatever the zero-initialised
    /// buffer holds, which decodes to `A`. Reading past the end here would
    /// panic, so the same bytes are produced explicitly.
    ///
    /// Only the requested interval is decoded: the whole-read unpack of this
    /// function is on the hot path of every alignment attempt, and decoding
    /// the entire read for a few hundred bases dominates the correction time.
    pub fn sub_region(&self, start_position: i64, length: i64, reverse_strand: bool) -> Vec<u8> {
        let tables = decoding_tables();
        let requested = length.max(0) as usize;
        let mut region = vec![b'A'; requested];
        let read_length = self.length as i64;
        if requested == 0 || read_length == 0 {
            return region;
        }

        if reverse_strand {
            // Output position `p` is the complement of read position
            // `read_length - 1 - (start + p)`; positions outside the read stay
            // `A`, exactly as the whole-read unpack followed by a slice does.
            let original_start = read_length - 1 - start_position.max(0);
            if original_start >= 0 {
                let in_range = requested.min(original_start as usize + 1);
                let mut output = 0usize;
                while output < in_range {
                    let original = original_start - output as i64;
                    let byte_index = (original >> 2) as usize;
                    let offset_in_byte = (original & 3) as usize;
                    let packed = self.packed_bases.get(byte_index).copied().unwrap_or(0);
                    let chunk = &tables.reverse_complement[packed as usize];
                    // `reverse_complement[k]` is the base at byte position
                    // `3 - k`, and positions run downwards here.
                    let chunk_start = 3 - offset_in_byte;
                    let take = (in_range - output).min(4 - chunk_start);
                    region[output..output + take]
                        .copy_from_slice(&chunk[chunk_start..chunk_start + take]);
                    output += take;
                }

                // `N` sites are stored in ascending read order; their output
                // positions therefore run downwards.
                let last_site = original_start;
                let first_site = original_start - (requested as i64 - 1);
                let first = self
                    .ambiguous_positions
                    .partition_point(|&site| (site as i64) < first_site);
                for &site in &self.ambiguous_positions[first..] {
                    let site = site as i64;
                    if site > last_site {
                        break;
                    }
                    region[(original_start - site) as usize] = b'N';
                }
            }
        } else {
            let start = start_position.max(0);
            let end = (start + requested as i64).min(read_length);
            let mut position = start;
            while position < end {
                let byte_index = (position >> 2) as usize;
                let offset_in_byte = (position & 3) as usize;
                let packed = self.packed_bases.get(byte_index).copied().unwrap_or(0);
                let chunk = &tables.forward[packed as usize];
                let take = ((end - position) as usize).min(4 - offset_in_byte);
                let output = (position - start) as usize;
                region[output..output + take]
                    .copy_from_slice(&chunk[offset_in_byte..offset_in_byte + take]);
                position += take as i64;
            }

            let first = self
                .ambiguous_positions
                .partition_point(|&site| (site as i64) < start);
            for &site in &self.ambiguous_positions[first..] {
                let site = site as i64;
                if site >= end {
                    break;
                }
                region[(site - start) as usize] = b'N';
            }
        }

        region
    }
}

fn encode_base(base: u8, position: usize, ambiguous_positions: &mut Vec<u64>) -> u8 {
    let code = SEQUENCE_NT6_TABLE[base as usize];
    if code >= 4 {
        ambiguous_positions.push(position as u64);
        0
    } else {
        code
    }
}

/// `reverse_complement`: in place reverse complement of an ASCII sequence.
pub fn reverse_complement(sequence: &mut [u8]) {
    let length = sequence.len();
    let half = length / 2;
    for index in 0..half {
        let opposite = length - index - 1;
        let temporary = sequence[opposite];
        sequence[opposite] = complement_character(sequence[index]);
        sequence[index] = complement_character(temporary);
    }
    if length % 2 == 1 {
        sequence[half] = complement_character(sequence[half]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nt6_table_matches_the_original_layout() {
        assert_eq!(SEQUENCE_NT6_TABLE[b'A' as usize], 0);
        assert_eq!(SEQUENCE_NT6_TABLE[b'C' as usize], 1);
        assert_eq!(SEQUENCE_NT6_TABLE[b'G' as usize], 2);
        assert_eq!(SEQUENCE_NT6_TABLE[b'T' as usize], 3);
        assert_eq!(SEQUENCE_NT6_TABLE[b'N' as usize], 5);
        assert_eq!(SEQUENCE_NT6_TABLE[b'a' as usize], 0);
    }

    #[test]
    fn packing_follows_the_most_significant_pair_first_layout() {
        let packed = PackedSequence::from_ascii(b"ACGT");
        assert_eq!(packed.packed_bases, vec![0b00_01_10_11]);
        assert_eq!(packed.length, 4);
        assert!(packed.ambiguous_positions.is_empty());
    }

    #[test]
    fn round_trip_restores_ambiguous_bases() {
        for sequence in [
            b"ACGT".to_vec(),
            b"ACGTN".to_vec(),
            b"NACGTACGTN".to_vec(),
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec(),
            b"GN".to_vec(),
            b"N".to_vec(),
        ] {
            let packed = PackedSequence::from_ascii(&sequence);
            assert_eq!(packed.unpack(), sequence);
        }
    }

    #[test]
    fn reverse_complement_matches_the_reference_implementation() {
        let mut sequence = b"ACGTN".to_vec();
        reverse_complement(&mut sequence);
        assert_eq!(sequence, b"NACGT");

        let mut even_sequence = b"AACC".to_vec();
        reverse_complement(&mut even_sequence);
        assert_eq!(even_sequence, b"GGTT");
    }

    #[test]
    fn sub_region_on_both_strands_matches_unpacking() {
        let sequence = b"ACGTTGCANNGTCA";
        let packed = PackedSequence::from_ascii(sequence);

        assert_eq!(packed.sub_region(2, 6, false), b"GTTGCA");

        let reverse_complemented = packed.unpack_reverse_complement();
        assert_eq!(
            packed.sub_region(0, sequence.len() as i64, true),
            reverse_complemented
        );
    }

    /// The direct decoder must agree with the whole-read unpack plus a slice
    /// for every window, including the ones that run past the read end (which
    /// the reference fills with `A` from its zero-initialised buffer).
    #[test]
    fn sub_region_matches_the_whole_read_unpack_at_every_offset() {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut sequence = Vec::new();
        for _ in 0..37 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            // Deliberately include `N` so the ambiguous-position path is used.
            sequence.push(b"ACGTN"[((state >> 3) & 3) as usize]);
        }
        let packed = PackedSequence::from_ascii(&sequence);

        let forward = packed.unpack();
        let reverse = packed.unpack_reverse_complement();
        let read_length = sequence.len() as i64;
        for start in -3..=read_length + 2 {
            for length in 0..=read_length + 4 {
                for strand in [false, true] {
                    let reference = if strand { &reverse } else { &forward };
                    let begin = start.max(0).min(read_length);
                    let end = (begin + length).min(read_length);
                    let mut expected = reference[begin as usize..end as usize].to_vec();
                    expected.resize(length.max(0) as usize, b'A');
                    assert_eq!(
                        packed.sub_region(start, length, strand),
                        expected,
                        "start {start}, length {length}, strand {strand}"
                    );
                }
            }
        }
    }
}

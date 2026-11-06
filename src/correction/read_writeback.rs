//! Applying a correction cigar to a read and writing the corrected reads,
//! the counterparts of `worker_sl_ec` and `Output_corrected_reads`.

use std::io::{self, Write};

use crate::correction::cigar::{
    CigarRecord, OPERATION_DELETION, OPERATION_INSERTION, OPERATION_MATCH, OPERATION_MISMATCH,
};
use crate::reads::read_store::ReadStore;
use crate::sequence::encoding::BASE_CHARACTERS;

/// Rebuilds a corrected read from the original read and its correction cigar.
///
/// The cigar cells are visited in order; `lost_bases` carries the replacement
/// bases for mismatches (packed as `new << 3 | old`), insertions and deletions.
/// This is the same reconstruction `worker_sl_ec` performs on the compressed
/// cigar stored per read.
pub fn corrected_read_from_cigar(original: &[u8], cigar: &CigarRecord) -> Vec<u8> {
    let mut corrected = Vec::with_capacity(cigar.new_read_length as usize);
    let mut original_index = 0usize;
    let mut lost_index = 0usize;

    for (operation, length) in cigar.cells() {
        match operation {
            OPERATION_MATCH => {
                for _ in 0..length {
                    corrected.push(
                        original
                            .get(original_index)
                            .copied()
                            .unwrap_or(b'N'),
                    );
                    original_index += 1;
                }
            }
            OPERATION_MISMATCH => {
                for _ in 0..length {
                    let merged = cigar.lost_bases.get(lost_index).copied().unwrap_or(0);
                    lost_index += 1;
                    // Packed as `new_base << 3 | old_base`; `s_H` maps the
                    // 2-bit code back to the character.
                    let new_base = BASE_CHARACTERS[(merged >> 3) as usize & 7];
                    corrected.push(new_base);
                    original_index += 1;
                }
            }
            OPERATION_INSERTION => {
                for _ in 0..length {
                    corrected.push(cigar.lost_bases.get(lost_index).copied().unwrap_or(b'N'));
                    lost_index += 1;
                }
            }
            OPERATION_DELETION => {
                for _ in 0..length {
                    // The reference records the deleted bases but does not
                    // emit them.
                    lost_index += 1;
                }
                original_index += length as usize;
            }
            other => panic!("unknown cigar operation {other}"),
        }
    }

    corrected
}

/// `sl_ec_r` for one read: replaces the stored sequence with its corrected
/// version and updates the length and total base counters.
pub fn write_corrected_read_back(
    reads: &mut ReadStore,
    read_identifier: usize,
    corrected: &[u8],
) {
    let previous_length = reads.read_length(read_identifier);
    reads.total_bases = reads.total_bases - previous_length + corrected.len() as u64;
    reads.read_lengths[read_identifier] = corrected.len() as u64;
    reads.sequences[read_identifier] =
        crate::sequence::encoding::PackedSequence::from_ascii(corrected);
}

/// `worker_sl_ec`'s reconstruction: rebuilds a read from its **consensus
/// cigar** (`scc`), the cells produced by correction rather than the final
/// overlap cigar.
///
/// Match cells copy the original base, mismatch and insertion cells emit the
/// target base stored in bits 12–13 of a cell, and deletion cells only consume
/// the original. The reference's first pass adds the output length over every
/// cell except deletions and the error count over every cell except matches;
/// when the read has no error (`tot_e == 0`) it is left untouched, which this
/// function reports as `None`.
pub fn corrected_read_from_consensus_cigar(original: &[u8], cells: &[u16]) -> Option<Vec<u8>> {
    use crate::overlap::stored_overlap::pop_trace_bp_f_full;

    let mut cell = 0usize;
    let mut output_length = 0usize;
    let mut error_length = 0usize;
    while cell < cells.len() {
        let (next, operation, _, _, length) = pop_trace_bp_f_full(cells, cell);
        if operation != 3 {
            output_length += length as usize;
        }
        if operation != 0 {
            error_length += length as usize;
        }
        cell = next;
    }
    if error_length == 0 {
        return None;
    }

    let mut corrected = vec![0u8; output_length];
    let mut original_index = 0usize;
    let mut corrected_index = 0usize;
    cell = 0;
    while cell < cells.len() {
        let (next, operation, _, target_base, length) = pop_trace_bp_f_full(cells, cell);
        let length = length as usize;
        let original_window_start = original_index;
        let corrected_window_start = corrected_index;
        if operation != 2 {
            original_index += length;
        }
        if operation != 3 {
            corrected_index += length;
        }
        match operation {
            0 => {
                corrected[corrected_window_start..corrected_index].copy_from_slice(
                    &original[original_window_start..original_index],
                );
            }
            1 | 2 => {
                let base = BASE_CHARACTERS[target_base as usize & 3];
                for position in corrected_window_start..corrected_index {
                    corrected[position] = base;
                }
            }
            _ => {}
        }
        cell = next;
    }

    Some(corrected)
}

/// `worker_sl_ec` for one read: applies the read's consensus cigar to the
/// stored sequence. Returns whether the read changed (`false` when the cigar
/// has no error and the reference leaves it alone).
pub fn apply_consensus_to_read(reads: &mut ReadStore, read_identifier: usize, cells: &[u16]) -> bool {
    let original = reads.read_sequence(read_identifier);
    match corrected_read_from_consensus_cigar(&original, cells) {
        Some(corrected) => {
            write_corrected_read_back(reads, read_identifier, &corrected);
            true
        }
        None => false,
    }
}

/// `Output_corrected_reads`: writes `>name\\nsequence\\n` for every read.
pub fn write_corrected_reads<W: Write>(writer: &mut W, reads: &ReadStore) -> io::Result<()> {
    for identifier in 0..reads.total_reads as usize {
        writer.write_all(b">")?;
        writer.write_all(reads.read_name(identifier).as_bytes())?;
        writer.write_all(b"\n")?;
        writer.write_all(&reads.read_sequence(identifier))?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

/// `Output_corrected_fastq`: writes a FASTQ record with placeholder quality
/// (`!`, the lowest quality) for every read.
pub fn write_corrected_fastq<W: Write>(writer: &mut W, reads: &ReadStore) -> io::Result<()> {
    for identifier in 0..reads.total_reads as usize {
        let sequence = reads.read_sequence(identifier);
        writer.write_all(b"@")?;
        writer.write_all(reads.read_name(identifier).as_bytes())?;
        writer.write_all(b"\n")?;
        writer.write_all(&sequence)?;
        writer.write_all(b"\n+\n")?;
        writer.write_all(&vec![b'!'; sequence.len()])?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn store() -> ReadStore {
        let mut reads = ReadStore::default();
        reads.push_read("read_a", b"ACGTACGT");
        reads.push_read("read_b", b"TTTTGGGG");
        reads
    }

    #[test]
    fn a_match_cigar_reproduces_the_original() {
        let mut cigar = CigarRecord::new();
        cigar.add(&[0u8; 8], 8, OPERATION_MATCH);
        assert_eq!(corrected_read_from_cigar(b"ACGTACGT", &cigar), b"ACGTACGT");
    }

    #[test]
    fn a_substitution_uses_the_new_base() {
        let mut cigar = CigarRecord::new();
        cigar.add(&[0u8; 3], 3, OPERATION_MATCH);
        // new base C (code 1), old base G (code 2)
        cigar.add(&[(1 << 3) | 2], 1, OPERATION_MISMATCH);
        cigar.add(&[0u8; 4], 4, OPERATION_MATCH);
        assert_eq!(corrected_read_from_cigar(b"ACGTACGT", &cigar), b"ACGCACGT");
    }

    #[test]
    fn insertions_and_deletions_change_the_length() {
        let mut cigar = CigarRecord::new();
        cigar.add(&[0u8; 2], 2, OPERATION_MATCH);
        cigar.add(b"TT", 2, OPERATION_INSERTION);
        cigar.add(&[0u8; 2], 2, OPERATION_MATCH);
        cigar.add(b"AA", 2, OPERATION_DELETION);
        cigar.add(&[0u8; 2], 2, OPERATION_MATCH);
        assert_eq!(
            corrected_read_from_cigar(b"ACGTACGT", &cigar),
            b"ACTTGTGT"
        );
    }

    #[test]
    fn writeback_updates_length_and_total_bases() {
        let mut reads = store();
        assert_eq!(reads.total_bases, 16);
        write_corrected_read_back(&mut reads, 0, b"ACG");
        assert_eq!(reads.read_length(0), 3);
        assert_eq!(reads.read_sequence(0), b"ACG");
        assert_eq!(reads.total_bases, 11);
    }

    #[test]
    fn corrected_fasta_matches_the_reference_layout() {
        let reads = store();
        let mut output = Vec::new();
        write_corrected_reads(&mut output, &reads).unwrap();
        assert_eq!(
            String::from_utf8(output).unwrap(),
            ">read_a\nACGTACGT\n>read_b\nTTTTGGGG\n"
        );
    }

    #[test]
    fn corrected_fastq_has_matching_quality_length() {
        let reads = store();
        let mut output = Vec::new();
        write_corrected_fastq(&mut output, &reads).unwrap();
        assert_eq!(
            String::from_utf8(output).unwrap(),
            "@read_a\nACGTACGT\n+\n!!!!!!!!\n@read_b\nTTTTGGGG\n+\n!!!!!!!!\n"
        );
    }

    /// Golden values printed by the real `get_corrected_read_from_cigar`.
    #[test]
    fn corrected_read_matches_the_reference_cigar_expansion() {
        struct Case {
            name: &'static str,
            cells: Vec<(u32, u32, u8)>,
            expected: &'static [u8],
        }
        // The third tuple member is the byte stored in `lost_bases`.
        let cases = [
            Case {
                name: "match",
                cells: vec![(0, 8, 0)],
                expected: b"ACGTACGT",
            },
            Case {
                name: "substitution",
                cells: vec![(0, 3, 0), (1, 1, 10), (0, 4, 0)],
                expected: b"ACGCACGT",
            },
            Case {
                name: "indel",
                cells: vec![
                    (0, 2, 0),
                    (2, 2, b'T'),
                    (0, 2, 0),
                    (3, 2, b'A'),
                    (0, 2, 0),
                ],
                expected: b"ACTTGTGT",
            },
        ];

        for case in cases {
            let mut cigar = CigarRecord::new();
            for (operation, length, lost) in &case.cells {
                let bases = vec![*lost; *length as usize];
                cigar.add(&bases, *length, *operation);
            }
            assert_eq!(
                corrected_read_from_cigar(b"ACGTACGT", &cigar),
                case.expected,
                "{}",
                case.name
            );
        }
    }
}

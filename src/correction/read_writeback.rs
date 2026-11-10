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
/// version and updates the length counter.
///
/// `total_reads_bases` deliberately stays at its load-time value, matching
/// the reference: `ha_insert_read_len` counts the original bases and no
/// correction path updates `All_reads::total_reads_bases`.
pub fn write_corrected_read_back(
    reads: &mut ReadStore,
    read_identifier: usize,
    corrected: &[u8],
) {
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

/// `get_wqual` over an already-unpacked quality track: the minimum mean of a
/// `window`-wide sliding sum around `position`.
///
/// `length` is the reference's `Get_READ_LENGTH`, which `worker_sl_ec` has
/// already updated to the corrected length before it rebuilds the quality
/// track. The source levels still describe the pre-correction read, so the
/// window can stop before the end of `levels`.
fn window_quality_levels(levels: &[u8], position: usize, window: usize, length: usize) -> u8 {
    let start = if position + 1 >= window {
        position + 1 - window
    } else {
        0
    };
    let end = (position + window).min(length);
    if start >= end {
        return 0;
    }
    let end = end.min(levels.len());
    if start >= end {
        return 0;
    }
    let track = &levels[start..end];
    let count = window.min(track.len());
    let mut sum: u64 = track[..count].iter().map(|&value| u64::from(value)).sum();
    let mut minimum = sum;
    for index in count..track.len() {
        sum -= u64::from(track[index - window]);
        sum += u64::from(track[index]);
        if sum < minimum {
            minimum = sum;
        }
    }
    (minimum as f64 / count as f64).round() as u8
}

/// The quality half of `worker_sl_ec`: rebuilds the corrected quality track
/// from the consensus cigar.
///
/// Match and mismatch cells copy the original quality levels; insertion
/// cells take the minimum window quality around the insertion point. The
/// rebuilt levels are packed with the reference's 2-bit layout.
fn corrected_quality_from_consensus_cigar(
    reads: &ReadStore,
    read_identifier: usize,
    cells: &[u16],
    corrected_length: usize,
) -> Option<Vec<u8>> {
    use crate::overlap::stored_overlap::pop_trace_bp_f_full;

    let original_length = reads.read_lengths[read_identifier] as i64;
    let source = reads.read_quality(read_identifier, 0, original_length, false);
    if source.is_empty() {
        return None;
    }

    let mut levels = vec![0u8; corrected_length];
    let mut cell = 0usize;
    let mut original_index = 0usize;
    let mut corrected_index = 0usize;
    while cell < cells.len() {
        let (next, operation, _, _, length) = pop_trace_bp_f_full(cells, cell);
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
            0 | 1 => {
                let available = source.len().saturating_sub(original_window_start);
                let copy_length = length.min(available).min(corrected_index - corrected_window_start);
                levels[corrected_window_start..corrected_window_start + copy_length]
                    .copy_from_slice(&source[original_window_start..original_window_start + copy_length]);
            }
            2 => {
                let quality =
                    window_quality_levels(&source, original_window_start, 5, corrected_length);
                for position in corrected_window_start..corrected_index {
                    levels[position] = quality;
                }
            }
            _ => {}
        }
        cell = next;
    }

    Some(crate::reads::read_store::pack_quality_levels(&levels))
}

/// `worker_sl_ec` for one read: applies the read's consensus cigar to the
/// stored sequence. Returns whether the read changed (`false` when the cigar
/// has no error and the reference leaves it alone).
pub fn apply_consensus_to_read(reads: &mut ReadStore, read_identifier: usize, cells: &[u16]) -> bool {
    let original = reads.read_sequence(read_identifier);
    match corrected_read_from_consensus_cigar(&original, cells) {
        Some(corrected) => {
            let quality =
                corrected_quality_from_consensus_cigar(reads, read_identifier, cells, corrected.len());
            write_corrected_read_back(reads, read_identifier, &corrected);
            if let Some(quality) = quality {
                reads.qualities[read_identifier] = quality;
            }
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
        // `sc_tb[level] + 33 - 1`: the 2-bit level expanded back to a Phred
        // character. Reads without qualities fall back to the lowest level.
        let quality = reads.read_quality(
            identifier,
            0,
            sequence.len() as i64,
            false,
        );
        let encoded: Vec<u8> = quality
            .iter()
            .map(|&level| {
                crate::reads::read_store::QUALITY_THRESHOLDS[level as usize]
                    .wrapping_add(33)
                    .wrapping_sub(1)
            })
            .collect();
        writer.write_all(&encoded)?;
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
    fn writeback_updates_length_and_keeps_the_total_frozen() {
        let mut reads = store();
        assert_eq!(reads.total_bases, 16);
        write_corrected_read_back(&mut reads, 0, b"ACG");
        assert_eq!(reads.read_length(0), 3);
        assert_eq!(reads.read_sequence(0), b"ACG");
        assert_eq!(reads.total_bases, 16);
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
        // No quality track: every level decodes to `sc_tb[0] + 32` (`*`).
        assert_eq!(
            String::from_utf8(output).unwrap(),
            "@read_a\nACGTACGT\n+\n********\n@read_b\nTTTTGGGG\n+\n********\n"
        );
    }

    #[test]
    fn corrected_fastq_expands_the_two_bit_levels() {
        let mut reads = ReadStore::default();
        reads.push_read_with_quality("read_a", b"ACGT", Some(b"5I5I"));
        let mut output = Vec::new();
        write_corrected_fastq(&mut output, &reads).unwrap();
        // `5` (Phred 20) is level 1 (`4`), `I` (Phred 40) is level 3 (`H`).
        assert_eq!(
            String::from_utf8(output).unwrap(),
            "@read_a\nACGT\n+\n4H4H\n"
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

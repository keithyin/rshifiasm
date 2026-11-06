//! The traceback variant of the banded bit-parallel edit distance,
//! `Reserve_Banded_BPM_PATH` in `Levenshtein_distance.h`.
//!
//! It runs the same Myers band as [`crate::alignment::banded_bpm`], but keeps
//! the five bit vectors of every column so that the alignment can be walked
//! back afterwards. The resulting `path` uses the reference numbering:
//! `0` match, `1` mismatch, `2` the pattern has an extra base, `3` the text
//! has an extra base. `generate_cigar` turns that path into packed cells.

use crate::alignment::banded_bpm::NO_ALIGNMENT_ERROR;

/// One column of the bit-parallel matrix.
#[derive(Debug, Clone, Copy, Default)]
struct MatrixColumn {
    d0: u64,
    positive_vertical: u64,
    negative_vertical: u64,
    horizontal_positive: u64,
    horizontal_negative: u64,
}

/// Result of `Reserve_Banded_BPM_PATH`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BandedPath {
    /// End position in the pattern, `-1` when no alignment was found.
    pub end_site: i32,
    /// Start position in the pattern.
    pub start_site: i32,
    /// Edit distance of the alignment.
    pub error: u32,
    /// Traceback operations, most recent first (`path[0]` is the last step).
    pub path: Vec<u8>,
}

/// The four bases whose equality masks survive the initial `memset`.
fn is_shifted_base(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

/// Computes one Myers column and returns all five bit vectors.
fn advance_column(
    equality_masks: &[u64; 256],
    text_base: u8,
    positive_vertical: u64,
    negative_vertical: u64,
) -> MatrixColumn {
    let equality = equality_masks[text_base as usize];
    let horizontal = equality | negative_vertical;
    let d0 = ((positive_vertical.wrapping_add(horizontal & positive_vertical))
        ^ positive_vertical)
        | horizontal;
    let horizontal_negative = positive_vertical & d0;
    let horizontal_positive = negative_vertical | !(positive_vertical | d0);
    let shifted = d0 >> 1;
    let next_negative = shifted & horizontal_positive;
    let next_positive = horizontal_negative | !(shifted | horizontal_positive);
    MatrixColumn {
        d0,
        positive_vertical: next_positive,
        negative_vertical: next_negative,
        horizontal_positive,
        horizontal_negative,
    }
}

/// `try_cigar`: an alignment with only substitutions needs no traceback.
///
/// Returns `(start_site, path)` when the ungapped alignment really has exactly
/// `old_error` mismatches.
fn try_ungapped_cigar(
    pattern: &[u8],
    text: &[u8],
    text_length: usize,
    old_end_site: i32,
    old_error: i32,
) -> Option<(i32, Vec<u8>)> {
    let start_site = old_end_site - text_length as i32 + 1;
    if start_site < 0 {
        return None;
    }

    let mut path = vec![0u8; text_length];
    let mut mismatches = 0i32;
    for (index, base) in text[..text_length].iter().enumerate() {
        let pattern_index = index + start_site as usize;
        if pattern.get(pattern_index).copied().unwrap_or(b'N') != *base {
            // The reference fills `path` backwards.
            path[text_length - index - 1] = 1;
            mismatches += 1;
            if mismatches > old_error {
                return None;
            }
        }
    }

    if mismatches == old_error {
        Some((start_site, path))
    } else {
        None
    }
}

/// `Reserve_Banded_BPM_PATH`.
///
/// `old_error` and `old_end_site` come from the earlier scalar alignment of
/// the same window; when they describe a valid ungapped alignment the
/// traceback is skipped entirely. `None` means no alignment within the
/// threshold exists.
pub fn reserve_banded_bpm_path(
    pattern: &[u8],
    pattern_length: usize,
    text: &[u8],
    text_length: usize,
    error_threshold: usize,
    old_error: i32,
    old_end_site: i32,
) -> Option<BandedPath> {
    if text_length == 0 || pattern_length == 0 {
        return None;
    }

    if old_error != -1 && old_end_site != -1 {
        if old_error == 0 {
            return Some(BandedPath {
                end_site: old_end_site,
                start_site: old_end_site - text_length as i32 + 1,
                error: 0,
                // `generate_cigar` special-cases error zero and never reads
                // the path, but a full match path is the faithful value.
                path: vec![0; text_length],
            });
        }
        if let Some((start_site, path)) =
            try_ungapped_cigar(pattern, text, text_length, old_end_site, old_error)
        {
            return Some(BandedPath {
                end_site: old_end_site,
                start_site,
                error: old_error as u32,
                path,
            });
        }
    }

    let band_length = (error_threshold << 1) + 1;
    let mut equality_masks = [0u64; 256];
    let mut band_bit = 1u64;
    for index in 0..band_length {
        if let Some(base) = pattern.get(index) {
            if is_shifted_base(*base) {
                equality_masks[*base as usize] |= band_bit;
            }
        }
        band_bit <<= 1;
    }

    let band_top_bit = 1u64 << (error_threshold << 1);
    // In the reference each column occupies eight words at `matrix_bit[i<<3]`;
    // one struct here carries the same five values.
    let mut columns = vec![MatrixColumn::default(); text_length + 1];
    let mut positive_vertical = 0u64;
    let mut negative_vertical = 0u64;
    let mut error: i32 = 0;
    let threshold = error_threshold as i32;
    let last_high = (error_threshold << 1) as i32;
    let mut band_position = error_threshold << 1;

    let mut text_index = 0usize;
    while text_index < text_length - 1 {
        let column = advance_column(
            &equality_masks,
            text.get(text_index).copied().unwrap_or(b'N'),
            positive_vertical,
            negative_vertical,
        );
        positive_vertical = column.positive_vertical;
        negative_vertical = column.negative_vertical;
        if column.d0 & 1 == 0 {
            error += 1;
            if error - last_high > threshold {
                return None;
            }
        }

        for slot in [b'A', b'C', b'G', b'T'] {
            equality_masks[slot as usize] >>= 1;
        }
        text_index += 1;
        band_position += 1;
        if let Some(base) = pattern.get(band_position) {
            equality_masks[*base as usize] |= band_top_bit;
        }

        // The matrix column is written after the band has slid, using the
        // vectors computed for the current text base.
        columns[text_index] = column;
    }

    let final_column = advance_column(
        &equality_masks,
        text.get(text_index).copied().unwrap_or(b'N'),
        positive_vertical,
        negative_vertical,
    );
    if final_column.d0 & 1 == 0 {
        error += 1;
        if error - last_high > threshold {
            return None;
        }
    }
    columns[text_index + 1] = final_column;

    let site = text_length as i32 - 1;
    let mut return_site = -1i32;
    let mut return_error = NO_ALIGNMENT_ERROR;
    let available_positions = pattern_length as i32 - text_length as i32;

    if error <= threshold && (error as u32) <= return_error {
        return_error = error as u32;
        return_site = site;
    }

    let mut ungap_error = NO_ALIGNMENT_ERROR;
    let mut position = 0i32;
    while position < available_positions {
        error += ((final_column
            .positive_vertical
            .checked_shr(position as u32)
            .unwrap_or(0))
            & 1) as i32;
        error -= ((final_column
            .negative_vertical
            .checked_shr(position as u32)
            .unwrap_or(0))
            & 1) as i32;
        position += 1;
        if error <= threshold && (error as u32) <= return_error {
            return_error = error as u32;
            return_site = site + position;
        }
        if position == threshold {
            ungap_error = error as u32;
        }
    }

    if return_error == NO_ALIGNMENT_ERROR {
        return None;
    }
    if ungap_error <= threshold as u32 && ungap_error == return_error {
        return_site = site + threshold;
    }

    Some(trace_back(
        &columns,
        text_length,
        error_threshold,
        return_site,
        return_error,
    ))
}

/// Walks the stored matrix columns back from the chosen end site.
fn trace_back(
    columns: &[MatrixColumn],
    text_length: usize,
    error_threshold: usize,
    end_site: i32,
    error: u32,
) -> BandedPath {
    let band_length = (error_threshold << 1) + 1;
    let corrected_pattern_length = text_length as i32 + 2 * error_threshold as i32;
    let mut start_site = end_site;
    let mut back_track_site = band_length as i32 - (corrected_pattern_length - end_site);
    let low_bound = band_length as i32 - 1;

    let mut path: Vec<u8> = Vec::new();
    let mut current_value = error as i64;
    let mut text_index = text_length as i32;
    let mut direction: i64 = 0;

    while text_index > 0 {
        if current_value == 0 {
            break;
        }
        let column = columns[text_index as usize];
        let bit = |value: u64, shift: i32| -> i64 {
            ((value.checked_shr(shift as u32).unwrap_or(0)) & 1) as i64
        };
        let delta_value = current_value - (1 - bit(column.d0, back_track_site));
        let mut min_value = delta_value;
        direction = 0;

        if back_track_site == 0 {
            let h_value = current_value - bit(column.horizontal_positive, 0)
                + bit(column.horizontal_negative, 0);
            if h_value < min_value {
                min_value = h_value;
                direction = 3;
            }
        } else if back_track_site == low_bound {
            let v_value = current_value - bit(column.positive_vertical, back_track_site - 1)
                + bit(column.negative_vertical, back_track_site - 1);
            if v_value < min_value {
                min_value = v_value;
                direction = 2;
            }
        } else {
            let h_value = current_value - bit(column.horizontal_positive, back_track_site)
                + bit(column.horizontal_negative, back_track_site);
            let v_value = current_value - bit(column.positive_vertical, back_track_site - 1)
                + bit(column.negative_vertical, back_track_site - 1);
            if v_value < min_value {
                min_value = v_value;
                direction = 2;
            }
            if h_value < min_value {
                min_value = h_value;
                direction = 3;
            }
        }

        if direction == 0 {
            if delta_value != current_value {
                direction = 1;
            }
            text_index -= 1;
            start_site -= 1;
        }
        if direction == 2 {
            back_track_site -= 1;
            start_site -= 1;
        } else if direction == 3 {
            text_index -= 1;
            back_track_site += 1;
        }

        path.push(direction as u8);
        current_value = min_value;
    }

    if text_index > 0 {
        path.resize(path.len() + text_index as usize, 0);
        start_site -= text_index;
        direction = 0;
    }
    if direction != 3 {
        start_site += 1;
    }

    BandedPath {
        end_site,
        start_site,
        error,
        path,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn padded(text: &str, length: usize) -> Vec<u8> {
        let mut buffer = vec![b'A'; length + 8];
        buffer[..text.len()].copy_from_slice(text.as_bytes());
        buffer
    }

    #[test]
    fn a_perfect_old_alignment_skips_traceback() {
        let pattern = padded("ACGTACGTACGG", 14);
        let text = padded("GTACGTACGG", 10);
        let result = reserve_banded_bpm_path(&pattern, 14, &text, 10, 2, 0, 11).unwrap();

        assert_eq!(result.end_site, 11);
        assert_eq!(result.start_site, 2);
        assert_eq!(result.error, 0);
        assert_eq!(result.path, vec![0; 10]);
    }

    #[test]
    fn an_ungapped_old_alignment_is_reused() {
        let pattern = padded("ACGTACGTACGG", 14);
        let text = padded("GTACGTACAG", 10);
        let result = reserve_banded_bpm_path(&pattern, 14, &text, 10, 2, 1, 11).unwrap();

        assert_eq!(result.end_site, 11);
        assert_eq!(result.error, 1);
        assert_eq!(result.path.len(), 10);
        assert_eq!(result.path.iter().filter(|code| **code == 1).count(), 1);
    }

    /// Golden values printed by the real `Reserve_Banded_BPM_PATH`, including
    /// mismatches, an indel (`300000000`) and two alignments that fail.
    #[test]
    fn matches_the_reference_traceback() {
        let cases: &[(&str, &str, usize, usize, usize, i32, i32, i32, u32, i32, &str)] = &[
            ("ACGTACGTACGG", "GTACGTACGG", 14usize, 10usize, 2usize, -1i32, -1i32, 11i32, 0u32, 2i32, "0000000000"),
            ("ACGTACGTACGG", "GTACGTACAG", 14usize, 10usize, 2usize, -1i32, -1i32, 11i32, 1u32, 2i32, "0100000000"),
            ("ACGTACGTACGG", "GTACGTACAA", 14usize, 10usize, 2usize, -1i32, -1i32, 11i32, 2u32, 2i32, "1100000000"),
            ("GCTATCCATAATTGGACCCCCTGCGATT", "CTATCCATAATTAGACCCCCTG", 28usize, 22usize, 3usize, -1i32, -1i32, 22i32, 1u32, 1i32, "0000000001000000000000"),
            ("TCTAGTAGGAGAGTTCGCTATTCAGGTTTAAC", "CTACTAGGAGAGTTCGCTATTC", 32usize, 22usize, 5usize, -1i32, -1i32, 22i32, 1u32, 1i32, "0000000000000000001000"),
            ("TTTTTGCGTGAGAGTACGT", "TGCGTGACAGTACGT", 19usize, 15usize, 2usize, -1i32, -1i32, 18i32, 1u32, 4i32, "000000010000000"),
            ("TCCTCTTTATCATACTCTCC", "TCTTTATCATACTCTC", 20usize, 16usize, 2usize, -1i32, -1i32, 18i32, 0u32, 3i32, "0000000000000000"),
            ("GTCGAGCAACACCGCTA", "GTCGAGGAGCACCGC", 17usize, 15usize, 1usize, -1i32, -1i32, -1i32, 4294967295u32, -7i32, ""),
            ("GGGAGCGAACCGGTTAAAACG", "GGGACCGAACCGGTT", 21usize, 15usize, 3usize, -1i32, -1i32, 14i32, 1u32, 0i32, "000000000010000"),
            ("CGACGCTTCGAGTCGC", "GACTCTTCGAGT", 16usize, 12usize, 2usize, -1i32, -1i32, 12i32, 1u32, 1i32, "000000001000"),
            ("GCGCTTGCATAGCCTGATGTCTCGGA", "GCGCTTGCCTAGCCTGAT", 26usize, 18usize, 4usize, -1i32, -1i32, 17i32, 1u32, 0i32, "000000000100000000"),
            ("CCTCTAATGGTACGACTGCTGAA", "TCTCATGGTAAGACTGCTGAA", 23usize, 21usize, 1usize, -1i32, -1i32, -1i32, 4294967295u32, -7i32, ""),
            ("TAGGAATTAGGGTGTCCCCAG", "AATCAGGGTGTCCCCAG", 21usize, 17usize, 2usize, -1i32, -1i32, 20i32, 1u32, 4i32, "00000000000001000"),
            ("TAATGATGTCGCCGTAAAT", "TGTCGCCGA", 19usize, 9usize, 5usize, -1i32, -1i32, 13i32, 1u32, 6i32, "300000000"),
            ("CGAGACATAGTGGTTTTCAACCCGAA", "GACATAGTGGTATTCAAC", 26usize, 18usize, 4usize, -1i32, -1i32, 20i32, 1u32, 3i32, "000000100000000000"),
        ];

        for (pattern_text, text_text, pattern_length, text_length, threshold, old_error,
             old_end_site, site, error, start, path_text) in cases
        {
            let pattern = padded(pattern_text, *pattern_length);
            let text = padded(text_text, *text_length);
            let result = reserve_banded_bpm_path(
                &pattern,
                *pattern_length,
                &text,
                *text_length,
                *threshold,
                *old_error,
                *old_end_site,
            );

            if *site == -1 {
                assert!(result.is_none(), "{text_text} must fail");
                continue;
            }
            let result = result.unwrap_or_else(|| panic!("{text_text} must align"));
            assert_eq!(result.end_site, *site, "end site of {text_text}");
            assert_eq!(result.error, *error, "error of {text_text}");
            assert_eq!(result.start_site, *start, "start site of {text_text}");
            let path: String = result.path.iter().map(|code| char::from(b'0' + *code)).collect();
            assert_eq!(path, *path_text, "path of {text_text}");
        }
    }

    /// Golden values printed by the real `Reserve_Banded_BPM_PATH` for a 200
    /// base text against the 232 base window that a boundary alignment uses.
    #[test]
    fn matches_the_reference_for_a_wide_perfect_window() {
        let pattern = vec![b'A'; 232];
        let text = vec![b'A'; 200];
        let result = reserve_banded_bpm_path(&pattern, 232, &text, 200, 16, -1, -1).unwrap();
        assert_eq!(result.end_site, 215);
        assert_eq!(result.start_site, 16);
        assert_eq!(result.error, 0);
        assert_eq!(result.path.len(), 200);

        let mut mutated = pattern.clone();
        mutated[96..106].fill(b'C');
        let result = reserve_banded_bpm_path(&mutated, 232, &text, 200, 16, -1, -1).unwrap();
        assert_eq!(result.end_site, 215);
        assert_eq!(result.start_site, 16);
        assert_eq!(result.error, 10);
    }
}


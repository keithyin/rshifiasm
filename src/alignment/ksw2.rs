//! Banded extension alignment, the counterpart of `ksw_extz2_sse` in
//! `ksw2_extz2_sse.c`.
//!
//! The original processes sixteen anti-diagonals per SIMD instruction; this
//! port computes the same recurrence one cell at a time while keeping the
//! eight-bit wrapping arithmetic, the `max_sc` clamp and the tie-breaking rules
//! that decide the traceback direction. Only cells inside the band are
//! computed: cells outside it are never read by cells inside it.

/// Negative infinity used for 32-bit scores.
pub const NEGATIVE_INFINITY: i32 = -0x4000_0000;

/// Do not record the alignment path.
pub const EXTZ_SCORE_ONLY: u32 = 0x01;
/// Right-align gaps.
pub const EXTZ_RIGHT: u32 = 0x02;
/// Use a full scoring matrix instead of match/mismatch.
pub const EXTZ_GENERIC_SCORE: u32 = 0x04;
/// Approximate the maximum score, faster but less precise.
pub const EXTZ_APPROXIMATE_MAXIMUM: u32 = 0x08;
/// Approximate Z-drop.
pub const EXTZ_APPROXIMATE_DROP: u32 = 0x10;
/// Only perform the extension, not the global alignment.
pub const EXTZ_EXTENSION_ONLY: u32 = 0x40;
/// Reverse the CIGAR in the output.
pub const EXTZ_REVERSE_CIGAR: u32 = 0x80;

/// Result of an extension alignment, `ksw_extz_t`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ExtensionResult {
    /// Whether the alignment stopped because of a Z-drop.
    pub zdropped: bool,
    /// Maximum score reached, `ez->max`.
    pub maximum: i32,
    /// Query coordinate of the maximum.
    pub maximum_query_position: i32,
    /// Target coordinate of the maximum.
    pub maximum_target_position: i32,
    /// Best score reaching the end of the query, `ez->mqe`.
    pub query_end_score: i32,
    /// Target coordinate of `query_end_score`.
    pub query_end_target_position: i32,
    /// Best score reaching the end of the target, `ez->mte`.
    pub target_end_score: i32,
    /// Query coordinate of `target_end_score`.
    pub target_end_query_position: i32,
    /// Score of the global alignment, may be [`NEGATIVE_INFINITY`].
    pub score: i32,
    /// Whether the query end was reached during an extension.
    pub reach_end: bool,
    /// BAM encoded CIGAR operations.
    pub cigar: Vec<u32>,
}

fn reset_result(result: &mut ExtensionResult) {
    result.maximum_query_position = -1;
    result.maximum_target_position = -1;
    result.query_end_target_position = -1;
    result.target_end_query_position = -1;
    result.maximum = 0;
    result.score = NEGATIVE_INFINITY;
    result.query_end_score = NEGATIVE_INFINITY;
    result.target_end_score = NEGATIVE_INFINITY;
    result.cigar.clear();
    result.zdropped = false;
    result.reach_end = false;
}

fn push_cigar(cigar: &mut Vec<u32>, operation: u32, length: i32) {
    match cigar.last_mut() {
        Some(last) if (*last & 0xf) == operation => *last += (length as u32) << 4,
        _ => cigar.push(((length as u32) << 4) | operation),
    }
}

/// `ksw_apply_zdrop`
fn apply_zdrop(
    result: &mut ExtensionResult,
    rotated: bool,
    score: i32,
    first: i32,
    second: i32,
    zdrop: i32,
    gap_extension: i8,
) -> bool {
    let (r, t) = if rotated {
        (first, second)
    } else {
        (first + second, first)
    };

    if score > result.maximum {
        result.maximum = score;
        result.maximum_target_position = t;
        result.maximum_query_position = r - t;
    } else if t >= result.maximum_target_position && r - t >= result.maximum_query_position {
        let target_shift = t - result.maximum_target_position;
        let query_shift = (r - t) - result.maximum_query_position;
        let diagonal_shift = if target_shift > query_shift {
            target_shift - query_shift
        } else {
            query_shift - target_shift
        };
        if zdrop >= 0 && result.maximum - score > zdrop + diagonal_shift * gap_extension as i32 {
            result.zdropped = true;
            return true;
        }
    }
    false
}

/// `ksw_backtrack` for the rotated (anti-diagonal) traceback matrix.
fn backtrack(
    traceback: &[u8],
    column_count: i32,
    row_offsets: &[i32],
    row_end_offsets: &[i32],
    start_row: i32,
    start_column: i32,
    reverse_cigar: bool,
    cigar: &mut Vec<u32>,
) {
    let mut row = start_row;
    let mut column = start_column;
    let mut state = 0i32;

    while row >= 0 && column >= 0 {
        let mut forced_state = -1i32;
        let diagonal = row + column;
        if row < row_offsets[diagonal as usize] {
            forced_state = 2;
        }
        if row > row_end_offsets[diagonal as usize] {
            forced_state = 1;
        }
        let cell = if forced_state < 0 {
            traceback[(diagonal as i64 * column_count as i64
                + (row - row_offsets[diagonal as usize]) as i64) as usize]
        } else {
            0
        };

        if state == 0 {
            state = (cell & 7) as i32;
        } else if (cell >> (state + 2)) & 1 == 0 {
            state = 0;
        }
        if state == 0 {
            state = (cell & 7) as i32;
        }
        if forced_state >= 0 {
            state = forced_state;
        }

        if state == 0 {
            push_cigar(cigar, 0, 1);
            row -= 1;
            column -= 1;
        } else if state == 1 {
            push_cigar(cigar, 2, 1);
            row -= 1;
        } else {
            push_cigar(cigar, 1, 1);
            column -= 1;
        }
    }

    if row >= 0 {
        push_cigar(cigar, 2, row + 1);
    }
    if column >= 0 {
        push_cigar(cigar, 1, column + 1);
    }

    if !reverse_cigar {
        cigar.reverse();
    }
}

/// `ksw_extz2_sse`: banded extension alignment.
///
/// `query` and `target` hold codes in `0..m`; `scoring_matrix` is a row-major
/// `m * m` matrix. A negative `band_width` disables the band.
#[allow(clippy::too_many_arguments)]
pub fn extz2(
    query_length: i32,
    query: &[u8],
    target_length: i32,
    target: &[u8],
    residue_type_count: i8,
    scoring_matrix: &[i8],
    gap_open: i8,
    gap_extension: i8,
    band_width: i32,
    zdrop: i32,
    end_bonus: i32,
    flags: u32,
    result: &mut ExtensionResult,
) {
    reset_result(result);

    let m = residue_type_count as i32;
    let qlen = query_length;
    let tlen = target_length;
    if m <= 0 || qlen <= 0 || tlen <= 0 {
        return;
    }

    let with_cigar = flags & EXTZ_SCORE_ONLY == 0;
    let approximate_maximum = flags & EXTZ_APPROXIMATE_MAXIMUM != 0;
    let gap_cost = gap_open as i32 + gap_extension as i32;
    let q = gap_open;
    let e = gap_extension;

    let w = if band_width < 0 {
        tlen.max(qlen)
    } else {
        band_width
    };
    let wl = w;
    let wr = w;

    let tlen_blocks = (tlen + 15) / 16;
    let n_col = {
        let shorter = qlen.min(tlen);
        let banded = shorter.min(w + 1);
        (banded + 15) / 16 + 1
    };

    let mut maximum_score = scoring_matrix[0] as i32;
    let mut minimum_score = scoring_matrix[1] as i32;
    for entry in scoring_matrix.iter().take((m * m) as usize).skip(1) {
        let value = *entry as i32;
        maximum_score = maximum_score.max(value);
        minimum_score = minimum_score.min(value);
    }
    if -minimum_score > 2 * gap_cost {
        // No mismatch can ever be observed, the original gives up immediately.
        return;
    }

    let match_score = scoring_matrix[0];
    let mismatch_score = scoring_matrix[1];
    let ambiguous_score = if scoring_matrix[(m * m - 1) as usize] == 0 {
        -gap_extension
    } else {
        scoring_matrix[(m * m - 1) as usize]
    };
    let wildcard = (m - 1) as u8;
    let score_ceiling = (scoring_matrix[0] as i32 + gap_cost * 2) as u8;

    // `qr` is the reversed query, `sf` the target.
    let mut reversed_query = vec![0u8; qlen as usize];
    for index in 0..qlen as usize {
        reversed_query[index] = query[qlen as usize - 1 - index];
    }
    let target_sequence = target;

    let storage = (tlen_blocks * 16 + 1) as usize;
    let mut u = vec![0u8; storage];
    let mut v = vec![0u8; storage];
    let mut x = vec![0u8; storage];
    let mut y = vec![0u8; storage];
    let mut s = vec![0u8; storage];

    let mut h_scores = if approximate_maximum {
        Vec::new()
    } else {
        vec![NEGATIVE_INFINITY; storage]
    };

    let traceback_columns = n_col as usize * 16;
    let row_count = (qlen + tlen - 1) as usize;
    let mut traceback = if with_cigar {
        vec![0u8; row_count * traceback_columns]
    } else {
        Vec::new()
    };
    let mut row_offsets = vec![0i32; row_count];
    let mut row_end_offsets = vec![0i32; row_count];

    let mut last_start: i32 = -1;
    let mut last_end: i32 = -1;
    let mut maximum_h0: i32 = 0;
    let mut last_h0_target: i32 = 0;

    for r in 0..(qlen + tlen - 1) {
        // Band boundaries for this anti-diagonal.
        let mut start = 0i32;
        let mut end = tlen - 1;
        if start < r - qlen + 1 {
            start = r - qlen + 1;
        }
        if end > r {
            end = r;
        }
        if start < (r - wr + 1) >> 1 {
            start = (r - wr + 1) >> 1;
        }
        if end > (r + wl) >> 1 {
            end = (r + wl) >> 1;
        }
        if start > end {
            result.zdropped = true;
            break;
        }

        let start_true = start;
        let end_true = end;

        // Boundary values carried in from the previous anti-diagonal.
        let (boundary_x, boundary_v) = if start_true > 0 {
            if start_true - 1 >= last_start && start_true - 1 <= last_end {
                (x[start_true as usize - 1], v[start_true as usize - 1])
            } else {
                (0u8, 0u8)
            }
        } else {
            (0u8, if r > 0 { q as u8 } else { 0 })
        };
        if end_true >= r {
            y[r as usize] = 0;
            u[r as usize] = if r > 0 { q as u8 } else { 0 };
        }

        // Substitution scores for the band. The query index of cell (r, t)
        // is r - t, which is always inside the query for cells in the band.
        for t in start_true..=end_true {
            let query_index = r - t;
            let query_code = if query_index >= 0 && query_index < qlen {
                reversed_query[(qlen - 1 - r + t) as usize]
            } else {
                0
            };
            let target_code = target_sequence[t as usize];
            s[t as usize] = if flags & EXTZ_GENERIC_SCORE != 0 {
                scoring_matrix[(target_code as i32 * m + query_code as i32) as usize] as u8
            } else if query_code == wildcard || target_code == wildcard {
                ambiguous_score as u8
            } else if query_code == target_code {
                match_score as u8
            } else {
                mismatch_score as u8
            };
        }

        let mut carry_x = boundary_x;
        let mut carry_v = boundary_v;

        for t in start_true..=end_true {
            let diagonal_x = carry_x;
            let diagonal_v = carry_v;
            carry_x = x[t as usize];
            carry_v = v[t as usize];

            let mut z = s[t as usize].wrapping_add((gap_cost * 2) as u8);
            let a = diagonal_x.wrapping_add(diagonal_v);
            let previous_u = u[t as usize];
            let b = y[t as usize].wrapping_add(previous_u);

            // The original compares the eight-bit lanes as *signed* bytes, so
            // the comparisons below go through `as i8` on purpose.
            let signed_a = a as i8;
            let signed_b = b as i8;
            let signed_z = z as i8;

            let mut direction = 0u8;
            if with_cigar {
                if flags & EXTZ_RIGHT == 0 {
                    direction = u8::from(signed_a > signed_z);
                } else {
                    direction = u8::from(signed_z <= signed_a);
                }
            }

            z = z.max(a);
            if with_cigar {
                let signed_updated_z = z as i8;
                if flags & EXTZ_RIGHT == 0 {
                    if signed_b > signed_updated_z {
                        direction = 2;
                    }
                } else if signed_updated_z <= signed_b {
                    direction = 2;
                }
            }

            z = z.max(b);
            z = z.min(score_ceiling);

            u[t as usize] = z.wrapping_sub(diagonal_v);
            v[t as usize] = z.wrapping_sub(previous_u);
            let z_adjusted = z.wrapping_sub(q as u8);
            let a_adjusted = a.wrapping_sub(z_adjusted);
            let b_adjusted = b.wrapping_sub(z_adjusted);

            if flags & EXTZ_RIGHT == 0 {
                // d = a > 0 ? 8 : 0, keeping a positive `a`.
                let a_positive = (a_adjusted as i8) > 0;
                x[t as usize] = if a_positive { a_adjusted } else { 0 };
                let b_positive = (b_adjusted as i8) > 0;
                y[t as usize] = if b_positive { b_adjusted } else { 0 };
                if with_cigar {
                    if a_positive {
                        direction |= 0x08;
                    }
                    if b_positive {
                        direction |= 0x10;
                    }
                }
            } else {
                // d = 0 > a ? 0 : 8, keeping `a` unless it is negative.
                let a_negative = (a_adjusted as i8) < 0;
                x[t as usize] = if a_negative { 0 } else { a_adjusted };
                let b_negative = (b_adjusted as i8) < 0;
                y[t as usize] = if b_negative { 0 } else { b_adjusted };
                if with_cigar {
                    if !a_negative {
                        direction |= 0x08;
                    }
                    if !b_negative {
                        direction |= 0x10;
                    }
                }
            }
            if with_cigar {
                let offset =
                    r as i64 * traceback_columns as i64 + (t - start_true) as i64;
                traceback[offset as usize] = direction;
            }
        }

        row_offsets[r as usize] = start_true;

        if !approximate_maximum {
            let (max_h, max_t);
            if r > 0 {
                // The last cell of the diagonal needs special casing because
                // its vertical difference is not available yet.
                let mut running_max = if end_true > 0 {
                    h_scores[end_true as usize - 1] + u[end_true as usize] as i32 - gap_cost
                } else {
                    h_scores[end_true as usize] + v[end_true as usize] as i32 - gap_cost
                };
                h_scores[end_true as usize] = running_max;
                let mut running_index = end_true;

                for t in start_true..end_true {
                    h_scores[t as usize] += v[t as usize] as i32 - gap_cost;
                    if h_scores[t as usize] > running_max {
                        running_max = h_scores[t as usize];
                        running_index = t;
                    }
                }
                max_h = running_max;
                max_t = running_index;
            } else {
                h_scores[0] = v[0] as i32 - gap_cost - gap_cost;
                max_h = h_scores[0];
                max_t = 0;
            }

            if end_true == tlen - 1 && h_scores[end_true as usize] > result.target_end_score {
                result.target_end_score = h_scores[end_true as usize];
                result.target_end_query_position = r - end_true;
            }
            if r - start_true == qlen - 1 && h_scores[start_true as usize] > result.query_end_score
            {
                result.query_end_score = h_scores[start_true as usize];
                result.query_end_target_position = start_true;
            }
            if apply_zdrop(result, true, max_h, r, max_t, zdrop, e) {
                break;
            }
            if r == qlen + tlen - 2 && end_true == tlen - 1 {
                result.score = h_scores[(tlen - 1) as usize];
            }
        } else {
            if r > 0 {
                if last_h0_target >= start_true
                    && last_h0_target <= end_true
                    && last_h0_target + 1 >= start_true
                    && last_h0_target + 1 <= end_true
                {
                    let first_delta = v[last_h0_target as usize] as i32 - gap_cost;
                    let second_delta = u[last_h0_target as usize + 1] as i32 - gap_cost;
                    if first_delta > second_delta {
                        maximum_h0 += first_delta;
                    } else {
                        maximum_h0 += second_delta;
                        last_h0_target += 1;
                    }
                } else if last_h0_target >= start_true && last_h0_target <= end_true {
                    maximum_h0 += v[last_h0_target as usize] as i32 - gap_cost;
                } else {
                    last_h0_target += 1;
                    maximum_h0 += u[last_h0_target as usize] as i32 - gap_cost;
                }
                if flags & EXTZ_APPROXIMATE_DROP != 0
                    && apply_zdrop(result, true, maximum_h0, r, last_h0_target, zdrop, e)
                {
                    break;
                }
            } else {
                maximum_h0 = v[0] as i32 - gap_cost - gap_cost;
                last_h0_target = 0;
            }
            if r == qlen + tlen - 2 && end_true == tlen - 1 {
                result.score = maximum_h0;
            }
        }

        // `off_end` records the padded end of the traceback row.
        row_end_offsets[r as usize] = (end_true + 16) / 16 * 16 - 1;
        last_start = start_true;
        last_end = end_true;
    }

    if with_cigar {
        let reverse_cigar = flags & EXTZ_REVERSE_CIGAR != 0;
        if !result.zdropped && flags & EXTZ_EXTENSION_ONLY == 0 {
            backtrack(
                &traceback,
                traceback_columns as i32,
                &row_offsets,
                &row_end_offsets,
                tlen - 1,
                qlen - 1,
                reverse_cigar,
                &mut result.cigar,
            );
        } else if !result.zdropped
            && flags & EXTZ_EXTENSION_ONLY != 0
            && result.query_end_score + end_bonus > result.maximum
        {
            result.reach_end = true;
            backtrack(
                &traceback,
                traceback_columns as i32,
                &row_offsets,
                &row_end_offsets,
                result.query_end_target_position,
                qlen - 1,
                reverse_cigar,
                &mut result.cigar,
            );
        } else if result.maximum_target_position >= 0 && result.maximum_query_position >= 0 {
            backtrack(
                &traceback,
                traceback_columns as i32,
                &row_offsets,
                &row_end_offsets,
                result.maximum_target_position,
                result.maximum_query_position,
                reverse_cigar,
                &mut result.cigar,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Simple 5x5 matrix: match 2, mismatch -4, ambiguous 0.
    const SIMPLE_MATRIX: [i8; 25] = [
        2, -4, -4, -4, 0, //
        -4, 2, -4, -4, 0, //
        -4, -4, 2, -4, 0, //
        -4, -4, -4, 2, 0, //
        0, 0, 0, 0, 0,
    ];

    fn align(query: &[u8], target: &[u8]) -> ExtensionResult {
        let mut result = ExtensionResult::default();
        extz2(
            query.len() as i32,
            query,
            target.len() as i32,
            target,
            5,
            &SIMPLE_MATRIX,
            4,
            2,
            -1,
            -1,
            0,
            0,
            &mut result,
        );
        result
    }

    #[test]
    fn identical_sequences_align_end_to_end() {
        let result = align(&[0, 1, 2, 3], &[0, 1, 2, 3]);
        assert_eq!(result.score, 8);
        assert!(!result.cigar.is_empty());
        let total: u32 = result.cigar.iter().map(|entry| entry >> 4).sum();
        assert_eq!(total, 4);
    }

    #[test]
    fn single_mismatch_is_reported() {
        let result = align(&[0, 1, 2, 3], &[0, 1, 3, 3]);
        assert_eq!(result.score, 2 + 2 - 4 + 2);
    }

    #[allow(dead_code)]
    fn cigar_text(cigar: &[u32]) -> String {
        let operations = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B'];
        cigar
            .iter()
            .map(|entry| format!("{}{}", entry >> 4, operations[(entry & 0xf) as usize]))
            .collect()
    }

    fn align_with(query: &[u8], target: &[u8], band_width: i32) -> ExtensionResult {
        let mut result = ExtensionResult::default();
        extz2(
            query.len() as i32,
            query,
            target.len() as i32,
            target,
            5,
            &SIMPLE_MATRIX,
            4,
            2,
            band_width,
            -1,
            0,
            0,
            &mut result,
        );
        result
    }

    /// Straightforward global affine alignment used as an independent oracle
    /// for the score of unbounded (`band_width < 0`) alignments.
    fn reference_global_score(query: &[u8], target: &[u8]) -> i32 {
        const GAP_OPEN: i32 = 4;
        const GAP_EXTENSION: i32 = 2;
        const NEGATIVE: i32 = -1_000_000;

        let rows = query.len() + 1;
        let columns = target.len() + 1;
        let mut h = vec![vec![NEGATIVE; columns]; rows];
        let mut e = vec![vec![NEGATIVE; columns]; rows];
        let mut f = vec![vec![NEGATIVE; columns]; rows];
        h[0][0] = 0;
        for row in 1..rows {
            e[row][0] = -GAP_OPEN as i32 - row as i32 * GAP_EXTENSION as i32;
            h[row][0] = e[row][0];
        }
        for column in 1..columns {
            f[0][column] = -GAP_OPEN as i32 - column as i32 * GAP_EXTENSION as i32;
            h[0][column] = f[0][column];
        }

        for row in 1..rows {
            for column in 1..columns {
                // A wildcard takes precedence over equality and is scored as
                // `-e` when `mat[m*m-1] == 0`, which is the case here.
                let score = if query[row - 1] == 4 || target[column - 1] == 4 {
                    -GAP_EXTENSION
                } else if query[row - 1] == target[column - 1] {
                    2
                } else {
                    -4
                };
                // ksw2 charges `gapo + l * gape` for a gap of length `l`, so
                // the opening step already includes one extension.
                e[row][column] = (h[row - 1][column] - GAP_OPEN - GAP_EXTENSION).max(
                    e[row - 1][column] - GAP_EXTENSION,
                );
                f[row][column] = (h[row][column - 1] - GAP_OPEN - GAP_EXTENSION).max(
                    f[row][column - 1] - GAP_EXTENSION,
                );
                h[row][column] = (h[row - 1][column - 1] + score)
                    .max(e[row][column])
                    .max(f[row][column]);
            }
        }

        h[rows - 1][columns - 1]
    }

    #[test]
    fn score_agrees_with_a_full_matrix_alignment() {
        let sequences: Vec<Vec<u8>> = vec![
            vec![0, 1, 2, 3],
            vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1],
            vec![0, 1, 2, 0, 1, 2, 3, 0, 1, 2],
            vec![3, 3, 1, 1, 0, 2],
            vec![0, 4, 2, 3, 4, 1, 2, 3, 4, 1],
        ];

        for query in &sequences {
            for target in &sequences {
                let result = align_with(query, target, -1);
                assert_eq!(
                    result.score,
                    reference_global_score(query, target),
                    "score of {:?} against {:?}",
                    query,
                    target
                );
            }
        }
    }

    #[test]
    fn empty_input_is_ignored() {
        let mut result = ExtensionResult::default();
        extz2(0, &[], 0, &[], 5, &SIMPLE_MATRIX, 4, 2, -1, -1, 0, 0, &mut result);
        assert_eq!(result.cigar.len(), 0);
        assert_eq!(result.score, NEGATIVE_INFINITY);
    }
}

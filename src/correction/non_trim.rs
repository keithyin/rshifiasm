//! `verify_sub_window` and `non_trim_error_rate` from `Correct.cpp`.
//!
//! A window that the first pass could not align leaves a hole in the overlap.
//! Before the overlap is classified, every hole is re-aligned from both of its
//! sides: the left attempt starts at the end of the previous aligned window and
//! walks forwards, the right attempt starts at the beginning of the next one
//! and walks backwards. A hole that neither attempt covers contributes its full
//! length to the error budget, which is what makes a divergent region visible
//! in the final error rate.

use crate::alignment::banded_bpm::{NO_ALIGNMENT_ERROR, alignment_extension};
use crate::correction::recalculate::{double_error_threshold, initial_error_threshold};
use crate::correction::region::{THRESHOLD, determine_overlap_region, fill_subregion};
use crate::correction::window::{number_of_windows, window_id_by_end, window_id_by_start};
use crate::overlap::region::OverlapRegion;
use crate::correction::target::TargetSource;
use crate::reads::read_store::ReadStore;

/// Outcome of `verify_sub_window`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SubWindowAlignment {
    /// Edit distance, [`NO_ALIGNMENT_ERROR`] when nothing aligned.
    pub error: u32,
    /// Target end of the alignment, `-1` on failure.
    pub y_end: i32,
    /// Query end of the alignment, `-1` on failure.
    pub x_end: i32,
    /// Query bases the alignment covers, `0` on failure.
    pub aligned_length: i64,
}

impl SubWindowAlignment {
    fn failure() -> Self {
        Self {
            error: NO_ALIGNMENT_ERROR,
            y_end: -1,
            x_end: -1,
            aligned_length: 0,
        }
    }
}

/// `verify_sub_window`: maps the query window onto the target and aligns it
/// from one of its ends.
///
/// `alignment_strand` picks the direction: `0` walks forwards from the start of
/// both windows, anything else walks backwards from their ends. This is the
/// HiFi branch of the reference; the ultra-long variant reads the unitig
/// instead of a read.
#[allow(clippy::too_many_arguments)]
pub fn verify_sub_window(
    reads: &ReadStore,
    query_sequence: &[u8],
    x_begin: i64,
    x_length: i64,
    y_begin: i64,
    y_length: i64,
    y_identifier: usize,
    y_strand: u8,
    threshold: i64,
    alignment_strand: i32,
) -> SubWindowAlignment {
    let target_length = reads.read_length(y_identifier) as i64;
    let Some(region) = determine_overlap_region(threshold, y_begin, target_length, y_length) else {
        return SubWindowAlignment::failure();
    };

    let target = fill_subregion(
        &TargetSource::Reads(reads),
        y_identifier,
        region.start,
        region.length,
        y_strand,
        region.extra_begin,
        region.extra_end,
    );

    let query_start = x_begin.max(0) as usize;
    let query_end = (query_start + x_length as usize).min(query_sequence.len());
    if query_start >= query_end {
        return SubWindowAlignment::failure();
    }
    let query = &query_sequence[query_start..query_end];

    let alignment = alignment_extension(
        &target,
        y_length as usize,
        query,
        x_length as usize,
        threshold as usize,
        alignment_strand,
    );

    SubWindowAlignment {
        error: alignment.error,
        y_end: alignment.pattern_end,
        x_end: alignment.text_end,
        aligned_length: alignment.aligned_text_length as i64,
    }
}

/// The parameters `get_init_paras` derives for one hole of an overlap.
struct HoleParameters {
    /// Target position the forward attempt starts at.
    y_start: i64,
    /// `N` bases padded in front of the target window.
    extra_begin: i32,
    /// Single-window edit budget before it is doubled.
    threshold: i64,
}

/// `get_init_paras`: where an unaligned window should be looked up on the
/// target, and how many errors it may have.
///
/// The target position comes from the chain gap structure (`f_cigar`), not from
/// the neighbouring windows, so it is only a starting guess.
fn initial_parameters(
    region: &OverlapRegion,
    reads: &ReadStore,
    x_start: i64,
    x_end: i64,
    error_rate: f64,
    block_length: i64,
) -> Option<HoleParameters> {
    let length = x_end + 1 - x_start;
    let threshold = initial_error_threshold(length, error_rate, block_length, THRESHOLD);

    let mut y_start = (x_start - region.x_position_start as i64) + region.y_position_start as i64;
    y_start += region.fake_cigar.y_start_offset(x_start)?;

    let window_length = length + (threshold << 1);
    let target_length = reads.read_length(region.y_identifier as usize) as i64;
    let mapped = determine_overlap_region(threshold, y_start, target_length, window_length)?;

    Some(HoleParameters {
        y_start: mapped.start,
        extra_begin: mapped.extra_begin,
        threshold,
    })
}

/// The error a hole contributes when both of its ends were attempted.
///
/// When both attempts cover it, their errors are scaled by how much of the hole
/// they cover together; the reference does that scaling in single precision, so
/// the same width is used here.
fn hole_error(left: &SubWindowAlignment, right: &SubWindowAlignment, x_length: i64) -> i64 {
    if left.aligned_length != 0 && right.aligned_length != 0 {
        let covered = left.aligned_length + right.aligned_length;
        if covered <= x_length {
            left.error as i64 + right.error as i64 + (x_length - covered)
        } else {
            let rate = x_length as f32 / covered as f32;
            ((left.error + right.error) as f32 * rate) as i64
        }
    } else if left.aligned_length == 0 && right.aligned_length == 0 {
        x_length
    } else if left.aligned_length != 0 {
        left.error as i64 + (x_length - left.aligned_length)
    } else {
        right.error as i64 + (x_length - right.aligned_length)
    }
}

/// `non_trim_error_rate`: the error rate of an overlap after adding the holes
/// that the window pass could not align.
///
/// Unlike the per-window rate this counts the whole overlap, including the
/// regions before the first and after the last aligned window. It is the value
/// compared with `max_ov_diff_final` to decide whether an overlap is a match.
pub fn non_trim_error_rate(
    region: &OverlapRegion,
    reads: &ReadStore,
    query_sequence: &[u8],
    error_rate: f64,
    block_length: i64,
) -> f64 {
    let windows = &region.window_cigar.windows;
    let aligned_windows = windows.len() as i64;
    let window_count = number_of_windows(
        region.x_position_start as i64,
        region.x_position_end as i64 + 1,
        block_length,
    );
    debug_assert!(window_count >= aligned_windows && aligned_windows > 0);

    let overlap_start = region.x_position_start as i32;
    let overlap_end = region.x_position_end as i32;
    let y_identifier = region.y_identifier as usize;
    let y_strand = region.y_strand as u8;

    let mut total_error: i64 = 0;
    let mut total_length: i64 = 0;

    // `index_end` is the identifier of the next aligned window to the right. A
    // hole is the range of windows between two aligned ones, so the inner loop
    // only runs when their identifiers are not consecutive.
    let mut index_end = window_count;
    let mut index = aligned_windows - 1;
    while index >= 0 {
        let window = windows[index as usize];
        let (window_identifier, _) =
            window_id_by_end(overlap_start, window.x_end as i64, block_length);

        total_length += window.x_end as i64 + 1 - window.x_start as i64;
        total_error += window.error as i64;

        let mut hole = window_identifier + 1;
        let mut hole_end = window.x_end as i64;
        while hole < index_end {
            let hole_start = hole_end + 1;
            let (hole_identifier, computed_end) =
                window_id_by_start(overlap_start, overlap_end, hole_start, block_length);
            hole_end = computed_end;
            debug_assert_eq!(hole_identifier, hole);
            let x_length = hole_end + 1 - hole_start;
            total_length += x_length;

            let Some(parameters) =
                initial_parameters(region, reads, hole_start, hole_end, error_rate, block_length)
            else {
                total_error += x_length;
                hole += 1;
                continue;
            };

            let single_threshold = parameters.threshold;
            let threshold = double_error_threshold(single_threshold, x_length);
            let window_length = x_length + (threshold << 1);

            // Prefer continuing from the aligned window on either side; the
            // side that is missing falls back to the chain estimate.
            let mut y_begin_left = -1i64;
            let mut y_begin_right = -1i64;
            if hole == window_identifier + 1 {
                y_begin_left = windows[index as usize].y_end as i64 + 1;
            }
            if hole + 1 == index_end && index + 1 < aligned_windows {
                y_begin_right = windows[(index + 1) as usize].y_start as i64 - x_length;
            }
            if y_begin_left == -1 && y_begin_right == -1 {
                y_begin_left = parameters.y_start;
                if parameters.extra_begin >= 0 {
                    y_begin_left += single_threshold - parameters.extra_begin as i64;
                }
                y_begin_right = y_begin_left;
            }
            if y_begin_left == -1 && y_begin_right != -1 {
                y_begin_left = y_begin_right;
            }
            if y_begin_right == -1 && y_begin_left != -1 {
                y_begin_right = y_begin_left;
            }

            let left = attempt(
                reads,
                query_sequence,
                hole_start,
                x_length,
                y_begin_left,
                window_length,
                y_identifier,
                y_strand,
                threshold,
                0,
            );
            let right = attempt(
                reads,
                query_sequence,
                hole_start,
                x_length,
                y_begin_right,
                window_length,
                y_identifier,
                y_strand,
                threshold,
                1,
            );

            total_error += hole_error(&left, &right, x_length);
            hole += 1;
        }

        index_end = window_identifier;
        index -= 1;
    }

    // Windows before the first aligned one: there is no left neighbour, so the
    // only anchor is the first aligned window on the right.
    if index_end > 0 {
        let mut hole = 0i64;
        let mut hole_end = region.x_position_start as i64 - 1;
        while hole < index_end {
            let hole_start = hole_end + 1;
            let (hole_identifier, computed_end) =
                window_id_by_start(overlap_start, overlap_end, hole_start, block_length);
            hole_end = computed_end;
            debug_assert_eq!(hole_identifier, hole);
            let x_length = hole_end + 1 - hole_start;
            total_length += x_length;

            let Some(parameters) =
                initial_parameters(region, reads, hole_start, hole_end, error_rate, block_length)
            else {
                total_error += x_length;
                hole += 1;
                continue;
            };

            let single_threshold = parameters.threshold;
            let threshold = double_error_threshold(single_threshold, x_length);
            let window_length = x_length + (threshold << 1);

            let mut y_begin_left = -1i64;
            let mut y_begin_right = -1i64;
            // The reference still holds `index == -1` here, so `windows[index
            // + 1]` is the first aligned window; that anchor always exists.
            if hole + 1 == index_end && 0 < aligned_windows {
                y_begin_right = windows[0].y_start as i64 - x_length;
            }
            if y_begin_left == -1 && y_begin_right == -1 {
                y_begin_left = parameters.y_start;
                if parameters.extra_begin >= 0 {
                    y_begin_left += single_threshold - parameters.extra_begin as i64;
                }
                y_begin_right = y_begin_left;
            }
            if y_begin_left == -1 && y_begin_right != -1 {
                y_begin_left = y_begin_right;
            }
            if y_begin_right == -1 && y_begin_left != -1 {
                y_begin_right = y_begin_left;
            }

            let left = attempt(
                reads,
                query_sequence,
                hole_start,
                x_length,
                y_begin_left,
                window_length,
                y_identifier,
                y_strand,
                threshold,
                0,
            );
            let right = attempt(
                reads,
                query_sequence,
                hole_start,
                x_length,
                y_begin_right,
                window_length,
                y_identifier,
                y_strand,
                threshold,
                1,
            );

            total_error += hole_error(&left, &right, x_length);
            hole += 1;
        }
    }

    debug_assert_eq!(
        total_length,
        region.x_position_end as i64 + 1 - region.x_position_start as i64
    );
    total_error as f64 / total_length as f64
}

/// Calls [`verify_sub_window`] for one side, or reports the failure the
/// reference would see when that side has no anchor.
#[allow(clippy::too_many_arguments)]
fn attempt(
    reads: &ReadStore,
    query_sequence: &[u8],
    x_begin: i64,
    x_length: i64,
    y_begin: i64,
    y_length: i64,
    y_identifier: usize,
    y_strand: u8,
    threshold: i64,
    alignment_strand: i32,
) -> SubWindowAlignment {
    if y_begin == -1 {
        return SubWindowAlignment::failure();
    }
    verify_sub_window(
        reads,
        query_sequence,
        x_begin,
        x_length,
        y_begin,
        y_length,
        y_identifier,
        y_strand,
        threshold,
        alignment_strand,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::correction::window::WindowCigar;
    use crate::overlap::region::FakeCigar;

    /// Two reads: the query is a thousand `A`s, the target is supplied by the
    /// caller and shares that length.
    fn reads_with_target(target: &[u8]) -> ReadStore {
        let mut reads = ReadStore::default();
        reads.push_read("query", &vec![b'A'; 1000]);
        reads.push_read("target", target);
        reads
    }

    fn region_of(windows: &[(i32, i32, i32, i32)], target_end: i32) -> OverlapRegion {
        let mut window_cigar = WindowCigar::default();
        for (x_start, x_end, y_start, y_end) in windows {
            window_cigar.append_window(
                *x_start,
                *x_end,
                *y_start,
                *y_end,
                0,
                0,
                0,
                THRESHOLD as i16,
            );
        }
        let mut fake_cigar = FakeCigar::default();
        fake_cigar.push(0, 0);
        fake_cigar.push(999, 0);
        OverlapRegion {
            x_identifier: 0,
            x_position_start: 0,
            x_position_end: 999,
            y_identifier: 1,
            y_position_start: windows.first().map(|window| window.2 as u32).unwrap_or(0),
            y_position_end: target_end as u32,
            y_strand: 0,
            window_cigar,
            fake_cigar,
            ..OverlapRegion::default()
        }
    }

    /// Golden values printed by the real `non_trim_error_rate`.
    #[test]
    fn matches_the_reference_error_rate() {
        let query = vec![b'A'; 1000];

        // A: three contiguous aligned windows; nothing is re-verified.
        let reads = reads_with_target(&vec![b'A'; 1000]);
        let region = region_of(
            &[
                (0, 374, 0, 374),
                (375, 749, 375, 749),
                (750, 999, 750, 999),
            ],
            999,
        );
        assert_eq!(non_trim_error_rate(&region, &reads, &query, 0.04, 375), 0.0);

        // B: the target ends inside the second window, so the chain estimate of
        // the tail window cannot be mapped and the whole hole is charged.
        let reads = reads_with_target(&vec![b'A'; 600]);
        let region = region_of(&[(0, 374, 0, 374), (375, 749, 500, 874)], 599);
        assert_eq!(non_trim_error_rate(&region, &reads, &query, 0.04, 375), 0.25);

        // D: one window missing in the middle, re-verified from both sides
        // against a target with twenty mismatching bases.
        let mut target = vec![b'A'; 1000];
        target[500..520].fill(b'C');
        let reads = reads_with_target(&target);
        let region = region_of(&[(0, 374, 0, 374), (750, 999, 750, 999)], 999);
        assert_eq!(
            non_trim_error_rate(&region, &reads, &query, 0.04, 375),
            0.02
        );

        // E: block size 100 with the sixth window missing; every base matches,
        // so the hole costs nothing even though it is re-verified.
        let reads = reads_with_target(&vec![b'A'; 1000]);
        let windows: Vec<(i32, i32, i32, i32)> = (0..10)
            .filter(|index| *index != 5)
            .map(|index| {
                let start = index * 100;
                (start, start + 99, start, start + 99)
            })
            .collect();
        let region = region_of(&windows, 999);
        assert_eq!(non_trim_error_rate(&region, &reads, &query, 0.04, 100), 0.0);

        // F: the same hole, but it covers sixty mismatches, so neither attempt
        // can align through it.
        let mut target = vec![b'A'; 1000];
        target[500..560].fill(b'C');
        let reads = reads_with_target(&target);
        let windows: Vec<(i32, i32, i32, i32)> = (0..10)
            .filter(|index| *index != 5)
            .map(|index| {
                let start = index * 100;
                (start, start + 99, start, start + 99)
            })
            .collect();
        let region = region_of(&windows, 999);
        assert_eq!(non_trim_error_rate(&region, &reads, &query, 0.04, 100), 0.03);
    }

    #[test]
    fn a_window_is_verified_from_both_ends() {
        let mut target = vec![b'A'; 1000];
        target[500..520].fill(b'C');
        let reads = reads_with_target(&target);
        let query = vec![b'A'; 1000];

        // Golden values from the reference: forwards ends at pattern index 405
        // and reports query end 374, backwards ends at index 32 and reports
        // query end 1 (both are translated from the reversed coordinates).
        let forward = verify_sub_window(&reads, &query, 375, 375, 375, 437, 1, 0, 31, 0);
        assert_eq!(forward.error, 20);
        assert_eq!((forward.y_end, forward.x_end, forward.aligned_length), (405, 374, 375));

        let backward = verify_sub_window(&reads, &query, 375, 375, 375, 437, 1, 0, 31, 1);
        assert_eq!(backward.error, 20);
        assert_eq!((backward.y_end, backward.x_end, backward.aligned_length), (32, 1, 375));
    }

    #[test]
    fn a_window_outside_the_target_is_rejected() {
        let reads = reads_with_target(&vec![b'A'; 300]);
        let query = vec![b'A'; 1000];
        let result = verify_sub_window(&reads, &query, 375, 375, 800, 437, 1, 0, 31, 0);
        assert_eq!(result.error, NO_ALIGNMENT_ERROR);
        assert_eq!(result.aligned_length, 0);
        assert_eq!(result.y_end, -1);
    }
}

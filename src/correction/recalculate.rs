//! Threshold helpers and the final overlap classification of
//! `recalcate_window_advance`, plus `append_unmatched_wins`.
//!
//! The window cigar is a sparse array indexed by window identifier. Once an
//! overlap is accepted, the missing windows are filled with sentinel entries so
//! that the consensus builder can walk the whole overlap by index.

use crate::correction::region::{THRESHOLD_MAXIMUM, adjust_threshold};
use crate::correction::window::{WindowCigar, WindowList, number_of_windows, window_id_by_end,
                                window_id_by_start};
use crate::overlap::region::OverlapRegion;

/// `OVERLAP_THRESHOLD_HIFI_FILTER`: fraction of an overlap that must have been
/// aligned before it is even considered for an exact cigar.
pub const OVERLAP_HIFI_ALIGNMENT_FRACTION: f64 = 0.9;

/// `MIN_UL_ALIN_RATE` / `MIN_UL_ALIN_LEN`.
pub const MINIMUM_ULTRA_LONG_ALIGNMENT_RATE: f64 = 0.5;
pub const MINIMUM_ULTRA_LONG_ALIGNMENT_LENGTH: i64 = 6 * 75;

/// `get_init_err_thres`: the initial edit budget of a window.
///
/// A full-length window gets the fixed `block_error`; a shorter one derives its
/// budget from the error rate.
pub fn initial_error_threshold(
    length: i64,
    error_rate: f64,
    block_size: i64,
    block_error: i64,
) -> i64 {
    if length >= block_size {
        return block_error;
    }
    let mut threshold = (length as f64 * error_rate) as i64;
    threshold = adjust_threshold(threshold, length);
    if threshold > THRESHOLD_MAXIMUM {
        threshold = THRESHOLD_MAXIMUM;
    }
    threshold
}

/// `double_error_threshold`: doubles the budget for a second alignment
/// attempt, with a floor for long windows and a ceiling for all windows.
pub fn double_error_threshold(previous_threshold: i64, length: i64) -> i64 {
    let adjusted = adjust_threshold(previous_threshold, length);
    let mut threshold = adjusted * 2;
    if length >= 300 && threshold < THRESHOLD_MAXIMUM {
        threshold = THRESHOLD_MAXIMUM;
    }
    if threshold > THRESHOLD_MAXIMUM {
        threshold = THRESHOLD_MAXIMUM;
    }
    threshold
}

/// `simi_pass`: whether the aligned part of an overlap is long enough.
///
/// `error_rate` is `Some` for the ultra-long branch, `None` for HiFi; in the
/// latter case `alignment_fraction` decides, unless `second_check` is set.
pub fn similarity_pass(
    overlap_length: i64,
    aligned_length: i64,
    second_check: bool,
    alignment_fraction: f64,
    error_rate: Option<f64>,
) -> bool {
    if aligned_length == 0 || overlap_length == 0 {
        return false;
    }
    if !second_check && error_rate.is_none() {
        overlap_length as f64 * alignment_fraction <= aligned_length as f64
    } else if let Some(rate) = error_rate {
        overlap_length as f64 * (1.0 - rate) <= aligned_length as f64
    } else if second_check {
        overlap_length as f64 * MINIMUM_ULTRA_LONG_ALIGNMENT_RATE <= aligned_length as f64
            && aligned_length >= MINIMUM_ULTRA_LONG_ALIGNMENT_LENGTH
    } else {
        false
    }
}

/// One sentinel window, the values `append_unmatched_wins` writes for a window
/// that no alignment covered.
fn unmatched_window(x_start: i32, x_end: i32) -> WindowList {
    WindowList {
        x_start,
        x_end,
        y_start: -1,
        y_end: -1,
        extra_begin: 0,
        extra_end: 0,
        error: -1,
        error_threshold: -1,
        cigar_start: 0,
        cigar_length: 0,
    }
}

/// `append_unmatched_wins`: expands a sparse window cigar to one entry per
/// window identifier, keeping the verified windows at their own identifier and
/// filling the gaps with sentinels.
///
/// `region_start` and `region_end` are the inclusive query bounds of the
/// overlap, and `block_length` is the window length.
pub fn append_unmatched_windows(
    region_start: i32,
    region_end: i32,
    window_cigar: &mut WindowCigar,
    block_length: i64,
) {
    let expected = number_of_windows(
        region_start as i64,
        region_end as i64 + 1,
        block_length,
    ) as usize;
    let existing = window_cigar.windows.len();
    if expected <= existing {
        return;
    }

    // Keep a copy: the algorithm scatters the verified windows into their
    // identifier slots while reading the original list from the back.
    let original = window_cigar.windows.clone();
    let mut slots = vec![unmatched_window(0, 0); expected];
    let mut next_empty_end = expected;

    for index in (0..existing).rev() {
        let window = original[index];
        let (identifier, _) = window_id_by_end(
            region_start,
            window.x_end as i64,
            block_length,
        );
        let identifier = identifier as usize;
        if identifier >= expected {
            continue;
        }

        slots[identifier] = window;
        let mut window_end = window.x_end as i64;
        for slot in identifier + 1..next_empty_end {
            let window_start = window_end + 1;
            let (_, computed_end) = window_id_by_start(
                region_start,
                region_end,
                window_start,
                block_length,
            );
            slots[slot] = unmatched_window(window_start as i32, computed_end as i32);
            window_end = computed_end;
        }
        next_empty_end = identifier;
    }

    if next_empty_end > 0 {
        let mut window_end = region_start as i64 - 1;
        for slot in 0..next_empty_end {
            let window_start = window_end + 1;
            let (_, computed_end) = window_id_by_start(
                region_start,
                region_end,
                window_start,
                block_length,
            );
            slots[slot] = unmatched_window(window_start as i32, computed_end as i32);
            window_end = computed_end;
        }
    }

    window_cigar.windows = slots;
}

/// The final `is_match` decision of `recalcate_window_advance`, once the
/// non-trimmed error rate of the overlap is known.
///
/// `1` means the overlap is a real match, `3` means it is kept as a weak match,
/// `0` means it is discarded.
pub fn classify_overlap(
    region: &mut OverlapRegion,
    error_rate: f64,
    final_error_rate: f64,
    block_length: i64,
) -> u8 {
    let overlap_length = region.x_position_end as i64 + 1 - region.x_position_start as i64;
    let aligned_length = region.alignment_length as i64;

    if overlap_length as f64 * OVERLAP_HIFI_ALIGNMENT_FRACTION > aligned_length as f64 {
        region.is_match = 0;
        return region.is_match;
    }

    region.is_match = if error_rate <= final_error_rate {
        1
    } else if error_rate <= final_error_rate * 1.5 {
        3
    } else {
        0
    };

    if region.is_match == 1 {
        append_unmatched_windows(
            region.x_position_start as i32,
            region.x_position_end as i32,
            &mut region.window_cigar,
            block_length,
        );
    }
    region.is_match
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Golden values printed by the real `get_init_err_thres` and
    /// `double_error_threshold` with `block_s=375, block_err=15, rate=0.04`.
    #[test]
    fn thresholds_match_the_reference() {
        let cases = [
            (375i64, 15i64, 31i64),
            (100, 4, 8),
            (10, 1, 2),
            (3, 0, 0),
            (4, 1, 2),
            (300, 12, 31),
            (1, 0, 0),
        ];
        for (length, initial, doubled) in cases {
            assert_eq!(
                initial_error_threshold(length, 0.04, 375, 15),
                initial,
                "initial threshold at length {length}"
            );
            assert_eq!(
                double_error_threshold(initial, length),
                doubled,
                "doubled threshold at length {length}"
            );
        }

        // A full window with the ultra-long block size.
        assert_eq!(initial_error_threshold(500, 0.04, 500, 31), 31);
        assert_eq!(double_error_threshold(31, 500), 31);
    }

    #[test]
    fn similarity_pass_requires_enough_aligned_bases() {
        assert!(!similarity_pass(0, 100, false, 0.9, None));
        assert!(!similarity_pass(100, 0, false, 0.9, None));
        assert!(similarity_pass(100, 90, false, 0.9, None));
        assert!(!similarity_pass(100, 89, false, 0.9, None));
        assert!(similarity_pass(100, 99, false, 0.9, Some(0.02)));
        assert!(!similarity_pass(100, 97, false, 0.9, Some(0.02)));
    }

    #[test]
    fn unmatched_windows_fill_the_gaps_by_window_identifier() {
        let mut cigar = WindowCigar::default();
        // Three windows: [0,99], [100,199], [200,249]. Only the middle one was
        // verified, so it must end up in slot 1.
        let index = cigar.append_window(100, 199, 100, 199, 0, 0, 0, 15);
        assert_eq!(index, 0);

        append_unmatched_windows(0, 249, &mut cigar, 100);

        assert_eq!(cigar.windows.len(), 3);
        assert_eq!((cigar.windows[0].x_start, cigar.windows[0].x_end), (0, 99));
        assert_eq!(cigar.windows[0].y_end, -1, "sentinel target end");
        assert_eq!(cigar.windows[0].error, -1);
        assert_eq!((cigar.windows[1].x_start, cigar.windows[1].x_end), (100, 199));
        assert_eq!(cigar.windows[1].y_end, 199, "the verified window survived");
        assert_eq!((cigar.windows[2].x_start, cigar.windows[2].x_end), (200, 249));
        assert_eq!(cigar.windows[2].y_end, -1);
    }

    #[test]
    fn unmatched_windows_do_nothing_when_every_window_is_present() {
        let mut cigar = WindowCigar::default();
        for start in [0, 100, 200] {
            cigar.append_window(start, start + 99, start, start + 99, 0, 0, 0, 15);
        }
        append_unmatched_windows(0, 299, &mut cigar, 100);
        assert_eq!(cigar.windows.len(), 3);
        assert_eq!(cigar.windows[1].y_end, 199);
    }

    #[test]
    fn classification_uses_the_error_rate_and_expands_accepted_overlaps() {
        let mut region = OverlapRegion {
            x_position_start: 0,
            x_position_end: 249,
            alignment_length: 250,
            ..OverlapRegion::default()
        };
        region
            .window_cigar
            .append_window(100, 199, 100, 199, 0, 0, 0, 15);
        assert_eq!(classify_overlap(&mut region, 0.01, 0.03, 100), 1);
        assert_eq!(region.window_cigar.windows.len(), 3);
        assert_eq!(region.window_cigar.windows[0].y_end, -1);
        assert_eq!(region.window_cigar.windows[1].y_end, 199);
        assert_eq!(region.window_cigar.windows[2].y_end, -1);

        let mut weak = OverlapRegion {
            x_position_end: 249,
            alignment_length: 250,
            ..OverlapRegion::default()
        };
        assert_eq!(classify_overlap(&mut weak, 0.04, 0.03, 100), 3);
        assert!(weak.window_cigar.windows.is_empty());

        let mut rejected = OverlapRegion {
            x_position_end: 249,
            alignment_length: 100,
            ..OverlapRegion::default()
        };
        assert_eq!(classify_overlap(&mut rejected, 0.01, 0.03, 100), 0);
    }

    /// Golden values printed by the real `append_unmatched_wins` for the four
    /// ways verified windows can be scattered across a three-window overlap.
    #[test]
    fn unmatched_windows_match_the_reference_scatter() {
        // Each case lists the verified (start, end) pairs; the expected rows
        // are (x_start, x_end, y_end, error), where `y_end == -1` means the
        // window is a sentinel.
        let cases = [
            (
                "middle_and_last",
                vec![(100, 199), (200, 249)],
                vec![(0, 99, -1, -1), (100, 199, 200, 0), (200, 249, 201, 1)],
            ),
            (
                "only_middle",
                vec![(100, 199)],
                vec![(0, 99, -1, -1), (100, 199, 200, 0), (200, 249, -1, -1)],
            ),
            (
                "all_present",
                vec![(0, 99), (100, 199), (200, 249)],
                vec![(0, 99, 200, 0), (100, 199, 201, 1), (200, 249, 202, 2)],
            ),
            (
                "first_and_last",
                vec![(0, 99), (200, 249)],
                vec![(0, 99, 200, 0), (100, 199, -1, -1), (200, 249, 201, 1)],
            ),
        ];

        for (name, verified, expected) in cases {
            let mut cigar = WindowCigar::default();
            for (index, (start, end)) in verified.iter().enumerate() {
                cigar.append_window(
                    *start,
                    *end,
                    100 + index as i32,
                    200 + index as i32,
                    index as i16,
                    0,
                    0,
                    15,
                );
            }

            append_unmatched_windows(0, 249, &mut cigar, 100);

            let rows: Vec<(i32, i32, i32, i16)> = cigar
                .windows
                .iter()
                .map(|window| {
                    (
                        window.x_start,
                        window.x_end,
                        window.y_end,
                        window.error,
                    )
                })
                .collect();
            assert_eq!(rows, expected, "{name}");
        }
    }
}

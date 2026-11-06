//! `get_available_fully_covered_interval` and `check_if_fully_covered` from
//! `Correct.cpp`.
//!
//! After every overlap of a read has been classified, the read is walked once
//! more to decide whether the correction covered all of it. A window is covered
//! when enough full-length overlaps that passed the refinement reach it; a
//! window with no such overlap at all marks the read as abnormal, which the
//! caller uses to keep the uncorrected sequence.

use crate::correction::interval::WindowPool;
use crate::correction::region::MINIMUM_COVERAGE_THRESHOLD;
use crate::overlap::region::OverlapRegion;

/// Outcome of [`check_if_fully_covered`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CoverageSummary {
    /// Every window had enough usable overlaps, `return_flag`.
    pub fully_covered: bool,
    /// `0` when nothing suspicious was seen, `1` when a middle window had no
    /// usable overlap and `2` when only an end window did.
    pub abnormal: i32,
}

/// `get_available_fully_covered_interval`: counts the overlaps that touch the
/// window and, separately, those that cover all of it and are real matches.
///
/// `start_index` is the persistent scan cursor; it is updated exactly like the
/// reference so that each window only rescans the overlaps that can still
/// reach it. The returned flag is `1` when some overlap touched the window,
/// `0` when none did and `-2` when the cursor ran past the last overlap.
fn available_fully_covered_interval(
    window_start: i64,
    window_end: i64,
    overlaps: &[OverlapRegion],
    start_index: &mut u64,
    real_length: &mut i64,
    real_length_100: &mut i64,
) -> i32 {
    if window_start == 0 {
        *start_index = 0;
    }

    let mut index = *start_index as usize;
    while index < overlaps.len() {
        if window_end < overlaps[index].x_position_start as i64 {
            *start_index = 0;
            return 0;
        }
        *start_index = index as u64;
        break;
    }

    if index >= overlaps.len() {
        *start_index = overlaps.len() as u64;
        return -2;
    }

    let overlap_length = window_end - window_start + 1;
    *real_length = 0;
    let mut touched = 0i64;

    while index < overlaps.len() {
        let overlap = &overlaps[index];
        let covered = (overlap.x_position_end.min(window_end as u32) as i64)
            - (overlap.x_position_start.max(window_start as u32) as i64)
            + 1;

        if covered > 0 {
            touched += 1;
            if overlap_length == covered && overlap.is_match == 1 {
                *real_length += 1;
            }
            if overlap_length == covered && overlap.is_match == 100 {
                *real_length_100 += 1;
            }
        }

        if overlap.x_position_start as i64 > window_end {
            break;
        }
        index += 1;
    }

    if touched == 0 {
        0
    } else {
        1
    }
}

/// `check_if_fully_covered`: whether every window of the read has enough
/// verified overlaps, plus the abnormality flag.
///
/// `window_length` is `WINDOW` and `maximum_overlap_divergence` is
/// `asm_opt.max_ov_diff_ec`, which only sets how many windows the pool keeps as
/// context for the tail.
pub fn check_if_fully_covered(
    query_length: i64,
    overlaps: &[OverlapRegion],
    window_length: i64,
    maximum_overlap_divergence: f64,
) -> CoverageSummary {
    let mut pool = WindowPool::new(
        query_length,
        window_length,
        (1.0 / maximum_overlap_divergence) as i64,
    );

    let mut fully_covered = true;
    let mut abnormal = 0i32;
    let mut start_index = 0u64;
    // The reference keeps these across windows: the early `-2` and `0` returns
    // do not reset them.
    let mut real_length = 0i64;
    let mut real_length_100 = 0i64;
    let mut flag = 0i32;

    while flag != -2 {
        let Some((window_start, window_end)) = pool.next_window() else {
            break;
        };

        flag = available_fully_covered_interval(
            window_start,
            window_end,
            overlaps,
            &mut start_index,
            &mut real_length,
            &mut real_length_100,
        );

        if real_length < (MINIMUM_COVERAGE_THRESHOLD * 2) as i64 {
            fully_covered = false;
        }

        if real_length == 0 {
            if window_start != 0 && window_end != query_length - 1 {
                abnormal = 1;
            } else if abnormal == 0 {
                abnormal = 2;
            }
        }
    }

    CoverageSummary {
        fully_covered,
        abnormal,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn overlaps_of(ranges: &[(u32, u32)], copies: usize) -> Vec<OverlapRegion> {
        let mut overlaps = Vec::new();
        for (start, end) in ranges {
            for _ in 0..copies {
                overlaps.push(OverlapRegion {
                    x_position_start: *start,
                    x_position_end: *end,
                    is_match: 1,
                    ..OverlapRegion::default()
                });
            }
        }
        overlaps
    }

    /// Golden values printed by the real `check_if_fully_covered` with
    /// `max_ov_diff_ec = 0.04` and a 375 base window.
    #[test]
    fn matches_the_reference_coverage_decision() {
        let cases: &[(&str, &[(u32, u32)], usize, i64, bool, i32)] = &[
            ("six_whole", &[(0, 999)], 6, 1000, true, 0),
            ("five_whole", &[(0, 999)], 5, 1000, false, 0),
            ("none", &[], 0, 1000, false, 2),
            ("two_windows", &[(0, 749)], 6, 1000, false, 2),
            (
                "middle_gap",
                &[(0, 374), (750, 999)],
                6,
                1000,
                false,
                1,
            ),
            ("short_single_window", &[(0, 249)], 6, 250, true, 0),
        ];

        for (name, ranges, copies, length, covered, abnormal) in cases {
            let overlaps = overlaps_of(ranges, *copies);
            let summary = check_if_fully_covered(*length, &overlaps, 375, 0.04);
            assert_eq!(summary.fully_covered, *covered, "coverage of {name}");
            assert_eq!(summary.abnormal, *abnormal, "abnormality of {name}");
        }
    }

    #[test]
    fn a_weak_match_does_not_count_as_coverage() {
        // `is_match == 3` overlaps touch the window but are not counted.
        let mut overlaps = overlaps_of(&[(0, 999)], 6);
        for overlap in &mut overlaps {
            overlap.is_match = 3;
        }
        let summary = check_if_fully_covered(1000, &overlaps, 375, 0.04);
        assert!(!summary.fully_covered);
        // The first window sees no usable overlap, which only raises the flag
        // to 2; the second window is a middle window and raises it to 1.
        assert_eq!(summary.abnormal, 1);
    }
}

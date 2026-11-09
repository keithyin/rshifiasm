//! Correction constants and the window-to-target mapping, the counterpart of
//! the constants at the top of `Correct.h`, `determine_overlap_region` and
//! `fill_subregion` in `Correct.cpp`.

use crate::correction::target::TargetSource;

/// Correction window length, `WINDOW`.
pub const WINDOW_LENGTH: i64 = 375;

/// Edit budget of a window, `THRESHOLD`.
pub const THRESHOLD: i64 = 15;

/// Largest edit budget `determine_overlap_region` accepts, `THRESHOLD_MAX_SIZE`.
pub const THRESHOLD_MAXIMUM: i64 = 31;

/// Number of overlaps verified at once, `GROUP_SIZE`.
pub const GROUP_SIZE: usize = 4;

/// Fraction of the overlap length a window may differ by, `CORRECT_THRESHOLD`.
pub const CORRECTION_THRESHOLD: f64 = 0.60;

/// Same, for homopolymer runs, `CORRECT_THRESHOLD_HOMOPOLYMER`.
pub const HOMOPOLYMER_CORRECTION_THRESHOLD: f64 = 0.515;

/// Minimum number of overlapping reads before a window is corrected,
/// `MIN_COVERAGE_THRESHOLD`.
pub const MINIMUM_COVERAGE_THRESHOLD: i32 = 3;

/// Longest gap the corrector inserts or removes, `CORRECT_INDEL_LENGTH`.
pub const CORRECT_INDEL_LENGTH: i32 = 2;

/// `Adjust_Threshold`: a window that would allow no edit still gets one when it
/// is long enough to make a difference.
pub fn adjust_threshold(threshold: i64, length: i64) -> i64 {
    if threshold == 0 && length >= 4 {
        1
    } else {
        threshold
    }
}

/// The part of the target read a window is compared against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TargetRegion {
    /// `N` bases prepended because the window starts before the target.
    pub extra_begin: i32,
    /// `N` bases appended because the window ends past the target.
    pub extra_end: i32,
    /// First base of the target that is compared.
    pub start: i64,
    /// Number of target bases that are compared.
    pub length: i64,
}

/// `determine_overlap_region`: maps a window position on the query to the
/// region of the target it must be compared with.
///
/// Returns `None` when the window falls outside the target or too close to its
/// end to be worth aligning.
pub fn determine_overlap_region(
    threshold: i64,
    target_start: i64,
    target_length: i64,
    window_length: i64,
) -> Option<TargetRegion> {
    if target_start < 0
        || target_length <= target_start
        || target_length - target_start + 2 * threshold + THRESHOLD_MAXIMUM < window_length
    {
        return None;
    }

    let mut extra_begin = 0i32;
    let mut extended_start = target_start - threshold;
    let mut length = window_length.min(target_length - extended_start);
    let extra_end = window_length - length;

    if extended_start < 0 {
        extra_begin = -extended_start as i32;
        extended_start = 0;
        length -= extra_begin as i64;
    }

    Some(TargetRegion {
        extra_begin,
        extra_end: extra_end as i32,
        start: extended_start,
        length,
    })
}

/// `fill_subregion`: extracts a slice of the target and pads it with `N`.
///
/// `strand` selects the reverse complement, and `extra_begin`/`extra_end` are
/// the padding lengths `determine_overlap_region` asked for. This is the
/// non-HPC half of `fill_subregion`/`return_str_seq` in the reference: the
/// core bases come from the read store or, through
/// [`TargetSource::sub_region`], from `retrieve_u_seq`.
pub fn fill_subregion(
    targets: &TargetSource,
    read_identifier: usize,
    start: i64,
    length: i64,
    strand: u8,
    extra_begin: i32,
    extra_end: i32,
) -> Vec<u8> {
    let core = targets.sub_region(read_identifier, start, length, strand);
    let mut filled = vec![b'N'; extra_begin as usize];
    filled.extend_from_slice(&core);
    filled.resize(filled.len() + extra_end as usize, b'N');
    filled
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_short_window_still_gets_one_edit() {
        assert_eq!(adjust_threshold(0, 4), 1);
        assert_eq!(adjust_threshold(0, 3), 0);
        assert_eq!(adjust_threshold(2, 100), 2);
    }

    #[test]
    fn a_window_before_the_target_is_rejected() {
        assert_eq!(determine_overlap_region(-1, -5, 1000, 400), None);
    }

    #[test]
    fn a_window_past_the_target_is_rejected() {
        assert_eq!(determine_overlap_region(-1, 1000, 1000, 400), None);
        assert_eq!(determine_overlap_region(-1, 1200, 1000, 400), None);
    }

    #[test]
    fn a_window_too_close_to_the_target_end_is_rejected() {
        // 1000 - 990 + 2*(-1) + 31 = 39 < 400.
        assert_eq!(determine_overlap_region(-1, 990, 1000, 400), None);
    }

    #[test]
    fn a_window_inside_the_target_maps_directly() {
        let region = determine_overlap_region(15, 500, 2000, 405).unwrap();
        assert_eq!(region.extra_begin, 0);
        assert_eq!(region.extra_end, 0);
        // The window starts 15 bases earlier than the anchor.
        assert_eq!(region.start, 485);
        assert_eq!(region.length, 405);
    }

    #[test]
    fn a_window_near_the_target_start_is_padded() {
        let region = determine_overlap_region(15, 5, 2000, 405).unwrap();
        // start = 5 - 15 = -10, so ten Ns are prepended.
        assert_eq!(region.extra_begin, 10);
        assert_eq!(region.start, 0);
        // 405 window bases minus the ten that fell before the target.
        assert_eq!(region.length, 395);
        assert_eq!(region.extra_end, 0);
    }

    #[test]
    fn a_window_near_the_target_end_is_padded() {
        let region = determine_overlap_region(15, 1640, 2000, 405).unwrap();
        assert_eq!(region.extra_begin, 0);
        assert_eq!(region.start, 1625);
        // Only 375 target bases remain, so 30 bases are padded.
        assert_eq!(region.length, 375);
        assert_eq!(region.extra_end, 30);
    }

    /// Golden values printed by the real `determine_overlap_region`.
    #[test]
    fn matches_the_reference_mapping() {
        let cases = [
            // start, target length, window, threshold, extra begin, extra end, target start, length
            (5i64, 2000i64, 405i64, 15i64, 10i32, 0i32, 0i64, 395i64),
            (1640, 2000, 405, 15, 0, 30, 1625, 375),
            (1595, 2000, 405, 15, 0, 0, 1580, 405),
            (50, 2000, 405, 15, 0, 0, 35, 405),
            (0, 2000, 405, 15, 15, 0, 0, 390),
            (500, 2000, 405, 15, 0, 0, 485, 405),
        ];

        for (start, target_length, window, threshold, begin, end, mapped_start, mapped_length) in
            cases
        {
            let region = determine_overlap_region(threshold, start, target_length, window)
                .unwrap_or_else(|| panic!("start {start} must map"));
            assert_eq!(region.extra_begin, begin, "extra begin of start {start}");
            assert_eq!(region.extra_end, end, "extra end of start {start}");
            assert_eq!(region.start, mapped_start, "start of start {start}");
            assert_eq!(region.length, mapped_length, "length of start {start}");
        }

        // The three rejected inputs.
        assert_eq!(determine_overlap_region(-1, 990, 1000, 400), None);
        assert_eq!(determine_overlap_region(-1, -5, 1000, 400), None);
        assert_eq!(determine_overlap_region(-1, 1000, 1000, 400), None);
    }

    #[test]
    fn padding_keeps_the_window_length() {
        for (start, length, window) in [
            (5i64, 2000i64, 405i64),
            (1640, 2000, 405),
            (1595, 2000, 405),
            (50, 2000, 405),
            (0, 2000, 405),
        ] {
            let region = determine_overlap_region(15, start, length, window).unwrap();
            assert_eq!(
                region.extra_begin as i64 + region.length + region.extra_end as i64,
                window,
                "start {start}"
            );
        }
    }
}

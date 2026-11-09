//! Boundary correction of a window alignment, `fix_boundary` in `Correct.cpp`.
//!
//! When a traceback ends exactly at one edge of the padded target window, the
//! alignment may have been pushed against the padding rather than finding its
//! best position. The reference shifts the target window by up to `threshold`
//! bases and retries; the new alignment is accepted only when it has strictly
//! fewer errors.

use crate::alignment::banded_bpm::NO_ALIGNMENT_ERROR;
use crate::alignment::banded_path::{BandedPath, reserve_banded_bpm_path};
use crate::correction::region::{TargetRegion, determine_overlap_region, fill_subregion};
use crate::correction::target::TargetSource;
use crate::reads::read_store::ReadStore;

/// Inputs that do not change while fixing one window's boundary.
pub struct BoundaryContext<'a> {
    /// The query window, `x_string`.
    pub query: &'a [u8],
    /// Length of the query window.
    pub query_length: i64,
    /// Edit budget, `threshold`.
    pub threshold: i64,
    /// Padded target window length, `Window_Len`.
    pub target_window_length: i64,
    /// Original (global) start of the target window, `y_start`.
    pub total_target_start: i64,
    /// Traceback start inside the padded window, `real_y_start`.
    pub local_start: i64,
    /// Traceback end inside the padded window, `end_site`.
    pub local_end: i64,
    /// Padding of the original alignment.
    pub old_extra_begin: i32,
    pub old_extra_end: i32,
    /// Edit distance of the original alignment.
    pub old_error: u32,
    /// Target read and its orientation.
    pub reads: &'a ReadStore,
    pub target_identifier: usize,
    pub target_strand: u8,
}

/// Result of an accepted boundary correction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundaryFix {
    /// New global target start.
    pub total_target_start: i64,
    /// New traceback start inside the padded window.
    pub local_start: i32,
    /// New traceback end inside the padded window.
    pub local_end: i32,
    pub extra_begin: i32,
    pub extra_end: i32,
    pub error: u32,
    /// Padded target substring of the accepted alignment, needed by
    /// `generate_cigar` as its `y` argument.
    pub target_sequence: Vec<u8>,
    /// Traceback path of the accepted alignment.
    pub path: Vec<u8>,
}

/// `fix_boundary`: tries a shifted target window and returns the improved
/// alignment, or `None` when the original one is already the best available.
pub fn fix_boundary(context: &BoundaryContext<'_>) -> Option<BoundaryFix> {
    let target_length = context.reads.read_length(context.target_identifier) as i64;

    let (shifted_start, comparison_start) = if context.local_start == 0 {
        if context.old_extra_begin != 0 {
            return None;
        }
        (context.total_target_start, context.total_target_start)
    } else if context.local_end == context.target_window_length - 1 {
        if context.old_extra_end != 0 {
            return None;
        }
        let total_end = context.total_target_start + context.local_end;
        (
            total_end - context.query_length + 1,
            total_end - context.local_end,
        )
    } else {
        return None;
    };

    let TargetRegion {
        extra_begin,
        extra_end,
        start,
        length,
    } = determine_overlap_region(
        context.threshold,
        shifted_start,
        target_length,
        context.target_window_length,
    )?;

    // A shift that leaves the window where it was cannot improve anything.
    if start == comparison_start {
        return None;
    }

    let target_sequence = fill_subregion(
        &TargetSource::Reads(context.reads),
        context.target_identifier,
        start,
        length,
        context.target_strand,
        extra_begin,
        extra_end,
    );

    let BandedPath {
        end_site,
        start_site,
        error,
        path,
    } = reserve_banded_bpm_path(
        &target_sequence,
        context.target_window_length as usize,
        context.query,
        context.query_length as usize,
        context.threshold as usize,
        -1,
        -1,
    )?;

    if error == NO_ALIGNMENT_ERROR || error >= context.old_error {
        return None;
    }

    Some(BoundaryFix {
        total_target_start: start,
        local_start: start_site,
        local_end: end_site,
        extra_begin,
        extra_end,
        error,
        target_sequence,
        path,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reads() -> ReadStore {
        let mut store = ReadStore::default();
        store.push_read("query", &vec![b'A'; 2000]);
        store.push_read("target", &vec![b'A'; 2000]);
        store
    }

    #[test]
    fn an_alignment_away_from_both_boundaries_is_left_alone() {
        let store = reads();
        let context = BoundaryContext {
            query: &vec![b'A'; 10],
            query_length: 10,
            threshold: 2,
            target_window_length: 14,
            total_target_start: 500,
            local_start: 3,
            local_end: 10,
            old_extra_begin: 0,
            old_extra_end: 0,
            old_error: 2,
            reads: &store,
            target_identifier: 1,
            target_strand: 0,
        };
        assert!(fix_boundary(&context).is_none());
    }

    #[test]
    fn padding_at_a_boundary_prevents_a_shift() {
        let store = reads();
        let mut context = BoundaryContext {
            query: &vec![b'A'; 10],
            query_length: 10,
            threshold: 2,
            target_window_length: 14,
            total_target_start: 500,
            local_start: 0,
            local_end: 9,
            old_extra_begin: 5,
            old_extra_end: 0,
            old_error: 2,
            reads: &store,
            target_identifier: 1,
            target_strand: 0,
        };
        assert!(fix_boundary(&context).is_none());

        context.old_extra_begin = 0;
        context.local_start = 3;
        context.local_end = 13;
        context.old_extra_end = 5;
        assert!(fix_boundary(&context).is_none());
    }

    #[test]
    fn a_better_shifted_alignment_is_accepted() {
        let store = reads();
        // Start at the left edge of the target: the shifted window is clamped
        // to zero, so determine_overlap_region returns a different start.
        let context = BoundaryContext {
            query: &vec![b'A'; 10],
            query_length: 10,
            threshold: 2,
            target_window_length: 14,
            total_target_start: 5,
            local_start: 0,
            local_end: 9,
            old_extra_begin: 0,
            old_extra_end: 0,
            old_error: 9,
            reads: &store,
            target_identifier: 1,
            target_strand: 0,
        };

        let fix = fix_boundary(&context).expect("the padded window should align perfectly");
        // Shifting by `threshold` moves the window start from 5 to 3.
        assert_eq!(fix.total_target_start, 3);
        assert_eq!(fix.extra_begin, 0);
        assert_eq!(fix.error, 0);
    }

    /// Golden values printed by the real `fix_boundary` for a 2000 base target
    /// of all `A`s, a 10 base query window and threshold 2.
    #[test]
    fn matches_the_reference_boundary_correction() {
        let store = reads();
        struct Case {
            name: &'static str,
            total_start: i64,
            local_start: i64,
            local_end: i64,
            extra_begin: i32,
            extra_end: i32,
            old_error: u32,
            accepted: bool,
            new_total_start: i64,
            new_local_start: i32,
            new_local_end: i32,
            new_error: u32,
        }
        let cases = [
            Case { name: "shift_left", total_start: 5, local_start: 0, local_end: 9,
                   extra_begin: 0, extra_end: 0, old_error: 9, accepted: true,
                   new_total_start: 3, new_local_start: 2, new_local_end: 11, new_error: 0 },
            Case { name: "shift_left_no_gain", total_start: 5, local_start: 0, local_end: 9,
                   extra_begin: 0, extra_end: 0, old_error: 0, accepted: false,
                   new_total_start: 0, new_local_start: 0, new_local_end: 0, new_error: 0 },
            Case { name: "left_edge_no_shift", total_start: 0, local_start: 0, local_end: 9,
                   extra_begin: 0, extra_end: 0, old_error: 9, accepted: false,
                   new_total_start: 0, new_local_start: 0, new_local_end: 0, new_error: 0 },
            Case { name: "with_extra_begin", total_start: 5, local_start: 0, local_end: 9,
                   extra_begin: 3, extra_end: 0, old_error: 9, accepted: false,
                   new_total_start: 0, new_local_start: 0, new_local_end: 0, new_error: 0 },
            Case { name: "away_from_boundary", total_start: 500, local_start: 3, local_end: 10,
                   extra_begin: 0, extra_end: 0, old_error: 9, accepted: false,
                   new_total_start: 0, new_local_start: 0, new_local_end: 0, new_error: 0 },
            Case { name: "right_boundary", total_start: 500, local_start: 3, local_end: 13,
                   extra_begin: 0, extra_end: 0, old_error: 9, accepted: true,
                   new_total_start: 502, new_local_start: 2, new_local_end: 11, new_error: 0 },
        ];

        let query = vec![b'A'; 10];
        for case in cases {
            let context = BoundaryContext {
                query: &query,
                query_length: 10,
                threshold: 2,
                target_window_length: 14,
                total_target_start: case.total_start,
                local_start: case.local_start,
                local_end: case.local_end,
                old_extra_begin: case.extra_begin,
                old_extra_end: case.extra_end,
                old_error: case.old_error,
                reads: &store,
                target_identifier: 1,
                target_strand: 0,
            };

            match fix_boundary(&context) {
                None => assert!(!case.accepted, "{} must be accepted", case.name),
                Some(fix) => {
                    assert!(case.accepted, "{} must be rejected", case.name);
                    assert_eq!(fix.total_target_start, case.new_total_start, "total of {}", case.name);
                    assert_eq!(fix.local_start, case.new_local_start, "start of {}", case.name);
                    assert_eq!(fix.local_end, case.new_local_end, "end of {}", case.name);
                    assert_eq!(fix.error, case.new_error, "error of {}", case.name);
                    assert_eq!(fix.path.len(), 10, "path length of {}", case.name);
                    assert_eq!(fix.extra_begin, 0);
                    assert_eq!(fix.extra_end, 0);
                }
            }
        }
    }
}

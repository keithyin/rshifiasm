//! `correct_overlap`, the per-read correction driver of `Correct.cpp`.
//!
//! One read is corrected in four steps:
//!
//! 1. every window is verified against the overlaps that reach it
//!    (`verify_window`),
//! 2. the windows that no overlap reached are retried and every window gets an
//!    exact cigar, then each overlap is classified (`recalcate_window_advance`),
//! 3. the corrected sequence is built window by window from the classified
//!    overlaps (`generate_consensus`),
//! 4. the read is walked once more to report whether correction covered it
//!    (`check_if_fully_covered`).
//!
//! `partition_overlaps_advance`, the haplotype evidence pass that the reference
//! runs between steps 2 and 3, only feeds phasing and is not needed for the
//! corrected sequence, so it is not called here yet.

use crate::correction::advance::recalcate_window_advance;
use crate::correction::coverage::{CoverageSummary, check_if_fully_covered};
use crate::correction::interval::{
    INTERVAL_PAST_THE_END, IntervalBuffer, WindowPool, get_interval,
};
use crate::correction::region::WINDOW_LENGTH;
use crate::correction::verify::{VerificationContext, verify_window};
use crate::correction::window_consensus::generate_consensus;
use crate::overlap::region::OverlapRegion;
use crate::reads::read_store::ReadStore;

/// Result of correcting one read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CorrectionOutcome {
    /// Corrected sequence of the read.
    pub sequence: Vec<u8>,
    /// Number of bases the correction changed.
    pub corrected_base: u64,
    /// Whether enough overlaps covered every window.
    pub fully_covered: bool,
    /// `0` when nothing suspicious was seen, see [`CoverageSummary`].
    pub abnormal: i32,
}

impl CorrectionOutcome {
    /// The summary of a read that was not corrected at all.
    pub fn unchanged(sequence: &[u8]) -> Self {
        Self {
            sequence: sequence.to_vec(),
            corrected_base: 0,
            fully_covered: false,
            abnormal: 0,
        }
    }
}

/// `correct_overlap`: corrects one read against its overlaps.
///
/// `overlaps` must be sorted by query start, as `ha_get_candidates_interface`
/// leaves them. `maximum_overlap_divergence` is `asm_opt.max_ov_diff_ec` and
/// `final_error_rate` is `asm_opt.max_ov_diff_final`.
pub fn correct_overlap(
    overlaps: &mut [OverlapRegion],
    reads: &ReadStore,
    query_sequence: &[u8],
    maximum_overlap_divergence: f64,
    final_error_rate: f64,
) -> CorrectionOutcome {
    // `ha_get_new_candidates` rebuilds the overlap list from scratch, so each
    // overlap starts without windows.
    for region in overlaps.iter_mut() {
        region.window_cigar.clear();
        region.boundary_cigars.clear();
        region.alignment_length = 0;
        region.is_match = 0;
    }

    verify_all_windows(
        overlaps,
        reads,
        query_sequence,
        maximum_overlap_divergence,
    );

    recalcate_window_advance(
        overlaps,
        reads,
        query_sequence,
        WINDOW_LENGTH,
        maximum_overlap_divergence,
        final_error_rate,
    );

    // `generate_consensus` is only called for the correction rounds that
    // produce a sequence; the reference guards it with `is_consensus`.
    let consensus = generate_consensus(
        query_sequence,
        overlaps,
        reads,
        maximum_overlap_divergence,
        &mut crate::correction::poa::PoaGraph::new(),
        &mut crate::correction::dag_con::DagCon::default(),
        None,
    );

    let coverage: CoverageSummary =
        check_if_fully_covered(query_sequence.len() as i64, overlaps, WINDOW_LENGTH, maximum_overlap_divergence);

    CorrectionOutcome {
        sequence: consensus.sequence,
        corrected_base: consensus.corrected_base,
        fully_covered: coverage.fully_covered,
        abnormal: coverage.abnormal,
    }
}

/// The `get_Window` / `get_interval` / `verify_window` loop of
/// `correct_overlap`.
fn verify_all_windows(
    overlaps: &mut [OverlapRegion],
    reads: &ReadStore,
    query_sequence: &[u8],
    maximum_overlap_divergence: f64,
) {
    if overlaps.is_empty() {
        return;
    }

    let mut pool = WindowPool::new(
        query_sequence.len() as i64,
        WINDOW_LENGTH,
        (1.0 / maximum_overlap_divergence) as i64,
    );
    let mut buffer = IntervalBuffer::with_capacity(overlaps.len());
    let context = VerificationContext {
        reads,
        query_sequence,
        maximum_overlap_divergence,
    };

    let mut flag = 0i32;
    while flag != INTERVAL_PAST_THE_END {
        let Some((window_start, window_end)) = pool.next_window() else {
            break;
        };
        buffer.length = 0;
        buffer.partial_length = 0;
        flag = get_interval(window_start, window_end, overlaps, &mut buffer, pool.window_length);
        verify_window(window_start, window_end, overlaps, &buffer, &context);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::correction::region::THRESHOLD;
    use crate::overlap::region::FakeCigar;

    fn reads() -> ReadStore {
        let mut reads = ReadStore::default();
        reads.push_read("query", &vec![b'A'; 1000]);
        reads.push_read("target", &vec![b'A'; 1000]);
        reads
    }

    fn full_overlap() -> OverlapRegion {
        let mut fake_cigar = FakeCigar::default();
        fake_cigar.push(0, 0);
        fake_cigar.push(999, 0);
        OverlapRegion {
            x_identifier: 0,
            x_position_start: 0,
            x_position_end: 999,
            y_identifier: 1,
            y_position_start: 0,
            y_position_end: 999,
            fake_cigar,
            ..OverlapRegion::default()
        }
    }

    #[test]
    fn a_read_without_overlaps_is_copied_verbatim() {
        let reads = reads();
        let query = b"ACGTACGTACGT".to_vec();
        let outcome = correct_overlap(&mut [], &reads, &query, 0.04, 0.05);

        assert_eq!(outcome.sequence, query);
        assert_eq!(outcome.corrected_base, 0);
        assert!(!outcome.fully_covered);
        assert_eq!(outcome.abnormal, 2);
    }

    #[test]
    fn a_matching_overlap_keeps_the_query_unchanged() {
        let reads = reads();
        let query = vec![b'A'; 1000];
        let mut overlaps = vec![full_overlap()];

        let outcome = correct_overlap(&mut overlaps, &reads, &query, 0.04, 0.05);

        assert_eq!(outcome.sequence, query);
        assert_eq!(outcome.corrected_base, 0);
        // One overlap is below the coverage the reference requires per window,
        // but the window it does reach is aligned.
        assert!(!outcome.fully_covered);
        assert_eq!(outcome.abnormal, 0);
        // The overlap was verified, classified and given one cigar per window.
        assert_eq!(overlaps[0].is_match, 1);
        assert_eq!(overlaps[0].window_cigar.windows.len(), 3);
        assert_eq!(overlaps[0].alignment_length, 1000);
        assert_eq!(overlaps[0].window_cigar.windows[0].error_threshold, THRESHOLD as i16);
    }

    #[test]
    fn six_matching_overlaps_cover_the_read() {
        let reads = reads();
        let query = vec![b'A'; 1000];
        let mut overlaps: Vec<OverlapRegion> = (0..6).map(|_| full_overlap()).collect();

        let outcome = correct_overlap(&mut overlaps, &reads, &query, 0.04, 0.05);

        assert_eq!(outcome.sequence, query);
        assert_eq!(outcome.corrected_base, 0);
        assert!(outcome.fully_covered, "six overlaps cover every window");
        assert_eq!(outcome.abnormal, 0);
        for region in &overlaps {
            assert_eq!(region.is_match, 1);
            assert_eq!(region.window_cigar.windows.len(), 3);
        }
    }

    #[test]
    fn a_single_overlap_is_not_enough_coverage() {
        let reads = reads();
        let query = vec![b'A'; 1000];
        let mut overlaps = vec![full_overlap(), full_overlap()];

        let outcome = correct_overlap(&mut overlaps, &reads, &query, 0.04, 0.05);

        // Two overlaps are below the six the reference requires per window.
        assert!(!outcome.fully_covered);
        assert_eq!(outcome.abnormal, 0);
        assert_eq!(outcome.sequence, query);
    }
}

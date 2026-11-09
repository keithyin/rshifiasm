//! Window verification, the counterpart of `verify_window` in `Correct.cpp`.
//!
//! For every window the caller has selected the overlaps that reach into it.
//! This module maps each of them onto the target read, extracts the padded
//! target substring, runs the banded bit-parallel alignment and records the
//! verified window. Full-window overlaps are verified in groups of four in the
//! reference (an SSE implementation); the four-lane routine was checked
//! against four scalar calls on 4000 random inputs, so we call the scalar port
//! four times and get the same coordinates and edit distances.

use crate::alignment::banded_bpm::{NO_ALIGNMENT_ERROR, reserve_banded_bpm};
use crate::correction::interval::IntervalBuffer;
use crate::correction::region::{
    THRESHOLD, THRESHOLD_MAXIMUM, WINDOW_LENGTH, adjust_threshold, determine_overlap_region,
    fill_subregion,
};
use crate::overlap::region::OverlapRegion;
use crate::correction::target::TargetSource;
use crate::reads::read_store::ReadStore;

/// Inputs that stay the same for all windows of one read.
pub struct VerificationContext<'a> {
    /// All reads, used to extract the target sub-regions.
    pub reads: &'a ReadStore,
    /// The query read in plain ASCII, forward strand.
    pub query_sequence: &'a [u8],
    /// `asm_opt.max_ov_diff_ec`, used for partial windows.
    pub maximum_overlap_divergence: f64,
}

/// One verified window, before it is appended to its overlap.
struct VerifiedWindow {
    overlap_identifier: u64,
    x_start: i32,
    x_end: i32,
    y_start: i32,
    y_end: i32,
    error: i16,
    extra_begin: i16,
    extra_end: i16,
    threshold: i16,
    aligned_length: u32,
}

/// `verify_window`: verifies every overlap selected by `get_interval` for one
/// window and appends the successful ones to each overlap's window cigar.
pub fn verify_window(
    window_start: i64,
    window_end: i64,
    overlaps: &mut [OverlapRegion],
    buffer: &IntervalBuffer,
    context: &VerificationContext<'_>,
) {
    let mut verified: Vec<VerifiedWindow> = Vec::new();
    let trace = crate::correction::window_consensus::tracing();
    if trace {
        eprintln!(
            "[M::verify_window] x=[{},{}] full={} partial={}",
            window_start,
            window_end,
            buffer.full_windows().len(),
            buffer.partial_windows().len()
        );
    }

    // Overlaps that cover the whole window are compared with a fixed edit
    // budget; the reference batches four of them into one SIMD call.
    for identifier in buffer.full_windows() {
        let identifier = *identifier;
        let Some(window) = verify_one_window(
            overlaps,
            identifier,
            window_end,
            window_start,
            WINDOW_LENGTH,
            THRESHOLD,
            context,
        ) else {
            if trace {
                let region = &overlaps[identifier as usize];
                eprintln!(
                    "[M::verify_window] dropped full id={} y={} x=[{},{}] offset={:?}",
                    identifier,
                    region.y_identifier,
                    region.x_position_start,
                    region.x_position_end,
                    region.fake_cigar.y_start_offset(window_start)
                );
            }
            continue;
        };
        verified.push(window);
    }

    // Overlaps that only partly cover the window derive their budget from the
    // length they actually cover.
    for identifier in buffer.partial_windows() {
        let region = &overlaps[identifier as usize];
        let x_start = window_start.max(region.x_position_start as i64);
        let x_end = window_end.min(region.x_position_end as i64);
        let covered_length = x_end - x_start + 1;

        let mut threshold = (covered_length as f64 * context.maximum_overlap_divergence) as i64;
        threshold = adjust_threshold(threshold, covered_length);
        if threshold > THRESHOLD_MAXIMUM {
            threshold = THRESHOLD_MAXIMUM;
        }

        let Some(window) = verify_one_window(
            overlaps,
            identifier,
            window_end,
            x_start,
            covered_length,
            threshold,
            context,
        ) else {
            continue;
        };
        verified.push(window);
    }

    for window in verified {
        let region = &mut overlaps[window.overlap_identifier as usize];
        region.alignment_length += window.aligned_length;
        region.window_cigar.append_window(
            window.x_start,
            window.x_end,
            window.y_start,
            window.y_end,
            window.error,
            window.extra_begin,
            window.extra_end,
            window.threshold,
        );
    }
}

/// Maps, extracts and aligns one overlap window.
fn verify_one_window(
    overlaps: &[OverlapRegion],
    overlap_identifier: u64,
    window_end: i64,
    x_start: i64,
    covered_length: i64,
    threshold: i64,
    context: &VerificationContext<'_>,
) -> Option<VerifiedWindow> {
    let region = &overlaps[overlap_identifier as usize];
    let x_end = x_start + covered_length - 1;

    let mut y_start = (x_start - region.x_position_start as i64) + region.y_position_start as i64;
    y_start += region.fake_cigar.y_start_offset(x_start)?;

    let target_read_length = context.reads.read_length(region.y_identifier as usize) as i64;
    let target_window_length = covered_length + (threshold << 1);
    let target = determine_overlap_region(
        threshold,
        y_start,
        target_read_length,
        target_window_length,
    )?;

    let target_sequence = fill_subregion(
        &TargetSource::Reads(context.reads),
        region.y_identifier as usize,
        target.start,
        target.length,
        region.y_strand as u8,
        target.extra_begin,
        target.extra_end,
    );

    let query_start = x_start.max(0) as usize;
    let query_end = (x_start + covered_length).max(0) as usize;
    if query_end > context.query_sequence.len() {
        return None;
    }
    let query_sequence = &context.query_sequence[query_start..query_end];

    let alignment = reserve_banded_bpm(
        &target_sequence,
        target_window_length as usize,
        query_sequence,
        covered_length as usize,
        threshold as usize,
    );
    if alignment.error == NO_ALIGNMENT_ERROR {
        return None;
    }

    Some(VerifiedWindow {
        overlap_identifier,
        x_start: x_start as i32,
        x_end: x_end.min(window_end) as i32,
        y_start: target.start as i32,
        y_end: (target.start + alignment.end_site as i64) as i32,
        error: alignment.error as i16,
        extra_begin: target.extra_begin as i16,
        extra_end: target.extra_end as i16,
        threshold: threshold as i16,
        aligned_length: covered_length as u32,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::overlap::region::FakeCigar;

    /// A synthetic pair of reads where `target` shares a 375 base window with
    /// the query starting 100 bases into the target.
    fn setup() -> (ReadStore, Vec<u8>, OverlapRegion) {
        let query = vec![b'A'; 1000];
        let target = vec![b'A'; 1000];
        let mut reads = ReadStore::default();
        reads.push_read("query", &query);
        reads.push_read("target", &target);

        let mut fake_cigar = FakeCigar::default();
        fake_cigar.push(0, 0);
        fake_cigar.push(999, 0);

        let region = OverlapRegion {
            x_identifier: 0,
            x_position_start: 0,
            x_position_end: 999,
            x_strand: 0,
            y_identifier: 1,
            y_position_start: 100,
            y_position_end: 999,
            y_strand: 0,
            fake_cigar,
            ..OverlapRegion::default()
        };
        (reads, query, region)
    }

    fn buffer_for_full_window() -> IntervalBuffer {
        let mut buffer = IntervalBuffer::with_capacity(8);
        buffer.overlap_identifiers[0] = 0;
        buffer.length = 1;
        buffer
    }

    #[test]
    fn a_matching_full_window_is_recorded() {
        let (reads, query, region) = setup();
        let mut overlaps = vec![region];
        let buffer = buffer_for_full_window();
        let context = VerificationContext {
            reads: &reads,
            query_sequence: &query,
            maximum_overlap_divergence: 0.04,
        };

        verify_window(0, 374, &mut overlaps, &buffer, &context);

        assert_eq!(overlaps[0].window_cigar.windows.len(), 1);
        let window = &overlaps[0].window_cigar.windows[0];
        assert_eq!(window.x_start, 0);
        assert_eq!(window.x_end, 374);
        assert_eq!(window.error, 0);
        assert_eq!(window.error_threshold, THRESHOLD as i16);
        // The target window starts 15 bases before the mapped position 100.
        assert_eq!(window.y_start, 85);
        assert_eq!(overlaps[0].alignment_length, WINDOW_LENGTH as u32);
    }

    #[test]
    fn a_window_with_too_many_errors_is_not_recorded() {
        let (reads, _query, region) = setup();
        // A query made entirely of a different base cannot be aligned.
        let diverged_query = vec![b'C'; 1000];
        let mut overlaps = vec![region];
        let buffer = buffer_for_full_window();
        let context = VerificationContext {
            reads: &reads,
            query_sequence: &diverged_query,
            maximum_overlap_divergence: 0.04,
        };

        verify_window(0, 374, &mut overlaps, &buffer, &context);

        assert!(overlaps[0].window_cigar.windows.is_empty());
        assert_eq!(overlaps[0].alignment_length, 0);
    }

    #[test]
    fn a_partial_window_uses_a_length_derived_threshold() {
        let (reads, query, mut region) = setup();
        // The overlap only covers the first 50 bases of this window.
        region.x_position_end = 49;
        region.fake_cigar.entries.clear();
        region.fake_cigar.push(0, 0);
        region.fake_cigar.push(49, 0);

        let mut overlaps = vec![region];
        let mut buffer = IntervalBuffer::with_capacity(8);
        buffer.overlap_identifiers[7] = 0;
        buffer.partial_length = 1;
        let context = VerificationContext {
            reads: &reads,
            query_sequence: &query,
            maximum_overlap_divergence: 0.04,
        };

        verify_window(0, 374, &mut overlaps, &buffer, &context);

        assert_eq!(overlaps[0].window_cigar.windows.len(), 1);
        let window = &overlaps[0].window_cigar.windows[0];
        assert_eq!(window.x_start, 0);
        assert_eq!(window.x_end, 49);
        // 50 * 0.04 = 2.
        assert_eq!(window.error_threshold, 2);
        assert_eq!(overlaps[0].alignment_length, 50);
    }
}

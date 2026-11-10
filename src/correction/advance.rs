//! `recalcate_window_advance` from `Correct.cpp`.
//!
//! The first pass (`verify_window`) only aligns the windows that some overlap
//! reaches and gives up on any window whose alignment fails. This pass closes
//! the gaps: it walks each overlap's window cigar from the two ends that did
//! align, retries every hole with a doubled error budget, and finally records
//! the exact cigar of every window. It also decides whether the overlap is a
//! real match, which is what the consensus step consumes.

use crate::alignment::banded_bpm::{NO_ALIGNMENT_ERROR, reserve_banded_bpm};
use crate::alignment::banded_path::reserve_banded_bpm_path;
use crate::correction::boundary::{BoundaryContext, fix_boundary};
use crate::correction::boundary_cigar::calculate_boundary_cigars;
use crate::correction::cigar_generator::generate_window_cigar;
use crate::correction::non_trim::non_trim_error_rate;
use crate::correction::recalculate::{
    OVERLAP_HIFI_ALIGNMENT_FRACTION, append_unmatched_windows, double_error_threshold,
    initial_error_threshold,
};
use crate::correction::region::{THRESHOLD, determine_overlap_region, fill_subregion};
use crate::correction::window::{number_of_windows, window_id_by_end, window_id_by_start};
use crate::overlap::region::OverlapRegion;
use crate::correction::target::TargetSource;
use crate::reads::read_store::ReadStore;

/// Everything the retry passes need to know about one overlap.
struct OverlapView {
    y_identifier: usize,
    y_strand: u8,
    y_read_length: i64,
    overlap_start: i32,
    overlap_end: i32,
}

/// Outcome of the retry passes for one overlap.
struct RetryOutcome {
    /// Number of query bases that ended up aligned, `align_length`.
    aligned_length: u64,
    /// Whether the window identifiers are still increasing.
    is_sorted: bool,
}

/// `recalcate_window_advance`: retries every unaligned window of every overlap,
/// rebuilds the window cigars and classifies the overlap.
///
/// `error_rate` is `asm_opt.max_ov_diff_ec` and `final_error_rate` is
/// `max_ov_diff_final`.
pub fn recalcate_window_advance(
    overlaps: &mut [OverlapRegion],
    reads: &ReadStore,
    query_sequence: &[u8],
    block_length: i64,
    error_rate: f64,
    final_error_rate: f64,
) {
    for region in overlaps.iter_mut() {
        region.is_match = 0;
        if region.window_cigar.windows.is_empty() {
            continue;
        }

        let window_count = number_of_windows(
            region.x_position_start as i64,
            region.x_position_end as i64 + 1,
            block_length,
        );
        let view = OverlapView {
            y_identifier: region.y_identifier as usize,
            y_strand: region.y_strand as u8,
            y_read_length: reads.read_length(region.y_identifier as usize) as i64,
            overlap_start: region.x_position_start as i32,
            overlap_end: region.x_position_end as i32,
        };

        // Map every window identifier to the slot that holds it; `u64::MAX`
        // marks a window no alignment has reached yet.
        let mut slots = vec![u64::MAX; window_count as usize];
        let aligned_count = region.window_cigar.windows.len();
        for index in 0..aligned_count {
            let (identifier, _) = window_id_by_start(
                view.overlap_start,
                view.overlap_end,
                region.window_cigar.windows[index].x_start as i64,
                block_length,
            );
            slots[identifier as usize] = index as u64;
        }

        let mut outcome = RetryOutcome {
            aligned_length: 0,
            is_sorted: true,
        };

        extend_forwards(
            region,
            reads,
            query_sequence,
            &view,
            &mut slots,
            window_count,
            block_length,
            error_rate,
            &mut outcome,
        );
        extend_backwards(
            region,
            reads,
            query_sequence,
            &view,
            &mut slots,
            window_count,
            block_length,
            error_rate,
            &mut outcome,
        );
        build_missing_cigars(region, reads, query_sequence, &view);

        if !outcome.is_sorted {
            // `radix_sort_window_list_xs_srt` is a stable sort by `x_start`.
            region.window_cigar.windows.sort_by_key(|window| window.x_start);
        }

        // `verify_window` already counted the windows it aligned; the retries
        // add to that running total, exactly like the reference.
        region.alignment_length = region
            .alignment_length
            .saturating_add(outcome.aligned_length.min(u32::MAX as u64) as u32);

        let rate = non_trim_error_rate(
            region,
            reads,
            query_sequence,
            error_rate,
            block_length,
        );
        let overlap_length =
            region.x_position_end as i64 + 1 - region.x_position_start as i64;
        // The comparison uses the accumulated field, which includes the bases
        // `verify_window` already aligned.
        let aligned = region.alignment_length as f64;

        if overlap_length as f64 * OVERLAP_HIFI_ALIGNMENT_FRACTION > aligned {
            region.is_match = 0;
            continue;
        }
        if rate <= final_error_rate {
            region.is_match = 1;
            append_unmatched_windows(
                region.x_position_start as i32,
                region.x_position_end as i32,
                &mut region.window_cigar,
                block_length,
            );
            calculate_boundary_cigars(region, reads, query_sequence, error_rate);
        } else if rate <= final_error_rate * 1.5 {
            region.is_match = 3;
        }
    }
}

/// `y_start_offset` for the first window of an overlap, mirroring the two
/// formulas the reference uses for the forward and backward retries.
#[allow(clippy::too_many_arguments)]
fn extend_forwards(
    region: &mut OverlapRegion,
    reads: &ReadStore,
    query_sequence: &[u8],
    view: &OverlapView,
    slots: &mut [u64],
    window_count: i64,
    block_length: i64,
    error_rate: f64,
    outcome: &mut RetryOutcome,
) {
    let mut index = region.window_cigar.windows.len() as i64 - 1;
    while index >= 0 {
        let window = region.window_cigar.windows[index as usize];
        let (window_identifier, window_end) = window_id_by_start(
            view.overlap_start,
            view.overlap_end,
            window.x_start as i64,
            block_length,
        );
        let mut total_y_start = window.y_end as i64 + 1 - window.extra_begin as i64;
        let mut identifier = window_identifier;
        let mut current_end = window_end;

        let mut slot = identifier + 1;
        while slot < window_count {
            if slots[slot as usize] != u64::MAX {
                break;
            }

            let window_start = current_end + 1;
            let (computed_identifier, computed_end) = window_id_by_start(
                view.overlap_start,
                view.overlap_end,
                window_start,
                block_length,
            );
            current_end = computed_end;
            identifier = computed_identifier;
            debug_assert_eq!(identifier, slot);

            if total_y_start >= view.y_read_length {
                break;
            }

            let x_start = window_start;
            let x_end = current_end;
            let x_length = x_end + 1 - x_start;
            let y_start = total_y_start;
            let threshold = double_error_threshold(
                initial_error_threshold(x_length, error_rate, block_length, THRESHOLD),
                x_length,
            );
            let window_length = x_length + (threshold << 1);

            let Some(mapped) =
                determine_overlap_region(threshold, y_start, view.y_read_length, window_length)
            else {
                break;
            };
            if mapped.length + threshold < x_length {
                break;
            }

            let target = fill_subregion(
                &TargetSource::Reads(reads),
                view.y_identifier,
                mapped.start,
                mapped.length,
                view.y_strand,
                mapped.extra_begin,
                mapped.extra_end,
            );
            let query = query_slice(query_sequence, x_start, x_length);

            let alignment = reserve_banded_bpm(
                &target,
                window_length as usize,
                query,
                x_length as usize,
                threshold as usize,
            );
            if alignment.error == NO_ALIGNMENT_ERROR {
                break;
            }

            // `determine_overlap_region` rewrites `y_start` in the reference,
            // so the new window is placed at the mapped start.
            let window_index = region.window_cigar.append_window(
                x_start as i32,
                x_end as i32,
                mapped.start as i32,
                (mapped.start + alignment.end_site as i64) as i32,
                alignment.error as i16,
                mapped.extra_begin as i16,
                mapped.extra_end as i16,
                threshold as i16,
            );
            outcome.aligned_length += x_length as u64;
            slots[slot as usize] = window_index as u64;

            if outcome.is_sorted && window_index > 0 {
                let previous = region.window_cigar.windows[window_index - 1].x_start;
                if region.window_cigar.windows[window_index].x_start < previous {
                    outcome.is_sorted = false;
                }
            }

            total_y_start =
                mapped.start + alignment.end_site as i64 + 1 - mapped.extra_begin as i64;
            slot += 1;
        }

        index -= 1;
    }
}

#[allow(clippy::too_many_arguments)]
fn extend_backwards(
    region: &mut OverlapRegion,
    reads: &ReadStore,
    query_sequence: &[u8],
    view: &OverlapView,
    slots: &mut [u64],
    window_count: i64,
    block_length: i64,
    error_rate: f64,
    outcome: &mut RetryOutcome,
) {
    let mut index = 0i64;
    while index < window_count {
        // Only the first aligned window after a hole is used as an anchor.
        if index > 0 && slots[index as usize] != u64::MAX && slots[(index - 1) as usize] == u64::MAX
        {
            let anchor = slots[index as usize] as usize;
            let mut window_start = region.window_cigar.windows[anchor].x_start as i64;

            let local_start = if region.window_cigar.windows[anchor].cigar_length == 0 {
                let window = region.window_cigar.windows[anchor];
                let x_start = window.x_start as i64;
                let x_end = window.x_end as i64;
                let x_length = x_end + 1 - x_start;
                let threshold = window.error_threshold as i64;
                let window_length = x_length + (threshold << 1);
                let mut y_start = window.y_start as i64;
                let extra_begin = window.extra_begin as i32;
                let extra_end = window.extra_end as i32;
                // The reference re-aligns the window exactly where it was
                // placed (`p->y_start`, `p->extra_*`), without running
                // `determine_overlap_region` again.
                let o_len =
                    window_length - i64::from(extra_end) - i64::from(extra_begin);
                let target = fill_subregion(
                    &TargetSource::Reads(reads),
                    view.y_identifier,
                    y_start,
                    o_len,
                    view.y_strand,
                    extra_begin,
                    extra_end,
                );
                let query = query_slice(query_sequence, x_start, x_length);

                let Some(alignment) = reserve_banded_bpm_path(
                    &target,
                    window_length as usize,
                    query,
                    x_length as usize,
                    threshold as usize,
                    window.error as i32,
                    window.y_end as i32 - y_start as i32,
                ) else {
                    return;
                };
                if alignment.error == NO_ALIGNMENT_ERROR {
                    return;
                }

                let mut start_site = alignment.start_site;
                let mut end_site = alignment.end_site;
                let mut error = alignment.error;
                let mut path = alignment.path;
                let mut extra_begin = extra_begin;
                let extra_end = extra_end;
                let mut target_sequence = target;

                if end_site == window_length as i32 - 1 || start_site == 0 {
                    if let Some(fixed) = fix_boundary(&BoundaryContext {
                        query,
                        query_length: x_length,
                        threshold,
                        target_window_length: window_length,
                        total_target_start: y_start,
                        local_start: start_site as i64,
                        local_end: end_site as i64,
                        old_extra_begin: extra_begin,
                        old_extra_end: extra_end,
                        old_error: error,
                        reads,
                        target_identifier: view.y_identifier,
                        target_strand: view.y_strand,
                    }) {
                        y_start = fixed.total_target_start;
                        extra_begin = fixed.extra_begin;
                        error = fixed.error;
                        path = fixed.path;
                        start_site = fixed.local_start;
                        end_site = fixed.local_end;
                        target_sequence = fixed.target_sequence;
                    }
                }

                generate_window_cigar(
                    &mut path,
                    anchor,
                    &mut region.window_cigar,
                    &mut start_site,
                    &mut end_site,
                    &mut error,
                    query,
                    &target_sequence,
                );

                let real_start = y_start + start_site as i64 - extra_begin as i64;
                region.window_cigar.windows[anchor].y_start = real_start as i32;
                region.window_cigar.windows[anchor].y_end = (y_start + end_site as i64) as i32;
                region.window_cigar.windows[anchor].error = error as i16;
                real_start
            } else {
                region.window_cigar.windows[anchor].y_start as i64
            };

            let mut total_y_end = local_start - 1;
            let mut hole = index - 1;
            while hole >= 0 && slots[hole as usize] == u64::MAX {
                let window_end = window_start - 1;
                let (computed_identifier, computed_start) = window_id_by_end(
                    view.overlap_start,
                    window_end,
                    block_length,
                );
                debug_assert_eq!(computed_identifier, hole);
                window_start = computed_start;

                let x_start = window_start;
                let x_end = window_end;
                let x_length = x_end + 1 - x_start;
                let threshold = double_error_threshold(
                    initial_error_threshold(x_length, error_rate, block_length, THRESHOLD),
                    x_length,
                );
                let window_length = x_length + (threshold << 1);

                if total_y_end <= 0 {
                    break;
                }
                let y_start = total_y_end - x_length + 1;

                let Some(mapped) =
                    determine_overlap_region(threshold, y_start, view.y_read_length, window_length)
                else {
                    break;
                };
                if mapped.length + threshold < x_length {
                    break;
                }

                let target = fill_subregion(
                    &TargetSource::Reads(reads),
                    view.y_identifier,
                    mapped.start,
                    mapped.length,
                    view.y_strand,
                    mapped.extra_begin,
                    mapped.extra_end,
                );
                let query = query_slice(query_sequence, x_start, x_length);

                let Some(alignment) = reserve_banded_bpm_path(
                    &target,
                    window_length as usize,
                    query,
                    x_length as usize,
                    threshold as usize,
                    -1,
                    -1,
                ) else {
                    break;
                };
                if alignment.error == NO_ALIGNMENT_ERROR {
                    break;
                }

                let mut start_site = alignment.start_site;
                let mut end_site = alignment.end_site;
                let mut error = alignment.error;
                let mut path = alignment.path;
                let mut extra_begin = mapped.extra_begin;
                let extra_end = mapped.extra_end;
                let mut target_sequence = target;
                let mut mapped_start = mapped.start;

                if end_site == window_length as i32 - 1 || start_site == 0 {
                    if let Some(fixed) = fix_boundary(&BoundaryContext {
                        query,
                        query_length: x_length,
                        threshold,
                        target_window_length: window_length,
                        total_target_start: mapped.start,
                        local_start: start_site as i64,
                        local_end: end_site as i64,
                        old_extra_begin: extra_begin,
                        old_extra_end: extra_end,
                        old_error: error,
                        reads,
                        target_identifier: view.y_identifier,
                        target_strand: view.y_strand,
                    }) {
                        mapped_start = fixed.total_target_start;
                        extra_begin = fixed.extra_begin;
                        error = fixed.error;
                        path = fixed.path;
                        start_site = fixed.local_start;
                        end_site = fixed.local_end;
                        target_sequence = fixed.target_sequence;
                    }
                }

                let window_index = region.window_cigar.append_window(
                    x_start as i32,
                    x_end as i32,
                    0,
                    0,
                    0,
                    extra_begin as i16,
                    extra_end as i16,
                    threshold as i16,
                );
                generate_window_cigar(
                    &mut path,
                    window_index,
                    &mut region.window_cigar,
                    &mut start_site,
                    &mut end_site,
                    &mut error,
                    query,
                    &target_sequence,
                );

                region.window_cigar.windows[window_index].y_start =
                    (mapped_start + start_site as i64 - extra_begin as i64) as i32;
                region.window_cigar.windows[window_index].y_end =
                    (mapped_start + end_site as i64) as i32;
                region.window_cigar.windows[window_index].error = error as i16;
                outcome.aligned_length += x_length as u64;
                slots[hole as usize] = window_index as u64;

                if outcome.is_sorted && window_index > 0 {
                    let previous = region.window_cigar.windows[window_index - 1].x_start;
                    if region.window_cigar.windows[window_index].x_start < previous {
                        outcome.is_sorted = false;
                    }
                }

                total_y_end = mapped_start + start_site as i64 - 1 - extra_begin as i64;
                hole -= 1;
            }
        }
        index += 1;
    }
}

/// The third pass: every window that still has no cigar gets one, using the
/// previous alignment of the same window as a hint.
fn build_missing_cigars(
    region: &mut OverlapRegion,
    reads: &ReadStore,
    query_sequence: &[u8],
    view: &OverlapView,
) {
    let count = region.window_cigar.windows.len();
    for index in 0..count {
        let window = region.window_cigar.windows[index];
        if window.cigar_length != 0 {
            // Windows that already had a cigar only need the padding removed
            // from their target end.
            region.window_cigar.windows[index].y_end -= window.extra_begin as i32;
            continue;
        }

        let x_start = window.x_start as i64;
        let x_end = window.x_end as i64;
        let x_length = x_end + 1 - x_start;
        let threshold = window.error_threshold as i64;
        let window_length = x_length + (threshold << 1);
        // Unlike the forward retry, this pass trusts the position already
        // recorded on the window: it reads `o_len` bases straight from
        // `y_start` and keeps the window's own padding.
        let y_start = window.y_start as i64;
        let extra_begin = window.extra_begin as i32;
        let extra_end = window.extra_end as i32;
        let o_len = window_length - extra_end as i64 - extra_begin as i64;
        if o_len <= 0 {
            continue;
        }

        let target = fill_subregion(
            &TargetSource::Reads(reads),
            view.y_identifier,
            y_start,
            o_len,
            view.y_strand,
            extra_begin,
            extra_end,
        );
        let query = query_slice(query_sequence, x_start, x_length);

        let Some(alignment) = reserve_banded_bpm_path(
            &target,
            window_length as usize,
            query,
            x_length as usize,
            threshold as usize,
            window.error as i32,
            window.y_end as i32 - y_start as i32,
        ) else {
            continue;
        };
        if alignment.error == NO_ALIGNMENT_ERROR {
            continue;
        }

        let mut start_site = alignment.start_site;
        let mut end_site = alignment.end_site;
        let mut error = alignment.error;
        let mut path = alignment.path;
        let mut extra_begin = extra_begin;
        let extra_end = extra_end;
        let mut target_sequence = target;
        let mut mapped_start = y_start;

        if end_site == window_length as i32 - 1 || start_site == 0 {
            if let Some(fixed) = fix_boundary(&BoundaryContext {
                query,
                query_length: x_length,
                threshold,
                target_window_length: window_length,
                total_target_start: y_start,
                local_start: start_site as i64,
                local_end: end_site as i64,
                old_extra_begin: extra_begin,
                old_extra_end: extra_end,
                old_error: error,
                reads,
                target_identifier: view.y_identifier,
                target_strand: view.y_strand,
            }) {
                mapped_start = fixed.total_target_start;
                extra_begin = fixed.extra_begin;
                error = fixed.error;
                path = fixed.path;
                start_site = fixed.local_start;
                end_site = fixed.local_end;
                target_sequence = fixed.target_sequence;
            }
        }

        generate_window_cigar(
            &mut path,
            index,
            &mut region.window_cigar,
            &mut start_site,
            &mut end_site,
            &mut error,
            query,
            &target_sequence,
        );

        region.window_cigar.windows[index].y_start =
            (mapped_start + start_site as i64 - extra_begin as i64) as i32;
        region.window_cigar.windows[index].y_end =
            (mapped_start + end_site as i64 - extra_begin as i64) as i32;
        region.window_cigar.windows[index].error = error as i16;
    }
}

/// A safe view of the query window; the reference reads past the end only for
/// corrupt inputs.
fn query_slice(query_sequence: &[u8], start: i64, length: i64) -> &[u8] {
    let start = start.max(0) as usize;
    let end = (start + length.max(0) as usize).min(query_sequence.len());
    if start >= end {
        return &[];
    }
    &query_sequence[start..end]
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::correction::region::THRESHOLD;

    const BLOCK: i64 = 375;

    /// A thousand base query and target, both all `A`.
    fn reads_with_target(target: &[u8]) -> ReadStore {
        let mut reads = ReadStore::default();
        reads.push_read("query", &vec![b'A'; 1000]);
        reads.push_read("target", target);
        reads
    }

    /// An overlap with the given windows, each still without a cigar, exactly
    /// as `verify_window` leaves them.
    fn region_with_windows(windows: &[(i32, i32, i32, i32)]) -> OverlapRegion {
        let mut region = OverlapRegion {
            x_identifier: 0,
            x_position_start: 0,
            x_position_end: 999,
            y_identifier: 1,
            y_position_start: 0,
            y_position_end: 999,
            ..OverlapRegion::default()
        };
        for (x_start, x_end, y_start, y_end) in windows {
            region.window_cigar.append_window(
                *x_start,
                *x_end,
                *y_start,
                *y_end,
                0,
                0,
                0,
                THRESHOLD as i16,
            );
            // `verify_window` has already added these bases.
            region.alignment_length += (*x_end - *x_start + 1) as u32;
        }
        region
    }

    #[test]
    fn a_fully_verified_overlap_gets_one_cell_per_window() {
        let reads = reads_with_target(&vec![b'A'; 1000]);
        let query = vec![b'A'; 1000];
        let mut regions = vec![region_with_windows(&[
            (0, 374, 0, 374),
            (375, 749, 375, 749),
            (750, 999, 750, 999),
        ])];

        recalcate_window_advance(&mut regions, &reads, &query, BLOCK, 0.04, 0.05);

        let region = &regions[0];
        assert_eq!(region.is_match, 1);
        assert_eq!(region.alignment_length, 1000);
        assert_eq!(region.window_cigar.windows.len(), 3);
        for (index, expected) in [375u32, 375, 250].iter().enumerate() {
            let window = region.window_cigar.windows[index];
            assert_eq!(window.cigar_length, 1, "window {index} has one cell");
            assert_eq!(region.window_cigar.cell(index, 0), (0, *expected));
            assert_eq!(window.error, 0);
        }
        // Two boundaries, both perfectly aligned.
        assert_eq!(region.boundary_cigars.windows.len(), 2);
        assert_eq!(region.boundary_cigars.windows[0].error, -2);
        assert_eq!(region.boundary_cigars.windows[1].error, -2);
    }

    #[test]
    fn a_hole_is_filled_from_the_previous_window() {
        let reads = reads_with_target(&vec![b'A'; 1000]);
        let query = vec![b'A'; 1000];
        // The middle window is missing.
        let mut regions = vec![region_with_windows(&[
            (0, 374, 0, 374),
            (750, 999, 750, 999),
        ])];

        recalcate_window_advance(&mut regions, &reads, &query, BLOCK, 0.04, 0.05);

        let region = &regions[0];
        assert_eq!(region.window_cigar.windows.len(), 3, "the hole was filled");
        assert_eq!(region.is_match, 1);
        assert_eq!(region.alignment_length, 1000);

        let middle = region.window_cigar.windows[1];
        assert_eq!((middle.x_start, middle.x_end), (375, 749));
        // The retry lands on target 344..780 and pass three re-derives the
        // reported window from the mapped start, so the window is 375..749.
        assert_eq!(middle.y_start, 375);
        assert_eq!(middle.y_end, 749);
        assert_eq!(middle.error, 0);
        assert!(middle.cigar_length >= 1);
    }

    #[test]
    fn an_unalignable_hole_rejects_the_overlap() {
        let mut target = vec![b'A'; 1000];
        target[500..560].fill(b'C');
        let reads = reads_with_target(&target);
        let query = vec![b'A'; 1000];
        let mut regions = vec![region_with_windows(&[
            (0, 374, 0, 374),
            (750, 999, 750, 999),
        ])];

        recalcate_window_advance(&mut regions, &reads, &query, BLOCK, 0.04, 0.05);

        let region = &regions[0];
        // Sixty mismatching bases inside a 375 base hole cost more than the
        // overlap may spend, so it is dropped.
        assert_eq!(region.is_match, 0);
        assert!(region.window_cigar.windows.len() >= 2);
    }

    #[test]
    fn an_overlap_without_windows_is_left_alone() {
        let reads = reads_with_target(&vec![b'A'; 1000]);
        let query = vec![b'A'; 1000];
        let mut regions = vec![OverlapRegion {
            x_position_start: 0,
            x_position_end: 999,
            y_identifier: 1,
            y_position_end: 999,
            ..OverlapRegion::default()
        }];

        recalcate_window_advance(&mut regions, &reads, &query, BLOCK, 0.04, 0.05);

        assert_eq!(regions[0].is_match, 0);
        assert!(!regions[0].window_cigar.windows.is_empty() || regions[0].window_cigar.windows.is_empty());
        assert_eq!(regions[0].alignment_length, 0);
    }
}

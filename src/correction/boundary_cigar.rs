//! `calculate_boundary_cigars` from `Correct.cpp`.
//!
//! Two neighbouring windows are aligned independently, so the bases around
//! their boundary are covered by neither cigar. This pass aligns a window that
//! straddles the boundary (half of it on each side) and keeps the result only
//! if it explains the boundary better than the two window cigars do. The second
//! round of correction uses these cigars to change bases that the first round
//! could not decide.

use crate::alignment::banded_bpm::NO_ALIGNMENT_ERROR;
use crate::alignment::banded_path::reserve_banded_bpm_path;
use crate::correction::cigar_generator::generate_window_cigar;
use crate::correction::recalculate::double_error_threshold;
use crate::correction::region::{adjust_threshold, determine_overlap_region, fill_subregion};
use crate::overlap::region::OverlapRegion;
use crate::correction::target::TargetSource;
use crate::reads::read_store::ReadStore;

/// `scanLen`: query bases inspected on each side of a boundary when deciding
/// whether it needs a cigar at all.
const SCAN_LENGTH: i64 = 10;

/// `boundaryLen`: total width of the straddling window.
const BOUNDARY_LENGTH: i64 = 200;

/// `alpha`: slack that the straddling alignment must beat.
const ALPHA: i64 = 1;

/// `calculate_boundary_cigars`: fills `region.boundary_cigars` with one entry
/// per pair of neighbouring windows.
///
/// An entry whose `error` is `-1` was rejected (or never needed), one with
/// `error == -2` means both window cigars already agree exactly at the
/// boundary.
pub fn calculate_boundary_cigars(
    region: &mut OverlapRegion,
    reads: &ReadStore,
    query_sequence: &[u8],
    error_rate: f64,
) {
    let window_cigar = std::mem::take(&mut region.window_cigar);
    let mut boundary = std::mem::take(&mut region.boundary_cigars);

    let windows = &window_cigar.windows;
    let window_count = windows.len();
    if window_count == 0 {
        region.window_cigar = window_cigar;
        region.boundary_cigars = boundary;
        return;
    }

    boundary.reset_windows(window_count - 1);

    let y_identifier = region.y_identifier as usize;
    let y_strand = region.y_strand as u8;
    let y_read_length = reads.read_length(y_identifier) as i64;
    let single_side_length = BOUNDARY_LENGTH / 2;
    let force_useless_side = single_side_length / 2;

    for index in 0..window_count - 1 {
        let current = windows[index];
        let next = windows[index + 1];

        // A boundary next to an unaligned window has nothing to connect.
        if current.y_end == -1 || next.y_end == -1 {
            boundary.windows[index].error = -1;
            boundary.windows[index].y_end = -1;
            continue;
        }

        let mut y_distance = next.y_start as i64 - current.y_end as i64 - 1;

        // Both windows align exactly at the boundary and neither has an error
        // within the scanned bases: the two cigars already join up.
        if y_distance == 0 {
            let (_, backward_error) = window_cigar.scan_error(index, SCAN_LENGTH, 1);
            let (_, forward_error) = window_cigar.scan_error(index + 1, SCAN_LENGTH, 0);
            if backward_error == 0 && forward_error == 0 {
                boundary.windows[index].error = -2;
                boundary.windows[index].y_end = -1;
                continue;
            }
        }

        let (mut y_start, mut x_start) = if current.y_end != -1 {
            (current.y_end as i64, current.x_end as i64)
        } else {
            (next.y_start as i64, next.x_start as i64)
        };

        // The reference stores the original positions here before they are
        // overwritten; nothing reads them back.
        boundary.windows[index].extra_begin = x_start as i16;
        boundary.windows[index].extra_end = y_start as i16;

        let left_length = (x_start - current.x_start as i64)
            .min(y_start)
            .min(single_side_length);
        let right_length = (next.x_end as i64 + 1 - x_start)
            .min(y_read_length - y_start)
            .min(single_side_length);

        let x_length = left_length + right_length;
        x_start -= left_length;
        let x_end = x_start + x_length - 1;
        y_start -= left_length;

        let mut threshold = (x_length as f64 * error_rate) as i64;
        threshold = adjust_threshold(threshold, x_length);
        threshold = double_error_threshold(threshold, x_length);

        let y_length = x_length + (threshold << 1);
        let Some(mapped) =
            determine_overlap_region(threshold, y_start, y_read_length, y_length)
        else {
            boundary.windows[index].error = -1;
            boundary.windows[index].y_end = -1;
            continue;
        };
        if mapped.length < x_length {
            boundary.windows[index].error = -1;
            boundary.windows[index].y_end = -1;
            continue;
        }

        let target_sequence = fill_subregion(
            &TargetSource::Reads(reads),
            y_identifier,
            mapped.start,
            mapped.length,
            y_strand,
            mapped.extra_begin,
            mapped.extra_end,
        );
        let query = &query_sequence[x_start.max(0) as usize
            ..(x_start.max(0) as usize + x_length as usize).min(query_sequence.len())];

        let Some(mut alignment) = reserve_banded_bpm_path(
            &target_sequence,
            y_length as usize,
            query,
            x_length as usize,
            threshold as usize,
            -1,
            -1,
        ) else {
            boundary.windows[index].error = -1;
            boundary.windows[index].y_end = -1;
            continue;
        };

        if alignment.error == NO_ALIGNMENT_ERROR {
            boundary.windows[index].error = -1;
            boundary.windows[index].y_end = -1;
            continue;
        }

        boundary.windows[index].x_start = x_start as i32;
        boundary.windows[index].x_end = x_end as i32;

        let mut start_site = alignment.start_site;
        let mut end_site = alignment.end_site;
        let mut error = alignment.error;
        generate_window_cigar(
            &mut alignment.path,
            index,
            &mut boundary,
            &mut start_site,
            &mut end_site,
            &mut error,
            query,
            &target_sequence,
        );

        if y_distance < 0 {
            y_distance = -y_distance;
        }

        // The outer halves of the first and last boundary overlap a single
        // window, so they carry no new information.
        let mut left_useless_side = force_useless_side;
        let mut right_useless_side = force_useless_side;
        if index == 0 && x_start == windows[0].x_start as i64 {
            left_useless_side = 0;
        }
        if index == window_count - 2 && x_end == windows[window_count - 1].x_end as i64 {
            right_useless_side = 0;
        }

        if left_length <= left_useless_side || right_length <= right_useless_side {
            boundary.windows[index].error = -1;
            boundary.windows[index].y_end = -1;
            boundary.cells.truncate(boundary.windows[index].cigar_start as usize);
            continue;
        }

        let (_, boundary_error) = boundary.scan_error_interval(
            index,
            left_useless_side,
            x_length - right_useless_side - 1,
        );
        let (_, backward_error) =
            window_cigar.scan_error(index, left_length - left_useless_side, 1);
        let (_, forward_error) =
            window_cigar.scan_error(index + 1, right_length - right_useless_side, 0);

        // Keep the new cigar only when it explains the boundary at least as
        // well as the two window cigars already do.
        if forward_error + backward_error + y_distance + ALPHA < boundary_error {
            boundary.windows[index].error = -1;
            boundary.windows[index].y_end = -1;
            boundary.cells.truncate(boundary.windows[index].cigar_start as usize);
            continue;
        }

        boundary.windows[index].error = error as i16;
        // The reference overwrites its `y_start` with the mapped target start,
        // so the final positions are measured from there.
        boundary.windows[index].y_start =
            (mapped.start + start_site as i64 - mapped.extra_begin as i64) as i32;
        boundary.windows[index].y_end =
            (mapped.start + end_site as i64 - mapped.extra_begin as i64) as i32;
        boundary.windows[index].x_start = x_start as i32;
        boundary.windows[index].x_end = x_end as i32;
        boundary.windows[index].extra_begin = left_useless_side as i16;
        boundary.windows[index].extra_end = right_useless_side as i16;
    }

    region.window_cigar = window_cigar;
    region.boundary_cigars = boundary;
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::correction::region::THRESHOLD;

    /// Two reads of a thousand bases; the target is given by the caller.
    fn reads_with_target(target: &[u8]) -> ReadStore {
        let mut reads = ReadStore::default();
        reads.push_read("query", &vec![b'A'; 1000]);
        reads.push_read("target", target);
        reads
    }

    /// The window layout of the probe: the first window matches for 495 bases
    /// and then has five mismatches, the second one mirrors that.
    fn two_windows(boundary_errors: bool) -> OverlapRegion {
        let mut region = OverlapRegion {
            x_identifier: 0,
            x_position_start: 0,
            x_position_end: 999,
            y_identifier: 1,
            y_position_start: 0,
            y_position_end: 999,
            ..OverlapRegion::default()
        };

        let first = region.window_cigar.append_window(0, 499, 0, 499, 5, 0, 0, THRESHOLD as i16);
        let cells = region.window_cigar.cells.len();
        region.window_cigar.push_cell(0, 495);
        region.window_cigar.push_cell(1, 5);
        region.window_cigar.assign_cells(first, cells);

        let second = if boundary_errors {
            region
                .window_cigar
                .append_window(500, 999, 500, 999, 5, 0, 0, THRESHOLD as i16)
        } else {
            region.window_cigar.append_window(
                500,
                999,
                -1,
                -1,
                -1,
                0,
                0,
                THRESHOLD as i16,
            )
        };
        if boundary_errors {
            let cells = region.window_cigar.cells.len();
            region.window_cigar.push_cell(1, 5);
            region.window_cigar.push_cell(0, 495);
            region.window_cigar.assign_cells(second, cells);
        }
        region
    }

    /// Golden values printed by the real `calculate_boundary_cigars` for a one
    /// thousand base overlap split into two windows of five hundred bases.
    #[test]
    fn matches_the_reference_boundary_cigars() {
        let query = vec![b'A'; 1000];

        // A perfectly matching boundary: the straddling window is a single
        // 200 base match, with fifty "useless" bases kept on each side.
        let reads = reads_with_target(&vec![b'A'; 1000]);
        let mut region = two_windows(true);
        calculate_boundary_cigars(&mut region, &reads, &query, 0.04);

        let window = region.boundary_cigars.windows[0];
        assert_eq!(window.error, 0);
        assert_eq!((window.y_start, window.y_end), (399, 598));
        assert_eq!((window.x_start, window.x_end), (399, 598));
        assert_eq!((window.extra_begin, window.extra_end), (50, 50));
        assert_eq!(window.cigar_length, 1);
        assert_eq!(region.boundary_cigars.cells, vec![200]);

        // Ten mismatching bases right at the boundary: the trimmed cigar keeps
        // the run in the middle.
        let mut target = vec![b'A'; 1000];
        target[495..505].fill(b'C');
        let reads = reads_with_target(&target);
        let mut region = two_windows(true);
        calculate_boundary_cigars(&mut region, &reads, &query, 0.04);

        let window = region.boundary_cigars.windows[0];
        assert_eq!(window.error, 10);
        assert_eq!((window.y_start, window.y_end), (399, 598));
        assert_eq!((window.x_start, window.x_end), (399, 598));
        let cells: Vec<u16> = region.boundary_cigars.cells.clone();
        assert_eq!(cells, vec![96, (1 << 14) | 10, 94]);

        // The second window is unaligned, so the boundary is skipped and the
        // entry keeps the "not computed" marker.
        let reads = reads_with_target(&vec![b'A'; 1000]);
        let mut region = two_windows(false);
        calculate_boundary_cigars(&mut region, &reads, &query, 0.04);

        let window = region.boundary_cigars.windows[0];
        assert_eq!(window.error, -1);
        assert_eq!(window.y_end, -1);
        assert_eq!((window.x_start, window.x_end), (0, 0));
        assert_eq!((window.extra_begin, window.extra_end), (0, 0));
        assert!(region.boundary_cigars.cells.is_empty());
    }

    #[test]
    fn a_perfect_boundary_needs_no_cigar() {
        // Both window cigars are a single match cell, so the ten bases at each
        // side scan clean and the boundary is marked `-2`.
        let reads = reads_with_target(&vec![b'A'; 1000]);
        let query = vec![b'A'; 1000];
        let mut region = OverlapRegion {
            x_identifier: 0,
            x_position_start: 0,
            x_position_end: 999,
            y_identifier: 1,
            y_position_start: 0,
            y_position_end: 999,
            ..OverlapRegion::default()
        };
        for start in [0i32, 500] {
            let index = region.window_cigar.append_window(
                start,
                start + 499,
                start,
                start + 499,
                0,
                0,
                0,
                THRESHOLD as i16,
            );
            let cells = region.window_cigar.cells.len();
            region.window_cigar.push_cell(0, 500);
            region.window_cigar.assign_cells(index, cells);
        }

        calculate_boundary_cigars(&mut region, &reads, &query, 0.04);

        assert_eq!(region.boundary_cigars.windows.len(), 1);
        assert_eq!(region.boundary_cigars.windows[0].error, -2);
        assert_eq!(region.boundary_cigars.windows[0].y_end, -1);
        assert!(region.boundary_cigars.cells.is_empty());
    }

    #[test]
    fn the_window_list_is_restored_after_the_pass() {
        let reads = reads_with_target(&vec![b'A'; 1000]);
        let query = vec![b'A'; 1000];
        let mut region = two_windows(true);
        let before = region.window_cigar.clone();

        calculate_boundary_cigars(&mut region, &reads, &query, 0.04);

        assert_eq!(region.window_cigar, before, "the window cigars are untouched");
    }
}

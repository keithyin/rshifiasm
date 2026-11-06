//! Chimeric-read detection, the counterpart of `collect_sides`,
//! `collect_contain`, `boundary_verify`, `verify_single_window`,
//! `intersection_check_by_base` and `detect_chimeric_reads` in
//! `Overlaps.cpp`/`Correct.cpp`.
//!
//! A read is chimeric when its left- and right-anchored overlaps cannot be
//! joined by any overlap spanning the gap. The gap is first checked with the
//! stored overlap intervals and, when an overlap does span it, with a banded
//! edit-distance alignment of the raw bases.

use crate::alignment::banded_bpm::{reserve_banded_bpm, NO_ALIGNMENT_ERROR};
use crate::assembly::hit::{delete_all_edges, CoverageCut, PRIMARY_LABEL};
use crate::command_line::options::AssemblyOptions;
use crate::correction::region::{adjust_threshold, determine_overlap_region, fill_subregion};
use crate::overlap::paf::PafHit;
use crate::reads::ReadStore;

/// `WINDOW`, the longest base-level region verified in one alignment.
const WINDOW: u32 = 375;

/// The left- or right-extent accumulator of `detect_chimeric_reads`, the
/// `s`/`e` pair of a `ma_sub_t`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct Side {
    start: u32,
    end: u32,
}

/// `collect_sides`: the union of the left- and right-anchored overlaps.
///
/// A left-anchored overlap starts at query base 0; a right-anchored one ends
/// at the last query base. An overlap spanning the whole read counts for both.
fn collect_sides(
    read_identifier: usize,
    sources: &[Vec<PafHit>],
    read_length: i64,
    max_left: &mut Side,
    max_right: &mut Side,
) {
    for hit in &sources[read_identifier] {
        if hit.deleted {
            continue;
        }

        let query_start = hit.query_start;
        let query_end = hit.query_end;

        if query_start == 0 {
            if query_start < max_left.start {
                max_left.start = query_start;
            }
            if query_end > max_left.end {
                max_left.end = query_end;
            }
        }

        if query_end as u64 == read_length as u64 {
            if query_start < max_right.start {
                max_right.start = query_start;
            }
            if query_end > max_right.end {
                max_right.end = query_end;
            }
        }
    }
}

/// `collect_contain`: extends the left and right extents with contained
/// overlaps that cross the current boundary by more than `overlap_rate`.
fn collect_contain(
    sources: &[Vec<PafHit>],
    read_identifier: usize,
    read_length: i64,
    max_left: &mut Side,
    max_right: &mut Side,
    overlap_rate: f32,
) {
    let mut new_left_end = max_left.end as i64;
    let mut new_right_start = max_right.start as i64;

    for hit in &sources[read_identifier] {
        if hit.deleted {
            continue;
        }

        let query_start = hit.query_start;
        let query_end = hit.query_end;
        if query_start == 0 || query_end as u64 == read_length as u64 {
            continue;
        }

        if query_start < max_left.end
            && query_end > max_left.end
            && max_left.end.wrapping_sub(query_start) as f32
                > overlap_rate * (query_end.wrapping_sub(query_start) as f32)
        {
            if query_end > max_left.end && (query_end as i64) > new_left_end {
                new_left_end = query_end as i64;
            }
        }

        if query_start < max_right.start
            && query_end > max_right.start
            && query_end.wrapping_sub(max_right.start) as f32
                > overlap_rate * (query_end.wrapping_sub(query_start) as f32)
        {
            if query_start < max_right.start && (query_start as i64) < new_right_start {
                new_right_start = query_start as i64;
            }
        }
    }

    max_left.end = new_left_end as u32;
    max_right.start = new_right_start as u32;
}

/// `boundary_verify`: whether the raw bases around `[start, end)` on the query
/// match the target through the overlap `hit`.
///
/// `true` corresponds to the reference's return value 1 (the alignment
/// succeeded).
pub fn boundary_verify(
    x_interval_start: u32,
    x_interval_end: u32,
    hit: &PafHit,
    reads: &ReadStore,
    options: &AssemblyOptions,
) -> bool {
    let direction = hit.reverse as u32;
    let query_start = hit.query_start;
    let query_identifier = hit.query_identifier as usize;
    let target_identifier = hit.target_identifier as usize;
    let target_length = reads.read_length(target_identifier) as i64;

    let target_start = if direction == 1 {
        // The target is on the opposite strand, so its coordinates are
        // mirrored. `yLen - (te - 1) - 1` is computed in the reference's
        // signed 64-bit arithmetic.
        (target_length - (hit.target_end as i64 - 1) - 1) as u32
    } else {
        hit.target_start
    };

    let interval_length = x_interval_end.wrapping_sub(x_interval_start);
    let y_interval_start = x_interval_start
        .wrapping_sub(query_start)
        .wrapping_add(target_start);
    if y_interval_start as i64 >= target_length {
        return false;
    }

    let mut y_interval_end = y_interval_start
        .wrapping_add(interval_length)
        .wrapping_sub(1);
    if y_interval_end as i64 >= target_length {
        y_interval_end = (target_length - 1) as u32;
    }
    if y_interval_end < y_interval_start {
        return false;
    }

    let y_interval_length = y_interval_end
        .wrapping_sub(y_interval_start)
        .wrapping_add(1);

    if y_interval_length <= WINDOW {
        return verify_single_window(
            y_interval_start as i64,
            y_interval_end as i64,
            target_start as i64,
            query_start as i64,
            target_identifier,
            query_identifier,
            direction as u8,
            reads,
            options,
        );
    }

    if !verify_single_window(
        y_interval_start as i64,
        y_interval_start.wrapping_add(WINDOW).wrapping_sub(1) as i64,
        target_start as i64,
        query_start as i64,
        target_identifier,
        query_identifier,
        direction as u8,
        reads,
        options,
    ) {
        return false;
    }

    verify_single_window(
        y_interval_end.wrapping_sub(WINDOW).wrapping_add(1) as i64,
        y_interval_end as i64,
        target_start as i64,
        query_start as i64,
        target_identifier,
        query_identifier,
        direction as u8,
        reads,
        options,
    )
}

/// `verify_single_window`: one banded alignment of the query window against
/// the target read.
///
/// The query read (`y_id` here, mirroring the reference's swapped arguments)
/// is always read on the forward strand; the target (`x_id`) is read on
/// `x_strand`.
fn verify_single_window(
    x_start: i64,
    x_end: i64,
    overlap_x_start: i64,
    overlap_y_start: i64,
    x_identifier: usize,
    y_identifier: usize,
    x_strand: u8,
    reads: &ReadStore,
    options: &AssemblyOptions,
) -> bool {
    let x_length = (x_end - x_start + 1) as i32;
    let threshold = (x_length as f64
        * options.maximum_overlap_divergence_error_correction) as i32;
    let threshold = adjust_threshold(threshold as i64, x_length as i64) as i32;

    let y_start = (x_start - overlap_x_start) + overlap_y_start;
    let window_length = x_length as i64 + ((threshold << 1) as i64);

    let y_length = reads.read_length(y_identifier) as i64;
    let Some(region) = determine_overlap_region(
        threshold as i64,
        y_start,
        y_length,
        window_length,
    ) else {
        return false;
    };

    let y_string = fill_subregion(
        reads,
        y_identifier,
        region.start,
        region.length,
        0,
        region.extra_begin,
        region.extra_end,
    );
    let x_string = reads.read_sub_region(
        x_identifier,
        x_start,
        x_length as i64,
        x_strand,
    );

    let alignment = reserve_banded_bpm(
        &y_string,
        window_length as usize,
        &x_string,
        x_length as usize,
        threshold as usize,
    );

    alignment.error != NO_ALIGNMENT_ERROR
}

/// `intersection_check_by_base`: whether some live overlap that spans the gap
/// proves the read is a chimera.
///
/// A hit spanning `[interval_s, interval_e)` that is flagged as an exact
/// region, or whose raw bases fail [`boundary_verify`], is evidence against
/// the junction. The reference returns 1 in that case and 0 otherwise.
fn intersection_check_by_base(
    sources: &[Vec<PafHit>],
    read_identifier: usize,
    interval_start: u32,
    interval_end: u32,
    reads: &ReadStore,
    options: &AssemblyOptions,
) -> bool {
    for hit in &sources[read_identifier] {
        if hit.deleted {
            continue;
        }

        if hit.query_start <= interval_start && hit.query_end >= interval_end {
            if hit.exact_region
                || !boundary_verify(interval_start, interval_end, hit, reads, options)
            {
                return true;
            }
        }
    }

    false
}

/// `detect_chimeric_reads`: removes every edge of the reads whose two ends
/// are not joined by a supported overlap.
///
/// `shift_rate` is `maximum_overlap_divergence_final * 2.0`. The ultra-long
/// overlap list of the reference (`all_ul_t *x`) is only used by the `--ar`
/// path, which is not part of the default configuration and is therefore not
/// passed here.
pub fn detect_chimeric_reads(
    sources: &mut [Vec<PafHit>],
    read_lengths: &[u64],
    coverage_cut: &mut [CoverageCut],
    shift_rate: f32,
    reads: &ReadStore,
    options: &AssemblyOptions,
) {
    let start_time = crate::logging::realtime_seconds();
    let read_count = sources.len();
    let mut simple_removed = 0i64;
    let mut complex_removed = 0i64;
    let mut complex_removed_real = 0i64;

    for i in 0..read_count {
        coverage_cut[i].class = PRIMARY_LABEL;
        let read_length = read_lengths[i] as i64;

        let mut max_left = Side {
            start: read_length as u32,
            end: 0,
        };
        let mut max_right = Side {
            start: read_length as u32,
            end: 0,
        };

        collect_sides(i, sources, read_length, &mut max_left, &mut max_right);

        // No overlap anchored at one of the two ends: an end node, not a
        // chimera.
        if max_left.start as i64 == read_length || max_right.start as i64 == read_length {
            continue;
        }

        collect_contain(
            sources,
            i,
            read_length,
            &mut max_left,
            &mut max_right,
            0.1,
        );

        // The two sides overlap by more than the allowed divergence: a normal
        // read.
        if max_left.end > max_right.start
            && max_left.end.wrapping_sub(max_right.start) as f32
                >= read_length as f32 * shift_rate
        {
            continue;
        }

        // The two sides do not reach each other: a simple chimera.
        if max_left.end <= max_right.start {
            delete_all_edges(sources, coverage_cut, i as u32);
            simple_removed += 1;
            continue;
        }

        let interval_start = max_right.start;
        let interval_end = max_left.end;

        if intersection_check_by_base(
            sources,
            i,
            interval_start,
            interval_end,
            reads,
            options,
        ) {
            delete_all_edges(sources, coverage_cut, i as u32);
            complex_removed_real += 1;
        }

        complex_removed += 1;
    }

    eprintln!(
        "[M::detect_chimeric_reads] takes {:.2} s, n_simple_remove: {}, n_complex_remove: {}/{}\n",
        crate::logging::realtime_seconds() - start_time,
        simple_removed,
        complex_removed_real,
        complex_removed
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(
        query: u32,
        query_start: u32,
        query_end: u32,
        target: u32,
        target_start: u32,
        target_end: u32,
    ) -> PafHit {
        PafHit {
            query_identifier: query,
            query_start,
            query_end,
            target_identifier: target,
            target_start,
            target_end,
            matching_length: 1,
            block_length: 100,
            reverse: false,
            without_large_indel: true,
            exact_region: false,
            deleted: false,
        }
    }

    fn cut() -> CoverageCut {
        CoverageCut {
            start: 0,
            end: 100,
            deleted: false,
            class: 0,
        }
    }

    #[test]
    fn sides_take_both_anchored_overlaps() {
        let sources = vec![vec![
            hit(0, 0, 40, 1, 0, 40),
            hit(0, 60, 100, 2, 0, 40),
            hit(0, 10, 30, 3, 0, 20),
        ]];
        let mut max_left = Side { start: 100, end: 0 };
        let mut max_right = Side { start: 100, end: 0 };
        collect_sides(0, &sources, 100, &mut max_left, &mut max_right);
        assert_eq!(max_left, Side { start: 0, end: 40 });
        assert_eq!(max_right, Side { start: 60, end: 100 });
    }

    #[test]
    fn contain_extends_the_sides() {
        let sources = vec![vec![
            hit(0, 0, 40, 1, 0, 40),
            hit(0, 60, 100, 2, 0, 40),
            hit(0, 20, 70, 3, 0, 50),
            hit(0, 30, 80, 4, 0, 50),
        ]];
        let mut max_left = Side { start: 0, end: 40 };
        let mut max_right = Side { start: 60, end: 100 };
        collect_contain(&sources, 0, 100, &mut max_left, &mut max_right, 0.1);
        // [20, 70) crosses both boundaries; [30, 80) then extends the left
        // end further, while the right start keeps the smallest crossing
        // start seen.
        assert_eq!(max_left.end, 80);
        assert_eq!(max_right.start, 20);
    }

    #[test]
    fn a_simple_chimera_loses_its_edges() {
        let mut sources = vec![
            vec![
                hit(0, 0, 40, 1, 0, 40),
                hit(0, 60, 100, 2, 0, 40),
            ],
            vec![hit(1, 0, 40, 0, 0, 40)],
            vec![hit(2, 0, 40, 0, 60, 100)],
        ];
        let mut cuts = vec![cut(), cut(), cut()];
        let reads = ReadStore::default();
        let options = AssemblyOptions::default();

        detect_chimeric_reads(
            &mut sources,
            &[100, 40, 40],
            &mut cuts,
            0.06,
            &reads,
            &options,
        );

        assert!(cuts[0].deleted);
        assert!(sources[0].iter().all(|hit| hit.deleted));
        assert!(sources[1][0].deleted);
        assert!(sources[2][0].deleted);
        assert_eq!(cuts[0].class, PRIMARY_LABEL);
    }

    #[test]
    fn a_joined_read_is_not_chimeric() {
        let mut sources = vec![
            vec![
                hit(0, 0, 40, 1, 0, 40),
                hit(0, 60, 100, 2, 0, 40),
                hit(0, 0, 100, 3, 0, 100),
            ],
            vec![hit(1, 0, 40, 0, 0, 40)],
            vec![hit(2, 0, 40, 0, 60, 100)],
            vec![hit(3, 0, 100, 0, 0, 100)],
        ];
        let mut cuts = vec![cut(), cut(), cut(), cut()];
        let reads = ReadStore::default();
        let options = AssemblyOptions::default();

        detect_chimeric_reads(
            &mut sources,
            &[100, 40, 40, 100],
            &mut cuts,
            0.06,
            &reads,
            &options,
        );

        // The whole-read overlap joins both sides, so read 0 stays.
        assert!(!cuts[0].deleted);
        assert!(!sources[0][0].deleted);
    }
}

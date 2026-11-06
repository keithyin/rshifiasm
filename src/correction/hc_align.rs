//! The post-chain refinement of the new HiFi correction path, the counterpart
//! of `align_hc_ed_post_extz` and its window helpers in `Correct.cpp`.
//!
//! The older path (`verify_window` in [`super::verify`]) verifies one fixed
//! sliding window at a time and stores the coordinates only. This path walks
//! the overlap in blocks of the same length, aligns each block against the
//! mapped target region with the absent-diagonal semi-global variant, fills the
//! blocks the direct pass skipped by walking backwards from the next aligned
//! block, and records a packed CIGAR for every block it reconstructs. The
//! direct blocks carry no CIGAR: their cells are produced later by the
//! consensus stage from the stored coordinates.
//!
//! Coordinates follow the reference exactly. `y_start` of a direct block is the
//! start of the recovered target region (which is `threshold` bases before the
//! mapped position), while `y_end` is the end of the aligned part. The comment
//! in the reference, "ts do not have aux_beg, while te has", refers to that
//! asymmetry.

use crate::alignment::banded_bpm::reserve_banded_bpm_extension;
use crate::alignment::edit_distance::{
    EditDistanceResult, TraceCigar, semi_global_64_absent_diagonal,
    semi_global_64_absent_diagonal_trace,
};
use crate::correction::recalculate::{
    double_error_threshold, initial_error_threshold, similarity_pass,
};
use crate::correction::region::{THRESHOLD_MAXIMUM, TargetRegion, adjust_threshold,
                                determine_overlap_region, fill_subregion};
use crate::correction::window::{
    WindowCigar, WindowList, number_of_windows, window_id_by_end, window_id_by_start,
};
use crate::overlap::region::OverlapRegion;
use crate::reads::read_store::ReadStore;

/// `bit_extz_t` as far as the post-chain alignment uses it: the result of the
/// last alignment plus the packed traceback it produced.
///
/// The reference keeps the DP state (bit vectors, path) inside the same object
/// and only clears the error between calls; every call here re-runs the DP from
/// scratch, so the Rust object only needs the public result and the CIGAR.
#[derive(Debug, Clone, Default)]
pub struct BitExtz {
    /// Error, threshold and boundaries of the last alignment.
    pub result: EditDistanceResult,
    /// Packed traceback of the last trace alignment, `operation << 14 | length`.
    pub cigar: TraceCigar,
}

impl BitExtz {
    pub fn new() -> Self {
        Self::default()
    }

    /// `is_align`
    pub fn is_aligned(&self) -> bool {
        self.result.error <= self.result.threshold
    }

    /// `clear_align`
    pub(crate) fn clear_align(&mut self) {
        self.result.error = i32::MAX;
    }

    /// `ed_band_cal_semi_64_w_absent_diag`
    fn semi_global_absent_diagonal(
        &mut self,
        pattern: &[u8],
        text: &[u8],
        threshold: i64,
        absolute_diagonal: i64,
    ) {
        semi_global_64_absent_diagonal(
            pattern,
            pattern.len() as i32,
            text,
            text.len() as i32,
            threshold as i32,
            absolute_diagonal as i32,
            &mut self.result,
        );
    }

    /// `ed_band_cal_semi_64_w_absent_diag_trace`
    fn semi_global_absent_diagonal_trace(
        &mut self,
        pattern: &[u8],
        text: &[u8],
        threshold: i64,
        absolute_diagonal: i64,
    ) {
        semi_global_64_absent_diagonal_trace(
            pattern,
            pattern.len() as i32,
            text,
            text.len() as i32,
            threshold as i32,
            absolute_diagonal as i32,
            &mut self.result,
            &mut self.cigar,
        );
    }
}

/// Parameters shared by every overlap of one read, `align_hc_ed_post_extz`'s
/// scalar arguments.
#[derive(Debug, Clone, Copy)]
pub struct HcAlignParameters {
    /// Divergence allowed per window, `e_rate`.
    pub error_rate: f64,
    /// Window length, `w_l`.
    pub window_length: i64,
    /// Fraction of the overlap that must be aligned, `ovlp_cut`.
    pub overlap_cutoff: f64,
    /// Accept even when the overlap fraction check fails, `force_aln`.
    pub force_alignment: bool,
}

/// `align_hc_ed_post_extz`: refines one overlap window by window.
///
/// Returns `false` when the overlap fails the final fraction check; the caller
/// then drops it. On success the overlap's window list and `alignment_length`
/// are updated; `is_match` stays `0` until the caller sets it.
pub fn align_hc_ed_post_extz(
    region: &mut OverlapRegion,
    extz: &mut BitExtz,
    reads: &ReadStore,
    query: &[u8],
    parameters: &HcAlignParameters,
) -> bool {
    let window_length = parameters.window_length;
    let x_start = region.x_position_start as i64;
    let x_end = region.x_position_end as i64;
    let y_start = region.y_position_start as i64;
    let target_identifier = region.y_identifier as usize;
    let target_strand = region.y_strand as u8;
    let target_total_length = reads.read_length(target_identifier) as i64;

    region.window_cigar.windows.clear();
    region.is_match = 0;
    region.alignment_length = 0;

    let window_count = number_of_windows(x_start, x_end + 1, window_length);
    let origin = (x_start / window_length) * window_length;
    let (mut query_start, mut query_end) =
        window_se_by_normalized_start(region, origin, window_length);

    for _ in 0..window_count {
        let query_length = 1 + query_end - query_start;
        let mut threshold =
            adjust_threshold((query_length as f64 * parameters.error_rate) as i64, query_length);
        if threshold > THRESHOLD_MAXIMUM {
            threshold = THRESHOLD_MAXIMUM;
        }

        // The mapped target position, shifted by the gap the chain already
        // knows about.
        let mut target_start = (query_start - x_start) + y_start;
        target_start += region
            .fake_cigar
            .y_start_offset(query_start)
            .expect("the fake cigar must cover every query window start");

        let alignment_length = query_length + (threshold << 1);
        if let Some(target) = determine_overlap_region(
            threshold,
            target_start,
            target_total_length,
            alignment_length,
        ) {
            let target_sequence = reads.read_sub_region(
                target_identifier,
                target.start,
                target.length,
                target_strand,
            );
            let query_sequence = &query[query_start as usize..(query_start + query_length) as usize];

            extz.semi_global_absent_diagonal(
                &target_sequence,
                query_sequence,
                threshold,
                target.extra_begin as i64,
            );

            if extz.is_aligned()
                && !push_hc_wlst_exz(
                    region,
                    extz,
                    reads,
                    query,
                    THRESHOLD_MAXIMUM,
                    query_start,
                    query_end,
                    target.start,
                    target.start + extz.result.pattern_end as i64,
                    target_total_length,
                    target.extra_begin as i64,
                    target.extra_end as i64,
                    parameters,
                )
            {
                return false;
            }
        }

        query_start = query_end + 1;
        query_end = query_start + window_length - 1;
        if query_end >= x_end {
            query_end = x_end;
        }
    }

    if !parameters.force_alignment
        && !pass_qovlp(
            x_end + 1 - x_start,
            region.alignment_length as i64,
            parameters.overlap_cutoff,
        )
    {
        return false;
    }
    true
}

/// `push_hc_wlst_exz`: records the directly aligned window `p` after filling
/// the blocks between it and the previously accepted window.
#[allow(clippy::too_many_arguments)]
fn push_hc_wlst_exz(
    region: &mut OverlapRegion,
    extz: &mut BitExtz,
    reads: &ReadStore,
    query: &[u8],
    maximum_errors: i64,
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    target_total_length: i64,
    extra_begin: i64,
    extra_end: i64,
    parameters: &HcAlignParameters,
) -> bool {
    let x_start = region.x_position_start as i64;
    let x_end = region.x_position_end as i64;
    let target_identifier = region.y_identifier;
    let target_strand = region.y_strand as u8;

    let mut window = WindowList {
        x_start: query_start as i32,
        x_end: query_end as i32,
        y_start: target_start as i32,
        y_end: target_end as i32,
        extra_begin: extra_begin as i16,
        extra_end: extra_end as i16,
        error: extz.result.error as i16,
        error_threshold: extz.result.threshold as i16,
        cigar_start: 0,
        cigar_length: 0,
    };

    // Fill the gap between the previous accepted window and this one by
    // extending forwards from the previous window's end.
    let mut covered_start = x_start;
    if let Some(last) = region.window_cigar.windows.last().copied() {
        let mut window_end = last.x_end as i64;
        let mut target_position = last.y_end as i64 + 1;
        while window_end < query_start - 1 && target_position < target_total_length {
            let window_start = window_end + 1;
            let (_, end) = window_id_by_start(
                region.x_position_start as i32,
                region.x_position_end as i32,
                window_start,
                parameters.window_length,
            );
            window_end = end;
            if aln_wlst_adv_exz(
                region,
                extz,
                reads,
                query,
                maximum_errors,
                target_strand,
                target_identifier,
                window_start,
                window_end,
                target_position,
                parameters,
                false,
            ) {
                target_position = region
                    .window_cigar
                    .windows
                    .last()
                    .expect("a window was just appended")
                    .y_end as i64
                    + 1;
            } else {
                break;
            }
        }
        covered_start = region
            .window_cigar
            .windows
            .last()
            .expect("the previous window is still there")
            .x_end as i64
            + 1;
    }

    // Fill the gap before this window by walking backwards from it.
    let windows_before = region.window_cigar.windows.len();
    let mut window_start = query_start;
    if window_start > covered_start {
        gen_backtrace_adv_exz(
            &mut window,
            region,
            reads,
            query,
            extz,
            target_strand,
            target_identifier,
        );
        let mut target_position = window.y_start as i64 - 1;
        while window_start > covered_start {
            let window_end = window_start - 1;
            let (_, start) = window_id_by_end(
                region.x_position_start as i32,
                window_end,
                parameters.window_length,
            );
            window_start = start;
            let y_start = target_position + 1 - (window_end + 1 - window_start);
            if y_start >= 0
                && aln_wlst_adv_exz(
                    region,
                    extz,
                    reads,
                    query,
                    maximum_errors,
                    target_strand,
                    target_identifier,
                    window_start,
                    window_end,
                    y_start,
                    parameters,
                    true,
                )
            {
                target_position = region
                    .window_cigar
                    .windows
                    .last()
                    .expect("a window was just appended")
                    .y_start as i64
                    - 1;
            } else {
                break;
            }
        }
    }

    region.alignment_length += (query_end + 1 - query_start) as u32;
    let overlap_length = x_end + 1 - x_start;
    let unaligned = (query_end + 1 - x_start) - region.alignment_length as i64;
    let aligned = overlap_length - unaligned;
    if !parameters.force_alignment
        && !pass_qovlp(overlap_length, aligned, parameters.overlap_cutoff)
    {
        region.window_cigar.windows.push(window);
        return false;
    }

    // The blocks filled backwards were appended in reverse chronological
    // order; the reference restores their order before appending this window.
    if region.window_cigar.windows.len() > windows_before {
        region.window_cigar.windows[windows_before..].reverse();
    }
    region.window_cigar.windows.push(window);
    true
}

/// `aln_wlst_adv_exz`: aligns one skipped block, optionally with its CIGAR.
#[allow(clippy::too_many_arguments)]
fn aln_wlst_adv_exz(
    region: &mut OverlapRegion,
    extz: &mut BitExtz,
    reads: &ReadStore,
    query: &[u8],
    maximum_errors: i64,
    reverse_strand: u8,
    target_identifier: u32,
    query_start: i64,
    query_end: i64,
    target_start: i64,
    parameters: &HcAlignParameters,
    is_cigar: bool,
) -> bool {
    let query_length = query_end + 1 - query_start;
    let threshold = double_error_threshold(
        initial_error_threshold(
            query_length,
            parameters.error_rate,
            parameters.window_length,
            maximum_errors,
        ),
        query_length,
    );
    let alignment_length = query_length + (threshold << 1);
    let target_total_length = reads.read_length(target_identifier as usize) as i64;

    let Some(target) =
        determine_overlap_region(threshold, target_start, target_total_length, alignment_length)
    else {
        return false;
    };
    if target.length + threshold < query_length {
        return false;
    }

    let target_sequence =
        reads.read_sub_region(target_identifier as usize, target.start, target.length, reverse_strand);
    let query_sequence = &query[query_start as usize..(query_start + query_length) as usize];

    if is_cigar {
        extz.clear_align();
        extz.semi_global_absent_diagonal_trace(
            &target_sequence,
            query_sequence,
            threshold,
            target.extra_begin as i64,
        );
    } else {
        extz.semi_global_absent_diagonal(
            &target_sequence,
            query_sequence,
            threshold,
            target.extra_begin as i64,
        );
        extz.result.pattern_start = 0;
    }

    if !extz.is_aligned() {
        return false;
    }

    let mut window = WindowList {
        x_start: query_start as i32,
        x_end: query_end as i32,
        y_start: (target.start + extz.result.pattern_start as i64) as i32,
        y_end: (target.start + extz.result.pattern_end as i64) as i32,
        error: extz.result.error as i16,
        error_threshold: threshold as i16,
        extra_begin: target.extra_begin as i16,
        extra_end: target.extra_end as i16,
        cigar_start: 0,
        cigar_length: 0,
    };
    let mut extra_begin = target.extra_begin as i64;
    let mut extra_end = target.extra_end as i64;

    if is_cigar {
        push_wcigar(&mut window, &mut region.window_cigar, extz);
        let at_boundary = (extz.result.pattern_end as i64 + 1 == target.length
            || extz.result.pattern_start == 0)
            && extz.result.error > 0;
        if at_boundary {
            if let Some((new_start, new_begin, new_end)) = recal_boundary_exz(
                query_sequence,
                reads,
                query_length,
                target.length,
                threshold,
                target.start,
                extz.result.pattern_start,
                extz.result.pattern_end,
                extz.result.error,
                target_identifier,
                reverse_strand,
                extz,
            ) {
                region.window_cigar.cells.truncate(window.cigar_start as usize);
                push_wcigar(&mut window, &mut region.window_cigar, extz);
                window.y_start = (new_start + extz.result.pattern_start as i64) as i32;
                window.y_end = (new_start + extz.result.pattern_end as i64) as i32;
                window.error = extz.result.error as i16;
                extra_begin = new_begin;
                extra_end = new_end;
            }
        }
    }

    window.extra_begin = extra_begin as i16;
    window.extra_end = extra_end as i16;

    region.alignment_length += query_length as u32;
    region.window_cigar.windows.push(window);
    true
}

/// `gen_backtrace_adv_exz`: re-runs the DP in trace mode from the seeded state
/// of a directly aligned window and stores its CIGAR.
///
/// The entered error and end position make the reference reuse the previous
/// result when it is already exact; otherwise the trace macro keeps the seeded
/// bounds and backtracks from them.
#[allow(clippy::too_many_arguments)]
fn gen_backtrace_adv_exz(
    window: &mut WindowList,
    region: &mut OverlapRegion,
    reads: &ReadStore,
    query: &[u8],
    extz: &mut BitExtz,
    reverse_strand: u8,
    target_identifier: u32,
) -> bool {
    if window.error < 0 || window.y_end < 0 {
        return false;
    }

    let query_start = window.x_start as i64;
    let query_end = window.x_end as i64;
    let query_length = query_end + 1 - query_start;
    let threshold = window.error_threshold as i64;
    let alignment_length = query_length + (threshold << 1);

    let mut target_start = window.y_start as i64;
    let mut extra_begin = window.extra_begin as i64;
    let mut extra_end = window.extra_end as i64;

    let target_length = if extra_end >= 0 {
        alignment_length - extra_begin - extra_end
    } else {
        let target_total_length = reads.read_length(target_identifier as usize) as i64;
        let mut length = target_start + alignment_length - extra_begin;
        if length > target_total_length {
            length = target_total_length;
        }
        length - target_start
    };

    let target_sequence =
        reads.read_sub_region(target_identifier as usize, target_start, target_length, reverse_strand);
    let query_sequence = &query[query_start as usize..(query_start + query_length) as usize];

    // `gen_backtrace_adv_exz` seeds the result and does not call `clear_align`.
    extz.result.text_start = 0;
    extz.result.text_end = (window.x_end - window.x_start) as i32;
    extz.result.text_length = query_length as i32;
    extz.result.pattern_start = -1;
    extz.result.pattern_end = (window.y_end - window.y_start) as i32;
    extz.result.pattern_length = target_length as i32;
    extz.result.error = window.error as i32;
    extz.result.threshold = window.error_threshold as i32;

    extz.semi_global_absent_diagonal_trace(
        &target_sequence,
        query_sequence,
        threshold,
        extra_begin,
    );

    if !extz.is_aligned() {
        window.error = -1;
        return false;
    }

    window.y_start = (target_start + extz.result.pattern_start as i64) as i32;
    window.y_end = (target_start + extz.result.pattern_end as i64) as i32;
    window.error = extz.result.error as i16;
    push_wcigar(window, &mut region.window_cigar, extz);

    let at_boundary = (extz.result.pattern_end as i64 + 1 == target_length
        || extz.result.pattern_start == 0)
        && extz.result.error > 0;
    if at_boundary {
        if let Some((new_start, new_begin, new_end)) = recal_boundary_exz(
            query_sequence,
            reads,
            query_length,
            target_length,
            threshold,
            target_start,
            extz.result.pattern_start,
            extz.result.pattern_end,
            extz.result.error,
            target_identifier,
            reverse_strand,
            extz,
        ) {
            region.window_cigar.cells.truncate(window.cigar_start as usize);
            push_wcigar(window, &mut region.window_cigar, extz);
            target_start = new_start;
            window.y_start = (target_start + extz.result.pattern_start as i64) as i32;
            window.y_end = (target_start + extz.result.pattern_end as i64) as i32;
            window.error = extz.result.error as i16;
            extra_begin = new_begin;
            extra_end = new_end;
        }
    }

    window.extra_begin = extra_begin as i16;
    window.extra_end = extra_end as i16;
    true
}

/// `recal_boundary_exz`: shifts a window that touched a boundary and accepts
/// the new alignment only when it has fewer errors.
///
/// Returns the new target start together with the new padding on success.
#[allow(clippy::too_many_arguments)]
fn recal_boundary_exz(
    query_sequence: &[u8],
    reads: &ReadStore,
    query_length: i64,
    pattern_length: i64,
    threshold: i64,
    target_offset: i64,
    pattern_start: i32,
    pattern_end: i32,
    previous_error: i32,
    target_identifier: u32,
    reverse_strand: u8,
    extz: &mut BitExtz,
) -> Option<(i64, i64, i64)> {
    let position = if pattern_start == 0 {
        target_offset
    } else if pattern_end as i64 + 1 == pattern_length {
        target_offset + pattern_end as i64 - query_length + 1
    } else {
        return None;
    };

    let alignment_length = query_length + (threshold << 1);
    let target_total_length = reads.read_length(target_identifier as usize) as i64;
    let target =
        determine_overlap_region(threshold, position, target_total_length, alignment_length)?;
    if target.start == target_offset && pattern_length == target.length {
        return None;
    }

    let target_sequence = reads.read_sub_region(
        target_identifier as usize,
        target.start,
        target.length,
        reverse_strand,
    );
    extz.clear_align();
    extz.semi_global_absent_diagonal_trace(
        &target_sequence,
        query_sequence,
        threshold,
        target.extra_begin as i64,
    );

    if extz.is_aligned() && extz.result.error < previous_error {
        Some((target.start, target.extra_begin as i64, target.extra_end as i64))
    } else {
        None
    }
}

/// `push_wcigar`: appends the packed CIGAR of `extz` to the shared cell buffer
/// and points the window at it.
fn push_wcigar(window: &mut WindowList, cigar: &mut WindowCigar, extz: &BitExtz) {
    window.cigar_start = cigar.cells.len() as u32;
    window.cigar_length = extz.cigar.packed.len() as u32;
    cigar.cells.extend_from_slice(&extz.cigar.packed);
}

/// `pass_qovlp`
fn pass_qovlp(overlap_length: i64, aligned_length: i64, ratio: f64) -> bool {
    aligned_length > 0 && (overlap_length as f64 * ratio) <= aligned_length as f64
}

/// `get_win_se_by_normalize_xs`: the block [start, end] that contains the
/// aligned origin, clipped to the overlap.
fn window_se_by_normalized_start(
    region: &OverlapRegion,
    normalized_start: i64,
    block_length: i64,
) -> (i64, i64) {
    let x_start = region.x_position_start as i64;
    let x_end = region.x_position_end as i64;
    let origin = (x_start / block_length) * block_length;
    let identifier = (normalized_start - origin) / block_length;
    let mut start = origin + identifier * block_length;
    if start < x_start {
        start = x_start;
    }
    let mut end = origin + (identifier + 1) * block_length - 1;
    if end > x_end {
        end = x_end;
    }
    (start, end)
}

/// `Reserve_Banded_BPM_Extension_REV`: the backward variant of
/// [`reserve_banded_bpm_extension`].
///
/// The reference reverses both strings, runs the forward extension and maps
/// the reported ends back with `length - 1 - end`. That is *not* the same
/// convention as C's `alignment_extension` wrapper, which reports
/// `length - end`; `ed_cut` uses this function directly.
fn reserve_banded_bpm_extension_reverse(
    pattern: &[u8],
    pattern_length: usize,
    text: &[u8],
    text_length: usize,
    error_threshold: usize,
) -> (u32, i32, i32) {
    let reversed_pattern: Vec<u8> = pattern[..pattern_length].iter().rev().copied().collect();
    let reversed_text: Vec<u8> = text[..text_length].iter().rev().copied().collect();
    let (error, pattern_end, text_end) = reserve_banded_bpm_extension(
        &reversed_pattern,
        pattern_length,
        &reversed_text,
        text_length,
        error_threshold,
    );
    if pattern_end == -1 || text_end == -1 {
        return (error, -1, -1);
    }
    (
        error,
        pattern_length as i32 - 1 - pattern_end,
        text_length as i32 - 1 - text_end,
    )
}

/// One accepted `ed_cut` alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CutAlignment {
    /// `r_err`
    error: i64,
    /// `q_end`, in query coordinates; only meaningful when the call succeeds.
    query_end: i64,
    /// `t_end`, in target coordinates; only meaningful when the call succeeds.
    target_end: i64,
    /// `aln_qlen`
    aligned_query_length: i64,
}

/// `ed_cut`: proves that a window can be aligned within the doubled budget and
/// reports how much of it was consumed.
#[allow(clippy::too_many_arguments)]
fn ed_cut(
    reads: &ReadStore,
    query: &[u8],
    target_identifier: usize,
    reverse_strand: u8,
    target_total_length: i64,
    query_start: i64,
    query_end: i64,
    target_start: i64,
    block_length: i64,
    error_rate: f64,
    maximum_errors: i64,
    reverse: bool,
) -> Option<CutAlignment> {
    let query_length = query_end + 1 - query_start;
    let threshold = double_error_threshold(
        initial_error_threshold(query_length, error_rate, block_length, maximum_errors),
        query_length,
    );
    let alignment_length = query_length + (threshold << 1);
    let target = determine_overlap_region(
        threshold,
        target_start,
        target_total_length,
        alignment_length,
    )?;
    let target_sequence = fill_subregion(
        reads,
        target_identifier,
        target.start,
        target.length,
        reverse_strand,
        target.extra_begin,
        target.extra_end,
    );
    let query_sequence = &query[query_start as usize..(query_start + query_length) as usize];
    let (error, pattern_end, text_end) = if reverse {
        reserve_banded_bpm_extension_reverse(
            &target_sequence,
            alignment_length as usize,
            query_sequence,
            query_length as usize,
            threshold as usize,
        )
    } else {
        reserve_banded_bpm_extension(
            &target_sequence,
            alignment_length as usize,
            query_sequence,
            query_length as usize,
            threshold as usize,
        )
    };

    let mut aligned_query_length = 0;
    if pattern_end != -1 && text_end != -1 {
        aligned_query_length = if reverse {
            query_length - text_end as i64
        } else {
            text_end as i64 + 1
        };
    }
    if aligned_query_length == 0 {
        return None;
    }
    Some(CutAlignment {
        error: error as i64,
        query_end: text_end as i64,
        target_end: pattern_end as i64,
        aligned_query_length,
    })
}

/// Re-runs the traceback of the window at `index` in place, as the reference
/// does when a gap fill needs the CIGAR of a directly aligned neighbour.
pub(crate) fn backtrace_window_in_place(
    region: &mut OverlapRegion,
    index: i64,
    reads: &ReadStore,
    query: &[u8],
    extz: &mut BitExtz,
    reverse_strand: u8,
) {
    let target_identifier = region.y_identifier;
    let mut window = region.window_cigar.windows[index as usize].clone();
    gen_backtrace_adv_exz(
        &mut window,
        region,
        reads,
        query,
        extz,
        reverse_strand,
        target_identifier,
    );
    region.window_cigar.windows[index as usize] = window;
}

/// `gen_extend_err_0_exz`: the errors in one block that no window covered.
///
/// Both neighbours of the block are used when they touch it, the left one
/// anchored at its end and the right one at `y_start - query_length`, so a
/// block is extended from whichever side the chain makes available.
#[allow(clippy::too_many_arguments)]
pub fn gen_extend_err_0_exz(
    region: &mut OverlapRegion,
    reads: &ReadStore,
    query: &[u8],
    extz: &mut BitExtz,
    block_length: i64,
    error_rate: f64,
    maximum_errors: i64,
    query_start: i64,
    query_end: i64,
    previous_window: i64,
) -> i64 {
    let query_length = query_end + 1 - query_start;
    let target_identifier = region.y_identifier as usize;
    let reverse_strand = region.y_strand as u8;
    let target_total_length = reads.read_length(target_identifier) as i64;
    let mut estimated_target =
        (query_start - region.x_position_start as i64) + region.y_position_start as i64;
    estimated_target += region
        .fake_cigar
        .y_start_offset(query_start)
        .expect("the fake cigar must cover every block start");

    let mut error = [0i64; 2];
    let mut aligned = [0i64; 2];
    let mut boundary = [-1i64; 2];

    if previous_window > 0 {
        let window = region.window_cigar.windows[previous_window as usize];
        if query_start == window.x_end as i64 + 1 {
            if window.cigar_length == 0 {
                backtrace_window_in_place(
                    region,
                    previous_window,
                    reads,
                    query,
                    extz,
                    reverse_strand,
                );
            }
            boundary[0] = region.window_cigar.windows[previous_window as usize].y_end as i64 + 1;
        }
    }

    let next_window = previous_window + 1;
    if next_window < region.window_cigar.windows.len() as i64 {
        let window = region.window_cigar.windows[next_window as usize];
        if query_end + 1 == window.x_start as i64 {
            if window.cigar_length == 0 {
                backtrace_window_in_place(region, next_window, reads, query, extz, reverse_strand);
            }
            boundary[1] = region.window_cigar.windows[next_window as usize].y_start as i64
                - query_length;
        }
    }

    if boundary[0] == -1 && boundary[1] == -1 {
        boundary[0] = estimated_target;
        boundary[1] = estimated_target;
    } else if boundary[0] == -1 && boundary[1] != -1 {
        boundary[0] = boundary[1];
    } else if boundary[1] == -1 && boundary[0] != -1 {
        boundary[1] = boundary[0];
    }

    if boundary[0] != -1 {
        if let Some(cut) = ed_cut(
            reads,
            query,
            target_identifier,
            reverse_strand,
            target_total_length,
            query_start,
            query_end,
            boundary[0],
            block_length,
            error_rate,
            maximum_errors,
            false,
        ) {
            error[0] = cut.error;
            aligned[0] = cut.aligned_query_length;
        } else {
            error[0] = query_length;
            aligned[0] = 0;
        }
    }

    if boundary[1] != -1 {
        if let Some(cut) = ed_cut(
            reads,
            query,
            target_identifier,
            reverse_strand,
            target_total_length,
            query_start,
            query_end,
            boundary[1],
            block_length,
            error_rate,
            maximum_errors,
            true,
        ) {
            error[1] = cut.error;
            aligned[1] = cut.aligned_query_length;
        } else {
            error[1] = query_length;
            aligned[1] = 0;
        }
    }

    if aligned[0] != 0 && aligned[1] != 0 {
        if aligned[0] + aligned[1] <= query_length {
            error[0] + error[1] + query_length - (aligned[0] + aligned[1])
        } else {
            let ratio = query_length as f64 / (aligned[0] + aligned[1]) as f64;
            ((error[0] + error[1]) as f64 * ratio) as i64
        }
    } else if aligned[0] == 0 && aligned[1] == 0 {
        query_length
    } else if aligned[0] != 0 {
        error[0] + (query_length - aligned[0])
    } else {
        error[1] + (query_length - aligned[1])
    }
}

/// Result of [`gen_extend_err_exz`]: the mean error rate over the whole
/// overlap and the error count, `rr` and `re` in the reference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExtensionError {
    pub error_rate: f64,
    pub error_bases: i64,
}

/// `gen_extend_err_exz`: sums the alignment errors of one overlap, including
/// the blocks between the windows that were never aligned, and returns the
/// mean error rate.
///
/// The caller aborts at `rr > error_rate`, so a partially walked overlap can
/// still exceed `maximum_error_rate` and come back as `f64::MAX`.
#[allow(clippy::too_many_arguments)]
pub fn gen_extend_err_exz(
    region: &mut OverlapRegion,
    reads: &ReadStore,
    query: &[u8],
    extz: &mut BitExtz,
    block_length: i64,
    overlap_cutoff: f64,
    error_rate: f64,
    maximum_error_rate: f64,
    maximum_errors: i64,
    secondary_check: bool,
) -> ExtensionError {
    let overlap_length = region.x_position_end as i64 + 1 - region.x_position_start as i64;
    let mut result = ExtensionError {
        error_rate: f64::MAX,
        error_bases: i64::MAX,
    };
    if secondary_check
        && !similarity_pass(
            overlap_length,
            region.alignment_length as i64,
            false,
            overlap_cutoff,
            Some(error_rate),
        )
    {
        return result;
    }

    let mut total_length = 0i64;
    let mut total_error = 0i64;
    let mut current_end = region.x_position_end as i64;
    let mut window_identifier = region.window_cigar.windows.len() as i64 - 1;
    while window_identifier >= 0 {
        let window = region.window_cigar.windows[window_identifier as usize];
        total_length += window.x_end as i64 + 1 - window.x_start as i64;
        total_error += window.error as i64;

        let mut window_end = window.x_end as i64;
        while window_end < current_end {
            let window_start = window_end + 1;
            let (_, end) = window_id_by_start(
                region.x_position_start as i32,
                region.x_position_end as i32,
                window_start,
                block_length,
            );
            window_end = end;
            let length = window_end + 1 - window_start;
            total_length += length;
            total_error += gen_extend_err_0_exz(
                region,
                reads,
                query,
                extz,
                block_length,
                error_rate,
                maximum_errors,
                window_start,
                window_end,
                window_identifier,
            );
            if maximum_error_rate > 0.0
                && total_error as f64 > overlap_length as f64 * maximum_error_rate
            {
                return result;
            }
        }

        current_end = region.window_cigar.windows[window_identifier as usize].x_start as i64 - 1;
        if maximum_error_rate > 0.0
            && total_error as f64 > overlap_length as f64 * maximum_error_rate
        {
            return result;
        }
        window_identifier -= 1;
    }

    if current_end >= region.x_position_start as i64 {
        let mut window_end = region.x_position_start as i64 - 1;
        while window_end < current_end {
            let window_start = window_end + 1;
            let (_, end) = window_id_by_start(
                region.x_position_start as i32,
                region.x_position_end as i32,
                window_start,
                block_length,
            );
            window_end = end;
            let length = window_end + 1 - window_start;
            total_length += length;
            total_error += gen_extend_err_0_exz(
                region,
                reads,
                query,
                extz,
                block_length,
                error_rate,
                maximum_errors,
                window_start,
                window_end,
                window_identifier,
            );
            if maximum_error_rate > 0.0
                && total_error as f64 > overlap_length as f64 * maximum_error_rate
            {
                return result;
            }
        }
    }

    debug_assert_eq!(total_length, overlap_length);
    result.error_bases = total_error;
    result.error_rate = total_error as f64 / total_length as f64;
    result
}

/// The parameters `init_waln` derives from a window: kept here because the
/// window mapping is shared with [`determine_overlap_region`].
pub type WindowRegion = TargetRegion;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::overlap::region::FakeCigar;

    /// The linear congruential generator the golden probe uses, so the
    /// sequences match the reference byte for byte.
    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> u32 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 33) as u32
        }
    }

    fn random_sequence(state: &mut Lcg, length: usize) -> Vec<u8> {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        (0..length)
            .map(|_| BASES[(state.next() & 3) as usize])
            .collect()
    }

    fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect()
    }

    fn derive_target(
        left: &[u8],
        core: &[u8],
        right: &[u8],
        substitution_step: usize,
        insertion: Option<(usize, &[u8])>,
    ) -> Vec<u8> {
        let mut target = Vec::new();
        target.extend_from_slice(left);
        target.extend_from_slice(core);
        target.extend_from_slice(right);
        if substitution_step > 0 {
            for offset in (0..core.len()).step_by(substitution_step) {
                let base = &mut target[left.len() + offset];
                *base = if *base == b'A' { b'C' } else { b'A' };
            }
        }
        if let Some((site, bases)) = insertion {
            target.splice(left.len() + site..left.len() + site, bases.iter().copied());
        }
        target
    }

    fn mutate(target: &mut [u8], offset: usize, range: std::ops::Range<usize>) {
        for index in range {
            let base = &mut target[offset + index];
            *base = if *base == b'A' { b'C' } else { b'A' };
        }
    }

    /// One window of the golden output.
    #[derive(Debug, PartialEq, Eq)]
    struct WindowExpectation {
        x_start: i32,
        x_end: i32,
        y_start: i32,
        y_end: i32,
        extra_begin: i16,
        extra_end: i16,
        error: i16,
        threshold: i16,
        cells: Vec<u16>,
    }

    /// The two runs of [`gen_extend_err_exz`] the golden probe performs: the
    /// caller's pass with the early abort, and a second pass with the abort
    /// disabled so the whole-overlap arithmetic is checked as well.
    #[derive(Debug, Clone, Copy, PartialEq)]
    struct ExtensionExpectation {
        early: ExtensionError,
        full: ExtensionError,
    }

    struct Case {
        accepted: bool,
        alignment_length: u32,
        windows: Vec<WindowExpectation>,
        extension: Option<ExtensionExpectation>,
        /// The window list after the extension passes; a gap fill may backtrace
        /// a direct window and rewrite its coordinates and CIGAR.
        final_windows: Vec<WindowExpectation>,
    }

    fn snapshot_windows(region: &OverlapRegion) -> Vec<WindowExpectation> {
        region
            .window_cigar
            .windows
            .iter()
            .map(|window| {
                let cells = (0..window.cigar_length)
                    .map(|offset| (window.cigar_start + offset) as usize)
                    .map(|index| region.window_cigar.cells[index])
                    .collect();
                WindowExpectation {
                    x_start: window.x_start,
                    x_end: window.x_end,
                    y_start: window.y_start,
                    y_end: window.y_end,
                    extra_begin: window.extra_begin,
                    extra_end: window.extra_end,
                    error: window.error,
                    threshold: window.error_threshold,
                    cells,
                }
            })
            .collect()
    }

    #[allow(clippy::too_many_arguments)]
    fn run_case(
        reads: &[Vec<u8>],
        query_identifier: u32,
        target_identifier: u32,
        x_start: u32,
        x_end: u32,
        y_start: u32,
        y_end: u32,
        target_strand: u32,
        gaps: &[(u32, i32)],
        error_rate: f64,
        window_length: i64,
        overlap_cutoff: f64,
        force_alignment: bool,
    ) -> Case {
        let mut store = ReadStore::default();
        for sequence in reads {
            store.push_read("read", sequence);
        }

        let mut fake_cigar = FakeCigar::default();
        for (site, shift) in gaps {
            fake_cigar.push(*site, *shift);
        }

        let mut region = OverlapRegion {
            x_identifier: query_identifier,
            x_position_start: x_start,
            x_position_end: x_end,
            y_identifier: target_identifier,
            y_position_start: y_start,
            y_position_end: y_end,
            y_strand: target_strand,
            fake_cigar,
            ..OverlapRegion::default()
        };

        let mut extz = BitExtz::new();
        let parameters = HcAlignParameters {
            error_rate,
            window_length,
            overlap_cutoff,
            force_alignment,
        };
        let accepted =
            align_hc_ed_post_extz(&mut region, &mut extz, &store, &reads[query_identifier as usize], &parameters);

        let windows = snapshot_windows(&region);

        // The reference only runs the extension pass on an accepted overlap,
        // but the gap arithmetic is worth checking on rejected ones too (they
        // still carry windows). The probe calls it whenever either holds.
        let extension = if accepted || !windows.is_empty() {
            let maximum_error_rate = error_rate * 1.5 + 0.000001;
            let early = gen_extend_err_exz(
                &mut region,
                &store,
                &reads[query_identifier as usize],
                &mut extz,
                window_length,
                -1.0,
                error_rate,
                maximum_error_rate,
                THRESHOLD_MAXIMUM,
                false,
            );
            let full = gen_extend_err_exz(
                &mut region,
                &store,
                &reads[query_identifier as usize],
                &mut extz,
                window_length,
                -1.0,
                error_rate,
                0.0,
                THRESHOLD_MAXIMUM,
                false,
            );
            Some(ExtensionExpectation { early, full })
        } else {
            None
        };
        let final_windows = snapshot_windows(&region);

        Case {
            accepted,
            alignment_length: region.alignment_length,
            windows,
            extension,
            final_windows,
        }
    }

    fn window(
        x_start: i32,
        x_end: i32,
        y_start: i32,
        y_end: i32,
        extra_begin: i16,
        extra_end: i16,
        error: i16,
        threshold: i16,
        cells: Vec<u16>,
    ) -> WindowExpectation {
        WindowExpectation {
            x_start,
            x_end,
            y_start,
            y_end,
            extra_begin,
            extra_end,
            error,
            threshold,
            cells,
        }
    }

    fn expect_extension(case: &Case, early: (f64, i64), full: (f64, i64)) {
        assert_eq!(
            case.extension,
            Some(ExtensionExpectation {
                early: ExtensionError {
                    error_rate: early.0,
                    error_bases: early.1,
                },
                full: ExtensionError {
                    error_rate: full.0,
                    error_bases: full.1,
                },
            })
        );
    }

    /// Golden values printed by the real `align_hc_ed_post_extz` linked
    /// against `Correct.o`; the probe lives in `/tmp/ec_align_probe.cpp`.
    #[test]
    fn matches_the_reference_window_records() {
        let mut state = Lcg(123_456_789);
        let query = random_sequence(&mut state, 2000);
        let prefix = random_sequence(&mut state, 600);
        let suffix = random_sequence(&mut state, 300);

        // A: exact, full two-window overlap.
        {
            let target = derive_target(&prefix, &query, &suffix, 0, None);
            let case = run_case(
                &[query.clone(), target],
                0,
                1,
                200,
                1799,
                800,
                2399,
                0,
                &[(200, 0), (1799, 0)],
                0.04,
                775,
                0.9,
                false,
            );
            assert!(case.accepted);
            assert_eq!(case.alignment_length, 1600);
            assert_eq!(
                case.windows,
                vec![
                    window(200, 774, 777, 1374, 0, 0, 0, 23, vec![]),
                    window(775, 1549, 1344, 2149, 0, 0, 0, 31, vec![]),
                    window(1550, 1799, 2140, 2399, 0, 0, 0, 10, vec![]),
                ]
            );
            expect_extension(&case, (0.0, 0), (0.0, 0));
            assert_eq!(case.final_windows, case.windows);
        }

        // B: one substitution every 100 bases of the target.
        {
            let target = derive_target(&prefix, &query, &suffix, 100, None);
            let case = run_case(
                &[query.clone(), target],
                0,
                1,
                200,
                1799,
                800,
                2399,
                0,
                &[(200, 0), (1799, 0)],
                0.04,
                775,
                0.9,
                false,
            );
            assert!(case.accepted);
            assert_eq!(case.alignment_length, 1600);
            assert_eq!(
                case.windows,
                vec![
                    window(200, 774, 777, 1374, 0, 0, 6, 23, vec![]),
                    window(775, 1549, 1344, 2149, 0, 0, 8, 31, vec![]),
                    window(1550, 1799, 2140, 2399, 0, 0, 2, 10, vec![]),
                ]
            );
            expect_extension(&case, (16.0 / 1600.0, 16), (16.0 / 1600.0, 16));
            assert_eq!(case.final_windows, case.windows);
        }

        // C: a three base insertion in the target.
        {
            let target = derive_target(&prefix, &query, &suffix, 0, Some((900, b"ACG")));
            let case = run_case(
                &[query.clone(), target],
                0,
                1,
                200,
                1799,
                800,
                2402,
                0,
                &[(200, 0), (1799, 0)],
                0.04,
                775,
                0.9,
                false,
            );
            assert!(case.accepted);
            assert_eq!(
                case.windows,
                vec![
                    window(200, 774, 777, 1374, 0, 0, 0, 23, vec![]),
                    window(775, 1549, 1344, 2152, 0, 0, 3, 31, vec![]),
                    window(1550, 1799, 2140, 2402, 0, 0, 0, 10, vec![]),
                ]
            );
            expect_extension(&case, (3.0 / 1600.0, 3), (3.0 / 1600.0, 3));
            assert_eq!(case.final_windows, case.windows);
        }

        // D: the target starts ten bases before the overlap, so the first
        // window is padded on the left.
        {
            let mut target = random_sequence(&mut state, 10);
            target.extend_from_slice(&query[200..]);
            let case = run_case(
                &[query.clone(), target],
                0,
                1,
                200,
                1799,
                10,
                1809,
                0,
                &[(200, 0), (1799, 0)],
                0.04,
                775,
                0.9,
                false,
            );
            assert!(case.accepted);
            assert_eq!(
                case.windows,
                vec![
                    window(200, 774, 0, 584, 13, 0, 0, 23, vec![]),
                    window(775, 1549, 554, 1359, 0, 0, 0, 31, vec![]),
                    window(1550, 1799, 1350, 1609, 0, 0, 0, 10, vec![]),
                ]
            );
            expect_extension(&case, (0.0, 0), (0.0, 0));
            assert_eq!(case.final_windows, case.windows);
        }

        // E: the overlap ends at the target end, so the last window is padded
        // on the right.
        {
            let mut target = prefix.clone();
            target.extend_from_slice(&query[200..900]);
            let case = run_case(
                &[query.clone(), target],
                0,
                1,
                200,
                899,
                600,
                1299,
                0,
                &[(200, 0), (899, 0)],
                0.04,
                775,
                0.9,
                false,
            );
            assert!(case.accepted);
            assert_eq!(case.alignment_length, 700);
            assert_eq!(
                case.windows,
                vec![
                    window(200, 774, 577, 1174, 0, 0, 0, 23, vec![]),
                    window(775, 899, 1170, 1299, 0, 5, 0, 5, vec![]),
                ]
            );
            expect_extension(&case, (0.0, 0), (0.0, 0));
            assert_eq!(case.final_windows, case.windows);
        }

        // F: the target overlaps the reverse strand.
        {
            let target = reverse_complement(&derive_target(&prefix, &query, &suffix, 0, None));
            let case = run_case(
                &[query.clone(), target],
                0,
                1,
                200,
                1799,
                800,
                2399,
                1,
                &[(200, 0), (1799, 0)],
                0.04,
                775,
                0.9,
                false,
            );
            assert!(case.accepted);
            assert_eq!(case.alignment_length, 1600);
            assert_eq!(
                case.windows,
                vec![
                    window(200, 774, 777, 1374, 0, 0, 0, 23, vec![]),
                    window(775, 1549, 1344, 2149, 0, 0, 0, 31, vec![]),
                    window(1550, 1799, 2140, 2399, 0, 0, 0, 10, vec![]),
                ]
            );
            expect_extension(&case, (0.0, 0), (0.0, 0));
            assert_eq!(case.final_windows, case.windows);
        }

        // G: a 60 base deletion inside every window makes the overlap fail.
        {
            let mut target = derive_target(&prefix, &query, &suffix, 0, None);
            target.drain(1100..1160);
            let case = run_case(
                &[query.clone(), target],
                0,
                1,
                200,
                1799,
                800,
                2339,
                0,
                &[(200, 0), (1799, 0)],
                0.04,
                775,
                0.9,
                false,
            );
            assert!(!case.accepted);
            assert_eq!(case.alignment_length, 0);
            assert!(case.windows.is_empty());
            assert_eq!(case.extension, None);
        }

        // H: heavily diverged, rejected by the final fraction check.
        {
            let target = derive_target(&prefix, &query, &suffix, 8, None);
            let case = run_case(
                &[query.clone(), target],
                0,
                1,
                200,
                1799,
                800,
                2399,
                0,
                &[(200, 0), (1799, 0)],
                0.04,
                775,
                0.9,
                false,
            );
            assert!(!case.accepted);
            assert!(case.windows.is_empty());
            assert_eq!(case.extension, None);
        }

        // I: the same overlap is accepted when the fraction check is forced,
        // but no window passes its edit budget.
        {
            let target = derive_target(&prefix, &query, &suffix, 8, None);
            let case = run_case(
                &[query.clone(), target],
                0,
                1,
                200,
                1799,
                800,
                2399,
                0,
                &[(200, 0), (1799, 0)],
                0.04,
                775,
                0.9,
                true,
            );
            assert!(case.accepted);
            assert!(case.windows.is_empty());
            expect_extension(&case, (f64::MAX, i64::MAX), (500.0 / 1600.0, 500));
            assert_eq!(case.final_windows, case.windows);
        }

        // K: the first window is too diverged, so the second one is reached
        // across a gap; its CIGAR comes from the backtrace.
        {
            let mut target = derive_target(&prefix, &query, &suffix, 0, None);
            mutate(&mut target, 600, 200..600);
            mutate(&mut target, 600, 1000..1005);
            let case = run_case(
                &[query.clone(), target],
                0,
                1,
                200,
                1799,
                800,
                2399,
                0,
                &[(200, 0), (1799, 0)],
                0.04,
                775,
                0.9,
                false,
            );
            assert!(!case.accepted);
            assert_eq!(case.alignment_length, 775);
            assert_eq!(
                case.windows,
                vec![window(
                    775,
                    1549,
                    1375,
                    2149,
                    0,
                    0,
                    5,
                    31,
                    vec![(0 << 14) | 225, (1 << 14) | 3, (3 << 14) | 1, 2, (2 << 14) | 1, (0 << 14) | 544],
                )]
            );
            expect_extension(&case, (f64::MAX, i64::MAX), (350.0 / 1600.0, 350));
            assert_eq!(case.final_windows, case.windows);
        }

        // L: the middle window is too diverged; the third window is reached
        // backwards and its CIGAR is an exact run.
        {
            let mut target = derive_target(&prefix, &query, &suffix, 0, None);
            mutate(&mut target, 600, 975..1550);
            let case = run_case(
                &[query.clone(), target],
                0,
                1,
                200,
                1799,
                800,
                2399,
                0,
                &[(200, 0), (1799, 0)],
                0.04,
                775,
                0.9,
                false,
            );
            assert!(!case.accepted);
            assert_eq!(case.alignment_length, 825);
            assert_eq!(
                case.windows,
                vec![
                    window(200, 774, 777, 1374, 0, 0, 0, 23, vec![]),
                    window(1550, 1799, 2150, 2399, 0, 0, 0, 10, vec![(0 << 14) | 250]),
                ]
            );
            expect_extension(&case, (f64::MAX, i64::MAX), (527.0 / 1600.0, 527));
            assert_eq!(case.final_windows, case.windows);
        }

        // M: the last window is too diverged, so the gap ahead of the second
        // window is filled by backtracing it; the backtrace rewrites its
        // y_start from the recovered region start to the aligned start and
        // gives it the exact-match CIGAR.
        {
            let mut target = derive_target(&prefix, &query, &suffix, 0, None);
            mutate(&mut target, 600, 1600..1799);
            let case = run_case(
                &[query.clone(), target],
                0,
                1,
                200,
                1799,
                800,
                2399,
                0,
                &[(200, 0), (1799, 0)],
                0.04,
                775,
                0.9,
                false,
            );
            assert!(!case.accepted);
            assert_eq!(case.alignment_length, 1350);
            assert_eq!(
                case.windows,
                vec![
                    window(200, 774, 777, 1374, 0, 0, 0, 23, vec![]),
                    window(775, 1549, 1344, 2149, 0, 0, 0, 31, vec![]),
                ]
            );
            expect_extension(&case, (f64::MAX, i64::MAX), (163.0 / 1600.0, 163));
            assert_eq!(
                case.final_windows,
                vec![
                    window(200, 774, 777, 1374, 0, 0, 0, 23, vec![]),
                    window(775, 1549, 1375, 2149, 0, 0, 0, 31, vec![(0 << 14) | 775]),
                ]
            );
        }

        // J: a gap shift in the fake CIGAR moves the third window.
        {
            let target = derive_target(&prefix, &query, &suffix, 0, Some((900, b"ACG")));
            let case = run_case(
                &[query.clone(), target],
                0,
                1,
                200,
                1799,
                800,
                2402,
                0,
                &[(200, 0), (1100, 3), (1799, 3)],
                0.04,
                775,
                0.9,
                false,
            );
            assert!(case.accepted);
            assert_eq!(
                case.windows,
                vec![
                    window(200, 774, 777, 1374, 0, 0, 0, 23, vec![]),
                    window(775, 1549, 1344, 2152, 0, 0, 3, 31, vec![]),
                    window(1550, 1799, 2143, 2402, 0, 0, 0, 10, vec![]),
                ]
            );
            expect_extension(&case, (3.0 / 1600.0, 3), (3.0 / 1600.0, 3));
            assert_eq!(case.final_windows, case.windows);
        }
    }
}

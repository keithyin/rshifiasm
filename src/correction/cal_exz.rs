//! The `cal_exz_*` alignment layer of the new HiFi correction path, the
//! counterpart of the dispatch wrappers, `cal_exact_exz`, `cal_exz_infi_adv`,
//! `cal_estimate_err*` and the window bookkeeping helpers in `Correct.cpp`.
//!
//! The reference passes `uref`/`hpc_g` for the unitig and Hi-C paths and a
//! scratch `UC_Read` for the target; the HiFi correction path this port
//! follows always passes NULL for the first two and uses `rref`, so those
//! parameters are gone and the target always comes from
//! [`ReadStore::read_sub_region`]. The prefix reuse inside
//! `retrieve_str_seq_exz` only keeps `tu->seq` alive between calls; the bytes
//! it exposes are exactly `[start, start + length)`, so always retrieving the
//! full interval is content-equivalent and the scratch buffer disappears.
//!
//! Every clamp below is applied in the order the reference applies it,
//! because the band a call selects depends on the threshold after clamping,
//! and the band is part of the output.

use crate::alignment::edit_distance::{
    band_word_count, extension_backward_64, extension_backward_no_trace,
    extension_backward_trace, extension_forward_64, extension_forward_no_trace,
    extension_forward_trace, global_64, global_no_trace, global_trace,
    semi_global_64_absent_diagonal, semi_global_absent_diagonal_no_trace,
    semi_global_absent_diagonal_trace,
};
use crate::correction::hc_align::BitExtz;
use crate::correction::region::THRESHOLD_MAXIMUM;
use crate::correction::target::TargetSource;
use crate::correction::window::{window_id_by_start, WindowCigar};
use crate::overlap::region::OverlapRegion;

/// Largest single-pass error the new correction path allows, `MAX_SIN_E`.
pub const MAXIMUM_SINGLE_ERROR: i64 = 2047;
/// Longest span aligned in one pass, `MAX_SIN_L`.
pub const MAXIMUM_SINGLE_LENGTH: i64 = 10000;
/// Spans up to this length are forced through every threshold, `FORCE_SIN_L`.
pub const FORCED_SINGLE_LENGTH: i64 = 512;

/// `mode` of `cal_exact_exz`/`cal_exz_infi_adv`: global alignment.
pub const MODE_GLOBAL: i64 = 0;
/// Forward extension from a seeded target start.
pub const MODE_FORWARD_EXTENSION: i64 = 1;
/// Backward extension from a seeded target end.
pub const MODE_BACKWARD_EXTENSION: i64 = 2;
/// Semi-global alignment against the mapped target window.
pub const MODE_SEMI_GLOBAL: i64 = 3;

/// The result of `init_waln`: the clipped target interval plus the number of
/// bases the reference skipped before it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WindowAlignment {
    /// First target base of the interval, after clipping.
    pub start: i64,
    /// Number of target bases in the interval.
    pub length: i64,
    /// Bases before `start` that the mapped start pointed at, `aux_beg`.
    pub auxiliary_begin: i64,
    /// Unused target bases after the interval, `aux_end`.
    ///
    /// The reference computes it before `aux_beg` shortens the interval, so it
    /// can differ from `alignment_length - length`.
    pub auxiliary_end: i64,
}

/// `pts`/`pte`/`pthre` as threaded through the `cal_exz_infi_adv` calls of one
/// alignment attempt: the cached target boundaries and the highest threshold
/// already tried.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExzAdvState {
    pub pattern_start: i64,
    pub pattern_end: i64,
    pub threshold: i64,
}

impl Default for ExzAdvState {
    fn default() -> Self {
        Self {
            pattern_start: -1,
            pattern_end: -1,
            threshold: -1,
        }
    }
}

/// `scale_ed_thre`: rounds the band width an error count needs up to a whole
/// number of 64-bit words, converts the rounded width back to the error it can
/// hold, and caps it at `maximum_error`.
///
/// The `err << 1` is 32-bit arithmetic in the reference and is converted to
/// 64 bits only afterwards, so the shift and the addition wrap here as well.
pub fn scale_ed_thre(error: u32, maximum_error: u32) -> u64 {
    let band_width = u64::from(error.wrapping_shl(1).wrapping_add(1));
    let mut rounded = band_width >> 6;
    rounded <<= 6;
    if rounded < band_width {
        rounded += 64;
    }
    let mut scaled = ((rounded - 1) >> 1) as u32;
    if scaled > maximum_error {
        scaled = maximum_error;
    }
    u64::from(scaled)
}

/// `cal_exz_global`: global alignment with traceback.
///
/// The width dispatch of the reference (`ed_band_cal_global_64/128/infi`) is
/// inside [`global_trace`], which is width-generic.
pub fn cal_exz_global(extz: &mut BitExtz, pattern: &[u8], text: &[u8], threshold: i32) {
    global_trace(
        pattern,
        pattern.len() as i32,
        text,
        text.len() as i32,
        threshold,
        &mut extz.result,
        &mut extz.cigar,
    );
}

/// `cal_exz_extension_0`: forward extension from the seeded target start.
pub fn cal_exz_extension_0(extz: &mut BitExtz, pattern: &[u8], text: &[u8], threshold: i32) {
    extension_forward_trace(
        pattern,
        pattern.len() as i32,
        text,
        text.len() as i32,
        threshold,
        &mut extz.result,
        &mut extz.cigar,
    );
}

/// `cal_exz_extension_1`: backward extension from the seeded target end.
pub fn cal_exz_extension_1(extz: &mut BitExtz, pattern: &[u8], text: &[u8], threshold: i32) {
    extension_backward_trace(
        pattern,
        pattern.len() as i32,
        text,
        text.len() as i32,
        threshold,
        &mut extz.result,
        &mut extz.cigar,
    );
}

/// `cal_exz_semi`: semi-global alignment on an absent anti-diagonal.
pub fn cal_exz_semi(
    extz: &mut BitExtz,
    pattern: &[u8],
    text: &[u8],
    threshold: i32,
    auxiliary_begin: i32,
) {
    semi_global_absent_diagonal_trace(
        pattern,
        pattern.len() as i32,
        text,
        text.len() as i32,
        threshold,
        auxiliary_begin,
        &mut extz.result,
        &mut extz.cigar,
    );
}

/// `init_waln`: clips the mapped target interval to the target, padding the
/// right side with the unaligned remainder.
///
/// Returns `None` when the mapped start lies outside the target or when the
/// target cannot supply enough bases for the requested alignment length.
pub fn initialise_window_alignment(
    threshold: i64,
    start: i64,
    target_total_length: i64,
    alignment_length: i64,
) -> Option<WindowAlignment> {
    if start < 0
        || start >= target_total_length
        || (target_total_length - start + (2 * threshold) + THRESHOLD_MAXIMUM) < alignment_length
    {
        return None;
    }

    let mut auxiliary_begin = 0;
    let mut result_start = start - threshold;
    let mut result_length = target_total_length - result_start;
    if result_length > alignment_length {
        result_length = alignment_length;
    }
    let auxiliary_end = alignment_length - result_length;
    if result_start < 0 {
        auxiliary_begin = -result_start;
        result_start = 0;
        result_length -= auxiliary_begin;
    }

    Some(WindowAlignment {
        start: result_start,
        length: result_length,
        auxiliary_begin,
        auxiliary_end,
    })
}

/// `update_semi_coord`: maps `[query_start, query_end)` through the fake cigar
/// and clips it against the target read.
///
/// Returns `(target_start, target_end, auxiliary_begin)`; every value is `-1`
/// when the mapping fails, which is how the caller learns that the
/// semi-global alignment is impossible.
pub fn update_semi_coord(
    region: &OverlapRegion,
    targets: &TargetSource,
    query_start: i64,
    query_end: i64,
    threshold: i64,
) -> (i64, i64, i64) {
    let query_length = query_end - query_start;
    let mut start = (query_start - i64::from(region.x_position_start))
        + i64::from(region.y_position_start);
    start += region
        .fake_cigar
        .y_start_offset(query_start)
        .expect("the fake cigar must cover the query start");
    let alignment_length = query_length + (threshold << 1);
    let target_total_length = targets.length(region.y_identifier as usize);

    match initialise_window_alignment(threshold, start, target_total_length, alignment_length) {
        Some(window) => (
            window.start,
            window.start + window.length,
            window.auxiliary_begin,
        ),
        None => (-1, -1, -1),
    }
}

/// `adjust_ext_offset`: shrinks a forward or backward extension to the target
/// bases that exist, then clamps the four coordinates.
pub fn adjust_ext_offset(
    query_start: &mut i64,
    query_end: &mut i64,
    target_start: &mut i64,
    target_end: &mut i64,
    query_length: i64,
    target_length: i64,
    threshold: i64,
    mode: i64,
) {
    if mode == MODE_FORWARD_EXTENSION {
        let query_offset = query_length - *query_start;
        let target_offset = target_length - *target_start;
        if query_offset <= target_offset {
            *query_end = query_length;
            *target_end = *target_start + query_offset + threshold;
        } else {
            *target_end = target_length;
            *query_end = *query_start + target_offset + threshold;
        }
    } else if mode == MODE_BACKWARD_EXTENSION {
        let query_offset = *query_end;
        let target_offset = *target_end;
        if query_offset <= target_offset {
            *query_start = 0;
            *target_start = *target_end - query_offset - threshold;
        } else {
            *target_start = 0;
            *query_start = *query_end - target_offset - threshold;
        }
    }
    if *query_start < 0 {
        *query_start = 0;
    }
    if *target_start < 0 {
        *target_start = 0;
    }
    if *query_end > query_length {
        *query_end = query_length;
    }
    if *target_end > target_length {
        *target_end = target_length;
    }
}

/// `cal_exact_exz`: checks whether the mapped target interval is identical to
/// the query interval without running the DP.
///
/// On success the result describes a single match run; on failure the result
/// is cleared but the cached boundaries in `state` are still updated, exactly
/// like the reference (it stores them before comparing the sequences).
#[allow(clippy::too_many_arguments)]
pub fn cal_exact_exz(
    region: &OverlapRegion,
    targets: &TargetSource,
    extz: &mut BitExtz,
    query: &[u8],
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    state: &mut ExzAdvState,
    mode: i64,
) -> bool {
    extz.clear_align();
    extz.result.threshold = 0;
    extz.cigar.packed.clear();

    let query_length = query_end - query_start;
    let target_total_length = targets.length(region.y_identifier as usize);
    let mut target_start = target_start;
    let mut target_end = target_end;

    if mode == MODE_SEMI_GLOBAL {
        target_start = (query_start - i64::from(region.x_position_start))
            + i64::from(region.y_position_start);
        target_start += region
            .fake_cigar
            .y_start_offset(query_start)
            .expect("the fake cigar must cover the query start");
        target_end = target_start + query_length;
    } else if mode == MODE_FORWARD_EXTENSION {
        target_end = target_start + query_length;
    } else if mode == MODE_BACKWARD_EXTENSION {
        target_start = target_end - query_length;
    }

    if target_start < 0 {
        target_start = 0;
    }
    if target_start > target_total_length {
        target_start = target_total_length;
    }
    if target_end > target_total_length {
        target_end = target_total_length;
    }
    let target_length = target_end - target_start;
    if query_length != target_length {
        return false;
    }

    let query_sequence = &query[query_start as usize..query_end as usize];
    let target_sequence = targets.sub_region(
        region.y_identifier as usize,
        target_start,
        target_length,
        region.y_strand as u8,
    );
    // `retrieve_str_seq_exz` stores the boundaries before the comparison, so
    // they are updated even when the sequences differ.
    state.pattern_start = target_start;
    state.pattern_end = target_end;

    if query_sequence != target_sequence.as_slice() {
        return false;
    }

    extz.result.error = 0;
    extz.cigar.push(0, query_length as u32);
    extz.result.pattern_length = target_length as i32;
    extz.result.pattern_start = 0;
    extz.result.pattern_end = target_length as i32 - 1;
    extz.result.text_length = query_length as i32;
    extz.result.text_start = 0;
    extz.result.text_end = query_length as i32 - 1;
    extz.result.pattern_start += target_start as i32;
    extz.result.pattern_end += target_start as i32;
    extz.result.text_start += query_start as i32;
    extz.result.text_end += query_start as i32;
    true
}

/// `cal_exz_infi_adv`: aligns `[query_start, query_end)` against the target
/// interval, mapping it first when `mode` is 3, and returns whether the
/// alignment stayed inside `threshold`.
///
/// `state.threshold` is the `pthre` of the reference: once a threshold has
/// been tried, any later call with a threshold no larger is skipped, because
/// a band that wide already failed.
#[allow(clippy::too_many_arguments)]
pub fn cal_exz_infi_adv(
    region: &OverlapRegion,
    targets: &TargetSource,
    extz: &mut BitExtz,
    query: &[u8],
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    state: &mut ExzAdvState,
    threshold: i64,
    query_total_length: i64,
    mode: i64,
) -> bool {
    extz.clear_align();

    let mut query_start = query_start;
    let mut query_end = query_end;
    let mut target_start = target_start;
    let mut target_end = target_end;
    let mut threshold = threshold;
    let mut auxiliary_begin = 0i64;

    let target_total_length = targets.length(region.y_identifier as usize);
    let mut query_length = query_end - query_start;
    let mut target_length = target_end - target_start;
    let mut diagonal_distance = query_length.max(target_length);

    // The clamp for the coordinate mapping uses the distance before the
    // mapping adjusts the interval.
    let mapping_threshold = threshold.min(diagonal_distance);
    if mode == MODE_SEMI_GLOBAL {
        let (start, end, auxiliary) =
            update_semi_coord(region, targets, query_start, query_end, mapping_threshold);
        target_start = start;
        target_end = end;
        auxiliary_begin = auxiliary;
    } else if mode == MODE_FORWARD_EXTENSION || mode == MODE_BACKWARD_EXTENSION {
        adjust_ext_offset(
            &mut query_start,
            &mut query_end,
            &mut target_start,
            &mut target_end,
            query_total_length,
            target_total_length,
            mapping_threshold,
            mode,
        );
    }

    if !(query_end > query_start
        && target_end > target_start
        && target_start != -1
        && target_end != -1)
    {
        return false;
    }

    query_length = query_end - query_start;
    target_length = target_end - target_start;
    diagonal_distance = query_length.max(target_length);
    if threshold > diagonal_distance {
        threshold = diagonal_distance;
    }
    if threshold <= state.threshold {
        return false;
    }
    state.threshold = threshold;

    let query_sequence = &query[query_start as usize..query_end as usize];
    let target_sequence = targets.sub_region(
        region.y_identifier as usize,
        target_start,
        target_length,
        region.y_strand as u8,
    );
    state.pattern_start = target_start;
    state.pattern_end = target_end;

    let threshold = threshold as i32;
    match mode {
        MODE_GLOBAL => cal_exz_global(extz, &target_sequence, query_sequence, threshold),
        MODE_FORWARD_EXTENSION => {
            cal_exz_extension_0(extz, &target_sequence, query_sequence, threshold)
        }
        MODE_BACKWARD_EXTENSION => {
            cal_exz_extension_1(extz, &target_sequence, query_sequence, threshold)
        }
        _ => cal_exz_semi(
            extz,
            &target_sequence,
            query_sequence,
            threshold,
            auxiliary_begin as i32,
        ),
    }

    if extz.is_aligned() {
        extz.result.pattern_start += target_start as i32;
        extz.result.pattern_end += target_start as i32;
        extz.result.text_start += query_start as i32;
        extz.result.text_end += query_start as i32;
        return true;
    }
    false
}

// ---------------------------------------------------------------------------
// The `_simi` variants. The UL second-pass alignment (`hc_aln_exz_simi_adv`)
// calls the no-traceback kernels; the dispatch below mirrors
// `cal_exz_global_simi` and friends (`Correct.cpp` 15496-15607). The 128-bit
// and multi-word `infi` instances are word-level equivalent, so one
// width-generic kernel covers both, exactly like the `_trace` layer.
// ---------------------------------------------------------------------------

/// `cal_exz_global_simi`: global alignment without traceback.
pub fn cal_exz_global_simi(extz: &mut BitExtz, pattern: &[u8], text: &[u8], threshold: i32) {
    if band_word_count(threshold) <= 1 {
        global_64(
            pattern,
            pattern.len() as i32,
            text,
            text.len() as i32,
            threshold,
            &mut extz.result,
        );
    } else {
        global_no_trace(
            pattern,
            pattern.len() as i32,
            text,
            text.len() as i32,
            threshold,
            &mut extz.result,
        );
    }
}

/// `cal_exz_extension_0_simi`: forward extension without traceback.
pub fn cal_exz_extension_0_simi(extz: &mut BitExtz, pattern: &[u8], text: &[u8], threshold: i32) {
    if band_word_count(threshold) <= 1 {
        extension_forward_64(
            pattern,
            pattern.len() as i32,
            text,
            text.len() as i32,
            threshold,
            &mut extz.result,
        );
    } else {
        extension_forward_no_trace(
            pattern,
            pattern.len() as i32,
            text,
            text.len() as i32,
            threshold,
            &mut extz.result,
        );
    }
}

/// `cal_exz_extension_1_simi`: backward extension without traceback.
pub fn cal_exz_extension_1_simi(extz: &mut BitExtz, pattern: &[u8], text: &[u8], threshold: i32) {
    if band_word_count(threshold) <= 1 {
        extension_backward_64(
            pattern,
            pattern.len() as i32,
            text,
            text.len() as i32,
            threshold,
            &mut extz.result,
        );
    } else {
        extension_backward_no_trace(
            pattern,
            pattern.len() as i32,
            text,
            text.len() as i32,
            threshold,
            &mut extz.result,
        );
    }
}

/// `cal_exz_semi_simi`: semi-global alignment on an absent anti-diagonal
/// without traceback.
pub fn cal_exz_semi_simi(
    extz: &mut BitExtz,
    pattern: &[u8],
    text: &[u8],
    threshold: i32,
    auxiliary_begin: i32,
) {
    if band_word_count(threshold) <= 1 {
        semi_global_64_absent_diagonal(
            pattern,
            pattern.len() as i32,
            text,
            text.len() as i32,
            threshold,
            auxiliary_begin,
            &mut extz.result,
        );
    } else {
        semi_global_absent_diagonal_no_trace(
            pattern,
            pattern.len() as i32,
            text,
            text.len() as i32,
            threshold,
            auxiliary_begin,
            &mut extz.result,
        );
    }
}

/// `adjust_ext_offset_fixed_t` (`Correct.cpp` 15811): the extension shrink of
/// [`adjust_ext_offset`] with caller-supplied clamping bounds instead of the
/// sequence lengths.
#[allow(clippy::too_many_arguments)]
pub fn adjust_ext_offset_fixed_t(
    query_start: &mut i64,
    query_end: &mut i64,
    target_start: &mut i64,
    target_end: &mut i64,
    query_minimum: i64,
    query_maximum: i64,
    target_minimum: i64,
    target_maximum: i64,
    threshold: i64,
    mode: i64,
) {
    if mode == MODE_FORWARD_EXTENSION {
        let query_offset = query_maximum - *query_start;
        let target_offset = target_maximum - *target_start;
        if query_offset <= target_offset {
            *query_end = query_maximum;
            *target_end = *target_start + query_offset + threshold;
        } else {
            *target_end = target_maximum;
            *query_end = *query_start + target_offset + threshold;
        }
    } else if mode == MODE_BACKWARD_EXTENSION {
        let query_offset = *query_end - query_minimum;
        let target_offset = *target_end - target_minimum;
        if query_offset <= target_offset {
            *query_start = query_minimum;
            *target_start = *target_end - query_offset - threshold;
        } else {
            *target_start = target_minimum;
            *query_start = *query_end - target_offset - threshold;
        }
    }
    if *query_start < query_minimum {
        *query_start = query_minimum;
    }
    if *target_start < target_minimum {
        *target_start = target_minimum;
    }
    if *query_end > query_maximum {
        *query_end = query_maximum;
    }
    if *target_end > target_maximum {
        *target_end = target_maximum;
    }
}

/// `cal_exact_simi_exz` (`Correct.cpp` 15918): the exact-match shortcut of
/// [`cal_exact_exz`] with the target identity and strand taken from the
/// caller (the UL second pass may call without an overlap region).
///
/// `state` carries only `pattern_start`/`pattern_end` (`pts`/`pte`); the
/// `pthre` half of the state is not touched.
#[allow(clippy::too_many_arguments)]
pub fn cal_exact_simi_exz(
    targets: &TargetSource,
    extz: &mut BitExtz,
    query: &[u8],
    identifier: usize,
    strand: u8,
    query_start: i64,
    query_end: i64,
    mut target_start: i64,
    mut target_end: i64,
    state: &mut ExzAdvState,
    mode: i64,
    region: Option<&OverlapRegion>,
) -> bool {
    extz.clear_align();
    extz.result.threshold = 0;
    extz.cigar.packed.clear();

    let target_total_length = targets.length(identifier);
    let mut query_length = query_end - query_start;

    if mode == MODE_SEMI_GLOBAL {
        if let Some(region) = region {
            target_start = (query_start - i64::from(region.x_position_start))
                + i64::from(region.y_position_start);
            target_start += region
                .fake_cigar
                .y_start_offset(query_start)
                .expect("the fake cigar must cover the query start");
            target_end = target_start + query_length;
        }
    } else if mode == MODE_FORWARD_EXTENSION {
        target_end = target_start + query_length;
    } else if mode == MODE_BACKWARD_EXTENSION {
        target_start = target_end - query_length;
    }
    if target_start < 0 {
        target_start = 0;
    }
    if target_start > target_total_length {
        target_start = target_total_length;
    }
    if target_end > target_total_length {
        target_end = target_total_length;
    }
    query_length = query_end - query_start;
    let target_length = target_end - target_start;
    if query_length != target_length {
        return false;
    }

    let query_sequence = &query[query_start as usize..query_end as usize];
    let target_sequence = targets.sub_region(identifier, target_start, target_length, strand);
    state.pattern_start = target_start;
    state.pattern_end = target_end;

    if query_sequence != target_sequence.as_slice() {
        return false;
    }

    extz.result.error = 0;
    extz.cigar.push(0, query_length as u32);
    extz.result.pattern_length = target_length as i32;
    extz.result.pattern_start = 0;
    extz.result.pattern_end = target_length as i32 - 1;
    extz.result.text_length = query_length as i32;
    extz.result.text_start = 0;
    extz.result.text_end = query_length as i32 - 1;

    extz.result.pattern_start += target_start as i32;
    extz.result.pattern_end += target_start as i32;
    extz.result.text_start += query_start as i32;
    extz.result.text_end += query_start as i32;
    true
}

/// `cal_exz_infi_simi_adv` (`Correct.cpp` 15838): the no-traceback sibling of
/// [`cal_exz_infi_adv`]. The target identity and strand come from the caller;
/// `region` is only consulted for the semi-global mapping and is `None` on
/// every UL second-pass call site.
#[allow(clippy::too_many_arguments)]
pub fn cal_exz_infi_simi_adv(
    targets: &TargetSource,
    extz: &mut BitExtz,
    query: &[u8],
    identifier: usize,
    strand: u8,
    mut query_start: i64,
    mut query_end: i64,
    mut target_start: i64,
    mut target_end: i64,
    state: &mut ExzAdvState,
    mut threshold: i64,
    query_minimum: i64,
    query_maximum: i64,
    target_minimum: i64,
    target_maximum: i64,
    mode: i64,
    region: Option<&OverlapRegion>,
) -> bool {
    extz.clear_align();
    let mut auxiliary_begin = 0i64;

    let mut query_length = query_end - query_start;
    let mut target_length = target_end - target_start;
    let mut diagonal_distance = query_length.max(target_length);
    if mode == MODE_SEMI_GLOBAL {
        if let Some(region) = region {
            let (start, end, auxiliary) = update_semi_coord(
                region,
                targets,
                query_start,
                query_end,
                threshold.min(diagonal_distance),
            );
            target_start = start;
            target_end = end;
            auxiliary_begin = auxiliary;
        }
    } else if mode == MODE_FORWARD_EXTENSION || mode == MODE_BACKWARD_EXTENSION {
        adjust_ext_offset_fixed_t(
            &mut query_start,
            &mut query_end,
            &mut target_start,
            &mut target_end,
            query_minimum,
            query_maximum,
            target_minimum,
            target_maximum,
            threshold.min(diagonal_distance),
            mode,
        );
    }

    if !(query_end > query_start
        && target_end > target_start
        && target_start != -1
        && target_end != -1)
    {
        return false;
    }

    query_length = query_end - query_start;
    target_length = target_end - target_start;
    diagonal_distance = query_length.max(target_length);
    if threshold > diagonal_distance {
        threshold = diagonal_distance;
    }
    if threshold <= state.threshold {
        return false;
    }
    state.threshold = threshold;

    let query_sequence = &query[query_start as usize..query_end as usize];
    let target_sequence = targets.sub_region(identifier, target_start, target_length, strand);
    state.pattern_start = target_start;
    state.pattern_end = target_end;

    let threshold = threshold as i32;
    match mode {
        MODE_GLOBAL => cal_exz_global_simi(extz, &target_sequence, query_sequence, threshold),
        MODE_FORWARD_EXTENSION => {
            cal_exz_extension_0_simi(extz, &target_sequence, query_sequence, threshold)
        }
        MODE_BACKWARD_EXTENSION => {
            cal_exz_extension_1_simi(extz, &target_sequence, query_sequence, threshold)
        }
        _ => cal_exz_semi_simi(
            extz,
            &target_sequence,
            query_sequence,
            threshold,
            auxiliary_begin as i32,
        ),
    }

    if extz.is_aligned() {
        extz.result.pattern_start += target_start as i32;
        extz.result.pattern_end += target_start as i32;
        extz.result.text_start += query_start as i32;
        extz.result.text_end += query_start as i32;
        return true;
    }
    false
}

/// `cal_estimate_err`: estimates the number of errors between `[query_start,
/// query_end)` and its target by summing the stored window errors, weighting
/// the partially covered windows by the fraction of the window covered, and
/// charging `error_rate` for the query bases outside every window.
pub fn cal_estimate_err(
    region: &OverlapRegion,
    window_length: i64,
    query_start: i64,
    query_end: i64,
    error_rate: f64,
) -> i64 {
    // The estimate uses the requested interval, before clipping.
    let estimated = ((query_end - query_start) as f64 * error_rate) as i64;
    let windows = &region.window_cigar.windows;
    let window_count = windows.len() as i64;
    if window_count == 0 {
        return estimated;
    }

    let mut query_start = query_start;
    let mut query_end = query_end;
    if query_start < i64::from(region.x_position_start) {
        query_start = i64::from(region.x_position_start);
    }
    if query_end > i64::from(region.x_position_end) + 1 {
        query_end = i64::from(region.x_position_end) + 1;
    }

    let window_start = query_start / window_length * window_length;
    let (mut identifier, _) = window_id_by_start(
        region.x_position_start as i32,
        region.x_position_end as i32,
        window_start,
        window_length,
    );
    if identifier >= window_count {
        identifier = window_count - 1;
    }

    let mut k = identifier;
    while k < window_count && query_start > i64::from(windows[k as usize].x_end) {
        k += 1;
    }
    if k == window_count {
        return estimated;
    }
    while k >= 0 && query_start < i64::from(windows[k as usize].x_start) {
        k -= 1;
    }
    if k < 0 {
        k = 0;
    }

    let mut total: i64 = 0;
    let mut covered: i64 = 0;
    while k < window_count && i64::from(windows[k as usize].x_start) < query_end {
        let window = &windows[k as usize];
        if window.y_end == -1 {
            k += 1;
            continue;
        }
        let window_start = i64::from(window.x_start);
        let window_end = i64::from(window.x_end) + 1;
        let overlap_start = query_start.max(window_start);
        let overlap_end = query_end.min(window_end);
        let overlap = if overlap_end > overlap_start {
            overlap_end - overlap_start
        } else {
            0
        };
        if overlap == 0 {
            k += 1;
            continue;
        }
        covered += overlap;
        if overlap == window_end - window_start {
            total += i64::from(window.error);
        } else {
            // The reference performs this addition in `double` and truncates.
            total = (total as f64
                + (f64::from(window.error)) * (overlap as f64) / ((window_end - window_start) as f64))
                as i64;
        }
        k += 1;
    }

    total = (total as f64 + ((query_end - query_start) - covered) as f64 * error_rate) as i64;
    total
}

/// `cal_estimate_err_hc`: [`cal_estimate_err`] that additionally reports
/// whether every covered window had zero error and the window coordinates
/// chain into exactly `[target_start, target_end)`.
pub fn cal_estimate_err_hc(
    region: &OverlapRegion,
    window_length: i64,
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    error_rate: f64,
) -> (i64, bool) {
    let estimated = ((query_end - query_start) as f64 * error_rate) as i64;
    let windows = &region.window_cigar.windows;
    let window_count = windows.len() as i64;
    if window_count == 0 {
        return (estimated, false);
    }

    let mut query_start = query_start;
    let mut query_end = query_end;
    if query_start < i64::from(region.x_position_start) {
        query_start = i64::from(region.x_position_start);
    }
    if query_end > i64::from(region.x_position_end) + 1 {
        query_end = i64::from(region.x_position_end) + 1;
    }

    let window_start = query_start / window_length * window_length;
    let (mut identifier, _) = window_id_by_start(
        region.x_position_start as i32,
        region.x_position_end as i32,
        window_start,
        window_length,
    );
    if identifier >= window_count {
        identifier = window_count - 1;
    }

    let mut k = identifier;
    while k < window_count && query_start > i64::from(windows[k as usize].x_end) {
        k += 1;
    }
    if k == window_count {
        return (estimated, false);
    }
    while k >= 0 && query_start < i64::from(windows[k as usize].x_start) {
        k -= 1;
    }
    if k < 0 {
        k = 0;
    }

    let mut total: i64 = 0;
    let mut covered: i64 = 0;
    let mut exact = true;
    let mut exact_start: i64 = -1;
    let mut exact_end: i64 = -1;
    while k < window_count && i64::from(windows[k as usize].x_start) < query_end {
        let window = &windows[k as usize];
        if window.y_end == -1 {
            k += 1;
            continue;
        }
        let window_start = i64::from(window.x_start);
        let window_end = i64::from(window.x_end) + 1;
        let overlap_start = query_start.max(window_start);
        let overlap_end = query_end.min(window_end);
        let overlap = if overlap_end > overlap_start {
            overlap_end - overlap_start
        } else {
            0
        };
        if overlap == 0 {
            k += 1;
            continue;
        }
        covered += overlap;
        if overlap == window_end - window_start {
            total += i64::from(window.error);
        } else {
            total = (total as f64
                + (f64::from(window.error)) * (overlap as f64) / ((window_end - window_start) as f64))
                as i64;
        }

        if window.error > 0 {
            exact = false;
        }
        if exact {
            // The query-side margins of the covered part must match the
            // target-side margins. `q[1]`/`t[1]` of the reference are computed
            // and never read, and `t[0] == t[0]` is always true, so only the
            // left margins matter.
            let query_leading = overlap_start - window_start;
            let target_window_start =
                i64::from(window.y_end) + 1 - (i64::from(window.x_end) + 1 - window_start);
            let target_window_end = i64::from(window.y_end) + 1;
            let target_overlap_start = target_start.max(target_window_start);
            let target_overlap_end = target_end.min(target_window_end);
            let target_overlap = if target_overlap_end > target_overlap_start {
                target_overlap_end - target_overlap_start
            } else {
                0
            };
            let target_leading = target_overlap_start - target_window_start;
            if target_overlap != 0 && query_leading == target_leading {
                if exact_end == -1 {
                    exact_start = target_overlap_start;
                    exact_end = target_overlap_end;
                } else if exact_end == target_overlap_start {
                    exact_end = target_overlap_end;
                } else {
                    exact = false;
                }
            } else {
                exact = false;
            }
        }
        k += 1;
    }

    total = (total as f64 + ((query_end - query_start) - covered) as f64 * error_rate) as i64;
    if exact && (query_end - query_start) == covered {
        if (query_end - query_start) == (exact_end - exact_start)
            && exact_start == target_start
            && exact_end == target_end
        {
            return (total, true);
        }
    }
    (total, false)
}

/// `push_wcigar`: attaches the cells of `extz` to the window at `index`.
fn push_window_cigar(aux_o: &mut WindowCigar, index: usize, extz: &BitExtz) {
    let cigar_start = aux_o.cells.len() as u32;
    let cigar_length = extz.cigar.packed.len() as u32;
    aux_o.windows[index].cigar_start = cigar_start;
    aux_o.windows[index].cigar_length = cigar_length;
    aux_o.cells.extend_from_slice(&extz.cigar.packed);
}

/// `pop_trace` on a standalone packed cigar: the first operation with its
/// consecutive equal-code cells merged, plus the index just past them.
fn pop_trace_first(cells: &[u16], mut index: usize) -> (u16, u32, usize) {
    let code = cells[index] >> 14;
    let mut length = u32::from(cells[index] & 0x3fff);
    index += 1;
    while index < cells.len() && cells[index] >> 14 == code {
        length += u32::from(cells[index] & 0x3fff);
        index += 1;
    }
    (code, length, index)
}

/// `push_trace` on the shared cell buffer: appends a run, splitting it when it
/// does not fit in the fourteen bits available for the length.
fn push_packed_run(cells: &mut Vec<u16>, code: u16, mut length: u32) {
    const MAXIMUM_RUN: u32 = 0x3fff;
    let shifted = code << 14;
    while length >= MAXIMUM_RUN {
        cells.push(shifted | MAXIMUM_RUN as u16);
        length -= MAXIMUM_RUN;
    }
    if length > 0 {
        cells.push(shifted | length as u16);
    }
}

/// `append_wcigar`: merges the cigar of `extz` into the window at `index`,
/// joining the last run of the stored cigar with the first run of the new one
/// when they share an operation.
fn append_window_cigar(aux_o: &mut WindowCigar, index: usize, extz: &BitExtz) {
    if aux_o.windows[index].cigar_length > 0 {
        assert_eq!(
            aux_o.windows[index].cigar_start as usize + aux_o.windows[index].cigar_length as usize,
            aux_o.cells.len()
        );
        if !extz.cigar.packed.is_empty() {
            let last = *aux_o
                .cells
                .last()
                .expect("a window with a cigar has at least one cell");
            let last_code = last >> 14;
            let last_length = u32::from(last & 0x3fff);
            let (code, mut length, next) = pop_trace_first(&extz.cigar.packed, 0);
            if last_code == code {
                length += last_length;
                aux_o.cells.pop();
                aux_o.windows[index].cigar_length -= 1;
            }
            push_packed_run(&mut aux_o.cells, code, length);
            aux_o.windows[index].cigar_length =
                aux_o.cells.len() as u32 - aux_o.windows[index].cigar_start;

            let remaining = extz.cigar.packed.len() - next;
            if remaining > 0 {
                aux_o.cells.extend_from_slice(&extz.cigar.packed[next..]);
                aux_o.windows[index].cigar_length += remaining as u32;
            }
        }
    } else {
        // `exz` is the first cigar of the window.
        push_window_cigar(aux_o, index, extz);
    }
}

/// `push_alnw`: records the alignment in `extz` as a window of `aux_o`,
/// merging it with the previous window when the two are adjacent in both
/// reads and their combined error still fits in an `i16`.
pub fn push_alignment_window(aux_o: &mut WindowCigar, extz: &BitExtz) {
    let window_count = aux_o.windows.len();
    if window_count > 0 {
        let last = &aux_o.windows[window_count - 1];
        let merged_error = i64::from(last.error) + i64::from(extz.result.error);
        if last.cigar_length > 0
            && last.x_end + 1 == extz.result.text_start
            && last.y_end + 1 == extz.result.pattern_start
            && merged_error < i64::from(i16::MAX)
        {
            let last = &mut aux_o.windows[window_count - 1];
            last.x_end = extz.result.text_end;
            last.y_end = extz.result.pattern_end;
            last.error += extz.result.error as i16;
            append_window_cigar(aux_o, window_count - 1, extz);
            return;
        }
    }

    // The reference leaves `extra_begin`/`extra_end` unwritten here; no
    // consumer reads them for an aligned window because `is_ualn_win` already
    // fails on `error != INT16_MAX`. They are recorded as zero.
    let index = aux_o.append_window(
        extz.result.text_start,
        extz.result.text_end,
        extz.result.pattern_start,
        extz.result.pattern_end,
        extz.result.error as i16,
        0,
        0,
        0,
    );
    push_window_cigar(aux_o, index, extz);
}

/// `push_unmap_alnw`: records a region no alignment covers, `[qs, qe]` on the
/// query against `[ts, te]` on the target.
pub fn push_unmapped_window(
    aux_o: &mut WindowCigar,
    query_start: i32,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    mode: i64,
) {
    aux_o.append_window(
        query_start,
        query_end as i32,
        target_start as i32,
        target_end as i32,
        i16::MAX,
        -1,
        -1,
        mode as i16,
    );
}

/// `push_replace_alnw`: records a region replaced by a fake cigar without an
/// alignment, `[qs, qe]` on the query against `[ts, te]` on the target.
pub fn push_replace_window(
    aux_o: &mut WindowCigar,
    query_start: i32,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    mode: i64,
) {
    aux_o.append_window(
        query_start,
        query_end as i32,
        target_start as i32,
        target_end as i32,
        0,
        0,
        0,
        mode as i16,
    );
}

/// `push_replace_alnw_adv`: like [`push_replace_window`], but the target
/// coordinates are taken from the existing windows of `region` between the
/// two query boundaries, splitting the replacement when the accumulated error
/// would reach `INT16_MAX`.
///
/// The reference asserts that both boundaries sit on window boundaries; the
/// asserts are kept because the reference binary is built with them enabled.
pub fn push_replace_window_adv(
    region: &OverlapRegion,
    window_length: i64,
    aux_o: &mut WindowCigar,
    query_start: i32,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    mode: i64,
) {
    let (start_id, _) = window_id_by_start(
        region.x_position_start as i32,
        region.x_position_end as i32,
        i64::from(query_start),
        window_length,
    );
    let (end_id, _) = window_id_by_start(
        region.x_position_start as i32,
        region.x_position_end as i32,
        query_end + 1,
        window_length,
    );
    assert_eq!(
        region.window_cigar.windows[start_id as usize].x_start,
        query_start
    );
    assert_eq!(
        i64::from(region.window_cigar.windows[end_id as usize].x_end),
        query_end + 1
    );

    let mut index = aux_o.append_window(
        query_start,
        query_start,
        target_start as i32,
        target_start as i32,
        0,
        0,
        0,
        mode as i16,
    );
    let mut error: i64 = 0;
    for k in start_id..=end_id {
        let source = &region.window_cigar.windows[k as usize];
        assert_ne!(source.y_end, -1);
        let summed = error + i64::from(source.error);
        if summed < i64::from(i16::MAX) {
            error += i64::from(source.error);
            aux_o.windows[index].x_end = source.x_end;
            aux_o.windows[index].y_end = source.y_end;
        } else {
            error = i64::from(source.error);
            index = aux_o.append_window(
                source.x_start,
                source.x_end,
                source.y_start,
                source.y_end,
                0,
                0,
                0,
                mode as i16,
            );
        }
        aux_o.windows[index].error = error as i16;
    }
    aux_o.windows[index].x_end = query_end as i32;
    aux_o.windows[index].y_end = target_end as i32;
}

/// `set_exact_exz`: fills `extz` with a single match run covering
/// `[query_start, query_end)`, used when the error estimate already says the
/// region is identical.
pub fn set_exact_exz(
    extz: &mut BitExtz,
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
) {
    extz.clear_align();
    extz.result.threshold = 0;
    extz.cigar.packed.clear();
    extz.result.error = 0;
    extz.cigar.push(0, (query_end - query_start) as u32);
    extz.result.pattern_length = (target_end - target_start) as i32;
    extz.result.pattern_start = 0;
    extz.result.pattern_end = extz.result.pattern_length - 1;
    extz.result.text_length = (query_end - query_start) as i32;
    extz.result.text_start = 0;
    extz.result.text_end = extz.result.text_length - 1;
    extz.result.pattern_start += target_start as i32;
    extz.result.pattern_end += target_start as i32;
    extz.result.text_start += query_start as i32;
    extz.result.text_end += query_start as i32;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::correction::window::WindowCigar;
    use crate::overlap::region::{FakeCigar, OverlapRegion};
    use crate::reads::read_store::ReadStore;

    fn region_with_windows(
        x_start: u32,
        x_end: u32,
        windows: &[(i32, i32, i32, i32, i16)],
    ) -> OverlapRegion {
        let mut window_cigar = WindowCigar::default();
        for &(x_start, x_end, y_start, y_end, error) in windows {
            window_cigar.append_window(x_start, x_end, y_start, y_end, error, 0, 0, 0);
        }
        OverlapRegion {
            x_position_start: x_start,
            x_position_end: x_end,
            y_identifier: 0,
            window_cigar,
            ..OverlapRegion::default()
        }
    }

    #[test]
    fn scale_ed_thre_rounds_to_whole_words() {
        // 2*err+1 fits in one word up to err == 31.
        assert_eq!(scale_ed_thre(0, 2047), 31);
        assert_eq!(scale_ed_thre(15, 2047), 31);
        assert_eq!(scale_ed_thre(31, 2047), 31);
        // err == 32 needs a band 65 bits wide, which rounds to two words.
        assert_eq!(scale_ed_thre(32, 2047), 63);
        assert_eq!(scale_ed_thre(63, 2047), 63);
        assert_eq!(scale_ed_thre(64, 2047), 95);
        // The result is capped at the maximum error.
        assert_eq!(scale_ed_thre(64, 50), 50);
    }

    #[test]
    fn initialise_window_alignment_pads_and_clips() {
        // The mapped start is before the target: the extra bases are reported
        // as `auxiliary_begin` and the interval starts at zero.
        let window = initialise_window_alignment(15, 5, 1000, 130).unwrap();
        assert_eq!(window.start, 0);
        assert_eq!(window.length, 120);
        assert_eq!(window.auxiliary_begin, 10);

        // The alignment length is clipped to what the target can supply.
        let window = initialise_window_alignment(15, 950, 1000, 70).unwrap();
        assert_eq!(window.start, 935);
        assert_eq!(window.length, 65);
        assert_eq!(window.auxiliary_begin, 0);

        // A start outside the target fails.
        assert!(initialise_window_alignment(15, -1, 1000, 130).is_none());
        assert!(initialise_window_alignment(15, 1000, 1000, 130).is_none());
        // So does a target too short to hold the alignment plus the margin.
        assert!(initialise_window_alignment(15, 900, 1000, 500).is_none());
    }

    #[test]
    fn adjust_ext_offset_clips_to_the_available_bases() {
        let (mut query_start, mut query_end, mut target_start, mut target_end) = (20, 80, 500, 600);
        adjust_ext_offset(
            &mut query_start,
            &mut query_end,
            &mut target_start,
            &mut target_end,
            100,
            1000,
            10,
            MODE_FORWARD_EXTENSION,
        );
        assert_eq!(
            (query_start, query_end, target_start, target_end),
            (20, 100, 500, 590)
        );

        let (mut query_start, mut query_end, mut target_start, mut target_end) = (20, 80, 500, 600);
        adjust_ext_offset(
            &mut query_start,
            &mut query_end,
            &mut target_start,
            &mut target_end,
            100,
            1000,
            10,
            MODE_BACKWARD_EXTENSION,
        );
        assert_eq!(
            (query_start, query_end, target_start, target_end),
            (0, 80, 510, 600)
        );
    }

    #[test]
    fn update_semi_coord_applies_the_fake_cigar_shift() {
        let mut reads = ReadStore::default();
        reads.push_read("target", &vec![b'A'; 5000]);
        let mut region = OverlapRegion {
            x_position_start: 100,
            x_position_end: 1999,
            y_position_start: 3000,
            y_identifier: 0,
            ..OverlapRegion::default()
        };
        region.fake_cigar = FakeCigar::default();
        region.fake_cigar.push(100, 0);
        region.fake_cigar.push(1999, 0);

        // qs == x_pos_s maps to y_pos_s.
        let (start, end, auxiliary) = update_semi_coord(&region, &TargetSource::Reads(&reads), 100, 900, 15);
        assert_eq!((start, end, auxiliary), (2985, 3815, 0));

        // A query position shifted by the fake cigar moves the target too.
        region.fake_cigar.clear();
        region.fake_cigar.push(100, 0);
        region.fake_cigar.push(1100, 40);
        region.fake_cigar.push(1999, 40);
        let (start, _, _) = update_semi_coord(&region, &TargetSource::Reads(&reads), 1000, 1500, 15);
        assert_eq!(start, 3885);
        let (start, _, _) = update_semi_coord(&region, &TargetSource::Reads(&reads), 1200, 1500, 15);
        assert_eq!(start, 4140 - 15);
    }

    #[test]
    fn cal_estimate_err_weights_partial_windows() {
        let region = region_with_windows(
            0,
            99,
            &[(0, 49, 100, 149, 10), (50, 99, 150, 199, 20)],
        );
        // The whole overlap, both windows fully covered.
        assert_eq!(cal_estimate_err(&region, 50, 0, 100, 0.1), 30);
        // Half of the first window: ten errors over fifty bases becomes five.
        assert_eq!(cal_estimate_err(&region, 50, 0, 25, 0.1), 5);
        // A query beyond every window falls back to the raw estimate.
        assert_eq!(cal_estimate_err(&region, 50, 100, 120, 0.1), 2);
    }

    #[test]
    fn cal_estimate_err_hc_detects_exact_tiling() {
        let region = region_with_windows(
            0,
            99,
            &[(0, 49, 100, 149, 0), (50, 99, 150, 199, 0)],
        );
        assert_eq!(
            cal_estimate_err_hc(&region, 50, 0, 100, 100, 200, 0.1),
            (0, true)
        );
        // One error unmasks the exactness.
        let region = region_with_windows(
            0,
            99,
            &[(0, 49, 100, 149, 0), (50, 99, 150, 199, 3)],
        );
        assert_eq!(
            cal_estimate_err_hc(&region, 50, 0, 100, 100, 200, 0.1),
            (3, false)
        );
        // Window coordinates that do not chain into the target fail too.
        let region = region_with_windows(
            0,
            99,
            &[(0, 49, 100, 149, 0), (50, 99, 150, 199, 0)],
        );
        assert_eq!(
            cal_estimate_err_hc(&region, 50, 0, 100, 100, 201, 0.1),
            (0, false)
        );
    }

    fn alignment(
        text_start: i32,
        text_end: i32,
        pattern_start: i32,
        pattern_end: i32,
        error: i32,
        runs: &[(u16, u32)],
    ) -> BitExtz {
        let mut extz = BitExtz::new();
        extz.result.threshold = error;
        extz.result.error = error;
        extz.result.text_start = text_start;
        extz.result.text_end = text_end;
        extz.result.text_length = text_end + 1 - text_start;
        extz.result.pattern_start = pattern_start;
        extz.result.pattern_end = pattern_end;
        extz.result.pattern_length = pattern_end + 1 - pattern_start;
        for &(code, length) in runs {
            extz.cigar.push(code, length);
        }
        extz
    }

    #[test]
    fn push_alignment_window_merges_adjacent_windows_and_cigars() {
        let mut windows = WindowCigar::default();
        let first = alignment(0, 99, 100, 199, 5, &[(0, 50), (1, 3), (0, 47)]);
        push_alignment_window(&mut windows, &first);
        assert_eq!(windows.windows.len(), 1);
        assert_eq!(windows.cells.len(), 3);

        // Adjacent in both reads: the windows merge and the shared operation
        // at the junction is folded into one run.
        let second = alignment(100, 199, 200, 299, 7, &[(1, 2), (0, 98)]);
        push_alignment_window(&mut windows, &second);
        assert_eq!(windows.windows.len(), 1);
        let window = &windows.windows[0];
        assert_eq!(
            (window.x_start, window.x_end, window.y_start, window.y_end, window.error),
            (0, 199, 100, 299, 12)
        );
        assert_eq!(
            windows.cells,
            vec![(0 << 14) | 50, (1 << 14) | 3, (0 << 14) | 47, (1 << 14) | 2, (0 << 14) | 98]
        );

        // A third window continues the merge; now the first new operation is
        // a match and is added to the trailing match run.
        let third = alignment(200, 299, 300, 399, 3, &[(0, 5), (0, 10)]);
        push_alignment_window(&mut windows, &third);
        assert_eq!(windows.windows.len(), 1);
        assert_eq!(windows.windows[0].error, 15);
        assert_eq!(
            windows.cells,
            vec![(0 << 14) | 50, (1 << 14) | 3, (0 << 14) | 47, (1 << 14) | 2, (0 << 14) | 113]
        );
    }

    #[test]
    fn push_alignment_window_starts_a_new_window_on_a_gap() {
        let mut windows = WindowCigar::default();
        push_alignment_window(&mut windows, &alignment(0, 99, 100, 199, 5, &[(0, 100)]));
        // One query base skipped: the windows stay separate.
        push_alignment_window(&mut windows, &alignment(101, 199, 200, 298, 3, &[(0, 99)]));
        assert_eq!(windows.windows.len(), 2);
        assert_eq!(windows.windows[1].cigar_start, 1);
        assert_eq!(windows.windows[1].cigar_length, 1);
        assert_eq!(windows.cells, vec![(0 << 14) | 100, (0 << 14) | 99]);
    }

    #[test]
    fn push_unmapped_and_replace_windows_record_the_mode() {
        let mut windows = WindowCigar::default();
        push_unmapped_window(&mut windows, 10, 19, 100, 119, 3);
        push_replace_window(&mut windows, 20, 29, 120, 139, 1);
        let unmapped = &windows.windows[0];
        assert_eq!(unmapped.error, i16::MAX);
        assert_eq!((unmapped.extra_begin, unmapped.extra_end), (-1, -1));
        assert_eq!(unmapped.error_threshold, 3);
        let replaced = &windows.windows[1];
        assert_eq!(replaced.error, 0);
        assert_eq!((replaced.extra_begin, replaced.extra_end), (0, 0));
        assert_eq!(replaced.error_threshold, 1);
        assert!(windows.cells.is_empty());
    }

    #[test]
    fn push_replace_window_adv_splits_on_error_overflow() {
        let region = region_with_windows(
            0,
            149,
            &[(0, 49, 100, 149, 10), (50, 99, 150, 199, 20), (100, 149, 200, 249, 30)],
        );
        let mut windows = WindowCigar::default();
        // From the first window to the end of the second: the interval must
        // start on a window start and end just before a window end.
        push_replace_window_adv(&region, 50, &mut windows, 0, 98, 1000, 1198, 0);
        assert_eq!(windows.windows.len(), 1);
        let window = &windows.windows[0];
        assert_eq!(
            (window.x_start, window.x_end, window.y_start, window.y_end, window.error),
            (0, 98, 1000, 1198, 30)
        );
    }

    #[test]
    fn set_exact_exz_describes_one_match_run() {
        let mut extz = BitExtz::new();
        set_exact_exz(&mut extz, 500, 600, 1000, 1100);
        assert!(extz.is_aligned());
        assert_eq!(extz.result.error, 0);
        assert_eq!(extz.cigar.packed, vec![(0 << 14) | 100]);
        assert_eq!(
            (
                extz.result.pattern_start,
                extz.result.pattern_end,
                extz.result.text_start,
                extz.result.text_end
            ),
            (1000, 1099, 500, 599)
        );
    }
}

/// Replay of the reference golden probe.
///
/// `/tmp/ec_calexz_probe.cpp` links against the unmodified hifiasm objects and
/// walks `init_waln`, `adjust_ext_offset`, `update_semi_coord`, both error
/// estimators, `cal_exact_exz` and `cal_exz_infi_adv` over exact, substituted,
/// inserted and reverse-strand targets. The output is stored under
/// `reference/`; the test rebuilds the same reads, windows and call sequence
/// and compares every printed line, including the stale fields the reference
/// leaves behind after a failed call.
#[cfg(test)]
mod probe_goldens {
    use super::*;
    use crate::overlap::region::{FakeCigar, OverlapRegion};
    use crate::reads::read_store::ReadStore;

    struct WindowSpec {
        x_start: i32,
        x_end: i32,
        y_start: i32,
        y_end: i32,
        error: i16,
    }

    struct SemiSpec {
        query_start: i64,
        query_end: i64,
        threshold: i64,
    }

    struct EstimateSpec {
        query_start: i64,
        query_end: i64,
    }

    struct HcEstimateSpec {
        query_start: i64,
        query_end: i64,
        target_start: i64,
        target_end: i64,
    }

    struct ExactSpec {
        query_start: i64,
        query_end: i64,
        target_start: i64,
        target_end: i64,
        mode: i64,
    }

    struct AdvSpec {
        query_start: i64,
        query_end: i64,
        target_start: i64,
        target_end: i64,
        threshold: i64,
        mode: i64,
        reuse_state: bool,
    }

    struct Scenario {
        name: &'static str,
        reads: Vec<Vec<u8>>,
        x_start: u32,
        x_end: u32,
        y_start: u32,
        y_end: u32,
        y_strand: u32,
        gaps: Vec<(u32, i32)>,
        windows: Vec<WindowSpec>,
        semis: Vec<SemiSpec>,
        estimates: Vec<EstimateSpec>,
        hc_estimates: Vec<HcEstimateSpec>,
        exacts: Vec<ExactSpec>,
        advs: Vec<AdvSpec>,
    }

    struct Lcg(u64);

    impl Lcg {
        fn next_u32(&mut self) -> u32 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 33) as u32
        }

        fn sequence(&mut self, length: usize) -> Vec<u8> {
            const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
            (0..length)
                .map(|_| BASES[(self.next_u32() & 3) as usize])
                .collect()
        }
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
                other => *other,
            })
            .collect()
    }

    fn exz_suffix(extz: &BitExtz) -> String {
        let result = &extz.result;
        let mut out = format!(
            " err={} thre={} ps={} pe={} ts={} te={} pl={} tl={} cn={} C",
            result.error,
            result.threshold,
            result.pattern_start,
            result.pattern_end,
            result.text_start,
            result.text_end,
            result.pattern_length,
            result.text_length,
            extz.cigar.packed.len()
        );
        for cell in &extz.cigar.packed {
            out.push_str(&format!(
                " {}{}",
                cell & 0x3fff,
                b"EMDI"[(cell >> 14) as usize] as char
            ));
        }
        out
    }

    fn run_scenario(scenario: &Scenario) -> Vec<String> {
        let mut reads = ReadStore::default();
        for (index, sequence) in scenario.reads.iter().enumerate() {
            reads.push_read(&format!("read{}", index), sequence);
        }

        let mut region = OverlapRegion {
            x_identifier: 0,
            x_position_start: scenario.x_start,
            x_position_end: scenario.x_end,
            y_identifier: 1,
            y_position_start: scenario.y_start,
            y_position_end: scenario.y_end,
            y_strand: scenario.y_strand,
            overlap_length: scenario.x_end + 1 - scenario.x_start,
            ..OverlapRegion::default()
        };
        region.fake_cigar = FakeCigar::default();
        for &(site, shift) in &scenario.gaps {
            region.fake_cigar.push(site, shift);
        }
        for window in &scenario.windows {
            region.window_cigar.append_window(
                window.x_start,
                window.x_end,
                window.y_start,
                window.y_end,
                window.error,
                0,
                0,
                0,
            );
        }

        let query = &scenario.reads[0];
        let query_total_length = query.len() as i64;
        let mut extz = BitExtz::new();
        let mut lines = Vec::new();

        for spec in &scenario.semis {
            let (start, end, auxiliary) = update_semi_coord(
                &region,
                &TargetSource::Reads(&reads),
                spec.query_start,
                spec.query_end,
                spec.threshold,
            );
            let aligned = i64::from(start != -1 && end != -1);
            lines.push(format!(
                "SEMI {} qs={} qe={} thre={} -> ok={} ts={} te={} ab={}",
                scenario.name,
                spec.query_start,
                spec.query_end,
                spec.threshold,
                aligned,
                start,
                end,
                auxiliary
            ));
        }

        for spec in &scenario.estimates {
            let value = cal_estimate_err(&region, 400, spec.query_start, spec.query_end, 0.04);
            lines.push(format!(
                "EST {} qs={} qe={} -> {}",
                scenario.name, spec.query_start, spec.query_end, value
            ));
        }

        for spec in &scenario.hc_estimates {
            let (value, exact) = cal_estimate_err_hc(
                &region,
                400,
                spec.query_start,
                spec.query_end,
                spec.target_start,
                spec.target_end,
                0.04,
            );
            lines.push(format!(
                "ESTHC {} qs={} qe={} ts={} te={} -> {} exact={}",
                scenario.name,
                spec.query_start,
                spec.query_end,
                spec.target_start,
                spec.target_end,
                value,
                i32::from(exact)
            ));
        }

        for spec in &scenario.exacts {
            let mut state = ExzAdvState {
                pattern_start: -1,
                pattern_end: -1,
                threshold: 0,
            };
            let aligned = cal_exact_exz(
                &region,
                &TargetSource::Reads(&reads),
                &mut extz,
                query,
                spec.query_start,
                spec.query_end,
                spec.target_start,
                spec.target_end,
                &mut state,
                spec.mode,
            );
            lines.push(format!(
                "EXACT {} qs={} qe={} ts={} te={} mode={} -> ret={} pts={} pte={}{}",
                scenario.name,
                spec.query_start,
                spec.query_end,
                spec.target_start,
                spec.target_end,
                spec.mode,
                i64::from(aligned),
                state.pattern_start,
                state.pattern_end,
                exz_suffix(&extz)
            ));
        }

        let mut state = ExzAdvState::default();
        for spec in &scenario.advs {
            if !spec.reuse_state {
                state = ExzAdvState::default();
            }
            let aligned = cal_exz_infi_adv(
                &region,
                &TargetSource::Reads(&reads),
                &mut extz,
                query,
                spec.query_start,
                spec.query_end,
                spec.target_start,
                spec.target_end,
                &mut state,
                spec.threshold,
                query_total_length,
                spec.mode,
            );
            lines.push(format!(
                "INFIA {} qs={} qe={} ts={} te={} thre={} mode={} reuse={} -> ret={} pts={} pte={} pthre={}{}",
                scenario.name,
                spec.query_start,
                spec.query_end,
                spec.target_start,
                spec.target_end,
                spec.threshold,
                spec.mode,
                i32::from(spec.reuse_state),
                i64::from(aligned),
                state.pattern_start,
                state.pattern_end,
                state.threshold,
                exz_suffix(&extz)
            ));
        }

        lines
    }

    fn scenarios() -> Vec<Scenario> {
        let mut state = Lcg(987654321);
        let query = state.sequence(2000);
        let prefix = state.sequence(600);
        let suffix = state.sequence(400);

        let mut exact_target = prefix.clone();
        exact_target.extend_from_slice(&query);
        exact_target.extend_from_slice(&suffix);

        let mut substituted = exact_target.clone();
        for index in (0..query.len()).step_by(97) {
            let position = prefix.len() + index;
            substituted[position] = if substituted[position] == b'A' { b'C' } else { b'A' };
        }

        let mut inserted = prefix.clone();
        inserted.extend_from_slice(&query[..700]);
        inserted.extend_from_slice(b"ACGTAC");
        inserted.extend_from_slice(&query[700..]);
        inserted.extend_from_slice(&suffix);

        let reversed = reverse_complement(&exact_target);

        let exact_windows = [
            WindowSpec { x_start: 200, x_end: 599, y_start: 800, y_end: 1199, error: 0 },
            WindowSpec { x_start: 600, x_end: 999, y_start: 1200, y_end: 1599, error: 0 },
            WindowSpec { x_start: 1000, x_end: 1399, y_start: 1600, y_end: 1999, error: 0 },
            WindowSpec { x_start: 1400, x_end: 1799, y_start: 2000, y_end: 2399, error: 0 },
        ];

        vec![
            Scenario {
                name: "A_exact",
                reads: vec![query.clone(), exact_target.clone()],
                x_start: 200,
                x_end: 1799,
                y_start: 800,
                y_end: 2399,
                y_strand: 0,
                gaps: vec![(200, 0), (1799, 0)],
                windows: exact_windows
                    .iter()
                    .map(|window| WindowSpec {
                        x_start: window.x_start,
                        x_end: window.x_end,
                        y_start: window.y_start,
                        y_end: window.y_end,
                        error: window.error,
                    })
                    .collect(),
                semis: vec![
                    SemiSpec { query_start: 200, query_end: 1800, threshold: 15 },
                    SemiSpec { query_start: 200, query_end: 400, threshold: 15 },
                    SemiSpec { query_start: 600, query_end: 1200, threshold: 15 },
                    SemiSpec { query_start: 200, query_end: 1800, threshold: 100 },
                    SemiSpec { query_start: 200, query_end: 1800, threshold: 0 },
                ],
                estimates: vec![
                    EstimateSpec { query_start: 200, query_end: 1800 },
                    EstimateSpec { query_start: 200, query_end: 400 },
                    EstimateSpec { query_start: 250, query_end: 900 },
                    EstimateSpec { query_start: 250, query_end: 850 },
                    EstimateSpec { query_start: 1000, query_end: 1000 },
                    EstimateSpec { query_start: 1900, query_end: 2100 },
                ],
                hc_estimates: vec![
                    HcEstimateSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400 },
                    HcEstimateSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2401 },
                    HcEstimateSpec { query_start: 200, query_end: 400, target_start: 800, target_end: 1000 },
                    HcEstimateSpec { query_start: 200, query_end: 1800, target_start: 799, target_end: 2399 },
                    HcEstimateSpec { query_start: 1900, query_end: 2100, target_start: 2500, target_end: 2700 },
                ],
                exacts: vec![
                    ExactSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, mode: 0 },
                    ExactSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, mode: 1 },
                    ExactSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, mode: 2 },
                    ExactSpec { query_start: 200, query_end: 1800, target_start: -1, target_end: -1, mode: 3 },
                    ExactSpec { query_start: 700, query_end: 1300, target_start: -1, target_end: -1, mode: 3 },
                    ExactSpec { query_start: 700, query_end: 1300, target_start: 1300, target_end: 1900, mode: 0 },
                    ExactSpec { query_start: 200, query_end: 1799, target_start: 800, target_end: 2400, mode: 0 },
                    ExactSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2399, mode: 0 },
                ],
                advs: vec![
                    AdvSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, threshold: 40, mode: 0, reuse_state: false },
                    AdvSpec { query_start: 200, query_end: 1800, target_start: -1, target_end: -1, threshold: 40, mode: 3, reuse_state: false },
                    AdvSpec { query_start: 200, query_end: 1000, target_start: 800, target_end: 1600, threshold: 40, mode: 1, reuse_state: false },
                    AdvSpec { query_start: 1000, query_end: 1800, target_start: 1600, target_end: 2400, threshold: 40, mode: 2, reuse_state: false },
                    AdvSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, threshold: 3, mode: 0, reuse_state: false },
                    AdvSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, threshold: 3, mode: 0, reuse_state: true },
                    AdvSpec { query_start: 700, query_end: 1300, target_start: -1, target_end: -1, threshold: 40, mode: 3, reuse_state: false },
                    AdvSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, threshold: 40, mode: 0, reuse_state: false },
                ],
            },
            Scenario {
                name: "B_subs",
                reads: vec![query.clone(), substituted],
                x_start: 200,
                x_end: 1799,
                y_start: 800,
                y_end: 2399,
                y_strand: 0,
                gaps: vec![(200, 0), (1799, 0)],
                windows: vec![
                    WindowSpec { x_start: 200, x_end: 599, y_start: 800, y_end: 1199, error: 4 },
                    WindowSpec { x_start: 600, x_end: 999, y_start: 1200, y_end: 1599, error: 4 },
                    WindowSpec { x_start: 1000, x_end: 1399, y_start: 1600, y_end: 1999, error: 4 },
                    WindowSpec { x_start: 1400, x_end: 1799, y_start: 2000, y_end: 2399, error: 5 },
                ],
                semis: vec![
                    SemiSpec { query_start: 200, query_end: 1800, threshold: 15 },
                    SemiSpec { query_start: 200, query_end: 1800, threshold: 100 },
                ],
                estimates: vec![
                    EstimateSpec { query_start: 200, query_end: 1800 },
                    EstimateSpec { query_start: 200, query_end: 1000 },
                    EstimateSpec { query_start: 250, query_end: 900 },
                    EstimateSpec { query_start: 250, query_end: 850 },
                ],
                hc_estimates: vec![
                    HcEstimateSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400 },
                    HcEstimateSpec { query_start: 200, query_end: 400, target_start: 800, target_end: 1000 },
                    HcEstimateSpec { query_start: 200, query_end: 800, target_start: 800, target_end: 1400 },
                ],
                exacts: vec![
                    ExactSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, mode: 0 },
                    ExactSpec { query_start: 200, query_end: 1800, target_start: -1, target_end: -1, mode: 3 },
                ],
                advs: vec![
                    AdvSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, threshold: 8, mode: 0, reuse_state: false },
                    AdvSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, threshold: 40, mode: 0, reuse_state: true },
                    AdvSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, threshold: 30, mode: 0, reuse_state: true },
                    AdvSpec { query_start: 200, query_end: 1800, target_start: -1, target_end: -1, threshold: 40, mode: 3, reuse_state: false },
                    AdvSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, threshold: 5, mode: 0, reuse_state: false },
                    AdvSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, threshold: 17, mode: 0, reuse_state: false },
                    AdvSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, threshold: 18, mode: 0, reuse_state: false },
                    AdvSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, threshold: 19, mode: 0, reuse_state: false },
                ],
            },
            Scenario {
                name: "C_ins6",
                reads: vec![query.clone(), inserted],
                x_start: 200,
                x_end: 1799,
                y_start: 800,
                y_end: 2405,
                y_strand: 0,
                gaps: vec![(200, 0), (700, 6), (1799, 6)],
                windows: vec![
                    WindowSpec { x_start: 200, x_end: 699, y_start: 800, y_end: 1299, error: 1 },
                    WindowSpec { x_start: 700, x_end: 1099, y_start: 1306, y_end: 1705, error: 1 },
                    WindowSpec { x_start: 1100, x_end: 1399, y_start: 1706, y_end: 2005, error: 1 },
                    WindowSpec { x_start: 1400, x_end: 1799, y_start: 2006, y_end: 2405, error: 1 },
                ],
                semis: vec![
                    SemiSpec { query_start: 200, query_end: 1800, threshold: 15 },
                    SemiSpec { query_start: 700, query_end: 1300, threshold: 15 },
                    SemiSpec { query_start: 200, query_end: 1800, threshold: 100 },
                ],
                estimates: vec![
                    EstimateSpec { query_start: 200, query_end: 1800 },
                    EstimateSpec { query_start: 700, query_end: 1300 },
                    EstimateSpec { query_start: 650, query_end: 750 },
                ],
                hc_estimates: vec![
                    HcEstimateSpec { query_start: 700, query_end: 1300, target_start: 1306, target_end: 1906 },
                    HcEstimateSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400 },
                    HcEstimateSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2406 },
                ],
                exacts: vec![
                    ExactSpec { query_start: 700, query_end: 1300, target_start: -1, target_end: -1, mode: 3 },
                    ExactSpec { query_start: 700, query_end: 1300, target_start: 1306, target_end: 1906, mode: 0 },
                    ExactSpec { query_start: 700, query_end: 1300, target_start: 1300, target_end: 1900, mode: 0 },
                    ExactSpec { query_start: 200, query_end: 1800, target_start: -1, target_end: -1, mode: 3 },
                ],
                advs: vec![
                    AdvSpec { query_start: 700, query_end: 1300, target_start: -1, target_end: -1, threshold: 20, mode: 3, reuse_state: false },
                    AdvSpec { query_start: 200, query_end: 1800, target_start: -1, target_end: -1, threshold: 20, mode: 3, reuse_state: false },
                    AdvSpec { query_start: 700, query_end: 1300, target_start: 1306, target_end: 1906, threshold: 8, mode: 0, reuse_state: false },
                    AdvSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, threshold: 20, mode: 0, reuse_state: false },
                    AdvSpec { query_start: 700, query_end: 1300, target_start: 1306, target_end: 1906, threshold: 5, mode: 0, reuse_state: false },
                    AdvSpec { query_start: 700, query_end: 1300, target_start: 1306, target_end: 1906, threshold: 6, mode: 0, reuse_state: false },
                ],
            },
            Scenario {
                name: "D_reverse",
                reads: vec![query.clone(), reversed],
                x_start: 200,
                x_end: 1799,
                y_start: 800,
                y_end: 2399,
                y_strand: 1,
                gaps: vec![(200, 0), (1799, 0)],
                windows: vec![
                    WindowSpec { x_start: 200, x_end: 599, y_start: 800, y_end: 1199, error: 0 },
                    WindowSpec { x_start: 600, x_end: 999, y_start: 1200, y_end: 1599, error: 0 },
                    WindowSpec { x_start: 1000, x_end: 1399, y_start: 1600, y_end: 1999, error: 0 },
                    WindowSpec { x_start: 1400, x_end: 1799, y_start: 2000, y_end: 2399, error: 0 },
                ],
                semis: vec![
                    SemiSpec { query_start: 200, query_end: 1800, threshold: 15 },
                    SemiSpec { query_start: 200, query_end: 400, threshold: 15 },
                    SemiSpec { query_start: 200, query_end: 1800, threshold: 100 },
                ],
                estimates: vec![
                    EstimateSpec { query_start: 200, query_end: 1800 },
                    EstimateSpec { query_start: 250, query_end: 900 },
                ],
                hc_estimates: vec![
                    HcEstimateSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400 },
                    HcEstimateSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2401 },
                ],
                exacts: vec![
                    ExactSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, mode: 0 },
                    ExactSpec { query_start: 200, query_end: 1800, target_start: -1, target_end: -1, mode: 3 },
                    ExactSpec { query_start: 700, query_end: 1300, target_start: -1, target_end: -1, mode: 3 },
                ],
                advs: vec![
                    AdvSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, threshold: 40, mode: 0, reuse_state: false },
                    AdvSpec { query_start: 200, query_end: 1800, target_start: -1, target_end: -1, threshold: 40, mode: 3, reuse_state: false },
                    AdvSpec { query_start: 700, query_end: 1300, target_start: -1, target_end: -1, threshold: 40, mode: 3, reuse_state: false },
                    AdvSpec { query_start: 200, query_end: 1000, target_start: 800, target_end: 1600, threshold: 40, mode: 1, reuse_state: false },
                    AdvSpec { query_start: 1000, query_end: 1800, target_start: 1600, target_end: 2400, threshold: 40, mode: 2, reuse_state: false },
                ],
            },
        ]
    }

    fn generated_lines() -> Vec<String> {
        let mut lines = Vec::new();

        const ERRORS: [i64; 4] = [1, 15, 31, 100];
        const STARTS: [i64; 9] = [-5, 0, 1, 10, 100, 500, 999, 1000, 1500];
        const TOTALS: [i64; 2] = [1000, 10000];
        const LENGTHS: [i64; 4] = [100, 130, 500, 2000];
        for &error in &ERRORS {
            for &start in &STARTS {
                for &total in &TOTALS {
                    for &length in &LENGTHS {
                        let mut auxiliary_begin = 0;
                        let mut auxiliary_end = 0;
                        let mut result_start = 0;
                        let mut result_length = 0;
                        let ok = match initialise_window_alignment_from_parts(
                            error,
                            start,
                            total,
                            length,
                            &mut auxiliary_begin,
                            &mut auxiliary_end,
                            &mut result_start,
                            &mut result_length,
                        ) {
                            true => 1,
                            false => 0,
                        };
                        lines.push(format!(
                            "INITWALN err={} s={} l={} wl={} -> ok={} ab={} ae={} rs={} rl={}",
                            error, start, total, length, ok, auxiliary_begin, auxiliary_end,
                            result_start, result_length
                        ));
                    }
                }
            }
        }

        let adjusts: [(i64, i64, i64, i64, i64, i64, i64, i64); 8] = [
            (20, 80, 500, 600, 100, 1000, 10, 1),
            (20, 80, 500, 600, 100, 1000, 10, 2),
            (0, 100, 900, 1000, 100, 1000, 5, 1),
            (0, 100, 900, 1000, 100, 1000, 5, 2),
            (10, 90, 5, 95, 100, 200, 0, 1),
            (10, 90, 5, 95, 100, 200, 0, 2),
            (90, 100, 100, 200, 100, 200, 30, 1),
            (0, 10, 0, 100, 100, 200, 30, 2),
        ];
        for &(original_query_start, original_query_end, original_target_start,
               original_target_end, query_length, target_length, threshold, mode) in &adjusts
        {
            let mut query_start = original_query_start;
            let mut query_end = original_query_end;
            let mut target_start = original_target_start;
            let mut target_end = original_target_end;
            adjust_ext_offset(
                &mut query_start,
                &mut query_end,
                &mut target_start,
                &mut target_end,
                query_length,
                target_length,
                threshold,
                mode,
            );
            lines.push(format!(
                "ADJUST mode={} qs={} qe={} ts={} te={} ql={} tl={} thre={} -> qs={} qe={} ts={} te={}",
                mode, original_query_start, original_query_end, original_target_start,
                original_target_end, query_length, target_length, threshold, query_start,
                query_end, target_start, target_end
            ));
        }

        for scenario in scenarios() {
            lines.extend(run_scenario(&scenario));
        }
        lines
    }

    /// The C test needs the `aux_end` output of `init_waln`, which no Rust
    /// caller uses, so the replay computes it from the public function.
    #[allow(clippy::too_many_arguments)]
    fn initialise_window_alignment_from_parts(
        threshold: i64,
        start: i64,
        target_total_length: i64,
        alignment_length: i64,
        auxiliary_begin: &mut i64,
        auxiliary_end: &mut i64,
        result_start: &mut i64,
        result_length: &mut i64,
    ) -> bool {
        *auxiliary_begin = -1;
        *auxiliary_end = -1;
        *result_start = -1;
        *result_length = -1;
        match initialise_window_alignment(threshold, start, target_total_length, alignment_length) {
            Some(window) => {
                *auxiliary_begin = window.auxiliary_begin;
                *auxiliary_end = window.auxiliary_end;
                *result_start = window.start;
                *result_length = window.length;
                true
            }
            None => false,
        }
    }

    #[test]
    #[ignore = "requires reference/ec_calexz_probe.out from the C++ probe"]
    fn cal_exz_layer_matches_the_reference() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("reference")
            .join("ec_calexz_probe.out");
        let contents = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("cannot read {}: {}", path.display(), error));
        let golden: Vec<&str> = contents.lines().collect();
        let generated = generated_lines();

        assert_eq!(
            generated.len(),
            golden.len(),
            "the probe printed {} lines, the replay generated {}",
            golden.len(),
            generated.len()
        );
        for (index, (expected, actual)) in golden.iter().zip(generated.iter()).enumerate() {
            assert_eq!(
                expected, actual,
                "line {} differs\nexpected: {}\nactual:   {}",
                index + 1, expected, actual
            );
        }
    }
}

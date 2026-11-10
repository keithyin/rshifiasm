//! The unitig-level refinement of the `trans_base_infer` path, the counterpart
//! of `ovlp_base_direct`, `rechain_aln`, `ul_raw_aln`, `ug_lalign` and
//! `split_ug_lalign` in `Correct.cpp` / `inter.cpp`.
//!
//! `gen_trans_adaptive_aln` (`inter.cpp` 10297) chains every unitig of the
//! polished graph against the global minimizer index and then aligns the
//! chains in two rounds through [`split_ug_lalign`]: the first pass aligns at
//! the optimistic error rate, the second at the fallback rate. With
//! `auxiliary` present the second pass re-aligns from scratch through
//! [`ul_raw_aln`], which aligns the chain blocks directly
//! ([`ovlp_base_direct`]), retries the blocks that failed by re-chaining the
//! covered windows ([`rechain_alignment_ul`]) and finally replaces the
//! overlap by its assembled windows.
//!
//! The overlap list is a fixed array in the reference; every function here
//! takes the slice plus the logical length (`ol->length`) and returns the new
//! length, so the elements past the logical end survive the shuffling of
//! [`split_ug_lalign`] exactly as they do in the reference.

use crate::correction::aln_adv::hc_aln_exz_adv;
use crate::correction::cal_exz::{
    adjust_ext_offset, push_unmapped_window, FORCED_SINGLE_LENGTH, MAXIMUM_SINGLE_ERROR,
    MAXIMUM_SINGLE_LENGTH, MODE_BACKWARD_EXTENSION, MODE_FORWARD_EXTENSION, MODE_GLOBAL,
};
use crate::correction::fast_cigar::return_target_chain;
use crate::correction::hc_align::{
    align_ul_ed_post_extz, gen_extend_err_exz, BitExtz, HcAlignParameters,
};
use crate::correction::recalculate::similarity_pass;
use crate::correction::region::{THRESHOLD_MAXIMUM, WINDOW_LENGTH};
use crate::correction::target::TargetSource;
use crate::correction::window_chain::{generate_window_chain, is_unaligned_window};
use crate::overlap::region::{CandidatesList, KmerHit, OverlapRegion};

/// `ul_ov_t` as far as the direct block pass consumes it: the bracket of the
/// chain plus the fixed interval ends.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ChainBracket {
    /// `qn`: the first chain hit of the bracket, `(uint32_t)-1` for an open
    /// left end.
    pub qn: u32,
    /// `tn`: the end of the bracket; `ch_n` (past the hits) extends to
    /// `qs` of the bracket.
    pub tn: u32,
    pub qs: u32,
    pub qe: u32,
    pub ts: u32,
    pub te: u32,
}

/// `ovlp_base_direct` (`Correct.cpp` 17367): aligns every block of the
/// bracket in one pass, without the primary-hit gate and the zero-error
/// shortcut of the Hi-C variant, and records unaligned windows for the retry.
///
/// The reference passes `rid` for debug logging only.
#[allow(clippy::too_many_arguments)]
pub(crate) fn ovlp_base_direct(
    region: &OverlapRegion,
    hits: &[KmerHit],
    hit_count: i64,
    bracket: &ChainBracket,
    window_length: i64,
    targets: &TargetSource,
    query: &[u8],
    extz: &mut BitExtz,
    output: &mut crate::correction::window::WindowCigar,
    error_rate: f64,
    query_total_length: i64,
    target_total_length: i64,
) {
    let first_block: i64 = if bracket.qn == u32::MAX {
        -1
    } else {
        i64::from(bracket.qn)
    };
    let last_block = i64::from(bracket.tn);
    debug_assert!(last_block >= first_block + 1);

    let mut previous = first_block;
    let mut index = first_block + 1;
    while index <= last_block {
        let mut query_start;
        let mut query_end;
        let mut target_start;
        let mut target_end;
        let mode;

        if previous >= 0 {
            query_start = i64::from(hits[previous as usize].self_offset);
            target_start = i64::from(hits[previous as usize].offset);
        } else {
            query_start = i64::from(bracket.qs);
            target_start = -1;
        }

        if index < hit_count {
            query_end = i64::from(hits[index as usize].self_offset);
            target_end = i64::from(hits[index as usize].offset);
        } else {
            query_end = i64::from(bracket.qe);
            target_end = -1;
        }

        if target_start != -1 && target_end != -1 {
            mode = MODE_GLOBAL; // global
        } else if target_start != -1 {
            mode = MODE_FORWARD_EXTENSION; // forward extension
        } else if target_end != -1 {
            mode = MODE_BACKWARD_EXTENSION; // backward extension
        } else {
            mode = 3; // no primary hit within the bracket
        }

        if mode == MODE_FORWARD_EXTENSION || mode == MODE_BACKWARD_EXTENSION {
            adjust_ext_offset(
                &mut query_start,
                &mut query_end,
                &mut target_start,
                &mut target_end,
                query_total_length,
                target_total_length,
                0,
                mode,
            );
        }

        let aligned = hc_aln_exz_adv(
            region,
            targets,
            extz,
            query,
            query_start,
            query_end,
            target_start,
            target_end,
            mode,
            window_length,
            query_total_length,
            error_rate,
            MAXIMUM_SINGLE_LENGTH,
            MAXIMUM_SINGLE_ERROR,
            FORCED_SINGLE_LENGTH,
            -1,
            output,
        );

        if !aligned {
            // `push_unmap_alnw` takes `[q0, q1]` and `[t0, t1]`, so an
            // exclusive end loses one base; an unset end stays `-1`/`-2`.
            push_unmapped_window(
                output,
                query_start as i32,
                query_end - 1,
                target_start,
                target_end - 1,
                mode,
            );
        }

        previous = index;
        index += 1;
    }
}

/// `rechain_aln` (`Correct.cpp` 17145): re-chains the covered windows of one
/// unaligned window and retries it once, at the `error_threshold` mode the
/// window records.
///
/// The retry is accepted only when the fresh alignment produced exactly one
/// new window whose fixed boundaries match the old unaligned one; the new
/// window then replaces the old entry in place. Unlike the Hi-C twin this
/// builds the chain bracket from the chain ends and calls
/// [`ovlp_base_direct`].
#[allow(clippy::too_many_arguments)]
pub(crate) fn rechain_alignment_ul(
    region: &mut OverlapRegion,
    candidates: &mut CandidatesList,
    output: &mut crate::correction::window::WindowCigar,
    auxiliary_index: usize,
    window_length: i64,
    targets: &TargetSource,
    query: &[u8],
    extz: &mut BitExtz,
    error_rate: f64,
    query_total_length: i64,
    target_total_length: i64,
    h_kmer_hit: i64,
) {
    let initial_length = candidates.length;
    let window = output.windows[auxiliary_index];
    let query_start = i64::from(window.x_start);
    let query_end = i64::from(window.x_end) + 1;
    let target_start = i64::from(window.y_start);
    let target_end = i64::from(window.y_end) + 1;
    if query_end - query_start < FORCED_SINGLE_LENGTH
        || target_end - target_start < FORCED_SINGLE_LENGTH
    {
        return;
    }

    let mode = i64::from(window.error_threshold);
    // The ultra-long caller passes `is_accurate = 0` to `gen_win_chain`.
    let chain_count = generate_window_chain(
        region,
        candidates,
        query_start,
        query_end,
        target_start,
        target_end,
        window_length,
        targets,
        query,
        extz,
        query_total_length,
        target_total_length,
        error_rate,
        h_kmer_hit,
        mode,
        false,
    );

    if chain_count != 0 {
        let base = initial_length as usize;
        let last = base + chain_count as usize - 1;
        let mut retry = true;
        let mut bracket = ChainBracket {
            qn: 0,
            tn: 0,
            qs: 0,
            qe: 0,
            ts: u32::MAX,
            te: u32::MAX,
        };
        match mode {
            MODE_GLOBAL => {
                debug_assert_eq!(i64::from(candidates.list[base].self_offset), query_start);
                debug_assert_eq!(i64::from(candidates.list[base].offset), target_start);
                debug_assert_eq!(i64::from(candidates.list[last].self_offset), query_end);
                debug_assert_eq!(i64::from(candidates.list[last].offset), target_end);
                bracket.qn = 0;
                bracket.tn = (chain_count - 1) as u32;
                bracket.qs = candidates.list[base].self_offset;
                bracket.ts = candidates.list[base].offset;
                bracket.qe = candidates.list[last].self_offset;
                bracket.te = candidates.list[last].offset;
                if chain_count <= 2 {
                    retry = false;
                }
            }
            MODE_FORWARD_EXTENSION => {
                debug_assert_eq!(i64::from(candidates.list[base].self_offset), query_start);
                debug_assert_eq!(i64::from(candidates.list[base].offset), target_start);
                bracket.qn = 0;
                bracket.tn = chain_count as u32;
                bracket.qs = candidates.list[base].self_offset;
                bracket.ts = candidates.list[base].offset;
                bracket.qe = query_total_length as u32;
                if chain_count <= 1 {
                    retry = false;
                }
            }
            MODE_BACKWARD_EXTENSION => {
                debug_assert_eq!(i64::from(candidates.list[last].self_offset), query_end);
                debug_assert_eq!(i64::from(candidates.list[last].offset), target_end);
                bracket.qn = u32::MAX;
                bracket.tn = (chain_count - 1) as u32;
                bracket.qs = 0;
                bracket.qe = candidates.list[last].self_offset;
                bracket.te = candidates.list[last].offset;
                if chain_count <= 1 {
                    retry = false;
                }
            }
            _ => {}
        }

        if retry {
            let previous_count = output.windows.len();
            ovlp_base_direct(
                region,
                &candidates.list[base..base + chain_count as usize],
                chain_count,
                &bracket,
                window_length,
                targets,
                query,
                extz,
                output,
                error_rate,
                query_total_length,
                target_total_length,
            );
            let new_count = output.windows.len();

            // The old unaligned window can be replaced only when the retry
            // produced exactly one aligned window.
            if new_count == previous_count + 1
                && !is_unaligned_window(&output.windows[new_count - 1])
            {
                let old_window = output.windows[auxiliary_index];
                let new_window = output.windows[new_count - 1];
                let query_start_matches = old_window.x_start == new_window.x_start;
                let query_end_matches = old_window.x_end == new_window.x_end;
                let target_start_matches = old_window.y_start == new_window.y_start;
                let target_end_matches = old_window.y_end == new_window.y_end;
                let replace = match mode {
                    MODE_GLOBAL => {
                        query_start_matches
                            && query_end_matches
                            && target_start_matches
                            && target_end_matches
                    }
                    MODE_FORWARD_EXTENSION => query_start_matches && target_start_matches,
                    MODE_BACKWARD_EXTENSION => query_end_matches && target_end_matches,
                    _ => false,
                };
                if replace {
                    output.windows[auxiliary_index] = new_window;
                    output.windows.pop();
                }
            }
        }
    }

    candidates.length = initial_length;
}

/// `ul_raw_aln` (`Correct.cpp` 24913): re-aligns one overlap from scratch
/// through its first candidate chain, retries the unaligned blocks and
/// replaces the overlap by the assembled windows.
///
/// Returns whether the assembled overlap is worth keeping.
#[allow(clippy::too_many_arguments)]
pub(crate) fn ul_raw_aln(
    region: &mut OverlapRegion,
    candidates: &mut CandidatesList,
    targets: &TargetSource,
    query: &[u8],
    extz: &mut BitExtz,
    error_rate: f64,
    window_length: i64,
    query_total_length: i64,
    h_kmer_hit: i64,
    auxiliary: &mut OverlapRegion,
) -> bool {
    return_target_chain(region, candidates);
    let chain_index = i64::from(region.shared_seed);
    if chain_index < 0 || chain_index >= candidates.length {
        return false;
    }
    let chain_start = chain_index as usize;
    let chain_read = candidates.list[chain_start].read_identifier;
    let mut scan = chain_start;
    while scan < candidates.length as usize && candidates.list[scan].read_identifier == chain_read
    {
        scan += 1;
    }
    let chain_count = scan as i64 - chain_index;
    let target_total_length = targets.length(region.y_identifier as usize);

    auxiliary.window_cigar.clear();
    auxiliary.y_identifier = region.y_identifier;
    auxiliary.y_strand = region.y_strand;
    auxiliary.x_position_start = region.x_position_start;
    auxiliary.x_position_end = region.x_position_end;
    auxiliary.y_position_start = region.y_position_start;
    auxiliary.y_position_end = region.y_position_end;

    // `ov.sec = 3; ov.ts = ov.te = -1; ov.qs = 0; ov.qn = -1; ov.qe = ql;
    // ov.tn = ch_n;` — one bracket spanning the whole overlap, open to the
    // left.
    let bracket = ChainBracket {
        qn: u32::MAX,
        tn: chain_count as u32,
        qs: 0,
        qe: query_total_length as u32,
        ts: u32::MAX,
        te: u32::MAX,
    };
    ovlp_base_direct(
        region,
        &candidates.list[chain_start..scan],
        chain_count,
        &bracket,
        window_length,
        targets,
        query,
        extz,
        &mut auxiliary.window_cigar,
        error_rate,
        query_total_length,
        target_total_length,
    );

    let initial_window_count = auxiliary.window_cigar.windows.len();
    for index in 0..initial_window_count {
        if !is_unaligned_window(&auxiliary.window_cigar.windows[index]) {
            continue;
        }
        // The retry may overwrite the candidate entries the first pass used;
        // `rechain_alignment_ul` restores the list length on return.
        rechain_alignment_ul(
            region,
            candidates,
            &mut auxiliary.window_cigar,
            index,
            window_length,
            targets,
            query,
            extz,
            error_rate,
            query_total_length,
            target_total_length,
            h_kmer_hit,
        );
    }
    // The retry appended windows at the end. Drop the initial windows that
    // are still unaligned, keep everything else, and restore the query order.
    if auxiliary.window_cigar.windows.len() > initial_window_count {
        let mut write = 0usize;
        for read in 0..auxiliary.window_cigar.windows.len() {
            if read < initial_window_count
                && is_unaligned_window(&auxiliary.window_cigar.windows[read])
            {
                continue;
            }
            auxiliary.window_cigar.windows[write] = auxiliary.window_cigar.windows[read];
            write += 1;
        }
        auxiliary.window_cigar.windows.truncate(write);
        // `radix_sort_window_list_xs_srt` is a stable 32-bit radix sort on
        // the signed `x_start`; window starts are never negative.
        crate::container::ksort::radix_sort_by_key_32(
            &mut auxiliary.window_cigar.windows,
            &|window: &crate::correction::window::WindowList| window.x_start as u32,
        );
    }

    /// `update_overlap_region` (`Correct.cpp` 17249): replaces the overlap by
    /// its assembled windows and widens both ends along the diagonal.
    region.update_from_windows(
        &auxiliary.window_cigar,
        query_total_length,
        target_total_length,
    );
    debug_assert!(region.x_position_end >= region.x_position_start);
    debug_assert!(region.y_position_end >= region.y_position_start);

    let mut total_error = 0i64;
    let overlap_length = i64::from(region.x_position_end) + 1 - i64::from(region.x_position_start);
    for window in &region.window_cigar.windows {
        if is_unaligned_window(window) {
            let query_extent =
                i64::from(window.x_end) + 1 - i64::from(window.x_start);
            let target_extent =
                i64::from(window.y_end) + 1 - i64::from(window.y_start);
            total_error += query_extent.max(target_extent);
        } else {
            total_error += i64::from(window.error);
        }
    }
    region.non_homopolymer_errors = total_error as u32;

    if total_error >= overlap_length {
        return false;
    }
    if total_error <= 0 && overlap_length > 0 {
        return true;
    }
    if total_error as f64 > overlap_length as f64 * error_rate {
        return false;
    }
    true
}

/// One overlap of the reference's `overlap_region_alloc` list.
pub type RegionList = Vec<OverlapRegion>;

/// `ug_lalign` (`Correct.cpp` 25850): refines every overlap of the list at
/// `error_rate`.
///
/// `auxiliary` selects the mode: `None` aligns the chain windows directly and
/// then measures the whole-overlap error (the first `split_ug_lalign` pass);
/// `Some` re-aligns each overlap from scratch through [`ul_raw_aln`] (the
/// second pass). Returns the new logical length; the elements past it are
/// kept for [`split_ug_lalign`]'s shuffling. The reference's
/// `mapped_overlaps_length` accumulator is not consumed downstream of this
/// call and is not carried over.
#[allow(clippy::too_many_arguments)]
fn ug_lalign(
    overlaps: &mut RegionList,
    length: usize,
    candidates: &mut CandidatesList,
    targets: &TargetSource,
    query: &[u8],
    query_total_length: i64,
    extz: &mut BitExtz,
    auxiliary: Option<&mut OverlapRegion>,
    error_rate: f64,
    window_length: i64,
    h_kmer_hit: i64,
    chain_cut: u64,
) -> usize {
    if length == 0 {
        return length;
    }

    let Some(auxiliary) = auxiliary else {
        // Direct pass: refine every chain's windows block by block.
        for region in overlaps.iter_mut().take(length) {
            let overlap_length =
                i64::from(region.x_position_end) + 1 - i64::from(region.x_position_start);
            // `z->shared_seed = z->non_homopolymer_errors` marks the chain
            // index inside the candidate list.
            region.shared_seed = region.non_homopolymer_errors as i32;
            let chain_index = region.shared_seed as usize;
            let chain_read = candidates.list[chain_index].read_identifier;
            let mut chain_hits = 0u64;
            let mut scan = chain_index;
            while scan < candidates.length as usize
                && candidates.list[scan].read_identifier == chain_read
                && chain_hits < chain_cut
            {
                scan += 1;
                chain_hits += 1;
            }
            if chain_hits < chain_cut {
                continue;
            }

            let accepted = align_ul_ed_post_extz(
                region,
                extz,
                targets,
                query,
                &HcAlignParameters {
                    error_rate,
                    window_length,
                    overlap_cutoff: -1.0,
                    force_alignment: false,
                },
            );
            if !accepted {
                continue;
            }
            // `uref` is always set on this path, so the secondary check is
            // on and the alignment fraction is ignored (`-1`).
            if similarity_pass(overlap_length, i64::from(region.alignment_length), true, -1.0, None)
            {
                region.is_match = 3;
            }
        }

        // Extension pass: measure the whole-overlap error including the
        // blocks between the windows.
        let maximum_error_rate = error_rate * 1.5 + 0.000001;
        let mut kept = 0usize;
        for index in 0..length {
            let extension = gen_extend_err_exz(
                &mut overlaps[index],
                targets,
                query,
                extz,
                window_length,
                -1.0,
                error_rate,
                maximum_error_rate,
                THRESHOLD_MAXIMUM,
                true,
            );
            overlaps[index].is_match = 0; // must be reset here
            if extension.error_rate <= error_rate {
                if kept != index {
                    overlaps.swap(kept, index);
                }
                overlaps[kept].is_match = 1;
                overlaps[kept].non_homopolymer_errors = extension.error_bases as u32;
                kept += 1;
            }
        }
        return kept;
    };

    // Second pass: re-align from scratch through the chain.
    let mut kept = 0usize;
    for index in 0..length {
        let accepted = ul_raw_aln(
            &mut overlaps[index],
            candidates,
            targets,
            query,
            extz,
            error_rate,
            window_length,
            query_total_length,
            h_kmer_hit,
            auxiliary,
        );
        if !accepted {
            overlaps[index].is_match = 0;
            continue;
        }
        if kept != index {
            overlaps.swap(kept, index);
        }
        overlaps[kept].is_match = 1;
        kept += 1;
    }
    kept
}

/// `split_ug_lalign` (`inter.cpp` 9967): aligns the first `head` overlaps at
/// `error_high` and the rest at `error_low`, restoring the reference's
/// element shuffling so both passes see exactly their own slice.
///
/// The block length is `min(THRESHOLD_MAX_SIZE / error_rate, WINDOW)`.
#[allow(clippy::too_many_arguments)]
pub fn split_ug_lalign(
    head: usize,
    overlaps: &mut RegionList,
    length: usize,
    error_high: f64,
    error_low: f64,
    candidates: &mut CandidatesList,
    targets: &TargetSource,
    query: &[u8],
    query_total_length: i64,
    extz: &mut BitExtz,
    mut auxiliary: Option<&mut OverlapRegion>,
    h_kmer_hit: i64,
    chain_cut: u64,
) -> usize {
    let mut tail = length - head;
    let mut head = head;
    if head > 0 {
        let original_length = length;
        let window_length = block_length(error_high);
        let new_head = ug_lalign(
            overlaps,
            head,
            candidates,
            targets,
            query,
            query_total_length,
            extz,
            auxiliary.as_deref_mut(),
            error_high,
            window_length,
            h_kmer_hit,
            chain_cut,
        );
        // Move the untouched tail entries behind the compacted prefix.
        let mut write = new_head;
        for index in head..original_length {
            if index != write {
                overlaps.swap(index, write);
            }
            write += 1;
        }
        head = new_head;
        let new_length = write;
        tail = new_length - head;
    }

    if tail > 0 {
        let original_length = head + tail;
        let window_length = block_length(error_low);
        // Swap the head and tail slices so the tail pass sees [0, tail).
        if head > 0 {
            let mut write = 0usize;
            for index in head..original_length {
                if index != write {
                    overlaps.swap(index, write);
                }
                write += 1;
            }
        }
        let new_tail = ug_lalign(
            overlaps,
            tail,
            candidates,
            targets,
            query,
            query_total_length,
            extz,
            auxiliary.as_deref_mut(),
            error_low,
            window_length,
            h_kmer_hit,
            chain_cut,
        );
        // Move the head entries behind the compacted tail prefix.
        let mut write = new_tail;
        for index in tail..original_length {
            if index != write {
                overlaps.swap(index, write);
            }
            write += 1;
        }
        let new_length = write;
        tail = new_tail;
        head = new_length - tail;
        // Swap back so the head occupies the front.
        if tail > 0 {
            let mut write = 0usize;
            for index in tail..new_length {
                if index != write {
                    overlaps.swap(index, write);
                }
                write += 1;
            }
        }
    }
    head
}

/// `wl = MIN(THRESHOLD_MAX_SIZE / erate, WINDOW)` of the reference: a double
/// quotient truncated into the `int64_t` block length.
fn block_length(error_rate: f64) -> i64 {
    ((THRESHOLD_MAXIMUM as f64 / error_rate).min(WINDOW_LENGTH as f64)) as i64
}

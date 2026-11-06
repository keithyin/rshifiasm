//! `gen_hc_r_alin`: the whole-overlap refinement pass of the new HiFi
//! correction path, plus the two terminal filters it applies before the gap
//! reassignment, `ff_lunalign` and `ff_tend`.
//!
//! The pass walks the overlap list of one query read. Every overlap is
//! refined window by window ([`align_hc_ed_post_extz`]), its residual error is
//! summed over the whole span ([`gen_extend_err_exz`]), its candidate chain is
//! rechecked and turned into a window cigar ([`generate_fast_cigar`]), and the
//! large-gap and terminal-error filters run. Survivors are compacted to the
//! front of the list and get `is_match = 1`; the failed overlaps stay behind
//! them in the permutation the reference's in-place swap leaves, and the list
//! is reported shrunk by the returned length.
//!
//! `move_gap` is the reference's twelve-argument `buf` slot: the function is
//! declared with it but neither reads it nor forwards it, so it is not part of
//! the Rust signature.

use crate::container::ksort::radix_sort_u64;
use crate::correction::cal_exz::{push_alignment_window, set_exact_exz};
use crate::correction::fast_cigar::generate_fast_cigar;
use crate::correction::gap_reassign::{reassign_gaps, TargetSource};
use crate::correction::hc_align::{
    BitExtz, HcAlignParameters, align_hc_ed_post_extz, gen_extend_err_exz,
};
use crate::correction::region::THRESHOLD_MAXIMUM;
use crate::correction::window::{CELL_LENGTH_MASK, WindowCigar};
use crate::correction::window_chain::is_unaligned_window;
use crate::overlap::paf::PafHit;
use crate::overlap::region::{CandidatesList, OverlapRegion};
use crate::reads::read_store::ReadStore;

/// Fraction of the overlap that must align, `OVERLAP_THRESHOLD_HIFI_FILTER`.
pub const OVERLAP_THRESHOLD_HIFI_FILTER: f64 = 0.9;

/// Window size the terminal filter measures errors in, `ff_tend`'s `wn`.
pub const TERMINAL_WINDOW_LENGTH: i64 = 384;

/// Distance from the read end the terminal filter inspects, `ff_tend`'s `dn`.
pub const TERMINAL_DISTANCE: i64 = 2000;

/// Fraction of the overlap [`TERMINAL_DISTANCE`] is clamped to, `ff_tend`'s
/// `dr`.
pub const TERMINAL_DISTANCE_RATE: f64 = 0.1;

/// Upper bound on the terminal error rate, the `0.36` of `ff_tend`'s `er`.
pub const TERMINAL_ERROR_RATE_CAP: f64 = 0.36;

/// Error count below which the terminal filter never triggers, `ff_tend`'s
/// `min_err`.
pub const TERMINAL_MINIMUM_ERRORS: i64 = 128;

/// Cigar operation for a matching base, `EAC_M`.
const EXACT_MATCH: u8 = 0;

/// Cigar operation for a target base the query does not have, `MOR_YP`.
const MORE_TARGET: u8 = 2;

/// `gen_hc_r_alin`: refines every overlap of one read.
///
/// `overlaps` is the read's whole overlap list, the reference's `ol->list` up
/// to `ol->length`. The returned value is the new `ol->length`: the first that
/// many entries are the survivors of the pass (`is_match = 1`, `strong` and
/// `without_large_indel` cleared) and the rest are the dropped overlaps in the
/// order the in-place swap leaves them.
///
/// The list is not truncated here because the reference does not clear the
/// dropped entries either; callers shrink it to the returned length.
#[allow(clippy::too_many_arguments)]
pub fn gen_hc_r_alin(
    overlaps: &mut [OverlapRegion],
    candidates: &mut CandidatesList,
    reads: &ReadStore,
    query: &[u8],
    extz: &mut BitExtz,
    auxiliary: &mut OverlapRegion,
    error_rate: f64,
    window_length: i64,
    h_kmer_hit: i64,
    scratch: &mut Vec<u16>,
    chemical_drop: bool,
    alignment_gap_rate: f64,
    alignment_gap_max: i64,
) -> usize {
    refine_overlaps(
        overlaps,
        candidates,
        reads,
        query,
        extz,
        auxiliary,
        error_rate,
        window_length,
        h_kmer_hit,
        scratch,
        chemical_drop,
        alignment_gap_rate,
        alignment_gap_max,
        false,
    )
}

/// `gen_hc_r_alin_nec`: the same refinement, but overlaps already marked
/// `is_match == 1` are kept as they are and only the rest is aligned.
///
/// The reference sets the flag itself for the exact hits it recovered from the
/// stored overlap list, so the exact overlaps never pay for a realignment.
/// Every entry still reaches the swap-and-mark tail of the loop, matching the
/// reference's unconditional `k++`.
#[allow(clippy::too_many_arguments)]
pub fn gen_hc_r_alin_nec(
    overlaps: &mut [OverlapRegion],
    candidates: &mut CandidatesList,
    reads: &ReadStore,
    query: &[u8],
    extz: &mut BitExtz,
    auxiliary: &mut OverlapRegion,
    error_rate: f64,
    window_length: i64,
    h_kmer_hit: i64,
    scratch: &mut Vec<u16>,
    chemical_drop: bool,
    alignment_gap_rate: f64,
    alignment_gap_max: i64,
) -> usize {
    refine_overlaps(
        overlaps,
        candidates,
        reads,
        query,
        extz,
        auxiliary,
        error_rate,
        window_length,
        h_kmer_hit,
        scratch,
        chemical_drop,
        alignment_gap_rate,
        alignment_gap_max,
        true,
    )
}

/// The shared loop of `gen_hc_r_alin` and `gen_hc_r_alin_nec`.
///
/// `skip_matched` is the only difference between the two: when set, an
/// overlap whose `is_match` is already one skips alignment, the filters and
/// the gap reassignment, and goes straight to the compaction tail.
#[allow(clippy::too_many_arguments)]
fn refine_overlaps(
    overlaps: &mut [OverlapRegion],
    candidates: &mut CandidatesList,
    reads: &ReadStore,
    query: &[u8],
    extz: &mut BitExtz,
    auxiliary: &mut OverlapRegion,
    error_rate: f64,
    window_length: i64,
    h_kmer_hit: i64,
    scratch: &mut Vec<u16>,
    chemical_drop: bool,
    alignment_gap_rate: f64,
    alignment_gap_max: i64,
    skip_matched: bool,
) -> usize {
    if overlaps.is_empty() {
        return 0;
    }

    let maximum_error_rate = error_rate * 1.5 + 0.000001;
    let parameters = HcAlignParameters {
        error_rate,
        window_length,
        overlap_cutoff: OVERLAP_THRESHOLD_HIFI_FILTER,
        force_alignment: false,
    };
    let query_length = query.len() as i64;

    let mut kept = 0usize;
    for index in 0..overlaps.len() {
        if !(skip_matched && overlaps[index].is_match == 1) {
            // The chain index of this overlap is parked in the error field so
            // the post-chain passes, which only receive the overlap, can find
            // it. It is overwritten by the measured error count below.
            overlaps[index].shared_seed = overlaps[index].non_homopolymer_errors as i32;

            if !align_hc_ed_post_extz(&mut overlaps[index], extz, reads, query, &parameters) {
                continue;
            }

            let extension = gen_extend_err_exz(
                &mut overlaps[index],
                reads,
                query,
                extz,
                window_length,
                -1.0,
                error_rate,
                maximum_error_rate,
                THRESHOLD_MAXIMUM,
                false,
            );
            overlaps[index].is_match = 0;
            if extension.error_rate > error_rate {
                continue;
            }
            overlaps[index].non_homopolymer_errors = extension.error_bases as u32;

            // The reference ignores the alignment's return value (it is always
            // one) and the error count it reports here; the count that
            // survives is the one `gen_extend_err_exz` measured above.
            generate_fast_cigar(
                &mut overlaps[index],
                candidates,
                window_length,
                reads,
                query,
                extz,
                auxiliary,
                error_rate,
                query_length,
                h_kmer_hit,
            );

            if alignment_gap_max >= 0
                && !large_gap_pass(
                    &overlaps[index],
                    error_rate,
                    alignment_gap_rate,
                    alignment_gap_max,
                )
            {
                continue;
            }

            if chemical_drop
                && terminal_error_trigger(
                    &overlaps[index],
                    TERMINAL_WINDOW_LENGTH,
                    TERMINAL_DISTANCE,
                    TERMINAL_DISTANCE_RATE,
                    if error_rate * 10.0 < TERMINAL_ERROR_RATE_CAP {
                        error_rate * 10.0
                    } else {
                        TERMINAL_ERROR_RATE_CAP
                    },
                    TERMINAL_MINIMUM_ERRORS,
                )
            {
                continue;
            }

            reassign_gaps(
                &mut overlaps[index],
                auxiliary,
                query,
                query_length,
                TargetSource::Reads(reads),
                -1,
                scratch,
            );
        }

        if kept != index {
            overlaps.swap(kept, index);
        }
        overlaps[kept].is_match = 1;
        overlaps[kept].strong = 0;
        overlaps[kept].without_large_indel = 0;
        kept += 1;
    }
    kept
}

/// `gen_hc_r_alin_ea`: refines the overlaps, but first recovers every overlap
/// that the previous correction round already stored as an exact match.
///
/// `exact_overlaps` is the read's stored overlap list (`R_INF.paf[i]`). The
/// entries whose `exact_region` is set are sorted by `(target, strand,
/// position-in-list)` and each overlap of the current chain list is looked up
/// in them. When the current overlap covers exactly the stored exact region
/// and the two reads agree base for base over that region, the overlap is
/// marked matched with a single exact window without any alignment. Everything
/// not matched that way goes through [`gen_hc_r_alin_nec`], which only aligns
/// the entries still missing the flag.
///
/// `exact_order` is the reference's `srt` scratch: on return it holds the
/// sorted lookup keys, exactly the first `srt->n` entries the reference leaves
/// behind. The returned value is the new overlap count, as in
/// [`gen_hc_r_alin`].
#[allow(clippy::too_many_arguments)]
pub fn gen_hc_r_alin_ea(
    overlaps: &mut [OverlapRegion],
    candidates: &mut CandidatesList,
    reads: &ReadStore,
    query: &[u8],
    extz: &mut BitExtz,
    auxiliary: &mut OverlapRegion,
    exact_overlaps: &[PafHit],
    error_rate: f64,
    window_length: i64,
    h_kmer_hit: i64,
    exact_order: &mut Vec<u64>,
    scratch: &mut Vec<u16>,
    chemical_drop: bool,
    alignment_gap_rate: f64,
    alignment_gap_max: i64,
) -> usize {
    if overlaps.is_empty() {
        return 0;
    }

    exact_order.clear();
    for (position, hit) in exact_overlaps.iter().enumerate() {
        if hit.exact_region {
            exact_order.push(
                ((u64::from(hit.target_identifier) << 1) | u64::from(hit.reverse as u8)) << 32
                    | position as u64,
            );
        }
    }

    if exact_order.is_empty() {
        return gen_hc_r_alin(
            overlaps,
            candidates,
            reads,
            query,
            extz,
            auxiliary,
            error_rate,
            window_length,
            h_kmer_hit,
            scratch,
            chemical_drop,
            alignment_gap_rate,
            alignment_gap_max,
        );
    }

    // The current overlaps, keyed the same way. The low half keeps the index
    // into `overlaps`, so the sort only orders the lookup, never the list.
    let mut ordered_overlaps: Vec<u64> = (0..overlaps.len())
        .map(|index| {
            let region = &mut overlaps[index];
            region.is_match = 0;
            region.strong = 0;
            region.without_large_indel = 0;
            ((u64::from(region.y_identifier) << 1) | u64::from(region.y_strand)) << 32
                | index as u64
        })
        .collect();

    radix_sort_u64(exact_order);
    radix_sort_u64(&mut ordered_overlaps);

    let mut exact_count = 0usize;
    let mut lookup = 0usize;
    for entry in &ordered_overlaps {
        let index = (*entry & 0xffff_ffff) as usize;
        let region = &mut overlaps[index];
        let target_key =
            (u64::from(region.y_identifier) << 1) | u64::from(region.y_strand);
        while lookup < exact_order.len() && (exact_order[lookup] >> 32) < target_key {
            lookup += 1;
        }
        if lookup >= exact_order.len() || (exact_order[lookup] >> 32) != target_key {
            continue;
        }

        let hit = &exact_overlaps[(exact_order[lookup] & 0xffff_ffff) as usize];
        if region.x_position_start != hit.query_start
            || region.x_position_end.wrapping_add(1) != hit.query_end
            || region.y_position_start != hit.target_start
            || region.y_position_end.wrapping_add(1) != hit.target_end
        {
            continue;
        }

        let target_length = hit.target_end.wrapping_sub(hit.target_start) as i64;
        let target = reads.read_sub_region(
            hit.target_identifier as usize,
            i64::from(hit.target_start),
            target_length,
            region.y_strand as u8,
        );
        if !exact_ec_check(query, hit.query_start, hit.query_end, &target) {
            continue;
        }

        region.is_match = 1;
        // The chain index moves to the shared-seed field for the post-chain
        // passes, exactly as the alignment path does.
        region.shared_seed = region.non_homopolymer_errors as i32;
        region.non_homopolymer_errors = 0;
        region.strong = 0;
        region.without_large_indel = 0;
        set_exact_exz(
            extz,
            i64::from(region.x_position_start),
            i64::from(region.x_position_end) + 1,
            i64::from(region.y_position_start),
            i64::from(region.y_position_end) + 1,
        );
        push_alignment_window(&mut region.window_cigar, extz);
        exact_count += 1;
    }

    if overlaps.len() > exact_count {
        return gen_hc_r_alin_nec(
            overlaps,
            candidates,
            reads,
            query,
            extz,
            auxiliary,
            error_rate,
            window_length,
            h_kmer_hit,
            scratch,
            chemical_drop,
            alignment_gap_rate,
            alignment_gap_max,
        );
    }

    overlaps.len()
}

/// `exact_ec_check`: whether the stored exact region of a hit matches the
/// query base for base. The reference's `ql`/`tl` arguments are unused.
pub(crate) fn exact_ec_check(query: &[u8], query_start: u32, query_end: u32, target: &[u8]) -> bool {
    let length = query_end.wrapping_sub(query_start) as usize;
    if length != target.len() {
        return false;
    }
    if length == 0 {
        return true;
    }
    let start = query_start as usize;
    if start + length > query.len() {
        return false;
    }
    query[start..start + length] == target[..]
}

/// `ff_lunalign`: rejects an overlap whose windows leave a large unaligned
/// gap.
///
/// The error and the gap are measured on the query and on the target side;
/// either side exceeding `maximum_gap` or `gap_rate` of the overlap drops it.
/// The final error check adds the gap to the window errors, because an
/// unaligned block counts as error too.
///
/// Returns whether the overlap passes the filter, matching the reference's
/// `1`/`0`.
pub fn large_gap_pass(
    region: &OverlapRegion,
    error_rate: f64,
    gap_rate: f64,
    maximum_gap: i64,
) -> bool {
    let windows = &region.window_cigar.windows;
    let x_start = i64::from(region.x_position_start);
    let x_end = i64::from(region.x_position_end);
    let y_start = i64::from(region.y_position_start);
    let y_end = i64::from(region.y_position_end);

    // One aligned window covering the whole overlap exactly: nothing to check.
    if windows.len() == 1
        && !is_unaligned_window(&windows[0])
        && i64::from(windows[0].x_start) == x_start
        && i64::from(windows[0].x_end) == x_end
        && i64::from(windows[0].y_start) == y_start
        && i64::from(windows[0].y_end) == y_end
    {
        return true;
    }

    let mut total_error = 0i64;
    let mut total_gap = 0i64;
    let mut query_position = x_start;
    let mut target_position = y_start;
    for window in windows {
        if is_unaligned_window(window) {
            continue;
        }
        let window_query_start = i64::from(window.x_start);
        let window_target_start = i64::from(window.y_start);
        if window_query_start != query_position {
            total_gap += (window_query_start - query_position).abs();
        }
        if window_target_start != target_position {
            total_gap += (window_target_start - target_position).abs();
        }
        query_position = i64::from(window.x_end) + 1;
        target_position = i64::from(window.y_end) + 1;
        total_error += i64::from(window.error);
    }

    let query_end = x_end + 1;
    let target_end = y_end + 1;
    if query_end != query_position {
        total_gap += (query_end - query_position).abs();
    }
    if target_end != target_position {
        total_gap += (target_end - target_position).abs();
    }

    if total_gap == 0 {
        return true;
    }
    if total_gap > maximum_gap {
        return false;
    }

    let query_length = x_end + 1 - x_start;
    let target_length = y_end + 1 - y_start;
    if total_gap as f64 > query_length as f64 * gap_rate
        || total_gap as f64 > target_length as f64 * gap_rate
    {
        return false;
    }

    total_error += total_gap;
    if total_error as f64 > query_length as f64 * error_rate
        || total_error as f64 > target_length as f64 * error_rate
    {
        return false;
    }
    true
}

/// `ff_tend`: whether either end of the overlap is mostly unaligned.
///
/// Both ends are scanned over the windows that start (or end) inside the
/// terminal region of length `distance`. Aligned windows are walked through
/// their cigar cells: matching bases count as covering the region, every other
/// operation and every unaligned window counts as error. The filter triggers
/// as soon as the accumulated error passes `minimum_errors`, is wider than
/// `window_threshold`, and dominates the bases seen in the region by
/// `error_rate`.
pub fn terminal_error_trigger(
    region: &OverlapRegion,
    window_threshold: i64,
    distance: i64,
    distance_rate: f64,
    error_rate: f64,
    minimum_errors: i64,
) -> bool {
    let windows = &region.window_cigar.windows;
    let window_count = windows.len() as i64;
    let query_start = i64::from(region.x_position_start);
    let query_end = i64::from(region.x_position_end.wrapping_add(1));
    let overlap_length = query_end - query_start;

    if overlap_length < window_threshold {
        return false;
    }
    let mut distance = distance;
    if distance as f64 > overlap_length as f64 * distance_rate {
        distance = (overlap_length as f64 * distance_rate) as i64;
    }
    if distance < window_threshold {
        return false;
    }
    if overlap_length < distance {
        return false;
    }

    // Forward scan of the left terminal region.
    let region_start = query_start;
    let region_end = query_start + distance;
    let mut error = 0i64;
    let mut matched = 0i64;
    let mut window_index = 0i64;
    while window_index < window_count
        && i64::from(windows[window_index as usize].x_start) < region_end
    {
        let window = windows[window_index as usize];
        if !is_unaligned_window(&window) {
            let mut cell_index = 0usize;
            let mut query_position = i64::from(window.x_start);
            while cell_index < window.cigar_length as usize && query_position < region_end {
                let run_start = query_position;
                let (next, operation, length) = pop_trace(
                    &region.window_cigar,
                    window_index as usize,
                    cell_index,
                );
                cell_index = next;
                if operation != MORE_TARGET {
                    query_position += i64::from(length);
                }
                let run_end = query_position;

                if operation == EXACT_MATCH {
                    let overlap_start = region_start.max(run_start);
                    let overlap_end = region_end.min(run_end);
                    if overlap_end > overlap_start {
                        matched += overlap_end - overlap_start;
                    }
                } else {
                    error += i64::from(length);
                }
                if terminal_triggered(
                    error,
                    matched,
                    window_threshold,
                    error_rate,
                    minimum_errors,
                ) {
                    return true;
                }
            }
        } else {
            let run_start = i64::from(window.x_start);
            let run_end = i64::from(window.x_end) + 1;
            error += run_end - run_start;
            if terminal_triggered(
                error,
                matched,
                window_threshold,
                error_rate,
                minimum_errors,
            ) {
                return true;
            }
        }
        window_index += 1;
    }

    // Backward scan of the right terminal region.
    let region_start = query_end - distance;
    let region_end = query_end;
    error = 0;
    matched = 0;
    let mut window_index = window_count - 1;
    while window_index >= 0
        && i64::from(windows[window_index as usize].x_end) + 1 > region_start
    {
        let window = windows[window_index as usize];
        if !is_unaligned_window(&window) {
            let mut cell_index = window.cigar_length as isize - 1;
            let mut query_position = i64::from(window.x_end) + 1;
            while cell_index >= 0 && query_position > region_start {
                let run_end = query_position;
                let (next, operation, length) = pop_trace_back(
                    &region.window_cigar,
                    window_index as usize,
                    cell_index,
                );
                cell_index = next;
                if operation != MORE_TARGET {
                    query_position -= i64::from(length);
                }
                let run_start = query_position;

                if operation == EXACT_MATCH {
                    let overlap_start = region_start.max(run_start);
                    let overlap_end = region_end.min(run_end);
                    if overlap_end > overlap_start {
                        matched += overlap_end - overlap_start;
                    }
                } else {
                    error += i64::from(length);
                }
                if terminal_triggered(
                    error,
                    matched,
                    window_threshold,
                    error_rate,
                    minimum_errors,
                ) {
                    return true;
                }
            }
        } else {
            let run_start = i64::from(window.x_start);
            let run_end = i64::from(window.x_end) + 1;
            error += run_end - run_start;
            if terminal_triggered(
                error,
                matched,
                window_threshold,
                error_rate,
                minimum_errors,
            ) {
                return true;
            }
        }
        window_index -= 1;
    }

    false
}

/// The trigger condition shared by the four scans of `ff_tend`.
fn terminal_triggered(
    error: i64,
    matched: i64,
    window_threshold: i64,
    error_rate: f64,
    minimum_errors: i64,
) -> bool {
    error > minimum_errors
        && matched + error > window_threshold
        && (error as f64) > ((matched + error) as f64 * error_rate)
}

/// `pop_trace`: coalesces the cell at `offset` with the following cells that
/// carry the same operation. Returns the offset just past the run, together
/// with the operation and the summed length.
fn pop_trace(cigar: &WindowCigar, window_index: usize, mut offset: usize) -> (usize, u8, u32) {
    let window = cigar.windows[window_index];
    let start = window.cigar_start as usize;
    let operation = (cigar.cells[start + offset] >> 14) as u8;
    let mut length = u32::from(cigar.cells[start + offset] & CELL_LENGTH_MASK);
    offset += 1;
    while offset < window.cigar_length as usize {
        let cell = cigar.cells[start + offset];
        if (cell >> 14) as u8 != operation {
            break;
        }
        length += u32::from(cell & CELL_LENGTH_MASK);
        offset += 1;
    }
    (offset, operation, length)
}

/// `pop_trace_back`: the backward counterpart of [`pop_trace`]. Returns the
/// offset just before the run, or `-1` when the run starts the window.
fn pop_trace_back(
    cigar: &WindowCigar,
    window_index: usize,
    mut offset: isize,
) -> (isize, u8, u32) {
    let window = cigar.windows[window_index];
    let start = window.cigar_start as isize;
    let operation = (cigar.cells[(start + offset) as usize] >> 14) as u8;
    let mut length = u32::from(cigar.cells[(start + offset) as usize] & CELL_LENGTH_MASK);
    offset -= 1;
    while offset >= 0 {
        let cell = cigar.cells[(start + offset) as usize];
        if (cell >> 14) as u8 != operation {
            break;
        }
        length += u32::from(cell & CELL_LENGTH_MASK);
        offset -= 1;
    }
    (offset, operation, length)
}

#[cfg(test)]
mod probe_goldens {
    use super::*;
    use crate::overlap::region::KmerHit;

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

    /// `print_windows`: the reference zeroes the extra boundaries of every
    /// aligned window before printing, so the golden carries those zeros.
    fn push_windows(lines: &mut Vec<String>, tag: &str, cigar: &WindowCigar) {
        lines.push(format!(
            "{} n={} cn={}",
            tag,
            cigar.windows.len(),
            cigar.cells.len()
        ));
        for (index, window) in cigar.windows.iter().enumerate() {
            let (extra_begin, extra_end) = if window.error != i16::MAX {
                (0, 0)
            } else {
                (window.extra_begin, window.extra_end)
            };
            let mut line = format!(
                " {}W {} x=[{},{}] y=[{},{}] eb={} ee={} err={} thre={} cidx={} clen={} :",
                tag,
                index,
                window.x_start,
                window.x_end,
                window.y_start,
                window.y_end,
                extra_begin,
                extra_end,
                window.error,
                window.error_threshold,
                window.cigar_start,
                window.cigar_length
            );
            for offset in 0..window.cigar_length {
                let cell = cigar.cells[(window.cigar_start + offset) as usize];
                line.push_str(&format!(
                    " {}{}",
                    cell & CELL_LENGTH_MASK,
                    b"EMDI"[(cell >> 14) as usize] as char
                ));
            }
            lines.push(line);
        }
    }

    struct OverlapSpec {
        y_identifier: u32,
        x_start: u32,
        x_end: u32,
        y_start: u32,
        y_end: u32,
        strand: u32,
        hits: Vec<(u32, u32)>,
    }

    struct BatchSpec {
        name: &'static str,
        reads: Vec<Vec<u8>>,
        overlaps: Vec<OverlapSpec>,
        filler_hits: u32,
        error_rate: f64,
        window_length: i64,
        h_kmer_hit: i64,
        chemical_drop: bool,
        alignment_gap_rate: f64,
        alignment_gap_max: i64,
    }

    /// `setup_overlap`: an overlap exactly as the C++ probe builds it, with a
    /// fake cigar spanning the whole query interval.
    fn build_overlap(spec: &OverlapSpec, non_homopolymer_errors: u32) -> OverlapRegion {
        let mut region = OverlapRegion {
            x_identifier: 0,
            x_position_start: spec.x_start,
            x_position_end: spec.x_end,
            x_strand: 0,
            y_identifier: spec.y_identifier,
            y_position_start: spec.y_start,
            y_position_end: spec.y_end,
            y_strand: spec.strand,
            overlap_length: spec.x_end + 1 - spec.x_start,
            alignment_length: spec.x_end + 1 - spec.x_start,
            non_homopolymer_errors,
            ..OverlapRegion::default()
        };
        region.fake_cigar.push(200, 0);
        region.fake_cigar.push(1799, 0);
        region
    }

    fn run_batch(batch: &BatchSpec, lines: &mut Vec<String>) {
        let mut reads = ReadStore::default();
        for (index, sequence) in batch.reads.iter().enumerate() {
            reads.push_read(&format!("read{}", index), sequence);
        }

        let mut candidates = CandidatesList::default();
        for index in 0..batch.filler_hits {
            candidates.push(KmerHit {
                read_identifier: 4096,
                reverse_strand: false,
                offset: 5000 + index,
                self_offset: 5000 + index,
                count: 1,
            });
        }

        let mut overlaps = Vec::new();
        for spec in &batch.overlaps {
            let region = build_overlap(spec, candidates.length as u32);
            for &(self_offset, offset) in &spec.hits {
                candidates.push(KmerHit {
                    read_identifier: spec.y_identifier,
                    reverse_strand: false,
                    offset,
                    self_offset,
                    count: 1,
                });
            }
            overlaps.push(region);
        }
        let original_length = overlaps.len();

        let mut extz = BitExtz::new();
        let mut auxiliary = OverlapRegion::default();
        let mut scratch = Vec::new();
        let query = batch.reads[0].clone();
        let kept = gen_hc_r_alin(
            &mut overlaps,
            &mut candidates,
            &reads,
            &query,
            &mut extz,
            &mut auxiliary,
            batch.error_rate,
            batch.window_length,
            batch.h_kmer_hit,
            &mut scratch,
            batch.chemical_drop,
            batch.alignment_gap_rate,
            batch.alignment_gap_max,
        );

        lines.push(format!(
            "RB name={} e_rate={:.3} wl={} khit={} drop={} graw={:.4} gmax={} fill={} on={} len={}",
            batch.name,
            batch.error_rate,
            batch.window_length,
            batch.h_kmer_hit,
            u32::from(batch.chemical_drop),
            batch.alignment_gap_rate,
            batch.alignment_gap_max,
            batch.filler_hits,
            original_length,
            kept
        ));
        for (index, region) in overlaps.iter().enumerate() {
            lines.push(format!(
                "SL i={} {} y={} x=[{},{}] y=[{},{}] strand={} ss={} nhpe={} aim={} alen={} ovl={} strong={} wli={}",
                index,
                if index < kept { "KEEP" } else { "TAIL" },
                region.y_identifier,
                region.x_position_start,
                region.x_position_end,
                region.y_position_start,
                region.y_position_end,
                region.y_strand,
                region.shared_seed,
                region.non_homopolymer_errors,
                region.is_match,
                region.alignment_length,
                region.overlap_length,
                region.strong,
                region.without_large_indel
            ));
            push_windows(lines, "  SW", &region.window_cigar);
        }
    }

    #[derive(Clone)]
    struct ManualWindow {
        x_start: i32,
        x_end: i32,
        y_start: i32,
        y_end: i32,
        error: i16,
        threshold: i16,
        cells: Vec<(u8, u32)>,
        extra_end: i16,
    }

    fn build_manual(
        x_end: u32,
        y_end: u32,
        strand: u32,
        windows: &[ManualWindow],
    ) -> OverlapRegion {
        let mut region = OverlapRegion {
            x_identifier: 0,
            x_position_start: 0,
            x_position_end: x_end,
            y_identifier: 1,
            y_position_start: 0,
            y_position_end: y_end,
            y_strand: strand,
            overlap_length: x_end + 1,
            alignment_length: x_end + 1,
            ..OverlapRegion::default()
        };
        region.fake_cigar.push(0, 0);
        region.fake_cigar.push(x_end, 0);
        for window in windows {
            let index = region.window_cigar.append_window(
                window.x_start,
                window.x_end,
                window.y_start,
                window.y_end,
                window.error,
                0,
                window.extra_end,
                window.threshold,
            );
            let start = region.window_cigar.cells.len();
            for &(operation, length) in &window.cells {
                region.window_cigar.push_cell(operation, length);
            }
            region.window_cigar.assign_cells(index, start);
        }
        region
    }

    #[allow(clippy::too_many_arguments)]
    fn run_ff_tend(
        name: &str,
        x_end: u32,
        y_end: u32,
        windows: &[ManualWindow],
        window_threshold: i64,
        distance: i64,
        distance_rate: f64,
        error_rate: f64,
        minimum_errors: i64,
        lines: &mut Vec<String>,
    ) {
        let region = build_manual(x_end, y_end, 0, windows);
        let triggered = terminal_error_trigger(
            &region,
            window_threshold,
            distance,
            distance_rate,
            error_rate,
            minimum_errors,
        );
        lines.push(format!(
            "FFT name={} ret={}",
            name,
            u32::from(triggered)
        ));
    }

    fn generate_lines() -> Vec<String> {
        let mut lines = Vec::new();
        let mut generator = Lcg(20241106);
        let query = generator.sequence(2000);
        let prefix = generator.sequence(600);
        let suffix = generator.sequence(300);

        let mut exact = prefix.clone();
        exact.extend_from_slice(&query);
        exact.extend_from_slice(&suffix);

        let mut subs = exact.clone();
        let mut position = 0;
        while position < query.len() {
            let base = &mut subs[prefix.len() + position];
            *base = if *base == b'A' { b'C' } else { b'A' };
            position += 100;
        }
        let mut ins = exact.clone();
        ins.splice(
            prefix.len() + 900..prefix.len() + 900,
            b"ACG".iter().copied(),
        );
        let mut del = exact.clone();
        del.drain(prefix.len() + 500..prefix.len() + 560);
        let mut diverged = exact.clone();
        let mut position = 0;
        while position < query.len() {
            let base = &mut diverged[prefix.len() + position];
            *base = if *base == b'A' { b'C' } else { b'A' };
            position += 8;
        }
        let mut middled = exact.clone();
        for position in 975..1550 {
            let base = &mut middled[prefix.len() + position];
            *base = if *base == b'A' { b'C' } else { b'A' };
        }
        let mut high = exact.clone();
        let mut position = 0;
        while position < query.len() {
            let base = &mut high[prefix.len() + position];
            *base = if *base == b'A' { b'C' } else { b'A' };
            position += 15;
        }
        let reversed = reverse_complement(&exact);
        let filler = generator.sequence(400);
        let mut gapblock = exact.clone();
        for position in 0..100 {
            gapblock[prefix.len() + 600 + position] = filler[position];
        }
        let mut gapblock2 = exact.clone();
        for position in 0..150 {
            gapblock2[prefix.len() + 500 + position] = filler[100 + position];
        }

        let mut homopolymer_query = query.clone();
        for base in &mut homopolymer_query[700..730] {
            *base = b'A';
        }
        let homopolymer_prefix = generator.sequence(600);
        let homopolymer_suffix = generator.sequence(300);
        let mut homopolymer_insert = homopolymer_prefix.clone();
        homopolymer_insert.extend_from_slice(&homopolymer_query);
        homopolymer_insert.extend_from_slice(&homopolymer_suffix);
        homopolymer_insert.insert(600 + 715, b'A');
        let mut homopolymer_delete = homopolymer_prefix;
        homopolymer_delete.extend_from_slice(&homopolymer_query);
        homopolymer_delete.extend_from_slice(&homopolymer_suffix);
        homopolymer_delete.remove(600 + 715);

        let mut del2 = exact.clone();
        del2.drain(prefix.len() + 1300..prefix.len() + 1330);

        let hits: Vec<(u32, u32)> = vec![(200, 800), (600, 1200), (1000, 1600), (1400, 2000)];
        let shifted_hits: Vec<(u32, u32)> =
            vec![(200, 800), (600, 1200), (1000, 1601), (1400, 2001)];
        let deleted_hits: Vec<(u32, u32)> =
            vec![(200, 800), (600, 1200), (1000, 1599), (1400, 1999)];

        let spec = |y_identifier: u32,
                    x_start: u32,
                    x_end: u32,
                    y_start: u32,
                    y_end: u32,
                    strand: u32,
                    hits: &[(u32, u32)]| OverlapSpec {
            y_identifier,
            x_start,
            x_end,
            y_start,
            y_end,
            strand,
            hits: hits.to_vec(),
        };

        let standard = |y_end: u32, strand: u32, hits: &[(u32, u32)]| {
            spec(1, 200, 1799, 800, y_end, strand, hits)
        };

        let batches = vec![
            BatchSpec {
                name: "exact",
                reads: vec![query.clone(), exact.clone()],
                overlaps: vec![standard(2399, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
            BatchSpec {
                name: "subs",
                reads: vec![query.clone(), subs.clone()],
                overlaps: vec![standard(2399, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
            BatchSpec {
                name: "ins",
                reads: vec![query.clone(), ins.clone()],
                overlaps: vec![standard(2402, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
            BatchSpec {
                name: "del",
                reads: vec![query.clone(), del.clone()],
                overlaps: vec![standard(2339, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
            BatchSpec {
                name: "del2",
                reads: vec![query.clone(), del2.clone()],
                overlaps: vec![standard(2369, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
            BatchSpec {
                name: "gapblock",
                reads: vec![query.clone(), gapblock.clone()],
                overlaps: vec![standard(2399, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
            BatchSpec {
                name: "gapblock_g0",
                reads: vec![query.clone(), gapblock.clone()],
                overlaps: vec![standard(2399, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: 0,
            },
            BatchSpec {
                name: "gapblock_g1000",
                reads: vec![query.clone(), gapblock.clone()],
                overlaps: vec![standard(2399, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: 1000,
            },
            BatchSpec {
                name: "gapblock2",
                reads: vec![query.clone(), gapblock2.clone()],
                overlaps: vec![standard(2399, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
            BatchSpec {
                name: "gapblock2_g0",
                reads: vec![query.clone(), gapblock2.clone()],
                overlaps: vec![standard(2399, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: 0,
            },
            BatchSpec {
                name: "gapblock2_g1000",
                reads: vec![query.clone(), gapblock2.clone()],
                overlaps: vec![standard(2399, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: 1000,
            },
            BatchSpec {
                name: "homopolymer_insert",
                reads: vec![homopolymer_query.clone(), homopolymer_insert.clone()],
                overlaps: vec![standard(2400, 0, &shifted_hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
            BatchSpec {
                name: "homopolymer_delete",
                reads: vec![homopolymer_query.clone(), homopolymer_delete.clone()],
                overlaps: vec![standard(2398, 0, &deleted_hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
            BatchSpec {
                name: "div",
                reads: vec![query.clone(), diverged.clone()],
                overlaps: vec![standard(2399, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
            BatchSpec {
                name: "mid",
                reads: vec![query.clone(), middled.clone()],
                overlaps: vec![standard(2399, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
            BatchSpec {
                name: "high",
                reads: vec![query.clone(), high.clone()],
                overlaps: vec![standard(2399, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
            BatchSpec {
                name: "rev",
                reads: vec![query.clone(), reversed.clone()],
                overlaps: vec![standard(2399, 1, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
            BatchSpec {
                name: "wl400",
                reads: vec![query.clone(), subs.clone()],
                overlaps: vec![standard(2399, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 400,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
            BatchSpec {
                name: "prefix",
                reads: vec![query.clone(), subs.clone()],
                overlaps: vec![standard(2399, 0, &hits)],
                filler_hits: 3,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
            BatchSpec {
                name: "ont",
                reads: vec![query.clone(), exact.clone()],
                overlaps: vec![standard(2399, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.07,
                window_length: 375,
                h_kmer_hit: 31,
                chemical_drop: true,
                alignment_gap_rate: 0.006,
                alignment_gap_max: 64,
            },
            BatchSpec {
                name: "ont_ins",
                reads: vec![query.clone(), ins.clone()],
                overlaps: vec![standard(2402, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.07,
                window_length: 375,
                h_kmer_hit: 31,
                chemical_drop: true,
                alignment_gap_rate: 0.006,
                alignment_gap_max: 64,
            },
            BatchSpec {
                name: "gapmax_exact",
                reads: vec![query.clone(), exact.clone()],
                overlaps: vec![standard(2399, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: 0,
            },
            BatchSpec {
                name: "gapmax_ins",
                reads: vec![query.clone(), ins.clone()],
                overlaps: vec![standard(2402, 0, &hits)],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: 0,
            },
            BatchSpec {
                name: "mix",
                reads: vec![
                    query.clone(),
                    exact.clone(),
                    diverged.clone(),
                    subs.clone(),
                    high.clone(),
                    reversed.clone(),
                ],
                overlaps: vec![
                    spec(1, 200, 1799, 800, 2399, 0, &hits),
                    spec(2, 200, 1799, 800, 2399, 0, &hits),
                    spec(3, 200, 1799, 800, 2399, 0, &hits),
                    spec(4, 200, 1799, 800, 2399, 0, &hits),
                    spec(5, 200, 1799, 800, 2399, 1, &hits),
                ],
                filler_hits: 0,
                error_rate: 0.04,
                window_length: 775,
                h_kmer_hit: 31,
                chemical_drop: false,
                alignment_gap_rate: -1.0,
                alignment_gap_max: -1,
            },
        ];
        for batch in &batches {
            run_batch(batch, &mut lines);
        }

        // ff_tend: a 5000 base overlap cut into four windows of 1250. The
        // window length is 384, the inspected distance shrinks to 500 (ten
        // percent of the overlap) and the rate is 0.36, so the terminal
        // regions are [0, 500) and [4500, 5000).
        let clean = vec![
            ManualWindow {
                x_start: 0,
                x_end: 1249,
                y_start: 0,
                y_end: 1249,
                error: 0,
                threshold: 0,
                cells: vec![(0, 1250)],
                extra_end: 0,
            },
            ManualWindow {
                x_start: 1250,
                x_end: 2499,
                y_start: 1250,
                y_end: 2499,
                error: 0,
                threshold: 0,
                cells: vec![(0, 1250)],
                extra_end: 0,
            },
            ManualWindow {
                x_start: 2500,
                x_end: 3749,
                y_start: 2500,
                y_end: 3749,
                error: 0,
                threshold: 0,
                cells: vec![(0, 1250)],
                extra_end: 0,
            },
            ManualWindow {
                x_start: 3750,
                x_end: 4999,
                y_start: 3750,
                y_end: 4999,
                error: 0,
                threshold: 0,
                cells: vec![(0, 1250)],
                extra_end: 0,
            },
        ];
        run_ff_tend("clean", 4999, 4999, &clean, 384, 2000, 0.1, 0.36, 128, &mut lines);

        let mut left_bad = clone_windows(&clean);
        left_bad[0].cells = vec![(0, 100), (1, 200), (0, 950)];
        run_ff_tend(
            "left_bad", 4999, 4999, &left_bad, 384, 2000, 0.1, 0.36, 128, &mut lines,
        );

        let mut left_bad_low = clone_windows(&clean);
        left_bad_low[0].cells = vec![(0, 100), (1, 200), (0, 950)];
        run_ff_tend(
            "left_bad_low", 4999, 4999, &left_bad_low, 384, 2000, 0.1, 0.9, 128, &mut lines,
        );

        let mut right_bad = clone_windows(&clean);
        right_bad[3].cells = vec![(0, 900), (1, 250), (0, 100)];
        run_ff_tend(
            "right_bad", 4999, 4999, &right_bad, 384, 2000, 0.1, 0.36, 128, &mut lines,
        );

        let mut middle_bad = clone_windows(&clean);
        middle_bad[1].cells = vec![(0, 100), (1, 1050)];
        run_ff_tend(
            "middle_bad", 4999, 4999, &middle_bad, 384, 2000, 0.1, 0.36, 128, &mut lines,
        );

        let mut unaligned_first = clone_windows(&clean);
        unaligned_first[0].error = i16::MAX;
        unaligned_first[0].cells.clear();
        unaligned_first[0].extra_end = -1;
        run_ff_tend(
            "unaligned_first",
            4999,
            4999,
            &unaligned_first,
            384,
            2000,
            0.1,
            0.36,
            128,
            &mut lines,
        );

        lines
    }

    // ------------------------------------------------- exact-recovery probe

    #[derive(Clone)]
    struct EaOverlapSpec {
        y_identifier: u32,
        x_start: u32,
        x_end: u32,
        y_start: u32,
        y_end: u32,
        strand: u32,
        hits: Vec<(u32, u32)>,
        pre_match: u8,
        windows: Vec<ManualWindow>,
    }

    struct ExactHitSpec {
        target_identifier: u32,
        reverse: u32,
        query_start: u32,
        query_end: u32,
        target_start: u32,
        target_end: u32,
        exact_region: u8,
    }

    struct EaBatchSpec {
        name: &'static str,
        direct_nec: bool,
        reads: Vec<Vec<u8>>,
        overlaps: Vec<EaOverlapSpec>,
        exact_hits: Vec<ExactHitSpec>,
        error_rate: f64,
        window_length: i64,
        h_kmer_hit: i64,
        chemical_drop: bool,
        alignment_gap_rate: f64,
        alignment_gap_max: i64,
        filler_hits: u32,
    }

    /// `setup_ea_overlap`: the overlap of the recovery probe, whose fake cigar
    /// spans the overlap's own query interval.
    fn build_ea_overlap(spec: &EaOverlapSpec, non_homopolymer_errors: u32) -> OverlapRegion {
        let mut region = OverlapRegion {
            x_identifier: 0,
            x_position_start: spec.x_start,
            x_position_end: spec.x_end,
            x_strand: 0,
            y_identifier: spec.y_identifier,
            y_position_start: spec.y_start,
            y_position_end: spec.y_end,
            y_strand: spec.strand,
            overlap_length: spec.x_end + 1 - spec.x_start,
            alignment_length: spec.x_end + 1 - spec.x_start,
            is_match: spec.pre_match,
            non_homopolymer_errors,
            ..OverlapRegion::default()
        };
        region.fake_cigar.push(spec.x_start, 0);
        region.fake_cigar.push(spec.x_end, 0);
        for window in &spec.windows {
            let index = region.window_cigar.append_window(
                window.x_start,
                window.x_end,
                window.y_start,
                window.y_end,
                window.error,
                0,
                window.extra_end,
                window.threshold,
            );
            let start = region.window_cigar.cells.len();
            for &(operation, length) in &window.cells {
                region.window_cigar.push_cell(operation, length);
            }
            region.window_cigar.assign_cells(index, start);
        }
        region
    }

    fn run_ea_batch(batch: &EaBatchSpec, lines: &mut Vec<String>) {
        let mut reads = ReadStore::default();
        for (index, sequence) in batch.reads.iter().enumerate() {
            reads.push_read(&format!("read{}", index), sequence);
        }

        let mut candidates = CandidatesList::default();
        for index in 0..batch.filler_hits {
            candidates.push(KmerHit {
                read_identifier: 4096,
                reverse_strand: false,
                offset: 5000 + index,
                self_offset: 5000 + index,
                count: 1,
            });
        }

        let mut overlaps = Vec::new();
        for spec in &batch.overlaps {
            let region = build_ea_overlap(spec, candidates.length as u32);
            for &(self_offset, offset) in &spec.hits {
                candidates.push(KmerHit {
                    read_identifier: spec.y_identifier,
                    reverse_strand: false,
                    offset,
                    self_offset,
                    count: 1,
                });
            }
            overlaps.push(region);
        }
        let original_length = overlaps.len();

        let exact_overlaps: Vec<PafHit> = batch
            .exact_hits
            .iter()
            .map(|hit| PafHit {
                query_identifier: 0,
                query_start: hit.query_start,
                query_end: hit.query_end,
                target_identifier: hit.target_identifier,
                target_start: hit.target_start,
                target_end: hit.target_end,
                matching_length: 0,
                block_length: 0,
                reverse: hit.reverse != 0,
                without_large_indel: false,
                exact_region: hit.exact_region != 0,
                deleted: false,
            })
            .collect();

        let mut extz = BitExtz::new();
        let mut auxiliary = OverlapRegion::default();
        let mut scratch = Vec::new();
        let mut exact_order = Vec::new();
        let query = batch.reads[0].clone();
        let kept = if batch.direct_nec {
            gen_hc_r_alin_nec(
                &mut overlaps,
                &mut candidates,
                &reads,
                &query,
                &mut extz,
                &mut auxiliary,
                batch.error_rate,
                batch.window_length,
                batch.h_kmer_hit,
                &mut scratch,
                batch.chemical_drop,
                batch.alignment_gap_rate,
                batch.alignment_gap_max,
            )
        } else {
            gen_hc_r_alin_ea(
                &mut overlaps,
                &mut candidates,
                &reads,
                &query,
                &mut extz,
                &mut auxiliary,
                &exact_overlaps,
                batch.error_rate,
                batch.window_length,
                batch.h_kmer_hit,
                &mut exact_order,
                &mut scratch,
                batch.chemical_drop,
                batch.alignment_gap_rate,
                batch.alignment_gap_max,
            )
        };

        lines.push(format!(
            "RB name={} mode={} e_rate={:.3} wl={} khit={} drop={} graw={:.4} gmax={} fill={} hits={} on={} len={} srtn={}",
            batch.name,
            if batch.direct_nec { "nec" } else { "ea" },
            batch.error_rate,
            batch.window_length,
            batch.h_kmer_hit,
            u32::from(batch.chemical_drop),
            batch.alignment_gap_rate,
            batch.alignment_gap_max,
            batch.filler_hits,
            batch.exact_hits.len(),
            original_length,
            kept,
            exact_order.len()
        ));
        for (index, key) in exact_order.iter().enumerate() {
            lines.push(format!(" SRT i={} key={:016x}", index, key));
        }
        for (index, region) in overlaps.iter().enumerate() {
            lines.push(format!(
                "SL i={} {} y={} x=[{},{}] y=[{},{}] strand={} ss={} nhpe={} aim={} alen={} ovl={} strong={} wli={}",
                index,
                if index < kept { "KEEP" } else { "TAIL" },
                region.y_identifier,
                region.x_position_start,
                region.x_position_end,
                region.y_position_start,
                region.y_position_end,
                region.y_strand,
                region.shared_seed,
                region.non_homopolymer_errors,
                region.is_match,
                region.alignment_length,
                region.overlap_length,
                region.strong,
                region.without_large_indel
            ));
            push_windows(lines, "  SW", &region.window_cigar);
        }
    }

    fn witness_window() -> ManualWindow {
        ManualWindow {
            x_start: 0,
            x_end: 99,
            y_start: 600,
            y_end: 699,
            error: 0,
            threshold: 0,
            cells: vec![(0, 100)],
            extra_end: 0,
        }
    }

    fn ea_standard(y_identifier: u32, y_end: u32, strand: u32) -> EaOverlapSpec {
        EaOverlapSpec {
            y_identifier,
            x_start: 0,
            x_end: 1999,
            y_start: 600,
            y_end,
            strand,
            hits: vec![(200, 800), (600, 1200), (1000, 1600), (1400, 2000)],
            pre_match: 0,
            windows: Vec::new(),
        }
    }

    fn with_witness(mut spec: EaOverlapSpec) -> EaOverlapSpec {
        spec.windows.push(witness_window());
        spec
    }

    fn exact_hit(
        target_identifier: u32,
        reverse: u32,
        target_start: u32,
        target_end: u32,
        exact_region: u8,
    ) -> ExactHitSpec {
        ExactHitSpec {
            target_identifier,
            reverse,
            query_start: 0,
            query_end: 2000,
            target_start,
            target_end,
            exact_region,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn ea_batch(
        name: &'static str,
        direct_nec: bool,
        reads: Vec<Vec<u8>>,
        overlaps: Vec<EaOverlapSpec>,
        exact_hits: Vec<ExactHitSpec>,
        error_rate: f64,
        window_length: i64,
        h_kmer_hit: i64,
        chemical_drop: bool,
        alignment_gap_rate: f64,
        alignment_gap_max: i64,
        filler_hits: u32,
    ) -> EaBatchSpec {
        EaBatchSpec {
            name,
            direct_nec,
            reads,
            overlaps,
            exact_hits,
            error_rate,
            window_length,
            h_kmer_hit,
            chemical_drop,
            alignment_gap_rate,
            alignment_gap_max,
            filler_hits,
        }
    }

    fn generate_ea_lines() -> Vec<String> {
        let mut lines = Vec::new();
        let mut generator = Lcg(20241106);
        let query = generator.sequence(2000);
        let prefix = generator.sequence(600);
        let suffix = generator.sequence(300);

        let mut exact = prefix.clone();
        exact.extend_from_slice(&query);
        exact.extend_from_slice(&suffix);

        let mut subs = exact.clone();
        let mut position = 0;
        while position < query.len() {
            let base = &mut subs[prefix.len() + position];
            *base = if *base == b'A' { b'C' } else { b'A' };
            position += 100;
        }
        let mut high = exact.clone();
        let mut position = 0;
        while position < query.len() {
            let base = &mut high[prefix.len() + position];
            *base = if *base == b'A' { b'C' } else { b'A' };
            position += 15;
        }
        let reversed = reverse_complement(&exact);

        let mut batches = Vec::new();

        // A filler group shifts the chain index away from zero: the exact
        // branch must park that index in shared_seed and zero the error
        // field.
        batches.push(ea_batch(
            "ea_exact_prefix",
            false,
            vec![query.clone(), exact.clone()],
            vec![with_witness(ea_standard(1, 2599, 0))],
            vec![exact_hit(1, 0, 600, 2600, 1)],
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            3,
        ));
        // The same for an overlap the shortcut rejects.
        batches.push(ea_batch(
            "ea_fallback_prefix",
            false,
            vec![query.clone(), subs.clone()],
            vec![with_witness(ea_standard(1, 2599, 0))],
            vec![exact_hit(1, 0, 600, 2600, 1)],
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            3,
        ));
        // An empty chain list returns before anything else, even with stored
        // hits.
        batches.push(ea_batch(
            "ea_empty",
            false,
            vec![query.clone(), exact.clone()],
            Vec::new(),
            vec![exact_hit(1, 0, 600, 2600, 1)],
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            0,
        ));
        // No exact hit at all: the plain refinement runs.
        batches.push(ea_batch(
            "ea_no_el",
            false,
            vec![query.clone(), exact.clone()],
            vec![ea_standard(1, 2599, 0)],
            vec![exact_hit(1, 0, 600, 2600, 0)],
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            0,
        ));
        // One stored exact region that covers the whole overlap and matches.
        batches.push(ea_batch(
            "ea_exact",
            false,
            vec![query.clone(), exact.clone()],
            vec![with_witness(ea_standard(1, 2599, 0))],
            vec![exact_hit(1, 0, 600, 2600, 1)],
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            0,
        ));
        // Same coordinates, one base differs inside the exact region.
        batches.push(ea_batch(
            "ea_exact_wrong_seq",
            false,
            vec![query.clone(), subs.clone()],
            vec![with_witness(ea_standard(1, 2599, 0))],
            vec![exact_hit(1, 0, 600, 2600, 1)],
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            0,
        ));
        // Stored coordinates do not match the new overlap.
        batches.push(ea_batch(
            "ea_coord_mismatch",
            false,
            vec![query.clone(), exact.clone()],
            vec![with_witness(ea_standard(1, 2599, 0))],
            vec![exact_hit(1, 0, 601, 2601, 1)],
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            0,
        ));
        // A reverse-strand exact region; the stored y coordinates are in the
        // RC frame.
        batches.push(ea_batch(
            "ea_rev_exact",
            false,
            vec![query.clone(), reversed.clone()],
            vec![with_witness(ea_standard(1, 2599, 1))],
            vec![exact_hit(1, 1, 600, 2600, 1)],
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            0,
        ));
        // The only stored hit has the other strand: the lookup key misses.
        batches.push(ea_batch(
            "ea_strand_mismatch",
            false,
            vec![query.clone(), exact.clone()],
            vec![with_witness(ea_standard(1, 2599, 0))],
            vec![exact_hit(1, 1, 600, 2600, 1)],
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            0,
        ));
        // Two stored hits for the same target: only the first sorted one is
        // ever consulted.
        batches.push(ea_batch(
            "ea_two_hits_first_mismatch",
            false,
            vec![query.clone(), exact.clone()],
            vec![with_witness(ea_standard(1, 2599, 0))],
            vec![exact_hit(1, 0, 601, 2601, 1), exact_hit(1, 0, 600, 2600, 1)],
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            0,
        ));
        batches.push(ea_batch(
            "ea_two_hits_first_match",
            false,
            vec![query.clone(), exact.clone()],
            vec![with_witness(ea_standard(1, 2599, 0))],
            vec![exact_hit(1, 0, 600, 2600, 1), exact_hit(1, 0, 601, 2601, 1)],
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            0,
        ));
        // Two targets on different strands: the sort orders them by
        // (target, strand) before the merge.
        batches.push(ea_batch(
            "ea_two_targets",
            false,
            vec![query.clone(), exact.clone(), reversed.clone()],
            vec![
                with_witness(ea_standard(2, 2599, 1)),
                with_witness(ea_standard(1, 2599, 0)),
            ],
            vec![exact_hit(2, 1, 600, 2600, 1), exact_hit(1, 0, 600, 2600, 1)],
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            0,
        ));
        // Mixed list: exact recovery, a dropped overlap, a surviving one.
        batches.push(ea_batch(
            "ea_mixed",
            false,
            vec![query.clone(), exact.clone(), high.clone(), subs.clone()],
            vec![
                with_witness(ea_standard(1, 2599, 0)),
                ea_standard(2, 2599, 0),
                ea_standard(3, 2599, 0),
            ],
            vec![exact_hit(1, 0, 600, 2600, 1)],
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            0,
        ));
        // Fallback with the ONT filters enabled.
        batches.push(ea_batch(
            "ea_nec_filters",
            false,
            vec![query.clone(), high.clone()],
            vec![ea_standard(1, 2599, 0)],
            vec![exact_hit(1, 0, 600, 2600, 0)],
            0.07,
            375,
            31,
            true,
            0.006,
            64,
            0,
        ));
        // A stored exact hit whose window list already holds an aligned
        // window.
        let mut merge = ea_standard(1, 2599, 0);
        merge.windows.push(witness_window());
        batches.push(ea_batch(
            "ea_exact_merge",
            false,
            vec![query.clone(), exact.clone()],
            vec![merge.clone()],
            vec![exact_hit(1, 0, 600, 2600, 1)],
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            0,
        ));
        merge.windows[0].x_end = 98;
        batches.push(ea_batch(
            "ea_exact_no_merge",
            false,
            vec![query.clone(), exact.clone()],
            vec![merge.clone()],
            vec![exact_hit(1, 0, 600, 2600, 1)],
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            0,
        ));
        // Direct gen_hc_r_alin_nec calls: a pre-matched entry passes through
        // without being aligned, keeping its window list and taking part in
        // the compaction.
        let mut premixed = ea_standard(4, 2599, 0);
        premixed.pre_match = 1;
        premixed.windows.push(ManualWindow {
            x_start: 0,
            x_end: 1999,
            y_start: 600,
            y_end: 2599,
            error: 0,
            threshold: 0,
            cells: vec![(0, 2000)],
            extra_end: 0,
        });
        batches.push(ea_batch(
            "nec_mixed",
            true,
            vec![query.clone(), subs.clone(), high.clone(), exact.clone()],
            vec![
                ea_standard(1, 2599, 0),
                ea_standard(2, 2599, 0),
                premixed.clone(),
            ],
            Vec::new(),
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            0,
        ));
        let mut prem = ea_standard(1, 2599, 0);
        prem.pre_match = 1;
        batches.push(ea_batch(
            "nec_all_prematched",
            true,
            vec![query.clone(), exact.clone()],
            vec![prem],
            Vec::new(),
            0.04,
            775,
            31,
            false,
            -1.0,
            -1,
            0,
        ));

        for batch in &batches {
            run_ea_batch(batch, &mut lines);
        }
        lines
    }

    #[test]
    #[ignore = "requires reference/ec_hc_ralin_ea_probe.out from the C++ probe"]
    fn exact_overlap_recovery_matches_the_reference() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("reference")
            .join("ec_hc_ralin_ea_probe.out");
        let contents = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("cannot read {}: {}", path.display(), error));
        let golden: Vec<&str> = contents.lines().collect();
        let generated = generate_ea_lines();

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
                index + 1,
                expected,
                actual
            );
        }
    }

    fn clone_windows(windows: &[ManualWindow]) -> Vec<ManualWindow> {
        windows
            .iter()
            .map(|window| ManualWindow {
                x_start: window.x_start,
                x_end: window.x_end,
                y_start: window.y_start,
                y_end: window.y_end,
                error: window.error,
                threshold: window.threshold,
                cells: window.cells.clone(),
                extra_end: window.extra_end,
            })
            .collect()
    }

    #[test]
    #[ignore = "requires reference/ec_hc_ralin_probe.out from the C++ probe"]
    fn whole_overlap_refinement_matches_the_reference() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("reference")
            .join("ec_hc_ralin_probe.out");
        let contents = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("cannot read {}: {}", path.display(), error));
        let golden: Vec<&str> = contents.lines().collect();
        let generated = generate_lines();

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
                index + 1,
                expected,
                actual
            );
        }
    }
}

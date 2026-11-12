//! `gen_shared_intervals` (`inter.cpp` 7169) and its consensus core
//! `gen_shared_trace` (`inter.cpp` 7106), the shared-interval pass of the UL
//! re-alignment pipeline. The sweep line collects the query intervals covered
//! by at least two overlaps; each such region is then reduced to its shared
//! pieces.
//!
//! `convert_ul_ov_t` (`Correct.cpp` 19160), `govlp_check` (`Correct.cpp`
//! 19216) and `check_connect_ug` / `check_connect_rg` (`Correct.cpp` 19174 /
//! 19190) support the containment check inside the trace.

use crate::assembly::asg::{Asg, MaUtg};
use crate::assembly::hit::ma_hit_to_arc;
use crate::assembly::ul::{infer_rovlp_unitigs, UlOv};
use crate::assembly::unitig::UnitigGraph;
use crate::container::ksort::{radix_sort_by_key_32, radix_sort_u64};
use crate::correction::hc_align::{
    align_ul_ed_post_extz, gen_extend_err_exz, window_se_by_normalized_start, BitExtz,
    HcAlignParameters,
};
use crate::correction::interval::{get_interval, IntervalBuffer, WindowPool};
use crate::correction::recalculate::similarity_pass;
use crate::correction::region::THRESHOLD_MAXIMUM;
use crate::correction::target::TargetSource;
use crate::correction::ug_lalign::{rechain_alignment_ul, ChainBracket};
use crate::correction::ul_align::{
    fusion_chain_ovlp, is_alnw, ovlp_base_aln, ALNW_SENTINEL,
};
use crate::correction::window_chain::is_unaligned_window;
use crate::correction::window::number_of_windows;
use crate::overlap::paf::PafHit;
use crate::overlap::region::{CandidatesList, ChainData, KmerHit, OverlapRegion};
use crate::reads::read_store::ReadStore;

/// `G_CHAIN_BW` (`inter.h` 7).
pub const G_CHAIN_BW: i64 = 16;
/// `G_CHAIN_INDEL` (`inter.h` 27).
pub const G_CHAIN_INDEL: u32 = 128;
/// `N_GCHAIN_RATE` (`inter.h` 29).
pub const N_GCHAIN_RATE: f64 = 0.04;

/// The `ug_opt_t` slice the shared-interval pass reads.
pub struct UlRescallContext<'a> {
    pub unitig_graph: &'a UnitigGraph,
    pub sources: &'a [Vec<PafHit>],
    /// The `ma_hit_t::cc` counters parallel to `sources`.
    pub cc: &'a [Vec<u32>],
    pub reads: &'a ReadStore,
    pub maximum_hang: i32,
    pub maximum_hang_rate: f32,
    pub minimum_overlap: i32,
}

/// `convert_ul_ov_t` (`Correct.cpp` 19160): an overlap region on the unitig's
/// own orientation.
fn convert_ul_ov_t(region: &OverlapRegion, unitigs: &[MaUtg]) -> UlOv {
    let tn = region.y_identifier;
    let rev = region.y_strand != 0;
    let (ts, te) = if rev {
        let length = unitigs[tn as usize].length;
        (
            length.wrapping_sub(region.y_position_end + 1),
            length.wrapping_sub(region.y_position_start),
        )
    } else {
        (region.y_position_start, region.y_position_end + 1)
    };
    UlOv {
        qn: u32::MAX,
        qs: region.x_position_start,
        qe: region.x_position_end + 1,
        tn,
        ts,
        te,
        sec: region.non_homopolymer_errors,
        el: true,
        rev,
    }
}

/// `diff <= mm` with `mm = max(dq, dt) * diff_ec_ul`, floored at `bw`.
fn within_band(dq: i64, dt: i64, bw: i64, diff_ec_ul: f64) -> bool {
    let difference = (dq - dt).abs();
    let mut allowed = (dq.max(dt) as f64 * diff_ec_ul) as i64;
    if allowed < bw {
        allowed = bw;
    }
    difference <= allowed
}

/// `check_connect_ug` (`Correct.cpp` 19174).
pub(crate) fn check_connect_ug(graph: &Asg, v: u32, w: u32, bw: i64, diff_ec_ul: f64, dq: i64) -> bool {
    for arc in graph.arcs_of(v) {
        if arc.deleted || arc.v != w {
            continue;
        }
        return within_band(dq, i64::from(arc.overlap), bw, diff_ec_ul);
    }
    false
}

/// `check_connect_rg` (`Correct.cpp` 19190).
pub(crate) fn check_connect_rg(
    context: &UlRescallContext<'_>,
    uv: u32,
    uw: u32,
    bw: i64,
    diff_ec_ul: f64,
    dq: i64,
) -> bool {
    let uv_unitig = &context.unitig_graph.unitigs[(uv >> 1) as usize];
    let uw_unitig = &context.unitig_graph.unitigs[(uw >> 1) as usize];
    if uv_unitig.circular || uw_unitig.circular {
        return false;
    }
    let end_vertex = |unitig: &MaUtg, strand: u32| -> u32 {
        let index = if strand == 1 { 0 } else { unitig.reads.len() - 1 };
        ((unitig.reads[index] >> 32) as u32) ^ strand
    };
    let rv = end_vertex(uv_unitig, uv & 1);
    let rw = end_vertex(uw_unitig, uw & 1);
    let query = (rv >> 1) as usize;
    for hit in &context.sources[query] {
        let qn = hit.query_identifier;
        let tn = hit.target_identifier;
        if tn != rw >> 1 {
            continue;
        }
        let (result, arc) = ma_hit_to_arc(
            hit,
            context.reads.read_length(qn as usize) as i32,
            context.reads.read_length(tn as usize) as i32,
            context.maximum_hang,
            context.maximum_hang_rate,
            context.minimum_overlap,
        );
        if result < 0 {
            continue;
        }
        if (arc.ul >> 32) as u32 != rv || arc.v != rw {
            continue;
        }
        return within_band(dq, i64::from(arc.overlap), bw, diff_ec_ul);
    }
    false
}

/// `govlp_check` (`Correct.cpp` 19216).
fn govlp_check(
    context: &UlRescallContext<'_>,
    bw: i64,
    diff_ec_ul: f64,
    li: &UlOv,
    lj: &UlOv,
) -> bool {
    let query_overlap = infer_rovlp_unitigs(li, lj, &context.unitig_graph.unitigs);
    let li_v = (((li.tn) << 1) | u32::from(li.rev)) ^ 1;
    let lj_v = (((lj.tn) << 1) | u32::from(lj.rev)) ^ 1;
    check_connect_ug(&context.unitig_graph.graph, li_v, lj_v, bw, diff_ec_ul, query_overlap)
        || check_connect_rg(context, li_v, lj_v, bw, diff_ec_ul, query_overlap)
}

/// `gen_shared_trace` (`inter.cpp` 7106): reduces the overlaps active in
/// `[s, e]` to their shared query pieces, appending them to `res`. Returns the
/// number of `id_a` entries that reach past `e` (the survivors, compacted to
/// the front).
#[allow(clippy::too_many_arguments)]
fn gen_shared_trace(
    regions: &mut [OverlapRegion],
    id_a: &mut Vec<u64>,
    s: u64,
    e: i64,
    context: &UlRescallContext<'_>,
    bw: i64,
    diff_ec_ul: f64,
    is_sorted: &mut bool,
    res: &mut Vec<UlOv>,
) -> usize {
    id_a.retain(|&id| (regions[id as usize].x_position_end as u64) >= s);
    let id_n = id_a.len();
    if id_n == 0 {
        return 0;
    }
    let mut remove_count = 0usize;
    for k in 0..id_n {
        let p_index = id_a[k];
        let mut p = convert_ul_ov_t(&regions[p_index as usize], &context.unitig_graph.unitigs);
        p.qn = p_index as u32;
        if (regions[p_index as usize].x_position_end as i64) <= e {
            remove_count += 1;
        }
        let p_end = regions[p_index as usize].x_position_end;
        let mut i = 0usize;
        while i < id_n
            && (regions[id_a[i] as usize].x_position_start as u32) <= p_end
        {
            if i != k {
                let q_index = id_a[i];
                let mut q = convert_ul_ov_t(&regions[q_index as usize], &context.unitig_graph.unitigs);
                q.qn = q_index as u32;
                let (li, lj) = if p.qe > q.qe || (p.qe == q.qe && p.qs >= q.qs) {
                    (&p, &q)
                } else {
                    (&q, &p)
                };
                let overlap_start = li.qs.max(lj.qs);
                let overlap_end = li.qe.min(lj.qe);
                let overlap = if overlap_end > overlap_start {
                    overlap_end - overlap_start
                } else {
                    0
                };
                let same_as_last = res.last().map(|t| t.qn == p.qn).unwrap_or(false);
                if overlap != 0 {
                    if same_as_last {
                        let last = res.last().unwrap();
                        if last.qs <= overlap_start && last.qe >= overlap_end {
                            i += 1;
                            continue;
                        }
                    }
                    let mut skip = false;
                    if lj.qs <= li.qs + G_CHAIN_INDEL {
                        skip = govlp_check(context, bw, diff_ec_ul, li, lj);
                    } else if lj.qe + G_CHAIN_INDEL >= li.qe
                        && lj.qs + G_CHAIN_INDEL >= li.qs
                    {
                        skip = govlp_check(context, bw, diff_ec_ul, lj, li);
                    }
                    if !skip {
                        let extend = if same_as_last {
                            let last = res.last().unwrap();
                            last.qs <= overlap_start && last.qe >= overlap_start
                        } else {
                            false
                        };
                        if extend {
                            let last = res.last_mut().unwrap();
                            if overlap_end > last.qe {
                                last.qe = overlap_end;
                            }
                        } else {
                            res.push(UlOv {
                                qn: p.qn,
                                qs: overlap_start,
                                qe: overlap_end,
                                ..UlOv::default()
                            });
                        }
                    }
                }
            }
            i += 1;
        }
        if let Some(last) = res.last() {
            if last.qn == p.qn {
                let region = &mut regions[p_index as usize];
                region.alignment_length += 1;
                if region.alignment_length > 1 {
                    *is_sorted = false;
                }
            }
        }
    }

    if remove_count != 0 {
        id_a.retain(|&id| (regions[id as usize].x_position_end as i64) > e);
    }
    id_a.len()
}

/// `gen_shared_intervals` (`inter.cpp` 7169): the sweep line over every
/// overlap's query interval, calling [`gen_shared_trace`] for each region
/// covered by at least two overlaps, then sorting the result by `(qn, qs)`.
///
/// `res` is the caller's `aln` buffer; it is cleared and repopulated. Returns
/// the number of shared intervals.
pub fn gen_shared_intervals(
    regions: &mut [OverlapRegion],
    context: &UlRescallContext<'_>,
    _window_length: i64,
    res: &mut Vec<UlOv>,
) -> usize {
    res.clear();
    let overlap_count = regions.len();
    let endpoint_count = overlap_count << 1;
    let mut endpoints: Vec<u64> = Vec::with_capacity(endpoint_count);
    for index in 0..overlap_count {
        endpoints.push((u64::from(regions[index].x_position_start) << 1) << 32 | index as u64);
        endpoints.push((((u64::from(regions[index].x_position_end) << 1) + 1) << 32) | index as u64);
        regions[index].alignment_length = 0;
    }
    radix_sort_u64(&mut endpoints[..endpoint_count]);

    let mut active: Vec<u64> = Vec::new();
    let mut depth: i64 = 0;
    let mut begin: u64 = 0;
    let mut is_sorted = true;
    for position in 0..endpoint_count {
        let previous_depth = depth;
        let entry = endpoints[position];
        if ((entry >> 32) & 1) != 0 {
            depth -= 1;
        } else {
            depth += 1;
            active.push(entry & 0xffff_ffff);
        }
        if previous_depth < 2 && depth >= 2 {
            begin = entry >> 33;
        } else if previous_depth >= 2 && depth < 2 {
            let end = (entry >> 33) as i64;
            let survivors = gen_shared_trace(
                regions,
                &mut active,
                begin,
                end,
                context,
                G_CHAIN_BW,
                N_GCHAIN_RATE,
                &mut is_sorted,
                res,
            );
            active.truncate(survivors);
        }
    }

    if res.is_empty() {
        return 0;
    }
    let result_count = res.len();
    radix_sort_by_key_32(res, &|entry: &UlOv| entry.qn);
    let mut run_start = 0usize;
    for index in 1..=result_count {
        if index == result_count || res[index].qn != res[run_start].qn {
            if index - run_start > 1 {
                radix_sort_by_key_32(&mut res[run_start..index], &|entry: &UlOv| entry.qs);
            }
            run_start = index;
        }
    }
    let _ = is_sorted;
    res.len()
}

/// `get_adjust_winid` (`Correct.cpp` 3079): the window index in `region` whose
/// `x_start` equals `win_beg`, or `-1`.
fn get_adjust_winid(region: &OverlapRegion, win_beg: i64, win_len: i64) -> i64 {
    let x_start = region.x_position_start as i64;
    let mut win_id = (win_beg - (x_start / win_len) * win_len) / win_len;
    let count = region.window_cigar.windows.len();
    if win_id >= 0
        && (win_id as usize) < count
        && region.window_cigar.windows[win_id as usize].x_start as i64 == win_beg
    {
        return win_id;
    }
    if count == 0 {
        return -1;
    }
    if win_id as u64 > count as u64 {
        win_id = count as i64;
    }
    let mut index = win_id - 1;
    while index >= 0 {
        let window = &region.window_cigar.windows[index as usize];
        if window.x_start as i64 == win_beg {
            return index;
        }
        if (window.x_start as i64) < win_beg {
            return -1;
        }
        index -= 1;
    }
    -1
}

/// `set_herror_win` (`Correct.cpp` 3099): marks overlaps whose window has no
/// supporting alignment as `is_match = 4` and records the freshly matched
/// windows in `v_idx` (`(offset << 32) | count`, followed by the
/// `(window_start << 32) | window_end` runs).
#[allow(clippy::too_many_arguments)]
pub fn set_herror_win(
    overlaps: &mut [OverlapRegion],
    dummy: &mut IntervalBuffer,
    v_idx: &mut Vec<u64>,
    maximum_overlap_divergence: f64,
    read_length: i64,
    block_length: i64,
) {
    let mut pool = WindowPool::new(
        read_length,
        block_length,
        (1.0 / maximum_overlap_divergence) as i64,
    );
    let window_length = pool.window_length;
    let index_count = number_of_windows(0, read_length, block_length) as usize;
    v_idx.clear();
    v_idx.resize(index_count, 0);
    let mut index_i = 0usize;
    let mut flag = 0i32;
    while let Some((window_start, window_end)) = pool.next_window() {
        if flag == -2 {
            break;
        }
        dummy.length = 0;
        dummy.partial_length = 0;
        flag = get_interval(window_start, window_end, overlaps, dummy, window_length);
        v_idx[index_i] = (v_idx.len() as u64) << 32;
        index_i += 1;

        // Does a matched window of the full-window group already cover it?
        let mut has_matched = false;
        for i in 0..dummy.length as usize {
            let identifier = dummy.overlap_identifiers[i] as u32 as usize;
            let is_match = overlaps[identifier].is_match;
            if is_match != 3 && is_match != 4 {
                continue;
            }
            if get_adjust_winid(&overlaps[identifier], window_start, window_length) >= 0 {
                has_matched = true;
                break;
            }
        }
        if has_matched {
            continue;
        }

        // Every full window is unmatched: mark the overlaps it covers.
        let mut matched = 0usize;
        for i in 0..dummy.length as usize {
            let identifier = dummy.overlap_identifiers[i] as u32 as usize;
            let is_match = overlaps[identifier].is_match;
            if is_match != 3 && is_match != 4 {
                continue;
            }
            overlaps[identifier].is_match = 4;
            overlaps[identifier].alignment_length += (window_end + 1 - window_start) as u32;
            matched += 1;
        }
        if matched > 0 {
            v_idx.push(((window_start as u64) << 32) | window_end as u64);
            v_idx[index_i - 1] += 1;
        }

        // The partial (short) windows at the end of the buffer.
        let size = dummy.overlap_identifiers.len();
        let partial = dummy.partial_length as usize;
        let start_short = size - partial;
        let mut mlen = start_short;
        let mut unmatched = 0usize;
        for i in start_short..size {
            let identifier = dummy.overlap_identifiers[i] as u32 as usize;
            let is_match = overlaps[identifier].is_match;
            if is_match != 3 && is_match != 4 {
                continue;
            }
            let (window_se, window_ee) = window_se_by_normalized_start(
                &overlaps[identifier],
                window_start,
                block_length,
            );
            let window_id = get_adjust_winid(&overlaps[identifier], window_se, block_length);
            if window_id >= 0 {
                dummy.overlap_identifiers[i] =
                    ((window_id as u64) << 32) | (dummy.overlap_identifiers[i] as u32 as u64);
                if mlen != i {
                    dummy.overlap_identifiers.swap(mlen, i);
                }
                mlen += 1;
            } else {
                dummy.overlap_identifiers[i] =
                    (0xffff_ffffu64 << 32) | (dummy.overlap_identifiers[i] as u32 as u64);
                unmatched += 1;
                let _ = window_ee;
            }
        }
        if unmatched == 0 {
            continue;
        }
        for i in mlen..size {
            let identifier = dummy.overlap_identifiers[i] as u32 as usize;
            let is_match = overlaps[identifier].is_match;
            if is_match != 3 && is_match != 4 {
                continue;
            }
            let (window_se, window_ee) =
                window_se_by_normalized_start(&overlaps[identifier], window_start, block_length);
            let mut covered = false;
            for k in start_short..mlen {
                let matched_id = dummy.overlap_identifiers[k] as u32 as usize;
                let window_id = (dummy.overlap_identifiers[k] >> 32) as usize;
                let window = &overlaps[matched_id].window_cigar.windows[window_id];
                if window.x_start as i64 <= window_se && window.x_end as i64 >= window_ee {
                    covered = true;
                    break;
                }
            }
            if !covered {
                let region = &mut overlaps[identifier];
                region.is_match = 4;
                region.alignment_length += (window_ee + 1 - window_se) as u32;
                v_idx.push(((window_se as u64) << 32) | window_ee as u64);
                v_idx[index_i - 1] += 1;
            }
        }
    }
}

/// `ul_lalign` first pass (`Correct.cpp` 22443, the `aln == NULL` branch):
/// base-aligns every overlap against the unitig index, marks the ones whose
/// windowed alignment passes the similarity check, charges the coverage gaps
/// via `set_herror_win`/`gen_extend_err_exz`, then compacts the survivors to
/// the front (`is_match = 1`, `non_homopolymer_errors` = extension error).
///
/// Returns the number of surviving overlaps; `regions` is compacted in place.
#[allow(clippy::too_many_arguments)]
pub fn ul_lalign_base(
    regions: &mut [OverlapRegion],
    targets: &TargetSource<'_>,
    query: &[u8],
    extz: &mut BitExtz,
    dummy: &mut IntervalBuffer,
    v_idx: &mut Vec<u64>,
    error_rate: f64,
    window_length: i64,
    query_length: i64,
) -> usize {
    // `clear_Correct_dumy(dumy, ol, km)`: the identifier buffer holds one slot
    // per overlap and the window scan restarts from the front.
    if dummy.overlap_identifiers.len() < regions.len() {
        dummy.overlap_identifiers.resize(regions.len(), 0);
    }
    dummy.length = 0;
    dummy.partial_length = 0;
    dummy.start_index = 0;

    let overlap_count = regions.len();
    let parameters = HcAlignParameters {
        error_rate,
        window_length,
        overlap_cutoff: -1.0,
        force_alignment: false,
    };

    let mut mapped_overlaps_length = 0i64;
    for index in 0..overlap_count {
        let overlap_length =
            regions[index].x_position_end as i64 + 1 - regions[index].x_position_start as i64;
        regions[index].shared_seed = regions[index].non_homopolymer_errors as i32;
        if !align_ul_ed_post_extz(&mut regions[index], extz, targets, query, &parameters) {
            continue;
        }
        if similarity_pass(
            overlap_length,
            regions[index].alignment_length as i64,
            true,
            -1.0,
            None,
        ) {
            regions[index].is_match = 3;
            mapped_overlaps_length += regions[index].alignment_length as i64;
        }
    }

    if mapped_overlaps_length > 0 {
        set_herror_win(
            regions,
            dummy,
            v_idx,
            error_rate,
            query_length,
            window_length,
        );
    }

    let maximum_error_rate = error_rate * 1.5;
    let mut kept = 0usize;
    for index in 0..overlap_count {
        let extension = gen_extend_err_exz(
            &mut regions[index],
            targets,
            query,
            extz,
            window_length,
            -1.0,
            error_rate,
            maximum_error_rate + 0.000001,
            THRESHOLD_MAXIMUM,
            true,
            Some(v_idx),
        );
        let error_rate_result = extension.error_rate;
        let error_bases = extension.error_bases;
        regions[index].is_match = 0;
        if error_rate_result <= error_rate {
            if kept != index {
                regions.swap(kept, index);
            }
            regions[kept].is_match = 1;
            regions[kept].non_homopolymer_errors = error_bases as u32;
            kept += 1;
        }
    }
    kept
}

/// `get_win_aln` (`Correct.cpp` 14917): the target coordinates and error of an
/// aligned window, only when it is adjacent (in target space) to the previous
/// window.
pub(crate) fn get_win_aln(region: &OverlapRegion, wid: usize) -> Option<(i64, i64, i64)> {
    let windows = &region.window_cigar.windows;
    if wid > 0
        && windows[wid].y_end != -1
        && windows[wid - 1].y_end != -1
        && windows[wid].y_end > windows[wid - 1].y_end
    {
        return Some((
            windows[wid - 1].y_end as i64 + 1,
            windows[wid].y_end as i64 + 1,
            windows[wid].error as i64,
        ));
    }
    None
}

/// `ovlp_win_check` (`Correct.cpp` 22190): whether two windows of the same
/// overlap are collinear enough to belong to one consensus chain.
pub(crate) fn ovlp_win_check(
    region: &OverlapRegion,
    id0: u32,
    id1: u32,
    max_lgap: i64,
    small_bw_rate: f64,
    min_small_bw: i64,
) -> bool {
    if id0 == u32::MAX || id1 == u32::MAX {
        return true;
    }
    let Some((ts0, te0, _)) = get_win_aln(region, id0 as usize) else {
        return true;
    };
    let windows = &region.window_cigar.windows;
    let qs0 = windows[id0 as usize].x_start as i64;
    let qe0 = windows[id0 as usize].x_end as i64 + 1;
    let Some((ts1, te1, _)) = get_win_aln(region, id1 as usize) else {
        return true;
    };
    let qs1 = windows[id1 as usize].x_start as i64;
    let qe1 = windows[id1 as usize].x_end as i64 + 1;

    if qs1 < qs0 || qe1 < qe0 {
        return false;
    }
    if ts1 < ts0 || te1 < te0 {
        return false;
    }
    let dq = qe1 - qs0;
    let dr = te1 - ts0;
    let dd = if dq >= dr { dq - dr } else { dr - dq };
    if ts1 < te0 && qe0 == qs1 {
        let dm = dq.min(dr) as f64;
        if (dd as f64) > dm * small_bw_rate && dd > min_small_bw {
            return false;
        }
    } else if dd > max_lgap {
        return false;
    }
    true
}

/// `dp_commen_sketch` (`Correct.cpp` 22211): longest chain of the shared query
/// windows `win_a` under the `ovlp_win_check` collinearity, then marks the
/// chain windows' `extra_end`.
#[allow(clippy::too_many_arguments)]
pub fn dp_commen_sketch(
    aln: &[UlOv],
    regions: &mut [OverlapRegion],
    id_a: &[u64],
    id_n: i64,
    win_a: &mut [u64],
    mut win_n: i64,
    dp: &mut [u64],
    n_skip: i64,
    wl: i64,
    cov: i64,
    max_lgap: i64,
    sgap_rate: f64,
    sgap: i64,
) -> u64 {
    if win_n == 0 {
        return 0;
    }
    let high_bit = 0x8000_0000_0000_0000u64;
    if n_skip < win_n {
        let mut long_sc: i64 = -1;
        let mut long_idx: i64 = -1;
        let mut k = 0i64;
        while k < n_skip {
            let parent = if k > 0 { (k - 1) as u64 } else { 0xffff_ffffu64 };
            dp[k as usize] = (parent << 32) | (k + 1) as u64;
            let max_sc = k + 1;
            if long_sc < max_sc {
                long_sc = max_sc;
                long_idx = k;
            }
            k += 1;
        }
        while k < win_n {
            let ws = win_a[k as usize];
            let sc: i64 = 1;
            let mut p: i64 = -1;
            let mut max_sc: i64 = sc;
            let mut m = k - 1;
            while m >= 0 {
                let ws0 = win_a[m as usize];
                let mut i = 0i64;
                while i < id_n {
                    let region = &regions[aln[id_a[i as usize] as usize].qn as usize];
                    let wid = crate::correction::window::window_id_by_start(
                        region.x_position_start as i32,
                        region.x_position_end as i32,
                        ws as i64,
                        wl,
                    )
                    .0;
                    let wid0 = crate::correction::window::window_id_by_start(
                        region.x_position_start as i32,
                        region.x_position_end as i32,
                        ws0 as i64,
                        wl,
                    )
                    .0;
                    if !ovlp_win_check(
                        region,
                        wid0 as u32,
                        wid as u32,
                        max_lgap,
                        sgap_rate,
                        sgap,
                    ) {
                        break;
                    }
                    i += 1;
                }
                if i >= id_n && sc + (dp[m as usize] as u32 as i64) > max_sc {
                    max_sc = sc + (dp[m as usize] as u32 as i64);
                    p = m;
                }
                m -= 1;
            }
            let parent = if p >= 0 { p as u64 } else { 0xffff_ffffu64 };
            dp[k as usize] = (parent << 32) | (max_sc as u32 as u64);
            if long_sc < max_sc {
                long_sc = max_sc;
                long_idx = k;
            }
            k += 1;
        }

        let mut k = long_idx;
        while k >= 0 {
            win_a[k as usize] |= high_bit;
            let parent = dp[k as usize] >> 32;
            k = if parent != 0xffff_ffff {
                parent as i64
            } else {
                -1
            };
        }
        let mut write = 0usize;
        for read in 0..win_n as usize {
            if win_a[read] & high_bit == 0 {
                continue;
            }
            win_a[write] = (win_a[read] << 1) >> 1;
            write += 1;
        }
        win_n = write as i64;
    }

    let cov = -cov;
    for k in 0..win_n as usize {
        let ws = win_a[k];
        for i in 0..id_n as usize {
            let qn = aln[id_a[i] as usize].qn as usize;
            let wid = crate::correction::window::window_id_by_start(
                regions[qn].x_position_start as i32,
                regions[qn].x_position_end as i32,
                ws as i64,
                wl,
            )
            .0;
            let value = if cov > i64::from(i16::MIN) {
                cov
            } else {
                i64::from(i16::MIN)
            };
            regions[qn].window_cigar.windows[wid as usize].extra_end = value as i16;
            regions[qn].alignment_length = wid as u32;
        }
    }
    win_n as u64
}

/// `kv_pushp_cl(k_mer_hit, *cl, &p)`: write at the logical length
/// (`cl.length`), overwriting any stale scratch beyond it, then extend.
fn push_kmer_hit_at_length(cl: &mut CandidatesList, hit: KmerHit) {
    let index = cl.length as usize;
    if index < cl.list.len() {
        cl.list[index] = hit;
    } else {
        cl.list.push(hit);
    }
    cl.length = index as i64 + 1;
}

/// `push_adp_k_hits` (`Correct.cpp` 14758): push the anchor hits of one
/// `[qs, qs)` / `[ts, te)` bracket plus its two `readID = 0x7fffffff`
/// sentinels onto `cl`, returning the next scan position.
#[allow(clippy::too_many_arguments)]
fn push_adp_k_hits(
    cl: &mut CandidatesList,
    cln: i64,
    qs: u64,
    qe: u64,
    ts: u64,
    te: u64,
    read_id: u32,
    ci: i64,
    _debug_id: i64,
) -> i64 {
    let mut k = ci;
    let mut previous: i64 = -1;
    if cl.length > cln {
        previous = cl.length - 1;
        debug_assert!(is_alnw(&cl.list[previous as usize]));
    }
    while k < cln
        && cl.list[k as usize].read_identifier == read_id
        && (cl.list[k as usize].self_offset as u64) < qs
    {
        let hit = cl.list[k as usize];
        k += 1;
        if hit.offset as u64 >= ts {
            continue;
        }
        if previous >= 0 {
            let last = &cl.list[previous as usize];
            if last.offset >= hit.offset || last.self_offset >= hit.self_offset {
                continue;
            }
        }
        let mut pushed = hit;
        pushed.read_identifier = hit.count;
        push_kmer_hit_at_length(cl, pushed);
    }
    if cln == 6008 {
        eprintln!(
            "[PSH] ci_in={} k_out={} len={} m={}",
            ci,
            k,
            cl.length,
            previous
        );
    }
    if qs == u64::MAX || qe == u64::MAX {
        return k;
    }
    push_kmer_hit_at_length(
        cl,
        KmerHit {
            read_identifier: ALNW_SENTINEL,
            count: u32::MAX,
            reverse_strand: false,
            self_offset: qs as u32,
            offset: ts as u32,
        },
    );
    while k < cln
        && cl.list[k as usize].read_identifier == read_id
        && (cl.list[k as usize].self_offset as u64) < qe
    {
        k += 1;
    }
    push_kmer_hit_at_length(
        cl,
        KmerHit {
            read_identifier: ALNW_SENTINEL,
            count: u32::MAX,
            reverse_strand: true,
            self_offset: (qe - 1) as u32,
            offset: (te - 1) as u32,
        },
    );
    k
}

/// `gen_weight_khits0` (`Correct.cpp` 14799): label the hits of `[qs, qe)` with
/// the chain weight `dp`.
fn gen_weight_khits0(
    qs: u32,
    qe: u32,
    hits: &mut [KmerHit],
    count: i64,
    mut k: i64,
    dp: i64,
) -> i64 {
    // The reference scans `a[k]` with `k == an` after a full window; the
    // out-of-bounds read either stops immediately or walks back to the last
    // hit, and both branches leave `k` unchanged in effect — guard it here.
    while k >= 0
        && (k as usize) < hits.len()
        && hits[k as usize].self_offset >= qs
    {
        k -= 1;
    }
    if k < 0 {
        k = 0;
    }
    while k < count && hits[k as usize].self_offset < qe {
        if hits[k as usize].self_offset >= qs
            && hits[k as usize].self_offset < qe
            && !is_alnw(&hits[k as usize])
        {
            hits[k as usize].read_identifier = dp as u32;
        }
        k += 1;
    }
    k
}

/// `gen_weight_khits` (`Correct.cpp` 14810): sweep the `iidx` endpoints and
/// weight each anchor run by its coverage depth.
fn gen_weight_khits(idx: &[u64], hits: &mut [KmerHit]) {
    let count = hits.len() as i64;
    let mut k = 0i64;
    let mut depth = 0i64;
    let mut begin = 0i64;
    for entry in idx {
        let previous_depth = depth;
        let end;
        if ((entry >> 32) & 1) != 0 {
            depth -= 1;
            end = ((entry >> 33) + 1) as i64;
        } else {
            depth += 1;
            end = (entry >> 33) as i64;
        }
        if end > begin {
            k = gen_weight_khits0(begin as u32, end as u32, hits, count, k, previous_depth);
        }
        begin = end;
    }
}

/// `gen_cns_chain` (`Correct.cpp` 14847): turns one overlap's aligned windows
/// into a candidate anchor chain in `cl` (appended from the current length) and
/// returns the number of hits. `need_filter_khit` must be false (the UL
/// gap-filling driver passes 0), so `lchain_dp_trace` is not run.
pub fn gen_cns_chain(
    regions: &mut [OverlapRegion],
    region_index: usize,
    cl: &mut CandidatesList,
    iidx: Option<&[u64]>,
    _maximum_length_gap: i64,
    _small_gap_rate: f64,
    need_filter_khit: bool,
) -> i64 {
    let window_count = regions[region_index].window_cigar.windows.len() as i64;
    if window_count <= 0 {
        return 0;
    }
    {
        use std::sync::atomic::{AtomicU32, Ordering};
        static CALL: AtomicU32 = AtomicU32::new(0);
        let call = CALL.fetch_add(1, Ordering::Relaxed);
        if call == 0 {
            let region = &regions[region_index];
            eprintln!(
                "[GCE] y_id={} wn={} rcn={} seed={}",
                region.y_identifier,
                window_count,
                cl.length,
                region.shared_seed
            );
            for (k, m) in region
                .window_cigar
                .windows
                .iter()
                .enumerate()
                .take(300)
            {
                eprintln!(
                    "[GCEW] w={} x=[{},{}] y=[{},{}] anchor={}",
                    k,
                    m.x_start,
                    m.x_end,
                    m.y_start,
                    m.y_end,
                    (m.extra_end < 0 && m.y_end != -1) as i32
                );
            }
        }
    }
    let target_identifier = regions[region_index].y_identifier as i64;
    let mut query_start: i64 = -1;
    let mut query_end: i64 = -1;
    let mut target_start: i64 = -1;
    let mut target_end: i64 = -1;
    let mut chain_index = regions[region_index].shared_seed as i64;
    let read_identifier = cl.list[chain_index as usize].read_identifier;
    let chain_start = cl.length;
    let mut aligned = 0i64;

    for k in 0..window_count {
        let window = regions[region_index].window_cigar.windows[k as usize];
        if window.extra_end < 0 && window.y_end != -1 {
            if query_start == -1 {
                query_start = window.x_start as i64;
                target_start =
                    regions[region_index].window_cigar.windows[(k - 1) as usize].y_end as i64 + 1;
            }
            query_end = window.x_end as i64 + 1;
            target_end = window.y_end as i64 + 1;
        } else {
            if query_start != -1 {
                chain_index = push_adp_k_hits(
                    cl,
                    chain_start,
                    query_start as u64,
                    query_end as u64,
                    target_start as u64,
                    target_end as u64,
                    read_identifier,
                    chain_index,
                    target_identifier,
                );
                aligned += 1;
            }
            query_start = -1;
            query_end = -1;
            target_start = -1;
            target_end = -1;
        }
    }
    if query_start != -1 {
        chain_index = push_adp_k_hits(
            cl,
            chain_start,
            query_start as u64,
            query_end as u64,
            target_start as u64,
            target_end as u64,
            read_identifier,
            chain_index,
            target_identifier,
        );
        aligned += 1;
    }
    push_adp_k_hits(
        cl,
        chain_start,
        u64::MAX,
        u64::MAX,
        u64::MAX,
        u64::MAX,
        read_identifier,
        chain_index,
        target_identifier,
    );
    let hit_count = cl.length - chain_start;
    if let Some(iidx) = iidx {
        // The reference bounds the chain hits: `gen_weight_khits(iidx, ka, kn)`.
        gen_weight_khits(
            iidx,
            &mut cl.list[chain_start as usize..(chain_start + hit_count) as usize],
        );
    }
    let _ = (aligned, need_filter_khit);
    cl.length = chain_start;
    hit_count
}

/// `khit_long_gap` (`Correct.cpp` 14966): whether two anchors stay inside the
/// small-gap band. Returns `true` when the gap is acceptable.
#[allow(dead_code)]
fn khit_long_gap(a: &KmerHit, b: &KmerHit, small_bw_rate: f64, min_small_bw: i64) -> bool {
    let dq = b.self_offset as i64 - a.self_offset as i64;
    let dr = b.offset as i64 - a.offset as i64;
    let dd = if dq >= dr { dq - dr } else { dr - dq };
    let dm = if dq >= dr { dr } else { dq };
    !((dd as f64) > dm as f64 * small_bw_rate && dd > min_small_bw)
}

/// `filter_bad_khits` (`Correct.cpp` 14978): drop the leading/trailing anchors
/// that sit too far from the neighbouring bracket and compact the survivors.
/// `hits` is exactly the `a_n` entries of the reference.
#[allow(dead_code)]
fn filter_bad_khits(
    start: Option<&KmerHit>,
    end: Option<&KmerHit>,
    hits: &mut [KmerHit],
    small_bw_rate: f64,
    min_small_bw: i64,
) -> usize {
    let count = hits.len();
    let mut k = 0usize;
    let mut removed = 0usize;
    if let Some(start) = start {
        while k < count {
            let hit = hits[k];
            let bad = !khit_long_gap(start, &hit, small_bw_rate, min_small_bw)
                || end
                    .map(|end| !khit_long_gap(&hit, end, small_bw_rate, min_small_bw))
                    .unwrap_or(false);
            if !bad {
                break;
            }
            hits[k].offset = u32::MAX;
            hits[k].self_offset = u32::MAX;
            removed += 1;
            k += 1;
        }
    }
    if end.is_some() && k < count {
        let mut j = count as i64 - 1;
        while j >= 0 {
            let hit = hits[j as usize];
            if hit.count < hit.read_identifier {
                j -= 1;
                continue;
            }
            let bad = start
                .map(|start| !khit_long_gap(start, &hit, small_bw_rate, min_small_bw))
                .unwrap_or(false)
                || !khit_long_gap(&hit, end.unwrap(), small_bw_rate, min_small_bw);
            if !bad {
                break;
            }
            hits[j as usize].offset = u32::MAX;
            hits[j as usize].self_offset = u32::MAX;
            removed += 1;
            j -= 1;
        }
    }
    if removed != 0 {
        let mut write = 0usize;
        for read in 0..count {
            if hits[read].offset == u32::MAX {
                continue;
            }
            hits[write] = hits[read];
            write += 1;
        }
        return write;
    }
    count
}

/// `k_hits_bcheck` (`Correct.cpp` 14937): verify that the anchors of one run
/// share the same `khit`-base target prefix.
#[allow(dead_code)]
fn k_hits_bcheck(
    targets: &TargetSource<'_>,
    regions: &[OverlapRegion],
    candidates: &CandidatesList,
    khit: u64,
    a: &[u64],
    a_n: i64,
) -> bool {
    if a_n < 2 {
        return true;
    }
    let first = candidates.list[a[0] as u32 as usize];
    let Some(e) = (first.offset as u64).checked_sub(khit) else {
        return false;
    };
    let region0 = &regions[first.read_identifier as usize];
    let reference0 = targets.sub_region(
        region0.y_identifier as usize,
        e as i64,
        khit as i64,
        region0.y_strand as u8,
    );
    for k in 1..a_n as usize {
        let hit = candidates.list[a[k] as u32 as usize];
        let Some(e) = (hit.offset as u64).checked_sub(khit) else {
            return false;
        };
        let region = &regions[hit.read_identifier as usize];
        let reference = targets.sub_region(
            region.y_identifier as usize,
            e as i64,
            khit as i64,
            region.y_strand as u8,
        );
        if reference0 != reference {
            return false;
        }
    }
    true
}

/// `filter_adp_k_hits` (`Correct.cpp` 15017): copy the anchors of one bracket
/// into the compacted array at `cmi`, clip the ones too far from the
/// neighbouring brackets, refine them with [`crate::overlap::chaining::chain_refine`],
/// and return the next scan position.
#[allow(clippy::too_many_arguments)]
fn filter_adp_k_hits(
    anchors: &mut [KmerHit],
    anchor_count: i64,
    chain_data: &mut ChainData,
    query_start: u64,
    query_end: u64,
    target_start: u64,
    _target_end: u64,
    read_identifier: u32,
    cursor: i64,
    copy_index: &mut i64,
    previous_end: Option<&KmerHit>,
    small_bw_rate: f64,
    min_small_bw: i64,
) -> i64 {
    let mut k = cursor;
    let copy_start = *copy_index;
    while k < anchor_count
        && anchors[k as usize].read_identifier == read_identifier
        && (anchors[k as usize].self_offset as u64) < query_start
    {
        let hit = anchors[k as usize];
        k += 1;
        if hit.offset as u64 >= target_start {
            continue;
        }
        if let Some(previous) = previous_end {
            if previous.offset >= hit.offset || previous.self_offset >= hit.self_offset {
                continue;
            }
        }
        anchors[*copy_index as usize] = hit;
        *copy_index += 1;
    }
    let boundary = KmerHit {
        self_offset: query_start as u32,
        offset: target_start as u32,
        ..KmerHit::default()
    };
    let end_anchor: Option<&KmerHit> = if query_start == u64::MAX || query_end == u64::MAX {
        None
    } else {
        Some(&boundary)
    };
    if *copy_index > copy_start {
        let kept = filter_bad_khits(
            previous_end,
            end_anchor,
            &mut anchors[copy_start as usize..*copy_index as usize],
            small_bw_rate,
            min_small_bw,
        );
        *copy_index = copy_start + kept as i64;
        if *copy_index > copy_start {
            *copy_index = copy_start
                + crate::overlap::chaining::chain_refine(
                    &mut anchors[copy_start as usize..*copy_index as usize],
                    *copy_index - copy_start,
                    chain_data,
                    50,
                    5000,
                    512,
                    16,
                );
        }
    }
    while k < anchor_count
        && anchors[k as usize].read_identifier == read_identifier
        && (anchors[k as usize].self_offset as u64) < query_end
    {
        k += 1;
    }
    k
}

/// `refine_khits` (`Correct.cpp` 15046): rebuild one overlap's anchor runs.
pub fn refine_khits(
    region: &mut OverlapRegion,
    candidates: &mut CandidatesList,
    small_gap_rate: f64,
) {
    let window_count = region.window_cigar.windows.len() as i64;
    if window_count <= 0 {
        return;
    }
    let mut query_start: i64 = -1;
    let mut query_end: i64 = -1;
    let mut target_start: i64 = -1;
    let mut target_end: i64 = -1;
    let mut cursor = region.shared_seed as i64;
    let mut copy_index = cursor;
    let read_identifier = candidates.list[cursor as usize].read_identifier;
    let mut previous = KmerHit {
        self_offset: u32::MAX,
        offset: u32::MAX,
        ..KmerHit::default()
    };
    let mut aligned = 0i64;
    let minimum_gap = 64i64;

    for k in 0..window_count {
        let window = region.window_cigar.windows[k as usize];
        if window.extra_end < 0 && window.y_end != -1 {
            if query_start == -1 {
                query_start = window.x_start as i64;
                target_start =
                    region.window_cigar.windows[(k - 1) as usize].y_end as i64 + 1;
            }
            query_end = window.x_end as i64 + 1;
            target_end = window.y_end as i64 + 1;
        } else {
            if query_start != -1 {
                let previous_end = if aligned != 0 { Some(&previous) } else { None };
                cursor = filter_adp_k_hits(
                    &mut candidates.list,
                    candidates.length,
                    &mut candidates.chain_data,
                    query_start as u64,
                    query_end as u64,
                    target_start as u64,
                    target_end as u64,
                    read_identifier,
                    cursor,
                    &mut copy_index,
                    previous_end,
                    small_gap_rate,
                    minimum_gap,
                );
                previous.self_offset = (query_end - 1) as u32;
                previous.offset = (target_end - 1) as u32;
                aligned += 1;
            }
            query_start = -1;
            query_end = -1;
            target_start = -1;
            target_end = -1;
        }
    }
    if query_start != -1 {
        let previous_end = if aligned != 0 { Some(&previous) } else { None };
        cursor = filter_adp_k_hits(
            &mut candidates.list,
            candidates.length,
            &mut candidates.chain_data,
            query_start as u64,
            query_end as u64,
            target_start as u64,
            target_end as u64,
            read_identifier,
            cursor,
            &mut copy_index,
            previous_end,
            small_gap_rate,
            minimum_gap,
        );
        previous.self_offset = (query_end - 1) as u32;
        previous.offset = (target_end - 1) as u32;
        aligned += 1;
    }
    let _ = aligned;
    let previous_end = if aligned != 0 { Some(&previous) } else { None };
    cursor = filter_adp_k_hits(
        &mut candidates.list,
        candidates.length,
        &mut candidates.chain_data,
        u64::MAX,
        u64::MAX,
        u64::MAX,
        u64::MAX,
        read_identifier,
        cursor,
        &mut copy_index,
        previous_end,
        small_gap_rate,
        minimum_gap,
    );
    eprintln!(
        "[RK] seed={} cmi={} kept={}",
        region.shared_seed,
        copy_index,
        copy_index - region.shared_seed as i64
    );
    for k in copy_index..cursor {
        candidates.list[k as usize].read_identifier = ALNW_SENTINEL;
    }
}

/// `count_k_hits_adv` (`Correct.cpp` 15118): refine every overlap's anchors,
/// group them by query position, and label each position's support count
/// (`cnt`), with the `khit`-base check when `base_check` is set.
pub fn count_k_hits_adv(
    targets: &TargetSource<'_>,
    regions: &mut [OverlapRegion],
    candidates: &mut CandidatesList,
    index: &mut Vec<u64>,
    small_gap_rate: f64,
    khit: u64,
    base_check: bool,
) {
    let overlap_count = regions.len();
    let mut candidate_count = candidates.length;
    index.clear();
    for k in 0..overlap_count {
        index.push(((regions[k].shared_seed as u64) << 32) | k as u64);
        refine_khits(&mut regions[k], candidates, small_gap_rate);
    }
    radix_sort_u64(index);
    {
        use std::sync::atomic::{AtomicU32, Ordering};
        static CHECK: AtomicU32 = AtomicU32::new(0);
        let call = CHECK.fetch_add(1, Ordering::Relaxed);
        if call == 0 && overlap_count == 5 {
            for (k, region) in regions.iter().enumerate() {
                let mut cursor = region.shared_seed as i64;
                let pid = candidates.list[cursor as usize].read_identifier;
                let mut count = 0i64;
                while cursor < candidates.length
                    && candidates.list[cursor as usize].read_identifier == pid
                    && candidates.list[cursor as usize].read_identifier != ALNW_SENTINEL
                {
                    count += 1;
                    cursor += 1;
                }
                let first = &candidates.list[region.shared_seed as usize];
                let last = &candidates.list[(cursor - 1) as usize];
                eprintln!(
                    "[SLC] k={} seed={} pid={} n={} first=[{},{}] last=[{},{}]",
                    k,
                    region.shared_seed,
                    pid,
                    count,
                    first.self_offset,
                    first.offset,
                    last.self_offset,
                    last.offset
                );
            }
        }
    }

    let mut write = 0usize;
    for k in 0..overlap_count {
        let overlap = index[k] as u32 as usize;
        let mut cursor = regions[overlap].shared_seed as i64;
        let target_identifier = candidates.list[cursor as usize].read_identifier;
        regions[overlap].shared_seed = write as i32;
        while cursor < candidate_count
            && candidates.list[cursor as usize].read_identifier == target_identifier
            && candidates.list[cursor as usize].read_identifier != ALNW_SENTINEL
        {
            let mut hit = candidates.list[cursor as usize];
            hit.read_identifier = index[k] as u32;
            hit.count = 0;
            candidates.list[write] = hit;
            index.push(((hit.self_offset as u64) << 32) | write as u64);
            write += 1;
            cursor += 1;
        }
    }
    candidate_count = write as i64;
    candidates.length = candidate_count;

    let sorted_start = overlap_count;
    radix_sort_u64(&mut index[sorted_start..]);
    let sorted_count = index.len() - overlap_count;
    // Group the same self_offset together and set each hit's support count.
    let mut run_start = 0usize;
    for run_end in 1..=sorted_count {
        if run_end == sorted_count
            || (index[sorted_start + run_start] >> 32) != (index[sorted_start + run_end] >> 32)
        {
            let mut support = (run_end - run_start) as u32;
            if base_check && support > 1 {
                let ids = &index[sorted_start + run_start..sorted_start + run_end];
                let checked = k_hits_bcheck(targets, regions, candidates, khit, ids, ids.len() as i64);
                if !checked {
                    support = 1;
                }
            }
            for i in run_start..run_end {
                candidates.list[index[sorted_start + i] as u32 as usize].count = support;
            }
            run_start = run_end;
        }
    }
}

/// `update_overlap_region` (`Correct.cpp` 17249): copy `source`'s window cigar
/// into `target` and re-derive the overlap coordinates (including the same
/// left/right padding rule as `push_ovlp_chain_qgen`).
fn update_overlap_region(target: &mut OverlapRegion, source: &OverlapRegion, xl: i64, yl: i64) {
    target.window_cigar.cells = source.window_cigar.cells.clone();
    target.window_cigar.windows = source.window_cigar.windows.clone();
    let count = source.window_cigar.windows.len();
    if count != 0 {
        target.x_position_start = source.window_cigar.windows[0].x_start as u32;
        target.x_position_end = source.window_cigar.windows[count - 1].x_end as u32;
        target.y_position_start = source.window_cigar.windows[0].y_start as u32;
        target.y_position_end = source.window_cigar.windows[count - 1].y_end as u32;
    }
    if target.x_position_start <= target.y_position_start {
        target.y_position_start -= target.x_position_start;
        target.x_position_start = 0;
    } else {
        target.x_position_start -= target.y_position_start;
        target.y_position_start = 0;
    }
    let x_remainder = xl - target.x_position_end as i64 - 1;
    let y_remainder = yl - target.y_position_end as i64 - 1;
    if x_remainder <= y_remainder {
        target.x_position_end = (xl - 1) as u32;
        target.y_position_end += x_remainder as u32;
    } else {
        target.y_position_end = (yl - 1) as u32;
        target.x_position_end += y_remainder as u32;
    }
}

/// `cigar_gen_by_chain_adv` (`Correct.cpp` 17279): turn one overlap's
/// consensus chain (`ch_a`) into per-piece alignments, retry the unmatched
/// pieces with [`rechain_alignment_ul`], and write the final window cigar back
/// into `region` (`z`).
#[allow(clippy::too_many_arguments)]
pub fn cigar_gen_by_chain_adv(
    region: &mut OverlapRegion,
    candidates: &mut CandidatesList,
    chain_index: i64,
    chain_count: i64,
    overlaps: &mut [UlOv],
    window_length: i64,
    targets: &TargetSource<'_>,
    query: &[u8],
    extz: &mut BitExtz,
    auxiliary: &mut OverlapRegion,
    error_rate: f64,
    query_length: i64,
    h_kmer_hit: i64,
) {
    if overlaps.is_empty() {
        return;
    }
    let overlap_count = overlaps.len() as i64;
    let target_total_length = targets.length(region.y_identifier as usize);

    auxiliary.window_cigar.windows.clear();
    auxiliary.window_cigar.cells.clear();
    auxiliary.y_identifier = region.y_identifier;
    auxiliary.y_strand = region.y_strand;
    auxiliary.x_position_start = region.x_position_start;
    auxiliary.x_position_end = region.x_position_end;
    auxiliary.y_position_start = region.y_position_start;
    auxiliary.y_position_end = region.y_position_end;

    let entry_count;
    {
        let hits = &candidates.list[chain_index as usize..];
        entry_count = fusion_chain_ovlp(hits, chain_count, overlaps, overlap_count, query_length);
        for index in 0..entry_count as usize {
            let bracket = ChainBracket {
                qn: overlaps[index].qn,
                tn: overlaps[index].tn,
                qs: overlaps[index].qs,
                qe: overlaps[index].qe,
                ts: overlaps[index].ts,
                te: overlaps[index].te,
            };
            ovlp_base_aln(
                region,
                hits,
                chain_count,
                &bracket,
                window_length,
                targets,
                query,
                extz,
                &mut auxiliary.window_cigar,
                error_rate,
                query_length,
                target_total_length,
            );
        }
    }

    let auxiliary_count = auxiliary.window_cigar.windows.len();
    for index in 0..auxiliary_count {
        if !is_unaligned_window(&auxiliary.window_cigar.windows[index]) {
            continue;
        }
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
            query_length,
            target_total_length,
            h_kmer_hit,
        );
    }

    if auxiliary.window_cigar.windows.len() > auxiliary_count {
        let mut write = 0usize;
        for index in 0..auxiliary.window_cigar.windows.len() {
            if index < auxiliary_count && is_unaligned_window(&auxiliary.window_cigar.windows[index]) {
                continue;
            }
            auxiliary.window_cigar.windows[write] = auxiliary.window_cigar.windows[index];
            write += 1;
        }
        auxiliary.window_cigar.windows.truncate(write);
        auxiliary.window_cigar.windows.sort_by_key(|window| window.x_start);
    }

    update_overlap_region(region, auxiliary, query_length, target_total_length);
}

/// `ul_gap_filling_adv` (`Correct.cpp` 22139): the UL second-pass consensus
/// driver. Recomputes the anchor support, then for each overlap (group of
/// `cigar` entries sharing one `qn`) rebuilds the consensus chain and aligns
/// it through [`cigar_gen_by_chain_adv`].
#[allow(clippy::too_many_arguments)]
pub fn ul_gap_filling_adv(
    regions: &mut [OverlapRegion],
    candidates: &mut CandidatesList,
    cigar: &mut [UlOv],
    window_length: i64,
    targets: &TargetSource<'_>,
    query: &[u8],
    extz: &mut BitExtz,
    auxiliary: &mut OverlapRegion,
    buffer: &mut Vec<u64>,
    iidx: &[u64],
    error_rate: f64,
    query_length: i64,
    khit: u64,
    base_check: i64,
    max_length_gap: i64,
) {
    count_k_hits_adv(
        targets,
        regions,
        candidates,
        buffer,
        error_rate,
        khit,
        base_check != 0,
    );

    let entry_count = cigar.len();
    let mut group_start = 0usize;
    let mut cleared = 0usize;
    let mut k = 1usize;
    while k <= entry_count {
        if k == entry_count || cigar[group_start].qn != cigar[k].qn {
            let qn = cigar[group_start].qn as usize;
            for pk in cleared..qn {
                regions[pk].window_cigar.windows.clear();
            }
            cleared = qn + 1;
            let chain_count = gen_cns_chain(
                regions,
                qn,
                candidates,
                Some(iidx),
                max_length_gap,
                error_rate,
                false,
            );
            if query_length == 73261 {
                eprintln!(
                    "[CNS] qn={} ch_n={} rcn={} len={}",
                    qn,
                    chain_count,
                    candidates.length - chain_count,
                    candidates.length
                );
                for c in (candidates.length - chain_count)..candidates.length {
                    let hit = &candidates.list[c as usize];
                    eprintln!(
                        "[CNSC] c={} readID={} self_off={} off={}",
                        c,
                        hit.read_identifier,
                        hit.self_offset,
                        hit.offset
                    );
                }
            }
            if chain_count != 0 {
                let chain_index = candidates.length;
                cigar_gen_by_chain_adv(
                    &mut regions[qn],
                    candidates,
                    chain_index,
                    chain_count,
                    &mut cigar[group_start..k],
                    window_length,
                    targets,
                    query,
                    extz,
                    auxiliary,
                    error_rate,
                    query_length,
                    khit as i64,
                );
            }

            group_start = k;
        }
        k += 1;
    }
    for pk in cleared..regions.len() {
        regions[pk].window_cigar.windows.clear();
    }
}

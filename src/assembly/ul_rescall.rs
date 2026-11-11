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
use crate::correction::ec_driver::PafList;
use crate::correction::hc_align::{
    align_ul_ed_post_extz, gen_extend_err_exz, window_se_by_normalized_start, BitExtz,
    HcAlignParameters,
};
use crate::correction::interval::{get_interval, IntervalBuffer, WindowPool};
use crate::correction::recalculate::similarity_pass;
use crate::correction::region::THRESHOLD_MAXIMUM;
use crate::correction::target::TargetSource;
use crate::correction::window::number_of_windows;
use crate::overlap::region::OverlapRegion;
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
    pub sources: &'a [PafList],
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
fn check_connect_ug(graph: &Asg, v: u32, w: u32, bw: i64, diff_ec_ul: f64, dq: i64) -> bool {
    for arc in graph.arcs_of(v) {
        if arc.deleted || arc.v != w {
            continue;
        }
        return within_band(dq, i64::from(arc.overlap), bw, diff_ec_ul);
    }
    false
}

/// `check_connect_rg` (`Correct.cpp` 19190).
fn check_connect_rg(
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
    for hit in &context.sources[query].hits {
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

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

//! `break_ug_contig` and its helpers (`Overlaps.cpp` 11503–12565): breaks
//! unitigs at within-unitig read-coverage cliffs.
//!
//! `--b-cov` breaks where the depth of same-unitig overlaps drops below the
//! threshold, `--h-cov` where it rises above it; a candidate break must
//! additionally carry a fraction of exact (`el`) overlaps of at most
//! `--m-rate`. The two switches share the interval machinery but the
//! reference never enables both at once — every call site passes one and
//! `NULL` for the other.
//!
//! The 64-bit depth keys pack `<unitig offset>:rest | <not-circle>:1 |
//! <is-begin>:1`, so an ascending sort visits ends before begins at equal
//! offsets (`is-begin = 1` raises the key).

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::hit::ReadToUnitig;
use crate::assembly::unitig::{fill_unitig, renew_utg, UnitigGraph};
use crate::container::ksort::radix_sort_u64;
use crate::correction::ec_driver::PafList;
use crate::overlap::paf::PafHit;
use crate::reads::read_store::ReadStore;

/// `in_sub_t`: one depth interval of the sweep, `[k_beg, k_end)`.
#[derive(Debug, Clone, Copy, Default)]
struct InSubInterval {
    k_beg: u64,
    k_end: u64,
    dp: u32,
}

/// `push_cov_interval_direct`: one `[x_beg, x_end]` span as a begin/end key
/// pair; `None` when the span is outside the unitig.
fn push_cov_interval_direct(
    buffer: &mut Vec<u64>,
    x_beg: i64,
    x_end: i64,
    unitig_length: i64,
    is_circle: bool,
) -> bool {
    if x_beg <= x_end && x_beg >= 0 && x_end >= 0 && x_beg < unitig_length && x_end < unitig_length
    {
        let mut key = x_beg as u64;
        key <<= 1;
        key |= u64::from(!is_circle);
        key <<= 1;
        key |= 1;
        buffer.push(key);
        let mut key = (x_end + 1) as u64;
        key <<= 1;
        key |= u64::from(!is_circle);
        key <<= 1;
        buffer.push(key);
        true
    } else {
        false
    }
}

/// `push_cov_interval_advance`: places a span that may under/overflow the
/// unitig, wrapping it around for circular unitigs.
fn push_cov_interval_advance(
    buffer: &mut Vec<u64>,
    x_beg: i64,
    x_end: i64,
    unitig_length: i64,
    is_circle: bool,
) -> bool {
    if x_beg > x_end {
        return false;
    }

    if push_cov_interval_direct(buffer, x_beg, x_end, unitig_length, false) {
        return true;
    }

    if x_beg < 0 && x_end < 0 {
        return push_cov_interval_direct(
            buffer,
            unitig_length + x_beg,
            unitig_length + x_end,
            unitig_length,
            false,
        );
    }

    if x_beg >= unitig_length && x_end >= unitig_length {
        return push_cov_interval_direct(
            buffer,
            x_beg - unitig_length,
            x_end - unitig_length,
            unitig_length,
            false,
        );
    }

    if x_beg < 0 && x_end >= 0 {
        let x_beg = unitig_length + x_beg;
        if push_cov_interval_direct(buffer, x_beg, unitig_length - 1, unitig_length, is_circle)
            || push_cov_interval_direct(buffer, 0, x_end, unitig_length, is_circle)
        {
            return true;
        }
    }

    false
}

/// `get_break_point_cov_advance`: pushes the unitig-coordinate span of one
/// read.
///
/// The first pass (`first_pass`) tries every member read's own span once;
/// the second pass derives partner spans from the stored overlaps. The
/// `r_flag` bit 0 marks unitig membership and bit 2 marks spans already
/// placed. After the membership check the target id is deliberately reset
/// to the raw overlap target (`///must!!!!` in the reference), so the
/// `&2` dedup marks and reads raw target ids.
fn get_break_point_cov_advance(
    vertex: u32,
    unitig_offset: i64,
    read_graph: &Asg,
    reads: &ReadStore,
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    r_flag: &mut [u8],
    depth: &mut Vec<u64>,
    unitig_length: i64,
    is_circle: bool,
    first_pass: bool,
) {
    let read = (vertex >> 1) as usize;
    if first_pass
        && r_flag[read] & 1 != 0
        && r_flag[read] & 2 == 0
    {
        let v_beg = 0i64;
        let v_end = reads.read_length(read) as i64 - 1;
        if push_cov_interval_advance(
            depth,
            v_beg + unitig_offset,
            v_end + unitig_offset,
            unitig_length,
            is_circle,
        ) {
            r_flag[read] |= 2;
        }
        return;
    }

    for hit in &sources[read].hits {
        let resolved = resolve_target(read_graph, read_to_unitig, hit.target_identifier);
        let tn = match resolved {
            Some(tn) => tn,
            None => continue,
        };
        if r_flag[tn as usize] & 1 == 0 {
            continue;
        }
        // `tn = Get_tn((*h));///must!!!!` — back to the raw target.
        let tn = hit.target_identifier;
        if r_flag[tn as usize] & 2 != 0 {
            continue;
        }

        let query_length = reads.read_length(read) as i64;
        let target_length = reads.read_length(tn as usize) as i64;
        let mut ts = hit.target_start as i64;
        let mut te = hit.target_end as i64 - 1;
        if hit.reverse {
            ts = target_length - (hit.target_end as i64 - 1) - 1;
            te = target_length - hit.target_start as i64 - 1;
        }
        let qs = hit.query_start as i64;
        let qe = hit.query_end as i64 - 1;
        let ts = qs - ts;
        let te = qe + (target_length - te - 1);
        let mut w_beg = ts;
        let mut w_end = te;

        if vertex & 1 != 0 {
            let mirrored_begin = query_length - ts - 1;
            let mirrored_end = query_length - te - 1;
            w_beg = mirrored_end;
            w_end = mirrored_begin;
        }

        if push_cov_interval_advance(
            depth,
            w_beg + unitig_offset,
            w_end + unitig_offset,
            unitig_length,
            is_circle,
        ) {
            r_flag[tn as usize] |= 2;
        }
    }
}

/// The shared target check of the coverage/`m_rate` scans: deleted targets
/// resolve through the read-to-unitig map.
fn resolve_target(
    read_graph: &Asg,
    read_to_unitig: &ReadToUnitig,
    target: u32,
) -> Option<u32> {
    let mut tn = target;
    if read_graph.sequences[tn as usize].deleted {
        let (resolved, is_unitig) = read_to_unitig.get(tn);
        tn = resolved;
        if tn == u32::MAX || is_unitig == 1 || read_graph.sequences[tn as usize].deleted {
            return None;
        }
    }
    if read_graph.sequences[tn as usize].deleted {
        return None;
    }
    Some(tn)
}

/// `get_overlap_contig_dir`: which side of `[p_beg, p_end]` the overlap of
/// read vertex `v` (packed in `unitig.reads[k]`) reaches: 0 = before,
/// 1 = after, 2 = spanning.
fn get_overlap_contig_dir(
    vertex: u32,
    target: u32,
    hit: &PafHit,
    reads: &ReadStore,
    unitig_offset: i64,
    p_beg: u32,
    p_end: u32,
) -> usize {
    let target_length = reads.read_length(target as usize) as i64;
    let mut ts = hit.target_start as i64;
    let mut te = hit.target_end as i64 - 1;
    if hit.reverse {
        ts = target_length - (hit.target_end as i64 - 1) - 1;
        te = target_length - hit.target_start as i64 - 1;
    }
    let qs = hit.query_start as i64;
    let qe = hit.query_end as i64 - 1;
    let ts = qs - ts;
    let te = qe + (target_length - te - 1);
    let mut c_beg = ts;
    let mut c_end = te;

    if vertex & 1 != 0 {
        let query_length = reads.read_length((vertex >> 1) as usize) as i64;
        let mirrored_begin = query_length - ts - 1;
        let mirrored_end = query_length - te - 1;
        c_beg = mirrored_end;
        c_end = mirrored_begin;
    }

    let c_end = c_end + 1;
    let c_beg = c_beg + unitig_offset;
    let c_end = c_end + unitig_offset;

    if (p_beg != p_end && c_beg <= p_beg as i64 && c_end >= p_end as i64)
        || (p_beg == p_end && c_beg < p_beg as i64 && c_end > p_end as i64)
    {
        return 2;
    }

    if c_beg < p_beg as i64 {
        return 0;
    }
    if c_end > p_end as i64 {
        return 1;
    }

    2
}

/// `get_break_point_idx`: the read index a depth interval's break belongs
/// to, validated by the exact-overlap ratio. `u32::MAX` rejects the break.
fn get_break_point_idx(
    unitig: &MaUtg,
    reads: &ReadStore,
    r_flag: &[u8],
    sources: &[PafList],
    read_graph: &Asg,
    read_to_unitig: &ReadToUnitig,
    p_beg: u32,
    p_end: u32,
    exact_match_rate: f64,
) -> u32 {
    let mut min_k = u32::MAX;
    let mut min_l = 0u32;
    let mut e_occ: u32 = 0;
    let mut ne_occ: u32 = 0;

    let mut offset: u32 = 0;
    for (k, packed) in unitig.reads.iter().enumerate() {
        let c_beg = offset;
        let c_end = c_beg + reads.read_length((packed >> 33) as usize) as u32;
        offset = offset.wrapping_add(*packed as u32);

        let covers = (p_beg != p_end && c_beg <= p_beg && c_end >= p_end)
            || (p_beg == p_end && c_beg < p_beg && c_end > p_end);
        if covers {
            if min_k == u32::MAX {
                min_k = k as u32;
                min_l = c_beg;
            }
            for hit in &sources[(packed >> 33) as usize].hits {
                let tn = match resolve_target(read_graph, read_to_unitig, hit.target_identifier) {
                    Some(tn) => tn,
                    None => continue,
                };
                if r_flag[tn as usize] & 1 == 0 {
                    continue;
                }
                if hit.exact_region {
                    e_occ += 1;
                } else {
                    ne_occ += 1;
                }
            }
        } else if min_k != u32::MAX {
            break;
        }
    }

    if min_k == u32::MAX {
        return u32::MAX;
    }
    if e_occ as f64 <= (e_occ + ne_occ) as f64 * exact_match_rate {
        return c_break(unitig, reads, min_k, min_l, p_beg, p_end);
    }

    // Direction-aware second pass: per covered read, the exact/inexact
    // overlaps are split into "reaches before" / "reaches after" / "spans",
    // the spanning ones count for both sides, and the side with the higher
    // inexact rate accumulates into the totals.
    e_occ = 0;
    ne_occ = 0;
    let mut offset = min_l;
    for &packed in &unitig.reads[min_k as usize..] {
        let c_beg = offset;
        let c_end = c_beg + reads.read_length((packed >> 33) as usize) as u32;
        offset = offset.wrapping_add(packed as u32);

        let covers = (p_beg != p_end && c_beg <= p_beg && c_end >= p_end)
            || (p_beg == p_end && c_beg < p_beg && c_end > p_end);
        if covers {
            let mut e_occ_dir = [0.0f64; 3];
            let mut ne_occ_dir = [0.0f64; 3];

            for hit in &sources[(packed >> 33) as usize].hits {
                let tn = match resolve_target(read_graph, read_to_unitig, hit.target_identifier) {
                    Some(tn) => tn,
                    None => continue,
                };
                if r_flag[tn as usize] & 1 == 0 {
                    continue;
                }
                // `tn = Get_tn((*h));///must!!!!`
                let tn = hit.target_identifier;

                let index = get_overlap_contig_dir(
                    (packed >> 32) as u32,
                    tn,
                    hit,
                    reads,
                    c_beg as i64,
                    p_beg,
                    p_end,
                );
                if hit.exact_region {
                    e_occ_dir[index] += 1.0;
                } else {
                    ne_occ_dir[index] += 1.0;
                }
            }

            e_occ_dir[0] += e_occ_dir[2];
            e_occ_dir[1] += e_occ_dir[2];
            ne_occ_dir[0] += ne_occ_dir[2];
            ne_occ_dir[1] += ne_occ_dir[2];
            let rate_left = ne_occ_dir[0] / (ne_occ_dir[0] + e_occ_dir[0]);
            let rate_right = ne_occ_dir[1] / (ne_occ_dir[1] + e_occ_dir[1]);
            if rate_left >= rate_right {
                e_occ += e_occ_dir[0] as u32;
                ne_occ += ne_occ_dir[0] as u32;
            } else {
                e_occ += e_occ_dir[1] as u32;
                ne_occ += ne_occ_dir[1] as u32;
            }
        } else {
            break;
        }
    }

    if e_occ as f64 <= (e_occ + ne_occ) as f64 * exact_match_rate {
        return c_break(unitig, reads, min_k, min_l, p_beg, p_end);
    }
    u32::MAX
}

/// The `c_break` label of `get_break_point_idx`: pick the read nearer to the
/// interval and reject cuts at the ends of non-circular unitigs.
fn c_break(
    unitig: &MaUtg,
    reads: &ReadStore,
    min_k: u32,
    min_l: u32,
    p_beg: u32,
    p_end: u32,
) -> u32 {
    let k = min_k as usize;
    let c_beg = min_l;
    let c_end = c_beg + reads.read_length((unitig.reads[k] >> 33) as usize) as u32;
    // Unsigned subtraction, as in the reference: a negative underflow wraps
    // and fails the comparison.
    if p_beg.wrapping_sub(c_beg) <= c_end.wrapping_sub(p_end) {
        if k == 0 && !unitig.circular {
            return u32::MAX;
        }
        return k as u32;
    }
    let k = k as u32;
    if k + 1 < unitig.reads.len() as u32 {
        return k + 1;
    }
    if k + 1 == unitig.reads.len() as u32 && unitig.circular {
        return 0;
    }
    u32::MAX
}

/// `detect_break_point_advance`: builds the depth sweep of one unitig and
/// turns the transitions across the threshold into break-point candidates
/// `(unitig_id << 32) | read index`.
#[allow(clippy::too_many_arguments)]
fn detect_break_point_advance(
    unitig: &MaUtg,
    unitig_id: u32,
    read_graph: &Asg,
    reads: &ReadStore,
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    r_flag: &mut [u8],
    depth: &mut Vec<u64>,
    depth_intervals: &mut Vec<InSubInterval>,
    break_points: &mut Vec<u64>,
    low_coverage: Option<i32>,
    high_coverage: Option<i32>,
    exact_match_rate: f64,
) {
    if unitig.reads.len() < 2 {
        return;
    }
    r_flag.fill(0);
    for &packed in &unitig.reads {
        r_flag[(packed >> 33) as usize] = 1;
    }

    depth.clear();
    let unitig_length = unitig.length as i64;
    for first_pass in [true, false] {
        let mut offset = 0i64;
        for packed in &unitig.reads {
            let unitig_offset = offset;
            offset += (*packed as u32) as i64;
            get_break_point_cov_advance(
                (packed >> 32) as u32,
                unitig_offset,
                read_graph,
                reads,
                sources,
                read_to_unitig,
                r_flag,
                depth,
                unitig_length,
                unitig.circular,
                first_pass,
            );
        }
    }

    radix_sort_u64(depth);

    depth_intervals.clear();
    let mut dp: i64 = 0;
    let mut o_idx: u32 = 0;
    let mut o_dir: u32 = 1;
    for &key in depth.iter() {
        let o_dp = dp;
        // At equal offsets the end key sorts before the begin key.
        if key & 1 != 0 {
            dp += 1;
        } else {
            dp -= 1;
        }

        let dir = (key & 1) as u32;
        let idx = (key >> 2) as u32;
        if idx.wrapping_sub(o_idx) > 0 || (idx == o_idx && o_dir != dir) {
            if let Some(low) = low_coverage {
                if o_dp >= i64::from(low)
                    && !depth_intervals.is_empty()
                    && depth_intervals.last().unwrap().dp as i32 >= low
                {
                    let p = depth_intervals.last_mut().unwrap();
                    p.k_end = u64::from(idx);
                    p.dp = o_dp as u32;
                } else {
                    depth_intervals.push(InSubInterval {
                        k_beg: u64::from(o_idx),
                        k_end: u64::from(idx),
                        dp: o_dp as u32,
                    });
                }
            }

            if let Some(high) = high_coverage {
                if o_dp <= i64::from(high)
                    && !depth_intervals.is_empty()
                    && depth_intervals.last().unwrap().dp as i32 <= high
                {
                    let p = depth_intervals.last_mut().unwrap();
                    p.k_end = u64::from(idx);
                    p.dp = o_dp as u32;
                } else {
                    depth_intervals.push(InSubInterval {
                        k_beg: u64::from(o_idx),
                        k_end: u64::from(idx),
                        dp: o_dp as u32,
                    });
                }
            }
        }
        o_idx = idx;
        o_dir = dir;
    }

    if o_idx as i64 != unitig_length {
        depth_intervals.push(InSubInterval {
            k_beg: u64::from(o_idx),
            k_end: unitig_length as u64,
            dp: 0,
        });
    }

    // Valley/hill transitions: a run of intervals across the threshold, the
    // lowest-depth one of which decides the break. Intervals that reach the
    // array ends wrap around circular unitigs.
    let mut beg_idx = u32::MAX;
    let mut end_idx = u32::MAX;
    let mut cir_beg_idx = u32::MAX;
    let mut cir_end_idx = u32::MAX;
    let interval_count = depth_intervals.len() as u32;

    for k in 0..interval_count {
        if k > 0 {
            let previous = depth_intervals[k as usize - 1].dp as i32;
            let current = depth_intervals[k as usize].dp as i32;
            let crossing = match low_coverage {
                Some(low) => previous >= low && current < low,
                None => match high_coverage {
                    Some(high) => previous <= high && current > high,
                    None => false,
                },
            };
            if crossing {
                beg_idx = k;
            }
        }

        if k < interval_count - 1 {
            let current = depth_intervals[k as usize].dp as i32;
            let next = depth_intervals[k as usize + 1].dp as i32;
            let crossing = match low_coverage {
                Some(low) => current < low && next >= low,
                None => match high_coverage {
                    Some(high) => current > high && next <= high,
                    None => false,
                },
            };
            if crossing {
                end_idx = k;
                if beg_idx == u32::MAX {
                    cir_end_idx = k;
                }

                if beg_idx != u32::MAX && end_idx >= beg_idx {
                    let mut min = u32::MAX;
                    let mut min_idx = u32::MAX;
                    for i in beg_idx..=end_idx {
                        if depth_intervals[i as usize].dp < min {
                            let cur_idx = get_break_point_idx(
                                unitig,
                                reads,
                                r_flag,
                                sources,
                                read_graph,
                                read_to_unitig,
                                depth_intervals[i as usize].k_beg as u32,
                                depth_intervals[i as usize].k_end as u32,
                                exact_match_rate,
                            );
                            if cur_idx == u32::MAX {
                                continue;
                            }
                            min = depth_intervals[i as usize].dp;
                            min_idx = cur_idx;
                        }
                    }
                    if min_idx != u32::MAX {
                        break_points.push((u64::from(unitig_id) << 32) | u64::from(min_idx));
                    }
                }
                beg_idx = u32::MAX;
                end_idx = u32::MAX;
            }
        }
    }

    if beg_idx != u32::MAX && end_idx == u32::MAX {
        cir_beg_idx = beg_idx;
    }

    if unitig.circular && (cir_beg_idx != u32::MAX || cir_end_idx != u32::MAX) {
        let beg_idx = cir_beg_idx;
        let end_idx = cir_end_idx;

        let mut min = u32::MAX;
        let mut min_idx = u32::MAX;
        if beg_idx != u32::MAX {
            for i in beg_idx..interval_count {
                if depth_intervals[i as usize].dp < min {
                    let cur_idx = get_break_point_idx(
                        unitig,
                        reads,
                        r_flag,
                        sources,
                        read_graph,
                        read_to_unitig,
                        depth_intervals[i as usize].k_beg as u32,
                        depth_intervals[i as usize].k_end as u32,
                        exact_match_rate,
                    );
                    if cur_idx == u32::MAX {
                        continue;
                    }
                    min = depth_intervals[i as usize].dp;
                    min_idx = cur_idx;
                }
            }
        }

        if end_idx != u32::MAX {
            for i in 0..=end_idx {
                if depth_intervals[i as usize].dp < min {
                    let cur_idx = get_break_point_idx(
                        unitig,
                        reads,
                        r_flag,
                        sources,
                        read_graph,
                        read_to_unitig,
                        depth_intervals[i as usize].k_beg as u32,
                        depth_intervals[i as usize].k_end as u32,
                        exact_match_rate,
                    );
                    if cur_idx == u32::MAX {
                        continue;
                    }
                    min = depth_intervals[i as usize].dp;
                    min_idx = cur_idx;
                }
            }
        }

        if min_idx != u32::MAX {
            break_points.push((u64::from(unitig_id) << 32) | u64::from(min_idx));
        }
    }
}

/// `push_sub_unitig`: one split piece as a fresh unitig.
fn push_sub_unitig(
    new_unitigs: &mut Vec<MaUtg>,
    source: &MaUtg,
    read_graph: &Asg,
    edges: &[AsgArc],
    beg_idx: usize,
    occ: usize,
) {
    let mut unitig = MaUtg::default();
    unitig.circular = beg_idx == 0 && occ == source.reads.len() && source.circular;
    unitig.reads = source.reads[beg_idx..beg_idx + occ].to_vec();

    let total_length = fill_unitig(&mut unitig.reads, read_graph, edges, unitig.circular);
    unitig.length = total_length as u32;
    if !unitig.circular {
        unitig.start = (unitig.reads[0] >> 32) as u32;
        unitig.end = ((unitig.reads[unitig.reads.len() - 1] >> 32) as u32) ^ 1;
    } else {
        unitig.start = u32::MAX;
        unitig.end = u32::MAX;
    }
    new_unitigs.push(unitig);
}

/// `break_all_contigs`: rebuilds the unitig graph with the broken pieces and
/// remaps every arc through the old→new vertex table.
fn break_all_contigs(
    unitig_graph: &mut UnitigGraph,
    read_graph: &Asg,
    edges: &[AsgArc],
    break_points: &mut Vec<u64>,
) {
    unitig_graph.graph.cleanup();
    radix_sort_u64(break_points);
    break_points.dedup();

    let mut utg_idx = vec![u32::MAX; unitig_graph.unitigs.len() * 2];
    let mut new_unitigs: Vec<MaUtg> = Vec::new();
    let mut new_graph = Asg::default();

    let unitig_count = unitig_graph.unitigs.len() as u32;
    let mut group_start = 0usize;
    let mut next_untouched = 0u32;
    let mut k = 1usize;
    while k <= break_points.len() {
        if k == break_points.len() || (break_points[k] >> 32) != (break_points[group_start] >> 32)
        {
            let occ = k - group_start;
            let u_id = (break_points[group_start] >> 32) as u32;
            group_start = k;

            // The reference defers the carry-over when the group is empty or
            // the unitig is too short to break: those unitigs are copied by
            // a later group (or by the tail loop) instead.
            if occ == 0 || unitig_graph.unitigs[u_id as usize].reads.len() < 2 {
                k += 1;
                continue;
            }

            for s_i in next_untouched..u_id {
                let carried = unitig_graph.unitigs[s_i as usize].clone();
                let source = &mut unitig_graph.unitigs[s_i as usize];
                source.length = 0;
                source.circular = false;
                source.reads = Vec::new();

                utg_idx[(s_i << 1) as usize] = ((new_unitigs.len() as u32) << 1) as u32;
                utg_idx[((s_i << 1) + 1) as usize] = ((new_unitigs.len() as u32) << 1) as u32 + 1;
                new_unitigs.push(carried);
            }
            next_untouched = u_id + 1;

            utg_idx[((u_id << 1) + 1) as usize] = ((new_unitigs.len() as u32) << 1) as u32 + 1;

            let mut previous_index = 0usize;
            for piece in &break_points[k - occ..k] {
                let idx = *piece as u32 as usize;
                if !unitig_graph.unitigs[u_id as usize].circular && idx == 0 {
                    eprintln!("ERROR 1");
                    continue;
                }

                if unitig_graph.unitigs[u_id as usize].circular && idx == 0 {
                    // Un-circularise: drop the closing arcs of the reverse end.
                    let (arc_start, arc_count) = unitig_graph.graph.arc_range((u_id << 1) + 1);
                    let mut targets = Vec::new();
                    for arc_index in arc_start..arc_start + arc_count as usize {
                        if unitig_graph.graph.arcs[arc_index].deleted {
                            continue;
                        }
                        targets.push(
                            (unitig_graph.graph.arcs[arc_index].v, unitig_graph.graph.arcs[arc_index].ul >> 32),
                        );
                        unitig_graph.graph.arcs[arc_index].deleted = true;
                    }
                    for (target, query) in targets {
                        unitig_graph.graph.delete_arc(target ^ 1, (query as u32) ^ 1, true);
                    }
                    let source = &mut unitig_graph.unitigs[u_id as usize];
                    source.circular = false;
                    source.start = (source.reads[0] >> 32) as u32;
                    source.end = ((source.reads[source.reads.len() - 1] >> 32) as u32) ^ 1;
                    continue;
                }

                if idx.wrapping_sub(previous_index) <= 0 {
                    eprintln!(
                        "ERROR2: uID: {}, u->circ: {}, idx: {}, p_u_idx: {}",
                        u_id,
                        unitig_graph.unitigs[u_id as usize].circular as u32,
                        idx,
                        previous_index
                    );
                    continue;
                }

                push_sub_unitig(
                    &mut new_unitigs,
                    &unitig_graph.unitigs[u_id as usize],
                    read_graph,
                    edges,
                    previous_index,
                    idx - previous_index,
                );
                previous_index = idx;
            }
            push_sub_unitig(
                &mut new_unitigs,
                &unitig_graph.unitigs[u_id as usize],
                read_graph,
                edges,
                previous_index,
                unitig_graph.unitigs[u_id as usize].reads.len() - previous_index,
            );
            utg_idx[(u_id << 1) as usize] = ((new_unitigs.len() as u32 - 1) << 1) as u32;
        }
        k += 1;
    }

    for s_i in next_untouched..unitig_count {
        let carried = unitig_graph.unitigs[s_i as usize].clone();
        let source = &mut unitig_graph.unitigs[s_i as usize];
        source.length = 0;
        source.circular = false;
        source.reads = Vec::new();

        utg_idx[(s_i << 1) as usize] = ((new_unitigs.len() as u32) << 1) as u32;
        utg_idx[((s_i << 1) + 1) as usize] = ((new_unitigs.len() as u32) << 1) as u32 + 1;
        new_unitigs.push(carried);
    }

    // Remap the live arcs of the old graph into the new one.
    for arc in &unitig_graph.graph.arcs {
        if arc.deleted {
            continue;
        }
        let mut remapped = *arc;
        remapped.v = utg_idx[(arc.v ^ 1) as usize] ^ 1;
        let query_vertex = utg_idx[(arc.ul >> 32) as usize];
        remapped.ul = (remapped.ul << 32) >> 32;
        remapped.ul |= u64::from(query_vertex) << 32;
        new_graph.arcs.push(remapped);
    }

    for (k, unitig) in new_unitigs.iter().enumerate() {
        new_graph.set_sequence(k, unitig.length, false);
    }
    new_graph.cleanup();

    unitig_graph.graph = new_graph;
    unitig_graph.unitigs = new_unitigs;
    renew_utg(unitig_graph, read_graph, edges);
}

/// `break_ug_contig`: breaks every unitig at the coverage cliffs of its own
/// reads. Only one of the two switches is ever active, mirroring the
/// reference call sites.
pub fn break_ug_contig(
    unitig_graph: &mut UnitigGraph,
    read_graph: &Asg,
    reads: &ReadStore,
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    edges: &[AsgArc],
    low_coverage: Option<i32>,
    high_coverage: Option<i32>,
    exact_match_rate: f64,
) {
    if let Some(low) = low_coverage {
        eprintln!(
            "[M::break_ug_contig] break potential misassemblies with <{}-fold coverage",
            low
        );
    }
    if let Some(high) = high_coverage {
        eprintln!(
            "[M::break_ug_contig] break potential misassemblies with >{}-fold coverage",
            high
        );
    }

    let mut depth: Vec<u64> = Vec::new();
    let mut depth_intervals: Vec<InSubInterval> = Vec::new();
    let mut break_points: Vec<u64> = Vec::new();
    let mut primary_flag = vec![0u8; read_graph.sequences.len()];

    for i in 0..unitig_graph.unitigs.len() {
        let unitig = &unitig_graph.unitigs[i];
        if unitig.reads.is_empty() || unitig.reads.len() < 2 {
            continue;
        }
        detect_break_point_advance(
            unitig,
            i as u32,
            read_graph,
            reads,
            sources,
            read_to_unitig,
            &mut primary_flag,
            &mut depth,
            &mut depth_intervals,
            &mut break_points,
            low_coverage,
            high_coverage,
            exact_match_rate,
        );
    }

    break_all_contigs(unitig_graph, read_graph, edges, &mut break_points);
}

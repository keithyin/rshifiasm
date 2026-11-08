//! The read-to-unitig coverage chaining of the primary cleaning, the
//! counterparts of `set_utg_offset`, `chain_trans_ovlp`, `quick_LIS`,
//! `GetCeilIndex`, `get_xy_pos_by_pos`, `cmp_hap_alignment_chaining` and
//! `collect_trans_cov` in `Overlaps.cpp`/`Purge_Dups.cpp`.
//!
//! [`collect_trans_cov`] transfers coverage from an auxiliary unitig chain to
//! the primary one: [`set_utg_offset`] lays the auxiliary chain out as
//! `<offset>:32 | <read vertex>` positions, [`chain_trans_ovlp`] chains the
//! overlaps that cross from the primary chain into it, and the average
//! coverage of the aligned part is handed back to the primary reads. Every
//! read on both aligned parts is flagged `P_HET`.
//!
//! `chain_origin_trans_uid_by_distance` appends the read-to-read transitions
//! to the trans chain for the trio output stage. `collect_trans_ovlp` only
//! fills a different trio transition buffer and is still deferred.

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::hap_cov::{
    chain_origin_trans_uid_by_distance, ArcOffset, HapCov, P_HET, RC_1,
};
use crate::assembly::hit::{ma_hit_to_arc, CoverageCut};
use crate::correction::ec_driver::PafList;

/// `Get_xOff`: the high half of an [`ArcOffset`] offset.
pub fn get_x_offset(offset: u64) -> i64 {
    i64::from((offset >> 32) as u32)
}

/// `Get_yOff`: the low half of an [`ArcOffset`] offset.
pub fn get_y_offset(offset: u64) -> i64 {
    i64::from(offset as u32)
}

/// The vertex of the read at `index` when a unitig is walked in `orientation`,
/// the `v` computation shared by `set_utg_offset` and `chain_trans_ovlp`.
///
/// The reference flips the vertex's orientation bit with
/// `a[u->n - r_i - 1] ^ 0x100000000` on the reverse strand.
pub fn chain_read_vertex(unitig: &MaUtg, orientation: u32, index: usize) -> u32 {
    if orientation == 1 {
        ((unitig.reads[unitig.reads.len() - index - 1] ^ 0x1_0000_0000) >> 32) as u32
    } else {
        (unitig.reads[index] >> 32) as u32
    }
}

/// The read id at `index` when a unitig is walked in `orientation`, the
/// `rid` computation of `collect_trans_cov`. The orientation only reorders
/// the reads, it does not flip the id.
pub fn chain_read_id(unitig: &MaUtg, orientation: u32, index: usize) -> u32 {
    if orientation == 1 {
        (unitig.reads[unitig.reads.len() - index - 1] >> 33) as u32
    } else {
        (unitig.reads[index] >> 33) as u32
    }
}

/// `asg_arc_len` of the live arc from `from` to `to`, or zero and the
/// reference's error message when no such arc exists.
fn live_arc_length(read_graph: &Asg, from: u32, to: u32) -> u32 {
    let (start, count) = read_graph.arc_range(from);
    for i in 0..count {
        let arc = read_graph.arcs[start + i];
        if arc.deleted {
            continue;
        }
        if arc.v == to {
            return arc.node_length();
        }
    }
    eprintln!("ERROR-set_utg_offset");
    0
}

/// `set_utg_offset`: lays a unitig chain out along the read graph.
///
/// Every read on the chain stores `<cumulative length>:32 | <read vertex>` in
/// `pos_idx`, then `len` finally grows by the last read's length. With
/// `is_clear` set the entries are reset to `u64::MAX` instead and no arc is
/// looked up; with `only_len` set nothing is stored.
pub fn set_utg_offset(
    nodes: &[u32],
    unitigs: &[MaUtg],
    read_graph: &Asg,
    pos_idx: &mut [u64],
    is_clear: u32,
    only_len: u32,
) -> u32 {
    let mut len = 0u32;
    let mut previous = u32::MAX;

    for &node in nodes {
        let unitig = &unitigs[(node >> 1) as usize];
        let orientation = node & 1;
        if unitig.reads.is_empty() {
            continue;
        }

        for r_i in 0..unitig.reads.len() {
            let v = chain_read_vertex(unitig, orientation, r_i);
            let mut l = 0u32;
            if previous != u32::MAX && is_clear == 0 {
                l = live_arc_length(read_graph, previous, v);
            }
            previous = v;
            len = len.wrapping_add(l);
            if only_len != 0 {
                continue;
            }
            if is_clear == 1 {
                pos_idx[(v >> 1) as usize] = u64::MAX;
            } else {
                pos_idx[(v >> 1) as usize] = (u64::from(len) << 32) | u64::from(v);
            }
        }
    }

    if previous != u32::MAX {
        len = len.wrapping_add(read_graph.sequences[(previous >> 1) as usize].length);
    }

    len
}

/// `GetCeilIndex`: the first tail entry whose y position is `>= key`.
fn get_ceil_index(arr: &[ArcOffset], tail: &[i32], mut l: i32, mut r: i32, key: u32) -> i32 {
    while r - l > 1 {
        let m = l + (r - l) / 2;
        if get_y_offset(arr[tail[m as usize] as usize].offset) >= i64::from(key) {
            r = m;
        } else {
            l = m;
        }
    }
    r
}

/// `quick_LIS`: the longest increasing subsequence of the y positions.
///
/// `x` must already be sorted by x position. On return `tail_index` holds the
/// chain of indices into `x`, in increasing order; `prev_index` keeps the
/// predecessor of every element, sized like `x`.
pub fn quick_lis(x: &[ArcOffset], tail_index: &mut Vec<i32>, prev_index: &mut Vec<i32>) {
    tail_index.clear();
    prev_index.clear();
    let n = x.len();
    if n == 0 {
        return;
    }

    tail_index.resize(n, 0);
    prev_index.resize(n, 0);

    let mut len: i64 = 1;
    tail_index[0] = 0;
    prev_index[0] = -1;

    for i in 1..n as i64 {
        let index = i as usize;
        let y = get_y_offset(x[index].offset);
        if y < get_y_offset(x[tail_index[0] as usize].offset) {
            // A new smallest value; its predecessor is never read.
            tail_index[0] = i as i32;
        } else if y > get_y_offset(x[tail_index[(len - 1) as usize] as usize].offset) {
            prev_index[index] = tail_index[(len - 1) as usize];
            tail_index[len as usize] = i as i32;
            len += 1;
        } else {
            let pos = get_ceil_index(x, tail_index, -1, (len - 1) as i32, y as u32);
            prev_index[index] = if pos > 0 {
                tail_index[(pos - 1) as usize]
            } else {
                -1
            };
            tail_index[pos as usize] = i as i32;
        }
    }

    let mut i = tail_index[(len - 1) as usize];
    for m in 0..len {
        tail_index[(len - m - 1) as usize] = i;
        if m + 1 < len {
            i = prev_index[i as usize];
        }
    }
    tail_index.truncate(len as usize);
}

/// `get_xy_pos_by_pos`: the packed `<x position>:32 | <y position>:32` of one
/// crossing overlap, with the relative orientation written back into `rev`
/// (the reference stores it in the arc's `el` field).
#[allow(clippy::too_many_arguments)]
fn get_xy_pos_by_pos(
    read_graph: &Asg,
    arc: AsgArc,
    v_in_unitig: u32,
    w_in_unitig: u32,
    v_in_pos: u32,
    w_in_pos: u32,
    x_unitig_len: u32,
    y_unitig_len: u32,
    rev: &mut u8,
) -> u64 {
    let mut x_pos = u32::MAX;
    let mut y_pos = u32::MAX;
    let mut x_dir = 0u32;
    let mut y_dir = 0u32;

    let query_vertex = (arc.ul >> 32) as u32;
    if query_vertex == v_in_unitig {
        // End position of the read.
        x_pos = v_in_pos
            .wrapping_add(read_graph.sequences[(v_in_unitig >> 1) as usize].length)
            .wrapping_sub(1);
        x_dir = 0;
    } else if query_vertex == (v_in_unitig ^ 1) {
        // Start position of the read.
        x_pos = v_in_pos;
        x_dir = 1;
    } else {
        eprintln!("ERROR");
    }

    if arc.v == w_in_unitig {
        y_pos = w_in_pos.wrapping_add(arc.overlap).wrapping_sub(1);
        y_dir = 0;
    } else if arc.v == (w_in_unitig ^ 1) {
        y_pos = w_in_pos
            .wrapping_add(read_graph.sequences[(w_in_unitig >> 1) as usize].length)
            .wrapping_sub(arc.overlap);
        y_dir = 1;
    } else {
        eprintln!("ERROR");
    }

    *rev = (x_dir ^ y_dir) as u8;
    if *rev != 0 {
        if y_unitig_len <= y_pos {
            y_pos = u32::MAX;
        } else {
            y_pos = y_unitig_len.wrapping_sub(y_pos).wrapping_sub(1);
        }
    }

    if x_pos >= x_unitig_len {
        x_pos = u32::MAX;
    }
    if y_pos >= y_unitig_len {
        y_pos = u32::MAX;
    }

    (u64::from(x_pos) << 32) | u64::from(y_pos)
}

/// `cmp_hap_alignment_chaining`: orders by `x.el`, then x position, then y
/// position. The reference calls `qsort`, which is a stable merge sort in
/// glibc, so the port sorts stably.
pub fn compare_arc_offsets(left: &ArcOffset, right: &ArcOffset) -> std::cmp::Ordering {
    left.arc
        .exact_region
        .cmp(&right.arc.exact_region)
        .then_with(|| (left.offset >> 32).cmp(&(right.offset >> 32)))
        .then_with(|| (left.offset as u32).cmp(&(right.offset as u32)))
}

/// `chain_trans_ovlp`: chains the overlaps that cross from `x_reads` into the
/// unitig chain [`set_utg_offset`] has laid out in `cov.pos_idx`.
///
/// Returns the occurrence in `x_reads` where the alignment reaches its end,
/// zero when no chain was found. The reference `utg_trans_t *o` variant only
/// fills the trio transition buffers, so this port implements the
/// `o == NULL` path the primary cleaning uses.
#[allow(clippy::too_many_arguments)]
pub fn chain_trans_ovlp(
    cov: &mut HapCov,
    unitigs: &[MaUtg],
    read_graph: &Asg,
    reverse_sources: &[PafList],
    coverage_cut: &[CoverageCut],
    x_reads: &[u32],
    target_base_len: u32,
    maximum_hang_rate: f32,
) -> u32 {
    let max_hang = cov.max_hang;
    let min_ovlp = cov.min_ovlp;
    let pos_idx = &mut cov.pos_idx;
    let u_buffer = &mut cov.u_buffer;
    let tail_index = &mut cov.tail_index;
    let prev_index = &mut cov.prev_index;

    u_buffer.clear();
    let mut x_end = 0u32;

    // First pass: how many reads of the primary chain the auxiliary chain
    // covers before its own length is exhausted.
    let mut a_occ = 0u32;
    let mut len = 0u32;
    let mut previous = u32::MAX;
    let mut x_occ = u32::MAX;
    'first: for &node in x_reads {
        let unitig = &unitigs[(node >> 1) as usize];
        let orientation = node & 1;
        if unitig.reads.is_empty() {
            continue;
        }
        for r_i in 0..unitig.reads.len() {
            let v = chain_read_vertex(unitig, orientation, r_i);
            let mut l = 0u32;
            if previous != u32::MAX {
                l = live_arc_length(read_graph, previous, v);
            }
            previous = v;
            len = len.wrapping_add(l);
            if len >= target_base_len {
                x_occ = a_occ;
                break 'first;
            }
            a_occ += 1;
        }
    }
    if x_occ == u32::MAX {
        x_occ = a_occ;
    }
    if x_occ == 0 {
        x_occ = 1;
    }

    // Second pass: collect every crossing overlap of the covered reads.
    let mut a_occ = 0u32;
    let mut len = 0u32;
    let mut previous = u32::MAX;
    for &node in x_reads {
        let unitig = &unitigs[(node >> 1) as usize];
        let orientation = node & 1;
        if unitig.reads.is_empty() {
            continue;
        }
        for r_i in 0..unitig.reads.len() {
            if a_occ >= x_occ {
                break;
            }
            let v = chain_read_vertex(unitig, orientation, r_i);
            let mut l = 0u32;
            if previous != u32::MAX {
                l = live_arc_length(read_graph, previous, v);
            }
            previous = v;
            len = len.wrapping_add(l);

            for hit in &reverse_sources[(v >> 1) as usize].hits {
                let query = hit.query_identifier as usize;
                let target = hit.target_identifier as usize;
                if coverage_cut[target].deleted || read_graph.sequences[target].deleted {
                    continue;
                }
                let query_length = coverage_cut[query].end.wrapping_sub(coverage_cut[query].start) as i32;
                let target_length = coverage_cut[target].end.wrapping_sub(coverage_cut[target].start) as i32;
                let (result, mut arc) = ma_hit_to_arc(
                    hit,
                    query_length,
                    target_length,
                    max_hang,
                    maximum_hang_rate,
                    min_ovlp,
                );
                // A contained or short overlap cannot be chained.
                if result < 0 {
                    continue;
                }

                let target_read = (arc.v >> 1) as usize;
                if read_graph.sequences[target_read].deleted {
                    continue;
                }
                if pos_idx[target_read] == u64::MAX {
                    continue;
                }
                let w = pos_idx[target_read] as u32;
                if target_read != (w >> 1) as usize {
                    continue;
                }

                let tmp = get_xy_pos_by_pos(
                    read_graph,
                    arc,
                    v,
                    w,
                    len,
                    (pos_idx[target_read] >> 32) as u32,
                    u32::MAX,
                    target_base_len,
                    &mut arc.exact_region,
                );
                if (tmp >> 32) == u64::from(u32::MAX) || (tmp as u32) == u32::MAX {
                    continue;
                }
                // `t.el` now holds the relative orientation; reverse matches
                // cannot join the chain.
                if arc.exact_region != 0 {
                    continue;
                }

                u_buffer.push(ArcOffset {
                    arc,
                    offset: tmp,
                    weight: 1,
                });
            }

            a_occ += 1;
        }
        if a_occ >= x_occ {
            break;
        }
    }

    if u_buffer.is_empty() {
        return x_end;
    }

    u_buffer.sort_by(compare_arc_offsets);

    // Merge the equal `(el, Off)` entries, keeping the longest overlap.
    let count = u_buffer.len();
    let mut k = 1usize;
    let mut l = 0usize;
    let mut m = 0usize;
    while k <= count {
        if k == count
            || u_buffer[k].arc.exact_region != u_buffer[l].arc.exact_region
            || u_buffer[k].offset != u_buffer[l].offset
        {
            u_buffer[m] = u_buffer[l];
            let mut merge = l + 1;
            while merge < k {
                u_buffer[m].weight += u_buffer[merge].weight;
                if u_buffer[merge].arc.overlap > u_buffer[m].arc.overlap {
                    u_buffer[m].arc = u_buffer[merge].arc;
                }
                merge += 1;
            }
            l = k;
            m += 1;
        }
        k += 1;
    }
    u_buffer.truncate(m);

    quick_lis(u_buffer, tail_index, prev_index);
    if tail_index.is_empty() {
        return x_end;
    }

    // Third pass: walk the primary chain until the alignment's end.
    let best = u_buffer[*tail_index.last().expect("the chain is not empty") as usize];
    let mut x_len_thres = u32::MAX;
    let mut a_occ = 0u32;
    let mut len = 0u32;
    let mut previous = u32::MAX;
    for &node in x_reads {
        let unitig = &unitigs[(node >> 1) as usize];
        let orientation = node & 1;
        if unitig.reads.is_empty() {
            continue;
        }
        for r_i in 0..unitig.reads.len() {
            let v = chain_read_vertex(unitig, orientation, r_i);
            let mut l = 0u32;
            if previous != u32::MAX {
                l = live_arc_length(read_graph, previous, v);
            }
            previous = v;
            len = len.wrapping_add(l);

            if (v >> 1) == (best.arc.ul >> 33) as u32 && x_len_thres == u32::MAX {
                for hit in &reverse_sources[(v >> 1) as usize].hits {
                    let query = hit.query_identifier as usize;
                    let target = hit.target_identifier as usize;
                    if coverage_cut[target].deleted || read_graph.sequences[target].deleted {
                        continue;
                    }
                    let query_length =
                        coverage_cut[query].end.wrapping_sub(coverage_cut[query].start) as i32;
                    let target_length =
                        coverage_cut[target].end.wrapping_sub(coverage_cut[target].start) as i32;
                    let (result, mut arc) = ma_hit_to_arc(
                        hit,
                        query_length,
                        target_length,
                        max_hang,
                        maximum_hang_rate,
                        min_ovlp,
                    );
                    if result < 0 {
                        continue;
                    }

                    let target_read = (arc.v >> 1) as usize;
                    if read_graph.sequences[target_read].deleted {
                        continue;
                    }
                    if pos_idx[target_read] == u64::MAX {
                        continue;
                    }
                    let w = pos_idx[target_read] as u32;
                    if target_read != (w >> 1) as usize {
                        continue;
                    }

                    let tmp = get_xy_pos_by_pos(
                        read_graph,
                        arc,
                        v,
                        w,
                        len,
                        (pos_idx[target_read] >> 32) as u32,
                        u32::MAX,
                        target_base_len,
                        &mut arc.exact_region,
                    );
                    if (tmp >> 32) == u64::from(u32::MAX) || (tmp as u32) == u32::MAX {
                        continue;
                    }
                    if arc.exact_region != 0 {
                        continue;
                    }

                    if tmp == best.offset
                        && arc.v == best.arc.v
                        && arc.ul == best.arc.ul
                    {
                        x_len_thres = (get_x_offset(best.offset) + i64::from(target_base_len)
                            - get_y_offset(best.offset))
                            as u32;
                    }
                }
            }

            if len >= x_len_thres {
                return a_occ;
            }
            a_occ += 1;
        }
    }

    x_end = a_occ;
    x_end
}

/// `collect_trans_cov`: moves the coverage of the auxiliary chain onto the
/// primary one and flags both sides `P_HET`.
///
/// `primary`/`auxiliary` are unitig vertex chains; `primary_offset` and
/// `auxiliary_offset` are the overlap lengths the reference passes through to
/// [`chain_origin_trans_uid_by_distance`].
#[allow(clippy::too_many_arguments)]
pub fn collect_trans_cov(
    primary: &[u32],
    primary_offset: u64,
    auxiliary: &[u32],
    auxiliary_offset: u64,
    unitigs: &[MaUtg],
    read_graph: &Asg,
    coverage_cut: &[CoverageCut],
    reverse_sources: &[PafList],
    maximum_hang_rate: f32,
    cov: &mut HapCov,
) {
    let _ = (primary_offset, auxiliary_offset);
    if primary.is_empty() || auxiliary.is_empty() {
        return;
    }

    let len_aux = set_utg_offset(auxiliary, unitigs, read_graph, &mut cov.pos_idx, 0, 0);
    let thre_pri = chain_trans_ovlp(
        cov,
        unitigs,
        read_graph,
        reverse_sources,
        coverage_cut,
        primary,
        len_aux,
        maximum_hang_rate,
    );

    if thre_pri > 0 {
        // Append the read-to-read transitions of this alignment; they feed
        // the trio/haplotype output stage.
        if let Some(t_ch) = cov.t_ch.as_mut() {
            chain_origin_trans_uid_by_distance(
                t_ch,
                unitigs,
                &mut cov.pos_idx,
                &mut cov.u_buffer,
                &mut cov.tail_index,
                primary,
                primary_offset as u32,
                None,
                auxiliary,
                auxiliary_offset as u32,
                Some(len_aux),
                RC_1,
                -1024.0,
            );
        }

        let has_trans_chain = cov.t_ch.is_some();

        // The auxiliary side contributes all of its coverage.
        let mut u_cov = 0u64;
        for &node in auxiliary {
            let unitig = &unitigs[(node >> 1) as usize];
            let orientation = node & 1;
            if unitig.reads.is_empty() {
                continue;
            }
            for k in 0..unitig.reads.len() {
                let read = chain_read_id(unitig, orientation, k) as usize;
                u_cov += u64::from(cov.cov[read]);
                if has_trans_chain {
                    cov.is_r_het[read] |= P_HET;
                }
            }
        }

        // The primary side only up to the alignment end.
        let mut u_len = 0u64;
        let mut occ = 0u32;
        'primary_len: for &node in primary {
            let unitig = &unitigs[(node >> 1) as usize];
            let orientation = node & 1;
            if unitig.reads.is_empty() {
                continue;
            }
            for k in 0..unitig.reads.len() {
                if occ >= thre_pri {
                    break;
                }
                let read = chain_read_id(unitig, orientation, k) as usize;
                u_len += u64::from(read_graph.sequences[read].length);
                if has_trans_chain {
                    cov.is_r_het[read] |= P_HET;
                }
                occ += 1;
            }
            if occ >= thre_pri {
                break 'primary_len;
            }
        }

        let u_cov = if u_len == 0 { 0 } else { u_cov / u_len };

        let mut occ = 0u32;
        'primary_add: for &node in primary {
            let unitig = &unitigs[(node >> 1) as usize];
            let orientation = node & 1;
            if unitig.reads.is_empty() {
                continue;
            }
            for k in 0..unitig.reads.len() {
                if occ >= thre_pri {
                    break;
                }
                let read = chain_read_id(unitig, orientation, k) as usize;
                let add = u_cov.wrapping_mul(u64::from(read_graph.sequences[read].length)) as u32;
                cov.cov[read] = cov.cov[read].wrapping_add(add);
                occ += 1;
            }
            if occ >= thre_pri {
                break 'primary_add;
            }
        }
    }

    set_utg_offset(auxiliary, unitigs, read_graph, &mut cov.pos_idx, 1, 0);
}

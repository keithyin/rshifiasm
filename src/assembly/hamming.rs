//! Hamming-error bubble reduction and phase-bubble refinement, the counterpart
//! of `rd_hm_bub`, `rd_hm_drop`, `rd_hamming_symm_simple0`,
//! `gen_switch_phasing`, `gen_fg`, `recall_arcs`, `renew_phase_bubble`,
//! `refine_bubble_popping` and `reduce_hamming_error_adv` in `Overlaps.cpp`.
//!
//! The C++ code keeps the phase graph of the fly-simplification in
//! `rd_hamming_fly_simp_t`: the cached unitig graph `fg`, the bubble-node list
//! `srt`, the arc marks `vs`, the inserted-edge counter and the parameters of
//! the overlap scan. Here that becomes [`FlySimp`], which owns a copy of the
//! graph and borrows the overlap lists and the caller's edge buffer.

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::bubble::{BubbleBuffer, BubbleInfo};
use crate::assembly::bubble_identify::{asg_bub_backtrack_primary, asg_bub_pop1_primary_trio};
use crate::assembly::gfa_clean::get_real_length;
use crate::assembly::graph_build::asg_arc_del_trans;
use crate::assembly::hap_cov::{get_num_trio_flag, HapCov};
use crate::assembly::hit::{ma_hit_to_arc, CoverageCut, ALTER_LABEL, HAP_LABEL, PRIMARY_LABEL};
use crate::assembly::trio_clean::DROP;
use crate::assembly::unitig::{
    copy_untig_graph, get_edge_from_source, ma_ug_gen_primary, UnitigGraph, AMBIGU, FATHER, MOTHER,
};
use crate::container::ksort::radix_sort_by_key;
use crate::correction::ec_driver::PafList;
use crate::logging;

/// `Uc_beg`: the read vertex at the beginning of a unitig.
fn uc_beg(unitig: &MaUtg) -> u32 {
    (unitig.reads[0] >> 32) as u32
}

/// `Uc_end`: the read vertex at the end of a unitig.
fn uc_end(unitig: &MaUtg) -> u32 {
    ((unitig.reads[unitig.reads.len() - 1] >> 32) as u32) ^ 1
}

/// The read vertex leaving unitig vertex `z`.
fn read_out(unitigs: &[MaUtg], z: u32) -> u32 {
    let unitig = &unitigs[(z >> 1) as usize];
    if z & 1 != 0 {
        uc_beg(unitig) ^ 1
    } else {
        uc_end(unitig) ^ 1
    }
}

/// The read vertex entering unitig vertex `z`.
fn read_in(unitigs: &[MaUtg], z: u32) -> u32 {
    let unitig = &unitigs[(z >> 1) as usize];
    if z & 1 != 0 {
        uc_end(unitig)
    } else {
        uc_beg(unitig)
    }
}

/// `rd_hamming_fly_simp_t`: the cached phase graph and its scratch state.
///
/// `graph` is the reference's `fg`: a copy of the unitig graph that receives
/// the edges recovered by the phase refinement. `srt` is `bc`, the bubble-node
/// list of [`renew_phase_bubble`]; `vs` is `s`, the per-vertex mark array;
/// `edges` is the `new_rtg_edges` buffer every recovered edge is appended to.
pub struct FlySimp<'a> {
    /// `fg`.
    pub graph: UnitigGraph,
    /// `srt`/`bc`.
    pub srt: Vec<u32>,
    /// `vs`.
    pub vs: Vec<u8>,
    /// `n_insert`.
    pub n_insert: u32,
    /// `src`, the forward overlap lists.
    pub sources: &'a [PafList],
    /// `cov`.
    pub coverage_cut: &'a [CoverageCut],
    /// `max_hang`.
    pub max_hang: i32,
    /// `asm_opt.max_hang_rate`, read by `ma_hit2arc`.
    pub max_hang_rate: f32,
    /// `min_ovlp`.
    pub min_ovlp: i32,
    /// `gap_fuzz`.
    pub gap_fuzz: i32,
    /// `ae`, the polished-edge buffer of the caller.
    pub edges: &'a mut Vec<AsgArc>,
}

/// `gen_rd_hamming_fly_simp_t`: builds the cached phase graph.
#[allow(clippy::too_many_arguments)]
pub fn gen_rd_hamming_fly_simp_t<'a>(
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
    sources: &'a [PafList],
    coverage_cut: &'a [CoverageCut],
    max_hang: i32,
    max_hang_rate: f32,
    min_ovlp: i32,
    gap_fuzz: i32,
    edges: &'a mut Vec<AsgArc>,
) -> FlySimp<'a> {
    let graph = gen_fg(
        unitig_graph,
        read_graph,
        sources,
        coverage_cut,
        max_hang,
        max_hang_rate,
        min_ovlp,
        gap_fuzz,
    );
    let vertex_count = unitig_graph.graph.sequences.len() << 1;
    FlySimp {
        graph,
        srt: Vec::new(),
        vs: vec![0u8; vertex_count],
        n_insert: 0,
        sources,
        coverage_cut,
        max_hang,
        max_hang_rate,
        min_ovlp,
        gap_fuzz,
        edges,
    }
}

/// `recall_arcs`: appends to `des` every live arc of `src` that `des` lost.
fn recall_arcs(des: &mut Asg, src: &Asg) {
    let vertex_count = src.sequences.len() * 2;
    let mut missing: Vec<AsgArc> = Vec::new();

    for v in 0..vertex_count as u32 {
        if src.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        let (source_start, source_count) = src.arc_range(v);
        let (target_start, target_count) = des.arc_range(v);
        for zi in 0..source_count {
            let za = src.arcs[source_start + zi];
            if za.deleted {
                continue;
            }
            let w = za.v;
            let mut ai = 0usize;
            while ai < target_count {
                let av = des.arcs[target_start + ai];
                if !av.deleted && av.v == w {
                    break;
                }
                ai += 1;
            }
            if ai >= target_count {
                missing.push(za);
            }
        }
    }

    let transitive = missing.len() as u32;
    if !missing.is_empty() {
        for arc in missing {
            des.push_arc(arc);
        }
        des.index.clear();
        des.is_sorted = false;
        des.cleanup();
        // The reference does not symmetrise here.
    }
    eprintln!("[M::recall_arcs] # transitive arcs::{}", transitive);
    eprintln!(
        "[M::recall_arcs] # new arcs::{}, # old arcs::{}",
        des.arcs.len(),
        src.arcs.len()
    );
}

/// `gen_fg`: copies the unitig graph and adds the phase edges implied by the
/// reads shared between neighbouring unitigs.
///
/// `srt` packs `<unitig vertex>:32 | <read vertex>:32` for both orientations
/// of every non-circular unitig and is sorted by read vertex; `idx` maps a read
/// to the first entry of its run.
#[allow(clippy::too_many_arguments)]
fn gen_fg(
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    max_hang: i32,
    max_hang_rate: f32,
    min_ovlp: i32,
    gap_fuzz: i32,
) -> UnitigGraph {
    let mut index = vec![u32::MAX; read_graph.sequences.len()];
    let mut fg = copy_untig_graph(unitig_graph);
    fg.graph.cleanup();

    let mut srt: Vec<u64> = Vec::new();
    for k in 0..fg.unitigs.len() {
        fg.graph.sequences[k].coverage_class = PRIMARY_LABEL;
        let unitig = &unitig_graph.unitigs[k];
        if unitig.circular || unitig.reads.is_empty() {
            continue;
        }
        let m = ((k as u64) << 1) | (u64::from(uc_beg(unitig)) << 32);
        srt.push(m);
        let m = (((k as u64) << 1) + 1) | (u64::from(uc_end(unitig)) << 32);
        srt.push(m);
    }

    radix_sort_by_key(&mut srt, &|value: &u64| *value);
    let mut l = 0usize;
    for k in 1..=srt.len() {
        if k == srt.len() || (srt[k] >> 33) != (srt[l] >> 33) {
            index[(srt[l] >> 33) as usize] = l as u32;
            l = k;
        }
    }

    let mut n_insert = 0u64;
    for k in 0..fg.unitigs.len() {
        let unitig = &unitig_graph.unitigs[k];
        if unitig.circular || unitig.reads.is_empty() {
            continue;
        }

        let uv = (k as u32) << 1;
        let rv = uc_end(unitig) ^ 1;
        n_insert += gen_fg_orientation(
            &mut fg,
            unitig_graph,
            read_graph,
            sources,
            coverage_cut,
            &srt,
            &index,
            uv,
            rv,
            max_hang,
            max_hang_rate,
            min_ovlp,
        );

        let uv = ((k as u32) << 1) + 1;
        let rv = uc_beg(unitig) ^ 1;
        n_insert += gen_fg_orientation(
            &mut fg,
            unitig_graph,
            read_graph,
            sources,
            coverage_cut,
            &srt,
            &index,
            uv,
            rv,
            max_hang,
            max_hang_rate,
            min_ovlp,
        );
    }

    if n_insert > 0 {
        fg.graph.index.clear();
        fg.graph.is_sorted = false;
        fg.graph.cleanup();
        fg.graph.symmetrize();
        asg_arc_del_trans(&mut fg.graph, gap_fuzz);
        // Some old edges may be lost due to the transitive reduction.
        recall_arcs(&mut fg.graph, &unitig_graph.graph);
    }

    fg
}

/// One orientation of the [`gen_fg`] scan, the body of the reference's two
/// identical loops.
#[allow(clippy::too_many_arguments)]
fn gen_fg_orientation(
    fg: &mut UnitigGraph,
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    srt: &[u64],
    index: &[u32],
    uv: u32,
    rv: u32,
    max_hang: i32,
    max_hang_rate: f32,
    min_ovlp: i32,
) -> u64 {
    let mut n_insert = 0u64;
    let (za_start, za_count) = unitig_graph.graph.arc_range(uv);

    for hit in &sources[(rv >> 1) as usize].hits {
        let sq = &coverage_cut[hit.query_identifier as usize];
        let st = &coverage_cut[hit.target_identifier as usize];
        if st.deleted || read_graph.sequences[hit.target_identifier as usize].deleted {
            continue;
        }
        let (result, t0) = ma_hit_to_arc(
            hit,
            sq.end.wrapping_sub(sq.start) as i32,
            st.end.wrapping_sub(st.start) as i32,
            max_hang,
            max_hang_rate,
            min_ovlp,
        );
        if result < 0 {
            continue;
        }
        if t0.query_vertex() != rv {
            continue;
        }
        let rw = t0.v;
        let first = index[(rw >> 1) as usize];
        if first == u32::MAX {
            continue;
        }
        let t0 = match get_edge_from_source(
            sources,
            coverage_cut,
            max_hang,
            max_hang_rate,
            min_ovlp,
            t0.query_vertex(),
            t0.v,
        ) {
            Some(arc) => arc,
            None => continue,
        };

        let mut m = first as usize;
        assert_eq!((srt[m] >> 33) as u32, rw >> 1);
        while m < srt.len() && (srt[m] >> 33) as u32 == (rw >> 1) {
            if rw == (srt[m] >> 32) as u32 {
                let uw = srt[m] as u32;
                if uv == uw {
                    m += 1;
                    continue;
                }
                let mut z = 0usize;
                while z < za_count {
                    let arc = unitig_graph.graph.arcs[za_start + z];
                    if !arc.deleted && arc.v == uw {
                        break;
                    }
                    z += 1;
                }
                if z < za_count {
                    m += 1;
                    continue;
                }
                if let Some(t1) = get_edge_from_source(
                    sources,
                    coverage_cut,
                    max_hang,
                    max_hang_rate,
                    min_ovlp,
                    t0.v ^ 1,
                    t0.query_vertex() ^ 1,
                ) {
                    let mut forward = t0;
                    forward.ul = (forward.ul & 0xffff_ffff) | (u64::from(uv) << 32);
                    forward.v = uw;
                    fg.graph.push_arc(forward);

                    let mut reverse = t1;
                    reverse.ul = (reverse.ul & 0xffff_ffff) | (u64::from(uw ^ 1) << 32);
                    reverse.v = uv ^ 1;
                    fg.graph.push_arc(reverse);
                    n_insert += 1;
                }
            }
            m += 1;
        }
    }

    n_insert
}

/// `rd_hm_bub`: the bubble search used by the hamming-error reduction.
///
/// Unlike [`asg_bub_pop1_primary_trio`], the search may follow arcs of a second
/// `reference` graph, and the first vertex of the search sees zero-length
/// edges. The sink is left in `buffer.sources[0]` on success and the visited
/// nodes in `buffer.visited`.
fn rd_hm_bub(
    graph: &Asg,
    reference: Option<&Asg>,
    v0: u32,
    max_distance: u64,
    buffer: &mut BubbleBuffer,
) -> u64 {
    if graph.sequences[(v0 >> 1) as usize].deleted {
        return 0;
    }
    let outgoing0 = get_real_length(graph, v0);
    let outgoing1 = reference.map_or(0, |reference| get_real_length(reference, v0));
    if outgoing0 + outgoing1 < 2 {
        return 0;
    }

    buffer.sources.clear();
    buffer.visited_edges.clear();
    buffer.visited.clear();
    {
        let info = &mut buffer.info[v0 as usize];
        info.positive_count = 0;
        info.distance = 0;
        info.negative_count = 0;
        info.node_count = 0;
        info.non_positive_count = 0;
    }
    buffer.sources.push(v0);

    let mut n_pending = 0u32;
    let mut is_first = true;
    let mut n_tips = 0u32;
    let mut tip_end = u32::MAX;
    let mut success = false;

    'search: loop {
        let v = buffer
            .sources
            .pop()
            .expect("the loop condition keeps the stack non-empty");
        let distance = buffer.info[v as usize].distance;
        let mut loops_complete = true;

        let (start, count) = graph.arc_range(v);
        let mut i = 0usize;
        while i < count {
            let arc = graph.arcs[start + i];
            i += 1;
            if arc.deleted {
                continue;
            }
            let w = arc.v;
            let mut length = arc.node_length();
            if (w >> 1) == (v0 >> 1) {
                break 'search;
            }
            if is_first {
                length = 0;
            }
            let candidate = distance.wrapping_add(length);
            if u64::from(candidate) > max_distance {
                loops_complete = false;
                break;
            }

            if !buffer.info[w as usize].state {
                buffer.visited.push(w);
                let remaining = get_real_length(graph, w ^ 1)
                    + reference.map_or(0, |reference| get_real_length(reference, w ^ 1));
                let info = &mut buffer.info[w as usize];
                info.parent = v;
                info.state = true;
                info.distance = candidate;
                info.remaining = remaining;
                n_pending += 1;
            } else {
                let info = &mut buffer.info[w as usize];
                if candidate < info.distance {
                    info.distance = candidate;
                }
            }

            let remaining = {
                let info = &mut buffer.info[w as usize];
                info.remaining = info.remaining.wrapping_sub(1);
                info.remaining
            };
            if remaining == 0 {
                let outgoing = get_real_length(graph, w)
                    + reference.map_or(0, |reference| get_real_length(reference, w));
                if outgoing > 0 {
                    buffer.sources.push(w);
                } else if n_tips != 0 {
                    break 'search;
                } else {
                    n_tips += 1;
                    tip_end = w;
                }
                n_pending = n_pending.wrapping_sub(1);
            }
        }

        if loops_complete {
            if let Some(reference) = reference {
                let (start, count) = reference.arc_range(v);
                let mut i = 0usize;
                while i < count {
                    let arc = reference.arcs[start + i];
                    i += 1;
                    if arc.deleted {
                        continue;
                    }
                    let w = arc.v;
                    let mut length = arc.node_length();
                    if (w >> 1) == (v0 >> 1) {
                        break 'search;
                    }
                    if is_first {
                        length = 0;
                    }
                    let candidate = distance.wrapping_add(length);
                    if u64::from(candidate) > max_distance {
                        loops_complete = false;
                        break;
                    }

                    if !buffer.info[w as usize].state {
                        buffer.visited.push(w);
                        let remaining = get_real_length(graph, w ^ 1)
                            + get_real_length(reference, w ^ 1);
                        let info = &mut buffer.info[w as usize];
                        info.parent = v;
                        info.state = true;
                        info.distance = candidate;
                        info.remaining = remaining;
                        n_pending += 1;
                    } else {
                        let info = &mut buffer.info[w as usize];
                        if candidate < info.distance {
                            info.distance = candidate;
                        }
                    }

                    let remaining = {
                        let info = &mut buffer.info[w as usize];
                        info.remaining = info.remaining.wrapping_sub(1);
                        info.remaining
                    };
                    if remaining == 0 {
                        let outgoing =
                            get_real_length(graph, w) + get_real_length(reference, w);
                        if outgoing > 0 {
                            buffer.sources.push(w);
                        } else if n_tips != 0 {
                            break 'search;
                        } else {
                            n_tips += 1;
                            tip_end = w;
                        }
                        n_pending = n_pending.wrapping_sub(1);
                    }
                }
            }
        }

        is_first = false;
        if n_tips == 1 {
            if tip_end != u32::MAX && n_pending == 0 && buffer.sources.is_empty() {
                buffer.sources.push(tip_end);
                success = true;
                break 'search;
            }
            break 'search;
        }

        if !loops_complete || buffer.sources.is_empty() {
            break 'search;
        }
        if !(buffer.sources.len() > 1 || n_pending > 0) {
            success = true;
            break 'search;
        }
    }

    for i in 0..buffer.visited.len() {
        let node = buffer.visited[i];
        let info = &mut buffer.info[node as usize];
        info.state = false;
        info.positive_count = 0;
        info.distance = 0;
        info.negative_count = 0;
        info.node_count = 0;
        info.non_positive_count = 0;
    }

    if success {
        1
    } else {
        0
    }
}

/// `rd_hm_drop0`: drops the arcs of `v` that cannot be part of the hamming
/// bubble.
fn rd_hm_drop0(
    graph: &mut Asg,
    reference: &Asg,
    v: u32,
    cutoff: f64,
    drop_inexact: bool,
) -> u32 {
    let (start, count) = graph.arc_range(v);

    if drop_inexact {
        let mut n_cut = 0u32;
        let mut i = 0usize;
        while i < count {
            let arc = graph.arcs[start + i];
            i += 1;
            if arc.deleted || arc.exact_region != 0 {
                continue;
            }
            graph.arcs[start + i - 1].deleted = true;
            graph.delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, true);
            n_cut += 1;
        }
        return n_cut;
    }

    if cutoff < 1.0 {
        let mut maximum_overlap = 0u32;
        for k in 0..count {
            let arc = graph.arcs[start + k];
            if !arc.deleted && maximum_overlap < arc.overlap {
                maximum_overlap = arc.overlap;
            }
        }
        let (reference_start, reference_count) = reference.arc_range(v);
        for k in 0..reference_count {
            let arc = reference.arcs[reference_start + k];
            if !arc.deleted && maximum_overlap < arc.overlap {
                maximum_overlap = arc.overlap;
            }
        }

        let mut n_cut = 0u32;
        if maximum_overlap > 0 {
            let threshold = f64::from(maximum_overlap) * cutoff;
            let mut i = 0usize;
            while i < count {
                let arc = graph.arcs[start + i];
                i += 1;
                if arc.deleted {
                    continue;
                }
                if f64::from(arc.overlap) < threshold {
                    graph.arcs[start + i - 1].deleted = true;
                    graph.delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, true);
                    n_cut += 1;
                }
            }
        }
        return n_cut;
    }

    let mut n_cut = 0u32;
    let mut i = 0usize;
    while i < count {
        let arc = graph.arcs[start + i];
        i += 1;
        if arc.deleted {
            continue;
        }
        graph.arcs[start + i - 1].deleted = true;
        graph.delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, true);
        n_cut += 1;
    }
    n_cut
}

/// `rd_hm_drop`: walks `reference` from `v0` and drops the arcs of every
/// visited vertex, including the two ends of the bubble.
///
/// `buffer.visited` is *appended* to, never cleared: the reference relies on
/// the union of the last failed bubble search and this walk.
fn rd_hm_drop(
    graph: &mut Asg,
    reference: &Asg,
    v0: u32,
    v1: u32,
    cutoff: f64,
    drop_inexact: bool,
    buffer: &mut BubbleBuffer,
) -> u32 {
    if graph.sequences[(v0 >> 1) as usize].deleted {
        return 0;
    }

    buffer.sources.clear();
    buffer.sources.push(v0);
    while let Some(v) = buffer.sources.pop() {
        if buffer.info[v as usize].state {
            continue;
        }
        buffer.info[v as usize].state = true;
        buffer.visited.push(v);

        let (start, count) = reference.arc_range(v);
        for k in 0..count {
            let arc = reference.arcs[start + k];
            if arc.deleted {
                continue;
            }
            let w = arc.v;
            if buffer.info[w as usize].state || w == v1 {
                continue;
            }
            buffer.sources.push(w);
        }
    }

    let mut n_cut = 0u32;
    for i in 0..buffer.visited.len() {
        let v = buffer.visited[i];
        buffer.info[v as usize].state = false;
        if v == v0 || v == v1 {
            continue;
        }
        n_cut += rd_hm_drop0(graph, reference, v, cutoff, drop_inexact);
        n_cut += rd_hm_drop0(graph, reference, v ^ 1, cutoff, drop_inexact);
    }
    n_cut += rd_hm_drop0(graph, reference, v0, cutoff, drop_inexact);
    n_cut += rd_hm_drop0(graph, reference, v1 ^ 1, cutoff, drop_inexact);
    n_cut
}

/// `rd_hamming_symm_simple0`: fixes the bubble between `st` and `ed^1`.
///
/// Returns whether the bubble was closed and the distance limit the reference
/// reports back through `r_max_dist` (the recomputed one when the first search
/// failed).
fn rd_hamming_symm_simple0(
    buffer: &mut BubbleBuffer,
    reference: &Asg,
    graph: &mut Asg,
    st: u32,
    ed: u32,
    max_distance: u64,
) -> (bool, u64) {
    let mut max_dist = max_distance;

    let mut popped = rd_hm_bub(graph, Some(reference), st, max_dist, buffer);
    if popped != 0 {
        assert_eq!(buffer.sources[0], ed ^ 1);
        return (true, max_dist);
    }

    // Recompute the distance limit from the reference-only bubble.
    popped = rd_hm_bub(reference, None, st, max_dist, buffer);
    assert!(popped != 0);
    assert_eq!(buffer.sources[0], ed ^ 1);
    max_dist = 0;
    for k in 0..buffer.visited.len() {
        let node = buffer.visited[k];
        if node == st || node == buffer.sources[0] {
            continue;
        }
        max_dist += u64::from(reference.sequences[(node >> 1) as usize].length);
    }
    max_dist += u64::from(reference.sequences[(st >> 1) as usize].length);
    max_dist += u64::from(reference.sequences[(buffer.sources[0] >> 1) as usize].length);

    popped = rd_hm_bub(graph, Some(reference), st, max_dist, buffer);
    if popped != 0 {
        assert_eq!(buffer.sources[0], ed ^ 1);
        return (true, max_dist);
    }

    // Drop the inexact edges first.
    let _ = rd_hm_drop(graph, reference, st, ed ^ 1, -1.0, true, buffer);
    popped = rd_hm_bub(graph, Some(reference), st, max_dist, buffer);
    if popped != 0 {
        assert_eq!(buffer.sources[0], ed ^ 1);
        return (true, max_dist);
    }

    const STEP: f64 = 0.2;
    let mut cutoff = STEP;
    while cutoff < 1.0 {
        let n_cut = rd_hm_drop(graph, reference, st, ed ^ 1, cutoff, false, buffer);
        popped = rd_hm_bub(graph, Some(reference), st, max_dist, buffer);
        if popped != 0 {
            assert_eq!(buffer.sources[0], ed ^ 1);
            return (true, max_dist);
        }
        if n_cut == 0 {
            break;
        }
        cutoff += STEP;
    }

    rd_hm_drop(graph, reference, st, ed ^ 1, 1024.0, false, buffer);
    popped = rd_hm_bub(graph, Some(reference), st, max_dist, buffer);
    assert!(popped != 0);
    assert_eq!(buffer.sources[0], ed ^ 1);
    (false, max_dist)
}

/// `asg_arc_del_trans_aux`: transitive reduction of `graph` against the
/// `auxiliary` graph.
///
/// The mark array deliberately keeps the stale marks of the auxiliary graph
/// between vertices, exactly like the reference.
fn asg_arc_del_trans_aux(
    graph: &mut Asg,
    auxiliary: &Asg,
    mark: &mut [u8],
    fuzz: i32,
) -> u32 {
    let vertex_count = graph.sequences.len() * 2;
    for value in mark.iter_mut().take(vertex_count) {
        *value = 0;
    }

    let mut n_reduced = 0u32;
    for v in 0..vertex_count as u32 {
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        let (graph_start, graph_count) = graph.arc_range(v);
        let (auxiliary_start, auxiliary_count) = auxiliary.arc_range(v);
        if graph_count + auxiliary_count == 0 {
            continue;
        }

        let mut live = 0u32;
        for i in 0..graph_count {
            let arc = graph.arcs[graph_start + i];
            if arc.deleted {
                continue;
            }
            mark[arc.v as usize] = 1;
            live += 1;
        }
        for i in 0..auxiliary_count {
            let arc = auxiliary.arcs[auxiliary_start + i];
            if arc.deleted {
                continue;
            }
            mark[arc.v as usize] = 1;
            live += 1;
        }
        if live == 0 {
            continue;
        }

        // The last slot of each span is the longest outgoing edge.
        let mut limit = 0u32;
        if graph_count > 0 {
            limit = graph.arcs[graph_start + graph_count - 1].node_length();
        }
        if auxiliary_count > 0 {
            let candidate = auxiliary.arcs[auxiliary_start + auxiliary_count - 1].node_length();
            if limit < candidate {
                limit = candidate;
            }
        }
        let limit = limit.wrapping_add(fuzz as u32);

        for i in 0..graph_count {
            let arc = graph.arcs[graph_start + i];
            let w = arc.v;
            if mark[w as usize] != 1 {
                continue;
            }
            let (w_start, w_count) = graph.arc_range(w);
            let (auxiliary_w_start, auxiliary_w_count) = auxiliary.arc_range(w);
            let arc_length = arc.node_length();

            let mut j = 0usize;
            while j < w_count
                && graph.arcs[w_start + j].node_length().wrapping_add(arc_length) <= limit
            {
                let target = graph.arcs[w_start + j].v;
                if mark[target as usize] != 0 {
                    mark[target as usize] = 2;
                }
                j += 1;
            }
            let mut j = 0usize;
            while j < auxiliary_w_count
                && auxiliary.arcs[auxiliary_w_start + j]
                    .node_length()
                    .wrapping_add(arc_length)
                    <= limit
            {
                let target = auxiliary.arcs[auxiliary_w_start + j].v;
                if mark[target as usize] != 0 {
                    mark[target as usize] = 2;
                }
                j += 1;
            }
        }

        for i in 0..graph_count {
            let target = graph.arcs[graph_start + i].v as usize;
            if mark[target] == 2 {
                graph.arcs[graph_start + i].deleted = true;
                n_reduced += 1;
            }
            mark[target] = 0;
        }
    }

    if n_reduced != 0 {
        graph.cleanup();
        graph.symmetrize();
    }
    n_reduced
}

/// `asg_bub_backtrack_check_switch`: marks the bubble reads that had no
/// assignment as `DROP` and reports whether both parents appear on the path.
fn asg_bub_backtrack_check_switch(
    unitigs: &[MaUtg],
    v0: u32,
    buffer: &BubbleBuffer,
    trio_labels: &mut [u8],
) -> u32 {
    let sink = buffer.sources[0];
    let mut father_occurrence = 0u32;
    let mut mother_occurrence = 0u32;

    let mut v = sink;
    loop {
        let u = buffer.info[v as usize].parent;
        if v != sink {
            let unitig = &unitigs[(v >> 1) as usize];
            for &packed in &unitig.reads {
                let read = (packed >> 33) as usize;
                if trio_labels[read] == FATHER {
                    father_occurrence += 1;
                }
                if trio_labels[read] == MOTHER {
                    mother_occurrence += 1;
                }
                if trio_labels[read] != AMBIGU {
                    continue;
                }
                trio_labels[read] = DROP as u8;
            }
        }
        v = u;
        if v == v0 {
            break;
        }
    }

    if father_occurrence > 0 && mother_occurrence > 0 {
        1
    } else {
        0
    }
}

/// `asg_bub_pop1_primary_trio_switch_check`: [`asg_bub_pop1_primary_trio`]
/// without the zero-length first edges and with the switch test on success.
///
/// The second return value is `Some(is_switch)` only when a bubble was found;
/// the reference then leaves the caller's `is_switch` untouched.
pub fn asg_bub_pop1_primary_trio_switch_check(
    graph: &mut Asg,
    unitigs: &[MaUtg],
    v0: u32,
    max_distance: u64,
    buffer: &mut BubbleBuffer,
    positive_flag: u32,
    negative_flag: u32,
    is_pop: bool,
    trio_labels: &mut [u8],
) -> (u64, Option<u32>) {
    let sequence = graph.sequences[(v0 >> 1) as usize];
    if sequence.deleted || sequence.coverage_class == ALTER_LABEL {
        return (0, None);
    }
    if get_real_length(graph, v0) < 2 {
        return (0, None);
    }

    buffer.sources.clear();
    buffer.visited.clear();
    buffer.visited_edges.clear();
    {
        let info = &mut buffer.info[v0 as usize];
        info.positive_count = 0;
        info.distance = 0;
        info.negative_count = 0;
        info.node_count = 0;
        info.non_positive_count = 0;
    }
    buffer.sources.push(v0);

    let mut n_pending = 0u32;
    let mut n_tips = 0u32;
    let mut tip_end = u32::MAX;
    let mut non_positive_flag = u32::MAX;
    if positive_flag == u32::from(FATHER) {
        non_positive_flag = u32::from(MOTHER);
    }
    if positive_flag == u32::from(MOTHER) {
        non_positive_flag = u32::from(FATHER);
    }

    let mut success = false;
    'search: loop {
        let v = buffer
            .sources
            .pop()
            .expect("the loop condition keeps the stack non-empty");
        let parent = buffer.info[v as usize];
        let distance = parent.distance;
        let positive = parent.positive_count;
        let negative = parent.negative_count;
        let non_positive = parent.non_positive_count;
        let (start, count) = graph.arc_range(v);
        let mut too_far = false;

        let mut i = 0usize;
        while i < count {
            let arc_index = start + i;
            let arc = graph.arcs[arc_index];
            i += 1;
            if arc.deleted {
                continue;
            }

            let w = arc.v;
            let length = arc.node_length();
            if (w >> 1) == (v0 >> 1) {
                break 'search;
            }

            buffer.visited_edges.push(arc_index as u32);
            if u64::from(distance.wrapping_add(length)) > max_distance {
                too_far = true;
                break;
            }

            let cur_positive = get_num_trio_flag(unitigs, graph, w >> 1, positive_flag, trio_labels);
            let cur_negative = get_num_trio_flag(unitigs, graph, w >> 1, negative_flag, trio_labels);
            let cur_non_positive = if non_positive_flag != u32::MAX {
                get_num_trio_flag(unitigs, graph, w >> 1, non_positive_flag, trio_labels)
            } else {
                0
            };
            let cur_node_count = unitigs[(w >> 1) as usize].reads.len() as u32;

            if !buffer.info[w as usize].state {
                buffer.visited.push(w);
                let info = &mut buffer.info[w as usize];
                info.parent = v;
                info.state = true;
                info.distance = distance.wrapping_add(length);
                info.positive_count = positive.wrapping_add(cur_positive);
                info.negative_count = negative.wrapping_add(cur_negative);
                info.non_positive_count = non_positive.wrapping_add(cur_non_positive);
                info.node_count = parent.node_count.wrapping_add(cur_node_count);
                info.remaining = get_real_length(graph, w ^ 1);
                n_pending += 1;
            } else {
                let info = &mut buffer.info[w as usize];
                let current_weight = i64::from(positive.wrapping_add(cur_positive))
                    - (i64::from(negative.wrapping_add(cur_negative))
                        + i64::from(non_positive.wrapping_add(cur_non_positive)));
                let maximum_weight = i64::from(info.positive_count)
                    - (i64::from(info.negative_count) + i64::from(info.non_positive_count));
                let candidate_count = parent.node_count.wrapping_add(cur_node_count);
                let candidate_distance = distance.wrapping_add(length);

                let replace = if current_weight > maximum_weight {
                    true
                } else if current_weight == maximum_weight {
                    if candidate_count > info.node_count {
                        true
                    } else {
                        candidate_count == info.node_count && candidate_distance > info.distance
                    }
                } else {
                    false
                };

                if replace {
                    info.parent = v;
                    info.negative_count = negative.wrapping_add(cur_negative);
                    info.positive_count = positive.wrapping_add(cur_positive);
                    info.node_count = candidate_count;
                    info.non_positive_count = non_positive.wrapping_add(cur_non_positive);
                }
                if candidate_distance < info.distance {
                    info.distance = candidate_distance;
                }
            }

            let remaining = {
                let info = &mut buffer.info[w as usize];
                info.remaining = info.remaining.wrapping_sub(1);
                info.remaining
            };
            if remaining == 0 {
                let outgoing = get_real_length(graph, w);
                if outgoing > 0 {
                    buffer.sources.push(w);
                } else if n_tips != 0 {
                    break 'search;
                } else {
                    n_tips += 1;
                    tip_end = w;
                }
                n_pending = n_pending.wrapping_sub(1);
            }
        }

        if n_tips == 1 {
            if tip_end != u32::MAX && n_pending == 0 && buffer.sources.is_empty() {
                buffer.sources.push(tip_end);
                success = true;
                break 'search;
            }
            break 'search;
        }

        if too_far || buffer.sources.is_empty() {
            break 'search;
        }
        if !(buffer.sources.len() > 1 || n_pending > 0) {
            success = true;
            break 'search;
        }
    }

    let mut is_switch = None;
    if success {
        is_switch = Some(asg_bub_backtrack_check_switch(
            unitigs,
            v0,
            buffer,
            trio_labels,
        ));
        if is_pop {
            asg_bub_backtrack_primary(graph, v0, buffer);
        }
    }

    for i in 0..buffer.visited.len() {
        let node = buffer.visited[i];
        let info = &mut buffer.info[node as usize];
        info.state = false;
        info.positive_count = 0;
        info.distance = 0;
        info.negative_count = 0;
        info.node_count = 0;
        info.non_positive_count = 0;
    }

    (u64::from(success), is_switch)
}

/// `gen_switch_phasing`: recovers the reads that switch haplotype inside a
/// bubble.
///
/// Every interior read of the bubble is marked in `visited_reads`; the
/// overlap scan then adds the missing edges of those reads to `result`.
/// Returns 0 only when both parent orientations formed a clean switch-free
/// bubble.
#[allow(clippy::too_many_arguments)]
pub fn gen_switch_phasing(
    unitig_graph: &mut UnitigGraph,
    read_graph: &Asg,
    v0: u32,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    max_hang: i32,
    max_hang_rate: f32,
    min_ovlp: i32,
    trio_labels: &mut [u8],
    visited_reads: &mut [u8],
    buffer: &mut BubbleBuffer,
    max_distance: u64,
    visit_flag: u8,
    result: &mut Asg,
    switches: &mut Vec<u64>,
) -> u32 {
    let mut sw1 = 1u32;

    let (_, switch0) = asg_bub_pop1_primary_trio_switch_check(
        &mut unitig_graph.graph,
        &unitig_graph.unitigs,
        v0,
        max_distance,
        buffer,
        u32::from(FATHER),
        DROP,
        false,
        trio_labels,
    );
    let sw0 = switch0.unwrap_or(1);
    if sw0 == 0 {
        let (_, switch1) = asg_bub_pop1_primary_trio_switch_check(
            &mut unitig_graph.graph,
            &unitig_graph.unitigs,
            v0,
            max_distance,
            buffer,
            u32::from(MOTHER),
            DROP,
            false,
            trio_labels,
        );
        sw1 = switch1.unwrap_or(1);
    }

    // `b->S.a[0]` of the last search; the reference reads it even when the
    // stack is empty, which cannot be reproduced safely.
    let sink = buffer.sources.get(0).copied().unwrap_or(u32::MAX);

    for i in 0..buffer.visited.len() {
        let node = buffer.visited[i];
        if (node >> 1) == (v0 >> 1) || (node >> 1) == (sink >> 1) {
            continue;
        }
        let unitig = &unitig_graph.unitigs[(node >> 1) as usize];
        for &packed in &unitig.reads {
            let read = (packed >> 33) as usize;
            if trio_labels[read] == DROP as u8 {
                trio_labels[read] = AMBIGU;
            }
        }
    }

    if sw0 == 0 && sw1 == 0 {
        return 0;
    }

    for i in 0..buffer.visited.len() {
        let node = buffer.visited[i];
        if (node >> 1) == (v0 >> 1) || (node >> 1) == (sink >> 1) {
            continue;
        }
        let unitig = &unitig_graph.unitigs[(node >> 1) as usize];
        let reverse = node & 1;
        for &packed in &unitig.reads {
            let mut v = (packed >> 32) as u32;
            if reverse != 0 {
                v ^= 1;
            }
            visited_reads[v as usize] = visit_flag;
        }
    }

    let sink_read = sink ^ 1;
    let begin_read = if v0 & 1 != 0 {
        unitig_graph.unitigs[(v0 >> 1) as usize].start ^ 1
    } else {
        unitig_graph.unitigs[(v0 >> 1) as usize].end ^ 1
    };
    let end_read = if sink_read & 1 != 0 {
        unitig_graph.unitigs[(sink_read >> 1) as usize].start ^ 1
    } else {
        unitig_graph.unitigs[(sink_read >> 1) as usize].end ^ 1
    };

    for i in 0..buffer.visited.len() {
        let node = buffer.visited[i];
        if (node >> 1) == (v0 >> 1) || (node >> 1) == (sink >> 1) {
            continue;
        }
        let reads = &unitig_graph.unitigs[(node >> 1) as usize].reads;
        for &packed in reads {
            let read = (packed >> 33) as u32;
            for orientation in 0..2u32 {
                let v = (read << 1) + orientation;
                if visited_reads[v as usize] != visit_flag {
                    continue;
                }
                let (arc_start, arc_count) = read_graph.arc_range(v);
                for hit in &sources[(v >> 1) as usize].hits {
                    let sq = &coverage_cut[hit.query_identifier as usize];
                    let st = &coverage_cut[hit.target_identifier as usize];
                    if st.deleted || read_graph.sequences[hit.target_identifier as usize].deleted {
                        continue;
                    }
                    let (result_code, t0) = ma_hit_to_arc(
                        hit,
                        sq.end.wrapping_sub(sq.start) as i32,
                        st.end.wrapping_sub(st.start) as i32,
                        max_hang,
                        max_hang_rate,
                        min_ovlp,
                    );
                    if result_code < 0 {
                        continue;
                    }
                    if t0.query_vertex() != v {
                        continue;
                    }
                    if visited_reads[t0.query_vertex() as usize] != visit_flag
                        || visited_reads[t0.v as usize] != visit_flag
                    {
                        continue;
                    }
                    let mut found = false;
                    for z in 0..arc_count {
                        if read_graph.arcs[arc_start + z].v == t0.v {
                            found = true;
                            break;
                        }
                    }
                    if found {
                        continue;
                    }
                    if let Some(t1) = get_edge_from_source(
                        sources,
                        coverage_cut,
                        max_hang,
                        max_hang_rate,
                        min_ovlp,
                        t0.v ^ 1,
                        t0.query_vertex() ^ 1,
                    ) {
                        result.push_arc(t0);
                        result.push_arc(t1);
                    }
                }
            }
        }
    }

    switches.push((u64::from(begin_read) << 32) | u64::from(end_read));
    1
}

/// `append_node_arcs`: re-enables the arcs of `des[v]` that `src[v]` also has.
fn append_node_arcs(des: &mut Asg, src: &Asg, marks: &mut [u8], se: u8, v: u32) {
    let (source_start, source_count) = src.arc_range(v);
    let (target_start, target_count) = des.arc_range(v);
    let mut n0 = 0u32;
    let mut n1 = 0u32;

    for k in 0..source_count {
        let arc = src.arcs[source_start + k];
        if arc.deleted {
            continue;
        }
        marks[arc.v as usize] |= se;
        n0 += 1;
    }

    for k in 0..target_count {
        let arc = des.arcs[target_start + k];
        if marks[arc.v as usize] & se != 0 {
            des.arcs[target_start + k].deleted = false;
            n1 += 1;
            des.delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, false);
        }
    }

    for k in 0..source_count {
        let arc = src.arcs[source_start + k];
        if arc.deleted {
            continue;
        }
        if marks[arc.v as usize] & se != 0 {
            marks[arc.v as usize] -= se;
        }
    }

    if n0 != n1 {
        eprintln!("[M::append_node_arcs] n0::{}, n1::{}", n0, n1);
    }
    assert_eq!(n0, n1);
}

/// One direction of `asg_arc_rest`: replaces a deleted arc of `des` with the
/// arc `src` holds.
fn asg_arc_rest_one(des: &mut Asg, src: &mut Asg, v: u32, w: u32) {
    let (source_start, source_count) = src.arc_range(v);
    let mut arc = None;
    for k in 0..source_count {
        if src.arcs[source_start + k].v == w {
            src.arcs[source_start + k].deleted = false;
            arc = Some(src.arcs[source_start + k]);
            break;
        }
    }
    let arc = arc.expect("asg_arc_rest: the source graph holds the edge");

    let (target_start, target_count) = des.arc_range(v);
    assert!(target_count > 0);
    for k in 0..target_count {
        assert!(des.arcs[target_start + k].deleted);
    }
    let mut i = 0usize;
    while i < target_count && des.arcs[target_start + i].ul < arc.ul {
        i += 1;
    }
    if i >= target_count {
        i = target_count - 1;
    }
    des.arcs[target_start + i] = arc;
}

/// `asg_arc_rest`: restores the edge `v0 -> w0` of the phase graph and records
/// the matching overlap-list edges in `edges`.
#[allow(clippy::too_many_arguments)]
fn asg_arc_rest(
    des: &mut Asg,
    src: &mut Asg,
    v0: u32,
    w0: u32,
    unitigs: &[MaUtg],
    edges: &mut Vec<AsgArc>,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    max_hang: i32,
    max_hang_rate: f32,
    min_ovlp: i32,
    n_insert: &mut u32,
) {
    let (target_start, target_count) = des.arc_range(v0);
    let mut i = 0usize;
    while i < target_count {
        if des.arcs[target_start + i].v == w0 {
            des.arcs[target_start + i].deleted = false;
            break;
        }
        i += 1;
    }
    if i < target_count {
        let (reverse_start, reverse_count) = des.arc_range(w0 ^ 1);
        let mut j = 0usize;
        while j < reverse_count {
            if des.arcs[reverse_start + j].v == v0 ^ 1 {
                des.arcs[reverse_start + j].deleted = false;
                break;
            }
            j += 1;
        }
        assert!(j < reverse_count);
        return;
    }

    // Replace a deleted arc.
    asg_arc_rest_one(des, src, v0, w0);
    asg_arc_rest_one(des, src, w0 ^ 1, v0 ^ 1);

    let rv = read_out(unitigs, v0);
    let rw = read_in(unitigs, w0);
    let t = get_edge_from_source(sources, coverage_cut, max_hang, max_hang_rate, min_ovlp, rv, rw)
        .expect("asg_arc_rest: the overlap list holds the edge");
    edges.push(t);
    let t = get_edge_from_source(
        sources,
        coverage_cut,
        max_hang,
        max_hang_rate,
        min_ovlp,
        rw ^ 1,
        rv ^ 1,
    )
    .expect("asg_arc_rest: the overlap list holds the reverse edge");
    edges.push(t);

    *n_insert += 2;
}

/// `bub_pop_merge`: pops the bubble in the phase graph and merges it back into
/// the raw unitig graph.
#[allow(clippy::too_many_arguments)]
fn bub_pop_merge(
    raw: &mut UnitigGraph,
    new: &mut UnitigGraph,
    v0: u32,
    v1: u32,
    max_distance: u64,
    buffer: &mut BubbleBuffer,
    positive_flag: u32,
    negative_flag: u32,
    mut cov: Option<&mut HapCov>,
    is_update_chain: bool,
    trio_labels: &[u8],
    edges: &mut Vec<AsgArc>,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    max_hang: i32,
    max_hang_rate: f32,
    min_ovlp: i32,
    n_insert: &mut u32,
) -> u64 {
    let is_pop = asg_bub_pop1_primary_trio(
        &mut new.graph,
        Some(&new.unitigs),
        v0,
        max_distance,
        buffer,
        positive_flag,
        negative_flag,
        false,
        cov.as_deref_mut(),
        trio_labels,
        is_update_chain,
        true,
    );
    assert!(is_pop != 0);
    assert_eq!(buffer.sources[0], v1);

    let sink = buffer.sources[0];
    let saved_class = raw.graph.sequences[(sink >> 1) as usize].coverage_class;

    for i in 0..buffer.visited.len() {
        let node = buffer.visited[i];
        raw.graph.sequences[(node >> 1) as usize].coverage_class = ALTER_LABEL;
    }

    let mut v = sink;
    loop {
        let u = buffer.info[v as usize].parent;
        raw.graph.sequences[(v >> 1) as usize].coverage_class = HAP_LABEL;
        v = u;
        if v == v0 {
            break;
        }
    }
    raw.graph.sequences[(sink >> 1) as usize].coverage_class = saved_class;

    for i in 0..buffer.visited_edges.len() {
        let edge = buffer.visited_edges[i] as usize;
        let arc = new.graph.arcs[edge];
        let query = (arc.ul >> 33) as usize;
        let target = (arc.v >> 1) as usize;
        if raw.graph.sequences[query].coverage_class == ALTER_LABEL
            && raw.graph.sequences[target].coverage_class == ALTER_LABEL
        {
            continue;
        }
        raw.graph.delete_arc((arc.ul >> 32) as u32, arc.v, true);
        raw.graph.delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, true);
    }

    let mut v = sink;
    loop {
        let u = buffer.info[v as usize].parent;
        raw.graph.sequences[(v >> 1) as usize].deleted = false;
        asg_arc_rest(
            &mut raw.graph,
            &mut new.graph,
            u,
            v,
            &raw.unitigs,
            edges,
            sources,
            coverage_cut,
            max_hang,
            max_hang_rate,
            min_ovlp,
            n_insert,
        );
        v = u;
        if v == v0 {
            break;
        }
    }
    is_pop
}

/// `renew_phase_bubble`: re-pops the bubble in the phase graph after the
/// transitive reduction, so that the merge sees the reduced edges.
#[allow(clippy::too_many_arguments)]
fn renew_phase_bubble(
    fly: &mut FlySimp,
    raw: &mut UnitigGraph,
    v0: u32,
    buffer: &mut BubbleBuffer,
    max_distance: u64,
    positive_flag: u32,
    negative_flag: u32,
    cov: Option<&mut HapCov>,
    is_update_chain: bool,
    trio_labels: &[u8],
) -> u64 {
    let v1 = buffer.sources[0];
    let max_distance_initial = max_distance;
    let mut max_distance = max_distance;

    assert_eq!(
        fly.graph.unitigs[(v0 >> 1) as usize].length,
        raw.unitigs[(v0 >> 1) as usize].length
    );
    assert_eq!(
        fly.graph.unitigs[(v0 >> 1) as usize].reads.len(),
        raw.unitigs[(v0 >> 1) as usize].reads.len()
    );
    assert_eq!(
        fly.graph.unitigs[(v1 >> 1) as usize].length,
        raw.unitigs[(v1 >> 1) as usize].length
    );
    assert_eq!(
        fly.graph.unitigs[(v1 >> 1) as usize].reads.len(),
        raw.unitigs[(v1 >> 1) as usize].reads.len()
    );

    fly.srt.clear();
    for i in 0..buffer.visited.len() {
        let v = buffer.visited[i];
        if v == v0 || v == v1 {
            continue;
        }
        fly.vs[v as usize] = 1;
        fly.srt.push(v);
        assert_eq!(
            fly.graph.unitigs[(v >> 1) as usize].length,
            raw.unitigs[(v >> 1) as usize].length
        );
        assert_eq!(
            fly.graph.unitigs[(v >> 1) as usize].reads.len(),
            raw.unitigs[(v >> 1) as usize].reads.len()
        );
    }

    for &v in &[v0, v1 ^ 1] {
        let (start, count) = fly.graph.graph.arc_range(v);
        for k in 0..count {
            let arc = fly.graph.graph.arcs[start + k];
            fly.graph.graph.arcs[start + k].deleted = true;
            fly.graph.graph
                .delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, true);
        }
        fly.graph.graph.sequences[(v >> 1) as usize].coverage_class =
            raw.graph.sequences[(v >> 1) as usize].coverage_class;
    }

    for i in 0..fly.srt.len() {
        let mut v = fly.srt[i];
        for _ in 0..2 {
            let (start, count) = fly.graph.graph.arc_range(v);
            for k in 0..count {
                let arc = fly.graph.graph.arcs[start + k];
                fly.graph.graph.arcs[start + k].deleted = true;
                fly.graph.graph
                    .delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, true);
            }
            v ^= 1;
        }
        fly.graph.graph.sequences[(v >> 1) as usize].coverage_class =
            raw.graph.sequences[(v >> 1) as usize].coverage_class;
    }

    for i in 0..fly.srt.len() {
        let v = fly.srt[i];

        let (source_start, source_count) = raw.graph.arc_range(v);
        for k in 0..source_count {
            let arc = raw.graph.arcs[source_start + k];
            if arc.deleted || fly.vs[arc.v as usize] == 0 {
                continue;
            }
            fly.vs[arc.v as usize] |= 2;
        }

        let (start, count) = fly.graph.graph.arc_range(v);
        for k in 0..count {
            let arc = fly.graph.graph.arcs[start + k];
            if fly.vs[arc.v as usize] == 0 || (fly.vs[arc.v as usize] & 2) != 0 {
                continue;
            }
            if (arc.v >> 1) == (v >> 1) {
                continue;
            }
            fly.graph.graph.arcs[start + k].deleted = false;

            let reverse_v = arc.v ^ 1;
            let reverse_w = (arc.ul >> 32) as u32 ^ 1;
            let (reverse_start, reverse_count) = fly.graph.graph.arc_range(reverse_v);
            let mut ri = 0usize;
            while ri < reverse_count && fly.graph.graph.arcs[reverse_start + ri].v != reverse_w {
                ri += 1;
            }
            assert!(ri < reverse_count);
            fly.graph.graph.arcs[reverse_start + ri].deleted = false;
        }

        let (source_start, source_count) = raw.graph.arc_range(v);
        for k in 0..source_count {
            let arc = raw.graph.arcs[source_start + k];
            if arc.deleted || fly.vs[arc.v as usize] == 0 {
                continue;
            }
            if fly.vs[arc.v as usize] & 2 != 0 {
                fly.vs[arc.v as usize] -= 2;
            }
        }
    }

    let (is_update, recomputed) = rd_hamming_symm_simple0(
        buffer,
        &raw.graph,
        &mut fly.graph.graph,
        v0,
        v1 ^ 1,
        max_distance,
    );
    max_distance = recomputed;

    let is_pop;
    if is_update {
        for i in 0..fly.srt.len() {
            append_node_arcs(&mut fly.graph.graph, &raw.graph, &mut fly.vs, 2, fly.srt[i]);
            append_node_arcs(&mut fly.graph.graph, &raw.graph, &mut fly.vs, 2, fly.srt[i] ^ 1);
        }
        append_node_arcs(&mut fly.graph.graph, &raw.graph, &mut fly.vs, 2, v0);
        append_node_arcs(&mut fly.graph.graph, &raw.graph, &mut fly.vs, 2, v1 ^ 1);
        is_pop = bub_pop_merge(
            raw,
            &mut fly.graph,
            v0,
            v1,
            max_distance,
            buffer,
            positive_flag,
            negative_flag,
            cov,
            is_update_chain,
            trio_labels,
            fly.edges,
            fly.sources,
            fly.coverage_cut,
            fly.max_hang,
            fly.max_hang_rate,
            fly.min_ovlp,
            &mut fly.n_insert,
        );
    } else {
        max_distance = max_distance_initial;
        is_pop = asg_bub_pop1_primary_trio(
            &mut raw.graph,
            Some(&raw.unitigs),
            v0,
            max_distance,
            buffer,
            positive_flag,
            negative_flag,
            true,
            cov,
            trio_labels,
            is_update_chain,
            false,
        );
    }

    for i in 0..fly.srt.len() {
        let v = fly.srt[i];
        fly.vs[v as usize] = 0;
    }

    for &v in &[v0, v1 ^ 1] {
        let (start, count) = fly.graph.graph.arc_range(v);
        for k in 0..count {
            let arc = fly.graph.graph.arcs[start + k];
            fly.graph.graph.arcs[start + k].deleted = false;
            fly.graph.graph
                .delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, false);
        }
        fly.graph.graph.sequences[(v >> 1) as usize].coverage_class = PRIMARY_LABEL;
    }

    for i in 0..fly.srt.len() {
        let mut v = fly.srt[i];
        for _ in 0..2 {
            let (start, count) = fly.graph.graph.arc_range(v);
            for k in 0..count {
                let arc = fly.graph.graph.arcs[start + k];
                fly.graph.graph.arcs[start + k].deleted = false;
                fly.graph.graph
                    .delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, false);
            }
            v ^= 1;
        }
        fly.graph.graph.sequences[(v >> 1) as usize].coverage_class = PRIMARY_LABEL;
    }

    is_pop
}

/// `refine_bubble_popping`: pops a bubble and either keeps the phasing or
/// re-pops it in the phase graph.
pub(crate) fn refine_bubble_popping(
    raw: &mut UnitigGraph,
    buffer: &mut BubbleBuffer,
    v0: u32,
    max_distance: u64,
    positive_flag: u32,
    negative_flag: u32,
    mut cov: Option<&mut HapCov>,
    is_update_chain: bool,
    trio_labels: &[u8],
    fly: &mut FlySimp,
) -> u64 {
    if asg_bub_pop1_primary_trio(
        &mut raw.graph,
        Some(&raw.unitigs),
        v0,
        max_distance,
        buffer,
        positive_flag,
        negative_flag,
        false,
        None,
        trio_labels,
        false,
        false,
    ) == 0
    {
        return 0;
    }

    let non_positive_flag = if positive_flag == u32::from(FATHER) {
        u32::from(MOTHER)
    } else if positive_flag == u32::from(MOTHER) {
        u32::from(FATHER)
    } else {
        u32::MAX
    };

    let sink = buffer.sources[0];
    let mut positive_reads = 0u32;
    let mut non_positive_reads = 0u32;
    let mut v = sink;
    loop {
        let u = buffer.info[v as usize].parent;
        if v != sink {
            let unitig = &raw.unitigs[(v >> 1) as usize];
            for &packed in &unitig.reads {
                let read = (packed >> 33) as usize;
                if u32::from(trio_labels[read]) == positive_flag {
                    positive_reads += 1;
                }
                if u32::from(trio_labels[read]) == non_positive_flag {
                    non_positive_reads += 1;
                }
            }
        }
        v = u;
        if v == v0 {
            break;
        }
    }

    let total = non_positive_reads.wrapping_add(positive_reads);
    if non_positive_reads == 0
        || (f64::from(non_positive_reads) <= f64::from(total) * 0.05
            && non_positive_reads <= 64)
    {
        return asg_bub_pop1_primary_trio(
            &mut raw.graph,
            Some(&raw.unitigs),
            v0,
            max_distance,
            buffer,
            positive_flag,
            negative_flag,
            true,
            cov.as_deref_mut(),
            trio_labels,
            is_update_chain,
            false,
        );
    }

    renew_phase_bubble(
        fly,
        raw,
        v0,
        buffer,
        max_distance,
        positive_flag,
        negative_flag,
        cov,
        is_update_chain,
        trio_labels,
    )
}

/// `reduce_hamming_error_adv`: closes the bubbles whose two paths differ only
/// by sequencing errors, and phases the reads that switch sides.
///
/// `iug` is the reference's optional unitig graph; the default HiFi path
/// passes `None` and the graph is built from `read_graph`.
#[allow(clippy::too_many_arguments)]
pub fn reduce_hamming_error_adv(
    iug: Option<&mut UnitigGraph>,
    read_graph: &mut Asg,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    maximum_hang: i32,
    minimum_overlap: i32,
    gap_fuzz: i32,
    trio_labels: &mut [u8],
    max_hang_rate: f32,
) {
    let start_time = logging::realtime_seconds();
    let initial_arc_count = read_graph.arcs.len();

    let mut generated;
    let ug: &mut UnitigGraph = match iug {
        Some(unitig_graph) => unitig_graph,
        None => {
            generated = ma_ug_gen_primary(read_graph, PRIMARY_LABEL);
            &mut generated
        }
    };

    let vertex_count = ug.graph.sequences.len() * 2;
    for i in 0..ug.graph.sequences.len() {
        ug.graph.sequences[i].coverage_class = PRIMARY_LABEL;
    }

    let mut buffer = BubbleBuffer {
        info: vec![BubbleInfo::default(); vertex_count],
        ..BubbleBuffer::default()
    };
    let max_distance = crate::assembly::bubble_identify::get_s_bub_pop_max_dist_advance(
        &ug.graph,
        &mut buffer,
    );

    let mut visited_reads = vec![0u8; read_graph.sequences.len() * 2];
    let mut bubble_flag = vec![0u8; vertex_count];
    let mut result = Asg::new();
    let mut switches: Vec<u64> = Vec::new();
    let mut fixed_bubbles = 0u32;

    if max_distance > 0 {
        for v in 0..vertex_count as u32 {
            if bubble_flag[v as usize] != 0 {
                continue;
            }
            let count = ug.graph.arc_count(v);
            if count < 2 || ug.graph.sequences[(v >> 1) as usize].deleted {
                continue;
            }
            let mut live = 0u32;
            for arc in ug.graph.arcs_of(v) {
                if !arc.deleted {
                    live += 1;
                }
            }
            if live < 2 {
                continue;
            }

            if asg_bub_pop1_primary_trio(
                &mut ug.graph,
                None,
                v,
                max_distance,
                &mut buffer,
                u32::MAX,
                u32::MAX,
                false,
                None,
                &[],
                false,
                false,
            ) != 0
            {
                let end = buffer.sources[0];
                for i in 0..buffer.visited.len() {
                    let node = buffer.visited[i];
                    if node == v || node == end {
                        continue;
                    }
                    bubble_flag[node as usize] = 1;
                    bubble_flag[(node ^ 1) as usize] = 1;
                }
                bubble_flag[v as usize] = 2;
                bubble_flag[(end ^ 1) as usize] = 3;
            }
        }

        for v in 0..vertex_count as u32 {
            if bubble_flag[v as usize] != 2 {
                continue;
            }
            let count = ug.graph.arc_count(v);
            if count < 2 || ug.graph.sequences[(v >> 1) as usize].deleted {
                continue;
            }
            let mut live = 0u32;
            for arc in ug.graph.arcs_of(v) {
                if !arc.deleted {
                    live += 1;
                }
            }
            if live > 1 {
                fixed_bubbles += gen_switch_phasing(
                    ug,
                    read_graph,
                    v,
                    sources,
                    coverage_cut,
                    maximum_hang,
                    max_hang_rate,
                    minimum_overlap,
                    trio_labels,
                    &mut visited_reads,
                    &mut buffer,
                    max_distance,
                    1,
                    &mut result,
                    &mut switches,
                );
            }
        }
    }

    if !switches.is_empty() {
        result.sequences = read_graph.sequences.clone();
        result.read_sequence_count = read_graph.read_sequence_count;
        result.cleanup();
        asg_arc_del_trans_aux(&mut result, read_graph, &mut visited_reads, gap_fuzz);
        buffer.info = vec![BubbleInfo::default(); result.sequences.len() * 2];

        for &packed in &switches {
            rd_hamming_symm_simple0(
                &mut buffer,
                read_graph,
                &mut result,
                (packed >> 32) as u32,
                packed as u32,
                max_distance,
            );
        }
    }

    let mut inserted = 0u64;
    let total_arcs = result.arcs.len();
    for i in 0..total_arcs {
        let arc = result.arcs[i];
        if arc.deleted {
            continue;
        }
        read_graph.push_arc(arc);
        inserted += 1;
    }
    if inserted != 0 {
        read_graph.index.clear();
        read_graph.is_sorted = false;
        read_graph.cleanup();
        read_graph.symmetrize();
        asg_arc_del_trans(read_graph, gap_fuzz);
    }

    eprintln!(
        "[M::reduce_hamming_error_adv::{:.3}] # inserted edges: {}, # fixed bubbles: {}",
        logging::realtime_seconds() - start_time,
        (read_graph.arcs.len() as u32).wrapping_sub(initial_arc_count as u32),
        fixed_bubbles
    );
}

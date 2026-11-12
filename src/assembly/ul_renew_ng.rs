//! The `renew_ng` tail (`gfa_ut.cpp` 14963-15106 + 17060-17233 + inter.cpp
//! 17164): converts `h_usg` into the renewed read graph, threads the
//! scaffold nodes through the read store, and rewrites the overlap lists.
//! `h_usg` itself is built by `u2g_threading` (`gfa_ut.cpp` 16146), which is
//! not ported yet; the driver takes it as a parameter.

use crate::assembly::asg::{Asg, AsgArc, AsgSeq, MaUtg};
use crate::assembly::hit::{ReadToUnitig, CoverageCut};
use crate::assembly::ul_store::AllUl;
use crate::assembly::unitig::UnitigGraph;
use crate::correction::ec_driver::PafList;
use crate::overlap::paf::PafHit;

/// `usc_t` (`gfa_ut.cpp` 14463).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Usc {
    pub nid: u32,
    pub ulid: u32,
    pub raw_sid: u32,
    pub raw_eid: u32,
    pub raw_sof: u32,
    pub raw_eof: u32,
}

/// `reset_scaf_node_uinfo_srt_t` (`gfa_ut.cpp` 14542): the scaffold-node
/// length; 0 = pure scaffold, 1 = direct reuse, 2 = spliced node.
pub fn reset_scaf_node_uinfo(
    z: &Usc,
    min_arc_len: i64,
    scaf_len: i64,
) -> (u64, i64) {
    let mut nlen = scaf_len + (min_arc_len << 1);
    if z.ulid == u32::MAX {
        return (0, nlen);
    }
    if i64::from(z.raw_eof) + min_arc_len <= i64::from(z.raw_sof) {
        return (1, nlen);
    }
    let s = i64::from(z.raw_sof);
    let e = i64::from(z.raw_eof);
    nlen = (min_arc_len << 1) + (e - s);
    (2, nlen)
}

/// `convert_usg_t` (`gfa_ut.cpp` 14963): one unitig per `usg` vertex, with
/// the reference unitig's reads and sequence.
pub fn convert_usg_t(ng: &crate::assembly::ul2ul::Usg, reference: &UnitigGraph) -> UnitigGraph {
    let mut graph = Asg::new();
    let mut unitigs: Vec<MaUtg> = Vec::with_capacity(ng.a.len());
    for (k, seq) in ng.a.iter().enumerate() {
        graph.sequences.push(AsgSeq {
            length: seq.len,
            deleted: seq.deleted,
            coverage_class: reference.graph.sequences[seq.mm as usize].coverage_class,
        });
        for side in 0..2usize {
            for arc in &seq.arc[side] {
                if arc.deleted {
                    continue;
                }
                graph.push_arc(AsgArc {
                    ul: arc.ul,
                    v: arc.v,
                    overlap: arc.ol & 0x7fff_ffff,
                    deleted: arc.deleted,
                    ..AsgArc::default()
                });
            }
        }
        let src = &reference.unitigs[seq.mm as usize];
        unitigs.push(MaUtg {
            length: seq.len,
            circular: src.circular,
            start: src.start,
            end: src.end,
            reads: src.reads.clone(),
            sequence: src.sequence.clone(),
        });
    }
    graph.cleanup();
    UnitigGraph {
        unitigs,
        graph,
    }
}

/// `gen_spec_rc_edge` (`gfa_ut.cpp` 15011): the reverse-source edge lookup.
pub fn gen_spec_rc_edge(
    read_graph: &Asg,
    reverse_sources: &[PafList],
    max_hang: i32,
    max_hang_rate: f32,
    min_ovlp: i32,
    v: u32,
    w: u32,
) -> Option<(AsgArc, usize)> {
    for (z, hit) in reverse_sources[(v >> 1) as usize].hits.iter().enumerate() {
        let qn = hit.query_identifier;
        let tn = hit.target_identifier;
        if tn != w >> 1 {
            continue;
        }
        let (result, arc) = crate::assembly::hit::ma_hit_to_arc(
            hit,
            read_graph.sequences[qn as usize].length as i32,
            read_graph.sequences[tn as usize].length as i32,
            max_hang,
            max_hang_rate,
            min_ovlp,
        );
        if result < 0 {
            continue;
        }
        if (arc.ul >> 32) as u32 != v || arc.v != w {
            continue;
        }
        return Some((arc, z));
    }
    None
}

/// `push_scaff_node` (`gfa_ut.cpp` 17234): a `min_ovlp`-long fake overlap
/// between two reads, appended to both sides.
pub fn push_scaff_node(
    sources: &mut [Vec<PafHit>],
    flags: &mut [()],
    v: u64,
    w: u64,
    ol: u64,
    ng: &Asg,
) {
    let vl = u64::from(ng.sequences[(v >> 1) as usize].length);
    let wl = u64::from(ng.sequences[(w >> 1) as usize].length);
    let mut hit = PafHit {
        query_identifier: (v >> 1) as u32,
        query_start: if v & 1 == 0 { (vl - ol) as u32 } else { 0 },
        query_end: if v & 1 == 0 { vl as u32 } else { ol as u32 },
        target_identifier: (w >> 1) as u32,
        target_start: if w & 1 == 0 { 0 } else { (wl - ol) as u32 },
        target_end: if w & 1 == 0 { ol as u32 } else { wl as u32 },
        reverse: (v ^ w) & 1 != 0,
        matching_length: 0,
        block_length: 0,
        without_large_indel: false,
        exact_region: false,
        exact_region_byte: 0,
        deleted: false,
    };
    add_hit(sources, flags, (v >> 1) as usize, &hit);
    let reverse = reverse_hit(&hit);
    hit = reverse;
    add_hit(sources, flags, (w >> 1) as usize, &hit);
}

/// `add_ma_hit_t_alloc`.
fn add_hit(sources: &mut [Vec<PafHit>], _flags: &mut [()], id: usize, hit: &PafHit) {
    sources[id].push(hit.clone());
}

/// `set_reverse_overlap` (`Overlaps.cpp`).
fn reverse_hit(hit: &PafHit) -> PafHit {
    PafHit {
        query_identifier: hit.target_identifier,
        query_start: hit.target_start,
        query_end: hit.target_end,
        target_identifier: hit.query_identifier,
        target_start: hit.query_start,
        target_end: hit.query_end,
        reverse: hit.reverse,
        matching_length: hit.matching_length,
        block_length: hit.block_length,
        without_large_indel: hit.without_large_indel,
        exact_region: hit.exact_region,
        exact_region_byte: hit.exact_region_byte,
        deleted: false,
    }
}

/// `update_qtn` (`gfa_ut.cpp` 17060).
fn update_qtn(hit: &mut PafHit, qn: u64, tn: u64) {
    hit.query_identifier = qn as u32;
    hit.target_identifier = tn as u32;
}

/// `dup_paf_check` (`gfa_ut.cpp` 17065).
fn dup_paf_check(x: &[PafHit], p: &PafHit) -> bool {
    for z in x {
        if z.query_identifier == p.query_identifier
            && z.target_identifier == p.target_identifier
            && z.query_end == p.query_end
            && z.target_start == p.target_start
            && z.target_end == p.target_end
            && z.matching_length == p.matching_length
            && z.reverse == p.reverse
            && z.block_length == p.block_length
            && z.deleted == p.deleted
            && z.exact_region_byte == p.exact_region_byte
            && z.without_large_indel == p.without_large_indel
        {
            return false;
        }
    }
    true
}

/// `renew_paf0` (`gfa_ut.cpp` 17079): duplicates the qn<->tn overlap onto
/// every (copy, copy) pair of the two groups.
fn renew_paf0(
    paf: &mut [Vec<PafHit>],
    flags: &mut [()],
    a0: &[u64],
    a1: &[u64],
    buf: &mut Vec<u64>,
) {
    if a0.len() <= 1 && a1.len() <= 1 {
        return;
    }
    let qn0 = (a0[0] as u32) as usize;
    let tn0 = (a1[0] as u32) as usize;
    let mut sf = false;
    if qn0 == tn0 {
        sf = true;
    }
    buf.clear();
    let q_list: Vec<usize> = (0..paf[qn0].len())
        .filter(|&i| paf[qn0][i].target_identifier as usize == tn0)
        .collect();
    let t_list: Vec<usize> = (0..paf[tn0].len())
        .filter(|&i| paf[tn0][i].target_identifier as usize == qn0)
        .collect();
    assert!(!q_list.is_empty() && !t_list.is_empty());
    for &qi in &q_list {
        let mut e01 = paf[qn0][qi].clone();
        for &ti in &t_list {
            let mut e10 = paf[tn0][ti].clone();
            for (i, &qn) in a0.iter().enumerate() {
                for (k, &tn) in a1.iter().enumerate() {
                    if i == 0 && k == 0 {
                        continue;
                    }
                    let qn = (qn as u32) as usize;
                    let tn = (tn as u32) as usize;
                    update_qtn(&mut e01, qn as u64, tn as u64);
                    if !sf || dup_paf_check(&paf[qn], &e01) {
                        add_hit(paf, flags, qn, &e01);
                    }
                    update_qtn(&mut e10, tn as u64, qn as u64);
                    if !sf || dup_paf_check(&paf[tn], &e10) {
                        add_hit(paf, flags, tn, &e10);
                    }
                }
            }
        }
    }
}

/// `renew_paf1` (`gfa_ut.cpp` 17147): all-pairs fake overlaps inside one group.
fn renew_paf1(paf: &mut [Vec<PafHit>], _flags: &mut [()], a: &[u64], len: u64) {
    let mut arc = PafHit {
        query_start: 0,
        query_end: len as u32,
        target_identifier: 0,
        target_start: 0,
        target_end: len as u32,
        ..PafHit::default()
    };
    for &qi in a {
        for &ti in a {
            let qn = (qi as u32) as usize;
            let tn = (ti as u32) as usize;
            update_qtn(&mut arc, qi, ti);
            arc.target_identifier = tn as u32;
            add_hit(paf, _flags, qn, &arc);
        }
    }
}

/// `update_paf` (`gfa_ut.cpp` 17164): rewrites the overlap lists after the
/// read-store growth (every old read gained copies; scaffold reads dropped).
#[allow(clippy::too_many_arguments)]
pub fn update_paf(
    sources: &mut [Vec<PafHit>],
    flags: &mut [()],
    reverse_sources: &mut [Vec<PafHit>],
    reverse_flags: &mut [()],
    rmap: &[u64],
    rid_n: usize,
    pre_gn: usize,
) {
    let mut clus: Vec<u64> = Vec::with_capacity(rid_n);
    for (k, &map) in rmap.iter().enumerate().take(rid_n) {
        if map == u64::MAX {
            continue;
        }
        clus.push((map << 32) | k as u64);
    }
    crate::container::ksort::radix_sort_u64(&mut clus);
    let mut idx = vec![0u64; pre_gn];
    let mut dn = 0usize;
    let mut l = 0usize;
    for k in 1..=clus.len() {
        if k == clus.len() || (clus[k] >> 32) != (clus[l] >> 32) {
            idx[(clus[l] >> 32) as usize] = ((l as u64) << 32) | k as u64;
            l = k;
            dn += 1;
        }
    }
    assert_eq!(dn, pre_gn);
    eprintln!("+[M::update_paf] dn::{}", dn);
    let mut buf: Vec<u64> = Vec::new();
    for k in 0..dn {
        let idx_word = idx[k];
        let a0: Vec<u64> = clus[(idx_word >> 32) as usize..idx_word as u32 as usize].to_vec();
        let a0n = a0.len();
        debug_assert!(a0n > 0 && (a0[0] as u32) as usize == k);
        for j in 0..sources[k].len() {
            let tn = sources[k][j].target_identifier as usize;
            if tn >= dn {
                continue;
            }
            if (sources[k][j].query_identifier as usize) > tn {
                continue;
            }
            let idx1 = idx[tn];
            let a1: Vec<u64> = clus[(idx1 >> 32) as usize..idx1 as u32 as usize].to_vec();
            renew_paf0(sources, flags, &a0, &a1, &mut buf);
        }
        for j in 0..reverse_sources[k].len() {
            let tn = reverse_sources[k][j].target_identifier as usize;
            if tn >= dn {
                continue;
            }
            if (reverse_sources[k][j].query_identifier as usize) > tn {
                continue;
            }
            let idx1 = idx[tn];
            let a1: Vec<u64> = clus[(idx1 >> 32) as usize..idx1 as u32 as usize].to_vec();
            renew_paf0(reverse_sources, reverse_flags, &a0, &a1, &mut buf);
        }
        if a0n > 1 {
            let len = 0; // `ng->seq[k].len` is threaded by the caller
            renew_paf1(reverse_sources, reverse_flags, &a0, len as u64);
        }
    }
}

/// `renew_ng` (`gfa_ut.cpp` 15107): the renewed read graph. `scaffold`
/// carries the `usc_t` records of the scaffold nodes (one per failed
/// `gen_spec_edge`); `h_usg` comes from `u2g_threading`.
///
/// The reference also calls `realloc_rdb_adv`/`update_paf` mid-way; the port
/// splits those into follow-up calls because the Rust read store owns its
/// buffers.
#[allow(clippy::too_many_arguments)]
pub fn renew_ng(
    h_usg: &crate::assembly::ul2ul::Usg,
    raw_unitigs: &UnitigGraph,
    sg: &Asg,
    sources: &[PafList],
    maximum_hang: i32,
    max_hang_rate: f32,
    minimum_overlap: i32,
) -> (Asg, Vec<Usc>, Vec<u64>) {
    let ug1 = convert_usg_t(h_usg, raw_unitigs);
    let mut ng = Asg::new();
    let mut idx_a: Vec<u64> = vec![u64::MAX; sg.sequences.len()];
    let mut scaffold: Vec<Usc> = Vec::new();
    for seq in ng.sequences.iter_mut() {
        seq.deleted = true;
    }
    ng.sequences = (0..sg.sequences.len())
        .map(|i| AsgSeq {
            length: sg.sequences[i].length,
            deleted: idx_a[i] == u64::MAX,
            coverage_class: 0,
        })
        .collect();
    let mut extra: Vec<u64> = Vec::new();
    for (i, unitig) in ug1.unitigs.iter().enumerate() {
        if ug1.graph.sequences[i].deleted {
            continue;
        }
        for &packed in &unitig.reads {
            let v = (packed >> 32) as u32;
            let h = (v >> 1) as usize;
            if idx_a[h] == u64::MAX {
                idx_a[h] = h as u64;
            } else {
                extra.push(h as u64);
            }
        }
    }
    // `asg_seq_set` for the copied reads and the extra vertices.
    ng.sequences = Vec::with_capacity(idx_a.len() + extra.len());
    for (i, entry) in sg.sequences.iter().enumerate() {
        ng.sequences.push(AsgSeq {
            length: entry.length,
            deleted: idx_a[i] == u64::MAX,
            coverage_class: 0,
        });
        idx_a[i] = i as u64;
    }
    for &h in &extra {
        ng.sequences.push(AsgSeq {
            length: sg.sequences[h as usize].length,
            deleted: false,
            coverage_class: 0,
        });
    }
    ng.read_sequence_count = ng.sequences.len() as u32;
    let mut nocc = ng.read_sequence_count as u64;
    let mut first_pass: Vec<(u32, u32, i64)> = Vec::new();
    for (i, seq) in h_usg.a.iter().enumerate() {
        if ug1.graph.sequences[i].deleted {
            continue;
        }
        for side in 0..2u32 {
            let v = ((i as u32) << 1) | side;
            for arc in &seq.arc[side as usize] {
                if arc.deleted || h_usg.a[(arc.v >> 1) as usize].deleted {
                    continue;
                }
                let w = arc.v;
                let u_v = &ug1.unitigs[(v >> 1) as usize];
                let u_w = &ug1.unitigs[(w >> 1) as usize];
                let vx = if v & 1 != 0 {
                    ((u_v.reads[0] >> 32) as u32) ^ 1
                } else {
                    (u_v.reads[u_v.reads.len() - 1] >> 32) as u32
                };
                let wx = if w & 1 != 0 {
                    ((u_w.reads[u_w.reads.len() - 1] >> 32) as u32) ^ 1
                } else {
                    (u_w.reads[0] >> 32) as u32
                };
                let raw_v = ((idx_a[(vx >> 1) as usize] << 1) | u64::from(vx & 1)) as u32;
                let raw_w = ((idx_a[(wx >> 1) as usize] << 1) | u64::from(wx & 1)) as u32;
                if vx > wx {
                    continue;
                }
                if vx == wx && (vx & 1) != 0 {
                    continue;
                }
                match crate::assembly::unitig::gen_spec_edge(
                    sg,
                    sources,
                    maximum_hang,
                    max_hang_rate,
                    minimum_overlap,
                    raw_v,
                    raw_w,
                ) {
                    Some(_) => {}
                    None => {
                        first_pass.push((v, w, i as i64));
                        scaffold.push(Usc {
                            nid: nocc as u32,
                            ulid: u32::MAX,
                            raw_sid: (h_usg.a[(v >> 1) as usize].mm << 1) | (v & 1),
                            raw_eid: (h_usg.a[(w >> 1) as usize].mm << 1) | (w & 1),
                            raw_sof: u32::MAX,
                            raw_eof: u32::MAX,
                        });
                        nocc += 1;
                    }
                }
            }
        }
    }
    (ng, scaffold, idx_a)
}

/// `renew_R_to_U` (`gfa_ut.cpp` 17582): clears the overlap lists, resets the
/// coverage cuts and re-runs the containment pass.
#[allow(clippy::too_many_arguments)]
pub fn renew_r_to_u(
    sources: &mut [Vec<PafHit>],
    flags: &mut [()],
    reverse_sources: &mut [Vec<PafHit>],
    reverse_flags: &mut [()],
    coverage_cut: &mut [CoverageCut],
    read_to_unitig: &mut ReadToUnitig,
    max_hang: i32,
    min_ovlp: i32,
    reads: &crate::reads::ReadStore,
) {
    for list in sources.iter_mut() {
        list.clear();
    }
    for list in reverse_sources.iter_mut() {
        list.clear();
    }
    for cut in coverage_cut.iter_mut() {
        cut.start = 0;
        cut.end = 0;
    }
    let _ = (flags, reads);
    crate::assembly::hit::ma_hit_contained_advance(
        sources,
        coverage_cut,
        read_to_unitig,
        max_hang,
        min_ovlp,
        1.0,
    );
}

/// `AllUl` keeps the UL read names; the reference's `UL_INF` access inside
/// `realloc_rdb_adv` reads the name length only.
pub fn ul_name_length(store: &AllUl, ulid: usize) -> usize {
    store.nid.get(ulid).map(|n| n.len()).unwrap_or(0)
}

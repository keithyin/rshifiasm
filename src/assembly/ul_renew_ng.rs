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
    sources: &mut [PafList],
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
    sources[(v >> 1) as usize].hits.push(hit);
    let reverse = reverse_hit(&hit);
    sources[(w >> 1) as usize].hits.push(reverse);
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
    paf: &mut [PafList],
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
    let q_list: Vec<usize> = (0..paf[qn0].hits.len())
        .filter(|&i| paf[qn0].hits[i].target_identifier as usize == tn0)
        .collect();
    let t_list: Vec<usize> = (0..paf[tn0].hits.len())
        .filter(|&i| paf[tn0].hits[i].target_identifier as usize == qn0)
        .collect();
    assert!(!q_list.is_empty() && !t_list.is_empty());
    for &qi in &q_list {
        let mut e01 = paf[qn0].hits[qi].clone();
        for &ti in &t_list {
            let mut e10 = paf[tn0].hits[ti].clone();
            for (i, &qn) in a0.iter().enumerate() {
                for (k, &tn) in a1.iter().enumerate() {
                    if i == 0 && k == 0 {
                        continue;
                    }
                    let qn = (qn as u32) as usize;
                    let tn = (tn as u32) as usize;
                    update_qtn(&mut e01, qn as u64, tn as u64);
                    if !sf || dup_paf_check(&paf[qn].hits, &e01) {
                        paf[qn].hits.push(e01.clone());
                    }
                    update_qtn(&mut e10, tn as u64, qn as u64);
                    if !sf || dup_paf_check(&paf[tn].hits, &e10) {
                        paf[tn].hits.push(e10.clone());
                    }
                }
            }
        }
    }
}

/// `renew_paf1` (`gfa_ut.cpp` 17147): all-pairs fake overlaps inside one group.
fn renew_paf1(paf: &mut [PafList], a: &[u64], len: u64) {
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
            paf[qn].hits.push(arc.clone());
        }
    }
}

/// `update_paf` (`gfa_ut.cpp` 17164): rewrites the overlap lists after the
/// read-store growth (every old read gained copies; scaffold reads dropped).
#[allow(clippy::too_many_arguments)]
pub fn update_paf(
    sources: &mut [PafList],
    reverse_sources: &mut [PafList],
    rmap: &[u64],
    rid_n: usize,
    pre_gn: usize,
    ng: &Asg,
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
        for j in 0..sources[k].hits.len() {
            let tn = sources[k].hits[j].target_identifier as usize;
            if tn >= dn {
                continue;
            }
            if (sources[k].hits[j].query_identifier as usize) > tn {
                continue;
            }
            let idx1 = idx[tn];
            let a1: Vec<u64> = clus[(idx1 >> 32) as usize..idx1 as u32 as usize].to_vec();
            renew_paf0(sources, &a0, &a1, &mut buf);
        }
        for j in 0..reverse_sources[k].hits.len() {
            let tn = reverse_sources[k].hits[j].target_identifier as usize;
            if tn >= dn {
                continue;
            }
            if (reverse_sources[k].hits[j].query_identifier as usize) > tn {
                continue;
            }
            let idx1 = idx[tn];
            let a1: Vec<u64> = clus[(idx1 >> 32) as usize..idx1 as u32 as usize].to_vec();
            renew_paf0(reverse_sources, &a0, &a1, &mut buf);
        }
        if a0n > 1 {
            renew_paf1(reverse_sources, &a0, u64::from(ng.sequences[k].length));
        }
    }
}

/// `renew_ng` (`gfa_ut.cpp` 15107): the renewed read graph. `scaffold`
/// carries the `usc_t` records of the scaffold nodes (one per failed
/// `gen_spec_edge`); `h_usg` comes from `u2g_threading`.
///
/// The scaffold/materialisation phases live in this module as well
/// (`fill_scaffolds`, `realloc_rdb_adv`, `update_paf`); the whole reference
/// body runs inside this call.
#[allow(clippy::too_many_arguments)]
pub fn renew_ng(
    h_usg: &crate::assembly::ul2ul::Usg,
    raw_unitigs: &UnitigGraph,
    sg: &Asg,
    sources: &mut [PafList],
    reverse_sources: &mut [PafList],
    coverage_cut: &mut Vec<crate::assembly::hit::CoverageCut>,
    read_to_unitig: &mut crate::assembly::hit::ReadToUnitig,
    reads: &mut crate::reads::ReadStore,
    ul_store: &mut AllUl,
    maximum_hang: i32,
    max_hang_rate: f32,
    minimum_overlap: i64,
    scaffold_len: u64,
    output_file_name: &str,
    bin_file: &str,
    b_mask: &mut crate::assembly::bubble::BubbleLabel,
) -> Asg {
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
                    minimum_overlap as i32,
                    raw_v,
                    raw_w,
                ) {
                    Some(_) => {}
                    None => {
                        first_pass.push((v, w, i as i64));
                        // `asg_seq_set(ng, nocc, 0, 0)`: a scaffold node
                        ng.sequences.push(AsgSeq {
                            length: 0,
                            deleted: false,
                            coverage_class: 0,
                        });
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
    // `ng->seq_vis`
    ng.sequence_visited = vec![0u8; ng.sequences.len() * 2];

    // `if(sa.n > 0) fill_scaffolds(...)`
    let hifi_read_count = reads.total_reads as usize;
    if !scaffold.is_empty() {
        fill_scaffolds(
            &mut scaffold,
            ul_store,
            reads,
            raw_unitigs,
            minimum_overlap,
            output_file_name,
            bin_file,
            hifi_read_count,
        );
    }
    // `realloc_rdb_adv(&(R_INF), &UL_INF, cov, ruI, ext.idx_a.a, ...)`
    realloc_rdb_adv(
        reads,
        ul_store,
        coverage_cut,
        read_to_unitig,
        &idx_a,
        scaffold_len,
        "scaf",
        &mut ng,
        minimum_overlap,
        raw_unitigs,
        &scaffold,
        b_mask,
    );
    // `update_paf(R_INF.paf, R_INF.reverse_paf, ext.idx_a.a, ext.idx_a.n, sg->n_seq, ng)`
    update_paf(
        sources,
        reverse_sources,
        &idx_a,
        idx_a.len(),
        sg.sequences.len(),
        &ng,
    );

    // the arc-building pass
    let mut nocc = ng.read_sequence_count as u64;
    for (i, seq) in h_usg.a.iter().enumerate() {
        if ug1.graph.sequences[i].deleted {
            continue;
        }
        let u = &ug1.unitigs[i];
        for m in 1..u.reads.len() {
            let v = (u.reads[m - 1] >> 32) as u32;
            let w = (u.reads[m] >> 32) as u32;
            let t = crate::assembly::unitig::gen_spec_edge(
                &ng,
                sources,
                maximum_hang,
                max_hang_rate,
                minimum_overlap as i32,
                v,
                w,
            )
            .expect("the unitig-internal edge must exist");
            ng.arcs.push(t);
            let t = crate::assembly::unitig::gen_spec_edge(
                &ng,
                sources,
                maximum_hang,
                max_hang_rate,
                minimum_overlap as i32,
                w ^ 1,
                v ^ 1,
            )
            .expect("the unitig-internal edge must exist");
            ng.arcs.push(t);
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
                if vx > wx {
                    continue;
                }
                if vx == wx && (vx & 1) != 0 {
                    continue;
                }
                match crate::assembly::unitig::gen_spec_edge(
                    &ng,
                    sources,
                    maximum_hang,
                    max_hang_rate,
                    minimum_overlap as i32,
                    vx,
                    wx,
                ) {
                    Some(t) => {
                        ng.arcs.push(t);
                        let t = crate::assembly::unitig::gen_spec_edge(
                            &ng,
                            sources,
                            maximum_hang,
                            max_hang_rate,
                            minimum_overlap as i32,
                            wx ^ 1,
                            vx ^ 1,
                        )
                        .expect("the mirrored edge must exist");
                        ng.arcs.push(t);
                    }
                    None => {
                        let psa = &scaffold[(nocc - ng.read_sequence_count as u64) as usize];
                        let z = nocc << 1;
                        debug_assert_eq!(psa.nid as u64, nocc);
                        debug_assert_eq!(
                            psa.raw_sid,
                            (h_usg.a[(v >> 1) as usize].mm << 1) | (v & 1)
                        );
                        debug_assert_eq!(
                            psa.raw_eid,
                            (h_usg.a[(w >> 1) as usize].mm << 1) | (w & 1)
                        );
                        if ng.sequences[nocc as usize].deleted {
                            push_direct_scaff_node(
                                psa,
                                sources,
                                reverse_sources,
                                &mut ng,
                                maximum_hang,
                                max_hang_rate,
                                minimum_overlap,
                                vx as u64,
                                wx as u64,
                            );
                        } else {
                            push_scaff_node(sources, vx as u64, z as u64, minimum_overlap as u64, &ng);
                            push_scaff_node(sources, z as u64, wx as u64, minimum_overlap as u64, &ng);
                            let t = crate::assembly::unitig::gen_spec_edge(
                                &ng, sources, maximum_hang, max_hang_rate,
                                minimum_overlap as i32, vx, z as u32,
                            )
                            .expect("the scaffold edge must exist");
                            ng.arcs.push(t);
                            let t = crate::assembly::unitig::gen_spec_edge(
                                &ng, sources, maximum_hang, max_hang_rate,
                                minimum_overlap as i32, z as u32 ^ 1, vx ^ 1,
                            )
                            .expect("the scaffold edge must exist");
                            ng.arcs.push(t);
                            let t = crate::assembly::unitig::gen_spec_edge(
                                &ng, sources, maximum_hang, max_hang_rate,
                                minimum_overlap as i32, z as u32, wx,
                            )
                            .expect("the scaffold edge must exist");
                            ng.arcs.push(t);
                            let t = crate::assembly::unitig::gen_spec_edge(
                                &ng, sources, maximum_hang, max_hang_rate,
                                minimum_overlap as i32, wx ^ 1, z as u32 ^ 1,
                            )
                            .expect("the scaffold edge must exist");
                            ng.arcs.push(t);
                        }
                        nocc += 1;
                    }
                }
            }
        }
    }

    ng.cleanup();
    ng.read_sequence_count = ng.sequences.len() as u32;
    eprintln!(
        "[M::renew_ng] nocc::{}, ng->n_seq::{}, sg->n_seq::{}",
        nocc,
        ng.sequences.len(),
        sg.sequences.len()
    );
    assert_eq!(nocc as usize, ng.sequences.len());
    ng
}

/// `renew_R_to_U` (`gfa_ut.cpp` 17582): clears the overlap lists, resets the
/// coverage cuts and re-runs the containment pass.
#[allow(clippy::too_many_arguments)]
pub fn renew_r_to_u(
    sources: &mut [PafList],
    flags: &mut [()],
    reverse_sources: &mut [PafList],
    reverse_flags: &mut [()],
    coverage_cut: &mut [CoverageCut],
    read_to_unitig: &mut ReadToUnitig,
    max_hang: i32,
    min_ovlp: i32,
    reads: &crate::reads::ReadStore,
) {
    // `clear_ma_hit_t_alloc`: clears the `del` flag of every stored hit
    // (the lists survive); `reset_ma_sub_t` with `g == NULL` clears the cut
    // deletions.
    for list in sources.iter_mut() {
        for hit in list.hits.iter_mut() {
            hit.deleted = false;
        }
    }
    for list in reverse_sources.iter_mut() {
        for hit in list.hits.iter_mut() {
            hit.deleted = false;
        }
    }
    for cut in coverage_cut.iter_mut() {
        cut.deleted = false;
    }
    let _ = (flags, reads);
    let mut forward: Vec<Vec<PafHit>> = sources
        .iter()
        .map(|list| list.hits.clone())
        .collect();
    crate::assembly::hit::ma_hit_contained_advance(
        &mut forward,
        coverage_cut,
        read_to_unitig,
        max_hang,
        min_ovlp,
        1.0,
    );
    for (list, hits) in sources.iter_mut().zip(forward) {
        list.hits = hits;
    }
}

/// `AllUl` keeps the UL read names; the reference's `UL_INF` access inside
/// `realloc_rdb_adv` reads the name length only.
pub fn ul_name_length(store: &AllUl, ulid: usize) -> usize {
    store.nid.get(ulid).map(|n| n.len()).unwrap_or(0)
}

/// `get_end_hifi` (`gfa_ut.cpp` 14400): the terminal `len` bases of the read
/// behind vertex `v`, in the vertex orientation.
pub fn get_end_hifi(reads: &crate::reads::ReadStore, v: u32, len: i64) -> Vec<u8> {
    let rid = (v >> 1) as usize;
    let rlen = reads.read_length(rid) as i64;
    assert!(rlen >= len);
    let seq = reads.sequences[rid].unpack();
    let start = if len < rlen {
        if v & 1 == 0 { (rlen - len) as usize } else { 0 }
    } else {
        0
    };
    let mut out: Vec<u8> = seq[start..start + len as usize].to_vec();
    if v & 1 != 0 {
        out.reverse();
        for base in out.iter_mut() {
            *base = crate::correction::target::complement_base(*base);
        }
    }
    out
}

/// `get_ul_subregion` (`gfa_ut.cpp` 14563): the ASCII bases of
/// `[ssp, sep)` of UL read `id`, with the N sites restored.
pub fn get_ul_subregion(
    store: &AllUl,
    id: u32,
    rev: bool,
    ssp: i64,
    sep: i64,
    reset_ns: bool,
) -> Vec<u8> {
    let z = &store.reads[id as usize];
    let sl = (sep - ssp) as usize;
    let mut des = vec![0u8; sl];
    // `r_base` packs four bases per byte, most significant pair first.
    for (i, out) in des.iter_mut().enumerate() {
        let src_i = ssp as usize + i;
        let code = (z.r_base[src_i >> 2] >> (6 - 2 * (src_i & 3))) & 3;
        *out = match code {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
    }
    for &site in &z.n_site {
        if (site as i64) >= ssp && (site as i64) < sep {
            des[(site as i64 - ssp) as usize] = b'N';
        } else if (site as i64) >= sep {
            break;
        }
    }
    if reset_ns {
        // `seq_nt4_table >= 4` is impossible in the two-bit packing; the
        // reference resets any leftover ambiguity to `A`.
        for out in des.iter_mut() {
            if !matches!(*out, b'A' | b'C' | b'G' | b'T') {
                *out = b'A';
            }
        }
    }
    if rev {
        des.reverse();
        for base in des.iter_mut() {
            *base = crate::correction::target::complement_base(*base);
        }
    }
    des
}

/// `push_scaf_bases` (`gfa_ut.cpp` 14620): materialises the scaffold node
/// bases; returns `(packed, ff)`.
#[allow(clippy::too_many_arguments)]
pub fn push_scaf_bases(
    reads: &crate::reads::ReadStore,
    store: &AllUl,
    z: &Usc,
    min_arc_len: i64,
    scaf_len: i64,
    nlen: i64,
    raw_g: &UnitigGraph,
    n_site: &mut Vec<u32>,
) -> (Vec<u8>, i64) {
    let (ff, nlen0) = reset_scaf_node_uinfo(z, min_arc_len, scaf_len);
    debug_assert_eq!(nlen0, nlen);
    if ff == 1 {
        return (vec![0u8; (nlen / 4 + 1) as usize], ff as i64);
    }
    let uv = z.raw_sid;
    let uw = z.raw_eid;
    let u_v = &raw_g.unitigs[(uv >> 1) as usize];
    let u_w = &raw_g.unitigs[(uw >> 1) as usize];
    let sv = if uv & 1 != 0 {
        ((u_v.reads[0] >> 32) as u32) ^ 1
    } else {
        (u_v.reads[u_v.reads.len() - 1] >> 32) as u32
    };
    let sw = if uw & 1 != 0 {
        ((u_w.reads[u_w.reads.len() - 1] >> 32) as u32) ^ 1
    } else {
        (u_w.reads[0] >> 32) as u32
    };
    let mut da: Vec<u8> = vec![0u8; nlen as usize];
    assert!(reads.read_length((sv >> 1) as usize) as i64 > min_arc_len);
    assert!(reads.read_length((sw >> 1) as usize) as i64 > min_arc_len);
    if ff == 2 {
        let s = i64::from(z.raw_sof);
        let e = i64::from(z.raw_eof);
        let ul = e - s;
        assert!(ul >= -min_arc_len);
        da[..min_arc_len as usize].copy_from_slice(&get_end_hifi(reads, sv, min_arc_len));
        da[(min_arc_len + ul) as usize..(min_arc_len + ul + min_arc_len) as usize]
            .copy_from_slice(&get_end_hifi(reads, sw, min_arc_len));
        if ul > 0 {
            let mid = get_ul_subregion(store, z.ulid >> 1, z.ulid & 1 != 0, s, e, true);
            da[min_arc_len as usize..(min_arc_len + ul) as usize].copy_from_slice(&mid);
        }
    } else {
        da[..min_arc_len as usize].copy_from_slice(&get_end_hifi(reads, sv, min_arc_len));
        for base in da[min_arc_len as usize..(min_arc_len + scaf_len) as usize].iter_mut() {
            *base = b'A';
        }
        da[(min_arc_len + scaf_len) as usize..]
            .copy_from_slice(&get_end_hifi(reads, sw, min_arc_len));
    }
    let mut packed = vec![0u8; (nlen / 4 + 1) as usize];
    n_site.clear();
    crate::assembly::ul_store::encode_base_public(&mut packed, &da, n_site);
    (packed, ff as i64)
}

/// `load_scaf_base` (`gfa_ut.cpp` 14330): restores the `dd == 4` UL reads'
/// bases from the `.uidx.ucr.bin` cache.
pub fn load_scaf_base(store: &mut AllUl, output_file_name: &str, bin_file: &str) -> bool {
    use std::io::Read;
    let path = format!("{}.{}.uidx.ucr.bin", output_file_name, bin_file);
    let Ok(mut file) = std::fs::File::open(&path) else {
        return false;
    };
    let mut ucr_path = path.clone();
    let _ = &mut ucr_path;
    let mut len_buf = [0u8; 8];
    loop {
        if file.read_exact(&mut len_buf).is_err() {
            break;
        }
        let rid = u64::from_le_bytes(len_buf);
        let mut word4 = [0u8; 4];
        if file.read_exact(&mut word4).is_err() {
            break;
        }
        let len = u32::from_le_bytes(word4);
        if (rid as usize) >= store.reads.len() {
            break;
        }
        let z = &mut store.reads[rid as usize];
        debug_assert_eq!(z.rlen, len);
        let mut n_word = [0u8; 8];
        if file.read_exact(&mut n_word).is_err() {
            break;
        }
        let n_site_n = u64::from_le_bytes(n_word) as usize;
        let mut sites = vec![0u32; n_site_n];
        for site in sites.iter_mut() {
            let mut w = [0u8; 4];
            if file.read_exact(&mut w).is_err() {
                return false;
            }
            *site = u32::from_le_bytes(w);
        }
        let base_n = ((len as usize) / 4) + (((len as usize) & 3 != 0) as usize);
        let mut packed = vec![0u8; base_n];
        if file.read_exact(&mut packed).is_err() {
            return false;
        }
        if z.dd != 4 {
            continue;
        }
        z.n_site = sites;
        z.r_base = packed;
    }
    true
}

/// `reload_uu` (`gfa_ut.cpp` 14921): reloads the UL store from the round
/// cache against the renewed raw graph.
pub fn reload_uu(
    store: &mut AllUl,
    raw_g: &UnitigGraph,
    output_file_name: &str,
    bin_file: &str,
    hifi_read_count: usize,
) {
    let gfa_name = format!("{}.{}", output_file_name, bin_file);
    crate::assembly::ul_renew::clear_all_ul_t(store);
    if let Ok(Some(cached)) = crate::assembly::ul_store::load_all_ul_t(&gfa_name) {
        *store = cached;
    }
    crate::assembly::ul_store::gen_ul_vec_rid_t(store, hifi_read_count);
    let _ = raw_g;
}

/// `gen_scaffold_id` (`gfa_ut.cpp` 14837): assigns `(ulid, raw_sof,
/// raw_eof)` to one scaffold record from the UL chains.
#[allow(clippy::too_many_arguments)]
pub fn gen_scaffold_id(
    store: &AllUl,
    reads: &crate::reads::ReadStore,
    raw_g: &UnitigGraph,
    z: &mut Usc,
) {
    let v = z.raw_sid;
    let w = z.raw_eid;
    let u_v = &raw_g.unitigs[(v >> 1) as usize];
    let u_w = &raw_g.unitigs[(w >> 1) as usize];
    let sv = if v & 1 != 0 {
        ((u_v.reads[0] >> 32) as u32) ^ 1
    } else {
        (u_v.reads[u_v.reads.len() - 1] >> 32) as u32
    };
    let sw = if w & 1 != 0 {
        ((u_w.reads[u_w.reads.len() - 1] >> 32) as u32) ^ 1
    } else {
        (u_w.reads[0] >> 32) as u32
    };
    let sv_len = reads.read_length((sv >> 1) as usize) as i64;
    let sw_len = reads.read_length((sw >> 1) as usize) as i64;
    let max_ql = -sv_len.min(sw_len);

    let start = store.ridx.idx[(v >> 1) as usize] as usize;
    let end = store.ridx.idx[((v >> 1) + 1) as usize] as usize;
    let mut mmql = i64::MAX;
    let mut rev = 0u32;
    let mut id = u32::MAX;
    let mut m0: Option<crate::assembly::ul_store::UcBlock> = None;
    let mut m1: Option<crate::assembly::ul_store::UcBlock> = None;
    for &word in &store.ridx.occ[start..end] {
        let read_id = (word >> 32) as usize;
        let slot = (word & 0xffff_ffff) as usize;
        let p = &store.reads[read_id].bb[slot];
        debug_assert_eq!(p.hid, v >> 1);
        if store.reads[read_id].dd != 3 {
            continue;
        }
        if p.base || !p.el || p.pchain == 0 {
            continue;
        }
        let uv = (p.hid << 1) | u32::from(p.rev);
        if uv == v && p.aidx == u32::MAX && slot + 1 < store.reads[read_id].bb.len() {
            let n = &store.reads[read_id].bb[slot + 1];
            if n.base || !n.el || n.pchain == 0 || n.pidx != u32::MAX {
                continue;
            }
            let uw2 = (n.hid << 1) | u32::from(n.rev);
            if uw2 == w {
                let ql = i64::from(n.qs) - i64::from(p.qe);
                if ql > max_ql && n.qe > p.qe && n.qs > p.qs && ql < mmql {
                    mmql = ql;
                    rev = 0;
                    m0 = Some(*p);
                    m1 = Some(*n);
                    id = read_id as u32;
                }
            }
        }
        if (uv ^ 1) == v && p.pidx == u32::MAX && slot > 0 {
            let n = &store.reads[read_id].bb[slot - 1];
            if n.base || !n.el || n.pchain == 0 || n.aidx != u32::MAX {
                continue;
            }
            let uw2 = ((n.hid << 1) | u32::from(n.rev)) ^ 1;
            if uw2 == w {
                let ql = i64::from(p.qs) - i64::from(n.qe);
                if ql > max_ql && n.qe < p.qe && n.qs < p.qs && ql < mmql {
                    mmql = ql;
                    rev = 1;
                    m0 = Some(*n);
                    m1 = Some(*p);
                    id = read_id as u32;
                }
            }
        }
    }
    if let (Some(m0), Some(m1)) = (m0, m1) {
        if id != u32::MAX {
            z.ulid = (id << 1) | rev;
            z.raw_sof = m0.qe;
            z.raw_eof = m1.qs;
        }
    }
}

/// `fill_scaffolds` (`gfa_ut.cpp` 14932).
pub fn fill_scaffolds(
    a: &mut [Usc],
    store: &mut AllUl,
    reads: &crate::reads::ReadStore,
    raw_g: &UnitigGraph,
    min_arc_len: i64,
    output_file_name: &str,
    bin_file: &str,
    hifi_read_count: usize,
) {
    reload_uu(store, raw_g, output_file_name, bin_file, hifi_read_count);
    for z in a.iter_mut() {
        gen_scaffold_id(store, reads, raw_g, z);
    }
    let mut occ = 0usize;
    for z in a.iter() {
        if z.ulid == u32::MAX {
            continue;
        }
        if z.raw_eof <= z.raw_sof {
            continue;
        }
        store.reads[(z.ulid >> 1) as usize].dd = 4;
        occ += 1;
    }
    if occ > 0 {
        load_scaf_base(store, output_file_name, bin_file);
    }
}

/// `push_direct_scaff_node` (`gfa_ut.cpp` 15019): the direct link between
/// two vertices bridged by a scaffold record whose UL gap is short enough.
#[allow(clippy::too_many_arguments)]
pub fn push_direct_scaff_node(
    psa: &Usc,
    sources: &mut [PafList],
    reverse_sources: &mut [PafList],
    ng: &mut Asg,
    maximum_hang: i32,
    max_hang_rate: f32,
    minimum_overlap: i64,
    vx: u64,
    wx: u64,
) {
    let uol = i64::from(psa.raw_sof) - i64::from(psa.raw_eof);
    assert!(uol >= minimum_overlap);
    let mut r = 1i32;
    let mut t0 = None;
    let mut rc0: Option<(AsgArc, usize)> = None;
    let mut t1 = None;
    let mut rc1: Option<(AsgArc, usize)> = None;
    if let Some(hit) = gen_spec_rc_edge(
        ng,
        reverse_sources,
        maximum_hang,
        max_hang_rate,
        minimum_overlap as i32,
        vx as u32,
        wx as u32,
    ) {
        let t = hit.0;
        let dif = if uol > i64::from(t.overlap) {
            uol - i64::from(t.overlap)
        } else {
            i64::from(t.overlap) - uol
        };
        if (dif as f64 > uol as f64 * 0.05) || (dif as f64 > t.overlap as f64 * 0.05) {
            if dif > 16 {
                r = -1;
            }
        }
        if r >= 0 {
            t0 = Some(t);
            rc0 = Some(hit);
        }
    } else {
        r = -1;
    }
    if r >= 0 {
        if let Some(hit) = gen_spec_rc_edge(
            ng,
            reverse_sources,
            maximum_hang,
            max_hang_rate,
            minimum_overlap as i32,
            (wx ^ 1) as u32,
            (vx ^ 1) as u32,
        ) {
            let t = hit.0;
            let dif = if uol > i64::from(t.overlap) {
                uol - i64::from(t.overlap)
            } else {
                i64::from(t.overlap) - uol
            };
            if (dif as f64 > uol as f64 * 0.05) || (dif as f64 > t.overlap as f64 * 0.05) {
                if dif > 16 {
                    r = -1;
                }
            }
            if r >= 0 {
                t1 = Some(t);
                rc1 = Some(hit);
            }
        } else {
            r = -1;
        }
    }

    if r >= 0 {
        let (t0, rc0) = t0.zip(rc0).expect("the first edge must exist");
        ng.arcs.push(t0);
        let (t1, rc1) = t1.zip(rc1).expect("the second edge must exist");
        ng.arcs.push(t1);
        // `add_ma_hit_t_alloc(&(R_INF.paf[Get_qn((*rc0))]), rc0)`: rc0 is
        // the hit found in the reverse list, copied onto the forward list.
        let hit0 = reverse_sources[(vx >> 1) as usize].hits[rc0.1].clone();
        sources[hit0.query_identifier as usize].hits.push(hit0);
        let hit1 = reverse_sources[((wx ^ 1) >> 1) as usize].hits[rc1.1].clone();
        sources[hit1.query_identifier as usize].hits.push(hit1);
    } else {
        push_scaff_node(sources, vx, wx, uol as u64, ng);
        let t = crate::assembly::unitig::gen_spec_edge(
            ng,
            sources,
            maximum_hang,
            max_hang_rate,
            minimum_overlap as i32,
            vx as u32,
            wx as u32,
        )
        .expect("the pushed edge must exist");
        ng.arcs.push(t);
        let t = crate::assembly::unitig::gen_spec_edge(
            ng,
            sources,
            maximum_hang,
            max_hang_rate,
            minimum_overlap as i32,
            (wx ^ 1) as u32,
            (vx ^ 1) as u32,
        )
        .expect("the pushed edge must exist");
        ng.arcs.push(t);
    }
}

/// `realloc_rdb_adv` (`gfa_ut.cpp` 14663): extends the read store with the
/// scaffold pseudo-reads and resets the coverage/label tables.
#[allow(clippy::too_many_arguments)]
pub fn realloc_rdb_adv(
    reads: &mut crate::reads::ReadStore,
    ul_store: &AllUl,
    coverage_cut: &mut Vec<crate::assembly::hit::CoverageCut>,
    read_to_unitig: &mut crate::assembly::hit::ReadToUnitig,
    rmap: &[u64],
    scaf_len: u64,
    scaf_id: &str,
    ng: &mut Asg,
    min_ovlp: i64,
    raw_g: &UnitigGraph,
    a: &[Usc],
    b_mask: &mut crate::assembly::bubble::BubbleLabel,
) {
    let rid_n0 = reads.total_reads as usize;
    let rid_n = rmap.len();
    reads.total_reads = rid_n as u64;
    eprintln!("+[M::realloc_rdb_adv] rid_n0::{}, rid_n::{}", rid_n0, rid_n);
    for i in rid_n0..rid_n {
        if rmap[i] != u64::MAX {
            let src = rmap[i] as usize;
            let name = reads.read_names[src].clone();
            let length = reads.read_lengths[src];
            let sequence = reads.sequences[src].clone();
            let trio = reads.trio_labels[src];
            let cut = coverage_cut[src];
            reads.read_names.push(name);
            reads.read_lengths.push(length);
            reads.sequences.push(sequence);
            reads.trio_labels.push(trio);
            coverage_cut.push(cut);
        } else {
            let z = &a[i - ng.read_sequence_count as usize];
            debug_assert_eq!(z.nid as usize, i);
            let (ff, nlen) = reset_scaf_node_uinfo(z, min_ovlp, scaf_len as i64);
            let name = if ff == 2 {
                let ulid = (z.ulid >> 1) as usize;
                String::from_utf8_lossy(&ul_store.nid[ulid]).to_string()
            } else {
                scaf_id.to_string()
            };
            reads.read_names.push(name);
            reads.read_lengths.push(nlen as u64);
            reads.trio_labels.push(crate::assembly::unitig::AMBIGU);
            coverage_cut.push(crate::assembly::hit::CoverageCut {
                start: 0,
                end: nlen as u32,
                deleted: false,
                class: 0,
            });
            ng.sequences[i].length = nlen as u32;
            let mut n_site: Vec<u32> = Vec::new();
            let (packed, ff) = push_scaf_bases(
                reads,
                ul_store,
                z,
                min_ovlp,
                scaf_len as i64,
                nlen,
                raw_g,
                &mut n_site,
            );
            reads.sequences.push(crate::sequence::encoding::PackedSequence {
                packed_bases: packed,
                length: nlen as u64,
                ambiguous_positions: n_site.iter().map(|&s| u64::from(s)).collect(),
            });
            if ff == 1 {
                ng.sequences[i].deleted = true;
            }
        }
        reads.total_bases += reads.read_lengths[i];
        ng.sequences[i].length = reads.read_lengths[i] as u32;
    }
    *read_to_unitig = crate::assembly::hit::ReadToUnitig::from_entries(vec![u32::MAX; rid_n]);
    b_mask.reset(ng.sequences.len(), 0, 0);
    eprintln!("-[M::realloc_rdb_adv] rid_n0::{}, rid_n::{}", rid_n0, rid_n);
}

/// `shrink_1b` (`gfa_ut.cpp` 7837): trims one base off a full-length end
/// block when the neighbour still dominates it.
pub fn shrink_1b(
    ug: &UnitigGraph,
    reads: &crate::reads::ReadStore,
    z: &mut crate::assembly::ul::UcBlock,
    lim: Option<&crate::assembly::ul::UcBlock>,
    is_forward: bool,
) {
    if z.ts != 0 || z.te != ug.graph.sequences[z.hid as usize].length {
        return;
    }
    let bc = *z;
    if z.ts + 1 < z.te && z.qs + 1 < z.qe {
        let off = crate::assembly::ul::get_offset_adjust(1, z.te - z.ts, z.qe - z.qs);
        if is_forward {
            if !z.rev {
                z.ts += 1;
                z.qs += off;
            } else {
                z.te -= 1;
                z.qs += off;
            }
            if let Some(lim) = lim {
                if !(z.qs <= lim.qs && z.qe <= lim.qe) {
                    *z = bc;
                }
            }
        } else {
            if !z.rev {
                z.te -= 1;
                z.qe -= off;
            } else {
                z.ts += 1;
                z.qe -= off;
            }
            if let Some(lim) = lim {
                if !(z.qs >= lim.qs && z.qe >= lim.qe) {
                    *z = bc;
                }
            }
        }
        let cover_before = crate::assembly::ul_renew::ugl_cover_check(
            u64::from(bc.ts),
            u64::from(bc.te),
            &ug.unitigs[bc.hid as usize],
            reads,
        );
        let cover_after = crate::assembly::ul_renew::ugl_cover_check(
            u64::from(z.ts),
            u64::from(z.te),
            &ug.unitigs[z.hid as usize],
            reads,
        );
        if cover_before && !cover_after {
            *z = bc;
        }
    }
}

/// `renew_ul_vec_t` (`gfa_ut.cpp` 7848): re-derives the query coordinates of
/// a rewritten block chain from the `pidx`/`pdis` links.
pub fn renew_ul_vec_t(x: &mut crate::assembly::ul_store::UlReadVec, ug: &UnitigGraph) {
    if x.bb.is_empty() {
        return;
    }
    let bn = x.bb.len();
    let mut qs: i64 = 0;
    let mut l = 0usize;
    let mut k = 1usize;
    loop {
        let at_end = k == bn;
        if at_end || x.bb[k].pidx == u32::MAX {
            if k > l {
                let mut dd: i64 = 0;
                let mut z = k as i64 - 1;
                while z >= l as i64 {
                    let zu = z as usize;
                    if zu > l {
                        debug_assert_eq!(x.bb[zu].pidx as usize, zu - 1);
                        dd += i64::from(x.bb[zu].pdis);
                    } else {
                        debug_assert_eq!(x.bb[zu].pidx, u32::MAX);
                        dd += i64::from(ug.graph.sequences[x.bb[zu].hid as usize].length);
                    }
                    dd -= i64::from(
                        ug.graph.sequences[x.bb[zu].hid as usize].length
                            - (x.bb[zu].te - x.bb[zu].ts),
                    );
                    z -= 1;
                }
                if dd < 0 {
                    dd = 0;
                }
                let mut qe = qs + dd;
                let mut dd_next = 0i64;
                if k < bn {
                    dd_next = qe + i64::from(x.bb[k].qs) - i64::from(x.bb[k - 1].qe);
                    if dd_next < 0 {
                        dd_next = 0;
                    }
                }
                let mut z = k as i64 - 1;
                while z >= l as i64 {
                    let zu = z as usize;
                    qs = qe - i64::from(x.bb[zu].te - x.bb[zu].ts);
                    if qs < 0 {
                        qs = 0;
                    }
                    x.bb[zu].qs = qs as u32;
                    x.bb[zu].qe = qe as u32;
                    if zu > l {
                        qe += i64::from(
                            ug.graph.sequences[x.bb[zu].hid as usize].length
                                - (x.bb[zu].te - x.bb[zu].ts),
                        );
                        qe -= i64::from(x.bb[zu].pdis);
                        if qe < 0 {
                            qe = 0;
                        }
                    }
                    z -= 1;
                }
                if k < bn {
                    qs = dd_next;
                }
            }
            l = k;
        }
        if at_end {
            break;
        }
        k += 1;
    }
    x.rlen = x.bb[bn - 1].qe;
}

/// `shrink_ul0` (`gfa_ut.cpp` 7742): splits one UL read's block chain at the
/// disconnected `pidx`/`pdis` links, appending the earlier pieces as new
/// reads and rewriting the last piece in place.
pub fn shrink_ul0(
    uls: &mut AllUl,
    str_: &mut crate::assembly::ul_integer::UlStr,
    id: usize,
    buf: &mut Vec<u64>,
    ug: &UnitigGraph,
    reads: &crate::reads::ReadStore,
) {
    if (str_.cn as usize) < 2 {
        return;
    }
    buf.clear();
    let mut bl = 0u32;
    let mut p_k = u32::MAX;
    let mut pv = u32::MAX;
    let mut k = 0usize;
    while k < str_.cn as usize {
        let xi = uls.reads[id].bb[(str_.a[k] >> 32) as usize];
        let mut is_conn = false;
        let c_k = k as u32;
        let cv = (str_.a[k] as u32) ^ 1;
        debug_assert_eq!(
            (xi.hid << 1) + u32::from(xi.rev),
            str_.a[k] as u32
        );
        if p_k != u32::MAX {
            if xi.pidx == (str_.a[p_k as usize] >> 32) as u32 {
                let arcs = ug.graph.arcs_of(cv);
                let mut i = 0usize;
                while i < arcs.len() {
                    if arcs[i].deleted {
                        i += 1;
                        continue;
                    }
                    if arcs[i].v == pv {
                        break;
                    }
                    i += 1;
                }
                if i < arcs.len() {
                    let s = arcs[i].ul as u32;
                    if s == xi.pdis {
                        is_conn = true;
                    } else {
                        let d = if s >= xi.pdis { s - xi.pdis } else { xi.pdis - s };
                        let mm = s.max(xi.pdis);
                        if (d as f64) <= (mm as f64 * 0.08) || d <= 512 {
                            is_conn = true;
                        }
                    }
                }
            } else if xi.pidx == u32::MAX {
                is_conn = true;
            }
        }
        if is_conn {
            is_conn = false;
            assert!(k > 0);
            let yi = uls.reads[id].bb[(str_.a[k - 1] >> 32) as usize];
            if xi.qs >= yi.qs && xi.qe >= yi.qe {
                is_conn = true;
            }
        }
        if is_conn {
            bl += 1;
        } else {
            if bl > 0 {
                let word = (((k as u64) - 1 - bl as u64) << 32) + bl as u64 + 1;
                buf.push(word);
            }
            bl = 0;
        }
        p_k = c_k;
        pv = cv;
        k += 1;
    }
    if bl > 0 {
        let word = (((k as u64) - 1 - bl as u64) << 32) + bl as u64 + 1;
        buf.push(word);
    }

    if buf.is_empty() {
        str_.cn = 0;
        return;
    }
    let copy_group = |src: &[crate::assembly::ul::UcBlock],
                      str_a: &[u64],
                      slots: &[u64],
                      out: &mut crate::assembly::ul_store::UlReadVec,
                      base_nid: usize| {
        let mut m = 0usize;
        for (j, &kk) in slots.iter().enumerate() {
            let blk = src[(str_a[kk as usize] >> 32) as usize];
            let pidx = blk.pidx;
            let aidx = blk.aidx;
            let pdis = blk.pdis;
            out.bb[m] = blk;
            if j == 0 {
                out.bb[m].pidx = u32::MAX;
                out.bb[m].pdis = u32::MAX;
            } else if pidx != u32::MAX {
                out.bb[m].pidx = m as u32 - 1;
                out.bb[m].pdis = pdis;
            }
            if j + 1 == slots.len() {
                out.bb[m].aidx = u32::MAX;
            } else if aidx != u32::MAX {
                out.bb[m].aidx = m as u32 + 1;
            }
            m += 1;
        }
        out.bb.truncate(m);
        assert!(out.bb.len() > 1);
        let lim1 = if out.bb.len() >= 2 {
            Some(out.bb[1])
        } else {
            None
        };
        shrink_1b(ug, reads, &mut out.bb[0], lim1.as_ref(), true);
        let lim2 = if out.bb.len() >= 2 {
            Some(out.bb[out.bb.len() - 2])
        } else {
            None
        };
        let last = out.bb.len() - 1;
        shrink_1b(ug, reads, &mut out.bb[last], lim2.as_ref(), false);
        let d = out.bb[0].qs;
        for blk in out.bb.iter_mut() {
            blk.qs -= d;
            blk.qe -= d;
        }
        out.rlen = out.bb[out.bb.len() - 1].qe;
        let _ = base_nid;
        renew_ul_vec_t(out, ug);
    };

    if buf.len() > 1 {
        for i in 0..buf.len() - 1 {
            let s = (buf[i] >> 32) as usize;
            let e = s + (buf[i] & 0xffff_ffff) as usize;
            // new name
            let name = uls.nid[id].clone();
            uls.nid.push(name);
            // new overlap list
            let mut x = crate::assembly::ul_store::UlReadVec::default();
            x.bb = vec![crate::assembly::ul::UcBlock::default(); e - s];
            let slots: Vec<u64> = (s..e).map(|kk| kk as u64).collect();
            // The source blocks are read from the untouched read; the new
            // group is written into `x`.
            let src = uls.reads[id].bb.clone();
            copy_group(&src, &str_.a, &slots, &mut x, uls.nid.len());
            uls.reads.push(x);
        }
    }
    // the last group updates in place
    let i = buf.len() - 1;
    let s = (buf[i] >> 32) as usize;
    let e = s + (buf[i] & 0xffff_ffff) as usize;
    let slots: Vec<u64> = (s..e).map(|kk| kk as u64).collect();
    let mut x = std::mem::take(&mut uls.reads[id]);
    // The reference copies block `m` from `x->bb.a[z]` in place; `m <= z`
    // always holds, so a cloned source reproduces it safely here.
    let src = x.bb.clone();
    x.bb.resize(e - s, crate::assembly::ul::UcBlock::default());
    copy_group(&src, &str_.a, &slots, &mut x, uls.nid.len());
    uls.reads[id] = x;
}

/// `shrink_uls` (`gfa_ut.cpp` 7755) plus `rebuid_idx` (7653): the whole
/// post-correction compaction.
pub fn shrink_uls(
    uls: &mut AllUl,
    str_idx: &mut crate::assembly::ul_integer::UlStrIdx,
    ug: &mut UnitigGraph,
    reads: &crate::reads::ReadStore,
    hifi_read_count: usize,
) {
    let ul_n = uls.reads.len();
    let mut buf: Vec<u64> = Vec::new();
    for k in 0..ul_n {
        let mut str_ = std::mem::take(&mut str_idx.str[k]);
        shrink_ul0(uls, &mut str_, k, &mut buf, ug, reads);
        str_idx.str[k] = str_;
    }
    // `rebuid_idx`
    crate::assembly::ul_renew::filter_ul_ug(uls, ug, reads);
    crate::assembly::ul_store::gen_ul_vec_rid_t(uls, hifi_read_count);
    crate::assembly::ul_renew::update_ug_arch_ul_mul(ug, uls);
    *str_idx = crate::assembly::ul_integer::init_ul_str_idx(&uls.reads, ug.unitigs.len());
}

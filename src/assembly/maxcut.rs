//! The max-cut haplotype partitioner of `rcut.cpp`.
//!
//! The reference solves the same problem twice: a block-level pass
//! (`mb_solve_core`) that groups unitigs into diploid blocks and balances the
//! spins inside each block, then a unitig-level pass (`mc_solve_core_adv`)
//! that optimises the spins of every connected component, optionally
//! rearranging bubble chains from `bubble_type`. Only the paths reached from
//! `mc_solve(ovlp=NULL, t_ch=NULL, ...)` are ported, which is exactly the
//! trio-dual call `mc_solve(NULL, NULL, &in, ug, sg, 0.8, NULL, 0, s->s, 1,
//! bub, ref, 0, 0)`.
//!
//! Every loop is written for a single thread: `kt_for` with `n_threads == 1`
//! invokes the worker once per item in ascending order, which is what the
//! byte-identical `-t1` gate compares against.

use crate::assembly::asg::MaUtg;
use crate::assembly::bubble_chain::{get_bubbles, BubbleType};
use crate::assembly::hap_cov::{UTrans, N_HET, P_HET};
use crate::assembly::purge_dups::{
    get_specific_hap_overlap, set_reverse_hap_overlap, HapOverlap, SET_TRIO,
};
use crate::assembly::unitig::{UnitigGraph, FATHER, MOTHER};
use crate::container::ksort::{radix_sort_by_key, radix_sort_by_key_32};

/// `mc_node_t`: a haplotype spin, `+1`/`-1` for the two haplotypes and `0`
/// when unassigned.
type McNode = i8;

#[derive(Debug, Clone, Copy)]
struct McOpt {
    max_iter: i32,
    n_perturb: i32,
    n_b_perturb: i32,
    n_s_perturb: i32,
    f_perturb: f64,
    seed: u64,
}

/// `mc_opt_init`: the perturbation parameters, seeded from `--seed`,
/// `--n-perturb` and `--f-perturb`.
fn mc_opt_init(n_perturb: i32, f_perturb: f64, seed: u64) -> McOpt {
    McOpt {
        max_iter: 1000,
        n_perturb,
        f_perturb,
        seed,
        n_s_perturb: n_perturb,
        n_b_perturb: n_perturb / 2,
    }
}

/// `kr_splitmix64`.
fn kr_splitmix64(x: u64) -> u64 {
    let z = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    let z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// `kr_drand_r`: a double in `[0, 1)` built from the top 52 RNG bits.
fn kr_drand_r(x: &mut u64) -> f64 {
    *x = kr_splitmix64(*x);
    f64::from_bits(0x3FFu64 << 52 | (*x >> 12)) - 1.0
}

/// `mc_edge_t`.
#[derive(Debug, Clone, Copy, Default)]
struct McEdge {
    x: u64,
    w: f64,
}

/// `mc_match_t`.
#[derive(Debug, Clone, Default)]
struct McMatch {
    idx: Vec<u64>,
    ma: Vec<McEdge>,
    cc: Vec<u64>,
    n_seq: u32,
}

/// `mb_node_t`.
#[derive(Debug, Clone, Copy, Default)]
struct MbNode {
    a: [u32; 2],
    occ: [u32; 2],
    s: [McNode; 2],
    z: [f64; 4],
}

/// `mb_nodes_t`.
#[derive(Debug, Clone, Default)]
struct MbNodes {
    bid: Vec<u32>,
    idx: Vec<u32>,
    u: Vec<MbNode>,
}

/// `mb_edge_t`.
#[derive(Debug, Clone, Copy, Default)]
struct MbEdge {
    x: u64,
    w: [f64; 4],
}

/// `mb_match_t`.
#[derive(Debug, Clone, Default)]
struct MbMatch {
    idx: Vec<u64>,
    ma: Vec<MbEdge>,
    cc: Vec<u64>,
    n_seq: u32,
}

/// `mb_g_t`.
#[derive(Debug, Clone, Default)]
struct MbG {
    u: MbNodes,
    e: MbMatch,
}

/// `mc_g_t`.
#[derive(Debug, Clone, Default)]
struct McG {
    s: Vec<McNode>,
    e: McMatch,
}

fn ma_x(x: u64) -> u32 {
    (x >> 32) as u32
}

fn ma_y(x: u64) -> u32 {
    x as u32
}

fn mc_span(ma: &McMatch, sid: u32) -> (usize, usize) {
    let packed = ma.idx[sid as usize];
    ((packed >> 32) as usize, packed as u32 as usize)
}

fn mb_span(ma: &MbMatch, sid: u32) -> (usize, usize) {
    let packed = ma.idx[sid as usize];
    ((packed >> 32) as usize, packed as u32 as usize)
}

fn u_trans_span(index: &[u64], qn: u32) -> (usize, usize) {
    let packed = index[qn as usize];
    ((packed >> 32) as usize, packed as u32 as usize)
}

fn get_mc_edge(ma: &McMatch, sid1: u32, sid2: u32) -> Option<usize> {
    let (off, count) = mc_span(ma, sid1);
    (off..off + count).find(|&k| ma_y(ma.ma[k].x) == sid2)
}

fn get_mb_edge(ma: &MbMatch, sid1: u32, sid2: u32) -> Option<usize> {
    let (off, count) = mb_span(ma, sid1);
    (off..off + count).find(|&k| ma_y(ma.ma[k].x) == sid2)
}

/// `mc_merge_dup`: sums the weights of the edges sharing an endpoint pair.
/// The vector must already be sorted by `x`.
fn mc_merge_dup(ma: &mut Vec<McEdge>) {
    let n = ma.len();
    let mut out: Vec<McEdge> = Vec::with_capacity(n);
    let mut st = 0usize;
    let mut i = 1usize;
    while i <= n {
        if i == n || ma[i].x != ma[st].x {
            if i - st > 1 {
                let mut w = 0.0f64;
                for edge in &ma[st..i] {
                    w += edge.w;
                }
                let mut entry = ma[st];
                entry.w = w;
                out.push(entry);
            } else {
                out.push(ma[st]);
            }
            st = i;
        }
        i += 1;
    }
    *ma = out;
}

/// `mb_merge_dup`: the four-weight variant of [`mc_merge_dup`].
fn mb_merge_dup(ma: &mut Vec<MbEdge>) {
    let n = ma.len();
    let mut out: Vec<MbEdge> = Vec::with_capacity(n);
    let mut st = 0usize;
    let mut i = 1usize;
    while i <= n {
        if i == n || ma[i].x != ma[st].x {
            if i - st > 1 {
                let mut w = [0.0f64; 4];
                for edge in &ma[st..i] {
                    w[0] += edge.w[0];
                    w[1] += edge.w[1];
                    w[2] += edge.w[2];
                    w[3] += edge.w[3];
                }
                let mut entry = ma[st];
                entry.w = w;
                out.push(entry);
            } else {
                out.push(ma[st]);
            }
            st = i;
        }
        i += 1;
    }
    *ma = out;
}

/// `mc_edges_idx`/`mb_edges_idx`: the `<start>:32 | <count>` span of every
/// query vertex.
fn mc_edges_idx(idx: &mut Vec<u64>, ma: &[McEdge], n_seq: u32) {
    *idx = vec![0u64; n_seq as usize];
    let n = ma.len();
    let mut st = 0usize;
    for i in 1..=n {
        if i == n || ma_x(ma[i].x) != ma_x(ma[st].x) {
            idx[ma_x(ma[st].x) as usize] = ((st as u64) << 32) | (i - st) as u64;
            st = i;
        }
    }
}

fn mb_edges_idx(idx: &mut Vec<u64>, ma: &[MbEdge], n_seq: u32) {
    *idx = vec![0u64; n_seq as usize];
    let n = ma.len();
    let mut st = 0usize;
    for i in 1..=n {
        if i == n || ma_x(ma[i].x) != ma_x(ma[st].x) {
            idx[ma_x(ma[st].x) as usize] = ((st as u64) << 32) | (i - st) as u64;
            st = i;
        }
    }
}

/// `mc_edges_symm`: every edge keeps the weight of the stronger direction;
/// self edges and edges without a mirror are dropped.
fn mc_edges_symm(ma: &mut McMatch) -> u32 {
    let mut del = vec![false; ma.ma.len()];
    let mut removed = 0u32;
    for i in 0..ma.ma.len() {
        let m = ma.ma[i];
        if ma_x(m.x) == ma_y(m.x) {
            del[i] = true;
            removed += 1;
            continue;
        }
        match get_mc_edge(ma, ma_y(m.x), ma_x(m.x)) {
            None => {
                del[i] = true;
                removed += 1;
            }
            Some(t) => {
                if m.w >= ma.ma[t].w {
                    ma.ma[t].w = m.w;
                } else {
                    ma.ma[i].w = ma.ma[t].w;
                }
            }
        }
    }
    if removed > 0 {
        let mut k = 0usize;
        for i in 0..ma.ma.len() {
            if !del[i] {
                ma.ma[k] = ma.ma[i];
                k += 1;
            }
        }
        ma.ma.truncate(k);
        let mut idx = std::mem::take(&mut ma.idx);
        mc_edges_idx(&mut idx, &ma.ma, ma.n_seq);
        ma.idx = idx;
    }
    removed
}

/// `normalize_mb_edge`: mirrors the partner's four weights when it is the
/// stronger direction.
fn mb_edges_symm(ma: &mut MbMatch) -> u32 {
    let mut del = vec![false; ma.ma.len()];
    let mut removed = 0u32;
    for i in 0..ma.ma.len() {
        let m = ma.ma[i];
        if ma_x(m.x) == ma_y(m.x) {
            del[i] = true;
            removed += 1;
            continue;
        }
        match get_mb_edge(ma, ma_y(m.x), ma_x(m.x)) {
            None => {
                del[i] = true;
                removed += 1;
            }
            Some(t) => {
                if m.w[0] >= ma.ma[t].w[0] {
                    ma.ma[t].w = [m.w[0], m.w[2], m.w[1], m.w[3]];
                } else {
                    let other = ma.ma[t];
                    ma.ma[i].w = [other.w[0], other.w[2], other.w[1], other.w[3]];
                }
            }
        }
    }
    if removed > 0 {
        let mut k = 0usize;
        for i in 0..ma.ma.len() {
            if !del[i] {
                ma.ma[k] = ma.ma[i];
                k += 1;
            }
        }
        ma.ma.truncate(k);
        let mut idx = std::mem::take(&mut ma.idx);
        mb_edges_idx(&mut idx, &ma.ma, ma.n_seq);
        ma.idx = idx;
    }
    removed
}

/// `init_mc_g_t`: the spin array and an empty edge set.
fn init_mc_g_t(n_seq: u32, s: &[McNode]) -> McG {
    let mut spins = s.to_vec();
    if spins.len() < n_seq as usize {
        spins.resize(n_seq as usize, 0);
    }
    McG {
        s: spins,
        e: McMatch::default(),
    }
}

/// `update_mc_edges` with no `hap_overlaps_list` and no `trans_chain`.
///
/// Without a `t_ch` the heterozygosity interval filter is skipped, so every
/// live transition becomes an edge weighted by `nw`; zero-weight edges are
/// dropped before the radix sort, exactly like the reference.
fn update_mc_edges(mg: &mut McG, ta: &[UTrans], is_sys: bool) {
    let n_seq = mg.e.n_seq;
    let mut edges: Vec<McEdge> = Vec::with_capacity(ta.len());
    for transition in ta {
        if transition.deleted {
            continue;
        }
        edges.push(McEdge {
            x: ((transition.query_identifier as u64) << 32)
                | transition.target_identifier as u64,
            w: transition.network_weight,
        });
    }
    edges.retain(|edge| edge.w != 0.0);
    radix_sort_by_key(&mut edges, &|edge: &McEdge| edge.x);
    mc_merge_dup(&mut edges);
    let mut idx = Vec::new();
    mc_edges_idx(&mut idx, &edges, n_seq);
    mg.e.ma = edges;
    mg.e.idx = idx;
    if is_sys {
        mc_edges_symm(&mut mg.e);
    }
}

/// `mc_interval_t`: one run of reads of a unitig with the same heterozygosity
/// flag, in both node and base coordinates.
#[derive(Debug, Clone, Copy, Default)]
struct McInterval {
    b_s: u32,
    b_e: u32,
    n_s: u32,
    n_e: u32,
    u_id: u32,
    hs: u8,
}

/// The `p`/`p_idx` interval tables `update_mc_edges` builds from `t_ch`.
fn build_mc_intervals(
    unitig_reads: &[MaUtg],
    read_lengths: &[u32],
    het: &[u8],
) -> (Vec<McInterval>, Vec<u32>) {
    let mut p: Vec<McInterval> = Vec::new();
    let mut p_idx: Vec<u32> = Vec::with_capacity(unitig_reads.len() + 1);
    p_idx.push(0);
    for (v, u) in unitig_reads.iter().enumerate() {
        let reads = &u.reads;
        let mut l = 0usize;
        let mut offset: u32 = 0;
        let mut l_pos: u32 = 0;
        for k in 1..=reads.len() {
            let end = k == reads.len()
                || het[(reads[k] >> 33) as usize] != het[(reads[l] >> 33) as usize];
            if end {
                p.push(McInterval {
                    u_id: v as u32,
                    hs: het[(reads[l] >> 33) as usize],
                    b_s: l_pos,
                    b_e: offset
                        .wrapping_add(read_lengths[(reads[k - 1] >> 33) as usize])
                        .wrapping_sub(1),
                    n_s: l as u32,
                    n_e: (k - 1) as u32,
                });
                l = k;
                l_pos = offset.wrapping_add(reads[k - 1] as u32);
            }
            offset = offset.wrapping_add(reads[k - 1] as u32);
        }
        p_idx.push(p.len() as u32);
    }
    (p, p_idx)
}

/// The per-endpoint het/hom run-length filter shared by the `ha` and `ta`
/// loops of `update_mc_edges`.
fn mc_interval_filter(
    p: &[McInterval],
    p_idx: &[u32],
    polyploidy: i32,
    qn: u32,
    qs: u32,
    qe: u32,
    f_rate: f64,
) -> bool {
    if p.is_empty() || p_idx.len() < 2 {
        return true;
    }
    let beg = p_idx[qn as usize] as usize;
    let occ = p_idx[qn as usize + 1] as usize - beg;
    let mut het_len: u64 = 0;
    let mut hom_len: u64 = 0;
    for interval in &p[beg..beg + occ] {
        let a_s = interval.b_s;
        let a_e = interval.b_e;
        let o_len = if qe.min(a_e) >= qs.max(a_s) {
            u64::from(qe.min(a_e) - qs.max(a_s) + 1)
        } else {
            0
        };
        if hom_len + het_len > 0 && o_len == 0 {
            break;
        }
        if o_len == 0 {
            continue;
        }
        if interval.hs == N_HET {
            hom_len += o_len;
        } else if polyploidy <= 2 && (interval.hs & P_HET) != 0 {
            hom_len += o_len;
        } else {
            het_len += o_len;
        }
    }
    het_len > ((het_len + hom_len) as f64 * f_rate) as u64
}

/// `update_mc_edges` with a `hap_overlaps_list` and an optional `trans_chain`.
///
/// Only the hap-overlap edges are built here: the purge call passes no
/// transition table. The het/hom interval filter applies to both endpoints and
/// is skipped entirely without a `t_ch`.
fn update_mc_edges_ovlp(
    mg: &mut McG,
    ha: &[Vec<HapOverlap>],
    het: Option<&[u8]>,
    read_lengths: &[u32],
    unitig_reads: &[MaUtg],
    f_rate: f64,
    polyploidy: i32,
    is_sys: bool,
) {
    let (p, p_idx) = match het {
        Some(het) => build_mc_intervals(unitig_reads, read_lengths, het),
        None => (Vec::new(), Vec::new()),
    };
    let n_seq = mg.e.n_seq;
    let mut edges: Vec<McEdge> = Vec::new();
    for list in ha {
        for ov in list {
            if ov.score <= 0 || ov.x_uid == ov.y_uid {
                continue;
            }
            if !p.is_empty() && !p_idx.is_empty() {
                if !mc_interval_filter(
                    &p,
                    &p_idx,
                    polyploidy,
                    ov.x_uid,
                    ov.x_beg_pos,
                    ov.x_end_pos.wrapping_sub(1),
                    f_rate,
                ) {
                    continue;
                }
                if !mc_interval_filter(
                    &p,
                    &p_idx,
                    polyploidy,
                    ov.y_uid,
                    ov.y_beg_pos,
                    ov.y_end_pos.wrapping_sub(1),
                    f_rate,
                ) {
                    continue;
                }
            }
            edges.push(McEdge {
                x: ((ov.x_uid as u64) << 32) | ov.y_uid as u64,
                w: ov.score as f64,
            });
        }
    }
    edges.retain(|edge| edge.w != 0.0);
    radix_sort_by_key(&mut edges, &|edge: &McEdge| edge.x);
    mc_merge_dup(&mut edges);
    let mut idx = Vec::new();
    mc_edges_idx(&mut idx, &edges, n_seq);
    mg.e.ma = edges;
    mg.e.idx = idx;
    if is_sys {
        mc_edges_symm(&mut mg.e);
    }
}

/// `filter_ovlp_by_mc`: marks the overlaps that cross the partition.
fn filter_ovlp_by_mc(ha: &mut [Vec<HapOverlap>], ma: &McMatch, s: &[McNode], u_id: u32) {
    let (off, count) = mc_span(ma, u_id);
    for k in off..off + count {
        let qn = ma_x(ma.ma[k].x);
        let tn = ma_y(ma.ma[k].x);
        if (s[qn as usize] as i32) * (s[tn as usize] as i32) != -1 {
            continue;
        }
        let score = ma.ma[k].w as i64;
        let q_index = get_specific_hap_overlap(&ha[qn as usize], qn, tn)
            .filter(|&index| ha[qn as usize][index].score == score);
        let target = match q_index {
            Some(index) => Some((qn, index)),
            None => get_specific_hap_overlap(&ha[tn as usize], tn, qn)
                .filter(|&index| ha[tn as usize][index].score == score)
                .map(|index| (tn, index)),
        };
        match target {
            Some((uid, index)) => ha[uid as usize][index].status = crate::assembly::purge_dups::FLIP,
            None => eprintln!("ERROR"),
        }
    }
}

/// `clean_ovlp_by_mc`: keeps only the overlaps that conflict with the
/// partition, mirroring each kept overlap onto its partner.
pub fn clean_ovlp_by_mc(mg: &McG, ha: &mut [Vec<HapOverlap>]) {
    let types = [1u8, 0, 3, 2];
    for i in 0..mg.e.n_seq {
        filter_ovlp_by_mc(ha, &mg.e, &mg.s, i);
    }

    let num = ha.len();
    for v in 0..num {
        for i in 0..ha[v].len() {
            let qn = ha[v][i].x_uid;
            let tn = ha[v][i].y_uid;
            if ha[v][i].status != crate::assembly::purge_dups::FLIP {
                continue;
            }
            let x = ha[v][i];
            if let Some(index) = get_specific_hap_overlap(&ha[tn as usize], tn, qn) {
                let mut y = ha[tn as usize][index];
                set_reverse_hap_overlap(&mut y, &x, &types);
                y.status = crate::assembly::purge_dups::FLIP;
                ha[tn as usize][index] = y;
            } else {
                eprintln!("ERROR");
            }
        }
    }

    for list in ha.iter_mut() {
        list.retain(|overlap| overlap.status == crate::assembly::purge_dups::FLIP);
    }
}

/// `set_p_flag`: records the partition result on every read of one unitig.
fn set_p_flag(unitig: &MaUtg, het: &[u8], trio_flags: &mut [u8], spin: McNode) {
    for &packed in &unitig.reads {
        let read = (packed >> 33) as usize;
        trio_flags[read] |= SET_TRIO;
        let flag = het[read];
        if flag == N_HET {
            continue;
        }
        if spin == 0 {
            if flag & P_HET != 0 {
                trio_flags[read] |= FATHER;
            }
            continue;
        }
        trio_flags[read] |= if spin > 0 { FATHER } else { MOTHER };
    }
}

/// `p_nodes`: writes the spin of every unitig onto its reads.
fn p_nodes(mg: &McG, unitig_reads: &[MaUtg], het: &[u8], trio_flags: &mut [u8]) {
    for (uid, unitig) in unitig_reads.iter().enumerate() {
        set_p_flag(unitig, het, trio_flags, mg.s[uid]);
    }
}

/// `mc_g_cc_core`/`mb_g_cc_core`: connected components of the edge graph,
/// numbered in vertex order.
fn graph_cc(n_seq: u32, edges: &[u64], idx: &[u64]) -> Vec<u64> {
    let n = n_seq as usize;
    let mut flag = vec![u32::MAX; n];
    let mut stack: Vec<u32> = Vec::new();
    for i in 0..n {
        if flag[i] != u32::MAX {
            continue;
        }
        stack.clear();
        stack.push(i as u32);
        while let Some(k) = stack.pop() {
            flag[k as usize] = i as u32;
            let packed = idx[k as usize];
            let off = (packed >> 32) as usize;
            let count = packed as u32 as usize;
            for j in 0..count {
                let t = ma_y(edges[off + j]);
                if flag[t as usize] != u32::MAX {
                    continue;
                }
                stack.push(t);
            }
        }
    }
    let mut group: Vec<u64> = (0..n).map(|i| ((flag[i] as u64) << 32) | i as u64).collect();
    radix_sort_by_key(&mut group, &|value: &u64| *value);
    let mut x = 0usize;
    let mut y = 0u64;
    for i in 1..=n {
        if i == n || group[i] >> 32 != group[x] >> 32 {
            for j in x..i {
                group[j] = (y << 32) | (group[j] as u32 as u64);
            }
            y += 1;
            x = i;
        }
    }
    group
}

fn mc_g_cc_core(ma: &McMatch) -> Vec<u64> {
    let edges: Vec<u64> = ma.ma.iter().map(|edge| edge.x).collect();
    graph_cc(ma.n_seq, &edges, &ma.idx)
}

fn mb_g_cc_core(ma: &MbMatch) -> Vec<u64> {
    let edges: Vec<u64> = ma.ma.iter().map(|edge| edge.x).collect();
    graph_cc(ma.n_seq, &edges, &ma.idx)
}

/// `u_trans_a`/`u_trans_n` restricted to one query.
fn transitions_of<'a>(ta: &'a [UTrans], idx: &[u64], qn: u32) -> &'a [UTrans] {
    let (off, count) = u_trans_span(idx, qn);
    &ta[off..off + count]
}

// ---------------------------------------------------------------------------
// Block level
// ---------------------------------------------------------------------------

/// `mb_nodes_core`: groups the unitigs into the blocks the `ref` transitions
/// alternate between, then numbers the groups.
fn mb_nodes_core(
    refs: &[UTrans],
    ref_idx: &[u64],
    ma: &McMatch,
    occ: u32,
    fg: &mut Vec<u32>,
    sk: &mut Vec<u32>,
) -> Vec<u64> {
    *fg = vec![u32::MAX; occ as usize];
    for i in 0..occ {
        if fg[i as usize] != u32::MAX {
            continue;
        }
        if mc_span(ma, i).1 == 0 {
            fg[i as usize] = i;
            continue;
        }
        sk.clear();
        sk.push(i);
        while let Some(k) = sk.pop() {
            fg[k as usize] = i;
            let entries = transitions_of(refs, ref_idx, k);
            let mut st = 0usize;
            for j in 1..=entries.len() {
                if j == entries.len() || entries[j].target_identifier != entries[st].target_identifier
                {
                    let t = entries[st].target_identifier;
                    if mc_span(ma, t).1 != 0 && fg[t as usize] == u32::MAX {
                        sk.push(t);
                    }
                    st = j;
                }
            }
        }
    }

    let mut group: Vec<u64> = (0..occ as usize)
        .map(|i| ((fg[i] as u64) << 32) | i as u64)
        .collect();
    radix_sort_by_key(&mut group, &|value: &u64| *value);
    let mut x = 0usize;
    let mut y = 0u64;
    let n = occ as usize;
    for i in 1..=n {
        if i == n || group[i] >> 32 != group[x] >> 32 {
            for j in x..i {
                group[j] = (y << 32) | (group[j] as u32 as u64);
            }
            y += 1;
            x = i;
        }
    }
    group
}

/// `assgin_mb_node`: two-colours one connected component; a conflict gives
/// every vertex its own single-strand block.
fn assgin_mb_node(
    x: &mut MbNodes,
    refs: &[UTrans],
    ref_idx: &[u64],
    ma: &McMatch,
    fg: &mut [u32],
    sk: &mut Vec<u32>,
    cc: &[u64],
    cc_off: usize,
    cc_size: usize,
) {
    for i in 0..cc_size {
        let v = cc[cc_off + i] as u32;
        fg[v as usize] = u32::MAX;
    }
    sk.clear();
    let first = cc[cc_off] as u32;
    sk.push(first);
    fg[first as usize] = 0;
    let mut pass = true;
    while let Some(v) = sk.pop() {
        if mc_span(ma, v).1 == 0 {
            continue;
        }
        let entries = transitions_of(refs, ref_idx, v);
        let mut st = 0usize;
        for j in 1..=entries.len() {
            if j == entries.len() || entries[j].target_identifier != entries[st].target_identifier
            {
                let t = entries[st].target_identifier;
                if mc_span(ma, t).1 == 0 {
                    st = j;
                    continue;
                }
                if fg[t as usize] == u32::MAX {
                    fg[t as usize] = 1 - fg[v as usize];
                    sk.push(t);
                } else if fg[t as usize] == fg[v as usize] {
                    pass = false;
                    break;
                }
                st = j;
            }
        }
        if !pass {
            break;
        }
    }

    if pass {
        let uid = x.u.len() as u32;
        let mut node = MbNode {
            a: [0, 0],
            occ: [0, 0],
            s: [0, 0],
            z: [0.0; 4],
        };
        node.a[0] = x.bid.len() as u32;
        for i in 0..cc_size {
            let v = cc[cc_off + i] as u32;
            if fg[v as usize] == 0 {
                x.bid.push(v);
                node.occ[0] += 1;
                x.idx[v as usize] = uid << 1;
            }
        }
        node.a[1] = x.bid.len() as u32;
        for i in 0..cc_size {
            let v = cc[cc_off + i] as u32;
            if fg[v as usize] == 1 {
                x.bid.push(v);
                node.occ[1] += 1;
                x.idx[v as usize] = (uid << 1) + 1;
            }
        }
        x.u.push(node);
    } else {
        for i in 0..cc_size {
            let v = cc[cc_off + i] as u32;
            let uid = x.u.len() as u32;
            let mut node = MbNode {
                a: [0, 0],
                occ: [0, 0],
                s: [0, 0],
                z: [0.0; 4],
            };
            node.a[0] = x.bid.len() as u32;
            x.bid.push(v);
            node.occ[0] = 1;
            x.idx[v as usize] = uid << 1;
            node.a[1] = x.bid.len() as u32;
            x.u.push(node);
        }
    }
}

/// `update_mb_nodes_t`.
fn update_mb_nodes_t(refs: &[UTrans], ref_idx: &[u64], ma: &McMatch) -> MbNodes {
    let occ = ma.n_seq;
    let mut x = MbNodes {
        bid: Vec::new(),
        idx: vec![u32::MAX; occ as usize],
        u: Vec::new(),
    };
    let mut sk: Vec<u32> = Vec::new();
    let mut fg: Vec<u32> = Vec::new();
    let cc = mb_nodes_core(refs, ref_idx, ma, occ, &mut fg, &mut sk);
    let n = occ as usize;
    let mut st = 0usize;
    for i in 1..=n {
        if i == n || cc[st] >> 32 != cc[i] >> 32 {
            assgin_mb_node(&mut x, refs, ref_idx, ma, &mut fg, &mut sk, &cc, st, i - st);
            st = i;
        }
    }
    x
}

/// `init_mb_g_t`: the block nodes and the block-to-block weights.
fn init_mb_g_t(e: &McMatch, refs: &[UTrans], ref_idx: &[u64], is_sys: bool) -> MbG {
    let nodes = update_mb_nodes_t(refs, ref_idx, e);
    let mut edges: Vec<MbEdge> = Vec::new();
    for qb in 0..nodes.u.len() {
        for strand in 0..2 {
            let node = nodes.u[qb];
            let a = &nodes.bid[node.a[strand] as usize
                ..node.a[strand] as usize + node.occ[strand] as usize];
            for &qn in a {
                let (off, count) = mc_span(e, qn);
                for m in off..off + count {
                    let tn = ma_y(e.ma[m].x);
                    let tb = (nodes.idx[tn as usize] >> 1) as usize;
                    if tb == qb {
                        continue;
                    }
                    let mut edge = MbEdge {
                        x: ((qb as u64) << 32) | tb as u64,
                        w: [0.0; 4],
                    };
                    edge.w[(nodes.idx[tn as usize] & 1) as usize + strand * 2] = e.ma[m].w;
                    edges.push(edge);
                }
            }
        }
    }
    radix_sort_by_key(&mut edges, &|edge: &MbEdge| edge.x);
    mb_merge_dup(&mut edges);
    let mut match_graph = MbMatch {
        idx: Vec::new(),
        ma: edges,
        cc: Vec::new(),
        n_seq: nodes.u.len() as u32,
    };
    mb_edges_idx(&mut match_graph.idx, &match_graph.ma, match_graph.n_seq);
    if is_sys {
        mb_edges_symm(&mut match_graph);
    }
    MbG {
        u: nodes,
        e: match_graph,
    }
}

// ---------------------------------------------------------------------------
// Scores, spins and the local search
// ---------------------------------------------------------------------------

/// `mc_svaux_t`.
#[derive(Debug, Clone, Default)]
struct McSvaux {
    x: u64,
    cc_off: u32,
    cc_size: u32,
    cc_node: Vec<u32>,
    cc_edge: Vec<u64>,
    bfs: Vec<u32>,
    bfs_mark: Vec<u32>,
    z: Vec<[f64; 2]>,
    z_opt: Vec<[f64; 2]>,
    s: Vec<McNode>,
    s_opt: Vec<McNode>,
    f: Vec<u8>,
}

/// `mb_svaux_t`.
#[derive(Debug, Clone, Default)]
struct MbSvaux {
    x: u64,
    cc_off: u32,
    cc_size: u32,
    cc_node: Vec<u32>,
    bfs: Vec<u32>,
    bfs_mark: Vec<u32>,
    u_opt: Vec<MbNode>,
}

/// `mc_svaux_init`: the scratch arrays sized for the largest component.
///
/// The reference aliases `b->s` with `mg->s.a`; the port keeps a copy and
/// copies it back at the phase boundaries (see [`mc_solve`]).
fn mc_svaux_init(mg: &McG, x: u64) -> McSvaux {
    let ma = &mg.e;
    let n = ma.n_seq as usize;
    let mut max_cc = 0usize;
    let mut st = 0usize;
    for i in 1..=n {
        if i == n || ma.cc[st] >> 32 != ma.cc[i] >> 32 {
            max_cc = max_cc.max(i - st);
            st = i;
        }
    }
    McSvaux {
        x,
        cc_off: 0,
        cc_size: 0,
        cc_node: vec![0u32; max_cc],
        cc_edge: Vec::new(),
        bfs: vec![0u32; n],
        bfs_mark: vec![u32::MAX; n],
        z: vec![[0.0; 2]; n],
        z_opt: vec![[0.0; 2]; n],
        s: mg.s.clone(),
        s_opt: vec![0; n],
        f: vec![0u8; n],
    }
}

/// `mb_svaux_init`.
fn mb_svaux_init(mg: &MbG, x: u64) -> MbSvaux {
    let ma = &mg.e;
    let n = ma.n_seq as usize;
    let mut max_cc = 0usize;
    let mut st = 0usize;
    for i in 1..=n {
        if i == n || ma.cc[st] >> 32 != ma.cc[i] >> 32 {
            max_cc = max_cc.max(i - st);
            st = i;
        }
    }
    MbSvaux {
        x,
        cc_off: 0,
        cc_size: 0,
        cc_node: vec![0u32; max_cc],
        bfs: vec![0u32; n],
        bfs_mark: vec![u32::MAX; n],
        u_opt: vec![MbNode::default(); n],
    }
}

/// `mc_best`: the node whose flip would help most.
///
/// The reference's update is `w = max_w; max_i = k;`, so `max_w` never
/// changes and the *last* candidate wins. The port keeps that behaviour.
fn mc_best(ma: &McMatch, b: &McSvaux) -> u32 {
    let mut max_i = u32::MAX;
    let max_w = -1.0f64;
    for i in 0..b.cc_size {
        let k = ma.cc[b.cc_off as usize + i as usize] as u32;
        if b.f[k as usize] != 0 || b.s[k as usize] == 0 {
            continue;
        }
        let w = b.s[k as usize] as f64
            * (b.z[k as usize][0] - b.z[k as usize][1])
            * 2.0;
        if w <= 0.0 {
            continue;
        }
        if w > max_w {
            max_i = k;
        }
    }
    max_i
}

/// `mb_score`.
fn mb_score(mbg: &MbG, b: &MbSvaux) -> f64 {
    let ma = &mbg.e;
    let mut z = 0.0f64;
    for i in 0..b.cc_size {
        let k = ma.cc[b.cc_off as usize + i as usize] as u32 as usize;
        let node = mbg.u.u[k];
        z += -(node.s[0] as f64) * (node.z[0] - node.z[1]);
        z += -(node.s[1] as f64) * (node.z[2] - node.z[3]);
    }
    z
}

/// `mc_score`.
fn mc_score(ma: &McMatch, b: &McSvaux) -> f64 {
    let mut z = 0.0f64;
    for i in 0..b.cc_size {
        let k = ma.cc[b.cc_off as usize + i as usize] as u32 as usize;
        z += -(b.s[k] as f64) * (b.z[k][0] - b.z[k][1]);
    }
    z
}

/// `mc_reset_z`.
fn mc_reset_z(ma: &McMatch, b: &mut McSvaux) {
    for i in 0..b.cc_size {
        let k = ma.cc[b.cc_off as usize + i as usize] as u32;
        let (o, n) = mc_span(ma, k);
        b.z[k as usize] = [0.0, 0.0];
        for j in 0..n {
            let edge = ma.ma[o + j];
            let t = ma_y(edge.x) as usize;
            if b.s[t] > 0 {
                b.z[k as usize][0] += edge.w;
            } else if b.s[t] < 0 {
                b.z[k as usize][1] += edge.w;
            }
        }
    }
}

/// `mb_reset_z`: the four-weight counterpart of [`mc_reset_z`].
fn reset_mb_g_t_z(mbg: &mut MbG) {
    let ma = &mbg.e;
    for k in 0..ma.n_seq as usize {
        let (o, n) = mb_span(&mbg.e, k as u32);
        mbg.u.u[k].z = [0.0; 4];
        for j in 0..n {
            let edge = ma.ma[o + j];
            let t = ma_y(edge.x) as usize;
            if mbg.u.u[t].s[0] > 0 {
                mbg.u.u[k].z[0] += edge.w[0];
            } else if mbg.u.u[t].s[0] < 0 {
                mbg.u.u[k].z[1] += edge.w[0];
            }
            if mbg.u.u[t].s[1] > 0 {
                mbg.u.u[k].z[0] += edge.w[1];
            } else if mbg.u.u[t].s[1] < 0 {
                mbg.u.u[k].z[1] += edge.w[1];
            }
            if mbg.u.u[t].s[0] > 0 {
                mbg.u.u[k].z[2] += edge.w[2];
            } else if mbg.u.u[t].s[0] < 0 {
                mbg.u.u[k].z[3] += edge.w[2];
            }
            if mbg.u.u[t].s[1] > 0 {
                mbg.u.u[k].z[2] += edge.w[3];
            } else if mbg.u.u[t].s[1] < 0 {
                mbg.u.u[k].z[3] += edge.w[3];
            }
        }
    }
}

/// `ks_shuffle_uint32_t`: a Fisher-Yates shuffle driven by the local RNG.
fn ks_shuffle_uint32_t(n: usize, a: &mut [u32], x: &mut u64) {
    let mut i = n;
    while i > 1 {
        let j = (kr_drand_r(x) * i as f64) as usize;
        a.swap(j, i - 1);
        i -= 1;
    }
}

/// `mc_init_spin`: the greedy initial spins of one component.
fn mc_init_spin(ma: &McMatch, b: &mut McSvaux) -> f64 {
    b.cc_edge.clear();
    if b.cc_size <= 2 {
        for i in 0..b.cc_size {
            let k = ma.cc[b.cc_off as usize + i as usize] as u32 as usize;
            b.s[k] = 0;
        }
    }

    let mut first = 0usize;
    while first < b.cc_size as usize {
        let k = ma.cc[b.cc_off as usize + first] as u32;
        b.cc_node[first] = k;
        if b.s[k as usize] == 0 {
            break;
        }
        first += 1;
    }
    if first >= b.cc_size as usize {
        mc_reset_z(ma, b);
        return mc_score(ma, b);
    }

    for i in 0..b.cc_size as usize {
        let k = ma.cc[b.cc_off as usize + i] as u32;
        let (o, n) = mc_span(ma, k);
        b.cc_node[i] = k;
        for j in 0..n {
            let w = ma.ma[o + j].w.abs();
            let w = w as u32;
            b.cc_edge
                .push(((u32::MAX - w) as u64) << 32 | (o + j) as u64);
        }
    }
    radix_sort_by_key(&mut b.cc_edge, &|entry: &u64| *entry);
    for i in 0..b.cc_edge.len() {
        let entry = b.cc_edge[i];
        let edge = ma.ma[entry as u32 as usize];
        let n1 = ma_x(edge.x) as usize;
        let n2 = ma_y(edge.x) as usize;
        if b.s[n1] == 0 && b.s[n2] == 0 {
            b.x = kr_splitmix64(b.x);
            b.s[n1] = if b.x & 1 != 0 { 1 } else { -1 };
            b.s[n2] = if edge.w > 0.0 { -b.s[n1] } else { b.s[n1] };
        } else if b.s[n1] == 0 {
            b.s[n1] = if edge.w > 0.0 { -b.s[n2] } else { b.s[n2] };
        } else if b.s[n2] == 0 {
            b.s[n2] = if edge.w > 0.0 { -b.s[n1] } else { b.s[n1] };
        }
    }

    mc_reset_z(ma, b);
    mc_score(ma, b)
}

/// `mc_set_spin`: moves the weight of `k`'s edges between the two score
/// buckets of every neighbour.
fn mc_set_spin(ma: &McMatch, b: &mut McSvaux, k: u32, s: McNode) {
    let s0 = b.s[k as usize];
    if s0 == s {
        return;
    }
    let (o, n) = mc_span(ma, k);
    for j in 0..n {
        let edge = ma.ma[o + j];
        let t = ma_y(edge.x) as usize;
        if s0 != 0 {
            b.z[t][(s0 < 0) as usize] -= edge.w;
        }
        if s != 0 {
            b.z[t][(s < 0) as usize] += edge.w;
        }
    }
    b.s[k as usize] = s;
}

/// `mb_flip_spin`: flips both spins of a block and repairs its neighbours.
fn mb_flip_spin(mbg: &mut MbG, _b: &mut MbSvaux, k: u32) {
    if mbg.u.u[k as usize].s[0] == 0 && mbg.u.u[k as usize].s[1] == 0 {
        return;
    }
    let (o, n) = mb_span(&mbg.e, k);
    mbg.u.u[k as usize].s[0] = -mbg.u.u[k as usize].s[0];
    mbg.u.u[k as usize].s[1] = -mbg.u.u[k as usize].s[1];
    let s0 = mbg.u.u[k as usize].s[0];
    let s1 = mbg.u.u[k as usize].s[1];
    for j in 0..n {
        let edge = mbg.e.ma[o + j];
        let t = ma_y(edge.x) as usize;
        if s0 != 0 {
            let l = (s0 < 0) as usize;
            mbg.u.u[t].z[l] += edge.w[0];
            mbg.u.u[t].z[1 - l] -= edge.w[0];
            mbg.u.u[t].z[2 + l] += edge.w[1];
            mbg.u.u[t].z[3 - l] -= edge.w[1];
        }
        if s1 != 0 {
            let l = (s1 < 0) as usize;
            mbg.u.u[t].z[l] += edge.w[2];
            mbg.u.u[t].z[1 - l] -= edge.w[2];
            mbg.u.u[t].z[2 + l] += edge.w[3];
            mbg.u.u[t].z[3 - l] -= edge.w[3];
        }
    }
}

/// `mc_best_flip`: flips every improving node once.
fn mc_best_flip(ma: &McMatch, b: &mut McSvaux) {
    for idx in 0..b.cc_size {
        let k = ma.cc[b.cc_off as usize + idx as usize] as u32 as usize;
        b.f[k] = 0;
    }
    loop {
        let idx = mc_best(ma, b);
        if idx == u32::MAX {
            break;
        }
        let spin = -b.s[idx as usize];
        mc_set_spin(ma, b, idx, spin);
        b.f[idx as usize] = 1;
    }
}

/// `mc_optimize_local`.
fn mc_optimize_local(opt: &McOpt, ma: &McMatch, b: &mut McSvaux, n_iter: &mut u32) -> f64 {
    let mut n_iter_local = 0i32;
    while n_iter_local < opt.max_iter {
        *n_iter += 1;
        ks_shuffle_uint32_t(b.cc_size as usize, &mut b.cc_node, &mut b.x);
        let mut n_flip = 0u32;
        for i in 0..b.cc_size as usize {
            let k = b.cc_node[i];
            let z = b.z[k as usize];
            if z[0] == z[1] {
                continue;
            }
            let s: McNode = if z[0] > z[1] { -1 } else { 1 };
            if b.s[k as usize] != s {
                mc_set_spin(ma, b, k, s);
                n_flip += 1;
            }
        }
        n_iter_local += 1;
        if n_flip == 0 {
            break;
        }
    }
    mc_score(ma, b)
}

/// `mb_optimize_local`.
fn mb_optimize_local(opt: &McOpt, mbg: &mut MbG, b: &mut MbSvaux, n_iter: &mut u32) -> f64 {
    let mut n_iter_local = 0i32;
    while n_iter_local < opt.max_iter {
        *n_iter += 1;
        ks_shuffle_uint32_t(b.cc_size as usize, &mut b.cc_node, &mut b.x);
        let mut n_flip = 0u32;
        for i in 0..b.cc_size as usize {
            let k = b.cc_node[i];
            let node = mbg.u.u[k as usize];
            if node.z[0] == node.z[1] && node.z[2] == node.z[3] {
                continue;
            }
            let z = -(node.s[0] as f64) * (node.z[0] - node.z[1])
                - (node.s[1] as f64) * (node.z[2] - node.z[3]);
            if z >= 0.0 {
                continue;
            }
            mb_flip_spin(mbg, b, k);
            n_flip += 1;
        }
        n_iter_local += 1;
        if n_flip == 0 {
            break;
        }
    }
    mb_score(mbg, b)
}

/// `mc_perturb`.
fn mc_perturb(opt: &McOpt, ma: &McMatch, b: &mut McSvaux) {
    for i in 0..b.cc_size {
        let k = ma.cc[b.cc_off as usize + i as usize] as u32;
        let y = kr_drand_r(&mut b.x);
        if y < opt.f_perturb {
            let spin = -b.s[k as usize];
            mc_set_spin(ma, b, k, spin);
        }
    }
}

/// `mb_perturb`.
fn mb_perturb(opt: &McOpt, mbg: &mut MbG, b: &mut MbSvaux) {
    for i in 0..b.cc_size {
        let k = mbg.e.cc[b.cc_off as usize + i as usize] as u32;
        let y = kr_drand_r(&mut b.x);
        if y < opt.f_perturb {
            mb_flip_spin(mbg, b, k);
        }
    }
}

/// `mc_bfs`.
fn mc_bfs(ma: &McMatch, b: &mut McSvaux, k0: u32, bfs_round: u32, max_size: u32) -> u32 {
    b.bfs[0] = k0;
    b.bfs_mark[k0 as usize] = k0;
    let mut n_bfs = 1usize;
    let mut st = 0usize;
    let mut en = n_bfs;
    for _ in 0..bfs_round {
        for i in st..en {
            let k = b.bfs[i];
            let (o, n) = mc_span(ma, k);
            for j in 0..n {
                let t = ma_y(ma.ma[o + j].x);
                if b.bfs_mark[t as usize] != k0 {
                    b.bfs[n_bfs] = t;
                    n_bfs += 1;
                    b.bfs_mark[t as usize] = k0;
                }
            }
        }
        st = en;
        en = n_bfs;
        if max_size > 0 && n_bfs as u32 > max_size {
            break;
        }
    }
    n_bfs as u32
}

/// `mc_perturb_node`: flips the BFS neighbourhood of a random node.
fn mc_perturb_node(opt: &McOpt, ma: &McMatch, b: &mut McSvaux, bfs_round: i32) {
    let mut k = (kr_drand_r(&mut b.x) * b.cc_size as f64 + 0.499) as u32;
    if k >= b.cc_size {
        k = b.cc_size - 1;
    }
    let k = ma.cc[b.cc_off as usize + k as usize] as u32;
    let max_size = (b.cc_size as f64 * opt.f_perturb) as i32 as u32;
    let n_bfs = mc_bfs(ma, b, k, bfs_round as u32, max_size);
    for i in 0..n_bfs as usize {
        let node = b.bfs[i];
        let spin = -b.s[node as usize];
        mc_set_spin(ma, b, node, spin);
    }
}

/// `mb_bfs`.
fn mb_bfs(ma: &MbMatch, b: &mut MbSvaux, k0: u32, bfs_round: u32, max_size: u32) -> u32 {
    b.bfs[0] = k0;
    b.bfs_mark[k0 as usize] = k0;
    let mut n_bfs = 1usize;
    let mut st = 0usize;
    let mut en = n_bfs;
    for _ in 0..bfs_round {
        for i in st..en {
            let k = b.bfs[i];
            let (o, n) = mb_span(ma, k);
            for j in 0..n {
                let t = ma_y(ma.ma[o + j].x);
                if b.bfs_mark[t as usize] != k0 {
                    b.bfs[n_bfs] = t;
                    n_bfs += 1;
                    b.bfs_mark[t as usize] = k0;
                }
            }
        }
        st = en;
        en = n_bfs;
        if max_size > 0 && n_bfs as u32 > max_size {
            break;
        }
    }
    n_bfs as u32
}

/// `mb_perturb_node`.
fn mb_perturb_node(opt: &McOpt, mbg: &mut MbG, b: &mut MbSvaux, bfs_round: i32) {
    let mut k = (kr_drand_r(&mut b.x) * b.cc_size as f64 + 0.499) as u32;
    if k >= b.cc_size {
        k = b.cc_size - 1;
    }
    let k = mbg.e.cc[b.cc_off as usize + k as usize] as u32;
    let max_size = (b.cc_size as f64 * opt.f_perturb) as i32 as u32;
    let n_bfs = mb_bfs(&mbg.e, b, k, bfs_round as u32, max_size);
    for i in 0..n_bfs as usize {
        let node = b.bfs[i];
        mb_flip_spin(mbg, b, node);
    }
}

/// `mb_solve_cc`: the block-level local search of one component.
fn mb_solve_cc(opt: &McOpt, mbg: &mut MbG, b: &mut MbSvaux, cc_off: u32, cc_size: u32) -> u32 {
    b.cc_off = cc_off;
    b.cc_size = cc_size;
    if cc_size <= 2 {
        return 0;
    }
    for j in 0..cc_size {
        let k = mbg.e.cc[(cc_off + j) as usize] as u32;
        b.cc_node[j as usize] = k;
    }
    let mut sc_opt = mb_score(mbg, b);
    for j in 0..cc_size {
        let k = b.cc_node[j as usize] as usize;
        b.u_opt[k] = mbg.u.u[k];
    }
    let mut n_iter = 0u32;
    let mut sc = mb_optimize_local(opt, mbg, b, &mut n_iter);
    if sc > sc_opt {
        for j in 0..cc_size {
            let k = b.cc_node[j as usize] as usize;
            b.u_opt[k] = mbg.u.u[k];
        }
        sc_opt = sc;
    } else {
        for j in 0..cc_size {
            let k = b.cc_node[j as usize] as usize;
            mbg.u.u[k] = b.u_opt[k];
        }
    }

    let flush = (opt.max_iter * 50) as u32;
    for k in 0..opt.n_perturb as u32 {
        if k & 1 != 0 {
            mb_perturb(opt, mbg, b);
        } else {
            mb_perturb_node(opt, mbg, b, 3);
        }
        sc = mb_optimize_local(opt, mbg, b, &mut n_iter);
        if sc > sc_opt {
            for j in 0..cc_size {
                let node = b.cc_node[j as usize] as usize;
                b.u_opt[node] = mbg.u.u[node];
            }
            sc_opt = sc;
        } else {
            for j in 0..cc_size {
                let node = b.cc_node[j as usize] as usize;
                mbg.u.u[node] = b.u_opt[node];
            }
        }
        if n_iter % flush == 0 {
            reset_mb_g_t_z(mbg);
            sc = mb_score(mbg, b);
            for j in 0..cc_size {
                let node = b.cc_node[j as usize] as usize;
                b.u_opt[node].z = mbg.u.u[node].z;
            }
            sc_opt = sc;
        }
    }
    for j in 0..cc_size {
        let node = b.cc_node[j as usize] as usize;
        mbg.u.u[node] = b.u_opt[node];
    }
    n_iter
}

/// `incre_weight`: the marginal gain of adding `uid` to the current cluster.
fn incre_weight(b_aux: &McSvaux, ma: &McMatch, vis: &[u8], uid: u32) -> f64 {
    let w = b_aux.s[uid as usize] as f64
        * (b_aux.z[uid as usize][0] - b_aux.z[uid as usize][1])
        * 2.0;
    let mut w_off = 0.0f64;
    let (o, n) = mc_span(ma, uid);
    for j in 0..n {
        let t = ma_y(ma.ma[o + j].x);
        if vis[t as usize] == 0 {
            continue;
        }
        if t == uid {
            continue;
        }
        w_off += b_aux.s[uid as usize] as f64
            * b_aux.s[t as usize] as f64
            * ma.ma[o + j].w;
    }
    w - w_off * 4.0
}

/// `mc_init_spin_all`: greedy spins for every component, then an adjustment
/// from the block structure.
fn mc_init_spin_all(_opt: &McOpt, mg: &mut McG, mbg: Option<&mut MbG>, b: &mut McSvaux) {
    let n = mg.e.n_seq as usize;
    let mut st = 0usize;
    for i in 1..=n {
        if i == n || mg.e.cc[st] >> 32 != mg.e.cc[i] >> 32 {
            b.cc_off = st as u32;
            b.cc_size = (i - st) as u32;
            if b.cc_size >= 2 {
                mc_init_spin(&mg.e, b);
            }
            st = i;
        }
    }

    let Some(mbg) = mbg else {
        return;
    };
    let mut s: Vec<u32> = Vec::new();
    let mut vis = vec![0u8; n];
    let mut max_w = 0.0f64;
    for i in 0..mbg.u.u.len() {
        let node = mbg.u.u[i];
        let a0 = &mbg.u.bid[node.a[0] as usize..node.a[0] as usize + node.occ[0] as usize];
        let a1 = &mbg.u.bid[node.a[1] as usize..node.a[1] as usize + node.occ[1] as usize];
        if a0.len() + a1.len() <= 1 {
            continue;
        }
        s.clear();
        for &qn in a0 {
            if b.s[qn as usize] == 0 {
                continue;
            }
            let mut qn = qn;
            if b.s[qn as usize] < 0 {
                qn += 1 << 31;
            }
            s.push(qn);
        }
        for &qn in a1 {
            if b.s[qn as usize] == 0 {
                continue;
            }
            let mut qn = qn;
            if b.s[qn as usize] > 0 {
                qn += 1 << 31;
            }
            s.push(qn);
        }
        if s.is_empty() {
            continue;
        }
        radix_sort_by_key_32(&mut s, &|value: &u32| *value);

        let mut ws = u32::MAX;
        let mut we = u32::MAX;
        let mut n_flip = 0u32;
        let mut group_start = 0usize;
        for k in 1..=s.len() {
            if k == s.len() || (s[k] >> 31) != (s[group_start] >> 31) {
                let mut w = 0.0f64;
                for m in group_start..k {
                    vis[(s[m] & 0x7fff_ffff) as usize] = 1;
                }
                for m in group_start..k {
                    w += incre_weight(b, &mg.e, &vis, s[m] & 0x7fff_ffff);
                }
                for m in group_start..k {
                    vis[(s[m] & 0x7fff_ffff) as usize] = 0;
                }
                if ws == u32::MAX || max_w < w {
                    max_w = w;
                    ws = group_start as u32;
                    we = k as u32;
                }
                group_start = k;
                n_flip += 1;
            }
        }
        if n_flip <= 1 {
            continue;
        }
        for m in ws as usize..we as usize {
            let node = s[m] & 0x7fff_ffff;
            let spin = -b.s[node as usize];
            mc_set_spin(&mg.e, b, node, spin);
        }
    }

    for i in 0..mbg.u.u.len() {
        let node = mbg.u.u[i];
        let a0 = &mbg.u.bid[node.a[0] as usize..node.a[0] as usize + node.occ[0] as usize];
        let a1 = &mbg.u.bid[node.a[1] as usize..node.a[1] as usize + node.occ[1] as usize];
        mbg.u.u[i].s = [0, 0];
        mbg.u.u[i].z = [0.0; 4];
        for &qn in a0 {
            if b.s[qn as usize] == 0 {
                continue;
            }
            mbg.u.u[i].s[0] = b.s[qn as usize];
            break;
        }
        for &qn in a1 {
            if b.s[qn as usize] == 0 {
                continue;
            }
            mbg.u.u[i].s[1] = b.s[qn as usize];
            break;
        }
    }
    reset_mb_g_t_z(mbg);
}

/// `mc_set_by_mbg`: projects the block spins back onto the unitigs.
fn mc_set_by_mbg(mg: &mut McG, mbg: &MbG) {
    for i in 0..mbg.u.u.len() {
        let node = mbg.u.u[i];
        let a0 = &mbg.u.bid[node.a[0] as usize..node.a[0] as usize + node.occ[0] as usize];
        let a1 = &mbg.u.bid[node.a[1] as usize..node.a[1] as usize + node.occ[1] as usize];
        for &qn in a0 {
            if mg.s[qn as usize] == 0 {
                continue;
            }
            mg.s[qn as usize] = node.s[0];
        }
        for &qn in a1 {
            if mg.s[qn as usize] == 0 {
                continue;
            }
            mg.s[qn as usize] = node.s[1];
        }
    }
}

/// `mb_solve_core`: the block-level pass.
fn mb_solve_core(opt: &mut McOpt, mg: &mut McG, refs: &[UTrans], ref_idx: &[u64], is_sys: bool) {
    let mut mbg = init_mb_g_t(&mg.e, refs, ref_idx, is_sys);

    mg.e.cc = mc_g_cc_core(&mg.e);
    let mut b = mc_svaux_init(mg, opt.seed);
    mc_init_spin_all(opt, mg, Some(&mut mbg), &mut b);
    mg.s = b.s.clone();
    mg.e.cc.clear();

    mbg.e.cc = mb_g_cc_core(&mbg.e);
    let mut bb = mb_svaux_init(&mbg, opt.seed);

    opt.n_perturb = opt.n_b_perturb;
    let n = mbg.e.n_seq as usize;
    let mut st = 0usize;
    for i in 1..=n {
        if i == n || mbg.e.cc[st] >> 32 != mbg.e.cc[i] >> 32 {
            mb_solve_cc(opt, &mut mbg, &mut bb, st as u32, (i - st) as u32);
            st = i;
        }
    }
    opt.n_perturb = opt.n_s_perturb - opt.n_b_perturb;

    mc_set_by_mbg(mg, &mbg);
}

// ---------------------------------------------------------------------------
// The cluster layout (`layout_mc_clus_t` and its `osg_t` helpers)
// ---------------------------------------------------------------------------

/// `osg_arc_t`.
#[derive(Debug, Clone, Copy, Default)]
struct OsgArc {
    v: u32,
    u: u32,
    occ: u32,
    del: bool,
    w: f64,
    nw: f64,
}

/// `osg_t`: a string graph over the chain boxes of one connected component.
#[derive(Debug, Clone, Default)]
struct Osg {
    arcs: Vec<OsgArc>,
    /// `seq[].del`, one flag per box.
    seqs: Vec<bool>,
    idx: Vec<u64>,
    is_sorted: bool,
}

impl Osg {
    fn seq_set(&mut self, sid: usize, del: bool) {
        if sid >= self.seqs.len() {
            self.seqs.resize(sid + 1, false);
        }
        self.seqs[sid] = del;
    }

    fn arc_push(&mut self) -> &mut OsgArc {
        self.arcs.push(OsgArc::default());
        self.arcs.last_mut().unwrap()
    }

    fn span(&self, v: u32) -> (usize, usize) {
        let packed = self.idx[v as usize];
        ((packed >> 32) as usize, packed as u32 as usize)
    }

    fn arc_count(&self, v: u32) -> u32 {
        self.idx[v as usize] as u32
    }

    fn arc_first(&self, v: u32) -> OsgArc {
        self.arcs[(self.idx[v as usize] >> 32) as usize]
    }

    fn arc_rm(&mut self) {
        let n = self.arcs.len();
        let mut kept: Vec<OsgArc> = Vec::with_capacity(n);
        for arc in &self.arcs {
            let u = (arc.u >> 1) as usize;
            let v = (arc.v >> 1) as usize;
            if !arc.del && !self.seqs[u] && !self.seqs[v] {
                kept.push(*arc);
            }
        }
        if kept.len() < n {
            self.idx.clear();
        }
        self.arcs = kept;
    }

    fn arc_index(&mut self) {
        let mut index = vec![0u64; self.seqs.len() * 2];
        let n = self.arcs.len();
        let mut last = 0usize;
        for i in 1..=n {
            if i == n || self.arcs[i - 1].u != self.arcs[i].u {
                index[self.arcs[i - 1].u as usize] = ((last as u64) << 32) | (i - last) as u64;
                last = i;
            }
        }
        self.idx = index;
    }

    fn cleanup(&mut self) {
        self.arc_rm();
        if !self.is_sorted {
            radix_sort_by_key_32(&mut self.arcs, &|arc: &OsgArc| arc.u);
            self.is_sorted = true;
        }
        if self.idx.is_empty() {
            self.arc_index();
        }
    }
}

/// `cal_chain_arch`: the box-to-box arcs of one chain side.
fn cal_chain_arch(
    sg: &mut Osg,
    ma: &McMatch,
    v: u32,
    nodes: &[u32],
    idx: &[u64],
    srt: &mut Vec<u64>,
) {
    srt.clear();
    for &node in nodes {
        let packed = idx[node as usize];
        debug_assert!(v == (packed >> 32) as u32 || v == packed as u32);
        let (o, n) = mc_span(ma, node);
        for j in 0..n {
            let t = ma_y(ma.ma[o + j].x);
            let w = (idx[t as usize] >> 32) as u32;
            if w != u32::MAX && (w >> 1) > (v >> 1) {
                srt.push(((w as u64) << 32) | (o + j) as u64);
            }
            let w = idx[t as usize] as u32;
            if w != u32::MAX && (w >> 1) > (v >> 1) {
                srt.push(((w as u64) << 32) | (o + j) as u64);
            }
        }
    }
    radix_sort_by_key(srt, &|value: &u64| *value);
    let mut l = 0usize;
    for k in 1..=srt.len() {
        if k == srt.len() || (srt[l] >> 32) != (srt[k] >> 32) {
            let w = (srt[l] >> 32) as u32;
            let mut mw = 0.0f64;
            for &entry in &srt[l..k] {
                mw += ma.ma[entry as u32 as usize].w.abs();
            }
            let forward = OsgArc {
                v: w,
                u: v,
                occ: 0,
                del: false,
                w: mw,
                nw: mw,
            };
            sg.arcs.push(forward);
            let reverse = OsgArc {
                v,
                u: w,
                occ: 0,
                del: false,
                w: mw,
                nw: mw,
            };
            sg.arcs.push(reverse);
            l = k;
        }
    }
}

/// `mc_clus_cut`: drops the weakest of `sa` when it is dominated on both
/// ends.
fn mc_clus_cut(sg: &mut Osg, sa: &mut [McEdge], cut_rate: f64, force_cut: bool) {
    for k in 0..sa.len() {
        if sa[k].x == u64::MAX {
            continue;
        }
        let xi = sa[k].x as usize;
        if sg.arcs[xi].del {
            sa[k].x = u64::MAX;
            continue;
        }
        let v = sg.arcs[xi].u;
        let w = sg.arcs[xi].v;
        let nv = sg.arc_count(v);
        let nw = sg.arc_count(w);
        if nv <= 1 && nw <= 1 {
            sa[k].x = u64::MAX;
            continue;
        }
        let (wo, wn) = sg.span(w);
        let mut we = None;
        for i in 0..wn {
            if sg.arcs[wo + i].v == v {
                we = Some(wo + i);
                break;
            }
        }
        let we = we.expect("the reverse arc of a live candidate exists");
        let mm_ol = sg.arcs[xi].nw;

        let (vo, vn) = sg.span(v);
        let mut kv = 0u32;
        let mut ol_max = 0.0f64;
        for i in 0..vn {
            if sg.arcs[vo + i].del {
                continue;
            }
            kv += 1;
            if ol_max < sg.arcs[vo + i].nw {
                ol_max = sg.arcs[vo + i].nw;
            }
        }
        if kv < 1 {
            sa[k].x = u64::MAX;
            continue;
        }
        if kv >= 2 && !force_cut && mm_ol > ol_max * cut_rate {
            continue;
        }

        let (wo, wn) = sg.span(w);
        let mut kw = 0u32;
        let mut ol_max = 0.0f64;
        for i in 0..wn {
            if sg.arcs[wo + i].del {
                continue;
            }
            kw += 1;
            if ol_max < sg.arcs[wo + i].nw {
                ol_max = sg.arcs[wo + i].nw;
            }
        }
        if kw < 1 {
            sa[k].x = u64::MAX;
            continue;
        }
        if kw >= 2 && !force_cut && mm_ol > ol_max * cut_rate {
            continue;
        }
        if kv <= 1 && kw <= 1 {
            sa[k].x = u64::MAX;
            continue;
        }

        sg.arcs[xi].del = true;
        sg.arcs[we].del = true;
        sa[k].x = u64::MAX;
    }
}

/// `gen_mc_clus_backbone_layout`: rewrites the cluster node order along the
/// chain backbones.
fn gen_mc_clus_backbone_layout(sg: &Osg, res: &mut Vec<u64>, out: &mut [u32], buf: &mut [u32]) {
    let n_seq = sg.seqs.len();
    res.clear();
    res.resize(n_seq, 0);
    let rn = res.len();

    for k in 0..n_seq {
        let n0 = sg.arc_count((k as u32) << 1);
        let n1 = sg.arc_count(((k as u32) << 1) + 1);
        debug_assert!(n0 <= 1 && n1 <= 1);
        if (n0 > 0) != (n1 > 0) || (n0 == 0 && n1 == 0) {
            let mut v = if n0 > 0 {
                (k as u32) << 1
            } else {
                ((k as u32) << 1) + 1
            };
            if res[k] != 0 {
                continue;
            }
            res.push((((v ^ 1) as u64) << 32) | v as u64);
            res[(v >> 1) as usize] = 1;
            while sg.arc_count(v) > 0 {
                v = sg.arc_first(v).v ^ 1;
                res.push((((v ^ 1) as u64) << 32) | v as u64);
                res[(v >> 1) as usize] = 1;
            }
        }
    }

    for k in 0..n_seq {
        if res[k] != 0 {
            continue;
        }
        let start = (k as u32) << 1;
        debug_assert_eq!(sg.arc_count(start), 1);
        debug_assert_eq!(sg.arc_count(start + 1), 1);
        let mut v = start;
        let mut weakest: Option<OsgArc> = None;
        while sg.arc_count(v) > 0 {
            let first = sg.arc_first(v);
            if weakest.map_or(true, |t| t.nw > first.nw) {
                weakest = Some(first);
            }
            v = first.v ^ 1;
            if v == start {
                break;
            }
        }
        let mut v = weakest.expect("a circle has an arc").v ^ 1;
        res.push((((v ^ 1) as u64) << 32) | v as u64);
        res[(v >> 1) as usize] = 1;
        loop {
            v = sg.arc_first(v).v ^ 1;
            if res[(v >> 1) as usize] != 0 {
                break;
            }
            res.push((((v ^ 1) as u64) << 32) | v as u64);
            res[(v >> 1) as usize] = 1;
        }
    }
    debug_assert_eq!(res.len() - rn, n_seq);

    let out_n = out.len();
    let mut i = 0usize;
    let mut l = 0usize;
    for k in 1..=out_n {
        if k == out_n || out[k] == u32::MAX {
            if k > l {
                res[i] = ((l as u64) << 32) | k as u64;
                i += 1;
            } else {
                debug_assert_eq!(k, out_n);
            }
            l = k + 1;
        }
    }
    debug_assert_eq!(i, n_seq);

    let mut z = 0usize;
    for k in rn..res.len() {
        let v = res[k] >> 32;
        let packed = res[(v >> 1) as usize];
        let s = (packed >> 32) as u32;
        let e = packed as u32;
        if v & 1 == 0 {
            for m in s..e {
                buf[z] = out[m as usize];
                z += 1;
            }
        } else {
            let mut m = e;
            while m > s {
                m -= 1;
                buf[z] = out[m as usize];
                z += 1;
            }
        }
        buf[z] = u32::MAX;
        z += 1;
    }
    debug_assert_eq!(z, out_n);
    out.copy_from_slice(&buf[..out_n]);
}

/// `layout_mc_clus_t`: orders the boxes and their chains, then rewrites the
/// cluster node list.
fn layout_mc_clus_t(
    ma: &McMatch,
    a: &mut [u32],
    sg: &mut Osg,
    buf: &mut [u32],
    idx: &mut [u64],
    unitig_graph: &UnitigGraph,
) {
    *sg = Osg::default();
    idx.fill(u64::MAX);

    let an = a.len();
    let mut l = 0usize;
    let mut i = 0u32;
    for k in 1..=an {
        if k == an || a[k] == u32::MAX {
            if k > l {
                sg.seq_set(i as usize, false);
                let mut len = 0u64;
                for z in l..k {
                    len += u64::from(unitig_graph.graph.sequences[a[z] as usize].length);
                }
                let cutoff = len >> 1;
                let mut offset = 0u64;
                for z in l..k {
                    let s = offset;
                    offset += u64::from(unitig_graph.graph.sequences[a[z] as usize].length);
                    let e = offset;
                    let node = a[z] as usize;
                    if s <= cutoff {
                        idx[node] <<= 32;
                        idx[node] |= (i as u64) << 1;
                    }
                    if e >= cutoff {
                        idx[node] <<= 32;
                        idx[node] |= ((i as u64) << 1) + 1;
                    }
                }
                i += 1;
            } else {
                debug_assert_eq!(k, an);
            }
            l = k + 1;
        }
    }

    let mut srt: Vec<u64> = Vec::new();
    l = 0;
    i = 0;
    for k in 1..=an {
        if k == an || a[k] == u32::MAX {
            if k > l {
                for side in 0..2u32 {
                    let v = (i << 1) + side;
                    let mut nodes: Vec<u32> = Vec::new();
                    for z in l..k {
                        let node = a[z];
                        let packed = idx[node as usize];
                        if v == (packed >> 32) as u32 || v == packed as u32 {
                            nodes.push(node);
                        }
                    }
                    cal_chain_arch(sg, ma, v, &nodes, idx, &mut srt);
                }
                i += 1;
            } else {
                debug_assert_eq!(k, an);
            }
            l = k + 1;
        }
    }
    sg.cleanup();

    let mut sm: Vec<McEdge> = Vec::new();
    for k in 0..sg.arcs.len() {
        let arc = sg.arcs[k];
        if sg.arc_count(arc.u) <= 1 && sg.arc_count(arc.v) <= 1 {
            continue;
        }
        sm.push(McEdge {
            x: k as u64,
            w: arc.nw,
        });
    }
    sm.sort_by(|x, y| x.w.partial_cmp(&y.w).unwrap_or(std::cmp::Ordering::Equal));

    let min_cut = 0.199999f64;
    let max_cut = 0.800001f64;
    let cut_round = 4usize;
    let step = if cut_round == 1 {
        max_cut
    } else {
        (max_cut - min_cut) / (cut_round as f64 - 1.0)
    };
    let mut drop = min_cut;
    for _ in 0..cut_round {
        if drop > max_cut {
            drop = max_cut;
        }
        mc_clus_cut(sg, &mut sm, drop, false);
        drop += step;
    }
    mc_clus_cut(sg, &mut sm, 1.1, true);
    sg.cleanup();

    gen_mc_clus_backbone_layout(sg, &mut srt, a, buf);
}

// ---------------------------------------------------------------------------
// Cluster-level optimisation (`mc_solve_cc_adv`)
// ---------------------------------------------------------------------------

/// `bits_p`.
#[derive(Debug, Clone, Default)]
struct BitsP {
    a: Vec<u8>,
    n: usize,
}

fn is_set_bits_p(bits: &BitsP, i: u32) -> bool {
    (bits.a[(i >> 3) as usize] >> (i & 7)) & 1 != 0
}

fn set_bits_p(bits: &mut BitsP, i: u32) {
    bits.a[(i >> 3) as usize] |= 1u8 << (i & 7);
}

/// `clus_flip_aux`.
#[derive(Debug, Clone, Default)]
struct ClusFlipAux {
    vis: BitsP,
    w: f64,
    off: u32,
    occ: u32,
}

/// `nn_clus_t`.
#[derive(Debug, Clone, Default)]
struct NnClus {
    nn: Vec<u32>,
    ng: Vec<u64>,
}

/// `mc_clus_t`.
struct McClus<'a> {
    bub: &'a BubbleType,
    mg_e: &'a McMatch,
    lock: Vec<u8>,
    lock_max: u8,
    n: u32,
    n_thread: usize,
    aux: Vec<ClusFlipAux>,
    asn: Option<Vec<u64>>,
    sg: Osg,
    cc: NnClus,
}

/// `init_mc_clus_t`.
fn init_mc_clus_t<'a>(
    mg_e: &'a McMatch,
    bub: &'a BubbleType,
    n_thread: u64,
    n_seq: u32,
    flp_max: u32,
) -> McClus<'a> {
    let n_thread = n_thread.min(64) as usize;
    let ss = (n_seq as usize >> 3) + usize::from(n_seq & 7 != 0);
    let mut aux = Vec::with_capacity(n_thread);
    for _ in 0..n_thread {
        aux.push(ClusFlipAux {
            vis: BitsP {
                a: vec![0u8; ss],
                n: ss,
            },
            w: 0.0,
            off: 0,
            occ: 0,
        });
    }
    let lock_max = if flp_max <= 255 { flp_max as u8 } else { 255 };
    McClus {
        bub,
        mg_e,
        lock: vec![0u8; n_seq as usize],
        lock_max: lock_max.max(1),
        n: n_seq,
        n_thread,
        aux,
        asn: None,
        sg: Osg::default(),
        cc: NnClus::default(),
    }
}

/// `clean_clus_flip_aux`.
fn clean_clus_flip_aux(z: &mut ClusFlipAux) {
    z.w = -1.0;
    z.occ = u32::MAX;
    z.off = u32::MAX;
    z.vis.a.fill(0);
}

/// `reorder_bub`: orders one bubble chain by the strongest outside links.
#[allow(clippy::too_many_arguments)]
fn reorder_bub(
    ma: &McMatch,
    a: &mut [u32],
    ff: &mut [u8],
    lf: u8,
    rf: u8,
    cf: u8,
    cuf: u8,
    sc_l: &mut [f64],
    sc_r: &mut [f64],
    sc_m: &mut [f64],
    res: &mut [u32],
) {
    let a_n = a.len();
    for k in 0..a_n {
        ff[a[k] as usize] = cf;
    }
    let mut m_inner = 0.0f64;
    let mut m_inner_k = -1i64;
    let mut mm_l = 0.0f64;
    let mut mm_r = 0.0f64;
    let mut mmlk = -1i64;
    let mut mmrk = -1i64;
    let rev;
    for k in 0..a_n {
        sc_l[k] = 0.0;
        sc_r[k] = 0.0;
        sc_m[k] = 0.0;
        let (o, n) = mc_span(ma, a[k]);
        for j in 0..n {
            let edge = ma.ma[o + j];
            let t = ma_y(edge.x) as usize;
            if ff[t] == lf {
                sc_l[k] += edge.w.abs();
            }
            if ff[t] == rf {
                sc_r[k] += edge.w.abs();
            }
            if ff[t] == cf {
                sc_m[k] += edge.w.abs();
            }
        }
        if sc_l[k] > 0.0 {
            if sc_l[k] > mm_l {
                mm_l = sc_l[k];
                mmlk = k as i64;
            }
        } else if sc_r[k] > 0.0 {
            if sc_r[k] > mm_r {
                mm_r = sc_r[k];
                mmrk = k as i64;
            }
        } else if sc_m[k] > 0.0 && sc_m[k] > m_inner {
            m_inner = sc_m[k];
            m_inner_k = k as i64;
        }
    }
    if mmlk == -1 && mmrk == -1 && m_inner_k == -1 {
        return;
    }
    let mmk = if mmlk != -1 {
        rev = 0i64;
        mmlk
    } else if mmrk != -1 {
        rev = 1i64;
        mmrk
    } else {
        rev = 0i64;
        m_inner_k
    };

    let mut z = 0usize;
    res[z] = a[mmk as usize];
    z += 1;
    ff[a[mmk as usize] as usize] = cuf;
    while z < a_n {
        let mut mm = 0.0f64;
        let mut best = -1i64;
        for k in 0..a_n {
            if ff[a[k] as usize] == cuf {
                continue;
            }
            let (o, n) = mc_span(ma, a[k]);
            for j in 0..n {
                let edge = ma.ma[o + j];
                let t = ma_y(edge.x) as usize;
                if ff[t] == cuf {
                    if rev == 0 {
                        sc_l[k] += edge.w.abs();
                    } else {
                        sc_r[k] += edge.w.abs();
                    }
                }
            }
            let score = if rev == 0 { sc_l[k] } else { sc_r[k] };
            if score > 0.0 && score > mm {
                mm = score;
                best = k as i64;
            }
        }
        if best == -1 {
            break;
        }
        res[z] = a[best as usize];
        z += 1;
        ff[a[best as usize] as usize] = cuf;
    }
    if z < a_n {
        for k in 0..a_n {
            if ff[a[k] as usize] == cuf {
                continue;
            }
            res[z] = a[k];
            z += 1;
            ff[a[k] as usize] = cuf;
        }
    }
    debug_assert_eq!(z, a_n);
    if rev == 0 {
        for k in 0..a_n {
            a[k] = res[k];
            ff[a[k] as usize] = lf;
        }
    } else {
        for k in 0..a_n {
            a[k] = res[a_n - k - 1];
            ff[a[k] as usize] = lf;
        }
    }
}

/// The `clus_a`/`clus_n` accessors: `(start, count)` of chain `id`.
fn clus_span(bc: &McClus, id: usize) -> (usize, usize) {
    let packed = bc.cc.ng[id];
    ((packed >> 32) as usize, packed as u32 as usize)
}

/// `renew_mc_clus_t_adv`: rebuilds the bubble chains of one component, with
/// the reference-cluster nodes from `asn`.
fn renew_mc_clus_t_adv(bc: &mut McClus, a: &[u32], a_n: usize, unitig_graph: &UnitigGraph) {
    bc.cc.ng.clear();
    bc.cc.nn.clear();
    bc.lock.fill(0);
    for &node in &a[..a_n] {
        bc.lock[node as usize] = 1;
    }
    bc.cc.nn.reserve(a_n);

    let empty: &[MaUtg] = &[];
    let unitigs = bc.bub.bubble_unitigs.as_deref().unwrap_or(empty);
    let mut sc_l: Vec<f64> = Vec::new();
    let mut sc_r: Vec<f64> = Vec::new();
    let mut sc_m: Vec<f64> = Vec::new();
    let mut tmp: Vec<u32> = Vec::new();
    let mut bid = 0usize;
    let mut bk = 0usize;
    for unitig in unitigs {
        if unitig.reads.is_empty() {
            continue;
        }
        let mut cocc = 0usize;
        let iin = bc.cc.nn.len();
        for k in 0..unitig.reads.len() {
            let (_, _, body, _) = get_bubbles(bc.bub, unitig.reads[k] >> 33);
            let start = bc.cc.nn.len() as u64;
            bc.cc.ng.push(start);
            for &node in body {
                if bc.lock[(node >> 1) as usize] == 0 {
                    continue;
                }
                bc.cc.nn.push(node >> 1);
                bc.lock[(node >> 1) as usize] = 2;
            }
            if let Some(asn) = &bc.asn {
                while bk < asn.len() && (asn[bk] >> 32) < bid as u64 {
                    bk += 1;
                }
                while bk < asn.len() && (asn[bk] >> 32) == bid as u64 {
                    let node = asn[bk] as u32;
                    if bc.lock[node as usize] != 0 {
                        bc.cc.nn.push(node);
                        bc.lock[node as usize] = 2;
                    }
                    bk += 1;
                }
            }
            bid += 1;
            if bc.cc.nn.len() as u64 <= start {
                bc.cc.ng.pop();
                continue;
            }
            let last = bc.cc.ng.last_mut().unwrap();
            *last = (*last << 32) | (bc.cc.nn.len() as u64 - start);
            cocc += 1;
        }
        if cocc > 0 {
            let ng_len = bc.cc.ng.len();
            for k in ng_len - cocc..ng_len {
                let (off, occ) = clus_span(bc, k);
                sc_l.resize(occ, 0.0);
                sc_r.resize(occ, 0.0);
                sc_m.resize(occ, 0.0);
                tmp.resize(occ, 0);
                {
                    let (nn, lock) = (&mut bc.cc.nn, &mut bc.lock);
                    let nodes = &mut nn[off..off + occ];
                    reorder_bub(
                        bc.mg_e,
                        nodes,
                        lock,
                        3,
                        2,
                        4,
                        5,
                        &mut sc_l[..occ],
                        &mut sc_r[..occ],
                        &mut sc_m[..occ],
                        &mut tmp[..occ],
                    );
                }
            }
            for k in iin..bc.cc.nn.len() {
                bc.lock[bc.cc.nn[k] as usize] = 1;
            }
            bc.cc.ng.push(u64::MAX);
            bc.cc.nn.push(u32::MAX);
        }
    }

    if let Some(asn) = &bc.asn {
        let n_asn = asn.len();
        let mut bl = bk;
        let mut cur = bk + 1;
        while cur <= n_asn {
            if cur == n_asn || (asn[cur] >> 32) != (asn[bl] >> 32) {
                let mut cocc = 0usize;
                let iin = bc.cc.nn.len();
                let start = bc.cc.nn.len() as u64;
                bc.cc.ng.push(start);
                for m in bl..cur {
                    let node = asn[m] as u32;
                    if bc.lock[node as usize] == 0 {
                        continue;
                    }
                    bc.cc.nn.push(node);
                    bc.lock[node as usize] = 2;
                }
                if bc.cc.nn.len() as u64 <= start {
                    bc.cc.ng.pop();
                    bl = cur;
                    cur += 1;
                    continue;
                }
                let last = bc.cc.ng.last_mut().unwrap();
                *last = (*last << 32) | (bc.cc.nn.len() as u64 - start);
                cocc += 1;
                if cocc > 0 {
                    let ng_len = bc.cc.ng.len();
                    for k in ng_len - cocc..ng_len {
                        let (off, occ) = clus_span(bc, k);
                        sc_l.resize(occ, 0.0);
                        sc_r.resize(occ, 0.0);
                        sc_m.resize(occ, 0.0);
                        tmp.resize(occ, 0);
                        let (nn, lock) = (&mut bc.cc.nn, &mut bc.lock);
                        let nodes = &mut nn[off..off + occ];
                        reorder_bub(
                            bc.mg_e,
                            nodes,
                            lock,
                            3,
                            2,
                            4,
                            5,
                            &mut sc_l[..occ],
                            &mut sc_r[..occ],
                            &mut sc_m[..occ],
                            &mut tmp[..occ],
                        );
                    }
                    for k in iin..bc.cc.nn.len() {
                        bc.lock[bc.cc.nn[k] as usize] = 1;
                    }
                    bc.cc.ng.push(u64::MAX);
                    bc.cc.nn.push(u32::MAX);
                }
                bl = cur;
            }
            cur += 1;
        }
    }

    let mut tmp2 = vec![0u32; bc.cc.nn.len()];
    let mut idx = vec![0u64; unitig_graph.graph.sequences.len()];
    layout_mc_clus_t(
        bc.mg_e,
        &mut bc.cc.nn,
        &mut bc.sg,
        &mut tmp2,
        &mut idx,
        unitig_graph,
    );

    for &node in &a[..a_n] {
        bc.lock[node as usize] = 0;
    }
}

/// `clus_weight`: the marginal gain of `uid` inside a visited cluster.
fn clus_weight(b_aux: &McSvaux, ma: &McMatch, vis: &BitsP, uid: u32) -> f64 {
    let w = b_aux.s[uid as usize] as f64
        * (b_aux.z[uid as usize][0] - b_aux.z[uid as usize][1])
        * 2.0;
    let mut w_off = 0.0f64;
    let (o, n) = mc_span(ma, uid);
    for j in 0..n {
        let t = ma_y(ma.ma[o + j].x);
        if !is_set_bits_p(vis, t) {
            continue;
        }
        if t == uid {
            continue;
        }
        w_off += b_aux.s[uid as usize] as f64
            * b_aux.s[t as usize] as f64
            * ma.ma[o + j].w;
    }
    w - w_off * 4.0
}

/// `cal_clus_sc0`: scores every prefix of one chain.
#[allow(clippy::too_many_arguments)]
fn cal_clus_sc0(
    ma: &McMatch,
    b_aux: &McSvaux,
    lock: &[u8],
    lock_max: u8,
    a: &[u32],
    a_n: usize,
    id: usize,
    r: &mut ClusFlipAux,
) {
    r.vis.a.fill(0);
    let mut max_occ = u32::MAX;
    let mut w = 0.0f64;
    let mut max_w = -1.0f64;
    let mut i = id;
    while i < a_n && a[i] != u32::MAX {
        let node = a[i];
        if lock[node as usize] >= lock_max {
            i += 1;
            continue;
        }
        if is_set_bits_p(&r.vis, node) {
            i += 1;
            continue;
        }
        w += clus_weight(b_aux, ma, &r.vis, node);
        set_bits_p(&mut r.vis, node);
        if max_w < w {
            max_w = w;
            max_occ = (i + 1 - id) as u32;
        }
        i += 1;
    }
    let mut i = id;
    while i < a_n && a[i] != u32::MAX {
        r.vis.a[(a[i] >> 3) as usize] = 0;
        i += 1;
    }
    if max_w <= 0.000001 || max_occ == u32::MAX {
        return;
    }
    if max_w > r.w || (max_w == r.w && (id as u32) < r.off) {
        r.w = max_w;
        r.off = id as u32;
        r.occ = max_occ;
    }
}

/// `gen_best_clus`: the best chain prefix among all start positions.
fn gen_best_clus(bc: &mut McClus, b_aux: &McSvaux) -> Option<(u32, u32, f64)> {
    for aux in bc.aux.iter_mut() {
        clean_clus_flip_aux(aux);
    }
    let a_n = bc.cc.nn.len();
    for i in 0..a_n {
        if bc.cc.nn[i] == u32::MAX {
            continue;
        }
        let mut aux = std::mem::take(&mut bc.aux[0]);
        cal_clus_sc0(
            bc.mg_e,
            b_aux,
            &bc.lock,
            bc.lock_max,
            &bc.cc.nn,
            a_n,
            i,
            &mut aux,
        );
        bc.aux[0] = aux;
    }
    let mut off = u32::MAX;
    let mut occ = u32::MAX;
    let mut rw = -1.0f64;
    for aux in &bc.aux {
        if aux.off == u32::MAX {
            continue;
        }
        if aux.w < 0.0 {
            continue;
        }
        if aux.w > rw || (aux.w == rw && aux.off < off) {
            off = aux.off;
            occ = aux.occ;
            rw = aux.w;
        }
    }
    if off != u32::MAX {
        Some((off, occ, rw))
    } else {
        None
    }
}

/// `flip_chain`: flips the chosen prefix and returns the new score.
fn flip_chain(
    ma: &McMatch,
    b: &mut McSvaux,
    a: &[u32],
    off: u32,
    occ: u32,
    vis: &mut BitsP,
    lock: &mut [u8],
    lock_max: u8,
) -> f64 {
    let kn = off + occ;
    for k in off..kn {
        vis.a[(a[k as usize] >> 3) as usize] = 0;
    }
    for k in off..kn {
        let node = a[k as usize];
        if lock[node as usize] >= lock_max {
            continue;
        }
        if is_set_bits_p(vis, node) {
            continue;
        }
        set_bits_p(vis, node);
        lock[node as usize] += 1;
        let spin = -b.s[node as usize];
        mc_set_spin(ma, b, node, spin);
    }
    mc_score(ma, b)
}

/// `mc_solve_clus`.
fn mc_solve_clus(bc: &mut McClus, b: &mut McSvaux) -> f64 {
    let mut sc_opt = mc_score(bc.mg_e, b);
    bc.lock.fill(0);
    loop {
        let Some((off, occ, _)) = gen_best_clus(bc, b) else {
            break;
        };
        // The reference reuses `aux[0].vis` for both the chain scoring and the
        // flip, so it must be live while `gen_best_clus` runs; take it only
        // around the flip.
        let mut vis = std::mem::take(&mut bc.aux[0].vis);
        let sc = flip_chain(
            bc.mg_e,
            b,
            &bc.cc.nn,
            off,
            occ,
            &mut vis,
            &mut bc.lock,
            bc.lock_max,
        );
        bc.aux[0].vis = vis;
        if sc < sc_opt {
            eprintln!(
                "\nwrong::[M::mc_solve_clus::] sc_opt::{} sc::{} off::{} occ::{}",
                sc_opt, sc, off, occ
            );
        }
        sc_opt = sc;
    }
    sc_opt
}

/// `mc_clus_cc`.
fn mc_clus_cc(bc: &mut McClus, b: &mut McSvaux) -> f64 {
    mc_reset_z(bc.mg_e, b);
    let mut sc_opt = mc_score(bc.mg_e, b);
    let mut sc;
    loop {
        sc = mc_solve_clus(bc, b);
        if sc <= sc_opt + 0.0000001 {
            break;
        }
        sc_opt = sc;
    }
    sc
}

/// `mc_solve_cc_adv`: component optimisation with occasional cluster
/// refinement.
fn mc_solve_cc_adv(
    opt: &McOpt,
    mg: &McG,
    b: &mut McSvaux,
    cc_off: u32,
    cc_size: u32,
    mut bc: Option<&mut McClus>,
    unitig_graph: Option<&UnitigGraph>,
) -> u32 {
    b.cc_off = cc_off;
    b.cc_size = cc_size;
    if cc_size < 2 {
        return 0;
    }
    let mut sc_opt = mc_init_spin(&mg.e, b);
    if cc_size == 2 {
        return 0;
    }
    for j in 0..cc_size {
        let k = b.cc_node[j as usize] as usize;
        b.s_opt[k] = b.s[k];
        b.z_opt[k] = b.z[k];
    }

    if let Some(cluster) = bc.as_deref_mut() {
        let unitig_graph =
            unitig_graph.expect("the bubble-chain layout requires the unitig graph");
        renew_mc_clus_t_adv(cluster, &b.cc_node[..cc_size as usize], cc_size as usize, unitig_graph);
    }

    let mut n_iter = 0u32;
    let mut sc = mc_optimize_local(opt, &mg.e, b, &mut n_iter);
    if sc > sc_opt {
        for j in 0..cc_size {
            let k = b.cc_node[j as usize] as usize;
            b.s_opt[k] = b.s[k];
            b.z_opt[k] = b.z[k];
        }
        sc_opt = sc;
    } else {
        for j in 0..cc_size {
            let k = b.cc_node[j as usize] as usize;
            b.s[k] = b.s_opt[k];
            b.z[k] = b.z_opt[k];
        }
    }

    if let Some(cluster) = bc.as_deref_mut() {
        let sc2 = mc_clus_cc(cluster, b);
        if sc2 > sc_opt {
            for j in 0..cc_size {
                let k = b.cc_node[j as usize] as usize;
                b.s_opt[k] = b.s[k];
                b.z_opt[k] = b.z[k];
            }
            sc_opt = sc2;
        } else {
            for j in 0..cc_size {
                let k = b.cc_node[j as usize] as usize;
                b.s[k] = b.s_opt[k];
                b.z[k] = b.z_opt[k];
            }
        }
    }

    let flush = (opt.max_iter * 50) as u32;
    let n_skip_flush = (opt.n_perturb / 16) as u32;
    let mut n_skip = 0u32;
    for k in 0..opt.n_perturb as u32 {
        if k & 1 != 0 {
            mc_perturb(opt, &mg.e, b);
        } else {
            mc_perturb_node(opt, &mg.e, b, 3);
        }
        sc = mc_optimize_local(opt, &mg.e, b, &mut n_iter);
        if sc > sc_opt {
            for j in 0..cc_size {
                let node = b.cc_node[j as usize] as usize;
                b.s_opt[node] = b.s[node];
                b.z_opt[node] = b.z[node];
            }
            sc_opt = sc;
            n_skip = 0;
        } else {
            for j in 0..cc_size {
                let node = b.cc_node[j as usize] as usize;
                b.s[node] = b.s_opt[node];
                b.z[node] = b.z_opt[node];
            }
            n_skip += 1;
        }
        if n_skip >= n_skip_flush && bc.is_some() {
            if let Some(cluster) = bc.as_deref_mut() {
                let sc2 = mc_clus_cc(cluster, b);
                if sc2 > sc_opt {
                    for j in 0..cc_size {
                        let node = b.cc_node[j as usize] as usize;
                        b.s_opt[node] = b.s[node];
                        b.z_opt[node] = b.z[node];
                    }
                    sc_opt = sc2;
                } else {
                    for j in 0..cc_size {
                        let node = b.cc_node[j as usize] as usize;
                        b.s[node] = b.s_opt[node];
                        b.z[node] = b.z_opt[node];
                    }
                }
            }
            n_skip = 0;
        }
        if n_iter % flush == 0 {
            mc_reset_z(&mg.e, b);
            sc = mc_score(&mg.e, b);
            for j in 0..cc_size {
                let node = b.cc_node[j as usize] as usize;
                b.z_opt[node] = b.z[node];
            }
            sc_opt = sc;
        }
    }

    if let Some(cluster) = bc.as_deref_mut() {
        let sc2 = mc_clus_cc(cluster, b);
        if sc2 > sc_opt {
            for j in 0..cc_size {
                let node = b.cc_node[j as usize] as usize;
                b.s_opt[node] = b.s[node];
                b.z_opt[node] = b.z[node];
            }
        } else {
            for j in 0..cc_size {
                let node = b.cc_node[j as usize] as usize;
                b.s[node] = b.s_opt[node];
                b.z[node] = b.z_opt[node];
            }
        }
    }

    for j in 0..cc_size {
        let node = b.cc_node[j as usize] as usize;
        b.s[node] = b.s_opt[node];
        b.z[node] = b.z_opt[node];
    }
    n_iter
}

/// `gen_ref_bub`: the cluster id of every non-bubble node from the reference
/// transitions.
fn gen_ref_bub(bc: &McClus, refs: &[UTrans], ref_idx: &[u64]) -> Vec<u64> {
    let n_seq = bc.n as usize;
    let mut aux = vec![u64::MAX; bc.n as usize];
    let empty: &[MaUtg] = &[];
    let unitigs = bc.bub.bubble_unitigs.as_deref().unwrap_or(empty);
    let mut bid = 0usize;
    for unitig in unitigs {
        for &packed in &unitig.reads {
            let (_, _, body, _) = get_bubbles(bc.bub, packed >> 33);
            for &node in body {
                aux[(node >> 1) as usize] = bid as u64;
            }
            bid += 1;
        }
    }

    for i in 0..n_seq {
        if aux[i] != u64::MAX {
            continue;
        }
        let entries = transitions_of(refs, ref_idx, i as u32);
        if entries.is_empty() {
            continue;
        }
        let mut buf: Vec<u64> = Vec::with_capacity(entries.len());
        for (k, entry) in entries.iter().enumerate() {
            let p = aux[entry.target_identifier as usize] << 32 | k as u64;
            buf.push(p);
        }
        radix_sort_by_key(&mut buf, &|value: &u64| *value);
        let mut mm = u32::MAX;
        let mut mmw = -1.0f64;
        let mut l = 0usize;
        for k in 1..=buf.len() {
            if k == buf.len()
                || (buf[l] >> 32) == u32::MAX as u64
                || (buf[l] >> 32) != (buf[k] >> 32)
            {
                let mut w = 0.0f64;
                for &entry in &buf[l..k] {
                    w += entries[entry as u32 as usize].network_weight.abs();
                }
                if mmw < w {
                    mmw = w;
                    mm = (buf[l] >> 32) as u32;
                }
                l = k;
            }
        }
        if mm != u32::MAX {
            aux[i] = mm as u64;
        }
    }

    let mut stack: Vec<u64> = Vec::new();
    for i in 0..n_seq {
        if aux[i] != u64::MAX {
            continue;
        }
        stack.clear();
        stack.push(i as u64);
        while let Some(k) = stack.pop() {
            if aux[k as usize] != u64::MAX {
                continue;
            }
            aux[k as usize] = bid as u64;
            let entries = transitions_of(refs, ref_idx, k as u32);
            let mut st = 0usize;
            for j in 1..=entries.len() {
                if j == entries.len()
                    || entries[j].target_identifier != entries[st].target_identifier
                {
                    if aux[entries[st].target_identifier as usize] == u64::MAX {
                        stack.push(entries[st].target_identifier as u64);
                    }
                    st = j;
                }
            }
        }
        bid += 1;
    }

    let empty: &[MaUtg] = &[];
    let unitigs = bc.bub.bubble_unitigs.as_deref().unwrap_or(empty);
    for unitig in unitigs {
        for &packed in &unitig.reads {
            let (_, _, body, _) = get_bubbles(bc.bub, packed >> 33);
            for &node in body {
                aux[(node >> 1) as usize] = u64::MAX;
            }
        }
    }
    let mut out: Vec<u64> = Vec::new();
    for (i, &cluster) in aux.iter().enumerate() {
        if cluster == u64::MAX {
            continue;
        }
        out.push((cluster << 32) | i as u64);
    }
    radix_sort_by_key(&mut out, &|value: &u64| *value);
    out
}

/// `mc_solve_core_adv`: the unitig-level pass with the bubble-chain layout.
fn mc_solve_core_adv(
    opt: &McOpt,
    mg: &mut McG,
    unitig_graph: Option<&UnitigGraph>,
    bub: Option<&BubbleType>,
    refs: &[UTrans],
    ref_idx: &[u64],
    n_thread: u64,
) {
    mg.e.cc = mc_g_cc_core(&mg.e);
    let mut b = mc_svaux_init(mg, opt.seed);
    let mut bc = match (bub, unitig_graph) {
        (Some(bub), Some(unitig_graph)) => Some(init_mc_clus_t(
            &mg.e,
            bub,
            n_thread,
            unitig_graph.graph.sequences.len() as u32,
            16,
        )),
        _ => None,
    };
    if let Some(cluster) = bc.as_mut() {
        if !refs.is_empty() && !ref_idx.is_empty() {
            cluster.asn = Some(gen_ref_bub(cluster, refs, ref_idx));
        }
    }

    let n = mg.e.n_seq as usize;
    let mut st = 0usize;
    for i in 1..=n {
        if i == n || mg.e.cc[st] >> 32 != mg.e.cc[i] >> 32 {
            mc_solve_cc_adv(
                opt,
                mg,
                &mut b,
                st as u32,
                (i - st) as u32,
                bc.as_mut(),
                unitig_graph,
            );
            st = i;
        }
    }
    mg.s = b.s.clone();
}

/// The `mc_solve` entry used by `trio_phasing_refine`.
///
/// `ta` is the transition table built by the trio stage, `units`/`units_idx`
/// are the reference transitions (`ref`), and `s` is `ps_t::s`, which the
/// solver updates in place.
#[allow(clippy::too_many_arguments)]
pub fn mc_solve(
    unitig_graph: &UnitigGraph,
    ta: &[UTrans],
    refs: &[UTrans],
    ref_idx: &[u64],
    s: &mut [McNode],
    is_sys: bool,
    bub: &BubbleType,
    n_thread: u64,
    n_perturb: i32,
    f_perturb: f64,
    seed: u64,
) {
    let mut opt = mc_opt_init(n_perturb, f_perturb, seed);
    let n_seq = unitig_graph.graph.sequences.len() as u32;
    let mut mg = init_mc_g_t(n_seq, s);
    mg.e.n_seq = n_seq;
    update_mc_edges(&mut mg, ta, is_sys);
    eprintln!("[M::mc_solve:: # edges: {}]", mg.e.ma.len());

    if !refs.is_empty() && !ref_idx.is_empty() {
        mb_solve_core(&mut opt, &mut mg, refs, ref_idx, is_sys);
    }
    mc_solve_core_adv(
        &opt,
        &mut mg,
        Some(unitig_graph),
        Some(bub),
        refs,
        ref_idx,
        n_thread,
    );
    let n = s.len().min(mg.s.len());
    s[..n].copy_from_slice(&mg.s[..n]);
}

/// The `mc_solve` entry used by `purge_dups`:
/// `mc_solve(ovlp, t_ch, NULL, ug, read_g, 0.8, trio_flag, 1, NULL, 1, NULL,
/// NULL, 1, 0)`.
///
/// The purge call has no reference transitions (so the block pass returns
/// immediately) and no bubble chain (so the component pass runs without the
/// cluster refinement). The fresh spin array starts at zero (`renew_s = 1`
/// with `s = NULL`). On return the conflicting hap overlaps are compacted in
/// place by [`clean_ovlp_by_mc`], as `clean_ov = 1` requests.
#[allow(clippy::too_many_arguments)]
pub fn mc_solve_purge(
    unitig_graph: &UnitigGraph,
    ha: &mut [Vec<HapOverlap>],
    het: &[u8],
    read_lengths: &[u32],
    f_rate: f64,
    polyploidy: i32,
    n_perturb: i32,
    f_perturb: f64,
    seed: u64,
    mut trio_flags: Option<&mut [u8]>,
) {
    let opt = mc_opt_init(n_perturb, f_perturb, seed);
    let n_seq = unitig_graph.graph.sequences.len() as u32;
    let mut mg = init_mc_g_t(n_seq, &[]);
    mg.e.n_seq = n_seq;
    update_mc_edges_ovlp(
        &mut mg,
        ha,
        Some(het),
        read_lengths,
        &unitig_graph.unitigs,
        f_rate,
        polyploidy,
        true,
    );
    eprintln!("[M::mc_solve:: # edges: {}]", mg.e.ma.len());

    mc_solve_core_adv(&opt, &mut mg, None, None, &[], &[], 1);

    if let Some(trio_flags) = trio_flags.as_deref_mut() {
        p_nodes(&mg, &unitig_graph.unitigs, het, trio_flags);
    }
    clean_ovlp_by_mc(&mg, ha);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splitmix_matches_the_reference_output() {
        // Standard SplitMix64: the first output of seed 0.
        assert_eq!(kr_splitmix64(0), 0xE220_A839_7B1D_CDAF);
        let mut x = 11u64;
        let first = kr_drand_r(&mut x);
        assert!((0.0..1.0).contains(&first));
        assert_eq!(x, kr_splitmix64(11));
    }

    #[test]
    fn connected_components_follow_vertex_order() {
        // 0 - 1 - 2 and an isolated 3.
        let mut ma = McMatch {
            n_seq: 4,
            ..Default::default()
        };
        ma.ma = vec![
            McEdge { x: (0u64 << 32) | 1, w: 1.0 },
            McEdge { x: (1u64 << 32) | 0, w: 1.0 },
            McEdge { x: (1u64 << 32) | 2, w: 1.0 },
            McEdge { x: (2u64 << 32) | 1, w: 1.0 },
        ];
        mc_edges_idx(&mut ma.idx, &ma.ma, 4);
        let cc = mc_g_cc_core(&ma);
        // Group ids are assigned in the order of the first vertex of the
        // component; vertices are sorted inside each group.
        assert_eq!(cc, vec![(0u64 << 32) | 0, (0 << 32) | 1, (0 << 32) | 2, (1 << 32) | 3]);
    }
}

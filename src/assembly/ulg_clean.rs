//! The `usg` graph-cleaning cluster of `u2g_hybrid_clean` (`gfa_ut.cpp`
//! 10519-11150 + 15332 + 15616-16050): tip cuts, bubble identification and
//! the bridge marking over the hybrid UL/unitig graph (`usg_t`).

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::bubble::BubbleBuffer;
use crate::assembly::bubble_chain::BubbleType;
use crate::assembly::ul2ul::{Ul2ulIdx, UlBg, UlCov, Usg, UsgArc};
use crate::container::ksort::radix_sort_u64;

/// `bg_*` flags (`gfa_ut.cpp` 8617-8620).
pub const BG_UNAVAILABLE: u32 = u32::MAX;
pub const BG_CORRECT: u32 = 1;
pub const BG_WRONG: u32 = 2;
pub const BG_AMBIGUOUS: u32 = 3;

/// `get_arcs` (`gfa_ut.cpp` 463): the live out-degree, optionally recording
/// the arc indices.
pub fn get_arcs(g: &Asg, v: u32, idx: Option<&mut Vec<usize>>, idx_n: usize) -> usize {
    let arcs = g.arcs_of(v);
    let mut kv = 0usize;
    if let Some(list) = idx {
        for (i, arc) in arcs.iter().enumerate() {
            if arc.deleted {
                continue;
            }
            if kv < idx_n {
                list.push(i);
            }
            kv += 1;
        }
    } else {
        kv = arcs.iter().filter(|a| !a.deleted).count();
    }
    kv
}

/// `usg_end` (`gfa_ut.cpp` 10698).
fn usg_end(g: &Usg, v: u32, lw: Option<&mut u64>) -> i32 {
    let arcs = g.arcs_of(v ^ 1);
    let mut i0 = -1i32;
    let mut nv = 0usize;
    for (i, arc) in arcs.iter().enumerate() {
        if !arc.deleted {
            i0 = i as i32;
            nv += 1;
        }
    }
    if nv == 0 {
        return 1; // ASG_ET_TIP
    }
    if nv > 1 {
        return 3; // ASG_ET_MULTI_OUT
    }
    if let Some(out) = lw {
        *out = ((u64::from(v ^ 1)) << 32) | u64::from(arcs[i0 as usize].v);
    }
    let w = arcs[i0 as usize].v ^ 1;
    let nw = g.arcs_of(w).iter().filter(|a| !a.deleted).count();
    if nw != 1 {
        return 4; // ASG_ET_MULTI_NEI
    }
    0 // ASG_ET_MERGEABLE
}

/// `usg_arc_occ` (`gfa_ut.h`): the live out-degree of `v`.
fn usg_arc_occ(g: &Usg, v: u32, lw: Option<&mut u64>) -> usize {
    let arcs = g.arcs_of(v);
    let mut count = 0usize;
    let mut chosen = 0usize;
    for (i, arc) in arcs.iter().enumerate() {
        if arc.deleted {
            continue;
        }
        count += 1;
        chosen = i;
    }
    if count == 1 {
        if let Some(out) = lw {
            *out = ((u64::from(v)) << 32) | chosen as u64;
        }
    }
    count
}

/// `get_usg_arc` (`gfa_ut.cpp` 10551).
pub fn get_usg_arc<'a>(g: &'a Usg, v: u32, w: u32) -> Option<&'a UsgArc> {
    g.arcs_of(v).iter().find(|arc| arc.v == w)
}

/// `usg_seq_del` (`gfa_ut.cpp` 10569): deletes one `usg` vertex with both
/// arc sides; returns the vertex's occupancy for the caller's counters.
pub fn usg_seq_del(g: &mut Usg, s: u32) -> u32 {
    g.a[s as usize].deleted = true;
    for side in 0..2usize {
        for arc in g.a[s as usize].arc[side].iter_mut() {
            arc.deleted = true;
        }
    }
    g.a[s as usize].occ
}

/// `usg_tip_del` (`gfa_ut.cpp` 15782): the walk-and-delete of one tip.
pub fn usg_tip_del(g: &mut Usg, v: u32, z_a: &mut Vec<u32>, z_b: &mut Vec<u32>) -> i32 {
    let a_n = z_a.len();
    let b_n = z_b.len();
    let mut n_ext = 0i32;
    let kv = g
        .arcs_of(v ^ 1)
        .iter()
        .filter(|a| !a.deleted && !g.a[(a.v >> 1) as usize].deleted)
        .count();
    if kv > 1 {
        z_a.truncate(a_n);
        z_b.truncate(b_n);
        return 0;
    }
    n_ext += g.a[v as usize].occ as i32;
    g.a[v as usize].deleted = true;
    z_b.push(v >> 1);
    for arc in g.arcs_of(v) {
        if !arc.deleted && !g.a[(arc.v >> 1) as usize].deleted {
            z_a.push(arc.v);
        }
    }
    while z_a.len() > a_n {
        let v = z_a.pop().unwrap();
        if g.a[(v >> 1) as usize].deleted {
            continue;
        }
        let kv = g
            .arcs_of(v ^ 1)
            .iter()
            .filter(|a| !a.deleted && !g.a[(a.v >> 1) as usize].deleted)
            .count();
        if kv > 0 {
            continue;
        }
        n_ext += g.a[(v >> 1) as usize].occ as i32;
        g.a[(v >> 1) as usize].deleted = true;
        z_b.push(v >> 1);
        for arc in g.arcs_of(v) {
            if !arc.deleted && !g.a[(arc.v >> 1) as usize].deleted {
                z_a.push(arc.v);
            }
        }
    }
    z_b.truncate(b_n);
    n_ext
}

/// `usg_tip_detect` (`gfa_ut.cpp` 10871 region): the read-only variant of
/// [`usg_tip_del`]; `f` marks the visited vertices.
pub fn usg_tip_detect(
    g: &Usg,
    v0: u32,
    max_ext: i32,
    f: &mut [u8],
    z_a: &mut Vec<u32>,
    z_b: &mut Vec<u32>,
) -> i32 {
    let a_n = z_a.len();
    let b_n = z_b.len();
    let mut v = v0;
    let mut n_ext = 0i32;
    let kv = g
        .arcs_of(v ^ 1)
        .iter()
        .filter(|a| !a.deleted && f[(a.v >> 1) as usize] == 0)
        .count();
    if kv > 1 {
        z_a.truncate(a_n);
        z_b.truncate(b_n);
        return 0;
    }
    n_ext += g.a[(v >> 1) as usize].occ as i32;
    f[(v >> 1) as usize] = 1;
    z_b.push(v >> 1);
    for arc in g.arcs_of(v) {
        if !arc.deleted && f[(arc.v >> 1) as usize] == 0 {
            z_a.push(arc.v);
        }
    }
    while z_a.len() > a_n && n_ext < max_ext {
        let v = z_a.pop().unwrap();
        if f[(v >> 1) as usize] != 0 {
            continue;
        }
        let kv = g
            .arcs_of(v ^ 1)
            .iter()
            .filter(|a| !a.deleted && f[(a.v >> 1) as usize] == 0)
            .count();
        if kv > 0 {
            continue;
        }
        n_ext += g.a[(v >> 1) as usize].occ as i32;
        f[(v >> 1) as usize] = 1;
        z_b.push(v >> 1);
        for arc in g.arcs_of(v) {
            if !arc.deleted && f[(arc.v >> 1) as usize] == 0 {
                z_a.push(arc.v);
            }
        }
    }
    for i in b_n..z_b.len() {
        f[z_b[i] as usize] = 0;
    }
    z_b.truncate(b_n);
    z_a.truncate(a_n);
    n_ext
}

/// `usg_naive_topocut_aux` (`gfa_ut.cpp` 10871).
pub fn usg_naive_topocut_aux(
    g: &Usg,
    v: u32,
    max_ext: i32,
    f: &mut [u8],
    b0: &mut Vec<u32>,
    b1: &mut Vec<u32>,
) -> i32 {
    let mut v = v;
    let v0 = v;
    let mut w;
    let mut n_ext = 0i32;
    let mut tip = false;
    while n_ext < max_ext {
        tip = false;
        let kv = g
            .arcs_of(v ^ 1)
            .iter()
            .filter(|a| !a.deleted)
            .count();
        if kv != 1 {
            break;
        }
        n_ext += g.a[(v >> 1) as usize].occ as i32;
        let av: Vec<u32> = g
            .arcs_of(v)
            .iter()
            .filter(|a| !a.deleted)
            .map(|a| a.v)
            .collect();
        let kv = av.len();
        if kv != 1 {
            if kv > 1 {
                tip = true;
            }
            break;
        }
        w = av[0];
        v = w;
    }
    if n_ext < max_ext && tip {
        n_ext = usg_tip_detect(g, v0, max_ext, f, b0, b1);
    }
    n_ext
}

/// `usg_arc_cut_tips` (`gfa_ut.cpp` 10717).
pub fn usg_arc_cut_tips(g: &mut Usg, max_ext: u32, ignore_ul: bool, b: &mut Vec<u64>) -> u32 {
    let n_vtx = (g.a.len() << 1) as u32;
    b.clear();
    for v in 0..n_vtx {
        if g.a[(v >> 1) as usize].deleted {
            continue;
        }
        let kv = g.arcs_of(v ^ 1).iter().filter(|a| !a.deleted).count();
        if kv > 0 {
            continue;
        }
        let mut is_telo = g.a[(v >> 1) as usize].telo;
        let mut w = v;
        let mut kv = g.a[(v >> 1) as usize].occ;
        for _ in 0..max_ext {
            let mut lw = 0u64;
            if usg_end(g, w ^ 1, Some(&mut lw)) != 0 {
                break;
            }
            w = (lw & 0xffff_ffff) as u32;
            kv += g.a[(w >> 1) as usize].occ;
            if g.a[(w >> 1) as usize].telo {
                is_telo = true;
            }
        }
        if kv <= max_ext && !is_telo {
            b.push(((kv as u64) << 32) | u64::from(v));
        }
    }
    crate::container::ksort::radix_sort_u64(b);
    let mut cnt = 0u32;
    for k in 0..b.len() {
        let v = (b[k] & 0xffff_ffff) as u32;
        if g.a[(v >> 1) as usize].deleted {
            continue;
        }
        let kv = g.arcs_of(v ^ 1).iter().filter(|a| !a.deleted).count();
        if kv > 0 {
            continue;
        }
        let pb = b.len();
        b.push(u64::from(v));
        let mut is_telo = g.a[(v >> 1) as usize].telo;
        let mut w = v;
        let mut kv = g.a[(v >> 1) as usize].occ;
        for _ in 0..max_ext {
            let mut lw = 0u64;
            if usg_end(g, w ^ 1, Some(&mut lw)) != 0 {
                break;
            }
            w = (lw & 0xffff_ffff) as u32;
            kv += g.a[(w >> 1) as usize].occ;
            b.push(lw);
            if g.a[(w >> 1) as usize].telo {
                is_telo = true;
            }
        }
        if kv <= max_ext && !is_telo {
            let mut ff = false;
            if !ignore_ul {
                let mut i = pb;
                while i + 1 < b.len() {
                    if let Some(p) = get_usg_arc(g, (b[i] & 0xffff_ffff) as u32, (b[i + 1] & 0xffff_ffff) as u32)
                    {
                        if p.ou > 1 {
                            ff = true;
                            break;
                        }
                    }
                    i += 1;
                }
                if !ff && i < b.len() {
                    let last = (b[i] & 0xffff_ffff) as u32;
                    for arc in g.arcs_of(last) {
                        if !arc.deleted && arc.ou > 1 {
                            ff = true;
                            break;
                        }
                    }
                }
            }
            if !ff
                && pb < b.len()
                && !usg_real_tip(g, (b[b.len() - 1] & 0xffff_ffff) as u32, 0.75)
            {
                ff = true;
            }
            if !ff {
                for i in pb..b.len() {
                    usg_seq_del(g, ((b[i] & 0xffff_ffff) as u32) >> 1);
                }
                cnt += 1;
            }
        }
        b.truncate(pb);
    }
    if cnt > 0 {
        usg_cleanup(g);
    }
    cnt
}

/// `usg_real_tip` (`gfa_ut.cpp` 10702): whether the tip is a real (read-end)
/// tip by the occupancy ratio.
pub fn usg_real_tip(g: &Usg, v0: u32, rate: f64) -> bool {
    let arcs_v = g.arcs_of(v0);
    let mut ov = 0u64;
    for arc in arcs_v {
        if !arc.deleted {
            ov += u64::from(arc.ol);
        }
    }
    let arcs_w = g.arcs_of(v0 ^ 1);
    let mut ow = 0u64;
    for arc in arcs_w {
        if !arc.deleted {
            ow += u64::from(arc.ol);
        }
    }
    let total = ov + ow;
    if total == 0 {
        return false;
    }
    let small = ov.min(ow) as f64;
    small / total as f64 <= rate
}

/// `usg_cleanup` (`gfa_ut.cpp` 10663): compacts every arc list, sorts it and
/// merges the `arc_mm` entries into the `idx` bookkeeping.
pub fn usg_cleanup(g: &mut Usg) {
    for seq in g.a.iter_mut() {
        if seq.deleted {
            seq.arc[0].clear();
            seq.arc[1].clear();
        }
        for side in 0..2usize {
            let kept: Vec<UsgArc> = seq.arc[side]
                .iter()
                .copied()
                .filter(|a| !a.deleted)
                .map(|mut a| {
                    a.idx = 0;
                    a
                })
                .collect();
            seq.arc[side] = kept;
        }
        for side in 0..2usize {
            crate::container::ksort::radix_sort_by_key_32(
                &mut seq.arc[side],
                &|arc: &UsgArc| arc.v,
            );
        }
    }
}

/// `extracr_clean_arc` (`gfa_ut.cpp` 15757): collects the zero-occupancy
/// arcs of `v` whose sibling arcs carry bridged occupancy.
pub fn extracr_clean_arc(g: &Usg, v: u32, b: &mut Vec<u64>, ub: &mut Vec<u64>) {
    let arcs = g.arcs_of(v);
    let nv = arcs.len();
    if nv < 2 {
        return;
    }
    let mut kv = 0usize;
    let mut kocc = [0usize; 2];
    for arc in arcs {
        if arc.deleted {
            continue;
        }
        kv += 1;
        if (arc.ou >> 1) > 0 {
            kocc[1] = kocc[1].max((arc.ou >> 1) as usize);
        } else if arc.ou == 0 {
            kocc[0] += 1;
        }
    }
    if kv < 2 || kocc[0] == 0 || kocc[1] == 0 {
        return;
    }
    for (i, arc) in arcs.iter().enumerate() {
        if arc.deleted || arc.ou != 0 {
            continue;
        }
        let x = (u64::from(arc.ol)) << 32;
        b.push(x | ub.len() as u64);
        ub.push(((u64::from(v)) << 32) | i as u64);
    }
}

/// `usg_bub_clean0` (`gfa_ut.cpp` 15824): the arc-level cleanup inside one
/// bubble.
#[allow(clippy::too_many_arguments)]
pub fn usg_bub_clean0(
    g: &mut Usg,
    a: &[u32],
    sid: u32,
    eid: u32,
    is_topo: bool,
    max_ext: i32,
    len_rat: f32,
    f: &mut [u8],
    b: &mut Vec<u64>,
    ub: &mut Vec<u64>,
) -> u64 {
    let mut z_a: Vec<u32> = Vec::new();
    let mut z_b: Vec<u32> = Vec::new();
    b.clear();
    ub.clear();
    for &vertex in a {
        extracr_clean_arc(g, vertex, b, ub);
        extracr_clean_arc(g, vertex ^ 1, b, ub);
    }
    extracr_clean_arc(g, sid, b, ub);
    extracr_clean_arc(g, eid, b, ub);
    crate::container::ksort::radix_sort_u64(b);
    let mut cnt = 0u64;
    for k in 0..b.len() {
        let slot = (b[k] & 0xffff_ffff) as usize;
        let v = (ub[slot] >> 32) as u32;
        let arc_index = (ub[slot] & 0xffff_ffff) as usize;
        let ve = g.a[(v >> 1) as usize].arc[(v & 1) as usize][arc_index];
        let w = ve.v ^ 1;
        if ve.deleted
            || g.a[(v >> 1) as usize].deleted
            || g.a[(w >> 1) as usize].deleted
            || ve.ou != 0
        {
            continue;
        }
        let nv = g.arcs_of(v).len();
        let nw = g.arcs_of(w).len();
        if nv <= 1 && nw <= 1 {
            continue;
        }
        let we_opt = g.arcs_of(w).iter().position(|a| a.v == v ^ 1);
        let Some(we_index) = we_opt else {
            continue;
        };
        let we = g.a[(w >> 1) as usize].arc[(w & 1) as usize][we_index];
        let mm_ol = ve.ol.min(we.ol);
        let mut kocc = [0usize; 2];
        let mut kv = 0usize;
        let mut ol_max = 0u32;
        for arc in g.arcs_of(v) {
            if arc.deleted {
                continue;
            }
            kv += 1;
            if arc.ou != 1 {
                kocc[usize::from(arc.ou != 0)] += 1;
            }
            if ol_max < arc.ol {
                ol_max = arc.ol;
            }
        }
        if kv < 1 || kocc[0] < 1 || kocc[1] < 1 {
            continue;
        }
        if kv >= 2 && (mm_ol as f32) > ol_max as f32 * len_rat {
            continue;
        }
        let mut kw = 0usize;
        let mut ol_max = 0u32;
        for arc in g.arcs_of(w) {
            if arc.deleted {
                continue;
            }
            kw += 1;
            if ol_max < arc.ol {
                ol_max = arc.ol;
            }
        }
        if kw < 1 {
            continue;
        }
        if kw >= 2 && (mm_ol as f32) > ol_max as f32 * len_rat {
            continue;
        }
        if kv <= 1 && kw <= 1 {
            continue;
        }
        let mut to_del = true;
        let mut n_tip = 0i32;
        let mut tip_v = u32::MAX;
        if is_topo {
            to_del = false;
            if kv > 1 && kw > 1 {
                to_del = true;
            } else if kw == 1 {
                n_tip = usg_naive_topocut_aux(g, w ^ 1, max_ext, f, &mut z_a, &mut z_b);
                if n_tip < max_ext {
                    to_del = true;
                    tip_v = w ^ 1;
                }
            } else if kv == 1 {
                n_tip = usg_naive_topocut_aux(g, v ^ 1, max_ext, f, &mut z_a, &mut z_b);
                if n_tip < max_ext {
                    to_del = true;
                    tip_v = v ^ 1;
                }
            }
        }
        if to_del {
            if n_tip > 0 && tip_v != u32::MAX {
                let r_tip = usg_tip_del(g, tip_v, &mut z_a, &mut z_b);
                debug_assert_eq!(n_tip, r_tip);
            }
            cnt += 1;
        }
    }
    cnt
}

fn b0buf(b: &mut Vec<u64>) -> &mut Vec<u64> {
    b
}

fn ubbuf(ub: &mut Vec<u64>) -> &mut Vec<u64> {
    ub
}

/// `usg_bub_clean` (`gfa_ut.cpp` 15917): identifies bubbles and cleans the
/// interior of each.
pub fn usg_bub_clean(
    g: &mut Usg,
    max_ext: i32,
    len_rat: f32,
    is_topo: bool,
    bs: &mut [u8],
    f: &mut [u8],
    b: &mut Vec<u64>,
    ub: &mut Vec<u64>,
) -> u64 {
    let mut n_cut = 0u64;
    let max_dist = usg_bub_dis(g);
    for slot in bs.iter_mut() {
        *slot = 0;
    }
    if max_dist > 0 {
        let n_vtx = (g.a.len() << 1) as u32;
        for v in 0..n_vtx {
            if bs[v as usize] != 0 {
                continue;
            }
            if g.a[(v >> 1) as usize].deleted {
                continue;
            }
            let nv = g.arcs_of(v).len();
            if nv < 2 {
                continue;
            }
            let n_arc = g.arcs_of(v).iter().filter(|a| !a.deleted).count();
            if n_arc < 2 {
                continue;
            }
            let Some((sink, body)) = usg_bub_identify(g, v, max_dist, u64::MAX) else {
                continue;
            };
            let mut compacted: Vec<u32> = Vec::new();
            for &node in &body {
                if node == v || node == sink {
                    continue;
                }
                bs[node as usize] = 1;
                bs[(node ^ 1) as usize] = 1;
                compacted.push(node);
            }
            bs[v as usize] = 2;
            bs[(sink ^ 1) as usize] = 3;
            n_cut += usg_bub_clean0(
                g, &compacted, v, sink ^ 1, is_topo, max_ext, len_rat, f, b, ub,
            );
        }
    }
    if n_cut > 0 {
        usg_cleanup(g);
    }
    n_cut
}

/// `usg_bub_dis` (`gfa_ut.cpp` 15616).
pub fn usg_bub_dis(g: &mut Usg) -> u64 {
    let n_vtx = (g.a.len() << 1) as u32;
    let mut seen = vec![false; g.a.len()];
    let mut m_len = 0u64;
    for v in 0..n_vtx {
        if seen[(v >> 1) as usize] || g.a[(v >> 1) as usize].deleted {
            continue;
        }
        let mut w = 0u32;
        let arcs: Vec<u32> = g
            .arcs_of(v)
            .iter()
            .filter(|a| !a.deleted)
            .map(|a| a.v ^ 1)
            .collect();
        let kv = arcs.len();
        if kv == 1 {
            w = arcs[0];
            let arcs2: Vec<u32> = g
                .arcs_of(w)
                .iter()
                .filter(|a| !a.deleted)
                .map(|a| a.v ^ 1)
                .collect();
            if arcs2.len() == 1 {
                continue;
            }
        }
        let (c_len, p_bub) = dfs_usg_t_dis(g, &mut seen, v ^ 1);
        if !p_bub {
            continue;
        }
        if c_len > m_len {
            m_len = c_len;
        }
    }
    for slot in seen.iter_mut() {
        *slot = false;
    }
    m_len
}

/// `dfs_usg_t_dis` (`gfa_ut.cpp` 15616's helper): the connected component
/// length with a bubble flag.
fn dfs_usg_t_dis(g: &mut Usg, seen: &mut [bool], x: u32) -> (u64, bool) {
    let mut p_bub = false;
    if seen[(x >> 1) as usize] || g.a[(x >> 1) as usize].deleted {
        return (0, false);
    }
    let mut total = 0u64;
    let mut flag_0 = 0u32;
    let mut flag_1 = 0u32;
    let mut stack: Vec<u32> = vec![x];
    while let Some(c_v) = stack.pop() {
        if seen[(c_v >> 1) as usize] {
            continue;
        }
        seen[(c_v >> 1) as usize] = true;
        let mut convex = u32::MAX;
        let arcs: Vec<u32> = g
            .arcs_of(c_v)
            .iter()
            .filter(|a| !a.deleted)
            .map(|a| a.v)
            .collect();
        if arcs.len() == 1 {
            convex = arcs[0];
        }
        let mut u_len = 0u64;
        if convex != u32::MAX {
            u_len = u64::from(g.a[(convex >> 1) as usize].len);
        }
        total += u_len;
        let mut kv = 0u32;
        for arc in g.arcs_of(c_v ^ 1) {
            if arc.deleted {
                continue;
            }
            kv += 1;
            if !seen[(arc.v >> 1) as usize] {
                stack.push(arc.v);
            }
        }
        let mut kv1 = 0u32;
        if convex != u32::MAX {
            for arc in g.arcs_of(convex ^ 1) {
                if arc.deleted {
                    continue;
                }
                kv1 += 1;
                if !seen[(arc.v >> 1) as usize] {
                    stack.push(arc.v);
                }
            }
        }
        if kv > 0 && kv1 > 0 {
            flag_0 += 1;
        }
        if kv > 1 {
            flag_1 += 1;
        }
        if kv1 > 1 {
            flag_1 += 1;
        }
    }
    if flag_0 > 0 && flag_1 > 1 {
        p_bub = true;
    }
    (total, p_bub)
}

/// `usg_bub_identify` (`gfa_ut.cpp` 15647): the single-sink bubble walk;
/// returns `(sink, body)` when the component closes.
pub fn usg_bub_identify(
    g: &mut Usg,
    v0: u32,
    max_dist: u64,
    max_occ: u64,
) -> Option<(u32, Vec<u32>)> {
    if g.a[(v0 >> 1) as usize].deleted {
        return None;
    }
    if usg_arc_occ(g, v0, None) < 2 {
        return None;
    }
    let mut dist = vec![0u32; g.a.len()];
    let mut parent = vec![u32::MAX; g.a.len()];
    let mut visited = vec![false; g.a.len()];
    let mut remaining = vec![0u32; g.a.len()];
    let mut body: Vec<u32> = Vec::new();
    let mut stack: Vec<u32> = vec![v0];
    let mut n_pending = 0usize;
    let mut is_first = true;
    let mut n_tips = 0u32;
    let mut tip_end = u32::MAX;
    dist[v0 as usize] = 0;
    while !stack.is_empty() {
        let v = stack.pop().unwrap();
        let d = dist[v as usize];
        for arc in g.arcs_of(v) {
            if arc.deleted {
                continue;
            }
            let mut w = arc.v;
            let mut l = (arc.ul & 0xffff_ffff) as u32;
            if (w >> 1) == (v0 >> 1) {
                return None;
            }
            if is_first {
                l = 0;
            }
            if u64::from(d) + u64::from(l) > max_dist {
                break;
            }
            if body.len() as u64 > max_occ {
                break;
            }
            if !visited[w as usize] {
                body.push(w);
                visited[w as usize] = true;
                parent[w as usize] = v;
                dist[w as usize] = d + l;
                remaining[w as usize] = usg_arc_occ(g, w ^ 1, None) as u32;
                n_pending += 1;
            } else if d + l < dist[w as usize] {
                dist[w as usize] = d + l;
            }
            remaining[w as usize] -= 1;
            if remaining[w as usize] == 0 {
                if usg_arc_occ(g, w, None) > 0 {
                    stack.push(w);
                } else {
                    if n_tips != 0 {
                        return None;
                    }
                    n_tips += 1;
                    tip_end = w;
                }
                n_pending -= 1;
            }
        }
        is_first = false;
        if n_tips == 1 && tip_end != u32::MAX && n_pending == 0 && stack.is_empty() {
            return Some((tip_end, body));
        }
    }
    None
}

/// `usg_arc_cut_srt_length` (`gfa_ut.cpp` 11104): the length-sorted weak-arc
/// cut with the double-ended topological check.
pub fn usg_arc_cut_srt_length(
    g: &mut Usg,
    max_ext: i32,
    len_rat: f32,
    is_topo: bool,
    max_drop_len: Option<&mut u32>,
    b: &mut Vec<u64>,
    ub: &mut Vec<u64>,
) -> u32 {
    let n_vtx = (g.a.len() << 1) as u32;
    let mut f = vec![0u8; g.a.len()];
    b.clear();
    ub.clear();
    for v in 0..n_vtx {
        if g.a[(v >> 1) as usize].deleted || f[(v >> 1) as usize] != 0 {
            continue;
        }
        let arcs: Vec<UsgArc> = g.arcs_of(v).to_vec();
        if arcs.len() < 2 {
            continue;
        }
        let mut kv = 0usize;
        let mut kocc = [0usize; 2];
        for arc in &arcs {
            if arc.deleted {
                continue;
            }
            kv += 1;
            if (arc.ou >> 1) > 0 {
                kocc[1] = kocc[1].max((arc.ou >> 1) as usize);
            } else if arc.ou == 0 {
                kocc[0] += 1;
            }
        }
        if kv < 2 || kocc[0] == 0 || kocc[1] == 0 {
            continue;
        }
        for (i, arc) in arcs.iter().enumerate() {
            if arc.deleted || arc.ou != 0 {
                continue;
            }
            if let Some(limit) = max_drop_len.as_deref() {
                if arc.ol >= *limit {
                    continue;
                }
            }
            let x = (u64::from(arc.ol)) << 32;
            b.push(x | ub.len() as u64);
            ub.push(((u64::from(v)) << 32) | i as u64);
        }
    }
    crate::container::ksort::radix_sort_u64(b);
    let mut cnt = 0u32;
    for k in 0..b.len() {
        let slot = (b[k] & 0xffff_ffff) as usize;
        let v = (ub[slot] >> 32) as u32;
        let arc_index = (ub[slot] & 0xffff_ffff) as usize;
        let ve = g.a[(v >> 1) as usize].arc[(v & 1) as usize][arc_index];
        let w = ve.v ^ 1;
        if ve.deleted
            || g.a[(v >> 1) as usize].deleted
            || g.a[(w >> 1) as usize].deleted
            || ve.ou != 0
        {
            continue;
        }
        let nv = g.arcs_of(v).len();
        let nw = g.arcs_of(w).len();
        if nv <= 1 && nw <= 1 {
            continue;
        }
        let Some(we_index) = g.arcs_of(w).iter().position(|a| a.v == v ^ 1) else {
            continue;
        };
        let we = g.a[(w >> 1) as usize].arc[(w & 1) as usize][we_index];
        let mm_ol = ve.ol.min(we.ol);
        let mut kocc = [0usize; 2];
        let mut kv = 0usize;
        let mut ol_max = 0u32;
        for arc in g.arcs_of(v) {
            if arc.deleted {
                continue;
            }
            kv += 1;
            if arc.ou != 1 {
                kocc[usize::from(arc.ou != 0)] += 1;
            }
            if ol_max < arc.ol {
                ol_max = arc.ol;
            }
        }
        if kv < 1 || kocc[0] < 1 || kocc[1] < 1 {
            continue;
        }
        if kv >= 2 && (mm_ol as f32) > ol_max as f32 * len_rat {
            continue;
        }
        let mut kw = 0usize;
        let mut ol_max = 0u32;
        for arc in g.arcs_of(w) {
            if arc.deleted {
                continue;
            }
            kw += 1;
            if ol_max < arc.ol {
                ol_max = arc.ol;
            }
        }
        if kw < 1 {
            continue;
        }
        if kw >= 2 && (mm_ol as f32) > ol_max as f32 * len_rat {
            continue;
        }
        if kv <= 1 && kw <= 1 {
            continue;
        }
        let mut to_del = true;
        if is_topo {
            to_del = false;
            if kv > 1 && kw > 1 {
                to_del = true;
            } else if kw == 1 {
                let mut z_a: Vec<u32> = Vec::new();
                let mut z_b: Vec<u32> = Vec::new();
                if usg_naive_topocut_aux(g, w ^ 1, max_ext, &mut f, &mut z_a, &mut z_b) < max_ext {
                    to_del = true;
                }
            } else if kv == 1 {
                let mut z_a: Vec<u32> = Vec::new();
                let mut z_b: Vec<u32> = Vec::new();
                if usg_naive_topocut_aux(g, v ^ 1, max_ext, &mut f, &mut z_a, &mut z_b) < max_ext {
                    to_del = true;
                }
            }
        }
        if to_del {
            set_arc_deleted(g, v, ve.v, true);
            set_arc_deleted(g, w, v ^ 1, true);
            if usg_naive_topocut_aux_sec(g, v, max_ext) < max_ext
                && usg_naive_topocut_aux_sec(g, w, max_ext) < max_ext
            {
                cnt += 1;
            } else {
                set_arc_deleted(g, v, ve.v, false);
                set_arc_deleted(g, w, v ^ 1, false);
            }
        }
    }
    if cnt > 0 {
        usg_cleanup(g);
    }
    cnt
}

fn set_arc_deleted(g: &mut Usg, v: u32, w: u32, del: bool) {
    if let Some(arc) = g
        .a[(v >> 1) as usize]
        .arc[(v & 1) as usize]
        .iter_mut()
        .find(|a| a.v == w)
    {
        arc.deleted = del;
    }
}

/// `usg_naive_topocut_aux_sec` (`gfa_ut.cpp` 10902): the both-direction walk.
pub fn usg_naive_topocut_aux_sec(g: &Usg, v0: u32, max_ext: i32) -> i32 {
    let mut n_ext = 0i32;
    let mut tip = [false; 2];
    let mut v = v0;
    loop {
        let arcs: Vec<u32> = g
            .arcs_of(v)
            .iter()
            .filter(|a| !a.deleted)
            .map(|a| a.v)
            .collect();
        let kv = arcs.len();
        n_ext += g.a[(v >> 1) as usize].occ as i32;
        if kv != 1 {
            if kv == 0 {
                tip[0] = true;
            }
            break;
        }
        v = arcs[0];
        let kv = g.arcs_of(v ^ 1).iter().filter(|a| !a.deleted).count();
        if kv != 1 {
            break;
        }
        if v == v0 {
            return 0; // circle, it is ok to remove it
        }
    }
    let mut v = v0 ^ 1;
    loop {
        let arcs: Vec<u32> = g
            .arcs_of(v)
            .iter()
            .filter(|a| !a.deleted)
            .map(|a| a.v)
            .collect();
        let kv = arcs.len();
        n_ext += g.a[(v >> 1) as usize].occ as i32;
        if kv != 1 {
            if kv == 0 {
                tip[1] = true;
            }
            break;
        }
        v = arcs[0];
        let kv = g.arcs_of(v ^ 1).iter().filter(|a| !a.deleted).count();
        if kv != 1 {
            break;
        }
        if v == v0 ^ 1 {
            return 0;
        }
    }
    if n_ext < max_ext && tip[0] && tip[1] {
        return n_ext;
    }
    max_ext
}

/// `get_ul_path_info` (`gfa_ut.cpp` 9004): walks the linear path from `s`,
/// returning the terminal kind; records the visited vertices and the
/// UL/bridge occupancies.
#[allow(clippy::too_many_arguments)]
pub fn get_ul_path_info(
    g: &Asg,
    cov: &crate::assembly::ul2ul::UlCov,
    bg: Option<&crate::assembly::ul2ul::UlBg>,
    s: u32,
    mut e: Option<&mut u32>,
    mut occ: Option<&mut u32>,
    mut ul_cnt: Option<&mut u32>,
    mut bridge_w: Option<&mut u32>,
    mut bridge_am: Option<&mut u32>,
    b: &mut Vec<u64>,
) -> u32 {
    let mut v = s;
    let mut uv = u32::MAX;
    if let Some(out) = occ.as_deref_mut() {
        *out = 0;
    }
    if let Some(out) = ul_cnt.as_deref_mut() {
        *out = 0;
    }
    if let Some(out) = bridge_w.as_deref_mut() {
        *out = 0;
    }
    if let Some(out) = bridge_am.as_deref_mut() {
        *out = 0;
    }
    loop {
        if let Some(out) = occ.as_deref_mut() {
            *out += 1;
        }
        let mut arc_idx: Vec<usize> = Vec::new();
        let kv = get_arcs(g, v, Some(&mut arc_idx), 1);
        if let Some(out) = e.as_deref_mut() {
            *out = v;
        }
        b.push(u64::from(v));
        if let Some(out) = ul_cnt.as_deref_mut() {
            *out += cov.uc[(v >> 1) as usize];
        }
        if bridge_w.as_deref_mut().is_some() || bridge_am.as_deref_mut().is_some() {
            let uw = v;
            if let Some(out) = bridge_w.as_deref_mut() {
                *out += bg.map_or(0, |t| t.w_n.get((v >> 1) as usize).copied().unwrap_or(0));
            }
            if let Some(out) = bridge_am.as_deref_mut() {
                *out += bg.map_or(0, |t| t.a_n.get((v >> 1) as usize).copied().unwrap_or(0));
            }
            if uv != u32::MAX {
                if let Some(t) = bg {
                    // The reference compares the two boundary raw vertices of
                    // the previous and current integer unitigs through
                    // `cc.iug_a`; the bridge table stores the flag on the
                    // same pair.
                    let ia = &cov.iug_a[(uv >> 1) as usize];
                    let mut bv = ia[if uv & 1 != 0 { 0 } else { ia.len() - 1 }].v;
                    if uv & 1 != 0 {
                        bv ^= 1;
                    }
                    let iw = &cov.iug_a[(uw >> 1) as usize];
                    let mut bw = iw[if uw & 1 != 0 { iw.len() - 1 } else { 0 }].v;
                    if uw & 1 != 0 {
                        bw ^= 1;
                    }
                    if crate::assembly::ul2ul::ulg_type(cov.uln, (bv >> 1) as usize) != 0
                        && crate::assembly::ul2ul::ulg_type(cov.uln, (bw >> 1) as usize) != 0
                    {
                        let l = cov.get_bg_flag_pub(&t.bg, bv, bw);
                        if l == BG_WRONG {
                            if let Some(out) = bridge_w.as_deref_mut() {
                                *out += 1;
                            }
                        }
                        if l == BG_AMBIGUOUS {
                            if let Some(out) = bridge_am.as_deref_mut() {
                                *out += 1;
                            }
                        }
                    }
                }
            }
            uv = uw;
        }
        if kv == 0 {
            return u32::from(crate::assembly::gfa_clean::END_TIPS);
        }
        if kv == 2 {
            return 2; // TWO_OUTPUT
        }
        if kv > 2 {
            return 4; // MUL_OUTPUT
        }
        let w = g.arcs_of(v)[arc_idx[0]].v;
        let kw = get_arcs(g, w ^ 1, None, 1);
        v = w;
        if kw == 2 {
            return 1; // TWO_INPUT
        }
        if kw > 2 {
            return 3; // MUL_INPUT
        }
        if v == s {
            return u32::from(crate::assembly::gfa_clean::LOOP);
        }
    }
}

/// `is_het_bridge` (`gfa_ut.cpp` 9554): whether the boundary node at
/// `match_bound` sits inside a het region.
pub fn is_het_bridge(bub: &BubbleType, p_a: &[u64], match_bound: i64) -> bool {
    let mut z = match_bound;
    while z >= 0 && crate::assembly::ul_renew_gfa::if_hom((p_a[z as usize] >> 1) as u32, bub) {
        z -= 1;
    }
    if z < 0 {
        return false;
    }
    let mut z = match_bound + 1;
    while (z as usize) < p_a.len()
        && crate::assembly::ul_renew_gfa::if_hom((p_a[z as usize] >> 1) as u32, bub)
    {
        z += 1;
    }
    (z as usize) < p_a.len()
}

/// `idx_check_rate` (`gfa_ut.cpp` 10228): the sliding-window match-rate
/// check of the ordered positions.
pub fn idx_check_rate(a: &[u64], tot: u64, winlen: u64, match_rate: f32) -> bool {
    if a.is_empty() {
        return false;
    }
    let mut winlen = winlen;
    if tot < winlen {
        winlen = tot;
    }
    let mut s = 0u32;
    let mut e = winlen as u32;
    let mut k = 0usize;
    while k < a.len() {
        debug_assert_ne!(a[k], u64::MAX);
        if (a[k] as u32) >= s && (a[k] as u32) < e {
            k += 1;
            continue;
        }
        break;
    }
    let mut p = k as i64;
    if p > 0 && p >= (winlen as f32 * match_rate) as i64 {
        return true;
    }
    let mut ks = 0usize;
    while k < a.len() {
        debug_assert_ne!(a[k], u64::MAX);
        e = (a[k] as u32) + 1;
        p += 1;
        while ks < k {
            if i64::from(e) - i64::from(a[ks] as u32) <= i64::from(winlen as u32) {
                break;
            }
            p -= 1;
            ks += 1;
        }
        debug_assert!(p >= 0);
        if p > 0 && p >= (winlen as f32 * match_rate) as i64 {
            return true;
        }
        k += 1;
    }
    false
}

/// `gen_ug_integer_seq_on_fly` (`gfa_ut.cpp` 9203): the raw-read sequence of
/// an integer-unitig path, appending to `res`.
pub fn gen_ug_integer_seq_on_fly(
    cov: &crate::assembly::ul2ul::UlCov,
    z_items: &crate::assembly::ul2ul::Ul2ulIdx,
    raw: &crate::assembly::unitig::UnitigGraph,
    u_a: &[u64],
    res: &mut Vec<u64>,
) {
    let uln = z_items.uln as usize;
    for k in 0..u_a.len() {
        let seq = &cov.iug_a[(u_a[k] >> 1) as usize];
        let u_rev = u_a[k] & 1 != 0;
        if k > 0 {
            let uv = u_a[k - 1];
            let uw = u_a[k];
            let mut bv = endpoint_of(z_items, iug_end(&cov.iug_a, uv));
            if uv & 1 != 0 {
                bv ^= 1;
            }
            let mut bw = endpoint_of(z_items, iug_end(&cov.iug_a, uw));
            if uw & 1 != 0 {
                bw ^= 1;
            }
            let is_bw_ul = ulg_type_of(z_items, bw);
            if is_bw_ul {
                let m = crate::assembly::ul2ul::get_ul_spec_ovlp(
                    z_items,
                    (bv >> 1) as usize,
                    (bw >> 1) as usize,
                )
                .expect("the string link must exist");
                assert!(!m.is_del);
                let nn = (m.te_k - m.ts_k) as usize;
                let zn = seq.len();
                if !u_rev {
                    for z in 0..nn {
                        debug_assert_eq!(res[res.len() - nn + z], u64::from(seq[z].v));
                    }
                    for z in nn..zn {
                        res.push(u64::from(seq[z].v));
                    }
                } else {
                    for z in (zn - nn..zn).rev() {
                        debug_assert_eq!(res[res.len() - nn + (zn - 1 - z)], u64::from(seq[z].v ^ 1));
                    }
                    for z in (0..zn - nn).rev() {
                        res.push(u64::from(seq[z].v ^ 1));
                    }
                }
            } else {
                let mut x = bw >> 1;
                x <<= 1;
                x += bw & 1;
                let is_bv_ul = ulg_type_of(z_items, bv);
                if !is_bv_ul {
                    assert!(
                        raw.graph
                            .arcs_of(res[res.len() - 1] as u32)
                            .iter()
                            .any(|a| !a.deleted && a.v == x)
                    );
                } else if res[res.len() - 1] == u64::from(x) {
                    res.pop();
                } else {
                    assert!(
                        raw.graph
                            .arcs_of(res[res.len() - 1] as u32)
                            .iter()
                            .any(|a| !a.deleted && a.v == x)
                    );
                }
                if !u_rev {
                    for entry in seq {
                        res.push(u64::from(entry.v));
                    }
                } else {
                    for entry in seq.iter().rev() {
                        res.push(u64::from(entry.v ^ 1));
                    }
                }
            }
        } else if !u_rev {
            for entry in seq {
                res.push(u64::from(entry.v));
            }
        } else {
            for entry in seq.iter().rev() {
                res.push(u64::from(entry.v ^ 1));
            }
        }
    }
}

/// The raw unitig id of an integer-unitig path's boundary vertex.
fn iug_end(iug_a: &[Vec<crate::assembly::ul2ul::UinfoSrt>], v: u64) -> u32 {
    let walk = &iug_a[(v >> 1) as usize];
    if v & 1 != 0 {
        walk[0].v
    } else {
        walk[walk.len() - 1].v
    }
}

fn endpoint_of(z_items: &crate::assembly::ul2ul::Ul2ulIdx, raw_vertex: u32) -> u32 {
    let _ = z_items;
    raw_vertex
}

fn ulg_type_of(z_items: &crate::assembly::ul2ul::Ul2ulIdx, raw_vertex: u32) -> bool {
    let raw_id = (raw_vertex >> 1) as usize;
    raw_id < z_items.uln as usize
}

/// `get_integer_seq_ovlps` (`gfa_ut.cpp` 9301): the string-supported
/// overhangs around `match_bound`; returns the support count.
#[allow(clippy::too_many_arguments)]
pub fn get_integer_seq_ovlps(
    str_idx: &crate::assembly::ul_integer::UlStrIdx,
    store: &[crate::assembly::ul_store::UlReadVec],
    raw_unitigs: &[crate::assembly::asg::MaUtg],
    read_lengths: &[u64],
    unitig_reads: &[Vec<u64>],
    bub: &BubbleType,
    p_a: &[u64],
    match_bound: i64,
    skip_hom: bool,
    mut res: Option<&mut Vec<u64>>,
    r_w: &mut u64,
) -> u32 {
    let v = p_a[match_bound as usize] as u32;
    let vertex = (v >> 1) as usize;
    let start = str_idx.idx[vertex] as usize;
    let end = str_idx.idx[vertex + 1] as usize;
    let mut ww = [0u64; 2];
    let mut occ = 0u32;
    *r_w = 0;
    for z in start..end {
        let word = str_idx.occ[z];
        let read_id = (word >> 32) as usize;
        let pos = (word & 0xffff_ffff) as usize;
        let str_ = &str_idx.str[read_id];
        let s_n = str_.cn as usize;
        if s_n < 2 {
            continue;
        }
        let vz = str_.a[pos] as u32;
        debug_assert_eq!(v >> 1, vz >> 1);
        let mut sw = [0u64; 2];
        let mut ps = u64::MAX;
        let mut pe = u64::MAX;
        if v == vz {
            let mut s = pos as i64 + 1;
            let mut p = match_bound + 1;
            while s < s_n as i64
                && (p as usize) < p_a.len()
                && str_.a[s as usize] as u32 == p_a[p as usize] as u32
            {
                let xi = &store[read_id].bb[str_.a[s as usize] as usize >> 32];
                debug_assert_eq!(
                    (xi.hid << 1) + u32::from(xi.rev),
                    str_.a[s as usize] as u32
                );
                if !(skip_hom && crate::assembly::ul_renew_gfa::if_hom(xi.hid, bub)) {
                    sw[1] += crate::assembly::ul_integer::ug_occ_w(
                        u64::from(xi.ts),
                        u64::from(xi.te),
                        read_lengths,
                        &unitig_reads[xi.hid as usize],
                    );
                }
                s += 1;
                p += 1;
            }
            if s < s_n as i64 && (p as usize) < p_a.len() {
                continue;
            }
            if p <= match_bound + 1 {
                continue;
            }
            if sw[1] == 0 {
                continue;
            }
            pe = p as u64;
            let mut s = pos as i64;
            let mut p = match_bound;
            while s >= 0 && p >= 0 && str_.a[s as usize] as u32 == p_a[p as usize] as u32 {
                let xi = &store[read_id].bb[str_.a[s as usize] as usize >> 32];
                debug_assert_eq!(
                    (xi.hid << 1) + u32::from(xi.rev),
                    str_.a[s as usize] as u32
                );
                if !(skip_hom && crate::assembly::ul_renew_gfa::if_hom(xi.hid, bub)) {
                    sw[0] += crate::assembly::ul_integer::ug_occ_w(
                        u64::from(xi.ts),
                        u64::from(xi.te),
                        read_lengths,
                        &unitig_reads[xi.hid as usize],
                    );
                }
                s -= 1;
                p -= 1;
            }
            debug_assert!(p < match_bound);
            if s >= 0 && p >= 0 {
                continue;
            }
            if sw[0] == 0 {
                continue;
            }
            ps = (p + 1) as u64;
        } else {
            let mut s = pos as i64 - 1;
            let mut p = match_bound + 1;
            while s >= 0
                && (p as usize) < p_a.len()
                && str_.a[s as usize] as u32 == (p_a[p as usize] as u32 ^ 1)
            {
                let xi = &store[read_id].bb[str_.a[s as usize] as usize >> 32];
                debug_assert_eq!(
                    (xi.hid << 1) + u32::from(xi.rev),
                    str_.a[s as usize] as u32
                );
                if !(skip_hom && crate::assembly::ul_renew_gfa::if_hom(xi.hid, bub)) {
                    sw[1] += crate::assembly::ul_integer::ug_occ_w(
                        u64::from(xi.ts),
                        u64::from(xi.te),
                        read_lengths,
                        &unitig_reads[xi.hid as usize],
                    );
                }
                s -= 1;
                p += 1;
            }
            if s >= 0 && (p as usize) < p_a.len() {
                continue;
            }
            if p <= match_bound + 1 {
                continue;
            }
            if sw[1] == 0 {
                continue;
            }
            pe = p as u64;
            let mut s = pos as i64;
            let mut p = match_bound;
            while s < s_n as i64 && p >= 0 && str_.a[s as usize] as u32 == (p_a[p as usize] as u32 ^ 1)
            {
                let xi = &store[read_id].bb[str_.a[s as usize] as usize >> 32];
                debug_assert_eq!(
                    (xi.hid << 1) + u32::from(xi.rev),
                    str_.a[s as usize] as u32
                );
                if !(skip_hom && crate::assembly::ul_renew_gfa::if_hom(xi.hid, bub)) {
                    sw[0] += crate::assembly::ul_integer::ug_occ_w(
                        u64::from(xi.ts),
                        u64::from(xi.te),
                        read_lengths,
                        &unitig_reads[xi.hid as usize],
                    );
                }
                s += 1;
                p -= 1;
            }
            debug_assert!(p < match_bound);
            if s < s_n as i64 && p >= 0 {
                continue;
            }
            if sw[0] == 0 {
                continue;
            }
            ps = (p + 1) as u64;
        }
        if let Some(out) = res.as_deref_mut() {
            out.push((ps << 32) | pe);
        }
        occ += 1;
        ww[0] += sw[0];
        ww[1] += sw[1];
        let _ = raw_unitigs;
    }
    *r_w = ww[0].min(ww[1]);
    occ
}

/// The `ul_resolve_t` slices shared by the ulg cleanup cluster. The C
/// driver threads `uidx` through every helper; the port passes this context
/// instead.
pub struct UlResolve<'a> {
    /// `uidx->uovl.cc`.
    pub cov: &'a mut UlCov,
    /// `uidx->uovl`, the ul2ul index.
    pub z_items: &'a Ul2ulIdx,
    /// `uidx->uovl.i_ug`, the integer unitig graph being cleaned.
    pub iug: &'a mut crate::assembly::unitig::UnitigGraph,
    /// `uidx->l1_ug`, the raw unitig graph.
    pub raw: &'a crate::assembly::unitig::UnitigGraph,
    /// The UL read store (`uidx->UL_INF.a`).
    pub store: &'a [crate::assembly::ul_store::UlReadVec],
    /// `raw->u.a`.
    pub raw_unitigs: &'a [MaUtg],
    /// `R_INF.read_lengths`.
    pub read_lengths: &'a [u64],
    /// The per-unitig read lists consumed by `get_integer_seq_ovlps`.
    pub unitig_reads: &'a [Vec<u64>],
    /// `uidx->pstr`.
    pub str_idx: &'a crate::assembly::ul_integer::UlStrIdx,
    /// `uidx->bub`.
    pub bub: &'a BubbleType,
    /// `uidx->uovl.bg`.
    pub bg: &'a UlBg,
    /// `uidx->uovl.telo` (`NULL` without `--telo-motif`).
    pub telo: Option<&'a [u8]>,
}

impl UlResolve<'_> {
    /// `get_iug_u_raw_occ` (`gfa_ut.cpp` 8041).
    pub fn get_iug_u_raw_occ(
        &self,
        id: u32,
        ul_occ: Option<&mut u32>,
        raw_ug_occ: Option<&mut u32>,
    ) {
        if let Some(out) = ul_occ {
            *out = self.cov.raw_uc[id as usize];
        }
        if let Some(out) = raw_ug_occ {
            *out = self.iug.unitigs[id as usize].reads.len() as u32
                - self.cov.raw_uc[id as usize];
        }
    }

    /// `get_bridges` (`gfa_ut.cpp` 8636): the wrong/ambiguous bridge counts
    /// along a vertex list.
    pub fn get_bridges(&self, v_a: &[u64], w_occ: Option<&mut u32>, am_occ: Option<&mut u32>) {
        let mut w = 0u64;
        let mut am = 0u64;
        let mut uv = u32::MAX;
        for &k in v_a {
            w += u64::from(self.bg.w_n[(k >> 1) as usize]);
            am += u64::from(self.bg.a_n[(k >> 1) as usize]);
            let uw = k as u32;
            if uv != u32::MAX {
                let ia = &self.cov.iug_a[(uv >> 1) as usize];
                let mut bv = ia[if uv & 1 != 0 { 0 } else { ia.len() - 1 }].v;
                if uv & 1 != 0 {
                    bv ^= 1;
                }
                let iw = &self.cov.iug_a[(uw >> 1) as usize];
                let mut bw = iw[if uw & 1 != 0 { iw.len() - 1 } else { 0 }].v;
                if uw & 1 != 0 {
                    bw ^= 1;
                }
                if crate::assembly::ul2ul::ulg_type(self.cov.uln, (bv >> 1) as usize) != 0
                    && crate::assembly::ul2ul::ulg_type(self.cov.uln, (bw >> 1) as usize) != 0
                {
                    let l = self.cov.get_bg_flag_pub(&self.bg.bg, bv, bw);
                    if l == BG_WRONG {
                        w += 1;
                    }
                    if l == BG_AMBIGUOUS {
                        am += 1;
                    }
                }
            }
            uv = uw;
        }
        if let Some(out) = w_occ {
            *out = w as u32;
        }
        if let Some(out) = am_occ {
            *out = am as u32;
        }
    }
}

/// `get_remove_hifi_occ` (`gfa_ut.cpp` 8499): masks the raw-read words of
/// the candidate vertices in the CSR and reports whether the deletion is
/// acceptable — the reference returns 1 when `del_n <= thres`. `del_occ`
/// carries the deleted base count.
pub fn get_remove_hifi_occ(
    ctx: &mut UlResolve,
    thres: u64,
    v_a: &[u64],
    ub: &mut Vec<u64>,
    del_occ: Option<&mut u32>,
) -> bool {
    let cov = &mut *ctx.cov;
    let iug_g = &ctx.iug.graph;
    let mut del_n: u64 = 0;
    let mut buf: Vec<u64> = Vec::new();
    let mut fn_count = 0usize;
    let mut dup_raw = false;
    for &vk in v_a {
        let iu = cov.iug_a[(vk >> 1) as usize].clone();
        for (z, item) in iu.iter().enumerate() {
            let raw_id = (item.v >> 1) as usize;
            let start = cov.iug_idx[raw_id] as usize;
            let end = cov.iug_idx[raw_id + 1] as usize;
            let mut keep_ns = [0u32; 2];
            let pn = buf.len();
            let mut n_mask = 0usize;
            for m in start..end {
                let word = cov.iug_b[m];
                let iug_id = ((word << 1) >> 33) as usize;
                let iug_off = (word & 0x7fff_ffff) as usize;
                if iug_g.sequences[iug_id].deleted {
                    continue;
                }
                if (word >> 63) != 0 {
                    n_mask += 1;
                    continue;
                }
                if iug_id == (vk >> 1) as usize && iug_off == z {
                    cov.iug_b[m] |= 0x8000_0000_0000_0000;
                    buf.push((u64::MAX << 32) | m as u64);
                } else {
                    let entry = &cov.iug_a[iug_id][iug_off];
                    if entry.v & 1 == 0 {
                        keep_ns[0] = keep_ns[0].max(entry.n);
                    } else {
                        keep_ns[1] = keep_ns[1].max(entry.n);
                    }
                }
            }
            debug_assert_eq!(buf.len(), pn + 1);
            let raw_len = cov.raw_len(raw_id);
            if u64::from(keep_ns[0]) + u64::from(keep_ns[1]) >= raw_len {
                keep_ns[0] = raw_len as u32;
                keep_ns[1] = raw_len as u32;
            }
            let is = u64::from(keep_ns[0]);
            let ie = raw_len - u64::from(keep_ns[1]);
            if ie > is {
                let (ms, me) = if item.v & 1 == 0 {
                    (0u64, u64::from(item.n))
                } else {
                    (raw_len - u64::from(item.n), raw_len)
                };
                if is.max(ms) < ie.min(me) {
                    del_n += ie.min(me) - is.max(ms);
                    buf[pn] = ((buf[pn] << 32) >> 32) | (raw_id as u64) << 32;
                    fn_count += 1;
                    if n_mask > 0 {
                        dup_raw = true;
                    }
                }
            }
        }
    }
    if del_n > thres && dup_raw {
        radix_sort_u64(&mut buf);
        del_n = 0;
        let mut l = 0usize;
        let mut k = 1usize;
        while k <= fn_count {
            if k == fn_count || (buf[k] >> 32) != (buf[l] >> 32) {
                let raw_id = (buf[l] >> 32) as usize;
                let start = cov.iug_idx[raw_id] as usize;
                let end = cov.iug_idx[raw_id + 1] as usize;
                let mut del_ns = [0u32; 2];
                let mut keep_ns = [0u32; 2];
                for m in start..end {
                    let word = cov.iug_b[m];
                    let iug_id = ((word << 1) >> 33) as usize;
                    let iug_off = (word & 0x7fff_ffff) as usize;
                    if iug_g.sequences[iug_id].deleted {
                        continue;
                    }
                    let entry = &cov.iug_a[iug_id][iug_off];
                    if (word >> 63) != 0 {
                        if entry.v & 1 == 0 {
                            del_ns[0] = del_ns[0].max(entry.n);
                        } else {
                            del_ns[1] = del_ns[1].max(entry.n);
                        }
                    } else if entry.v & 1 == 0 {
                        keep_ns[0] = keep_ns[0].max(entry.n);
                    } else {
                        keep_ns[1] = keep_ns[1].max(entry.n);
                    }
                }
                debug_assert!(del_ns[0] + del_ns[1] > 0);
                let raw_len = cov.raw_len(raw_id);
                if u64::from(keep_ns[0]) + u64::from(keep_ns[1]) >= raw_len {
                    keep_ns[0] = raw_len as u32;
                    keep_ns[1] = raw_len as u32;
                }
                let is = u64::from(keep_ns[0]);
                let ie = raw_len - u64::from(keep_ns[1]);
                if ie > is {
                    if u64::from(del_ns[0]) + u64::from(del_ns[1]) >= raw_len {
                        del_n += ie - is;
                    } else {
                        let mut ms = 0u64;
                        let me = u64::from(del_ns[0]);
                        if ms < me && is.max(ms) < ie.min(me) {
                            del_n += ie.min(me) - is.max(ms);
                        }
                        ms = raw_len - u64::from(del_ns[1]);
                        let me = raw_len;
                        if ms < me && is.max(ms) < ie.min(me) {
                            del_n += ie.min(me) - is.max(ms);
                        }
                    }
                }
                if del_n > thres {
                    break;
                }
                l = k;
            }
            k += 1;
        }
    }
    for &word in &buf {
        cov.iug_b[(word & 0xffff_ffff) as usize] <<= 1;
        cov.iug_b[(word & 0xffff_ffff) as usize] >>= 1;
    }
    if let Some(out) = del_occ {
        *out = del_n as u32;
    }
    ub.clear();
    del_n <= thres
}

/// `usg_topocut_aux` (`gfa_ut.cpp` 8784): the naive topological cut on the
/// integer unitig graph.
pub fn usg_topocut_aux(
    ctx: &mut UlResolve,
    v: u32,
    max_ext: u32,
    max_ext_hifi: u32,
    b: &mut Vec<u64>,
    ub: &mut Vec<u64>,
) -> bool {
    let mut v = v;
    let pn = b.len();
    let mut n_ext = 0u32;
    while n_ext < max_ext {
        let kv = get_arcs(&ctx.iug.graph, v ^ 1, None, 1);
        if kv != 1 {
            break;
        }
        let mut ul_occ = 0u32;
        ctx.get_iug_u_raw_occ(v >> 1, Some(&mut ul_occ), None);
        n_ext += ul_occ;
        b.push(u64::from(v));
        let av: Vec<u32> = ctx
            .iug
            .graph
            .arcs_of(v)
            .iter()
            .filter(|a| !a.deleted)
            .map(|a| a.v)
            .collect();
        if av.len() != 1 {
            break;
        }
        v = av[0];
    }
    if n_ext < max_ext {
        let tail: Vec<u64> = b[pn..].to_vec();
        if get_remove_hifi_occ(ctx, u64::from(max_ext_hifi), &tail, ub, None) {
            b.truncate(pn);
            return true;
        }
    }
    b.truncate(pn);
    false
}

/// The shared walk of `get_ul_arc_supports` / `check_ul_contain_arc_supports`
/// (`gfa_ut.cpp` 9577 / 9677): the query/target support split of one arc
/// against its siblings; returns the collapse flag. The two reference
/// variants differ only in zeroing `min_w_v` on collapse.
#[allow(clippy::too_many_arguments)]
fn ul_arc_supports_walk(
    ctx: &mut UlResolve,
    ve: &AsgArc,
    b_int: &mut Vec<u64>,
    b_raw: &mut Vec<u64>,
    skip_hom: bool,
    return_w_v: &mut u64,
    return_w_r: &mut u64,
    zero_on_collapse: bool,
) -> u32 {
    let cov = &mut *ctx.cov;
    let iug = &*ctx.iug;
    let b_int_s = b_int.len();
    let b_raw_s = b_raw.len();
    *return_w_v = u64::MAX;
    *return_w_r = u64::MAX;
    let mut is_collapse = false;
    let v = (ve.ul >> 32) as u32;
    get_ul_path_info(&iug.graph, cov, None, v ^ 1, None, None, None, None, None, b_int);
    let nv = b_int.len() - b_int_s;
    for k in 0..nv / 2 {
        let front = b_int[b_int_s + k];
        let back = b_int[b_int_s + nv - k - 1];
        b_int[b_int_s + k] = back ^ 1;
        b_int[b_int_s + nv - k - 1] = front ^ 1;
    }
    if nv & 1 != 0 {
        b_int[b_int_s + nv / 2] ^= 1;
    }
    // `n_pre` is captured right after the reversed `v^1` path: the per-sibling
    // reset keeps that prefix and drops the `ve->v` path with it.
    let n_pre = b_int.len();
    let mut v_occ = 0u32;
    get_ul_path_info(
        &iug.graph, cov, None, ve.v, None, Some(&mut v_occ), None, None, None, b_int,
    );
    let v_path: Vec<u64> = b_int[b_int_s..].to_vec();
    gen_ug_integer_seq_on_fly(cov, ctx.z_items, ctx.raw, &v_path, b_raw);
    let l_v = b_raw.len() - b_raw_s;
    let av: Vec<AsgArc> = iug.graph.arcs_of(v).to_vec();
    let mut min_w_v = u64::MAX;
    let mut min_w_r = u64::MAX;
    for arc in &av {
        if arc.deleted || arc.v == ve.v {
            continue;
        }
        b_int.truncate(n_pre);
        b_raw.truncate(b_raw_s + l_v);
        let mut r_occ = 0u32;
        get_ul_path_info(
            &iug.graph, cov, None, arc.v, None, Some(&mut r_occ), None, None, None, b_int,
        );
        let r_path: Vec<u64> = b_int[b_int_s..].to_vec();
        gen_ug_integer_seq_on_fly(cov, ctx.z_items, ctx.raw, &r_path, b_raw);
        let l_r = b_raw.len() - b_raw_s - l_v;
        let zn = l_v.min(l_r);
        let mut z = 0usize;
        while z < zn && b_raw[b_raw_s + z] == b_raw[b_raw_s + l_v + z] {
            z += 1;
        }
        let mut skip_hom_local = skip_hom;
        if skip_hom_local && z < l_v {
            skip_hom_local = is_het_bridge(ctx.bub, &b_raw[b_raw_s..b_raw_s + l_v], z as i64 - 1);
        }
        if skip_hom_local && z < l_r {
            skip_hom_local =
                is_het_bridge(ctx.bub, &b_raw[b_raw_s + l_v..b_raw_s + l_v + l_r], z as i64 - 1);
        }
        assert!(z > 0);
        let mut w_v = u64::MAX;
        let mut w_r = u64::MAX;
        if z < l_v {
            get_integer_seq_ovlps(
                ctx.str_idx,
                ctx.store,
                ctx.raw_unitigs,
                ctx.read_lengths,
                ctx.unitig_reads,
                ctx.bub,
                &b_raw[b_raw_s..b_raw_s + l_v],
                z as i64 - 1,
                skip_hom_local,
                None,
                &mut w_v,
            );
        }
        if z < l_r {
            get_integer_seq_ovlps(
                ctx.str_idx,
                ctx.store,
                ctx.raw_unitigs,
                ctx.read_lengths,
                ctx.unitig_reads,
                ctx.bub,
                &b_raw[b_raw_s + l_v..b_raw_s + l_v + l_r],
                z as i64 - 1,
                skip_hom_local,
                None,
                &mut w_r,
            );
        }
        if w_v == u64::MAX {
            w_v = 0;
        }
        if w_r == u64::MAX {
            w_r = 0;
        }
        if min_w_v == u64::MAX || z == zn || min_w_v > w_v || (min_w_v == w_v && min_w_r < w_r) {
            min_w_v = w_v;
            min_w_r = w_r;
            if z == zn && v_occ < r_occ {
                is_collapse = true;
            }
        }
    }
    if zero_on_collapse && is_collapse {
        min_w_v = 0;
    }
    *return_w_v = min_w_v;
    *return_w_r = min_w_r;
    b_int.truncate(b_int_s);
    b_raw.truncate(b_raw_s);
    if is_collapse {
        1
    } else {
        0
    }
}

/// `get_ul_arc_supports` (`gfa_ut.cpp` 9577).
#[allow(clippy::too_many_arguments)]
pub fn get_ul_arc_supports(
    ctx: &mut UlResolve,
    ve: &AsgArc,
    b_int: &mut Vec<u64>,
    b_raw: &mut Vec<u64>,
    skip_hom: bool,
    return_w_v: &mut u64,
    return_w_r: &mut u64,
) -> u32 {
    ul_arc_supports_walk(
        ctx, ve, b_int, b_raw, skip_hom, return_w_v, return_w_r, true,
    )
}

/// `check_ul_contain_arc_supports` (`gfa_ut.cpp` 9677): the identical walk
/// without the collapse zeroing of `min_w_v`.
#[allow(clippy::too_many_arguments)]
pub fn check_ul_contain_arc_supports(
    ctx: &mut UlResolve,
    ve: &AsgArc,
    b_int: &mut Vec<u64>,
    b_raw: &mut Vec<u64>,
    skip_hom: bool,
    return_w_v: &mut u64,
    return_w_r: &mut u64,
) -> u32 {
    ul_arc_supports_walk(
        ctx, ve, b_int, b_raw, skip_hom, return_w_v, return_w_r, false,
    )
}

/// `check_ulg_to_del` (`gfa_ut.cpp` 9735).
#[allow(clippy::too_many_arguments)]
pub fn check_ulg_to_del(
    ctx: &mut UlResolve,
    v: u32,
    w: u32,
    kv: u32,
    kw: u32,
    max_ext: u32,
    max_ext_hifi: u64,
    mut topo_level: u32,
    collapse: bool,
    b: &mut Vec<u64>,
    ub: &mut Vec<u64>,
) -> bool {
    let mut to_del = false;
    if collapse {
        topo_level = 3;
    }
    if topo_level == 0 {
        to_del = true;
    } else if topo_level == 2 {
        if kv > 1 && kw > 1 {
            to_del = true;
        }
    } else if kv > 1 && kw > 1 {
        to_del = true;
    } else if kw == 1 {
        if usg_topocut_aux(ctx, w ^ 1, max_ext, max_ext_hifi as u32, b, ub) {
            to_del = true;
        }
    } else if kv == 1
        && usg_topocut_aux(ctx, v ^ 1, max_ext, max_ext_hifi as u32, b, ub)
    {
        to_del = true;
    }
    to_del
}

/// `ulg_arc_cut_supports` (`gfa_ut.cpp` 9757): the support-driven weak-arc
/// cut over the integer unitig graph. `is_trio`/`max_drop_len` are unused in
/// the reference body and dropped here.
#[allow(clippy::too_many_arguments)]
pub fn ulg_arc_cut_supports(
    ctx: &mut UlResolve,
    max_ext: i32,
    max_ext_hifi: u32,
    len_rat: f32,
    topo_level: u32,
    skip_hom: bool,
    collapse_check: bool,
    b: &mut Vec<u64>,
    ub: &mut Vec<u64>,
) -> u32 {
    let n_vtx = (ctx.iug.graph.sequences.len() << 1) as u32;
    b.clear();
    for v in 0..n_vtx {
        if ctx.iug.graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        // `g->seq_vis[v] == 0`: the bubble vertices marked by
        // `asg_arc_identify_simple_bubbles_multi` are skipped.
        if ctx.iug.graph.sequence_visited[v as usize] != 0 {
            continue;
        }
        let nv = ctx.iug.graph.arcs_of(v).len();
        if nv < 2 {
            continue;
        }
        let kv = ctx
            .iug
            .graph
            .arcs_of(v)
            .iter()
            .filter(|a| !a.deleted)
            .take(2)
            .count();
        if kv < 2 {
            continue;
        }
        for (i, arc) in ctx.iug.graph.arcs_of(v).iter().enumerate() {
            if arc.deleted {
                continue;
            }
            b.push(i as u64);
        }
    }
    // `worker_update_ul_arc_supports`: the sequential form of the kt_for
    // pre-pass; every arc word grows its query support in the high half so
    // the radix sort orders the cuts weakest-support first.
    let mut scratch_int: Vec<u64> = Vec::new();
    let mut scratch_raw: Vec<u64> = Vec::new();
    for entry in b.iter_mut() {
        let arc_index = *entry as usize;
        let ve = ctx.iug.graph.arcs[arc_index];
        let mut w_v = 0u64;
        let mut w_r = 0u64;
        get_ul_arc_supports(
            ctx, &ve, &mut scratch_int, &mut scratch_raw, true, &mut w_v, &mut w_r,
        );
        *entry |= w_v << 32;
    }
    radix_sort_u64(b);
    let mut cnt = 0u32;
    for k in 0..b.len() {
        let arc_index = (b[k] & 0xffff_ffff) as usize;
        if ctx.iug.graph.arcs[arc_index].deleted {
            continue;
        }
        let v = (ctx.iug.graph.arcs[arc_index].ul >> 32) as u32;
        let w = ctx.iug.graph.arcs[arc_index].v ^ 1;
        if ctx.iug.graph.sequences[(v >> 1) as usize].deleted
            || ctx.iug.graph.sequences[(w >> 1) as usize].deleted
        {
            continue;
        }
        let nv = ctx.iug.graph.arcs_of(v).len();
        let nw = ctx.iug.graph.arcs_of(w).len();
        if nv <= 1 && nw <= 1 {
            continue;
        }
        let ve = ctx.iug.graph.arcs[arc_index];
        let kv = ctx
            .iug
            .graph
            .arcs_of(v)
            .iter()
            .filter(|a| !a.deleted)
            .count() as u32;
        let we = ctx
            .iug
            .graph
            .arcs_of(w)
            .iter()
            .find(|a| !a.deleted && a.v == v ^ 1)
            .copied();
        let kw = ctx
            .iug
            .graph
            .arcs_of(w)
            .iter()
            .filter(|a| !a.deleted)
            .count() as u32;
        if kv <= 1 && kw <= 1 {
            continue;
        }
        let mut collapse = false;
        if collapse_check {
            if kv > 1 {
                let mut w_q = 0u64;
                let mut w_t = 0u64;
                collapse = check_ul_contain_arc_supports(
                    ctx, &ve, b, ub, skip_hom, &mut w_q, &mut w_t,
                ) != 0;
            }
            if !collapse && kw > 1 {
                if let Some(we) = we {
                    let mut w_q = 0u64;
                    let mut w_t = 0u64;
                    collapse = check_ul_contain_arc_supports(
                        ctx, &we, b, ub, skip_hom, &mut w_q, &mut w_t,
                    ) != 0;
                }
            }
        }
        if !collapse {
            if kv > 1 {
                let mut w_q = 0u64;
                let mut w_t = 0u64;
                get_ul_arc_supports(ctx, &ve, b, ub, skip_hom, &mut w_q, &mut w_t);
                if w_q == u64::MAX {
                    continue;
                }
                if (w_q as f32) > (w_t as f32 * len_rat) {
                    continue;
                }
            }
            if kw > 1 {
                let Some(we) = we else {
                    continue;
                };
                let mut w_q = 0u64;
                let mut w_t = 0u64;
                get_ul_arc_supports(ctx, &we, b, ub, skip_hom, &mut w_q, &mut w_t);
                if w_q == u64::MAX {
                    continue;
                }
                if (w_q as f32) > (w_t as f32 * len_rat) {
                    continue;
                }
            }
        }
        let to_del = check_ulg_to_del(
            ctx,
            v,
            w,
            kv,
            kw,
            max_ext as u32,
            u64::from(max_ext_hifi),
            topo_level,
            collapse,
            b,
            ub,
        );
        if to_del {
            // `ve->del = we->del = 1` (the mirror arcs are left intact,
            // exactly like the reference).
            ctx.iug.graph.arcs[arc_index].deleted = true;
            ctx.iug.graph.delete_first_arc(w, v ^ 1);
            cnt += 1;
        }
    }
    if cnt > 0 {
        ctx.iug.graph.cleanup();
    }
    cnt
}

/// `ulg_seq_del` (`gfa_ut.cpp` 8662): deletes one integer unitig vertex with
/// both arc directions and empties the unitig.
pub fn ulg_seq_del(iug: &mut crate::assembly::unitig::UnitigGraph, s: u32) {
    iug.graph.sequences[s as usize].deleted = true;
    for k in 0..2u32 {
        let v = (s << 1) | k;
        let targets: Vec<u32> = iug.graph.arcs_of(v).iter().map(|a| a.v).collect();
        for t in targets {
            iug.graph.delete_arc(v, t, true);
            iug.graph.delete_arc(t ^ 1, v ^ 1, true);
        }
    }
    iug.unitigs[s as usize] = MaUtg::default();
}

/// `ulg_arc_cut_tips` (`gfa_ut.cpp` 8673).
#[allow(clippy::too_many_arguments)]
pub fn ulg_arc_cut_tips(
    ctx: &mut UlResolve,
    max_ext: u32,
    max_ext_hifi: u32,
    is_double_check: bool,
    b: &mut Vec<u64>,
    ub: &mut Vec<u64>,
) -> u32 {
    let n_vtx = (ctx.iug.graph.sequences.len() << 1) as u32;
    b.clear();
    let mut cnt = 0u32;
    for v in 0..n_vtx {
        if ctx.iug.graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if get_arcs(&ctx.iug.graph, v ^ 1, None, 1) != 0 {
            continue;
        }
        let mut ul_occ = 0u32;
        ctx.get_iug_u_raw_occ(v >> 1, Some(&mut ul_occ), None);
        let mut kv = ul_occ;
        let mut is_telo = ctx
            .telo
            .is_some_and(|t| t[(v >> 1) as usize] != 0);
        let mut w = v;
        for _ in 0..max_ext {
            let (kind, lw) = ctx.iug.graph.is_unitig_end(w ^ 1);
            if kind != crate::assembly::asg::ASG_ET_MERGEABLE {
                break;
            }
            w = lw.expect("a mergeable edge always reports its endpoint") as u32;
            ctx.get_iug_u_raw_occ(w >> 1, Some(&mut ul_occ), None);
            kv += ul_occ;
            if let Some(t) = ctx.telo {
                if t[(w >> 1) as usize] != 0 {
                    is_telo = true;
                }
            }
        }
        if kv <= max_ext && !is_telo {
            b.push((u64::from(kv) << 32) | u64::from(v));
        }
    }
    radix_sort_u64(b);
    for k in 0..b.len() {
        let v = b[k] as u32;
        if ctx.iug.graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if get_arcs(&ctx.iug.graph, v ^ 1, None, 1) != 0 {
            continue;
        }
        let pb = b.len();
        b.push(u64::from(v));
        let mut ul_occ = 0u32;
        ctx.get_iug_u_raw_occ(v >> 1, Some(&mut ul_occ), None);
        let mut kv = ul_occ;
        let mut is_telo = ctx
            .telo
            .is_some_and(|t| t[(v >> 1) as usize] != 0);
        let mut w = v;
        for _ in 0..max_ext {
            let (kind, lw) = ctx.iug.graph.is_unitig_end(w ^ 1);
            if kind != crate::assembly::asg::ASG_ET_MERGEABLE {
                break;
            }
            w = lw.expect("a mergeable edge always reports its endpoint") as u32;
            b.push(u64::from(w));
            ctx.get_iug_u_raw_occ(w >> 1, Some(&mut ul_occ), None);
            kv += ul_occ;
            if let Some(t) = ctx.telo {
                if t[(w >> 1) as usize] != 0 {
                    is_telo = true;
                }
            }
        }
        if !is_telo && kv <= max_ext {
            let mut del_occ = 0u32;
            if get_remove_hifi_occ(
                ctx,
                u64::from(max_ext_hifi),
                &b[pb..],
                ub,
                Some(&mut del_occ),
            ) {
                let (mut w_occ, mut a_occ) = (0u32, 0u32);
                if is_double_check {
                    ctx.get_bridges(&b[pb..], Some(&mut w_occ), Some(&mut a_occ));
                }
                if !is_double_check || w_occ > 0 || a_occ > 0 || kv == 0 || del_occ == 0 {
                    for i in pb..b.len() {
                        ulg_seq_del(ctx.iug, (b[i] >> 1) as u32);
                    }
                    cnt += (b.len() - pb) as u32;
                }
            }
        }
        b.truncate(pb);
    }
    if cnt > 0 {
        ctx.iug.graph.cleanup();
    }
    cnt
}

/// `get_ug_integer_seq_occ` (`gfa_ut.cpp` 10162): the raw read count of the
/// integer sequence of a path.
fn get_ug_integer_seq_occ(ctx: &UlResolve, u_a: &[u64], b: &mut Vec<u64>) -> u32 {
    let bn = b.len();
    gen_ug_integer_seq_on_fly(&*ctx.cov, ctx.z_items, ctx.raw, u_a, b);
    let mut occ = 0u32;
    for k in bn..b.len() {
        occ += ctx.raw.unitigs[(b[k] >> 1) as usize].reads.len() as u32;
    }
    b.truncate(bn);
    occ
}

/// `ul_occ_check` (`gfa_ut.cpp` 10173).
#[allow(clippy::too_many_arguments)]
fn ul_occ_check(
    ctx: &mut UlResolve,
    qocc_ul: u32,
    qocc_hifi: u32,
    tv: u32,
    occ_rate: f32,
    b: &mut Vec<u64>,
    ub: &mut Vec<u64>,
) -> bool {
    let bn = b.len();
    get_ul_path_info(
        &ctx.iug.graph, &*ctx.cov, None, tv, None, None, None, None, None, b,
    );
    let mut tocc_ul = 0u32;
    for z in bn..b.len() {
        let mut ul = 0u32;
        ctx.get_iug_u_raw_occ((b[z] >> 1) as u32, Some(&mut ul), None);
        tocc_ul += ul;
    }
    let path: Vec<u64> = b[bn..].to_vec();
    let tocc_hifi = get_ug_integer_seq_occ(ctx, &path, ub);
    b.truncate(bn);
    (qocc_ul as f32) <= (tocc_ul as f32 * occ_rate)
        && (qocc_hifi as f32) <= (tocc_hifi as f32 * occ_rate)
}

/// `ul_homo_path_check` (`gfa_ut.cpp` 10198).
#[allow(clippy::too_many_arguments)]
fn ul_homo_path_check(
    ctx: &mut UlResolve,
    v: u32,
    w: u32,
    raw_ug_occ: u32,
    match_rate: f32,
    b: &mut Vec<u64>,
    rb: &mut Vec<u64>,
) -> bool {
    let bn = b.len();
    let rbn = rb.len();
    get_ul_path_info(&ctx.iug.graph, &*ctx.cov, None, v, None, None, None, None, None, b);
    let vn = b.len() - bn;
    get_ul_path_info(&ctx.iug.graph, &*ctx.cov, None, w, None, None, None, None, None, b);
    let wn = b.len() - bn - vn;
    let va: Vec<u64> = b[bn..bn + vn].to_vec();
    let wa: Vec<u64> = b[bn + vn..].to_vec();
    gen_ug_integer_seq_on_fly(&*ctx.cov, ctx.z_items, ctx.raw, &va, rb);
    let rvn = rb.len() - rbn;
    let rva: Vec<u64> = rb[rbn..].to_vec();
    gen_ug_integer_seq_on_fly(&*ctx.cov, ctx.z_items, ctx.raw, &wa, rb);
    let rwn = rb.len() - rbn - rvn;
    let rwa: Vec<u64> = rb[rbn + rvn..].to_vec();
    b.truncate(bn);
    for (k, &word) in rva.iter().enumerate() {
        let rid = (word >> 1) as u32;
        let mut x: u64;
        if ctx.bub.is_bubble(rid) {
            x = (u64::from(ctx.bub.index[rid as usize]) | 0x8000_0000) << 32;
        } else {
            x = u64::from(rid) << 32;
        }
        x |= (k as u64) << 1;
        b.push(x);
    }
    for (k, &word) in rwa.iter().enumerate() {
        let rid = (word >> 1) as u32;
        let mut x: u64;
        if ctx.bub.is_bubble(rid) {
            x = (u64::from(ctx.bub.index[rid as usize]) | 0x8000_0000) << 32;
        } else {
            x = u64::from(rid) << 32;
        }
        x |= 1;
        x |= (k as u64) << 1;
        b.push(x);
    }
    radix_sort_u64(&mut b[bn..]);
    let mut o = [0u64; 2];
    let mut l = bn;
    let mut k = bn + 1;
    loop {
        let at_end = k == b.len();
        if at_end || (b[k] >> 32) != (b[l] >> 32) {
            if k - l > 1 {
                let mut c = [0u64; 2];
                let mut z = l;
                while z < k {
                    c[(b[z] & 1) as usize] += 1;
                    if c[0] > 0 && c[1] > 0 {
                        break;
                    }
                    z += 1;
                }
                if c[0] > 0 && c[1] > 0 {
                    for z in l..k {
                        o[(b[z] & 1) as usize] += 1;
                        let mut x = (b[z] << 32) >> 32;
                        x >>= 1;
                        if b[z] & 1 != 0 {
                            x |= 1 << 32;
                        }
                        b[z] = x;
                    }
                } else {
                    for z in l..k {
                        b[z] = u64::MAX;
                    }
                }
            } else {
                for z in l..k {
                    b[z] = u64::MAX;
                }
            }
            l = k;
        }
        if at_end {
            break;
        }
        k += 1;
    }
    radix_sort_u64(&mut b[bn..]);
    // The reference shrinks `b->n` but keeps reading the surviving entries
    // past it through raw pointers; copy them out before the truncation.
    let survivors: Vec<u64> = b[bn..].to_vec();
    b.truncate(bn);
    rb.truncate(rbn);
    if o[0] > 0
        && idx_check_rate(
            &survivors[..o[0] as usize],
            o[0],
            u64::from(raw_ug_occ),
            match_rate,
        )
    {
        return true;
    }
    if o[1] > 0
        && idx_check_rate(
            &survivors[o[0] as usize..o[0] as usize + o[1] as usize],
            o[1],
            u64::from(raw_ug_occ),
            match_rate,
        )
    {
        return true;
    }
    false
}

/// `ulg_arc_cut_z` (`gfa_ut.cpp` 10337). `is_trio`/`max_drop_len` are unused
/// in the reference body and dropped here.
#[allow(clippy::too_many_arguments)]
pub fn ulg_arc_cut_z(
    ctx: &mut UlResolve,
    max_ext: u32,
    max_ext_hifi: u32,
    len_rat: f32,
    small_occ_rate: f32,
    raw_ug_occ: u32,
    raw_match_rate: f32,
    skip_hom: bool,
    b: &mut Vec<u64>,
    ub: &mut Vec<u64>,
) -> u32 {
    let n_vtx = (ctx.iug.graph.sequences.len() << 1) as u32;
    b.clear();
    for v in 0..n_vtx {
        if ctx.iug.graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        let av: Vec<AsgArc> = ctx.iug.graph.arcs_of(v).to_vec();
        let nv = av.len();
        if nv < 2 {
            continue;
        }
        let mut kv = 0u32;
        for arc in &av {
            if arc.deleted {
                continue;
            }
            kv += 1;
            if kv > 2 {
                break;
            }
        }
        if kv != 2 {
            continue;
        }
        for (i, arc) in av.iter().enumerate() {
            if arc.deleted {
                continue;
            }
            let kw = get_arcs(&ctx.iug.graph, arc.v ^ 1, None, 1) as u32;
            if kw == 2 {
                b.push(i as u64);
            } else if kw == 1 {
                ub.clear();
                if get_ul_path_info(
                    &ctx.iug.graph, &*ctx.cov, None, arc.v, None, None, None, None, None, ub,
                ) == 1
                {
                    let mut raw_ul = 0u32;
                    for &word in ub.iter() {
                        let mut ul = 0u32;
                        ctx.get_iug_u_raw_occ((word >> 1) as u32, Some(&mut ul), None);
                        raw_ul += ul;
                    }
                    b.push((u64::from(raw_ul) << 32) | i as u64);
                }
            }
        }
    }
    radix_sort_u64(b);
    let mut cnt = 0u32;
    for k in 0..b.len() {
        let arc_index = (b[k] & 0xffff_ffff) as usize;
        if ctx.iug.graph.arcs[arc_index].deleted {
            continue;
        }
        let v = (ctx.iug.graph.arcs[arc_index].ul >> 32) as u32;
        let w = ctx.iug.graph.arcs[arc_index].v ^ 1;
        if ctx.iug.graph.sequences[(v >> 1) as usize].deleted
            || ctx.iug.graph.sequences[(w >> 1) as usize].deleted
        {
            continue;
        }
        let av: Vec<AsgArc> = ctx.iug.graph.arcs_of(v).to_vec();
        let nv = av.len();
        let aw: Vec<AsgArc> = ctx.iug.graph.arcs_of(w).to_vec();
        let nw = aw.len();
        if nv <= 1 && nw <= 1 {
            continue;
        }
        let mut vp = [v ^ 1, u32::MAX];
        let ve = ctx.iug.graph.arcs[arc_index];
        let mut kv = 0u32;
        for arc in &av {
            if arc.deleted {
                continue;
            }
            if arc.v != w ^ 1 {
                vp[1] = arc.v;
            }
            kv += 1;
        }
        let mut wp = [w ^ 1, u32::MAX];
        let mut we: Option<AsgArc> = None;
        let mut kw = 0u32;
        for arc in &aw {
            if arc.deleted {
                continue;
            }
            if arc.v == v ^ 1 {
                we = Some(*arc);
            } else {
                wp[1] = arc.v;
            }
            kw += 1;
        }
        if kv <= 1 && kw <= 1 {
            continue;
        }
        if kv != 2 || kw > 2 {
            continue;
        }
        let mut raw_ul = 0u32;
        let mut wt = w;
        let mut wte: Option<AsgArc> = we;
        let mut kwt = kw;
        let pb = b.len();
        ub.clear();
        if kw == 1 {
            let kind = get_ul_path_info(
                &ctx.iug.graph,
                &*ctx.cov,
                None,
                w ^ 1,
                Some(&mut wt),
                None,
                None,
                None,
                None,
                b,
            );
            if kind == 1 {
                // TWO_INPUT
                for z in pb..b.len() {
                    let mut ul = 0u32;
                    ctx.get_iug_u_raw_occ((b[z] >> 1) as u32, Some(&mut ul), None);
                    raw_ul += ul;
                }
                let ul_saved = wt;
                // `get_arcs(ug->g, wt, &wt, 1); wt = ug->g->arc[wt].v^1;`
                let mut idx1: Vec<usize> = Vec::new();
                get_arcs(&ctx.iug.graph, wt, Some(&mut idx1), 1);
                wt = ctx.iug.graph.arcs[idx1[0]].v ^ 1;
                wp[0] = wt ^ 1;
                wp[1] = u32::MAX;
                let wtw: Vec<AsgArc> = ctx.iug.graph.arcs_of(wt).to_vec();
                kwt = 0;
                for arc in &wtw {
                    if arc.deleted {
                        continue;
                    }
                    if arc.v == ul_saved ^ 1 {
                        wte = Some(*arc);
                    } else {
                        wp[1] = arc.v;
                    }
                    kwt += 1;
                }
            } else {
                b.truncate(pb);
                continue;
            }
        }
        assert!(kv == 2 && kwt == 2);
        if raw_ul <= max_ext {
            let tail: Vec<u64> = b[pb..].to_vec();
            if get_remove_hifi_occ(ctx, u64::from(max_ext_hifi), &tail, ub, None) {
                let raw_hifi = get_ug_integer_seq_occ(ctx, &tail, ub);
                b.truncate(pb);
                if raw_ul > 0 {
                    let mut ok = true;
                    for tv in [vp[0], vp[1], wp[0], wp[1]] {
                        if !ul_occ_check(ctx, raw_ul, raw_hifi, tv, small_occ_rate, b, ub) {
                            ok = false;
                            break;
                        }
                    }
                    if !ok {
                        continue;
                    }
                }
                let mut w_q = 0u64;
                let mut w_t = 0u64;
                get_ul_arc_supports(ctx, &ve, b, ub, skip_hom, &mut w_q, &mut w_t);
                b.truncate(pb);
                ub.clear();
                if w_q == u64::MAX || (w_q as f32) > (w_t as f32 * len_rat) {
                    continue;
                }
                let wte = wte.expect("the extended target edge must exist");
                let mut w_q = 0u64;
                let mut w_t = 0u64;
                get_ul_arc_supports(ctx, &wte, b, ub, skip_hom, &mut w_q, &mut w_t);
                b.truncate(pb);
                ub.clear();
                if w_q == u64::MAX || (w_q as f32) > (w_t as f32 * len_rat) {
                    continue;
                }
                if !ul_homo_path_check(ctx, vp[0], wp[1], raw_ug_occ, raw_match_rate, b, ub) {
                    continue;
                }
                if !ul_homo_path_check(ctx, wp[0], vp[1], raw_ug_occ, raw_match_rate, b, ub) {
                    continue;
                }
                if kw == 1 {
                    let start = b.len();
                    let mut wt2 = wt;
                    get_ul_path_info(
                        &ctx.iug.graph,
                        &*ctx.cov,
                        None,
                        w ^ 1,
                        Some(&mut wt2),
                        None,
                        None,
                        None,
                        None,
                        b,
                    );
                    for z in start..b.len() {
                        ulg_seq_del(ctx.iug, (b[z] >> 1) as u32);
                    }
                }
                ctx.iug.graph.arcs[arc_index].deleted = true;
                ctx.iug.graph.delete_first_arc(w, v ^ 1);
                cnt += 1;
            }
        }
        b.truncate(pb);
    }
    if cnt > 0 {
        ctx.iug.graph.cleanup();
    }
    cnt
}

/// `ulg_bub_pop_cut_aux` (`gfa_ut.cpp` 9865): the cut form of the bubble pop
/// — deletes the off-path bubble vertices when they carry little coverage.
#[allow(clippy::too_many_arguments)]
fn ulg_bub_pop_cut_aux(
    ctx: &mut UlResolve,
    v0: u32,
    x: &mut crate::assembly::bubble::BubbleBuffer,
    max_ext: u32,
    max_ext_hifi: u32,
    b: &mut Vec<u64>,
    ub: &mut Vec<u64>,
) -> u64 {
    let bn = b.len();
    let mut v = x.sources[0];
    loop {
        let u = x.info[v as usize].parent;
        x.info[v as usize].distance = u32::MAX;
        v = u;
        if v == v0 {
            break;
        }
    }
    let mut n_ext = 0u32;
    for &node in &x.visited {
        if x.info[node as usize].distance == u32::MAX {
            continue;
        }
        let mut ul_occ = 0u32;
        ctx.get_iug_u_raw_occ(node >> 1, Some(&mut ul_occ), None);
        n_ext += ul_occ;
        b.push(u64::from(node));
    }
    let mut r = 0u64;
    if n_ext <= max_ext {
        let tail: Vec<u64> = b[bn..].to_vec();
        if get_remove_hifi_occ(ctx, u64::from(max_ext_hifi), &tail, ub, None) {
            r = 1;
        }
    }
    b.truncate(bn);
    r
}

/// `ulg_bub_pop_backtrack` (`gfa_ut.cpp` 9887): the pop form — deletes the
/// bubble and restores the optimal path.
fn ulg_bub_pop_backtrack(
    iug: &mut crate::assembly::unitig::UnitigGraph,
    v0: u32,
    x: &crate::assembly::bubble::BubbleBuffer,
) -> u64 {
    let mut cnt = x.visited_edges.len() as u64;
    for &node in &x.visited {
        iug.graph.sequences[(node >> 1) as usize].deleted = true;
    }
    for &edge in &x.visited_edges {
        let arc = iug.graph.arcs[edge as usize];
        iug.graph.arcs[edge as usize].deleted = true;
        iug.graph.delete_arc(arc.v ^ 1, ((arc.ul >> 32) as u32) ^ 1, true);
    }
    let mut v = x.sources[0];
    loop {
        let u = x.info[v as usize].parent;
        iug.graph.sequences[(v >> 1) as usize].deleted = false;
        iug.graph.delete_arc(u, v, false);
        iug.graph.delete_arc(v ^ 1, u ^ 1, false);
        cnt -= 1;
        v = u;
        if v == v0 {
            break;
        }
    }
    let visited: Vec<u32> = x.visited.clone();
    for &node in &visited {
        if !iug.graph.sequences[(node >> 1) as usize].deleted {
            continue;
        }
        ulg_seq_del(iug, node >> 1);
    }
    cnt
}

/// `ulg_bub_pop1` (`gfa_ut.cpp` 9941).
#[allow(clippy::too_many_arguments)]
pub fn ulg_bub_pop1(
    ctx: &mut UlResolve,
    v0: u32,
    max_dist: u64,
    x: &mut crate::assembly::bubble::BubbleBuffer,
    max_ext: u32,
    max_ext_hifi: u32,
    check_bubble_only: bool,
    is_pop: bool,
    skip_hom: bool,
    b: &mut Vec<u64>,
    ub: &mut Vec<u64>,
    r_w_c: &mut u32,
    r_w_m: &mut u32,
) -> u64 {
    *r_w_c = u32::MAX;
    *r_w_m = u32::MAX;
    if ctx.iug.graph.sequences[(v0 >> 1) as usize].deleted
        || get_arcs(&ctx.iug.graph, v0, None, 1) < 2
    {
        return 0;
    }
    x.sources.clear();
    x.visited.clear();
    x.visited_edges.clear();
    let head = &mut x.info[v0 as usize];
    head.positive_count = 0;
    head.distance = 0;
    head.negative_count = 0;
    head.node_count = 0;
    head.non_positive_count = 0;
    x.sources.push(v0);
    let mut cnt = 0u64;
    let mut fail_b = false;
    let mut n_tips = 0u32;
    let mut tip_end = u32::MAX;
    let mut n_pending = 0u32;
    loop {
        let v = x.sources.pop().expect("the search stack must not be empty");
        let d = x.info[v as usize].distance;
        let c = x.info[v as usize].positive_count;
        let m = x.info[v as usize].negative_count;
        let av: Vec<AsgArc> = ctx.iug.graph.arcs_of(v).to_vec();
        let nv = av.len();
        let kv = get_arcs(&ctx.iug.graph, v, None, 1) as u32;
        let mut i = 0usize;
        while i < nv {
            let arc = av[i];
            if arc.deleted {
                i += 1;
                continue;
            }
            let w = arc.v;
            // `l = (uint32_t)av[i].ul`: the node length in the low half.
            let l = if v == v0 { 0u32 } else { arc.ul as u32 };
            if (w >> 1) == (v0 >> 1) {
                fail_b = true;
                break;
            }
            let (arc_start, _) = ctx.iug.graph.arc_range(v);
            x.visited_edges.push((arc_start + i) as u32);
            if u64::from(d.wrapping_add(l)) > max_dist {
                fail_b = true;
                break;
            }
            let kw = get_arcs(&ctx.iug.graph, w ^ 1, None, 1) as u32;
            let mut wc = 0u64;
            let wm = ctx.cov.get_ul_occ((w >> 1) as usize);
            if !check_bubble_only {
                if kv > 1 {
                    let mut w_q = 0u64;
                    let mut w_t = 0u64;
                    get_ul_arc_supports(ctx, &arc, b, ub, skip_hom, &mut w_q, &mut w_t);
                    if w_q < w_t {
                        let ww = if w_q == 0 {
                            10u64
                        } else if w_t != 0 {
                            ((w_t - w_q) * 10) / w_t
                        } else {
                            0
                        };
                        if ww > wc {
                            wc = ww;
                        }
                    }
                }
                if kw > 1 {
                    let we = ctx
                        .iug
                        .graph
                        .get_arc(w ^ 1, v ^ 1)
                        .expect("the mirror edge must exist");
                    let mut w_q = 0u64;
                    let mut w_t = 0u64;
                    get_ul_arc_supports(ctx, &we, b, ub, skip_hom, &mut w_q, &mut w_t);
                    if w_q < w_t {
                        let ww = if w_q == 0 {
                            10u64
                        } else if w_t != 0 {
                            ((w_t - w_q) * 10) / w_t
                        } else {
                            0
                        };
                        if ww > wc {
                            wc = ww;
                        }
                    }
                }
            }
            if !x.info[w as usize].state {
                x.visited.push(w);
                let t = &mut x.info[w as usize];
                t.parent = v;
                t.state = true;
                t.distance = d.wrapping_add(l);
                t.positive_count = (u64::from(c) + wc) as u32;
                t.negative_count = (u64::from(m) + wm) as u32;
                t.remaining = kw;
                n_pending += 1;
            } else {
                let mut to_replace = false;
                if u64::from(c) + wc < u64::from(x.info[w as usize].positive_count) {
                    to_replace = true;
                } else if u64::from(c) + wc == u64::from(x.info[w as usize].positive_count)
                    && u64::from(m) + wm > u64::from(x.info[w as usize].negative_count)
                {
                    to_replace = true;
                } else if u64::from(c) + wc == u64::from(x.info[w as usize].positive_count)
                    && u64::from(m) + wm == u64::from(x.info[w as usize].negative_count)
                    && d.wrapping_add(l) > x.info[w as usize].distance
                {
                    to_replace = true;
                }
                if to_replace {
                    let t = &mut x.info[w as usize];
                    t.parent = v;
                    t.positive_count = (u64::from(c) + wc) as u32;
                    t.negative_count = (u64::from(m) + wm) as u32;
                }
                if d.wrapping_add(l) < x.info[w as usize].distance {
                    x.info[w as usize].distance = d.wrapping_add(l);
                }
            }
            x.info[w as usize].remaining -= 1;
            if x.info[w as usize].remaining == 0 {
                let z = get_arcs(&ctx.iug.graph, w, None, 1) as u32;
                if z > 0 {
                    x.sources.push(w);
                } else {
                    if n_tips != 0 {
                        fail_b = true;
                        break;
                    }
                    n_tips += 1;
                    tip_end = w;
                }
                n_pending -= 1;
            }
            i += 1;
        }
        if fail_b {
            break;
        }
        if n_tips == 1 {
            if tip_end != u32::MAX && n_pending == 0 && x.sources.is_empty() {
                x.sources.push(tip_end);
                break;
            }
            fail_b = true;
            break;
        }
        if i < nv || x.sources.is_empty() {
            fail_b = true;
            break;
        }
        if !(x.sources.len() > 1 || n_pending > 0) {
            break;
        }
    }
    if !fail_b {
        cnt = 1;
        *r_w_c = x.info[x.sources[0] as usize].positive_count;
        *r_w_m = x.info[x.sources[0] as usize].negative_count;
        if !check_bubble_only {
            if is_pop {
                cnt = ulg_bub_pop_backtrack(ctx.iug, v0, x);
            } else {
                cnt = ulg_bub_pop_cut_aux(ctx, v0, x, max_ext, max_ext_hifi, b, ub);
            }
        }
    }
    for &node in &x.visited {
        let t = &mut x.info[node as usize];
        t.state = false;
        t.positive_count = 0;
        t.distance = 0;
        t.negative_count = 0;
        t.node_count = 0;
        t.non_positive_count = 0;
    }
    if cnt == 0 {
        *r_w_c = u32::MAX;
        *r_w_m = u32::MAX;
    }
    cnt
}

/// `ulg_pop_bubble` (`gfa_ut.cpp` 10084).
#[allow(clippy::too_many_arguments)]
pub fn ulg_pop_bubble(
    ctx: &mut UlResolve,
    i_max_dist: Option<&mut u64>,
    max_ext: u32,
    max_ext_hifi: u32,
    skip_hom: bool,
    ob: &mut Vec<u64>,
    ub: &mut Vec<u64>,
) -> u64 {
    if !ctx.iug.graph.is_symmetric {
        ctx.iug.graph.symmetrize();
    }
    let n_vtx = (ctx.iug.graph.sequences.len() << 1) as u32;
    let mut x = crate::assembly::bubble::BubbleBuffer::default();
    x.info = vec![crate::assembly::bubble::BubbleInfo::default(); n_vtx as usize];
    let max_dist = match i_max_dist {
        Some(dist) => *dist,
        None => crate::assembly::ul_renew_gfa::get_bub_pop_max_dist_advance(
            &mut ctx.iug.graph,
        ),
    };
    ob.clear();
    ub.clear();
    let mut n_pop = 0u64;
    if max_dist > 0 {
        for v in 0..n_vtx {
            let nv = ctx.iug.graph.arcs_of(v).len();
            if nv < 2 || ctx.iug.graph.sequences[(v >> 1) as usize].deleted {
                continue;
            }
            let n_arc = ctx
                .iug
                .graph
                .arcs_of(v)
                .iter()
                .filter(|a| !a.deleted)
                .count();
            if n_arc < 2 {
                continue;
            }
            ob.clear();
            ub.clear();
            let mut wc = [0u32; 2];
            let mut wm = [0u32; 2];
            if ulg_bub_pop1(
                ctx, v, max_dist, &mut x, max_ext, max_ext_hifi, true, false, skip_hom, ob, ub,
                &mut wc[0], &mut wm[0],
            ) != 0
            {
                let w = x.sources[0] ^ 1;
                let mut mm_c = u32::MAX;
                let mut mm_m = u32::MAX;
                let mut mm_v = u32::MAX;
                ob.clear();
                ub.clear();
                if ulg_bub_pop1(
                    ctx, v, max_dist, &mut x, max_ext, max_ext_hifi, false, false, skip_hom, ob,
                    ub, &mut wc[0], &mut wm[0],
                ) != 0
                {
                    mm_c = wc[0];
                    mm_m = wm[0];
                    mm_v = v;
                }
                ob.clear();
                ub.clear();
                if ulg_bub_pop1(
                    ctx, w, max_dist, &mut x, max_ext, max_ext_hifi, false, false, skip_hom, ob,
                    ub, &mut wc[1], &mut wm[1],
                ) != 0
                {
                    // The reference's condition requires both `wc[1] < mm_c`
                    // and `wc[1] == mm_c`, so it only fires on the unset
                    // sentinel; replicated verbatim.
                    if (wc[1] < mm_c) && (wc[1] == mm_c && wm[1] > mm_m) {
                        mm_c = wc[1];
                        mm_m = wm[1];
                        mm_v = w;
                    }
                }
                if mm_v != u32::MAX {
                    ob.clear();
                    ub.clear();
                    let popped = ulg_bub_pop1(
                        ctx, mm_v, max_dist, &mut x, max_ext, max_ext_hifi, false, true,
                        skip_hom, ob, ub, &mut wc[0], &mut wm[0],
                    );
                    assert!(popped != 0);
                    n_pop += popped;
                }
            }
        }
    }
    if n_pop > 0 {
        ctx.iug.graph.cleanup();
    }
    n_pop
}
/// `ulg_opt_t`, reduced to the fields `u2g_clean` consumes (`gfa_ut.cpp`
/// 16906).
#[derive(Debug, Clone)]
pub struct UlgOptions {
    /// `clean_round`.
    pub clean_round: i64,
    /// `min_path_drop_ratio`.
    pub min_path_drop_ratio: f64,
    /// `max_path_drop_ratio`.
    pub max_path_drop_ratio: f64,
    /// `max_tip`.
    pub max_tip: i64,
    /// `max_tip_hifi`.
    pub max_tip_hifi: i64,
    /// `is_trio`.
    pub is_trio: bool,
}

/// `u2g_clean` (`gfa_ut.cpp` 16281): the whole integer-unitig cleanup —
/// repeated support/tip cuts, bubble pops and the z-shaped chimeric cuts —
/// then the `renew_ul2_utg` rebuild and the `u2g_threading` of `h_usg`.
///
/// `cov`/`bg` are the caller's `mem::take` of `z.cov`/`z.bg`; the rebuild
/// refills both.
#[allow(clippy::too_many_arguments)]
pub fn u2g_clean(
    cov: &mut UlCov,
    bg: &mut UlBg,
    z_items: &Ul2ulIdx,
    iug: &mut crate::assembly::unitig::UnitigGraph,
    i_g: &Asg,
    raw: &crate::assembly::unitig::UnitigGraph,
    sg: &Asg,
    store: &[crate::assembly::ul_store::UlReadVec],
    raw_unitigs: &[MaUtg],
    read_lengths: &[u64],
    unitig_reads: &[Vec<u64>],
    str_idx: &crate::assembly::ul_integer::UlStrIdx,
    bub: &BubbleType,
    sources: &[Vec<crate::overlap::paf::PafHit>],
    reads: &crate::reads::ReadStore,
    maximum_hang: i32,
    maximum_hang_rate: f32,
    minimum_overlap: i64,
    b_mask: &mut crate::assembly::bubble::BubbleLabel,
    options: &UlgOptions,
    keep_raw_utg: bool,
    is_bridg: bool,
) -> Usg {
    // `CALLOC(iug->g->seq_vis, iug->g->n_seq*2)`
    iug.graph.sequence_visited = vec![0u8; iug.graph.sequences.len() * 2];
    let mut bu: Vec<u64> = Vec::new();
    let mut uu: Vec<u64> = Vec::new();
    let step = if options.clean_round == 1 {
        options.max_path_drop_ratio
    } else {
        (options.max_path_drop_ratio - options.min_path_drop_ratio)
            / ((options.clean_round - 1) as f64)
    };
    let mut max_tip_hifi = options.max_tip_hifi;
    for _ss in 0..2i64 {
        for i in 0..options.clean_round {
            let mut drop = options.min_path_drop_ratio + step * (i as f64);
            if drop > options.max_path_drop_ratio {
                drop = options.max_path_drop_ratio;
            }
            for (topo_level, mm_tip, double_check) in [
                (2u32, options.max_tip, true),
                (1, options.max_tip, true),
                (1, 0x7fff_ffffi64, false),
            ] {
                let mut cnt = 1u32;
                while cnt != 0 {
                    cnt = 0;
                    crate::assembly::bubble_identify::asg_arc_identify_simple_bubbles_multi(
                        &mut iug.graph, b_mask, false, 0,
                    );
                    let mut ctx = UlResolve {
                        cov: &mut *cov,
                        z_items,
                        iug,
                        raw,
                        store,
                        raw_unitigs,
                        read_lengths,
                        unitig_reads,
                        str_idx,
                        bub,
                        bg: &*bg,
                        // `uidx->uovl.telo`: only populated with
                        // `--telo-motif`, which is not ported yet.
                        telo: None,
                    };
                    cnt += ulg_arc_cut_supports(
                        &mut ctx,
                        mm_tip as i32,
                        max_tip_hifi as u32,
                        drop as f32,
                        topo_level,
                        true,
                        keep_raw_utg,
                        &mut bu,
                        &mut uu,
                    );
                    cnt += ulg_arc_cut_tips(
                        &mut ctx,
                        mm_tip as u32,
                        max_tip_hifi as u32,
                        double_check,
                        &mut bu,
                        &mut uu,
                    );
                }
            }
        }
        {
            let mut ctx = UlResolve {
                cov: &mut *cov,
                z_items,
                iug,
                raw,
                store,
                raw_unitigs,
                read_lengths,
                unitig_reads,
                str_idx,
                bub,
                bg: &*bg,
                telo: None,
            };
            ulg_pop_bubble(
                &mut ctx,
                None,
                0x7fff_ffffu32,
                max_tip_hifi as u32,
                true,
                &mut bu,
                &mut uu,
            );
            while ulg_arc_cut_z(
                &mut ctx,
                0x7fff_ffffu32,
                max_tip_hifi as u32,
                1.5,
                0.15,
                100,
                0.8,
                true,
                &mut bu,
                &mut uu,
            ) != 0
            {}
        }
        max_tip_hifi <<= 1;
    }
    // `u2g_threading(uidx, ulopt, 3, is_bridg, &bu, &uu)` — with the
    // `renew_ul2_utg` rebuild of `i_ug` and the coverage/bridge tables first.
    crate::assembly::unitig::renew_utg(iug, i_g, &[]);
    iug.graph.sequence_visited = vec![0u8; iug.graph.sequences.len() * 2];
    crate::assembly::ul2ul::renew_u2g_cov(
        z_items,
        &iug.unitigs,
        raw,
        sources,
        reads,
        maximum_hang,
        maximum_hang_rate,
        minimum_overlap as i32,
        cov,
        str_idx,
    );
    *bg = crate::assembly::ul2ul::renew_u2g_bg(cov, str_idx, &iug.graph, raw, bub);
    crate::assembly::ul2ul::u2g_threading(
        iug,
        raw,
        sg,
        cov,
        str_idx,
        store,
        bub,
        minimum_overlap,
        3,
        is_bridg,
        &mut bu,
        &mut uu,
    )
}

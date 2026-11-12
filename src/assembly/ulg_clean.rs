//! The `usg` graph-cleaning cluster of `u2g_hybrid_clean` (`gfa_ut.cpp`
//! 10519-11150 + 15332 + 15616-16050): tip cuts, bubble identification and
//! the bridge marking over the hybrid UL/unitig graph (`usg_t`).

use crate::assembly::ul2ul::{Usg, UsgArc};

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

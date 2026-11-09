//! The tangle search behind `is_small_tangle` (`Overlaps.cpp` 22547).
//!
//! `gen_double_scaffold_gfa` links two contigs across the gap between them
//! when the reference unitig graph has a short tangle path joining their
//! flanks. This module ports that search: flank extraction
//! (`get_ref_nid_flank`), the two BFS variants and the trace-back.

use crate::assembly::asg::Asg;
use crate::assembly::ul::ScafRes;
use crate::assembly::hap_cov::UTrans;
use crate::assembly::unitig::UnitigGraph;

/// `scaf_hit2arc` (`Overlaps.cpp` 21876).
pub fn scaf_hit2arc(
    a0: &UTrans,
    a1: &UTrans,
    rev: bool,
    el: bool,
    glen: i64,
    ul_out: &mut u64,
    v_out: &mut u32,
    ol_out: &mut u32,
    el_out: &mut bool,
    strong_out: &mut bool,
) {
    let mut v = (a0.target_identifier << 1) | u32::from(a0.reverse);
    v ^= u32::from(rev);
    let mut w = (a1.target_identifier << 1) | u32::from(a1.reverse);
    w ^= u32::from(rev);
    *ul_out = u64::from(v) << 32;
    *v_out = w;
    *ol_out = glen.unsigned_abs() as u32;
    *el_out = el;
    *strong_out = glen >= 0;
}

/// `is_noisy_ov` (`Overlaps.cpp` 21888).
pub fn is_noisy_ov(
    mut s: u64,
    mut e: u64,
    len: u64,
    bd: u64,
    cov_rate: f64,
    a: &[UTrans],
) -> bool {
    if e - s < (bd << 1) {
        let w = (e + s) >> 1;
        s = w.saturating_sub(bd);
        e = (w + bd).min(len);
    }
    if s >= e {
        return false;
    }
    let mut tot = 0u64;
    for item in a {
        if item.query_start as u64 >= e {
            break;
        }
        if item.query_end as u64 <= s {
            continue;
        }
        let zs = s.max(item.query_start as u64);
        let ze = e.min(item.query_end as u64);
        if ze > zs {
            tot += ze - zs;
        }
    }
    tot > 0 && tot as f64 >= (e - s) as f64 * cov_rate
}

/// `set_fsk` (`Overlaps.cpp` 22648): marks the unitigs of one contig path.
pub fn set_fsk(fsk: &mut [u8], value: u8, scaf: &crate::assembly::ul::UlVec) {
    for k in 0..scaf.n {
        let header = scaf.bb[k];
        for block in &scaf.bb[header.ts as usize..header.te as usize] {
            fsk[block.hid as usize] = value;
        }
    }
}

/// `get_ref_nid_flank` (`Overlaps.cpp` 22040): the reference path flanking a
/// transition end, plus the distance of every node to that end.
#[allow(clippy::too_many_arguments)]
pub fn get_ref_nid_flank(
    z: &UTrans,
    hu1: &UnitigGraph,
    hu2: &UnitigGraph,
    cp1: &ScafRes,
    cp2: &ScafRes,
    is_beg: bool,
    gap_bd: u64,
    flank_len: u64,
    res: &mut Vec<u32>,
) -> usize {
    let (sc, id) = if (z.target_identifier as usize) < hu1.unitigs.len() {
        (cp1, z.target_identifier as usize)
    } else {
        (cp2, z.target_identifier as usize - hu1.unitigs.len())
    };
    let _ = hu2;
    let (mut s_end, mut e_end) = if is_beg { (false, true) } else { (true, false) };
    if z.reverse {
        s_end ^= true;
        e_end ^= true;
    }

    let idx = &sc.a[id];
    let mut best: Option<(usize, usize)> = None; // (first block, count)
    let mut max_len = 0u64;
    for k in 0..idx.n {
        let header = idx.bb[k];
        let a = &idx.bb[header.ts as usize..header.te as usize];
        let a_n = a.len();
        let os = u64::from(z.target_start).max(u64::from(header.qs));
        let oe = u64::from(z.target_end).min(u64::from(header.qe));
        if os >= oe {
            continue;
        }
        let mut m = 0usize;
        while m < a_n && u64::from(a[m].qe) <= os {
            m += 1;
        }
        if m >= a_n {
            continue;
        }
        let m0 = m;
        while m < a_n && u64::from(a[m].qs) < oe {
            m += 1;
        }
        if m <= m0 {
            continue;
        }
        let os = u64::from(z.target_start).max(u64::from(a[m0].qs));
        let oe = u64::from(z.target_end).min(u64::from(a[m - 1].qe));
        if oe <= os {
            continue;
        }
        let clen = oe - os;
        let cgap = if s_end {
            os - u64::from(z.target_start)
        } else if e_end {
            u64::from(z.target_end) - oe
        } else {
            0
        };
        if cgap > gap_bd {
            continue;
        }
        if clen > max_len {
            max_len = clen;
            best = Some((header.ts as usize + m0, m - m0));
        }
    }
    let Some((start, count)) = best else {
        return 0;
    };
    let am: Vec<crate::assembly::ul::UcBlock> = idx.bb[start..start + count].to_vec();
    let amn = am.len();
    if amn == 0 {
        return 0;
    }

    let (mut as_i, mut ae_i) = (0usize, amn);
    if s_end {
        let mut ami = 0usize;
        while ami < amn && u64::from(z.target_start) + flank_len >= u64::from(am[ami].qe) {
            ami += 1;
        }
        if ami < amn && u64::from(z.target_start) + flank_len >= u64::from(am[ami].qs) {
            ami += 1;
        }
        as_i = 0;
        ae_i = ami;
    } else if e_end {
        let mut ami = amn as i64 - 1;
        while ami >= 0 && u64::from(z.target_end) <= u64::from(am[ami as usize].qs) + flank_len {
            ami -= 1;
        }
        if ami >= 0 && u64::from(z.target_end) <= u64::from(am[ami as usize].qe) + flank_len {
            ami -= 1;
        }
        ami += 1;
        as_i = ami as usize;
        ae_i = amn;
    }
    if as_i >= ae_i {
        return 0;
    }

    // Drop duplicated unitigs.
    if s_end {
        let mut ami = as_i;
        while ami < ae_i {
            let mut amz = ami + 1;
            while amz < ae_i && am[ami].hid != am[amz].hid {
                amz += 1;
            }
            if amz < ae_i {
                ae_i = ami + 1;
                break;
            }
            ami += 1;
        }
    } else if e_end {
        let mut ami = ae_i as i64 - 1;
        while ami >= as_i as i64 {
            let mut amz = ami - 1;
            while amz >= as_i as i64 && am[ami as usize].hid != am[amz as usize].hid {
                amz -= 1;
            }
            if amz >= as_i as i64 {
                as_i = ami as usize;
                break;
            }
            ami -= 1;
        }
    }
    if as_i >= ae_i {
        return 0;
    }

    let rn0 = res.len();
    for ami in as_i..ae_i {
        let mut m = (am[ami].hid << 1) | u32::from(am[ami].rev);
        m ^= u32::from(z.reverse);
        res.push(m);
    }
    let rn = res.len() - rn0;
    if z.reverse {
        res[rn0..rn0 + rn].reverse();
    }
    if s_end {
        for ami in as_i..ae_i {
            let m = am[ami].qs.saturating_sub(z.target_start);
            res.push(m);
        }
    } else if e_end {
        for ami in as_i..ae_i {
            let m = z.target_end.saturating_sub(am[ami].qe);
            res.push(m);
        }
    }
    if z.reverse {
        // The reference re-derives the pointer as `res->n - rn`, so the second
        // reversal covers the last `rn` words: the distance block, not the
        // node block reversed just above.
        let end = res.len();
        res[end - rn..end].reverse();
    }
    rn
}

/// `dedup_tang_flank` (`Overlaps.cpp` 22314): trims the flanks at the first
/// shared unitig. `buf` is `[va(vn), vd(vn), wa(wn), wd(wn)]`.
fn dedup_tang_flank(buf: &mut [u32], vn0: usize, wn0: usize) -> (usize, usize) {
    let mut vn = vn0;
    let mut wn = wn0;
    // The reference re-evaluates `k < vn && k < wn` on every step, so the
    // bounds move as `vn`/`wn` shrink.
    let mut k = 0usize;
    while k < vn && k < wn {
        if k < vn {
            let vp = buf[vn - k - 1] >> 1;
            let mut z = k;
            while z < wn && vp != (buf[2 * vn0 + z] >> 1) {
                z += 1;
            }
            if z < wn {
                wn = z;
            }
        }
        if k < wn {
            let wp = buf[2 * vn0 + k] >> 1;
            let mut z = k;
            while z < vn && wp != (buf[vn - z - 1] >> 1) {
                z += 1;
            }
            if z < vn {
                vn = z;
            }
        }
        k += 1;
    }
    if vn < vn0 {
        let shift = vn0 - vn;
        for k in 0..vn {
            buf[k] = buf[k + shift];
            buf[vn0 + k] = buf[vn0 + k + shift];
        }
    }
    (vn, wn)
}

/// `cal_self_close_tang_flank0` (`Overlaps.cpp` 22348).
#[allow(clippy::too_many_arguments)]
fn cal_self_close_tang_flank0(
    g: &Asg,
    v0: u32,
    s0: u32,
    e0: u32,
    sf: u8,
    ef: u8,
    b1: &mut Vec<u32>,
    b2: &mut Vec<u32>,
    f: &mut [u8],
    fs: &mut [u8],
    max_path: u64,
) -> u64 {
    if f[v0 as usize] != 0 {
        return u64::MAX;
    }
    let mut se = false;
    let mut tot = 0u64;
    b1.clear();
    b2.clear();
    b1.push(v0);
    while let Some(v) = b1.pop() {
        if f[v as usize] != 0 && f[(v ^ 1) as usize] != 0 {
            continue;
        }
        let mut v = v;
        if f[v as usize] == 0 && v != (s0 ^ 1) && v != e0 {
            if fs[(v ^ 1) as usize] != sf && fs[v as usize] != ef {
                f[v as usize] = 1;
                b2.push(v);
                for arc in g.arcs_of(v) {
                    if arc.deleted {
                        continue;
                    }
                    if arc.v == s0 || arc.v == (e0 ^ 1) {
                        se = true;
                        break;
                    }
                    if fs[arc.v as usize] == sf {
                        fs[arc.v as usize] = 0;
                        if f[(arc.v ^ 1) as usize] == 0 {
                            b1.push(arc.v ^ 1);
                        }
                    }
                    if fs[(arc.v ^ 1) as usize] == ef {
                        fs[(arc.v ^ 1) as usize] = 0;
                        if f[arc.v as usize] == 0 {
                            b1.push(arc.v);
                        }
                    }
                    if f[arc.v as usize] == 0 {
                        b1.push(arc.v);
                    }
                }
                if se {
                    break;
                }
            }
        }

        v ^= 1;
        if f[v as usize] == 0 && v != (s0 ^ 1) && v != e0 {
            if fs[(v ^ 1) as usize] != sf && fs[v as usize] != ef {
                f[v as usize] = 1;
                b2.push(v);
                for arc in g.arcs_of(v) {
                    if arc.deleted {
                        continue;
                    }
                    if arc.v == s0 || arc.v == (e0 ^ 1) {
                        se = true;
                        break;
                    }
                    if fs[arc.v as usize] == sf {
                        fs[arc.v as usize] = 0;
                        if f[(arc.v ^ 1) as usize] == 0 {
                            b1.push(arc.v ^ 1);
                        }
                    }
                    if fs[(arc.v ^ 1) as usize] == ef {
                        // The reverse half only clears the flag; the node is
                        // reached later through its own forward arcs.
                        fs[(arc.v ^ 1) as usize] = 0;
                    }
                    if f[arc.v as usize] == 0 {
                        b1.push(arc.v);
                    }
                }
                if se {
                    break;
                }
            }
        }

        if f[v as usize] != 0 && f[(v ^ 1) as usize] != 0 {
            let mut o = 0u32;
            for arc in g.arcs_of(v) {
                if arc.deleted {
                    continue;
                }
                if f[(arc.v ^ 1) as usize] != 0 && f[arc.v as usize] != 0 && o < arc.overlap {
                    o = arc.overlap;
                }
            }
            let o0 = o;
            let mut o = 0u32;
            for arc in g.arcs_of(v ^ 1) {
                if arc.deleted {
                    continue;
                }
                if f[(arc.v ^ 1) as usize] != 0 && f[arc.v as usize] != 0 && o < arc.overlap {
                    o = arc.overlap;
                }
            }
            let o1 = o;
            let length = g.sequences[(v >> 1) as usize].length;
            if o0 + o1 < length {
                tot += u64::from(length - o0 - o1);
            }
        }
        if tot > max_path {
            break;
        }
    }
    if !se && tot <= max_path {
        tot
    } else {
        u64::MAX
    }
}

/// `cal_self_close_tang_flank1` (`Overlaps.cpp` 22441).
fn cal_self_close_tang_flank1(
    g: &Asg,
    s0: u32,
    e0: u32,
    b1: &mut Vec<u32>,
    b2: &mut Vec<u32>,
    f: &mut [u8],
    fsk: &[u8],
    max_path: u64,
) -> u64 {
    let mut se = false;
    let mut tot = 0u64;
    b1.clear();
    b2.clear();
    b1.push(s0);
    while let Some(v) = b1.pop() {
        if f[v as usize] != 0 && f[(v ^ 1) as usize] != 0 {
            continue;
        }
        let mut v = v;
        if f[v as usize] == 0 && v != (s0 ^ 1) && v != e0 {
            f[v as usize] = 1;
            b2.push(v);
            for arc in g.arcs_of(v) {
                if arc.deleted {
                    continue;
                }
                if arc.v == s0 || arc.v == (e0 ^ 1) {
                    se = true;
                    break;
                }
                if f[arc.v as usize] == 0 {
                    b1.push(arc.v);
                }
            }
            if se {
                break;
            }
        }
        v ^= 1;
        if f[v as usize] == 0 && v != (s0 ^ 1) && v != e0 {
            f[v as usize] = 1;
            b2.push(v);
            for arc in g.arcs_of(v) {
                if arc.deleted {
                    continue;
                }
                if arc.v == s0 || arc.v == (e0 ^ 1) {
                    se = true;
                    break;
                }
                if f[arc.v as usize] == 0 {
                    b1.push(arc.v);
                }
            }
            if se {
                break;
            }
        }
        if f[v as usize] != 0 && f[(v ^ 1) as usize] != 0 && fsk[(v >> 1) as usize] == 0 {
            let mut o = 0u32;
            for arc in g.arcs_of(v) {
                if arc.deleted {
                    continue;
                }
                if f[(arc.v ^ 1) as usize] != 0
                    && f[arc.v as usize] != 0
                    && fsk[(arc.v >> 1) as usize] == 0
                    && o < arc.overlap
                {
                    o = arc.overlap;
                }
            }
            let o0 = o;
            let mut o = 0u32;
            for arc in g.arcs_of(v ^ 1) {
                if arc.deleted {
                    continue;
                }
                if f[(arc.v ^ 1) as usize] != 0
                    && f[arc.v as usize] != 0
                    && fsk[(arc.v >> 1) as usize] == 0
                    && o < arc.overlap
                {
                    o = arc.overlap;
                }
            }
            let o1 = o;
            let length = g.sequences[(v >> 1) as usize].length;
            if o0 + o1 < length {
                tot += u64::from(length - o0 - o1);
            }
        }
        if tot > max_path {
            break;
        }
    }
    if !se && tot <= max_path {
        tot
    } else {
        u64::MAX
    }
}

/// `trace_back_tangle_dis` (`Overlaps.cpp` 22497).
#[allow(clippy::too_many_arguments)]
fn trace_back_tangle_dis(
    g: &Asg,
    s: &[u32],
    sd: &[u32],
    e: &[u32],
    ed: &[u32],
    b: &mut Vec<u32>,
    b1: &mut Vec<u32>,
    b2: &mut Vec<u32>,
    f: &mut [u8],
    fs: &[u8],
    fsk: &[u8],
    sf: u8,
    ef: u8,
    max_path: u64,
    is_dip: &mut u64,
) -> u64 {
    b.clear();
    *is_dip = 0;
    let mut sm = 0usize;
    for sk in 0..s.len() {
        let idx = s.len() - sk - 1;
        let v = s[idx];
        if fs[v as usize] == sf && f[v as usize] != 0 && f[(v ^ 1) as usize] == 0 {
            b.push(idx as u32);
            sm += 1;
        }
    }
    if sm == 0 {
        return u64::MAX;
    }
    let mut em = 0usize;
    for ek in 0..e.len() {
        let v = e[ek];
        if fs[v as usize] == ef && f[v as usize] == 0 && f[(v ^ 1) as usize] != 0 {
            b.push(ek as u32);
            em += 1;
        }
    }
    if em == 0 {
        return u64::MAX;
    }
    let si: Vec<u32> = b[..sm].to_vec();
    let ei: Vec<u32> = b[sm..sm + em].to_vec();
    for &sidx in &si {
        let s0 = s[sidx as usize];
        let slen = u64::from(sd[sidx as usize]);
        for &eidx in &ei {
            let e0 = e[eidx as usize];
            let elen = u64::from(ed[eidx as usize]);
            for &z in b2.iter() {
                f[z as usize] = 0;
            }
            let tot = cal_self_close_tang_flank1(
                g,
                s0,
                e0,
                b1,
                b2,
                f,
                fsk,
                max_path + slen + elen,
            );
            if tot != u64::MAX && tot <= max_path + slen + elen {
                if fsk[(s0 >> 1) as usize] != 0 && fsk[(e0 >> 1) as usize] != 0 {
                    *is_dip = 1;
                    return 0;
                }
                return if tot >= slen + elen {
                    tot - slen - elen
                } else {
                    0
                };
            }
        }
    }
    u64::MAX
}

/// `cal_self_close_tang_flank` (`Overlaps.cpp` 22540).
#[allow(clippy::too_many_arguments)]
fn cal_self_close_tang_flank(
    g: &Asg,
    s: &[u32],
    sd: &[u32],
    e: &[u32],
    ed: &[u32],
    b: &mut Vec<u32>,
    b1: &mut Vec<u32>,
    b2: &mut Vec<u32>,
    f: &mut [u8],
    fs: &mut [u8],
    fsk: &[u8],
    max_path: u64,
    is_dip: &mut u64,
) -> u64 {
    *is_dip = 0;
    if s.is_empty() || e.is_empty() {
        return u64::MAX;
    }
    let sf = 1u8;
    let ef = 2u8;
    for &v in s {
        fs[v as usize] = sf;
    }
    let s0 = s[0];
    for &v in e {
        fs[v as usize] = ef;
    }
    let e0 = e[e.len() - 1];

    b1.clear();
    b2.clear();
    let mut tot = u64::MAX;
    for k in 0..s.len() {
        let v = s[s.len() - k - 1];
        if fs[v as usize] != sf {
            continue;
        }
        tot = cal_self_close_tang_flank0(
            g,
            v,
            s0,
            e0,
            sf,
            ef,
            b1,
            b2,
            f,
            fs,
            max_path,
        );
        if tot == u64::MAX || tot > max_path {
            tot = u64::MAX;
            break;
        }
        let found = trace_back_tangle_dis(
            g, s, sd, e, ed, b, b1, b2, f, fs, fsk, sf, ef, max_path, is_dip,
        );
        if found != u64::MAX {
            tot = found;
            break;
        }
        for &v in s {
            fs[v as usize] = sf;
        }
        for &v in e {
            fs[v as usize] = ef;
        }
        for &z in b2.iter() {
            f[z as usize] = 0;
        }
    }

    for &v in s {
        fs[v as usize] = 0;
    }
    for &v in e {
        fs[v as usize] = 0;
    }
    for &z in b2.iter() {
        f[z as usize] = 0;
    }
    tot
}

/// `is_small_tangle` (`Overlaps.cpp` 22547): the gap between two contig
/// transitions, or `u64::MAX` when no short tangle joins them.
#[allow(clippy::too_many_arguments)]
pub fn is_small_tangle(
    s: &UTrans,
    e: &UTrans,
    reference: &UnitigGraph,
    hu1: &UnitigGraph,
    hu2: &UnitigGraph,
    cp1: &ScafRes,
    cp2: &ScafRes,
    f: &mut [u8],
    fs: &mut [u8],
    fsk: &[u8],
    max_len: u64,
    is_dip: &mut u64,
) -> u64 {
    *is_dip = 0;
    let mut flank = Vec::new();
    let vn0 = get_ref_nid_flank(s, hu1, hu2, cp1, cp2, true, 1000, max_len, &mut flank);
    if vn0 == 0 {
        return u64::MAX;
    }
    let wn0 = get_ref_nid_flank(e, hu1, hu2, cp1, cp2, false, 1000, max_len, &mut flank);
    if wn0 == 0 {
        return u64::MAX;
    }
    let (vn, wn) = dedup_tang_flank(&mut flank, vn0, wn0);
    if vn == 0 || wn == 0 {
        return u64::MAX;
    }
    // The reference computes `va`, `vd`, `wa` and `wd` before the dedup, so
    // the distance and second-node blocks keep their original offsets.
    let va = flank[..vn].to_vec();
    let vd = flank[vn0..vn0 + vn].to_vec();
    let wa = flank[2 * vn0..2 * vn0 + wn].to_vec();
    let wd = flank[2 * vn0 + wn0..2 * vn0 + wn0 + wn].to_vec();
    let mut b1 = Vec::new();
    let mut b2 = Vec::new();
    let mut b3 = Vec::new();
    let mut tot = cal_self_close_tang_flank(
        &reference.graph,
        &va,
        &vd,
        &wa,
        &wd,
        &mut b3,
        &mut b1,
        &mut b2,
        f,
        fs,
        fsk,
        max_len,
        is_dip,
    );
    if tot > max_len {
        tot = u64::MAX;
    }
    tot
}

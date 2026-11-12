//! The `gl_chain` cluster (`inter.cpp` 2412-7290): turns one UL read's
//! overlap list into graph chains over the unitig graph.
//!
//! Two drivers consume it: `gl_chain` (the rescall worker, with
//! `gl_chain_graph_adv`) and `gl_chain_refine_advance` (the recorrect worker,
//! with `gl_chain_advance`). Dead variants skipped: the `trav_rate`
//! `gl_chain_advance` (commented out in the reference), `gl_chain_refine`,
//! `gl_chain_advance_back` + `determine_containment_chain` (only reachable
//! from the dead `_with_trans` combine), `infer_read_ovlp`, `get_het_site`,
//! `gl_exact_chain`, `dump_all_chain`, `ins_merge_ul_ov`,
//! `debug_*`/`prt_*` printers.

use crate::assembly::asg::{Asg, MaUtg};
use crate::assembly::hit::ma_hit_to_arc;
use crate::assembly::ul::{
    get_offset_adjust, hc_shortest_k, infer_rovlp_unitigs, max_ovlp, find_mg_lchain_max,
    MgLchain, MgPathDst, MgPathv, UlOv,
};
use crate::assembly::ul_phase::SEC_MODE;
use crate::assembly::ul_rescall::UlRescallContext;
use crate::container::ksort::{
    introsort_klib_by, radix_sort_by_key_32, radix_sort_i64, radix_sort_u64,
};
use crate::correction::haplotype_evidence::{HaplotypeEvidence, HaplotypeEvidenceAlloc};
use crate::correction::region::THRESHOLD_MAXIMUM;
use crate::overlap::region::{ChainData, OverlapRegion};
use crate::reads::read_store::ReadStore;

/// `UG_TRANS_W` (`inter.h` 22).
pub const UG_TRANS_W: i64 = 2;
/// `UG_TRANS_ERR_W` (`inter.h` 24).
pub const UG_TRANS_ERR_W: i64 = 64;
/// `W_CHN_PEN_GAP` (`inter.h` 28).
pub const W_CHN_PEN_GAP: f32 = 0.1;
/// `ERROR_RATE` (`Correct.h` 22).
pub const ERROR_RATE: f64 = 1.5;
/// `UG_SKIP` (`inter.h` 15).
pub const UG_SKIP: i64 = 5;
/// `G_CHAIN_TRANS_RATE` (`inter.h` 25).
pub const G_CHAIN_TRANS_RATE: f32 = 0.25;
/// `G_CHAIN_TRANS_WEIGHT` (`inter.h` 26).
pub const G_CHAIN_TRANS_WEIGHT: i64 = -1;
/// `P_CHAIN_COV` (`inter.h` 9).
pub const P_CHAIN_COV: f32 = 0.985;
/// `P_FRAGEMENT_CHAIN_COV` (`inter.h` 10).
pub const P_FRAGEMENT_CHAIN_COV: f32 = 0.20;
/// `P_FRAGEMENT_PRIMARY_CHAIN_COV` (`inter.h` 11).
pub const P_FRAGEMENT_PRIMARY_CHAIN_COV: f32 = 0.70;
/// `PRIMARY_UL_CHAIN_MIN` (`inter.h` 30).
pub const PRIMARY_UL_CHAIN_MIN: u64 = 75000;
/// `OU_MASK` (`Overlaps.h` 189).
pub const OU_MASK: u32 = 0x3fff;
/// `MG_MAX_SHORT_K` (`inter.cpp` 50).
pub const MG_MAX_SHORT_K: usize = 15;

/// `aln_sc` (`inter.cpp` 5221).
#[inline]
pub fn aln_sc(region: &OverlapRegion, trans_w: i64, err_w: i64) -> i64 {
    i64::from(region.alignment_length)
        - i64::from(region.overlap_length as i64 - i64::from(region.alignment_length) * trans_w)
        - i64::from(region.non_homopolymer_errors) * err_w
}

/// `push_sc_pre` (`inter.cpp` 3697).
#[inline]
pub fn push_sc_pre(mm_sc: i64, mm_idx: i64) -> u64 {
    let sc = if mm_sc >= 0 {
        (mm_sc as u32) | 0x8000_0000
    } else {
        (-(mm_sc)) as u32
    };
    let mut x = u64::from(sc);
    x <<= 32;
    x |= if mm_idx >= 0 { mm_idx as u64 } else { 0x7fff_ffff };
    x
}

/// `pop_sc` (`inter.cpp` 3704).
#[inline]
pub fn pop_sc(x: u64) -> i64 {
    let x = x >> 32;
    if x & 0x8000_0000 != 0 {
        (x - u64::from(0x8000_0000u32)) as i64
    } else {
        (x as i64) * -1
    }
}

/// `pop_pre` (`inter.cpp` 3714).
#[inline]
pub fn pop_pre(x: u64) -> i64 {
    if x & 0x7fff_ffff == 0x7fff_ffff {
        -1
    } else {
        (x & 0x7fff_ffff) as i64
    }
}

/// `find_ul_ov_max` (`inter.cpp` 2504): the last chain whose `qe` is `< x`.
pub fn find_ul_ov_max(a: &[UlOv], x: u32) -> i64 {
    let n = a.len() as i64;
    if n == 0 {
        return -1;
    }
    if a[(n - 1) as usize].qe < x {
        return n - 1;
    }
    if a[0].qe >= x {
        return -1;
    }
    let (mut s, mut e) = (0i64, n);
    while e > s {
        let m = s + (e - s) / 2;
        if a[m as usize].qe >= x {
            e = m;
        } else {
            s = m + 1;
        }
    }
    debug_assert_eq!(s, e);
    s
}

/// `max_ovlp_src` (`inter.cpp` 1505).
pub fn max_ovlp_src(context: &UlRescallContext<'_>, v: u32) -> i64 {
    let mut o = 0u32;
    let x = (v >> 1) as usize;
    for hit in &context.sources[x] {
        let qn = hit.query_identifier;
        let tn = hit.target_identifier;
        let (result, arc) = ma_hit_to_arc(
            hit,
            context.reads.read_length(qn as usize) as i32,
            context.reads.read_length(tn as usize) as i32,
            context.maximum_hang as i32,
            context.maximum_hang_rate,
            context.minimum_overlap as i32,
        );
        if result < 0 {
            continue;
        }
        if (arc.ul >> 32) as u32 != v {
            continue;
        }
        if o < arc.overlap {
            o = arc.overlap;
        }
    }
    i64::from(o)
}

/// `get_ecov_adv` (`inter.cpp` 3405) with the full `mode`/`contain_off`
/// surface; `mode = 0` reads the unitig graph, `mode = 1` the read overlaps.
pub fn get_ecov_adv_full(
    context: &UlRescallContext<'_>,
    v: u32,
    w: u32,
    bw: i64,
    diff_ec_ul: f64,
    dq: i64,
    mode: u64,
    mut contain_off: Option<&mut i64>,
) -> bool {
    let mut dt: i64 = -1;
    if let Some(out) = contain_off.as_deref_mut() {
        *out = 0;
    }
    if mode == 0 {
        let g = &context.unitig_graph.graph;
        for arc in g.arcs_of(v) {
            if arc.deleted || arc.v != w {
                continue;
            }
            dt = i64::from(arc.overlap);
            if let Some(out) = contain_off.as_deref_mut() {
                *out = i64::from(arc.occurrence);
                if u32::from(arc.occurrence) >= OU_MASK {
                    if let Some(hit) = get_ug_edge_src(context, (arc.ul >> 32) as u32, arc.v) {
                        let x = (hit.query_identifier) as usize;
                        if let Some(pos) = context.sources[x].iter().position(|h| std::ptr::eq(h, hit)) {
                            *out = i64::from(context.cc[x][pos]);
                        }
                    }
                }
            }
            break;
        }
    } else {
        let x = (v >> 1) as usize;
        for hit in &context.sources[x] {
            let qn = hit.query_identifier;
            let tn = hit.target_identifier;
            if tn != w >> 1 {
                continue;
            }
            let (result, arc) = ma_hit_to_arc(
                hit,
                context.reads.read_length(qn as usize) as i32,
                context.reads.read_length(tn as usize) as i32,
                context.maximum_hang as i32,
                context.maximum_hang_rate,
                context.minimum_overlap as i32,
            );
            if result < 0 {
                continue;
            }
            if (arc.ul >> 32) as u32 != v || arc.v != w {
                continue;
            }
            dt = i64::from(arc.overlap);
            if let Some(out) = contain_off.as_deref_mut() {
                let pos = context.sources[x]
                    .iter()
                    .position(|h| std::ptr::eq(h, hit))
                    .unwrap_or(0);
                *out = i64::from(context.cc[x][pos]);
            }
            break;
        }
    }
    if dt < 0 {
        return false;
    }
    let dif = (dq - dt).abs();
    let mut mm = (dq.max(dt) as f64 * diff_ec_ul) as i64;
    if mm < bw {
        mm = bw;
    }
    dif <= mm
}

/// `get_ug_edge_src` (`inter.cpp` 3381).
fn get_ug_edge_src<'a>(context: &'a UlRescallContext<'_>, uv: u32, uw: u32) -> Option<&'a crate::overlap::paf::PafHit> {
    let unitigs = &context.unitig_graph.unitigs;
    let uu = &unitigs[(uv >> 1) as usize];
    let wu = &unitigs[(uw >> 1) as usize];
    if uu.circular || wu.circular {
        return None;
    }
    let v = if uv & 1 != 0 { uu.start ^ 1 } else { uu.end ^ 1 };
    let w = if uw & 1 != 0 { wu.end } else { wu.start };
    for hit in &context.sources[(v >> 1) as usize] {
        let qn = hit.query_identifier;
        let tn = hit.target_identifier;
        if qn != v >> 1 || tn != w >> 1 {
            continue;
        }
        let (result, arc) = ma_hit_to_arc(
            hit,
            context.reads.read_length((v >> 1) as usize) as i32,
            context.reads.read_length((w >> 1) as usize) as i32,
            context.maximum_hang as i32,
            context.maximum_hang_rate,
            context.minimum_overlap as i32,
        );
        if result < 0 {
            continue;
        }
        if (arc.ul >> 32) as u32 != v || arc.v != w {
            continue;
        }
        return Some(hit);
    }
    None
}

/// `update_ava_het_site` (`inter.cpp` 2412): the het-site records of overlap
/// `oid`, compacted in place; `beg`/`end` delimit the scan.
pub fn update_ava_het_site(
    hap: &mut HaplotypeEvidenceAlloc,
    oid: u32,
    mut beg: Option<&mut usize>,
    mut end: Option<&mut usize>,
    is_srt: bool,
) -> usize {
    let n = hap.length;
    let mut start = match beg.as_deref() {
        Some(b) => *b,
        None => 0,
    };
    if let Some(e) = end.as_deref_mut() {
        *e = n;
    }
    if let Some(b) = beg.as_deref_mut() {
        *b = n;
    }
    if start < n && hap.list[start].overlap_identifier > oid {
        if let Some(e) = end {
            *e = start;
        }
        return 0;
    }
    let mut occ = 0usize;
    let mut need_srt = false;
    let mut k = start + 1;
    let mut l = start;
    while k <= n {
        if hap.list[l].overlap_identifier > oid {
            if let Some(e) = end.as_deref_mut() {
                *e = l;
            }
            break;
        }
        if k == n || hap.list[k].overlap_identifier != hap.list[l].overlap_identifier {
            if hap.list[l].overlap_identifier == oid {
                for i in l..k {
                    if hap.list[i].evidence_type != 1 {
                        continue;
                    }
                    let site = hap.list[i].overlap_site as usize;
                    let s = &hap.snp_stat[site];
                    if s.score == 1 && !(s.occ_0 < 2 || s.occ_1 < 2) {
                        if start + occ != i {
                            hap.list.swap(start + occ, i);
                        }
                        if occ > 0
                            && hap.list[start + occ].coverage < hap.list[start + occ - 1].coverage
                        {
                            need_srt = true;
                        }
                        occ += 1;
                    }
                }
                if let Some(b) = beg.as_deref_mut() {
                    *b = l;
                }
                if let Some(e) = end.as_deref_mut() {
                    *e = k;
                }
                break;
            }
            l = k;
        }
        k += 1;
    }
    let _ = start;
    if occ > 0 && is_srt && need_srt {
        hap.list[l..l + occ].sort_by_key(|e| e.coverage);
    }
    occ
}

/// `gl_chain_gen` (`inter.cpp` 2462): the `ul_ov_t` list of the current
/// overlaps; `rec_trans` selects whether `is_match == 2` overlaps join in.
pub fn gl_chain_gen(
    olist: &[OverlapRegion],
    unitigs: &[MaUtg],
    res: &mut Vec<UlOv>,
    rec_trans: bool,
    hap: &mut HaplotypeEvidenceAlloc,
) -> usize {
    let mut o2 = 0usize;
    let mut si = 0usize;
    let mut ei = 0usize;
    res.clear();
    for (k, region) in olist.iter().enumerate() {
        if region.is_match == 2 {
            o2 += 1;
        }
        if !rec_trans && region.is_match != 1 {
            continue;
        }
        if rec_trans && region.is_match != 1 && region.is_match != 2 {
            continue;
        }
        let rev = region.y_strand != 0;
        let tn = region.y_identifier;
        let (ts, te) = if rev {
            let len = unitigs[tn as usize].length;
            (len.wrapping_sub(region.y_position_end + 1), len.wrapping_sub(region.y_position_start))
        } else {
            (region.y_position_start, region.y_position_end + 1)
        };
        let mut sec = 0u32;
        if region.is_match == 2 {
            sec = update_ava_het_site(hap, k as u32, Some(&mut si), Some(&mut ei), true) as u32;
            assert!(sec > 0);
            si = ei;
        }
        res.push(UlOv {
            qn: k as u32,
            qs: region.x_position_start,
            qe: region.x_position_end + 1,
            tn,
            ts,
            te,
            sec,
            el: region.is_match == 1,
            rev,
        });
    }
    o2
}

/// `get_chain_x` (`inter.cpp` 2578): the query coordinate of target position
/// `q` along the chain's fake cigar.
fn get_chain_x(region: &OverlapRegion, q: i64) -> i64 {
    let o = &region.fake_cigar;
    let get_pos = |i: usize| -> i64 { o.gap_position(i) };
    let get_shift = |i: usize| -> i64 { o.gap_shift(i) };
    let mut x = get_pos(o.entries.len() - 1);
    let off = get_shift(o.entries.len() - 1);
    let mut y = x - region.x_position_start as i64 + region.y_position_start as i64 + off;
    if y == q {
        return x;
    }
    let mut lx = -1i64;
    let mut ly = -1i64;
    let mut hit = usize::MAX;
    for i in 0..o.entries.len() {
        x = get_pos(i);
        let off = get_shift(i);
        y = x - region.x_position_start as i64 + region.y_position_start as i64 + off;
        if q < y {
            lx = x;
            ly = y;
            hit = i;
            break;
        }
    }
    assert!(hit != 0 && hit != o.entries.len());
    x = get_pos(hit - 1);
    let off = get_shift(hit - 1);
    y = x - region.x_position_start as i64 + region.y_position_start as i64 + off;
    let mut y = ((((q - y) as f64) / ((ly - y) as f64)) * ((lx - x) as f64) + x as f64) as i64;
    if y < region.x_position_start as i64 {
        y = region.x_position_start as i64;
    }
    if y > region.x_position_end as i64 {
        y = region.x_position_end as i64;
    }
    y
}

/// `es_win_err` (`inter.cpp` 2630).
fn es_win_err(o: &OverlapRegion, win_len: i64, s: i64, e: i64) -> f64 {
    let windows = &o.window_cigar.windows;
    let os = (o.x_position_start as i64 / win_len) * win_len;
    let si = (s - os) / win_len;
    let ei = (e - os) / win_len;
    let mut t_err = 0i64;
    let mut t_len = 0i64;
    for k in (si + 1)..=(ei - 1) {
        let w = &windows[k as usize];
        t_len += i64::from(w.x_end) + 1 - i64::from(w.x_start);
        if w.y_end != -1 {
            t_err += i64::from(w.error);
        } else {
            t_err += i64::from(w.x_end) + 1 - i64::from(w.x_start);
        }
    }
    for k in [si, ei] {
        let w = &windows[k as usize];
        let max_s = s.max(i64::from(w.x_start));
        let min_e = (e.min(i64::from(w.x_end))) + 1;
        let ov = (min_e - max_s).max(0);
        t_len += ov;
        if w.y_end != -1 {
            t_err += (ov * i64::from(w.error)) / (i64::from(w.x_end) + 1 - i64::from(w.x_start));
        } else {
            t_err += ov;
        }
    }
    (t_err as f64) / (t_len as f64)
}

/// `utg_ct_t` (`Overlaps.h` 278).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct UtgCt {
    pub x: u32,
    pub s: u32,
    pub e: u32,
}

/// `gen_contain_chain` (`inter.cpp` 2662).
fn gen_contain_chain(
    unitigs: &[MaUtg],
    reads: &ReadStore,
    p: &UtgCt,
    o: &OverlapRegion,
    chains: &mut Vec<UlOv>,
    diff_ec_ul: f64,
    win_len: i64,
) -> bool {
    let tl = u64::from(unitigs[o.y_identifier as usize].length);
    let (mut y_s, mut y_e) = if o.y_strand != 0 {
        (tl as i64 - p.e as i64, tl as i64 - p.s as i64 - 1)
    } else {
        (i64::from(p.s), i64::from(p.e) - 1)
    };
    y_s = y_s.max(i64::from(o.y_position_start));
    y_e = y_e.min(i64::from(o.y_position_end));
    if y_s > y_e {
        return false;
    }
    let x_s = get_chain_x(o, y_s);
    let x_e = get_chain_x(o, y_e) + 1;
    if win_len > 0 && diff_ec_ul > 0.0 && es_win_err(o, win_len, x_s, x_e - 1) > diff_ec_ul {
        return false;
    }
    let (y_bs, y_be) = if o.y_strand != 0 {
        (tl as i64 - (y_e + 1), tl as i64 - y_s)
    } else {
        (y_s, y_e + 1)
    };
    let mut q_s = 0i64;
    let mut q_e = i64::from(p.e) - i64::from(p.s);
    if p.x & 1 != 0 {
        q_s += p.e as i64 - y_be;
        q_e -= y_bs - p.s as i64;
    } else {
        q_s += y_bs - p.s as i64;
        q_e -= p.e as i64 - y_be;
    }
    chains.push(UlOv {
        qn: o.x_identifier,
        qs: x_s as u32,
        qe: x_e as u32,
        tn: 0x8000_0000 | (p.x >> 1),
        ts: q_s as u32,
        te: q_e as u32,
        el: true,
        rev: (o.y_strand != 0) != (p.x & 1 != 0),
        ..UlOv::default()
    });
    true
}

/// `get_het_occ` (`inter.cpp` 2746).
fn get_het_occ(
    he_a: &[HaplotypeEvidence],
    c_k: usize,
    ylen: i64,
    p: &UtgCt,
    rev: bool,
) -> i64 {
    let mut occ = 0i64;
    if !rev {
        for k in (0..=c_k).rev() {
            let cov = i64::from(he_a[k].overlap_site);
            if cov >= i64::from(p.s) && cov < i64::from(p.e) {
                occ += 1;
            } else {
                break;
            }
        }
        for k in (c_k + 1)..he_a.len() {
            let cov = i64::from(he_a[k].overlap_site);
            if cov >= i64::from(p.s) && cov < i64::from(p.e) {
                occ += 1;
            } else {
                break;
            }
        }
    } else {
        for k in (0..=c_k).rev() {
            let ss = ylen - i64::from(he_a[k].overlap_site) - 1;
            if ss >= i64::from(p.s) && ss < i64::from(p.e) {
                occ += 1;
            } else {
                break;
            }
        }
        for k in (c_k + 1)..he_a.len() {
            let ss = ylen - i64::from(he_a[k].overlap_site) - 1;
            if ss >= i64::from(p.s) && ss < i64::from(p.e) {
                occ += 1;
            } else {
                break;
            }
        }
    }
    assert!(occ != 0);
    occ
}

/// `rescue_contain_ul_chains` (`inter.cpp` 2774).
#[allow(clippy::too_many_arguments)]
fn rescue_contain_ul_chains(
    unitigs: &[MaUtg],
    reads: &ReadStore,
    o: &OverlapRegion,
    he_a: &[HaplotypeEvidence],
    ct_a: &[UtgCt],
    chains: &mut Vec<UlOv>,
    diff_ec_ul: f64,
    win_len: i64,
    rescue_trans: bool,
) -> i64 {
    let mut t0 = 0i64;
    let tl = u64::from(unitigs[o.y_identifier as usize].length);
    let (ys, ye) = if o.y_strand == 0 {
        (i64::from(o.y_position_start), i64::from(o.y_position_end) + 1)
    } else {
        (
            tl as i64 - i64::from(o.y_position_end) - 1,
            tl as i64 - i64::from(o.y_position_start),
        )
    };
    let mut k = if o.y_strand == 0 { 0usize } else { he_a.len() };
    for p in ct_a.iter() {
        if (p.e as i64) <= ys {
            continue;
        }
        if (p.s as i64) >= ye {
            break;
        }
        let mut ff = true;
        if !he_a.is_empty() {
            if o.y_strand == 0 {
                while k < he_a.len() {
                    let cov = i64::from(he_a[k].overlap_site);
                    if cov >= i64::from(p.s) && cov < i64::from(p.e) {
                        ff = false;
                        break;
                    }
                    if cov >= i64::from(p.e) {
                        break;
                    }
                    k += 1;
                }
            } else {
                while k > 0 {
                    k -= 1;
                    let ss = tl as i64 - i64::from(he_a[k].overlap_site) - 1;
                    if ss >= i64::from(p.s) && ss < i64::from(p.e) {
                        ff = false;
                        break;
                    }
                    if ss >= i64::from(p.e) {
                        break;
                    }
                }
            }
        }
        if ff {
            if gen_contain_chain(unitigs, reads, p, o, chains, diff_ec_ul, win_len) {
                t0 += 1;
            }
        } else if rescue_trans {
            if gen_contain_chain(unitigs, reads, p, o, chains, diff_ec_ul, win_len) {
                t0 += 1;
                let last = chains.len() - 1;
                chains[last].el = false;
                chains[last].sec =
                    get_het_occ(he_a, k, tl as i64, p, o.y_strand != 0) as u32;
            }
        }
    }
    t0
}

/// `rescue_trans_ul_chains` (`inter.cpp` 2936).
#[allow(clippy::too_many_arguments)]
fn rescue_trans_ul_chains(
    unitigs: &[MaUtg],
    reads: &ReadStore,
    o: &OverlapRegion,
    he_a: &[HaplotypeEvidence],
    u: &MaUtg,
    chains: &mut Vec<UlOv>,
    diff_ec_ul: f64,
    win_len: i64,
    rescue_trans: bool,
    cis_occ: Option<&mut u64>,
) -> i64 {
    let mut t0 = 0i64;
    let mut cis = 0u64;
    let tl = u64::from(unitigs[o.y_identifier as usize].length);
    let (ys, ye) = if o.y_strand == 0 {
        (i64::from(o.y_position_start), i64::from(o.y_position_end) + 1)
    } else {
        (
            tl as i64 - i64::from(o.y_position_end) - 1,
            tl as i64 - i64::from(o.y_position_start),
        )
    };
    let mut k = if o.y_strand == 0 { 0usize } else { he_a.len() };
    let mut l = 0u32;
    for &packed in &u.reads {
        let p = UtgCt {
            x: (packed >> 32) as u32,
            s: l,
            e: l + reads.read_length((packed >> 33) as usize) as u32,
        };
        l = l.wrapping_add(packed as u32);
        if (p.e as i64) <= ys {
            continue;
        }
        if (p.s as i64) >= ye {
            break;
        }
        let mut ff = true;
        if !he_a.is_empty() {
            if o.y_strand == 0 {
                while k < he_a.len() {
                    let cov = i64::from(he_a[k].overlap_site);
                    if cov >= i64::from(p.s) && cov < i64::from(p.e) {
                        ff = false;
                        break;
                    }
                    if cov >= i64::from(p.e) {
                        break;
                    }
                    k += 1;
                }
            } else {
                while k > 0 {
                    k -= 1;
                    let ss = tl as i64 - i64::from(he_a[k].overlap_site) - 1;
                    if ss >= i64::from(p.s) && ss < i64::from(p.e) {
                        ff = false;
                        break;
                    }
                    if ss >= i64::from(p.e) {
                        break;
                    }
                }
            }
        }
        if ff {
            if gen_contain_chain(unitigs, reads, &p, o, chains, diff_ec_ul, win_len) {
                t0 += 1;
            }
        } else if rescue_trans {
            if gen_contain_chain(unitigs, reads, &p, o, chains, diff_ec_ul, win_len) {
                t0 += 1;
                cis += 1;
                let last = chains.len() - 1;
                chains[last].el = false;
                chains[last].sec =
                    get_het_occ(he_a, k, tl as i64, &p, o.y_strand != 0) as u32;
            }
        }
    }
    if let Some(out) = cis_occ {
        *out = cis;
    }
    t0
}

/// `dedup_sort_ul_ov_t` (`inter.cpp` 3036).
fn dedup_sort_ul_ov_t(a: &mut [UlOv]) -> usize {
    let a_n = a.len();
    let rr = 0.9f32;
    let mut i = 0usize;
    let mut l = 0usize;
    for k in 1..=a_n {
        if k == a_n || a[k].tn != a[l].tn {
            for z in l..k {
                let mut r = i as i64 - 1;
                while r >= 0 && a[r as usize].tn == a[z].tn {
                    let zr = r as usize;
                    if a[z].qn == a[zr].qn && a[z].tn == a[zr].tn && a[z].rev == a[zr].rev {
                        let qo = if a[z].qe.min(a[zr].qe) > a[z].qs.max(a[zr].qs) {
                            i64::from(a[z].qe.min(a[zr].qe) - a[z].qs.max(a[zr].qs))
                        } else {
                            0
                        };
                        let to = if a[z].te.min(a[zr].te) > a[z].ts.max(a[zr].ts) {
                            i64::from(a[z].te.min(a[zr].te) - a[z].ts.max(a[zr].ts))
                        } else {
                            0
                        };
                        if (qo as f32) >= ((a[zr].qe - a[zr].qs) as f32) * rr
                            && (qo as f32) >= ((a[z].qe - a[z].qs) as f32) * rr
                            && (to as f32) >= ((a[zr].te - a[zr].ts) as f32) * rr
                            && (to as f32) >= ((a[z].te - a[z].ts) as f32) * rr
                        {
                            break;
                        }
                    }
                    r -= 1;
                }
                if r >= 0 && a[r as usize].tn == a[z].tn {
                    if a[z].el {
                        a[r as usize].el = true;
                    }
                    continue;
                }
                a[i] = a[z];
                i += 1;
            }
            l = k;
        }
    }
    i
}

/// `check_contain_pair` (`inter.cpp` 3082).
fn check_contain_pair(context: &UlRescallContext<'_>, x: u32, y: u32, check_el: bool) -> bool {
    for hit in &context.sources[x as usize] {
        if check_el && hit.exact_region_byte == 0 {
            continue;
        }
        let qn = hit.query_identifier;
        let tn = hit.target_identifier;
        if tn != y {
            continue;
        }
        let (result, _arc) = ma_hit_to_arc(
            hit,
            context.reads.read_length(qn as usize) as i32,
            context.reads.read_length(tn as usize) as i32,
            context.maximum_hang as i32,
            context.maximum_hang_rate,
            context.minimum_overlap as i32,
        );
        // `r == MA_HT_QCONT || r == MA_HT_TCONT`
        if result == crate::assembly::hit::MA_HT_QCONT || result == crate::assembly::hit::MA_HT_TCONT {
            return true;
        }
    }
    false
}

/// `get_overlap_region_sub_err` (`inter.cpp` 5295): the error of overlap
/// `o`'s cigar right of query position `qs`, with the running cursor `it`.
pub fn get_overlap_region_sub_err(
    o: &OverlapRegion,
    it: &mut crate::correction::ul_align::RtraceIter,
    qs: i64,
    sec_err: &mut i64,
) -> i64 {
    let windows = &o.window_cigar.windows;
    if windows.is_empty() {
        *sec_err = 0;
        return 0;
    }
    if it.k == i32::MAX as i64 {
        it.k = windows.len() as i64 - 1;
        it.cur_qoff = i64::from(o.x_position_end) + 1;
        it.qoff = i64::from(windows[it.k as usize].x_end);
        it.werr = 0;
        it.werr0 = 0;
    }
    assert!(qs <= it.cur_qoff);
    it.cur_qoff = qs;
    *sec_err = it.werr0;
    if qs == it.qoff {
        *sec_err = it.werr0;
        return it.werr;
    }
    while it.k >= 0 && qs < i64::from(windows[it.k as usize].x_end) {
        let w = &windows[it.k as usize];
        if qs >= i64::from(w.x_start) && qs < i64::from(w.x_end) {
            *sec_err = 0;
            if w.cigar_length > 0 {
                let rr = (i64::from(w.x_end) - qs) as f64
                    / (i64::from(w.x_end) - i64::from(w.x_start)) as f64;
                *sec_err = (rr * w.cigar_length as f64) as i64;
            }
            *sec_err += it.werr0;
            let terr = if w.cigar_length > 0 {
                i64::from(w.x_end) - qs
            } else {
                0
            };
            return it.werr + terr;
        }
        if w.cigar_length > 0 {
            it.werr += i64::from(w.x_end) - i64::from(w.x_start);
            it.werr0 += i64::from(w.cigar_length);
        }
        it.qoff = i64::from(w.x_start);
        it.k -= 1;
    }
    *sec_err = it.werr0;
    it.werr
}

/// `get_ecov_el` (`inter.cpp` 5341): like [`get_ecov_adv_full`] but with the
/// graph arc's `el` bit (or the hit's `el` byte) on the side channel.
#[allow(clippy::too_many_arguments)]
pub fn get_ecov_el(
    context: &UlRescallContext<'_>,
    v: u32,
    w: u32,
    bw: i64,
    diff_ec_ul: f64,
    dq: i64,
    mode: u64,
    el: &mut i64,
) -> bool {
    let mut dt: i64 = -1;
    *el = 0;
    if mode == 0 {
        let g = &context.unitig_graph.graph;
        for arc in g.arcs_of(v) {
            if arc.deleted || arc.v != w {
                continue;
            }
            dt = i64::from(arc.overlap);
            *el = i64::from(arc.strong as i32);
            break;
        }
    } else {
        let x = (v >> 1) as usize;
        for (z, hit) in context.sources[x].iter().enumerate() {
            let qn = hit.query_identifier;
            let tn = hit.target_identifier;
            if tn != w >> 1 {
                continue;
            }
            let (result, arc) = ma_hit_to_arc(
                hit,
                context.reads.read_length(qn as usize) as i32,
                context.reads.read_length(tn as usize) as i32,
                context.maximum_hang as i32,
                context.maximum_hang_rate,
                context.minimum_overlap as i32,
            );
            if result < 0 {
                continue;
            }
            if (arc.ul >> 32) as u32 != v || arc.v != w {
                continue;
            }
            dt = i64::from(arc.overlap);
            *el = i64::from(context.cc[x][z] & 1);
            break;
        }
    }
    if dt < 0 {
        return false;
    }
    let dif = (dq - dt).abs();
    let mut mm = (dq.max(dt) as f64 * diff_ec_ul) as i64;
    if mm < bw {
        mm = bw;
    }
    dif <= mm
}

/// `cal_gl_chain_lin_sc` (`inter.cpp` 5391).
#[allow(clippy::too_many_arguments)]
fn cal_gl_chain_lin_sc(
    a: &[UlOv],
    ai: usize,
    aj: usize,
    tc: &mut crate::correction::ul_align::RtraceIter,
    ol: &[OverlapRegion],
    context: &UlRescallContext<'_>,
    bw: i64,
    diff_ec_ul: f64,
    mode: u64,
    trans_sc: i64,
    sec_sec: i64,
    f: &[i32],
) -> i64 {
    let li = &a[ai];
    let lj = &a[aj];
    if lj.qs >= li.qs {
        return i32::MIN as i64;
    }
    let li_v = (li.tn << 1) | u32::from(li.rev);
    let lj_v = (lj.tn << 1) | u32::from(lj.rev);
    let qo = infer_rovlp_unitigs(li, lj, &context.unitig_graph.unitigs) as i64;
    let mut el = 0i64;
    let mut sec_err = 0i64;
    if get_ecov_el(context, li_v ^ 1, lj_v ^ 1, bw, diff_ec_ul, qo, mode, &mut el) {
        let trans_l = get_overlap_region_sub_err(&ol[li.qn as usize], tc, i64::from(lj.qe), &mut sec_err);
        let mut sc = i64::from(f[aj]) + i64::from(li.qe) - i64::from(lj.qe)
            - trans_l * trans_sc
            - sec_err * sec_sec;
        if el == 0 && (trans_l > 0 || sec_err > 0) && lj.qe > li.qs {
            let mut tr = crate::correction::ul_align::RtraceIter::default();
            let sc0 = sc;
            sc = i64::from(f[aj]) + aln_sc(&ol[li.qn as usize], trans_sc, sec_sec);
            let mut sec_err2 = 0i64;
            let trans_l2 = get_overlap_region_sub_err(&ol[lj.qn as usize], &mut tr, i64::from(li.qs), &mut sec_err2);
            sc -= (i64::from(lj.qe) - i64::from(li.qs)) - trans_l2 * trans_sc - sec_err2 * sec_sec;
            if sc < sc0 {
                sc = sc0;
            }
        }
        return sc;
    }
    i32::MIN as i64
}

/// `gl_chain_lin` (`inter.cpp` 5431): the linear DP over the overlap list;
/// chains are dumped into `ex` with `res` rewritten to one entry per chain.
#[allow(clippy::too_many_arguments)]
pub fn gl_chain_lin(
    res: &mut Vec<UlOv>,
    ol: &[OverlapRegion],
    ex: &mut Vec<UlOv>,
    context: &UlRescallContext<'_>,
    bw: i64,
    diff_ec_ul: f64,
    qlen: i64,
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
    dp: &mut ChainData,
    trans_sc: i64,
    mode: u64,
    need_srt: bool,
) -> i64 {
    if res.is_empty() {
        return 0;
    }
    let res_n = res.len();
    dp.resize(res_n);
    let mut tc = crate::correction::ul_align::RtraceIter::default();
    let f: &mut Vec<i32> = &mut dp.score;
    let p: &mut Vec<i64> = &mut dp.predecessor;
    let t: &mut Vec<i64> = &mut dp.scratch;
    let c_n: &mut Vec<i32> = &mut dp.occurrence;
    let c_sc: &mut Vec<i32> = &mut dp.self_length;

    if need_srt {
        radix_sort_by_key_32(res, &|entry: &UlOv| entry.qe);
        let mut j = 0usize;
        for i in 1..=res_n {
            if i == res_n || res[i].qe != res[j].qe {
                if i - j > 1 {
                    radix_sort_by_key_32(&mut res[j..i], &|entry: &UlOv| entry.qs);
                }
                j = i;
            }
        }
    }
    for entry in t.iter_mut().take(res_n) {
        *entry = 0;
    }
    let g = &context.unitig_graph.graph;
    let mut st = 0i64;
    let mut plus = 0i64;
    let mut max_ii = -1i64;
    for i in 0..res_n {
        let li = res[i];
        let mm_ovlp = if mode != 0 {
            max_ovlp_src(context, ((li.tn << 1) | u32::from(li.rev)) ^ 1)
        } else {
            max_ovlp(g, ((li.tn << 1) | u32::from(li.rev)) ^ 1)
        };
        let mut x = ((i64::from(li.qs) + mm_ovlp) as f64 * diff_ec_ul) as i64;
        if x < bw {
            x = bw;
        }
        x += i64::from(li.qs) + mm_ovlp;
        if x > qlen + 1 {
            x = qlen + 1;
        }
        let x = find_ul_ov_max(&res[..i], (x + i64::from(crate::assembly::ul_rescall::G_CHAIN_INDEL)) as u32);
        let csc = aln_sc(&ol[li.qn as usize], trans_sc, UG_TRANS_ERR_W);
        let mut mm_sc = csc;
        let mut mm_idx = -1i64;
        let mut n_skip = 0i64;
        let mut end_j = -1i64;
        tc = crate::correction::ul_align::RtraceIter::default();
        if (x - st) > max_iter {
            st = x - max_iter;
        }
        let mut j = x;
        while j >= st {
            let lj = res[j as usize];
            if i64::from(lj.qe) + i64::from(crate::assembly::ul_rescall::G_CHAIN_INDEL) <= i64::from(li.qs) {
                break;
            }
            let sc = cal_gl_chain_lin_sc(res, i, j as usize, &mut tc, ol, context, bw, diff_ec_ul, mode, trans_sc, UG_TRANS_ERR_W, f);
            if sc == i32::MIN as i64 {
                j -= 1;
                continue;
            }
            if sc > mm_sc {
                mm_sc = sc;
                mm_idx = j;
                if n_skip > 0 {
                    n_skip -= 1;
                }
            } else if t[j as usize] == i as i64 {
                n_skip += 1;
                if n_skip > max_skip {
                    break;
                }
            }
            let pj = p[j as usize];
            if pj >= 0 {
                t[pj as usize] = i as i64;
            }
            j -= 1;
        }
        end_j = j;
        if max_ii < 0 || i64::from(res[i].qe) > i64::from(res[max_ii as usize].qe) + max_dis {
            let mut max = i32::MIN as i64;
            max_ii = -1;
            let mut j = i as i64 - 1;
            while j >= st && i64::from(res[i].qe) <= max_dis + i64::from(res[j as usize].qe) {
                if max < i64::from(f[j as usize]) {
                    max = i64::from(f[j as usize]);
                    max_ii = j;
                }
                j -= 1;
            }
        }
        if max_ii >= 0 && max_ii < end_j {
            let lj = res[max_ii as usize];
            if i64::from(lj.qe) + i64::from(crate::assembly::ul_rescall::G_CHAIN_INDEL) > i64::from(li.qs)
                && lj.qs < li.qs
            {
                let sc = cal_gl_chain_lin_sc(res, i, max_ii as usize, &mut tc, ol, context, bw, diff_ec_ul, mode, trans_sc, UG_TRANS_ERR_W, f);
                if sc != i32::MIN as i64 && sc > mm_sc {
                    mm_sc = sc;
                    mm_idx = max_ii;
                }
            }
        }
        if mm_sc < 0 {
            mm_sc = csc;
            mm_idx = -1;
        }
        f[i] = mm_sc as i32;
        p[i] = mm_idx;
        if max_ii < 0
            || (i64::from(res[i].qe) <= max_dis + i64::from(res[max_ii as usize].qe)
                && i64::from(f[max_ii as usize]) < i64::from(f[i]))
        {
            max_ii = i as i64;
        }
        if mm_sc < plus {
            plus = mm_sc;
        }
    }
    for i in 0..res_n {
        f[i] -= plus as i32;
        t[i] = ((f[i] as i64) << 32) + (i as i64) << 1;
    }
    let mut n_v = 0usize;
    let mut n_u = 0usize;
    let mut plus2 = 0i64;
    radix_sort_i64(&mut t[..res_n]);
    for k in (0..res_n).rev() {
        let n_v0 = n_v;
        let mut i = ((t[k] as u32) >> 1) as i64;
        while i >= 0 && (t[i as usize] & 1) == 0 {
            ex.push(res[i as usize]);
            n_v += 1;
            t[i as usize] |= 1;
            i = p[i as usize];
        }
        if n_v0 == n_v {
            continue;
        }
        let sc = if i < 0 {
            t[k] >> 32
        } else {
            (t[k] >> 32) - i64::from(f[i as usize])
        };
        c_n[n_u] = (n_v - n_v0) as i32;
        c_sc[n_u] = sc as i32;
        n_u += 1;
        if sc < plus2 {
            plus2 = sc;
        }
    }
    n_v = 0;
    let mut n_v0;
    for k in 0..n_u {
        n_v0 = n_v;
        n_v += c_n[k] as usize;
        res[k].qn = (c_sc[k] - plus2 as i32) as u32;
        res[k].ts = n_v0 as u32;
        res[k].te = n_v as u32;
        let rev_n = (c_n[k] as u32) >> 1;
        res[k].qs = u32::MAX;
        res[k].qe = ex[n_v0].qe;
        for i in 0..rev_n as usize {
            ex.swap(n_v0 + i, n_v - i - 1);
            if res[k].qs > ex[n_v0 + i].qs {
                res[k].qs = ex[n_v0 + i].qs;
            }
            if res[k].qs > ex[n_v - i - 1].qs {
                res[k].qs = ex[n_v - i - 1].qs;
            }
            ex[n_v0 + i].sec = SEC_MODE;
            ex[n_v - i - 1].sec = SEC_MODE;
        }
        if c_n[k] & 1 != 0 {
            let i = rev_n as usize;
            if res[k].qs > ex[n_v0 + i].qs {
                res[k].qs = ex[n_v0 + i].qs;
            }
            ex[n_v0 + i].sec = SEC_MODE;
        }
    }
    res.truncate(n_u);
    radix_sort_by_key_32(res, &|entry: &UlOv| entry.qn);
    n_v as i64
}

/// `gen_gl_aln` (`inter.cpp` 5223): the `ul_ov_t` list used by `gl_chain`,
/// `sec` = the alignment length.
pub fn gen_gl_aln(olist: &[OverlapRegion], unitigs: &[MaUtg], res: &mut Vec<UlOv>) {
    res.clear();
    res.reserve(olist.len());
    for (k, region) in olist.iter().enumerate() {
        assert!(region.overlap_length >= region.alignment_length);
        let rev = region.y_strand != 0;
        let tn = region.y_identifier;
        let (ts, te) = if rev {
            let len = unitigs[tn as usize].length;
            (len.wrapping_sub(region.y_position_end + 1), len.wrapping_sub(region.y_position_start))
        } else {
            (region.y_position_start, region.y_position_end + 1)
        };
        res.push(UlOv {
            qn: k as u32,
            qs: region.x_position_start,
            qe: region.x_position_end + 1,
            tn,
            ts,
            te,
            sec: region.alignment_length,
            el: region.is_match == 1,
            rev,
        });
    }
}

/// `gen_gg_aln` (`inter.cpp` 5246): the seed `mg_lchain_t` list.
pub fn gen_gg_aln(olist: &[OverlapRegion], unitigs: &[MaUtg], trans_sc: i64, res: &mut Vec<MgLchain>) {
    res.clear();
    res.reserve(olist.len());
    for (k, region) in olist.iter().enumerate() {
        let v = (region.y_identifier << 1) | region.y_strand;
        let (rs, re) = if v & 1 != 0 {
            let len = unitigs[(v >> 1) as usize].length;
            (len.wrapping_sub(region.y_position_end + 1) as i32, len.wrapping_sub(region.y_position_start) as i32)
        } else {
            (region.y_position_start as i32, (region.y_position_end + 1) as i32)
        };
        res.push(MgLchain {
            v,
            off: k as i32,
            score: aln_sc(region, trans_sc, UG_TRANS_ERR_W) as i32,
            qs: region.x_position_start as i32,
            qe: (region.x_position_end + 1) as i32,
            rs,
            re,
            ..MgLchain::default()
        });
    }
}

/// `regen_ul_ov_t_lst` (`inter.cpp` 7000): the `ul_ov_t` list used by
/// `gl_chain_flter`, `sec` = the substitution count.
pub fn regen_ul_ov_t_lst(unitigs: &[MaUtg], olist: &[OverlapRegion], idx: &mut Vec<UlOv>) {
    idx.clear();
    idx.reserve(olist.len());
    for (k, region) in olist.iter().enumerate() {
        let rev = region.y_strand != 0;
        let tn = region.y_identifier;
        let (ts, te) = if rev {
            let len = unitigs[tn as usize].length;
            (len.wrapping_sub(region.y_position_end + 1), len.wrapping_sub(region.y_position_start))
        } else {
            (region.y_position_start, region.y_position_end + 1)
        };
        idx.push(UlOv {
            qn: k as u32,
            qs: region.x_position_start,
            qe: region.x_position_end + 1,
            tn,
            ts,
            te,
            sec: region.non_homopolymer_errors,
            el: true,
            rev,
        });
    }
}

/// `g_adjacent_dis` (`inter.cpp` 4575).
fn g_adjacent_dis(g: &Asg, v: u32, w: u32) -> i32 {
    for arc in g.arcs_of(v) {
        if arc.deleted || arc.v != w {
            continue;
        }
        return arc.overlap as i32;
    }
    -1
}

/// `l2g_res_chain_sc` (`inter.cpp` 4580): converts the best chain into the
/// `mg_lchain_t` group layout (a `(uint32_t)-1` header plus the members).
pub fn l2g_res_chain_sc(
    unitigs: &[MaUtg],
    g: &Asg,
    a: &[UlOv],
    gchains: &mut Vec<MgLchain>,
) -> i64 {
    let mut a_n = a.len();
    if a_n == 0 {
        return 0;
    }
    a_n += 1;
    gchains.clear();
    gchains.resize(a_n, MgLchain::default());
    gchains[0] = MgLchain {
        cnt: (a_n - 1) as i32,
        v: u32::MAX,
        ..MgLchain::default()
    };
    let mut l = 0i64;
    for k in 1..a_n {
        let m = k - 1;
        gchains[k] = MgLchain {
            v: (a[m].tn << 1) | u32::from(a[m].rev),
            dist_pre: -1,
            off: a[m].qn as i32,
            score: a[m].sec as i32,
            qs: a[m].qs as i32,
            qe: a[m].qe as i32,
            rs: a[m].ts as i32,
            re: a[m].te as i32,
            ..MgLchain::default()
        };
        if k > 1 {
            gchains[k - 1].dist_pre = g_adjacent_dis(g, gchains[k].v ^ 1, gchains[k - 1].v ^ 1);
            assert!(gchains[k - 1].dist_pre >= 0);
            l += i64::from(g.sequences[(gchains[k - 1].v >> 1) as usize].length)
                + i64::from(gchains[k - 1].dist_pre);
        }
    }
    let _ = unitigs;
    1
}

/// `hc_target_len` (`inter.cpp` 5588).
fn hc_target_len(g: &Asg, s: &MgLchain, e: &MgLchain) -> i64 {
    let ql = i64::from(s.qs) - i64::from(e.qe);
    let sts = if s.v & 1 != 0 {
        i64::from(g.sequences[(s.v >> 1) as usize].length) - i64::from(s.re)
    } else {
        i64::from(s.rs)
    };
    let tp = i64::from(g.sequences[(s.v >> 1) as usize].length) - sts;
    let ete = if e.v & 1 != 0 {
        i64::from(g.sequences[(e.v >> 1) as usize].length) - i64::from(e.rs)
    } else {
        i64::from(e.re)
    };
    let tm = i64::from(g.sequences[(e.v >> 1) as usize].length) - ete;
    ql + tp - tm
}

/// `set_ul_ov_t_by_mg_lchain_t` (`inter.cpp` 5646).
pub fn set_ul_ov_t_by_mg_lchain_t(u: &mut UlOv, l: &MgLchain) {
    u.tn = l.v >> 1;
    u.rev = l.v & 1 != 0;
    u.ts = l.rs as u32;
    u.te = l.re as u32;
    u.qs = l.qs as u32;
    u.qe = l.qe as u32;
}

/// `mg_log2` (`inter.cpp` 654): the bit-twiddling approximation of log2.
fn mg_log2(x: f32) -> f32 {
    let z = x.to_bits();
    let mut log_2 = (((z >> 23) & 255) as i32 - 128) as f32;
    let z2 = (z & !(255u32 << 23)) + (127u32 << 23);
    let zf = f32::from_bits(z2);
    log_2 += (-0.34484843f32 * zf + 2.02466578f32) * zf - 0.67487759f32;
    log_2
}

/// `cal_gchain_sc_adv` (`inter.cpp` 6036): scores one shortest-path candidate.
#[allow(clippy::too_many_arguments)]
fn cal_gchain_sc_adv(
    context: &UlRescallContext<'_>,
    ol: &[OverlapRegion],
    dj: &MgPathDst,
    li: &MgLchain,
    ui: &UlOv,
    lc: &[MgLchain],
    f: &[i64],
    b_w: i64,
    diff_thre: f64,
    chn_pen_gap: f32,
    tc: &mut crate::correction::ul_align::RtraceIter,
    trans_sc: i64,
    sec_sec: i64,
) -> i32 {
    if dj.n_path == 0 {
        return i32::MIN;
    }
    let mut gap = i64::from(dj.dist) - i64::from(dj.target_dist);
    if gap < 0 {
        gap = -gap;
    }
    if gap > ((dj.target_dist as f64) * diff_thre) as i64 && gap > b_w {
        return i32::MIN;
    }
    let lj = &lc[dj.meta as usize];
    let mut sec_err = 0i64;
    let mut trans_l = get_overlap_region_sub_err(&ol[li.off as usize], tc, i64::from(lj.qe), &mut sec_err);
    let mut sc = i64::from(li.qe) - i64::from(lj.qe) - trans_l * trans_sc - sec_err * sec_sec;
    sc += f[dj.meta as usize];
    if (trans_l > 0 || sec_err > 0) && lj.qe > li.qs {
        let mut uj = UlOv::default();
        set_ul_ov_t_by_mg_lchain_t(&mut uj, lj);
        let qo = infer_rovlp_unitigs(ui, &uj, &context.unitig_graph.unitigs);
        let mut el = 0i64;
        if !get_ecov_el(
            context,
            li.v ^ 1,
            lj.v ^ 1,
            b_w,
            crate::assembly::ul_rescall::N_GCHAIN_RATE,
            i64::from(qo),
            0,
            &mut el,
        ) || el == 0
        {
            let mut tr = crate::correction::ul_align::RtraceIter::default();
            let sc0 = sc;
            sc = f[dj.meta as usize] + i64::from(li.score);
            let mut sec_err2 = 0i64;
            trans_l = get_overlap_region_sub_err(&ol[lj.off as usize], &mut tr, i64::from(li.qs), &mut sec_err2);
            sc -= (i64::from(lj.qe) - i64::from(li.qs)) - trans_l * trans_sc - sec_err2 * sec_sec;
            if sc < sc0 {
                sc = sc0;
            }
        }
    }
    let lin_pen = chn_pen_gap * gap as f32;
    let log_pen = if gap >= 2 { mg_log2(gap as f32) } else { 0.0 };
    sc -= i64::from((lin_pen + log_pen) as i32);
    sc as i32
}

/// `gl_chain_graph_adv` (`inter.cpp` 6086): the graph DP over the seed
/// chains, with `hc_shortest_k` resolving multi-edge connections.
#[allow(clippy::too_many_arguments)]
pub fn gl_chain_graph_adv(
    ol: &[OverlapRegion],
    context: &UlRescallContext<'_>,
    lc: &mut Vec<MgLchain>,
    sw: &mut Vec<MgLchain>,
    dst: &mut Vec<MgPathDst>,
    out: &mut Vec<usize>,
    path: &mut Vec<MgPathv>,
    qlen: i64,
    bw: i64,
    diff_thre: f64,
    bf: &mut Vec<u64>,
    dp: &mut ChainData,
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
    trans_sc: i64,
    sec_sec: i64,
    need_srt: bool,
) -> i64 {
    bf.clear();
    if lc.is_empty() {
        return 0;
    }
    let lc_n = lc.len();
    let g = &context.unitig_graph.graph;
    let mut srt: Vec<u64> = Vec::with_capacity(lc_n);
    let mut n_ext = 0usize;
    if !need_srt {
        srt.clear();
        for (i, r) in lc.iter_mut().enumerate() {
            r.dist_pre = -1;
            let isolated = (r.re as i64) < i64::from(g.sequences[(r.v >> 1) as usize].length) && r.rs > 0;
            if !isolated {
                n_ext += 1;
                srt.push(((r.qe as u64) << 32) | i as u64);
            }
        }
        if n_ext < lc_n {
            for (i, r) in lc.iter().enumerate() {
                let isolated = (r.re as i64) < i64::from(g.sequences[(r.v >> 1) as usize].length) && r.rs > 0;
                if isolated {
                    srt.push(((r.qe as u64) << 32) | i as u64 | (1u64 << 63));
                }
            }
        }
        debug_assert_eq!(srt.len(), lc_n);
    } else {
        for (i, r) in lc.iter_mut().enumerate() {
            r.dist_pre = -1;
            let isolated = (r.re as i64) < i64::from(g.sequences[(r.v >> 1) as usize].length) && r.rs > 0;
            if !isolated {
                n_ext += 1;
            }
            srt.push(((r.qe as u64) << 32) | i as u64 | (u64::from(isolated) << 63));
        }
        radix_sort_u64(&mut srt);
        let mut j = 0usize;
        for i in 1..=lc_n {
            if i == lc_n || (srt[i] >> 32) != (srt[j] >> 32) {
                if i - j > 1 {
                    for x in j..i {
                        srt[x] = (srt[x] << 32) >> 32;
                        srt[x] |= (lc[(srt[x] as u32) as usize].qs as u64) << 32;
                    }
                    radix_sort_u64(&mut srt[j..i]);
                }
                j = i;
            }
        }
    }
    if n_ext != lc_n || need_srt {
        let reordered: Vec<MgLchain> = srt[..lc_n]
            .iter()
            .map(|entry| lc[(*entry as u32) as usize])
            .collect();
        *sw = reordered.clone();
        lc.copy_from_slice(&reordered);
    }
    dp.resize(lc_n);
    let mut tc = crate::correction::ul_align::RtraceIter::default();
    let p: &mut Vec<i32> = &mut dp.score;
    let f: &mut Vec<i64> = &mut dp.predecessor;
    let t: &mut Vec<i64> = &mut dp.scratch;
    for entry in t.iter_mut().take(n_ext) {
        *entry = 0;
    }
    let mut plus = 0i64;
    for i in 0..n_ext {
        let li = lc[i];
        let mut ui = UlOv::default();
        set_ul_ov_t_by_mg_lchain_t(&mut ui, &li);
        let mm_ovlp = max_ovlp(g, li.v ^ 1);
        let mut x = ((i64::from(li.qs) + mm_ovlp) as f64 * diff_thre) as i64;
        if x < bw {
            x = bw;
        }
        x += i64::from(li.qs) + mm_ovlp;
        if x > qlen + 1 {
            x = qlen + 1;
        }
        let x = find_mg_lchain_max(&lc[..i], (x + i64::from(crate::assembly::ul_rescall::G_CHAIN_INDEL)) as i32);
        let mut n_skip = 0i64;
        let mut is_f = false;
        dst.clear();
        let mut max_target_dist = -1i64;
        let mut j = x;
        while j >= 0 {
            let lj = lc[j as usize];
            if lj.qs >= li.qs {
                j -= 1;
                continue;
            }
            let target_dist = hc_target_len(g, &li, &lj);
            if target_dist < 0 {
                j -= 1;
                continue;
            }
            dst.push(MgPathDst {
                v: lj.v ^ 1,
                meta: j as u32,
                qlen: li.qs - lj.qe,
                target_dist: target_dist as i32,
                ..MgPathDst::default()
            });
            if max_target_dist < target_dist {
                max_target_dist = target_dist;
            }
            if t[j as usize] == i as i64 {
                n_skip += 1;
                if n_skip > max_skip {
                    break;
                }
            }
            let pj = p[j as usize];
            if pj >= 0 {
                t[pj as usize] = i as i64;
            }
            if !is_f && i64::from(lj.qe) + i64::from(crate::assembly::ul_rescall::G_CHAIN_INDEL) > i64::from(li.qs) {
                let mut uj = UlOv::default();
                set_ul_ov_t_by_mg_lchain_t(&mut uj, &lj);
                let qo = infer_rovlp_unitigs(&ui, &uj, &context.unitig_graph.unitigs);
                if get_ecov_adv_full(context, li.v ^ 1, lj.v ^ 1, bw, crate::assembly::ul_rescall::N_GCHAIN_RATE, i64::from(qo), 0, None) {
                    is_f = true;
                    if n_skip > 0 {
                        n_skip -= 1;
                    }
                    if n_skip < max_skip >> 1 {
                        n_skip = max_skip >> 1;
                    }
                }
            }
            j -= 1;
        }
        let mut max_f = i64::from(li.score);
        let mut max_j = -1i64;
        let mut max_d = -1i64;
        let mut max_inner = false;
        let mut max_hash = 0u32;
        if !dst.is_empty() {
            max_target_dist = (max_target_dist as f64 * (1.0 + diff_thre)) as i64;
            if max_target_dist < bw {
                max_target_dist = bw;
            }
            hc_shortest_k(
                g,
                li.v ^ 1,
                dst,
                max_target_dist as i32,
                MG_MAX_SHORT_K as i32,
                bf,
                out,
                path,
                true,
                false,
                0.0,
            );
            let dst_n = dst.len();
            tc = crate::correction::ul_align::RtraceIter::default();
            for dj in dst.iter().take(dst_n) {
                if dj.n_path == 0 {
                    continue;
                }
                let sc = cal_gchain_sc_adv(context, ol, dj, &li, &ui, lc, f, bw, diff_thre, W_CHN_PEN_GAP, &mut tc, trans_sc, sec_sec);
                if sc == i32::MIN {
                    continue;
                }
                if i64::from(sc) > max_f {
                    max_f = i64::from(sc);
                    max_j = dj.meta as i64;
                    max_d = i64::from(dj.dist);
                    max_hash = dj.hash;
                    max_inner = dj.inner;
                }
            }
        }
        if max_f < 0 {
            max_f = i64::from(li.score);
            max_j = -1;
        }
        f[i] = max_f;
        p[i] = max_j as i32;
        lc[i].dist_pre = max_d as i32;
        lc[i].hash_pre = max_hash;
        lc[i].inner_pre = max_inner;
        if max_f < plus {
            plus = max_f;
        }
    }
    for i in n_ext..lc_n {
        let li = lc[i];
        let max_f = i64::from(li.score);
        f[i] = max_f;
        p[i] = -1;
        lc[i].dist_pre = -1;
        lc[i].hash_pre = 0;
        lc[i].inner_pre = false;
        if max_f < plus {
            plus = max_f;
        }
    }
    for i in 0..lc_n {
        f[i] -= plus;
        t[i] = (f[i] << 32) + (i as i64) << 1;
    }
    sw.clear();
    sw.reserve(lc_n);
    bf.clear();
    bf.reserve(lc_n);
    let mut n_v = 0usize;
    let mut n_u = 0usize;
    radix_sort_i64(&mut t[..lc_n]);
    let mut plus2 = 0i64;
    let mut u: Vec<u64> = Vec::new();
    for k in (0..lc_n).rev() {
        let n_v0 = n_v;
        let mut i = ((t[k] as u32) >> 1) as i64;
        while i >= 0 && (t[i as usize] & 1) == 0 {
            sw.push(lc[i as usize]);
            n_v += 1;
            t[i as usize] |= 1;
            i = p[i as usize] as i64;
        }
        if n_v0 == n_v {
            continue;
        }
        let mut sc = if i < 0 {
            t[k] >> 32
        } else {
            (t[k] >> 32) - f[i as usize]
        };
        if sc < plus2 {
            plus2 = sc;
        }
        let ff;
        if sc >= 0 {
            ff = 0x8000_0000_0000_0000u64;
        } else {
            ff = 0;
            sc = -sc;
        }
        u.push(((sc as u64) << 32) | ((n_v - n_v0) as u64) | ff);
    }
    let mut m_idx = -1i64;
    let mut m_sc = -1i64;
    let mut k = 0usize;
    for (i, entry) in u.iter_mut().enumerate() {
        let mut sc;
        if *entry & 0x8000_0000_0000_0000 != 0 {
            *entry -= 0x8000_0000_0000_0000;
            sc = (*entry >> 32) as i64;
        } else {
            sc = (*entry >> 32) as i64;
            sc = -sc;
        }
        sc -= plus2;
        *entry <<= 32;
        *entry >>= 32;
        *entry |= (sc as u64) << 32;
        let k0 = k;
        let ni = (*entry as u32) as usize;
        for j in 0..ni {
            lc[k] = sw[k0 + (ni - j - 1)];
            k += 1;
        }
        if m_idx < 0 || m_sc < (*entry >> 32) as i64 {
            m_idx = i as i64;
            m_sc = (*entry >> 32) as i64;
        }
    }
    debug_assert_eq!(k, n_v);
    *bf = u;
    m_idx
}

/// `comput_sc_partial_cigar` (`inter.cpp` 6637).
fn comput_sc_partial_cigar(
    sc: i64,
    ol: i64,
    err_sc_r: f64,
    z: &OverlapRegion,
    wi: &mut i64,
    werr: &mut i64,
) -> i64 {
    let mut k = *wi;
    let wn = z.window_cigar.windows.len() as i64;
    let mut err = *werr;
    let e = i64::from(z.x_position_start) + ol;
    if ol == 0 {
        return sc;
    }
    if i64::from(z.x_position_end) + 1 - i64::from(z.x_position_start) <= ol {
        return 0;
    }
    let windows = &z.window_cigar.windows;
    let get_err = |i: i64| -> i64 {
        if windows[i as usize].y_end != -1 {
            i64::from(windows[i as usize].error)
        } else {
            0
        }
    };
    if k == wn {
        k -= 1;
        err -= get_err(k);
    }
    if i64::from(windows[k as usize].x_start) >= e {
        err += get_err(k);
        while k >= 0 && i64::from(windows[k as usize].x_start) >= e {
            err -= get_err(k);
            k -= 1;
        }
    } else {
        while k < wn && i64::from(windows[k as usize].x_start) < e {
            err += get_err(k);
            k += 1;
        }
        k -= 1;
    }
    if k < 0 {
        k = 0;
        err = 0;
        *wi = k;
        *werr = err;
    } else {
        err -= get_err(k);
        *wi = k;
        *werr = err;
        if i64::from(windows[k as usize].x_end) + 1 <= e {
            err += get_err(k);
        } else {
            err += (((e - i64::from(windows[k as usize].x_start)) as f64
                / (i64::from(windows[k as usize].x_end) + 1 - i64::from(windows[k as usize].x_start)) as f64)
                * get_err(k) as f64) as i64;
        }
    }
    let mut sc = sc - (ol - ((err as f64) * err_sc_r) as i64);
    if sc <= 0 {
        sc = 1;
    }
    sc
}

/// `ed_dp_c` (`inter.cpp` 6762).
#[allow(clippy::too_many_arguments)]
fn ed_dp_c(
    o: &[OverlapRegion],
    res: &mut Vec<UlOv>,
    ex: &mut Vec<UlOv>,
    context: &UlRescallContext<'_>,
    bw: i64,
    diff_ec_ul: f64,
    qlen: i64,
    max_skip: i64,
    srt: &mut Vec<u64>,
    idx: &mut Vec<u64>,
    track: &mut Vec<u64>,
    err_sc: f64,
    mode: u64,
    need_srt: bool,
) -> i64 {
    if res.is_empty() {
        return 0;
    }
    let res_n = res.len();
    let g = &context.unitig_graph.graph;
    if need_srt {
        radix_sort_by_key_32(res, &|entry: &UlOv| entry.qe);
        let mut j = 0usize;
        for i in 1..=res_n {
            if i == res_n || res[i].qe != res[j].qe {
                if i - j > 1 {
                    radix_sort_by_key_32(&mut res[j..i], &|entry: &UlOv| entry.qs);
                }
                j = i;
            }
        }
    }
    if res[0].qe == qlen as u32 && res[res_n - 1].qs == 0 {
        for i in 0..res_n {
            let li = res[i];
            assert!(li.qs == 0 && li.qe == qlen as u32);
            let mut csc = i64::from(li.qe - li.qs);
            let minus_sc = (o[li.qn as usize].non_homopolymer_errors as f64 * err_sc) as i64;
            csc -= minus_sc;
            if csc <= 0 {
                csc = 1;
            }
            let mut mm_sc = csc;
            let mm_idx = -1i64;
            if mm_sc > 0x7fffffff {
                mm_sc = 0x7fffffff;
            }
            track[i] = push_sc_pre(mm_sc, mm_idx);
            srt[i] = (track[i] >> 32) << 32 | i as u64;
        }
    } else {
        idx.clear();
        idx.resize(res_n, 0);
        for i in 0..res_n {
            let li = res[i];
            let li_v = (li.tn << 1) | u32::from(li.rev);
            let mm_ovlp = if mode != 0 {
                max_ovlp_src(context, li_v ^ 1)
            } else {
                max_ovlp(g, li_v ^ 1)
            };
            let mut x = ((i64::from(li.qs) + mm_ovlp) as f64 * diff_ec_ul) as i64;
            if x < bw {
                x = bw;
            }
            x += i64::from(li.qs) + mm_ovlp;
            if x > qlen + 1 {
                x = qlen + 1;
            }
            let x = find_ul_ov_max(&res[..i], (x + i64::from(crate::assembly::ul_rescall::G_CHAIN_INDEL)) as u32);
            let mut csc = i64::from(li.qe - li.qs);
            let minus_sc = (o[li.qn as usize].non_homopolymer_errors as f64 * err_sc) as i64;
            csc -= minus_sc;
            if csc <= 0 {
                csc = 1;
            }
            let mut mm_sc = csc;
            let mut mm_idx = -1i64;
            let mut n_skip = 0i64;
            let mut wi = 0i64;
            let mut werr = 0i64;
            let mut j = x;
            while j >= 0 {
                let lj = res[j as usize];
                let lj_v = (lj.tn << 1) | u32::from(lj.rev);
                if i64::from(lj.qe) + i64::from(crate::assembly::ul_rescall::G_CHAIN_INDEL) <= i64::from(li.qs) {
                    break;
                }
                if lj.qs > li.qs + crate::assembly::ul_rescall::G_CHAIN_INDEL {
                    j -= 1;
                    continue;
                }
                let qo = infer_rovlp_unitigs(&li, &lj, &context.unitig_graph.unitigs);
                if li_v != lj_v
                    && get_ecov_adv_full(context, li_v ^ 1, lj_v ^ 1, bw, diff_ec_ul, qo, mode, None)
                {
                    let qovl = if li.qe.min(lj.qe) > li.qs.max(lj.qs) {
                        i64::from(li.qe.min(lj.qe) - li.qs.max(lj.qs))
                    } else {
                        0
                    };
                    let sc = comput_sc_partial_cigar(csc, qovl, err_sc, &o[li.qn as usize], &mut wi, &mut werr)
                        + pop_sc(track[j as usize]);
                    if sc > mm_sc {
                        mm_sc = sc;
                        mm_idx = j;
                        if n_skip > 0 {
                            n_skip -= 1;
                        }
                    } else if idx[j as usize] == i as u64 {
                        n_skip += 1;
                        if n_skip > max_skip {
                            break;
                        }
                    }
                    let pj = pop_pre(track[j as usize]);
                    if pj >= 0 {
                        idx[pj as usize] = i as u64;
                    }
                }
                j -= 1;
            }
            if mm_sc > 0x7fffffff {
                mm_sc = 0x7fffffff;
            }
            track[i] = push_sc_pre(mm_sc, mm_idx);
            srt[i] = (track[i] >> 32) << 32 | i as u64;
        }
    }
    let mut n_v = 0usize;
    let mut n_u = 0usize;
    radix_sort_u64(&mut srt[..res_n]);
    for k in (0..res_n).rev() {
        let n_v0 = n_v;
        let mut i = srt[k] as u32 as i64;
        while i >= 0 && (track[i as usize] & 0x8000_0000) == 0 {
            ex.push(res[i as usize]);
            n_v += 1;
            track[i as usize] |= 0x8000_0000;
            i = pop_pre(track[i as usize]);
        }
        if n_v0 == n_v {
            continue;
        }
        let sc = if i < 0 {
            pop_sc(srt[k])
        } else {
            pop_sc(srt[k]) - pop_sc(track[i as usize])
        };
        if sc < 0 {
            n_v = n_v0;
            continue;
        }
        idx[n_u] = ((sc as u64) << 32) | (n_v - n_v0) as u64;
        n_u += 1;
    }
    n_v = 0;
    let mut n_v0;
    for k in 0..n_u {
        n_v0 = n_v;
        n_v += idx[k] as u32 as usize;
        res[k].qn = (idx[k] >> 32) as u32;
        res[k].ts = n_v0 as u32;
        res[k].te = n_v as u32;
        let rev_n = (idx[k] as u32) >> 1;
        res[k].qs = u32::MAX;
        res[k].qe = ex[n_v0].qe;
        for i in 0..rev_n as usize {
            ex.swap(n_v0 + i, n_v - i - 1);
            if res[k].qs > ex[n_v0 + i].qs {
                res[k].qs = ex[n_v0 + i].qs;
            }
            if res[k].qs > ex[n_v - i - 1].qs {
                res[k].qs = ex[n_v - i - 1].qs;
            }
        }
        if idx[k] as u32 & 1 != 0 {
            let i = rev_n as usize;
            if res[k].qs > ex[n_v0 + i].qs {
                res[k].qs = ex[n_v0 + i].qs;
            }
        }
    }
    res.truncate(n_u);
    radix_sort_by_key_32(res, &|entry: &UlOv| entry.qn);
    n_v as i64
}

/// `get_num_wins` (`Correct.cpp` 782).
pub fn get_num_wins(s: i64, e: i64, block_s: i64) -> i64 {
    let nl = e - (s / block_s) * block_s;
    let mut nw = nl / block_s;
    if nl % block_s > 0 {
        nw += 1;
    }
    nw
}

/// `set_sec_e_min` (`inter.cpp` 6906): the minimum window error per block.
fn set_sec_e_min(z: &OverlapRegion, w_idx: &mut [u64], wl: i64, ql: i64) {
    let windows = &z.window_cigar.windows;
    for w in windows.iter() {
        let wid = (i64::from(w.x_start) / wl) as usize;
        let ws = wid as i64 * wl;
        let mut we = ws + wl;
        if we > ql {
            we = ql;
        }
        we -= 1;
        if ws == i64::from(w.x_start) && we == i64::from(w.x_end) && w.y_end != -1 {
            if w_idx[wid] == u64::MAX || w_idx[wid] > w.error as u64 {
                w_idx[wid] = w.error as u64;
            }
        }
    }
}

/// `cal_sec_e_min` (`inter.cpp` 6922).
fn cal_sec_e_min(z: &OverlapRegion, w_idx: &[u64], wl: i64, ql: i64, sec_err: Option<&mut i64>) -> i64 {
    let windows = &z.window_cigar.windows;
    let mut o = [0i64; 2];
    let mut tot_e = 0i64;
    for w in windows.iter() {
        let wid = (i64::from(w.x_start) / wl) as usize;
        if w_idx[wid] == u64::MAX {
            continue;
        }
        let ws = wid as i64 * wl;
        let mut we = ws + wl;
        if we > ql {
            we = ql;
        }
        we -= 1;
        let self_err: u64 = if w.y_end != -1 {
            w.error as u64
        } else {
            ((THRESHOLD_MAXIMUM + 1) as i64) as u64
        };
        if ws == i64::from(w.x_start) && we == i64::from(w.x_end) {
            if self_err <= w_idx[wid] {
                o[0] += we + 1 - ws;
            } else {
                o[1] += we + 1 - ws;
                tot_e += (self_err - w_idx[wid]) as i64;
            }
        }
    }
    if let Some(out) = sec_err {
        *out = tot_e;
    }
    o[0] - ((o[1] as f64 * ERROR_RATE) as i64) - tot_e * 5
}

/// `set_w_e` (`inter.cpp` 6949): the maximum window error per block.
fn set_w_e(z: &OverlapRegion, w_idx: &mut [u64], wl: i64, ql: i64) {
    let windows = &z.window_cigar.windows;
    for w in windows.iter() {
        let wid = (i64::from(w.x_start) / wl) as usize;
        let ws = wid as i64 * wl;
        let mut we = ws + wl;
        if we > ql {
            we = ql;
        }
        we -= 1;
        if ws == i64::from(w.x_start) && we == i64::from(w.x_end) && w.y_end != -1 {
            if w_idx[wid] == u64::MAX || w_idx[wid] < w.error as u64 {
                w_idx[wid] = w.error as u64;
            }
        }
    }
}

/// `ck_w_err` (`inter.cpp` 6962).
fn ck_w_err(z: &mut OverlapRegion, w_idx: &[u64], wl: i64, ql: i64) -> u32 {
    let windows = &z.window_cigar.windows;
    let mut ol = 0i64;
    let mut e = [0i64; 2];
    for w in windows.iter() {
        let wid = (i64::from(w.x_start) / wl) as usize;
        if w_idx[wid] == u64::MAX {
            continue;
        }
        let ws = wid as i64 * wl;
        let mut we = ws + wl;
        if we > ql {
            we = ql;
        }
        we -= 1;
        if ws == i64::from(w.x_start) && we == i64::from(w.x_end) {
            ol += we + 1 - ws;
            e[0] += w_idx[wid] as i64;
            if w.y_end != -1 {
                e[1] += i64::from(w.error);
            } else {
                e[1] += THRESHOLD_MAXIMUM + 1;
            }
        }
    }
    if e[1] > e[0] + 64 && e[1] as f64 > e[0] as f64 + (ol as f64) * 0.01 {
        z.non_homopolymer_errors = (e[1] - e[0]) as u32;
        return 0;
    }
    1
}

/// `filter_sec` (`inter.cpp` 6983): keeps only the overlaps that belong to
/// the best chains or beat their block's minimum error; compacts `ol`.
fn filter_sec(
    ol: &mut Vec<OverlapRegion>,
    idx: &[UlOv],
    a: &[UlOv],
    w_idx: &mut [u64],
    nw: usize,
    wl: u64,
    ql: u64,
) -> bool {
    if idx.is_empty() {
        return true;
    }
    let on = ol.len();
    for entry in w_idx.iter_mut().take(nw) {
        *entry = u64::MAX;
    }
    for region in ol.iter_mut().take(on) {
        region.is_match = 0;
    }
    let mut on_contain = 0i64;
    for entry in idx {
        for z in entry.ts as usize..entry.te as usize {
            ol[a[z].qn as usize].is_match = 2;
            let qn = a[z].qn as usize;
            let (region_a, region_b) = (&ol[qn], &ol[qn]);
            let _ = region_b;
            set_w_e(region_a, w_idx, wl as i64, ql as i64);
        }
        on_contain += i64::from((entry.te - entry.ts) == 1);
    }
    if on_contain == idx.len() as i64 {
        on_contain = 0;
    } else {
        on_contain = -(on as i64) - 1;
    }
    let mut max_i = a[idx[idx.len() - 1].ts as usize].qn as i64;
    let mut max_k = -1i64;
    let mut z = 0usize;
    for k in 0..on {
        if ol[k].is_match == 0 {
            ol[k].is_match = ck_w_err(&mut ol[k], w_idx, wl as i64, ql as i64) as u8;
        }
        if ol[k].is_match == 0 {
            continue;
        }
        if z != k {
            ol.swap(k, z);
        }
        if ol[z].x_position_start == 0 && ol[z].x_position_end == ql as u32 - 1 {
            on_contain += 1;
            if max_i == k as i64 {
                max_k = z as i64;
            }
        }
        ol[z].is_match = 1;
        z += 1;
    }
    if on_contain == z as i64 {
        let k = max_k;
        if z != 0 && k >= 0 {
            ol.swap(k as usize, 0);
        }
        z = 1;
    }
    ol.truncate(z);
    z == 1
}

/// `gl_chain_flter` (`inter.cpp` 7054): the fast pre-filter of the rescall
/// worker; `need_phase` reports whether the phase pass must run.
#[allow(clippy::too_many_arguments)]
pub fn gl_chain_flter(
    context: &UlRescallContext<'_>,
    olist: &mut Vec<OverlapRegion>,
    dumy_srt: &mut Vec<u64>,
    ll_srt: &mut Vec<u64>,
    sps: &mut Vec<u64>,
    ll_tk: &mut Vec<UlOv>,
    ll_lo: &mut Vec<UlOv>,
    diff_ec_ul: f64,
    wl: i64,
    ql: i64,
    need_phase: &mut u32,
) -> i64 {
    *need_phase = 1;
    ll_tk.clear();
    ll_lo.clear();
    let idx = ll_lo;
    introsort_overlaps_by_xe(olist);
    regen_ul_ov_t_lst(&context.unitig_graph.unitigs, olist, idx);
    if idx.is_empty() {
        return 0;
    }
    ll_srt.reserve(idx.len());
    sps.reserve(idx.len());
    ll_tk.reserve(idx.len());
    let mut ex = std::mem::take(ll_tk);
    let mut idx_scratch = std::mem::take(ll_srt);
    let mut track = std::mem::take(sps);
    let occ = ed_dp_c(
        olist,
        idx,
        &mut ex,
        context,
        crate::assembly::ul_rescall::G_CHAIN_BW,
        crate::assembly::ul_rescall::N_GCHAIN_RATE,
        ql,
        75,
        dumy_srt,
        &mut idx_scratch,
        &mut track,
        ERROR_RATE,
        0,
        false,
    );
    *sps = track;
    *ll_srt = idx_scratch;
    *ll_tk = ex;
    if occ == 0 || idx.is_empty() {
        return 0;
    }
    let idx_n = idx.len();
    let p_last = idx[idx_n - 1];
    if idx_n <= 1 {
        if p_last.te - p_last.ts <= 1 {
            *need_phase = 0;
        }
        return 0;
    }
    if (p_last.qe - p_last.qs) as i64 <= (ql as f64 * 0.333333) as i64 {
        return 0;
    }
    let mut i = idx_n as i64 - 1;
    let mut occ2 = i64::from(p_last.te - p_last.ts);
    if ((p_last.qe - p_last.qs) as i64) < ql && idx_n > 1 {
        occ2 = 0;
        while i >= 0 {
            let p = idx[i as usize];
            let mut k = i as usize + 1;
            while k < idx.len() {
                let m = idx[k];
                let ovlp = if m.qe.min(p.qe) > m.qs.max(p.qs) {
                    i64::from(m.qe.min(p.qe) - m.qs.max(p.qs))
                } else {
                    0
                };
                if ((ovlp as f64) > ((m.qe - m.qs) as f64 * 0.005)
                    || (ovlp as f64) > ((p.qe - p.qs) as f64 * 0.015))
                    && ovlp > 32
                {
                    break;
                }
                if ovlp == i64::from(m.qe - m.qs) || ovlp == i64::from(p.qe - p.qs) {
                    break;
                }
                k += 1;
            }
            if k < idx.len() {
                break;
            }
            occ2 += i64::from(p.te - p.ts);
            i -= 1;
        }
        i += 1;
    }
    if occ2 == olist.len() as i64 {
        return 1;
    }
    let nw = get_num_wins(0, ql, wl) as usize;
    ll_srt.clear();
    ll_srt.resize(nw, u64::MAX);
    if filter_sec(
        olist,
        &idx[i as usize..],
        ll_tk,
        ll_srt,
        nw,
        wl as u64,
        ql as u64,
    ) {
        *need_phase = 0;
    }
    1
}

/// `ks_introsort_or_xe` (`inter.cpp` 23): orders overlaps by
/// `(x_pos_e << 32 | x_pos_s)`.
pub fn introsort_overlaps_by_xe(overlaps: &mut [OverlapRegion]) {
    introsort_klib_by(overlaps, &|a: &OverlapRegion, b: &OverlapRegion| {
        ((u64::from(a.x_position_end) << 32) | u64::from(a.x_position_start))
            < ((u64::from(b.x_position_end) << 32) | u64::from(b.x_position_start))
    });
}

/// `check_trans_rate` (`inter.cpp` 3981).
pub fn check_trans_rate(a: &[UlOv], trans_thres: f32) -> bool {
    let mut sp = u32::MAX;
    let mut ep = u32::MAX;
    let mut tts = u32::MAX;
    let mut tte = 0u32;
    let mut el = 0u32;
    for entry in a.iter().rev() {
        if entry.qs < tts {
            tts = entry.qs;
        }
        if entry.qe > tte {
            tte = entry.qe;
        }
        if !entry.el {
            continue;
        }
        if sp == u32::MAX || entry.qe <= sp {
            if sp != u32::MAX {
                el += ep - sp;
            }
            sp = entry.qs;
            ep = entry.qe;
        } else {
            sp = sp.min(entry.qs);
        }
    }
    if sp != u32::MAX {
        el += ep - sp;
    }
    let iel = (tte - tts) - el;
    iel == 0 || (iel as f32) <= ((tte - tts) as f32 * trans_thres)
}

/// `ff_chain` (`inter.cpp` 4004).
#[allow(clippy::too_many_arguments)]
pub fn ff_chain(
    context: &UlRescallContext<'_>,
    idx: &mut Vec<UlOv>,
    qlen: i64,
    cov_rate: f32,
    trans_thres: f32,
    a: &[UlOv],
    olist: Option<&mut Vec<OverlapRegion>>,
    hap: Option<&mut HaplotypeEvidenceAlloc>,
    diff_ec_ul: f64,
    win_len: i64,
) -> bool {
    if idx.is_empty() {
        return false;
    }
    let m = idx[idx.len() - 1];
    if i64::from(m.qe - m.qs) <= (qlen as f32 * cov_rate) as i64 {
        return false;
    }
    if trans_thres >= 0.0 {
        if check_trans_rate(&a[m.ts as usize..m.te as usize], trans_thres) {
            return true;
        }
        if let (Some(olist), Some(hap)) = (olist, hap) {
            let idx_n = idx.len();
            let mut resc_tk = 0i64;
            let mut f = false;
            let mut si = 0usize;
            for z in m.ts as usize..m.te as usize {
                if a[z].el {
                    idx.push(a[z]);
                } else {
                    let i = a[z].qn;
                    let het_n = update_ava_het_site(hap, i, Some(&mut si), None, true);
                    assert!(het_n > 0 && olist[i as usize].is_match == 2);
                    let u = &context.unitig_graph.unitigs[olist[i as usize].y_identifier as usize];
                    if u.reads.len() > 1 {
                        let he: Vec<HaplotypeEvidence> = hap.list[si..si + het_n].to_vec();
                        resc_tk += rescue_trans_ul_chains(
                            &context.unitig_graph.unitigs,
                            context.reads,
                            &olist[i as usize],
                            &he,
                            u,
                            idx,
                            diff_ec_ul,
                            win_len,
                            false,
                            None,
                        );
                    }
                }
            }
            if resc_tk > 0 {
                radix_sort_by_key_32(&mut idx[idx_n..], &|entry: &UlOv| entry.qe);
                f = check_trans_rate(&idx[idx_n..], trans_thres);
            }
            idx.truncate(idx_n);
            return f;
        }
        return false;
    }
    true
}

/// `fill_unaligned_alignments` (`inter.cpp` 4944): the middle links of a
/// chain get `qs = qe = -1` before `update_uovlp_chain_qse`; here the
/// reference only fills them when the chain is longer than two.
fn fill_unaligned_alignments(a: &mut [MgLchain], offset: i64) {
    if a.len() > 2 {
        let n = a.len() - 1;
        for entry in a.iter_mut().take(n).skip(1) {
            entry.qs = -1;
            entry.qe = -1;
        }
    }
    let _ = offset;
}

/// `update_ul_vec_t_ug` (`inter.cpp` 4974): rewrites the read's block list
/// with the winning chains.
pub fn update_ul_vec_t_ug(
    unitigs: &[MaUtg],
    g: &Asg,
    rch: &mut crate::assembly::ul_store::UlReadVec,
    uc: &mut Vec<MgLchain>,
) {
    let ucn = uc.len();
    let mut a_n = 0usize;
    let mut k = 0usize;
    while k < ucn {
        let cnt = uc[k].cnt as usize;
        fill_unaligned_alignments(&mut uc[k + 1..k + 1 + cnt], (k + 1) as i64);
        uc[k].hash_pre = u32::MAX;
        uc[k].off = -1;
        a_n += cnt;
        k += cnt + 1;
    }
    rch.bb.clear();
    rch.bb.reserve(a_n);
    let mut k = 0usize;
    while k < ucn {
        let cnt = uc[k].cnt as usize;
        uc[k].hash_pre = u32::MAX;
        uc[k].off = -1;
        let a = &uc[k + 1..k + 1 + cnt];
        for (m, entry) in a.iter().enumerate() {
            rch.bb.push(crate::assembly::ul_store::UcBlock {
                hid: entry.v >> 1,
                rev: entry.v & 1 != 0,
                pchain: 1,
                base: false,
                el: true,
                qs: entry.qs as u32,
                qe: entry.qe as u32,
                te: entry.re as u32,
                ts: entry.rs as u32,
                pidx: (k + 1 + m) as u32,
                pdis: u32::MAX,
                aidx: u32::MAX,
            });
        }
        k += cnt + 1;
    }
    let a_n = rch.bb.len();
    crate::container::ksort::radix_sort_by_key_32(&mut rch.bb, &|b: &crate::assembly::ul_store::UcBlock| b.qe);
    let mut lk = 0usize;
    let mut m: i64 = -1;
    let mut l: i64 = -1;
    for k in 1..=a_n {
        uc[rch.bb[k - 1].pidx as usize].off = (k - 1) as i32;
        if m < rch.bb[k - 1].pidx as i64 {
            m = rch.bb[k - 1].pidx as i64;
            l = k as i64 - 1;
        }
        if k == a_n || rch.bb[k].qe != rch.bb[lk].qe {
            if k - lk > 1 {
                crate::container::ksort::radix_sort_by_key_32(&mut rch.bb[lk..k], &|b: &crate::assembly::ul_store::UcBlock| b.qs);
            }
            lk = k;
        }
    }
    for k in 0..a_n {
        if uc[rch.bb[k].pidx as usize].hash_pre == u32::MAX {
            rch.bb[k].pidx = u32::MAX;
            rch.bb[k].pdis = u32::MAX;
            rch.bb[k].aidx = u32::MAX;
            continue;
        }
        m = rch.bb[k].pidx as i64;
        let hp = uc[m as usize].hash_pre as usize;
        rch.bb[k].pidx = uc[hp].off as u32;
        rch.bb[k].pdis = uc[hp].dist_pre as u32;
        let pidx = rch.bb[k].pidx as usize;
        rch.bb[pidx].aidx = k as u32;
    }
    let mut sp = u32::MAX;
    let mut ep = u32::MAX;
    let mut ch_n = 0u64;
    let mut k = l;
    let mut l_sum = 0i64;
    while k >= 0 {
        let entry = &rch.bb[k as usize];
        if sp == u32::MAX || entry.qe <= sp {
            if sp != u32::MAX {
                l_sum += i64::from(ep - sp);
            }
            sp = entry.qs;
            ep = entry.qe;
        } else {
            sp = sp.min(entry.qs);
        }
        if entry.pidx == u32::MAX {
            k = -1;
        } else {
            k = i64::from(entry.pidx);
        }
        ch_n += 1;
    }
    rch.dd = 0;
    if sp != u32::MAX {
        l_sum += i64::from(ep - sp);
    }
    let l_rest = i64::from(rch.rlen) - l_sum;
    if l_rest == 0 {
        rch.dd = 1;
    } else if (l_rest as f64) < (rch.rlen as f64) * 0.001 {
        rch.dd = 2;
    } else if ch_n < rch.bb.len() as u64 {
        rch.dd = 3;
    }
    let _ = unitigs;
    let _ = g;
}

/// `filter_topN` (`inter.cpp` 7223): keeps at most `max_cov` query bases of
/// chains, ranked by the per-block minimum error.
#[allow(clippy::too_many_arguments)]
pub fn filter_topN(
    ol: &mut Vec<OverlapRegion>,
    aln: &mut Vec<UlOv>,
    ql: u64,
    wl: u64,
    max_cov: u64,
    ll_srt: &mut Vec<u64>,
) {
    let mut cc = 0u64;
    for entry in aln.iter() {
        cc += u64::from(entry.qe - entry.qs);
    }
    let cc_max = ql * max_cov;
    if cc <= cc_max {
        return;
    }
    let nw = get_num_wins(0, ql as i64, wl as i64) as usize;
    ll_srt.clear();
    ll_srt.resize(nw + ol.len(), u64::MAX);
    let (w_idx, srt) = ll_srt.split_at_mut(nw);
    for entry in w_idx.iter_mut() {
        *entry = u64::MAX;
    }
    for region in ol.iter_mut() {
        region.is_match = 0;
        region.alignment_length = 0;
        crate::correction::recalculate::append_unmatched_windows(
            region.x_position_start as i32,
            region.x_position_end as i32,
            &mut region.window_cigar,
            wl as i64,
        );
        set_sec_e_min(region, w_idx, wl as i64, ql as i64);
    }
    let mut j = 0usize;
    for i in 1..=aln.len() {
        if i == aln.len() || aln[i].qn != aln[j].qn {
            let z = &mut ol[aln[j].qn as usize];
            z.alignment_length = 0;
            for k in j..i {
                z.alignment_length += aln[k].qe - aln[k].qs;
                assert!(k <= j || aln[k].qs >= aln[k - 1].qe);
            }
            j = i;
        }
    }
    let mut m = 0usize;
    let mut h0 = 0usize;
    let mut h1 = 0usize;
    for (i, region) in ol.iter().enumerate() {
        let sc = cal_sec_e_min(region, w_idx, wl as i64, ql as i64, None);
        if sc >= 0 {
            srt[m] = (sc as u64) << 32 | i as u64 | 0x8000_0000_0000_0000;
            h0 += 1;
        } else {
            srt[m] = ((-sc) as u64) << 32 | i as u64;
            h1 += 1;
        }
        m += 1;
    }
    radix_sort_u64(&mut srt[..ol.len()]);
    if h1 > 0 {
        let zn = h1 >> 1;
        for i in 0..zn {
            srt.swap(i, h1 - i - 1);
        }
    }
    let mut m = ol.len() as i64 - 1;
    let mut cc = 0u64;
    while m >= 0 && cc <= cc_max {
        cc += u64::from(ol[(srt[m as usize] as u32) as usize].alignment_length);
        ol[(srt[m as usize] as u32) as usize].is_match = 1;
        m -= 1;
    }
    let mut m = 0usize;
    for k in 0..ol.len() {
        ol[k].alignment_length = u32::MAX;
        if ol[k].is_match == 0 {
            continue;
        }
        ol[k].alignment_length = m as u32;
        m += 1;
    }
    let mut m = 0usize;
    for i in 0..aln.len() {
        if ol[aln[i].qn as usize].is_match == 0 {
            continue;
        }
        aln[m] = aln[i];
        aln[m].qn = ol[aln[m].qn as usize].alignment_length;
        m += 1;
    }
    aln.truncate(m);
    let mut m = 0usize;
    let mut sec_err = 0i64;
    for k in 0..ol.len() {
        if ol[k].is_match == 0 {
            let _ = cal_sec_e_min(&ol[k], w_idx, wl as i64, ql as i64, Some(&mut sec_err));
            ol[k].non_homopolymer_errors = sec_err as u32;
            continue;
        }
        if m != k {
            ol.swap(k, m);
        }
        m += 1;
    }
    ol.truncate(m);
}

/// `gdpchain_t` (`Correct.h`): the graph-chain scratch buffers of one worker.
#[derive(Default)]
pub struct GdpChain {
    /// `gdp->l`.
    pub l: Vec<MgLchain>,
    /// `gdp->swap`: the winning chains.
    pub swap: Vec<MgLchain>,
    /// `gdp->dst`.
    pub dst: Vec<MgPathDst>,
    /// `gdp->out` (`vec_sp_node_t` in the reference; folded into
    /// `hc_shortest_k`, so this carries the `dst_group` array instead).
    pub out: Vec<usize>,
    /// `gdp->path`.
    pub path: Vec<MgPathv>,
    /// `gdp->dst_done`.
    pub dst_done: Vec<u64>,
}

/// `gl_chain` (`inter.cpp` 6605): the rescall worker's chain driver.
#[allow(clippy::too_many_arguments)]
pub fn gl_chain(
    context: &UlRescallContext<'_>,
    rch: &mut crate::assembly::ul_store::UlReadVec,
    olist: &mut Vec<OverlapRegion>,
    dp: &mut ChainData,
    hap: &mut HaplotypeEvidenceAlloc,
    sps: &mut Vec<u64>,
    ll_tk: &mut Vec<UlOv>,
    ll_lo: &mut Vec<UlOv>,
    ll_srt: &mut Vec<u64>,
    gdp: &mut GdpChain,
    diff_ec_ul: f64,
    win_len: i64,
    qlen: i64,
) -> i64 {
    ll_tk.clear();
    ll_lo.clear();
    let idx = ll_lo;
    introsort_overlaps_by_xe(olist);
    gen_gl_aln(olist, &context.unitig_graph.unitigs, idx);
    if idx.is_empty() {
        return 0;
    }
    let mut occ = 0i64;
    let mut f = false;
    ll_tk.reserve(idx.len());
    let mut ex = std::mem::take(ll_tk);
    occ = gl_chain_lin(
        idx,
        olist,
        &mut ex,
        context,
        crate::assembly::ul_rescall::G_CHAIN_BW,
        crate::assembly::ul_rescall::N_GCHAIN_RATE,
        qlen,
        crate::assembly::ul::UG_SKIP_N,
        crate::assembly::ul::UG_ITER_N,
        crate::assembly::ul::UG_DIS_N,
        dp,
        UG_TRANS_W,
        0,
        false,
    );
    *ll_tk = ex;
    if occ != 0 {
        if ff_chain(
            context,
            idx,
            qlen,
            0.99,
            -1.0,
            ll_tk,
            None,
            None,
            diff_ec_ul,
            win_len,
        ) {
            let last = idx.last().unwrap();
            f = l2g_res_chain_sc(
                &context.unitig_graph.unitigs,
                &context.unitig_graph.graph,
                &ll_tk[last.ts as usize..last.te as usize],
                &mut gdp.swap,
            ) != 0;
        }
    }
    if !f {
        gen_gg_aln(
            olist,
            &context.unitig_graph.unitigs,
            UG_TRANS_W,
            &mut gdp.l,
        );
        ll_tk.clear();
        let max_idx = gl_chain_graph_adv(
            olist,
            context,
            &mut gdp.l,
            &mut gdp.swap,
            &mut gdp.dst,
            &mut gdp.out,
            &mut gdp.path,
            rch.rlen as i64,
            crate::assembly::ul_rescall::G_CHAIN_BW,
            diff_ec_ul,
            sps,
            dp,
            crate::assembly::ul::UG_SKIP_GRAPH_N,
            crate::assembly::ul::UG_ITER_N,
            crate::assembly::ul::UG_DIS_N,
            0,
            UG_TRANS_W,
            UG_TRANS_ERR_W != 0,
        );
        if max_idx >= 0
            && crate::assembly::ul_gchain::select_max_gchain(
                &context.unitig_graph.graph,
                sps,
                &gdp.l,
                ll_tk,
                &mut gdp.dst_done,
                &mut gdp.out,
                &mut gdp.path,
                &mut gdp.swap,
            ) != 0
        {
            f = true;
        }
    }
    if f {
        update_ul_vec_t_ug(
            &context.unitig_graph.unitigs,
            &context.unitig_graph.graph,
            rch,
            &mut gdp.swap,
        );
    }
    1
}

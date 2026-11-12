//! The graph-binning mask engine (`inter.cpp` 9500-12300): aligns every
//! unitig against the unitig minimizer index and records, per unitig pair,
//! the masked end pieces (`emask_t`) that the rescall phase pass later
//! honours. `graph_ovlp_binning` (`inter.cpp` 22860) builds the
//! `idx_emask_t` table; `push_graph_bin` fills one unitig's entry.

use crate::assembly::asg::{Asg, MaUtg};
use crate::assembly::ul::UlOv;
use crate::assembly::ul_phase::{cal_s_interval, cal_x_interval, Emask, MaskUlOv};
pub use crate::assembly::ul_phase::Emask as BinningEmask;
use crate::correction::ul_align::{gen_err_unaligned, WindowExtz};
use crate::correction::window::WindowList;
use crate::container::ksort::{radix_sort_by_key_32, radix_sort_u64};
use crate::correction::region::THRESHOLD_MAXIMUM;
use crate::overlap::region::OverlapRegion;

/// `GBIN_K` (`inter.cpp` 61).
pub const GBIN_K: i32 = 31;
/// `GBIN_W` (`inter.cpp` 62).
pub const GBIN_W: i32 = 31;
/// `GBIN_E` (`inter.cpp` 63).
pub const GBIN_E: f64 = 0.04;
/// `GBIN_BE` (`inter.cpp` 64).
pub const GBIN_BE: f64 = 0.005;
/// `GBIN_ME` (`inter.cpp` 65).
pub const GBIN_ME: u64 = 4;
/// `GBIN_L` (`inter.cpp` 68).
pub const GBIN_L: u64 = 256;
/// `REA_ALIGN_CUTOFF` (`inter.cpp` 58).
pub const REA_ALIGN_CUTOFF: i32 = 32;

/// `is_ualn_win` (`Correct.h` 1362).
#[inline]
fn is_ualn_win(window: &WindowList) -> bool {
    window.error == i16::MAX && window.cigar_length == 0 && window.extra_end < 0
}

/// `is_est_aln` (`Correct.h` 1364).
#[inline]
fn is_est_aln(window: &WindowList) -> bool {
    window.error < i16::MAX && window.cigar_length == 0
}

/// `gen_src_shared_interval_adv` (`inter.cpp` 11417): the graph-arc intervals
/// of `src`'s two- and three-hop neighbourhood, merged per target.
fn gen_src_shared_interval_adv(src: u32, graph: &Asg, res: &mut Vec<UlOv>) {
    let rn = res.len();
    let mut s1 = UlOv::default();
    let mut s2 = UlOv::default();
    let mut s3 = UlOv::default();
    let mut s4 = UlOv::default();
    let mut s5 = UlOv::default();

    for half in 0..2u32 {
        let st = (src << 1) ^ half;
        for arc in graph.arcs_of(st) {
            if arc.deleted {
                continue;
            }
            if (arc.v >> 1) == (st >> 1) {
                continue;
            }
            let w = arc.v;
            cal_s_interval(graph, arc, &mut s1);
            if s1.qn == src && s1.tn != src {
                res.push(s1);
            }
            for aw in graph.arcs_of(w ^ 1) {
                if aw.deleted {
                    continue;
                }
                if (aw.v >> 1) == (st >> 1) {
                    continue;
                }
                cal_s_interval(graph, aw, &mut s2);
                if !cal_x_interval(graph, &s1, &s2, &mut s3) {
                    continue;
                }
                if s3.qn == src && s3.tn != src {
                    res.push(s3);
                }
                for az in graph.arcs_of(aw.v ^ 1) {
                    if az.deleted {
                        continue;
                    }
                    cal_s_interval(graph, az, &mut s4);
                    if !cal_x_interval(graph, &s3, &s4, &mut s5) {
                        continue;
                    }
                    if s5.qn == src && s5.tn != src {
                        res.push(s5);
                    }
                }
            }
        }
    }

    radix_sort_by_key_32(&mut res[rn..], &|entry: &UlOv| entry.tn);
    if res.len() > rn {
        let mut i = rn;
        for k in rn + 1..=res.len() {
            if k == res.len() || res[k].tn != res[i].tn {
                if k - i > 1 {
                    for v in i..k {
                        if res[v].tn == u32::MAX {
                            continue;
                        }
                        for z in i..k {
                            if v == z {
                                continue;
                            }
                            if res[z].tn == u32::MAX {
                                continue;
                            }
                            if res[v].rev != res[z].rev {
                                continue;
                            }
                            let os = res[v].qs.max(res[z].qs);
                            let oe = res[v].qe.min(res[z].qe);
                            if oe <= os {
                                continue;
                            }
                            let ovlpq = u64::from(oe - os);
                            let os = res[v].ts.max(res[z].ts);
                            let oe = res[v].te.min(res[z].te);
                            if oe <= os {
                                continue;
                            }
                            let ovlpt = u64::from(oe - os);
                            let mut is_cov = false;
                            if (ovlpq == u64::from(res[v].qe - res[v].qs)
                                && ovlpt == u64::from(res[v].te - res[v].ts))
                                || (ovlpq == u64::from(res[z].qe - res[z].qs)
                                    && ovlpt == u64::from(res[z].te - res[z].ts))
                            {
                                is_cov = true;
                            }
                            if !is_cov {
                                if (ovlpq as f64) <= (u64::from(res[v].qe - res[v].qs) as f64) * 0.95 {
                                    continue;
                                }
                                if (ovlpq as f64) <= (u64::from(res[z].qe - res[z].qs) as f64) * 0.95 {
                                    continue;
                                }
                                if (ovlpt as f64) <= (u64::from(res[v].te - res[v].ts) as f64) * 0.95 {
                                    continue;
                                }
                                if (ovlpt as f64) <= (u64::from(res[z].te - res[z].ts) as f64) * 0.95 {
                                    continue;
                                }
                            }
                            if res[v].qs > res[z].qs {
                                res[v].qs = res[z].qs;
                            }
                            if res[v].qe < res[z].qe {
                                res[v].qe = res[z].qe;
                            }
                            if res[v].ts > res[z].ts {
                                res[v].ts = res[z].ts;
                            }
                            if res[v].te < res[z].te {
                                res[v].te = res[z].te;
                            }
                            res[z].tn = u32::MAX;
                        }
                    }
                    radix_sort_by_key_32(&mut res[i..k], &|entry: &UlOv| entry.qs);
                }
                i = k;
            }
        }
        let mut z = rn;
        for k in rn..res.len() {
            if res[k].tn == u32::MAX {
                continue;
            }
            res[z] = res[k];
            z += 1;
        }
        res.truncate(z);
    }
}

/// `trans_ovlp_connect` (`inter.cpp` 9537): whether the overlap is consistent
/// with the graph arc in either direction.
pub fn trans_ovlp_connect(p: &UlOv, graph: &Asg) -> bool {
    let check = |v: u32, w: u32| -> bool {
        let mut dg: i64 = -1;
        for arc in graph.arcs_of(v) {
            if arc.deleted || arc.v != w {
                continue;
            }
            dg = i64::from(arc.overlap);
            break;
        }
        if dg < 0 {
            return false;
        }
        for da in [i64::from(p.qe - p.qs), i64::from(p.te - p.ts)] {
            let dif = (dg - da).abs();
            let mut mm = (dg.max(da) as f64 * 0.06) as i64;
            if mm < 8 {
                mm = 8;
            }
            if dif <= mm {
                return true;
            }
        }
        false
    };
    let v0 = p.qn << 1;
    let w0 = (p.tn << 1) + u32::from(p.rev);
    if check(v0, w0) {
        return false;
    }
    if check(v0 ^ 1, w0 ^ 1) {
        return false;
    }
    true
}

/// `is_mask_err_full` (`inter.cpp` 11800): whether one masked interval fully
/// covers `[qs, qe) x [ts, te)`.
fn is_mask_err_full(sa: &[UlOv], qs: u64, qe: u64, ts: u64, te: u64) -> bool {
    for entry in sa {
        if u64::from(entry.qs) <= qs
            && u64::from(entry.qe) >= qe
            && u64::from(entry.ts) <= ts
            && u64::from(entry.te) >= te
        {
            return true;
        }
    }
    false
}

/// `cal_exact_len` (`inter.cpp` 11822): walks the window's cigar and reports
/// the longest masked-error prefix (forward) or suffix (backward).
#[allow(clippy::too_many_arguments)]
fn cal_exact_len(
    ez: &WindowExtz<'_>,
    rev: bool,
    sa: &[UlOv],
    ylen: i64,
    xoff: &mut i64,
    yoff: &mut i64,
    exac: &mut i64,
    err: &mut i64,
    backward: bool,
) -> bool {
    *xoff = 0;
    *yoff = 0;
    *exac = 0;
    *err = 0;
    if ez.cells.is_empty() {
        return false;
    }
    let cn = ez.cells.len() as i64;
    if !backward {
        let mut xk = i64::from(ez.text_start);
        let mut yk = i64::from(ez.pattern_start);
        for ck in 0..cn {
            let xs = xk;
            let ys = yk;
            let cell = ez.cells[ck as usize];
            let op = cell >> 14;
            let se = i64::from(cell & 0x3fff);
            if op != 2 {
                xk += se;
            }
            if op != 3 {
                yk += se;
            }
            let xe = xk;
            let ye = yk;
            if op == 0 {
                *exac += se;
                continue;
            }
            let ts = if rev { ylen - ye } else { ys };
            let te = if rev { ylen - ys } else { ye };
            if sa.is_empty() || !is_mask_err_full(sa, xs as u64, xe as u64, ts as u64, te as u64) {
                *xoff = xs - i64::from(ez.text_start);
                *yoff = ys - i64::from(ez.pattern_start);
                return false;
            }
            *exac = 0;
            *err += se;
        }
    } else {
        let mut xk = i64::from(ez.text_end) + 1;
        let mut yk = i64::from(ez.pattern_end) + 1;
        for ck in (0..cn).rev() {
            let xe = xk;
            let ye = yk;
            let cell = ez.cells[ck as usize];
            let op = cell >> 14;
            let se = i64::from(cell & 0x3fff);
            if op != 2 {
                xk -= se;
            }
            if op != 3 {
                yk -= se;
            }
            let xs = xk;
            let ys = yk;
            if op == 0 {
                *exac += se;
                continue;
            }
            let ts = if rev { ylen - ye } else { ys };
            let te = if rev { ylen - ys } else { ye };
            if sa.is_empty() || !is_mask_err_full(sa, xs as u64, xe as u64, ts as u64, te as u64) {
                *xoff = i64::from(ez.text_end) + 1 - xe;
                *yoff = i64::from(ez.pattern_end) + 1 - ye;
                return false;
            }
            *exac = 0;
            *err += se;
        }
    }
    *xoff = i64::from(ez.text_end) + 1 - i64::from(ez.text_start);
    *yoff = i64::from(ez.pattern_end) + 1 - i64::from(ez.pattern_start);
    true
}

/// `cal_exact_batch` (`inter.cpp` 11881): the unmasked error of overlap `z`
/// after the masked pieces `sa`, split into the left (`q0l/t0l/e0l`) and
/// right (`q1l/t1l/e1l`) exact runs.
#[allow(clippy::too_many_arguments)]
fn cal_exact_batch(
    z: &OverlapRegion,
    sa: &[UlOv],
    ylen: i64,
    q0l: &mut u64,
    t0l: &mut u64,
    e0l: &mut u64,
    q1l: &mut u64,
    t1l: &mut u64,
    e1l: &mut u64,
) -> i64 {
    let windows = &z.window_cigar.windows;
    let cells = &z.window_cigar.cells;
    let wn = windows.len() as i64;
    let rev = z.y_strand != 0;
    let tot = z.non_homopolymer_errors as i64;
    if tot == 0 {
        *q0l = u64::from(z.x_position_end + 1 - z.x_position_start);
        *q1l = *q0l;
        *t0l = u64::from(z.y_position_end + 1 - z.y_position_start);
        *t1l = *t0l;
        *e0l = u64::from(z.x_position_end + 1 - z.x_position_start);
        *e1l = *e0l;
        return 0;
    }
    *q0l = u64::MAX;
    *t0l = u64::MAX;
    *e0l = u64::MAX;
    *q1l = u64::MAX;
    *t1l = u64::MAX;
    *e1l = u64::MAX;
    if wn == 0 {
        return i32::MAX as i64;
    }
    *q0l = 0;
    *t0l = 0;
    *e0l = 0;
    let mut err0 = 0i64;
    let mut mask_win = false;
    let mut wk = 0i64;
    let mut xk = i64::from(windows[0].x_start);
    let mut yk = i64::from(windows[0].y_start);
    let mut is_t = false;
    if xk == i64::from(z.x_position_start) && yk == i64::from(z.y_position_start) {
        is_t = true;
    }
    if !is_t && !sa.is_empty() {
        let qs = i64::from(z.x_position_start);
        let qe = xk;
        let (ts, te) = if !rev {
            (i64::from(z.y_position_start), yk)
        } else {
            (ylen - yk, ylen - i64::from(z.y_position_start))
        };
        is_t = is_mask_err_full(sa, qs as u64, qe as u64, ts as u64, te as u64);
        if is_t {
            *q0l += (qe - qs) as u64;
            *t0l += (te - ts) as u64;
            *e0l = 0;
        }
    }
    if is_t {
        while wk < wn {
            let m = &windows[wk as usize];
            let qs = i64::from(m.x_start);
            let qe = i64::from(m.x_end) + 1;
            let (ts, te) = if !rev {
                (i64::from(m.y_start), i64::from(m.y_end) + 1)
            } else {
                (ylen - (i64::from(m.y_end) + 1), ylen - i64::from(m.y_start))
            };
            if is_ualn_win(m) || (is_est_aln(m) && m.error > 0) {
                if is_mask_err_full(sa, qs as u64, qe as u64, ts as u64, te as u64) {
                    *q0l += (qe - qs) as u64;
                    *t0l += (te - ts) as u64;
                    *e0l = 0;
                    mask_win = true;
                    wk += 1;
                    continue;
                }
                break;
            }
            if m.error == 0 {
                *q0l += (qe - qs) as u64;
                *t0l += (te - ts) as u64;
                *e0l += (qe - qs) as u64;
                wk += 1;
                continue;
            }
            let ez = WindowExtz::new(m, cells);
            let mut qoff = 0i64;
            let mut toff = 0i64;
            let mut elen = 0i64;
            let mut sube = 0i64;
            let is_f = cal_exact_len(&ez, rev, sa, ylen, &mut qoff, &mut toff, &mut elen, &mut sube, false);
            *q0l += qoff as u64;
            *t0l += toff as u64;
            err0 += sube;
            if elen == qoff && elen == toff {
                *e0l += elen as u64;
            } else {
                *e0l = elen as u64;
            }
            if !is_f {
                break;
            }
            wk += 1;
        }
    }
    *q1l = 0;
    *t1l = 0;
    *e1l = 0;
    let mut err1 = 0i64;
    let mut wk = wn - 1;
    let mut xk = i64::from(windows[wk as usize].x_end) + 1;
    let mut yk = i64::from(windows[wk as usize].y_end) + 1;
    let mut is_t = false;
    if xk == i64::from(z.x_position_end) + 1 && yk == i64::from(z.y_position_end) + 1 {
        is_t = true;
    }
    if !is_t && !sa.is_empty() {
        let qs = xk;
        let qe = i64::from(z.x_position_end) + 1;
        let (ts, te) = if !rev {
            (yk, i64::from(z.y_position_end) + 1)
        } else {
            (ylen - (i64::from(z.y_position_end) + 1), ylen - yk)
        };
        is_t = is_mask_err_full(sa, qs as u64, qe as u64, ts as u64, te as u64);
        if is_t {
            *q1l += (qe - qs) as u64;
            *t1l += (te - ts) as u64;
            *e1l = 0;
        }
    }
    if is_t {
        loop {
            let m = &windows[wk as usize];
            let qs = i64::from(m.x_start);
            let qe = i64::from(m.x_end) + 1;
            let (ts, te) = if !rev {
                (i64::from(m.y_start), i64::from(m.y_end) + 1)
            } else {
                (ylen - (i64::from(m.y_end) + 1), ylen - i64::from(m.y_start))
            };
            if is_ualn_win(m) || (is_est_aln(m) && m.error > 0) {
                if is_mask_err_full(sa, qs as u64, qe as u64, ts as u64, te as u64) {
                    *q1l += (qe - qs) as u64;
                    *t1l += (te - ts) as u64;
                    *e1l = 0;
                    mask_win = true;
                    if wk == 0 {
                        break;
                    }
                    wk -= 1;
                    continue;
                }
                break;
            }
            if m.error == 0 {
                *q1l += (qe - qs) as u64;
                *t1l += (te - ts) as u64;
                *e1l += (qe - qs) as u64;
                if wk == 0 {
                    break;
                }
                wk -= 1;
                continue;
            }
            let ez = WindowExtz::new(m, cells);
            let mut qoff = 0i64;
            let mut toff = 0i64;
            let mut elen = 0i64;
            let mut sube = 0i64;
            let is_f = cal_exact_len(&ez, rev, sa, ylen, &mut qoff, &mut toff, &mut elen, &mut sube, true);
            *q1l += qoff as u64;
            *t1l += toff as u64;
            err1 += sube;
            if elen == qoff && elen == toff {
                *e1l += elen as u64;
            } else {
                *e1l = elen as u64;
            }
            if !is_f {
                break;
            }
            if wk == 0 {
                break;
            }
            wk -= 1;
        }
    }
    if *q0l >= u64::from(z.x_position_end + 1 - z.x_position_start)
        && *t0l >= u64::from(z.y_position_end + 1 - z.y_position_start)
    {
        assert!(*q0l == *q1l && *t0l == *t1l && err0 == err1);
        assert!(mask_win || err0 == tot);
        return 0;
    }
    assert!(err0 + err1 < tot);
    tot - (err0 + err1)
}

/// `get_mask_interval` (`inter.cpp` 11766): the masked intervals of one
/// target; `ii` is the persistent cursor (`usize::MAX` means -1), `idx` is
/// permuted in place so consistent entries compact to the group head.
/// Returns `(start, n, n_skip)` into `idx`.
fn get_mask_interval(
    inside: &UlOv,
    idx: &mut [UlOv],
    ii: &mut usize,
    ql: i64,
    tl: i64,
) -> Option<(usize, usize, usize)> {
    let idx_n = idx.len() as i64;
    let mut n = 0usize;
    let mut n_skip = 0usize;
    let mut i = *ii as i64;
    if i >= idx_n {
        i = idx_n - 1;
    }
    if i < 0 {
        i = 0;
    }
    let mut found = None;
    if i >= 0 && i < idx_n {
        while i >= 0 && inside.tn <= idx[i as usize].tn {
            i -= 1;
        }
        if i < 0 && idx_n > 0 && inside.tn == idx[0].tn {
            i = 0;
        }
        if i >= 0 {
            while i < idx_n && inside.tn > idx[i as usize].tn {
                i += 1;
            }
            if i < idx_n && inside.tn == idx[i as usize].tn {
                let mut m = i;
                let start = i as usize;
                while i < idx_n && inside.tn == idx[i as usize].tn {
                    if !crate::assembly::ul_phase::check_ul_ov_t_consist(
                        inside,
                        &idx[i as usize],
                        ql,
                        tl,
                        0.06,
                    ) {
                        n_skip += 1;
                        i += 1;
                        continue;
                    }
                    if m != i {
                        idx.swap(i as usize, m as usize);
                    }
                    m += 1;
                    n += 1;
                    i += 1;
                }
                *ii = i as usize;
                if n > 0 {
                    found = Some((start, n, n_skip));
                }
            } else {
                *ii = i.max(0) as usize;
            }
        } else {
            *ii = (-1i64) as usize;
        }
    } else {
        *ii = i.max(0) as usize;
    }
    found
}

/// `full_lst_set` (`inter.cpp` 12014).
fn full_lst_set(a: &mut [Emask], sa: &[UlOv], skip_n: usize) {
    let mut ff;
    if sa.is_empty() && skip_n == 0 {
        ff = true;
    } else if skip_n > 0 {
        ff = false;
    } else {
        let mut k = 0usize;
        while k < sa.len() {
            let mut sql = i64::from(sa[k].qe - sa[k].qs);
            let mut stl = i64::from(sa[k].te - sa[k].ts);
            for entry in a.iter() {
                if sql <= 0 || stl <= 0 {
                    break;
                }
                let os = sa[k].qs.max(entry.qs);
                let oe = sa[k].qe.min(entry.qe);
                let ovlp = if oe > os { i64::from(oe - os) } else { 0 };
                sql -= ovlp;
                let os = sa[k].ts.max(entry.ts);
                let oe = sa[k].te.min(entry.te);
                let ovlp = if oe > os { i64::from(oe - os) } else { 0 };
                stl -= ovlp;
            }
            let q_len = i64::from(sa[k].qe - sa[k].qs);
            let t_len = i64::from(sa[k].te - sa[k].ts);
            if sql > 256 && (sql as f64) > q_len as f64 * 0.06 {
                break;
            }
            if sql <= 256 && (sql as f64) > q_len as f64 * 0.6 {
                break;
            }
            if stl > 256 && (stl as f64) > t_len as f64 * 0.06 {
                break;
            }
            if stl <= 256 && (stl as f64) > t_len as f64 * 0.6 {
                break;
            }
            k += 1;
        }
        ff = k >= sa.len();
    }
    for entry in a.iter_mut() {
        entry.full = ff;
    }
}

/// `push_emask0` (`inter.cpp` 12043).
#[allow(clippy::too_many_arguments)]
fn push_emask0(
    res: &mut Vec<Emask>,
    q0l: u64,
    t0l: u64,
    e0l: u64,
    q1l: u64,
    t1l: u64,
    e1l: u64,
    end: u64,
    is_pe: bool,
    minov: u64,
    p: &UlOv,
) {
    let mut m = Emask {
        tn: p.tn,
        rev: p.rev,
        rn: u32::MAX,
        ..Emask::default()
    };
    m.dir = end as u8;
    m.pe = is_pe;
    if end == 0 {
        assert!(q0l >= minov && t0l >= minov);
        m.el = e0l as u32;
        m.qs = p.qs;
        m.qe = p.qs + q0l as u32;
        if !m.rev {
            m.ts = p.ts;
            m.te = p.ts + t0l as u32;
        } else {
            m.te = p.te;
            m.ts = p.te - t0l as u32;
        }
    } else {
        assert!(q1l >= minov && t1l >= minov);
        m.el = e1l as u32;
        m.qe = p.qe;
        m.qs = p.qe - q1l as u32;
        if !m.rev {
            m.te = p.te;
            m.ts = p.te - t1l as u32;
        } else {
            m.ts = p.ts;
            m.te = p.ts + t1l as u32;
        }
    }
    res.push(m);
}

/// `overlap_region2ul_ov_t` (`inter.cpp` 11533).
fn overlap_region2ul_ov_t(region: &OverlapRegion, tl: u32) -> UlOv {
    let rev = region.y_strand != 0;
    let (ts, te) = if rev {
        (tl - (region.y_position_end + 1), tl - region.y_position_start)
    } else {
        (region.y_position_start, region.y_position_end + 1)
    };
    UlOv {
        qn: region.x_identifier,
        tn: region.y_identifier,
        el: true,
        sec: region.non_homopolymer_errors,
        rev,
        qs: region.x_position_start,
        qe: region.x_position_end + 1,
        ts,
        te,
    }
}

/// `push_emask_lst` (`inter.cpp` 12189).
#[allow(clippy::too_many_arguments)]
fn push_emask_lst(
    mk_srt: &mut Vec<UlOv>,
    mk_idx: &[u64],
    ol: &[OverlapRegion],
    minov: u64,
    unitigs: &[MaUtg],
    a: &mut [u64],
    res: &mut Vec<Emask>,
    is_pe: bool,
) {
    if is_pe {
        for entry in a.iter_mut() {
            *entry = (*entry << 32) >> 32;
        }
        radix_sort_u64(a);
        let mut l = 0usize;
        let mut si = 0usize;
        for i in 1..=a.len() {
            if i == a.len() || (a[i] >> 1) != (a[l] >> 1) {
                let z = &ol[(mk_idx[(a[l] >> 1) as usize] as u32) as usize];
                let p = overlap_region2ul_ov_t(z, unitigs[z.y_identifier as usize].length);
                assert!(i <= l + 2);
                assert!(u64::from(p.qe - p.qs) >= minov);
                assert!(u64::from(p.te - p.ts) >= minov);
                assert!(p.sec > 0);
                let ql = i64::from(unitigs[p.qn as usize].length);
                let tl = i64::from(unitigs[p.tn as usize].length);
                let (sa, _sn, skip_n) = match get_mask_interval(&p, mk_srt, &mut si, ql, tl) {
                    Some((start, n, skip)) => {
                        let slice = mk_srt[start..start + n].to_vec();
                        (slice, n, skip)
                    }
                    None => (Vec::new(), 0, 0),
                };
                let (q0l, t0l, e0l, q1l, t1l, e1l, re) =
                    cal_exact_batch_out(z, &sa, tl);
                assert!(re != 0);
                let cn = i - l;
                push_emask0(res, q0l, t0l, e0l, q1l, t1l, e1l, a[l] & 1, i == l + 2, minov, &p);
                l += 1;
                if l < i {
                    push_emask0(res, q0l, t0l, e0l, q1l, t1l, e1l, a[l] & 1, i == l + 2, minov, &p);
                }
                let begin = res.len() - cn;
                full_lst_set(&mut res[begin..], &sa, skip_n);
                l = i;
            }
        }
    } else {
        for entry in a.iter_mut() {
            *entry = (*entry << 32) >> 32;
            let z = &ol[(mk_idx[(*entry >> 1) as usize] as u32) as usize];
            let p = overlap_region2ul_ov_t(z, unitigs[z.y_identifier as usize].length);
            assert!(u64::from(p.qe - p.qs) >= minov);
            assert!(u64::from(p.te - p.ts) >= minov);
            if p.sec > 0 {
                let ql = i64::from(unitigs[p.qn as usize].length);
                let tl = i64::from(unitigs[p.tn as usize].length);
                let mut si = 0usize;
                let (sa, _sn, _skip) = match get_mask_interval(&p, mk_srt, &mut si, ql, tl) {
                    Some((start, n, skip)) => {
                        let slice = mk_srt[start..start + n].to_vec();
                        (slice, n, skip)
                    }
                    None => (Vec::new(), 0, 0),
                };
                let (_q0l, _t0l, _e0l, _q1l, _t1l, _e1l, re) =
                    cal_exact_batch_out(z, &sa, tl);
                assert!(re == 0);
            }
            let m = Emask {
                dir: 2,
                pe: false,
                tn: p.tn,
                rev: p.rev,
                full: true,
                el: u32::MAX,
                rn: u32::MAX,
                qs: p.qs,
                qe: p.qe,
                ts: p.ts,
                te: p.te,
                ..Emask::default()
            };
            res.push(m);
        }
    }
}

/// `cal_exact_batch` returning all six outputs plus the residual error, to
/// mirror the C six-out-param call.
fn cal_exact_batch_out(
    z: &OverlapRegion,
    sa: &[UlOv],
    tl: i64,
) -> (u64, u64, u64, u64, u64, u64, i64) {
    let mut q0l = 0u64;
    let mut t0l = 0u64;
    let mut e0l = 0u64;
    let mut q1l = 0u64;
    let mut t1l = 0u64;
    let mut e1l = 0u64;
    let re = cal_exact_batch(z, sa, tl, &mut q0l, &mut t0l, &mut e0l, &mut q1l, &mut t1l, &mut e1l);
    (q0l, t0l, e0l, q1l, t1l, e1l, re)
}

/// `srt_kv_emask_t` (`inter.cpp` 12239): pairs the `pe` entries and rebuilds
/// both `res` and `mk.srt` sorted by target.
fn srt_kv_emask_t(res: &mut Vec<Emask>, mk: &mut MaskUlOv) {
    mk.srt.clear();
    mk.idx.clear();
    mk.idx.reserve(res.len());
    for k in 0..res.len() {
        let m = k + 1;
        let mut z = (u64::from(res[k].tn) << 32) | (k << 1) as u64;
        if m < res.len()
            && res[k].tn == res[m].tn
            && res[k].rev == res[m].rev
            && res[k].pe
            && res[m].pe
            && res[k].dir == 0
            && res[m].dir == 1
        {
            z |= 1;
        }
        mk.idx.push(z);
    }
    radix_sort_u64(&mut mk.idx);
    mk.srt.reserve(res.len());
    for z in 0..mk.idx.len() {
        let k = ((mk.idx[z] as u32) >> 1) as usize;
        let m = ((res[k].dir as u32) << 1) | u32::from(res[k].pe);
        mk.srt.push(UlOv {
            tn: res[k].tn,
            qn: res[k].el,
            qs: res[k].qs,
            qe: res[k].qe,
            ts: res[k].ts,
            te: res[k].te,
            rev: res[k].rev,
            el: res[k].full,
            sec: m,
        });
        if mk.idx[z] & 1 != 0 {
            let k = k + 1;
            let m = ((res[k].dir as u32) << 1) | u32::from(res[k].pe);
            mk.srt.push(UlOv {
                tn: res[k].tn,
                qn: res[k].el,
                qs: res[k].qs,
                qe: res[k].qe,
                ts: res[k].ts,
                te: res[k].te,
                rev: res[k].rev,
                el: res[k].full,
                sec: m,
            });
        }
    }
    for (k, p) in mk.srt.iter().enumerate() {
        res[k].tn = p.tn;
        res[k].el = p.qn;
        res[k].qs = p.qs;
        res[k].qe = p.qe;
        res[k].ts = p.ts;
        res[k].te = p.te;
        res[k].rev = p.rev;
        res[k].full = p.el;
        res[k].rn = u32::MAX;
        res[k].pe = p.sec & 1 != 0;
        res[k].dir = (p.sec >> 1) as u8;
    }
}

/// `push_graph_bin` (`inter.cpp` 12297): fills unitig `rid`'s emask entry
/// from its unitig-vs-unitig overlaps.
#[allow(clippy::too_many_arguments)]
pub fn push_graph_bin(
    ol: &mut [OverlapRegion],
    unitigs: &[MaUtg],
    graph: &Asg,
    erate: f64,
    maxe: u64,
    minov: u64,
    rid: u32,
    b0: &mut Vec<u64>,
    mk: &mut MaskUlOv,
    res: &mut Vec<Emask>,
) {
    b0.clear();
    let tl = u64::from(unitigs[rid as usize].length);
    let mocc = (tl / 8).max(16);
    mk.srt.clear();
    mk.idx.clear();
    gen_src_shared_interval_adv(rid, graph, &mut mk.srt);
    for (k, region) in ol.iter().enumerate() {
        mk.idx.push((u64::from(region.y_identifier) << 32) | k as u64);
    }
    radix_sort_u64(&mut mk.idx);
    let mut occ = [0u64; 3];
    let mut si = 0usize;
    for i in 0..mk.idx.len() {
        let k = (mk.idx[i] as u32) as usize;
        let z = &ol[k];
        let p = overlap_region2ul_ov_t(z, unitigs[z.y_identifier as usize].length);
        if u64::from(p.qe - p.qs) >= minov && u64::from(p.te - p.ts) >= minov {
            if !trans_ovlp_connect(&p, graph) {
                continue;
            }
            if p.sec == 0 {
                let zt = (u64::from(u32::MAX) << 32) | ((i << 1) as u64);
                occ[2] += 1;
                b0.push(zt);
                continue;
            }
            let (sa, _sn, _skip_n) = {
                match get_mask_interval(
                    &p,
                    &mut mk.srt,
                    &mut si,
                    i64::from(unitigs[p.qn as usize].length),
                    i64::from(unitigs[p.tn as usize].length),
                ) {
                    Some((start, n, skip)) => {
                        let slice = mk.srt[start..start + n].to_vec();
                        (slice, n, skip)
                    }
                    None => (Vec::new(), 0, 0),
                }
            };
            let (q0l, t0l, _e0l, q1l, t1l, _e1l, re) =
                cal_exact_batch_out(z, &sa, i64::from(unitigs[p.tn as usize].length));
            if re == 0 {
                let zt = (u64::from(u32::MAX) << 32) | ((i << 1) as u64);
                occ[2] += 1;
                b0.push(zt);
                continue;
            }
            let mm: u64 = if (re as u64) <= maxe {
                0x8000_0000_0000_0000
            } else {
                0
            };
            if q0l >= minov && t0l >= minov {
                let zt = (q0l.min(t0l) << 32) | (((i << 1) as u64)) | mm;
                occ[usize::from(mm != 0)] += 1;
                b0.push(zt);
            }
            if q1l >= minov && t1l >= minov {
                let zt = (q1l.min(t1l) << 32) | (((i << 1) as u64) + 1) | mm;
                occ[usize::from(mm != 0)] += 1;
                b0.push(zt);
            }
        }
    }
    radix_sort_u64(b0);
    let a2n = occ[0] as usize;
    let a1n = occ[1] as usize;
    let a0n = occ[2] as usize;
    let mut a2: Vec<u64> = b0[..a2n].to_vec();
    let mut a1: Vec<u64> = b0[a2n..a2n + a1n].to_vec();
    let mut a0: Vec<u64> = b0[a2n + a1n..a2n + a1n + a0n].to_vec();
    let mut mn = [0u64; 3];
    mn[2] = mocc >> 2;
    mn[1] = mocc >> 1;
    mn[0] = mocc >> 2;
    for group in [&mut a2, &mut a1] {
        let sn = group.len();
        for i in 0..sn >> 1 {
            group.swap(i, sn - i - 1);
        }
    }
    for entry in a0.iter_mut() {
        let k = ((*entry as u32) >> 1) as usize;
        let z = &ol[k];
        let q0 = u64::from(z.x_position_end + 1 - z.x_position_start);
        let t0 = u64::from(z.y_position_end + 1 - z.y_position_start);
        let m = u64::from(u32::MAX - (q0.min(t0) as u32));
        *entry = (m << 32) | u64::from(*entry as u32);
    }
    radix_sort_u64(&mut a0);
    if mn[1] > a1n as u64 {
        mn[1] = a1n as u64;
    }
    if mn[2] > a2n as u64 {
        mn[2] = a2n as u64;
    }
    if mn[0] > a0n as u64 {
        mn[0] = a0n as u64;
    }
    mn[1] = mocc - mn[0] - mn[2];
    if mn[1] > a1n as u64 {
        mn[1] = a1n as u64;
    }
    mn[2] = mocc - mn[1] - mn[0];
    if mn[2] > a2n as u64 {
        mn[2] = a2n as u64;
    }
    mn[0] = mocc - mn[1] - mn[2];
    if mn[0] > a0n as u64 {
        mn[0] = a0n as u64;
    }
    res.clear();
    res.reserve((mn[0] + mn[1] + mn[2]) as usize);
    let MaskUlOvParts { idx, mut srt } = MaskUlOvParts::take(mk);
    push_emask_lst(&mut srt, &idx, ol, minov, unitigs, &mut a2, res, true);
    push_emask_lst(&mut srt, &idx, ol, minov, unitigs, &mut a1, res, true);
    push_emask_lst(&mut srt, &idx, ol, minov, unitigs, &mut a0, res, false);
    MaskUlOvParts::give(mk, idx, srt);
    srt_kv_emask_t(res, mk);
    let _ = erate;
}

/// Split borrow of [`MaskUlOv`]'s two fields.
struct MaskUlOvParts {
    idx: Vec<u64>,
    srt: Vec<UlOv>,
}

impl MaskUlOvParts {
    fn take(mk: &mut MaskUlOv) -> Self {
        Self {
            idx: std::mem::take(&mut mk.idx),
            srt: std::mem::take(&mut mk.srt),
        }
    }
    fn give(mk: &mut MaskUlOv, idx: Vec<u64>, srt: Vec<UlOv>) {
        mk.idx = idx;
        mk.srt = srt;
    }
}

/// `write_emask_t` (`inter.cpp` 21515): the `.ul.msk.bin` cache — a unitig
/// fingerprint, the table length, then one `(u64 n, n records)` per unitig.
/// `emask_t` is 32 bytes in the reference (28 payload + one flag byte +
/// three zero padding bytes the `memset` of `push_emask0` guarantees).
pub fn write_emask_table(
    table: &crate::assembly::ul_phase::EmaskTable,
    file_name: &str,
    unitigs: &[MaUtg],
) -> std::io::Result<()> {
    use std::io::Write;
    let path = format!("{}.ul.msk.bin", file_name);
    let mut file = std::fs::File::create(&path)?;
    crate::assembly::ul_index::write_unitig_fingerprint(unitigs, &mut file)?;
    file.write_all(&(table.len() as u32).to_le_bytes())?;
    for entry in table {
        file.write_all(&(entry.len() as u64).to_le_bytes())?;
        for mask in entry {
            let mut record = [0u8; 32];
            record[0..4].copy_from_slice(&mask.tn.to_le_bytes());
            record[4..8].copy_from_slice(&mask.rn.to_le_bytes());
            record[8..12].copy_from_slice(&mask.el.to_le_bytes());
            record[12..16].copy_from_slice(&mask.qs.to_le_bytes());
            record[16..20].copy_from_slice(&mask.qe.to_le_bytes());
            record[20..24].copy_from_slice(&mask.ts.to_le_bytes());
            record[24..28].copy_from_slice(&mask.te.to_le_bytes());
            let mut flags = mask.dir & 0x1f;
            if mask.pe {
                flags |= 0x20;
            }
            if mask.full {
                flags |= 0x40;
            }
            if mask.rev {
                flags |= 0x80;
            }
            record[28] = flags;
            file.write_all(&record)?;
        }
    }
    eprintln!("[M::write_emask_t] Index has been written.");
    Ok(())
}

/// `load_emask_t` (`inter.cpp` 21488).
pub fn load_emask_table(
    file_name: &str,
    unitigs: &[MaUtg],
) -> Option<crate::assembly::ul_phase::EmaskTable> {
    use std::io::Read;
    let path = format!("{}.ul.msk.bin", file_name);
    let mut file = std::fs::File::open(&path).ok()?;
    match crate::assembly::ul_index::test_unitig_fingerprint(unitigs, &mut file) {
        Ok(true) => {}
        _ => {
            eprintln!("[M::load_emask_t] Renew UL Index");
            return None;
        }
    }
    let mut length_buffer = [0u8; 4];
    file.read_exact(&mut length_buffer).ok()?;
    let n = u32::from_le_bytes(length_buffer) as usize;
    let mut table: crate::assembly::ul_phase::EmaskTable = Vec::with_capacity(n);
    for _ in 0..n {
        let mut count_buffer = [0u8; 8];
        file.read_exact(&mut count_buffer).ok()?;
        let count = u64::from_le_bytes(count_buffer) as usize;
        let mut entry = Vec::with_capacity(count);
        for _ in 0..count {
            let mut record = [0u8; 32];
            file.read_exact(&mut record).ok()?;
            let flags = record[28];
            entry.push(Emask {
                tn: u32::from_le_bytes(record[0..4].try_into().unwrap()),
                rn: u32::from_le_bytes(record[4..8].try_into().unwrap()),
                el: u32::from_le_bytes(record[8..12].try_into().unwrap()),
                qs: u32::from_le_bytes(record[12..16].try_into().unwrap()),
                qe: u32::from_le_bytes(record[16..20].try_into().unwrap()),
                ts: u32::from_le_bytes(record[20..24].try_into().unwrap()),
                te: u32::from_le_bytes(record[24..28].try_into().unwrap()),
                dir: flags & 0x1f,
                pe: flags & 0x20 != 0,
                full: flags & 0x40 != 0,
                rev: flags & 0x80 != 0,
            });
        }
        table.push(entry);
    }
    eprintln!("[M::load_emask_t] Index has been loaded.");
    Some(table)
}

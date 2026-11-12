//! The UL masking/phase subsystem of the rescall pass (`Correct.cpp`
//! 18050-22426): `rphase_rr` / `rphase_hl` decide, per piece of every
//! candidate overlap, whether it is a primary alignment or a shadowed
//! secondary one, using the graph-arc masks built by `mask_ovlps`.
//!
//! The reference stores the per-piece windows past `w_list.n` (counted by
//! `align_length`) and compacts them over the primary windows at the end of
//! `rphase_hl`. The port keeps them in a separate [`SecondaryWindows`] slot
//! per overlap, which is swapped together with the overlap list whenever the
//! reference swaps `overlap_region` structs.

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::ul::{get_offset_adjust, UlOv};
use crate::assembly::ul_rescall::{
    check_connect_ug, check_connect_rg, dp_commen_sketch, get_win_aln, ovlp_win_check,
    G_CHAIN_INDEL, UlRescallContext,
};
use crate::correction::recalculate::append_unmatched_windows;
use crate::correction::target::TargetSource;
use crate::correction::ul_align::{gen_err_unaligned, WindowExtz};
use crate::correction::window::{window_id_by_start, WindowList};
use crate::container::ksort::{radix_sort_by_key_32, radix_sort_u64};
use crate::overlap::paf::PafHit;
use crate::overlap::region::OverlapRegion;
use crate::reads::read_store::ReadStore;

/// `MAX_SEC_ERR` / `SEC_MODE` (`Correct.cpp` 75, `Overlaps.h` 54).
pub const SEC_MODE: u32 = 0x3fff_ffff;
/// `UL_TOPN` (`Correct.h` 23).
pub const UL_TOPN: u64 = 50;
/// `SGAP` (`Correct.h` 24).
pub const SGAP: i64 = 16;

/// `MAX_LGAP(ql)` (`Correct.h` 25).
pub fn max_lgap(query_length: i64) -> i64 {
    let scaled = query_length as f64 * 0.2;
    if scaled < 256.0 {
        scaled as i64
    } else {
        256
    }
}

/// `emask_t` (`Overlaps.h` 67): one mask entry of the graph-binning table.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Emask {
    pub tn: u32,
    pub rn: u32,
    /// `el == u32::MAX` marks "no measured difference".
    pub el: u32,
    pub qs: u32,
    pub qe: u32,
    pub ts: u32,
    pub te: u32,
    /// `dir:5`.
    pub dir: u8,
    /// `pe:1`.
    pub pe: bool,
    /// `full:1`.
    pub full: bool,
    /// `rev:1`.
    pub rev: bool,
}

/// `idx_emask_t` (`Overlaps.h` 80): one [`Emask`] list per unitig.
pub type EmaskTable = Vec<Vec<Emask>>;

/// `mask_ul_ov_t` (`Overlaps.h` 265).
#[derive(Clone, Debug, Default)]
pub struct MaskUlOv {
    /// `idx[k] = (begin << 32) | end` into `srt`, one entry per overlap.
    pub idx: Vec<u64>,
    pub srt: Vec<UlOv>,
}

/// The secondary window slots, `w_list.a[w_list.n .. w_list.n + align_length)`.
pub type SecondaryWindows = Vec<Vec<WindowList>>;

/// `ovlp_id(x)` (`Correct.h` 1404).
#[inline]
pub fn ovlp_id(x: &UlOv) -> u32 {
    x.tn
}

/// `extend_batch` (`inter.cpp` 11709).
pub fn extend_batch(z: &UlOv, ql: i64, tl: i64) -> (i64, i64, i64, i64) {
    let (mut rqs, mut rqe, mut rts, mut rte) = (-1i64, -1i64, -1i64, -1i64);
    if !z.rev {
        let q0 = i64::from(z.qs);
        let t0 = i64::from(z.ts);
        if q0 <= t0 {
            rqs = i64::from(z.qs) - q0;
            rts = i64::from(z.ts) - q0;
        } else {
            rqs = i64::from(z.qs) - t0;
            rts = i64::from(z.ts) - t0;
        }
        let q1 = ql - i64::from(z.qe);
        let t1 = tl - i64::from(z.te);
        if q1 <= t1 {
            rqe = i64::from(z.qe) + q1;
            rte = i64::from(z.te) + q1;
        } else {
            rqe = i64::from(z.qe) + t1;
            rte = i64::from(z.te) + t1;
        }
    } else {
        let q0 = i64::from(z.qs);
        let t0 = tl - i64::from(z.te);
        if q0 <= t0 {
            rqs = i64::from(z.qs) - q0;
            rte = i64::from(z.te) + q0;
        } else {
            rqs = i64::from(z.qs) - t0;
            rte = i64::from(z.te) + t0;
        }
        let q1 = ql - i64::from(z.qe);
        let t1 = i64::from(z.ts);
        if q1 <= t1 {
            rqe = i64::from(z.qe) + q1;
            rts = i64::from(z.ts) - q1;
        } else {
            rqe = i64::from(z.qe) + t1;
            rts = i64::from(z.ts) - t1;
        }
    }
    (rqs, rqe, rts, rte)
}

/// `check_ul_ov_t_consist` (`inter.cpp` 11744): `ql`/`tl` are the graph
/// vertex lengths `ug->g->seq[qn].len` / `seq[tn].len`.
pub fn check_ul_ov_t_consist(x: &UlOv, y: &UlOv, ql: i64, tl: i64, diff: f64) -> bool {
    if !(x.qn == y.qn && x.tn == y.tn && x.rev == y.rev) {
        return false;
    }
    let (xqs, xqe, xts, xte) = extend_batch(x, ql, tl);
    let (yqs, yqe, yts, yte) = extend_batch(y, ql, tl);
    let os = xqs.max(yqs);
    let oe = xqe.min(yqe);
    if oe <= os {
        return false;
    }
    let ovlp = oe - os;
    if (ovlp as f64) <= ((xqe - xqs) as f64 * (1.0 - diff)) {
        return false;
    }
    if (ovlp as f64) <= ((yqe - yqs) as f64 * (1.0 - diff)) {
        return false;
    }
    let os = xts.max(yts);
    let oe = xte.min(yte);
    if oe <= os {
        return false;
    }
    if (ovlp as f64) <= ((xte - xts) as f64 * (1.0 - diff)) {
        return false;
    }
    if (ovlp as f64) <= ((yte - yts) as f64 * (1.0 - diff)) {
        return false;
    }
    true
}

/// `infer_se` (`inter.cpp` 11171).
pub fn infer_se(
    qs: u32,
    qe: u32,
    ts: u32,
    te: u32,
    rev: bool,
    rqs: u32,
    rqe: u32,
) -> Option<(u32, u32)> {
    if rqs <= rqe && rqs >= qs && rqe <= qe {
        let mut s_shift = get_offset_adjust(rqs - qs, qe - qs, te - ts);
        let mut e_shift = get_offset_adjust(qe - rqe, qe - qs, te - ts);
        if rev {
            std::mem::swap(&mut s_shift, &mut e_shift);
        }
        let rts = ts + s_shift;
        let rte = te - e_shift;
        if rts <= rte && rts >= ts && rte <= te {
            return Some((rts, rte));
        }
    }
    None
}

/// `cal_s_interval` (`inter.cpp` 11162).
pub(crate) fn cal_s_interval(graph: &Asg, arc: &AsgArc, res: &mut UlOv) {
    *res = UlOv::default();
    res.rev = (((arc.ul >> 32) as u32 ^ arc.v) & 1) != 0;
    res.qn = (arc.ul >> 33) as u32;
    res.tn = arc.v >> 1;
    let ql = graph.sequences[res.qn as usize].length;
    let tl = graph.sequences[res.tn as usize].length;
    let overlap = arc.overlap;
    if (arc.ul >> 32) & 1 != 0 {
        res.qs = 0;
        res.qe = overlap.min(ql);
    } else {
        res.qe = ql;
        res.qs = if overlap <= res.qe { res.qe - overlap } else { 0 };
    }
    if arc.v & 1 != 0 {
        res.te = tl;
        res.ts = if overlap <= res.te { res.te - overlap } else { 0 };
    } else {
        res.ts = 0;
        res.te = overlap.min(tl);
    }
}

/// `cal_x_interval` (`inter.cpp` 11202).
pub(crate) fn cal_x_interval(graph: &Asg, v: &UlOv, w: &UlOv, res: &mut UlOv) -> bool {
    *res = UlOv::default();
    if v.tn == w.qn {
        let os = v.ts.max(w.qs);
        let oe = v.te.min(w.qe);
        if oe <= os {
            return false;
        }
        if let (Some((rqs, rqe)), Some((rts, rte))) = (
            infer_se(v.ts, v.te, v.qs, v.qe, v.rev, os, oe),
            infer_se(w.qs, w.qe, w.ts, w.te, w.rev, os, oe),
        ) {
            if rqe > rqs && rte > rts {
                res.qn = v.qn;
                res.tn = w.tn;
                res.rev = v.rev != w.rev;
                res.qs = rqs;
                res.qe = rqe;
                res.ts = rts;
                res.te = rte;
                return true;
            }
        }
    }
    false
}

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

/// `retrieve_cigar_err` (`Correct.cpp` 18050).
fn retrieve_cigar_err(ez: &WindowExtz<'_>, s: i64, e: i64, xk: &mut i64, ck: &mut i64) -> i64 {
    if ez.cells.is_empty() {
        return 0;
    }
    let cn = ez.cells.len() as i64;
    let mut err = 0i64;
    if *ck < 0 || *ck > cn {
        *ck = 0;
        *xk = i64::from(ez.text_start);
    }
    while *ck > 0 && *xk > s {
        *ck -= 1;
        let cell = ez.cells[*ck as usize];
        let op = cell >> 14;
        if op != 2 {
            *xk -= i64::from(cell & 0x3fff);
        }
    }
    while *ck < cn && *xk < e {
        let ws = *xk;
        let cell = ez.cells[*ck as usize];
        let op = cell >> 14;
        if op != 2 {
            *xk += i64::from(cell & 0x3fff);
        }
        let we = *xk;
        let os = s.max(ws);
        let oe = e.min(we);
        let mut ovlp = (oe - os).max(0);
        if op == 2 && ws >= s && ws < e {
            ovlp = i64::from(cell & 0x3fff);
        }
        *ck += 1;
        if ovlp == 0 || op == 0 {
            continue;
        }
        err += ovlp;
    }
    err
}

/// `extract_sub_cigar_err` (`Correct.cpp` 18093); `bd = Some(..)` selects the
/// `_rr` variant (`Correct.cpp` 18465), which clamps `[s, e)` by the boundary
/// padding first. The cursor triple lives in `p`: `ovlp_cur_wid` -> `qn`,
/// `ovlp_cur_xoff` -> `qs`, `ovlp_cur_coff` -> `qe`.
fn extract_sub_cigar_err_impl(
    region: &OverlapRegion,
    mut s: i64,
    mut e: i64,
    p: &mut UlOv,
    bd: Option<i64>,
) -> i64 {
    let mut wk = i64::from(p.qn);
    let mut xk = i64::from(p.qs);
    let mut ck = i64::from(p.qe);
    let min_w = i64::from(p.ts);
    let max_w = i64::from(p.te);
    if let Some(bd) = bd {
        let s0 = i64::from(region.window_cigar.windows[min_w as usize].x_start) + bd;
        let e0 = i64::from(region.window_cigar.windows[max_w as usize].x_end) + 1 - bd;
        if s < s0 {
            s = s0;
        }
        if e > e0 {
            e = e0;
        }
        if s >= e {
            return -1;
        }
    }
    let windows = &region.window_cigar.windows;
    let cells = &region.window_cigar.cells;
    let mut err = 0i64;
    let mut tot = e - s;
    if wk < min_w || wk > max_w {
        wk = min_w;
    }
    while wk >= min_w && i64::from(windows[wk as usize].x_start) > s {
        wk -= 1;
    }
    if wk < min_w || wk > max_w {
        return -1;
    }
    while wk <= max_w && i64::from(windows[wk as usize].x_end) < s {
        wk += 1;
    }
    if wk < min_w || wk > max_w {
        return -1;
    }
    if wk != i64::from(p.qn) {
        xk = i64::from(windows[wk as usize].x_start);
        ck = 0;
    }
    while wk <= max_w && i64::from(windows[wk as usize].x_start) < e {
        let m = &windows[wk as usize];
        let ws = i64::from(m.x_start);
        let we = i64::from(m.x_end) + 1;
        let os = s.max(ws);
        let oe = e.min(we);
        let ovlp = (oe - os).max(0);
        if ovlp > 0 {
            let xl = i64::from(m.x_end) + 1 - i64::from(m.x_start);
            let yl = i64::from(m.y_end) + 1 - i64::from(m.y_start);
            if is_ualn_win(m) || is_est_aln(m) {
                let mut werr = if is_ualn_win(m) {
                    gen_err_unaligned(xl, yl)
                } else {
                    i64::from(m.error)
                };
                if ovlp < xl {
                    werr = (((ovlp as f64) / (xl as f64)) * (werr as f64)) as i64;
                }
                err += werr;
                xk = we;
                ck = i64::from(m.cigar_length);
            } else if ovlp == xl {
                err += i64::from(m.error);
                xk = we;
                ck = i64::from(m.cigar_length);
            } else {
                let ez = WindowExtz::new(m, cells);
                err += retrieve_cigar_err(&ez, os, oe, &mut xk, &mut ck);
            }
        }
        tot -= ovlp;
        if xk >= e {
            break;
        }
        wk += 1;
        if wk > max_w {
            break;
        }
        xk = i64::from(windows[wk as usize].x_start);
        ck = 0;
    }
    if tot != 0 {
        eprintln!(
            "[ETOT] tot={} s={} e={} y_id={} x=[{},{}] p=[ts={},te={},sec={},rev={}] minw={} maxw={} curwid={}",
            tot,
            s,
            e,
            region.y_identifier,
            region.x_position_start,
            region.x_position_end,
            p.ts,
            p.te,
            p.sec,
            p.rev as u32,
            min_w,
            max_w,
            wk
        );
        for w in min_w..=max_w.min(windows.len() as i64 - 1) {
            let m = &windows[w as usize];
            eprintln!(
                "[ETOTW] w={} x=[{},{}] y=[{},{}] err={} clen={} ualn={} est={}",
                w,
                m.x_start,
                m.x_end,
                m.y_start,
                m.y_end,
                m.error,
                m.cigar_length,
                is_ualn_win(m) as i32,
                is_est_aln(m) as i32
            );
        }
    }
    assert_eq!(tot, 0);
    p.qn = wk as u32;
    p.qs = xk as u32;
    p.qe = ck as u32;
    err
}

/// `extract_sub_cigar_err` (`Correct.cpp` 18093).
pub fn extract_sub_cigar_err(region: &OverlapRegion, s: i64, e: i64, p: &mut UlOv) -> i64 {
    extract_sub_cigar_err_impl(region, s, e, p, None)
}

/// `extract_sub_cigar_err_rr` (`Correct.cpp` 18465).
fn extract_sub_cigar_err_rr(region: &OverlapRegion, s: i64, e: i64, p: &mut UlOv, bd: i64) -> i64 {
    extract_sub_cigar_err_impl(region, s, e, p, Some(bd))
}

/// `is_mask_ov` (`Correct.cpp` 18948): whether the best candidate `sec_id` is
/// itself a masked overlap of one of the equally-best `bes_id` overlaps.
fn is_mask_ov(mk: &MaskUlOv, bes_id: &[u64], sec_id: u32) -> bool {
    let entry = mk.idx[sec_id as usize];
    let s = (entry >> 32) as usize;
    let e = entry as u32 as usize;
    let mut bk = 0usize;
    let mut si = s;
    while bk < bes_id.len() && si < e {
        let tn = mk.srt[si].tn;
        if bes_id[bk] < u64::from(tn) {
            bk += 1;
        } else if u64::from(tn) < bes_id[bk] {
            si += 1;
        } else {
            return mk.srt[si].qs == 0;
        }
    }
    false
}

/// `push_sec_aln_robust` (`Correct.cpp` 18201): appends one secondary piece to
/// `slots[oid]`, merging into the previous piece when it ends exactly at `s`,
/// carries no cigar and the new piece has no error either.
fn push_sec_aln_robust(slots: &mut [Vec<WindowList>], oid: usize, s: i64, e: i64, sec_err: i64) {
    let slot = &mut slots[oid];
    if let Some(last) = slot.last_mut() {
        if i64::from(last.x_end) == s && last.cigar_length == 0 && sec_err == 0 {
            last.x_end = e as i32;
            return;
        }
    }
    slot.push(WindowList {
        x_start: s as i32,
        x_end: e as i32,
        y_start: 0,
        y_end: 0,
        extra_begin: 0,
        extra_end: 0,
        error: 0,
        error_threshold: 0,
        cigar_start: 0,
        cigar_length: sec_err as u32,
    });
}

/// `gen_region_phase_robust_rr` (`Correct.cpp` 18966): picks the lowest-error
/// candidate of one shared region, pushes the primary pieces and the
/// (error-charged) secondary pieces, and compacts `id_a`.
fn gen_region_phase_robust_rr(
    regions: &mut [OverlapRegion],
    id_a: &mut [u64],
    s: u64,
    e: u64,
    dp: u64,
    c_idx: &mut [UlOv],
    buf: &mut Vec<u64>,
    mk: Option<&MaskUlOv>,
    slots: &mut [Vec<WindowList>],
) -> usize {
    let mut id_n = id_a.len();
    if id_n == 0 {
        return id_n;
    }
    buf.clear();
    buf.reserve(dp as usize);
    let mut rm_n = 0usize;
    for k in 0..id_n {
        let p = c_idx[id_a[k] as usize];
        let q0 =
            i64::from(regions[ovlp_id(&p) as usize].window_cigar.windows[p.ts as usize].x_start)
                + i64::from(p.sec);
        let q1 = i64::from(
            regions[ovlp_id(&p) as usize].window_cigar.windows[p.te as usize].x_end,
        ) + 1
            - i64::from(p.sec);
        if q0 <= s as i64 && q1 >= e as i64 {
            buf.push(id_a[k]);
        }
        if q1 < e as i64 {
            rm_n += 1;
        }
    }
    let buf_n = buf.len();
    assert_eq!(buf_n, dp as usize);
    let mut mn = 0usize;
    if buf_n > 0 {
        let mut msc = i32::MAX as i64;
        let mut msc_k = -1i64;
        let mut msc_n = 0i64;
        for k in 0..buf_n {
            let cid = buf[k] as usize;
            let mut p = c_idx[cid];
            let oid = ovlp_id(&p) as usize;
            let bd = p.sec as i64;
            let err = extract_sub_cigar_err_rr(&regions[oid], s as i64, e as i64, &mut p, bd);
            c_idx[cid] = p;
            assert!(err >= 0);
            if err < msc {
                msc = err;
                msc_k = k as i64;
                msc_n = 1;
            } else if err == msc {
                msc_n += 1;
            }
            buf[k] |= (err as u64) << 32;
        }
        if msc_n == 1 {
            mn_select_single(buf, msc_k as usize);
            mn = 1;
        } else {
            mn = mn_select_many(buf, buf_n, msc as u64, msc_n as usize);
        }
        let _ = mn;
        for k in 0..mn {
            let oid = ovlp_id(&c_idx[buf[k] as u32 as usize]) as usize;
            buf[k] = u64::from(buf[k] as u32);
            push_sec_aln_robust(slots, oid, s as i64, e as i64, 0);
        }
        match mk {
            None => {
                for k in mn..buf_n {
                    let err = ((buf[k] & 0x7fff_ffff_ffff_ffff) >> 32) as i64 - msc;
                    let oid = ovlp_id(&c_idx[buf[k] as u32 as usize]) as usize;
                    push_sec_aln_robust(slots, oid, s as i64, e as i64, err);
                }
            }
            Some(mk) => {
                if mn > 1 {
                    radix_sort_u64(&mut buf[..mn]);
                }
                for k in mn..buf_n {
                    let mut err = ((buf[k] & 0x7fff_ffff_ffff_ffff) >> 32) as i64 - msc;
                    let cid = buf[k] as u32;
                    if is_mask_ov(mk, &buf[..mn], cid) {
                        err = 0;
                    }
                    let oid = ovlp_id(&c_idx[cid as usize]) as usize;
                    push_sec_aln_robust(slots, oid, s as i64, e as i64, err);
                }
            }
        }
    }
    if rm_n > 0 {
        let mut m = 0usize;
        for k in 0..id_n {
            let p = &c_idx[id_a[k] as usize];
            let q1 = i64::from(
                regions[ovlp_id(p) as usize].window_cigar.windows[p.te as usize].x_end,
            ) + 1
                - i64::from(p.sec);
            if q1 < e as i64 {
                continue;
            }
            id_a[m] = id_a[k];
            m += 1;
        }
        id_n = m;
    }
    id_n
}

/// The `msc_n == 1` branch of `gen_region_phase_robust_rr`: bring the unique
/// best candidate to the front.
fn mn_select_single(buf: &mut [u64], msc_k: usize) {
    if msc_k != 0 {
        buf.swap(msc_k, 0);
    }
}

/// The `msc_n > 1` branch: gather every equally-best candidate to the front.
fn mn_select_many(buf: &mut [u64], buf_n: usize, msc: u64, msc_n: usize) -> usize {
    let mut mn = 0usize;
    let mut k = 0usize;
    while k < buf_n && mn < msc_n {
        if (buf[k] >> 32) == msc {
            if mn != k {
                buf.swap(k, mn);
            }
            mn += 1;
        }
        k += 1;
    }
    mn
}

/// `rphase_rr` (`Correct.cpp` 19457): the sliding-window phase sweep; the
/// `idx`/`buf`/`c_idx` of the reference are scratch here.
pub fn rphase_rr(
    regions: &mut [OverlapRegion],
    bd: i64,
    mk: Option<&MaskUlOv>,
    slots: &mut [Vec<WindowList>],
) {
    let on = regions.len();
    let mut idx: Vec<u64> = Vec::with_capacity(on << 1);
    let mut c_idx: Vec<UlOv> = Vec::with_capacity(on);
    for k in 0..on {
        let zwn = regions[k].window_cigar.windows.len();
        regions[k].alignment_length = 0;
        regions[k].overlap_length = u32::MAX;
        regions[k].non_homopolymer_errors = 0;
        // The reference resets `align_length` here, which makes the secondary
        // windows pushed by the previous pass invisible; the slot accumulates
        // from scratch.
        slots[k].clear();
        if zwn == 0 {
            continue;
        }
        let windows = &regions[k].window_cigar.windows;
        let mut q = [i32::MIN as i64; 2];
        let mut t = [i32::MIN as i64; 2];
        let mut w = [0usize; 2];
        let mut started = false;
        for i in 0..zwn {
            if i64::from(windows[i].x_start) == q[1] + 1
                && i64::from(windows[i].y_start) == t[1] + 1
            {
                q[1] = i64::from(windows[i].x_end);
                t[1] = i64::from(windows[i].y_end);
                w[1] = i;
            } else {
                if started {
                    let q0 = q[0] + bd;
                    let q1 = q[1] - bd;
                    if q1 >= q0 {
                        let mut m = (q0 as u64) << 1;
                        m <<= 32;
                        m += c_idx.len() as u64;
                        idx.push(m);
                        let mut m = (((q1 as u64) << 1) + 1) << 32;
                        m += c_idx.len() as u64;
                        idx.push(m);
                        c_idx.push(UlOv {
                            tn: k as u32,
                            ts: w[0] as u32,
                            te: w[1] as u32,
                            qn: w[0] as u32,
                            qs: windows[w[0]].x_start as u32,
                            qe: 0,
                            sec: bd as u32,
                            ..UlOv::default()
                        });
                    }
                }
                q[0] = i64::from(windows[i].x_start);
                q[1] = i64::from(windows[i].x_end);
                t[0] = i64::from(windows[i].y_start);
                t[1] = i64::from(windows[i].y_end);
                w[0] = i;
                w[1] = i;
                started = true;
            }
        }
        if started {
            let q0 = q[0] + bd;
            let q1 = q[1] - bd;
            if q1 >= q0 {
                let mut m = (q0 as u64) << 1;
                m <<= 32;
                m += c_idx.len() as u64;
                idx.push(m);
                let mut m = (((q1 as u64) << 1) + 1) << 32;
                m += c_idx.len() as u64;
                idx.push(m);
                c_idx.push(UlOv {
                    tn: k as u32,
                    ts: w[0] as u32,
                    te: w[1] as u32,
                    qn: w[0] as u32,
                    qs: windows[w[0]].x_start as u32,
                    qe: 0,
                    sec: bd as u32,
                    ..UlOv::default()
                });
            }
        }
    }
    radix_sort_u64(&mut idx);

    let srt_n = idx.len();
    let mut buf: Vec<u64> = Vec::new();
    let mut dp = 0i64;
    let mut beg: i64 = 0;
    let mut end: i64 = -1;
    for i in 0..srt_n {
        let old_dp = dp;
        if (idx[i] >> 32) & 1 != 0 {
            dp -= 1;
            end = ((idx[i] >> 33) + 1) as i64;
        } else {
            dp += 1;
            end = (idx[i] >> 33) as i64;
            idx.push(idx[i] & 0xffff_ffff);
        }
        if end > beg && old_dp >= 2 {
            let survivors = gen_region_phase_robust_rr(
                regions,
                &mut idx[srt_n..],
                beg as u64,
                end as u64,
                old_dp as u64,
                &mut c_idx,
                &mut buf,
                mk,
                slots,
            );
            idx.truncate(srt_n + survivors);
        }
        beg = end;
    }
}

/// `convert_ul_ov_t` (`Correct.cpp` 19160): the overlap on the unitig's own
/// orientation, `qn` left as the sentinel `(uint32_t)-1`.
fn convert_ul_ov_t(region: &OverlapRegion, unitigs: &[MaUtg]) -> UlOv {
    let tn = region.y_identifier;
    let rev = region.y_strand != 0;
    let (ts, te) = if rev {
        let length = unitigs[tn as usize].length;
        (
            length.wrapping_sub(region.y_position_end + 1),
            length.wrapping_sub(region.y_position_start),
        )
    } else {
        (region.y_position_start, region.y_position_end + 1)
    };
    UlOv {
        qn: u32::MAX,
        qs: region.x_position_start,
        qe: region.x_position_end + 1,
        tn,
        ts,
        te,
        sec: region.non_homopolymer_errors,
        el: true,
        rev,
    }
}

/// `gen_aln_ul_ov_t` (`Correct.cpp` 20329): the window-span of one overlap as
/// an `ul_ov_t`, `qn` = the overlap id.
fn gen_aln_ul_ov_t(in_id: usize, tl: u32, region: &OverlapRegion) -> Option<UlOv> {
    let zwn = region.window_cigar.windows.len();
    if zwn == 0 {
        return None;
    }
    let windows = &region.window_cigar.windows;
    let rev = region.y_strand != 0;
    let (ts, te) = if !rev {
        (
            windows[0].y_start as u32,
            windows[zwn - 1].y_end as u32 + 1,
        )
    } else {
        (
            tl - (windows[zwn - 1].y_end as u32 + 1),
            tl - windows[0].y_start as u32,
        )
    };
    Some(UlOv {
        qn: in_id as u32,
        tn: region.y_identifier,
        rev,
        qs: windows[0].x_start as u32,
        qe: windows[zwn - 1].x_end as u32 + 1,
        ts,
        te,
        ..UlOv::default()
    })
}

/// `push_emask_flt` (`Correct.cpp` 20350): copies the mask entries of `qn`
/// that also appear in the sorted `flt` filter into `res`; the entry's `sec`
/// records its index inside the mask table.
fn push_emask_flt(
    masks: &[Emask],
    flt: &[u64],
    qn: u64,
    res: &mut Vec<UlOv>,
) -> usize {
    let rn0 = res.len();
    let mut k = 0usize;
    let mut z = 0usize;
    while k < masks.len() && z < flt.len() {
        if masks[k].tn < (flt[z] >> 32) as u32 {
            k += 1;
        } else if ((flt[z] >> 32) as u32) < masks[k].tn {
            z += 1;
        } else {
            let tn = masks[k].tn;
            while k < masks.len() && masks[k].tn == tn {
                res.push(UlOv {
                    qn: qn as u32,
                    qs: masks[k].qs,
                    qe: masks[k].qe,
                    tn: masks[k].tn,
                    ts: masks[k].ts,
                    te: masks[k].te,
                    rev: masks[k].rev,
                    el: masks[k].full,
                    sec: k as u32,
                });
                k += 1;
            }
            while z < flt.len() && (flt[z] >> 32) as u32 == tn {
                z += 1;
            }
        }
    }
    res.len() - rn0
}

/// `is_ul_ov_pe` (`Correct.cpp` 20362): `z[i]` and `z[i+1]` are the two halves
/// of one paired mask record.
fn is_ul_ov_pe(z: &[UlOv], zn: usize, i: usize, mm: &EmaskTable) -> bool {
    if i + 1 >= zn {
        return false;
    }
    let a = &z[i];
    let b = &z[i + 1];
    b.qn == a.qn
        && b.tn == a.tn
        && b.rev == a.rev
        && b.sec == a.sec + 1
        && (a.sec as usize) < mm[a.qn as usize].len()
        && (b.sec as usize) < mm[b.qn as usize].len()
        && mm[a.qn as usize][a.sec as usize].pe
        && mm[b.qn as usize][b.sec as usize].pe
        && mm[a.qn as usize][a.sec as usize].dir == 0
        && mm[b.qn as usize][b.sec as usize].dir == 1
}

/// `is_exact_ov` (`Correct.cpp` 20366).
fn is_exact_ov(z: &UlOv, mm: &EmaskTable) -> bool {
    (z.sec as usize) < mm[z.qn as usize].len()
        && !mm[z.qn as usize][z.sec as usize].pe
        && mm[z.qn as usize][z.sec as usize].el == u32::MAX
}

/// `is_cover_ul_ov_t` (`Correct.cpp` 20381); `ref_n <= 2`.
fn is_cover_ul_ov_t(
    ql: i64,
    tl: i64,
    diff: f64,
    reference: &[UlOv],
    inside: &UlOv,
) -> bool {
    if !check_ul_ov_t_consist(inside, &reference[0], ql, tl, diff) {
        return false;
    }
    if reference.len() > 1 && !check_ul_ov_t_consist(inside, &reference[reference.len() - 1], ql, tl, 0.06) {
        return false;
    }
    let mut sql = i64::from(inside.qe) - i64::from(inside.qs);
    let mut stl = i64::from(inside.te) - i64::from(inside.ts);
    for r in reference {
        if sql <= 0 || stl <= 0 {
            break;
        }
        let os = inside.qs.max(r.qs);
        let oe = inside.qe.min(r.qe);
        let ovlp = if oe > os { i64::from(oe - os) } else { 0 };
        sql -= ovlp;
        let os = inside.ts.max(r.ts);
        let oe = inside.te.min(r.te);
        let ovlp = if oe > os { i64::from(oe - os) } else { 0 };
        stl -= ovlp;
    }
    let q_len = i64::from(inside.qe) - i64::from(inside.qs);
    let t_len = i64::from(inside.te) - i64::from(inside.ts);
    if sql > 256 && (sql as f64) > q_len as f64 * 0.06 {
        return false;
    }
    if sql <= 256 && (sql as f64) > q_len as f64 * 0.6 {
        return false;
    }
    if stl > 256 && (stl as f64) > t_len as f64 * 0.06 {
        return false;
    }
    if stl <= 256 && (stl as f64) > t_len as f64 * 0.6 {
        return false;
    }
    true
}

/// `dedup_src_shared1` (`Correct.cpp` 20432): drops the entries of the query
/// segment that are covered by a mask record of the target segment.
fn dedup_src_shared1(
    graph: &Asg,
    res: &mut Vec<UlOv>,
    tocc: usize,
    qocc: usize,
    mm: &EmaskTable,
) {
    if tocc == 0 || qocc == 0 {
        return;
    }
    let rn = res.len();
    let ti = rn - tocc - qocc;
    let qi = rn - qocc;
    let tn = ti + tocc;
    let qn = qi + qocc;
    let t_seg: Vec<UlOv> = res[ti..tn].to_vec();
    let q_seg: Vec<UlOv> = res[qi..qn].to_vec();
    res.truncate(ti);

    let vertex_len = |id: u32| i64::from(graph.sequences[id as usize].length);
    let mut ff = false;
    let mut a = 0usize; // cursor into t_seg
    let mut b = 0usize; // cursor into q_seg
    while a < t_seg.len() && b < q_seg.len() {
        if t_seg[a].tn < q_seg[b].tn {
            let entry = t_seg[a];
            res.push(entry);
            a += 1;
        } else if q_seg[b].tn < t_seg[a].tn {
            let entry = q_seg[b];
            res.push(entry);
            b += 1;
        } else {
            let id = t_seg[a].tn;
            let ts = a;
            while a < t_seg.len() && t_seg[a].tn == id {
                let entry = t_seg[a];
                res.push(entry);
                a += 1;
            }
            let te = a; // relative end of the t run inside t_seg
            while b < q_seg.len() && q_seg[b].tn == id {
                let mut k = ts;
                let mut covered = false;
                while k < te {
                    let l = k;
                    let mut kk = k;
                    if is_ul_ov_pe(&t_seg[..te], te, kk, mm) {
                        kk += 1;
                    }
                    k = kk;
                    if is_cover_ul_ov_t(
                        vertex_len(t_seg[l].qn),
                        vertex_len(t_seg[l].tn),
                        0.06,
                        &t_seg[l..=kk],
                        &q_seg[b],
                    ) {
                        covered = true;
                        break;
                    }
                    k += 1;
                }
                if !covered {
                    let entry = q_seg[b];
                    res.push(entry);
                } else {
                    ff = true;
                }
                b += 1;
            }
        }
    }
    while a < t_seg.len() {
        let entry = t_seg[a];
        res.push(entry);
        a += 1;
    }
    while b < q_seg.len() {
        let entry = q_seg[b];
        res.push(entry);
        b += 1;
    }

    if !ff {
        res.truncate(rn);
    } else if res.len() > rn {
        let l = res.len();
        let mut write = rn - tocc - qocc;
        for k in rn..l {
            res[write] = res[k];
            write += 1;
        }
        res.truncate(write);
    }
}

/// The `aln2ov` macro (`Correct.cpp` 20492).
fn aln2ov(region: &OverlapRegion, tl: u32) -> UlOv {
    let rev = region.y_strand != 0;
    let (ts, te) = if rev {
        (tl - (region.y_position_end + 1), tl - region.y_position_start)
    } else {
        (region.y_position_start, region.y_position_end + 1)
    };
    UlOv {
        qn: region.x_identifier,
        tn: region.y_identifier,
        rev,
        qs: region.x_position_start,
        qe: region.x_position_end + 1,
        ts,
        te,
        ..UlOv::default()
    }
}

/// `cal_x_ul_ovlp` (`Correct.cpp` 20498).
fn cal_x_ul_ovlp(graph: &Asg, q: &OverlapRegion, t: &OverlapRegion, res: &mut UlOv) -> bool {
    let ql = graph.sequences[q.y_identifier as usize].length;
    let tl = graph.sequences[t.y_identifier as usize].length;
    let qo = aln2ov(q, ql);
    let to = aln2ov(t, tl);
    *res = UlOv::default();
    let os = qo.qs.max(to.qs);
    let oe = qo.qe.min(to.qe);
    if oe <= os {
        return false;
    }
    if let (Some((rqs, rqe)), Some((rts, rte))) = (
        infer_se(qo.qs, qo.qe, qo.ts, qo.te, qo.rev, os, oe),
        infer_se(to.qs, to.qe, to.ts, to.te, to.rev, os, oe),
    ) {
        if rqe > rqs && rte > rts {
            res.qn = qo.tn;
            res.tn = to.tn;
            res.rev = qo.rev != to.rev;
            res.qs = rqs;
            res.qe = rqe;
            res.ts = rts;
            res.te = rte;
            return true;
        }
    }
    false
}

/// `check_mask_exist` (`Correct.cpp` 20516).
#[allow(clippy::too_many_arguments)]
fn check_mask_exist(
    graph: &Asg,
    q: &OverlapRegion,
    t: &OverlapRegion,
    ov_db: &[UlOv],
    mm: &EmaskTable,
    len_diff: f64,
    is_exact: &mut u64,
    ovdb_idx: &mut u64,
) -> bool {
    *is_exact = 0;
    *ovdb_idx = u64::MAX;
    let mut rr = UlOv::default();
    if !cal_x_ul_ovlp(graph, q, t, &mut rr) {
        return false;
    }
    let vertex_len = |id: u32| i64::from(graph.sequences[id as usize].length);
    let mut ss = q.overlap_length as usize;
    if q.x_strand != 0 {
        let eid = rr.tn;
        let a = &ov_db[ss..];
        let an = q.x_strand as usize;
        let mut k = 0usize;
        while k < an && a[k].tn < eid {
            k += 1;
        }
        while k < an && a[k].tn == eid {
            if !check_ul_ov_t_consist(&rr, &a[k], vertex_len(rr.qn), vertex_len(rr.tn), len_diff) {
                k += 1;
                continue;
            }
            let mut kk = k;
            if is_ul_ov_pe(a, an, kk, mm) {
                kk += 1;
                if !check_ul_ov_t_consist(&rr, &a[kk], vertex_len(rr.qn), vertex_len(rr.tn), len_diff) {
                    k = kk + 1;
                    continue;
                }
            }
            if is_exact_ov(&a[kk], mm) {
                *is_exact = 1;
            }
            *ovdb_idx = (ss + kk) as u64;
            return true;
        }
    }

    std::mem::swap(&mut rr.qn, &mut rr.tn);
    std::mem::swap(&mut rr.qs, &mut rr.ts);
    std::mem::swap(&mut rr.qe, &mut rr.te);
    ss = t.overlap_length as usize;
    if t.x_strand != 0 {
        let eid = rr.tn;
        let a = &ov_db[ss..];
        let an = t.x_strand as usize;
        let mut k = 0usize;
        while k < an && a[k].tn < eid {
            k += 1;
        }
        while k < an && a[k].tn == eid {
            if !check_ul_ov_t_consist(&rr, &a[k], vertex_len(rr.qn), vertex_len(rr.tn), len_diff) {
                k += 1;
                continue;
            }
            let mut kk = k;
            if is_ul_ov_pe(a, an, kk, mm) {
                kk += 1;
                if !check_ul_ov_t_consist(&rr, &a[kk], vertex_len(rr.qn), vertex_len(rr.tn), len_diff) {
                    k = kk + 1;
                    continue;
                }
            }
            if is_exact_ov(&a[kk], mm) {
                *is_exact = 1;
            }
            *ovdb_idx = (ss + kk) as u64;
            return true;
        }
    }
    false
}

/// `push_consist_ul_ovlps` (`Correct.cpp` 20555).
#[allow(clippy::too_many_arguments)]
fn push_consist_ul_ovlps(
    graph: &Asg,
    ov_db: &[UlOv],
    s: usize,
    e: usize,
    len_diff: f64,
    reference: &UlOv,
    mm: &EmaskTable,
    out: &mut Vec<UlOv>,
    flip_res: bool,
    is_exact: &mut u64,
) {
    *is_exact = 0;
    if s >= e {
        return;
    }
    let ql = i64::from(graph.sequences[reference.qn as usize].length);
    let tl = i64::from(graph.sequences[reference.tn as usize].length);
    let mut k = s;
    while k < e && ov_db[k].tn < reference.tn {
        k += 1;
    }
    while k < e && ov_db[k].tn == reference.tn {
        if !check_ul_ov_t_consist(reference, &ov_db[k], ql, tl, len_diff) {
            k += 1;
            continue;
        }
        let mut z = k;
        let mut kk = k;
        if is_ul_ov_pe(ov_db, ov_db.len(), kk, mm) {
            kk += 1;
            if !check_ul_ov_t_consist(reference, &ov_db[kk], ql, tl, len_diff) {
                k = kk + 1;
                continue;
            }
        }
        if is_exact_ov(&ov_db[kk], mm) {
            *is_exact = 1;
            return;
        }
        while z <= kk {
            let mut p = ov_db[kk];
            if flip_res {
                std::mem::swap(&mut p.qn, &mut p.tn);
                std::mem::swap(&mut p.qs, &mut p.ts);
                std::mem::swap(&mut p.qe, &mut p.te);
            }
            out.push(p);
            z += 1;
        }
        k += 1;
    }
}

/// `cal_no_bd_coor` (`Correct.cpp` 20583).
fn cal_no_bd_coor(inside: &UlOv, out: &mut UlOv, mm: &EmaskTable, bd: i64) -> bool {
    let mut qs = i64::from(inside.qs);
    let mut qe = i64::from(inside.qe);
    let mut ts = i64::from(inside.ts);
    let mut te = i64::from(inside.te);
    if inside.sec != SEC_MODE {
        let entry = mm[inside.qn as usize].get(inside.sec as usize);
        let matched = match entry {
            Some(z) => {
                inside.tn == z.tn
                    && inside.rev == z.rev
                    && inside.qs == z.qs
                    && inside.qe == z.qe
                    && inside.ts == z.ts
                    && inside.te == z.te
            }
            None => false,
        };
        if matched {
            let z = &mm[inside.qn as usize][inside.sec as usize];
            assert!(z.el != u32::MAX);
            if z.dir == 0 {
                qe -= i64::from(z.el >> 1);
                if inside.rev {
                    ts += i64::from(z.el >> 1);
                } else {
                    te -= i64::from(z.el >> 1);
                }
            } else {
                qs += i64::from(z.el >> 1);
                if inside.rev {
                    te -= i64::from(z.el >> 1);
                } else {
                    ts += i64::from(z.el >> 1);
                }
            }
        } else {
            let z = &mm[inside.tn as usize][inside.sec as usize];
            assert!(
                (inside.sec as usize) < mm[inside.tn as usize].len()
                    && inside.qn == z.tn
                    && inside.rev == z.rev
                    && inside.qs == z.ts
                    && inside.qe == z.te
                    && inside.ts == z.qs
                    && inside.te == z.qe
            );
            assert!(z.el != u32::MAX);
            if z.dir == 1 {
                qe -= i64::from(z.el >> 1);
                if inside.rev {
                    ts += i64::from(z.el >> 1);
                } else {
                    te -= i64::from(z.el >> 1);
                }
            } else {
                qs += i64::from(z.el >> 1);
                if inside.rev {
                    te -= i64::from(z.el >> 1);
                } else {
                    ts += i64::from(z.el >> 1);
                }
            }
        }
    } else {
        qs += bd;
        qe -= bd;
        ts += bd;
        te -= bd;
    }
    if qe > qs && te > ts {
        *out = *inside;
        out.qs = qs as u32;
        out.qe = qe as u32;
        out.ts = ts as u32;
        out.te = te as u32;
        return true;
    }
    false
}

/// `retrieve_cigar_xcoord` (`Correct.cpp` 20742).
fn retrieve_cigar_xcoord(
    ez: &WindowExtz<'_>,
    is: i64,
    ie: i64,
    yk: &mut i64,
    xk: &mut i64,
    ck: &mut i64,
    rxs: &mut i64,
    rxe: &mut i64,
) {
    if ez.cells.is_empty() {
        return;
    }
    let cn = ez.cells.len() as i64;
    if *ck < 0 || *ck > cn {
        *ck = 0;
        *yk = i64::from(ez.pattern_start);
        *xk = i64::from(ez.text_start);
    }
    while *ck > 0 && *yk > is {
        *ck -= 1;
        let cell = ez.cells[*ck as usize];
        let op = cell >> 14;
        if op != 3 {
            *yk -= i64::from(cell & 0x3fff);
        }
        if op != 2 {
            *xk -= i64::from(cell & 0x3fff);
        }
    }
    while *ck < cn && *yk < ie {
        let ws = *yk;
        let mut ws0 = *xk;
        let cell = ez.cells[*ck as usize];
        let op = cell >> 14;
        if op != 3 {
            *yk += i64::from(cell & 0x3fff);
        }
        if op != 2 {
            *xk += i64::from(cell & 0x3fff);
        }
        let we = *yk;
        let mut we0 = *xk;
        if ws >= is && we <= ie {
            if *rxs > ws0 {
                *rxs = ws0;
            }
            if *rxe < we0 {
                *rxe = we0;
            }
        } else if op == 0 || op == 1 {
            let os = is.max(ws);
            let oe = ie.min(we);
            if oe > os {
                ws0 += os - ws;
                we0 = ws0 + (oe - os);
                if *rxs > ws0 {
                    *rxs = ws0;
                }
                if *rxe < we0 {
                    *rxe = we0;
                }
            }
        }
        *ck += 1;
    }
}

/// `extract_xcoordates0` (`Correct.cpp` 20780): the target-coordinate cursor
/// variant; `ovlp_cur_xoff` carries `yk` and `ovlp_id` carries `xk`.
fn extract_xcoordates0(
    region: &OverlapRegion,
    ys: i64,
    ye: i64,
    rxs: &mut i64,
    rxe: &mut i64,
    p: &mut UlOv,
) -> bool {
    let mut wk = i64::from(p.qn);
    let mut yk = i64::from(p.qs);
    let mut xk = i64::from(p.tn);
    let mut ck = i64::from(p.qe);
    let min_w = i64::from(p.ts);
    let max_w = i64::from(p.te);
    *rxs = i32::MAX as i64;
    *rxe = -1;
    if ys >= ye {
        return false;
    }
    let windows = &region.window_cigar.windows;
    let cells = &region.window_cigar.cells;
    let mut tot = ye - ys;
    if wk < min_w || wk > max_w {
        wk = min_w;
    }
    while wk >= min_w && i64::from(windows[wk as usize].y_start) > ys {
        wk -= 1;
    }
    if wk < min_w || wk > max_w {
        return false;
    }
    while wk <= max_w && i64::from(windows[wk as usize].y_end) < ys {
        wk += 1;
    }
    if wk < min_w || wk > max_w {
        return false;
    }
    if wk != i64::from(p.qn) {
        yk = i64::from(windows[wk as usize].y_start);
        xk = i64::from(windows[wk as usize].x_start);
        ck = 0;
    }
    while wk <= max_w && i64::from(windows[wk as usize].y_start) < ye {
        let m = &windows[wk as usize];
        let ws = i64::from(m.y_start);
        let we = i64::from(m.y_end) + 1;
        let os = ys.max(ws);
        let oe = ye.min(we);
        let ovlp = (oe - os).max(0);
        if ovlp > 0 {
            let yl = i64::from(m.y_end) + 1 - i64::from(m.y_start);
            if is_ualn_win(m) || is_est_aln(m) {
                if ovlp >= yl {
                    if *rxs > i64::from(m.x_start) {
                        *rxs = i64::from(m.x_start);
                    }
                    if *rxe < i64::from(m.x_end) + 1 {
                        *rxe = i64::from(m.x_end) + 1;
                    }
                }
                yk = we;
                xk = we;
                ck = i64::from(m.cigar_length);
            } else if ovlp == yl {
                yk = we;
                xk = we;
                ck = i64::from(m.cigar_length);
                if *rxs > i64::from(m.x_start) {
                    *rxs = i64::from(m.x_start);
                }
                if *rxe < i64::from(m.x_end) + 1 {
                    *rxe = i64::from(m.x_end) + 1;
                }
            } else {
                let ez = WindowExtz::new(m, cells);
                retrieve_cigar_xcoord(&ez, os, oe, &mut yk, &mut xk, &mut ck, rxs, rxe);
            }
        }
        tot -= ovlp;
        if yk >= ye {
            break;
        }
        wk += 1;
        if wk > max_w {
            break;
        }
        yk = i64::from(windows[wk as usize].y_start);
        xk = i64::from(windows[wk as usize].x_start);
        ck = 0;
    }
    assert_eq!(tot, 0);
    p.qn = wk as u32;
    p.qs = yk as u32;
    p.tn = xk as u32;
    p.qe = ck as u32;
    *rxe > *rxs
}

/// `extract_xcoordates` (`Correct.cpp` 20839): the query intervals of `z` that
/// overlap the sorted interval list `a`, expressed in target coordinates.
fn extract_xcoordates(region: &OverlapRegion, a: &[u64], b: &mut Vec<u64>) -> usize {
    let bn = b.len();
    let zwn = region.window_cigar.windows.len();
    if zwn == 0 {
        return 0;
    }
    let windows = &region.window_cigar.windows;
    let mut q = [i32::MIN as i64; 2];
    let mut t = [i32::MIN as i64; 2];
    let mut w = [0usize; 2];
    let mut started = false;
    let mut m = UlOv::default();
    let mut a_i = 0usize;
    let mut a_z;
    for i in 0..zwn {
        if i64::from(windows[i].x_start) == q[1] + 1
            && i64::from(windows[i].y_start) == t[1] + 1
        {
            q[1] = i64::from(windows[i].x_end);
            t[1] = i64::from(windows[i].y_end);
            w[1] = i;
        } else {
            if started {
                m.tn = windows[w[0]].x_start as u32;
                m.ts = w[0] as u32;
                m.te = w[1] as u32;
                m.qn = w[0] as u32;
                m.qs = windows[w[0]].y_start as u32;
                m.qe = 0;

                let is = t[0];
                let ie = t[1] + 1;
                a_z = a_i as i64;
                while a_z >= 0 {
                    let ae = a[a_z as usize] as u32 as i64;
                    if ae <= is {
                        break;
                    }
                    a_z -= 1;
                }
                if a_z < 0 {
                    a_z = 0;
                }
                while (a_z as usize) < a.len() {
                    let as_ = (a[a_z as usize] >> 32) as i64;
                    let ae = a[a_z as usize] as u32 as i64;
                    let os = is.max(as_);
                    let oe = ie.min(ae);
                    let ovlp = (oe - os).max(0);
                    if ovlp > 0 {
                        let mut qs = 0i64;
                        let mut qe = 0i64;
                        if extract_xcoordates0(region, os, oe, &mut qs, &mut qe, &mut m) {
                            b.push(((qs as u64) << 32) | qe as u64);
                        }
                    }
                    if as_ >= ie {
                        break;
                    }
                    a_z += 1;
                }
                a_i = a_z.max(0) as usize;
            }
            q[0] = i64::from(windows[i].x_start);
            q[1] = i64::from(windows[i].x_end);
            t[0] = i64::from(windows[i].y_start);
            t[1] = i64::from(windows[i].y_end);
            w[0] = i;
            w[1] = i;
            started = true;
        }
    }
    if started {
        m.tn = windows[w[0]].x_start as u32;
        m.ts = w[0] as u32;
        m.te = w[1] as u32;
        m.qn = w[0] as u32;
        m.qs = windows[w[0]].y_start as u32;
        m.qe = 0;

        let is = t[0];
        let ie = t[1] + 1;
        a_z = a_i as i64;
        while a_z >= 0 {
            let ae = a[a_z as usize] as u32 as i64;
            if ae <= is {
                break;
            }
            a_z -= 1;
        }
        if a_z < 0 {
            a_z = 0;
        }
        while (a_z as usize) < a.len() {
            let as_ = (a[a_z as usize] >> 32) as i64;
            let ae = a[a_z as usize] as u32 as i64;
            let os = is.max(as_);
            let oe = ie.min(ae);
            let ovlp = (oe - os).max(0);
            if ovlp > 0 {
                let mut qs = 0i64;
                let mut qe = 0i64;
                if extract_xcoordates0(region, os, oe, &mut qs, &mut qe, &mut m) {
                    b.push(((qs as u64) << 32) | qe as u64);
                }
            }
            if as_ >= ie {
                break;
            }
            a_z += 1;
        }
        let _ = a_i;
    }
    b.len() - bn
}

/// `cal_xerr` (`Correct.cpp` 20914): the total error of `z` inside the sorted
/// interval list `a`.
fn cal_xerr(region: &OverlapRegion, a: &[u64]) -> i64 {
    let zwn = region.window_cigar.windows.len();
    if zwn == 0 {
        return 0;
    }
    let windows = &region.window_cigar.windows;
    let mut q = [i32::MIN as i64; 2];
    let mut t = [i32::MIN as i64; 2];
    let mut w = [0usize; 2];
    let mut started = false;
    let mut err = 0i64;
    let mut a_i = 0usize;
    let mut a_z;
    for i in 0..zwn {
        if i64::from(windows[i].x_start) == q[1] + 1
            && i64::from(windows[i].y_start) == t[1] + 1
        {
            q[1] = i64::from(windows[i].x_end);
            t[1] = i64::from(windows[i].y_end);
            w[1] = i;
        } else {
            if started {
                let mut m = UlOv::default();
                m.tn = 0;
                m.ts = w[0] as u32;
                m.te = w[1] as u32;
                m.qn = w[0] as u32;
                m.qs = windows[w[0]].x_start as u32;
                m.qe = 0;
                m.sec = 0;

                let is = q[0];
                let ie = q[1] + 1;
                a_z = a_i as i64;
                while a_z >= 0 {
                    let ae = a[a_z as usize] as u32 as i64;
                    if ae <= is {
                        break;
                    }
                    a_z -= 1;
                }
                if a_z < 0 {
                    a_z = 0;
                }
                while (a_z as usize) < a.len() {
                    let as_ = (a[a_z as usize] >> 32) as i64;
                    let ae = a[a_z as usize] as u32 as i64;
                    let os = is.max(as_);
                    let oe = ie.min(ae);
                    let ovlp = (oe - os).max(0);
                    if ovlp > 0 {
                        err += extract_sub_cigar_err(region, os, oe, &mut m);
                    }
                    if as_ >= ie {
                        break;
                    }
                    a_z += 1;
                }
                a_i = a_z.max(0) as usize;
            }
            q[0] = i64::from(windows[i].x_start);
            q[1] = i64::from(windows[i].x_end);
            t[0] = i64::from(windows[i].y_start);
            t[1] = i64::from(windows[i].y_end);
            w[0] = i;
            w[1] = i;
            started = true;
        }
    }
    if started {
        let mut m = UlOv::default();
        m.tn = 0;
        m.ts = w[0] as u32;
        m.te = w[1] as u32;
        m.qn = w[0] as u32;
        m.qs = windows[w[0]].x_start as u32;
        m.qe = 0;
        m.sec = 0;

        let is = q[0];
        let ie = q[1] + 1;
        a_z = a_i as i64;
        while a_z >= 0 {
            let ae = a[a_z as usize] as u32 as i64;
            if ae <= is {
                break;
            }
            a_z -= 1;
        }
        if a_z < 0 {
            a_z = 0;
        }
        while (a_z as usize) < a.len() {
            let as_ = (a[a_z as usize] >> 32) as i64;
            let ae = a[a_z as usize] as u32 as i64;
            let os = is.max(as_);
            let oe = ie.min(ae);
            let ovlp = (oe - os).max(0);
            if ovlp > 0 {
                err += extract_sub_cigar_err(region, os, oe, &mut m);
            }
            if as_ >= ie {
                break;
            }
            a_z += 1;
        }
    }
    err
}

/// `update_masks` (`Correct.cpp` 20993): computes the mutually un-masked
/// intervals of the two overlaps' unaligned regions; `in` carries the pair
/// intervals in and the masked regions out.
fn update_masks(input: &mut Vec<u64>, q: &OverlapRegion, t: &OverlapRegion, bd: i64) {
    let in_n0 = input.len();
    for k in 0..in_n0 {
        let s = input[k] >> 32;
        let e = input[k] as u32 as u64;
        input.push(s << 1);
        input.push((((e - 1) << 1) | 1));
    }
    let in_n1 = input.len();

    radix_sort_u64(&mut input[in_n0..]);
    let mut start = 0u64;
    let mut dp = 0i64;
    for k in in_n0..in_n1 {
        let old_dp = dp;
        if input[k] & 1 != 0 {
            dp -= 1;
        } else {
            dp += 1;
        }
        if old_dp < 2 && dp >= 2 {
            start = input[k] >> 1;
        } else if old_dp >= 2 && dp < 2 {
            let end = input[k] >> 1;
            if end >= start {
                input.push((start << 32) | end);
            }
        }
    }
    let compacted: Vec<u64> = input[in_n1..].to_vec();
    input.clear();
    input.extend_from_slice(&compacted);
    let in_n0 = input.len();

    for z in [&q, &t] {
        let zwn = z.window_cigar.windows.len();
        if zwn == 0 {
            continue;
        }
        let windows = &z.window_cigar.windows;
        let mut qv = [i32::MIN as i64; 2];
        let mut tv = [i32::MIN as i64; 2];
        let mut started = false;
        for i in 0..zwn {
            if i64::from(windows[i].x_start) == qv[1] + 1
                && i64::from(windows[i].y_start) == tv[1] + 1
            {
                qv[1] = i64::from(windows[i].x_end);
                tv[1] = i64::from(windows[i].y_end);
            } else {
                if started {
                    let q0 = qv[0] + bd;
                    let q1 = qv[1] - bd;
                    if q1 >= q0 {
                        let m = ((q0 as u64) << 1) << 32;
                        input.push(m);
                        let m = ((((q1 as u64) << 1) + 1) << 32);
                        input.push(m);
                    }
                }
                qv[0] = i64::from(windows[i].x_start);
                qv[1] = i64::from(windows[i].x_end);
                tv[0] = i64::from(windows[i].y_start);
                tv[1] = i64::from(windows[i].y_end);
                started = true;
            }
        }
        if started {
            let q0 = qv[0] + bd;
            let q1 = qv[1] - bd;
            if q1 >= q0 {
                let m = ((q0 as u64) << 1) << 32;
                input.push(m);
                let m = ((((q1 as u64) << 1) + 1) << 32);
                input.push(m);
            }
        }
    }

    for k in 0..in_n0 {
        let s = input[k] >> 32;
        let e = input[k] as u32 as u64;
        if e >= s {
            let m = (((s as u64) << 1) << 32) | 1;
            input.push(m);
            let m = ((((e as u64) << 1) + 1) << 32) | 1;
            input.push(m);
        }
    }
    let in_n1 = input.len();

    radix_sort_u64(&mut input[in_n0..]);
    let mut dp = 0i64;
    let mut dp_mask = 0i64;
    let mut old_dp_mask = 0i64;
    let mut start = 0u64;
    let mut end: i64 = -1;
    for k in in_n0..in_n1 {
        let old_dp = dp;
        old_dp_mask = dp_mask;
        if (input[k] >> 32) & 1 != 0 {
            dp -= 1;
            end = ((input[k] >> 33) + 1) as i64;
            if input[k] as u32 != 0 {
                dp_mask -= 1;
            }
        } else {
            dp += 1;
            end = (input[k] >> 33) as i64;
            if input[k] as u32 != 0 {
                dp_mask += 1;
            }
        }
        if end > start as i64 && old_dp >= 2 && old_dp_mask <= 0 {
            input.push((start << 32) | end as u64);
        }
        start = end as u64;
    }
    let compacted: Vec<u64> = input[in_n1..].to_vec();
    input.clear();
    input.extend_from_slice(&compacted);
}

/// `cal_paired_distance` (`Correct.cpp` 21132).
#[allow(clippy::too_many_arguments)]
fn cal_paired_distance(
    graph: &Asg,
    q: &OverlapRegion,
    t: &OverlapRegion,
    a: &[UlOv],
    coor_srt: &mut Vec<u64>,
    buf1: &mut Vec<u64>,
    mm: &EmaskTable,
    bd: i64,
) -> i64 {
    let a_n = a.len();
    coor_srt.clear();
    coor_srt.reserve(a_n << 1);
    let mut ou = UlOv::default();
    for entry in a {
        if !cal_no_bd_coor(entry, &mut ou, mm, bd) {
            continue;
        }
        coor_srt.push(((ou.qs as u64) << 32) | ou.qe as u64);
        coor_srt.push(((ou.ts as u64) << 32) | ou.te as u64);
    }
    // The reference splits `srt` into `qa = srt->a` and `ta = qa + a_n`.
    let cn = coor_srt.len() / 2;
    let mut qa: Vec<u64> = coor_srt[..cn].to_vec();
    let mut ta: Vec<u64> = coor_srt[cn..].to_vec();

    radix_sort_u64(&mut qa);
    let mut m = merge_compact(&mut qa);
    let mut rev = q.y_strand != 0;
    let mut tl = u64::from(graph.sequences[q.y_identifier as usize].length);
    if rev {
        mirror_intervals(&mut qa, m, tl);
    }
    let qn = m;

    radix_sort_u64(&mut ta);
    m = merge_compact(&mut ta);
    rev = t.y_strand != 0;
    tl = u64::from(graph.sequences[t.y_identifier as usize].length);
    if rev {
        mirror_intervals(&mut ta, m, tl);
    }
    let tn = m;

    buf1.clear();
    let mut qocc = 0usize;
    let mut tocc = 0usize;
    if qn > 0 && tn > 0 {
        qocc = extract_xcoordates(q, &qa, buf1);
        tocc = extract_xcoordates(t, &ta, buf1);
        if qocc == 0 || tocc == 0 {
            buf1.clear();
        }
    }
    if buf1.is_empty() {
        return i32::MIN as i64;
    }
    update_masks(buf1, q, t, bd);

    let mut eq = 0i64;
    let mut et = 0i64;
    if !buf1.is_empty() {
        eq = cal_xerr(q, buf1);
        et = cal_xerr(t, buf1);
    }
    eq - et
}

/// The sort-and-merge compaction of `cal_paired_distance` / the mirror halves.
fn merge_compact(ca: &mut Vec<u64>) -> usize {
    let mut m = 0usize;
    for k in 0..ca.len() {
        if m == 0 || (ca[m - 1] as u32) < (ca[k] >> 32) as u32 {
            ca[m] = ca[k];
            m += 1;
        } else if (ca[m - 1] as u32) < (ca[k] as u32) {
            ca[m - 1] += u64::from((ca[k] as u32) - (ca[m - 1] as u32));
        }
    }
    ca.truncate(m);
    m
}

/// The reverse-strand mirror of `cal_paired_distance`.
fn mirror_intervals(ca: &mut [u64], m: usize, tl: u64) {
    let rev_n = m >> 1;
    for k in 0..rev_n {
        let mt = ca[k];
        ca[k] = ca[m - k - 1];
        ca[m - k - 1] = mt;

        let s = tl - (ca[k] as u32) as u64;
        let e = tl - (ca[k] >> 32);
        ca[k] = (s << 32) | e;
        let s = tl - (ca[m - k - 1] as u32) as u64;
        let e = tl - (ca[m - k - 1] >> 32);
        ca[m - k - 1] = (s << 32) | e;
    }
    if m & 1 != 0 {
        let k = rev_n;
        let s = tl - (ca[k] as u32) as u64;
        let e = tl - (ca[k] >> 32);
        ca[k] = (s << 32) | e;
    }
}

/// `gen_mask_ovlp0` (`Correct.cpp` 21270).
#[allow(clippy::too_many_arguments)]
fn gen_mask_ovlp0(
    graph: &Asg,
    regions: &[OverlapRegion],
    qi: usize,
    ti: usize,
    ov_db: &[UlOv],
    buf: &mut Vec<UlOv>,
    mm: &EmaskTable,
    len_diff: f64,
    rlen: u64,
    ovdb_idx: usize,
    coor_srt: &mut Vec<u64>,
    buf1: &mut Vec<u64>,
    bd: i64,
    res: &mut UlOv,
) -> u32 {
    let bn = buf.len();
    let mut is_exact = 0u64;
    let q = &regions[qi];
    let t = &regions[ti];
    assert!(q.y_identifier < t.y_identifier);
    assert!(ovdb_idx < ov_db.len());
    let mut r0 = UlOv::default();
    cal_x_ul_ovlp(graph, q, t, &mut r0);
    let mut r1 = r0;
    std::mem::swap(&mut r1.qn, &mut r1.tn);
    std::mem::swap(&mut r1.qs, &mut r1.ts);
    std::mem::swap(&mut r1.qe, &mut r1.te);

    let mut dd = 0i64;
    if ov_db[ovdb_idx].qn == q.y_identifier && ov_db[ovdb_idx].tn == t.y_identifier {
        let s = ovdb_idx;
        let e = q.overlap_length as usize + q.x_strand as usize;
        if s as u64 != 0xffff_ffff && e > s {
            push_consist_ul_ovlps(graph, ov_db, s, e, len_diff, &r0, mm, buf, false, &mut is_exact);
        }
        if is_exact == 0 {
            let s = t.overlap_length as usize;
            let e = t.overlap_length as usize + t.x_strand as usize;
            if s as u64 != 0xffff_ffff && e > s {
                push_consist_ul_ovlps(graph, ov_db, s, e, len_diff, &r1, mm, buf, true, &mut is_exact);
            }
        }
    } else {
        assert!(ov_db[ovdb_idx].tn == q.y_identifier && ov_db[ovdb_idx].qn == t.y_identifier);
        let s = ovdb_idx;
        let e = t.overlap_length as usize + t.x_strand as usize;
        if s as u64 != 0xffff_ffff && e > s {
            push_consist_ul_ovlps(graph, ov_db, s, e, len_diff, &r1, mm, buf, true, &mut is_exact);
        }
        if is_exact == 0 {
            let s = q.overlap_length as usize;
            let e = q.overlap_length as usize + q.x_strand as usize;
            if s as u64 != 0xffff_ffff && e > s {
                push_consist_ul_ovlps(graph, ov_db, s, e, len_diff, &r0, mm, buf, false, &mut is_exact);
            }
        }
    }
    *res = UlOv {
        qn: qi as u32,
        tn: ti as u32,
        el: true,
        ..UlOv::default()
    };
    if is_exact == 0 {
        assert!(buf.len() > bn);
        let pieces: Vec<UlOv> = buf[bn..].to_vec();
        dd = cal_paired_distance(graph, q, t, &pieces, coor_srt, buf1, mm, bd);
        if dd != 0 {
            res.el = false;
            if dd > 0 {
                res.qs = dd as u32;
                res.ts = 0;
            } else {
                res.qs = 0;
                res.ts = (-(dd)) as u32;
            }
        }
    }
    buf.truncate(bn);
    let _ = rlen;
    if dd != i32::MIN as i64 {
        1
    } else {
        0
    }
}

/// `if_direct_mask` (`Correct.cpp` 21326).
fn if_direct_mask(
    context: &UlRescallContext<'_>,
    p: &UlOv,
    q: &UlOv,
    bw: u64,
    len_diff: f64,
) -> bool {
    use crate::assembly::ul_rescall::{check_connect_ug, check_connect_rg};
    let (li, lj) = if p.qe > q.qe || (p.qe == q.qe && p.qs >= q.qs) {
        (p, q)
    } else {
        (q, p)
    };
    let os = li.qs.max(lj.qs);
    let oe = li.qe.min(lj.qe);
    let ovlp = if oe > os { i64::from(oe - os) } else { 0 };
    if ovlp == 0 {
        return false;
    }
    if i64::from(lj.qs) <= i64::from(li.qs) + i64::from(crate::assembly::ul_rescall::G_CHAIN_INDEL) {
        if govlp_check(context, li, lj, bw, len_diff) {
            return true;
        }
    } else if i64::from(lj.qe) + i64::from(crate::assembly::ul_rescall::G_CHAIN_INDEL)
        >= i64::from(li.qe)
        && i64::from(lj.qs) + i64::from(crate::assembly::ul_rescall::G_CHAIN_INDEL)
            >= i64::from(li.qs)
    {
        if govlp_check(context, lj, li, bw, len_diff) {
            return true;
        }
    }
    false
}

/// `govlp_check` (`Correct.cpp` 19216) re-exported for this module.
fn govlp_check(context: &UlRescallContext<'_>, li: &UlOv, lj: &UlOv, bw: u64, diff_ec_ul: f64) -> bool {
    use crate::assembly::ul_rescall::{check_connect_ug, check_connect_rg};
    let query_overlap = crate::assembly::ul::infer_rovlp_unitigs(li, lj, &context.unitig_graph.unitigs);
    let li_v = (((li.tn) << 1) | u32::from(li.rev)) ^ 1;
    let lj_v = (((lj.tn) << 1) | u32::from(lj.rev)) ^ 1;
    check_connect_ug(
        &context.unitig_graph.graph,
        li_v,
        lj_v,
        bw as i64,
        diff_ec_ul,
        i64::from(query_overlap),
    ) || check_connect_rg(context, li_v, lj_v, bw as i64, diff_ec_ul, i64::from(query_overlap))
}

/// `mask_ovlps` (`Correct.cpp` 21344).
#[allow(clippy::too_many_arguments)]
fn mask_ovlps(
    context: &UlRescallContext<'_>,
    regions: &mut [OverlapRegion],
    mk: &mut MaskUlOv,
    ov_db: &mut Vec<UlOv>,
    coor_srt: &mut Vec<u64>,
    buf1: &mut Vec<u64>,
    mm: &EmaskTable,
    len_diff: f64,
    rlen: u64,
    bd: i64,
) -> bool {
    mk.srt.clear();
    mk.idx.clear();
    let unitigs = &context.unitig_graph.unitigs;
    let graph = &context.unitig_graph.graph;
    let on = regions.len();
    let mut osrt: Vec<UlOv> = Vec::new();
    for k in 0..on {
        if regions[k].x_strand == 0 {
            continue;
        }
        let tl = unitigs[regions[k].y_identifier as usize].length;
        if let Some(ou) = gen_aln_ul_ov_t(k, tl, &regions[k]) {
            osrt.push(ou);
        }
    }
    if osrt.is_empty() {
        return false;
    }
    let osrt_n = osrt.len();
    radix_sort_by_key_32(osrt.as_mut_slice(), &|entry: &UlOv| entry.qs);
    let mut tot = 0u64;
    let mut tol = 0u64;
    for k in 0..osrt_n {
        tol += u64::from(osrt[k].qe - osrt[k].qs);
        let mut p = convert_ul_ov_t(&regions[osrt[k].qn as usize], unitigs);
        p.qn = osrt[k].qn;
        let mut i = k + 1;
        while i < osrt_n && osrt[i].qs < osrt[k].qe {
            if regions[osrt[i].qn as usize].y_identifier >= regions[osrt[k].qn as usize].y_identifier {
                i += 1;
                continue;
            }
            let mut is_exact = 0u64;
            let mut sid = u64::MAX;
            let mut qv = convert_ul_ov_t(&regions[osrt[i].qn as usize], unitigs);
            qv.qn = osrt[i].qn;
            if if_direct_mask(context, &p, &qv, 16, len_diff) {
                is_exact = 1;
            }
            if is_exact == 0
                && !check_mask_exist(
                    graph,
                    &regions[osrt[i].qn as usize],
                    &regions[osrt[k].qn as usize],
                    ov_db,
                    mm,
                    len_diff,
                    &mut is_exact,
                    &mut sid,
                )
            {
                i += 1;
                continue;
            }
            let mut ou = UlOv::default();
            if is_exact != 0 {
                ou.qs = 0;
                ou.qe = sid as u32;
                ou.ts = 0;
                ou.te = 0;
                ou.el = true;
            } else {
                let w0 = regions[osrt[k].qn as usize].non_homopolymer_errors;
                let w1 = regions[osrt[i].qn as usize].non_homopolymer_errors;
                let mut wt = w0.wrapping_add(w1).wrapping_add(w0.wrapping_sub(w1));
                if wt > u32::MAX {
                    wt = u32::MAX;
                }
                ou.qs = wt;
                ou.qe = sid as u32;
                ou.el = false;
                ou.ts = osrt[i].qe - osrt[i].qs;
                ou.te = osrt[k].qe - osrt[k].qs;
            }
            ou.qn = osrt[i].qn;
            ou.tn = osrt[k].qn;
            tot += u64::from(ou.ts) + u64::from(ou.te);
            osrt.push(ou);
            i += 1;
        }
        let mut ii = k as i64 - 1;
        while ii >= 0 {
            let i = ii as usize;
            if regions[osrt[i].qn as usize].y_identifier >= regions[osrt[k].qn as usize].y_identifier {
                ii -= 1;
                continue;
            }
            let mut is_exact = 0u64;
            let mut sid = u64::MAX;
            let mut qv = convert_ul_ov_t(&regions[osrt[i].qn as usize], unitigs);
            qv.qn = osrt[i].qn;
            if if_direct_mask(context, &p, &qv, 16, len_diff) {
                is_exact = 1;
            }
            if is_exact == 0
                && !check_mask_exist(
                    graph,
                    &regions[osrt[i].qn as usize],
                    &regions[osrt[k].qn as usize],
                    ov_db,
                    mm,
                    len_diff,
                    &mut is_exact,
                    &mut sid,
                )
            {
                ii -= 1;
                continue;
            }
            let mut ou = UlOv::default();
            if is_exact != 0 {
                ou.qs = 0;
                ou.qe = sid as u32;
                ou.ts = 0;
                ou.te = 0;
                ou.el = true;
            } else {
                let w0 = regions[osrt[k].qn as usize].non_homopolymer_errors;
                let w1 = regions[osrt[i].qn as usize].non_homopolymer_errors;
                let mut wt = w0.wrapping_add(w1).wrapping_add(w0.wrapping_sub(w1));
                if wt > u32::MAX {
                    wt = u32::MAX;
                }
                ou.qs = wt;
                ou.qe = sid as u32;
                ou.el = false;
                ou.ts = osrt[i].qe - osrt[i].qs;
                ou.te = osrt[k].qe - osrt[k].qs;
            }
            ou.qn = osrt[i].qn;
            ou.tn = osrt[k].qn;
            tot += u64::from(ou.ts) + u64::from(ou.te);
            osrt.push(ou);
            ii -= 1;
        }
    }

    if tot > rlen * 256 && tot > tol * 32 {
        tot = (rlen * 256).max(tol * 32);
        let osrt_n1 = osrt.len();
        radix_sort_by_key_32(&mut osrt[osrt_n..], &|entry: &UlOv| entry.qs);
        let mut k = osrt_n;
        let mut wt = 0u64;
        while k < osrt_n1 && wt <= tot {
            if osrt[k].el {
                k += 1;
                continue;
            }
            wt += u64::from(osrt[k].ts) + u64::from(osrt[k].te);
            k += 1;
        }
        osrt.truncate(k);
    }

    let osrt_n1 = osrt.len();
    let mut ou_out: Vec<UlOv> = Vec::new();
    for k in osrt_n..osrt_n1 {
        let qi = osrt[k].qn as usize;
        let ti = osrt[k].tn as usize;
        let sid = osrt[k].qe as usize;
        let is_exact_entry = osrt[k].el;
        let mut wt;
        let mut ou = UlOv::default();
        if is_exact_entry {
            ou.qn = qi as u32;
            ou.tn = ti as u32;
            ou.el = true;
            ou.qs = 0;
            ou.ts = 0;
            wt = 1;
        } else {
            wt = gen_mask_ovlp0(
                graph,
                regions,
                qi,
                ti,
                ov_db,
                &mut osrt,
                mm,
                len_diff,
                rlen,
                sid,
                coor_srt,
                buf1,
                bd,
                &mut ou,
            );
        }
        if wt != 0 {
            ou_out.push(ou);
        }
    }
    osrt.clear();
    osrt.extend_from_slice(&ou_out);
    if osrt.is_empty() {
        return false;
    }

    let osrt_n = osrt.len();
    let mirrored: Vec<UlOv> = osrt
        .iter()
        .map(|entry| {
            let mut m = *entry;
            std::mem::swap(&mut m.qn, &mut m.tn);
            std::mem::swap(&mut m.qs, &mut m.ts);
            m
        })
        .collect();
    osrt.extend_from_slice(&mirrored);
    let osrt_n1 = osrt.len();

    radix_sort_by_key_32(osrt.as_mut_slice(), &|entry: &UlOv| entry.qn);
    let mut idx: Vec<u64> = vec![0; regions.len()];
    let mut i = 0usize;
    for k in 1..=osrt_n1 {
        if k == osrt_n1 || osrt[i].qn != osrt[k].qn {
            if k - i > 1 {
                radix_sort_by_key_32(&mut osrt[i..k], &|entry: &UlOv| entry.tn);
            }
            idx[osrt[i].qn as usize] = ((i as u64) << 32) | k as u64;
            i = k;
        }
    }
    mk.srt = osrt;
    mk.idx = idx;
    true
}

/// `push_rphase` (`Correct.cpp` 21676): appends one piece window, merging and
/// overriding shadowed pieces; `is_pri` rides in `cidx` until the final reset.
fn push_rphase(
    region: &mut OverlapRegion,
    zs: u64,
    ze: u64,
    zerr: u64,
    is_pri: u64,
) {
    let windows = &mut region.window_cigar.windows;
    let pe;
    let p_pri;
    if let Some(last) = windows.last() {
        pe = i64::from(last.x_end);
        p_pri = u64::from(last.cigar_start);
    } else {
        pe = i64::from(region.x_position_start);
        p_pri = 0;
    }
    if pe <= zs as i64 {
        if pe < zs as i64 {
            windows.push(WindowList {
                x_start: pe as i32,
                x_end: zs as i32,
                y_start: 0,
                y_end: 0,
                extra_begin: 0,
                extra_end: 0,
                error: 0,
                error_threshold: 0,
                cigar_start: 0,
                cigar_length: 0,
            });
        }
        windows.push(WindowList {
            x_start: zs as i32,
            x_end: ze as i32,
            y_start: 0,
            y_end: 0,
            extra_begin: 0,
            extra_end: 0,
            error: 0,
            error_threshold: 0,
            cigar_start: is_pri as u32,
            cigar_length: zerr as u32,
        });
    } else if pe <= ze as i64 {
        if is_pri != 0 {
            if p_pri != 0 {
                let n = windows.len();
                windows.push(WindowList {
                    x_start: zs as i32,
                    x_end: ze as i32,
                    y_start: 0,
                    y_end: 0,
                    extra_begin: 0,
                    extra_end: 0,
                    error: 0,
                    error_threshold: 0,
                    cigar_start: is_pri as u32,
                    cigar_length: zerr as u32,
                });
                windows[n - 2].x_end = zs as i32;
            } else {
                let last = windows.last_mut().unwrap();
                last.cigar_length = zerr as u32;
                last.cigar_start = is_pri as u32;
                last.x_end = ze as i32;
            }
        } else {
            // If current is not primary: extend the last window.
            let last = windows.last_mut().unwrap();
            last.x_end = ze as i32;
        }
    } else if is_pri != 0 {
        if p_pri == 0 {
            let last = windows.last_mut().unwrap();
            last.cigar_length = zerr as u32;
            last.cigar_start = is_pri as u32;
        } else {
            let n = windows.len();
            let prev_end = windows[n - 2].x_end;
            windows.push(WindowList {
                x_start: zs as i32,
                x_end: prev_end,
                y_start: 0,
                y_end: 0,
                extra_begin: 0,
                extra_end: 0,
                error: 0,
                error_threshold: 0,
                cigar_start: is_pri as u32,
                cigar_length: zerr as u32,
            });
            windows[n - 2].x_end = zs as i32;
        }
    }
}

/// `refine_rphase` (`Correct.cpp` 21648): rebuilds overlap `zid`'s windows
/// from the mask pieces (`a`, the primary reference) and its secondary
/// pieces; `ref` pieces never override an existing primary piece, `qry`
/// pieces do.
pub fn refine_rphase(
    regions: &mut [OverlapRegion],
    slots: &mut SecondaryWindows,
    zid: usize,
    a: &[UlOv],
    buf: &mut Vec<u64>,
) {
    let zid_s = regions[zid].x_position_start;
    let zid_e = regions[zid].x_position_end;
    buf.clear();
    buf.reserve(a.len());
    for k in 0..a.len() {
        if a[k].qs > 0 {
            let tn = a[k].tn as usize;
            let packed = (u64::from(regions[tn].x_position_start) << 32) | k as u64;
            buf.push(packed);
        }
    }
    radix_sort_u64(buf);

    let bn = buf.len();
    let (mut qs, mut qe, mut err) = (0u64, 0u64, 0u64);
    for k in 0..bn {
        if buf.len() > bn {
            qs = buf[buf.len() - 2] >> 32;
            qe = buf[buf.len() - 2] as u32 as u64;
            err = buf[buf.len() - 1];
        }
        let tn = a[buf[k] as u32 as usize].tn as usize;
        let zs = regions[tn].x_position_start.max(zid_s) as u64;
        let ze = regions[tn].x_position_end.min(zid_e) as u64 + 1;
        if ze <= zs {
            continue;
        }
        let zerr = u64::from(a[buf[k] as u32 as usize].qs);
        if qe <= zs {
            buf.push((zs << 32) | ze);
            buf.push(zerr);
        } else {
            if ze > qe {
                qe = ze;
            }
            if zerr > err {
                err = zerr;
            }
            let n = buf.len();
            buf[n - 2] = (qs << 32) | qe;
            buf[n - 1] = err;
        }
    }

    let pushed = buf.len() - bn;
    for k in 0..pushed {
        buf[k] = buf[bn + k];
    }
    buf.truncate(pushed);
    if buf.is_empty() {
        return;
    }
    let ref_entries: Vec<u64> = buf[..].to_vec();
    let ref_n = ref_entries.len();

    let sec_count = regions[zid].alignment_length as usize;
    let mut qry_entries: Vec<u64> = Vec::new();
    for k in 0..sec_count {
        let p = slots[zid][k];
        if p.cigar_length == 0 {
            continue;
        }
        qry_entries.push(((p.x_start as u64) << 32) | p.x_end as u64);
        qry_entries.push(u64::from(p.cigar_length));
    }

    let region = &mut regions[zid];
    region.window_cigar.windows.clear();
    let mut ref_i = 0usize;
    let mut qry_i = 0usize;
    while ref_i < ref_n && qry_i < qry_entries.len() {
        if (ref_entries[ref_i] >> 32) < (qry_entries[qry_i] >> 32) {
            push_rphase(region, ref_entries[ref_i] >> 32, ref_entries[ref_i] as u32 as u64, ref_entries[ref_i + 1], 0);
            ref_i += 2;
        } else {
            push_rphase(region, qry_entries[qry_i] >> 32, qry_entries[qry_i] as u32 as u64, qry_entries[qry_i + 1], 1);
            qry_i += 2;
        }
    }
    while qry_i < qry_entries.len() {
        push_rphase(region, qry_entries[qry_i] >> 32, qry_entries[qry_i] as u32 as u64, qry_entries[qry_i + 1], 1);
        qry_i += 2;
    }
    while ref_i < ref_n {
        push_rphase(region, ref_entries[ref_i] >> 32, ref_entries[ref_i] as u32 as u64, ref_entries[ref_i + 1], 0);
        ref_i += 2;
    }

    let z_pos_s = u64::from(region.x_position_start);
    let z_pos_e1 = u64::from(region.x_position_end) + 1;
    let (mut qs, mut qe) = (z_pos_s, z_pos_s);
    if let Some(last) = region.window_cigar.windows.last() {
        qs = last.x_start as u64;
        qe = last.x_end as u64;
    }
    if qe < z_pos_e1 {
        region.window_cigar.windows.push(WindowList {
            x_start: qe as i32,
            x_end: z_pos_e1 as i32,
            y_start: 0,
            y_end: 0,
            extra_begin: 0,
            extra_end: 0,
            error: 0,
            error_threshold: 0,
            cigar_start: 0,
            cigar_length: 0,
        });
    }

    region.x_strand = u32::MAX;
    region.overlap_length = region.x_position_end + 1 - region.x_position_start;
    region.non_homopolymer_errors = 0;
    let mut zwn = 0u64;
    let mut total = 0u64;
    for k in 0..region.window_cigar.windows.len() {
        let window = &mut region.window_cigar.windows[k];
        if window.cigar_length > 0 {
            region.non_homopolymer_errors += window.cigar_length;
            zwn += (window.x_end - window.x_start) as u64;
        }
        window.cigar_start = 0;
        total += (window.x_end - window.x_start) as u64;
    }
    assert_eq!(total, u64::from(region.overlap_length));
    region.alignment_length = region.overlap_length - zwn as u32;
}


/// `rphase_hl` (`Correct.cpp` 21764): the UL phase pass over one read's
/// candidate overlaps; `len_diff` is the caller's `err` (the UL error rate),
/// `bd` the boundary padding (64), `rlen` the query length.
#[allow(clippy::too_many_arguments)]
pub fn rphase_hl(
    context: &UlRescallContext<'_>,
    regions: &mut [OverlapRegion],
    slots: &mut SecondaryWindows,
    aln: &mut Vec<UlOv>,
    mk: &mut MaskUlOv,
    mm: &EmaskTable,
    len_diff: f64,
    bd: i64,
    rlen: i64,
) {
    let on0 = regions.len();
    let mut active = 0usize;
    {
        use std::sync::atomic::{AtomicU32, Ordering};
        static CALL: AtomicU32 = AtomicU32::new(0);
        let call = CALL.fetch_add(1, Ordering::Relaxed);
        if call == 0 {
            for (k, region) in regions.iter().enumerate() {
                eprintln!(
                    "[RHL] k={} y_id={} x=[{},{}] y=[{},{}] n={}",
                    k,
                    region.y_identifier,
                    region.x_position_start,
                    region.x_position_end,
                    region.y_position_start,
                    region.y_position_end,
                    region.window_cigar.windows.len()
                );
                for (w, m) in region.window_cigar.windows.iter().enumerate() {
                    eprintln!(
                        "[RHLW] k={} w={} x=[{},{}] y=[{},{}] err={} clen={} ualn={} est={} eb={} ee={}",
                        k,
                        w,
                        m.x_start,
                        m.x_end,
                        m.y_start,
                        m.y_end,
                        m.error,
                        m.cigar_length,
                        is_ualn_win(m) as i32,
                        is_est_aln(m) as i32,
                        m.extra_begin,
                        m.extra_end
                    );
                }
            }
        }
    }
    for k in 0..on0 {
        {
            let z = &mut regions[k];
            z.alignment_length = 0;
            z.non_homopolymer_errors = 0;
            z.overlap_length = u32::MAX;
            z.x_strand = 0;
        }
        if regions[k].window_cigar.windows.is_empty() {
            continue;
        }
        if k != active {
            regions.swap(k, active);
            slots.swap(k, active);
        }
        active += 1;
    }

    let mut buf: Vec<u64> = Vec::new();
    let mut idx: Vec<u64> = Vec::with_capacity(active);
    let mut buf1: Vec<u64> = Vec::new();
    let mut coor_srt: Vec<u64> = Vec::new();
    {
        let active_regions = &mut regions[..active];
        let active_slots = &mut slots[..active];
        rphase_rr(active_regions, bd, None, active_slots);
    }
    for k in 0..active {
        let mut sec = 0u64;
        for entry in &slots[k] {
            sec += u64::from(entry.cigar_length);
        }
        regions[k].non_homopolymer_errors = sec as u32;
        let sec_norm = sec - (sec & 3);
        buf.push((sec_norm << 32) | k as u64);
        idx.push((u64::from(regions[k].y_identifier) << 32) | k as u64);
    }
    radix_sort_u64(&mut idx);
    radix_sort_u64(&mut buf);

    let mut l = 0usize;
    for k in 1..=buf.len() {
        if k == buf.len() || (buf[l] >> 32) == (buf[k] >> 32) {
            if k - l > 1 && k < buf.len() {
                for i in l..k {
                    let oid = buf[i] as u32 as usize;
                    let m0 = u64::from(regions[oid].x_position_end + 1 - regions[oid].x_position_start);
                    let m1 = u64::from(regions[oid].y_position_end + 1 - regions[oid].y_position_start);
                    let sec = m0.min(m1);
                    let sec = u32::MAX - sec as u32;
                    buf[i] = ((sec as u64) << 32) | (buf[i] as u32 as u64);
                }
                radix_sort_u64(&mut buf[l..k]);
            }
            l = k;
        }
    }

    mk.idx.clear();
    mk.srt.clear();
    aln.clear();
    let mut mmov = 1_000_000usize;
    if mmov > active * 32 {
        mmov = active * 32;
    }
    for k in 0..buf.len() {
        if aln.len() > mmov {
            break;
        }
        let oid = buf[k] as u32 as usize;
        regions[oid].overlap_length = aln.len() as u32;
        let m0 = push_emask_flt(&mm[regions[oid].y_identifier as usize], &idx, u64::from(regions[oid].y_identifier), aln);
        let m1 = gen_src_shared_interval_simple(
            regions[oid].y_identifier,
            &context.unitig_graph.graph,
            &idx,
            aln,
        );
        dedup_src_shared1(&context.unitig_graph.graph, aln, m0, m1, mm);
        regions[oid].x_strand = (aln.len() as u32).wrapping_sub(regions[oid].overlap_length);
        if regions[oid].x_strand == 0 {
            regions[oid].overlap_length = u32::MAX;
        }
    }

    let mut masked = false;
    if !aln.is_empty() || active > 1 {
        masked = mask_ovlps(
            context,
            &mut regions[..active],
            mk,
            aln,
            &mut coor_srt,
            &mut buf1,
            mm,
            len_diff,
            rlen as u64,
            bd,
        );
    }
    if masked {
        {
            let active_regions = &mut regions[..active];
            let active_slots = &mut slots[..active];
            rphase_rr(active_regions, bd, Some(mk), active_slots);
        }
        for k in 0..active {
            let entry = mk.idx[k];
            if (entry as u32 as u64) > (entry >> 32) {
                let begin = (entry >> 32) as usize;
                let end = entry as u32 as usize;
                let masks: Vec<UlOv> = mk.srt[begin..end].to_vec();
                refine_rphase(regions, slots, k, &masks, &mut coor_srt);
            }
        }
    }

    for k in 0..on0 {
        let z = &mut regions[k];
        if z.x_strand == u32::MAX {
            z.x_strand = 0;
            continue;
        }
        z.overlap_length = z.x_position_end + 1 - z.x_position_start;
        z.non_homopolymer_errors = 0;
        let mut zwn = 0u64;
        let windows = &mut z.window_cigar.windows;
        for i in 0..slots[k].len() {
            if i < windows.len() {
                windows[i] = slots[k][i];
            } else {
                windows.push(slots[k][i]);
            }
            if windows[i].cigar_length > 0 {
                z.non_homopolymer_errors += windows[i].cigar_length;
                zwn += (windows[i].x_end - windows[i].x_start) as u64;
            }
        }
        windows.truncate(slots[k].len());
        z.x_strand = 0;
        eprintln!(
            "[RPH] k={} y_id={} ol={} zwn={} al={} n={}",
            k,
            z.y_identifier,
            z.overlap_length,
            zwn,
            slots[k].len(),
            windows.len()
        );
        assert!(zwn <= u64::from(z.overlap_length));
        z.alignment_length = z.overlap_length - zwn as u32;
    }
}

/// `gen_src_shared_interval_simple` (`inter.cpp` 11565): the graph-arc
/// intervals of `src`'s neighbourhood, kept only when the target appears in
/// the sorted filter.
fn gen_src_shared_interval_simple(
    src: u32,
    graph: &Asg,
    flt: &[u64],
    res: &mut Vec<UlOv>,
) -> usize {
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
        let mut fi = 0usize;
        let mut i = rn;
        for k in rn + 1..=res.len() {
            if k == res.len() || res[k].tn != res[i].tn {
                while fi < flt.len() && (flt[fi] >> 32) < u64::from(res[i].tn) {
                    fi += 1;
                }
                if fi < flt.len() && (flt[fi] >> 32) == u64::from(res[i].tn) {
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
                    }
                } else {
                    for v in i..k {
                        res[v].tn = u32::MAX;
                    }
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
            res[z].sec = SEC_MODE;
            z += 1;
        }
        res.truncate(z);
    }
    res.len() - rn
}

/// `gen_commen_sketch` (`Correct.cpp` 22282): walks `[qs, qe)` window by
/// window, keeps the windows where every candidate overlap has an identical
/// target piece, and chains the kept starts with `dp_commen_sketch`.
#[allow(clippy::too_many_arguments)]
fn gen_commen_sketch(
    targets: &TargetSource<'_>,
    regions: &mut [OverlapRegion],
    id_a: &mut [u64],
    s: u64,
    e: u64,
    ql: u64,
    wl: u64,
    buf: &mut [u64],
    dp: u64,
    str0: &mut [u8],
    str1: &mut [u8],
    aln: &mut [UlOv],
    trace: &mut Vec<u64>,
    max_lgap: i64,
    sgap_rate: f64,
    sgap: i64,
) -> usize {
    let mut id_n = id_a.len();
    if id_n == 0 {
        return id_n;
    }
    let mut rm_n = 0usize;
    let mut qs = (s / wl) * wl;
    if qs < s {
        qs += wl;
    }
    if qs >= ql {
        return id_n;
    }
    let mut qe = (e / wl) * wl;
    if qe >= ql {
        qe = ql;
    }
    if qs >= qe {
        return id_n;
    }
    let mut buf_n = 0usize;
    for k in 0..id_n {
        let entry = &aln[id_a[k] as usize];
        if u64::from(entry.qs) <= qs && u64::from(entry.qe) >= qe {
            buf[buf_n] = u64::from(regions[aln[id_a[k] as usize].qn as usize].alignment_length);
            buf[buf_n] <<= 32;
            buf[buf_n] |= id_a[k];
            buf_n += 1;
        }
        if u64::from(entry.qe) < e {
            rm_n += 1;
        }
    }
    assert!(buf_n == dp as usize && buf_n > 1);

    if buf_n > 0 {
        let mut co = 1u64;
        trace.clear();
        let mut occ = 0u64;
        let mut k = qs;
        while k < qe {
            let ws = k;
            let mut we = ws + wl;
            if we > qe {
                we = qe;
            }
            let first = aln[id_a[buf[0] as u32 as usize] as usize].qn as usize;
            let (wid, _) = window_id_by_start(
                regions[first].x_position_start as i32,
                regions[first].x_position_end as i32,
                ws as i64,
                wl as i64,
            );
            let Some((r_y0, r_y1, r_err)) =
                get_win_aln(&regions[first], wid as usize)
            else {
                k += wl;
                continue;
            };
            if !ovlp_win_check(
                &regions[first],
                regions[first].alignment_length,
                wid as u32,
                max_lgap,
                sgap_rate,
                sgap,
            ) {
                k += wl;
                continue;
            }
            let rxl = i64::from(regions[first].window_cigar.windows[wid as usize].x_end) + 1
                - i64::from(regions[first].window_cigar.windows[wid as usize].x_start);
            let ryl = r_y1 - r_y0;
            if rxl > 1 && ryl > 1 {
                k += wl;
                continue;
            }
            let mut i = 1usize;
            while i < buf_n {
                let other = aln[id_a[buf[i] as u32 as usize] as usize].qn as usize;
                let (wid2, _) = window_id_by_start(
                    regions[other].x_position_start as i32,
                    regions[other].x_position_end as i32,
                    ws as i64,
                    wl as i64,
                );
                let Some((p_y0, p_y1, p_err)) = get_win_aln(&regions[other], wid2 as usize) else {
                    break;
                };
                if !ovlp_win_check(
                    &regions[other],
                    regions[other].alignment_length,
                    wid2 as u32,
                    max_lgap,
                    sgap_rate,
                    sgap,
                ) {
                    break;
                }
                let pxl = i64::from(regions[other].window_cigar.windows[wid2 as usize].x_end) + 1
                    - i64::from(regions[other].window_cigar.windows[wid2 as usize].x_start);
                if pxl != rxl {
                    break;
                }
                if (r_y1 - r_y0) != (p_y1 - p_y0) || r_err != p_err {
                    break;
                }
                if r_err == 0 {
                    i += 1;
                    continue;
                }
                let qname = targets.sub_region(
                    regions[first].y_identifier as usize,
                    r_y0,
                    r_y1 - r_y0,
                    regions[first].y_strand as u8,
                );
                str0[..qname.len()].copy_from_slice(&qname);
                let tname = targets.sub_region(
                    regions[other].y_identifier as usize,
                    p_y0,
                    p_y1 - p_y0,
                    regions[other].y_strand as u8,
                );
                str1[..tname.len()].copy_from_slice(&tname);
                if str0[..(p_y1 - p_y0) as usize] != str1[..(p_y1 - p_y0) as usize] {
                    break;
                }
                i += 1;
            }
            if i < buf_n {
                k += wl;
                continue;
            }
            if co != 0 {
                for i in 0..buf_n {
                    let other = aln[id_a[buf[i] as u32 as usize] as usize].qn as usize;
                    let (wid3, _) = window_id_by_start(
                        regions[other].x_position_start as i32,
                        regions[other].x_position_end as i32,
                        ws as i64,
                        wl as i64,
                    );
                    co = u64::from(ovlp_win_check(
                        &regions[other],
                        (buf[i] >> 32) as u32,
                        wid3 as u32,
                        max_lgap,
                        sgap_rate,
                        sgap,
                    ));
                    let m = ((wid3 as u64) << 32) | (buf[i] as u32 as u64);
                    buf[i] = m;
                }
                if co != 0 {
                    occ += 1;
                }
            }
            trace.push(ws);
            k += wl;
        }

        if !trace.is_empty() {
            if co == 0 {
                trace.reserve(trace.len());
            }
            let dp_needed = trace.len();
            let mut dp_buf = vec![0u64; dp_needed];
            let kept = {
                let trace_slice = trace.as_mut_slice();
                dp_commen_sketch(
                    aln,
                    regions,
                    &id_a[..buf_n],
                    buf_n as i64,
                    trace_slice,
                    dp_needed as i64,
                    &mut dp_buf,
                    occ as i64,
                    wl as i64,
                    dp as i64,
                    max_lgap,
                    sgap_rate,
                    sgap,
                )
            };
            trace.truncate(kept as usize);
        }
    }
    if rm_n > 0 {
        let mut m = 0usize;
        for i in 0..id_n {
            if u64::from(aln[id_a[i] as usize].qe) < e {
                continue;
            }
            id_a[m] = id_a[i];
            m += 1;
        }
        id_n = m;
    }
    id_n
}

/// `update_sketch_trace` (`Correct.cpp` 22372): marks the co-linear shared
/// windows of `[qs, qe)` stretches covered by at least two overlaps.
#[allow(clippy::too_many_arguments)]
pub fn update_sketch_trace(
    targets: &TargetSource<'_>,
    regions: &mut [OverlapRegion],
    aln: &mut [UlOv],
    ql: i64,
    wl: i64,
    idx: &mut Vec<u64>,
    b0: &mut Vec<u64>,
    b1: &mut Vec<u64>,
    max_lgap: i64,
    sgap_rate: f64,
) {
    idx.clear();
    if aln.is_empty() {
        return;
    }
    for z in regions.iter_mut() {
        append_unmatched_windows(z.x_position_start as i32, z.x_position_end as i32, &mut z.window_cigar, wl);
        for window in z.window_cigar.windows.iter_mut() {
            if window.extra_end < 0 {
                window.extra_end = 0;
            }
        }
        z.alignment_length = u32::MAX;
    }

    idx.reserve(aln.len() << 1);
    b0.reserve(aln.len());
    let mut srt_n = 0usize;
    for (i, entry) in aln.iter_mut().enumerate() {
        regions[entry.qn as usize].alignment_length = u32::MAX;
        idx.push((((u64::from(entry.qs)) << 1) << 32) + i as u64);
        srt_n += 1;
        idx.push(((((u64::from(entry.qe) - 1) << 1) + 1) << 32) + i as u64);
        srt_n += 1;
        entry.el = true;
    }

    radix_sort_u64(&mut idx[..srt_n]);
    idx.truncate(srt_n);
    let mut str0 = vec![0u8; (wl << 1) as usize];
    let mut str1 = vec![0u8; (wl << 1) as usize];
    let mut dp = 0i64;
    let mut beg: i64 = 0;
    let mut end: i64 = -1;
    for i in 0..srt_n {
        let old_dp = dp;
        if (idx[i] >> 32) & 1 != 0 {
            dp -= 1;
            end = ((idx[i] >> 33) + 1) as i64;
        } else {
            dp += 1;
            end = (idx[i] >> 33) as i64;
            idx.push(idx[i] & 0xffff_ffff);
        }
        if end > beg && end - beg > wl && old_dp >= 2 {
            if b0.len() < old_dp as usize {
                b0.resize(old_dp as usize, 0);
            }
            let survivors = {
                let dp_buf = &mut b0[..old_dp as usize];
                gen_commen_sketch(
                    targets,
                    regions,
                    &mut idx[srt_n..],
                    beg as u64,
                    end as u64,
                    ql as u64,
                    wl as u64,
                    dp_buf,
                    old_dp as u64,
                    &mut str0,
                    &mut str1,
                    aln,
                    b1,
                    max_lgap,
                    sgap_rate,
                    SGAP,
                )
            };
            idx.truncate(srt_n + survivors);
        }
        beg = end;
    }

    for i in 0..aln.len() {
        if i == 0 || aln[i].qn != aln[i - 1].qn {
            regions[aln[i].qn as usize].alignment_length = i as u32;
        }
    }
    idx.truncate(srt_n);
}

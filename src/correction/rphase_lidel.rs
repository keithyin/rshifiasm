//! The ONT large-indel phasing of `rphase_hc`, the counterparts of
//! `rphase_lidel`, `rphase_lidel_cc`, `rcall_lidel_variant`,
//! `gen_ov_lidel_variant`, `push_idel_info`, `extract_sub_cigar_sv`,
//! `iter_sub_cigar_sv`, `extract_sub_err`, `cal_lindel_dd`, `set_cgid` and
//! `is_get_group` in `Correct.cpp`.
//!
//! The reference repurposes `ul_ov_t` twice: as interval variants (`qs/qe`
//! query span, `ts/te` target span, `qn` error, `tn` overlap id, `sec` cluster
//! degree then id, `el` homopolymer flag) and as window cursors (`tn` overlap
//! id, `qn` window id, `qs/ts` x/y offsets, `qe` cigar offset, `sec` boundary).
//! [`UlOv`] is reused for both, with the meaning fixed by the function.

use crate::assembly::ul::UlOv;
use crate::container::ksort::{radix_sort_by_key_32, radix_sort_u64};
use crate::correction::haplotype_evidence::{HaplotypeEvidence, HaplotypeEvidenceAlloc, SnpStats};
use crate::correction::rphase::hpc_mask_ff_region;
use crate::overlap::region::OverlapRegion;
use crate::reads::read_store::ReadStore;

/// `cal_lindel_dd` (`Correct.cpp` 19740).
fn cal_lindel_dd(a: &UlOv, b: &UlOv, ol_r: f64, ol_w: u64, err_dif: f64) -> f64 {
    if a.tn == b.tn {
        return -1.0;
    }
    let os = a.qs.max(b.qs) as i64;
    let oe = a.qe.min(b.qe) as i64;
    if oe <= os + ol_w as i64 {
        return -1.0;
    }
    let ol = (oe - os) as u64;
    if (ol as f64) < ((a.qe - a.qs) as f64 * ol_r)
        && (ol as f64) < ((b.qe - b.qs) as f64 * ol_r)
    {
        return -1.0;
    }
    let ea = a.qn as i64;
    let eb = b.qn as i64;
    let ed = (ea - eb).abs();
    if (ed as f64) > (ea as f64 * err_dif) || (ed as f64) > (eb as f64 * err_dif) {
        return -1.0;
    }
    let mut sc = ((ea - ed + eb - ed) as f64) / ((ea + eb) as f64);
    sc += ((ol + ol) as f64) / ((a.qe - a.qs + b.qe - b.qs) as f64);
    sc
}

/// `set_cgid` (`Correct.cpp` 19760): assigns the cluster id, one member per
/// overlap.
fn set_cgid(
    z: &mut UlOv,
    ia: &mut Vec<u64>,
    ian: usize,
    ca: &mut [u64],
    ci: u64,
    gni: &mut i64,
) -> bool {
    assert!(ca[z.tn as usize] <= 1);
    if z.ts != u32::MAX {
        return false;
    }
    if ca[z.tn as usize] != 0 {
        return false;
    }
    z.ts = ci as u32;
    ca[z.tn as usize] += 1;
    ia.push(z.tn as u64);
    assert!(ia.len() <= ian);
    *gni += 1;
    true
}

/// `is_get_group` (`Correct.cpp` 19779).
fn is_get_group(a: &[UlOv], ga: &[u64], ia: &[u64], gi: u64, tn: u32) -> bool {
    let mut k = ga[gi as usize] >> 32;
    while k != u64::MAX && a[k as usize].tn != tn {
        k = ia[k as usize];
    }
    !(k != u64::MAX && a[k as usize].tn == tn)
}

/// `rphase_lidel_cc` (`Correct.cpp` 19796): clusters the indel variants.
#[allow(clippy::too_many_arguments)]
fn rphase_lidel_cc(
    a: &mut Vec<UlOv>,
    overlap_count: usize,
    len_st: f64,
    len_w: u64,
    err_dif: f64,
    c_sz: u64,
    buf: &mut Vec<u64>,
    idx: &mut Vec<u64>,
) -> i64 {
    let mut an = a.len();
    let ian = overlap_count;
    // `idx` is the reference's `ia`/`ca` scratch; the cluster pass below keeps
    // its own vectors with the same semantics and `rcall` clears `idx` again.
    idx.resize(ian * 2, 0);

    for k in 0..an {
        let mut ol = ((a[k].qe - a[k].qs) as f64 * len_st) as u64;
        if ol < len_w {
            ol = len_w;
        }
        if (a[k].qe - a[k].qs) < ol as u32 {
            continue;
        }
        let mut z = k + 1;
        while z < an && a[z].qs < a[k].qe {
            if cal_lindel_dd(&a[k], &a[z], len_st, len_w, err_dif) >= 0.0 {
                a[k].sec += 1;
                a[z].sec += 1;
            }
            z += 1;
        }
    }

    let mut write = 0usize;
    for k in 0..an {
        if a[k].sec == 0 {
            continue;
        }
        a[k].ts = u32::MAX;
        a[k].te = u32::MAX;
        a[write] = a[k];
        write += 1;
    }
    an = write;
    a.truncate(an);
    if an == 0 {
        return 0;
    }

    let mut ck = 0u64;
    let mut gni = 0i64;
    let mut ia: Vec<u64> = Vec::new();
    let mut ca = vec![0u64; ian];
    for k in 0..an {
        if (a[k].sec as u64) < c_sz || a[k].ts != u32::MAX {
            continue;
        }
        buf.clear();
        buf.push(k as u64);
        ia.clear();
        while let Some(v) = buf.pop() {
            let v = v as usize;
            let mut ol = ((a[v].qe - a[v].qs) as f64 * len_st) as u64;
            if ol < len_w {
                ol = len_w;
            }
            if (a[v].qe - a[v].qs) < ol as u32 {
                continue;
            }
            let assigned = set_cgid(&mut a[v], &mut ia, ian, &mut ca, ck, &mut gni);
            if assigned || a[v].ts == ck as u32 {
                let mut z = 0usize;
                while z < an && a[z].qs < a[v].qe {
                    if (a[z].sec as u64) < c_sz || a[z].ts != u32::MAX || v == z {
                        z += 1;
                        continue;
                    }
                    if cal_lindel_dd(&a[v], &a[z], len_st, len_w, err_dif) >= 0.0
                        && set_cgid(&mut a[z], &mut ia, ian, &mut ca, ck, &mut gni)
                    {
                        buf.push(z as u64);
                    }
                    z += 1;
                }
            }
        }
        ck += 1;
    }

    if ck == 0 || gni <= 0 {
        return 0;
    }

    if gni < an as i64 {
        // Build the per-cluster linked index, then attach the singletons.
        buf.resize(an + ck as usize, u64::MAX);
        let (ca2, ia2) = buf.split_at_mut(ck as usize);
        ia2.fill(u64::MAX);
        for k in 0..an {
            if a[k].ts == u32::MAX {
                continue;
            }
            let g = a[k].ts as usize;
            if ca2[g] != u64::MAX {
                let prev = ca2[g] as u32 as usize;
                assert!(ia2[prev] == u64::MAX);
                ia2[prev] = k as u64;
                ca2[g] = (ca2[g] >> 32 << 32) | k as u64;
            } else {
                ca2[g] = ((k as u64) << 32) | k as u64;
            }
        }
        let mut ck2 = 0u64;
        for k in 0..an {
            if a[k].ts != u32::MAX {
                ck2 += 1;
                continue;
            }
            let mut ol = ((a[k].qe - a[k].qs) as f64 * len_st) as u64;
            if ol < len_w {
                ol = len_w;
            }
            if (a[k].qe - a[k].qs) < ol as u32 {
                continue;
            }
            let mut mk = -1i64;
            let mut msw = -1.0f64;
            let mut z = 0usize;
            while z < an && a[z].qs < a[k].qe {
                if (a[z].sec as u64) < c_sz || a[z].ts == u32::MAX || z == k {
                    z += 1;
                    continue;
                }
                let sw = cal_lindel_dd(&a[k], &a[z], len_st, len_w, err_dif);
                if sw < 0.0 {
                    z += 1;
                    continue;
                }
                if sw > msw && is_get_group(a, ca2, ia2, a[z].ts as u64, a[k].tn) {
                    mk = a[z].ts as i64;
                    msw = sw;
                }
                z += 1;
            }
            if mk != -1 {
                ck2 += 1;
                a[k].ts = mk as u32;
                let g = a[k].ts as usize;
                let prev = ca2[g] as u32 as usize;
                assert!(ia2[prev] == u64::MAX);
                ia2[prev] = k as u64;
                ca2[g] = (ca2[g] >> 32 << 32) | k as u64;
            }
        }
        ck = ck2;
    } else {
        ck = gni as u64;
    }
    radix_sort_by_key_32(a.as_mut_slice(), &|v: &UlOv| v.ts);
    ck as i64
}

/// The cells of one window, `set_bit_extz_t(ez, z, wk)`.
fn window_cells<'a>(z: &'a OverlapRegion, wk: usize) -> &'a [u16] {
    let window = &z.window_cigar.windows[wk];
    let start = window.cigar_start as usize;
    &z.window_cigar.cells[start..start + window.cigar_length as usize]
}

/// `extract_sub_err` (`Correct.cpp` 18674).
#[allow(clippy::too_many_arguments)]
fn extract_sub_err(
    z: &OverlapRegion,
    s: i64,
    e: i64,
    os0: i64,
    oe0: i64,
    err0: i64,
    err_sec_rate: f64,
    p: &mut UlOv,
) -> i64 {
    let wk = p.qn as usize;
    let mut xk = p.qs as i64;
    let mut yk = p.ts as i64;
    let mut ck = p.qe as i64;
    let bd = p.sec as i64;
    let window = &z.window_cigar.windows[wk];
    let cells = window_cells(z, wk);
    let cn = cells.len() as i64;
    let s0 = window.x_start as i64 + bd;
    let e0 = window.x_end as i64 + 1 - bd;
    let mut s = s;
    let mut e = e;
    if s < s0 {
        s = s0;
    }
    if e > e0 {
        e = e0;
    }
    if s >= e {
        return -1;
    }
    let os = s.max(s0);
    let oe = e.min(e0);
    if oe <= os {
        return -1;
    }
    if cn == 0 {
        return -1;
    }
    if ck < 0 || ck > cn {
        ck = 0;
        xk = window.x_start as i64;
        yk = window.y_start as i64;
    }

    while ck > 0 && xk >= s {
        ck -= 1;
        let cell = cells[ck as usize];
        let op = cell >> 14;
        let len = (cell & 0x3fff) as i64;
        if op != 2 {
            xk -= len;
        }
        if op != 3 {
            yk -= len;
        }
    }

    let mut l = [0i64; 2];
    while ck < cn && xk < e {
        let ws = xk;
        let cell = cells[ck as usize];
        let op = cell >> 14;
        let ol = (cell & 0x3fff) as i64;
        if op != 2 {
            xk += ol;
        }
        if op != 3 {
            yk += ol;
        }
        ck += 1;
        let we = xk;
        let mut t = false;
        let os = s.max(ws);
        let oe = e.min(we);
        let ovlp = if oe > os { oe - os } else { 0 };
        if op == 2 && ws == we && ws >= s && ws <= e {
            t = true;
        }
        if ovlp > 0 {
            t = true;
        }
        if !t {
            continue;
        }
        if op != 2 {
            l[usize::from(op != 0)] += ovlp;
        } else {
            l[1] += ol;
        }
    }

    while ck < cn && xk <= e {
        let ws = xk;
        let cell = cells[ck as usize];
        let op = cell >> 14;
        let mut ol = (cell & 0x3fff) as i64;
        if op != 2 {
            break;
        }
        ck += 1;
        while ck < cn && op == cells[ck as usize] >> 14 {
            ol += (cells[ck as usize] & 0x3fff) as i64;
            ck += 1;
        }
        yk += ol;
        if ws >= s && ws <= e {
            l[1] += ol;
        }
    }

    p.qs = xk as u32;
    p.ts = yk as u32;
    p.qe = ck as u32;
    p.te = 0;

    let mut err = err0 as f64;
    if (oe0 - os0) < (e - s) {
        err = ((oe0 - os0) as f64 / (e - s) as f64) * err0 as f64;
    }
    let err = (err * err_sec_rate) as i64;
    if l[1] <= err {
        l[1]
    } else {
        -1
    }
}

/// The walk cursor of `iter_sub_cigar_sv`.
struct SvCursor {
    ck: i64,
    xk: i64,
    yk: i64,
}

/// `iter_sub_cigar_sv` (`Correct.cpp` 18760).
#[allow(clippy::too_many_arguments)]
fn iter_sub_cigar_sv(
    zid: u32,
    cp: &mut UlOv,
    cells: &[u16],
    cur: &mut SvCursor,
    qstr: Option<&[u8]>,
    tstr: Option<&[u8]>,
    ts: i64,
    te: i64,
    hpc_len: i64,
    hpc_rr: i64,
    exd_err_bd: i64,
    ql: i64,
) -> (i64, i64, i64) {
    let cn = cells.len() as i64;
    let mut err = 0i64;
    let mut herr = 0i64;
    let mut ex = [-1i64; 2];
    let mut ey = [-1i64; 2];
    let s = cp.qs as i64;
    let e = cp.qe as i64;

    while cur.ck > 0 && cur.xk >= s {
        cur.ck -= 1;
        let cell = cells[cur.ck as usize];
        let op = cell >> 14;
        let len = (cell & 0x3fff) as i64;
        if op != 2 {
            cur.xk -= len;
        }
        if op != 3 {
            cur.yk -= len;
        }
    }

    while cur.ck < cn && cur.xk < e {
        let wx0 = cur.xk;
        let wy0 = cur.yk;
        let cell = cells[cur.ck as usize];
        let op = cell >> 14;
        let cl = (cell & 0x3fff) as i64;
        cur.ck += 1;
        if op != 2 {
            cur.xk += cl;
        }
        if op != 3 {
            cur.yk += cl;
        }
        let wx1 = cur.xk;
        let wy1 = cur.yk;
        if op == 0 {
            continue;
        }
        if wx0 >= s && wx1 <= e {
            if ex[0] < 0 {
                ex[0] = wx0;
            }
            if ey[0] < 0 {
                ey[0] = wy0;
            }
            ex[1] = wx1;
            ey[1] = wy1;
            err += cl;
            if hpc_len > 0 {
                assert!(wy0 >= ts && wy1 <= te);
                let q = qstr.unwrap();
                let t = tstr.unwrap();
                if hpc_mask_ff_region(q, wx0, wx1, hpc_len, hpc_rr, 6, 0.51)
                    || hpc_mask_ff_region(t, wy0 - ts, wy1 - ts, hpc_len, hpc_rr, 6, 0.51)
                {
                    herr += cl;
                }
            }
        }
    }

    while cur.ck < cn && cur.xk <= e {
        let wx0 = cur.xk;
        let wy0 = cur.yk;
        let cell = cells[cur.ck as usize];
        let op = cell >> 14;
        let cl = (cell & 0x3fff) as i64;
        if op != 2 {
            break;
        }
        cur.yk += cl;
        cur.ck += 1;
        let wx1 = cur.xk;
        let wy1 = cur.yk;
        if wx0 >= s && wx1 <= e {
            if ex[0] < 0 {
                ex[0] = wx0;
            }
            if ey[0] < 0 {
                ey[0] = wy0;
            }
            ex[1] = wx1;
            ey[1] = wy1;
            err += cl;
            if hpc_len > 0 {
                assert!(wy0 >= ts && wy1 <= te);
                let q = qstr.unwrap();
                let t = tstr.unwrap();
                if hpc_mask_ff_region(q, wx0, wx1, hpc_len, hpc_rr, 6, 0.51)
                    || hpc_mask_ff_region(t, wy0 - ts, wy1 - ts, hpc_len, hpc_rr, 6, 0.51)
                {
                    herr += cl;
                }
            }
        }
    }

    if hpc_len > 0 {
        assert!(cp.qn as i64 == err);
        cp.qs = s as u32;
        cp.qe = e as u32;
        cp.ts = ex[0] as u32;
        cp.te = ex[1] as u32;
        cp.tn = zid;
        cp.qn = err as u32;
        cp.el = !(herr > 0 && (herr as f64) > (err as f64 * 0.66));
        if err <= exd_err_bd {
            let mut p = (s + e) / 2 - exd_err_bd;
            if p < 0 {
                p = 0;
            }
            if p < cp.qs as i64 {
                cp.qs = p as u32;
            }
            let mut p = (s + e) / 2 + exd_err_bd;
            if p > ql {
                p = ql;
            }
            if p > cp.qe as i64 {
                cp.qe = p as u32;
            }
        }
    }
    (err, ey[0], ey[1])
}

/// `extract_sub_cigar_sv` (`Correct.cpp` 18839).
#[allow(clippy::too_many_arguments)]
fn extract_sub_cigar_sv(
    z: &OverlapRegion,
    zid: usize,
    wk: usize,
    reads: &ReadStore,
    qstr: &[u8],
    rr: &mut Vec<UlOv>,
    min_err: u64,
    hpc_len: i64,
    hpc_rr: i64,
) -> i64 {
    let window = &z.window_cigar.windows[wk];
    let mut xk = window.x_start as i64;
    let mut yk = window.y_start as i64;
    let cells = window_cells(z, wk);
    let cn = cells.len() as i64;
    if cn == 0 {
        return -1;
    }
    let ql = qstr.len() as i64;
    let tl = reads.read_length(z.y_identifier as usize) as i64;
    let rr_n0 = rr.len();
    let mut ck = 0i64;
    let mut ex = [-1i64; 2];
    let mut ey = [-1i64; 2];
    let mut is_srt = true;

    let push_interval = |rr: &mut Vec<UlOv>, ex: [i64; 2], ey: [i64; 2], is_srt: &mut bool| {
        let cp = UlOv {
            qs: ex[0] as u32,
            qe: ex[1] as u32,
            ts: ey[0] as u32,
            te: ey[1] as u32,
            sec: 0,
            el: false,
            rev: false,
            qn: 0,
            tn: 0,
        };
        if rr.len() > rr_n0 + 1 && cp.qs < rr[rr.len() - 1].qs {
            *is_srt = false;
        }
        rr.push(cp);
    };

    while ck < cn {
        let wx0 = xk;
        let wy0 = yk;
        let cell = cells[ck as usize];
        let mut op = cell >> 14;
        let mut cl = (cell & 0x3fff) as i64;
        if op != 2 {
            xk += cl;
        }
        if op != 3 {
            yk += cl;
        }
        let op0 = op;
        let mut el0 = cl;
        ck += 1;
        while ck < cn {
            let cell = cells[ck as usize];
            op = cell >> 14;
            if (op != 0) != (op0 != 0) {
                break;
            }
            cl = (cell & 0x3fff) as i64;
            if op != 2 {
                xk += cl;
            }
            if op != 3 {
                yk += cl;
            }
            el0 += cl;
            ck += 1;
        }
        let wx1 = xk;
        let wy1 = yk;
        if op0 == 0 {
            continue;
        }
        let mut ff = true;
        let ox0 = if wx0 >= el0 { wx0 - el0 } else { 0 };
        let ox1 = if wx1 + el0 <= ql { wx1 + el0 } else { ql };
        let oy0 = if wy0 >= el0 { wy0 - el0 } else { 0 };
        let oy1 = if wy1 + el0 <= tl { wy1 + el0 } else { tl };
        if ex[0] < 0 || ex[1] < 0 || ey[0] < 0 || ey[1] < 0 {
            ff = false;
        }
        if ff {
            if ox1.min(ex[1]) <= ox0.max(ex[0]) {
                ff = false;
            }
            if oy1.min(ey[1]) <= oy0.max(ey[0]) {
                ff = false;
            }
        }
        if ff {
            ex[0] = ox0.min(ex[0]);
            ex[1] = ox1.max(ex[1]);
            ey[0] = oy0.min(ey[0]);
            ey[1] = oy1.max(ey[1]);
        } else {
            if ex[0] >= 0 && ex[1] >= 0 && ey[0] >= 0 && ey[1] >= 0 {
                push_interval(rr, ex, ey, &mut is_srt);
            }
            ex = [ox0, ox1];
            ey = [oy0, oy1];
        }
    }
    if ex[0] >= 0 && ex[1] >= 0 && ey[0] >= 0 && ey[1] >= 0 {
        push_interval(rr, ex, ey, &mut is_srt);
    }
    if rr.len() <= rr_n0 {
        return 0;
    }
    if !is_srt {
        radix_sort_by_key_32(&mut rr[rr_n0..], &|v: &UlOv| v.qs);
    }

    // Merge overlapping intervals and drop the short ones.
    let rn = rr.len() - rr_n0;
    let mut ck2 = 0usize;
    for k in 0..rn {
        let entry = rr[rr_n0 + k];
        if ck2 > 0 && entry.qs < rr[rr_n0 + ck2 - 1].qe {
            let prev = &mut rr[rr_n0 + ck2 - 1];
            if entry.qe > prev.qe {
                prev.qe = entry.qe;
            }
            if entry.ts < prev.ts {
                prev.ts = entry.ts;
            }
            if entry.te > prev.te {
                prev.te = entry.te;
            }
        } else {
            if ck2 > 0 {
                let prev = rr[rr_n0 + ck2 - 1];
                if (prev.qe - prev.qs) < min_err as u32 && (prev.te - prev.ts) < min_err as u32 {
                    ck2 -= 1;
                }
            }
            rr[rr_n0 + ck2] = entry;
            ck2 += 1;
        }
    }
    if ck2 > 0 {
        let prev = rr[rr_n0 + ck2 - 1];
        if (prev.qe - prev.qs) < min_err as u32 && (prev.te - prev.ts) < min_err as u32 {
            ck2 -= 1;
        }
    }
    if ck2 == 0 {
        rr.truncate(rr_n0);
        return 0;
    }
    let rn = ck2;
    let mut cur = SvCursor {
        ck: 0,
        xk: window.x_start as i64,
        yk: window.y_start as i64,
    };
    let mut rn1 = 0usize;
    for k in 0..rn {
        let mut cp = rr[rr_n0 + k];
        let (err, mut ey0, mut ey1) = iter_sub_cigar_sv(
            zid as u32, &mut cp, cells, &mut cur, None, None, -1, -1, -1, -1, -1, -1,
        );
        cp.qn = err as u32;
        if cp.qn < min_err as u32 {
            continue;
        }
        assert!(ey0 >= 0 && ey1 >= 0 && ey1 >= ey0);
        ey0 -= hpc_len;
        if ey0 < 0 {
            ey0 = 0;
        }
        ey1 += hpc_len;
        if ey1 > tl {
            ey1 = tl;
        }
        let tstr = reads.read_sub_region(z.y_identifier as usize, ey0, ey1 - ey0, z.y_strand as u8);
        iter_sub_cigar_sv(
            zid as u32,
            &mut cp,
            cells,
            &mut cur,
            Some(qstr),
            Some(&tstr),
            ey0,
            ey1,
            hpc_len,
            hpc_rr,
            (min_err << 1) as i64,
            ql,
        );
        rr[rr_n0 + rn1] = cp;
        rn1 += 1;
    }
    rr.truncate(rr_n0 + rn1);
    if rn1 > 0 {
        1
    } else {
        0
    }
}

/// `gen_ov_lidel_variant` (`Correct.cpp` 19986): the window cursors of the
/// match overlaps, sorted by query start then end.
fn gen_ov_lidel_variant(cz: &mut Vec<UlOv>, overlaps: &[OverlapRegion], idx: &mut Vec<u64>, bd: i64) {
    let cbn = cz.len();
    idx.clear();
    for (k, z) in overlaps.iter().enumerate() {
        let zwn = z.window_cigar.windows.len();
        if zwn == 0 || z.is_match != 1 {
            continue;
        }
        for (i, window) in z.window_cigar.windows.iter().enumerate() {
            if window.y_end == -1 {
                continue;
            }
            let q0 = window.x_start as i64 + bd;
            let q1 = window.x_end as i64 - bd;
            if q1 >= q0 {
                idx.push(((q0 as u64) << 32) | (cz.len() - cbn) as u64);
                cz.push(UlOv {
                    tn: k as u32,
                    qn: i as u32,
                    qs: window.x_start as u32,
                    ts: window.y_start as u32,
                    te: 0,
                    qe: 0,
                    sec: bd as u32,
                    el: false,
                    rev: false,
                });
            }
        }
    }
    radix_sort_u64(idx);
    let n = idx.len();
    let mut i = 0usize;
    let mut k = 1usize;
    while k <= n {
        if k == n || (idx[k] >> 32) != (idx[i] >> 32) {
            if k - i > 1 {
                for t in i..k {
                    let cursor = &cz[cbn + idx[t] as u32 as usize];
                    let window = &overlaps[cursor.tn as usize].window_cigar.windows[cursor.qn as usize];
                    let m = ((window.x_end as i64 + 1 - bd) as u64) << 32 | (idx[t] as u32 as u64);
                    idx[t] = m;
                }
                radix_sort_u64(&mut idx[i..k]);
            }
            i = k;
        }
        k += 1;
    }
}

/// `push_idel_info` (`Correct.cpp` 19929).
#[allow(clippy::too_many_arguments)]
fn push_idel_info(
    hap: &mut HaplotypeEvidenceAlloc,
    s: u64,
    e: u64,
    err: u64,
    ov: &mut [UlOv],
    sv: &[UlOv],
    overlaps: &[OverlapRegion],
    idx: &[u64],
    het_cov: i64,
    hom_cov: i64,
    n_hap: i64,
    cut_rate: f64,
    cut_bd: i64,
) {
    let bl0 = hap.length;
    let mut occ0 = 0i64;
    let mut occ2 = 0i64;
    for &key in idx {
        let p = &mut ov[key as usize];
        let z = &overlaps[p.tn as usize];
        if z.is_match != 1 {
            continue;
        }
        let window = &z.window_cigar.windows[p.qn as usize];
        let q0 = window.x_start as i64 + p.sec as i64;
        let q1 = window.x_end as i64 + 1 - p.sec as i64;
        let os = q0.max(s as i64);
        let oe = q1.min(e as i64);
        if oe > os && (oe - os) >= ((e as i64 - s as i64) - (oe - os)) {
            let ec = extract_sub_err(z, s as i64, e as i64, os, oe, err as i64, 0.2, p);
            if ec >= 0 {
                occ0 += 1;
                hap.add(HaplotypeEvidence {
                    site: s as u32,
                    overlap_identifier: p.tn,
                    overlap_site: hap.snp_stat.len() as u32,
                    evidence_type: 0,
                    coverage: ec as u32,
                    mismatch_base: 0,
                });
            } else {
                occ2 += 1;
            }
        }
    }
    let occ1 = sv.len() as i64;
    for entry in sv {
        hap.add(HaplotypeEvidence {
            site: s as u32,
            overlap_identifier: entry.tn,
            overlap_site: hap.snp_stat.len() as u32,
            evidence_type: 1,
            coverage: entry.qn,
            mismatch_base: 0,
        });
    }
    occ0 += 1;
    let mut cc = if het_cov > 0 { het_cov } else { hom_cov / n_hap } as f64;
    cc *= cut_rate;
    let cc = (cc as i64).max(cut_bd);
    if occ0 < cc || occ1 < cc {
        hap.length = bl0;
        return;
    }
    hap.snp_stat.push(SnpStats {
        identifier: hap.snp_stat.len() as u32,
        overlap_num: 0,
        occ_0: occ0 as u32,
        occ_1: occ1 as u32,
        occ_2: occ2 as u32,
        homopolymer_num: 0,
        non_homopolymer_num: 0,
        score: -1,
        site: s as u32,
        is_homopolymer: 0,
    });
}

/// `rcall_lidel_variant` (`Correct.cpp` 20028).
#[allow(clippy::too_many_arguments)]
fn rcall_lidel_variant(
    cz: &mut Vec<UlOv>,
    overlaps: &mut [OverlapRegion],
    hap: &mut HaplotypeEvidenceAlloc,
    idx: &mut Vec<u64>,
    idz: &mut Vec<u64>,
    het_cov: i64,
    hom_cov: i64,
    n_hap: i64,
) -> u64 {
    hap.length = 0;
    hap.snp_stat.clear();
    idx.clear();
    idz.clear();
    let sv_n = cz.len();

    // Consensus interval/error per cluster.
    let mut i = 0usize;
    let mut k = 1usize;
    while k <= sv_n {
        if k == sv_n || cz[k].ts != cz[i].ts {
            let mut nec = [0u64; 2];
            // `bn0` stays zero in the reference: `idx->n = bn0` resets the
            // scratch at the start of every sub-loop.
            let bn0 = 0usize;
            idx.clear();
            let mut pmm = u64::MAX;
            let mut pmn = 0u64;
            for zi in i..k {
                nec[cz[zi].el as usize] += 1;
                let p = ((cz[zi].qs as u64) << 32) | cz[zi].qe as u64;
                idx.push(p);
                if pmm == u64::MAX {
                    pmm = p;
                }
                if pmm == p {
                    pmn += 1;
                }
            }
            if nec[0] < nec[1] {
                if pmn == 0 || pmn < ((idx.len() - bn0) as u64 - pmn) {
                    let ia = &mut idx[bn0..];
                    if k - i > 1 {
                        radix_sort_u64(ia);
                    }
                    let in_len = ia.len();
                    let mut zi = 0usize;
                    let mut zk = 1usize;
                    pmm = u64::MAX;
                    pmn = 0;
                    while zk <= in_len {
                        if zk == in_len || ia[zi] != ia[zk] {
                            if (zk - zi) as u64 > pmn {
                                pmn = (zk - zi) as u64;
                                pmm = ia[zi];
                            }
                            zi = zk;
                        }
                        zk += 1;
                    }
                }
                let mut ez = UlOv {
                    qs: (pmm >> 32) as u32,
                    qe: pmm as u32,
                    ts: i as u32,
                    te: k as u32,
                    qn: 0,
                    tn: cz[i].ts,
                    sec: 0,
                    el: false,
                    rev: false,
                };
                ez.qn = 0;
                idx.truncate(bn0);
                let mut pmm = u64::MAX;
                let mut pmn = 0u64;
                for zi in i..k {
                    if cz[zi].qs != ez.qs || cz[zi].qe != ez.qe {
                        continue;
                    }
                    idx.push(cz[zi].qn as u64);
                    if pmm == u64::MAX {
                        pmm = cz[zi].qn as u64;
                    }
                    if pmm == cz[zi].qn as u64 {
                        pmn += 1;
                    }
                }
                if pmn == 0 || pmn < ((idx.len() - bn0) as u64 - pmn) {
                    let ia = &mut idx[bn0..];
                    if k - i > 1 {
                        radix_sort_u64(ia);
                    }
                    let in_len = ia.len();
                    let mut zi = 0usize;
                    let mut zk = 1usize;
                    pmm = u64::MAX;
                    pmn = 0;
                    while zk <= in_len {
                        if zk == in_len || ia[zi] != ia[zk] {
                            if (zk - zi) as u64 > pmn {
                                pmn = (zk - zi) as u64;
                                pmm = ia[zi];
                            }
                            zi = zk;
                        }
                        zk += 1;
                    }
                }
                ez.qn = pmm as u32;
                cz.push(ez);
            }
            i = k;
        }
        k += 1;
    }

    let svi_n = cz.len() - sv_n;
    idx.clear();
    gen_ov_lidel_variant(cz, overlaps, idx, 0);
    let ovn = cz.len() - sv_n - svi_n;

    let svi_start = sv_n;
    let svi_end = sv_n + svi_n;
    let ov_start = svi_end;
    radix_sort_by_key_32(&mut cz[svi_start..svi_end], &|v: &UlOv| v.qs);

    let mut cursor = 0usize;
    for k in 0..svi_n {
        let s = cz[svi_start + k].qs as i64;
        let e = cz[svi_start + k].qe as i64;
        let cluster_begin = cz[svi_start + k].ts as usize;
        let cluster_end = cz[svi_start + k].te as usize;
        for p in cluster_begin..cluster_end {
            assert!(cz[p].ts == cz[svi_start + k].tn);
            let oid = cz[p].tn as usize;
            assert!(overlaps[oid].is_match == 1);
            overlaps[oid].is_match = 2;
        }

        // Drop the cursors that end before the variant.
        let mut rm_n = ovn;
        for m in ovn..idx.len() {
            let cp = &cz[ov_start + idx[m] as usize];
            let z = &overlaps[cp.tn as usize];
            let window = &z.window_cigar.windows[cp.qn as usize];
            let q1 = window.x_end as i64 + 1 - cp.sec as i64;
            if q1 <= s {
                continue;
            }
            idx[rm_n] = idx[m];
            rm_n += 1;
        }
        idx.truncate(rm_n);

        while cursor < ovn {
            let cp = &cz[ov_start + idx[cursor] as u32 as usize];
            let z = &overlaps[cp.tn as usize];
            let window = &z.window_cigar.windows[cp.qn as usize];
            let q0 = window.x_start as i64 + cp.sec as i64;
            if q0 >= e {
                break;
            }
            let q1 = window.x_end as i64 + 1 - cp.sec as i64;
            let os = q0.max(s);
            let oe = q1.min(e);
            if oe > os {
                idx.push(idx[cursor] as u32 as u64);
            }
            cursor += 1;
        }

        let sv_copy = cz[cluster_begin..cluster_end].to_vec();
        let variant_err = cz[svi_start + k].qn as u64;
        let ov_slice = &mut cz[ov_start..ov_start + ovn];
        push_idel_info(
            hap,
            s as u64,
            e as u64,
            variant_err,
            ov_slice,
            &sv_copy,
            overlaps,
            &idx[ovn..],
            het_cov,
            hom_cov,
            n_hap,
            0.333333,
            5,
        );

        for p in cluster_begin..cluster_end {
            let oid = cz[p].tn as usize;
            assert!(overlaps[oid].is_match == 2);
            overlaps[oid].is_match = 1;
        }
    }
    hap.length as u64
}

/// `rphase_lidel` (`Correct.cpp` 20155): the ONT indel pass. Returns the number
/// of evidence records, zero when no large indel was found.
#[allow(clippy::too_many_arguments)]
pub fn rphase_lidel(
    overlaps: &mut [OverlapRegion],
    reads: &ReadStore,
    hap: &mut HaplotypeEvidenceAlloc,
    query: &[u8],
    c_idx: &mut Vec<UlOv>,
    idx: &mut Vec<u64>,
    buf: &mut Vec<u64>,
    bd: i64,
    hpc_len: u64,
    het_cov: i64,
    hom_cov: i64,
    n_hap: i64,
) -> u64 {
    hap.length = 0;
    hap.snp_stat.clear();
    c_idx.clear();
    idx.clear();
    for (k, z) in overlaps.iter().enumerate() {
        let zwn = z.window_cigar.windows.len();
        if zwn == 0 || z.is_match != 1 {
            continue;
        }
        for i in 0..zwn {
            let window = &z.window_cigar.windows[i];
            if window.y_end == -1 {
                continue;
            }
            let q0 = window.x_start as i64;
            let q1 = window.x_end as i64 + 1;
            let t0 = window.y_start as i64;
            let t1 = window.y_end as i64 + 1;
            if q1 > q0 && t1 > t0 {
                extract_sub_cigar_sv(z, k, i, reads, query, c_idx, 16, hpc_len as i64, 2);
            }
        }
    }
    if c_idx.is_empty() {
        return 0;
    }
    radix_sort_by_key_32(c_idx.as_mut_slice(), &|v: &UlOv| v.qs);
    // The reference's secondary `qe` sort is guarded by `k > i + k`, which is
    // never true, so the stable `qs` order is what reaches the cluster pass.
    let on = rphase_lidel_cc(c_idx, overlaps.len(), 0.500001, 3, 0.25, 3, buf, idx);
    c_idx.truncate(on as usize);
    rcall_lidel_variant(c_idx, overlaps, hap, idx, buf, het_cov, hom_cov, n_hap)
}

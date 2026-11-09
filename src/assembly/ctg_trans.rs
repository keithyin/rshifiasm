//! The UL contig-transition chaining used by `--dual-scaf`
//! (`gen_contig_self` in `inter.cpp` 21941 and its callees).
//!
//! A contig (scaffold) is mapped onto the merged haplotype unitig graph by
//! matching its blocks against the reference scaffold results: first bubble
//! matches, then exact unitig matches, then the read-to-unitig transitions
//! collected by the shallow/deep passes. The matches are chained by
//! `linear_ctg_trans_chain_dp` and filtered by `filter_ctg_trans`.

use crate::assembly::asg::MaUtg;
use crate::assembly::bubble_chain::{get_bubbles, BubbleType};
use crate::assembly::dual_scaf::TransTable;
use crate::assembly::hap_cov::{UTrans, CHAIN_MATCH, RC_0, RC_1};
use crate::assembly::trans_clean::kt_u_trans_t_idx;
use crate::assembly::unitig::UnitigGraph;
use crate::assembly::ul::{
    get_r_ug_region, UlOv, UlVec, UtgRid, UtgRidDt, MAX_SEC_ERR,
};
use crate::assembly::ul::ScafRes;

/// `G_CHAIN_BW` / `UG_SKIP_N` / `UG_ITER_N` / `UG_DIS_N` (`inter.cpp`).
pub const G_CHAIN_BW: i64 = 16;
pub const UG_SKIP_N: i64 = 100;
pub const UG_ITER_N: i64 = 5000;
pub const UG_DIS_N: i64 = 50000;
/// `asm_opt.ul_error_rate` default (`CommandLines.cpp` 335).
pub const UL_ERROR_RATE: f64 = 0.2;

const INT32_MIN_I64: i64 = i32::MIN as i64;

/// `Chain_Data` (`Hash_Table.h` 122): the DP scratch buffers.
#[derive(Default)]
pub struct ChainData {
    pub score: Vec<i32>,
    pub pre: Vec<i64>,
    pub indels: Vec<i32>,
    pub tmp: Vec<i64>,
}

impl ChainData {
    fn resize(&mut self, n: usize) {
        self.score.resize(n, 0);
        self.pre.resize(n, 0);
        self.indels.resize(n, 0);
        self.tmp.resize(n, 0);
    }
}

/// `gen_r_ug_sc_idx` (`inter.cpp` 21352): the scaffold index
/// unitig -> (scaffold, occurrence).
pub fn gen_r_ug_sc_idx(scaffolds: &ScafRes, unitig_count: usize) -> UtgRid {
    let mut idx = vec![0u32; unitig_count + 1];
    for idx_ul in &scaffolds.a {
        for z in 0..idx_ul.n {
            let header = idx_ul.bb[z];
            for k in header.ts as usize..header.te as usize {
                idx[idx_ul.bb[k].hid as usize] += 1;
            }
        }
    }
    let mut l = 0u32;
    for k in 0..=unitig_count {
        let count = idx[k];
        idx[k] = l;
        l += count;
    }
    let mut p = vec![UtgRidDt::default(); l as usize];
    for (i, idx_ul) in scaffolds.a.iter().enumerate() {
        for z in 0..idx_ul.n {
            let header = idx_ul.bb[z];
            for k in header.ts as usize..header.te as usize {
                let rid = idx_ul.bb[k].hid as usize;
                let base = idx[rid] as usize;
                let a_n = (idx[rid + 1] - idx[rid]) as usize;
                if a_n == 0 {
                    continue;
                }
                let last = base + a_n - 1;
                let slot = if p[last].off as usize == a_n - 1 {
                    last
                } else {
                    let cursor = base + p[last].off as usize;
                    p[last].off += 1;
                    cursor
                };
                p[slot].u = ((i as u32) << 1) | u32::from(idx_ul.bb[k].rev);
                p[slot].pos = idx_ul.bb[k].qs;
                p[slot].off = k as u32;
            }
        }
    }
    UtgRid { idx, p }
}

/// `get_seq_bub_id` (`Overlaps.cpp` 16620): the bubble whose end a vertex is.
fn get_seq_bub_id(v: u32, bubble: &BubbleType) -> Option<(u64, u32)> {
    let (id0, id1) = bubble.get_bub_id(v >> 1, 1);
    let f_bub = bubble.full_bubble_count;
    for id in [id0, id1] {
        if id == u64::MAX || id >= f_bub {
            continue;
        }
        let (beg, sink, _, _) = get_bubbles(bubble, id);
        if v == beg {
            return Some((id, sink ^ 1));
        }
        if v == sink {
            return Some((id, beg ^ 1));
        }
    }
    None
}

/// `pick_bubble` (`Overlaps.cpp` 16613): how many blocks of a path form one
/// small bubble (0 when the head is not a bubble beg).
fn pick_bubble(blocks: &[crate::assembly::ul::UcBlock], bubble: &BubbleType) -> (usize, u64) {
    if blocks.len() <= 1 {
        return (0, u64::MAX);
    }
    let v = (blocks[0].hid << 1) | u32::from(blocks[0].rev);
    let Some((bid, w)) = get_seq_bub_id(v, bubble) else {
        return (0, u64::MAX);
    };
    for k in 1..blocks.len() {
        let v = (blocks[k].hid << 1) | u32::from(blocks[k].rev);
        if v == w {
            return (k, bid);
        }
        if bubble.index[(v >> 1) as usize] as u64 != bid {
            return (0, u64::MAX);
        }
    }
    (0, bid)
}

/// `qry_bub_sc` (`Overlaps.cpp` 16629): the two reference occurrences that
/// bracket the query bubble.
#[allow(clippy::too_many_arguments)]
fn qry_bub_sc(
    index: &UtgRid,
    reference: &UnitigGraph,
    reference_scaf: &ScafRes,
    bubble: &BubbleType,
    bid: u64,
    qsidx: usize,
    qeidx: usize,
    qstr: &UlVec,
    avoid_id: u32,
    res: &mut Vec<UlOv>,
) -> bool {
    let (beg, sink, _, _) = get_bubbles(bubble, bid);
    if beg == u32::MAX || sink == u32::MAX {
        return false;
    }
    let qv = (qstr.bb[qsidx].hid << 1) | u32::from(qstr.bb[qsidx].rev);
    let qw = (qstr.bb[qeidx].hid << 1) | u32::from(qstr.bb[qeidx].rev);
    let qrev = if qv == beg && qw == sink ^ 1 {
        0
    } else if qv == sink && qw == beg ^ 1 {
        1
    } else {
        panic!("qry_bub_sc: query does not match the bubble");
    };
    let mut count = 0usize;
    for (v, w, trev) in [(beg, sink ^ 1, 0u32), (sink, beg ^ 1, 1u32)] {
        let region = get_r_ug_region(index, (v >> 1) as usize);
        for entry in region {
            if entry.u & 1 != v & 1 {
                continue;
            }
            if entry.u >> 1 == avoid_id {
                continue;
            }
            let q = &reference_scaf.a[(entry.u >> 1) as usize];
            let si = entry.off as usize;
            let mut ei = si;
            debug_assert_eq!(q.bb[si].hid, v >> 1);
            debug_assert_eq!(u32::from(q.bb[si].rev), v & 1);
            let mut tn = 0usize;
            for k in 0..q.n {
                tn += (q.bb[k].te - q.bb[k].ts) as usize;
            }
            let mut k = si + 1;
            while k < tn && q.bb[k].aidx == q.bb[si].aidx {
                let z = (q.bb[k].hid << 1) | u32::from(q.bb[k].rev);
                if z == w {
                    ei = k;
                    break;
                }
                if bubble.index[(z >> 1) as usize] as u64 != bid {
                    break;
                }
                k += 1;
            }
            if ei > si {
                res.push(UlOv {
                    qn: q.bb[si].aidx,
                    qs: qsidx as u32,
                    qe: qeidx as u32,
                    tn: entry.u >> 1,
                    ts: si as u32,
                    te: ei as u32,
                    rev: (qrev == trev),
                    ..UlOv::default()
                });
                let _ = reference;
                count += 1;
            }
        }
    }
    count > 0
}

/// `cc_bub_match` (`inter.cpp` 16718): merges consecutive bubble matches on
/// the same target into one alignment.
fn cc_bub_match(ch: &mut [UlOv]) -> usize {
    if ch.is_empty() {
        return 0;
    }
    ch.sort_by_key(|item| item.qs);
    let n = ch.len();
    let mut j = 0usize;
    let mut i = 1usize;
    while i <= n {
        if i == n || ch[i].qs != ch[j].qs {
            if i - j > 1 {
                ch[j..i].sort_by_key(|item| item.qe);
            }
            j = i;
        }
        i += 1;
    }

    let mut k = 0usize;
    for i in 0..n {
        if ch[i].tn == u32::MAX {
            continue;
        }
        let mut j = i + 1;
        while j < n && ch[i].qe >= ch[j].qs {
            if ch[j].tn != u32::MAX
                && ch[i].qe == ch[j].qs
                && ch[i].rev == ch[j].rev
                && ch[i].qn == ch[j].qn
            {
                if !ch[i].rev {
                    if ch[i].te == ch[j].ts {
                        ch[i].qe = ch[j].qe;
                        ch[i].te = ch[j].te;
                        ch[j].tn = u32::MAX;
                    }
                } else if ch[j].te == ch[i].ts {
                    ch[i].qe = ch[j].qe;
                    ch[i].ts = ch[j].ts;
                    ch[j].tn = u32::MAX;
                }
            }
            j += 1;
        }
        ch[k] = ch[i];
        k += 1;
    }
    k
}

/// `extract_scaf_res_t` (`inter.cpp` 16760): exact unitig matches.
fn extract_scaf_res_t(
    block: &crate::assembly::ul::UcBlock,
    index: &UtgRid,
    reference_scaf: &ScafRes,
    avid: u32,
    res: &mut Vec<UlOv>,
) -> usize {
    let region = get_r_ug_region(index, block.hid as usize);
    let mut count = 0usize;
    for entry in region {
        if entry.u >> 1 == avid {
            continue;
        }
        let ou = &reference_scaf.a[(entry.u >> 1) as usize].bb[entry.off as usize];
        let nw = f64::from((ou.qe - ou.qs).max(block.qe - block.qs)) * CHAIN_MATCH;
        res.push(UlOv {
            qn: u32::MAX,
            qs: block.qs,
            qe: block.qe,
            tn: entry.u >> 1,
            ts: ou.qs,
            te: ou.qe,
            sec: if nw >= f64::from(MAX_SEC_ERR) {
                MAX_SEC_ERR
            } else {
                nw as u32
            },
            el: true,
            rev: block.rev != ou.rev,
        });
        count += 1;
    }
    count
}

/// `gen_trans_ovlp_scaf` (`inter.cpp` 16778): maps one transition onto the
/// query/reference blocks.
fn gen_trans_ovlp_scaf(
    gfa: &UnitigGraph,
    map: &UTrans,
    qin: &crate::assembly::ul::UcBlock,
    rin: &crate::assembly::ul::UcBlock,
    res: &mut UlOv,
) -> bool {
    let qlen = gfa.unitigs[qin.hid as usize].length;
    let (s0, e0) = (map.query_start, map.query_end);
    if s0 > qlen || e0 > qlen {
        return false;
    }
    let (s, e) = if !qin.rev {
        (s0, e0)
    } else {
        (qlen - e0, qlen - s0)
    };
    let (s, e) = (s + qin.qs, e + qin.qs);
    if s >= e {
        return false;
    }
    res.qs = s;
    res.qe = e;

    let tlen = gfa.unitigs[rin.hid as usize].length;
    let (s0, e0) = (map.target_start, map.target_end);
    if s0 > tlen || e0 > tlen {
        return false;
    }
    let (s, e) = if !rin.rev {
        (s0, e0)
    } else {
        (tlen - e0, tlen - s0)
    };
    let (s, e) = (s + rin.qs, e + rin.qs);
    if s >= e {
        return false;
    }
    res.ts = s;
    res.te = e;
    res.rev = map.reverse ^ (qin.rev ^ rin.rev);
    true
}

/// `extract_trans_scaf_res_t` (`inter.cpp` 16800): transition matches.
#[allow(clippy::too_many_arguments)]
fn extract_trans_scaf_res_t(
    qin: &crate::assembly::ul::UcBlock,
    index: &UtgRid,
    reference_scaf: &ScafRes,
    gfa: &UnitigGraph,
    ta: &TransTable,
    avid: u32,
    is_exact: bool,
    res: &mut Vec<UlOv>,
) {
    for u in ta.span(qin.hid) {
        if is_exact && u.flag != RC_0 && u.flag != RC_1 {
            continue;
        }
        let region = get_r_ug_region(index, u.target_identifier as usize);
        for entry in region {
            if entry.u >> 1 == avid {
                continue;
            }
            let rin = &reference_scaf.a[(entry.u >> 1) as usize].bb[entry.off as usize];
            let mut p = UlOv::default();
            if gen_trans_ovlp_scaf(gfa, u, qin, rin, &mut p) {
                p.tn = entry.u >> 1;
                p.qn = u32::from(u.flag);
                p.el = false;
                p.sec = u.network_weight as u32;
                res.push(p);
            }
        }
    }
}

/// `cl_trans_gen` (`inter.cpp` 16836): the initial query alignments.
#[allow(clippy::too_many_arguments)]
pub fn cl_trans_gen(
    qid: u32,
    qlen: u32,
    qstr: &UlVec,
    index: &UtgRid,
    reference: &UnitigGraph,
    reference_scaf: &ScafRes,
    gfa: &UnitigGraph,
    bubble: &BubbleType,
    ta: &TransTable,
    is_self: bool,
    is_exact: bool,
) -> Vec<UlOv> {
    let _ = qlen;
    let avid = if is_self { qid } else { u32::MAX };
    let mut res: Vec<UlOv> = Vec::new();

    for k in 0..qstr.n {
        let header = qstr.bb[k];
        let mut z = header.ts as usize;
        let a_n = header.te as usize;
        while z < a_n {
            let (zn, bid) = pick_bubble(&qstr.bb[z..a_n], bubble);
            if zn == 0
                || !qry_bub_sc(
                    index,
                    reference,
                    reference_scaf,
                    bubble,
                    bid,
                    z,
                    z + zn,
                    qstr,
                    avid,
                    &mut res,
                )
            {
                z += 1;
            } else {
                z += zn;
            }
        }
    }

    let mut m = 0usize;
    if !res.is_empty() {
        for item in &mut res {
            item.tn = (item.tn << 1) | u32::from(item.rev);
        }
        res.sort_by_key(|item| item.tn);

        let mut k = 1usize;
        let mut l = 0usize;
        while k <= res.len() {
            if k == res.len() || (res[k].tn >> 1) != (res[l].tn >> 1) {
                for item in &mut res[l..k] {
                    item.tn >>= 1;
                }
                let an = if k - l > 1 {
                    l + cc_bub_match(&mut res[l..k])
                } else {
                    k
                };
                for z in l..an {
                    res[m] = res[z];
                    m += 1;
                }
                l = k;
            }
            k += 1;
        }
        res.truncate(m);
        res.sort_by_key(|item| item.qe);
    }

    for k in 0..qstr.n {
        let header = qstr.bb[k];
        for z in header.ts as usize..header.te as usize {
            let mut mi = 0usize;
            while mi < m && z as u32 <= res[mi].qe {
                if z as u32 >= res[mi].qs && z as u32 <= res[mi].qe {
                    break;
                }
                mi += 1;
            }
            if mi < m && z as u32 >= res[mi].qs && z as u32 <= res[mi].qe {
                continue;
            }
            if extract_scaf_res_t(&qstr.bb[z], index, reference_scaf, avid, &mut res) > 0 {
                continue;
            }
            extract_trans_scaf_res_t(
                &qstr.bb[z],
                index,
                reference_scaf,
                gfa,
                ta,
                avid,
                is_exact,
                &mut res,
            );
        }
    }

    for k in 0..m {
        res[k].el = true;
        res[k].qn = res[k].qs;
        res[k].qs = qstr.bb[res[k].qs as usize].qs;
        res[k].qe = qstr.bb[res[k].qe as usize].qe;
        let target = &reference_scaf.a[res[k].tn as usize];
        res[k].ts = target.bb[res[k].ts as usize].qs;
        res[k].te = target.bb[res[k].te as usize].qe;
        let nw = f64::from((res[k].te - res[k].ts).max(res[k].qe - res[k].qs)) * CHAIN_MATCH;
        res[k].sec = if nw >= f64::from(MAX_SEC_ERR) {
            MAX_SEC_ERR
        } else {
            nw as u32
        };
    }
    res
}

/// `comput_linear_ctg_trans_sc` (`inter.cpp` 16883): the DP step between two
/// alignments, where `li` is the suffix of `lj`.
fn comput_linear_ctg_trans_sc(
    li: &UlOv,
    lj: &UlOv,
    diff_ec_ul: f64,
    ovlp_max: f64,
    bw: i64,
) -> i64 {
    if lj.rev != li.rev {
        return INT32_MIN_I64;
    }
    if i64::from(lj.te) > i64::from(li.te) {
        return INT32_MIN_I64;
    }

    let dq = i64::from(li.qs) - i64::from(lj.qe);
    let dt = i64::from(li.ts) - i64::from(lj.te);
    let dd = (dq - dt).abs();

    let dq = i64::from(li.qe) - i64::from(lj.qs);
    let dt = i64::from(li.te) - i64::from(lj.ts);
    let mut mm = (dq.max(dt) as f64 * diff_ec_ul) as i64;
    if mm < bw {
        mm = bw;
    }
    if dd > mm {
        return INT32_MIN_I64;
    }

    let mut dd = i64::from(li.qe - li.qs).min(i64::from(lj.qe - lj.qs));
    dd = (dd as f64 * ovlp_max) as i64;
    if dd < bw {
        dd = bw;
    }
    let os = i64::from(li.qs.max(lj.qs));
    let oe = i64::from(li.qe.min(lj.qe));
    let mm = if oe > os { oe - os } else { 0 };
    if mm > 0 && mm > dd {
        return INT32_MIN_I64;
    }

    let mut dd = i64::from(li.te - li.ts).min(i64::from(lj.te - lj.ts));
    dd = (dd as f64 * ovlp_max) as i64;
    if dd < bw {
        dd = bw;
    }
    let os = i64::from(li.ts.max(lj.ts));
    let oe = i64::from(li.te.min(lj.te));
    let mm = if oe > os { oe - os } else { 0 };
    if mm > dd {
        return INT32_MIN_I64;
    }

    i64::from(li.sec)
}

/// `adjust_rev_tse` (`inter.cpp` 16364) on a [`UlOv`].
fn adjust_rev_tse(x: &UlOv, tlen: u32) -> (u32, u32) {
    if x.rev {
        (tlen.wrapping_sub(x.te), tlen.wrapping_sub(x.ts))
    } else {
        (x.ts, x.te)
    }
}

/// `extend_end_coord` (`inter.cpp` 13793) on a [`UlOv`].
fn extend_end_coord_ul(li: &UlOv, qlen: i64, rlen: i64) -> (i64, i64, i64, i64) {
    let mut qs = i64::from(li.qs);
    let qe = i64::from(li.qe);
    let mut rs = i64::from(li.ts);
    let mut re = i64::from(li.te);
    let rev = li.rev;
    if rev {
        rs = rlen - i64::from(li.te);
        re = rlen - i64::from(li.ts);
    }
    if qs <= rs {
        rs -= qs;
        qs = 0;
    } else {
        qs -= rs;
        rs = 0;
    }
    let qtail = qlen - qe;
    let rtail = rlen - re;
    let mut qe = qe;
    if qtail <= rtail {
        qe = qlen;
        re += qtail;
    } else {
        re = rlen;
        qe += rtail;
    }
    if rev {
        (qs, qe, rlen - re, rlen - rs)
    } else {
        (qs, qe, rs, re)
    }
}

/// `linear_ctg_trans_chain_dp` (`inter.cpp` 16910): the trans DP. Returns
/// `false` when no chain survives and leaves `ea` untouched on failure.
#[allow(clippy::too_many_arguments)]
fn linear_ctg_trans_chain_dp(
    ch: &mut [UlOv],
    reference: &UnitigGraph,
    diff_ec_ul: f64,
    ovlp_max: f64,
    unmatch_max: f64,
    qlen: i64,
    tlen: i64,
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
    dp: &mut ChainData,
    res: &mut UlOv,
    ea: &mut Vec<u64>,
) -> bool {
    *res = UlOv::default();
    let ch_n = ch.len();
    if ch_n == 0 {
        return false;
    }
    dp.resize(ch_n);
    let ea_n = ea.len();

    for item in ch.iter_mut() {
        item.tn >>= 1;
        if item.rev {
            let (rts, rte) = adjust_rev_tse(item, tlen as u32);
            item.ts = rts;
            item.te = rte;
        }
    }
    let mut j = 0usize;
    let mut i = 1usize;
    while i <= ch_n {
        if i == ch_n || ch[i].rev != ch[j].rev {
            if i - j > 1 {
                ch[j..i].sort_by_key(|item| item.qe);
                let mut k = j + 1;
                let mut l = j;
                while k <= i {
                    if k == i || ch[k].qe != ch[l].qe {
                        if k - l > 1 {
                            ch[l..k].sort_by_key(|item| item.qs);
                        }
                        l = k;
                    }
                    k += 1;
                }
            }
            j = i;
        }
        i += 1;
    }

    let n = ch_n;
    dp.indels[..n].fill(0);
    let mut peak_sc = INT32_MIN_I64;
    let mut peak_e = INT32_MIN_I64;
    let mut peak_i: i64 = -1;
    let mut max_ii: i64 = -1;
    let mut f = vec![0i64; n];
    let mut e = vec![0i64; n];
    let mut p = vec![-1i64; n];
    let mut t = vec![0i32; n];

    for i in 0..n {
        let li = ch[i];
        let csc = i64::from(li.sec);
        let cse = if li.el || li.qn == u32::from(RC_0) || li.qn == u32::from(RC_1) {
            i64::from(li.qe - li.qs)
        } else {
            0
        };
        let mut mm_sc = csc;
        let mut mm_se = cse;
        let mut mm_idx: i64 = -1;
        let mut n_skip = 0i64;
        let mut st = if (i as i64) < max_iter {
            0i64
        } else {
            i as i64 - max_iter
        };
        while ch[i].rev != ch[st as usize].rev {
            st += 1;
        }
        let mut j = i as i64 - 1;
        while j >= st {
            let lj = &ch[j as usize];
            let sc = comput_linear_ctg_trans_sc(&li, lj, diff_ec_ul, ovlp_max, G_CHAIN_BW);
            if sc == INT32_MIN_I64 {
                j -= 1;
                continue;
            }
            let sc = sc + f[j as usize];
            let se = cse + e[j as usize];
            if sc > mm_sc || (sc == mm_sc && se > mm_se) {
                mm_sc = sc;
                mm_se = se;
                mm_idx = j;
                if n_skip > 0 {
                    n_skip -= 1;
                }
            } else if t[j as usize] as i64 == i as i64 {
                n_skip += 1;
                if n_skip > max_skip {
                    break;
                }
            }
            if p[j as usize] >= 0 {
                t[p[j as usize] as usize] = i as i32;
            }
            j -= 1;
        }
        let end_j = j;

        if max_ii < 0 || i64::from(ch[i].qe) > i64::from(ch[max_ii as usize].qe) + max_dis {
            let mut max = INT32_MIN_I64;
            let mut max_e = INT32_MIN_I64;
            max_ii = -1;
            let mut j = i as i64 - 1;
            while j >= st && i64::from(ch[i].qe) <= max_dis + i64::from(ch[j as usize].qe) {
                if max < f[j as usize] || (max == f[j as usize] && max_e < e[j as usize]) {
                    max = f[j as usize];
                    max_e = e[j as usize];
                    max_ii = j;
                }
                j -= 1;
            }
        }

        if max_ii >= 0 && max_ii < end_j {
            let lj = &ch[max_ii as usize];
            let sc = comput_linear_ctg_trans_sc(&li, lj, diff_ec_ul, ovlp_max, G_CHAIN_BW);
            if sc != INT32_MIN_I64 {
                let sc = sc + f[max_ii as usize];
                let se = cse + e[max_ii as usize];
                if sc > mm_sc || (sc == mm_sc && se > mm_se) {
                    mm_sc = sc;
                    mm_se = se;
                    mm_idx = max_ii;
                }
            }
        }

        f[i] = mm_sc;
        e[i] = mm_se;
        p[i] = mm_idx;
        if max_ii < 0
            || (i64::from(ch[i].qe) <= max_dis + i64::from(ch[max_ii as usize].qe)
                && (f[max_ii as usize] < f[i]
                    || (f[max_ii as usize] == f[i] && e[max_ii as usize] < e[i])))
        {
            max_ii = i as i64;
        }

        if f[i] > peak_sc || (f[i] == peak_sc || e[i] > peak_e) {
            peak_sc = f[i];
            peak_e = e[i];
            peak_i = i as i64;
        }
    }

    if peak_i < 0 {
        return false;
    }

    let mut qs = u64::from(ch[peak_i as usize].qs);
    let mut qgap = 0u64;
    let mut ts = u64::from(ch[peak_i as usize].ts);
    let mut tgap = 0u64;
    let mut i = peak_i;
    *res = ch[peak_i as usize];
    res.qn = if peak_sc > i64::from(u32::MAX) {
        u32::MAX
    } else {
        peak_sc as u32
    };
    res.sec = 0;

    while i >= 0 {
        let item = &ch[i as usize];
        if u64::from(item.qe) >= qs {
            if u64::from(item.qs) < qs {
                qs = u64::from(item.qs);
            }
        } else {
            qgap += qs - u64::from(item.qe);
            qs = u64::from(item.qs);
        }
        if u64::from(item.te) >= ts {
            if u64::from(item.ts) < ts {
                ts = u64::from(item.ts);
            }
        } else {
            tgap += ts - u64::from(item.te);
            ts = u64::from(item.ts);
        }
        if item.el || item.qn == u32::from(RC_0) || item.qn == u32::from(RC_1) {
            if ea.len() > ea_n && u64::from(item.qe) >= (ea[ea.len() - 1] >> 32) {
                if u64::from(item.qs) < (ea[ea.len() - 1] >> 32) {
                    let eqs = u64::from(item.qs);
                    let eqe = ea[ea.len() - 1] as u32;
                    let last = ea.len() - 1;
                    ea[last] = eqs << 32 | u64::from(eqe);
                }
            } else {
                ea.push(u64::from(item.qs) << 32 | u64::from(item.qe));
            }
        }
        if item.qs < res.qs {
            res.qs = item.qs;
        }
        if item.qe > res.qe {
            res.qe = item.qe;
        }
        if item.ts < res.ts {
            res.ts = item.ts;
        }
        if item.te > res.te {
            res.te = item.te;
        }
        i = p[i as usize];
    }

    let (rts, rte) = adjust_rev_tse(res, tlen as u32);
    res.ts = rts;
    res.te = rte;
    let (rqs, rqe, rts, rte) = extend_end_coord_ul(res, qlen, tlen);
    qgap = qgap.wrapping_add((rqe - rqs) as u64).wrapping_sub(u64::from(res.qe - res.qs));
    tgap = tgap.wrapping_add((rte - rts) as u64).wrapping_sub(u64::from(res.te - res.ts));

    if (qgap as f64) > ((rqe - rqs) as f64 * unmatch_max)
        || (tgap as f64) > ((rte - rts) as f64 * unmatch_max)
    {
        ea.truncate(ea_n);
        return false;
    }
    let qrate = ((rqe - rqs) as f64 - qgap as f64) / (rqe - rqs) as f64;
    let trate = ((rte - rts) as f64 - tgap as f64) / (rte - rts) as f64;
    let mut mrate = qrate.min(trate);
    if mrate < 0.0 {
        mrate = 0.0;
    }
    let scaled = f64::from(res.qn) * mrate;
    res.qn = scaled as u32;
    if res.qn == 0 {
        ea.truncate(ea_n);
        return false;
    }
    res.qs = rqs as u32;
    res.qe = rqe as u32;
    res.ts = rts as u32;
    res.te = rte as u32;
    res.sec = MAX_SEC_ERR;
    if ea.len() > ea_n {
        let added = ea.len() - ea_n;
        let half = added >> 1;
        for ai in 0..half {
            ea.swap(ea_n + ai, ea_n + added - ai - 1);
        }
        res.sec = ea.len() as u32;
        ea.push(added as u64);
    }
    true
}

/// `ctg_trans_gp_chain` (`inter.cpp` 17168): chains per target and keeps the
/// best path, collecting reliable intervals into `b`.
#[allow(clippy::too_many_arguments)]
fn ctg_trans_gp_chain(
    qid: u32,
    res: &mut Vec<UlOv>,
    reference: &UnitigGraph,
    qlen: i64,
    dp: &mut ChainData,
    is_self: bool,
    b: &mut Vec<u64>,
) {
    let rn = res.len();
    for item in res.iter_mut() {
        item.tn = (item.tn << 1) | u32::from(item.rev);
    }
    res.sort_by_key(|item| item.tn);

    let avid = if is_self { qid } else { u32::MAX };
    let mut m = 0usize;
    b.clear();
    let mut k = 1usize;
    let mut l = 0usize;
    while k <= rn {
        if k == rn || (res[l].tn >> 1) != (res[k].tn >> 1) {
            let tlen = i64::from(reference.unitigs[(res[l].tn >> 1) as usize].length);
            if (res[l].tn >> 1) != avid {
                let mut rr = UlOv::default();
                if linear_ctg_trans_chain_dp(
                    &mut res[l..k],
                    reference,
                    UL_ERROR_RATE,
                    0.4,
                    0.666666,
                    qlen,
                    tlen,
                    UG_SKIP_N,
                    UG_ITER_N,
                    UG_DIS_N,
                    dp,
                    &mut rr,
                    b,
                ) {
                    res[m] = rr;
                    m += 1;
                }
            }
            l = k;
        }
        k += 1;
    }
    res.truncate(m);
}

/// `cal_reliable_ovlp` (`inter.cpp` 17208).
fn cal_reliable_ovlp(z: &UlOv, ea: &[u64], s: u64, e: u64) -> u64 {
    if u64::from(z.sec) >= ea.len() as u64 {
        return 0;
    }
    let mut length = 0u64;
    let start = z.sec as usize - ea[z.sec as usize] as usize;
    for &item in &ea[start..z.sec as usize] {
        let zs = item >> 32;
        let ze = item as u32 as u64;
        if zs >= e {
            break;
        }
        if ze <= s {
            continue;
        }
        let os = s.max(zs);
        let oe = e.min(ze);
        if oe > os {
            length += oe - os;
        }
    }
    length
}

/// `filter_ctg_trans` (`inter.cpp` 17219): reliability and containment
/// filtering. `is_self` is always true on the dual-scaf path.
#[allow(clippy::too_many_arguments)]
fn filter_ctg_trans(
    aln: &mut Vec<UlOv>,
    aln_e: &[u64],
    reliable_rate: f64,
    reliable_len: u64,
    sec_rate: f64,
    sec_bd: u64,
    sec_sc: f64,
    qid: u32,
    soff: u32,
    is_self: bool,
) {
    if aln.is_empty() {
        return;
    }
    let sce = [0u32, 1u32];
    for item in aln.iter_mut() {
        item.el = false;
        let mut clen = 0u64;
        if !is_self
            || sce[usize::from(qid >= soff)] != sce[usize::from(item.tn >= soff)]
        {
            if (item.sec as usize) < aln_e.len() {
                let start = item.sec as usize - aln_e[item.sec as usize] as usize;
                for &z in &aln_e[start..item.sec as usize] {
                    clen = clen.wrapping_add(u64::from(z as u32)).wrapping_sub(z >> 32);
                }
            }
            if clen > reliable_len
                && clen as f64 > f64::from(item.qe - item.qs) * reliable_rate
            {
                item.el = true;
            }
        }
        item.qn = u32::MAX - item.qn;
    }
    if aln.len() > 1 {
        aln.sort_by_key(|item| item.qn);
    }
    for item in aln.iter_mut() {
        item.qn = u32::MAX - item.qn;
    }

    let n = aln.len();
    let mut sorted: Vec<(u32, usize)> = aln.iter().enumerate().map(|(k, z)| (z.qs, k)).collect();
    sorted.sort();

    for k in 0..n {
        if aln[k].el {
            continue;
        }
        let cz = aln[k];
        let clen = f64::from(cz.qe - cz.qs) * sec_rate;
        let mut bu: Vec<(u32, u32)> = Vec::new();
        let mut wu: Vec<(u32, u32)> = Vec::new();
        let mut l = 0u64;
        let mut wl = 0u64;
        let mut le = 0u64;
        for &(_, z) in &sorted {
            if z == k {
                continue;
            }
            // `if(ss[z] > k) continue;` — only earlier alignments count.
            if z > k {
                continue;
            }
            let mz = aln[z];
            if mz.qs >= cz.qe {
                break;
            }
            if mz.tn == u32::MAX {
                continue;
            }
            if !mz.el && f64::from(cz.qn) > f64::from(mz.qn) * 0.9 {
                continue;
            }
            let os = cz.qs.max(mz.qs);
            let oe = cz.qe.min(mz.qe);
            if oe <= os {
                continue;
            }
            if let Some(last) = bu.last_mut() {
                if last.1 >= os {
                    if oe > last.1 {
                        le += cal_reliable_ovlp(&cz, aln_e, u64::from(last.1), u64::from(oe));
                        l += u64::from(oe - last.1);
                        last.1 = oe;
                    }
                } else {
                    le += cal_reliable_ovlp(&cz, aln_e, u64::from(os), u64::from(oe));
                    l += u64::from(oe - os);
                    bu.push((os, oe));
                }
            } else {
                le += cal_reliable_ovlp(&cz, aln_e, u64::from(os), u64::from(oe));
                l += u64::from(oe - os);
                bu.push((os, oe));
            }

            if let Some(last) = wu.last_mut() {
                if last.1 >= mz.qs {
                    if mz.qe > last.1 {
                        wl += u64::from(mz.qe - last.1);
                        last.1 = mz.qe;
                    }
                } else {
                    wl += u64::from(mz.qe - mz.qs);
                    wu.push((mz.qs, mz.qe));
                }
            } else {
                wl += u64::from(mz.qe - mz.qs);
                wu.push((mz.qs, mz.qe));
            }

            if l > 0 && l > sec_bd && (le as f64) <= l as f64 * 0.5 {
                if l as f64 >= clen || (l as f64) >= (wl as f64 * sec_rate) {
                    if f64::from(cz.qn) <= f64::from(mz.qn) * sec_sc {
                        aln[k].tn = u32::MAX;
                        break;
                    }
                }
            }
        }
    }

    let mut m = 0usize;
    for k in 0..n {
        if aln[k].tn == u32::MAX {
            continue;
        }
        aln[m] = aln[k];
        m += 1;
    }
    aln.truncate(m);
}

/// `push_ctg_trans_res` (`inter.cpp` 17343).
fn push_ctg_trans_res(id: u32, input: &mut [UlOv], out: &mut Vec<UlOv>) {
    for item in input {
        item.sec = if item.qn >= MAX_SEC_ERR {
            MAX_SEC_ERR
        } else {
            item.qn
        };
        item.qn = id;
        out.push(*item);
    }
}

/// `direct_ctg_trans_chain` (`inter.cpp` 17353) with `qry == NULL`.
#[allow(clippy::too_many_arguments)]
pub fn direct_ctg_trans_chain(
    id: u32,
    query_scaf: &UlVec,
    index: &UtgRid,
    reference: &UnitigGraph,
    reference_scaf: &ScafRes,
    gfa: &UnitigGraph,
    bubble: &BubbleType,
    ta: &TransTable,
    soff: u32,
    is_exact: bool,
    dp: &mut ChainData,
    out: &mut Vec<UlOv>,
) -> bool {
    let qlen = reference.unitigs[id as usize].length;
    let mut idx = cl_trans_gen(
        id, qlen, query_scaf, index, reference, reference_scaf, gfa, bubble, ta, true, is_exact,
    );
    if idx.is_empty() {
        return false;
    }
    let mut b0: Vec<u64> = Vec::new();
    ctg_trans_gp_chain(id, &mut idx, reference, i64::from(qlen), dp, true, &mut b0);
    if idx.is_empty() {
        return false;
    }
    filter_ctg_trans(
        &mut idx,
        &b0,
        0.5,
        10000,
        0.85,
        10000,
        0.3,
        id,
        soff,
        true,
    );
    if idx.is_empty() {
        return false;
    }
    push_ctg_trans_res(id, &mut idx, out);
    true
}

/// `work_ctg_path_trans_self` (`inter.cpp` 18983): chains every contig of
/// the merged haplotype against the reference.
pub fn work_ctg_path_trans_self(
    reference: &UnitigGraph,
    reference_scaf: &ScafRes,
    gfa: &UnitigGraph,
    ta: &TransTable,
    soff: u32,
    is_exact: bool,
    bubble: &BubbleType,
) -> TransTable {
    let index = gen_r_ug_sc_idx(reference_scaf, gfa.unitigs.len());
    let mut per_thread: Vec<Vec<UlOv>> = Vec::new();
    let mut dp = ChainData::default();
    for i in 0..reference.unitigs.len() {
        let mut out = Vec::new();
        direct_ctg_trans_chain(
            i as u32,
            &reference_scaf.a[i],
            &index,
            reference,
            reference_scaf,
            gfa,
            bubble,
            ta,
            soff,
            is_exact,
            &mut dp,
            &mut out,
        );
        per_thread.push(out);
    }
    fill_ctg_path_trans_idx(&per_thread, reference.unitigs.len() as u32)
}

/// `fill_ctg_path_trans_idx` (`inter.cpp` 18909): flattens the per-thread
/// lists and indexes them by query unitig, emitting both directions.
pub fn fill_ctg_path_trans_idx(lists: &[Vec<UlOv>], qn: u32) -> TransTable {
    let mut index = vec![u64::MAX; qn as usize];
    let mut total = 0usize;
    for list in lists {
        total += list.len();
    }
    let mut transitions: Vec<UTrans> = Vec::with_capacity(total * 2);
    let mut starts = vec![0usize; lists.len() + 1];
    for (i, list) in lists.iter().enumerate() {
        let mut l = 0usize;
        for k in 1..=list.len() {
            if k == list.len() || list[l].qn != list[k].qn {
                index[list[l].qn as usize] = ((i as u64) << 32) | l as u64;
                l = k;
            }
        }
        starts[i + 1] = starts[i] + list.len();
    }
    let _ = starts;

    for q in 0..qn as usize {
        if index[q] == u64::MAX {
            index[q] = 0;
            continue;
        }
        let thread = (index[q] >> 32) as usize;
        let list = &lists[thread];
        let mut k = (index[q] as u32) as usize;
        while k < list.len() && list[k].qn == q as u32 {
            let z = list[k];
            transitions.push(UTrans {
                query_identifier: z.qn,
                query_start: z.qs,
                query_end: z.qe,
                target_identifier: z.tn,
                target_start: z.ts,
                target_end: z.te,
                reverse: z.rev,
                flag: u8::from(z.el),
                network_weight: f64::from(z.sec),
                occurrence: z.sec,
                deleted: false,
            });
            transitions.push(UTrans {
                query_identifier: z.tn,
                query_start: z.ts,
                query_end: z.te,
                target_identifier: z.qn,
                target_start: z.qs,
                target_end: z.qe,
                reverse: z.rev,
                flag: u8::from(z.el),
                network_weight: f64::from(z.sec),
                occurrence: z.sec,
                deleted: false,
            });
            k += 1;
        }
    }
    TransTable {
        transitions,
        index,
    }
}

/// `get_u_trans_spec` (`Overlaps.cpp` 17000): the mirrored span of a
/// `qn -> tn` group.
fn get_u_trans_spec(ta: &TransTable, qn: u32, tn: u32) -> &[UTrans] {
    let span = ta.span(qn);
    let mut st = 0usize;
    for i in 1..=span.len() {
        if i == span.len() || span[i].target_identifier != span[st].target_identifier {
            if span[st].target_identifier == tn {
                return &span[st..i];
            }
            st = i;
        }
    }
    &[]
}

/// `update_ctg_trans_cc` (`Overlaps.cpp` 21810): keeps one copy of every
/// mirrored transition pair, appends the reverse directions and reindexes.
pub fn update_ctg_trans_cc(ta: &mut TransTable) {
    let n = ta.index.len();
    let mut deleted = vec![false; ta.transitions.len()];
    for id in 0..n {
        let (start, count) = {
            let entry = ta.index[id];
            ((entry >> 32) as usize, entry as u32 as usize)
        };
        let mut l = 0usize;
        let mut k = 1usize;
        while k <= count {
            let target = ta.transitions[start + l].target_identifier;
            if k == count || ta.transitions[start + k].target_identifier != target {
                if k > l {
                    if id as u32 > target {
                        l = k;
                        k += 1;
                        continue;
                    }
                    let (r_start, r_count) = {
                        let entry = ta.index[target as usize];
                        ((entry >> 32) as usize, entry as u32 as usize)
                    };
                    let (mirror_start, mirror_count) =
                        find_mirror_span(&ta.transitions[r_start..r_start + r_count], id as u32);
                    // (is_mirror, index) of the current best entry.
                    let mut best: Option<(bool, usize)> = None;
                    for z in l..k {
                        let cz = ta.transitions[start + z];
                        deleted[start + z] = true;
                        let keep = match best {
                            None => true,
                            Some((is_mirror, idx)) => {
                                let mz = if is_mirror {
                                    ta.transitions[r_start + idx]
                                } else {
                                    ta.transitions[start + idx]
                                };
                                if mz.flag == 1 && cz.flag != 1 {
                                    false
                                } else {
                                    (cz.flag == 1 && mz.flag != 1)
                                        || (cz.flag == mz.flag
                                            && (mz.network_weight < cz.network_weight
                                                || (mz.network_weight == cz.network_weight
                                                    && (mz.query_end - mz.query_start)
                                                        < (cz.query_end - cz.query_start))))
                                }
                            }
                        };
                        if keep {
                            best = Some((false, z));
                        }
                    }
                    for z in 0..mirror_count {
                        let cz = ta.transitions[r_start + z];
                        deleted[r_start + z] = true;
                        let keep = match best {
                            None => true,
                            Some((is_mirror, idx)) => {
                                let mz = if is_mirror {
                                    ta.transitions[r_start + idx]
                                } else {
                                    ta.transitions[start + idx]
                                };
                                if mz.flag == 1 && cz.flag != 1 {
                                    false
                                } else {
                                    (cz.flag == 1 && mz.flag != 1)
                                        || (cz.flag == mz.flag
                                            && (mz.network_weight < cz.network_weight
                                                || (mz.network_weight == cz.network_weight
                                                    && (mz.query_end - mz.query_start)
                                                        < (cz.query_end - cz.query_start))))
                                }
                            }
                        };
                        if keep {
                            best = Some((true, z));
                        }
                    }
                    if let Some((is_mirror, idx)) = best {
                        if is_mirror {
                            deleted[r_start + idx] = false;
                        } else {
                            deleted[start + idx] = false;
                        }
                    }
                }
                l = k;
            }
            k += 1;
        }
    }

    let mut kept: Vec<UTrans> = Vec::new();
    for (k, item) in ta.transitions.iter().enumerate() {
        if !deleted[k] {
            kept.push(*item);
        }
    }
    let m = kept.len();
    for k in 0..m {
        let item = kept[k];
        let mut mirror = item;
        mirror.query_identifier = item.target_identifier;
        mirror.query_start = item.target_start;
        mirror.query_end = item.target_end;
        mirror.target_identifier = item.query_identifier;
        mirror.target_start = item.query_start;
        mirror.target_end = item.query_end;
        kept.push(mirror);
    }
    ta.transitions = kept;
    kt_u_trans_t_idx(&mut ta.transitions, &mut ta.index, n as u32);
}

fn find_mirror_span(span: &[UTrans], tn: u32) -> (usize, usize) {
    let mut st = 0usize;
    for i in 1..=span.len() {
        if i == span.len() || span[i].target_identifier != span[st].target_identifier {
            if span[st].target_identifier == tn {
                return (st, i - st);
            }
            st = i;
        }
    }
    (0, 0)
}

/// The public entry: `gen_contig_self` (`inter.cpp` 21941) restricted to the
/// `is_sep == 0` branch used by `--dual-scaf`.
#[allow(clippy::too_many_arguments)]
pub fn gen_contig_self(
    reference: &UnitigGraph,
    reference_scaf: &ScafRes,
    gfa: &UnitigGraph,
    ta: &TransTable,
    soff: u32,
    bubble: &BubbleType,
    query_unitigs: usize,
) -> TransTable {
    let mut result = work_ctg_path_trans_self(
        reference,
        reference_scaf,
        gfa,
        ta,
        soff,
        false,
        bubble,
    );
    crate::assembly::trans_clean::kt_u_trans_t_idx(
        &mut result.transitions,
        &mut result.index,
        query_unitigs as u32,
    );
    update_ctg_trans_cc(&mut result);
    result
}

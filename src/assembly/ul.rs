//! The ultra-long (`inter.cpp`) scaffolding foundation.
//!
//! Only the slice needed by `--dual-scaf` is ported: the read-to-unitig
//! index (`gen_r_ug_idx`), the contig-to-unitig alignment bridge
//! (`ctg_rg2ug_gen`) and the linear chaining DP (`gen_linear_chains_ctg`).
//! The minimizer-based graph chaining and the UL correction pipeline are not
//! ported yet.

use crate::assembly::asg::MaUtg;

/// `MAX_SEC_ERR`: the `sec` field saturation value (`Correct.cpp` 68).
pub const MAX_SEC_ERR: u32 = 0x3fff_ffff;

/// `asm_opt.ul_error_rate` default (`CommandLines.cpp` 335).
const UL_ERROR_RATE: f64 = 0.2;
/// `G_CHAIN_BW` (`inter.cpp`).
const G_CHAIN_BW: i64 = 16;
/// `N_GCHAIN_RATE` (`inter.cpp`).
const N_GCHAIN_RATE: f64 = 0.04;
/// `UG_SKIP_N` / `UG_ITER_N` / `UG_DIS_N` (`inter.cpp` 18830).
const UG_SKIP_N: i64 = 100;
const UG_ITER_N: i64 = 5000;
const UG_DIS_N: i64 = 50000;
/// `UG_SKIP_GRAPH_N` (`inter.cpp`).
const UG_SKIP_GRAPH_N: i64 = 72;

/// `ul_ov_t` (`Overlaps.h` 58): one unitig alignment of a query read/contig.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct UlOv {
    /// Query id: read id or contig position, depending on the caller.
    pub qn: u32,
    pub qs: u32,
    pub qe: u32,
    /// Target id; `(unitig << 1) | rev` until `gen_linear_chains_ctg` demotes it.
    pub tn: u32,
    pub ts: u32,
    pub te: u32,
    /// `sec:30`.
    pub sec: u32,
    /// `el:1`.
    pub el: bool,
    /// `rev:1`.
    pub rev: bool,
}

/// `uc_block_t` (`Process_Read.h` 168): one block of a contig path.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct UcBlock {
    pub hid: u32,
    pub qs: u32,
    pub qe: u32,
    pub ts: u32,
    pub te: u32,
    pub pidx: u32,
    pub pdis: u32,
    pub aidx: u32,
    /// `pchain:5`.
    pub pchain: u8,
    /// `rev:1`, `base:1`, `el:1`.
    pub rev: bool,
    pub base: bool,
    pub el: bool,
}

/// `ul_vec_t` (`Process_Read.h` 202), reduced to the fields the scaffolding
/// path consumes.
#[derive(Clone, Debug, Default)]
pub struct UlVec {
    pub rlen: u32,
    /// `bb.n`: the number of chain headers. The block data of header `k`
    /// lives at `bb[bb[k].ts..bb[k].te]`, past `n` in the reference's
    /// capacity-only region; `bb` keeps everything in one array here.
    pub n: usize,
    pub bb: Vec<UcBlock>,
}

/// `scaf_res_t` (`Process_Read.h` 216): one [`UlVec`] per input contig.
#[derive(Clone, Debug, Default)]
pub struct ScafRes {
    pub a: Vec<UlVec>,
}

impl ScafRes {
    /// `init_scaf_res_t`.
    pub fn new(n: usize) -> Self {
        Self {
            a: vec![UlVec::default(); n],
        }
    }
}

/// `utg_rid_dt` (`Overlaps.h` 290): one occurrence of a read in a unitig.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct UtgRidDt {
    /// `(unitig << 1) | orientation`.
    pub u: u32,
    /// Occurrence cursor / read index inside the unitig.
    pub off: u32,
    /// Read offset along the unitig when `keep_offset` is off.
    pub pos: u32,
}

/// `utg_rid_t` (`Overlaps.h` 293): CSR read-to-unitig index.
#[derive(Clone, Debug, Default)]
pub struct UtgRid {
    pub idx: Vec<u32>,
    pub p: Vec<UtgRidDt>,
}

/// `gen_r_ug_idx` (`inter.cpp` 21278) with `ct == NULL`.
///
/// `keep_offset` mirrors the reference: the scaffolding index keeps the read
/// index inside the unitig (`off = k`), the UL index keeps the base offset
/// (`off = l`).
pub fn gen_r_ug_idx(unitigs: &[MaUtg], read_count: usize, keep_offset: bool) -> UtgRid {
    let mut idx = vec![0u32; read_count + 1];
    for unitig in unitigs {
        for &packed in &unitig.reads {
            idx[(packed >> 33) as usize] += 1;
        }
    }
    let mut l = 0u32;
    for k in 0..=read_count {
        let count = idx[k];
        idx[k] = l;
        l += count;
    }

    let total = idx[read_count] as usize;
    let mut p = vec![UtgRidDt::default(); total];
    for (i, unitig) in unitigs.iter().enumerate() {
        let mut l = 0u32;
        for (k, &packed) in unitig.reads.iter().enumerate() {
            let rid = (packed >> 33) as usize;
            let base = idx[rid] as usize;
            let a_n = (idx[rid + 1] - idx[rid]) as usize;
            if a_n > 0 {
                let last = base + a_n - 1;
                let slot = if p[last].off as usize == a_n - 1 {
                    last
                } else {
                    let cursor = base + p[last].off as usize;
                    p[last].off += 1;
                    cursor
                };
                p[slot].u = ((i as u32) << 1) | (((packed >> 32) & 1) as u32);
                p[slot].pos = l;
                p[slot].off = if keep_offset { k as u32 } else { l };
            }
            l = l.wrapping_add(packed as u32);
        }
    }

    UtgRid { idx, p }
}

/// `get_r_ug_region` (`inter.cpp` 13268).
pub fn get_r_ug_region<'a>(index: &'a UtgRid, rid: usize) -> &'a [UtgRidDt] {
    let start = index.idx[rid] as usize;
    let end = index.idx[rid + 1] as usize;
    &index.p[start..end]
}

/// `rov2uov_ctg` (`inter.cpp` 16223): convert one read-graph block to a
/// unitig coordinate alignment.
#[allow(clippy::too_many_arguments)]
pub fn rov2uov_ctg(
    read_length: u64,
    unitig_length: u32,
    ru_map: &UtgRidDt,
    read_block: &UcBlock,
    qidx: u32,
    adjust_rev: bool,
) -> Option<UlOv> {
    let ori = ru_map.u & 1;
    let (mut ts, mut te) = if ori == 0 {
        (u64::from(read_block.ts), u64::from(read_block.te))
    } else {
        (
            read_length.wrapping_sub(u64::from(read_block.te)),
            read_length.wrapping_sub(u64::from(read_block.ts)),
        )
    };
    ts = ts.wrapping_add(u64::from(ru_map.pos));
    te = te.wrapping_add(u64::from(ru_map.pos));
    if te > u64::from(unitig_length) {
        return None;
    }
    let mut res = UlOv {
        qn: qidx,
        qs: read_block.qs,
        qe: read_block.qe,
        tn: ru_map.u >> 1,
        ts: ts as u32,
        te: te as u32,
        sec: ru_map.off,
        el: read_block.el,
        rev: read_block.rev != (ori != 0),
    };
    if adjust_rev && res.rev {
        res.ts = unitig_length.wrapping_sub(res.te);
        res.te = unitig_length.wrapping_sub(ts as u32);
    }
    Some(res)
}

/// `get_unique_rctg_aln` (`inter.cpp` 16269). Scaffold (mate-separated)
/// reads never appear in the HiFi graph, so the `IS_SCAF_READ` guard is a
/// no-op.
pub fn get_unique_rctg_aln(
    unitigs: &[MaUtg],
    index: &UtgRid,
    read_lengths: &[u64],
    v: u32,
    vi: u32,
    vs: u64,
    ulen: u64,
    ures: Option<&mut UlOv>,
    mut rr: Option<&mut Vec<UlOv>>,
    sec: u32,
) -> usize {
    let rid = (v >> 1) as usize;
    let read_length = read_lengths[rid];
    let region = get_r_ug_region(index, rid);
    let mut matched = 0usize;
    let mut ures_value = UlOv {
        tn: u32::MAX,
        qn: u32::MAX,
        ..UlOv::default()
    };

    let mut z = UcBlock {
        hid: v >> 1,
        rev: v & 1 != 0,
        el: true,
        qs: vs.min(ulen) as u32,
        qe: (vs + read_length).min(ulen) as u32,
        ts: 0,
        te: read_length as u32,
        ..UcBlock::default()
    };
    if vs + read_length > ulen {
        let over_bd = (z.qe - z.qs).min(read_length as u32);
        if !z.rev {
            z.ts = 0;
            z.te = over_bd;
        } else {
            z.ts = read_length as u32 - over_bd;
            z.te = read_length as u32;
        }
    }
    if z.qe <= z.qs || z.te <= z.ts {
        return 0;
    }

    for entry in region {
        let unitig_length = unitigs[(entry.u >> 1) as usize].length;
        let Some(mut p) = rov2uov_ctg(read_length, unitig_length, entry, &z, vi, true) else {
            continue;
        };
        p.el = true;
        p.tn = (p.tn << 1) | u32::from(p.rev);
        if sec != u32::MAX {
            p.sec = sec;
        }
        if let Some(list) = rr.as_deref_mut() {
            list.push(p);
        }
        ures_value = p;
        matched += 1;
    }
    if matched == 1 {
        if let Some(out) = ures {
            *out = ures_value;
        }
    }
    matched
}

/// `cmp_shrink_ul_ov_t` (`inter.cpp` 16315).
fn cmp_shrink_ul_ov_t(a0: &UlOv, a1: &UlOv, off: u32) -> bool {
    if a0.tn == u32::MAX || a1.tn == u32::MAX {
        return true;
    }
    !(a0.tn == a1.tn && a0.qn.wrapping_add(off) == a1.qn && a0.sec.wrapping_add(off) == a1.sec)
}

/// `ctg_rg2ug_gen` (`inter.cpp` 16323): map one contig's reads to unitig
/// alignments, merging runs of consecutive reads that stay on one unitig.
pub fn ctg_rg2ug_gen(
    contig: &MaUtg,
    u_cl: &mut Vec<UlOv>,
    unitigs: &[MaUtg],
    index: &UtgRid,
    read_lengths: &[u64],
) {
    u_cl.clear();
    if contig.reads.is_empty() {
        return;
    }
    let mut lp = UlOv {
        tn: u32::MAX,
        qn: u32::MAX,
        ..UlOv::default()
    };
    let mut kp = lp;
    let mut kp0: UlOv;

    let mut k = 1usize;
    let mut l = 0usize;
    let mut ls = 0u64;
    let mut vs = 0u64;
    get_unique_rctg_aln(
        unitigs,
        index,
        read_lengths,
        (contig.reads[0] >> 32) as u32,
        0,
        ls,
        u64::from(contig.length),
        Some(&mut lp),
        None,
        u32::MAX,
    );
    vs += contig.reads[0] as u32 as u64;

    while k <= contig.reads.len() {
        let mut is_push = false;
        if k < contig.reads.len() {
            get_unique_rctg_aln(
                unitigs,
                index,
                read_lengths,
                (contig.reads[k] >> 32) as u32,
                k as u32,
                vs,
                u64::from(contig.length),
                Some(&mut kp),
                None,
                u32::MAX,
            );
            if cmp_shrink_ul_ov_t(&lp, &kp, (k - l) as u32) {
                is_push = true;
            }
        } else {
            is_push = true;
        }

        if is_push {
            kp0 = kp;
            if k - l > 1 {
                for i in l..k {
                    let mut merged = UlOv::default();
                    let m = get_unique_rctg_aln(
                        unitigs,
                        index,
                        read_lengths,
                        (contig.reads[i] >> 32) as u32,
                        i as u32,
                        ls,
                        u64::from(contig.length),
                        Some(&mut merged),
                        None,
                        u32::MAX,
                    );
                    debug_assert_eq!(m, 1);
                    debug_assert_eq!(merged.tn, lp.tn);
                    debug_assert_eq!(merged.qn, lp.qn + (i - l) as u32);
                    debug_assert_eq!(merged.sec, lp.sec + (i - l) as u32);
                    if merged.qs < lp.qs {
                        lp.qs = merged.qs;
                    }
                    if merged.qe > lp.qe {
                        lp.qe = merged.qe;
                    }
                    if merged.ts < lp.ts {
                        lp.ts = merged.ts;
                    }
                    if merged.te > lp.te {
                        lp.te = merged.te;
                    }
                    ls += contig.reads[i] as u32 as u64;
                }
                lp.sec = (k - l) as u32;
                u_cl.push(lp);
            } else if k - l == 1 {
                get_unique_rctg_aln(
                    unitigs,
                    index,
                    read_lengths,
                    (contig.reads[l] >> 32) as u32,
                    l as u32,
                    ls,
                    u64::from(contig.length),
                    None,
                    Some(u_cl),
                    1,
                );
            }
            l = k;
            lp = kp0;
            ls = vs;
        }

        if k < contig.reads.len() {
            vs += contig.reads[k] as u32 as u64;
        }
        k += 1;
    }
}

/// `comput_linear_ctg_sc` (`inter.cpp` 16369).
fn comput_linear_ctg_sc(li: &UlOv, lj: &UlOv, diff_ec_ul: f64, bw: i64) -> i64 {
    if lj.te > li.te {
        return i64::from(i32::MIN);
    }
    let dq = i64::from(li.qe.wrapping_sub(lj.qs));
    let dt = i64::from(li.te.wrapping_sub(lj.ts));
    let dd = (dq - dt).abs();
    let mut mm = (dq.max(dt) as f64 * diff_ec_ul) as i64;
    if mm < bw {
        mm = bw;
    }
    if dd > mm {
        i64::from(i32::MIN)
    } else {
        i64::from(li.sec)
    }
}

/// `is_available_ctg_aln` (`inter.cpp` 16380).
fn is_available_ctg_aln(unitig: &MaUtg, z: &UlOv) -> bool {
    if z.te <= z.ts {
        return false;
    }
    let span = z.te - z.ts;
    unitig.length >= span
        && unitig.length - span <= 16
        && f64::from(z.sec) >= f64::from(unitig.reads.len() as u32) * 0.333333
}

/// `adjust_rev_tse` (`inter.cpp` 13349).
fn adjust_rev_tse(x: &UlOv, tlen: u32) -> (u32, u32) {
    if x.rev {
        (tlen.wrapping_sub(x.te), tlen.wrapping_sub(x.ts))
    } else {
        (x.ts, x.te)
    }
}

/// `linear_ctg_chain_dp_adv` (`inter.cpp` 16391): the per-unitig linear
/// chaining DP. `ch` is `res[l..k]`, `sv` the matching scratch slice; the
/// chains are compacted into the front of `ch` and their count is returned.
#[allow(clippy::too_many_arguments)]
fn linear_ctg_chain_dp_adv(
    ch: &mut [UlOv],
    sv: &mut [UlOv],
    unitigs: &[MaUtg],
    bw: i64,
    diff_ec_ul: f64,
    qlen: i64,
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
    chain_offset: i64,
) -> usize {
    let ch_n = ch.len();
    if ch_n == 0 {
        return 0;
    }

    ch.sort_by_key(|item| item.qe);
    let mut j = 0usize;
    for i in 1..=ch_n {
        if i == ch_n || ch[i].qe != ch[j].qe {
            if i - j > 1 {
                ch[j..i].sort_by_key(|item| item.qs);
            }
            j = i;
        }
    }

    let mut t = vec![0i64; ch_n];
    let mut f = vec![0i32; ch_n];
    let mut p = vec![-1i64; ch_n];

    let int32_min = i64::from(i32::MIN);
    let mut max_ii: i64 = -1;
    for i in 0..ch_n {
        let li = ch[i];
        let csc = li.sec as i64;
        let mut mm_sc = csc;
        let mut mm_idx: i64 = -1;
        let mut n_skip: i64 = 0;
        let st = if (i as i64) < max_iter { 0 } else { i as i64 - max_iter };
        let mut j = i as i64 - 1;
        while j >= st {
            let lj = ch[j as usize];
            let mut sc = comput_linear_ctg_sc(&li, &lj, diff_ec_ul, bw);
            if sc != int32_min {
                sc += i64::from(f[j as usize]);
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
            }
            j -= 1;
        }
        let end_j = j;

        if max_ii < 0 || i64::from(ch[i].qe) > i64::from(ch[max_ii as usize].qe) + max_dis {
            let mut max = i32::MIN;
            max_ii = -1;
            let mut j = i as i64 - 1;
            while j >= st && i64::from(ch[i].qe) <= max_dis + i64::from(ch[j as usize].qe) {
                if max < f[j as usize] {
                    max = f[j as usize];
                    max_ii = j;
                }
                j -= 1;
            }
        }

        if max_ii >= 0 && max_ii < end_j {
            let lj = ch[max_ii as usize];
            let mut sc = comput_linear_ctg_sc(&li, &lj, diff_ec_ul, bw);
            if sc != int32_min {
                sc += i64::from(f[max_ii as usize]);
                if sc > mm_sc {
                    mm_sc = sc;
                    mm_idx = max_ii;
                }
            }
        }

        f[i] = mm_sc as i32;
        p[i] = mm_idx;
        if max_ii < 0
            || (i64::from(ch[i].qe) <= max_dis + i64::from(ch[max_ii as usize].qe)
                && f[max_ii as usize] < f[i])
        {
            max_ii = i as i64;
        }
        sv[i] = li;
    }

    for value in t.iter_mut() {
        *value = 0;
    }
    let mut n_u = 0usize;
    for k in (0..ch_n).rev() {
        if t[k] != 0 {
            continue;
        }
        let mut i = k as i64;
        let mut chain = sv[k];
        let mut sc = i64::from(f[k]);
        while i >= 0 {
            let item = sv[i as usize];
            if item.qs < chain.qs {
                chain.qs = item.qs;
            }
            if item.ts < chain.ts {
                chain.ts = item.ts;
            }
            if item.qe > chain.qe {
                chain.qe = item.qe;
            }
            if item.te > chain.te {
                chain.te = item.te;
            }
            t[i as usize] = 1;
            i = p[i as usize];
        }
        let tlen = unitigs[chain.tn as usize].length;
        let (its, ite) = adjust_rev_tse(&chain, tlen);
        chain.ts = its;
        chain.te = ite;
        if i >= 0 {
            sc -= i64::from(f[i as usize]);
        }
        chain.sec = if sc > i64::from(MAX_SEC_ERR) {
            MAX_SEC_ERR
        } else {
            sc as u32
        };
        chain.qn = (k as i64 + chain_offset) as u32;
        ch[n_u] = chain;
        n_u += 1;
    }

    for i in 0..ch_n {
        let tlen = unitigs[sv[i].tn as usize].length;
        let (its, ite) = adjust_rev_tse(&sv[i], tlen);
        sv[i].ts = its;
        sv[i].te = ite;
        sv[i].tn = if p[i] >= 0 {
            (p[i] + chain_offset) as u32
        } else {
            u32::MAX
        };
    }
    let _ = qlen;
    n_u
}

/// `gen_linear_chains_ctg` (`inter.cpp` 16497): demote the packed target
/// ids, chain every unitig's alignments and keep the full-length ones.
#[allow(clippy::too_many_arguments)]
pub fn gen_linear_chains_ctg(
    res: &mut Vec<UlOv>,
    unitigs: &[MaUtg],
    bw: i64,
    diff_ec_ul: f64,
    qlen: i64,
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
) {
    res.sort_by_key(|item| item.tn);
    let mut sv = vec![UlOv::default(); res.len()];
    let mut m = 0usize;
    let n = res.len();
    let mut l = 0usize;
    while l < n {
        let mut k = l + 1;
        while k < n && res[k].tn == res[l].tn {
            k += 1;
        }
        for z in l..k {
            res[z].tn >>= 1;
        }
        let mut an = k;
        if k - l > 1 {
            an = l + linear_ctg_chain_dp_adv(
                &mut res[l..k],
                &mut sv[l..k],
                unitigs,
                bw,
                diff_ec_ul,
                qlen,
                max_skip,
                max_iter,
                max_dis,
                l as i64,
            );
        }
        for z in l..an {
            if is_available_ctg_aln(&unitigs[res[z].tn as usize], &res[z]) {
                res[z].ts = 0;
                res[z].te = unitigs[res[z].tn as usize].length;
                res[m] = res[z];
                m += 1;
            }
        }
        l = k;
    }
    res.truncate(m);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn packed_read(rid: u64, length: u32) -> u64 {
        (rid << 33) | u64::from(length)
    }

    #[test]
    fn gen_r_ug_idx_keeps_the_read_offset_and_base_position() {
        let unitigs = vec![MaUtg {
            length: 150,
            reads: vec![packed_read(0, 100), packed_read(1, 100)],
            ..MaUtg::default()
        }];
        let index = gen_r_ug_idx(&unitigs, 2, true);
        assert_eq!(index.idx, vec![0, 1, 2]);
        assert_eq!(index.p[0], UtgRidDt { u: 0, off: 0, pos: 0 });
        assert_eq!(index.p[1], UtgRidDt { u: 0, off: 1, pos: 100 });
    }

    #[test]
    fn ctg_rg2ug_merges_a_run_and_keeps_full_length_alignments() {
        let reads = [packed_read(0, 100), packed_read(1, 100)];
        let unitigs = vec![MaUtg {
            length: 150,
            reads: reads.to_vec(),
            ..MaUtg::default()
        }];
        let index = gen_r_ug_idx(&unitigs, 2, true);
        let contig = MaUtg {
            length: 200,
            reads: reads.to_vec(),
            ..MaUtg::default()
        };
        let mut u_cl = Vec::new();
        ctg_rg2ug_gen(&contig, &mut u_cl, &unitigs, &index, &[50, 50]);
        assert_eq!(u_cl.len(), 1);
        let z = u_cl[0];
        assert_eq!(
            (z.qn, z.qs, z.qe, z.tn, z.ts, z.te, z.sec, z.el, z.rev),
            (0, 0, 150, 0, 0, 150, 2, true, false)
        );

        let mut res = u_cl.clone();
        gen_linear_chains_ctg(&mut res, &unitigs, 1000, 0.02, 200, 25, 10_000, 1000);
        assert_eq!(res.len(), 1);
        assert_eq!((res[0].ts, res[0].te), (0, 150));
    }
}

// ---------------------------------------------------------------------------
// Graph chaining (`inter.cpp` 5080-6456, 13833-16570) used by
// `gen_contig_path` / `work_ctg_path_gchains`.
// ---------------------------------------------------------------------------

use std::collections::BTreeSet;
use std::collections::HashMap;

use crate::assembly::asg::Asg;
use crate::assembly::unitig::UnitigGraph;

/// `mg_lchain_t` (`Overlaps.h` 97).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MgLchain {
    pub off: i32,
    pub cnt: i32,
    pub inner_pre: bool,
    /// `ref_id|rev`.
    pub v: u32,
    pub rs: i32,
    pub re: i32,
    pub qs: i32,
    pub qe: i32,
    pub score: i32,
    pub dist_pre: i32,
    pub hash_pre: u32,
}

/// `mg_path_dst_t` (`inter.cpp` 288).
#[derive(Clone, Copy, Debug, Default)]
pub struct MgPathDst {
    pub v: u32,
    pub target_dist: i32,
    pub target_hash: u32,
    /// `meta:30 | check_hash:1 | inner:1`.
    pub meta: u32,
    pub check_hash: bool,
    pub inner: bool,
    pub qlen: i32,
    pub n_path: u32,
    pub is_0: bool,
    pub path_end: i32,
    pub dist: i32,
    pub mlen: i32,
    pub hash: u32,
}

/// `mg_pathv_t` (`inter.cpp` 191).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MgPathv {
    pub v: u32,
    pub d: u32,
    pub pre: i32,
}

const MG_MAX_SHORT_K: usize = 15;
const MG_SHORT_K_EXT: i32 = 10000;
const SEC_LEN_DIF: f64 = 0.03;
const GC_OFFSET_POS: i64 = 8;

/// `__ac_Wang_hash` (`khash.h` 406).
pub fn wang_hash(mut key: u32) -> u32 {
    key = key.wrapping_add(!(key << 15));
    key ^= key >> 10;
    key = key.wrapping_add(key << 3);
    key ^= key >> 6;
    key = key.wrapping_add(!(key << 11));
    key ^= key >> 16;
    key
}

/// `max_ovlp` (`inter.cpp` 1494).
pub fn max_ovlp(g: &Asg, v: u32) -> i64 {
    let mut o = 0u32;
    for arc in g.arcs_of(v) {
        if arc.deleted {
            continue;
        }
        if o < arc.overlap {
            o = arc.overlap;
        }
    }
    i64::from(o)
}

/// `g_adjacent_dis` (`inter.cpp` 4490).
pub fn g_adjacent_dis(g: &Asg, v: u32, w: u32) -> i64 {
    for arc in g.arcs_of(v) {
        if arc.deleted || arc.v != w {
            continue;
        }
        return i64::from(arc.node_length());
    }
    -1
}

/// `set_ul_ov_t_by_mg_lchain_t` (`inter.cpp` 5568).
pub fn set_ul_ov_t_by_mg_lchain(l: &MgLchain) -> UlOv {
    UlOv {
        tn: l.v >> 1,
        rev: (l.v & 1) != 0,
        ts: l.rs as u32,
        te: l.re as u32,
        qs: l.qs as u32,
        qe: l.qe as u32,
        ..UlOv::default()
    }
}

/// `get_offset_adjust` (`Overlaps.h` 1169).
pub(crate) fn get_offset_adjust(offset: u32, offset_len: u32, target_len: u32) -> u32 {
    (f64::from(offset) / f64::from(offset_len) * f64::from(target_len)) as u32
}

/// `infer_rovlp` (`Correct.cpp` 19083) with `bi = bj = NULL` and a unitig
/// graph (`ug != NULL`): the lengths come from the unitigs.
pub fn infer_rovlp(li: &UlOv, lj: &UlOv, unitigs: &[MaUtg]) -> i64 {
    infer_rovlp_with_lengths(
        li,
        lj,
        &|id| i64::from(unitigs[id as usize].length),
    )
}

/// `infer_rovlp` with `ug == NULL`: the lengths come from `R_INF` (the entry
/// `tn` is a read id).
pub fn infer_rovlp_reads(li: &UlOv, lj: &UlOv, read_lengths: &[u64]) -> i64 {
    infer_rovlp_with_lengths(li, lj, &|id| read_lengths[id as usize] as i64)
}

/// `infer_rovlp(li, lj, NULL, NULL, NULL, ug)`: the query overlap length with
/// the unitig-graph target lengths.
pub fn infer_rovlp_unitigs(li: &UlOv, lj: &UlOv, unitigs: &[crate::assembly::asg::MaUtg]) -> i64 {
    infer_rovlp_with_lengths(li, lj, &|id| unitigs[id as usize].length as i64)
}

fn infer_rovlp_with_lengths(li: &UlOv, lj: &UlOv, length: &dyn Fn(u32) -> i64) -> i64 {
    let in_len = length(li.tn);
    let mut is = i64::from(li.ts);
    let mut ie = i64::from(li.te);
    let irev = li.rev;
    let iqs = i64::from(li.qs);
    let iqe = i64::from(li.qe);

    let jn = length(lj.tn);
    let mut js = i64::from(lj.ts);
    let mut je = i64::from(lj.te);
    let jrev = lj.rev;
    let jqs = i64::from(lj.qs);
    let jqe = i64::from(lj.qe);

    let max_s = iqs.max(jqs);
    let min_e = iqe.min(jqe);
    if min_e <= max_s {
        return 0;
    }
    let mut s_shift = get_offset_adjust((max_s - iqs) as u32, (iqe - iqs) as u32, (ie - is) as u32);
    let mut e_shift = get_offset_adjust((iqe - min_e) as u32, (iqe - iqs) as u32, (ie - is) as u32);
    if irev {
        std::mem::swap(&mut s_shift, &mut e_shift);
    }
    is += i64::from(s_shift);
    ie -= i64::from(e_shift);

    let mut s_shift = get_offset_adjust((max_s - jqs) as u32, (jqe - jqs) as u32, (je - js) as u32);
    let mut e_shift = get_offset_adjust((jqe - min_e) as u32, (jqe - jqs) as u32, (je - js) as u32);
    if jrev {
        std::mem::swap(&mut s_shift, &mut e_shift);
    }
    js += i64::from(s_shift);
    je -= i64::from(e_shift);

    if irev {
        let ts = in_len - ie;
        let te = in_len - is;
        is = ts;
        ie = te;
    }
    if jrev {
        let ts = jn - je;
        let te = jn - js;
        js = ts;
        je = te;
    }

    if is <= js {
        js -= is;
        is = 0;
    } else {
        is -= js;
        js = 0;
    }

    let mut ir = in_len - ie;
    let jr = jn - je;
    if ir <= jr {
        ie = in_len;
        je += ir;
    } else {
        je = jn;
        ie += jr;
    }
    ir = ie - is;
    let jr = je - js;
    ir.max(jr)
}

/// `get_ecov_adv` (`inter.cpp` 3405) with `mode = 0` (unitig graph only).
pub fn get_ecov_adv(g: &Asg, v: u32, w: u32, bw: i64, diff_ec_ul: f64, dq: i64) -> bool {
    let mut dt: i64 = -1;
    for arc in g.arcs_of(v) {
        if arc.deleted || arc.v != w {
            continue;
        }
        dt = i64::from(arc.overlap);
        break;
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

/// `find_mg_lchain_max` (`inter.cpp` 5578).
fn find_mg_lchain_max(a: &[MgLchain], x: i32) -> i64 {
    let n = a.len();
    if n == 0 {
        return -1;
    }
    if a[n - 1].qe < x {
        return n as i64 - 1;
    }
    if a[0].qe >= x {
        return -1;
    }
    let mut s = 0usize;
    let mut e = n;
    while e > s {
        let m = s + (e - s) / 2;
        if a[m].qe >= x {
            e = m;
        } else {
            s = m + 1;
        }
    }
    s as i64
}

/// `extend_end_coord` (`inter.cpp` 13793) on a `mg_lchain_t`.
pub fn extend_end_coord_lchain(
    li: &MgLchain,
    qlen: i64,
    rlen: i64,
) -> (i64, i64, i64, i64) {
    let mut qs = i64::from(li.qs);
    let qe = i64::from(li.qe);
    let mut rs = i64::from(li.rs);
    let mut re = i64::from(li.re);
    let rev = (li.v & 1) != 0;
    if rev {
        rs = rlen - i64::from(li.re);
        re = rlen - i64::from(li.rs);
    }
    let (mut qs, qe, mut rs, mut re) = (qs, qe, rs, re);
    if qs <= rs {
        rs -= qs;
        qs = 0;
    } else {
        qs -= rs;
        rs = 0;
    }
    let qtail = qlen - qe;
    let rtail = rlen - re;
    let (mut qe, mut re) = (qe, re);
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

/// One node of the shortest-path search (`sp_node_t`, `inter.cpp` 310).
#[derive(Clone, Copy, Debug, Default)]
struct SpNode {
    v: u32,
    /// `dist << 32 | id`; afterwards the low half is the index in `out`.
    di: u64,
    pre: i64,
    hash: u32,
    is_0: bool,
}

/// `sp_topk_t` (`inter.cpp` 314): the max-heap of paths reaching one vertex.
#[derive(Clone, Debug, Default)]
struct SpTopk {
    p: Vec<usize>,
    qs: i32,
    qe: i32,
    mlen: i32,
}

fn heap_less(nodes: &[SpNode], a: usize, b: usize) -> bool {
    nodes[a].di < nodes[b].di
}

fn heap_up(nodes: &[SpNode], heap: &mut [usize], mut i: usize) {
    while i > 0 {
        let parent = (i - 1) / 2;
        if heap_less(nodes, heap[i], heap[parent]) {
            heap.swap(i, parent);
            i = parent;
        } else {
            break;
        }
    }
}

fn heap_down(nodes: &[SpNode], heap: &mut [usize], mut i: usize) {
    let n = heap.len();
    loop {
        let left = i * 2 + 1;
        if left >= n {
            break;
        }
        let right = left + 1;
        let mut best = left;
        if right < n && heap_less(nodes, heap[right], heap[left]) {
            best = right;
        }
        if heap_less(nodes, heap[best], heap[i]) {
            heap.swap(i, best);
            i = best;
        } else {
            break;
        }
    }
}

fn push_sp_node(
    nodes: &mut Vec<SpNode>,
    queue: &mut BTreeSet<(u64, usize)>,
    next_id: &mut u64,
    v: u32,
    dist: i32,
    pre: i64,
    hash: u32,
    is_0: bool,
) -> usize {
    let id = *next_id;
    *next_id += 1;
    let di = ((dist as u32 as u64) << 32) | id;
    nodes.push(SpNode { v, di, pre, hash, is_0 });
    let index = nodes.len() - 1;
    queue.insert((di, index));
    index
}

/// `hc_shortest_k` (`inter.cpp` 14036): the bounded top-k shortest-path search
/// used by the contig chainer. Returns nothing; `dst`, `out`, `res` and
/// `dst_done` carry the results.
#[allow(clippy::too_many_arguments)]
pub fn hc_shortest_k(
    g: &Asg,
    src: u32,
    dst: &mut [MgPathDst],
    max_dist: i32,
    max_k: i32,
    dst_done: &mut Vec<u64>,
    out: &mut Vec<usize>,
    res: &mut Vec<MgPathv>,
    mut first_src_ban: bool,
    detect_mul_way: bool,
    len_dif: f64,
) {
    if dst.is_empty() {
        return;
    }
    let n_dst = dst.len();
    for t in dst.iter_mut() {
        if t.inner {
            t.dist = 0;
            t.n_path = 1;
            t.path_end = -1;
        } else {
            t.dist = -1;
            t.n_path = 0;
            t.path_end = -1;
        }
    }
    let mut max_k = max_k as usize;
    if max_k > MG_MAX_SHORT_K {
        max_k = MG_MAX_SHORT_K;
    }

    dst_done.clear();
    dst_done.resize(n_dst, 0);
    let mut dst_group: Vec<u64> = (0..n_dst)
        .map(|i| (u64::from(dst[i].v) << 32) | i as u64)
        .collect();
    dst_group.sort();

    let mut h2: HashMap<u32, Vec<usize>> = HashMap::new();
    for &entry in &dst_group {
        h2.entry((entry >> 32) as u32).or_default().push((entry as u32) as usize);
    }

    out.clear();
    let mut nodes: Vec<SpNode> = Vec::new();
    let mut queue: BTreeSet<(u64, usize)> = BTreeSet::new();
    let mut topk: HashMap<u32, SpTopk> = HashMap::new();
    let mut next_id: u64 = 0;

    let src_node = push_sp_node(
        &mut nodes,
        &mut queue,
        &mut next_id,
        src,
        0,
        -1,
        wang_hash(src),
        true,
    );
    topk.insert(
        src,
        SpTopk {
            p: vec![src_node],
            qs: -1,
            qe: -1,
            mlen: 0,
        },
    );

    while let Some(&(_, first)) = queue.iter().next() {
        queue.remove(&(nodes[first].di, first));
        let out_index = out.len() as u64;
        nodes[first].di = (nodes[first].di >> 32 << 32) | out_index;
        out.push(first);

        let vertex = nodes[first].v;
        let dist = (nodes[first].di >> 32) as i32;
        if first_src_ban {
            if let Some(group) = h2.get(&vertex) {
                for &dst_index in group {
                    let t = &mut dst[dst_index];
                    let mut done = false;
                    if t.inner {
                        done = true;
                    } else {
                        let mlen = 0i32;
                        let mut copy = false;
                        if t.n_path == 0 {
                            copy = true;
                        } else if t.target_dist >= 0 {
                            if dist == t.target_dist
                                && t.check_hash
                                && nodes[first].hash == t.target_hash
                            {
                                copy = true;
                                done = true;
                            } else {
                                let d0 = (t.dist - t.target_dist).abs();
                                let d1 = (dist - t.target_dist).abs();
                                if d1 - mlen / 2 < d0 - t.mlen / 2 {
                                    copy = true;
                                }
                            }
                        }
                        if copy {
                            t.path_end = out.len() as i32 - 1;
                            t.dist = dist;
                            t.hash = nodes[first].hash;
                            t.mlen = mlen;
                            t.is_0 = nodes[first].is_0;
                            if t.target_dist >= 0 {
                                if dist == t.target_dist
                                    && t.check_hash
                                    && nodes[first].hash == t.target_hash
                                {
                                    done = true;
                                } else if dist > t.target_dist + MG_SHORT_K_EXT
                                    && dist > (t.target_dist >> 4)
                                    && f64::from(dist) > f64::from(t.target_dist) * 1.25
                                {
                                    done = true;
                                }
                            }
                        }
                        t.n_path += 1;
                        if t.n_path as usize >= max_k {
                            done = true;
                        }
                    }
                    if !detect_mul_way && dst_done[dst_index] == 0 && done {
                        dst_done[dst_index] = 1;
                    }
                }
            }
        }
        first_src_ban = true;

        for arc in g.arcs_of(vertex) {
            let d = dist.wrapping_add(arc.node_length() as i32);
            if d > max_dist {
                continue;
            }
            let entry = topk.entry(arc.v).or_default();
            if entry.p.len() < max_k {
                let child_hash = nodes[first].hash.wrapping_add(wang_hash(arc.v));
                let child_is_0 = nodes[first].is_0;
                let node = push_sp_node(
                    &mut nodes,
                    &mut queue,
                    &mut next_id,
                    arc.v,
                    d,
                    out.len() as i64 - 1,
                    child_hash,
                    child_is_0,
                );
                entry.p.push(node);
                let len = entry.p.len();
                heap_up(&nodes, &mut entry.p, len - 1);
            } else if (nodes[entry.p[0]].di >> 32) as i32 > d {
                let victim = entry.p[0];
                queue.remove(&(nodes[victim].di, victim));
                let id = next_id;
                next_id += 1;
                let child_hash = nodes[first].hash.wrapping_add(wang_hash(arc.v));
                let child_is_0 = nodes[first].is_0;
                nodes[victim].di = ((d as u32 as u64) << 32) | id;
                nodes[victim].pre = out.len() as i64 - 1;
                nodes[victim].hash = child_hash;
                nodes[victim].is_0 = child_is_0;
                queue.insert((nodes[victim].di, victim));
                let len = entry.p.len();
                heap_down(&nodes, &mut entry.p, 0);
                let _ = len;
            }
        }
    }

    let n_found = dst.iter().filter(|t| t.n_path > 0).count();
    if n_found == 0 {
        return;
    }

    let mut trans = vec![0u64; out.len()];
    let res_start = res.len();
    let mut n_mpath = 1usize;
    if detect_mul_way {
        n_mpath = phase_mul_ways(
            res, &mut trans, &nodes, out, dst, max_k, len_dif, res_start,
        );
    }
    if n_mpath == 1 {
        for t in dst.iter() {
            if t.n_path > 0 && t.target_dist >= 0 && t.path_end >= 0 {
                trans[(nodes[out[t.path_end as usize]].di & 0xffff_ffff) as usize] = 1;
            }
        }
        for i in (0..out.len()).rev() {
            if trans[i] != 0 {
                let pre = nodes[out[i]].pre;
                if pre >= 0 {
                    trans[pre as usize] = 1;
                }
            }
        }
        let mut n = 0u64;
        for value in trans.iter_mut() {
            if *value != 0 {
                *value = n;
                n += 1;
            } else {
                *value = u64::MAX;
            }
        }
        res.resize(res_start + n as usize, MgPathv::default());
        for (i, &node) in out.iter().enumerate() {
            if trans[i] == u64::MAX {
                continue;
            }
            let index = res_start + trans[i] as usize;
            let pre = nodes[node].pre;
            res[index] = MgPathv {
                v: nodes[node].v,
                d: (nodes[node].di >> 32) as u32,
                pre: if pre < 0 {
                    pre as i32
                } else {
                    trans[pre as usize] as i32
                },
            };
        }
        for t in dst.iter_mut() {
            if t.path_end >= 0 {
                t.path_end = trans[t.path_end as usize] as i32;
            }
        }
    }
}

/// `select_mul_way_nodes` (`inter.cpp` 13931).
fn select_mul_way_nodes(
    a: &mut [MgPathv],
    nodes: &[SpNode],
    out: &[usize],
    len_dif: f64,
    m_pathn: usize,
    flag: &mut [u64],
) -> usize {
    let a_n = a.len();
    if a_n > 1 {
        a.sort_by_key(|item| item.d);
    }
    let mut tt = 0;
    for i in 0..a_n {
        let p = a[i];
        let pp = (p.d << 1) >> 1;
        let pd = (nodes[out[p.pre as usize]].di >> 32) as i64;
        if p.d & 0x8000_0000 != 0 {
            break;
        }
        let mut occ = 1usize;
        for (k, item) in a.iter().enumerate() {
            if k == i {
                continue;
            }
            if ((item.d << 1) >> 1) == pp {
                continue;
            }
            if item.v != p.v {
                continue;
            }
            let kd = (nodes[out[item.pre as usize]].di >> 32) as i64;
            if (kd as f64) >= pd as f64 * (1.0 - len_dif)
                && (kd as f64) <= pd as f64 * (1.0 + len_dif)
            {
                occ += 1;
            }
        }
        if occ >= m_pathn {
            flag[p.pre as usize] = 1;
            tt += 1;
        }
    }
    tt
}

/// `phase_mul_ways` (`inter.cpp` 13911).
#[allow(clippy::too_many_arguments)]
fn phase_mul_ways(
    res: &mut Vec<MgPathv>,
    trans: &mut [u64],
    nodes: &[SpNode],
    out: &[usize],
    dst: &mut [MgPathDst],
    max_k: usize,
    len_dif: f64,
    res_start: usize,
) -> usize {
    let mut n = 0usize;
    let mut n_mpath = 0usize;

    for t in dst.iter_mut() {
        if t.n_path > 0 && t.target_dist >= 0 && t.path_end >= 0 {
            if t.n_path as usize >= max_k {
                t.n_path = 0;
                let mut z = t.path_end;
                let mut zl;
                loop {
                    zl = z;
                    if z < 0 {
                        break;
                    }
                    z = nodes[out[z as usize]].pre as i32;
                }
                n += 2;
                trans[t.path_end as usize] = 1;
                if zl >= 0 {
                    trans[zl as usize] = 1;
                }
            } else {
                let target = f64::from(t.target_dist);
                let lo = target * (1.0 - len_dif);
                let hi = target * (1.0 + len_dif);
                let mut kk_p = 0i64;
                let mut j = t.path_end as i64;
                while j < out.len() as i64 {
                    let od = (nodes[out[j as usize]].di >> 32) as i64 as f64;
                    if od < lo || od > hi {
                        break;
                    }
                    if nodes[out[j as usize]].v == t.v {
                        kk_p += 1;
                    }
                    j += 1;
                }
                let mut j = t.path_end as i64 - 1;
                while j >= 0 {
                    let od = (nodes[out[j as usize]].di >> 32) as i64 as f64;
                    if od < lo || od > hi {
                        break;
                    }
                    if nodes[out[j as usize]].v == t.v {
                        kk_p += 1;
                    }
                    j -= 1;
                }
                n_mpath += kk_p as usize;
            }
        }
    }

    if n_mpath > 1 {
        let mut pid = 0u32;
        for t in dst.iter() {
            if !(t.n_path > 0 && t.target_dist >= 0 && t.path_end >= 0) {
                continue;
            }
            let target = f64::from(t.target_dist);
            let lo = target * (1.0 - len_dif);
            let hi = target * (1.0 + len_dif);
            let mut j = t.path_end as i64;
            while j < out.len() as i64 {
                let od = (nodes[out[j as usize]].di >> 32) as i64 as f64;
                if od < lo || od > hi {
                    break;
                }
                if nodes[out[j as usize]].v == t.v {
                    let mut z = j;
                    while z >= 0 {
                        let mut d = pid;
                        if j != t.path_end as i64 {
                            d |= 0x8000_0000;
                        }
                        res.push(MgPathv {
                            v: nodes[out[z as usize]].v,
                            pre: z as i32,
                            d,
                        });
                        z = nodes[out[z as usize]].pre;
                    }
                    pid += 1;
                }
                j += 1;
            }
            let mut j = t.path_end as i64 - 1;
            while j >= 0 {
                let od = (nodes[out[j as usize]].di >> 32) as i64 as f64;
                if od < lo || od > hi {
                    break;
                }
                if nodes[out[j as usize]].v == t.v {
                    let mut z = j;
                    while z >= 0 {
                        let mut d = pid;
                        if j != t.path_end as i64 {
                            d |= 0x8000_0000;
                        }
                        res.push(MgPathv {
                            v: nodes[out[z as usize]].v,
                            pre: z as i32,
                            d,
                        });
                        z = nodes[out[z as usize]].pre;
                    }
                    pid += 1;
                }
                j -= 1;
            }
        }

        res[res_start..].sort_by_key(|item| item.v);
        let res_end = res.len();
        let mut j = res_start;
        let mut i = res_start + 1;
        while i <= res_end {
            if i == res_end || res[i].v != res[j].v {
                n += select_mul_way_nodes(
                    &mut res[j..i],
                    nodes,
                    out,
                    len_dif,
                    n_mpath,
                    trans,
                );
                j = i;
            }
            i += 1;
        }
        res.truncate(res_start);
    }

    if n > 0 {
        let mut count = 0u64;
        for value in trans.iter_mut() {
            if *value != 0 {
                *value = count;
                count += 1;
            } else {
                *value = u64::MAX;
            }
        }
        res.resize(res_start + count as usize, MgPathv::default());
        for (i, &node) in out.iter().enumerate() {
            if trans[i] == u64::MAX {
                continue;
            }
            let index = res_start + trans[i] as usize;
            let pre = nodes[node].pre;
            res[index] = MgPathv {
                v: nodes[node].v,
                d: (nodes[node].di >> 32) as u32,
                pre: if pre < 0 {
                    pre as i32
                } else if trans[pre as usize] == u64::MAX {
                    -2
                } else {
                    trans[pre as usize] as i32
                },
            };
        }
        for t in dst.iter_mut() {
            if t.path_end >= 0 {
                t.path_end = trans[t.path_end as usize] as i32;
            }
        }
    }
    n_mpath
}

/// `reverse_track` (`inter.cpp` 14640).
fn reverse_track(a: &mut [MgPathv]) {
    let a_n = a.len();
    if a_n <= 1 {
        return;
    }
    for k in (1..a_n).rev() {
        a[k].d -= a[k - 1].d;
    }
    let hn = a_n / 2;
    for k in 0..hn {
        a.swap(k, a_n - k - 1);
        a[k].v ^= 1;
        a[a_n - k - 1].v ^= 1;
    }
    if a_n & 1 != 0 {
        a[hn].v ^= 1;
    }
}

/// `gen_gchain_track` (`inter.cpp` 14670).
pub fn gen_gchain_track(
    g: &Asg,
    a: &[MgLchain],
    dst_done: &mut Vec<u64>,
    out: &mut Vec<usize>,
    res: &mut Vec<MgPathv>,
) -> usize {
    res.clear();
    res.push(MgPathv {
        v: u32::MAX,
        d: 0,
        pre: 0,
    });
    for k in 1..a.len() {
        let l0 = a[k - 1];
        let l1 = a[k];
        debug_assert!(!l1.inner_pre);
        debug_assert!(l1.dist_pre >= 0);
        let mut dst = MgPathDst {
            v: l0.v ^ 1,
            target_dist: l1.dist_pre,
            target_hash: l1.hash_pre,
            check_hash: true,
            ..MgPathDst::default()
        };
        let p_n = res.len();
        if dst.target_hash != wang_hash(l1.v ^ 1).wrapping_add(wang_hash(dst.v))
            || g_adjacent_dis(g, l1.v ^ 1, dst.v) != i64::from(dst.target_dist)
        {
            let max_dist = (f64::from(dst.target_dist) * (1.0 + SEC_LEN_DIF)) as i32;
            hc_shortest_k(
                g,
                l1.v ^ 1,
                std::slice::from_mut(&mut dst),
                max_dist,
                MG_MAX_SHORT_K as i32,
                dst_done,
                out,
                res,
                false,
                true,
                SEC_LEN_DIF,
            );
            debug_assert!(res.len() - p_n > 1);
            let last = res.len() - 1;
            res[last].d = res[last].d.wrapping_sub(res[last - 1].d);
            res.pop();
            reverse_track(&mut res[p_n..]);
            res.pop();
        } else {
            res[p_n - 1].d = dst.target_dist as u32;
        }
        res.push(MgPathv {
            v: u32::MAX,
            pre: k as i32,
            d: 0,
        });
    }
    res.len()
}

/// `gl_chain_linear` (`inter.cpp` 5654): linear chaining over the
/// linear-chain list, returning the best chain index.
#[allow(clippy::too_many_arguments)]
pub fn gl_chain_linear(
    g: &Asg,
    lc: &mut [MgLchain],
    unitigs: &[MaUtg],
    qlen: i64,
    bw: i64,
    ng_diff_thre: f64,
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
    n_ext: usize,
) -> (i64, Vec<u64>) {
    let lc_n = lc.len();
    let mut t = vec![0i64; lc_n];
    let mut f = vec![0i64; lc_n];
    let mut p = vec![-1i64; lc_n];

    let mut st: i64 = 0;
    let mut max_ii: i64 = -1;
    for i in 0..n_ext {
        let li = lc[i];
        let ui = set_ul_ov_t_by_mg_lchain(&li);
        let mm_ovlp = max_ovlp(g, li.v ^ 1);
        let mut x = ((i64::from(li.qs) + mm_ovlp) as f64 * ng_diff_thre) as i64;
        if x < bw {
            x = bw;
        }
        x += i64::from(li.qs) + mm_ovlp;
        if x > qlen + 1 {
            x = qlen + 1;
        }
        let x = find_mg_lchain_max(&lc[..i], (x + 128) as i32);
        let csc = i64::from(li.score);
        let mut mm_sc = csc;
        let mut mm_idx: i64 = -1;
        let mut n_skip: i64 = 0;
        if (x - st) > max_iter {
            st = x - max_iter;
        }
        let mut j = x;
        while j >= st {
            let lj = lc[j as usize];
            if i64::from(lj.qe) + 128 <= i64::from(li.qs) {
                break;
            }
            if lj.qs >= li.qs {
                j -= 1;
                continue;
            }
            let uj = set_ul_ov_t_by_mg_lchain(&lj);
            let qo = infer_rovlp(&ui, &uj, unitigs);
            if get_ecov_adv(g, li.v ^ 1, lj.v ^ 1, bw, ng_diff_thre, qo) {
                let sc = csc + f[j as usize];
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
                if p[j as usize] >= 0 {
                    t[p[j as usize] as usize] = i as i64;
                }
            }
            j -= 1;
        }
        let end_j = j;

        if max_ii < 0
            || i64::from(lc[i].qe) > i64::from(lc[max_ii as usize].qe) + max_dis
        {
            let mut max = i64::from(i32::MIN);
            max_ii = -1;
            let mut j = i as i64 - 1;
            while j >= st
                && i64::from(lc[i].qe) <= max_dis + i64::from(lc[j as usize].qe)
            {
                if max < f[j as usize] {
                    max = f[j as usize];
                    max_ii = j;
                }
                j -= 1;
            }
        }

        if max_ii >= 0 && max_ii < end_j {
            let lj = lc[max_ii as usize];
            if i64::from(lj.qe) + 128 > i64::from(li.qs) && lj.qs < li.qs {
                let uj = set_ul_ov_t_by_mg_lchain(&lj);
                let qo = infer_rovlp(&ui, &uj, unitigs);
                if get_ecov_adv(g, li.v ^ 1, lj.v ^ 1, bw, ng_diff_thre, qo) {
                    let sc = csc + f[max_ii as usize];
                    if sc > mm_sc {
                        mm_sc = sc;
                        mm_idx = max_ii;
                    }
                }
            }
        }

        f[i] = mm_sc;
        p[i] = mm_idx;
        if max_ii < 0
            || (i64::from(lc[i].qe) <= max_dis + i64::from(lc[max_ii as usize].qe)
                && f[max_ii as usize] < f[i])
        {
            max_ii = i as i64;
        }
        lc[i].dist_pre = if mm_idx < 0 {
            -1
        } else {
            g_adjacent_dis(g, li.v ^ 1, lc[mm_idx as usize].v ^ 1) as i32
        };
        lc[i].inner_pre = false;
        lc[i].hash_pre = if mm_idx < 0 {
            0
        } else {
            wang_hash(li.v ^ 1).wrapping_add(wang_hash(lc[mm_idx as usize].v ^ 1))
        };
    }

    for i in n_ext..lc_n {
        let max_f = i64::from(lc[i].score);
        f[i] = max_f;
        p[i] = -1;
        lc[i].dist_pre = -1;
        lc[i].hash_pre = 0;
        lc[i].inner_pre = false;
    }

    let mut srt: Vec<u64> = (0..lc_n)
        .map(|i| ((f[i] as u64) << 32) | ((i as u64) << 1))
        .collect();
    srt.sort();

    let mut sw: Vec<MgLchain> = Vec::with_capacity(lc_n);
    let mut u: Vec<u64> = Vec::with_capacity(lc_n);
    let mut n_v = 0usize;
    for k in (0..lc_n).rev() {
        let n_v0 = n_v;
        let mut i = ((srt[k] as u32) >> 1) as i64;
        while i >= 0 && (srt[i as usize] & 1) == 0 {
            sw.push(lc[i as usize]);
            n_v += 1;
            srt[i as usize] |= 1;
            i = p[i as usize];
        }
        if n_v == n_v0 {
            continue;
        }
        let mut sc = if i < 0 {
            (srt[k] >> 32) as i64
        } else {
            (srt[k] >> 32) as i64 - f[i as usize]
        };
        let ff = if sc >= 0 {
            0x8000_0000_0000_0000u64
        } else {
            sc = -sc;
            0
        };
        u.push((((sc as u64) << 32) | (n_v - n_v0) as u64) | ff);
    }

    let mut m_idx: i64 = -1;
    let mut m_sc: i64 = -1;
    let mut k = 0usize;
    for i in 0..u.len() {
        let saved = u[i];
        let sc;
        if saved & 0x8000_0000_0000_0000 != 0 {
            u[i] -= 0x8000_0000_0000_0000;
            sc = (u[i] >> 32) as i64;
        } else {
            sc = -((saved >> 32) as i64);
        }
        u[i] = (u[i] << 32) >> 32;
        u[i] |= (sc as u64) << 32;

        let k0 = k;
        let ni = u[i] as u32 as usize;
        for j in 0..ni {
            lc[k] = sw[k0 + (ni - j - 1)];
            k += 1;
        }
        if m_idx < 0 || m_sc < (u[i] >> 32) as i64 {
            m_idx = i as i64;
            m_sc = (u[i] >> 32) as i64;
        }
    }
    debug_assert_eq!(k, n_v);
    (m_idx, u)
}

/// `ctg_chain_graph` (`inter.cpp` 5989).
#[allow(clippy::too_many_arguments)]
pub fn ctg_chain_graph(
    g: &Asg,
    unitigs: &[MaUtg],
    lc: &mut Vec<MgLchain>,
    qlen: i64,
    bw: i64,
    ng_diff_thre: f64,
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
) -> (i64, Vec<u64>) {
    if lc.is_empty() {
        return (0, Vec::new());
    }
    let lc_n = lc.len();
    let mut srt: Vec<u64> = Vec::with_capacity(lc_n);
    for i in 0..lc_n {
        lc[i].dist_pre = -1;
        srt.push(((lc[i].qe as u32 as u64) << 32) | i as u64);
    }
    srt.sort();
    let mut j = 0usize;
    for i in 1..=lc_n {
        if i == lc_n || (srt[i] >> 32) != (srt[j] >> 32) {
            if i - j > 1 {
                for x in j..i {
                    srt[x] = (srt[x] << 32 >> 32) | ((lc[(srt[x] as u32) as usize].qs as u32 as u64) << 32);
                }
                srt[j..i].sort();
            }
            j = i;
        }
    }
    let mut sw: Vec<MgLchain> = (0..lc_n).map(|i| lc[(srt[i] as u32) as usize]).collect();
    lc.clear();
    lc.extend_from_slice(&sw);
    let _ = &mut sw;

    let n_ext = lc_n;
    gl_chain_linear(
        g, lc, unitigs, qlen, bw, ng_diff_thre, max_skip, max_iter, max_dis, n_ext,
    )
}

/// `is_compat_contig_chain` (`inter.cpp` 6438).
fn is_compat_contig_chain(g: &Asg, p: &UlOv, m: &UlOv, a: &[MgLchain]) -> bool {
    let ovlp = if m.qe.min(p.qe) > m.qs.max(p.qs) {
        u64::from(m.qe.min(p.qe) - m.qs.max(p.qs))
    } else {
        0
    };
    if ovlp == 0 {
        return true;
    }
    if ovlp == u64::from(m.qe - m.qs) {
        return false;
    }
    if ovlp == u64::from(p.qe - p.qs) {
        return false;
    }
    let (i0, i1) = if p.qs <= m.qs { (p, m) } else { (m, p) };
    let o0 = max_ovlp(g, a[i0.te as usize - 1].v);
    let o1 = max_ovlp(g, a[i1.ts as usize].v ^ 1);
    let cut = (o0.min(o1) as f64 * 1.333333) as u64;
    ovlp <= cut
}

/// `select_max_ctg_chain` (`inter.cpp` 6456). `raw_idx` is debug-only in the
/// reference and is therefore omitted.
#[allow(clippy::too_many_arguments)]
pub fn select_max_ctg_chain(
    g: &Asg,
    idx: &[u64],
    e: &[MgLchain],
    dst_done: &mut Vec<u64>,
    out: &mut Vec<usize>,
    res: &mut Vec<MgPathv>,
    gchains: &mut Vec<MgLchain>,
) -> usize {
    gchains.clear();
    if idx.is_empty() {
        return 0;
    }
    let mut gb: Vec<UlOv> = Vec::with_capacity(idx.len());
    let mut a_n = 0usize;
    for &entry in idx {
        let count = (entry as u32) as usize;
        gb.push(UlOv {
            qn: (entry >> 32) as u32,
            ts: a_n as u32,
            te: (a_n + count) as u32,
            qs: e[a_n].qs as u32,
            qe: e[a_n + count - 1].qe as u32,
            tn: 1,
            ..UlOv::default()
        });
        a_n += count;
    }
    gb.sort_by_key(|item| item.qn);
    let gb_n = gb.len();
    for k in 0..gb_n / 2 {
        gb.swap(k, gb_n - k - 1);
    }

    for k in 0..gb_n {
        let om = u64::from(gb[k].qe - gb[k].qs);
        if gb[k].tn == 0 {
            continue;
        }
        let mut i = k as i64 - 1;
        while i >= 0 {
            let p = gb[i as usize];
            let m = gb[k];
            let ovlp = if m.qe.min(p.qe) > m.qs.max(p.qs) {
                u64::from(m.qe.min(p.qe) - m.qs.max(p.qs))
            } else {
                0
            };
            if ovlp == 0 {
                i -= 1;
                continue;
            }
            let min_sc = u64::from(p.qn.min(m.qn));
            let max_sc = u64::from(p.qn.max(m.qn));
            let mut ok = u64::from(p.qe - p.qs);
            ok = ok.max(om);
            if (min_sc as f64) < max_sc as f64 * 0.98 {
                break;
            }
            if ovlp > GC_OFFSET_POS as u64
                && (min_sc as f64) > max_sc as f64 * 0.98
                && (ovlp as f64) > ok as f64 * 0.8
            {
                gb[k].tn = 0;
                gb[i as usize].tn = 0;
            }
            i -= 1;
        }
        let mut i = k + 1;
        while i < gb_n {
            let p = gb[i];
            let m = gb[k];
            let ovlp = if m.qe.min(p.qe) > m.qs.max(p.qs) {
                u64::from(m.qe.min(p.qe) - m.qs.max(p.qs))
            } else {
                0
            };
            if ovlp == 0 {
                i += 1;
                continue;
            }
            let min_sc = u64::from(p.qn.min(m.qn));
            let max_sc = u64::from(p.qn.max(m.qn));
            let mut ok = u64::from(p.qe - p.qs);
            ok = ok.max(om);
            if (min_sc as f64) < max_sc as f64 * 0.98 {
                break;
            }
            if ovlp > GC_OFFSET_POS as u64
                && (min_sc as f64) > max_sc as f64 * 0.98
                && (ovlp as f64) > ok as f64 * 0.8
            {
                gb[k].tn = 0;
                gb[i].tn = 0;
            }
            i += 1;
        }
    }
    gb.retain(|item| item.tn != 0);

    let mut accepted: Vec<UlOv> = Vec::new();
    for m in gb {
        let mut compat = true;
        for p in &accepted {
            if !is_compat_contig_chain(g, p, &m, e) {
                compat = false;
                break;
            }
        }
        if compat {
            accepted.push(m);
        }
    }

    let mut n_mchain = 0usize;
    if !accepted.is_empty() {
        gchains.clear();
        for m in accepted {
            // `kv_pushp`: the chain header is appended before the path is
            // generated and later carries the path length in `cnt`.
            let header = gchains.len();
            gchains.push(MgLchain {
                v: u32::MAX,
                qs: m.qs as i32,
                qe: m.qe as i32,
                rs: m.ts as i32,
                re: m.te as i32,
                cnt: 0,
                off: 0,
                ..MgLchain::default()
            });
            let rs = m.ts as usize;
            let re = m.te as usize;
            gen_gchain_track(g, &e[rs..re], dst_done, out, res);
            let res_n = res.len();
            let base = gchains.len();
            gchains.resize(base + res_n, MgLchain::default());
            for (i, p) in res.iter().enumerate() {
                let target = &mut gchains[base + i];
                if p.v == u32::MAX {
                    *target = e[p.pre as usize + rs];
                    target.dist_pre = p.d as i32;
                } else {
                    target.v = p.v;
                    target.off = -1;
                    target.dist_pre = p.d as i32;
                    target.rs = 0;
                    target.re = unitig_length_for(g, p.v >> 1);
                }
            }
            gchains[header].cnt = res_n as i32;
            n_mchain += 1;
        }
    }
    n_mchain
}

fn unitig_length_for(g: &Asg, id: u32) -> i32 {
    g.sequences[id as usize].length as i32
}

/// `dump_linear_chain` (`inter.cpp` 13833).
pub fn dump_linear_chain(
    unitigs: &[MaUtg],
    lidx: &[UlOv],
    qlen: i64,
) -> Vec<MgLchain> {
    let mut res = Vec::with_capacity(lidx.len());
    for (i, item) in lidx.iter().enumerate() {
        let mut p = MgLchain {
            v: (item.tn << 1) | u32::from(item.rev),
            off: i as i32,
            score: item.sec as i32,
            qs: item.qs as i32,
            qe: item.qe as i32,
            rs: item.ts as i32,
            re: item.te as i32,
            ..MgLchain::default()
        };
        let rlen = i64::from(unitigs[(p.v >> 1) as usize].length);
        let (qs, qe, rs, re) = extend_end_coord_lchain(&p, qlen, rlen);
        p.qs = qs as i32;
        p.qe = qe as i32;
        p.rs = rs as i32;
        p.re = re as i32;
        res.push(p);
    }
    res
}

/// `push_ctg_res` (`inter.cpp` 16055).
pub fn push_ctg_res(res: &mut UlVec, uc: &[MgLchain]) {
    res.bb.clear();
    let ucn = uc.len();
    let mut tt = 0usize;
    let mut k = 0usize;
    while k < ucn {
        let ix = uc[k];
        debug_assert_eq!(ix.v, u32::MAX);
        res.bb.push(UcBlock {
            hid: k as u32,
            qs: ix.qs as u32,
            qe: ix.qe as u32,
            ts: ix.cnt as u32,
            ..UcBlock::default()
        });
        tt += ix.cnt as usize;
        k += ix.cnt as usize + 1;
    }
    res.bb.sort_by_key(|item| item.qe);

    let header_n = res.bb.len();
    res.bb.resize(header_n + tt, UcBlock::default());
    let mut tt = header_n;
    for k in 0..header_n {
        let s = tt;
        let e = tt + res.bb[k].ts as usize;
        res.bb[k].ts = s as u32;
        res.bb[k].te = e as u32;
        tt = e;
        let hid = res.bb[k].hid as usize;
        let ix = uc[hid];
        debug_assert_eq!(ix.v, u32::MAX);
        let rrn = (res.bb[k].te - res.bb[k].ts) as usize;
        debug_assert_eq!(ix.cnt as usize, rrn);
        for z in 0..rrn {
            let src = uc[hid + 1 + z];
            res.bb[s + z] = UcBlock {
                hid: src.v >> 1,
                rev: (src.v & 1) != 0,
                qs: src.qs as u32,
                qe: src.qe as u32,
                ts: src.rs as u32,
                te: src.re as u32,
                aidx: hid as u32,
                ..UcBlock::default()
            };
        }
        res.bb[k].hid = u32::MAX;
    }
    debug_assert!(tt <= res.bb.len());
    // `rch->bb.n` stays at the header count in the reference; the data
    // blocks live in the capacity-only tail.
    res.n = header_n;
}

/// `direct_gchain_scaf` (`inter.cpp` 16528) followed by
/// `work_ctg_path_gchains`'s per-contig body (`inter.cpp` 18880).
pub fn gen_contig_path_for_contig(
    contig: &MaUtg,
    unitigs: &[MaUtg],
    index: &UtgRid,
    read_lengths: &[u64],
    ug_graph: &Asg,
) -> UlVec {
    let mut res = UlVec::default();
    let mut u_cl: Vec<UlOv> = Vec::new();
    ctg_rg2ug_gen(contig, &mut u_cl, unitigs, index, read_lengths);
    if u_cl.is_empty() {
        return res;
    }
    gen_linear_chains_ctg(
        &mut u_cl,
        unitigs,
        16,
        UL_ERROR_RATE,
        i64::from(contig.length),
        UG_SKIP_N,
        UG_ITER_N,
        UG_DIS_N,
    );
    if u_cl.is_empty() {
        return res;
    }
    let mut lc = dump_linear_chain(unitigs, &u_cl, i64::from(contig.length));
    if lc.is_empty() {
        return res;
    }
    let (max_idx, idx) = ctg_chain_graph(
        ug_graph,
        unitigs,
        &mut lc,
        i64::from(contig.length),
        G_CHAIN_BW,
        N_GCHAIN_RATE,
        UG_SKIP_GRAPH_N,
        UG_ITER_N,
        UG_DIS_N * 100,
    );
    if max_idx >= 0 {
        let mut dst_done: Vec<u64> = Vec::new();
        let mut out: Vec<usize> = Vec::new();
        let mut pathv: Vec<MgPathv> = Vec::new();
        let mut gchains: Vec<MgLchain> = Vec::new();
        let selected = select_max_ctg_chain(
            ug_graph,
            &idx,
            &lc,
            &mut dst_done,
            &mut out,
            &mut pathv,
            &mut gchains,
        );
        if selected > 0 {
            push_ctg_res(&mut res, &gchains);
        }
    }
    res
}

/// `gen_contig_path` (`inter.cpp` 21914) single-threaded: build the read-to-
/// unitig index of `ref` and map every `ctg` unitig onto it.
pub fn gen_contig_path(
    ctg: &UnitigGraph,
    reference: &UnitigGraph,
    read_lengths: &[u64],
) -> ScafRes {
    let read_count = read_lengths.len();
    let index = gen_r_ug_idx(&reference.unitigs, read_count, true);
    let mut res = ScafRes::new(ctg.unitigs.len());
    for (i, contig) in ctg.unitigs.iter().enumerate() {
        res.a[i] = gen_contig_path_for_contig(
            contig,
            &reference.unitigs,
            &index,
            read_lengths,
            &reference.graph,
        );
    }
    res
}

/// Text dump used to compare `gen_contig_path` against the reference under
/// gdb (one line per chain header and per block).
pub fn dump_scaf_res_text(path: &std::path::Path, res: &ScafRes) -> std::io::Result<()> {
    let mut text = String::new();
    for (i, contig) in res.a.iter().enumerate() {
        text.push_str(&format!("CTG {}\t{}\n", i, contig.n));
        for (k, header) in contig.bb.iter().take(contig.n).enumerate() {
            text.push_str(&format!(
                "HDR {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                k,
                header.hid,
                header.qs,
                header.qe,
                header.ts,
                header.te,
                header.aidx,
                u8::from(header.rev),
                u8::from(header.el),
                header.pchain
            ));
            for (z, block) in contig.bb[header.ts as usize..header.te as usize]
                .iter()
                .enumerate()
            {
                text.push_str(&format!(
                    "BLK {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    k,
                    z,
                    block.hid,
                    block.qs,
                    block.qe,
                    block.ts,
                    block.te,
                    block.aidx,
                    u8::from(block.rev),
                    u8::from(block.el),
                    block.pchain
                ));
            }
        }
    }
    std::fs::write(path, text)
}

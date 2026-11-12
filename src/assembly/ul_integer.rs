//! The `gen_ul2ul` integer-overlap engine core (`gfa_ut.cpp` 4368-5520):
//! the per-UL-read unitig strings (`ul_str_idx_t`), the anchored chaining of
//! string positions (`integer_chain`/`integer_chain_dp`) and the exact-overlap
//! materialisation (`update_exact_ul_ovlps`). The POA-consensus half of
//! `integer_normalize_ovlp` (`poa_cns_chain`/`update_raw_integer_seq`),
//! `integer_align_extention` and `gap_chain_check` are not ported yet; see
//! docs/consistency.md for the inventory.

use crate::assembly::bubble_chain::BubbleType;
use crate::assembly::ul_store::{UcBlock, UlReadVec};
use crate::container::ksort::{radix_sort_by_key_32, radix_sort_i64};

/// `ug_occ_w` (`gfa_ut.cpp` 4099): the number of the unitig's reads inside
/// `[is, ie)`; the full span returns the whole read count.
pub fn ug_occ_w(is: u64, ie: u64, read_lengths: &[u64], unitig: &[u64]) -> u64 {
    let unitig_len = unitig.iter().map(|b| (b & 0xffff_ffff) as u64).sum::<u64>();
    let unitig_len = unitig_len.max(
        unitig
            .iter()
            .map(|b| u64::from((b & 0xffff_ffff) as u32))
            .sum::<u64>(),
    );
    let _ = unitig_len;
    if is == 0 && ie == unitig_last_len(unitig, read_lengths) {
        return unitig.len() as u64;
    }
    let mut l = 0u64;
    let mut occ = 0u64;
    for &packed in unitig {
        let us = l;
        let ue = l + read_lengths[(packed >> 33) as usize];
        if is <= us && ie >= ue {
            occ += 1;
        }
        if us >= ie {
            break;
        }
        l += (packed & 0xffff_ffff) as u64;
    }
    occ
}

/// The unitig length as `u->len` (`ma_utg_t::len`).
fn unitig_last_len(unitig: &[u64], read_lengths: &[u64]) -> u64 {
    let _ = read_lengths;
    unitig.iter().map(|b| u64::from((b & 0xffff_ffff) as u32)).sum()
}

/// `ul_str_t` (`gfa_ut.cpp` 225): one UL read's unitig string.
#[derive(Clone, Debug, Default)]
pub struct UlStr {
    /// Packed `(block_index << 32) | vertex`; `cn` live entries.
    pub a: Vec<u64>,
    /// `cn` (and `n` — the reference keeps both equal).
    pub cn: u32,
    /// `is_cir`.
    pub is_cir: bool,
}

/// `ul_str_idx_t` (`gfa_ut.cpp` 236).
#[derive(Clone, Debug, Default)]
pub struct UlStrIdx {
    /// CSR: per unitig, the `(read, string position)` words.
    pub idx: Vec<u64>,
    pub occ: Vec<u64>,
    /// One string per UL read.
    pub str: Vec<UlStr>,
}

/// `integer_aln_t` (`gfa_ut.cpp` 251).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct IntegerAln {
    /// The target string position (`tk` in the reference).
    pub tk: u32,
    /// The query vertex `(hid<<1)|rev`.
    pub vq: u32,
    pub sc: u32,
    /// `(tid << 34) | (rev << 32) | qk`.
    pub tn_rev_qk: u64,
}

/// `ul_chain_t` (`gfa_ut.cpp` 265).
#[derive(Clone, Copy, Debug, Default)]
pub struct UlChain {
    /// The target `(hid<<1)|rev`, `u32::MAX` when unset.
    pub v: u32,
    /// `[s, e)` into the anchor array.
    pub s: u32,
    pub e: u32,
    pub sc: u64,
    pub q_sidx: u32,
    pub q_eidx: u32,
    pub t_sidx: u32,
    pub t_eidx: u32,
}

/// `ul2ul_t` (`gfa_ut.cpp` 106).
#[derive(Clone, Copy, Debug, Default)]
pub struct Ul2ulT {
    pub hid: u32,
    pub qs: u32,
    pub qe: u32,
    pub ts: u32,
    pub te: u32,
    pub qs_k: u32,
    pub qe_k: u32,
    pub ts_k: u32,
    pub te_k: u32,
    pub is_rev: bool,
    pub is_del: bool,
    pub is_ct: bool,
}

/// `ul2ul_item_t` (`gfa_ut.cpp` 119).
#[derive(Clone, Debug, Default)]
pub struct Ul2ulItem {
    pub a: Vec<Ul2ulT>,
    /// `id:31 | is_del:1`.
    pub id: u32,
    pub is_del: bool,
    /// `cn`: the live prefix of `a`.
    pub cn: usize,
}

/// `normlize_gdis` (`gfa_ut.cpp` 4368).
pub fn normlize_gdis(
    unitig_lengths: &[u32],
    unitig_reads: &[Vec<u64>],
    read_lengths: &[u64],
    i: &UcBlock,
    k: &UcBlock,
    is_i2k_forward: bool,
) -> i64 {
    let mut i_len = if !i.rev {
        i.qe as i64 + (unitig_lengths[i.hid as usize] as i64 - i.te as i64)
    } else {
        i.qe as i64 + i.ts as i64
    };
    let mut k_len = if !k.rev {
        k.qe as i64 + (unitig_lengths[k.hid as usize] as i64 - k.te as i64)
    } else {
        k.qe as i64 + k.ts as i64
    };
    if is_i2k_forward {
        i_len -= unitig_lengths[i.hid as usize] as i64;
        k_len -= unitig_lengths[k.hid as usize] as i64;
    }
    let _ = (unitig_reads, read_lengths);
    if i_len >= k_len {
        return i_len - k_len;
    }
    0
}

/// `normlize_gdis_exact` (`gfa_ut.cpp` 4424): the exact path distance along
/// the block chain between block indices `i > k`.
pub fn normlize_gdis_exact(
    unitig_lengths: &[u32],
    blocks: &[UcBlock],
    mut i: usize,
    k: usize,
    is_i2k_forward: bool,
) -> i64 {
    assert!(i > k);
    let bi = i;
    let mut li = i;
    let mut l = 0i64;
    while i != usize::MAX && i >= k {
        li = i;
        if blocks[i].pidx != u32::MAX && i > k {
            l += blocks[i].pdis as i64;
        }
        i = if blocks[i].pidx == u32::MAX {
            usize::MAX
        } else {
            blocks[i].pidx as usize
        };
    }
    if is_i2k_forward {
        l += unitig_lengths[blocks[li].hid as usize] as i64;
        l -= unitig_lengths[blocks[bi].hid as usize] as i64;
    }
    i = li;
    if i == k {
        return l;
    }
    assert!(i > k);
    let pk = k;
    let mut lk = k;
    let mut kk = k;
    while kk != usize::MAX && kk <= i {
        lk = kk;
        kk = if blocks[kk].aidx == u32::MAX {
            usize::MAX
        } else {
            blocks[kk].aidx as usize
        };
    }
    let mut k2 = lk;
    while k2 != usize::MAX && k2 != pk {
        l += blocks[k2].pdis as i64;
        k2 = if blocks[k2].pidx == u32::MAX {
            usize::MAX
        } else {
            blocks[k2].pidx as usize
        };
    }
    if is_i2k_forward {
        l += unitig_lengths[blocks[pk].hid as usize] as i64;
        l -= unitig_lengths[blocks[lk].hid as usize] as i64;
    }
    if l < 0 {
        l = 0;
    }
    let kk = lk;
    assert!(i > kk);
    l + normlize_gdis(
        unitig_lengths,
        &[],
        &[],
        &blocks[i],
        &blocks[kk],
        is_i2k_forward,
    )
}

/// `dis_check_integer_aln_t` (`gfa_ut.cpp` 4453).
#[allow(clippy::too_many_arguments)]
pub fn dis_check_integer_aln_t(
    store: &[UlReadVec],
    str_idx: &UlStrIdx,
    unitig_lengths: &[u32],
    li: &IntegerAln,
    lk: &IntegerAln,
    qid: usize,
    tid: usize,
    is_rev: bool,
    diff_rate: f64,
    hard_thres: i64,
) -> bool {
    debug_assert!((li.tn_rev_qk as u32) >= (lk.tn_rev_qk as u32));
    if (li.tn_rev_qk as u32) == (lk.tn_rev_qk as u32) {
        return false;
    }
    if li.tk <= lk.tk {
        return false;
    }
    let q = &str_idx.str[qid];
    let t = &str_idx.str[tid];
    let i_qk = li.tn_rev_qk as u32 as usize;
    let k_qk = lk.tn_rev_qk as u32 as usize;
    let (i_tk, k_tk) = if !is_rev {
        (li.tk as usize, lk.tk as usize)
    } else {
        (t.cn as usize - lk.tk as usize - 1, t.cn as usize - li.tk as usize - 1)
    };
    let iq = &store[qid].bb[q.a[i_qk] as usize >> 32];
    let it = &store[tid].bb[t.a[i_tk] as usize >> 32];
    let kq = &store[qid].bb[q.a[k_qk] as usize >> 32];
    let kt = &store[tid].bb[t.a[k_tk] as usize >> 32];
    let mut qlen = normlize_gdis(unitig_lengths, &[], &[], iq, kq, false);
    let mut tlen = normlize_gdis(unitig_lengths, &[], &[], it, kt, is_rev);
    let (mm, dd) = if qlen < tlen {
        (qlen, tlen - qlen)
    } else {
        (tlen, qlen - tlen)
    };
    if (dd as f64) <= (mm as f64 * diff_rate) || dd < hard_thres {
        return true;
    }
    qlen = normlize_gdis_exact(
        unitig_lengths,
        &store[qid].bb,
        q.a[i_qk] as usize >> 32,
        q.a[k_qk] as usize >> 32,
        false,
    );
    tlen = normlize_gdis_exact(
        unitig_lengths,
        &store[tid].bb,
        t.a[i_tk] as usize >> 32,
        t.a[k_tk] as usize >> 32,
        is_rev,
    );
    let (mm, dd) = if qlen < tlen {
        (qlen, tlen - qlen)
    } else {
        (tlen, qlen - tlen)
    };
    if (dd as f64) <= (mm as f64 * diff_rate) || dd < hard_thres {
        return true;
    }
    false
}

/// `integer_chain` (`gfa_ut.cpp` 4506): chains one target group of anchors;
/// on success the anchors are compacted into the chain order.
pub fn integer_chain(
    qid: u32,
    a: &mut [IntegerAln],
    offset: usize,
    str_idx: Option<&UlStrIdx>,
    store: Option<&[UlReadVec]>,
    unitig_lengths: &[u32],
    f_buf: &mut Vec<i64>,
    p_buf: &mut Vec<i64>,
    res: &mut UlChain,
) -> bool {
    res.v = u32::MAX;
    res.s = u32::MAX;
    res.e = u32::MAX;
    res.sc = u64::MAX;
    res.q_sidx = u32::MAX;
    res.q_eidx = u32::MAX;
    res.t_sidx = u32::MAX;
    res.t_eidx = u32::MAX;
    let a_n = a.len();
    if a_n == 0 {
        return false;
    }
    let tid = (a[0].tn_rev_qk >> 33) as u32;
    let is_rev = ((a[0].tn_rev_qk >> 32) & 1) != 0;
    let mut sc = 0i64;
    let mut i = 1usize;
    while i < a_n {
        sc += a[i].sc as i64;
        if a[i].tk <= a[i - 1].tk {
            break;
        }
        if (a[i].tn_rev_qk as u32) <= (a[i - 1].tn_rev_qk as u32) {
            break;
        }
        if let (Some(str_idx), Some(store)) = (str_idx, store) {
            if !dis_check_integer_aln_t(
                store,
                str_idx,
                unitig_lengths,
                &a[i],
                &a[i - 1],
                qid as usize,
                tid as usize,
                is_rev,
                0.08,
                2000,
            ) {
                break;
            }
        }
        i += 1;
    }
    if i >= a_n {
        sc += a[0].sc as i64;
        res.v = (a[0].tn_rev_qk >> 32) as u32 & 0x7fff_ffff;
        res.v = (a[0].tn_rev_qk >> 33) as u32;
        res.s = offset as u32;
        res.e = (offset + a_n) as u32;
        res.sc = sc as u64;
        return true;
    }
    f_buf.clear();
    p_buf.clear();
    f_buf.resize(a_n, 0);
    p_buf.resize(a_n, 0);
    let f: &mut Vec<i64> = f_buf;
    let p: &mut Vec<i64> = p_buf;
    let mut tf = -1i64;
    let mut ti = -1i64;
    for i in 0..a_n {
        let li = a[i];
        let csc = li.sc as i64;
        let mut max_f = csc;
        let mut max_k = -1i64;
        let mut k = i as i64 - 1;
        while k >= 0 {
            let lk = a[k as usize];
            if lk.tk >= li.tk || (lk.tn_rev_qk as u32) >= (li.tn_rev_qk as u32) {
                k -= 1;
                continue;
            }
            if let (Some(str_idx), Some(store)) = (str_idx, store) {
                if !dis_check_integer_aln_t(
                    store,
                    str_idx,
                    unitig_lengths,
                    &li,
                    &lk,
                    qid as usize,
                    tid as usize,
                    is_rev,
                    0.08,
                    2000,
                ) {
                    k -= 1;
                    continue;
                }
            }
            let sc = csc + f[k as usize];
            if sc > max_f {
                max_f = sc;
                max_k = k;
            }
            k -= 1;
        }
        f[i] = max_f;
        p[i] = max_k;
        if tf < max_f {
            tf = max_f;
            ti = i as i64;
        }
    }
    if ti < 0 {
        return false;
    }
    let mut order: Vec<usize> = Vec::new();
    let mut cur = ti as usize;
    while cur != usize::MAX {
        order.push(cur);
        cur = if p[cur] < 0 {
            usize::MAX
        } else {
            p[cur] as usize
        };
    }
    order.reverse();
    let mut sc = 0i64;
    let picked: Vec<IntegerAln> = order.iter().map(|&idx| a[idx]).collect();
    for (i, entry) in picked.iter().enumerate() {
        a[i] = *entry;
        sc += entry.sc as i64;
    }
    res.v = (a[0].tn_rev_qk >> 33) as u32;
    res.s = offset as u32;
    res.e = (offset + picked.len()) as u32;
    res.sc = sc as u64;
    true
}

/// `append_connective` (`gfa_ut.cpp` 4670).
fn append_connective(
    aln: &[IntegerAln],
    idx: &mut UlChain,
    str_i: usize,
    occ_thres: u64,
    res: &mut [u64],
) -> i64 {
    let kl = (idx.e - idx.s) as usize;
    let base = idx.s as usize;
    let mut fp = 0i64;
    let mut k = idx.sc as usize;
    let mut str_k = -1i64;
    while k < kl {
        str_k = (aln[base + k].tn_rev_qk) as u32 as i64;
        if str_k >= str_i as i64 {
            break;
        }
        k += 1;
    }
    idx.sc = k as u64;
    if k >= kl || str_k != str_i as i64 {
        return 0;
    }
    let mut k = k as i64 - 1;
    while k >= 0 {
        let str_k = (aln[base + k as usize].tn_rev_qk) as u32 as usize;
        if res[str_k] < occ_thres {
            res[str_k] += 1;
            if res[str_k] == occ_thres {
                fp += 1;
            }
        }
        k -= 1;
    }
    fp
}

/// `connective_conform` (`gfa_ut.cpp` 4691).
fn connective_conform(
    aln: &[IntegerAln],
    idx_a: &mut [UlChain],
    str_i0: usize,
    occ_thres: u64,
    is_cov_check: bool,
) -> bool {
    if str_i0 == 0 {
        return true;
    }
    let mut matched = 0i64;
    let mut exact = 0i64;
    for x in idx_a.iter_mut() {
        let kl = (x.e - x.s) as usize;
        let base = x.s as usize;
        let mut k = x.sc as usize;
        let mut str_k = usize::MAX;
        while k < kl {
            str_k = (aln[base + k].tn_rev_qk) as u32 as usize;
            if str_k >= str_i0 {
                break;
            }
            k += 1;
        }
        x.sc = k as u64;
        if k >= kl || str_k != str_i0 {
            continue;
        }
        if k > 0 {
            let str_k_prev = (aln[base + k - 1].tn_rev_qk) as u32 as usize;
            if str_k_prev + 1 == str_i0 {
                matched += 1;
                if aln[base + k - 1].tk + 1 == aln[base + k].tk {
                    exact += 1;
                }
            }
        }
    }
    if is_cov_check {
        if (matched as u64) < occ_thres {
            return false;
        }
    } else {
        assert!(matched as u64 >= occ_thres);
    }
    if exact == matched {
        return true;
    }
    if (exact as f64) > (matched as f64 * 0.51) && exact > matched / 2 {
        return true;
    }
    false
}

/// `integer_chain_dp` (`gfa_ut.cpp` 4727): the query-side DP over the unitig
/// string; returns the consensus positions in `o_buf`.
#[allow(clippy::too_many_arguments)]
pub fn integer_chain_dp(
    bub: &BubbleType,
    f_buf: &mut Vec<i64>,
    p_buf: &mut Vec<i64>,
    o_buf: &mut Vec<u64>,
    aln: &[IntegerAln],
    idx: &mut [UlChain],
    qstr: &UlStr,
    weights: &[u64],
    is_hom: bool,
    occ_thres: u64,
    mut corrected: Option<&mut u64>,
) -> usize {
    let q_n = qstr.cn as usize;
    if let Some(out) = corrected.as_deref_mut() {
        *out = 0;
    }
    for entry in idx.iter_mut() {
        entry.sc = 0;
    }
    f_buf.clear();
    p_buf.clear();
    o_buf.clear();
    f_buf.resize(q_n, 0);
    p_buf.resize(q_n, 0);
    o_buf.resize(q_n, 0);
    let mut n_skip = 0u64;
    let mut tf = -1i64;
    let mut tk = -1i64;
    for k in 0..q_n {
        let mut csc = weights[k] as i64;
        if !is_hom && crate::assembly::ul_renew_gfa::if_hom((qstr.a[k] >> 1) as u32, bub) {
            csc = -1;
            n_skip += 1;
        }
        let mut max_p = -1i64;
        let mut max_f = csc;
        if k > 0 && max_f >= 0 {
            let mut done_z = 0usize;
            if is_hom {
                for entry in o_buf[..k].iter_mut() {
                    *entry = 0;
                }
            } else {
                for z in 0..k {
                    o_buf[z] = 0;
                    if f_buf[z] == -1 {
                        o_buf[z] = occ_thres;
                        done_z += 1;
                    }
                }
            }
            for entry in idx.iter_mut() {
                if done_z >= k {
                    break;
                }
                done_z += append_connective(aln, entry, k, occ_thres, o_buf) as usize;
            }
            debug_assert!(done_z <= k);
            for z in (0..k).rev() {
                if o_buf[z] < occ_thres {
                    continue;
                }
                if f_buf[z] == -1 {
                    continue;
                }
                let sc = csc + f_buf[z];
                if sc > max_f {
                    max_f = sc;
                    max_p = z as i64;
                }
            }
        }
        f_buf[k] = max_f;
        p_buf[k] = max_p;
        if tf < max_f && max_f >= 0 {
            tf = max_f;
            tk = k as i64;
        }
    }
    if tk < 0 {
        return 0;
    }
    let mut done_z = 0usize;
    let mut k = tk as usize;
    while {
        done_z += 1;
        if p_buf[k] < 0 {
            false
        } else {
            k = p_buf[k] as usize;
            true
        }
    } {}
    if let Some(out) = corrected.as_deref_mut() {
        if n_skip as usize + done_z == q_n {
            for entry in idx.iter_mut() {
                entry.sc = 0;
            }
            *out = 0;
            let mut k = 1usize;
            while k < q_n {
                let cov_check = !(f_buf[k] >= 0 && f_buf[k - 1] >= 0);
                if !connective_conform(aln, idx, k, occ_thres, cov_check) {
                    break;
                }
                k += 1;
            }
            if k >= q_n {
                *out = 1;
            }
        }
    }
    let mut sc = done_z;
    let mut k = tk as usize;
    while {
        done_z -= 1;
        o_buf[done_z] = k as u64;
        if p_buf[k] < 0 {
            false
        } else {
            k = p_buf[k] as usize;
            true
        }
    } {}
    let _ = sc;
    sc
}

/// `uc_block_qse_cutoff` (`gfa_ut.cpp` 5220).
pub fn uc_block_qse_cutoff(c_ts: u32, c_te: u32, x: &UcBlock) -> (u32, u32) {
    assert!(c_ts >= x.ts && c_te <= x.te);
    let (l, r) = if !x.rev {
        (c_ts - x.ts, x.te - c_te)
    } else {
        (c_ts - x.ts, x.te - c_te)
    };
    let span = x.te - x.ts;
    let qspan = x.qe - x.qs;
    let nqs = x.qs + crate::assembly::ul::get_offset_adjust(l, span, qspan);
    let nqe = x.qe - crate::assembly::ul::get_offset_adjust(r, span, qspan);
    assert!(nqs >= x.qs && nqe <= x.qe);
    (nqs, nqe)
}

/// `update_exact_ul_ovlps` (`gfa_ut.cpp` 5235): the exact-overlap assembly
/// for one chain; the overlap coordinates land in `res`.
#[allow(clippy::too_many_arguments)]
pub fn update_exact_ul_ovlps(
    store: &[UlReadVec],
    str_idx: &[UlStr],
    aln: &[IntegerAln],
    idx: &mut UlChain,
    qid: usize,
    p_buf: &mut Vec<i64>,
    f_buf: &mut Vec<i64>,
    o_buf: &mut Vec<u64>,
    u_buf: &mut Vec<u64>,
    res: &mut Ul2ulT,
) -> bool {
    if idx.e <= idx.s {
        return false;
    }
    let q_str = &str_idx[qid];
    let tid = (aln[idx.s as usize].tn_rev_qk >> 33) as usize;
    let is_rev = (aln[idx.s as usize].tn_rev_qk >> 32) & 1 != 0;
    let t_str = &str_idx[tid];
    let x = &aln[idx.s as usize];
    let qk = x.tn_rev_qk as u32 as usize;
    let tk = if !is_rev {
        x.tk as usize
    } else {
        t_str.cn as usize - x.tk as usize - 1
    };
    if qk > 0 && x.tk > 0 {
        return false;
    }
    idx.q_sidx = x.tn_rev_qk as u32;
    idx.t_sidx = x.tk;
    let x = &aln[idx.e as usize - 1];
    let qk = x.tn_rev_qk as u32 as usize;
    let tk = if !is_rev {
        x.tk as usize
    } else {
        t_str.cn as usize - x.tk as usize - 1
    };
    if (qk + 1) < q_str.cn as usize && (x.tk as usize + 1) < t_str.cn as usize {
        return false;
    }
    idx.q_eidx = (x.tn_rev_qk as u32) + 1;
    idx.t_eidx = x.tk + 1;
    if (idx.q_eidx - idx.q_sidx) as usize != (idx.t_eidx - idx.t_sidx) as usize {
        return false;
    }
    if idx.q_eidx <= idx.q_sidx {
        return false;
    }
    {
        let mut qk = idx.q_sidx as usize;
        let q_end = idx.q_eidx as usize;
        let mut tk = idx.t_sidx as usize;
        let t_end = idx.t_eidx as usize;
        let mut broke = false;
        while qk < q_end && tk < t_end {
            let v = q_str.a[qk] as u32;
            let w = if !is_rev {
                t_str.a[tk] as u32
            } else {
                t_str.a[t_str.cn as usize - tk - 1] as u32 ^ 1
            };
            if v != w {
                broke = true;
                break;
            }
            qk += 1;
            tk += 1;
        }
        if broke || qk != q_end || tk != t_end {
            return false;
        }
    }
    p_buf.clear();
    f_buf.clear();
    o_buf.clear();
    u_buf.clear();
    let span = (idx.q_eidx - idx.q_sidx) as usize;
    p_buf.resize(span, 0);
    f_buf.resize(span, 0);
    o_buf.resize(span, 0);
    u_buf.resize(span, 0);
    let mut acc = 0usize;
    let mut chained = true;
    let mut tf = -1i64;
    let mut tz = -1i64;
    let mut ps_flag = false;
    let mut qk = idx.q_sidx as usize;
    let q_end = idx.q_eidx as usize;
    let mut tk = idx.t_sidx as usize;
    let t_end = idx.t_eidx as usize;
    while qk < q_end && tk < t_end {
        let qi = &store[qid].bb[q_str.a[qk] as usize >> 32];
        let ti = if !is_rev {
            &store[tid].bb[t_str.a[tk] as usize >> 32]
        } else {
            &store[tid].bb[t_str.a[t_str.cn as usize - tk - 1] as usize >> 32]
        };
        assert_eq!(qi.hid, ti.hid);
        let os = qi.ts.max(ti.ts);
        let oe = qi.te.min(ti.te);
        if oe <= os {
            qk += 1;
            tk += 1;
            continue;
        }
        let (q_ns, q_ne) = uc_block_qse_cutoff(os, oe, qi);
        let (t_ns, t_ne) = uc_block_qse_cutoff(os, oe, ti);
        o_buf[acc] = ((q_ns as u64) << 32) | q_ne as u64;
        u_buf[acc] = ((t_ns as u64) << 32) | t_ne as u64;
        let mut z = acc as i64 - 1;
        let mut csc = (q_ne - q_ns) as i64;
        if csc <= 0 {
            csc = 1;
        }
        let mut max_f = csc;
        let mut max_z = -1i64;
        let mut ps = false;
        if chained && z >= 0 {
            let qbs = (o_buf[z as usize] >> 32) as u32;
            let qbe = o_buf[z as usize] as u32;
            let tbs = (u_buf[z as usize] >> 32) as u32;
            let tbe = u_buf[z as usize] as u32;
            ps = if !is_rev {
                qbs <= q_ns && qbe <= q_ne && tbs <= t_ns && tbe <= t_ne
            } else {
                qbs <= q_ns && qbe <= q_ne && tbs >= t_ns && tbe >= t_ne
            };
            if ps {
                let sc = csc + f_buf[z as usize];
                if sc > max_f {
                    max_f = sc;
                    max_z = z;
                }
            } else {
                chained = false;
            }
        }
        if !ps || !chained {
            while z >= 0 {
                let qbs = (o_buf[z as usize] >> 32) as u32;
                let qbe = o_buf[z as usize] as u32;
                let tbs = (u_buf[z as usize] >> 32) as u32;
                let tbe = u_buf[z as usize] as u32;
                ps = if !is_rev {
                    qbs <= q_ns && qbe <= q_ne && tbs <= t_ns && tbe <= t_ne
                } else {
                    qbs <= q_ns && qbe <= q_ne && tbs >= t_ns && tbe >= t_ne
                };
                if !ps {
                    z -= 1;
                    continue;
                }
                let sc = csc + f_buf[z as usize];
                if sc > max_f {
                    max_f = sc;
                    max_z = z;
                }
                z -= 1;
            }
        }
        f_buf[acc] = max_f;
        p_buf[acc] = max_z;
        if tf < max_f {
            tf = max_f;
            tz = max_z;
        }
        ps_flag = ps;
        acc += 1;
        qk += 1;
        tk += 1;
    }
    let _ = ps_flag;
    if acc == 0 {
        return false;
    }
    assert!(tz >= 0);
    *res = Ul2ulT {
        hid: tid as u32,
        qs: u32::MAX,
        qe: u32::MAX,
        ts: u32::MAX,
        te: u32::MAX,
        is_rev,
        is_del: false,
        is_ct: false,
        ..Ul2ulT::default()
    };
    res.qs_k = idx.q_sidx;
    res.qe_k = idx.q_eidx;
    if !is_rev {
        res.ts_k = idx.t_sidx;
        res.te_k = idx.t_eidx;
    } else {
        res.ts_k = t_str.cn - idx.t_eidx;
        res.te_k = t_str.cn - idx.t_sidx;
    }
    res.qe = o_buf[tz as usize] as u32;
    if !is_rev {
        res.te = u_buf[tz as usize] as u32;
    } else {
        res.ts = (u_buf[tz as usize] >> 32) as u32;
    }
    let mut z = tz as usize;
    loop {
        res.qs = (o_buf[z] >> 32) as u32;
        if !is_rev {
            res.ts = (u_buf[z] >> 32) as u32;
        } else {
            res.te = u_buf[z] as u32;
        }
        if p_buf[z] < 0 {
            break;
        }
        z = p_buf[z] as usize;
    }
    assert!(res.qs <= res.qe && res.ts <= res.te);
    let mut qs = res.qs as i64;
    let mut qe = res.qe as i64;
    let mut rs = res.ts as i64;
    let mut re = res.te as i64;
    if is_rev {
        rs = store[tid].rlen as i64 - res.te as i64;
        re = store[tid].rlen as i64 - res.ts as i64;
    }
    if qs <= rs {
        rs -= qs;
        qs = 0;
    } else {
        qs -= rs;
        rs = 0;
    }
    let qtail = store[qid].rlen as i64 - qe;
    let rtail = store[tid].rlen as i64 - re;
    if qtail <= rtail {
        qe = store[qid].rlen as i64;
        re += qtail;
    } else {
        re = store[tid].rlen as i64;
        qe += rtail;
    }
    res.qs = qs as u32;
    res.qe = qe as u32;
    res.ts = rs as u32;
    res.te = re as u32;
    if is_rev {
        res.ts = store[tid].rlen as u32 - re as u32;
        res.te = store[tid].rlen as u32 - rs as u32;
    }
    true
}

/// The CSR string index of `init_ul_str_idx_t` (`gfa_ut.cpp` 4137).
pub fn init_ul_str_idx(store: &[UlReadVec], unitig_count: usize) -> UlStrIdx {
    let mut str_idx = UlStrIdx::default();
    str_idx.str = vec![UlStr::default(); store.len()];
    str_idx.idx = vec![0u64; unitig_count + 1];
    for (i, read) in store.iter().enumerate() {
        for (k, block) in read.bb.iter().enumerate() {
            if block.base || !block.el || block.pchain == 0 {
                continue;
            }
            let packed = ((k as u64) << 32) + ((block.hid << 1) + u32::from(block.rev)) as u64;
            str_idx.idx[block.hid as usize] += 1;
            str_idx.str[i].a.push(packed);
        }
        str_idx.str[i].cn = str_idx.str[i].a.len() as u32;
    }
    let mut l = 0u64;
    for entry in str_idx.idx.iter_mut() {
        let m = *entry;
        *entry = l;
        l += m;
    }
    str_idx.occ = vec![0u64; l as usize];
    for k in 0..unitig_count {
        let start = str_idx.idx[k] as usize;
        let end = str_idx.idx[k + 1] as usize;
        if end > start {
            str_idx.occ[end - 1] = 0;
        }
    }
    for (k, z) in str_idx.str.iter_mut().enumerate() {
        for i in 0..z.cn as usize {
            let vertex = (z.a[i] as u32) as usize >> 1;
            let start = str_idx.idx[vertex] as usize;
            let end = str_idx.idx[vertex + 1] as usize;
            if end > start {
                let m;
                if str_idx.occ[end - 1] == (end - start - 1) as u64 {
                    m = end - start - 1;
                    str_idx.occ[end - 1] = ((k as u64) << 32) | i as u64;
                } else {
                    m = str_idx.occ[end - 1] as usize;
                    let next = str_idx.occ[m] as usize;
                    str_idx.occ[next] = ((k as u64) << 32) | i as u64;
                    str_idx.occ[m] = next as u64 + 1;
                }
                let v0 = (z.a[(str_idx.occ[m]) as u32 as usize] as u32) as u32;
                let mut m = m as i64;
                while m > 0 && !z.is_cir {
                    m -= 1;
                    if (str_idx.occ[m as usize] >> 32) != k as u64 {
                        break;
                    }
                    let v1 = (z.a[(str_idx.occ[m as usize]) as u32 as usize] as u32) as u32;
                    if v0 == v1 {
                        z.is_cir = true;
                    }
                }
            }
        }
    }
    str_idx
}

/// `utg_cover_read_occ_by_qs` (`gfa_ut.cpp` 4890).
pub fn utg_cover_read_occ_by_qs(
    oqs: i64,
    oqe: i64,
    x: &UcBlock,
    read_lengths: &[u64],
    unitig_reads: &[u64],
) -> i64 {
    assert!(oqs >= i64::from(x.qs) && oqs <= i64::from(x.qe));
    assert!(oqe >= i64::from(x.qs) && oqe <= i64::from(x.qe));
    assert!(oqe > oqs);
    let mut s_off = crate::assembly::ul::get_offset_adjust(
        (oqs - i64::from(x.qs)) as u32,
        x.qe - x.qs,
        x.te - x.ts,
    ) as i64;
    let mut e_off = crate::assembly::ul::get_offset_adjust(
        (i64::from(x.qe) - oqe) as u32,
        x.qe - x.qs,
        x.te - x.ts,
    ) as i64;
    if x.rev {
        std::mem::swap(&mut s_off, &mut e_off);
    }
    ug_occ_w(
        u64::from(x.ts + s_off as u32),
        u64::from(x.te - e_off as u32),
        read_lengths,
        unitig_reads,
    ) as i64
}

/// `estimate_ul_len` (`gfa_ut.cpp` 4905): the string position at the end of
/// a query-space extension, writing per-node coverages into `cov_buf`.
pub fn estimate_ul_len(
    raw: &UlReadVec,
    str_: &UlStr,
    idx: usize,
    ext_len: i64,
    cov_buf: &mut [u64],
    read_lengths: &[u64],
    unitig_reads: impl Fn(u32) -> Vec<u64>,
) -> usize {
    if ext_len == 0 {
        return idx;
    }
    let str_n = str_.cn as i64;
    let mut k = idx as i64;
    if ext_len < 0 {
        let mut aim_s = i64::from(raw.bb[str_.a[idx] as usize >> 32].qs) + ext_len;
        if aim_s < 0 {
            aim_s = 0;
        }
        let aim_e = i64::from(raw.bb[str_.a[idx] as usize >> 32].qe);
        k = idx as i64 - 1;
        while k >= 0 {
            let block = &raw.bb[str_.a[k as usize] as usize >> 32];
            let qs = i64::from(block.qs);
            let qe = i64::from(block.qe);
            let ovlp_s = aim_s.max(qs);
            let ovlp_e = aim_e.min(qe);
            if ovlp_e <= ovlp_s {
                k += 1;
                break;
            }
            let cov_occ = utg_cover_read_occ_by_qs(
                ovlp_s,
                ovlp_e,
                block,
                read_lengths,
                &unitig_reads(block.hid),
            );
            if cov_occ == 0 {
                k += 1;
                break;
            }
            cov_buf[k as usize] = cov_occ as u64;
            k -= 1;
        }
        if k < 0 {
            k += 1;
        }
    } else {
        let aim_s = i64::from(raw.bb[str_.a[idx] as usize >> 32].qs);
        let mut aim_e = i64::from(raw.bb[str_.a[idx] as usize >> 32].qe) + ext_len;
        if aim_e > i64::from(raw.rlen) {
            aim_e = i64::from(raw.rlen);
        }
        k = idx as i64 + 1;
        while k < str_n {
            let block = &raw.bb[str_.a[k as usize] as usize >> 32];
            let qs = i64::from(block.qs);
            let qe = i64::from(block.qe);
            let ovlp_s = aim_s.max(qs);
            let ovlp_e = aim_e.min(qe);
            if ovlp_e <= ovlp_s {
                k -= 1;
                break;
            }
            let cov_occ = utg_cover_read_occ_by_qs(
                ovlp_s,
                ovlp_e,
                block,
                read_lengths,
                &unitig_reads(block.hid),
            );
            if cov_occ == 0 {
                k -= 1;
                break;
            }
            cov_buf[k as usize] = cov_occ as u64;
            k += 1;
        }
        if k >= str_n {
            k -= 1;
        }
    }
    k as usize
}

/// `integer_align_extention` (`gfa_ut.cpp` 4977).
#[allow(clippy::too_many_arguments)]
pub fn integer_align_extention(
    store: &[UlReadVec],
    qid: usize,
    q_str: &UlStr,
    tid: usize,
    t_str: &UlStr,
    is_rev: bool,
    q_cov_buf: &mut [u64],
    t_cov_buf: &mut [u64],
    aln_pair: &IntegerAln,
    is_backward: bool,
    read_lengths: &[u64],
    unitig_reads: impl Fn(u32) -> Vec<u64>,
) -> (usize, usize) {
    let qk = aln_pair.tn_rev_qk as u32 as usize;
    let tk = if !is_rev {
        aln_pair.tk as usize
    } else {
        t_str.cn as usize - aln_pair.tk as usize - 1
    };
    let q_b = &store[qid].bb[q_str.a[qk] as usize >> 32];
    let t_b = &store[tid].bb[t_str.a[tk] as usize >> 32];
    let qlen = i64::from(store[qid].rlen);
    let tlen = i64::from(store[tid].rlen);
    let (qoff, toff) = if is_backward {
        (
            i64::from(q_b.qs),
            if is_rev {
                tlen - i64::from(t_b.qe)
            } else {
                i64::from(t_b.qs)
            },
        )
    } else {
        (
            qlen - i64::from(q_b.qe),
            if is_rev {
                i64::from(t_b.qs)
            } else {
                tlen - i64::from(t_b.qe)
            },
        )
    };
    let ext = qoff.min(toff);
    let (q_end, t_end) = if is_backward {
        (
            estimate_ul_len(&store[qid], q_str, qk, -ext, q_cov_buf, read_lengths, &unitig_reads),
            estimate_ul_len(
                &store[tid],
                t_str,
                tk,
                if is_rev { ext } else { -ext },
                t_cov_buf,
                read_lengths,
                &unitig_reads,
            ),
        )
    } else {
        (
            estimate_ul_len(&store[qid], q_str, qk, ext, q_cov_buf, read_lengths, &unitig_reads),
            estimate_ul_len(
                &store[tid],
                t_str,
                tk,
                if is_rev { -ext } else { ext },
                t_cov_buf,
                read_lengths,
                &unitig_reads,
            ),
        )
    };
    (q_end, t_end)
}

/// `gap_chain_check` (`gfa_ut.cpp` 5005).
#[allow(clippy::too_many_arguments)]
pub fn gap_chain_check(
    unitig_lengths: &[u32],
    i0: &IntegerAln,
    i1: &IntegerAln,
    is_rev: bool,
    q_str: &UlStr,
    t_str: &UlStr,
    q_block: &[UcBlock],
    t_block: &[UcBlock],
    diff_rate: f64,
    hard_thres: i64,
) -> bool {
    let i0_qk_str = i0.tn_rev_qk as u32 as usize;
    let i1_qk_str = i1.tn_rev_qk as u32 as usize;
    let (i0_tk_str, i1_tk_str) = if !is_rev {
        (i0.tk as usize, i1.tk as usize)
    } else {
        (
            t_str.cn as usize - i0.tk as usize - 1,
            t_str.cn as usize - i1.tk as usize - 1,
        )
    };
    debug_assert!(i0_qk_str < i1_qk_str);
    debug_assert!(i0_tk_str < i1_tk_str);
    let mut q_near = 0u32;
    let mut t_near = 0u32;
    if i0_qk_str + 1 == i1_qk_str {
        q_near = 1;
    }
    if i0_tk_str + 1 == i1_tk_str {
        t_near = 1;
    }
    if (q_near + t_near) != 1 {
        return false;
    }
    let i0_qk = q_str.a[i0_qk_str] as usize >> 32;
    let i1_qk = q_str.a[i1_qk_str] as usize >> 32;
    let i0_tk = t_str.a[i0_tk_str] as usize >> 32;
    let i1_tk = t_str.a[i1_tk_str] as usize >> 32;
    let walk = |blocks: &[UcBlock], mut k: usize, e: usize| -> Option<i64> {
        let mut l = 0i64;
        while k != usize::MAX && k != e {
            if blocks[k].pidx != u32::MAX {
                l += blocks[k].pdis as i64;
            }
            k = if blocks[k].pidx == u32::MAX {
                usize::MAX
            } else {
                blocks[k].pidx as usize
            };
        }
        if k != e {
            None
        } else {
            Some(l)
        }
    };
    let Some(mut lq) = walk(q_block, i1_qk, i0_qk) else {
        return false;
    };
    let Some(mut lt) = walk(t_block, i1_tk, i0_tk) else {
        return false;
    };
    if is_rev {
        lt += unitig_lengths[t_block[i0_tk].hid as usize] as i64;
        lt -= unitig_lengths[t_block[i1_tk].hid as usize] as i64;
    }
    if lt < 0 {
        lt = 0;
    }
    let check = |min: i64, max: i64, diff_rate: f64, hard_thres: i64| -> bool {
        let dif = max - min;
        !((dif as f64) >= (min as f64 * diff_rate))
    };
    let (min, max) = if lq <= lt { (lq, lt) } else { (lt, lq) };
    if !check(min, max, diff_rate, hard_thres) {
        return true;
    }
    lq = normlize_gdis(unitig_lengths, &[], &[], &q_block[i1_qk], &q_block[i0_qk], false);
    lt = normlize_gdis(unitig_lengths, &[], &[], &t_block[i1_tk], &t_block[i0_tk], is_rev);
    let (min, max) = if lq <= lt { (lq, lt) } else { (lt, lq) };
    let dif = max - min;
    if (dif as f64) >= (min as f64 * diff_rate) || dif < hard_thres {
        return false;
    }
    true
}

//! The UL gchain cluster (`inter.cpp` 5785–16660, 18852, 21892): re-derives
//! every UL read's block chain from the current read graph — `gl_rg2ug_gen`
//! projects the blocks onto unitigs, `gen_linear_chains` + `gl_chain_graph`
//! chain them across the unitig graph, `select_max_gchain` picks the best
//! chains and `update_ul_vec_t` rewrites the blocks with the new paths
//! (`ul_refine_alignment`).
//!
//! The DP buffers of the reference (`Chain_Data`, `st_mt_t`, scratch vectors)
//! are locals here; the shared helpers (`gl_chain_linear`, `hc_shortest_k`,
//! `gen_gchain_track`, `dump_linear_chain`, `max_ovlp`, `get_ecov_adv`) come
//! from [`crate::assembly::ul`].

use crate::assembly::asg::{Asg, MaUtg};
use crate::assembly::hit::{ma_hit_to_arc, ReadToUnitig};
use crate::assembly::ul::{
    dump_linear_chain, find_mg_lchain_max, gen_gchain_track, get_ecov_adv, gl_chain_linear,
    hc_shortest_k, infer_rovlp, max_ovlp, set_ul_ov_t_by_mg_lchain, MgLchain, MgPathDst, MgPathv,
    UlOv, UtgRidDt,
};
use crate::assembly::ul_index::{UlIndex, UtgCt};
use crate::assembly::ul::UcBlock;
use crate::assembly::ul_store::{gen_ul_vec_rid_t, AllUl, UlReadVec};
use crate::assembly::ul_src::{update_ovlp_src, update_ovlp_src_bl, UlSourceOptions};
use crate::command_line::options::AssemblyOptions;
use crate::logging;
use crate::overlap::paf::PafHit;
use crate::reads::ReadStore;

/// `UG_SKIP_N` (`inter.h`).
const UG_SKIP_N: i64 = 100;
/// `UG_ITER_N` (`inter.h`).
const UG_ITER_N: i64 = 5000;
/// `UG_DIS_N` (`inter.h`).
const UG_DIS_N: i64 = 50000;
/// `UG_SKIP_GRAPH_N` (`inter.h`).
const UG_SKIP_GRAPH_N: i64 = 72;
/// `G_CHAIN_BW` (`inter.h`).
const G_CHAIN_BW: i64 = 16;
/// `G_CHAIN_INDEL` (`inter.h`).
const G_CHAIN_INDEL: i64 = 100;
/// `P_CHAIN_COV` (`inter.h`).
const P_CHAIN_COV: f64 = 0.985;
/// `N_GCHAIN_RATE` (`inter.h`).
const N_GCHAIN_RATE: f64 = 0.04;
/// `W_CHN_PEN_GAP` (`inter.h`).
const W_CHN_PEN_GAP: f64 = 0.1;
/// `GC_OFFSET_POS` (`inter.h`).
const GC_OFFSET_POS: i64 = 8;
/// `GC_OFFSET_RATE` (`inter.cpp` 54).
const GC_OFFSET_RATE: f64 = 0.0001;

/// `mg_log2` (`inter.cpp` 654): the bit-level approximate log2.
fn mg_log2(x: f32) -> f32 {
    let mut z = f32::to_bits(x);
    let mut log_2 = (((z >> 23) & 255) - 128) as f32;
    z &= !(255u32 << 23);
    z += 127u32 << 23;
    let mantissa = f32::from_bits(z);
    log_2 + (-0.34484843 * mantissa + 2.02466578) * mantissa - 0.67487759
}

/// `retrieve_u_cov_region` (`Process_Read.cpp` 1827): the summed depth over
/// `[s, e)` of unitig `id`'s coverage track. `cursor` is the caller's cached
/// region index (`pi`).
pub fn retrieve_u_cov_region(
    index: &UlIndex,
    id: usize,
    strand: bool,
    s: u64,
    e: u64,
    cursor: &mut i64,
) -> u64 {
    let track = &index.coverage;
    let start = track.idx[id] as usize;
    let end = track.idx[id + 1] as usize;
    let a_n = end - start;
    if a_n == 0 {
        return 0;
    }
    let mut s = s;
    let mut e = e;
    let mut cc_i = if *cursor + 1 >= a_n as i64 || *cursor < 0 {
        0
    } else {
        *cursor
    };
    if strand {
        let total = u64::from(index.unitig.unitigs[id].length);
        let old_s = s;
        s = total - e;
        e = total - old_s;
    }
    let mut k = cc_i as usize;
    let tk = s;
    if tk < track.interval[start + k] >> 32 {
        while k > 0 {
            if tk >= track.interval[start + k] >> 32 && tk < track.interval[start + k + 1] >> 32 {
                break;
            }
            k -= 1;
        }
    } else if tk >= track.interval[start + k + 1] >> 32 {
        while k + 1 < a_n {
            if tk >= track.interval[start + k] >> 32 && tk < track.interval[start + k + 1] >> 32 {
                break;
            }
            k += 1;
        }
    }
    cc_i = k as i64;
    *cursor = cc_i;

    let mut total = 0u64;
    while k + 1 < a_n {
        let ts = track.interval[start + k] >> 32;
        let te = track.interval[start + k + 1] >> 32;
        let cc = track.interval[start + k] & 0xFFFF_FFFF;
        let overlap = if e.min(te) > s.max(ts) {
            e.min(te) - s.max(ts)
        } else {
            0
        };
        total += overlap * cc;
        if e >= track.interval[start + k] >> 32 && e < track.interval[start + k + 1] >> 32 {
            break;
        }
        k += 1;
    }
    total
}

/// `rov2uov` (`inter.cpp` 13274): projects one block of a HiFi read onto the
/// unitig occurrence `map`.
#[allow(clippy::too_many_arguments)]
fn rov2uov(
    read_id: usize,
    index: &UlIndex,
    read_lengths: &[u64],
    map: &UtgRidDt,
    block: &UcBlock,
    adjust_rev: bool,
) -> Option<UlOv> {
    let ori = map.u & 1;
    let (ts, te) = if ori == 0 {
        (u64::from(block.ts), u64::from(block.te))
    } else {
        let len = read_lengths[read_id];
        (len - u64::from(block.te), len - u64::from(block.ts))
    };
    let ts = ts + u64::from(map.pos);
    let te = te + u64::from(map.pos);
    let unitig_length = u64::from(index.unitig.unitigs[(map.u >> 1) as usize].length);
    // `if(ts >= 0 && te <= len)`: the unsigned `ts` never underflows the
    // check, so only the upper bound is tested
    if te <= unitig_length {
        let mut ts = ts as u32;
        let mut te = te as u32;
        let rev = (block.rev as u32) != ori;
        if adjust_rev && rev {
            // for linear chaining; both expressions read the pre-flip values
            let (orig_ts, orig_te) = (ts, te);
            ts = unitig_length as u32 - orig_te;
            te = unitig_length as u32 - orig_ts;
        }
        Some(UlOv {
            qn: 0,
            qs: block.qs,
            qe: block.qe,
            tn: map.u >> 1,
            ts,
            te,
            sec: 0,
            el: block.el,
            rev,
        })
    } else {
        None
    }
}

/// `gl_rg2ug_gen` (`inter.cpp` 13312): the unitig projections of the kept
/// blocks of `r_cl`.
fn gl_rg2ug_gen(
    r_cl: &UlReadVec,
    u_cl: &mut Vec<UlOv>,
    index: &UlIndex,
    read_lengths: &[u64],
    is_el: bool,
    n_pchain: u8,
) {
    u_cl.clear();
    for (k, block) in r_cl.bb.iter().enumerate() {
        if block.base {
            continue;
        }
        if is_el && !block.el {
            continue;
        }
        if block.pchain == n_pchain {
            continue;
        }
        let read_id = block.hid as usize;
        let start = index.read_unitig.idx[read_id] as usize;
        let end = index.read_unitig.idx[read_id + 1] as usize;
        for map in &index.read_unitig.p[start..end] {
            let Some(mut p) = rov2uov(read_id, index, read_lengths, map, block, true) else {
                continue;
            };
            p.el = true;
            p.tn <<= 1;
            p.tn |= u32::from(p.rev);
            p.qn = k as u32; // for linear chain
            u_cl.push(p);
        }
    }
}

/// `comput_linear_sc` (`inter.cpp` 13489).
fn comput_linear_sc(li: &UlOv, lj: &UlOv, diff_ec_ul: f64, bw: i64) -> i64 {
    // li is the suffix of lj
    if lj.te > li.te {
        return i32::MIN as i64;
    }
    let dq = i64::from(li.qe - lj.qs);
    let dt = i64::from(li.te - lj.ts);
    let dd = (dq - dt).abs();
    let mut mm = dq.max(dt);
    mm = (mm as f64 * diff_ec_ul) as i64;
    if mm < bw {
        mm = bw;
    }
    if dd > mm {
        return i32::MIN as i64;
    }
    i64::from(li.qe - li.qs)
}

/// `linear_chain_dp_adv` (`inter.cpp` 13506): the forward linear-chain DP
/// over one unitig's projections; the surviving chains keep `off`/`tn` links.
#[allow(clippy::too_many_arguments)]
fn linear_chain_dp_adv(
    ch: &mut [UlOv],
    sv: &mut [UlOv],
    bw: i64,
    diff_ec_ul: f64,
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
    unitig_lengths: &[u32],
    chain_offset: i64,
) -> usize {
    let ch_n = ch.len();
    if ch_n == 0 {
        return 0;
    }
    let mut f = vec![0i64; ch_n];
    let mut p = vec![0i64; ch_n];
    let mut t: Vec<u64> = vec![0; ch_n];

    crate::container::ksort::radix_sort_by_key_32(ch, &|entry: &UlOv| entry.qe);
    let mut j = 0usize;
    for i in 1..=ch_n {
        if i == ch_n || ch[i].qe != ch[j].qe {
            if i - j > 1 {
                crate::container::ksort::radix_sort_by_key_32(&mut ch[j..i], &|e: &UlOv| e.qs);
            }
            j = i;
        }
    }

    let mut st;
    let mut max_ii: i64 = -1;
    for i in 0..ch_n {
        let li = ch[i];
        let csc = i64::from(li.qe - li.qs);
        let mut mm_sc = csc;
        let mut mm_idx: i64 = -1;
        let mut n_skip = 0i64;
        st = if (i as i64) < max_iter { 0 } else { i as i64 - max_iter };
        let mut end_j: i64;
        let mut jj = i as i64 - 1;
        while jj >= st {
            let lj = ch[jj as usize];
            let sc = comput_linear_sc(&li, &lj, diff_ec_ul, bw);
            if sc == i32::MIN as i64 {
                jj -= 1;
                continue;
            }
            let sc = sc + f[jj as usize];
            if sc > mm_sc {
                mm_sc = sc;
                mm_idx = jj;
                if n_skip > 0 {
                    n_skip -= 1;
                }
            } else if t[jj as usize] == i as u64 {
                n_skip += 1;
                if n_skip > max_skip {
                    break;
                }
            }
            if p[jj as usize] >= 0 {
                t[p[jj as usize] as usize] = i as u64;
            }
            jj -= 1;
        }
        end_j = jj;
        if max_ii < 0 || (i64::from(ch[i].qe) > i64::from(ch[max_ii as usize].qe) + max_dis) {
            let mut max = i32::MIN as i64;
            max_ii = -1;
            let mut j2 = i as i64 - 1;
            while j2 >= st && i64::from(ch[i].qe) <= max_dis + i64::from(ch[j2 as usize].qe) {
                if max < f[j2 as usize] {
                    max = f[j2 as usize];
                    max_ii = j2;
                }
                j2 -= 1;
            }
        }

        if max_ii >= 0 && max_ii < end_j {
            let lj = ch[max_ii as usize];
            let sc = comput_linear_sc(&li, &lj, diff_ec_ul, bw);
            if sc != i32::MIN as i64 {
                let sc = sc + f[max_ii as usize];
                if sc > mm_sc {
                    mm_sc = sc;
                    mm_idx = max_ii;
                }
            }
        }

        f[i] = mm_sc;
        p[i] = mm_idx;
        if max_ii < 0
            || (i64::from(ch[i].qe) <= max_dis + i64::from(ch[max_ii as usize].qe)
                && f[max_ii as usize] < f[i])
        {
            max_ii = i as i64;
        }
        ch[i].sec = if mm_idx < 0 {
            0x3FFF_FFFF
        } else {
            (i as i64 - mm_idx) as u32
        };
        sv[i] = ch[i];
    }

    let mut n_u = 0usize;
    for k in (0..ch_n).rev() {
        if t[k] != 0 {
            continue;
        }
        let mut i = k;
        ch[n_u] = sv[i];
        let sc = f[i];
        loop {
            if sv[i].qs < ch[n_u].qs {
                ch[n_u].qs = sv[i].qs;
            }
            if sv[i].ts < ch[n_u].ts {
                ch[n_u].ts = sv[i].ts;
            }
            if sv[i].qe > ch[n_u].qe {
                ch[n_u].qe = sv[i].qe;
            }
            if sv[i].te > ch[n_u].te {
                ch[n_u].te = sv[i].te;
            }
            t[i] = 1;
            if p[i] < 0 {
                break;
            }
            i = p[i] as usize;
        }
        let len = i64::from(unitig_lengths[ch[n_u].tn as usize]);
        let (its, ite) = adjust_rev_tse(&ch[n_u], len);
        ch[n_u].ts = its as u32;
        ch[n_u].te = ite as u32;
        ch[n_u].sec = if sc > 0x3FFF_FFFF { 0x3FFF_FFFF } else { sc as u32 };
        ch[n_u].qn = (k as i64 + chain_offset) as u32;
        n_u += 1;
    }
    for i in 0..ch_n {
        let len = i64::from(unitig_lengths[sv[i].tn as usize]);
        let (its, ite) = adjust_rev_tse(&sv[i], len);
        sv[i].ts = its as u32;
        sv[i].te = ite as u32;
        let k = p[i];
        sv[i].tn = if k >= 0 {
            (k + chain_offset) as u32
        } else {
            u32::MAX
        };
    }
    n_u
}

/// `adjust_rev_tse` (`inter.cpp` 13379).
fn adjust_rev_tse(x: &UlOv, tlen: i64) -> (i64, i64) {
    let mut ts = i64::from(x.ts);
    let mut te = i64::from(x.te);
    if x.rev {
        ts = tlen - i64::from(x.te);
        te = tlen - i64::from(x.ts);
    }
    (ts, te)
}

/// `gen_linear_chains` (`inter.cpp` 13571). `buf` receives the DP's per-entry
/// scratch copies (`init`/`ll->tk` in the reference): one entry per original
/// position, with `tn` = the predecessor's chain index and `qn` kept as the
/// block index — exactly what `update_existing_anchors` walks.
fn gen_linear_chains(
    res: &mut Vec<UlOv>,
    buf: &mut Vec<UlOv>,
    index: &UlIndex,
    bw: i64,
    diff_ec_ul: f64,
    qlen: i64,
) {
    let unitig_lengths: Vec<u32> = index
        .unitig
        .unitigs
        .iter()
        .map(|unitig| unitig.length)
        .collect();
    crate::container::ksort::radix_sort_by_key_32(res, &|entry: &UlOv| entry.tn);
    let n = res.len();
    buf.clear();
    buf.resize(n, UlOv::default());
    let mut l = 0usize;
    let mut m = 0usize;
    for k in 1..=n {
        if k == n || res[k].tn != res[l].tn {
            for z in l..k {
                res[z].tn >>= 1;
            }
            let mut sv: Vec<UlOv> = vec![UlOv::default(); k - l];
            let an = l
                + linear_chain_dp_adv(
                    &mut res[l..k],
                    &mut sv,
                    bw,
                    diff_ec_ul,
                    UG_SKIP_N,
                    UG_ITER_N,
                    UG_DIS_N,
                    &unitig_lengths,
                    l as i64,
                );
            buf[l..k].copy_from_slice(&sv);
            for z in l..an {
                res[m] = res[z];
                m += 1;
            }
            l = k;
        }
    }
    res.truncate(m);
}

/// `hc_target_len` (`inter.cpp` 5591).
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

/// `cal_gchain_sc` (`inter.cpp` 5607).
fn cal_gchain_sc(
    dj: &MgPathDst,
    li: &MgLchain,
    f: &[i64],
    b_w: i64,
    diff_thre: f64,
    chn_pen_gap: f64,
) -> i64 {
    if dj.n_path == 0 {
        return i32::MIN as i64;
    }
    let mut gap = i64::from(dj.dist - dj.target_dist);
    if gap < 0 {
        gap = -gap;
    }
    if (gap as f64 > dj.target_dist as f64 * diff_thre) && gap > b_w {
        return i32::MIN as i64;
    }
    // the reference computes the penalties in `f32` and truncates once
    let mut sc = i32::try_from(li.score).unwrap_or(i32::MAX);
    let lin_pen = (chn_pen_gap as f32) * (gap as f32);
    let log_pen = if gap >= 2 { mg_log2(gap as f32) } else { 0.0 };
    sc -= (lin_pen + log_pen) as i32;
    i64::from(sc) + f[dj.meta as usize]
}

/// `primary_chain_check` (`inter.cpp` 5632).
fn primary_chain_check(idx: &[u64], a: &[MgLchain]) -> bool {
    let idx_n = idx.len();
    if idx_n == 0 {
        return false;
    }
    let mut m_sc: i64 = -1;
    let mut m_qn = 0usize;
    let mut m = (0usize, 0usize, 0u32, 0u32); // (ts, te, qs, qe)
    let mut a_n = 0usize;
    for (i, &entry) in idx.iter().enumerate() {
        if (entry >> 32) as i64 > m_sc {
            m_sc = (entry >> 32) as i64;
            m_qn = i;
            let ts = a_n;
            let te = a_n + (entry as u32) as usize;
            m = (ts, te, a[ts].qs as u32, a[te - 1].qe as u32);
        }
        a_n += (entry as u32) as usize;
    }
    let (m_ts, m_te, m_qs, m_qe) = m;
    debug_assert!(a[m_ts].qs as u32 <= a[m_te - 1].qs as u32 && a[m_te - 1].qe as u32 >= m_qe);

    let mut a_n = 0usize;
    for (i, &entry) in idx.iter().enumerate() {
        let s_idx = a[a_n].qs as i64;
        a_n += (entry as u32) as usize;
        let e_idx = a[a_n - 1].qe as i64;
        if i == m_qn {
            continue;
        }
        let ovlp = if m_qe.min(e_idx as u32) > m_qs.max(s_idx as u32) {
            m_qe.min(e_idx as u32) - m_qs.max(s_idx as u32)
        } else {
            0
        };
        let novlp = (e_idx - s_idx) - i64::from(ovlp);
        if novlp as f64 > (m_qe - m_qs) as f64 * GC_OFFSET_RATE && novlp > GC_OFFSET_POS {
            return false;
        }
    }
    true
}

/// `gl_chain_graph` (`inter.cpp` 5785): the full graph-chain DP. `need_srt`
/// is `1` from `direct_gchain`, so the chains always arrive pre-sorted by the
/// caller and are sorted here.
#[allow(clippy::too_many_arguments)]
fn gl_chain_graph(
    index: &UlIndex,
    g: &Asg,
    lc: &mut Vec<MgLchain>,
    qlen: i64,
    bw: i64,
    diff_thre: f64,
    ng_diff_thre: f64,
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
) -> (i64, Vec<u64>) {
    if lc.is_empty() {
        return (0, Vec::new());
    }
    let lc_n = lc.len();
    let mut srt: Vec<u64> = vec![0; lc_n];
    let mut n_ext = 0usize;
    for (i, r) in lc.iter_mut().enumerate() {
        r.dist_pre = -1;
        // dist_pre -> parent in graph chain
        let isolated = (r.re as i64) < i64::from(g.sequences[(r.v >> 1) as usize].length)
            && r.rs > 0;
        if !isolated {
            n_ext += 1;
        }
        srt[i] = ((r.qe as u32 as u64) << 32) | i as u64 | ((u64::from(isolated)) << 63);
    }
    // `radix_sort_gfa64` then per-`qe` group by `qs`.
    srt.sort();
    let mut j = 0usize;
    for i in 1..=lc_n {
        if i == lc_n || (srt[i] >> 32) != (srt[j] >> 32) {
            if i - j > 1 {
                for x in j..i {
                    srt[x] = (srt[x] << 32 >> 32)
                        | (((lc[(srt[x] as u32) as usize].qs as u32 as u64) << 32));
                }
                srt[j..i].sort();
            }
            j = i;
        }
    }
    // `need_srt = 1`: the reorder always happens
    let sorted: Vec<MgLchain> = (0..lc_n).map(|i| lc[(srt[i] as u32) as usize]).collect();
    lc.copy_from_slice(&sorted);

    let mut m_idx: i64 = -1;
    let mut idx_out: Vec<u64> = Vec::new();
    // `if(ng_diff_thre >= 0)`: the linear pre-pass.
    if ng_diff_thre >= 0.0 {
        let (linear_idx, linear_out) = gl_chain_linear(
            g,
            lc,
            &index.unitig.unitigs,
            qlen,
            bw,
            ng_diff_thre,
            max_skip,
            max_iter,
            max_dis,
            n_ext,
        );
        if diff_thre < 0.0 {
            return (linear_idx, linear_out);
        }
        if linear_idx >= 0 && primary_chain_check(&linear_out, lc) {
            return (linear_idx, linear_out);
        }
    }

    let mut f = vec![0i64; lc_n];
    let mut p = vec![-1i64; lc_n];
    let mut t = vec![0u64; lc_n];
    let mut dst: Vec<MgPathDst> = Vec::new();
    let mut dst_done: Vec<u64> = Vec::new();
    let mut out: Vec<usize> = Vec::new();
    let mut path: Vec<MgPathv> = Vec::new();

    let mut plus: i64 = 0;
    let mut i = 0usize;
    while i < n_ext {
        let li = lc[i];
        let ui = set_ul_ov_t_by_mg_lchain(&li);
        let mm_ovlp = max_ovlp(g, li.v ^ 1);
        let mut x = ((i64::from(li.qs) + mm_ovlp) as f64 * diff_thre) as i64;
        if x < bw {
            x = bw;
        }
        x += i64::from(li.qs) + mm_ovlp;
        if x > qlen + 1 {
            x = qlen + 1;
        }
        let x = i64::from(find_mg_lchain_max(&lc[..i], (x + G_CHAIN_INDEL) as i32));

        let mut n_skip = 0i64;
        let mut is_f = false;
        let mut max_target_dist: i64 = -1;
        dst.clear();
        let mut jj = x;
        while jj >= 0 {
            let lj = lc[jj as usize];
            if lj.qs >= li.qs {
                jj -= 1;
                continue;
            }
            let target_dist = hc_target_len(g, &li, &lj);
            if target_dist < 0 {
                jj -= 1;
                continue;
            }
            let mut q = MgPathDst {
                v: lj.v ^ 1, // must be v^1 instead of v
                meta: jj as u32,
                qlen: li.qs - lj.qe, // might be negative
                target_dist: target_dist as i32,
                ..MgPathDst::default()
            };
            q.inner = false; // allow circles
            q.target_hash = 0;
            q.check_hash = false;
            if max_target_dist < target_dist {
                max_target_dist = target_dist;
            }
            dst.push(q);
            if t[jj as usize] == i as u64 {
                n_skip += 1;
                if n_skip > max_skip {
                    break;
                }
            }
            if p[jj as usize] >= 0 {
                t[p[jj as usize] as usize] = i as u64;
            }
            if !is_f && i64::from(lj.qe) + G_CHAIN_INDEL > i64::from(li.qs) {
                let uj = set_ul_ov_t_by_mg_lchain(&lj);
                let qo = infer_rovlp(&ui, &uj, &index.unitig.unitigs);
                if get_ecov_adv(g, li.v ^ 1, lj.v ^ 1, bw, N_GCHAIN_RATE, qo) {
                    is_f = true;
                    if n_skip > 0 {
                        n_skip -= 1;
                    }
                    if n_skip < (max_skip >> 1) {
                        n_skip = max_skip >> 1;
                    }
                }
            }
            jj -= 1;
        }

        // confirm reach-ability
        let mut max_f = i64::from(li.score);
        let mut max_j: i64 = -1;
        let mut max_d: i64 = -1;
        let mut max_hash: u32 = 0;
        let mut max_inner = false;
        if !dst.is_empty() {
            let mut max_target_dist = (max_target_dist as f64 * (1.0 + diff_thre)) as i64;
            if max_target_dist < bw {
                max_target_dist = bw;
            }
            hc_shortest_k(
                g,
                li.v ^ 1,
                &mut dst,
                max_target_dist as i32,
                15,
                &mut dst_done,
                &mut out,
                &mut path,
                true,
                false,
                0.0,
            );
            for dj in &dst {
                if dj.n_path == 0 {
                    continue; // unreachable
                }
                let sc = cal_gchain_sc(dj, &li, &f, bw, diff_thre, W_CHN_PEN_GAP);
                if sc == i32::MIN as i64 {
                    continue; // out of band
                }
                if sc > max_f {
                    max_f = sc;
                    max_j = dj.meta as i64;
                    max_d = dj.dist as i64;
                    max_hash = dj.hash;
                    max_inner = dj.inner;
                }
            }
        }
        if max_f < 0 {
            max_f = i64::from(li.score);
            max_j = -1;
            max_d = -1;
            max_hash = 0;
            max_inner = false;
        }

        f[i] = max_f;
        p[i] = max_j;
        // same time for gchain
        lc[i].dist_pre = max_d as i32;
        lc[i].hash_pre = max_hash;
        lc[i].inner_pre = max_inner;
        if max_f < plus {
            plus = max_f;
        }
        i += 1;
    }
    while i < lc_n {
        let li = lc[i];
        f[i] = i64::from(li.score);
        p[i] = -1;
        lc[i].dist_pre = -1;
        lc[i].hash_pre = 0;
        lc[i].inner_pre = false;
        if f[i] < plus {
            plus = f[i];
        }
        i += 1;
    }

    for i in 0..lc_n {
        f[i] -= plus;
        t[i] = ((f[i] as u64) << 32) + ((i as u64) << 1);
    }

    let mut sw: Vec<MgLchain> = Vec::with_capacity(lc_n);
    let mut u: Vec<u64> = Vec::with_capacity(lc_n);
    let mut n_v = 0usize;
    t.sort(); // `radix_sort_gfa64i` ascending, ties impossible
    let mut plus2: i64 = 0;
    for k in (0..lc_n).rev() {
        let n_v0 = n_v;
        let mut i = ((t[k]) & 0xFFFF_FFFF) as usize >> 1;
        while i != usize::MAX && t[i] & 1 == 0 {
            sw.push(lc[i]);
            n_v += 1;
            t[i] |= 1;
            i = if p[i] >= 0 { p[i] as usize } else { usize::MAX };
        }
        if n_v0 == n_v {
            continue;
        }
        let mut sc = if i == usize::MAX {
            (t[k] >> 32) as i64
        } else {
            (t[k] >> 32) as i64 - f[i]
        };
        if sc < plus2 {
            plus2 = sc;
        }
        let flag: u64 = if sc >= 0 {
            0x8000_0000_0000_0000
        } else {
            sc = -sc;
            0
        };
        u.push((((sc as u64) << 32) | (n_v - n_v0) as u64) | flag);
    }

    let mut m_sc: i64 = -1;
    let mut k = 0usize;
    for i in 0..u.len() {
        let entry = u[i];
        let mut sc = (entry >> 32) as i64;
        if entry & 0x8000_0000_0000_0000 == 0 {
            sc = -sc;
        }
        sc -= plus2;
        let count = (entry & 0xFFFF_FFFF) as usize;
        u[i] = (u[i] << 32 >> 32) | ((sc as u64) << 32);

        let k0 = k;
        for j in 0..count {
            lc[k] = sw[k0 + count - j - 1];
            k += 1;
        }
        if m_idx < 0 || m_sc < (u[i] >> 32) as i64 {
            m_idx = i as i64;
            m_sc = (u[i] >> 32) as i64;
        }
    }
    debug_assert_eq!(k, n_v);
    idx_out = u;
    (m_idx, idx_out)
}

/// `select_max_gchain` (`inter.cpp` 6310): filters overlapping chains and
/// materialises the survivors via `gen_gchain_track`.
#[allow(clippy::too_many_arguments)]
fn select_max_gchain(
    g: &Asg,
    idx: &[u64],
    e: &[MgLchain],
    raw_idx: &mut Vec<UlOv>,
    dst_done: &mut Vec<u64>,
    out: &mut Vec<usize>,
    res: &mut Vec<MgPathv>,
    gchains: &mut Vec<MgLchain>,
) -> usize {
    gchains.clear();
    if idx.is_empty() {
        return 0;
    }
    let mut a_n = 0usize;
    let raw_idx_n = raw_idx.len();
    for &entry in idx {
        let p = UlOv {
            qn: (entry >> 32) as u32,
            ts: a_n as u32,
            te: (a_n + entry as u32 as usize) as u32,
            qs: e[a_n].qs as u32,
            qe: e[a_n + (entry as u32 as usize) - 1].qe as u32,
            tn: 1, // (tn = 1) -> normal; (tn = 0) -> duplicated chain
            ..UlOv::default()
        };
        raw_idx.push(p);
        a_n += entry as u32 as usize;
    }
    let mut gb: Vec<UlOv> = raw_idx[raw_idx_n..].to_vec();
    // sort by scores (`radix_sort_ul_ov_srt_qn`), then reverse
    crate::container::ksort::radix_sort_by_key_32(&mut gb, &|entry: &UlOv| entry.qn);
    gb.reverse();
    let mut gb_n = gb.len();

    for k in 0..gb.len() {
        let om = u64::from(gb[k].qe - gb[k].qs);
        if gb[k].tn == 0 {
            continue;
        }
        let mut i = k as i64 - 1;
        while i >= 0 {
            let p = gb[i as usize];
            let m = gb[k];
            let ovlp = if m.qe.min(p.qe) > m.qs.max(p.qs) {
                m.qe.min(p.qe) - m.qs.max(p.qs)
            } else {
                0
            };
            if ovlp == 0 {
                i -= 1;
                continue;
            }
            let min_sc = p.qn.min(m.qn);
            let max_sc = p.qn.max(m.qn);
            let ok = (p.qe - p.qs).max(om as u32) as u64;
            if u64::from(min_sc) < u64::from(max_sc) * 98 / 100 {
                break;
            }
            if u64::from(ovlp) > GC_OFFSET_POS as u64
                && u64::from(min_sc) > u64::from(max_sc) * 98 / 100
                && ovlp as f64 > ok as f64 * 0.8
            {
                gb[k].tn = 0;
                gb[i as usize].tn = 0;
            }
            i -= 1;
        }
        let mut i = k + 1;
        while i < gb.len() {
            let p = gb[i];
            let m = gb[k];
            let ovlp = if m.qe.min(p.qe) > m.qs.max(p.qs) {
                m.qe.min(p.qe) - m.qs.max(p.qs)
            } else {
                0
            };
            if ovlp == 0 {
                i += 1;
                continue;
            }
            let min_sc = p.qn.min(m.qn);
            let max_sc = p.qn.max(m.qn);
            let ok = (p.qe - p.qs).max(om as u32) as u64;
            if u64::from(min_sc) < u64::from(max_sc) * 98 / 100 {
                break;
            }
            if u64::from(ovlp) > GC_OFFSET_POS as u64
                && u64::from(min_sc) > u64::from(max_sc) * 98 / 100
                && ovlp as f64 > ok as f64 * 0.8
            {
                gb[k].tn = 0;
                gb[i].tn = 0;
            }
            i += 1;
        }
    }

    gb.retain(|entry| entry.tn != 0);
    gb_n = gb.len();

    let mut n_mchain = 0usize;
    for k in 0..gb_n {
        let om = u64::from(gb[k].qe - gb[k].qs);
        let mut i = 0usize;
        while i < n_mchain {
            let p = gb[i];
            let m = gb[k];
            let ovlp = if m.qe.min(p.qe) > m.qs.max(p.qs) {
                m.qe.min(p.qe) - m.qs.max(p.qs)
            } else {
                0
            };
            if ovlp == 0 {
                i += 1;
                continue;
            }
            let ok = u64::from(p.qe - p.qs);
            if u64::from(ovlp) > ok / 10 || u64::from(ovlp) > om / 10 {
                break;
            }
            i += 1;
        }
        if i < n_mchain {
            continue;
        }
        gb[n_mchain] = gb[k];
        n_mchain += 1;
    }
    gb.truncate(n_mchain);

    if !gb.is_empty() {
        for item in &gb {
            let rs = item.ts as usize;
            let re = item.te as usize;
            // the group header is pushed first; the materialised entries go
            // right behind it
            gchains.push(MgLchain {
                v: u32::MAX,
                qs: item.qs as i32,
                qe: item.qe as i32,
                rs: item.ts as i32,
                re: item.te as i32,
                ..MgLchain::default()
            });
            gen_gchain_track(g, &e[rs..re], dst_done, out, res);
            let res_n = res.len();
            let base = gchains.len() - 1;
            gchains.resize(base + 1 + res_n, MgLchain::default());
            for (i, path) in res.iter().enumerate() {
                if path.v == u32::MAX {
                    gchains[base + 1 + i] = e[path.pre as usize + rs];
                    gchains[base + 1 + i].dist_pre = path.d as i32;
                } else {
                    gchains[base + 1 + i].v = path.v;
                    gchains[base + 1 + i].off = -1;
                    gchains[base + 1 + i].dist_pre = path.d as i32;
                    // the graph-chained nodes are fully covered
                    gchains[base + 1 + i].rs = 0;
                    gchains[base + 1 + i].re =
                        g.sequences[(path.v >> 1) as usize].length as i32;
                }
            }
            gchains[base].cnt = res_n as i32;
        }
    }

    raw_idx.truncate(raw_idx_n);
    n_mchain
}

/// `get_r_offset` (`inter.cpp` 4501).
fn get_r_offset(g: &Asg, x: &MgLchain) -> (i64, i64, i64, i64) {
    let qs = i64::from(x.qs);
    let qe = i64::from(x.qe);
    let (rs, re) = if x.score & 1 == 0 {
        (i64::from(x.rs) + i64::from(x.off), i64::from(x.re) + i64::from(x.off))
    } else {
        let len = i64::from(g.sequences[(x.score >> 1) as usize].length);
        (
            i64::from(x.off) + len - i64::from(x.re),
            i64::from(x.off) + len - i64::from(x.rs),
        )
    };
    (rs, re, qs, qe)
}

/// `cal_qext_coor` (`inter.cpp` 4779).
fn cal_qext_coor(pr: i64, ar: i64, pq: i64, aq: i64, r_off: i64) -> i64 {
    let q_off;
    if r_off >= pr && r_off <= ar && ar >= pr && aq >= pq {
        q_off = pq
            + i64::from(crate::assembly::ul::get_offset_adjust(
                (r_off - pr) as u32,
                (ar - pr) as u32,
                (aq - pq) as u32,
            ));
    } else {
        let pd = if r_off >= pr { r_off - pr } else { pr - r_off };
        let ad = if r_off >= ar { r_off - ar } else { ar - r_off };
        q_off = if ad <= pd { aq } else { pq };
    }
    q_off
}

/// `extract_rovlp_by_ug` (`inter.cpp` 15272).
fn extract_rovlp_by_ug(p: &UtgCt, o: &MgLchain, chains: &mut Vec<MgLchain>, t_off: i64) -> bool {
    let rs = i64::from(p.s).max(i64::from(o.rs));
    let re = i64::from(p.e).min(i64::from(o.re));
    if rs > re {
        return false;
    }
    let mut v = (p.x >> 1) << 1;
    v += u32::from((o.v & 1) != (p.x & 1));
    chains.push(MgLchain {
        v,
        rs: rs as i32,
        re: re as i32,
        off: t_off as i32,
        hash_pre: u32::MAX,
        dist_pre: -1,
        qs: -1,
        qe: -1,
        ..MgLchain::default()
    });
    true
}

/// `gl_ug2rg_gen` (`inter.cpp` 15369).
fn gl_ug2rg_gen(
    rg: &Asg,
    ug_g: &Asg,
    r_cl: &UlReadVec,
    raw_idx: &[UlOv],
    raw_chn: &[UlOv],
    unitigs: &[MaUtg],
    read_lengths: &[u64],
    uo: &MgLchain,
    res: &mut Vec<MgLchain>,
    t_off: i64,
) {
    let res_n0 = res.len();
    let unitig = &unitigs[(uo.v >> 1) as usize];
    let rs = uo.rs;
    let re = uo.re;
    let mut l = 0i64;
    for (i, read) in unitig.reads.iter().enumerate() {
        let p = UtgCt {
            x: (read >> 32) as u32,
            s: l as u32,
            e: (l + read_lengths[(read >> 33) as usize] as i64) as u32,
        };
        l += (read & 0xFFFF_FFFF) as i64;
        if p.e as i64 <= i64::from(rs) {
            continue;
        }
        if p.s as i64 >= i64::from(re) {
            break;
        }
        if !extract_rovlp_by_ug(&p, uo, res, t_off) {
            panic!(
                "[gl_ug2rg] empty overlap: uo v={} rs={} re={} off={} | read {} p.x={} p.s={} p.e={}",
                uo.v, uo.rs, uo.re, uo.off, i, p.x, p.s, p.e
            );
        }
        let last = res.len() - 1;
        res[last].score = uo.v as i32;
        res[last].cnt = i as i32;
    }

    let a_n = res.len() - res_n0;
    if uo.v & 1 != 0 {
        res[res_n0..].reverse();
    }

    let unitig = &unitigs[(uo.v >> 1) as usize];
    update_existing_anchors(rg, ug_g, r_cl, raw_idx, raw_chn, uo, unitig, res, res_n0);
}

/// `update_existing_anchors` (`inter.cpp` 15287).
#[allow(clippy::too_many_arguments)]
fn update_existing_anchors(
    rg: &Asg,
    ug_g: &Asg,
    r_cl: &UlReadVec,
    raw_idx: &[UlOv],
    raw_chn: &[UlOv],
    uo: &MgLchain,
    unitig: &MaUtg,
    res: &mut Vec<MgLchain>,
    res_n0: usize,
) {
    let mut z: i64 = -1;
    if uo.off >= 0 {
        z = raw_idx[uo.off as usize].qn as i64;
    }
    let mut m = res.len() as i64 - 1;
    while z >= 0 {
        let block = &r_cl.bb[raw_chn[z as usize].qn as usize];
        let zv = (block.hid << 1) | u32::from(block.rev);
        if !rg.sequences[(zv >> 1) as usize].deleted {
            let mut midx: i64 = -1;
            let mut mdif: i32 = -1;
            let mbeg = m;
            while m >= res_n0 as i64 {
                if zv == res[m as usize].v {
                    if raw_chn[z as usize].ts as i32 == res[m as usize].rs
                        && raw_chn[z as usize].te as i32 == res[m as usize].re
                    {
                        midx = m;
                        mdif = 0;
                        break;
                    } else {
                        let ovlp = if (raw_chn[z as usize].te as i32).min(res[m as usize].re)
                            > (raw_chn[z as usize].ts as i32).max(res[m as usize].rs)
                        {
                            (raw_chn[z as usize].te as i32).min(res[m as usize].re)
                                - (raw_chn[z as usize].ts as i32).max(res[m as usize].rs)
                        } else {
                            0
                        };
                        let novlp = (raw_chn[z as usize].te as i32 - raw_chn[z as usize].ts as i32
                            - ovlp)
                            + (res[m as usize].re - res[m as usize].rs - ovlp);
                        if midx == -1 || mdif > novlp {
                            midx = m;
                            mdif = novlp;
                        }
                    }
                }
                m -= 1;
            }
            if mdif != 0 {
                m = res.len() as i64 - 1;
                while m > mbeg {
                    if zv == res[m as usize].v {
                        if raw_chn[z as usize].ts as i32 == res[m as usize].rs
                            && raw_chn[z as usize].te as i32 == res[m as usize].re
                        {
                            midx = m;
                            mdif = 0;
                            break;
                        } else {
                            let ovlp = if (raw_chn[z as usize].te as i32).min(res[m as usize].re)
                                > (raw_chn[z as usize].ts as i32).max(res[m as usize].rs)
                            {
                                (raw_chn[z as usize].te as i32).min(res[m as usize].re)
                                    - (raw_chn[z as usize].ts as i32).max(res[m as usize].rs)
                            } else {
                                0
                            };
                            let novlp = (raw_chn[z as usize].te as i32
                                - raw_chn[z as usize].ts as i32
                                - ovlp)
                                + (res[m as usize].re - res[m as usize].rs - ovlp);
                            if midx == -1 || mdif > novlp {
                                midx = m;
                                mdif = novlp;
                            }
                        }
                    }
                    m -= 1;
                }
            }
            m = midx;
            debug_assert!(m >= res_n0 as i64);
            res[m as usize].qs = raw_chn[z as usize].qs as i32;
            res[m as usize].qe = raw_chn[z as usize].qe as i32;
            res[m as usize].dist_pre = raw_chn[z as usize].qn as i32;
            // the co-linearity check against the neighbours of `m`
            let left: (i64, i64) = if m > 0 {
                let (rs, re, _, _) = get_r_offset(ug_g, &res[(m - 1) as usize]);
                (rs, re)
            } else {
                (-1, -1)
            };
            let right: (i64, i64) = if m + 1 < res.len() as i64 {
                let (rs, re, _, _) = get_r_offset(ug_g, &res[(m + 1) as usize]);
                (rs, re)
            } else {
                (
                    i64::from(unitig.length) + 1,
                    i64::from(unitig.length) + 1,
                )
            };
            if (raw_chn[z as usize].ts as i64) >= left.0
                && (raw_chn[z as usize].ts as i64) <= right.0
                && (raw_chn[z as usize].te as i64) >= left.1
                && (raw_chn[z as usize].te as i64) <= right.1
            {
                res[m as usize].rs = raw_chn[z as usize].ts as i32;
                res[m as usize].re = raw_chn[z as usize].te as i32;
            }
        }
        z = if raw_chn[z as usize].tn == u32::MAX {
            -1
        } else {
            raw_chn[z as usize].tn as i64
        };
    }
}

/// `flat_rovlp_chain_sc` (`inter.cpp` 15575).
fn flat_rovlp_chain_sc(g: &Asg, li: &MgLchain, lj: &MgLchain) -> i64 {
    let (its, ite, iqs, iqe) = get_r_offset(g, li);
    let (jts, jte, jqs, jqe) = get_r_offset(g, lj);
    if jqs <= iqs && jqe <= iqe && jts <= its && jte <= ite {
        if jqs == iqs && jqe == iqe && jts == its && jte == ite {
            return i32::MIN as i64;
        }
        let ovlp = if jqe.min(iqe) > jqs.max(iqs) {
            jqe.min(iqe) - jqs.max(iqs)
        } else {
            0
        };
        return iqe - iqs - ovlp;
    }
    i32::MIN as i64
}

/// `flat_rovlp_chain` (`inter.cpp` 15588).
fn flat_rovlp_chain(g: &Asg, x: &mut [MgLchain], max_skip: i64, max_iter: i64, max_dis: i64) -> i64 {
    let x_n = x.len();
    if x_n == 0 {
        return 0;
    }
    let mut id: Vec<i64> = Vec::with_capacity(x_n);
    let mut cl = 1i64;
    for (i, item) in x.iter().enumerate() {
        if item.qs >= 0 {
            if cl == 1 && !id.is_empty() {
                let li = item;
                let lj = &x[id[id.len() - 1] as usize];
                let sc = flat_rovlp_chain_sc(g, li, lj);
                if sc == i32::MIN as i64 {
                    cl = 0;
                }
            }
            id.push(i as i64);
        }
    }
    let a_n = id.len();
    if a_n == 0 {
        return 0;
    }
    if cl == 1 {
        return a_n as i64;
    }

    let mut f = vec![0i64; a_n];
    let mut p = vec![-1i64; a_n];
    let mut t: Vec<u64> = vec![0; a_n];
    let mut c_n = vec![0i32; a_n];
    let mut st;
    let mut max_ii: i64 = -1;
    let mut tot_sc: i64 = i32::MIN as i64;
    let mut tot_n: i64 = i32::MIN as i64;
    let mut tot_i: i64 = -1;
    for i in 0..a_n {
        let li = x[id[i] as usize];
        let mut mm_sc = i64::from(li.qe - li.qs);
        let mut mm_n: i64 = 1;
        let mut mm_idx: i64 = -1;
        let mut n_skip = 0i64;
        st = if (i as i64) < max_iter { 0 } else { i as i64 - max_iter };
        let mut end_j: i64;
        let mut jj = i as i64 - 1;
        while jj >= st {
            let lj = &x[id[jj as usize] as usize];
            let sc = flat_rovlp_chain_sc(g, &li, lj);
            if sc == i32::MIN as i64 {
                jj -= 1;
                continue;
            }
            let sc = sc + f[jj as usize];
            let sn = i64::from(c_n[jj as usize]) + 1;
            if (sc > mm_sc) || (sc == mm_sc && sn > mm_n) {
                mm_sc = sc;
                mm_idx = jj;
                mm_n = sn;
                if n_skip > 0 {
                    n_skip -= 1;
                }
            } else if t[jj as usize] == i as u64 {
                n_skip += 1;
                if n_skip > max_skip {
                    break;
                }
            }
            if p[jj as usize] >= 0 {
                t[p[jj as usize] as usize] = i as u64;
            }
            jj -= 1;
        }
        end_j = jj;
        if max_ii < 0
            || i64::from(x[id[i] as usize].qe)
                > i64::from(x[id[max_ii as usize] as usize].qe) + max_dis
        {
            let mut max = i32::MIN as i64;
            let mut max_n: i64 = i32::MIN as i64;
            max_ii = -1;
            let mut j2 = i as i64 - 1;
            while j2 >= st
                && i64::from(x[id[i] as usize].qe)
                    <= max_dis + i64::from(x[id[j2 as usize] as usize].qe)
            {
                if max < f[j2 as usize]
                    || (max == f[j2 as usize] && max_n < i64::from(c_n[j2 as usize]))
                {
                    max = f[j2 as usize];
                    max_n = i64::from(c_n[j2 as usize]);
                    max_ii = j2;
                }
                j2 -= 1;
            }
        }

        if max_ii >= 0 && max_ii < end_j {
            let lj = &x[id[max_ii as usize] as usize];
            let sc = flat_rovlp_chain_sc(g, &li, lj);
            if sc != i32::MIN as i64 {
                let sc = sc + f[max_ii as usize];
                let sn = i64::from(c_n[max_ii as usize]) + 1;
                if (sc > mm_sc) || (sc == mm_sc && sn > mm_n) {
                    mm_sc = sc;
                    mm_idx = max_ii;
                    mm_n = sn;
                }
            }
        }

        f[i] = mm_sc;
        p[i] = mm_idx;
        c_n[i] = mm_n as i32;
        if max_ii < 0
            || (i64::from(x[id[i] as usize].qe)
                <= max_dis + i64::from(x[id[max_ii as usize] as usize].qe)
                && f[max_ii as usize] < f[i])
        {
            max_ii = i as i64;
        }
        if tot_sc < f[i] || (tot_sc == f[i] && tot_n < i64::from(c_n[i])) {
            tot_sc = f[i];
            tot_n = i64::from(c_n[i]);
            tot_i = i as i64;
        }
    }

    let mut cl: i64 = 0;
    let mut i = tot_i;
    while i >= 0 {
        id[i as usize] = -1;
        i = p[i as usize];
        cl += 1;
    }

    if cl < a_n as i64 {
        for &index in id.iter() {
            if index < 0 {
                continue;
            }
            let li = &mut x[index as usize];
            li.dist_pre = -1;
            li.qs = -1;
            li.qe = -1;
        }
    }
    cl
}

/// `update_rovlp_chain_qse` (`inter.cpp` 15495).
fn update_rovlp_chain_qse(
    g: &Asg,
    sidx: i64,
    eidx: i64,
    a: &mut [MgLchain],
    qlen: i64,
) {
    if eidx - sidx <= 1 {
        return;
    }
    debug_assert!(sidx >= 0 || (eidx as usize) < a.len());
    let (left_r0, left_r1, left_q0, left_q1) = if sidx >= 0 {
        get_r_offset(g, &a[sidx as usize])
    } else {
        get_r_offset(g, &a[0])
    };
    let (right_r0, right_r1, right_q0, right_q1) = if (eidx as usize) < a.len() {
        get_r_offset(g, &a[eidx as usize])
    } else {
        get_r_offset(g, &a[a.len() - 1])
    };
    debug_assert!(
        (left_q0 >= 0 && left_q1 >= 0) || (right_q0 >= 0 && right_q1 >= 0)
    );

    let mut left_r0 = left_r0;
    let mut left_r1 = left_r1;
    let mut left_q0 = left_q0;
    let mut left_q1 = left_q1;

    if left_q0 >= 0 && left_q1 >= 0 && right_q0 >= 0 && right_q1 >= 0 {
        for i in sidx + 1..eidx {
            let (rs, re, _, _) = get_r_offset(g, &a[i as usize]);
            a[i as usize].qe =
                cal_qext_coor(left_r1, right_r1, left_q1, right_q1, re) as i32;
            debug_assert!(a[i as usize].qe >= 0 && i64::from(a[i as usize].qe) <= qlen);
            a[i as usize].qs = if a[i as usize].qe as i64 <= right_q0 {
                cal_qext_coor(left_r0, re, left_q0, i64::from(a[i as usize].qe), rs)
            } else {
                cal_qext_coor(left_r0, right_r0, left_q0, right_q0, rs)
            } as i32;
            debug_assert!(a[i as usize].qs >= 0 && i64::from(a[i as usize].qs) <= qlen);
            if a[i as usize].qs > a[i as usize].qe {
                std::mem::swap(&mut a[i as usize].qs, &mut a[i as usize].qe);
            }
            left_q0 = i64::from(a[i as usize].qs);
            left_q1 = i64::from(a[i as usize].qe);
            left_r0 = rs;
            left_r1 = re;
        }
    }

    if right_q0 < 0 || right_q1 < 0 {
        for i in sidx + 1..eidx {
            let (rs, re, _, _) = get_r_offset(g, &a[i as usize]);
            a[i as usize].qe =
                cal_qext_coor(left_r1, re, left_q1, left_q1 + re - left_r1, re) as i32;
            if a[i as usize].qe < 0 {
                a[i as usize].qe = 0;
            }
            if i64::from(a[i as usize].qe) > qlen {
                a[i as usize].qe = qlen as i32;
            }
            a[i as usize].qs = if a[i as usize].qe as i64 <= left_q1 {
                cal_qext_coor(left_r0, re, left_q0, i64::from(a[i as usize].qe), rs)
            } else {
                cal_qext_coor(left_r0, left_r1, left_q0, left_q1, rs)
            } as i32;
            debug_assert!(a[i as usize].qs >= 0 && i64::from(a[i as usize].qs) <= qlen);
            if a[i as usize].qs > a[i as usize].qe {
                std::mem::swap(&mut a[i as usize].qs, &mut a[i as usize].qe);
            }
            left_q0 = i64::from(a[i as usize].qs);
            left_q1 = i64::from(a[i as usize].qe);
            left_r0 = rs;
            left_r1 = re;
        }
    }

    if left_q0 < 0 || left_q1 < 0 {
        let mut right_r0 = right_r0;
        let mut right_r1 = right_r1;
        let mut right_q0 = right_q0;
        let mut right_q1 = right_q1;
        let mut i = eidx - 1;
        while i > sidx {
            let (rs, re, _, _) = get_r_offset(g, &a[i as usize]);
            a[i as usize].qe = cal_qext_coor(right_r0, right_r1, right_q0, right_q1, re) as i32;
            debug_assert!(a[i as usize].qe >= 0 && i64::from(a[i as usize].qe) <= qlen);
            a[i as usize].qs =
                cal_qext_coor(rs, right_r0, right_q0 - (right_r0 - rs), right_q0, rs) as i32;
            if a[i as usize].qs < 0 {
                a[i as usize].qs = 0;
            }
            if i64::from(a[i as usize].qs) > qlen {
                a[i as usize].qs = qlen as i32;
            }
            if a[i as usize].qs > a[i as usize].qe {
                std::mem::swap(&mut a[i as usize].qs, &mut a[i as usize].qe);
            }
            right_q0 = i64::from(a[i as usize].qs);
            right_q1 = i64::from(a[i as usize].qe);
            right_r0 = rs;
            right_r1 = re;
            i -= 1;
        }
    }
}

/// `gen_rovlp_chain_by_ul` (`inter.cpp` 15675).
#[allow(clippy::too_many_arguments)]
fn gen_rovlp_chain_by_ul(
    rg: &Asg,
    r_cl: &UlReadVec,
    index: &UlIndex,
    read_lengths: &[u64],
    raw_idx: &[UlOv],
    raw_chn: &[UlOv],
    a: &[MgLchain],
    res: &mut Vec<MgLchain>,
    dp_max_skip: i64,
    dp_max_iter: i64,
    dp_max_dis: i64,
) {
    if a.is_empty() {
        return;
    }
    let res_n0 = res.len();
    let ug_g = &index.unitig.graph;
    let mut l = i64::from(ug_g.sequences[(a[0].v >> 1) as usize].length);
    for item in a {
        l -= i64::from(ug_g.sequences[(item.v >> 1) as usize].length);
        gl_ug2rg_gen(
            rg,
            &index.unitig.graph,
            r_cl,
            raw_idx,
            raw_chn,
            &index.unitig.unitigs,
            read_lengths,
            item,
            res,
            l,
        );
        l += i64::from(ug_g.sequences[(item.v >> 1) as usize].length) + i64::from(item.dist_pre);
    }
    let x_n = res.len() - res_n0;
    let fn_flat = flat_rovlp_chain(ug_g, &mut res[res_n0..], dp_max_skip, dp_max_iter, dp_max_dis);
    if fn_flat > 0 {
        let mut tt = 0usize;
        for k in 0..x_n {
            let idx = res_n0 + k;
            if k > 0 {
                res[idx].hash_pre = (k - 1 + res_n0) as u32;
            } else {
                res[idx].hash_pre = u32::MAX;
            }
            if res[idx].qs >= 0 {
                tt += 1;
            }
        }
        let mut l: i64 = -1;
        for k in 0..=x_n {
            if k == x_n || res[res_n0 + k].qs >= 0 {
                if k as i64 - l > 1 {
                    update_rovlp_chain_qse(
                        ug_g,
                        res_n0 as i64 + l,
                        res_n0 as i64 + k as i64,
                        res,
                        i64::from(r_cl.rlen),
                    );
                }
                l = k as i64;
            }
        }
        debug_assert!(tt > 0);
    } else {
        res.truncate(res_n0);
    }
}

/// `convert_mg_lchain_t` (`inter.cpp` 15703).
fn convert_mg_lchain_t(p: &UtgCt, o: &mut MgLchain) {
    let mut rs = i64::from(p.s).max(i64::from(o.rs));
    let mut re = i64::from(p.e).min(i64::from(o.re));
    debug_assert!(rs < re);
    if p.x & 1 == 0 {
        o.rs = (rs - i64::from(p.s)) as i32;
        o.re = (re - i64::from(p.s)) as i32;
    } else {
        // both expressions read the clipped values, as in the reference
        let (orig_rs, orig_re) = (rs, re);
        rs = i64::from(p.e) - orig_re;
        re = i64::from(p.e) - orig_rs;
        o.rs = rs as i32;
        o.re = re as i32;
    }
}

/// `renew_mg_lchains` (`inter.cpp` 15716).
fn renew_mg_lchains(unitigs: &[MaUtg], read_lengths: &[u64], a: &mut [MgLchain]) {
    if a.is_empty() {
        return;
    }
    let rev = a[0].score & 1;
    let u = &unitigs[(a[0].score >> 1) as usize];
    if rev == 0 {
        let mut k = 0usize;
        while k < a.len() {
            let mut l = 0i64;
            let mut i = 0usize;
            while i < u.reads.len() {
                if i as i64 == i64::from(a[k].cnt) {
                    break;
                }
                l += (u.reads[i] & 0xFFFF_FFFF) as i64;
                i += 1;
            }
            debug_assert!(i < u.reads.len());
            while k < a.len() {
                if i64::from(a[k].cnt) != i as i64 {
                    break;
                }
                debug_assert_eq!((a[k].v >> 1) as usize, (u.reads[i] >> 33) as usize);
                let p = UtgCt {
                    x: (u.reads[i] >> 32) as u32,
                    s: l as u32,
                    e: (l + read_lengths[(u.reads[i] >> 33) as usize] as i64) as u32,
                };
                convert_mg_lchain_t(&p, &mut a[k]);
                l += (u.reads[i] & 0xFFFF_FFFF) as i64;
                k += 1;
                i += 1;
            }
        }
    } else {
        let mut k = a.len() as i64 - 1;
        while k >= 0 {
            let mut l = 0i64;
            let mut i = 0usize;
            while i < u.reads.len() {
                if i as i64 == i64::from(a[k as usize].cnt) {
                    break;
                }
                l += (u.reads[i] & 0xFFFF_FFFF) as i64;
                i += 1;
            }
            debug_assert!(i < u.reads.len());
            while k >= 0 {
                if i64::from(a[k as usize].cnt) != i as i64 {
                    break;
                }
                debug_assert_eq!((a[k as usize].v >> 1) as usize, (u.reads[i] >> 33) as usize);
                let p = UtgCt {
                    x: (u.reads[i] >> 32) as u32,
                    s: l as u32,
                    e: (l + read_lengths[(u.reads[i] >> 33) as usize] as i64) as u32,
                };
                convert_mg_lchain_t(&p, &mut a[k as usize]);
                l += (u.reads[i] & 0xFFFF_FFFF) as i64;
                k -= 1;
                i += 1;
            }
        }
    }
}

/// `dd_ul_vec_t` (`inter.cpp` 15812).
fn dd_ul_vec_t(
    rg: &Asg,
    index: &UlIndex,
    read_lengths: &[u64],
    a: &mut [MgLchain],
    rch: &mut UlReadVec,
) {
    let mut l = 0usize;
    for k in 1..=a.len() {
        if k == a.len() || a[k].score != a[l].score {
            renew_mg_lchains(&index.unitig.unitigs, read_lengths, &mut a[l..k]);
            l = k;
        }
    }

    for block in rch.bb.iter_mut() {
        block.pidx = 0x0fff_fffe;
        block.aidx = u32::MAX;
        block.pdis = u32::MAX;
    }

    for (k, item) in a.iter().enumerate() {
        if item.dist_pre >= 0 {
            // not a new alignment
            let z = &mut rch.bb[item.dist_pre as usize];
            debug_assert_eq!(z.hid, item.v >> 1);
            debug_assert_eq!(z.rev, item.v & 1 != 0);
            if z.qs as i32 == item.qs
                && z.qe as i32 == item.qe
                && z.ts as i32 == item.rs
                && z.te as i32 == item.re
            {
                z.pidx = k as u32;
                z.pchain = 1;
            } else {
                let mut tt: i64 = 0;
                let mut novlp: i64 = 0;
                let ovlp = if (z.qe as i32).min(item.qe) > (z.qs as i32).max(item.qs) {
                    (z.qe as i32).min(item.qe) - (z.qs as i32).max(item.qs)
                } else {
                    0
                };
                tt += i64::from(ovlp);
                novlp += i64::from(item.qe) - i64::from(item.qs) - i64::from(ovlp)
                    + i64::from(z.qe)
                    - i64::from(z.qs)
                    - i64::from(ovlp);

                let ovlp = if (z.te as i32).min(item.re) > (z.ts as i32).max(item.rs) {
                    (z.te as i32).min(item.re) - (z.ts as i32).max(item.rs)
                } else {
                    0
                };
                tt += i64::from(ovlp);
                novlp += i64::from(item.re) - i64::from(item.rs) - i64::from(ovlp)
                    + i64::from(z.te)
                    - i64::from(z.ts)
                    - i64::from(ovlp);
                // `novlp > (tt*0.01)` compares against a double in the reference
                if novlp > 8 || novlp as f64 > tt as f64 * 0.01 {
                    rch.bb.push(UcBlock {
                        hid: item.v >> 1,
                        rev: item.v & 1 != 0,
                        pchain: 2,
                        base: false,
                        el: true,
                        qs: item.qs as u32,
                        qe: item.qe as u32,
                        te: item.re as u32,
                        ts: item.rs as u32,
                        pidx: k as u32,
                        aidx: u32::MAX,
                        pdis: u32::MAX,
                    });
                } else {
                    z.qs = item.qs as u32;
                    z.qe = item.qe as u32;
                    z.te = item.re as u32;
                    z.ts = item.rs as u32;
                    z.pidx = k as u32;
                    z.pchain = 1;
                }
            }
        } else {
            rch.bb.push(UcBlock {
                hid: item.v >> 1,
                rev: item.v & 1 != 0,
                pchain: 2,
                base: false,
                el: true,
                qs: item.qs as u32,
                qe: item.qe as u32,
                te: item.re as u32,
                ts: item.rs as u32,
                pidx: k as u32,
                aidx: u32::MAX,
                pdis: u32::MAX,
            });
        }
    }

    // drop the stale blocks
    let mut k_out = 0usize;
    for i in 0..rch.bb.len() {
        if rch.bb[i].pidx == 0x0fff_fffe {
            continue;
        }
        rch.bb[k_out] = rch.bb[i];
        k_out += 1;
    }
    rch.bb.truncate(k_out);

    crate::container::ksort::radix_sort_by_key_32(&mut rch.bb, &|block: &UcBlock| block.qe);
    let mut l2 = 0usize;
    for k in 1..=rch.bb.len() {
        if k == rch.bb.len() || rch.bb[k].qe != rch.bb[l2].qe {
            if k - l2 > 1 {
                crate::container::ksort::radix_sort_by_key_32(
                    &mut rch.bb[l2..k],
                    &|block: &UcBlock| block.qs,
                );
            }
            while l2 < k {
                if rch.bb[l2].pidx == 0x0fff_fffe {
                    rch.bb[l2].pidx = u32::MAX;
                } else {
                    a[rch.bb[l2].pidx as usize].dist_pre = l2 as i32;
                }
                l2 += 1;
            }
        }
    }

    // rebuild the pidx/pdis/aidx links and the full-coverage flag
    let mut k: i64 = -1;
    for i in 0..rch.bb.len() {
        if rch.bb[i].pidx == u32::MAX {
            continue;
        }
        if k < 0 {
            k = i as i64; // in case there is only one UL-to-HiFi alignment
        }
        if a[rch.bb[i].pidx as usize].hash_pre == u32::MAX {
            rch.bb[i].pidx = u32::MAX;
            continue;
        }
        let new_pidx = a[a[rch.bb[i].pidx as usize].hash_pre as usize].dist_pre as u32;
        rch.bb[i].pidx = new_pidx;
        rch.bb[new_pidx as usize].aidx = i as u32;
        let tt = crate::assembly::ul::g_adjacent_dis(
            rg,
            ((rch.bb[i].hid << 1) | u32::from(rch.bb[i].rev)) ^ 1,
            ((rch.bb[rch.bb[i].pidx as usize].hid << 1)
                | u32::from(rch.bb[rch.bb[i].pidx as usize].rev))
                ^ 1,
        );
        if tt >= 0 {
            rch.bb[i].pdis = tt as u32;
        } else {
            let p = &a[a[rch.bb[i].pidx as usize].hash_pre as usize];
            let c = &a[rch.bb[i].pidx as usize];
            let ug_g = &index.unitig.graph;
            let rs = i64::from(p.off) + i64::from(ug_g.sequences[(p.score >> 1) as usize].length);
            let re = i64::from(c.off) + i64::from(ug_g.sequences[(c.score >> 1) as usize].length);
            rch.bb[i].pdis = if re >= rs {
                (re - rs) as u32
            } else {
                u32::MAX
            };
        }
        k = i as i64;
    }

    // make sure if this UL read has been done
    let mut sp = u32::MAX;
    let mut ep = u32::MAX;
    let mut l3: i64 = 0;
    let mut kk = k;
    while kk >= 0 {
        if sp == u32::MAX || rch.bb[kk as usize].qe <= sp {
            if sp != u32::MAX {
                l3 += i64::from(ep - sp);
            }
            sp = rch.bb[kk as usize].qs;
            ep = rch.bb[kk as usize].qe;
        } else {
            sp = sp.min(rch.bb[kk as usize].qs);
        }
        kk = if rch.bb[kk as usize].pidx == u32::MAX {
            -1
        } else {
            rch.bb[kk as usize].pidx as i64
        };
    }
    if sp != u32::MAX {
        l3 += i64::from(ep - sp);
    }
    if l3 == i64::from(rch.rlen) {
        rch.dd = 1;
    }
}

/// `update_ul_vec_t` (`inter.cpp` 16028).
#[allow(clippy::too_many_arguments)]
fn update_ul_vec_t(
    rg: &Asg,
    index: &UlIndex,
    read_lengths: &[u64],
    raw_idx: &[UlOv],
    raw_chn: &[UlOv],
    uc: &[MgLchain],
    rch: &mut UlReadVec,
    dp_max_skip: i64,
    dp_max_iter: i64,
    dp_max_dis: i64,
) {
    let mut ucn = 0usize;
    let mut swap: Vec<MgLchain> = Vec::new();
    while ucn < uc.len() {
        let ix = &uc[ucn];
        debug_assert_eq!(ix.v, u32::MAX);
        let count = ix.cnt as usize;
        gen_rovlp_chain_by_ul(
            rg,
            rch,
            index,
            read_lengths,
            raw_idx,
            raw_chn,
            &uc[ucn + 1..ucn + 1 + count],
            &mut swap,
            dp_max_skip,
            dp_max_iter,
            dp_max_dis,
        );
        ucn += count + 1;
    }
    dd_ul_vec_t(rg, index, read_lengths, &mut swap, rch);
}

/// `refine_rid_chain` (`inter.cpp` 17420).
fn refine_rid_chain(rg: &Asg, rch: &UlReadVec) -> u32 {
    if rch.bb.len() == 1 && rch.bb[0].base {
        return 1; // no alignment
    }
    if rch.bb.is_empty() {
        return 1; // no alignment
    }
    let mut m = 0u64;
    let mut nc = 0u64;
    let mut c0 = u32::MAX;
    let mut i = rch.bb.len() as u32 - 1;
    while i != u32::MAX {
        m += 1;
        c0 = i;
        if !rg.sequences[rch.bb[i as usize].hid as usize].deleted {
            nc += 1;
        }
        i = rch.bb[i as usize].pidx;
    }
    if m == rch.bb.len() as u64 {
        // only one chain
        if nc == m || nc == 0 {
            return 1; // all alignments are non-contained/contained
        }
        let mut cc = 0u64;
        let mut i = rch.bb.len() as u32 - 1;
        while i != u32::MAX {
            if rg.sequences[rch.bb[i as usize].hid as usize].deleted {
                cc += 1;
            } else {
                break;
            }
            i = rch.bb[i as usize].pidx;
        }
        let mut i = c0;
        while i != u32::MAX {
            if rg.sequences[rch.bb[i as usize].hid as usize].deleted {
                cc += 1;
            } else {
                break;
            }
            i = rch.bb[i as usize].aidx;
        }
        if nc + cc == m {
            return 1;
        }
        return 2;
    }
    0
}

/// `direct_gchain` (`inter.cpp` 16151).
#[allow(clippy::too_many_arguments)]
fn direct_gchain(
    rg: &Asg,
    index: &UlIndex,
    read_lengths: &[u64],
    rch: &mut UlReadVec,
    diff_ec_ul: f64,
    linear_only: bool,
) -> u32 {
    let mut u_cl: Vec<UlOv> = Vec::new();
    gl_rg2ug_gen(rch, &mut u_cl, index, read_lengths, true, 2);
    if u_cl.is_empty() {
        return 0;
    }
    // generate linear chains; `raw_chn` collects the DP scratch the anchor
    // walk consumes
    let mut raw_chn: Vec<UlOv> = Vec::new();
    gen_linear_chains(
        &mut u_cl,
        &mut raw_chn,
        index,
        G_CHAIN_BW,
        diff_ec_ul,
        i64::from(rch.rlen),
    );
    if u_cl.is_empty() {
        return 0;
    }
    let mut lc = dump_linear_chain(&index.unitig.unitigs, &u_cl, i64::from(rch.rlen));
    if lc.is_empty() {
        return 0;
    }

    let (max_idx, chain_idx) = gl_chain_graph(
        index,
        &index.unitig.graph,
        &mut lc,
        i64::from(rch.rlen),
        G_CHAIN_BW,
        if linear_only { -1.0 } else { diff_ec_ul },
        N_GCHAIN_RATE,
        UG_SKIP_GRAPH_N,
        UG_ITER_N,
        UG_DIS_N,
    );
    if max_idx >= 0 {
        let mut dst_done: Vec<u64> = Vec::new();
        let mut out: Vec<usize> = Vec::new();
        let mut path: Vec<MgPathv> = Vec::new();
        let mut gchains: Vec<MgLchain> = Vec::new();
        let mut raw_idx: Vec<UlOv> = Vec::new();
        let selected = select_max_gchain(
            &index.unitig.graph,
            &chain_idx,
            &lc,
            &mut raw_idx,
            &mut dst_done,
            &mut out,
            &mut path,
            &mut gchains,
        );
        if selected > 0 {
            if !linear_only
                || (gchains[0].qe - gchains[0].qs) as f64
                    >= rch.rlen as f64 * P_CHAIN_COV
            {
                update_ul_vec_t(
                    rg,
                    index,
                    read_lengths,
                    &u_cl,
                    &raw_chn,
                    &gchains,
                    rch,
                    UG_SKIP_GRAPH_N,
                    UG_ITER_N,
                    UG_DIS_N,
                );
            }
            return u32::from(rch.dd == 1);
        }
    }
    0
}

/// `ul_refine_alignment` (`inter.cpp` 21892): rebuilds every UL read's chain
/// against the current read graph.
pub fn ul_refine_alignment(
    options: &mut AssemblyOptions,
    reads: &ReadStore,
    sources: &mut [Vec<PafHit>],
    read_to_unitig: &ReadToUnitig,
    store: &mut AllUl,
    sg: &mut Asg,
) -> bool {
    eprintln!("[M::ul_refine_alignment::] ==> UL refinement...");
    options.update_coverage(options.homozygous_coverage);
    let index = crate::assembly::ul_index::gen_ul_index_t(
        options,
        reads,
        sources,
        sg,
        false,
        false,
    );
    let start_time = logging::realtime_seconds();
    let mut tried = 0u64;
    let mut done = 0u64;
    for read_id in 0..store.reads.len() {
        let ff = refine_rid_chain(sg, &store.reads[read_id]);
        if ff == 1 {
            continue;
        }
        tried += 1;
        let linear_only = ff == 2;
        done += u64::from(direct_gchain(
            sg,
            &index,
            &reads.read_lengths,
            &mut store.reads[read_id],
            options.ultra_long_error_rate,
            linear_only,
        ) != 0);
    }
    eprintln!(
        "[M::work_ul_gchains::] # try:{}, # done:{}",
        tried, done
    );
    let _ = start_time;
    if done == 0 {
        return false;
    }

    store.ridx = Default::default();
    gen_ul_vec_rid_t(store, reads.total_reads as usize);
    let source_options = UlSourceOptions {
        maximum_hang: options.maximum_hang_length,
        maximum_hang_rate: options.maximum_hang_rate,
        minimum_overlap: options.minimum_overlap_length,
    };
    update_ovlp_src(store, sources, reads, &source_options);
    update_ovlp_src_bl(sources, reads);
    true
}

// silence the unused warnings for helpers kept for the later renewal chain
#[allow(unused_imports)]
use crate::assembly::hit::ma_hit_to_arc as _ma_hit_to_arc;

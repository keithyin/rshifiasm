//! The contained-read consensus cluster of `clean_contain_g`
//! (`inter.cpp` 17640–19990): for every UL read whose primary chain contains
//! HiFi reads, re-chain the covered regions across the read graph and rebuild
//! the block chain without them, then drop the covered contained reads from
//! the graph.
//!
//! Skipped (dead in the reference): `dedup_contain_g`/
//! `work_rg_contain_dedup` and its `gen_dup_path`/`collect_pp_ovlps`/
//! `renew_consensus_chain`/`gen_cns_chain_linear_hard` helpers — every call
//! site in the reference is commented out.

use crate::assembly::asg::{Asg, AsgArc};
use crate::assembly::hit::{ma_hit_to_arc, ReadToUnitig};
use crate::assembly::ul::get_offset_adjust;
use crate::assembly::ul_raw_chain::find_ul_ov_max;
use crate::assembly::ul::UcBlock;
use crate::assembly::ul_store::{gen_ul_vec_rid_t, AllUl, UlReadVec};
use crate::assembly::ul_src::{update_ovlp_src, update_ovlp_src_bl, UlSourceOptions};
use crate::command_line::options::AssemblyOptions;
use crate::logging;
use crate::overlap::paf::PafHit;
use crate::reads::ReadStore;

/// `is_contain_r` (`gfa_ut.h` 6): the read maps to a unitig *read* slot
/// (`index != -1` and the unitig flag is clear).
fn is_contain_r(read_to_unitig: &ReadToUnitig, read: usize) -> bool {
    let (unitig, is_unitig) = read_to_unitig.get(read as u32);
    read < read_to_unitig.len() && unitig != u32::MAX && is_unitig == 0
}

/// `G_CHAIN_BW` (`inter.h`).
const G_CHAIN_BW: i64 = 16;
/// `UG_SKIP_N` (`inter.h`).
const UG_SKIP_N: i64 = 100;
/// `UG_ITER_N` (`inter.h`).
const UG_ITER_N: i64 = 5000;
/// `UG_DIS_N` (`inter.h`).
const UG_DIS_N: i64 = 50000;
/// `G_CHAIN_INDEL` (`inter.h`).
const G_CHAIN_INDEL: i64 = 128;

fn is_rg_connect(graph: &Asg, v: u32, w: u32) -> bool {
    let (start, count) = graph.arc_range(v);
    graph.arcs[start..start + count]
        .iter()
        .any(|arc| !arc.deleted && arc.v == w)
}

/// `refine_contain_consensus_chain` (`inter.cpp` 17640): finds the primary
/// chains of `rch` that contain reads. On return `heads[0..result]` holds the
/// block indices of those chains' tails.
pub(crate) fn refine_contain_consensus_chain(
    graph: &Asg,
    rch: &UlReadVec,
    read_to_unitig: &ReadToUnitig,
    heads: &mut Vec<u64>,
) -> u32 {
    heads.clear();
    if rch.bb.len() == 1 && rch.bb[0].base {
        return 0; // no alignment
    }
    if rch.bb.is_empty() {
        return 0; // no alignment
    }
    let block_count = rch.bb.len();
    heads.resize(block_count, u64::MAX);

    // The chain walk continues through `pidx` while the anchor link
    // (`aidx`) stays inside the connected component of the read graph.
    let walk_next = |i: usize| -> usize {
        let next = rch.bb[i].pidx;
        if next == u32::MAX {
            u32::MAX as usize
        } else {
            let anchor = rch.bb[next as usize].aidx as usize;
            let connected = is_rg_connect(
                graph,
                (rch.bb[next as usize].hid << 1) | u32::from(rch.bb[next as usize].rev),
                (rch.bb[anchor].hid << 1) | u32::from(rch.bb[anchor].rev),
            );
            if connected {
                next as usize
            } else {
                u32::MAX as usize
            }
        }
    };

    for k in (0..block_count).rev() {
        if heads[k] != u64::MAX {
            continue;
        }
        let mut is_end = false;
        if rch.bb[k].aidx == u32::MAX && rch.bb[k].pidx != u32::MAX {
            is_end = true;
        } else if rch.bb[k].aidx != u32::MAX && rch.bb[k].pidx != u32::MAX {
            let anchor = rch.bb[k].aidx as usize;
            if !is_rg_connect(
                graph,
                (rch.bb[k].hid << 1) | u32::from(rch.bb[k].rev),
                (rch.bb[anchor].hid << 1) | u32::from(rch.bb[anchor].rev),
            ) {
                is_end = true;
            }
        }
        if !is_end {
            continue;
        }
        let mut i = k;
        let mut m = k;
        let mut occ = 0u64;
        while i != u32::MAX as usize {
            if graph.sequences[rch.bb[i].hid as usize].deleted {
                if occ > 1 {
                    heads[m] = occ;
                }
                heads[i] = 0;
                m = u32::MAX as usize;
                occ = 0;
            } else {
                if m == u32::MAX as usize {
                    m = i;
                }
                heads[i] = 0;
                occ += 1;
            }
            i = walk_next(i);
        }
        if occ > 1 {
            heads[m] = occ;
        }
    }

    let mut result = 0u32;
    for k in 0..block_count {
        if heads[k] == u64::MAX || heads[k] == 0 {
            continue;
        }
        let mut i = k;
        let mut occ = 0u64;
        let mut nc = 0u64;
        while i != u32::MAX as usize {
            if graph.sequences[rch.bb[i].hid as usize].deleted {
                break;
            }
            occ += 1;
            if is_contain_r(read_to_unitig, rch.bb[i].hid as usize) {
                nc += 1;
            }
            i = walk_next(i);
        }
        debug_assert_eq!(occ, heads[k]);
        if nc > 0 {
            heads[result as usize] = k as u64;
            result += 1;
        }
    }
    result
}

/// `extract_nccov` (`inter.cpp` 17784): projects the covered stretch of one
/// HiFi overlap onto the UL block's target-read coordinates. `min_ovlp` and
/// `max_hang` are dead parameters of the reference and stay unused here.
#[allow(clippy::too_many_arguments)]
fn extract_nccov(hit: &PafHit, rovlp: &UcBlock, graph: &Asg, adjust_rev: bool) -> Option<(u32, u32, u32, u32, u32, bool)> {
    let tn = hit.target_identifier as usize;
    if graph.sequences[tn].deleted {
        return None;
    }
    let os = i64::from(rovlp.ts).max(i64::from(hit.query_start));
    let oe = i64::from(rovlp.te).min(i64::from(hit.query_end));
    if oe <= os {
        return None;
    }

    // `[os, oe) -> rovlp->t*` (the block's coordinates on its own read)
    let mut s_shift = get_offset_adjust(
        (os - i64::from(rovlp.ts)) as u32,
        rovlp.te - rovlp.ts,
        rovlp.qe - rovlp.qs,
    ) as i64;
    let mut e_shift = get_offset_adjust(
        (rovlp.te as i64 - oe) as u32,
        rovlp.te - rovlp.ts,
        rovlp.qe - rovlp.qs,
    ) as i64;
    if rovlp.rev {
        std::mem::swap(&mut s_shift, &mut e_shift);
    }
    let qs = i64::from(rovlp.qs) + s_shift;
    let qe = i64::from(rovlp.qe) - e_shift;
    if qs >= qe {
        return None;
    }

    // `[os, oe) -> in->q*` (the overlap's coordinates on the target read)
    let query_span = hit.query_end - hit.query_start;
    let target_span = hit.target_end - hit.target_start;
    let mut s_shift = get_offset_adjust(
        (os - i64::from(hit.query_start)) as u32,
        query_span,
        target_span,
    ) as i64;
    let mut e_shift = get_offset_adjust(
        (i64::from(hit.query_end) - oe) as u32,
        query_span,
        target_span,
    ) as i64;
    if hit.reverse {
        std::mem::swap(&mut s_shift, &mut e_shift);
    }
    let ts = i64::from(hit.target_start) + s_shift;
    let te = i64::from(hit.target_end) - e_shift;
    if ts >= te {
        return None;
    }

    let rev = rovlp.rev != hit.reverse;
    let mut ts = ts;
    let mut te = te;
    if adjust_rev && rev {
        // for linear chaining; both expressions read the original values
        let len = i64::from(graph.sequences[tn].length);
        let (orig_ts, orig_te) = (ts, te);
        ts = len - orig_te;
        te = len - orig_ts;
    }
    Some((qs as u32, qe as u32, tn as u32, ts as u32, te as u32, rev))
}

/// `collect_nc_ovlps` (`inter.cpp` 17819): gathers the read-graph neighbours
/// of every contained block of the chain plus the chain blocks themselves.
pub(crate) fn collect_nc_ovlps(
    rch: &UlReadVec,
    rch_i: usize,
    res: &mut Vec<crate::assembly::ul::UlOv>,
    sources: &[Vec<PafHit>],
    graph: &Asg,
    read_to_unitig: &ReadToUnitig,
) -> (u64, u64) {
    use crate::assembly::ul::UlOv;
    res.clear();
    let mut mqs = u64::MAX;
    let mut mqe = u64::MAX;
    let mut occ = 0u64;
    let mut nc = 0u64;

    let walk_next = |i: usize| -> usize {
        let next = rch.bb[i].pidx;
        if next == u32::MAX {
            u32::MAX as usize
        } else {
            let anchor = rch.bb[next as usize].aidx as usize;
            let connected = is_rg_connect(
                graph,
                (rch.bb[next as usize].hid << 1) | u32::from(rch.bb[next as usize].rev),
                (rch.bb[anchor].hid << 1) | u32::from(rch.bb[anchor].rev),
            );
            if connected {
                next as usize
            } else {
                u32::MAX as usize
            }
        }
    };

    let mut i = rch_i;
    while i != u32::MAX as usize {
        if graph.sequences[rch.bb[i].hid as usize].deleted {
            break;
        }
        if is_contain_r(read_to_unitig, rch.bb[i].hid as usize) {
            for hit in &sources[rch.bb[i].hid as usize] {
                if graph.sequences[hit.target_identifier as usize].deleted {
                    continue; // tn must exist
                }
                let Some((qs, qe, tn, ts, te, rev)) =
                    extract_nccov(hit, &rch.bb[i], graph, true)
                else {
                    continue;
                };
                let tn = (tn << 1) | u32::from(rev);
                res.push(UlOv {
                    qn: i as u32,
                    qs,
                    qe,
                    tn,
                    ts,
                    te,
                    sec: 0,
                    el: false,
                    rev,
                });
            }
            nc += 1;
        }
        // push itself into the chain
        let mut ts = rch.bb[i].ts;
        let mut te = rch.bb[i].te;
        let rev = rch.bb[i].rev;
        if rev {
            let len = i64::from(graph.sequences[rch.bb[i].hid as usize].length);
            let (orig_ts, orig_te) = (ts, te);
            ts = (len - i64::from(orig_te)) as u32;
            te = (len - i64::from(orig_ts)) as u32;
        }
        res.push(UlOv {
            qn: i as u32,
            qs: rch.bb[i].qs,
            qe: rch.bb[i].qe,
            tn: (rch.bb[i].hid << 1) | u32::from(rev),
            ts,
            te,
            sec: 0,
            el: true,
            rev,
        });
        if mqs == u64::MAX || mqs > u64::from(rch.bb[i].qs) {
            mqs = u64::from(rch.bb[i].qs);
        }
        if mqe == u64::MAX || mqe < u64::from(rch.bb[i].qe) {
            mqe = u64::from(rch.bb[i].qe);
        }
        occ += 1;

        i = walk_next(i);
    }
    debug_assert!(occ > 1 && nc > 0);
    (mqs, mqe)
}

/// `comput_rlinear_sc` (`inter.cpp` 17867).
fn comput_rlinear_sc(
    li: &crate::assembly::ul::UlOv,
    lj: &crate::assembly::ul::UlOv,
    jidx: usize,
    bq: &[i32],
    bt: &[i32],
    diff_ec_ul: f64,
    bw: i64,
) -> i64 {
    use crate::assembly::ul::UlOv;
    // li is the suffix of lj; sorted by qe, so li->qe >= lj->qe
    if li.te < lj.te {
        return i32::MIN as i64;
    }
    let mut oe = li.qe.min(lj.qe);
    let mut os = li.qs.max(lj.qs);
    if oe <= os {
        return i32::MIN as i64;
    }
    oe = li.te.min(lj.te);
    os = li.ts.max(lj.ts);
    if oe <= os {
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
    let mut dq = i64::from(li.qe) - i64::from(bq[jidx]);
    let mut dt = i64::from(li.te) - i64::from(bt[jidx]);
    if dq < i64::from(li.qe - li.qs) {
        dq = i64::from(li.qe - li.qs);
    }
    if dt < i64::from(li.te - li.ts) {
        dt = i64::from(li.te - li.ts);
    }
    dq.min(dt)
}

/// `linear_rchain_dp_adv` (`inter.cpp` 17890): the backwards linear-chain DP
/// over one target read's entries.
#[allow(clippy::too_many_arguments)]
fn linear_rchain_dp_adv(
    ch: &mut [crate::assembly::ul::UlOv],
    bw: i64,
    diff_ec_ul: f64,
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
    graph: &Asg,
) -> usize {
    use crate::assembly::ul::UlOv;
    let ch_n = ch.len();
    if ch_n == 0 {
        return 0;
    }
    let mut sv = vec![UlOv::default(); ch_n];
    let mut f = vec![0i32; ch_n];
    let mut p = vec![0i64; ch_n];
    let mut t: Vec<u64> = vec![0; ch_n];
    let mut bq = vec![0i32; ch_n];
    let mut bt = vec![0i32; ch_n];

    // `radix_sort_ul_ov_srt_qe`, then per-`qe` group by `qs`.
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

    let mut st: i64 = 0;
    let mut max_ii: i64 = -1;
    for i in 0..ch_n {
        let li = ch[i];
        let csc = i64::from((li.qe - li.qs).min(li.te - li.ts));
        let mut mm_sc: i64 = -1;
        let mut mm_idx: i64 = -1;
        let mut n_skip = 0i64;
        st = if (i as i64) < max_iter { 0 } else { i as i64 - max_iter };
        let mut jj = i as i64 - 1;
        let mut end_j: i64;
        while jj >= st {
            let lj = ch[jj as usize];
            if lj.qe <= li.qs {
                break;
            }
            let sc = comput_rlinear_sc(&li, &lj, jj as usize, &bq, &bt, diff_ec_ul, bw);
            if sc == i32::MIN as i64 {
                jj -= 1;
                continue;
            }
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
                if max < i64::from(f[j2 as usize]) {
                    max = i64::from(f[j2 as usize]);
                    max_ii = j2;
                }
                j2 -= 1;
            }
        }

        if max_ii >= 0 && max_ii < end_j {
            let lj = ch[max_ii as usize];
            if lj.qe > li.qs {
                let sc = comput_rlinear_sc(&li, &lj, max_ii as usize, &bq, &bt, diff_ec_ul, bw);
                if sc != i32::MIN as i64 && sc > mm_sc {
                    mm_sc = sc;
                    mm_idx = max_ii;
                }
            }
        }

        if mm_idx < 0 {
            mm_sc = csc;
        }
        f[i] = mm_sc as i32;
        p[i] = mm_idx;
        bq[i] = li.qs as i32;
        bt[i] = li.ts as i32;
        if mm_idx >= 0 {
            if bq[i] > bq[mm_idx as usize] {
                bq[i] = bq[mm_idx as usize];
            }
            if bt[i] > bt[mm_idx as usize] {
                bt[i] = bt[mm_idx as usize];
            }
        }

        if max_ii < 0
            || (i64::from(ch[i].qe) <= max_dis + i64::from(ch[max_ii as usize].qe)
                && i64::from(f[max_ii as usize]) < mm_sc)
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

    // `t[i] = f[i]<<32 | i`, then the descending score order of
    // `radix_sort_gfa64i` walked from the back.
    for i in 0..ch_n {
        t[i] = ((f[i] as u32 as u64) << 32) + i as u64;
        bt[i] = 0;
    }
    t.sort(); // ascending, like the stable radix sort; ties impossible

    let mut n_u = 0usize;
    for z in (0..ch_n).rev() {
        let k = (t[z] & 0xFFFF_FFFF) as usize;
        if bt[k] != 0 {
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
            bt[i] = 1;
            if p[i] < 0 {
                break;
            }
            i = p[i] as usize;
        }
        let len = i64::from(graph.sequences[ch[n_u].tn as usize].length);
        let (its, ite) = adjust_rev_tse(&ch[n_u], len);
        ch[n_u].ts = its as u32;
        ch[n_u].te = ite as u32;
        ch[n_u].sec = if sc > 0x3FFF_FFFF { 0x3FFF_FFFF } else { sc as u32 };
        ch[n_u].qn = 0;
        n_u += 1;
    }
    n_u
}

/// `adjust_rev_tse` (`inter.cpp` 13379).
fn adjust_rev_tse(x: &crate::assembly::ul::UlOv, tlen: i64) -> (i64, i64) {
    let mut ts = i64::from(x.ts);
    let mut te = i64::from(x.te);
    if x.rev {
        ts = tlen - i64::from(x.te);
        te = tlen - i64::from(x.ts);
    }
    (ts, te)
}

/// `gen_linear_rchains` (`inter.cpp` 17960): chains the collected overlaps
/// per target read. The reference's `buf` only carries DP scratch space and
/// is not part of the result.
pub(crate) fn gen_linear_rchains(
    res: &mut Vec<crate::assembly::ul::UlOv>,
    graph: &Asg,
    bw: i64,
    diff_ec_ul: f64,
    qlen: i64,
) {
    use crate::assembly::ul::UlOv;
    let _ = qlen;
    crate::container::ksort::radix_sort_by_key_32(res, &|entry: &UlOv| entry.tn);
    let n = res.len();
    let mut l = 0usize;
    let mut m = 0usize;
    for k in 1..=n {
        if k == n || res[k].tn != res[l].tn {
            for z in l..k {
                res[z].tn >>= 1;
            }
            let an = l
                + linear_rchain_dp_adv(
                    &mut res[l..k],
                    bw,
                    diff_ec_ul,
                    UG_SKIP_N,
                    UG_ITER_N,
                    UG_DIS_N,
                    graph,
                );
            for z in l..an {
                res[m] = res[z];
                m += 1;
            }
            l = k;
        }
    }
    res.truncate(m);
}

/// `gconnect_test` (`inter.cpp` 18011).
fn gconnect_test(graph: &Asg, v: u32, w: u32, bw: i64, diff_ec_ul: f64, dq: i64) -> bool {
    let (start, count) = graph.arc_range(v);
    let mut dt: i64 = -1;
    for arc in &graph.arcs[start..start + count] {
        if arc.deleted || arc.v != w {
            continue;
        }
        // the reference reads the  field here, not 
        dt = i64::from(arc.overlap);
        break;
    }
    if dt < 0 {
        return false;
    }
    let dif = (dq - dt).abs();
    let mut mm = dq.max(dt);
    mm = (mm as f64 * diff_ec_ul) as i64;
    if mm < bw {
        mm = bw;
    }
    dif <= mm
}

/// `infer_rovlp` (`Correct.cpp` 19083) with `ug = NULL`: the lengths are the
/// **full** `R_INF` read lengths, not the coverage-trimmed graph lengths.
fn infer_rovlp_read_graph(
    li: &crate::assembly::ul::UlOv,
    lj: &crate::assembly::ul::UlOv,
    read_lengths: &[u64],
) -> i64 {
    let len = |id: u32| read_lengths[id as usize] as i64;
    let i_len = len(li.tn);
    let j_len = len(lj.tn);
    let is = i64::from(li.ts);
    let ie = i64::from(li.te);
    let irev = li.rev;
    let iqs = i64::from(li.qs);
    let iqe = i64::from(li.qe);
    let js = i64::from(lj.ts);
    let je = i64::from(lj.te);
    let jrev = lj.rev;
    let jqs = i64::from(lj.qs);
    let jqe = i64::from(lj.qe);

    let max_s = iqs.max(jqs);
    let min_e = iqe.min(jqe);
    if min_e <= max_s {
        return 0;
    }
    let mut s_shift = get_offset_adjust(
        (max_s - iqs) as u32,
        (iqe - iqs) as u32,
        (ie - is) as u32,
    ) as i64;
    let mut e_shift = get_offset_adjust(
        (iqe - min_e) as u32,
        (iqe - iqs) as u32,
        (ie - is) as u32,
    ) as i64;
    if irev {
        std::mem::swap(&mut s_shift, &mut e_shift);
    }
    let mut is = is + s_shift;
    let mut ie = ie - e_shift;

    let mut s_shift = get_offset_adjust(
        (max_s - jqs) as u32,
        (jqe - jqs) as u32,
        (je - js) as u32,
    ) as i64;
    let mut e_shift = get_offset_adjust(
        (jqe - min_e) as u32,
        (jqe - jqs) as u32,
        (je - js) as u32,
    ) as i64;
    if jrev {
        std::mem::swap(&mut s_shift, &mut e_shift);
    }
    let mut js = js + s_shift;
    let mut je = je - e_shift;

    if irev {
        let ts = i_len - ie;
        let te = i_len - is;
        is = ts;
        ie = te;
    }
    if jrev {
        let ts = j_len - je;
        let te = j_len - js;
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
    let ir = i_len - ie;
    let jr = j_len - je;
    if ir <= jr {
        ie = i_len;
        je += ir;
    } else {
        je = j_len;
        ie += jr;
    }
    (ie - is).max(je - js)
}

/// `gen_cns_chain_linear` (`inter.cpp` 18037): picks the single best
/// left-to-right chain that spans `[mqs, mqe]`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn gen_cns_chain_linear(
    a: &mut [crate::assembly::ul::UlOv],
    graph: &Asg,
    read_lengths: &[u64],
    read_to_unitig: &ReadToUnitig,
    qlen: i64,
    bw: i64,
    diff_thre: f64,
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
    mqs: u64,
    mqe: u64,
) -> usize {
    use crate::assembly::ul::UlOv;
    let a_n = a.len();
    if a_n == 0 {
        return 0;
    }
    let mut f = vec![0i32; a_n];
    let mut p = vec![0i64; a_n];
    let mut t = vec![0i64; a_n];
    let mut c_n = vec![0i32; a_n];
    let mut len = vec![0i32; a_n];

    crate::container::ksort::radix_sort_by_key_32(a, &|entry: &UlOv| entry.qe);
    let mut j = 0usize;
    for i in 1..=a_n {
        if i == a_n || a[i].qe != a[j].qe {
            if i - j > 1 {
                crate::container::ksort::radix_sort_by_key_32(&mut a[j..i], &|e: &UlOv| e.qs);
            }
            j = i;
        }
    }

    let int_min = i32::MIN as i64;
    let mut ch_sc: i64 = int_min;
    let mut ch_i: i64 = int_min;
    let mut cn_sn: i64 = int_min;
    let mut ch_ln: i64 = i32::MAX as i64;
    let mut st: i64 = 0;
    let mut max_ii: i64 = -1;
    for i in 0..a_n {
        let li = a[i];
        let mm_ovlp = crate::assembly::ul::max_ovlp(graph, ((li.tn << 1) | u32::from(li.rev)) ^ 1);
        let mut x = ((i64::from(li.qs) + mm_ovlp) as f64 * diff_thre) as i64;
        if x < bw {
            x = bw;
        }
        x += i64::from(li.qs) + mm_ovlp;
        if x > qlen + 1 {
            x = qlen + 1;
        }
        let x = i64::from(find_ul_ov_max(i as i32, a, (x + G_CHAIN_INDEL) as u32));

        let mut csc: i64 = 0;
        let mut csn: i64 = 1;
        let mut cln: i64 = 0;
        if is_contain_r(read_to_unitig, li.tn as usize) {
            csc = -1;
            csn = 0;
            cln = i64::from(graph.sequences[li.tn as usize].length);
        }
        let mut mm_sc: i64 = int_min + 1;
        let mut mm_sn: i64 = int_min + 1;
        let mut mm_ln: i64 = i32::MAX as i64;
        let mut mm_idx: i64 = -1;
        let mut n_skip = 0i64;
        let mut end_j: i64;

        let li_v = ((li.tn << 1) | u32::from(li.rev)) ^ 1;
        if x - st > max_iter {
            st = x - max_iter;
        }
        let mut jj = x;
        while jj >= st {
            let lj = a[jj as usize];
            if i64::from(lj.qe) + G_CHAIN_INDEL <= i64::from(li.qs) {
                break;
            }
            if f[jj as usize] == i32::MIN {
                jj -= 1;
                continue;
            }
            let mut sc: i64 = int_min;
            let mut sn: i64 = int_min;
            let mut ln: i64 = int_min;
            if lj.qs <= li.qs && lj.qe <= li.qe {
                // not contained
                let lj_v = ((lj.tn << 1) | u32::from(lj.rev)) ^ 1;
                let qo = infer_rovlp_read_graph(&li, &lj, read_lengths);
                if gconnect_test(graph, li_v, lj_v, bw, diff_thre, qo) {
                    sc = i64::from(f[jj as usize]) + csc;
                    sn = i64::from(c_n[jj as usize]) + csn;
                    ln = i64::from(len[jj as usize]) + cln;
                }
            }
            if sc == int_min {
                jj -= 1;
                continue;
            }
            if (sc > mm_sc)
                || (sc == mm_sc && sn > mm_sn)
                || (sc == mm_sc && sn == mm_sn && ln < mm_ln)
            {
                mm_sc = sc;
                mm_idx = jj;
                mm_sn = sn;
                mm_ln = ln;
                if n_skip > 0 {
                    n_skip -= 1;
                }
            } else if t[jj as usize] == i as i64 {
                n_skip += 1;
                if n_skip > max_skip {
                    break;
                }
            }
            if p[jj as usize] >= 0 {
                t[p[jj as usize] as usize] = i as i64;
            }
            jj -= 1;
        }

        end_j = jj;
        if max_ii < 0 || (i64::from(a[i].qe) > i64::from(a[max_ii as usize].qe) + max_dis) {
            let mut max = int_min;
            max_ii = -1;
            let mut j2 = i as i64 - 1;
            while j2 >= st && i64::from(a[i].qe) <= max_dis + i64::from(a[j2 as usize].qe) {
                if max < i64::from(f[j2 as usize]) {
                    max = i64::from(f[j2 as usize]);
                    max_ii = j2;
                }
                j2 -= 1;
            }
        }

        if max_ii >= 0 && max_ii < end_j {
            let lj = a[max_ii as usize];
            if i64::from(lj.qe) + G_CHAIN_INDEL > i64::from(li.qs)
                && lj.qs <= li.qs
                && lj.qe <= li.qe
                && f[max_ii as usize] != i32::MIN
            {
                let lj_v = ((lj.tn << 1) | u32::from(lj.rev)) ^ 1;
                let mut sc: i64 = int_min;
                let mut sn: i64 = int_min;
                let mut ln: i64 = int_min;
                let qo = infer_rovlp_read_graph(&li, &lj, read_lengths);
                if gconnect_test(graph, li_v, lj_v, bw, diff_thre, qo) {
                    sc = i64::from(f[max_ii as usize]) + csc;
                    sn = i64::from(c_n[max_ii as usize]) + csn;
                    ln = i64::from(len[max_ii as usize]) + cln;
                }
                if sc != int_min
                    && ((sc > mm_sc)
                        || (sc == mm_sc && sn > mm_sn)
                        || (sc == mm_sc && sn == mm_sn && ln < mm_ln))
                {
                    mm_sc = sc;
                    mm_idx = max_ii;
                    mm_sn = sn;
                    mm_ln = ln;
                }
            }
        }

        if mm_idx == -1 && i64::from(li.qs) > mqs as i64 + bw {
            mm_sc = int_min;
            mm_idx = int_min;
            mm_sn = int_min;
            mm_ln = i32::MAX as i64;
        }
        if mm_sc == int_min + 1 {
            mm_sc = csc;
            mm_sn = csn;
            mm_ln = cln;
        }

        f[i] = mm_sc as i32;
        p[i] = mm_idx;
        c_n[i] = mm_sn as i32;
        len[i] = mm_ln as i32;
        if max_ii < 0
            || (i64::from(a[i].qe) <= max_dis + i64::from(a[max_ii as usize].qe)
                && i64::from(f[max_ii as usize]) < mm_sc)
        {
            max_ii = i as i64;
        }

        if mm_sc != int_min && (i64::from(a[i].qe) + bw) >= mqe as i64 {
            if (mm_sc > ch_sc)
                || (mm_sc == ch_sc && mm_sn > cn_sn)
                || (mm_sc == ch_sc && mm_sn == cn_sn && mm_ln < ch_ln)
            {
                ch_sc = mm_sc;
                ch_i = i as i64;
                cn_sn = mm_sn;
                ch_ln = mm_ln;
            }
        }
    }

    if ch_i == int_min {
        return 0;
    }

    let mut path = Vec::new();
    let mut i = ch_i;
    while i >= 0 {
        path.push(i);
        i = p[i as usize];
    }
    let cl = path.len();
    for (slot, _) in path.iter().enumerate() {
        a[slot] = a[path[cl - 1 - slot] as usize];
    }
    cl
}

/// `gen_contain_consensus_chain` (`inter.cpp` 18184).
#[allow(clippy::too_many_arguments)]
fn gen_contain_consensus_chain(
    rch: &mut UlReadVec,
    rch_i: usize,
    idx: &mut Vec<crate::assembly::ul::UlOv>,
    dump: &mut Vec<crate::assembly::ul::UlOv>,
    sources: &[Vec<PafHit>],
    graph: &Asg,
    read_lengths: &[u64],
    read_to_unitig: &ReadToUnitig,
    bw: i64,
    diff_ec_ul: f64,
) {
    use crate::assembly::ul::UlOv;
    let dn = dump.len();
    idx.clear();
    let (mqs, mqe) = collect_nc_ovlps(rch, rch_i, idx, sources, graph, read_to_unitig);
    if idx.len() <= 1 || mqs == u64::MAX || mqe == u64::MAX || mqe <= mqs {
        return;
    }

    gen_linear_rchains(idx, graph, bw, diff_ec_ul, i64::from(rch.rlen));
    debug_assert!(!idx.is_empty());
    let chained = gen_cns_chain_linear(
        idx,
        graph,
        read_lengths,
        read_to_unitig,
        i64::from(rch.rlen),
        bw,
        diff_ec_ul,
        UG_SKIP_N,
        UG_ITER_N,
        UG_DIS_N,
        mqs,
        mqe,
    );
    idx.truncate(chained);
    if !idx.is_empty() {
        dump.resize(dn + idx.len(), UlOv::default());
        // mask the existing overlaps of the chain
        let mut i = rch_i;
        while i != u32::MAX as usize {
            if graph.sequences[rch.bb[i].hid as usize].deleted {
                break;
            }
            dump[i].tn = u32::MAX;
            dump[i].qn = u32::MAX;
            i = if rch.bb[i].pidx == u32::MAX {
                u32::MAX as usize
            } else {
                rch.bb[i].pidx as usize
            };
        }
        for (slot, entry) in idx.iter().enumerate() {
            let mut record = *entry;
            let tlen = i64::from(graph.sequences[record.tn as usize].length);
            let (iqs, iqe, its, ite) = extend_end_coord_ul_ov(&record, i64::from(rch.rlen), tlen);
            record.qs = iqs as u32;
            record.qe = iqe as u32;
            record.ts = its as u32;
            record.te = ite as u32;
            record.qn = if slot > 0 { (dn + slot - 1) as u32 } else { u32::MAX };
            dump[dn + slot] = record;
        }
    }
}

/// `extend_end_coord` (`inter.cpp` 13607) for a `ul_ov_t`.
pub(crate) fn extend_end_coord_ul_ov(
    z: &crate::assembly::ul::UlOv,
    qlen: i64,
    tlen: i64,
) -> (i64, i64, i64, i64) {
    let mut qs = i64::from(z.qs);
    let mut qe = i64::from(z.qe);
    let mut ts = i64::from(z.ts);
    let mut te = i64::from(z.te);
    if z.rev {
        ts = tlen - i64::from(z.te);
        te = tlen - i64::from(z.ts);
    }
    if qs <= ts {
        ts -= qs;
        qs = 0;
    } else {
        qs -= ts;
        ts = 0;
    }
    let qtail = qlen - qe;
    let ttail = tlen - te;
    if qtail <= ttail {
        qe = qlen;
        te += qtail;
    } else {
        te = tlen;
        qe += ttail;
    }
    let (r_qs, r_qe, r_ts, r_te) = (qs, qe, ts, te);
    if z.rev {
        return (r_qs, r_qe, tlen - te, tlen - ts);
    }
    (r_qs, r_qe, r_ts, r_te)
}

/// `update_consensus_chain` (`inter.cpp` 18212).
fn update_consensus_chain(
    sources: &[Vec<PafHit>],
    reads: &ReadStore,
    graph: &Asg,
    idx: &mut Vec<crate::assembly::ul::UlOv>,
    dump: &mut Vec<crate::assembly::ul::UlOv>,
    rch: &mut UlReadVec,
    idm: &mut Vec<u64>,
    read_to_unitig: &ReadToUnitig,
    maximum_hang: i32,
    maximum_hang_rate: f32,
    minimum_overlap: i32,
) -> u32 {
    use crate::assembly::ul::UlOv;
    debug_assert!(dump.len() >= rch.bb.len());
    let mut kn = 0usize;
    idm.clear();
    idm.resize(dump.len(), u64::MAX);
    for k in 0..rch.bb.len() {
        idm[k] = u64::MAX;
        if dump[k].qn == u32::MAX && dump[k].tn == u32::MAX {
            continue;
        }
        dump[k].qn = rch.bb[dump[k].qn as usize].pidx;
        idm[k] = kn as u64;
        kn += 1;
    }
    for k in rch.bb.len()..dump.len() {
        idm[k] = kn as u64;
        kn += 1;
    }

    let dn = dump.len();
    let mut kept = 0usize;
    idx.clear();
    for k in 0..dn {
        if idm[k] == u64::MAX {
            continue;
        }
        let target = idm[k] as usize;
        dump[target] = dump[k];
        kn -= 1;
        if dump[target].qn != u32::MAX {
            dump[target].qn = idm[dump[target].qn as usize] as u32;
        }
        let mut entry = dump[target];
        entry.qn = target as u32;
        idx.push(entry);
        kept += 1;
    }
    debug_assert_eq!(kn, 0);
    dump.truncate(kept);

    crate::container::ksort::radix_sort_by_key_32(idx, &|entry: &UlOv| entry.qe);
    let n = idx.len();
    let mut l = 0usize;
    for k in 1..=n {
        if k == n || idx[l].qe != idx[k].qe {
            if k - l > 1 {
                crate::container::ksort::radix_sort_by_key_32(&mut idx[l..k], &|e: &UlOv| e.qs);
            }
            l = k;
        }
    }
    for k in 0..n {
        dump[idx[k].qn as usize].tn = k as u32;
    }
    for k in 0..n {
        let dump_index = idx[k].qn as usize;
        idx[k].qn = if dump[dump_index].qn != u32::MAX {
            dump[dump[dump_index].qn as usize].tn
        } else {
            u32::MAX
        };
        idx[k].el = true;
    }

    rch.bb = Vec::with_capacity(n);
    for entry in idx.iter() {
        rch.bb.push(UcBlock {
            hid: entry.tn,
            rev: entry.rev,
            pchain: 1,
            base: false,
            el: true,
            qs: entry.qs,
            qe: entry.qe,
            te: entry.te,
            ts: entry.ts,
            pidx: entry.qn,
            pdis: u32::MAX,
            aidx: u32::MAX,
        });
    }
    let mut cc = 0u32;
    for k in 0..rch.bb.len() {
        if is_contain_r(read_to_unitig, rch.bb[k].hid as usize) {
            cc += 1;
        }
        if rch.bb[k].pidx == u32::MAX {
            continue;
        }
        let z = rch.bb[k];
        let p = &mut rch.bb[z.pidx as usize];
        debug_assert_eq!(p.aidx, u32::MAX);
        p.aidx = k as u32;
        let tt = g_adjacent_dis_mul(
            sources,
            reads,
            maximum_hang,
            maximum_hang_rate,
            minimum_overlap,
            ((z.hid << 1) | u32::from(z.rev)) ^ 1,
            ((p.hid << 1) | u32::from(p.rev)) ^ 1,
        );
        if tt >= 0 {
            rch.bb[k].pdis = tt as u32;
        }
    }
    if cc == 0 {
        1
    } else {
        0
    }
}

/// `g_adjacent_dis_mul` (`inter.cpp` 15837) with `g == NULL`: the arc length
/// between two read vertices recovered from the source overlaps.
fn g_adjacent_dis_mul(
    sources: &[Vec<PafHit>],
    reads: &ReadStore,
    maximum_hang: i32,
    maximum_hang_rate: f32,
    minimum_overlap: i32,
    v: u32,
    w: u32,
) -> i64 {
    for hit in &sources[(v >> 1) as usize] {
        let qn = hit.query_identifier;
        let tn = hit.target_identifier;
        if qn == (v >> 1) && tn == (w >> 1) {
            let (result, arc) = ma_hit_to_arc(
                hit,
                reads.read_length(qn as usize) as i32,
                reads.read_length(tn as usize) as i32,
                maximum_hang,
                maximum_hang_rate,
                minimum_overlap,
            );
            if result < 0 {
                continue;
            }
            if (arc.ul >> 32) != u64::from(v) || u64::from(arc.v) != u64::from(w) {
                continue;
            }
            return i64::from(arc.node_length());
        }
    }
    -1
}

/// The per-read body of `worker_for_contain_consensus` (`inter.cpp` 18319).
/// Returns `(tried, corrected)` for the `# try/# done` statistics.
#[allow(clippy::too_many_arguments)]
fn contain_consensus_for_read(
    store: &mut AllUl,
    read_id: usize,
    sources: &[Vec<PafHit>],
    reads: &ReadStore,
    graph: &Asg,
    read_to_unitig: &ReadToUnitig,
    diff_ec_ul: f64,
    maximum_hang: i32,
    maximum_hang_rate: f32,
    minimum_overlap: i32,
    scratch_heads: &mut Vec<u64>,
    idx: &mut Vec<crate::assembly::ul::UlOv>,
    dump: &mut Vec<crate::assembly::ul::UlOv>,
    idm: &mut Vec<u64>,
) -> (u32, u32) {
    let chain_count = refine_contain_consensus_chain(
        graph,
        &store.reads[read_id],
        read_to_unitig,
        scratch_heads,
    );
    if chain_count == 0 {
        return (0, 0);
    }

    let tried = 1u32;
    dump.clear();
    dump.resize(store.reads[read_id].bb.len(), crate::assembly::ul::UlOv::default());
    for (k, entry) in dump.iter_mut().enumerate() {
        let block = &store.reads[read_id].bb[k];
        *entry = crate::assembly::ul::UlOv {
            qn: k as u32,
            qs: block.qs,
            qe: block.qe,
            tn: block.hid,
            ts: block.ts,
            te: block.te,
            sec: 0,
            el: true,
            rev: block.rev,
        };
    }

    for k in 0..chain_count as usize {
        gen_contain_consensus_chain(
            &mut store.reads[read_id],
            scratch_heads[k] as usize,
            idx,
            dump,
            sources,
            graph,
            &reads.read_lengths,
            read_to_unitig,
            G_CHAIN_BW,
            diff_ec_ul,
        );
    }

    let corrected = update_consensus_chain(
        sources,
        reads,
        graph,
        idx,
        dump,
        &mut store.reads[read_id],
        idm,
        read_to_unitig,
        maximum_hang,
        maximum_hang_rate,
        minimum_overlap,
    );
    (tried, corrected)
}

/// `clean_contained_chg` (`inter.cpp` 19937) plus the deletion pass of
/// The probe counterpart of the C++ `trace_read` helper: replays one read's
/// consensus body stage by stage and dumps the same format.
#[allow(clippy::too_many_arguments)]
pub fn trace_read_for_probe(
    store: &AllUl,
    read_id: usize,
    sources: &[Vec<PafHit>],
    graph: &Asg,
    read_lengths: &[u64],
    read_to_unitig: &ReadToUnitig,
    diff_ec_ul: f64,
    tag: &str,
) {
    use crate::assembly::ul::UlOv;
    let rch = &store.reads[read_id];
    let mut heads = Vec::new();
    let chain_count = refine_contain_consensus_chain(graph, rch, read_to_unitig, &mut heads);
    let mut out = String::new();
    out.push_str(&format!("ff {} heads {}:", chain_count, chain_count));
    for k in 0..chain_count as usize {
        out.push_str(&format!(" {}", heads[k]));
    }
    out.push('\n');

    if chain_count > 0 {
        for k in 0..chain_count as usize {
            let rch_i = heads[k] as usize;
            let mut idx: Vec<UlOv> = Vec::new();
            let (mqs, mqe) = collect_nc_ovlps(rch, rch_i, &mut idx, sources, graph, read_to_unitig);
            out.push_str(&format!(
                "chain {} collected {} mqs {} mqe {}:",
                rch_i, idx.len(), mqs, mqe
            ));
            for e in &idx {
                out.push_str(&format!(
                    " {}:{}:{}:{}:{}:{}:{}:{}",
                    e.qn, e.qs, e.qe, e.tn, e.ts, e.te, u8::from(e.el), u8::from(e.rev)
                ));
            }
            out.push('\n');
            if idx.len() <= 1 || mqs == u64::MAX || mqe == u64::MAX || mqe <= mqs {
                continue;
            }
            gen_linear_rchains(&mut idx, graph, G_CHAIN_BW, diff_ec_ul, i64::from(rch.rlen));
            out.push_str(&format!("chain {} linear {}:", rch_i, idx.len()));
            for e in &idx {
                out.push_str(&format!(
                    " {}:{}:{}:{}:{}:{}:{}:{}",
                    e.qn, e.qs, e.qe, e.tn, e.ts, e.te, u8::from(e.el), u8::from(e.rev)
                ));
            }
            out.push('\n');
            let cn = gen_cns_chain_linear(
                &mut idx,
                graph,
                read_lengths,
                read_to_unitig,
                i64::from(rch.rlen),
                G_CHAIN_BW,
                diff_ec_ul,
                UG_SKIP_N,
                UG_ITER_N,
                UG_DIS_N,
                mqs,
                mqe,
            );
            out.push_str(&format!("chain {} cns {}:", rch_i, cn));
            for e in idx.iter().take(cn) {
                out.push_str(&format!(
                    " {}:{}:{}:{}:{}:{}:{}:{}",
                    e.qn, e.qs, e.qe, e.tn, e.ts, e.te, u8::from(e.el), u8::from(e.rev)
                ));
            }
            out.push('\n');
        }
    }
    std::fs::write(tag, out).unwrap();
}

/// `clean_contain_g0` (`inter.cpp` 21724): deletes the contained reads whose
/// source overlaps carry no UL support.
fn clean_contain_g0(graph: &mut Asg, sources: &[Vec<PafHit>], read_to_unitig: &ReadToUnitig) -> u32 {
    let read_count = graph.sequences.len();
    let mut marked = vec![false; read_count];
    for read in 0..read_count {
        if graph.sequences[read].deleted {
            continue;
        }
        if !is_contain_r(read_to_unitig, read) {
            continue;
        }
        if sources[read].iter().any(|hit| hit.block_length != 0) {
            continue;
        }
        marked[read] = true;
    }
    let mut cnt = 0u32;
    for (read, mark) in marked.iter().enumerate() {
        if *mark {
            graph.delete_sequence(read as u32);
            cnt += 1;
        }
    }
    if cnt > 0 {
        graph.cleanup();
    }
    eprintln!("[M::clean_contain_g0::] # discard cread::{}", cnt);
    cnt
}

/// The structured outcome of [`clean_contain_g`]: the counters behind the
/// reference's log lines.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CleanContainOutcome {
    /// `clean_contain_g`'s return value.
    pub changed: bool,
    /// `work_ul_gchains_consensus`'s `# try`.
    pub tried: u32,
    /// `work_ul_gchains_consensus`'s `# done`.
    pub corrected: u32,
    /// `clean_contain_g0`'s `# discard cread`.
    pub discarded: u32,
}

/// `clean_contain_g` (`inter.cpp` 21857): the whole contained-read cleanup.
///
/// `push_trans` mirrors the reference parameter: `gen_init_sg` calls this with
/// `0`, the `ul_clean_gfa` rounds with `1`.
pub fn clean_contain_g(
    options: &mut AssemblyOptions,
    reads: &ReadStore,
    sources: &mut [Vec<PafHit>],
    read_to_unitig: &ReadToUnitig,
    store: &mut AllUl,
    graph: &mut Asg,
    push_trans: bool,
) -> CleanContainOutcome {
    // `init_aux_table(); ha_opt_update_cov(&asm_opt, asm_opt.hom_cov);`
    options.update_coverage(options.homozygous_coverage);
    if push_trans {
        repush_contain_trans_archs(options, sources, read_to_unitig, graph);
    }

    let start_time = logging::realtime_seconds();
    let mut tried = 0u32;
    let mut corrected = 0u32;
    let diff_ec_ul = options.ultra_long_error_rate;
    let maximum_hang = options.maximum_hang_length;
    let maximum_hang_rate = options.maximum_hang_rate;
    let minimum_overlap = options.minimum_overlap_length;
    let mut scratch_heads = Vec::new();
    let mut idx = Vec::new();
    let mut dump = Vec::new();
    let mut idm = Vec::new();
    for read_id in 0..store.reads.len() {
        let (t, c) = contain_consensus_for_read(
            store,
            read_id,
            sources,
            reads,
            graph,
            read_to_unitig,
            diff_ec_ul,
            maximum_hang,
            maximum_hang_rate,
            minimum_overlap,
            &mut scratch_heads,
            &mut idx,
            &mut dump,
            &mut idm,
        );
        tried += t;
        corrected += c;
    }
    eprintln!(
        "[M::work_ul_gchains_consensus::] # try:{}, # done:{}",
        tried, corrected
    );
    let _ = start_time;

    let mut outcome = CleanContainOutcome {
        changed: false,
        tried,
        corrected,
        discarded: 0,
    };
    if tried > 0 {
        // `free(UL_INF.ridx...)` + `gen_ul_vec_rid_t(&UL_INF, &R_INF, NULL)`
        store.ridx = Default::default();
        gen_ul_vec_rid_t(store, reads.total_reads as usize);
        let source_options = UlSourceOptions {
            maximum_hang,
            maximum_hang_rate,
            minimum_overlap,
        };
        update_ovlp_src(store, sources, reads, &source_options);
        update_ovlp_src_bl(sources, reads);
        outcome.discarded = clean_contain_g0(graph, sources, read_to_unitig);
        outcome.changed = true;
    }
    if push_trans {
        crate::assembly::ul_graph::asg_arc_del_trans_ul(graph, options.gap_fuzz);
    }
    outcome
}

/// `repush_contain_trans_archs` (`inter.cpp` 21844): `label_contained_chg`
/// (`inter.cpp` 19951) followed by `asg_arc_push_contain_trans`
/// (`inter.cpp` 21742).
fn repush_contain_trans_archs(
    options: &AssemblyOptions,
    sources: &mut [Vec<PafHit>],
    read_to_unitig: &ReadToUnitig,
    graph: &mut Asg,
) {
    let vertex_count = graph.sequences.len() * 2;
    let mut mark = vec![false; vertex_count];
    for v in 0..vertex_count {
        if graph.sequences[v >> 1].deleted {
            continue;
        }
        if is_contain_r(read_to_unitig, v >> 1) {
            mark[v] = true;
            continue;
        }
        let (start, count) = graph.arc_range(v as u32);
        if graph.arcs[start..start + count]
            .iter()
            .any(|arc| !arc.deleted && is_contain_r(read_to_unitig, (arc.v >> 1) as usize))
        {
            mark[v] = true;
        }
    }
    asg_arc_push_contain_trans(
        graph,
        sources,
        i64::from(options.minimum_overlap_length),
        i64::from(options.maximum_hang_length),
        f64::from(options.maximum_hang_rate),
        i64::from(options.gap_fuzz),
        read_to_unitig,
        &mark,
    );
}

/// `asg_arc_push_contain_trans` (`inter.cpp` 21742).
///
/// The reference juggles raw array windows here: newly pushed arcs whose
/// target is contained stay at the tail (they seed the next round), the others
/// are swapped to the absolute front position `nc`. No arc is ever lost — the
/// final `asg_cleanup` rebuilds the graph from the arc multiset — but the
/// scan windows alias the swapped slots, so the exact absolute-index semantics
/// are reproduced here.
#[allow(clippy::too_many_arguments)]
fn asg_arc_push_contain_trans(
    graph: &mut Asg,
    overlaps: &[Vec<PafHit>],
    min_ovlp: i64,
    max_hang: i64,
    max_hang_rate: f64,
    fuzz: i64,
    read_to_unitig: &ReadToUnitig,
    mark: &[bool],
) {
    let n_vtx = (graph.sequences.len() * 2) as u32;
    let mut dis = vec![u32::MAX; n_vtx as usize];
    let mut idx = vec![0u32; n_vtx as usize];
    let mut cnt = 0u32;

    for v in 0..n_vtx {
        if graph.sequences[(v >> 1) as usize].deleted || !mark[v as usize] {
            continue;
        }
        let (avi, nv_count) = graph.arc_range(v);
        let mut avi = avi as u32;
        let mut nv = nv_count as u32;
        if nv == 0 {
            continue; // no hits
        }
        let has_contained_target = graph.arcs[avi as usize..avi as usize + nv as usize]
            .iter()
            .any(|arc| !arc.deleted && is_contain_r(read_to_unitig, (arc.v >> 1) as usize));
        if !has_contained_target {
            continue;
        }

        for (i, hit) in overlaps[(v >> 1) as usize].iter().enumerate() {
            if hit.deleted {
                continue;
            }
            let Some(p) = hit_to_arc(graph, hit, max_hang, max_hang_rate, min_ovlp) else {
                continue;
            };
            if (p.ul >> 32) != u64::from(v) {
                continue;
            }
            if graph.sequences[(p.v >> 1) as usize].deleted {
                continue;
            }
            if mark[(p.v ^ 1) as usize] {
                dis[p.v as usize] = p.node_length().wrapping_add(fuzz as u32);
                idx[p.v as usize] = i as u32;
            }
        }

        // `nv`/`avi` drive the scan window; new arcs whose target is contained
        // accumulate at the tail and seed the next round.
        while nv > 0 {
            let mut nc = 0u32;
            let mut cc = 0u32;
            let mut i = 0u32;
            while i < nv {
                let slot = (avi + i) as usize;
                if graph.arcs[slot as usize].deleted {
                    i += 1;
                    continue;
                }
                let w = graph.arcs[slot as usize].v;
                if !is_contain_r(read_to_unitig, (w >> 1) as usize) {
                    i += 1;
                    continue; // new arcs must be bridged by contained reads
                }
                debug_assert!(!graph.sequences[(w >> 1) as usize].deleted);
                let (awi, nw) = graph.arc_range(w);
                for k in 0..nw {
                    let aw = graph.arcs[awi + k];
                    if aw.deleted {
                        continue;
                    }
                    let target = aw.v as usize;
                    if dis[target] == u32::MAX || dis[target] == 0 {
                        continue;
                    }
                    debug_assert!(!graph.sequences[(aw.v >> 1) as usize].deleted);
                    let av_len = graph.arcs[slot as usize].node_length() as u64;
                    if av_len + (aw.node_length() as u64) <= u64::from(dis[target]) {
                        dis[target] = u32::MAX;
                        let hit = &overlaps[(v >> 1) as usize][idx[target] as usize];
                        debug_assert!(!hit.deleted);
                        let mut p =
                            hit_to_arc(graph, hit, max_hang, max_hang_rate, min_ovlp)
                                .expect("the recorded hit must form an arc");
                        debug_assert_eq!(p.ul >> 32, u64::from(v));
                        debug_assert_eq!(p.v, aw.v);
                        cnt += 1;
                        p.occurrence = 0;
                        graph.push_arc(p);
                        // refresh windows after the potential reallocation
                        if is_contain_r(read_to_unitig, (p.v >> 1) as usize) {
                            cc += 1;
                        } else {
                            // `if(g->n_arc != nc+1)` with the absolute `nc`.
                            if graph.arcs.len() != (nc as usize) + 1 {
                                let last = graph.arcs.len() - 1;
                                graph.arcs.swap(last, nc as usize);
                            }
                            nc += 1;
                        }
                    }
                }
                i += 1;
            }

            if cc > 0 {
                // new arcs must be bridged by contained reads
                nv = cc;
                avi = graph.arcs.len() as u32 - cc;
            } else {
                nv = 0;
                avi = u32::MAX;
            }
        }

        for hit in &overlaps[(v >> 1) as usize] {
            if hit.deleted {
                continue;
            }
            let Some(p) = hit_to_arc(graph, hit, max_hang, max_hang_rate, min_ovlp) else {
                continue;
            };
            if (p.ul >> 32) != u64::from(v) {
                continue;
            }
            dis[p.v as usize] = 0;
        }
    }

    if cnt > 0 {
        graph.index.clear();
        graph.is_sorted = false;
        graph.cleanup();
        graph.symmetrize();
    }
}

/// `ma_hit2arc` against the graph's read lengths.
fn hit_to_arc(
    graph: &Asg,
    hit: &PafHit,
    max_hang: i64,
    max_hang_rate: f64,
    min_ovlp: i64,
) -> Option<AsgArc> {
    let qlen = graph.sequences[hit.query_identifier as usize].length as i32;
    let tlen = graph.sequences[hit.target_identifier as usize].length as i32;
    let (result, arc) = ma_hit_to_arc(
        hit,
        qlen,
        tlen,
        max_hang as i32,
        max_hang_rate as f32,
        min_ovlp as i32,
    );
    if result < 0 {
        None
    } else {
        Some(arc)
    }
}

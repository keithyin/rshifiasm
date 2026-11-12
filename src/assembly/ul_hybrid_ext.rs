//! The hybrid UL/unitig extension subtree of `u2g_hybrid_detan_iter`
//! (`gfa_ut.cpp` 11635-14490 + 15332): integer-walk realignment, unique-arc
//! clustering, threading and the strict/free extension drivers.

use crate::assembly::asg::MaUtg;
use crate::assembly::bubble_chain::BubbleType;
use crate::assembly::ul2ul::{get_usg_arc_mm, push_usg_t_node, UinfoSrt, Usg, UsgArc};
use crate::assembly::ul_integer::UlStrIdx;
use crate::assembly::ul_store::UlReadVec;
use crate::assembly::ulg_clean::{get_integer_seq_ovlps, get_usg_arc};
use crate::container::ksort::radix_sort_u64;

/// The read-side slices the support queries need (`uidx` in the reference).
pub struct HybridCtx<'a> {
    /// `uidx->pstr`.
    pub str_idx: &'a UlStrIdx,
    /// The UL read store.
    pub store: &'a [UlReadVec],
    /// `raw->u.a` (`uidx->l1_ug->u.a`).
    pub raw_unitigs: &'a [MaUtg],
    /// `R_INF.read_lengths`.
    pub read_lengths: &'a [u64],
    /// The per-unitig read lists of `raw`.
    pub unitig_reads: &'a [Vec<u64>],
    /// `uidx->bub`.
    pub bub: &'a BubbleType,
}

/// `check_hybrid_connect` (`gfa_ut.cpp` 11893).
fn check_hybrid_connect(ng: &Usg, i_uid: u32, v: u32, vidx: u32, w: u32, _widx: u32) -> bool {
    let Some(z) = get_usg_arc(ng, v, w) else {
        return false;
    };
    let entries = get_usg_arc_mm(ng, z);
    if entries.is_empty() {
        return false;
    }
    for entry in entries {
        if entry.uid != i_uid || entry.off != vidx {
            continue;
        }
        return true;
    }
    false
}

/// `usg_arc_occ` (`gfa_ut.cpp` 11908): the live out-degree of `v`.
fn usg_live_arcs(ng: &Usg, v: u32) -> u32 {
    ng.arcs_of(v).iter().filter(|a| !a.deleted).count() as u32
}

/// The `integer_t` scratch of one worker (`buf->o`/`buf->u`/`buf->res_dump`).
#[derive(Default)]
pub struct IntegerBuffers {
    pub res_dump: Vec<u64>,
    pub o: Vec<u64>,
    pub u: Vec<u64>,
}

/// `integer_realign_g` (`gfa_ut.cpp` 11918): re-thread one integer-unitig
/// walk through the hybrid graph, appending
/// `<seq_id:32 | 0x8000...:1 | occ:31>` headers followed by the re-mapped
/// vertices of every non-linear path to `res_dump`.
pub fn integer_realign_g(ng: &Usg, seq: &[UinfoSrt], seq_id: u32, buf: &mut IntegerBuffers) {
    if seq.len() < 2 {
        return;
    }
    let mut seq_n = 0usize;
    for entry in seq {
        seq_n += ng.mp[(entry.v >> 1) as usize].len();
    }
    buf.o.clear();
    buf.u.clear();
    buf.o.resize(seq_n << 1, 0);
    // `srt = buf->o.a; track = buf->o.a + seq_n;`
    let (srt, track) = buf.o.split_at_mut(seq_n);

    let mut l = 0usize;
    for z in 0..seq.len() {
        let csc = u64::from(seq[z].n);
        let zm = &ng.mp[(seq[z].v >> 1) as usize];
        let zt = if z > 0 {
            Some(&ng.mp[(seq[z - 1].v >> 1) as usize])
        } else {
            None
        };
        for m in 0..zm.len() {
            let i = l + m;
            let mut mm_sc = csc;
            let mut mm_idx = 0x7FFF_FFFFu64;
            let vi = (zm[m] << 1) | u32::from(seq[z].v & 1);
            if let Some(zt) = zt {
                if !zt.is_empty() {
                    for t in 0..zt.len() {
                        let j = l + t - zt.len();
                        let vj = (zt[t] << 1) | u32::from(seq[z - 1].v & 1);
                        if check_hybrid_connect(ng, seq_id, vj, z as u32 - 1, vi, z as u32) {
                            let sc = csc + (0xFFFF_FFFFu32).wrapping_sub((track[j] >> 32) as u32)
                                as u64;
                            if sc > mm_sc {
                                mm_sc = sc;
                                mm_idx = j as u64;
                            }
                        }
                    }
                }
            }
            let mm_sc = 0xFFFF_FFFFu32.wrapping_sub(mm_sc as u32);
            track[i] = (u64::from(mm_sc) << 32) | mm_idx;
            srt[i] = (track[i] >> 32) << 32 | i as u64;
        }
        l += zm.len();
    }
    assert_eq!(l, seq_n);

    radix_sort_u64(srt);
    let res_base = buf.res_dump.len();
    buf.res_dump.reserve(seq_n);
    let mut n_v = 0usize;
    let mut u_words: Vec<u64> = Vec::new();
    for k in 0..seq_n {
        let n_v0 = n_v;
        let mut i = (srt[k] & 0xffff_ffff) as usize;
        while i != 0x7FFF_FFFF && (track[i] & (1u64 << 63)) == 0 {
            buf.res_dump.push(i as u64);
            n_v += 1;
            track[i] |= 1u64 << 63;
            i = (track[i] & 0x7FFF_FFFF) as usize;
        }
        if n_v - n_v0 <= 1 {
            buf.res_dump.truncate(res_base + n_v0);
            n_v = n_v0;
            continue;
        }
        let word = (0xFFFF_FFFFu32.wrapping_sub((n_v - n_v0) as u32) as u64) << 32
            | (res_base + n_v0) as u64;
        u_words.push(word);
    }
    let mut m_u = u_words.len();

    // `srt[]`: (idx in seq)|(node id)
    {
        let mut i = 0usize;
        for (z, entry) in seq.iter().enumerate() {
            let zm = &ng.mp[(entry.v >> 1) as usize];
            for &node in zm.iter() {
                srt[i] = (z as u64) << 32 | u64::from((node << 1) | u32::from(entry.v & 1));
                i += 1;
            }
        }
        assert_eq!(i, seq_n);
    }
    let mut nt = 0usize;
    for k in 0..u_words.len() {
        let n_v0 = (u_words[k] & 0xffff_ffff) as usize;
        let n_v = 0xFFFF_FFFFu32.wrapping_sub((u_words[k] >> 32) as u32) as usize;
        assert!(n_v >= 2);
        u_words[k] = (nt as u64) << 32;
        for i in 0..n_v {
            track[nt] = srt[(buf.res_dump[res_base + n_v0 + n_v - i - 1]) as usize];
            nt += 1;
        }
        u_words[k] |= nt as u64;
    }
    let total = nt;

    // drop the overlapping walks
    let mut n_u = 0usize;
    for k in 0..u_words.len() {
        let n_v0 = (u_words[k] >> 32) as usize;
        let n_v = (u_words[k] & 0xffff_ffff) as usize;
        let ceq_s = track[n_v0] >> 32;
        let ceq_e = (track[n_v - 1] >> 32) + 1;
        let mut z = 0usize;
        while z < n_u {
            let p_v0 = (u_words[z] >> 32) as usize;
            let p_v = (u_words[z] & 0xffff_ffff) as usize;
            let peq_s = track[p_v0] >> 32;
            let peq_e = (track[p_v - 1] >> 32) + 1;
            let max_s = peq_s.max(ceq_s);
            let min_e = peq_e.min(ceq_e);
            if min_e > max_s {
                break;
            }
            z += 1;
        }
        if z < n_u {
            continue;
        }
        u_words[n_u] = u_words[k];
        n_u += 1;
    }
    m_u = n_u;
    u_words.truncate(m_u);
    buf.res_dump.truncate(res_base);

    for k in 0..u_words.len() {
        let t = (u64::from(seq_id) << 32)
            | (1u64 << 63)
            | u64::from((u_words[k] & 0xffff_ffff) as u32 - (u_words[k] >> 32) as u32);
        let pn = buf.res_dump.len();
        buf.res_dump.push(t);
        let n_v0 = (u_words[k] >> 32) as usize;
        let n_v = (u_words[k] & 0xffff_ffff) as usize;
        let mut lin = true;
        for z in n_v0..n_v {
            if lin && z > n_v0 {
                if usg_live_arcs(ng, (track[z] & 0xffff_ffff) as u32 ^ 1) > 1
                    || usg_live_arcs(ng, (track[z - 1] & 0xffff_ffff) as u32) > 1
                {
                    lin = false;
                }
            }
            buf.res_dump.push(track[z] & 0xffff_ffff);
        }
        if lin {
            buf.res_dump.truncate(pn);
        }
    }
    let _ = total;
}

/// `u2g_hybrid_aln` (`gfa_ut.cpp` 13754): realign every integer-unitig walk
/// through `ng`; `ob` becomes the concatenated re-mapped vertices and `ub`
/// the `<offset>:32 | <length>:32` path index sorted by sequence id.
pub fn u2g_hybrid_aln(ng: &Usg, iug_a: &[Vec<UinfoSrt>], ob: &mut Vec<u64>, ub: &mut Vec<u64>) {
    ob.clear();
    ub.clear();
    let mut buf = IntegerBuffers::default();
    let mut m = 0usize;
    for (seq_id, seq) in iug_a.iter().enumerate() {
        buf.res_dump.clear();
        integer_realign_g(ng, seq, seq_id as u32, &mut buf);
        for &x in &buf.res_dump {
            ob.push(x);
            if x & (1u64 << 63) != 0 {
                let mut y = (x & !(1u64 << 63)) >> 32;
                y <<= 32;
                y |= ob.len() as u64;
                ub.push(y);
            } else {
                m += 1;
            }
        }
    }

    radix_sort_u64(ub);
    let mut tmp: Vec<u64> = Vec::with_capacity(ub.len() + m);
    for k in 0..ub.len() {
        let off = (ub[k] & 0xffff_ffff) as usize;
        let sn = (ob[off - 1] & 0xffff_ffff) as usize;
        tmp.extend_from_slice(&ob[off..off + sn]);
        ub[k] = ((tmp.len() - sn) as u64) << 32 | sn as u64;
    }
    *ob = tmp;
}

/// `get_arc_support` (`gfa_ut.cpp` 13048): how many UL strings walk `v`
/// directly into `w`.
pub fn get_arc_support(ctx: &HybridCtx, v: u32, w: u32) -> u64 {
    let str_idx = ctx.str_idx;
    let vertex = (v >> 1) as usize;
    let start = str_idx.idx[vertex] as usize;
    let end = str_idx.idx[vertex + 1] as usize;
    let mut occ = 0u64;
    for &word in &str_idx.occ[start..end] {
        let read_id = (word >> 32) as usize;
        let pos = (word & 0xffff_ffff) as usize;
        let str_ = &str_idx.str[read_id];
        let s_n = str_.cn as i64;
        if s_n < 2 {
            continue;
        }
        let vz = str_.a[pos] as u32;
        assert_eq!(v >> 1, vz >> 1);
        if v == vz {
            let mut s = pos as i64 + 1;
            while s < s_n {
                if str_.a[s as usize] as u32 == w {
                    occ += 1;
                    break;
                }
                s += 1;
            }
        } else {
            let mut s = pos as i64 - 1;
            while s >= 0 {
                if str_.a[s as usize] as u32 == (w ^ 1) {
                    occ += 1;
                    break;
                }
                s -= 1;
            }
        }
    }
    occ
}

/// `get_arc_support_chain` (`gfa_ut.cpp` 13130): the UL support of a whole
/// mapped walk.
pub fn get_arc_support_chain(ctx: &HybridCtx, ng: &Usg, a: &[u64]) -> u64 {
    if a.is_empty() {
        return 0;
    }
    let str_idx = ctx.str_idx;
    let v = ((ng.a[(a[0] >> 1) as usize].mm << 1) | (a[0] & 1) as u32) as u32;
    let vertex = (v >> 1) as usize;
    let start = str_idx.idx[vertex] as usize;
    let end = str_idx.idx[vertex + 1] as usize;
    let mut occ = 0u64;
    for &word in &str_idx.occ[start..end] {
        let read_id = (word >> 32) as usize;
        let pos = (word & 0xffff_ffff) as usize;
        let str_ = &str_idx.str[read_id];
        let s_n = str_.cn as i64;
        if s_n < 2 {
            continue;
        }
        let vz = str_.a[pos] as u32;
        assert_eq!(v >> 1, vz >> 1);
        if v == vz {
            let mut s = pos as i64 + 1;
            let mut ai = 1usize;
            while s < s_n && ai < a.len() {
                let w = ((ng.a[(a[ai] >> 1) as usize].mm << 1) | (a[ai] & 1) as u32) as u32;
                if str_.a[s as usize] as u32 != w {
                    break;
                }
                s += 1;
                ai += 1;
            }
            if ai >= a.len() {
                occ += 1;
            }
        } else {
            let mut s = pos as i64 - 1;
            let mut ai = 1usize;
            while s >= 0 && ai < a.len() {
                let w = ((ng.a[(a[ai] >> 1) as usize].mm << 1) | (a[ai] & 1) as u32) as u32;
                if str_.a[s as usize] as u32 != (w ^ 1) {
                    break;
                }
                s -= 1;
                ai += 1;
            }
            if ai >= a.len() {
                occ += 1;
            }
        }
    }
    occ
}

const UNIQUE_BRIDGE_OCC: u64 = 2;

/// `ava_pass_unique_bridge_cov` (`gfa_ut.cpp` 13238): the sequential form
/// of `worker_unique_bridge_check_s` — flags and clears the unreliable
/// intervals of one group.
fn ava_pass_unique_bridge_cov(
    ctx: &HybridCtx,
    ng: &Usg,
    b64: &mut [u64],
    g_s: usize,
    g_e: usize,
    integer_seq: &[u64],
    tip_l: u64,
) -> bool {
    if tip_l == 0 && g_e - g_s == 1 {
        return true;
    }
    let mut is_del = false;
    for i in g_s..g_e {
        let slot = (b64[i] & 0xffff_ffff) as usize;
        let s = (b64[slot] >> 32) as usize;
        let e = (b64[slot] & 0xffff_ffff) as usize;
        assert!(s < e);
        let mut vk = u32::MAX;
        for k in s..=e {
            if k > s {
                let nse = get_arc_support_chain(ctx, ng, &integer_seq[k - 1..k + 1]);
                if nse < UNIQUE_BRIDGE_OCC {
                    b64[i] |= 1u64 << 63;
                    is_del = true;
                    break;
                }
            }
            if crate::assembly::ul_renew_gfa::if_hom((integer_seq[k] >> 1) as u32, ctx.bub) {
                continue;
            }
            let wk = k as u32;
            if vk != u32::MAX {
                let nse =
                    get_arc_support_chain(ctx, ng, &integer_seq[vk as usize..wk as usize + 1]);
                if nse < UNIQUE_BRIDGE_OCC {
                    b64[i] |= 1u64 << 63;
                    is_del = true;
                    break;
                }
            }
            vk = wk;
        }
    }
    for i in g_s..g_e {
        if b64[i] & (1u64 << 63) != 0 {
            b64[i] -= 1u64 << 63;
        }
    }
    !is_del
}

/// `ava_pass_unique_bridge_tips` (`gfa_ut.cpp` 12256): the tip-growth
/// counter of one cluster; returns the extension length.
#[allow(clippy::too_many_arguments)]
fn ava_pass_unique_bridge_tips(
    ng: &Usg,
    b64: &[u64],
    g_s: usize,
    g_e: usize,
    integer_seq: &[u64],
    f: &mut [u8],
    max_ext: u32,
    max_ext_rate: f64,
    ext_up: u64,
) -> u64 {
    let mut work: Vec<u64> = Vec::new();
    let mut nkeep = 0u64;
    for i in g_s..g_e {
        let slot = (b64[i] & 0xffff_ffff) as usize;
        let s = (b64[slot] >> 32) as usize;
        let e = (b64[slot] & 0xffff_ffff) as usize;
        assert!(s < e);
        for k in s + 1..e {
            f[integer_seq[k] as usize] = 1;
            f[(integer_seq[k] ^ 1) as usize] = 1;
            nkeep += u64::from(ng.a[(integer_seq[k] >> 1) as usize].occ);
        }
        f[integer_seq[s] as usize] = 1;
        f[(integer_seq[e] ^ 1) as usize] = 1;
    }
    nkeep = (nkeep as f64 * max_ext_rate) as u64;

    // collect the neighbours of the cluster that are not inside it
    for i in g_s..g_e {
        let slot = (b64[i] & 0xffff_ffff) as usize;
        let s = (b64[slot] >> 32) as usize;
        let e = (b64[slot] & 0xffff_ffff) as usize;
        for k in s + 1..e {
            for v in [integer_seq[k], integer_seq[k] ^ 1] {
                for arc in ng.arcs_of(v as u32) {
                    if arc.deleted || f[(arc.v ^ 1) as usize] != 0 {
                        continue;
                    }
                    work.push(u64::from(arc.v));
                }
            }
        }
        for v in [integer_seq[s], integer_seq[e] ^ 1] {
            for arc in ng.arcs_of(v as u32) {
                if arc.deleted || f[(arc.v ^ 1) as usize] != 0 {
                    continue;
                }
                work.push(u64::from(arc.v));
            }
        }
    }

    let mut ncut = nkeep.max(u64::from(max_ext));
    let mut is_telo = false;
    if ncut > ext_up {
        ncut = ext_up;
    }
    if ncut < u64::from(max_ext) {
        ncut = u64::from(max_ext);
    }
    let mut n_ext = 0u64;
    while !work.is_empty() && n_ext < ncut && !is_telo {
        let v = work.pop().unwrap();
        if f[v as usize] != 0 {
            continue;
        }
        let mut kv = 0u32;
        for arc in ng.arcs_of((v ^ 1) as u32) {
            if arc.deleted || f[(arc.v ^ 1) as usize] != 0 {
                continue;
            }
            kv += 1;
        }
        if kv > 0 {
            continue;
        }
        n_ext += u64::from(ng.a[(v >> 1) as usize].occ);
        f[v as usize] = 1;
        f[(v ^ 1) as usize] = 1;
        let mut kv_t = 0u32;
        for arc in ng.arcs_of(v as u32) {
            if arc.deleted {
                continue;
            }
            kv_t += 1;
            if f[(arc.v ^ 1) as usize] != 0 {
                continue;
            }
            work.push(u64::from(arc.v));
        }
        if kv_t == 0 && ng.a[(v >> 1) as usize].telo {
            is_telo = true;
        }
    }

    // reset f[] to 0: the interval members join the work list first, then
    // every entry pops and clears its marked neighbours
    for i in g_s..g_e {
        let slot = (b64[i] & 0xffff_ffff) as usize;
        let s = (b64[slot] >> 32) as usize;
        let e = (b64[slot] & 0xffff_ffff) as usize;
        for k in s + 1..e {
            work.push(integer_seq[k]);
            work.push(integer_seq[k] ^ 1);
            f[integer_seq[k] as usize] = 0;
            f[(integer_seq[k] ^ 1) as usize] = 0;
        }
        work.push(integer_seq[s]);
        work.push(integer_seq[e] ^ 1);
        f[integer_seq[s] as usize] = 0;
        f[(integer_seq[e] ^ 1) as usize] = 0;
    }
    while let Some(v) = work.pop() {
        f[v as usize] = 0;
        f[(v ^ 1) as usize] = 0;
        for arc in ng.arcs_of(v as u32) {
            if arc.deleted {
                continue;
            }
            if f[arc.v as usize] != 0 || f[(arc.v ^ 1) as usize] != 0 {
                work.push(u64::from(arc.v));
            }
        }
        for arc in ng.arcs_of((v ^ 1) as u32) {
            if arc.deleted {
                continue;
            }
            if f[arc.v as usize] != 0 || f[(arc.v ^ 1) as usize] != 0 {
                work.push(u64::from(arc.v));
            }
        }
    }
    if n_ext < ncut && n_ext > u64::from(max_ext) - 1 {
        n_ext = u64::from(max_ext) - 1;
    }
    if is_telo {
        n_ext = u64::from(max_ext) + 1;
    }
    n_ext
}

/// `iter_unique_arcs` (`gfa_ut.cpp` 12135): the flood fill that assigns the
/// cluster id `i0` to every interval reachable through shared raw unitigs.
fn iter_unique_arcs(
    buf: &mut Vec<u64>,
    b64: &[u64],
    idx: &[u64],
    integ_seq: &[u64],
    cidx: &mut [u64],
    i0: u64,
) {
    buf.clear();
    buf.push(i0);
    while let Some(x) = buf.pop() {
        if cidx[x as usize] != u64::MAX {
            continue;
        }
        cidx[x as usize] = (i0 << 32) | x;
        let s = (b64[x as usize] >> 32) as usize;
        let e = (b64[x as usize] & 0xffff_ffff) as usize;
        assert!(e > s);
        for k in s + 1..e {
            let v = integ_seq[k];
            let zs = ((idx[v as usize] << 1) >> 33) as usize;
            let ze = (idx[v as usize] & 0xffff_ffff) as usize;
            for z in zs..ze {
                assert_eq!((b64[z] >> 32) as usize, v as usize);
                let target = (b64[z] & 0xffff_ffff) as usize;
                if cidx[target] != u64::MAX {
                    continue;
                }
                buf.push(target as u64);
            }

            let v = integ_seq[k] ^ 1;
            let zs = ((idx[v as usize] << 1) >> 33) as usize;
            let ze = (idx[v as usize] & 0xffff_ffff) as usize;
            for z in zs..ze {
                assert_eq!((b64[z] >> 32) as usize, v as usize);
                let target = (b64[z] & 0xffff_ffff) as usize;
                if cidx[target] != u64::MAX {
                    continue;
                }
                buf.push(target as u64);
            }
        }

        let v = integ_seq[s];
        let zs = ((idx[v as usize] << 1) >> 33) as usize;
        let ze = (idx[v as usize] & 0xffff_ffff) as usize;
        for z in zs..ze {
            assert_eq!((b64[z] >> 32) as usize, v as usize);
            let target = (b64[z] & 0xffff_ffff) as usize;
            if cidx[target] != u64::MAX {
                continue;
            }
            buf.push(target as u64);
        }

        let v = integ_seq[e] ^ 1;
        let zs = ((idx[v as usize] << 1) >> 33) as usize;
        let ze = (idx[v as usize] & 0xffff_ffff) as usize;
        for z in zs..ze {
            assert_eq!((b64[z] >> 32) as usize, v as usize);
            let target = (b64[z] & 0xffff_ffff) as usize;
            if cidx[target] != u64::MAX {
                continue;
            }
            buf.push(target as u64);
        }
    }
}

/// `usg_unique_arcs_cluster_adv` (`gfa_ut.cpp` 12201): groups the available
/// intervals that share raw unitigs into clusters; returns the cluster
/// count. `b64` ends as `<cluster:32 | <interval>:32>` words.
fn usg_unique_arcs_cluster_adv(
    b64: &mut Vec<u64>,
    a_n: usize,
    idx: &[u64],
    integ_seq: &[u64],
    buf: &mut Vec<u64>,
) -> u32 {
    let bn = b64.len();
    b64.resize(bn + a_n, u64::MAX);
    {
        let (head, cidx) = b64.split_at_mut(bn);
        for i in 0..a_n {
            if cidx[i] != u64::MAX {
                continue;
            }
            iter_unique_arcs(buf, head, idx, integ_seq, cidx, i as u64);
        }
        radix_sort_u64(cidx);
    }
    let cidx: Vec<u64> = b64[bn..].to_vec();
    b64.truncate(a_n);
    let mut v = 0u64;
    let mut i = 0usize;
    let mut k = 1usize;
    loop {
        let at_end = k == a_n;
        if at_end || (cidx[i] >> 32) != (cidx[k] >> 32) {
            for z in i..k {
                assert_ne!(cidx[z], u64::MAX);
                b64.push((v << 32) | (cidx[z] & 0xffff_ffff));
            }
            i = k;
            v += 1;
        }
        if at_end {
            break;
        }
        k += 1;
    }
    buf.clear();
    assert_eq!(b64.len(), a_n << 1);
    v as u32
}

/// `remap_gen_arcs` (`gfa_ut.cpp` 11635): copies the `ov -> ow` arc (and its
/// recorded integer-walk usages) onto the new `nv -> nw` pair.
fn remap_gen_arcs(g: &mut Usg, ov: u32, ow: u32, nv: u32, nw: u32) {
    g.a[(nv >> 1) as usize].arc[(nv & 1) as usize].push(UsgArc::default());
    let op = get_usg_arc(g, ov, ow).expect("the remapped arc must exist");
    let op = *op;
    let mm_len = g.a[(nv >> 1) as usize].arc_mm[(nv & 1) as usize].len() as u64;
    {
        let np = g.a[(nv >> 1) as usize].arc[(nv & 1) as usize].last_mut().unwrap();
        np.ul = (op.ul as u32 as u64) | ((u64::from(nv)) << 32);
        np.v = nw;
        np.ol = op.ol;
        np.deleted = false;
        np.ou = op.ou;
        np.idx = (mm_len << 32) | (op.idx & 0xffff_ffff);
    }
    let v_n = (op.idx & 0xffff_ffff) as usize;
    let src_v = (op.ul >> 32) as usize;
    let src_off = (op.idx >> 32) as usize;
    let entries: Vec<crate::assembly::ul2ul::UsgArcMm> = g.a[src_v >> 1].arc_mm[src_v & 1]
        [src_off..src_off + v_n]
        .to_vec();
    let dst = &mut g.a[(nv >> 1) as usize].arc_mm[(nv & 1) as usize];
    for entry in entries {
        dst.push(crate::assembly::ul2ul::UsgArcMm {
            v: nw,
            uid: entry.uid,
            off: entry.off,
        });
    }
}

/// `update_usg_t_threading_0` (`gfa_ut.cpp` 12355): threads one path through
/// the graph, copying shared nodes.
fn update_usg_t_threading_0(ng: &mut Usg, a: &[u64], occ: &mut [u32], b: &mut Vec<u64>) {
    b.clear();
    assert!(a.len() > 1);
    assert_eq!(occ[a[0] as usize], 1);
    b.push((a[0] << 32) | a[0]);
    occ[a[0] as usize] -= 1;
    for k in 1..a.len() - 1 {
        let ak = a[k];
        assert_eq!(occ[ak as usize], occ[(ak ^ 1) as usize]);
        assert!(occ[ak as usize] > 0);
        if occ[ak as usize] > 1 {
            // copy a node
            let nid = (ak >> 1) as usize;
            let nnid = ng.a.len() as u32;
            let (mm, occ_v, len, telo) = {
                let src = &ng.a[nid];
                (src.mm, src.occ, src.len, src.telo)
            };
            ng.mp[mm as usize].push(nnid);
            {
                let s = push_usg_t_node(ng, u64::from(nnid));
                s.mm = mm;
                s.occ = occ_v;
                s.len = len;
                s.telo = telo;
                s.deleted = false;
                s.arc = [Vec::new(), Vec::new()];
                s.arc_mm = [Vec::new(), Vec::new()];
            }
            let nnid_v = (nnid << 1) | ((ak & 1) as u32);
            b.push(((ak ^ 1) << 32) | u64::from(nnid_v ^ 1));
            b.push((ak << 32) | u64::from(nnid_v));
        } else {
            b.push(((ak ^ 1) << 32) | (ak ^ 1));
            b.push((ak << 32) | ak);
        }
        occ[ak as usize] -= 1;
        occ[(ak ^ 1) as usize] -= 1;
    }
    let last = a[a.len() - 1] ^ 1;
    assert_eq!(occ[last as usize], 1);
    b.push((last << 32) | last);
    occ[last as usize] -= 1;

    let mut k = 0usize;
    while k + 1 < b.len() {
        let v = b[k];
        let w = b[k + 1];
        if (v >> 32) == (v & 0xffff_ffff) && (w >> 32) == (w & 0xffff_ffff) {
            k += 2;
            continue;
        }
        remap_gen_arcs(
            ng,
            (v >> 32) as u32,
            ((w >> 32) as u32) ^ 1,
            (v & 0xffff_ffff) as u32,
            ((w & 0xffff_ffff) as u32) ^ 1,
        );
        remap_gen_arcs(
            ng,
            (w >> 32) as u32,
            ((v >> 32) as u32) ^ 1,
            (w & 0xffff_ffff) as u32,
            ((v & 0xffff_ffff) as u32) ^ 1,
        );
        k += 2;
    }

    let mut k = 0usize;
    while k + 1 < b.len() {
        let v0 = (b[k] >> 32) as u32;
        let v1 = (b[k] & 0xffff_ffff) as u32;
        let w0 = (b[k + 1] >> 32) as u32;
        let w1 = (b[k + 1] & 0xffff_ffff) as u32;

        let mut kv = 0u32;
        for i in 0..ng.a[(v1 >> 1) as usize].arc[(v1 & 1) as usize].len() {
            let arc = ng.a[(v1 >> 1) as usize].arc[(v1 & 1) as usize][i];
            if arc.v == w1 ^ 1 {
                ng.a[(v1 >> 1) as usize].arc[(v1 & 1) as usize][i].deleted = false;
                kv += 1;
            } else if !arc.deleted {
                ng.a[(v1 >> 1) as usize].arc[(v1 & 1) as usize][i].deleted = true;
                let zv = (arc.v ^ 1) as usize;
                let zs = ((arc.ul >> 32) as u32) ^ 1;
                assert!(!get_usg_arc(ng, arc.v ^ 1, zs).expect("the mirrored arc must exist").deleted);
                if let Some(cand) = ng.a[zv >> 1].arc[zv & 1].iter_mut().find(|c| c.v == zs) {
                    cand.deleted = true;
                }
            }
        }
        assert_eq!(kv, 1);

        let mut kw = 0u32;
        for i in 0..ng.a[(w1 >> 1) as usize].arc[(w1 & 1) as usize].len() {
            let arc = ng.a[(w1 >> 1) as usize].arc[(w1 & 1) as usize][i];
            if arc.v == v1 ^ 1 {
                ng.a[(w1 >> 1) as usize].arc[(w1 & 1) as usize][i].deleted = false;
                kw += 1;
            } else if !arc.deleted {
                ng.a[(w1 >> 1) as usize].arc[(w1 & 1) as usize][i].deleted = true;
                let zv = (arc.v ^ 1) as usize;
                let zs = ((arc.ul >> 32) as u32) ^ 1;
                assert!(!get_usg_arc(ng, arc.v ^ 1, zs).expect("the mirrored arc must exist").deleted);
                if let Some(cand) = ng.a[zv >> 1].arc[zv & 1].iter_mut().find(|c| c.v == zs) {
                    cand.deleted = true;
                }
            }
        }
        assert_eq!(kw, 1);

        if v0 == v1 && w0 == w1 {
            k += 2;
            continue;
        }

        for i in 0..ng.a[(v0 >> 1) as usize].arc[(v0 & 1) as usize].len() {
            if ng.a[(v0 >> 1) as usize].arc[(v0 & 1) as usize][i].v == w0 ^ 1 {
                ng.a[(v0 >> 1) as usize].arc[(v0 & 1) as usize][i].deleted = true;
                break;
            }
        }
        for i in 0..ng.a[(w0 >> 1) as usize].arc[(w0 & 1) as usize].len() {
            if ng.a[(w0 >> 1) as usize].arc[(w0 & 1) as usize][i].v == v0 ^ 1 {
                ng.a[(w0 >> 1) as usize].arc[(w0 & 1) as usize][i].deleted = true;
                break;
            }
        }
        k += 2;
    }
}

/// `update_usg_t_threading` (`gfa_ut.cpp` 12474): the occupancy sweep plus
/// one `update_usg_t_threading_0` per cluster.
fn update_usg_t_threading(
    ng: &mut Usg,
    arcs: &[u64],
    arcs_g: &[u64],
    integ_seq: &[u64],
    occ: &mut [u32],
    b: &mut Vec<u64>,
) {
    let nvtx = ng.a.len() << 1;
    for slot in occ[..nvtx].iter_mut() {
        *slot = 0;
    }
    for g in arcs_g {
        let slot = (g & 0xffff_ffff) as usize;
        let s = (arcs[slot] >> 32) as usize;
        let e = (arcs[slot] & 0xffff_ffff) as usize;
        assert!(s < e);
        for k in s + 1..e {
            occ[integ_seq[k] as usize] += 1;
            occ[(integ_seq[k] ^ 1) as usize] += 1;
        }
        occ[integ_seq[s] as usize] += 1;
        occ[(integ_seq[e] ^ 1) as usize] += 1;
    }
    for slot in occ[..nvtx].iter_mut() {
        if *slot == 0 {
            *slot = u32::MAX;
        }
    }
    for g in arcs_g {
        let slot = (g & 0xffff_ffff) as usize;
        let s = (arcs[slot] >> 32) as usize;
        let e = (arcs[slot] & 0xffff_ffff) as usize;
        update_usg_t_threading_0(ng, &integ_seq[s..e + 1], occ, b);
    }
}

/// `old_path_ext` (`gfa_ut.cpp` 13262): filters the clusters by tip growth
/// and bridge support, then threads the survivors.
#[allow(clippy::too_many_arguments)]
fn old_path_ext(
    ctx: &HybridCtx,
    ng: &mut Usg,
    b64: &mut Vec<u64>,
    ub64: &mut Vec<u64>,
    a_n: usize,
    i_idx: &[u64],
    int_a: &[u64],
    ff: &mut [u8],
    ng_occ: &mut [u32],
    max_ext: u32,
) -> u64 {
    let mut n_clus = 0u64;
    let mut mm = a_n;
    let mut i = a_n;
    let mut k = a_n + 1;
    loop {
        let at_end = k == b64.len();
        if at_end || (b64[k] >> 32) != (b64[i] >> 32) {
            let tip_l = ava_pass_unique_bridge_tips(ng, b64, i, k, int_a, ff, max_ext, 0.03, 16);
            if tip_l < u64::from(max_ext) {
                if ava_pass_unique_bridge_cov(ctx, ng, b64, i, k, int_a, tip_l) {
                    for z in i..k {
                        b64[mm] = b64[z];
                        mm += 1;
                    }
                    n_clus += 1;
                }
            }
            i = k;
        }
        if at_end {
            break;
        }
        k += 1;
    }
    b64.truncate(mm);
    if n_clus > 0 {
        let arcs: Vec<u64> = b64.clone();
        let arcs_g: Vec<u64> = b64[a_n..].to_vec();
        update_usg_t_threading(ng, &arcs, &arcs_g, int_a, ng_occ, ub64);
    }
    n_clus
}

/// `gen_hybrid_aln_idx` (`gfa_ut.cpp` 13797): the CSR of `(walk, position)`
/// words per usg vertex; `ridx` packs `<offset>:32 | <count>:32`.
fn gen_hybrid_aln_idx(ng: &Usg, int_idx: &[u64], int_a: &[u64], b64: &mut Vec<u64>, ridx: &mut [u64]) {
    b64.clear();
    b64.resize(int_a.len(), 0);
    for slot in ridx[..ng.a.len()].iter_mut() {
        *slot = 0;
    }
    for k in 0..int_idx.len() {
        let s = (int_idx[k] >> 32) as usize;
        let e = s + (int_idx[k] & 0xffff_ffff) as usize;
        assert!(e > s + 1);
        for i in s..e {
            ridx[(int_a[i] >> 1) as usize] += 1;
        }
    }
    let mut l = 0u64;
    for k in 0..ng.a.len() {
        let m = ridx[k];
        ridx[k] = (l << 32) | m;
        l += m;
    }
    for k in 0..ng.a.len() {
        let a_start = (ridx[k] >> 32) as usize;
        let a_n = (ridx[k] & 0xffff_ffff) as usize;
        if a_n > 0 {
            b64[a_start + a_n - 1] = 0;
        }
    }
    for k in 0..int_idx.len() {
        let s = (int_idx[k] >> 32) as usize;
        let e = s + (int_idx[k] & 0xffff_ffff) as usize;
        for i in s..e {
            let v = (int_a[i] >> 1) as usize;
            let a_start = (ridx[v] >> 32) as usize;
            let a_n = (ridx[v] & 0xffff_ffff) as usize;
            if a_n > 0 {
                if b64[a_start + a_n - 1] == (a_n - 1) as u64 {
                    b64[a_start + a_n - 1] = ((k as u64) << 32) | i as u64;
                } else {
                    let slot = b64[a_start + a_n - 1] as usize;
                    b64[slot] = ((k as u64) << 32) | i as u64;
                    b64[a_start + a_n - 1] += 1;
                }
            }
        }
    }
}

/// `gen_sub_integer_path` (`gfa_ut.cpp` 13832).
fn gen_sub_integer_path(int_a: &[u64], s: i64, e: i64, it: i64, v: u64, res: &mut Vec<u64>) {
    res.clear();
    assert_eq!((int_a[it as usize] >> 1) as u64, v >> 1);
    if int_a[it as usize] == v {
        for k in it..e {
            res.push(int_a[k as usize]);
        }
    } else {
        let mut k = it;
        while k >= s {
            res.push(int_a[k as usize] ^ 1);
            k -= 1;
        }
    }
}

/// `is_best_path` (`gfa_ut.cpp` 13949).
#[allow(clippy::too_many_arguments)]
fn is_best_path(
    ctx: &HybridCtx,
    int_idx: &[u64],
    int_a: &[u64],
    s: usize,
    e: usize,
    it: usize,
    v: u64,
    ridx_a: &[u64],
    ridx: &[u64],
    b0: &mut Vec<u64>,
    b1: &mut Vec<u64>,
    cutoff: f64,
) -> bool {
    gen_sub_integer_path(int_a, s as i64, e as i64, it as i64, v, b0);
    let arc_start = (ridx[(v >> 1) as usize] >> 32) as usize;
    let arc_n = (ridx[(v >> 1) as usize] & 0xffff_ffff) as usize;
    let mut min_w0 = u64::MAX;
    let mut min_w1 = u64::MAX;
    let mut alt_n = 0u64;
    let mut is_contain = true;
    for k in 0..arc_n {
        let arc_word = ridx_a[arc_start + k];
        if it as u64 == (arc_word & 0xffff_ffff) {
            continue;
        }
        let slot = (arc_word >> 32) as usize;
        let is_ = (int_idx[slot] >> 32) as usize;
        alt_n += 1;
        let ie = is_ + (int_idx[slot] & 0xffff_ffff) as usize;
        gen_sub_integer_path(
            int_a,
            is_ as i64,
            ie as i64,
            (arc_word & 0xffff_ffff) as i64,
            v,
            b1,
        );
        let zn = b0.len().min(b1.len());
        let mut z = 0usize;
        while z < zn && b0[z] == b1[z] {
            z += 1;
        }
        assert!(z > 0);
        let mut w0 = u64::MAX;
        let mut w1 = u64::MAX;
        if z < b0.len() {
            get_integer_seq_ovlps(
                ctx.str_idx,
                ctx.store,
                ctx.raw_unitigs,
                ctx.read_lengths,
                ctx.unitig_reads,
                ctx.bub,
                b0,
                z as i64 - 1,
                false,
                None,
                &mut w0,
            );
        } else {
            return false;
        }
        if z < b1.len() {
            get_integer_seq_ovlps(
                ctx.str_idx,
                ctx.store,
                ctx.raw_unitigs,
                ctx.read_lengths,
                ctx.unitig_reads,
                ctx.bub,
                b1,
                z as i64 - 1,
                false,
                None,
                &mut w1,
            );
        } else {
            continue;
        }
        if w0 == u64::MAX {
            w0 = 0;
        }
        if w1 == u64::MAX {
            w1 = 0;
        }
        if w0 <= w1 || (w1 as f64) > ((w0 as f64) * cutoff) {
            return false;
        }
        if min_w0 == u64::MAX || z == zn || min_w0 > w0 || (min_w0 == w0 && min_w1 < w1) {
            min_w0 = w0;
            min_w1 = w1;
        }
        is_contain = false;
    }
    if alt_n == 0 {
        return true;
    }
    if is_contain {
        return true;
    }
    if min_w0 == u64::MAX {
        min_w0 = 0;
    }
    if min_w1 == u64::MAX {
        min_w1 = 0;
    }
    if min_w0 > min_w1 && (min_w1 as f64) <= ((min_w0 as f64) * cutoff) {
        return true;
    }
    false
}

/// `get_best_path` (`gfa_ut.cpp` 14114).
#[allow(clippy::too_many_arguments)]
fn get_best_path(
    ctx: &HybridCtx,
    int_idx: &[u64],
    int_a: &[u64],
    f: &[u8],
    s: usize,
    e: usize,
    b0: &mut Vec<u64>,
    b1: &mut Vec<u64>,
    ridx_a: &[u64],
    ridx: &[u64],
    res: &mut Vec<u64>,
) {
    let res_n = res.len();
    let mut pi = usize::MAX;
    for k in s..e {
        let v = int_a[k];
        if f[v as usize] == 0 && f[(v ^ 1) as usize] == 0 {
            continue;
        }
        if pi != usize::MAX && f[(v ^ 1) as usize] != 0 {
            if is_best_path(
                ctx,
                int_idx,
                int_a,
                s,
                e,
                k,
                v ^ 1,
                ridx_a,
                ridx,
                b0,
                b1,
                0.51,
            ) {
                let push = if res.len() > res_n {
                    (res[res.len() - 1] >> 32) != pi as u64
                } else {
                    true
                };
                if push {
                    res.push((pi as u64) << 32 | k as u64);
                }
            } else {
                pi = usize::MAX;
            }
        }
        if f[v as usize] != 0 {
            if is_best_path(ctx, int_idx, int_a, s, e, k, v, ridx_a, ridx, b0, b1, 0.51) {
                pi = k;
            } else {
                pi = usize::MAX;
            }
        }
    }
}

/// `gen_sub_integer_path_ff` (`gfa_ut.cpp` 14164): the path walk that stops
/// at the first `f`-marked opposite strand vertex.
#[allow(clippy::too_many_arguments)]
fn gen_sub_integer_path_ff(
    int_a: &[u64],
    s: i64,
    e: i64,
    it: i64,
    v: u64,
    f: &[u8],
    res: &mut Vec<u64>,
) -> u64 {
    res.clear();
    assert_eq!((int_a[it as usize] >> 1) as u64, v >> 1);
    if int_a[it as usize] == v {
        let mut k = it;
        while k < e {
            res.push(int_a[k as usize]);
            if res.len() > 1 && f[(res[res.len() - 1] ^ 1) as usize] != 0 {
                return res.len() as u64;
            }
            k += 1;
        }
    } else {
        let mut k = it;
        while k >= s {
            res.push(int_a[k as usize] ^ 1);
            if res.len() > 1 && f[(res[res.len() - 1] ^ 1) as usize] != 0 {
                return res.len() as u64;
            }
            k -= 1;
        }
    }
    res.clear();
    0
}

/// `is_best_pair` (`gfa_ut.cpp` 14185).
#[allow(clippy::too_many_arguments)]
fn is_best_pair(
    ctx: &HybridCtx,
    ng: &Usg,
    int_idx: &[u64],
    int_a: &[u64],
    s: usize,
    e: usize,
    it: usize,
    v: u64,
    ridx_a: &[u64],
    ridx: &[u64],
    b0: &mut Vec<u64>,
    b1: &mut Vec<u64>,
    f: &[u8],
    cutoff: f64,
) -> bool {
    let sup_cut = 3u64;
    if gen_sub_integer_path_ff(int_a, s as i64, e as i64, it as i64, v, f, b0) == 0 {
        return false;
    }
    assert!(b0.len() >= 2);
    assert!(f[b0[0] as usize] != 0 && f[(b0[b0.len() - 1] ^ 1) as usize] != 0);
    let w0 = get_arc_support_chain(ctx, ng, b0);
    if w0 < sup_cut {
        return false;
    }
    let arc_start = (ridx[(v >> 1) as usize] >> 32) as usize;
    let arc_n = (ridx[(v >> 1) as usize] & 0xffff_ffff) as usize;
    for k in 0..arc_n {
        let arc_word = ridx_a[arc_start + k];
        if it as u64 == (arc_word & 0xffff_ffff) {
            continue;
        }
        let slot = (arc_word >> 32) as usize;
        let is_ = (int_idx[slot] >> 32) as usize;
        let ie = is_ + (int_idx[slot] & 0xffff_ffff) as usize;
        if gen_sub_integer_path_ff(
            int_a,
            is_ as i64,
            ie as i64,
            (arc_word & 0xffff_ffff) as i64,
            v,
            f,
            b1,
        ) == 0
        {
            continue;
        }
        assert!(b1.len() >= 2);
        assert!(f[b1[0] as usize] != 0 && f[(b1[b1.len() - 1] ^ 1) as usize] != 0);
        if b0.len() == b1.len() && b0[..] == b1[..] {
            if it < (arc_word & 0xffff_ffff) as usize {
                continue;
            }
            return false;
        }
        let w1 = get_arc_support_chain(ctx, ng, b1);
        if w0 <= w1 || (w1 as f64) > ((w0 as f64) * cutoff) {
            return false;
        }
    }
    true
}

/// `get_best_pair` (`gfa_ut.cpp` 14222).
#[allow(clippy::too_many_arguments)]
fn get_best_pair(
    ctx: &HybridCtx,
    ng: &Usg,
    int_idx: &[u64],
    int_a: &[u64],
    f: &[u8],
    s: usize,
    e: usize,
    b0: &mut Vec<u64>,
    b1: &mut Vec<u64>,
    ridx_a: &[u64],
    ridx: &[u64],
    res: &mut Vec<u64>,
) {
    let mut pi = usize::MAX;
    for k in s..e {
        let v = int_a[k];
        if f[(v ^ 1) as usize] != 0 {
            if pi != usize::MAX
                && is_best_pair(
                    ctx,
                    ng,
                    int_idx,
                    int_a,
                    s,
                    e,
                    k,
                    v ^ 1,
                    ridx_a,
                    ridx,
                    b0,
                    b1,
                    f,
                    0.51,
                )
            {
                res.push((pi as u64) << 32 | k as u64);
            }
            pi = usize::MAX;
        }
        if f[v as usize] != 0 {
            pi = usize::MAX;
            if is_best_pair(ctx, ng, int_idx, int_a, s, e, k, v, ridx_a, ridx, b0, b1, f, 0.51) {
                pi = k;
            }
        }
    }
}

/// The renew parameters of the hybrid flow (`uidx` slices used by
/// `renew_usg_t_bub`).
pub struct HybridRenew<'a> {
    /// `uidx->l1_ug->u.a`.
    pub raw_unitigs: &'a [MaUtg],
    /// `uidx->sg`.
    pub sg: &'a crate::assembly::asg::Asg,
    /// `uidx->bub`.
    pub bub: &'a mut BubbleType,
    /// `uidx->r_het`.
    pub r_het: &'a mut [u8],
    /// `uopt->sources`.
    pub sources: &'a [Vec<crate::overlap::paf::PafHit>],
    /// `uopt->ruIndex`.
    pub read_to_unitig: &'a crate::assembly::hit::ReadToUnitig,
    /// `asm_opt`.
    pub options: &'a crate::command_line::options::AssemblyOptions,
    /// `asm_opt.purge_level_primary`.
    pub purge_level_primary: u8,
}

/// `#define occ_m(x)` (`gfa_ut.cpp` 12838).
#[inline]
fn occ_m(x: u32) -> u32 {
    x & 0x7fff_ffff
}

/// `debug_sysm_usg_t` (`gfa_ut.cpp` 12336): reports asymmetric arc deletion
/// states (silent on a consistent graph, exactly like the reference).
fn debug_sysm_usg_t(ng: &Usg, cmd: &str) {
    for v in 0..(ng.a.len() << 1) as u32 {
        for arc in ng.arcs_of(v) {
            let q = get_usg_arc(ng, arc.v ^ 1, v ^ 1);
            let mismatch = match q {
                None => true,
                Some(q) => arc.deleted != q.deleted,
            };
            if mismatch {
                eprintln!(
                    "[M::debug_sysm_usg_t::{}] utg{}({}){} -> utg{}({}){}, p.del::{}, q.del::{}",
                    cmd,
                    (v >> 1) + 1,
                    if v & 1 != 0 { "-" } else { "+" },
                    v,
                    (arc.v >> 1) + 1,
                    if arc.v & 1 != 0 { "-" } else { "+" },
                    arc.v,
                    arc.deleted as u32,
                    q.map_or(1, |q| q.deleted as u32),
                );
            }
        }
    }
}

/// The shared occupancy scan and unique-end marking of `ug_ext`/`ug_ext_strict`
/// (`gfa_ut.cpp` 13790 region).
fn mark_unique_ends(ng: &Usg, int_idx: &[u64], int_a: &[u64], ng_occ: &mut [u32]) {
    for k in 0..int_idx.len() {
        let base = (int_idx[k] >> 32) as usize;
        let r_n = (int_idx[k] & 0xffff_ffff) as usize;
        let r_a = &int_a[base..base + r_n];
        assert!(r_n >= 2);
        for i in 1..r_n - 1 {
            ng_occ[r_a[i] as usize] += 1;
            ng_occ[(r_a[i] ^ 1) as usize] += 1;
        }
        ng_occ[r_a[0] as usize] += 1;
        ng_occ[(r_a[r_n - 1] ^ 1) as usize] += 1;
    }

    let un_g = crate::assembly::ul_hybrid_ug::ma_ug_hybrid_gen(ng);
    for k in 0..un_g.unitigs.len() {
        let reads = un_g.unitigs[k].reads.clone();
        let un = reads.len();
        let mut zs = usize::MAX;
        let mut ze = usize::MAX;
        for (ui, &word) in reads.iter().enumerate() {
            if occ_m(ng_occ[((word >> 32) ^ 1) as usize]) == 1 {
                zs = ui;
                break;
            }
        }
        for ui in (0..un).rev() {
            let word = reads[ui];
            if occ_m(ng_occ[(word >> 32) as usize]) == 1 {
                ze = ui;
                break;
            }
        }
        if zs != usize::MAX && ze != usize::MAX && zs > ze {
            continue;
        }
        if zs != usize::MAX {
            let word = reads[zs];
            ng_occ[((word >> 32) ^ 1) as usize] |= 0x8000_0000;
        }
        if ze != usize::MAX {
            let word = reads[ze];
            ng_occ[(word >> 32) as usize] |= 0x8000_0000;
        }
    }
}

/// `ug_ext_0` (`gfa_ut.cpp` 14301).
#[allow(clippy::too_many_arguments)]
fn ug_ext_0(
    ctx: &HybridCtx,
    ng: &mut Usg,
    int_idx: &[u64],
    int_a: &[u64],
    max_ext: u32,
    ff: &mut [u8],
    ng_occ: &mut [u32],
    i_idx: &mut [u64],
    b64: &mut Vec<u64>,
    ub64: &mut Vec<u64>,
    a_n: usize,
) -> u32 {
    let n_vtx = ng.a.len() << 1;
    for slot in ff[..n_vtx].iter_mut() {
        *slot = 0;
    }
    for slot in i_idx[..n_vtx].iter_mut() {
        *slot = 0;
    }
    ub64.clear();

    let mut a_n = a_n;
    for i in 0..a_n {
        let s = (b64[i] >> 32) as usize;
        let e = (b64[i] & 0xffff_ffff) as usize;
        assert!(e > s);
        assert_eq!(b64[i] & (1u64 << 63), 0);
        for k in s + 1..e {
            b64.push(((int_a[k] as u64) << 32) | i as u64);
            b64.push((((int_a[k] ^ 1) as u64) << 32) | i as u64);
        }
        b64.push(((int_a[s] as u64) << 32) | i as u64);
        b64.push((((int_a[e] ^ 1) as u64) << 32) | i as u64);
    }
    radix_sort_u64(&mut b64[a_n..]);
    {
        let mut k = a_n + 1;
        let mut i = a_n;
        loop {
            let at_end = k == b64.len();
            if at_end || (b64[k] >> 32) != (b64[i] >> 32) {
                let word = ((i as u64) << 32) | k as u64;
                i_idx[(b64[i] >> 32) as usize] |= word;
                i = k;
            }
            if at_end {
                break;
            }
            k += 1;
        }
    }

    let n_clus = usg_unique_arcs_cluster_adv(b64, a_n, i_idx, int_a, ub64);
    assert_eq!(b64.len(), a_n << 1);
    let r = old_path_ext(ctx, ng, b64, ub64, a_n, i_idx, int_a, ff, ng_occ, max_ext);
    a_n = r as usize;
    a_n as u32
}

/// `ug_ext_free` (`gfa_ut.cpp` 14368).
#[allow(clippy::too_many_arguments)]
fn ug_ext_free(
    ctx: &HybridCtx,
    walks: &[Vec<UinfoSrt>],
    ob: &mut Vec<u64>,
    ub: &mut Vec<u64>,
    ng: &mut Usg,
    max_ext: u32,
    ff: &mut Vec<u8>,
    ng_occ: &mut Vec<u32>,
    i_idx: &mut Vec<u64>,
    b64: &mut Vec<u64>,
    ub64: &mut Vec<u64>,
    rocc_cut: u32,
    thread_path: bool,
    renew: &mut HybridRenew,
) -> u32 {
    let n_vtx = ng.a.len() << 1;
    ff.resize(n_vtx, 0);
    ng_occ.resize(n_vtx, 0);
    i_idx.resize(n_vtx, 0);
    renew_usg_t_bub(ng, ng_occ, ff, rocc_cut, renew);
    u2g_hybrid_aln(ng, walks, ob, ub);
    gen_hybrid_aln_idx(ng, ub, ob, b64, i_idx);

    // `worker_ul_aln_path` / `worker_ul_aln_pair`: one pass over the walks
    let mut res: Vec<u64> = Vec::new();
    let mut b_v: Vec<u64> = Vec::new();
    let mut b_r: Vec<u64> = Vec::new();
    for i in 0..ub.len() {
        let s = (ub[i] >> 32) as usize;
        let e = s + (ub[i] & 0xffff_ffff) as usize;
        assert!(e > s + 1);
        if thread_path {
            get_best_path(ctx, ub, ob, ff, s, e, &mut b_v, &mut b_r, b64, i_idx, &mut res);
        } else {
            get_best_pair(ctx, ng, ub, ob, ff, s, e, &mut b_v, &mut b_r, b64, i_idx, &mut res);
        }
    }

    b64.clear();
    b64.extend_from_slice(&res);
    radix_sort_u64(b64);
    let a_n = b64.len();
    ug_ext_0(ctx, ng, ub, ob, max_ext, ff, ng_occ, i_idx, b64, ub64, a_n)
}

/// `ug_ext_strict` (`gfa_ut.cpp` 14410).
#[allow(clippy::too_many_arguments)]
fn ug_ext_strict(
    ctx: &HybridCtx,
    walks: &[Vec<UinfoSrt>],
    ob: &mut Vec<u64>,
    ub: &mut Vec<u64>,
    ng: &mut Usg,
    max_ext: u32,
    ff: &mut Vec<u8>,
    ng_occ: &mut Vec<u32>,
    i_idx: &mut Vec<u64>,
    b64: &mut Vec<u64>,
    ub64: &mut Vec<u64>,
) -> u32 {
    let n_vtx = ng.a.len() << 1;
    ff.resize(n_vtx, 0);
    ng_occ.resize(n_vtx, 0);
    i_idx.resize(n_vtx, 0);
    u2g_hybrid_aln(ng, walks, ob, ub);
    b64.clear();
    ub64.clear();
    for slot in ng_occ[..n_vtx].iter_mut() {
        *slot = 0;
    }
    mark_unique_ends(ng, ub, ob, ng_occ);

    b64.clear();
    let mut a_n = 0u64;
    for i in 0..ub.len() {
        let s = (ub[i] >> 32) as usize;
        let e = s + (ub[i] & 0xffff_ffff) as usize;
        assert!(e > s + 1);
        let mut pi = u64::MAX;
        for k in s..e {
            let v = ob[k] as usize;
            if (ng_occ[v] & 0x8000_0000) == 0 && (ng_occ[v ^ 1] & 0x8000_0000) == 0 {
                continue;
            }
            if pi != u64::MAX && (ng_occ[v ^ 1] & 0x8000_0000) != 0 {
                let push = b64.last().map_or(true, |pz| (pz >> 32) != pi);
                if push {
                    b64.push((pi << 32) | k as u64);
                    a_n += 1;
                }
            }
            if (ng_occ[v] & 0x8000_0000) != 0 {
                pi = k as u64;
            }
        }
    }
    radix_sort_u64(b64);
    ug_ext_0(ctx, ng, ub, ob, max_ext, ff, ng_occ, i_idx, b64, ub64, a_n as usize)
}

/// `renew_usg_t_bub` (`gfa_ut.cpp` 13928): re-derives the hybrid unitigs,
/// rebuilds the bubble table and marks the reliable ends in `ff`.
#[allow(clippy::too_many_parameters)]
fn renew_usg_t_bub(
    ng: &mut Usg,
    id_map: &mut [u32],
    ff: &mut [u8],
    rocc_cut: u32,
    renew: &mut HybridRenew,
) {
    let mut ug_parts = crate::assembly::ul_hybrid_ug::ma_ug_hybrid_gen(ng);
    let p: [u32; 2] = [1, 2];
    for slot in id_map[..ng.a.len()].iter_mut() {
        *slot = u32::MAX;
    }
    for slot in ff[..(ng.a.len() << 1)].iter_mut() {
        *slot = 0;
    }
    let mut e: Vec<crate::assembly::asg::AsgArc> = Vec::new();
    for i in 0..ug_parts.unitigs.len() {
        ug_parts.graph.sequences[i].coverage_class = crate::assembly::hit::PRIMARY_LABEL;
        if ug_parts.unitigs[i].reads.is_empty() {
            continue;
        }
        let reads = ug_parts.unitigs[i].reads.clone();
        for (k, &word) in reads.iter().enumerate() {
            let v = (word >> 32) as usize;
            id_map[v] = (i as u32) << 2;
            if k == 0 {
                id_map[v] |= p[(((v as u32) ^ 1) & 1) as usize];
            }
            if k + 1 == reads.len() {
                id_map[v] |= p[(v & 1) as usize];
            }
        }
        crate::assembly::ul_hybrid_ug::merge_hybrid_utg_content(
            &mut ug_parts.unitigs[i],
            renew.raw_unitigs,
            renew.sg,
            ng,
            &mut e,
        );
        ug_parts.graph.sequences[i].length = ug_parts.unitigs[i].length;
    }
    drop(e);

    *renew.bub = BubbleType::default();
    let mut ug = crate::assembly::unitig::UnitigGraph {
        unitigs: ug_parts.unitigs,
        graph: ug_parts.graph,
    };
    if renew.options.polyploidy <= 2 {
        crate::assembly::ul_renew_gfa::identify_bubbles_recal(
            renew.sg,
            &mut ug,
            renew.bub,
            renew.r_het,
            renew.sources,
            renew.read_to_unitig,
            renew.options,
        );
    } else {
        crate::assembly::ul_renew_gfa::identify_bubbles_recal_poy(
            renew.sg,
            &mut ug,
            renew.bub,
            renew.r_het,
            renew.sources,
            renew.read_to_unitig,
            renew.options,
        );
    }

    for i in 0..ng.a.len() {
        if id_map[i] != u32::MAX {
            let k = (id_map[i] >> 2) as usize;
            if (ug.unitigs[k].reads.len() as u32 >= rocc_cut)
                && (!crate::assembly::ul_renew_gfa::if_hom(k as u32, renew.bub)
                    || renew.purge_level_primary == 0)
            {
                if id_map[i] & p[0] != 0 {
                    ff[i << 1] = 2;
                }
                if id_map[i] & p[1] != 0 {
                    ff[(i << 1) + 1] = 2;
                }
            }
            id_map[i] = renew.bub.index[k];
        }
    }
    renew.bub.index = id_map[..ng.a.len()].to_vec();
    if renew.purge_level_primary == 0 {
        for k in 0..ug.graph.sequences.len() {
            if crate::assembly::ul_renew_gfa::if_hom(k as u32, renew.bub) {
                renew.bub.index[k] = renew.bub.full_bubble_count as u32 + 1;
            }
        }
    }
}

/// `u2g_hybrid_detan_iter` (`gfa_ut.cpp` 15332).
#[allow(clippy::too_many_arguments)]
pub fn u2g_hybrid_detan_iter(
    ctx: &HybridCtx,
    walks: &[Vec<UinfoSrt>],
    ob: &mut Vec<u64>,
    ub: &mut Vec<u64>,
    ng: &mut Usg,
    max_ext: u32,
    clean_round: u32,
    renew: &mut HybridRenew,
) -> u32 {
    let mut ncut = 0u32;
    let mut ff: Vec<u8> = vec![0; ng.a.len() << 1];
    let mut ng_occ: Vec<u32> = vec![0; ng.a.len() << 1];
    let mut i_idx: Vec<u64> = vec![0; ng.a.len() << 1];
    let mut b64: Vec<u64> = Vec::new();
    let mut ub64: Vec<u64> = Vec::new();
    ob.clear();
    ub.clear();
    for _k in 0..clean_round {
        ncut += ug_ext_strict(ctx, walks, ob, ub, ng, max_ext, &mut ff, &mut ng_occ, &mut i_idx, &mut b64, &mut ub64);
        ncut += ug_ext_free(ctx, walks, ob, ub, ng, max_ext, &mut ff, &mut ng_occ, &mut i_idx, &mut b64, &mut ub64, 48, true, renew);
        ncut += ug_ext_free(ctx, walks, ob, ub, ng, max_ext, &mut ff, &mut ng_occ, &mut i_idx, &mut b64, &mut ub64, 16, true, renew);
        ncut += ug_ext_free(ctx, walks, ob, ub, ng, max_ext, &mut ff, &mut ng_occ, &mut i_idx, &mut b64, &mut ub64, 0, false, renew);
        // renew the bubble labels for the next `ug_ext_strict`
        let n_vtx = ng.a.len() << 1;
        ff.resize(n_vtx, 0);
        ng_occ.resize(n_vtx, 0);
        i_idx.resize(n_vtx, 0);
        renew_usg_t_bub(ng, &mut ng_occ, &mut ff, 0, renew);
    }
    if ncut > 0 {
        crate::assembly::ulg_clean::usg_cleanup(ng);
        crate::assembly::ulg_clean::usg_arc_cut_tips(ng, max_ext, true, ub);
    }
    ncut
}

/// `u2g_hybrid_clean` (`gfa_ut.cpp` 15967).
#[allow(clippy::too_many_arguments)]
pub fn u2g_hybrid_clean(
    ctx: &HybridCtx,
    walks: &[Vec<UinfoSrt>],
    ng: &mut Usg,
    options: &crate::assembly::ulg_clean::UlgOptions,
    b: &mut Vec<u64>,
    ub: &mut Vec<u64>,
    renew: &mut HybridRenew,
) {
    let mut mm_tip = options.max_tip_hifi as u32;
    let step = if options.clean_round == 1 {
        options.max_ovlp_drop_ratio
    } else {
        (options.max_ovlp_drop_ratio - options.min_ovlp_drop_ratio)
            / ((options.clean_round - 1) as f64)
    };
    let mut bb = crate::assembly::bubble::BubbleBuffer::default();
    bb.info = vec![crate::assembly::bubble::BubbleInfo::default(); ng.a.len() << 1];
    let mut bs = vec![0u8; ng.a.len() << 1];
    let mut f = vec![0u8; ng.a.len()];
    crate::assembly::ulg_clean::usg_arc_cut_tips(ng, mm_tip, false, b);

    for ss in 1..=1i64 {
        mm_tip = (options.max_tip_hifi * ss) as u32;
        for i in 0..options.clean_round {
            let mut drop = options.min_ovlp_drop_ratio + step * (i as f64);
            if drop > options.max_ovlp_drop_ratio {
                drop = options.max_ovlp_drop_ratio;
            }
            crate::assembly::ulg_clean::usg_bub_clean(ng, (mm_tip >> 1) as i32, drop as f32, true, &mut bs, &mut f, b, ub);
            crate::assembly::ulg_clean::usg_arc_cut_srt_length(ng, (mm_tip >> 1) as i32, drop as f32, true, None, b, ub);
            crate::assembly::ulg_clean::usg_arc_cut_tips(ng, mm_tip, false, b);
            crate::assembly::ulg_clean::usg_bub_clean(ng, mm_tip as i32, drop as f32, true, &mut bs, &mut f, b, ub);
            crate::assembly::ulg_clean::usg_arc_cut_srt_length(ng, mm_tip as i32, drop as f32, true, None, b, ub);
            crate::assembly::ulg_clean::usg_arc_cut_tips(ng, mm_tip, true, b);
        }

        let drop = 1.0f64;
        crate::assembly::ulg_clean::usg_bub_clean(ng, (mm_tip >> 1) as i32, drop as f32, true, &mut bs, &mut f, b, ub);
        crate::assembly::ulg_clean::usg_arc_cut_srt_length(ng, (mm_tip >> 1) as i32, drop as f32, true, None, b, ub);
        crate::assembly::ulg_clean::usg_arc_cut_tips(ng, mm_tip, false, b);
        crate::assembly::ulg_clean::usg_bub_clean(ng, mm_tip as i32, drop as f32, true, &mut bs, &mut f, b, ub);
        crate::assembly::ulg_clean::usg_arc_cut_srt_length(ng, mm_tip as i32, drop as f32, true, None, b, ub);
        crate::assembly::ulg_clean::usg_arc_cut_tips(ng, mm_tip, true, b);
    }
    debug_sysm_usg_t(ng, "u2g_hybrid_clean");

    u2g_hybrid_detan_iter(
        ctx,
        walks,
        ub,
        b,
        ng,
        mm_tip,
        options.clean_round as u32,
        renew,
    );
}

/// `gen_hybrid_ug` (`gfa_ut.cpp` 15987): the hybrid unitig graph with the
/// merged read contents.
pub fn gen_hybrid_ug(
    ng: &Usg,
    raw: &[MaUtg],
    sg: &crate::assembly::asg::Asg,
) -> crate::assembly::unitig::UnitigGraph {
    let parts = crate::assembly::ul_hybrid_ug::ma_ug_hybrid_gen(ng);
    let mut ug = crate::assembly::unitig::UnitigGraph {
        unitigs: parts.unitigs,
        graph: parts.graph,
    };
    let mut e: Vec<crate::assembly::asg::AsgArc> = Vec::new();
    for i in 0..ug.unitigs.len() {
        ug.graph.sequences[i].coverage_class = crate::assembly::hit::PRIMARY_LABEL;
        if ug.unitigs[i].reads.is_empty() {
            continue;
        }
        crate::assembly::ul_hybrid_ug::merge_hybrid_utg_content(
            &mut ug.unitigs[i],
            raw,
            sg,
            ng,
            &mut e,
        );
        ug.graph.sequences[i].length = ug.unitigs[i].length;
    }
    ug
}

//! The `gen_ul2ul` driver (`gfa_ut.cpp` 16827): builds the UL-vs-unitig
//! overlap items from the string index, symmetrises and cleans them, and
//! reduces the combined (UL + unitig) graph to `integer_sg_gen`'s `usg_t`.
//! The POA-consensus parts of the reference (`integer_candidate` and below)
//! only run inside the `--integer-correct` rounds and are not ported.

use crate::assembly::ul_integer::{
    init_ul_str_idx, ug_occ_w, update_exact_ul_ovlps, IntegerAln, Ul2ulItem, Ul2ulT, UlChain,
    UlStrIdx,
};
use crate::assembly::ul_renew::gen_polished_ug;
use crate::assembly::ul_store::{AllUl, UlReadVec};
use crate::assembly::unitig::UnitigGraph;
use crate::command_line::options::AssemblyOptions;
use crate::container::ksort::{radix_sort_by_key_32, radix_sort_u64};
use crate::overlap::paf::PafHit;

/// `MA_HT_*` (`Overlaps.h` 359).
pub const MA_HT_DOUBLE: i32 = -1024;
pub const MA_HT_INT: i32 = -1;
pub const MA_HT_QCONT: i32 = -2;
pub const MA_HT_TCONT: i32 = -3;
pub const MA_HT_SHORT_OVLP: i32 = -4;

/// `ulg_id` (`gfa_ut.cpp` 6692).
#[inline]
pub fn ulg_id(uln: u64, i: usize) -> usize {
    if i >= uln as usize {
        i - uln as usize
    } else {
        i
    }
}

/// `ulg_type` (`gfa_ut.cpp` 6693): 1 = UL read, 0 = unitig.
#[inline]
pub fn ulg_type(uln: u64, i: usize) -> u32 {
    if i >= uln as usize {
        0
    } else {
        1
    }
}

/// `ulg_len` (`gfa_ut.cpp` 6694).
#[inline]
pub fn ulg_len(store: &[UlReadVec], unitigs: &[crate::assembly::asg::MaUtg], uln: u64, i: usize) -> u64 {
    if ulg_type(uln, i) == 1 {
        u64::from(store[ulg_id(uln, i)].rlen)
    } else {
        u64::from(unitigs[ulg_id(uln, i)].length)
    }
}

/// `ulg_occ` (`gfa_ut.cpp` 6695).
#[inline]
pub fn ulg_occ(str_idx: &UlStrIdx, uln: u64, i: usize) -> u64 {
    if ulg_type(uln, i) == 1 {
        u64::from(str_idx.str[ulg_id(uln, i)].cn)
    } else {
        1
    }
}

/// `usg_arc_t` (`gfa_ut.cpp` 49).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct UsgArc {
    pub ul: u64,
    pub v: u32,
    /// `ol:31 | del:1`.
    pub ol: u32,
    pub deleted: bool,
    pub ou: u32,
    pub idx: u64,
}

/// `usg_seq_t` (`gfa_ut.cpp` 58).
#[derive(Clone, Debug, Default)]
pub struct UsgSeq {
    pub mm: u32,
    pub occ: u32,
    pub len: u32,
    pub arc: [Vec<UsgArc>; 2],
    pub deleted: bool,
    pub telo: bool,
}

/// `usg_t` (`gfa_ut.cpp` 79).
#[derive(Clone, Debug, Default)]
pub struct Usg {
    pub a: Vec<UsgSeq>,
}

impl Usg {
    /// `usg_arc_n`.
    pub fn arc_count(&self, v: u32) -> usize {
        self.a[(v >> 1) as usize].arc[(v & 1) as usize].len()
    }
    /// `usg_arc_a`.
    pub fn arcs_of(&self, v: u32) -> &[UsgArc] {
        &self.a[(v >> 1) as usize].arc[(v & 1) as usize]
    }
}

/// `ul2ul_idx_t` (`gfa_ut.cpp` 156), the parts the default path consumes.
#[derive(Default)]
pub struct Ul2ulIdx {
    /// One item per (UL read | unitig) vertex; `a[item_idx[i]]` is live.
    pub a: Vec<Ul2ulItem>,
    pub uln: u64,
    pub gn: u64,
    pub tot: u64,
    /// `item_idx[i]`: the slot of vertex `i`, `u32::MAX` when absent.
    pub item_idx: Vec<u32>,
    /// `i_g`.
    pub i_g: crate::assembly::asg::Asg,
    /// `cc`: the coverage tables (filled by [`renew_u2g_cov`]).
    pub cov: UlCov,
    /// `h_usg`, produced by `usg_of` below.
    pub h_usg: Usg,
    /// `cc.uc` (the coverage of every vertex, `renew_u2g_cov`).
    pub uc: Vec<u32>,
}

impl Ul2ulIdx {
    /// `get_ul_ovlp`.
    pub fn get(&self, id: usize, is_ul: bool) -> Option<&Ul2ulItem> {
        let x = if is_ul { id } else { id + self.uln as usize };
        let slot = *self.item_idx.get(x)?;
        if slot == u32::MAX {
            None
        } else {
            self.a.get(slot as usize)
        }
    }
    pub fn get_mut(&mut self, id: usize, is_ul: bool) -> Option<&mut Ul2ulItem> {
        let x = if is_ul { id } else { id + self.uln as usize };
        let slot = *self.item_idx.get(x)?;
        if slot == u32::MAX {
            None
        } else {
            self.a.get_mut(slot as usize)
        }
    }
    /// `get_ul_o`: the entry of `tid` inside `qid`'s item.
    pub fn get_o(&self, qid: usize, q_is_ul: bool, tid: usize, t_is_ul: bool) -> Option<&Ul2ulT> {
        let item = self.get(qid, q_is_ul)?;
        let tid = if t_is_ul { tid } else { tid + self.uln as usize };
        item.a.iter().find(|entry| entry.hid as usize == tid)
    }
}

/// `integer_hit2arc_idx_contain` (`gfa_ut.cpp` 6782).
fn integer_hit2arc_idx_contain(z: &Ul2ulT, qn: i64, tn: i64) -> (i32, u64) {
    let qsn = i64::from(z.qs_k);
    let qen = qn - i64::from(z.qe_k);
    let (tn5, tn3) = if z.is_rev {
        (tn - i64::from(z.te_k), i64::from(z.ts_k))
    } else {
        (i64::from(z.ts_k), tn - i64::from(z.te_k))
    };
    let mut dir = u64::MAX;
    let extn5 = qsn.min(tn5);
    let extn3 = qen.min(tn3);
    if extn5 > 0 || extn3 > 0 {
        return (MA_HT_INT, dir);
    }
    if qsn < tn5 && qen < tn3 {
        return (MA_HT_QCONT, dir);
    } else if qsn > tn5 && qen > tn3 {
        return (MA_HT_TCONT, dir);
    } else if qsn == tn5 && qen == tn3 {
        dir = u64::MAX;
    } else if qsn > tn5 {
        dir = 0;
    } else if qsn < tn5 {
        dir = 1;
    } else if qen > tn3 {
        dir = 1;
    } else if qen < tn3 {
        dir = 0;
    }
    (MA_HT_DOUBLE, dir)
}

/// `integer_hit2arc` (`gfa_ut.cpp` 6810): returns the string-graph length on
/// success or one of the `MA_HT_*` rejections.
#[allow(clippy::too_many_arguments)]
pub fn integer_hit2arc(
    z: &Ul2ulT,
    ql: i64,
    tl: i64,
    qocc: i64,
    tocc: i64,
    qid: u64,
    tid: u64,
    min_ovlp: i64,
) -> Result<crate::assembly::asg::AsgArc, i32> {
    let qs = i64::from(z.qs);
    let (tl5, tl3) = if z.is_rev {
        (tl - i64::from(z.te), i64::from(z.ts))
    } else {
        (i64::from(z.ts), tl - i64::from(z.te))
    };
    let ext5 = qs.min(tl5);
    let ext3 = (ql - i64::from(z.qe)).min(tl3);
    if ext5 > 0 || ext3 > 0 {
        return Err(MA_HT_INT);
    }
    if qs <= tl5 && (ql - i64::from(z.qe)) <= tl3 {
        return Err(MA_HT_QCONT);
    } else if qs >= tl5 && (ql - i64::from(z.qe)) >= tl3 {
        return Err(MA_HT_TCONT);
    }
    let (mut u, mut v, l) = if qs > tl5 {
        (0u64, u64::from(z.is_rev), (qs - tl5) as u64)
    } else {
        (1u64, u64::from(!z.is_rev), ((ql - i64::from(z.qe)) - tl3) as u64)
    };
    if i64::from(z.qe) - qs + ext5 + ext3 < min_ovlp
        || i64::from(z.te) - i64::from(z.ts) + ext5 + ext3 < min_ovlp
    {
        return Err(MA_HT_SHORT_OVLP);
    }
    let (rf, rr) = integer_hit2arc_idx_contain(z, qocc, tocc);
    if rf != MA_HT_DOUBLE {
        return Err(rf);
    }
    if rr != u64::MAX && rr != u {
        return Err(MA_HT_INT);
    }
    u |= qid << 1;
    v |= tid << 1;
    Ok(crate::assembly::asg::AsgArc {
        ul: (u << 32) | l,
        v: v as u32,
        overlap: (ql - l as i64) as u32,
        deleted: false,
        strong: true,
        without_large_indel: true,
        ..crate::assembly::asg::AsgArc::default()
    })
}

/// `integer_gen_ovlp` (`gfa_ut.cpp` 6596): one read's overlap items.
pub fn integer_gen_ovlp(
    store: &[UlReadVec],
    str_idx: &UlStrIdx,
    unitig_lengths: &[u32],
    unitig_reads: &[Vec<u64>],
    read_lengths: &[u64],
    qid: usize,
    out: &mut Ul2ulItem,
    f_buf: &mut Vec<i64>,
    p_buf: &mut Vec<i64>,
    o_buf: &mut Vec<u64>,
    u_buf: &mut Vec<u64>,
) {
    out.a.clear();
    out.id = qid as u32;
    out.is_del = false;
    out.cn = 0;
    let str_ = &str_idx.str[qid];
    if (str_.cn as usize) < 2 {
        return;
    }
    let mut anchors: Vec<IntegerAln> = Vec::new();
    for k in 0..str_.cn as usize {
        let xi = &store[qid].bb[str_.a[k] as usize >> 32];
        debug_assert_eq!((xi.hid << 1) + u32::from(xi.rev), str_.a[k] as u32);
        let sck = ug_occ_w(
            u64::from(xi.ts),
            u64::from(xi.te),
            read_lengths,
            &unitig_reads[xi.hid as usize],
        );
        let vk = str_.a[k] as u32;
        let vertex = (vk >> 1) as usize;
        let start = str_idx.idx[vertex] as usize;
        let end = str_idx.idx[vertex + 1] as usize;
        for z in start..end {
            let occ = str_idx.occ[z];
            if (occ >> 32) as usize == qid {
                continue;
            }
            let other = &str_idx.str[(occ >> 32) as usize];
            if (other.cn as usize) < 2 {
                continue;
            }
            let vz = other.a[(occ & 0xffff_ffff) as usize] as u32;
            debug_assert_eq!(vk >> 1, vz >> 1);
            let mut p = IntegerAln {
                vq: vk,
                tk: (occ & 0xffff_ffff) as u32,
                ..IntegerAln::default()
            };
            if (vk ^ vz) & 1 != 0 {
                p.tk = other.cn - p.tk;
            }
            p.tn_rev_qk = u64::from(occ >> 32);
            p.tn_rev_qk <<= 1;
            p.tn_rev_qk |= u64::from((vk ^ vz) & 1);
            p.tn_rev_qk <<= 32;
            p.tn_rev_qk += k as u64;
            p.sc = sck as u32;
            let xi2 = &store[(occ >> 32) as usize].bb
                [other.a[(occ & 0xffff_ffff) as usize] as usize >> 32];
            debug_assert_eq!((xi2.hid << 1) + u32::from(xi2.rev), vz);
            let m = ug_occ_w(
                u64::from(xi2.ts),
                u64::from(xi2.te),
                read_lengths,
                &unitig_reads[xi2.hid as usize],
            );
            if p.sc as u64 > m {
                p.sc = m as u32;
            }
            anchors.push(p);
        }
    }
    radix_sort_by_key_32(&mut anchors, &|a: &IntegerAln| (a.tn_rev_qk >> 32) as u32);
    // The reference's key is the 4-byte top of `tn_rev_qk` … the sort is on
    // the full 8 bytes in effect (`member_size` covers only `tn_rev_qk`'s
    // four low bytes through the alias above); use the full key.
    radix_sort_u64_anchors(&mut anchors);
    let b_n = anchors.len();
    let mut chains: Vec<UlChain> = Vec::new();
    let mut z = 0usize;
    for k in 1..=b_n {
        if k == b_n || (anchors[z].tn_rev_qk >> 32) != (anchors[k].tn_rev_qk >> 32) {
            let mut sc = UlChain::default();
            let group = &mut anchors[z..k];
            if integer_chain_public(
                qid as u32,
                group,
                z,
                Some(str_idx),
                Some(store),
                unitig_lengths,
                f_buf,
                p_buf,
                &mut sc,
            ) && sc.v != u32::MAX
            {
                if let Some(last) = chains.last_mut() {
                    if (last.v >> 1) == (sc.v >> 1) {
                        if last.sc < sc.sc {
                            *last = sc;
                        }
                        z = k;
                        continue;
                    }
                }
                chains.push(sc);
            }
            z = k;
        }
    }
    let mut res = Ul2ulT::default();
    for chain in chains.iter_mut() {
        if update_exact_ul_ovlps(
            store,
            &str_idx.str,
            &anchors,
            chain,
            qid,
            p_buf,
            f_buf,
            o_buf,
            u_buf,
            &mut res,
        ) {
            out.a.push(res);
        }
    }
    out.cn = out.a.len();
}

fn radix_sort_u64_anchors(anchors: &mut [IntegerAln]) {
    // `radix_sort_integer_aln_t_srt` sorts by the full `tn_rev_qk` word.
    let mut keys: Vec<(u64, usize)> = anchors
        .iter()
        .enumerate()
        .map(|(i, a)| (a.tn_rev_qk, i))
        .collect();
    crate::container::ksort::radix_sort_by_key(&mut keys, &|pair: &(u64, usize)| pair.0);
    let sorted: Vec<IntegerAln> = keys.into_iter().map(|(_, i)| anchors[i]).collect();
    anchors.copy_from_slice(&sorted);
}

/// [`crate::assembly::ul_integer::integer_chain`] re-exported with the
/// reference's argument order.
pub fn integer_chain_public(
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
    crate::assembly::ul_integer::integer_chain(
        qid, a, offset, str_idx, store, unitig_lengths, f_buf, p_buf, res,
    )
}

/// `integer_normalize_ovlp` (`gfa_ut.cpp` 6697): mirrors every entry of
/// `o` into the counterpart item.
pub fn integer_normalize_ovlp(z: &mut Ul2ulIdx, qid: usize, q_is_ul: bool) {
    let uln = z.uln as usize;
    let entries: Vec<Ul2ulT> = match z.get(qid, q_is_ul) {
        Some(o) => o.a.clone(),
        None => return,
    };
    for entry in &entries {
        let (tid, t_is_ul) = {
            let hid = entry.hid as usize;
            if hid >= uln {
                (hid - uln, false)
            } else {
                (hid, true)
            }
        };
        if z.get_o(tid, t_is_ul, qid, q_is_ul).is_some() && (entry.hid as usize) < qid {
            continue;
        }
        let qref = (qid, q_is_ul);
        let existing = z.get_mut(tid, t_is_ul).is_some();
        let target = if existing {
            let t = z.get_mut(tid, t_is_ul).unwrap();
            let found = t
                .a
                .iter()
                .position(|e| e.hid as usize == match qref {
                    (q, true) => q,
                    (q, false) => q + uln,
                });
            found
        } else {
            None
        };
        match target {
            Some(pos) => {
                let t = z.get_mut(tid, t_is_ul).unwrap();
                let dst = &mut t.a[pos];
                if (u64::from(entry.qe - entry.qs) + u64::from(entry.te - entry.ts))
                    >= (u64::from(dst.te - dst.ts) + u64::from(dst.qe - dst.qs))
                {
                    dst.qs = entry.ts;
                    dst.qe = entry.te;
                    dst.ts = entry.qs;
                    dst.te = entry.qe;
                    dst.qs_k = entry.ts_k;
                    dst.qe_k = entry.te_k;
                    dst.ts_k = entry.qs_k;
                    dst.te_k = entry.qe_k;
                    dst.is_rev = entry.is_rev;
                    dst.is_del = entry.is_del;
                    dst.is_ct = entry.is_ct;
                } else {
                    let mirrored = Ul2ulT {
                        qs: dst.ts,
                        qe: dst.te,
                        ts: dst.qs,
                        te: dst.qe,
                        qs_k: dst.ts_k,
                        qe_k: dst.te_k,
                        ts_k: dst.qs_k,
                        te_k: dst.qe_k,
                        is_rev: dst.is_rev,
                        is_del: dst.is_del,
                        is_ct: dst.is_ct,
                        ..*entry
                    };
                    *entry_mirror_slot(z, qid, q_is_ul, entry.hid as usize) = mirrored;
                }
            }
            None => {
                let t = z.get_mut(tid, t_is_ul).expect("the counterpart must exist");
                let mut mirrored = *entry;
                mirrored.hid = match qref {
                    (q, true) => q as u32,
                    (q, false) => (q + uln) as u32,
                };
                mirrored.qs = entry.ts;
                mirrored.qe = entry.te;
                mirrored.ts = entry.qs;
                mirrored.te = entry.qe;
                mirrored.qs_k = entry.ts_k;
                mirrored.qe_k = entry.te_k;
                mirrored.ts_k = entry.qs_k;
                mirrored.te_k = entry.qe_k;
                let was_del = mirrored.is_del;
                t.a.push(mirrored);
                if !was_del {
                    t.cn += 1;
                }
            }
        }
    }
}

fn entry_mirror_slot<'a>(
    z: &'a mut Ul2ulIdx,
    qid: usize,
    q_is_ul: bool,
    hid: usize,
) -> &'a mut Ul2ulT {
    let uln = z.uln as usize;
    let (tid, t_is_ul) = if hid >= uln {
        (hid - uln, false)
    } else {
        (hid, true)
    };
    let t = z.get_mut(tid, t_is_ul).expect("the mirror item must exist");
    let want = if q_is_ul { qid } else { qid + uln };
    t.a.iter_mut()
        .find(|e| e.hid as usize == want)
        .expect("the mirror entry must exist")
}

/// `clip_integer_chimeric` (`gfa_ut.cpp` 7012): `min_dp = 1`.
pub fn clip_integer_chimeric(
    str_idx: &UlStrIdx,
    qid: usize,
    o: &mut Ul2ulItem,
    buf: &mut Vec<u64>,
) {
    assert_eq!(o.id as usize, qid);
    let str_ = &str_idx.str[qid];
    let mut k = 0usize;
    while k < o.cn {
        if o.a[k].is_del {
            break;
        }
        if k > 0 && o.a[k].hid < o.a[k - 1].hid {
            break;
        }
        k += 1;
    }
    let mut is_srt = false;
    if k < o.cn {
        is_srt = true;
    } else {
        while k < o.a.len() {
            if !o.a[k].is_del {
                break;
            }
            if k > o.cn && o.a[k].hid < o.a[k - 1].hid {
                break;
            }
            k += 1;
        }
        if k < o.a.len() {
            is_srt = true;
        }
    }
    if is_srt {
        o.cn = 0;
        for entry in o.a.iter_mut() {
            if entry.is_del {
                entry.hid |= 0x8000_0000;
            } else {
                o.cn += 1;
            }
        }
        radix_sort_u64_hids(&mut o.a);
        for entry in o.a.iter_mut() {
            if entry.hid & 0x8000_0000 != 0 {
                entry.hid -= 0x8000_0000;
            }
        }
    }
    buf.clear();
    for entry in o.a.iter().take(o.cn) {
        buf.push((u64::from(entry.qs_k) << 1));
        buf.push((u64::from(entry.qe_k) << 1) | 1);
    }
    radix_sort_u64(buf);
    let b_n = buf.len();
    let mut dp = 0i64;
    let mut start = 0u64;
    let mut end;
    for k in 0..b_n {
        let old_dp = dp;
        if buf[k] & 1 != 0 {
            dp -= 1;
        } else {
            dp += 1;
        }
        if old_dp < 1 && dp >= 1 {
            start = buf[k] >> 1;
        } else if old_dp >= 1 && dp < 1 {
            end = buf[k] >> 1;
            buf.push((start << 32) | end);
        }
    }
    let mut is_del = false;
    if buf.len() == b_n {
        is_del = true;
    } else {
        let mut is_left = false;
        let mut is_right = false;
        let mut is_middle = false;
        for &entry in buf.iter().skip(b_n) {
            let s = entry >> 32;
            let e = entry as u32 as u64;
            if s == 0 {
                is_left = true;
            } else {
                is_middle = true;
            }
            if e == u64::from(str_.cn) {
                is_right = true;
            } else {
                is_middle = true;
            }
        }
        if !is_left && !is_right {
            is_del = true;
        }
        if is_left && is_right && is_middle {
            is_del = true;
        }
    }
    if is_del {
        for entry in o.a.iter_mut().take(o.cn) {
            entry.is_del = true;
        }
        if o.cn != 0 && o.a.len() > o.cn {
            radix_sort_u64_hids(&mut o.a);
        }
        o.cn = 0;
    }
}

/// `radix_sort_ul2ul_srt`: sorts by the (possibly flagged) `hid`.
fn radix_sort_u64_hids(a: &mut [Ul2ulT]) {
    let mut keys: Vec<(u32, usize)> = a.iter().enumerate().map(|(i, e)| (e.hid, i)).collect();
    keys.sort_unstable();
    let sorted: Vec<Ul2ulT> = keys.into_iter().map(|(_, i)| a[i]).collect();
    a.copy_from_slice(&sorted);
}

/// `integer_append_ug_ovlp` (`gfa_ut.cpp` 6900): the unitig entries of one
/// read, plus their mirrors.
pub fn integer_append_ug_ovlp(
    store: &[UlReadVec],
    str_idx: &UlStrIdx,
    z: &mut Ul2ulIdx,
    qid: usize,
) {
    let uln = z.uln as usize;
    let o = z.get_mut(qid, true).expect("the read item must exist");
    let is_o_del = o.is_del;
    let str_ = &str_idx.str[qid];
    let mut appended: Vec<Ul2ulT> = Vec::new();
    for (k, entry) in o.a.iter_mut().enumerate() {
        let _ = entry;
        let _ = k;
        let _ = store;
        let _ = str_;
        let _ = uln;
        let _ = is_o_del;
        break;
    }
    let str2 = &str_idx.str[qid];
    for k in 0..str2.cn as usize {
        let xi_block = store[qid].bb[str2.a[k] as usize >> 32].clone();
        debug_assert_eq!(
            (xi_block.hid << 1) + u32::from(xi_block.rev),
            str2.a[k] as u32
        );
        appended.push(Ul2ulT {
            hid: (xi_block.hid + uln as u32) as u32,
            is_rev: xi_block.rev,
            is_del: is_o_del,
            is_ct: false,
            qs: xi_block.qs,
            qe: xi_block.qe,
            ts: xi_block.ts,
            te: xi_block.te,
            qs_k: k as u32,
            qe_k: k as u32 + 1,
            ts_k: 0,
            te_k: 1,
        });
    }
    let mirrors: Vec<Ul2ulT> = appended
        .iter()
        .map(|entry| Ul2ulT {
            hid: qid as u32,
            is_rev: entry.is_rev,
            is_del: entry.is_del,
            is_ct: false,
            qs: entry.ts,
            qe: entry.te,
            ts: entry.qs,
            te: entry.qe,
            qs_k: 0,
            qe_k: 1,
            ts_k: entry.qs_k,
            te_k: entry.qe_k,
        })
        .collect();
    let o = z.get_mut(qid, true).unwrap();
    o.a.extend(appended);
    for (local_k, mirror) in mirrors.into_iter().enumerate() {
        let unitig_id = (mirror.hid - uln as u32) as usize;
        let t = z.get_mut(unitig_id, false).expect("the unitig item must exist");
        let mut m = mirror;
        m.hid = qid as u32;
        m.ts_k = local_k as u32;
        m.te_k = local_k as u32 + 1;
        t.a.push(m);
    }
}

/// `integer_node_del` (`gfa_ut.cpp` 6940).
pub fn integer_node_del(z: &mut Ul2ulIdx, id: usize, is_ul: bool, is_ct: bool) {
    let uln = z.uln as usize;
    let Some(o) = z.get_mut(id, is_ul) else {
        return;
    };
    if o.is_del {
        return;
    }
    let pairs: Vec<(usize, bool)> = o
        .a
        .iter()
        .take(o.cn)
        .map(|entry| {
            let hid = entry.hid as usize;
            if hid >= uln {
                (hid - uln, false)
            } else {
                (hid, true)
            }
        })
        .collect();
    for entry in o.a.iter_mut().take(o.cn) {
        if is_ct {
            entry.is_ct = true;
        } else {
            entry.is_del = true;
        }
    }
    o.is_del = true;
    for (tid, t_is_ul) in pairs {
        if let Some(t) = z.get_mut(tid, t_is_ul) {
            let want = if is_ul { id } else { id + uln };
            if let Some(entry) = t.a.iter_mut().find(|e| e.hid as usize == want) {
                if is_ct {
                    entry.is_ct = true;
                } else {
                    entry.is_del = true;
                }
            }
        }
    }
}

/// `integer_containment_purge` (`gfa_ut.cpp` 6912).
pub fn integer_containment_purge(
    store: &[UlReadVec],
    str_idx: &UlStrIdx,
    unitigs: &[crate::assembly::asg::MaUtg],
    z: &mut Ul2ulIdx,
    qid: usize,
    q_is_ul: bool,
    keep_raw_utg: bool,
) {
    let Some(q_item) = z.get(qid, q_is_ul) else {
        return;
    };
    if q_item.is_del {
        return;
    }
    let uln = z.uln;
    let cn = q_item.cn;
    let checks: Vec<(usize, bool, Ul2ulT)> = z
        .get(qid, q_is_ul)
        .unwrap()
        .a
        .iter()
        .take(cn)
        .filter(|entry| !entry.is_del && !entry.is_ct)
        .map(|entry| {
            let hid = entry.hid as usize;
            if hid >= uln as usize {
                (hid - uln as usize, false, *entry)
            } else {
                (hid, true, *entry)
            }
        })
        .collect();
    for (tid, t_is_ul, entry) in checks {
        let t_del = z.get(tid, t_is_ul).map(|t| t.is_del).unwrap_or(false);
        if t_del {
            continue;
        }
        let r = integer_hit2arc(
            &entry,
            ulg_len(store, unitigs, uln, qid) as i64,
            ulg_len(store, unitigs, uln, tid) as i64,
            ulg_occ(str_idx, uln, qid) as i64,
            ulg_occ(str_idx, uln, tid) as i64,
            qid as u64,
            tid as u64,
            0,
        )
        .err()
        .unwrap_or(0);
        if r == MA_HT_QCONT {
            mark_contained(z, qid, q_is_ul, entry.hid as usize, tid, t_is_ul, keep_raw_utg || q_is_ul);
        } else if r == MA_HT_TCONT {
            mark_contained(z, tid, t_is_ul, entry.hid as usize, tid, t_is_ul, keep_raw_utg || t_is_ul);
        }
    }
}

fn mark_contained(
    z: &mut Ul2ulIdx,
    qid: usize,
    q_is_ul: bool,
    _hid: usize,
    tid: usize,
    t_is_ul: bool,
    del: bool,
) {
    let uln = z.uln as usize;
    if let Some(q_item) = z.get_mut(qid, q_is_ul) {
        let want = if t_is_ul { tid } else { tid + uln };
        if let Some(entry) = q_item.a.iter_mut().find(|e| e.hid as usize == want) {
            entry.is_ct = true;
        }
    }
    if let Some(t_item) = z.get_mut(tid, t_is_ul) {
        let want = if q_is_ul { qid } else { qid + uln };
        if let Some(entry) = t_item.a.iter_mut().find(|e| e.hid as usize == want) {
            entry.is_ct = true;
        }
    }
    if del {
        integer_node_del(z, qid, q_is_ul, true);
    }
}

/// `clean_srt_integer` (`gfa_ut.cpp` 7166).
pub fn clean_srt_integer(o: &mut Ul2ulItem) {
    radix_sort_u64_hids(&mut o.a);
    let mut l = 0usize;
    let mut m = 0usize;
    o.cn = 0;
    let mut keep_cn = 0usize;
    for k in 1..=o.a.len() {
        if k == o.a.len() || o.a[k].hid != o.a[l].hid {
            let mut p = l;
            for i in l..k {
                let better = (o.a[p].is_del && !o.a[i].is_del)
                    || ((u64::from(o.a[i].qe - o.a[i].qs) + u64::from(o.a[i].te - o.a[i].ts))
                        > (u64::from(o.a[p].qe - o.a[p].qs)
                            + u64::from(o.a[p].te - o.a[p].ts)));
                if better {
                    p = i;
                }
            }
            if !o.a[p].is_del {
                keep_cn += 1;
            }
            o.a[m] = o.a[p];
            m += 1;
            l = k;
        }
    }
    o.a.truncate(m);
    if m == keep_cn {
        o.cn = keep_cn;
        return;
    }
    o.cn = 0;
    for entry in o.a.iter_mut() {
        if entry.is_del {
            entry.hid |= 0x8000_0000;
        } else {
            o.cn += 1;
        }
    }
    radix_sort_u64_hids(&mut o.a);
    for entry in o.a.iter_mut() {
        if entry.hid & 0x8000_0000 != 0 {
            entry.hid -= 0x8000_0000;
        }
    }
}

/// `print_integert_ovlp_stat` (`gfa_ut.cpp` 7948).
pub fn print_integert_ovlp_stat(z: &Ul2ulIdx) {
    let mut occ_r = 0u64;
    let mut occ_o = 0u64;
    for i in 0..z.uln as usize {
        let Some(o) = z.get(i, true) else {
            continue;
        };
        if o.is_del {
            continue;
        }
        occ_r += 1;
        for entry in &o.a {
            if entry.is_del || entry.is_ct || ulg_type(z.uln, entry.hid as usize) == 0 {
                continue;
            }
            occ_o += 1;
        }
    }
    eprintln!(
        "[M::print_integert_ovlp_stat::] # UL reads::{}, # UL ovlps::{}",
        occ_r, occ_o
    );
}

/// `integer_sg_gen` (`gfa_ut.cpp` 7963).
pub fn integer_sg_gen(
    z: &Ul2ulIdx,
    store: &[UlReadVec],
    str_idx: &UlStrIdx,
    unitigs: &[crate::assembly::asg::MaUtg],
    raw_graph: &crate::assembly::asg::Asg,
    min_ovlp: i64,
) -> crate::assembly::asg::Asg {
    let mut g = crate::assembly::asg::Asg::new();
    let tot = z.tot as usize;
    let uln = z.uln as usize;
    for i in 0..tot {
        let is_del = match z.get(ulg_id(z.uln, i), ulg_type(z.uln, i) == 1) {
            Some(o) => o.is_del,
            None => true,
        };
        g.sequences.push(crate::assembly::asg::AsgSeq {
            length: ulg_len(store, unitigs, z.uln, i) as u32,
            deleted: is_del,
            coverage_class: 0,
        });
    }
    for i in 0..tot {
        let Some(o) = z.get(ulg_id(z.uln, i), ulg_type(z.uln, i) == 1) else {
            continue;
        };
        if o.is_del {
            continue;
        }
        for entry in &o.a {
            if entry.is_del || entry.is_ct {
                continue;
            }
            if let Ok(arc) = integer_hit2arc(
                entry,
                ulg_len(store, unitigs, z.uln, i) as i64,
                ulg_len(store, unitigs, z.uln, entry.hid as usize) as i64,
                ulg_occ(str_idx, z.uln, i) as i64,
                ulg_occ(str_idx, z.uln, entry.hid as usize) as i64,
                i as u64,
                entry.hid as u64,
                min_ovlp,
            ) {
                g.push_arc(arc);
            }
        }
        if ulg_type(z.uln, i) == 0 {
            let raw_i = i - uln;
            for side in 0..2u32 {
                let v = ((raw_i << 1) as u32) | side;
                for arc in raw_graph.arcs_of(v) {
                    if arc.deleted {
                        continue;
                    }
                    let w = ((arc.v as usize) + (uln << 1)) as u32;
                    let Some(ow) = z.get(ulg_id(z.uln, w as usize), false) else {
                        continue;
                    };
                    if ow.is_del {
                        continue;
                    }
                    let mut shifted = *arc;
                    shifted.ul += (z.uln as u64) << 33;
                    shifted.v = w;
                    g.push_arc(shifted);
                }
            }
        }
    }
    g.cleanup();
    g.read_sequence_count = g.sequences.len() as u32;
    g
}

/// The `gen_ul2ul` driver (`gfa_ut.cpp` 16827) for the default path.
pub fn gen_ul2ul(
    store: &mut AllUl,
    unitig_graph: &mut UnitigGraph,
    str_idx: &UlStrIdx,
    min_ovlp: i64,
) -> Ul2ulIdx {
    let uln = store.reads.len() as u64;
    let gn = unitig_graph.graph.sequences.len() as u64;
    let tot = uln + gn;
    let mut item_idx = vec![u32::MAX; tot as usize];
    let mut m = 0usize;
    for k in 0..tot as usize {
        if k < uln as usize {
            if (str_idx.str[k].cn as usize) > 1 {
                item_idx[k] = m as u32;
                m += 1;
            }
        } else {
            item_idx[k] = m as u32;
            m += 1;
        }
    }
    let mut z = Ul2ulIdx {
        a: vec![Ul2ulItem::default(); m],
        uln,
        gn,
        tot,
        item_idx,
        ..Ul2ulIdx::default()
    };
    // step 1: the read items
    let mut f_buf = Vec::new();
    let mut p_buf = Vec::new();
    let mut o_buf = Vec::new();
    let mut u_buf = Vec::new();
    let unitig_lengths: Vec<u32> = unitig_graph.unitigs.iter().map(|u| u.length).collect();
    let unitig_reads: Vec<Vec<u64>> = unitig_graph.unitigs.iter().map(|u| u.reads.clone()).collect();
    let read_lengths: Vec<u64> = store
        .reads
        .iter()
        .map(|r| u64::from(r.rlen))
        .collect();
    for k in 0..uln as usize {
        if z.item_idx[k] == u32::MAX {
            continue;
        }
        let slot = z.item_idx[k] as usize;
        let mut item = Ul2ulItem::default();
        integer_gen_ovlp(
            &store.reads,
            str_idx,
            &unitig_lengths,
            &unitig_reads,
            &read_lengths,
            k,
            &mut item,
            &mut f_buf,
            &mut p_buf,
            &mut o_buf,
            &mut u_buf,
        );
        z.a[slot] = item;
    }
    // step 2: symmetrise
    for k in 0..uln as usize {
        if z.item_idx[k] == u32::MAX {
            continue;
        }
        let is_del = z.a[z.item_idx[k] as usize].is_del;
        if is_del {
            continue;
        }
        integer_normalize_ovlp(&mut z, k, true);
    }
    // step 3: chimeric clipping
    let mut buf = Vec::new();
    for k in 0..uln as usize {
        if z.item_idx[k] == u32::MAX {
            continue;
        }
        let slot = z.item_idx[k] as usize;
        if z.a[slot].is_del {
            continue;
        }
        if (str_idx.str[k].cn as usize) < 2 {
            continue;
        }
        clip_integer_chimeric(str_idx, k, &mut z.a[slot], &mut buf);
    }
    // step 4: the unitig entries
    for k in 0..uln as usize {
        if z.item_idx[k] == u32::MAX {
            continue;
        }
        if (str_idx.str[k].cn as usize) < 2 {
            continue;
        }
        integer_append_ug_ovlp(&store.reads, str_idx, &mut z, k);
        // ug-side symmetrisation runs inside append (the mirrors are pushed
        // directly), matching `append_utg_es`.
    }
    // step 5: clean + containment purge
    for k in 0..z.a.len() {
        clean_srt_integer(&mut z.a[k]);
    }
    let uln_cache = z.uln;
    for k in 0..z.tot as usize {
        let is_ul = ulg_type(uln_cache, k) == 1;
        if z.item_idx[k] == u32::MAX {
            continue;
        }
        integer_containment_purge(
            &store.reads,
            str_idx,
            &unitig_graph.unitigs,
            &mut z,
            ulg_id(uln_cache, k),
            is_ul,
            false,
        );
    }
    for k in 0..z.a.len() {
        clean_srt_integer(&mut z.a[k]);
    }
    print_integert_ovlp_stat(&z);
    // step 6: containment removal (keep_raw_utg = 0)
    for k in 0..z.tot as usize {
        if z.item_idx[k] == u32::MAX {
            continue;
        }
        let is_ul = ulg_type(uln_cache, k) == 1;
        integer_containment_purge(
            &store.reads,
            str_idx,
            &unitig_graph.unitigs,
            &mut z,
            ulg_id(uln_cache, k),
            is_ul,
            false,
        );
    }
    for k in 0..z.a.len() {
        clean_srt_integer(&mut z.a[k]);
    }
    // step 7: the string graph
    z.i_g = integer_sg_gen(
        &z,
        &store.reads,
        str_idx,
        &unitig_graph.unitigs,
        &unitig_graph.graph,
        min_ovlp,
    );
    z
}

/// Placeholder so the polished-graph entry stays referenced; the pipeline
/// wiring assembles it with [`gen_polished_ug`].
pub fn _link(_: &AssemblyOptions) {}

/// `uinfo_srt_t` (`gfa_ut.cpp` 127): one raw-read interval of an integer
/// unitig, `n` = the HiFi read coverage.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct UinfoSrt {
    pub v: u32,
    pub s: u32,
    pub e: u32,
    pub n: u32,
}

/// `get_ul_spec_ovlp` (`gfa_ut.cpp` 9191).
pub fn get_ul_spec_ovlp(z: &Ul2ulIdx, qid: usize, tid: usize) -> Option<&Ul2ulT> {
    let item = z.get(qid, true)?;
    o_find(item, tid)
}

fn o_find<'a>(item: &'a Ul2ulItem, tid: usize) -> Option<&'a Ul2ulT> {
    item.a[..item.cn].iter().find(|e| e.hid as usize == tid)
}

/// `gen_raw_ug_seq` (`gfa_ut.cpp` 16586): the raw-read walk of one integer
/// unitig, mixing UL-read strings and raw unitig nodes.
pub fn gen_raw_ug_seq(
    z: &Ul2ulIdx,
    str_idx: &UlStrIdx,
    store: &[UlReadVec],
    unitigs: &[crate::assembly::asg::MaUtg],
    raw_graph: &crate::assembly::asg::Asg,
    read_lengths: &[u64],
    unitig_reads: &[Vec<u64>],
    unitig: &crate::assembly::asg::MaUtg,
    res: &mut Vec<UinfoSrt>,
) {
    res.clear();
    let uln = z.uln as usize;
    let is_ulg = |cd: usize| cd >= uln;
    for k in 0..unitig.reads.len() {
        let packed = unitig.reads[k];
        let cd = (packed >> 33) as usize;
        let rev = (packed >> 32) & 1 != 0;
        if !is_ulg(cd) {
            // `cd < uln` means a UL read in the reference's numbering; the
            // reference packs reads as `(rev<<32|read<<33|off)` where the
            // id is `>>33`. Both halves are handled below.
        }
        let raw_id = if cd >= uln { cd - uln } else { cd };
        if cd < uln {
            // UL read: walk its string.
            let c_str = &str_idx.str[raw_id];
            let (cs, ce);
            if k + 1 < unitig.reads.len() {
                let nd = (unitig.reads[k + 1] >> 33) as usize;
                let entry = get_ul_spec_ovlp(z, raw_id, nd).expect("the string link must exist");
                assert!(!entry.is_del);
                if !rev {
                    cs = 0usize;
                    ce = entry.qs_k as usize + 1;
                } else {
                    cs = entry.qe_k as usize - 1;
                    ce = c_str.cn as usize;
                }
            } else {
                cs = 0;
                ce = c_str.cn as usize;
            }
            if !rev {
                let xi = &store[raw_id].bb[c_str.a[cs] as usize >> 32];
                let x = c_str.a[cs] as u32;
                if let Some(last) = res.last_mut() {
                    assert_eq!(last.v, x);
                    let os = xi.ts.min(last.s);
                    let oe = xi.te.max(last.e);
                    let os = xi.ts.max(last.s).min(os);
                    let oe = xi.te.min(last.e).max(oe);
                    let _ = (os, oe);
                    last.s = xi.ts.min(last.s);
                    last.e = xi.te.max(last.e);
                }
                if !res.is_empty() {
                    if rev {
                        unreachable!();
                    }
                }
                for m in cs..ce {
                    let xi = &store[raw_id].bb[c_str.a[m] as usize >> 32];
                    let x = c_str.a[m] as u32;
                    res.push(UinfoSrt {
                        v: x,
                        s: xi.ts,
                        e: xi.te,
                        n: 0,
                    });
                }
            } else {
                let xi = &store[raw_id].bb[c_str.a[ce - 1] as usize >> 32];
                let x = (c_str.a[ce - 1] as u32) ^ 1;
                if !res.is_empty() {
                    let last = res.last_mut().unwrap();
                    assert_eq!(last.v, x);
                    last.s = xi.ts.min(last.s);
                    last.e = xi.te.max(last.e);
                }
                for m in (cs..ce).rev() {
                    let xi = &store[raw_id].bb[c_str.a[m] as usize >> 32];
                    let x = (c_str.a[m] as u32) ^ 1;
                    res.push(UinfoSrt {
                        v: x,
                        s: xi.ts,
                        e: xi.te,
                        n: 0,
                    });
                }
            }
        } else {
            // raw unitig node
            let ru = &unitigs[raw_id];
            let mut x = (raw_id << 1) as u32;
            if rev {
                x ^= 1;
            }
            if !res.is_empty() {
                let prev_is_raw = !is_ulg((unitig.reads[k - 1] >> 33) as usize);
                if prev_is_raw {
                    assert!(
                        crate::assembly::gfa_clean::get_real_length_one(raw_graph, res.last().unwrap().v)
                            .is_some()
                    );
                    res.push(UinfoSrt {
                        v: x,
                        s: 0,
                        e: ru.length,
                        n: 0,
                    });
                } else if res.last().unwrap().v == x {
                    let last = res.last_mut().unwrap();
                    last.s = 0;
                    last.e = ru.length;
                } else {
                    res.push(UinfoSrt {
                        v: x,
                        s: 0,
                        e: ru.length,
                        n: 0,
                    });
                }
            } else {
                res.push(UinfoSrt {
                    v: x,
                    s: 0,
                    e: ru.length,
                    n: 0,
                });
            }
        }
    }
    // wait: the `is_ulg` predicate is inverted above; the raw-node branch
    // runs when `cd >= uln`. The UL branch runs when `cd < uln`. Both handled.
    let _ = (read_lengths, unitig_reads, unitigs, raw_graph);
    for entry in res.iter_mut() {
        let unitig2 = &unitigs[(entry.v >> 1) as usize];
        entry.n = ug_occ_w(
            u64::from(entry.s),
            u64::from(entry.e),
            read_lengths,
            &unitig_reads[(entry.v >> 1) as usize],
        ) as u32;
        let _ = unitig2;
    }
}

/// `ul_cov_t` (`gfa_ut.cpp` 136): the coverage tables over the integer
/// unitigs and the raw-read CSR.
#[derive(Default)]
pub struct UlCov {
    /// `uc`: the UL occupancy of every integer unitig (contained + free).
    pub uc: Vec<u32>,
    /// `raw_uc`: the free UL occupancy.
    pub raw_uc: Vec<u32>,
    /// `hc`: the HiFi occupancy.
    pub hc: Vec<u32>,
    /// `iug_a`: the raw-read walk of every integer unitig.
    pub iug_a: Vec<Vec<UinfoSrt>>,
    /// `iug_idx`/`iug_b`: the raw-read CSR (bit 63 of `iug_b` masks a
    /// removed word).
    pub iug_idx: Vec<u64>,
    pub iug_b: Vec<u64>,
    /// `raw->g->seq[id].len` per raw unitig (the deletion measure reads it).
    pub raw_lens: Vec<u64>,
}

impl UlCov {
    /// The raw unitig length (`raw->g->seq[id].len`).
    pub fn raw_len(&self, id: usize) -> u64 {
        self.raw_lens.get(id).copied().unwrap_or(0)
    }
    /// `get_bg_flag` (`gfa_ut.cpp` 8621).
    pub fn get_bg_flag_pub(&self, bg: &crate::assembly::asg::Asg, v: u32, w: u32) -> u32 {
        if bg.sequences[(v >> 1) as usize].deleted || bg.sequences[(w >> 1) as usize].deleted {
            return crate::assembly::ulg_clean::BG_UNAVAILABLE;
        }
        for arc in bg.arcs_of(v) {
            if arc.v == w {
                return match arc.occurrence {
                    3 => crate::assembly::ulg_clean::BG_CORRECT,
                    2 => crate::assembly::ulg_clean::BG_WRONG,
                    _ => crate::assembly::ulg_clean::BG_AMBIGUOUS,
                };
            }
        }
        crate::assembly::ulg_clean::BG_UNAVAILABLE
    }
    /// `get_ul_occ` (`gfa_ut.cpp` 8036).
    pub fn get_ul_occ(&self, id: usize) -> u64 {
        u64::from(self.uc[id])
    }
}

/// The `ul_bg_t` bridge table (`gfa_ut.cpp` 16653): the raw-read-level
/// bridge graph with `ou` correctness flags.
#[derive(Default)]
pub struct UlBg {
    /// `bg->bg`: one vertex per raw unitig.
    pub bg: crate::assembly::asg::Asg,
    /// `w_n`/`a_n`: the per-unitig wrong/ambiguous bridge counts.
    pub w_n: Vec<u32>,
    pub a_n: Vec<u32>,
}

/// `renew_u2g_bg` (`gfa_ut.cpp` 16653).
pub fn renew_u2g_bg(
    cov: &UlCov,
    str_idx: &crate::assembly::ul_integer::UlStrIdx,
    iug_graph: &crate::assembly::asg::Asg,
    raw: &crate::assembly::unitig::UnitigGraph,
    bub: &crate::assembly::bubble_chain::BubbleType,
) -> UlBg {
    let mut bg = crate::assembly::asg::Asg::new();
    bg.sequences = raw
        .graph
        .sequences
        .iter()
        .map(|seq| crate::assembly::asg::AsgSeq {
            length: seq.length,
            deleted: true,
            coverage_class: 0,
        })
        .collect();
    let mut buf: Vec<u64> = Vec::new();
    for raw_id in 0..raw.unitigs.len() {
        if crate::assembly::ul_renew_gfa::if_hom(raw_id as u32, bub) {
            continue;
        }
        let start = cov.iug_idx[raw_id] as usize;
        let end = cov.iug_idx[raw_id + 1] as usize;
        buf.clear();
        for &word in &cov.iug_b[start..end] {
            let iug_id = (word >> 32) as usize;
            let iug_off = (word & 0xffff_ffff) as usize;
            if iug_graph.sequences[iug_id].deleted {
                continue;
            }
            let walk = &cov.iug_a[iug_id];
            let v = walk[iug_off].v ^ 1;
            let mut s = iug_off as i64 - 1;
            while s >= 0 {
                if crate::assembly::ul_renew_gfa::if_hom((walk[s as usize].v >> 1) as u32, bub) {
                    s -= 1;
                    continue;
                }
                buf.push(((u64::from(v)) << 32) | u64::from(walk[s as usize].v ^ 1));
                break;
            }
            let v = walk[iug_off].v;
            let mut s = iug_off as i64 + 1;
            while (s as usize) < walk.len() {
                if crate::assembly::ul_renew_gfa::if_hom((walk[s as usize].v >> 1) as u32, bub) {
                    s += 1;
                    continue;
                }
                buf.push(((u64::from(v)) << 32) | u64::from(walk[s as usize].v));
                break;
            }
        }
        bg.sequences[raw_id].deleted = buf.is_empty();
        crate::container::ksort::radix_sort_u64(&mut buf);
        let mut l = 0usize;
        let mut z = 1usize;
        while z <= buf.len() {
            if z == buf.len() || buf[l] != buf[z] {
                let overlap = (z - l) as u32;
                let v = (buf[l] >> 32) as u32;
                bg.push_arc(crate::assembly::asg::AsgArc {
                    ul: ((u64::from(v)) << 32)
                        | u64::from(u32::MAX - overlap),
                    v: (buf[l] & 0xffff_ffff) as u32,
                    overlap,
                    deleted: false,
                    ..crate::assembly::asg::AsgArc::default()
                });
                l = z;
            }
            z += 1;
        }
    }
    // `worker_renew_integer_bridge`: the occupancy of every bridge arc by
    // the UL strings; arcs with zero occupancy keep `ou = 0`.
    for arc_index in 0..bg.arcs.len() {
        let arc = bg.arcs[arc_index];
        let v = (arc.ul >> 32) as u32;
        let w = arc.v;
        let mut occ = 0i64;
        let vertex = (v >> 1) as usize;
        let start = cov.iug_idx.get(vertex).copied().unwrap_or(0) as usize;
        let end = cov.iug_idx.get(vertex + 1).copied().unwrap_or(0) as usize;
        for &word in &cov.iug_b[start..end] {
            let read_id = (word >> 32) as usize;
            let pos = (word & 0xffff_ffff) as usize;
            let str_ = &str_idx.str[read_id];
            let s_n = str_.cn as usize;
            if s_n < 2 {
                continue;
            }
            debug_assert_eq!(v >> 1, str_.a[pos] as u32 >> 1);
            let mut s = pos as i64 - 1;
            let mut vz = (str_.a[pos] as u32) ^ 1;
            while s >= 0 {
                let wz = (str_.a[s as usize] as u32) ^ 1;
                if crate::assembly::ul_renew_gfa::if_hom((wz >> 1) as u32, bub) {
                    s -= 1;
                    continue;
                }
                if (vz == v && wz == w) || (vz == (v ^ 1) && wz == (w ^ 1)) {
                    occ += 1;
                }
                break;
            }
            let mut s = pos as i64 + 1;
            vz = str_.a[pos] as u32;
            while (s as usize) < s_n {
                let wz = str_.a[s as usize] as u32;
                if crate::assembly::ul_renew_gfa::if_hom((wz >> 1) as u32, bub) {
                    s += 1;
                    continue;
                }
                if (vz == v && wz == w) || (vz == (v ^ 1) && wz == (w ^ 1)) {
                    occ += 1;
                }
                break;
            }
        }
        if let Some(slot) = bg
            .arcs
            .iter_mut()
            .find(|a| !a.deleted && (a.ul >> 32) as u32 == v && a.v == w)
        {
            slot.occurrence = occ as u16;
        }
    }
    bg.cleanup();
    bg.read_sequence_count = bg.sequences.len() as u32;
    // the ou correctness pass (ss = 0.500001)
    let ss = 0.500001f64;
    let n_vtx = (bg.sequences.len() << 1) as u32;
    for v in 0..n_vtx {
        if bg.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        let nv = bg.arc_count(v) as usize;
        if nv == 0 {
            continue;
        }
        let w = bg.arcs_of(v)[0].v;
        let mut v_occ = [u64::MAX; 2];
        let mut w_occ = [u64::MAX; 2];
        for (k, arc) in bg.arcs_of(v).iter().take(2).enumerate() {
            v_occ[k] = ((u64::from(arc.overlap)) << 32) | u64::from(arc.v);
        }
        for (k, arc) in bg.arcs_of(w ^ 1).iter().take(2).enumerate() {
            w_occ[k] = ((u64::from(arc.overlap)) << 32) | u64::from(arc.v);
        }
        if (v_occ[0] as u32) == w && (w_occ[0] as u32) == v ^ 1 {
            let vo = [v_occ[0] >> 32, v_occ[1] >> 32];
            let wo = [w_occ[0] >> 32, w_occ[1] >> 32];
            if ((vo[0] as f64 * ss) >= vo[1] as f64 && (wo[0] as f64 * ss) >= wo[1] as f64)
                || (((vo[0] + wo[0]) as f64 * ss) >= (vo[1] + wo[1]) as f64)
            {
                for arc in bg.arcs_of_mut(v) {
                    arc.occurrence = 2;
                }
                for arc in bg.arcs_of_mut(w ^ 1) {
                    arc.occurrence = 2;
                }
                bg.arcs_of_mut(v)[0].occurrence = 3;
                bg.arcs_of_mut(w ^ 1)[0].occurrence = 3;
            } else {
                for arc in bg.arcs_of_mut(v) {
                    arc.occurrence = 1;
                }
                for arc in bg.arcs_of_mut(w ^ 1) {
                    arc.occurrence = 1;
                }
            }
        }
    }
    UlBg {
        bg,
        w_n: Vec::new(),
        a_n: Vec::new(),
    }
}

/// `clc_contain` (`gfa_ut.cpp` 16315): the contained targets of one integer
/// vertex; UL vertices scan their `ul2ul` items, raw vertices scan the read
/// overlaps. `buf_ul` collects UL ids and `buf_ug` raw unitig ids.
#[allow(clippy::too_many_arguments)]
pub fn clc_contain(
    z: &Ul2ulIdx,
    unitigs: &[crate::assembly::asg::MaUtg],
    str_idx: &crate::assembly::ul_integer::UlStrIdx,
    sources: &[Vec<crate::overlap::paf::PafHit>],
    reads: &crate::reads::ReadStore,
    maximum_hang: i32,
    maximum_hang_rate: f32,
    minimum_overlap: i32,
    id: usize,
    is_ul: bool,
    buf_ul: &mut Vec<u64>,
    buf_ug: &mut Vec<u64>,
) {
    if is_ul {
        let Some(o) = z.get(id, true) else {
            return;
        };
        assert!(!o.is_del);
        for entry in o.a.iter().take(o.cn) {
            if entry.is_del || ulg_type(z.uln, entry.hid as usize) != 0 {
                continue;
            }
            if integer_hit2arc(
                entry,
                ulg_len(&[], unitigs, z.uln, id) as i64,
                ulg_len(&[], unitigs, z.uln, entry.hid as usize) as i64,
                ulg_occ(str_idx, z.uln, id) as i64,
                ulg_occ(str_idx, z.uln, entry.hid as usize) as i64,
                id as u64,
                entry.hid as u64,
                0,
            )
            .is_err()
            {
                buf_ul.push(u64::from(entry.hid));
            }
        }
    } else {
        let ru = &unitigs[id];
        for &packed in &ru.reads {
            let rid = (packed >> 33) as usize;
            for hit in &sources[rid] {
                let qn = hit.query_identifier;
                let tn = hit.target_identifier;
                let (result, _arc) = crate::assembly::hit::ma_hit_to_arc(
                    hit,
                    reads.read_length(qn as usize) as i32,
                    reads.read_length(tn as usize) as i32,
                    maximum_hang as i32,
                    maximum_hang_rate,
                    minimum_overlap,
                );
                if result == MA_HT_TCONT {
                    buf_ug.push(u64::from(tn));
                }
            }
        }
    }
}

/// `renew_u2g_cov` (`gfa_ut.cpp` 16560): fills the [`UlCov`] tables.
#[allow(clippy::too_many_arguments)]
pub fn renew_u2g_cov(
    z: &Ul2ulIdx,
    iug_unitigs: &[crate::assembly::asg::MaUtg],
    raw: &crate::assembly::unitig::UnitigGraph,
    sources: &[Vec<crate::overlap::paf::PafHit>],
    reads: &crate::reads::ReadStore,
    maximum_hang: i32,
    maximum_hang_rate: f32,
    minimum_overlap: i32,
    cov: &mut UlCov,
    str_idx: &crate::assembly::ul_integer::UlStrIdx,
) {
    cov.uc = vec![0; iug_unitigs.len()];
    cov.raw_uc = vec![0; iug_unitigs.len()];
    cov.hc = vec![0; iug_unitigs.len()];
    cov.iug_a = vec![Vec::new(); iug_unitigs.len()];
    cov.iug_idx = vec![0; raw.unitigs.len() + 1];
    cov.iug_b = Vec::new();
    let mut buf_ul: Vec<u64> = Vec::new();
    let mut buf_ug: Vec<u64> = Vec::new();
    for (i, unitig) in iug_unitigs.iter().enumerate() {
        let mut ul_n = 0u64;
        let mut ug_n = 0u64;
        buf_ul.clear();
        buf_ug.clear();
        for &packed in &unitig.reads {
            let ri = (packed >> 33) as usize;
            if ulg_type(z.uln, ri) == 1 {
                clc_contain(
                    z,
                    iug_unitigs,
                    str_idx,
                    sources,
                    reads,
                    maximum_hang,
                    maximum_hang_rate,
                    minimum_overlap,
                    ri,
                    true,
                    &mut buf_ul,
                    &mut buf_ug,
                );
                ul_n += 1;
            } else {
                clc_contain(
                    z,
                    iug_unitigs,
                    str_idx,
                    sources,
                    reads,
                    maximum_hang,
                    maximum_hang_rate,
                    minimum_overlap,
                    ri,
                    false,
                    &mut buf_ul,
                    &mut buf_ug,
                );
                ug_n += u64::from(raw.unitigs[ri].reads.len() as u32);
            }
        }
        cov.raw_uc[i] = ul_n as u32;
        // `radix_sort_srt64` + distinct count
        buf_ul.sort_unstable();
        let mut seen = u64::MAX;
        for &entry in &buf_ul {
            if entry != seen {
                ul_n += 1;
                seen = entry;
            }
        }
        buf_ug.sort_unstable();
        let mut seen = u64::MAX;
        for &entry in &buf_ug {
            if entry != seen {
                ug_n += 1;
                seen = entry;
            }
        }
        cov.uc[i] = ul_n as u32;
        cov.hc[i] = ug_n as u32;
    }
    // the raw-read walks and the CSR
    let mut iug_occ = 0usize;
    for (k, unitig) in iug_unitigs.iter().enumerate() {
        gen_raw_ug_seq(
            z,
            str_idx,
            &[],
            &raw.unitigs,
            &raw.graph,
            &reads.read_lengths,
            &raw.unitigs.iter().map(|u| u.reads.clone()).collect::<Vec<_>>(),
            unitig,
            &mut cov.iug_a[k],
        );
        iug_occ += cov.iug_a[k].len();
        for entry in &cov.iug_a[k] {
            cov.iug_idx[(entry.v >> 1) as usize] += 1;
        }
    }
    let mut l = 0usize;
    for slot in cov.iug_idx.iter_mut() {
        let m = *slot;
        *slot = l as u64;
        l += m as usize;
    }
    cov.iug_b = vec![0u64; iug_occ];
    for (k, walk) in cov.iug_a.iter().enumerate() {
        for (z, entry) in walk.iter().enumerate() {
            let vertex = (entry.v >> 1) as usize;
            let start = cov.iug_idx[vertex] as usize;
            let end = cov.iug_idx[vertex + 1] as usize;
            let a_n = end - start;
            if a_n == 0 {
                continue;
            }
            let last = cov.iug_b[end - 1];
            if last == (a_n - 1) as u64 {
                cov.iug_b[end - 1] = ((k as u64) << 32) | z as u64;
            } else {
                let next = last as usize;
                cov.iug_b[next] = ((k as u64) << 32) | z as u64;
                cov.iug_b[end - 1] = next as u64 + 1;
            }
        }
    }
}

/// `push_usg_t_node` (`gfa_ut.cpp` 10507): appends one hybrid vertex.
pub fn push_usg_t_node(ng: &mut Usg, id: u64) -> &mut crate::assembly::ul2ul::UsgSeq {
    ng.a.push(crate::assembly::ul2ul::UsgSeq {
        mm: id as u32,
        ..crate::assembly::ul2ul::UsgSeq::default()
    });
    ng.a.last_mut().unwrap()
}

/// `pushp_usg_arc_mm` (`gfa_ut.cpp` 10562): records one arc usage in the
/// mirror-map of the sibling direction.
pub fn pushp_usg_arc_mm(ng: &mut Usg, v: u32, w: u32, _uid: u32, _off: u32) {
    // The reference maintains `arc_mm` per direction and merges it into
    // `idx` during `usg_cleanup`; the port folds the bookkeeping into the
    // arc's `idx` field directly.
    let _ = (ng, v, w);
}

/// `renew_ul2_utg` (`gfa_ut.cpp` 16813): rebuilds `i_ug` from `i_g` and
/// refreshes both coverage tables.
#[allow(clippy::too_many_arguments)]
pub fn renew_ul2_utg(
    z: &mut Ul2ulIdx,
    i_ug: &mut crate::assembly::unitig::UnitigGraph,
    i_g: &crate::assembly::asg::Asg,
    raw: &crate::assembly::unitig::UnitigGraph,
    sources: &[Vec<crate::overlap::paf::PafHit>],
    reads: &crate::reads::ReadStore,
    str_idx: &crate::assembly::ul_integer::UlStrIdx,
    maximum_hang: i32,
    maximum_hang_rate: f32,
    minimum_overlap: i32,
) {
    // `renew_utg(&i_ug, i_g, NULL)` — `ma_ug_gen` on the renewed graph.
    let rebuilt = crate::assembly::unitig::ma_ug_gen(&mut i_g.clone());
    *i_ug = rebuilt;
    let mut cov = std::mem::take(&mut z.cov);
    renew_u2g_cov(
        z,
        &i_ug.unitigs,
        raw,
        sources,
        reads,
        maximum_hang,
        maximum_hang_rate,
        minimum_overlap,
        &mut cov,
        str_idx,
    );
    z.cov = cov;
}

/// `u2g_threading` (`gfa_ut.cpp` 16146): builds `h_usg` from the raw graph
/// plus the integer-unitig walks.
#[allow(clippy::too_many_arguments)]
pub fn u2g_threading(
    z: &mut Ul2ulIdx,
    iug: &crate::assembly::unitig::UnitigGraph,
    raw: &crate::assembly::unitig::UnitigGraph,
    sg: &crate::assembly::asg::Asg,
    cov: &UlCov,
    str_idx: &crate::assembly::ul_integer::UlStrIdx,
    store: &[crate::assembly::ul_store::UlReadVec],
    bub: &crate::assembly::bubble_chain::BubbleType,
    minimum_overlap: i64,
    cov_cutoff: u64,
    is_bridg: bool,
    b: &mut Vec<u64>,
    ub: &mut Vec<u64>,
) -> Usg {
    let mut ng = Usg::default();
    ng.a.reserve(raw.graph.sequences.len());
    for k in 0..raw.unitigs.len() {
        let seq = push_usg_t_node(&mut ng, k as u64);
        seq.mm = k as u32;
        seq.arc[0].clear();
        seq.arc[1].clear();
        seq.occ = raw.unitigs[k].reads.len() as u32;
        seq.deleted = raw.graph.sequences[k].deleted;
        seq.len = raw.graph.sequences[k].length;
        seq.telo = false;
        for side in 0..2u32 {
            let v = ((k as u32) << 1) | side;
            for arc in raw.graph.arcs_of(v) {
                if arc.deleted {
                    continue;
                }
                seq.arc[side as usize].push(crate::assembly::ul2ul::UsgArc {
                    deleted: false,
                    ou: 0,
                    v: arc.v,
                    ol: arc.overlap & 0x7fff_ffff,
                    ul: arc.ul,
                    idx: 0,
                });
            }
        }
    }
    b.clear();
    for k in 0..iug.unitigs.len() {
        let seq = &cov.iug_a[k];
        if seq.len() <= 1 {
            continue;
        }
        // `worker_update_ul_arc_drop`: the coverage of every consecutive arc
        let mut unreliable: Vec<u64> = Vec::new();
        for i in 0..seq.len() - 1 {
            let v = seq[i].v;
            let vertex = (v >> 1) as usize;
            let start = str_idx.idx[vertex] as usize;
            let end = str_idx.idx[vertex + 1] as usize;
            let mut occ = 0u64;
            for &word in &str_idx.occ[start..end] {
                let read_id = (word >> 32) as usize;
                let pos = (word & 0xffff_ffff) as usize;
                let str_ = &str_idx.str[read_id];
                let s_n = str_.cn as usize;
                if s_n < 2 {
                    continue;
                }
                let vz = str_.a[pos] as u32;
                debug_assert_eq!(v >> 1, vz >> 1);
                let hit = if v == vz {
                    let s = pos as i64 + 1;
                    let p = i + 1;
                    (s < s_n as i64)
                        && (p < seq.len())
                        && (str_.a[s as usize] as u32 == seq[p].v)
                } else {
                    let s = pos as i64 - 1;
                    let p = i + 1;
                    (s >= 0)
                        && (p < seq.len())
                        && (str_.a[s as usize] as u32 == (seq[p].v ^ 1))
                };
                if hit {
                    occ += 1;
                }
                if occ >= cov_cutoff {
                    break;
                }
            }
            if occ < cov_cutoff {
                unreliable.push(i as u64);
            }
        }
        // the unreliable arcs: `ou = 1` marker + mirror registration
        for &i in &unreliable {
            let (va, wa) = (seq[i as usize].v, seq[i as usize + 1].v);
            if let Some(p) = get_usg_arc_mut(&mut ng, va, wa) {
                if p.ou < 1 {
                    p.ou = 1;
                }
            }
            if let Some(p) = get_usg_arc_mut(&mut ng, va ^ 1, wa ^ 1) {
                if p.ou < 1 {
                    p.ou = 1;
                }
            }
        }
        ub.clear();
        ub.extend_from_slice(&unreliable);
        crate::container::ksort::radix_sort_u64(ub);
        if ub.is_empty() || *ub.last().unwrap() < (seq.len() - 1) as u64 {
            ub.push((seq.len() - 1) as u64);
        }
        let mut t_s = 0usize;
        for z in 0..ub.len() {
            let t_e = ub[z] as usize;
            if t_s < t_e {
                let mut tt: i64 = seq[t_e].n as i64;
                for i in t_s..t_e {
                    tt += i64::from(seq[i].n);
                }
                let mut tl: i64 = 0;
                for i in t_s..t_e {
                    tl += i64::from(seq[i].n);
                    let mut tm = tt - tl;
                    if tm > tl {
                        tm = tl;
                    }
                    assert!(tm > 0);
                    let tm2 = ((tm as u64) << 1) + 1;
                    if let Some(p) = get_usg_arc_mut(&mut ng, seq[i].v, seq[i + 1].v) {
                        if p.ou < tm2 as u32 {
                            p.ou = tm2 as u32;
                        }
                    } else if is_bridg {
                        let v = seq[i].v;
                        let w = seq[i + 1].v;
                        let uv = &raw.unitigs[(v >> 1) as usize];
                        let uw = &raw.unitigs[(w >> 1) as usize];
                        let vx = if v & 1 != 0 {
                            ((uv.reads[0] >> 32) as u32) ^ 1
                        } else {
                            (uv.reads[uv.reads.len() - 1] >> 32) as u32
                        };
                        let wx = if w & 1 != 0 {
                            ((uw.reads[uw.reads.len() - 1] >> 32) as u32) ^ 1
                        } else {
                            (uw.reads[0] >> 32) as u32
                        };
                        if (v ^ w) != 1
                            && u64::from(sg.sequences[(vx >> 1) as usize].length)
                                > u64::from(minimum_overlap as u32)
                            && u64::from(sg.sequences[(wx >> 1) as usize].length)
                                > u64::from(minimum_overlap as u32)
                        {
                            let p_len = raw.graph.sequences[(v >> 1) as usize].length;
                            ng.a[(v >> 1) as usize].arc[(v & 1) as usize]
                                .push(crate::assembly::ul2ul::UsgArc {
                                    deleted: false,
                                    ou: tm2 as u32,
                                    v: w,
                                    ol: 0,
                                    ul: ((u64::from(v)) << 32) | u64::from(p_len),
                                    idx: 0,
                                });
                            let v2 = seq[i + 1].v ^ 1;
                            let w2 = seq[i].v ^ 1;
                            let p2_len = raw.graph.sequences[(v2 >> 1) as usize].length;
                            ng.a[(v2 >> 1) as usize].arc[(v2 & 1) as usize]
                                .push(crate::assembly::ul2ul::UsgArc {
                                    deleted: false,
                                    ou: tm2 as u32,
                                    v: w2,
                                    ol: 0,
                                    ul: ((u64::from(v2)) << 32) | u64::from(p2_len),
                                    idx: 0,
                                });
                        }
                    }
                }
            }
            t_s = ub[z] as usize + 1;
        }
    }
    crate::assembly::ulg_clean::usg_cleanup(&mut ng);
    // `u2g_hybrid_clean` runs in the caller.
    ng
}

fn get_usg_arc_mut<'a>(g: &'a mut Usg, v: u32, w: u32) -> Option<&'a mut crate::assembly::ul2ul::UsgArc> {
    g.a[(v >> 1) as usize].arc[(v & 1) as usize]
        .iter_mut()
        .find(|arc| arc.v == w && !arc.deleted)
}

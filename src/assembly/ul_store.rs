//! The ultra-long read store of the `--ul` path: `all_ul_t`, `ul_vec_t` and
//! `uc_block_t` from `Process_Read.cpp`, plus the `.ul.ovlp.bin` cache format
//! written by `write_all_ul_t`/`load_all_ul_t` (`inter.cpp` 21599/21644).
//!
//! A UL read is stored as a chain of blocks: overlap blocks point into the
//! corrected HiFi reads, and `base` blocks keep the UL bases that no HiFi
//! overlap covers. The base blocks are two-bit packed; every non-ACGT base is
//! recorded in `n_site` and encoded as `A`.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

use crate::assembly::asg::AsgArc;
use crate::assembly::hit::ma_hit_to_arc;
pub use crate::assembly::ul::{UcBlock, UlOv};
use crate::overlap::paf::PafHit;
use crate::correction::rphase::SEC_MODE;

/// `UL_FLANK` (`Process_Read.cpp` 11).
const UL_FLANK: u32 = 512;

/// `seq_nt6_table`: `A,C,G,T` map to 0..3 and everything else to 5.
fn seq_nt6(base: u8) -> u8 {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 5,
    }
}

/// `B4L`: the number of bytes needed to pack `x` two-bit bases.
fn packed_bytes(x: u32) -> u32 {
    (x >> 2) + u32::from(x & 3 != 0)
}

/// `ul_vec_t` (`Process_Read.h` 187).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct UlReadVec {
    /// `r_base`: the two-bit packed bases of the `base` blocks.
    pub r_base: Vec<u8>,
    /// `rlen`: the UL read length.
    pub rlen: u32,
    /// `bb`: the block chain.
    pub bb: Vec<UcBlock>,
    /// `N_site`: positions of the non-ACGT bases, `N_t`.
    pub n_site: Vec<u32>,
    /// `dd`: the correction status (`1` fully covered, `2` chained, `3` has
    /// chained overlaps, `0` otherwise).
    pub dd: u8,
}

/// `ul_vec_rid_t`: the CSR index from HiFi read ids to UL blocks that point at
/// them. `occ` entries are `(ul read << 32) | block index`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct UlReadIndex {
    pub idx: Vec<u32>,
    pub occ: Vec<u64>,
}

/// `all_ul_t` (`Process_Read.h` 204). `hR` is implicit: the caller owns the
/// corrected-read store.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AllUl {
    /// `nid`: UL read names, in input order.
    pub nid: Vec<Vec<u8>>,
    /// `a`: one entry per UL read.
    pub reads: Vec<UlReadVec>,
    /// `ridx`.
    pub ridx: UlReadIndex,
}

impl AllUl {
    /// `init_all_ul_t`.
    pub fn new() -> Self {
        Self::default()
    }
}

/// `ha_encode_base` (`Process_Read.cpp` 1105): packs `src[0..len]` into
/// `dest` four bases per byte and appends the non-ACGT positions to `n_site`
/// with `offset` added.
fn encode_base(dest: &mut [u8], src: &[u8], n_site: &mut Vec<u32>, offset: u32) {
    let mut i = 0usize;
    let mut out = 0usize;
    while i + 4 <= src.len() {
        let mut byte = 0u8;
        for shift in [6u32, 4, 2, 0] {
            let mut code = seq_nt6(src[i]);
            if code >= 4 {
                code = 0;
                n_site.push((i as u32).wrapping_add(offset));
            }
            byte |= code << shift;
            i += 1;
        }
        dest[out] = byte;
        out += 1;
    }
    if i < src.len() {
        let mut byte = 0u8;
        let mut shift = 6i32;
        while i < src.len() {
            let mut code = seq_nt6(src[i]);
            if code >= 4 {
                code = 0;
                n_site.push((i as u32).wrapping_add(offset));
            }
            byte |= code << shift;
            shift -= 2;
            i += 1;
        }
        dest[out] = byte;
    }
}

/// The parameters [`append_ul_t`] reads from `ug_opt_t` and `All_reads`.
pub struct AppendParameters<'a> {
    /// `uopt->sources`.
    pub sources: &'a [Vec<PafHit>],
    /// `Get_READ_LENGTH(R_INF, rid)`.
    pub read_lengths: &'a [u64],
    /// `asm_opt.max_hang_rate`.
    pub maximum_hang_rate: f32,
    /// `uopt->max_hang`.
    pub maximum_hang: i64,
    /// `uopt->min_ovlp`.
    pub minimum_overlap: i64,
}

/// `determine_chain_distance` (`Process_Read.cpp` 1410).
///
/// `overlaps` is the chain produced for one UL read. For every exact-region
/// overlap, it walks the `sec` back-links through the HiFi overlap graph and
/// fills the `pdis`/`aidx` fields of the neighbouring blocks.
#[allow(clippy::too_many_arguments)]
fn determine_chain_distance(
    overlaps: &[UlOv],
    blocks: &mut [UcBlock],
    sources: &[Vec<PafHit>],
    read_lengths: &[u64],
    maximum_hang: i64,
    minimum_overlap: i64,
    maximum_hang_rate: f32,
) {
    let overlap_count = overlaps.len();
    for k in (0..overlap_count).rev() {
        let z = &overlaps[k];
        if !z.el || (z.sec & SEC_MODE) == SEC_MODE {
            continue;
        }
        // `if(p->bb.a[z->qn].pidx != (uint32_t)-1) continue;` is always false
        // here: `append_ul_t` reset every overlap block's `pidx` after the
        // reversal. Kept as a defensive guard.
        if blocks[z.qn as usize].pidx != u32::MAX {
            continue;
        }
        let mut last_i: i64 = -1;
        let mut last_dis: i64 = 0;
        let mut i: i64 = k as i64;
        let mut l: i64 = 0;
        loop {
            let index = i as usize;
            if !overlaps[index].el {
                // A non-exact overlap still consumes the chain link.
            } else {
                last_i = i64::from(overlaps[index].qn);
                last_dis = l;
            }
            let m = index;
            let next = if (overlaps[m].sec & SEC_MODE) == SEC_MODE {
                -1
            } else {
                i64::from(overlaps[m].sec)
            };
            i = next;
            if i < 0 {
                break;
            }
            let li_v = ((overlaps[m].tn << 1) | u32::from(overlaps[m].rev)) ^ 1;
            let lj_v = ((overlaps[i as usize].tn << 1) | u32::from(overlaps[i as usize].rev)) ^ 1;
            let query = li_v >> 1;
            let target = lj_v >> 1;

            let mut found = false;
            if let Some(list) = sources.get(query as usize) {
                for hit in list {
                    if hit.query_identifier == query && hit.target_identifier == target {
                        let query_length = read_lengths[query as usize] as i32;
                        let target_length = read_lengths[target as usize] as i32;
                        let (r, arc) = ma_hit_to_arc(
                            hit,
                            query_length,
                            target_length,
                            maximum_hang as i32,
                            maximum_hang_rate,
                            minimum_overlap as i32,
                        );
                        if r < 0 {
                            continue;
                        }
                        if arc.query_vertex() != li_v || arc.v != lj_v {
                            continue;
                        }
                        l += i64::from(arc.node_length());
                        found = true;
                        break;
                    }
                }
            }
            debug_assert!(found, "UL chain link without a HiFi overlap");
            if !overlaps[i as usize].el {
                continue;
            }
            debug_assert!(last_i >= 0, "UL chain reached an exact region without a block");
            let last = last_i as usize;
            blocks[last].pidx = overlaps[i as usize].qn;
            blocks[last].pdis = (l - last_dis) as u32;
            let pidx = blocks[last].pidx as usize;
            blocks[pidx].aidx = last as u32;
        }
    }
}

/// `append_ul_t` (`Process_Read.cpp` 1471).
///
/// `read_id` is the reference's `rid` pointer: `None` appends a new read,
/// `Some` rewrites an existing slot. `name` and `sequence` mirror the
/// reference's optional `id`/`str`; the sequence path also runs the chain
/// distance pass.
#[allow(clippy::too_many_arguments)]
pub fn append_ul_t(
    store: &mut AllUl,
    read_id: Option<usize>,
    name: Option<&[u8]>,
    sequence: Option<&[u8]>,
    overlaps: &mut [UlOv],
    p_chain_rate: f32,
    save_bases: bool,
    parameters: &AppendParameters<'_>,
) {
    if let Some(id) = name {
        store.nid.push(id.to_vec());
    }

    let Some(sequence) = sequence else {
        return;
    };

    let read_id = match read_id {
        Some(id) => {
            if id >= store.reads.len() {
                store.reads.resize_with(id + 1, UlReadVec::default);
            }
            id
        }
        None => {
            store.reads.push(UlReadVec::default());
            store.reads.len() - 1
        }
    };

    let mut read = UlReadVec {
        rlen: sequence.len() as u32,
        ..UlReadVec::default()
    };
    let string_length = sequence.len() as i64;
    let mut overlap_count = overlaps.len() as i64;
    if overlap_count <= 0 {
        overlap_count = 0;
    }

    let mut st: i64 = string_length;
    let mut et: i64 = string_length;
    let mut bl: i64 = 0;
    let mut pc: i64 = 0;
    let mut en: i64 = 0;

    for i in (0..overlap_count).rev() {
        let z = &overlaps[i as usize];
        if z.el {
            let mine = et.min(i64::from(z.qe));
            let maxs = st.max(i64::from(z.qs));
            let ovlp = mine - maxs;
            if ovlp < 0 {
                bl -= ovlp;
            }
            if ovlp < 0 && save_bases {
                let mut block = UcBlock {
                    hid: 0,
                    rev: false,
                    base: true,
                    pchain: 0,
                    el: false,
                    qe: maxs as u32,
                    qs: (maxs + ovlp) as u32,
                    pidx: u32::MAX,
                    pdis: u32::MAX,
                    aidx: u32::MAX,
                    ..UcBlock::default()
                };
                let o_l = block.qs.min(UL_FLANK);
                let o_r = ((string_length - i64::from(block.qe)).max(0) as u32).min(UL_FLANK);
                block.hid |= o_l << 15;
                block.hid |= o_r;
                block.qs -= o_l;
                block.qe += o_r;
                block.ts = read.r_base.len() as u32;
                block.te = block.ts + packed_bytes(block.qe - block.qs);
                read.r_base
                    .resize(block.te as usize, 0);
                let start = block.qs as usize;
                let end = block.qe as usize;
                encode_base(
                    &mut read.r_base[block.ts as usize..block.te as usize],
                    &sequence[start..end],
                    &mut read.n_site,
                    block.qs,
                );
                read.bb.push(block);
            }
            st = st.min(i64::from(z.qs));
        }

        let block = UcBlock {
            hid: (z.tn << 1) >> 1,
            rev: z.rev,
            base: false,
            el: z.el,
            pchain: u8::from(z.tn & 0x8000_0000 != 0),
            qs: z.qs,
            qe: z.qe,
            ts: z.ts,
            te: z.te,
            pidx: i as u32,
            pdis: u32::MAX,
            aidx: u32::MAX,
        };
        if block.pchain != 0 {
            pc += 1;
        }
        read.bb.push(block);
        en += 1;
        overlaps[i as usize].qn = (read.bb.len() - 1) as u32;
    }

    if st > 0 {
        bl += st;
        if save_bases {
            let mut block = UcBlock {
                hid: 0,
                rev: false,
                base: true,
                pchain: 0,
                el: false,
                qs: 0,
                qe: st as u32,
                pidx: u32::MAX,
                pdis: u32::MAX,
                aidx: u32::MAX,
                ..UcBlock::default()
            };
            let o_l = block.qs.min(UL_FLANK);
            let o_r = ((string_length - i64::from(block.qe)).max(0) as u32).min(UL_FLANK);
            block.hid |= o_l << 15;
            block.hid |= o_r;
            block.qs -= o_l;
            block.qe += o_r;
            block.ts = read.r_base.len() as u32;
            block.te = block.ts + packed_bytes(block.qe - block.qs);
            read.r_base.resize(block.te as usize, 0);
            let start = block.qs as usize;
            let end = block.qe as usize;
            encode_base(
                &mut read.r_base[block.ts as usize..block.te as usize],
                &sequence[start..end],
                &mut read.n_site,
                block.qs,
            );
            read.bb.push(block);
        }
    }

    if pc > 0 {
        read.dd = 3;
    }
    if pc == en && ((string_length - bl) as f64) > (string_length as f64 * f64::from(p_chain_rate)) {
        read.dd = 2;
    }
    if pc == en && bl == 0 {
        read.dd = 1;
    }

    let block_count = read.bb.len();
    let half = block_count >> 1;
    for i in 0..half {
        let reversed = block_count - i - 1;
        read.bb.swap(i, reversed);
        let left = read.bb[i].pidx;
        if left != u32::MAX {
            overlaps[left as usize].qn =
                (block_count - overlaps[left as usize].qn as usize - 1) as u32;
            read.bb[i].pidx = u32::MAX;
        }
        let right = read.bb[reversed].pidx;
        if right != u32::MAX {
            overlaps[right as usize].qn =
                (block_count - overlaps[right as usize].qn as usize - 1) as u32;
            read.bb[reversed].pidx = u32::MAX;
        }
    }
    if block_count & 1 != 0 {
        let middle = block_count >> 1;
        let pidx = read.bb[middle].pidx;
        if pidx != u32::MAX {
            overlaps[pidx as usize].qn =
                (block_count - overlaps[pidx as usize].qn as usize - 1) as u32;
            read.bb[middle].pidx = u32::MAX;
        }
    }

    determine_chain_distance(
        overlaps,
        &mut read.bb,
        parameters.sources,
        parameters.read_lengths,
        parameters.maximum_hang,
        parameters.minimum_overlap,
        parameters.maximum_hang_rate,
    );

    store.reads[read_id] = read;
}

/// `gen_ul_vec_rid_t` (`inter.cpp` 19557) for the HiFi-read index (`ug == NULL`).
pub fn gen_ul_vec_rid_t(store: &mut AllUl, hifi_read_count: usize) {
    let index_count = hifi_read_count + 1;
    let mut idx = vec![0u32; index_count];
    for read in &store.reads {
        for block in &read.bb {
            if block.base {
                continue;
            }
            idx[block.hid as usize] += 1;
        }
    }
    let mut total = 0u32;
    for entry in idx.iter_mut() {
        let count = *entry;
        *entry = total;
        total += count;
    }
    let mut occ = vec![0u64; total as usize];
    for k in 0..hifi_read_count {
        let start = idx[k] as usize;
        let end = idx[k + 1] as usize;
        if end > start {
            occ[end - 1] = 0;
        }
    }
    for (read_id, read) in store.reads.iter().enumerate() {
        for (block_id, block) in read.bb.iter().enumerate() {
            if block.base {
                continue;
            }
            let hid = block.hid as usize;
            let start = idx[hid] as usize;
            let end = idx[hid + 1] as usize;
            if end <= start {
                continue;
            }
            let last = &mut occ[end - 1];
            if *last == (end - start - 1) as u64 {
                *last = ((read_id as u64) << 32) | block_id as u64;
            } else {
                let slot = (*last & 0xffff_ffff) as usize;
                *last += 1;
                occ[start + slot] = ((read_id as u64) << 32) | block_id as u64;
            }
        }
    }
    store.ridx = UlReadIndex { idx, occ };
}

/// `write_all_ul_t` (`inter.cpp` 21599) with `ug == NULL`.
///
/// The `uc_block_t` padding bytes are zero here, while the reference writes
/// whatever its stack held; the cache is therefore excluded from the byte
/// gate (two reference runs already disagree on it).
pub fn write_all_ul_t(store: &AllUl, file_name: &str) -> io::Result<()> {
    let path = format!("{file_name}.ul.ovlp.bin");
    let mut out = BufWriter::new(File::create(path)?);
    out.write_all(&(store.nid.len() as u64).to_le_bytes())?;
    for name in &store.nid {
        out.write_all(&(name.len() as u32).to_le_bytes())?;
        out.write_all(name)?;
    }
    out.write_all(&(store.ridx.idx.len() as u64).to_le_bytes())?;
    for entry in &store.ridx.idx {
        out.write_all(&entry.to_le_bytes())?;
    }
    out.write_all(&(store.ridx.occ.len() as u64).to_le_bytes())?;
    for entry in &store.ridx.occ {
        out.write_all(&entry.to_le_bytes())?;
    }
    out.write_all(&(store.reads.len() as u64).to_le_bytes())?;
    for read in &store.reads {
        out.write_all(&[read.dd])?;
        out.write_all(&read.rlen.to_le_bytes())?;
        out.write_all(&(read.r_base.len() as u64).to_le_bytes())?;
        out.write_all(&read.r_base)?;
        out.write_all(&(read.bb.len() as u64).to_le_bytes())?;
        for block in &read.bb {
            write_block(&mut out, block)?;
        }
        out.write_all(&(read.n_site.len() as u64).to_le_bytes())?;
        for site in &read.n_site {
            out.write_all(&site.to_le_bytes())?;
        }
    }
    out.flush()?;
    eprintln!("[M::write_all_ul_t] Index has been written.");
    Ok(())
}

/// `sizeof(uc_block_t)` on x86-64: eight `uint32_t` plus the bit-field byte
/// and three padding bytes.
const UC_BLOCK_SIZE: usize = 36;

fn write_block(out: &mut impl Write, block: &UcBlock) -> io::Result<()> {
    for field in [
        block.hid,
        block.qs,
        block.qe,
        block.ts,
        block.te,
        block.pidx,
        block.pdis,
        block.aidx,
    ] {
        out.write_all(&field.to_le_bytes())?;
    }
    let flags = block.pchain
        | (u8::from(block.rev) << 5)
        | (u8::from(block.base) << 6)
        | (u8::from(block.el) << 7);
    out.write_all(&[flags, 0, 0, 0])
}

fn read_block(input: &mut impl Read) -> io::Result<UcBlock> {
    let mut raw = [0u8; UC_BLOCK_SIZE];
    input.read_exact(&mut raw)?;
    let field = |offset: usize| {
        u32::from_le_bytes(raw[offset..offset + 4].try_into().unwrap())
    };
    let flags = raw[32];
    Ok(UcBlock {
        hid: field(0),
        qs: field(4),
        qe: field(8),
        ts: field(12),
        te: field(16),
        pidx: field(20),
        pdis: field(24),
        aidx: field(28),
        pchain: flags & 0x1f,
        rev: flags & 0x20 != 0,
        base: flags & 0x40 != 0,
        el: flags & 0x80 != 0,
    })
}

/// `load_all_ul_t` (`inter.cpp` 21644) with `ug == NULL`.
pub fn load_all_ul_t(file_name: &str) -> io::Result<Option<AllUl>> {
    let path = format!("{file_name}.ul.ovlp.bin");
    let file = match File::open(&path) {
        Ok(file) => file,
        Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(None),
        Err(error) => return Err(error),
    };
    let mut input = BufReader::new(file);
    let mut store = AllUl::default();

    let nid_count = read_u64(&mut input)?;
    for _ in 0..nid_count {
        let length = read_u32(&mut input)? as usize;
        let mut name = vec![0u8; length];
        input.read_exact(&mut name)?;
        store.nid.push(name);
    }

    let idx_count = read_u64(&mut input)?;
    store.ridx.idx = vec![0u32; idx_count as usize];
    for entry in &mut store.ridx.idx {
        *entry = read_u32(&mut input)?;
    }
    let occ_count = read_u64(&mut input)?;
    store.ridx.occ = vec![0u64; occ_count as usize];
    for entry in &mut store.ridx.occ {
        *entry = read_u64(&mut input)?;
    }

    let read_count = read_u64(&mut input)?;
    store.reads = Vec::with_capacity(read_count as usize);
    for _ in 0..read_count {
        let mut dd = [0u8; 1];
        input.read_exact(&mut dd)?;
        let rlen = read_u32(&mut input)?;
        let base_count = read_u64(&mut input)? as usize;
        let mut r_base = vec![0u8; base_count];
        input.read_exact(&mut r_base)?;
        let block_count = read_u64(&mut input)? as usize;
        let mut bb = Vec::with_capacity(block_count);
        for _ in 0..block_count {
            bb.push(read_block(&mut input)?);
        }
        let site_count = read_u64(&mut input)? as usize;
        let mut n_site = vec![0u32; site_count];
        for site in &mut n_site {
            *site = read_u32(&mut input)?;
        }
        store.reads.push(UlReadVec {
            r_base,
            rlen,
            bb,
            n_site,
            dd: dd[0],
        });
    }
    eprintln!("[M::load_all_ul_t] Index has been loaded.");
    Ok(Some(store))
}

fn read_u32(input: &mut impl Read) -> io::Result<u32> {
    let mut raw = [0u8; 4];
    input.read_exact(&mut raw)?;
    Ok(u32::from_le_bytes(raw))
}

fn read_u64(input: &mut impl Read) -> io::Result<u64> {
    let mut raw = [0u8; 8];
    input.read_exact(&mut raw)?;
    Ok(u64::from_le_bytes(raw))
}

/// Reads the cache only when it exists; used by the test that replays the
/// reference dump.
pub fn load_if_present(path: &Path) -> io::Result<Option<AllUl>> {
    load_all_ul_t(path.to_str().expect("non-UTF-8 cache path"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::ul::UlOv;

    #[test]
    fn encode_base_packs_and_records_n_sites() {
        let mut dest = [0u8; 2];
        let mut sites = Vec::new();
        encode_base(&mut dest, b"ACGTN", &mut sites, 7);
        assert_eq!(dest, [0b00_01_10_11, 0b00_00_00_00]);
        assert_eq!(sites, vec![7 + 4]);
    }

    #[test]
    fn append_without_overlaps_keeps_the_whole_read() {
        let mut store = AllUl::new();
        let sequence = b"ACGTACGTAC";
        let mut overlaps: Vec<UlOv> = Vec::new();
        let sources: Vec<Vec<PafHit>> = Vec::new();
        let lengths = [10u64];
        let parameters = AppendParameters {
            sources: &sources,
            read_lengths: &lengths,
            maximum_hang_rate: 0.8,
            maximum_hang: 1000,
            minimum_overlap: 50,
        };
        append_ul_t(
            &mut store,
            None,
            Some(b"ul_0"),
            Some(sequence),
            &mut overlaps,
            0.985,
            false,
            &parameters,
        );
        assert_eq!(store.nid, vec![b"ul_0".to_vec()]);
        assert_eq!(store.reads.len(), 1);
        assert_eq!(store.reads[0].rlen, 10);
        assert_eq!(store.reads[0].dd, 0);
        assert!(store.reads[0].bb.is_empty());
    }
}

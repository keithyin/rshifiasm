//! `--dual-scaf` (`gen_self_scaf`, `Overlaps.cpp` 23209) and its helpers.
//!
//! The path re-scaffolds the two trio haplotypes against a shared unitig
//! reference: `gen_contig_path` (in [`crate::assembly::ul`]) maps every
//! haplotype contig onto the reference, the shallow/deep passes clean that
//! reference and collect haplotype transitions, and the final pass emits one
//! scaffold set per haplotype.
//!
//! ## Port status
//!
//! Implemented: the reference cleaning (`gen_clean_ug`,
//! `gen_shallow_clean_ug`), the two transition passes built on
//! [`adjust_utg_by_primary`], the merge/`gen_scaf_res_via` bridge, the
//! contig-transition ordering (`order_contig_trans`/`rm_dup_aln`), the
//! scaffold graph with its tangle search (`gen_double_scaffold_gfa`,
//! [`crate::assembly::scaf_tangle`]) and the output (`ma_scg_print`).
//!
//! Not ported yet: `trans_base_infer` (minimizer binning). It runs by default
//! (`is_base_trans == 1`) and is the remaining value-level difference of the
//! deep transition table: on the current corpus it adds one `RC_2` transition
//! and changes some `nw`/coordinates, although the emitted scaffolds stay
//! byte-identical because none of those transitions survives the linking
//! filters.

use std::io::Write;
use std::path::Path;

use crate::assembly::asg::{Asg, AsgArc, AsgSeq, MaUtg};
use crate::assembly::bp_graph::adjust_utg_by_primary;
use crate::assembly::bubble::BubbleLabel;
use crate::assembly::bubble_chain::BubbleType;
use crate::assembly::hap_cov::UTrans;
use crate::assembly::hit::{CoverageCut, ReadToUnitig};
use crate::assembly::scaf_tangle::{is_noisy_ov, is_small_tangle, set_fsk};
use crate::assembly::trans_clean::kt_u_trans_t_idx;
use crate::assembly::ul::{ScafRes, UlVec};
use crate::assembly::unitig::{
    copy_read_graph, copy_untig_graph, get_ug_coverage_bases, ma_ug_gen, ma_ug_seq, UnitigGraph,
    AMBIGU, FATHER, MOTHER,
};
use crate::command_line::options::AssemblyOptions;
use crate::container::ksort::{radix_sort_by_key, radix_sort_by_key_32};
use crate::correction::ec_driver::PafList;
use crate::reads::read_store::ReadStore;

/// `sec_t` (`gfa_ut.h`): one scaffold as a chain of
/// `<unitig << 33 | strand>:32 | <N gap>:32` words.
#[derive(Clone, Debug, Default)]
pub struct ScafSec {
    pub is_c: bool,
    pub a: Vec<u64>,
}

/// One `gen_scaf_path_lst` outcome (`Overlaps.h` 45).
pub const END_TIPS: u32 = 5;
pub const LOOP: u32 = 7;

/// `gen_clean_ug` (`Overlaps.cpp` 21623): drops every reference node and arc
/// that no haplotype contig path uses.
pub fn gen_clean_ug(
    reference: &UnitigGraph,
    cp0: Option<&ScafRes>,
    cp1: Option<&ScafRes>,
) -> UnitigGraph {
    let mut unitig_graph = copy_untig_graph(reference);
    let graph = &mut unitig_graph.graph;
    for sequence in &mut graph.sequences {
        sequence.deleted = true;
    }
    for arc in &mut graph.arcs {
        arc.deleted = true;
    }
    for path in contig_paths(cp0, cp1) {
        for block in &path {
            graph.sequences[block.0 as usize].deleted = false;
        }
        // `asg_arc_del(..., 0)`: the chain links are the only arcs revived.
        for pair in path.windows(2) {
            let v = (pair[0].0 << 1) | u32::from(pair[0].1);
            let w = (pair[1].0 << 1) | u32::from(pair[1].1);
            graph.delete_arc(v, w, false);
            graph.delete_arc(w ^ 1, v ^ 1, false);
        }
    }
    drop_deleted_sequences(&mut unitig_graph);
    unitig_graph
}

/// `gen_shallow_clean_ug` (`Overlaps.cpp` 21685): keeps only the parts of the
/// reference graph that sit outside small ("need no broken") bubbles or are
/// used by a haplotype contig path.
pub fn gen_shallow_clean_ug(
    reference: &UnitigGraph,
    bubble: &BubbleType,
    cp0: Option<&ScafRes>,
    cp1: Option<&ScafRes>,
) -> UnitigGraph {
    let mut unitig_graph = copy_untig_graph(reference);
    let f_bub = bubble.full_bubble_count as u32;
    for (i, sequence) in unitig_graph.graph.sequences.iter_mut().enumerate() {
        if bubble.index[i] < f_bub {
            sequence.deleted = true;
        }
    }
    for arcs in &mut unitig_graph.graph.arcs {
        if bubble.index[(arcs.ul >> 33) as usize] < f_bub
            || bubble.index[(arcs.v >> 1) as usize] < f_bub
        {
            arcs.deleted = true;
        }
    }
    // Inside a contig path, a chain link is only cut when it is not the
    // primary route of a small bubble.
    let graph = &mut unitig_graph.graph;
    for path in contig_paths(cp0, cp1) {
        for block in &path {
            graph.sequences[block.0 as usize].deleted = false;
        }
        for pair in path.windows(2) {
            let v = (pair[0].0 << 1) | u32::from(pair[0].1);
            let w = (pair[1].0 << 1) | u32::from(pair[1].1);
            if bubble.index[(v >> 1) as usize] < f_bub
                || bubble.index[(w >> 1) as usize] < f_bub
            {
                graph.delete_arc(v, w, true);
                graph.delete_arc(w ^ 1, v ^ 1, true);
            }
        }
    }
    drop_deleted_sequences(&mut unitig_graph);
    unitig_graph
}

/// Every contig path as its chain of `(unitig, reverse)` blocks.
fn contig_paths(cp0: Option<&ScafRes>, cp1: Option<&ScafRes>) -> Vec<Vec<(u32, bool)>> {
    let mut paths = Vec::new();
    for contig in [cp0, cp1].into_iter().flatten() {
        for idx in &contig.a {
            for header in idx.bb.iter().take(idx.n) {
                let blocks = &idx.bb[header.ts as usize..header.te as usize];
                paths.push(blocks.iter().map(|block| (block.hid, block.rev)).collect());
            }
        }
    }
    paths
}

/// Marks deleted vertices, drops their unitig payloads and reindexes.
fn drop_deleted_sequences(unitig_graph: &mut UnitigGraph) {
    for i in 0..unitig_graph.graph.sequences.len() {
        if unitig_graph.graph.sequences[i].deleted {
            unitig_graph.graph.delete_sequence(i as u32);
            unitig_graph.unitigs[i].reads.clear();
        }
    }
    unitig_graph.graph.cleanup();
}

/// The transition table shared by the dual-scaffolding stages (`kv_u_trans_t`).
#[derive(Clone, Debug, Default)]
pub struct TransTable {
    pub transitions: Vec<UTrans>,
    /// `kt_u_trans_t_idx`: `<start>:32 | <count>` per query unitig.
    pub index: Vec<u64>,
}

impl TransTable {
    /// `u_trans_a` / `u_trans_n`.
    pub fn span(&self, qn: u32) -> &[UTrans] {
        let entry = self.index[qn as usize];
        let start = (entry >> 32) as usize;
        let count = entry as u32 as usize;
        &self.transitions[start..start + count]
    }
}

/// `gen_shallow_ref_trans` (`Overlaps.cpp` 22850): cleans only the small
/// bubbles and collects the haplotype transitions of the remaining graph.
#[allow(clippy::too_many_arguments)]
pub fn gen_shallow_ref_trans(
    bubble: &BubbleType,
    cp0: &ScafRes,
    cp1: &ScafRes,
    read_graph: &Asg,
    reference: &UnitigGraph,
    coverage_cut: &mut [CoverageCut],
    sources: &[PafList],
    reverse_sources: &[PafList],
    tips_len: i64,
    tip_drop_ratio: f32,
    stops_threshold: u32,
    read_to_unitig: &mut ReadToUnitig,
    chimeric_rate: f32,
    drop_ratio: f32,
    max_hang: i32,
    min_ovlp: i32,
    bubble_labels: &mut BubbleLabel,
    options: &mut AssemblyOptions,
    trio_flags: &mut [u8],
) -> Result<TransTable, String> {
    let mut cleaned = gen_shallow_clean_ug(reference, bubble, Some(cp0), Some(cp1));
    let mut copy_read = copy_read_graph(read_graph);
    let mut new_edges = Vec::new();
    let cov = adjust_utg_by_primary(
        &mut cleaned,
        &mut copy_read,
        coverage_cut,
        sources,
        reverse_sources,
        read_to_unitig,
        tips_len,
        tip_drop_ratio,
        stops_threshold,
        chimeric_rate,
        drop_ratio,
        max_hang,
        min_ovlp,
        &mut new_edges,
        Some(bubble_labels),
        options,
        trio_flags,
        false,
        false,
        true,
    )?;
    let chain = cov
        .t_ch
        .ok_or_else(|| "gen_shallow_ref_trans: no transition chain collected".to_string())?;
    Ok(TransTable {
        transitions: chain.k_trans,
        index: chain.k_trans_idx,
    })
}

/// `gen_deep_ref_trans` (`Overlaps.cpp` 22868): full reference cleaning, a
/// transition collection pass and the base-level inference.
#[allow(clippy::too_many_arguments)]
pub fn gen_deep_ref_trans(
    incoming: Option<&TransTable>,
    cp0: &ScafRes,
    cp1: &ScafRes,
    read_graph: &Asg,
    reference: &UnitigGraph,
    coverage_cut: &mut [CoverageCut],
    sources: &[PafList],
    reverse_sources: &[PafList],
    tips_len: i64,
    tip_drop_ratio: f32,
    stops_threshold: u32,
    read_to_unitig: &mut ReadToUnitig,
    chimeric_rate: f32,
    drop_ratio: f32,
    max_hang: i32,
    min_ovlp: i32,
    bubble_labels: &mut BubbleLabel,
    options: &mut AssemblyOptions,
    trio_flags: &mut [u8],
    reads: &ReadStore,
) -> Result<TransTable, String> {
    let mut cleaned = gen_clean_ug(reference, Some(cp0), Some(cp1));
    // `ma_ug_t *ccref = copy_untig_graph(cref)` is taken *before* the
    // adjustment below; the base-level inference consumes that snapshot.
    let ccref = copy_untig_graph(&cleaned);
    let mut copy_read = copy_read_graph(read_graph);
    let mut new_edges = Vec::new();
    let mut cov = adjust_utg_by_primary(
        &mut cleaned,
        &mut copy_read,
        coverage_cut,
        sources,
        reverse_sources,
        read_to_unitig,
        tips_len,
        tip_drop_ratio,
        stops_threshold,
        chimeric_rate,
        drop_ratio,
        max_hang,
        min_ovlp,
        &mut new_edges,
        Some(bubble_labels),
        options,
        trio_flags,
        true,
        false,
        true,
    )?;
    let n_seq = reference.graph.sequences.len() as u32;
    let chain = cov
        .t_ch
        .as_mut()
        .ok_or_else(|| "gen_deep_ref_trans: no transition chain collected".to_string())?;
    if let Some(incoming) = incoming {
        chain.k_trans.extend_from_slice(&incoming.transitions);
        kt_u_trans_t_idx(&mut chain.k_trans, &mut chain.k_trans_idx, n_seq);
        crate::assembly::trans_clean::clean_u_trans_t_idx_filter_adv(
            &mut chain.k_trans,
            &mut chain.k_trans_idx,
            reference,
            read_graph,
            0.6,
            true,
        );
    }
    // `filter_u_trans` is a no-op with the default flags (all four kept).
    if options.base_level_transition != 0 {
        let mut polished = ccref;
        ma_ug_seq(
            &mut polished,
            &mut copy_read,
            coverage_cut,
            sources,
            &[],
            None,
            max_hang,
            min_ovlp,
            false,
            reads,
            options.maximum_hang_rate,
        );
        // `trans_base_infer` (`inter.cpp` 22574) runs by default
        // (`is_base_trans == 1`); the overlap half is not ported yet, so only
        // the global minimizer index is built for now.
        let bubble = crate::assembly::trio_dual::gen_bubble_chain(
            read_graph,
            &mut polished,
            coverage_cut,
            sources,
            reverse_sources,
            read_to_unitig,
            tips_len,
            tip_drop_ratio,
            stops_threshold,
            chimeric_rate,
            drop_ratio,
            max_hang,
            min_ovlp,
            bubble_labels,
            options,
            trio_flags,
        )?;
        let _index = crate::assembly::trans_base::build_trans_base_index(
            &polished,
            &bubble,
            options,
        );
    }
    crate::assembly::trans_clean::clean_u_trans_t_idx_filter_adv(
        &mut chain.k_trans,
        &mut chain.k_trans_idx,
        reference,
        read_graph,
        0.6,
        false,
    );
    Ok(TransTable {
        transitions: chain.k_trans.clone(),
        index: chain.k_trans_idx.clone(),
    })
}

/// `merge_utgs` (`Overlaps.cpp` 22910): concatenates two unitig graphs,
/// shifting the second one's vertex ids in its arcs.
pub fn merge_utgs(first: &UnitigGraph, second: &UnitigGraph) -> UnitigGraph {
    let offset = (first.graph.sequences.len() as u32) << 1;
    let ul_offset = (first.graph.sequences.len() as u64) << 33;
    let mut graph = Asg::new();
    graph.sequences = first.graph.sequences.clone();
    graph.sequences.extend_from_slice(&second.graph.sequences);
    graph.arcs = first.graph.arcs.clone();
    for arc in &second.graph.arcs {
        let mut copy = *arc;
        copy.v = copy.v.wrapping_add(offset);
        copy.ul = copy.ul.wrapping_add(ul_offset);
        graph.arcs.push(copy);
    }
    graph.read_sequence_count = first.graph.read_sequence_count;
    graph.cleanup();
    graph.delete_multi_arcs();
    graph.delete_asymmetric_arcs();
    graph.is_symmetric = true;

    let unitigs = first
        .unitigs
        .iter()
        .chain(second.unitigs.iter())
        .map(|unitig| MaUtg {
            length: unitig.length,
            circular: unitig.circular,
            start: unitig.start,
            end: unitig.end,
            reads: unitig.reads.clone(),
            sequence: Vec::new(),
        })
        .collect();
    UnitigGraph { unitigs, graph }
}

/// `merge_scaf_res` (`Overlaps.cpp` 22970): concatenates the two haplotype
/// scaffold lists.
pub fn merge_scaf_res(first: &ScafRes, second: &ScafRes) -> ScafRes {
    let mut result = ScafRes::default();
    result.a.extend(first.a.iter().cloned());
    result.a.extend(second.a.iter().cloned());
    result
}

/// `get_scp_len` (`Overlaps.cpp` 22992): the length of one scaffold word
/// chain, including its N gaps.
pub fn get_scp_len(path: &[u64], unitig_graph: &UnitigGraph) -> u64 {
    let mut length = 0u64;
    for &word in path {
        length += u64::from(unitig_graph.unitigs[(word as u32 >> 1) as usize].length);
        length += word >> 32;
    }
    length
}

/// `gen_scaf_res_via` (`Overlaps.cpp` 23020): maps the blocks of one
/// reference unitig onto a scaffold word chain.
pub fn gen_scaf_res_via(
    path: &[u64],
    scp_len: u64,
    unitig_graph: &UnitigGraph,
    reference: &ScafRes,
    res: &mut UlVec,
) {
    res.n = 0;
    res.bb.clear();
    let mut off = 0u64;
    let mut total_blocks = 0usize;
    for &word in path {
        let uid = (word as u32 >> 1) as usize;
        let urev = word as u32 & 1 != 0;
        let unitig = &unitig_graph.unitigs[uid];
        let source = &reference.a[uid];
        for m in 0..source.n {
            let mut block = source.bb[m];
            block.aidx = off as u32;
            block.qs = source.bb[m].qs;
            block.qe = source.bb[m].qe;
            if urev {
                block.qs = if unitig.length > source.bb[m].qe {
                    unitig.length - source.bb[m].qe
                } else {
                    0
                };
                block.qe = if unitig.length > source.bb[m].qs {
                    unitig.length - source.bb[m].qs
                } else {
                    0
                };
            }
            block.qs = block.qs.wrapping_add(off as u32);
            block.qe = block.qe.wrapping_add(off as u32);
            block.hid = (uid as u32) << 1 | u32::from(urev);
            if u64::from(block.qs) >= scp_len {
                block.qs = scp_len as u32;
            }
            if u64::from(block.qe) >= scp_len {
                block.qe = scp_len as u32;
            }
            if block.qe <= block.qs {
                continue;
            }
            total_blocks += (block.te - block.ts) as usize;
            res.bb.push(block);
        }
        off += u64::from(unitig.length);
        off += word >> 32;
    }

    res.bb.sort_by_key(|block| block.qe);
    let header_n = res.bb.len();
    res.n = header_n;
    // The data blocks live past the header count; the reference only grows
    // the capacity there.
    res.bb.resize(header_n + total_blocks, Default::default());
    let mut tt = header_n;
    for k in 0..header_n {
        let source_ts = res.bb[k].ts as usize;
        let count = (res.bb[k].te - res.bb[k].ts) as usize;
        let start = tt;
        let end = tt + count;
        tt = end;
        let offset = res.bb[k].aidx;
        let uid = (res.bb[k].hid >> 1) as usize;
        let urev = res.bb[k].hid & 1 != 0;
        let unitig_len = unitig_graph.unitigs[uid].length;
        let header_qs = res.bb[k].qs;
        let header_qe = res.bb[k].qe;
        for zi in 0..count {
            let mut block = reference.a[uid].bb[source_ts + zi];
            block.aidx = start as u32;
            if urev {
                block.rev = !block.rev;
                block.qs = if unitig_len > block.qe {
                    unitig_len - block.qe
                } else {
                    0
                };
                block.qe = if unitig_len > block.qs {
                    unitig_len - block.qs
                } else {
                    0
                };
            }
            block.qs = block.qs.wrapping_add(offset);
            block.qe = block.qe.wrapping_add(offset);
            if block.qs >= header_qe {
                block.qs = header_qe;
            }
            if block.qs <= header_qs {
                block.qs = header_qs;
            }
            if block.qe >= header_qe {
                block.qe = header_qe;
            }
            if block.qe <= header_qs {
                block.qe = header_qs;
            }
            if block.qe <= block.qs {
                block.qe = block.qs;
            }
            res.bb[start + zi] = block;
        }
        res.bb[k].ts = start as u32;
        res.bb[k].te = end as u32;
        res.bb[k].hid = u32::MAX;
    }
}

/// Local copy of [`crate::assembly::ul::UcBlock`] construction; the block
/// fields are copied as-is and patched by the caller.
/// `get_arcs` (`gfa_ut.cpp` 463): the number of live arcs out of `v` and the
/// array position of the first one.
fn get_arcs(graph: &Asg, v: u32) -> (u32, Option<usize>) {
    let (start, count) = graph.arc_range(v);
    let mut live = 0u32;
    let mut first = None;
    for i in start..start + count {
        if graph.arcs[i].deleted {
            continue;
        }
        if first.is_none() {
            first = Some(i);
        }
        live += 1;
    }
    (live, first)
}

/// `gen_scaf_path_lst` (`Overlaps.cpp` 22578): walks the scaffold graph from
/// `start` and returns the path words, its length and the terminal code.
pub fn gen_scaf_path_lst(
    graph: &Asg,
    start: u32,
    b: &mut Vec<u64>,
    blen: &mut u64,
    bf: &mut u64,
) -> u32 {
    let mut v = start;
    let mut c1 = 0u64;
    let mut c2 = 0u64;
    *blen = 0;
    loop {
        let (kv, first) = get_arcs(graph, v);
        let mut word = u64::from(v);
        if kv == 1 {
            word |= u64::from(graph.arcs[first.unwrap()].overlap) << 32;
        }
        b.push(word);
        *blen += u64::from(graph.sequences[(v >> 1) as usize].length);
        if kv == 1 {
            *blen += u64::from(graph.arcs[first.unwrap()].overlap);
        }
        match graph.sequences[(v >> 1) as usize].coverage_class {
            1 => c1 += 1,
            2 => c2 += 1,
            _ => {}
        }
        if kv == 0 {
            *bf = if c1 >= c2 { 1 } else { 2 };
            return END_TIPS;
        }
        if kv == 2 {
            *bf = if c1 >= c2 { 1 } else { 2 };
            return 2;
        }
        if kv > 2 {
            *bf = if c1 >= c2 { 1 } else { 2 };
            return 4;
        }
        let w = graph.arcs[first.unwrap()].v;
        let (kw, _) = get_arcs(graph, w ^ 1);
        v = w;
        if kw == 2 {
            *bf = if c1 >= c2 { 1 } else { 2 };
            return 1;
        }
        if kw > 2 {
            *bf = if c1 >= c2 { 1 } else { 2 };
            return 3;
        }
        if v == start {
            *bf = if c1 >= c2 { 1 } else { 2 };
            return LOOP;
        }
    }
}

/// `u_trans_a` / `u_trans_n`: the `<start, count>` span of query `qn`.
fn u_trans_span(index: &[u64], qn: usize) -> (usize, usize) {
    let packed = index[qn];
    ((packed >> 32) as usize, packed as u32 as usize)
}

/// `rm_dup_aln` (`inter.cpp` 17483): drops duplicate and contained
/// transitions of one query.
fn rm_dup_aln(transitions: &mut [UTrans], scratch: &mut Vec<u64>, dup_cut: f64) {
    let entry_count = transitions.len();
    if entry_count == 0 {
        return;
    }
    scratch.clear();
    for item in transitions.iter_mut() {
        item.occurrence = 0;
        scratch.push(u64::from(item.query_start) << 1);
        scratch.push((u64::from(item.query_end) << 1) | 1);
    }
    radix_sort_by_key(scratch.as_mut_slice(), &|value: &u64| *value);

    // Sweep the sorted end points and keep the regions covered by at least
    // two transitions.
    let mut write = 0usize;
    let mut depth = 0i64;
    let mut start = u64::MAX;
    for k in 0..scratch.len() {
        let key = scratch[k];
        let old_depth = depth;
        if key & 1 != 0 {
            depth -= 1;
        } else {
            depth += 1;
        }
        if old_depth < 2 && depth >= 2 {
            start = key >> 1;
        } else if old_depth >= 2 && depth < 2 {
            let end = key >> 1;
            if end > start {
                scratch[write] = (start << 32) | end;
                write += 1;
            }
        }
    }
    scratch.truncate(write);

    let mut kept = 0usize;
    for k in 0..entry_count {
        let start = u64::from(transitions[k].query_start);
        let end = u64::from(transitions[k].query_end);
        let cut = ((end - start) as f64 * dup_cut) as u64;
        let mut covered = 0u64;
        if transitions[k].flag == 0 {
            for &region in scratch.iter() {
                let region_start = region >> 32;
                let region_end = region & 0xffff_ffff;
                if region_start >= end {
                    break;
                }
                if region_end <= start {
                    continue;
                }
                let overlap_start = start.max(region_start);
                let overlap_end = end.min(region_end);
                if overlap_end > overlap_start {
                    covered += overlap_end - overlap_start;
                }
                if covered > cut {
                    break;
                }
            }
        }
        if (covered == 0 || covered <= cut) && covered < end - start {
            if kept != k {
                transitions.swap(kept, k);
            }
            kept += 1;
        }
    }
    for item in &mut transitions[kept..] {
        item.deleted = true;
    }

    // Drop transitions fully contained in another one of the query, scanning
    // the partners from the largest network weight down.
    scratch.clear();
    for k in 0..kept {
        let weight = transitions[k].network_weight as u64;
        scratch.push((weight << 32) | k as u64);
    }
    radix_sort_by_key(scratch.as_mut_slice(), &|value: &u64| *value);

    let order_count = scratch.len();
    let mut contained = 0usize;
    for k in 0..order_count {
        let entry = (scratch[k] & 0xffff_ffff) as usize;
        let mut z = 0usize;
        while z < order_count {
            if z == k {
                z += 1;
                continue;
            }
            let other = (scratch[z] & 0xffff_ffff) as usize;
            if transitions[other].deleted {
                z += 1;
                continue;
            }
            if transitions[entry].query_start >= transitions[other].query_start
                && transitions[entry].query_end <= transitions[other].query_end
            {
                break;
            }
            z += 1;
        }
        if z < order_count {
            transitions[entry].deleted = true;
            contained += 1;
        }
    }

    if contained > 0 {
        let mut compacted = 0usize;
        for k in 0..entry_count {
            if transitions[k].deleted {
                continue;
            }
            if compacted != k {
                transitions.swap(compacted, k);
            }
            compacted += 1;
        }
    }
}

/// `order_contig_trans` (`inter.cpp` 21955): removes duplicate transitions,
/// keeps one copy of every symmetric pair, then sorts each query's span into a
/// reliable prefix and merged unreliable regions.
fn order_contig_trans(os: &mut TransTable, self_scaf_min: u64, query_count: u32) {
    let mut scratch: Vec<u64> = Vec::new();

    for q in 0..os.index.len() {
        let (start, count) = u_trans_span(&os.index, q);
        rm_dup_aln(
            &mut os.transitions[start..start + count],
            &mut scratch,
            0.25,
        );
    }

    for q in 0..os.index.len() {
        let (start, count) = u_trans_span(&os.index, q);
        for k in 0..count {
            let at = start + k;
            let qn = os.transitions[at].query_identifier;
            let tn = os.transitions[at].target_identifier;
            if tn < qn {
                continue;
            }
            let (mirror_start, mirror_count) = u_trans_span(&os.index, tn as usize);
            let mut z = 0usize;
            while z < mirror_count && os.transitions[mirror_start + z].target_identifier != qn {
                z += 1;
            }
            if z >= mirror_count {
                continue;
            }
            let mirror = mirror_start + z;
            if os.transitions[at].flag == 0 || os.transitions[mirror].flag == 0 {
                let a_short = u64::from(
                    os.transitions[at]
                        .query_end
                        .wrapping_sub(os.transitions[at].query_start),
                ) <= self_scaf_min
                    || u64::from(
                        os.transitions[at]
                            .target_end
                            .wrapping_sub(os.transitions[at].target_start),
                    ) <= self_scaf_min;
                let mirror_short = u64::from(
                    os.transitions[mirror]
                        .query_end
                        .wrapping_sub(os.transitions[mirror].query_start),
                ) <= self_scaf_min
                    || u64::from(
                        os.transitions[mirror]
                            .target_end
                            .wrapping_sub(os.transitions[mirror].target_start),
                    ) <= self_scaf_min;
                if a_short || mirror_short {
                    os.transitions[at].deleted = true;
                    os.transitions[mirror].deleted = true;
                }
            }
            if os.transitions[at].deleted != os.transitions[mirror].deleted {
                os.transitions[at].deleted = false;
                os.transitions[mirror].deleted = false;
                os.transitions[at].occurrence = u32::MAX;
                os.transitions[mirror].occurrence = u32::MAX;
            }
        }
    }

    os.transitions.retain(|transition| !transition.deleted);
    kt_u_trans_t_idx(&mut os.transitions, &mut os.index, query_count);

    for q in 0..os.index.len() {
        let (start, count) = u_trans_span(&os.index, q);
        let a = &mut os.transitions[start..start + count];
        let mut reliable = 0usize;
        for k in 0..count {
            if a[k].occurrence != u32::MAX {
                if reliable != k {
                    a.swap(reliable, k);
                }
                reliable += 1;
            }
        }
        radix_sort_by_key_32(&mut a[..reliable], &|transition: &UTrans| {
            transition.query_start
        });
        radix_sort_by_key_32(&mut a[reliable..], &|transition: &UTrans| {
            transition.query_start
        });
        let unstable = count - reliable;
        if unstable == 0 {
            continue;
        }
        for item in &mut a[reliable..] {
            item.deleted = true;
        }
        let mut write = 0usize;
        for k in 0..unstable {
            if write > 0 && a[reliable + write - 1].query_end >= a[reliable + k].query_start {
                if a[reliable + k].query_end > a[reliable + write - 1].query_end {
                    a[reliable + write - 1].query_end = a[reliable + k].query_end;
                }
            } else {
                a[reliable + write] = a[reliable + k];
                write += 1;
            }
        }
        for k in write..unstable {
            let item = &mut a[reliable + k];
            item.query_identifier = u32::MAX;
            item.query_start = u32::MAX;
            item.query_end = u32::MAX;
            item.target_identifier = u32::MAX;
            item.target_start = u32::MAX;
            item.target_end = u32::MAX;
        }
    }
}

/// `scaf_hit2arc` (`Overlaps.cpp` 21876): turns one pair of transitions into
/// a scaffold-graph arc.
fn push_scaf_arc(graph: &mut Asg, a0: &UTrans, a1: &UTrans, reverse: bool, el: bool, glen: i64) {
    let mut v = (a0.target_identifier << 1) | u32::from(a0.reverse);
    v ^= u32::from(reverse);
    let mut w = (a1.target_identifier << 1) | u32::from(a1.reverse);
    w ^= u32::from(reverse);
    graph.push_arc(AsgArc {
        ul: u64::from(v) << 32,
        v: w,
        overlap: glen.unsigned_abs() as u32,
        deleted: false,
        occurrence: 0,
        strong: glen >= 0,
        without_large_indel: false,
        exact_region: u8::from(el),
    });
}

/// `gen_double_scaffold_gfa` (`Overlaps.cpp` 22689): turns the contig
/// transitions into two scaffold lists, one per haplotype.
#[allow(clippy::too_many_arguments)]
pub fn gen_double_scaffold_gfa(
    os: &mut TransTable,
    hu1: &UnitigGraph,
    hu2: &UnitigGraph,
    cp1: &ScafRes,
    cp2: &ScafRes,
    reference: &UnitigGraph,
    min_ns: i64,
    options: &AssemblyOptions,
) -> Result<(Vec<ScafSec>, Vec<ScafSec>), String> {
    let query_count = os.index.len() as u32;
    order_contig_trans(os, options.self_scaffolding_minimum_length, query_count);

    if let Some(directory) = crate::assembly::probe::dump_directory() {
        let mut text = String::new();
        text.push_str(&format!("N\t{}\t{}\n", os.transitions.len(), os.index.len()));
        for (i, entry) in os.index.iter().enumerate() {
            text.push_str(&format!("I\t{}\t{}\t{}\n", i, entry >> 32, *entry as u32));
        }
        for (i, t) in os.transitions.iter().enumerate() {
            text.push_str(&format!(
                "T\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.0}\n",
                i,
                t.query_identifier,
                t.query_start,
                t.query_end,
                t.target_identifier,
                t.target_start,
                t.target_end,
                u8::from(t.reverse),
                t.flag,
                u8::from(t.deleted),
                t.occurrence,
                t.network_weight
            ));
        }
        let _ = std::fs::write(directory.join("dual_os_ordered.txt"), text);
    }

    let reliable_min = options.self_scaffolding_reliable_minimum_length;
    let gap_max = options.self_scaffolding_maximum_gap_length;
    let gap_min = -gap_max;
    let hu1_count = hu1.unitigs.len();

    let mut graph = Asg::new();
    graph.sequences = Vec::with_capacity(os.index.len());
    for k in 0..os.index.len() {
        let length = if k < hu1_count {
            hu1.unitigs[k].length
        } else {
            hu2.unitigs[k - hu1_count].length
        };
        let mut sequence = AsgSeq::default();
        sequence.length = length;
        sequence.coverage_class = if k < hu1_count { 1 } else { 2 };
        graph.sequences.push(sequence);
    }

    let reference_count = reference.graph.sequences.len();
    let mut f = vec![0u8; (reference_count * 2).max(graph.sequences.len())];
    let mut fs = vec![0u8; reference_count * 2];
    let mut fsk = vec![0u8; reference_count];

    for k in 0..os.index.len() {
        let (start, count) = u_trans_span(&os.index, k);
        if count <= 1 {
            continue;
        }
        let transitions = &os.transitions[start..start + count];
        let mut pn = 0usize;
        while pn < count && !transitions[pn].deleted {
            pn += 1;
        }
        let mut dn = 0usize;
        while pn + dn < count && transitions[pn + dn].query_start != u32::MAX {
            dn += 1;
        }
        if pn <= 1 {
            continue;
        }

        let reliable = |t: &UTrans| {
            t.flag != 0
                || (u64::from(t.query_end - t.query_start) >= reliable_min
                    && u64::from(t.target_end - t.target_start) >= reliable_min)
        };
        let mut span_start: i64 = -1;
        let mut span_end: i64 = -1;
        let mut z = 0usize;
        while z < pn {
            if reliable(&transitions[z]) {
                span_start = z as i64;
                span_end = z as i64;
                break;
            }
            z += 1;
        }
        while z < pn {
            if reliable(&transitions[z]) {
                span_end = z as i64;
            }
            z += 1;
        }

        let query_length = graph.sequences[k].length;
        let scaf_path = if k < hu1_count {
            &cp1.a[k]
        } else {
            &cp2.a[k - hu1_count]
        };
        set_fsk(&mut fsk, 1, scaf_path);
        let noisy = &transitions[pn..pn + dn];
        for z in 1..pn {
            let previous = &transitions[z - 1];
            let current = &transitions[z];
            let previous_target = if (previous.target_identifier as usize) < hu1_count {
                &hu1.unitigs[previous.target_identifier as usize]
            } else {
                &hu2.unitigs[previous.target_identifier as usize - hu1_count]
            };
            if previous_target.circular {
                continue;
            }
            let current_target = if (current.target_identifier as usize) < hu1_count {
                &hu1.unitigs[current.target_identifier as usize]
            } else {
                &hu2.unitigs[current.target_identifier as usize - hu1_count]
            };
            if current_target.circular {
                continue;
            }

            let mut el = false;
            let mut gl: i64 = 0;
            let previous_cross = ((previous.target_identifier as usize) < hu1_count)
                != ((previous.query_identifier as usize) < hu1_count);
            let current_cross = ((current.target_identifier as usize) < hu1_count)
                != ((current.query_identifier as usize) < hu1_count);
            if previous_cross && current_cross {
                if (z as i64 - 1) >= span_start
                    && (z as i64 - 1) <= span_end
                    && (z as i64) >= span_start
                    && (z as i64) <= span_end
                {
                    el = true;
                }
                if el {
                    let zs = i64::from(previous.query_start.max(current.query_start));
                    let ze = i64::from(previous.query_end.min(current.query_end));
                    if is_noisy_ov(
                        zs.min(ze) as u64,
                        zs.max(ze) as u64,
                        u64::from(query_length),
                        2_500_000,
                        0.5,
                        noisy,
                    ) {
                        el = false;
                    }
                    if el {
                        gl = zs - ze;
                        if previous.flag == 0 || current.flag == 0 {
                            el = gl >= gap_min && gl <= gap_max;
                        } else {
                            el = gl <= gap_max;
                        }
                    }
                }
                let mut is_dip = 0u64;
                let glen = is_small_tangle(
                    previous,
                    current,
                    reference,
                    hu1,
                    hu2,
                    cp1,
                    cp2,
                    &mut f,
                    &mut fs,
                    &fsk,
                    5_000_000,
                    &mut is_dip,
                );
                if glen != u64::MAX && glen as i64 <= gap_max {
                    el = true;
                    if is_dip == 0 {
                        gl = glen as i64;
                    }
                }
            }
            if gl < min_ns {
                gl = min_ns;
            }
            if el && gl > gap_max {
                el = false;
            }
            push_scaf_arc(&mut graph, previous, current, false, el, gl);
            push_scaf_arc(&mut graph, current, previous, true, el, gl);
        }
        set_fsk(&mut fsk, 0, scaf_path);
    }

    graph.cleanup();
    graph.symmetrize();

    if let Some(directory) = crate::assembly::probe::dump_directory() {
        let mut text = String::new();
        for (i, arc) in graph.arcs.iter().enumerate() {
            text.push_str(&format!(
                "A\t{}\t{}\t{}\t{}\t{}\t{}\n",
                i,
                arc.ul >> 32,
                arc.v,
                arc.overlap,
                arc.exact_region,
                u8::from(arc.strong)
            ));
        }
        let _ = std::fs::write(directory.join("dual_raw_arcs.txt"), text);
    }

    // Drop every arc leaving a branching vertex.
    let vertex_count = (graph.sequences.len() * 2) as u32;
    for v in 0..vertex_count {
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        let (nv, _) = get_arcs(&graph, v);
        if nv <= 1 {
            continue;
        }
        let (start, count) = graph.arc_range(v);
        for k in start..start + count {
            if graph.arcs[k].deleted {
                continue;
            }
            let target = graph.arcs[k].v;
            let source = (graph.arcs[k].ul >> 32) as u32;
            graph.arcs[k].deleted = true;
            graph.delete_arc(target ^ 1, source ^ 1, true);
        }
    }
    for arc in &mut graph.arcs {
        if arc.exact_region == 0 {
            arc.deleted = true;
        }
    }
    graph.cleanup();

    let mut seen = vec![0u8; graph.sequences.len()];
    for value in &mut f[..graph.sequences.len()] {
        *value = 0;
    }
    let mut first = Vec::new();
    let mut second = Vec::new();
    let mut words: Vec<u64> = Vec::new();
    for v in 0..vertex_count {
        let (nv, _) = get_arcs(&graph, v);
        assert!(nv <= 1);
        if seen[(v >> 1) as usize] != 0 {
            continue;
        }
        if get_arcs(&graph, v ^ 1).0 != 0 {
            continue;
        }
        words.clear();
        let mut blen = 0u64;
        let mut bf = 0u64;
        let code = gen_scaf_path_lst(&graph, v, &mut words, &mut blen, &mut bf);
        assert_eq!(code, END_TIPS);
        let mut len = [0u64; 2];
        for &word in &words {
            let id = (word as u32 >> 1) as usize;
            seen[id] = 1;
            let class = graph.sequences[id].coverage_class as usize;
            len[class - 1] += u64::from(graph.sequences[id].length);
        }
        let target = if len[0] >= len[1] {
            &mut first
        } else {
            &mut second
        };
        target.push(ScafSec {
            is_c: false,
            a: words.clone(),
        });
    }
    for v in 0..vertex_count {
        let (nv, _) = get_arcs(&graph, v);
        if seen[(v >> 1) as usize] != 0 {
            continue;
        }
        assert_eq!(nv, 1);
        assert_eq!(get_arcs(&graph, v ^ 1).0, 1);
        words.clear();
        let mut blen = 0u64;
        let mut bf = 0u64;
        let code = gen_scaf_path_lst(&graph, v, &mut words, &mut blen, &mut bf);
        assert_eq!(code, LOOP);
        let mut len = [0u64; 2];
        for &word in &words {
            let id = (word as u32 >> 1) as usize;
            seen[id] = 1;
            let class = graph.sequences[id].coverage_class as usize;
            len[class - 1] += u64::from(graph.sequences[id].length);
        }
        let target = if len[0] >= len[1] {
            &mut first
        } else {
            &mut second
        };
        target.push(ScafSec {
            is_c: true,
            a: words.clone(),
        });
    }
    Ok((first, second))
}

/// `gen_scaffold_bases` (`Overlaps.cpp` 23162) with the `p == NULL` branch:
/// materialise the scaffold unitig sequences.
#[allow(clippy::too_many_arguments)]
pub fn gen_scaffold_bases(
    unitig_graph: &mut UnitigGraph,
    read_graph: &mut Asg,
    coverage_cut: &[CoverageCut],
    sources: &[PafList],
    arcs: &[crate::assembly::asg::AsgArc],
    max_hang: i32,
    min_ovlp: i32,
    reads: &ReadStore,
    options: &AssemblyOptions,
) {
    ma_ug_seq(
        unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        arcs,
        None,
        max_hang,
        min_ovlp,
        true,
        reads,
        options.maximum_hang_rate,
    );
}

/// `ma_scg_print` (`Overlaps.cpp` 11025): the segment/alignment lines of one
/// haplotype's scaffold list.
#[allow(clippy::too_many_arguments)]
pub fn ma_scg_print(
    scaffolds: &[ScafSec],
    reads: &ReadStore,
    read_graph: &Asg,
    coverage_cut: &[CoverageCut],
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    unitig_graph: &UnitigGraph,
    prefix: &str,
    path: &Path,
) -> std::io::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut out = std::io::BufWriter::new(file);
    let mut read_flags = vec![0u8; read_graph.sequences.len()];
    for (i, scaffold) in scaffolds.iter().enumerate() {
        let mut total_length = 0u64;
        let mut read_bases = 0u64;
        let mut covered_bases = 0u64;
        for &word in &scaffold.a {
            let unitig = &unitig_graph.unitigs[(word as u32 >> 1) as usize];
            total_length += u64::from(unitig.length);
            total_length += word >> 32;
            let (rb, cb) = get_ug_coverage_bases(
                unitig,
                read_graph,
                coverage_cut,
                sources,
                read_to_unitig,
                &mut read_flags,
            );
            read_bases += rb;
            covered_bases += cb;
        }
        let coverage = if read_bases == 0 {
            0
        } else {
            covered_bases / read_bases
        };
        let name = format!(
            "{}{:06}{}",
            prefix,
            i + 1,
            if scaffold.is_c { 'c' } else { 'l' }
        );
        write!(
            out,
            "S\t{}\t{}\tLN:i:{}\trd:i:{}\n",
            name,
            segment_sequence(scaffold, unitig_graph),
            total_length,
            coverage
        )?;

        let mut offset = 0u64;
        for &word in &scaffold.a {
            let unitig = &unitig_graph.unitigs[(word as u32 >> 1) as usize];
            let mut position = offset;
            if word as u32 & 1 == 0 {
                for &element in &unitig.reads {
                    let read_id = (element >> 33) as u32;
                    write_alignment_line(
                        &mut out,
                        &name,
                        position,
                        (element >> 32) & 1 != 0,
                        read_id,
                        reads,
                        coverage_cut,
                    )?;
                    position += element as u32 as u64;
                }
            } else {
                for element in unitig.reads.iter().rev() {
                    let read_id = (*element >> 33) as u32;
                    write_alignment_line(
                        &mut out,
                        &name,
                        position,
                        (element >> 32) & 1 == 0,
                        read_id,
                        reads,
                        coverage_cut,
                    )?;
                    position += *element as u32 as u64;
                }
            }
            offset = position;
            let gap = word >> 32;
            if gap > 0 {
                write!(
                    out,
                    "A\t{}\t{}\t+\tNs\t0\t{}\tid:i:-1\tHG:A:a\n",
                    name, position, gap
                )?;
            }
            offset += gap;
        }
    }
    out.flush()
}

/// The `S` line sequence of one scaffold: every unitig in its path
/// orientation with `N` gap runs between them.
fn segment_sequence(scaffold: &ScafSec, unitig_graph: &UnitigGraph) -> String {
    let mut sequence = Vec::new();
    for &word in &scaffold.a {
        let unitig = &unitig_graph.unitigs[(word as u32 >> 1) as usize];
        if word as u32 & 1 == 0 {
            sequence.extend_from_slice(&unitig.sequence);
        } else {
            for base in unitig.sequence.iter().rev() {
                sequence.push(if *base >= 128 {
                    b'N'
                } else {
                    complement(*base)
                });
            }
        }
        for _ in 0..(word >> 32) {
            sequence.push(b'N');
        }
    }
    String::from_utf8_lossy(&sequence).into_owned()
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        other => other,
    }
}

#[allow(clippy::too_many_arguments)]
fn write_alignment_line(
    out: &mut impl Write,
    name: &str,
    position: u64,
    reverse: bool,
    read_id: u32,
    reads: &ReadStore,
    coverage_cut: &[CoverageCut],
) -> std::io::Result<()> {
    if (read_id as usize) < reads.total_reads as usize {
        let flag = reads.trio_labels[read_id as usize];
        let flag = if flag != FATHER && flag != MOTHER {
            AMBIGU
        } else {
            flag
        };
        let hg = match flag {
            1 => 'p',
            2 => 'm',
            _ => 'a',
        };
        write!(
            out,
            "A\t{}\t{}\t{}\t{}\t{}\t{}\tid:i:{}\tHG:A:{}\n",
            name,
            position,
            if reverse { '-' } else { '+' },
            reads.read_name(read_id as usize),
            coverage_cut[read_id as usize].start,
            coverage_cut[read_id as usize].end,
            read_id,
            hg
        )
    } else {
        write!(
            out,
            "A\t{}\t{}\t{}\tFAKE\t{}\t{}\tid:i:{}\tHG:A:*\n",
            name,
            position,
            if reverse { '-' } else { '+' },
            coverage_cut[read_id as usize].start,
            coverage_cut[read_id as usize].end,
            read_id
        )
    }
}

/// `gen_self_scaf` (`Overlaps.cpp` 23209): the `--dual-scaf` driver.
#[allow(clippy::too_many_arguments)]
pub fn gen_self_scaf(
    read_graph: &mut Asg,
    coverage_cut: &mut [CoverageCut],
    output_file_name: &str,
    hu0: &mut UnitigGraph,
    hu1: &mut UnitigGraph,
    sources: &[PafList],
    reverse_sources: &[PafList],
    tips_len: i64,
    tip_drop_ratio: f32,
    stops_threshold: u32,
    read_to_unitig: &mut ReadToUnitig,
    chimeric_rate: f32,
    drop_ratio: f32,
    max_hang: i32,
    min_ovlp: i32,
    bubble_labels: &mut BubbleLabel,
    options: &mut AssemblyOptions,
    trio_flags: &mut [u8],
    reads: &ReadStore,
    arcs: &[crate::assembly::asg::AsgArc],
) -> Result<(), String> {
    let mut reference = ma_ug_gen(read_graph);
    let cp0 = crate::assembly::ul::gen_contig_path(hu0, &reference, &reads.read_lengths);
    let cp1 = crate::assembly::ul::gen_contig_path(hu1, &reference, &reads.read_lengths);

    // `gen_bubble_chain(sg, ref, opt, &bf, 0)`.
    let bubble = crate::assembly::trio_dual::gen_bubble_chain(
        read_graph,
        &mut reference,
        coverage_cut,
        sources,
        reverse_sources,
        read_to_unitig,
        tips_len,
        tip_drop_ratio,
        stops_threshold,
        chimeric_rate,
        drop_ratio,
        max_hang,
        min_ovlp,
        bubble_labels,
        options,
        trio_flags,
    )?;
    let shallow = gen_shallow_ref_trans(
        &bubble,
        &cp0,
        &cp1,
        read_graph,
        &reference,
        coverage_cut,
        sources,
        reverse_sources,
        tips_len,
        tip_drop_ratio,
        stops_threshold,
        read_to_unitig,
        chimeric_rate,
        drop_ratio,
        max_hang,
        min_ovlp,
        bubble_labels,
        options,
        trio_flags,
    )?;
    let deep = gen_deep_ref_trans(
        Some(&shallow),
        &cp0,
        &cp1,
        read_graph,
        &reference,
        coverage_cut,
        sources,
        reverse_sources,
        tips_len,
        tip_drop_ratio,
        stops_threshold,
        read_to_unitig,
        chimeric_rate,
        drop_ratio,
        max_hang,
        min_ovlp,
        bubble_labels,
        options,
        trio_flags,
        reads,
    )?;

    // `gen_contig_trans_func` with `is_sep == 0`: merge the two haplotypes.
    let mut merged = merge_utgs(hu0, hu1);
    let merged_scaf = merge_scaf_res(&cp0, &cp1);
    let mut transitions = crate::assembly::ctg_trans::gen_contig_self(
        &merged,
        &merged_scaf,
        &reference,
        &deep,
        hu0.unitigs.len() as u32,
        &bubble,
        merged.unitigs.len(),
    );

    if let Some(directory) = crate::assembly::probe::dump_directory() {
        let mut text = String::new();
        for (i, t) in transitions.transitions.iter().enumerate() {
            text.push_str(&format!(
                "T {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                i,
                t.query_identifier,
                t.query_start,
                t.query_end,
                t.target_identifier,
                t.target_start,
                t.target_end,
                u8::from(t.reverse),
                t.flag,
                t.network_weight,
                t.occurrence
            ));
        }
        let _ = std::fs::write(directory.join("dual_trans.txt"), text);
        let mut index = String::new();
        for (i, entry) in transitions.index.iter().enumerate() {
            index.push_str(&format!("I {}\t{}\t{}\n", i, entry >> 32, *entry as u32));
        }
        let _ = std::fs::write(directory.join("dual_trans_idx.txt"), index);
    }
    let (sc0, sc1) =
        gen_double_scaffold_gfa(&mut transitions, hu0, hu1, &cp0, &cp1, &reference, 100, options)?;
    gen_scaffold_bases(
        &mut merged,
        read_graph,
        coverage_cut,
        sources,
        arcs,
        max_hang,
        min_ovlp,
        reads,
        options,
    );

    let f0 = format!("{}.hap1.p_ctg.gfa", output_file_name);
    eprintln!("Writing {} to disk... ", f0);
    ma_scg_print(
        &sc0,
        reads,
        read_graph,
        coverage_cut,
        sources,
        read_to_unitig,
        &merged,
        "h1tg",
        Path::new(&f0),
    )
    .map_err(|error| format!("failed to write {}: {}", f0, error))?;
    let f1 = format!("{}.hap2.p_ctg.gfa", output_file_name);
    eprintln!("Writing {} to disk... ", f1);
    ma_scg_print(
        &sc1,
        reads,
        read_graph,
        coverage_cut,
        sources,
        read_to_unitig,
        &merged,
        "h2tg",
        Path::new(&f1),
    )
    .map_err(|error| format!("failed to write {}: {}", f1, error))?;

    let _ = kt_u_trans_t_idx;
    Ok(())
}

//! `trans_base_infer` (`inter.cpp` 22574): base-level transitions inferred
//! by mapping every unitig against the unitig graph itself with a global
//! minimizer index.
//!
//! The index is built in two sketch passes (`gen_trans_base_count_comp`,
//! `inter.cpp` 22181): the first counts the `mz2` minimizers of every
//! non-homozygous unitig, the second fills the flat array. The array is then
//! sorted by hash, every entry's `x` is remapped to its sorted position, and
//! the occupancy quantile tightens `mini_cut`. The sorted array doubles as
//! the global lookup index of `minimizers_qgen_input`.

use crate::assembly::bubble_chain::BubbleType;
use crate::assembly::unitig::UnitigGraph;
use crate::command_line::options::AssemblyOptions;
use crate::container::ksort::radix_sort_by_key;
use crate::counting::minimizer::sketch;

/// One global minimizer, `ha_mzl_t`: `x` is the hash key before the remap
/// and the sorted position afterwards; `rid`/`pos`/`rev`/`span` locate the
/// k-mer on its unitig.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MzlEntry {
    pub x: u64,
    pub rid: u32,
    pub pos: u32,
    pub reverse: bool,
    pub span: u8,
}

/// The global minimizer index produced by the two count passes.
pub struct TransBaseIndex {
    /// `idx_a`: per-unitig minimizers, `x` remapped to the sorted position.
    pub entries: Vec<MzlEntry>,
    /// `idx_n`: prefix offsets, `entries[offsets[i]..offsets[i+1]]` belongs
    /// to unitig `i`.
    pub offsets: Vec<u64>,
    /// `srt_a`: the hash-sorted global index (original hash in `x`).
    pub sorted: Vec<MzlEntry>,
    /// `mini_cut` after the occupancy-quantile adjustment.
    pub mini_cut: u64,
}

/// `mz2_ha_sketch` through the shared `HA_SC_INIT` selection: the `mz2`
/// instantiation differs from `mz1` only in the field widths, and this call
/// path passes no filter table and no position index, so the `mz1` port
/// reproduces it exactly.
fn sketch_unitig(
    unitig_graph: &UnitigGraph,
    id: usize,
    options: &AssemblyOptions,
) -> Vec<MzlEntry> {
    let result = sketch(
        &unitig_graph.unitigs[id].sequence,
        options.transition_window_size,
        options.transition_k_mer_length,
        id as u32,
        false,
        None,
        options.minimizer_sample_distance,
        options.minimizer_rewindow_size,
    );
    result
        .minimizers
        .iter()
        .map(|m| MzlEntry {
            x: m.hash,
            rid: id as u32,
            pos: m.position,
            reverse: m.reverse,
            span: m.span as u8,
        })
        .collect()
}

/// Whether `worker_for_trans_ovlp_adv` skips the unitig: homozygous bubble
/// nodes and deleted sequences take no part in the base transitions.
fn skipped(unitig_graph: &UnitigGraph, bubble: &BubbleType, id: usize) -> bool {
    bubble.is_hom(id as u32) || unitig_graph.graph.sequences[id].deleted
}

/// The count/fill/sort/remap half of `gen_trans_base_count_comp`
/// (`inter.cpp` 22181): builds the global index and adjusts `mini_cut`.
pub fn build_trans_base_index(
    unitig_graph: &UnitigGraph,
    bubble: &BubbleType,
    options: &AssemblyOptions,
) -> TransBaseIndex {
    let n = unitig_graph.unitigs.len();
    let mut mini_cut = 512u64;

    // First pass: count the minimizers of every active unitig.
    let mut offsets = vec![0u64; n + 1];
    let mut total = 0u64;
    for i in 0..n {
        if !skipped(unitig_graph, bubble, i) {
            total += sketch_unitig(unitig_graph, i, options).len() as u64;
        }
        offsets[i + 1] = total;
    }

    // Second pass: fill the flat array in unitig order.
    let mut entries: Vec<MzlEntry> = Vec::with_capacity(total as usize);
    for i in 0..n {
        if skipped(unitig_graph, bubble, i) {
            continue;
        }
        entries.extend(sketch_unitig(unitig_graph, i, options));
    }
    assert_eq!(entries.len() as u64, total);

    // srt_a keeps the original index packed into rid/pos while sorting.
    let mut sorted: Vec<MzlEntry> = entries
        .iter()
        .enumerate()
        .map(|(i, e)| MzlEntry {
            x: e.x,
            rid: (i as u64 >> 32) as u32,
            pos: i as u32,
            reverse: e.reverse,
            span: e.span,
        })
        .collect();
    radix_sort_by_key(&mut sorted, &|e: &MzlEntry| e.x);

    // Remap every entry to its sorted position and collect group sizes.
    let mut cut: Vec<u64> = Vec::new();
    let mut l = 0usize;
    for k in 1..=sorted.len() {
        if k == sorted.len() || sorted[l].x != sorted[k].x {
            for i in l..k {
                let m = ((sorted[i].rid as u64) << 32) | sorted[i].pos as u64;
                assert_eq!(sorted[i].x, entries[m as usize].x);
                sorted[i] = entries[m as usize];
                entries[m as usize].x = i as u64;
            }
            cut.push((k - l) as u64);
            l = k;
        }
    }

    if !cut.is_empty() {
        radix_sort_by_key(&mut cut, &|&c: &u64| c);
        let m = (cut.len() as f64 * 0.0002) as usize;
        let mut cc = cut[cut.len() - 1] + 1;
        if m > 0 && m <= cut.len() {
            cc = cut[cut.len() - m] + 1;
        }
        if cc < mini_cut {
            mini_cut = cc;
        }
    }

    if let Some(directory) = crate::assembly::probe::dump_directory() {
        let mut out = format!("N {} {}\n", sorted.len(), mini_cut);
        for (i, e) in sorted.iter().enumerate() {
            out.push_str(&format!(
                "S {} {} {} {} {} {}\n",
                i, e.x, e.rid, e.pos, u8::from(e.reverse), e.span
            ));
        }
        for (i, &off) in offsets.iter().enumerate() {
            out.push_str(&format!("O {} {}\n", i, off));
        }
        for (i, e) in entries.iter().enumerate() {
            out.push_str(&format!(
                "I {} {} {} {} {} {}\n",
                i, e.x, e.rid, e.pos, u8::from(e.reverse), e.span
            ));
        }
        let _ = std::fs::write(directory.join("trans_base_idx.txt"), out);
    }

    TransBaseIndex {
        entries,
        offsets,
        sorted,
        mini_cut,
    }
}
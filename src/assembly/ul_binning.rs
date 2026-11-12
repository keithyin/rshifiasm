//! The graph-binning driver (`inter.cpp` 22599-22875): builds the per-unitig
//! minimizer index (`ug_bin_t::idx_a/idx_n/srt_a`), maps every unitig against
//! it (`ug_map_lchain_simple` = sketch + `minimizers_qgen_input` +
//! `lchain_qgen_mcopy`), refines the chains with `ug_lalign` and fills the
//! `idx_emask_t` table through `push_graph_bin`.

use crate::assembly::trans_base::{MzlEntry, TransBaseIndex};
use crate::assembly::ul_mask::{
    push_graph_bin, GBIN_BE, GBIN_E, GBIN_K, GBIN_L, GBIN_ME, GBIN_W,
};
use crate::assembly::ul_phase::Emask;
use crate::assembly::ul_phase::{EmaskTable, MaskUlOv};
use crate::assembly::ul_rescall::UlRescallContext;
use crate::command_line::options::AssemblyOptions;
use crate::container::ksort::{radix_sort_by_key, radix_sort_u64};
use crate::correction::hc_align::BitExtz;
use crate::correction::target::TargetSource;
use crate::correction::ug_lalign::ug_lalign_public;
use crate::correction::ul_align::gen_aux_ovlp;
use crate::overlap::lchain_parameters::set_lchain_parameters;
use crate::overlap::qgen_input::{InputChainParameters, SeedInput, lchain_qgen_mcopy_input, minimizers_qgen_input};
use crate::overlap::region::CandidatesList;

/// The `mg_idxopt_t` of the binning pass: `k = GBIN_K`, `w = GBIN_W`.
const BINNING_KMER: u32 = GBIN_K as u32;
const BINNING_WINDOW: u32 = GBIN_W as u32;

/// The unitig sequence, materialised through `retrieve_u_seq` when the
/// graph does not carry it.
fn materialize_unitig_sequence(
    unitig_graph: &crate::assembly::unitig::UnitigGraph,
    reads: &crate::reads::ReadStore,
    id: usize,
    cached: Vec<u8>,
) -> Vec<u8> {
    if !cached.is_empty() {
        return cached;
    }
    let length = unitig_graph.unitigs[id].length as i64;
    crate::correction::target::unitig_sub_region(
        &unitig_graph.unitigs,
        reads,
        id,
        0,
        length,
        0,
    )
}

/// `ug_map_lchain_simple`'s sketch half (`mz2_ha_sketch` with
/// `is_HPC = 0`), via the shared `mz1` port the trans-base path uses.
fn sketch_unitig_for_binning(
    sequence: &[u8],
    id: usize,
    options: &AssemblyOptions,
) -> Vec<MzlEntry> {
    let result = crate::counting::minimizer::sketch(
        sequence,
        BINNING_WINDOW as i32,
        BINNING_KMER as i32,
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

/// The count/fill/sort/remap half of `cal_graph_ovlp_binning` (22789): the
/// global minimizer index over the unitig sequences (k = 31, w = 31).
pub fn build_binning_index(
    unitig_graph: &crate::assembly::unitig::UnitigGraph,
    reads: &crate::reads::ReadStore,
    options: &AssemblyOptions,
) -> TransBaseIndex {
    let n = unitig_graph.unitigs.len();
    let mut mini_cut = 512u64;

    let sequences: Vec<Vec<u8>> = unitig_graph
        .unitigs
        .iter()
        .enumerate()
        .map(|(id, u)| materialize_unitig_sequence(unitig_graph, reads, id, u.sequence.clone()))
        .collect();

    // First pass: count the minimizers of every non-deleted unitig.
    let mut offsets = vec![0u64; n + 1];
    let mut total = 0u64;
    for i in 0..n {
        if !unitig_graph.graph.sequences[i].deleted {
            total += sketch_unitig_for_binning(&sequences[i], i, options).len() as u64;
        }
        offsets[i + 1] = total;
    }

    // Second pass: fill the flat array in unitig order.
    let mut entries: Vec<MzlEntry> = Vec::with_capacity(total as usize);
    for i in 0..n {
        if unitig_graph.graph.sequences[i].deleted {
            continue;
        }
        entries.extend(sketch_unitig_for_binning(&sequences[i], i, options));
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
        radix_sort_u64(&mut cut);
        let m = (cut.len() as f64 * 0.0002) as usize;
        let mut cc = cut[cut.len() - 1] + 1;
        if m > 0 && m <= cut.len() {
            cc = cut[cut.len() - m] + 1;
        }
        if cc < mini_cut {
            mini_cut = cc;
        }
    }

    TransBaseIndex {
        entries,
        offsets,
        sorted,
        mini_cut,
    }
}

/// The per-unitig overlap pass of `worker_for_graph_bin` (12442): maps the
/// unitig, refines twice and pushes the emask entry.
#[allow(clippy::too_many_arguments)]
fn bin_one_unitig(
    context: &UlRescallContext<'_>,
    index: &TransBaseIndex,
    options: &AssemblyOptions,
    sequences: &[Vec<u8>],
    target_lengths: &[u32],
    rid: usize,
    extz: &mut BitExtz,
    b0: &mut Vec<u64>,
    mk: &mut MaskUlOv,
    entry: &mut Vec<Emask>,
) {
    let unitig = &context.unitig_graph.unitigs[rid];
    let sequence = &sequences[rid];
    let length = sequence.len() as u64;

    let query_minimizers = &index.entries[index.offsets[rid] as usize..index.offsets[rid + 1] as usize];
    let high_occurrence = (options.polyploidy + 1) as u32;
    let seed = SeedInput {
        query_identifier: rid as u32,
        query_minimizers,
        index: &index.sorted,
        minimizer_cutoff: index.mini_cut,
        chain_cutoff: 3,
        known_transitions: &[],
        target_lengths,
        high_occurrence: Some(high_occurrence),
        low_occurrence: None,
    };
    let anchors = minimizers_qgen_input(&seed);

    let dp = set_lchain_parameters(true, BINNING_KMER);
    let parameters = InputChainParameters {
        maximum_chain_count: options.maximum_chain_count as u64,
        maximum_skip: dp.maximum_skip,
        maximum_iterations: dp.maximum_iterations,
        maximum_distance: dp.maximum_distance,
        gap_penalty: dp.gap_penalty,
        skip_penalty: dp.skip_penalty,
        band_width_rate: GBIN_E,
        band_width_rate_sec: GBIN_E,
        append_both_ends: true,
        generate_offset: true,
    };
    let (mut overlaps, candidate_list) = lchain_qgen_mcopy_input(
        anchors.candidates,
        rid as u32,
        length,
        target_lengths,
        &anchors.group_starts,
        &parameters,
    );
    let mut candidates = CandidatesList {
        list: candidate_list,
        length: 0,
        size: 0,
        chain_data: Default::default(),
    };
    candidates.length = candidates.list.len() as i64;
    candidates.size = candidates.length;

    // `erate = s->diff_ov; wl = MIN(THRESHOLD_MAX_SIZE / erate, WINDOW)`.
    let erate = GBIN_E;
    let wl = ((31.0 / erate) as i64).min(375);
    ug_lalign_public(
        &mut overlaps,
        &mut candidates,
        &TargetSource::Unitigs {
            unitigs: &context.unitig_graph.unitigs,
            reads: context.reads,
        },
        sequence,
        length as i64,
        extz,
        None,
        erate,
        wl,
        BINNING_KMER as i64,
        3,
    );
    let mut auxiliary = gen_aux_ovlp(&overlaps, overlaps.len());
    ug_lalign_public(
        &mut overlaps,
        &mut candidates,
        &TargetSource::Unitigs {
            unitigs: &context.unitig_graph.unitigs,
            reads: context.reads,
        },
        sequence,
        length as i64,
        extz,
        Some(&mut auxiliary),
        erate,
        wl,
        BINNING_KMER as i64,
        3,
    );

    push_graph_bin(
        &mut overlaps,
        &context.unitig_graph.unitigs,
        &context.unitig_graph.graph,
        GBIN_BE,
        GBIN_ME,
        GBIN_L,
        rid as u32,
        b0,
        mk,
        entry,
    );
}

/// `graph_ovlp_binning` (`inter.cpp` 22860): the whole `idx_emask_t`
/// construction; `sg` is the read graph (unused after the seq check).
pub fn graph_ovlp_binning(
    context: &UlRescallContext<'_>,
    options: &AssemblyOptions,
) -> EmaskTable {
    let mut mm: EmaskTable = vec![Vec::new(); context.unitig_graph.graph.sequences.len()];
    let index = build_binning_index(context.unitig_graph, context.reads, options);
    let sequences: Vec<Vec<u8>> = context
        .unitig_graph
        .unitigs
        .iter()
        .enumerate()
        .map(|(id, u)| {
            materialize_unitig_sequence(context.unitig_graph, context.reads, id, u.sequence.clone())
        })
        .collect();
    let target_lengths: Vec<u32> = context
        .unitig_graph
        .unitigs
        .iter()
        .map(|u| u.length)
        .collect();
    let mut extz = BitExtz::new();
    let mut b0: Vec<u64> = Vec::new();
    let mut mk = MaskUlOv::default();
    for rid in 0..context.unitig_graph.unitigs.len() {
        if context.unitig_graph.graph.sequences[rid].deleted {
            continue;
        }
        bin_one_unitig(
            context,
            &index,
            options,
            &sequences,
            &target_lengths,
            rid,
            &mut extz,
            &mut b0,
            &mut mk,
            &mut mm[rid],
        );
    }
    mm
}

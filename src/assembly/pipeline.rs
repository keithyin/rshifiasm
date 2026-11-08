//! Top-level assembly-graph driver, the counterpart of
//! `build_string_graph_without_clean` and `clean_graph` in `Overlaps.cpp`,
//! plus their `gfa_ut.cpp` callees `ul_clean_gfa` and `post_rescue`.
//!
//! The correction stage hands over its final overlap lists; everything below
//! is the reference's graph construction, cleaning and output sequence.

use crate::assembly::asg::Asg;
use crate::assembly::bubble::{BubbleBuffer, BubbleInfo, BubbleLabel};
use crate::assembly::bubble_identify::{
    asg_arc_identify_simple_bubbles_multi, flat_soma_v,
};
use crate::assembly::gfa_clean::{
    asg_arc_cut_bub_links, asg_arc_cut_chimeric, asg_arc_cut_complex_bub_links,
    asg_arc_cut_inexact, asg_arc_cut_length, asg_arc_cut_tips, asg_cut_large_indel,
    asg_cut_semi_circ, asg_iterative_semi_circ, TeloEnd, LIM_LEN,
};
use crate::assembly::graph_build::gen_init_sg;
use crate::assembly::hit::{
    clean_weak_ma_hit_t, normalize_ma_hit_t_single_side_advance, renew_graph_init,
    try_rescue_overlaps, CoverageCut, ReadToUnitig,
};
use crate::assembly::probe;
use crate::assembly::purge_dups::set_hom_global_coverage;
use crate::assembly::rescue::{
    rescue_contained_reads_aggressive, rescue_missing_overlaps_aggressive,
    rescue_missing_overlaps_backward, RescueParameters,
};
use crate::assembly::unitig::{output_contig_graph_primary_pre, output_unitig_graph, ug_ext_gfa};
use crate::command_line::options::{
    AssemblyOptions, FLAG_BAN_ASSEMBLY, FLAG_PARTITION, MINIMUM_CHAIN_COUNT,
};
use crate::correction::ec_driver::{CorrectionState, PafList};
use crate::overlap::paf::PafHit;
use crate::reads::read_store::ReadStore;

/// `HARD_OU_DROP` in `ul_clean_gfa`.
const HARD_OU_DROP: f32 = 0.75;
/// `HARD_OL_DROP` in `ul_clean_gfa`.
const HARD_OL_DROP: f32 = 0.6;
/// `HARD_OL_SEC_DROP` in `ul_clean_gfa`.
const HARD_OL_SEC_DROP: f32 = 0.85;
/// `HARD_ORTHOLOGY_DROP` in `ul_clean_gfa`.
const HARD_ORTHOLOGY_DROP: f32 = 0.4;
/// `ug_ext_len` in `Overlaps.h`.
const UG_EXT_LEN: u32 = 75_000;

/// `get_outfile_name`: the prefix shared by every assembly output.
pub fn get_outfile_name(options: &AssemblyOptions) -> String {
    if options.trio_binning_requested() && options.hic_partition_requested() {
        format!("{}.hic.bench", options.output_file_name)
    } else if options.trio_binning_requested() {
        format!("{}.dip", options.output_file_name)
    } else if options.hic_partition_requested() {
        format!("{}.hic", options.output_file_name)
    } else if options.flags & FLAG_PARTITION != 0 {
        format!("{}.bp", options.output_file_name)
    } else {
        options.output_file_name.clone()
    }
}

/// The per-read flags attached to a [`PafList`] and preserved across the
/// raw-list detour.
type PafFlags = (bool, u8);

fn split_paf_lists(lists: Vec<PafList>) -> (Vec<Vec<PafHit>>, Vec<PafFlags>) {
    let mut raw = Vec::with_capacity(lists.len());
    let mut flags = Vec::with_capacity(lists.len());
    for list in lists {
        raw.push(list.hits);
        flags.push((list.is_fully_corrected, list.is_abnormal));
    }
    (raw, flags)
}

fn join_paf_lists(raw: Vec<Vec<PafHit>>, flags: Vec<PafFlags>) -> Vec<PafList> {
    raw.into_iter()
        .zip(flags)
        .map(|(hits, (is_fully_corrected, is_abnormal))| PafList {
            hits,
            is_fully_corrected,
            is_abnormal,
        })
        .collect()
}

/// `build_string_graph_without_clean` followed by `clean_graph`, stopping at
/// the output dispatch.
///
/// The default HiFi path (`HA_F_PARTITION`, `purge_level_primary > 0`) ends in
/// `output_bp_graph`, which is not ported yet; the function therefore returns
/// an error after the read-graph and primary-unitig outputs have been written.
pub fn build_string_graph_without_clean(
    options: &mut AssemblyOptions,
    reads: &ReadStore,
    state: CorrectionState,
    write_caches: bool,
) -> Result<(), String> {
    // `ha_opt_update_cov_min(&asm_opt, asm_opt.hom_cov, MIN_N_CHAIN)`.
    options.update_coverage_with_minimum(options.homozygous_coverage, MINIMUM_CHAIN_COUNT);

    let min_dp = options.minimum_overlap_coverage;
    let minimum_overlap = options.minimum_overlap_length;
    let maximum_hang = options.maximum_hang_length;
    let minimum_overlap_length = i64::from(options.minimum_overlap_length);
    let gap_fuzz = options.gap_fuzz;
    let is_ont = options.ont_mode != 0;
    let max_short_tip = options.maximum_short_tip_read_count;
    let min_thres = i64::from(max_short_tip + 1);

    let mut read_to_unitig = ReadToUnitig::new(reads.total_reads as usize);
    // `R_INF.trio_flag`: one live array threaded through the cleaning, the
    // bubble rescue and the output dispatch (the reference mutates the
    // global in place).
    let mut trio_flags = reads.trio_labels.clone();

    // `if (asm_opt.write_index_to_disk && write) write_all_data_to_disk(...)`
    // (`Overlaps.cpp` 39707), before `try_rescue_overlaps`: the caches hold
    // the pre-rescue overlaps. A run that loaded them passes `write = 0`.
    if write_caches && options.write_index_to_disk != 0 {
        crate::assembly::read_cache::write_all_data_to_disk(
            options,
            reads,
            &state.forward,
            &state.reverse,
        )
        .map_err(|error| format!("failed to write the restart caches: {}", error))?;
    }

    let (mut forward, forward_flags) = split_paf_lists(state.forward);
    let (mut reverse, reverse_flags) = split_paf_lists(state.reverse);

    // `if (asm_opt.flag & HA_F_BAN_ASSEMBLY)` (`Overlaps.cpp` 39714): with
    // `-e` the reference stops after the caches are written and exits 0.
    if options.flags & FLAG_BAN_ASSEMBLY != 0 {
        return Ok(());
    }

    // Stage dumps for local debugging; see `assembly::probe`.
    let dump_directory = probe::prepare_directory().ok().flatten();

    // `try_rescue_overlaps(sources, reverse_sources, n_read, 4, is_ont)`.
    try_rescue_overlaps(
        &mut forward,
        &mut reverse,
        &reads.read_lengths,
        4,
        is_ont,
    );
    probe::try_dump_stage_raw(
        dump_directory.as_deref(),
        "d00_after_try_rescue",
        None,
        &forward,
        &forward_flags,
        &reverse,
        &reverse_flags,
        None,
        &read_to_unitig,
        &trio_flags,
    );

    // `clean_graph` head: `if(!cmk) renew_graph_init(...)`. On the HiFi path
    // `cmk` is NULL, so every deletion flag is cleared and the read-to-unitig
    // map is reset to the unmapped sentinel.
    renew_graph_init(&mut forward, &mut reverse, &mut read_to_unitig);

    // `clean_graph`: mirror both sides, then drop the weakly supported hits.
    normalize_ma_hit_t_single_side_advance(&mut forward, is_ont, None, 0);
    normalize_ma_hit_t_single_side_advance(&mut reverse, false, None, 0);
    probe::try_dump_stage_raw(
        dump_directory.as_deref(),
        "d00a_after_normalize",
        None,
        &forward,
        &forward_flags,
        &reverse,
        &reverse_flags,
        None,
        &read_to_unitig,
        &trio_flags,
    );
    // `if (ha_opt_triobin(&asm_opt)) drop_edges_by_trio(sources, n_read);
    // else memset(R_INF.trio_flag, AMBIGU, ...)` — the else half is a no-op
    // here because the store already holds `TRIO_UNKNOWN` (0) everywhere.
    if options.trio_binning_requested() {
        crate::assembly::hit::drop_edges_by_trio(&mut forward, &trio_flags);
    }
    probe::try_dump_stage_raw(
        dump_directory.as_deref(),
        "d00b_after_drop_edges",
        None,
        &forward,
        &forward_flags,
        &reverse,
        &reverse_flags,
        None,
        &read_to_unitig,
        &trio_flags,
    );
    clean_weak_ma_hit_t(&mut forward, &reverse, &reads.read_lengths, u32::MAX);
    probe::try_dump_stage_raw(
        dump_directory.as_deref(),
        "d01_after_weak_clean",
        None,
        &forward,
        &forward_flags,
        &reverse,
        &reverse_flags,
        None,
        &read_to_unitig,
        &trio_flags,
    );

    let initial = gen_init_sg(
        min_dp,
        &mut forward,
        &reads.read_lengths,
        &mut read_to_unitig,
        maximum_hang,
        minimum_overlap,
        gap_fuzz,
        minimum_overlap_length,
        options.thread_count,
        reads,
        options,
    );

    // The cast of the corrected overlaps is over; the later stages only read
    // the lists, so they go back to the flagged representation.
    let sources = join_paf_lists(forward, forward_flags);
    let reverse_sources = join_paf_lists(reverse, reverse_flags);

    let mut graph = initial.graph;
    let mut bubble_labels = initial.bubble_label;
    let mut coverage_cut = initial.coverage_cut;
    let output_file_name = get_outfile_name(options);

    probe::try_dump_stage(
        dump_directory.as_deref(),
        "d02_after_gen_init_sg",
        Some(&graph),
        &sources,
        &reverse_sources,
        Some(&coverage_cut),
        &read_to_unitig,
        &trio_flags,
    );

    // `if(asm_opt.telo_motif) te = gen_telo_end_t(...)` (`Overlaps.cpp` 39347).
    let telo_end = options.telomere_motif.as_ref().map(|motif| {
        crate::assembly::gfa_clean::gen_telo_ends(
            reads,
            motif,
            options.telomere_minimum_score,
            options.telomere_penalty,
            options.telomere_maximum_drop,
        )
    });

    ul_clean_gfa(
        options,
        reads,
        &mut trio_flags,
        &mut graph,
        &sources,
        &reverse_sources,
        &mut read_to_unitig,
        &mut coverage_cut,
        &mut bubble_labels,
        &output_file_name,
        min_thres,
        telo_end.as_ref(),
        dump_directory.as_deref(),
    )?;

    // `clean_graph`, non-`ar` tail.
    ug_ext_gfa(
        &mut graph,
        &sources,
        &read_to_unitig,
        (max_short_tip * 2) as u32,
        maximum_hang,
        minimum_overlap,
        options.maximum_hang_rate,
        UG_EXT_LEN,
    );
    // `hic_clean_adv(sg, &uopt)` follows here in the reference. It needs the
    // Hi-C bubble-chain machinery that is not ported yet; the stage does not
    // change the graph on the chr11-2M corpus (probe d19 == d20).
    probe::try_dump_stage(
        dump_directory.as_deref(),
        "d20_after_hic_clean",
        Some(&graph),
        &sources,
        &reverse_sources,
        Some(&coverage_cut),
        &read_to_unitig,
        &trio_flags,
    );

    output_contig_graph_primary_pre(
        &mut graph,
        &coverage_cut,
        &output_file_name,
        &sources,
        &read_to_unitig,
        options.small_bubble_population_size as u64,
        maximum_hang,
        minimum_overlap,
        reads,
        options.maximum_hang_rate,
        options.bed_inconsistency_rate as u32,
    )
    .map_err(|error| format!("failed to write the primary unitig GFA: {}", error))?;

    // Output dispatch of `clean_graph` (`Overlaps.cpp` 39605-39651), in the
    // reference's order: polyploid, trio+Hi-C benchmark, trio, Hi-C, default.
    if options.trio_polyploid_file_name.is_some() {
        return Err("output_poly_trio (-5) is not implemented yet".to_string());
    }
    if options.trio_binning_requested() && options.hic_partition_requested() {
        return Err("benchmark_hic_graph is not implemented yet".to_string());
    }
    if options.trio_binning_requested() {
        if options.trio_heterozygous_overlap_coverage > 0 {
            return Err("output_bp_trio_graph (--trio-dual) is not implemented yet".to_string());
        }
        // `if(asm_opt.flag & HA_F_PARTITION) asm_opt.flag -= HA_F_PARTITION;`
        options.flags &= !FLAG_PARTITION;
        // `output_trio_graph_joint(sg, coverage_cut, o_file, sources,
        // reverse_sources, (asm_opt.max_short_tip*2), 0.15, 3, ruIndex,
        // 0.05, 0.9, max_hang_length, mini_overlap_length, gap_fuzz,
        // &b_mask_t, NULL, NULL, &uopt)`.
        let joint_status = crate::assembly::bp_graph::output_trio_graph_joint(
            &mut graph,
            &mut coverage_cut,
            &output_file_name,
            &sources,
            &reverse_sources,
            i64::from(options.maximum_short_tip_read_count) * 2,
            0.15,
            3,
            &mut read_to_unitig,
            0.05,
            0.9,
            maximum_hang,
            minimum_overlap,
            options.gap_fuzz,
            &mut bubble_labels,
            options,
            &mut trio_flags,
            reads,
        );
        // The reference prints the BED threshold at the end of `clean_graph`,
        // after every dispatch branch.
        eprintln!(
            "Inconsistency threshold for low-quality regions in BED files: {}%",
            options.bed_inconsistency_rate
        );
        return joint_status;
    }
    if options.hic_partition_requested() {
        return Err("output_hic_graph is not implemented yet".to_string());
    }
    if options.flags & FLAG_PARTITION != 0 && options.purge_level_primary > 0 {
        // The reference's `R_INF.trio_flag` starts out without any trio
        // evidence and is filled by the purge/trio stages; keep a local copy
        // so the output drivers can mutate it.
        return crate::assembly::bp_graph::output_bp_graph(
            &mut graph,
            &mut coverage_cut,
            &output_file_name,
            &sources,
            &reverse_sources,
            &mut read_to_unitig,
            &mut bubble_labels,
            options,
            reads,
            &mut trio_flags,
        );
    }

    // The reference's final `else` clears `HA_F_PARTITION` defensively
    // (`Overlaps.cpp` 39660) and writes the primary and alternate contig
    // graphs from the same, now relabelled, read graph.
    options.flags &= !FLAG_PARTITION;
    crate::assembly::bp_graph::output_contig_graph_primary(
        &mut graph,
        &mut coverage_cut,
        &output_file_name,
        &sources,
        &reverse_sources,
        &mut read_to_unitig,
        &mut bubble_labels,
        options,
        reads,
        &mut trio_flags,
    )?;
    crate::assembly::bp_graph::output_contig_graph_alternative(
        &mut graph,
        &coverage_cut,
        &output_file_name,
        &sources,
        &read_to_unitig,
        reads,
        &trio_flags,
        options,
    )?;

    // The last statement of `build_string_graph_without_clean`
    // (`Overlaps.cpp` 39673), reached by every dispatch branch.
    eprintln!(
        "Inconsistency threshold for low-quality regions in BED files: {}%",
        options.bed_inconsistency_rate
    );

    Ok(())
}

/// `ul_clean_gfa`: the cleaning loop, the post-rescue pass and the read
/// unitig output.
#[allow(clippy::too_many_arguments)]
fn ul_clean_gfa(
    options: &mut AssemblyOptions,
    reads: &ReadStore,
    trio_flags: &mut Vec<u8>,
    graph: &mut Asg,
    sources: &[PafList],
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    coverage_cut: &mut [CoverageCut],
    bubble_labels: &mut BubbleLabel,
    output_file_name: &str,
    min_thres: i64,
    telo_end: Option<&TeloEnd>,
    dump_directory: Option<&std::path::Path>,
) -> Result<(), String> {
    let maximum_hang = options.maximum_hang_length;
    let minimum_overlap = options.minimum_overlap_length;
    let maximum_hang_rate = options.maximum_hang_rate;
    let bed_inconsistency_rate = options.bed_inconsistency_rate as u32;
    let clean_round = options.cleaning_round_count;
    let min_drop = f64::from(options.minimum_drop_rate);
    let max_drop = f64::from(options.maximum_drop_rate);
    let ou_drop_rate = 0.6f32;
    let max_tip = options.maximum_short_tip_read_count as u32;
    let max_tip_i32 = max_tip as i32;
    let large_bubble_population_size = options.large_bubble_population_size as u64;
    let maximum_short_tip_read_count = options.maximum_short_tip_read_count;
    // `clean_graph` passes `ha_opt_triobin(&asm_opt)` into `ul_clean_gfa`
    // (`Overlaps.cpp` 39447); the round-loop cutters consume it.
    let is_trio = options.trio_binning_requested();

    let step = if clean_round == 1 {
        max_drop
    } else {
        (max_drop - min_drop) / f64::from(clean_round - 1)
    };
    let mut drop = min_drop;

    let mut buffer: Vec<u64> = Vec::new();
    let vertex_count = graph.sequences.len() * 2;

    // `asg_arc_cut_tips(sg, max_tip, &bu, is_ou, is_ou?rI:NULL, uopt->te)`.
    asg_arc_cut_tips(graph, max_tip, &mut buffer, false, None, telo_end);
    probe::try_dump_stage(
        dump_directory,
        "d03_after_cut_tips",
        Some(graph),
        sources,
        reverse_sources,
        Some(coverage_cut),
        read_to_unitig,
        &trio_flags,
    );

    for round in 0..clean_round {
        if drop > max_drop {
            drop = max_drop;
        }

        asg_iterative_semi_circ(graph, sources, &mut buffer, max_tip, true, telo_end);
        asg_arc_identify_simple_bubbles_multi(graph, bubble_labels, true, min_thres);
        asg_arc_cut_chimeric(graph, sources, &mut buffer, u32::MAX, telo_end);

        asg_arc_cut_tips(graph, max_tip, &mut buffer, false, None, telo_end);
        asg_arc_identify_simple_bubbles_multi(graph, bubble_labels, false, min_thres);
        asg_arc_cut_inexact(
            graph,
            sources,
            &mut buffer,
            max_tip_i32,
            false,
            is_trio,
            0,
            ou_drop_rate,
            &trio_flags,
        );
        asg_arc_cut_tips(graph, max_tip, &mut buffer, false, None, telo_end);

        asg_arc_identify_simple_bubbles_multi(graph, bubble_labels, true, min_thres);
        asg_arc_cut_length(
            graph,
            &mut buffer,
            max_tip_i32,
            drop as f32,
            ou_drop_rate,
            false,
            is_trio,
            true,
            0,
            1,
            None,
            None,
            None,
            &trio_flags,
        );
        asg_arc_cut_tips(graph, max_tip, &mut buffer, false, None, telo_end);

        asg_arc_identify_simple_bubbles_multi(graph, bubble_labels, true, min_thres);
        // The reference allocates a fresh `buf_t` per call, so the workspace
        // is not carried between rounds.
        let mut workspace = BubbleBuffer {
            info: vec![BubbleInfo::default(); vertex_count],
            ..BubbleBuffer::default()
        };
        asg_arc_cut_bub_links(
            graph,
            &mut buffer,
            HARD_OL_DROP,
            HARD_OL_SEC_DROP,
            HARD_OU_DROP,
            false,
            large_bubble_population_size,
            Some(reverse_sources),
            Some(read_to_unitig),
            max_tip_i32,
            &mut workspace,
        );

        asg_arc_identify_simple_bubbles_multi(graph, bubble_labels, true, min_thres);
        asg_arc_cut_complex_bub_links(
            graph,
            &mut buffer,
            HARD_OL_DROP,
            HARD_OU_DROP,
            false,
            bubble_labels,
        );
        asg_arc_cut_tips(graph, max_tip, &mut buffer, false, None, telo_end);

        probe::try_dump_stage(
            dump_directory,
            &format!("d{:02}_after_round_{}", 4 + round, round),
            Some(graph),
            sources,
            reverse_sources,
            Some(coverage_cut),
            read_to_unitig,
            &trio_flags,
        );

        drop += step;
    }

    asg_iterative_semi_circ(graph, sources, &mut buffer, max_tip, true, telo_end);
    asg_arc_identify_simple_bubbles_multi(graph, bubble_labels, false, min_thres);
    asg_cut_large_indel(
        graph,
        &mut buffer,
        max_tip_i32,
        HARD_OU_DROP,
        false,
        0,
    );
    asg_arc_cut_tips(graph, max_tip, &mut buffer, false, None, telo_end);

    probe::try_dump_stage(
        dump_directory,
        "d08_after_large_indel",
        Some(graph),
        sources,
        reverse_sources,
        Some(coverage_cut),
        read_to_unitig,
        &trio_flags,
    );

    asg_arc_identify_simple_bubbles_multi(graph, bubble_labels, false, min_thres);
    // The two post-loop passes pass `0/**is_trio**/` in the reference
    // (`gfa_ut.cpp` 3157/3163); only the round loop is trio aware.
    asg_arc_cut_length(
        graph,
        &mut buffer,
        max_tip_i32,
        HARD_ORTHOLOGY_DROP,
        ou_drop_rate,
        false,
        false,
        false,
        0,
        1,
        Some(reverse_sources),
        Some(read_to_unitig),
        None,
        &trio_flags,
    );
    asg_arc_cut_tips(graph, max_tip, &mut buffer, false, None, telo_end);

    asg_arc_identify_simple_bubbles_multi(graph, bubble_labels, false, min_thres);
    asg_arc_cut_length(
        graph,
        &mut buffer,
        max_tip_i32,
        options.minimum_drop_rate,
        ou_drop_rate,
        false,
        false,
        false,
        0,
        1,
        Some(reverse_sources),
        Some(read_to_unitig),
        Some(2000),
        &trio_flags,
    );
    asg_arc_cut_tips(graph, max_tip, &mut buffer, false, None, telo_end);

    probe::try_dump_stage(
        dump_directory,
        "d09_after_cut_length_final",
        Some(graph),
        sources,
        reverse_sources,
        Some(coverage_cut),
        read_to_unitig,
        &trio_flags,
    );

    asg_cut_semi_circ(graph, LIM_LEN, true);
    probe::try_dump_stage(
        dump_directory,
        "d10_after_semi_circ",
        Some(graph),
        sources,
        reverse_sources,
        Some(coverage_cut),
        read_to_unitig,
        &trio_flags,
    );

    post_rescue(
        options,
        trio_flags,
        graph,
        sources,
        reverse_sources,
        read_to_unitig,
        coverage_cut,
        bubble_labels,
        min_thres,
        dump_directory,
    );

    ug_ext_gfa(
        graph,
        sources,
        read_to_unitig,
        (maximum_short_tip_read_count * 2) as u32,
        maximum_hang,
        minimum_overlap,
        maximum_hang_rate,
        UG_EXT_LEN,
    );
    probe::try_dump_stage(
        dump_directory,
        "d17_after_ug_ext",
        Some(graph),
        sources,
        reverse_sources,
        Some(coverage_cut),
        read_to_unitig,
        &trio_flags,
    );

    output_unitig_graph(
        graph,
        coverage_cut,
        output_file_name,
        sources,
        read_to_unitig,
        maximum_hang,
        minimum_overlap,
        reads,
        maximum_hang_rate,
        bed_inconsistency_rate,
    )
    .map_err(|error| format!("failed to write the read unitig GFA: {}", error))?;

    probe::try_dump_stage(
        dump_directory,
        "d18_after_output_unitig",
        Some(graph),
        sources,
        reverse_sources,
        Some(coverage_cut),
        read_to_unitig,
        &trio_flags,
    );

    // `flat_soma_v(sg, src, rI)` closes `ul_clean_gfa`: it collapses the
    // heterozygous bubbles of the read graph, which is what turns the raw
    // unitigs into the primary ones that `output_contig_graph_primary_pre`
    // later writes.
    flat_soma_v(options, graph, sources, read_to_unitig);
    probe::try_dump_stage(
        dump_directory,
        "d19_after_flat_soma",
        Some(graph),
        sources,
        reverse_sources,
        Some(coverage_cut),
        read_to_unitig,
        &trio_flags,
    );

    Ok(())
}

/// `post_rescue`: the rescue passes that run at the end of `ul_clean_gfa`.
#[allow(clippy::too_many_arguments)]
fn post_rescue(
    options: &mut AssemblyOptions,
    trio_flags: &mut Vec<u8>,
    graph: &mut Asg,
    sources: &[PafList],
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    coverage_cut: &mut [CoverageCut],
    bubble_labels: &mut BubbleLabel,
    min_thres: i64,
    dump_directory: Option<&std::path::Path>,
) {
    let parameters = RescueParameters {
        maximum_hang: options.maximum_hang_length,
        maximum_hang_rate: options.maximum_hang_rate,
        minimum_overlap: options.minimum_overlap_length,
        minimum_length_threshold: min_thres,
    };
    let trio_labels: &[u8] = trio_flags;

    rescue_contained_reads_aggressive(
        None,
        graph,
        sources,
        coverage_cut,
        read_to_unitig,
        trio_labels,
        parameters,
        10,
        true,
        false,
        None,
        None,
        bubble_labels,
    );
    probe::try_dump_stage(
        dump_directory,
        "d12_after_rescue_contained",
        Some(graph),
        sources,
        reverse_sources,
        Some(coverage_cut),
        read_to_unitig,
        trio_labels,
    );

    rescue_missing_overlaps_aggressive(
        None,
        graph,
        sources,
        coverage_cut,
        read_to_unitig,
        trio_labels,
        parameters,
        true,
        false,
        None,
        bubble_labels,
    );
    probe::try_dump_stage(
        dump_directory,
        "d13_after_rescue_missing",
        Some(graph),
        sources,
        reverse_sources,
        Some(coverage_cut),
        read_to_unitig,
        trio_labels,
    );

    rescue_missing_overlaps_backward(
        None,
        graph,
        sources,
        coverage_cut,
        read_to_unitig,
        trio_labels,
        parameters,
        10,
        true,
        false,
        bubble_labels,
    );
    probe::try_dump_stage(
        dump_directory,
        "d14_after_rescue_backward",
        Some(graph),
        sources,
        reverse_sources,
        Some(coverage_cut),
        read_to_unitig,
        trio_labels,
    );

    // `cmk` is NULL on the HiFi path, so `rescue_chimeric_reads_aggressive`
    // is not called.

    set_hom_global_coverage(
        options,
        graph,
        coverage_cut,
        sources,
        read_to_unitig,
        parameters.maximum_hang,
        parameters.minimum_overlap,
    );
    probe::try_dump_stage(
        dump_directory,
        "d15_after_set_hom_cov",
        Some(graph),
        sources,
        reverse_sources,
        Some(coverage_cut),
        read_to_unitig,
        trio_labels,
    );

    // `rescue_bubble_by_chain(sg, uopt->coverage_cut, src, rev,
    // (asm_opt.max_short_tip*2), 0.15, 3, rI, 0.05, 0.9, uopt->max_hang,
    // uopt->min_ovlp, 10, uopt->gap_fuzz, b_mask_t, 0, NULL)` — the final
    // pass of `post_rescue`. The reference mutates `R_INF.trio_flag` through
    // it; the label writes of the live trio path are read-only, so a local
    // copy stands in for the global (the trio gate validates this).
    // The reference runs this pass unconditionally; on a non-trio corpus it
    // is a no-op (probe d15 == d16), so the port only enters it when the
    // trio front end is active until its remaining dependencies
    // (the Purge_Dups hap-alignment candidate scan) are ported.
    if !options.trio_binning_requested() {
        return;
    }
    if let Err(error) = crate::assembly::bubble_chain::rescue_bubble_by_chain(
        graph,
        coverage_cut,
        sources,
        reverse_sources,
        read_to_unitig,
        bubble_labels,
        options,
        trio_flags,
        i64::from(options.maximum_short_tip_read_count) * 2,
        0.15,
        3,
        0.05,
        0.9,
        options.maximum_hang_length,
        options.minimum_overlap_length,
        10,
        min_thres,
    ) {
        eprintln!("[ERROR] rescue_bubble_by_chain failed: {}", error);
    }
    probe::try_dump_stage(
        dump_directory,
        "d16_after_bubble_chain",
        Some(graph),
        sources,
        reverse_sources,
        Some(coverage_cut),
        read_to_unitig,
        trio_flags,
    );
}

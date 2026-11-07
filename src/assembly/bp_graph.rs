//! `output_bp_graph` and its p_ctg-side callees `adjust_utg_by_primary` and
//! `print_utg` from `Overlaps.cpp`.
//!
//! The default HiFi pipeline reaches this driver through the
//! `HA_F_PARTITION && purge_level_primary > 0` branch of `clean_graph`. It
//! clones the read graph and the primary unitig graph, cleans the clone with
//! the coverage-based purge pass, and writes `bp.p_ctg.*`; the trio output
//! stage that follows writes the haplotype graphs.
//!
//! The haplotype half of `output_bp_graph` is `output_trio_graph_joint`
//! (`Overlaps.cpp` 23244): it builds a phase unitig graph per haplotype with
//! `ma_ug_gen_phase`, cleans each with [`adjust_utg_by_trio`], deduplicates
//! exact duplicates across the two phase graphs, appends the missed unitigs
//! and writes `bp.hap1/hap2.p_ctg.*` through [`output_hap_graph`].

use std::fs::File;
use std::io::{BufWriter, Write};

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::bubble::BubbleLabel;
use crate::assembly::hamming::reduce_hamming_error_adv;
use crate::assembly::hap_cov::{init_hap_cov_t, HapCov};
use crate::assembly::hit::{CoverageCut, ReadToUnitig, ALTER_LABEL, PRIMARY_LABEL};
use crate::assembly::purge_dups::{purge_dups, HOM_PEAK_RATE};
use crate::assembly::rescue::{
    rescue_contained_reads_aggressive, rescue_missing_overlaps_aggressive, RescueParameters,
};
use crate::assembly::trio_clean::{
    asg_seq_drop, clean_primary_untig_graph, clean_trio_untig_graph, delete_useless_trio_nodes,
    force_trio_clean, if_primary_unitig, set_drop_trio_flag, update_hap_label,
    update_unitig_graph, DOUBLE_CHECK_THRES, TRIO_THRES,
};
use crate::assembly::trans_clean::{clean_u_trans_t_idx, set_trio_flag_by_cov};
use crate::assembly::unitig::{
    adjust_utg_advance, copy_read_graph, copy_untig_graph, delete_useless_nodes,
    discard_small_ctg, hic_clean, ma_ug_gen_phase, ma_ug_gen_primary, ma_ug_print,
    ma_ug_print_bed, ma_ug_print_simple, ma_ug_seq, recover_utg_by_coverage, renew_utg,
    UnitigGraph, FATHER, MOTHER, PHASE_SEF, PHASE_SEP, PHASE_SEP_RATE,
};
use crate::command_line::options::{AssemblyOptions, FLAG_BAN_POST_JOIN};
use crate::correction::ec_driver::PafList;
use crate::reads::read_store::ReadStore;

/// `adjust_utg_primary`'s `tipsLen`: `clean_graph` passes
/// `asm_opt.max_short_tip * 2` to `output_bp_graph`.
fn primary_tips_length(options: &AssemblyOptions) -> i64 {
    i64::from(options.maximum_short_tip_read_count * 2)
}

/// `output_bp_graph`: the default partition output driver.
///
/// `read_graph` is the cleaned read graph `sg`, `coverage_cut` the per-read
/// coverage cuts and `read_to_unitig` the reference's `ruIndex`; the function
/// proceeds through `hic_clean`, `adjust_utg_by_primary` and `print_utg` and
/// then reaches the haplotype stage, which is not ported yet. The primary
/// contig files are complete when this returns `Err`.
#[allow(clippy::too_many_arguments)]
pub fn output_bp_graph(
    read_graph: &mut Asg,
    coverage_cut: &mut [CoverageCut],
    output_file_name: &str,
    sources: &[PafList],
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    bubble_labels: &mut BubbleLabel,
    options: &mut AssemblyOptions,
    reads: &ReadStore,
    trio_flags: &mut [u8],
) -> Result<(), String> {
    hic_clean(read_graph);

    let unitig_graph = ma_ug_gen_primary(read_graph, PRIMARY_LABEL);
    let mut copy_read_graph = copy_read_graph(read_graph);
    let mut copy_unitig_graph = copy_untig_graph(&unitig_graph);

    // The reference allocates `new_rtg_edges` empty here; only the polish
    // passes can append to it, and both output drivers pass `E = NULL`.
    let new_rtg_edges: Vec<AsgArc> = Vec::new();

    let mut cov = adjust_utg_by_primary(
        &mut copy_unitig_graph,
        &mut copy_read_graph,
        coverage_cut,
        sources,
        reverse_sources,
        read_to_unitig,
        primary_tips_length(options),
        0.15,
        3,
        0.05,
        0.9,
        options.maximum_hang_length,
        options.minimum_overlap_length,
        &new_rtg_edges,
        Some(bubble_labels),
        options,
        trio_flags,
        true,
        true,
    )?;

    print_utg(
        &mut copy_unitig_graph,
        &mut copy_read_graph,
        coverage_cut,
        output_file_name,
        sources,
        read_to_unitig,
        options.maximum_hang_length,
        options.minimum_overlap_length,
        &new_rtg_edges,
        reads,
        trio_flags,
        options,
    )?;

    // The trio half of the driver: clean the collected transitions against the
    // original primary unitig graph and derive the coverage-based trio flags.
    let t_ch = cov
        .t_ch
        .as_mut()
        .ok_or_else(|| "output_bp_graph requires the transition chain".to_string())?;
    clean_u_trans_t_idx(
        &mut t_ch.k_trans,
        &mut t_ch.k_trans_idx,
        &unitig_graph.unitigs,
        &unitig_graph.graph,
        read_graph,
    );
    set_trio_flag_by_cov(
        &unitig_graph.unitigs,
        read_graph,
        &cov.is_r_het,
        &t_ch.k_trans,
        &t_ch.k_trans_idx,
        trio_flags,
    );

    // `output_trio_graph_joint(sg, coverage_cut, output_file_name, sources,
    // reverse_sources, (asm_opt.max_short_tip*2), 0.15, 3, ruIndex, 0.05, 0.9,
    // max_hang, min_ovlp, gap_fuzz, b_mask_t, NULL, NULL, opt)`.
    output_trio_graph_joint(
        read_graph,
        coverage_cut,
        output_file_name,
        sources,
        reverse_sources,
        primary_tips_length(options),
        0.15,
        3,
        read_to_unitig,
        0.05,
        0.9,
        options.maximum_hang_length,
        options.minimum_overlap_length,
        options.gap_fuzz,
        bubble_labels,
        options,
        trio_flags,
        reads,
    )
}

/// `adjust_utg_by_primary`: cleans and purges the copied primary unitig graph.
///
/// The function returns the `hap_cov_t` the reference hands back through its
/// `i_cov` out-parameter; `print_utg` and the later trio stages read it.
#[allow(clippy::too_many_arguments)]
pub fn adjust_utg_by_primary(
    unitig_graph: &mut UnitigGraph,
    read_graph: &mut Asg,
    coverage_cut: &mut [CoverageCut],
    sources: &[PafList],
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    tips_len: i64,
    tip_drop_ratio: f32,
    stops_threshold: u32,
    chimeric_rate: f32,
    drop_ratio: f32,
    max_hang: i32,
    min_ovlp: i32,
    new_rtg_edges: &[AsgArc],
    mut bubble_labels: Option<&mut BubbleLabel>,
    options: &mut AssemblyOptions,
    trio_flags: &mut [u8],
    collect_p_trans: bool,
    collect_p_trans_f: bool,
) -> Result<HapCov, String> {
    // `(asm_opt.purge_level_primary > 0 || i_cov) ? 1 : 0`: every caller of
    // this port asks for the coverage object, so `i_cov` is non-null.
    let mut cov = init_hap_cov_t(
        unitig_graph,
        read_graph,
        sources,
        read_to_unitig,
        coverage_cut,
        options,
        true,
    );

    // The reference aliases `cov->t_ch->ir_het` with `cov->is_r_het`; the
    // Rust `TransChain` has no separate buffer.

    adjust_utg_advance(
        read_graph,
        unitig_graph,
        reverse_sources,
        read_to_unitig,
        bubble_labels.as_deref_mut(),
        &cov.is_r_het,
        options.purge_similarity_threshold,
    );

    for v in 0..unitig_graph.graph.sequences.len() {
        if unitig_graph.graph.sequences[v].deleted {
            continue;
        }
        unitig_graph.graph.sequences[v].coverage_class = PRIMARY_LABEL;
        unitig_graph.unitigs[v].start = unitig_graph.unitigs[v].reads.len() as u32;
    }

    clean_primary_untig_graph(
        unitig_graph,
        read_graph,
        reverse_sources,
        coverage_cut,
        tips_len,
        tip_drop_ratio,
        stops_threshold,
        read_to_unitig,
        chimeric_rate,
        false,
        drop_ratio,
        &mut cov,
        trio_flags,
        options.maximum_short_tip_read_count,
        options.purge_similarity_threshold,
        options.maximum_hang_rate,
    );
    delete_useless_nodes(unitig_graph);
    renew_utg(unitig_graph, read_graph, new_rtg_edges);

    if options.purge_level_primary > 0 {
        // `just_contain = (purge_level_primary == 1)`.
        let just_contain = options.purge_level_primary == 1;
        let collect = cov.t_ch.is_some() && collect_p_trans;
        purge_dups(
            unitig_graph,
            read_graph,
            coverage_cut,
            sources,
            reverse_sources,
            read_to_unitig,
            options.purge_similarity_threshold,
            options.purge_overlap_length as u32,
            max_hang,
            min_ovlp,
            drop_ratio,
            just_contain,
            false,
            &mut cov,
            collect,
            collect_p_trans_f,
            options,
            trio_flags,
        )?;
        delete_useless_nodes(unitig_graph);
        renew_utg(unitig_graph, read_graph, new_rtg_edges);
    }

    if options.flags & FLAG_BAN_POST_JOIN == 0 {
        let parameters = RescueParameters {
            maximum_hang: max_hang,
            maximum_hang_rate: options.maximum_hang_rate,
            minimum_overlap: min_ovlp,
            minimum_length_threshold: i64::from(options.maximum_short_tip_read_count + 1),
        };
        let Some(labels) = bubble_labels.as_deref_mut() else {
            return Err(
                "adjust_utg_by_primary requires the bubble labels for its rescue passes"
                    .to_string(),
            );
        };

        rescue_missing_overlaps_aggressive(
            Some(unitig_graph),
            read_graph,
            sources,
            coverage_cut,
            read_to_unitig,
            trio_flags,
            parameters,
            false,
            true,
            None,
            &mut *labels,
        );
        renew_utg(unitig_graph, read_graph, new_rtg_edges);
        rescue_contained_reads_aggressive(
            Some(unitig_graph),
            read_graph,
            sources,
            coverage_cut,
            read_to_unitig,
            trio_flags,
            parameters,
            10,
            false,
            true,
            None,
            None,
            &mut *labels,
        );
        renew_utg(unitig_graph, read_graph, new_rtg_edges);
    }

    if options.maximum_contig_tip_read_count > 0 {
        discard_small_ctg(
            unitig_graph,
            read_graph,
            new_rtg_edges,
            options.maximum_contig_tip_read_count,
            false,
            true,
        );
    }

    // Every read is alternate unless a surviving primary unitig claims it;
    // the unclaimed reads are then dropped from the copied read graph.
    for sequence in &mut read_graph.sequences {
        sequence.coverage_class = ALTER_LABEL;
    }
    for v in 0..unitig_graph.graph.sequences.len() {
        if unitig_graph.graph.sequences[v].deleted {
            continue;
        }
        if unitig_graph.graph.sequences[v].coverage_class == ALTER_LABEL {
            continue;
        }
        let unitig = &unitig_graph.unitigs[v];
        if unitig.reads.is_empty() {
            continue;
        }
        let class = unitig_graph.graph.sequences[v].coverage_class;
        for &packed in &unitig.reads {
            read_graph.sequences[(packed >> 33) as usize].coverage_class = class;
        }
    }
    for v in 0..read_graph.sequences.len() {
        if read_graph.sequences[v].coverage_class == ALTER_LABEL {
            asg_seq_drop(read_graph, v as u32);
        }
    }

    if options.recover_alternate_contig_coverage_minimum == -1024 {
        let maximum = if options.global_homozygous_coverage_is_set != 0 {
            f64::from(options.global_homozygous_coverage)
        } else {
            f64::from(options.global_homozygous_coverage) / HOM_PEAK_RATE
        };
        options.recover_alternate_contig_coverage_maximum = maximum as i64;
        options.recover_alternate_contig_coverage_minimum = (maximum * 0.85) as i64;
        options.recover_alternate_contig_coverage_maximum = i64::from(i32::MAX);
    }

    if options.recover_alternate_contig_coverage_maximum != i64::from(i32::MAX) {
        eprintln!(
            "[M::adjust_utg_by_primary] primary contig coverage range: [{}, {}]",
            options.recover_alternate_contig_coverage_minimum,
            options.recover_alternate_contig_coverage_maximum
        );
    } else {
        eprintln!(
            "[M::adjust_utg_by_primary] primary contig coverage range: [{}, infinity]",
            options.recover_alternate_contig_coverage_minimum
        );
    }

    let is_r_het = if cov.t_ch.is_some() {
        Some(&mut cov.is_r_het[..])
    } else {
        None
    };
    recover_utg_by_coverage(
        unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        read_to_unitig,
        options.recover_alternate_contig_coverage_minimum as i32,
        options.recover_alternate_contig_coverage_maximum as i32,
        is_r_het,
    );

    Ok(cov)
}

/// `print_utg`: materialises the unitig sequences and writes `p_ctg.*`.
#[allow(clippy::too_many_arguments)]
pub fn print_utg(
    unitig_graph: &mut UnitigGraph,
    read_graph: &mut Asg,
    coverage_cut: &[CoverageCut],
    output_file_name: &str,
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    max_hang: i32,
    min_ovlp: i32,
    new_rtg_edges: &[AsgArc],
    reads: &ReadStore,
    trio_labels: &[u8],
    options: &AssemblyOptions,
) -> Result<(), String> {
    // `break_ug_contig` for `b_low_cov`/`b_high_cov` is not ported; both
    // switches are off in the default configuration.
    if options.breakpoint_low_coverage > 0 || options.breakpoint_high_coverage > 0 {
        return Err(
            "print_utg needs break_ug_contig (b_low_cov/b_high_cov), which is not ported yet"
                .to_string(),
        );
    }

    ma_ug_seq(
        unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        new_rtg_edges,
        None,
        max_hang,
        min_ovlp,
        true,
        reads,
        options.maximum_hang_rate,
    );

    let gfa_name = format!("{}.p_ctg.gfa", output_file_name);
    eprintln!("Writing {} to disk... ", gfa_name);
    let mut output = BufWriter::new(
        File::create(&gfa_name).map_err(|error| format!("failed to create {}: {}", gfa_name, error))?,
    );
    ma_ug_print(
        unitig_graph,
        read_graph,
        Some(coverage_cut),
        Some(sources),
        Some(read_to_unitig),
        "ptg",
        &mut output,
        reads,
        trio_labels,
    )
    .map_err(|error| format!("failed to write {}: {}", gfa_name, error))?;
    output
        .flush()
        .map_err(|error| format!("failed to write {}: {}", gfa_name, error))?;

    let noseq_name = format!("{}.p_ctg.noseq.gfa", output_file_name);
    let mut output = BufWriter::new(
        File::create(&noseq_name)
            .map_err(|error| format!("failed to create {}: {}", noseq_name, error))?,
    );
    ma_ug_print_simple(
        unitig_graph,
        read_graph,
        Some(coverage_cut),
        Some(sources),
        Some(read_to_unitig),
        "ptg",
        &mut output,
        reads,
        trio_labels,
    )
    .map_err(|error| format!("failed to write {}: {}", noseq_name, error))?;
    output
        .flush()
        .map_err(|error| format!("failed to write {}: {}", noseq_name, error))?;

    if options.bed_inconsistency_rate != 0 {
        let bed_name = format!("{}.p_ctg.lowQ.bed", output_file_name);
        let mut output = BufWriter::new(
            File::create(&bed_name)
                .map_err(|error| format!("failed to create {}: {}", bed_name, error))?,
        );
        ma_ug_print_bed(
            unitig_graph,
            read_graph,
            coverage_cut,
            sources,
            new_rtg_edges,
            max_hang,
            min_ovlp,
            options.bed_inconsistency_rate as u32,
            "ptg",
            &mut output,
            reads,
            options.maximum_hang_rate,
        )
        .map_err(|error| format!("failed to write {}: {}", bed_name, error))?;
        output
            .flush()
            .map_err(|error| format!("failed to write {}: {}", bed_name, error))?;
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// `output_trio_graph_joint` (`Overlaps.cpp` 23244) and its callees.
// ---------------------------------------------------------------------------

/// `PHASE_MISS_LEN` / `PHASE_MISS_N`: the length and read-count floors a
/// missed unitig must clear to be appended to a haplotype graph.
const PHASE_MISS_LEN: u64 = 1_000_000;
const PHASE_MISS_N: u64 = 8;

/// `FINAL_DOUBLE_CHECK_THRES`: the `double_check_rate` of the final
/// [`update_unitig_graph`] pass (`Overlaps.h` 31).
const FINAL_DOUBLE_CHECK_THRES: f64 = 0.2;

/// `ma_ug_gen_phase`'s `pscut`: the phase read-occurrence floor.
fn phase_cut_threshold(options: &AssemblyOptions) -> u64 {
    let base = if options.global_homozygous_coverage_is_set != 0 {
        f64::from(options.global_homozygous_coverage)
    } else {
        f64::from(options.global_homozygous_coverage) / HOM_PEAK_RATE
    };
    ((base * PHASE_SEF) as u64).max(PHASE_SEP)
}

/// `update_recover_atg_cov`: fills the recover-coverage window from the
/// homozygous peak when it is still unset.
fn update_recover_atg_cov(options: &mut AssemblyOptions) {
    if options.recover_alternate_contig_coverage_minimum == -1024 {
        let maximum = if options.global_homozygous_coverage_is_set != 0 {
            f64::from(options.global_homozygous_coverage)
        } else {
            f64::from(options.global_homozygous_coverage) / HOM_PEAK_RATE
        };
        options.recover_alternate_contig_coverage_maximum = maximum as i64;
        options.recover_alternate_contig_coverage_minimum = (maximum * 0.85) as i64;
        options.recover_alternate_contig_coverage_maximum = i64::from(i32::MAX);
    }
}

/// `delete_ug_node`: marks the unitig vertex deleted and frees its reads.
///
/// `length` is deliberately preserved, matching the reference: the callers
/// read it after the deletion to accumulate the removed base count.
fn delete_ug_node(unitig_graph: &mut UnitigGraph, node: usize) {
    unitig_graph.graph.delete_sequence(node as u32);
    if !unitig_graph.unitigs[node].reads.is_empty() {
        unitig_graph.unitigs[node].reads = Vec::new();
    }
}

/// `adjust_utg_by_trio`: the per-haplotype counterpart of
/// [`adjust_utg_by_primary`] (`Overlaps.cpp` 20995).
///
/// The phase unitig graph is cleaned against the trio flags, the tangles and
/// the primary coverage; every read of a surviving unitig that has no trio
/// assignment yet is finally marked [`DROP`].
#[allow(clippy::too_many_arguments)]
pub fn adjust_utg_by_trio(
    unitig_graph: &mut UnitigGraph,
    read_graph: &mut Asg,
    flag: u32,
    drop_rate: f32,
    sources: &[PafList],
    reverse_sources: &[PafList],
    coverage_cut: &mut [CoverageCut],
    tips_len: i64,
    tip_drop_ratio: f32,
    stops_threshold: u32,
    read_to_unitig: &mut ReadToUnitig,
    chimeric_rate: f32,
    drop_ratio: f32,
    max_hang: i32,
    min_ovlp: i32,
    gap_fuzz: i32,
    new_rtg_edges: &mut Vec<AsgArc>,
    bubble_labels: &mut BubbleLabel,
    options: &mut AssemblyOptions,
    trio_flags: &mut [u8],
) -> Result<HapCov, String> {
    // `init_hap_cov_t(..., asm_opt.purge_level_trio>0?1:0)`.
    let mut cov = init_hap_cov_t(
        unitig_graph,
        read_graph,
        sources,
        read_to_unitig,
        coverage_cut,
        options,
        options.purge_level_trio > 0,
    );

    update_recover_atg_cov(options);
    if options.recover_alternate_contig_coverage_maximum != i64::from(i32::MAX) {
        eprintln!(
            "[M::adjust_utg_by_trio] primary contig coverage range: [{}, {}]",
            options.recover_alternate_contig_coverage_minimum,
            options.recover_alternate_contig_coverage_maximum
        );
    } else {
        eprintln!(
            "[M::adjust_utg_by_trio] primary contig coverage range: [{}, infinity]",
            options.recover_alternate_contig_coverage_minimum
        );
    }

    adjust_utg_advance(
        read_graph,
        unitig_graph,
        reverse_sources,
        read_to_unitig,
        Some(&mut *bubble_labels),
        &cov.is_r_het,
        options.purge_similarity_threshold,
    );

    update_unitig_graph(
        unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        reverse_sources,
        read_to_unitig,
        &cov.is_r_het,
        false,
        DOUBLE_CHECK_THRES,
        flag,
        drop_rate,
        trio_flags,
        options.purge_similarity_threshold,
        options,
    );

    for v in 0..unitig_graph.graph.sequences.len() {
        if unitig_graph.graph.sequences[v].deleted {
            continue;
        }
        unitig_graph.graph.sequences[v].coverage_class = PRIMARY_LABEL;
        unitig_graph.unitigs[v].start = unitig_graph.unitigs[v].reads.len() as u32;
    }

    clean_trio_untig_graph(
        unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        reverse_sources,
        tips_len,
        tip_drop_ratio,
        stops_threshold,
        read_to_unitig,
        chimeric_rate,
        false,
        drop_ratio,
        flag,
        drop_rate,
        max_hang,
        min_ovlp,
        gap_fuzz,
        &mut cov,
        new_rtg_edges,
        trio_flags,
        options.purge_similarity_threshold,
        options.maximum_hang_rate,
        options,
    );

    delete_useless_trio_nodes(
        unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        read_to_unitig,
        flag,
        0.8,
        0.15,
        16,
        trio_flags,
        options,
    );

    update_hap_label(Some(unitig_graph), read_graph);

    update_unitig_graph(
        unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        reverse_sources,
        read_to_unitig,
        &cov.is_r_het,
        false,
        DOUBLE_CHECK_THRES,
        flag,
        drop_rate,
        trio_flags,
        options.purge_similarity_threshold,
        options,
    );

    force_trio_clean(
        unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        read_to_unitig,
        flag,
        0.55,
        0.01,
        5,
        trio_flags,
        options,
    );

    renew_utg(unitig_graph, read_graph, new_rtg_edges);

    if options.flags & FLAG_BAN_POST_JOIN == 0 {
        let parameters = RescueParameters {
            maximum_hang: max_hang,
            maximum_hang_rate: options.maximum_hang_rate,
            minimum_overlap: min_ovlp,
            minimum_length_threshold: i64::from(options.maximum_short_tip_read_count + 1),
        };
        rescue_missing_overlaps_aggressive(
            Some(unitig_graph),
            read_graph,
            sources,
            coverage_cut,
            read_to_unitig,
            trio_flags,
            parameters,
            false,
            true,
            None,
            &mut *bubble_labels,
        );
        renew_utg(unitig_graph, read_graph, new_rtg_edges);

        rescue_contained_reads_aggressive(
            Some(unitig_graph),
            read_graph,
            sources,
            coverage_cut,
            read_to_unitig,
            trio_flags,
            parameters,
            10,
            false,
            true,
            None,
            None,
            &mut *bubble_labels,
        );
        renew_utg(unitig_graph, read_graph, new_rtg_edges);
    }

    update_unitig_graph(
        unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        reverse_sources,
        read_to_unitig,
        &cov.is_r_het,
        true,
        FINAL_DOUBLE_CHECK_THRES,
        flag,
        drop_rate,
        trio_flags,
        options.purge_similarity_threshold,
        options,
    );

    update_hap_label(None, read_graph);

    renew_utg(unitig_graph, read_graph, new_rtg_edges);

    delete_useless_trio_nodes(
        unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        read_to_unitig,
        flag,
        0.8,
        0.15,
        16,
        trio_flags,
        options,
    );

    if options.purge_level_trio == 1 {
        // `purge_dups(*ug, read_g, ..., drop_ratio, 1, 0, cov, 0, 0)`.
        purge_dups(
            unitig_graph,
            read_graph,
            coverage_cut,
            sources,
            reverse_sources,
            read_to_unitig,
            options.purge_similarity_threshold,
            options.purge_overlap_length as u32,
            max_hang,
            min_ovlp,
            drop_ratio,
            true,
            false,
            &mut cov,
            false,
            false,
            options,
            trio_flags,
        )?;
        delete_useless_trio_nodes(
            unitig_graph,
            read_graph,
            coverage_cut,
            sources,
            read_to_unitig,
            flag,
            0.8,
            0.15,
            16,
            trio_flags,
            options,
        );
    }

    if options.maximum_contig_tip_read_count > 0 {
        discard_small_ctg(
            unitig_graph,
            read_graph,
            new_rtg_edges,
            options.maximum_contig_tip_read_count,
            true,
            false,
        );
    }

    set_drop_trio_flag(unitig_graph, trio_flags);
    renew_utg(unitig_graph, read_graph, new_rtg_edges);

    Ok(cov)
}

/// `output_trio_unitig_graph`: builds and cleans one haplotype's phase unitig
/// graph (`Overlaps.cpp` 21159).
///
/// Only the benchmark path (`is_bench`) is ported: the graph is returned to
/// the caller together with the re-linked edges collected during cleaning,
/// and the GFA writing happens in [`output_hap_graph`]. The `break_ug_contig`
/// and `update_dump_trio` branches are off in the default configuration.
#[allow(clippy::too_many_arguments)]
pub fn output_trio_unitig_graph(
    read_graph: &mut Asg,
    coverage_cut: &mut [CoverageCut],
    output_file_name: &str,
    flag: u32,
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
    gap_fuzz: i32,
    is_bench: bool,
    bubble_labels: &mut BubbleLabel,
    r_edges: &mut Vec<AsgArc>,
    options: &mut AssemblyOptions,
    trio_flags: &mut [u8],
    _reads: &ReadStore,
) -> Result<UnitigGraph, String> {
    let _ = output_file_name;
    if !is_bench {
        return Err(
            "output_trio_unitig_graph is only ported for the benchmark path (is_bench)"
                .to_string(),
        );
    }

    let phase_threshold = phase_cut_threshold(options);
    let mut unitig_graph = ma_ug_gen_phase(read_graph, phase_threshold as u32, PHASE_SEP_RATE, trio_flags);

    let mut new_rtg_edges: Vec<AsgArc> = Vec::new();
    adjust_utg_by_trio(
        &mut unitig_graph,
        read_graph,
        flag,
        TRIO_THRES,
        sources,
        reverse_sources,
        coverage_cut,
        tips_len,
        tip_drop_ratio,
        stops_threshold,
        read_to_unitig,
        chimeric_rate,
        drop_ratio,
        max_hang,
        min_ovlp,
        gap_fuzz,
        &mut new_rtg_edges,
        bubble_labels,
        options,
        trio_flags,
    )?;

    if options.breakpoint_low_coverage > 0 || options.breakpoint_high_coverage > 0 {
        return Err(
            "output_trio_unitig_graph needs break_ug_contig (b_low_cov/b_high_cov), which is not ported yet"
                .to_string(),
        );
    }

    // `is_bench`: hand the re-linked edges to the caller and return the graph.
    r_edges.extend_from_slice(&new_rtg_edges);
    Ok(unitig_graph)
}

/// `dedup_idx_t`: the per-read occurrence lists of a unitig graph.
///
/// `ridx[r]` is the prefix sum of the occurrences of read `r` and `ra` holds
/// one `<unitig>:32 | <position>:32` entry per occurrence. The last slot of
/// every read's span doubles as the insertion cursor of its intrusive list,
/// exactly as in the reference.
struct DedupIndex {
    ridx: Vec<u64>,
    ra: Vec<u64>,
}

impl DedupIndex {
    /// `gen_dedup_idx_t`.
    fn new(unitig_graph: &UnitigGraph, read_count: usize) -> Self {
        let mut ridx = vec![0u64; read_count + 1];
        let mut ra_count = 0usize;
        for unitig in &unitig_graph.unitigs {
            ra_count += unitig.reads.len();
            if unitig.reads.is_empty() {
                continue;
            }
            for &packed in &unitig.reads {
                ridx[(packed >> 33) as usize] += 1;
            }
        }
        let mut next = 0u64;
        for slot in ridx.iter_mut() {
            let count = *slot;
            *slot = next;
            next += count;
        }

        let mut ra = vec![u64::MAX; ra_count];
        for k in 1..ridx.len() {
            let base = ridx[k - 1] as usize;
            let len = (ridx[k] - ridx[k - 1]) as usize;
            if len > 0 {
                ra[base + len - 1] = 0;
            }
        }

        for (k, unitig) in unitig_graph.unitigs.iter().enumerate() {
            if unitig.reads.is_empty() {
                continue;
            }
            for (z, &packed) in unitig.reads.iter().enumerate() {
                let read = (packed >> 33) as usize;
                let base = ridx[read] as usize;
                let len = (ridx[read + 1] - ridx[read]) as usize;
                if len == 0 {
                    continue;
                }
                if ra[base + len - 1] == (len - 1) as u64 {
                    ra[base + len - 1] = ((k as u64) << 32) | z as u64;
                } else {
                    let cursor = ra[base + len - 1] as usize;
                    ra[base + cursor] = ((k as u64) << 32) | z as u64;
                    ra[base + len - 1] = cursor as u64 + 1;
                }
            }
        }

        Self { ridx, ra }
    }
}

/// `cal_exact_ug_o`: whether `unitig` has an exact duplicate in the graph
/// behind `index`.
///
/// A duplicate is a unitig that contains the same start and end reads with
/// the same orientation and exactly the query's read count. Returns
/// `i32::MIN` when there is none, `-1` when the query is the copy to delete
/// and the duplicate's unitig id when the duplicate should go instead.
fn cal_exact_ug_o(
    index: &DedupIndex,
    unitig_graph: &UnitigGraph,
    unitig: &MaUtg,
    flag: u32,
    trio_flags: &[u8],
) -> i64 {
    if unitig.reads.is_empty() {
        return i32::MIN as i64;
    }
    let opposite_flag = if flag == u32::from(FATHER) {
        u32::from(MOTHER)
    } else if flag == u32::from(MOTHER) {
        u32::from(FATHER)
    } else {
        u32::MAX
    };

    let start_vertex = unitig.reads[0] >> 32;
    let end_vertex = unitig.reads[unitig.reads.len() - 1] >> 32;
    let s_base = index.ridx[(start_vertex >> 1) as usize] as usize;
    let s_len = (index.ridx[(start_vertex >> 1) as usize + 1] - index.ridx[(start_vertex >> 1) as usize]) as usize;
    let e_base = index.ridx[(end_vertex >> 1) as usize] as usize;
    let e_len = (index.ridx[(end_vertex >> 1) as usize + 1] - index.ridx[(end_vertex >> 1) as usize]) as usize;

    for si in 0..s_len {
        let s_entry = index.ra[s_base + si];
        let s_uid = (s_entry >> 32) as usize;
        let s_unitig = &unitig_graph.unitigs[s_uid];
        if s_unitig.reads.is_empty() || unitig_graph.graph.sequences[s_uid].deleted {
            continue;
        }
        let s_position = s_entry as u32 as usize;
        let reverse = (s_unitig.reads[s_position] >> 32) != start_vertex;

        for ei in 0..e_len {
            let e_entry = index.ra[e_base + ei];
            let e_uid = (e_entry >> 32) as usize;
            if e_uid != s_uid {
                continue;
            }
            let e_unitig = &unitig_graph.unitigs[e_uid];
            if e_unitig.reads.is_empty() || unitig_graph.graph.sequences[e_uid].deleted {
                continue;
            }
            let e_position = e_entry as u32 as usize;

            let (window_length, window_start) = if s_position >= e_position {
                (s_position - e_position + 1, e_position)
            } else {
                (e_position - s_position + 1, s_position)
            };
            if unitig.reads.len() != window_length {
                continue;
            }
            let window = &s_unitig.reads[window_start..window_start + window_length];

            if !reverse {
                let mut k = 0;
                while k < window_length && (unitig.reads[k] >> 32) == (window[k] >> 32) {
                    k += 1;
                }
                if k < window_length {
                    continue;
                }
            } else {
                let mut k = 0;
                while k < window_length
                    && (unitig.reads[k] >> 32) == ((window[window_length - k - 1] >> 32) ^ 1)
                {
                    k += 1;
                }
                if k < window_length {
                    continue;
                }
            }

            if unitig.reads.len() < s_unitig.reads.len() {
                return -1;
            } else if unitig.reads.len() == s_unitig.reads.len() {
                let mut flag_occ = 0u32;
                let mut opposite_occ = 0u32;
                for &packed in &unitig.reads {
                    let label = trio_flags[(packed >> 33) as usize];
                    if u32::from(label) == flag {
                        flag_occ += 1;
                    }
                    if u32::from(label) == opposite_flag {
                        opposite_occ += 1;
                    }
                }
                if flag_occ < opposite_occ {
                    return -1;
                }
                return s_uid as i64;
            }
        }
    }

    i32::MIN as i64
}

/// `dedup_exact_ug`: deletes the exact duplicates between the reference and
/// query phase graphs; returns the removed base total.
#[allow(clippy::too_many_arguments)]
fn dedup_exact_ug(
    ref_index: &DedupIndex,
    ref_graph: &mut UnitigGraph,
    // The reference also receives the query index but only walks `qry->ug`.
    _qry_index: &DedupIndex,
    qry_graph: &mut UnitigGraph,
    read_graph: &Asg,
    coverage_cut: &[CoverageCut],
    sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    primary_flag: &mut [u8],
    trio_flag: u32,
    trio_flags: &[u8],
    options: &AssemblyOptions,
) -> u64 {
    let mut n_base = 0u64;
    for k in 0..qry_graph.unitigs.len() {
        if qry_graph.unitigs[k].reads.is_empty() || qry_graph.graph.sequences[k].deleted {
            continue;
        }
        if if_primary_unitig(
            &qry_graph.unitigs[k],
            read_graph,
            coverage_cut,
            sources,
            read_to_unitig,
            primary_flag,
            options,
        ) {
            continue;
        }
        let found = cal_exact_ug_o(ref_index, ref_graph, &qry_graph.unitigs[k], trio_flag, trio_flags);
        if found == i32::MIN as i64 {
            continue;
        }
        if found == -1 {
            let length = qry_graph.unitigs[k].length;
            delete_ug_node(qry_graph, k);
            n_base += u64::from(length);
        } else {
            let target = found as usize;
            let length = ref_graph.unitigs[target].length;
            delete_ug_node(ref_graph, target);
            n_base += u64::from(length);
        }
    }
    n_base
}

/// `push_ma_utg_t`: appends a copy of `unitig` to the haplotype graph.
fn push_ma_utg(unitig_graph: &mut UnitigGraph, unitig: &MaUtg) {
    let sid = unitig_graph.unitigs.len();
    unitig_graph.graph.set_sequence(sid, unitig.length, false);
    unitig_graph.graph.sequences[sid].coverage_class = 0;

    let mut copy = MaUtg::default();
    copy.length = unitig.length;
    copy.reads = unitig.reads.clone();
    copy.sequence = unitig.sequence.clone();
    if !copy.reads.is_empty() {
        copy.start = (copy.reads[0] >> 32) as u32;
        copy.end = ((copy.reads[copy.reads.len() - 1] >> 32) as u32) ^ 1;
    }
    unitig_graph.unitigs.push(copy);
}

/// The `if(hapNn)` epilogue of `append_miss_nid`: rebuilds the graph index
/// after the appended unitigs.
fn finalize_appended_hap(unitig_graph: &mut UnitigGraph) {
    unitig_graph.graph.is_sorted = false;
    unitig_graph.graph.cleanup();
    unitig_graph.graph.symmetrize();
    if !unitig_graph.graph.sequence_visited.is_empty() {
        let size = unitig_graph.graph.sequences.len() * 2;
        unitig_graph.graph.sequence_visited.clear();
        unitig_graph.graph.sequence_visited.resize(size, 0);
    }
}

/// `append_miss_nid`: re-derives the phase graph over both haplotypes' reads
/// and appends the unitigs neither haplotype took (`Overlaps.cpp` 21514).
#[allow(clippy::too_many_arguments)]
fn append_miss_nid(
    read_graph: &mut Asg,
    hap0: &mut UnitigGraph,
    hap1: &mut UnitigGraph,
    ff: &mut [u8],
    len_cut: u64,
    occ_cut: u64,
    options: &AssemblyOptions,
    trio_flags: &[u8],
) -> u64 {
    ff.fill(0);
    for unitig_graph in [&mut *hap0, &mut *hap1] {
        for v in 0..unitig_graph.unitigs.len() {
            let unitig = &unitig_graph.unitigs[v];
            if unitig.reads.is_empty() || unitig_graph.graph.sequences[v].deleted {
                continue;
            }
            for &packed in &unitig.reads {
                ff[(packed >> 33) as usize] = 1;
            }
        }
    }

    let phase_threshold = phase_cut_threshold(options);
    let mut ug = ma_ug_gen_phase(read_graph, phase_threshold as u32, PHASE_SEP_RATE, trio_flags);

    for k in 0..ug.unitigs.len() {
        if ug.unitigs[k].reads.is_empty() || ug.graph.sequences[k].deleted {
            continue;
        }
        let mut n_set = 0usize;
        for &packed in &ug.unitigs[k].reads {
            if ff[(packed >> 33) as usize] != 0 {
                n_set += 1;
            }
        }
        if n_set > 0 && n_set as f64 >= ug.unitigs[k].reads.len() as f64 * 0.2 {
            delete_ug_node(&mut ug, k);
        }
    }
    let fe: Vec<AsgArc> = Vec::new();
    renew_utg(&mut ug, read_graph, &fe);

    let mut n_base = 0u64;
    let mut hap0_count = 0u64;
    let mut hap1_count = 0u64;
    for k in 0..ug.unitigs.len() {
        let unitig = &ug.unitigs[k];
        if unitig.reads.is_empty() || ug.graph.sequences[k].deleted {
            continue;
        }
        if u64::from(unitig.length) < len_cut || (unitig.reads.len() as u64) < occ_cut {
            continue;
        }
        let mut n_set = 0usize;
        for &packed in &unitig.reads {
            if ff[(packed >> 33) as usize] != 0 {
                n_set += 1;
            }
        }
        if n_set > 0 && n_set as f64 >= unitig.reads.len() as f64 * 0.2 {
            continue;
        }

        let mut flag_occ = 0u32;
        let mut opposite_occ = 0u32;
        for &packed in &unitig.reads {
            let label = trio_flags[(packed >> 33) as usize];
            if u32::from(label) == u32::from(FATHER) {
                flag_occ += 1;
            }
            if u32::from(label) == u32::from(MOTHER) {
                opposite_occ += 1;
            }
        }

        if flag_occ > opposite_occ {
            push_ma_utg(hap0, unitig);
            hap0_count += 1;
        } else {
            push_ma_utg(hap1, unitig);
            hap1_count += 1;
        }
        n_base += u64::from(unitig.length);
    }

    if hap0_count > 0 {
        finalize_appended_hap(hap0);
    }
    if hap1_count > 0 {
        finalize_appended_hap(hap1);
    }
    n_base
}

/// `output_hap_graph`: writes one haplotype's `p_ctg` triple
/// (`Overlaps.cpp` 21237).
#[allow(clippy::too_many_arguments)]
fn output_hap_graph(
    unitig_graph: &mut UnitigGraph,
    read_graph: &mut Asg,
    arcs: &[AsgArc],
    coverage_cut: &[CoverageCut],
    output_file_name: &str,
    flag: u32,
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    max_hang: i32,
    min_ovlp: i32,
    reads: &ReadStore,
    trio_labels: &[u8],
    options: &AssemblyOptions,
) -> Result<(), String> {
    let (file_tag, gfa_tag) = if flag == u32::from(FATHER) {
        ("hap1", "h1tg")
    } else {
        ("hap2", "h2tg")
    };

    let gfa_name = format!("{}.{}.p_ctg.gfa", output_file_name, file_tag);
    eprintln!("Writing {} to disk... ", gfa_name);
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
    let mut output = BufWriter::new(
        File::create(&gfa_name).map_err(|error| format!("failed to create {}: {}", gfa_name, error))?,
    );
    ma_ug_print(
        unitig_graph,
        read_graph,
        Some(coverage_cut),
        Some(sources),
        Some(read_to_unitig),
        gfa_tag,
        &mut output,
        reads,
        trio_labels,
    )
    .map_err(|error| format!("failed to write {}: {}", gfa_name, error))?;
    output
        .flush()
        .map_err(|error| format!("failed to write {}: {}", gfa_name, error))?;

    let noseq_name = format!("{}.{}.p_ctg.noseq.gfa", output_file_name, file_tag);
    let mut output = BufWriter::new(
        File::create(&noseq_name)
            .map_err(|error| format!("failed to create {}: {}", noseq_name, error))?,
    );
    ma_ug_print_simple(
        unitig_graph,
        read_graph,
        Some(coverage_cut),
        Some(sources),
        Some(read_to_unitig),
        gfa_tag,
        &mut output,
        reads,
        trio_labels,
    )
    .map_err(|error| format!("failed to write {}: {}", noseq_name, error))?;
    output
        .flush()
        .map_err(|error| format!("failed to write {}: {}", noseq_name, error))?;

    if options.bed_inconsistency_rate != 0 {
        let bed_name = format!("{}.{}.p_ctg.lowQ.bed", output_file_name, file_tag);
        let mut output = BufWriter::new(
            File::create(&bed_name)
                .map_err(|error| format!("failed to create {}: {}", bed_name, error))?,
        );
        ma_ug_print_bed(
            unitig_graph,
            read_graph,
            coverage_cut,
            sources,
            arcs,
            max_hang,
            min_ovlp,
            options.bed_inconsistency_rate as u32,
            gfa_tag,
            &mut output,
            reads,
            options.maximum_hang_rate,
        )
        .map_err(|error| format!("failed to write {}: {}", bed_name, error))?;
        output
            .flush()
            .map_err(|error| format!("failed to write {}: {}", bed_name, error))?;
    }

    Ok(())
}

/// `output_trio_graph_joint`: the haplotype output driver
/// (`Overlaps.cpp` 23244).
///
/// Both phase graphs are built from the same read graph in sequence (the
/// father pass marks its reads as [`DROP`], which the mother pass sees),
/// deduplicated against each other, topped up with the missed unitigs and
/// written as `bp.hap1.p_ctg.*` / `bp.hap2.p_ctg.*`.
#[allow(clippy::too_many_arguments)]
pub fn output_trio_graph_joint(
    read_graph: &mut Asg,
    coverage_cut: &mut [CoverageCut],
    output_file_name: &str,
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
    gap_fuzz: i32,
    bubble_labels: &mut BubbleLabel,
    options: &mut AssemblyOptions,
    trio_flags: &mut [u8],
    reads: &ReadStore,
) -> Result<(), String> {
    reduce_hamming_error_adv(
        None,
        read_graph,
        sources,
        coverage_cut,
        max_hang,
        min_ovlp,
        gap_fuzz,
        trio_flags,
        options.maximum_hang_rate,
    );

    let mut arcs0: Vec<AsgArc> = Vec::new();
    let mut arcs1: Vec<AsgArc> = Vec::new();
    let mut hu0 = output_trio_unitig_graph(
        read_graph,
        coverage_cut,
        output_file_name,
        u32::from(FATHER),
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
        gap_fuzz,
        true,
        bubble_labels,
        &mut arcs0,
        options,
        trio_flags,
        reads,
    )?;
    let mut hu1 = output_trio_unitig_graph(
        read_graph,
        coverage_cut,
        output_file_name,
        u32::from(MOTHER),
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
        gap_fuzz,
        true,
        bubble_labels,
        &mut arcs1,
        options,
        trio_flags,
        reads,
    )?;

    let mut ff = vec![0u8; read_graph.sequences.len()];
    let hidx0 = DedupIndex::new(&hu0, read_graph.sequences.len());
    let hidx1 = DedupIndex::new(&hu1, read_graph.sequences.len());
    update_recover_atg_cov(options);

    let mut primary_flag = vec![0u8; read_graph.sequences.len()];
    let mut dedup_base = 0u64;
    dedup_base += dedup_exact_ug(
        &hidx1,
        &mut hu1,
        &hidx0,
        &mut hu0,
        read_graph,
        coverage_cut,
        sources,
        read_to_unitig,
        &mut primary_flag,
        u32::from(FATHER),
        trio_flags,
        options,
    );
    dedup_base += dedup_exact_ug(
        &hidx0,
        &mut hu0,
        &hidx1,
        &mut hu1,
        read_graph,
        coverage_cut,
        sources,
        read_to_unitig,
        &mut primary_flag,
        u32::from(MOTHER),
        trio_flags,
        options,
    );

    let miss_base = append_miss_nid(
        read_graph,
        &mut hu0,
        &mut hu1,
        &mut ff,
        PHASE_MISS_LEN,
        PHASE_MISS_N,
        options,
        trio_flags,
    );

    renew_utg(&mut hu0, read_graph, &arcs0);
    renew_utg(&mut hu1, read_graph, &arcs1);
    eprintln!(
        "[M::output_trio_graph_joint] dedup_base::{}, miss_base::{}",
        dedup_base, miss_base
    );

    if options.self_scaffolding != 0 {
        return Err(
            "output_trio_graph_joint: the --dual-scaf scaffolding path is not ported yet"
                .to_string(),
        );
    }

    output_hap_graph(
        &mut hu0,
        read_graph,
        &arcs0,
        coverage_cut,
        output_file_name,
        u32::from(FATHER),
        sources,
        read_to_unitig,
        max_hang,
        min_ovlp,
        reads,
        trio_flags,
        options,
    )?;
    output_hap_graph(
        &mut hu1,
        read_graph,
        &arcs1,
        coverage_cut,
        output_file_name,
        u32::from(MOTHER),
        sources,
        read_to_unitig,
        max_hang,
        min_ovlp,
        reads,
        trio_flags,
        options,
    )?;

    Ok(())
}


// ---------------------------------------------------------------------------
// `output_contig_graph_primary` / `output_contig_graph_alternative`
// (`Overlaps.cpp` 32644 / 32699): the `--primary` tail of `clean_graph`.
// ---------------------------------------------------------------------------

/// The `ma_ug_seq` + GFA triple shared by the two `--primary` drivers.
///
/// The C++ duplicates this block per driver with a fixed log line and tag;
/// the log text is part of the consistency scope, so it stays parameterised.
#[allow(clippy::too_many_arguments)]
fn print_primary_contig_graph(
    unitig_graph: &mut UnitigGraph,
    read_graph: &mut Asg,
    coverage_cut: &[CoverageCut],
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    new_rtg_edges: &[AsgArc],
    output_file_name: &str,
    file_tag: &str,
    gfa_tag: &str,
    log_line: &str,
    reads: &ReadStore,
    trio_labels: &[u8],
    options: &AssemblyOptions,
) -> Result<(), String> {
    ma_ug_seq(
        unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        new_rtg_edges,
        None,
        options.maximum_hang_length,
        options.minimum_overlap_length,
        true,
        reads,
        options.maximum_hang_rate,
    );

    eprintln!("{}", log_line);
    let gfa_name = format!("{}.{}.gfa", output_file_name, file_tag);
    let mut output = BufWriter::new(
        File::create(&gfa_name).map_err(|error| format!("failed to create {}: {}", gfa_name, error))?,
    );
    ma_ug_print(
        unitig_graph,
        read_graph,
        Some(coverage_cut),
        Some(sources),
        Some(read_to_unitig),
        gfa_tag,
        &mut output,
        reads,
        trio_labels,
    )
    .map_err(|error| format!("failed to write {}: {}", gfa_name, error))?;
    output
        .flush()
        .map_err(|error| format!("failed to write {}: {}", gfa_name, error))?;

    let noseq_name = format!("{}.{}.noseq.gfa", output_file_name, file_tag);
    let mut output = BufWriter::new(
        File::create(&noseq_name)
            .map_err(|error| format!("failed to create {}: {}", noseq_name, error))?,
    );
    ma_ug_print_simple(
        unitig_graph,
        read_graph,
        Some(coverage_cut),
        Some(sources),
        Some(read_to_unitig),
        gfa_tag,
        &mut output,
        reads,
        trio_labels,
    )
    .map_err(|error| format!("failed to write {}: {}", noseq_name, error))?;
    output
        .flush()
        .map_err(|error| format!("failed to write {}: {}", noseq_name, error))?;

    if options.bed_inconsistency_rate != 0 {
        let bed_name = format!("{}.{}.lowQ.bed", output_file_name, file_tag);
        let mut output = BufWriter::new(
            File::create(&bed_name)
                .map_err(|error| format!("failed to create {}: {}", bed_name, error))?,
        );
        ma_ug_print_bed(
            unitig_graph,
            read_graph,
            coverage_cut,
            sources,
            new_rtg_edges,
            options.maximum_hang_length,
            options.minimum_overlap_length,
            options.bed_inconsistency_rate as u32,
            gfa_tag,
            &mut output,
            reads,
            options.maximum_hang_rate,
        )
        .map_err(|error| format!("failed to write {}: {}", bed_name, error))?;
        output
            .flush()
            .map_err(|error| format!("failed to write {}: {}", bed_name, error))?;
    }

    Ok(())
}

/// `output_contig_graph_primary`: the primary contig output of the
/// non-partition (`--primary`) dispatch.
///
/// Runs [`adjust_utg_by_primary`] on the live read graph — the reference
/// passes `sg` itself here, so the read classes are relabelled and the
/// unclaimed reads dropped before [`output_contig_graph_alternative`] runs —
/// and writes the `p_ctg` triple. `break_ug_contig` for
/// `b_low_cov`/`b_high_cov` is not ported; both switches are off in the
/// default configuration.
#[allow(clippy::too_many_arguments)]
pub fn output_contig_graph_primary(
    read_graph: &mut Asg,
    coverage_cut: &mut [CoverageCut],
    output_file_name: &str,
    sources: &[PafList],
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    bubble_labels: &mut BubbleLabel,
    options: &mut AssemblyOptions,
    reads: &ReadStore,
    trio_flags: &mut [u8],
) -> Result<(), String> {
    let mut unitig_graph = ma_ug_gen_primary(read_graph, PRIMARY_LABEL);
    let new_rtg_edges: Vec<AsgArc> = Vec::new();

    adjust_utg_by_primary(
        &mut unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        reverse_sources,
        read_to_unitig,
        primary_tips_length(options),
        0.15,
        3,
        0.05,
        0.9,
        options.maximum_hang_length,
        options.minimum_overlap_length,
        &new_rtg_edges,
        Some(bubble_labels),
        options,
        trio_flags,
        false,
        false,
    )?;

    if options.breakpoint_low_coverage > 0 || options.breakpoint_high_coverage > 0 {
        return Err(
            "output_contig_graph_primary needs break_ug_contig (b_low_cov/b_high_cov), which is not ported yet"
                .to_string(),
        );
    }

    print_primary_contig_graph(
        &mut unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        read_to_unitig,
        &new_rtg_edges,
        output_file_name,
        "p_ctg",
        "ptg",
        "Writing primary contig GFA to disk... ",
        reads,
        trio_flags,
        options,
    )
}

/// `output_contig_graph_alternative`: the alternate contig output of the
/// non-partition dispatch. The read graph has already been relabelled by
/// [`output_contig_graph_primary`], so the alternate unitig graph is built
/// from the dropped reads.
#[allow(clippy::too_many_arguments)]
pub fn output_contig_graph_alternative(
    read_graph: &mut Asg,
    coverage_cut: &[CoverageCut],
    output_file_name: &str,
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    reads: &ReadStore,
    trio_labels: &[u8],
    options: &AssemblyOptions,
) -> Result<(), String> {
    let mut unitig_graph = ma_ug_gen_primary(read_graph, ALTER_LABEL);
    let new_rtg_edges: Vec<AsgArc> = Vec::new();

    print_primary_contig_graph(
        &mut unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        read_to_unitig,
        &new_rtg_edges,
        output_file_name,
        "a_ctg",
        "atg",
        "Writing alternate contig GFA to disk... ",
        reads,
        trio_labels,
        options,
    )
}

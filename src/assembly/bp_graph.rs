//! `output_bp_graph` and its p_ctg-side callees `adjust_utg_by_primary` and
//! `print_utg` from `Overlaps.cpp`.
//!
//! The default HiFi pipeline reaches this driver through the
//! `HA_F_PARTITION && purge_level_primary > 0` branch of `clean_graph`. It
//! clones the read graph and the primary unitig graph, cleans the clone with
//! the coverage-based purge pass, and writes `bp.p_ctg.*`; the trio output
//! stage that follows writes the haplotype graphs.
//!
//! The haplotype half of `output_bp_graph` (`clean_u_trans_t_idx`,
//! `set_trio_flag_by_cov` and `output_trio_graph_joint`) is not ported yet, so
//! [`output_bp_graph`] returns an error after the primary contig files are
//! written.

use std::fs::File;
use std::io::{BufWriter, Write};

use crate::assembly::asg::{Asg, AsgArc};
use crate::assembly::bubble::BubbleLabel;
use crate::assembly::hap_cov::{init_hap_cov_t, HapCov};
use crate::assembly::hit::{CoverageCut, ReadToUnitig, ALTER_LABEL, PRIMARY_LABEL};
use crate::assembly::purge_dups::{purge_dups, HOM_PEAK_RATE};
use crate::assembly::rescue::{
    rescue_contained_reads_aggressive, rescue_missing_overlaps_aggressive, RescueParameters,
};
use crate::assembly::trio_clean::{asg_seq_drop, clean_primary_untig_graph};
use crate::assembly::trans_clean::{clean_u_trans_t_idx, set_trio_flag_by_cov};
use crate::assembly::unitig::{
    adjust_utg_advance, copy_read_graph, copy_untig_graph, delete_useless_nodes,
    discard_small_ctg, hic_clean, ma_ug_gen_primary, ma_ug_print, ma_ug_print_bed,
    ma_ug_print_simple, ma_ug_seq, recover_utg_by_coverage, renew_utg, UnitigGraph,
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

    // `output_trio_graph_joint` starts by closing the hamming-error bubbles of
    // the read graph; the phase graphs below are built from the result.
    crate::assembly::hamming::reduce_hamming_error_adv(
        None,
        read_graph,
        sources,
        coverage_cut,
        options.maximum_hang_length,
        options.minimum_overlap_length,
        options.gap_fuzz,
        trio_flags,
        options.maximum_hang_rate,
    );

    Err(
        "output_bp_graph: the haplotype output stage (output_trio_graph_joint) is not \
         implemented yet"
            .to_string(),
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

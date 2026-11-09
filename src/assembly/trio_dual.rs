//! The trio-dual output driver (`--trio-dual`), the counterpart of
//! `output_bp_trio_graph` (`Overlaps.cpp` 17680), `trio_phasing_refine`
//! (`hic.cpp` 17741) and `gen_bubble_chain` (`gfa_ut.cpp` 3306).
//!
//! The stage re-phases the primary unitig graph with the max-cut partitioner
//! (`rcut.cpp`) and rewrites `R_INF.trio_flag` before the regular trio output
//! writes the two haplotype contig graphs. `r_utg`/`p_utg` are written earlier
//! and are identical between the two trio modes.

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::bp_graph::{adjust_utg_by_primary, output_trio_graph_joint};
use crate::assembly::bubble::BubbleLabel;
use crate::assembly::bubble_chain::{identify_bubbles, BubbleType};
use crate::assembly::hap_cov::{UTrans, RC_2};
use crate::assembly::hit::{CoverageCut, ReadToUnitig, PRIMARY_LABEL};
use crate::assembly::maxcut::mc_solve;
use crate::assembly::trans_clean::{
    clean_u_trans_t_idx_filter_adv, kt_u_trans_t_idx,
};
use crate::assembly::unitig::{
    copy_read_graph, copy_untig_graph, hic_clean, ma_ug_gen_primary, UnitigGraph, AMBIGU, FATHER,
    MOTHER,
};
use crate::command_line::options::AssemblyOptions;
use crate::correction::ec_driver::PafList;
use crate::reads::read_store::ReadStore;

/// `ug_occ_w`: the number of reads of `u` contained in `[is, ie)`.
fn ug_occ_w(is: u64, ie: u64, unitig: &MaUtg, read_lengths: &[u32]) -> u64 {
    if is == 0 && ie == u64::from(unitig.length) {
        return unitig.reads.len() as u64;
    }
    let mut offset = 0u64;
    let mut occ = 0u64;
    for &packed in &unitig.reads {
        let read_start = offset;
        let read_end = offset + u64::from(read_lengths[(packed >> 33) as usize]);
        if is <= read_start && ie >= read_end {
            occ += 1;
        }
        if read_start >= ie {
            break;
        }
        offset += u64::from(packed as u32);
    }
    occ
}

/// `ug_occ_hap_w`: the number of trio-marked reads of `u` contained in
/// `[is, ie)`.
fn ug_occ_hap_w(
    is: u64,
    ie: u64,
    unitig: &MaUtg,
    read_lengths: &[u32],
    trio_flags: &[u8],
) -> u64 {
    let mut offset = 0u64;
    let mut occ = 0u64;
    for &packed in &unitig.reads {
        let read_start = offset;
        let read_end = offset + u64::from(read_lengths[(packed >> 33) as usize]);
        if is <= read_start && ie >= read_end {
            let flag = trio_flags[(packed >> 33) as usize];
            if flag == FATHER || flag == MOTHER {
                occ += 1;
            }
        }
        if read_start >= ie {
            break;
        }
        offset += u64::from(packed as u32);
    }
    occ
}

/// `gen_bubble_chain`: the bubble chains of `unitig_graph`, plus the
/// heterozygosity flags computed by a throw-away `adjust_utg_by_primary`.
///
/// The reference runs `adjust_utg_by_primary` on copies of the read and
/// unitig graphs only to obtain `cov->t_ch->ir_het`; the bubble detection
/// itself runs on the original unitig graph.
#[allow(clippy::too_many_arguments)]
/// `gen_bubble_chain` (`gfa_ut.cpp` 3306): the shallow bubble classification
/// used by the trio phasing and by `--dual-scaf`.
pub(crate) fn gen_bubble_chain(
    read_graph: &Asg,
    unitig_graph: &mut UnitigGraph,
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
    bubble_labels: &mut BubbleLabel,
    options: &mut AssemblyOptions,
    trio_flags: &mut [u8],
) -> Result<BubbleType, String> {
    let mut copy_read = copy_read_graph(read_graph);
    let mut copy_unitig = copy_untig_graph(unitig_graph);
    let mut new_edges: Vec<AsgArc> = Vec::new();
    let cov = adjust_utg_by_primary(
        &mut copy_unitig,
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
    let ir_het = cov.is_r_het.clone();

    let mut bub = BubbleType::default();
    identify_bubbles(
        &mut unitig_graph.graph,
        &unitig_graph.unitigs,
        &mut bub,
        &ir_het,
        None,
    );

    if options.purge_level_primary == 0 {
        for k in 0..unitig_graph.graph.sequences.len() {
            if bub.is_hom(k as u32) {
                bub.index[k] = bub.full_bubble_count as u32 + 1;
            }
        }
    }
    Ok(bub)
}

/// `trio_phasing_refine`: phases the unitig graph with the transition table
/// collected by `adjust_utg_by_primary`, then overwrites the trio labels of
/// every read of the assigned unitigs.
#[allow(clippy::too_many_arguments)]
fn trio_phasing_refine(
    input_unitig_graph: &UnitigGraph,
    read_graph: &mut Asg,
    transitions: &mut Vec<UTrans>,
    transition_index: &mut Vec<u64>,
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
    bubble_labels: &mut BubbleLabel,
    options: &mut AssemblyOptions,
    trio_flags: &mut [u8],
) -> Result<(), String> {
    let mut unitig_graph = copy_untig_graph(input_unitig_graph);
    let ug_n0 = unitig_graph.graph.sequences.len() as u32;
    let read_lengths: Vec<u32> = read_graph
        .sequences
        .iter()
        .map(|sequence| sequence.length)
        .collect();

    let mut bub = gen_bubble_chain(
        read_graph,
        &mut unitig_graph,
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

    let logical_n = clean_u_trans_t_idx_filter_adv(
        transitions,
        transition_index,
        &unitig_graph,
        read_graph,
        0.95,
        false,
    );
    transitions.truncate(logical_n);

    // The two pseudo unitigs that carry the father and mother anchors.
    let frid = ug_n0;
    let mrid = ug_n0 + 1;
    for (id, length) in [(frid, u32::from(FATHER)), (mrid, u32::from(MOTHER))] {
        unitig_graph
            .graph
            .set_sequence(id as usize, length, false);
        unitig_graph.unitigs.push(MaUtg {
            length,
            circular: true,
            start: u32::MAX,
            end: u32::MAX,
            reads: vec![(u64::from(u32::MAX) << 32) | u64::from(length)],
            sequence: Vec::new(),
        });
    }
    unitig_graph.graph.index.clear();
    unitig_graph.graph.is_sorted = false;
    unitig_graph.graph.cleanup();

    let new_count = unitig_graph.graph.sequences.len();
    bub.index
        .resize(new_count, bub.full_bubble_count as u32 + 1);
    bub.b_s_idx.resize(new_count, u64::MAX);

    // Bubbles / trio coverage of the original unitigs.
    let mut cis_n = 0u64;
    let mut tot_hap = 0u64;
    let mut tot_r = 0u64;
    for k in 0..ug_n0 as usize {
        if unitig_graph.graph.sequences[k].deleted {
            continue;
        }
        let unitig = &unitig_graph.unitigs[k];
        tot_r += unitig.reads.len() as u64;
        if unitig.reads.is_empty() {
            continue;
        }
        let mut fh = false;
        let mut mh = false;
        for &packed in &unitig.reads {
            let flag = trio_flags[(packed >> 33) as usize];
            if flag == FATHER {
                tot_hap += 1;
                fh = true;
            } else if flag == MOTHER {
                tot_hap += 1;
                mh = true;
            }
            if flag == FATHER || flag == MOTHER {
                tot_hap += 1;
            }
        }
        if fh {
            cis_n += 2;
        }
        if mh {
            cis_n += 2;
        }
    }

    // The anchor transitions, then the re-index.
    for (qn, tn) in [(frid, mrid), (mrid, frid)] {
        transitions.push(UTrans {
            query_start: 0,
            query_end: unitig_graph.graph.sequences[qn as usize].length,
            query_identifier: qn,
            target_start: 0,
            target_end: unitig_graph.graph.sequences[tn as usize].length,
            target_identifier: tn,
            occurrence: 0,
            network_weight: f64::MAX / 2.0,
            flag: RC_2,
            reverse: false,
            deleted: false,
        });
    }
    kt_u_trans_t_idx(transitions, transition_index, unitig_graph.graph.sequences.len() as u32);

    // `trans_n + cis_n` links between the original unitigs.
    let mut links: Vec<UTrans> = Vec::new();
    let mut trans_n = 0usize;
    for entry in transitions.iter() {
        if entry.deleted {
            continue;
        }
        let hom_q = bub.is_hom(entry.query_identifier);
        let hom_t = bub.is_hom(entry.target_identifier);
        if hom_q || hom_t {
            continue;
        }
        let mut weight = 0u64;
        if entry.query_identifier < ug_n0 && entry.target_identifier < ug_n0 {
            let q = &unitig_graph.unitigs[entry.query_identifier as usize];
            let t = &unitig_graph.unitigs[entry.target_identifier as usize];
            weight += ug_occ_hap_w(
                u64::from(entry.query_start),
                u64::from(entry.query_end),
                q,
                &read_lengths,
                trio_flags,
            );
            weight += ug_occ_hap_w(
                u64::from(entry.target_start),
                u64::from(entry.target_end),
                t,
                &read_lengths,
                trio_flags,
            );
            weight += ug_occ_w(
                u64::from(entry.query_start),
                u64::from(entry.query_end),
                q,
                &read_lengths,
            );
            weight += ug_occ_w(
                u64::from(entry.target_start),
                u64::from(entry.target_end),
                t,
                &read_lengths,
            );
        } else {
            weight += tot_hap + tot_r;
        }
        weight >>= 1;
        let mut link = *entry;
        link.network_weight = if weight != 0 { weight as f64 } else { 1.0 };
        links.push(link);
        trans_n += 1;
    }

    // The anchor links of every unitig that carries trio evidence.
    for k in 0..ug_n0 as usize {
        let unitig = &unitig_graph.unitigs[k];
        if unitig.reads.is_empty() {
            continue;
        }
        let mut fh = 0i64;
        let mut mh = 0i64;
        for &packed in &unitig.reads {
            let flag = trio_flags[(packed >> 33) as usize];
            if flag == FATHER {
                fh -= 1;
            }
            if flag == MOTHER {
                mh -= 1;
            }
        }
        for (anchor, count) in [(frid, fh), (mrid, mh)] {
            if count == 0 {
                continue;
            }
            let bound = u64::from(
                unitig_graph.graph.sequences[anchor as usize]
                    .length
                    .min(unitig_graph.graph.sequences[k].length),
            );
            let forward = UTrans {
                query_start: 0,
                query_end: bound as u32,
                query_identifier: anchor,
                target_start: 0,
                target_end: bound as u32,
                target_identifier: k as u32,
                occurrence: 0,
                network_weight: count as f64,
                flag: RC_2,
                reverse: false,
                deleted: false,
            };
            let mut reverse = forward;
            reverse.query_identifier = k as u32;
            reverse.target_identifier = anchor;
            links.push(forward);
            links.push(reverse);
        }
    }

    debug_assert_eq!(links.len() as u64, trans_n as u64 + cis_n);

    let mut link_index = Vec::new();
    kt_u_trans_t_idx(&mut links, &mut link_index, unitig_graph.graph.sequences.len() as u32);

    let mut s = vec![0i8; unitig_graph.graph.sequences.len()];
    s[frid as usize] = 1;
    s[mrid as usize] = -1;
    mc_solve(
        &unitig_graph,
        &links,
        transitions,
        transition_index,
        &mut s,
        true,
        &bub,
        options.thread_count.max(1) as u64,
        options.perturbation_round_count,
        options.perturbation_fraction,
        options.random_seed,
    );

    if s[frid as usize] != 0
        && s[mrid as usize] != 0
        && s[frid as usize] != s[mrid as usize]
    {
        trio_flags.fill(AMBIGU);
        for i in 0..ug_n0 as usize {
            if unitig_graph.graph.sequences[i].deleted {
                continue;
            }
            if s[i] == 0 {
                continue;
            }
            let flag = if s[i] == s[frid as usize] {
                FATHER
            } else {
                MOTHER
            };
            for &packed in &unitig_graph.unitigs[i].reads {
                trio_flags[(packed >> 33) as usize] = flag;
            }
        }
    }
    Ok(())
}

/// `output_bp_trio_graph`: the `--trio-dual` output driver.
#[allow(clippy::too_many_arguments)]
pub fn output_bp_trio_graph(
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
    hic_clean(read_graph);

    let mut unitig_graph = ma_ug_gen_primary(read_graph, PRIMARY_LABEL);
    let mut copy_read = copy_read_graph(read_graph);
    let mut copy_unitig = copy_untig_graph(&unitig_graph);
    let mut new_edges: Vec<AsgArc> = Vec::new();
    let cov = adjust_utg_by_primary(
        &mut copy_unitig,
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

    let mut transitions = cov
        .t_ch
        .as_ref()
        .map(|chain| chain.k_trans.clone())
        .unwrap_or_default();
    let mut transition_index = cov
        .t_ch
        .as_ref()
        .map(|chain| chain.k_trans_idx.clone())
        .unwrap_or_default();

    trio_phasing_refine(
        &unitig_graph,
        read_graph,
        &mut transitions,
        &mut transition_index,
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

    output_trio_graph_joint(
        read_graph,
        coverage_cut,
        output_file_name,
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
        bubble_labels,
        options,
        trio_flags,
        reads,
    )
}

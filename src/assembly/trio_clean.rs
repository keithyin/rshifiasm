//! The trio-aware tip cleaning of `clean_primary_untig_graph`, the
//! counterparts of `check_tip`, `cal_trio_vec`, `check_different_haps`,
//! `cut_trio_tip_primary` and `asg_seq_drop` in `Overlaps.cpp`/`Overlaps.h`.
//!
//! [`check_different_haps`] decides whether two chains belong to different
//! haplotypes: it labels every read of the longer chain in the read-to-unitig
//! map and counts how many reads of the shorter chain, or of their reverse
//! overlaps, run into a labelled unitig. The map and the `is_r_het` flags are
//! shared scratch state, exactly as in the reference.

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::bubble_identify::{asg_pop_bubble_primary_trio, get_unitig};
use crate::assembly::gfa_clean::{
    get_real_length, get_real_length_one, END_TIPS, LONG_TIPS, LOOP, MUL_INPUT, MUL_OUTPUT,
};
use crate::assembly::hap_cov::{HapCov, C_HET, P_HET};
use crate::assembly::hit::{CoverageCut, ReadToUnitig, ALTER_LABEL, HAP_LABEL};
use crate::assembly::tangle::resolve_tangles;
use crate::assembly::trans_cov::collect_trans_cov;
use crate::assembly::unitig::{UnitigGraph, AMBIGU, FATHER, MOTHER};
use crate::container::ksort::radix_sort_u64;
use crate::correction::ec_driver::PafList;

/// `UNAVAILABLE`: the two chains cannot be compared.
pub const UNAVAILABLE: u32 = u32::MAX;
/// `PLOID`: the two chains agree, so they are not from different haplotypes.
pub const PLOID: u32 = 0;
/// `NON_PLOID`: the two chains disagree.
pub const NON_PLOID: u32 = 1;

/// `TRIM`: the tip reads are dropped and their edges removed.
pub const TRIM: u8 = 10;
/// `CUT`: the tip is cut because it comes from another haplotype.
pub const CUT: u8 = 11;
/// `CUT_DIF_HAP`: the tip is cut because its haplotypes differ.
pub const CUT_DIF_HAP: u8 = 12;

/// `TRIO_THRES`: the majority a trio label needs inside a tip.
pub const TRIO_THRES: f32 = 0.9;
/// `HET_HOM_RATE`: the share of heterozygous reads a ploid decision needs.
pub const HET_HOM_RATE: f64 = 0.7;
/// `DROP`: the label of a read the trio data discards.
pub const DROP: u32 = 5;
/// `DOUBLE_CHECK_THRES`: the second trio-label threshold of [`trio_check`].
pub const DOUBLE_CHECK_THRES: f64 = 0.1;
/// `TRIO_DROP_THRES`: the share of foreign reads a tip may carry.
pub const TRIO_DROP_THRES: f64 = 0.9;
/// `TRIO_DROP_LENGTH_THRES`: how long an alternative tip must be to compete.
pub const TRIO_DROP_LENGTH_THRES: f64 = 0.8;
/// `MAX_STOP_RATE`: how much of a chain a stop must cover to be accepted.
pub const MAX_STOP_RATE: f64 = 0.6;
/// `CHIMERIC_TRIM_THRES`: the read count below which a chimeric chain is cut.
pub const CHIMERIC_TRIM_THRES: u32 = 4;

/// `Trio_counter`: the trio labels of one tip chain.
#[derive(Debug, Clone, Copy, Default)]
pub struct TrioCounter {
    pub father_occ: u32,
    pub mother_occ: u32,
    pub ambig_occ: u32,
    pub drop_occ: u32,
    pub total: u32,
}

/// `check_tip`: walks the unbranched chain from `beg_node`, collecting the
/// visited vertices in `buffer`.
///
/// The chain must branch inwards at every step to continue. The returned flag
/// is one of [`END_TIPS`], [`MUL_OUTPUT`], [`MUL_INPUT`], [`LOOP`] or
/// [`LONG_TIPS`]; `end_node` is the last visited vertex.
pub fn check_tip(
    graph: &Asg,
    beg_node: u32,
    end_node: &mut u32,
    buffer: &mut Vec<u32>,
    max_ext: u32,
) -> u8 {
    let mut v = beg_node;
    let mut chain_len = 0u32;
    *end_node = u32::MAX;
    buffer.clear();
    loop {
        let live_arcs = get_real_length(graph, v);
        *end_node = v;
        chain_len = chain_len.wrapping_add(1);
        buffer.push(v);
        if live_arcs == 0 {
            return END_TIPS;
        }
        if live_arcs > 1 {
            return MUL_OUTPUT;
        }
        let w = get_real_length_one(graph, v).expect("a single live arc has one target");
        // `kv` is one here; the chain continues only while the next vertex
        // has exactly one incoming live arc.
        if get_real_length(graph, w ^ 1) != 1 {
            return MUL_INPUT;
        }
        v = w;
        if v == beg_node {
            return LOOP;
        }
        if chain_len >= max_ext {
            return LONG_TIPS;
        }
    }
}

/// `cal_trio_vec`: the majority trio label of the reads on `nodes`.
///
/// With `unitigs == None` every node is a read; otherwise every read stored
/// in each node's unitig is counted. The reference skips a unitig whose `m`
/// is zero, which only happens when it stores no reads at all.
pub fn cal_trio_vec(
    nodes: &[u32],
    unitigs: Option<&[MaUtg]>,
    trio_labels: &[u8],
    threshold: f32,
) -> u8 {
    let mut father_occ = 0u32;
    let mut mother_occ = 0u32;
    if let Some(unitigs) = unitigs {
        for &node in nodes {
            let unitig = &unitigs[(node >> 1) as usize];
            for &packed in &unitig.reads {
                let read = (packed >> 33) as usize;
                if trio_labels[read] == FATHER {
                    father_occ += 1;
                } else if trio_labels[read] == MOTHER {
                    mother_occ += 1;
                }
            }
        }
    } else {
        for &node in nodes {
            let read = (node >> 1) as usize;
            if trio_labels[read] == FATHER {
                father_occ += 1;
            } else if trio_labels[read] == MOTHER {
                mother_occ += 1;
            }
        }
    }

    let total = father_occ.wrapping_add(mother_occ) as f32;
    if father_occ as f32 >= threshold * total {
        return FATHER;
    }
    if mother_occ as f32 >= threshold * total {
        return MOTHER;
    }
    AMBIGU
}

/// `check_different_haps`: whether two chains carry different haplotypes.
///
/// Both chains are collected first; the longer one's reads are labelled in
/// `read_to_unitig`, and every read of the shorter one is counted as
/// heterozygous (`min_count`) or as reaching a labelled unitig through itself
/// or one of its reverse overlaps (`max_count`). The labels are removed again
/// before returning.
///
/// Returns [`UNAVAILABLE`] when a chain loops or is too short, [`PLOID`] when
/// the overlap counts and the heterozygosity share pass the thresholds, and
/// [`NON_PLOID`] otherwise.
#[allow(clippy::too_many_arguments)]
pub fn check_different_haps(
    graph: &Asg,
    unitigs: Option<&[MaUtg]>,
    read_graph: &Asg,
    v_0: u32,
    v_1: u32,
    reverse_sources: &[PafList],
    b_0: &mut Vec<u32>,
    b_1: &mut Vec<u32>,
    read_to_unitig: &mut ReadToUnitig,
    is_r_het: &[u8],
    min_edge_length: u32,
    stops_threshold: u32,
    purge_similarity_threshold: f32,
) -> u32 {
    b_0.clear();
    b_1.clear();
    let mut end_node = 0u32;
    let mut node_len_0 = 0i64;
    let mut node_len_1 = 0i64;
    let mut base_len = 0i64;
    let mut max_stop_node_len = 0i64;
    let mut max_stop_base_len = 0i64;

    if get_unitig(
        graph,
        unitigs,
        v_0,
        &mut end_node,
        &mut node_len_0,
        &mut base_len,
        &mut max_stop_node_len,
        &mut max_stop_base_len,
        stops_threshold,
        Some(b_0),
    ) == LOOP
    {
        return UNAVAILABLE;
    }
    if get_unitig(
        graph,
        unitigs,
        v_1,
        &mut end_node,
        &mut node_len_1,
        &mut base_len,
        &mut max_stop_node_len,
        &mut max_stop_base_len,
        stops_threshold,
        Some(b_1),
    ) == LOOP
    {
        return UNAVAILABLE;
    }
    if node_len_0 <= i64::from(min_edge_length) || node_len_1 <= i64::from(min_edge_length) {
        return UNAVAILABLE;
    }

    let (b_min, b_max): (&[u32], &[u32]) = if node_len_0 <= node_len_1 {
        (b_0.as_slice(), b_1.as_slice())
    } else {
        (b_1.as_slice(), b_0.as_slice())
    };

    let mut max_count = 0u32;
    let mut min_count = 0u32;
    let mut n_het = 0u32;
    let mut n_hom = 0u32;

    if let Some(unitigs) = unitigs {
        // Label every read of the longer chain with its unitig.
        for &node in b_max {
            let unitig = &unitigs[(node >> 1) as usize];
            for &packed in &unitig.reads {
                let qn = (packed >> 33) as u32;
                read_to_unitig.set(
                    qn,
                    node >> 1,
                    true,
                    Some(read_graph.sequences[qn as usize].coverage_class),
                );
            }
        }

        for &node in b_min {
            let unitig = &unitigs[(node >> 1) as usize];
            for &packed in &unitig.reads {
                let qn = (packed >> 33) as u32;
                if !reverse_sources[qn as usize].hits.is_empty() {
                    min_count += 1;
                }
                if is_r_het[qn as usize] & (C_HET | P_HET) != 0 {
                    n_het += 1;
                }
                n_hom += 1;

                for hit in &reverse_sources[qn as usize].hits {
                    let mut tn = hit.target_identifier;
                    if read_graph.sequences[tn as usize].deleted {
                        let (mapped, is_unitig) = read_to_unitig.get(tn);
                        tn = mapped;
                        if tn == u32::MAX
                            || is_unitig == 1
                            || read_graph.sequences[tn as usize].deleted
                        {
                            continue;
                        }
                    }

                    let (unitig_id, is_unitig) = read_to_unitig.get(tn);
                    if unitig_id != u32::MAX && is_unitig == 1 {
                        max_count += 1;
                        break;
                    }
                }
            }
        }

        for &node in b_max {
            let unitig = &unitigs[(node >> 1) as usize];
            for &packed in &unitig.reads {
                read_to_unitig.set_raw((packed >> 33) as u32, u32::MAX);
            }
        }
    } else {
        // Without a unitig graph every node is a read of the string graph.
        for &node in b_max {
            let qn = node >> 1;
            read_to_unitig.set(
                qn,
                1,
                true,
                Some(read_graph.sequences[qn as usize].coverage_class),
            );
        }

        for &node in b_min {
            let qn = node >> 1;
            if !reverse_sources[qn as usize].hits.is_empty() {
                min_count += 1;
            }
            if is_r_het[qn as usize] & (C_HET | P_HET) != 0 {
                n_het += 1;
            }
            n_hom += 1;

            for hit in &reverse_sources[qn as usize].hits {
                let mut tn = hit.target_identifier;
                if graph.sequences[tn as usize].deleted {
                    let (mapped, is_unitig) = read_to_unitig.get(tn);
                    tn = mapped;
                    if tn == u32::MAX || is_unitig == 1 || graph.sequences[tn as usize].deleted {
                        continue;
                    }
                }

                let (unitig_id, is_unitig) = read_to_unitig.get(tn);
                if unitig_id != u32::MAX && is_unitig == 1 {
                    max_count += 1;
                    break;
                }
            }
        }

        for &node in b_max {
            read_to_unitig.set_raw(node >> 1, u32::MAX);
        }
    }

    if min_count == 0 {
        return UNAVAILABLE;
    }
    if max_count as f32 > min_count as f32 * purge_similarity_threshold
        && n_het as f64 >= n_hom as f64 * HET_HOM_RATE
    {
        return PLOID;
    }
    NON_PLOID
}

/// `asg_seq_drop`: removes the edges of `s` that leave the alternate label.
///
/// Only sequences already marked [`ALTER_LABEL`] are processed. For both
/// orientations every live outgoing arc whose target is not [`ALTER_LABEL`]
/// is deleted, together with its mirror arc. The sequence itself is not
/// marked deleted.
pub fn asg_seq_drop(graph: &mut Asg, s: u32) {
    if graph.sequences[s as usize].coverage_class != ALTER_LABEL {
        return;
    }

    for k in 0..2u32 {
        let v = (s << 1) | k;
        let (start, count) = graph.arc_range(v);
        for i in 0..count {
            let arc = graph.arcs[start + i];
            if arc.deleted {
                continue;
            }
            if graph.sequences[(arc.v >> 1) as usize].coverage_class == ALTER_LABEL {
                continue;
            }
            graph.arcs[start + i].deleted = true;
            graph.delete_arc(arc.v ^ 1, v ^ 1, true);
        }
    }
}

/// `cut_trio_tip_primary`: trims the short tips and cuts the ones that come
/// from another haplotype.
///
/// A vertex starts a tip when its opposite orientation has no live outgoing
/// arc. The tip chain is collected by [`check_tip`]; chains that are too long
/// or leave the graph at an output node are skipped. An [`END_TIPS`] chain is
/// cut when [`cal_trio_vec`] assigns it a trio label other than
/// `non_trio_flag`; a [`MUL_INPUT`] chain is cut when any of the opposite
/// branches is [`NON_PLOID`]. The tip's reads are marked [`ALTER_LABEL`] and
/// their edges dropped by [`asg_seq_drop`]; cut chains are finally relabelled
/// [`CUT`].
///
/// `trio_flag` is the flag the reads were labelled with (`u32::MAX` when the
/// sample has no trio data). `max_ext` is the tip length limit in nodes, and
/// `trio_labels` must hold one label per read.
#[allow(clippy::too_many_arguments)]
pub fn cut_trio_tip_primary(
    graph: &mut Asg,
    unitigs: Option<&[MaUtg]>,
    max_ext: u32,
    trio_flag: u32,
    keep_out_node: u32,
    read_graph: &Asg,
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    is_r_het: &[u8],
    min_edge_length: u32,
    trio_labels: &[u8],
    purge_similarity_threshold: f32,
) -> u32 {
    let mut non_trio_flag = u32::MAX;
    if trio_flag == u32::from(FATHER) {
        non_trio_flag = u32::from(MOTHER);
    }
    if trio_flag == u32::from(MOTHER) {
        non_trio_flag = u32::from(FATHER);
    }

    let mut buffer: Vec<u32> = Vec::new();
    let mut b_0: Vec<u32> = Vec::new();
    let mut b_1: Vec<u32> = Vec::new();
    let mut count = 0u32;
    let vertex_count = graph.sequences.len() * 2;

    for v in 0..vertex_count as u32 {
        let sequence = graph.sequences[(v >> 1) as usize];
        if sequence.deleted || sequence.coverage_class == ALTER_LABEL {
            continue;
        }
        if sequence.coverage_class == CUT
            || sequence.coverage_class == CUT_DIF_HAP
            || sequence.coverage_class == TRIM
        {
            continue;
        }
        if get_real_length(graph, v ^ 1) != 0 {
            continue;
        }

        let mut end_node = 0u32;
        let flag = check_tip(graph, v, &mut end_node, &mut buffer, max_ext);
        if flag == LOOP || flag == LONG_TIPS {
            continue;
        }
        if keep_out_node != 0 && flag == MUL_OUTPUT {
            continue;
        }
        if let Some(unitigs) = unitigs {
            let mut tip_evaluate_len = 0u32;
            for &node in &buffer {
                tip_evaluate_len =
                    tip_evaluate_len.wrapping_add(unitigs[(node >> 1) as usize].start);
            }
            if tip_evaluate_len > max_ext {
                continue;
            }
        }

        // Two operations are possible for a tip: trimming the reads only, or
        // cutting the chain altogether.
        let mut operation = TRIM;
        if flag == END_TIPS {
            let tip_trio_flag = cal_trio_vec(&buffer, unitigs, trio_labels, TRIO_THRES);
            if (tip_trio_flag == FATHER || tip_trio_flag == MOTHER)
                && u32::from(tip_trio_flag) != non_trio_flag
            {
                operation = CUT;
            }
        }

        if flag == MUL_INPUT {
            // The tip runs into a branching node; cut it when the opposite
            // chain belongs to another haplotype.
            let last = *buffer
                .last()
                .expect("check_tip pushes at least one node");
            let opposite = get_real_length_one(graph, last)
                .expect("the tip end has exactly one live arc")
                ^ 1;
            let (start, arc_count) = graph.arc_range(opposite);
            for i in 0..arc_count {
                if operation == CUT {
                    break;
                }
                let arc = graph.arcs[start + i];
                if arc.deleted {
                    continue;
                }
                if arc.v == (last ^ 1) {
                    continue;
                }
                let inner_flag = check_different_haps(
                    graph,
                    unitigs,
                    read_graph,
                    last ^ 1,
                    arc.v,
                    reverse_sources,
                    &mut b_0,
                    &mut b_1,
                    read_to_unitig,
                    is_r_het,
                    min_edge_length,
                    1,
                    purge_similarity_threshold,
                );
                if inner_flag == NON_PLOID {
                    operation = CUT;
                }
            }
        }

        for &node in &buffer {
            graph.sequences[(node >> 1) as usize].coverage_class = ALTER_LABEL;
        }

        for &node in &buffer {
            asg_seq_drop(graph, node >> 1);
        }

        if operation == CUT {
            for &node in &buffer {
                graph.sequences[(node >> 1) as usize].coverage_class = CUT;
            }
        }

        count += 1;
    }

    if count > 0 {
        graph.cleanup();
    }

    count
}

/// `get_unitig_trio_flag`: how many reads of one unitig carry `flag`, carry
/// another label, or are ambiguous/dropped.
pub fn get_unitig_trio_flag(unitig: &MaUtg, flag: u32, trio_labels: &[u8]) -> (u32, u32, u32) {
    let mut require = 0u32;
    let mut non_require = 0u32;
    let mut ambigious = 0u32;

    for &packed in &unitig.reads {
        let read = (packed >> 33) as usize;
        let label = u32::from(trio_labels[read]);
        if label == u32::from(AMBIGU) || label == DROP {
            ambigious = ambigious.wrapping_add(1);
            continue;
        }
        if label == flag {
            require = require.wrapping_add(1);
            continue;
        }
        if label != flag {
            non_require = non_require.wrapping_add(1);
            continue;
        }
    }

    (require, non_require, ambigious)
}

/// `trio_check`: whether a tip chain is dominated by `flag` reads.
///
/// The check rejects chains that are too small, lack a clear majority, or
/// carry too few `flag` reads for their size.
pub fn trio_check(unitigs: &[MaUtg], nodes: &[u32], flag: u32, trio_labels: &[u8]) -> bool {
    if flag != u32::from(FATHER) && flag != u32::from(MOTHER) {
        return false;
    }

    let mut flag_occ = 0u32;
    let mut non_flag_occ = 0u32;
    let mut ambigious = 0u32;
    let mut u_n = 0u32;
    for &node in nodes {
        let unitig = &unitigs[(node >> 1) as usize];
        let (require, non_require, amb) = get_unitig_trio_flag(unitig, flag, trio_labels);
        flag_occ = flag_occ.wrapping_add(require);
        non_flag_occ = non_flag_occ.wrapping_add(non_require);
        ambigious = ambigious.wrapping_add(amb);
        u_n = u_n.wrapping_add(unitig.reads.len() as u32);
    }
    let _ = ambigious;

    if flag_occ.wrapping_add(non_flag_occ) == 0 {
        return false;
    }
    if f64::from(flag_occ) <= f64::from(non_flag_occ.wrapping_add(flag_occ)) * 0.75 {
        return false;
    }
    if non_flag_occ == 0 && flag_occ >= 20 {
        return true;
    }
    if u_n >= 100 {
        if f64::from(flag_occ) < f64::from(u_n) * DOUBLE_CHECK_THRES {
            return false;
        }
    } else if u_n >= 50 {
        if f64::from(flag_occ) < f64::from(u_n) * DOUBLE_CHECK_THRES * 0.5 {
            return false;
        }
    } else if f64::from(flag_occ) < f64::from(u_n) * DOUBLE_CHECK_THRES * 0.25 {
        return false;
    }

    true
}

/// `get_trio_labs`: the trio label counts of a tip chain.
pub fn get_trio_labs(nodes: &[u32], unitigs: &[MaUtg], trio_labels: &[u8]) -> TrioCounter {
    let mut counter = TrioCounter::default();
    for &node in nodes {
        let unitig = &unitigs[(node >> 1) as usize];
        for &packed in &unitig.reads {
            let read = (packed >> 33) as usize;
            let label = u32::from(trio_labels[read]);
            if label == u32::from(FATHER) {
                counter.father_occ = counter.father_occ.wrapping_add(1);
            } else if label == u32::from(MOTHER) {
                counter.mother_occ = counter.mother_occ.wrapping_add(1);
            } else if label == u32::from(AMBIGU) {
                counter.ambig_occ = counter.ambig_occ.wrapping_add(1);
            } else if label == DROP {
                counter.drop_occ = counter.drop_occ.wrapping_add(1);
            }
        }
    }
    counter.total = counter
        .father_occ
        .wrapping_add(counter.mother_occ)
        .wrapping_add(counter.ambig_occ)
        .wrapping_add(counter.drop_occ);
    counter
}

/// `renew_longest_tip_by_drop`: picks a better longest tip among the tips of a
/// branching vertex when trio labels are available.
///
/// The renewal only runs for a father/mother trio flag. It aborts as soon as
/// any candidate is not an [`END_TIPS`] chain.
#[allow(clippy::too_many_arguments)]
pub fn renew_longest_tip_by_drop(
    graph: &Asg,
    unitigs: &[MaUtg],
    arcs: &[AsgArc],
    base_max_len: &mut i64,
    base_max_len_i: &mut i64,
    stops_threshold: u32,
    buffer: &mut Vec<u32>,
    trio_flag: u32,
    trio_labels: &[u8],
) {
    if trio_flag != u32::from(FATHER) && trio_flag != u32::from(MOTHER) {
        return;
    }

    let mut max;
    let mut cur;
    let mut max_weight = 0i64;
    let mut cur_weight = 0i64;
    let mut node_len = 0i64;
    let mut base_len = 0i64;
    let mut max_stop_node_len = 0i64;
    let mut max_stop_base_len = 0i64;
    let mut end_node = 0u32;

    buffer.clear();
    let return_flag = get_unitig(
        graph,
        Some(unitigs),
        arcs[*base_max_len_i as usize].v,
        &mut end_node,
        &mut node_len,
        &mut base_len,
        &mut max_stop_node_len,
        &mut max_stop_base_len,
        stops_threshold,
        Some(buffer),
    );
    if return_flag != END_TIPS {
        return;
    }

    max = get_trio_labs(buffer, unitigs, trio_labels);

    let mut best_tip_i = *base_max_len_i as u32;
    let mut best_tip_len = *base_max_len as u32;

    for i in 0..arcs.len() {
        let index = i as u32;
        if arcs[i].deleted {
            continue;
        }
        if index == *base_max_len_i as u32 {
            continue;
        }

        buffer.clear();
        let return_flag = get_unitig(
            graph,
            Some(unitigs),
            arcs[i].v,
            &mut end_node,
            &mut node_len,
            &mut base_len,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            stops_threshold,
            Some(buffer),
        );
        if return_flag != END_TIPS {
            return;
        }
        // The alternative tip must be long enough to compete.
        if (base_len as f64) < TRIO_DROP_LENGTH_THRES * (*base_max_len as f64) {
            continue;
        }

        cur = get_trio_labs(buffer, unitigs, trio_labels);

        if trio_flag == u32::from(FATHER) {
            max_weight = i64::from(max.father_occ)
                - i64::from(max.drop_occ)
                - i64::from(max.mother_occ);
            cur_weight = i64::from(cur.father_occ)
                - i64::from(cur.drop_occ)
                - i64::from(cur.mother_occ);
            // This unitig most likely belongs to the other haplotype.
            if f64::from(cur.drop_occ.wrapping_add(cur.mother_occ))
                >= f64::from(cur.total) * TRIO_DROP_THRES
            {
                continue;
            }
        }

        if trio_flag == u32::from(MOTHER) {
            max_weight = i64::from(max.mother_occ)
                - i64::from(max.drop_occ)
                - i64::from(max.father_occ);
            cur_weight = i64::from(cur.mother_occ)
                - i64::from(cur.drop_occ)
                - i64::from(cur.father_occ);
            if f64::from(cur.drop_occ.wrapping_add(cur.father_occ))
                >= f64::from(cur.total) * TRIO_DROP_THRES
            {
                continue;
            }
        }

        if cur_weight > max_weight {
            max = cur;
            best_tip_i = index;
            best_tip_len = base_len as u32;
        } else if cur_weight == max_weight && (base_len as u32) > best_tip_len {
            max = cur;
            best_tip_i = index;
            best_tip_len = base_len as u32;
        }
    }

    *base_max_len_i = i64::from(best_tip_i);
    *base_max_len = i64::from(best_tip_len);
}

/// `asg_arc_cut_trio_long_equal_tips_assembly`: of several equally long tips at
/// one vertex, keeps the longest one and cuts the others when the two chains
/// belong to the same haplotype.
///
/// The kept tip's reads are finally relabelled [`HAP_LABEL`]. The reference
/// `utg_trans_t *o` variant is deferred, so this port implements the
/// `o == NULL` path of the primary cleaning.
#[allow(clippy::too_many_arguments)]
pub fn asg_arc_cut_trio_long_equal_tips_assembly(
    graph: &mut Asg,
    unitigs: &[MaUtg],
    read_graph: &Asg,
    reverse_sources: &[PafList],
    minimum_edge_length: i64,
    read_to_unitig: &mut ReadToUnitig,
    trio_flag: u32,
    trio_labels: &[u8],
    coverage_cut: &[CoverageCut],
    purge_similarity_threshold: f32,
    maximum_hang_rate: f32,
    cov: &mut HapCov,
) -> u32 {
    let vertex_count = graph.sequences.len() * 2;
    let mut n_reduced = 0u32;
    let mut b: Vec<u32> = Vec::new();
    let mut b_0: Vec<u32> = Vec::new();
    let mut b_1: Vec<u32> = Vec::new();
    let mut convex = 0u32;
    let mut ll = 0i64;
    let mut tmp = 0i64;
    let mut max_stop_node_len = 0i64;
    let mut max_stop_base_len = 0i64;

    for v in 0..vertex_count as u32 {
        let (arc_start, nv) = graph.arc_range(v);
        let sequence = graph.sequences[(v >> 1) as usize];
        if nv < 2 || sequence.deleted || sequence.coverage_class == ALTER_LABEL {
            continue;
        }
        let n_arc = get_real_length(graph, v);
        if n_arc < 2 {
            continue;
        }

        let mut base_max_len = -1i64;
        let mut base_max_len_i = -1i64;
        let mut n_tips = 0u32;
        let mut is_hap = 0u32;

        for i in 0..nv {
            let arc = graph.arcs[arc_start + i];
            if arc.deleted {
                continue;
            }
            let return_flag = get_unitig(
                graph,
                Some(unitigs),
                arc.v,
                &mut convex,
                &mut tmp,
                &mut ll,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                1,
                None,
            );
            if return_flag == LOOP {
                continue;
            }
            if return_flag == END_TIPS {
                n_tips += 1;
            }
            if base_max_len < ll {
                base_max_len = ll;
                base_max_len_i = i as i64;
            }
        }

        if n_tips == 0 {
            continue;
        }
        // Every live arc is a tip.
        if n_arc == n_tips {
            renew_longest_tip_by_drop(
                graph,
                unitigs,
                &graph.arcs[arc_start..arc_start + nv],
                &mut base_max_len,
                &mut base_max_len_i,
                1,
                &mut b,
                trio_flag,
                trio_labels,
            );
        }

        for i in 0..nv {
            let index = i as u32;
            if i as i64 == base_max_len_i {
                continue;
            }
            let arc = graph.arcs[arc_start + i];
            if arc.deleted {
                continue;
            }

            b.clear();
            let return_flag = get_unitig(
                graph,
                Some(unitigs),
                arc.v,
                &mut convex,
                &mut tmp,
                &mut ll,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                1,
                Some(&mut b),
            );
            if return_flag != END_TIPS {
                continue;
            }
            if trio_check(unitigs, &b, trio_flag, trio_labels) {
                continue;
            }

            let base_arc = graph.arcs[arc_start + base_max_len_i as usize];
            let flag = check_different_haps(
                graph,
                Some(unitigs),
                read_graph,
                base_arc.v,
                arc.v,
                reverse_sources,
                &mut b_0,
                &mut b_1,
                read_to_unitig,
                &cov.is_r_het,
                minimum_edge_length as u32,
                1,
                purge_similarity_threshold,
            );
            // `UNAVAILABLE`, `PLOID` and `NON_PLOID`: only `PLOID` cuts.
            if flag != PLOID {
                continue;
            }
            n_reduced = n_reduced.wrapping_add(1);

            for &node in &b {
                graph.sequences[(node >> 1) as usize].coverage_class = ALTER_LABEL;
            }
            for k in 0..b.len() {
                asg_seq_drop(graph, b[k] >> 1);
            }

            collect_trans_cov(
                &b_0,
                u64::from(base_arc.overlap),
                &b_1,
                u64::from(arc.overlap),
                unitigs,
                read_graph,
                coverage_cut,
                reverse_sources,
                maximum_hang_rate,
                cov,
            );

            is_hap = is_hap.wrapping_add(1);
            let _ = index;
        }

        if is_hap > 0 {
            let base_arc = graph.arcs[arc_start + base_max_len_i as usize];
            b.clear();
            get_unitig(
                graph,
                Some(unitigs),
                base_arc.v,
                &mut convex,
                &mut tmp,
                &mut ll,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                1,
                Some(&mut b),
            );
            for &node in &b {
                graph.sequences[(node >> 1) as usize].coverage_class = HAP_LABEL;
            }
        }
    }

    graph.cleanup();
    graph.symmetrize();

    n_reduced
}

/// `get_graph_statistic`: the number of live arcs over every live, non-alternate
/// sequence.
///
/// The reference counts `get_real_length` for both orientations of each
/// sequence, so every live arc is counted once from its query vertex.
pub fn get_graph_statistic(graph: &Asg) -> i64 {
    let mut num_arc = 0i64;
    let vertex_count = graph.sequences.len() * 2;
    for v in 0..vertex_count as u32 {
        let sequence = graph.sequences[(v >> 1) as usize];
        if sequence.deleted || sequence.coverage_class == ALTER_LABEL {
            continue;
        }
        num_arc += i64::from(get_real_length(graph, v));
    }
    num_arc
}

/// `asg_arc_cut_trio_long_tip_primary`: trims the short side branches of a
/// branching vertex when a clearly longer tip exists.
///
/// The longest chain among the live arcs is kept; every other [`END_TIPS`]
/// chain shorter than `drop_ratio` of it is dropped, and cut instead of trimmed
/// when [`check_different_haps`] reports [`NON_PLOID`]. The removed chain's
/// coverage is transferred to the kept one through [`collect_trans_cov`].
#[allow(clippy::too_many_arguments)]
pub fn asg_arc_cut_trio_long_tip_primary(
    graph: &mut Asg,
    unitigs: &[MaUtg],
    read_graph: &Asg,
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    min_edge_length: u32,
    drop_ratio: f32,
    trio_flag: u32,
    cov: &mut HapCov,
    trio_labels: &[u8],
    purge_similarity_threshold: f32,
    maximum_hang_rate: f32,
    coverage_cut: &[CoverageCut],
) -> u32 {
    let vertex_count = graph.sequences.len() * 2;
    let mut n_reduced = 0u32;
    let mut b: Vec<u32> = Vec::new();
    let mut b_0: Vec<u32> = Vec::new();
    let mut b_1: Vec<u32> = Vec::new();
    let mut convex = 0u32;
    let mut ll = 0i64;
    let mut tmp = 0i64;
    let mut max_stop_node_len = 0i64;
    let mut max_stop_base_len = 0i64;

    for v in 0..vertex_count as u32 {
        let (arc_start, nv) = graph.arc_range(v);
        let sequence = graph.sequences[(v >> 1) as usize];
        if nv < 2 || sequence.deleted || sequence.coverage_class == ALTER_LABEL {
            continue;
        }
        if get_real_length(graph, v) < 2 {
            continue;
        }

        let mut v_max_len = -1i64;
        let mut v_max_len_i = usize::MAX;
        let mut n_tips = 0u32;

        for i in 0..nv {
            let arc = graph.arcs[arc_start + i as usize];
            if arc.deleted {
                continue;
            }
            let return_flag = get_unitig(
                graph,
                Some(unitigs),
                arc.v,
                &mut convex,
                &mut ll,
                &mut tmp,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                1,
                None,
            );
            if return_flag == LOOP {
                continue;
            }
            if return_flag == END_TIPS {
                n_tips += 1;
            }
            if v_max_len < ll {
                v_max_len = ll;
                v_max_len_i = i;
            }
        }

        if n_tips == 0 {
            continue;
        }

        for i in 0..nv {
            let arc = graph.arcs[arc_start + i as usize];
            if arc.deleted {
                continue;
            }
            // The longest chain stays.
            if v_max_len_i == i {
                continue;
            }

            b.clear();
            let return_flag = get_unitig(
                graph,
                Some(unitigs),
                arc.v,
                &mut convex,
                &mut ll,
                &mut tmp,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                1,
                Some(&mut b),
            );
            if return_flag != END_TIPS {
                continue;
            }
            if (ll as f32) >= (v_max_len as f32) * drop_ratio {
                continue;
            }
            if trio_check(unitigs, &b, trio_flag, trio_labels) {
                continue;
            }
            n_reduced = n_reduced.wrapping_add(1);

            let longest_arc = graph.arcs[arc_start + v_max_len_i];
            let mut operation = TRIM;
            let flag = check_different_haps(
                graph,
                Some(unitigs),
                read_graph,
                longest_arc.v,
                arc.v,
                reverse_sources,
                &mut b_0,
                &mut b_1,
                read_to_unitig,
                &cov.is_r_het,
                min_edge_length,
                1,
                purge_similarity_threshold,
            );
            if flag == NON_PLOID {
                operation = CUT;
            }

            for k in 0..b.len() {
                graph.sequences[(b[k] >> 1) as usize].coverage_class = ALTER_LABEL;
            }
            for k in 0..b.len() {
                asg_seq_drop(graph, b[k] >> 1);
            }
            if operation == CUT {
                for k in 0..b.len() {
                    graph.sequences[(b[k] >> 1) as usize].coverage_class = CUT;
                }
            }

            if operation != CUT {
                collect_trans_cov(
                    &b_0,
                    u64::from(longest_arc.overlap),
                    &b_1,
                    u64::from(arc.overlap),
                    unitigs,
                    read_graph,
                    coverage_cut,
                    reverse_sources,
                    maximum_hang_rate,
                    cov,
                );
            }
        }
    }

    graph.cleanup();
    graph.symmetrize();

    n_reduced
}

/// `asg_arc_cut_trio_long_tip_primary_complex`: drops a short tip that enters a
/// branching node when the competing branch is much longer and has a long
/// stop-free run.
///
/// The tip chain collected from `v ^ 1` is discarded as a whole; the
/// comparison is between the tip's node length (`ll` truncated to `u32`) and
/// the candidate branch, and the stop must cover `MAX_STOP_RATE` of the
/// candidate's node length.
#[allow(clippy::too_many_arguments)]
pub fn asg_arc_cut_trio_long_tip_primary_complex(
    graph: &mut Asg,
    unitigs: &[MaUtg],
    read_graph: &Asg,
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    min_edge_length: u32,
    drop_ratio: f32,
    stops_threshold: u32,
    trio_flag: u32,
    cov: &mut HapCov,
    trio_labels: &[u8],
    purge_similarity_threshold: f32,
    maximum_hang_rate: f32,
    coverage_cut: &[CoverageCut],
) -> u32 {
    let vertex_count = graph.sequences.len() * 2;
    let mut n_reduced = 0u32;
    let mut b: Vec<u32> = Vec::new();
    let mut b_0: Vec<u32> = Vec::new();
    let mut b_1: Vec<u32> = Vec::new();
    let mut convex = 0u32;
    let mut ll = 0i64;
    let mut tmp = 0i64;
    let mut max_stop_node_len = 0i64;
    let mut max_stop_base_len = 0i64;

    for v in 0..vertex_count as u32 {
        let sequence = graph.sequences[(v >> 1) as usize];
        if sequence.deleted || sequence.coverage_class == ALTER_LABEL {
            continue;
        }
        if get_real_length(graph, v) != 0 {
            continue;
        }
        if get_real_length(graph, v ^ 1) != 1 {
            continue;
        }

        b.clear();
        let return_flag = get_unitig(
            graph,
            Some(unitigs),
            v ^ 1,
            &mut convex,
            &mut ll,
            &mut tmp,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            1,
            Some(&mut b),
        );
        if return_flag != MUL_INPUT {
            continue;
        }
        if trio_check(unitigs, &b, trio_flag, trio_labels) {
            continue;
        }

        let in_node = convex ^ 1;
        if let Some(target) = get_real_length_one(graph, convex) {
            convex = target;
        }
        convex ^= 1;
        let (convex_start, n_convex) = graph.arc_range(convex);
        let convex_len = ll as u32;

        let mut convex_i = n_convex;
        for i in 0..n_convex {
            let arc = graph.arcs[convex_start + i as usize];
            if arc.deleted {
                continue;
            }
            if arc.v == in_node {
                convex_i = i;
                break;
            }
        }
        if convex_i >= n_convex {
            continue;
        }

        for i in 0..n_convex {
            let arc = graph.arcs[convex_start + i as usize];
            if arc.deleted {
                continue;
            }
            if i == convex_i {
                continue;
            }

            let _return_flag = get_unitig(
                graph,
                Some(unitigs),
                arc.v,
                &mut convex,
                &mut ll,
                &mut tmp,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                stops_threshold,
                None,
            );

            if (convex_len as f32) < (ll as f32) * drop_ratio
                && (max_stop_node_len as f64) >= (ll as f64) * MAX_STOP_RATE
            {
                n_reduced = n_reduced.wrapping_add(1);
                let mut operation = TRIM;
                let convex_arc = graph.arcs[convex_start + convex_i as usize];
                let flag = check_different_haps(
                    graph,
                    Some(unitigs),
                    read_graph,
                    convex_arc.v,
                    arc.v,
                    reverse_sources,
                    &mut b_0,
                    &mut b_1,
                    read_to_unitig,
                    &cov.is_r_het,
                    min_edge_length,
                    stops_threshold,
                    purge_similarity_threshold,
                );
                if flag == NON_PLOID {
                    operation = CUT;
                }

                for k in 0..b.len() {
                    graph.sequences[(b[k] >> 1) as usize].coverage_class = ALTER_LABEL;
                }
                for k in 0..b.len() {
                    asg_seq_drop(graph, b[k] >> 1);
                }
                if operation == CUT {
                    for k in 0..b.len() {
                        graph.sequences[(b[k] >> 1) as usize].coverage_class = CUT;
                    }
                }

                // The reference removes the tip chain's coverage, so the
                // auxiliary side is `b_0` and the primary side `b_1`.
                if operation != CUT {
                    collect_trans_cov(
                        &b_1,
                        u64::from(arc.overlap),
                        &b_0,
                        u64::from(convex_arc.overlap),
                        unitigs,
                        read_graph,
                        coverage_cut,
                        reverse_sources,
                        maximum_hang_rate,
                        cov,
                    );
                }

                break;
            }
        }
    }

    graph.cleanup();
    graph.symmetrize();

    n_reduced
}

/// `unitig_arc_del_short_diploid_by_length_topo`: cuts short overlaps at
/// branching vertices, using the chain length on the single-input side.
///
/// Every arc whose overlap is short relative to the longest arc at both ends is
/// removed; when one end has a single live arc the decision uses the length of
/// the chain starting at the opposite vertex instead.
pub fn unitig_arc_del_short_diploid_by_length_topo(
    graph: &mut Asg,
    unitigs: &[MaUtg],
    drop_ratio: f32,
    max_ext: i32,
    if_skip_bubble: bool,
    if_drop: bool,
) -> i64 {
    let mut keys: Vec<u64> = Vec::new();
    let vertex_count = graph.sequences.len() * 2;
    let mut convex = 0u32;
    let mut tmp = 0i64;
    let mut node_len = 0i64;
    let mut max_stop_node_len = 0i64;
    let mut max_stop_base_len = 0i64;

    for v in 0..vertex_count as u32 {
        if if_skip_bubble && graph.sequence_visited[v as usize] != 0 {
            continue;
        }
        if if_drop && graph.sequences[(v >> 1) as usize].coverage_class == ALTER_LABEL {
            continue;
        }
        let (arc_start, nv) = graph.arc_range(v);
        if nv < 2 {
            continue;
        }
        if get_real_length(graph, v) < 2 {
            continue;
        }
        for i in 0..nv {
            let index = (arc_start + i as usize) as u64;
            keys.push((u64::from(graph.arcs[arc_start + i as usize].overlap) << 32) | index);
        }
    }

    radix_sort_u64(&mut keys);

    let mut n_cut = 0i64;
    for &key in &keys {
        let arc_index = key as u32 as usize;
        let a = graph.arcs[arc_index];
        // `v` is the query vertex, `w` the reverse of the target end.
        let v = a.query_vertex();
        let w = a.v ^ 1;
        let (v_start, nv) = graph.arc_range(v);
        let (w_start, nw) = graph.arc_range(w);
        if nv == 1 && nw == 1 {
            continue;
        }

        let mut ov_max = 0u32;
        let mut kv = 0u32;
        for i in 0..nv {
            let arc = graph.arcs[v_start + i as usize];
            if arc.deleted {
                continue;
            }
            if ov_max < arc.overlap {
                ov_max = arc.overlap;
            }
            kv += 1;
        }
        if kv >= 2 && (a.overlap as f32) > (ov_max as f32) * drop_ratio {
            continue;
        }

        let mut ow_max = 0u32;
        let mut kw = 0u32;
        for i in 0..nw {
            let arc = graph.arcs[w_start + i as usize];
            if arc.deleted {
                continue;
            }
            if ow_max < arc.overlap {
                ow_max = arc.overlap;
            }
            kw += 1;
        }
        if kw >= 2 && (a.overlap as f32) > (ow_max as f32) * drop_ratio {
            continue;
        }
        if kv <= 1 && kw <= 1 {
            continue;
        }

        // Locate the two directions of the same edge; both arms are always
        // present in a symmetric graph.
        let mut iv = nv;
        for i in 0..nv {
            if graph.arcs[v_start + i as usize].v == (w ^ 1) {
                iv = i;
                break;
            }
        }
        let mut iw = nw;
        for i in 0..nw {
            if graph.arcs[w_start + i as usize].v == (v ^ 1) {
                iw = i;
                break;
            }
        }
        if iv >= nv || iw >= nw {
            continue;
        }
        if graph.arcs[v_start + iv as usize].deleted
            && graph.arcs[w_start + iw as usize].deleted
        {
            continue;
        }

        let mut to_del = false;
        if kv > 1 && kw > 1 {
            if (a.overlap as f32) < (ov_max as f32) * drop_ratio
                && (a.overlap as f32) < (ow_max as f32) * drop_ratio
            {
                to_del = true;
            }
        } else if kw == 1 {
            get_unitig(
                graph,
                Some(unitigs),
                w ^ 1,
                &mut convex,
                &mut node_len,
                &mut tmp,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                1,
                None,
            );
            if node_len < i64::from(max_ext) {
                to_del = true;
            }
        } else if kv == 1 {
            get_unitig(
                graph,
                Some(unitigs),
                v ^ 1,
                &mut convex,
                &mut node_len,
                &mut tmp,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                1,
                None,
            );
            if node_len < i64::from(max_ext) {
                to_del = true;
            }
        }

        if to_del {
            graph.arcs[v_start + iv as usize].deleted = true;
            graph.arcs[w_start + iw as usize].deleted = true;
            n_cut += 1;
        }
    }

    if n_cut != 0 {
        graph.cleanup();
        graph.symmetrize();
    }

    n_cut
}

/// `unitig_arc_del_short_diploid_by_length`: removes every arc shorter than
/// `drop_ratio` of the longest live arc of its query vertex.
///
/// The mirror arc is removed through [`Asg::delete_arc`].
pub fn unitig_arc_del_short_diploid_by_length(graph: &mut Asg, drop_ratio: f32) -> i64 {
    let mut keys: Vec<u64> = Vec::new();
    let vertex_count = graph.sequences.len() * 2;

    for v in 0..vertex_count as u32 {
        let sequence = graph.sequences[(v >> 1) as usize];
        if sequence.coverage_class == ALTER_LABEL || sequence.deleted {
            continue;
        }
        let (arc_start, nv) = graph.arc_range(v);
        if nv < 2 {
            continue;
        }
        for i in 0..nv {
            let index = (arc_start + i as usize) as u64;
            keys.push((u64::from(graph.arcs[arc_start + i as usize].overlap) << 32) | index);
        }
    }

    radix_sort_u64(&mut keys);

    let mut n_cut = 0i64;
    for &key in &keys {
        let arc_index = key as u32 as usize;
        let a = graph.arcs[arc_index];
        let v = a.query_vertex();
        let (v_start, nv) = graph.arc_range(v);
        if nv <= 1 {
            continue;
        }

        let mut ov_max = 0u32;
        let mut kv = 0u32;
        for i in 0..nv {
            let arc = graph.arcs[v_start + i as usize];
            if arc.deleted {
                continue;
            }
            if ov_max < arc.overlap {
                ov_max = arc.overlap;
            }
            kv += 1;
        }
        if kv <= 1 {
            continue;
        }
        if kv >= 2 && (a.overlap as f32) > (ov_max as f32) * drop_ratio {
            continue;
        }

        graph.arcs[arc_index].deleted = true;
        graph.delete_arc(a.v ^ 1, v ^ 1, true);
        n_cut += 1;
    }

    if n_cut != 0 {
        graph.cleanup();
        graph.symmetrize();
    }

    n_cut
}

/// `drop_semi_circle`: cuts one arm of a circular unitig that returns to the
/// same branching vertex.
///
/// When both `v` and its reverse are part of the graph, the chain from `v ^ 1`
/// ends at a multi-input node and the chain from each live arc of `v` reaches
/// the same node, the arc is deleted if the two chains are [`PLOID`].
#[allow(clippy::too_many_arguments)]
pub fn drop_semi_circle(
    unitig_graph: &mut UnitigGraph,
    read_graph: &Asg,
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    is_r_het: &[u8],
    purge_similarity_threshold: f32,
) {
    let UnitigGraph { unitigs, graph } = unitig_graph;
    let vertex_count = graph.sequences.len() * 2;
    let mut b_0: Vec<u32> = Vec::new();
    let mut b_1: Vec<u32> = Vec::new();
    let mut convex_b = 0u32;
    let mut convex_f = 0u32;
    let mut ll = 0i64;
    let mut tmp = 0i64;
    let mut max_stop_node_len = 0i64;
    let mut max_stop_base_len = 0i64;

    for v in 0..vertex_count as u32 {
        if get_real_length(graph, v) == 0 {
            continue;
        }
        if get_real_length(graph, v ^ 1) == 0 {
            continue;
        }
        if get_unitig(
            graph,
            Some(unitigs.as_slice()),
            v ^ 1,
            &mut convex_b,
            &mut ll,
            &mut tmp,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            1,
            None,
        ) != MUL_INPUT
        {
            continue;
        }
        // `get_real_length(nsg, convex_b, &convex_b)` leaves the variable
        // untouched when the vertex has no live arc.
        if let Some(target) = get_real_length_one(graph, convex_b) {
            convex_b = target;
        }

        let (arc_start, nv) = graph.arc_range(v);
        for i in 0..nv {
            let arc = graph.arcs[arc_start + i as usize];
            if arc.deleted {
                continue;
            }
            if get_unitig(
                graph,
                Some(unitigs.as_slice()),
                arc.v,
                &mut convex_f,
                &mut ll,
                &mut tmp,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                1,
                None,
            ) != MUL_INPUT
            {
                continue;
            }
            if let Some(target) = get_real_length_one(graph, convex_f) {
                convex_f = target;
            }
            if convex_f != convex_b {
                continue;
            }
            if check_different_haps(
                graph,
                Some(unitigs.as_slice()),
                read_graph,
                v ^ 1,
                arc.v,
                reverse_sources,
                &mut b_0,
                &mut b_1,
                read_to_unitig,
                is_r_het,
                2,
                1,
                purge_similarity_threshold,
            ) == PLOID
            {
                graph.arcs[arc_start + i as usize].deleted = true;
                graph.delete_arc(arc.v ^ 1, v ^ 1, true);
            }
        }
    }
}

/// `detect_chimeric_by_topo`: removes unitigs that form a chimeric connection
/// between two branching vertices.
///
/// The chain from `v_beg ^ 1` must end at a single-input node on both sides,
/// and every alternative arm must point back at the same branch while still
/// belonging to another haplotype. Chains with at most
/// [`CHIMERIC_TRIM_THRES`] reads are relabelled [`CUT`]; the others are only
/// dropped from the graph.
#[allow(clippy::too_many_arguments)]
pub fn detect_chimeric_by_topo(
    graph: &mut Asg,
    unitigs: &[MaUtg],
    read_graph: &Asg,
    reverse_sources: &[PafList],
    min_edge_length: i64,
    stops_threshold: u32,
    drop_rate: f32,
    read_to_unitig: &mut ReadToUnitig,
    is_r_het: &[u8],
    purge_similarity_threshold: f32,
) -> u32 {
    let vertex_count = graph.sequences.len() * 2;
    let mut n_reduced = 0u32;
    let mut b_0: Vec<u32> = Vec::new();
    let mut b_1: Vec<u32> = Vec::new();
    let mut convex = 0u32;
    let mut ll = 0i64;
    let mut tmp = 0i64;
    let mut max_stop_node_len = 0i64;
    let mut max_stop_base_len = 0i64;
    // The reference increments its `stops_threshold` parameter, so the value
    // grows across the outer loop.
    let mut stops = stops_threshold;

    for v_beg in 0..vertex_count as u32 {
        let sequence = graph.sequences[(v_beg >> 1) as usize];
        if sequence.deleted || sequence.coverage_class == ALTER_LABEL {
            continue;
        }
        if get_real_length(graph, v_beg) != 1 {
            continue;
        }
        let w1 = match get_real_length_one(graph, v_beg) {
            Some(target) => target,
            None => continue,
        };
        if get_real_length(graph, w1 ^ 1) <= 1 {
            continue;
        }

        let mut v_end = 0u32;
        if get_unitig(
            graph,
            Some(unitigs),
            v_beg ^ 1,
            &mut v_end,
            &mut ll,
            &mut tmp,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            1,
            None,
        ) == LOOP
        {
            continue;
        }
        // `selfLen` is a `uint32_t` in the reference.
        let self_len = ll as u32;
        if get_real_length(graph, v_end) != 1 {
            continue;
        }
        let w2 = match get_real_length_one(graph, v_end) {
            Some(target) => target,
            None => continue,
        };
        if get_real_length(graph, w2 ^ 1) <= 1 {
            continue;
        }

        get_unitig(
            graph,
            Some(unitigs),
            w1,
            &mut convex,
            &mut ll,
            &mut tmp,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            stops,
            None,
        );
        if (ll as f32) * drop_rate < (self_len as f32) {
            continue;
        }
        if (ll as f64) * 0.4 > (max_stop_node_len as f64) {
            continue;
        }

        get_unitig(
            graph,
            Some(unitigs),
            w2,
            &mut convex,
            &mut ll,
            &mut tmp,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            stops,
            None,
        );
        if (ll as f32) * drop_rate < (self_len as f32) {
            continue;
        }
        if (ll as f64) * 0.4 > (max_stop_node_len as f64) {
            continue;
        }

        stops += 1;

        // Side A: every live arm of `w1 ^ 1` must point back at `v_beg`.
        let w1_rev = w1 ^ 1;
        let wv = v_beg ^ 1;
        let (aw_start, nw) = graph.arc_range(w1_rev);
        let mut convex_t = u32::MAX;
        let mut aborted = false;
        for i in 0..nw {
            let arc = graph.arcs[aw_start + i as usize];
            if arc.deleted {
                continue;
            }
            get_unitig(
                graph,
                Some(unitigs),
                arc.v,
                &mut convex,
                &mut ll,
                &mut tmp,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                stops,
                None,
            );
            if (ll as f32) * drop_rate < (self_len as f32) {
                aborted = true;
                break;
            }
            if (ll as f64) * 0.4 > (max_stop_node_len as f64) {
                continue;
            }
            if (arc.v >> 1) == (v_beg >> 1) {
                if convex_t != u32::MAX {
                    aborted = true;
                    break;
                }
                convex_t = convex;
            }
        }
        if aborted {
            continue;
        }

        let mut aborted = false;
        for i in 0..nw {
            let arc = graph.arcs[aw_start + i as usize];
            if arc.deleted {
                continue;
            }
            if (arc.v >> 1) == (v_beg >> 1) {
                continue;
            }

            b_0.clear();
            get_unitig(
                graph,
                Some(unitigs),
                arc.v,
                &mut convex,
                &mut ll,
                &mut tmp,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                stops,
                Some(&mut b_0),
            );
            let mut reached_branch = false;
            for &node in &b_0 {
                if (node >> 1) == (convex_t >> 1) {
                    reached_branch = true;
                    break;
                }
            }
            if reached_branch {
                aborted = true;
                break;
            }

            if check_different_haps(
                graph,
                Some(unitigs),
                read_graph,
                wv,
                arc.v,
                reverse_sources,
                &mut b_0,
                &mut b_1,
                read_to_unitig,
                is_r_het,
                min_edge_length as u32,
                stops,
                purge_similarity_threshold,
            ) == PLOID
            {
                aborted = true;
                break;
            }
        }
        if aborted {
            continue;
        }

        // Side B: the same checks around `v_end`.
        let w2_rev = w2 ^ 1;
        let wv = v_end ^ 1;
        let (aw_start, nw) = graph.arc_range(w2_rev);
        let mut convex_t = u32::MAX;
        let mut aborted = false;
        for i in 0..nw {
            let arc = graph.arcs[aw_start + i as usize];
            if arc.deleted {
                continue;
            }
            get_unitig(
                graph,
                Some(unitigs),
                arc.v,
                &mut convex,
                &mut ll,
                &mut tmp,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                stops,
                None,
            );
            if (ll as f32) * drop_rate < (self_len as f32) {
                aborted = true;
                break;
            }
            if (ll as f64) * 0.4 > (max_stop_node_len as f64) {
                continue;
            }
            if (arc.v >> 1) == (v_end >> 1) {
                if convex_t != u32::MAX {
                    aborted = true;
                    break;
                }
                convex_t = convex;
            }
        }
        if aborted {
            continue;
        }

        let mut aborted = false;
        for i in 0..nw {
            let arc = graph.arcs[aw_start + i as usize];
            if arc.deleted {
                continue;
            }
            if (arc.v >> 1) == (v_end >> 1) {
                continue;
            }

            b_0.clear();
            get_unitig(
                graph,
                Some(unitigs),
                arc.v,
                &mut convex,
                &mut ll,
                &mut tmp,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                stops,
                Some(&mut b_0),
            );
            let mut reached_branch = false;
            for &node in &b_0 {
                if (node >> 1) == (convex_t >> 1) {
                    reached_branch = true;
                    break;
                }
            }
            if reached_branch {
                aborted = true;
                break;
            }

            if check_different_haps(
                graph,
                Some(unitigs),
                read_graph,
                wv,
                arc.v,
                reverse_sources,
                &mut b_0,
                &mut b_1,
                read_to_unitig,
                is_r_het,
                min_edge_length as u32,
                stops,
                purge_similarity_threshold,
            ) == PLOID
            {
                aborted = true;
                break;
            }
        }
        if aborted {
            continue;
        }

        n_reduced = n_reduced.wrapping_add(1);
        b_0.clear();
        let mut read_num = 0u32;
        get_unitig(
            graph,
            Some(unitigs),
            v_beg ^ 1,
            &mut v_end,
            &mut ll,
            &mut tmp,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            1,
            Some(&mut b_0),
        );
        for k in 0..b_0.len() {
            graph.sequences[(b_0[k] >> 1) as usize].coverage_class = ALTER_LABEL;
            read_num = read_num
                .wrapping_add(unitigs[(b_0[k] >> 1) as usize].reads.len() as u32);
        }
        for k in 0..b_0.len() {
            asg_seq_drop(graph, b_0[k] >> 1);
        }
        if read_num <= CHIMERIC_TRIM_THRES {
            for k in 0..b_0.len() {
                graph.sequences[(b_0[k] >> 1) as usize].coverage_class = CUT;
            }
        }
    }

    graph.cleanup();

    n_reduced
}

/// `asg_arc_cut_trio_long_equal_tips_assembly_complex`: of two tips of nearly
/// equal length that meet at one branching node, drops the shorter one when the
/// two chains belong to the same haplotype.
///
/// The comparison uses the base length (`ll`), and the kept chain is finally
/// relabelled [`HAP_LABEL`].
#[allow(clippy::too_many_arguments)]
pub fn asg_arc_cut_trio_long_equal_tips_assembly_complex(
    graph: &mut Asg,
    unitigs: &[MaUtg],
    read_graph: &Asg,
    reverse_sources: &[PafList],
    min_edge_length: i64,
    read_to_unitig: &mut ReadToUnitig,
    stops_threshold: u32,
    trio_flag: u32,
    cov: &mut HapCov,
    trio_labels: &[u8],
    purge_similarity_threshold: f32,
    maximum_hang_rate: f32,
    coverage_cut: &[CoverageCut],
) -> u32 {
    let vertex_count = graph.sequences.len() * 2;
    let mut n_reduced = 0u32;
    let mut b: Vec<u32> = Vec::new();
    let mut b_0: Vec<u32> = Vec::new();
    let mut b_1: Vec<u32> = Vec::new();
    let mut convex = 0u32;
    let mut ll = 0i64;
    let mut tmp = 0i64;
    let mut max_stop_node_len = 0i64;
    let mut max_stop_base_len = 0i64;

    for v in 0..vertex_count as u32 {
        let sequence = graph.sequences[(v >> 1) as usize];
        if sequence.deleted || sequence.coverage_class == ALTER_LABEL {
            continue;
        }
        if get_real_length(graph, v) != 0 {
            continue;
        }
        if get_real_length(graph, v ^ 1) != 1 {
            continue;
        }

        b.clear();
        let return_flag = get_unitig(
            graph,
            Some(unitigs),
            v ^ 1,
            &mut convex,
            &mut tmp,
            &mut ll,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            1,
            Some(&mut b),
        );
        if return_flag != MUL_INPUT {
            continue;
        }
        if trio_check(unitigs, &b, trio_flag, trio_labels) {
            continue;
        }

        let in_node = convex ^ 1;
        if let Some(target) = get_real_length_one(graph, convex) {
            convex = target;
        }
        convex ^= 1;
        let (convex_start, n_convex) = graph.arc_range(convex);
        let convex_len = ll as u32;

        let mut convex_i = n_convex;
        for i in 0..n_convex {
            let arc = graph.arcs[convex_start + i as usize];
            if arc.deleted {
                continue;
            }
            if arc.v == in_node {
                convex_i = i;
                break;
            }
        }
        if convex_i >= n_convex {
            continue;
        }

        for i in 0..n_convex {
            let arc = graph.arcs[convex_start + i as usize];
            if arc.deleted {
                continue;
            }
            if i == convex_i {
                continue;
            }

            let _return_flag = get_unitig(
                graph,
                Some(unitigs),
                arc.v,
                &mut convex,
                &mut tmp,
                &mut ll,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                stops_threshold,
                None,
            );

            if ll > i64::from(convex_len)
                && (max_stop_base_len as f64) >= (ll as f64) * MAX_STOP_RATE
            {
                let convex_arc = graph.arcs[convex_start + convex_i as usize];
                let flag = check_different_haps(
                    graph,
                    Some(unitigs),
                    read_graph,
                    convex_arc.v,
                    arc.v,
                    reverse_sources,
                    &mut b_0,
                    &mut b_1,
                    read_to_unitig,
                    &cov.is_r_het,
                    min_edge_length as u32,
                    stops_threshold,
                    purge_similarity_threshold,
                );
                if flag != PLOID {
                    continue;
                }
                n_reduced = n_reduced.wrapping_add(1);

                for k in 0..b.len() {
                    graph.sequences[(b[k] >> 1) as usize].coverage_class = ALTER_LABEL;
                }
                for k in 0..b.len() {
                    asg_seq_drop(graph, b[k] >> 1);
                }

                collect_trans_cov(
                    &b_1,
                    u64::from(arc.overlap),
                    &b_0,
                    u64::from(convex_arc.overlap),
                    unitigs,
                    read_graph,
                    coverage_cut,
                    reverse_sources,
                    maximum_hang_rate,
                    cov,
                );

                // Label the primary side.
                b_0.clear();
                get_unitig(
                    graph,
                    Some(unitigs),
                    arc.v,
                    &mut convex,
                    &mut tmp,
                    &mut ll,
                    &mut max_stop_node_len,
                    &mut max_stop_base_len,
                    stops_threshold,
                    Some(&mut b_0),
                );
                for k in 0..b_0.len() {
                    graph.sequences[(b_0[k] >> 1) as usize].coverage_class = HAP_LABEL;
                }

                break;
            }
        }
    }

    graph.cleanup();
    graph.symmetrize();

    n_reduced
}

/// `clean_primary_untig_graph`: the iterative topology clean of the primary
/// unitig graph.
///
/// The body is repeated up to three times (`T_ROUND = 2`); the first pass uses
/// the simple bubble pop and tip trims, while the later passes also cut chimeric
/// connections and, between rounds, overlapping short diploid arcs.
#[allow(clippy::too_many_arguments)]
pub fn clean_primary_untig_graph(
    unitig_graph: &mut UnitigGraph,
    read_graph: &Asg,
    reverse_sources: &[PafList],
    coverage_cut: &[CoverageCut],
    tips_len: i64,
    tip_drop_ratio: f32,
    stops_threshold: u32,
    read_to_unitig: &mut ReadToUnitig,
    chimeric_rate: f32,
    just_bubble_pop: bool,
    drop_ratio: f32,
    cov: &mut HapCov,
    trio_labels: &[u8],
    max_short_tip: i32,
    purge_similarity_threshold: f32,
    maximum_hang_rate: f32,
) {
    const T_ROUND: i32 = 2;
    let mut round = T_ROUND;

    loop {
        asg_pop_bubble_primary_trio(
            unitig_graph,
            None,
            u32::MAX,
            DROP,
            Some(&mut *cov),
            trio_labels,
            true,
            None,
        );
        if !just_bubble_pop {
            cut_trio_tip_primary(
                &mut unitig_graph.graph,
                Some(&unitig_graph.unitigs),
                tips_len as u32,
                u32::MAX,
                0,
                read_graph,
                reverse_sources,
                read_to_unitig,
                &cov.is_r_het,
                2,
                trio_labels,
                purge_similarity_threshold,
            );
        }

        let mut pre_cons = get_graph_statistic(&unitig_graph.graph);
        let mut cur_cons = 0i64;
        while pre_cons != cur_cons {
            pre_cons = get_graph_statistic(&unitig_graph.graph);
            asg_pop_bubble_primary_trio(
                unitig_graph,
                None,
                u32::MAX,
                DROP,
                Some(&mut *cov),
                trio_labels,
                true,
                None,
            );
            if !just_bubble_pop {
                // Tangles have to be considered here.
                asg_arc_cut_trio_long_tip_primary(
                    &mut unitig_graph.graph,
                    &unitig_graph.unitigs,
                    read_graph,
                    reverse_sources,
                    read_to_unitig,
                    2,
                    tip_drop_ratio,
                    u32::MAX,
                    cov,
                    trio_labels,
                    purge_similarity_threshold,
                    maximum_hang_rate,
                    coverage_cut,
                );
                asg_arc_cut_trio_long_equal_tips_assembly(
                    &mut unitig_graph.graph,
                    &unitig_graph.unitigs,
                    read_graph,
                    reverse_sources,
                    2,
                    read_to_unitig,
                    u32::MAX,
                    trio_labels,
                    coverage_cut,
                    purge_similarity_threshold,
                    maximum_hang_rate,
                    cov,
                );
                asg_arc_cut_trio_long_tip_primary_complex(
                    &mut unitig_graph.graph,
                    &unitig_graph.unitigs,
                    read_graph,
                    reverse_sources,
                    read_to_unitig,
                    2,
                    tip_drop_ratio,
                    stops_threshold,
                    u32::MAX,
                    cov,
                    trio_labels,
                    purge_similarity_threshold,
                    maximum_hang_rate,
                    coverage_cut,
                );
                asg_arc_cut_trio_long_equal_tips_assembly_complex(
                    &mut unitig_graph.graph,
                    &unitig_graph.unitigs,
                    read_graph,
                    reverse_sources,
                    2,
                    read_to_unitig,
                    stops_threshold,
                    u32::MAX,
                    cov,
                    trio_labels,
                    purge_similarity_threshold,
                    maximum_hang_rate,
                    coverage_cut,
                );
                detect_chimeric_by_topo(
                    &mut unitig_graph.graph,
                    &unitig_graph.unitigs,
                    read_graph,
                    reverse_sources,
                    2,
                    stops_threshold,
                    chimeric_rate,
                    read_to_unitig,
                    &cov.is_r_het,
                    purge_similarity_threshold,
                );
                if round != T_ROUND {
                    unitig_arc_del_short_diploid_by_length_topo(
                        &mut unitig_graph.graph,
                        &unitig_graph.unitigs,
                        drop_ratio,
                        max_short_tip,
                        false,
                        true,
                    );
                }
            }
            cur_cons = get_graph_statistic(&unitig_graph.graph);
        }

        if !just_bubble_pop {
            cut_trio_tip_primary(
                &mut unitig_graph.graph,
                Some(&unitig_graph.unitigs),
                tips_len as u32,
                u32::MAX,
                0,
                read_graph,
                reverse_sources,
                read_to_unitig,
                &cov.is_r_het,
                2,
                trio_labels,
                purge_similarity_threshold,
            );
        }
        resolve_tangles(
            unitig_graph,
            read_graph,
            reverse_sources,
            20,
            100,
            0.05,
            0.2,
            read_to_unitig,
            &cov.is_r_het,
            u32::MAX,
            drop_ratio,
            trio_labels,
            purge_similarity_threshold,
        );
        drop_semi_circle(
            unitig_graph,
            read_graph,
            reverse_sources,
            read_to_unitig,
            &cov.is_r_het,
            purge_similarity_threshold,
        );
        unitig_arc_del_short_diploid_by_length_topo(
            &mut unitig_graph.graph,
            &unitig_graph.unitigs,
            drop_ratio,
            max_short_tip,
            false,
            true,
        );

        if round > 0 {
            if round != T_ROUND {
                unitig_arc_del_short_diploid_by_length(&mut unitig_graph.graph, drop_ratio);
            }
            round -= 1;
            continue;
        }
        break;
    }
}

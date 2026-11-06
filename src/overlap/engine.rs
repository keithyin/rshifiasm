//! Per-read overlap computation, the counterpart of `ha_get_new_candidates`
//! plus the final sort of `ha_get_candidates_interface`.
//!
//! The pipeline is: sketch the read, look every minimizer up in the position
//! index, turn the hits into anchors, chain them per candidate read and strand,
//! then keep the strongest chains. Nothing here refines the chains yet; that is
//! the job of the alignment stage.

use crate::container::ksort::introsort_by;
use crate::counting::filter_table::KmerFilterTable;
use crate::counting::minimizer::sketch;
use crate::counting::position_index::PositionIndex;
use crate::overlap::candidates::collect_seeds;
use crate::overlap::chaining::chain_candidate_groups;
use crate::overlap::region::ChainData;
use crate::overlap::region::OverlapRegion;

/// Everything the overlap computation needs besides the read itself.
pub struct CandidateParameters<'a> {
    /// Minimizer window, `asm_opt.mz_win`.
    pub minimizer_window: i32,
    /// K-mer length, `asm_opt.k_mer_length`.
    pub kmer_length: i32,
    /// `asm_opt.mz_sample_dist`.
    pub sample_distance: i32,
    /// `asm_opt.mz_rewin`.
    pub rewindow: i32,
    /// Whether homopolymer runs are collapsed before sketching.
    pub use_homopolymer_compression: bool,
    /// `asm_opt.hom_cov`, used to weight seed occurrences.
    pub homozygous_coverage: i32,
    /// `asm_opt.max_n_chain`: maximum number of chains per containment class.
    pub maximum_chain_count: i32,
    /// Band width of the chaining step.
    pub band_width_threshold: f64,
    /// Lengths of every read, indexed by identifier.
    pub read_lengths: &'a [u64],
    /// Position index of all minimizers.
    pub index: &'a PositionIndex,
    /// High-occurrence filter, when one was built.
    pub filter_table: Option<&'a KmerFilterTable>,
}

/// `ha_ov_type`: distinguishes the four ways a chain can sit on its query.
///
/// * `0` starts at the query start,
/// * `1` ends at the query end,
/// * `2` is fully contained in the query,
/// * `3` fully contains the query.
pub fn overlap_type(region: &OverlapRegion, read_length: u32) -> usize {
    if region.x_position_start == 0 && region.x_position_end == read_length - 1 {
        2
    } else if region.x_position_start > 0 && region.x_position_end < read_length - 1 {
        3
    } else if region.x_position_start == 0 {
        0
    } else {
        1
    }
}

/// Computes the candidate overlaps of one read.
///
/// The result is sorted by query start and end, which is the order
/// `ha_get_candidates_interface` leaves behind.
pub fn compute_candidate_overlaps(
    sequence: &[u8],
    read_identifier: u64,
    parameters: &CandidateParameters,
    chain_data: &mut ChainData,
) -> Vec<OverlapRegion> {
    let minimizers = sketch(
        sequence,
        parameters.minimizer_window,
        parameters.kmer_length,
        read_identifier as u32,
        parameters.use_homopolymer_compression,
        parameters.filter_table,
        parameters.sample_distance,
        parameters.rewindow,
    )
    .minimizers;

    let seeds = collect_seeds(
        &minimizers,
        parameters.index,
        sequence.len() as u32,
        parameters.homozygous_coverage,
    );

    let mut overlaps: Vec<OverlapRegion> = Vec::new();
    chain_candidate_groups(
        &seeds,
        read_identifier,
        sequence.len() as u64,
        parameters.read_lengths,
        parameters.band_width_threshold,
        parameters.kmer_length as i64,
        true,
        chain_data,
        &mut overlaps,
    );

    apply_chain_limit(&mut overlaps, sequence.len() as u32, parameters.maximum_chain_count);

    // `ha_get_candidates_interface` ends with `ks_introsort_or_xs`: the overlap
    // list must be ordered by query start, because `get_interval` walks it with
    // a monotone cursor.
    //
    // `RSHIFIASM_OVL_SORT=yid` additionally breaks ties by target identifier.
    // It exists to test whether the remaining correction differences come from
    // the order in which equally positioned overlaps are folded into the POA
    // graph; it must not be used for the reference comparison itself.
    let tie_break_by_target = std::env::var("RSHIFIASM_OVL_SORT").as_deref() == Ok("yid");
    if tie_break_by_target {
        introsort_by(&mut overlaps, &|left: &OverlapRegion, right: &OverlapRegion| {
            (left.x_position_start, left.y_identifier)
                < (right.x_position_start, right.y_identifier)
        });
    } else {
        introsort_by(&mut overlaps, &|left: &OverlapRegion, right: &OverlapRegion| {
            left.x_position_start < right.x_position_start
        });
    }
    overlaps
}

/// `ha_get_new_candidates` tail: keeps the strongest
/// `maximum_chain_count` chains of every containment class.
fn apply_chain_limit(overlaps: &mut Vec<OverlapRegion>, read_length: u32, maximum_chain_count: i32) {
    if maximum_chain_count <= 0 || overlaps.len() <= maximum_chain_count as usize {
        return;
    }

    // Chains are visited by decreasing seed support, so the first
    // `maximum_chain_count` of a class are the ones worth keeping.
    introsort_by(overlaps, &|left: &OverlapRegion, right: &OverlapRegion| {
        left.shared_seed > right.shared_seed
    });

    let mut class_counts = [0i32; 4];
    let mut class_thresholds = [0i32; 4];
    for region in overlaps.iter() {
        let class = overlap_type(region, read_length);
        class_counts[class] += 1;
        if class_counts[class] == maximum_chain_count {
            class_thresholds[class] = region.shared_seed;
        }
    }

    if class_thresholds.iter().all(|threshold| *threshold == 0) {
        return;
    }

    let mut kept = 0usize;
    for index in 0..overlaps.len() {
        let class = overlap_type(&overlaps[index], read_length);
        if overlaps[index].shared_seed >= class_thresholds[class] {
            overlaps.swap(kept, index);
            kept += 1;
        }
    }
    overlaps.truncate(kept);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region(start: u32, end: u32, seed: i32) -> OverlapRegion {
        OverlapRegion {
            x_position_start: start,
            x_position_end: end,
            shared_seed: seed,
            ..OverlapRegion::default()
        }
    }

    #[test]
    fn overlap_types_follow_the_containment_classes() {
        let read_length = 100u32;
        assert_eq!(overlap_type(&region(0, 99, 1), read_length), 2);
        assert_eq!(overlap_type(&region(10, 80, 1), read_length), 3);
        assert_eq!(overlap_type(&region(0, 50, 1), read_length), 0);
        assert_eq!(overlap_type(&region(20, 99, 1), read_length), 1);
    }

    #[test]
    fn chain_limit_keeps_the_strongest_chains_of_each_class() {
        // Class 0 (starting at the read start) has three members, the others
        // one each; the limit of two must drop the weakest class 0 chain.
        let mut overlaps = vec![
            region(0, 50, 30),
            region(0, 60, 40),
            region(0, 70, 10),
            region(20, 99, 25),
            region(10, 80, 5),
        ];
        apply_chain_limit(&mut overlaps, 100, 2);
        assert_eq!(overlaps.len(), 4);
        assert!(overlaps.iter().all(|item| item.shared_seed != 10));
    }

    #[test]
    fn chain_limit_is_a_no_op_below_the_threshold() {
        let mut overlaps = vec![region(0, 50, 30), region(20, 99, 25)];
        apply_chain_limit(&mut overlaps, 100, 5);
        assert_eq!(overlaps.len(), 2);
    }

    /// `ha_get_candidates_interface` sorts by query start before returning, so
    /// the window cursor of `get_interval` can walk the list once.
    #[test]
    fn candidates_are_returned_in_query_order() {
        // Chains are produced grouped by target read, so the raw order is by
        // target identifier; the final sort must restore query order.
        let mut overlaps: Vec<OverlapRegion> = vec![
            OverlapRegion { x_position_start: 900, x_position_end: 999, y_identifier: 1, ..OverlapRegion::default() },
            OverlapRegion { x_position_start: 0, x_position_end: 999, y_identifier: 2, ..OverlapRegion::default() },
            OverlapRegion { x_position_start: 450, x_position_end: 999, y_identifier: 3, ..OverlapRegion::default() },
        ];
        introsort_by(&mut overlaps, &|left: &OverlapRegion, right: &OverlapRegion| {
            left.x_position_start < right.x_position_start
        });
        let starts: Vec<u32> = overlaps.iter().map(|region| region.x_position_start).collect();
        assert_eq!(starts, vec![0, 450, 900]);
    }
}

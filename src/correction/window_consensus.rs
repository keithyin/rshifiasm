//! Building one window's consensus and walking the whole read, the
//! counterparts of `window_consensus` and the main loop of
//! `generate_consensus` in `Correct.cpp`.
//!
//! The second-round boundary correction is applied at every window boundary
//! through [`crate::correction::second_round::process_boundary`];
//! [`generate_consensus`] walks the whole read and produces the corrected
//! sequence from the verified window cigars.

use crate::correction::consensus::{CorrectedRead, get_available_interval};
use crate::correction::consensus_graph::{GraphConsensusContext, get_sequence_from_graph};
use crate::correction::dag_con::DagCon;
use crate::correction::interval::{INTERVAL_PAST_THE_END, IntervalBuffer, WindowPool};
use crate::correction::poa::PoaGraph;
use crate::correction::region::{MINIMUM_COVERAGE_THRESHOLD, WINDOW_LENGTH};
use crate::overlap::region::OverlapRegion;
use crate::reads::read_store::ReadStore;
use std::sync::atomic::{AtomicBool, Ordering};

/// Set by the pipeline while the read selected by `RSHIFIASM_EC_READ` is being
/// corrected, so that the window layout can be printed.
static TRACED_READ: AtomicBool = AtomicBool::new(false);

/// Turns the per-window trace on or off.
pub fn set_traced_read(enabled: bool) {
    TRACED_READ.store(enabled, Ordering::Relaxed);
}

/// Whether the per-window trace is on.
pub(crate) fn tracing() -> bool {
    TRACED_READ.load(Ordering::Relaxed)
}

/// `window_consensus`: folds every usable overlap of one window into a POA
/// backbone and walks it to produce the corrected bases.
pub fn window_consensus(
    query_sequence: &[u8],
    window_start: i64,
    window_end: i64,
    overlaps: &[OverlapRegion],
    buffer: &IntervalBuffer,
    reads: &ReadStore,
    backbone_graph: &mut PoaGraph,
    dag: &mut DagCon,
    corrected: &mut CorrectedRead,
) {
    backbone_graph.clear();
    dag.graph.clear();
    dag.node_queue.clear();

    let backbone = &query_sequence[window_start as usize..=window_end as usize];
    backbone_graph.add_unmatched_sequence(backbone);

    for identifier in buffer.full_windows() {
        let region = &overlaps[*identifier as usize];
        let aligned_origin =
            (region.x_position_start as i64 / WINDOW_LENGTH) * WINDOW_LENGTH;
        let window_id = (window_start - aligned_origin) / WINDOW_LENGTH;
        if window_id < 0 || window_id as usize >= region.window_cigar.windows.len() {
            continue;
        }

        let window = region.window_cigar.windows[window_id as usize];
        if window.y_end == -1 {
            continue;
        }

        let target_start = window.y_start as i64;
        let target_length = (window.y_end + 1 - window.y_start) as i64;
        let target = reads.read_sub_region(
            region.y_identifier as usize,
            target_start,
            target_length,
            region.y_strand as u8,
        );

        let cells: Vec<(u32, u32)> = region
            .window_cigar
            .window_cells(window_id as usize)
            .map(|(operation, length)| (operation as u32, length))
            .collect();
        // `addmatchedSeqToGraph` receives the query window in C but never
        // reads it; only the target sequence and the cigar cells are used.
        let current_node = (window.x_start as i64 - window_start) as u64;
        backbone_graph.add_matched_sequence(current_node, &target, &cells);
    }

    let before = corrected.sequence.len();
    get_sequence_from_graph(
        backbone_graph,
        dag,
        corrected,
        &GraphConsensusContext {
            query_sequence,
            query_window: backbone,
            window_start,
        },
    );
    if tracing() {
        eprintln!(
            "[M::window_consensus] x=[{},{}] backbone={} emitted={} delta={}",
            window_start,
            window_end,
            backbone.len(),
            corrected.sequence.len() - before,
            (corrected.sequence.len() - before) as i64 - window_end + window_start - 1
        );
    }
}

/// Result of the first-round consensus pass.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConsensusResult {
    /// Corrected sequence of the whole read.
    pub sequence: Vec<u8>,
    /// Number of bases changed by correction.
    pub corrected_base: u64,
    /// Sum of usable overlaps seen across all windows, `num_availiable_win`.
    pub total_available_overlaps: u64,
}

/// `generate_consensus`: walks the read window by window, correcting the
/// windows with enough coverage and copying the others verbatim.
///
/// `maximum_overlap_divergence` is `asm_opt.max_ov_diff_ec`. The second-round
/// alignment is filled when a corrected window has a usable boundary; pass
/// `None` to skip that pass (for example in tests of the first round alone).
pub fn generate_consensus(
    query_sequence: &[u8],
    overlaps: &[OverlapRegion],
    reads: &ReadStore,
    maximum_overlap_divergence: f64,
    backbone_graph: &mut PoaGraph,
    dag: &mut DagCon,
    mut second_round: Option<&mut crate::correction::second_round::SecondRoundAlignment>,
) -> ConsensusResult {
    let mut corrected = CorrectedRead::new();
    let mut buffer = IntervalBuffer::with_capacity(overlaps.len().max(1));
    let mut pool = WindowPool::new(
        query_sequence.len() as i64,
        WINDOW_LENGTH,
        (1.0 / maximum_overlap_divergence) as i64,
    );
    let mut total_available_overlaps = 0u64;
    let mut flag = 0i32;
    let mut last_window_end = query_sequence.len() as i64;

    while let Some((window_start, window_end)) = pool.next_window() {
        if flag == INTERVAL_PAST_THE_END {
            break;
        }
        last_window_end = window_end;

        buffer.length = 0;
        buffer.partial_length = 0;
        flag = get_available_interval(window_start, window_end, overlaps, &mut buffer);
        total_available_overlaps += buffer.length;

        let previous_boundary = corrected.last_boundary_length;
        if buffer.length as i64 >= MINIMUM_COVERAGE_THRESHOLD as i64 {
            window_consensus(
                query_sequence,
                window_start,
                window_end,
                overlaps,
                &buffer,
                reads,
                backbone_graph,
                dag,
                &mut corrected,
            );

            if let Some(round2) = second_round.as_deref_mut() {
                if previous_boundary != 0 {
                    crate::correction::second_round::process_boundary(
                        &corrected.sequence,
                        overlaps,
                        &buffer,
                        reads,
                        window_start,
                        previous_boundary as i64,
                        maximum_overlap_divergence,
                        backbone_graph,
                        dag,
                        round2,
                    );
                }
            }
        } else {
            let segment = &query_sequence[window_start as usize..=window_end as usize];
            corrected.add_operation(segment, 0);
        }
        corrected.last_boundary_length = corrected.cigar.new_read_length as u64;
    }

    // The reference appends the tail after the last window when the final
    // `window_start` still falls inside the read.
    let tail_start = last_window_end + 1;
    if tail_start < query_sequence.len() as i64 {
        let segment = &query_sequence[tail_start as usize..];
        corrected.add_operation(segment, 0);
    }

    // The accumulated second-round cigar ends with the bases that no boundary
    // pass has claimed yet.
    if let Some(round2) = second_round {
        let remaining = corrected.cigar.new_read_length as i64 - round2.obtained_cigar_length;
        if remaining != 0 {
            let bases = vec![b'A'; remaining.max(0) as usize];
            round2.cigar.add(&bases, remaining.max(0) as u32, 0);
        }
    }

    ConsensusResult {
        sequence: corrected.sequence,
        corrected_base: corrected.corrected_base,
        total_available_overlaps,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reads() -> ReadStore {
        let mut store = ReadStore::default();
        store.push_read("query", &vec![b'A'; 1000]);
        store.push_read("target", &vec![b'A'; 1000]);
        store
    }

    #[test]
    fn a_read_without_overlaps_is_copied_verbatim() {
        let store = reads();
        let query = b"ACGTACGTACGTACGT".to_vec();
        let mut backbone = PoaGraph::new();
        let mut dag = DagCon::default();

        let result = generate_consensus(&query, &[], &store, 0.04, &mut backbone, &mut dag, None);
        assert_eq!(result.sequence, query);
        assert_eq!(result.corrected_base, 0);
        assert_eq!(result.total_available_overlaps, 0);
        assert_eq!(result.sequence.len(), query.len());
    }

    #[test]
    fn a_single_matching_overlap_does_not_change_the_sequence() {
        let store = reads();
        let query = b"ACGTACGTACGTACGT".to_vec();

        let mut region = OverlapRegion {
            x_identifier: 0,
            x_position_start: 0,
            x_position_end: 15,
            y_identifier: 1,
            y_position_start: 0,
            y_position_end: 15,
            is_match: 1,
            ..OverlapRegion::default()
        };
        region
            .window_cigar
            .append_window(0, 15, 0, 15, 0, 0, 0, 15);
        region.window_cigar.push_cell(0, 16);
        region.window_cigar.assign_cells(0, 0);

        let mut backbone = PoaGraph::new();
        let mut dag = DagCon::default();
        let mut corrected = CorrectedRead::new();
        let mut buffer = IntervalBuffer::with_capacity(4);
        buffer.overlap_identifiers[0] = 0;
        buffer.length = 1;

        window_consensus(
            &query,
            0,
            15,
            &[region],
            &buffer,
            &store,
            &mut backbone,
            &mut dag,
            &mut corrected,
        );

        assert_eq!(corrected.sequence, query);
        assert_eq!(corrected.corrected_base, 0);
    }

    /// Golden values printed by the real `generate_consensus` for reads with
    /// no overlaps: the corrected read is the input, byte for byte.
    #[test]
    fn matches_the_reference_without_overlaps() {
        let store = reads();
        for length in [1usize, 100, 375, 376, 750, 1000] {
            let query: Vec<u8> = (0..length).map(|index| b"ACGT"[index % 4]).collect();
            let mut backbone = PoaGraph::new();
            let mut dag = DagCon::default();
            let result =
                generate_consensus(&query, &[], &store, 0.04, &mut backbone, &mut dag, None);
            assert_eq!(result.sequence, query, "length {length}");
            assert_eq!(result.corrected_base, 0, "length {length}");
            assert_eq!(result.total_available_overlaps, 0);
        }
    }
}

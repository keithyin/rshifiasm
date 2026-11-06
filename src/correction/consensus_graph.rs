//! Walking the POA backbone to produce a corrected window, the counterpart of
//! `get_seq_from_Graph` and `add_path_to_correct_read_new` in `Correct.cpp`.

use crate::correction::cigar::{
    OPERATION_DELETION, OPERATION_INSERTION, OPERATION_MATCH, OPERATION_MISMATCH,
};
use crate::correction::consensus::{CorrectedRead, is_strict_homopolymer};
use crate::correction::dag_con::{DagCon, build_dag_con};
use crate::correction::poa::{Edge, PoaGraph};
use crate::correction::region::{
    CORRECTION_THRESHOLD, HOMOPOLYMER_CORRECTION_THRESHOLD,
};
use crate::sequence::encoding::SEQUENCE_NT6_TABLE;

/// Inputs that stay fixed while one window is walked.
pub struct GraphConsensusContext<'a> {
    /// The whole query read, `r_string`, used for the homopolymer test.
    pub query_sequence: &'a [u8],
    /// The query window itself, `self_string = r_string + window_start`.
    pub query_window: &'a [u8],
    /// Offset of the window in the read, `r_string_site`.
    pub window_start: i64,
}

/// `add_path_to_correct_read_new`: applies the chosen graph operation to the
/// corrected read and returns the next backbone node.
fn add_path_to_correct_read(
    backbone: &mut PoaGraph,
    dag: &mut DagCon,
    corrected: &mut CorrectedRead,
    current_node: u64,
    operation_type: i64,
    edge_index: usize,
    query_window: &[u8],
) -> u64 {
    match operation_type {
        1 => {
            // MISMATCH: a zero-length edge is the backbone match, a
            // length-one edge is an alternative base introduced by a read.
            let edge = backbone.nodes[current_node as usize].mismatch_edges.edges[edge_index];
            let alternative = edge.out_node;
            let base = backbone.nodes[alternative as usize].base;
            corrected.add_base_directly(base);

            if edge.length == 0 {
                corrected.cigar.add(&[base], 1, OPERATION_MATCH);
                return alternative;
            }

            corrected.corrected_base += 1;
            let next: Edge = backbone.nodes[alternative as usize].mismatch_edges.edges[0];
            let next_base = backbone.nodes[next.out_node as usize].base;
            let merged_base = (SEQUENCE_NT6_TABLE[base as usize] << 3)
                | SEQUENCE_NT6_TABLE[next_base as usize];
            corrected.cigar.add(&[merged_base], 1, OPERATION_MISMATCH);
            next.out_node
        }
        3 => {
            // DELETION: the backbone has bases that the consensus omits.
            let edge = backbone.nodes[current_node as usize].deletion_edges.edges[edge_index];
            let next = edge.out_node;
            let removed = (next - current_node) as usize;
            corrected.corrected_base += removed as u64;
            let start = current_node as usize;
            let mut bases = vec![b'N'; removed];
            if start + removed <= query_window.len() {
                bases.copy_from_slice(&query_window[start..start + removed]);
            }
            corrected
                .cigar
                .add(&bases, removed as u32, OPERATION_DELETION);
            next
        }
        2 => {
            // INSERTION: emit the best path of the insertion DAG.
            backbone.nodes[current_node as usize].num_insertions = 0;
            while let Some(base) = dag.node_queue.pop_front() {
                corrected.add_base_directly(base);
                corrected.cigar.add(&[base], 1, OPERATION_INSERTION);
                corrected.corrected_base += 1;
            }
            current_node
        }
        other => panic!("unknown graph operation {other}"),
    }
}

/// `get_seq_from_Graph`: walks the backbone from the start node to the end
/// node, choosing the best supported base at every position.
pub fn get_sequence_from_graph(
    backbone: &mut PoaGraph,
    dag: &mut DagCon,
    corrected: &mut CorrectedRead,
    context: &GraphConsensusContext<'_>,
) {
    let mut current = backbone.start_node_identifier;

    while current != backbone.end_node_identifier {
        let mut total_count: i64 = 0;
        let mut maximum_count: i64 = -1;
        let mut maximum_type: i64 = -1;
        let mut maximum_edge: usize = 0;

        let node = &backbone.nodes[current as usize];
        let has_insertions = node.num_insertions != 0;
        let mismatch_edges = node.mismatch_edges.edges.clone();

        for (index, edge) in mismatch_edges.iter().enumerate() {
            let weight = if has_insertions {
                edge.weight as i64 - edge.num_insertions as i64
            } else {
                edge.weight as i64
            };
            total_count += weight;
            if weight > maximum_count {
                maximum_count = weight;
                maximum_edge = index;
                maximum_type = 1;
            }
        }

        if has_insertions {
            total_count += backbone.nodes[current as usize].num_insertions as i64;
            let insertion_support = build_dag_con(dag, backbone, current) as i64;
            if insertion_support > maximum_count {
                maximum_count = insertion_support;
                maximum_type = 2;
            }
        }

        let deletion_edges = backbone.nodes[current as usize].deletion_edges.edges.clone();
        for (index, edge) in deletion_edges.iter().enumerate() {
            let weight = edge.weight as i64;
            total_count += weight;
            if weight > maximum_count {
                maximum_count = weight;
                maximum_edge = index;
                maximum_type = 3;
            }
        }

        let corrected_here = maximum_count as f64 >= total_count as f64 * CORRECTION_THRESHOLD
            || (current > 0
                && is_strict_homopolymer(
                    context.window_start + current as i64 - 1,
                    context.query_sequence,
                )
                && maximum_count as f64
                    >= total_count as f64 * HOMOPOLYMER_CORRECTION_THRESHOLD);

        if corrected_here && maximum_type != -1 {
            current = add_path_to_correct_read(
                backbone,
                dag,
                corrected,
                current,
                maximum_type,
                maximum_edge,
                context.query_window,
            );
        } else {
            current += 1;
            let base = backbone.nodes[current as usize].base;
            corrected.add_base_directly(base);
            corrected.cigar.add(&[base], 1, OPERATION_MATCH);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_backbone() -> PoaGraph {
        let mut graph = PoaGraph::new();
        graph.add_unmatched_sequence(b"ACGT");
        graph
    }

    #[test]
    fn an_uncorrected_backbone_reproduces_the_window() {
        let mut backbone = simple_backbone();
        let mut dag = DagCon::default();
        let mut corrected = CorrectedRead::new();
        let query = b"ACGT".to_vec();

        // No alternative edges: every step has one match edge of weight one,
        // so the 0.60 threshold is met and the backbone base is emitted.
        get_sequence_from_graph(
            &mut backbone,
            &mut dag,
            &mut corrected,
            &GraphConsensusContext {
                query_sequence: &query,
                query_window: &query,
                window_start: 0,
            },
        );

        // The loop visits backbone nodes 0..end-1 and emits the base of each
        // edge target, so a four-base window produces four bases.
        assert_eq!(corrected.sequence, b"ACGT");
        assert_eq!(corrected.cigar.new_read_length, 4);
    }

    #[test]
    fn a_supported_alternative_base_replaces_the_backbone_base() {
        let mut backbone = simple_backbone();
        // Add an alternative 'T' at the first backbone position with weight 9
        // against the backbone edge of weight 1.
        backbone.add_mismatch_edge_weight(1, b'T', -1);
        backbone.nodes[1].mismatch_edges.edges[1].weight = 9;

        let mut dag = DagCon::default();
        let mut corrected = CorrectedRead::new();
        let query = b"ACGT".to_vec();
        get_sequence_from_graph(
            &mut backbone,
            &mut dag,
            &mut corrected,
            &GraphConsensusContext {
                query_sequence: &query,
                query_window: &query,
                window_start: 0,
            },
        );

        // The alternative edge leaves backbone node 1 and therefore replaces
        // the second base of the window, not the first.
        assert_eq!(corrected.sequence, b"ATGT", "the alternative base wins");
        assert_eq!(corrected.corrected_base, 1);
        assert!(
            corrected
                .cigar
                .cells()
                .any(|(operation, _)| operation == OPERATION_MISMATCH),
            "a mismatch cell is recorded"
        );
    }

    /// Golden values printed by the real `get_seq_from_Graph`.
    #[test]
    fn matches_the_reference_consensus_walk() {
        struct Case {
            name: &'static str,
            alternative: bool,
            insertion: bool,
            weight: u64,
            sequence: &'static [u8],
            corrected_base: u64,
            cigar_cells: Vec<(u32, u32)>,
        }
        let cases = [
            Case {
                name: "plain",
                alternative: false,
                insertion: false,
                weight: 0,
                sequence: b"ACGT",
                corrected_base: 0,
                cigar_cells: vec![(0, 4)],
            },
            Case {
                name: "alternative",
                alternative: true,
                insertion: false,
                weight: 9,
                sequence: b"ATGT",
                corrected_base: 1,
                cigar_cells: vec![(0, 1), (1, 1), (0, 2)],
            },
            Case {
                name: "weak_alternative",
                alternative: true,
                insertion: false,
                weight: 1,
                sequence: b"ACGT",
                corrected_base: 0,
                cigar_cells: vec![(0, 4)],
            },
            Case {
                name: "insertion",
                alternative: false,
                insertion: true,
                weight: 9,
                sequence: b"AGGCGT",
                corrected_base: 2,
                cigar_cells: vec![(0, 1), (2, 2), (0, 3)],
            },
        ];

        for case in cases {
            let mut backbone = simple_backbone();
            if case.alternative {
                backbone.add_mismatch_edge_weight(1, b'T', -1);
                backbone.nodes[1].mismatch_edges.edges[1].weight = case.weight;
            }
            if case.insertion {
                backbone.add_insertion_edge_weight(1, b"GG");
                backbone.nodes[1].insertion_edges.edges[0].weight = case.weight;
                backbone.nodes[1].num_insertions = case.weight;
            }

            let mut dag = DagCon::default();
            let mut corrected = CorrectedRead::new();
            let query = b"ACGT".to_vec();
            get_sequence_from_graph(
                &mut backbone,
                &mut dag,
                &mut corrected,
                &GraphConsensusContext {
                    query_sequence: &query,
                    query_window: &query,
                    window_start: 0,
                },
            );

            assert_eq!(corrected.sequence, case.sequence, "sequence of {}", case.name);
            assert_eq!(
                corrected.corrected_base, case.corrected_base,
                "corrected base count of {}",
                case.name
            );
            let cells: Vec<(u32, u32)> = corrected.cigar.cells().collect();
            assert_eq!(cells, case.cigar_cells, "cigar of {}", case.name);
        }
    }
}

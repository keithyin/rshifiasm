//! The insertion consensus DAG, the counterpart of `build_DAGCon`,
//! `Merge_DAGCon`, `Merge_In_Nodes`, `Merge_Out_Nodes`,
//! `generate_best_seq_from_nodes` and `generate_seq_from_node`.
//!
//! When several reads insert different bases at the same query position, the
//! reference builds a small DAG out of those insertion paths, merges nodes
//! that carry the same base, and picks the highest-weight path. The chosen
//! bases are left in `node_queue` for `add_path_to_correct_read_new`.

use std::collections::VecDeque;

use crate::correction::poa::{DELETED_BASE, EdgeListKind, PoaGraph};

/// The DAG plus the queue that carries the selected insertion bases.
#[derive(Debug, Clone, Default)]
pub struct DagCon {
    pub graph: PoaGraph,
    /// `Graph::node_q`, filled by the best-path walk.
    pub node_queue: VecDeque<u8>,
}

/// Guards against unbounded recursion in the merging pass; the reference
/// relies on the graph acyclicity instead.
const MAXIMUM_MERGE_DEPTH: usize = 4096;

fn mark_edge(graph: &mut PoaGraph, in_node: u64, out_node: u64) -> Option<(usize, usize)> {
    let (forward_slot, backward_slot) = graph.bidirectional_edge_slots(in_node, out_node)?;
    graph.set_edge_visited_in_list(in_node, EdgeListKind::Deletion, forward_slot, true);
    graph.set_edge_visited_in_list(out_node, EdgeListKind::Insertion, backward_slot, true);
    Some((forward_slot, backward_slot))
}

/// Increases the weight of `in_node -> out_node` by `delta`, or creates the
/// edge (already marked visited) when it does not exist yet.
fn add_edge_weight(graph: &mut PoaGraph, in_node: u64, out_node: u64, delta: u64) {
    if let Some((forward_slot, backward_slot)) =
        graph.bidirectional_edge_slots(in_node, out_node)
    {
        mark_edge(graph, in_node, out_node);
        let weight = graph.edge_at(in_node, EdgeListKind::Deletion, forward_slot).weight + delta;
        graph.set_edge_weight_in_list(in_node, EdgeListKind::Deletion, forward_slot, weight);
        graph.set_edge_weight_in_list(out_node, EdgeListKind::Insertion, backward_slot, weight);
    } else {
        graph.add_bidirectional_edge(in_node, out_node, delta, 1);
    }
}

/// `Merge_Out_Nodes`: merges the out-nodes of `current` that carry the same
/// base and have exactly one in-edge, summing their weights.
fn merge_out_nodes(graph: &mut PoaGraph, current: u64, depth: usize) {
    if depth > MAXIMUM_MERGE_DEPTH || graph.live_edge_count(current, EdgeListKind::Deletion) == 0 {
        return;
    }

    for base in [b'A', b'C', b'G', b'T'] {
        let mut flag = 0i32;
        let mut weight = 0u64;
        let mut consensus_node = 0u64;
        let candidates = graph.live_output_nodes(current);

        for candidate in candidates {
            if graph.nodes[candidate as usize].base != base
                || graph.live_edge_count(candidate, EdgeListKind::Insertion) != 1
            {
                continue;
            }

            if flag == 0 {
                flag = 1;
                consensus_node = candidate;
                if let Some((forward_slot, _)) = mark_edge(graph, current, consensus_node) {
                    weight =
                        graph.edge_at(current, EdgeListKind::Deletion, forward_slot).weight;
                }
            } else {
                flag += 1;
                if let Some((forward_slot, _)) =
                    graph.bidirectional_edge_slots(current, candidate)
                {
                    weight +=
                        graph.edge_at(current, EdgeListKind::Deletion, forward_slot).weight;
                }

                for successor in graph.live_output_nodes(candidate) {
                    let successor_weight = graph.edge_weight(candidate, successor).unwrap_or(0);
                    add_edge_weight(graph, consensus_node, successor, successor_weight);
                }
                graph.delete_node(candidate);
            }
        }

        if flag > 1 {
            if let Some((forward_slot, backward_slot)) =
                graph.bidirectional_edge_slots(current, consensus_node)
            {
                graph.set_edge_weight_in_list(
                    current,
                    EdgeListKind::Deletion,
                    forward_slot,
                    weight,
                );
                graph.set_edge_weight_in_list(
                    consensus_node,
                    EdgeListKind::Insertion,
                    backward_slot,
                    weight,
                );
            }
        }
        if flag > 0 {
            merge_out_nodes(graph, consensus_node, depth + 1);
        }
    }
}

/// `Merge_In_Nodes`: the symmetric pass over the in-nodes of `current`.
fn merge_in_nodes(graph: &mut PoaGraph, current: u64, depth: usize) {
    if depth > MAXIMUM_MERGE_DEPTH || graph.live_edge_count(current, EdgeListKind::Insertion) == 0 {
        return;
    }

    for base in [b'A', b'C', b'G', b'T'] {
        let mut flag = 0i32;
        let mut weight = 0u64;
        let mut consensus_node = 0u64;
        let candidates = graph.live_input_nodes(current);

        for candidate in candidates {
            if graph.nodes[candidate as usize].base != base
                || graph.live_edge_count(candidate, EdgeListKind::Deletion) != 1
            {
                continue;
            }

            if flag == 0 {
                flag = 1;
                consensus_node = candidate;
                if let Some((forward_slot, _)) = mark_edge(graph, consensus_node, current) {
                    weight =
                        graph.edge_at(consensus_node, EdgeListKind::Deletion, forward_slot).weight;
                }
            } else {
                flag += 1;
                if let Some((forward_slot, _)) =
                    graph.bidirectional_edge_slots(candidate, current)
                {
                    weight +=
                        graph.edge_at(candidate, EdgeListKind::Deletion, forward_slot).weight;
                }

                for predecessor in graph.live_input_nodes(candidate) {
                    let predecessor_weight = graph.edge_weight(predecessor, candidate).unwrap_or(0);
                    add_edge_weight(graph, predecessor, consensus_node, predecessor_weight);
                }
                graph.delete_node(candidate);
            }
        }

        if flag > 1 {
            if let Some((forward_slot, backward_slot)) =
                graph.bidirectional_edge_slots(consensus_node, current)
            {
                graph.set_edge_weight_in_list(
                    consensus_node,
                    EdgeListKind::Deletion,
                    forward_slot,
                    weight,
                );
                graph.set_edge_weight_in_list(
                    current,
                    EdgeListKind::Insertion,
                    backward_slot,
                    weight,
                );
            }
        }
        if flag > 0 {
            merge_in_nodes(graph, consensus_node, depth + 1);
        }
    }
}

/// `Merge_DAGCon`: a topological walk from the start node that merges equal
/// bases and only advances to a node once all of its in-edges are visited.
pub fn merge_dag_con(dag: &mut DagCon) {
    let mut queue: VecDeque<u64> = VecDeque::new();
    queue.push_back(dag.graph.start_node_identifier);

    while let Some(current) = queue.pop_front() {
        merge_in_nodes(&mut dag.graph, current, 0);
        merge_out_nodes(&mut dag.graph, current, 0);

        for slot in dag.graph.live_output_edge_slots(current) {
            let edge = dag.graph.edge_at(current, EdgeListKind::Deletion, slot);
            mark_edge(&mut dag.graph, edge.in_node, edge.out_node);
        }

        for out_node in dag.graph.live_output_nodes(current) {
            let all_visited = dag
                .graph
                .live_input_edge_slots(out_node)
                .iter()
                .all(|slot| dag.graph.edge_visited_in_list(out_node, EdgeListKind::Insertion, *slot));
            if all_visited {
                queue.push_back(out_node);
            }
        }
    }
}

/// `generate_seq_from_node`: walks the highest-weight node path into
/// `node_queue`, reversing the result when walking backwards from the end.
fn generate_sequence_from_node(dag: &mut DagCon, mut node: u64, backwards: bool) {
    dag.node_queue.clear();

    if !backwards {
        while node != dag.graph.end_node_identifier {
            dag.node_queue.push_back(dag.graph.nodes[node as usize].base);
            let mut maximum = 0u64;
            let mut best = node;
            for out_node in dag.graph.live_output_nodes(node) {
                let weight = dag.graph.nodes[out_node as usize].weight;
                if weight > maximum {
                    maximum = weight;
                    best = out_node;
                }
            }
            if best == node {
                break;
            }
            node = best;
        }
    } else {
        while node != dag.graph.start_node_identifier {
            dag.node_queue.push_back(dag.graph.nodes[node as usize].base);
            let mut maximum = 0u64;
            let mut best = node;
            for in_node in dag.graph.live_input_nodes(node) {
                let weight = dag.graph.nodes[in_node as usize].weight;
                if weight > maximum {
                    maximum = weight;
                    best = in_node;
                }
            }
            if best == node {
                break;
            }
            node = best;
        }

        // The reference reverses the live part of the queue in place.
        let reversed: Vec<u8> = dag.node_queue.iter().rev().copied().collect();
        dag.node_queue = reversed.into();
    }
}

/// `generate_best_seq_from_nodes`: recomputes node weights, picks the best
/// start-side or end-side path and returns its support count.
pub fn generate_best_sequence(dag: &mut DagCon) -> u64 {
    for index in 0..dag.graph.nodes.len() {
        if dag.graph.nodes[index].base == DELETED_BASE {
            continue;
        }
        let weight: u64 = dag
            .graph
            .live_output_edge_slots(index as u64)
            .iter()
            .map(|slot| dag.graph.edge_at(index as u64, EdgeListKind::Deletion, *slot).weight)
            .sum();
        dag.graph.nodes[index].weight = weight;
    }

    let start = dag.graph.start_node_identifier;
    let start_weight: u64 = dag
        .graph
        .live_output_edge_slots(start)
        .iter()
        .map(|slot| dag.graph.edge_at(start, EdgeListKind::Deletion, *slot).weight)
        .sum();
    dag.graph.nodes[start as usize].weight = start_weight;

    let end = dag.graph.end_node_identifier;
    let end_weight: u64 = dag
        .graph
        .live_input_edge_slots(end)
        .iter()
        .map(|slot| dag.graph.edge_at(end, EdgeListKind::Insertion, *slot).weight)
        .sum();
    dag.graph.nodes[end as usize].weight = end_weight;

    let mut maximum_start = 0u64;
    let mut best_start = start;
    for out_node in dag.graph.live_output_nodes(start) {
        let weight = dag.graph.nodes[out_node as usize].weight;
        if weight > maximum_start {
            maximum_start = weight;
            best_start = out_node;
        }
    }

    let mut maximum_end = 0u64;
    let mut best_end = end;
    for in_node in dag.graph.live_input_nodes(end) {
        let weight = dag.graph.nodes[in_node as usize].weight;
        if weight > maximum_end {
            maximum_end = weight;
            best_end = in_node;
        }
    }

    if maximum_start == 0 && maximum_end == 0 {
        // The reference only reaches this function when the anchor has at
        // least one insertion bubble; without a path it would dereference a
        // null node. Returning no sequence is the safe equivalent.
        dag.node_queue.clear();
        return 0;
    }

    if maximum_start >= maximum_end {
        generate_sequence_from_node(dag, best_start, false);
        maximum_start
    } else {
        generate_sequence_from_node(dag, best_end, true);
        maximum_end
    }
}

/// `build_DAGCon`: turns the insertion bubbles anchored at `current_node` of
/// `backbone` into a merged DAG and leaves the best insertion sequence in
/// `node_queue`.
///
/// Returns the support count of that sequence. The reference also clears the
/// anchor's insertion counter, which is done here as well.
pub fn build_dag_con(
    dag: &mut DagCon,
    backbone: &mut PoaGraph,
    current_node: u64,
) -> u64 {
    dag.graph.clear();
    dag.node_queue.clear();

    let start = dag.graph.add_node(b'S');
    dag.graph.start_node_identifier = start;
    let end = dag.graph.add_node(b'E');
    dag.graph.end_node_identifier = end;

    let bubbles: Vec<(u64, u64, u64)> = backbone.nodes[current_node as usize]
        .insertion_edges
        .edges
        .iter()
        .filter(|edge| edge.exists())
        .map(|edge| (edge.weight, edge.length, edge.out_node))
        .collect();

    for (path_weight, step, mut source_node) in bubbles {
        if step == 0 {
            continue;
        }
        let mut last_node = start;
        for _ in 0..step {
            let base = backbone.nodes[source_node as usize].base;
            let new_node = dag.graph.add_node(base);
            dag.graph.add_bidirectional_edge(last_node, new_node, path_weight, 0);
            let next = backbone.nodes[source_node as usize]
                .insertion_edges
                .edges
                .first()
                .map(|edge| edge.out_node)
                .unwrap_or(source_node);
            source_node = next;
            last_node = new_node;
        }
        if last_node != start {
            dag.graph.add_bidirectional_edge(last_node, end, path_weight, 0);
        }
    }

    merge_dag_con(dag);
    let maximum_count = generate_best_sequence(dag);
    backbone.nodes[current_node as usize].num_insertions = 0;
    maximum_count
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds a backbone with two insertion bubbles at node 1: "GG" with
    /// weight 3 and "TT" with weight 1.
    fn backbone_with_bubbles() -> PoaGraph {
        let mut graph = PoaGraph::new();
        graph.add_unmatched_sequence(b"ACGT");
        graph.add_insertion_edge_weight(1, b"GG");
        graph.add_insertion_edge_weight(1, b"TT");
        // Raise the support of the first bubble ("GG") to three reads.
        // The anchor holds exactly one outgoing edge per bubble.
        let edges = &mut graph.nodes[1].insertion_edges.edges;
        edges[0].weight = 3;
        edges[1].weight = 1;
        graph.nodes[1].num_insertions = 4;
        graph
    }

    #[test]
    fn the_dag_picks_the_best_supported_insertion() {
        let mut backbone = backbone_with_bubbles();
        let mut dag = DagCon::default();

        let support = build_dag_con(&mut dag, &mut backbone, 1);

        assert_eq!(support, 3, "the three-read bubble wins");
        assert_eq!(dag.node_queue.iter().copied().collect::<Vec<u8>>(), b"GG");
        assert_eq!(backbone.nodes[1].num_insertions, 0, "the anchor is cleared");
    }

    #[test]
    fn the_dag_contains_start_and_end_nodes() {
        let mut backbone = backbone_with_bubbles();
        let mut dag = DagCon::default();
        build_dag_con(&mut dag, &mut backbone, 1);

        assert_eq!(dag.graph.nodes[dag.graph.start_node_identifier as usize].base, b'S');
        assert_eq!(dag.graph.nodes[dag.graph.end_node_identifier as usize].base, b'E');
        assert!(dag.graph.nodes.len() >= 4);
    }

    #[test]
    fn an_anchor_without_bubbles_produces_no_sequence() {
        let mut backbone = PoaGraph::new();
        backbone.add_unmatched_sequence(b"ACGT");
        let mut dag = DagCon::default();

        let support = build_dag_con(&mut dag, &mut backbone, 1);
        assert_eq!(support, 0);
        assert!(dag.node_queue.is_empty());
    }

    /// Golden values printed by the real `build_DAGCon` for four insertion
    /// configurations anchored at node 1 of an `ACGT` backbone.
    #[test]
    fn matches_the_reference_insertion_consensus() {
        struct Case {
            name: &'static str,
            bubbles: Vec<&'static [u8]>,
            weights: Vec<u64>,
            support: u64,
            sequence: &'static [u8],
            node_count: usize,
        }
        let cases = [
            Case {
                name: "two_bubbles",
                bubbles: vec![b"GG".as_slice(), b"TT".as_slice()],
                weights: vec![3, 1],
                support: 3,
                sequence: b"GG",
                node_count: 6,
            },
            Case {
                name: "one_bubble",
                bubbles: vec![b"GG".as_slice()],
                weights: vec![3],
                support: 3,
                sequence: b"GG",
                node_count: 4,
            },
            Case {
                name: "three_single",
                bubbles: vec![b"G".as_slice(), b"T".as_slice(), b"C".as_slice()],
                weights: vec![5, 2, 1],
                support: 5,
                sequence: b"G",
                node_count: 5,
            },
            Case {
                name: "shared_prefix",
                bubbles: vec![b"GG".as_slice(), b"GA".as_slice()],
                weights: vec![4, 1],
                support: 5,
                sequence: b"GG",
                node_count: 6,
            },
        ];

        for case in cases {
            let mut backbone = PoaGraph::new();
            backbone.add_unmatched_sequence(b"ACGT");
            for bubble in &case.bubbles {
                backbone.add_insertion_edge_weight(1, *bubble);
            }
            let edges = &mut backbone.nodes[1].insertion_edges.edges;
            for (index, weight) in case.weights.iter().enumerate() {
                edges[index].weight = *weight;
            }
            backbone.nodes[1].num_insertions = case.weights.iter().sum();

            let mut dag = DagCon::default();
            let support = build_dag_con(&mut dag, &mut backbone, 1);
            let sequence: Vec<u8> = dag.node_queue.iter().copied().collect();

            assert_eq!(support, case.support, "support of {}", case.name);
            assert_eq!(sequence, case.sequence, "sequence of {}", case.name);
            assert_eq!(dag.graph.nodes.len(), case.node_count, "nodes of {}", case.name);
            assert_eq!(backbone.nodes[1].num_insertions, 0);
        }
    }
}

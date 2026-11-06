//! Partial order alignment graph, the counterpart of `POA.h` and `POA.cpp`.
//!
//! The graph keeps three edge lists per node. Their names follow the
//! reference macros rather than the direction they intuitively carry:
//! `deletion_edges` holds the edges leaving a node, `insertion_edges` holds
//! both the edges entering a node and the insertion bubbles anchored at it,
//! and `mismatch_edges` holds the alternative bases leaving a node. Deleted
//! edges keep their slot with every field set to `DELETED`.

/// Sentinel written into every field of a removed edge, `(uint64_t)-1`.
pub const DELETED: u64 = u64::MAX;

/// Sentinel base of a removed node.
pub const DELETED_BASE: u8 = b'D';

/// Start node base, the one that precedes the first base of a read.
pub const START_BASE: u8 = b'S';

/// One graph edge, `Edge`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Edge {
    pub in_node: u64,
    pub out_node: u64,
    pub weight: u64,
    pub num_insertions: u64,
    pub length: u64,
    pub self_edge_id: u64,
    pub reverse_edge_id: u64,
}

impl Edge {
    /// A removed edge.
    pub const REMOVED: Self = Self {
        in_node: DELETED,
        out_node: DELETED,
        weight: DELETED,
        num_insertions: DELETED,
        length: DELETED,
        self_edge_id: DELETED,
        reverse_edge_id: DELETED,
    };

    /// `If_Edge_Exist`.
    pub fn exists(&self) -> bool {
        self.out_node != DELETED
    }
}

/// One edge list, `Edge_alloc`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct EdgeList {
    pub edges: Vec<Edge>,
    /// Number of removed edges still occupying a slot.
    pub delete_length: u64,
}

impl EdgeList {
    /// `append_Edge_alloc`.
    pub fn push(&mut self, in_node: u64, out_node: u64, weight: u64, length: u64) {
        let self_edge_id = self.edges.len() as u64;
        self.edges.push(Edge {
            in_node,
            out_node,
            weight,
            num_insertions: 0,
            length,
            self_edge_id,
            reverse_edge_id: 0,
        });
    }

    /// `Real_Length`: number of live edges.
    pub fn real_length(&self) -> u64 {
        self.edges.len() as u64 - self.delete_length
    }

    /// `clear_Edge_alloc`.
    pub fn clear(&mut self) {
        self.edges.clear();
        self.delete_length = 0;
    }

    /// Index of the first live edge from `in_node` to `out_node`.
    fn find(&self, in_node: u64, out_node: u64) -> Option<usize> {
        self.edges
            .iter()
            .position(|edge| edge.exists() && edge.in_node == in_node && edge.out_node == out_node)
    }
}

/// One graph node, `Node`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Node {
    pub identifier: u64,
    pub weight: u64,
    pub num_insertions: u64,
    pub base: u8,
    pub mismatch_edges: EdgeList,
    pub deletion_edges: EdgeList,
    pub insertion_edges: EdgeList,
}

impl Node {
    fn new(identifier: u64, base: u8) -> Self {
        Self {
            identifier,
            weight: 1,
            num_insertions: 0,
            base,
            mismatch_edges: EdgeList::default(),
            deletion_edges: EdgeList::default(),
            insertion_edges: EdgeList::default(),
        }
    }

    /// `If_Node_Exist`.
    pub fn exists(&self) -> bool {
        self.base != DELETED_BASE
    }
}

/// The whole graph, `Graph`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PoaGraph {
    pub nodes: Vec<Node>,
    /// Number of removed nodes.
    pub delete_length: u64,
    pub sequence_identifier: u64,
    pub start_node_identifier: u64,
    pub end_node_identifier: u64,
}

impl PoaGraph {
    /// `init_Graph` / `clear_Graph`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `add_Node_Graph`.
    pub fn add_node(&mut self, base: u8) -> u64 {
        let identifier = self.nodes.len() as u64;
        self.nodes.push(Node::new(identifier, base));
        identifier
    }

    /// Number of live nodes.
    pub fn live_node_count(&self) -> u64 {
        self.nodes.len() as u64 - self.delete_length
    }

    /// `get_bi_Edge`: finds the forward edge `in_node -> out_node`.
    fn find_forward_edge(&self, in_node: u64, out_node: u64) -> Option<usize> {
        if !self.nodes[in_node as usize].exists() || !self.nodes[out_node as usize].exists() {
            return None;
        }
        self.nodes[out_node as usize]
            .insertion_edges
            .edges
            .iter()
            .position(|edge| edge.exists() && edge.in_node == in_node)
    }

    /// `add_bi_direction_edge`: appends the edge to both endpoints and links
    /// the two copies with `reverse_edge_id`.
    pub fn add_bidirectional_edge(
        &mut self,
        in_node: u64,
        out_node: u64,
        weight: u64,
        length: u64,
    ) {
        self.nodes[in_node as usize]
            .deletion_edges
            .push(in_node, out_node, weight, length);
        self.nodes[out_node as usize]
            .insertion_edges
            .push(in_node, out_node, weight, length);

        let forward_index = self.nodes[in_node as usize].deletion_edges.edges.len() - 1;
        let backward_index = self.nodes[out_node as usize].insertion_edges.edges.len() - 1;

        self.nodes[in_node as usize].deletion_edges.edges[forward_index].reverse_edge_id =
            backward_index as u64;
        self.nodes[out_node as usize].insertion_edges.edges[backward_index].reverse_edge_id =
            forward_index as u64;
    }

    /// `get_bi_direction_edges`: locates both copies of an edge.
    ///
    /// Returns `(forward_slot, backward_slot)`: the index in `in_node`'s
    /// deletion-edge list and the index in `out_node`'s insertion-edge list.
    pub fn bidirectional_edge_slots(
        &self,
        in_node: u64,
        out_node: u64,
    ) -> Option<(usize, usize)> {
        let backward_slot = self.find_forward_edge(in_node, out_node)?;
        let backward = &self.nodes[out_node as usize].insertion_edges.edges[backward_slot];
        let forward_slot = backward.reverse_edge_id as usize;
        Some((forward_slot, backward_slot))
    }

    /// `get_Edge_Weight`: weight of the live edge `in_node -> out_node`.
    pub fn edge_weight(&self, in_node: u64, out_node: u64) -> Option<u64> {
        let (forward_slot, _) = self.bidirectional_edge_slots(in_node, out_node)?;
        Some(self.nodes[in_node as usize].deletion_edges.edges[forward_slot].weight)
    }

    /// Marks an edge as visited. The reference reuses `Edge::length` for this
    /// (`Visit(E)`), so the same field is used here.
    pub fn set_edge_visited(&mut self, in_node: u64, out_node: u64, visited: bool) -> bool {
        let Some((forward_slot, backward_slot)) = self.bidirectional_edge_slots(in_node, out_node)
        else {
            return false;
        };
        let value = u64::from(visited);
        self.nodes[in_node as usize].deletion_edges.edges[forward_slot].length = value;
        self.nodes[out_node as usize].insertion_edges.edges[backward_slot].length = value;
        true
    }

    /// `Visit(E)` of the edge copy at `slot` in `node`'s list.
    pub fn edge_visited_in_list(&self, node: u64, kind: EdgeListKind, slot: usize) -> bool {
        self.edge_list(node, kind).edges[slot].length != 0
    }

    /// Sets `Visit(E)` of one edge copy; used by the DAG merging pass which
    /// marks edges as visited through either endpoint.
    pub fn set_edge_visited_in_list(
        &mut self,
        node: u64,
        kind: EdgeListKind,
        slot: usize,
        visited: bool,
    ) {
        self.edge_list_mut(node, kind).edges[slot].length = u64::from(visited);
    }

    /// A copy of the edge at `slot` in `node`'s list.
    pub fn edge_at(&self, node: u64, kind: EdgeListKind, slot: usize) -> Edge {
        self.edge_list(node, kind).edges[slot]
    }

    /// Adds `delta` to the weight of one edge copy.
    pub fn add_edge_weight_in_list(
        &mut self,
        node: u64,
        kind: EdgeListKind,
        slot: usize,
        delta: u64,
    ) {
        self.edge_list_mut(node, kind).edges[slot].weight += delta;
    }

    /// Sets the weight of one edge copy.
    pub fn set_edge_weight_in_list(
        &mut self,
        node: u64,
        kind: EdgeListKind,
        slot: usize,
        weight: u64,
    ) {
        self.edge_list_mut(node, kind).edges[slot].weight = weight;
    }

    /// Number of live edges in one of a node's lists, `Real_Length`.
    pub fn live_edge_count(&self, node: u64, kind: EdgeListKind) -> u64 {
        self.edge_list(node, kind).real_length()
    }

    /// `clear_Graph`: drops every node and counter but keeps allocations.
    pub fn clear(&mut self) {
        self.nodes.clear();
        self.delete_length = 0;
        self.sequence_identifier = 0;
        self.start_node_identifier = 0;
        self.end_node_identifier = 0;
    }

    /// `getOutputEdges`: indices of the live edges leaving `node`.
    pub fn live_output_edge_slots(&self, node: u64) -> Vec<usize> {
        self.nodes[node as usize]
            .deletion_edges
            .edges
            .iter()
            .enumerate()
            .filter(|(_, edge)| edge.exists())
            .map(|(slot, _)| slot)
            .collect()
    }

    /// `getOutputNodes`: targets of the live edges leaving `node`.
    pub fn live_output_nodes(&self, node: u64) -> Vec<u64> {
        self.nodes[node as usize]
            .deletion_edges
            .edges
            .iter()
            .filter(|edge| edge.exists())
            .map(|edge| edge.out_node)
            .collect()
    }

    /// `getInputEdges`: indices of the live edges entering `node`.
    pub fn live_input_edge_slots(&self, node: u64) -> Vec<usize> {
        self.nodes[node as usize]
            .insertion_edges
            .edges
            .iter()
            .enumerate()
            .filter(|(_, edge)| edge.exists())
            .map(|(slot, _)| slot)
            .collect()
    }

    /// `getInputNodes`: sources of the live edges entering `node`.
    pub fn live_input_nodes(&self, node: u64) -> Vec<u64> {
        self.nodes[node as usize]
            .insertion_edges
            .edges
            .iter()
            .filter(|edge| edge.exists())
            .map(|edge| edge.in_node)
            .collect()
    }

    /// `remove_and_check_bi_direction_edge_from_nodes`.
    pub fn remove_bidirectional_edge(&mut self, in_node: u64, out_node: u64) -> bool {
        let Some((forward_slot, backward_slot)) = self.bidirectional_edge_slots(in_node, out_node)
        else {
            return false;
        };

        self.nodes[in_node as usize].deletion_edges.edges[forward_slot] = Edge::REMOVED;
        self.nodes[in_node as usize].deletion_edges.delete_length += 1;
        self.nodes[out_node as usize].insertion_edges.edges[backward_slot] = Edge::REMOVED;
        self.nodes[out_node as usize].insertion_edges.delete_length += 1;
        true
    }

    /// `delete_Node_DAGCon`: marks a node as deleted, removes every edge that
    /// touches it and empties its three edge lists.
    pub fn delete_node(&mut self, node: u64) -> bool {
        if node as usize >= self.nodes.len() || !self.nodes[node as usize].exists() {
            return false;
        }

        // The edge removal must happen while both endpoints still exist:
        // `get_bi_Edge` in the reference refuses deleted nodes, whereas the
        // deletion pass uses `get_bi_direction_edges` which does not. Doing the
        // removals first gives the same final graph.
        let (outgoing, incoming) = {
            let target = &self.nodes[node as usize];
            let outgoing: Vec<u64> = target
                .deletion_edges
                .edges
                .iter()
                .filter(|edge| edge.exists())
                .map(|edge| edge.out_node)
                .collect();
            let incoming: Vec<u64> = target
                .insertion_edges
                .edges
                .iter()
                .filter(|edge| edge.exists())
                .map(|edge| edge.in_node)
                .collect();
            (outgoing, incoming)
        };

        for neighbour in outgoing {
            self.remove_bidirectional_edge(node, neighbour);
        }
        for neighbour in incoming {
            self.remove_bidirectional_edge(neighbour, node);
        }

        self.delete_length += 1;
        let target = &mut self.nodes[node as usize];
        target.base = DELETED_BASE;
        target.num_insertions = DELETED;
        target.weight = DELETED;
        target.insertion_edges.clear();
        target.mismatch_edges.clear();
        target.deletion_edges.clear();
        true
    }

    /// `addUnmatchedSeqToGraph`: lays a whole read down as a fresh path.
    ///
    /// Returns the first and last node the path visits.
    pub fn add_unmatched_sequence(&mut self, sequence: &[u8]) -> (u64, u64) {
        let start = self.add_node(START_BASE);
        let mut last = start;

        for base in sequence {
            let node = self.add_node(*base);
            // Mismatch edges carry the sequence forward; a match edge would
            // have length zero, a mismatch edge has length one.
            self.nodes[last as usize]
                .mismatch_edges
                .push(last, node, 1, 0);
            last = node;
        }

        self.start_node_identifier = start;
        self.end_node_identifier = last;
        (start, last)
    }

    /// `add_mismatchEdge_weight`: records that `base` follows `in_node`.
    pub fn add_mismatch_edge_weight(&mut self, in_node: u64, base: u8, last_operation: i32) {
        let existing = self.nodes[in_node as usize]
            .mismatch_edges
            .edges
            .iter()
            .position(|edge| {
                edge.exists() && self.nodes[edge.out_node as usize].base == base
            });

        if let Some(index) = existing {
            let edge = &mut self.nodes[in_node as usize].mismatch_edges.edges[index];
            edge.weight += 1;
            if last_operation == 2 {
                edge.num_insertions += 1;
            }
            return;
        }

        let node = self.add_node(base);
        self.nodes[in_node as usize]
            .mismatch_edges
            .push(in_node, node, 1, 1);
        if last_operation == 2 {
            let edges = &mut self.nodes[in_node as usize].mismatch_edges.edges;
            let last = edges.len() - 1;
            edges[last].num_insertions += 1;
        }

        // The new node gets its own mismatch edge towards the next backbone
        // node, which is what `in_node + 1` means in the reference.
        self.nodes[node as usize]
            .mismatch_edges
            .push(node, in_node + 1, 1, 0);
    }

    /// `add_single_deletionEdge_weight`.
    fn add_single_deletion_edge_weight(
        &mut self,
        align_node: u64,
        next_node: u64,
        edge_length: u64,
    ) {
        let existing = self.nodes[align_node as usize]
            .deletion_edges
            .find(align_node, next_node);
        match existing {
            Some(index) => self.nodes[align_node as usize].deletion_edges.edges[index].weight += 1,
            None => self.nodes[align_node as usize]
                .deletion_edges
                .push(align_node, next_node, 1, edge_length),
        }
    }

    /// `add_deletionEdge_weight`: a deletion removes `deletion_length` bases of
    /// the backbone, so each of them gets a zero length edge to the next node.
    pub fn add_deletion_edge_weight(&mut self, align_node: u64, deletion_length: u64) {
        for offset in 0..deletion_length {
            self.add_single_deletion_edge_weight(align_node + offset, align_node + offset + 1, 0);
        }
    }

    /// Walks an existing single base insertion path, `getEdge`.
    fn find_edge(
        &self,
        in_node: u64,
        list_kind: EdgeListKind,
        edge_length: u64,
        base: u8,
    ) -> Option<usize> {
        let list = self.edge_list(in_node, list_kind);
        list.edges.iter().position(|edge| {
            edge.exists()
                && edge.length == edge_length
                && self.nodes[edge.out_node as usize].base == base
        })
    }

    fn edge_list(&self, node: u64, kind: EdgeListKind) -> &EdgeList {
        match kind {
            EdgeListKind::Insertion => &self.nodes[node as usize].insertion_edges,
            EdgeListKind::Mismatch => &self.nodes[node as usize].mismatch_edges,
            EdgeListKind::Deletion => &self.nodes[node as usize].deletion_edges,
        }
    }

    fn edge_list_mut(&mut self, node: u64, kind: EdgeListKind) -> &mut EdgeList {
        match kind {
            EdgeListKind::Insertion => &mut self.nodes[node as usize].insertion_edges,
            EdgeListKind::Mismatch => &mut self.nodes[node as usize].mismatch_edges,
            EdgeListKind::Deletion => &mut self.nodes[node as usize].deletion_edges,
        }
    }

    /// `get_insertion_Edges_new`: matches an insertion path by walking it.
    fn find_insertion_path(&self, anchor: u64, bases: &[u8]) -> Option<usize> {
        let edges = &self.nodes[anchor as usize].insertion_edges.edges;
        for (path_index, edge) in edges.iter().enumerate() {
            if !edge.exists() || edge.length != bases.len() as u64 {
                continue;
            }
            let mut node = edge.out_node;
            let mut matched = 0usize;
            for base in bases {
                if self.nodes[node as usize].base != *base {
                    break;
                }
                matched += 1;
                let next = &self.nodes[node as usize].insertion_edges.edges;
                if next.is_empty() {
                    break;
                }
                node = next[0].out_node;
            }
            if matched == bases.len() {
                return Some(path_index);
            }
        }
        None
    }

    /// `create_insertion_Edges`: builds a new bubble for an unseen insertion.
    fn create_insertion_edges(&mut self, anchor: u64, bases: &[u8]) -> bool {
        if bases.is_empty() {
            return false;
        }

        let mut previous = anchor;
        for (offset, base) in bases.iter().enumerate() {
            let node = self.add_node(*base);
            let remaining = bases.len() as u64 - offset as u64;
            self.edge_list_mut(previous, EdgeListKind::Insertion)
                .push(previous, node, 1, remaining);
            previous = node;
        }
        // The path returns to the anchor.
        self.edge_list_mut(previous, EdgeListKind::Insertion)
            .push(previous, anchor, 1, 0);
        true
    }

    /// `add_insertionEdge_weight`.
    pub fn add_insertion_edge_weight(&mut self, anchor: u64, bases: &[u8]) {
        if bases.is_empty() {
            return;
        }

        if bases.len() == 1 {
            match self.find_edge(anchor, EdgeListKind::Insertion, 1, bases[0]) {
                Some(index) => {
                    self.edge_list_mut(anchor, EdgeListKind::Insertion).edges[index].weight += 1
                }
                None => {
                    let node = self.add_node(bases[0]);
                    self.edge_list_mut(anchor, EdgeListKind::Insertion)
                        .push(anchor, node, 1, 1);
                    self.edge_list_mut(node, EdgeListKind::Insertion)
                        .push(node, anchor, 1, 0);
                }
            }
            return;
        }

        match self.find_insertion_path(anchor, bases) {
            Some(index) => {
                self.edge_list_mut(anchor, EdgeListKind::Insertion).edges[index].weight += 1
            }
            None => {
                self.create_insertion_edges(anchor, bases);
            }
        }
    }

    /// `addmatchedSeqToGraph`: folds a read whose cigar is known into the graph.
    ///
    /// `cells` are the decoded cigar cells `(operation, length)`; `target` is
    /// the read being folded in. Operations follow the reference numbering:
    /// `0` match, `1` mismatch, `2` the folded read has extra bases, `3` the
    /// backbone has extra bases.
    pub fn add_matched_sequence(
        &mut self,
        mut current_node: u64,
        target: &[u8],
        cells: &[(u32, u32)],
    ) {
        let mut target_index = 0usize;
        let mut last_operation: i32 = -1;

        for (operation, operation_length) in cells {
            match operation {
                0 | 1 => {
                    for _ in 0..*operation_length {
                        self.add_mismatch_edge_weight(
                            current_node,
                            target[target_index],
                            last_operation,
                        );
                        target_index += 1;
                        current_node += 1;
                    }
                }
                2 => {
                    self.add_insertion_edge_weight(
                        current_node,
                        &target[target_index..target_index + *operation_length as usize],
                    );
                    self.nodes[current_node as usize].num_insertions += 1;
                    target_index += *operation_length as usize;
                }
                3 => {
                    self.add_deletion_edge_weight(current_node, *operation_length as u64);
                    current_node += *operation_length as u64;
                }
                other => panic!("unknown cigar operation {other}"),
            }
            last_operation = *operation as i32;
        }
    }
}

/// Which edge list of a node to read or write.
///
/// `Mismatch` and `Deletion` are not queried yet: only the insertion path
/// lookup needs a list selector at the moment, the other two are addressed
/// directly by the edge-weight functions.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeListKind {
    Insertion,
    Mismatch,
    Deletion,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_unmatched_sequence_becomes_a_path_of_nodes() {
        let mut graph = PoaGraph::new();
        let (start, end) = graph.add_unmatched_sequence(b"ACG");

        assert_eq!(graph.nodes.len(), 4, "one start node plus three bases");
        assert_eq!(graph.nodes[start as usize].base, START_BASE);
        assert_eq!(graph.nodes[end as usize].base, b'G');
        assert_eq!(graph.graph_path(), b"ACG");
    }

    /// `render` prints the graph in the same shape as the C probe, so the
    /// expected text below is copied straight from the reference run.
    fn render(graph: &PoaGraph) -> String {
        fn edges(list: &EdgeList) -> String {
            let rendered: Vec<String> = list
                .edges
                .iter()
                .map(|edge| {
                    format!(
                        "({}->{},w{},l{},ins{})",
                        edge.in_node, edge.out_node, edge.weight, edge.length, edge.num_insertions
                    )
                })
                .collect();
            format!("n={}: {}", list.edges.len(), rendered.join(" "))
        }

        let mut text = format!(
            "nodes={} deleted={} start={} end={}\n",
            graph.nodes.len(),
            graph.delete_length,
            graph.start_node_identifier,
            graph.end_node_identifier
        );
        for node in &graph.nodes {
            text.push_str(&format!(
                "  node {} base={} weight={} ins={}\n",
                node.identifier, node.base as char, node.weight, node.num_insertions
            ));
            text.push_str(&format!("    mismatch {}\n", edges(&node.mismatch_edges)));
            text.push_str(&format!("    deletion {}\n", edges(&node.deletion_edges)));
            text.push_str(&format!("    insertion {}\n", edges(&node.insertion_edges)));
        }
        text
    }

    /// Exercises the same five steps as the C probe and compares the whole
    /// graph, node by node and edge by edge.
    #[test]
    fn matches_the_reference_graph_layout() {
        let mut graph = PoaGraph::new();
        graph.add_unmatched_sequence(b"ACGT");
        assert_eq!(
            render(&graph),
            "nodes=5 deleted=0 start=0 end=4\n\
             \u{20}\u{20}node 0 base=S weight=1 ins=0\n\
             \u{20}\u{20}\u{20}\u{20}mismatch n=1: (0->1,w1,l0,ins0)\n\
             \u{20}\u{20}\u{20}\u{20}deletion n=0: \n\
             \u{20}\u{20}\u{20}\u{20}insertion n=0: \n\
             \u{20}\u{20}node 1 base=A weight=1 ins=0\n\
             \u{20}\u{20}\u{20}\u{20}mismatch n=1: (1->2,w1,l0,ins0)\n\
             \u{20}\u{20}\u{20}\u{20}deletion n=0: \n\
             \u{20}\u{20}\u{20}\u{20}insertion n=0: \n\
             \u{20}\u{20}node 2 base=C weight=1 ins=0\n\
             \u{20}\u{20}\u{20}\u{20}mismatch n=1: (2->3,w1,l0,ins0)\n\
             \u{20}\u{20}\u{20}\u{20}deletion n=0: \n\
             \u{20}\u{20}\u{20}\u{20}insertion n=0: \n\
             \u{20}\u{20}node 3 base=G weight=1 ins=0\n\
             \u{20}\u{20}\u{20}\u{20}mismatch n=1: (3->4,w1,l0,ins0)\n\
             \u{20}\u{20}\u{20}\u{20}deletion n=0: \n\
             \u{20}\u{20}\u{20}\u{20}insertion n=0: \n\
             \u{20}\u{20}node 4 base=T weight=1 ins=0\n\
             \u{20}\u{20}\u{20}\u{20}mismatch n=0: \n\
             \u{20}\u{20}\u{20}\u{20}deletion n=0: \n\
             \u{20}\u{20}\u{20}\u{20}insertion n=0: \n"
        );

        // Folding the same read in again only raises the edge weights.
        let start = graph.start_node_identifier;
        graph.add_matched_sequence(start, b"ACGT", &[(0, 4)]);
        assert_eq!(graph.nodes.len(), 5, "no node is created");
        assert_eq!(graph.nodes[0].mismatch_edges.edges[0].weight, 2);
        assert_eq!(graph.nodes[3].mismatch_edges.edges[0].weight, 2);

        // A mismatch at the second base creates one alternative node, and that
        // node points back at the backbone node after the mismatch.
        graph.add_mismatch_edge_weight(start, b'A', -1);
        graph.add_mismatch_edge_weight(start + 1, b'T', 0);
        graph.add_mismatch_edge_weight(start + 2, b'G', 1);
        graph.add_mismatch_edge_weight(start + 3, b'T', 0);

        assert_eq!(graph.nodes.len(), 6);
        let branches = &graph.nodes[1].mismatch_edges.edges;
        assert_eq!(branches.len(), 2);
        assert_eq!(graph.nodes[branches[1].out_node as usize].base, b'T');
        assert_eq!(branches[1].weight, 1);
        assert_eq!(branches[1].length, 1);
        assert_eq!(graph.nodes[5].mismatch_edges.edges[0].in_node, 5);
        assert_eq!(graph.nodes[5].mismatch_edges.edges[0].out_node, 2);

        // An insertion hangs a bubble off the anchor and returns to it.
        graph.add_insertion_edge_weight(start + 1, b"G");
        assert_eq!(graph.nodes.len(), 7);
        let bubble = &graph.nodes[1].insertion_edges.edges[0];
        assert_eq!(bubble.out_node, 6);
        assert_eq!(bubble.length, 1);
        assert_eq!(graph.nodes[6].base, b'G');
        let return_edge = &graph.nodes[6].insertion_edges.edges[0];
        assert_eq!(return_edge.out_node, 1);
        assert_eq!(return_edge.length, 0);

        // A deletion links the node before the removed base to the node after.
        graph.add_deletion_edge_weight(start + 1, 1);
        let deletion = &graph.nodes[1].deletion_edges.edges[0];
        assert_eq!(deletion.out_node, 2);
        assert_eq!(deletion.length, 0);
        assert_eq!(deletion.weight, 1);
        assert_eq!(graph.nodes.len(), 7, "a deletion adds no node");
    }

    #[test]
    fn a_removed_edge_keeps_its_slot() {
        let mut graph = PoaGraph::new();
        graph.add_unmatched_sequence(b"AC");
        graph.add_bidirectional_edge(0, 1, 1, 0);

        assert!(graph.remove_bidirectional_edge(0, 1));
        assert!(!graph.nodes[0].deletion_edges.edges[0].exists());
        assert_eq!(graph.nodes[0].deletion_edges.edges.len(), 1);
        assert_eq!(graph.nodes[0].deletion_edges.delete_length, 1);
        assert_eq!(graph.nodes[0].deletion_edges.real_length(), 0);
        assert!(!graph.remove_bidirectional_edge(0, 1), "already removed");
    }

    impl PoaGraph {
        /// Walks the backbone path, for tests only.
        fn graph_path(&self) -> Vec<u8> {
            let mut path = Vec::new();
            let mut node = self.start_node_identifier;
            while let Some(edge) = self.nodes[node as usize]
                .mismatch_edges
                .edges
                .iter()
                .find(|edge| edge.exists())
            {
                node = edge.out_node;
                path.push(self.nodes[node as usize].base);
            }
            path
        }
    }
}

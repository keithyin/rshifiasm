//! Window consensus generation, the counterpart of [`cns_gfa`] and its helpers
//! in `ecovlp.cpp`.
//!
//! [`ConsensusScratch::wcns_gen`] is the reference's `wcns_gen`: it walks the
//! query in windows, votes base by base on the covered overlaps, and fills the
//! uncovered stretches with a consensus built from a small partial order graph
//! (`cns_gfa`). The result is appended to one auxiliary window as a packed,
//! base-aware correction cigar.
//!
//! The reference keeps a lot of state in flat vectors with hand-maintained
//! cursors. Here the same state lives in named structs, but every index and
//! every comparison keeps the reference's exact arithmetic, including the
//! places where it uses unsigned subtraction and 32-bit truncation.
//!
//! # Differences from the reference
//!
//! * `get_sub_seq` may return bytes past the end of a read when the requested
//!   sub-region starts beyond it: the reference clamps the recovered length and
//!   then hands back its scratch buffer anyway, so the caller reads whatever
//!   was left there. This port fills such bytes with `A` (`0` after encoding),
//!   which only differs on inputs that would make the reference read
//!   uninitialised memory.

use std::collections::VecDeque;

use crate::correction::cal_exz::cal_exz_global;
use crate::correction::hc_align::BitExtz;
use crate::correction::rphase::{UlOverlap, SEC_MODE};
use crate::correction::window::WindowCigar;
use crate::correction::window_chain::is_unaligned_window;
use crate::container::ksort::{radix_sort_by_key_32, radix_sort_u64};
use crate::overlap::region::OverlapRegion;
use crate::reads::read_store::ReadStore;
use crate::sequence::encoding::{BASE_CHARACTERS, SEQUENCE_NT6_TABLE};

/// `CNS_DEL_E`.
const DELETED_EDGE: u32 = 0x7fff_ffff;
/// `CNS_DEL_V`.
const DELETED_NODE: u32 = 0x1fff_ffff;
/// `simp_vote_len`.
const SIMPLE_VOTE_LENGTH: u64 = 6;
/// `(uint32_t)-1`.
const UNSET: u32 = u32::MAX;

/// `cns_arc`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct CnsArc {
    /// Target node, `v`.
    target: u32,
    /// Visited flag, `f`.
    visited: bool,
    /// Arc weight, `sc`.
    score: u32,
}

/// `cns_arc_v`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct CnsEdges {
    /// Out-edges come first, `nou` of them, then the in-edges.
    arcs: Vec<CnsArc>,
    /// Number of leading out-edges, `nou`.
    out_count: u32,
}

/// `cns_t`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct CnsNode {
    /// Base, two bits, `c`.
    base: u8,
    /// Visited flag used by [`CnsGfa::generate_sequence`], `f`.
    visited: bool,
    /// Weight, `sc`; also the in-degree during the sequence pass.
    score: u32,
    /// The arcs of this node.
    edges: CnsEdges,
}

/// `cns_gfa`.
#[derive(Debug, Clone, Default)]
struct CnsGfa {
    nodes: Vec<CnsNode>,
    /// Start node, `si`.
    start: u32,
    /// End node, `ei`.
    end: u32,
    /// Node id of the first backbone base, `off`.
    offset: u32,
    /// Backbone length, `bn`.
    backbone_length: u32,
    /// First backbone node, `bb0`.
    backbone_start: u32,
    /// One past the last backbone node, `bb1`.
    backbone_end: u32,
    /// Window length used to split long consensus gaps, `cns_g_wl`.
    consensus_window_length: u32,
    /// `q`.
    queue: VecDeque<u32>,
}

impl CnsGfa {
    fn is_deleted_edge(&self, node: u32, index: usize) -> bool {
        self.nodes[node as usize].edges.arcs[index].target == DELETED_EDGE
    }

    fn is_deleted_node(&self, node: u32) -> bool {
        self.nodes[node as usize].score == DELETED_NODE
    }

    /// `insert_cns_arc`: appends one arc; an out-edge is moved in front of the
    /// in-edge block.
    fn insert_arc_with_flag(&mut self, source: u32, target: u32, is_out: bool, weight: u32, flag: bool) {
        let edges = &mut self.nodes[source as usize].edges;
        edges.arcs.push(CnsArc {
            target,
            visited: flag,
            score: weight,
        });
        if is_out {
            edges.out_count += 1;
            let last = edges.arcs.len() - 1;
            let position = edges.out_count as usize - 1;
            if position < last {
                edges.arcs.swap(position, last);
            }
        }
    }

    fn insert_arc(&mut self, source: u32, target: u32, is_out: bool, weight: u32) {
        self.insert_arc_with_flag(source, target, is_out, weight, false);
    }

    /// `insert_cns_node`: appends a zeroed node and returns its id.
    fn insert_node(&mut self) -> u32 {
        self.nodes.push(CnsNode::default());
        (self.nodes.len() - 1) as u32
    }

    /// `add_cns_arc`: adds `plus` to the weight of an existing arc.
    fn add_arc(&mut self, source: u32, target: u32, is_out: bool, plus: u32) -> bool {
        let edges = &mut self.nodes[source as usize].edges;
        let (begin, end) = if is_out {
            (0, edges.out_count as usize)
        } else {
            (edges.out_count as usize, edges.arcs.len())
        };
        for k in begin..end {
            if edges.arcs[k].target == target {
                edges.arcs[k].score += plus;
                return true;
            }
        }
        false
    }

    /// `get_cns_arc_bp`: the index of the out-arc of `source` that leads to a
    /// node carrying `base`, or `UNSET`.
    fn get_arc_base(&self, source: u32, base: u8, avoid_backbone: bool) -> u32 {
        let edges = &self.nodes[source as usize].edges;
        for k in 0..edges.out_count as usize {
            let target = edges.arcs[k].target;
            if target == 0 || target == 1 {
                continue;
            }
            if avoid_backbone && target >= self.backbone_start && target < self.backbone_end {
                continue;
            }
            // A deleted arc points at `DELETED_EDGE`; the reference would read
            // out of bounds there, which cannot happen before `refine`.
            if target as usize >= self.nodes.len() {
                continue;
            }
            if self.nodes[target as usize].base == base {
                return k as u32;
            }
        }
        UNSET
    }

    /// `add_cns_arc_bp`: follows or creates the arc carrying `base`.
    fn add_arc_base(&mut self, source: u32, base: u8, plus0: u32, avoid_backbone: bool) -> u32 {
        let found = self.get_arc_base(source, base, avoid_backbone);
        if found != UNSET {
            let target = self.nodes[source as usize].edges.arcs[found as usize].target;
            self.nodes[target as usize].score += 1;
            self.nodes[source as usize].edges.arcs[found as usize].score += plus0;
            assert!(self.add_arc(target, source, false, plus0));
            return target;
        }
        let target = self.insert_node();
        self.nodes[target as usize].score += 1;
        self.nodes[target as usize].base = base & 3;
        self.insert_arc(source, target, true, plus0);
        self.insert_arc(target, source, false, plus0);
        target
    }

    /// `init_cns_g`: resets the graph to the backbone `sequence`.
    fn init(&mut self, sequence: &[u8]) {
        self.nodes.clear();
        self.queue.clear();
        self.start = 0;
        self.end = 1;
        self.offset = 2;
        self.nodes.push(CnsNode::default());
        self.nodes.push(CnsNode::default());
        self.backbone_start = 2;
        let length = sequence.len();
        for (k, &base) in sequence.iter().enumerate() {
            self.nodes.push(CnsNode {
                base: SEQUENCE_NT6_TABLE[base as usize] & 3,
                visited: false,
                score: 1,
                edges: CnsEdges::default(),
            });
            if k + 1 < length {
                self.insert_arc(k as u32 + 2, k as u32 + 3, true, 1);
            }
            if k > 0 {
                self.insert_arc(k as u32 + 2, k as u32 + 1, false, 1);
            }
        }
        if length > 0 {
            let last = length as u32 - 1 + 2;
            self.insert_arc(self.start, 2, true, 1);
            self.insert_arc(2, self.start, false, 1);
            self.insert_arc(last, self.end, true, 1);
            self.insert_arc(self.end, last, false, 1);
        } else {
            self.insert_arc(self.start, self.end, true, 1);
            self.insert_arc(self.end, self.start, false, 1);
        }
        self.backbone_length = self.nodes.len() as u32;
        self.backbone_end = self.nodes.len() as u32;
    }

    /// `del_cns_g_nn`: removes a node and every arc that touches it.
    fn delete_node(&mut self, node: u32) {
        let out_count = self.nodes[node as usize].edges.out_count as usize;
        let edge_count = self.nodes[node as usize].edges.arcs.len();
        for k in 0..out_count {
            if self.is_deleted_edge(node, k) {
                continue;
            }
            let target = self.nodes[node as usize].edges.arcs[k].target;
            self.nodes[node as usize].edges.arcs[k].target = DELETED_EDGE;
            let target_out = self.nodes[target as usize].edges.out_count as usize;
            let target_count = self.nodes[target as usize].edges.arcs.len();
            for i in target_out..target_count {
                if self.nodes[target as usize].edges.arcs[i].target == node
                    && !self.is_deleted_edge(target, i)
                {
                    self.nodes[target as usize].edges.arcs[i].target = DELETED_EDGE;
                    break;
                }
            }
        }
        for k in out_count..edge_count {
            if self.is_deleted_edge(node, k) {
                continue;
            }
            let source = self.nodes[node as usize].edges.arcs[k].target;
            self.nodes[node as usize].edges.arcs[k].target = DELETED_EDGE;
            let source_out = self.nodes[source as usize].edges.out_count as usize;
            let source_count = self.nodes[source as usize].edges.arcs.len();
            for i in source_out..source_count {
                if self.nodes[source as usize].edges.arcs[i].target == node
                    && !self.is_deleted_edge(source, i)
                {
                    self.nodes[source as usize].edges.arcs[i].target = DELETED_EDGE;
                    break;
                }
            }
        }
        self.nodes[node as usize].edges.arcs.clear();
        self.nodes[node as usize].edges.out_count = 0;
        self.nodes[node as usize].base = 0;
        self.nodes[node as usize].visited = false;
        self.nodes[node as usize].score = DELETED_NODE;
    }

    /// `gen_mm_cns_arc`: adds `score` to the `source -> target` arc, creating
    /// it when absent, and marks both directions with `flag`.
    fn merge_arc(&mut self, source: u32, target: u32, score: u32, flag: bool) {
        let out_count = self.nodes[source as usize].edges.out_count as usize;
        for k in 0..out_count {
            if self.nodes[source as usize].edges.arcs[k].target != target
                || self.is_deleted_edge(source, k)
            {
                continue;
            }
            self.nodes[source as usize].edges.arcs[k].visited = flag;
            self.nodes[source as usize].edges.arcs[k].score += score;
            let target_out = self.nodes[target as usize].edges.out_count as usize;
            let target_count = self.nodes[target as usize].edges.arcs.len();
            let mut found = false;
            for i in target_out..target_count {
                if self.nodes[target as usize].edges.arcs[i].target == source
                    && !self.is_deleted_edge(target, i)
                {
                    self.nodes[target as usize].edges.arcs[i].visited = flag;
                    self.nodes[target as usize].edges.arcs[i].score += score;
                    found = true;
                    break;
                }
            }
            assert!(found);
            return;
        }
        self.insert_arc_with_flag(source, target, true, score, flag);
        self.insert_arc_with_flag(target, source, false, score, flag);
    }

    /// `merge_cns_g_in`: merges the in-neighbours of `v` that share a base.
    fn merge_in(&mut self, v: u32, scratch: &mut Vec<u32>) {
        scratch.clear();
        scratch.push(v);
        while let Some(current) = scratch.pop() {
            if self.is_deleted_node(current) {
                continue;
            }
            for base in 0..4u8 {
                let mut count = 0u32;
                let mut weight = 0u32;
                let mut merge_node = 0u32;
                let mut merge_arc = 0u32;
                let mut merge_reverse_arc = 0u32;
                let mut vk = self.nodes[current as usize].edges.out_count as usize;
                while vk < self.nodes[current as usize].edges.arcs.len() {
                    if self.is_deleted_edge(current, vk) {
                        vk += 1;
                        continue;
                    }
                    let w = self.nodes[current as usize].edges.arcs[vk].target;
                    if self.nodes[w as usize].base != base || w == self.start || w == self.end {
                        vk += 1;
                        continue;
                    }
                    let mut reverse_count = 0usize;
                    let mut reverse_arc = UNSET;
                    let w_out = self.nodes[w as usize].edges.out_count as usize;
                    let mut wk = 0usize;
                    while wk < w_out {
                        if !self.is_deleted_edge(w, wk) {
                            reverse_count += 1;
                            reverse_arc = wk as u32;
                            if reverse_count > 1 {
                                break;
                            }
                        }
                        wk += 1;
                    }
                    if reverse_count != 1 {
                        vk += 1;
                        continue;
                    }
                    assert_eq!(self.nodes[w as usize].edges.arcs[reverse_arc as usize].target, current);
                    if count == 0 {
                        merge_node = w;
                        merge_arc = vk as u32;
                        merge_reverse_arc = reverse_arc;
                        weight = self.nodes[current as usize].edges.arcs[vk].score;
                    } else {
                        weight += self.nodes[w as usize].edges.arcs[reverse_arc as usize].score;
                    }
                    if count > 0 {
                        let w_out = self.nodes[w as usize].edges.out_count as usize;
                        let mut wk = w_out;
                        while wk < self.nodes[w as usize].edges.arcs.len() {
                            if self.is_deleted_edge(w, wk) {
                                wk += 1;
                                continue;
                            }
                            let source = self.nodes[w as usize].edges.arcs[wk].target;
                            let arc_weight = self.nodes[w as usize].edges.arcs[wk].score;
                            let flag = self.nodes[w as usize].edges.arcs[wk].visited;
                            self.merge_arc(source, merge_node, arc_weight, flag);
                            wk += 1;
                        }
                        self.delete_node(w);
                    }
                    count += 1;
                    vk += 1;
                }
                if count > 0 {
                    self.nodes[current as usize].edges.arcs[merge_arc as usize].score = weight;
                    self.nodes[merge_node as usize].edges.arcs[merge_reverse_arc as usize].score = weight;
                    scratch.push(merge_node);
                }
            }
        }
    }

    /// `merge_cns_g_out`: merges the out-neighbours of `v` that share a base.
    fn merge_out(&mut self, v: u32, scratch: &mut Vec<u32>) {
        scratch.clear();
        scratch.push(v);
        while let Some(current) = scratch.pop() {
            if self.is_deleted_node(current) {
                continue;
            }
            for base in 0..4u8 {
                let mut count = 0u32;
                let mut weight = 0u32;
                let mut merge_node = 0u32;
                let mut merge_arc = 0u32;
                let mut merge_reverse_arc = 0u32;
                let mut vk = 0usize;
                while vk < self.nodes[current as usize].edges.out_count as usize {
                    if self.is_deleted_edge(current, vk) {
                        vk += 1;
                        continue;
                    }
                    let w = self.nodes[current as usize].edges.arcs[vk].target;
                    if self.nodes[w as usize].base != base || w == self.start || w == self.end {
                        vk += 1;
                        continue;
                    }
                    let in_count = self.nodes[w as usize].edges.arcs.len();
                    let mut in_edges = 0usize;
                    let mut reverse_arc = UNSET;
                    let mut wk = self.nodes[w as usize].edges.out_count as usize;
                    while wk < in_count {
                        if !self.is_deleted_edge(w, wk) {
                            in_edges += 1;
                            reverse_arc = wk as u32;
                            if in_edges > 1 {
                                break;
                            }
                        }
                        wk += 1;
                    }
                    if in_edges != 1 {
                        vk += 1;
                        continue;
                    }
                    assert_eq!(self.nodes[w as usize].edges.arcs[reverse_arc as usize].target, current);
                    if count == 0 {
                        merge_node = w;
                        merge_arc = vk as u32;
                        merge_reverse_arc = reverse_arc;
                        weight = self.nodes[current as usize].edges.arcs[vk].score;
                    } else {
                        weight += self.nodes[w as usize].edges.arcs[reverse_arc as usize].score;
                    }
                    if count > 0 {
                        let w_out = self.nodes[w as usize].edges.out_count as usize;
                        let mut wk = 0usize;
                        while wk < w_out {
                            if self.is_deleted_edge(w, wk) {
                                wk += 1;
                                continue;
                            }
                            let target = self.nodes[w as usize].edges.arcs[wk].target;
                            let arc_weight = self.nodes[w as usize].edges.arcs[wk].score;
                            let flag = self.nodes[w as usize].edges.arcs[wk].visited;
                            self.merge_arc(merge_node, target, arc_weight, flag);
                            wk += 1;
                        }
                        self.delete_node(w);
                    }
                    count += 1;
                    vk += 1;
                }
                if count > 0 {
                    self.nodes[current as usize].edges.arcs[merge_arc as usize].score = weight;
                    self.nodes[merge_node as usize].edges.arcs[merge_reverse_arc as usize].score = weight;
                    scratch.push(merge_node);
                }
            }
        }
    }

    /// `refine_cns_g`: repeatedly merges the neighbours of every reachable node
    /// until the graph is a simple path per node.
    fn refine(&mut self, scratch: &mut Vec<u32>) {
        self.queue.clear();
        self.queue.push_back(self.start);
        while let Some(v) = self.queue.pop_front() {
            if self.is_deleted_node(v) {
                continue;
            }
            self.merge_in(v, scratch);
            self.merge_out(v, scratch);
            let out_count = self.nodes[v as usize].edges.out_count as usize;
            for vk in 0..out_count {
                if self.is_deleted_edge(v, vk) || !self.nodes[v as usize].edges.arcs[vk].visited {
                    continue;
                }
                self.nodes[v as usize].edges.arcs[vk].visited = true;
                let w = self.nodes[v as usize].edges.arcs[vk].target;
                let mut wk = self.nodes[w as usize].edges.out_count as usize;
                while wk < self.nodes[w as usize].edges.arcs.len() {
                    if self.nodes[w as usize].edges.arcs[wk].target == v
                        && !self.is_deleted_edge(w, wk)
                    {
                        self.nodes[w as usize].edges.arcs[wk].visited = true;
                        break;
                    }
                    wk += 1;
                }
            }
            self.nodes[v as usize].visited = true;
            let out_count = self.nodes[v as usize].edges.out_count as usize;
            for vk in 0..out_count {
                if self.is_deleted_edge(v, vk) {
                    continue;
                }
                let w = self.nodes[v as usize].edges.arcs[vk].target;
                let mut wk = self.nodes[w as usize].edges.out_count as usize;
                let w_count = self.nodes[w as usize].edges.arcs.len();
                while wk < w_count {
                    if self.is_deleted_edge(w, wk) {
                        wk += 1;
                        continue;
                    }
                    if self.nodes[w as usize].edges.arcs[wk].visited {
                        wk += 1;
                        continue;
                    }
                    let source = self.nodes[w as usize].edges.arcs[wk].target;
                    if self.nodes[source as usize].visited {
                        wk += 1;
                        continue;
                    }
                    break;
                }
                if wk >= w_count {
                    self.queue.push_back(w);
                }
            }
        }
    }

    /// `gseq_cns_g`: picks the heaviest path through the refined graph and
    /// leaves it in `path` (start and end nodes excluded).
    ///
    /// The reference's `be = bl + cns->off` is [`CnsGfa::backbone_end`], which
    /// [`CnsGfa::init`] always sets to the same value.
    fn generate_sequence(&mut self, path: &mut Vec<u32>) {
        path.clear();
        path.resize(self.nodes.len(), 0);
        for v in 0..self.nodes.len() {
            if self.is_deleted_node(v as u32) {
                continue;
            }
            self.nodes[v].score = 0;
            self.nodes[v].visited = false;
            let out_count = self.nodes[v].edges.out_count as usize;
            for vk in out_count..self.nodes[v].edges.arcs.len() {
                if self.is_deleted_edge(v as u32, vk) {
                    continue;
                }
                self.nodes[v].score += 1;
            }
        }
        self.queue.clear();
        self.queue.push_back(self.start);
        assert_eq!(self.nodes[self.start as usize].score, 0);
        assert!(!self.is_deleted_node(self.start));
        while let Some(v) = self.queue.pop_front() {
            if self.is_deleted_node(v) {
                continue;
            }
            assert_eq!(self.nodes[v as usize].score, 0);
            let mut best = UNSET;
            let mut best_edge = 0u32;
            let mut best_path = 0u32;
            let mut best_backbone = 0u32;
            let out_count = self.nodes[v as usize].edges.out_count as usize;
            for vk in out_count..self.nodes[v as usize].edges.arcs.len() {
                if self.is_deleted_edge(v, vk) {
                    continue;
                }
                let w = self.nodes[v as usize].edges.arcs[vk].target;
                assert!(self.nodes[w as usize].visited);
                let is_backbone =
                    u32::from(w >= self.backbone_start && w < self.backbone_end);
                let edge_score = self.nodes[v as usize].edges.arcs[vk].score;
                let path_score = path[w as usize];
                if best == UNSET
                    || edge_score > best_edge
                    || (edge_score == best_edge && path_score > best_path)
                    || (edge_score == best_edge
                        && path_score == best_path
                        && is_backbone == 1
                        && best_backbone == 0)
                {
                    best = vk as u32;
                    best_edge = edge_score;
                    best_path = path_score;
                    best_backbone = is_backbone;
                }
            }
            path[v as usize] = best_edge + best_path;
            self.nodes[v as usize].visited = true;
            if best != UNSET {
                self.nodes[v as usize].score = self.nodes[v as usize].edges.arcs[best as usize].target;
            } else {
                self.nodes[v as usize].score = v;
            }
            let out_count = self.nodes[v as usize].edges.out_count as usize;
            for vk in 0..out_count {
                if self.is_deleted_edge(v, vk) {
                    continue;
                }
                let w = self.nodes[v as usize].edges.arcs[vk].target;
                assert!(!self.nodes[w as usize].visited);
                assert_ne!(self.nodes[w as usize].score, 0);
                self.nodes[w as usize].score -= 1;
                if self.nodes[w as usize].score == 0 {
                    self.queue.push_back(w);
                }
            }
        }
        let mut v = self.nodes[self.end as usize].score;
        path.clear();
        while v != self.start {
            path.push(v);
            v = self.nodes[v as usize].score;
        }
        assert_eq!(v, self.start);
        let length = path.len();
        let half = length / 2;
        for k in 0..half {
            path.swap(k, length - 1 - k);
        }
    }
}

/// `pop_trac_bpc_f`: decodes one cell into operation, base and length. The
/// base fields the cell does not carry keep their previous value, exactly like
/// the reference macro.
fn pop_cell_bases(cell: u16, base_query: &mut u16, base_target: &mut u16, length: &mut u32) -> u16 {
    let operation = cell >> 14;
    match operation {
        2 | 3 => {
            *base_target = (cell >> 12) & 3;
            *length = u32::from(cell & 0xfff);
        }
        1 => {
            *base_target = (cell >> 12) & 3;
            *base_query = (cell >> 10) & 3;
            *length = u32::from(cell & 0x3ff);
        }
        _ => {
            *length = u32::from(cell & 0x3fff);
        }
    }
    operation
}

/// `push_trace_bp_f`: appends one base-aware run, merging with the previous
/// cell when it carries the same operation and bases.
fn push_trace_bp_f(
    cells: &mut Vec<u16>,
    operation: u16,
    mut base_query: u16,
    mut base_target: u16,
    mut length: u32,
    is_append: bool,
) {
    if operation == 3 {
        base_target = base_query;
        base_query = u16::MAX;
    }
    if is_append && !cells.is_empty() {
        let mut last_query = base_query;
        let mut last_target = base_target;
        let mut last_length = 0u32;
        let last_operation = pop_cell_bases(
            *cells.last().unwrap(),
            &mut last_query,
            &mut last_target,
            &mut last_length,
        );
        if operation == last_operation && base_query == last_query && base_target == last_target {
            cells.pop();
            length += last_length;
        }
    }
    let operation_shifted = u32::from(operation) << 14;
    let (mask, encoded) = match operation {
        2 | 3 => (
            0xfff_u32,
            operation_shifted + (u32::from(base_target & 3) << 12),
        ),
        1 => (
            0x3ff,
            operation_shifted
                + (u32::from(base_target & 3) << 12)
                + (u32::from(base_query & 3) << 10),
        ),
        _ => (0x3fff, operation_shifted),
    };
    while length >= mask {
        cells.push((encoded + mask) as u16);
        length -= mask;
    }
    if length > 0 {
        cells.push((encoded + length) as u16);
    }
}

/// `pop_trace_bp_f`: decodes a base-aware run starting at `index`, coalescing
/// following cells that carry the same operation and bases.
fn pop_trace_bp_f_full(cells: &[u16], mut index: usize) -> (usize, u16, u16, u16, u32) {
    let operation = cells[index] >> 14;
    let mut base_query = u16::MAX;
    let mut base_target = u16::MAX;
    let mut length;
    match operation {
        2 | 3 => {
            base_target = (cells[index] >> 12) & 3;
            length = u32::from(cells[index] & 0xfff);
        }
        1 => {
            base_target = (cells[index] >> 12) & 3;
            base_query = (cells[index] >> 10) & 3;
            length = u32::from(cells[index] & 0x3ff);
        }
        _ => {
            length = u32::from(cells[index] & 0x3fff);
        }
    }
    index += 1;
    while index < cells.len() && cells[index] >> 14 == operation {
        let mut next_query = u16::MAX;
        let mut next_target = u16::MAX;
        let next_length = match operation {
            2 | 3 => {
                next_target = (cells[index] >> 12) & 3;
                u32::from(cells[index] & 0xfff)
            }
            1 => {
                next_target = (cells[index] >> 12) & 3;
                next_query = (cells[index] >> 10) & 3;
                u32::from(cells[index] & 0x3ff)
            }
            _ => u32::from(cells[index] & 0x3fff),
        };
        if base_query != next_query || base_target != next_target {
            break;
        }
        length += next_length;
        index += 1;
    }
    if operation == 3 {
        base_query = base_target;
        base_target = u16::MAX;
    }
    (index, operation, base_query, base_target, length)
}

/// `pop_trace`: decodes a plain run starting at `index`, coalescing following
/// cells with the same operation.
fn pop_trace(cells: &[u16], mut index: usize) -> (usize, u16, u32) {
    let operation = cells[index] >> 14;
    let mut length = u32::from(cells[index] & 0x3fff);
    index += 1;
    while index < cells.len() && cells[index] >> 14 == operation {
        length += u32::from(cells[index] & 0x3fff);
        index += 1;
    }
    (index, operation, length)
}

/// `push_cns_c0`: adds the backbone path `[s, e)` after `s0` and returns the
/// last node it touched.
fn push_cns_c0(graph: &mut CnsGfa, s0: u64, s: u64, e: u64, plus0: u32) -> u64 {
    if s > e {
        return s0;
    }
    if !graph.add_arc(s0 as u32, s as u32, true, plus0) {
        graph.insert_arc(s0 as u32, s as u32, true, plus0);
        graph.insert_arc(s as u32, s0 as u32, false, plus0);
    } else {
        assert!(graph.add_arc(s as u32, s0 as u32, false, plus0));
    }
    graph.nodes[s as usize].score += 1;
    let mut re = s;
    let mut k = s + 1;
    while k < e {
        assert!(graph.add_arc((k - 1) as u32, k as u32, true, 1));
        assert!(graph.add_arc(k as u32, (k - 1) as u32, false, 1));
        graph.nodes[k as usize].score += 1;
        re = k;
        k += 1;
    }
    re
}

/// `trace_cns_bp`: walks the graph along `tstr`, reusing every node already on
/// the path, and returns how many bases it consumed plus the last node.
fn trace_cns_bp(
    graph: &mut CnsGfa,
    s0: u64,
    tstr: &[u8],
    tl: u64,
    values: &mut Vec<u32>,
    plus0: u32,
    max_trace: u64,
    avoid_backbone: bool,
) -> (u32, u32) {
    let mut reach = s0 as u32;
    if tl == 0 {
        return (0, reach);
    }
    values.clear();
    values.push(s0 as u32);
    values.push(UNSET);
    let mut nm = 2u32;
    let mut i = 0u32;
    let mut finished = false;
    while (i as u64) < tl && !finished {
        let base = SEQUENCE_NT6_TABLE[tstr[i as usize] as usize];
        let before = values.len() as u32;
        let mut bi = (before - nm) as usize;
        while bi < before as usize {
            let m = values[bi];
            let out_count = graph.nodes[m as usize].edges.out_count as usize;
            let mut k = 0usize;
            while k < out_count {
                let target = graph.nodes[m as usize].edges.arcs[k].target;
                if target == 0 || target == 1 {
                    k += 1;
                    continue;
                }
                if avoid_backbone && target >= graph.backbone_start && target < graph.backbone_end {
                    k += 1;
                    continue;
                }
                if target as usize >= graph.nodes.len() {
                    k += 1;
                    continue;
                }
                if graph.nodes[target as usize].base == base {
                    values.push(target);
                    values.push(bi as u32);
                    if (i as u64 + 1) == tl || values.len() as u64 > max_trace {
                        finished = true;
                        break;
                    }
                }
                k += 1;
            }
            if finished {
                break;
            }
            bi += 2;
        }
        if values.len() as u32 <= before {
            break;
        }
        nm = values.len() as u32 - before;
        i += 1;
    }
    if i > 0 && nm > 0 {
        let start = values.len() - nm as usize;
        reach = values[start];
        let mut bi = start;
        while values[bi + 1] != UNSET {
            let destination = values[bi];
            let source = values[values[bi + 1] as usize];
            let weight = if u64::from(source) != s0 { 1 } else { plus0 };
            assert!(graph.add_arc(source, destination, true, weight));
            assert!(graph.add_arc(destination, source, false, weight));
            graph.nodes[destination as usize].score += 1;
            bi = values[bi + 1] as usize;
        }
    } else {
        i = 0;
    }
    (i, reach)
}

/// `push_cns_c1`: follows `tstr` from `s0`, creating nodes where the graph has
/// no matching arc yet.
fn push_cns_c1(
    graph: &mut CnsGfa,
    s0: u64,
    tstr: &[u8],
    tl: u64,
    plus0: u32,
    values: &mut Vec<u32>,
    max_trace: u64,
) -> u64 {
    if tl == 0 {
        return s0;
    }
    let mut weight = plus0;
    let (mut k, mut reach) = trace_cns_bp(graph, s0, tstr, tl, values, plus0, max_trace, true);
    if k > 0 {
        weight = 1;
    }
    while (k as u64) < tl {
        reach = graph.add_arc_base(
            reach,
            SEQUENCE_NT6_TABLE[tstr[k as usize] as usize],
            weight,
            true,
        );
        weight = 1;
        k += 1;
    }
    u64::from(reach)
}

/// `append_cns_g`: dispatches one cigar cell of a candidate to the graph and
/// returns the node the path currently ends at.
#[allow(clippy::too_many_arguments)]
fn append_cns_g(
    graph: &mut CnsGfa,
    tstr: &[u8],
    tl: u64,
    qs: u64,
    qe: u64,
    cp: u16,
    cl: u64,
    pe: u64,
    values: &mut Vec<u32>,
    max_trace: u64,
    insert_position: i64,
) -> u64 {
    let mut s0 = pe;
    if pe == u64::MAX {
        if qs > 0 {
            s0 = qs - 1 + u64::from(graph.offset);
        } else {
            s0 = 0;
        }
    }
    let mut ns = qs + u64::from(graph.offset);
    let mut ne = qe + u64::from(graph.offset);
    if cp == 0 {
        if cl == 0 && qs == qe && (qe as i64) == insert_position {
            s0 = 0;
            ns = 1;
            ne = 1;
        }
        return push_cns_c0(graph, s0, ns, ne, 1);
    }
    if cp == 1 || cp == 2 {
        return push_cns_c1(graph, s0, tstr, tl, 1, values, max_trace);
    }
    s0
}

/// `rr_seq_t`: the recovered target sequence of the candidate being extracted,
/// cached between cells.
#[derive(Debug, Clone, Default)]
struct RrSequence {
    sequence: Vec<u8>,
    filled: usize,
    start: i64,
    end: i64,
    n0: i64,
    n1: i64,
    read_identifier: usize,
    strand: u8,
}

impl RrSequence {
    /// `get_sub_seq`: the target bases `[s, e)`, or nothing when the range is
    /// empty.
    fn get(&mut self, reads: &ReadStore, s: i64, e: i64) -> Option<Vec<u8>> {
        if s >= e {
            return None;
        }
        if s >= self.start && e <= self.end {
            let offset = (s - self.start) as usize;
            return Some(self.materialize(offset, (e - s) as usize));
        }
        let mut length = e - s;
        if self.start >= self.end {
            if length < self.n0 {
                length = self.n0;
            }
        } else if length < self.n1 {
            length = self.n1;
        }
        self.start = s;
        self.end = s + length;
        let read_length = reads.read_length(self.read_identifier) as i64;
        if self.end > read_length {
            self.end = read_length;
            length = self.end - self.start;
        }
        self.sequence = reads.read_sub_region(self.read_identifier, self.start, length, self.strand);
        self.filled = self.sequence.len();
        Some(self.materialize(0, (e - s) as usize))
    }

    fn materialize(&mut self, offset: usize, length: usize) -> Vec<u8> {
        let end = offset + length;
        if self.sequence.len() < end {
            self.sequence.resize(end, b'A');
        }
        self.sequence[offset..end].to_vec()
    }
}

/// `extract_sub_cigar_cns`: adds one candidate's path over the window
/// `[s, e)` to the consensus graph.
#[allow(clippy::too_many_arguments)]
fn extract_sub_cigar_cns(
    overlap: &OverlapRegion,
    reads: &ReadStore,
    s_input: i64,
    e_input: i64,
    insert_window_start: i64,
    insert_window_end: i64,
    s_end_input: bool,
    cache: &mut RrSequence,
    cursor: &mut UlOverlap,
    graph: &mut CnsGfa,
    values: &mut Vec<u32>,
    max_trace: u64,
) -> bool {
    let window_index = cursor.current_window as usize;
    let boundary = cursor.boundary as i64;
    let window = overlap.window_cigar.windows[window_index];
    let window_start = window.x_start as i64 + boundary;
    let window_end = window.x_end as i64 + 1 - boundary;
    let mut s = s_input;
    let mut e = e_input;
    if s < window_start {
        s = window_start;
    }
    if e > window_end {
        e = window_end;
    }
    if s > e {
        return false;
    }
    let os = s.max(window_start);
    let oe = e.min(window_end);
    if oe < os {
        return false;
    }
    let cigar_length = window.cigar_length as i64;
    if cigar_length == 0 {
        return false;
    }
    let mut xk = cursor.current_x_offset as i64;
    let mut yk = cursor.current_y_offset as i64;
    let mut ck = cursor.current_cigar_offset as i64;
    if ck < 0 || ck > cigar_length {
        ck = 0;
        xk = window.x_start as i64;
        yk = window.y_start as i64;
    }
    let cell = |index: i64| -> u16 { overlap.window_cigar.cells[(window.cigar_start as i64 + index) as usize] };
    while ck > 0 && xk >= s {
        ck -= 1;
        let operation = cell(ck) >> 14;
        if operation != 2 {
            xk -= i64::from(cell(ck) & 0x3fff);
        }
        if operation != 3 {
            yk -= i64::from(cell(ck) & 0x3fff);
        }
    }
    let s_end = !(!s_end_input && s == insert_window_start);
    let mut ii = [-1i64; 2];
    let mut it = [-1i64; 2];
    cache.start = 0;
    cache.end = 0;
    cache.n0 = e - s;
    cache.read_identifier = overlap.y_identifier as usize;
    cache.strand = overlap.y_strand as u8;
    if cache.n0 == 0 {
        cache.n0 = cache.n1;
    }
    let insert_position = if insert_window_start == insert_window_end {
        0
    } else {
        -1
    };
    let mut pe = u64::MAX;
    while ck < cigar_length && xk < e {
        let ws = xk;
        let wts = yk;
        let operation = cell(ck) >> 14;
        let mut length = i64::from(cell(ck) & 0x3fff);
        ck += 1;
        while ck < cigar_length && operation == cell(ck) >> 14 {
            length += i64::from(cell(ck) & 0x3fff);
            ck += 1;
        }
        if operation != 2 {
            xk += length;
        }
        if operation != 3 {
            yk += length;
        }
        let we = xk;
        let wte = yk;
        let os = s.max(ws);
        let oe = e.min(we);
        let overlap_length = if oe > os { oe - os } else { 0 };
        if s == e {
            if operation != 0 || ws >= s || we <= e || e != insert_window_end || s != insert_window_start
            {
                continue;
            }
        } else if operation != 2 {
            if overlap_length == 0 {
                continue;
            }
        } else if ws < s || ws >= e {
            continue;
        }
        if !s_end && operation == 2 && ws == s {
            continue;
        }
        let (ots, ote) = if operation < 2 {
            (os - ws + wts, oe - ws + wts)
        } else {
            (wts, wte)
        };
        if ii[0] == -1 {
            ii[0] = os;
            it[0] = ots;
        }
        ii[1] = oe;
        it[1] = ote;
        let length = if operation != 2 { oe - os } else { length };
        let target = if operation != 0 {
            cache.get(reads, ots, ote).unwrap_or_default()
        } else {
            Vec::new()
        };
        pe = append_cns_g(
            graph,
            &target,
            (ote - ots) as u64,
            (os - insert_window_start) as u64,
            (oe - insert_window_start) as u64,
            operation,
            length as u64,
            pe,
            values,
            max_trace,
            insert_position,
        );
    }
    while ck < cigar_length && xk <= e {
        let ws = xk;
        let wts = yk;
        let operation = cell(ck) >> 14;
        let mut length = i64::from(cell(ck) & 0x3fff);
        if operation != 2 {
            break;
        }
        ck += 1;
        while ck < cigar_length && operation == cell(ck) >> 14 {
            length += i64::from(cell(ck) & 0x3fff);
            ck += 1;
        }
        yk += length;
        let we = xk;
        let wte = yk;
        if ws >= s && ws <= e {
            let (ots, ote) = (wts, wte);
            if ii[0] == -1 {
                ii[0] = ws;
                it[0] = ots;
            }
            ii[1] = we;
            it[1] = ote;
            let target = cache.get(reads, ots, ote).unwrap_or_default();
            pe = append_cns_g(
                graph,
                &target,
                (ote - ots) as u64,
                (ws - insert_window_start) as u64,
                (we - insert_window_start) as u64,
                operation,
                length as u64,
                pe,
                values,
                max_trace,
                insert_position,
            );
        }
    }
    if ii[1] == -1 {
        return false;
    }
    let ae = if ii[1] == insert_window_end {
        1u64
    } else {
        (ii[1] + i64::from(graph.offset) - insert_window_start) as u64
    };
    if pe == u64::MAX {
        pe = 0;
    }
    if pe != ae {
        if !graph.add_arc(pe as u32, ae as u32, true, 1) {
            graph.insert_arc(pe as u32, ae as u32, true, 1);
            graph.insert_arc(ae as u32, pe as u32, false, 1);
        } else {
            graph.add_arc(ae as u32, pe as u32, false, 1);
        }
    }
    cursor.current_x_offset = xk as u32;
    cursor.current_y_offset = yk as u32;
    cursor.current_cigar_offset = ck as u32;
    cursor.current_y_length = 0;
    let _ = &mut it;
    true
}

/// `extract_sub_cigar_mm`: adds one candidate's vote for the query bases
/// `[s, e)` to the per-position counters.
///
/// The counters live at `counters[counter_base..]`; entry `2 * (t - s)` holds
/// the number of candidates that align query base `t` (the high half counts
/// matches), entry `2 * (t - s) + 1` the number of candidates with a target
/// base there.
fn extract_sub_cigar_mm(
    overlap: &OverlapRegion,
    s_input: i64,
    e_input: i64,
    cursor: &mut UlOverlap,
    counters: &mut [u64],
    counter_base: usize,
) -> i64 {
    let window_index = cursor.current_window as usize;
    let boundary = cursor.boundary as i64;
    let window = overlap.window_cigar.windows[window_index];
    let window_start = window.x_start as i64 + boundary;
    let window_end = window.x_end as i64 + 1 - boundary;
    let mut s = s_input;
    let mut e = e_input;
    if s < window_start {
        s = window_start;
    }
    if e > window_end {
        e = window_end;
    }
    if s >= e {
        return -1;
    }
    let os = s.max(window_start);
    let oe = e.min(window_end);
    if oe <= os {
        return -1;
    }
    let cigar_length = window.cigar_length as i64;
    if cigar_length == 0 {
        return -1;
    }
    let mut xk = cursor.current_x_offset as i64;
    let mut yk = cursor.current_y_offset as i64;
    let mut ck = cursor.current_cigar_offset as i64;
    if ck < 0 || ck > cigar_length {
        ck = 0;
        xk = window.x_start as i64;
        yk = window.y_start as i64;
    }
    let cell = |index: i64| -> u16 { overlap.window_cigar.cells[(window.cigar_start as i64 + index) as usize] };
    while ck > 0 && xk >= s {
        ck -= 1;
        let operation = cell(ck) >> 14;
        if operation != 2 {
            xk -= i64::from(cell(ck) & 0x3fff);
        }
        if operation != 3 {
            yk -= i64::from(cell(ck) & 0x3fff);
        }
    }
    while ck < cigar_length && xk < e {
        let ws = xk;
        let operation = cell(ck) >> 14;
        if operation != 2 {
            xk += i64::from(cell(ck) & 0x3fff);
        }
        if operation != 3 {
            yk += i64::from(cell(ck) & 0x3fff);
        }
        ck += 1;
        let we = xk;
        let os = s.max(ws);
        let oe = e.min(we);
        let overlap_length = if oe > os { oe - os } else { 0 };
        if operation != 2 {
            if overlap_length == 0 {
                continue;
            }
        } else if ws < s || ws >= e {
            continue;
        }
        if operation == 0 {
            let mut t = os + 1;
            while t < oe {
                let index = counter_base + ((t - s) as usize) * 2;
                counters[index] = counters[index].wrapping_add(1 + 0x1_0000_0000);
                counters[index + 1] = counters[index + 1].wrapping_add(1 + 0x1_0000_0000);
                t += 1;
            }
            let t = os;
            if t < oe {
                let index = counter_base + ((t - s) as usize) * 2;
                counters[index] = counters[index].wrapping_add(1 + 0x1_0000_0000);
                if os > ws {
                    counters[index + 1] = counters[index + 1].wrapping_add(1 + 0x1_0000_0000);
                }
            }
        } else if operation != 2 {
            let mut t = os + 1;
            while t < oe {
                let index = counter_base + ((t - s) as usize) * 2;
                counters[index] += 1;
                counters[index + 1] += 1;
                t += 1;
            }
            let t = os;
            if t < oe {
                let index = counter_base + ((t - s) as usize) * 2;
                counters[index] += 1;
                if os > ws {
                    counters[index + 1] += 1;
                }
            }
        } else {
            let index = counter_base + ((ws - s) as usize) * 2;
            counters[index + 1] += 1;
        }
    }
    cursor.current_x_offset = xk as u32;
    cursor.current_y_offset = yk as u32;
    cursor.current_cigar_offset = ck as u32;
    cursor.current_y_length = 0;
    1
}

/// `extract_sub_cigar_ii`: one candidate's packed vote over the gap `[s, e)`.
///
/// Returns `UNSET` when the candidate cannot vote.
fn extract_sub_cigar_ii(
    overlap: &OverlapRegion,
    reads: &ReadStore,
    query_length: i64,
    s_input: i64,
    e_input: i64,
    insert_window_start: i64,
    insert_window_end: i64,
    cursor: &mut UlOverlap,
) -> u32 {
    let window_index = cursor.current_window as usize;
    let boundary = cursor.boundary as i64;
    let window = overlap.window_cigar.windows[window_index];
    let window_start = window.x_start as i64 + boundary;
    let window_end = window.x_end as i64 + 1 - boundary;
    let mut s = s_input;
    let mut e = e_input;
    if s < window_start {
        s = window_start;
    }
    if e > window_end {
        e = window_end;
    }
    if s > e {
        return UNSET;
    }
    let os = s.max(window_start);
    let oe = e.min(window_end);
    if oe < os {
        return UNSET;
    }
    if !((window_start < insert_window_start || window_start == 0)
        && (window_end > insert_window_end || window_end == query_length))
    {
        return UNSET;
    }
    let cigar_length = window.cigar_length as i64;
    if cigar_length == 0 {
        return UNSET;
    }
    let mut xk = cursor.current_x_offset as i64;
    let mut yk = cursor.current_y_offset as i64;
    let mut ck = cursor.current_cigar_offset as i64;
    if ck < 0 || ck > cigar_length {
        ck = 0;
        xk = window.x_start as i64;
        yk = window.y_start as i64;
    }
    let cell = |index: i64| -> u16 { overlap.window_cigar.cells[(window.cigar_start as i64 + index) as usize] };
    while ck > 0 && xk >= s {
        ck -= 1;
        let operation = cell(ck) >> 14;
        if operation != 2 {
            xk -= i64::from(cell(ck) & 0x3fff);
        }
        if operation != 3 {
            yk -= i64::from(cell(ck) & 0x3fff);
        }
    }
    let mut ii = [-1i64; 2];
    let mut it = [-1i64; 2];
    let mut result = 0u32;
    let mut count = 0i64;
    while ck < cigar_length && xk < e {
        let ws = xk;
        let wts = yk;
        let operation = cell(ck) >> 14;
        let mut length = i64::from(cell(ck) & 0x3fff);
        if operation != 2 {
            xk += length;
        }
        if operation != 3 {
            yk += length;
        }
        ck += 1;
        let we = xk;
        let wte = yk;
        let os = s.max(ws);
        let oe = e.min(we);
        let overlap_length = if oe > os { oe - os } else { 0 };
        if s == e {
            if operation != 0 || ws >= s || we <= e || e != insert_window_end || s != insert_window_start
            {
                continue;
            }
        } else if operation != 2 {
            if overlap_length == 0 {
                continue;
            }
        } else if ws < s || ws >= e {
            continue;
        }
        if ii[0] == -1 {
            ii[0] = os;
            it[0] = if operation < 2 { os - ws + wts } else { wts };
        }
        ii[1] = oe;
        it[1] = if operation < 2 { oe - ws + wts } else { wte };
        if operation != 2 {
            length = oe - os;
        }
        count += length;
        if count <= SIMPLE_VOTE_LENGTH as i64 {
            for _ in 0..length {
                result <<= 2;
                result |= u32::from(operation);
            }
        }
    }
    while ck < cigar_length && xk <= e {
        let ws = xk;
        let wts = yk;
        let operation = cell(ck) >> 14;
        let length = i64::from(cell(ck) & 0x3fff);
        if operation != 2 {
            break;
        }
        yk += length;
        ck += 1;
        let we = xk;
        let wte = yk;
        if ws >= s && ws <= e {
            if ii[0] == -1 {
                ii[0] = ws;
                it[0] = wts;
            }
            ii[1] = we;
            it[1] = wte;
            count += length;
            if count <= SIMPLE_VOTE_LENGTH as i64 {
                for _ in 0..length {
                    result <<= 2;
                    result |= u32::from(operation);
                }
            }
        }
    }
    let mut packed = UNSET;
    if count <= SIMPLE_VOTE_LENGTH as i64
        && ii[1] >= ii[0]
        && ii[1] - ii[0] <= SIMPLE_VOTE_LENGTH as i64
        && it[1] >= it[0]
        && it[1] - it[0] <= SIMPLE_VOTE_LENGTH as i64
    {
        if ii[0] == insert_window_start && ii[1] == insert_window_end {
            let mut value = u16::from(count as u8);
            value <<= 12;
            result |= u32::from(value);
            result <<= 16;
            let base_count = it[1] - it[0];
            let mut bases = 0u16;
            if base_count > 0 {
                let target = reads.read_sub_region(
                    overlap.y_identifier as usize,
                    it[0],
                    base_count,
                    overlap.y_strand as u8,
                );
                for &base in target.iter() {
                    bases <<= 2;
                    bases |= u16::from(SEQUENCE_NT6_TABLE[base as usize]);
                }
            }
            result |= u32::from(bases);
            let length_bits = u16::from(base_count as u8);
            result |= u32::from(length_bits << 12);
            packed = result;
        }
    }
    cursor.current_x_offset = xk as u32;
    cursor.current_y_offset = yk as u32;
    cursor.current_cigar_offset = ck as u32;
    cursor.current_y_length = 0;
    let _ = &mut it;
    packed
}

/// `cc_idx_t`: the cursor state of one candidate index.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct CcIndex {
    /// Scan position, `i`.
    scan: i64,
    /// Number of sorted candidates, `srt_n`.
    sorted_count: i64,
    /// Set when a candidate ends inside the window, `rr`.
    reduce_next: i64,
    /// Set when [`ConsensusScratch::cns_gen0`] failed, `ru`.
    skip_next: i64,
    /// Start of the trailing anchor, `mms`.
    anchor_start: u64,
    /// End of the trailing anchor, `mme`.
    anchor_end: u64,
}

/// `iter_cc_idx_t`: returns the candidate ids that overlap `[s, e)`, appending
/// them to `index`.
fn iter_cc_idx_t(
    overlaps: &[OverlapRegion],
    candidates: &[UlOverlap],
    index: &mut CcIndex,
    entries: &mut Vec<u64>,
    s: i64,
    e: i64,
    is_reduce: bool,
    is_insert: bool,
) -> Vec<u32> {
    if index.skip_next == 0 {
        if is_reduce {
            let mut destination = index.sorted_count as usize;
            let mut m = index.sorted_count as usize;
            while m < entries.len() {
                let candidate = entries[m] as u32 as usize;
                let cursor = candidates[candidate];
                let window = overlaps[cursor.overlap_identifier as usize].window_cigar.windows
                    [cursor.current_window as usize];
                let q0 = window.x_start as i64 + cursor.boundary as i64;
                let q1 = window.x_end as i64 + 1 - cursor.boundary as i64;
                let os = q0.max(s);
                let oe = q1.min(e);
                if oe > os || (is_insert && s == e && s >= q0 && s <= q1) {
                    entries[destination] = entries[m];
                    destination += 1;
                }
                m += 1;
            }
            entries.truncate(destination);
        }
        let sorted_count = index.sorted_count as usize;
        let mut i = index.scan as usize;
        while i < sorted_count {
            let candidate = entries[i] as u32 as usize;
            let cursor = candidates[candidate];
            let window = overlaps[cursor.overlap_identifier as usize].window_cigar.windows
                [cursor.current_window as usize];
            let q0 = window.x_start as i64 + cursor.boundary as i64;
            let q1 = window.x_end as i64 + 1 - cursor.boundary as i64;
            if q0 > e {
                break;
            }
            if !is_insert && q0 >= e {
                break;
            }
            let os = q0.max(s);
            let oe = q1.min(e);
            if oe > os || (is_insert && s == e && s >= q0 && s <= q1) {
                let value = entries[i];
                entries.push(u64::from(value as u32));
            }
            i += 1;
        }
        index.scan = i as i64;
    } else {
        index.skip_next = 0;
    }
    entries[index.sorted_count as usize..]
        .iter()
        .map(|value| *value as u32)
        .collect()
}

/// The cell buffer writer of one window, `window_list` plus its shared cell
/// vector.
///
/// The reference keeps the window length in `window_list::clen` and refreshes
/// it by hand after every push, never after a truncation. Truncating without
/// refreshing is observable: the next push still merges with the cell that
/// preceded the window. [`Self::length`] therefore tracks the reference's
/// `clen` instead of deriving it from the cell vector.
struct WindowWriter<'a> {
    cigar: &'a mut WindowCigar,
    window_index: usize,
    length: u32,
}

impl<'a> WindowWriter<'a> {
    fn new(cigar: &'a mut WindowCigar, window_index: usize) -> Self {
        let length = cigar.windows[window_index].cigar_length;
        Self {
            cigar,
            window_index,
            length,
        }
    }

    fn window_length(&self) -> u32 {
        self.length
    }

    fn len(&self) -> usize {
        self.cigar.cells.len()
    }

    fn sync(&mut self) {
        let start = self.cigar.windows[self.window_index].cigar_start as usize;
        self.length = (self.cigar.cells.len() - start) as u32;
        self.cigar.windows[self.window_index].cigar_length = self.length;
    }

    fn push(&mut self, operation: u16, base_query: u16, base_target: u16, length: u32) {
        let append = self.window_length() > 0;
        push_trace_bp_f(
            &mut self.cigar.cells,
            operation,
            base_query,
            base_target,
            length,
            append,
        );
        self.sync();
    }

    /// `sc->n = cn0`: the reference shortens the vector but leaves `clen`
    /// alone, so [`Self::length`] is deliberately not refreshed here.
    fn truncate(&mut self, length: usize) {
        self.cigar.cells.truncate(length);
    }

    fn cell(&self, index: usize) -> u16 {
        self.cigar.cells[index]
    }
}

/// `push_correct0_fhc`: appends one match run or one decoded vote.
fn push_correct0_fhc(writer: &mut WindowWriter, length0: u32, packed: u32, query: &[u8]) -> u64 {
    let mut nec = 0u64;
    if length0 != u32::MAX {
        writer.push(0, u16::MAX, u16::MAX, length0);
    } else if packed != u32::MAX {
        let cigar_ops = (packed << 4) >> 20;
        let cigar_count = packed >> 28;
        let base_ops = (packed << 20) >> 20;
        let base_count = (packed << 16) >> 28;
        let mut target_index = 0u32;
        let mut query_index = 0u32;
        for k in 0..cigar_count {
            let shift = (cigar_count - 1 - k) << 1;
            let operation = ((cigar_ops >> shift) & 3) as u16;
            let mut base_query = u16::MAX;
            let mut base_target = u16::MAX;
            if operation != 3 {
                let shift = (base_count - 1 - target_index) << 1;
                base_target = ((base_ops >> shift) & 3) as u16;
                target_index += 1;
            }
            if operation != 2 {
                base_query = u16::from(SEQUENCE_NT6_TABLE[query[query_index as usize] as usize]);
                query_index += 1;
            }
            writer.push(operation, base_query, base_target, 1);
            if operation != 0 {
                nec += 1;
            }
        }
    }
    nec
}

/// `push_correct1_fhc_indel_exz`: rewrites a deletion that follows an
/// insertion into the cheapest alignment of the two sequences, when the
/// reference alignment allows it.
#[allow(clippy::too_many_arguments)]
fn push_correct1_fhc_indel_exz(
    writer: &mut WindowWriter,
    cell_start: usize,
    graph: &CnsGfa,
    text: &[u8],
    exz: &mut BitExtz,
    genome_begin: u32,
    c0: u16,
    cl0: i64,
    ok0_input: i64,
) -> u64 {
    assert_eq!(c0, 3);
    let mut ok = 0i64;
    let mut nk = 0i64;
    let mut e0 = 0i64;
    ok += cl0;
    e0 += cl0;
    let mut ck = writer.len() as i64 - 1;
    let mut f = 0u32;
    while ck >= cell_start as i64 {
        let packed = writer.cell(ck as usize);
        let operation = packed >> 14;
        if operation == 0 {
            break;
        }
        let length = match operation {
            2 | 3 => u32::from(packed & 0xfff),
            1 => u32::from(packed & 0x3ff),
            _ => u32::from(packed & 0x3fff),
        };
        let length = i64::from(length);
        if operation != 2 {
            ok += length;
        }
        if operation != 3 {
            nk += length;
        }
        if operation != 0 {
            e0 += length;
        }
        if c0 != operation {
            f = 1;
        }
        ck -= 1;
    }
    let cn0 = ck + 1;
    let cn = writer.len() as i64;
    let mut nec = 0u64;
    if f == 0 || ok == 0 || nk == 0 {
        let mut index = ok0_input + i64::from(genome_begin);
        for _ in 0..cl0 {
            writer.push(c0, u16::from(graph.nodes[index as usize].base), u16::MAX, 1);
            nec += 1;
            index += 1;
        }
        return nec;
    }
    let mut ok0 = ok0_input + cl0;
    ok0 -= ok;
    assert!(ok0 >= 0);
    let oseq = &text[ok0 as usize..(ok0 + ok) as usize];
    let ol = ok;
    let nl = nk;
    let mut nseq = vec![0u8; nl as usize];
    let mut index = cn0;
    let mut x = 0i64;
    let mut y = 0i64;
    while index < cn {
        let (next, operation, _base_query, base_target, length) =
            pop_trace_bp_f_full(&writer.cigar.cells, index as usize);
        let x0 = x;
        let y0 = y;
        index = next as i64;
        if operation != 2 {
            x += i64::from(length);
        }
        if operation != 3 {
            y += i64::from(length);
        }
        if operation == 0 {
            for k in 0..(x - x0) {
                nseq[(y0 + k) as usize] = oseq[(x0 + k) as usize];
            }
        } else if operation == 1 || operation == 2 {
            for k in y0..y {
                nseq[k as usize] = BASE_CHARACTERS[base_target as usize];
            }
        }
    }
    assert_eq!(x + cl0, ol);
    assert_eq!(y, nl);
    if nl == ol && nl == 1 {
        writer.truncate(cn0 as usize);
        if oseq[0] == nseq[0] {
            writer.push(0, u16::MAX, u16::MAX, 1);
        } else {
            writer.push(
                1,
                u16::from(SEQUENCE_NT6_TABLE[oseq[0] as usize]),
                u16::from(SEQUENCE_NT6_TABLE[nseq[0] as usize]),
                1,
            );
            nec += 1;
        }
        return nec;
    }
    let ml = ol.max(nl);
    let mut diff = 31i64.min(ml);
    let diff0 = diff;
    exz.clear_align();
    cal_exz_global(exz, &nseq, oseq, diff as i32);
    let mut aligned = exz.is_aligned();
    if !aligned {
        diff = 63i64.min(ml);
        if diff > diff0 {
            exz.clear_align();
            cal_exz_global(exz, &nseq, oseq, diff as i32);
            aligned = exz.is_aligned();
        }
    }
    if aligned && i64::from(exz.result.error) < e0 {
        writer.truncate(cn0 as usize);
        let cn = exz.cigar.packed.len();
        let mut x = 0i64;
        let mut y = 0i64;
        let mut ck = 0usize;
        while ck < cn {
            let x0 = x;
            let y0 = y;
            let (next, operation, length) = pop_trace(&exz.cigar.packed, ck);
            ck = next;
            if operation != 2 {
                x += i64::from(length);
            }
            if operation != 3 {
                y += i64::from(length);
            }
            match operation {
                0 => {
                    writer.push(0, u16::MAX, u16::MAX, length);
                }
                1 => {
                    for k in 0..u32::from(length) {
                        writer.push(
                            1,
                            u16::from(SEQUENCE_NT6_TABLE[oseq[(x0 + i64::from(k)) as usize] as usize]),
                            u16::from(SEQUENCE_NT6_TABLE[nseq[(y0 + i64::from(k)) as usize] as usize]),
                            1,
                        );
                        nec += 1;
                    }
                }
                2 => {
                    for k in 0..u32::from(length) {
                        writer.push(
                            2,
                            u16::MAX,
                            u16::from(SEQUENCE_NT6_TABLE[nseq[(y0 + i64::from(k)) as usize] as usize]),
                            1,
                        );
                        nec += 1;
                    }
                }
                3 => {
                    for k in 0..u32::from(length) {
                        writer.push(
                            3,
                            u16::from(SEQUENCE_NT6_TABLE[oseq[(x0 + i64::from(k)) as usize] as usize]),
                            u16::MAX,
                            1,
                        );
                        nec += 1;
                    }
                }
                _ => {}
            }
        }
    } else if ml < e0 {
        writer.truncate(cn0 as usize);
        let ml0 = ol.min(nl);
        for k in 0..ml0 {
            writer.push(
                1,
                u16::from(SEQUENCE_NT6_TABLE[oseq[k as usize] as usize]),
                u16::from(SEQUENCE_NT6_TABLE[nseq[k as usize] as usize]),
                1,
            );
            nec += 1;
        }
        if ol > ml0 {
            for k in ml0..ol {
                writer.push(
                    3,
                    u16::from(SEQUENCE_NT6_TABLE[oseq[k as usize] as usize]),
                    u16::MAX,
                    1,
                );
                nec += 1;
            }
        } else if nl > ml0 {
            for k in ml0..nl {
                writer.push(
                    2,
                    u16::MAX,
                    u16::from(SEQUENCE_NT6_TABLE[nseq[k as usize] as usize]),
                    1,
                );
                nec += 1;
            }
        }
    } else {
        let mut index = ok0 + i64::from(genome_begin);
        for _ in 0..cl0 {
            writer.push(c0, u16::from(graph.nodes[index as usize].base), u16::MAX, 1);
            nec += 1;
            index += 1;
        }
    }
    nec
}

/// `push_correct1_fhc`: writes the consensus path `path` as cigar cells,
/// squeezing the deletions against the backbone.
fn push_correct1_fhc(
    writer: &mut WindowWriter,
    graph: &CnsGfa,
    query: &[u8],
    exz: &mut BitExtz,
    path: &[u32],
    block_length: u32,
) -> u64 {
    let mut nec = 0u64;
    let begin = graph.offset;
    let end = block_length + graph.offset;
    let mut bend = graph.offset;
    let mut is_insertion = false;
    let cell_start = writer.len();
    if !path.is_empty() {
        let count = path.len();
        let mut l = 0usize;
        for k in 1..=count {
            let mut split = false;
            let mut is_match = false;
            if k == count {
                if l < count {
                    is_match = path[l] >= begin && path[l] < end;
                }
                split = true;
            } else {
                is_match = path[l] >= begin && path[l] < end;
                let next_match = path[k] >= begin && path[k] < end;
                if is_match != next_match {
                    split = true;
                } else if is_match && path[k] - path[l] != (k - l) as u32 {
                    split = true;
                }
            }
            if !split {
                continue;
            }
            if is_match {
                if path[l] > bend {
                    if is_insertion {
                        nec += push_correct1_fhc_indel_exz(
                            writer,
                            cell_start,
                            graph,
                            query,
                            exz,
                            graph.offset,
                            3,
                            i64::from(path[l] - bend),
                            i64::from(bend - graph.offset),
                        );
                    } else {
                        for i in bend..path[l] {
                            writer.push(3, u16::from(graph.nodes[i as usize].base), u16::MAX, 1);
                            nec += 1;
                        }
                    }
                }
                writer.push(0, u16::MAX, u16::MAX, path[k - 1] + 1 - path[l]);
                bend = path[k - 1] + 1;
                is_insertion = false;
            } else {
                for i in l..k {
                    writer.push(2, u16::MAX, u16::from(graph.nodes[path[i] as usize].base), 1);
                    nec += 1;
                }
                is_insertion = true;
            }
            l = k;
        }
    }
    if end > bend {
        if is_insertion {
            nec += push_correct1_fhc_indel_exz(
                writer,
                cell_start,
                graph,
                query,
                exz,
                graph.offset,
                3,
                i64::from(end - bend),
                i64::from(bend - graph.offset),
            );
        } else {
            for i in bend..end {
                writer.push(3, u16::from(graph.nodes[i as usize].base), u16::MAX, 1);
                nec += 1;
            }
        }
    }
    nec
}

/// The reference's `x_length` of a base-aware cigar, `cal_cigar_xlen_fhc`.
pub fn cigar_x_length(cells: &[u16]) -> u32 {
    let mut index = 0usize;
    let mut x = 0u32;
    while index < cells.len() {
        let (next, operation, _base_query, _base_target, length) = pop_trace_bp_f_full(cells, index);
        index = next;
        if operation != 2 {
            x += length;
        }
    }
    x
}

/// The scratch state of one [`ConsensusScratch::wcns_gen`] run.
///
/// The reference keeps two independent candidate indices over the same cursor
/// array: `ii_a` sorts and scans the vote candidates (`idx`), while `ii_b`
/// scans the per-position counters and their own copy of the sorted keys
/// (`buf`). They advance separately, so they are separate vectors here too.
#[derive(Debug, Clone, Default)]
pub struct ConsensusScratch {
    /// Candidate windows, `c_idx`.
    candidates: Vec<UlOverlap>,
    /// `ii_a`: the vote candidate cursor.
    index_a: CcIndex,
    /// The sorted vote keys and their working set, `idx->idx`.
    entries: Vec<u64>,
    /// `ii_b`'s per-position counters, the first `2 * wl` entries of `buf`.
    counters: Vec<u64>,
    /// `ii_b`: the gap candidate cursor.
    index_b: CcIndex,
    /// `ii_b`'s copy of the sorted keys and its working set, `buf` past the
    /// counters.
    counter_entries: Vec<u64>,
    /// The packed votes and the graph path, `b32`.
    values: Vec<u32>,
    /// The consensus graph, `cns`.
    graph: CnsGfa,
    /// The auxiliary overlap the cigar is written to, `aux_o`.
    pub auxiliary: WindowCigar,
}

impl ConsensusScratch {
    pub fn new() -> Self {
        Self::default()
    }

    /// `cns_gen0`: votes over the short gap `[s, e)` and returns the packed
    /// winner.
    fn cns_gen0(
        &mut self,
        overlaps: &[OverlapRegion],
        reads: &ReadStore,
        s: u64,
        e: u64,
        query_length: u64,
        occurrence_total: u64,
        occurrence_exact: f64,
    ) -> Option<u32> {
        if e > s + SIMPLE_VOTE_LENGTH {
            return None;
        }
        self.values.clear();
        let is_reduce = self.index_b.reduce_next != 0;
        let ids = iter_cc_idx_t(
            overlaps,
            &self.candidates,
            &mut self.index_b,
            &mut self.counter_entries,
            s as i64,
            e as i64,
            is_reduce,
            s == e,
        );
        self.index_b.reduce_next = 0;
        let mut aligned = 0u64;
        for &id in &ids {
            let cursor = self.candidates[id as usize];
            let overlap = &overlaps[cursor.overlap_identifier as usize];
            let window = overlap.window_cigar.windows[cursor.current_window as usize];
            let q0 = window.x_start as i64 + cursor.boundary as i64;
            let q1 = window.x_end as i64 + 1 - cursor.boundary as i64;
            if q1 <= e as i64 {
                self.index_b.reduce_next = 1;
            }
            let os = q0.max(s as i64);
            let oe = q1.min(e as i64);
            if oe > os || (s == e && (s as i64) > q0 && (s as i64) < q1) {
                aligned += 1;
                let packed = extract_sub_cigar_ii(
                    overlap,
                    reads,
                    query_length as i64,
                    os,
                    oe,
                    s as i64,
                    e as i64,
                    &mut self.candidates[id as usize],
                );
                if packed != UNSET {
                    self.values.push(packed);
                }
            }
        }
        let mut count = [self.values.len() as u64, aligned + 1];
        if count[0] as f64 > (count[1] as f64) * occurrence_exact
            && count[0] > count[1].wrapping_sub(count[0])
            && count[1] >= occurrence_total
            && count[0] > 1
        {
            radix_sort_by_key_32(&mut self.values, &|value: &u32| *value);
            let mut best = 0usize;
            let mut best_start = 0usize;
            let mut l = 0usize;
            for k in 1..=self.values.len() {
                if k == self.values.len() || self.values[k] != self.values[l] {
                    if k - l > best {
                        best = k - l;
                        best_start = l;
                    }
                    l = k;
                }
            }
            count[0] = best as u64;
            if count[0] as f64 > (count[1] as f64) * occurrence_exact
                && count[0] > count[1].wrapping_sub(count[0])
                && count[1] >= occurrence_total
                && count[0] > 1
            {
                return Some(self.values[best_start]);
            }
        }
        self.index_b.skip_next = 1;
        None
    }

    /// `cns_gen_full0`: builds and walks the consensus graph over `[s, e)`.
    #[allow(clippy::too_many_arguments)]
    fn cns_gen_full0(
        &mut self,
        overlaps: &[OverlapRegion],
        reads: &ReadStore,
        query: &[u8],
        exz: &mut BitExtz,
        s: u64,
        e: u64,
        s_end: bool,
        _occurrence_total: u64,
        _occurrence_exact: f64,
        max_trace: u64,
        window_index: usize,
    ) -> u64 {
        self.graph
            .init(&query[s as usize..e as usize]);
        self.values.clear();
        let mut cache = RrSequence {
            n0: 32,
            n1: 32,
            ..RrSequence::default()
        };
        let is_reduce = self.index_b.reduce_next != 0;
        let ids = iter_cc_idx_t(
            overlaps,
            &self.candidates,
            &mut self.index_b,
            &mut self.counter_entries,
            s as i64,
            e as i64,
            is_reduce,
            s == e,
        );
        self.index_b.reduce_next = 0;
        for &id in &ids {
            let cursor = self.candidates[id as usize];
            let overlap = &overlaps[cursor.overlap_identifier as usize];
            let window = overlap.window_cigar.windows[cursor.current_window as usize];
            let q0 = window.x_start as i64 + cursor.boundary as i64;
            let q1 = window.x_end as i64 + 1 - cursor.boundary as i64;
            if q1 <= e as i64 {
                self.index_b.reduce_next = 1;
            }
            let os = q0.max(s as i64);
            let oe = q1.min(e as i64);
            if oe > os || (s == e && (s as i64) > q0 && (s as i64) < q1) {
                extract_sub_cigar_cns(
                    overlap,
                    reads,
                    os,
                    oe,
                    s as i64,
                    e as i64,
                    s_end,
                    &mut cache,
                    &mut self.candidates[id as usize],
                    &mut self.graph,
                    &mut self.values,
                    max_trace,
                );
            }
        }
        let mut scratch = Vec::new();
        self.graph.refine(&mut scratch);
        self.graph.generate_sequence(&mut self.values);
        let mut writer = WindowWriter::new(&mut self.auxiliary, window_index);
        push_correct1_fhc(
            &mut writer,
            &self.graph,
            &query[s as usize..],
            exz,
            &self.values,
            (e - s) as u32,
        )
    }

    /// `cns_gen_full`: fills the long gap `[s0, e0)` window by window.
    #[allow(clippy::too_many_arguments)]
    fn cns_gen_full(
        &mut self,
        overlaps: &[OverlapRegion],
        reads: &ReadStore,
        s0: u64,
        e0: u64,
        window_length: u64,
        query: &[u8],
        exz: &mut BitExtz,
        occurrence_total: u64,
        occurrence_exact: f64,
        max_trace: u64,
        window_index: usize,
    ) -> u64 {
        let mut nec = 0u64;
        if e0 - s0 <= window_length {
            nec += self.cns_gen_full0(
                overlaps,
                reads,
                query,
                exz,
                s0,
                e0,
                true,
                occurrence_total,
                occurrence_exact,
                max_trace,
                window_index,
            );
            return nec;
        }
        let mut s = s0;
        let mut e = (s0 + window_length).min(e0);
        while s < e0 {
            nec += self.cns_gen_full0(
                overlaps,
                reads,
                query,
                exz,
                s,
                e,
                s == s0,
                occurrence_total,
                occurrence_exact,
                max_trace,
                window_index,
            );
            s += window_length;
            e += window_length;
            e = e.min(e0);
        }
        nec
    }

    /// `push_cns_anchor`: emits the anchor `[s, e)` and, before it, the gap
    /// left by the previous anchor.
    #[allow(clippy::too_many_arguments)]
    fn push_cns_anchor(
        &mut self,
        overlaps: &[OverlapRegion],
        reads: &ReadStore,
        s: u64,
        e: u64,
        query: &[u8],
        exz: &mut BitExtz,
        is_tail: bool,
        occurrence_total: u64,
        occurrence_exact: f64,
    ) -> u64 {
        if !is_tail && s >= e {
            return 0;
        }
        let mut e0 = 0u64;
        if !self.auxiliary.windows.is_empty() {
            e0 = self.auxiliary.windows[self.auxiliary.windows.len() - 1].x_end as i64 as u64 + 1;
        }
        assert!(s >= e0);
        if s == e && is_tail && s == e0 {
            return 0;
        }
        if self.auxiliary.windows.is_empty() {
            self.auxiliary.append_window(-1, -1, 0, 0, 0, 0, 0, 0);
        }
        let window_index = self.auxiliary.windows.len() - 1;
        let mut nec = 0u64;
        let window_length = self.graph.consensus_window_length as u64;
        let query_length = query.len() as u64;
        if (!is_tail && s > 0) || (is_tail && s > e0) {
            let packed = self.cns_gen0(
                overlaps,
                reads,
                e0,
                s,
                query_length,
                occurrence_total,
                occurrence_exact,
            );
            if let Some(packed) = packed {
                if self.auxiliary.windows[window_index].x_start == -1
                    || self.auxiliary.windows[window_index].x_end == -1
                {
                    self.auxiliary.windows[window_index].x_start = e0 as i32;
                    self.auxiliary.windows[window_index].x_end = (s - 1) as i32;
                }
                let mut writer = WindowWriter::new(&mut self.auxiliary, window_index);
                nec += push_correct0_fhc(&mut writer, u32::MAX, packed, &query[e0 as usize..]);
            } else {
                nec += self.cns_gen_full(
                    overlaps,
                    reads,
                    e0,
                    s,
                    window_length,
                    query,
                    exz,
                    occurrence_total,
                    occurrence_exact,
                    query_length,
                    window_index,
                );
            }
            self.auxiliary.windows[window_index].x_end = (s - 1) as i32;
        }
        if self.auxiliary.windows[window_index].x_start == -1
            || self.auxiliary.windows[window_index].x_end == -1
        {
            self.auxiliary.windows[window_index].x_start = s as i32;
            self.auxiliary.windows[window_index].x_end = (e - 1) as i32;
        }
        let mut writer = WindowWriter::new(&mut self.auxiliary, window_index);
        nec += push_correct0_fhc(&mut writer, (e - s) as u32, u32::MAX, &[]);
        self.auxiliary.windows[window_index].x_end = (e - 1) as i32;
        nec
    }

    /// `wcns_vote`: votes over one window and flushes anchors into
    /// [`Self::push_cns_anchor`].
    #[allow(clippy::too_many_arguments)]
    fn wcns_vote(
        &mut self,
        overlaps: &[OverlapRegion],
        reads: &ReadStore,
        query: &[u8],
        exz: &mut BitExtz,
        candidate_ids: &[u32],
        s: i64,
        e: i64,
        occurrence_total: u64,
        occurrence_exact: f64,
    ) -> (i64, u64) {
        let mut rr = 0i64;
        let mut nec = 0u64;
        for &candidate in candidate_ids {
            let cursor = self.candidates[candidate as usize];
            let overlap = &overlaps[cursor.overlap_identifier as usize];
            let window = overlap.window_cigar.windows[cursor.current_window as usize];
            let q0 = window.x_start as i64 + cursor.boundary as i64;
            let q1 = window.x_end as i64 + 1 - cursor.boundary as i64;
            if q1 <= e {
                rr = 1;
            }
            let os = q0.max(s);
            let oe = q1.min(e);
            if oe > os {
                extract_sub_cigar_mm(
                    overlap,
                    os,
                    oe,
                    &mut self.candidates[candidate as usize],
                    &mut self.counters,
                    (os - s) as usize,
                );
            }
        }
        let window_length = (e - s) as usize;
        let mut os = self.index_b.anchor_start as i64;
        let mut oe = self.index_b.anchor_end as i64;
        for k in 0..window_length {
            let mut oc0 = (self.counters[k << 1] >> 32) + 1;
            let mut oc1 = (self.counters[k << 1] as u32 as u64) + 1;
            if vote_passes(oc0, oc1, occurrence_total, occurrence_exact)
                || oc1 < occurrence_total
            {
                let mut insertion = true;
                oc0 = (self.counters[(k << 1) + 1] >> 32) + 1;
                oc1 = (self.counters[(k << 1) + 1] as u32 as u64) + 1;
                if vote_passes(oc0, oc1, occurrence_total, occurrence_exact)
                    || oc1 < occurrence_total
                {
                    insertion = false;
                }
                if insertion {
                    if oe > os && os != -1 {
                        nec += self.push_cns_anchor(
                            overlaps,
                            reads,
                            os as u64,
                            oe as u64,
                            query,
                            exz,
                            false,
                            occurrence_total,
                            occurrence_exact,
                        );
                    }
                    os = -1;
                    oe = -1;
                }
                if s + k as i64 == oe {
                    oe += 1;
                } else {
                    if oe > os && os != -1 {
                        nec += self.push_cns_anchor(
                            overlaps,
                            reads,
                            os as u64,
                            oe as u64,
                            query,
                            exz,
                            false,
                            occurrence_total,
                            occurrence_exact,
                        );
                    }
                    os = s + k as i64;
                    oe = s + k as i64 + 1;
                }
            } else {
                if oe > os && os != -1 {
                    nec += self.push_cns_anchor(
                        overlaps,
                        reads,
                        os as u64,
                        oe as u64,
                        query,
                        exz,
                        false,
                        occurrence_total,
                        occurrence_exact,
                    );
                }
                os = -1;
                oe = -1;
            }
            self.counters[k << 1] = 0;
            self.counters[(k << 1) + 1] = 0;
        }
        self.index_b.anchor_start = u64::MAX;
        self.index_b.anchor_end = u64::MAX;
        if oe > os && os != -1 {
            self.index_b.anchor_start = os as u64;
            self.index_b.anchor_end = oe as u64;
        }
        (rr, nec)
    }

    /// `wcns_gen`: the windowed consensus over the whole query.
    #[allow(clippy::too_many_arguments)]
    pub fn wcns_gen(
        &mut self,
        overlaps: &mut [OverlapRegion],
        reads: &ReadStore,
        query: &[u8],
        exz: &mut BitExtz,
        boundary: i64,
        window_length: u64,
        query_length: i64,
        occurrence_total: u64,
        occurrence_exact: f64,
        consensus_window: u64,
    ) -> u64 {
        let overlap_count = overlaps.len();
        self.candidates.clear();
        self.entries.clear();
        for k in 0..overlap_count {
            let window_count = overlaps[k].window_cigar.windows.len();
            let mut large_indel = false;
            overlaps[k].without_large_indel = 0;
            if window_count == 0 || overlaps[k].is_match != 1 {
                continue;
            }
            let mut last_window: i64 = -1;
            for i in 0..window_count {
                let window = overlaps[k].window_cigar.windows[i];
                if is_unaligned_window(&window) {
                    let length = window.x_end as i64 + 1 - window.x_start as i64;
                    if length >= 6 {
                        large_indel = true;
                    }
                    continue;
                }
                let q0 = window.x_start as i64 + boundary;
                let q1 = window.x_end as i64 - boundary;
                if q1 >= q0 {
                    let key = ((q0 as u64) << 32) | self.candidates.len() as u64;
                    self.entries.push(key);
                    self.candidates.push(UlOverlap {
                        overlap_identifier: k as u32,
                        current_window: i as u32,
                        current_x_offset: window.x_start as u32,
                        current_y_offset: window.y_start as u32,
                        current_y_length: 0,
                        current_cigar_offset: 0,
                        boundary: (boundary as u32) & SEC_MODE,
                    });
                }
                if !large_indel {
                    if i == 0 {
                        let difference =
                            (window.x_start as i64 - overlaps[k].x_position_start as i64).abs();
                        if difference >= 6 {
                            large_indel = true;
                        }
                    }
                    if last_window != -1 {
                        let previous = overlaps[k].window_cigar.windows[last_window as usize];
                        let difference =
                            (window.x_start as i64 - (previous.x_end as i64 + 1)).abs();
                        if difference >= 6 {
                            large_indel = true;
                        }
                        let difference =
                            (window.y_start as i64 - (previous.y_end as i64 + 1)).abs();
                        if difference >= 6 {
                            large_indel = true;
                        }
                    }
                    if i + 1 == window_count {
                        let difference =
                            (window.x_end as i64 - overlaps[k].x_position_end as i64).abs();
                        if difference >= 6 {
                            large_indel = true;
                        }
                    }
                    if !large_indel {
                        let start = window.cigar_start as usize;
                        let length = window.cigar_length as usize;
                        let cells = &overlaps[k].window_cigar.cells[start..start + length];
                        let mut index = 0usize;
                        while index < cells.len() && !large_indel {
                            let (next, operation, run) = pop_trace(cells, index);
                            index = next;
                            if operation >= 2 && run >= 6 {
                                large_indel = true;
                            }
                        }
                    }
                }
                last_window = i as i64;
            }
            overlaps[k].without_large_indel = u8::from(!large_indel);
        }
        let sorted_count = self.entries.len();
        radix_sort_u64(&mut self.entries);
        let mut group_start = 0usize;
        for k in 1..sorted_count {
            if (self.entries[k] >> 32) != (self.entries[group_start] >> 32) {
                if k - group_start > 1 {
                    for t in group_start..k {
                        let candidate = self.entries[t] as u32 as usize;
                        let cursor = self.candidates[candidate];
                        let window = overlaps[cursor.overlap_identifier as usize].window_cigar.windows
                            [cursor.current_window as usize];
                        let q1 = window.x_end as i64 + 1 - cursor.boundary as i64;
                        self.entries[t] = ((q1 as u64) << 32) | (self.entries[t] as u32 as u64);
                    }
                    radix_sort_u64(&mut self.entries[group_start..k]);
                }
                group_start = k;
            }
        }
        let counter_base = (window_length as usize) * 2;
        self.counters.clear();
        self.counters.resize(counter_base, 0);
        self.counter_entries.clear();
        self.counter_entries.extend_from_slice(&self.entries);
        self.index_a.scan = 0;
        self.index_a.sorted_count = sorted_count as i64;
        self.index_a.reduce_next = 0;
        self.index_a.skip_next = 0;
        self.index_a.anchor_start = u64::MAX;
        self.index_a.anchor_end = u64::MAX;
        self.index_b.scan = 0;
        self.index_b.sorted_count = self.counter_entries.len() as i64;
        self.index_b.reduce_next = 0;
        self.index_b.skip_next = 0;
        self.index_b.anchor_start = u64::MAX;
        self.index_b.anchor_end = u64::MAX;
        self.graph.consensus_window_length = consensus_window as u32;
        self.auxiliary.clear();
        let mut nec = 0u64;
        let mut s = 0i64;
        let mut e = (window_length as i64).min(query_length);
        let mut rr = 0i64;
        while s < query_length {
            let ids = iter_cc_idx_t(
                &*overlaps,
                &self.candidates,
                &mut self.index_a,
                &mut self.entries,
                s,
                e,
                rr != 0,
                false,
            );
            let (next_rr, added) = self.wcns_vote(
                &*overlaps,
                reads,
                query,
                exz,
                &ids,
                s,
                e,
                occurrence_total,
                occurrence_exact,
            );
            rr = next_rr;
            nec += added;
            s += window_length as i64;
            e += window_length as i64;
            e = e.min(query_length);
        }
        if self.index_b.anchor_end > self.index_b.anchor_start
            && self.index_b.anchor_start != u64::MAX
        {
            nec += self.push_cns_anchor(
                &*overlaps,
                reads,
                self.index_b.anchor_start,
                self.index_b.anchor_end,
                query,
                exz,
                false,
                occurrence_total,
                occurrence_exact,
            );
        }
        nec += self.push_cns_anchor(
            &*overlaps,
            reads,
            query_length as u64,
            query_length as u64,
            query,
            exz,
            true,
            occurrence_total,
            occurrence_exact,
        );
        nec
    }
}

/// The site predicate of `wcns_vote` and `cns_gen0`, reference read included.
fn vote_passes(oc0: u64, oc1: u64, occurrence_total: u64, occurrence_exact: f64) -> bool {
    (oc0 as f64) > (oc1 as f64) * occurrence_exact
        && oc0 > oc1.wrapping_sub(oc0)
        && oc1 >= occurrence_total
        && oc0 > 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn push_trace_bp_f_merges_equal_cells() {
        let mut cells = Vec::new();
        push_trace_bp_f(&mut cells, 0, u16::MAX, u16::MAX, 3, false);
        assert_eq!(cells, vec![3]);
        push_trace_bp_f(&mut cells, 0, u16::MAX, u16::MAX, 4, true);
        assert_eq!(cells, vec![7]);
        push_trace_bp_f(&mut cells, 1, 2, 1, 1, true);
        assert_eq!(cells, vec![7, 0x4000 | (1 << 12) | (2 << 10) | 1]);
    }

    #[test]
    fn pop_trace_bp_f_round_trips_bases() {
        let operation_two = (2 << 14) | (2 << 12);
        let cells = vec![operation_two | 5, operation_two | 3];
        let (next, operation, base_query, base_target, length) = pop_trace_bp_f_full(&cells, 0);
        assert_eq!((next, operation, base_query, base_target, length), (2, 2, u16::MAX, 2, 8));

        let deletion = (3 << 14) | (3 << 12) | 2;
        let (next, operation, base_query, base_target, length) =
            pop_trace_bp_f_full(&[deletion], 0);
        assert_eq!((next, operation, base_query, base_target, length), (1, 3, 3, u16::MAX, 2));
    }
}

/// Golden probes against the C++ reference. Ignored by default; run with
/// `cargo test --lib -- --ignored`.
#[cfg(test)]
mod probe_goldens {
    use super::*;

    /// The parameters `worker_hap_ec` passes to `wcns_gen` on the HiFi path.
    const BOUNDARY: i64 = 0;
    const WINDOW_LENGTH: u64 = 512;
    const OCCURRENCE_TOTAL: u64 = 3;
    const OCCURRENCE_EXACT: f64 = 0.500001;
    const CONSENSUS_WINDOW: u64 = 256;

    /// One window of the golden input, before it is turned into a
    /// [`WindowCigar`].
    #[derive(Debug)]
    struct GoldenWindow {
        x_start: i32,
        x_end: i32,
        y_start: i32,
        y_end: i32,
        error: i16,
        extra_begin: i16,
        extra_end: i16,
        error_threshold: i16,
        cells: Vec<(u8, u32)>,
    }

    /// One overlap of the golden input.
    #[derive(Debug)]
    struct GoldenOverlap {
        y_identifier: u32,
        strand: u32,
        x_start: u32,
        x_end: u32,
        y_start: u32,
        y_end: u32,
        alignment_length: u32,
        is_match: u8,
        windows: Vec<GoldenWindow>,
    }

    /// One complete scenario: the input the probe fed to the reference and the
    /// output it recorded.
    #[derive(Debug)]
    struct GoldenScenario {
        name: String,
        query_length: i64,
        query: Vec<u8>,
        reads: Vec<Vec<u8>>,
        overlaps: Vec<GoldenOverlap>,
        nec: u64,
        without_large_indel: Vec<u8>,
        /// The `AUXW` lines, kept verbatim for an exact comparison.
        auxiliary: Vec<String>,
    }

    /// Reads the `value` out of a `key=value` token.
    fn value<'a>(token: &'a str, key: &str) -> &'a str {
        token
            .strip_prefix(key)
            .unwrap_or_else(|| panic!("expected `{key}` in `{token}`"))
    }

    /// Parses a `key=value` token.
    fn parse_number<T: std::str::FromStr>(token: &str, key: &str) -> T {
        value(token, key)
            .parse()
            .unwrap_or_else(|_| panic!("cannot parse `{token}` as `{key}`"))
    }

    /// Parses a `<label> <number>` line.
    fn parse_labeled_number(line: &str, label: &str) -> u64 {
        let mut tokens = line.split_whitespace();
        assert_eq!(tokens.next(), Some(label), "expected `{label}` in `{line}`");
        tokens
            .next()
            .and_then(|token| token.parse().ok())
            .unwrap_or_else(|| panic!("cannot parse `{line}`"))
    }

    /// A cursor over the lines of one scenario block.
    struct Lines<'a> {
        lines: &'a [&'a str],
        at: usize,
    }

    impl<'a> Lines<'a> {
        fn next(&mut self) -> &'a str {
            let line = self.lines[self.at];
            self.at += 1;
            line
        }
    }

    fn parse_window(line: &str) -> GoldenWindow {
        let tokens: Vec<&str> = line.split_whitespace().collect();
        assert_eq!(tokens[0], "WINDOW", "expected a window line: `{line}`");
        let cell_count: usize = tokens[11]
            .parse()
            .unwrap_or_else(|_| panic!("bad cell count in `{line}`"));
        assert_eq!(
            tokens.len(),
            12 + cell_count,
            "cell count does not match the cells in `{line}`"
        );
        let mut cells = Vec::with_capacity(cell_count);
        for token in &tokens[12..] {
            let (operation, length) = token
                .split_once(':')
                .unwrap_or_else(|| panic!("bad cell `{token}` in `{line}`"));
            cells.push((operation.parse().unwrap(), length.parse().unwrap()));
        }
        GoldenWindow {
            x_start: parse_number(tokens[3], "xs="),
            x_end: parse_number(tokens[4], "xe="),
            y_start: parse_number(tokens[5], "ys="),
            y_end: parse_number(tokens[6], "ye="),
            error: parse_number(tokens[7], "err="),
            extra_begin: parse_number(tokens[8], "eb="),
            extra_end: parse_number(tokens[9], "ee="),
            error_threshold: parse_number(tokens[10], "et="),
            cells,
        }
    }

    fn parse_overlap(line: &str, expected_index: usize) -> GoldenOverlap {
        let tokens: Vec<&str> = line.split_whitespace().collect();
        assert_eq!(tokens[0], "OVERLAP", "expected an overlap line: `{line}`");
        assert_eq!(tokens[1].parse::<usize>().unwrap(), expected_index);
        GoldenOverlap {
            y_identifier: parse_number(tokens[2], "y="),
            strand: parse_number(tokens[3], "strand="),
            x_start: parse_number(tokens[4], "xps="),
            x_end: parse_number(tokens[5], "xpe="),
            y_start: parse_number(tokens[6], "yps="),
            y_end: parse_number(tokens[7], "ype="),
            alignment_length: parse_number(tokens[8], "alen="),
            is_match: parse_number(tokens[9], "match="),
            windows: Vec::new(),
        }
    }

    fn parse_scenario(lines: &[&str]) -> GoldenScenario {
        let mut input = Lines { lines, at: 0 };

        let header = input.next();
        let tokens: Vec<&str> = header.split_whitespace().collect();
        assert_eq!(tokens[0], "SCENARIO", "expected a scenario header");
        let name = tokens[1].to_string();
        let query_length: i64 = parse_number(tokens[2], "ql=");
        assert_eq!(parse_number::<i64>(tokens[3], "bd="), BOUNDARY);
        assert_eq!(parse_number::<u64>(tokens[4], "wl="), WINDOW_LENGTH);
        assert_eq!(
            parse_number::<u64>(tokens[5], "occ_tot="),
            OCCURRENCE_TOTAL
        );
        assert_eq!(
            parse_number::<f64>(tokens[6], "occ_exact="),
            OCCURRENCE_EXACT
        );
        assert_eq!(
            parse_number::<u64>(tokens[7], "cns_wl="),
            CONSENSUS_WINDOW
        );

        let query_line = input.next();
        let query = query_line
            .strip_prefix("QUERY ")
            .unwrap_or_else(|| panic!("expected a query line: `{query_line}`"))
            .as_bytes()
            .to_vec();

        let read_count = parse_labeled_number(input.next(), "NREADS") as usize;
        let mut reads = Vec::with_capacity(read_count);
        for _ in 0..read_count {
            let line = input.next();
            let tokens: Vec<&str> = line.split_whitespace().collect();
            assert_eq!(tokens[0], "READ", "expected a read line: `{line}`");
            reads.push(tokens[2].as_bytes().to_vec());
        }

        let overlap_count = parse_labeled_number(input.next(), "NOVERLAPS") as usize;
        let mut overlaps = Vec::with_capacity(overlap_count);
        for k in 0..overlap_count {
            let mut overlap = parse_overlap(input.next(), k);
            let window_header = input.next();
            let tokens: Vec<&str> = window_header.split_whitespace().collect();
            assert_eq!(tokens[0], "NWINDOWS");
            assert_eq!(tokens[1].parse::<usize>().unwrap(), k);
            let window_count: usize = tokens[2].parse().unwrap();
            for _ in 0..window_count {
                overlap.windows.push(parse_window(input.next()));
            }
            overlaps.push(overlap);
        }

        let result_line = input.next();
        assert_eq!(
            result_line.strip_prefix("RESULT "),
            Some(name.as_str()),
            "RESULT does not repeat the scenario name"
        );
        let nec = parse_labeled_number(input.next(), "NEC");
        let output_count = parse_labeled_number(input.next(), "NOUTOVERLAPS") as usize;
        let mut without_large_indel = Vec::with_capacity(output_count);
        for k in 0..output_count {
            let line = input.next();
            let tokens: Vec<&str> = line.split_whitespace().collect();
            assert_eq!(tokens[0], "OVL", "expected an overlap result: `{line}`");
            assert_eq!(tokens[1].parse::<usize>().unwrap(), k);
            without_large_indel.push(parse_number(tokens[2], "wli="));
        }
        let auxiliary_count = parse_labeled_number(input.next(), "AUX") as usize;
        let mut auxiliary = Vec::with_capacity(auxiliary_count);
        for k in 0..auxiliary_count {
            let line = input.next();
            let tokens: Vec<&str> = line.split_whitespace().collect();
            assert_eq!(tokens[0], "AUXW", "expected an auxiliary window: `{line}`");
            assert_eq!(tokens[1].parse::<usize>().unwrap(), k);
            auxiliary.push(line.to_string());
        }
        assert_eq!(input.at, lines.len(), "trailing lines in `{name}`");

        GoldenScenario {
            name,
            query_length,
            query,
            reads,
            overlaps,
            nec,
            without_large_indel,
            auxiliary,
        }
    }

    /// Splits the golden into its scenario blocks, one per `SCENARIO` line.
    fn parse_golden(text: &str) -> Vec<GoldenScenario> {
        let mut blocks: Vec<Vec<&str>> = Vec::new();
        for line in text.lines() {
            if line.starts_with("SCENARIO ") {
                blocks.push(Vec::new());
            }
            blocks
                .last_mut()
                .expect("the golden must begin with a SCENARIO line")
                .push(line);
        }
        blocks.iter().map(|block| parse_scenario(block)).collect()
    }

    /// Rebuilds the input side of one scenario.
    fn build_overlaps(scenario: &GoldenScenario) -> Vec<OverlapRegion> {
        scenario
            .overlaps
            .iter()
            .map(|overlap| {
                let mut region = OverlapRegion {
                    x_identifier: 0,
                    x_position_start: overlap.x_start,
                    x_position_end: overlap.x_end,
                    y_identifier: overlap.y_identifier,
                    y_position_start: overlap.y_start,
                    y_position_end: overlap.y_end,
                    y_strand: overlap.strand,
                    overlap_length: overlap.x_end - overlap.x_start + 1,
                    alignment_length: overlap.alignment_length,
                    is_match: overlap.is_match,
                    ..OverlapRegion::default()
                };
                for window in &overlap.windows {
                    let index = region.window_cigar.append_window(
                        window.x_start,
                        window.x_end,
                        window.y_start,
                        window.y_end,
                        window.error,
                        window.extra_begin,
                        window.extra_end,
                        window.error_threshold,
                    );
                    let start = region.window_cigar.cells.len();
                    for &(operation, length) in &window.cells {
                        region.window_cigar.push_cell(operation, length);
                    }
                    region.window_cigar.assign_cells(index, start);
                }
                region
            })
            .collect()
    }

    /// Renders the auxiliary window exactly like the probe's `AUXW` line.
    fn render_auxiliary_window(scratch: &ConsensusScratch, index: usize) -> String {
        let window = scratch.auxiliary.windows[index];
        let mut line = format!(
            "AUXW {index} xs={} xe={} clen={} cidx={}",
            window.x_start, window.x_end, window.cigar_length, window.cigar_start
        );
        line.push_str(&format!(" {}", window.cigar_length));
        for offset in 0..window.cigar_length as usize {
            let cell = scratch.auxiliary.cells[window.cigar_start as usize + offset];
            line.push_str(&format!(" {cell:04x}"));
        }
        line
    }

    #[test]
    #[ignore = "requires reference/ec_wcns_probe.out from the C++ probe"]
    fn replays_the_wcns_gen_golden() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/reference/ec_wcns_probe.out");
        let text = std::fs::read_to_string(path)
            .unwrap_or_else(|error| panic!("cannot read `{path}`: {error}"));
        let scenarios = parse_golden(&text);
        assert!(!scenarios.is_empty(), "the golden has no scenarios");

        for scenario in &scenarios {
            let mut reads = ReadStore::default();
            for (k, sequence) in scenario.reads.iter().enumerate() {
                reads.push_read(&format!("read{k}"), sequence);
            }
            let mut overlaps = build_overlaps(scenario);
            let mut scratch = ConsensusScratch::new();
            let mut exz = BitExtz::new();
            let nec = scratch.wcns_gen(
                &mut overlaps,
                &reads,
                &scenario.query,
                &mut exz,
                BOUNDARY,
                WINDOW_LENGTH,
                scenario.query_length,
                OCCURRENCE_TOTAL,
                OCCURRENCE_EXACT,
                CONSENSUS_WINDOW,
            );

            assert_eq!(nec, scenario.nec, "{}: NEC", scenario.name);
            assert_eq!(
                overlaps.len(),
                scenario.without_large_indel.len(),
                "{}: overlap count",
                scenario.name
            );
            for (k, &expected) in scenario.without_large_indel.iter().enumerate() {
                assert_eq!(
                    overlaps[k].without_large_indel, expected,
                    "{}: OVL {k} wli",
                    scenario.name
                );
            }
            assert_eq!(
                scratch.auxiliary.windows.len(),
                scenario.auxiliary.len(),
                "{}: AUX count",
                scenario.name
            );
            for (k, expected) in scenario.auxiliary.iter().enumerate() {
                assert_eq!(
                    &render_auxiliary_window(&scratch, k),
                    expected,
                    "{}: AUXW {k}",
                    scenario.name
                );
            }
        }
    }
}

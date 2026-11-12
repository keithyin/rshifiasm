//! The string graph used by the assembly stage, the counterpart of the
//! `asg_t` family in `Overlaps.h` and `Overlaps.cpp`.
//!
//! Every vertex is a read end: read `i` contributes vertices `2*i` (forward)
//! and `2*i+1` (reverse complement). An arc stores the query vertex in its
//! high 32 bits (with the overlap flank as the low 32 bits) and the target
//! vertex in `v`. All pointers of the reference become indices here.

use crate::container::ksort::radix_sort_by_key;

/// `ASG_ET_MERGEABLE` and friends, the edge types classified by `asg_extend`.
pub const ASG_ET_MERGEABLE: u8 = 0;
pub const ASG_ET_TIP: u8 = 1;
pub const ASG_ET_MULTI_OUT: u8 = 2;
pub const ASG_ET_MULTI_NEI: u8 = 3;

/// `OU_MASK`: the 14-bit `ou` field of `asg_arc_t`.
pub const OU_MASK: u32 = 0x3fff;

/// One arc, `asg_arc_t`.
///
/// `ul` packs `<query vertex>:32 | <node length>:32`; `v` is the target
/// vertex. The C structure packs `ol`/`del` and `ou`/`strong`/`no_l_indel`
/// into bit fields, which are plain fields here because every value assigned
/// already fits the field width.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct AsgArc {
    /// Query vertex in the high half, node length in the low half.
    pub ul: u64,
    /// Target vertex.
    pub v: u32,
    /// Overlap length, `ol` (31 bits in the reference).
    pub overlap: u32,
    /// `del`.
    pub deleted: bool,
    /// `ou`, a 14-bit occurrence count.
    pub occurrence: u16,
    /// `strong`.
    pub strong: bool,
    /// `no_l_indel`.
    pub without_large_indel: bool,
    /// `el`.
    pub exact_region: u8,
}

impl AsgArc {
    /// `asg_arc_len`: the node length stored in the low half of `ul`.
    pub fn node_length(&self) -> u32 {
        self.ul as u32
    }

    /// The query vertex stored in the high half of `ul`.
    pub fn query_vertex(&self) -> u32 {
        (self.ul >> 32) as u32
    }
}

/// One vertex, `asg_seq_t`: the trimmed read length and its deletion flag.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct AsgSeq {
    /// `len` (31 bits in the reference).
    pub length: u32,
    /// `del`.
    pub deleted: bool,
    /// `c`, the coverage class assigned by `ma_hit_sub`.
    pub coverage_class: u8,
}

/// One unitig of a unitig graph, `ma_utg_t`.
#[derive(Debug, Clone, Default)]
pub struct MaUtg {
    /// `len` (31 bits in the reference) and `circ`.
    pub length: u32,
    pub circular: bool,
    /// Starting vertex of the unitig in the read graph.
    pub start: u32,
    /// Ending vertex of the unitig in the read graph.
    pub end: u32,
    /// Reads contained in this unitig.
    pub reads: Vec<u64>,
    /// Unitig sequence, empty when it has not been materialised.
    pub sequence: Vec<u8>,
}

/// The string graph, `asg_t`.
#[derive(Debug, Clone, Default)]
pub struct Asg {
    pub arcs: Vec<AsgArc>,
    /// `idx[v] = start << 32 | count`; rebuilt by `arc_index`.
    pub index: Vec<u64>,
    /// Whether `arcs` is sorted by `ul`.
    pub is_sorted: bool,
    pub is_symmetric: bool,
    pub sequences: Vec<AsgSeq>,
    /// `r_seq`: the number of read vertices before any unitig vertex.
    pub read_sequence_count: u32,
    /// Vertices whose sequence has been materialised (`seq_vis`).
    pub sequence_visited: Vec<u8>,
}

impl Asg {
    /// `asg_init`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `asg_arc_n`: the number of outgoing arcs of vertex `v`.
    pub fn arc_count(&self, v: u32) -> u32 {
        self.index[v as usize] as u32
    }

    /// `asg_arc_a`: the start index and length of vertex `v`'s arc span.
    pub fn arc_range(&self, v: u32) -> (usize, usize) {
        let packed = self.index[v as usize];
        ((packed >> 32) as usize, packed as u32 as usize)
    }

    /// The outgoing arcs of vertex `v`.
    pub fn arcs_of(&self, v: u32) -> &[AsgArc] {
        let (start, count) = self.arc_range(v);
        &self.arcs[start..start + count]
    }

    /// The mutable outgoing arcs of vertex `v`.
    pub fn arcs_of_mut(&mut self, v: u32) -> &mut [AsgArc] {
        let (start, count) = self.arc_range(v);
        &mut self.arcs[start..start + count]
    }

    /// `asg_seq_set`: grows the vertex array and assigns a length.
    pub fn set_sequence(&mut self, sid: usize, length: u32, deleted: bool) {
        if sid >= self.sequences.len() {
            self.sequences.resize(sid + 1, AsgSeq::default());
        }
        let entry = &mut self.sequences[sid];
        entry.deleted = deleted;
        entry.length = length;
    }

    /// `asg_arc_pushp`: appends an arc.
    pub fn push_arc(&mut self, arc: AsgArc) {
        self.arcs.push(arc);
    }

    /// `asg_arc_sort`: sorts the arcs by `ul` with the reference radix sort.
    pub fn sort_arcs(&mut self) {
        radix_sort_by_key(&mut self.arcs, &|arc: &AsgArc| arc.ul);
    }

    /// `asg_arc_index_core`/`asg_arc_index`: rebuilds the per-vertex index.
    pub fn arc_index(&mut self) {
        let vertex_count = self.sequences.len() * 2;
        let mut index = vec![0u64; vertex_count];
        let mut last = 0usize;
        for i in 1..=self.arcs.len() {
            if i == self.arcs.len() || self.arcs[i - 1].ul >> 32 != self.arcs[i].ul >> 32 {
                index[(self.arcs[i - 1].ul >> 32) as usize] =
                    ((last as u64) << 32) | (i - last) as u64;
                last = i;
            }
        }
        self.index = index;
    }

    /// `asg_arc_rm`: drops deleted arcs and the arcs of deleted vertices.
    pub fn remove_deleted_arcs(&mut self) {
        let mut kept: Vec<AsgArc> = Vec::with_capacity(self.arcs.len());
        for arc in &self.arcs {
            let u = (arc.ul >> 32) as usize;
            let v = arc.v as usize;
            if !arc.deleted
                && !self.sequences[u >> 1].deleted
                && !self.sequences[v >> 1].deleted
            {
                kept.push(*arc);
            }
        }
        if kept.len() < self.arcs.len() {
            self.arcs = kept;
            self.index.clear();
        }
    }

    /// `asg_cleanup`: removes deleted arcs, sorts and reindexes.
    pub fn cleanup(&mut self) {
        self.remove_deleted_arcs();
        if !self.is_sorted {
            self.sort_arcs();
            self.is_sorted = true;
        }
        if self.index.is_empty() {
            self.arc_index();
        }
    }

    /// `insert_index`: the arc-array position where a new arc of `v` belongs.
    ///
    /// The reference keeps the arcs sorted by query vertex, so the insertion
    /// point is the end of `v`'s own span, or the start of the next vertex that
    /// owns arcs when `v` has none.
    pub fn insert_index(&self, v: u32) -> u32 {
        let vertex_count = self.sequences.len() as u32 * 2;
        if v >= vertex_count {
            return u32::MAX;
        }
        let count = self.arc_count(v);
        if count != 0 {
            return ((self.index[v as usize] >> 32) as u32).wrapping_add(count);
        }
        let mut w = v;
        while w < vertex_count && self.arc_count(w) == 0 {
            w += 1;
        }
        if w >= vertex_count {
            return self.arcs.len() as u32;
        }
        (self.index[w as usize] >> 32) as u32
    }

    /// `insert_arc`: inserts one arc at `index` and keeps the graph sorted.
    ///
    /// `m_distance` is always one on the paths ported here, so the tail of the
    /// arc array simply moves up by one slot. The new arc walks down through
    /// the arcs of the same query vertex until the sort order is restored; if
    /// it meets an identical live arc on the way it is marked deleted instead.
    /// The per-vertex index is updated in the same pass.
    pub fn insert_arc(
        &mut self,
        index: usize,
        src: u32,
        source_length: u32,
        dest: u32,
        overlap_length: u32,
        strong: bool,
        exact_region: u8,
        without_large_indel: bool,
    ) {
        // `insert_index_p` with `m_distance == 1`.
        self.arcs.insert(index, AsgArc::default());

        let mut arc = AsgArc {
            ul: (u64::from(src) << 32) | u64::from(source_length.wrapping_sub(overlap_length)),
            v: dest,
            // `ol` is a 31-bit bit field in the reference.
            overlap: overlap_length & 0x7fff_ffff,
            deleted: false,
            occurrence: 0,
            strong,
            without_large_indel,
            exact_region,
        };

        // The new arc bubbles down through the arcs of this query vertex while
        // it sorts before them; every displaced arc moves one slot up.
        let new_ul = arc.ul;
        let mut position = index;
        while position > 0 && new_ul <= self.arcs[position - 1].ul {
            let previous = self.arcs[position - 1];
            if !previous.deleted
                && previous.v == arc.v
                && (previous.ul >> 32) == (arc.ul >> 32)
            {
                arc.deleted = true;
                break;
            }
            self.arcs[position] = previous;
            position -= 1;
        }
        self.arcs[position] = arc;

        // `asg_arc_n` is still the pre-insertion count here.
        if self.arc_count(src) == 0 {
            self.index[src as usize] = ((index as u64) << 32) | 1;
        } else {
            self.index[src as usize] += 1;
        }

        let vertex_count = self.sequences.len() as u32 * 2;
        for v in (src + 1)..vertex_count {
            if self.arc_count(v) == 0 {
                continue;
            }
            self.index[v as usize] = self.index[v as usize].wrapping_add(1u64 << 32);
        }
    }

    /// `asg_append_edges_to_srt`: appends the forward arc `src -> dest`.
    ///
    /// The reference only links `src -> dest`; the reverse arc is commented
    /// out. Returns false, without touching the graph, when either endpoint is
    /// past the vertex array.
    pub fn append_edge_to_sorted(
        &mut self,
        src: u32,
        source_length: u32,
        dest: u32,
        overlap_length: u32,
        strong: bool,
        exact_region: u8,
        without_large_indel: bool,
    ) -> bool {
        let vertex_count = self.sequences.len() as u32 * 2;
        if src >= vertex_count || dest >= vertex_count {
            return false;
        }
        let index = self.insert_index(src) as usize;
        self.insert_arc(
            index,
            src,
            source_length,
            dest,
            overlap_length,
            strong,
            exact_region,
            without_large_indel,
        );
        true
    }

    /// `asg_seq_del`: marks a vertex deleted and deletes its arcs.
    pub fn delete_sequence(&mut self, s: u32) {
        self.sequences[s as usize].deleted = true;
        for d in 0..2u32 {
            let v = (s << 1) | d;
            let targets: Vec<u32> = self.arcs_of(v).iter().map(|arc| arc.v).collect();
            for arc in self.arcs_of_mut(v) {
                arc.deleted = true;
            }
            for target in targets {
                self.delete_arc(target ^ 1, v ^ 1, true);
            }
        }
    }

    /// `asg_arc_del`: clears or sets the deletion flag of every `v -> w` arc.
    pub fn delete_arc(&mut self, v: u32, w: u32, delete: bool) {
        let (start, count) = self.arc_range(v);
        for arc in &mut self.arcs[start..start + count] {
            if arc.v == w {
                arc.deleted = delete;
            }
        }
    }

    /// `we->del = 1`: marks the first live `v -> w` arc deleted.
    pub fn delete_first_arc(&mut self, v: u32, w: u32) -> bool {
        let (start, count) = self.arc_range(v);
        for arc in &mut self.arcs[start..start + count] {
            if !arc.deleted && arc.v == w {
                arc.deleted = true;
                return true;
            }
        }
        false
    }

    /// `asg_arc_del_multi`: removes duplicate outgoing arcs to the same target.
    pub fn delete_multi_arcs(&mut self) -> u32 {
        let vertex_count = self.sequences.len() * 2;
        let mut counts = vec![0u32; vertex_count];
        let mut removed = 0u32;
        for v in 0..vertex_count {
            let nv = self.arc_count(v as u32);
            if nv < 2 {
                continue;
            }
            let (start, count) = self.arc_range(v as u32);
            for arc in &self.arcs[start..start + count] {
                counts[arc.v as usize] += 1;
            }
            // The reference walks the span backwards, so among duplicates the
            // arc at the lowest index (shortest node length, since the span is
            // sorted by `ul`) survives.
            for i in (start..start + count).rev() {
                let target = self.arcs[i].v as usize;
                counts[target] -= 1;
                if counts[target] != 0 {
                    self.arcs[i].deleted = true;
                    removed += 1;
                }
            }
        }
        if removed > 0 {
            self.cleanup();
        }
        removed
    }

    /// `asg_arc_del_asymm`: removes arcs whose reverse complement is absent.
    pub fn delete_asymmetric_arcs(&mut self) -> u32 {
        let mut removed = 0u32;
        for e in 0..self.arcs.len() {
            let v = self.arcs[e].v ^ 1;
            let u = ((self.arcs[e].ul >> 32) as u32) ^ 1;
            let mut found = false;
            for arc in self.arcs_of(v) {
                if arc.v == u {
                    found = true;
                    break;
                }
            }
            if !found {
                self.arcs[e].deleted = true;
                removed += 1;
            }
        }
        if removed > 0 {
            self.cleanup();
        }
        removed
    }

    /// `asg_symm`: makes the graph symmetric.
    pub fn symmetrize(&mut self) {
        self.delete_multi_arcs();
        self.delete_asymmetric_arcs();
        self.is_symmetric = true;
    }

    /// `asg_get_arc`: finds the live `v -> w` arc.
    pub fn get_arc(&self, v: u32, w: u32) -> Option<AsgArc> {
        self.arcs_of(v).iter().find(|arc| !arc.deleted && arc.v == w).copied()
    }

    /// `asg_is_utg_end`: classifies the edge leaving `v`.
    ///
    /// Returns the edge type and, when mergeable, the packed
    /// `<node length>:32 | <target vertex>:32` word (`lw`).
    pub fn is_unitig_end(&self, v: u32) -> (u8, Option<u64>) {
        let (kind, packed, _) = self.is_unitig_end_with_occurrence(v);
        (kind, packed)
    }

    /// `asg_end`: [`Asg::is_unitig_end`] plus the `ou` field of the arc.
    ///
    /// The reference fills `ou` whenever the vertex has exactly one live arc,
    /// even when the neighbour check then fails, but every caller reads it
    /// only for a mergeable edge.
    pub fn is_unitig_end_with_occurrence(&self, v: u32) -> (u8, Option<u64>, u32) {
        let mut live0 = 0u32;
        let mut i0 = 0usize;
        for (i, arc) in self.arcs_of(v ^ 1).iter().enumerate() {
            if !arc.deleted {
                i0 = i;
                live0 += 1;
            }
        }
        if live0 == 0 {
            return (ASG_ET_TIP, None, 0);
        }
        if live0 > 1 {
            return (ASG_ET_MULTI_OUT, None, 0);
        }
        let arc = self.arcs_of(v ^ 1)[i0];
        let lw = (arc.ul << 32) | u64::from(arc.v);
        let w = arc.v ^ 1;
        let mut nw = 0u32;
        for out in self.arcs_of(w) {
            if !out.deleted {
                nw += 1;
            }
        }
        if nw != 1 {
            return (ASG_ET_MULTI_NEI, None, 0);
        }
        (ASG_ET_MERGEABLE, Some(lw), u32::from(arc.occurrence))
    }

    /// `asg_extend`: walks the unambiguous prefix chain from `v`.
    ///
    /// `buffer` receives the starting vertex first and then one packed
    /// `<node length>:32 | <target vertex>:32` word per merged step; the return
    /// value is the classification that stopped the walk.
    pub fn extend(&self, v: u32, max_ext: i32, buffer: &mut Vec<u64>) -> u8 {
        buffer.clear();
        let mut v = v;
        let mut max_ext = max_ext;
        let mut result;
        buffer.push(u64::from(v));
        loop {
            let (kind, lw) = self.is_unitig_end(v ^ 1);
            result = kind;
            if kind != ASG_ET_MERGEABLE {
                break;
            }
            let lw = lw.expect("a mergeable edge always reports its endpoint");
            buffer.push(lw);
            v = lw as u32;
            max_ext -= 1;
            if max_ext <= 0 {
                break;
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph_with_arcs(arcs: &[(u32, u32, u32)]) -> Asg {
        let mut graph = Asg::new();
        graph.set_sequence(0, 100, false);
        graph.set_sequence(1, 200, false);
        for &(u, v, length) in arcs {
            graph.push_arc(AsgArc {
                ul: ((u as u64) << 32) | length as u64,
                v,
                overlap: 50,
                ..AsgArc::default()
            });
        }
        graph.cleanup();
        graph
    }

    #[test]
    fn vertex_arcs_are_indexed_by_query_vertex() {
        let graph = graph_with_arcs(&[(0, 2, 100), (0, 3, 90), (1, 3, 80)]);
        assert_eq!(graph.arc_count(0), 2);
        assert_eq!(graph.arc_count(1), 1);
        // `cleanup` sorts each span by `ul`, so the shorter node (90) comes
        // first even though it was pushed second.
        assert_eq!(graph.arcs_of(0)[0].v, 3);
        assert_eq!(graph.arcs_of(0)[1].v, 2);
        assert_eq!(graph.arcs_of(3).len(), 0);
    }

    #[test]
    fn symmetric_graphs_lose_their_asymmetric_arcs() {
        let mut graph = graph_with_arcs(&[(0, 2, 100), (2, 1, 100)]);
        // 0->2 implies 3->1; 2->1 implies 0->3 (missing), so both go away.
        assert_eq!(graph.delete_asymmetric_arcs(), 2);
        assert!(graph.arcs.is_empty());
    }

    #[test]
    fn duplicate_arcs_keep_the_first_in_sorted_order() {
        // `cleanup` sorts the span by `ul`, so the first arc is the one with
        // the shorter node length; the reference keeps that one.
        let mut graph = graph_with_arcs(&[(0, 2, 100), (0, 2, 90)]);
        assert_eq!(graph.arcs_of(0)[0].node_length(), 90);
        assert_eq!(graph.delete_multi_arcs(), 1);
        assert_eq!(graph.arcs.len(), 1);
        assert_eq!(graph.arcs[0].node_length(), 90);
    }

    #[test]
    fn arc_sort_follows_the_unsigned_key() {
        let mut graph = Asg::new();
        graph.set_sequence(0, 10, false);
        graph.push_arc(AsgArc { ul: 2, ..AsgArc::default() });
        graph.push_arc(AsgArc { ul: 1, ..AsgArc::default() });
        graph.sort_arcs();
        assert_eq!(graph.arcs[0].ul, 1);
        assert_eq!(graph.arcs[1].ul, 2);
    }

    #[test]
    fn deleted_vertices_take_their_arcs_with_them() {
        let mut graph = graph_with_arcs(&[(0, 2, 100), (2, 1, 100), (3, 1, 80)]);
        graph.delete_sequence(0);
        assert!(graph.sequences[0].deleted);
        assert!(graph.arcs.iter().all(|arc| arc.deleted || arc.v != 2));
    }
}

//! Bubble-label workspace, the counterpart of `bub_label_t`,
//! `init_bub_label_t`, `reset_bub_label_t` and `destory_bub_label_t` in
//! `Overlaps.cpp`.
//!
//! The reference stores a pointer to the graph being labelled. The Rust
//! version keeps only the per-threshold scratch buffers; every function that
//! consumes them also receives the graph, so no borrow is stored.

/// `binfo_s_t`: the label of one vertex during a bubble search.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BubbleInfo {
    /// `p`, the optimal parent vertex.
    pub parent: u32,
    /// `d`, the shortest distance from the initial vertex.
    pub distance: u32,
    /// `r`, the number of remaining incoming arcs (31 bits in the reference).
    pub remaining: u32,
    /// `s`, the search state.
    pub state: bool,
    /// `c`, the accumulated positive-read count of the primary-trio pop.
    pub positive_count: u32,
    /// `m`, the accumulated negative-read count.
    pub negative_count: u32,
    /// `np`, the accumulated non-positive-read count.
    pub non_positive_count: u32,
    /// `nc`, the accumulated read count regardless of side.
    pub node_count: u32,
}

/// `buf_s_t`: one threshold's scratch buffers.
#[derive(Debug, Clone, Default)]
pub struct BubbleBuffer {
    /// `a`, one label per vertex (two per read).
    pub info: Vec<BubbleInfo>,
    /// `S`, the vertices whose incoming arcs have all been visited.
    pub sources: Vec<u32>,
    /// `b`, the visited vertices.
    pub visited: Vec<u32>,
    /// `e`, the visited edges.
    pub visited_edges: Vec<u32>,
}

/// `bub_label_t`: the label workspace shared by the cleaning passes.
#[derive(Debug, Clone, Default)]
pub struct BubbleLabel {
    /// `b`, `n_thres` independent label arrays.
    pub buffers: Vec<BubbleBuffer>,
    /// `n_thres`.
    pub threshold_count: u32,
    /// `n_reads`, the number of vertices (twice the sequence count).
    pub read_count: u32,
    /// `check_cross`.
    pub check_cross: u32,
    /// `bub_dist`.
    pub bubble_distance: u64,
}

impl BubbleLabel {
    /// `init_bub_label_t`.
    pub fn new(threshold_count: u32, read_count: u32) -> Self {
        let buffers = (0..threshold_count)
            .map(|_| BubbleBuffer {
                info: vec![BubbleInfo::default(); (read_count as usize) << 1],
                ..BubbleBuffer::default()
            })
            .collect();

        Self {
            buffers,
            threshold_count,
            read_count,
            check_cross: 0,
            bubble_distance: 0,
        }
    }

    /// `reset_bub_label_t`.
    ///
    /// `sequence_count` is the graph's `n_seq`, which may have grown past the
    /// read count once unitig vertices exist.
    pub fn reset(&mut self, sequence_count: usize, bubble_distance: u64, check_cross: u32) {
        self.bubble_distance = bubble_distance;
        self.check_cross = check_cross;

        if (self.read_count as usize) < sequence_count {
            self.read_count = sequence_count as u32;
            for buffer in &mut self.buffers {
                buffer.info.resize(sequence_count << 1, BubbleInfo::default());
            }
        }

        for buffer in &mut self.buffers {
            buffer.sources.clear();
            buffer.visited.clear();
            buffer.visited_edges.clear();
            buffer.info.fill(BubbleInfo::default());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn init_allocates_two_labels_per_read() {
        let labels = BubbleLabel::new(3, 5);
        assert_eq!(labels.threshold_count, 3);
        assert_eq!(labels.read_count, 5);
        assert_eq!(labels.buffers.len(), 3);
        assert_eq!(labels.buffers[0].info.len(), 10);
        assert_eq!(labels.check_cross, 0);
        assert_eq!(labels.bubble_distance, 0);
    }

    #[test]
    fn reset_grows_the_label_arrays_and_clears_the_buffers() {
        let mut labels = BubbleLabel::new(2, 1);
        labels.buffers[0].sources.push(7);
        labels.buffers[0].info[0].parent = 9;
        labels.reset(4, 75000, 1);
        assert_eq!(labels.read_count, 4);
        assert_eq!(labels.buffers[0].info.len(), 8);
        assert_eq!(labels.buffers[0].info[0], BubbleInfo::default());
        assert!(labels.buffers[0].sources.is_empty());
        assert_eq!(labels.bubble_distance, 75000);
        assert_eq!(labels.check_cross, 1);
    }
}

//! The bubble-chain machinery of `hic.cpp`, shared by the trio bubble
//! rescues (`rescue_bubble_by_chain`) and the Hi-C phasing front end.
//!
//! Everything here is a literal transcription of the reference; the `ug`
//! pointer of `bubble_type` is dropped (the Rust port passes the unitig
//! graph to the functions that need it) and the Hi-C-only weight helpers
//! (`hc_links` users) are not ported.

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::hit::{CoverageCut, ReadToUnitig};
use crate::overlap::paf::PafHit;
use crate::correction::ec_driver::PafList;
use crate::assembly::bubble::{BubbleBuffer, BubbleLabel};
use crate::assembly::gfa_clean::get_real_length;
use crate::assembly::rescue::{lable_all_bubbles, pre_clean, RescueParameters};
use crate::assembly::unitig::get_edge_from_source;
use crate::assembly::bubble_identify::{asg_bub_pop1_primary_trio, get_s_bub_pop_max_dist_advance};
use crate::container::ksort::radix_sort_by_key_32;
use crate::assembly::hap_cov::TransChain;

/// `C_HET` (`Overlaps.h` 962).
pub const C_HET: u8 = 1;
/// `P_HET` (`Overlaps.h` 963).
pub const P_HET: u8 = 2;

/// `chain_w_type` (`hic.h` 30-35): one entry of the chain weight list. The
/// trio path never fills it (the weight helpers are Hi-C-only); the field
/// set is kept for the structure's sake.
#[derive(Debug, Clone, Copy, Default)]
pub struct ChainWeightEntry {
    /// `g_occ`.
    pub graph_occurrence: i64,
    /// `b_occ`.
    pub bubble_occurrence: i64,
    /// `id`.
    pub id: u64,
    /// `del`.
    pub deleted: u8,
}

/// `bubble_type` (`hic.h` 37-53).
#[derive(Debug, Clone, Default)]
pub struct BubbleType {
    /// `index`, one classification word per unitig (two per read vertex).
    pub index: Vec<u32>,
    /// `round_id`.
    pub round_id: u32,
    /// `n_round`.
    pub n_round: u32,
    /// `list`, the flattened `[beg, sink, body...]` of every bubble.
    pub list: Vec<u32>,
    /// `num`, the per-bubble offset into `list` (`f_bub`+1 entries, the last
    /// being `list.len()`).
    pub num: Vec<u32>,
    /// `pathLen`, one entry per bubble.
    pub path_len: Vec<u64>,
    /// `b_s_idx`, per unitig: the (up to two) bubble ids whose beg/sink it
    /// is, packed two per word (`0xFFFFFFFF` in a half means empty).
    pub b_s_idx: Vec<u64>,
    /// `s_bub`.
    pub simple_bubble_count: u64,
    /// `f_bub`.
    pub full_bubble_count: u64,
    /// `b_bub`.
    pub broken_bubble_count: u64,
    /// `b_end_bub`.
    pub broken_end_bubble_count: u64,
    /// `tangle_bub`.
    pub tangle_bubble_count: u64,
    /// `cross_bub`.
    pub cross_bubble_count: u64,
    /// `mess_bub`.
    pub mess_bubble_count: u64,
    /// `check_het`.
    pub check_het: u32,
    /// `b_g`, the bubble graph built by `get_bub_graph`.
    pub bubble_graph: Option<Asg>,
    /// `b_ug`, the bubble-chain unitig graph.
    pub bubble_unitigs: Option<Vec<MaUtg>>,
    /// `b_ug->g`, the compact graph of the bubble chains (vertices are
    /// chains).
    pub bubble_chain_graph: Option<Asg>,
    /// `chain_weight` (empty on the trio path).
    pub chain_weight: Vec<ChainWeightEntry>,
}

impl BubbleType {
    /// `M_het(B)` (`hic.h` 83).
    // NOTE: the C names are counterintuitive — `P_het(B)` is `num.n` and
    // `M_het(B)` is `num.n + 1` (`hic.h` 82-83), the opposite of what the
    // letters suggest.
    pub fn m_het(&self) -> u32 {
        self.num.len() as u32 + 1
    }

    /// `P_het(B)` (`hic.h` 82).
    pub fn p_het(&self) -> u32 {
        self.num.len() as u32
    }

    /// `IF_HET(ID, B)` (`hic.h` 88).
    pub fn is_het(&self, id: u32) -> bool {
        self.index[id as usize] == self.full_bubble_count as u32 + 1
    }
}

/// `pdq` (`hic.h` 66-78): the Dijkstra workspace of `get_shortest_path`.
#[derive(Debug, Clone, Default)]
pub struct Pdq {
    /// `vis`.
    pub visited: Vec<u8>,
    /// `x`, the heap (the reference pushes a leading `0`).
    pub heap: Vec<u64>,
    /// `dis`.
    pub distance: Vec<u64>,
    /// `uID_mode`.
    pub id_mode: u64,
    /// `uID_shift`.
    pub id_shift: u64,
    /// `tmp_v`.
    pub tmp_v: u64,
    /// `tmp_d`.
    pub tmp_d: u64,
}

impl Pdq {
    /// `init_pdq` (`hic.cpp` 3387).
    pub fn new(unitig_number: u64) -> Self {
        let mut queue = Pdq {
            visited: vec![0; unitig_number as usize],
            distance: vec![0; unitig_number as usize],
            heap: vec![0],
            id_mode: 0,
            id_shift: 0,
            tmp_v: 0,
            tmp_d: 0,
        };
        let mut shift = 1u64;
        while (1u64 << shift) < unitig_number {
            shift += 1;
        }
        queue.id_shift = shift;
        // The reference computes `((uint64_t)-1) >> (64-i)`, i.e. the low
        // `shift` bits — `u64::MAX >> shift` would keep the distance bits in
        // the mask and corrupt every unpacked vertex.
        queue.id_mode = u64::MAX >> (64 - shift);
        queue
    }
}

/// `get_bubbles` (`hic.cpp` 2217): the `(beg, sink, body, path_base)` of
/// bubble `id`; the body excludes beg and sink.
///
/// The reference only reads `num[id+1]` when the caller requests the body;
/// a bubble freshly pushed by `update_bubble_graph` has no sentinel entry
/// yet, and the body length the reference computes there is uninitialized
/// memory that no caller consumes. The port falls back to the list end.
pub fn get_bubbles<'a>(bub: &'a BubbleType, id: u64) -> (u32, u32, &'a [u32], u64) {
    let entry = bub.num[id as usize] as usize;
    let next = if id as usize + 1 < bub.num.len() {
        bub.num[id as usize + 1] as usize
    } else {
        bub.list.len()
    };
    let beg = bub.list[entry];
    let sink = bub.list[entry + 1];
    let body = &bub.list[entry + 2..next];
    let path_base = bub.path_len[id as usize];
    (beg, sink, body, path_base)
}

/// `dfs_bubble` (`hic.cpp` 2286): collects every read reachable from `v`
/// (both strands) without passing through the reads `beg`/`sink`.
pub fn dfs_bubble(
    graph: &Asg,
    stack: &mut Vec<u32>,
    result: &mut Vec<u32>,
    v: u32,
    beg: u32,
    sink: u32,
) {
    stack.clear();
    result.clear();
    let start = v << 1;
    stack.push(start);
    while let Some(cur) = stack.pop() {
        if (start >> 1) != (cur >> 1) {
            result.push(cur >> 1);
        }
        for arc in graph.arcs_of(cur) {
            if arc.deleted {
                continue;
            }
            if (arc.v >> 1) == beg || (arc.v >> 1) == sink {
                continue;
            }
            stack.push(arc.v);
        }
    }

    stack.push(start + 1);
    while let Some(cur) = stack.pop() {
        if (start >> 1) != (cur >> 1) {
            result.push(cur >> 1);
        }
        for arc in graph.arcs_of(cur) {
            if arc.deleted {
                continue;
            }
            if (arc.v >> 1) == beg || (arc.v >> 1) == sink {
                continue;
            }
            stack.push(arc.v);
        }
    }
}

/// `dfs_bubble_broken` (`hic.cpp` 2232): like [`dfs_bubble`] but with a
/// visited flag and without collecting the start read.
pub fn dfs_bubble_broken(
    graph: &Asg,
    stack: &mut Vec<u32>,
    result: &mut Vec<u32>,
    visited: &mut [u8],
    v_d: u32,
    beg_d: u32,
    sink_d: u32,
) {
    visited.iter_mut().for_each(|flag| *flag = 0);
    stack.clear();
    result.clear();
    let mut v = v_d;
    if v != (beg_d ^ 1) && v != (sink_d ^ 1) {
        stack.push(v);
    }
    while let Some(cur) = stack.pop() {
        visited[cur as usize] = 1;
        if (v >> 1) != (cur >> 1) {
            result.push(cur >> 1);
        }
        for arc in graph.arcs_of(cur) {
            if arc.deleted || visited[arc.v as usize] != 0 {
                continue;
            }
            if (arc.v >> 1) == (beg_d >> 1) {
                continue;
            }
            if (arc.v >> 1) == (sink_d >> 1) {
                continue;
            }
            stack.push(arc.v);
        }
    }

    visited.iter_mut().for_each(|flag| *flag = 0);
    v ^= 1;
    if v != (beg_d ^ 1) && v != (sink_d ^ 1) {
        stack.push(v);
    }
    while let Some(cur) = stack.pop() {
        visited[cur as usize] = 1;
        if (v >> 1) != (cur >> 1) {
            result.push(cur >> 1);
        }
        for arc in graph.arcs_of(cur) {
            if arc.deleted || visited[arc.v as usize] != 0 {
                continue;
            }
            if (arc.v >> 1) == (beg_d >> 1) {
                continue;
            }
            if (arc.v >> 1) == (sink_d >> 1) {
                continue;
            }
            stack.push(arc.v);
        }
    }
}

/// `get_unitig_het_arb` (`hic.cpp` 2327): classifies a unitig by the
/// heterozygosity flags of its reads.
pub fn get_unitig_het_arb(
    unitigs: &[MaUtg],
    uid: u32,
    read_het_flag: &[u8],
    transitions: Option<&TransChain>,
    m_het_occ: u32,
    m_het_label: u32,
    p_het_label: u32,
    n_het_label: u32,
) -> u32 {
    if let Some(trans) = transitions {
        // `ref->idx.n > 0 && u_trans_n(*ref, uid) > 0`.
        if !trans.k_trans_idx.is_empty()
            && trans.k_trans_idx[uid as usize] & 0xffffffff != 0
        {
            return m_het_label;
        }
    }
    let unitig = &unitigs[uid as usize];
    let mut het_occurrence = 0u32;
    let mut hom_occurrence = 0u32;
    for read in &unitig.reads {
        let read_id = (read >> 33) as usize;
        if read_het_flag[read_id] & P_HET != 0 {
            return m_het_label;
        }
        if read_het_flag[read_id] & C_HET != 0 || read_het_flag[read_id] & P_HET != 0 {
            het_occurrence += 1;
        } else {
            hom_occurrence += 1;
        }
    }

    if het_occurrence + hom_occurrence == 0 {
        return n_het_label;
    }
    if (het_occurrence as f64) > ((het_occurrence + hom_occurrence) as f64) * 0.8
        && het_occurrence + hom_occurrence > m_het_occ
    {
        return m_het_label;
    }
    if het_occurrence >= hom_occurrence {
        return p_het_label;
    }
    n_het_label
}

/// `asg_bub_backtrack_primary_length` (`Overlaps.cpp` 24135): walks the
/// parent chain of a found bubble from its sink back to `v0`, accumulating
/// the path length and node count.
pub fn asg_bub_backtrack_primary_length(
    graph: &Asg,
    unitigs: Option<&[MaUtg]>,
    v0: u32,
    buffer: &BubbleBuffer,
    path_base_len: &mut u64,
    path_nodes: &mut u64,
) {
    let sink = buffer.sources[0];
    let mut len = 0u64;
    let mut node = 0u64;
    let mut v = sink;
    loop {
        let u = buffer.info[v as usize].parent;
        if u == v0 {
            break;
        }
        if v == sink {
            len += u64::from(graph.sequences[(u >> 1) as usize].length);
        } else {
            let mut arc_length = 0u32;
            for arc in graph.arcs_of(u) {
                if arc.deleted {
                    continue;
                }
                if arc.v == v {
                    arc_length = (arc.ul & 0xffffffff) as u32;
                    break;
                }
            }
            len += u64::from(arc_length);
        }

        match unitigs {
            Some(unitigs) => node += unitigs[(u >> 1) as usize].reads.len() as u64,
            None => node += 1,
        }

        v = u;
    }
    *path_base_len = len;
    *path_nodes = node;
}

/// `update_bub_b_s_idx` (`hic.cpp` 7663): rebuilds the beg/sink → bubble
/// index from the current `num`/`list` (two bubble ids per word).
pub fn update_bub_b_s_idx(bub: &mut BubbleType) {
    bub.b_s_idx.iter_mut().for_each(|entry| *entry = u64::MAX);
    let n_bub = bub.num.len() - 1;
    for i in 0..n_bub {
        let (beg, sink, _, _) = get_bubbles(bub, i as u64);

        if beg != u32::MAX {
            let v = (beg >> 1) as usize;
            if bub.b_s_idx[v] == u64::MAX {
                bub.b_s_idx[v] <<= 32;
                bub.b_s_idx[v] |= i as u64;
            } else if (bub.b_s_idx[v] & 0xffffffff00000000) == 0xffffffff00000000 {
                bub.b_s_idx[v] <<= 32;
                bub.b_s_idx[v] |= i as u64;
            }
        }

        if sink != u32::MAX {
            let v = (sink >> 1) as usize;
            if bub.b_s_idx[v] == u64::MAX {
                bub.b_s_idx[v] <<= 32;
                bub.b_s_idx[v] |= i as u64;
            } else if (bub.b_s_idx[v] & 0xffffffff00000000) == 0xffffffff00000000 {
                bub.b_s_idx[v] <<= 32;
                bub.b_s_idx[v] |= i as u64;
            }
        }
    }
}

/// `identify_bubbles` (`hic.cpp` 2354): finds every bubble of the unitig
/// graph, classifies the unitigs by heterozygosity and builds the bubble
/// graph and its chains.
///
/// `read_het_flag` is the reference's `r_het_flag` (`cov->t_ch->ir_het`,
/// aliased with `cov->is_r_het`); `transitions` is the `kv_u_trans_t* ref`
/// the reference passes as NULL on every trio/recovery path.
pub fn identify_bubbles(
    graph: &mut Asg,
    unitigs: &[MaUtg],
    bub: &mut BubbleType,
    read_het_flag: &[u8],
    transitions: Option<&TransChain>,
) {
    graph.cleanup();
    if !graph.is_symmetric {
        graph.symmetrize();
    }
    let n_vtx = (graph.sequences.len() * 2) as u32;
    let mode: u32 = u32::MAX << 2;
    bub.broken_bubble_count = 0;
    bub.broken_end_bubble_count = 0;
    bub.tangle_bubble_count = 0;
    bub.cross_bubble_count = 0;
    bub.mess_bubble_count = 0;

    if bub.round_id == 0 {
        let mut buffer = BubbleBuffer {
            info: vec![Default::default(); n_vtx as usize],
            ..Default::default()
        };
        let t_len = get_s_bub_pop_max_dist_advance(graph, &mut buffer);
        bub.list = Vec::new();
        bub.num = Vec::new();
        bub.path_len = Vec::new();
        bub.b_s_idx = vec![u64::MAX; graph.sequences.len()];
        bub.bubble_unitigs = None;
        bub.chain_weight = Vec::new();
        bub.index = vec![0; n_vtx as usize];
        for sequence in graph.sequences.iter_mut() {
            sequence.coverage_class = 0;
        }

        // First pass: mark the interior, beg and sink of every bubble.
        for v in 0..n_vtx {
            if graph.sequences[(v >> 1) as usize].deleted {
                continue;
            }
            if graph.arc_count(v) < 2 {
                continue;
            }
            if bub.index[v as usize] & 3 != 0 {
                continue;
            }
            if asg_bub_pop1_primary_trio(
                graph,
                None,
                v,
                t_len,
                &mut buffer,
                u32::MAX,
                u32::MAX,
                false,
                None,
                &[],
                false,
                false,
            ) != 0
            {
                let sink = buffer.sources[0];
                for &node in &buffer.visited {
                    if node == v || node == sink {
                        continue;
                    }
                    bub.index[node as usize] &= mode;
                    bub.index[node as usize] += 1;
                    bub.index[(node ^ 1) as usize] &= mode;
                    bub.index[(node ^ 1) as usize] += 1;
                }
                bub.index[v as usize] &= mode;
                bub.index[v as usize] += 2;
                bub.index[(sink ^ 1) as usize] &= mode;
                bub.index[(sink ^ 1) as usize] += 3;
            }
        }

        // Second pass: collect the bubbles rooted at every marked beg.
        let mut stack: Vec<u32> = Vec::new();
        let mut result: Vec<u32> = Vec::new();
        for v in 0..n_vtx {
            if bub.index[v as usize] & 3 != 2 {
                continue;
            }
            if asg_bub_pop1_primary_trio(
                graph,
                Some(unitigs),
                v,
                t_len,
                &mut buffer,
                u32::MAX,
                u32::MAX,
                false,
                None,
                &[],
                false,
                false,
            ) != 0
            {
                let sink = buffer.sources[0];
                // `i = b.b.n + 1` and the validation loop resets `i = 0`,
                // walking every body node; acceptance is `i == b.b.n` — the
                // loop completed without break. Bodies outside {2,3,5} never
                // run the loop and stay rejected (`i = b.b.n + 1`).
                let body_len = buffer.visited.len() as u32;
                let mut accepted = false;
                if body_len == 2 || body_len == 3 || body_len == 5 {
                    let mut completed = true;
                    for &node in &buffer.visited {
                        if node == v || node == sink {
                            continue;
                        }
                        dfs_bubble(graph, &mut stack, &mut result, node >> 1, v >> 1, sink >> 1);
                        if (result.len() as u32) + 3 != body_len
                            && (result.len() as u32) + 2 != body_len
                        {
                            completed = false;
                            break;
                        }
                    }
                    accepted = completed;
                }

                if accepted {
                    bub.num.push(v);
                } else {
                    bub.num.push(v | (1 << 31));
                }
            }
        }
        drop((stack, result));

        radix_sort_by_key_32(&mut bub.num, &|entry: &u32| *entry);
        bub.simple_bubble_count = 0;
        for k in 0..bub.num.len() {
            if bub.num[k] >> 31 == 0 {
                bub.simple_bubble_count += 1;
            }
            let v = (bub.num[k] << 1) >> 1;
            bub.num[k] = bub.list.len() as u32;
            if asg_bub_pop1_primary_trio(
                graph,
                Some(unitigs),
                v,
                t_len,
                &mut buffer,
                u32::MAX,
                u32::MAX,
                false,
                None,
                &[],
                false,
                false,
            ) != 0
            {
                let sink = buffer.sources[0];
                let mut path_len = 0u64;
                let mut path_nodes = 0u64;
                asg_bub_backtrack_primary_length(
                    graph,
                    Some(unitigs),
                    v,
                    &buffer,
                    &mut path_len,
                    &mut path_nodes,
                );
                bub.path_len.push(path_len);
                bub.list.push(v);
                bub.list.push(sink ^ 1);
                for &node in &buffer.visited {
                    if node == v || node == sink {
                        continue;
                    }
                    bub.list.push(node);
                }
            }
        }
        bub.num.push(bub.list.len() as u32);
        bub.full_bubble_count = bub.num.len() as u64 - 1;

        for i in 0..graph.sequences.len() as u32 {
            bub.index[i as usize] = get_unitig_het_arb(
                unitigs,
                i,
                read_het_flag,
                transitions,
                20,
                bub.m_het(),
                bub.p_het(),
                u32::MAX,
            );
        }
        for i in 0..bub.full_bubble_count as usize {
            let (beg, sink, body, _) = get_bubbles(bub, i as u64);
            let body: Vec<u32> = body.to_vec();
            let mut n_occ = 0u32;
            for &node in &body {
                bub.index[(node >> 1) as usize] = i as u32;
                n_occ += unitigs[(node >> 1) as usize].reads.len() as u32;
            }

            if n_occ > 3 {
                if bub.index[(beg >> 1) as usize] != bub.m_het() {
                    bub.index[(beg >> 1) as usize] = u32::MAX;
                }
                if bub.index[(sink >> 1) as usize] != bub.m_het() {
                    bub.index[(sink >> 1) as usize] = u32::MAX;
                }
            }

            let v = (beg >> 1) as usize;
            if bub.b_s_idx[v] == u64::MAX {
                bub.b_s_idx[v] <<= 32;
                bub.b_s_idx[v] |= i as u64;
            } else if (bub.b_s_idx[v] & 0xffffffff00000000) == 0xffffffff00000000 {
                bub.b_s_idx[v] <<= 32;
                bub.b_s_idx[v] |= i as u64;
            }

            let v = (sink >> 1) as usize;
            if bub.b_s_idx[v] == u64::MAX {
                bub.b_s_idx[v] <<= 32;
                bub.b_s_idx[v] |= i as u64;
            } else if (bub.b_s_idx[v] & 0xffffffff00000000) == 0xffffffff00000000 {
                bub.b_s_idx[v] <<= 32;
                bub.b_s_idx[v] |= i as u64;
            }
        }
        for i in 0..graph.sequences.len() as u32 {
            if bub.index[i as usize] == bub.m_het() {
                bub.index[i as usize] = bub.p_het();
            }
        }
    } else {
        bub.num.truncate(bub.full_bubble_count as usize + 1);
        bub.path_len.truncate(bub.full_bubble_count as usize);
        bub.list.truncate(bub.num[bub.full_bubble_count as usize] as usize);
        update_bub_b_s_idx(bub);
        bub.check_het = 0;
        bub.bubble_graph = None;
        bub.bubble_unitigs = None;
        bub.bubble_chain_graph = None;
        bub.chain_weight = Vec::new();
    }
    bub.bubble_graph = None;
    bub.bubble_unitigs = None;
    bub.bubble_chain_graph = None;
    build_bub_graph(graph, unitigs, bub);
}

impl BubbleType {
    /// `IF_HOM(ID, B)` (`hic.h` 89).
    pub fn is_hom(&self, id: u32) -> bool {
        self.index[id as usize] > self.full_bubble_count as u32 + 1
    }

    /// `get_bub_id` (`hic.cpp` 6924): the (up to two) bubble ids a unitig
    /// connects, honouring `check_het` and normalising the empty half.
    pub fn get_bub_id(&self, root: u32, check_het: u32) -> (u64, u64) {
        let mut id0 = u64::MAX;
        let mut id1 = u64::MAX;
        let word = self.b_s_idx[root as usize];

        if word & 0xffffffff != 0xffffffff {
            id0 = word & 0xffffffff;
            if check_het != 0 {
                let (beg, sink, _, _) = get_bubbles(self, id0);
                if (beg == u32::MAX || self.is_het(beg >> 1))
                    && (sink == u32::MAX || self.is_het(sink >> 1))
                {
                    id0 = u64::MAX;
                }
            }
        }

        if word & 0xffffffff00000000 != 0xffffffff00000000 {
            id1 = (word & 0xffffffff00000000) >> 32;
            if check_het != 0 {
                let (beg, sink, _, _) = get_bubbles(self, id1);
                if (beg == u32::MAX || self.is_het(beg >> 1))
                    && (sink == u32::MAX || self.is_het(sink >> 1))
                {
                    id1 = u64::MAX;
                }
            }
        }

        if id0 == u64::MAX && id1 != u64::MAX {
            id0 = id1;
            id1 = u64::MAX;
        }
        (id0, id1)
    }

    /// `connect_bub_occ` (`hic.cpp` 6962).
    pub fn connect_bub_occ(&self, root_id: u32, check_het: u32) -> u32 {
        let (id0, id1) = self.get_bub_id(root_id, check_het);
        let mut occurrence = 2u32;
        if id0 == u64::MAX {
            occurrence -= 1;
        }
        if id1 == u64::MAX {
            occurrence -= 1;
        }
        occurrence
    }
}

impl Pdq {
    fn weight(&self, index: usize) -> (u64, u64) {
        let entry = self.heap[index];
        (entry & self.id_mode, entry >> self.id_shift)
    }

    fn set_entry(v: u64, dis: u64, v_mode: u64, v_shift: u64) -> u64 {
        (dis << v_shift) | (v & v_mode)
    }

    /// `reset_pdq` (`hic.cpp` 3406).
    pub fn reset(&mut self) {
        self.heap.truncate(1);
        self.heap[0] = 0;
        self.distance.iter_mut().for_each(|entry| *entry = u64::MAX);
        self.visited.iter_mut().for_each(|flag| *flag = 0);
    }

    /// `push_pdq` (`hic.cpp` 3436).
    pub fn push(&mut self, v: u64, dis: u64) {
        self.heap
            .push(Self::set_entry(v, dis, self.id_mode, self.id_shift));
        self.heap[0] += 1;
        let mut child = self.heap[0] as usize;
        let mut parent = child >> 1;
        while parent > 0 && self.weight(child).1 < self.weight(parent).1 {
            self.heap.swap(child, parent);
            child = parent;
            parent = child >> 1;
        }
    }

    /// `pop_pdq` (`hic.cpp` 3450).
    pub fn pop(&mut self) -> Option<(u64, u64)> {
        if self.heap[0] == 0 {
            return None;
        }
        let (min_v, min_dis) = self.weight(1);
        let count = self.heap[0] as usize;
        self.heap.swap(1, count);
        self.heap[0] -= 1;
        self.heap.truncate(count);
        let count = self.heap[0] as usize;

        let mut current = 1usize;
        loop {
            let left = current << 1;
            let right = left + 1;
            if left > self.heap[0] as usize {
                break;
            }
            let min_child = if right > self.heap[0] as usize {
                left
            } else if self.weight(left).1 < self.weight(right).1 {
                left
            } else {
                right
            };
            if self.weight(current).1 > self.weight(min_child).1 {
                self.heap.swap(current, min_child);
                current = min_child;
            } else {
                break;
            }
        }
        Some((min_v, min_dis))
    }

    /// `get_shortest_path` (`hic.cpp` 3495): Dijkstra over the read graph
    /// with arc node-lengths as weights.
    pub fn shortest_path(&mut self, src: u32, graph: &Asg, pre: &mut [u32]) {
        self.reset();
        self.distance[src as usize] = 0;
        pre[src as usize] = u32::MAX;
        self.push(src as u64, 0);
        while let Some((v, _)) = self.pop() {
            self.visited[v as usize] = 1;
            for arc in graph.arcs_of(v as u32) {
                if arc.deleted {
                    continue;
                }
                let u = arc.v as u64;
                let w = (arc.ul & 0xffffffff) as u64;
                if self.visited[u as usize] == 0
                    && self.distance[u as usize] > self.distance[v as usize] + w
                {
                    self.distance[u as usize] = self.distance[v as usize] + w;
                    self.push(u, self.distance[u as usize]);
                    pre[u as usize] = v as u32;
                }
            }
        }
    }
}


/// `check_bubble_tangle` (`hic.cpp` 8092): flood-fills from `beg` on both
/// strands of `graph` and decides whether the region up to `sink` is a
/// tangle. `bub`/`chain_unitigs` enable the occurrence-based rejection.
#[allow(clippy::too_many_arguments)]
fn check_bubble_tangle(
    bub: Option<&BubbleType>,
    graph: &Asg,
    chain_unitigs: Option<&[MaUtg]>,
    unitigs: Option<&[MaUtg]>,
    beg: u32,
    sink: u32,
    side_rate: f64,
    total_rate: f64,
    beg_occ: u32,
    sink_occ: u32,
    is_vis: &mut [u8],
    stack: &mut Vec<u32>,
    res: Option<&mut Vec<u32>>,
    chain_flag: Option<&[u8]>,
    mut extra_check: Option<&mut u32>,
) -> bool {
    if let Some(extra) = extra_check.as_deref_mut() {
        *extra = 1;
    }
    let mut tan_occ = 0u32;
    let mut no_first = false;
    is_vis.iter_mut().for_each(|flag| *flag = 0);
    stack.clear();
    stack.push(beg);
    let mut res = res;
    if let Some(res) = res.as_deref_mut() {
        res.clear();
    }
    while let Some(cur) = stack.pop() {
        if no_first && cur == beg {
            return false;
        }
        if sink != u32::MAX && cur == sink {
            return false;
        }

        if is_vis[cur as usize] == 0 && is_vis[(cur ^ 1) as usize] == 0 {
            if (cur >> 1) != (beg >> 1) && (sink == u32::MAX || (cur >> 1) != (sink >> 1)) {
                if let Some(res) = res.as_deref_mut() {
                    res.push(cur);
                }
                if let Some(flags) = chain_flag {
                    if flags[(cur >> 1) as usize] != 0 {
                        if let Some(extra) = extra_check.as_deref_mut() {
                            *extra = 0;
                        }
                    }
                }

                if let (Some(bub), Some(chain_unitigs)) = (bub, chain_unitigs) {
                    tan_occ += get_btg_occ(bub, chain_unitigs, unitigs.unwrap_or(&[]), cur >> 1);
                    if (tan_occ as f64) * side_rate >= beg_occ as f64 {
                        return false;
                    }
                    if sink != u32::MAX && (tan_occ as f64) * side_rate >= sink_occ as f64 {
                        return false;
                    }
                    if (tan_occ as f64) * total_rate
                        >= (beg_occ + if sink != u32::MAX { sink_occ } else { 0 }) as f64
                    {
                        return false;
                    }
                }
            }
        }

        if is_vis[cur as usize] == 0 && cur != beg ^ 1 && (sink == u32::MAX || cur != sink ^ 1) {
            is_vis[cur as usize] = 1;
            for arc in graph.arcs_of(cur) {
                if arc.deleted {
                    continue;
                }
                if arc.v == beg {
                    return false;
                }
                if sink != u32::MAX && arc.v == sink {
                    return false;
                }
                if is_vis[arc.v as usize] != 0 {
                    continue;
                }
                stack.push(arc.v);
            }
        }
        is_vis[cur as usize] = 1;

        let cur = cur ^ 1;
        if is_vis[cur as usize] == 0 && cur != beg ^ 1 && (sink == u32::MAX || cur != sink ^ 1) {
            is_vis[cur as usize] = 1;
            for arc in graph.arcs_of(cur) {
                if arc.deleted {
                    continue;
                }
                if arc.v == beg {
                    return false;
                }
                if sink != u32::MAX && arc.v == sink {
                    return false;
                }
                if is_vis[arc.v as usize] != 0 {
                    continue;
                }
                stack.push(arc.v);
            }
        }
        is_vis[cur as usize] = 1;
        no_first = true;
    }

    if let (Some(bub), Some(chain_unitigs)) = (bub, chain_unitigs) {
        if (tan_occ as f64) * side_rate >= beg_occ as f64 {
            return false;
        }
        if sink != u32::MAX && (tan_occ as f64) * side_rate >= sink_occ as f64 {
            return false;
        }
        if (tan_occ as f64) * total_rate
            >= (beg_occ + if sink != u32::MAX { sink_occ } else { 0 }) as f64
        {
            return false;
        }
    }

    true
}

/// `get_utg_end_from_btg` (`hic.cpp` 8363).
fn get_utg_end_from_btg(chain_unitigs: &[MaUtg], v: u32) -> u32 {
    let chain = &chain_unitigs[(v >> 1) as usize];
    if v & 1 == 1 {
        (chain.reads[0] >> 32) as u32 ^ 1
    } else {
        (chain.reads[chain.reads.len() - 1] >> 32) as u32
    }
}

/// `drop_g_edges_by_utg` (`hic.cpp` 8381): deletes every bubble-graph arc
/// around the tangle's boundary chains.
fn drop_g_edges_by_utg(
    chain_unitigs: &[MaUtg],
    bsg: &mut Asg,
    res_btg: &[u32],
    b_utg_beg: u32,
    b_utg_sink: u32,
) {
    let mut drop_at = |root: u32| {
        let v = get_utg_end_from_btg(chain_unitigs, root);
        let (start, count) = bsg.arc_range(v);
        for i in 0..count {
            let arc = bsg.arcs[start + i];
            if arc.deleted {
                continue;
            }
            bsg.arcs[start + i].deleted = true;
            bsg.delete_arc(arc.v ^ 1, ((arc.ul >> 32) as u32) ^ 1, true);
        }
    };

    if b_utg_beg != u32::MAX {
        drop_at(b_utg_beg);
    }
    if b_utg_sink != u32::MAX {
        drop_at(b_utg_sink);
    }

    for &root in res_btg {
        drop_at(root);
        drop_at(root ^ 1);
    }
}

/// `update_bsg` (`hic.cpp` 8513): flushes the pending edges into the bubble
/// graph and keeps only the longest chain per source vertex.
fn update_bsg(bsg: &mut Asg, edges: &mut Vec<AsgArc>) {
    for arc in edges.iter() {
        bsg.push_arc(*arc);
    }
    // `bsg->is_srt = 0; free(bsg->idx); bsg->idx = 0; asg_cleanup(bsg);` —
    // the node growth of `update_bubble_graph` leaves the index sized for the
    // old vertex count; freeing forces the rebuild.
    bsg.is_sorted = false;
    bsg.index.clear();
    bsg.cleanup();

    crate::container::ksort::radix_sort_by_key_32(edges, &|arc: &AsgArc| (arc.ul >> 32) as u32);
    let mut l = 0usize;
    for k in 1..=edges.len() {
        if k == edges.len() || (edges[k].ul >> 32) != (edges[l].ul >> 32) {
            if k - l > 1 {
                let mut max = -1i64;
                let mut max_i = usize::MAX;
                for i in l..k {
                    let mut convex = 0u32;
                    let mut node_len = 0i64;
                    let mut base_len = 0i64;
                    let mut max_stop_node_len = 0i64;
                    let mut max_stop_base_len = 0i64;
                    crate::assembly::bubble_identify::get_unitig(
                        bsg,
                        None,
                        edges[i].v,
                        &mut convex,
                        &mut node_len,
                        &mut base_len,
                        &mut max_stop_node_len,
                        &mut max_stop_base_len,
                        1,
                        None,
                    );
                    if max < node_len {
                        max = node_len;
                        max_i = i;
                    }
                }

                for i in l..k {
                    if i == max_i {
                        continue;
                    }
                    bsg.delete_arc((edges[i].ul >> 32) as u32, edges[i].v, true);
                    bsg.delete_arc(edges[i].v ^ 1, ((edges[i].ul >> 32) as u32) ^ 1, true);
                }
            }
            l = k;
        }
    }

    bsg.cleanup();
}

/// Extracts `(bubble-graph vertex, read-unitig vertex)` from the start or
/// end entry of a bubble chain — the `get_bubbles(bub, (u->a[k]>>32)>>1, ...)`
/// + `utg = w^1` pattern of `find_bubble_tangle`.
fn chain_end_to_bubble(bub: &BubbleType, chain_unitigs: &[MaUtg], x: u32) -> (u32, u32) {
    let chain = &chain_unitigs[(x >> 1) as usize];
    let entry = if x & 1 == 1 {
        chain.reads[0]
    } else {
        chain.reads[chain.reads.len() - 1]
    };
    let bubble_vertex = (entry >> 32) as u32;
    let bubble_id = (entry >> 33) as u64;
    let (beg, sink, _, _) = get_bubbles(bub, bubble_id);
    let ori = (entry >> 32) & 1;
    let w = if ori == 0 { beg } else { sink };
    (bubble_vertex, w ^ 1)
}

/// The five outputs of [`find_bubble_tangle`] (`r_b_utg_beg` …
/// `r_utg_sink`).
#[derive(Debug, Clone, Copy, Default)]
pub struct TangleBounds {
    /// `r_b_utg_beg`.
    pub chain_beg: u32,
    /// `r_b_utg_sink`.
    pub chain_sink: u32,
    /// `r_b_tg_beg`.
    pub chain_bubble_beg: u32,
    /// `r_b_tg_sink`.
    pub chain_bubble_sink: u32,
    /// `r_utg_beg`.
    pub unitig_beg: u32,
    /// `r_utg_sink`.
    pub unitig_sink: u32,
}

/// `find_bubble_tangle` (`hic.cpp` 8185): walks the bubble chain from `v`
/// and decides whether it enters a tangle; on success reports the tangle's
/// chain-level and unitig-level boundaries.
#[allow(clippy::too_many_arguments)]
fn find_bubble_tangle(
    bub: &BubbleType,
    chain_unitigs: &[MaUtg],
    chain_graph: &Asg,
    unitig_graph: &Asg,
    unitigs: &[MaUtg],
    is_vis: &mut [u8],
    is_vis2: &mut [u8],
    v: u32,
    side_rate: f64,
    total_rate: f64,
    stack: &mut Vec<u32>,
    stack2: &mut Vec<u32>,
    res_btg: &mut Vec<u32>,
    res_utg: &mut Vec<u32>,
    chain_flag: &mut [u8],
) -> (bool, TangleBounds) {
    let mut bounds = TangleBounds::default();
    let mut tan_occ = 0u32;
    let mut no_root = false;
    stack.clear();
    is_vis.iter_mut().for_each(|flag| *flag = 0);
    stack.push(v);
    let v_occ = get_btg_occ(bub, chain_unitigs, unitigs, v >> 1);
    while let Some(cur) = stack.pop() {
        if is_vis[cur as usize] != 0 {
            continue;
        }

        let mut c_occ = 0u32;
        if no_root && cur == v {
            return (false, bounds);
        }
        if no_root && is_vis[cur as usize] == 0 && is_vis[(cur ^ 1) as usize] == 0 {
            c_occ = get_btg_occ(bub, chain_unitigs, unitigs, cur >> 1);
            if (tan_occ as f64) * side_rate < c_occ as f64
                && (tan_occ as f64) * total_rate < (c_occ + v_occ) as f64
            {
                if check_bubble_tangle(
                    Some(bub),
                    chain_graph,
                    Some(unitigs),
                    Some(unitigs),
                    v,
                    cur ^ 1,
                    side_rate,
                    total_rate,
                    v_occ,
                    c_occ,
                    is_vis2,
                    stack2,
                    None,
                    None,
                    None,
                ) {
                    let mut extra_check = 0u32;
                    check_bubble_tangle(
                        Some(bub),
                        chain_graph,
                        Some(unitigs),
                        Some(unitigs),
                        v,
                        cur ^ 1,
                        side_rate,
                        total_rate,
                        v_occ,
                        c_occ,
                        is_vis2,
                        stack2,
                        Some(res_btg),
                        None,
                        Some(&mut extra_check),
                    );

                    for &node in res_btg.iter() {
                        set_b_utg_weight_flag(bub, chain_unitigs, unitigs, node, Some(chain_flag), 0, None);
                    }
                    let btg_beg = v;
                    let btg_sink = cur ^ 1;

                    let (chain_bubble_beg, utg_beg) =
                        chain_end_to_bubble(bub, chain_unitigs, btg_beg);
                    let (chain_bubble_sink, utg_sink) =
                        chain_end_to_bubble(bub, chain_unitigs, btg_sink);
                    bounds.chain_bubble_beg = chain_bubble_beg;
                    bounds.chain_bubble_sink = chain_bubble_sink;

                    let mut is_t = check_bubble_tangle(
                        None,
                        unitig_graph,
                        None,
                        Some(unitigs),
                        utg_beg,
                        utg_sink,
                        side_rate,
                        total_rate,
                        u32::MAX,
                        u32::MAX,
                        is_vis2,
                        stack2,
                        Some(res_utg),
                        Some(chain_flag),
                        Some(&mut extra_check),
                    );

                    if is_t && extra_check == 0 {
                        let mut utg_occ = 0u32;
                        for &node in res_utg.iter() {
                            if bub.is_hom(node >> 1) {
                                continue;
                            }
                            utg_occ += unitigs[(node >> 1) as usize].reads.len() as u32;
                        }
                        if (utg_occ as f64) * total_rate >= (v_occ + c_occ) as f64 {
                            is_t = false;
                        }
                    }

                    for &node in res_btg.iter() {
                        set_b_utg_weight_flag(bub, chain_unitigs, unitigs, node, Some(chain_flag), 1, None);
                    }

                    if is_t {
                        bounds.chain_beg = btg_beg;
                        bounds.chain_sink = btg_sink;
                        bounds.unitig_beg = utg_beg;
                        bounds.unitig_sink = utg_sink;
                        return (true, bounds);
                    }
                }
            }
        }

        is_vis[cur as usize] = 1;
        if cur != v ^ 1 {
            for arc in chain_graph.arcs_of(cur) {
                if arc.deleted {
                    continue;
                }
                if arc.v == v {
                    return (false, bounds);
                }
                if is_vis[arc.v as usize] != 0 {
                    continue;
                }
                stack.push(arc.v);
            }
        }

        if no_root {
            tan_occ += c_occ;
        }
        if (tan_occ as f64) * side_rate >= v_occ as f64 {
            return (false, bounds);
        }
        no_root = true;
    }

    if (tan_occ as f64) * side_rate >= v_occ as f64 {
        return (false, bounds);
    }
    if (tan_occ as f64) * total_rate >= v_occ as f64 {
        return (false, bounds);
    }
    if !check_bubble_tangle(
        Some(bub),
        chain_graph,
        Some(unitigs),
        Some(unitigs),
        v,
        u32::MAX,
        side_rate,
        total_rate,
        v_occ,
        u32::MAX,
        is_vis2,
        stack2,
        Some(res_btg),
        None,
        None,
    ) {
        return (false, bounds);
    }

    for &node in res_btg.iter() {
        set_b_utg_weight_flag(bub, chain_unitigs, unitigs, node, Some(chain_flag), 0, None);
    }

    let btg_beg = v;
    let (chain_bubble_beg, utg_beg) = chain_end_to_bubble(bub, chain_unitigs, btg_beg);
    bounds.chain_bubble_beg = chain_bubble_beg;

    let mut extra_check = 0u32;
    let mut is_t = check_bubble_tangle(
        None,
        unitig_graph,
        None,
        Some(unitigs),
        utg_beg,
        u32::MAX,
        side_rate,
        total_rate,
        u32::MAX,
        u32::MAX,
        is_vis2,
        stack2,
        Some(res_utg),
        Some(chain_flag),
        Some(&mut extra_check),
    );

    if is_t && extra_check == 0 {
        let mut utg_occ = 0u32;
        for &node in res_utg.iter() {
            if bub.is_hom(node >> 1) {
                continue;
            }
            utg_occ += unitigs[(node >> 1) as usize].reads.len() as u32;
        }
        if (utg_occ as f64) * total_rate >= v_occ as f64 {
            is_t = false;
        }
    }

    for &node in res_btg.iter() {
        set_b_utg_weight_flag(bub, chain_unitigs, unitigs, node, Some(chain_flag), 1, None);
    }

    if is_t {
        bounds.chain_beg = btg_beg;
        bounds.unitig_beg = utg_beg;
    }
    (is_t, bounds)
}

/// `ma_2_bub_arc` (`hic.cpp` 6972): builds the bubble-graph arc between the
/// beg/sink vertices of two bubbles. Returns `false` when the vertices do
/// not bound the bubbles.
#[allow(clippy::too_many_arguments)]
pub fn ma_2_bub_arc(
    bub: &BubbleType,
    x_0: u32,
    mut x_0_b_id: Option<u64>,
    x_1: u32,
    mut x_1_b_id: Option<u64>,
    arc: &mut AsgArc,
    check_het: u32,
) -> bool {
    let mut id0 = x_0_b_id.unwrap_or(u64::MAX);
    let mut id1 = x_1_b_id.unwrap_or(u64::MAX);

    if x_0 != u32::MAX && x_1 != u32::MAX {
        if (x_0 >> 1) == (x_1 >> 1) {
            if x_0_b_id.is_none() && x_1_b_id.is_none() {
                let (a, b) = bub.get_bub_id(x_0 >> 1, check_het);
                id0 = a;
                id1 = b;
            }

            let (beg_0, sink_0, _, _) = get_bubbles(bub, id0);
            let (beg_1, sink_1, _, _) = get_bubbles(bub, id1);

            let mut ori_0 = u64::MAX;
            if x_0 == beg_0 ^ 1 {
                ori_0 = 1;
            } else if x_0 == sink_0 ^ 1 {
                ori_0 = 0;
            } else if x_0 == beg_1 ^ 1 {
                ori_0 = 3;
            } else if x_0 == sink_1 ^ 1 {
                ori_0 = 2;
            } else {
                eprintln!("error 0");
                return false;
            }

            let mut ori_1 = u64::MAX;
            if x_1 == beg_0 ^ 1 {
                ori_1 = 1;
            } else if x_1 == sink_0 ^ 1 {
                ori_1 = 0;
            } else if x_1 == beg_1 ^ 1 {
                ori_1 = 3;
            } else if x_1 == sink_1 ^ 1 {
                ori_1 = 2;
            } else {
                eprintln!("error 1");
                return false;
            }

            if ((ori_0 >> 1) ^ (ori_1 >> 1)) & 1 != 1 {
                eprintln!("error 10");
                return false;
            }

            if ori_0 & 2 != 0 {
                std::mem::swap(&mut id0, &mut id1);
            }

            ori_0 &= 1;
            ori_1 &= 1;
            ori_1 ^= 1;

            arc.ul = ((id0 << 1) | ori_0) << 32;
            arc.v = ((id1 << 1) | ori_1) as u32;
            arc.overlap = 0;
            arc.deleted = false;
            arc.exact_region = 1;
            arc.without_large_indel = true;
            arc.strong = true;
        } else {
            if x_0_b_id.is_none() {
                let (a, _) = bub.get_bub_id(x_0 >> 1, check_het);
                id0 = a;
            }
            if x_1_b_id.is_none() {
                let (a, _) = bub.get_bub_id(x_1 >> 1, check_het);
                id1 = a;
            }

            let (beg, sink, _, _) = get_bubbles(bub, id0);
            let mut ori_0 = u64::MAX;
            if x_0 == beg ^ 1 {
                ori_0 = 1;
            } else if x_0 == sink ^ 1 {
                ori_0 = 0;
            } else {
                eprintln!("error 0");
                return false;
            }

            let (beg, sink, _, _) = get_bubbles(bub, id1);
            let mut ori_1 = u64::MAX;
            if x_1 == beg ^ 1 {
                ori_1 = 1;
            } else if x_1 == sink ^ 1 {
                ori_1 = 0;
            } else {
                eprintln!("error 1");
                return false;
            }
            ori_0 &= 1;
            ori_1 &= 1;
            ori_1 ^= 1;

            arc.ul = ((id0 << 1) | ori_0) << 32;
            arc.v = ((id1 << 1) | ori_1) as u32;
            arc.overlap = 0;
            arc.deleted = false;
            arc.exact_region = 1;
            arc.without_large_indel = true;
            arc.strong = true;
        }
    } else {
        let mut x = u32::MAX;
        if x_0 != u32::MAX {
            x = x_0;
        }
        if x_1 != u32::MAX {
            x = x_1;
        }
        if x == u32::MAX {
            return false;
        }
        if x_0_b_id.is_none() && x_1_b_id.is_none() {
            let (a, b) = bub.get_bub_id(x >> 1, check_het);
            id0 = a;
            id1 = b;
        }

        if id0 != u64::MAX {
            let (beg, sink, _, _) = get_bubbles(bub, id0);
            if x == beg ^ 1 || x == sink ^ 1 {
                return true;
            }
            return false;
        }

        if id1 != u64::MAX {
            let (beg, sink, _, _) = get_bubbles(bub, id1);
            if x == beg ^ 1 || x == sink ^ 1 {
                return true;
            }
            return false;
        }
    }

    true
}

/// `set_bub_idx` (`hic.cpp` 7220): flips the hap status of the bubble
/// beg/sink vertices along one bubble chain segment.
fn set_bub_idx(
    bub: &mut BubbleType,
    chain: &[u64],
    unitig_graph: &Asg,
    beg_idx: i64,
    end_idx: i64,
    is_to_hom: u32,
    check_het: u32,
) {
    if beg_idx > end_idx {
        return;
    }
    for k in beg_idx..=end_idx {
        let entry = chain[k as usize];
        let read_id = (entry >> 33) as u32;
        let ori = (entry >> 32) & 1;
        let (beg, sink, _, _) = get_bubbles(bub, read_id as u64);
        let root = if ori == 1 { beg } else { sink };

        if is_to_hom != 0 && bub.is_hom(root >> 1) {
            continue;
        }
        if is_to_hom == 0 && bub.is_het(root >> 1) {
            continue;
        }

        let (id0, id1) = bub.get_bub_id(root >> 1, check_het);
        if id0 == u64::MAX || id1 == u64::MAX {
            continue;
        }
        let (_, _, _, len0) = get_bubbles(bub, id0);
        let (_, _, _, len1) = get_bubbles(bub, id1);

        if is_to_hom != 0 {
            if u64::from(unitig_graph.sequences[(root >> 1) as usize].length)
                > len0.min(len1) * 3
            {
                continue;
            }
            bub.index[(root >> 1) as usize] = u32::MAX;
        } else {
            bub.index[(root >> 1) as usize] = bub.full_bubble_count as u32 + 1;
        }
    }
}

/// `determine_bub_idx` (`hic.cpp` 7254).
#[allow(clippy::too_many_arguments)]
fn determine_bub_idx(
    bub: &mut BubbleType,
    chain: &[u64],
    unitig_graph: &Asg,
    p_len: u64,
    r_len: u64,
    r_het_len: u64,
    beg_idx: i64,
    end_idx: i64,
    check_het: u32,
) {
    if beg_idx > end_idx {
        return;
    }

    let r_hom_len = r_len - r_het_len;
    if p_len > 0
        && r_len > 0
        && r_het_len > 0
        && (r_len as f64) < (p_len as f64) * 0.5
        && (r_het_len as f64) < (r_len as f64) * 0.2
    {
        set_bub_idx(bub, chain, unitig_graph, beg_idx, end_idx, 1, check_het);
    } else if p_len > 0
        && r_len > 0
        && r_hom_len > 0
        && (r_len as f64) > (p_len as f64) * 0.9
        && (r_hom_len as f64) < (r_len as f64) * 0.1
    {
        set_bub_idx(bub, chain, unitig_graph, beg_idx, end_idx, 0, check_het);
    }
}

/// `detect_bub_graph` (`hic.cpp` 7271): walks every bubble chain of the
/// bubble graph and re-derives the het/hom classification of its beg/sink
/// unitigs from the chain path lengths.
fn detect_bub_graph(bub: &mut BubbleType, unitig_graph: &Asg) {
    let mut bsg = bub.bubble_graph.take().expect("bubble graph built");
    let chain_graph = crate::assembly::unitig::ma_ug_gen(&mut bsg);

    for chain in &chain_graph.unitigs {
        if chain.reads.is_empty() {
            continue;
        }
        let mut p_len = 0u64;
        let mut r_len = 0u64;
        let mut r_het_len = 0u64;
        let mut beg_idx = 0i64;
        let mut end_idx = -1i64;
        let mut root = 0u32;
        for k in 0..chain.reads.len() {
            let entry = chain.reads[k];
            let read_id = (entry >> 33) as u32;
            let ori = (entry >> 32) & 1;
            let (beg, sink, _, _) = get_bubbles(bub, read_id as u64);
            root = if ori == 1 { beg } else { sink };

            let t = if k + 1 < chain.reads.len() {
                let (start, count) = bsg.arc_range((entry >> 32) as u32);
                if count > 0 {
                    Some(bsg.arcs[start])
                } else {
                    None
                }
            } else {
                None
            };

            p_len += u64::from(bsg.sequences[read_id as usize].length);
            if end_idx < beg_idx {
                let length = unitig_graph.sequences[(root >> 1) as usize].length;
                p_len += u64::from(length);
                r_len += u64::from(length);
                if bub.is_het(root >> 1) {
                    r_het_len += u64::from(length);
                }
            }

            if let Some(t) = t {
                if t.exact_region == 0 {
                    if end_idx >= beg_idx {
                        let length = unitig_graph.sequences[(root >> 1) as usize].length;
                        p_len += u64::from(length);
                        r_len += u64::from(length);
                        if bub.is_het(root >> 1) {
                            r_het_len += u64::from(length);
                        }
                        determine_bub_idx(
                            bub,
                            &chain.reads,
                            unitig_graph,
                            p_len,
                            r_len,
                            r_het_len,
                            beg_idx,
                            end_idx,
                            bub.check_het,
                        );
                    }

                    p_len = 0;
                    r_len = 0;
                    r_het_len = 0;
                    beg_idx = k as i64 + 1;
                    end_idx = k as i64;
                } else {
                    p_len += u64::from(t.overlap);
                    r_len += u64::from(t.overlap);
                    if bub.is_het(root >> 1) {
                        r_het_len += u64::from(t.overlap);
                    }
                    end_idx = k as i64;
                }
            }
        }

        if end_idx >= beg_idx {
            let length = unitig_graph.sequences[(root >> 1) as usize].length;
            p_len += u64::from(length);
            r_len += u64::from(length);
            if bub.is_het(root >> 1) {
                r_het_len += u64::from(length);
            }
            determine_bub_idx(
                bub,
                &chain.reads,
                unitig_graph,
                p_len,
                r_len,
                r_het_len,
                beg_idx,
                end_idx,
                bub.check_het,
            );
        }
    }

    bub.bubble_graph = Some(bsg);
}

/// `get_bub_graph` (`hic.cpp` 7341): builds the bubble graph — one node per
/// bubble (length = path base), arcs between beg/sink vertices that are not
/// separated by another bubble.
fn get_bub_graph(unitig_graph: &Asg, bub: &mut BubbleType) -> Asg {
    let sg = unitig_graph;
    let mut pq = Pdq::new((sg.sequences.len() * 2) as u64);
    let n_vtx = (sg.sequences.len() * 2) as u32;
    let mut pre = vec![0u32; n_vtx as usize];
    let mut bub_g = Asg::new();
    for v in 0..bub.full_bubble_count as u32 {
        let (beg, sink, _, path_base) = get_bubbles(bub, v as u64);
        let het = bub.check_het != 0 && bub.is_het(beg >> 1) && bub.is_het(sink >> 1);
        bub_g.set_sequence(v as usize, path_base as u32, het);
        bub_g.sequences[v as usize].coverage_class = crate::assembly::hit::PRIMARY_LABEL;
    }

    for v in 0..n_vtx {
        if sg.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if bub.b_s_idx[(v >> 1) as usize] == u64::MAX {
            continue;
        }
        let bub_occ = bub.connect_bub_occ(v >> 1, bub.check_het);
        if bub_occ == 0 {
            continue;
        }
        if bub_occ == 2 {
            let mut t = AsgArc::default();
            if ma_2_bub_arc(bub, v, None, v ^ 1, None, &mut t, bub.check_het) {
                t.overlap = sg.sequences[(v >> 1) as usize].length;
                bub_g.push_arc(t);
            }
            continue;
        }
        let mut t = AsgArc::default();
        if !ma_2_bub_arc(bub, v, None, u32::MAX, None, &mut t, bub.check_het) {
            continue;
        }

        pq.shortest_path(v, sg, &mut pre);
        for k in 0..pq.distance.len() as u32 {
            if pq.distance[k as usize] == u64::MAX {
                continue;
            }
            if bub.b_s_idx[(k >> 1) as usize] == u64::MAX {
                continue;
            }
            if bub.connect_bub_occ(k >> 1, bub.check_het) == 0 {
                continue;
            }
            if (k >> 1) == (v >> 1) {
                continue;
            }
            let mut pre_id = pre[k as usize];
            let mut adjacent = false;
            while pre_id != v {
                if bub.connect_bub_occ(pre_id >> 1, bub.check_het) > 0 {
                    adjacent = true;
                    break;
                }
                pre_id = pre[pre_id as usize];
            }

            if !adjacent {
                let mut t = AsgArc::default();
                if ma_2_bub_arc(bub, v, None, k ^ 1, None, &mut t, bub.check_het) {
                    t.exact_region = 0;
                    t.overlap = (pq.distance[k as usize]
                        + u64::from(sg.sequences[(k >> 1) as usize].length))
                        as u32;
                    bub_g.push_arc(t);
                }
            }
        }
    }

    bub_g.cleanup();
    bub_g.read_sequence_count = bub_g.sequences.len() as u32;
    bub_g
}

/// `is_simple_broken_bubble` (`hic.cpp` 7476).
fn is_simple_broken_bubble(
    chain_graph: &Asg,
    x: u32,
    beg: u32,
    sink: u32,
    new_het: &mut u32,
) -> bool {
    *new_het = u32::MAX;
    let mut v = u32::MAX;
    if chain_graph.arc_count(x) == 1 && chain_graph.arc_count(x ^ 1) == 0 {
        v = x;
    }
    if chain_graph.arc_count(x ^ 1) == 1 && chain_graph.arc_count(x) == 0 {
        v = x ^ 1;
    }
    if v == u32::MAX {
        return false;
    }

    let mut u_s = u32::MAX;
    let mut u_e = u32::MAX;
    for arc in chain_graph.arcs_of(v) {
        if arc.deleted {
            continue;
        }
        if (arc.v >> 1) == (beg >> 1) {
            u_s = beg;
            u_e = sink;
        }
        if (arc.v >> 1) == (sink >> 1) {
            u_s = sink;
            u_e = beg;
        }
    }
    if u_s == u32::MAX || u_e == u32::MAX {
        return false;
    }

    let arcs = chain_graph.arcs_of(u_s);
    if chain_graph.arc_count(u_s) != 2 {
        return false;
    }
    for arc in arcs {
        if arc.deleted {
            continue;
        }
        if arc.v == v ^ 1 {
            continue;
        }
        if arc.v == u_e ^ 1 {
            *new_het = u_e;
            return true;
        }
    }
    false
}

/// `double_check_broken_bubble` (`hic.cpp` 7525): counts the vertices a
/// flood fill from `beg` reaches on both strands and rejects the merge when
/// it exceeds the collected `broken` set.
fn double_check_broken_bubble(
    graph: &Asg,
    broken: &mut Vec<u32>,
    beg: u32,
    sink: u32,
    visited: &mut [u8],
    stack: &mut Vec<u32>,
    bubble_graph: Option<&mut Asg>,
    previous_arc: Option<&mut AsgArc>,
) -> bool {
    crate::container::ksort::radix_sort_by_key_32(broken, &|entry: &u32| *entry);
    let mut n = 0usize;
    let mut pre = u32::MAX;
    for i in 0..broken.len() {
        if (broken[i] >> 1) == (pre >> 1) {
            continue;
        }
        pre = broken[i];
        broken[n] = pre;
        n += 1;
    }
    broken.truncate(n);

    visited.iter_mut().for_each(|flag| *flag = 0);
    stack.clear();
    stack.push(beg);
    let mut occ = 0u32;
    while let Some(mut cur) = stack.pop() {
        if visited[cur as usize] == 0 && visited[(cur ^ 1) as usize] == 0 {
            occ += 1;
        }
        if visited[cur as usize] == 0 && cur != beg ^ 1 && (sink == u32::MAX || cur != sink ^ 1) {
            visited[cur as usize] = 1;
            for arc in graph.arcs_of(cur) {
                if arc.deleted || visited[arc.v as usize] != 0 {
                    continue;
                }
                stack.push(arc.v);
            }
        }
        visited[cur as usize] = 1;

        cur ^= 1;
        if visited[cur as usize] == 0 && cur != beg ^ 1 && (sink == u32::MAX || cur != sink ^ 1) {
            visited[cur as usize] = 1;
            for arc in graph.arcs_of(cur) {
                if arc.deleted || visited[arc.v as usize] != 0 {
                    continue;
                }
                stack.push(arc.v);
            }
        }
        visited[cur as usize] = 1;
    }

    let mut n = broken.len();
    if beg != u32::MAX {
        n += 1;
    }
    if sink != u32::MAX {
        n += 1;
    }

    if occ > n as u32 {
        if let (Some(bsg), Some(p_t)) = (bubble_graph, previous_arc) {
            if beg != u32::MAX && sink != u32::MAX {
                p_t.deleted = true;
                bsg.delete_arc(p_t.v ^ 1, ((p_t.ul >> 32) as u32) ^ 1, true);
            }
        }
        return false;
    }

    true
}

/// `is_local_simple_circle` (`hic.cpp` 7641).
fn is_local_simple_circle(graph: &Asg, v: u32) -> bool {
    if graph.arc_count(v) != graph.arc_count(v ^ 1) {
        return false;
    }
    let mut v = v;
    if graph.arc_count(v) == 1 {
        let (start, _) = graph.arc_range(v);
        v = graph.arcs[start].v;
    }
    if graph.arc_count(v) != graph.arc_count(v ^ 1) {
        return false;
    }
    if graph.arc_count(v) != 2 {
        return false;
    }
    for arc in graph.arcs_of(v) {
        if arc.deleted {
            continue;
        }
        let u = arc.v;
        if graph.arc_count(u) != 1 || graph.arc_count(u ^ 1) != 1 {
            continue;
        }
        let (start, _) = graph.arc_range(u);
        if graph.arcs[start].v != v {
            continue;
        }
        return true;
    }
    false
}

/// `update_bubble_graph` (`hic.cpp` 7705): turns a collected `broken` set
/// into a new bubble node of the bubble graph and links it to the beg/sink
/// bubbles.
#[allow(clippy::too_many_arguments)]
fn update_bubble_graph(
    broken: &mut Vec<u32>,
    beg: u32,
    beg_bub_id: u32,
    sink: u32,
    sink_bub_id: u32,
    bub: &mut BubbleType,
    edges: &mut Vec<AsgArc>,
    bsg: &mut Asg,
    // `p_t`: the slot index of the arc between the two chain entries; the
    // reference holds a pointer into `bsg->arc` and marks it deleted in
    // place.
    previous_arc: Option<usize>,
    mut bsg_idx: Option<&mut [u8]>,
    chain_unitigs: &[MaUtg],
    chain_graph: &Asg,
    occurrence_threshold: Option<u64>,
    is_b_bub: bool,
) {
    crate::container::ksort::radix_sort_by_key_32(broken, &|entry: &u32| *entry);
    let mut n = 0usize;
    let mut occ = 0u64;
    let mut pre = u32::MAX;
    for i in 0..broken.len() {
        if (broken[i] >> 1) == (pre >> 1) {
            continue;
        }
        if bub.is_hom(broken[i] >> 1) {
            if is_local_simple_circle(chain_graph, broken[i]) {
                bub.index[(broken[i] >> 1) as usize] = bub.full_bubble_count as u32 + 1;
            } else {
                continue;
            }
        }

        pre = broken[i];
        broken[n] = pre;
        occ += chain_unitigs[(broken[n] >> 1) as usize].reads.len() as u64;
        n += 1;
    }
    broken.truncate(n);

    if broken.len() == 1 {
        let mut v = 0u32;
        if beg != u32::MAX
            && sink != u32::MAX
            && is_simple_broken_bubble(chain_graph, broken[0], beg, sink, &mut v)
        {
            if (v >> 1) != (broken[0] >> 1) {
                broken.push(v);
                occ += chain_unitigs[(v >> 1) as usize].reads.len() as u64;
                if let Some(thres) = occurrence_threshold {
                    if occ > thres {
                        return;
                    }
                }
                bub.index[(v >> 1) as usize] = bub.full_bubble_count as u32 + 1;
            }
        }
    }

    if let Some(thres) = occurrence_threshold {
        if occ > thres {
            return;
        }
    }

    let bub_id = bsg.sequences.len() as u32;
    bsg.set_sequence(bub_id as usize, 0, false);
    bsg.sequences[bub_id as usize].coverage_class = crate::assembly::hit::HAP_LABEL;
    if is_b_bub {
        bub.broken_bubble_count += 1;
    }

    bub.num.push(bub.list.len() as u32);
    bub.path_len.push(0);
    bub.list.push(beg);
    bub.list.push(sink);

    for &entry in broken.iter() {
        bub.list.push(entry);
        if let Some(idx) = bsg_idx.as_deref_mut() {
            idx[(entry >> 1) as usize] = 1;
        }
    }

    if beg != u32::MAX {
        let mut t_f = AsgArc::default();
        let mut t_r = AsgArc::default();
        if ma_2_bub_arc(bub, beg, Some(beg_bub_id as u64), beg ^ 1, Some(bub_id as u64), &mut t_f, bub.check_het)
            && ma_2_bub_arc(bub, beg ^ 1, Some(bub_id as u64), beg, Some(beg_bub_id as u64), &mut t_r, bub.check_het)
        {
            t_f.exact_region = 0;
            t_f.without_large_indel = false;
            t_f.deleted = false;
            edges.push(t_f);
            t_r.exact_region = 0;
            t_r.without_large_indel = false;
            t_r.deleted = false;
            edges.push(t_r);
        }
    }

    if sink != u32::MAX {
        let mut t_f = AsgArc::default();
        let mut t_r = AsgArc::default();
        if ma_2_bub_arc(bub, sink ^ 1, Some(bub_id as u64), sink, Some(sink_bub_id as u64), &mut t_f, bub.check_het)
            && ma_2_bub_arc(bub, sink, Some(sink_bub_id as u64), sink ^ 1, Some(bub_id as u64), &mut t_r, bub.check_het)
        {
            t_f.exact_region = 0;
            t_f.without_large_indel = false;
            t_f.deleted = false;
            edges.push(t_f);
            t_r.exact_region = 0;
            t_r.without_large_indel = false;
            t_r.deleted = false;
            edges.push(t_r);
        }
    }

    if beg != u32::MAX && sink != u32::MAX {
        if let Some(slot) = previous_arc {
            // `p_t->del = 1; asg_arc_del(bsg, p_t->v^1, (p_t->ul>>32)^1, 1);`
            bsg.arcs[slot].deleted = true;
            let v = bsg.arcs[slot].v;
            bsg.delete_arc(v ^ 1, ((bsg.arcs[slot].ul >> 32) as u32) ^ 1, true);
        }
    }
}

/// `get_related_bub_nodes` (`hic.cpp` 7807).
fn get_related_bub_nodes(
    broken: &mut Vec<u32>,
    bub: &BubbleType,
    pq: &mut Pdq,
    unitig_graph: &Asg,
    pre: &mut [u32],
    src: u32,
    dest: u32,
    bsg_idx: Option<&[u8]>,
) {
    let src = src ^ 1;
    pq.shortest_path(src, unitig_graph, pre);

    for j_i in 0..pq.distance.len() as u32 {
        if pq.distance[j_i as usize] == u64::MAX {
            continue;
        }
        if (j_i >> 1) == (src >> 1) {
            continue;
        }
        if dest != u32::MAX && (j_i >> 1) == (dest >> 1) {
            continue;
        }

        let mut pre_id = pre[j_i as usize];
        let mut adjacent = false;
        while pre_id != src {
            if (dest != u32::MAX && (pre_id >> 1) == (dest >> 1))
                || (pre_id >> 1) == (src >> 1)
            {
                adjacent = true;
                break;
            }
            pre_id = pre[pre_id as usize];
        }

        if !adjacent {
            if broken.is_empty() || (j_i >> 1) != (broken[broken.len() - 1] >> 1) {
                if let Some(idx) = bsg_idx {
                    if idx[(j_i >> 1) as usize] != 0 {
                        broken.clear();
                        return;
                    }
                }
                broken.push(j_i);
            }
        }
    }
}

/// `calculate_chain_weight` (`hic.cpp` 7900).
fn calculate_chain_weight(
    chain: &MaUtg,
    bub: &BubbleType,
    bubble_graph: &Asg,
    chain_unitigs: &[MaUtg],
    entry: &mut ChainWeightEntry,
) -> u64 {
    entry.bubble_occurrence = 0;
    entry.graph_occurrence = 0;
    let mut occ = 0u64;
    let mut occ_n = 0u64;
    for read in &chain.reads {
        let bubble_id = (read >> 33) as usize;
        if bubble_graph.sequences[bubble_id].coverage_class != crate::assembly::hit::HAP_LABEL {
            let (_, _, body, _) = get_bubbles(bub, bubble_id as u64);
            for &node in body {
                occ += chain_unitigs[(node >> 1) as usize].reads.len() as u64;
            }
            occ_n += 1;
        }
    }

    let mut thres = u64::MAX;
    if occ_n > 0 {
        thres = (occ * 6) / occ_n;
    }

    for read in &chain.reads {
        let bubble_id = (read >> 33) as usize;
        let mut occ = 0u64;
        let (_, _, body, _) = get_bubbles(bub, bubble_id as u64);
        for &node in body {
            occ += chain_unitigs[(node >> 1) as usize].reads.len() as u64;
        }

        if bubble_graph.sequences[bubble_id].coverage_class != crate::assembly::hit::HAP_LABEL
            || occ < thres
        {
            entry.graph_occurrence += occ as i64;
        } else {
            entry.bubble_occurrence += occ as i64;
        }
    }

    thres
}

/// `get_btg_occ` (`hic.cpp` 8074): the summed read count of every unitig in
/// the bubbles of chain `v`. `unitigs` is the reference's `bub->ug->u`.
fn get_btg_occ(bub: &BubbleType, chain_unitigs: &[MaUtg], unitigs: &[MaUtg], v: u32) -> u32 {
    let mut tan_occ = 0u32;
    let chain = &chain_unitigs[v as usize];
    for read in &chain.reads {
        let (_, _, body, _) = get_bubbles(bub, (read >> 33) as u64);
        for &node in body {
            tan_occ += unitigs[(node >> 1) as usize].reads.len() as u32;
        }
    }
    tan_occ
}

/// `set_b_utg_weight_flag` (`hic.cpp` 8960), the `b == NULL` shape used by
/// every trio-path call site (the `buf_t` variant is Hi-C only).
fn set_b_utg_weight_flag(
    bub: &BubbleType,
    chain_unitigs: &[MaUtg],
    unitigs: &[MaUtg],
    v: u32,
    mut vis_flag: Option<&mut [u8]>,
    flag: u32,
    mut occurrence: Option<&mut u32>,
) {
    if let Some(occ) = occurrence.as_deref_mut() {
        *occ = 0;
    }
    let chain = &chain_unitigs[(v >> 1) as usize];
    if chain.reads.is_empty() {
        return;
    }
    for read in &chain.reads {
        let (beg, sink, body, _) = get_bubbles(bub, (read >> 33) as u64);
        for &node in body {
            if let Some(vis) = vis_flag.as_deref_mut() {
                vis[(node >> 1) as usize] = flag as u8;
            }
            if let Some(occ) = occurrence.as_deref_mut() {
                *occ += unitigs[(node >> 1) as usize].reads.len() as u32;
            }
        }
        if beg != u32::MAX {
            if let Some(vis) = vis_flag.as_deref_mut() {
                vis[(beg >> 1) as usize] = flag as u8;
            }
        }
        if sink != u32::MAX {
            if let Some(vis) = vis_flag.as_deref_mut() {
                vis[(sink >> 1) as usize] = flag as u8;
            }
        }
    }
}

/// `resolve_bubble_chain_tangle` (`hic.cpp` 8562): repeatedly finds tangles
/// between bubble chains, converts them into new bubbles and rebuilds the
/// chain graph until nothing changes; finishes by weighing the chains.
fn resolve_bubble_chain_tangle(unitig_graph: &Asg, unitigs: &[MaUtg], bsg: &mut Asg, bub: &mut BubbleType) {
    let mut is_vis: Vec<u8>;
    let mut n_vx = {
        let chain_graph = bub.bubble_chain_graph.as_ref().expect("chain graph built");
        (unitig_graph.sequences.len() * 2)
            .max(chain_graph.sequences.len() * 2)
            .max(bsg.sequences.len() * 2)
    };
    is_vis = vec![0; n_vx];
    let mut is_vis2 = vec![0; n_vx];
    let mut is_used = vec![0; n_vx];
    let mut chain_flag = vec![0; n_vx];
    let mut stack: Vec<u32> = Vec::new();
    let mut stack2: Vec<u32> = Vec::new();
    let mut res_btg: Vec<u32> = Vec::new();
    let mut res_utg: Vec<u32> = Vec::new();
    let side_rate = 2.5f64;
    let total_rate = 8f64;

    loop {
        res_btg.clear();
        res_utg.clear();
        let mut occ_idx: Vec<u64> = Vec::new();
        let mut edges: Vec<AsgArc> = Vec::new();
        {
            let chain_graph = bub.bubble_chain_graph.as_ref().expect("chain graph built");
            let needed = (unitig_graph.sequences.len() * 2)
                .max(chain_graph.sequences.len() * 2)
                .max(bsg.sequences.len() * 2);
            if n_vx < needed {
                n_vx = needed;
                is_vis.resize(n_vx, 0);
                is_vis2.resize(n_vx, 0);
                is_used.resize(n_vx, 0);
                chain_flag.resize(n_vx, 0);
            }
        }
        is_vis.iter_mut().for_each(|flag| *flag = 0);
        is_vis2.iter_mut().for_each(|flag| *flag = 0);
        is_used.iter_mut().for_each(|flag| *flag = 0);
        chain_flag.iter_mut().for_each(|flag| *flag = 0);
        if !bub.num.is_empty() {
            bub.num.pop();
        }
        let new_bub_base = bsg.sequences.len() as u32;

        let chain_unitigs = bub.bubble_unitigs.take().expect("chain unitigs built");
        let chain_graph = bub.bubble_chain_graph.take().expect("chain graph built");
        for k in 0..chain_unitigs.len() {
            let occurrence = get_btg_occ(bub, &chain_unitigs, unitigs, k as u32);
            occ_idx.push(((occurrence as u64) << 32) + k as u64);
            set_b_utg_weight_flag(bub, &chain_unitigs, unitigs, (k as u32) << 1, Some(&mut chain_flag), 1, None);
        }
        crate::container::ksort::radix_sort_u64(&mut occ_idx);
        occ_idx.reverse();

        for k in 0..chain_unitigs.len() {
            let v = ((occ_idx[k] as u32) & 0xffffffff) << 1;
            if is_used[v as usize] == 0 && chain_graph.arc_count(v) > 0 {
                let (found, bounds) = find_bubble_tangle(
                    bub,
                    &chain_unitigs,
                    &chain_graph,
                    unitig_graph,
                    unitigs,
                    &mut is_vis,
                    &mut is_vis2,
                    v,
                    side_rate,
                    total_rate,
                    &mut stack,
                    &mut stack2,
                    &mut res_btg,
                    &mut res_utg,
                    &mut chain_flag,
                );
                if found {
                    if bounds.unitig_beg != u32::MAX && !bub.is_hom(bounds.unitig_beg >> 1) {
                        res_utg.push(bounds.unitig_beg);
                    }
                    if bounds.unitig_sink != u32::MAX && !bub.is_hom(bounds.unitig_sink >> 1) {
                        res_utg.push(bounds.unitig_sink);
                    }

                    for &node in res_btg.iter() {
                        is_used[node as usize] = 1;
                        is_used[(node ^ 1) as usize] = 1;
                    }
                    if bounds.chain_beg != u32::MAX {
                        is_used[bounds.chain_beg as usize] = 1;
                    }
                    if bounds.chain_sink != u32::MAX {
                        is_used[bounds.chain_sink as usize] = 1;
                    }
                    let mut b_tg_beg = bounds.chain_bubble_beg;
                    let mut b_tg_sink = bounds.chain_bubble_sink;
                    if b_tg_beg != u32::MAX {
                        b_tg_beg >>= 1;
                    }
                    if b_tg_sink != u32::MAX {
                        b_tg_sink >>= 1;
                    }
                    let mut broken = res_utg.clone();
                    update_bubble_graph(
                        &mut broken,
                        bounds.unitig_beg,
                        b_tg_beg,
                        bounds.unitig_sink,
                        b_tg_sink,
                        bub,
                        &mut edges,
                        bsg,
                        None,
                        None,
                        unitigs,
                        unitig_graph,
                        None,
                        false,
                    );
                    drop_g_edges_by_utg(&chain_unitigs, bsg, &res_btg, bounds.chain_beg, bounds.chain_sink);
                }
            }

            let v = v ^ 1;
            if is_used[v as usize] == 0 && chain_graph.arc_count(v) > 0 {
                let (found, bounds) = find_bubble_tangle(
                    bub,
                    &chain_unitigs,
                    &chain_graph,
                    unitig_graph,
                    unitigs,
                    &mut is_vis,
                    &mut is_vis2,
                    v,
                    side_rate,
                    total_rate,
                    &mut stack,
                    &mut stack2,
                    &mut res_btg,
                    &mut res_utg,
                    &mut chain_flag,
                );
                if found {
                    if bounds.unitig_beg != u32::MAX && !bub.is_hom(bounds.unitig_beg >> 1) {
                        res_utg.push(bounds.unitig_beg);
                    }
                    if bounds.unitig_sink != u32::MAX && !bub.is_hom(bounds.unitig_sink >> 1) {
                        res_utg.push(bounds.unitig_sink);
                    }

                    for &node in res_btg.iter() {
                        is_used[node as usize] = 1;
                        is_used[(node ^ 1) as usize] = 1;
                    }
                    if bounds.chain_beg != u32::MAX {
                        is_used[bounds.chain_beg as usize] = 1;
                    }
                    if bounds.chain_sink != u32::MAX {
                        is_used[bounds.chain_sink as usize] = 1;
                    }
                    let mut b_tg_beg = bounds.chain_bubble_beg;
                    let mut b_tg_sink = bounds.chain_bubble_sink;
                    if b_tg_beg != u32::MAX {
                        b_tg_beg >>= 1;
                    }
                    if b_tg_sink != u32::MAX {
                        b_tg_sink >>= 1;
                    }
                    let mut broken = res_utg.clone();
                    update_bubble_graph(
                        &mut broken,
                        bounds.unitig_beg,
                        b_tg_beg,
                        bounds.unitig_sink,
                        b_tg_sink,
                        bub,
                        &mut edges,
                        bsg,
                        None,
                        None,
                        unitigs,
                        unitig_graph,
                        None,
                        false,
                    );
                    drop_g_edges_by_utg(&chain_unitigs, bsg, &res_btg, bounds.chain_beg, bounds.chain_sink);
                }
            }
        }
        bub.num.push(bub.list.len() as u32);
        let new_bub = bsg.sequences.len() as u32 - new_bub_base;
        bub.tangle_bubble_count += new_bub as u64;
        if new_bub != 0 {
            update_bub_b_s_idx(bub);
        }
        update_bsg(bsg, &mut edges);

        let chain_graph = crate::assembly::unitig::ma_ug_gen(bsg);
        bub.bubble_unitigs = Some(chain_graph.unitigs);
        bub.bubble_chain_graph = Some(chain_graph.graph);
        if new_bub == 0 {
            break;
        }
    }

    bub.chain_weight.clear();
    let chain_unitigs = bub.bubble_unitigs.as_ref().expect("chain unitigs built");
    for i in 0..chain_unitigs.len() {
        let mut entry = ChainWeightEntry {
            id: i as u64,
            deleted: 0,
            ..Default::default()
        };
        calculate_chain_weight(&chain_unitigs[i], bub, bsg, unitigs, &mut entry);
        bub.chain_weight.push(entry);
    }
    bub.chain_weight.sort_by(|a, b| {
        if a.deleted != b.deleted {
            return if a.deleted > b.deleted { std::cmp::Ordering::Greater } else { std::cmp::Ordering::Less };
        }
        let a_occ = a.graph_occurrence - a.bubble_occurrence;
        let b_occ = b.graph_occurrence - b.bubble_occurrence;
        b_occ.cmp(&a_occ)
    });
}

/// `update_bubble_chain` (`hic.cpp` 8736): extends broken bubbles at chain
/// joints (`is_middle`) and at chain ends (`is_end`).
fn update_bubble_chain(
    unitig_graph: &Asg,
    unitigs: &[MaUtg],
    bsg: &mut Asg,
    bub: &mut BubbleType,
    is_middle: bool,
    is_end: bool,
) {
    bub.bubble_unitigs = None;
    bub.bubble_chain_graph = None;
    bub.chain_weight.clear();
    let mut broken: Vec<u32> = Vec::new();
    let mut edges: Vec<AsgArc> = Vec::new();
    let mut pq = Pdq::new((unitig_graph.sequences.len() * 2) as u64);
    let mut bub_ug = crate::assembly::unitig::ma_ug_gen(bsg);
    let n_vtx = (unitig_graph.sequences.len() * 2) as usize;
    let mut pre = vec![0u32; n_vtx];
    let mut vis_flag = vec![0u8; n_vtx];
    let mut stack: Vec<u32> = Vec::new();
    let mut bsg_idx = vec![0u8; n_vtx >> 1];

    for chain in &bub_ug.unitigs {
        if chain.reads.is_empty() {
            continue;
        }
        for entry in &chain.reads {
            let (root_0, root_1, body, _) = get_bubbles(bub, (entry >> 33) as u64);
            for &node in body {
                bsg_idx[(node >> 1) as usize] = 1;
            }
            bsg_idx[(root_0 >> 1) as usize] = 1;
            bsg_idx[(root_1 >> 1) as usize] = 1;
        }
    }

    if !bub.num.is_empty() {
        bub.num.pop();
    }
    let new_bub_base = bsg.sequences.len() as u32;
    for i in 0..bub_ug.unitigs.len() {
        let chain = &bub_ug.unitigs[i];
        if chain.reads.is_empty() {
            continue;
        }
        if is_middle {
            for k_i in 0..chain.reads.len() {
                if k_i + 1 >= chain.reads.len() {
                    continue;
                }
                // `arc_first(g, v) = g->arc[g->idx[v]>>32]`: the slot at the
                // offset, regardless of the arc count.
                let (start, _count) = bsg.arc_range((chain.reads[k_i] >> 32) as u32);
                if bsg.arcs[start].exact_region == 1 {
                    continue;
                }

                let entry_0 = chain.reads[k_i];
                let r_id_0 = (entry_0 >> 33) as u32;
                let ori_0 = (entry_0 >> 32) & 1;
                let (beg_0, sink_0, _, _) = get_bubbles(bub, r_id_0 as u64);
                let root_0 = if ori_0 == 1 { beg_0 } else { sink_0 };

                let entry_1 = chain.reads[k_i + 1];
                let r_id_1 = (entry_1 >> 33) as u32;
                let ori_1 = ((entry_1 >> 32) & 1) ^ 1;
                let (beg_1, sink_1, _, _) = get_bubbles(bub, r_id_1 as u64);
                let root_1 = if ori_1 == 1 { beg_1 } else { sink_1 };

                broken.clear();
                get_related_bub_nodes(&mut broken, bub, &mut pq, unitig_graph, &mut pre, root_0, root_1, None);
                get_related_bub_nodes(&mut broken, bub, &mut pq, unitig_graph, &mut pre, root_1, root_0, None);
                if !double_check_broken_bubble(
                    unitig_graph,
                    &mut broken,
                    root_0 ^ 1,
                    root_1 ^ 1,
                    &mut vis_flag,
                    &mut stack,
                    None,
                    None,
                ) {
                    continue;
                }
                if !bub.is_hom(root_0 >> 1) {
                    broken.push(root_0);
                }
                if !bub.is_hom(root_1 >> 1) {
                    broken.push(root_1);
                }
                if !broken.is_empty() {
                    update_bubble_graph(
                        &mut broken,
                        root_0 ^ 1,
                        r_id_0,
                        root_1 ^ 1,
                        r_id_1,
                        bub,
                        &mut edges,
                        bsg,
                        Some(start),
                        Some(&mut bsg_idx),
                        unitigs,
                        unitig_graph,
                        None,
                        true,
                    );
                }
            }
        }

        if is_end {
            if !chain.reads.is_empty() && bub_ug.graph.arc_count(((i as u32) << 1) + 1) == 0 {
                let entry_0 = chain.reads[0];
                let r_id_0 = (entry_0 >> 33) as u32;
                let ori_0 = ((entry_0 >> 32) & 1) ^ 1;
                let (beg_0, sink_0, _, _) = get_bubbles(bub, r_id_0 as u64);
                let root_0 = if ori_0 == 1 { beg_0 } else { sink_0 };
                broken.clear();
                get_related_bub_nodes(&mut broken, bub, &mut pq, unitig_graph, &mut pre, root_0, u32::MAX, None);

                if double_check_broken_bubble(
                    unitig_graph,
                    &mut broken,
                    root_0 ^ 1,
                    u32::MAX,
                    &mut vis_flag,
                    &mut stack,
                    None,
                    None,
                ) {
                    if !bub.is_hom(root_0 >> 1) {
                        broken.push(root_0);
                    }
                    if !broken.is_empty() {
                        update_bubble_graph(
                            &mut broken,
                            root_0 ^ 1,
                            r_id_0,
                            u32::MAX,
                            u32::MAX,
                            bub,
                            &mut edges,
                            bsg,
                            None,
                            Some(&mut bsg_idx),
                            unitigs,
                            unitig_graph,
                            None,
                            false,
                        );
                    }
                }
            }

            if !chain.reads.is_empty() && bub_ug.graph.arc_count((i as u32) << 1) == 0 {
                let entry_1 = chain.reads[chain.reads.len() - 1];
                let r_id_1 = (entry_1 >> 33) as u32;
                let ori_1 = (entry_1 >> 32) & 1;
                let (beg_1, sink_1, _, _) = get_bubbles(bub, r_id_1 as u64);
                let root_1 = if ori_1 == 1 { beg_1 } else { sink_1 };
                broken.clear();
                get_related_bub_nodes(&mut broken, bub, &mut pq, unitig_graph, &mut pre, root_1, u32::MAX, Some(&bsg_idx));

                if double_check_broken_bubble(
                    unitig_graph,
                    &mut broken,
                    root_1 ^ 1,
                    u32::MAX,
                    &mut vis_flag,
                    &mut stack,
                    None,
                    None,
                ) {
                    if !bub.is_hom(root_1 >> 1) {
                        broken.push(root_1);
                    }
                    if !broken.is_empty() {
                        update_bubble_graph(
                            &mut broken,
                            u32::MAX,
                            u32::MAX,
                            root_1 ^ 1,
                            r_id_1,
                            bub,
                            &mut edges,
                            bsg,
                            None,
                            Some(&mut bsg_idx),
                            unitigs,
                            unitig_graph,
                            None,
                            false,
                        );
                    }
                }
            }
        }
    }
    bub.num.push(bub.list.len() as u32);
    let new_bub = bsg.sequences.len() as u32 - new_bub_base;
    if is_end {
        bub.broken_end_bubble_count += new_bub as u64;
    }
    if new_bub != 0 {
        update_bub_b_s_idx(bub);
    }

    for arc in edges.iter() {
        bsg.push_arc(*arc);
    }
    // `bsg->is_srt = 0; free(bsg->idx); bsg->idx = 0; asg_cleanup(bsg);` —
    // the node growth of `update_bubble_graph` leaves the index sized for the
    // old vertex count; freeing forces the rebuild.
    bsg.is_sorted = false;
    bsg.index.clear();
    bsg.cleanup();

    let chain_graph = crate::assembly::unitig::ma_ug_gen(bsg);
    bub.bubble_unitigs = Some(chain_graph.unitigs);
    bub.bubble_chain_graph = Some(chain_graph.graph);
    bub.chain_weight.clear();
    let chain_unitigs = bub.bubble_unitigs.take().expect("chain unitigs built");
    for i in 0..chain_unitigs.len() {
        let mut entry = ChainWeightEntry {
            id: i as u64,
            deleted: 0,
            ..Default::default()
        };
        calculate_chain_weight(&chain_unitigs[i], bub, bsg, unitigs, &mut entry);
        bub.chain_weight.push(entry);
    }
    bub.bubble_unitigs = Some(chain_unitigs);
    bub.chain_weight.sort_by(|a, b| {
        if a.deleted != b.deleted {
            return if a.deleted > b.deleted { std::cmp::Ordering::Greater } else { std::cmp::Ordering::Less };
        }
        let a_occ = a.graph_occurrence - a.bubble_occurrence;
        let b_occ = b.graph_occurrence - b.bubble_occurrence;
        b_occ.cmp(&a_occ)
    });
}

/// `build_bub_graph` (`hic.cpp` 9910): builds the bubble graph twice
/// (`check_het` 0 then 1) and resolves its chains and tangles.
pub fn build_bub_graph(unitig_graph: &mut Asg, unitigs: &[MaUtg], bub: &mut BubbleType) {
    bub.check_het = 0;
    let bsg = get_bub_graph(unitig_graph, bub);
    bub.bubble_graph = Some(bsg);
    detect_bub_graph(bub, unitig_graph);
    let _ = bub.bubble_graph.take();
    bub.check_het = 1;
    let mut bsg = get_bub_graph(unitig_graph, bub);
    // The reference assigns `bub->b_g = bub_g` inside `get_bub_graph`, so the
    // field is live during `update_bubble_chain`; the port hands the graph to
    // `calculate_chain_weight` through a parameter instead, because Rust
    // cannot alias the `&mut bsg` and the field.
    update_bubble_chain(unitig_graph, unitigs, &mut bsg, bub, true, false);
    update_bubble_chain(unitig_graph, unitigs, &mut bsg, bub, false, true);
    resolve_bubble_chain_tangle(unitig_graph, unitigs, &mut bsg, bub);
    bub.bubble_graph = Some(bsg);
}

/// `get_num_edges2existing_nodes_advance_by_broken_bub` (`Overlaps.cpp`
/// 33258): counts the distinct existing unitigs the query's deleted-overlap
/// hits reach; `expect_vis` filters the targets, `utg_vis` deduplicates
/// within one call (cleared at the first valid hit, like the reference).
#[allow(clippy::too_many_arguments)]
fn get_num_edges2existing_nodes_advance_by_broken_bub(
    unitigs: &[MaUtg],
    unitig_graph: &Asg,
    read_graph: &Asg,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    read_to_unitig: &ReadToUnitig,
    trio_labels: &[u8],
    parameters: crate::assembly::rescue::RescueParameters,
    expect_vis: &[u8],
    query: u32,
    mut utg_vis: Option<&mut Vec<u8>>,
    ignore_trio_flag: u32,
) -> (u32, u32) {
    let mut overlap_length = 0u32;
    let query_read = (query >> 1) as usize;
    let hits = &sources[query_read].hits;
    let trio_flag = trio_labels[query_read];
    let mut non_trio_flag = u32::MAX;
    if ignore_trio_flag == 0 {
        if trio_flag == crate::command_line::options::TRIO_FATHER {
            non_trio_flag = u32::from(crate::command_line::options::TRIO_MOTHER);
        }
        if trio_flag == crate::command_line::options::TRIO_MOTHER {
            non_trio_flag = u32::from(crate::command_line::options::TRIO_FATHER);
        }
    }

    let mut occ = 0u32;
    let mut is_first = true;
    for h in hits {
        let query_cut = &coverage_cut[h.query_identifier as usize];
        let target_cut = &coverage_cut[h.target_identifier as usize];
        if target_cut.deleted || read_graph.sequences[h.target_identifier as usize].deleted {
            continue;
        }
        if u32::from(trio_labels[h.target_identifier as usize]) == non_trio_flag {
            continue;
        }

        let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
        let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
        let (result, arc) = crate::assembly::hit::ma_hit_to_arc(
            h,
            query_length,
            target_length,
            parameters.maximum_hang,
            parameters.maximum_hang_rate,
            parameters.minimum_overlap,
        );
        if result < 0 {
            continue;
        }
        if (arc.ul >> 32) as u32 != query {
            continue;
        }
        let (u_id, is_unitig) = read_to_unitig.get(arc.v >> 1);
        if u_id == u32::MAX || is_unitig == 0 || unitig_graph.sequences[u_id as usize].deleted {
            continue;
        }
        if expect_vis[(arc.v >> 1) as usize] == 0 {
            continue;
        }

        if is_first {
            if let Some(vis) = utg_vis.as_deref_mut() {
                vis.iter_mut().for_each(|flag| *flag = 0);
            }
        }
        let counted = match utg_vis.as_deref_mut() {
            None => true,
            Some(vis) => vis[u_id as usize] == 0,
        };
        if counted {
            occ += 1;
        }
        if let Some(vis) = utg_vis.as_deref_mut() {
            vis[u_id as usize] = 1;
        }

        overlap_length = overlap_length.wrapping_add(arc.overlap);
        is_first = false;
    }

    (occ, overlap_length)
}

/// `get_edge2existing_node_advance_by_broken_bub` (`Overlaps.cpp` 33409).
#[allow(clippy::too_many_arguments)]
fn get_edge2existing_node_advance_by_broken_bub(
    unitigs: &[MaUtg],
    unitig_graph: &Asg,
    read_graph: &Asg,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    read_to_unitig: &ReadToUnitig,
    trio_labels: &[u8],
    parameters: crate::assembly::rescue::RescueParameters,
    expect_vis: &[u8],
    query: u32,
    index: &mut usize,
    arc_out: &mut AsgArc,
    ignore_trio_flag: u32,
) -> bool {
    let _ = unitigs;
    let query_read = (query >> 1) as usize;
    let hits = &sources[query_read].hits;
    let trio_flag = trio_labels[query_read];
    let mut non_trio_flag = u32::MAX;
    if ignore_trio_flag == 0 {
        if trio_flag == crate::command_line::options::TRIO_FATHER {
            non_trio_flag = u32::from(crate::command_line::options::TRIO_MOTHER);
        }
        if trio_flag == crate::command_line::options::TRIO_MOTHER {
            non_trio_flag = u32::from(crate::command_line::options::TRIO_FATHER);
        }
    }

    while *index < hits.len() {
        let h = &hits[*index];
        *index += 1;
        let query_cut = &coverage_cut[h.query_identifier as usize];
        let target_cut = &coverage_cut[h.target_identifier as usize];
        if target_cut.deleted || read_graph.sequences[h.target_identifier as usize].deleted {
            continue;
        }
        if u32::from(trio_labels[h.target_identifier as usize]) == non_trio_flag {
            continue;
        }

        let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
        let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
        let (result, arc) = crate::assembly::hit::ma_hit_to_arc(
            h,
            query_length,
            target_length,
            parameters.maximum_hang,
            parameters.maximum_hang_rate,
            parameters.minimum_overlap,
        );
        if result < 0 {
            continue;
        }
        *arc_out = arc;

        if (arc_out.ul >> 32) as u32 != query {
            continue;
        }
        let (u_id, is_unitig) = read_to_unitig.get(arc_out.v >> 1);
        if u_id == u32::MAX || is_unitig == 0 || unitig_graph.sequences[u_id as usize].deleted {
            continue;
        }
        if expect_vis[(arc_out.v >> 1) as usize] == 0 {
            continue;
        }

        return true;
    }

    false
}

/// `get_contained_reads_chain_by_broken_bub` (`Overlaps.cpp` 34276): walks
/// the chain of deleted (contained) overlaps starting at `first_hit` until
/// it either reaches an existing unitig (`true`) or dead-ends.
#[allow(clippy::too_many_arguments)]
fn get_contained_reads_chain_by_broken_bub(
    unitigs: &[MaUtg],
    unitig_graph: &Asg,
    read_graph: &Asg,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    read_to_unitig: &ReadToUnitig,
    trio_labels: &[u8],
    parameters: crate::assembly::rescue::RescueParameters,
    mut end_rid: u32,
    first_hit: &PafHit,
    chain_edges: &mut Vec<AsgArc>,
    return_ava_cur: &mut u32,
    return_ava_ol: &mut u32,
    return_chain_len: &mut u32,
    expect_vis: &[u8],
    circle_vis: &mut [u8],
    mut utg_vis: Option<&mut Vec<u8>>,
    ignore_trio_flag: u32,
) -> bool {
    *return_chain_len = u32::MAX;
    *return_ava_cur = u32::MAX;
    *return_ava_ol = u32::MAX;
    chain_edges.clear();
    let mut chain_len = 0u32;

    let mut h: Option<&PafHit> = Some(first_hit);
    while let Some(hit) = h {
        let query_cut = &coverage_cut[hit.query_identifier as usize];
        let target_cut = &coverage_cut[hit.target_identifier as usize];
        let trio_flag = trio_labels[hit.query_identifier as usize];

        let mut non_trio_flag = u32::MAX;
        if ignore_trio_flag == 0 {
            if trio_flag == crate::command_line::options::TRIO_FATHER {
                non_trio_flag = u32::from(crate::command_line::options::TRIO_MOTHER);
            }
            if trio_flag == crate::command_line::options::TRIO_MOTHER {
                non_trio_flag = u32::from(crate::command_line::options::TRIO_FATHER);
            }
            if u32::from(trio_labels[hit.target_identifier as usize]) == non_trio_flag {
                break;
            }
        }

        // Only deleted edges of deleted targets describe contained reads.
        if !target_cut.deleted {
            break;
        }
        if !hit.deleted {
            break;
        }

        let (contain_rid, is_unitig) = read_to_unitig.get(hit.target_identifier);
        if contain_rid == u32::MAX || is_unitig == 1 {
            break;
        }
        if read_graph.sequences[contain_rid as usize].deleted {
            break;
        }

        let (contain_uid, is_unitig) = read_to_unitig.get(contain_rid);
        if contain_uid == u32::MAX || is_unitig != 1 {
            break;
        }
        if unitig_graph.sequences[contain_uid as usize].deleted {
            break;
        }

        let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
        let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
        let (result, t) = crate::assembly::hit::ma_hit_to_arc(
            hit,
            query_length,
            target_length,
            parameters.maximum_hang,
            parameters.maximum_hang_rate,
            parameters.minimum_overlap,
        );
        if result < 0 {
            break;
        }
        if (t.ul >> 32) as u32 != end_rid {
            break;
        }

        chain_edges.push(t);
        chain_len += 1;

        if circle_vis[(t.ul >> 33) as usize] != 0 || circle_vis[(t.v >> 1) as usize] != 0 {
            break;
        }
        circle_vis[(t.ul >> 33) as usize] = 1;
        circle_vis[(t.v >> 1) as usize] = 1;

        let (ava_cur, test_o_len) = get_num_edges2existing_nodes_advance_by_broken_bub(
            unitigs,
            unitig_graph,
            read_graph,
            sources,
            coverage_cut,
            read_to_unitig,
            trio_labels,
            parameters,
            expect_vis,
            t.v,
            utg_vis.as_deref_mut(),
            ignore_trio_flag,
        );
        if ava_cur > 0 {
            *return_ava_cur = ava_cur;
            *return_ava_ol = test_o_len;
            *return_chain_len = chain_len;
            for edge in chain_edges.iter() {
                circle_vis[(edge.ul >> 33) as usize] = 0;
                circle_vis[(edge.v >> 1) as usize] = 0;
            }
            return true;
        }

        end_rid = t.v;
        h = crate::assembly::rescue::get_best_contained_read(
            unitig_graph,
            read_graph,
            sources,
            coverage_cut,
            read_to_unitig,
            trio_labels,
            parameters,
            end_rid,
            ignore_trio_flag,
        );
    }

    for edge in chain_edges.iter() {
        circle_vis[(edge.ul >> 33) as usize] = 0;
        circle_vis[(edge.v >> 1) as usize] = 0;
    }
    false
}

/// `set_rtg_flag_by_bubble` (`Overlaps.cpp` 35321): flags the reads of
/// bubble `v` (body, beg and sink) in `vis_flag`.
fn set_rtg_flag_by_bubble(
    bub: &BubbleType,
    unitigs: &[MaUtg],
    read_graph: &Asg,
    v: u32,
    vis_flag: &mut [u8],
    flag: u32,
) {
    let (beg, sink, body, _) = get_bubbles(bub, v as u64);
    let mut flag_reads = |unitig_id: u32| {
        let chain = &unitigs[(unitig_id >> 1) as usize];
        if chain.reads.is_empty() {
            return;
        }
        for read in &chain.reads {
            let read_id = (read >> 33) as usize;
            if read_graph.sequences[read_id].deleted {
                continue;
            }
            vis_flag[read_id] = flag as u8;
        }
    };
    for &node in body {
        flag_reads(node);
    }
    if beg != u32::MAX {
        flag_reads(beg);
    }
    if sink != u32::MAX {
        flag_reads(sink);
    }
}

/// `minor_transitive_reduction_r_g` (`Overlaps.cpp` 35509): marks every
/// other pending edge redundant when its target is directly reachable.
fn minor_transitive_reduction_r_g(read_graph: &Asg, rbub_edges: &mut [AsgArc]) {
    for i in 0..rbub_edges.len() {
        let r_id = rbub_edges[i].v;
        for arc in read_graph.arcs_of(r_id) {
            if arc.deleted {
                continue;
            }
            let w = arc.v;
            for k in 0..rbub_edges.len() {
                if k == i {
                    continue;
                }
                if rbub_edges[k].v == w {
                    rbub_edges[k].deleted = true;
                }
            }
        }
    }
}

/// The shared tail of the three bubble rescues (`Overlaps.cpp` 35807-35853
/// and the copies in the missing-overlap rescues): flushes `new_edges` into
/// the read graph, re-labels the bubbles and deletes the new edges that sit
/// inside a bubble.
#[allow(clippy::too_many_arguments)]
fn finish_bubble_rescue(
    read_graph: &mut Asg,
    sources: &[PafList],
    coverage_cut: &mut [CoverageCut],
    bubble_labels: &mut BubbleLabel,
    trio_labels: &mut [u8],
    new_edges: &mut Vec<AsgArc>,
    minimum_length_threshold: i64,
) {
    for arc in new_edges.iter() {
        read_graph.push_arc(*arc);
    }
    if !new_edges.is_empty() {
        read_graph.cleanup();
        read_graph.symmetrize();
    }

    pre_clean(read_graph);
    lable_all_bubbles(read_graph, bubble_labels, minimum_length_threshold);

    for k in 0..new_edges.len() {
        let v = (new_edges[k].ul >> 32) as u32;
        let w = new_edges[k].v;
        if read_graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if read_graph.sequences[(w >> 1) as usize].deleted {
            continue;
        }
        // If this edge is inside a bubble, keep it.
        if read_graph.sequence_visited[v as usize] != 0
            && read_graph.sequence_visited[(w ^ 1) as usize] != 0
        {
            continue;
        }
        read_graph.delete_arc(v, w, true);
        read_graph.delete_arc(w ^ 1, v ^ 1, true);
    }
    read_graph.cleanup();
    read_graph.symmetrize();
    let _ = trio_labels;
}

/// `rescue_bubbles_by_contained_reads` (`Overlaps.cpp` 35620).
#[allow(clippy::too_many_arguments)]
pub fn rescue_bubbles_by_contained_reads(
    unitig_graph: &mut Asg,
    unitigs: &[MaUtg],
    read_graph: &mut Asg,
    sources: &[PafList],
    coverage_cut: &mut [CoverageCut],
    read_to_unitig: &mut ReadToUnitig,
    trio_labels: &mut [u8],
    bub: &BubbleType,
    bubble_labels: &mut BubbleLabel,
    parameters: crate::assembly::rescue::RescueParameters,
    chain_len_thres: u32,
    beg_idx: u32,
    occ: u32,
) {
    let chain_unitigs = bub.bubble_unitigs.as_ref().expect("chain unitigs built");
    let chain_graph = bub.bubble_chain_graph.as_ref().expect("chain graph built");
    for sequence in unitig_graph.sequences.iter_mut() {
        sequence.coverage_class = crate::assembly::hit::PRIMARY_LABEL;
    }

    let mut new_edges: Vec<AsgArc> = Vec::new();
    let mut rbub_edges: Vec<AsgArc> = Vec::new();
    let mut chain_edges: Vec<AsgArc> = Vec::new();
    let mut u_vecs: Vec<u64> = Vec::new();
    let mut expect_vis = vec![0u8; read_graph.sequences.len()];
    let mut circle_vis = vec![0u8; read_graph.sequences.len()];
    let mut utg_vis = vec![0u8; unitig_graph.sequences.len()];

    for v in 0..unitig_graph.sequences.len() as u32 {
        let unitig = &unitigs[v as usize];
        if unitig.reads.is_empty() {
            continue;
        }
        if unitig_graph.sequences[v as usize].deleted {
            continue;
        }
        for read in &unitig.reads {
            let read_id = (read >> 33) as u32;
            read_to_unitig.set(
                read_id,
                v,
                true,
                Some(read_graph.sequences[read_id as usize].coverage_class),
            );
        }
    }

    for chain in chain_unitigs.iter() {
        if chain.reads.len() < 3 {
            continue;
        }
        for bub_i in 1..chain.reads.len() - 1 {
            let m_bub = (chain.reads[bub_i] >> 33) as u64;
            if m_bub < beg_idx as u64 || m_bub >= beg_idx as u64 + occ as u64 {
                continue;
            }
            let l_bub = (chain.reads[bub_i - 1] >> 33) as u64;
            let r_bub = (chain.reads[bub_i + 1] >> 33) as u64;

            set_rtg_flag_by_bubble(bub, unitigs, read_graph, l_bub as u32, &mut expect_vis, 1);
            set_rtg_flag_by_bubble(bub, unitigs, read_graph, r_bub as u32, &mut expect_vis, 1);
            set_rtg_flag_by_bubble(bub, unitigs, read_graph, m_bub as u32, &mut expect_vis, 1);

            let (_, _, body, _) = get_bubbles(bub, m_bub);
            for &node in body {
                let u_id = (node >> 1) as usize;
                let unitig = &unitigs[u_id];
                if unitig.reads.is_empty() {
                    continue;
                }
                for k_v in 0..2u32 {
                    let v = ((u_id as u32) << 1) + k_v;
                    if get_real_length(unitig_graph, v) != 0 {
                        continue;
                    }

                    let end_rid = if v & 1 == 1 {
                        unitig.start ^ 1
                    } else {
                        unitig.end ^ 1
                    };

                    let hits = &sources[(end_rid >> 1) as usize].hits;
                    let mut ava_max = 0u32;
                    let mut ava_ol_max = 0u32;
                    let mut ava_min_chain = u32::MAX;
                    let mut best: Option<usize> = None;
                    let mut ava_cur = 0u32;
                    let mut test_o_len = 0u32;
                    let mut ava_chain_len = 0u32;
                    for k in 0..hits.len() {
                        if get_contained_reads_chain_by_broken_bub(
                            unitigs,
                            unitig_graph,
                            read_graph,
                            sources,
                            coverage_cut,
                            read_to_unitig,
                            trio_labels,
                            parameters,
                            end_rid,
                            &hits[k],
                            &mut chain_edges,
                            &mut ava_cur,
                            &mut test_o_len,
                            &mut ava_chain_len,
                            &expect_vis,
                            &mut circle_vis,
                            Some(&mut utg_vis),
                            1,
                        ) {
                            let mut is_update = false;
                            if ava_cur > ava_max {
                                is_update = true;
                            } else if ava_cur == ava_max {
                                if ava_chain_len < ava_min_chain {
                                    is_update = true;
                                } else if ava_chain_len == ava_min_chain && test_o_len > ava_ol_max {
                                    is_update = true;
                                }
                            }

                            if is_update {
                                ava_min_chain = ava_chain_len;
                                ava_max = ava_cur;
                                ava_ol_max = test_o_len;
                                best = Some(k);
                            }
                        }
                    }

                    if ava_max > 0 {
                        let k = best.expect("a best hit exists");
                        if !get_contained_reads_chain_by_broken_bub(
                            unitigs,
                            unitig_graph,
                            read_graph,
                            sources,
                            coverage_cut,
                            read_to_unitig,
                            trio_labels,
                            parameters,
                            end_rid,
                            &hits[k],
                            &mut chain_edges,
                            &mut ava_cur,
                            &mut test_o_len,
                            &mut ava_chain_len,
                            &expect_vis,
                            &mut circle_vis,
                            None,
                            1,
                        ) || chain_edges.is_empty()
                        {
                            continue;
                        }
                        let t_max = *chain_edges.last().expect("non-empty chain");

                        let mut rbub_edges: Vec<AsgArc> = Vec::new();
                        let mut scan_index = 0usize;
                        let mut r_edge = AsgArc::default();
                        while get_edge2existing_node_advance_by_broken_bub(
                            unitigs,
                            unitig_graph,
                            read_graph,
                            sources,
                            coverage_cut,
                            read_to_unitig,
                            trio_labels,
                            parameters,
                            &expect_vis,
                            t_max.v,
                            &mut scan_index,
                            &mut r_edge,
                            1,
                        ) {
                            rbub_edges.push(r_edge);
                        }

                        minor_transitive_reduction_r_g(read_graph, &mut rbub_edges);

                        for edge in chain_edges.iter() {
                            let t_max = *edge;
                            let (contain_rid, _) = read_to_unitig.get(t_max.v >> 1);
                            let a_nodes = ((contain_rid as u64) << 32) | (t_max.v >> 1) as u64;
                            u_vecs.push(a_nodes);

                            read_graph.sequences[(t_max.v >> 1) as usize].deleted = false;
                            coverage_cut[(t_max.v >> 1) as usize].deleted = false;
                            coverage_cut[(t_max.v >> 1) as usize].class =
                                crate::assembly::hit::PRIMARY_LABEL;

                            if let Some(t) = get_edge_from_source(
                                sources,
                                coverage_cut,
                                parameters.maximum_hang,
                                parameters.maximum_hang_rate,
                                parameters.minimum_overlap,
                                (t_max.ul >> 32) as u32,
                                t_max.v,
                            ) {
                                new_edges.push(t);
                            }
                            if let Some(t) = get_edge_from_source(
                                sources,
                                coverage_cut,
                                parameters.maximum_hang,
                                parameters.maximum_hang_rate,
                                parameters.minimum_overlap,
                                t_max.v ^ 1,
                                ((t_max.ul >> 32) as u32) ^ 1,
                            ) {
                                new_edges.push(t);
                            }
                        }

                        for edge in rbub_edges.iter() {
                            let t = *edge;
                            if t.deleted {
                                continue;
                            }
                            new_edges.push(t);
                            if let Some(t) = get_edge_from_source(
                                sources,
                                coverage_cut,
                                parameters.maximum_hang,
                                parameters.maximum_hang_rate,
                                parameters.minimum_overlap,
                                t.v ^ 1,
                                ((t.ul >> 32) as u32) ^ 1,
                            ) {
                                new_edges.push(t);
                            }
                        }
                    }
                }
            }

            set_rtg_flag_by_bubble(bub, unitigs, read_graph, l_bub as u32, &mut expect_vis, 0);
            set_rtg_flag_by_bubble(bub, unitigs, read_graph, r_bub as u32, &mut expect_vis, 0);
            set_rtg_flag_by_bubble(bub, unitigs, read_graph, m_bub as u32, &mut expect_vis, 0);
        }
    }
    let _ = chain_graph;

    finish_bubble_rescue(
        read_graph,
        sources,
        coverage_cut,
        bubble_labels,
        trio_labels,
        &mut new_edges,
        parameters.minimum_length_threshold,
    );

    for a_nodes in u_vecs.iter() {
        let mut w = (*a_nodes & 0xffffffff) as u32;
        w <<= 1;
        if !read_graph.sequences[(w >> 1) as usize].deleted
            && (get_real_length(read_graph, w) != 0 || get_real_length(read_graph, w ^ 1) != 0)
        {
            read_to_unitig.set_raw(w >> 1, u32::MAX);
        } else {
            read_graph.delete_sequence(w >> 1);
            coverage_cut[(w >> 1) as usize].deleted = true;
        }
    }

    for v in 0..read_to_unitig.len() {
        let (u_id, is_unitig) = read_to_unitig.get(v as u32);
        if is_unitig == 1 {
            read_to_unitig.set_raw(v as u32, u32::MAX);
        }
        let _ = u_id;
    }
}

/// `rescue_bubbles_by_missing_ovlp` (`Overlaps.cpp` 36092).
#[allow(clippy::too_many_arguments)]
pub fn rescue_bubbles_by_missing_ovlp(
    unitig_graph: &mut Asg,
    unitigs: &[MaUtg],
    read_graph: &mut Asg,
    sources: &[PafList],
    coverage_cut: &mut [CoverageCut],
    read_to_unitig: &mut ReadToUnitig,
    trio_labels: &mut [u8],
    bub: &BubbleType,
    bubble_labels: &mut BubbleLabel,
    parameters: RescueParameters,
    beg_idx: u32,
    occ: u32,
) {
    let chain_unitigs = bub.bubble_unitigs.as_ref().expect("chain unitigs built");
    for sequence in unitig_graph.sequences.iter_mut() {
        sequence.coverage_class = crate::assembly::hit::PRIMARY_LABEL;
    }

    let mut new_edges: Vec<AsgArc> = Vec::new();
    let mut expect_vis = vec![0u8; read_graph.sequences.len()];

    for v in 0..unitig_graph.sequences.len() as u32 {
        let unitig = &unitigs[v as usize];
        if unitig.reads.is_empty() {
            continue;
        }
        if unitig_graph.sequences[v as usize].deleted {
            continue;
        }
        for read in &unitig.reads {
            let read_id = (read >> 33) as u32;
            read_to_unitig.set(
                read_id,
                v,
                true,
                Some(read_graph.sequences[read_id as usize].coverage_class),
            );
        }
    }

    for chain in chain_unitigs.iter() {
        if chain.reads.len() < 3 {
            continue;
        }
        for bub_i in 1..chain.reads.len() - 1 {
            let m_bub = (chain.reads[bub_i] >> 33) as u64;
            if m_bub < beg_idx as u64 || m_bub >= beg_idx as u64 + occ as u64 {
                continue;
            }
            let l_bub = (chain.reads[bub_i - 1] >> 33) as u64;
            let r_bub = (chain.reads[bub_i + 1] >> 33) as u64;

            set_rtg_flag_by_bubble(bub, unitigs, read_graph, l_bub as u32, &mut expect_vis, 1);
            set_rtg_flag_by_bubble(bub, unitigs, read_graph, r_bub as u32, &mut expect_vis, 1);
            set_rtg_flag_by_bubble(bub, unitigs, read_graph, m_bub as u32, &mut expect_vis, 1);

            let (_, _, body, _) = get_bubbles(bub, m_bub);
            for &node in body {
                let u_id = (node >> 1) as usize;
                let unitig = &unitigs[u_id];
                if unitig.reads.is_empty() {
                    continue;
                }
                for k_v in 0..2u32 {
                    let v = ((u_id as u32) << 1) + k_v;
                    if get_real_length(unitig_graph, v) != 0 {
                        continue;
                    }

                    let end_rid = if v & 1 == 1 {
                        unitig.start ^ 1
                    } else {
                        unitig.end ^ 1
                    };

                    let mut rbub_edges: Vec<AsgArc> = Vec::new();
                    let mut scan_index = 0usize;
                    let mut r_edge = AsgArc::default();
                    while get_edge2existing_node_advance_by_broken_bub(
                        unitigs,
                        unitig_graph,
                        read_graph,
                        sources,
                        coverage_cut,
                        read_to_unitig,
                        trio_labels,
                        parameters,
                        &expect_vis,
                        end_rid,
                        &mut scan_index,
                        &mut r_edge,
                        1,
                    ) {
                        rbub_edges.push(r_edge);
                    }

                    if !rbub_edges.is_empty() {
                        minor_transitive_reduction_r_g(read_graph, &mut rbub_edges);
                        for edge in rbub_edges.iter() {
                            let t = *edge;
                            if t.deleted {
                                continue;
                            }
                            new_edges.push(t);
                            if let Some(t) = get_edge_from_source(
                                sources,
                                coverage_cut,
                                parameters.maximum_hang,
                                parameters.maximum_hang_rate,
                                parameters.minimum_overlap,
                                t.v ^ 1,
                                ((t.ul >> 32) as u32) ^ 1,
                            ) {
                                new_edges.push(t);
                            }
                        }
                    }
                }
            }

            set_rtg_flag_by_bubble(bub, unitigs, read_graph, l_bub as u32, &mut expect_vis, 0);
            set_rtg_flag_by_bubble(bub, unitigs, read_graph, r_bub as u32, &mut expect_vis, 0);
            set_rtg_flag_by_bubble(bub, unitigs, read_graph, m_bub as u32, &mut expect_vis, 0);
        }
    }

    finish_bubble_rescue(
        read_graph,
        sources,
        coverage_cut,
        bubble_labels,
        trio_labels,
        &mut new_edges,
        parameters.minimum_length_threshold,
    );

    for v in 0..read_to_unitig.len() {
        let (_, is_unitig) = read_to_unitig.get(v as u32);
        if is_unitig == 1 {
            read_to_unitig.set_raw(v as u32, u32::MAX);
        }
    }
}

/// `rescue_bubbles_by_missing_ovlp_backward` (`Overlaps.cpp` 36247).
#[allow(clippy::too_many_arguments)]
pub fn rescue_bubbles_by_missing_ovlp_backward(
    unitig_graph: &mut Asg,
    unitigs: &mut [MaUtg],
    read_graph: &mut Asg,
    sources: &[PafList],
    coverage_cut: &mut [CoverageCut],
    read_to_unitig: &mut ReadToUnitig,
    trio_labels: &mut [u8],
    bub: &BubbleType,
    bubble_labels: &mut BubbleLabel,
    parameters: RescueParameters,
    backward_steps: u32,
    beg_idx: u32,
    occ: u32,
) {
    let chain_unitigs = bub.bubble_unitigs.as_ref().expect("chain unitigs built");
    for sequence in unitig_graph.sequences.iter_mut() {
        sequence.coverage_class = crate::assembly::hit::PRIMARY_LABEL;
    }

    let mut new_edges: Vec<AsgArc> = Vec::new();
    let mut recover_edges: Vec<AsgArc> = Vec::new();
    let mut u_vecs: Vec<u64> = Vec::new();
    let mut expect_vis = vec![0u8; read_graph.sequences.len()];

    for v in 0..unitig_graph.sequences.len() as u32 {
        let unitig = &unitigs[v as usize];
        if unitig.reads.is_empty() {
            continue;
        }
        if unitig_graph.sequences[v as usize].deleted {
            continue;
        }
        for read in &unitig.reads {
            let read_id = (read >> 33) as u32;
            read_to_unitig.set(
                read_id,
                v,
                true,
                Some(read_graph.sequences[read_id as usize].coverage_class),
            );
        }
    }

    for chain in chain_unitigs.iter() {
        if chain.reads.len() < 3 {
            continue;
        }
        for bub_i in 1..chain.reads.len() - 1 {
            let m_bub = (chain.reads[bub_i] >> 33) as u64;
            if m_bub < beg_idx as u64 || m_bub >= beg_idx as u64 + occ as u64 {
                continue;
            }
            let l_bub = (chain.reads[bub_i - 1] >> 33) as u64;
            let r_bub = (chain.reads[bub_i + 1] >> 33) as u64;

            set_rtg_flag_by_bubble(bub, unitigs, read_graph, l_bub as u32, &mut expect_vis, 1);
            set_rtg_flag_by_bubble(bub, unitigs, read_graph, r_bub as u32, &mut expect_vis, 1);
            set_rtg_flag_by_bubble(bub, unitigs, read_graph, m_bub as u32, &mut expect_vis, 1);
            let (_, _, body, _) = get_bubbles(bub, m_bub);
            for &node in body {
                let u_id = (node >> 1) as usize;
                let unitig = &unitigs[u_id];
                if unitig.reads.is_empty() {
                    continue;
                }

                for k_v in 0..2u32 {
                    let v = ((u_id as u32) << 1) + k_v;
                    if get_real_length(unitig_graph, v) != 0 {
                        continue;
                    }

                    let (mut init, step, mode) = if v & 1 == 1 {
                        (0i64, 1i64, 1u32)
                    } else {
                        (unitig.reads.len() as i64 - 1, -1i64, 0u32)
                    };

                    let mut rbub_edges: Vec<AsgArc> = Vec::new();
                    let mut cur_backward_steps = backward_steps;
                    let mut round = 0u32;
                    while round < cur_backward_steps && init >= 0 && (init as usize) < unitig.reads.len() {
                        let mut end_rid = (unitig.reads[init as usize] >> 32) as u32;
                        end_rid ^= mode;

                        rbub_edges.clear();
                        let mut scan_index = 0usize;
                        let mut r_edge = AsgArc::default();
                        while get_edge2existing_node_advance_by_broken_bub(
                            unitigs,
                            unitig_graph,
                            read_graph,
                            sources,
                            coverage_cut,
                            read_to_unitig,
                            trio_labels,
                            parameters,
                            &expect_vis,
                            end_rid,
                            &mut scan_index,
                            &mut r_edge,
                            1,
                        ) {
                            rbub_edges.push(r_edge);
                        }

                        if !rbub_edges.is_empty() {
                            break;
                        }

                        init += step;
                        round += 1;
                        if rbub_edges.is_empty() && round < backward_steps {
                            let remaining = unitig.reads.len() as i64 - round as i64 - 1;
                            cur_backward_steps = if remaining > backward_steps as i64 {
                                backward_steps
                            } else {
                                remaining.max(0) as u32
                            };
                            let _ = cur_backward_steps;
                        }
                    }

                    if !rbub_edges.is_empty() {
                        // Save for the revert; the reference packs
                        // `mode<<31` then shifts the whole word left by 32,
                        // so `mode` lands on bit 63 — packing it at bit 31
                        // would corrupt the `uId` field.
                        let tmp = ((mode as u64) << 63)
                            | ((init as u64) & 0x7fffffff) << 32
                            | u_id as u64;
                        u_vecs.push(tmp);
                        minor_transitive_reduction_r_g(read_graph, &mut rbub_edges);

                        // Delete the reads between the hit and the unitig end.
                        let mut walk = init - step;
                        while walk >= 0 && (walk as usize) < unitig.reads.len() {
                            let w = (unitig.reads[walk as usize] >> 32) as u32;
                            for arc in read_graph.arcs_of(w) {
                                if arc.deleted {
                                    continue;
                                }
                                recover_edges.push(*arc);
                                match read_graph.get_arc(arc.v ^ 1, ((arc.ul >> 32) as u32) ^ 1) {
                                    Some(t) => recover_edges.push(t),
                                    None => eprintln!("error"),
                                }
                            }
                            for arc in read_graph.arcs_of(w ^ 1) {
                                if arc.deleted {
                                    continue;
                                }
                                recover_edges.push(*arc);
                                match read_graph.get_arc(arc.v ^ 1, ((arc.ul >> 32) as u32) ^ 1) {
                                    Some(t) => recover_edges.push(t),
                                    None => eprintln!("error"),
                                }
                            }

                            read_graph.delete_sequence(w >> 1);
                            expect_vis[(w >> 1) as usize] = 0;
                            walk -= step;
                        }

                        for edge in rbub_edges.iter() {
                            let t = *edge;
                            if t.deleted {
                                continue;
                            }
                            new_edges.push(t);
                            if let Some(t) = get_edge_from_source(
                                sources,
                                coverage_cut,
                                parameters.maximum_hang,
                                parameters.maximum_hang_rate,
                                parameters.minimum_overlap,
                                t.v ^ 1,
                                ((t.ul >> 32) as u32) ^ 1,
                            ) {
                                new_edges.push(t);
                            }
                        }
                    }
                }
            }

            set_rtg_flag_by_bubble(bub, unitigs, read_graph, l_bub as u32, &mut expect_vis, 0);
            set_rtg_flag_by_bubble(bub, unitigs, read_graph, r_bub as u32, &mut expect_vis, 0);
            set_rtg_flag_by_bubble(bub, unitigs, read_graph, m_bub as u32, &mut expect_vis, 0);
        }
    }

    finish_bubble_rescue(
        read_graph,
        sources,
        coverage_cut,
        bubble_labels,
        trio_labels,
        &mut new_edges,
        parameters.minimum_length_threshold,
    );

    for arc in recover_edges.iter_mut() {
        arc.deleted = true;
    }

    for entry in u_vecs.iter() {
        let mode = (entry >> 63) as u32;
        let step: i64 = if mode == 1 { 1 } else { -1 };
        let init = ((entry >> 32) & 0x7fffffff) as i64;
        let u_id = (entry & 0xffffffff) as usize;
        let unitig = &mut unitigs[u_id];

        let mut end_rid = (unitig.reads[init as usize] >> 32) as u32;
        end_rid ^= mode;
        if get_real_length(read_graph, end_rid) > 0 {
            crate::assembly::rescue::update_unitig(
                step,
                init,
                unitig,
                read_graph,
                Some(&mut recover_edges),
                1,
            );
        } else {
            crate::assembly::rescue::update_unitig(
                step,
                init,
                unitig,
                read_graph,
                Some(&mut recover_edges),
                0,
            );
        }
    }

    let mut recov_occ = 0u64;
    for arc in recover_edges.iter() {
        if arc.deleted {
            continue;
        }
        read_graph.push_arc(*arc);
        recov_occ += 1;
    }

    if recov_occ != 0 {
        read_graph.cleanup();
        read_graph.symmetrize();
    }

    for v in 0..read_to_unitig.len() {
        let (_, is_unitig) = read_to_unitig.get(v as u32);
        if is_unitig == 1 {
            read_to_unitig.set_raw(v as u32, u32::MAX);
        }
    }
}

/// `reset_bub` (`Overlaps.cpp` 36534): destroys and rebuilds the bubble
/// bookkeeping from the current unitig graph.
fn reset_bub(
    bub: &mut BubbleType,
    unitig_graph: &mut Asg,
    unitigs: &[MaUtg],
    ir_het: &[u8],
    new_rtg_edges: &mut Vec<AsgArc>,
) {
    *bub = BubbleType::default();
    new_rtg_edges.clear();
    identify_bubbles(unitig_graph, unitigs, bub, ir_het, None);
}

/// `rescue_bubble_by_chain` (`Overlaps.cpp` 36832): the bubble rescue
/// driver — classifies the bubbles through a primary-adjustment on graph
/// copies, then runs the three bubble rescues against the live read graph.
#[allow(clippy::too_many_arguments)]
pub fn rescue_bubble_by_chain(
    read_graph: &mut Asg,
    coverage_cut: &mut [CoverageCut],
    sources: &[PafList],
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    bubble_labels: &mut BubbleLabel,
    options: &mut crate::command_line::options::AssemblyOptions,
    trio_labels: &mut [u8],
    tips_len: i64,
    tip_drop_ratio: f32,
    stops_threshold: u32,
    chimeric_rate: f32,
    drop_ratio: f32,
    max_hang: i32,
    min_ovlp: i32,
    chain_len_thres: u32,
    minimum_length_threshold: i64,
) -> Result<(), String> {
    let parameters = RescueParameters {
        maximum_hang: max_hang,
        maximum_hang_rate: options.maximum_hang_rate,
        minimum_overlap: min_ovlp,
        minimum_length_threshold,
    };
    let mut new_rtg_edges: Vec<AsgArc> = Vec::new();
    let mut unitig_graph = crate::assembly::unitig::ma_ug_gen_primary(
        read_graph,
        crate::assembly::hit::PRIMARY_LABEL,
    );

    // `adjust_utg_by_primary` on copies produces the coverage object whose
    // `is_r_het` classifies the unitigs; the live coverage/labels are still
    // the ones mutated.
    let mut copy_sg = crate::assembly::unitig::copy_read_graph(read_graph);
    let mut copy_ug = crate::assembly::unitig::copy_untig_graph(&unitig_graph);
    let cov = crate::assembly::bp_graph::adjust_utg_by_primary(
        &mut copy_ug,
        &mut copy_sg,
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
        &new_rtg_edges,
        Some(&mut *bubble_labels),
        options,
        trio_labels,
        false,
        false,
        true,
    )?;
    let ir_het = cov.is_r_het;

    let mut bub = BubbleType::default();
    reset_bub(&mut bub, &mut unitig_graph.graph, &unitig_graph.unitigs, &ir_het, &mut new_rtg_edges);
    eprintln!(
        "[M::rescue_bubble_by_chain] f_bub: {}, s_bub: {}, b_bub: {}, b_end: {}, tangle: {}",
        bub.full_bubble_count,
        bub.simple_bubble_count,
        bub.broken_bubble_count,
        bub.broken_end_bubble_count,
        bub.tangle_bubble_count
    );
    rescue_bubbles_by_contained_reads(
        &mut unitig_graph.graph,
        &unitig_graph.unitigs,
        read_graph,
        sources,
        coverage_cut,
        read_to_unitig,
        trio_labels,
        &bub,
        bubble_labels,
        parameters,
        chain_len_thres,
        bub.full_bubble_count as u32,
        (bub.broken_bubble_count + bub.broken_end_bubble_count + bub.tangle_bubble_count) as u32,
    );

    unitig_graph = crate::assembly::unitig::ma_ug_gen_primary(
        read_graph,
        crate::assembly::hit::PRIMARY_LABEL,
    );
    reset_bub(&mut bub, &mut unitig_graph.graph, &unitig_graph.unitigs, &ir_het, &mut new_rtg_edges);
    rescue_bubbles_by_missing_ovlp(
        &mut unitig_graph.graph,
        &unitig_graph.unitigs,
        read_graph,
        sources,
        coverage_cut,
        read_to_unitig,
        trio_labels,
        &bub,
        bubble_labels,
        parameters,
        bub.full_bubble_count as u32,
        (bub.broken_bubble_count + bub.broken_end_bubble_count + bub.tangle_bubble_count) as u32,
    );

    unitig_graph = crate::assembly::unitig::ma_ug_gen_primary(
        read_graph,
        crate::assembly::hit::PRIMARY_LABEL,
    );
    reset_bub(&mut bub, &mut unitig_graph.graph, &unitig_graph.unitigs, &ir_het, &mut new_rtg_edges);
    rescue_bubbles_by_missing_ovlp_backward(
        &mut unitig_graph.graph,
        &mut unitig_graph.unitigs,
        read_graph,
        sources,
        coverage_cut,
        read_to_unitig,
        trio_labels,
        &bub,
        bubble_labels,
        parameters,
        10,
        bub.simple_bubble_count as u32,
        (bub.broken_bubble_count + bub.broken_end_bubble_count + bub.tangle_bubble_count) as u32,
    );

    Ok(())
}

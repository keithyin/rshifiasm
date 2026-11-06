//! The haplotype coverage state of the primary assembly, the counterpart of
//! `hap_cov_t`, `trans_chain`, `init_hap_cov_t`, `set_r_het_flag`,
//! `set_ug_coverage_aggressive`, `init_trans_chain`,
//! `asg_bub_backtrack_primary_cov`, `topologicalSortUtil` and
//! `dfs_trans_chain_bub` in `Overlaps.cpp`.
//!
//! `hap_cov_t::ir_het` and `hap_cov_t::is_r_het` are the same buffer in the
//! reference (`cov->t_ch->ir_het = cov->is_r_het`); here the flags live only in
//! [`HapCov::is_r_het`] and every writer takes that array directly.

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::bubble::BubbleBuffer;
use crate::assembly::bubble_identify::get_unitig;
use crate::assembly::gfa_clean::get_real_length;
use crate::assembly::hit::{CoverageCut, ReadToUnitig};
use crate::assembly::purge_dups::HOM_PEAK_RATE;
use crate::assembly::unitig::{UnitigGraph, FATHER, MOTHER};
use crate::command_line::options::AssemblyOptions;
use crate::correction::ec_driver::PafList;

/// `N_HET`: the read carries no heterozygosity evidence.
pub const N_HET: u8 = 0;
/// `C_HET`: the read looks heterozygous from the coverage layout.
pub const C_HET: u8 = 1;
/// `P_HET`: the read lies on the alternate side of a resolved bubble.
pub const P_HET: u8 = 2;

/// One interval of `bed_in`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BedInterval {
    pub beg: u64,
    pub end: u64,
}

/// `u_trans_t`: one read-to-read transition collected for the trio stage.
#[derive(Debug, Clone, Copy, Default)]
pub struct UTrans {
    pub query_start: u32,
    pub query_end: u32,
    pub query_identifier: u32,
    pub target_start: u32,
    pub target_end: u32,
    pub target_identifier: u32,
    pub occurrence: u32,
    pub network_weight: f64,
    /// `f` (6 bits).
    pub flag: u8,
    /// `rev`.
    pub reverse: bool,
    /// `del`.
    pub deleted: bool,
}

/// `u_trans_hit_t`: the four candidate intervals of one transition.
#[derive(Debug, Clone, Copy, Default)]
pub struct UTransHit {
    pub query_start_pre: u32,
    pub query_end_pre: u32,
    pub query_start_cur: u32,
    pub query_end_cur: u32,
    pub query_identifier: u32,
    pub target_start_pre: u32,
    pub target_end_pre: u32,
    pub target_start_cur: u32,
    pub target_end_cur: u32,
    pub target_identifier: u32,
}

/// `ca_buf_t`: one candidate transition pair of the chain builder.
#[derive(Debug, Clone, Copy, Default)]
pub struct CABuf {
    pub p_x_p: u32,
    pub p_y_p: u32,
    pub p_x: u32,
    pub p_y: u32,
    pub c_x_p: u32,
    pub c_y_p: u32,
    pub c_reverse: bool,
}

/// `sub_tran_t`: the chain buffer header.
#[derive(Debug, Clone, Default)]
pub struct SubTran {
    pub u_ids: Vec<u32>,
    pub i_dxs: Vec<u32>,
    pub chain_num: u32,
}

/// `trans_chain`: the read-to-read transition state accumulated while the
/// primary graph is cleaned and consumed by the trio output stage.
///
/// `ir_het` aliases [`HapCov::is_r_het`] in the reference and is therefore not
/// duplicated here. `k_trans`/`k_t_b`/`c_buf` are filled by
/// `chain_origin_trans_uid_*`; those writers are deferred until the trio stage
/// is ported, so the vectors stay empty for now.
#[derive(Debug, Clone, Default)]
pub struct TransChain {
    /// `rUidx`: the vertex of the unitig holding each read.
    pub r_u_idx: Vec<u32>,
    /// `rUpos`: the packed `<offset>:32 | read index>` position of each read.
    pub r_u_pos: Vec<u64>,
    /// `r_num`.
    pub r_num: u32,
    /// `u_num`.
    pub u_num: u32,
    /// `bed`: the per-unitig interval lists.
    pub bed: Vec<Vec<BedInterval>>,
    /// `topo_buf`.
    pub topo_buf: Vec<u32>,
    /// `topo_res`: the topological order of the bubble interior.
    pub topo_res: Vec<u32>,
    /// `b_buf_0`.
    pub b_buf_0: BubbleBuffer,
    /// `b_buf_1`.
    pub b_buf_1: BubbleBuffer,
    /// `k_trans`.
    pub k_trans: Vec<UTrans>,
    /// `k_t_b`.
    pub k_t_b: Vec<UTransHit>,
    /// `c_buf`.
    pub c_buf: Vec<CABuf>,
    /// `st`.
    pub st: SubTran,
}

/// `asg_arc_t_offset`: one read-to-unitig transition of the haplotype
/// alignment chaining, `asg_arc_t_offset` in `Overlaps.h`.
///
/// `offset` is `<x position>:32 | <y position>:32` along the two chains the
/// transition connects.
#[derive(Debug, Clone, Copy, Default)]
pub struct ArcOffset {
    pub arc: AsgArc,
    /// `Off`.
    pub offset: u64,
    /// `weight`.
    pub weight: u64,
}

/// `hap_cov_t`: per-read coverage, the heterozygosity flags and the transition
/// chain of the primary unitig graph.
#[derive(Debug, Clone, Default)]
pub struct HapCov {
    /// `n`: the number of reads.
    pub n: u32,
    /// `cov[rId]`: the read's coverage inside its unitig.
    pub cov: Vec<u32>,
    /// `pos_idx`: `topologicalSortUtil`'s visit flags, `u64::MAX` when unset.
    pub pos_idx: Vec<u64>,
    /// `read_g`'s `seq[].len`, which the bubble backtracking reads.
    pub read_lengths: Vec<u32>,
    /// `max_hang`, kept for the callers that mirror the reference signature.
    pub max_hang: i32,
    /// `min_ovlp`, kept for the callers that mirror the reference signature.
    pub min_ovlp: i32,
    /// `u_buffer`: the chaining scratch buffer.
    pub u_buffer: Vec<ArcOffset>,
    /// `tailIndex`.
    pub tail_index: Vec<i32>,
    /// `prevIndex`.
    pub prev_index: Vec<i32>,
    /// `is_r_het` (also `trans_chain::ir_het`): the `C_HET`/`P_HET` flags.
    pub is_r_het: Vec<u8>,
    /// `t_ch`: present when the pipeline collects read-to-read transitions.
    pub t_ch: Option<TransChain>,
}

/// `get_num_trio_flag`: how many reads of unitig `v` carry `flag`.
///
/// `trio_labels` is `R_INF.trio_flag`; `flag == (uint32_t)-1` means "no flag".
pub fn get_num_trio_flag(
    unitigs: &[MaUtg],
    graph: &Asg,
    v: u32,
    flag: u32,
    trio_labels: &[u8],
) -> u32 {
    if flag == u32::MAX {
        return 0;
    }
    if graph.sequences[v as usize].deleted {
        return 0;
    }
    let unitig = &unitigs[v as usize];
    if unitig.reads.is_empty() {
        return 0;
    }

    let mut occurrence = 0u32;
    for &packed in &unitig.reads {
        let read = (packed >> 33) as usize;
        if u32::from(trio_labels[read]) == flag {
            occurrence += 1;
        }
    }
    occurrence
}

/// `init_hap_cov_t`: builds the coverage array and the heterozygosity flags.
///
/// The `reverse_sources` argument of the reference is only stored and first
/// read by later stages, so it stays with the functions that consume it.
#[allow(clippy::too_many_arguments)]
pub fn init_hap_cov_t(
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    coverage_cut: &[CoverageCut],
    options: &AssemblyOptions,
    is_collect_trans: bool,
) -> HapCov {
    let read_count = read_graph.sequences.len();
    let position_index = vec![u64::MAX; read_count];
    let mut set = vec![0u8; read_count << 1];
    let mut cov = vec![0u32; read_count];

    let graph = &unitig_graph.graph;
    let unitigs = &unitig_graph.unitigs;
    for i in 0..graph.sequences.len() {
        if graph.sequences[i].deleted {
            continue;
        }
        let unitig = &unitigs[i];
        for &packed in &unitig.reads {
            set[(packed >> 33) as usize] = 1;
        }

        for orientation in 0..2u32 {
            let v = ((i as u32) << 1) + orientation;
            for arc in graph.arcs_of(v) {
                if arc.deleted {
                    continue;
                }
                let w = (arc.v >> 1) as usize;
                if graph.sequences[w].deleted {
                    continue;
                }
                for &packed in &unitigs[w].reads {
                    set[(packed >> 33) as usize] = 1;
                }
            }
        }

        let unitig = &unitigs[i];
        for &packed in &unitig.reads {
            let read = (packed >> 33) as usize;
            let mut c_bases = 0u32;
            for hit in &sources[read].hits {
                let mut target = hit.target_identifier;
                // `get_R_to_U` rewrites a deleted target to its containing
                // read; the reference skips targets that stay deleted.
                if read_graph.sequences[target as usize].deleted {
                    let (mapped, is_unitig) = read_to_unitig.get(target);
                    target = mapped;
                    if target == u32::MAX
                        || is_unitig == 1
                        || read_graph.sequences[target as usize].deleted
                    {
                        continue;
                    }
                }
                if read_graph.sequences[target as usize].deleted {
                    continue;
                }
                if set[target as usize] == 0 {
                    continue;
                }
                c_bases = c_bases.wrapping_add(hit.query_end.wrapping_sub(hit.query_start));
            }
            cov[read] = cov[read].max(c_bases);
        }

        for &packed in &unitig.reads {
            set[(packed >> 33) as usize] = 0;
        }
        for orientation in 0..2u32 {
            let v = ((i as u32) << 1) + orientation;
            for arc in graph.arcs_of(v) {
                if arc.deleted {
                    continue;
                }
                let w = (arc.v >> 1) as usize;
                if graph.sequences[w].deleted {
                    continue;
                }
                for &packed in &unitigs[w].reads {
                    set[(packed >> 33) as usize] = 0;
                }
            }
        }
    }

    let mut is_r_het = vec![0u8; read_count];
    set_r_het_flag(
        unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        read_to_unitig,
        &mut is_r_het,
        options,
    );

    let trans_chain = if is_collect_trans {
        Some(init_trans_chain(unitig_graph, read_count))
    } else {
        None
    };

    HapCov {
        n: read_count as u32,
        cov,
        pos_idx: position_index,
        read_lengths: read_graph.sequences.iter().map(|seq| seq.length).collect(),
        max_hang: options.maximum_hang_length,
        min_ovlp: options.minimum_overlap_length,
        u_buffer: Vec::new(),
        tail_index: Vec::new(),
        prev_index: Vec::new(),
        is_r_het,
        t_ch: trans_chain,
    }
}

/// `set_r_het_flag`: marks every read whose coverage looks heterozygous.
#[allow(clippy::too_many_arguments)]
pub fn set_r_het_flag(
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
    coverage_cut: &[CoverageCut],
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    is_r_het: &mut [u8],
    options: &AssemblyOptions,
) {
    let mut read_flag = vec![0u8; read_graph.sequences.len()];

    let mut dip_thre_max = if options.global_homozygous_coverage_is_set != 0 {
        options.global_homozygous_coverage as u64
    } else {
        // The reference converts the double to `uint64_t`; the negative
        // default (`-1 / HOM_PEAK_RATE`) therefore truncates to zero.
        (f64::from(options.global_homozygous_coverage) / HOM_PEAK_RATE) as u64
    };
    dip_thre_max = ((dip_thre_max as f64)
        - (((dip_thre_max as f64) * 0.5) / f64::from(options.polyploidy)))
        as u64;

    for m in 0..unitig_graph.graph.sequences.len() {
        set_ug_coverage_aggressive(
            unitig_graph,
            m as u32,
            read_graph,
            coverage_cut,
            sources,
            read_to_unitig,
            &mut read_flag,
            is_r_het,
            dip_thre_max as i64,
        );
    }
}

/// `set_ug_coverage_aggressive`: flags the reads of one unitig that are
/// heterozygous by their in-unitig versus out-of-unitig exact overlaps.
#[allow(clippy::too_many_arguments)]
fn set_ug_coverage_aggressive(
    unitig_graph: &UnitigGraph,
    unitig_id: u32,
    read_graph: &Asg,
    coverage_cut: &[CoverageCut],
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    read_flag: &mut [u8],
    is_r_het: &mut [u8],
    het_coverage_threshold: i64,
) {
    let unitigs = &unitig_graph.unitigs;
    let graph = &unitig_graph.graph;
    let unitig = &unitigs[unitig_id as usize];
    if unitig.reads.is_empty() {
        return;
    }

    for &packed in &unitig.reads {
        read_flag[(packed >> 33) as usize] = 1;
    }
    // The reference walks both orientations and does not test the arc's
    // deletion flag here.
    for orientation in 0..2u32 {
        let v = (unitig_id << 1) + orientation;
        for arc in graph.arcs_of(v) {
            let neighbour = &unitigs[(arc.v >> 1) as usize];
            for &packed in &neighbour.reads {
                read_flag[(packed >> 33) as usize] = 2;
            }
        }
    }

    let unitig = &unitigs[unitig_id as usize];
    for &packed in &unitig.reads {
        let read = (packed >> 33) as usize;
        let mut cov_in = 0i64;
        let mut cov_out = 0i64;
        let mut cov_in_start = 0i64;
        let mut cov_in_end = 0i64;
        let mut cov_out_start = 0i64;
        let mut cov_out_end = 0i64;

        for hit in &sources[read].hits {
            if !hit.exact_region {
                continue;
            }
            let mut target = hit.target_identifier;
            if read_graph.sequences[target as usize].deleted {
                let (mapped, is_unitig) = read_to_unitig.get(target);
                target = mapped;
                if target == u32::MAX
                    || is_unitig == 1
                    || read_graph.sequences[target as usize].deleted
                {
                    continue;
                }
            }
            if read_graph.sequences[target as usize].deleted {
                continue;
            }
            if read_flag[target as usize] == 0 {
                continue;
            }
            let overlap_bases = i64::from(hit.query_end.wrapping_sub(hit.query_start));
            if read_flag[target as usize] == 1 {
                cov_in += overlap_bases;
                if hit.query_start <= coverage_cut[read].start {
                    cov_in_start += 1;
                }
                if hit.query_end >= coverage_cut[read].end {
                    cov_in_end += 1;
                }
            }
            if read_flag[target as usize] == 2 {
                cov_out += overlap_bases;
                if hit.query_start <= coverage_cut[read].start {
                    cov_out_start += 1;
                }
                if hit.query_end >= coverage_cut[read].end {
                    cov_out_end += 1;
                }
            }
        }

        let (central_cov, try_cov) = if cov_out_start >= cov_out_end {
            (cov_in_end, cov_in_start + cov_out_start)
        } else {
            (cov_in_start, cov_in_end + cov_out_end)
        };

        let (c_bases, r_bases) = if (cov_out as f64) <= (cov_in as f64) * 0.1 {
            (
                cov_in + cov_out,
                i64::from(coverage_cut[read].end.wrapping_sub(coverage_cut[read].start)),
            )
        } else if (try_cov as f64) <= (central_cov as f64) * 1.1
            && cov_out_start.min(cov_out_end) as f64
                <= (cov_out_start.max(cov_out_end) as f64) * 0.2
        {
            (
                cov_in + cov_out,
                i64::from(coverage_cut[read].end.wrapping_sub(coverage_cut[read].start)),
            )
        } else {
            (central_cov, 1)
        };

        if r_bases <= 0 || (c_bases / r_bases) <= het_coverage_threshold {
            is_r_het[read] |= C_HET;
        }
    }

    for &packed in &unitig.reads {
        read_flag[(packed >> 33) as usize] = 0;
    }
    for orientation in 0..2u32 {
        let v = (unitig_id << 1) + orientation;
        for arc in graph.arcs_of(v) {
            let neighbour = &unitigs[(arc.v >> 1) as usize];
            for &packed in &neighbour.reads {
                read_flag[(packed >> 33) as usize] = 0;
            }
        }
    }
}

/// `init_trans_chain`: allocates the transition chain and indexes every read
/// inside its unitig.
pub fn init_trans_chain(unitig_graph: &UnitigGraph, read_count: usize) -> TransChain {
    let graph = &unitig_graph.graph;
    let mut r_u_idx = vec![u32::MAX; read_count];
    let mut r_u_pos = vec![u64::MAX; read_count];
    let mut is_duplicate = false;

    for v in 0..graph.sequences.len() {
        if graph.sequences[v].deleted {
            continue;
        }
        let unitig = &unitig_graph.unitigs[v];
        if unitig.reads.is_empty() {
            continue;
        }
        let mut offset = 0u64;
        for (k, &packed) in unitig.reads.iter().enumerate() {
            let read = (packed >> 33) as usize;

            let vertex = ((v as u32) << 1) | ((packed >> 32) & 1) as u32;
            if r_u_idx[read] != u32::MAX && r_u_idx[read] != vertex {
                is_duplicate = true;
            }
            r_u_idx[read] = vertex;

            let position = (offset << 32) | k as u64;
            if r_u_pos[read] != u64::MAX && r_u_pos[read] != position {
                is_duplicate = true;
            }
            r_u_pos[read] = position;

            offset += u64::from(packed as u32);
        }
    }

    if is_duplicate {
        for v in 0..graph.sequences.len() {
            if graph.sequences[v].deleted {
                continue;
            }
            let unitig = &unitig_graph.unitigs[v];
            if unitig.reads.is_empty() {
                continue;
            }
            let mut offset = 0u64;
            for (k, &packed) in unitig.reads.iter().enumerate() {
                let read = (packed >> 33) as usize;
                let mut duplicate = false;

                let vertex = ((v as u32) << 1) | ((packed >> 32) & 1) as u32;
                if r_u_idx[read] != vertex {
                    duplicate = true;
                }
                let position = (offset << 32) | k as u64;
                if r_u_pos[read] != position {
                    duplicate = true;
                }

                if duplicate {
                    r_u_idx[read] = u32::MAX;
                    r_u_pos[read] = u64::MAX;
                }
                offset += u64::from(packed as u32);
            }
        }
    }

    let unitig_count = graph.sequences.len();
    let mut bed = Vec::with_capacity(unitig_count);
    for _ in 0..unitig_count {
        bed.push(Vec::new());
    }
    let mut state = SubTran::default();
    state.i_dxs.push(0);

    TransChain {
        r_u_idx,
        r_u_pos,
        r_num: read_count as u32,
        u_num: unitig_count as u32,
        bed,
        topo_buf: Vec::new(),
        topo_res: Vec::new(),
        b_buf_0: BubbleBuffer::default(),
        b_buf_1: BubbleBuffer::default(),
        k_trans: Vec::new(),
        k_t_b: Vec::new(),
        c_buf: Vec::new(),
        st: state,
    }
}

/// `reset_trans_chain`: clears the heterozygosity flag of every read of `u`.
///
/// The reference's `t_ch->ir_het` aliases [`HapCov::is_r_het`], so the caller
/// passes that array directly. `u->n == 0 || u->m == 0` is modelled by an
/// empty read list: every unitig built by `ma_ug_gen` or copied by
/// `copy_untig_graph` has a nonzero capacity exactly when it has reads.
pub fn reset_trans_chain(is_r_het: &mut [u8], unitig: &MaUtg) {
    if unitig.reads.is_empty() {
        return;
    }
    for &packed in &unitig.reads {
        is_r_het[(packed >> 33) as usize] = N_HET;
    }
}

/// `topologicalSortUtil`: fills `topo_res` with the bubble interior in
/// topological order, beginning at `beg` and ending before `sink`.
fn topological_sort_util_impl(
    graph: &Asg,
    pos_idx: &mut [u64],
    topo_res: &mut Vec<u32>,
    b_buf_0: &mut BubbleBuffer,
    beg: u32,
    sink: u32,
) {
    b_buf_0.visited.clear();
    topo_res.clear();
    b_buf_0.visited.push(beg);

    while !b_buf_0.visited.is_empty() {
        let v = *b_buf_0.visited.last().expect("the loop keeps the stack non-empty");
        if pos_idx[(v >> 1) as usize] == u64::MAX {
            pos_idx[(v >> 1) as usize] = 0;
        }

        let mut kv = 0u32;
        for arc in graph.arcs_of(v) {
            if arc.deleted {
                continue;
            }
            if (arc.v >> 1) == (beg >> 1) || (arc.v >> 1) == (sink >> 1) {
                continue;
            }
            if pos_idx[(arc.v >> 1) as usize] != u64::MAX {
                continue;
            }
            b_buf_0.visited.push(arc.v);
            kv += 1;
        }

        if kv != 0 {
            continue;
        }
        b_buf_0.visited.pop();
        if pos_idx[(v >> 1) as usize] != 1 {
            topo_res.push(v);
            pos_idx[(v >> 1) as usize] = 1;
        }
    }

    for &node in topo_res.iter() {
        pos_idx[(node >> 1) as usize] = u64::MAX;
    }

    // Drop `beg`, then reverse the remaining order.
    topo_res.pop();
    let count = topo_res.len();
    for i in 0..(count >> 1) {
        topo_res.swap(i, count - i - 1);
    }
}

/// `topologicalSortUtil` over a [`HapCov`]'s trans chain.
pub fn topological_sort_util(graph: &Asg, cov: &mut HapCov, beg: u32, sink: u32) {
    let HapCov { pos_idx, t_ch, .. } = cov;
    let t_ch = t_ch
        .as_mut()
        .expect("topologicalSortUtil requires a transition chain");
    let TransChain {
        topo_res, b_buf_0, ..
    } = t_ch;
    topological_sort_util_impl(graph, pos_idx, topo_res, b_buf_0, beg, sink);
}

/// `dfs_trans_chain_bub`: collects the topologically sorted bubble nodes that
/// are *not* reachable from either orientation of `v` into `b_buf_0`.
fn dfs_trans_chain_bub_impl(
    graph: &Asg,
    pos_idx: &mut [u64],
    topo_res: &[u32],
    b_buf_0: &mut BubbleBuffer,
    v: u32,
    beg: u32,
    sink: u32,
) {
    b_buf_0.visited.clear();
    if v == beg || v == sink {
        return;
    }

    let mut start = v << 1;
    b_buf_0.visited.push(start);
    while let Some(current) = b_buf_0.visited.pop() {
        if pos_idx[(current >> 1) as usize] == 0 && (current >> 1) != (start >> 1) {
            continue;
        }
        pos_idx[(current >> 1) as usize] = 0;

        for arc in graph.arcs_of(current) {
            if arc.deleted {
                continue;
            }
            if (arc.v >> 1) == beg || (arc.v >> 1) == sink {
                continue;
            }
            if pos_idx[(arc.v >> 1) as usize] == 0 {
                continue;
            }
            b_buf_0.visited.push(arc.v);
        }
    }

    start += 1;
    b_buf_0.visited.push(start);
    while let Some(current) = b_buf_0.visited.pop() {
        if pos_idx[(current >> 1) as usize] == 0 && (current >> 1) != (start >> 1) {
            continue;
        }
        pos_idx[(current >> 1) as usize] = 0;

        for arc in graph.arcs_of(current) {
            if arc.deleted {
                continue;
            }
            if (arc.v >> 1) == beg || (arc.v >> 1) == sink {
                continue;
            }
            if pos_idx[(arc.v >> 1) as usize] == 0 {
                continue;
            }
            b_buf_0.visited.push(arc.v);
        }
    }

    b_buf_0.visited.clear();
    for &node in topo_res {
        if pos_idx[(node >> 1) as usize] == 0 {
            pos_idx[(node >> 1) as usize] = u64::MAX;
        } else {
            b_buf_0.visited.push(node);
        }
    }
}

/// `dfs_trans_chain_bub` over a [`HapCov`]'s trans chain.
pub fn dfs_trans_chain_bub(graph: &Asg, cov: &mut HapCov, v: u32, beg: u32, sink: u32) {
    let HapCov { pos_idx, t_ch, .. } = cov;
    let t_ch = t_ch
        .as_mut()
        .expect("dfs_trans_chain_bub requires a transition chain");
    let TransChain {
        topo_res, b_buf_0, ..
    } = t_ch;
    dfs_trans_chain_bub_impl(graph, pos_idx, topo_res, b_buf_0, v, beg, sink);
}

/// `asg_bub_backtrack_primary_cov`: redistributes the bubble coverage and
/// marks the alternate side's reads with `P_HET`.
///
/// The two `chain_origin_trans_uid_*` calls of the reference only append to
/// `t_ch->k_trans`/`k_t_b`/`c_buf` and read `pos_idx` without writing it
/// (`set_utg_offset(..., only_len = 1)`); those transitions are consumed by
/// the trio output stage, so their construction is deferred here and does not
/// affect any primary output.
///
/// The caller must pass the bubble buffer of the successful pop (its
/// `sources[0]` is the sink and its `visited` the bubble nodes), before
/// `asg_bub_backtrack_primary` deletes the bubble's arcs.
pub fn asg_bub_backtrack_primary_cov(
    graph: &Asg,
    unitigs: &[MaUtg],
    v0: u32,
    buffer: &BubbleBuffer,
    cov: &mut HapCov,
    is_update_chain: bool,
) {
    let sink = buffer.sources[0];

    // Sum the coverage of every visited unitig except the sink's.
    let mut unitig_cov = 0u32;
    for &node in &buffer.visited {
        let unitig_id = (node >> 1) as usize;
        if unitig_id == (sink >> 1) as usize {
            continue;
        }
        for &packed in &unitigs[unitig_id].reads {
            unitig_cov = unitig_cov.wrapping_add(cov.cov[(packed >> 33) as usize]);
        }
    }

    // Subtract the coverage of the reads on the sink-to-beg path and measure
    // the path length.
    let mut unitig_len = 0u32;
    let mut v = sink;
    loop {
        let u = buffer.info[v as usize].parent;
        if v != sink {
            let unitig_id = (v >> 1) as usize;
            for &packed in &unitigs[unitig_id].reads {
                let read = (packed >> 33) as usize;
                unitig_cov = unitig_cov.wrapping_sub(cov.cov[read]);
                unitig_len = unitig_len.wrapping_add(cov.read_lengths[read]);
            }
        }
        v = u;
        if v == v0 {
            break;
        }
    }

    let unitig_cov = if unitig_len == 0 {
        0
    } else {
        unitig_cov / unitig_len
    };

    // Hand that average coverage back to every read on the path.
    let mut v = sink;
    loop {
        let u = buffer.info[v as usize].parent;
        if v != sink {
            let unitig_id = (v >> 1) as usize;
            for &packed in &unitigs[unitig_id].reads {
                let read = (packed >> 33) as usize;
                cov.cov[read] = cov.cov[read]
                    .wrapping_add(unitig_cov.wrapping_mul(cov.read_lengths[read]));
            }
        }
        v = u;
        if v == v0 {
            break;
        }
    }

    if !is_update_chain || cov.t_ch.is_none() {
        return;
    }

    if asg_bub_backtrack_primary_cov_simple_bubble(graph, unitigs, v0, buffer, cov) {
        return;
    }

    let HapCov {
        pos_idx,
        is_r_het,
        t_ch,
        ..
    } = cov;
    let t_ch = t_ch.as_mut().expect("checked above");
    let TransChain {
        topo_res,
        b_buf_0,
        ..
    } = t_ch;

    topological_sort_util_impl(graph, pos_idx, topo_res, b_buf_0, v0, sink);
    if topo_res.is_empty() {
        return;
    }

    for index in 0..topo_res.len() {
        let node = topo_res[index];
        let unitig_id = (node >> 1) as usize;
        if unitig_id == (sink >> 1) as usize {
            continue;
        }

        dfs_trans_chain_bub_impl(
            graph,
            pos_idx,
            topo_res,
            b_buf_0,
            node >> 1,
            v0 >> 1,
            sink >> 1,
        );

        if b_buf_0.visited.is_empty() {
            continue;
        }

        // Deferred: `chain_origin_trans_uid_c_bubble(node, &t_ch->b_buf_0,
        // b, ug, cov)` only appends to `k_trans`/`k_t_b`/`c_buf`.

        mark_unitig_reads_heterozygous(unitigs, node, is_r_het);
        for k_i in 0..b_buf_0.visited.len() {
            let other = b_buf_0.visited[k_i];
            mark_unitig_reads_heterozygous(unitigs, other, is_r_het);
        }
    }
}

/// The simple two-side branch of [`asg_bub_backtrack_primary_cov`].
///
/// Returns `true` when the bubble was recognised as a simple bubble; the
/// reference then marks both sides and returns before the generic
/// topological walk. The two `chain_origin_trans_uid_s_bubble` calls that
/// follow the marking are deferred with the other transition builders.
fn asg_bub_backtrack_primary_cov_simple_bubble(
    graph: &Asg,
    unitigs: &[MaUtg],
    v0: u32,
    buffer: &BubbleBuffer,
    cov: &mut HapCov,
) -> bool {
    let sink = buffer.sources[0];
    if get_real_length(graph, v0) != 2 || get_real_length(graph, sink ^ 1) != 2 {
        return false;
    }

    let mut branches = [0u32; 2];
    get_live_targets(graph, v0, &mut branches);
    if (branches[0] >> 1) == (sink >> 1) || (branches[1] >> 1) == (sink >> 1) {
        return false;
    }

    let (mut end_0, mut end_1) = (0u32, 0u32);
    let (mut occ_0, mut occ_1) = (0i64, 0i64);
    let (mut base_0, mut base_1) = (0i64, 0i64);
    let (mut max_stop_node, mut max_stop_base) = (0i64, 0i64);

    get_unitig(
        graph,
        None,
        branches[0],
        &mut end_0,
        &mut occ_0,
        &mut base_0,
        &mut max_stop_node,
        &mut max_stop_base,
        1,
        None,
    );
    get_unitig(
        graph,
        None,
        branches[1],
        &mut end_1,
        &mut occ_1,
        &mut base_1,
        &mut max_stop_node,
        &mut max_stop_base,
        1,
        None,
    );

    if occ_0 + occ_1 + 1 != buffer.visited.len() as i64
        || get_real_length(graph, end_0) != 1
        || get_real_length(graph, end_1) != 1
    {
        return false;
    }

    let mut targets = [0u32; 1];
    get_live_targets(graph, end_0, &mut targets);
    end_0 = targets[0];
    get_live_targets(graph, end_1, &mut targets);
    end_1 = targets[0];
    if end_0 != sink || end_1 != sink {
        return false;
    }

    let HapCov {
        is_r_het, t_ch, ..
    } = cov;
    let t_ch = t_ch.as_mut().expect("the caller checked the chain");
    let TransChain {
        b_buf_0, b_buf_1, ..
    } = t_ch;

    b_buf_0.visited.clear();
    let (mut end, mut occ, mut base) = (0u32, 0i64, 0i64);
    get_unitig(
        graph,
        None,
        branches[0],
        &mut end,
        &mut occ,
        &mut base,
        &mut max_stop_node,
        &mut max_stop_base,
        1,
        Some(&mut b_buf_0.visited),
    );
    for i in 0..b_buf_0.visited.len() {
        let node = b_buf_0.visited[i];
        mark_unitig_reads_heterozygous(unitigs, node, is_r_het);
    }

    b_buf_1.visited.clear();
    let (mut end, mut occ, mut base) = (0u32, 0i64, 0i64);
    get_unitig(
        graph,
        None,
        branches[1],
        &mut end,
        &mut occ,
        &mut base,
        &mut max_stop_node,
        &mut max_stop_base,
        1,
        Some(&mut b_buf_1.visited),
    );
    for i in 0..b_buf_1.visited.len() {
        let node = b_buf_1.visited[i];
        mark_unitig_reads_heterozygous(unitigs, node, is_r_het);
    }

    // Deferred: `chain_origin_trans_uid_s_bubble(&t_ch->b_buf_0,
    // &t_ch->b_buf_1, v0, b->S.a[0] ^ 1, ug, cov)`.

    true
}

/// Marks every read of the unitig behind `node` with `P_HET`, walking the
/// unitig's read list in the orientation the vertex encodes.
fn mark_unitig_reads_heterozygous(unitigs: &[MaUtg], node: u32, is_r_het: &mut [u8]) {
    let unitig = &unitigs[(node >> 1) as usize];
    if unitig.reads.is_empty() {
        return;
    }
    let orientation = node & 1;
    for k in 0..unitig.reads.len() {
        let packed = if orientation == 1 {
            unitig.reads[unitig.reads.len() - k - 1]
        } else {
            unitig.reads[k]
        };
        is_r_het[(packed >> 33) as usize] |= P_HET;
    }
}

/// `get_real_length(g, v, v_s)`: the live out-arcs of `v`, copied into `out`.
fn get_live_targets(graph: &Asg, v: u32, out: &mut [u32]) -> u32 {
    let mut count = 0usize;
    for arc in graph.arcs_of(v) {
        if arc.deleted {
            continue;
        }
        if count < out.len() {
            out[count] = arc.v;
        }
        count += 1;
    }
    count as u32
}

/// The `FATHER`/`MOTHER` trio labels in the `u32` form the pop flags use.
pub fn parent_flags() -> (u32, u32) {
    (u32::from(FATHER), u32::from(MOTHER))
}

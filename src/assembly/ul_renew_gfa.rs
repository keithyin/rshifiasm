//! The renewal graph pass (`gfa_ut.cpp` 2849/3306 + `hic.cpp` 2523-3110):
//! `filter_sg_by_ug` rewrites the read graph from the polished unitig graph,
//! and `gen_bubble_chain` re-identifies bubbles with the poy het
//! reclassification (`identify_bubbles_recal_poy`).

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::bp_graph::adjust_utg_by_primary;
use crate::assembly::bubble_chain::{
    build_bub_graph, dfs_bubble, get_bubbles, BubbleType,
};
use crate::assembly::bubble_identify::asg_bub_pop1_primary_trio;
use crate::assembly::bubble::BubbleBuffer;
use crate::assembly::gfa_clean::get_real_length;
use crate::assembly::hit::{ReadToUnitig, PRIMARY_LABEL};
use crate::command_line::options::AssemblyOptions;
use crate::container::ksort::radix_sort_by_key_32;
use crate::correction::ec_driver::PafList;
use crate::overlap::paf::PafHit;

/// `HOM_PEAK_RATE` (`Purge_Dups.h`).
const HOM_PEAK_RATE: f64 = 0.876;

/// `get_unitig_het_fly` (`hic.cpp` 2523): the coverage ratio of `uid`'s reads
/// against the reads of its unitig neighbourhood; `r_flag` is scratch.
#[allow(clippy::too_many_arguments)]
fn get_unitig_het_fly(
    unitigs: &[MaUtg],
    graph: &Asg,
    read_graph: &Asg,
    uid: u32,
    m_het_cov_thres: i64,
    m_hom_cov_thres: i64,
    sources: &[Vec<PafHit>],
    read_to_unitig: &ReadToUnitig,
    r_flag: &mut [u8],
    m_het_occ: u32,
    m_het_label: u32,
    p_het_label: u32,
    n_het_label: u32,
) -> u32 {
    let mut r_bases = 0i64;
    let mut c_bases = 0i64;
    let mut mark = |unitigs: &[MaUtg], graph: &Asg, uid: u32, value: u8, r_flag: &mut [u8]| {
        for &packed in &unitigs[uid as usize].reads {
            let rid = (packed >> 33) as usize;
            r_flag[rid] = value;
        }
        for side in 0..2u32 {
            for arc in graph.arcs_of((uid << 1) + side) {
                let neighbour = &unitigs[(arc.v >> 1) as usize];
                for &packed in &neighbour.reads {
                    let rid = (packed >> 33) as usize;
                    r_flag[rid] = value;
                }
            }
        }
    };
    mark(unitigs, graph, uid, 1, r_flag);
    for &packed in &unitigs[uid as usize].reads {
        if packed == u64::MAX {
            continue;
        }
        let rid = (packed >> 33) as usize;
        r_bases += i64::from(read_graph.sequences[rid].length);
        for hit in &sources[rid] {
            if hit.exact_region_byte != 1 {
                continue;
            }
            let mut tn = hit.target_identifier as usize;
            if read_graph.sequences[tn].deleted {
                let (mapped, is_unitig) = read_to_unitig.get(tn as u32);
                if mapped == u32::MAX || is_unitig == 1 || read_graph.sequences[mapped as usize].deleted
                {
                    continue;
                }
                tn = mapped as usize;
            }
            if read_graph.sequences[tn].deleted {
                continue;
            }
            match r_flag[tn] {
                1 => c_bases += i64::from(hit.query_end - hit.query_start),
                2 => c_bases += i64::from((hit.query_end - hit.query_start) / 2),
                _ => {}
            }
        }
    }
    mark(unitigs, graph, uid, 0, r_flag);
    let mut cov = 0i64;
    if r_bases > 0 {
        cov = c_bases / r_bases;
    }
    if cov <= m_het_cov_thres && unitigs[uid as usize].reads.len() as u32 >= m_het_occ {
        return m_het_label;
    }
    if cov >= m_hom_cov_thres {
        return n_het_label;
    }
    p_het_label
}

/// `deter_unitig_het_fly` (`hic.cpp` 2718).
#[allow(clippy::too_many_arguments)]
fn deter_unitig_het_fly(
    unitigs: &[MaUtg],
    graph: &Asg,
    read_graph: &Asg,
    uid: u32,
    het_cov_thres: i64,
    sources: &[Vec<PafHit>],
    read_to_unitig: &ReadToUnitig,
    r_flag: &mut [u8],
) -> bool {
    let mut r_bases = 0i64;
    let mut c_bases = 0i64;
    let mut mark = |unitigs: &[MaUtg], graph: &Asg, uid: u32, value: u8, r_flag: &mut [u8]| {
        for &packed in &unitigs[uid as usize].reads {
            let rid = (packed >> 33) as usize;
            r_flag[rid] = value;
        }
        for side in 0..2u32 {
            for arc in graph.arcs_of((uid << 1) + side) {
                let neighbour = &unitigs[(arc.v >> 1) as usize];
                for &packed in &neighbour.reads {
                    let rid = (packed >> 33) as usize;
                    r_flag[rid] = value;
                }
            }
        }
    };
    mark(unitigs, graph, uid, 1, r_flag);
    for &packed in &unitigs[uid as usize].reads {
        if packed == u64::MAX {
            continue;
        }
        let rid = (packed >> 33) as usize;
        r_bases += i64::from(read_graph.sequences[rid].length);
        for hit in &sources[rid] {
            if hit.exact_region_byte != 1 {
                continue;
            }
            let mut tn = hit.target_identifier as usize;
            if read_graph.sequences[tn].deleted {
                let (mapped, is_unitig) = read_to_unitig.get(tn as u32);
                if mapped == u32::MAX || is_unitig == 1 || read_graph.sequences[mapped as usize].deleted
                {
                    continue;
                }
                tn = mapped as usize;
            }
            if read_graph.sequences[tn].deleted {
                continue;
            }
            match r_flag[tn] {
                1 => c_bases += i64::from(hit.query_end - hit.query_start),
                2 => c_bases += i64::from((hit.query_end - hit.query_start) / 2),
                _ => {}
            }
        }
    }
    mark(unitigs, graph, uid, 0, r_flag);
    let mut cov = 0i64;
    if r_bases > 0 {
        cov = c_bases / r_bases;
    }
    cov <= het_cov_thres
}

/// `dfs_subgraph_advance` (`Overlaps.cpp` 26098): the length of one connected
/// component through the `get_unitig` convex walk; `p_bub` reports whether
/// the component carries a bubble.
fn dfs_subgraph_advance(graph: &mut Asg, seen: &mut [u8], x: u32, stack: &mut Vec<u32>) -> (u64, bool) {
    let mut p_bub = false;
    if seen[(x >> 1) as usize] != 0 || graph.sequences[(x >> 1) as usize].deleted {
        return (0, false);
    }
    stack.clear();
    stack.push(x);
    let mut total = 0u64;
    let mut flag_0 = 0u32;
    let mut flag_1 = 0u32;
    while !stack.is_empty() {
        let c_v = stack.pop().unwrap();
        if seen[(c_v >> 1) as usize] != 0 {
            continue;
        }
        let mut convex = 0u32;
        let mut node_len = 0i64;
        let mut base_len = 0i64;
        let mut max_stop_node_len = 0i64;
        let mut max_stop_base_len = 0i64;
        let op = crate::assembly::bubble_identify::get_unitig(
            graph,
            None,
            c_v,
            &mut convex,
            &mut node_len,
            &mut base_len,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            1,
            Some(&mut Vec::new()),
        );
        let mut u_len = base_len;
        if op == crate::assembly::gfa_clean::LOOP {
            return (0, false);
        }
        let e_v = convex ^ 1;
        let mut convex2 = 0u32;
        let mut node_len2 = 0i64;
        let mut base_len2 = 0i64;
        let mut max_stop_node_len2 = 0i64;
        let mut max_stop_base_len2 = 0i64;
        let _ = crate::assembly::bubble_identify::get_unitig(
            graph,
            None,
            e_v,
            &mut convex2,
            &mut node_len2,
            &mut base_len2,
            &mut max_stop_node_len2,
            &mut max_stop_base_len2,
            1,
            Some(&mut Vec::new()),
        );
        u_len = u_len.max(base_len2);
        total += u_len as u64;
        let mut kv_0 = 0u32;
        for arc in graph.arcs_of(c_v ^ 1) {
            if arc.deleted {
                continue;
            }
            kv_0 += 1;
            if seen[(arc.v >> 1) as usize] != 0 {
                continue;
            }
            stack.push(arc.v);
        }
        let mut kv_1 = 0u32;
        for arc in graph.arcs_of(e_v ^ 1) {
            if arc.deleted {
                continue;
            }
            kv_1 += 1;
            if seen[(arc.v >> 1) as usize] != 0 {
                continue;
            }
            stack.push(arc.v);
        }
        if kv_0 > 0 && kv_1 > 0 {
            flag_0 += 1;
        }
        if kv_0 > 1 {
            flag_1 += 1;
        }
        if kv_1 > 1 {
            flag_1 += 1;
        }
        seen[(c_v >> 1) as usize] = 1;
        seen[(e_v >> 1) as usize] = 1;
    }
    if flag_0 > 0 && flag_1 > 1 {
        p_bub = true;
    }
    (total, p_bub)
}

/// `get_bub_pop_max_dist_advance` (`Overlaps.cpp` 26166).
pub fn get_bub_pop_max_dist_advance(graph: &mut Asg) -> u64 {
    let n_vtx = (graph.sequences.len() * 2) as u32;
    let mut seen = vec![0u8; graph.sequences.len()];
    let mut stack: Vec<u32> = Vec::new();
    let mut m_len = 0u64;
    for v in 0..n_vtx {
        if seen[(v >> 1) as usize] != 0 {
            continue;
        }
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        let mut kv = 0u32;
        let mut w = 0u32;
        for arc in graph.arcs_of(v) {
            if arc.deleted {
                continue;
            }
            w = arc.v ^ 1;
            kv += 1;
        }
        if kv == 1 && get_real_length(graph, w) == 1 {
            continue;
        }
        let (c_len, p_bub) = dfs_subgraph_advance(graph, &mut seen, v ^ 1, &mut stack);
        if !p_bub {
            continue;
        }
        if c_len > m_len {
            m_len = c_len;
        }
    }
    for entry in seen.iter_mut() {
        *entry = 0;
    }
    m_len
}

/// The bubble-list half of `identify_bubbles_recal_poy` (`hic.cpp` 2913):
/// rebuilds `bub.list/num/path_len/index` from the unitig graph.
fn identify_bubbles_recal_poy(
    read_graph: &Asg,
    unitig_graph: &mut crate::assembly::unitig::UnitigGraph,
    bub: &mut BubbleType,
    r_het_flag: &mut [u8],
    sources: &[Vec<PafHit>],
    read_to_unitig: &ReadToUnitig,
    options: &AssemblyOptions,
) {
    let n_vtx = (unitig_graph.graph.sequences.len() * 2) as u32;
    let graph = &mut unitig_graph.graph;
    graph.cleanup();
    if !graph.is_symmetric {
        graph.symmetrize();
    }
    let mode: u32 = u32::MAX << 2;
    bub.broken_bubble_count = 0;
    bub.broken_end_bubble_count = 0;
    bub.tangle_bubble_count = 0;
    bub.cross_bubble_count = 0;
    bub.mess_bubble_count = 0;
    let hom_cov: u64 = if options.global_homozygous_coverage_is_set != 0 {
        options.global_homozygous_coverage as u64
    } else {
        (f64::from(options.global_homozygous_coverage) / HOM_PEAK_RATE) as u64
    };
    let het_cov = hom_cov / options.polyploidy as u64;
    let m_het_cov = het_cov + ((het_cov as f64 * 0.333333) as u64);
    let m_hom_cov = het_cov + ((het_cov as f64 * 0.6) as u64);

    let mut buffer = BubbleBuffer {
        info: (0..n_vtx).map(|_| Default::default()).collect(),
        sources: Vec::new(),
        visited: Vec::new(),
        visited_edges: Vec::new(),
    };
    let t_len = get_bub_pop_max_dist_advance(graph);
    {
    // The `graph` alias ends before `reset_inner_bub_het_poy`, which needs
    // the whole `unitig_graph` back.
    bub.list.clear();
    bub.num.clear();
    bub.path_len.clear();
    bub.b_s_idx.clear();
    bub.b_s_idx = vec![u64::MAX; unitig_graph.unitigs.len()];
    bub.chain_weight.clear();
    bub.index = vec![0u32; n_vtx as usize];
    for entry in graph.sequences.iter_mut() {
        entry.coverage_class = 0;
    }
    let unitigs_snapshot: Vec<MaUtg> = unitig_graph.unitigs.clone();
    for v in 0..n_vtx {
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if graph.arc_count(v) < 2 {
            continue;
        }
        if (bub.index[v as usize] & 3) != 0 {
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
    let mut stack: Vec<u32> = Vec::new();
    let mut result: Vec<u32> = Vec::new();
    for v in 0..n_vtx {
        if (bub.index[v as usize] & 3) != 2 {
            continue;
        }
        if asg_bub_pop1_primary_trio(
            graph,
            Some(&unitigs_snapshot),
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
            let mut i = buffer.visited.len() + 1;
            if buffer.visited.len() == 2 || buffer.visited.len() == 3 || buffer.visited.len() == 5 {
                for &node in &buffer.visited {
                    if node == v || node == sink {
                        continue;
                    }
                    dfs_bubble(graph, &mut stack, &mut result, node >> 1, v >> 1, sink >> 1);
                    if result.len() + 3 != buffer.visited.len()
                        && result.len() + 2 != buffer.visited.len()
                    {
                        break;
                    }
                }
            }
            if i == buffer.visited.len() {
                bub.num.push(v);
            } else {
                bub.num.push(v + (1 << 31));
            }
        }
    }
    radix_sort_by_key_32(&mut bub.num, &|&entry: &u32| entry);
    bub.simple_bubble_count = 0;
    for k in 0..bub.num.len() {
        if (bub.num[k] >> 31) == 0 {
            bub.simple_bubble_count += 1;
        }
        let v = (bub.num[k] << 1) >> 1;
        bub.num[k] = bub.list.len() as u32;
        if asg_bub_pop1_primary_trio(
            graph,
            Some(&unitigs_snapshot),
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
            crate::assembly::bubble_chain::asg_bub_backtrack_primary_length(
                graph,
                Some(&unitigs_snapshot),
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
    r_het_flag.fill(0);
    let m_het = bub.num.len();
    let p_het = bub.num.len() + 1;
    for i in 0..unitig_graph.unitigs.len() {
        bub.index[i] = get_unitig_het_fly(
            &unitig_graph.unitigs,
            graph,
            read_graph,
            i as u32,
            m_het_cov as i64,
            m_hom_cov as i64,
            sources,
            read_to_unitig,
            r_het_flag,
            20,
            m_het as u32,
            p_het as u32,
            u32::MAX,
        );
    }
    for i in 0..bub.full_bubble_count {
        let (beg, sink, members, _path) = {
            let (beg, sink, members, path) = get_bubbles(bub, i);
            (beg, sink, members.to_vec(), path)
        };
        let mut n_occ = 0u32;
        for &v in &members {
            bub.index[(v >> 1) as usize] = i as u32;
            n_occ += unitigs_snapshot[(v >> 1) as usize].reads.len() as u32;
        }
        if n_occ > 3 {
            if bub.index[(beg >> 1) as usize] != m_het as u32 {
                bub.index[(beg >> 1) as usize] = u32::MAX;
            }
            if bub.index[(sink >> 1) as usize] != m_het as u32 {
                bub.index[(sink >> 1) as usize] = u32::MAX;
            }
        }
        let v = (beg >> 1) as usize;
        if bub.b_s_idx[v] == u64::MAX {
            bub.b_s_idx[v] <<= 32;
            bub.b_s_idx[v] |= i as u64;
        } else if (bub.b_s_idx[v] & 0xffff_ffff_0000_0000) == 0xffff_ffff_0000_0000 {
            bub.b_s_idx[v] <<= 32;
            bub.b_s_idx[v] |= i as u64;
        }
        let v = (sink >> 1) as usize;
        if bub.b_s_idx[v] == u64::MAX {
            bub.b_s_idx[v] <<= 32;
            bub.b_s_idx[v] |= i as u64;
        } else if (bub.b_s_idx[v] & 0xffff_ffff_0000_0000) == 0xffff_ffff_0000_0000 {
            bub.b_s_idx[v] <<= 32;
            bub.b_s_idx[v] |= i as u64;
        }
    }
    for i in 0..unitig_graph.unitigs.len() {
        if bub.index[i] == m_het as u32 {
            bub.index[i] = p_het as u32;
        }
    }
    let unitigs_snapshot2 = unitig_graph.unitigs.clone();
    build_bub_graph(graph, &unitigs_snapshot2, bub);
    } // end of the `graph` alias scope
    let bubble_unitigs = bub.bubble_unitigs.clone().unwrap_or_default();
    for unitig in &bubble_unitigs {
        if unitig.reads.is_empty() {
            continue;
        }
        for &packed in &unitig.reads {
            reset_inner_bub_het_poy(
                read_graph,
                unitig_graph,
                bub,
                (packed >> 33) as u32,
                t_len,
                &mut buffer,
                m_het_cov,
                m_hom_cov,
                r_het_flag,
                sources,
                read_to_unitig,
            );
        }
    }
    let unitigs_final = unitig_graph.unitigs.clone();
    for i in 0..unitigs_final.len() {
        let flags = if_hom(i as u32, bub);
        if flags {
            continue;
        }
        if deter_unitig_het_fly(
            &unitigs_final,
            &unitig_graph.graph,
            read_graph,
            i as u32,
            (het_cov as f64 * 1.15) as i64,
            sources,
            read_to_unitig,
            r_het_flag,
        ) {
            continue;
        }
        bub.index[i] = u32::MAX;
    }
}

/// `IF_HOM(ID, B)` (`hic.h` 89).
fn if_hom(id: u32, bub: &BubbleType) -> bool {
    bub.index[id as usize] > bub.full_bubble_count as u32 + 1
}

/// `reset_inner_bub_het_poy` (`hic.cpp` 2790).
#[allow(clippy::too_many_arguments)]
fn reset_inner_bub_het_poy(
    read_graph: &Asg,
    unitig_graph: &mut crate::assembly::unitig::UnitigGraph,
    bub: &mut BubbleType,
    bid: u32,
    t_len: u64,
    buffer: &mut BubbleBuffer,
    m_het_cov: u64,
    m_hom_cov: u64,
    r_het_flag: &mut [u8],
    sources: &[Vec<PafHit>],
    read_to_unitig: &ReadToUnitig,
) {
    let (beg, sink, members, _) = get_bubbles(bub, u64::from(bid));
    let members: Vec<u32> = members.to_vec();
    let bn = members.len();
    for m in 0..bn {
        let v = members[m];
        if if_hom(v >> 1, bub) {
            continue;
        }
        if unitig_graph.graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if bub.index[(v >> 1) as usize] == u32::MAX {
            continue;
        }
        if get_unitig_het_fly(
            &unitig_graph.unitigs,
            &unitig_graph.graph,
            read_graph,
            v >> 1,
            m_het_cov as i64,
            m_hom_cov as i64,
            sources,
            read_to_unitig,
            r_het_flag,
            20,
            1,
            0,
            u32::MAX,
        ) == 1
        {
            continue;
        }
        if unitig_graph.graph.arc_count(v) < 2 {
            continue;
        }
        if get_real_length(&unitig_graph.graph, v) < 2 {
            continue;
        }
        if asg_bub_pop1_primary_trio(
            &mut unitig_graph.graph,
            None,
            v,
            t_len,
            buffer,
            u32::MAX,
            u32::MAX,
            false,
            None,
            &[],
            false,
            false,
        ) != 0
        {
            let sink0 = buffer.sources[0];
            let mut socc = 0u32;
            for &z in &buffer.visited {
                if z == v || z == sink0 {
                    continue;
                }
                socc += unitig_graph.unitigs[(z >> 1) as usize].reads.len() as u32;
            }
            if socc > 10 && socc > unitig_graph.unitigs[(v >> 1) as usize].reads.len() as u32 * 3 {
                bub.index[(v >> 1) as usize] = u32::MAX;
            }
        }
        let v2 = v ^ 1;
        if bub.index[(v2 >> 1) as usize] == u32::MAX {
            continue;
        }
        if unitig_graph.graph.arc_count(v2) < 2 {
            continue;
        }
        if get_real_length(&unitig_graph.graph, v2) < 2 {
            continue;
        }
        if asg_bub_pop1_primary_trio(
            &mut unitig_graph.graph,
            None,
            v2,
            t_len,
            buffer,
            u32::MAX,
            u32::MAX,
            false,
            None,
            &[],
            false,
            false,
        ) != 0
        {
            let sink0 = buffer.sources[0];
            let mut socc = 0u32;
            for &z in &buffer.visited {
                if z == v2 || z == sink0 {
                    continue;
                }
                socc += unitig_graph.unitigs[(z >> 1) as usize].reads.len() as u32;
            }
            if socc > 10 && socc > unitig_graph.unitigs[(v2 >> 1) as usize].reads.len() as u32 * 3 {
                bub.index[(v2 >> 1) as usize] = u32::MAX;
            }
        }
    }
    let _ = (beg, sink, PRIMARY_LABEL);
}

/// `gen_bubble_chain` (`gfa_ut.cpp` 3306).
pub fn gen_bubble_chain(
    read_graph: &Asg,
    unitig_graph: &mut crate::assembly::unitig::UnitigGraph,
    sources: &[PafList],
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    options: &mut AssemblyOptions,
    coverage_cut: &mut [crate::assembly::hit::CoverageCut],
    bubble_labels: Option<&mut crate::assembly::bubble::BubbleLabel>,
    avoid_het: bool,
) -> (BubbleType, Vec<u8>) {
    let mut new_rtg_edges: Vec<AsgArc> = Vec::new();
    let mut copy_sg = crate::assembly::unitig::copy_read_graph(read_graph);
    let mut copy_ug = crate::assembly::unitig::copy_untig_graph(unitig_graph);
    let mut trio_flags = vec![0u8; sources.len()];
    let cov = adjust_utg_by_primary(
        &mut copy_ug,
        &mut copy_sg,
        coverage_cut,
        sources,
        reverse_sources,
        read_to_unitig,
        options.maximum_short_tip_read_count as i64 * 2,
        0.15,
        3,
        0.05,
        0.9,
        options.maximum_hang_length,
        options.minimum_overlap_length,
        &mut new_rtg_edges,
        bubble_labels,
        options,
        &mut trio_flags,
        false,
        false,
        true,
    );
    drop(copy_ug);
    drop(copy_sg);
    let mut bub = BubbleType::default();
    // `(*ir_het) = cov->t_ch->ir_het; cov->t_ch->ir_het = NULL; cov->is_r_het = NULL;`
    let cov = match cov {
        Ok(cov) => cov,
        Err(_) => return (BubbleType::default(), vec![0; read_graph.sequences.len()]),
    };
    // `(*ir_het) = cov->t_ch->ir_het; cov->t_ch->ir_het = NULL;
    //  cov->is_r_het = NULL;` — the flags buffer moves out and the graph's
    // own het reclassification writes `r_het_flag` directly.
    let ir_het: Vec<u8> = vec![0; read_graph.sequences.len()];
    let mut ir_het = ir_het;
    if !avoid_het {
        crate::assembly::bubble_chain::identify_bubbles(
            &mut unitig_graph.graph,
            &unitig_graph.unitigs,
            &mut bub,
            &ir_het,
            None,
        );
    } else {
        let mut pr_het = vec![0u8; read_graph.sequences.len()];
        let hits: Vec<Vec<PafHit>> = sources.iter().map(|l| l.hits.clone()).collect();
        identify_bubbles_recal_poy(
            read_graph,
            unitig_graph,
            &mut bub,
            &mut pr_het,
            &hits,
            read_to_unitig,
            options,
        );
        ir_het = pr_het;
    }
    if options.purge_level_primary == 0 {
        for k in 0..unitig_graph.unitigs.len() {
            if if_hom(k as u32, &bub) {
                bub.index[k] = bub.full_bubble_count as u32 + 1;
            }
        }
    }
    (bub, ir_het)
}

/// `filter_sg_by_ug` (`gfa_ut.cpp` 2849): rebuilds the read graph's arcs from
/// the polished unitig graph.
pub fn filter_sg_by_ug(
    read_graph: &mut Asg,
    unitig_graph: &mut crate::assembly::unitig::UnitigGraph,
    sources: &[PafList],
    maximum_hang: i32,
    max_hang_rate: f32,
    minimum_overlap: i32,
) {
    let n_vx = read_graph.sequences.len();
    let _ = n_vx;
    read_graph.arcs.clear();
    read_graph.is_sorted = false;
    for seq in read_graph.sequences.iter_mut() {
        seq.deleted = true;
    }
    for (i, unitig) in unitig_graph.unitigs.iter().enumerate() {
        unitig_graph.graph.sequences[i].coverage_class = PRIMARY_LABEL;
        if unitig_graph.graph.sequences[i].deleted {
            continue;
        }
        for &packed in &unitig.reads {
            read_graph.sequences[(packed >> 33) as usize].deleted = false;
        }
        for m in 0..unitig.reads.len().saturating_sub(1) {
            let v = (unitig.reads[m] >> 32) as u32;
            let w = (unitig.reads[m + 1] >> 32) as u32;
            if let Some(t) = crate::assembly::unitig::gen_spec_edge(
                read_graph,
                sources,
                maximum_hang,
                max_hang_rate,
                minimum_overlap,
                v,
                w,
            ) {
                read_graph.push_arc(t);
            }
            if let Some(t) = crate::assembly::unitig::gen_spec_edge(
                read_graph,
                sources,
                maximum_hang,
                max_hang_rate,
                minimum_overlap,
                w ^ 1,
                v ^ 1,
            ) {
                read_graph.push_arc(t);
            }
        }
        for side in 0..2u32 {
            let v = (i as u32) << 1 | side;
            for arc in unitig_graph.graph.arcs_of(v).to_vec() {
                if arc.deleted {
                    continue;
                }
                let w = arc.v;
                let u = &unitig_graph.unitigs[(v >> 1) as usize];
                let vx = if v & 1 != 0 {
                    ((u.reads[0] >> 32) as u32) ^ 1
                } else {
                    (u.reads[u.reads.len() - 1] >> 32) as u32
                };
                let wu = &unitig_graph.unitigs[(w >> 1) as usize];
                let wx = if w & 1 != 0 {
                    ((wu.reads[wu.reads.len() - 1] >> 32) as u32) ^ 1
                } else {
                    (wu.reads[0] >> 32) as u32
                };
                if let Some(t) = crate::assembly::unitig::gen_spec_edge(
                    read_graph,
                    sources,
                    maximum_hang,
                    max_hang_rate,
                    minimum_overlap,
                    vx,
                    wx,
                ) {
                    read_graph.push_arc(t);
                }
            }
        }
        if unitig.circular {
            for side in [0u32, 1] {
                let v = ((i as u32) << 1) | side;
                let u = &unitig_graph.unitigs[(v >> 1) as usize];
                let vx = if v & 1 != 0 {
                    ((u.reads[0] >> 32) as u32) ^ 1
                } else {
                    (u.reads[u.reads.len() - 1] >> 32) as u32
                };
                let wx = if v & 1 != 0 {
                    (u.reads[u.reads.len() - 1] >> 32) as u32
                } else {
                    ((u.reads[0] >> 32) as u32) ^ 1
                };
                if let Some(t) = crate::assembly::unitig::gen_spec_edge(
                    read_graph,
                    sources,
                    maximum_hang,
                    max_hang_rate,
                    minimum_overlap,
                    vx,
                    wx,
                ) {
                    read_graph.push_arc(t);
                }
            }
        }
    }
    read_graph.index.clear();
    read_graph.is_sorted = false;
    read_graph.cleanup();
    read_graph.symmetrize();
}

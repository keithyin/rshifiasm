//! The tangle-resolution subsystem of the primary unitig cleaning, the
//! counterpart of `resolve_tangles`, `walk_through`, `explore_graph`,
//! `merge_nodes`, `unroll_tangle`, `process_tangles` and their helpers in
//! `Overlaps.cpp`.
//!
//! The reference works on two unitig graphs: `src`, which is modified in
//! place, and a copy that is only used to measure unitig lengths after the
//! merges. [`resolve_tangles`] pulls every merged node of the copy with
//! [`unroll_tangle`], which drops tips, pops bubbles, cuts edges and finally
//! restores the sequences and edges it hid from the rest of the graph.
//!
//! All scratch fields follow the reference macros: `EvaluateLen` is
//! [`MaUtg::start`] and `IsMerge` is [`MaUtg::end`].

use std::collections::VecDeque;

use crate::assembly::asg::{Asg, MaUtg};
use crate::assembly::bubble::{BubbleBuffer, BubbleInfo};
use crate::assembly::bubble_identify::{asg_bub_pop1_primary_trio, get_unitig};
use crate::assembly::gfa_clean::{
    get_real_length, get_real_length_one, LOOP, MUL_OUTPUT,
};
use crate::assembly::hit::{ReadToUnitig, ALTER_LABEL};
use crate::assembly::trio_clean::{asg_seq_drop, check_different_haps, DROP, PLOID};
use crate::assembly::unitig::{copy_untig_graph, UnitigGraph};
use crate::container::ksort::radix_sort_u64;
use crate::correction::ec_driver::PafList;

/// `UNROLL`: merge the tangle nodes into a line.
pub const UNROLL: u32 = 0;
/// `CONVEX`: merge the tangle nodes into a convex shape.
pub const CONVEX: u32 = 1;
/// `CONVEX_M`: [`walk_through`] merged a convex tangle.
pub const CONVEX_M: u32 = 1;
/// `UNROLL_M`: [`walk_through`] unrolled a tangle that continues.
pub const UNROLL_M: u32 = 2;
/// `UNROLL_E`: [`walk_through`] unrolled a tangle that ends the graph.
pub const UNROLL_E: u32 = 3;

/// `TANGLE_MISSED_THRES`: how much of the tangle may be lost.
pub const TANGLE_MISSED_THRES: f64 = 0.6;

/// `Get_vis`: the visit bit of orientation `d` of vertex `v`.
fn get_vis(visit: &[u8], v: u32, d: u32) -> u8 {
    let mask = ((((v << d) & 1) as u8)) + 1;
    visit[(v >> 1) as usize] & mask
}

/// `Set_vis`: raises the visit bit of orientation `d` of vertex `v`.
fn set_vis(visit: &mut [u8], v: u32, d: u32) {
    let mask = ((((v << d) & 1) as u8)) + 1;
    visit[(v >> 1) as usize] |= mask;
}

/// `if_long_tip_length`: whether the chain from `beg_node` is a long tip.
///
/// `unitig_len` is the reference's optional `untigLen` pointer: when it is
/// given the chain is not measured again. `main_len` is the length of the
/// chain the tip branches off.
#[allow(clippy::too_many_arguments)]
fn if_long_tip_length(
    graph: &Asg,
    unitigs: &[MaUtg],
    beg_node: u32,
    unitig_len: Option<u32>,
    min_long_untig: i64,
    max_short_untig: i64,
    short_untig_rate: f32,
    main_len: i64,
) -> u32 {
    let length = match unitig_len {
        Some(length) => length,
        None => {
            let mut end_node = 0u32;
            let mut node_len = 0i64;
            let mut base_len = 0i64;
            let mut max_stop_node_len = 0i64;
            let mut max_stop_base_len = 0i64;
            if get_unitig(
                graph,
                Some(unitigs),
                beg_node,
                &mut end_node,
                &mut node_len,
                &mut base_len,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                1,
                None,
            ) == LOOP
            {
                // A loop is infinitely long.
                return 1;
            }
            node_len as u32
        }
    };

    if length == 0 {
        return 0;
    }
    if (length as f32) < short_untig_rate * (main_len as f32) {
        return 0;
    }
    if i64::from(length) >= max_short_untig {
        return 1;
    }
    if i64::from(length) >= min_long_untig
        && (length as f32) >= short_untig_rate * (main_len as f32)
    {
        return 1;
    }
    0
}

/// `double_check_tangle`: whether every live arc around the tangle stays
/// inside it.
fn double_check_tangle(v_beg: u32, v_end: u32, u_vecs: &[u32], graph: &Asg) -> u32 {
    if v_beg != u32::MAX {
        for arc in graph.arcs_of(v_beg) {
            if arc.deleted {
                continue;
            }
            if !u_vecs.iter().any(|&u| (arc.v >> 1) == (u >> 1)) {
                return 0;
            }
        }
    }

    if v_end != u32::MAX {
        let v = v_end ^ 1;
        for arc in graph.arcs_of(v) {
            if arc.deleted {
                continue;
            }
            if !u_vecs.iter().any(|&u| (arc.v >> 1) == (u >> 1)) {
                return 0;
            }
        }
    }

    for &vec_v in u_vecs {
        for arc in graph.arcs_of(vec_v) {
            if arc.deleted {
                continue;
            }
            let inside = u_vecs.iter().any(|&u| (arc.v >> 1) == (u >> 1));
            if !inside && arc.v != (v_beg ^ 1) && arc.v != v_end {
                return 0;
            }
        }

        let v = vec_v ^ 1;
        for arc in graph.arcs_of(v) {
            if arc.deleted {
                continue;
            }
            let inside = u_vecs.iter().any(|&u| (arc.v >> 1) == (u >> 1));
            if !inside && arc.v != (v_beg ^ 1) && arc.v != v_end {
                return 0;
            }
        }
    }

    1
}

/// `explore_graph`: breadth-first search across the unitigs of a tangle.
///
/// The search stops when it meets a second long unitig (`r_id` is left in
/// `None` when the tangle ends), when it exceeds `single_threshold` or
/// `threshold` reads, or when it can return to `v_beg`. `result` receives
/// every visited vertex.
#[allow(clippy::too_many_arguments)]
fn explore_graph(
    graph: &Asg,
    unitigs: &[MaUtg],
    v_beg: u32,
    single_threshold: f32,
    l_untig_rate_threshold: f32,
    min_long_untig: i64,
    max_short_untig: i64,
    ignore_d: u32,
    buffer: &mut Vec<u32>,
    mut result: Option<&mut Vec<u32>>,
    visit: &mut [u8],
    r_id: &mut u32,
) -> u32 {
    *r_id = u32::MAX;
    let has_result = result.is_some();
    if let Some(result) = result.as_mut() {
        result.clear();
    }
    // `memset(visit, 0, nsg->n_seq)`: the visit bitmap has one byte per read.
    visit.fill(0);
    let mut queue: VecDeque<u32> = VecDeque::new();

    let mut end_node = 0u32;
    let mut node_len = 0i64;
    let mut base_len = 0i64;
    let mut max_stop_node_len = 0i64;
    let mut max_stop_base_len = 0i64;

    buffer.clear();
    if get_unitig(
        graph,
        Some(unitigs),
        v_beg,
        &mut end_node,
        &mut node_len,
        &mut base_len,
        &mut max_stop_node_len,
        &mut max_stop_base_len,
        1,
        Some(&mut *buffer),
    ) == LOOP
    {
        return 0;
    }
    let v_beg_end = end_node;
    let threshold = node_len as u32;
    for i in 0..buffer.len() {
        set_vis(visit, buffer[i], ignore_d);
        set_vis(visit, buffer[i] ^ 1, ignore_d);
    }

    let mut in_count = 0u32;
    let mut end = u32::MAX;
    let mut num_reads = 0u32;
    queue.push_back(v_beg);

    while let Some(mut v) = queue.pop_front() {
        in_count = in_count.wrapping_add(1);

        // `in_count == 1` is the start node itself, which carries no
        // information.
        if in_count != 1 {
            buffer.clear();
            if get_unitig(
                graph,
                Some(unitigs),
                v,
                &mut end_node,
                &mut node_len,
                &mut base_len,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                1,
                Some(&mut *buffer),
            ) == LOOP
            {
                return 0;
            }
            let v_e_len = node_len as u32;

            if if_long_tip_length(
                graph,
                unitigs,
                v,
                Some(v_e_len),
                min_long_untig,
                max_short_untig,
                l_untig_rate_threshold,
                i64::from(threshold),
            ) == 1
            {
                if end == u32::MAX {
                    end = v;
                    continue;
                } else {
                    in_count = u32::MAX;
                    break;
                }
            }

            if (v_e_len as f32) > (threshold as f32) * single_threshold {
                in_count = u32::MAX;
                break;
            }

            num_reads = num_reads.wrapping_add(v_e_len);
            if num_reads > threshold {
                in_count = u32::MAX;
                break;
            }

            if let Some(result) = result.as_mut() {
                for i in 0..buffer.len() {
                    result.push(buffer[i]);
                }
            }
        }

        let tmp = v ^ 1;
        v = end_node;
        let (start, nv) = graph.arc_range(v);
        for k in 0..nv {
            let arc = graph.arcs[start + k];
            if arc.deleted {
                continue;
            }

            if arc.v == v_beg {
                in_count = u32::MAX;
                break;
            }

            if get_vis(visit, arc.v, ignore_d) == 0 {
                queue.push_back(arc.v);

                buffer.clear();
                get_unitig(
                    graph,
                    Some(unitigs),
                    arc.v,
                    &mut end_node,
                    &mut node_len,
                    &mut base_len,
                    &mut max_stop_node_len,
                    &mut max_stop_base_len,
                    1,
                    Some(&mut *buffer),
                );

                for i in 0..buffer.len() {
                    set_vis(visit, buffer[i], ignore_d);
                }
            }
        }

        // For the start node only one direction is explored.
        if in_count != 1 && ignore_d != 0 {
            v = tmp;
            let (start, nv) = graph.arc_range(v);
            for k in 0..nv {
                let arc = graph.arcs[start + k];
                if arc.deleted {
                    continue;
                }

                if arc.v == v_beg {
                    in_count = u32::MAX;
                    break;
                }

                if get_vis(visit, arc.v, ignore_d) == 0 {
                    queue.push_back(arc.v);

                    buffer.clear();
                    get_unitig(
                        graph,
                        Some(unitigs),
                        arc.v,
                        &mut end_node,
                        &mut node_len,
                        &mut base_len,
                        &mut max_stop_node_len,
                        &mut max_stop_base_len,
                        1,
                        Some(&mut *buffer),
                    );

                    for i in 0..buffer.len() {
                        set_vis(visit, buffer[i], ignore_d);
                    }
                }
            }
        }

        // A `break` out of a `for` above must also leave the queue loop; the
        // reference jumps to `termi` with `in_count` already invalidated.
        if in_count == u32::MAX {
            break;
        }
    }

    *r_id = end;
    // `in_count == 1` means the tangle was the long unitig itself.
    if in_count == u32::MAX || in_count == 1 {
        return 0;
    }
    if has_result {
        let result = result.as_deref().expect("the result buffer is present");
        return double_check_tangle(v_beg_end, *r_id, result, graph);
    }
    1
}

/// `append_ma_utg_t`: moves the reads of `v_y` to the end of `v_x`.
///
/// The reference frees `v_y`'s read array and leaves every other field of
/// `v_y` untouched.
#[allow(dead_code)]
fn append_ma_utg_t(v_x: &mut MaUtg, v_y: &mut MaUtg) {
    v_x.reads.extend_from_slice(&v_y.reads);
    v_y.reads.clear();
}

/// `merge_nodes`: merges the tangle nodes into the first one.
///
/// Every node but the first is deleted from the graph and its reads moved to
/// the first. The first node's `EvaluateLen` becomes the largest of the
/// merged nodes, or half the total when that is larger; `IsMerge` is
/// incremented and the length becomes the longest base length. Depending on
/// `merge_type` the start and end vertices are reconnected around the new
/// unitig.
fn merge_nodes(
    unitig_graph: &mut UnitigGraph,
    start_id: u32,
    end_id: u32,
    nodes: &[u32],
    merge_type: u32,
) {
    let mut max_evaluate_len = 0u32;
    let mut max_base_len = 0u32;
    let mut total_evaluate_len = 0u32;
    let first = if nodes.is_empty() {
        usize::MAX
    } else {
        (nodes[0] >> 1) as usize
    };

    if !nodes.is_empty() {
        unitig_graph.graph.delete_sequence(nodes[0] >> 1);
        max_evaluate_len = unitig_graph.unitigs[first].start;
        max_base_len = unitig_graph.unitigs[first].length;
        total_evaluate_len = unitig_graph.unitigs[first].start;
    }

    for i in 1..nodes.len() {
        let y = (nodes[i] >> 1) as usize;
        if max_evaluate_len <= unitig_graph.unitigs[y].start {
            max_evaluate_len = unitig_graph.unitigs[y].start;
            max_base_len = unitig_graph.unitigs[y].length;
        }
        total_evaluate_len = total_evaluate_len.wrapping_add(unitig_graph.unitigs[y].start);

        let y_reads = std::mem::take(&mut unitig_graph.unitigs[y].reads);
        unitig_graph.unitigs[first].reads.extend_from_slice(&y_reads);
        unitig_graph.graph.delete_sequence(nodes[i] >> 1);
    }

    if !nodes.is_empty() {
        unitig_graph.graph.sequences[first].deleted = false;

        unitig_graph.unitigs[first].start = max_evaluate_len;
        total_evaluate_len /= 2;
        if total_evaluate_len > unitig_graph.unitigs[first].start {
            unitig_graph.unitigs[first].start = total_evaluate_len;
        }
        unitig_graph.unitigs[first].end = unitig_graph.unitigs[first].end.wrapping_add(1);

        unitig_graph.unitigs[first].length = max_base_len;
        unitig_graph.graph.sequences[first].length = max_base_len;
        let v = nodes[0];

        if merge_type == UNROLL {
            if start_id != u32::MAX {
                let length = unitig_graph.unitigs[(start_id >> 1) as usize].length;
                unitig_graph.graph.append_edge_to_sorted(
                    start_id, length, v, 0, false, 0, false,
                );
                let length = unitig_graph.unitigs[(v >> 1) as usize].length;
                unitig_graph.graph.append_edge_to_sorted(
                    v ^ 1, length, start_id ^ 1, 0, false, 0, false,
                );
            }

            if end_id != u32::MAX {
                let length = unitig_graph.unitigs[(v >> 1) as usize].length;
                unitig_graph.graph.append_edge_to_sorted(
                    v, length, end_id, 0, false, 0, false,
                );
                let length = unitig_graph.unitigs[(end_id >> 1) as usize].length;
                unitig_graph.graph.append_edge_to_sorted(
                    end_id ^ 1, length, v ^ 1, 0, false, 0, false,
                );
            }
        }

        if merge_type == CONVEX {
            if start_id != u32::MAX {
                let length = unitig_graph.unitigs[(v >> 1) as usize].length;
                unitig_graph.graph.append_edge_to_sorted(
                    v, length, start_id ^ 1, 0, false, 0, false,
                );
                let length = unitig_graph.unitigs[(start_id >> 1) as usize].length;
                unitig_graph.graph.append_edge_to_sorted(
                    start_id, length, v ^ 1, 0, false, 0, false,
                );
            }

            if end_id != u32::MAX {
                let length = unitig_graph.unitigs[(v >> 1) as usize].length;
                unitig_graph.graph.append_edge_to_sorted(
                    v, length, end_id, 0, false, 0, false,
                );
                let length = unitig_graph.unitigs[(end_id >> 1) as usize].length;
                unitig_graph.graph.append_edge_to_sorted(
                    end_id ^ 1, length, v ^ 1, 0, false, 0, false,
                );
            }
        }
    }
}

/// `walk_through`: follows one unambiguous chain and resolves the tangle it
/// runs into.
///
/// Returns [`CONVEX_M`], [`UNROLL_M`], [`UNROLL_E`] or 0; `r_beg`/`r_end`
/// hold the measured chain and `r_next_u_id` the vertex the next walk should
/// start from after an unroll.
#[allow(clippy::too_many_arguments)]
fn walk_through(
    read_graph: &Asg,
    unitig_graph: &mut UnitigGraph,
    reverse_sources: &[PafList],
    min_long_untig: i64,
    max_short_untig: i64,
    l_untig_rate: f32,
    max_node_threshold: f32,
    b_0: &mut Vec<u32>,
    b_1: &mut Vec<u32>,
    u_vecs: &mut Vec<u32>,
    visit: &mut [u8],
    v: u32,
    r_beg: &mut u32,
    r_end: &mut u32,
    r_next_u_id: &mut u32,
    read_to_unitig: &mut ReadToUnitig,
    is_r_het: &[u8],
    trio_flag: u32,
    purge_similarity_threshold: f32,
) -> u32 {
    let _ = trio_flag;
    *r_beg = u32::MAX;
    *r_end = u32::MAX;
    let beg = v;

    if unitig_graph.graph.sequences[(beg >> 1) as usize].deleted
        || unitig_graph.graph.arc_count(beg) == 0
        || get_real_length(&unitig_graph.graph, beg) <= 0
    {
        return 0;
    }
    if get_real_length(&unitig_graph.graph, beg ^ 1) == 1 {
        let end = get_real_length_one(&unitig_graph.graph, beg ^ 1)
            .expect("a single live arc has one target");
        if get_real_length(&unitig_graph.graph, end ^ 1) == 1 {
            return 0;
        }
    }

    let mut end_node = 0u32;
    let mut node_len = 0i64;
    let mut base_len = 0i64;
    let mut max_stop_node_len = 0i64;
    let mut max_stop_base_len = 0i64;
    let primary_flag = get_unitig(
        &unitig_graph.graph,
        Some(&unitig_graph.unitigs),
        beg,
        &mut end_node,
        &mut node_len,
        &mut base_len,
        &mut max_stop_node_len,
        &mut max_stop_base_len,
        1,
        None,
    );
    if primary_flag == LOOP || node_len < min_long_untig {
        return 0;
    }
    let primary_len = node_len;

    if get_real_length(&unitig_graph.graph, end_node) <= 0 {
        return 0;
    }
    *r_beg = beg;
    *r_end = end_node;

    let end = end_node;
    let (end_start, nw) = unitig_graph.graph.arc_range(end);
    let mut n_arc = 0u32;
    for i in 0..nw {
        let arc = unitig_graph.graph.arcs[end_start + i];
        if arc.deleted {
            continue;
        }
        n_arc += 1;
        if if_long_tip_length(
            &unitig_graph.graph,
            &unitig_graph.unitigs,
            arc.v,
            None,
            min_long_untig,
            max_short_untig,
            l_untig_rate,
            primary_len,
        ) == 1
        {
            n_arc = 0;
            break;
        }
    }
    if n_arc == 0 {
        return 0;
    }

    let mut next_u_id = u32::MAX;
    if explore_graph(
        &unitig_graph.graph,
        &unitig_graph.unitigs,
        beg,
        max_node_threshold,
        l_untig_rate,
        min_long_untig,
        max_short_untig,
        1,
        b_0,
        Some(&mut *u_vecs),
        visit,
        &mut next_u_id,
    ) == 1
    {
        *r_next_u_id = next_u_id;
        if next_u_id != u32::MAX
            && u_vecs.len() == 1
            && unitig_graph.unitigs[(u_vecs[0] >> 1) as usize].end > 0
        {
            return 0;
        }

        if next_u_id != u32::MAX {
            let mut next_u_id_verify = u32::MAX;
            let mut return_flag = explore_graph(
                &unitig_graph.graph,
                &unitig_graph.unitigs,
                beg,
                max_node_threshold,
                l_untig_rate,
                min_long_untig,
                max_short_untig,
                0,
                b_0,
                None,
                visit,
                &mut next_u_id_verify,
            );

            if (return_flag == 1 && next_u_id_verify != next_u_id) || return_flag == 0 {
                return_flag = 0;
            } else {
                return_flag = 1;
            }

            if return_flag == 1
                && check_different_haps(
                    &unitig_graph.graph,
                    Some(&unitig_graph.unitigs),
                    read_graph,
                    beg,
                    next_u_id,
                    reverse_sources,
                    b_0,
                    b_1,
                    read_to_unitig,
                    is_r_het,
                    (min_long_untig - 1) as u32,
                    1,
                    purge_similarity_threshold,
                ) == PLOID
            {
                return_flag = 0;
            }

            if return_flag == 0 {
                merge_nodes(unitig_graph, end, next_u_id, u_vecs, CONVEX);
                return CONVEX_M;
            }
        }

        if u_vecs.is_empty() || (next_u_id >> 1) == (beg >> 1) {
            return 0;
        }
        merge_nodes(unitig_graph, end, next_u_id, u_vecs, UNROLL);
        if next_u_id != u32::MAX {
            *r_next_u_id = next_u_id;
            return UNROLL_M;
        } else {
            if !u_vecs.is_empty() {
                *r_next_u_id = u_vecs[0];
            }
            return UNROLL_E;
        }
    }

    0
}

/// `save_all_edges`: hides the arcs of the tangle nodes and records the
/// original coverage classes.
///
/// Every entry of `nodes` is rewritten to `<class>:32 | <unitig id>:32`. The
/// returned pair is the original class of the two boundary unitigs.
fn save_all_edges(
    edges: &mut Vec<u64>,
    graph: &Asg,
    nodes: &mut [u64],
    beg: u32,
    end: u32,
) -> (u32, u32) {
    edges.clear();
    let beg_class = u32::from(graph.sequences[(beg >> 1) as usize].coverage_class);
    let end_class = u32::from(graph.sequences[(end >> 1) as usize].coverage_class);

    for i in 0..nodes.len() {
        let v = nodes[i] as u32;
        let class = u64::from(graph.sequences[v as usize].coverage_class);
        nodes[i] = (class << 32) | u64::from(v);

        if graph.sequences[v as usize].deleted {
            continue;
        }
        if graph.sequences[v as usize].coverage_class == ALTER_LABEL {
            continue;
        }

        let vertex = v << 1;
        for &v in &[vertex, vertex ^ 1] {
            let (start, nv) = graph.arc_range(v);
            for k in 0..nv {
                let arc = graph.arcs[start + k];
                if arc.deleted {
                    continue;
                }
                edges.push((u64::from(arc.overlap) << 32) | (start + k) as u64);
            }
        }
    }

    let v = beg;
    if !graph.sequences[(v >> 1) as usize].deleted
        && graph.sequences[(v >> 1) as usize].coverage_class != ALTER_LABEL
    {
        let (start, nv) = graph.arc_range(v);
        for k in 0..nv {
            let arc = graph.arcs[start + k];
            if arc.deleted {
                continue;
            }
            let mut i = 0usize;
            while i < nodes.len() {
                if (arc.v >> 1) as u64 == (nodes[i] & 0xffff_ffff) {
                    break;
                }
                i += 1;
            }
            if i == nodes.len() {
                continue;
            }
            edges.push((u64::from(arc.overlap) << 32) | (start + k) as u64);
        }
    }

    let v = end ^ 1;
    if !graph.sequences[(v >> 1) as usize].deleted
        && graph.sequences[(v >> 1) as usize].coverage_class != ALTER_LABEL
    {
        let (start, nv) = graph.arc_range(v);
        for k in 0..nv {
            let arc = graph.arcs[start + k];
            if arc.deleted {
                continue;
            }
            let mut i = 0usize;
            while i < nodes.len() {
                if (arc.v >> 1) as u64 == (nodes[i] & 0xffff_ffff) {
                    break;
                }
                i += 1;
            }
            if i == nodes.len() {
                continue;
            }
            edges.push((u64::from(arc.overlap) << 32) | (start + k) as u64);
        }
    }

    (beg_class, end_class)
}

/// `recover_edges`: restores the tangle sequences, their classes and, when
/// `is_recover_edges` is set, every hidden arc.
#[allow(clippy::too_many_arguments)]
fn recover_edges(
    graph: &mut Asg,
    edges: &[u64],
    nodes: &[u64],
    beg: u32,
    end: u32,
    beg_class: u32,
    end_class: u32,
    is_recover_edges: bool,
) {
    for &node in nodes {
        let v = node as u32;
        graph.sequences[v as usize].deleted = false;
        graph.sequences[v as usize].coverage_class = (node >> 32) as u8;
    }
    graph.sequences[(beg >> 1) as usize].coverage_class = beg_class as u8;
    graph.sequences[(beg >> 1) as usize].deleted = false;
    graph.sequences[(end >> 1) as usize].coverage_class = end_class as u8;
    graph.sequences[(end >> 1) as usize].deleted = false;

    if is_recover_edges {
        for &edge in edges {
            graph.arcs[edge as u32 as usize].deleted = false;
        }
    }
}

/// `drop_tips_at_tangle`: drops every tip that leaves the tangle until none is
/// left.
fn drop_tips_at_tangle(
    graph: &mut Asg,
    nodes: &[u64],
    beg: u32,
    end: u32,
    buffer: &mut Vec<u32>,
) -> u32 {
    let mut next_tips = 1u32;
    let mut reduce = 0u32;

    while next_tips > 0 {
        next_tips = 0;

        for i in 0..nodes.len() {
            let v = nodes[i] as u32;
            let sequence = graph.sequences[v as usize];
            if sequence.deleted || sequence.coverage_class == ALTER_LABEL {
                continue;
            }
            let vertex = v << 1;

            for k in 0..2u32 {
                let v = vertex | k;
                if get_real_length(graph, v ^ 1) != 0 {
                    continue;
                }

                buffer.clear();
                let mut convex = 0u32;
                let mut ll = 0i64;
                let mut tmp = 0i64;
                let mut max_stop_node_len = 0i64;
                let mut max_stop_base_len = 0i64;
                let return_flag = get_unitig(
                    graph,
                    None,
                    v,
                    &mut convex,
                    &mut ll,
                    &mut tmp,
                    &mut max_stop_node_len,
                    &mut max_stop_base_len,
                    1,
                    Some(&mut *buffer),
                );
                if return_flag == LOOP {
                    continue;
                }
                if return_flag == MUL_OUTPUT {
                    next_tips += 1;
                }

                let mut m = 0usize;
                while m < buffer.len() {
                    if (buffer[m] >> 1) == (beg >> 1) || (buffer[m] >> 1) == (end >> 1) {
                        break;
                    }
                    m += 1;
                }
                // The chain contains one of the two boundary unitigs.
                if m != buffer.len() {
                    continue;
                }

                for m in 0..buffer.len() {
                    graph.sequences[(buffer[m] >> 1) as usize].coverage_class = ALTER_LABEL;
                }
                for m in 0..buffer.len() {
                    asg_seq_drop(graph, buffer[m] >> 1);
                }

                reduce += 1;
            }
        }
    }

    reduce
}

/// `pop_bubble_at_tangle`: pops the bubbles starting at the two boundaries and
/// at every tangle node.
#[allow(clippy::too_many_arguments)]
fn pop_bubble_at_tangle(
    unitig_graph: &mut UnitigGraph,
    max_dist: u64,
    nodes: &[u64],
    beg: u32,
    end: u32,
    positive_flag: u32,
    negative_flag: u32,
    trio_labels: &[u8],
) -> u64 {
    if !unitig_graph.graph.is_symmetric {
        unitig_graph.graph.symmetrize();
    }
    let vertex_count = unitig_graph.graph.sequences.len() * 2;
    let mut buffer = BubbleBuffer {
        info: vec![BubbleInfo::default(); vertex_count],
        ..BubbleBuffer::default()
    };
    let mut n_pop = 0u64;

    let v = beg;
    {
        let sequence = unitig_graph.graph.sequences[(v >> 1) as usize];
        if !sequence.deleted
            && sequence.coverage_class != ALTER_LABEL
            && get_real_length(&unitig_graph.graph, v) >= 2
        {
            n_pop += asg_bub_pop1_primary_trio(
                &mut unitig_graph.graph,
                Some(&unitig_graph.unitigs),
                v,
                max_dist,
                &mut buffer,
                positive_flag,
                negative_flag,
                true,
                None,
                trio_labels,
                false,
                true,
            );
        }
    }

    let v = end ^ 1;
    {
        let sequence = unitig_graph.graph.sequences[(v >> 1) as usize];
        if !sequence.deleted
            && sequence.coverage_class != ALTER_LABEL
            && get_real_length(&unitig_graph.graph, v) >= 2
        {
            n_pop += asg_bub_pop1_primary_trio(
                &mut unitig_graph.graph,
                Some(&unitig_graph.unitigs),
                v,
                max_dist,
                &mut buffer,
                positive_flag,
                negative_flag,
                true,
                None,
                trio_labels,
                false,
                true,
            );
        }
    }

    for i in 0..nodes.len() {
        let v = nodes[i] as u32;
        let sequence = unitig_graph.graph.sequences[v as usize];
        if sequence.deleted || sequence.coverage_class == ALTER_LABEL {
            continue;
        }
        let vertex = v << 1;

        for k in 0..2u32 {
            let v = vertex | k;
            if get_real_length(&unitig_graph.graph, v) <= 1 {
                continue;
            }
            n_pop += asg_bub_pop1_primary_trio(
                &mut unitig_graph.graph,
                Some(&unitig_graph.unitigs),
                v,
                max_dist,
                &mut buffer,
                positive_flag,
                negative_flag,
                true,
                None,
                trio_labels,
                false,
                true,
            );
        }
    }

    n_pop
}

/// `find_spec_node`: whether `v_dest` is reachable from `v_beg` within
/// `max_nodes` of chain length.
///
/// `visit` is the shared visit bitmap, cleared on entry; the vertices of
/// `forbidden_nodes` are never entered.
#[allow(clippy::too_many_arguments)]
fn find_spec_node(
    graph: &Asg,
    unitigs: Option<&[MaUtg]>,
    v_beg: u32,
    v_dest: u32,
    forbidden_nodes: &[u32],
    max_nodes: i64,
    buffer: &mut Vec<u32>,
    visit: &mut [u8],
) -> u32 {
    // `memset(visit, 0, nsg->n_seq)`.
    visit.fill(0);
    let mut queue: VecDeque<u32> = VecDeque::new();
    let mut is_found = 0u32;
    let mut total_nodes = 0i64;

    let mut end_node = 0u32;
    let mut node_len = 0i64;
    let mut base_len = 0i64;
    let mut max_stop_node_len = 0i64;
    let mut max_stop_base_len = 0i64;

    buffer.clear();
    if get_unitig(
        graph,
        unitigs,
        v_beg,
        &mut end_node,
        &mut node_len,
        &mut base_len,
        &mut max_stop_node_len,
        &mut max_stop_base_len,
        1,
        Some(&mut *buffer),
    ) == LOOP
    {
        let mut i = 0usize;
        while i < buffer.len() {
            if buffer[i] == v_dest {
                break;
            }
            i += 1;
        }
        return if i != buffer.len() { 1 } else { 0 };
    }

    for i in 0..buffer.len() {
        set_vis(visit, buffer[i], 0);
    }
    for &node in forbidden_nodes {
        set_vis(visit, node, 0);
    }

    queue.push_back(v_beg);

    while let Some(v) = queue.pop_front() {
        buffer.clear();
        let return_flag = get_unitig(
            graph,
            unitigs,
            v,
            &mut end_node,
            &mut node_len,
            &mut base_len,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            1,
            Some(&mut *buffer),
        );

        let mut i = 0usize;
        while i < buffer.len() {
            if buffer[i] == v_dest {
                break;
            }
            i += 1;
        }
        if i != buffer.len() {
            is_found = 1;
            break;
        }

        if return_flag == LOOP {
            break;
        }

        total_nodes += node_len;
        if total_nodes > max_nodes {
            break;
        }

        let v = end_node;
        let (start, nv) = graph.arc_range(v);
        for k in 0..nv {
            let arc = graph.arcs[start + k];
            if arc.deleted {
                continue;
            }

            if get_vis(visit, arc.v, 0) != 0 {
                continue;
            }

            let mut un_visit = 0u32;
            buffer.clear();
            get_unitig(
                graph,
                unitigs,
                arc.v,
                &mut end_node,
                &mut node_len,
                &mut base_len,
                &mut max_stop_node_len,
                &mut max_stop_base_len,
                1,
                Some(&mut *buffer),
            );
            let mut i = 0usize;
            while i < buffer.len() {
                if get_vis(visit, buffer[i], 0) != 0 {
                    break;
                }
                i += 1;
            }
            if i == buffer.len() {
                un_visit = 1;
            }

            if un_visit != 0 {
                queue.push_back(arc.v);
                for i in 0..buffer.len() {
                    set_vis(visit, buffer[i], 0);
                }
            }
        }
    }

    is_found
}

/// `drop_useless_edges`: deletes every arc of the tangle that cannot reach the
/// other side.
fn drop_useless_edges(
    unitig_graph: &mut UnitigGraph,
    buffer: &mut Vec<u32>,
    visit: &mut [u8],
    beg: u32,
    end: u32,
    nodes: &[u64],
    nodes_threshold: u32,
) -> u32 {
    // `drop_useless_edges` receives the threshold as a `uint32_t` and hands it
    // to `find_spec_node` as a `long long`.
    let nodes_threshold = i64::from(nodes_threshold);
    let mut n_reduce = 0u32;

    let mut v = beg;
    let mut dest = end;
    if get_real_length(&unitig_graph.graph, v) > 1 {
        let w = v ^ 1;
        let (start, nv) = unitig_graph.graph.arc_range(v);
        for k in 0..nv {
            let arc = unitig_graph.graph.arcs[start + k];
            if arc.deleted {
                continue;
            }
            if find_spec_node(
                &unitig_graph.graph,
                Some(&unitig_graph.unitigs),
                arc.v,
                dest,
                &[w],
                nodes_threshold,
                buffer,
                visit,
            ) == 0
            {
                unitig_graph.graph.arcs[start + k].deleted = true;
                unitig_graph
                    .graph
                    .delete_arc(arc.v ^ 1, ((arc.ul >> 32) as u32) ^ 1, true);
                n_reduce += 1;
            }
        }
    }

    v = end ^ 1;
    dest = beg ^ 1;
    if get_real_length(&unitig_graph.graph, v) > 1 {
        let w = v ^ 1;
        let (start, nv) = unitig_graph.graph.arc_range(v);
        for k in 0..nv {
            let arc = unitig_graph.graph.arcs[start + k];
            if arc.deleted {
                continue;
            }
            if find_spec_node(
                &unitig_graph.graph,
                Some(&unitig_graph.unitigs),
                arc.v,
                dest,
                &[w],
                nodes_threshold,
                buffer,
                visit,
            ) == 0
            {
                unitig_graph.graph.arcs[start + k].deleted = true;
                unitig_graph
                    .graph
                    .delete_arc(arc.v ^ 1, ((arc.ul >> 32) as u32) ^ 1, true);
                n_reduce += 1;
            }
        }
    }

    let mut reach_beg = [0u32; 2];
    let mut reach_end = [0u32; 2];
    for i in 0..nodes.len() {
        let v = nodes[i] as u32;
        let sequence = unitig_graph.graph.sequences[v as usize];
        if sequence.deleted || sequence.coverage_class == ALTER_LABEL {
            continue;
        }
        let vertex = v << 1;

        for k in 0..2u32 {
            let v = vertex | k;

            reach_beg[k as usize] = find_spec_node(
                &unitig_graph.graph,
                Some(&unitig_graph.unitigs),
                v,
                beg ^ 1,
                &[end],
                nodes_threshold,
                buffer,
                visit,
            );
            reach_end[k as usize] = find_spec_node(
                &unitig_graph.graph,
                Some(&unitig_graph.unitigs),
                v,
                end,
                &[beg ^ 1],
                nodes_threshold,
                buffer,
                visit,
            );
        }

        // The reference sums the reachability of the two directions.
        let reach_sum_0 = reach_beg[0].wrapping_add(reach_end[0]);
        let reach_sum_1 = reach_beg[1].wrapping_add(reach_end[1]);
        if reach_sum_0 == 0 || reach_sum_1 == 0 {
            continue;
        }

        let mut k = u32::MAX;
        let mut drop_dest = u32::MAX;
        let mut forbid_dest = u32::MAX;
        if reach_sum_0 == 1 && reach_sum_1 > 1 {
            k = 0;
        }
        if reach_sum_0 > 1 && reach_sum_1 == 1 {
            k = 1;
        }
        if k == u32::MAX {
            continue;
        }

        let drop_node = ((nodes[i] as u32) << 1) | (1 - k);
        if reach_beg[k as usize] == 1 {
            drop_dest = end;
            forbid_dest = beg ^ 1;
        }
        if reach_end[k as usize] == 1 {
            drop_dest = beg ^ 1;
            forbid_dest = end;
        }
        if drop_dest == u32::MAX || forbid_dest == u32::MAX {
            continue;
        }

        let (start, nv) = unitig_graph.graph.arc_range(drop_node);
        for k in 0..nv {
            let arc = unitig_graph.graph.arcs[start + k];
            if arc.deleted {
                continue;
            }
            if find_spec_node(
                &unitig_graph.graph,
                Some(&unitig_graph.unitigs),
                arc.v,
                drop_dest,
                &[forbid_dest],
                nodes_threshold,
                buffer,
                visit,
            ) == 0
            {
                unitig_graph.graph.arcs[start + k].deleted = true;
                unitig_graph
                    .graph
                    .delete_arc(arc.v ^ 1, ((arc.ul >> 32) as u32) ^ 1, true);
                n_reduce += 1;
            }
        }
    }

    n_reduce
}

/// `process_tangles`: drop tips, pop bubbles and cut useless edges until the
/// two boundaries are connected.
#[allow(clippy::too_many_arguments)]
fn process_tangles(
    unitig_graph: &mut UnitigGraph,
    nodes: &[u64],
    beg: u32,
    end: u32,
    buffer: &mut Vec<u32>,
    visit: &mut [u8],
    trio_flag: u32,
    total_node_len: i64,
    total_base_len: i64,
    trio_labels: &[u8],
) -> u32 {
    let mut n_reduce = 1u32;

    while n_reduce != 0 {
        while n_reduce != 0 {
            n_reduce = 0;
            n_reduce = n_reduce.wrapping_add(drop_tips_at_tangle(
                &mut unitig_graph.graph,
                nodes,
                beg,
                end,
                buffer,
            ));
            n_reduce = n_reduce.wrapping_add(
                pop_bubble_at_tangle(
                    unitig_graph,
                    total_base_len as u64,
                    nodes,
                    beg,
                    end,
                    trio_flag,
                    DROP,
                    trio_labels,
                ) as u32,
            );
        }
        n_reduce = drop_useless_edges(
            unitig_graph,
            buffer,
            visit,
            beg,
            end,
            nodes,
            total_node_len as u32,
        );
    }

    let mut is_found = 0u32;
    let mut v = beg;
    loop {
        if v == end {
            is_found = 1;
        }
        if is_found == 1 {
            break;
        }

        let kv = get_real_length(&unitig_graph.graph, v);
        if kv != 1 {
            break;
        }
        let w = get_real_length_one(&unitig_graph.graph, v)
            .expect("a single live arc has one target");
        if get_real_length(&unitig_graph.graph, w ^ 1) != 1 {
            break;
        }

        v = w;
        if v == beg {
            break;
        }
    }

    is_found
}

/// `cut_edges_progressive`: cuts the tangle's weakest edges until the
/// boundaries are connected, restoring each cut that disconnects them.
#[allow(clippy::too_many_arguments)]
fn cut_edges_progressive(
    unitig_graph: &mut UnitigGraph,
    edges: &mut [u64],
    drop_ratio: f32,
    nodes: &[u64],
    beg: u32,
    end: u32,
    buffer: &mut Vec<u32>,
    visit: &mut [u8],
    trio_flag: u32,
    total_node_len: i64,
    total_base_len: i64,
    max_arc_n: u32,
    trio_labels: &[u8],
) -> u32 {
    radix_sort_u64(edges);
    let mut is_found = 0u32;
    let limit = edges.len().min(max_arc_n as usize);

    for k in 0..limit {
        let arc_index = edges[k] as u32 as usize;
        let arc = unitig_graph.graph.arcs[arc_index];
        if arc.deleted {
            continue;
        }
        let v = (arc.ul >> 32) as u32;
        let w = arc.v;
        let sequence = unitig_graph.graph.sequences[(v >> 1) as usize];
        if sequence.deleted || sequence.coverage_class == ALTER_LABEL {
            continue;
        }
        let sequence = unitig_graph.graph.sequences[(w >> 1) as usize];
        if sequence.deleted || sequence.coverage_class == ALTER_LABEL {
            continue;
        }

        let (start, nv) = unitig_graph.graph.arc_range(v);
        if nv <= 1 {
            continue;
        }

        let mut ov_max = 0u32;
        let mut kv = 0u32;
        for i in 0..nv {
            let a = unitig_graph.graph.arcs[start + i];
            if a.deleted {
                continue;
            }
            if ov_max < a.overlap {
                ov_max = a.overlap;
            }
            kv += 1;
        }
        if kv <= 1 {
            continue;
        }
        if (arc.overlap as f32) > ov_max as f32 * drop_ratio {
            continue;
        }

        unitig_graph.graph.delete_arc(v, w, true);
        unitig_graph.graph.delete_arc(w ^ 1, v ^ 1, true);

        let ban = beg ^ 1;
        if find_spec_node(
            &unitig_graph.graph,
            Some(&unitig_graph.unitigs),
            beg,
            end,
            &[ban],
            total_node_len,
            buffer,
            visit,
        ) == 0
        {
            is_found = 0;
            unitig_graph.graph.delete_arc(v, w, false);
            unitig_graph.graph.delete_arc(w ^ 1, v ^ 1, false);
            break;
        }

        is_found = process_tangles(
            unitig_graph,
            nodes,
            beg,
            end,
            buffer,
            visit,
            trio_flag,
            total_node_len,
            total_base_len,
            trio_labels,
        );
        if is_found == 1 {
            break;
        }
    }

    is_found
}

/// `unroll_tangle`: resolves one tangle of `target` (the reference's `src`).
///
/// `length_graph` is the copy the reference calls `length_ug`: it is only read
/// when measuring the two boundary unitigs. `nodes` holds the unitig ids of
/// the merged nodes and is rewritten by [`save_all_edges`].
#[allow(clippy::too_many_arguments)]
fn unroll_tangle(
    target: &mut UnitigGraph,
    length_graph: &UnitigGraph,
    nodes: &mut Vec<u64>,
    beg: u32,
    end: u32,
    edges: &mut Vec<u64>,
    trio_flag: u32,
    buffer: &mut Vec<u32>,
    visit: &mut [u8],
    drop_ratio: f32,
    trio_labels: &[u8],
) {
    edges.clear();
    let (beg_class, end_class) = save_all_edges(edges, &target.graph, nodes, beg, end);

    let mut total_node_len = 0i64;
    let mut total_base_len = 0i64;
    for i in 0..nodes.len() {
        let v = nodes[i] as u32;
        total_node_len += target.unitigs[v as usize].reads.len() as i64;
        total_base_len += target.graph.sequences[v as usize].length as i64;
    }
    total_node_len += target.unitigs[(beg >> 1) as usize].reads.len() as i64
        + target.unitigs[(end >> 1) as usize].reads.len() as i64
        + 1;
    total_base_len += target.graph.sequences[(beg >> 1) as usize].length as i64
        + target.graph.sequences[(end >> 1) as usize].length as i64
        + 1;
    // Each node might be visited twice.
    total_node_len *= 2;
    total_base_len *= 2;

    let mut is_found = process_tangles(
        target,
        nodes,
        beg,
        end,
        buffer,
        visit,
        trio_flag,
        total_node_len,
        total_base_len,
        trio_labels,
    );

    if is_found == 0 {
        is_found = cut_edges_progressive(
            target,
            edges,
            drop_ratio,
            nodes,
            beg,
            end,
            buffer,
            visit,
            trio_flag,
            total_node_len,
            total_base_len,
            48,
            trio_labels,
        );
    }

    if is_found == 1 {
        let mut convex = 0u32;
        let mut beg_len = 0i64;
        let mut end_len = 0i64;
        let mut tmp = 0i64;
        let mut max_stop_node_len = 0i64;
        let mut max_stop_base_len = 0i64;
        get_unitig(
            &length_graph.graph,
            Some(&length_graph.unitigs),
            beg ^ 1,
            &mut convex,
            &mut beg_len,
            &mut tmp,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            1,
            None,
        );
        get_unitig(
            &length_graph.graph,
            Some(&length_graph.unitigs),
            end,
            &mut convex,
            &mut end_len,
            &mut tmp,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            1,
            None,
        );

        let mut miss_len = 0i64;
        for i in 0..nodes.len() {
            let v = nodes[i] as u32;
            if target.graph.sequences[v as usize].deleted
                || target.graph.sequences[v as usize].coverage_class == ALTER_LABEL
            {
                miss_len += target.unitigs[v as usize].reads.len() as i64;
            }
        }
        // Resolved when little enough of the tangle was lost.
        if (miss_len as f64) <= (beg_len + end_len) as f64 * TANGLE_MISSED_THRES {
            let w = beg ^ 1;
            is_found = find_spec_node(
                &target.graph,
                Some(&target.unitigs),
                beg,
                end,
                &[w],
                total_node_len,
                buffer,
                visit,
            );
        } else {
            is_found = 0;
        }
    }

    recover_edges(
        &mut target.graph,
        edges,
        nodes,
        beg,
        end,
        beg_class,
        end_class,
        is_found == 0,
    );
}

/// `resolve_tangles`: unrolls every merged tangle of the primary unitig graph.
///
/// The graph is modified in place; the read map is first relabelled to the
/// current unitigs and every unitig entry is removed again at the end.
#[allow(clippy::too_many_arguments)]
pub fn resolve_tangles(
    unitig_graph: &mut UnitigGraph,
    read_graph: &Asg,
    reverse_sources: &[PafList],
    min_long_untig: i64,
    max_short_untig: i64,
    l_untig_rate: f32,
    max_node_threshold: f32,
    read_to_unitig: &mut ReadToUnitig,
    is_r_het: &[u8],
    trio_flag: u32,
    drop_ratio: f32,
    trio_labels: &[u8],
    purge_similarity_threshold: f32,
) {
    // `EvaluateLen` has to start from fresh read counts.
    for v in 0..unitig_graph.graph.sequences.len() {
        if unitig_graph.graph.sequences[v].deleted {
            continue;
        }
        if unitig_graph.unitigs[v].reads.is_empty() {
            continue;
        }
        unitig_graph.unitigs[v].start = unitig_graph.unitigs[v].reads.len() as u32;
    }

    let mut copy = copy_untig_graph(unitig_graph);
    for v in 0..copy.graph.sequences.len() {
        copy.unitigs[v].start = copy.unitigs[v].reads.len() as u32;
        copy.unitigs[v].end = 0;
    }

    let mut visit = vec![0u8; copy.graph.sequences.len()];
    let mut b_0: Vec<u32> = Vec::new();
    let mut b_1: Vec<u32> = Vec::new();
    let mut u_vecs: Vec<u32> = Vec::new();
    let mut next_u_id = u32::MAX;
    let vertex_count = copy.graph.sequences.len() * 2;

    loop {
        let mut n_reduce = 0u32;
        for v in 0..vertex_count as u32 {
            let sv = v;
            loop {
                let mut beg = u32::MAX;
                let mut end = u32::MAX;
                let flag = walk_through(
                    read_graph,
                    &mut copy,
                    reverse_sources,
                    min_long_untig,
                    max_short_untig,
                    l_untig_rate,
                    max_node_threshold,
                    &mut b_0,
                    &mut b_1,
                    &mut u_vecs,
                    &mut visit,
                    sv,
                    &mut beg,
                    &mut end,
                    &mut next_u_id,
                    read_to_unitig,
                    is_r_het,
                    trio_flag,
                    purge_similarity_threshold,
                );
                n_reduce = n_reduce.wrapping_add(flag);
                if flag != UNROLL_M {
                    break;
                }
            }
        }
        if n_reduce == 0 {
            break;
        }
    }

    copy.graph.cleanup();
    copy.graph.symmetrize();

    // Relabel the reads of every live source unitig to that unitig.
    for v in 0..unitig_graph.graph.sequences.len() {
        let unitig_id = v as u32;
        if unitig_graph.unitigs[v].reads.is_empty() {
            continue;
        }
        if unitig_graph.graph.sequences[v].deleted {
            continue;
        }
        for i in 0..unitig_graph.unitigs[v].reads.len() {
            let read = (unitig_graph.unitigs[v].reads[i] >> 33) as u32;
            read_to_unitig.set(
                read,
                unitig_id,
                true,
                Some(read_graph.sequences[read as usize].coverage_class),
            );
        }
    }

    // Rewrite the copy's read lists with the unitig ids they now belong to.
    for v in 0..copy.graph.sequences.len() {
        if copy.graph.sequences[v].deleted {
            continue;
        }
        let unitig = &mut copy.unitigs[v];
        let mut m = 0usize;
        let mut pre = u32::MAX;
        for i in 0..unitig.reads.len() {
            let read = (unitig.reads[i] >> 33) as u32;
            let (unitig_id, is_unitig) = read_to_unitig.get(read);
            if is_unitig != 1 || unitig_id == u32::MAX {
                continue;
            }
            if i > 0 && pre == unitig_id {
                continue;
            }
            pre = unitig_id;
            unitig.reads[m] = u64::from(unitig_id);
            m += 1;
        }
        unitig.reads.truncate(m);
    }

    let vertex_count = copy.graph.sequences.len();
    for i in 0..vertex_count {
        let v = i as u32;
        if copy.graph.sequences[i].deleted {
            continue;
        }
        if copy.unitigs[i].end == 0 {
            continue;
        }
        if copy.graph.sequences[i].coverage_class == ALTER_LABEL {
            continue;
        }

        let v = v << 1;
        if get_real_length(&copy.graph, v) != 1 {
            continue;
        }
        let w = get_real_length_one(&copy.graph, v).expect("a single live arc has one target");
        if get_real_length(&copy.graph, w ^ 1) != 1 {
            continue;
        }
        if copy.unitigs[(w >> 1) as usize].end != 0 {
            continue;
        }
        let beg = w ^ 1;

        let v = v ^ 1;
        if get_real_length(&copy.graph, v) != 1 {
            continue;
        }
        let w = get_real_length_one(&copy.graph, v).expect("a single live arc has one target");
        if get_real_length(&copy.graph, w ^ 1) != 1 {
            continue;
        }
        if copy.unitigs[(w >> 1) as usize].end != 0 {
            continue;
        }
        let end = w;

        let mut nodes = copy.unitigs[i].reads.clone();
        let mut edges: Vec<u64> = Vec::new();
        unroll_tangle(
            unitig_graph,
            &copy,
            &mut nodes,
            beg,
            end,
            &mut edges,
            trio_flag,
            &mut b_0,
            &mut visit,
            drop_ratio,
            trio_labels,
        );
    }

    read_to_unitig.clear_unitig_entries();

    // Restore `EvaluateLen` for the rest of the cleaning.
    for v in 0..unitig_graph.graph.sequences.len() {
        if unitig_graph.graph.sequences[v].deleted {
            continue;
        }
        if unitig_graph.unitigs[v].reads.is_empty() {
            continue;
        }
        unitig_graph.unitigs[v].start = unitig_graph.unitigs[v].reads.len() as u32;
    }
}

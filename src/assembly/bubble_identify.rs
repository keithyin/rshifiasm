//! Bubble identification over the read graph, the counterpart of
//! `asg_bub_pop1_label`, the `get_unitig_s`/`dfs_subgraph_s_advance` distance
//! search and `asg_arc_identify_simple_bubbles_multi` in `Overlaps.cpp` and
//! `gfa_ut.cpp`.
//!
//! The reference threads the worker over the graph; the port runs the same
//! workers sequentially, which is exactly what `kt_for` does for the
//! single-thread default configuration that produces the golden outputs.

use crate::assembly::asg::{Asg, MaUtg};
use crate::assembly::bubble::{BubbleBuffer, BubbleInfo, BubbleLabel};
use crate::assembly::hap_cov::{asg_bub_backtrack_primary_cov, get_num_trio_flag, HapCov};
use crate::assembly::hamming::FlySimp;
use crate::assembly::hit::{ALTER_LABEL, HAP_LABEL};
use crate::assembly::gfa_clean::{
    get_real_length, get_real_length_one, END_TIPS, LONG_TIPS_UNDER_MAX_EXT, LOOP, MUL_INPUT,
    MUL_OUTPUT, TWO_INPUT, TWO_OUTPUT,
};
use crate::assembly::unitig::{UnitigGraph, FATHER, MOTHER};
use crate::logging;

/// `get_unitig`: walks the unambiguous chain from `beg_node`, measuring it.
///
/// With `unitigs == None` this is the read-graph variant `get_unitig_s`; with
/// a unitig graph the node length accumulates the unitig's `EvaluateLen`
/// (`ma_utg_t::start`) instead of one per vertex, as the reference does.
#[allow(clippy::too_many_arguments)]
pub fn get_unitig(
    graph: &Asg,
    unitigs: Option<&[MaUtg]>,
    beg_node: u32,
    end_node: &mut u32,
    node_len: &mut i64,
    base_len: &mut i64,
    max_stop_node_len: &mut i64,
    max_stop_base_len: &mut i64,
    stops_threshold: u32,
    mut visited: Option<&mut Vec<u32>>,
) -> u8 {
    let mut v = beg_node;
    let mut n_stops = 0u32;
    let mut pre_base_len = 0i64;
    let mut pre_node_len = 0i64;
    *max_stop_node_len = 0;
    *max_stop_base_len = 0;
    *node_len = 0;
    *base_len = 0;
    *end_node = u32::MAX;

    let return_flag;
    loop {
        let kv = get_real_length(graph, v);
        *end_node = v;
        if let Some(unitigs) = unitigs {
            *node_len += i64::from(unitigs[(v >> 1) as usize].start);
        } else {
            *node_len += 1;
        }
        if let Some(visited) = visited.as_deref_mut() {
            visited.push(v);
        }
        if kv != 1 {
            *base_len += i64::from(graph.sequences[(v >> 1) as usize].length);
        }
        if kv == 0 {
            return_flag = END_TIPS;
            break;
        }
        if kv > 1 {
            return_flag = MUL_OUTPUT;
            break;
        }

        let w = get_real_length_one(graph, v).expect("kv == 1 has a live target");
        if get_real_length(graph, w ^ 1) != 1 {
            n_stops += 1;
            if n_stops >= stops_threshold {
                *base_len += i64::from(graph.sequences[(v >> 1) as usize].length);
                return_flag = MUL_INPUT;
                break;
            }
        }

        if let Some(arc) = graph.arcs_of(v).iter().find(|arc| !arc.deleted) {
            *base_len += i64::from(arc.node_length());
        }

        if get_real_length(graph, w ^ 1) != 1 {
            let current = *base_len - pre_base_len;
            pre_base_len = *base_len;
            if current > *max_stop_base_len {
                *max_stop_base_len = current;
            }

            let current = *node_len - pre_node_len;
            pre_node_len = *node_len;
            if current > *max_stop_node_len {
                *max_stop_node_len = current;
            }
        }

        v = w;
        if v == beg_node {
            return_flag = LOOP;
            break;
        }
    }

    let current = *base_len - pre_base_len;
    let _ = pre_base_len;
    if current > *max_stop_base_len {
        *max_stop_base_len = current;
    }

    let current = *node_len - pre_node_len;
    let _ = pre_node_len;
    if current > *max_stop_node_len {
        *max_stop_node_len = current;
    }

    return_flag
}

/// `get_unitig_s`: the read-graph variant, `get_unitig(g, NULL, ...)`.
#[allow(clippy::too_many_arguments)]
pub fn get_unitig_s(
    graph: &Asg,
    beg_node: u32,
    end_node: &mut u32,
    node_len: &mut i64,
    base_len: &mut i64,
    max_stop_node_len: &mut i64,
    max_stop_base_len: &mut i64,
    stops_threshold: u32,
    visited: Option<&mut Vec<u32>>,
) -> u8 {
    get_unitig(
        graph,
        None,
        beg_node,
        end_node,
        node_len,
        base_len,
        max_stop_node_len,
        max_stop_base_len,
        stops_threshold,
        visited,
    )
}

/// `dfs_subgraph_s_advance`: sums the maximum unitig length reachable from
/// `x`, and flags whether the region looks like a bubble (`p_bubble`).
pub fn dfs_subgraph_s_advance(
    graph: &Asg,
    buffer: &mut BubbleBuffer,
    x: u32,
    p_bubble: &mut u32,
) -> u64 {
    let mut length = 0u64;
    *p_bubble = 0;
    if buffer.info[(x >> 1) as usize].state || graph.sequences[(x >> 1) as usize].deleted {
        return 0;
    }

    buffer.sources.clear();
    buffer.sources.push(x);
    buffer.visited.clear();

    let mut flag_0 = 0u32;
    let mut flag_1 = 0u32;

    while !buffer.sources.is_empty() {
        let current_vertex = buffer.sources.pop().expect("the loop keeps it non-empty");
        if buffer.info[(current_vertex >> 1) as usize].state {
            continue;
        }

        buffer.visited.clear();
        let mut convex = 0u32;
        let mut node_len = 0i64;
        let mut base_len = 0i64;
        let mut max_stop_node_len = 0i64;
        let mut max_stop_base_len = 0i64;
        let operation = get_unitig_s(
            graph,
            current_vertex,
            &mut convex,
            &mut node_len,
            &mut base_len,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            1,
            Some(&mut buffer.visited),
        );

        let mut unitig_length = base_len;
        {
            let (visited, info) = (&buffer.visited, &mut buffer.info);
            for &node in visited {
                info[(node >> 1) as usize].state = true;
            }
        }
        if operation == LOOP {
            return 0;
        }

        let end_vertex = convex ^ 1;
        buffer.visited.clear();
        let operation = get_unitig_s(
            graph,
            end_vertex,
            &mut convex,
            &mut node_len,
            &mut base_len,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            1,
            Some(&mut buffer.visited),
        );
        let _ = operation;
        unitig_length = unitig_length.max(base_len);

        length += unitig_length as u64;

        let mut kv_0 = 0u32;
        buffer.visited.clear();
        for arc in graph.arcs_of(current_vertex ^ 1) {
            if arc.deleted {
                continue;
            }
            kv_0 += 1;
            if buffer.info[(arc.v >> 1) as usize].state {
                continue;
            }
            buffer.sources.push(arc.v);
        }

        let mut kv_1 = 0u32;
        for arc in graph.arcs_of(end_vertex ^ 1) {
            if arc.deleted {
                continue;
            }
            kv_1 += 1;
            if buffer.info[(arc.v >> 1) as usize].state {
                continue;
            }
            buffer.sources.push(arc.v);
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
    }

    if flag_0 > 0 && flag_1 > 1 {
        *p_bubble = 1;
    }
    length
}

/// `get_s_bub_pop_max_dist_advance`: the largest bubble subgraph, which sets
/// the distance limit used by [`asg_bub_pop1_label`].
pub fn get_s_bub_pop_max_dist_advance(graph: &Asg, buffer: &mut BubbleBuffer) -> u64 {
    let vertex_count = graph.sequences.len() * 2;
    let mut max_len = 0u64;

    for v in 0..vertex_count as u32 {
        if buffer.info[(v >> 1) as usize].state {
            continue;
        }
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }

        let mut w = 0u32;
        let mut kv = 0u32;
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

        let mut p_bubble = 0u32;
        let current_len = dfs_subgraph_s_advance(graph, buffer, v ^ 1, &mut p_bubble);
        if p_bubble == 0 {
            continue;
        }
        if current_len > max_len {
            max_len = current_len;
        }
    }

    for k in 0..graph.sequences.len() {
        buffer.info[k].state = false;
    }

    buffer.sources.clear();
    buffer.visited.clear();
    max_len
}

/// `asg_bub_pop1_label`: labels every vertex of the bubble starting at `v0`.
///
/// Returns 1 when a bubble was found; the end vertex is left in
/// `buffer.sources[0]`, as the reference leaves it in `b->S.a[0]`.
pub fn asg_bub_pop1_label(graph: &Asg, v0: u32, max_dist: u64, buffer: &mut BubbleBuffer) -> u64 {
    if graph.sequences[(v0 >> 1) as usize].deleted {
        return 0;
    }
    if get_real_length(graph, v0) < 2 {
        return 0;
    }

    buffer.sources.clear();
    buffer.visited.clear();
    buffer.visited_edges.clear();
    buffer.info[v0 as usize].distance = 0;
    buffer.sources.push(v0);

    let mut n_pending = 0u32;
    let mut n_tips = 0u32;
    let mut tip_end = u32::MAX;
    let mut success = false;

    'search: loop {
        let v = buffer.sources.pop().expect("the loop condition keeps it non-empty");
        let distance = buffer.info[v as usize].distance;
        let (start, count) = graph.arc_range(v);
        let mut early = false;

        let mut i = 0usize;
        while i < count {
            let arc = graph.arcs[start + i];
            i += 1;
            if arc.deleted {
                continue;
            }

            let w = arc.v;
            let length = arc.node_length();
            if (w >> 1) == (v0 >> 1) {
                early = true;
                break;
            }
            if u64::from(distance) + u64::from(length) > max_dist {
                early = true;
                break;
            }

            let info = &mut buffer.info[w as usize];
            if !info.state {
                buffer.visited.push(w);
                info.parent = v;
                info.state = true;
                info.distance = distance.wrapping_add(length);
                info.remaining = get_real_length(graph, w ^ 1);
                n_pending += 1;
            } else if distance.wrapping_add(length) < info.distance {
                info.distance = distance.wrapping_add(length);
                info.parent = v;
            }

            info.remaining = info.remaining.wrapping_sub(1);
            if info.remaining == 0 {
                let outgoing = get_real_length(graph, w);
                if outgoing > 0 {
                    buffer.sources.push(w);
                } else {
                    if n_tips != 0 {
                        early = true;
                        break;
                    }
                    n_tips += 1;
                    tip_end = w;
                }
                n_pending = n_pending.wrapping_sub(1);
            }
        }

        if early {
            break 'search;
        }

        if n_tips == 1 {
            if tip_end != u32::MAX && n_pending == 0 && buffer.sources.is_empty() {
                buffer.sources.push(tip_end);
                success = true;
                break 'search;
            }
            break 'search;
        }

        if buffer.sources.is_empty() {
            break 'search;
        }
        if !(buffer.sources.len() > 1 || n_pending > 0) {
            success = true;
            break 'search;
        }
    }

    for &node in &buffer.visited {
        let info = &mut buffer.info[node as usize];
        info.state = false;
        info.distance = 0;
    }

    if success {
        1
    } else {
        0
    }
}

/// `asg_bub_pop1_primary_trio`: the bubble pop used by the primary cleaning
/// passes.
///
/// The end vertex is left in `buffer.sources[0]`. When `unitigs` is given the
/// pop also accumulates the trio weights and read counts of every visited
/// unitig; when `cov` is given and the pop succeeds, the bubble coverage is
/// redistributed and the alternate side's reads flagged through
/// [`asg_bub_backtrack_primary_cov`].
///
/// `is_pop` is the reference's `is_pop` argument: when it is set and the pop
/// succeeds, the bubble is recorded in the graph through
/// [`asg_bub_backtrack_primary`]. `hic_clean` uses that to remove the popped
/// bubbles; the topology clean passes never do.
#[allow(clippy::too_many_arguments)]
pub fn asg_bub_pop1_primary_trio(
    graph: &mut Asg,
    unitigs: Option<&[MaUtg]>,
    v0: u32,
    max_distance: u64,
    buffer: &mut BubbleBuffer,
    positive_flag: u32,
    negative_flag: u32,
    is_pop: bool,
    mut cov: Option<&mut HapCov>,
    trio_labels: &[u8],
    is_update_chain: bool,
    keep_d: bool,
) -> u64 {
    let sequence = graph.sequences[(v0 >> 1) as usize];
    if sequence.deleted || sequence.coverage_class == ALTER_LABEL {
        return 0;
    }
    if get_real_length(graph, v0) < 2 {
        return 0;
    }

    buffer.sources.clear();
    buffer.visited.clear();
    buffer.visited_edges.clear();
    {
        let info = &mut buffer.info[v0 as usize];
        info.positive_count = 0;
        info.distance = 0;
        info.negative_count = 0;
        info.node_count = 0;
        info.non_positive_count = 0;
    }
    buffer.sources.push(v0);

    let mut n_pending = 0u32;
    let mut is_first = true;
    let mut n_tips = 0u32;
    let mut tip_end = u32::MAX;
    let mut non_positive_flag = u32::MAX;
    if positive_flag == u32::from(FATHER) {
        non_positive_flag = u32::from(MOTHER);
    }
    if positive_flag == u32::from(MOTHER) {
        non_positive_flag = u32::from(FATHER);
    }

    let mut success = false;
    'search: loop {
        let v = buffer
            .sources
            .pop()
            .expect("the loop condition keeps the stack non-empty");
        let parent = buffer.info[v as usize];
        let distance = parent.distance;
        let positive = parent.positive_count;
        let negative = parent.negative_count;
        let non_positive = parent.non_positive_count;
        let (start, count) = graph.arc_range(v);
        let mut too_far = false;

        let mut i = 0usize;
        while i < count {
            let arc_index = start + i;
            let arc = graph.arcs[arc_index];
            i += 1;
            if arc.deleted {
                continue;
            }

            let w = arc.v;
            let mut length = arc.node_length();
            if (w >> 1) == (v0 >> 1) {
                break 'search;
            }
            if is_first && keep_d {
                length = 0;
            }

            buffer.visited_edges.push(arc_index as u32);
            if u64::from(distance.wrapping_add(length)) > max_distance {
                too_far = true;
                break;
            }

            let (cur_positive, cur_negative, cur_non_positive, cur_node_count) =
                if let Some(unitigs) = unitigs {
                    let cur_positive =
                        get_num_trio_flag(unitigs, graph, w >> 1, positive_flag, trio_labels);
                    let cur_negative =
                        get_num_trio_flag(unitigs, graph, w >> 1, negative_flag, trio_labels);
                    let mut cur_non_positive = 0;
                    if non_positive_flag != u32::MAX {
                        cur_non_positive = get_num_trio_flag(
                            unitigs,
                            graph,
                            w >> 1,
                            non_positive_flag,
                            trio_labels,
                        );
                    }
                    let cur_node_count = unitigs[(w >> 1) as usize].reads.len() as u32;
                    (cur_positive, cur_negative, cur_non_positive, cur_node_count)
                } else {
                    (0, 0, 0, 1)
                };

            if !buffer.info[w as usize].state {
                buffer.visited.push(w);
                let info = &mut buffer.info[w as usize];
                info.parent = v;
                info.state = true;
                info.distance = distance.wrapping_add(length);
                info.positive_count = positive.wrapping_add(cur_positive);
                info.negative_count = negative.wrapping_add(cur_negative);
                info.non_positive_count = non_positive.wrapping_add(cur_non_positive);
                info.node_count = parent.node_count.wrapping_add(cur_node_count);
                info.remaining = get_real_length(graph, w ^ 1);
                n_pending += 1;
            } else {
                let info = &mut buffer.info[w as usize];
                let current_weight = i64::from(positive.wrapping_add(cur_positive))
                    - (i64::from(negative.wrapping_add(cur_negative))
                        + i64::from(non_positive.wrapping_add(cur_non_positive)));
                let maximum_weight = i64::from(info.positive_count)
                    - (i64::from(info.negative_count) + i64::from(info.non_positive_count));
                let candidate_count = parent.node_count.wrapping_add(cur_node_count);
                let candidate_distance = distance.wrapping_add(length);

                let replace = if current_weight > maximum_weight {
                    true
                } else if current_weight == maximum_weight {
                    if candidate_count > info.node_count {
                        true
                    } else {
                        candidate_count == info.node_count && candidate_distance > info.distance
                    }
                } else {
                    false
                };

                if replace {
                    info.parent = v;
                    info.negative_count = negative.wrapping_add(cur_negative);
                    info.positive_count = positive.wrapping_add(cur_positive);
                    info.node_count = candidate_count;
                    info.non_positive_count = non_positive.wrapping_add(cur_non_positive);
                }
                if candidate_distance < info.distance {
                    info.distance = candidate_distance;
                }
            }

            let remaining = {
                let info = &mut buffer.info[w as usize];
                info.remaining = info.remaining.wrapping_sub(1);
                info.remaining
            };
            if remaining == 0 {
                let outgoing = get_real_length(graph, w);
                if outgoing > 0 {
                    buffer.sources.push(w);
                } else if n_tips != 0 {
                    // `goto pop_reset`: a second tip aborts the whole search.
                    break 'search;
                } else {
                    n_tips += 1;
                    tip_end = w;
                }
                n_pending = n_pending.wrapping_sub(1);
            }
        }

        is_first = false;
        if n_tips == 1 {
            if tip_end != u32::MAX && n_pending == 0 && buffer.sources.is_empty() {
                buffer.sources.push(tip_end);
                success = true;
                break 'search;
            }
            break 'search;
        }

        // `i < nv` means the arc walk stopped at the distance limit.
        if too_far || buffer.sources.is_empty() {
            break 'search;
        }
        if !(buffer.sources.len() > 1 || n_pending > 0) {
            success = true;
            break 'search;
        }
    }

    if success {
        if let (Some(cov_ref), Some(unitigs)) = (cov.as_deref_mut(), unitigs) {
            asg_bub_backtrack_primary_cov(graph, unitigs, v0, buffer, cov_ref, is_update_chain);
        }
        if is_pop {
            asg_bub_backtrack_primary(graph, v0, buffer);
        }
    }

    for &node in &buffer.visited {
        let info = &mut buffer.info[node as usize];
        info.state = false;
        info.positive_count = 0;
        info.distance = 0;
        info.negative_count = 0;
        info.node_count = 0;
        info.non_positive_count = 0;
    }

    if success {
        1
    } else {
        0
    }
}

/// `asg_pop_bubble_primary_trio`: the two-pass bubble pop of the primary
/// cleaning.
///
/// The first pass searches every vertex's bubble with `utg`/`cov` unset and
/// records which unitigs belong to a bubble; the second pass pops the bubbles
/// again at their start vertices, this time with the trio weights and the
/// coverage backtracking enabled.
#[allow(clippy::too_many_arguments)]
pub fn asg_pop_bubble_primary_trio(
    unitig_graph: &mut UnitigGraph,
    max_distance_override: Option<u64>,
    positive_flag: u32,
    negative_flag: u32,
    mut cov: Option<&mut HapCov>,
    trio_labels: &[u8],
    is_update_chain: bool,
    mut fly: Option<&mut FlySimp>,
) -> u64 {
    if !unitig_graph.graph.is_symmetric {
        unitig_graph.graph.symmetrize();
    }

    let vertex_count = unitig_graph.graph.sequences.len() * 2;
    let mut buffer = BubbleBuffer {
        info: vec![BubbleInfo::default(); vertex_count],
        ..BubbleBuffer::default()
    };
    let max_distance = match max_distance_override {
        Some(distance) => distance,
        None => get_s_bub_pop_max_dist_advance(&unitig_graph.graph, &mut buffer),
    };
    let mut bubble_flag = vec![0u8; vertex_count];
    let mut n_pop = 0u64;

    if max_distance > 0 {
        // First pass: locate every bubble and flag the unitigs inside it.
        for v in 0..vertex_count as u32 {
            if bubble_flag[v as usize] != 0 {
                continue;
            }
            {
                let graph = &unitig_graph.graph;
                if graph.arc_count(v) < 2
                    || graph.sequences[(v >> 1) as usize].deleted
                    || graph.sequences[(v >> 1) as usize].coverage_class == ALTER_LABEL
                {
                    continue;
                }
                let mut live_arcs = 0u32;
                for arc in graph.arcs_of(v) {
                    if !arc.deleted {
                        live_arcs += 1;
                    }
                }
                if live_arcs < 2 {
                    continue;
                }
            }

            if asg_bub_pop1_primary_trio(
                &mut unitig_graph.graph,
                None,
                v,
                max_distance,
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
                let end = buffer.sources[0];
                for i in 0..buffer.visited.len() {
                    let node = buffer.visited[i];
                    if node == v || node == end {
                        continue;
                    }
                    bubble_flag[node as usize] = 1;
                    bubble_flag[(node ^ 1) as usize] = 1;
                }
                bubble_flag[v as usize] = 2;
                bubble_flag[(end ^ 1) as usize] = 3;
            }
        }

        // Second pass: pop the bubbles whose start vertex was flagged. When a
        // fly-simplification is attached the pop is refined through the cached
        // phase graph, as `refine_bubble_popping` does in the reference.
        for v in 0..vertex_count as u32 {
            if bubble_flag[v as usize] != 2 {
                continue;
            }
            {
                let graph = &unitig_graph.graph;
                if graph.arc_count(v) < 2
                    || graph.sequences[(v >> 1) as usize].deleted
                    || graph.sequences[(v >> 1) as usize].coverage_class == ALTER_LABEL
                {
                    continue;
                }
                let mut live_arcs = 0u32;
                for arc in graph.arcs_of(v) {
                    if !arc.deleted {
                        live_arcs += 1;
                    }
                }
                if live_arcs <= 1 {
                    continue;
                }
            }

            match fly.as_deref_mut() {
                Some(fly) => {
                    n_pop += crate::assembly::hamming::refine_bubble_popping(
                        unitig_graph,
                        &mut buffer,
                        v,
                        max_distance,
                        positive_flag,
                        negative_flag,
                        cov.as_deref_mut(),
                        is_update_chain,
                        trio_labels,
                        fly,
                    );
                }
                None => {
                    n_pop += asg_bub_pop1_primary_trio(
                        &mut unitig_graph.graph,
                        Some(unitig_graph.unitigs.as_slice()),
                        v,
                        max_distance,
                        &mut buffer,
                        positive_flag,
                        negative_flag,
                        true,
                        cov.as_deref_mut(),
                        trio_labels,
                        is_update_chain,
                        false,
                    );
                }
            }
        }
    }

    if n_pop != 0 {
        unitig_graph.graph.cleanup();
    }
    n_pop
}

/// `asg_bub_backtrack_primary`: records a successful primary pop in the graph
/// and removes the bubble's interior edges.
///
/// `buffer.sources[0]` must hold the sink found by
/// [`asg_bub_pop1_primary_trio`], `buffer.visited` its nodes and
/// `buffer.visited_edges` the edge slots the search walked. Every visited node
/// becomes `ALTER_LABLE`, the path from the sink back to `v0` becomes
/// `HAP_LABLE` (the sink itself keeps its old class, and `v0` is left alone),
/// and the visited edges are deleted in both directions unless both ends are
/// already alternate.
pub fn asg_bub_backtrack_primary(graph: &mut Asg, v0: u32, buffer: &BubbleBuffer) {
    // `b->S.a[0]` is the sink of this bubble.
    let sink = buffer.sources[0];
    let saved_class = graph.sequences[(sink >> 1) as usize].coverage_class;

    // First remove all nodes in this bubble.
    for &node in &buffer.visited {
        graph.sequences[(node >> 1) as usize].coverage_class = ALTER_LABEL;
    }

    // Recover the nodes on the sink-to-beg path. The sink is temporarily
    // labelled HAP as well; the reference restores it just below because a
    // unitig may carry many reads.
    let mut v = sink;
    loop {
        let u = buffer.info[v as usize].parent;
        graph.sequences[(v >> 1) as usize].coverage_class = HAP_LABEL;
        v = u;
        if v == v0 {
            break;
        }
    }
    graph.sequences[(sink >> 1) as usize].coverage_class = saved_class;

    // Remove all edges (self and reverse for each edge) in this bubble.
    for &edge in &buffer.visited_edges {
        let arc = graph.arcs[edge as usize];
        let query = (arc.ul >> 33) as u32;
        let target = arc.v >> 1;
        if graph.sequences[query as usize].coverage_class == ALTER_LABEL
            && graph.sequences[target as usize].coverage_class == ALTER_LABEL
        {
            continue;
        }
        graph.arcs[edge as usize].deleted = true;
        graph.delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, true);
    }

    // The path nodes survive the pop.
    let mut v = sink;
    loop {
        let u = buffer.info[v as usize].parent;
        graph.sequences[(v >> 1) as usize].deleted = false;
        graph.delete_arc(u, v, false);
        graph.delete_arc(v ^ 1, u ^ 1, false);
        v = u;
        if v == v0 {
            break;
        }
    }
}

/// `asg_is_single_edge`: the raw out-degree of `v^1`, with the reference's
/// diagnostic when `start_node` is not among its targets.
///
/// The reference deliberately does not skip deleted arcs here.
fn asg_is_single_edge(graph: &Asg, v: u32, start_node: u32) -> i32 {
    let arcs = graph.arcs_of(v ^ 1);
    let mut flag = false;
    for arc in arcs {
        if arc.v >> 1 == start_node {
            flag = true;
        }
    }
    if !flag {
        eprintln!("****ERROR");
    }
    arcs.len() as i32
}

/// `detect_single_path`: follows the raw single-arc chain from `beg_node`.
pub fn detect_single_path(
    graph: &Asg,
    beg_node: u32,
    end_node: &mut u32,
    length: &mut i64,
    mut visited: Option<&mut Vec<u32>>,
) -> u8 {
    let mut v = beg_node;
    *length = 0;

    loop {
        *length += 1;
        let nv = graph.arc_count(v);
        *end_node = v;

        if let Some(visited) = visited.as_deref_mut() {
            visited.push(v >> 1);
        }

        if nv == 0 {
            return END_TIPS;
        }
        if nv == 2 {
            return TWO_OUTPUT;
        }
        if nv > 2 {
            return MUL_OUTPUT;
        }

        // Up to here `nv == 1`; `rnv` is at least 1.
        let target = graph.arcs_of(v)[0].v;
        let rnv = asg_is_single_edge(graph, target, v >> 1);
        v = target;
        *end_node = v;
        if rnv == 2 {
            *length += 1;
            if let Some(visited) = visited.as_deref_mut() {
                visited.push(v >> 1);
            }
            return TWO_INPUT;
        }
        if rnv > 2 {
            *length += 1;
            if let Some(visited) = visited.as_deref_mut() {
                visited.push(v >> 1);
            }
            return MUL_INPUT;
        }

        if (v >> 1) == (beg_node >> 1) {
            return LOOP;
        }
    }
}

/// `detect_simple_bubble`: whether `beg_node` splits into two paths that meet
/// again at the same vertex.
pub fn detect_simple_bubble(
    graph: &Asg,
    beg_node: u32,
    end_node: &mut u32,
    min_len: &mut i64,
    mut visited: Option<&mut Vec<u32>>,
) -> bool {
    if graph.arc_count(beg_node) != 2 {
        return false;
    }

    let first = graph.arcs_of(beg_node)[0].v;
    let second = graph.arcs_of(beg_node)[1].v;
    if asg_is_single_edge(graph, first, beg_node >> 1) != 1
        || asg_is_single_edge(graph, second, beg_node >> 1) != 1
    {
        return false;
    }

    if let Some(visited) = visited.as_deref_mut() {
        visited.push(beg_node >> 1);
    }

    let mut e1 = 0u32;
    let mut l1 = 0i64;
    let mut e2 = 0u32;
    let mut l2 = 0i64;
    if detect_single_path(graph, first, &mut e1, &mut l1, visited.as_deref_mut()) == TWO_INPUT
        && detect_single_path(graph, second, &mut e2, &mut l2, visited.as_deref_mut()) == TWO_INPUT
        && e1 == e2
    {
        *end_node = e1;
        *min_len = l1.min(l2) + 1;
        return true;
    }

    false
}

/// `detect_single_path_with_single_bubbles`: a single-path walk that steps
/// through simple bubbles on the way.
pub fn detect_single_path_with_single_bubbles(
    graph: &Asg,
    beg_node: u32,
    end_node: &mut u32,
    length: &mut i64,
    mut visited: Option<&mut Vec<u32>>,
    max_ext: u32,
) -> u8 {
    let mut v = beg_node;
    *length = 0;
    let max_ext = i64::from(max_ext);

    loop {
        let nv = graph.arc_count(v);
        *end_node = v;
        *length += 1;

        if *length > max_ext {
            return LONG_TIPS_UNDER_MAX_EXT;
        }

        if let Some(visited) = visited.as_deref_mut() {
            visited.push(v >> 1);
        }

        if nv == 0 {
            return END_TIPS;
        }
        if nv == 2 {
            let previous_length = visited.as_ref().map(|visited| visited.len());
            let mut bubble_end = 0u32;
            let mut bubble_len = 0i64;
            if !detect_simple_bubble(
                graph,
                v,
                &mut bubble_end,
                &mut bubble_len,
                visited.as_deref_mut(),
            ) {
                if let (Some(visited), Some(previous_length)) =
                    (visited.as_deref_mut(), previous_length)
                {
                    visited.truncate(previous_length);
                }
                return TWO_OUTPUT;
            }
            *length = *length + bubble_len - 2;
            v = bubble_end;
            continue;
        }
        if nv > 2 {
            return MUL_OUTPUT;
        }

        // Up to here `nv == 1`.
        let target = graph.arcs_of(v)[0].v;
        let rnv = asg_is_single_edge(graph, target, v >> 1);
        v = target;
        *end_node = v;
        if rnv == 2 {
            if let Some(visited) = visited.as_deref_mut() {
                visited.push(v >> 1);
            }
            *length += 1;
            return TWO_INPUT;
        }
        if rnv > 2 {
            if let Some(visited) = visited.as_deref_mut() {
                visited.push(v >> 1);
            }
            *length += 1;
            return MUL_INPUT;
        }

        if (v >> 1) == (beg_node >> 1) {
            return LOOP;
        }
    }
}

/// `detect_bubble_end_with_bubbles`: whether two walks starting at `beg_node1`
/// and `beg_node2` reach the same merge vertex.
pub fn detect_bubble_end_with_bubbles(
    graph: &Asg,
    beg_node1: u32,
    beg_node2: u32,
    end_node: &mut u32,
    min_len: &mut i64,
    mut visited: Option<&mut Vec<u32>>,
) -> bool {
    let mut e1 = 0u32;
    let mut l1 = 0i64;
    let mut e2 = 0u32;
    let mut l2 = 0i64;

    if detect_single_path_with_single_bubbles(
        graph,
        beg_node1,
        &mut e1,
        &mut l1,
        visited.as_deref_mut(),
        u32::MAX,
    ) == TWO_INPUT
        && detect_single_path_with_single_bubbles(
            graph,
            beg_node2,
            &mut e2,
            &mut l2,
            visited.as_deref_mut(),
            u32::MAX,
        ) == TWO_INPUT
        && e1 == e2
    {
        *end_node = e1;
        *min_len = l1.min(l2);
        return true;
    }

    false
}

/// `check_if_cross`: whether `v` sits at the crossing of two bubbles.
///
/// `minimum_length_threshold` is the reference's global `min_thres`, set to
/// `maximum_short_tip + 1` by the graph builder.
pub fn check_if_cross(graph: &Asg, v: u32, minimum_length_threshold: i64) -> bool {
    let mut nodes = [0u32; 5];
    if graph.sequences[(v >> 1) as usize].deleted {
        return false;
    }

    if graph.arc_count(v) != 2 {
        return false;
    }
    let arcs = graph.arcs_of(v);
    if asg_is_single_edge(graph, arcs[0].v, v >> 1) != 2
        || asg_is_single_edge(graph, arcs[1].v, v >> 1) != 2
    {
        return false;
    }
    if arcs[0].v == arcs[1].v {
        return false;
    }

    nodes[0] = v;
    nodes[1] = arcs[0].v ^ 1;
    nodes[2] = arcs[1].v ^ 1;
    if graph.arc_count(nodes[0]) != 2
        || graph.arc_count(nodes[1]) != 2
        || graph.arc_count(nodes[2]) != 2
    {
        return false;
    }

    let first_arcs = graph.arcs_of(nodes[1]);
    let second_arcs = graph.arcs_of(nodes[2]);
    if first_arcs[0].v == first_arcs[1].v {
        return false;
    }
    if second_arcs[0].v == second_arcs[1].v {
        return false;
    }

    if first_arcs[0].v == (nodes[0] ^ 1) {
        nodes[3] = first_arcs[1].v ^ 1;
    } else if first_arcs[1].v == (nodes[0] ^ 1) {
        nodes[3] = first_arcs[0].v ^ 1;
    }

    if second_arcs[0].v == (nodes[0] ^ 1) {
        nodes[4] = second_arcs[1].v ^ 1;
    } else if second_arcs[1].v == (nodes[0] ^ 1) {
        nodes[4] = second_arcs[0].v ^ 1;
    }

    if nodes[3] != nodes[4] {
        return false;
    }

    if graph.arc_count(nodes[0]) != 2
        || graph.arc_count(nodes[1]) != 2
        || graph.arc_count(nodes[2]) != 2
        || graph.arc_count(nodes[3]) != 2
    {
        return false;
    }

    let mut convex1 = 0u32;
    let mut convex2 = 0u32;
    let mut length1 = 0i64;
    let mut length2 = 0i64;
    let f1 = detect_bubble_end_with_bubbles(
        graph,
        nodes[0] ^ 1,
        nodes[3] ^ 1,
        &mut convex1,
        &mut length1,
        None,
    );
    let f2 = detect_bubble_end_with_bubbles(
        graph,
        nodes[1] ^ 1,
        nodes[2] ^ 1,
        &mut convex2,
        &mut length2,
        None,
    );

    if f1 && f2 {
        length1 > minimum_length_threshold && length2 > minimum_length_threshold
    } else if f1 {
        length1 > minimum_length_threshold
    } else if f2 {
        length2 > minimum_length_threshold
    } else {
        false
    }
}

/// `asg_arc_identify_simple_bubbles_multi`: marks bubble vertices and, when
/// asked, bubble-crossing vertices in `sequence_visited`.
///
/// The reference dispatches the worker through `kt_for`; with the default
/// single thread this is the same sequential sweep over all vertices, using
/// the single label buffer.
pub fn asg_arc_identify_simple_bubbles_multi(
    graph: &mut Asg,
    labels: &mut BubbleLabel,
    check_cross: bool,
    minimum_length_threshold: i64,
) -> u64 {
    let start_time = logging::realtime_seconds();
    let vertex_count = graph.sequences.len() * 2;
    if graph.sequence_visited.len() < vertex_count {
        graph.sequence_visited.resize(vertex_count, 0);
    }
    graph.sequence_visited[..vertex_count].fill(0);

    let bubble_distance = get_s_bub_pop_max_dist_advance(graph, &mut labels.buffers[0]);
    labels.reset(graph.sequences.len(), bubble_distance, check_cross as u32);

    for v in 0..vertex_count as u32 {
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if graph.arc_count(v) < 2 || get_real_length(graph, v) < 2 {
            continue;
        }

        if graph.sequence_visited[v as usize] != 1 {
            let buffer = &mut labels.buffers[0];
            if asg_bub_pop1_label(graph, v, bubble_distance, buffer) != 0 {
                let end = buffer.sources[0];
                for i in 0..buffer.visited.len() {
                    let node = buffer.visited[i];
                    if node == v || node == end {
                        continue;
                    }
                    graph.sequence_visited[node as usize] = 1;
                    graph.sequence_visited[(node ^ 1) as usize] = 1;
                }
                graph.sequence_visited[v as usize] = 1;
                graph.sequence_visited[(end ^ 1) as usize] = 1;
            }
        }

        if check_cross
            && graph.sequence_visited[v as usize] == 0
            && check_if_cross(graph, v, minimum_length_threshold)
        {
            graph.sequence_visited[v as usize] = 2;
        }
    }

    let mut nodes = 0i64;
    let mut bubble_nodes = 0i64;
    let mut cross_nodes = 0i64;
    for v in 0..vertex_count {
        if graph.sequences[v >> 1].deleted {
            continue;
        }
        nodes += 1;
        if graph.sequence_visited[v] == 1 {
            bubble_nodes += 1;
        }
        if graph.sequence_visited[v] == 2 {
            cross_nodes += 1;
        }
    }

    let _ = nodes;
    eprintln!(
        "[M::asg_arc_identify_simple_bubbles_multi] takes {:.2} s\n\n",
        logging::realtime_seconds() - start_time
    );

    (bubble_nodes + cross_nodes) as u64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::asg::AsgArc;
    use crate::assembly::bubble::{BubbleBuffer, BubbleInfo, BubbleLabel};

    fn arc(u: u32, v: u32, length: u32) -> AsgArc {
        AsgArc {
            ul: (u64::from(u) << 32) | u64::from(length),
            v,
            overlap: 10,
            ..AsgArc::default()
        }
    }

    /// A symmetric four-read diamond: read 0 branches into reads 1 and 2 and
    /// both meet again at read 3.
    fn diamond() -> Asg {
        let mut graph = Asg::new();
        for read in 0..4 {
            graph.set_sequence(read, 100, false);
        }
        for &(u, v, length) in &[
            (0, 2, 10),
            (0, 4, 10),
            (2, 6, 10),
            (4, 6, 10),
            (3, 1, 10),
            (5, 1, 10),
            (7, 3, 10),
            (7, 5, 10),
        ] {
            graph.push_arc(arc(u, v, length));
        }
        graph.cleanup();
        graph
    }

    fn chain() -> Asg {
        let mut graph = Asg::new();
        for read in 0..4 {
            graph.set_sequence(read, 100, false);
        }
        for &(u, v, length) in &[(0, 2, 10), (2, 6, 10), (3, 1, 10), (7, 3, 10)] {
            graph.push_arc(arc(u, v, length));
        }
        graph.cleanup();
        graph
    }

    fn buffer(sequence_count: usize) -> BubbleBuffer {
        BubbleBuffer {
            info: vec![BubbleInfo::default(); sequence_count * 2],
            ..BubbleBuffer::default()
        }
    }

    #[test]
    fn unitig_walk_measures_a_chain() {
        let graph = chain();
        let mut end = 0u32;
        let mut node_len = 0i64;
        let mut base_len = 0i64;
        let mut max_stop_node_len = 0i64;
        let mut max_stop_base_len = 0i64;
        let mut visited = Vec::new();
        let flag = get_unitig_s(
            &graph,
            0,
            &mut end,
            &mut node_len,
            &mut base_len,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            1,
            Some(&mut visited),
        );
        assert_eq!(flag, END_TIPS);
        assert_eq!(end, 6);
        assert_eq!(node_len, 3);
        assert_eq!(base_len, 120);
        assert_eq!(max_stop_base_len, 120);
        assert_eq!(max_stop_node_len, 3);
        assert_eq!(visited, vec![0, 2, 6]);
    }

    #[test]
    fn a_simple_bubble_is_detected() {
        let graph = diamond();
        let mut end = 0u32;
        let mut min_len = 0i64;
        assert!(detect_simple_bubble(&graph, 0, &mut end, &mut min_len, None));
        assert_eq!(end, 6);
        assert_eq!(min_len, 3);
    }

    #[test]
    fn bubble_labeling_stops_at_the_merge_vertex() {
        let graph = diamond();
        let mut workspace = buffer(4);
        assert_eq!(asg_bub_pop1_label(&graph, 0, 100, &mut workspace), 1);
        assert_eq!(workspace.sources[0], 6);
        let mut visited = workspace.visited.clone();
        visited.sort_unstable();
        assert_eq!(visited, vec![2, 4, 6]);
        // Every label is reset on the way out.
        for info in &workspace.info {
            assert!(!info.state);
            assert_eq!(info.distance, 0);
        }
    }

    #[test]
    fn bubble_identification_marks_both_orientations() {
        let mut graph = diamond();
        let mut labels = BubbleLabel::new(1, 4);
        let marked = asg_arc_identify_simple_bubbles_multi(&mut graph, &mut labels, false, 4);
        assert_eq!(marked, 6);
        for vertex in [0u32, 2, 3, 4, 5, 7] {
            assert_eq!(graph.sequence_visited[vertex as usize], 1, "vertex {vertex}");
        }
        for vertex in [1u32, 6] {
            assert_eq!(graph.sequence_visited[vertex as usize], 0, "vertex {vertex}");
        }
    }

    #[test]
    fn subgraph_distance_search_finds_the_diamond() {
        let graph = diamond();
        let mut workspace = buffer(4);
        let distance = get_s_bub_pop_max_dist_advance(&graph, &mut workspace);
        assert!(distance > 0);
    }

    #[test]
    fn popping_a_bubble_relabels_its_nodes_and_removes_its_edges() {
        let mut graph = diamond();
        let mut workspace = buffer(4);

        assert_eq!(
            asg_bub_pop1_primary_trio(
                &mut graph, None, 0, 100, &mut workspace, u32::MAX, u32::MAX, true, None, &[], false,
                false,
            ),
            1
        );
        assert_eq!(workspace.sources[0], 6, "the sink is left in S");

        // Read 3 is the sink (its class is restored), read 2 lies on the
        // optimal parent path and becomes hap, read 1 is only visited and
        // becomes alternate.
        assert_eq!(graph.sequences[2].coverage_class, HAP_LABEL);
        assert_eq!(graph.sequences[1].coverage_class, ALTER_LABEL);
        assert_eq!(graph.sequences[3].coverage_class, 0);

        // The edge pass removes every visited edge, then the sink-to-beg walk
        // restores the optimal path and its reverse complements.
        for (u, v) in [(0u32, 4u32), (4, 6)] {
            assert!(graph.get_arc(u, v).is_some(), "{u}->{v} survives");
            assert!(graph.get_arc(v ^ 1, u ^ 1).is_some(), "the reverse of {u}->{v}");
        }
        for (u, v) in [(0u32, 2u32), (2, 6)] {
            assert!(graph.get_arc(u, v).is_none(), "{u}->{v} must be removed");
            assert!(
                graph.get_arc(v ^ 1, u ^ 1).is_none(),
                "the reverse of {u}->{v} must be removed"
            );
        }
    }

    #[test]
    fn popping_without_is_pop_only_finds_the_bubble() {
        let mut graph = diamond();
        let mut workspace = buffer(4);

        assert_eq!(
            asg_bub_pop1_primary_trio(
                &mut graph, None, 0, 100, &mut workspace, u32::MAX, u32::MAX, false, None, &[],
                false, false,
            ),
            1
        );
        assert_eq!(workspace.sources[0], 6);
        // The graph is untouched: the cleaning passes pop without recording.
        for sequence in &graph.sequences {
            assert_eq!(sequence.coverage_class, 0);
            assert!(!sequence.deleted);
        }
        assert!(graph.get_arc(0, 4).is_some());
    }
}

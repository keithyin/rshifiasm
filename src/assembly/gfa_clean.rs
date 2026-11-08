//! String-graph cleaning passes, the counterparts of the `asg_arc_cut_*`
//! family, the bubble-identification helpers and `ul_clean_gfa` in
//! `gfa_ut.cpp`/`Overlaps.cpp`.
//!
//! The reference passes a shared `asg64_v` scratch buffer into most of these
//! functions; here the caller owns a `Vec<u64>` and the functions clear or
//! truncate it exactly where the reference resets `b->n`.

use crate::assembly::asg::{Asg, ASG_ET_MERGEABLE};
use crate::assembly::bubble::{BubbleBuffer, BubbleLabel};
use crate::assembly::bubble_identify::{
    asg_arc_identify_simple_bubbles_multi, asg_bub_pop1_primary_trio,
};
use crate::assembly::hit::ReadToUnitig;
use crate::container::ksort::radix_sort_u64;
use crate::correction::ec_driver::PafList;
use crate::reads::read_store::ReadStore;

/// Path classifications shared by the topology walkers, from `Overlaps.h`.
pub const LONG_TIPS: u8 = 0;
pub const TWO_INPUT: u8 = 1;
pub const TWO_OUTPUT: u8 = 2;
pub const MUL_INPUT: u8 = 3;
pub const MUL_OUTPUT: u8 = 4;
pub const END_TIPS: u8 = 5;
pub const LONG_TIPS_UNDER_MAX_EXT: u8 = 6;
pub const LOOP: u8 = 7;

/// `LIM_LEN`, the limit `asg_iterative_semi_circ` passes to
/// [`asg_cut_semi_circ`].
pub const LIM_LEN: u32 = 100;

/// `TRIO_THRES` (`Overlaps.h` 29): a tip belongs to a haplotype only when the
/// dominant label covers more than this share of its flagged reads.
const TRIO_THRES: f64 = 0.9;
/// `DOUBLE_CHECK_THRES` (`Overlaps.h` 30): the same tip must also cover at
/// least this share of its whole walk.
const DOUBLE_CHECK_THRES: f64 = 0.1;

/// `telo_end_t`: which reads carry a telomere motif on either end.
///
/// Bit 0 marks a 5' telomere and bit 1 a 3' telomere; the cleaning passes only
/// test the byte for being nonzero. The reference generates this table from
/// `--telo-motif`; the default run passes `NULL`, which is `None` here.
#[derive(Debug, Clone, Default)]
pub struct TeloEnd {
    /// One byte per read, `hh`.
    pub hh: Vec<u8>,
}

impl TeloEnd {
    /// Whether read `read` carries a telomere motif on either end.
    pub fn is_telomere(&self, read: u32) -> bool {
        self.hh.get(read as usize).copied().unwrap_or(0) != 0
    }
}

/// `get_real_length`: the number of live outgoing arcs of `v`.
pub fn get_real_length(graph: &Asg, v: u32) -> u32 {
    graph.arcs_of(v).iter().filter(|arc| !arc.deleted).count() as u32
}

/// `get_real_length` with its optional target output.
///
/// Every call site either passes `NULL` or has already established that `v`
/// has exactly one live arc, so only the first live target is returned.
pub fn get_real_length_one(graph: &Asg, v: u32) -> Option<u32> {
    graph
        .arcs_of(v)
        .iter()
        .find(|arc| !arc.deleted)
        .map(|arc| arc.v)
}

/// `get_tip_trio_infor` (`Overlaps.cpp` 6890): walks the unitary chain that
/// starts at `begin_node` and classifies it by the trio labels of the reads
/// it covers.
///
/// The walk advances while the current vertex has exactly one live arc whose
/// target has exactly one live arc into the reverse complement, stopping at
/// the first branching vertex (or when the walk returns to `begin_node`). A
/// chain is labelled with the dominant label only when it covers more than
/// [`TRIO_THRES`] of the flagged reads *and* at least
/// [`DOUBLE_CHECK_THRES`] of the whole walk; otherwise it stays ambiguous.
pub fn get_tip_trio_infor(graph: &Asg, begin_node: u32, trio_labels: &[u8]) -> u8 {
    let mut v = begin_node;
    let mut walk_length = 0u32;
    let mut father_occurrence = 0u32;
    let mut mother_occurrence = 0u32;
    // The reference counts `ambigious_occ` but never reads it; keep the
    // counter for the sake of the mirror.
    #[allow(unused_assignments)]
    let mut ambiguous_occurrence = 0u32;

    loop {
        let live_arcs = get_real_length(graph, v);
        walk_length += 1;

        let label = trio_labels[(v >> 1) as usize];
        if label == crate::command_line::options::TRIO_FATHER {
            father_occurrence += 1;
        } else if label == crate::command_line::options::TRIO_MOTHER {
            mother_occurrence += 1;
        } else if label == crate::command_line::options::TRIO_UNKNOWN
            || u32::from(label) == crate::assembly::trio_clean::DROP
        {
            ambiguous_occurrence += 1;
        }

        if live_arcs != 1 {
            break;
        }
        // `live_arcs == 1`, so the single live arc exists; the reference
        // refetches it through the optional output of `get_real_length`.
        let w = get_real_length_one(graph, v).expect("the single live arc exists");
        if get_real_length(graph, w ^ 1) != 1 {
            break;
        }
        v = w;
        if v == begin_node {
            break;
        }
    }

    let flagged_length = father_occurrence + mother_occurrence;
    if flagged_length == 0 {
        return crate::command_line::options::TRIO_UNKNOWN;
    }
    if father_occurrence >= mother_occurrence {
        if (father_occurrence as f64) > TRIO_THRES * (flagged_length as f64)
            && (father_occurrence as f64) >= DOUBLE_CHECK_THRES * (walk_length as f64)
        {
            return crate::command_line::options::TRIO_FATHER;
        }
    } else if (mother_occurrence as f64) > TRIO_THRES * (flagged_length as f64)
        && (mother_occurrence as f64) >= DOUBLE_CHECK_THRES * (walk_length as f64)
    {
        return crate::command_line::options::TRIO_MOTHER;
    }
    crate::command_line::options::TRIO_UNKNOWN
}

/// The `for (i = kv = 0; ...) { if (av[i].del) continue; kv++; break; }` test
/// used throughout the cleaning passes.
fn has_live_arc(graph: &Asg, v: u32) -> bool {
    graph.arcs_of(v).iter().any(|arc| !arc.deleted)
}

/// `asg_arc_cut_tips`: removes the chains that end within `max_ext` merges.
///
/// A vertex is a candidate when it has no live outgoing arc on its opposite
/// orientation, i.e. it starts a chain. The chain is walked through
/// [`Asg::is_unitig_end`] up to `max_ext` merges; when it stops earlier (a
/// tip), every read on the chain is deleted. Candidates are processed
/// shortest-first, and the graph is re-checked before each deletion because
/// earlier cuts may already have removed part of the chain.
///
/// `is_ou`/`read_to_unitig` select the ultra-long branch (`ru && is_ou`),
/// `telo_end` the telomere-aware branch; the default HiFi run uses neither.
pub fn asg_arc_cut_tips(
    graph: &mut Asg,
    max_ext: u32,
    buffer: &mut Vec<u64>,
    is_ou: bool,
    read_to_unitig: Option<&ReadToUnitig>,
    telo_end: Option<&TeloEnd>,
) -> u32 {
    let vertex_count = graph.sequences.len() * 2;
    let mut count = 0u32;
    buffer.clear();

    for v in 0..vertex_count as u32 {
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if telo_end.is_some_and(|te| te.is_telomere(v >> 1)) {
            continue;
        }
        if has_live_arc(graph, v ^ 1) {
            continue;
        }

        // `kv` starts at 1 for `v` itself; `mm_ou` is the smallest occurrence
        // seen along the chain, but the reference accumulates it even when
        // occurrence is not tracked, in which case it is zero.
        let mut kv = 1u32;
        let mut max_occurrence = u32::MAX;
        let mut is_telo = false;
        let mut i = 0u32;
        let mut w = v;
        while i < max_ext {
            let (kind, packed, occurrence) = graph.is_unitig_end_with_occurrence(w ^ 1);
            if kind != ASG_ET_MERGEABLE {
                break;
            }
            let lw = packed.expect("a mergeable edge always reports its endpoint");
            w = lw as u32;
            kv = kv.wrapping_add(1);
            max_occurrence = if is_ou {
                max_occurrence.min(occurrence)
            } else {
                0
            };
            if telo_end.is_some_and(|te| te.is_telomere(w >> 1)) {
                is_telo = true;
            }
            i += 1;
        }

        if max_occurrence == u32::MAX {
            max_occurrence = 0;
        }
        kv = kv.wrapping_add(max_occurrence);
        i = i.wrapping_add(max_occurrence);
        if i < max_ext && !is_telo {
            buffer.push((u64::from(kv) << 32) | u64::from(v));
        }
    }

    radix_sort_u64(buffer);

    let sorted_length = buffer.len();
    for k in 0..sorted_length {
        let v = buffer[k] as u32;

        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if telo_end.is_some_and(|te| te.is_telomere(v >> 1)) {
            continue;
        }
        if has_live_arc(graph, v ^ 1) {
            continue;
        }

        let previous_length = buffer.len();
        buffer.push(u64::from(v));
        let mut max_occurrence = u32::MAX;
        let mut is_telo = false;
        let mut i = 0u32;
        let mut w = v;
        while i < max_ext {
            let (kind, packed, occurrence) = graph.is_unitig_end_with_occurrence(w ^ 1);
            if kind != ASG_ET_MERGEABLE {
                break;
            }
            let lw = packed.expect("a mergeable edge always reports its endpoint");
            w = lw as u32;
            buffer.push(lw);
            max_occurrence = if is_ou {
                max_occurrence.min(occurrence)
            } else {
                0
            };
            if telo_end.is_some_and(|te| te.is_telomere(w >> 1)) {
                is_telo = true;
            }
            i += 1;
        }

        if max_occurrence == u32::MAX {
            max_occurrence = 0;
        }
        i = i.wrapping_add(max_occurrence);

        if i < max_ext && !is_telo {
            for entry in &buffer[previous_length..] {
                let read = (*entry as u32) >> 1;
                graph.delete_sequence(read);
            }
            count += 1;
        }
        buffer.truncate(previous_length);
    }

    if let (Some(read_to_unitig), true) = (read_to_unitig, is_ou) {
        // The reference keeps `b->n` across the whole loop, so every deletion
        // pass re-walks the entries accumulated so far; reproduced here.
        buffer.clear();
        for v in 0..vertex_count as u32 {
            if graph.sequences[(v >> 1) as usize].deleted {
                continue;
            }
            if has_live_arc(graph, v ^ 1) {
                continue;
            }

            let (unitig, is_unitig) = read_to_unitig.get(v >> 1);
            if unitig == u32::MAX || is_unitig == 1 {
                continue;
            }
            if telo_end.is_some_and(|te| te.is_telomere(v >> 1)) {
                continue;
            }
            buffer.push(u64::from(v));

            let mut i = 0u32;
            let mut w = v;
            while i < max_ext {
                let (kind, packed, _) = graph.is_unitig_end_with_occurrence(w ^ 1);
                if kind != ASG_ET_MERGEABLE {
                    break;
                }
                let lw = packed.expect("a mergeable edge always reports its endpoint");
                w = lw as u32;

                let (unitig, is_unitig) = read_to_unitig.get(w >> 1);
                if unitig == u32::MAX || is_unitig == 1 {
                    break;
                }
                if telo_end.is_some_and(|te| te.is_telomere(w >> 1)) {
                    break;
                }

                buffer.push(lw);
                i += 1;
            }

            for entry in buffer.iter() {
                let read = (*entry as u32) >> 1;
                graph.delete_sequence(read);
            }
            if !buffer.is_empty() {
                count += 1;
            }
        }
    }

    if count > 0 {
        graph.cleanup();
    }

    count
}

/// `get_arcs`: the number of live arcs of `v`, plus their raw indices in
/// `graph.arcs` when a slot array is supplied.
///
/// The reference returns C array indices (`av - g->arc + i`); the same
/// indices address `graph.arcs` here.
pub fn get_arcs(graph: &Asg, v: u32, indices: Option<&mut [u32]>) -> u32 {
    let (start, count) = graph.arc_range(v);
    let mut slots = indices;
    let mut live = 0u32;
    for i in 0..count {
        if graph.arcs[start + i].deleted {
            continue;
        }
        if let Some(slots) = slots.as_deref_mut() {
            if (live as usize) < slots.len() {
                slots[live as usize] = (start + i) as u32;
            }
        }
        live += 1;
    }
    live
}

/// `follow_limit_path`: walks the unambiguous chain starting at `s`.
///
/// `end` receives the last vertex reached and `occurrence` the number of
/// vertices visited. The returned code describes what stopped the walk.
pub fn follow_limit_path(
    graph: &Asg,
    s: u32,
    end: &mut u32,
    occurrence: &mut u32,
    mut buffer: Option<&mut Vec<u64>>,
    limit: u32,
) -> u8 {
    let mut v = s;
    *occurrence = 0;

    loop {
        *occurrence = occurrence.wrapping_add(1);
        let mut arc_index = [0u32; 1];
        let kv = get_arcs(graph, v, Some(&mut arc_index));
        *end = v;

        if let Some(buffer) = buffer.as_deref_mut() {
            buffer.push(u64::from(v));
        }

        if kv == 0 {
            return END_TIPS;
        }
        if kv == 2 {
            return TWO_OUTPUT;
        }
        if kv > 2 {
            return MUL_OUTPUT;
        }
        if *occurrence > limit {
            return LONG_TIPS;
        }
        let w = graph.arcs[arc_index[0] as usize].v;
        // Up to here `kv == 1`, so the next vertex always has at least one
        // incoming arc.
        let kw = get_arcs(graph, w ^ 1, None);
        v = w;

        if kw == 2 {
            return TWO_INPUT;
        }
        if kw > 2 {
            return MUL_INPUT;
        }
        if v == s {
            return LOOP;
        }
    }
}

/// `if_sup_chimeric`: whether the left- and right-anchored overlaps of a read
/// fail to support its middle.
///
/// Returns the reference's 1 for a chimeric read and 0 otherwise. `buffer`
/// is the shared scratch vector; its length is restored before returning.
pub fn if_sup_chimeric(
    source: &PafList,
    read_length: u64,
    buffer: &mut Vec<u64>,
    if_exact: bool,
) -> i32 {
    let initial = read_length as u32;
    let mut left = [initial, 0u32];
    let mut right = [initial, 0u32];

    for hit in &source.hits {
        if hit.deleted {
            continue;
        }
        if if_exact && !hit.exact_region {
            continue;
        }

        let query_start = hit.query_start;
        let query_end = hit.query_end;

        if query_start == 0 {
            if query_start < left[0] {
                left[0] = query_start;
            }
            if query_end > left[1] {
                left[1] = query_end;
            }
        }

        if u64::from(query_end) == read_length {
            if query_start < right[0] {
                right[0] = query_start;
            }
            if query_end > right[1] {
                right[1] = query_end;
            }
        }
    }

    if left[1] > right[0] {
        return 0;
    }
    if left[1] <= left[0] || right[1] <= right[0] {
        return 1;
    }

    let previous_length = buffer.len();
    if left[1] > left[0] {
        buffer.push(u64::from(left[0] << 1));
        buffer.push(u64::from((left[1] << 1) | 1));
    }
    if right[1] > right[0] {
        buffer.push(u64::from(right[0] << 1));
        buffer.push(u64::from((right[1] << 1) | 1));
    }

    // The contained overlaps fill whatever gap the two anchored sides leave.
    for hit in &source.hits {
        if hit.deleted {
            continue;
        }
        if if_exact && !hit.exact_region {
            continue;
        }

        let query_start = hit.query_start;
        let query_end = hit.query_end;
        if query_start == 0 || u64::from(query_end) == read_length {
            continue;
        }

        buffer.push(u64::from(query_start << 1));
        buffer.push(u64::from((query_end << 1) | 1));
    }

    radix_sort_u64(&mut buffer[previous_length..]);
    left = [initial, 0u32];
    right = [initial, 0u32];

    let mut depth = 0i32;
    let mut state_start = 0u32;
    for k in previous_length..buffer.len() {
        let previous_depth = depth;
        if buffer[k] & 1 != 0 {
            depth -= 1;
        } else {
            depth += 1;
        }

        if previous_depth < 1 && depth >= 1 {
            state_start = (buffer[k] >> 1) as u32;
        } else if previous_depth >= 1 && depth < 1 {
            if state_start == 0 {
                left[0] = state_start;
                left[1] = (buffer[k] >> 1) as u32;
            }
            if buffer[k] >> 1 == read_length {
                right[0] = state_start;
                right[1] = (buffer[k] >> 1) as u32;
            }
        }
    }

    buffer.truncate(previous_length);
    if left[1] > right[0] {
        return 0;
    }

    1
}

/// `asg_check_unambi1`: the target of the only live arc of `v`, or
/// `u32::MAX` when `v` does not have exactly one.
fn asg_check_unambi1(graph: &Asg, v: u32) -> u32 {
    let mut target = u32::MAX;
    let mut live = 0u32;
    for arc in graph.arcs_of(v) {
        if !arc.deleted {
            live += 1;
            target = arc.v;
        }
    }
    if live != 1 {
        return u32::MAX;
    }
    target
}

/// `asg_topocut_aux`: how far the unambiguous chain from `v` extends.
pub fn asg_topocut_aux(graph: &Asg, v: u32, max_ext: i32) -> i32 {
    let mut v = v;
    let mut n_ext = 1i32;
    while n_ext < max_ext && v != u32::MAX {
        if asg_check_unambi1(graph, v ^ 1) == u32::MAX {
            n_ext -= 1;
            break;
        }
        v = asg_check_unambi1(graph, v);
        n_ext += 1;
    }
    n_ext
}

/// `asg_cut_semi_circ`: cuts the input edge of a semi-circular path.
///
/// A vertex with several inputs but a single output starts a "semi-circle"
/// when its forward path ends inside `lim_len` at a branch that is not a tip
/// and not a loop; the edge `v^1 -> e^1` (and its mirror) is then removed.
pub fn asg_cut_semi_circ(graph: &mut Asg, lim_len: u32, is_clean: bool) -> u32 {
    let vertex_count = graph.sequences.len() * 2;
    let mut count = 0u32;

    for v in 0..vertex_count as u32 {
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if get_real_length(graph, v ^ 1) <= 1 {
            continue;
        }
        if get_real_length(graph, v) != 1 {
            continue;
        }

        // The reference loops over the raw arcs and breaks after the first
        // live one; `get_real_length(v) == 1` guarantees a single iteration.
        let mut end = 0u32;
        let mut occurrence = 0u32;
        let kind = follow_limit_path(graph, v, &mut end, &mut occurrence, None, lim_len);
        if occurrence > lim_len || kind == LONG_TIPS || kind == LOOP || kind == END_TIPS {
            continue;
        }

        let target = end ^ 1;
        let (start, _) = graph.arc_range(v ^ 1);
        let mut deleted_indices = Vec::new();
        for i in 0..graph.arc_count(v ^ 1) {
            let index = start + i as usize;
            if !graph.arcs[index].deleted && graph.arcs[index].v == target {
                deleted_indices.push(index);
            }
        }
        for index in deleted_indices {
            graph.arcs[index].deleted = true;
            let arc = graph.arcs[index];
            graph.delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, true);
            count += 1;
        }
    }

    if count > 0 && is_clean {
        graph.cleanup();
    }
    count
}

/// `asg_cut_chimeric_bub`: removes a chimeric read that joins two branches.
///
/// The read must have one live arc in each direction, its two neighbours must
/// each branch, its overlaps must be unsupported, and the alternative path
/// must be long enough and end at the read's other branch.
pub fn asg_cut_chimeric_bub(
    graph: &mut Asg,
    sources: &[PafList],
    buffer: &mut Vec<u64>,
    normal_len: u32,
    is_clean: bool,
    telo_end: Option<&TeloEnd>,
) -> u32 {
    let vertex_count = graph.sequences.len() * 2;
    let mut count = 0u32;
    buffer.clear();

    for v in 0..vertex_count as u32 {
        let read = (v >> 1) as usize;
        if graph.sequences[read].deleted {
            continue;
        }
        if telo_end.is_some_and(|te| te.is_telomere(v >> 1)) {
            continue;
        }

        let mut first_arc = [0u32; 1];
        let mut reverse_arc = [0u32; 1];
        if get_arcs(graph, v, Some(&mut first_arc)) != 1
            || get_arcs(graph, v ^ 1, Some(&mut reverse_arc)) != 1
        {
            continue;
        }

        let first_target = graph.arcs[first_arc[0] as usize].v;
        let reverse_target = graph.arcs[reverse_arc[0] as usize].v;
        if get_arcs(graph, first_target ^ 1, None) != 2
            || get_arcs(graph, reverse_target ^ 1, None) != 2
        {
            continue;
        }

        if if_sup_chimeric(
            &sources[read],
            u64::from(graph.sequences[read].length),
            buffer,
            true,
        ) == 0
        {
            continue;
        }

        // The other branch of the first neighbour, i.e. the live arc that is
        // not its edge back to `v`.
        let branch = first_target ^ 1;
        let mut branch_target = u32::MAX;
        for arc in graph.arcs_of(branch) {
            if arc.deleted {
                continue;
            }
            if arc.v == (v ^ 1) {
                continue;
            }
            branch_target = arc.v;
            break;
        }
        let branch_target = branch_target;
        if branch_target == u32::MAX {
            continue;
        }

        let mut end = 0u32;
        let mut occurrence = 0u32;
        if follow_limit_path(graph, branch_target, &mut end, &mut occurrence, None, u32::MAX)
            != TWO_INPUT
        {
            continue;
        }

        if occurrence > normal_len {
            let mut far_target = u32::MAX;
            let mut live = 0u32;
            for arc in graph.arcs_of(end) {
                if arc.deleted {
                    continue;
                }
                live += 1;
                far_target = arc.v;
                if live > 1 {
                    break;
                }
            }
            if live == 1 && far_target == reverse_target {
                graph.delete_sequence((v >> 1) as u32);
                count += 1;
            }
        }
    }

    if is_clean && count > 0 {
        graph.cleanup();
    }
    count
}

/// `asg_iterative_semi_circ`: repeats the semi-circle and chimeric-bubble
/// cuts until neither removes anything.
pub fn asg_iterative_semi_circ(
    graph: &mut Asg,
    sources: &[PafList],
    buffer: &mut Vec<u64>,
    normal_len: u32,
    pop_chimeric: bool,
    telo_end: Option<&TeloEnd>,
) {
    let mut total = 0u64;
    let mut step = 1u64;
    while step != 0 {
        step = u64::from(asg_cut_semi_circ(graph, LIM_LEN, false));
        if pop_chimeric {
            step += u64::from(asg_cut_chimeric_bub(
                graph,
                sources,
                buffer,
                normal_len,
                false,
                telo_end,
            ));
        }
        total += step;
    }

    if total > 0 {
        graph.cleanup();
    }
}

/// `asg_arc_cut_chimeric`: removes single chimeric reads from the graph.
///
/// A read is removed when it has one arc per direction, both neighbours
/// branch, its first arc is not exact, and its overlaps are unsupported.
pub fn asg_arc_cut_chimeric(
    graph: &mut Asg,
    sources: &[PafList],
    buffer: &mut Vec<u64>,
    ou_threshold: u32,
    telo_end: Option<&TeloEnd>,
) -> u32 {
    let vertex_count = graph.sequences.len() * 2;
    let mut count = 0u32;
    buffer.clear();

    for v in 0..vertex_count as u32 {
        let read = (v >> 1) as usize;
        if graph.sequences[read].deleted {
            continue;
        }
        if telo_end.is_some_and(|te| te.is_telomere(v >> 1)) {
            continue;
        }
        if graph.sequence_visited[v as usize] != 0 {
            continue;
        }

        let mut first_arc = [0u32; 1];
        let mut reverse_arc = [0u32; 1];
        if get_arcs(graph, v, Some(&mut first_arc)) != 1
            || get_arcs(graph, v ^ 1, Some(&mut reverse_arc)) != 1
        {
            continue;
        }

        let first = graph.arcs[first_arc[0] as usize];
        let reverse = graph.arcs[reverse_arc[0] as usize];
        if get_arcs(graph, first.v ^ 1, None) < 2 || get_arcs(graph, reverse.v ^ 1, None) < 2 {
            continue;
        }
        if first.exact_region != 0 {
            continue;
        }
        if ou_threshold != u32::MAX
            && u32::from(first.occurrence) >= ou_threshold
            && u32::from(reverse.occurrence) >= ou_threshold
        {
            continue;
        }
        if if_sup_chimeric(
            &sources[read],
            u64::from(graph.sequences[read].length),
            buffer,
            true,
        ) == 0
        {
            continue;
        }

        buffer.push((u64::from(first.overlap) << 32) | u64::from(first_arc[0]));
    }

    radix_sort_u64(buffer);

    let sorted_length = buffer.len();
    for k in 0..sorted_length {
        let raw_index = buffer[k] as u32 as usize;
        if graph.arcs[raw_index].deleted {
            continue;
        }

        let v = (graph.arcs[raw_index].ul >> 32) as u32;
        let w = graph.arcs[raw_index].v ^ 1;
        if graph.sequences[(v >> 1) as usize].deleted
            || graph.sequences[(w >> 1) as usize].deleted
        {
            continue;
        }

        let mut first_arc = [0u32; 1];
        let mut reverse_arc = [0u32; 1];
        if get_arcs(graph, v, Some(&mut first_arc)) != 1
            || get_arcs(graph, v ^ 1, Some(&mut reverse_arc)) != 1
        {
            continue;
        }
        let first = graph.arcs[first_arc[0] as usize];
        let reverse = graph.arcs[reverse_arc[0] as usize];
        if get_arcs(graph, first.v ^ 1, None) < 2 || get_arcs(graph, reverse.v ^ 1, None) < 2 {
            continue;
        }

        let mut has_exact_other = false;
        for arc in graph.arcs_of(w) {
            if arc.deleted || arc.v == (v ^ 1) || arc.exact_region == 0 {
                continue;
            }
            has_exact_other = true;
            break;
        }
        if !has_exact_other {
            continue;
        }
        if telo_end.is_some_and(|te| te.is_telomere(v >> 1)) {
            continue;
        }

        graph.delete_sequence(v >> 1);
        count += 1;
    }

    if count > 0 {
        graph.cleanup();
    }
    count
}

/// `asg_arc_cut_inexact`: removes an inexact edge that competes with a better
/// supported path.
///
/// With `is_trio` set (the `-3/-4`/`-1/-2` front ends), an arc whose target
/// tip belongs to the opposite haplotype of the split never counts as the
/// better supported path ([`get_tip_trio_infor`]).
#[allow(clippy::too_many_arguments)]
pub fn asg_arc_cut_inexact(
    graph: &mut Asg,
    sources: &[PafList],
    buffer: &mut Vec<u64>,
    max_ext: i32,
    is_ou: bool,
    is_trio: bool,
    min_diff: u32,
    occurrence_ratio: f32,
    trio_labels: &[u8],
) -> u32 {

    let vertex_count = graph.sequences.len() * 2;
    let mut count = 0u32;
    buffer.clear();

    for v in 0..vertex_count as u32 {
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if graph.sequence_visited[v as usize] != 0 {
            continue;
        }

        let (start, arc_count) = graph.arc_range(v);
        if arc_count < 2 {
            continue;
        }
        let live = graph.arcs[start..start + arc_count]
            .iter()
            .filter(|arc| !arc.deleted)
            .count();
        if live < 2 {
            continue;
        }

        for i in 0..arc_count {
            let arc = graph.arcs[start + i as usize];
            if arc.deleted || arc.exact_region != 0 {
                continue;
            }
            buffer.push((u64::from(arc.overlap) << 32) | (start + i as usize) as u64);
        }
    }

    radix_sort_u64(buffer);

    let sorted_length = buffer.len();
    for k in 0..sorted_length {
        let raw_index = buffer[k] as u32 as usize;
        let ve = graph.arcs[raw_index];
        if ve.deleted {
            continue;
        }

        let v = (ve.ul >> 32) as u32;
        let w = ve.v ^ 1;
        if graph.sequences[(v >> 1) as usize].deleted
            || graph.sequences[(w >> 1) as usize].deleted
        {
            continue;
        }

        let nv = graph.arc_count(v);
        let nw = graph.arc_count(w);
        if nv <= 1 && nw <= 1 {
            continue;
        }

        // `if(is_trio) { speedup; trioF = get_tip_trio_infor(g, v^1);
        // ntrioF = opposite }`: the haplotype an arc must not be matched by.
        let mut opposite_flag = u32::MAX;
        if is_trio {
            if get_real_length(graph, v) <= 1 && get_real_length(graph, w) <= 1 {
                continue;
            }
            let trio_flag = get_tip_trio_infor(graph, v ^ 1, trio_labels);
            opposite_flag = match trio_flag {
                crate::command_line::options::TRIO_FATHER => {
                    u32::from(crate::command_line::options::TRIO_MOTHER)
                }
                crate::command_line::options::TRIO_MOTHER => {
                    u32::from(crate::command_line::options::TRIO_FATHER)
                }
                _ => u32::MAX,
            };
        }

        // The mirror arc `w -> v^1`; the graph is kept symmetric, so it is
        // present even when it has been deleted.
        let (w_start, _) = graph.arc_range(w);
        let we_index = graph.arcs[w_start..w_start + graph.arc_count(w) as usize]
            .iter()
            .position(|arc| arc.v == (v ^ 1))
            .map(|offset| w_start + offset)
            .expect("a symmetric graph stores the mirror arc");
        let we = graph.arcs[we_index];

        let mutual_overlap = ve.overlap.min(we.overlap);
        let mutual_occurrence = u32::from(ve.occurrence).min(u32::from(we.occurrence));

        let mut kv = 0u32;
        let mut v_max_overlap = 0u32;
        let mut v_max_occurrence = 0u32;
        let mut v_max: Option<(usize, crate::assembly::asg::AsgArc)> = None;
        let (v_start, _) = graph.arc_range(v);
        for i in 0..graph.arc_count(v) {
            let arc = graph.arcs[v_start + i as usize];
            if arc.deleted {
                continue;
            }
            kv += 1;
            if is_trio
                && u32::from(get_tip_trio_infor(graph, arc.v, trio_labels)) == opposite_flag
            {
                continue;
            }
            if v_max_overlap < arc.overlap {
                v_max_overlap = arc.overlap;
                v_max = Some((v_start + i as usize, arc));
            }
            if v_max_occurrence < u32::from(arc.occurrence) {
                v_max_occurrence = u32::from(arc.occurrence);
            }
        }
        if kv < 1 {
            continue;
        }
        if kv >= 2 {
            if mutual_overlap >= v_max_overlap {
                continue;
            }
            if is_ou && (mutual_occurrence as f32) > (v_max_occurrence as f32) * occurrence_ratio {
                continue;
            }
            if mutual_overlap.wrapping_add(min_diff) > v_max_overlap {
                continue;
            }
        }

        let mut kw = 0u32;
        let mut w_max_overlap = 0u32;
        let mut w_max_occurrence = 0u32;
        for i in 0..graph.arc_count(w) {
            let arc = graph.arcs[w_start + i as usize];
            if arc.deleted {
                continue;
            }
            kw += 1;
            if is_trio
                && u32::from(get_tip_trio_infor(graph, arc.v, trio_labels)) == opposite_flag
            {
                continue;
            }
            if w_max_overlap < arc.overlap {
                w_max_overlap = arc.overlap;
            }
            if w_max_occurrence < u32::from(arc.occurrence) {
                w_max_occurrence = u32::from(arc.occurrence);
            }
        }
        if kw < 1 {
            continue;
        }
        if kw >= 2 {
            if mutual_overlap >= w_max_overlap {
                continue;
            }
            if is_ou && (mutual_occurrence as f32) > (w_max_occurrence as f32) * occurrence_ratio {
                continue;
            }
            if mutual_overlap.wrapping_add(min_diff) > w_max_overlap {
                continue;
            }
        }
        if kv <= 1 && kw <= 1 {
            continue;
        }

        let mut to_delete = false;
        let v_corrected = sources[(v >> 1) as usize].is_fully_corrected;
        let w_corrected = sources[(w >> 1) as usize].is_fully_corrected;

        if v_corrected && w_corrected {
            if kv > 1 && kw > 1 {
                to_delete = true;
            } else if kw == 1 {
                if asg_topocut_aux(graph, w ^ 1, max_ext) < max_ext {
                    to_delete = true;
                }
            } else if kv == 1 && asg_topocut_aux(graph, v ^ 1, max_ext) < max_ext {
                to_delete = true;
            }
        }

        if v_corrected && !w_corrected {
            if let Some((_, v_max_arc)) = v_max {
                if v_max_arc.v != ve.v
                    && v_max_arc.exact_region != 0
                    && sources[(v_max_arc.v >> 1) as usize].is_fully_corrected
                {
                    if kv > 1 && kw > 1 {
                        to_delete = true;
                    } else if kw == 1 {
                        if asg_topocut_aux(graph, w ^ 1, max_ext) < max_ext {
                            to_delete = true;
                        }
                    } else if kv == 1 && asg_topocut_aux(graph, v ^ 1, max_ext) < max_ext {
                        to_delete = true;
                    }
                }
            }
        }

        if to_delete {
            graph.arcs[raw_index].deleted = true;
            graph.arcs[we_index].deleted = true;
            count += 1;
        }
    }

    if count > 0 {
        graph.cleanup();
    }
    count
}


/// `trans_path_check`: whether two edges leave towards the same downstream
/// read set, judged from the stored reverse overlaps.
///
/// Returns `u32::MAX` for the reference's -1.
#[allow(clippy::too_many_arguments)]
pub fn trans_path_check(
    graph: &mut Asg,
    first: u32,
    second: u32,
    reverse_sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    minimum_length: u32,
    buffer: &mut Vec<u64>,
) -> u32 {
    if first == second {
        return u32::MAX;
    }

    let initial_length = buffer.len();
    let mut end = 0u32;
    let mut first_length = 0u32;
    let mut second_length = 0u32;

    let kind = follow_limit_path(
        graph,
        first,
        &mut end,
        &mut first_length,
        Some(buffer),
        u32::MAX,
    );
    if kind == LOOP || first_length <= minimum_length {
        buffer.truncate(initial_length);
        return u32::MAX;
    }

    let kind = follow_limit_path(
        graph,
        second,
        &mut end,
        &mut second_length,
        Some(buffer),
        u32::MAX,
    );
    if kind == LOOP || second_length <= minimum_length {
        buffer.truncate(initial_length);
        return u32::MAX;
    }

    let first_start = initial_length;
    let second_start = initial_length + first_length as usize;
    let (short_start, short_length, long_start, long_length) = if first_length > second_length {
        (second_start, second_length, first_start, first_length)
    } else {
        (first_start, first_length, second_start, second_length)
    };

    for i in 0..long_length as usize {
        let vertex = buffer[long_start + i] as u32;
        graph.sequence_visited[(vertex >> 1) as usize] = 1;
    }

    let mut minimum_count = 0.0f64;
    let mut maximum_count = 0.0f64;
    for i in 0..short_length as usize {
        let query = (buffer[short_start + i] as u32) >> 1;
        for hit in &reverse_sources[query as usize].hits {
            let mut target = hit.target_identifier;
            if graph.sequences[target as usize].deleted {
                let (unitig, is_unitig) = read_to_unitig.get(target);
                target = unitig;
                if target == u32::MAX
                    || is_unitig == 1
                    || graph.sequences[target as usize].deleted
                {
                    continue;
                }
            }
            minimum_count += 1.0;
            maximum_count += f64::from(graph.sequence_visited[target as usize]);
        }
    }

    for i in 0..long_length as usize {
        let vertex = buffer[long_start + i] as u32;
        graph.sequence_visited[(vertex >> 1) as usize] = 0;
    }
    buffer.truncate(initial_length);

    if minimum_count == 0.0 {
        return u32::MAX;
    }
    if maximum_count == 0.0 {
        return 0;
    }
    if maximum_count / minimum_count > 0.3 {
        return 1;
    }
    0
}

/// `asg_arc_cut_length`: cuts the shorter side of a split whose overlaps are
/// too short relative to the longest alternative.
///
/// With `is_trio` set, an arc whose target tip belongs to the opposite
/// haplotype of the split never counts as the longest alternative
/// ([`get_tip_trio_infor`]).
#[allow(clippy::too_many_arguments)]
pub fn asg_arc_cut_length(
    graph: &mut Asg,
    buffer: &mut Vec<u64>,
    max_ext: i32,
    length_ratio: f32,
    occurrence_ratio: f32,
    is_occurrence: bool,
    is_trio: bool,
    is_topological: bool,
    minimum_difference: u32,
    minimum_occurrence: u32,
    reverse_sources: Option<&[PafList]>,
    read_to_unitig: Option<&ReadToUnitig>,
    maximum_drop_length: Option<u32>,
    trio_labels: &[u8],
) -> u32 {

    buffer.clear();
    let vertex_count = graph.sequences.len() * 2;
    let use_reverse = reverse_sources.is_some() && read_to_unitig.is_some();

    for v in 0..vertex_count as u32 {
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if graph.sequence_visited[v as usize] != 0 {
            continue;
        }
        let (start, count) = graph.arc_range(v);
        if count < 2 {
            continue;
        }
        let mut live = 0u32;
        for arc in &graph.arcs[start..start + count] {
            if !arc.deleted {
                live += 1;
            }
        }
        if live < 2 {
            continue;
        }
        for i in 0..count {
            let arc = graph.arcs[start + i];
            if arc.deleted {
                continue;
            }
            if let Some(limit) = maximum_drop_length {
                if arc.overlap >= limit {
                    continue;
                }
            }
            buffer.push((u64::from(arc.overlap) << 32) | (start + i) as u64);
        }
    }

    if use_reverse {
        graph.sequence_visited.fill(0);
    }
    radix_sort_u64(buffer);

    let mut count_deleted = 0u32;
    for k in 0..buffer.len() {
        let raw_index = buffer[k] as u32 as usize;
        if graph.arcs[raw_index].deleted {
            continue;
        }

        let ve = graph.arcs[raw_index];
        let v = (ve.ul >> 32) as u32;
        let w = ve.v ^ 1;
        if graph.sequences[(v >> 1) as usize].deleted
            || graph.sequences[(w >> 1) as usize].deleted
        {
            continue;
        }

        let nv = graph.arc_count(v);
        let nw = graph.arc_count(w);
        if nv <= 1 && nw <= 1 {
            continue;
        }

        // `if(is_trio) { speedup; trioF = get_tip_trio_infor(g, v^1);
        // ntrioF = opposite }`, as in [`asg_arc_cut_inexact`].
        let mut opposite_flag = u32::MAX;
        if is_trio {
            if get_real_length(graph, v) <= 1 && get_real_length(graph, w) <= 1 {
                continue;
            }
            let trio_flag = get_tip_trio_infor(graph, v ^ 1, trio_labels);
            opposite_flag = match trio_flag {
                crate::command_line::options::TRIO_FATHER => {
                    u32::from(crate::command_line::options::TRIO_MOTHER)
                }
                crate::command_line::options::TRIO_MOTHER => {
                    u32::from(crate::command_line::options::TRIO_FATHER)
                }
                _ => u32::MAX,
            };
        }

        let mirror_offset = graph
            .arcs_of(w)
            .iter()
            .position(|arc| arc.v == (v ^ 1))
            .expect("the reverse arc exists in a symmetric graph");
        let (mirror_start, _) = graph.arc_range(w);
        let we_index = mirror_start + mirror_offset;
        let we = graph.arcs[we_index];
        let mutual_overlap = ve.overlap.min(we.overlap);
        let mutual_occurrence = u32::from(ve.occurrence).min(u32::from(we.occurrence));

        let mut kv = 0u32;
        let mut overlap_max = 0u32;
        let mut occurrence_max = 0u32;
        let mut v_max: Option<usize> = None;
        let (start, count) = graph.arc_range(v);
        for i in 0..count {
            let arc = graph.arcs[start + i];
            if arc.deleted {
                continue;
            }
            kv += 1;
            if is_trio
                && u32::from(get_tip_trio_infor(graph, arc.v, trio_labels)) == opposite_flag
            {
                continue;
            }
            if overlap_max < arc.overlap {
                overlap_max = arc.overlap;
                v_max = Some(start + i);
            }
            if occurrence_max < u32::from(arc.occurrence) {
                occurrence_max = u32::from(arc.occurrence);
            }
        }
        if kv < 1 {
            continue;
        }
        if kv >= 2 {
            if (mutual_overlap as f32) > (overlap_max as f32) * length_ratio {
                continue;
            }
            if is_occurrence
                && (mutual_occurrence as f32) > (occurrence_max as f32) * occurrence_ratio
                && mutual_occurrence > minimum_occurrence
            {
                continue;
            }
            if mutual_overlap.wrapping_add(minimum_difference) > overlap_max {
                continue;
            }
        }

        let mut kw = 0u32;
        let mut overlap_max = 0u32;
        let mut occurrence_max = 0u32;
        let mut w_max: Option<usize> = None;
        let (start, count) = graph.arc_range(w);
        for i in 0..count {
            let arc = graph.arcs[start + i];
            if arc.deleted {
                continue;
            }
            kw += 1;
            if is_trio
                && u32::from(get_tip_trio_infor(graph, arc.v, trio_labels)) == opposite_flag
            {
                continue;
            }
            if overlap_max < arc.overlap {
                overlap_max = arc.overlap;
                w_max = Some(start + i);
            }
            if occurrence_max < u32::from(arc.occurrence) {
                occurrence_max = u32::from(arc.occurrence);
            }
        }
        if kw < 1 {
            continue;
        }
        if kw >= 2 {
            if (mutual_overlap as f32) > (overlap_max as f32) * length_ratio {
                continue;
            }
            if is_occurrence
                && (mutual_occurrence as f32) > (occurrence_max as f32) * occurrence_ratio
                && mutual_occurrence > minimum_occurrence
            {
                continue;
            }
            if mutual_overlap.wrapping_add(minimum_difference) > overlap_max {
                continue;
            }
        }

        if kv <= 1 && kw <= 1 {
            continue;
        }

        let mut to_delete = false;
        if is_topological {
            if kv > 1 && kw > 1 {
                to_delete = true;
            } else if kw == 1 {
                if asg_topocut_aux(graph, w ^ 1, max_ext) < max_ext {
                    to_delete = true;
                }
            } else if kv == 1 && asg_topocut_aux(graph, v ^ 1, max_ext) < max_ext {
                to_delete = true;
            }
        }

        if use_reverse {
            let reverse_sources = reverse_sources.expect("checked above");
            let read_to_unitig = read_to_unitig.expect("checked above");
            let v_max_target = v_max.map(|index| graph.arcs[index].v);
            let w_max_target = w_max.map(|index| graph.arcs[index].v);

            if !to_delete {
                if let Some(target) = v_max_target {
                    if ve.v != target
                        && trans_path_check(
                            graph,
                            ve.v,
                            target,
                            reverse_sources,
                            read_to_unitig,
                            max_ext as u32,
                            buffer,
                        ) == 0
                    {
                        to_delete = true;
                    }
                }
            }
            if !to_delete {
                if let Some(target) = w_max_target {
                    if we.v != target
                        && trans_path_check(
                            graph,
                            we.v,
                            target,
                            reverse_sources,
                            read_to_unitig,
                            max_ext as u32,
                            buffer,
                        ) == 0
                    {
                        to_delete = true;
                    }
                }
            }
            if v_max_target.is_some() && w_max_target.is_some() {
                assert!(ve.v != v_max_target.unwrap() || we.v != w_max_target.unwrap());
            }
        }

        if to_delete {
            graph.arcs[raw_index].deleted = true;
            graph.arcs[we_index].deleted = true;
            count_deleted += 1;
        }
    }

    if count_deleted > 0 {
        graph.cleanup();
    }
    count_deleted
}

/// `if_false_bub_links`: whether the candidate links really belong to a
/// false bubble.
///
/// The links `buffer[start_index..]` are deleted while the bubble is popped
/// and restored before returning; the workspace is the caller's `buf_t`
/// substitute.
#[allow(clippy::too_many_arguments)]
pub fn if_false_bub_links(
    graph: &mut Asg,
    v: u32,
    workspace: &mut BubbleBuffer,
    buffer: &mut Vec<u64>,
    start_index: usize,
    check_distance: u64,
) -> u32 {
    if graph.sequences[(v >> 1) as usize].deleted {
        return 0;
    }

    let mut result = 1u32;
    for i in start_index..buffer.len() {
        let index = buffer[i] as u32 as usize;
        let arc = graph.arcs[index];
        graph.arcs[index].deleted = true;
        graph.delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, true);
    }

    if graph.arc_count(v) < 2 || get_arcs(graph, v, None) < 2 {
        result = 0;
    }

    if result != 0 {
        result = 0;
        if asg_bub_pop1_primary_trio(
            graph,
            None,
            v,
            check_distance,
            workspace,
            u32::MAX,
            u32::MAX,
            false,
            None,
            &[],
            false,
            false,
        ) != 0 {
            for i in start_index..buffer.len() {
                let index = buffer[i] as u32 as usize;
                let arc = graph.arcs[index];
                graph.arcs[index].deleted = false;
                graph.delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, false);
            }

            let (start, count) = graph.arc_range(v);
            buffer.truncate(start_index);
            for i in 0..count {
                if graph.arcs[start + i].deleted {
                    continue;
                }
                let arc = graph.arcs[start + i];
                graph.arcs[start + i].deleted = true;
                graph.delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, true);
                buffer.push((start + i) as u64);
            }

            let end = workspace.sources[0];
            if asg_bub_pop1_primary_trio(
                graph,
                None,
                end ^ 1,
                check_distance,
                workspace,
                u32::MAX,
                u32::MAX,
                false,
                None,
                &[],
                false,
                false,
            ) != 0
            {
                result = 1;
            }
        }
    }

    for i in start_index..buffer.len() {
        let index = buffer[i] as u32 as usize;
        let arc = graph.arcs[index];
        graph.arcs[index].deleted = false;
        graph.delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, false);
    }

    result
}

/// `asg_arc_cut_bub_links`: cuts short links that a bubble search shows to be
/// supported by a false link structure.
#[allow(clippy::too_many_arguments)]
pub fn asg_arc_cut_bub_links(
    graph: &mut Asg,
    buffer: &mut Vec<u64>,
    length_ratio: f32,
    secondary_length_ratio: f32,
    occurrence_ratio: f32,
    is_occurrence: bool,
    check_distance: u64,
    reverse_sources: Option<&[PafList]>,
    read_to_unitig: Option<&ReadToUnitig>,
    max_ext: i32,
    workspace: &mut BubbleBuffer,
) -> u32 {
    buffer.clear();
    let vertex_count = graph.sequences.len() * 2;
    let use_reverse = reverse_sources.is_some() && read_to_unitig.is_some();

    for v in 0..vertex_count as u32 {
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if graph.sequence_visited[v as usize] != 0 {
            continue;
        }
        let (start, count) = graph.arc_range(v);
        if count < 2 {
            continue;
        }
        let mut live = 0u32;
        let mut total_overlap = 0u32;
        for arc in &graph.arcs[start..start + count] {
            if arc.deleted {
                continue;
            }
            live += 1;
            total_overlap = total_overlap.wrapping_add(arc.overlap);
        }
        if live < 2 {
            continue;
        }
        buffer.push((u64::from(total_overlap) << 32) | u64::from(v));
    }

    if use_reverse {
        graph.sequence_visited.fill(0);
    }
    radix_sort_u64(buffer);
    let candidate_count = buffer.len();

    let mut count = 0u32;
    for k in 0..candidate_count {
        let v = buffer[k] as u32;
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        let nv = graph.arc_count(v);
        if nv < 2 || get_arcs(graph, v, None) < 2 {
            continue;
        }

        buffer.truncate(candidate_count);
        let mut secondary_check = 0u32;
        let mut broke = false;
        let (start, count_arcs) = graph.arc_range(v);

        let mut i = 0usize;
        while i < count_arcs {
            let arc_index = start + i;
            let arc = graph.arcs[arc_index];
            i += 1;
            if arc.deleted {
                continue;
            }

            let w = arc.v ^ 1;
            let nw = graph.arc_count(w);
            if nw < 2 {
                broke = true;
                break;
            }
            let (w_start, w_count) = graph.arc_range(w);
            let mut kw = 0u32;
            let mut minimum_overlap = u32::MAX;
            let mut minimum_occurrence = u32::MAX;
            for t in 0..w_count {
                let neighbor = graph.arcs[w_start + t];
                if neighbor.deleted {
                    continue;
                }
                kw += 1;
                if neighbor.v == (v ^ 1) {
                    continue;
                }
                if neighbor.overlap < minimum_overlap {
                    minimum_overlap = neighbor.overlap;
                }
                if u32::from(neighbor.occurrence) < minimum_occurrence {
                    minimum_occurrence = u32::from(neighbor.occurrence);
                }
                buffer.push((w_start + t) as u64);
            }
            if kw < 2 {
                broke = true;
                break;
            }

            if (arc.overlap as f32) > (minimum_overlap as f32) * length_ratio
                && (arc.overlap as f32) > (minimum_overlap as f32) * secondary_length_ratio
            {
                broke = true;
                break;
            }
            if (arc.overlap as f32) > (minimum_overlap as f32) * length_ratio {
                secondary_check += 1;
            }
            if is_occurrence
                && (u32::from(arc.occurrence) as f32)
                    > (minimum_occurrence as f32) * occurrence_ratio
            {
                broke = true;
                break;
            }
        }

        if broke {
            continue;
        }

        if secondary_check != 0 {
            let mut reference: Option<u32> = None;
            let mut forward_broke = false;
            for i in 0..count_arcs {
                let arc = graph.arcs[start + i];
                if arc.deleted {
                    continue;
                }
                if reference.is_none() {
                    reference = Some(arc.v);
                } else if trans_path_check(
                    graph,
                    reference.unwrap(),
                    arc.v,
                    reverse_sources.expect("checked above"),
                    read_to_unitig.expect("checked above"),
                    max_ext as u32,
                    buffer,
                ) != 1
                {
                    forward_broke = true;
                    break;
                }
            }

            if forward_broke {
                if buffer.len() < candidate_count + 2 {
                    continue;
                }
                let mut reference: Option<u32> = None;
                let mut backward_broke = false;
                for i in candidate_count..buffer.len() {
                    let index = buffer[i] as u32 as usize;
                    let arc = graph.arcs[index];
                    if arc.deleted {
                        continue;
                    }
                    if get_arcs(graph, (arc.ul >> 32) as u32, None) != 2 {
                        backward_broke = true;
                        break;
                    }
                    if reference.is_none() {
                        reference = Some(arc.v);
                    } else if trans_path_check(
                        graph,
                        reference.unwrap(),
                        arc.v,
                        reverse_sources.expect("checked above"),
                        read_to_unitig.expect("checked above"),
                        max_ext as u32,
                        buffer,
                    ) != 1
                    {
                        backward_broke = true;
                        break;
                    }
                }
                if backward_broke {
                    continue;
                }
            }
        }

        if if_false_bub_links(graph, v, workspace, buffer, candidate_count, check_distance) != 0 {
            let (start, count_arcs) = graph.arc_range(v);
            for i in 0..count_arcs {
                if graph.arcs[start + i].deleted {
                    continue;
                }
                let arc = graph.arcs[start + i];
                graph.arcs[start + i].deleted = true;
                graph.delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, true);
            }
            count += 1;
        }
    }

    if count > 0 {
        graph.cleanup();
    }
    count
}

/// `asg_arc_cut_complex_bub_links`: cuts every link of a vertex that fails the
/// per-link length test, then restores those whose endpoints are not near a
/// bubble.
pub fn asg_arc_cut_complex_bub_links(
    graph: &mut Asg,
    buffer: &mut Vec<u64>,
    length_ratio: f32,
    occurrence_ratio: f32,
    is_occurrence: bool,
    labels: &mut BubbleLabel,
) -> u32 {
    buffer.clear();
    let vertex_count = graph.sequences.len() * 2;

    for v in 0..vertex_count as u32 {
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if graph.sequence_visited[v as usize] != 0 {
            continue;
        }
        let (start, count) = graph.arc_range(v);
        if count < 2 {
            continue;
        }
        let mut live = 0u32;
        let mut total_overlap = 0u32;
        for arc in &graph.arcs[start..start + count] {
            if arc.deleted {
                continue;
            }
            live += 1;
            total_overlap = total_overlap.wrapping_add(arc.overlap);
        }
        if live < 2 {
            continue;
        }
        buffer.push((u64::from(total_overlap) << 32) | u64::from(v));
    }

    radix_sort_u64(buffer);
    let candidate_count = buffer.len();

    for k in 0..candidate_count {
        let v = buffer[k] as u32;
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        let nv = graph.arc_count(v);
        if nv < 2 || get_arcs(graph, v, None) < 2 {
            continue;
        }

        let (start, count_arcs) = graph.arc_range(v);
        let mut broke = false;
        for i in 0..count_arcs {
            let arc = graph.arcs[start + i];
            if arc.deleted {
                continue;
            }
            let w = arc.v ^ 1;
            let (w_start, w_count) = graph.arc_range(w);
            if w_count < 2 {
                broke = true;
                break;
            }
            let mut kw = 0u32;
            let mut minimum_overlap = u32::MAX;
            let mut minimum_occurrence = u32::MAX;
            for t in 0..w_count {
                let neighbor = graph.arcs[w_start + t];
                if neighbor.deleted {
                    continue;
                }
                kw += 1;
                if neighbor.v == (v ^ 1) {
                    continue;
                }
                if neighbor.overlap < minimum_overlap {
                    minimum_overlap = neighbor.overlap;
                }
                if u32::from(neighbor.occurrence) < minimum_occurrence {
                    minimum_occurrence = u32::from(neighbor.occurrence);
                }
            }
            if kw < 2 {
                broke = true;
                break;
            }
            if (arc.overlap as f32) > (minimum_overlap as f32) * length_ratio {
                broke = true;
                break;
            }
            if is_occurrence
                && (u32::from(arc.occurrence) as f32)
                    > (minimum_occurrence as f32) * occurrence_ratio
            {
                broke = true;
                break;
            }
        }
        if broke {
            continue;
        }

        let (start, count_arcs) = graph.arc_range(v);
        for i in 0..count_arcs {
            if graph.arcs[start + i].deleted {
                continue;
            }
            let arc = graph.arcs[start + i];
            graph.arcs[start + i].deleted = true;
            graph.delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, true);
            buffer.push((start + i) as u64);
        }
    }

    let mut count = 0u32;
    if buffer.len() > candidate_count {
        asg_arc_identify_simple_bubbles_multi(graph, labels, false, 0);
        for k in candidate_count..buffer.len() {
            let index = buffer[k] as u32 as usize;
            let arc = graph.arcs[index];
            let v = (arc.ul >> 32) as u32;
            let w = arc.v;
            if graph.sequence_visited[v as usize] != 0
                || graph.sequence_visited[(v ^ 1) as usize] != 0
                || graph.sequence_visited[w as usize] != 0
                || graph.sequence_visited[(w ^ 1) as usize] != 0
            {
                count += 1;
                continue;
            }
            graph.arcs[index].deleted = false;
            graph.delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, false);
        }
    }

    if count > 0 {
        graph.cleanup();
    }
    count
}

/// `asg_cut_large_indel`: cuts the flank of an overlap carrying a large
/// insertion or deletion.
pub fn asg_cut_large_indel(
    graph: &mut Asg,
    buffer: &mut Vec<u64>,
    max_ext: i32,
    occurrence_ratio: f32,
    is_occurrence: bool,
    minimum_difference: u32,
) -> u32 {
    buffer.clear();
    let vertex_count = graph.sequences.len() * 2;

    for v in 0..vertex_count as u32 {
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if graph.sequence_visited[v as usize] != 0 {
            continue;
        }
        let (start, count) = graph.arc_range(v);
        if count < 2 {
            continue;
        }
        let mut live = 0u32;
        for arc in &graph.arcs[start..start + count] {
            if !arc.deleted {
                live += 1;
            }
        }
        if live < 2 {
            continue;
        }
        for i in 0..count {
            let arc = graph.arcs[start + i];
            if arc.deleted || arc.without_large_indel {
                continue;
            }
            buffer.push((u64::from(arc.overlap) << 32) | (start + i) as u64);
        }
    }

    radix_sort_u64(buffer);

    let mut count_deleted = 0u32;
    for k in 0..buffer.len() {
        let raw_index = buffer[k] as u32 as usize;
        if graph.arcs[raw_index].deleted {
            continue;
        }

        let ve = graph.arcs[raw_index];
        let v = (ve.ul >> 32) as u32;
        let w = ve.v ^ 1;
        if graph.sequences[(v >> 1) as usize].deleted
            || graph.sequences[(w >> 1) as usize].deleted
        {
            continue;
        }
        let nv = graph.arc_count(v);
        let nw = graph.arc_count(w);
        if nv <= 1 && nw <= 1 {
            continue;
        }

        let (start, count) = graph.arc_range(v);
        let mut kv = 0u32;
        let mut occurrence_max = 0u32;
        let mut overlap_max = 0u32;
        let mut v_edge: Option<usize> = None;
        for i in 0..count {
            let arc = graph.arcs[start + i];
            if arc.deleted {
                continue;
            }
            if arc.v == (w ^ 1) {
                v_edge = Some(start + i);
            }
            kv += 1;
            if occurrence_max < u32::from(arc.occurrence) {
                occurrence_max = u32::from(arc.occurrence);
            }
            if overlap_max < arc.overlap {
                overlap_max = arc.overlap;
            }
        }
        if kv < 1 {
            continue;
        }
        if kv >= 2 {
            let v_edge = graph.arcs[v_edge.expect("the candidate arc is live")];
            if is_occurrence
                && (u32::from(v_edge.occurrence) as f32)
                    > (occurrence_max as f32) * occurrence_ratio
            {
                continue;
            }
            if v_edge.overlap.wrapping_add(minimum_difference) > overlap_max {
                continue;
            }
        }

        let (start, count) = graph.arc_range(w);
        let mut kw = 0u32;
        let mut occurrence_max = 0u32;
        let mut overlap_max = 0u32;
        let mut w_edge: Option<usize> = None;
        for i in 0..count {
            let arc = graph.arcs[start + i];
            if arc.deleted {
                continue;
            }
            if arc.v == (v ^ 1) {
                w_edge = Some(start + i);
            }
            kw += 1;
            if occurrence_max < u32::from(arc.occurrence) {
                occurrence_max = u32::from(arc.occurrence);
            }
            if overlap_max < arc.overlap {
                overlap_max = arc.overlap;
            }
        }
        if kw < 1 {
            continue;
        }
        if kw >= 2 {
            let w_edge = graph.arcs[w_edge.expect("the mirror arc is live")];
            if is_occurrence
                && (u32::from(w_edge.occurrence) as f32)
                    > (occurrence_max as f32) * occurrence_ratio
            {
                continue;
            }
            if w_edge.overlap.wrapping_add(minimum_difference) > overlap_max {
                continue;
            }
        }

        if kv <= 1 && kw <= 1 {
            continue;
        }

        let mut to_delete = false;
        if kv > 1 && kw > 1 {
            to_delete = true;
        } else if kw == 1 {
            if asg_topocut_aux(graph, w ^ 1, max_ext) < max_ext {
                to_delete = true;
            }
        } else if kv == 1 && asg_topocut_aux(graph, v ^ 1, max_ext) < max_ext {
            to_delete = true;
        }

        if to_delete {
            let v_edge = v_edge.expect("the candidate arc is live");
            let w_edge = w_edge.expect("the mirror arc is live");
            graph.arcs[v_edge].deleted = true;
            graph.arcs[w_edge].deleted = true;
            count_deleted += 1;
        }
    }

    if count_deleted > 0 {
        graph.cleanup();
    }
    count_deleted
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::asg::AsgArc;

    fn chain_graph(node_count: usize) -> Asg {
        // A linear chain 0 -> 2 -> 4 -> ... plus its reverse-complement
        // mirror, so every merge is symmetric.
        let mut graph = Asg::new();
        for read in 0..node_count {
            graph.set_sequence(read, 100, false);
        }
        for read in 0..node_count - 1 {
            let from = (read as u32) << 1;
            let to = ((read as u32) + 1) << 1;
            graph.push_arc(AsgArc {
                ul: (u64::from(from) << 32) | 100,
                v: to,
                ..AsgArc::default()
            });
            graph.push_arc(AsgArc {
                ul: (u64::from(to + 1) << 32) | 100,
                v: from + 1,
                ..AsgArc::default()
            });
        }
        graph.cleanup();
        graph
    }

    #[test]
    fn real_length_counts_only_live_arcs() {
        let mut graph = chain_graph(2);
        assert_eq!(get_real_length(&graph, 0), 1);
        assert_eq!(get_real_length_one(&graph, 0), Some(2));
        assert_eq!(get_real_length(&graph, 2), 0);
        graph.arcs_of_mut(0)[0].deleted = true;
        assert_eq!(get_real_length(&graph, 0), 0);
        assert_eq!(get_real_length_one(&graph, 0), None);
    }

    #[test]
    fn tips_shorter_than_the_limit_are_cut() {
        // Three reads in a line: two merges, both below max_ext = 3.
        let mut graph = chain_graph(3);
        let mut buffer = Vec::new();
        let count = asg_arc_cut_tips(&mut graph, 3, &mut buffer, false, None, None);
        assert_eq!(count, 1);
        assert!(graph.sequences.iter().all(|sequence| sequence.deleted));
        assert!(graph.arcs.is_empty());
    }

    #[test]
    fn chains_as_long_as_the_limit_survive() {
        let mut graph = chain_graph(3);
        let mut buffer = Vec::new();
        let count = asg_arc_cut_tips(&mut graph, 2, &mut buffer, false, None, None);
        assert_eq!(count, 0);
        assert!(graph.sequences.iter().all(|sequence| !sequence.deleted));
        // Two merges fill the limit, so no candidate is even recorded.
        assert!(buffer.is_empty());
    }

    #[test]
    fn a_chain_carrying_a_telomere_is_preserved() {
        let mut graph = chain_graph(3);
        let telo_end = TeloEnd {
            hh: vec![0, 0, 1],
        };
        let mut buffer = Vec::new();
        let count = asg_arc_cut_tips(
            &mut graph,
            3,
            &mut buffer,
            false,
            None,
            Some(&telo_end),
        );
        // The walk reaches read 2, whose telomere flag suppresses the cut of
        // the whole chain.
        assert_eq!(count, 0);
        assert!(graph.sequences.iter().all(|sequence| !sequence.deleted));
    }

    fn overlap(
        query_start: u32,
        query_end: u32,
        exact: bool,
    ) -> crate::overlap::paf::PafHit {
        crate::overlap::paf::PafHit {
            query_identifier: 0,
            query_start,
            query_end,
            target_identifier: 1,
            target_start: 0,
            target_end: query_end - query_start,
            matching_length: 90,
            block_length: query_end - query_start,
            reverse: false,
            without_large_indel: true,
            exact_region: exact,
            deleted: false,
        }
    }

    fn overlap_list(hits: Vec<crate::overlap::paf::PafHit>) -> PafList {
        PafList {
            hits,
            is_fully_corrected: true,
            is_abnormal: 0,
        }
    }

    #[test]
    fn unsupported_middles_are_chimeric() {
        let mut buffer = Vec::new();

        // The two anchored sides do not reach each other.
        let source = overlap_list(vec![overlap(0, 40, true), overlap(60, 100, true)]);
        assert_eq!(if_sup_chimeric(&source, 100, &mut buffer, true), 1);
        assert!(buffer.is_empty());

        // A contained exact overlap joins them into one interval.
        let source = overlap_list(vec![
            overlap(0, 40, true),
            overlap(60, 100, true),
            overlap(20, 80, true),
        ]);
        assert_eq!(if_sup_chimeric(&source, 100, &mut buffer, true), 0);

        // Without exact hits both sides stay empty, which counts as chimeric.
        let source = overlap_list(vec![overlap(0, 40, false), overlap(60, 100, false)]);
        assert_eq!(if_sup_chimeric(&source, 100, &mut buffer, true), 1);
    }

    #[test]
    fn a_spanning_exact_overlap_is_support() {
        let source = overlap_list(vec![
            overlap(0, 40, true),
            overlap(60, 100, true),
            overlap(0, 100, true),
        ]);
        let mut buffer = Vec::new();
        assert_eq!(if_sup_chimeric(&source, 100, &mut buffer, true), 0);
    }

    #[test]
    fn follow_limit_path_stops_at_a_dead_end() {
        let graph = chain_graph(3);
        let mut end = 0u32;
        let mut occurrence = 0u32;
        let kind = follow_limit_path(&graph, 0, &mut end, &mut occurrence, None, 100);
        assert_eq!(kind, END_TIPS);
        assert_eq!(end, 4);
        assert_eq!(occurrence, 3);
    }

    #[test]
    fn topocut_walks_until_the_chain_ends() {
        let graph = chain_graph(3);
        assert_eq!(asg_topocut_aux(&graph, 0, 10), 0);
        assert_eq!(asg_topocut_aux(&graph, 1, 10), 2);
    }

    fn asymmetric_graph(reads: usize, arcs: &[(u32, u32)]) -> Asg {
        let mut graph = Asg::new();
        for read in 0..reads {
            graph.set_sequence(read, 100, false);
        }
        for &(from, to) in arcs {
            graph.push_arc(AsgArc {
                ul: (u64::from(from) << 32) | 100,
                v: to,
                overlap: 90,
                ..AsgArc::default()
            });
        }
        graph.cleanup();
        graph
    }

    #[test]
    fn semi_circles_lose_their_closing_input() {
        // Vertex 0 has three inputs (arcs of vertex 1) and one output. Its
        // path runs 0 -> 2 -> 0 and stops because vertex 1 has three out-arcs,
        // so the closing edge 1 -> 3 and its mirror 2 -> 0 are cut.
        let mut graph = asymmetric_graph(
            3,
            &[(0, 2), (1, 0), (1, 4), (3, 0), (1, 3), (2, 0)],
        );
        let count = asg_cut_semi_circ(&mut graph, LIM_LEN, true);
        assert_eq!(count, 1);
        assert!(graph.get_arc(1, 3).is_none());
        assert!(graph.get_arc(2, 0).is_none());
        assert!(graph.get_arc(1, 0).is_some());
    }

    fn split_graph(second_without_large_indel: bool) -> Asg {
        // A symmetric split from vertex 0: a short branch through read 1 and
        // a long branch through read 2, both mirrored from the end read 3.
        let mut graph = Asg::new();
        for read in 0..4 {
            graph.set_sequence(read, 100, false);
        }
        graph.sequence_visited = vec![0u8; 8];
        for &(from, to, overlap, without_large_indel) in &[
            (0u32, 2u32, 10u32, false),
            (0, 4, 100, second_without_large_indel),
            (3, 1, 10, false),
            (5, 1, 100, false),
        ] {
            graph.push_arc(AsgArc {
                ul: (u64::from(from) << 32) | u64::from(overlap),
                v: to,
                overlap,
                without_large_indel,
                ..AsgArc::default()
            });
        }
        graph.cleanup();
        graph
    }

    #[test]
    fn short_branches_are_cut_when_the_long_side_survives() {
        let mut graph = split_graph(false);
        let mut buffer = Vec::new();
        let count = asg_arc_cut_length(
            &mut graph,
            &mut buffer,
            3,
            0.2,
            0.6,
            false,
            false,
            true,
            0,
            1,
            None,
            None,
            None,
            &[],
        );
        assert_eq!(count, 1);
        assert!(graph.get_arc(0, 2).is_none());
        assert!(graph.get_arc(3, 1).is_none());
        assert!(graph.get_arc(0, 4).is_some());
    }

    #[test]
    fn large_indels_cut_the_supported_flank() {
        let mut graph = split_graph(true);
        let mut buffer = Vec::new();
        let count = asg_cut_large_indel(&mut graph, &mut buffer, 3, 0.75, false, 0);
        assert_eq!(count, 1);
        assert!(graph.get_arc(0, 2).is_none());
        assert!(graph.get_arc(3, 1).is_none());
        assert!(graph.get_arc(0, 4).is_some());
    }
}

/// `gen_telo_end_t` (`Overlaps.cpp` 300) with its `mark_telo_ends` worker
/// (`Overlaps.cpp` 191): scans both ends of every read for the telomere
/// motif or any of its rotations.
///
/// Each end runs two passes: the first accumulates `+1` per motif hit and
/// `-penalty` per non-hit (aborting early once the best score is reached),
/// and a read end is telomeric only when the best score reaches
/// `minimum_score`. The second pass repeats the scan with the rotation
/// continuity refinement of the reference. Single-threaded: the reads are
/// independent, so `-t1` parity holds for any thread count.
pub fn gen_telo_ends(
    reads: &ReadStore,
    motif: &str,
    minimum_score: i64,
    penalty: i64,
    maximum_drop: i64,
) -> TeloEnd {
    let motif_length = motif.len();
    let mask = (1u64 << (motif_length * 2)) - 1;
    let mut motif_code = 0u64;
    for byte in motif.bytes() {
        let code = base_code(byte).expect("the telomere motif must be A/C/G/T only");
        motif_code = (motif_code << 2) | code;
    }

    // The rotation set of the motif (the reference's `aux->h` hash).
    let mut rotations = std::collections::HashSet::new();
    let mut code = motif_code;
    for _ in 0..motif_length {
        rotations.insert(code);
        code = (((code >> ((motif_length - 1) << 1)) & 3) | (code << 2)) & mask;
    }
    let penalty = penalty.abs();

    let mut hh = vec![0u8; reads.total_reads as usize];
    for rid in 0..reads.total_reads as usize {
        // 5' end: forward scan; 3' end: the reverse complement scanned
        // forward, which is exactly the reference's backward `RC_CHAR` walk.
        let forward = reads.read_sequence(rid);
        if end_is_telomeric(
            &forward,
            motif_length,
            mask,
            &rotations,
            minimum_score,
            penalty,
            maximum_drop,
        ) && refine_end(
            &forward,
            motif_length,
            mask,
            motif_code,
            &rotations,
            minimum_score,
            penalty,
            maximum_drop,
        ) {
            hh[rid] |= 1;
        }

        let reverse = reads.read_sequence_reverse_complement(rid);
        if end_is_telomeric(
            &reverse,
            motif_length,
            mask,
            &rotations,
            minimum_score,
            penalty,
            maximum_drop,
        ) && refine_end(
            &reverse,
            motif_length,
            mask,
            motif_code,
            &rotations,
            minimum_score,
            penalty,
            maximum_drop,
        ) {
            hh[rid] |= 2;
        }
    }

    let t5 = hh.iter().filter(|flag| *flag & 1 != 0).count();
    let t3 = hh.iter().filter(|flag| *flag & 2 != 0).count();
    eprintln!(
        "[M::gen_telo_end_t::] ==> # 5'-telomeres::{}, # 3'-telomeres::{}, # tot::{}, motif::{}, motif_len::{}",
        t5, t3, reads.total_reads, motif, motif_length
    );

    TeloEnd { hh }
}

/// `seq_nt6_table` restricted to A/C/G/T; `None` for every other base.
fn base_code(byte: u8) -> Option<u64> {
    match byte {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// The first (detection) pass of `mark_telo_ends`: whether the end's best
/// score reaches `minimum_score`.
fn end_is_telomeric(
    sequence: &[u8],
    motif_length: usize,
    mask: u64,
    rotations: &std::collections::HashSet<u64>,
    minimum_score: i64,
    penalty: i64,
    maximum_drop: i64,
) -> bool {
    let mut code = 0u64;
    let mut consecutive = 0i64;
    let mut score = 0i64;
    let mut max_score = 0i64;

    for (k, byte) in sequence.iter().enumerate() {
        let mut hit = false;
        match base_code(*byte) {
            Some(c) => {
                code = ((code << 2) | c) & mask;
                consecutive += 1;
                if consecutive as usize >= motif_length && rotations.contains(&code) {
                    hit = true;
                }
            }
            None => {
                consecutive = 0;
                code = 0;
            }
        }
        if k >= motif_length {
            score += if hit { 1 } else { -penalty };
        }
        if score > max_score {
            max_score = score;
        } else if max_score - score > maximum_drop {
            break;
        }

        if max_score >= minimum_score {
            break;
        }
    }

    max_score >= minimum_score
}

/// The second (refinement) pass of `mark_telo_ends`, with the rotation
/// continuity tracking; returns whether the end still reaches
/// `minimum_score`.
#[allow(clippy::too_many_arguments)]
fn refine_end(
    sequence: &[u8],
    motif_length: usize,
    mask: u64,
    motif_code: u64,
    rotations: &std::collections::HashSet<u64>,
    minimum_score: i64,
    penalty: i64,
    maximum_drop: i64,
) -> bool {
    let mut code = 0u64;
    let mut consecutive = 0i64;
    let mut score = 0i64;
    let mut max_score = 0i64;
    let mut expected = u64::MAX;

    for (k, byte) in sequence.iter().enumerate() {
        let mut hit = false;
        match base_code(*byte) {
            Some(c) => {
                code = ((code << 2) | c) & mask;
                consecutive += 1;
                if consecutive as usize >= motif_length && rotations.contains(&code) {
                    if expected != u64::MAX && expected == code {
                        hit = true;
                    } else {
                        let mut z = 0usize;
                        let mut walk = motif_code;
                        while z < motif_length && walk != code {
                            walk = (((walk >> ((motif_length - 1) << 1)) & 3) | (walk << 2)) & mask;
                            z += 1;
                        }
                        if z < motif_length && walk == code {
                            hit = true;
                        }
                    }
                }
            }
            None => {
                consecutive = 0;
                code = 0;
            }
        }

        if !hit {
            expected = u64::MAX;
        } else {
            expected = (((code >> ((motif_length - 1) << 1)) & 3) | (code << 2)) & mask;
        }
        if k >= motif_length {
            score += if hit { 1 } else { -penalty };
        }
        if score > max_score {
            max_score = score;
        } else if max_score - score > maximum_drop {
            break;
        }

        if max_score >= minimum_score {
            break;
        }
    }

    max_score >= minimum_score
}

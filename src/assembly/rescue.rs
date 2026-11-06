//! The read-recovery passes run at the end of `ul_clean_gfa`, the counterpart
//! of `pre_clean` and the `rescue_*_aggressive` family in `Overlaps.cpp`.
//!
//! Every function here works on the read graph (`asg_t`) and on a temporary
//! unitig graph built by `ma_ug_gen`, exactly like the reference. The
//! reference threads raw pointers through the whole family; the Rust version
//! passes the owning slices instead, and keeps the reference's iteration and
//! deletion order so the outputs stay byte-identical.
//!
//! The reference has several `i_ug != NULL` entry points: the caller can hand
//! in an existing unitig graph instead of having one built. The default HiFi
//! pipeline always passes `NULL`, but the primary-contig stage uses the
//! explicit form, so the parameter is kept.

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::bubble::BubbleLabel;
use crate::assembly::bubble_identify::asg_arc_identify_simple_bubbles_multi;
use crate::assembly::gfa_clean::{
    get_real_length, get_real_length_one, END_TIPS, LONG_TIPS, LOOP, MUL_INPUT, MUL_OUTPUT,
    TWO_INPUT, TWO_OUTPUT,
};
use crate::assembly::hit::{
    ma_hit_to_arc, CoverageCut, ReadToUnitig, ALTER_LABEL, PRIMARY_LABEL,
};
use crate::assembly::unitig::{ma_ug_gen, UnitigGraph, FATHER, MOTHER};
use crate::correction::ec_driver::PafList;
use crate::overlap::paf::PafHit;

/// The overlap-conversion constants shared by the rescue family.
#[derive(Debug, Clone, Copy)]
pub struct RescueParameters {
    /// `ug_opt_t::max_hang`.
    pub maximum_hang: i32,
    /// `asm_opt.max_hang_rate`.
    pub maximum_hang_rate: f32,
    /// `ug_opt_t::min_ovlp`.
    pub minimum_overlap: i32,
    /// The reference's global `min_thres`, used by the bubble labeler.
    pub minimum_length_threshold: i64,
}

/// The three output parameters of `get_contained_reads_chain`.
#[derive(Debug, Clone, Copy)]
struct ChainOutcome {
    /// The function's return value.
    found: bool,
    /// `return_ava_cur`; `u32::MAX` when the chain did not link to a unitig.
    available_count: u32,
    /// `return_ava_ol`; `u32::MAX` when the chain did not link to a unitig.
    available_overlap: u32,
    /// `return_chainLen`; `u32::MAX` when the chain did not link to a unitig.
    chain_length: u32,
}

/// `non_trio_flag`: the label of the opposite haplotype, or `-1`.
fn non_trio_flag(trio_flag: u8, ignore_trio_flag: u32) -> u32 {
    if ignore_trio_flag != 0 {
        return u32::MAX;
    }
    if trio_flag == FATHER {
        u32::from(MOTHER)
    } else if trio_flag == MOTHER {
        u32::from(FATHER)
    } else {
        u32::MAX
    }
}

/// `ug->u.a[vertex>>1].len`: the length of the unitig that owns `vertex`.
fn unitig_length(unitigs: &[MaUtg], vertex: u32) -> u32 {
    unitigs[(vertex >> 1) as usize].length
}

/// `detect_single_path_with_dels_by_length`: walks a chain until it ends, forks
/// or exceeds `maximum_length`.
///
/// `end` receives the last vertex visited and `length` the number of nodes.
/// When `buffer` is given, every visited vertex is appended to it.
pub fn detect_single_path_with_dels_by_length(
    graph: &Asg,
    beginning: u32,
    end: &mut u32,
    length: &mut i64,
    mut buffer: Option<&mut Vec<u32>>,
    maximum_length: i64,
) -> u8 {
    let mut v = beginning;
    *length = 0;
    loop {
        *length += 1;
        let kv = get_real_length(graph, v);
        *end = v;
        if let Some(visited) = buffer.as_deref_mut() {
            visited.push(v);
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
        if *length > maximum_length {
            return LONG_TIPS;
        }

        // `kv` is one here, so `get_real_length` writes the only live target.
        let w = get_real_length_one(graph, v).expect("kv == 1 implies one live arc");
        let kw = get_real_length(graph, w ^ 1);
        v = w;
        *end = v;

        if kw == 2 {
            *length += 1;
            if let Some(visited) = buffer.as_deref_mut() {
                visited.push(v);
            }
            return TWO_INPUT;
        }
        if kw > 2 {
            *length += 1;
            if let Some(visited) = buffer.as_deref_mut() {
                visited.push(v);
            }
            return MUL_INPUT;
        }
        if (v >> 1) == (beginning >> 1) {
            return LOOP;
        }
    }
}

/// `asg_arc_del_simple_circle_untig`: removes the short simple circles.
///
/// A vertex whose forward side has exactly one live arc is walked through
/// `detect_single_path_with_dels_by_length`; when the walk ends at a fork
/// before `circle_length` nodes, the reverse arc closing the circle is
/// deleted. The reference's `sources`/`coverage_cut` arguments are unused on
/// this path and are omitted here.
pub fn asg_arc_del_simple_circle_untig(
    graph: &mut Asg,
    circle_length: i64,
    is_drop: bool,
) -> u32 {
    let vertex_count = graph.sequences.len() as u32 * 2;
    let mut reduced = 0u32;
    let mut buffer: Vec<u32> = Vec::new();

    for v in 0..vertex_count {
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if is_drop && graph.sequences[(v >> 1) as usize].coverage_class == ALTER_LABEL {
            continue;
        }
        if get_real_length(graph, v ^ 1) <= 1 {
            continue;
        }
        if get_real_length(graph, v) != 1 {
            continue;
        }

        // The reference captures the span once and only rewrites `del` flags
        // inside the loop, so the span stays valid.
        let (arc_start, arc_count) = graph.arc_range(v);
        for i in 0..arc_count as usize {
            if graph.arcs[arc_start + i].deleted {
                continue;
            }

            buffer.clear();
            let mut convex = 0u32;
            let mut length = 0i64;
            let flag = detect_single_path_with_dels_by_length(
                graph,
                v,
                &mut convex,
                &mut length,
                Some(&mut buffer),
                circle_length,
            );
            if length > circle_length || flag == LONG_TIPS {
                break;
            }
            if flag == LOOP {
                break;
            }
            if flag != END_TIPS && flag != LONG_TIPS {
                if v == convex && buffer.len() > 1 {
                    convex = buffer[buffer.len() - 2];
                    buffer.pop();
                }

                let w = v ^ 1;
                if get_real_length(graph, w) == 0 {
                    break;
                }
                let (reverse_start, reverse_count) = graph.arc_range(w);
                for k in 0..reverse_count as usize {
                    let arc = graph.arcs[reverse_start + k];
                    if !arc.deleted && arc.v == (convex ^ 1) {
                        graph.arcs[reverse_start + k].deleted = true;
                        graph.delete_arc(arc.v ^ 1, ((arc.ul >> 32) as u32) ^ 1, true);
                        reduced += 1;
                    }
                }
            }
        }
    }

    if reduced != 0 {
        graph.cleanup();
        graph.symmetrize();
    }

    reduced
}

/// `pre_clean`: alternates the simple-circle removal until it stops reducing.
///
/// The reference's `pop_s_node` branch (`asg_arc_del_single_node_directly`) is
/// only used by commented-out call sites, so every live caller passes zero and
/// the branch is not ported.
pub fn pre_clean(graph: &mut Asg) {
    loop {
        let reduced = asg_arc_del_simple_circle_untig(graph, 100, false);
        if reduced == 0 {
            break;
        }
    }
}

/// `get_corresponding_uId`: the unitig vertex matching read vertex `t`.
fn get_corresponding_u_id(
    unitigs: &[MaUtg],
    unitig_graph: &Asg,
    read_to_unitig: &ReadToUnitig,
    t: u32,
) -> u32 {
    let (u_id, is_unitig) = read_to_unitig.get(t >> 1);
    if u_id == u32::MAX || is_unitig == 0 {
        return u32::MAX;
    }
    if unitig_graph.sequences[u_id as usize].deleted {
        return u32::MAX;
    }
    let unitig = &unitigs[u_id as usize];
    let mut w = u32::MAX;
    if t == unitig.start {
        w = u_id << 1;
    }
    if t == unitig.end {
        w = (u_id << 1) ^ 1;
    }
    w
}

/// `get_num_edges2existing_nodes_advance`: counts the overlaps of `query` that
/// reach an existing unitig and returns their total overlap length.
#[allow(clippy::too_many_arguments)]
fn get_num_edges2existing_nodes_advance(
    unitigs: &[MaUtg],
    unitig_graph: &Asg,
    read_graph: &Asg,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    read_to_unitig: &ReadToUnitig,
    trio_labels: &[u8],
    parameters: RescueParameters,
    query: u32,
    skip_unitigs: &[u32],
    ignore_trio_flag: u32,
) -> (u32, u32) {
    let query_read = (query >> 1) as usize;
    let hits = &sources[query_read].hits;
    let opposite = non_trio_flag(trio_labels[query_read], ignore_trio_flag);
    let mut occurrences = 0u32;
    let mut overlap_length = 0u32;

    for h in hits {
        let query_cut = &coverage_cut[h.query_identifier as usize];
        let target_cut = &coverage_cut[h.target_identifier as usize];
        // The query side may already be deleted; the target may not.
        if target_cut.deleted || read_graph.sequences[h.target_identifier as usize].deleted {
            continue;
        }
        if u32::from(trio_labels[h.target_identifier as usize]) == opposite {
            continue;
        }

        let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
        let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
        let (result, arc) = ma_hit_to_arc(
            h,
            query_length,
            target_length,
            parameters.maximum_hang,
            parameters.maximum_hang_rate,
            parameters.minimum_overlap,
        );
        // A contained read cannot extend the unitig.
        if result < 0 {
            continue;
        }
        if arc.query_vertex() != query {
            continue;
        }
        let (u_id, is_unitig) = read_to_unitig.get(arc.v >> 1);
        if u_id == u32::MAX || is_unitig == 0 {
            continue;
        }
        if skip_unitigs.contains(&u_id) {
            continue;
        }
        if unitig_graph.sequences[u_id as usize].deleted {
            continue;
        }

        let unitig = &unitigs[u_id as usize];
        if arc.v != unitig.start && arc.v != unitig.end {
            continue;
        }
        // The reference runs both tests in sequence; the second one wins when
        // start and end coincide.
        let mut v = 0u32;
        if arc.v == unitig.start {
            v = (u_id << 1) ^ 1;
        }
        if arc.v == unitig.end {
            v = u_id << 1;
        }

        if get_real_length(unitig_graph, v) == 1 {
            let w = get_real_length_one(unitig_graph, v).expect("one live arc");
            if get_real_length(unitig_graph, w ^ 1) == 1 {
                continue;
            }
        }

        occurrences += 1;
        overlap_length = overlap_length.wrapping_add(arc.overlap);
    }

    (occurrences, overlap_length)
}

/// `get_edge2existing_node_advance`: the first remaining hit of `query` that
/// reaches an existing unitig.
///
/// `index` is the scan position shared by successive calls; it advances past
/// every examined hit, including the matching one.
#[allow(clippy::too_many_arguments)]
fn get_edge2existing_node_advance(
    unitigs: &[MaUtg],
    unitig_graph: &Asg,
    read_graph: &Asg,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    read_to_unitig: &ReadToUnitig,
    trio_labels: &[u8],
    parameters: RescueParameters,
    query: u32,
    skip_unitigs: &[u32],
    index: &mut usize,
    t: &mut AsgArc,
    ignore_trio_flag: u32,
) -> bool {
    let query_read = (query >> 1) as usize;
    let hits = &sources[query_read].hits;
    let opposite = non_trio_flag(trio_labels[query_read], ignore_trio_flag);

    while *index < hits.len() {
        let h = &hits[*index];
        *index += 1;

        let query_cut = &coverage_cut[h.query_identifier as usize];
        let target_cut = &coverage_cut[h.target_identifier as usize];
        if target_cut.deleted || read_graph.sequences[h.target_identifier as usize].deleted {
            continue;
        }
        if u32::from(trio_labels[h.target_identifier as usize]) == opposite {
            continue;
        }

        let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
        let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
        let (result, arc) = ma_hit_to_arc(
            h,
            query_length,
            target_length,
            parameters.maximum_hang,
            parameters.maximum_hang_rate,
            parameters.minimum_overlap,
        );
        // The reference writes `t` on every successful conversion and tests the
        // filters against it; a rejection leaves the previous value in place.
        if result < 0 {
            continue;
        }
        *t = arc;

        if t.query_vertex() != query {
            continue;
        }
        let (u_id, is_unitig) = read_to_unitig.get(t.v >> 1);
        if u_id == u32::MAX || is_unitig == 0 {
            continue;
        }
        if skip_unitigs.contains(&u_id) {
            continue;
        }
        if unitig_graph.sequences[u_id as usize].deleted {
            continue;
        }

        let unitig = &unitigs[u_id as usize];
        if t.v != unitig.start && t.v != unitig.end {
            continue;
        }
        let mut v = 0u32;
        if t.v == unitig.start {
            v = (u_id << 1) ^ 1;
        }
        if t.v == unitig.end {
            v = u_id << 1;
        }

        if get_real_length(unitig_graph, v) == 1 {
            let w = get_real_length_one(unitig_graph, v).expect("one live arc");
            if get_real_length(unitig_graph, w ^ 1) == 1 {
                continue;
            }
        }

        return true;
    }

    false
}

/// `get_edge_from_source`: the stored overlap of `query` that becomes the
/// `query -> target` arc.
///
/// The reference overwrites `t` on every successful `ma_hit2arc` call, even
/// when the arc does not match, so the caller's stale value semantics are
/// preserved here.
fn get_edge_from_source_into(
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    parameters: RescueParameters,
    query: u32,
    target: u32,
    t: &mut AsgArc,
) -> bool {
    let query_read = (query >> 1) as usize;
    for h in &sources[query_read].hits {
        let query_cut = &coverage_cut[h.query_identifier as usize];
        let target_cut = &coverage_cut[h.target_identifier as usize];
        let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
        let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
        let (result, arc) = ma_hit_to_arc(
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
        *t = arc;
        if t.query_vertex() != query {
            continue;
        }
        if t.v != target {
            continue;
        }
        return true;
    }
    false
}

/// `get_best_contained_read`: the deleted overlap whose target is a contained
/// read with the longest overlap to `query`.
#[allow(clippy::too_many_arguments)]
fn get_best_contained_read<'a>(
    unitig_graph: &Asg,
    read_graph: &Asg,
    sources: &'a [PafList],
    coverage_cut: &[CoverageCut],
    read_to_unitig: &ReadToUnitig,
    trio_labels: &[u8],
    parameters: RescueParameters,
    query: u32,
    ignore_trio_flag: u32,
) -> Option<&'a PafHit> {
    let query_read = (query >> 1) as usize;
    let hits = &sources[query_read].hits;
    let opposite = non_trio_flag(trio_labels[query_read], ignore_trio_flag);
    let mut best: Option<&PafHit> = None;
    let mut best_overlap = 0u32;

    for h in hits {
        let query_cut = &coverage_cut[h.query_identifier as usize];
        let target_cut = &coverage_cut[h.target_identifier as usize];
        // Only deleted edges may describe contained reads.
        if !target_cut.deleted {
            continue;
        }
        if !h.deleted {
            continue;
        }
        if u32::from(trio_labels[h.target_identifier as usize]) == opposite {
            continue;
        }

        // The target must be contained in another live read.
        let (contained_read, is_unitig) = read_to_unitig.get(h.target_identifier);
        if contained_read == u32::MAX || is_unitig == 1 {
            continue;
        }
        if read_graph.sequences[contained_read as usize].deleted {
            continue;
        }

        let (contained_unitig, is_unitig) = read_to_unitig.get(contained_read);
        if contained_unitig == u32::MAX || is_unitig != 1 {
            continue;
        }
        if unitig_graph.sequences[contained_unitig as usize].deleted {
            continue;
        }

        let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
        let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
        let (result, arc) = ma_hit_to_arc(
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
        if arc.query_vertex() != query {
            continue;
        }

        if arc.overlap > best_overlap {
            best_overlap = arc.overlap;
            best = Some(h);
        }
    }

    best
}

/// `get_contained_reads_chain`: follows the contained reads starting at `h`
/// until one of them reaches an existing unitig.
#[allow(clippy::too_many_arguments)]
fn get_contained_reads_chain<'a>(
    mut h: Option<&'a PafHit>,
    sources: &'a [PafList],
    coverage_cut: &[CoverageCut],
    read_to_unitig: &ReadToUnitig,
    unitigs: &[MaUtg],
    unitig_graph: &Asg,
    read_graph: &Asg,
    trio_labels: &[u8],
    parameters: RescueParameters,
    mut end_read: u32,
    unitig_vertex: u32,
    chain_buffer: &mut Vec<u32>,
    mut chain_edges: Option<&mut Vec<AsgArc>>,
    threshold_length: u32,
    ignore_trio_flag: u32,
) -> ChainOutcome {
    chain_buffer.clear();
    chain_buffer.push(unitig_vertex);
    if let Some(edges) = chain_edges.as_deref_mut() {
        edges.clear();
    }

    let mut chain_length = 0u32;
    while let Some(hit) = h {
        let query_cut = &coverage_cut[hit.query_identifier as usize];
        let target_cut = &coverage_cut[hit.target_identifier as usize];
        let trio_flag = trio_labels[hit.query_identifier as usize];

        // Do not link reads of different haplotypes.
        if ignore_trio_flag == 0 {
            let opposite = non_trio_flag(trio_flag, ignore_trio_flag);
            if u32::from(trio_labels[hit.target_identifier as usize]) == opposite {
                break;
            }
        }

        // Only deleted edges may describe contained reads.
        if !target_cut.deleted {
            break;
        }
        if !hit.deleted {
            break;
        }

        let (contained_read, is_unitig) = read_to_unitig.get(hit.target_identifier);
        if contained_read == u32::MAX || is_unitig == 1 {
            break;
        }
        if read_graph.sequences[contained_read as usize].deleted {
            break;
        }

        let (contained_unitig, is_unitig) = read_to_unitig.get(contained_read);
        if contained_unitig == u32::MAX || is_unitig != 1 {
            break;
        }
        if unitig_graph.sequences[contained_unitig as usize].deleted {
            break;
        }

        let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
        let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
        let (result, t) = ma_hit_to_arc(
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
        if t.query_vertex() != end_read {
            break;
        }

        if let Some(edges) = chain_edges.as_deref_mut() {
            edges.push(t);
        }
        chain_length += 1;
        chain_buffer.push(contained_unitig);

        // The contained read now looks for existing unitigs of its own.
        let (available, overlap) = get_num_edges2existing_nodes_advance(
            unitigs,
            unitig_graph,
            read_graph,
            sources,
            coverage_cut,
            read_to_unitig,
            trio_labels,
            parameters,
            t.v,
            chain_buffer,
            ignore_trio_flag,
        );
        if available > 0 {
            return ChainOutcome {
                found: true,
                available_count: available,
                available_overlap: overlap,
                chain_length,
            };
        }

        if chain_length >= threshold_length {
            break;
        }

        end_read = t.v;
        h = get_best_contained_read(
            unitig_graph,
            read_graph,
            sources,
            coverage_cut,
            read_to_unitig,
            trio_labels,
            parameters,
            end_read,
            ignore_trio_flag,
        );
    }

    ChainOutcome {
        found: false,
        available_count: u32::MAX,
        available_overlap: u32::MAX,
        chain_length: u32::MAX,
    }
}

/// `append_rId_to_Unitig`: appends read vertex `read_vertex` to one end of a
/// unitig.
///
/// `unitig_vertex` is the unitig graph vertex whose direction decides the end:
/// even extends the unitig's suffix, odd extends its prefix. The reference
/// recomputes the edge length between the new read and the unitig from the
/// stored overlaps; `t` starts at `-1` so a missing edge keeps that value.
#[allow(clippy::too_many_arguments)]
fn append_read_id_to_unitig(
    read_graph: &Asg,
    unitigs: &mut [MaUtg],
    unitig_vertex: u32,
    mut read_vertex: u32,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    read_to_unitig: &mut ReadToUnitig,
    parameters: RescueParameters,
) {
    let mut t = AsgArc {
        ul: u64::MAX,
        ..AsgArc::default()
    };

    let unitig_index = (unitig_vertex >> 1) as usize;
    let extend_suffix = unitig_vertex & 1 == 0;
    let (beginning, ending) = if extend_suffix {
        (unitigs[unitig_index].end ^ 1, read_vertex)
    } else {
        (read_vertex ^ 1, unitigs[unitig_index].start)
    };
    // `t` is the edge that goes from `beginning` to `ending`.
    get_edge_from_source_into(sources, coverage_cut, parameters, beginning, ending, &mut t);

    if extend_suffix {
        let unitig = &mut unitigs[unitig_index];
        let last = unitig.reads.len() - 1;
        unitig.length = unitig.length.wrapping_sub(unitig.reads[last] as u32) & 0x7fff_ffff;

        let read_length = read_graph.sequences[(read_vertex >> 1) as usize].length;
        unitig
            .reads
            .push((u64::from(read_vertex) << 32) | u64::from(read_length));

        let appended = unitig.reads.len() - 1;
        let edge_length = t.node_length();
        let previous = unitig.reads[appended - 1];
        unitig.reads[appended - 1] = (previous & 0xffff_ffff_0000_0000) | u64::from(edge_length);
        unitig.length = unitig
            .length
            .wrapping_add(unitig.reads[appended - 1] as u32)
            .wrapping_add(unitig.reads[appended] as u32)
            & 0x7fff_ffff;
        unitig.end = read_vertex ^ 1;
    } else {
        read_vertex ^= 1;
        let unitig = &mut unitigs[unitig_index];
        unitig.reads.insert(0, 0);
        let edge_length = t.node_length();
        unitig.reads[0] = (u64::from(read_vertex) << 32) | u64::from(edge_length);
        unitig.length = unitig.length.wrapping_add(unitig.reads[0] as u32) & 0x7fff_ffff;
        unitig.start = read_vertex;
    }

    let read = (read_vertex >> 1) as usize;
    read_to_unitig.set(
        read as u32,
        unitig_vertex >> 1,
        true,
        Some(read_graph.sequences[read].coverage_class),
    );
}

/// `minor_transitive_reduction`: the unitig-level transitive reduction used on
/// the candidate recovery edges.
fn minor_transitive_reduction(
    unitigs: &[MaUtg],
    unitig_graph: &Asg,
    read_to_unitig: &ReadToUnitig,
    edges: &mut [AsgArc],
) {
    for i in 0..edges.len() {
        let t = edges[i];
        let u_id = get_corresponding_u_id(unitigs, unitig_graph, read_to_unitig, t.v);
        if u_id == u32::MAX {
            continue;
        }

        let (arc_start, arc_count) = unitig_graph.arc_range(u_id);
        for j in 0..arc_count as usize {
            let arc = unitig_graph.arcs[arc_start + j];
            if arc.deleted {
                continue;
            }
            let w = arc.v;

            // `w` is now a unitig id while `p->v` is a read id.
            for k in 0..edges.len() {
                if k == i {
                    continue;
                }
                if get_corresponding_u_id(unitigs, unitig_graph, read_to_unitig, edges[k].v) == w {
                    edges[k].deleted = true;
                }
            }
        }
    }
}

/// `lable_all_bubbles`: a full cleanup followed by bubble labelling.
pub fn lable_all_bubbles(graph: &mut Asg, labels: &mut BubbleLabel, minimum_threshold: i64) {
    // Required by the reference: the labeler needs a cleaned, sorted graph.
    graph.cleanup();
    asg_arc_identify_simple_bubbles_multi(graph, labels, false, minimum_threshold);
}

/// `update_unitig`: restores or truncates one end of a unitig after the bubble
/// check deleted the recovered reads again.
pub fn update_unitig(
    step: i64,
    init: i64,
    unitig: &mut MaUtg,
    read_graph: &mut Asg,
    recover_edges: Option<&mut [AsgArc]>,
    update_mode: u32,
) {
    if update_mode == 0 {
        // recover
        if step == 1 {
            unitig.start = (unitig.reads[0] >> 32) as u32;
        }
        if step == -1 {
            let last = unitig.reads.len() - 1;
            unitig.end = ((unitig.reads[last] >> 32) as u32) ^ 1;
        }
        if let Some(edges) = recover_edges {
            // The first node has not been deleted; the others have been.
            let mut position = init - step;
            while position >= 0 && position < unitig.reads.len() as i64 {
                let v = (unitig.reads[position as usize] >> 32) as u32;
                read_graph.sequences[(v >> 1) as usize].deleted = false;
                for edge in edges.iter_mut() {
                    if (v >> 1) == ((edge.ul >> 33) as u32) || (v >> 1) == (edge.v >> 1) {
                        edge.deleted = false;
                    }
                }
                position -= step;
            }
        }
    } else {
        // update: truncate the unitig at the recovered end
        if step == -1 && init != unitig.reads.len() as i64 - 1 {
            let l = read_graph.sequences[(unitig.reads[init as usize] >> 33) as usize].length;
            let new_length = init as usize + 1;
            unitig.reads.truncate(new_length);
            let last = new_length - 1;
            unitig.reads[last] = (unitig.reads[last] & 0xffff_ffff_0000_0000) | u64::from(l);
        }
        if step == 1 && init != 0 {
            unitig.reads.drain(0..init as usize);
        }

        unitig.length = 0;
        for &element in &unitig.reads {
            unitig.length = unitig.length.wrapping_add(element as u32) & 0x7fff_ffff;
        }
    }
}

/// `rescue_contained_reads_aggressive`: reattaches chains of contained reads to
/// the unitig ends they overlap.
#[allow(clippy::too_many_arguments)]
pub fn rescue_contained_reads_aggressive(
    input_unitigs: Option<&mut UnitigGraph>,
    read_graph: &mut Asg,
    sources: &[PafList],
    coverage_cut: &mut [CoverageCut],
    read_to_unitig: &mut ReadToUnitig,
    trio_labels: &[u8],
    parameters: RescueParameters,
    chain_length_threshold: u32,
    is_bubble_check: bool,
    is_primary_check: bool,
    mut new_rtg_edges: Option<&mut Vec<AsgArc>>,
    mut new_rtg_nodes: Option<&mut Vec<u32>>,
    bubble_labels: &mut BubbleLabel,
) {
    let mut new_edges: Vec<AsgArc> = Vec::new();
    let mut rbub_edges: Vec<AsgArc> = Vec::new();
    let mut u_vecs: Vec<u64> = Vec::new();
    let mut chain_buffer: Vec<u32> = Vec::new();
    let mut chain_edges: Vec<AsgArc> = Vec::new();

    let mut owned_unitigs;
    let ug: &mut UnitigGraph = match input_unitigs {
        Some(unitigs) => unitigs,
        None => {
            owned_unitigs = ma_ug_gen(read_graph);
            for sequence in &mut owned_unitigs.graph.sequences {
                sequence.coverage_class = PRIMARY_LABEL;
            }
            &mut owned_unitigs
        }
    };

    // Make every read of a live unitig resolvable again.
    for v in 0..ug.graph.sequences.len() as u32 {
        let u_id = v;
        let unitig = &ug.unitigs[u_id as usize];
        if unitig.reads.is_empty() {
            continue;
        }
        if ug.graph.sequences[v as usize].deleted {
            continue;
        }
        if is_primary_check && ug.graph.sequences[v as usize].coverage_class == ALTER_LABEL {
            continue;
        }
        for k in 0..unitig.reads.len() {
            let read = (unitig.reads[k] >> 33) as u32;
            read_to_unitig.set(
                read,
                u_id,
                true,
                Some(read_graph.sequences[read as usize].coverage_class),
            );
        }
    }

    let vertex_count = ug.graph.sequences.len() as u32 * 2;
    for v in 0..vertex_count {
        let u_id = v >> 1;
        if ug.graph.sequences[u_id as usize].deleted {
            continue;
        }
        if is_primary_check && ug.graph.sequences[u_id as usize].coverage_class == ALTER_LABEL {
            continue;
        }
        {
            let unitig = &ug.unitigs[u_id as usize];
            if unitig.reads.is_empty() {
                continue;
            }
            if unitig.circular || unitig.start == u32::MAX || unitig.end == u32::MAX {
                continue;
            }
        }
        if get_real_length(&ug.graph, v) != 0 {
            continue;
        }

        // `x` is the end read of a tip; find every overlap of `x`.
        let end_read = if v & 1 != 0 {
            ug.unitigs[u_id as usize].start ^ 1
        } else {
            ug.unitigs[u_id as usize].end ^ 1
        };

        let mut available_max = 0u32;
        let mut available_overlap_max = 0u32;
        let mut available_min_chain = u32::MAX;
        let mut best_hit: Option<&PafHit> = None;

        for h in &sources[(end_read >> 1) as usize].hits {
            let outcome = get_contained_reads_chain(
                Some(h),
                sources,
                coverage_cut,
                read_to_unitig,
                &ug.unitigs,
                &ug.graph,
                read_graph,
                trio_labels,
                parameters,
                end_read,
                u_id,
                &mut chain_buffer,
                None,
                chain_length_threshold,
                1,
            );
            if !outcome.found {
                continue;
            }

            let mut is_update = false;
            if outcome.chain_length < available_min_chain {
                is_update = true;
            } else if outcome.chain_length == available_min_chain {
                if outcome.available_count > available_max {
                    is_update = true;
                } else if outcome.available_count == available_max
                    && outcome.available_overlap > available_overlap_max
                {
                    is_update = true;
                }
            }
            if is_update {
                available_min_chain = outcome.chain_length;
                available_max = outcome.available_count;
                available_overlap_max = outcome.available_overlap;
                best_hit = Some(h);
            }
        }

        if available_max == 0 {
            continue;
        }

        // Collect all edges between the contained reads.
        get_contained_reads_chain(
            best_hit,
            sources,
            coverage_cut,
            read_to_unitig,
            &ug.unitigs,
            &ug.graph,
            read_graph,
            trio_labels,
            parameters,
            end_read,
            u_id,
            &mut chain_buffer,
            Some(&mut chain_edges),
            chain_length_threshold,
            1,
        );
        if chain_edges.is_empty() {
            continue;
        }
        let last_contained = *chain_edges.last().unwrap();

        // Edges from the last contained read to other unitigs.
        let mut k = 0usize;
        rbub_edges.clear();
        let mut r_edge = AsgArc::default();
        while get_edge2existing_node_advance(
            &ug.unitigs,
            &ug.graph,
            read_graph,
            sources,
            coverage_cut,
            read_to_unitig,
            trio_labels,
            parameters,
            last_contained.v,
            &[u_id],
            &mut k,
            &mut r_edge,
            1,
        ) {
            rbub_edges.push(r_edge);
        }

        // The reduction here is not the standard transitive reduction.
        minor_transitive_reduction(&ug.unitigs, &ug.graph, read_to_unitig, &mut rbub_edges);

        let mut max_overlap_index = u32::MAX;
        if is_primary_check {
            let mut max_overlap = 0u32;
            for (i, edge) in rbub_edges.iter().enumerate() {
                if edge.deleted {
                    continue;
                }
                if edge.overlap > max_overlap {
                    max_overlap = edge.overlap;
                    max_overlap_index = i as u32;
                }
            }

            if max_overlap_index == u32::MAX {
                continue;
            }
            let t = rbub_edges[max_overlap_index as usize];
            let w = get_corresponding_u_id(&ug.unitigs, &ug.graph, read_to_unitig, t.v);
            if w == u32::MAX {
                continue;
            }
            if get_real_length(&ug.graph, w ^ 1) != 0 {
                continue;
            }
        }

        // Recover all contained reads.
        for edge_index in 0..chain_edges.len() {
            let edge = chain_edges[edge_index];
            let (contained_read, _is_unitig) = read_to_unitig.get(edge.v >> 1);
            u_vecs.push((u64::from(contained_read) << 32) | u64::from(edge.v >> 1));

            // The read and the edge that links it must both come back.
            read_graph.sequences[(edge.v >> 1) as usize].deleted = false;
            coverage_cut[(edge.v >> 1) as usize].deleted = false;
            coverage_cut[(edge.v >> 1) as usize].class = PRIMARY_LABEL;
            append_read_id_to_unitig(
                read_graph,
                &mut ug.unitigs,
                v,
                edge.v,
                sources,
                coverage_cut,
                read_to_unitig,
                parameters,
            );

            let mut forward = AsgArc::default();
            get_edge_from_source_into(
                sources,
                coverage_cut,
                parameters,
                (edge.ul >> 32) as u32,
                edge.v,
                &mut forward,
            );
            new_edges.push(forward);

            let mut backward = AsgArc::default();
            get_edge_from_source_into(
                sources,
                coverage_cut,
                parameters,
                edge.v ^ 1,
                ((edge.ul >> 32) as u32) ^ 1,
                &mut backward,
            );
            new_edges.push(backward);
        }

        if is_primary_check {
            let t = rbub_edges[max_overlap_index as usize];
            let w = get_corresponding_u_id(&ug.unitigs, &ug.graph, read_to_unitig, t.v);
            if w == u32::MAX {
                continue;
            }
            new_edges.push(t);
            ug.graph.append_edge_to_sorted(
                v,
                ug.unitigs[(v >> 1) as usize].length,
                w,
                t.overlap,
                t.strong,
                t.exact_region,
                t.without_large_indel,
            );

            let mut reverse = AsgArc::default();
            get_edge_from_source_into(
                sources,
                coverage_cut,
                parameters,
                t.v ^ 1,
                ((t.ul >> 32) as u32) ^ 1,
                &mut reverse,
            );
            new_edges.push(reverse);
            ug.graph.append_edge_to_sorted(
                w ^ 1,
                ug.unitigs[(w >> 1) as usize].length,
                v ^ 1,
                reverse.overlap,
                reverse.strong,
                reverse.exact_region,
                reverse.without_large_indel,
            );
        } else {
            for edge in &rbub_edges {
                let t = *edge;
                if t.deleted {
                    continue;
                }
                let w = get_corresponding_u_id(&ug.unitigs, &ug.graph, read_to_unitig, t.v);
                if w == u32::MAX {
                    continue;
                }

                new_edges.push(t);
                ug.graph.append_edge_to_sorted(
                    v,
                    ug.unitigs[(v >> 1) as usize].length,
                    w,
                    t.overlap,
                    t.strong,
                    t.exact_region,
                    t.without_large_indel,
                );

                let mut reverse = AsgArc::default();
                get_edge_from_source_into(
                    sources,
                    coverage_cut,
                    parameters,
                    t.v ^ 1,
                    ((t.ul >> 32) as u32) ^ 1,
                    &mut reverse,
                );
                new_edges.push(reverse);
                ug.graph.append_edge_to_sorted(
                    w ^ 1,
                    ug.unitigs[(w >> 1) as usize].length,
                    v ^ 1,
                    reverse.overlap,
                    reverse.strong,
                    reverse.exact_region,
                    reverse.without_large_indel,
                );
            }
        }
    }

    if is_bubble_check {
        for edge in &new_edges {
            read_graph.push_arc(*edge);
        }
        if !new_edges.is_empty() {
            read_graph.index.clear();
            read_graph.is_sorted = false;
            read_graph.cleanup();
        }

        pre_clean(read_graph);

        lable_all_bubbles(read_graph, bubble_labels, parameters.minimum_length_threshold);

        for edge in &new_edges {
            let v = (edge.ul >> 32) as u32;
            let w = edge.v;
            if read_graph.sequences[(v >> 1) as usize].deleted {
                continue;
            }
            if read_graph.sequences[(w >> 1) as usize].deleted {
                continue;
            }
            // Keep the edge when it sits at a bubble.
            if read_graph.sequence_visited[v as usize] != 0
                && read_graph.sequence_visited[(w ^ 1) as usize] != 0
            {
                continue;
            }

            read_graph.delete_arc(v, w, true);
            read_graph.delete_arc(w ^ 1, v ^ 1, true);
        }
        read_graph.cleanup();

        for entry in &u_vecs {
            let mut w = (*entry as u32) << 1;
            if !read_graph.sequences[(w >> 1) as usize].deleted
                && get_real_length(read_graph, w) != 0
            {
                continue;
            }
            if !read_graph.sequences[(w >> 1) as usize].deleted
                && get_real_length(read_graph, w ^ 1) != 0
            {
                continue;
            }
            w >>= 1;
            read_graph.sequences[w as usize].deleted = true;
            coverage_cut[w as usize].deleted = true;
            read_to_unitig.set(*entry as u32, (*entry >> 32) as u32, false, None);
        }
    }

    // The primary check must leave the graph untouched for the caller.
    if is_primary_check {
        if let Some(nodes) = new_rtg_nodes.as_deref_mut() {
            for entry in &u_vecs {
                nodes.push(*entry as u32);
            }
        }

        for edge in &new_edges {
            read_graph.push_arc(*edge);
            if let Some(rtg_edges) = new_rtg_edges.as_deref_mut() {
                rtg_edges.push(*edge);
            }
        }
        if !new_edges.is_empty() {
            read_graph.index.clear();
            read_graph.is_sorted = false;
            read_graph.cleanup();
        }
    }

    read_to_unitig.clear_unitig_entries();
}

/// `rescue_missing_overlaps_aggressive`: links unitig ends to unitigs that the
/// stored overlaps mention but the graph does not contain.
#[allow(clippy::too_many_arguments)]
pub fn rescue_missing_overlaps_aggressive(
    input_unitigs: Option<&mut UnitigGraph>,
    read_graph: &mut Asg,
    sources: &[PafList],
    coverage_cut: &mut [CoverageCut],
    read_to_unitig: &mut ReadToUnitig,
    trio_labels: &[u8],
    parameters: RescueParameters,
    is_bubble_check: bool,
    is_primary_check: bool,
    mut new_rtg_edges: Option<&mut Vec<AsgArc>>,
    bubble_labels: &mut BubbleLabel,
) {
    let mut new_edges: Vec<AsgArc> = Vec::new();
    let mut rbub_edges: Vec<AsgArc> = Vec::new();

    let mut owned_unitigs;
    let ug: &mut UnitigGraph = match input_unitigs {
        Some(unitigs) => unitigs,
        None => {
            owned_unitigs = ma_ug_gen(read_graph);
            for sequence in &mut owned_unitigs.graph.sequences {
                sequence.coverage_class = PRIMARY_LABEL;
            }
            &mut owned_unitigs
        }
    };

    for v in 0..ug.graph.sequences.len() as u32 {
        let u_id = v;
        let unitig = &ug.unitigs[u_id as usize];
        if unitig.reads.is_empty() {
            continue;
        }
        if ug.graph.sequences[v as usize].deleted {
            continue;
        }
        if is_primary_check && ug.graph.sequences[v as usize].coverage_class == ALTER_LABEL {
            continue;
        }
        for k in 0..unitig.reads.len() {
            let read = (unitig.reads[k] >> 33) as u32;
            read_to_unitig.set(
                read,
                u_id,
                true,
                Some(read_graph.sequences[read as usize].coverage_class),
            );
        }
    }

    let vertex_count = ug.graph.sequences.len() as u32 * 2;
    for v in 0..vertex_count {
        let u_id = v >> 1;
        if ug.graph.sequences[u_id as usize].deleted {
            continue;
        }
        if is_primary_check && ug.graph.sequences[u_id as usize].coverage_class == ALTER_LABEL {
            continue;
        }
        {
            let unitig = &ug.unitigs[u_id as usize];
            if unitig.reads.is_empty() {
                continue;
            }
            if unitig.circular || unitig.start == u32::MAX || unitig.end == u32::MAX {
                continue;
            }
        }
        if get_real_length(&ug.graph, v) != 0 {
            continue;
        }

        let end_read = if v & 1 != 0 {
            ug.unitigs[u_id as usize].start ^ 1
        } else {
            ug.unitigs[u_id as usize].end ^ 1
        };

        let mut k = 0usize;
        rbub_edges.clear();
        let mut r_edge = AsgArc::default();
        while get_edge2existing_node_advance(
            &ug.unitigs,
            &ug.graph,
            read_graph,
            sources,
            coverage_cut,
            read_to_unitig,
            trio_labels,
            parameters,
            end_read,
            &[u_id],
            &mut k,
            &mut r_edge,
            1,
        ) {
            rbub_edges.push(r_edge);
        }

        if rbub_edges.is_empty() {
            continue;
        }

        // The reduction here is not the standard transitive reduction.
        minor_transitive_reduction(&ug.unitigs, &ug.graph, read_to_unitig, &mut rbub_edges);

        if is_primary_check {
            let mut max_overlap = 0u32;
            let mut max_overlap_index = u32::MAX;
            for (i, edge) in rbub_edges.iter().enumerate() {
                if edge.deleted {
                    continue;
                }
                if edge.overlap > max_overlap {
                    max_overlap = edge.overlap;
                    max_overlap_index = i as u32;
                }
            }

            if max_overlap_index == u32::MAX {
                continue;
            }
            let t = rbub_edges[max_overlap_index as usize];
            let w = get_corresponding_u_id(&ug.unitigs, &ug.graph, read_to_unitig, t.v);
            if w == u32::MAX {
                continue;
            }
            if get_real_length(&ug.graph, w ^ 1) != 0 {
                continue;
            }

            new_edges.push(t);
            ug.graph.append_edge_to_sorted(
                v,
                unitig_length(&ug.unitigs, v),
                w,
                t.overlap,
                t.strong,
                t.exact_region,
                t.without_large_indel,
            );

            let mut reverse = AsgArc::default();
            get_edge_from_source_into(
                sources,
                coverage_cut,
                parameters,
                t.v ^ 1,
                ((t.ul >> 32) as u32) ^ 1,
                &mut reverse,
            );
            new_edges.push(reverse);
            ug.graph.append_edge_to_sorted(
                w ^ 1,
                unitig_length(&ug.unitigs, w),
                v ^ 1,
                reverse.overlap,
                reverse.strong,
                reverse.exact_region,
                reverse.without_large_indel,
            );
        } else {
            for edge in &rbub_edges {
                let t = *edge;
                if t.deleted {
                    continue;
                }
                let w = get_corresponding_u_id(&ug.unitigs, &ug.graph, read_to_unitig, t.v);
                if w == u32::MAX {
                    continue;
                }

                new_edges.push(t);
                ug.graph.append_edge_to_sorted(
                    v,
                    unitig_length(&ug.unitigs, v),
                    w,
                    t.overlap,
                    t.strong,
                    t.exact_region,
                    t.without_large_indel,
                );

                let mut reverse = AsgArc::default();
                get_edge_from_source_into(
                    sources,
                    coverage_cut,
                    parameters,
                    t.v ^ 1,
                    ((t.ul >> 32) as u32) ^ 1,
                    &mut reverse,
                );
                new_edges.push(reverse);
                ug.graph.append_edge_to_sorted(
                    w ^ 1,
                    unitig_length(&ug.unitigs, w),
                    v ^ 1,
                    reverse.overlap,
                    reverse.strong,
                    reverse.exact_region,
                    reverse.without_large_indel,
                );
            }
        }
    }

    if is_bubble_check {
        for edge in &new_edges {
            read_graph.push_arc(*edge);
        }
        if !new_edges.is_empty() {
            read_graph.index.clear();
            read_graph.is_sorted = false;
            read_graph.cleanup();
        }

        pre_clean(read_graph);

        lable_all_bubbles(read_graph, bubble_labels, parameters.minimum_length_threshold);

        for edge in &new_edges {
            let v = (edge.ul >> 32) as u32;
            let w = edge.v;
            if read_graph.sequences[(v >> 1) as usize].deleted {
                continue;
            }
            if read_graph.sequences[(w >> 1) as usize].deleted {
                continue;
            }
            if read_graph.sequence_visited[v as usize] != 0
                && read_graph.sequence_visited[(w ^ 1) as usize] != 0
            {
                continue;
            }

            read_graph.delete_arc(v, w, true);
            read_graph.delete_arc(w ^ 1, v ^ 1, true);
        }
        read_graph.cleanup();
    }

    if is_primary_check {
        for edge in &new_edges {
            read_graph.push_arc(*edge);
            if let Some(rtg_edges) = new_rtg_edges.as_deref_mut() {
                rtg_edges.push(*edge);
            }
        }
        if !new_edges.is_empty() {
            read_graph.index.clear();
            read_graph.is_sorted = false;
            read_graph.cleanup();
        }
    }

    read_to_unitig.clear_unitig_entries();
}

/// `rescue_missing_overlaps_backward`: walks from the far end of a unitig back
/// through its reads, deleting the reads whose overlaps are missing and
/// recovering them when the walk reaches a unitig.
#[allow(clippy::too_many_arguments)]
pub fn rescue_missing_overlaps_backward(
    input_unitigs: Option<&mut UnitigGraph>,
    read_graph: &mut Asg,
    sources: &[PafList],
    coverage_cut: &mut [CoverageCut],
    read_to_unitig: &mut ReadToUnitig,
    trio_labels: &[u8],
    parameters: RescueParameters,
    backward_steps: u32,
    is_bubble_check: bool,
    is_primary_check: bool,
    bubble_labels: &mut BubbleLabel,
) {
    let mut new_edges: Vec<AsgArc> = Vec::new();
    let mut recover_edges: Vec<AsgArc> = Vec::new();
    let mut rbub_edges: Vec<AsgArc> = Vec::new();
    let mut u_vecs: Vec<u64> = Vec::new();

    let mut owned_unitigs;
    let ug: &mut UnitigGraph = match input_unitigs {
        Some(unitigs) => unitigs,
        None => {
            owned_unitigs = ma_ug_gen(read_graph);
            for sequence in &mut owned_unitigs.graph.sequences {
                sequence.coverage_class = PRIMARY_LABEL;
            }
            &mut owned_unitigs
        }
    };

    for v in 0..ug.graph.sequences.len() as u32 {
        let u_id = v;
        let unitig = &ug.unitigs[u_id as usize];
        if unitig.reads.is_empty() {
            continue;
        }
        if ug.graph.sequences[v as usize].deleted {
            continue;
        }
        if is_primary_check && ug.graph.sequences[v as usize].coverage_class == ALTER_LABEL {
            continue;
        }
        for k in 0..unitig.reads.len() {
            let read = (unitig.reads[k] >> 33) as u32;
            read_to_unitig.set(
                read,
                u_id,
                true,
                Some(read_graph.sequences[read as usize].coverage_class),
            );
        }
    }

    for v_id in 0..ug.graph.sequences.len() as u32 {
        let mut round = 0u32;
        rbub_edges.clear();

        for direction in 0..2u32 {
            let current_backward_steps = if !rbub_edges.is_empty() {
                let unitig = &ug.unitigs[v_id as usize];
                let remaining = (unitig.reads.len() as u32)
                    .wrapping_sub(round)
                    .wrapping_sub(1);
                remaining.min(backward_steps)
            } else {
                backward_steps
            };

            let v = (v_id << 1) | direction;
            let u_id = v >> 1;
            if ug.graph.sequences[u_id as usize].deleted {
                continue;
            }
            if is_primary_check && ug.graph.sequences[u_id as usize].coverage_class == ALTER_LABEL {
                continue;
            }

            {
                let unitig = &ug.unitigs[u_id as usize];
                if unitig.reads.is_empty() {
                    continue;
                }
                if unitig.circular || unitig.start == u32::MAX || unitig.end == u32::MAX {
                    continue;
                }
            }
            if get_real_length(&ug.graph, v) != 0 {
                continue;
            }
            {
                // Both ends must still match the first and last read.
                let unitig = &ug.unitigs[u_id as usize];
                if unitig.start != (unitig.reads[0] >> 32) as u32 {
                    continue;
                }
                let last = unitig.reads.len() - 1;
                if (unitig.end ^ 1) != (unitig.reads[last] >> 32) as u32 {
                    continue;
                }
            }

            let (mut init, step, mode): (i64, i64, u32) = if v & 1 != 0 {
                (0, 1, 1)
            } else {
                (ug.unitigs[u_id as usize].reads.len() as i64 - 1, -1, 0)
            };

            rbub_edges.clear();
            round = 0;
            let mut end_read = 0u32;
            loop {
                let read_count = ug.unitigs[u_id as usize].reads.len() as i64;
                if !(round < current_backward_steps && init >= 0 && init < read_count) {
                    break;
                }

                end_read = (ug.unitigs[u_id as usize].reads[init as usize] >> 32) as u32;
                end_read ^= mode;

                let mut k = 0usize;
                rbub_edges.clear();
                let mut r_edge = AsgArc::default();
                while get_edge2existing_node_advance(
                    &ug.unitigs,
                    &ug.graph,
                    read_graph,
                    sources,
                    coverage_cut,
                    read_to_unitig,
                    trio_labels,
                    parameters,
                    end_read,
                    &[u_id],
                    &mut k,
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
            }

            if rbub_edges.is_empty() {
                continue;
            }

            {
                let unitig = &mut ug.unitigs[u_id as usize];
                if mode == 1 {
                    unitig.start = end_read ^ 1;
                }
                if mode == 0 {
                    unitig.end = end_read ^ 1;
                }
            }
            // Save everything needed to revert.
            u_vecs.push(
                ((((mode as u64) << 31) | init as u64) << 32) | u64::from(u_id),
            );

            // The reduction here is not the standard transitive reduction.
            minor_transitive_reduction(&ug.unitigs, &ug.graph, read_to_unitig, &mut rbub_edges);

            if is_primary_check {
                let mut max_overlap = 0u32;
                let mut max_overlap_index = u32::MAX;
                for (i, edge) in rbub_edges.iter().enumerate() {
                    if edge.deleted {
                        continue;
                    }
                    if edge.overlap > max_overlap {
                        max_overlap = edge.overlap;
                        max_overlap_index = i as u32;
                    }
                }

                if max_overlap_index == u32::MAX {
                    continue;
                }
                let t = rbub_edges[max_overlap_index as usize];
                let w = get_corresponding_u_id(&ug.unitigs, &ug.graph, read_to_unitig, t.v);
                if w == u32::MAX {
                    continue;
                }
                if get_real_length(&ug.graph, w ^ 1) != 0 {
                    continue;
                }

                new_edges.push(t);
                ug.graph.append_edge_to_sorted(
                    v,
                    unitig_length(&ug.unitigs, v),
                    w,
                    t.overlap,
                    t.strong,
                    t.exact_region,
                    t.without_large_indel,
                );

                let mut reverse = AsgArc::default();
                get_edge_from_source_into(
                    sources,
                    coverage_cut,
                    parameters,
                    t.v ^ 1,
                    ((t.ul >> 32) as u32) ^ 1,
                    &mut reverse,
                );
                new_edges.push(reverse);
                ug.graph.append_edge_to_sorted(
                    w ^ 1,
                    unitig_length(&ug.unitigs, w),
                    v ^ 1,
                    reverse.overlap,
                    reverse.strong,
                    reverse.exact_region,
                    reverse.without_large_indel,
                );
            } else {
                // Modify the read graph: every read between the found end and
                // the unitig end loses its arcs but they are kept for recovery.
                let mut position = init - step;
                while position >= 0 && position < ug.unitigs[u_id as usize].reads.len() as i64 {
                    let w = (ug.unitigs[u_id as usize].reads[position as usize] >> 32) as u32;

                    let (arc_start, arc_count) = read_graph.arc_range(w);
                    for k in 0..arc_count as usize {
                        let arc = read_graph.arcs[arc_start + k];
                        if arc.deleted {
                            continue;
                        }
                        recover_edges.push(arc);
                        let mut inverse = AsgArc::default();
                        // The reference prints an error when the reverse arc is
                        // missing and still pushes `t`; keep the default then.
                        if let Some(found) = read_graph.get_arc(arc.v ^ 1, ((arc.ul >> 32) as u32) ^ 1) {
                            inverse = found;
                        }
                        recover_edges.push(inverse);
                    }

                    let (arc_start, arc_count) = read_graph.arc_range(w ^ 1);
                    for k in 0..arc_count as usize {
                        let arc = read_graph.arcs[arc_start + k];
                        if arc.deleted {
                            continue;
                        }
                        recover_edges.push(arc);
                        let mut inverse = AsgArc::default();
                        read_graph.get_arc(arc.v ^ 1, ((arc.ul >> 32) as u32) ^ 1).map(|found| inverse = found);
                        recover_edges.push(inverse);
                    }

                    read_graph.delete_sequence(w >> 1);

                    position -= step;
                }

                for edge in &rbub_edges {
                    let t = *edge;
                    if t.deleted {
                        continue;
                    }
                    let w = get_corresponding_u_id(&ug.unitigs, &ug.graph, read_to_unitig, t.v);
                    if w == u32::MAX {
                        continue;
                    }

                    new_edges.push(t);
                    ug.graph.append_edge_to_sorted(
                        v,
                        unitig_length(&ug.unitigs, v),
                        w,
                        t.overlap,
                        t.strong,
                        t.exact_region,
                        t.without_large_indel,
                    );

                    let mut reverse = AsgArc::default();
                    get_edge_from_source_into(
                        sources,
                        coverage_cut,
                        parameters,
                        t.v ^ 1,
                        ((t.ul >> 32) as u32) ^ 1,
                        &mut reverse,
                    );
                    new_edges.push(reverse);
                    ug.graph.append_edge_to_sorted(
                        w ^ 1,
                        unitig_length(&ug.unitigs, w),
                        v ^ 1,
                        reverse.overlap,
                        reverse.strong,
                        reverse.exact_region,
                        reverse.without_large_indel,
                    );
                }
            }
        }
    }

    if is_bubble_check {
        for edge in &new_edges {
            read_graph.push_arc(*edge);
        }
        if !new_edges.is_empty() {
            read_graph.index.clear();
            read_graph.is_sorted = false;
            read_graph.cleanup();
        }

        pre_clean(read_graph);

        lable_all_bubbles(read_graph, bubble_labels, parameters.minimum_length_threshold);

        for edge in &new_edges {
            let v = (edge.ul >> 32) as u32;
            let w = edge.v;
            if read_graph.sequences[(v >> 1) as usize].deleted {
                continue;
            }
            if read_graph.sequences[(w >> 1) as usize].deleted {
                continue;
            }
            if read_graph.sequence_visited[v as usize] != 0
                && read_graph.sequence_visited[(w ^ 1) as usize] != 0
            {
                continue;
            }

            read_graph.delete_arc(v, w, true);
            read_graph.delete_arc(w ^ 1, v ^ 1, true);
        }

        read_graph.cleanup();

        for edge in &mut recover_edges {
            edge.deleted = true;
        }

        for entry in &u_vecs {
            let mode = (*entry >> 63) as u32;
            let step = if mode == 1 { 1 } else { -1 };
            let init = ((*entry >> 32) & 0x7fff_ffff) as i64;
            let u_id = *entry as u32;
            let end_read = ((ug.unitigs[u_id as usize].reads[init as usize] >> 32) as u32) ^ mode;

            if get_real_length(read_graph, end_read) > 0 {
                update_unitig(
                    step,
                    init,
                    &mut ug.unitigs[u_id as usize],
                    read_graph,
                    Some(&mut recover_edges),
                    1,
                );
            } else {
                update_unitig(
                    step,
                    init,
                    &mut ug.unitigs[u_id as usize],
                    read_graph,
                    Some(&mut recover_edges),
                    0,
                );
            }
        }

        let mut recovered = 0u32;
        for edge in &recover_edges {
            if edge.deleted {
                continue;
            }
            read_graph.push_arc(*edge);
            recovered += 1;
        }

        if recovered != 0 {
            read_graph.index.clear();
            read_graph.is_sorted = false;
        }

        read_graph.cleanup();
    }

    if is_primary_check {
        for entry in &u_vecs {
            let mode = (*entry >> 63) as u32;
            let step = if mode == 1 { 1 } else { -1 };
            let init = ((*entry >> 32) & 0x7fff_ffff) as i64;
            let u_id = *entry as u32;

            update_unitig(
                step,
                init,
                &mut ug.unitigs[u_id as usize],
                read_graph,
                Some(&mut recover_edges),
                1,
            );
        }
    }

    read_to_unitig.clear_unitig_entries();
}

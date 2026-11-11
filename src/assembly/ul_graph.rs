//! The UL string-graph primitives: `recover_contain_g` / `normalize_gou` /
//! `update_sg_uo` (`gfa_ut.cpp` 715 / 777 / 802) and `ma_sg_gen_ul`
//! (`Overlaps.cpp` 2982).
//!
//! These mirror the HiFi `ma_sg_gen` but keep the UL occurrence label (`ou`,
//! taken from `ma_hit_t::bl`) and revive contained reads whose overlaps carry
//! enough UL support.

use crate::assembly::asg::Asg;
use crate::assembly::hit::{ma_hit_to_arc, CoverageCut, ReadToUnitig};
use crate::overlap::paf::PafHit;

/// `OU_MASK`, the 14-bit `ou` field width.
pub const OU_MASK: u32 = 0x3fff;

/// `get_specific_overlap` (`Overlaps.cpp` 1075): the index of the hit in
/// `sources[owner]` whose target is `target` and query is `query`.
fn specific_overlap(sources: &[Vec<PafHit>], owner: u32, query: u32, target: u32) -> Option<usize> {
    sources[owner as usize]
        .iter()
        .position(|hit| hit.query_identifier == query && hit.target_identifier == target)
}

/// One entry of the `recover_contain_g` / `update_sg_uo` worker state.
struct SourceLabelOptions {
    maximum_hang: i32,
    maximum_hang_rate: f32,
    minimum_overlap: i32,
    ul_occ: u32,
}

/// `update_sg_contain` (`gfa_ut.cpp` 685): revive a contained read's overlaps
/// when both directions are strong UL links.
fn update_sg_contain(
    graph: &mut Asg,
    sources: &mut [Vec<PafHit>],
    read_to_unitig: &ReadToUnitig,
    read: usize,
    options: &SourceLabelOptions,
) {
    graph.sequence_visited[read] = 0;
    if !graph.sequences[read].deleted {
        return;
    }
    let (unitig, is_unitig) = read_to_unitig.get(read as u32);
    if unitig == u32::MAX || is_unitig == 1 {
        return;
    }
    let hit_count = sources[read].len();
    for k in 0..hit_count {
        let hit = sources[read][k];
        if hit.block_length < options.ul_occ {
            continue;
        }
        let query = hit.query_identifier;
        let target = hit.target_identifier;
        if graph.sequences[target as usize].deleted {
            let (target_unitig, target_is_unitig) = read_to_unitig.get(target);
            if target_unitig == u32::MAX || target_is_unitig == 1 {
                continue;
            }
        }
        let (result, _) = ma_hit_to_arc(
            &hit,
            graph.sequences[query as usize].length as i32,
            graph.sequences[target as usize].length as i32,
            options.maximum_hang,
            options.maximum_hang_rate,
            options.minimum_overlap,
        );
        if result < 0 {
            continue;
        }
        let Some(reverse) = specific_overlap(sources, target, target, query) else {
            continue;
        };
        let reverse_hit = sources[target as usize][reverse];
        debug_assert_eq!(reverse_hit.block_length, hit.block_length);
        let (reverse_result, _) = ma_hit_to_arc(
            &reverse_hit,
            graph.sequences[target as usize].length as i32,
            graph.sequences[query as usize].length as i32,
            options.maximum_hang,
            options.maximum_hang_rate,
            options.minimum_overlap,
        );
        if reverse_result < 0 {
            continue;
        }
        sources[read][k].deleted = false;
        if !graph.sequences[target as usize].deleted {
            sources[target as usize][reverse].deleted = false;
        }
        graph.sequence_visited[read] = 1;
    }
}

/// `recover_contain_g` (`gfa_ut.cpp` 715): revive the contained reads whose
/// overlaps reach `ul_occ`.
pub fn recover_contain_g(
    graph: &mut Asg,
    sources: &mut [Vec<PafHit>],
    read_to_unitig: &ReadToUnitig,
    maximum_hang: i32,
    maximum_hang_rate: f32,
    minimum_overlap: i32,
    ul_occurrence: u32,
) {
    let options = SourceLabelOptions {
        maximum_hang,
        maximum_hang_rate,
        minimum_overlap,
        ul_occ: ul_occurrence,
    };
    let read_count = graph.sequences.len();
    if graph.sequence_visited.len() < read_count * 2 {
        graph.sequence_visited.resize(read_count * 2, 0);
    }
    for read in 0..read_count {
        update_sg_contain(graph, sources, read_to_unitig, read, &options);
    }
    for read in 0..read_count {
        if graph.sequence_visited[read] != 0 {
            graph.sequences[read].deleted = false;
        }
    }
    for entry in &mut graph.sequence_visited[..read_count * 2] {
        *entry = 0;
    }
}

/// `normalize_gou` (`gfa_ut.cpp` 777): make the two directions of an arc carry
/// the same occurrence label.
pub fn normalize_gou(graph: &mut Asg) {
    let arc_count = graph.arcs.len();
    for arc_index in 0..arc_count {
        let arc = graph.arcs[arc_index];
        if arc.v > arc.query_vertex() {
            continue;
        }
        let v = arc.v ^ 1;
        let w = arc.query_vertex() ^ 1;
        let (start, count) = graph.arc_range(v);
        for neighbouring in start..start + count {
            if graph.arcs[neighbouring as usize].v == w {
                let occurrence = graph.arcs[neighbouring as usize]
                    .occurrence
                    .max(arc.occurrence);
                graph.arcs[neighbouring as usize].occurrence = occurrence;
                graph.arcs[arc_index].occurrence = occurrence;
                break;
            }
        }
    }
}

/// `update_sg_uo` (`gfa_ut.cpp` 802): copy every arc's `ou` from the matching
/// source overlap's `bl`, then print the coverage statistics.
pub fn update_sg_uo(graph: &mut Asg, sources: &[Vec<PafHit>]) {
    for arc_index in 0..graph.arcs.len() {
        let mut arc = graph.arcs[arc_index];
        arc.occurrence = 0;
        if !arc.deleted {
            let query = (arc.ul >> 33) as u32;
            let target = arc.v >> 1;
            for hit in &sources[query as usize] {
                if hit.query_identifier == query && hit.target_identifier == target {
                    arc.occurrence = hit.block_length.min(OU_MASK) as u16;
                    break;
                }
            }
        }
        graph.arcs[arc_index] = arc;
    }

    let mut covered = 0u32;
    let mut kept = 0u32;
    for sequence in 0..graph.sequences.len() {
        if graph.sequences[sequence].deleted {
            continue;
        }
        kept += 1;
        let has = |vertex: u32| {
            graph
                .arcs_of(vertex)
                .iter()
                .any(|arc| !arc.deleted && arc.occurrence != 0)
        };
        if has((sequence as u32) << 1) || has(((sequence as u32) << 1) | 1) {
            covered += 1;
        }
    }
    eprintln!(
        "[M::update_sg_uo::] ==> # gfa reads:{}, # covered gfa reads:{}",
        kept, covered
    );
}

/// `ma_sg_gen_ul` (`Overlaps.cpp` 2982): the UL variant of [`super::
/// graph_build::ma_sg_gen`]; it keeps the `bl` occurrence label and revives
/// contained reads through [`recover_contain_g`].
pub fn ma_sg_gen_ul(
    sources: &mut [Vec<PafHit>],
    coverage_cut: &[CoverageCut],
    read_to_unitig: &ReadToUnitig,
    maximum_hang: i32,
    minimum_overlap: i32,
    integration_fraction: f32,
    ul_occurrence: u32,
) -> Asg {
    let read_count = sources.len();
    let mut graph = Asg::new();
    for read in 0..read_count {
        let length = coverage_cut[read].end.wrapping_sub(coverage_cut[read].start) & 0x7fff_ffff;
        graph.set_sequence(read, length, coverage_cut[read].deleted);
        graph.sequences[read].coverage_class = coverage_cut[read].class;
    }
    graph.sequence_visited = vec![0u8; read_count * 2];

    recover_contain_g(
        &mut graph,
        sources,
        read_to_unitig,
        maximum_hang,
        integration_fraction,
        minimum_overlap,
        ul_occurrence,
    );

    for read in 0..read_count {
        if graph.sequences[read].deleted {
            continue;
        }
        for hit in &sources[read] {
            if hit.deleted {
                continue;
            }
            let query = hit.query_identifier as usize;
            let target = hit.target_identifier as usize;
            let (result, arc) = ma_hit_to_arc(
                hit,
                coverage_cut[query].end.wrapping_sub(coverage_cut[query].start) as i32,
                coverage_cut[target].end.wrapping_sub(coverage_cut[target].start) as i32,
                maximum_hang,
                integration_fraction,
                minimum_overlap,
            );
            debug_assert!(result >= 0);
            let mut arc = arc;
            arc.occurrence = hit.block_length.min(OU_MASK) as u16;
            graph.push_arc(arc);
        }
    }

    graph.cleanup();
    graph.symmetrize();
    graph.read_sequence_count = graph.sequences.len() as u32;
    graph
}

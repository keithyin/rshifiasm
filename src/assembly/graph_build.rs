//! Initial string-graph construction, the counterpart of `ma_sg_gen`,
//! `asg_arc_del_trans` and `gen_init_sg` in `Overlaps.cpp`.
//!
//! `gen_init_sg` runs the whole overlap-to-graph chain of `clean_graph`'s
//! first pass: coverage cut, chimeric-read removal, coordinate cutting,
//! filtering, contained-read folding, arc generation and transitive
//! reduction.

use crate::assembly::asg::Asg;
use crate::assembly::bubble::BubbleLabel;
use crate::assembly::chimeric::detect_chimeric_reads;
use crate::assembly::hit::{
    get_coverage, ma_hit_contained_advance, ma_hit_cut, ma_hit_flt, ma_hit_sub, ma_hit_to_arc,
    CoverageCut, ReadToUnitig,
};
use crate::command_line::options::AssemblyOptions;
use crate::correction::ec_driver::PafList;
use crate::logging;
use crate::overlap::paf::PafHit;
use crate::reads::ReadStore;
use std::io::Write;

/// `asg_arc_del_trans`: Myers' transitive reduction with a length tolerance.
///
/// The reference walks the raw arc arrays, including arcs already marked
/// deleted during this pass, so the same arrays are used here.
pub fn asg_arc_del_trans(graph: &mut Asg, fuzz: i32) -> u32 {
    let start_time = logging::realtime_seconds();
    let vertex_count = graph.sequences.len() * 2;
    let mut mark = vec![0u8; vertex_count];
    let mut reduced = 0u32;

    for v in 0..vertex_count as u32 {
        let (start, count) = graph.arc_range(v);
        if count == 0 {
            continue;
        }

        if graph.sequences[(v >> 1) as usize].deleted {
            for arc in &mut graph.arcs[start..start + count] {
                arc.deleted = true;
                reduced += 1;
            }
            continue;
        }

        // Every out-neighbour is a candidate for reduction.
        for arc in &graph.arcs[start..start + count] {
            mark[arc.v as usize] = 1;
        }

        // The span is sorted by `ul`, so its last arc has the longest node.
        let length_limit = graph.arcs[start + count - 1]
            .node_length()
            .wrapping_add(fuzz as u32);

        for i in 0..count {
            let target = graph.arcs[start + i].v;
            if mark[target as usize] != 1 {
                continue;
            }

            let (target_start, target_count) = graph.arc_range(target);
            let first_length = graph.arcs[start + i].node_length();
            for j in 0..target_count {
                let second = graph.arcs[target_start + j];
                if second.node_length().wrapping_add(first_length) > length_limit {
                    break;
                }
                if mark[second.v as usize] != 0 {
                    mark[second.v as usize] = 2;
                }
            }
        }

        for i in 0..count {
            let target = graph.arcs[start + i].v as usize;
            if mark[target] == 2 {
                graph.arcs[start + i].deleted = true;
                reduced += 1;
            }
            mark[target] = 0;
        }
    }

    eprintln!(
        "[M::asg_arc_del_trans] transitively reduced {} arcs",
        reduced
    );

    if reduced > 0 {
        graph.cleanup();
        graph.symmetrize();
    }

    eprintln!(
        "[M::asg_arc_del_trans] takes {:.2}s\n",
        logging::realtime_seconds() - start_time
    );

    reduced
}

/// `ma_sg_gen`: builds the initial string graph from the cleaned overlaps.
pub fn ma_sg_gen(
    sources: &[Vec<PafHit>],
    coverage_cut: &[CoverageCut],
    maximum_hang: i32,
    minimum_overlap: i32,
    integration_fraction: f32,
) -> Asg {
    let start_time = logging::realtime_seconds();
    let read_count = sources.len();
    let mut graph = Asg::new();

    for i in 0..read_count {
        // `asg_seq_t::len` is a 31-bit field.
        let length = coverage_cut[i]
            .end
            .wrapping_sub(coverage_cut[i].start)
            & 0x7fff_ffff;
        graph.set_sequence(i, length, coverage_cut[i].deleted);
        graph.sequences[i].coverage_class = coverage_cut[i].class;
    }

    graph.sequence_visited = vec![0u8; read_count * 2];

    for i in 0..read_count {
        for hit in &sources[i] {
            if hit.deleted {
                continue;
            }

            let query_cut = coverage_cut[hit.query_identifier as usize];
            let target_cut = coverage_cut[hit.target_identifier as usize];
            let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
            let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;

            let (result, arc) = ma_hit_to_arc(
                hit,
                query_length,
                target_length,
                maximum_hang,
                integration_fraction,
                minimum_overlap,
            );
            if result >= 0 {
                // `ma_sg_gen` does not set the arc's occurrence field: the
                // reference copies an uninitialized stack `asg_arc_t`, so the
                // value is garbage on the default HiFi path and is never read
                // while `is_ou == 0`.
                graph.push_arc(arc);
            } else {
                eprintln!("error");
            }
        }
    }

    graph.cleanup();
    graph.read_sequence_count = graph.sequences.len() as u32;

    eprintln!(
        "[M::ma_sg_gen] takes {:.2} s\n",
        logging::realtime_seconds() - start_time
    );

    graph
}

/// The graph, coverage cuts and label workspace produced by [`gen_init_sg`].
#[derive(Debug, Clone, Default)]
pub struct InitialGraph {
    pub graph: Asg,
    /// `coverage_cut`, the per-read covered interval and class. The reference
    /// hands it on to `clean_graph` and every output writer.
    pub coverage_cut: Vec<CoverageCut>,
    pub bubble_label: BubbleLabel,
    pub coverage: i64,
}

/// `gen_init_sg`: prepares the overlaps and builds the first string graph.
///
/// The reference's optional ultra-long (`all_ul_t`) input is `NULL` in the
/// default configuration, so it is not accepted here.
#[allow(clippy::too_many_arguments)]
pub fn gen_init_sg(
    minimum_depth: i32,
    sources: &mut [Vec<PafHit>],
    read_lengths: &[u64],
    read_to_unitig: &mut ReadToUnitig,
    maximum_hang: i32,
    minimum_overlap: i32,
    gap_fuzz: i32,
    minimum_overlap_length: i64,
    thread_count: i32,
    reads: &ReadStore,
    options: &AssemblyOptions,
) -> InitialGraph {
    let mut coverage_cut =
        ma_hit_sub(minimum_depth, sources, read_lengths, minimum_overlap_length);

    detect_chimeric_reads(
        sources,
        read_lengths,
        &mut coverage_cut,
        (options.maximum_overlap_divergence_final * 2.0) as f32,
        reads,
        options,
    );

    ma_hit_cut(sources, read_lengths, minimum_overlap_length, &mut coverage_cut);
    ma_hit_flt(
        sources,
        &mut coverage_cut,
        maximum_hang,
        minimum_overlap,
        options.maximum_hang_rate,
    );
    ma_hit_contained_advance(
        sources,
        &mut coverage_cut,
        read_to_unitig,
        maximum_hang,
        minimum_overlap,
        options.maximum_hang_rate,
    );

    let mut graph = ma_sg_gen(
        sources,
        &coverage_cut,
        maximum_hang,
        minimum_overlap,
        options.maximum_hang_rate,
    );
    if options.print_debug_gfa != 0 {
        // `if(asm_opt.prt_dbg_gfa) prt_dbg_gfa(sg, "raw", ...)`: a dump of the
        // raw unitig graph with sequences, taken before the transitive
        // reduction (`print_debug_gfa(..., is_seq = 1)`).
        write_raw_debug_gfa(
            &mut graph,
            sources,
            &coverage_cut,
            read_to_unitig,
            reads,
            options,
            maximum_hang,
            minimum_overlap,
        );
    }
    asg_arc_del_trans(&mut graph, gap_fuzz);

    let threshold_count = thread_count.clamp(0, 10) as u32;
    let bubble_label = BubbleLabel::new(threshold_count, graph.sequences.len() as u32);

    let coverage = get_coverage(sources, &coverage_cut, sources.len());

    InitialGraph {
        graph,
        coverage_cut,
        bubble_label,
        coverage,
    }
}

/// The `--prt-raw` dump of `gen_init_sg`: `print_debug_gfa` with `ug == NULL`
/// and `is_seq = 1`, i.e. `ma_ug_gen`, `ma_ug_seq(... is_polish = 0)` and the
/// full GFA writer.
#[allow(clippy::too_many_arguments)]
fn write_raw_debug_gfa(
    read_graph: &mut Asg,
    sources: &[Vec<PafHit>],
    coverage_cut: &[CoverageCut],
    read_to_unitig: &ReadToUnitig,
    reads: &ReadStore,
    options: &AssemblyOptions,
    maximum_hang: i32,
    minimum_overlap: i32,
) {
    let lists: Vec<PafList> = sources
        .iter()
        .map(|hits| PafList {
            hits: hits.clone(),
            ..PafList::default()
        })
        .collect();
    let mut unitig_graph = crate::assembly::unitig::ma_ug_gen(read_graph);
    let edges: Vec<crate::assembly::asg::AsgArc> = Vec::new();
    crate::assembly::unitig::ma_ug_seq(
        &mut unitig_graph,
        read_graph,
        coverage_cut,
        &lists,
        &edges,
        None,
        maximum_hang,
        minimum_overlap,
        false,
        reads,
        options.maximum_hang_rate,
    );
    let prefix = crate::assembly::pipeline::get_outfile_name(options);
    let path = format!("{}.raw.r_utg.gfa", prefix);
    let result = std::fs::File::create(&path).and_then(|file| {
        let mut output = std::io::BufWriter::new(file);
        crate::assembly::unitig::ma_ug_print(
            &unitig_graph,
            read_graph,
            Some(coverage_cut),
            Some(&lists),
            Some(read_to_unitig),
            "utg",
            &mut output,
            reads,
            &reads.trio_labels,
        )?;
        output.flush()
    });
    if let Err(error) = result {
        eprintln!("[ERROR] failed to write {}: {}", path, error);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::asg::AsgArc;
    use crate::assembly::hit::CoverageCut;

    fn hit(
        query: u32,
        query_start: u32,
        query_end: u32,
        target: u32,
        target_start: u32,
        target_end: u32,
    ) -> PafHit {
        PafHit {
            query_identifier: query,
            query_start,
            query_end,
            target_identifier: target,
            target_start,
            target_end,
            matching_length: 1,
            block_length: 200,
            reverse: false,
            without_large_indel: true,
            exact_region: false,
            exact_region_byte: 0,
            deleted: false,
        }
    }

    fn cut(length: u32) -> CoverageCut {
        CoverageCut {
            start: 0,
            end: length,
            deleted: false,
            class: 0,
        }
    }

    fn graph_from_arcs(arcs: &[(u32, u32, u32)]) -> Asg {
        let mut graph = Asg::new();
        for read in 0..4 {
            graph.set_sequence(read, 100, false);
        }
        for &(u, v, length) in arcs {
            graph.push_arc(AsgArc {
                ul: ((u as u64) << 32) | length as u64,
                v,
                overlap: 10,
                ..AsgArc::default()
            });
        }
        graph.cleanup();
        graph
    }

    #[test]
    fn transitive_edges_are_reduced() {
        // 0->4 (30) then 4->6 (60) stays within 0's longest edge (100), so
        // the direct 0->6 (90) is transitive. The mirror 7->1 is reduced the
        // same way.
        let mut graph = graph_from_arcs(&[
            (0, 4, 30),
            (0, 6, 90),
            (0, 2, 100),
            (4, 6, 60),
            (5, 1, 30),
            (7, 1, 90),
            (3, 1, 100),
            (7, 5, 60),
        ]);
        let reduced = asg_arc_del_trans(&mut graph, 0);
        assert_eq!(reduced, 2);
        assert!(graph.get_arc(0, 6).is_none());
        assert!(graph.get_arc(0, 4).is_some());
        assert!(graph.is_symmetric);
    }

    #[test]
    fn graph_generation_keeps_both_halves_of_an_overlap() {
        let sources = vec![
            vec![hit(0, 5, 195, 1, 0, 190)],
            vec![hit(1, 0, 190, 0, 5, 195)],
        ];
        let coverage_cut = vec![cut(200), cut(200)];
        let graph = ma_sg_gen(&sources, &coverage_cut, 1000, 50, 0.8);
        assert_eq!(graph.read_sequence_count, 2);
        assert_eq!(graph.arc_count(0), 1);
        assert_eq!(graph.arc_count(3), 1);
        let arc = graph.arcs_of(0)[0];
        assert_eq!(arc.v, 2);
        assert_eq!(arc.node_length(), 5);
    }

    #[test]
    fn init_graph_runs_the_whole_first_pass() {
        let mut sources = vec![
            vec![hit(0, 5, 195, 1, 0, 190)],
            vec![hit(1, 0, 190, 0, 5, 195)],
        ];
        let read_lengths = [200u64, 200];
        let mut read_to_unitig = ReadToUnitig::new(2);
        read_to_unitig.reset();
        let reads = ReadStore::default();
        let options = AssemblyOptions::default();

        let initial = gen_init_sg(
            options.minimum_overlap_coverage,
            &mut sources,
            &read_lengths,
            &mut read_to_unitig,
            1000,
            50,
            1000,
            50,
            options.thread_count,
            &reads,
            &options,
        );

        assert_eq!(initial.graph.sequences.len(), 2);
        assert_eq!(initial.graph.arc_count(0), 1);
        assert_eq!(initial.graph.arc_count(3), 1);
        assert_eq!(initial.bubble_label.read_count, 2);
        assert_eq!(initial.bubble_label.threshold_count, 1);
        assert_eq!(initial.coverage, 0);
    }
}

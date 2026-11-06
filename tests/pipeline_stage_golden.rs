//! Stage-by-stage validation of the graph pipeline against the reference
//! probe in `/tmp/clean_probe.cpp`.
//!
//! The reference probe replays the cleaning pipeline from the saved overlap
//! index and dumps every intermediate state; the Rust pipeline writes the same
//! dumps when `RSHIFIASM_PIPELINE_DUMP_DIR` is set. This test compares the two
//! trajectories and names the first stage that diverges.
//!
//! Generate the Rust side with:
//!
//! ```text
//! RSHIFIASM_PIPELINE_DUMP_DIR=/tmp/rustprobe \
//!     ./target/release/rshifiasm -t1 -o /tmp/probe_out /path/to/reads.fa.gz
//! ```

use std::path::{Path, PathBuf};

use rshifiasm::assembly::asg::Asg;
use rshifiasm::assembly::hit::{CoverageCut, ReadToUnitig};
use rshifiasm::assembly::probe::{
    load_asg, load_coverage_cut, load_read_to_unitig, load_sources,
};
use rshifiasm::correction::ec_driver::PafList;

/// Stages in pipeline order; `None` marks stages without a graph dump.
const STAGES: &[&str] = &[
    "d00_after_try_rescue",
    "d01_after_weak_clean",
    "d02_after_gen_init_sg",
    "d03_after_cut_tips",
    "d04_after_round_0",
    "d05_after_round_1",
    "d06_after_round_2",
    "d07_after_round_3",
    "d08_after_large_indel",
    "d09_after_cut_length_final",
    "d10_after_semi_circ",
    "d12_after_rescue_contained",
    "d13_after_rescue_missing",
    "d14_after_rescue_backward",
    "d15_after_set_hom_cov",
    "d17_after_ug_ext",
    "d18_after_output_unitig",
    "d20_after_hic_clean",
];

/// The last stage whose overlap lists the Rust pipeline still models exactly;
/// the rescue passes add overlaps on the reference side but only arcs on the
/// Rust side.
const LAST_SOURCE_STAGE: &str = "d10_after_semi_circ";

fn reference_directory() -> PathBuf {
    std::env::var_os("RSHIFIASM_CLEANPROBE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("/tmp/cleanprobe"))
}

fn candidate_directory() -> PathBuf {
    std::env::var_os("RSHIFIASM_PIPELINE_DUMP_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("/tmp/rustprobe"))
}

fn file(directory: &Path, tag: &str, suffix: &str) -> PathBuf {
    directory.join(format!("{}.all.debug.{}.bin", tag, suffix))
}

fn compare_graphs(produced: &Asg, golden: &Asg, label: &str) {
    assert_eq!(produced.sequences.len(), golden.sequences.len(), "{}: vertex count", label);
    assert_eq!(produced.arcs.len(), golden.arcs.len(), "{}: arc count", label);
    assert_eq!(produced.sequence_visited, golden.sequence_visited, "{}: sequence_visited", label);
    assert_eq!(produced.index, golden.index, "{}: index", label);
    assert_eq!(produced.is_sorted, golden.is_sorted, "{}: is_sorted", label);
    assert_eq!(produced.is_symmetric, golden.is_symmetric, "{}: is_symmetric", label);
    assert_eq!(
        produced.read_sequence_count, golden.read_sequence_count,
        "{}: read_sequence_count", label
    );

    for (i, (left, right)) in produced.sequences.iter().zip(&golden.sequences).enumerate() {
        assert_eq!(*left, *right, "{}: sequence {} differs", label, i);
    }
    for (i, (left, right)) in produced.arcs.iter().zip(&golden.arcs).enumerate() {
        // `ma_sg_gen` copies an uninitialized stack `asg_arc_t`, so the
        // reference's 14-bit `ou` field is garbage on the default HiFi path
        // (a few arcs carry leftover stack values, the rest happen to be 0).
        // Every consumer on this path is guarded by `is_ou`, so the field is
        // excluded from the stage comparison.
        assert_eq!(left.ul, right.ul, "{}: arc {} ul differs", label, i);
        assert_eq!(left.v, right.v, "{}: arc {} v differs", label, i);
        assert_eq!(left.overlap, right.overlap, "{}: arc {} ol differs", label, i);
        assert_eq!(left.deleted, right.deleted, "{}: arc {} del differs", label, i);
        assert_eq!(left.strong, right.strong, "{}: arc {} strong differs", label, i);
        assert_eq!(
            left.without_large_indel, right.without_large_indel,
            "{}: arc {} no_l_indel differs", label, i
        );
        assert_eq!(
            left.exact_region, right.exact_region,
            "{}: arc {} el differs", label, i
        );
    }
}

fn compare_coverage(produced: &[CoverageCut], golden: &[CoverageCut], label: &str) {
    assert_eq!(produced.len(), golden.len(), "{}: coverage length", label);
    for (i, (left, right)) in produced.iter().zip(golden).enumerate() {
        assert_eq!(*left, *right, "{}: coverage {} differs", label, i);
    }
}

fn compare_read_map(produced: &ReadToUnitig, golden: &ReadToUnitig, label: &str) {
    assert_eq!(produced.len(), golden.len(), "{}: read map length", label);
    for read in 0..golden.len() as u32 {
        assert_eq!(
            produced.get(read),
            golden.get(read),
            "{}: read map entry {} differs",
            label,
            read
        );
    }
}

fn compare_sources(produced: &[PafList], golden: &[PafList], label: &str) {
    assert_eq!(produced.len(), golden.len(), "{}: source list count", label);
    for (read, (left, right)) in produced.iter().zip(golden).enumerate() {
        assert_eq!(
            left.is_fully_corrected, right.is_fully_corrected,
            "{}: read {} is_fully_corrected", label, read
        );
        assert_eq!(left.is_abnormal, right.is_abnormal, "{}: read {} is_abnormal", label, read);
        assert_eq!(left.hits.len(), right.hits.len(), "{}: read {} hit count", label, read);
        for (index, (left, right)) in left.hits.iter().zip(&right.hits).enumerate() {
            assert_eq!(*left, *right, "{}: read {} hit {}", label, read, index);
        }
    }
}

/// Walks the stages in order and stops at the first divergence, so the panic
/// message names the stage where the Rust pipeline left the reference.
#[test]
fn the_pipeline_stages_reproduce_the_reference_probe() {
    let reference = reference_directory();
    let candidate = candidate_directory();

    if !file(&reference, "d02_after_gen_init_sg", "asg_t").exists() {
        eprintln!(
            "skipping: the reference probe dumps are missing from {}",
            reference.display()
        );
        return;
    }
    if !file(&candidate, "d02_after_gen_init_sg", "asg_t").exists() {
        eprintln!(
            "skipping: the Rust stage dumps are missing from {}; run the pipeline with \
             RSHIFIASM_PIPELINE_DUMP_DIR",
            candidate.display()
        );
        return;
    }

    let source_limit = STAGES
        .iter()
        .position(|stage| *stage == LAST_SOURCE_STAGE)
        .expect("the last source stage must be listed");
    for (index, tag) in STAGES.iter().enumerate() {
        let mut compared = false;

        // The graph dump is absent from the stages that run before the first
        // unitig graph exists.
        let golden_graph_path = file(&reference, tag, "asg_t");
        let produced_graph_path = file(&candidate, tag, "asg_t");
        match (produced_graph_path.exists(), golden_graph_path.exists()) {
            (true, true) => {
                compare_graphs(
                    &load_asg(&produced_graph_path),
                    &load_asg(&golden_graph_path),
                    tag,
                );
                compared = true;
            }
            (false, false) => {}
            (produced, _) => panic!(
                "{}: the graph dump exists on {} side only",
                tag,
                if produced { "the Rust" } else { "the reference" }
            ),
        }

        let golden_coverage_path = file(&reference, tag, "coverage_cut");
        let produced_coverage_path = file(&candidate, tag, "coverage_cut");
        match (produced_coverage_path.exists(), golden_coverage_path.exists()) {
            (true, true) => {
                compare_coverage(
                    &load_coverage_cut(&produced_coverage_path),
                    &load_coverage_cut(&golden_coverage_path),
                    tag,
                );
                compared = true;
            }
            (false, false) => {}
            (produced, _) => panic!(
                "{}: the coverage dump exists on {} side only",
                tag,
                if produced { "the Rust" } else { "the reference" }
            ),
        }

        let golden_map_path = file(&reference, tag, "ruIndex");
        let produced_map_path = file(&candidate, tag, "ruIndex");
        match (produced_map_path.exists(), golden_map_path.exists()) {
            (true, true) => {
                let (produced_map, produced_trio, _, _) = load_read_to_unitig(&produced_map_path);
                let (golden_map, golden_trio, _, _) = load_read_to_unitig(&golden_map_path);
                compare_read_map(&produced_map, &golden_map, tag);
                assert_eq!(produced_trio, golden_trio, "{}: trio labels", tag);
                compared = true;
            }
            (false, false) => {}
            (produced, _) => panic!(
                "{}: the read map dump exists on {} side only",
                tag,
                if produced { "the Rust" } else { "the reference" }
            ),
        }

        // From the rescue passes on, the reference appends the rescued
        // overlaps to the source lists while the Rust port only adds arcs, so
        // only the stages before that point can be compared hit by hit.
        if index <= source_limit {
            let golden_source_path = file(&reference, tag, "source");
            let produced_source_path = file(&candidate, tag, "source");
            match (produced_source_path.exists(), golden_source_path.exists()) {
                (true, true) => {
                    compare_sources(
                        &load_sources(&produced_source_path),
                        &load_sources(&golden_source_path),
                        tag,
                    );
                    compared = true;
                }
                (false, false) => {}
                (produced, _) => panic!(
                    "{}: the overlap dump exists on {} side only",
                    tag,
                    if produced { "the Rust" } else { "the reference" }
                ),
            }
        }

        if compared {
            eprintln!("stage {} matches", tag);
        } else {
            eprintln!("skipping stage {}: no dumps on either side", tag);
        }
    }
}

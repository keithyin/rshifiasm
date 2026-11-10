//! Byte-for-byte validation of the rescue family against the reference state
//! dumps.
//!
//! The stepwise probe (`/tmp/clean_probe.cpp`) replays `ul_clean_gfa` with a
//! state dump after every step. This test takes the dump written after
//! `asg_cut_semi_circ` (d10), runs the three rescue passes of `post_rescue`
//! with the same arguments as the probe, and compares the resulting graph,
//! coverage cuts and read map against the reference dumps d12, d13 and d14.
//!
//! The dumps live in a scratch directory and are not part of the repository,
//! so the test reports a skip when they are absent.

use std::fs;
use std::path::{Path, PathBuf};

use rshifiasm::assembly::asg::{Asg, AsgArc, AsgSeq};
use rshifiasm::assembly::bubble::BubbleLabel;
use rshifiasm::assembly::hit::{CoverageCut, ReadToUnitig};
use rshifiasm::assembly::rescue::{
    rescue_contained_reads_aggressive, rescue_missing_overlaps_aggressive,
    rescue_missing_overlaps_backward, RescueParameters,
};
use rshifiasm::correction::ec_driver::PafList;
use rshifiasm::overlap::paf::PafHit;

/// `asm_opt.max_hang_Len`.
const MAXIMUM_HANG: i32 = 1000;
/// `asm_opt.min_overlap_Len`.
const MINIMUM_OVERLAP: i32 = 50;
/// `asm_opt.max_hang_rate`.
const MAX_HANG_RATE: f32 = 0.8;
/// `min_thres = asm_opt.max_short_tip + 1`.
const MINIMUM_LENGTH_THRESHOLD: i64 = 4;

fn dump_directory() -> PathBuf {
    std::env::var_os("RSHIFIASM_CLEANPROBE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("/tmp/cleanprobe"))
}

fn read_u32(data: &[u8], position: &mut usize) -> u32 {
    let value = u32::from_le_bytes(data[*position..*position + 4].try_into().unwrap());
    *position += 4;
    value
}

fn read_u16(data: &[u8], position: &mut usize) -> u16 {
    let value = u16::from_le_bytes(data[*position..*position + 2].try_into().unwrap());
    *position += 2;
    value
}

fn read_u64(data: &[u8], position: &mut usize) -> u64 {
    let value = u64::from_le_bytes(data[*position..*position + 8].try_into().unwrap());
    *position += 8;
    value
}

fn read_u8(data: &[u8], position: &mut usize) -> u8 {
    let value = data[*position];
    *position += 1;
    value
}

/// `load_asg_t`.
fn load_asg(path: &Path) -> Asg {
    let data = fs::read(path).expect("the asg dump must be readable");
    let mut position = 0usize;
    let count_arc = read_u32(&data, &mut position) as usize;
    let _capacity_arc = read_u32(&data, &mut position);
    let is_sorted = read_u32(&data, &mut position) != 0;
    let count_sequence = read_u32(&data, &mut position) as usize;
    let _capacity_sequence = read_u32(&data, &mut position);
    let is_symmetric = read_u32(&data, &mut position) != 0;
    let read_sequence_count = read_u32(&data, &mut position);

    let mut sequence_visited = vec![0u8; count_sequence * 2];
    sequence_visited.copy_from_slice(&data[position..position + count_sequence * 2]);
    position += count_sequence * 2;

    let mut index = vec![0u64; count_sequence * 2];
    for entry in index.iter_mut() {
        *entry = read_u64(&data, &mut position);
    }

    let mut arcs = Vec::with_capacity(count_arc);
    for _ in 0..count_arc {
        let ul = read_u64(&data, &mut position);
        let v = read_u32(&data, &mut position);
        let overlap_word = read_u32(&data, &mut position);
        let field_word = read_u16(&data, &mut position);
        let exact_region = read_u8(&data, &mut position);
        position += 5;
        arcs.push(AsgArc {
            ul,
            v,
            overlap: overlap_word & 0x7fff_ffff,
            deleted: overlap_word >> 31 != 0,
            occurrence: field_word & 0x3fff,
            strong: (field_word >> 14) & 1 != 0,
            without_large_indel: (field_word >> 15) & 1 != 0,
            exact_region,
        });
    }

    let mut sequences = Vec::with_capacity(count_sequence);
    for _ in 0..count_sequence {
        let length_word = read_u32(&data, &mut position);
        let coverage_class = read_u8(&data, &mut position);
        position += 3;
        sequences.push(AsgSeq {
            length: length_word & 0x7fff_ffff,
            deleted: length_word >> 31 != 0,
            coverage_class,
        });
    }

    assert_eq!(position, data.len(), "the asg dump has trailing bytes");
    Asg {
        arcs,
        index,
        is_sorted,
        is_symmetric,
        sequences,
        read_sequence_count,
        sequence_visited,
    }
}

/// `load_ma_hit_ts`.
fn load_sources(path: &Path) -> Vec<PafList> {
    let data = fs::read(path).expect("the overlap dump must be readable");
    let mut position = 0usize;
    let read_count = read_u64(&data, &mut position) as usize;
    let mut lists = Vec::with_capacity(read_count);

    for _ in 0..read_count {
        let is_fully_corrected = read_u8(&data, &mut position) != 0;
        let is_abnormal = read_u8(&data, &mut position);
        let length = read_u32(&data, &mut position) as usize;
        let mut hits = Vec::with_capacity(length);
        for _ in 0..length {
            let qns = read_u64(&data, &mut position);
            let query_end = read_u32(&data, &mut position);
            let target_identifier = read_u32(&data, &mut position);
            let target_start = read_u32(&data, &mut position);
            let target_end = read_u32(&data, &mut position);
            let field_word = read_u32(&data, &mut position);
            let block_word = read_u32(&data, &mut position);
            let exact_region = read_u8(&data, &mut position);
            let without_large_indel = read_u8(&data, &mut position) != 0;
            position += 6;

            hits.push(PafHit {
                query_identifier: (qns >> 32) as u32,
                query_start: qns as u32,
                query_end,
                target_identifier,
                target_start,
                target_end,
                matching_length: (field_word >> 30) & 1,
                block_length: block_word & 0x7fff_ffff,
                reverse: field_word >> 31 != 0,
                without_large_indel,
                exact_region: exact_region != 0,
                exact_region_byte: exact_region,
                deleted: block_word >> 31 != 0,
            });
        }
        lists.push(PafList {
            hits,
            is_fully_corrected,
            is_abnormal,
        });
    }

    assert_eq!(position, data.len(), "the overlap dump has trailing bytes");
    lists
}

/// `load_coverage_cut`.
fn load_coverage_cut(path: &Path) -> Vec<CoverageCut> {
    let data = fs::read(path).expect("the coverage dump must be readable");
    let mut position = 0usize;
    let read_count = read_u64(&data, &mut position) as usize;
    let mut cuts = Vec::with_capacity(read_count);
    for _ in 0..read_count {
        let start_word = read_u32(&data, &mut position);
        let end = read_u32(&data, &mut position);
        let coverage_class = read_u8(&data, &mut position);
        position += 3;
        cuts.push(CoverageCut {
            start: start_word & 0x7fff_ffff,
            end,
            deleted: start_word >> 31 != 0,
            class: coverage_class,
        });
    }
    assert_eq!(position, data.len(), "the coverage dump has trailing bytes");
    cuts
}

/// `load_ruIndex`; the trio labels travel in the same file.
fn load_read_to_unitig(path: &Path) -> (ReadToUnitig, Vec<u8>, u32, u32) {
    let data = fs::read(path).expect("the read map dump must be readable");
    let mut position = 0usize;
    let length = read_u64(&data, &mut position) as usize;

    let mut values = Vec::with_capacity(length);
    for _ in 0..length {
        values.push(read_u32(&data, &mut position));
    }
    let mut trio_labels = vec![0u8; length];
    trio_labels.copy_from_slice(&data[position..position + length]);
    position += length;
    let hom_global_coverage_set = read_u32(&data, &mut position);
    let hom_global_coverage = read_u32(&data, &mut position);
    assert_eq!(position, data.len(), "the read map dump has trailing bytes");

    let mut read_to_unitig = ReadToUnitig::new(length);
    for (read, &value) in values.iter().enumerate() {
        read_to_unitig.set(read as u32, value & 0x7fff_ffff, value >> 31 != 0, None);
    }
    (
        read_to_unitig,
        trio_labels,
        hom_global_coverage_set,
        hom_global_coverage,
    )
}

fn parameters() -> RescueParameters {
    RescueParameters {
        maximum_hang: MAXIMUM_HANG,
        maximum_hang_rate: MAX_HANG_RATE,
        minimum_overlap: MINIMUM_OVERLAP,
        minimum_length_threshold: MINIMUM_LENGTH_THRESHOLD,
    }
}

fn compare_graphs(produced: &Asg, golden: &Asg, label: &str) {
    assert_eq!(
        produced.sequences.len(),
        golden.sequences.len(),
        "{}: vertex count",
        label
    );
    assert_eq!(produced.arcs.len(), golden.arcs.len(), "{}: arc count", label);
    assert_eq!(
        produced.sequence_visited, golden.sequence_visited,
        "{}: sequence_visited",
        label
    );
    assert_eq!(produced.index, golden.index, "{}: index", label);
    assert_eq!(
        produced.is_sorted, golden.is_sorted,
        "{}: is_sorted",
        label
    );
    assert_eq!(
        produced.is_symmetric, golden.is_symmetric,
        "{}: is_symmetric",
        label
    );
    assert_eq!(
        produced.read_sequence_count, golden.read_sequence_count,
        "{}: read_sequence_count",
        label
    );

    for (i, (left, right)) in produced.sequences.iter().zip(&golden.sequences).enumerate() {
        assert_eq!(*left, *right, "{}: sequence {} differs", label, i);
    }
    for (i, (left, right)) in produced.arcs.iter().zip(&golden.arcs).enumerate() {
        assert_eq!(*left, *right, "{}: arc {} differs", label, i);
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

#[test]
fn the_rescue_family_reproduces_the_reference_dumps() {
    let dump = dump_directory();
    let asg_path = dump.join("d10_after_semi_circ.all.debug.asg_t.bin");
    if !asg_path.exists() {
        eprintln!("skipping: the clean probe dumps are missing");
        return;
    }

    let mut graph = load_asg(&asg_path);
    let sources = load_sources(&dump.join("d10_after_semi_circ.all.debug.source.bin"));
    let mut coverage_cut = load_coverage_cut(&dump.join("d10_after_semi_circ.all.debug.coverage_cut.bin"));
    let (mut read_to_unitig, trio_labels, _hom_set, _hom_cov) =
        load_read_to_unitig(&dump.join("d10_after_semi_circ.all.debug.ruIndex.bin"));

    let read_count = graph.sequences.len() as u32;
    let mut bubble_labels = BubbleLabel::new(1, read_count);

    rescue_contained_reads_aggressive(
        None,
        &mut graph,
        &sources,
        &mut coverage_cut,
        &mut read_to_unitig,
        &trio_labels,
        parameters(),
        10,
        true,
        false,
        None,
        None,
        &mut bubble_labels,
    );

    let golden_graph = load_asg(&dump.join("d12_after_rescue_contained.all.debug.asg_t.bin"));
    let golden_coverage =
        load_coverage_cut(&dump.join("d12_after_rescue_contained.all.debug.coverage_cut.bin"));
    let (golden_map, golden_trio, _, _) =
        load_read_to_unitig(&dump.join("d12_after_rescue_contained.all.debug.ruIndex.bin"));
    compare_graphs(&graph, &golden_graph, "d12");
    compare_coverage(&coverage_cut, &golden_coverage, "d12");
    compare_read_map(&read_to_unitig, &golden_map, "d12");
    assert_eq!(trio_labels, golden_trio, "d12: trio labels");

    rescue_missing_overlaps_aggressive(
        None,
        &mut graph,
        &sources,
        &mut coverage_cut,
        &mut read_to_unitig,
        &trio_labels,
        parameters(),
        true,
        false,
        None,
        &mut bubble_labels,
    );

    let golden_graph = load_asg(&dump.join("d13_after_rescue_missing.all.debug.asg_t.bin"));
    let golden_coverage =
        load_coverage_cut(&dump.join("d13_after_rescue_missing.all.debug.coverage_cut.bin"));
    let (golden_map, _, _, _) =
        load_read_to_unitig(&dump.join("d13_after_rescue_missing.all.debug.ruIndex.bin"));
    compare_graphs(&graph, &golden_graph, "d13");
    compare_coverage(&coverage_cut, &golden_coverage, "d13");
    compare_read_map(&read_to_unitig, &golden_map, "d13");

    rescue_missing_overlaps_backward(
        None,
        &mut graph,
        &sources,
        &mut coverage_cut,
        &mut read_to_unitig,
        &trio_labels,
        parameters(),
        10,
        true,
        false,
        &mut bubble_labels,
    );

    let golden_graph = load_asg(&dump.join("d14_after_rescue_backward.all.debug.asg_t.bin"));
    let golden_coverage =
        load_coverage_cut(&dump.join("d14_after_rescue_backward.all.debug.coverage_cut.bin"));
    let (golden_map, _, _, _) =
        load_read_to_unitig(&dump.join("d14_after_rescue_backward.all.debug.ruIndex.bin"));
    compare_graphs(&graph, &golden_graph, "d14");
    compare_coverage(&coverage_cut, &golden_coverage, "d14");
    compare_read_map(&read_to_unitig, &golden_map, "d14");

    // The overlap lists must survive the rescue untouched.
    let golden_sources =
        load_sources(&dump.join("d14_after_rescue_backward.all.debug.source.bin"));
    assert_eq!(sources.len(), golden_sources.len(), "source list count");
    for (i, (left, right)) in sources.iter().zip(&golden_sources).enumerate() {
        assert_eq!(left.hits, right.hits, "source list {} differs", i);
        assert_eq!(left.is_fully_corrected, right.is_fully_corrected);
        assert_eq!(left.is_abnormal, right.is_abnormal);
    }
}

//! Byte-for-byte validation of the unitig stage against the reference state
//! dumped after `ul_clean_gfa`.
//!
//! The reference run with `HA_F_VERBOSE_GFA` writes
//! `<output>.all.debug.{source,reverse,coverage_cut,ruIndex,asg_t}` right after
//! `ul_clean_gfa` returns. Because `clean_sg_by_utg` leaves exactly the unitig
//! arcs live, running `ma_ug_gen` on that state reproduces the `r_utg` GFA that
//! the same run wrote a moment earlier; the test proves it for the chr11-2M
//! corpus.
//!
//! The dump lives in a scratch directory and is not part of the repository, so
//! the test reports a skip when it is absent.

use std::fs;
use std::path::{Path, PathBuf};

use rshifiasm::assembly::asg::{Asg, AsgArc, AsgSeq, MaUtg};
use rshifiasm::assembly::bubble::BubbleLabel;
use rshifiasm::assembly::bubble_identify::asg_pop_bubble_primary_trio;
use rshifiasm::assembly::hap_cov::{init_trans_chain, HapCov};
use rshifiasm::assembly::hit::{CoverageCut, ReadToUnitig, PRIMARY_LABEL};
use rshifiasm::assembly::tangle::resolve_tangles;
use rshifiasm::assembly::trio_clean::{
    asg_arc_cut_trio_long_equal_tips_assembly,
    asg_arc_cut_trio_long_equal_tips_assembly_complex, asg_arc_cut_trio_long_tip_primary,
    asg_arc_cut_trio_long_tip_primary_complex, clean_primary_untig_graph, cut_trio_tip_primary,
    detect_chimeric_by_topo, drop_semi_circle, unitig_arc_del_short_diploid_by_length,
    unitig_arc_del_short_diploid_by_length_topo,
};
use rshifiasm::assembly::unitig::{
    adjust_utg_advance, copy_read_graph, copy_untig_graph, delete_useless_nodes,
    discard_small_ctg, hic_clean, ma_ug_gen, ma_ug_gen_primary, ma_ug_print, ma_ug_print_bed,
    ma_ug_print_simple, ma_ug_seq, output_contig_graph_primary_pre, recover_utg_by_coverage,
    renew_utg, UnitigGraph,
};
use rshifiasm::command_line::options::AssemblyOptions;
use rshifiasm::correction::ec_driver::PafList;
use rshifiasm::overlap::paf::PafHit;
use rshifiasm::reads::ReadStore;
use rshifiasm::sequence::encoding::PackedSequence;

/// `asm_opt.max_hang_Len`.
const MAXIMUM_HANG: i32 = 1000;
/// `asm_opt.min_overlap_Len`.
const MINIMUM_OVERLAP: i32 = 50;
/// `asm_opt.max_hang_rate`.
const MAX_HANG_RATE: f32 = 0.8;
/// `DROP` from `Process_Read.h`: the flag of a read dropped by the trio logic.
const DROP: u32 = 5;

fn dump_directory() -> PathBuf {
    std::env::var_os("RSHIFIASM_UNITIG_DUMP_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("/tmp/dbgref"))
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
        // The C structure is 24 bytes; five padding bytes follow `el`.
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

/// `load_ma_hit_ts`, the raw `ma_hit_t` records behind `write_debug_ma_hit_ts`.
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

/// `load_coverage_cut`: raw `ma_sub_t` records.
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
fn load_read_to_unitig(path: &Path) -> (ReadToUnitig, Vec<u8>) {
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
    let _hom_global_coverage_set = read_u32(&data, &mut position);
    let _hom_global_coverage = read_u32(&data, &mut position);
    assert_eq!(position, data.len(), "the read map dump has trailing bytes");

    let mut read_to_unitig = ReadToUnitig::new(length);
    for (read, &value) in values.iter().enumerate() {
        read_to_unitig.set(
            read as u32,
            value & 0x7fff_ffff,
            value >> 31 != 0,
            None,
        );
    }
    (read_to_unitig, trio_labels)
}

/// The corrected reads, in `R_INF` order, written by `--write-ec` as one
/// sequence line per record.
fn load_corrected_reads(path: &Path) -> ReadStore {
    let data = fs::read(path).expect("the corrected FASTA must be readable");
    let mut store = ReadStore::default();
    let mut lines = data.split(|&byte| byte == b'\n');
    while let Some(header) = lines.next() {
        if header.is_empty() {
            continue;
        }
        assert_eq!(header[0], b'>', "every record must start with '>'");
        let name = String::from_utf8_lossy(&header[1..]).to_string();
        let sequence = lines.next().expect("every record must have a sequence");
        let sequence = sequence.strip_suffix(b"\r").unwrap_or(sequence);
        store.read_names.push(name);
        store.read_lengths.push(sequence.len() as u64);
        store
            .sequences
            .push(PackedSequence::from_ascii(sequence));
        store.trio_labels.push(0);
        store.total_reads += 1;
        store.total_bases += sequence.len() as u64;
    }
    store
}

fn repository_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

#[test]
fn the_unitig_stage_reproduces_the_reference_r_utg() {
    let dump = dump_directory();
    let asg_path = dump.join("reference.all.debug.asg_t.bin");
    let golden_path = repository_root().join("reference/default_repeat.bp.r_utg.gfa");
    let output_sequence_golden = repository_root().join("reference/default_repeat.bp.r_utg.noseq.gfa");
    let low_quality_golden = repository_root().join("reference/default_repeat.bp.r_utg.lowQ.bed");
    let corrected_reads = repository_root().join("reference/write_paf_ec.ec.fa");
    if !asg_path.exists() || !golden_path.exists() || !corrected_reads.exists() {
        eprintln!(
            "skipping: the reference dump ({}) or the goldens are missing",
            asg_path.display()
        );
        return;
    }

    let mut graph = load_asg(&asg_path);
    let sources = load_sources(&dump.join("reference.all.debug.source.bin"));
    let coverage_cut = load_coverage_cut(&dump.join("reference.all.debug.coverage_cut.bin"));
    let (read_to_unitig, trio_labels) =
        load_read_to_unitig(&dump.join("reference.all.debug.ruIndex.bin"));

    let mut reads = load_corrected_reads(&corrected_reads);
    assert_eq!(reads.total_reads as usize, trio_labels.len());
    reads.trio_labels = trio_labels;
    assert_eq!(
        reads.total_reads as usize,
        coverage_cut.len(),
        "every read must have a coverage cut"
    );

    let mut unitig_graph = ma_ug_gen(&mut graph);
    let edges = Vec::new();
    ma_ug_seq(
        &mut unitig_graph,
        &mut graph,
        &coverage_cut,
        &sources,
        &edges,
        None,
        MAXIMUM_HANG,
        MINIMUM_OVERLAP,
        true,
        &reads,
        MAX_HANG_RATE,
    );

    let mut produced: Vec<u8> = Vec::new();
    ma_ug_print(
        &unitig_graph,
        &graph,
        Some(&coverage_cut),
        Some(&sources),
        Some(&read_to_unitig),
        "utg",
        &mut produced,
        &reads,
        &reads.trio_labels,
    )
    .unwrap();

    let golden = fs::read(&golden_path).unwrap();
    assert_lines_equal(&produced, &golden, "r_utg.gfa");

    let mut produced_simple: Vec<u8> = Vec::new();
    ma_ug_print_simple(
        &unitig_graph,
        &graph,
        Some(&coverage_cut),
        Some(&sources),
        Some(&read_to_unitig),
        "utg",
        &mut produced_simple,
        &reads,
        &reads.trio_labels,
    )
    .unwrap();

    if output_sequence_golden.exists() {
        let golden_simple = fs::read(&output_sequence_golden).unwrap();
        assert_lines_equal(&produced_simple, &golden_simple, "r_utg.noseq.gfa");
    }

    if low_quality_golden.exists() {
        let mut produced_bed: Vec<u8> = Vec::new();
        ma_ug_print_bed(
            &unitig_graph,
            &graph,
            &coverage_cut,
            &sources,
            &edges,
            MAXIMUM_HANG,
            MINIMUM_OVERLAP,
            70,
            "utg",
            &mut produced_bed,
            &reads,
            MAX_HANG_RATE,
        )
        .unwrap();
        let golden_bed = fs::read(&low_quality_golden).unwrap();
        assert_lines_equal(&produced_bed, &golden_bed, "r_utg.lowQ.bed");
    }
}

#[test]
fn the_primary_unitig_writer_reproduces_the_reference_p_utg() {
    // The probe writes its post-`ul_clean_gfa` state (d15, after
    // `set_hom_global_coverage`); the later steps leave it untouched on this
    // corpus, so the `p_utg` writer can run directly on it.
    let dump = dump_directory();
    let asg_path = dump.join("reference.all.debug.asg_t.bin");
    let alternative_dump = PathBuf::from("/tmp/cleanprobe");
    let alternative_asg = alternative_dump.join("d15_after_set_hom_cov.all.debug.asg_t.bin");
    let (asg_path, sources_path, coverage_path, map_path) = if asg_path.exists() {
        (
            asg_path.clone(),
            dump.join("reference.all.debug.source.bin"),
            dump.join("reference.all.debug.coverage_cut.bin"),
            dump.join("reference.all.debug.ruIndex.bin"),
        )
    } else if alternative_asg.exists() {
        (
            alternative_asg,
            alternative_dump.join("d15_after_set_hom_cov.all.debug.source.bin"),
            alternative_dump.join("d15_after_set_hom_cov.all.debug.coverage_cut.bin"),
            alternative_dump.join("d15_after_set_hom_cov.all.debug.ruIndex.bin"),
        )
    } else {
        eprintln!("skipping: the post-clean dump is missing");
        return;
    };

    let golden_path = repository_root().join("reference/default_repeat.bp.p_utg.gfa");
    let output_sequence_golden =
        repository_root().join("reference/default_repeat.bp.p_utg.noseq.gfa");
    let low_quality_golden = repository_root().join("reference/default_repeat.bp.p_utg.lowQ.bed");
    let corrected_reads = repository_root().join("reference/write_paf_ec.ec.fa");
    if !golden_path.exists() || !corrected_reads.exists() {
        eprintln!("skipping: the p_utg goldens are missing");
        return;
    }

    let mut graph = load_asg(&asg_path);
    let sources = load_sources(&sources_path);
    let coverage_cut = load_coverage_cut(&coverage_path);
    let (read_to_unitig, trio_labels) = load_read_to_unitig(&map_path);

    let mut reads = load_corrected_reads(&corrected_reads);
    assert_eq!(reads.total_reads as usize, trio_labels.len());
    reads.trio_labels = trio_labels;
    assert_eq!(reads.total_reads as usize, coverage_cut.len());

    let output_directory = std::env::temp_dir().join("rshifiasm_p_utg_golden");
    fs::create_dir_all(&output_directory).unwrap();
    let prefix = output_directory.join("reference");
    let prefix = prefix.to_str().unwrap();

    output_contig_graph_primary_pre(
        &mut graph,
        &coverage_cut,
        prefix,
        &sources,
        &read_to_unitig,
        0,
        MAXIMUM_HANG,
        MINIMUM_OVERLAP,
        &reads,
        MAX_HANG_RATE,
        70,
    )
    .unwrap();

    let produced = fs::read(output_directory.join("reference.p_utg.gfa")).unwrap();
    let golden = fs::read(&golden_path).unwrap();
    assert_lines_equal(&produced, &golden, "p_utg.gfa");

    if output_sequence_golden.exists() {
        let produced = fs::read(output_directory.join("reference.p_utg.noseq.gfa")).unwrap();
        let golden = fs::read(&output_sequence_golden).unwrap();
        assert_lines_equal(&produced, &golden, "p_utg.noseq.gfa");
    }
    if low_quality_golden.exists() {
        let produced = fs::read(output_directory.join("reference.p_utg.lowQ.bed")).unwrap();
        let golden = fs::read(&low_quality_golden).unwrap();
        assert_lines_equal(&produced, &golden, "p_utg.lowQ.bed");
    }
}

/// Compares two graphs field by field, reporting the first difference.
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
    assert_eq!(produced.is_sorted, golden.is_sorted, "{}: is_sorted", label);
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

/// Compares two byte buffers and reports the first differing line.
fn assert_lines_equal(produced: &[u8], golden: &[u8], label: &str) {
    if produced == golden {
        return;
    }
    // Leave the artefact behind so the difference can be inspected with `cmp`.
    let artefact = std::env::temp_dir().join(format!("produced.{}", label));
    fs::write(&artefact, produced).unwrap();
    let produced_lines: Vec<&[u8]> = produced.split(|&byte| byte == b'\n').collect();
    let golden_lines: Vec<&[u8]> = golden.split(|&byte| byte == b'\n').collect();
    for (number, (left, right)) in produced_lines.iter().zip(golden_lines.iter()).enumerate() {
        if left != right {
            let offset = left
                .iter()
                .zip(right.iter())
                .position(|(a, b)| a != b)
                .unwrap_or(left.len().min(right.len()));
            let start = offset.saturating_sub(30);
            panic!(
                "{} differs at line {} byte {}:\n  produced ({} B): {:?}\n  golden   ({} B): {:?}",
                label,
                number + 1,
                offset,
                left.len(),
                String::from_utf8_lossy(&left[start..(offset + 30).min(left.len())]),
                right.len(),
                String::from_utf8_lossy(&right[start..(offset + 30).min(right.len())]),
            );
        }
    }
    panic!(
        "{}: the produced output has {} lines, the golden has {}",
        label,
        produced_lines.len(),
        golden_lines.len()
    );
}

#[test]
fn the_hic_clean_stage_reproduces_the_reference_dump() {
    // The probe (`/tmp/hic_probe.cpp`) loaded the cleaning state dumped after
    // `hic_clean_adv` (d20), ran the reference `hic_clean` on it and dumped
    // the result as e01. On this corpus the call is a no-op, but the test
    // still proves the whole unitig-build/pop/delete chain leaves the graph
    // untouched.
    let dump = PathBuf::from("/tmp/cleanprobe");
    let before = dump.join("d20_after_hic_clean.all.debug.asg_t.bin");
    let after = dump.join("e01_after_hic_clean.all.debug.asg_t.bin");
    if !before.exists() || !after.exists() {
        eprintln!("skipping: the hic_clean probe dumps are missing");
        return;
    }

    let mut graph = load_asg(&before);
    hic_clean(&mut graph);
    let golden = load_asg(&after);
    compare_graphs(&graph, &golden, "hic_clean");
}

/// The probe's unitig-list dump: `n`, then for every unitig its `len`,
/// `circ`, `start`, `end`, read count and packed read list.
fn load_unitig_list(path: &Path) -> Vec<MaUtg> {
    let data = fs::read(path).expect("the unitig dump must be readable");
    let mut position = 0usize;
    let count = read_u32(&data, &mut position) as usize;
    let mut unitigs = Vec::with_capacity(count);
    for _ in 0..count {
        let length = read_u32(&data, &mut position) & 0x7fff_ffff;
        let circular = read_u8(&data, &mut position) != 0;
        let start = read_u32(&data, &mut position);
        let end = read_u32(&data, &mut position);
        let read_count = read_u32(&data, &mut position) as usize;
        let mut reads = Vec::with_capacity(read_count);
        for _ in 0..read_count {
            reads.push(read_u64(&data, &mut position));
        }
        unitigs.push(MaUtg {
            length,
            circular,
            start,
            end,
            reads,
            sequence: Vec::new(),
        });
    }
    assert_eq!(position, data.len(), "the unitig dump has trailing bytes");
    unitigs
}

/// Loads the `<tag>.ug.bin` / `<tag>.asg_t.bin` pair the bubble probe wrote.
fn load_unitig_graph(directory: &Path, tag: &str) -> UnitigGraph {
    UnitigGraph {
        unitigs: load_unitig_list(&directory.join(format!("{}.ug.bin", tag))),
        graph: load_asg(&directory.join(format!("{}.asg_t.bin", tag))),
    }
}

/// Compares two unitig lists field by field.
fn compare_unitig_lists(produced: &[MaUtg], golden: &[MaUtg], label: &str) {
    assert_eq!(produced.len(), golden.len(), "{}: unitig count", label);
    for (i, (left, right)) in produced.iter().zip(golden).enumerate() {
        assert_eq!(left.length, right.length, "{}: unitig {} length", label, i);
        assert_eq!(left.circular, right.circular, "{}: unitig {} circ", label, i);
        assert_eq!(left.start, right.start, "{}: unitig {} start", label, i);
        assert_eq!(left.end, right.end, "{}: unitig {} end", label, i);
        assert_eq!(left.reads, right.reads, "{}: unitig {} reads", label, i);
    }
}

/// Compares two unitig graphs.
///
/// `strong` and `no_l_indel` are left out on purpose: `ma_ug_gen_primary`
/// never assigns them, so the reference graph keeps whatever the allocator
/// left in those bits, and no function on the assembly path reads them.
/// The coverage class byte is left out for the same reason: the reference's
/// `asg_seq_set` does not write `c`, so a freshly built unitig graph holds
/// uninitialised class bytes.
fn compare_unitig_graphs(produced: &UnitigGraph, golden: &UnitigGraph, label: &str) {
    compare_unitig_lists(&produced.unitigs, &golden.unitigs, label);
    assert_eq!(
        produced.graph.sequences.len(),
        golden.graph.sequences.len(),
        "{}: vertex count",
        label
    );
    assert_eq!(
        produced.graph.arcs.len(),
        golden.graph.arcs.len(),
        "{}: arc count",
        label
    );
    assert_eq!(produced.graph.is_sorted, golden.graph.is_sorted, "{}: is_sorted", label);
    assert_eq!(
        produced.graph.is_symmetric, golden.graph.is_symmetric,
        "{}: is_symmetric",
        label
    );
    assert_eq!(produced.graph.index, golden.graph.index, "{}: index", label);
    for (i, (left, right)) in produced.graph.sequences.iter().zip(&golden.graph.sequences).enumerate()
    {
        assert_eq!(left.length, right.length, "{}: sequence {} length", label, i);
        assert_eq!(left.deleted, right.deleted, "{}: sequence {} deleted", label, i);
    }
    for (i, (left, right)) in produced.graph.arcs.iter().zip(&golden.graph.arcs).enumerate() {
        assert_eq!(left.ul, right.ul, "{}: arc {} ul", label, i);
        assert_eq!(left.v, right.v, "{}: arc {} v", label, i);
        assert_eq!(left.overlap, right.overlap, "{}: arc {} overlap", label, i);
        assert_eq!(left.deleted, right.deleted, "{}: arc {} deleted", label, i);
        assert_eq!(left.occurrence, right.occurrence, "{}: arc {} occurrence", label, i);
        assert_eq!(left.exact_region, right.exact_region, "{}: arc {} el", label, i);
    }
}

/// The bubble probe states; the directory is not part of the repository.
fn bubble_probe_directory() -> PathBuf {
    PathBuf::from("/tmp/bpprobe")
}

#[test]
fn the_primary_unitig_build_reproduces_the_reference_probe_dump() {
    let dump = bubble_probe_directory();
    let input = dump.join("f00_after_hic_clean.all.debug.asg_t.bin");
    let expected = dump.join("f01_after_ma_ug_primary.ug.bin");
    if !input.exists() || !expected.exists() {
        eprintln!("skipping: the output_bp_graph probe dumps are missing");
        return;
    }

    let mut graph = load_asg(&input);
    let produced = ma_ug_gen_primary(&mut graph, PRIMARY_LABEL);
    let golden = load_unitig_graph(&dump, "f01_after_ma_ug_primary");
    compare_unitig_graphs(&produced, &golden, "ma_ug_gen_primary");
}

#[test]
fn copying_the_graphs_reproduces_the_reference_probe_dumps() {
    let dump = bubble_probe_directory();
    let read_source = dump.join("f00_after_hic_clean.all.debug.asg_t.bin");
    let unitig_source = dump.join("f01_after_ma_ug_primary.ug.bin");
    let read_expected = dump.join("f02_after_copy_sg.all.debug.asg_t.bin");
    if !read_source.exists() || !unitig_source.exists() || !read_expected.exists() {
        eprintln!("skipping: the output_bp_graph probe dumps are missing");
        return;
    }

    let source = load_asg(&read_source);
    let produced = copy_read_graph(&source);
    let golden = load_asg(&read_expected);
    compare_graphs(&produced, &golden, "copy_read_graph");

    let source = load_unitig_graph(&dump, "f01_after_ma_ug_primary");
    let produced = copy_untig_graph(&source);
    let golden = load_unitig_graph(&dump, "f02_after_copy_ug");
    compare_unitig_graphs(&produced, &golden, "copy_untig_graph");
}

#[test]
fn deleting_useless_nodes_reproduces_the_reference_probe_dump() {
    let dump = bubble_probe_directory();
    let source = dump.join("f04_after_clean_primary_untig.ug.bin");
    let expected = dump.join("f05_after_delete_useless_1.ug.bin");
    if !source.exists() || !expected.exists() {
        eprintln!("skipping: the output_bp_graph probe dumps are missing");
        return;
    }

    let mut unitig_graph = load_unitig_graph(&dump, "f04_after_clean_primary_untig");
    delete_useless_nodes(&mut unitig_graph);
    let golden = load_unitig_graph(&dump, "f05_after_delete_useless_1");
    compare_unitig_graphs(&unitig_graph, &golden, "delete_useless_nodes");
}

#[test]
fn renewing_the_unitig_graph_reproduces_the_reference_probe_dump() {
    let dump = bubble_probe_directory();
    let source = dump.join("f05_after_delete_useless_1.ug.bin");
    let expected = dump.join("f06_after_renew_1.ug.bin");
    let read_graph = dump.join("f02_after_copy_sg.all.debug.asg_t.bin");
    if !source.exists() || !expected.exists() || !read_graph.exists() {
        eprintln!("skipping: the output_bp_graph probe dumps are missing");
        return;
    }

    let mut unitig_graph = load_unitig_graph(&dump, "f05_after_delete_useless_1");
    let read_graph = load_asg(&read_graph);
    let edges: Vec<AsgArc> = Vec::new();
    renew_utg(&mut unitig_graph, &read_graph, &edges);
    let golden = load_unitig_graph(&dump, "f06_after_renew_1");
    compare_unitig_graphs(&unitig_graph, &golden, "renew_utg");
}

#[test]
fn advancing_the_unitig_graph_reproduces_the_reference_state() {
    // f02 -> f03: reference `adjust_utg_advance(copy_sg, copy_ug,
    // reverse_sources, ruIndex, &b_mask_t, cov->is_r_het)`.
    //
    // The heterozygosity flags are seeded from the clean_primary probe's h00
    // dump: it is written after `adjust_utg_advance` (which never changes the
    // flags) and was produced from the same `init_hap_cov_t` inputs.
    let dump = bubble_probe_directory();
    let source = dump.join("f02_after_copy_ug.ug.bin");
    let expected = dump.join("f03_after_adjust_utg_advance.ug.bin");
    let read_graph_path = dump.join("f02_after_copy_sg.all.debug.asg_t.bin");
    let reverse_path = dump.join("f02_after_copy_sg.all.debug.reverse.bin");
    let map_path = dump.join("f02_after_copy_sg.all.debug.ruIndex.bin");
    let cov_path = PathBuf::from("/tmp/cleanp/h00_before_clean_primary.cov.bin");
    if !source.exists()
        || !expected.exists()
        || !read_graph_path.exists()
        || !reverse_path.exists()
        || !map_path.exists()
        || !cov_path.exists()
    {
        eprintln!("skipping: the output_bp_graph probe dumps are missing");
        return;
    }

    let mut unitig_graph = load_unitig_graph(&dump, "f02_after_copy_ug");
    let read_graph = load_asg(&read_graph_path);
    // The `write_debug_graph` dumps use `write_debug_ma_hit_ts`, whose record
    // layout is the one `load_sources` parses.
    let reverse_sources = load_sources(&reverse_path);
    let (mut read_to_unitig, _trio_labels) = load_read_to_unitig(&map_path);

    let (cov_values, is_r_het) = load_hap_cov_fields(&cov_path);
    assert_eq!(cov_values.len(), read_graph.sequences.len());
    assert_eq!(is_r_het.len(), read_graph.sequences.len());

    // `init_bub_label_t(&b_mask_t, min(threads, 10), sg->n_seq)`: the reference
    // run uses one thread.
    let mut labels = BubbleLabel::new(1, read_graph.sequences.len() as u32);
    adjust_utg_advance(
        &read_graph,
        &mut unitig_graph,
        &reverse_sources,
        &mut read_to_unitig,
        Some(&mut labels),
        &is_r_het,
        PURGE_SIMILARITY_THRESHOLD,
    );

    let golden = load_unitig_graph(&dump, "f03_after_adjust_utg_advance");
    compare_unitig_graphs(&unitig_graph, &golden, "adjust_utg_advance");

    let (golden_map, _golden_labels) = load_read_to_unitig(
        &dump.join("f03_after_adjust_utg_advance_sg.all.debug.ruIndex.bin"),
    );
    assert_eq!(
        read_to_unitig.entries(),
        golden_map.entries(),
        "adjust_utg_advance: read map"
    );

    let golden_read_graph =
        load_asg(&dump.join("f03_after_adjust_utg_advance_sg.all.debug.asg_t.bin"));
    compare_graphs(&read_graph, &golden_read_graph, "adjust_utg_advance: read graph");
}

/// `discard_small_ctg` runs on the state after the contained-read rescue. The
/// reference call passes `asm_opt.max_contig_tip` (3) with `is_renew0 = 0` and
/// `is_renew1 = 1`, so the graph is only regenerated when something was cut.
#[test]
fn discarding_small_contigs_reproduces_the_reference_state() {
    let dump = bubble_probe_directory();
    let before = dump.join("f11_after_rescue_contained_renew.ug.bin");
    let after = dump.join("f12_after_discard_small_ctg.ug.bin");
    let read_graph_path = dump.join("f10_after_rescue_contained_sg.all.debug.asg_t.bin");
    if !before.exists() || !after.exists() || !read_graph_path.exists() {
        eprintln!("skipping: the output_bp_graph probe dumps are missing");
        return;
    }

    let mut unitig_graph = load_unitig_graph(&dump, "f11_after_rescue_contained_renew");
    let read_graph = load_asg(&read_graph_path);

    // The probe's `new_rtg_edges` is empty: the rescues ran with a NULL edge
    // vector, so nothing was recorded for the later renewal.
    discard_small_ctg(&mut unitig_graph, &read_graph, &[], 3, false, true);

    let golden = load_unitig_graph(&dump, "f12_after_discard_small_ctg");
    compare_unitig_graphs(&unitig_graph, &golden, "discard_small_ctg");
}

/// `recover_utg_by_coverage` runs after the primary label propagation dropped
/// every alternate read. On this corpus the read graph keeps no alternate
/// sequence, so `ma_ug_gen_primary(read_g, ALTER_LABLE)` starts empty and the
/// unitig graph must come back unchanged.
#[test]
fn recovering_alternate_unitigs_reproduces_the_reference_state() {
    let dump = bubble_probe_directory();
    let before = dump.join("f12_after_discard_small_ctg.ug.bin");
    let after = dump.join("f14_after_recover_by_coverage.ug.bin");
    let read_graph_path = dump.join("f13_after_label_and_drop.all.debug.asg_t.bin");
    let coverage_path = dump.join("f13_after_label_and_drop.all.debug.coverage_cut.bin");
    let source_path = dump.join("f13_after_label_and_drop.all.debug.source.bin");
    let map_path = dump.join("f13_after_label_and_drop.all.debug.ruIndex.bin");
    if !before.exists()
        || !after.exists()
        || !read_graph_path.exists()
        || !coverage_path.exists()
        || !source_path.exists()
        || !map_path.exists()
    {
        eprintln!("skipping: the output_bp_graph probe dumps are missing");
        return;
    }

    let mut unitig_graph = load_unitig_graph(&dump, "f12_after_discard_small_ctg");
    let mut read_graph = load_asg(&read_graph_path);
    let coverage_cut = load_coverage_cut(&coverage_path);
    let sources = load_sources(&source_path);
    let (read_to_unitig, _trio_labels) = load_read_to_unitig(&map_path);

    // `hom_global_coverage` is 45 with the "set" flag cleared, so the `-1024`
    // default becomes `[45/1.25*0.85, INT32_MAX]`: the integer stores
    // truncate to 30 and the reference logs "[30, infinity]".
    recover_utg_by_coverage(
        &mut unitig_graph,
        &mut read_graph,
        &coverage_cut,
        &sources,
        &read_to_unitig,
        30,
        i32::MAX,
        // `cov->t_ch->ir_het` only matters when a unitig is appended, which
        // cannot happen while the alternate read graph is empty.
        None,
    );

    let golden = load_unitig_graph(&dump, "f14_after_recover_by_coverage");
    compare_unitig_graphs(&unitig_graph, &golden, "recover_utg_by_coverage");
}

/// The clean_primary probe's `<tag>.cov.bin`: `n`, the per-read coverage and
/// the per-read heterozygosity flags.
fn load_hap_cov_fields(path: &Path) -> (Vec<u32>, Vec<u8>) {
    let data = fs::read(path).expect("the coverage dump must be readable");
    let mut position = 0usize;
    let count = read_u32(&data, &mut position) as usize;
    let mut cov = Vec::with_capacity(count);
    for _ in 0..count {
        cov.push(read_u32(&data, &mut position));
    }
    let mut is_r_het = vec![0u8; count];
    is_r_het.copy_from_slice(&data[position..position + count]);
    position += count;
    assert_eq!(position, data.len(), "the coverage dump has trailing bytes");
    (cov, is_r_het)
}

#[test]
fn popping_primary_bubbles_reproduces_the_reference_state() {
    // The clean_primary probe (`/tmp/cleanp_probe.cpp`) dumps the state right
    // before `clean_primary_untig_graph` as h00, then runs the reference
    // `asg_pop_bubble_primary_trio` with `i_max_dist = NULL`,
    // `positive_flag = -1`, `negative_flag = DROP` and
    // `is_update_chain = 1`, and dumps the result as h01.
    //
    // The seeded coverage comes from the same probe: `cov.cov` and the
    // `C_HET` flags were computed by the reference `init_hap_cov_t`, and the
    // read lengths are those of the read graph the probe loaded (`d20`).
    let dump = PathBuf::from("/tmp/cleanp");
    let before = dump.join("h00_before_clean_primary.ug.bin");
    let after = dump.join("h01_after_pop_bubble.ug.bin");
    let read_dump = PathBuf::from("/tmp/cleanprobe/d20_after_hic_clean.all.debug.asg_t.bin");
    if !before.exists() || !after.exists() || !read_dump.exists() {
        eprintln!("skipping: the clean_primary probe dumps are missing");
        return;
    }

    let mut unitig_graph = load_unitig_graph(&dump, "h00_before_clean_primary");
    let (cov_values, is_r_het) =
        load_hap_cov_fields(&dump.join("h00_before_clean_primary.cov.bin"));
    let read_graph = load_asg(&read_dump);
    let read_lengths: Vec<u32> = read_graph.sequences.iter().map(|seq| seq.length).collect();
    assert_eq!(
        cov_values.len(),
        read_lengths.len(),
        "the coverage and read-graph dumps cover different read sets"
    );
    let read_count = read_lengths.len();
    let mut cov = HapCov {
        n: read_count as u32,
        cov: cov_values,
        pos_idx: vec![u64::MAX; read_count],
        read_lengths,
        max_hang: MAXIMUM_HANG,
        min_ovlp: MINIMUM_OVERLAP,
        is_r_het,
        t_ch: Some(init_trans_chain(&unitig_graph, &read_graph, read_count)),
        ..HapCov::default()
    };

    // `R_INF.trio_flag` is filled with `AMBIGU` (0) for a HiFi-only run.
    let trio_labels = vec![0u8; read_count];
    asg_pop_bubble_primary_trio(
        &mut unitig_graph,
        None,
        u32::MAX,
        DROP,
        Some(&mut cov),
        &trio_labels,
        true,
        None,
    );

    let golden = load_unitig_graph(&dump, "h01_after_pop_bubble");
    compare_unitig_graphs(&unitig_graph, &golden, "asg_pop_bubble_primary_trio");

    let (golden_cov, golden_is_r_het) =
        load_hap_cov_fields(&dump.join("h01_after_pop_bubble.cov.bin"));
    assert_eq!(
        cov.cov, golden_cov,
        "asg_pop_bubble_primary_trio: read coverage"
    );
    assert_eq!(
        cov.is_r_het, golden_is_r_het,
        "asg_pop_bubble_primary_trio: read heterozygosity flags"
    );
}

/// The clean_primary probe's `<tag>.ruindex.bin`: the raw `R_to_U` array.
fn load_raw_read_to_unitig(path: &Path) -> ReadToUnitig {
    let data = fs::read(path).expect("the read map dump must be readable");
    let mut position = 0usize;
    let length = read_u64(&data, &mut position) as usize;
    let mut entries = Vec::with_capacity(length);
    for _ in 0..length {
        entries.push(read_u32(&data, &mut position));
    }
    assert_eq!(position, data.len(), "the read map dump has trailing bytes");
    ReadToUnitig::from_entries(entries)
}

/// The clean_primary probe's `<tag>.reverse.bin`: `R_INF.reverse_paf`, one
/// hit field at a time.
fn load_reverse_sources(path: &Path) -> Vec<PafList> {
    let data = fs::read(path).expect("the reverse overlap dump must be readable");
    let mut position = 0usize;
    let read_count = read_u32(&data, &mut position) as usize;
    let mut lists = Vec::with_capacity(read_count);
    for _ in 0..read_count {
        let hit_count = read_u32(&data, &mut position) as usize;
        let mut hits = Vec::with_capacity(hit_count);
        for _ in 0..hit_count {
            let qns = read_u64(&data, &mut position);
            let query_end = read_u32(&data, &mut position);
            let target_identifier = read_u32(&data, &mut position);
            let target_start = read_u32(&data, &mut position);
            let target_end = read_u32(&data, &mut position);
            let field_word = read_u32(&data, &mut position);
            let block_word = read_u32(&data, &mut position);
            let exact_region = read_u8(&data, &mut position) != 0;
            let without_large_indel = read_u8(&data, &mut position) != 0;
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
                exact_region,
                deleted: block_word >> 31 != 0,
            });
        }
        lists.push(PafList {
            hits,
            is_fully_corrected: false,
            is_abnormal: 0,
        });
    }
    assert_eq!(
        position,
        data.len(),
        "the reverse overlap dump has trailing bytes"
    );
    lists
}

#[test]
fn cutting_trio_tips_reproduces_the_reference_state() {
    // The clean_primary probe (`/tmp/cleanp_probe2.cpp`) dumps h01 after
    // `asg_pop_bubble_primary_trio`, then runs the reference
    // `cut_trio_tip_primary(g, ug, max_short_tip * 2, -1, 0, read_g,
    // reverse_sources, ruIndex, is_r_het, 2)` and dumps h02.
    //
    // `purge_simi_thres` is 0.55 because the probe's `CommandLine_process`
    // sees the default `-l3`.
    let dump = PathBuf::from("/tmp/cleanp");
    let before = dump.join("h01_after_pop_bubble.ug.bin");
    let after = dump.join("h02_after_cut_trio_tip.ug.bin");
    let read_dump = PathBuf::from("/tmp/cleanprobe/d20_after_hic_clean.all.debug.asg_t.bin");
    if !before.exists() || !after.exists() || !read_dump.exists() {
        eprintln!("skipping: the clean_primary probe dumps are missing");
        return;
    }

    let UnitigGraph { unitigs, mut graph } = load_unitig_graph(&dump, "h01_after_pop_bubble");
    let (_, is_r_het) = load_hap_cov_fields(&dump.join("h01_after_pop_bubble.cov.bin"));
    let read_graph = load_asg(&read_dump);
    let reverse_sources = load_reverse_sources(&dump.join("h01_after_pop_bubble.reverse.bin"));
    let mut read_to_unitig =
        load_raw_read_to_unitig(&dump.join("h01_after_pop_bubble.ruindex.bin"));
    // `R_INF.trio_flag` is filled with `AMBIGU` (0) for a HiFi-only run.
    let trio_labels = vec![0u8; read_graph.sequences.len()];

    let mut options = AssemblyOptions::default();
    options.purge_similarity_threshold = options.purge_similarity_rate_level_three;
    let max_ext = (options.maximum_short_tip_read_count * 2) as u32;

    let count = cut_trio_tip_primary(
        &mut graph,
        Some(&unitigs),
        max_ext,
        u32::MAX,
        0,
        &read_graph,
        &reverse_sources,
        &mut read_to_unitig,
        &is_r_het,
        2,
        &trio_labels,
        options.purge_similarity_threshold,
    );
    assert!(count > 0, "cut_trio_tip_primary: no tip was processed");

    let produced = UnitigGraph { unitigs, graph };
    let golden = load_unitig_graph(&dump, "h02_after_cut_trio_tip");
    compare_unitig_graphs(&produced, &golden, "cut_trio_tip_primary");
    for (i, (left, right)) in produced
        .graph
        .sequences
        .iter()
        .zip(&golden.graph.sequences)
        .enumerate()
    {
        assert_eq!(
            left.coverage_class, right.coverage_class,
            "cut_trio_tip_primary: sequence {} class",
            i
        );
    }

    let golden_map = load_raw_read_to_unitig(&dump.join("h02_after_cut_trio_tip.ruindex.bin"));
    assert_eq!(
        read_to_unitig.len(),
        golden_map.len(),
        "cut_trio_tip_primary: read map length"
    );
    assert_eq!(
        read_to_unitig.entries(),
        golden_map.entries(),
        "cut_trio_tip_primary: read map"
    );
}

/// The equal-tips cut of the primary cleaning, validated against the state
/// dumped after the reference `asg_arc_cut_trio_long_equal_tips_assembly`.
///
/// The probe replays `clean_primary_untig_graph` from h00 and dumps the state
/// as `h04_p0_i0_tip` before the call and `h05_p0_i0_eq` after it. The call
/// updates the unitig graph's class labels and the haplotype coverage
/// (`cov`/`is_r_het`); the read map comes back unchanged.
#[test]
fn cutting_equal_tips_reproduces_the_reference_state() {
    let dump = PathBuf::from("/tmp/cleanp");
    let before = dump.join("h04_p0_i0_tip.ug.bin");
    let after = dump.join("h05_p0_i0_eq.ug.bin");
    let read_dump = PathBuf::from("/tmp/cleanprobe/d20_after_hic_clean.all.debug.asg_t.bin");
    let coverage_dump =
        PathBuf::from("/tmp/cleanprobe/d20_after_hic_clean.all.debug.coverage_cut.bin");
    if !before.exists() || !after.exists() || !read_dump.exists() || !coverage_dump.exists() {
        eprintln!("skipping: the clean_primary probe dumps are missing");
        return;
    }

    let UnitigGraph { unitigs, mut graph } = load_unitig_graph(&dump, "h04_p0_i0_tip");
    let (cov_values, is_r_het) =
        load_hap_cov_fields(&dump.join("h04_p0_i0_tip.cov.bin"));
    let read_graph = load_asg(&read_dump);
    let coverage_cut = load_coverage_cut(&coverage_dump);
    let reverse_sources = load_reverse_sources(&dump.join("h04_p0_i0_tip.reverse.bin"));
    let mut read_to_unitig = load_raw_read_to_unitig(&dump.join("h04_p0_i0_tip.ruindex.bin"));

    let read_lengths: Vec<u32> = read_graph.sequences.iter().map(|seq| seq.length).collect();
    let read_count = read_lengths.len();
    assert_eq!(
        cov_values.len(),
        read_count,
        "the coverage and read-graph dumps cover different read sets"
    );
    let mut cov = HapCov {
        n: read_count as u32,
        cov: cov_values,
        pos_idx: vec![u64::MAX; read_count],
        read_lengths,
        max_hang: MAXIMUM_HANG,
        min_ovlp: MINIMUM_OVERLAP,
        is_r_het,
        t_ch: Some(init_trans_chain(&UnitigGraph {
            unitigs: unitigs.clone(),
            graph: graph.clone(),
        }, &read_graph, read_count)),
        ..HapCov::default()
    };

    // `R_INF.trio_flag` is filled with `AMBIGU` for a HiFi-only run.
    let trio_labels = vec![0u8; read_count];
    let mut options = AssemblyOptions::default();
    options.purge_similarity_threshold = options.purge_similarity_rate_level_three;

    let count = asg_arc_cut_trio_long_equal_tips_assembly(
        &mut graph,
        &unitigs,
        &read_graph,
        &reverse_sources,
        2,
        &mut read_to_unitig,
        u32::MAX,
        &trio_labels,
        &coverage_cut,
        options.purge_similarity_threshold,
        MAX_HANG_RATE,
        &mut cov,
    );
    assert!(count > 0, "asg_arc_cut_trio_long_equal_tips_assembly: no tip was cut");

    let produced = UnitigGraph { unitigs, graph };
    let golden = load_unitig_graph(&dump, "h05_p0_i0_eq");
    compare_unitig_graphs(&produced, &golden, "asg_arc_cut_trio_long_equal_tips_assembly");
    for (i, (left, right)) in produced
        .graph
        .sequences
        .iter()
        .zip(&golden.graph.sequences)
        .enumerate()
    {
        assert_eq!(
            left.coverage_class, right.coverage_class,
            "asg_arc_cut_trio_long_equal_tips_assembly: sequence {} class",
            i
        );
    }

    let (golden_cov, golden_is_r_het) = load_hap_cov_fields(&dump.join("h05_p0_i0_eq.cov.bin"));
    assert_eq!(
        cov.cov, golden_cov,
        "asg_arc_cut_trio_long_equal_tips_assembly: read coverage"
    );
    assert_eq!(
        cov.is_r_het, golden_is_r_het,
        "asg_arc_cut_trio_long_equal_tips_assembly: heterozygosity flags"
    );

    let golden_map = load_raw_read_to_unitig(&dump.join("h05_p0_i0_eq.ruindex.bin"));
    assert_eq!(
        read_to_unitig.entries(),
        golden_map.entries(),
        "asg_arc_cut_trio_long_equal_tips_assembly: read map"
    );
}

/// Every state the clean_primary probe records for one cleaning step: the
/// unitig graph, the seeded coverage and the read map.
struct CleanProbeState {
    unitig_graph: UnitigGraph,
    cov: HapCov,
    read_graph: Asg,
    coverage_cut: Vec<CoverageCut>,
    reverse_sources: Vec<PafList>,
    read_to_unitig: ReadToUnitig,
    trio_labels: Vec<u8>,
}

fn clean_probe_dump() -> PathBuf {
    PathBuf::from("/tmp/cleanp")
}

fn clean_probe_read_dump() -> PathBuf {
    PathBuf::from("/tmp/cleanprobe/d20_after_hic_clean.all.debug.asg_t.bin")
}

fn clean_probe_coverage_dump() -> PathBuf {
    PathBuf::from("/tmp/cleanprobe/d20_after_hic_clean.all.debug.coverage_cut.bin")
}

/// Loads the state dumped as `tag` by `/tmp/cleanp_probe2.cpp`.
///
/// The seeded coverage comes from the probe's own `cov.bin`: `cov.cov` and the
/// `C_HET`/`P_HET` flags were computed by the reference `init_hap_cov_t` and
/// updated by the reference steps that preceded `tag`, so the test does not
/// have to replay them.
fn load_clean_probe_state(tag: &str) -> Option<CleanProbeState> {
    let dump = clean_probe_dump();
    let cov_path = dump.join(format!("{}.cov.bin", tag));
    let map_path = dump.join(format!("{}.ruindex.bin", tag));
    let reverse_path = dump.join(format!("{}.reverse.bin", tag));
    let read_dump = clean_probe_read_dump();
    let coverage_dump = clean_probe_coverage_dump();
    if !cov_path.exists()
        || !map_path.exists()
        || !reverse_path.exists()
        || !read_dump.exists()
        || !coverage_dump.exists()
    {
        return None;
    }

    let unitig_graph = load_unitig_graph(&dump, tag);
    let (cov_values, is_r_het) = load_hap_cov_fields(&cov_path);
    let read_graph = load_asg(&read_dump);
    let read_lengths: Vec<u32> = read_graph.sequences.iter().map(|seq| seq.length).collect();
    let read_count = read_lengths.len();
    assert_eq!(
        cov_values.len(),
        read_count,
        "the coverage and read-graph dumps cover different read sets"
    );
    let cov = HapCov {
        n: read_count as u32,
        cov: cov_values,
        pos_idx: vec![u64::MAX; read_count],
        read_lengths,
        max_hang: MAXIMUM_HANG,
        min_ovlp: MINIMUM_OVERLAP,
        is_r_het,
        t_ch: Some(init_trans_chain(&unitig_graph, &read_graph, read_count)),
        ..HapCov::default()
    };

    Some(CleanProbeState {
        unitig_graph,
        cov,
        read_graph,
        coverage_cut: load_coverage_cut(&coverage_dump),
        reverse_sources: load_reverse_sources(&reverse_path),
        read_to_unitig: load_raw_read_to_unitig(&map_path),
        // `R_INF.trio_flag` is filled with `AMBIGU` (0) for a HiFi-only run.
        trio_labels: vec![0u8; read_count],
    })
}

/// Compares one cleaning result against the probe state dumped as `tag`.
fn compare_clean_probe_state(
    unitig_graph: &UnitigGraph,
    cov: &HapCov,
    read_to_unitig: &ReadToUnitig,
    tag: &str,
    label: &str,
) {
    let dump = clean_probe_dump();
    let golden = load_unitig_graph(&dump, tag);
    compare_unitig_graphs(unitig_graph, &golden, label);
    for (i, (left, right)) in unitig_graph
        .graph
        .sequences
        .iter()
        .zip(&golden.graph.sequences)
        .enumerate()
    {
        assert_eq!(
            left.coverage_class, right.coverage_class,
            "{}: sequence {} class",
            label, i
        );
    }

    let (golden_cov, golden_is_r_het) = load_hap_cov_fields(&dump.join(format!("{}.cov.bin", tag)));
    assert_eq!(cov.cov, golden_cov, "{}: read coverage", label);
    assert_eq!(cov.is_r_het, golden_is_r_het, "{}: heterozygosity flags", label);

    let golden_map = load_raw_read_to_unitig(&dump.join(format!("{}.ruindex.bin", tag)));
    assert_eq!(read_to_unitig.len(), golden_map.len(), "{}: read map length", label);
    assert_eq!(read_to_unitig.entries(), golden_map.entries(), "{}: read map", label);
}

/// `purge_simi_thres` as `CommandLine_process` sets it for the default `-l3`.
const PURGE_SIMILARITY_THRESHOLD: f32 = 0.55;

#[test]
fn cutting_long_tips_reproduces_the_reference_state() {
    // h03 -> h04: reference `asg_arc_cut_trio_long_tip_primary(g, ug, read_g,
    // reverse_sources, ruIndex, 2, tip_drop_ratio, -1, cov, NULL)`.
    let Some(state) = load_clean_probe_state("h03_p0_i0_pop") else {
        eprintln!("skipping: the clean_primary probe dumps are missing");
        return;
    };
    let CleanProbeState {
        mut unitig_graph,
        mut cov,
        read_graph,
        coverage_cut,
        reverse_sources,
        mut read_to_unitig,
        trio_labels,
    } = state;
    let UnitigGraph { unitigs, graph } = &mut unitig_graph;

    let _count = asg_arc_cut_trio_long_tip_primary(
        graph,
        unitigs,
        &read_graph,
        &reverse_sources,
        &mut read_to_unitig,
        2,
        0.15,
        u32::MAX,
        &mut cov,
        &trio_labels,
        PURGE_SIMILARITY_THRESHOLD,
        MAX_HANG_RATE,
        &coverage_cut,
    );

    compare_clean_probe_state(
        &unitig_graph,
        &cov,
        &read_to_unitig,
        "h04_p0_i0_tip",
        "asg_arc_cut_trio_long_tip_primary",
    );
}

#[test]
fn cutting_long_tip_complexes_reproduces_the_reference_state() {
    // h05 -> h06: reference `asg_arc_cut_trio_long_tip_primary_complex(g, ug,
    // read_g, reverse_sources, ruIndex, 2, tip_drop_ratio, stops_threshold,
    // cov, NULL, -1)`.
    let Some(state) = load_clean_probe_state("h05_p0_i0_eq") else {
        eprintln!("skipping: the clean_primary probe dumps are missing");
        return;
    };
    let CleanProbeState {
        mut unitig_graph,
        mut cov,
        read_graph,
        coverage_cut,
        reverse_sources,
        mut read_to_unitig,
        trio_labels,
    } = state;
    let UnitigGraph { unitigs, graph } = &mut unitig_graph;

    let _count = asg_arc_cut_trio_long_tip_primary_complex(
        graph,
        unitigs,
        &read_graph,
        &reverse_sources,
        &mut read_to_unitig,
        2,
        0.15,
        3,
        u32::MAX,
        &mut cov,
        &trio_labels,
        PURGE_SIMILARITY_THRESHOLD,
        MAX_HANG_RATE,
        &coverage_cut,
    );

    compare_clean_probe_state(
        &unitig_graph,
        &cov,
        &read_to_unitig,
        "h06_p0_i0_tipcx",
        "asg_arc_cut_trio_long_tip_primary_complex",
    );
}

#[test]
fn cutting_equal_tip_complexes_reproduces_the_reference_state() {
    // h06 -> h07: reference `asg_arc_cut_trio_long_equal_tips_assembly_complex(
    // g, ug, read_g, reverse_sources, 2, ruIndex, stops_threshold, cov, NULL,
    // -1)`.
    let Some(state) = load_clean_probe_state("h06_p0_i0_tipcx") else {
        eprintln!("skipping: the clean_primary probe dumps are missing");
        return;
    };
    let CleanProbeState {
        mut unitig_graph,
        mut cov,
        read_graph,
        coverage_cut,
        reverse_sources,
        mut read_to_unitig,
        trio_labels,
    } = state;
    let UnitigGraph { unitigs, graph } = &mut unitig_graph;

    let _count = asg_arc_cut_trio_long_equal_tips_assembly_complex(
        graph,
        unitigs,
        &read_graph,
        &reverse_sources,
        2,
        &mut read_to_unitig,
        3,
        u32::MAX,
        &mut cov,
        &trio_labels,
        PURGE_SIMILARITY_THRESHOLD,
        MAX_HANG_RATE,
        &coverage_cut,
    );

    compare_clean_probe_state(
        &unitig_graph,
        &cov,
        &read_to_unitig,
        "h07_p0_i0_eqcx",
        "asg_arc_cut_trio_long_equal_tips_assembly_complex",
    );
}

#[test]
fn detecting_chimeric_topologies_reproduces_the_reference_state() {
    // h07 -> h08: reference `detect_chimeric_by_topo(g, ug, read_g,
    // reverse_sources, 2, stops_threshold, chimeric_rate, ruIndex, NULL,
    // cov->is_r_het)`.
    let Some(state) = load_clean_probe_state("h07_p0_i0_eqcx") else {
        eprintln!("skipping: the clean_primary probe dumps are missing");
        return;
    };
    let CleanProbeState {
        mut unitig_graph,
        cov,
        read_graph,
        reverse_sources,
        mut read_to_unitig,
        ..
    } = state;
    let UnitigGraph { unitigs, graph } = &mut unitig_graph;

    let _count = detect_chimeric_by_topo(
        graph,
        unitigs,
        &read_graph,
        &reverse_sources,
        2,
        3,
        0.05,
        &mut read_to_unitig,
        &cov.is_r_het,
        PURGE_SIMILARITY_THRESHOLD,
    );

    compare_clean_probe_state(
        &unitig_graph,
        &cov,
        &read_to_unitig,
        "h08_p0_i0_chim",
        "detect_chimeric_by_topo",
    );
}

#[test]
fn cutting_late_trio_tips_reproduces_the_reference_state() {
    // h20 -> h21: the pass-final `cut_trio_tip_primary`.
    let Some(state) = load_clean_probe_state("h20_p0_i2_chim") else {
        eprintln!("skipping: the clean_primary probe dumps are missing");
        return;
    };
    let CleanProbeState {
        mut unitig_graph,
        cov,
        read_graph,
        reverse_sources,
        mut read_to_unitig,
        trio_labels,
        ..
    } = state;
    let UnitigGraph { unitigs, graph } = &mut unitig_graph;

    let _count = cut_trio_tip_primary(
        graph,
        Some(unitigs),
        6,
        u32::MAX,
        0,
        &read_graph,
        &reverse_sources,
        &mut read_to_unitig,
        &cov.is_r_het,
        2,
        &trio_labels,
        PURGE_SIMILARITY_THRESHOLD,
    );

    compare_clean_probe_state(
        &unitig_graph,
        &cov,
        &read_to_unitig,
        "h21_late_cut_trio_tip",
        "cut_trio_tip_primary (late)",
    );
}

#[test]
fn resolving_late_tangles_reproduces_the_reference_state() {
    // h21 -> h22: reference `resolve_tangles(ug, read_g, reverse_sources, 20,
    // 100, 0.05, 0.2, ruIndex, cov->is_r_het, -1, drop_ratio)`.
    let Some(state) = load_clean_probe_state("h21_late_cut_trio_tip") else {
        eprintln!("skipping: the clean_primary probe dumps are missing");
        return;
    };
    let CleanProbeState {
        mut unitig_graph,
        cov,
        read_graph,
        reverse_sources,
        mut read_to_unitig,
        trio_labels,
        ..
    } = state;

    resolve_tangles(
        &mut unitig_graph,
        &read_graph,
        &reverse_sources,
        20,
        100,
        0.05,
        0.2,
        &mut read_to_unitig,
        &cov.is_r_het,
        u32::MAX,
        0.9,
        &trio_labels,
        PURGE_SIMILARITY_THRESHOLD,
    );

    compare_clean_probe_state(
        &unitig_graph,
        &cov,
        &read_to_unitig,
        "h22_late_resolve_tangles",
        "resolve_tangles",
    );
}

#[test]
fn dropping_semi_circles_reproduces_the_reference_state() {
    // h22 -> h23: reference `drop_semi_circle(ug, g, read_g,
    // reverse_sources, ruIndex, cov->is_r_het)`.
    let Some(state) = load_clean_probe_state("h22_late_resolve_tangles") else {
        eprintln!("skipping: the clean_primary probe dumps are missing");
        return;
    };
    let CleanProbeState {
        mut unitig_graph,
        cov,
        read_graph,
        reverse_sources,
        mut read_to_unitig,
        ..
    } = state;

    drop_semi_circle(
        &mut unitig_graph,
        &read_graph,
        &reverse_sources,
        &mut read_to_unitig,
        &cov.is_r_het,
        PURGE_SIMILARITY_THRESHOLD,
    );

    compare_clean_probe_state(
        &unitig_graph,
        &cov,
        &read_to_unitig,
        "h23_late_drop_semi_circle",
        "drop_semi_circle",
    );
}

#[test]
fn dropping_short_diploid_topologies_reproduces_the_reference_state() {
    // h23 -> h24: reference `unitig_arc_del_short_diploid_by_length_topo(g, ug,
    // drop_ratio, max_short_tip, reverse_sources, 0, 1)`.
    let Some(state) = load_clean_probe_state("h23_late_drop_semi_circle") else {
        eprintln!("skipping: the clean_primary probe dumps are missing");
        return;
    };
    let CleanProbeState {
        mut unitig_graph,
        cov,
        read_to_unitig,
        ..
    } = state;
    let UnitigGraph { unitigs, graph } = &mut unitig_graph;

    let _count = unitig_arc_del_short_diploid_by_length_topo(graph, unitigs, 0.9, 3, false, true);

    compare_clean_probe_state(
        &unitig_graph,
        &cov,
        &read_to_unitig,
        "h24_late_short_diploid_topo",
        "unitig_arc_del_short_diploid_by_length_topo",
    );
}

#[test]
fn dropping_short_diploids_by_length_reproduces_the_reference_state() {
    // h37 -> h38: the pass-final `unitig_arc_del_short_diploid_by_length`.
    let Some(state) = load_clean_probe_state("h37_late_short_diploid_topo") else {
        eprintln!("skipping: the clean_primary probe dumps are missing");
        return;
    };
    let CleanProbeState {
        mut unitig_graph,
        cov,
        read_to_unitig,
        ..
    } = state;

    let _count = unitig_arc_del_short_diploid_by_length(&mut unitig_graph.graph, 0.9);

    compare_clean_probe_state(
        &unitig_graph,
        &cov,
        &read_to_unitig,
        "h38_pass_short_diploid",
        "unitig_arc_del_short_diploid_by_length",
    );
}

#[test]
fn cleaning_the_primary_unitig_graph_reproduces_the_reference_state() {
    // h00 -> h51: the whole `clean_primary_untig_graph` replay, with the
    // parameters `output_contig_graph_primary` passes for a HiFi-only run.
    let Some(mut state) = load_clean_probe_state("h00_before_clean_primary") else {
        eprintln!("skipping: the clean_primary probe dumps are missing");
        return;
    };

    clean_primary_untig_graph(
        &mut state.unitig_graph,
        &state.read_graph,
        &state.reverse_sources,
        &state.coverage_cut,
        6,
        0.15,
        3,
        &mut state.read_to_unitig,
        0.05,
        false,
        0.9,
        &mut state.cov,
        &state.trio_labels,
        3,
        PURGE_SIMILARITY_THRESHOLD,
        MAX_HANG_RATE,
    );

    compare_clean_probe_state(
        &state.unitig_graph,
        &state.cov,
        &state.read_to_unitig,
        "h51_late_short_diploid_topo",
        "clean_primary_untig_graph",
    );
}


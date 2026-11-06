//! Byte-for-byte validation of `output_bp_graph`'s primary-contig half.
//!
//! The reference probe `/tmp/bp_probe.cpp` replays `adjust_utg_by_primary` and
//! `print_utg` from the state dumped after `hic_clean_adv`
//! (`/tmp/cleanprobe/d20_after_hic_clean.all.debug.*`) and writes
//! `probe.bp.p_ctg.*`, which match the reference run's `bp.p_ctg.*` byte for
//! byte. This test replays the same state through the Rust port.
//!
//! The dumps live in a scratch directory and are not part of the repository,
//! so the test reports a skip when they are absent.

use std::fs;
use std::path::{Path, PathBuf};

use rshifiasm::assembly::asg::Asg;
use rshifiasm::assembly::bp_graph::{adjust_utg_by_primary, output_bp_graph};
use rshifiasm::assembly::bubble::BubbleLabel;
use rshifiasm::assembly::hap_cov::UTrans;
use rshifiasm::assembly::hamming::reduce_hamming_error_adv;
use rshifiasm::assembly::hit::PRIMARY_LABEL;
use rshifiasm::assembly::probe::{
    load_asg, load_coverage_cut, load_read_to_unitig, load_sources,
};
use rshifiasm::assembly::purge_dups::HOM_PEAK_RATE;
use rshifiasm::assembly::trans_clean::{clean_u_trans_t_idx, set_trio_flag_by_cov};
use rshifiasm::assembly::unitig::{
    copy_read_graph, copy_untig_graph, hic_clean, ma_ug_gen_phase, ma_ug_gen_primary,
    PHASE_SEF, PHASE_SEP, PHASE_SEP_RATE,
};
use rshifiasm::command_line::options::AssemblyOptions;
use rshifiasm::reads::ReadStore;
use rshifiasm::sequence::encoding::PackedSequence;

fn repository_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn probe_dump_directory() -> PathBuf {
    std::env::var_os("RSHIFIASM_BP_DUMP_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("/tmp/cleanprobe"))
}

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
        store.sequences.push(PackedSequence::from_ascii(sequence));
        store.trio_labels.push(0);
        store.total_reads += 1;
        store.total_bases += sequence.len() as u64;
    }
    store
}

fn assert_same_file(produced: &Path, golden: &Path, label: &str) {
    let produced_bytes = fs::read(produced)
        .unwrap_or_else(|error| panic!("cannot read {}: {}", produced.display(), error));
    let golden_bytes = fs::read(golden).expect("the golden must be readable");
    assert_eq!(
        produced_bytes.len(),
        golden_bytes.len(),
        "{}: length differs ({} vs {})",
        label,
        produced_bytes.len(),
        golden_bytes.len()
    );
    assert!(produced_bytes == golden_bytes, "{} differs from the golden", label);
}

#[test]
fn the_primary_contig_output_reproduces_the_reference() {
    let dump = probe_dump_directory();
    let config = "d20_after_hic_clean";
    let asg_path = dump.join(format!("{}.all.debug.asg_t.bin", config));
    let source_path = dump.join(format!("{}.all.debug.source.bin", config));
    let reverse_path = dump.join(format!("{}.all.debug.reverse.bin", config));
    let coverage_path = dump.join(format!("{}.all.debug.coverage_cut.bin", config));
    let map_path = dump.join(format!("{}.all.debug.ruIndex.bin", config));
    let corrected_reads = repository_root().join("reference/write_paf_ec.ec.fa");
    let golden_prefix = repository_root().join("reference/default_repeat.bp");

    if !asg_path.exists() || !corrected_reads.exists() {
        eprintln!(
            "skipping: the reference dump ({}) or the corrected reads are missing",
            asg_path.display()
        );
        return;
    }

    let mut graph = load_asg(&asg_path);
    let sources = load_sources(&source_path);
    let reverse_sources = load_sources(&reverse_path);
    let mut coverage_cut = load_coverage_cut(&coverage_path);
    let (mut read_to_unitig, trio_labels, hom_global_coverage_set, hom_global_coverage) =
        load_read_to_unitig(&map_path);
    let reads = load_corrected_reads(&corrected_reads);
    assert_eq!(reads.total_reads as usize, trio_labels.len());

    let mut options = AssemblyOptions::default();
    options.global_homozygous_coverage = hom_global_coverage as i32;
    options.global_homozygous_coverage_is_set = hom_global_coverage_set as i32;

    let mut bubble_labels = BubbleLabel::new(1, graph.sequences.len() as u32);
    let mut trio_flags = trio_labels.clone();

    let output_directory = std::env::temp_dir().join("bp_graph_replay");
    fs::create_dir_all(&output_directory).expect("the scratch directory must be creatable");
    let output_prefix = output_directory.join("replay");
    let output_name = output_prefix.to_str().expect("the path must be UTF-8");

    output_bp_graph(
        &mut graph,
        &mut coverage_cut,
        output_name,
        &sources,
        &reverse_sources,
        &mut read_to_unitig,
        &mut bubble_labels,
        &mut options,
        &reads,
        &mut trio_flags,
    )
    .expect("output_bp_graph must produce the primary and haplotype contigs");

    for suffix in [
        ".p_ctg.gfa",
        ".p_ctg.noseq.gfa",
        ".p_ctg.lowQ.bed",
        ".hap1.p_ctg.gfa",
        ".hap1.p_ctg.noseq.gfa",
        ".hap1.p_ctg.lowQ.bed",
        ".hap2.p_ctg.gfa",
        ".hap2.p_ctg.noseq.gfa",
        ".hap2.p_ctg.lowQ.bed",
    ] {
        let produced = format!("{}{}", output_name, suffix);
        let golden = format!("{}{}", golden_prefix.display(), suffix);
        assert_same_file(Path::new(&produced), Path::new(&golden), suffix);
    }
}

fn bpprobe_dump_directory() -> PathBuf {
    std::env::var_os("RSHIFIASM_BPPROBE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("/tmp/bpprobe"))
}

/// The on-disk `kv_u_trans_t` layout written by the probe: the entry count, the
/// index length, the index itself and then 48-byte `u_trans_t` records.
fn load_transition_dump(path: &Path) -> (Vec<UTrans>, Vec<u64>) {
    let data = fs::read(path)
        .unwrap_or_else(|error| panic!("cannot read {}: {}", path.display(), error));
    let count = u32::from_le_bytes(data[0..4].try_into().unwrap()) as usize;
    let index_count = u32::from_le_bytes(data[4..8].try_into().unwrap()) as usize;
    let mut position = 8usize;
    let mut index = Vec::with_capacity(index_count);
    for _ in 0..index_count {
        index.push(u64::from_le_bytes(
            data[position..position + 8].try_into().unwrap(),
        ));
        position += 8;
    }

    let mut entries = Vec::with_capacity(count);
    for _ in 0..count {
        let record = &data[position..position + 48];
        position += 48;
        // `occ` at offset 24 is never initialised by the reference and the
        // remaining bytes are padding, so both are ignored.
        let flags = record[40];
        entries.push(UTrans {
            query_start: u32::from_le_bytes(record[0..4].try_into().unwrap()),
            query_end: u32::from_le_bytes(record[4..8].try_into().unwrap()),
            query_identifier: u32::from_le_bytes(record[8..12].try_into().unwrap()),
            target_start: u32::from_le_bytes(record[12..16].try_into().unwrap()),
            target_end: u32::from_le_bytes(record[16..20].try_into().unwrap()),
            target_identifier: u32::from_le_bytes(record[20..24].try_into().unwrap()),
            occurrence: 0,
            network_weight: f64::from_le_bytes(record[32..40].try_into().unwrap()),
            flag: flags & 0x3f,
            reverse: flags & 0x40 != 0,
            deleted: flags & 0x80 != 0,
        });
    }
    assert_eq!(position, data.len(), "the transition dump has trailing bytes");
    (entries, index)
}

fn assert_transitions_match(produced: &[UTrans], golden: &[UTrans], label: &str) {
    assert_eq!(
        produced.len(),
        golden.len(),
        "{}: entry count differs ({} vs {})",
        label,
        produced.len(),
        golden.len()
    );
    for (i, (a, b)) in produced.iter().zip(golden.iter()).enumerate() {
        assert_eq!(a.query_start, b.query_start, "{}: entry {} qs", label, i);
        assert_eq!(a.query_end, b.query_end, "{}: entry {} qe", label, i);
        assert_eq!(
            a.query_identifier, b.query_identifier,
            "{}: entry {} qn",
            label, i
        );
        assert_eq!(a.target_start, b.target_start, "{}: entry {} ts", label, i);
        assert_eq!(a.target_end, b.target_end, "{}: entry {} te", label, i);
        assert_eq!(
            a.target_identifier, b.target_identifier,
            "{}: entry {} tn",
            label, i
        );
        assert_eq!(
            a.network_weight.to_bits(),
            b.network_weight.to_bits(),
            "{}: entry {} nw",
            label,
            i
        );
        assert_eq!(a.flag, b.flag, "{}: entry {} f", label, i);
        assert_eq!(a.reverse, b.reverse, "{}: entry {} rev", label, i);
        assert_eq!(a.deleted, b.deleted, "{}: entry {} del", label, i);
    }
}

/// Replays the trio half of `output_bp_graph` and compares the transition
/// table before and after `clean_u_trans_t_idx`, plus the trio flags assigned
/// by `set_trio_flag_by_cov`, against the reference probe dumps.
#[test]
fn the_trio_stage_reproduces_the_reference_flags() {
    let dump = probe_dump_directory();
    let config = "d20_after_hic_clean";
    let asg_path = dump.join(format!("{}.all.debug.asg_t.bin", config));
    let source_path = dump.join(format!("{}.all.debug.source.bin", config));
    let reverse_path = dump.join(format!("{}.all.debug.reverse.bin", config));
    let coverage_path = dump.join(format!("{}.all.debug.coverage_cut.bin", config));
    let map_path = dump.join(format!("{}.all.debug.ruIndex.bin", config));
    let corrected_reads = repository_root().join("reference/write_paf_ec.ec.fa");

    let probe_dumps = bpprobe_dump_directory();
    let pre_path = probe_dumps.join("f15_k_trans_pre.bin");
    let post_path = probe_dumps.join("f16_k_trans_post.bin");
    let flags_path = probe_dumps.join("f17_trio_flag.bin");

    if !asg_path.exists()
        || !corrected_reads.exists()
        || !pre_path.exists()
        || !post_path.exists()
        || !flags_path.exists()
    {
        eprintln!("skipping: the reference dumps are missing");
        return;
    }

    let mut read_graph = load_asg(&asg_path);
    let sources = load_sources(&source_path);
    let reverse_sources = load_sources(&reverse_path);
    let mut coverage_cut = load_coverage_cut(&coverage_path);
    let (mut read_to_unitig, trio_labels, hom_global_coverage_set, hom_global_coverage) =
        load_read_to_unitig(&map_path);
    let _reads = load_corrected_reads(&corrected_reads);

    let mut options = AssemblyOptions::default();
    options.global_homozygous_coverage = hom_global_coverage as i32;
    options.global_homozygous_coverage_is_set = hom_global_coverage_set as i32;

    hic_clean(&mut read_graph);
    let unitig_graph = ma_ug_gen_primary(&mut read_graph, PRIMARY_LABEL);
    let mut copy_read_graph = copy_read_graph(&read_graph);
    let mut copy_unitig_graph = copy_untig_graph(&unitig_graph);
    let mut bubble_labels = BubbleLabel::new(1, read_graph.sequences.len() as u32);
    let mut trio_flags = trio_labels.clone();

    let mut cov = adjust_utg_by_primary(
        &mut copy_unitig_graph,
        &mut copy_read_graph,
        &mut coverage_cut,
        &sources,
        &reverse_sources,
        &mut read_to_unitig,
        i64::from(options.maximum_short_tip_read_count * 2),
        0.15,
        3,
        0.05,
        0.9,
        options.maximum_hang_length,
        options.minimum_overlap_length,
        &[],
        Some(&mut bubble_labels),
        &mut options,
        &mut trio_flags,
        true,
        true,
    )
    .expect("adjust_utg_by_primary must succeed");

    let (golden_pre, _) = load_transition_dump(&pre_path);
    {
        let t_ch = cov.t_ch.as_mut().expect("the trans chain must exist");
        assert_transitions_match(&t_ch.k_trans, &golden_pre, "pre-clean");

        clean_u_trans_t_idx(
            &mut t_ch.k_trans,
            &mut t_ch.k_trans_idx,
            &unitig_graph.unitigs,
            &unitig_graph.graph,
            &read_graph,
        );
    }

    let (golden_post, golden_index) = load_transition_dump(&post_path);
    {
        let t_ch = cov.t_ch.as_ref().expect("the trans chain must exist");
        assert_transitions_match(&t_ch.k_trans, &golden_post, "post-clean");
        assert_eq!(
            t_ch.k_trans_idx, golden_index,
            "post-clean: the transition index differs"
        );
    }

    let t_ch = cov.t_ch.as_ref().expect("the trans chain must exist");
    set_trio_flag_by_cov(
        &unitig_graph.unitigs,
        &read_graph,
        &cov.is_r_het,
        &t_ch.k_trans,
        &t_ch.k_trans_idx,
        &mut trio_flags,
    );

    let golden_flags = fs::read(&flags_path).expect("the trio flag dump must be readable");
    let read_count = u32::from_le_bytes(golden_flags[0..4].try_into().unwrap()) as usize;
    assert_eq!(trio_flags.len(), read_count);
    assert_eq!(
        trio_flags.as_slice(),
        &golden_flags[4..4 + read_count],
        "the trio flags differ"
    );
}

/// The unitig list part of the probe dump: `u32 n` followed by, per unitig,
/// `u32 len`, `u8 circ`, `u32 start`, `u32 end`, `u32 n_reads` and the reads.
fn load_unitig_dump(path: &Path) -> Vec<(u32, bool, u32, u32, Vec<u64>)> {
    let data = fs::read(path)
        .unwrap_or_else(|error| panic!("cannot read {}: {}", path.display(), error));
    let mut position = 0usize;
    let count = u32::from_le_bytes(data[0..4].try_into().unwrap()) as usize;
    position += 4;
    let mut unitigs = Vec::with_capacity(count);
    for _ in 0..count {
        let length = u32::from_le_bytes(data[position..position + 4].try_into().unwrap());
        position += 4;
        let circular = data[position] != 0;
        position += 1;
        let start = u32::from_le_bytes(data[position..position + 4].try_into().unwrap());
        position += 4;
        let end = u32::from_le_bytes(data[position..position + 4].try_into().unwrap());
        position += 4;
        let read_count =
            u32::from_le_bytes(data[position..position + 4].try_into().unwrap()) as usize;
        position += 4;
        let mut reads = Vec::with_capacity(read_count);
        for _ in 0..read_count {
            reads.push(u64::from_le_bytes(
                data[position..position + 8].try_into().unwrap(),
            ));
            position += 8;
        }
        unitigs.push((length, circular, start, end, reads));
    }
    assert_eq!(position, data.len(), "the unitig dump has trailing bytes");
    unitigs
}

fn assert_graphs_match(
    produced: &Asg,
    golden: &Asg,
    label: &str,
    compare_classes: bool,
    compare_flags: bool,
) {
    assert_eq!(
        produced.sequences.len(),
        golden.sequences.len(),
        "{}: sequence count differs",
        label
    );
    for (i, (a, b)) in produced.sequences.iter().zip(&golden.sequences).enumerate() {
        assert_eq!(a.length, b.length, "{}: sequence {} length", label, i);
        assert_eq!(a.deleted, b.deleted, "{}: sequence {} deleted", label, i);
        if compare_classes {
            assert_eq!(
                a.coverage_class, b.coverage_class,
                "{}: sequence {} class",
                label, i
            );
        }
    }
    assert_eq!(produced.arcs.len(), golden.arcs.len(), "{}: arc count", label);
    for (i, (a, b)) in produced.arcs.iter().zip(&golden.arcs).enumerate() {
        assert_eq!(a.ul, b.ul, "{}: arc {} ul", label, i);
        assert_eq!(a.v, b.v, "{}: arc {} v", label, i);
        assert_eq!(a.overlap, b.overlap, "{}: arc {} overlap", label, i);
        assert_eq!(a.deleted, b.deleted, "{}: arc {} deleted", label, i);
        assert_eq!(a.occurrence, b.occurrence, "{}: arc {} occurrence", label, i);
        if compare_flags {
            assert_eq!(a.strong, b.strong, "{}: arc {} strong", label, i);
            assert_eq!(
                a.without_large_indel, b.without_large_indel,
                "{}: arc {} without_large_indel",
                label, i
            );
        }
        assert_eq!(
            a.exact_region, b.exact_region,
            "{}: arc {} exact_region",
            label, i
        );
    }
}

/// Replays `output_bp_graph` into `ma_ug_gen_phase` and compares both the read
/// graph handed to `reduce_hamming_error_adv` and the resulting phase unitigs
/// with the reference probe dumps.
#[test]
fn the_phase_unitig_graph_reproduces_the_reference() {
    let dump = probe_dump_directory();
    let config = "d20_after_hic_clean";
    let asg_path = dump.join(format!("{}.all.debug.asg_t.bin", config));
    let source_path = dump.join(format!("{}.all.debug.source.bin", config));
    let reverse_path = dump.join(format!("{}.all.debug.reverse.bin", config));
    let coverage_path = dump.join(format!("{}.all.debug.coverage_cut.bin", config));
    let map_path = dump.join(format!("{}.all.debug.ruIndex.bin", config));
    let corrected_reads = repository_root().join("reference/write_paf_ec.ec.fa");

    let probe_dumps = PathBuf::from(
        std::env::var_os("RSHIFIASM_TRIOPROBE_DIR")
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from("/tmp/trioprobe")),
    );
    let read_graph_path = probe_dumps.join("f20_after_reduce_hamming.asg_t.bin");
    let phase_ug_path = probe_dumps.join("f21_after_ma_ug_gen_phase.ug.bin");
    let phase_graph_path = probe_dumps.join("f21_after_ma_ug_gen_phase.asg_t.bin");

    if !asg_path.exists() || !corrected_reads.exists() || !read_graph_path.exists() {
        eprintln!("skipping: the reference dumps are missing");
        return;
    }

    let mut read_graph = load_asg(&asg_path);
    let sources = load_sources(&source_path);
    let reverse_sources = load_sources(&reverse_path);
    let mut coverage_cut = load_coverage_cut(&coverage_path);
    let (mut read_to_unitig, trio_labels, hom_global_coverage_set, hom_global_coverage) =
        load_read_to_unitig(&map_path);
    let reads = load_corrected_reads(&corrected_reads);

    let mut options = AssemblyOptions::default();
    options.global_homozygous_coverage = hom_global_coverage as i32;
    options.global_homozygous_coverage_is_set = hom_global_coverage_set as i32;

    let mut bubble_labels = BubbleLabel::new(1, read_graph.sequences.len() as u32);
    let mut trio_flags = trio_labels.clone();
    let output_directory = std::env::temp_dir().join("bp_graph_phase_replay");
    fs::create_dir_all(&output_directory).expect("the scratch directory must be creatable");
    let output_name = output_directory.join("replay");
    let output_name = output_name.to_str().expect("the path must be UTF-8");

    output_bp_graph(
        &mut read_graph,
        &mut coverage_cut,
        output_name,
        &sources,
        &reverse_sources,
        &mut read_to_unitig,
        &mut bubble_labels,
        &mut options,
        &reads,
        &mut trio_flags,
    )
    .expect("output_bp_graph must produce the primary and haplotype contigs");

    // Replay the `f20`/`f21` states on a fresh clone: the driver's
    // `reduce_hamming_error_adv` inserted no edge and fixed no bubble on this
    // corpus, so the read graph is unchanged when `ma_ug_gen_phase` runs.
    let mut replay_graph = load_asg(&asg_path);
    let mut replay_flags = trio_labels.clone();
    reduce_hamming_error_adv(
        None,
        &mut replay_graph,
        &sources,
        &coverage_cut,
        options.maximum_hang_length,
        options.minimum_overlap_length,
        options.gap_fuzz,
        &mut replay_flags,
        options.maximum_hang_rate,
    );
    let golden_read_graph = load_asg(&read_graph_path);
    assert_graphs_match(
        &replay_graph,
        &golden_read_graph,
        "before ma_ug_gen_phase",
        true,
        true,
    );

    let homozygous = f64::from(options.global_homozygous_coverage);
    let base = if options.global_homozygous_coverage_is_set != 0 {
        homozygous
    } else {
        homozygous / HOM_PEAK_RATE
    };
    let phase_threshold = ((base * PHASE_SEF) as u64).max(PHASE_SEP);
    let phase = ma_ug_gen_phase(
        &mut replay_graph,
        phase_threshold as u32,
        PHASE_SEP_RATE,
        &replay_flags,
    );

    let golden_unitigs = load_unitig_dump(&phase_ug_path);
    assert_eq!(
        phase.unitigs.len(),
        golden_unitigs.len(),
        "phase unitig count differs"
    );
    for (i, (unitig, golden)) in phase.unitigs.iter().zip(&golden_unitigs).enumerate() {
        assert_eq!(unitig.length, golden.0, "phase unitig {} length", i);
        assert_eq!(unitig.circular, golden.1, "phase unitig {} circular", i);
        assert_eq!(unitig.start, golden.2, "phase unitig {} start", i);
        assert_eq!(unitig.end, golden.3, "phase unitig {} end", i);
        assert_eq!(unitig.reads, golden.4, "phase unitig {} reads", i);
    }

    let golden_phase_graph = load_asg(&phase_graph_path);
    // The reference leaves the sequence classes and the `strong`/`no_l_indel`
    // arc flags uninitialised in the freshly allocated unitig graph; the next
    // stage overwrites the classes for every live unitig.
    assert_graphs_match(
        &phase.graph,
        &golden_phase_graph,
        "ma_ug_gen_phase graph",
        false,
        false,
    );
}

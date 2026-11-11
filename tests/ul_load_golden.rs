//! Golden replay of the `--ul` first pass against the reference run: the
//! `ul_load` orchestration (scall, `.uidx.bin` cache, `bl` labelling) and the
//! `gen_init_sg` UL branch (contained-read consensus cleanup).
//!
//! Requires `/tmp/ulref.uidx.bin` (the reference unitig index), `/tmp/ul.fq`
//! (the synthetic UL reads), `reference/ul.ul.ovlp.bin` (the reference's
//! post-scall `UL_INF` snapshot) and `/tmp/ul_consensus_after.txt` (the C++
//! probe's post-consensus `UL_INF` dump, written by
//! `g++ -O1 -I /home/keith/projects/hifiasm /tmp/ul_consensus_probe.cpp`).
//! Ignored by default.

use rshifiasm::assembly::hit::{
    normalize_ma_hit_t_single_side_advance, renew_graph_init, try_rescue_overlaps, ReadToUnitig,
};
use rshifiasm::assembly::read_cache;
use rshifiasm::assembly::ul_load::{ul_load, UL_COV_THRES};
use rshifiasm::command_line::options::{AssemblyOptions, FileList, MINIMUM_CHAIN_COUNT};

const REFERENCE_PREFIX: &str = "/tmp/ulref";
const SCRATCH_PREFIX: &str = "/tmp/ulout";

/// Compares two `.uidx.bin` dumps, ignoring the bytes the reference never
/// reads back (unused `khashl` buckets). Returns the first offending offset.
fn compare_uidx(actual: &[u8], expected: &[u8]) -> Option<usize> {
    fn read_u32(source: &[u8], cursor: &mut usize) -> u32 {
        let value = u32::from_le_bytes(source[*cursor..*cursor + 4].try_into().unwrap());
        *cursor += 4;
        value
    }
    fn read_u64(source: &[u8], cursor: &mut usize) -> u64 {
        let value = u64::from_le_bytes(source[*cursor..*cursor + 8].try_into().unwrap());
        *cursor += 8;
        value
    }
    // One `khashl` map dump: `{u64 key, value}` buckets; `value_bytes` is 2
    // for the filter table (`int16`) and 8 for the position partitions.
    fn compare_map(
        actual: &[u8],
        expected: &[u8],
        cursor: &mut usize,
        value_bytes: usize,
    ) -> Option<usize> {
        let bucket_count = read_u32(actual, cursor);
        debug_assert_eq!(bucket_count, read_u32(expected, cursor));
        *cursor += 8; // bits + count
        let mut used = Vec::new();
        let flag_size = if bucket_count < 32 { 1 } else { bucket_count >> 5 } as usize;
        for _ in 0..flag_size {
            let word = read_u32(actual, cursor);
            debug_assert_eq!(word, read_u32(expected, cursor));
            used.push(word);
        }
        *cursor += 1; // keys flag
        let bucket_bytes = 8 + value_bytes;
        let is_used = |bucket: u32| (used[(bucket >> 5) as usize] >> (bucket & 31)) & 1 != 0;
        for bucket in 0..bucket_count {
            let start = *cursor;
            *cursor += bucket_bytes;
            if is_used(bucket) && actual[start..*cursor] != expected[start..*cursor] {
                return Some(start);
            }
        }
        None
    }

    let mut cursor = 1; // the 'f' marker
    compare_map(actual, expected, &mut cursor, 2)?;
    cursor += 1; // the 'h' marker
    cursor += 4 + 4 + 8 + 8; // k, pre, tot, tot_pos
    let partitions = u32::from_le_bytes(actual[cursor - 20..cursor - 16].try_into().unwrap());
    for _ in 0..partitions {
        compare_map(actual, expected, &mut cursor, 8)?;
        let length = read_u64(actual, &mut cursor) as usize;
        debug_assert_eq!(length, read_u64(expected, &mut cursor) as usize);
        for _ in 0..length {
            let start = cursor;
            cursor += 16;
            if actual[start..cursor] != expected[start..cursor] {
                return Some(start);
            }
        }
    }
    None
}

/// Asserts `store` matches one `R` line of the C++ probe's dump.
fn assert_matches_dump(store: &rshifiasm::assembly::ul_store::AllUl, dump_path: &str) {
    let dump = std::fs::read_to_string(dump_path)
        .unwrap_or_else(|error| panic!("{dump_path} must be present: {error}"));
    let mut lines = dump.lines();
    for (read_id, read) in store.reads.iter().enumerate() {
        let line = lines.next().expect("the dump must cover every read");
        let fields: Vec<&str> = line.split_whitespace().collect();
        assert_eq!(fields[0], "R");
        assert_eq!(fields[1], read_id.to_string());
        assert_eq!(fields[2], read.rlen.to_string(), "rlen at read {}", read_id);
        assert_eq!(fields[3], read.dd.to_string(), "dd at read {}", read_id);
        let block_count: usize = fields[4].parse().unwrap();
        assert_eq!(
            block_count,
            read.bb.len(),
            "block count at read {}",
            read_id
        );
        for (block_index, block) in fields[5..5 + block_count].iter().enumerate() {
            let values: Vec<u64> = block.split(':').map(|v| v.parse().unwrap()).collect();
            let actual = &read.bb[block_index];
            let expected = [
                actual.hid as u64,
                actual.qs as u64,
                actual.qe as u64,
                actual.ts as u64,
                actual.te as u64,
                actual.pidx as u64,
                actual.pdis as u64,
                actual.aidx as u64,
                actual.pchain as u64,
                u64::from(actual.rev),
                u64::from(actual.base),
                u64::from(actual.el),
            ];
            assert_eq!(
                values.as_slice(),
                &expected[..],
                "block {} of read {}",
                block_index,
                read_id
            );
        }
    }
}

#[test]
#[ignore]
fn ul_load_reproduces_the_reference_first_pass() {
    let mut options = AssemblyOptions::default();
    // The restart caches live at the reference prefix; the run itself writes
    // to the scratch prefix so the reference files are never touched.
    options.output_file_name = REFERENCE_PREFIX.to_string();
    options.ultra_long_read_files = Some(FileList {
        file_names: vec!["/tmp/ul.fq".to_string()],
        name_lengths: vec!["/tmp/ul.fq".len()],
    });
    for suffix in [".uidx.bin", ".ul.ovlp.bin"] {
        let _ = std::fs::remove_file(format!("{}{}", SCRATCH_PREFIX, suffix));
    }

    let (reads, state) = read_cache::load_all_data_from_disk(&mut options)
        .expect("the reference caches must load");
    options.output_file_name = SCRATCH_PREFIX.to_string();
    // `ha_opt_update_cov_min` at the head of `build_string_graph_without_clean`.
    options.update_coverage_with_minimum(options.homozygous_coverage, MINIMUM_CHAIN_COUNT);

    let mut read_to_unitig = ReadToUnitig::new(reads.total_reads as usize);
    let (mut forward, forward_flags): (Vec<Vec<_>>, Vec<_>) = state
        .forward
        .iter()
        .map(|list| (list.hits.clone(), (list.is_fully_corrected, list.is_abnormal)))
        .unzip();
    let (mut reverse, _reverse_flags): (Vec<Vec<_>>, Vec<_>) = state
        .reverse
        .iter()
        .map(|list| (list.hits.clone(), (list.is_fully_corrected, list.is_abnormal)))
        .unzip();
    try_rescue_overlaps(&mut forward, &mut reverse, &reads.read_lengths, 4, false);
    renew_graph_init(&mut forward, &mut reverse, &mut read_to_unitig);
    normalize_ma_hit_t_single_side_advance(&mut forward, false, None, 0);
    normalize_ma_hit_t_single_side_advance(&mut reverse, false, None, 0);

    let mut lists = forward;
    let _ = forward_flags;
    let mut store = ul_load(&mut options, &reads, &mut lists);

    // `print_ovlp_src_bl_stat`'s two lines, recomputed from the labelled
    // sources: 5626/1184/1028 and 739/736/2/1 in the reference run.
    let mut covered = 0u64;
    let mut chained = 0u64;
    for read_id in 0..reads.total_reads as usize {
        if store.ridx.idx[read_id + 1] - store.ridx.idx[read_id] == 0 {
            continue;
        }
        covered += 1;
        if lists[read_id].iter().any(|hit| hit.block_length != 0) {
            chained += 1;
        }
    }
    assert_eq!(covered, 1184, "covered HiFi reads");
    assert_eq!(chained, 1028, "chained HiFi reads");
    let mut status = [0u64; 4];
    for read in &store.reads {
        status[read.dd as usize] += 1;
    }
    assert_eq!(status[0] + status[1] + status[2] + status[3], 739);
    assert_eq!(status[1], 736, "fully corrected UL reads");
    assert_eq!(status[2], 2, "almost fully corrected UL reads");
    assert_eq!(status[3], 1, "UL reads with primary chains");

    // The `.uidx.bin` cache: the reference writes its khashl tables raw, so
    // the *unused* bucket slots carry stale heap bytes and the file cannot be
    // compared byte for byte (same class as the `.ul.ovlp.bin` exemption).
    // Everything the reference actually reads back — the used buckets of both
    // tables and the occurrence lists — must match exactly.
    let actual = std::fs::read(format!("{}.uidx.bin", SCRATCH_PREFIX))
        .expect("the run must write the unitig index cache");
    let expected = std::fs::read(format!("{}.uidx.bin", REFERENCE_PREFIX))
        .expect("the reference unitig index must be present");
    assert_eq!(actual.len(), expected.len(), "uidx.bin size");
    assert_eq!(compare_uidx(&actual, &expected), None, "uidx.bin content");

    // The `UL_INF` snapshot must match the reference field by field (the
    // `.ul.ovlp.bin` bytes themselves carry uninitialized padding).
    let reference = rshifiasm::assembly::ul_store::load_all_ul_t("reference/ul")
        .expect("the reference UL cache must load")
        .expect("the reference UL cache must exist");
    assert_eq!(reference.nid.len(), store.nid.len(), "UL read count");
    for (index, (expected, actual)) in reference.nid.iter().zip(store.nid.iter()).enumerate() {
        assert_eq!(expected, actual, "UL read name mismatch at {}", index);
    }
    for (index, (expected, actual)) in
        reference.reads.iter().zip(store.reads.iter()).enumerate()
    {
        assert_eq!(expected.rlen, actual.rlen, "rlen mismatch at {}", index);
        assert_eq!(expected.dd, actual.dd, "dd mismatch at {}", index);
        assert_eq!(
            expected.bb.len(),
            actual.bb.len(),
            "block count mismatch at {}",
            index
        );
        for (block_index, (expected, actual)) in
            expected.bb.iter().zip(actual.bb.iter()).enumerate()
        {
            assert_eq!(
                expected, actual,
                "block mismatch at read {} block {}",
                index, block_index
            );
        }
    }

    // `gen_init_sg`'s UL branch (`Overlaps.cpp` 39394): the whole first graph
    // pass, whose contained cleanup must match the reference log anchors
    // (try:152/done:122 and discard:47) and the C++ probe's post-consensus
    // dump.
    let initial = rshifiasm::assembly::graph_build::gen_init_sg(
        options.minimum_overlap_coverage,
        &mut lists,
        &reads.read_lengths,
        &mut read_to_unitig,
        options.maximum_hang_length,
        options.minimum_overlap_length,
        options.gap_fuzz,
        i64::from(options.minimum_overlap_length),
        options.thread_count,
        &reads,
        &mut options,
        Some(&mut store),
    );
    let outcome = initial
        .ultra_long_outcome
        .expect("the UL branch must run the contained cleanup");
    assert_eq!(outcome.tried, 152, "work_ul_gchains_consensus # try");
    assert_eq!(outcome.corrected, 122, "work_ul_gchains_consensus # done");
    assert_eq!(outcome.discarded, 47, "clean_contain_g0 # discard cread");
    assert_matches_dump(&store, "/tmp/ul_consensus_after.txt");
    let _ = UL_COV_THRES;

    for suffix in [".uidx.bin", ".ul.ovlp.bin"] {
        let _ = std::fs::remove_file(format!("{}{}", SCRATCH_PREFIX, suffix));
    }
}

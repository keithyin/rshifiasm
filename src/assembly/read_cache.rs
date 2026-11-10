//! The restart caches, the counterparts of `write_all_data_to_disk` /
//! `load_all_data_from_disk` (`Overlaps.cpp` 23567 / 23590) and the
//! `write_All_reads` / `load_All_reads` / `write_ma_hit_ts` / `write_ma`
//! formats (`Process_Read.cpp` 69 / 127, `Overlaps.cpp` 23442 / 23420).
//!
//! Three files are written below the output prefix: `<out>.ec.bin` (the
//! corrected read store plus the homozygous/heterozygous peaks),
//! `<out>.ovlp.source.bin` and `<out>.ovlp.reverse.bin` (the per-read
//! overlap lists of the final round, before `try_rescue_overlaps`). A
//! following run loads them and skips the correction stage entirely.
//!
//! Two regions of `ec.bin` are *uninitialised memory* in the reference and
//! therefore not byte-reproducible: the padding bits of each read's final
//! packed-sequence byte and the `name_index` slots beyond the sentinel
//! (`index_size` is the allocation capacity, not the entry count). The
//! golden comparison for `ec.bin` has to mask both; the overlap caches are
//! fully deterministic.

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};

use crate::command_line::options::AssemblyOptions;
use crate::correction::ec_driver::{CorrectionState, PafList};
use crate::counting::minimizer_count::count_minimizers;
use crate::overlap::paf::PafHit;
use crate::reads::read_store::ReadStore;
use crate::sequence::encoding::PackedSequence;

/// `READ_INIT_NUMBER` and the growth step of `All_reads::index_size`
/// (`Process_Read.cpp` 421): start at 1000 and double-plus-two while below
/// `total_reads + 2`. The value is part of the `ec.bin` format because the
/// reference serialises the allocation capacity.
fn stored_capacity(read_count: u64) -> u64 {
    let mut size = 1000u64;
    while size < read_count + 2 {
        size = size * 2 + 2;
    }
    size
}

fn le_u64(value: u64) -> [u8; 8] {
    value.to_le_bytes()
}

fn le_u32(value: u32) -> [u8; 4] {
    value.to_le_bytes()
}

fn le_i32(value: i32) -> [u8; 4] {
    value.to_le_bytes()
}

/// `write_All_reads` (`Process_Read.cpp` 69): the corrected read store.
fn write_reads(
    output: &mut BufWriter<File>,
    options: &AssemblyOptions,
    reads: &ReadStore,
) -> std::io::Result<()> {
    eprintln!("Writing reads to disk... ");
    output.write_all(&le_i32(options.adapter_length))?;

    let read_count = reads.total_reads;
    let capacity = stored_capacity(read_count);
    output.write_all(&le_u64(capacity))?;
    output.write_all(&le_u64(capacity))?;
    output.write_all(&le_u64(read_count))?;
    output.write_all(&le_u64(reads.total_bases))?;

    let total_name_length: u64 = reads
        .read_names
        .iter()
        .map(|name| name.len() as u64)
        .sum();
    output.write_all(&le_u64(total_name_length))?;

    // Per read: the ambiguous-base count followed by the positions.
    for sequence in &reads.sequences {
        output.write_all(&le_u64(sequence.ambiguous_positions.len() as u64))?;
        for &position in &sequence.ambiguous_positions {
            output.write_all(&le_u64(position))?;
        }
    }

    for &length in &reads.read_lengths {
        output.write_all(&le_u64(length))?;
    }

    // `read_sperate`: `len/4+1` bytes per read. The Rust packing only holds
    // `ceil(len/4)` bytes, the reference's trailing padding byte is zeroed.
    for sequence in &reads.sequences {
        output.write_all(&sequence.packed_bases)?;
        let stored = sequence.length / 4 + 1;
        for _ in sequence.packed_bases.len() as u64..stored {
            output.write_all(&[0u8])?;
        }
    }

    for name in &reads.read_names {
        output.write_all(name.as_bytes())?;
    }

    // `name_index`: start offsets with a `total_name_length` sentinel; the
    // slots beyond the sentinel belong to the allocation capacity.
    let mut offset = 0u64;
    for name in &reads.read_names {
        output.write_all(&le_u64(offset))?;
        offset += name.len() as u64;
    }
    output.write_all(&le_u64(total_name_length))?;
    for _ in read_count + 1..capacity {
        output.write_all(&le_u64(0))?;
    }

    output.write_all(&reads.trio_labels)?;

    output.write_all(&le_i32(options.homozygous_coverage))?;
    output.write_all(&le_i32(options.heterozygous_coverage))?;

    // `write_All_reads`: in `is_sc` mode the `uint64_t` marker 1 is followed
    // by each read's packed quality, `read_length/sc_bn` rounded up bytes.
    // Only the first `ceil(len/4)` bytes carry levels; the reference writes
    // its uninitialised tail, the port zeroes it.
    if options.single_cell_mode != 0 {
        output.write_all(&le_u64(1))?;
        for (length, packed) in reads.read_lengths.iter().zip(reads.qualities.iter()) {
            output.write_all(packed)?;
            let stored = (length / 2 + (length % 2)) as usize;
            for _ in packed.len()..stored {
                output.write_all(&[0u8])?;
            }
        }
    }

    output.flush()?;
    eprintln!("Reads has been written.");
    Ok(())
}

/// `write_ma` (`Overlaps.cpp` 23420): one hit as a 42-byte field stream.
///
/// The bitfields are serialised unpacked: `qns` keeps its
/// `<read>:32 | <start>:32` packing while `ml`/`rev`/`bl`/`del` are written
/// as plain little-endian words.
fn write_hit(output: &mut impl Write, hit: &PafHit) -> std::io::Result<()> {
    output.write_all(&le_u64(
        ((hit.query_identifier as u64) << 32) | hit.query_start as u64,
    ))?;
    output.write_all(&le_u32(hit.query_end))?;
    output.write_all(&le_u32(hit.target_identifier))?;
    output.write_all(&le_u32(hit.target_start))?;
    output.write_all(&le_u32(hit.target_end))?;
    output.write_all(&[hit.exact_region_byte])?;
    output.write_all(&[u8::from(hit.without_large_indel)])?;
    output.write_all(&le_u32(hit.matching_length))?;
    output.write_all(&le_u32(u32::from(hit.reverse)))?;
    output.write_all(&le_u32(hit.block_length))?;
    output.write_all(&le_u32(u32::from(hit.deleted)))?;
    Ok(())
}

/// `write_ma_hit_ts` (`Overlaps.cpp` 23442): the per-read overlap lists.
fn write_hit_lists(
    output: &mut BufWriter<File>,
    lists: &[PafList],
) -> std::io::Result<()> {
    eprintln!("Writing ma_hit_ts to disk... ");
    output.write_all(&(lists.len() as i64).to_le_bytes())?;
    for list in lists {
        output.write_all(&[u8::from(list.is_fully_corrected)])?;
        output.write_all(&[list.is_abnormal])?;
        output.write_all(&le_u32(list.hits.len() as u32))?;
        for hit in &list.hits {
            write_hit(output, hit)?;
        }
    }
    output.flush()?;
    eprintln!("ma_hit_ts has been written.");
    Ok(())
}

/// `write_all_data_to_disk` (`Overlaps.cpp` 23567): the three cache files.
///
/// The trio yak-binning cache is only written in trio-binning mode.
pub fn write_all_data_to_disk(
    options: &AssemblyOptions,
    reads: &ReadStore,
    forward: &[PafList],
    reverse: &[PafList],
) -> std::io::Result<()> {
    let ec_name = format!("{}.ec.bin", options.output_file_name);
    let mut output = BufWriter::new(File::create(&ec_name)?);
    write_reads(&mut output, options, reads)?;
    drop(output);

    if options.trio_binning_requested()
        && options.trio_yak_file_names[0].is_some()
        && options.trio_yak_file_names[1].is_some()
    {
        crate::assembly::trio_bin::write_yak_binning(options, reads);
    }

    let source_name = format!("{}.ovlp.source.bin", options.output_file_name);
    let mut output = BufWriter::new(File::create(&source_name)?);
    write_hit_lists(&mut output, forward)?;

    let reverse_name = format!("{}.ovlp.reverse.bin", options.output_file_name);
    let mut output = BufWriter::new(File::create(&reverse_name)?);
    write_hit_lists(&mut output, reverse)?;
    drop(output);

    eprintln!("bin files have been written.");
    Ok(())
}

fn read_u64(input: &mut impl Read) -> std::io::Result<u64> {
    let mut buffer = [0u8; 8];
    input.read_exact(&mut buffer)?;
    Ok(u64::from_le_bytes(buffer))
}

fn read_u32(input: &mut impl Read) -> std::io::Result<u32> {
    let mut buffer = [0u8; 4];
    input.read_exact(&mut buffer)?;
    Ok(u32::from_le_bytes(buffer))
}

fn read_i32(input: &mut impl Read) -> std::io::Result<i32> {
    let mut buffer = [0u8; 4];
    input.read_exact(&mut buffer)?;
    Ok(i32::from_le_bytes(buffer))
}

/// Parses one read block (`<prefix>.ec[0-9]*.bin`) and appends it to
/// `store`; returns the stored homozygous/heterozygous peaks.
fn read_read_store(
    input: &mut impl Read,
    options: &AssemblyOptions,
    store: &mut ReadStore,
) -> std::io::Result<(i32, i32)> {
    let stored_adapter = read_i32(input)?;
    if stored_adapter != options.adapter_length {
        eprintln!(
            "the adapterLen of index is: {}, but the adapterLen set by user is: {}",
            stored_adapter, options.adapter_length
        );
        std::process::exit(1);
    }

    let _index_size = read_u64(input)?;
    let name_index_size = read_u64(input)?;
    let total_reads = read_u64(input)?;
    let total_bases = read_u64(input)?;
    let total_name_length = read_u64(input)?;

    let mut sequences = Vec::with_capacity(total_reads as usize);
    for _ in 0..total_reads {
        let count = read_u64(input)?;
        let mut ambiguous_positions = Vec::with_capacity(count as usize);
        for _ in 0..count {
            ambiguous_positions.push(read_u64(input)?);
        }
        sequences.push(PackedSequence {
            packed_bases: Vec::new(),
            length: 0,
            ambiguous_positions,
        });
    }

    let mut read_lengths = Vec::with_capacity(total_reads as usize);
    for _ in 0..total_reads {
        read_lengths.push(read_u64(input)?);
    }

    for (sequence, &length) in sequences.iter_mut().zip(&read_lengths) {
        // The reference stores `len/4+1` bytes; the Rust packing only keeps
        // the `ceil(len/4)` bytes that carry bases.
        let stored = length / 4 + 1;
        let mut packed = vec![0u8; stored as usize];
        input.read_exact(&mut packed)?;
        packed.truncate(((length + 3) / 4) as usize);
        sequence.packed_bases = packed;
        sequence.length = length;
    }

    let mut names = vec![0u8; total_name_length as usize];
    input.read_exact(&mut names)?;

    let mut name_index = Vec::with_capacity(name_index_size as usize);
    for _ in 0..name_index_size {
        name_index.push(read_u64(input)?);
    }
    let mut read_names = Vec::with_capacity(total_reads as usize);
    for read in 0..total_reads as usize {
        let start = name_index[read] as usize;
        let end = name_index[read + 1] as usize;
        read_names.push(
            String::from_utf8_lossy(&names[start..end]).into_owned(),
        );
    }

    let mut trio_labels = vec![0u8; total_reads as usize];
    input.read_exact(&mut trio_labels)?;

    let homozygous_coverage = read_i32(input)?;
    let heterozygous_coverage = read_i32(input)?;

    store.sequences.extend(sequences);
    store.read_lengths.extend(read_lengths);
    store.read_names.extend(read_names);
    store.trio_labels.extend(trio_labels);
    store.qualities.extend(vec![Vec::new(); total_reads as usize]);
    store.total_reads += total_reads;
    store.total_bases += total_bases;

    eprintln!("Reads has been loaded.");
    Ok((homozygous_coverage, heterozygous_coverage))
}

/// `load_All_reads` (`Process_Read.cpp` 127): rebuilds the read store and
/// overwrites the homozygous/heterozygous peaks.
fn load_reads(
    input: &mut BufReader<File>,
    options: &mut AssemblyOptions,
) -> std::io::Result<ReadStore> {
    let mut store = ReadStore::default();
    let (homozygous_coverage, heterozygous_coverage) = read_read_store(input, options, &mut store)?;
    options.homozygous_coverage = homozygous_coverage;
    options.heterozygous_coverage = heterozygous_coverage;

    // The `is_sc` tail is optional; marker 1 is followed by every read's
    // packed quality (`read_length/sc_bn` rounded up bytes).
    let mut tail = [0u8; 8];
    match input.read_exact(&mut tail) {
        Ok(()) if u64::from_le_bytes(tail) == 1 => {
            store.qualities.clear();
            store.qualities.reserve(store.total_reads as usize);
            for read in 0..store.total_reads as usize {
                let length = store.read_lengths[read];
                let stored = (length / 2 + (length % 2)) as usize;
                let mut packed = vec![0u8; stored];
                input.read_exact(&mut packed)?;
                packed.truncate(length.div_ceil(4) as usize);
                store.qualities.push(packed);
            }
        }
        _ => {}
    }

    Ok(store)
}

/// `append_All_reads` (`Process_Read.cpp` 260): appends `<prefix>.ec{id}.bin`
/// to `store` and consumes `<prefix>.ovlp{id}.source.bin`.
///
/// The reference reads every stored hit back only to discard it; the overlap
/// lists of the pair path are recomputed from the merged reads.
pub fn append_all_reads(
    store: &mut ReadStore,
    options: &mut AssemblyOptions,
    prefix: &str,
    id: u32,
) -> Option<Vec<(bool, u8)>> {
    let ec_name = format!("{}.ec{}.bin", prefix, id);
    let source_name = format!("{}.ovlp{}.source.bin", prefix, id);
    let ec_file = File::open(&ec_name).ok()?;
    let source_file = File::open(&source_name).ok()?;

    let mut input = BufReader::new(ec_file);
    match read_read_store(&mut input, options, store) {
        Ok((homozygous_coverage, heterozygous_coverage)) => {
            // `asm_opt.hom_cov += hom_cov0` (the reference accumulates).
            options.homozygous_coverage += homozygous_coverage;
            options.heterozygous_coverage += heterozygous_coverage;
        }
        Err(_) => return None,
    }

    read_source_overlap_flags(BufReader::new(source_file)).ok()
}

/// Reads the per-read `is_fully_corrected`/`is_abnormal` flags of one
/// `ma_hit_t_alloc` file and discards the hits: `append_All_reads` clears
/// every loaded list but keeps the two flags.
fn read_source_overlap_flags(mut input: impl Read) -> std::io::Result<Vec<(bool, u8)>> {
    let read_count = read_u64(&mut input)? as usize;
    let mut flags = [0u8; 2];
    let mut result = Vec::with_capacity(read_count);
    for _ in 0..read_count {
        input.read_exact(&mut flags)?;
        let hits = read_u32(&mut input)? as u64;
        result.push((flags[0] != 0, flags[1]));
        std::io::copy(
            &mut input.by_ref().take(hits * 42),
            &mut std::io::sink(),
        )?;
    }
    Ok(result)
}

fn read_hit(input: &mut impl Read) -> std::io::Result<PafHit> {
    let qns = read_u64(input)?;
    let query_end = read_u32(input)?;
    let target_identifier = read_u32(input)?;
    let target_start = read_u32(input)?;
    let target_end = read_u32(input)?;
    let mut flags = [0u8; 2];
    input.read_exact(&mut flags)?;
    let matching_length = read_u32(input)?;
    let reverse = read_u32(input)?;
    let block_length = read_u32(input)?;
    let deleted = read_u32(input)?;

    Ok(PafHit {
        query_identifier: (qns >> 32) as u32,
        query_start: qns as u32,
        query_end,
        target_identifier,
        target_start,
        target_end,
        matching_length,
        block_length,
        reverse: reverse != 0,
        without_large_indel: flags[1] != 0,
        exact_region: flags[0] != 0,
        exact_region_byte: flags[0],
        deleted: deleted != 0,
    })
}

/// `load_ma_hit_ts` (`Overlaps.cpp` 23328).
fn load_hit_lists(input: &mut BufReader<File>) -> std::io::Result<Vec<PafList>> {
    eprintln!("Loading ma_hit_ts from disk... ");
    let read_count = i64::from_le_bytes({
        let mut buffer = [0u8; 8];
        input.read_exact(&mut buffer)?;
        buffer
    });
    let mut lists = Vec::with_capacity(read_count as usize);
    for _ in 0..read_count {
        let mut flags = [0u8; 2];
        input.read_exact(&mut flags)?;
        let length = read_u32(input)?;
        let mut hits = Vec::with_capacity(length as usize);
        for _ in 0..length {
            hits.push(read_hit(input)?);
        }
        lists.push(PafList {
            hits,
            is_fully_corrected: flags[0] != 0,
            is_abnormal: flags[1],
        });
    }
    eprintln!("ma_hit_ts has been read.");
    Ok(lists)
}

/// `load_all_data_from_disk` (`Overlaps.cpp` 23590): `None` when the caches
/// are absent, so the caller falls back to a full run.
pub fn load_all_data_from_disk(
    options: &mut AssemblyOptions,
) -> Option<(ReadStore, CorrectionState)> {
    let ec_name = format!("{}.ec.bin", options.output_file_name);
    let ec_file = File::open(&ec_name).ok()?;
    let mut input = BufReader::new(ec_file);
    let mut reads = match load_reads(&mut input, options) {
        Ok(reads) => reads,
        Err(error) => {
            eprintln!("[WARN] failed to load {}: {}", ec_name, error);
            return None;
        }
    };
    drop(input);

    // `if (ha_opt_triobin(&asm_opt)) { if (!(yak0 && yak1 &&
    // load_yak_binning(...))) { ha_triobin(&asm_opt); write_yak_binning(...); } }`
    // — list-only runs replay the lists (and write their phase cache) too.
    if options.trio_binning_requested() {
        let dumps_ready = options.trio_yak_file_names[0].is_some()
            && options.trio_yak_file_names[1].is_some();
        if !(dumps_ready && crate::assembly::trio_bin::load_yak_binning(options, &mut reads)) {
            crate::assembly::trio_bin::run_trio_binning(options, &mut reads);
            crate::assembly::trio_bin::write_yak_binning(options, &reads);
        }
    }

    let source_name = format!("{}.ovlp.source.bin", options.output_file_name);
    let forward = match File::open(&source_name) {
        Ok(file) => {
            let mut input = BufReader::new(file);
            match load_hit_lists(&mut input) {
                Ok(lists) => lists,
                Err(error) => {
                    eprintln!("[WARN] failed to load {}: {}", source_name, error);
                    return None;
                }
            }
        }
        Err(_) => return None,
    };

    let reverse_name = format!("{}.ovlp.reverse.bin", options.output_file_name);
    let reverse = match File::open(&reverse_name) {
        Ok(file) => {
            let mut input = BufReader::new(file);
            match load_hit_lists(&mut input) {
                Ok(lists) => lists,
                Err(error) => {
                    eprintln!("[WARN] failed to load {}: {}", reverse_name, error);
                    return None;
                }
            }
        }
        Err(_) => return None,
    };

    Some((reads, CorrectionState { forward, reverse }))
}

/// `hap_recalculate_peaks` (`Assembly.cpp` 1881): with the corrected reads
/// already loaded, rebuild the k-mer filter and the minimizer statistics to
/// recover the homozygous/heterozygous peaks.
///
/// The reference frees and reloads the caches around the computation; the
/// net effect on the options is the same without it.
pub fn hap_recalculate_peaks(options: &mut AssemblyOptions, reads: &ReadStore) {
    let filter_table = if options.flags & crate::command_line::options::FLAG_NO_KMER_FILTER == 0 {
        let (filter_table, homozygous_coverage) =
            crate::counting::filter_table::generate_high_occurrence_filter(options, reads);
        options.update_coverage(homozygous_coverage);
        Some(filter_table)
    } else {
        None
    };

    let counts = count_minimizers(options, reads, filter_table.as_ref());
    options.homozygous_coverage = counts.homozygous_coverage;
    options.heterozygous_coverage = counts.heterozygous_coverage;
    eprintln!("M::hap_recalculate_peaks has done.");
}

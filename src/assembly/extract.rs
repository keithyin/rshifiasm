//! `--ex-list` / `--ex-iter`: dump the overlaps that connect a named set of
//! reads, the counterpart of `extract.cpp`.
//!
//! The reference runs this only when the restart caches were loaded: it takes
//! the names from a comma list or an `@file`, seeds those reads at round `-1`,
//! and then walks `n_rounds` of breadth-first layers over both overlap lists,
//! printing every edge that leaves a visited layer. The program exits right
//! after printing.

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

use crate::command_line::options::AssemblyOptions;
use crate::correction::ec_driver::{CorrectionState, PafList};
use crate::reads::read_store::ReadStore;

/// `gv_read_list`: `name1,name2,...`, or `@file` with one name per line
/// (gzip supported, the name being the first whitespace-delimited token).
fn read_list(spec: &str) -> Vec<String> {
    if let Some(path) = spec.strip_prefix('@') {
        return read_list_file(path);
    }
    spec.split(',').map(str::to_string).collect()
}

/// `gv_read_list`'s `@file` half.
fn read_list_file(path: &str) -> Vec<String> {
    let raw = match File::open(Path::new(path)) {
        Ok(file) => file,
        Err(_) => return Vec::new(),
    };
    let mut probe = [0u8; 2];
    let is_gzip = {
        let mut file = match File::open(Path::new(path)) {
            Ok(file) => file,
            Err(_) => return Vec::new(),
        };
        file.read(&mut probe).unwrap_or(0) == 2 && probe == [0x1f, 0x8b]
    };
    let reader: Box<dyn Read> = if is_gzip {
        Box::new(flate2::read::MultiGzDecoder::new(raw))
    } else {
        Box::new(raw)
    };

    let mut names = Vec::new();
    for line in BufReader::new(reader).lines() {
        let line = match line {
            Ok(line) => line,
            Err(_) => break,
        };
        let end = line
            .bytes()
            .position(|byte| matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r'))
            .unwrap_or(line.len());
        names.push(line[..end].to_string());
    }
    names
}

/// `ha_extract_print`: prints the overlap layers around the named reads.
pub fn print_overlap_layers<W: Write>(
    reads: &ReadStore,
    forward: &[PafList],
    reverse: &[PafList],
    rounds: i32,
    list_spec: &str,
    output: &mut W,
) -> io::Result<()> {
    let names = read_list(list_spec);
    if names.is_empty() {
        return Ok(());
    }

    let wanted: HashSet<&str> = names.iter().map(String::as_str).collect();
    // The reference's `h64_t`: the discovery round of every selected read.
    let mut layer: Vec<i32> = vec![i32::MAX; reads.total_reads as usize];
    for (read, name) in reads.read_names.iter().enumerate() {
        if wanted.contains(name.as_str()) {
            layer[read] = -1;
        }
    }

    // The layered passes only grow the discovered set; nothing is printed
    // until the final dump below.
    for round in 0..rounds {
        for read in 0..reads.total_reads as usize {
            for (reverse_list, lists) in [forward, reverse].into_iter().enumerate() {
                let hits = match lists.get(read) {
                    Some(list) => &list.hits,
                    None => continue,
                };
                for hit in hits {
                    let query = hit.query_identifier as usize;
                    let target = hit.target_identifier as usize;
                    let query_hit = layer[query] != i32::MAX && layer[query] < round;
                    let target_hit = layer[target] != i32::MAX && layer[target] < round;
                    if (!query_hit && !target_hit) || (query_hit && target_hit) {
                        continue;
                    }
                    if !query_hit && layer[query] == i32::MAX {
                        layer[query] = round;
                    }
                    if !target_hit && layer[target] == i32::MAX {
                        layer[target] = round;
                    }
                }
            }
        }
    }

    // The final dump: every edge that touches a discovered read, regardless of
    // its round. This runs after the layered passes and does not update the
    // layers.
    for read in 0..reads.total_reads as usize {
        for (reverse_list, lists) in [forward, reverse].into_iter().enumerate() {
            let hits = match lists.get(read) {
                Some(list) => &list.hits,
                None => continue,
            };
            for hit in hits {
                let query = hit.query_identifier as usize;
                let target = hit.target_identifier as usize;
                if layer[query] == i32::MAX && layer[target] == i32::MAX {
                    continue;
                }
                write!(
                    output,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    reads.read_name(query),
                    reads.read_length(query),
                    hit.query_start,
                    hit.query_end,
                    if hit.reverse { '-' } else { '+' },
                    reads.read_name(target),
                    reads.read_length(target),
                    hit.target_start,
                    hit.target_end,
                    hit.matching_length,
                    hit.block_length,
                    u8::from(reverse_list == 0),
                )?;
            }
        }
    }
    output.flush()
}

/// The `--ex-list` dispatch of `ha_assemble` (`Assembly.cpp` 2065): print the
/// layers to stdout. Only reachable when the restart caches were loaded.
pub fn print_requested_list(
    options: &AssemblyOptions,
    reads: &ReadStore,
    state: &CorrectionState,
) -> io::Result<()> {
    let spec = match options.extract_list_file_name.as_deref() {
        Some(list) => list,
        None => return Ok(()),
    };
    let stdout = io::stdout();
    let mut output = io::BufWriter::new(stdout.lock());
    print_overlap_layers(
        reads,
        &state.forward,
        &state.reverse,
        options.extract_iteration,
        spec,
        &mut output,
    )
}

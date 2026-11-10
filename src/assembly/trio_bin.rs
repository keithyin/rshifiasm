//! `ha_triobin` (`Trio.cpp`): the trio binning front end.
//!
//! Two input modes fill [`ReadStore::trio_labels`]: read-name lists
//! (`-3`/`-4`) and yak k-mer dumps (`-1`/`-2`). Both may run in one
//! invocation — the reference calls them in that order, so a read flagged by
//! both ends up with the yak verdict. The yak mode also has an on-disk
//! checkpoint (`write_yak_binning`/`load_yak_binning`), the `.hapN.phase.bin`
//! name lists keyed to the dump path.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};

use crate::assembly::trio_clean::DROP;
use crate::assembly::unitig::{AMBIGU, FATHER, MOTHER};
use crate::command_line::options::{AssemblyOptions, TRIO_FATHER, TRIO_MOTHER, TRIO_UNKNOWN};
use crate::container::khash::{hash_without_counter_bits, U64HashTable};
use crate::counting::hash::hash64_64;
use crate::logging;
use crate::reads::read_store::ReadStore;
use crate::sequence::encoding::SEQUENCE_NT4_TABLE;

/// `ha_triobin` (`Trio.cpp` 450): resets every label and runs the requested
/// binning modes.
pub fn run_trio_binning(options: &AssemblyOptions, reads: &mut ReadStore) {
    reads.trio_labels.fill(TRIO_UNKNOWN);
    // `if (opt->fn_bin_list[0] && opt->fn_bin_list[1]) ha_triobin_list(opt);`
    if options.trio_list_file_names[0].is_some() && options.trio_list_file_names[1].is_some() {
        ha_triobin_list(options, reads);
    }
    // `if (opt->fn_bin_yak[0] && opt->fn_bin_yak[1]) ha_triobin_yak(opt);`
    if options.trio_yak_file_names[0].is_some() && options.trio_yak_file_names[1].is_some() {
        ha_triobin_yak(options, reads);
    }
}

/// `ha_triobin_list` (`Trio.cpp` 328): indexes the read names and flags the
/// reads named in the paternal and maternal lists.
///
/// The reference keeps the names in a `cstr_ht_t` (a khashl string map);
/// nothing observable depends on its layout, so a plain hash map with the
/// same first-occurrence-wins rule stands in for it.
fn ha_triobin_list(options: &AssemblyOptions, reads: &mut ReadStore) {
    // `kh_val(h, k) = i` only when the key was absent: the first read with a
    // duplicated name owns it.
    let mut names: HashMap<String, u32> = HashMap::with_capacity(reads.total_reads as usize);
    for read_id in 0..reads.total_reads as u32 {
        let name = reads.read_name(read_id as usize);
        // `std` maps keep the first value for an existing key, matching the
        // reference's `if (absent) kh_val(h, k) = i;`.
        names.entry(name.to_string()).or_insert(read_id);
    }
    eprintln!(
        "[M::ha_triobin_list::{:.3}*{:.2}] created the hash table for read names",
        logging::realtime_seconds(),
        logging::cpu_usage()
    );
    ha_triobin_set_list(
        &names,
        options.trio_list_file_names[0].as_deref().expect("checked"),
        TRIO_FATHER,
        reads,
    );
    ha_triobin_set_list(
        &names,
        options.trio_list_file_names[1].as_deref().expect("checked"),
        TRIO_MOTHER,
        reads,
    );
    eprintln!(
        "[M::ha_triobin_list::{:.3}*{:.2}] ==> partitioned reads with external lists",
        logging::realtime_seconds(),
        logging::cpu_usage()
    );
}

/// `ha_triobin_set_list` (`Trio.cpp` 294): flags every read named in `path`.
///
/// Names are truncated at the first tab or space, and a name missing from
/// the store is silently skipped. `gzopen` transparently reads plain and
/// gzip files, so the same sniffing happens here.
fn ha_triobin_set_list(names: &HashMap<String, u32>, path: &str, flag: u8, reads: &mut ReadStore) {
    let Some(raw) = read_text(path) else {
        eprintln!("ERROR: failed to open file '{}'", path);
        // The reference returns -1 and the caller ignores it.
        return;
    };
    let text = raw;

    let mut total_lines = 0u64;
    let mut flagged_reads = 0u64;
    // `ks_getuntil(KS_SEP_LINE)` hands out every complete line including
    // empty ones and stops at the end of the file, so `split_inclusive`
    // over `\n` reproduces the count: a trailing newline adds no extra
    // entry, an unterminated last line still counts, and empty interior
    // lines count as lines that flag nothing.
    for chunk in text.split_inclusive(|byte| *byte == b'\n') {
        let line = chunk.strip_suffix(b"\n").unwrap_or(chunk);
        total_lines += 1;
        // `for (p = str.s; *p; ++p) if (*p == '\t' || *p == ' ') *p = 0;`
        let name = match line.iter().position(|byte| *byte == b'\t' || *byte == b' ') {
            Some(end) => &line[..end],
            None => {
                // A final `\r` stays part of the name, matching `kstream`,
                // which does not strip it either.
                line
            }
        };
        if let Ok(name) = std::str::from_utf8(name) {
            if let Some(read_id) = names.get(name) {
                reads.trio_labels[*read_id as usize] = flag;
                flagged_reads += 1;
            }
        }
    }
    eprintln!(
        "[M::ha_triobin_set_list::{:.3}*{:.2}] flagged {} reads, out of {} lines in file '{}'",
        logging::realtime_seconds(),
        logging::cpu_usage(),
        flagged_reads,
        total_lines,
        path
    );
}

/// Reads a file that may be plain or gzip compressed (`gzopen` + `kstream`).
fn read_text(path: &str) -> Option<Vec<u8>> {
    let raw = std::fs::File::open(path)
        .and_then(|mut file| {
            let mut raw = Vec::new();
            file.read_to_end(&mut raw)?;
            Ok(raw)
        })
        .ok()?;
    if raw.starts_with(&[0x1f, 0x8b]) {
        let mut plain = Vec::with_capacity(raw.len() * 4);
        flate2::read::MultiGzDecoder::new(raw.as_slice())
            .read_to_end(&mut plain)
            .ok()
            .map(|_| plain)
            .or(Some(raw))
    } else {
        Some(raw)
    }
}

/// `atoi` on a byte slice (leading sign and digits, trailing bytes ignored).
fn atoi(bytes: &[u8]) -> u32 {
    let mut index = 0usize;
    while index < bytes.len() && (bytes[index] == b' ' || bytes[index] == b'\t') {
        index += 1;
    }
    let mut sign = 1i64;
    if index < bytes.len() && (bytes[index] == b'+' || bytes[index] == b'-') {
        if bytes[index] == b'-' {
            sign = -1;
        }
        index += 1;
    }
    let mut value = 0i64;
    while index < bytes.len() && bytes[index].is_ascii_digit() {
        value = value * 10 + i64::from(bytes[index] - b'0');
        index += 1;
    }
    (value * sign) as i32 as u32
}

/// `phrase_hstatus` (`Trio.cpp` 373): split one line of the `-5` list into
/// the read name and the haplotype id encoded after the last `_` of the
/// first field.
fn phrase_hstatus(line: &[u8]) -> Option<(&[u8], u32)> {
    let first_sep = line
        .iter()
        .position(|byte| *byte == b'\t' || *byte == b' ')?;
    let rest = &line[first_sep + 1..];
    let second_sep = rest
        .iter()
        .position(|byte| *byte == b'\t' || *byte == b' ')
        .unwrap_or(rest.len());
    let read_name = &rest[..second_sep];
    if read_name.is_empty() {
        return None;
    }
    let first = &line[..first_sep];
    let underscore = first.iter().rposition(|byte| *byte == b'_')?;
    let hid = atoi(&first[underscore + 1..]);
    if hid == u32::MAX {
        return None;
    }
    Some((read_name, hid))
}

/// `ha_polybin_list` (`Trio.cpp` 391): the per-read haplotype bitmask of the
/// `-5` external list.
pub fn ha_polybin_list(path: &str, reads: &ReadStore) -> Result<Vec<u32>, String> {
    let mut names: HashMap<String, u32> = HashMap::with_capacity(reads.total_reads as usize);
    for read_id in 0..reads.total_reads as u32 {
        names
            .entry(reads.read_name(read_id as usize).to_string())
            .or_insert(read_id);
    }
    eprintln!(
        "[M::ha_polybin_list::{:.3}*{:.2}] created the hash table for read names",
        logging::realtime_seconds(),
        logging::cpu_usage()
    );

    let Some(text) = read_text(path) else {
        eprintln!("ERROR: failed to open file '{}'", path);
        return Err(format!("failed to open file '{}'", path));
    };

    let mut hap_s = vec![0u32; reads.total_reads as usize];
    for chunk in text.split_inclusive(|byte| *byte == b'\n') {
        let line = chunk.strip_suffix(b"\n").unwrap_or(chunk);
        let Some((read_name, hid)) = phrase_hstatus(line) else {
            eprintln!("ERROR: wrong hap status");
            continue;
        };
        let Ok(read_name) = std::str::from_utf8(read_name) else {
            continue;
        };
        if let Some(read_id) = names.get(read_name) {
            // `(((uint32_t)1)<<(hid-1))`; `wrapping_shl` reproduces the
            // hardware's shift-count masking for out-of-range ids.
            hap_s[*read_id as usize] |= 1u32.wrapping_shl(hid.wrapping_sub(1) & 31);
        }
    }

    eprintln!(
        "[M::ha_polybin_list::{:.3}*{:.2}] ==> partitioned reads with external lists",
        logging::realtime_seconds(),
        logging::cpu_usage()
    );
    Ok(hap_s)
}

/// `update_poly_trio` (`Overlaps.cpp` 14646): apply one haplotype bit of the
/// external list to the read labels. Reads outside every bin are labelled
/// [`MOTHER`]; reads already [`DROP`]ped keep that label.
pub fn update_poly_trio(mask: u32, hap_s: &[u32], trio_flags: &mut [u8]) {
    for (label, &haps) in trio_flags.iter_mut().zip(hap_s.iter()) {
        if *label == DROP as u8 {
            continue;
        }
        *label = AMBIGU;
        if haps == 0 {
            continue;
        }
        *label = if haps & mask != 0 { FATHER } else { MOTHER };
    }
}

// ---------------------------------------------------------------------------
// `ha_triobin_yak` (`Trio.cpp` 261): the yak k-mer dump front end.

/// `YAK_COUNTER_BITS` of `Trio.cpp` (yak writes 10, unlike hifiasm's 12).
const YAK_COUNTER_BITS: u32 = 10;
const YAK_MAX_COUNT: u64 = (1 << YAK_COUNTER_BITS) - 1;
const YAK_MAGIC: [u8; 4] = [b'Y', b'A', b'K', 2];

/// `yak_ch_hash`: the table hash ignores the packed count (the low bits).
fn yak_ch_hash(key: u64) -> u32 {
    hash_without_counter_bits(key, YAK_COUNTER_BITS)
}

/// `yak_hash64`: Thomas Wang's hash truncated into `mask` after every step,
/// exactly like the `Trio.cpp` copy (the masked variant differs from applying
/// `hash64_64` and masking once at the end).
fn yak_hash64(key: u64, mask: u64) -> u64 {
    let mut key = (!key).wrapping_add(key << 21) & mask;
    key ^= key >> 24;
    key = key.wrapping_add(key << 3).wrapping_add(key << 8) & mask;
    key ^= key >> 14;
    key = key.wrapping_add(key << 2).wrapping_add(key << 4) & mask;
    key ^= key >> 28;
    key.wrapping_add(key << 31) & mask
}

fn read_u32(input: &mut impl Read) -> Option<u32> {
    let mut buffer = [0u8; 4];
    input.read_exact(&mut buffer).ok()?;
    Some(u32::from_le_bytes(buffer))
}

fn read_u64(input: &mut impl Read) -> Option<u64> {
    let mut buffer = [0u8; 8];
    input.read_exact(&mut buffer).ok()?;
    Some(u64::from_le_bytes(buffer))
}

/// `yak_ch_restore_core` (`Trio.cpp` 66) in `YAK_LOAD_TRIOBIN1/2` mode:
/// merges one dump into the flat hash table, keeping only k-mers with
/// `count >= min_cnt` and labelling them 1 (or 2 for `count >= mid_cnt`) in
/// the class field selected by `shift` (0 paternal, 2 maternal).
///
/// The reference partitions by the low `pre` hash bits and stores
/// `(hash >> pre) << 10 | class`; the dump loses those `pre` bits, which are
/// recovered from the sub-table index here so the flat table can be keyed by
/// the full hash. A flat map is equivalent because the bucket and the stored
/// high bits together identify the hash.
fn restore_core(
    table: &mut U64HashTable,
    path: &str,
    shift: u32,
    min_cnt: u32,
    mid_cnt: u32,
) -> Option<(i32, i32)> {
    let mut input = std::io::BufReader::new(File::open(path).ok()?);
    let mut magic = [0u8; 4];
    if input.read_exact(&mut magic).is_err() {
        return None;
    }
    if magic != YAK_MAGIC {
        eprintln!("ERROR: wrong file magic.");
        return None;
    }
    let k = read_u32(&mut input)?;
    let pre = read_u32(&mut input)?;
    let counter_bits = read_u32(&mut input)?;
    if counter_bits != YAK_COUNTER_BITS {
        eprintln!(
            "ERROR: saved counter bits: {}; compile-time counter bits: {}",
            counter_bits, YAK_COUNTER_BITS
        );
        return None;
    }
    if pre < YAK_COUNTER_BITS {
        return None;
    }

    for index in 0..(1u32 << pre) {
        let _capacity = read_u32(&mut input)?;
        let size = read_u32(&mut input)?;
        for _ in 0..size {
            let key = read_u64(&mut input)?;
            let count = key & YAK_MAX_COUNT;
            let class = if count >= u64::from(mid_cnt) {
                2u64 << shift
            } else if count >= u64::from(min_cnt) {
                1u64 << shift
            } else {
                continue;
            };
            let hash = ((key >> YAK_COUNTER_BITS) << pre) | u64::from(index);
            let (bucket, absent) = table.put(hash);
            if absent {
                table.set_value_at(bucket, class);
            } else {
                let merged = table.value_at(bucket) | class;
                table.set_value_at(bucket, merged);
            }
        }
    }
    Some((k as i32, pre as i32))
}

/// `yak_ch_get`: the class bits of one hash, `-1` when absent.
fn yak_ch_get(table: &U64HashTable, hash: u64) -> i32 {
    match table.get(hash) {
        None => -1,
        Some(bucket) => (table.value_at(bucket) & 0xf) as i32,
    }
}

/// `tb_classify` (`Trio.cpp` 173): the haplotype verdict from the k-mer class
/// counts and the private-run spans.
fn tb_classify(spans: [i32; 2], counts: &[i32; 16], k: i32, ratio: f64) -> u8 {
    // `c[0<<2|2]` and `c[2<<2|0]`: paternal-only and maternal-only class 2.
    let paternal = counts[2];
    let maternal = counts[8];
    if spans[0] == 0 && spans[1] == 0 {
        if paternal == maternal {
            AMBIGU
        } else if paternal >= k - 4 + maternal
            && (maternal <= 1 || f64::from(paternal) * 0.05 > f64::from(maternal))
        {
            FATHER
        } else if maternal >= k - 4 + paternal
            && (paternal <= 1 || f64::from(maternal) * 0.05 > f64::from(paternal))
        {
            MOTHER
        } else {
            AMBIGU
        }
    } else if spans[0] > k && spans[1] > k {
        AMBIGU
    } else if spans[0] >= k - 4 + spans[1]
        && f64::from(spans[0]) * 0.05 >= f64::from(spans[1])
        && f64::from(paternal) * ratio > f64::from(maternal)
    {
        FATHER
    } else if spans[1] >= k - 4 + spans[0]
        && f64::from(spans[1]) * 0.05 >= f64::from(spans[0])
        && f64::from(maternal) * ratio > f64::from(paternal)
    {
        MOTHER
    } else {
        AMBIGU
    }
}

/// `tb_worker` (`Trio.cpp` 203): classifies one read by walking its k-mers.
fn classify_read(table: &U64HashTable, k: i32, sequence: &[u8]) -> u8 {
    let length = sequence.len();
    let mut states = vec![0u8; length];
    let mut counts = [0i32; 16];
    let (mask, shift) = if k < 32 {
        ((1u64 << (2 * k)) - 1, (2 * (k - 1)) as u32)
    } else {
        ((1u64 << k) - 1, (k - 1) as u32)
    };
    let mut x = [0u64; 4];
    let mut run = 0i32;
    for index in 0..length {
        let code = u64::from(SEQUENCE_NT4_TABLE[sequence[index] as usize]);
        if code < 4 {
            if k < 32 {
                x[0] = ((x[0] << 2) | code) & mask;
                x[1] = (x[1] >> 2) | ((3 - code) << shift);
            } else {
                x[0] = ((x[0] << 1) | (code & 1)) & mask;
                x[1] = ((x[1] << 1) | (code >> 1)) & mask;
                x[2] = (x[2] >> 1) | ((1 - (code & 1)) << shift);
                x[3] = (x[3] >> 1) | ((1 - (code >> 1)) << shift);
            }
            run += 1;
            if run >= k {
                let hash = if k < 32 {
                    yak_hash64(x[0].min(x[1]), mask)
                } else {
                    crate::counting::hash::hash_long_kmer(&x)
                };
                let flag = yak_ch_get(table, hash).max(0);
                counts[flag as usize] += 1;
                let paternal = flag & 3;
                let maternal = (flag >> 2) & 3;
                states[index] = if paternal == 2 && maternal == 0 {
                    1
                } else if maternal == 2 && paternal == 0 {
                    2
                } else {
                    0
                };
            }
        } else {
            run = 0;
            x = [0u64; 4];
        }
    }

    let mut spans = [0i32; 2];
    let mut start = 0usize;
    for end in 1..=length {
        if end == length || states[end] != states[start] {
            if states[start] > 0 && (end - start) as i32 >= k - 4 {
                spans[states[start] as usize - 1] += (end - start) as i32;
            }
            start = end;
        }
    }
    tb_classify(spans, &counts, k, 0.33)
}

/// `ha_triobin_yak` (`Trio.cpp` 261): loads both dumps and classifies every
/// read. A read's verdict overrides the external-list one.
fn ha_triobin_yak(options: &AssemblyOptions, reads: &mut ReadStore) {
    let paternal = options.trio_yak_file_names[0].as_deref().expect("checked");
    let maternal = options.trio_yak_file_names[1].as_deref().expect("checked");
    let min_cnt = options.minimum_binned_kmer_count.max(0) as u32;
    let mid_cnt = options.middle_binned_kmer_count.max(0) as u32;

    let mut table = U64HashTable::new(yak_ch_hash);
    let Some((k, _)) = restore_core(&mut table, paternal, 0, min_cnt, mid_cnt) else {
        return;
    };
    let Some((maternal_k, _)) = restore_core(&mut table, maternal, 2, min_cnt, mid_cnt) else {
        return;
    };
    assert_eq!(k, maternal_k, "yak dumps disagree on k");

    for read_id in 0..reads.total_reads as usize {
        let sequence = reads.read_sequence(read_id);
        reads.trio_labels[read_id] = classify_read(&table, k, &sequence);
    }
    eprintln!(
        "[M::ha_triobin_yak::{:.3}*{:.2}] ==> partitioned reads using yak dumps",
        logging::realtime_seconds(),
        logging::cpu_usage()
    );
}

/// `test_yak_binning` (`Trio.cpp` 352): whether the phase-list cache names
/// the same dump in its first line.
fn test_yak_binning(path: &str, command: &str) -> bool {
    let Some(raw) = read_text(path) else {
        return false;
    };
    let line_end = raw.iter().position(|&byte| byte == b'\n').unwrap_or(raw.len());
    let line = &raw[..line_end];
    let Some(tab) = line.iter().position(|&byte| byte == b'\t') else {
        return false;
    };
    tab + 1 < line.len() && &line[tab + 1..] == command.as_bytes()
}

/// `write_yak_binning` (`Overlaps.cpp` 23501): the `.hapN.phase.bin` name
/// lists, whose first line carries the dump path used to validate the cache.
pub fn write_yak_binning(options: &AssemblyOptions, reads: &ReadStore) {
    eprintln!("Writing binning to disk ...... ");
    // `%s` with a NULL argument prints "(null)" in glibc; list-only runs
    // reach here with no dump paths.
    let paternal_command = options.trio_yak_file_names[0].as_deref().unwrap_or("(null)");
    let maternal_command = options.trio_yak_file_names[1].as_deref().unwrap_or("(null)");

    for (haplotype, flag, command) in [
        (1, FATHER, paternal_command),
        (2, MOTHER, maternal_command),
    ] {
        let path = format!("{}.hap{}.phase.bin", options.output_file_name, haplotype);
        let result = File::create(&path).and_then(|file| {
            let mut output = BufWriter::new(file);
            let mut first = true;
            for read_id in 0..reads.total_reads as usize {
                if reads.trio_labels[read_id] != flag {
                    continue;
                }
                if first {
                    writeln!(output, "{}\t{}", reads.read_name(read_id), command)?;
                    first = false;
                } else {
                    writeln!(output, "{}", reads.read_name(read_id))?;
                }
            }
            output.flush()
        });
        if let Err(error) = result {
            eprintln!("[ERROR] failed to write {}: {}", path, error);
        }
    }
    eprintln!("Binning has been written.");
}

/// `load_yak_binning` (`Overlaps.cpp` 23534): replays a valid phase-list
/// cache. Returns whether both caches matched their dump paths.
pub fn load_yak_binning(options: &mut AssemblyOptions, reads: &mut ReadStore) -> bool {
    let paternal_list = format!("{}.hap1.phase.bin", options.output_file_name);
    let maternal_list = format!("{}.hap2.phase.bin", options.output_file_name);
    let paternal_dump = options.trio_yak_file_names[0].as_deref().unwrap_or("");
    let maternal_dump = options.trio_yak_file_names[1].as_deref().unwrap_or("");
    if !test_yak_binning(&paternal_list, paternal_dump)
        || !test_yak_binning(&maternal_list, maternal_dump)
    {
        return false;
    }

    // `ha_triobin` with the phase lists substituted for the dumps.
    let saved_lists = options.trio_list_file_names.clone();
    let saved_dumps = options.trio_yak_file_names.clone();
    options.trio_list_file_names = [Some(paternal_list), Some(maternal_list)];
    options.trio_yak_file_names = [None, None];
    run_trio_binning(options, reads);
    options.trio_list_file_names = saved_lists;
    options.trio_yak_file_names = saved_dumps;
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phrase_hstatus_reads_the_second_field_and_last_underscore() {
        let (name, hid) = phrase_hstatus(b"x_2\tread_A/ccs").unwrap();
        assert_eq!(name, b"read_A/ccs");
        assert_eq!(hid, 2);

        let (name, hid) = phrase_hstatus(b"foo_bar_3 read4 trailing").unwrap();
        assert_eq!(name, b"read4");
        assert_eq!(hid, 3);
    }

    #[test]
    fn phrase_hstatus_rejects_malformed_lines() {
        assert!(phrase_hstatus(b"single").is_none());
        assert!(phrase_hstatus(b"plain\tread").is_none());
        assert!(phrase_hstatus(b"x_1\t").is_none());
        assert!(phrase_hstatus(b"x_-1\tread").is_none());
        // A non-numeric suffix parses as `atoi` would: zero.
        assert_eq!(phrase_hstatus(b"no_underscore\tread").unwrap().1, 0);
    }

    #[test]
    fn update_poly_trio_follows_the_bitmask() {
        let hap_s = vec![0b01, 0b10, 0, 0b11];

        let mut flags = vec![DROP as u8, 0, 0, 0];
        update_poly_trio(1, &hap_s, &mut flags);
        assert_eq!(flags, vec![DROP as u8, MOTHER, AMBIGU, FATHER]);

        let mut flags = vec![DROP as u8, 0, 0, 0];
        update_poly_trio(2, &hap_s, &mut flags);
        assert_eq!(flags, vec![DROP as u8, FATHER, AMBIGU, FATHER]);
    }

    /// Writes one dump whose buckets are distributed by the low 10 hash bits,
    /// exactly like the files `yak count` produces.
    fn write_test_dump(path: &str, k: u32, hashes: &[u64]) {
        use std::io::Write;
        const PRE: u32 = 10;
        let mut buckets: Vec<Vec<u64>> = vec![Vec::new(); 1 << PRE];
        for &hash in hashes {
            let bucket = (hash & ((1 << PRE) - 1)) as usize;
            buckets[bucket].push(((hash >> PRE) << YAK_COUNTER_BITS) | 5);
        }
        let mut file = File::create(path).unwrap();
        file.write_all(&YAK_MAGIC).unwrap();
        file.write_all(&k.to_le_bytes()).unwrap();
        file.write_all(&PRE.to_le_bytes()).unwrap();
        file.write_all(&YAK_COUNTER_BITS.to_le_bytes()).unwrap();
        for keys in &buckets {
            file.write_all(&4u32.to_le_bytes()).unwrap();
            file.write_all(&(keys.len() as u32).to_le_bytes()).unwrap();
            for key in keys {
                file.write_all(&key.to_le_bytes()).unwrap();
            }
        }
    }

    fn read_hashes(sequence: &[u8], k: i32) -> Vec<u64> {
        let mask = (1u64 << (2 * k)) - 1;
        let shift = (2 * (k - 1)) as u32;
        let mut forward = 0u64;
        let mut reverse = 0u64;
        let mut run = 0i32;
        let mut hashes = Vec::new();
        for &base in sequence {
            let code = u64::from(SEQUENCE_NT4_TABLE[base as usize]);
            if code < 4 {
                forward = ((forward << 2) | code) & mask;
                reverse = (reverse >> 2) | ((3 - code) << shift);
                run += 1;
                if run >= k {
                    hashes.push(yak_hash64(forward.min(reverse), mask));
                }
            } else {
                run = 0;
                forward = 0;
                reverse = 0;
            }
        }
        hashes
    }

    #[test]
    fn yak_dump_round_trip_classifies_both_haplotypes() {
        let sequence = b"ACGTTGCAAGTCCTAGGATCCGATCGATGGCATTCAGGCATCGATCGTAGCTAGCATCGATCGGATC";
        let k = 15i32;
        let hashes = read_hashes(sequence, k);
        let directory = std::env::temp_dir().join(format!(
            "rshifiasm_yak_test_{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&directory).unwrap();
        let paternal = directory.join("paternal.yak");
        let maternal = directory.join("maternal.yak");

        write_test_dump(paternal.to_str().unwrap(), k as u32, &hashes);
        let mut table = U64HashTable::new(yak_ch_hash);
        let (stored_k, _) = restore_core(&mut table, paternal.to_str().unwrap(), 0, 2, 5)
            .expect("paternal dump must load");
        assert_eq!(stored_k, k);
        assert_eq!(classify_read(&table, k, sequence), FATHER);

        // The same k-mers in the maternal dump make every k-mer shared, which
        // leaves the read ambiguous.
        write_test_dump(maternal.to_str().unwrap(), k as u32, &hashes);
        restore_core(&mut table, maternal.to_str().unwrap(), 2, 2, 5)
            .expect("maternal dump must load");
        assert_eq!(classify_read(&table, k, sequence), AMBIGU);

        let _ = std::fs::remove_dir_all(&directory);
    }
}

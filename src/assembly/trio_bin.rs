//! `ha_triobin` (`Trio.cpp`): the trio binning front end.
//!
//! Two input modes fill [`ReadStore::trio_labels`]: read-name lists
//! (`-3`/`-4`) and yak k-mer dumps (`-1`/`-2`, a later stage). Both may run
//! in one invocation — the reference calls them in that order, so a read
//! flagged by both ends up with the yak verdict.

use std::collections::HashMap;
use std::io::Read;

use crate::command_line::options::{AssemblyOptions, TRIO_FATHER, TRIO_MOTHER, TRIO_UNKNOWN};
use crate::logging;
use crate::reads::read_store::ReadStore;

/// `ha_triobin` (`Trio.cpp` 450): resets every label and runs the requested
/// binning modes.
pub fn run_trio_binning(options: &AssemblyOptions, reads: &mut ReadStore) {
    reads.trio_labels.fill(TRIO_UNKNOWN);
    // `if (opt->fn_bin_list[0] && opt->fn_bin_list[1]) ha_triobin_list(opt);`
    if options.trio_list_file_names[0].is_some() && options.trio_list_file_names[1].is_some() {
        ha_triobin_list(options, reads);
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
    let raw = match std::fs::File::open(path).and_then(|mut file| {
        let mut raw = Vec::new();
        file.read_to_end(&mut raw)?;
        Ok(raw)
    }) {
        Ok(raw) => raw,
        Err(_) => {
            eprintln!("ERROR: failed to open file '{}'", path);
            // The reference returns -1 and the caller ignores it.
            return;
        }
    };
    // `kstream` over `gzopen`: decompress on the gzip magic, plain text
    // otherwise.
    let text = if raw.starts_with(&[0x1f, 0x8b]) {
        let mut plain = Vec::with_capacity(raw.len() * 4);
        flate2::read::MultiGzDecoder::new(raw.as_slice())
            .read_to_end(&mut plain)
            .ok()
            .map(|_| plain)
            .unwrap_or(raw)
    } else {
        raw
    };

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

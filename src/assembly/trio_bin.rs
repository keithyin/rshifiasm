//! `ha_triobin` (`Trio.cpp`): the trio binning front end.
//!
//! Two input modes fill [`ReadStore::trio_labels`]: read-name lists
//! (`-3`/`-4`) and yak k-mer dumps (`-1`/`-2`, a later stage). Both may run
//! in one invocation — the reference calls them in that order, so a read
//! flagged by both ends up with the yak verdict.

use std::collections::HashMap;
use std::io::Read;

use crate::assembly::trio_clean::DROP;
use crate::assembly::unitig::{AMBIGU, FATHER, MOTHER};
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
}

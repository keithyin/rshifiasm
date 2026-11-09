//! The assembly pipeline, the counterpart of `Assembly.cpp`.
//!
//! Only the top level orchestration lives here; every stage is implemented in
//! the module that owns its algorithm.

pub mod asg;
pub mod bp_graph;
pub mod bubble;
pub mod bubble_identify;
pub mod chimeric;
pub mod corrected_output;
pub mod extract;
pub mod gfa_clean;
pub mod graph_build;
pub mod hamming;
pub mod hap_cov;
pub mod hap_alignment;
pub mod hit;
pub mod maxcut;
pub mod pipeline;
pub mod probe;
pub mod purge_dups;
pub mod rescue;
pub mod tangle;
pub mod trans_clean;
pub mod bubble_chain;
pub mod trans_cov;
pub mod trio_dual;
pub mod trio_clean;
pub mod trio_bin;
pub mod break_contig;
pub mod read_cache;
pub mod unitig;

use crate::command_line::options::{AssemblyOptions, FLAG_NO_KMER_FILTER};
use crate::correction::ec_driver::run_error_correction;
use crate::logging;
use crate::reads::ReadStore;

/// `ha_assemble`: the default HiFi-only path.
pub fn assemble(options: &mut AssemblyOptions) -> i32 {
    // `ha_assemble`: with the default `load_index_from_disk` the restart
    // caches replace the whole correction stage (`Assembly.cpp` 2062).
    let loaded = if options.load_index_from_disk != 0 {
        read_cache::load_all_data_from_disk(options)
    } else {
        None
    };

    let (read_store, state, write_caches) = match loaded {
        Some((store, state)) => {
            eprintln!(
                "[M::ha_assemble::{:.3}*{:.2}] ==> loaded corrected reads and overlaps from disk",
                logging::realtime_seconds(),
                logging::cpu_usage()
            );
            // `if (asm_opt.extract_list) { ha_extract_print_list(...); exit(0); }`
            // (`Assembly.cpp` 2065), before any other output of the loaded run.
            if options.extract_list_file_name.is_some() {
                if let Err(error) = extract::print_requested_list(options, &store, &state) {
                    eprintln!("[ERROR] failed to print the extracted overlaps: {}", error);
                    return 1;
                }
                return 0;
            }
            // `if (asm_opt.het_cov == -1024) hap_recalculate_peaks(...)`.
            if options.heterozygous_coverage == -1024 {
                read_cache::hap_recalculate_peaks(options, &store);
            }
            // The loaded run still honours `--write-ec` / `--write-paf`
            // (`Assembly.cpp` 2070-2073).
            if let Err(error) = corrected_output::write_requested_outputs(
                options,
                &store,
                &state.forward,
                false,
            ) {
                eprintln!("[ERROR] failed to write the requested outputs: {}", error);
                return 1;
            }
            (store, state, false)
        }
        None => {
            let mut read_store = match ReadStore::load_from_files(options) {
                Ok(store) => store,
                Err(error) => {
                    eprintln!("[ERROR] failed to load the input reads: {}", error);
                    return 1;
                }
            };

            // High occurrence k-mers are filtered out before any overlap is
            // computed; the same pass infers the homozygous coverage.
            let filter_table = if options.flags & FLAG_NO_KMER_FILTER == 0 {
                let (filter_table, homozygous_coverage) =
                    crate::counting::filter_table::generate_high_occurrence_filter(
                        options,
                        &read_store,
                    );
                options.update_coverage(homozygous_coverage);
                Some(filter_table)
            } else {
                None
            };

            // `ha_ec`: correct every read before the assembly graph is built,
            // replacing the stored sequence with the corrected one. The call
            // also runs `ha_ec_ff`, which rebuilds the index and recomputes
            // every read's stored overlaps for the final round.
            let corrected =
                match run_error_correction(options, &mut read_store, filter_table.as_ref()) {
                    Ok(corrected) => corrected,
                    Err(error) => {
                        eprintln!("[ERROR] error correction failed: {}", error);
                        return 1;
                    }
                };
            eprintln!(
                "[M::ha_assemble::{:.3}*{:.2}@{:.3}GB] ==> found overlaps for the final round",
                logging::realtime_seconds(),
                logging::cpu_usage(),
                logging::peak_resident_memory_gigabytes()
            );

            if let Err(error) = corrected_output::write_requested_outputs(
                options,
                &read_store,
                &corrected.state.forward,
                false,
            ) {
                eprintln!("[ERROR] failed to write the requested outputs: {}", error);
                return 1;
            }
            // `Output_PAF(); ha_triobin(&asm_opt);` — binning runs on the
            // freshly computed overlaps only; a run that loaded the caches
            // replays it from the phase checkpoints instead (`Trio.cpp`,
            // `load_yak_binning`).
            if options.trio_binning_requested() {
                trio_bin::run_trio_binning(options, &mut read_store);
            }
            (read_store, corrected.state, true)
        }
    };

    // `ha_opt_update_cov_min`, `build_string_graph_without_clean` and
    // `clean_graph`: build, clean and write the assembly graphs. A run that
    // loaded the caches passes `write = 0` (`Assembly.cpp` 2190).
    if let Err(error) =
        pipeline::build_string_graph_without_clean(options, &read_store, state, write_caches)
    {
        eprintln!("[ERROR] {}", error);
        return 1;
    }

    0
}

/// `ha_assemble_pair`: the secondary-input variant, which loads two stores.
pub fn assemble_pair(options: &mut AssemblyOptions) -> i32 {
    eprintln!("[ERROR] ha_assemble_pair is not implemented yet");
    let _ = options;
    1
}

/// `ha_ec_dbg`: overlap calculation debug entry point.
pub fn error_correction_debug(options: &mut AssemblyOptions) -> i32 {
    match crate::correction::ec_driver::run_debug_overlap(options) {
        Ok(()) => 0,
        Err(error) => {
            eprintln!("[ERROR] {}", error);
            1
        }
    }
}

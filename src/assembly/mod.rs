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
pub mod ul;
pub mod ul_raw_chain;
pub mod ul_scall;
pub mod ul_index;
pub mod ul_consensus;
pub mod ul_load;
pub mod ul_store;
pub mod ul_src;
pub mod ul_graph;
pub mod ul_rescall;
pub mod dual_scaf;
pub mod ctg_trans;
pub mod scaf_tangle;
pub mod trans_base;

use crate::command_line::options::{AssemblyOptions, FLAG_NO_KMER_FILTER};
use crate::correction::ec_driver::{run_error_correction, run_single_correction_round};
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
                options.single_cell_mode != 0,
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

            // `-b`: `ha_ec` dumps the requested reads and exits after the
            // first correction round. The debug read ids are never filled in
            // this path, so all three files stay empty.
            if let Some(name_file) = options.required_read_name.clone() {
                eprintln!("Queried debugging reads at: {}", name_file);
                for suffix in [".debug.stdout", ".debug.r0.fa", ".debug.r1.fa"] {
                    let path = format!("{}{}", name_file, suffix);
                    if let Err(error) = std::fs::File::create(&path) {
                        eprintln!("[ERROR] failed to create {}: {}", path, error);
                        return 1;
                    }
                }
                eprintln!("Print debugging information to: {}.debug.stdout", name_file);
                eprintln!("Print raw reads to: {}.debug.r0.fa", name_file);
                eprintln!("Print corrected reads to: {}.debug.r1.fa", name_file);
                if let Err(error) =
                    run_single_correction_round(options, &mut read_store, filter_table.as_ref())
                {
                    eprintln!("[ERROR] error correction failed: {}", error);
                    return 1;
                }
                return 0;
            }

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
                options.single_cell_mode != 0,
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

/// `ha_assemble_pair` (`Assembly.cpp` 2128): the `--sec-in` variant, which
/// merges the two numbered caches `<prefix>.ec{0,1}.bin` and runs one legacy
/// overlap-and-correct pass before building the graph.
pub fn assemble_pair(options: &mut AssemblyOptions) -> i32 {
    let loaded = if options.load_index_from_disk != 0 {
        read_cache::load_all_data_from_disk(options)
    } else {
        None
    };

    let (read_store, state, write_caches) = match loaded {
        Some((store, state)) => {
            eprintln!(
                "[M::ha_assemble_pair::{:.3}*{:.2}] ==> loaded corrected reads and overlaps from disk",
                logging::realtime_seconds(),
                logging::cpu_usage()
            );
            if options.extract_list_file_name.is_some() {
                if let Err(error) = extract::print_requested_list(options, &store, &state) {
                    eprintln!("[ERROR] failed to print the extracted overlaps: {}", error);
                    return 1;
                }
                return 0;
            }
            if options.heterozygous_coverage == -1024 {
                read_cache::hap_recalculate_peaks(options, &store);
            }
            if let Err(error) = corrected_output::write_requested_outputs(
                options,
                &store,
                &state.forward,
                options.single_cell_mode != 0,
            ) {
                eprintln!("[ERROR] failed to write the requested outputs: {}", error);
                return 1;
            }
            (store, state, false)
        }
        None => {
            let mut read_store = ReadStore::default();
            let prefix = options.output_file_name.clone();
            let mut initial_flags: Vec<(bool, u8)> = Vec::new();
            for identifier in 0..2 {
                let Some(mut flags) =
                    read_cache::append_all_reads(&mut read_store, options, &prefix, identifier)
                else {
                    eprintln!(
                        "[M::ha_assemble_pair::] Cannot load {}.{}",
                        prefix, identifier
                    );
                    return 1;
                };
                initial_flags.append(&mut flags);
            }
            eprintln!("[M::ha_assemble_pair] size of buffer: 0.000GB");

            // `ha_ft_gen(..., read_from_store=1)`; the same pass updates the
            // chain count from the inferred homozygous coverage.
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

            // `r = asm_opt.number_of_round - 1; ha_opt_reset_to_round(&asm_opt, r);`
            options.reset_to_round(options.correction_round_count - 1);
            let state = match crate::correction::pair::run_overlap_correction(
                options,
                &read_store,
                filter_table.as_ref(),
                &initial_flags,
            ) {
                Ok(state) => state,
                Err(error) => {
                    eprintln!("[ERROR] the pair-mode overlap pass failed: {}", error);
                    return 1;
                }
            };
            eprintln!(
                "[M::ha_assemble_pair::{:.3}*{:.2}@{:.3}GB] ==> found overlaps for the final round",
                logging::realtime_seconds(),
                logging::cpu_usage(),
                logging::peak_resident_memory_gigabytes()
            );
            drop(filter_table);

            if let Err(error) = corrected_output::write_requested_outputs(
                options,
                &read_store,
                &state.forward,
                options.single_cell_mode != 0,
            ) {
                eprintln!("[ERROR] failed to write the requested outputs: {}", error);
                return 1;
            }
            // `Output_PAF(); ha_triobin(&asm_opt);`
            if options.trio_binning_requested() {
                trio_bin::run_trio_binning(options, &mut read_store);
            }
            (read_store, state, true)
        }
    };

    if let Err(error) =
        pipeline::build_string_graph_without_clean(options, &read_store, state, write_caches)
    {
        eprintln!("[ERROR] {}", error);
        return 1;
    }

    0
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

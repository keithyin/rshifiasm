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
pub mod gfa_clean;
pub mod graph_build;
pub mod hamming;
pub mod hap_cov;
pub mod hit;
pub mod pipeline;
pub mod probe;
pub mod purge_dups;
pub mod rescue;
pub mod tangle;
pub mod trans_clean;
pub mod trans_cov;
pub mod trio_clean;
pub mod unitig;

use crate::command_line::options::{AssemblyOptions, FLAG_NO_KMER_FILTER};
use crate::correction::ec_driver::run_error_correction;
use crate::logging;
use crate::reads::ReadStore;

/// `ha_assemble`: the default HiFi-only path.
pub fn assemble(options: &mut AssemblyOptions) -> i32 {
    let mut read_store = match ReadStore::load_from_files(options) {
        Ok(store) => store,
        Err(error) => {
            eprintln!("[ERROR] failed to load the input reads: {}", error);
            return 1;
        }
    };

    // High occurrence k-mers are filtered out before any overlap is computed;
    // the same pass infers the homozygous coverage.
    let filter_table = if options.flags & FLAG_NO_KMER_FILTER == 0 {
        let (filter_table, homozygous_coverage) =
            crate::counting::filter_table::generate_high_occurrence_filter(options, &read_store);
        options.update_coverage(homozygous_coverage);
        Some(filter_table)
    } else {
        None
    };

    // `ha_ec`: correct every read before the assembly graph is built,
    // replacing the stored sequence with the corrected one. The call also
    // runs `ha_ec_ff`, which rebuilds the index and recomputes every read's
    // stored overlaps for the final round.
    let corrected = match run_error_correction(options, &mut read_store, filter_table.as_ref()) {
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

    // `ha_opt_update_cov_min`, `build_string_graph_without_clean` and
    // `clean_graph`: build, clean and write the assembly graphs.
    if let Err(error) =
        pipeline::build_string_graph_without_clean(options, &read_store, corrected.state)
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
    eprintln!("[ERROR] ha_ec_dbg is not implemented yet");
    let _ = options;
    1
}

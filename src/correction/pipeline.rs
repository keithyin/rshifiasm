//! The read-level error-correction pass, the counterpart of `ha_ec` in
//! `Assembly.cpp`.
//!
//! For every read in input order the pass computes its candidate overlaps from
//! the shared position index, refines them and rebuilds the sequence, then
//! writes the corrected sequence back into the read store. The reference runs
//! several correction rounds; a HiFi assembly uses one, so the round parameter
//! is kept but always `0` for now.

use crate::command_line::options::AssemblyOptions;
use crate::correction::correct::correct_overlap;
use crate::correction::read_writeback::write_corrected_read_back;
use crate::counting::filter_table::KmerFilterTable;
use crate::counting::position_index::PositionIndex;
use crate::logging;
use crate::overlap::engine::{CandidateParameters, compute_candidate_overlaps};
use crate::overlap::region::ChainData;
use crate::reads::ReadStore;

/// Band width of the chaining step at the correction call site, `bw_thres`.
const BAND_WIDTH_THRESHOLD: f64 = 0.02;

/// Statistics of one correction pass.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ErrorCorrectionStats {
    /// Reads visited.
    pub total_reads: u64,
    /// Reads whose correction covered every window.
    pub fully_covered_reads: u64,
    /// Reads whose sequence changed.
    pub corrected_reads: u64,
    /// Bases changed over all reads.
    pub corrected_bases: u64,
}

/// `ha_ec`: corrects every read of the store in place.
///
/// Reads are visited in identifier order and the corrected sequence replaces
/// the stored one, which is what the downstream assembly stages consume and
/// what `--write-ec` prints.
pub fn correct_all_reads(
    options: &AssemblyOptions,
    reads: &mut ReadStore,
    index: &PositionIndex,
    filter_table: Option<&KmerFilterTable>,
) -> ErrorCorrectionStats {
    let read_lengths = reads.read_lengths.clone();
    let mut chain_data = ChainData::default();
    let mut stats = ErrorCorrectionStats::default();
    let total_reads = reads.total_reads as usize;

    // Debug aids: `RSHIFIASM_EC_STOP_AT` limits the pass to the first reads so
    // that one read can be inspected quickly, `RSHIFIASM_EC_READ` prints the
    // window layout of a single read.
    let stop_after = std::env::var("RSHIFIASM_EC_STOP_AT")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(total_reads)
        .min(total_reads);
    let traced_read = std::env::var("RSHIFIASM_EC_READ")
        .ok()
        .and_then(|value| value.parse::<usize>().ok());

    for read_identifier in 0..stop_after {
        let sequence = reads.read_sequence(read_identifier);
        if sequence.is_empty() {
            continue;
        }

        let parameters = CandidateParameters {
            minimizer_window: options.minimizer_window_size,
            kmer_length: options.k_mer_length,
            sample_distance: options.minimizer_sample_distance,
            rewindow: options.minimizer_rewindow_size,
            use_homopolymer_compression: options.flags & crate::command_line::options::FLAG_NO_HPC
                == 0,
            homozygous_coverage: options.homozygous_coverage,
            maximum_chain_count: options.maximum_chain_count,
            band_width_threshold: BAND_WIDTH_THRESHOLD,
            read_lengths: &read_lengths,
            index,
            filter_table,
        };

        let mut overlaps =
            compute_candidate_overlaps(&sequence, read_identifier as u64, &parameters, &mut chain_data);

        // The reference skips the refinement when no chain survived.
        if overlaps.is_empty() {
            stats.total_reads += 1;
            continue;
        }

        let tracing_read = traced_read == Some(read_identifier);
        if tracing_read {
            for (index, region) in overlaps.iter().enumerate() {
                eprintln!(
                    "[M::ha_ec::ovl {}] y={} x=[{},{}] shared={} ovl_len={}",
                    index,
                    region.y_identifier,
                    region.x_position_start,
                    region.x_position_end,
                    region.shared_seed,
                    region.overlap_length
                );
            }
        }
        crate::correction::window_consensus::set_traced_read(tracing_read);
        let outcome = correct_overlap(
            &mut overlaps,
            reads,
            &sequence,
            options.maximum_overlap_divergence_error_correction,
            options.maximum_overlap_divergence_final,
        );
        if tracing_read {
            for (index, region) in overlaps.iter().enumerate() {
                let mut ranges: Vec<(i32, i32, i32)> = Vec::new();
                for window in &region.window_cigar.windows {
                    match ranges.last_mut() {
                        Some(last) if last.1 + 1 == window.x_start => {
                            last.1 = window.x_end;
                        }
                        _ => ranges.push((window.x_start, window.x_end, window.y_end)),
                    }
                }
                eprintln!(
                    "[M::ha_ec::cov {}] match={} align={} x=[{},{}] ranges={:?}",
                    index,
                    region.is_match,
                    region.alignment_length,
                    region.x_position_start,
                    region.x_position_end,
                    ranges
                );
            }
            eprintln!(
                "[M::ha_ec::read {}] input={} corrected={} changed={} overlaps={}",
                read_identifier,
                sequence.len(),
                outcome.sequence.len(),
                outcome.corrected_base,
                overlaps.len()
            );
        }
        crate::correction::window_consensus::set_traced_read(false);

        stats.total_reads += 1;
        if outcome.fully_covered {
            stats.fully_covered_reads += 1;
        }
        if outcome.sequence.len() != sequence.len() || outcome.corrected_base != 0 {
            stats.corrected_reads += 1;
            stats.corrected_bases += outcome.corrected_base;
        }
        write_corrected_read_back(reads, read_identifier, &outcome.sequence);
    }

    eprintln!(
        "[M::ha_ec::{:.3}*{:.2}@ {:.3}GB] ==> corrected {} of {} reads, {} bases",
        logging::realtime_seconds(),
        logging::cpu_usage(),
        logging::peak_resident_memory_gigabytes(),
        stats.corrected_reads,
        stats.total_reads,
        stats.corrected_bases
    );

    stats
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::counting::kmer_counter::KmerCountTable;
    use crate::counting::position_index::PositionIndex;

    fn options() -> AssemblyOptions {
        AssemblyOptions::default()
    }

    #[test]
    fn an_empty_index_leaves_every_read_unchanged() {
        let mut options = options();
        options.homozygous_coverage = 10;
        let mut reads = ReadStore::default();
        reads.push_read("read_a", b"ACGTACGTACGT");
        reads.push_read("read_b", b"TTTTAAAACCCC");

        let counts = KmerCountTable::new(options.k_mer_length, 0, 0);
        let index = PositionIndex::build(&counts, options.k_mer_length);
        let stats = correct_all_reads(&options, &mut reads, &index, None);

        assert_eq!(stats.total_reads, 2);
        assert_eq!(stats.corrected_reads, 0);
        assert_eq!(reads.read_sequence(0), b"ACGTACGTACGT");
        assert_eq!(reads.read_sequence(1), b"TTTTAAAACCCC");
    }
}

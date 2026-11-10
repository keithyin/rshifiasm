//! The legacy overlap-and-correct pass of `ha_overlap_cal` (`Assembly.cpp`
//! 1142), used only by the `--sec-in` pair mode.
//!
//! The normal pipeline drives the newer `ha_ec` / `worker_hap_ec` chain; the
//! pair path stores no overlaps of its own and instead recomputes every read's
//! candidates from the position index, refines them with the legacy
//! `correct_overlap` driver, and writes the survivors back into the stored
//! PAF lists through `push_final_overlaps`.

use crate::command_line::options::{AssemblyOptions, FLAG_NO_HPC};
use crate::correction::correct::correct_overlap;
use crate::correction::ec_driver::{CorrectionState, PafList};
use crate::counting::filter_table::KmerFilterTable;
use crate::counting::minimizer_count::count_minimizers;
use crate::overlap::engine::{compute_candidate_overlaps, CandidateParameters};
use crate::overlap::final_overlap::push_final_overlaps;
use crate::overlap::region::ChainData;
use crate::reads::read_store::ReadStore;

/// `ha_overlap_cal(round, read_from_store=1)` with the legacy worker.
///
/// The position index is built from the already-loaded reads, which is the
/// `ha_pt_gen(..., read_from_store=1, ...)` mode; the coverage peaks it
/// infers overwrite the pair-mode sum, exactly like the reference.
pub fn run_overlap_correction(
    options: &mut AssemblyOptions,
    reads: &ReadStore,
    filter_table: Option<&KmerFilterTable>,
    initial_flags: &[(bool, u8)],
) -> Result<CorrectionState, String> {
    let counts = count_minimizers(options, reads, filter_table);
    let index = counts
        .index
        .ok_or_else(|| "the pair-mode position index is missing".to_string())?;
    // `ha_overlap_cal` only assigns the freshly inferred peaks after the worker
    // pass; the worker itself still sees the coverage from `ha_ft_gen`.
    let worker_homozygous_coverage = options.homozygous_coverage;
    if filter_table.is_none() {
        options.update_coverage(counts.homozygous_coverage);
    }

    let parameters = CandidateParameters {
        minimizer_window: options.minimizer_window_size,
        kmer_length: options.k_mer_length,
        sample_distance: options.minimizer_sample_distance,
        rewindow: options.minimizer_rewindow_size,
        use_homopolymer_compression: options.flags & FLAG_NO_HPC == 0,
        homozygous_coverage: worker_homozygous_coverage,
        maximum_chain_count: options.maximum_chain_count,
        // `ha_get_candidates_interface(..., 0.02, ...)`.
        band_width_threshold: 0.02,
        read_lengths: &reads.read_lengths,
        index: &index,
        filter_table,
    };
    let maximum_overlap_divergence =
        f64::from(options.maximum_overlap_divergence_error_correction);
    let final_error_rate = f64::from(options.maximum_overlap_divergence_final);

    let mut state = CorrectionState::new(reads.total_reads as usize);
    // `worker_ovec_cal0` reuses one overlap buffer per thread: fields the
    // append does not write keep the previous read's values at that slot.
    let mut slot_template: Vec<crate::overlap::region::OverlapRegion> = Vec::new();
    for read_identifier in 0..reads.total_reads as usize {
        let query = reads.read_sequence(read_identifier);
        let mut chain_data = ChainData::default();
        let mut overlaps = compute_candidate_overlaps(
            &query,
            read_identifier as u64,
            &parameters,
            &mut chain_data,
            Some(&mut slot_template),
        );
        // `is_consensus = 0` in `worker_ovec_cal0`.
        let _ = correct_overlap(
            &mut overlaps,
            reads,
            &query,
            maximum_overlap_divergence,
            final_error_rate,
            false,
            options.informative_coverage,
            options.single_haplotype_coverage,
        );
        let (forward, reverse, _) = push_final_overlaps(&overlaps, reads);
        // `clear_ma_hit_t_alloc` only clears the list; the flags loaded from
        // the `.ovlpN.source.bin` cache stay in place.
        let (is_fully_corrected, is_abnormal) = initial_flags
            .get(read_identifier)
            .copied()
            .unwrap_or((false, 0));
        state.forward[read_identifier] = PafList {
            hits: forward,
            is_fully_corrected,
            is_abnormal,
        };
        state.reverse[read_identifier] = PafList {
            hits: reverse,
            is_fully_corrected: false,
            is_abnormal: 0,
        };
        // The reference keeps the buffer's capacity across reads: slots past
        // the current length retain values from whichever earlier read filled
        // them, so the template is grown but never shrunk.
        if slot_template.len() < overlaps.len() {
            slot_template.resize(overlaps.len(), Default::default());
        }
        slot_template[..overlaps.len()].clone_from_slice(&overlaps);
    }

    // `asm_opt.hom_cov = hom_cov; asm_opt.het_cov = het_cov;` at the end of
    // `ha_overlap_cal`.
    options.homozygous_coverage = counts.homozygous_coverage;
    options.heterozygous_coverage = counts.heterozygous_coverage;

    Ok(state)
}

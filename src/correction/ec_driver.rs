//! The error-correction driver, the counterpart of `ha_ec`, `cal_ec_r`,
//! `cal_ec_multiple`, `cal_update_ec_multiple`, `sl_ec_r` and
//! `worker_hap_post_rev` in `ecovlp.cpp`.
//!
//! One round corrects every read: chaining against the round's position
//! index, refinement, haplotype voting and the windowed consensus produce the
//! read's consensus cigar (`scc`), which is written back into the read store;
//! the stored overlaps are then remeasured against the corrected read. Except
//! for the last round, reads and overlaps are reverse-complemented afterwards
//! so that the next round sees them in their other orientation, which the
//! reference re-flips on the following round.
//!
//! Only the HiFi branch is implemented. The ONT and single-cell branches
//! change almost every call site of the correction pass, so they are rejected
//! instead of silently producing different reads.

use std::io::Write;

use crate::command_line::options::{AssemblyOptions, FLAG_FAST, FLAG_NO_HPC};
use crate::correction::consensus_gfa::ConsensusScratch;
use crate::correction::dedup_chains::dedup_chains;
use crate::correction::haplotype_evidence::{HaplotypeEvidenceAlloc, WINDOW_MAXIMUM_SIZE};
use crate::correction::hc_align::BitExtz;
use crate::correction::hc_r_alin::gen_hc_r_alin_ea;
use crate::correction::read_writeback::{apply_consensus_to_read, write_corrected_read_back};
use crate::correction::rphase::{rphase_hc, RphaseScratch};
use crate::counting::filter_table::KmerFilterTable;
use crate::counting::minimizer::sketch;
use crate::counting::minimizer_count::count_minimizers;
use crate::counting::position_index::PositionIndex;
use crate::logging;
use crate::overlap::final_overlap::run_final_overlap_pass;
use crate::overlap::lchain_parameters::set_lchain_parameters;
use crate::overlap::paf::PafHit;
use crate::overlap::qgen_chain::{lchain_qgen_mcopy_fast, QgenChainParameters};
use crate::overlap::qgen_seeds::collect_qgen_seeds;
use crate::overlap::region::{CandidatesList, OverlapRegion};
use crate::overlap::stored_overlap::{
    check_well_cal, push_ne_ovlp, push_nec_re_from_window, quick_exact_match,
};
use crate::reads::read_store::{ReadStore, TRIO_LABEL_UNKNOWN};

/// `HA_KMER_GOOD_RATIO`, the share of the homozygous coverage below which a
/// k-mer count is considered low.
const KMER_GOOD_RATIO: f64 = 0.333;

/// `WINDOW_HC`: window length of the HiFi refinement.
const ALIGNMENT_WINDOW_LENGTH: i64 = 775;

/// `WINDOW_OHC`: window length of the ONT refinement.
const ONT_ALIGNMENT_WINDOW_LENGTH: i64 = 375;

/// `WINDOW_HC_FAST`: window length of the consensus vote.
const CONSENSUS_WINDOW_LENGTH: u64 = 512;

/// `WINDOW_MAX_SIZE`: window length of the haplotype pass.
const PHASING_WINDOW_LENGTH: i64 = WINDOW_MAXIMUM_SIZE as i64;

/// `E_KHIT`: k-mer hit threshold handed to the refinement.
const KMER_HIT_THRESHOLD: i64 = 31;

/// `COV_W`: window over which the copy-number filter measures coverage.
const COPY_NUMBER_COVERAGE_WINDOW: u64 = 3072;

/// Copy-number multiplier of the chaining filter, `mcopy_num`.
const COPY_NUMBER: i64 = 3;

/// Copy-number rate of the chaining filter, `mcopy_rate`.
const COPY_NUMBER_RATE: f64 = 0.7;

/// Chaining cutoff in alignment length, `chain_cutoff`.
const CHAIN_CUTOFF: u32 = 2;

/// K-mer hit cutoff of the copy-number filter, `mcopy_khit_cut`.
const COPY_NUMBER_KMER_HIT_CUTOFF: u32 = 32;

/// Band width of the chaining step, the HiFi branch of `bw_thres`.
const BAND_WIDTH_THRESHOLD: f64 = 0.02;

/// Band width of the ONT branch of `bw_thres`.
const ONT_BAND_WIDTH_THRESHOLD: f64 = 0.05;

/// `align_gap_rate` of the ONT refinement.
const ONT_ALIGNMENT_GAP_RATE: f64 = 0.006;

/// `align_gap_max` of the ONT refinement.
const ONT_ALIGNMENT_GAP_MAX: i64 = 64;

/// `HPC_PL`: homopolymer compression strength of the ONT phasing pass.
const ONT_HPC_LENGTH: u64 = 12;

/// `MIN_COVERAGE_THRESHOLD * 2`: depth every stretch of a read must reach to
/// be fully corrected.
const MINIMUM_CORRECTION_DEPTH: i64 = 6;

/// `occ_tot` of the consensus vote.
const CONSENSUS_VOTE_TOTAL_OCCURRENCE: u64 = 3;

/// `occ_exact` of the consensus vote.
const CONSENSUS_VOTE_EXACT_OCCURRENCE: f64 = 0.500001;

/// `cns_g_wl`: window length used to split long consensus gaps.
const CONSENSUS_WINDOW: u64 = 256;

/// `occ_thres` of the haplotype pass.
const PHASING_OCCURRENCE_THRESHOLD: u8 = 1;

/// One read's stored overlap list, `ma_hit_t_alloc` plus the two flags the
/// correction pipeline keeps next to it.
#[derive(Debug, Clone, Default)]
pub struct PafList {
    /// The stored overlaps, `buffer[0..length]`.
    pub hits: Vec<PafHit>,
    /// `is_fully_corrected`: every stretch of the read had enough coverage.
    pub is_fully_corrected: bool,
    /// `is_abnormal`: 1 when coverage drops inside the read, 2 when it drops
    /// at an end.
    pub is_abnormal: u8,
}

impl AsRef<[PafHit]> for PafList {
    /// The overlap list itself, for the writers that take `R_INF.paf`.
    fn as_ref(&self) -> &[PafHit] {
        &self.hits
    }
}

/// The per-read overlap lists the correction passes read and rewrite,
/// `R_INF.paf` and `R_INF.reverse_paf`.
#[derive(Debug, Clone, Default)]
pub struct CorrectionState {
    /// Overlaps in the read's current orientation, `paf`.
    pub forward: Vec<PafList>,
    /// Overlaps to reads on the other strand, `reverse_paf`.
    pub reverse: Vec<PafList>,
}

impl CorrectionState {
    /// Empty lists for `total_reads` reads.
    pub fn new(total_reads: usize) -> Self {
        Self {
            forward: vec![PafList::default(); total_reads],
            reverse: vec![PafList::default(); total_reads],
        }
    }
}

/// What the correction rounds leave behind for the following stages.
pub struct ErrorCorrectionResult {
    /// The overlap lists after the last round.
    pub state: CorrectionState,
    /// Bases seen in the last round, `tot_b`.
    pub total_bases: u64,
    /// Bases covered by the last round's consensus, `tot_e`.
    pub corrected_bases: u64,
    /// Always empty: the final overlap pass rebuilds and then destroys the
    /// position index, exactly like `ha_ec_ff`.
    pub index: Option<PositionIndex>,
}

/// `ha_ec`: runs every correction round and leaves the reads corrected in the
/// store.
///
/// `filter_table` is the high-occurrence filter, `None` when the run disabled
/// it; the reference distinguishes the two when it picks the shrink cutoff and
/// when it refreshes the coverage estimate after the first round.
pub fn run_error_correction(
    options: &mut AssemblyOptions,
    reads: &mut ReadStore,
    filter_table: Option<&KmerFilterTable>,
) -> Result<ErrorCorrectionResult, &'static str> {
    if options.single_cell_mode != 0 && options.ont_mode == 0 {
        return Err("the single-cell correction branch is not implemented yet");
    }
    if options.partial_correction_round_count != 0 {
        return Err("partial correction rounds are not implemented yet");
    }

    let number_of_rounds = options.correction_round_count.max(1) as usize;
    let total_reads = reads.total_reads as usize;
    let mut state = CorrectionState::new(total_reads);
    let mut consensus_cigars: Vec<Vec<u16>> = vec![Vec::new(); total_reads];
    let mut total_bases = 0u64;
    let mut corrected_bases = 0u64;

    for round in 0..number_of_rounds {
        options.reset_to_round(round as i32);
        let counts = count_minimizers(options, reads, filter_table);
        let homozygous_coverage = counts.homozygous_coverage;
        let heterozygous_coverage = counts.heterozygous_coverage;
        options.homozygous_coverage = homozygous_coverage;
        options.heterozygous_coverage = heterozygous_coverage;
        // `ha_opt_update_cov` only runs when the high-occurrence filter was
        // skipped; otherwise `ha_ft_gen` already updated the chain count.
        if round == 0 && filter_table.is_none() {
            options.update_coverage(homozygous_coverage);
        }
        let index = counts
            .index
            .ok_or("the correction rounds require the normal counting mode")?;

        let round_stats = correct_round(
            options,
            reads,
            &mut state,
            &mut consensus_cigars,
            &index,
            filter_table,
            round,
            number_of_rounds,
        );
        total_bases = round_stats.0;
        corrected_bases = round_stats.1;

        eprintln!(
            "[M::ha_assemble::{:.3}*{:.2}@{:.3}GB] ==> corrected reads for round {}",
            logging::realtime_seconds(),
            logging::cpu_usage(),
            logging::peak_resident_memory_gigabytes(),
            round + 1
        );
        eprintln!(
            "[M::ha_assemble] # bases: {}; # corrected bases: {}",
            total_bases, corrected_bases
        );

        // `des_idx`: every round but the last drops its index and rebuilds it
        // from the freshly corrected reads; the final overlap pass destroys
        // the last one as well.
    }

    // `ha_assemble` points the options at the round after the last one before
    // it computes the final overlaps.
    options.reset_to_round(number_of_rounds as i32);

    // `ha_ec_ff(1)`: rebuild the index from the corrected reads and recompute
    // every read's stored overlaps. The reference destroys the index at the
    // end of the pass, so `index` in the result is always empty.
    run_final_overlap_pass(options, reads, &mut state, filter_table)?;

    Ok(ErrorCorrectionResult {
        state,
        total_bases,
        corrected_bases,
        index: None,
    })
}

/// One correction round only, the `-b` debugging entry of `ha_ec`.
///
/// `ha_ec` exits right after `cal_ec_r` of round zero when a debug read list
/// was supplied, so the final overlap pass never runs.
pub fn run_single_correction_round(
    options: &mut AssemblyOptions,
    reads: &mut ReadStore,
    filter_table: Option<&KmerFilterTable>,
) -> Result<(), &'static str> {
    let total_reads = reads.total_reads as usize;
    let mut state = CorrectionState::new(total_reads);
    let mut consensus_cigars: Vec<Vec<u16>> = vec![Vec::new(); total_reads];

    options.reset_to_round(0);
    let counts = count_minimizers(options, reads, filter_table);
    options.homozygous_coverage = counts.homozygous_coverage;
    options.heterozygous_coverage = counts.heterozygous_coverage;
    if filter_table.is_none() {
        options.update_coverage(counts.homozygous_coverage);
    }
    let index = counts
        .index
        .ok_or("the correction rounds require the normal counting mode")?;
    let number_of_rounds = options.correction_round_count.max(1) as usize;
    correct_round(
        options,
        reads,
        &mut state,
        &mut consensus_cigars,
        &index,
        filter_table,
        0,
        number_of_rounds,
    );
    Ok(())
}

/// `cal_ec_r`: one correction round, from the worker pass to the orientation
/// flip.
#[allow(clippy::too_many_arguments)]
fn correct_round(
    options: &AssemblyOptions,
    reads: &mut ReadStore,
    state: &mut CorrectionState,
    consensus_cigars: &mut [Vec<u16>],
    index: &PositionIndex,
    filter_table: Option<&KmerFilterTable>,
    round: usize,
    number_of_rounds: usize,
) -> (u64, u64) {
    // `cal_ec_multiple`.
    let profiling = std::env::var_os("RSHIFIASM_EC_PROFILE").is_some();
    // Debug aid: `RSHIFIASM_EC_STOP_AT=<n>` limits the round to the first `n`
    // reads, so a single read's stages can be inspected without waiting for a
    // full pass. Never set during a golden run.
    let corrected_read_count = std::env::var("RSHIFIASM_EC_STOP_AT")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .map_or(reads.total_reads as usize, |limit| {
            limit.min(reads.total_reads as usize)
        });
    let started = logging::realtime_seconds();
    let mut total_bases = 0u64;
    let mut corrected_bases = 0u64;
    // `if(round == asm_opt.number_of_round-1 && asm_opt.is_dbg_het_cnt)
    // CALLOC(het_cnt, R_INF.total_reads);` — the live `cal_ec_r` worker never
    // fills the array (only the dead `worker_ovec*` drivers do), so the log
    // is all zeros, exactly like the reference.
    let is_final_round = round + 1 == number_of_rounds;
    for read_identifier in 0..corrected_read_count {
        let (bases, corrected) = correct_read(
            options,
            reads,
            state,
            consensus_cigars,
            read_identifier,
            index,
            filter_table,
            profiling,
        );
        total_bases += bases;
        corrected_bases += corrected;
    }
    eprintln!(
        "[M::pec::{:.3}] # bases: {}; # corrected bases: {}",
        logging::realtime_seconds() - started,
        total_bases,
        corrected_bases
    );
    if options.debug_heterozygosity_count != 0 && is_final_round {
        write_het_cnt_log(options, reads);
    }

    // `sl_ec_r`: the consensus cigar replaces the stored sequence.
    for read_identifier in 0..corrected_read_count {
        apply_consensus_to_read(reads, read_identifier, &consensus_cigars[read_identifier]);
    }

    // `cal_update_ec_multiple`.
    let started = logging::realtime_seconds();
    let (exact, non_exact) = update_exact_overlaps(state, reads, consensus_cigars);
    eprintln!(
        "[M::pec::{:.3}] # exact o: {}; # non-exact o: {}",
        logging::realtime_seconds() - started,
        exact,
        non_exact
    );

    // The reference flips both the reads and their overlap lists except on
    // the last round of an even number of rounds; `is_cr` alternates so that
    // two consecutive flips restore the original orientation.
    let is_reverse_round = round % 2 == 1;
    if !is_final_round || (is_final_round && is_reverse_round) {
        reverse_complement_round(state, reads);
    }

    (total_bases, corrected_bases)
}

/// `worker_hap_ec` for one read: chains, refines, votes and writes the read's
/// consensus cigar and overlap lists back into `state`.
///
/// Returns the read's length (the `# bases` contribution) and the number of
/// bases the consensus covered (the `# corrected bases` contribution).
#[allow(clippy::too_many_arguments)]
#[allow(unused_assignments)] // the last profile mark resets a timer nobody reads
fn correct_read(
    options: &AssemblyOptions,
    reads: &mut ReadStore,
    state: &mut CorrectionState,
    consensus_cigars: &mut [Vec<u16>],
    read_identifier: usize,
    index: &PositionIndex,
    filter_table: Option<&KmerFilterTable>,
    profiling: bool,
) -> (u64, u64) {
    let query = reads.read_sequence(read_identifier);
    let query_length = query.len() as u64;

    // `RSHIFIASM_EC_PROFILE=1` reports how long each stage of the per-read
    // correction takes.
    let mut stage_start = std::time::Instant::now();
    macro_rules! stage {
        ($name:literal) => {
            if profiling {
                eprintln!(
                    "[M::ec_profile] rid={} {}={:.3}ms",
                    read_identifier,
                    $name,
                    stage_start.elapsed().as_secs_f64() * 1000.0
                );
                stage_start = std::time::Instant::now();
            }
        };
    }

    // `worker_hap_ec` recomputes these from the round's coverage estimate.
    let high_occurrence = (f64::from(options.homozygous_coverage) * (2.0 - KMER_GOOD_RATIO)) as u32;
    let low_occurrence = (f64::from(options.homozygous_coverage) * KMER_GOOD_RATIO) as u32;
    let is_ont = options.ont_mode != 0;

    // `h_ec_lchain`: the HiFi parameters are the accurate chaining set.
    let chaining = set_lchain_parameters(true, options.k_mer_length as u32);
    let sketch_result = sketch(
        &query,
        options.minimizer_window_size,
        options.k_mer_length,
        read_identifier as u32,
        options.flags & FLAG_NO_HPC == 0,
        filter_table,
        options.minimizer_sample_distance,
        options.minimizer_rewindow_size,
    );
    let seeds = collect_qgen_seeds(
        &sketch_result.minimizers,
        index,
        &reads.read_lengths,
        high_occurrence,
        low_occurrence,
    );
    stage!("seeds");
    let chained = lchain_qgen_mcopy_fast(
        seeds,
        read_identifier as u32,
        query_length,
        &reads.read_lengths,
        &QgenChainParameters {
            maximum_chain_count: options.maximum_chain_count as u32,
            maximum_skip: chaining.maximum_skip,
            maximum_iterations: chaining.maximum_iterations,
            maximum_distance: chaining.maximum_distance,
            gap_penalty: chaining.gap_penalty,
            skip_penalty: chaining.skip_penalty,
            band_width_rate: if is_ont {
                ONT_BAND_WIDTH_THRESHOLD
            } else {
                BAND_WIDTH_THRESHOLD
            },
            quick_check: chaining.quick_check != 0,
            append_both_ends: true,
            generate_offset: true,
            copy_number: COPY_NUMBER,
            copy_rate: COPY_NUMBER_RATE,
            chain_cutoff: CHAIN_CUTOFF,
            copy_kmer_cutoff: COPY_NUMBER_KMER_HIT_CUTOFF,
            coverage_window: COPY_NUMBER_COVERAGE_WINDOW,
        },
    );
    stage!("chain");

    // `b->olist` and `b->clist`; the refinement appends single-hit candidates
    // at `length`, exactly like `gen_single_khit`.
    let mut candidates = CandidatesList::default();
    candidates.list = chained.candidates;
    candidates.length = chained.candidate_length as i64;
    let mut overlaps = chained.overlaps;

    let mut auxiliary = OverlapRegion::default();
    let mut extz = BitExtz::new();
    let mut exact_order: Vec<u64> = Vec::new();
    let mut cells: Vec<u16> = Vec::new();
    let kept = gen_hc_r_alin_ea(
        &mut overlaps,
        &mut candidates,
        reads,
        &query,
        &mut extz,
        &mut auxiliary,
        &state.forward[read_identifier].hits,
        options.maximum_overlap_divergence_error_correction,
        if is_ont {
            ONT_ALIGNMENT_WINDOW_LENGTH
        } else {
            ALIGNMENT_WINDOW_LENGTH
        },
        KMER_HIT_THRESHOLD,
        &mut exact_order,
        &mut cells,
        is_ont,
        if is_ont { ONT_ALIGNMENT_GAP_RATE } else { -1.0 },
        if is_ont { ONT_ALIGNMENT_GAP_MAX } else { -1 },
    );
    // The reference compacts `ol->length = k` in place; the port returns the
    // kept count instead.
    overlaps.truncate(kept);
    stage!("refine");

    let mut haplotype = HaplotypeEvidenceAlloc::new();
    let mut phasing_scratch = RphaseScratch::default();
    let mut mapped_overlaps_length = 0i64;
    rphase_hc(
        &mut overlaps,
        reads,
        &mut haplotype,
        &query,
        &mut phasing_scratch,
        0,
        PHASING_WINDOW_LENGTH,
        PHASING_OCCURRENCE_THRESHOLD,
        if is_ont { ONT_HPC_LENGTH } else { 0 },
        if is_ont { 1 } else { 0 },
        if is_ont { 1 } else { 0 },
        &mut mapped_overlaps_length,
        options.informative_coverage,
        options.single_haplotype_coverage,
        read_identifier as u64,
        i64::from(options.heterozygous_coverage),
        i64::from(options.homozygous_coverage),
        i64::from(options.polyploidy),
        if is_ont { Some(reads) } else { None },
        is_ont,
    );
    stage!("phase");

    dedup_chains(&mut overlaps);
    stage!("dedup");

    let mut consensus = ConsensusScratch::new();
    let corrected_bases = consensus.wcns_gen(
        &mut overlaps,
        reads,
        &query,
        &mut extz,
        0,
        CONSENSUS_WINDOW_LENGTH,
        query_length as i64,
        CONSENSUS_VOTE_TOTAL_OCCURRENCE,
        CONSENSUS_VOTE_EXACT_OCCURRENCE,
        CONSENSUS_WINDOW,
    );
    stage!("consensus");

    let consensus_cigar = push_nec_re_from_window(&consensus.auxiliary);
    let hits = push_ne_ovlp(&overlaps, 1, reads, Some(&consensus_cigar));
    let reverse_hits = push_ne_ovlp(&overlaps, 2, reads, None);
    state.reverse[read_identifier] = PafList {
        hits: reverse_hits,
        ..PafList::default()
    };

    let mut is_fully_corrected = false;
    let mut is_abnormal = 0u8;
    let mut indices: Vec<u64> = Vec::new();
    check_well_cal(
        &consensus_cigar,
        &mut indices,
        &mut is_fully_corrected,
        &mut is_abnormal,
        query_length as i64,
        MINIMUM_CORRECTION_DEPTH,
        &hits,
    );

    consensus_cigars[read_identifier] = consensus_cigar;
    state.forward[read_identifier] = PafList {
        hits,
        is_fully_corrected,
        is_abnormal,
    };
    reads.trio_labels[read_identifier] = TRIO_LABEL_UNKNOWN;

    (query_length, corrected_bases)
}

/// `print_het_cnt_log`: writes `<output>.het_cnt.log` with one
/// `>name<TAB>count` line per read. The counts are all zero because the live
/// correction worker never records them.
fn write_het_cnt_log(options: &AssemblyOptions, reads: &ReadStore) {
    use std::io::Write;
    let path = format!("{}.het_cnt.log", options.output_file_name);
    let file = match std::fs::File::create(&path) {
        Ok(file) => file,
        Err(error) => {
            eprintln!("[ERROR] failed to write {}: {}", path, error);
            return;
        }
    };
    let mut output = std::io::BufWriter::new(file);
    eprintln!(
        "[M::print_het_cnt_log::] ==> print cnt of het sites to {}...",
        path
    );
    for name in &reads.read_names {
        if let Err(error) = writeln!(output, ">{}\t{}", name, 0) {
            eprintln!("[ERROR] failed to write {}: {}", path, error);
            return;
        }
    }
    if let Err(error) = output.flush() {
        eprintln!("[ERROR] failed to write {}: {}", path, error);
    }
}

/// `cal_update_ec_multiple`: re-checks every stored overlap against the
/// corrected reads and updates its `el` flag.
fn update_exact_overlaps(
    state: &mut CorrectionState,
    reads: &ReadStore,
    consensus_cigars: &[Vec<u16>],
) -> (u64, u64) {
    let mut exact = 0u64;
    let mut non_exact = 0u64;
    for read_identifier in 0..reads.total_reads as usize {
        // `if(!(R_INF.paf[i].length)) return;`
        if state.forward[read_identifier].hits.is_empty() {
            continue;
        }
        let query = reads.read_sequence(read_identifier);
        for hit in state.forward[read_identifier].hits.iter_mut() {
            if hit.exact_region && quick_exact_match(hit, reads, &query, consensus_cigars) {
                hit.exact_region = true;
                hit.exact_region_byte = 1;
                exact += 1;
            } else {
                hit.exact_region = false;
                non_exact += 1;
            }
        }
    }
    (exact, non_exact)
}

/// `flip_paf_rc`: mirrors one overlap list into the reverse-complemented read.
///
/// Intervals are clamped to the current read lengths first — the correction
/// may have shortened a read — and overlaps that collapse are dropped while
/// the surviving order is preserved.
fn flip_paf_rc(read_identifier: usize, hits: &mut Vec<PafHit>, reads: &ReadStore) {
    let query_length = reads.read_length(read_identifier) as i64;
    let mut kept = 0usize;
    for index in 0..hits.len() {
        let hit = hits[index].clone();
        let target_length = reads.read_length(hit.target_identifier as usize) as i64;
        let query_start = i64::from(hit.query_start).clamp(0, query_length);
        let query_end = i64::from(hit.query_end).clamp(0, query_length);
        let target_start = i64::from(hit.target_start).clamp(0, target_length);
        let target_end = i64::from(hit.target_end).clamp(0, target_length);
        if query_end > query_start && target_end > target_start {
            let mut flipped = hit;
            flipped.query_start = (query_length - query_end) as u32;
            flipped.query_end = (query_length - query_start) as u32;
            flipped.target_start = (target_length - target_end) as u32;
            flipped.target_end = (target_length - target_start) as u32;
            hits[kept] = flipped;
            kept += 1;
        }
    }
    hits.truncate(kept);
}

/// `worker_hap_post_rev`: flips both overlap lists and reverse-complements the
/// read in place, so the next round reads it from the other end.
fn reverse_complement_round(state: &mut CorrectionState, reads: &mut ReadStore) {
    for read_identifier in 0..reads.total_reads as usize {
        flip_paf_rc(
            read_identifier,
            &mut state.forward[read_identifier].hits,
            reads,
        );
        flip_paf_rc(
            read_identifier,
            &mut state.reverse[read_identifier].hits,
            reads,
        );
        let sequence = reads.read_sequence_reverse_complement(read_identifier);
        // `worker_hap_post_rev` also reverses the 2-bit quality track.
        let quality = if reads.qualities[read_identifier].is_empty() {
            None
        } else {
            let mut levels = reads.read_quality(
                read_identifier,
                0,
                reads.read_lengths[read_identifier] as i64,
                false,
            );
            levels.reverse();
            Some(crate::reads::read_store::pack_quality_levels(&levels))
        };
        write_corrected_read_back(reads, read_identifier, &sequence);
        if let Some(quality) = quality {
            reads.qualities[read_identifier] = quality;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(query_start: u32, query_end: u32, target_start: u32, target_end: u32) -> PafHit {
        PafHit {
            query_identifier: 0,
            query_start,
            query_end,
            target_identifier: 1,
            target_start,
            target_end,
            matching_length: 1,
            block_length: 100,
            reverse: false,
            without_large_indel: true,
            exact_region: true,
            exact_region_byte: 0,
            deleted: false,
        }
    }

    #[test]
    fn flipping_mirrors_both_intervals() {
        let mut reads = ReadStore::default();
        reads.push_read("query", &vec![b'A'; 100]);
        reads.push_read("target", &vec![b'C'; 80]);
        let mut hits = vec![hit(10, 30, 5, 20)];

        flip_paf_rc(0, &mut hits, &reads);
        assert_eq!(hits[0].query_start, 70);
        assert_eq!(hits[0].query_end, 90);
        assert_eq!(hits[0].target_start, 60);
        assert_eq!(hits[0].target_end, 75);

        // A second flip restores the original intervals.
        flip_paf_rc(0, &mut hits, &reads);
        assert_eq!(hits[0].query_start, 10);
        assert_eq!(hits[0].query_end, 30);
        assert_eq!(hits[0].target_start, 5);
        assert_eq!(hits[0].target_end, 20);
    }

    #[test]
    fn flipping_drops_overlaps_that_collapse_after_shortening() {
        let mut reads = ReadStore::default();
        reads.push_read("query", &vec![b'A'; 20]);
        reads.push_read("target", &vec![b'C'; 10]);
        // The stored coordinates still point at the pre-correction lengths.
        let mut hits = vec![hit(30, 40, 20, 30), hit(0, 5, 0, 5)];

        flip_paf_rc(0, &mut hits, &reads);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].query_start, 15);
        assert_eq!(hits[0].query_end, 20);
        assert_eq!(hits[0].target_start, 5);
        assert_eq!(hits[0].target_end, 10);
    }
}

/// `ha_ec_dbg` + `cal_ec_r_dbg` + `worker_hap_ec_dbg_paf`
/// (`Assembly.cpp` 1053, `ecovlp.cpp` 6311–6375): dumps every read's aligned
/// overlap windows, cigars included, as PAF — the pre-consensus checkpoint
/// that isolates C++/Rust differences to chaining and alignment.
///
/// The index is built without a k-mer filter table; the reference processes
/// reads in chunks of 2000 through a pipeline, which only bounds memory —
/// the records are written in read order either way.
pub fn run_debug_overlap(options: &mut AssemblyOptions) -> Result<(), String> {
    let mut reads = ReadStore::load_from_files(options)
        .map_err(|error| format!("failed to load the input reads: {}", error))?;

    // `ha_pt_gen(&asm_opt, 0, 0, 0, ...)`: no filter table.
    let counts = count_minimizers(options, &reads, None);
    options.homozygous_coverage = counts.homozygous_coverage;
    options.heterozygous_coverage = counts.heterozygous_coverage;
    options.update_coverage(options.homozygous_coverage);
    let index = counts
        .index
        .ok_or_else(|| "the debug-overlap index is missing".to_string())?;

    let high_occurrence =
        (f64::from(options.homozygous_coverage) * (2.0 - KMER_GOOD_RATIO)) as u32;
    let low_occurrence = (f64::from(options.homozygous_coverage) * KMER_GOOD_RATIO) as u32;

    let paf_name = format!("{}.ovlp.paf", options.output_file_name);
    let mut output = std::io::BufWriter::new(
        std::fs::File::create(&paf_name)
            .map_err(|error| format!("failed to create {}: {}", paf_name, error))?,
    );

    let chaining = set_lchain_parameters(true, options.k_mer_length as u32);
    let chain_parameters = QgenChainParameters {
        maximum_chain_count: options.maximum_chain_count as u32,
        maximum_skip: chaining.maximum_skip,
        maximum_iterations: chaining.maximum_iterations,
        maximum_distance: chaining.maximum_distance,
        gap_penalty: chaining.gap_penalty,
        skip_penalty: chaining.skip_penalty,
        band_width_rate: BAND_WIDTH_THRESHOLD,
        quick_check: chaining.quick_check != 0,
        append_both_ends: true,
        generate_offset: true,
        copy_number: COPY_NUMBER,
        copy_rate: COPY_NUMBER_RATE,
        chain_cutoff: CHAIN_CUTOFF,
        copy_kmer_cutoff: COPY_NUMBER_KMER_HIT_CUTOFF,
        coverage_window: COPY_NUMBER_COVERAGE_WINDOW,
    };

    for read_identifier in 0..reads.total_reads as usize {
        let query = reads.read_sequence(read_identifier);
        let query_length = query.len() as u64;

        let sketch_result = sketch(
            &query,
            options.minimizer_window_size,
            options.k_mer_length,
            read_identifier as u32,
            options.flags & FLAG_NO_HPC == 0,
            None,
            options.minimizer_sample_distance,
            options.minimizer_rewindow_size,
        );
        let seeds = collect_qgen_seeds(
            &sketch_result.minimizers,
            &index,
            &reads.read_lengths,
            high_occurrence,
            low_occurrence,
        );
        let chained = lchain_qgen_mcopy_fast(
            seeds,
            read_identifier as u32,
            query_length,
            &reads.read_lengths,
            &chain_parameters,
        );

        let mut candidates = CandidatesList::default();
        candidates.list = chained.candidates;
        candidates.length = chained.candidate_length as i64;
        let mut overlaps = chained.overlaps;
        let mut auxiliary = OverlapRegion::default();
        let mut extz = BitExtz::new();
        let mut exact_order: Vec<u64> = Vec::new();
        let mut scratch: Vec<u16> = Vec::new();
        // The dbg call passes an empty `R_INF.paf[i]` and
        // `chem_drop = 0, align_gap_rate = -1, align_gap_max = -1`.
        let kept = gen_hc_r_alin_ea(
            &mut overlaps,
            &mut candidates,
            &reads,
            &query,
            &mut extz,
            &mut auxiliary,
            &[],
            options.maximum_overlap_divergence_error_correction,
            ALIGNMENT_WINDOW_LENGTH,
            KMER_HIT_THRESHOLD,
            &mut exact_order,
            &mut scratch,
            false,
            -1.0,
            -1,
        );
        overlaps.truncate(kept);

        // The record loop of `worker_hap_ec_dbg_paf`.
        for region in &overlaps {
            if region.window_cigar.windows.is_empty() {
                continue;
            }
            let target_length = reads.read_length(region.y_identifier as usize);
            for window in &region.window_cigar.windows {
                // `is_ualn_win(a)`: error == INT16_MAX && clen == 0 && extra_end < 0.
                if window.error == i16::MAX && window.cigar_length == 0 && window.extra_end < 0 {
                    continue;
                }

                let mut target_start = window.y_start as u64;
                let mut target_end = (window.y_end + 1) as u64;
                if region.y_strand != 0 {
                    target_start = target_length - (window.y_end as u64) - 1;
                    target_end = target_length - window.y_start as u64;
                }

                // `pop_trace` over the window's cells, op 0/1 → M, 2 → I, 3 → D.
                let mut cigar = String::new();
                let mut cell_index = window.cigar_start as usize;
                let cell_end = cell_index + window.cigar_length as usize;
                while cell_index < cell_end {
                    let operation = (region.window_cigar.cells[cell_index] >> 14) as usize;
                    let mut length =
                        u32::from(region.window_cigar.cells[cell_index] & 0x3fff);
                    cell_index += 1;
                    while cell_index < cell_end
                        && (region.window_cigar.cells[cell_index] >> 14) as usize == operation
                    {
                        length += u32::from(region.window_cigar.cells[cell_index] & 0x3fff);
                        cell_index += 1;
                    }
                    let operation_character = ['M', 'M', 'I', 'D'][operation];
                    cigar.push_str(&format!("{}{}", length, operation_character));
                }

                // The reference prints the TARGET (y) read first and the
                // QUERY (x) read second, with the strand character taken from
                // `"+-"[y_pos_strand]`.
                writeln!(
                    output,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\tcg:Z:{}",
                    reads.read_name(region.y_identifier as usize),
                    reads.read_length(region.y_identifier as usize),
                    target_start,
                    target_end,
                    ['+', '-'][region.y_strand as usize],
                    reads.read_name(region.x_identifier as usize),
                    reads.read_length(region.x_identifier as usize),
                    window.x_start,
                    window.x_end + 1,
                    cigar,
                )
                .map_err(|error| format!("failed to write {}: {}", paf_name, error))?;
            }
        }
    }

    output
        .flush()
        .map_err(|error| format!("failed to write {}: {}", paf_name, error))?;
    Ok(())
}

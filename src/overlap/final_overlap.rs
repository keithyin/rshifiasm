//! The final overlap pass, `ha_ec_ff` in `Assembly.cpp` and
//! `h_ec_lchain_fast_new`/`push_ff_ovlp` in `ecovlp.cpp`.
//!
//! Once the correction rounds have rewritten every read, the reads are
//! chained one last time against an index rebuilt from the corrected
//! sequences. The band width of this pass is much narrower than the rounds'
//! (`0.001` against `0.02`), so the chains and the overlaps stored by the last
//! round describe nearly the same alignments. `h_ec_lchain_fast_new` reconciles
//! the two:
//!
//! * a chain that reproduces a stored overlap inherits the stored `strong`
//!   and `exact` flags and survives;
//! * a chain that matches no stored overlap survives only when it is exactly
//!   equal to the target substring;
//! * a stored overlap flagged exact that no chain reproduced is appended as an
//!   overlap in its own right.
//!
//! Duplicated chains pointing at the same target read are collapsed to the one
//! with the fewest errors, and the survivors are split by strand into `paf`
//! and `reverse_paf` — the overlap lists the assembly graph is built from.

use crate::command_line::options::{AssemblyOptions, FLAG_NO_HPC};
use crate::correction::cal_exz::{push_alignment_window, set_exact_exz};
use crate::correction::ec_driver::CorrectionState;
use crate::correction::hc_align::BitExtz;
use crate::correction::hc_r_alin::exact_ec_check;
use crate::counting::filter_table::KmerFilterTable;
use crate::counting::minimizer::sketch;
use crate::counting::minimizer_count::count_minimizers;
use crate::counting::position_index::PositionIndex;
use crate::container::ksort::radix_sort_u64;
use crate::logging;
use crate::overlap::lchain_parameters::set_lchain_parameters;
use crate::overlap::paf::PafHit;
use crate::overlap::qgen_chain::{lchain_qgen_mcopy_fast, QgenChainParameters};
use crate::overlap::qgen_seeds::collect_qgen_seeds;
use crate::overlap::region::{sort_by_y_identifier, OverlapRegion};
use crate::reads::read_store::ReadStore;

/// `HA_KMER_GOOD_RATIO`, the share of the homozygous coverage below which a
/// k-mer count is considered low.
const KMER_GOOD_RATIO: f64 = 0.333;

/// `bw_thres` of the final chaining pass. The correction rounds chain at
/// `0.02`; the final pass only needs to rediscover the alignments the last
/// round already stored, so it uses a much narrower band.
const BAND_WIDTH_RATE: f64 = 0.001;

/// `sh`: the share of a stored overlap a chain must reproduce to match it.
const SHARED_SEED_THRESHOLD: f64 = 0.866666;

/// Copy-number multiplier of the chaining filter, `mcopy_num`.
const COPY_NUMBER: i64 = 3;

/// Copy-number rate of the chaining filter, `mcopy_rate`.
const COPY_NUMBER_RATE: f64 = 0.7;

/// Chaining cutoff in alignment length, `chain_cutoff`.
const CHAIN_CUTOFF: u32 = 2;

/// K-mer hit cutoff of the copy-number filter, `mcopy_khit_cut`.
const COPY_NUMBER_KMER_HIT_CUTOFF: u32 = 32;

/// `COV_W`: window over which the copy-number filter measures coverage.
const COPY_NUMBER_COVERAGE_WINDOW: u64 = 3072;

/// One worker's overlap counters, `cnt[0..6]` of `push_ff_ovlp`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct FinalOverlapCounters {
    /// `cnt[0]`: overlaps in the read's own orientation.
    pub forward: u64,
    /// `cnt[1]`: overlaps pointing at the other strand.
    pub reverse: u64,
    /// `cnt[2]`: forward overlaps with `ml == 1`.
    pub strong: u64,
    /// `cnt[3]`: forward overlaps with `ml == 0`.
    pub weak: u64,
    /// `cnt[4]`: forward overlaps flagged exact.
    pub exact: u64,
    /// `cnt[5]`: forward overlaps without a large indel.
    pub without_large_indel: u64,
}

impl FinalOverlapCounters {
    fn add_assign(&mut self, other: &Self) {
        self.forward += other.forward;
        self.reverse += other.reverse;
        self.strong += other.strong;
        self.weak += other.weak;
        self.exact += other.exact;
        self.without_large_indel += other.without_large_indel;
    }
}

/// `ha_ec_ff(1)`: rebuilds the position index from the corrected reads, runs
/// the final overlap pass over every read and prints the overlap statistics.
///
/// The reference destroys the index it built when the pass ends, so nothing is
/// returned; `state` carries the resulting overlap lists.
pub fn run_final_overlap_pass(
    options: &mut AssemblyOptions,
    reads: &ReadStore,
    state: &mut CorrectionState,
    filter_table: Option<&KmerFilterTable>,
) -> Result<(), &'static str> {
    // `ha_pt_gen(..., read_from_store = 1, ...)`: unlike the rounds, the
    // index is always rebuilt, because the last round has rewritten the reads
    // and possibly reverse-complemented them.
    let counts = count_minimizers(options, reads, filter_table);
    options.homozygous_coverage = counts.homozygous_coverage;
    options.heterozygous_coverage = counts.heterozygous_coverage;
    let index = counts
        .index
        .ok_or("the final overlap pass requires the normal counting mode")?;

    // `ha_print_ovlp_stat_0`.
    let started = logging::realtime_seconds();
    let mut counters = FinalOverlapCounters::default();
    for read_identifier in 0..reads.total_reads as usize {
        let read_counters =
            final_overlap_for_read(options, reads, state, read_identifier, &index, filter_table);
        counters.add_assign(&read_counters);
    }

    eprintln!(
        "[M::ha_print_ovlp_stat_0] # overlaps: {}",
        counters.forward
    );
    eprintln!(
        "[M::ha_print_ovlp_stat_0] # strong overlaps: {}",
        counters.strong
    );
    eprintln!(
        "[M::ha_print_ovlp_stat_0] # weak overlaps: {}",
        counters.weak
    );
    eprintln!(
        "[M::ha_print_ovlp_stat_0] # exact overlaps: {}",
        counters.exact
    );
    eprintln!(
        "[M::ha_print_ovlp_stat_0] # inexact overlaps: {}",
        counters.forward - counters.exact
    );
    eprintln!(
        "[M::ha_print_ovlp_stat_0] # overlaps without large indels: {}",
        counters.without_large_indel
    );
    eprintln!(
        "[M::ha_print_ovlp_stat_0] # reverse overlaps: {}",
        counters.reverse
    );
    eprintln!(
        "[M::ha_print_ovlp_stat_0] # running time: {:.3}",
        logging::realtime_seconds() - started
    );

    // `ha_ec_ff` destroys the index before returning (`ha_pt_destroy`).
    drop(index);
    Ok(())
}

/// `worker_hap_dc_ec_gen_new_idx`: chains one corrected read and rebuilds its
/// two stored overlap lists.
fn final_overlap_for_read(
    options: &AssemblyOptions,
    reads: &ReadStore,
    state: &mut CorrectionState,
    read_identifier: usize,
    index: &PositionIndex,
    filter_table: Option<&KmerFilterTable>,
) -> FinalOverlapCounters {
    let query = reads.read_sequence(read_identifier);
    let query_length = query.len() as u64;

    let high_occurrence = (f64::from(options.homozygous_coverage) * (2.0 - KMER_GOOD_RATIO)) as u32;
    let low_occurrence = (f64::from(options.homozygous_coverage) * KMER_GOOD_RATIO) as u32;

    // `h_ec_lchain`: the accurate chaining set, but with the narrow band width
    // of the final pass.
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
        Some(high_occurrence),
        Some(low_occurrence),
    );
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
            band_width_rate: BAND_WIDTH_RATE,
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
    let mut overlaps = chained.overlaps;

    // `overlap_region_sort_y_id` before `h_ec_lchain_fast_new`: the chain list
    // arrives ordered by query start, and the refinement walks it in target
    // order with a monotone cursor over the stored overlaps.
    sort_by_y_identifier(&mut overlaps);

    // The forward list is both an input and the destination of this pass; the
    // reference mutates it in place, clearing the `el` flags of the stored
    // overlaps a chain reproduced.
    let mut stored_forward = std::mem::take(&mut state.forward[read_identifier].hits);
    let mut extz = BitExtz::new();
    refine_chained_overlaps(
        &mut overlaps,
        &mut stored_forward,
        &state.reverse[read_identifier].hits,
        reads,
        &query,
        read_identifier as u32,
        &mut extz,
    );

    // `if(is_ont && is_uncorrected_read(..., 1600)) olist.length = 0`: ONT
    // reads with an uncovered stretch of at least 1600 bp keep no overlaps.
    if options.ont_mode != 0
        && is_uncorrected_read(&overlaps, query_length as i64, 1600)
    {
        overlaps.clear();
    }

    let (forward_hits, reverse_hits, counters) = push_final_overlaps(&overlaps, reads);
    state.forward[read_identifier].hits = forward_hits;
    state.reverse[read_identifier].hits = reverse_hits;
    counters
}

/// `is_uncorrected_read`: whether the chain runs leave an uncovered stretch
/// of at least `minimum_length` inside the read.
fn is_uncorrected_read(overlaps: &[OverlapRegion], read_length: i64, minimum_length: i64) -> bool {
    let mut events: Vec<u64> = Vec::with_capacity(overlaps.len() * 2);
    for overlap in overlaps {
        events.push(u64::from(overlap.x_position_start) << 1);
        events.push((u64::from(overlap.x_position_end) + 1) << 1 | 1);
    }
    radix_sort_u64(&mut events);

    let mut depth = 0i64;
    let mut start = 0i64;
    let mut end = 0i64;
    for &event in &events {
        let previous = depth;
        if event & 1 != 0 {
            depth -= 1;
        } else {
            depth += 1;
        }
        end = (event >> 1) as i64;
        if end > start && previous == 0 && end - start >= minimum_length {
            return true;
        }
        start = end;
    }
    end = read_length;
    if end > start && depth == 0 && (end - start >= minimum_length || end - start >= read_length) {
        return true;
    }
    false
}

/// `h_ec_lchain_fast_new`: reconciles the chains against the overlaps the last
/// round stored and rewrites `overlaps` in place.
///
/// `stored_forward` is `in0`, the read's own overlap list; it is mutated
/// because reproduced exact overlaps lose their `el` flag and are therefore
/// not appended again. `stored_reverse` is `in1`, the opposite-strand list,
/// which the reference only reads.
#[allow(clippy::too_many_arguments)]
fn refine_chained_overlaps(
    overlaps: &mut Vec<OverlapRegion>,
    stored_forward: &mut [PafHit],
    stored_reverse: &[PafHit],
    reads: &ReadStore,
    query: &[u8],
    read_identifier: u32,
    extz: &mut BitExtz,
) {
    // `srt_i`: the stored overlaps keyed by `(target << 1 | strand)` in the
    // high 32 bits and by their index in the low 32 bits. The low bit records
    // which list the entry came from (`1` for the reverse list).
    let mut sorted_index: Vec<u64> =
        Vec::with_capacity(stored_forward.len() + stored_reverse.len());
    for (stored_index, hit) in stored_forward.iter().enumerate() {
        let mut key = (u64::from(hit.target_identifier) << 1) | u64::from(hit.reverse);
        key <<= 32;
        key |= (stored_index as u64) << 1;
        sorted_index.push(key);
    }
    for (stored_index, hit) in stored_reverse.iter().enumerate() {
        let mut key = (u64::from(hit.target_identifier) << 1) | u64::from(hit.reverse);
        key <<= 32;
        key |= ((stored_index as u64) << 1) | 1;
        sorted_index.push(key);
    }
    radix_sort_u64(&mut sorted_index);

    // The chains are sorted by target, so the cursor only moves forward: each
    // distinct `(target, strand)` matches the first stored overlap of that
    // key, as the reference's `for` loop does.
    let mut cursor = 0usize;
    let mut kept = 0usize;
    for index in 0..overlaps.len() {
        let target = overlaps[index].y_identifier;
        let strand = overlaps[index].y_strand;
        let target_key = (u64::from(target) << 1) | u64::from(strand);
        overlaps[index].non_homopolymer_errors = 0;

        while cursor < sorted_index.len() && (sorted_index[cursor] >> 32) < target_key {
            cursor += 1;
        }

        let query_start = overlaps[index].x_position_start;
        let query_end = overlaps[index].x_position_end.wrapping_add(1);
        let target_start = overlaps[index].y_position_start;
        let target_end = overlaps[index].y_position_end.wrapping_add(1);
        // The reference widens both ends to `uint64_t` before subtracting.
        let chain_query_length = u64::from(query_end).wrapping_sub(u64::from(query_start));
        let chain_target_length = u64::from(target_end).wrapping_sub(u64::from(target_start));

        let mut keep = false;
        if cursor < sorted_index.len() && (sorted_index[cursor] >> 32) == target_key {
            let entry = sorted_index[cursor];
            let from_reverse = entry & 1 == 1;
            let stored_index = ((entry as u32) >> 1) as usize;
            let stored = if from_reverse {
                stored_reverse[stored_index]
            } else {
                stored_forward[stored_index]
            };

            let mut is_match: u8 = if from_reverse { 2 } else { 1 };
            overlaps[index].shared_seed = 0;
            overlaps[index].is_match = 0;
            overlaps[index].strong = stored.matching_length as i8;
            overlaps[index].without_large_indel = u8::from(stored.without_large_indel);

            let stored_query_start = u64::from(stored.query_start);
            let stored_query_end = u64::from(stored.query_end);
            let stored_target_start = u64::from(stored.target_start);
            let stored_target_end = u64::from(stored.target_end);

            let mut matched = true;

            let overlap_start = stored_query_start.max(u64::from(query_start));
            let overlap_end = stored_query_end.min(u64::from(query_end));
            let overlap = if overlap_end > overlap_start {
                overlap_end - overlap_start
            } else {
                0
            };
            if !overlap_matches(
                overlap,
                stored_query_end.wrapping_sub(stored_query_start),
                chain_query_length,
            ) {
                matched = false;
            }
            overlaps[index].non_homopolymer_errors = overlaps[index]
                .non_homopolymer_errors
                .wrapping_add(stored_query_end.wrapping_sub(stored_query_start).wrapping_sub(overlap) as u32);

            let overlap_start = stored_target_start.max(u64::from(target_start));
            let overlap_end = stored_target_end.min(u64::from(target_end));
            let overlap = if overlap_end > overlap_start {
                overlap_end - overlap_start
            } else {
                0
            };
            if !overlap_matches(
                overlap,
                stored_target_end.wrapping_sub(stored_target_start),
                chain_target_length,
            ) {
                matched = false;
            }
            overlaps[index].non_homopolymer_errors = overlaps[index]
                .non_homopolymer_errors
                .wrapping_add(stored_target_end.wrapping_sub(stored_target_start).wrapping_sub(overlap) as u32);

            if matched {
                // A reproduced forward overlap no longer needs to be appended
                // as a standalone exact region.
                if is_match == 1 && stored.exact_region {
                    stored_forward[stored_index].exact_region = false;
                    stored_forward[stored_index].exact_region_byte = 0;
                }

                // `exact_ec_check` against the target substring the chain
                // covers, in the chain's strand.
                let target_bases = reads.read_sub_region(
                    target as usize,
                    target_start as i64,
                    chain_target_length as i64,
                    strand as u8,
                );
                if exact_ec_check(query, query_start, query_end, &target_bases) {
                    if is_match == 2 {
                        overlaps[index].strong = 0;
                        overlaps[index].without_large_indel = 1;
                    }
                    is_match = 1;
                    overlaps[index].shared_seed = 1;
                }
                overlaps[index].is_match = is_match;
                keep = true;
            }
        } else {
            // No stored overlap for this target: the chain survives only when
            // it is exactly equal to the target substring.
            let target_bases = reads.read_sub_region(
                target as usize,
                target_start as i64,
                chain_target_length as i64,
                strand as u8,
            );
            if exact_ec_check(query, query_start, query_end, &target_bases) {
                overlaps[index].strong = 0;
                overlaps[index].without_large_indel = 1;
                overlaps[index].shared_seed = 1;
                overlaps[index].is_match = 1;
                keep = true;
            }
        }

        if keep {
            if kept != index {
                overlaps.swap(kept, index);
            }
            kept += 1;
        }
    }
    overlaps.truncate(kept);

    // Stored exact overlaps that no chain reproduced are kept as overlaps in
    // their own right.
    let mut appended = false;
    for stored in stored_forward.iter() {
        if !stored.exact_region {
            continue;
        }

        let x_position_start = stored.query_start;
        let x_position_end = stored.query_end.wrapping_sub(1);
        let y_position_start = stored.target_start;
        let y_position_end = stored.target_end.wrapping_sub(1);
        let mut region = OverlapRegion {
            x_identifier: read_identifier,
            x_position_start,
            x_position_end,
            x_strand: 0,
            y_identifier: stored.target_identifier,
            y_position_start,
            y_position_end,
            y_strand: u32::from(stored.reverse),
            non_homopolymer_errors: 0,
            is_match: 1,
            shared_seed: 1,
            strong: stored.matching_length as i8,
            without_large_indel: u8::from(stored.without_large_indel),
            ..OverlapRegion::default()
        };
        region.overlap_length = x_position_end
            .wrapping_add(1)
            .wrapping_sub(x_position_start);
        region.alignment_length = region.overlap_length;

        // `set_exact_exz` + `push_alnw`: one exact window spanning the whole
        // stored region.
        set_exact_exz(
            extz,
            i64::from(x_position_start),
            i64::from(x_position_end) + 1,
            i64::from(y_position_start),
            i64::from(y_position_end) + 1,
        );
        push_alignment_window(&mut region.window_cigar, extz);
        overlaps.push(region);
        appended = true;
    }

    if appended {
        sort_by_y_identifier(overlaps);
    }

    // `for duplicated chains`: one overlap per target read, the one with the
    // fewest errors; ties go to the longer query span.
    if overlaps.len() > 1 {
        let length = overlaps.len();
        let mut survivors = 0usize;
        let mut group_start = 0usize;
        for index in 1..=length {
            if index == length || overlaps[index].y_identifier != overlaps[group_start].y_identifier
            {
                let chosen = if index - group_start > 1 {
                    choose_best_region(overlaps, group_start, index)
                } else {
                    Some(group_start)
                };
                if let Some(chosen) = chosen {
                    if chosen != survivors {
                        overlaps.swap(chosen, survivors);
                    }
                    survivors += 1;
                }
                group_start = index;
            }
        }
        overlaps.truncate(survivors);
    }
}

/// Whether the shared part of the two intervals is long enough: non-empty and
/// at least [`SHARED_SEED_THRESHOLD`] of both.
fn overlap_matches(overlap: u64, stored_length: u64, chain_length: u64) -> bool {
    overlap != 0
        && (overlap as f64) >= (stored_length as f64) * SHARED_SEED_THRESHOLD
        && (overlap as f64) >= (chain_length as f64) * SHARED_SEED_THRESHOLD
}

/// Picks the region with the fewest non-homopolymer errors in
/// `overlaps[group_start..group_end]`.
///
/// The reference compares the negated error counts as `int64_t`; on a tie the
/// region with the longer query span wins.
fn choose_best_region(overlaps: &[OverlapRegion], group_start: usize, group_end: usize) -> Option<usize> {
    let mut best_score = i64::from(i32::MIN);
    let mut best: Option<usize> = None;
    for index in group_start..group_end {
        let region = &overlaps[index];
        let score = -i64::from(region.non_homopolymer_errors);
        let span = region
            .x_position_end
            .wrapping_add(1)
            .wrapping_sub(region.x_position_start);
        let replace = score > best_score
            || (score == best_score
                && best.map_or(true, |best| {
                    overlaps[best]
                        .x_position_end
                        .wrapping_add(1)
                        .wrapping_sub(overlaps[best].x_position_start)
                        < span
                }));
        if replace {
            best_score = score;
            best = Some(index);
        }
    }
    best
}

/// `push_ff_ovlp`: copies the surviving regions into the two stored overlap
/// lists, mirroring the target coordinates of reverse hits and updating the
/// per-thread counters.
pub(crate) fn push_final_overlaps(
    overlaps: &[OverlapRegion],
    reads: &ReadStore,
) -> (Vec<PafHit>, Vec<PafHit>, FinalOverlapCounters) {
    let forward_count = overlaps
        .iter()
        .filter(|region| region.is_match == 1)
        .count();
    let reverse_count = overlaps
        .iter()
        .filter(|region| region.is_match == 2)
        .count();
    let mut forward = Vec::with_capacity(forward_count);
    let mut reverse = Vec::with_capacity(reverse_count);
    let mut counters = FinalOverlapCounters::default();

    for region in overlaps {
        match region.is_match {
            1 => {
                let hit = final_hit(region, reads);
                // `push_ff_ovlp` only counts the forward list.
                if hit.matching_length == 1 {
                    counters.strong += 1;
                }
                if hit.matching_length == 0 {
                    counters.weak += 1;
                }
                if hit.exact_region {
                    counters.exact += 1;
                }
                if hit.without_large_indel {
                    counters.without_large_indel += 1;
                }
                forward.push(hit);
            }
            2 => reverse.push(final_hit(region, reads)),
            _ => {}
        }
    }

    counters.forward = forward.len() as u64;
    counters.reverse = reverse.len() as u64;
    (forward, reverse, counters)
}

/// One region in the `ma_hit_t` layout.
fn final_hit(region: &OverlapRegion, reads: &ReadStore) -> PafHit {
    let block_length = reads.read_length(region.y_identifier as usize) as u32;
    let reverse = region.y_strand != 0;
    let (target_start, target_end) = if reverse {
        // `z->ts = bl - y_pos_e - 1; z->te = bl - y_pos_s;`
        (
            block_length
                .wrapping_sub(region.y_position_end)
                .wrapping_sub(1),
            block_length.wrapping_sub(region.y_position_start),
        )
    } else {
        (
            region.y_position_start,
            region.y_position_end.wrapping_add(1),
        )
    };

    PafHit {
        query_identifier: region.x_identifier,
        query_start: region.x_position_start,
        query_end: region.x_position_end.wrapping_add(1),
        target_identifier: region.y_identifier,
        target_start,
        target_end,
        // `ml` is a one-bit field in the reference.
        matching_length: (region.strong as u32) & 1,
        block_length,
        reverse,
        without_large_indel: region.without_large_indel != 0,
        // `el` is the region's `shared_seed` byte.
        exact_region: region.shared_seed != 0,
        exact_region_byte: region.shared_seed as u8,
        deleted: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn store() -> ReadStore {
        let mut reads = ReadStore::default();
        reads.push_read("query", &vec![b'A'; 100]);
        reads.push_read("target", &vec![b'C'; 80]);
        reads
    }

    fn region(
        y_identifier: u32,
        y_strand: u32,
        x_start: u32,
        x_end: u32,
        y_start: u32,
        y_end: u32,
    ) -> OverlapRegion {
        OverlapRegion {
            x_identifier: 0,
            x_position_start: x_start,
            x_position_end: x_end,
            y_identifier,
            y_position_start: y_start,
            y_position_end: y_end,
            y_strand,
            is_match: 1,
            shared_seed: 1,
            strong: 1,
            ..OverlapRegion::default()
        }
    }

    #[test]
    fn forward_overlaps_keep_their_coordinates() {
        let reads = store();
        let hits = vec![region(1, 0, 10, 29, 5, 24)];
        let (forward, reverse, counters) = push_final_overlaps(&hits, &reads);

        assert_eq!(reverse.len(), 0);
        assert_eq!(forward.len(), 1);
        assert_eq!(forward[0].query_start, 10);
        assert_eq!(forward[0].query_end, 30);
        assert_eq!(forward[0].target_start, 5);
        assert_eq!(forward[0].target_end, 25);
        assert_eq!(forward[0].block_length, 80);
        assert!(!forward[0].reverse);
        assert_eq!(counters.forward, 1);
        assert_eq!(counters.strong, 1);
        assert_eq!(counters.weak, 0);
        assert_eq!(counters.exact, 1);
        assert_eq!(counters.without_large_indel, 0);
    }

    #[test]
    fn reverse_overlaps_mirror_the_target_interval() {
        let reads = store();
        // `is_match == 2` is what `push_ff_ovlp` selects for `reverse_paf`.
        let mut hits = vec![region(1, 1, 0, 9, 0, 9)];
        hits[0].is_match = 2;
        let (forward, reverse, counters) = push_final_overlaps(&hits, &reads);

        assert_eq!(forward.len(), 0);
        assert_eq!(reverse.len(), 1);
        assert!(reverse[0].reverse);
        // bl = 80: ts = 80 - 9 - 1 = 70, te = 80 - 0 = 80.
        assert_eq!(reverse[0].target_start, 70);
        assert_eq!(reverse[0].target_end, 80);
        // Query coordinates are not mirrored.
        assert_eq!(reverse[0].query_start, 0);
        assert_eq!(reverse[0].query_end, 10);
        // Only the forward list feeds the counters.
        assert_eq!(counters.forward, 0);
        assert_eq!(counters.reverse, 1);
    }

    #[test]
    fn duplicated_chains_keep_the_region_with_fewest_errors() {
        let mut overlaps = vec![
            region(3, 0, 0, 99, 0, 99),
            region(3, 0, 10, 60, 10, 60),
            region(5, 0, 0, 50, 0, 50),
        ];
        // Same error count: the longer span wins, which is the first region.
        overlaps[0].non_homopolymer_errors = 4;
        overlaps[1].non_homopolymer_errors = 4;
        // The second target has fewer errors and wins its group.
        let chosen = choose_best_region(&overlaps, 0, 2);
        assert_eq!(chosen, Some(0));
        assert_eq!(choose_best_region(&overlaps, 2, 3), Some(2));

        overlaps[1].non_homopolymer_errors = 0;
        assert_eq!(choose_best_region(&overlaps, 0, 2), Some(1));
    }

    #[test]
    fn overlap_match_requires_the_shared_share_on_both_sides() {
        assert!(overlap_matches(90, 100, 100));
        assert!(!overlap_matches(0, 100, 100));
        // 86 of 100 is below 0.866666 * 100 = 86.6666.
        assert!(!overlap_matches(86, 100, 100));
        assert!(overlap_matches(87, 100, 100));
        // The overlap must cover the shorter side too.
        assert!(!overlap_matches(90, 100, 200));
    }
}

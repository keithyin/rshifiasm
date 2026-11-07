//! The chaining pass of the new error-correction pipeline, the counterpart of
//! `lchain_qgen_mcopy_fast` and `lchain_qdp_mcopy_fast` in `Hash_Table.cpp`.
//!
//! `h_ec_lchain` (`anchor.cpp`) sketches the query, turns the minimizers into
//! the qgen anchor list, and hands that list to this module. Unlike the older
//! `chain_DP`, the qgen pass can emit several sub-chains per candidate read
//! (the `mcopy` branch), filters chains by occurrence-type cutoffs, and keeps
//! the anchor lists of every chain inside the candidate list so that the
//! haplotype and consensus stages can look them up through
//! `overlap_region::non_homopolymer_errors`.

use crate::container::ksort::{introsort_klib_by, radix_sort_by_key};
use crate::overlap::chain_score::{chain_band_width, chain_link_score, INVALID_SCORE};
use crate::overlap::engine::overlap_type;
use crate::overlap::region::{chain_length, KmerHit, OverlapRegion};

/// `OFL`: fraction of a low-support chain's query span that must be covered by
/// a stronger chain before the low-support one is dropped.
const OVERLAP_FRACTION: f64 = 0.95;

/// `CH_OCC`
const CHAIN_OCCURRENCE_SHIFT: u32 = 4;

/// `CH_SC`
const CHAIN_SCORE_FACTOR: i64 = 16;

/// Scratch state of `quick_ck_lchain`, returned to the caller.
#[derive(Debug, Clone, Copy)]
struct QuickChainState {
    /// Minimum score seen on a sorted run (`plus`).
    plus: i64,
    /// Score of the best chain (`msc`).
    maximum_score: i64,
    /// Index of the anchor ending the best chain (`msc_i`).
    maximum_score_index: i64,
    /// Query span of the best chain (`movl`).
    minimum_overlap_length: i64,
    /// First anchor the dynamic programme has to consider (`si`).
    start_index: i64,
    /// One past the last anchor the dynamic programme has to consider (`ei`).
    end_index: i64,
}

impl QuickChainState {
    fn empty(anchor_count: i64) -> Self {
        QuickChainState {
            plus: 0,
            maximum_score: i64::from(i32::MIN),
            maximum_score_index: i64::from(i32::MIN),
            minimum_overlap_length: i64::MAX,
            start_index: 0,
            end_index: anchor_count,
        }
    }
}

/// `quick_ck_lchain`: detects runs of anchors that are already collinear and
/// can skip the quadratic dynamic programme.
///
/// `predecessor` (`p`), `scratch` (`t`), `score` (`f`) and `occurrence` (`ii`)
/// are the scratch arrays of the reference `Chain_Data`.
#[allow(clippy::too_many_arguments)]
fn quick_chain_check(
    anchors: &[KmerHit],
    predecessor: &mut [i64],
    scratch: &mut [i64],
    score: &mut [i32],
    occurrence: &mut [i32],
    query_length: i64,
    target_length: i64,
    gap_penalty: f64,
    skip_penalty: f64,
    band_width_rate: f64,
) -> QuickChainState {
    let anchor_count = anchors.len() as i64;
    let mut state = QuickChainState::empty(anchor_count);
    if anchor_count <= 0 {
        return state;
    }

    let mut run_start = 0usize;
    let mut sorted = true;

    for position in 1..=anchors.len() {
        let position_i = position as i64;
        let end_of_run = position == anchors.len()
            || anchors[position].reverse_strand != anchors[run_start].reverse_strand;

        if end_of_run {
            let last = position - 1;
            scratch[last] = 0;
            occurrence[last] = 0;

            if sorted {
                let mut plus_run: i64 = 0;
                let mut best_score: i64 = i64::from(i32::MIN);
                let mut best_index: i64 = i64::from(i32::MIN);
                let mut total_diagonal_gap: i64 = 0;

                let first = run_start;
                predecessor[first] = -1;
                score[first] = (anchors[first].count & 0xff) as i32;
                if i64::from(score[first]) >= best_score {
                    best_score = i64::from(score[first]);
                    best_index = first as i64;
                }
                if i64::from(score[first]) < plus_run {
                    plus_run = i64::from(score[first]);
                }

                let mut cursor = run_start + 1;
                while cursor < position {
                    let current = &anchors[cursor];
                    let previous = &anchors[cursor - 1];
                    let query_distance = current.self_offset as i64 - previous.self_offset as i64;
                    if query_distance <= 0 {
                        break;
                    }
                    let target_distance = current.offset as i64 - previous.offset as i64;
                    if target_distance <= 0 {
                        break;
                    }
                    let diagonal_gap = (target_distance - query_distance).abs();
                    if diagonal_gap > 16
                        && diagonal_gap
                            > i64::from(chain_band_width(
                                previous,
                                current,
                                band_width_rate,
                                query_length,
                                target_length,
                            ))
                    {
                        break;
                    }

                    let shared_length = target_distance.min(query_distance);
                    let span = (current.count & 0xff) as i64;
                    let occurrence_weight = (current.count >> 8) as i64;
                    let mut link_score = if span < shared_length {
                        span
                    } else {
                        shared_length
                    };
                    link_score = if occurrence_weight > 0 && link_score >= occurrence_weight {
                        link_score / occurrence_weight
                    } else {
                        1
                    };

                    if diagonal_gap != 0 || (shared_length > span && shared_length > 0) {
                        let mut linear_penalty = gap_penalty * diagonal_gap as f64;
                        let adaptive_penalty = link_score as f64
                            * ((diagonal_gap as f64 / shared_length as f64) / band_width_rate);
                        if diagonal_gap < 4 {
                            linear_penalty = linear_penalty.min(adaptive_penalty);
                        } else {
                            linear_penalty = linear_penalty.max(adaptive_penalty);
                        }
                        linear_penalty += skip_penalty * shared_length as f64;
                        link_score -= linear_penalty as i32 as i64;
                    }

                    link_score += i64::from(score[cursor - 1]);
                    let span_score = (current.count & 0xff) as i64;
                    if link_score < span_score {
                        break;
                    }

                    predecessor[cursor] = (cursor - 1) as i64;
                    score[cursor] = link_score as i32;
                    total_diagonal_gap += diagonal_gap;

                    if i64::from(score[cursor]) >= best_score {
                        best_score = i64::from(score[cursor]);
                        best_index = cursor as i64;
                    }
                    if i64::from(score[cursor]) < plus_run {
                        plus_run = i64::from(score[cursor]);
                    }
                    cursor += 1;
                }

                if cursor >= position && best_index == (position - 1) as i64 {
                    if position - run_start >= 2
                        && total_diagonal_gap > 16
                        && total_diagonal_gap
                            > i64::from(chain_band_width(
                                &anchors[run_start],
                                &anchors[position - 1],
                                band_width_rate,
                                query_length,
                                target_length,
                            ))
                    {
                        best_index = i64::from(i32::MIN);
                    }

                    if best_index == (position - 1) as i64 {
                        if best_score >= state.maximum_score {
                            let overlap_length = chain_length(
                                anchors[best_index as usize].self_offset as i64,
                                anchors[best_index as usize].self_offset as i64,
                                query_length,
                                anchors[best_index as usize].offset as i64,
                                anchors[best_index as usize].offset as i64,
                                target_length,
                            );
                            if best_score > state.maximum_score
                                || overlap_length < state.minimum_overlap_length
                            {
                                state.maximum_score = best_score;
                                state.maximum_score_index = best_index;
                                state.minimum_overlap_length = overlap_length;
                            }
                        }
                        if plus_run < state.plus {
                            state.plus = plus_run;
                        }
                        if state.end_index > position_i {
                            state.start_index = position_i;
                        } else {
                            state.end_index = run_start as i64;
                        }
                    }
                }
            }

            run_start = position;
            sorted = true;
        } else {
            if anchors[position].self_offset <= anchors[position - 1].self_offset
                || anchors[position].offset <= anchors[position - 1].offset
            {
                sorted = false;
            }
            scratch[position - 1] = 0;
            occurrence[position - 1] = 0;
        }
    }

    state
}

/// `push_ovlp_chain_qgen`: turns the two end anchors of a chain into an
/// `overlap_region`.
///
/// The query interval is extended so that both reads start and end at the same
/// relative position; the qgen variant keeps the query coordinates
/// unmirrored (the chain anchors are already expressed on the query strand).
fn push_chain_overlap(
    overlaps: &mut Vec<OverlapRegion>,
    query_identifier: u32,
    query_length: i64,
    target_length: i64,
    score: i64,
    beginning: &KmerHit,
    end: &KmerHit,
) {
    let mut region = OverlapRegion::default();
    region.x_identifier = query_identifier;
    region.y_identifier = beginning.read_identifier;
    region.x_strand = 0;
    region.y_strand = u32::from(beginning.reverse_strand);
    region.x_position_start = beginning.self_offset;
    region.y_position_start = beginning.offset;
    region.x_position_end = end.self_offset;
    region.y_position_end = end.offset;

    if region.x_position_start <= region.y_position_start {
        region.y_position_start -= region.x_position_start;
        region.x_position_start = 0;
    } else {
        region.x_position_start -= region.y_position_start;
        region.y_position_start = 0;
    }

    let query_right = query_length - region.x_position_end as i64 - 1;
    let target_right = target_length - region.y_position_end as i64 - 1;
    if query_right <= target_right {
        region.x_position_end = (query_length - 1) as u32;
        region.y_position_end += query_right as u32;
    } else {
        region.y_position_end = (target_length - 1) as u32;
        region.x_position_end += target_right as u32;
    }

    region.shared_seed = score as i32;
    region.alignment_length = 0;
    region.is_match = 0;
    region.non_homopolymer_errors = 0;
    region.strong = 0;
    region.overlap_length = 0;
    overlaps.push(region);
}

/// `gen_fake_cigar`: records where the diagonal offset of the chain changes.
///
/// `append_both_ends` is the `apend_be` flag of the reference. The three
/// coordinate values describe the overlap the chain belongs to.
pub fn generate_fake_cigar(
    cigar: &mut crate::overlap::region::FakeCigar,
    x_position_start: u32,
    y_position_start: u32,
    x_position_end: u32,
    append_both_ends: bool,
    hits: &[KmerHit],
) {
    cigar.clear();
    if append_both_ends {
        cigar.push(x_position_start, 0);
    }

    let mut previous_delta = i64::from(i32::MAX);
    for hit in hits {
        let query_distance = hit.self_offset as i64 - i64::from(x_position_start);
        let target_distance = hit.offset as i64 - i64::from(y_position_start);
        let diagonal_delta = target_distance - query_distance;
        if diagonal_delta != previous_delta {
            previous_delta = diagonal_delta;
            cigar.push(hit.self_offset, diagonal_delta as i32);
        }
    }

    if append_both_ends && !cigar.is_empty() {
        let last = cigar.len() - 1;
        if cigar.gap_position(last) != i64::from(x_position_end) {
            let shift = cigar.gap_shift(last);
            cigar.push(x_position_end, shift as i32);
        }
    }
}

/// `lchain_qdp_mcopy_fast`: dynamic programme over one candidate read.
///
/// Returns the number of anchors copied into the candidate list for the
/// emitted chains. `candidates` mirrors `Candidates_list::list`; anchors of
/// every emitted chain are appended at `destination_index`.
#[allow(clippy::too_many_arguments)]
fn chain_candidate_read(
    candidates: &mut Vec<KmerHit>,
    anchor_index: usize,
    anchor_count: usize,
    destination_index: usize,
    overlaps: &mut Vec<OverlapRegion>,
    maximum_skip: i64,
    maximum_iterations: i64,
    maximum_distance: i64,
    gap_penalty: f64,
    skip_penalty: f64,
    band_width_rate: f64,
    query_identifier: u32,
    query_length: i64,
    target_length: i64,
    quick_check: bool,
    append_both_ends: bool,
    generate_cigar: bool,
    copy_number: i64,
    copy_rate: f64,
    copy_kmer_cutoff: i64,
    keep_anchor_count: bool,
) -> usize {
    if anchor_count == 0 {
        return 0;
    }

    let anchors: Vec<KmerHit> =
        candidates[anchor_index..anchor_index + anchor_count].to_vec();
    let mut predecessor = vec![-1i64; anchor_count];
    let mut scratch = vec![0i64; anchor_count];
    let mut score = vec![0i32; anchor_count];
    let mut occurrence = vec![0i32; anchor_count];

    let mut state = if quick_check {
        quick_chain_check(
            &anchors,
            &mut predecessor,
            &mut scratch,
            &mut score,
            &mut occurrence,
            query_length,
            target_length,
            gap_penalty,
            skip_penalty,
            band_width_rate,
        )
    } else {
        QuickChainState::empty(anchor_count as i64)
    };

    let mut start = state.start_index;
    let mut maximum_anchor_score_index: i64 = -1;

    let mut index = start;
    while index < state.end_index {
        let index_usize = index as usize;
        let mut maximum_score = i64::from(anchors[index_usize].count & 0xff);
        let mut skip_count: i64 = 0;
        let mut maximum_predecessor: i64 = -1;
        let end_predecessor: i64;

        if index - start > maximum_iterations {
            start = index - maximum_iterations;
        }
        while anchors[index_usize].reverse_strand
            != anchors[start as usize].reverse_strand
        {
            start += 1;
        }

        let mut predecessor_index = index - 1;
        while predecessor_index >= start {
            let link = chain_link_score(
                &anchors[predecessor_index as usize],
                &anchors[index_usize],
                band_width_rate,
                gap_penalty,
                skip_penalty,
                query_length,
                target_length,
            );
            if link == INVALID_SCORE {
                predecessor_index -= 1;
                continue;
            }
            let link_score = i64::from(link) + i64::from(score[predecessor_index as usize]);
            if link_score > maximum_score {
                maximum_score = link_score;
                maximum_predecessor = predecessor_index;
                if skip_count > 0 {
                    skip_count -= 1;
                }
            } else if scratch[predecessor_index as usize] == index {
                skip_count += 1;
                if skip_count > maximum_skip {
                    break;
                }
            }
            if predecessor[predecessor_index as usize] >= 0 {
                let target = predecessor[predecessor_index as usize] as usize;
                scratch[target] = index;
            }
            predecessor_index -= 1;
        }
        end_predecessor = predecessor_index;

        if maximum_anchor_score_index < 0
            || anchors[index_usize].self_offset
                > anchors[maximum_anchor_score_index as usize].self_offset as u32
                    + maximum_distance as u32
            || anchors[index_usize].reverse_strand
                != anchors[maximum_anchor_score_index as usize].reverse_strand
        {
            let mut best: i64 = i64::from(i32::MIN);
            maximum_anchor_score_index = -1;
            let mut candidate = index - 1;
            while candidate >= start
                && (anchors[index_usize].self_offset as i64)
                    <= maximum_distance + anchors[candidate as usize].self_offset as i64
                && anchors[index_usize].reverse_strand
                    == anchors[candidate as usize].reverse_strand
            {
                if best < i64::from(score[candidate as usize]) {
                    best = i64::from(score[candidate as usize]);
                    maximum_anchor_score_index = candidate;
                }
                candidate -= 1;
            }
        }

        if maximum_anchor_score_index >= 0
            && maximum_anchor_score_index < end_predecessor
            && anchors[index_usize].reverse_strand
                == anchors[maximum_anchor_score_index as usize].reverse_strand
        {
            let link = chain_link_score(
                &anchors[maximum_anchor_score_index as usize],
                &anchors[index_usize],
                band_width_rate,
                gap_penalty,
                skip_penalty,
                query_length,
                target_length,
            );
            if link != INVALID_SCORE {
                let candidate_score = i64::from(link)
                    + i64::from(score[maximum_anchor_score_index as usize]);
                if maximum_score < candidate_score {
                    maximum_score = candidate_score;
                    maximum_predecessor = maximum_anchor_score_index;
                }
            }
        }

        score[index_usize] = maximum_score as i32;
        predecessor[index_usize] = maximum_predecessor;

        if maximum_anchor_score_index < 0
            || ((anchors[index_usize].self_offset as i64)
                <= maximum_distance + anchors[maximum_anchor_score_index as usize].self_offset as i64
                && anchors[index_usize].reverse_strand
                    == anchors[maximum_anchor_score_index as usize].reverse_strand
                && score[maximum_anchor_score_index as usize] < score[index_usize])
        {
            maximum_anchor_score_index = index;
        }

        if i64::from(score[index_usize]) >= state.maximum_score {
            let overlap_length = chain_length(
                anchors[index_usize].self_offset as i64,
                anchors[index_usize].self_offset as i64,
                query_length,
                anchors[index_usize].offset as i64,
                anchors[index_usize].offset as i64,
                target_length,
            );
            if i64::from(score[index_usize]) > state.maximum_score
                || overlap_length < state.minimum_overlap_length
            {
                state.maximum_score = i64::from(score[index_usize]);
                state.maximum_score_index = index;
                state.minimum_overlap_length = overlap_length;
            }
        }
        if i64::from(score[index_usize]) < state.plus {
            state.plus = i64::from(score[index_usize]);
        }
        occurrence[index_usize] = 0;

        index += 1;
    }

    // Trace the best chain backwards, labelling its anchors in `occurrence`.
    let mut chain_length_count: usize = 0;
    let mut trace = state.maximum_score_index;
    while trace >= 0 {
        occurrence[trace as usize] = 1;
        if chain_length_count < scratch.len() {
            scratch[chain_length_count] = trace;
        }
        chain_length_count += 1;
        trace = predecessor[trace as usize];
    }

    let mut copied_anchors = anchor_count;

    if copy_number > 1 {
        if chain_length_count as i64 >= copy_kmer_cutoff {
            let best_score = state.maximum_score - state.plus;
            let minimum_score = (best_score as f64 * copy_rate) as i64;
            occurrence[state.maximum_score_index as usize] = 0;

            // Make every chain score positive and collect the candidate ends.
            let mut candidate_chain_count: usize = 0;
            for candidate in 0..anchor_count {
                score[candidate] = (i64::from(score[candidate]) - state.plus) as i32;
                if candidate >= candidate_chain_count {
                    scratch[candidate] = 0;
                }
                if occurrence[candidate] == 0
                    && i64::from(score[candidate]) >= minimum_score
                {
                    scratch[candidate_chain_count] =
                        ((u64::from(score[candidate] as u32) << 32)
                            + ((candidate as u64) << 1)) as i64;
                    candidate_chain_count += 1;
                }
            }

            if candidate_chain_count > 1 {
                radix_sort_by_key(&mut scratch[..candidate_chain_count], &|value: &i64| {
                    *value as u64
                });

                let initial_overlap_count = overlaps.len();
                let mut emitted: usize = 0;
                let mut visited_anchors: usize = 0;
                let mut candidate = candidate_chain_count as i64 - 1;

                while candidate >= 0 && emitted < copy_number as usize {
                    let candidate_usize = candidate as usize;
                    let anchor_visit_start = visited_anchors;
                    let mut anchor = ((scratch[candidate_usize] as u32) >> 1) as i64;

                    while anchor >= 0 && (scratch[anchor as usize] & 1) == 0 {
                        occurrence[visited_anchors] = anchor as i32;
                        visited_anchors += 1;
                        scratch[anchor as usize] |= 1;
                        anchor = predecessor[anchor as usize];
                    }

                    if anchor_visit_start == visited_anchors {
                        candidate -= 1;
                        continue;
                    }

                    let chain_score = if anchor < 0 {
                        scratch[candidate_usize] >> 32
                    } else {
                        (scratch[candidate_usize] >> 32) - i64::from(score[anchor as usize])
                    };

                    if chain_score >= minimum_score {
                        let beginning = anchors[occurrence[visited_anchors - 1] as usize];
                        let end = anchors[occurrence[anchor_visit_start] as usize];
                        let overlap_start = overlaps.len();
                        push_chain_overlap(
                            overlaps,
                            query_identifier,
                            query_length,
                            target_length,
                            chain_score + state.plus,
                            &beginning,
                            &end,
                        );
                        let new_anchor_count = visited_anchors - anchor_visit_start;
                        if emitted == 0 || new_anchor_count > 1 {
                            overlaps[overlap_start].alignment_length = new_anchor_count as u32;
                            overlaps[overlap_start].x_identifier = anchor_visit_start as u32;
                            emitted += 1;
                        } else {
                            overlaps.pop();
                            visited_anchors = anchor_visit_start;
                        }
                    } else {
                        visited_anchors = anchor_visit_start;
                    }
                    candidate -= 1;
                }

                let new_overlap_end = overlaps.len();
                if new_overlap_end > initial_overlap_count + 1 {
                    let total_new_anchors = visited_anchors;
                    if candidates.len() < destination_index + total_new_anchors {
                        candidates.resize(destination_index + total_new_anchors, KmerHit::default());
                    }
                    let mut copied = 0usize;
                    for overlap_index in initial_overlap_count..new_overlap_end {
                        let region = &mut overlaps[overlap_index];
                        region.non_homopolymer_errors = (destination_index + copied) as u32;
                        let first_anchor = region.x_identifier as usize;
                        let anchor_run = region.alignment_length as usize;
                        for offset in 0..anchor_run {
                            let source = occurrence[first_anchor + (anchor_run - offset - 1)];
                            let mut hit = anchors[source as usize];
                            hit.read_identifier = overlap_index as u32;
                            candidates[destination_index + copied] = hit;
                            copied += 1;
                        }
                        region.x_identifier = query_identifier;
                        if generate_cigar {
                            let start = destination_index + copied - anchor_run;
                            let x_start = region.x_position_start;
                            let y_start = region.y_position_start;
                            let x_end = region.x_position_end;
                            generate_fake_cigar(
                                &mut region.fake_cigar,
                                x_start,
                                y_start,
                                x_end,
                                append_both_ends,
                                &candidates[start..destination_index + copied],
                            );
                        }
                        if !keep_anchor_count {
                            region.alignment_length = 0;
                        }
                    }
                    copied_anchors = copied;
                } else if new_overlap_end == initial_overlap_count + 1 {
                    let region = &mut overlaps[initial_overlap_count];
                    let anchor_run = region.alignment_length as usize;
                    let first_anchor = region.x_identifier as usize;
                    if candidates.len() < destination_index + anchor_run {
                        candidates.resize(destination_index + anchor_run, KmerHit::default());
                    }
                    for offset in 0..anchor_run {
                        let source = occurrence[first_anchor + (anchor_run - offset - 1)];
                        let mut hit = anchors[source as usize];
                        hit.read_identifier = initial_overlap_count as u32;
                        candidates[destination_index + offset] = hit;
                    }
                    region.non_homopolymer_errors = destination_index as u32;
                    region.x_identifier = query_identifier;
                    if generate_cigar {
                        let x_start = region.x_position_start;
                        let y_start = region.y_position_start;
                        let x_end = region.x_position_end;
                        generate_fake_cigar(
                            &mut region.fake_cigar,
                            x_start,
                            y_start,
                            x_end,
                            append_both_ends,
                            &candidates[destination_index..destination_index + anchor_run],
                        );
                    }
                    if !keep_anchor_count {
                        region.alignment_length = 0;
                    }
                    copied_anchors = anchor_run;
                }
                return copied_anchors;
            } else {
                // Too few distinct sub-chains: fall back to the single best one.
                state.maximum_score += state.plus;
                let mut chain_end = state.maximum_score_index;
                let mut chain_index = 0usize;
                while chain_end >= 0 {
                    scratch[chain_index] = chain_end;
                    chain_index += 1;
                    chain_end = predecessor[chain_end as usize];
                }
                chain_length_count = chain_index;
            }
        }
    }

    // Emit the single best chain.
    let overlap_index = overlaps.len();
    let chain_start_anchor = anchors[scratch[chain_length_count - 1] as usize];
    let chain_end_anchor = anchors[scratch[0] as usize];
    push_chain_overlap(
        overlaps,
        query_identifier,
        query_length,
        target_length,
        state.maximum_score,
        &chain_start_anchor,
        &chain_end_anchor,
    );
    if candidates.len() < destination_index + chain_length_count {
        candidates.resize(destination_index + chain_length_count, KmerHit::default());
    }
    for offset in 0..chain_length_count {
        let source = scratch[chain_length_count - offset - 1] as usize;
        let mut hit = anchors[source];
        hit.read_identifier = overlap_index as u32;
        candidates[destination_index + offset] = hit;
    }
    let region = &mut overlaps[overlap_index];
    region.non_homopolymer_errors = destination_index as u32;
    if generate_cigar {
        let x_start = region.x_position_start;
        let y_start = region.y_position_start;
        let x_end = region.x_position_end;
        generate_fake_cigar(
            &mut region.fake_cigar,
            x_start,
            y_start,
            x_end,
            append_both_ends,
            &candidates[destination_index..destination_index + chain_length_count],
        );
    }
    if keep_anchor_count {
        region.alignment_length = chain_length_count as u32;
    }
    chain_length_count
}

/// Parameters of the qgen chain pass, mirroring the call site in
/// `worker_hap_ec`.
#[derive(Debug, Clone, Copy)]
pub struct QgenChainParameters {
    pub maximum_chain_count: u32,
    pub maximum_skip: i64,
    pub maximum_iterations: i64,
    pub maximum_distance: i64,
    pub gap_penalty: f64,
    pub skip_penalty: f64,
    pub band_width_rate: f64,
    pub quick_check: bool,
    pub append_both_ends: bool,
    pub generate_offset: bool,
    pub copy_number: i64,
    pub copy_rate: f64,
    pub chain_cutoff: u32,
    pub copy_kmer_cutoff: u32,
    pub coverage_window: u64,
}

/// Result of [`lchain_qgen_mcopy_fast`].
#[derive(Debug, Clone, Default)]
pub struct QgenChainResult {
    /// Chains emitted for this query, ordered by query start.
    pub overlaps: Vec<OverlapRegion>,
    /// The candidate list after the pass; chain anchors live in its first
    /// `candidate_length` entries.
    pub candidates: Vec<KmerHit>,
    /// Number of valid entries in `candidates` (`Candidates_list::length`).
    pub candidate_length: usize,
}

/// `lchain_qgen_mcopy_fast`: groups the sorted qgen anchors by candidate read,
/// chains every group and applies the chain-count / coverage filters.
pub fn lchain_qgen_mcopy_fast(
    candidates: Vec<KmerHit>,
    query_identifier: u32,
    query_length: u64,
    read_lengths: &[u64],
    parameters: &QgenChainParameters,
) -> QgenChainResult {
    let mut candidate_list = candidates;
    let mut overlaps: Vec<OverlapRegion> = Vec::new();
    let mut chain_cutoff_hit = false;
    let mut copied_anchors: usize = 0;

    let candidate_count = candidate_list.len();
    let mut group_start = 0usize;
    let mut position = 1usize;

    while position <= candidate_count {
        if position == candidate_count
            || candidate_list[position].read_identifier
                != candidate_list[group_start].read_identifier
        {
            let target_identifier = candidate_list[group_start].read_identifier;
            if target_identifier != query_identifier {
                let overlap_start = overlaps.len();
                let target_length = read_lengths
                    .get(target_identifier as usize)
                    .copied()
                    .unwrap_or(0);
                let copied = chain_candidate_read(
                    &mut candidate_list,
                    group_start,
                    position - group_start,
                    copied_anchors,
                    &mut overlaps,
                    parameters.maximum_skip,
                    parameters.maximum_iterations,
                    parameters.maximum_distance,
                    parameters.gap_penalty,
                    parameters.skip_penalty,
                    parameters.band_width_rate,
                    query_identifier,
                    query_length as i64,
                    target_length as i64,
                    parameters.quick_check,
                    parameters.append_both_ends,
                    parameters.generate_offset,
                    parameters.copy_number,
                    parameters.copy_rate,
                    i64::from(parameters.copy_kmer_cutoff),
                    true,
                );
                copied_anchors += copied;
                if parameters.chain_cutoff >= 2 && !chain_cutoff_hit {
                    for region in &overlaps[overlap_start..] {
                        if region.alignment_length < parameters.chain_cutoff {
                            chain_cutoff_hit = true;
                            break;
                        }
                    }
                }
            }
            group_start = position;
        }
        position += 1;
    }
    candidate_list.truncate(copied_anchors);

    // Limit the number of chains per overlap-type class.
    let maximum_chain_count = parameters.maximum_chain_count;
    if overlaps.len() as u64 > u64::from(maximum_chain_count) {
        introsort_klib_by(&mut overlaps, &|left: &OverlapRegion, right: &OverlapRegion| {
            // `oreg_ss_lt`: descending shared-seed order.
            left.shared_seed > right.shared_seed
        });

        let mut type_counts = [0u32; 4];
        let mut type_thresholds = [0i32; 4];
        for region in &overlaps {
            let overlap_class = overlap_type(region, query_length as u32);
            type_counts[overlap_class] += 1;
            if u64::from(type_counts[overlap_class]) == u64::from(maximum_chain_count) {
                type_thresholds[overlap_class] = region.shared_seed;
            }
        }

        if type_thresholds.iter().any(|threshold| *threshold > 0) {
            let mut coverage_counts: Vec<u64> = Vec::new();
            // The reference restarts the low-align flag in the `for` initialiser
            // (`for (i = 0, k = 0, lch = 0; ...)`), so whatever the per-read
            // pass saw is discarded once the threshold filter runs.
            chain_cutoff_hit = false;
            let mut kept: usize = 0;

            let coverage_window = parameters.coverage_window;
            if type_counts[3] as u64 >= u64::from(maximum_chain_count)
                && query_length >= coverage_window
            {
                let count = query_length / coverage_window
                    + u64::from(query_length % coverage_window != 0);
                coverage_counts.resize(count as usize, 0);
                let mut window_start = 0u64;
                for entry in coverage_counts.iter_mut() {
                    let window_end = (window_start + coverage_window).min(query_length);
                    let mut coverage =
                        (window_end - window_start) * u64::from(maximum_chain_count >> 1);
                    if coverage > u32::MAX as u64 {
                        coverage = u32::MAX as u64;
                    }
                    *entry = coverage << 32;
                    window_start += coverage_window;
                }
            }

            for index in 0..overlaps.len() {
                let region = overlaps[index].clone();
                let overlap_class = overlap_type(&region, query_length as u32);
                let threshold = type_thresholds[overlap_class];

                if region.shared_seed >= threshold {
                    if !coverage_counts.is_empty() {
                        let mut bucket =
                            (region.x_position_start as u64 / coverage_window) as usize;
                        let query_start = region.x_position_start as u64;
                        let query_end = region.x_position_end as u64 + 1;
                        let mut window_start = bucket as u64 * coverage_window;
                        while bucket < coverage_counts.len() {
                            let window_end = (window_start + coverage_window).min(query_length);
                            let overlap_start = query_start.max(window_start);
                            let overlap_end = query_end.min(window_end);
                            if overlap_end <= overlap_start {
                                break;
                            }
                            let covered = overlap_end - overlap_start;
                            let low = coverage_counts[bucket] as u32 as u64;
                            if low + covered < u32::MAX as u64 {
                                coverage_counts[bucket] += covered;
                            } else {
                                coverage_counts[bucket] =
                                    ((coverage_counts[bucket] >> 32) << 32) | u32::MAX as u64;
                            }
                            window_start += coverage_window;
                            bucket += 1;
                        }
                    }

                    if kept != index {
                        overlaps.swap(kept, index);
                    }
                    if overlaps[kept].alignment_length < parameters.chain_cutoff {
                        chain_cutoff_hit = true;
                    }
                    kept += 1;
                } else if overlap_class == 3 && !coverage_counts.is_empty() {
                    let mut bucket =
                        (region.x_position_start as u64 / coverage_window) as usize;
                    let query_start = region.x_position_start as u64;
                    let query_end = region.x_position_end as u64 + 1;
                    let mut window_start = bucket as u64 * coverage_window;
                    let mut uncovered: u64 = 0;
                    let mut covered: u64 = 0;
                    while bucket < coverage_counts.len() {
                        let window_end = (window_start + coverage_window).min(query_length);
                        let overlap_start = query_start.max(window_start);
                        let overlap_end = query_end.min(window_end);
                        if overlap_end <= overlap_start {
                            break;
                        }
                        let length = overlap_end - overlap_start;
                        let low = coverage_counts[bucket] as u32 as u64;
                        let capacity = coverage_counts[bucket] >> 32;
                        if length + low >= capacity {
                            covered += length;
                        } else {
                            uncovered += length;
                        }
                        window_start += coverage_window;
                        bucket += 1;
                    }

                    if (uncovered as f64) >= ((uncovered + covered) as f64 * 0.7) {
                        let mut bucket =
                            (region.x_position_start as u64 / coverage_window) as usize;
                        let mut window_start = bucket as u64 * coverage_window;
                        while bucket < coverage_counts.len() {
                            let window_end = (window_start + coverage_window).min(query_length);
                            let overlap_start = query_start.max(window_start);
                            let overlap_end = query_end.min(window_end);
                            if overlap_end <= overlap_start {
                                break;
                            }
                            let covered_length = overlap_end - overlap_start;
                            let low = coverage_counts[bucket] as u32 as u64;
                            if low + covered_length < u32::MAX as u64 {
                                coverage_counts[bucket] += covered_length;
                            } else {
                                coverage_counts[bucket] =
                                    ((coverage_counts[bucket] >> 32) << 32) | u32::MAX as u64;
                            }
                            window_start += coverage_window;
                            bucket += 1;
                        }
                        if kept != index {
                            overlaps.swap(kept, index);
                        }
                        if overlaps[kept].alignment_length < parameters.chain_cutoff {
                            chain_cutoff_hit = true;
                        }
                        kept += 1;
                    }
                }
            }
            overlaps.truncate(kept);
        }
    }

    // The window/pair logic of `get_interval` depends on this order.
    introsort_klib_by(&mut overlaps, &|left: &OverlapRegion, right: &OverlapRegion| {
        // `oreg_xs_lt`.
        ((u64::from(left.x_position_start) << 32) | u64::from(left.x_position_end))
            < ((u64::from(right.x_position_start) << 32) | u64::from(right.x_position_end))
    });

    if chain_cutoff_hit {
        // `@brief r485`: low-support chains are dropped when a stronger chain
        // covers most of their query span.
        let mut kept = 0usize;
        for index in 0..overlaps.len() {
            let region = overlaps[index].clone();
            if region.alignment_length < parameters.chain_cutoff {
                let query_start = region.x_position_start as u64;
                let query_end = region.x_position_end as u64 + 1;
                // `ob` is a uint64 in the reference: the fractional bound is
                // truncated before the integer comparison.
                let mut overlap_bound =
                    ((query_end - query_start) as f64 * OVERLAP_FRACTION) as u64;
                if overlap_bound < 16 {
                    overlap_bound = 16;
                }
                let minimum_seed_score = i64::from(region.shared_seed) * CHAIN_SCORE_FACTOR;
                let minimum_occurrence =
                    u64::from(region.alignment_length) << CHAIN_OCCURRENCE_SHIFT;

                let mut stronger = 0usize;
                while stronger < overlaps.len()
                    && query_end > overlaps[stronger].x_position_start as u64
                {
                    let other = overlaps[stronger].clone();
                    if other.alignment_length < parameters.chain_cutoff {
                        stronger += 1;
                        continue;
                    }
                    if u64::from(other.alignment_length) < minimum_occurrence {
                        stronger += 1;
                        continue;
                    }
                    if i64::from(other.shared_seed) < minimum_seed_score {
                        stronger += 1;
                        continue;
                    }
                    let other_start = other.x_position_start as u64;
                    let other_end = other.x_position_end as u64 + 1;
                    let overlap_start = other_start.max(query_start);
                    let overlap_end = other_end.min(query_end);
                    if overlap_end > overlap_start
                        && overlap_end - overlap_start >= overlap_bound
                    {
                        let anchor_index = other.non_homopolymer_errors as usize;
                        let target_identifier = candidate_list[anchor_index].read_identifier;
                        let mut contained = 0u64;
                        let mut cursor = anchor_index;
                        while cursor < candidate_list.len()
                            && candidate_list[cursor].read_identifier == target_identifier
                            && contained < minimum_occurrence
                        {
                            let end = candidate_list[cursor].self_offset;
                            let start = end
                                .wrapping_sub(candidate_list[cursor].count & 0xff);
                            if u64::from(start) >= overlap_start
                                && u64::from(end) <= overlap_end
                            {
                                contained += 1;
                            }
                            cursor += 1;
                        }
                        if contained >= minimum_occurrence {
                            break;
                        }
                    }
                    stronger += 1;
                }

                if stronger < overlaps.len()
                    && query_end > overlaps[stronger].x_position_start as u64
                {
                    continue;
                }
            }

            if kept != index {
                overlaps.swap(kept, index);
            }
            kept += 1;
        }
        overlaps.truncate(kept);
    }

    for region in &mut overlaps {
        region.alignment_length = 0;
    }

    let candidate_length = candidate_list.len();
    QgenChainResult {
        overlaps,
        candidates: candidate_list,
        candidate_length,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(read: u32, reverse: bool, self_offset: u32, offset: u32, weight: u32, span: u8) -> KmerHit {
        KmerHit {
            read_identifier: read,
            reverse_strand: reverse,
            offset,
            self_offset,
            count: (weight << 8) | span as u32,
        }
    }

    #[test]
    fn a_single_collinear_read_becomes_one_chain() {
        let anchors = vec![
            hit(1, false, 100, 100, 2, 51),
            hit(1, false, 150, 150, 2, 51),
            hit(1, false, 200, 200, 2, 51),
        ];
        let parameters = QgenChainParameters {
            maximum_chain_count: 500,
            maximum_skip: 25,
            maximum_iterations: 5000,
            maximum_distance: 5000,
            gap_penalty: 0.30024778842926025,
            skip_penalty: 0.00030024780269026496,
            band_width_rate: 0.02,
            quick_check: true,
            append_both_ends: true,
            generate_offset: true,
            copy_number: 3,
            copy_rate: 0.7,
            chain_cutoff: 2,
            copy_kmer_cutoff: 32,
            coverage_window: 3072,
        };
        let result = lchain_qgen_mcopy_fast(anchors, 0, 1000, &[1000, 1000], &parameters);
        assert_eq!(result.overlaps.len(), 1);
        let region = &result.overlaps[0];
        assert_eq!(region.y_identifier, 1);
        assert_eq!((region.x_position_start, region.x_position_end), (0, 999));
        assert_eq!((region.y_position_start, region.y_position_end), (0, 999));
        assert!(!region.fake_cigar.is_empty());
    }

    #[test]
    fn a_self_group_is_skipped() {
        let anchors = vec![
            hit(0, false, 100, 100, 2, 51),
            hit(0, false, 150, 150, 2, 51),
        ];
        let parameters = QgenChainParameters {
            maximum_chain_count: 500,
            maximum_skip: 25,
            maximum_iterations: 5000,
            maximum_distance: 5000,
            gap_penalty: 0.3,
            skip_penalty: 0.0003,
            band_width_rate: 0.02,
            quick_check: true,
            append_both_ends: true,
            generate_offset: true,
            copy_number: 3,
            copy_rate: 0.7,
            chain_cutoff: 2,
            copy_kmer_cutoff: 32,
            coverage_window: 3072,
        };
        let result = lchain_qgen_mcopy_fast(anchors, 0, 1000, &[1000], &parameters);
        assert!(result.overlaps.is_empty());
    }

    // The values below were printed by the unmodified `Hash_Table.o` through
    // `/tmp/qgen_probe.cpp`; see docs/consistency.md for the probe recipe.
    const GAP: f64 = 0.30024778842926025;
    const SKIP: f64 = 0.00030024780269026496;
    const RATE: f64 = 0.02;

    fn collinear_anchors() -> Vec<KmerHit> {
        (0..4)
            .map(|step| hit(1, false, 100 + step * 50, 100 + step * 50, 2, 51))
            .collect()
    }

    fn long_two_chain_anchors() -> Vec<KmerHit> {
        let mut anchors: Vec<KmerHit> = (0..40)
            .map(|step| hit(1, false, 100 + step * 10, 100 + step * 10, 2, 51))
            .collect();
        anchors.extend((0..40).map(|step| hit(1, false, 500 + step * 10, 540 + step * 10, 2, 51)));
        anchors
    }

    fn tiny_two_chain_anchors() -> Vec<KmerHit> {
        let mut anchors: Vec<KmerHit> = (0..4)
            .map(|step| hit(1, false, 100 + step * 50, 100 + step * 50, 2, 51))
            .collect();
        anchors.extend((0..4).map(|step| hit(1, false, 400 + step * 50, 440 + step * 50, 2, 51)));
        anchors
    }

    fn run_quick(anchors: &[KmerHit]) -> (QuickChainState, Vec<i32>, Vec<i64>) {
        let count = anchors.len();
        let mut predecessor = vec![0i64; count];
        let mut scratch = vec![0i64; count];
        let mut score = vec![0i32; count];
        let mut occurrence = vec![0i32; count];
        let state = quick_chain_check(
            anchors,
            &mut predecessor,
            &mut scratch,
            &mut score,
            &mut occurrence,
            1000,
            1000,
            GAP,
            SKIP,
            RATE,
        );
        (state, score, predecessor)
    }

    #[allow(clippy::too_many_arguments)]
    fn run_dp(
        anchors: &[KmerHit],
        quick_check: bool,
        copy_number: i64,
        copy_kmer_cutoff: i64,
        keep_anchor_count: bool,
        append_both_ends: bool,
    ) -> (usize, Vec<OverlapRegion>, Vec<KmerHit>) {
        let mut candidates = anchors.to_vec();
        candidates.resize(anchors.len() + 256, KmerHit::default());
        let mut overlaps: Vec<OverlapRegion> = Vec::new();
        let copied = chain_candidate_read(
            &mut candidates,
            0,
            anchors.len(),
            0,
            &mut overlaps,
            25,
            5000,
            5000,
            GAP,
            SKIP,
            RATE,
            0,
            1000,
            1000,
            quick_check,
            append_both_ends,
            true,
            copy_number,
            0.7,
            copy_kmer_cutoff,
            keep_anchor_count,
        );
        (copied, overlaps, candidates)
    }

    #[test]
    fn quick_check_matches_the_reference_on_a_sorted_run() {
        let (state, score, predecessor) = run_quick(&collinear_anchors());
        assert_eq!(state.plus, 0);
        assert_eq!(state.maximum_score, 126);
        assert_eq!(state.maximum_score_index, 3);
        assert_eq!(state.minimum_overlap_length, 1000);
        assert_eq!(state.start_index, 0);
        assert_eq!(state.end_index, 0);
        assert_eq!(score, vec![51, 76, 101, 126]);
        assert_eq!(predecessor, vec![-1, 0, 1, 2]);
    }

    #[test]
    fn quick_check_matches_the_reference_on_an_unsorted_run() {
        let anchors = vec![
            hit(1, false, 100, 100, 2, 51),
            hit(1, false, 150, 160, 2, 51),
            hit(1, false, 200, 155, 2, 51),
            hit(1, false, 250, 250, 2, 51),
        ];
        let (state, score, _) = run_quick(&anchors);
        assert_eq!(state.maximum_score, i64::from(i32::MIN));
        assert_eq!(state.maximum_score_index, i64::from(i32::MIN));
        assert_eq!(state.minimum_overlap_length, i64::MAX);
        assert_eq!(state.start_index, 0);
        assert_eq!(state.end_index, 4);
        assert_eq!(score, vec![0, 0, 0, 0]);
    }

    #[test]
    fn quick_check_only_reports_a_run_that_reaches_the_run_end() {
        let anchors = vec![
            hit(1, false, 100, 100, 2, 51),
            hit(1, false, 150, 150, 2, 51),
            hit(1, false, 200, 200, 2, 51),
            hit(1, false, 250, 250, 2, 51),
            hit(1, false, 500, 540, 2, 51),
        ];
        let (state, score, _) = run_quick(&anchors);
        assert_eq!(state.maximum_score, i64::from(i32::MIN));
        assert_eq!(state.maximum_score_index, i64::from(i32::MIN));
        assert_eq!(state.start_index, 0);
        assert_eq!(state.end_index, 5);
        assert_eq!(score[..4], [51, 76, 101, 126]);
    }

    #[test]
    fn dp_matches_the_reference_for_one_collinear_chain() {
        for quick_check in [true, false] {
            let (copied, overlaps, candidates) =
                run_dp(&collinear_anchors(), quick_check, 3, 32, true, true);
            assert_eq!(copied, 4);
            assert_eq!(overlaps.len(), 1);
            let region = &overlaps[0];
            assert_eq!(region.x_identifier, 0);
            assert_eq!((region.x_position_start, region.x_position_end), (0, 999));
            assert_eq!(region.y_identifier, 1);
            assert_eq!((region.y_position_start, region.y_position_end), (0, 999));
            assert_eq!(region.shared_seed, 126);
            assert_eq!(region.alignment_length, 4);
            assert_eq!(region.non_homopolymer_errors, 0);
            assert_eq!(
                region.fake_cigar.entries,
                vec![(0u64 << 32) | 0, (100u64 << 32) | 0, (999u64 << 32) | 0]
            );
            for offset in 0..4u32 {
                let candidate = candidates[offset as usize];
                assert_eq!(candidate.read_identifier, 0);
                assert_eq!(candidate.self_offset, 100 + offset * 50);
                assert_eq!(candidate.offset, 100 + offset * 50);
                assert_eq!(candidate.count, 563);
            }
        }
    }

    #[test]
    fn dp_matches_the_reference_for_two_long_chains() {
        let (copied, overlaps, candidates) =
            run_dp(&long_two_chain_anchors(), true, 3, 32, true, true);
        assert_eq!(copied, 80);
        assert_eq!(overlaps.len(), 2);

        let first = &overlaps[0];
        assert_eq!((first.x_position_start, first.x_position_end), (0, 959));
        assert_eq!((first.y_position_start, first.y_position_end), (40, 999));
        assert_eq!(first.shared_seed, 246);
        assert_eq!(first.alignment_length, 40);
        assert_eq!(first.non_homopolymer_errors, 0);
        assert_eq!(
            first.fake_cigar.entries,
            vec![0, (500u64 << 32), (959u64 << 32)]
        );
        for offset in 0..40u32 {
            let candidate = candidates[offset as usize];
            assert_eq!(candidate.read_identifier, 0);
            assert_eq!(candidate.self_offset, 500 + offset * 10);
            assert_eq!(candidate.offset, 540 + offset * 10);
        }

        let second = &overlaps[1];
        assert_eq!((second.x_position_start, second.x_position_end), (0, 999));
        assert_eq!((second.y_position_start, second.y_position_end), (0, 999));
        assert_eq!(second.shared_seed, 246);
        assert_eq!(second.alignment_length, 40);
        assert_eq!(second.non_homopolymer_errors, 40);
        assert_eq!(
            second.fake_cigar.entries,
            vec![0, (100u64 << 32), (999u64 << 32)]
        );
        for offset in 0..40u32 {
            let candidate = candidates[40 + offset as usize];
            assert_eq!(candidate.read_identifier, 1);
            assert_eq!(candidate.self_offset, 100 + offset * 10);
            assert_eq!(candidate.offset, 100 + offset * 10);
        }
    }

    #[test]
    fn dp_matches_the_reference_when_anchor_counts_are_not_kept() {
        let (copied, overlaps, _) = run_dp(&long_two_chain_anchors(), true, 3, 32, false, true);
        assert_eq!(copied, 80);
        assert_eq!(overlaps.len(), 2);
        assert!(overlaps.iter().all(|region| region.alignment_length == 0));
        assert_eq!(overlaps[0].non_homopolymer_errors, 0);
        assert_eq!(overlaps[1].non_homopolymer_errors, 40);
    }

    #[test]
    fn dp_matches_the_reference_for_tiny_mcopy_chains() {
        let (copied, overlaps, _) = run_dp(&tiny_two_chain_anchors(), true, 3, 1, true, true);
        assert_eq!(copied, 8);
        assert_eq!(overlaps.len(), 2);
        assert_eq!(
            (overlaps[0].x_position_start, overlaps[0].x_position_end),
            (0, 959)
        );
        assert_eq!(
            (overlaps[0].y_position_start, overlaps[0].y_position_end),
            (40, 999)
        );
        assert_eq!(overlaps[0].alignment_length, 4);
        assert_eq!(overlaps[0].non_homopolymer_errors, 0);
        assert_eq!(
            (overlaps[1].x_position_start, overlaps[1].x_position_end),
            (0, 999)
        );
        assert_eq!(
            (overlaps[1].y_position_start, overlaps[1].y_position_end),
            (0, 999)
        );
        assert_eq!(overlaps[1].alignment_length, 4);
        assert_eq!(overlaps[1].non_homopolymer_errors, 4);
    }

    #[test]
    fn dp_matches_the_reference_without_the_both_ends_cigar_entry() {
        let (_, overlaps, _) = run_dp(&long_two_chain_anchors(), true, 3, 32, true, false);
        assert_eq!(overlaps.len(), 2);
        assert_eq!(overlaps[0].fake_cigar.entries, vec![500u64 << 32]);
        assert_eq!(overlaps[1].fake_cigar.entries, vec![100u64 << 32]);
    }

    // Wrapper-level goldens.  The values come from the same probe, extended
    // with the `lchain_qgen_mcopy_fast` block that links the real `anchor.o`.
    fn wrapper_parameters(
        maximum_chain_count: u32,
        chain_cutoff: u32,
        coverage_window: u64,
    ) -> QgenChainParameters {
        QgenChainParameters {
            maximum_chain_count,
            maximum_skip: 25,
            maximum_iterations: 5000,
            maximum_distance: 5000,
            gap_penalty: GAP,
            skip_penalty: SKIP,
            band_width_rate: RATE,
            quick_check: true,
            append_both_ends: true,
            generate_offset: true,
            copy_number: 1,
            copy_rate: 0.7,
            chain_cutoff,
            copy_kmer_cutoff: 32,
            coverage_window,
        }
    }

    #[test]
    fn wrapper_matches_the_reference_on_grouping_and_self_skip() {
        let mut anchors: Vec<KmerHit> = (0..3)
            .map(|step| hit(0, false, 100 + step * 50, 100 + step * 50, 2, 51))
            .collect();
        anchors.extend(
            (0..4).map(|step| hit(1, false, 100 + step * 50, 100 + step * 50, 2, 51)),
        );
        anchors.extend(
            (0..3).map(|step| hit(2, false, 500 + step * 50, 500 + step * 50, 2, 51)),
        );
        let result = lchain_qgen_mcopy_fast(
            anchors,
            0,
            1000,
            &[1000, 1000, 800],
            &wrapper_parameters(100, 2, 3072),
        );

        assert_eq!(result.candidate_length, 7);
        assert_eq!(result.candidates.len(), 7);
        assert_eq!(result.overlaps.len(), 2);

        let first = &result.overlaps[0];
        assert_eq!(first.y_identifier, 2);
        assert_eq!((first.x_position_start, first.x_position_end), (0, 799));
        assert_eq!((first.y_position_start, first.y_position_end), (0, 799));
        assert_eq!(first.shared_seed, 101);
        assert_eq!(first.non_homopolymer_errors, 4);
        assert_eq!(
            first.fake_cigar.entries,
            vec![0, 500u64 << 32, 799u64 << 32]
        );

        let second = &result.overlaps[1];
        assert_eq!(second.y_identifier, 1);
        assert_eq!((second.x_position_start, second.x_position_end), (0, 999));
        assert_eq!((second.y_position_start, second.y_position_end), (0, 999));
        assert_eq!(second.shared_seed, 126);
        assert_eq!(second.non_homopolymer_errors, 0);
        assert_eq!(
            second.fake_cigar.entries,
            vec![0, 100u64 << 32, 999u64 << 32]
        );

        // The query's own group is dropped, so the copied anchors are the two
        // target chains, tagged with the overlap index they produced.
        for step in 0..4u32 {
            let candidate = result.candidates[step as usize];
            assert_eq!(candidate.read_identifier, 0);
            assert_eq!(candidate.self_offset, 100 + step * 50);
        }
        for step in 0..3u32 {
            let candidate = result.candidates[4 + step as usize];
            assert_eq!(candidate.read_identifier, 1);
            assert_eq!(candidate.self_offset, 500 + step * 50);
        }
        assert!(result
            .overlaps
            .iter()
            .all(|region| region.alignment_length == 0));
    }

    #[test]
    fn wrapper_applies_the_chain_limit_and_the_type_three_coverage_rule() {
        let mut anchors: Vec<KmerHit> = (0..5)  // target 1, sc 151, x [50, 449]
            .map(|step| hit(1, false, 100 + step * 50, 50 + step * 50, 2, 51))
            .collect();
        anchors.extend(
            (0..4).map(|step| hit(2, false, 400 + step * 50, 350 + step * 50, 2, 51)),
        );
        anchors.extend((0..3).map(|step| hit(3, false, 750 + step * 10, step * 10, 2, 51)));
        anchors.extend(
            (0..2).map(|step| hit(4, false, 200 + step * 50, 150 + step * 50, 2, 51)),
        );
        let result = lchain_qgen_mcopy_fast(
            anchors,
            0,
            1000,
            &[0, 400, 600, 90, 300],
            &wrapper_parameters(2, 2, 100),
        );

        assert_eq!(result.candidate_length, 14);
        assert_eq!(result.overlaps.len(), 3);

        let first = &result.overlaps[0];
        assert_eq!(first.y_identifier, 1);
        assert_eq!((first.x_position_start, first.x_position_end), (50, 449));
        assert_eq!((first.y_position_start, first.y_position_end), (0, 399));
        assert_eq!(first.shared_seed, 151);
        assert_eq!(first.non_homopolymer_errors, 0);
        assert_eq!(
            first.fake_cigar.entries,
            vec![50u64 << 32, 100u64 << 32, 449u64 << 32]
        );

        let second = &result.overlaps[1];
        assert_eq!(second.y_identifier, 2);
        assert_eq!((second.x_position_start, second.x_position_end), (50, 649));
        assert_eq!((second.y_position_start, second.y_position_end), (0, 599));
        assert_eq!(second.shared_seed, 126);
        assert_eq!(second.non_homopolymer_errors, 5);

        // Target 4 is fully covered by the two chains above and is dropped;
        // target 3 sits in an uncovered stretch shorter than a window, so the
        // 0.7 rule keeps it.
        let third = &result.overlaps[2];
        assert_eq!(third.y_identifier, 3);
        assert_eq!((third.x_position_start, third.x_position_end), (750, 839));
        assert_eq!((third.y_position_start, third.y_position_end), (0, 89));
        assert_eq!(third.shared_seed, 61);
        assert_eq!(third.non_homopolymer_errors, 9);
        assert_eq!(
            third.fake_cigar.entries,
            vec![750u64 << 32, 750u64 << 32, 839u64 << 32]
        );

        // The dropped chain's anchors stay in the candidate list, tagged with
        // the overlap index they had when they were copied.
        assert_eq!(result.candidates[12].read_identifier, 3);
        assert_eq!(result.candidates[12].self_offset, 200);
        assert_eq!(result.candidates[13].read_identifier, 3);
        assert_eq!(result.candidates[13].self_offset, 250);
    }

    #[test]
    fn wrapper_drops_a_low_support_chain_covered_by_a_stronger_one() {
        let mut anchors: Vec<KmerHit> = (0..32)
            .map(|step| hit(1, false, 100 + step * 50, 100 + step * 50, 2, 51))
            .collect();
        anchors.push(hit(2, false, 500, 100, 2, 51));
        let result = lchain_qgen_mcopy_fast(
            anchors,
            0,
            2000,
            &[0, 2000, 2000],
            &wrapper_parameters(100, 2, 3072),
        );

        assert_eq!(result.candidate_length, 33);
        assert_eq!(result.overlaps.len(), 1);
        let region = &result.overlaps[0];
        assert_eq!(region.y_identifier, 1);
        assert_eq!((region.x_position_start, region.x_position_end), (0, 1999));
        assert_eq!((region.y_position_start, region.y_position_end), (0, 1999));
        assert_eq!(region.shared_seed, 826);
        assert_eq!(region.non_homopolymer_errors, 0);

        assert_eq!(result.candidates[32].read_identifier, 1);
        assert_eq!(result.candidates[32].self_offset, 500);
        assert_eq!(result.candidates[32].offset, 100);
    }

    #[test]
    fn wrapper_does_not_filter_when_no_type_reaches_the_limit() {
        let mut anchors: Vec<KmerHit> = (0..3)
            .map(|step| hit(1, false, 100 + step * 50, 100 + step * 50, 2, 51))
            .collect();
        anchors.extend((0..3).map(|step| hit(2, false, 300 + step * 50, 100 + step * 50, 2, 51)));
        anchors.extend((0..3).map(|step| hit(3, false, 500 + step * 50, 500 + step * 50, 2, 51)));
        let result = lchain_qgen_mcopy_fast(
            anchors,
            0,
            1000,
            &[0, 800, 800, 1000],
            &wrapper_parameters(2, 2, 3072),
        );

        assert_eq!(result.candidate_length, 9);
        assert_eq!(result.overlaps.len(), 3);

        let expected = [
            (1u32, 0u32, 799u32, 0u32),
            (3, 0, 999, 6),
            (2, 200, 999, 3),
        ];
        for (region, (target, start, end, anchor_index)) in
            result.overlaps.iter().zip(expected)
        {
            assert_eq!(region.y_identifier, target);
            assert_eq!((region.x_position_start, region.x_position_end), (start, end));
            assert_eq!(region.shared_seed, 101);
            assert_eq!(region.non_homopolymer_errors, anchor_index);
        }
    }

    #[test]
    fn wrapper_keeps_types_whose_threshold_stays_zero() {
        let mut anchors: Vec<KmerHit> = (0..5)
            .map(|step| hit(1, false, 100 + step * 50, 100 + step * 50, 2, 51))
            .collect();
        anchors.extend((0..4).map(|step| hit(2, false, 350 + step * 50, 350 + step * 50, 2, 51)));
        anchors.extend((0..3).map(|step| hit(3, false, 550 + step * 50, 550 + step * 50, 2, 51)));
        anchors.extend((0..2).map(|step| hit(4, false, 700 + step * 50, 500 + step * 50, 2, 51)));
        let result = lchain_qgen_mcopy_fast(
            anchors,
            0,
            1000,
            &[0, 800, 700, 800, 700],
            &wrapper_parameters(2, 2, 3072),
        );

        assert_eq!(result.candidate_length, 14);
        assert_eq!(result.overlaps.len(), 3);

        let expected = [
            (2u32, 0u32, 699u32, 126i32, 5u32),
            (1, 0, 799, 151, 0),
            (4, 200, 899, 76, 12),
        ];
        for (region, (target, start, end, seed, anchor_index)) in
            result.overlaps.iter().zip(expected)
        {
            assert_eq!(region.y_identifier, target);
            assert_eq!((region.x_position_start, region.x_position_end), (start, end));
            assert_eq!(region.shared_seed, seed);
            assert_eq!(region.non_homopolymer_errors, anchor_index);
        }
    }
}

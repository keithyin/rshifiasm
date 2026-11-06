//! Minimizer chaining, the counterpart of `ha_chain_check` and `chain_DP` in
//! `Hash_Table.cpp`.
//!
//! Chaining turns the list of shared seeds between two reads into one
//! candidate overlap: a dynamic programme scores every anchor against its
//! predecessors, keeping the anchors inside the allowed diagonal band, and the
//! best scoring anchor is traced back into the overlap interval plus the gap
//! structure used later by the alignment stage.

use crate::overlap::chain_score::INVALID_SCORE;
use crate::overlap::region::{chain_length, normal_weight, ChainData, KmerHit, OverlapRegion};

/// Largest diagonal jump that is still considered a normal gap
/// (`THRESHOLD_MAX_SIZE`).
const MAXIMUM_DIAGONAL_STEP: i64 = 31;

/// `ha_chain_check`: accepts the anchor list only when it is a strictly
/// increasing collinear chain within the band.
///
/// Returns the number of usable anchors, or `-1` when the list must go through
/// the full dynamic programme instead.
pub fn chain_check(
    hits: &[KmerHit],
    hit_count: i32,
    data: &mut ChainData,
    minimum_score: i64,
    band_width_threshold: f64,
) -> i32 {
    if hit_count == 0 {
        return -1;
    }

    let mut index = 1i32;
    while index < hit_count {
        if hits[(index - 1) as usize].self_offset >= hits[index as usize].self_offset {
            break;
        }
        index += 1;
    }
    if index < hit_count {
        return -1;
    }

    let band_width_penalty = 1.0 / band_width_threshold;
    let mut total_indels: i32 = 0;
    let mut total_length: i32 = 0;

    data.score[0] = normal_weight(minimum_score, hits[0].count as i64) as i32;
    data.predecessor[0] = -1;
    data.indels[0] = 0;
    data.self_length[0] = 0;
    data.occurrence[0] = 1;

    let mut index = 1i32;
    while index < hit_count {
        let delta_x = hits[index as usize].offset as i32
            - hits[(index - 1) as usize].offset as i32;
        let delta_y = hits[index as usize].self_offset as i32
            - hits[(index - 1) as usize].self_offset as i32;
        let diagonal_delta = if delta_x > delta_y {
            delta_x - delta_y
        } else {
            delta_y - delta_x
        };

        total_indels += diagonal_delta;
        total_length += delta_y;
        if total_indels as f64 > total_length as f64 * band_width_threshold {
            break;
        }

        let distance_min = if delta_x < delta_y { delta_x } else { delta_y };
        if diagonal_delta > MAXIMUM_DIAGONAL_STEP as i32
            && diagonal_delta as f64 > distance_min as f64 * band_width_threshold
        {
            break;
        }

        let mut score = if distance_min < minimum_score as i32 {
            distance_min
        } else {
            minimum_score as i32
        };
        score = normal_weight(score as i64, hits[index as usize].count as i64) as i32;

        let gap_rate = total_indels as f64 / total_length as f64;
        score -= (gap_rate * score as f64 * band_width_penalty) as i32;

        data.score[index as usize] = data.score[(index - 1) as usize] + score;
        data.predecessor[index as usize] = (index - 1) as i64;
        data.indels[index as usize] = total_indels;
        data.self_length[index as usize] = total_length;
        data.occurrence[index as usize] = index + 1;
        index += 1;
    }

    if index < hit_count {
        return -1;
    }
    hit_count
}

/// `chain_DP`: scores the anchors and traces back the best chain.
///
/// `hits` must be sorted by `offset`. Returns the number of anchors in the
/// chosen chain, and fills `result` with the overlap interval and gap
/// structure.
#[allow(clippy::too_many_arguments)]
pub fn chain_dp(
    hits: &[KmerHit],
    hit_count: i64,
    data: &mut ChainData,
    result: &mut OverlapRegion,
    band_width_threshold: f64,
    maximum_skip: i32,
    x_read_length: i64,
    y_read_length: i64,
    k_mer_length: i64,
) -> i64 {
    if hit_count <= 0 {
        // The reference would index `a[-1]` here; no caller passes an empty
        // list, so the guard keeps the contract explicit.
        return 0;
    }

    data.resize(hit_count as usize);
    let minimum_score = k_mer_length;
    let band_width_penalty = 1.0 / band_width_threshold;
    let mut anchor_count = hit_count;

    let collinear = chain_check(
        hits,
        hit_count as i32,
        data,
        minimum_score,
        band_width_threshold,
    );
    if collinear > 0 {
        anchor_count = collinear as i64;
    } else {
        // Fill the score and backtrack arrays.
        for index in 0..anchor_count {
            data.scratch[index as usize] = -1;
        }

        for index in 0..anchor_count {
            let position = hits[index as usize].offset as i64;
            let self_position = hits[index as usize].self_offset as i64;
            let mut best_predecessor: i64 = -1;
            let mut best_score = normal_weight(minimum_score, hits[index as usize].count as i64);
            let mut best_indels: i64 = 0;
            let mut best_self_length: i64 = 0;
            let mut maximum_skip_count = 0i32;
            let mut chain_skip_count = 0i32;

            let mut predecessor = index - 1;
            while predecessor >= 0 {
                let distance_position = position - hits[predecessor as usize].offset as i64;
                let distance_self_position =
                    self_position - hits[predecessor as usize].self_offset as i64;
                if distance_position == 0 || distance_self_position <= 0 {
                    predecessor -= 1;
                    continue;
                }

                let distance_gap = if distance_position > distance_self_position {
                    distance_position - distance_self_position
                } else {
                    distance_self_position - distance_position
                };

                let total_indels = data.indels[predecessor as usize] as i64 + distance_gap;
                let total_self_length =
                    data.self_length[predecessor as usize] as i64 + distance_self_position;
                if total_indels as f64 > band_width_threshold * total_self_length as f64 {
                    predecessor -= 1;
                    continue;
                }

                let distance_min = if distance_position < distance_self_position {
                    distance_position
                } else {
                    distance_self_position
                };
                let mut score = if distance_min < minimum_score {
                    distance_min
                } else {
                    minimum_score
                };
                score = normal_weight(score, hits[predecessor as usize].count as i64);

                let gap_rate = total_indels as f64 / total_self_length as f64;
                score -= (gap_rate * score as f64 * band_width_penalty) as i64;
                score += data.score[predecessor as usize] as i64;

                // A strict comparison keeps the earliest predecessor on ties.
                if score > best_score {
                    best_score = score;
                    best_predecessor = predecessor;
                    best_indels = total_indels;
                    best_self_length = total_self_length;
                    maximum_skip_count = 0;
                    if chain_skip_count > 0 {
                        chain_skip_count -= 1;
                    }
                } else {
                    maximum_skip_count += 1;
                    if maximum_skip_count > maximum_skip {
                        break;
                    }
                    if data.scratch[predecessor as usize] == index {
                        chain_skip_count += 1;
                        if chain_skip_count > maximum_skip {
                            break;
                        }
                    }
                }

                let predecessor_of_predecessor = data.predecessor[predecessor as usize];
                if predecessor_of_predecessor >= 0 {
                    data.scratch[predecessor_of_predecessor as usize] = index;
                }
                predecessor -= 1;
            }

            data.score[index as usize] = best_score as i32;
            data.predecessor[index as usize] = best_predecessor;
            data.indels[index as usize] = best_indels as i32;
            data.self_length[index as usize] = best_self_length as i32;
            data.occurrence[index as usize] = 1;
            if best_predecessor != -1 {
                data.occurrence[index as usize] =
                    data.occurrence[best_predecessor as usize] + 1;
            }
        }
    }

    // Choose the best scoring anchor, preferring the one with the shortest
    // implied overlap when scores tie.
    let mut best_score: i64 = -1;
    let mut best_index: i64 = -1;
    let mut shortest_length = x_read_length * 2 + 2;
    for index in 0..anchor_count {
        let score = data.score[index as usize] as i64;
        let candidate_length = chain_length(
            hits[index as usize].self_offset as i64,
            hits[index as usize].self_offset as i64,
            x_read_length,
            hits[index as usize].offset as i64,
            hits[index as usize].offset as i64,
            y_read_length,
        );
        if score > best_score {
            best_score = score;
            best_index = index;
            shortest_length = candidate_length;
        } else if score == best_score && candidate_length < shortest_length {
            best_index = index;
            shortest_length = candidate_length;
        }
    }

    result.fake_cigar.clear();
    let mut index = best_index;
    result.x_position_end = hits[index as usize].self_offset;
    result.y_position_end = hits[index as usize].offset;
    result.shared_seed = best_score as i32;
    result.overlap_length = shortest_length as u32;

    let mut distance_self_position =
        result.x_position_end as i64 - hits[index as usize].self_offset as i64;
    let mut distance_position = result.y_position_end as i64 - hits[index as usize].offset as i64;
    let mut previous_gap = distance_position - distance_self_position;
    result
        .fake_cigar
        .push(hits[index as usize].self_offset, previous_gap as i32);

    let mut chain_size = 0i64;
    if result.x_strand == 1 {
        while index >= 0 {
            distance_self_position =
                result.x_position_end as i64 - hits[index as usize].self_offset as i64;
            distance_position = result.y_position_end as i64 - hits[index as usize].offset as i64;
            let distance_gap = distance_position - distance_self_position;
            if distance_gap != previous_gap {
                previous_gap = distance_gap;
                result
                    .fake_cigar
                    .push(hits[index as usize].self_offset, previous_gap as i32);
            }

            chain_size += 1;
            result.x_position_start = hits[index as usize].self_offset;
            result.y_position_start = hits[index as usize].offset;
            index = data.predecessor[index as usize];
        }
    } else {
        while index >= 0 {
            distance_self_position =
                result.x_position_end as i64 - hits[index as usize].self_offset as i64;
            distance_position = result.y_position_end as i64 - hits[index as usize].offset as i64;
            let distance_gap = distance_position - distance_self_position;
            if distance_gap == previous_gap {
                // Drop the entry that was just added for the same gap.
                result.fake_cigar.entries.pop();
                result
                    .fake_cigar
                    .push(hits[index as usize].self_offset, previous_gap as i32);
            } else {
                previous_gap = distance_gap;
                result
                    .fake_cigar
                    .push(hits[index as usize].self_offset, previous_gap as i32);
            }

            chain_size += 1;
            result.x_position_start = hits[index as usize].self_offset;
            result.y_position_start = hits[index as usize].offset;
            index = data.predecessor[index as usize];
        }
    }

    chain_size
}

/// `calculate_overlap_region_by_chaining`: chains every group of seeds that
/// points at the same read and strand, appending one overlap region per group.
///
/// `read_lengths` is indexed by read identifier and provides the length of the
/// candidate reads, which chaining needs to convert an anchor into an overlap
/// length. `add_begin_and_end` is the reference `keep_whole_chain` flag.
pub fn chain_candidate_groups(
    seeds: &[KmerHit],
    read_identifier: u64,
    read_length: u64,
    read_lengths: &[u64],
    band_width_threshold: f64,
    k_mer_length: i64,
    add_begin_and_end: bool,
    data: &mut ChainData,
    overlap_list: &mut Vec<OverlapRegion>,
) {
    /// The reference passes a maximum skip of 25 in this call site.
    const MAXIMUM_SKIP: i32 = 25;

    if seeds.is_empty() {
        return;
    }

    let mut index = 0usize;
    while index < seeds.len() {
        let current_identifier = seeds[index].read_identifier;
        let current_strand = seeds[index].strand_value();
        let group_start = index;
        index += 1;
        while index < seeds.len()
            && seeds[index].read_identifier == current_identifier
            && seeds[index].strand_value() == current_strand
        {
            index += 1;
        }

        if current_identifier as u64 == read_identifier {
            continue;
        }

        // The working region mirrors `f_cigar` in the reference: the query is
        // the read being processed and the target the candidate.
        let mut region = OverlapRegion {
            x_identifier: read_identifier as u32,
            x_strand: current_strand,
            y_identifier: current_identifier,
            y_strand: 0,
            ..OverlapRegion::default()
        };
        chain_dp(
            &seeds[group_start..index],
            (index - group_start) as i64,
            data,
            &mut region,
            band_width_threshold,
            MAXIMUM_SKIP,
            read_length as i64,
            read_lengths[current_identifier as usize] as i64,
            k_mer_length,
        );

        crate::overlap::append::append_inexact_overlap(
            overlap_list,
            &region,
            read_length as i64,
            read_lengths[current_identifier as usize] as i64,
            add_begin_and_end,
        );
    }
}

/// `lchain_qcheck`: linear-time acceptance test used when `quick_check` is on.
///
/// The anchor list is accepted only when it is already a strictly increasing
/// collinear chain inside the diagonal band; on success the score and
/// predecessor arrays are filled linearly so the caller can backtrack without
/// running the quadratic programme.
pub fn chain_check_quick(
    hits: &[KmerHit],
    hit_count: i32,
    data: &mut ChainData,
    band_width_threshold: f64,
) -> Option<i32> {
    if hit_count == 0 {
        return None;
    }

    if hit_count > 1 {
        let first = &hits[0];
        let last = &hits[(hit_count - 1) as usize];
        if first.self_offset >= last.self_offset || first.offset >= last.offset {
            return None;
        }
        let delta_query = (last.self_offset as i32).wrapping_sub(first.self_offset as i32);
        let delta_target = (last.offset as i32).wrapping_sub(first.offset as i32);
        let diagonal = if delta_query >= delta_target {
            delta_query.wrapping_sub(delta_target)
        } else {
            delta_target.wrapping_sub(delta_query)
        };
        let distance = if delta_query >= delta_target {
            delta_target
        } else {
            delta_query
        };
        if distance == 0 || f64::from(diagonal) > f64::from(distance) * band_width_threshold {
            return None;
        }
    }

    let mut index = 1i32;
    while index < hit_count {
        let previous = &hits[(index - 1) as usize];
        let current = &hits[index as usize];
        if previous.self_offset >= current.self_offset {
            break;
        }
        if previous.offset >= current.offset {
            break;
        }
        index += 1;
    }
    if index < hit_count {
        return None;
    }

    let band_width_penalty = 1.0 / band_width_threshold;
    data.score[0] = normal_weight(
        (hits[0].count & 0xff) as i64,
        (hits[0].count >> 8) as i64,
    ) as i32;
    data.predecessor[0] = -1;

    let mut total_gap: i32 = 0;
    let mut index = 1i32;
    while index < hit_count {
        let previous = &hits[(index - 1) as usize];
        let current = &hits[index as usize];
        let delta_query = (current.self_offset as i32).wrapping_sub(previous.self_offset as i32);
        let delta_target = (current.offset as i32).wrapping_sub(previous.offset as i32);
        let diagonal = if delta_query >= delta_target {
            delta_query.wrapping_sub(delta_target)
        } else {
            delta_target.wrapping_sub(delta_query)
        };
        let distance = if delta_query >= delta_target {
            delta_target
        } else {
            delta_query
        };
        if distance == 0 {
            break;
        }

        total_gap = total_gap.wrapping_add(diagonal);
        if i64::from(diagonal) > MAXIMUM_DIAGONAL_STEP
            && f64::from(diagonal) > f64::from(distance) * band_width_threshold
        {
            break;
        }

        let span = (current.count & 0xff) as i32;
        let mut score = if distance < span { distance } else { span };
        score = normal_weight(score as i64, (current.count >> 8) as i64) as i32;
        score = score.wrapping_sub(
            ((f64::from(diagonal) / f64::from(distance)) * band_width_penalty * f64::from(score))
                as i32,
        );

        data.score[index as usize] = data.score[(index - 1) as usize].wrapping_add(score);
        data.predecessor[index as usize] = (index - 1) as i64;
        index += 1;
    }
    if index < hit_count {
        return None;
    }

    if hit_count > 1 {
        let first = &hits[0];
        let last = &hits[(hit_count - 1) as usize];
        let delta_query = (last.self_offset as i32).wrapping_sub(first.self_offset as i32);
        let delta_target = (last.offset as i32).wrapping_sub(first.offset as i32);
        let distance = if delta_query >= delta_target {
            delta_target
        } else {
            delta_query
        };
        if f64::from(total_gap) > f64::from(distance) * band_width_threshold {
            return None;
        }
    }

    Some(hit_count)
}

/// `comput_sc_ff`: the link score used by [`chain_dp_fixed_ends`].
///
/// Unlike `comput_sc_ch` the band-width rejection is commented out in the
/// reference, so a link is rejected only when one coordinate does not advance.
fn fixed_end_link_score(
    later: &KmerHit,
    earlier: &KmerHit,
    band_width_rate: f64,
    gap_penalty: f64,
    skip_penalty: f64,
) -> i32 {
    // The reference casts each position to `int32_t` before subtracting.
    let delta_query = (i64::from(later.self_offset) - i64::from(earlier.self_offset)) as i32;
    if delta_query < 0 {
        return INVALID_SCORE;
    }
    let delta_target = (i64::from(later.offset) - i64::from(earlier.offset)) as i32;
    if delta_target < 0 {
        return INVALID_SCORE;
    }

    let diagonal = if delta_target > delta_query {
        delta_target - delta_query
    } else {
        delta_query - delta_target
    };
    let distance = if delta_target < delta_query {
        delta_target
    } else {
        delta_query
    };
    let span = (later.count & 0xff) as i32;
    let mut score = if span < distance { span } else { distance };
    score = normal_weight(score as i64, (later.count >> 8) as i64) as i32;

    if diagonal != 0 || (distance > span && distance > 0) {
        let linear_penalty = gap_penalty * f64::from(diagonal);
        let affine_penalty =
            f64::from(score) * ((f64::from(diagonal) / f64::from(distance)) / band_width_rate);
        let mut penalty = if linear_penalty > affine_penalty {
            affine_penalty
        } else {
            linear_penalty
        };
        penalty += skip_penalty * f64::from(distance);
        score = score.wrapping_sub(penalty as i32);
    }

    score
}

/// `rev_khit`: mirrors one anchor around the read ends.
fn reverse_kmer_hit(hit: &mut KmerHit, query_length: i64, target_length: i64) {
    let span = hit.count & 0xff;
    let query_position = hit.self_offset.wrapping_add(1).wrapping_sub(span);
    hit.self_offset = (query_length - 1 - i64::from(query_position)) as u32;
    let target_position = hit.offset.wrapping_add(1).wrapping_sub(span);
    hit.offset = (target_length - 1 - i64::from(target_position)) as u32;
}

/// `lchain_qdp_fix`: the window-rechaining dynamic programme.
///
/// When `right_fixed` is set without `left_fixed` the anchor array is mirrored
/// so that the fixed end becomes the left one, the chain is computed, and the
/// array and the traceback indices are mapped back. `left_fixed` and
/// `right_fixed` also make the chain start (respectively end) on the first
/// (last) anchor instead of the best scoring one.
#[allow(clippy::too_many_arguments)]
pub fn chain_dp_fixed_ends(
    hits: &mut [KmerHit],
    hit_count: i64,
    data: &mut ChainData,
    maximum_skip: i64,
    maximum_iterations: i64,
    maximum_distance: i64,
    gap_penalty: f64,
    skip_penalty: f64,
    band_width_rate: f64,
    query_length: i64,
    target_length: i64,
    quick_check: bool,
    left_fixed: bool,
    right_fixed: bool,
) -> i64 {
    let mut anchor_count = hit_count;
    data.resize(anchor_count as usize);

    let mut best_score: i64 = -1;
    let mut best_index: i64 = -1;
    let mut shortest_overlap: i64 = i64::from(i32::MAX);
    let mut quick_success = false;

    if quick_check {
        let accepted = chain_check_quick(
            &hits[..anchor_count as usize],
            anchor_count as i32,
            data,
            band_width_rate,
        );
        if let Some(accepted_count) = accepted {
            anchor_count = i64::from(accepted_count);
            best_index = anchor_count - 1;
            best_score = i64::from(data.score[best_index as usize]);
            quick_success = true;
        }
    }

    let mut is_reordered = false;
    if !quick_success {
        data.scratch[..anchor_count as usize].fill(0);

        if right_fixed && !left_fixed {
            let count = anchor_count;
            let half = count >> 1;
            for index in 0..half {
                let mirror = count - index - 1;
                hits.swap(index as usize, mirror as usize);
                reverse_kmer_hit(&mut hits[index as usize], query_length, target_length);
                reverse_kmer_hit(&mut hits[mirror as usize], query_length, target_length);
            }
            if count & 1 == 1 {
                reverse_kmer_hit(&mut hits[half as usize], query_length, target_length);
            }
            is_reordered = true;
        }

        // Without a fixed end the chain may start anywhere, so every anchor is
        // a valid head; with one, only the sentinel may start the chain.
        let must_use_prefix = left_fixed || right_fixed;

        let mut start = 0i64;
        let mut maximum_index: i64 = -1;
        for index in 0..anchor_count {
            let span = i64::from(hits[index as usize].count & 0xff);
            let mut best = if must_use_prefix {
                i64::from(i32::MIN)
            } else {
                span
            };
            let mut skip_count = 0i64;
            let mut best_predecessor: i64 = -1;

            if index - start > maximum_iterations {
                start = index - maximum_iterations;
            }

            let mut predecessor = index - 1;
            while predecessor >= 0 {
                let link = fixed_end_link_score(
                    &hits[index as usize],
                    &hits[predecessor as usize],
                    band_width_rate,
                    gap_penalty,
                    skip_penalty,
                );
                if link == INVALID_SCORE {
                    predecessor -= 1;
                    continue;
                }

                let candidate = i64::from(link) + i64::from(data.score[predecessor as usize]);
                if candidate > best {
                    best = candidate;
                    best_predecessor = predecessor;
                    if skip_count > 0 {
                        skip_count -= 1;
                    }
                } else if data.scratch[predecessor as usize] == i64::from(index as i32) {
                    skip_count += 1;
                    if skip_count > maximum_skip
                        && (best_predecessor != -1 || !must_use_prefix)
                    {
                        break;
                    }
                }

                let predecessor_of_predecessor = data.predecessor[predecessor as usize];
                if predecessor_of_predecessor >= 0 {
                    data.scratch[predecessor_of_predecessor as usize] = index;
                }

                if i64::from(hits[index as usize].self_offset)
                    > maximum_distance + i64::from(hits[predecessor as usize].self_offset)
                {
                    if best_predecessor != -1 {
                        break;
                    }
                }

                if predecessor < start {
                    if best_predecessor != -1 || !must_use_prefix {
                        break;
                    }
                }

                predecessor -= 1;
            }
            let end_of_scan = predecessor;

            // The best anchor within `max_dis` is tried as a fallback when the
            // band scan broke early.
            if maximum_index < 0
                || i64::from(hits[index as usize].self_offset)
                    - i64::from(hits[maximum_index as usize].self_offset)
                    > maximum_distance
            {
                let mut maximum_score = i32::MIN;
                maximum_index = -1;
                let mut candidate = index - 1;
                while candidate >= start
                    && i64::from(hits[index as usize].self_offset)
                        - i64::from(hits[candidate as usize].self_offset)
                        <= maximum_distance
                {
                    if maximum_score < data.score[candidate as usize] {
                        maximum_score = data.score[candidate as usize];
                        maximum_index = candidate;
                    }
                    candidate -= 1;
                }
            }

            if maximum_index >= 0 && maximum_index < end_of_scan {
                let link = fixed_end_link_score(
                    &hits[index as usize],
                    &hits[maximum_index as usize],
                    band_width_rate,
                    gap_penalty,
                    skip_penalty,
                );
                let candidate = link.wrapping_add(data.score[maximum_index as usize]);
                if link != INVALID_SCORE && best < i64::from(candidate) {
                    best = i64::from(candidate);
                    best_predecessor = maximum_index;
                }
            }

            if best_predecessor == -1 {
                data.score[index as usize] = 0;
            } else {
                data.score[index as usize] = best as i32;
            }
            data.predecessor[index as usize] = best_predecessor;

            if maximum_index < 0
                || (i64::from(hits[index as usize].self_offset)
                    - i64::from(hits[maximum_index as usize].self_offset)
                    <= maximum_distance
                    && data.score[maximum_index as usize] < data.score[index as usize])
            {
                maximum_index = index;
            }

            if i64::from(data.score[index as usize]) >= best_score {
                let overlap = chain_length(
                    i64::from(hits[index as usize].self_offset),
                    i64::from(hits[index as usize].self_offset),
                    query_length,
                    i64::from(hits[index as usize].offset),
                    i64::from(hits[index as usize].offset),
                    target_length,
                );
                if i64::from(data.score[index as usize]) > best_score || overlap < shortest_overlap
                {
                    best_score = i64::from(data.score[index as usize]);
                    best_index = index;
                    shortest_overlap = overlap;
                }
            }
        }
    }

    if right_fixed && left_fixed {
        best_index = anchor_count - 1;
    }

    // Trace the chosen chain. `scratch` doubles as the reference's `tmp`
    // array: it holds the anchor index of every chain step, in reverse order
    // until the final reversal below.
    let mut chain_length_count: i64 = 0;
    let mut index = best_index;
    while index >= 0 {
        data.scratch[chain_length_count as usize] = index;
        chain_length_count += 1;
        index = data.predecessor[index as usize];
    }

    if is_reordered {
        let count = anchor_count;
        let half = count >> 1;
        for index in 0..half {
            let mirror = count - index - 1;
            hits.swap(index as usize, mirror as usize);
            reverse_kmer_hit(&mut hits[index as usize], query_length, target_length);
            reverse_kmer_hit(&mut hits[mirror as usize], query_length, target_length);
        }
        if count & 1 == 1 {
            reverse_kmer_hit(&mut hits[half as usize], query_length, target_length);
        }
        for position in 0..chain_length_count {
            data.scratch[position as usize] = count - data.scratch[position as usize] - 1;
        }
    } else {
        let half = chain_length_count >> 1;
        for position in 0..half {
            let mirror = chain_length_count - position - 1;
            data.scratch.swap(position as usize, mirror as usize);
        }
    }

    chain_length_count
}

/// `lchain_simple`: the simple chaining used when the exact match runs are
/// already almost colinear.
///
/// The anchors are reordered in place, exactly like the reference, whose
/// `des` argument always aliases `a` at this call site.
pub fn chain_simple(
    hits: &mut [KmerHit],
    hit_count: i64,
    data: &mut ChainData,
    maximum_skip: i64,
    maximum_iterations: i64,
) -> i64 {
    if hit_count <= 0 {
        return 0;
    }

    data.resize(hit_count as usize);
    let mut best_index = hit_count - 1;

    data.score[0] = hits[0].count as i32;
    data.predecessor[0] = -1;
    let mut index = 1i64;
    while index < hit_count {
        let previous = index - 1;
        if hits[index as usize].self_offset > hits[previous as usize].self_offset
            && hits[index as usize].offset > hits[previous as usize].offset
        {
            data.predecessor[index as usize] = previous;
            // `f[j] + a[i].cnt` is an unsigned 32-bit addition in the
            // reference, because `cnt` is `uint32_t`.
            data.score[index as usize] = (data.score[previous as usize] as u32)
                .wrapping_add(hits[index as usize].count)
                as i32;
        } else {
            break;
        }
        index += 1;
    }

    if index < hit_count {
        data.scratch[..hit_count as usize].fill(0);
        data.score[0] = hits[0].count as i32;
        data.predecessor[0] = -1;
        let mut best_score = i64::from(data.score[0]);
        best_index = 0;
        let mut start = 0i64;

        for index in 1..hit_count {
            let mut best = i64::from(i32::MIN);
            let mut skip_count = 0i64;
            let mut best_predecessor: i64 = -1;

            if index - start > maximum_iterations {
                start = index - maximum_iterations;
            }

            let mut predecessor = index - 1;
            while predecessor >= start {
                if hits[index as usize].self_offset > hits[predecessor as usize].self_offset
                    && hits[index as usize].offset > hits[predecessor as usize].offset
                {
                    let candidate = i64::from(
                        (data.score[predecessor as usize] as u32)
                            .wrapping_add(hits[index as usize].count),
                    );
                    if candidate > best {
                        best = candidate;
                        best_predecessor = predecessor;
                        if skip_count > 0 {
                            skip_count -= 1;
                        }
                    } else if data.scratch[predecessor as usize] == i64::from(index as i32) {
                        skip_count += 1;
                        if skip_count > maximum_skip {
                            break;
                        }
                    }

                    let predecessor_of_predecessor = data.predecessor[predecessor as usize];
                    if predecessor_of_predecessor >= 0 {
                        data.scratch[predecessor_of_predecessor as usize] = index;
                    }
                }
                predecessor -= 1;
            }

            data.score[index as usize] = best as i32;
            data.predecessor[index as usize] = best_predecessor;
            if i64::from(data.score[index as usize]) > best_score {
                best_score = i64::from(data.score[index as usize]);
                best_index = index;
            }
        }
    }

    let mut chain_count: i64 = 0;
    let mut index = best_index;
    while index >= 0 {
        data.scratch[chain_count as usize] = index;
        chain_count += 1;
        index = data.predecessor[index as usize];
    }

    let half = chain_count >> 1;
    for position in 0..half {
        let mirror = chain_count - position - 1;
        data.scratch.swap(position as usize, mirror as usize);
    }

    // `des[i] = a[t[i]]` with `des == a`: the traceback indices are strictly
    // increasing, so a forward copy never reads a slot it has overwritten.
    let reordered: Vec<KmerHit> = (0..chain_count as usize)
        .map(|position| hits[data.scratch[position] as usize])
        .collect();
    hits[..chain_count as usize].copy_from_slice(&reordered);

    chain_count
}

/// `lchain_refine`: chains an already x-sorted hit list with the lightweight
/// DP used by `return_t_chain`.
///
/// Consecutive hits are linked by default; a diagonal gap up to `long_gap` is
/// free and a larger one costs its size. Only hits within `maximum_distance`
/// of the current one are considered by the backward scan. The best chain is
/// traced back and written in chain order to the front of `hits`; the return
/// value is its length.
pub fn chain_refine(
    hits: &mut [KmerHit],
    hit_count: i64,
    data: &mut ChainData,
    maximum_skip: i64,
    maximum_iterations: i64,
    maximum_distance: i64,
    long_gap: i64,
) -> i64 {
    if hit_count <= 0 {
        return 0;
    }

    data.resize(hit_count as usize);
    let mut best_index = hit_count - 1;

    // The fast path trusts the list while consecutive hits are linked by a
    // free diagonal gap or are too far apart in the query to be compared.
    data.score[0] = 0;
    data.predecessor[0] = -1;
    let mut index = 1i64;
    while index < hit_count {
        let previous = index - 1;
        let query_gap = i64::from(hits[index as usize].self_offset)
            - i64::from(hits[previous as usize].self_offset);
        let target_gap = i64::from(hits[index as usize].offset)
            - i64::from(hits[previous as usize].offset);
        let diagonal_gap = if target_gap > query_gap {
            target_gap - query_gap
        } else {
            query_gap - target_gap
        };
        if diagonal_gap <= long_gap || query_gap > maximum_distance {
            data.predecessor[index as usize] = previous;
            data.score[index as usize] = index as i32;
        } else {
            break;
        }
        index += 1;
    }

    if index < hit_count {
        data.scratch[..hit_count as usize].fill(0);
        data.score[0] = 0;
        data.predecessor[0] = -1;
        let mut start = 0i64;

        for index in 1..hit_count {
            let mut best = i64::from(i32::MIN);
            let mut skip_count = 0i64;
            let mut best_predecessor: i64 = -1;

            if index - start > maximum_iterations {
                start = index - maximum_iterations;
            }

            // The immediate predecessor is always scored; unlike the backward
            // scan it is not bounded by `maximum_distance`.
            let mut predecessor = index - 1;
            let query_gap = i64::from(hits[index as usize].self_offset)
                - i64::from(hits[predecessor as usize].self_offset);
            let target_gap = i64::from(hits[index as usize].offset)
                - i64::from(hits[predecessor as usize].offset);
            let mut diagonal_gap = if target_gap > query_gap {
                target_gap - query_gap
            } else {
                query_gap - target_gap
            };
            if diagonal_gap <= long_gap {
                diagonal_gap = 0;
            }
            if i64::from(data.score[predecessor as usize]) - diagonal_gap > best {
                best = i64::from(data.score[predecessor as usize]) - diagonal_gap;
                best_predecessor = predecessor;
            }
            let predecessor_of_predecessor = data.predecessor[predecessor as usize];
            if predecessor_of_predecessor >= 0 {
                data.scratch[predecessor_of_predecessor as usize] = index;
            }

            predecessor -= 1;
            while predecessor >= start
                && i64::from(hits[index as usize].self_offset)
                    <= maximum_distance + i64::from(hits[predecessor as usize].self_offset)
            {
                let query_gap = i64::from(hits[index as usize].self_offset)
                    - i64::from(hits[predecessor as usize].self_offset);
                let target_gap = i64::from(hits[index as usize].offset)
                    - i64::from(hits[predecessor as usize].offset);
                let mut diagonal_gap = if target_gap > query_gap {
                    target_gap - query_gap
                } else {
                    query_gap - target_gap
                };
                if diagonal_gap <= long_gap {
                    diagonal_gap = 0;
                }
                let candidate = i64::from(data.score[predecessor as usize]) - diagonal_gap;
                if candidate > best {
                    best = candidate;
                    best_predecessor = predecessor;
                    if skip_count > 0 {
                        skip_count -= 1;
                    }
                } else if data.scratch[predecessor as usize] == i64::from(index as i32) {
                    skip_count += 1;
                    if skip_count > maximum_skip {
                        break;
                    }
                }
                let predecessor_of_predecessor = data.predecessor[predecessor as usize];
                if predecessor_of_predecessor >= 0 {
                    data.scratch[predecessor_of_predecessor as usize] = index;
                }
                predecessor -= 1;
            }

            data.score[index as usize] = best as i32;
            data.predecessor[index as usize] = best_predecessor;
        }

        // The chain does not have to end at the final hit: any hit within
        // `maximum_distance` of it with a better score can close the chain.
        index = hit_count - 1;
        let mut best_score = i64::from(data.score[index as usize]);
        best_index = index;
        let mut predecessor = index - 1;
        while predecessor >= 0
            && i64::from(hits[index as usize].self_offset)
                <= maximum_distance + i64::from(hits[predecessor as usize].self_offset)
        {
            if best_score < i64::from(data.score[predecessor as usize])
                && data.predecessor[predecessor as usize] >= 0
            {
                best_score = i64::from(data.score[predecessor as usize]);
                best_index = predecessor;
            }
            predecessor -= 1;
        }
    }

    let mut chain_count: i64 = 0;
    let mut index = best_index;
    while index >= 0 {
        data.scratch[chain_count as usize] = index;
        chain_count += 1;
        index = data.predecessor[index as usize];
    }

    let half = chain_count >> 1;
    for position in 0..half {
        let mirror = chain_count - position - 1;
        data.scratch.swap(position as usize, mirror as usize);
    }

    // `des[i] = a[t[i]]` with `des == a`: the traceback indices are strictly
    // increasing, so a forward copy never reads a slot it has overwritten.
    let reordered: Vec<KmerHit> = (0..chain_count as usize)
        .map(|position| hits[data.scratch[position] as usize])
        .collect();
    hits[..chain_count as usize].copy_from_slice(&reordered);

    chain_count
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(offset: u32, self_offset: u32, count: u32) -> KmerHit {
        KmerHit {
            read_identifier: 1,
            reverse_strand: false,
            offset,
            self_offset,
            count,
        }
    }

    #[test]
    fn collinear_hits_take_the_fast_path() {
        let hits = vec![
            hit(100, 100, 1),
            hit(200, 200, 1),
            hit(300, 300, 1),
            hit(400, 400, 1),
            hit(500, 500, 1),
        ];
        let mut data = ChainData::default();
        let mut result = OverlapRegion::default();
        let chain_size = chain_dp(
            &hits,
            5,
            &mut data,
            &mut result,
            0.05,
            5,
            1000,
            1000,
            51,
        );

        assert_eq!(chain_size, 5);
        assert_eq!(result.x_position_start, 100);
        assert_eq!(result.x_position_end, 500);
        assert_eq!(result.y_position_start, 100);
        assert_eq!(result.y_position_end, 500);
    }

    /// Golden values printed by the real `chain_DP` object file, driven with
    /// the same anchors and `k_mer_length = 51`.
    #[test]
    fn chain_dp_matches_the_reference_implementation() {
        struct Case {
            name: &'static str,
            hits: Vec<KmerHit>,
            strand: u32,
            expected_chain: i64,
            expected_start: (u32, u32),
            expected_end: (u32, u32),
            expected_overlap: u32,
            expected_seed: i32,
            expected_cigar_length: usize,
            read_length: i64,
        }

        let cases = vec![
            Case {
                name: "collinear",
                hits: vec![
                    hit(100, 100, 1),
                    hit(200, 200, 1),
                    hit(300, 300, 1),
                    hit(400, 400, 1),
                    hit(500, 500, 1),
                ],
                strand: 0,
                expected_chain: 5,
                expected_start: (100, 100),
                expected_end: (500, 500),
                expected_overlap: 1000,
                expected_seed: 255,
                expected_cigar_length: 1,
                read_length: 1000,
            },
            Case {
                name: "indels",
                hits: vec![
                    hit(100, 100, 1),
                    hit(200, 200, 1),
                    hit(301, 300, 1),
                    hit(400, 400, 2),
                    hit(502, 500, 1),
                ],
                strand: 0,
                expected_chain: 5,
                expected_start: (100, 100),
                expected_end: (500, 502),
                expected_overlap: 998,
                expected_seed: 211,
                expected_cigar_length: 4,
                read_length: 1000,
            },
            Case {
                name: "zigzag",
                hits: vec![
                    hit(100, 500, 1),
                    hit(200, 200, 1),
                    hit(300, 300, 3),
                    hit(400, 400, 1),
                    hit(500, 100, 1),
                ],
                strand: 0,
                expected_chain: 3,
                expected_start: (200, 200),
                expected_end: (400, 400),
                expected_overlap: 1000,
                expected_seed: 119,
                expected_cigar_length: 1,
                read_length: 1000,
            },
            Case {
                name: "reverse",
                hits: vec![
                    hit(100, 100, 2),
                    hit(250, 260, 2),
                    hit(400, 405, 2),
                    hit(550, 560, 2),
                ],
                strand: 1,
                expected_chain: 3,
                expected_start: (100, 100),
                expected_end: (560, 550),
                expected_overlap: 790,
                expected_seed: 57,
                expected_cigar_length: 3,
                read_length: 800,
            },
        ];

        for case in cases {
            let mut data = ChainData::default();
            let mut result = OverlapRegion::default();
            result.x_strand = case.strand;

            let chain_size = chain_dp(
                &case.hits,
                case.hits.len() as i64,
                &mut data,
                &mut result,
                0.05,
                5,
                case.read_length,
                case.read_length,
                51,
            );

            assert_eq!(chain_size, case.expected_chain, "chain of {}", case.name);
            assert_eq!(
                (result.x_position_start, result.y_position_start),
                case.expected_start,
                "start of {}",
                case.name
            );
            assert_eq!(
                (result.x_position_end, result.y_position_end),
                case.expected_end,
                "end of {}",
                case.name
            );
            assert_eq!(
                result.overlap_length, case.expected_overlap,
                "overlap of {}",
                case.name
            );
            assert_eq!(
                result.shared_seed, case.expected_seed,
                "seed of {}",
                case.name
            );
            assert_eq!(
                result.fake_cigar.len(),
                case.expected_cigar_length,
                "fake cigar of {}",
                case.name
            );
        }
    }

    #[test]
    fn candidate_groups_are_chained_per_read_and_strand() {
        let seeds = vec![
            // Anchors pointing at the query itself must be skipped.
            KmerHit { read_identifier: 0, reverse_strand: false, offset: 100, self_offset: 100, count: 1 },
            KmerHit { read_identifier: 0, reverse_strand: false, offset: 200, self_offset: 200, count: 1 },
            // Two anchors for read 5, forward.
            KmerHit { read_identifier: 5, reverse_strand: false, offset: 100, self_offset: 100, count: 1 },
            KmerHit { read_identifier: 5, reverse_strand: false, offset: 200, self_offset: 200, count: 1 },
            // One anchor for read 5, reverse.
            KmerHit { read_identifier: 5, reverse_strand: true, offset: 300, self_offset: 300, count: 1 },
            // Two anchors for read 7.
            KmerHit { read_identifier: 7, reverse_strand: false, offset: 100, self_offset: 100, count: 1 },
            KmerHit { read_identifier: 7, reverse_strand: false, offset: 400, self_offset: 400, count: 1 },
        ];
        let read_lengths = vec![1000u64; 8];

        let mut data = ChainData::default();
        let mut overlap_list: Vec<OverlapRegion> = Vec::new();
        chain_candidate_groups(
            &seeds,
            0,
            1000,
            &read_lengths,
            0.05,
            51,
            false,
            &mut data,
            &mut overlap_list,
        );

        // The two read-5 chains share one target identifier, so the append
        // step keeps only the better one, exactly like the reference.
        assert_eq!(overlap_list.len(), 2);
        assert_eq!(overlap_list[0].y_identifier, 5);
        assert_eq!(overlap_list[1].y_identifier, 7);
        assert_eq!(overlap_list[0].shared_seed, 102);
    }

    #[test]
    fn chain_check_rejects_backwards_hits() {
        let hits = vec![hit(100, 200, 1), hit(200, 100, 1)];
        let mut data = ChainData::default();
        data.resize(2);
        assert_eq!(chain_check(&hits, 2, &mut data, 51, 0.05), -1);
    }

    #[test]
    fn chain_dp_handles_an_empty_list_gracefully() {
        let hits: Vec<KmerHit> = Vec::new();
        let mut data = ChainData::default();
        let mut result = OverlapRegion::default();
        let chain_size = chain_dp(&hits, 0, &mut data, &mut result, 0.05, 5, 100, 100, 51);
        assert_eq!(chain_size, 0);
    }
}

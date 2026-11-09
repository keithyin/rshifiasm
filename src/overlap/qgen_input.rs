//! Seed generation and chaining over a prebuilt global minimizer index, the
//! counterparts of `minimizers_qgen_input` / `gen_pair_chain` and
//! `lchain_qgen_mcopy_input` / `special_lchain` in `anchor.cpp`.
//!
//! `ug_map_lchain` (`anchor.cpp` 3602) drives this pass when the unitig graph
//! is mapped against itself (`trans_base_infer`): the query minimizers carry
//! their position in the hash-sorted global index, and every shared hash
//! group becomes a set of anchors against the other unitigs.

use crate::assembly::hap_cov::UTrans;
use crate::assembly::trans_base::MzlEntry;
use crate::container::ksort::{
    introsort_klib_by, radix_sort_by_key, radix_sort_by_key_32, radix_sort_u64,
};
use crate::overlap::engine::overlap_type;
use crate::overlap::qgen_chain::chain_candidate_read_ch;
use crate::overlap::region::{KmerHit, OverlapRegion};

/// One raw anchor, `anchor1_t`: `srt` packs `(rid << 33) | (rev << 32) |
/// self_off` and `count` packs `(group width << 8) | span`.
#[derive(Debug, Clone, Copy, Default)]
struct Anchor {
    srt: u64,
    self_off: u32,
    other_off: u32,
    count: u32,
}

/// Parameters of `lchain_qgen_mcopy_input` on the `ug_map_lchain` path.
///
/// `enable_mcopy` is hard-wired to zero and `khit_n` to zero at every call
/// site, so neither the copy-number nor the anchor-count flag appears here.
#[derive(Debug, Clone, Copy)]
pub(crate) struct InputChainParameters {
    pub maximum_chain_count: u64,
    pub maximum_skip: i64,
    pub maximum_iterations: i64,
    pub maximum_distance: i64,
    pub gap_penalty: f64,
    pub skip_penalty: f64,
    /// `bw_thres`, the band-width rate of the plain groups.
    pub band_width_rate: f64,
    /// `bw_thres_sec`, the rate `special_lchain` uses for both of its calls.
    pub band_width_rate_sec: f64,
    pub append_both_ends: bool,
    /// `gen_off`, forwarded as `gen_cigar`.
    pub generate_offset: bool,
}

/// The query-side minimizers and the global index of `minimizers_qgen_input`.
pub(crate) struct SeedInput<'a> {
    pub query_identifier: u32,
    /// `in`: the query's minimizers, `x` remapped to the sorted position.
    pub query_minimizers: &'a [MzlEntry],
    /// `idx`: the hash-sorted global index, `x` holding the original hash.
    pub index: &'a [MzlEntry],
    /// `mzl_cutoff` (`mini_cut` after the quantile adjustment).
    pub minimizer_cutoff: u64,
    /// `chain_cutoff`: minimum anchor count of a kept target group.
    pub chain_cutoff: u64,
    /// `kov`'s transitions of the query, deleted entries included.
    pub known_transitions: &'a [UTrans],
    /// `udb->ug->u.a[tid].len` for every target unitig.
    pub target_lengths: &'a [u32],
    /// `high_occ`; `None` keeps the reference default of `UINT32_MAX`.
    pub high_occurrence: Option<u32>,
    /// `low_occ`; `None` keeps the reference default of zero.
    pub low_occurrence: Option<u32>,
}

/// Counting scan of one query minimizer: walks the same-hash group right then
/// left, counting entries of other unitigs until the cutoff is reached.
/// Returns the count and the group boundaries `[ns, ne)` the scans stopped at.
fn scan_group(
    index: &[MzlEntry],
    position: u64,
    query_identifier: u32,
    cutoff: u64,
    known: Option<&[u64]>,
) -> (u64, usize, usize) {
    let hash = index[position as usize].x;
    let reverse = index[position as usize].reverse;
    let mut count = 0u64;
    let mut right = position as i64 + 1;
    while right < index.len() as i64 && index[right as usize].x == hash && count < cutoff {
        let entry = &index[right as usize];
        if entry.rid != query_identifier && group_member(entry, reverse, known) {
            count += 1;
        }
        right += 1;
    }
    let mut left = position as i64 - 1;
    while left >= 0 && index[left as usize].x == hash && count < cutoff {
        let entry = &index[left as usize];
        if entry.rid != query_identifier && group_member(entry, reverse, known) {
            count += 1;
        }
        left -= 1;
    }
    (count, (left + 1) as usize, right as usize)
}

/// Whether an index entry counts towards the known-target list of
/// `gen_pair_chain`; `None` accepts every entry (`minimizers_qgen_input`).
fn group_member(entry: &MzlEntry, query_reverse: bool, known: Option<&[u64]>) -> bool {
    match known {
        None => true,
        Some(targets) => {
            let key = (u64::from(entry.rid) << 1) | u64::from(entry.reverse != query_reverse);
            targets.contains(&key)
        }
    }
}

/// Emission scan of one query minimizer: appends an anchor for every
/// same-hash entry of another unitig (restricted to `known` when present),
/// walking right then left without the counting cutoff.
fn emit_group(
    anchors: &mut Vec<Anchor>,
    index: &[MzlEntry],
    position: u64,
    query_identifier: u32,
    group_width: u32,
    known: Option<&[u64]>,
) {
    let query = &index[position as usize];
    let hash = query.x;
    let mut emit = |entry: &MzlEntry| {
        if entry.rid == query_identifier || !group_member(entry, query.reverse, known) {
            return;
        }
        let reverse = query.reverse != entry.reverse;
        let other_off = if reverse {
            0xffff_fffeu32.wrapping_sub(
                entry
                    .pos
                    .wrapping_add(1)
                    .wrapping_sub(u32::from(entry.span)),
            )
        } else {
            entry.pos
        };
        anchors.push(Anchor {
            srt: (u64::from(entry.rid) << 33)
                | (u64::from(reverse) << 32)
                | u64::from(query.pos),
            self_off: query.pos,
            other_off,
            count: u32::from(query.span).min(0xff) | group_width,
        });
    };
    let mut right = position as usize + 1;
    while right < index.len() && index[right].x == hash {
        emit(&index[right]);
        right += 1;
    }
    let mut left = position as i64 - 1;
    while left >= 0 && index[left as usize].x == hash {
        emit(&index[left as usize]);
        left -= 1;
    }
}

/// `gen_pair_chain`: adds anchors against the known targets that the main
/// group filter dropped, and reports the start of every target group it
/// produced in `group_starts`.
fn gen_pair_chain(
    anchors: &mut Vec<Anchor>,
    query_identifier: u32,
    known: &[u64],
    query_minimizers: &[MzlEntry],
    index: &[MzlEntry],
    cutoff: u64,
    group_starts: &mut Vec<u64>,
) {
    if known.is_empty() {
        return;
    }
    let mut total = 0u64;
    for query in query_minimizers {
        let (count, _, _) = scan_group(
            index,
            query.x,
            query_identifier,
            cutoff,
            Some(known),
        );
        if count == 0 || count >= cutoff {
            continue;
        }
        total += count;
    }
    if total == 0 {
        return;
    }

    let base = anchors.len();
    anchors.reserve(total as usize);
    for query in query_minimizers {
        let (count, ns, ne) = scan_group(
            index,
            query.x,
            query_identifier,
            cutoff,
            Some(known),
        );
        if count == 0 || count >= cutoff {
            continue;
        }
        let width = ((ne - ns) as u64).min(0xff_ffff) as u32;
        emit_group(
            anchors,
            index,
            query.x,
            query_identifier,
            width << 8,
            Some(known),
        );
    }
    debug_assert_eq!(anchors.len() - base, total as usize);

    radix_sort_by_key(&mut anchors[base..], &|anchor: &Anchor| anchor.srt);
    let mut cursor = 0usize;
    let mut start = base;
    for end in base + 1..=anchors.len() {
        if end == anchors.len() || (anchors[end].srt >> 32) != (anchors[start].srt >> 32) {
            let key = anchors[start].srt >> 32;
            while cursor < known.len() && known[cursor] < key {
                cursor += 1;
            }
            if cursor < known.len() && known[cursor] == key {
                group_starts.push(start as u64);
            }
            while cursor < known.len() && known[cursor] == key {
                cursor += 1;
            }
            start = end;
        }
    }
}

/// Result of [`minimizers_qgen_input`].
pub(crate) struct SeedAnchors {
    /// `cl->list`: the anchors grouped by target and strand, weighted.
    pub candidates: Vec<KmerHit>,
    /// `sp`: anchor-array starts of the known-target groups, ascending.
    pub group_starts: Vec<u64>,
}

/// `minimizers_qgen_input` (`anchor.cpp` 1395): turns the query minimizers
/// into the weighted anchor list the chaining pass consumes.
pub(crate) fn minimizers_qgen_input(input: &SeedInput<'_>) -> SeedAnchors {
    let rid = input.query_identifier;
    let index = input.index;
    let cutoff = input.minimizer_cutoff;
    let max_count = u64::from(input.high_occurrence.unwrap_or(u32::MAX).max(2));
    let min_count = u64::from(input.low_occurrence.unwrap_or(0).max(2));

    // First pass: size the anchor array.
    let mut total = 0u64;
    for query in input.query_minimizers {
        debug_assert_eq!(query.pos, index[query.x as usize].pos);
        debug_assert_eq!(query.rid, index[query.x as usize].rid);
        let (count, _, _) = scan_group(index, query.x, rid, cutoff, None);
        if count == 0 || count >= cutoff {
            continue;
        }
        total += count;
    }

    // Second pass: emit the anchors, right side of the group then left.
    let mut anchors: Vec<Anchor> = Vec::with_capacity(total as usize);
    for query in input.query_minimizers {
        let (count, ns, ne) = scan_group(index, query.x, rid, cutoff, None);
        if count == 0 || count >= cutoff {
            continue;
        }
        let width = ((ne - ns) as u64).min(0xff_ffff) as u32;
        emit_group(&mut anchors, index, query.x, rid, width << 8, None);
    }
    debug_assert_eq!(anchors.len() as u64, total);

    // The known-target keys, sorted like `radix_sort_anc64`.
    let mut known: Vec<u64> = input
        .known_transitions
        .iter()
        .filter(|transition| !transition.deleted)
        .map(|transition| {
            (u64::from(transition.target_identifier) << 1) | u64::from(transition.reverse)
        })
        .collect();
    if known.len() > 1 {
        radix_sort_u64(&mut known);
    }
    let known_count = known.len();

    // Keep the groups with at least `chain_cutoff` anchors; the starts of the
    // kept known-target groups are collected in emission order.
    radix_sort_by_key(&mut anchors, &|anchor: &Anchor| anchor.srt);
    let mut group_starts: Vec<u64> = Vec::new();
    let mut cursor = 0usize;
    let mut kept = 0usize;
    let mut start = 0usize;
    for end in 1..=anchors.len() {
        if end == anchors.len() || (anchors[end].srt >> 32) != (anchors[start].srt >> 32) {
            if (end - start) as u64 >= input.chain_cutoff {
                let key = anchors[start].srt >> 32;
                while cursor < known_count && known[cursor] < key {
                    cursor += 1;
                }
                if cursor < known_count && known[cursor] == key {
                    group_starts.push(kept as u64);
                }
                while cursor < known_count && known[cursor] == key {
                    known[cursor] = u64::MAX;
                    cursor += 1;
                }
                for index in start..end {
                    anchors[kept] = anchors[index];
                    kept += 1;
                }
            }
            start = end;
        }
    }
    anchors.truncate(kept);

    // The keys no kept group matched go through `gen_pair_chain`, whose
    // group starts follow the ones collected above.
    let unmatched: Vec<u64> = known[..known_count]
        .iter()
        .copied()
        .filter(|&key| key != u64::MAX)
        .collect();
    if !unmatched.is_empty() {
        gen_pair_chain(
            &mut anchors,
            rid,
            &unmatched,
            input.query_minimizers,
            index,
            cutoff,
            &mut group_starts,
        );
    }
    if known_count == 0 {
        group_starts.clear();
    }

    // Copy over to the candidate list, weighting the group occupancy.
    let mut candidates: Vec<KmerHit> = Vec::with_capacity(anchors.len());
    let mut current_target = u64::MAX;
    let mut target_length = 0u64;
    let mut start = 0usize;
    for end in 1..=anchors.len() {
        if end == anchors.len() || anchors[end].srt != anchors[start].srt {
            if end - start > 1 {
                radix_sort_by_key_32(&mut anchors[start..end], &|anchor: &Anchor| {
                    anchor.other_off
                });
            }
            let target = anchors[start].srt >> 33;
            if target != current_target {
                current_target = target;
                target_length = u64::from(input.target_lengths[target as usize]);
            }
            for anchor in &anchors[start..end] {
                let strand = ((anchor.srt >> 32) & 1) != 0;
                let offset = if strand {
                    target_length.wrapping_sub(u64::from(
                        0xffff_ffffu32.wrapping_sub(anchor.other_off),
                    )) as u32
                } else {
                    anchor.other_off
                };
                let occupancy = u64::from(anchor.count >> 8);
                let mut weight: u32 = if occupancy < max_count && occupancy > min_count {
                    1
                } else if occupancy <= min_count {
                    2
                } else {
                    let scaled =
                        1 + ((occupancy + (max_count << 1) - 1) / (max_count << 1));
                    (scaled as f64).powf(1.1) as u32
                };
                if weight > 0xff_ffff {
                    weight = 0xff_ffff;
                }
                candidates.push(KmerHit {
                    read_identifier: target as u32,
                    reverse_strand: strand,
                    offset,
                    self_offset: anchor.self_off,
                    count: (weight << 8) | (anchor.count & 0xff),
                });
            }
            start = end;
        }
    }

    SeedAnchors {
        candidates,
        group_starts,
    }
}

/// `special_lchain` (`anchor.cpp` 2115): chains a known-target group in two
/// strand runs, marks the run that produced overlaps at a reported group
/// start, and keeps only the marked overlaps plus the best one.
#[allow(clippy::too_many_arguments)]
fn special_lchain(
    candidates: &mut Vec<KmerHit>,
    overlaps: &mut Vec<OverlapRegion>,
    query_identifier: u32,
    query_length: i64,
    target_length: i64,
    parameters: &InputChainParameters,
    group_starts: &[u64],
    start_index: &mut usize,
    mut copied: usize,
    group_begin: usize,
    group_end: usize,
) -> usize {
    let overlap_begin = overlaps.len();
    let mut write = copied;
    if group_begin >= group_end {
        return copied;
    }
    let mut split = group_begin;
    while split < group_end
        && candidates[split].reverse_strand == candidates[group_begin].reverse_strand
    {
        split += 1;
    }

    // Both runs are chained at the secondary band width, exactly like the
    // reference's `bw_thres_sec` argument.
    let mut chain_run = |begin: usize, end: usize, copied: usize, start_index: &mut usize| {
        let overlap_mark = overlaps.len();
        let added = chain_candidate_read_ch(
            candidates,
            begin,
            end - begin,
            copied,
            overlaps,
            parameters.maximum_skip,
            parameters.maximum_iterations,
            parameters.maximum_distance,
            parameters.gap_penalty,
            parameters.skip_penalty,
            parameters.band_width_rate_sec,
            query_identifier,
            query_length,
            target_length,
            parameters.append_both_ends,
            parameters.generate_offset,
        );
        if *start_index < group_starts.len()
            && group_starts[*start_index] == begin as u64
            && overlaps.len() > overlap_mark
        {
            overlaps.last_mut().unwrap().x_strand = 1;
        }
        while *start_index < group_starts.len()
            && group_starts[*start_index] == begin as u64
        {
            *start_index += 1;
        }
        added
    };

    if split > group_begin {
        copied += chain_run(group_begin, split, copied, start_index);
    }
    if split < group_end {
        copied += chain_run(split, group_end, copied, start_index);
    }

    if overlaps.len() > overlap_begin {
        let mut best = overlap_begin;
        for index in overlap_begin..overlaps.len() {
            if overlaps[best].shared_seed < overlaps[index].shared_seed
                || (overlaps[best].shared_seed == overlaps[index].shared_seed
                    && overlaps[index].x_strand == 1)
            {
                best = index;
            }
        }

        let mut kept = overlap_begin;
        for index in overlap_begin..overlaps.len() {
            if overlaps[index].x_strand == 0 && index != best {
                continue;
            }
            if kept != index {
                overlaps.swap(kept, index);
            }
            kept += 1;
        }

        if kept != overlaps.len() {
            overlaps.truncate(kept);
            for index in overlap_begin..overlaps.len() {
                let mut read = overlaps[index].non_homopolymer_errors as usize;
                let tag = candidates[read].read_identifier;
                overlaps[index].non_homopolymer_errors = write as u32;
                while read < copied && candidates[read].read_identifier == tag {
                    candidates[write] = candidates[read];
                    candidates[write].read_identifier = index as u32;
                    write += 1;
                    read += 1;
                }
            }
            copied = write;
        }
    }
    copied
}

/// `lchain_qgen_mcopy_input` (`anchor.cpp` 2197): groups the weighted anchors
/// by target, chains every group and applies the per-type chain limit.
///
/// Returns the overlaps sorted by query interval and the anchor lists of the
/// kept chains, tagged with the overlap index they belong to.
pub(crate) fn lchain_qgen_mcopy_input(
    candidates: Vec<KmerHit>,
    query_identifier: u32,
    query_length: u64,
    target_lengths: &[u32],
    group_starts: &[u64],
    parameters: &InputChainParameters,
) -> (Vec<OverlapRegion>, Vec<KmerHit>) {
    let mut candidates = candidates;
    let mut overlaps: Vec<OverlapRegion> = Vec::new();
    let candidate_count = candidates.len();
    let mut copied = 0usize;
    let mut start_index = 0usize;
    let mut start = 0usize;
    for end in 1..=candidate_count {
        if end == candidate_count
            || candidates[end].read_identifier != candidates[start].read_identifier
        {
            let target_identifier = candidates[start].read_identifier;
            if target_identifier != query_identifier {
                while start_index < group_starts.len()
                    && group_starts[start_index] < start as u64
                {
                    start_index += 1;
                }
                let target_length = i64::from(target_lengths[target_identifier as usize]);
                if start_index >= group_starts.len()
                    || group_starts[start_index] >= end as u64
                {
                    copied += chain_candidate_read_ch(
                        &mut candidates,
                        start,
                        end - start,
                        copied,
                        &mut overlaps,
                        parameters.maximum_skip,
                        parameters.maximum_iterations,
                        parameters.maximum_distance,
                        parameters.gap_penalty,
                        parameters.skip_penalty,
                        parameters.band_width_rate,
                        query_identifier,
                        query_length as i64,
                        target_length,
                        parameters.append_both_ends,
                        parameters.generate_offset,
                    );
                } else {
                    copied = special_lchain(
                        &mut candidates,
                        &mut overlaps,
                        query_identifier,
                        query_length as i64,
                        target_length,
                        parameters,
                        group_starts,
                        &mut start_index,
                        copied,
                        start,
                        end,
                    );
                }
            }
            start = end;
        }
    }
    candidates.truncate(copied);

    // Limit the number of chains per overlap-type class.
    let maximum_chain_count = parameters.maximum_chain_count;
    if overlaps.len() as u64 > maximum_chain_count {
        introsort_klib_by(&mut overlaps, &|left: &OverlapRegion, right: &OverlapRegion| {
            left.shared_seed > right.shared_seed
        });
        let mut type_counts = [0u32; 4];
        let mut type_thresholds = [0i32; 4];
        for region in &overlaps {
            let class = overlap_type(region, query_length as u32);
            type_counts[class] += 1;
            if u64::from(type_counts[class]) == maximum_chain_count {
                type_thresholds[class] = region.shared_seed;
            }
        }
        if type_thresholds.iter().any(|&threshold| threshold > 0) {
            let mut kept = 0usize;
            for index in 0..overlaps.len() {
                let class = overlap_type(&overlaps[index], query_length as u32);
                if overlaps[index].shared_seed >= type_thresholds[class]
                    || overlaps[index].x_strand == 1
                {
                    if kept != index {
                        overlaps.swap(kept, index);
                    }
                    kept += 1;
                }
            }
            overlaps.truncate(kept);
        }
    }
    introsort_klib_by(&mut overlaps, &|left: &OverlapRegion, right: &OverlapRegion| {
        ((u64::from(left.x_position_start) << 32) | u64::from(left.x_position_end))
            < ((u64::from(right.x_position_start) << 32) | u64::from(right.x_position_end))
    });

    (overlaps, candidates)
}
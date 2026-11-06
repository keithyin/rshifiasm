//! Cleanup of the collected read-to-read transitions, the counterparts of
//! `kt_u_trans_t_idx`, `get_u_trans_spec`, `merge_u_trans`,
//! `kt_u_trans_t_symm`, `filter_u_trans_t`, `clean_u_trans_t_idx` and
//! `set_trio_flag_by_cov` in `Overlaps.cpp`.
//!
//! The transition table is the reference's `kv_u_trans_t`: a flat
//! [`UTrans`] vector plus an index of `<start>:32 | <count>` spans, one per
//! query unitig. [`kt_u_trans_t_idx`] rebuilds the index after every sort or
//! compaction, exactly like the reference.
//!
//! [`clean_u_trans_t_idx`] makes the table symmetric, merges the transitions
//! of each unitig pair and drops the pairs that are not supported by at least
//! three fully covered reads. [`set_trio_flag_by_cov`] then assigns the
//! `FATHER`/`MOTHER` trio flags from the coverage of the partner unitigs.

use crate::assembly::asg::{Asg, MaUtg};
use crate::assembly::hap_cov::{get_offset_adjust, N_HET, UTrans};
use crate::assembly::purge_dups::SET_TRIO;
use crate::assembly::unitig::{AMBIGU, FATHER, MOTHER};
use crate::container::ksort::{radix_sort_by_key, radix_sort_by_key_32};

/// `u_trans_key`: the sort key `(qn << 32) | (tn << 1) | rev`.
fn u_trans_key(transition: &UTrans) -> u64 {
    (u64::from(transition.query_identifier) << 32)
        | (u64::from(transition.target_identifier) << 1)
        | u64::from(transition.reverse)
}

/// `u_trans_a` and `u_trans_n`: the `<start, count>` span of query `qn`.
fn u_trans_span(index: &[u64], qn: u32) -> (usize, usize) {
    let packed = index[qn as usize];
    ((packed >> 32) as usize, packed as u32 as usize)
}

/// `kt_u_trans_t_idx`: sorts the transitions by `u_trans_key` and rebuilds
/// the per-query spans. `n` is the number of unitig vertices.
pub fn kt_u_trans_t_idx(transitions: &mut Vec<UTrans>, index: &mut Vec<u64>, n: u32) {
    radix_sort_by_key(transitions.as_mut_slice(), &u_trans_key);

    index.clear();
    index.resize(n as usize, 0);

    let len = transitions.len();
    let mut st = 0usize;
    for i in 1..=len {
        if i == len || transitions[i].query_identifier != transitions[st].query_identifier {
            index[transitions[st].query_identifier as usize] =
                ((st as u64) << 32) | ((i - st) as u64);
            st = i;
        }
    }
}

/// `get_u_trans_spec`: the `<start, count>` of the `tn` subgroup inside the
/// `qn` span, or `None` when there is none.
///
/// Deleted entries are not skipped, exactly like the reference.
pub fn get_u_trans_spec(
    transitions: &[UTrans],
    index: &[u64],
    qn: u32,
    tn: u32,
) -> Option<(usize, usize)> {
    let (start, n) = u_trans_span(index, qn);
    let mut st = 0usize;
    for i in 1..=n {
        if i == n
            || transitions[start + i].target_identifier
                != transitions[start + st].target_identifier
        {
            if transitions[start + st].target_identifier == tn {
                return Some((start + st, i - st));
            }
            st = i;
        }
    }
    None
}

/// Removes the deleted entries of `transitions`, keeping their order.
fn compact_transitions(transitions: &mut Vec<UTrans>) {
    let mut write = 0usize;
    for read in 0..transitions.len() {
        if transitions[read].deleted {
            continue;
        }
        transitions[write] = transitions[read];
        write += 1;
    }
    transitions.truncate(write);
}

/// The flipped target interval of a reversed transition, `[len - te, len - ts)`.
fn flip_target(transition: &mut UTrans, length: u32) {
    let ts = transition.target_start;
    let te = transition.target_end.wrapping_sub(1);
    transition.target_start = length.wrapping_sub(te).wrapping_sub(1);
    transition.target_end = length
        .wrapping_sub(ts)
        .wrapping_sub(1)
        .wrapping_add(1);
}

/// `merge_u_trans`: merges the overlapping transitions of one direction and
/// returns the total weight of the survivors.
fn merge_u_trans(transitions: &mut Vec<UTrans>, unitig_graph: &Asg) -> f64 {
    radix_sort_by_key_32(transitions.as_mut_slice(), &|t: &UTrans| t.query_start);

    for i in 0..transitions.len() {
        if transitions[i].deleted || !transitions[i].reverse {
            continue;
        }
        let tn = transitions[i].target_identifier as usize;
        let length = unitig_graph.sequences[tn].length;
        flip_target(&mut transitions[i], length);
    }

    // Repeatedly merge every pair whose query intervals touch; each pass
    // compacts the table.
    loop {
        let mut is_found = false;
        let occ = transitions.len();
        for i in 0..occ {
            if transitions[i].deleted {
                continue;
            }
            for k in (i + 1)..occ {
                if transitions[k].deleted {
                    continue;
                }
                if transitions[i].query_end < transitions[k].query_start {
                    break;
                }
                if transitions[i].query_end >= transitions[k].query_start
                    && transitions[i].target_end >= transitions[k].target_start
                {
                    let overlap_query = {
                        let low = transitions[i].query_end.min(transitions[k].query_end);
                        let high = transitions[i].query_start.max(transitions[k].query_start);
                        if low > high {
                            low - high
                        } else {
                            0
                        }
                    };
                    let overlap_target = {
                        let low = transitions[i].target_end.min(transitions[k].target_end);
                        let high = transitions[i].target_start.max(transitions[k].target_start);
                        if low > high {
                            low - high
                        } else {
                            0
                        }
                    };

                    let query_len = transitions[k]
                        .query_end
                        .wrapping_sub(transitions[k].query_start);
                    let target_len = transitions[k]
                        .target_end
                        .wrapping_sub(transitions[k].target_start);
                    let ov_q = query_len.wrapping_sub(overlap_query);
                    let ov_t = target_len.wrapping_sub(overlap_target);

                    let weight_query =
                        (f64::from(ov_q) / f64::from(query_len)) * transitions[k].network_weight;
                    let weight_target =
                        (f64::from(ov_t) / f64::from(target_len)) * transitions[k].network_weight;

                    transitions[i].query_end =
                        transitions[i].query_end.max(transitions[k].query_end);
                    transitions[i].target_end =
                        transitions[i].target_end.max(transitions[k].target_end);
                    transitions[i].network_weight += weight_query.min(weight_target);
                    is_found = true;
                    transitions[k].deleted = true;
                }
            }
        }
        compact_transitions(transitions);
        if !is_found {
            break;
        }
    }

    // Trim the query overlaps so that no two survivors overlap.
    let occ = transitions.len();
    for i in 0..occ {
        if transitions[i].deleted {
            continue;
        }
        for k in (i + 1)..occ {
            if transitions[k].deleted {
                continue;
            }
            let overlap_query = {
                let low = transitions[i].query_end.min(transitions[k].query_end);
                let high = transitions[i].query_start.max(transitions[k].query_start);
                if low > high {
                    low - high
                } else {
                    0
                }
            };
            if overlap_query == 0 {
                break;
            }
            let query_len = transitions[k]
                .query_end
                .wrapping_sub(transitions[k].query_start);
            if overlap_query == query_len {
                transitions[k].deleted = true;
                continue;
            }
            let self_len = transitions[i]
                .query_end
                .wrapping_sub(transitions[i].query_start);
            if overlap_query == self_len {
                transitions[i].deleted = true;
                continue;
            }
            let overlap_target = get_offset_adjust(
                overlap_query,
                query_len,
                transitions[k]
                    .target_end
                    .wrapping_sub(transitions[k].target_start),
            );
            transitions[k].query_start =
                transitions[k].query_start.wrapping_add(overlap_query);
            transitions[k].target_start =
                transitions[k].target_start.wrapping_add(overlap_target);
            let remaining = transitions[k]
                .query_end
                .wrapping_sub(transitions[k].query_start);
            transitions[k].network_weight -=
                (f64::from(overlap_query) / f64::from(remaining)) * transitions[k].network_weight;
        }
    }
    compact_transitions(transitions);

    // Then trim the target overlaps.
    radix_sort_by_key_32(transitions.as_mut_slice(), &|t: &UTrans| t.target_start);
    let occ = transitions.len();
    for i in 0..occ {
        if transitions[i].deleted {
            continue;
        }
        for k in (i + 1)..occ {
            if transitions[k].deleted {
                continue;
            }
            let overlap_target = {
                let low = transitions[i].target_end.min(transitions[k].target_end);
                let high = transitions[i].target_start.max(transitions[k].target_start);
                if low > high {
                    low - high
                } else {
                    0
                }
            };
            if overlap_target == 0 {
                break;
            }
            let target_len = transitions[k]
                .target_end
                .wrapping_sub(transitions[k].target_start);
            if overlap_target == target_len {
                transitions[k].deleted = true;
                continue;
            }
            let self_len = transitions[i]
                .target_end
                .wrapping_sub(transitions[i].target_start);
            if overlap_target == self_len {
                transitions[i].deleted = true;
                continue;
            }
            let overlap_query = get_offset_adjust(
                overlap_target,
                target_len,
                transitions[k]
                    .query_end
                    .wrapping_sub(transitions[k].query_start),
            );
            transitions[k].target_start =
                transitions[k].target_start.wrapping_add(overlap_target);
            transitions[k].query_start =
                transitions[k].query_start.wrapping_add(overlap_query);
            let remaining = transitions[k]
                .target_end
                .wrapping_sub(transitions[k].target_start);
            transitions[k].network_weight -= (f64::from(overlap_target)
                / f64::from(remaining))
                * transitions[k].network_weight;
        }
    }

    let mut total = 0.0f64;
    let mut write = 0usize;
    for read in 0..transitions.len() {
        if transitions[read].deleted {
            continue;
        }
        total += transitions[read].network_weight;
        transitions[write] = transitions[read];
        write += 1;
    }
    transitions.truncate(write);

    for i in 0..transitions.len() {
        if transitions[i].deleted || !transitions[i].reverse {
            continue;
        }
        let tn = transitions[i].target_identifier as usize;
        let length = unitig_graph.sequences[tn].length;
        flip_target(&mut transitions[i], length);
    }

    total
}

/// `kt_u_trans_t_symm`: makes every transition pair symmetric, merges both
/// directions and appends the mirrored entries.
pub fn kt_u_trans_t_symm(
    transitions: &mut Vec<UTrans>,
    index: &mut Vec<u64>,
    unitig_graph: &Asg,
) {
    let query_count = index.len();
    let mut e0: Vec<UTrans> = Vec::new();
    let mut e1: Vec<UTrans> = Vec::new();

    for k in 0..query_count {
        let (span_start, n) = u_trans_span(index, k as u32);
        let mut st = 0usize;
        for i in 1..=n {
            let end_of_group = i == n
                || transitions[span_start + i].target_identifier
                    != transitions[span_start + st].target_identifier;
            if !end_of_group {
                continue;
            }

            let qn = transitions[span_start + st].query_identifier;
            let tn = transitions[span_start + st].target_identifier;
            let mirror = get_u_trans_spec(transitions, index, tn, qn);
            let mirror_count = mirror.map_or(0, |(_, count)| count);
            if i == st + 1 && mirror_count == 0 {
                st = i;
                continue;
            }

            e0.clear();
            e1.clear();
            for m in st..i {
                let entry = transitions[span_start + m];
                if entry.deleted {
                    continue;
                }
                if entry.reverse {
                    e1.push(entry);
                } else {
                    e0.push(entry);
                }
            }

            if let Some((mirror_start, _)) = mirror {
                for m in 0..mirror_count {
                    let entry = transitions[mirror_start + m];
                    if entry.deleted {
                        continue;
                    }
                    let mut copy = entry;
                    copy.query_identifier = entry.target_identifier;
                    copy.query_start = entry.target_start;
                    copy.query_end = entry.target_end;
                    copy.target_identifier = entry.query_identifier;
                    copy.target_start = entry.query_start;
                    copy.target_end = entry.query_end;
                    if copy.reverse {
                        e1.push(copy);
                    } else {
                        e0.push(copy);
                    }
                }
            }

            if e0.len() + e1.len() > 1 {
                for m in st..i {
                    transitions[span_start + m].deleted = true;
                }
                if let Some((mirror_start, _)) = mirror {
                    for m in 0..mirror_count {
                        transitions[mirror_start + m].deleted = true;
                    }
                }

                let weight0 = merge_u_trans(&mut e0, unitig_graph);
                let weight1 = merge_u_trans(&mut e1, unitig_graph);
                if weight0 >= weight1 {
                    transitions.extend_from_slice(&e0);
                } else {
                    transitions.extend_from_slice(&e1);
                }
            }

            st = i;
        }
    }

    compact_transitions(transitions);

    // Every surviving entry also gets its mirrored copy.
    let original_count = transitions.len();
    for i in 0..original_count {
        let entry = transitions[i];
        if entry.deleted {
            continue;
        }
        let mut copy = entry;
        copy.query_identifier = entry.target_identifier;
        copy.query_start = entry.target_start;
        copy.query_end = entry.target_end;
        copy.target_identifier = entry.query_identifier;
        copy.target_start = entry.query_start;
        copy.target_end = entry.query_end;
        transitions.push(copy);
    }

    kt_u_trans_t_idx(transitions, index, unitig_graph.sequences.len() as u32);
}

/// `filter_u_trans_t`: drops every transition pair not covered by `thres`
/// reads on both sides, then removes their mirrors.
pub fn filter_u_trans_t(
    transitions: &mut Vec<UTrans>,
    index: &mut Vec<u64>,
    unitigs: &[MaUtg],
    unitig_graph: &Asg,
    read_graph: &Asg,
    thres: u32,
) {
    // Count, for every query unitig, how many of its reads are fully covered
    // by each transition.
    for k in 0..index.len() {
        let (span_start, n) = u_trans_span(index, k as u32);
        if n == 0 {
            continue;
        }

        let mut counts = vec![0u32; n];
        let mut min = u64::from(transitions[span_start].query_start);
        let mut max = u64::from(transitions[span_start].query_end);
        for i in 0..n {
            counts[i] = 0;
            min = min.min(u64::from(transitions[span_start + i].query_start));
            max = max.max(u64::from(transitions[span_start + i].query_end));
        }

        let unitig = &unitigs[k];
        let mut pass = 0u64;
        let mut offset = 0u64;
        for j in 0..unitig.reads.len() {
            pass = 0;
            let packed = unitig.reads[j];
            let read = (packed >> 33) as usize;
            let r_beg = offset;
            let r_end = offset + u64::from(read_graph.sequences[read].length);
            offset += u64::from(packed as u32);
            if min >= r_end {
                continue;
            }
            if max <= r_beg {
                break;
            }
            for i in 0..n {
                let query_end = u64::from(transitions[span_start + i].query_end);
                if counts[i] >= thres || r_beg >= query_end {
                    pass += 1;
                    continue;
                }
                let query_start = u64::from(transitions[span_start + i].query_start);
                let ovlp = if r_end.min(query_end) > r_beg.max(query_start) {
                    r_end.min(query_end) - r_beg.max(query_start)
                } else {
                    0
                };
                if ovlp == r_end - r_beg {
                    counts[i] += 1;
                    if counts[i] >= thres {
                        pass += 1;
                    }
                }
            }

            if pass == n as u64 {
                for i in 0..n {
                    if counts[i] < thres {
                        transitions[span_start + i].deleted = true;
                    }
                }
                break;
            }
        }

        if pass < n as u64 {
            for i in 0..n {
                if counts[i] < thres {
                    transitions[span_start + i].deleted = true;
                }
            }
        }
    }

    // Delete the mirror of every deleted transition.
    for k in 0..index.len() {
        let (span_start, n) = u_trans_span(index, k as u32);
        let mut st = 0usize;
        for i in 1..=n {
            let end_of_group = i == n
                || transitions[span_start + i].target_identifier
                    != transitions[span_start + st].target_identifier;
            if !end_of_group {
                continue;
            }

            let qn = transitions[span_start + st].query_identifier;
            let tn = transitions[span_start + st].target_identifier;
            if let Some((mirror_start, mirror_count)) =
                get_u_trans_spec(transitions, index, tn, qn)
            {
                for m in st..i {
                    let entry = transitions[span_start + m];
                    if !entry.deleted {
                        continue;
                    }
                    for j in 0..mirror_count {
                        let mirror = transitions[mirror_start + j];
                        if mirror.target_identifier == entry.query_identifier
                            && mirror.query_identifier == entry.target_identifier
                            && mirror.target_start == entry.query_start
                            && mirror.target_end == entry.query_end
                            && mirror.query_start == entry.target_start
                            && mirror.query_end == entry.target_end
                            && mirror.network_weight == entry.network_weight
                            && mirror.reverse == entry.reverse
                        {
                            transitions[mirror_start + j].deleted = true;
                            break;
                        }
                    }
                }
            }

            st = i;
        }
    }

    compact_transitions(transitions);
    kt_u_trans_t_idx(transitions, index, unitig_graph.sequences.len() as u32);
}

/// `clean_u_trans_t_idx`: index, symmetrise and filter the transitions.
pub fn clean_u_trans_t_idx(
    transitions: &mut Vec<UTrans>,
    index: &mut Vec<u64>,
    unitigs: &[MaUtg],
    unitig_graph: &Asg,
    read_graph: &Asg,
) {
    kt_u_trans_t_idx(transitions, index, unitig_graph.sequences.len() as u32);
    kt_u_trans_t_symm(transitions, index, unitig_graph);
    filter_u_trans_t(transitions, index, unitigs, unitig_graph, read_graph, 3);
}

/// The `set_trio_flag_by_cov` coverage walk over one target unitig.
///
/// Returns `[father_or_mother, set_trio]` counts exactly like the reference.
fn count_target_reads(
    target: &MaUtg,
    transition: &UTrans,
    read_graph: &Asg,
    trio_flags: &[u8],
) -> [u32; 2] {
    let mut counts = [0u32, 0u32];
    let mut found = false;
    let mut offset = 0u64;
    for j in 0..target.reads.len() {
        let packed = target.reads[j];
        let read = (packed >> 33) as usize;
        let r_beg = offset;
        let r_end = offset + u64::from(read_graph.sequences[read].length);
        offset += u64::from(packed as u32);

        let ts = u64::from(transition.target_start);
        let te = u64::from(transition.target_end);
        let ovlp = if r_end.min(te) > r_beg.max(ts) {
            r_end.min(te) - r_beg.max(ts)
        } else {
            0
        };
        if found && ovlp == 0 {
            break;
        }
        if ovlp == 0 {
            continue;
        }
        found = true;
        if (ovlp as f64) <= (u64::from(read_graph.sequences[read].length) as f64) * 0.8 {
            continue;
        }

        let flag = trio_flags[read];
        if flag & (FATHER | MOTHER) != 0 {
            counts[0] += 1;
        }
        if flag & SET_TRIO != 0 {
            counts[1] += 1;
        }
    }
    counts
}

/// `set_trio_flag_by_cov`: assigns `FATHER`/`MOTHER` to the reads that lack a
/// trio label, using the labels of the reads on the partner unitigs.
pub fn set_trio_flag_by_cov(
    unitigs: &[MaUtg],
    read_graph: &Asg,
    is_r_het: &[u8],
    transitions: &[UTrans],
    index: &[u64],
    trio_flags: &mut [u8],
) {
    // Pass 1: collect the candidate unitigs, ordered by the number of reads
    // that carry a trio label.
    let mut candidates: Vec<u64> = Vec::new();
    for i in 0..index.len() {
        let qn = i;
        let unitig = &unitigs[qn];
        let mut k = 0usize;
        while k < unitig.reads.len() {
            if trio_flags[(unitig.reads[k] >> 33) as usize] & SET_TRIO == 0 {
                break;
            }
            k += 1;
        }
        if k >= unitig.reads.len() {
            continue;
        }

        let (span_start, n) = u_trans_span(index, qn as u32);
        let mut counts = [0u32, 0u32];
        for t in 0..n {
            let transition = transitions[span_start + t];
            let target = &unitigs[transition.target_identifier as usize];
            let target_counts = count_target_reads(target, &transition, read_graph, trio_flags);
            counts[0] += target_counts[0];
            counts[1] += target_counts[1];
        }

        if counts[1] > 0 {
            candidates.push((u64::from(u32::MAX - counts[0]) << 32) | qn as u64);
        }
    }

    // Pass 2: pick the parent with the larger count on the partner unitigs
    // and tag the unlabelled heterozygous reads.
    for &candidate in candidates.iter() {
        let qn = candidate as u32 as usize;
        let unitig = &unitigs[qn];
        let (span_start, n) = u_trans_span(index, qn as u32);

        let mut father = 0u32;
        let mut mother = 0u32;
        for t in 0..n {
            let transition = transitions[span_start + t];
            let target = &unitigs[transition.target_identifier as usize];
            let mut found = false;
            let mut offset = 0u64;
            for j in 0..target.reads.len() {
                let packed = target.reads[j];
                let read = (packed >> 33) as usize;
                let r_beg = offset;
                let r_end = offset + u64::from(read_graph.sequences[read].length);
                offset += u64::from(packed as u32);

                let ts = u64::from(transition.target_start);
                let te = u64::from(transition.target_end);
                let ovlp = if r_end.min(te) > r_beg.max(ts) {
                    r_end.min(te) - r_beg.max(ts)
                } else {
                    0
                };
                if found && ovlp == 0 {
                    break;
                }
                if ovlp == 0 {
                    continue;
                }
                found = true;
                if (ovlp as f64) <= (u64::from(read_graph.sequences[read].length) as f64) * 0.8
                {
                    continue;
                }
                if trio_flags[read] & FATHER != 0 {
                    father += 1;
                }
                if trio_flags[read] & MOTHER != 0 {
                    mother += 1;
                }
            }
        }

        let flag = if father >= mother { MOTHER } else { FATHER };
        for k in 0..unitig.reads.len() {
            let read = (unitig.reads[k] >> 33) as usize;
            if trio_flags[read] & SET_TRIO != 0 {
                continue;
            }
            if is_r_het[read] == N_HET {
                continue;
            }
            trio_flags[read] |= flag;
        }
    }

    // Pass 3: normalise every read of an indexed unitig to its trio label.
    for i in 0..index.len() {
        let unitig = &unitigs[i];
        for k in 0..unitig.reads.len() {
            let read = (unitig.reads[k] >> 33) as usize;
            let flag = trio_flags[read];
            if flag & FATHER != 0 {
                trio_flags[read] = FATHER;
            } else if flag & MOTHER != 0 {
                trio_flags[read] = MOTHER;
            } else {
                trio_flags[read] = AMBIGU;
            }
        }
    }
}

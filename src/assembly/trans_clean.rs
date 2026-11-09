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
use crate::assembly::hap_cov::{get_offset_adjust, N_HET, UTrans, RC_0, RC_1, RC_2, RC_3};
use crate::assembly::purge_dups::SET_TRIO;
use crate::assembly::unitig::{UnitigGraph, AMBIGU, FATHER, MOTHER};
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

// ---------------------------------------------------------------------------
// `clean_u_trans_t_idx_filter_adv` (`Overlaps.cpp` 15499)
// ---------------------------------------------------------------------------

/// `trans_ovlp_connect1`: whether the transition agrees with a direct unitig
/// arc within the reference's six-percent tolerance.
fn trans_ovlp_connect1(p: &UTrans, unitig_graph: &UnitigGraph) -> u32 {
    let mut v = p.query_identifier << 1;
    let mut w = (p.target_identifier << 1) | u32::from(p.reverse);

    let mut found = false;
    let mut dg = u32::MAX;
    for _pass in 0..2 {
        let (off, count) = unitig_graph.graph.arc_range(v);
        found = false;
        for i in 0..count {
            let arc = unitig_graph.graph.arcs[off + i];
            if arc.deleted || arc.v != w {
                continue;
            }
            dg = arc.overlap;
            found = true;
            break;
        }
        if found {
            let mut da = p.query_end - p.query_start;
            let mut dif = if dg > da { dg - da } else { da - dg };
            let mut mm = dg.max(da);
            mm = (f64::from(mm) * 0.06) as u32;
            if mm < 8 {
                mm = 8;
            }
            if dif <= mm {
                return 0;
            }

            da = p.target_end - p.target_start;
            dif = if dg > da { dg - da } else { da - dg };
            mm = dg.max(da);
            mm = (f64::from(mm) * 0.06) as u32;
            if mm < 8 {
                mm = 8;
            }
            if dif <= mm {
                return 0;
            }
        }
        v ^= 1;
        w ^= 1;
    }

    1
}

/// `ovlp_rocc`: whether the transition interval is covered by enough whole
/// reads on both unitigs.
fn ovlp_rocc(
    p: &UTrans,
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
    cov_rate: f64,
    cov_cutoff: u32,
) -> u32 {
    for (unitig_id, start, end) in [
        (p.query_identifier, p.query_start, p.query_end),
        (p.target_identifier, p.target_start, p.target_end),
    ] {
        let unitig = &unitig_graph.unitigs[unitig_id as usize];
        if f64::from(end - start) >= f64::from(unitig.length) * cov_rate {
            return 1;
        }
        let mut occ = 0u32;
        let mut offset = 0u32;
        for &packed in &unitig.reads {
            if occ >= cov_cutoff {
                break;
            }
            let read_start = offset;
            let read_end =
                offset + read_graph.sequences[(packed >> 33) as usize].length;
            if read_start >= end {
                break;
            }
            if start <= read_start && end >= read_end {
                occ += 1;
            }
            offset += packed as u32;
        }
        if occ >= cov_cutoff {
            return 1;
        }
    }
    0
}

/// `infer_utrans_ovlp_len`.
fn infer_utrans_ovlp_len(li: &UTrans, lj: &UTrans, unitig_graph: &UnitigGraph) -> i64 {
    let i_n = i64::from(unitig_graph.unitigs[li.target_identifier as usize].length);
    let mut i_s = i64::from(li.target_start);
    let mut i_e = i64::from(li.target_end);
    let i_rev = li.reverse;
    let i_qs = i64::from(li.query_start);
    let i_qe = i64::from(li.query_end);

    let j_n = i64::from(unitig_graph.unitigs[lj.target_identifier as usize].length);
    let mut j_s = i64::from(lj.target_start);
    let mut j_e = i64::from(lj.target_end);
    let j_rev = lj.reverse;
    let j_qs = i64::from(lj.query_start);
    let j_qe = i64::from(lj.query_end);

    let max_s = i_qs.max(j_qs);
    let min_e = i_qe.min(j_qe);
    if min_e <= max_s {
        return 0;
    }
    let mut s_shift = i64::from(get_offset_adjust(
        (max_s - i_qs) as u32,
        (i_qe - i_qs) as u32,
        (i_e - i_s) as u32,
    ));
    let mut e_shift = i64::from(get_offset_adjust(
        (i_qe - min_e) as u32,
        (i_qe - i_qs) as u32,
        (i_e - i_s) as u32,
    ));
    if i_rev {
        std::mem::swap(&mut s_shift, &mut e_shift);
    }
    i_s += s_shift;
    i_e -= e_shift;

    s_shift = i64::from(get_offset_adjust(
        (max_s - j_qs) as u32,
        (j_qe - j_qs) as u32,
        (j_e - j_s) as u32,
    ));
    e_shift = i64::from(get_offset_adjust(
        (j_qe - min_e) as u32,
        (j_qe - j_qs) as u32,
        (j_e - j_s) as u32,
    ));
    if j_rev {
        std::mem::swap(&mut s_shift, &mut e_shift);
    }
    j_s += s_shift;
    j_e -= e_shift;

    if i_rev {
        let ts = i_n - i_e;
        let te = i_n - i_s;
        i_s = ts;
        i_e = te;
    }
    if j_rev {
        let ts = j_n - j_e;
        let te = j_n - j_s;
        j_s = ts;
        j_e = te;
    }

    if i_s <= j_s {
        j_s -= i_s;
        i_s = 0;
    } else {
        i_s -= j_s;
        j_s = 0;
    }

    let i_r = i_n - i_e;
    let j_r = j_n - j_e;
    if i_r <= j_r {
        i_e = i_n;
        j_e += i_r;
    } else {
        j_e = j_n;
        i_e += j_r;
    }
    (i_e - i_s).max(j_e - j_s)
}

/// `is_connect_arc`: whether `li` (the suffix) is explained by the direct
/// unitig arc to `lj`.
fn is_connect_arc(li: &UTrans, lj: &UTrans, unitig_graph: &UnitigGraph, diff_ec_ul: f64) -> u32 {
    if li.query_start >= lj.query_start && li.query_end >= lj.query_end {
        let dq = infer_utrans_ovlp_len(li, lj, unitig_graph);
        let v = (li.target_identifier << 1) | u32::from(li.reverse);
        let w = (lj.target_identifier << 1) | u32::from(lj.reverse);
        let mut dt = -1i64;
        let (off, count) = unitig_graph.graph.arc_range(v);
        for i in 0..count {
            let arc = unitig_graph.graph.arcs[off + i];
            if arc.deleted || arc.v != w {
                continue;
            }
            dt = i64::from(arc.overlap);
            break;
        }
        if dt < 0 {
            return 0;
        }
        let dif = (dq - dt).abs();
        let mm = dq.max(dt);
        let mut mm = (mm as f64 * diff_ec_ul) as i64;
        if mm < 8 {
            mm = 8;
        }
        if dif <= mm {
            return 1;
        }
    }
    0
}

/// `test_arc_rm`: whether the transition `a[a_k]` can be dropped because the
/// remaining transitions of the same query already cover it.
fn test_arc_rm(
    a: &[UTrans],
    a_k: usize,
    length: u32,
    srt: &mut Vec<u64>,
    len_cut: u32,
    rate_cut: f64,
) -> bool {
    srt.clear();
    for (z, entry) in a.iter().enumerate() {
        if z == a_k {
            continue;
        }
        if entry.occurrence == u32::MAX {
            continue;
        }
        srt.push(((u64::from(entry.query_start)) << 32) | u64::from(entry.query_end));
    }
    radix_sort_by_key(srt, &|value: &u64| *value);
    let mut z = 0usize;
    for i in 0..srt.len() {
        let os = (srt[i] >> 32) as u32;
        let oe = srt[i] as u32;
        if z > 0 && os <= srt[z - 1] as u32 {
            if oe > srt[z - 1] as u32 {
                let hi = srt[z - 1] >> 32;
                srt[z - 1] = (hi << 32) | u64::from(oe);
            }
        } else {
            srt[z] = srt[i];
            z += 1;
        }
    }
    srt.truncate(z);

    let mut l = a[a_k].query_end - a[a_k].query_start;
    for &interval in srt.iter() {
        if l == 0 {
            break;
        }
        let os = a[a_k].query_start.max((interval >> 32) as u32);
        let oe = a[a_k].query_end.min(interval as u32);
        let ovq = if oe > os { oe - os } else { 0 };
        debug_assert!(l >= ovq);
        l -= ovq;
    }
    if f64::from(l) > f64::from(len_cut) && f64::from(l) > f64::from(length) * rate_cut {
        return false;
    }
    true
}

/// `deep_clean_u_trans`: retracts exactly one direction of an already-cut
/// pair when the rest of its query still covers the target.
fn deep_clean_u_trans(
    transitions: &mut [UTrans],
    index: &[u64],
    mz_idx: usize,
    unitig_graph: &UnitigGraph,
    srt: &mut Vec<u64>,
) {
    let mz = transitions[mz_idx];
    if mz.query_identifier > mz.target_identifier {
        return;
    }
    let (ro, rn) = u_trans_span(index, mz.target_identifier);
    let mut r_k = None;
    for i in 0..rn {
        if transitions[ro + i].target_identifier == mz.query_identifier {
            r_k = Some(i);
            break;
        }
    }
    let Some(r_k) = r_k else {
        return;
    };
    let cz = transitions[ro + r_k];
    if mz.occurrence == u32::MAX && cz.occurrence == u32::MAX {
        return;
    }
    if mz.occurrence != u32::MAX && cz.occurrence != u32::MAX {
        return;
    }

    let target_len = unitig_graph.unitigs[mz.target_identifier as usize].length;
    let can_drop_target = test_arc_rm(
        &transitions[ro..ro + rn],
        r_k,
        target_len,
        srt,
        3_000_000,
        0.5,
    );
    if !can_drop_target {
        transitions[mz_idx].occurrence = 0;
        transitions[ro + r_k].occurrence = 0;
        return;
    }
    let (qo, qn) = u_trans_span(index, mz.query_identifier);
    let mut k = None;
    for i in 0..qn {
        if transitions[qo + i].target_identifier == mz.target_identifier {
            k = Some(i);
            break;
        }
    }
    let Some(k) = k else {
        return;
    };
    let query_len = unitig_graph.unitigs[mz.query_identifier as usize].length;
    let can_drop_query =
        test_arc_rm(&transitions[qo..qo + qn], k, query_len, srt, 3_000_000, 0.5);
    if !can_drop_query {
        transitions[mz_idx].occurrence = 0;
        transitions[ro + r_k].occurrence = 0;
        return;
    }
    transitions[mz_idx].occurrence = u32::MAX;
    transitions[ro + r_k].occurrence = u32::MAX;
}

/// `worker_for_trans_clean_re`: one representative transition per unitig
/// pair, validated against the graph.
fn worker_for_trans_clean_re(
    transitions: &[UTrans],
    index: &[u64],
    id: u32,
    res: &mut Vec<UTrans>,
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
) {
    let (off, n) = u_trans_span(index, id);
    let a = &transitions[off..off + n];
    let mut l = 0usize;
    for k in 1..=n {
        if k == n || a[l].target_identifier != a[k].target_identifier {
            if k > l {
                // `get_u_trans_spec(ta, a[l].tn, a[l].qn, &r_a, &r_n)`: only
                // the partner's entries that point back at `id`.
                let (ro, rn) = u_trans_span(index, a[l].target_identifier);
                let mut r_beg = ro;
                let mut r_n = 0usize;
                let mut st = 0usize;
                for i in 1..=rn {
                    if i == rn
                        || transitions[ro + i].target_identifier
                            != transitions[ro + st].target_identifier
                    {
                        if transitions[ro + st].target_identifier == a[l].query_identifier {
                            r_beg = ro + st;
                            r_n = i - st;
                            break;
                        }
                        st = i;
                    }
                }
                let r_a = &transitions[r_beg..r_beg + r_n];
                if r_n > 0 && id > a[l].target_identifier {
                    l = k;
                    continue;
                }

                let mut mz: Option<UTrans> = None;
                let mut sz: Option<UTrans> = None;
                let mut consider = |cz: &UTrans, reverse_side: bool, mz: &mut Option<UTrans>| {
                    if cz.flag == RC_3 && sz.map_or(true, |s| s.network_weight < cz.network_weight) {
                        sz = Some(*cz);
                    }
                    if let Some(best) = *mz {
                        if best.flag == RC_0 && cz.flag != RC_0 {
                            return;
                        }
                        if best.flag == RC_1 && (cz.flag == RC_2 || cz.flag == RC_3) {
                            return;
                        }
                    }
                    let better = mz.is_none()
                        || (cz.flag < mz.unwrap().flag
                            && (cz.flag == RC_0 || cz.flag == RC_1))
                        || mz.unwrap().network_weight < cz.network_weight
                        || (mz.unwrap().network_weight == cz.network_weight
                            && if reverse_side {
                                (mz.unwrap().query_end - mz.unwrap().query_start)
                                    < (cz.target_end - cz.target_start)
                            } else {
                                (mz.unwrap().query_end - mz.unwrap().query_start)
                                    < (cz.query_end - cz.query_start)
                            });
                    if better {
                        *mz = Some(*cz);
                    }
                };
                for cz in &a[l..k] {
                    consider(cz, false, &mut mz);
                }
                for cz in r_a {
                    consider(cz, true, &mut mz);
                }

                if let Some(mut mz) = mz {
                    if mz.query_identifier != mz.target_identifier {
                        if let Some(sz) = sz {
                            if sz.reverse == mz.reverse {
                                if mz.query_identifier != sz.query_identifier
                                    || mz.target_identifier != sz.target_identifier
                                {
                                    // Re-orient the copy of `mz` to the other
                                    // direction before comparing with `sz`.
                                    let qn = mz.query_identifier;
                                    let qs = mz.query_start;
                                    let qe = mz.query_end;
                                    mz.query_identifier = mz.target_identifier;
                                    mz.query_start = mz.target_start;
                                    mz.query_end = mz.target_end;
                                    mz.target_identifier = qn;
                                    mz.target_start = qs;
                                    mz.target_end = qe;
                                }
                                if mz.query_identifier == sz.query_identifier
                                    && mz.target_identifier == sz.target_identifier
                                {
                                    let os = mz.query_start.max(sz.query_start);
                                    let oe = mz.query_end.min(sz.query_end);
                                    let ovq = if oe > os { oe - os } else { 0 };
                                    let os = mz.target_start.max(sz.target_start);
                                    let oe = mz.target_end.min(sz.target_end);
                                    let ovt = if oe > os { oe - os } else { 0 };
                                    if ovq > 0
                                        && f64::from(ovq)
                                            > f64::from(mz.query_end - mz.query_start) * 0.8
                                        && f64::from(ovq)
                                            > f64::from(sz.query_end - sz.query_start) * 0.8
                                        && ovt > 0
                                        && f64::from(ovt)
                                            > f64::from(mz.target_end - mz.target_start) * 0.8
                                        && f64::from(ovt)
                                            > f64::from(sz.target_end - sz.target_start) * 0.8
                                    {
                                        let flag = mz.flag;
                                        mz = sz;
                                        mz.flag = flag;
                                    }
                                }
                            }
                        }
                        if mz.query_identifier != id {
                            let qn = mz.query_identifier;
                            let tn = mz.target_identifier;
                            let qs = mz.query_start;
                            let qe = mz.query_end;
                            let ts = mz.target_start;
                            let te = mz.target_end;
                            mz.query_identifier = tn;
                            mz.query_start = ts;
                            mz.query_end = te;
                            mz.target_identifier = qn;
                            mz.target_start = qs;
                            mz.target_end = qe;
                        }

                        res.push(mz);
                        if mz.network_weight >= 0.0
                            && (trans_ovlp_connect1(&mz, unitig_graph) == 0
                                || ovlp_rocc(&mz, unitig_graph, read_graph, 0.8, 3) == 0)
                        {
                            res.pop();
                        }
                    }
                }
            }
            l = k;
        }
    }
}

/// `worker_for_trans_sec_cut`, first phase: removes the transitions whose
/// query interval is already covered by stronger ones.
fn trans_sec_cut_pass(
    transitions: &mut [UTrans],
    index: &[u64],
    id: u32,
    unitig_graph: &UnitigGraph,
    sec_rate: f64,
    sc_sec_rate: f64,
    srt: &mut Vec<u64>,
) {
    let (off, n) = u_trans_span(index, id);
    if n == 0 {
        return;
    }
    let mut a: Vec<UTrans> = transitions[off..off + n].to_vec();
    let mut r_n = 0usize;
    for k in 0..n {
        if a[k].flag == RC_0 || a[k].flag == RC_1 {
            a.swap(k, r_n);
            r_n += 1;
        }
    }
    if r_n >= n {
        transitions[off..off + n].copy_from_slice(&a);
        return;
    }
    if n - r_n > 1 {
        a[r_n..].sort_by(|x, y| {
            y.network_weight
                .partial_cmp(&x.network_weight)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    srt.clear();
    for k in 0..n {
        srt.push(((u64::from(a[k].query_start)) << 32) | k as u64);
    }
    radix_sort_by_key(srt, &|value: &u64| *value);
    let mut ss = vec![0usize; n];
    let mut hs = vec![0usize; n];
    for k in 0..n {
        ss[k] = srt[k] as u32 as usize;
        hs[ss[k]] = k;
    }

    let mut bu = vec![0u64; n];
    let mut wu = vec![0u64; n];
    for k in r_n..n {
        let cz = a[k];
        let len = (f64::from(cz.query_end - cz.query_start) * sec_rate) as u32;
        let mut m = 0usize;
        let mut wm = 0usize;
        let mut l = 0u32;
        let mut wl = 0u32;
        for z in 0..n {
            if z == hs[k] {
                continue;
            }
            if ss[z] > k {
                continue;
            }
            let mz = a[ss[z]];
            if mz.query_start >= cz.query_end {
                break;
            }
            if mz.occurrence == u32::MAX {
                continue;
            }
            if mz.flag != RC_0
                && mz.flag != RC_1
                && cz.network_weight > mz.network_weight * sc_sec_rate
            {
                continue;
            }
            let os = cz.query_start.max(mz.query_start);
            let oe = cz.query_end.min(mz.query_end);
            let ovq = if oe > os { oe - os } else { 0 };
            if ovq == 0 {
                continue;
            }
            if is_connect_arc(&cz, &mz, unitig_graph, 0.04) != 0
                || is_connect_arc(&mz, &cz, unitig_graph, 0.04) != 0
            {
                continue;
            }

            if m > 0 && (bu[m - 1] as u32) >= os {
                if oe > bu[m - 1] as u32 {
                    l += oe - (bu[m - 1] as u32);
                    bu[m - 1] += u64::from(oe - (bu[m - 1] as u32));
                }
            } else {
                l += oe - os;
                bu[m] = (u64::from(os) << 32) | u64::from(oe);
                m += 1;
            }

            if wm > 0 && (wu[wm - 1] as u32) >= mz.query_start {
                if mz.query_end > wu[wm - 1] as u32 {
                    wl += mz.query_end - (wu[wm - 1] as u32);
                    wu[wm - 1] += u64::from(mz.query_end - (wu[wm - 1] as u32));
                }
            } else {
                wl += mz.query_end - mz.query_start;
                wu[wm] = (u64::from(mz.query_start) << 32) | u64::from(mz.query_end);
                wm += 1;
            }

            if l > 0 && (l >= len || f64::from(l) >= f64::from(wl) * sec_rate) {
                a[k].occurrence = u32::MAX;
                break;
            }
        }
    }
    transitions[off..off + n].copy_from_slice(&a);
}

/// `gen_u_trans_cluster`: the connected components of the transition graph
/// and their spans.
fn gen_u_trans_cluster(transitions: &[UTrans], index: &[u64]) -> (Vec<u64>, Vec<u64>) {
    let n = index.len();
    let mut z = vec![u64::MAX; n];
    let mut stack: Vec<u32> = Vec::new();
    let mut cn = 0usize;
    for i in 0..n as u32 {
        if z[i as usize] != u64::MAX {
            continue;
        }
        stack.clear();
        stack.push(i);
        let mut m = 0usize;
        while let Some(v) = stack.pop() {
            if z[v as usize] != u64::MAX {
                continue;
            }
            z[v as usize] = ((i as u64) << 32) | u64::from(v);
            m += 1;
            let (off, count) = u_trans_span(index, v);
            for k in 0..count {
                let t = transitions[off + k].target_identifier;
                if z[t as usize] == u64::MAX {
                    stack.push(t);
                }
            }
        }
        if m > 1 {
            cn += 1;
        }
    }
    radix_sort_by_key(&mut z, &|value: &u64| *value);
    let mut idx = Vec::new();
    let mut l = 0usize;
    for k in 1..=n {
        if k == n || (z[k] >> 32) != (z[l] >> 32) {
            if k - l > 1 {
                idx.push(((l as u64) << 32) | k as u64);
            }
            l = k;
        }
    }
    debug_assert_eq!(idx.len(), cn);
    (z, idx)
}

/// `clean_u_trans_t_idx_filter_adv`: the transition filter used by the
/// trio-dual stage.
///
/// The reference leaves the table buffer longer than the logical length (the
/// compaction only moves the survivors forward) and keeps the stale index
/// entries of the queries that lost every transition; later span lookups read
/// those leftovers. The port returns the full buffer and the logical length
/// so the same reads are reproduced.
pub fn clean_u_trans_t_idx_filter_adv(
    transitions: &mut Vec<UTrans>,
    index: &mut Vec<u64>,
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
    sc_sec_rate: f64,
    uniform_only: bool,
) -> usize {
    let n_seq = unitig_graph.graph.sequences.len() as u32;
    kt_u_trans_t_idx(transitions, index, n_seq);

    let mut res: Vec<UTrans> = Vec::new();
    for i in 0..n_seq {
        worker_for_trans_clean_re(transitions, index, i, &mut res, unitig_graph, read_graph);
    }

    let mut rebuilt: Vec<UTrans> = Vec::with_capacity(res.len() * 2);
    for entry in &res {
        if entry.query_identifier == entry.target_identifier {
            continue;
        }
        let mut forward = *entry;
        forward.occurrence = 0;
        if forward.flag == RC_3 {
            forward.flag = RC_2;
        }
        rebuilt.push(forward);
        let mut reverse = forward;
        reverse.query_identifier = forward.target_identifier;
        reverse.query_start = forward.target_start;
        reverse.query_end = forward.target_end;
        reverse.target_identifier = forward.query_identifier;
        reverse.target_start = forward.query_start;
        reverse.target_end = forward.query_end;
        rebuilt.push(reverse);
    }
    *transitions = rebuilt;
    kt_u_trans_t_idx(transitions, index, n_seq);
    if uniform_only {
        return transitions.len();
    }

    let mut srt: Vec<u64> = Vec::new();
    for i in 0..n_seq {
        trans_sec_cut_pass(
            transitions,
            index,
            i,
            unitig_graph,
            0.5,
            sc_sec_rate,
            &mut srt,
        );
    }

    let (cluster_z, cluster_idx) = gen_u_trans_cluster(transitions, index);
    for id in 0..cluster_idx.len() as u32 {
        let span = cluster_idx[id as usize];
        let start = (span >> 32) as usize;
        let end = span as u32 as usize;
        for z in start..end {
            let node = cluster_z[z] as u32;
            let (off, count) = u_trans_span(index, node);
            for k in 0..count {
                deep_clean_u_trans(transitions, index, off + k, unitig_graph, &mut srt);
            }
        }
    }

    let mut st = 0usize;
    for k in 0..transitions.len() {
        if transitions[k].occurrence == u32::MAX {
            continue;
        }
        transitions[st] = transitions[k];
        st += 1;
    }
    let logical_n = st;

    let mut st = 0usize;
    for i in 1..=logical_n {
        if i == logical_n
            || transitions[i].query_identifier != transitions[st].query_identifier
        {
            index[transitions[st].query_identifier as usize] =
                ((st as u64) << 32) | ((i - st) as u64);
            st = i;
        }
    }
    logical_n
}

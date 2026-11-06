//! Window rechaining, the counterpart of `gen_win_chain` and
//! `rechain_aln_hc` in `Correct.cpp`.
//!
//! After the block-wise pass accepts an overlap, the blocks it could not align
//! become unaligned windows. Each of them is retried here: the exact-match runs
//! of the neighbouring windows are extracted from their stored CIGARs, chained
//! with the fixed-end dynamic programme, turned into seed anchors, and passed
//! back to [`hc_ovlp_base_direct`]. A successful retry replaces the unaligned
//! window with the newly produced aligned one.
//!
//! The anchor buffer is the reference `Candidates_list`: `length` is a fill
//! pointer that the caller rewinds after the retry, so the anchors produced
//! here are scratch data and never become part of the read's permanent chain.

use crate::container::ksort::radix_sort_by_key_32;
use crate::correction::cal_exz::{
    FORCED_SINGLE_LENGTH, MAXIMUM_SINGLE_LENGTH, MODE_BACKWARD_EXTENSION, MODE_FORWARD_EXTENSION,
    MODE_GLOBAL,
};
use crate::correction::hc_align::{backtrace_window_in_place, BitExtz};
use crate::correction::ovlp_base::hc_ovlp_base_direct;
use crate::correction::window::{window_id_by_start, WindowCigar, WindowList};
use crate::overlap::chaining::{chain_dp_fixed_ends, chain_simple};
use crate::overlap::lchain_parameters::set_lchain_parameters;
use crate::overlap::region::{CandidatesList, KmerHit, OverlapRegion};
use crate::reads::read_store::ReadStore;

/// `is_ualn_win`: an unaligned window has no CIGAR and no target end.
pub fn is_unaligned_window(window: &WindowList) -> bool {
    window.error == i16::MAX && window.cigar_length == 0 && window.extra_end < 0
}

/// `push_khit`: appends one run of `length` matching bases ending at
/// `(query_end, target_end)`.
///
/// The count packed into `cnt` is `coverage << 8 | span`; a run longer than
/// 255 bases is split into chunks of 255. An `initial_count` overrides the
/// coverage tag (the sentinels use `1`, the fallback match uses `3`). A zero
/// length still emits one entry, which is how the sentinels are recorded.
pub fn push_kmer_hit(
    candidates: &mut CandidatesList,
    query_end: i32,
    target_end: i32,
    length: u32,
    h_kmer_hit: u32,
    initial_count: Option<u32>,
) {
    let mut code = if length >= h_kmer_hit { 1 } else { 2 };
    if let Some(value) = initial_count {
        code = value;
    }
    code <<= 8;

    let mut remaining = length;
    if remaining > 0 {
        while remaining >= 0xff {
            candidates.push(KmerHit {
                read_identifier: 0,
                reverse_strand: false,
                offset: target_end as u32,
                self_offset: query_end as u32,
                count: code.wrapping_add(0xff),
            });
            remaining -= 0xff;
        }
        if remaining > 0 {
            candidates.push(KmerHit {
                read_identifier: 0,
                reverse_strand: false,
                offset: target_end as u32,
                self_offset: query_end as u32,
                count: code.wrapping_add(remaining),
            });
        }
    } else {
        candidates.push(KmerHit {
            read_identifier: 0,
            reverse_strand: false,
            offset: target_end as u32,
            self_offset: query_end as u32,
            count: code,
        });
    }
}

/// `extract_exact_cigar`: turns the match runs of one window CIGAR into seed
/// anchors clipped to the retried interval.
///
/// `pattern_start`/`text_start` are the target/query coordinates the CIGAR was
/// generated at; `pattern_min..pattern_max` and `text_min..text_max` are the
/// half-open retried intervals. A run is usable only when it is fully inside
/// both intervals. When no run reaches `minimum_run` but one is shorter, the
/// longest one is emitted with the special coverage tag `3` if it ends at
/// least `minimum_window_length` bases after the window start, which marks the
/// whole window as an exact block.
#[allow(clippy::too_many_arguments)]
pub fn extract_exact_cigar(
    cells: &[u16],
    pattern_start: i32,
    text_start: i32,
    pattern_min: i32,
    pattern_max: i32,
    text_min: i32,
    text_max: i32,
    candidates: &mut CandidatesList,
    minimum_run: i32,
    minimum_window_length: i64,
    h_kmer_hit: i64,
) -> u32 {
    let mut cell_index = 0usize;
    let mut pattern_index = pattern_start;
    let mut text_index = text_start;
    let mut pattern_offset = -1i32;
    let mut text_offset = -1i32;
    let mut maximum_length = -1i32;
    let mut occurrence = 0u32;

    while cell_index < cells.len() && pattern_index < pattern_max && text_index < text_max {
        // `pop_trace`: consecutive cells with the same operation form one run.
        let operation = (cells[cell_index] >> 14) as u8;
        let mut run_length = u32::from(cells[cell_index] & 0x3fff);
        cell_index += 1;
        while cell_index < cells.len() && (cells[cell_index] >> 14) as u8 == operation {
            run_length += u32::from(cells[cell_index] & 0x3fff);
            cell_index += 1;
        }
        let run_length = run_length as i32;

        match operation {
            0 => {
                let pattern_run_end = pattern_index.wrapping_add(run_length);
                let text_run_end = text_index.wrapping_add(run_length);
                let position = pattern_index.max(pattern_min);
                let pattern_overlap_end = pattern_run_end.min(pattern_max);
                let text_overlap_start = text_index.max(text_min);
                let text_overlap_end = text_run_end.min(text_max);
                if pattern_overlap_end > position && text_overlap_end > text_overlap_start {
                    let length = pattern_overlap_end - position;
                    if length == text_overlap_end - text_overlap_start {
                        let pattern_last = pattern_overlap_end - 1;
                        let text_last = text_overlap_end - 1;
                        if length > maximum_length {
                            pattern_offset = pattern_last;
                            text_offset = text_last;
                            maximum_length = length;
                        }
                        if length >= minimum_run {
                            push_kmer_hit(
                                candidates,
                                text_last,
                                pattern_last,
                                length as u32,
                                h_kmer_hit as u32,
                                None,
                            );
                            occurrence += 1;
                        }
                    }
                }
                pattern_index = pattern_run_end;
                text_index = text_run_end;
            }
            1 => {
                pattern_index = pattern_index.wrapping_add(run_length);
                text_index = text_index.wrapping_add(run_length);
            }
            // More target than query: only the target cursor advances.
            2 => pattern_index = pattern_index.wrapping_add(run_length),
            _ => text_index = text_index.wrapping_add(run_length),
        }
    }

    // A short exact run that covers the whole retried window is still usable.
    if maximum_length > 0
        && maximum_length < minimum_run
        && text_start >= text_min
        && i64::from(text_index) >= minimum_window_length + i64::from(text_start)
    {
        push_kmer_hit(
            candidates,
            text_offset,
            pattern_offset,
            maximum_length as u32,
            h_kmer_hit as u32,
            Some(3),
        );
        occurrence += 1;
    }

    occurrence
}

/// `gen_single_khit`: turns the chained exact-match runs into the seed list
/// consumed by [`hc_ovlp_base_direct`].
///
/// Every run becomes an end anchor and, when its span exceeds one base, a
/// start anchor as well; the sentinels mark the fixed ends of the window.
/// Exact anchors are deduplicated against the sentinels afterwards, and the
/// remaining list is chained again with the simple programme because the
/// run-level anchors already form a nearly colinear chain.
#[allow(clippy::too_many_arguments)]
pub fn generate_single_kmer_hit_chain(
    candidates: &mut CandidatesList,
    hit_count: i64,
    h_kmer_hit: i64,
    mode: i64,
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    maximum_skip: i64,
    maximum_iterations: i64,
) -> i64 {
    let base = candidates.length as usize;
    let prefix = i64::from(mode == MODE_GLOBAL || mode == MODE_FORWARD_EXTENSION);
    let suffix = i64::from(mode == MODE_GLOBAL || mode == MODE_BACKWARD_EXTENSION);

    // Drop the sentinels inserted by `push_kmer_hit`; every real anchor also
    // contributes a start anchor when its span exceeds one base.
    let mut occurrence = 0i64;
    let mut kept = 0usize;
    for index in 0..hit_count as usize {
        let hit = candidates.list[base + index];
        if hit.count & 0xff == 0 {
            continue;
        }
        occurrence += 1;
        if hit.count & 0xff > 1 {
            occurrence += 1;
        }
        candidates.list[base + kept] = hit;
        kept += 1;
    }
    let mut hit_count = kept as i64;
    if hit_count == 0 {
        return 0;
    }

    occurrence += prefix + suffix;
    candidates.reserve_at_length(occurrence);

    // Fill the array from the back so the anchors come out in ascending order.
    let mut write = occurrence;
    if suffix != 0 {
        write -= 1;
        candidates.list[base + write as usize] = KmerHit {
            read_identifier: 1,
            reverse_strand: false,
            offset: target_end as u32,
            self_offset: query_end as u32,
            count: 1,
        };
    }

    for index in (0..hit_count).rev() {
        let hit = candidates.list[base + index as usize];
        let span = hit.count & 0xff;
        if span == 0 {
            continue;
        }
        let mut read_identifier = span;
        if i64::from(span) < h_kmer_hit {
            read_identifier = span + 1;
        }

        write -= 1;
        candidates.list[base + write as usize] = KmerHit {
            read_identifier,
            reverse_strand: false,
            offset: hit.offset.wrapping_add(1),
            self_offset: hit.self_offset.wrapping_add(1),
            count: span,
        };

        if span > 1 {
            write -= 1;
            candidates.list[base + write as usize] = KmerHit {
                read_identifier,
                reverse_strand: false,
                offset: hit.offset.wrapping_add(1).wrapping_sub(span),
                self_offset: hit.self_offset.wrapping_add(1).wrapping_sub(span),
                count: span,
            };
        }
    }

    if prefix != 0 {
        write -= 1;
        candidates.list[base + write as usize] = KmerHit {
            read_identifier: 1,
            reverse_strand: false,
            offset: target_start as u32,
            self_offset: query_start as u32,
            count: 1,
        };
    }
    debug_assert_eq!(write, 0);
    hit_count = occurrence;

    // Deduplicate anchors sitting exactly on a sentinel. The first comparison
    // is an assignment in the reference — it overwrites the target coordinate
    // before dropping the anchor — and is reproduced as such.
    let mut query_boundary = [u64::MAX; 2];
    let mut target_boundary = [u64::MAX; 2];
    if prefix != 0 {
        query_boundary[0] = query_start as u64;
        target_boundary[0] = target_start as u64;
    }
    if suffix != 0 {
        query_boundary[1] = query_end as u64;
        target_boundary[1] = target_end as u64;
    }

    let mut sorted = true;
    let mut colinear = 0i64;
    let mut kept = 0usize;
    for index in 0..hit_count as usize {
        // The reference tests `(ch_a[k].offset = t[0])`, an assignment rather
        // than a comparison: the target coordinate is overwritten with the
        // prefix target before the entry may be dropped, and the condition is
        // true only when that stored value is not zero.
        if index > 0
            && u64::from(candidates.list[base + index].self_offset) == query_boundary[0]
        {
            let stored_target = target_boundary[0] as u32;
            candidates.list[base + index].offset = stored_target;
            if stored_target != 0 {
                continue;
            }
        }
        if index + 1 < hit_count as usize
            && u64::from(candidates.list[base + index].self_offset) == query_boundary[1]
            && u64::from(candidates.list[base + index].offset) == target_boundary[1]
        {
            continue;
        }

        if kept > 0 {
            let current = candidates.list[base + index];
            let previous = candidates.list[base + kept - 1];
            if current.self_offset > previous.self_offset && current.offset > previous.offset {
                colinear += 1;
            }
            if current.self_offset <= previous.self_offset {
                sorted = false;
            }
        } else {
            colinear += 1;
        }
        candidates.list[base + kept] = candidates.list[base + index];
        kept += 1;
    }
    hit_count = kept as i64;

    if colinear == hit_count {
        return hit_count;
    }

    if !sorted {
        radix_sort_by_key_32(
            &mut candidates.list[base..base + hit_count as usize],
            &|hit: &KmerHit| hit.self_offset,
        );
        let mut group_start = 0usize;
        let mut index = 1usize;
        while index <= hit_count as usize {
            if index == hit_count as usize
                || candidates.list[base + index].self_offset
                    != candidates.list[base + group_start].self_offset
            {
                if index - group_start > 1 {
                    radix_sort_by_key_32(
                        &mut candidates.list[base + group_start..base + index],
                        &|hit: &KmerHit| hit.offset,
                    );
                }
                group_start = index;
            }
            index += 1;
        }
    }

    let previous_count = hit_count;
    let sub_count = hit_count - prefix - suffix;
    let prefix_index = prefix as usize;
    let chain_count = chain_simple(
        &mut candidates.list
            [base + prefix_index..base + prefix_index + sub_count.max(0) as usize],
        sub_count,
        &mut candidates.chain_data,
        maximum_skip,
        maximum_iterations,
    );
    let final_count = chain_count + prefix + suffix;
    if suffix != 0 {
        candidates.list[base + (final_count - 1) as usize] =
            candidates.list[base + (previous_count - 1) as usize];
    }

    final_count
}

/// `gen_win_chain`: rebuilds one window of an overlap from the exact matches
/// its neighbouring windows already contain.
///
/// Only the `rref` variant (target reads in memory) is implemented, which is
/// the one `rechain_aln_hc` uses.
#[allow(clippy::too_many_arguments)]
pub fn generate_window_chain(
    region: &mut OverlapRegion,
    candidates: &mut CandidatesList,
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    window_length: i64,
    reads: &ReadStore,
    query: &[u8],
    extz: &mut BitExtz,
    query_total_length: i64,
    target_total_length: i64,
    error_rate: f64,
    h_kmer_hit: i64,
    mode: i64,
    accurate: bool,
) -> i64 {
    assert!(mode < MODE_SEMI_GLOBAL_LIMIT);
    let initial_length = candidates.length;
    let window_count = region.window_cigar.windows.len() as i64;

    let mut window_start = query_start;
    if window_start < i64::from(region.x_position_start) {
        window_start = i64::from(region.x_position_start);
    }
    let (mut window_index, _) = window_id_by_start(
        region.x_position_start as i32,
        region.x_position_end as i32,
        (window_start / window_length) * window_length,
        window_length,
    );

    // Locate the window that contains the start of the retried interval.
    if window_index >= window_count {
        window_index = window_count - 1;
    }
    while window_index < window_count
        && query_start
            > i64::from(region.window_cigar.windows[window_index as usize].x_end)
    {
        window_index += 1;
    }
    if window_index >= window_count {
        window_index = window_count - 1;
    }
    while window_index >= 0
        && query_start
            < i64::from(region.window_cigar.windows[window_index as usize].x_start)
    {
        window_index -= 1;
    }
    if window_index < 0 {
        window_index = 0;
    }

    // Sentinels marking the fixed ends, as in `gen_single_khit`'s output.
    if mode == MODE_GLOBAL || mode == MODE_FORWARD_EXTENSION {
        push_kmer_hit(
            candidates,
            query_start as i32,
            target_start as i32,
            0,
            0,
            Some(1),
        );
    }

    let mut occurrence = 0i64;
    let mut index = window_index;
    while index < window_count
        && i64::from(region.window_cigar.windows[index as usize].x_start) < query_end
    {
        if region.window_cigar.windows[index as usize].y_end != -1 {
            let window_x_start = i64::from(region.window_cigar.windows[index as usize].x_start);
            let window_x_end = i64::from(region.window_cigar.windows[index as usize].x_end) + 1;
            let overlap_start = query_start.max(window_x_start);
            let overlap_end = query_end.min(window_x_end);
            if overlap_end > overlap_start {
                // A directly aligned window stores no CIGAR yet; regenerate it
                // from its recorded boundaries.
                if region.window_cigar.windows[index as usize].cigar_length == 0 {
                    let reverse_strand = region.y_strand as u8;
                    backtrace_window_in_place(
                        region,
                        index,
                        reads,
                        query,
                        extz,
                        reverse_strand,
                    );
                }

                let window = region.window_cigar.windows[index as usize];
                let cell_start = window.cigar_start as usize;
                let cell_length = window.cigar_length as usize;
                let cells = &region.window_cigar.cells[cell_start..cell_start + cell_length];
                occurrence += i64::from(extract_exact_cigar(
                    cells,
                    window.y_start,
                    window.x_start,
                    target_start as i32,
                    target_end as i32,
                    query_start as i32,
                    query_end as i32,
                    candidates,
                    10,
                    window_length,
                    h_kmer_hit,
                ));
            }
        }
        index += 1;
    }

    if mode == MODE_GLOBAL || mode == MODE_BACKWARD_EXTENSION {
        push_kmer_hit(
            candidates,
            query_end as i32 - 1,
            target_end as i32 - 1,
            0,
            0,
            Some(1),
        );
    }

    if occurrence == 0 {
        candidates.length = initial_length;
        return 0;
    }

    let extended_count = candidates.length;
    candidates.length = initial_length;
    let chain_hit_count = extended_count - initial_length;

    let parameters = set_lchain_parameters(accurate, h_kmer_hit as u32);
    let maximum_distance = MAXIMUM_SINGLE_LENGTH >> 1;
    let left_fixed = mode == MODE_GLOBAL || mode == MODE_FORWARD_EXTENSION;
    let right_fixed = mode == MODE_GLOBAL || mode == MODE_BACKWARD_EXTENSION;

    let base = initial_length as usize;
    let chain_count = chain_dp_fixed_ends(
        &mut candidates.list[base..extended_count as usize],
        chain_hit_count,
        &mut candidates.chain_data,
        parameters.maximum_skip,
        parameters.maximum_iterations,
        maximum_distance,
        parameters.gap_penalty,
        parameters.skip_penalty,
        error_rate,
        query_total_length,
        target_total_length,
        true,
        left_fixed,
        right_fixed,
    );

    let mut occurrence = 0i64;
    for index in 0..chain_count {
        let source = candidates.chain_data.scratch[index as usize] as usize;
        let hit = candidates.list[base + source];
        candidates.list[base + index as usize] = hit;
        if hit.count & 0xff != 0 {
            occurrence += 1;
        }
    }
    if occurrence <= 0 {
        return 0;
    }

    generate_single_kmer_hit_chain(
        candidates,
        chain_count,
        h_kmer_hit,
        mode,
        query_start,
        query_end,
        target_start,
        target_end,
        parameters.maximum_skip,
        parameters.maximum_iterations,
    )
}

/// `gen_win_chain` only knows the three fixed-end modes; the semi-global
/// alignment has its own entry point.
const MODE_SEMI_GLOBAL_LIMIT: i64 = 3;

/// `rechain_aln_hc`: retries one unaligned window of `output`.
///
/// The retry is accepted only when the fresh alignment produced exactly one
/// new window whose fixed boundaries match the old unaligned one; the new
/// window then replaces the old entry in place.
#[allow(clippy::too_many_arguments)]
pub fn rechain_alignment_hc(
    region: &mut OverlapRegion,
    candidates: &mut CandidatesList,
    output: &mut WindowCigar,
    auxiliary_index: usize,
    window_length: i64,
    reads: &ReadStore,
    query: &[u8],
    extz: &mut BitExtz,
    error_rate: f64,
    query_total_length: i64,
    target_total_length: i64,
    h_kmer_hit: i64,
) {
    let initial_length = candidates.length;
    let window = output.windows[auxiliary_index];
    let query_start = i64::from(window.x_start);
    let query_end = i64::from(window.x_end) + 1;
    let target_start = i64::from(window.y_start);
    let target_end = i64::from(window.y_end) + 1;
    if query_end - query_start < FORCED_SINGLE_LENGTH
        || target_end - target_start < FORCED_SINGLE_LENGTH
    {
        return;
    }

    let mode = i64::from(window.error_threshold);
    let chain_count = generate_window_chain(
        region,
        candidates,
        query_start,
        query_end,
        target_start,
        target_end,
        window_length,
        reads,
        query,
        extz,
        query_total_length,
        target_total_length,
        error_rate,
        h_kmer_hit,
        mode,
        true,
    );

    if chain_count != 0 {
        let base = initial_length as usize;
        let mut retry = true;
        let last = base + chain_count as usize - 1;
        match mode {
            MODE_GLOBAL => {
                assert_eq!(i64::from(candidates.list[base].self_offset), query_start);
                assert_eq!(i64::from(candidates.list[base].offset), target_start);
                assert_eq!(i64::from(candidates.list[last].self_offset), query_end);
                assert_eq!(i64::from(candidates.list[last].offset), target_end);
                if chain_count <= 2 {
                    retry = false;
                }
            }
            MODE_FORWARD_EXTENSION => {
                assert_eq!(i64::from(candidates.list[base].self_offset), query_start);
                assert_eq!(i64::from(candidates.list[base].offset), target_start);
                if chain_count <= 1 {
                    retry = false;
                }
            }
            MODE_BACKWARD_EXTENSION => {
                assert_eq!(i64::from(candidates.list[last].self_offset), query_end);
                assert_eq!(i64::from(candidates.list[last].offset), target_end);
                if chain_count <= 1 {
                    retry = false;
                }
            }
            _ => {}
        }

        if retry {
            let previous_count = output.windows.len() as i64;
            hc_ovlp_base_direct(
                region,
                &candidates.list[base..base + chain_count as usize],
                chain_count,
                window_length,
                reads,
                query,
                extz,
                output,
                error_rate,
                query_total_length,
                target_total_length,
                mode,
            );
            let new_count = output.windows.len() as i64;

            // The old unaligned window can be replaced only when the retry
            // produced exactly one aligned window.
            if new_count == previous_count + 1
                && !is_unaligned_window(&output.windows[(new_count - 1) as usize])
            {
                let old_window = output.windows[auxiliary_index];
                let new_window = output.windows[(new_count - 1) as usize];
                let query_start_matches = old_window.x_start == new_window.x_start;
                let query_end_matches = old_window.x_end == new_window.x_end;
                let target_start_matches = old_window.y_start == new_window.y_start;
                let target_end_matches = old_window.y_end == new_window.y_end;
                let replace = match mode {
                    MODE_GLOBAL => {
                        query_start_matches
                            && query_end_matches
                            && target_start_matches
                            && target_end_matches
                    }
                    MODE_FORWARD_EXTENSION => query_start_matches && target_start_matches,
                    MODE_BACKWARD_EXTENSION => query_end_matches && target_end_matches,
                    _ => false,
                };
                if replace {
                    output.windows[auxiliary_index] = new_window;
                    output.windows.pop();
                }
            }
        }
    }

    candidates.length = initial_length;
}

#[cfg(test)]
mod probe_goldens {
    use super::*;
    use crate::correction::hc_align::BitExtz;
    use crate::correction::ovlp_base::hc_ovlp_base_direct;
    use crate::correction::window::WindowCigar;
    use crate::overlap::region::{CandidatesList, ChainData, KmerHit, OverlapRegion};
    use crate::reads::read_store::ReadStore;

    struct Lcg(u64);

    impl Lcg {
        fn next_u32(&mut self) -> u32 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 33) as u32
        }

        fn sequence(&mut self, length: usize) -> Vec<u8> {
            const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
            (0..length)
                .map(|_| BASES[(self.next_u32() & 3) as usize])
                .collect()
        }
    }

    fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                other => *other,
            })
            .collect()
    }

    fn push_hits(lines: &mut Vec<String>, tag: &str, hits: &[KmerHit], count: i64) {
        if count <= 0 {
            return;
        }
        for index in 0..count as usize {
            let hit = hits[index];
            lines.push(format!(
                "{} i={} x={} y={} cnt={} rid={} st={}",
                tag,
                index,
                hit.self_offset,
                hit.offset,
                hit.count,
                hit.read_identifier,
                hit.strand_value()
            ));
        }
    }

    /// `print_windows`: the reference zeroes the extra boundaries of every
    /// aligned window before printing, and that zeroing persists.
    fn push_windows(lines: &mut Vec<String>, tag: &str, cigar: &mut WindowCigar) {
        lines.push(format!(
            "{} n={} cn={}",
            tag,
            cigar.windows.len(),
            cigar.cells.len()
        ));
        for index in 0..cigar.windows.len() {
            if cigar.windows[index].error != i16::MAX {
                cigar.windows[index].extra_begin = 0;
                cigar.windows[index].extra_end = 0;
            }
            let window = cigar.windows[index];
            let mut line = format!(
                " {}W {} x=[{},{}] y=[{},{}] eb={} ee={} err={} thre={} cidx={} clen={} :",
                tag,
                index,
                window.x_start,
                window.x_end,
                window.y_start,
                window.y_end,
                window.extra_begin,
                window.extra_end,
                window.error,
                window.error_threshold,
                window.cigar_start,
                window.cigar_length
            );
            for offset in 0..window.cigar_length {
                let cell = cigar.cells[(window.cigar_start + offset) as usize];
                line.push_str(&format!(
                    " {}{}",
                    cell & 0x3fff,
                    b"EMDI"[(cell >> 14) as usize] as char
                ));
            }
            lines.push(line);
        }
    }

    // ------------------------------------------------------------- chaining

    struct ChainScenario {
        name: &'static str,
        count: usize,
        seed: u64,
        dx_min: i64,
        dx_span: u32,
        dy_min: i64,
        dy_span: u32,
        span_min: u32,
        span_span: u32,
        coverage_min: u32,
        coverage_span: u32,
        quick: bool,
        left: bool,
        right: bool,
        maximum_skip: i64,
        maximum_iterations: i64,
        maximum_distance: i64,
        query_length: i64,
        target_length: i64,
        band_width_rate: f64,
    }

    #[allow(clippy::too_many_arguments)]
    fn generate_chain_hits(
        count: usize,
        seed: u64,
        x0: i64,
        y0: i64,
        dx_min: i64,
        dx_span: u32,
        dy_min: i64,
        dy_span: u32,
        span_min: u32,
        span_span: u32,
        coverage_min: u32,
        coverage_span: u32,
    ) -> Vec<KmerHit> {
        let mut generator = Lcg(seed);
        let mut x = x0;
        let mut y = y0;
        let mut hits = Vec::with_capacity(count);
        for _ in 0..count {
            x += dx_min + i64::from(generator.next_u32() % dx_span);
            y += dy_min + i64::from(generator.next_u32() % dy_span);
            let span = span_min + generator.next_u32() % span_span;
            let coverage = coverage_min + generator.next_u32() % coverage_span;
            hits.push(KmerHit {
                read_identifier: 1,
                reverse_strand: false,
                offset: y as u32,
                self_offset: x as u32,
                count: (coverage << 8) | span,
            });
        }
        hits
    }

    fn run_chain_scenario(scenario: &ChainScenario, lines: &mut Vec<String>) {
        let mut hits = generate_chain_hits(
            scenario.count,
            scenario.seed,
            100,
            100,
            scenario.dx_min,
            scenario.dx_span,
            scenario.dy_min,
            scenario.dy_span,
            scenario.span_min,
            scenario.span_span,
            scenario.coverage_min,
            scenario.coverage_span,
        );
        let mut data = ChainData::default();
        let ret = chain_dp_fixed_ends(
            &mut hits,
            scenario.count as i64,
            &mut data,
            scenario.maximum_skip,
            scenario.maximum_iterations,
            scenario.maximum_distance,
            0.5,
            0.0005,
            scenario.band_width_rate,
            scenario.query_length,
            scenario.target_length,
            scenario.quick,
            scenario.left,
            scenario.right,
        );
        lines.push(format!(
            "CH name={} quick={} left={} right={} n={} ret={}",
            scenario.name,
            i32::from(scenario.quick),
            i32::from(scenario.left),
            i32::from(scenario.right),
            scenario.count,
            ret as u64
        ));
        push_hits(lines, "CHA", &hits, ret);
        for index in 0..ret as usize {
            lines.push(format!("CHT i={} t={}", index, data.scratch[index]));
        }
    }

    fn run_simple_scenario(
        lines: &mut Vec<String>,
        name: &str,
        mut hits: Vec<KmerHit>,
        maximum_skip: i64,
        maximum_iterations: i64,
    ) {
        let mut data = ChainData::default();
        let hit_count = hits.len() as i64;
        let ret = chain_simple(&mut hits, hit_count, &mut data, maximum_skip, maximum_iterations);
        lines.push(format!(
            "CS name={} n={} skip={} iter={} ret={}",
            name,
            hits.len(),
            maximum_skip,
            maximum_iterations,
            ret as u64
        ));
        push_hits(lines, "CSA", &hits, ret);
        for index in 0..ret as usize {
            lines.push(format!(
                "CSF i={} f={} pre={}",
                index, data.score[index], data.predecessor[index]
            ));
        }
    }

    // --------------------------------------------------------- exact cigar

    struct ExactScenario {
        name: &'static str,
        cells: Vec<(u8, u32)>,
        pattern_start: i32,
        text_start: i32,
        pattern_min: i32,
        pattern_max: i32,
        text_min: i32,
        text_max: i32,
        minimum_run: i32,
        minimum_window_length: i64,
        h_kmer_hit: i64,
    }

    fn run_exact_scenario(scenario: &ExactScenario, lines: &mut Vec<String>) {
        let packed: Vec<u16> = scenario
            .cells
            .iter()
            .map(|&(operation, length)| ((operation as u16) << 14) | (length as u16))
            .collect();
        let mut candidates = CandidatesList::default();
        let occurrence = extract_exact_cigar(
            &packed,
            scenario.pattern_start,
            scenario.text_start,
            scenario.pattern_min,
            scenario.pattern_max,
            scenario.text_min,
            scenario.text_max,
            &mut candidates,
            scenario.minimum_run,
            scenario.minimum_window_length,
            scenario.h_kmer_hit,
        );
        lines.push(format!(
            "EX name={} occ={} n={}",
            scenario.name, occurrence, candidates.length
        ));
        push_hits(lines, "EXH", &candidates.list, candidates.length);
    }

    // ------------------------------------------------------- single anchors

    struct SingleScenario {
        name: &'static str,
        mode: i64,
        h_kmer_hit: i64,
        query_start: i64,
        query_end: i64,
        target_start: i64,
        target_end: i64,
        maximum_skip: i64,
        maximum_iterations: i64,
        raw_hits: Vec<(u32, u32)>,
        spans: Vec<u32>,
    }

    fn run_single_scenario(scenario: &SingleScenario, lines: &mut Vec<String>) {
        let mut candidates = CandidatesList::default();
        for index in 0..3u32 {
            candidates.push(KmerHit {
                read_identifier: 1,
                reverse_strand: false,
                offset: 2000 + index,
                self_offset: 1000 + index,
                count: 1,
            });
        }
        let base = candidates.length;
        for (index, &(self_offset, offset)) in scenario.raw_hits.iter().enumerate() {
            candidates.push(KmerHit {
                read_identifier: 1,
                reverse_strand: false,
                offset,
                self_offset,
                count: scenario.spans[index],
            });
        }
        // The function reads its input at `candidates.length`, which the caller
        // rewinds to the start of the region to rebuild.
        candidates.length = base;
        let ret = generate_single_kmer_hit_chain(
            &mut candidates,
            scenario.raw_hits.len() as i64,
            scenario.h_kmer_hit,
            scenario.mode,
            scenario.query_start,
            scenario.query_end,
            scenario.target_start,
            scenario.target_end,
            scenario.maximum_skip,
            scenario.maximum_iterations,
        );
        lines.push(format!(
            "GS name={} mode={} hk={} ch_n={} skip={} iter={} qs={} qe={} ts={} te={} ret={} length={}",
            scenario.name,
            scenario.mode,
            scenario.h_kmer_hit,
            scenario.raw_hits.len(),
            scenario.maximum_skip,
            scenario.maximum_iterations,
            scenario.query_start,
            scenario.query_end,
            scenario.target_start,
            scenario.target_end,
            ret,
            candidates.length
        ));
        push_hits(lines, "GSH", &candidates.list[base as usize..], ret);
    }

    // --------------------------------------------------------------- reads

    struct WindowSpec {
        x_start: i32,
        x_end: i32,
        y_start: i32,
        y_end: i32,
        error: i16,
        threshold: i16,
    }

    struct ReadScenario {
        name: &'static str,
        reads: Vec<Vec<u8>>,
        x_start: u32,
        x_end: u32,
        y_start: u32,
        y_end: u32,
        y_strand: u32,
        non_homopolymer_errors: u32,
        gaps: Vec<(u32, i32)>,
        windows: Vec<WindowSpec>,
        hits: Vec<(u32, u32)>,
        modes: Vec<i64>,
        query_ranges: Vec<(i64, i64)>,
        target_ranges: Vec<(i64, i64)>,
        window_length: i64,
        h_kmer_hits: Vec<i64>,
    }

    fn build_base(scenario: &ReadScenario) -> OverlapRegion {
        let mut region = OverlapRegion {
            x_identifier: 0,
            x_position_start: scenario.x_start,
            x_position_end: scenario.x_end,
            y_identifier: 1,
            y_position_start: scenario.y_start,
            y_position_end: scenario.y_end,
            y_strand: scenario.y_strand,
            overlap_length: scenario.x_end + 1 - scenario.x_start,
            non_homopolymer_errors: scenario.non_homopolymer_errors,
            ..OverlapRegion::default()
        };
        for &(site, shift) in &scenario.gaps {
            region.fake_cigar.push(site, shift);
        }
        for window in &scenario.windows {
            region.window_cigar.append_window(
                window.x_start,
                window.x_end,
                window.y_start,
                window.y_end,
                window.error,
                0,
                0,
                window.threshold,
            );
        }
        region
    }

    fn build_reads(scenario: &ReadScenario) -> ReadStore {
        let mut reads = ReadStore::default();
        for (index, sequence) in scenario.reads.iter().enumerate() {
            reads.push_read(&format!("read{}", index), sequence);
        }
        reads
    }

    fn produce_alignment(
        base: &OverlapRegion,
        scenario: &ReadScenario,
        reads: &ReadStore,
        extz: &mut BitExtz,
    ) -> OverlapRegion {
        let hits: Vec<KmerHit> = scenario
            .hits
            .iter()
            .map(|&(self_offset, offset)| KmerHit {
                self_offset,
                offset,
                ..KmerHit::default()
            })
            .collect();
        let query = &scenario.reads[0];
        let query_total_length = query.len() as i64;
        let target_total_length = scenario.reads[1].len() as i64;
        let mut output = OverlapRegion::default();
        hc_ovlp_base_direct(
            base,
            &hits,
            hits.len() as i64,
            scenario.window_length,
            reads,
            query,
            extz,
            &mut output.window_cigar,
            0.04,
            query_total_length,
            target_total_length,
            -1,
        );
        output
    }

    fn run_read_scenario(scenario: &ReadScenario, lines: &mut Vec<String>) {
        let reads = build_reads(scenario);
        let base = build_base(scenario);
        let mut extz = BitExtz::new();
        let mut reference = produce_alignment(&base, scenario, &reads, &mut extz);

        let query = &scenario.reads[0];
        let query_total_length = query.len() as i64;
        let target_total_length = scenario.reads[1].len() as i64;

        for index in 0..scenario.modes.len() {
            let (query_start, query_end) = scenario.query_ranges[index];
            let (target_start, target_end) = scenario.target_ranges[index];
            let h_kmer_hit = scenario.h_kmer_hits[index];
            let mut candidates = CandidatesList::default();
            let ret = generate_window_chain(
                &mut reference,
                &mut candidates,
                query_start,
                query_end,
                target_start,
                target_end,
                scenario.window_length,
                &reads,
                query,
                &mut extz,
                query_total_length,
                target_total_length,
                0.04,
                h_kmer_hit,
                scenario.modes[index],
                true,
            );
            lines.push(format!(
                "GW name={} mode={} qs={} qe={} ts={} te={} wl={} hk={} ql={} tl={} ret={} length={}",
                scenario.name,
                scenario.modes[index],
                query_start,
                query_end,
                target_start,
                target_end,
                scenario.window_length,
                h_kmer_hit,
                query_total_length,
                target_total_length,
                ret,
                candidates.length
            ));
            push_hits(lines, "GWH", &candidates.list, ret);
            push_windows(lines, "GWW", &mut reference.window_cigar);
        }
    }

    fn run_rechain(
        scenario: &ReadScenario,
        run_name: &str,
        auxiliary_index: usize,
        mut shift: i64,
        mode: i64,
        lines: &mut Vec<String>,
    ) {
        let reads = build_reads(scenario);
        let mut base = build_base(scenario);
        let mut extz = BitExtz::new();
        let reference = produce_alignment(&base, scenario, &reads, &mut extz);

        // The retry pass starts from the window list the first alignment
        // produced; `gen_hc_fast_cigar0` leaves the windows it could not align
        // marked as unaligned and hands the whole list to `rechain_alignment_hc`.
        // The overlap region `z` that the retry consults is still the original
        // one, whose window list is the tiling the first pass consumed.
        let mut aux = OverlapRegion {
            x_identifier: base.x_identifier,
            x_position_start: base.x_position_start,
            x_position_end: base.x_position_end,
            x_strand: base.x_strand,
            y_identifier: base.y_identifier,
            y_position_start: base.y_position_start,
            y_position_end: base.y_position_end,
            y_strand: base.y_strand,
            window_cigar: reference.window_cigar.clone(),
            ..OverlapRegion::default()
        };

        let (query_start, query_end, target_start, target_end);
        {
            let window = &mut aux.window_cigar.windows[auxiliary_index];
            if i64::from(window.y_start) + shift < 0 {
                shift = 0;
            }
            query_start = window.x_start;
            query_end = window.x_end;
            target_start = window.y_start + shift as i32;
            target_end = window.y_end + shift as i32;
            window.y_start = target_start;
            window.y_end = target_end;
            window.error = i16::MAX;
            window.extra_begin = -1;
            window.extra_end = -1;
            window.cigar_start = 0;
            window.cigar_length = 0;
            window.error_threshold = mode as i16;
        }

        let query = &scenario.reads[0];
        let query_total_length = query.len() as i64;
        let target_total_length = scenario.reads[1].len() as i64;
        let mut candidates = CandidatesList::default();
        let before = aux.window_cigar.windows.len() as i64;
        rechain_alignment_hc(
            &mut base,
            &mut candidates,
            &mut aux.window_cigar,
            auxiliary_index,
            scenario.window_length,
            &reads,
            query,
            &mut extz,
            0.04,
            query_total_length,
            target_total_length,
            10,
        );
        lines.push(format!(
            "RC name={} aux_i={} mode={} qs={} qe={} ts={} te={} wl={} hk=10 before={} length={}",
            run_name,
            auxiliary_index,
            mode,
            query_start,
            query_end,
            target_start,
            target_end,
            scenario.window_length,
            before,
            candidates.length
        ));
        push_windows(lines, "RCW", &mut aux.window_cigar);
    }

    // ------------------------------------------------------------ scenarios

    fn chain_scenarios() -> Vec<ChainScenario> {
        vec![
            ChainScenario {
                name: "fix_global",
                count: 30,
                seed: 101,
                dx_min: 5,
                dx_span: 40,
                dy_min: 5,
                dy_span: 40,
                span_min: 1,
                span_span: 20,
                coverage_min: 1,
                coverage_span: 4,
                quick: false,
                left: true,
                right: true,
                maximum_skip: 5,
                maximum_iterations: 200,
                maximum_distance: 5000,
                query_length: 10000,
                target_length: 10000,
                band_width_rate: 0.01,
            },
            ChainScenario {
                name: "fix_global_quick",
                count: 30,
                seed: 101,
                dx_min: 5,
                dx_span: 40,
                dy_min: 5,
                dy_span: 40,
                span_min: 1,
                span_span: 20,
                coverage_min: 1,
                coverage_span: 4,
                quick: true,
                left: true,
                right: true,
                maximum_skip: 5,
                maximum_iterations: 200,
                maximum_distance: 5000,
                query_length: 10000,
                target_length: 10000,
                band_width_rate: 0.01,
            },
            ChainScenario {
                name: "fix_left",
                count: 30,
                seed: 102,
                dx_min: 5,
                dx_span: 40,
                dy_min: 5,
                dy_span: 40,
                span_min: 1,
                span_span: 20,
                coverage_min: 1,
                coverage_span: 4,
                quick: false,
                left: true,
                right: false,
                maximum_skip: 5,
                maximum_iterations: 200,
                maximum_distance: 5000,
                query_length: 10000,
                target_length: 10000,
                band_width_rate: 0.01,
            },
            ChainScenario {
                name: "fix_right",
                count: 30,
                seed: 103,
                dx_min: 5,
                dx_span: 40,
                dy_min: 5,
                dy_span: 40,
                span_min: 1,
                span_span: 20,
                coverage_min: 1,
                coverage_span: 4,
                quick: false,
                left: false,
                right: true,
                maximum_skip: 5,
                maximum_iterations: 200,
                maximum_distance: 5000,
                query_length: 10000,
                target_length: 10000,
                band_width_rate: 0.01,
            },
            ChainScenario {
                name: "fix_right_quick",
                count: 30,
                seed: 103,
                dx_min: 5,
                dx_span: 40,
                dy_min: 5,
                dy_span: 40,
                span_min: 1,
                span_span: 20,
                coverage_min: 1,
                coverage_span: 4,
                quick: true,
                left: false,
                right: true,
                maximum_skip: 5,
                maximum_iterations: 200,
                maximum_distance: 5000,
                query_length: 10000,
                target_length: 10000,
                band_width_rate: 0.01,
            },
            ChainScenario {
                name: "fix_right_colinear",
                count: 24,
                seed: 104,
                dx_min: 10,
                dx_span: 1,
                dy_min: 10,
                dy_span: 1,
                span_min: 5,
                span_span: 4,
                coverage_min: 1,
                coverage_span: 2,
                quick: true,
                left: false,
                right: true,
                maximum_skip: 5,
                maximum_iterations: 200,
                maximum_distance: 5000,
                query_length: 10000,
                target_length: 10000,
                band_width_rate: 0.01,
            },
            ChainScenario {
                name: "fix_global_colinear",
                count: 24,
                seed: 104,
                dx_min: 10,
                dx_span: 1,
                dy_min: 10,
                dy_span: 1,
                span_min: 5,
                span_span: 4,
                coverage_min: 1,
                coverage_span: 2,
                quick: true,
                left: true,
                right: true,
                maximum_skip: 5,
                maximum_iterations: 200,
                maximum_distance: 5000,
                query_length: 10000,
                target_length: 10000,
                band_width_rate: 0.01,
            },
            ChainScenario {
                name: "fix_narrow",
                count: 40,
                seed: 105,
                dx_min: 5,
                dx_span: 40,
                dy_min: 5,
                dy_span: 40,
                span_min: 1,
                span_span: 20,
                coverage_min: 1,
                coverage_span: 4,
                quick: false,
                left: true,
                right: true,
                maximum_skip: 0,
                maximum_iterations: 8,
                maximum_distance: 200,
                query_length: 10000,
                target_length: 10000,
                band_width_rate: 0.001,
            },
            ChainScenario {
                name: "fix_single",
                count: 1,
                seed: 106,
                dx_min: 5,
                dx_span: 40,
                dy_min: 5,
                dy_span: 40,
                span_min: 1,
                span_span: 20,
                coverage_min: 1,
                coverage_span: 4,
                quick: false,
                left: true,
                right: true,
                maximum_skip: 5,
                maximum_iterations: 200,
                maximum_distance: 5000,
                query_length: 10000,
                target_length: 10000,
                band_width_rate: 0.01,
            },
            ChainScenario {
                name: "fix_two_right",
                count: 2,
                seed: 107,
                dx_min: 5,
                dx_span: 40,
                dy_min: 5,
                dy_span: 40,
                span_min: 1,
                span_span: 20,
                coverage_min: 1,
                coverage_span: 4,
                quick: false,
                left: false,
                right: true,
                maximum_skip: 5,
                maximum_iterations: 200,
                maximum_distance: 5000,
                query_length: 10000,
                target_length: 10000,
                band_width_rate: 0.01,
            },
            ChainScenario {
                name: "fix_right_noquick_narrow",
                count: 40,
                seed: 108,
                dx_min: 5,
                dx_span: 40,
                dy_min: 5,
                dy_span: 40,
                span_min: 1,
                span_span: 20,
                coverage_min: 1,
                coverage_span: 4,
                quick: false,
                left: false,
                right: true,
                maximum_skip: 2,
                maximum_iterations: 16,
                maximum_distance: 300,
                query_length: 10000,
                target_length: 10000,
                band_width_rate: 0.01,
            },
        ]
    }

    fn exact_scenarios() -> Vec<ExactScenario> {
        vec![
            ExactScenario {
                name: "exact_simple",
                cells: vec![(0, 300)],
                pattern_start: 0,
                text_start: 0,
                pattern_min: 0,
                pattern_max: 1000,
                text_min: 0,
                text_max: 1000,
                minimum_run: 10,
                minimum_window_length: 100,
                h_kmer_hit: 20,
            },
            ExactScenario {
                name: "exact_clip",
                cells: vec![
                    (0, 100),
                    (1, 5),
                    (0, 80),
                    (2, 3),
                    (0, 120),
                    (3, 4),
                    (0, 90),
                ],
                pattern_start: 100,
                text_start: 200,
                pattern_min: 150,
                pattern_max: 350,
                text_min: 230,
                text_max: 480,
                minimum_run: 10,
                minimum_window_length: 64,
                h_kmer_hit: 20,
            },
            ExactScenario {
                name: "exact_short_fallback",
                cells: vec![(0, 8), (1, 40), (0, 5)],
                pattern_start: 0,
                text_start: 0,
                pattern_min: 0,
                pattern_max: 1000,
                text_min: 0,
                text_max: 1000,
                minimum_run: 10,
                minimum_window_length: 20,
                h_kmer_hit: 20,
            },
            ExactScenario {
                name: "exact_short_nofallback",
                cells: vec![(0, 8), (1, 40), (0, 5)],
                pattern_start: 0,
                text_start: 0,
                pattern_min: 0,
                pattern_max: 1000,
                text_min: 0,
                text_max: 1000,
                minimum_run: 10,
                minimum_window_length: 100,
                h_kmer_hit: 20,
            },
            ExactScenario {
                name: "exact_start_clip",
                cells: vec![(0, 40), (1, 2), (0, 60)],
                pattern_start: 0,
                text_start: 5,
                pattern_min: 0,
                pattern_max: 500,
                text_min: 10,
                text_max: 500,
                minimum_run: 10,
                minimum_window_length: 20,
                h_kmer_hit: 1000,
            },
            ExactScenario {
                name: "exact_ops",
                cells: vec![
                    (0, 50),
                    (2, 10),
                    (0, 60),
                    (3, 7),
                    (0, 70),
                    (1, 3),
                    (0, 80),
                ],
                pattern_start: 300,
                text_start: 400,
                pattern_min: 320,
                pattern_max: 700,
                text_min: 420,
                text_max: 800,
                minimum_run: 4,
                minimum_window_length: 32,
                h_kmer_hit: 40,
            },
            ExactScenario {
                name: "exact_split_runs",
                cells: vec![(0, 120), (1, 1), (0, 140)],
                pattern_start: 0,
                text_start: 0,
                pattern_min: 0,
                pattern_max: 400,
                text_min: 0,
                text_max: 400,
                minimum_run: 10,
                minimum_window_length: 64,
                h_kmer_hit: 20,
            },
        ]
    }

    fn single_scenarios() -> Vec<SingleScenario> {
        let shared_hits = vec![(199, 299), (259, 359), (319, 419), (420, 520)];
        let shared_spans = vec![11, 20, 5, 11];
        vec![
            SingleScenario {
                name: "gs_mode0",
                mode: 0,
                h_kmer_hit: 10,
                query_start: 100,
                query_end: 500,
                target_start: 200,
                target_end: 600,
                maximum_skip: 5,
                maximum_iterations: 200,
                raw_hits: shared_hits.clone(),
                spans: shared_spans.clone(),
            },
            SingleScenario {
                name: "gs_mode1",
                mode: 1,
                h_kmer_hit: 10,
                query_start: 100,
                query_end: 500,
                target_start: 200,
                target_end: 600,
                maximum_skip: 5,
                maximum_iterations: 200,
                raw_hits: shared_hits.clone(),
                spans: shared_spans.clone(),
            },
            SingleScenario {
                name: "gs_mode2",
                mode: 2,
                h_kmer_hit: 10,
                query_start: 100,
                query_end: 500,
                target_start: 200,
                target_end: 600,
                maximum_skip: 5,
                maximum_iterations: 200,
                raw_hits: shared_hits.clone(),
                spans: shared_spans.clone(),
            },
            SingleScenario {
                name: "gs_unsorted",
                mode: 0,
                h_kmer_hit: 10,
                query_start: 100,
                query_end: 500,
                target_start: 200,
                target_end: 600,
                maximum_skip: 5,
                maximum_iterations: 200,
                raw_hits: vec![
                    (300, 400),
                    (199, 299),
                    (260, 360),
                    (200, 300),
                    (420, 520),
                    (319, 419),
                ],
                spans: vec![11, 12, 20, 9, 5, 11],
            },
            SingleScenario {
                name: "gs_duplicates",
                mode: 0,
                h_kmer_hit: 10,
                query_start: 100,
                query_end: 500,
                target_start: 200,
                target_end: 600,
                maximum_skip: 5,
                maximum_iterations: 200,
                raw_hits: vec![(199, 299), (199, 299), (199, 299), (300, 400), (420, 520)],
                spans: vec![11, 11, 4, 15, 5],
            },
            SingleScenario {
                name: "gs_zero_spans",
                mode: 0,
                h_kmer_hit: 10,
                query_start: 100,
                query_end: 500,
                target_start: 200,
                target_end: 600,
                maximum_skip: 5,
                maximum_iterations: 200,
                raw_hits: vec![(100, 200), (199, 299), (250, 350), (500, 600)],
                spans: vec![0, 11, 0, 1],
            },
            SingleScenario {
                name: "gs_boundary",
                mode: 0,
                h_kmer_hit: 10,
                query_start: 100,
                query_end: 500,
                target_start: 200,
                target_end: 600,
                maximum_skip: 5,
                maximum_iterations: 200,
                raw_hits: vec![(100, 200), (250, 350), (499, 599), (500, 600)],
                spans: vec![11, 20, 4, 1],
            },
            SingleScenario {
                name: "gs_all_zero",
                mode: 0,
                h_kmer_hit: 10,
                query_start: 100,
                query_end: 500,
                target_start: 200,
                target_end: 600,
                maximum_skip: 5,
                maximum_iterations: 200,
                raw_hits: vec![(100, 200), (250, 350)],
                spans: vec![0, 0],
            },
            SingleScenario {
                name: "gs_large_hk",
                mode: 0,
                h_kmer_hit: 1000,
                query_start: 100,
                query_end: 500,
                target_start: 200,
                target_end: 600,
                maximum_skip: 5,
                maximum_iterations: 200,
                raw_hits: shared_hits,
                spans: shared_spans,
            },
            SingleScenario {
                name: "gs_negative_start",
                mode: 0,
                h_kmer_hit: 10,
                query_start: 100,
                query_end: 500,
                target_start: 200,
                target_end: 600,
                maximum_skip: 5,
                maximum_iterations: 200,
                raw_hits: vec![(105, 205), (259, 359), (319, 419)],
                spans: vec![11, 20, 5],
            },
        ]
    }

    fn read_scenarios() -> Vec<ReadScenario> {
        let mut generator = Lcg(987654321);
        let query = generator.sequence(2000);
        let prefix = generator.sequence(600);
        let suffix = generator.sequence(400);
        let tile_query = generator.sequence(1400);
        let mut exact_target = prefix.clone();
        exact_target.extend_from_slice(&query);
        exact_target.extend_from_slice(&suffix);
        let mut substituted = exact_target.clone();
        for index in (0..query.len()).step_by(97) {
            let position = prefix.len() + index;
            substituted[position] = if substituted[position] == b'A' {
                b'C'
            } else {
                b'A'
            };
        }
        let reversed = reverse_complement(&exact_target);

        let exact_tile: Vec<WindowSpec> = (0..4)
            .map(|index| {
                let start = index * 350;
                WindowSpec {
                    x_start: start,
                    x_end: start + 349,
                    y_start: start,
                    y_end: start + 349,
                    error: 0,
                    threshold: 0,
                }
            })
            .collect();
        let subs_windows = vec![
            WindowSpec {
                x_start: 200,
                x_end: 599,
                y_start: 800,
                y_end: 1199,
                error: 4,
                threshold: 0,
            },
            WindowSpec {
                x_start: 600,
                x_end: 999,
                y_start: 1200,
                y_end: 1599,
                error: 4,
                threshold: 0,
            },
            WindowSpec {
                x_start: 1000,
                x_end: 1399,
                y_start: 1600,
                y_end: 1999,
                error: 4,
                threshold: 0,
            },
            WindowSpec {
                x_start: 1400,
                x_end: 1799,
                y_start: 2000,
                y_end: 2399,
                error: 5,
                threshold: 0,
            },
        ];
        let reverse_windows = vec![
            WindowSpec {
                x_start: 200,
                x_end: 599,
                y_start: 800,
                y_end: 1199,
                error: 0,
                threshold: 0,
            },
            WindowSpec {
                x_start: 600,
                x_end: 999,
                y_start: 1200,
                y_end: 1599,
                error: 0,
                threshold: 0,
            },
            WindowSpec {
                x_start: 1000,
                x_end: 1399,
                y_start: 1600,
                y_end: 1999,
                error: 0,
                threshold: 0,
            },
            WindowSpec {
                x_start: 1400,
                x_end: 1799,
                y_start: 2000,
                y_end: 2399,
                error: 0,
                threshold: 0,
            },
        ];

        vec![
            ReadScenario {
                name: "W1_exact",
                reads: vec![tile_query.clone(), tile_query.clone()],
                x_start: 0,
                x_end: 1399,
                y_start: 0,
                y_end: 1399,
                y_strand: 0,
                non_homopolymer_errors: 0,
                gaps: vec![(0, 0), (1399, 0)],
                windows: exact_tile,
                hits: vec![],
                modes: vec![0, 1, 2],
                query_ranges: vec![(0, 1400), (350, 1050), (0, 1400)],
                target_ranges: vec![(0, 1400), (350, 1050), (0, 1400)],
                window_length: 350,
                h_kmer_hits: vec![10, 10, 1000],
            },
            ReadScenario {
                name: "W2_subs",
                reads: vec![query.clone(), substituted],
                x_start: 200,
                x_end: 1799,
                y_start: 800,
                y_end: 2399,
                y_strand: 0,
                non_homopolymer_errors: 16,
                gaps: vec![(200, 0), (1799, 0)],
                windows: subs_windows,
                hits: vec![(200, 800), (600, 1200), (1000, 1600), (1400, 2000)],
                modes: vec![0, 1, 2],
                query_ranges: vec![(200, 1800), (600, 1400), (200, 1800)],
                target_ranges: vec![(800, 2400), (1200, 2000), (800, 2400)],
                window_length: 400,
                h_kmer_hits: vec![10, 10, 1000],
            },
            ReadScenario {
                name: "W3_reverse",
                reads: vec![query, reversed],
                x_start: 200,
                x_end: 1799,
                y_start: 800,
                y_end: 2399,
                y_strand: 1,
                non_homopolymer_errors: 0,
                gaps: vec![(200, 0), (1799, 0)],
                windows: reverse_windows,
                hits: vec![(200, 800), (600, 1200), (1000, 1600), (1400, 2000)],
                modes: vec![0, 1],
                query_ranges: vec![(200, 1800), (600, 1400)],
                target_ranges: vec![(800, 2400), (1200, 2000)],
                window_length: 400,
                h_kmer_hits: vec![10, 1000],
            },
        ]
    }

    fn generate_lines() -> Vec<String> {
        let mut lines = Vec::new();

        for scenario in chain_scenarios() {
            run_chain_scenario(&scenario, &mut lines);
        }

        let mut simple_hits = generate_chain_hits(16, 201, 100, 100, 10, 1, 10, 1, 5, 6, 1, 3);
        run_simple_scenario(&mut lines, "simple_colinear", simple_hits.clone(), 5, 5000);
        simple_hits[8].offset = simple_hits[7].offset.wrapping_sub(5);
        run_simple_scenario(&mut lines, "simple_break", simple_hits.clone(), 5, 5000);
        run_simple_scenario(&mut lines, "simple_break_small", simple_hits, 0, 4);
        let single = vec![KmerHit {
            read_identifier: 1,
            reverse_strand: false,
            offset: 200,
            self_offset: 100,
            count: (3 << 8) | 7,
        }];
        run_simple_scenario(&mut lines, "simple_single", single, 5, 5000);

        for scenario in exact_scenarios() {
            run_exact_scenario(&scenario, &mut lines);
        }

        for scenario in single_scenarios() {
            run_single_scenario(&scenario, &mut lines);
        }

        let read_scenarios = read_scenarios();
        for scenario in &read_scenarios {
            run_read_scenario(scenario, &mut lines);
        }

        let rechain_specs: [(&str, usize, i64, i64); 4] = [
            ("keep_mode0", 0, 0, 0),
            ("shift_mode0", 0, -7, 0),
            ("keep_mode1", 0, 0, 1),
            ("keep_mode2", 0, 0, 2),
        ];
        for scenario in &read_scenarios {
            for (suffix, auxiliary_index, shift, mode) in rechain_specs {
                let run_name = format!("{}_{}", scenario.name, suffix);
                run_rechain(
                    scenario,
                    &run_name,
                    auxiliary_index,
                    shift,
                    mode,
                    &mut lines,
                );
            }
        }

        lines
    }

    #[test]
    #[ignore = "requires reference/ec_chain_probe.out from the C++ probe"]
    fn window_rechaining_matches_the_reference() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("reference")
            .join("ec_chain_probe.out");
        let contents = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("cannot read {}: {}", path.display(), error));
        let golden: Vec<&str> = contents.lines().collect();
        let generated = generate_lines();

        assert_eq!(
            generated.len(),
            golden.len(),
            "the probe printed {} lines, the replay generated {}",
            golden.len(),
            generated.len()
        );
        for (index, (expected, actual)) in golden.iter().zip(generated.iter()).enumerate() {
            assert_eq!(
                expected, actual,
                "line {} differs\nexpected: {}\nactual:   {}",
                index + 1,
                expected,
                actual
            );
        }
    }
}

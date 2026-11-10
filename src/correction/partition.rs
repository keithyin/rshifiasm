//! The legacy window-level phasing of the `--sec-in` pair mode:
//! `partition_overlaps_advance`, `cluster_advance` and the
//! `markSNP`/`addSNPtohaplotype` helpers of `Correct.cpp`.
//!
//! The pair pipeline calls `correct_overlap` after the overlap windows have
//! been verified and recalculated. Unlike the newer `rphase_hc` path, this
//! code marks mismatches into a per-window flag array, folds them into
//! haplotype evidence and finally re-labels large-indel overlaps.

use crate::correction::consensus::{get_available_interval, is_strict_homopolymer};
use crate::correction::coverage::available_fully_covered_interval;
use crate::correction::haplotype_evidence::{
    generate_haplotypes_naive_hifi, HaplotypeEvidence, HaplotypeEvidenceAlloc, SnpStats,
};
use crate::correction::interval::{IntervalBuffer, WindowPool, INTERVAL_PAST_THE_END};
use crate::correction::region::WINDOW_LENGTH;
use crate::correction::window::WindowCigar;
use crate::container::ksort::radix_sort_by_key_32;
use crate::overlap::region::OverlapRegion;
use crate::reads::read_store::ReadStore;
use crate::sequence::encoding::SEQUENCE_NT6_TABLE;

const SNP_THRESHOLD: i64 = 1;

/// `markSNP_detail`: adds one to the mismatch counter of every covered query
/// position of one window cigar.
#[allow(clippy::too_many_arguments)]
fn mark_snp_detail(
    cigar: &WindowCigar,
    window_index: usize,
    flag: &mut [u8],
    absolute_base: i64,
    x_begin: i64,
    x_end: i64,
    flag_offset: i64,
) {
    let mut x_position = 0i64;
    let mut y_position = 0i64;
    for (operation, length) in cigar.window_cells(window_index) {
        if x_position > x_end {
            break;
        }
        let length = i64::from(length);
        match operation {
            0 => {
                x_position += length;
                y_position += length;
            }
            1 => {
                for _ in 0..length {
                    if x_position >= x_begin && x_position <= x_end {
                        let index = absolute_base + x_position - flag_offset;
                        if index >= 0 && (index as usize) < flag.len() && flag[index as usize] < 127 {
                            flag[index as usize] += 1;
                        }
                    }
                    x_position += 1;
                    y_position += 1;
                }
            }
            2 => y_position += length,
            3 => x_position += length,
            _ => {}
        }
    }
}

/// `markSNP_advance`: marks the mismatches of one overlap window, including
/// the parts that spill into the neighbouring boundary cigars.
#[allow(clippy::too_many_arguments)]
fn mark_snp_advance(
    window_start: i64,
    x_total_start: i64,
    x_length: i64,
    _y_total_start: i64,
    _y_length: i64,
    current: &WindowCigar,
    current_window: usize,
    beg_cigar: Option<(&WindowCigar, usize)>,
    end_cigar: Option<(&WindowCigar, usize)>,
    hap: &mut HaplotypeEvidenceAlloc,
) {
    let x_total_end = x_total_start + x_length - 1;
    let inner_offset = x_total_start - window_start;
    let mut current_begin = 0i64;
    let mut current_end = x_length - 1;
    // The reference reads `hap->flag` while it appends to `hap->list`; the
    // flag array is small, so a local copy keeps the borrows apart.
    let flag = hap.flag;

    if let Some((cigar, index)) = beg_cigar {
        let window = &cigar.windows[index];
        if window.y_end != -1 {
            let left_useless = i64::from(window.extra_begin);
            let right_useless = i64::from(window.extra_end);
            let x_left = x_total_start - i64::from(window.x_start);
            let x_right = i64::from(window.x_end) - x_total_start + 1;
            if x_left > left_useless && x_right > right_useless {
                let interval_begin = x_left;
                let interval_end = interval_begin + (x_right - right_useless) - 1;
                current_begin = x_right - right_useless;
                mark_snp_detail(
                    cigar,
                    index,
                    &mut hap.flag,
                    inner_offset,
                    interval_begin,
                    interval_end,
                    interval_begin,
                );
            }
        }
    }

    if let Some((cigar, index)) = end_cigar {
        let window = &cigar.windows[index];
        if window.y_end != -1 {
            let left_useless = i64::from(window.extra_begin);
            let right_useless = i64::from(window.extra_end);
            let x_left = x_total_end - i64::from(window.x_start);
            let x_right = i64::from(window.x_end) - x_total_end + 1;
            if x_left > left_useless && x_right > right_useless {
                let interval_end = x_left;
                let interval_begin = interval_end - ((x_left + 1) - left_useless) + 1;
                current_end = (x_length - 1) - ((x_left + 1) - left_useless);
                mark_snp_detail(
                    cigar,
                    index,
                    &mut hap.flag,
                    i64::from(window.x_start) - window_start,
                    interval_begin,
                    interval_end,
                    0,
                );
            }
        }
    }

    mark_snp_detail(
        current,
        current_window,
        &mut hap.flag,
        inner_offset,
        current_begin,
        current_end,
        0,
    );
}

/// `addSNPtohaplotype_details`: records the evidence of the marked positions
/// of one window cigar.
#[allow(clippy::too_many_arguments)]
fn add_snp_to_haplotype_details(
    cigar: &WindowCigar,
    window_index: usize,
    flag: &[u8],
    absolute_base: i64,
    x_string: &[u8],
    y_string: &[u8],
    x_total_start: i64,
    y_total_start: i64,
    x_begin: i64,
    x_end: i64,
    overlap_identifier: u32,
    flag_offset: i64,
    hap: &mut HaplotypeEvidenceAlloc,
    snp_threshold: i64,
) {
    if x_begin > x_end {
        return;
    }
    let mut x_position = 0i64;
    let mut y_position = 0i64;

    for (operation, length) in cigar.window_cells(window_index) {
        if x_position > x_end {
            break;
        }
        let length = i64::from(length);
        let mut push = |hap: &mut HaplotypeEvidenceAlloc,
                        site: i64,
                        overlap_site: i64,
                        kind: u8,
                        base: u8| {
            hap.add(HaplotypeEvidence {
                site: site as u32,
                overlap_identifier,
                overlap_site: overlap_site as u32,
                evidence_type: kind,
                coverage: 1,
                mismatch_base: base,
            });
        };
        let index = absolute_base + x_position - flag_offset;
        let marked = index >= 0
            && (index as usize) < flag.len()
            && flag[index as usize] > snp_threshold as u8;
        match operation {
            0 => {
                for _ in 0..length {
                    let index = absolute_base + x_position - flag_offset;
                    let marked = index >= 0
                        && (index as usize) < flag.len()
                        && flag[index as usize] > snp_threshold as u8;
                    if x_position >= x_begin && x_position <= x_end && marked {
                        let base = y_string.get(y_position as usize).copied().unwrap_or(b'N');
                        push(hap, x_total_start + x_position, y_total_start + y_position, 0, base);
                    }
                    x_position += 1;
                    y_position += 1;
                }
            }
            1 => {
                for _ in 0..length {
                    let index = absolute_base + x_position - flag_offset;
                    let marked = index >= 0
                        && (index as usize) < flag.len()
                        && flag[index as usize] > snp_threshold as u8;
                    if x_position >= x_begin && x_position <= x_end && marked {
                        let base = y_string.get(y_position as usize).copied().unwrap_or(b'N');
                        push(hap, x_total_start + x_position, y_total_start + y_position, 1, base);
                    }
                    x_position += 1;
                    y_position += 1;
                }
            }
            2 => y_position += length,
            3 => {
                for _ in 0..length {
                    let index = absolute_base + x_position - flag_offset;
                    let marked = index >= 0
                        && (index as usize) < flag.len()
                        && flag[index as usize] > snp_threshold as u8;
                    if x_position >= x_begin && x_position <= x_end && marked {
                        push(hap, x_total_start + x_position, y_total_start + y_position, 2, b'N');
                    }
                    x_position += 1;
                }
            }
            _ => {}
        }
        let _ = marked;
    }
}

/// `addSNPtohaplotype_advance`: mirrors [`mark_snp_advance`] but records the
/// evidence instead of counting mismatches.
#[allow(clippy::too_many_arguments)]
fn add_snp_to_haplotype_advance(
    window_start: i64,
    overlap_identifier: u32,
    x_total_start: i64,
    x_length: i64,
    y_total_start: i64,
    _y_length: i64,
    current: &WindowCigar,
    current_window: usize,
    beg_cigar: Option<(&WindowCigar, usize)>,
    end_cigar: Option<(&WindowCigar, usize)>,
    hap: &mut HaplotypeEvidenceAlloc,
    x_string: &[u8],
    y_string: &[u8],
    snp_threshold: i64,
) {
    let x_total_end = x_total_start + x_length - 1;
    let inner_offset = x_total_start - window_start;
    let mut current_begin = 0i64;
    let mut current_end = x_length - 1;
    // The reference reads `hap->flag` while it appends to `hap->list`; the
    // flag array is small, so a local copy keeps the borrows apart.
    let flag = hap.flag;

    if let Some((cigar, index)) = beg_cigar {
        let window = &cigar.windows[index];
        if window.y_end != -1 {
            let left_useless = i64::from(window.extra_begin);
            let right_useless = i64::from(window.extra_end);
            let x_left = x_total_start - i64::from(window.x_start);
            let x_right = i64::from(window.x_end) - x_total_start + 1;
            if x_left > left_useless && x_right > right_useless {
                let interval_begin = x_left;
                let interval_end = interval_begin + (x_right - right_useless) - 1;
                current_begin = x_right - right_useless;
                let x_offset = window.x_start.max(0) as usize;
                let y_offset = window.y_start.max(0) as usize;
                add_snp_to_haplotype_details(
                    cigar,
                    index,
                    &flag,
                    inner_offset,
                    &x_string[x_offset.min(x_string.len())..],
                    &y_string[y_offset.min(y_string.len())..],
                    i64::from(window.x_start),
                    i64::from(window.y_start),
                    interval_begin,
                    interval_end,
                    overlap_identifier,
                    interval_begin,
                    hap,
                    snp_threshold,
                );
            }
        }
    }

    if let Some((cigar, index)) = end_cigar {
        let window = &cigar.windows[index];
        if window.y_end != -1 {
            let left_useless = i64::from(window.extra_begin);
            let right_useless = i64::from(window.extra_end);
            let x_left = x_total_end - i64::from(window.x_start);
            let x_right = i64::from(window.x_end) - x_total_end + 1;
            if x_left > left_useless && x_right > right_useless {
                let interval_end = x_left;
                let interval_begin = interval_end - ((x_left + 1) - left_useless) + 1;
                current_end = (x_length - 1) - ((x_left + 1) - left_useless);
                let x_offset = window.x_start.max(0) as usize;
                let y_offset = window.y_start.max(0) as usize;
                add_snp_to_haplotype_details(
                    cigar,
                    index,
                    &flag,
                    i64::from(window.x_start) - window_start,
                    &x_string[x_offset.min(x_string.len())..],
                    &y_string[y_offset.min(y_string.len())..],
                    i64::from(window.x_start),
                    i64::from(window.y_start),
                    interval_begin,
                    interval_end,
                    overlap_identifier,
                    0,
                    hap,
                    snp_threshold,
                );
            }
        }
    }

    let x_offset = x_total_start.max(0) as usize;
    let y_offset = y_total_start.max(0) as usize;
    add_snp_to_haplotype_details(
        current,
        current_window,
        &flag,
        inner_offset,
        &x_string[x_offset.min(x_string.len())..],
        &y_string[y_offset.min(y_string.len())..],
        x_total_start,
        y_total_start,
        current_begin,
        current_end,
        overlap_identifier,
        0,
        hap,
        snp_threshold,
    );
}

/// `cluster_advance`: marks the mismatches of every overlap of the current
/// window and records the haplotype evidence.
fn cluster_advance(
    query: &[u8],
    window_start: i64,
    window_end: i64,
    overlaps: &[OverlapRegion],
    buffer: &IntervalBuffer,
    reads: &ReadStore,
    hap: &mut HaplotypeEvidenceAlloc,
    snp_threshold: i64,
) {
    let useful_length = window_end - window_start + 1;
    let mut first_snp = -1i64;
    let mut last_snp = -1i64;
    let mut informative = 0i64;
    let evidence_before = hap.length;

    // Pass 1: mark mismatches.
    for &overlap_identifier in buffer.full_windows() {
        let overlap = &overlaps[overlap_identifier as usize];
        let correct_x_start = (i64::from(overlap.x_position_start) / WINDOW_LENGTH) * WINDOW_LENGTH;
        let window_id = (window_start - correct_x_start) / WINDOW_LENGTH;
        if window_id < 0 || window_id as usize >= overlap.window_cigar.windows.len() {
            continue;
        }
        let window = overlap.window_cigar.windows[window_id as usize];
        if window.y_end == -1 {
            continue;
        }
        let x_start = i64::from(window.x_start);
        let x_length = i64::from(window.x_end) + 1 - x_start;
        let y_start = i64::from(window.y_start);
        let y_length = i64::from(window.y_end) + 1 - y_start;

        let beg = if window_id >= 1 {
            Some((&overlap.boundary_cigars, (window_id - 1) as usize))
        } else {
            None
        };
        let end = if (window_id as usize) + 1 < overlap.window_cigar.windows.len() {
            Some((&overlap.boundary_cigars, window_id as usize))
        } else {
            None
        };
        mark_snp_advance(
            window_start,
            x_start,
            x_length,
            y_start,
            y_length,
            &overlap.window_cigar,
            window_id as usize,
            beg,
            end,
            hap,
        );
    }

    for index in 0..useful_length {
        let index = index as usize;
        if index >= hap.flag.len() {
            break;
        }
        if hap.flag[index] != 0 {
            last_snp = index as i64;
            if first_snp == -1 {
                first_snp = index as i64;
            }
        }
        if i64::from(hap.flag[index]) > snp_threshold {
            hap.nn_snp += 1;
            informative += 1;
        }
    }
    if first_snp == -1 || last_snp == -1 {
        first_snp = 0;
        last_snp = -1;
    }

    // Pass 2: record the evidence of the marked positions.
    for &overlap_identifier in buffer.full_windows() {
        let overlap = &overlaps[overlap_identifier as usize];
        let correct_x_start = (i64::from(overlap.x_position_start) / WINDOW_LENGTH) * WINDOW_LENGTH;
        let window_id = (window_start - correct_x_start) / WINDOW_LENGTH;
        if window_id < 0 || window_id as usize >= overlap.window_cigar.windows.len() {
            continue;
        }
        let window = overlap.window_cigar.windows[window_id as usize];
        if window.y_end == -1 {
            continue;
        }
        let x_start = i64::from(window.x_start);
        let x_length = i64::from(window.x_end) + 1 - x_start;
        let y_start = i64::from(window.y_start);
        let y_length = i64::from(window.y_end) + 1 - y_start;

        let target_length = reads.read_length(overlap.y_identifier as usize);
        let target = reads.read_sub_region(
            overlap.y_identifier as usize,
            0,
            target_length as i64,
            overlap.y_strand as u8,
        );

        let beg = if window_id >= 1 {
            Some((&overlap.boundary_cigars, (window_id - 1) as usize))
        } else {
            None
        };
        let end = if (window_id as usize) + 1 < overlap.window_cigar.windows.len() {
            Some((&overlap.boundary_cigars, window_id as usize))
        } else {
            None
        };
        add_snp_to_haplotype_advance(
            window_start,
            overlap_identifier as u32,
            x_start,
            x_length,
            y_start,
            y_length,
            &overlap.window_cigar,
            window_id as usize,
            beg,
            end,
            hap,
            query,
            &target,
            snp_threshold,
        );
    }

    // `RsetInitHaplotypeEvdienceFlag`.
    let begin = first_snp.max(0) as usize;
    let end = (last_snp + 1).max(0) as usize;
    if begin < hap.flag.len() {
        let end = end.min(hap.flag.len());
        for value in &mut hap.flag[begin..end] {
            *value = 0;
        }
    }

    if hap.length.saturating_sub(evidence_before) > 1 && informative > 1 {
        radix_sort_by_key_32(&mut hap.list[evidence_before..hap.length], &|record| {
            record.site
        });
    }
}

/// `insert_snp_ee`: folds the evidence of one query site into allele
/// statistics, keeping only alleles seen at least twice.
fn insert_snp_ee(
    hap: &mut HaplotypeEvidenceAlloc,
    start: usize,
    count: usize,
    output: usize,
    query: &[u8],
) -> usize {
    let mut occurrence_0 = 0u64;
    let mut occurrence_1 = [0u64; 6];
    let mut occurrence_2 = 0u64;
    let mut difference = 0u64;

    for index in start..start + count {
        let evidence = hap.list[index];
        let cov = u64::from(evidence.coverage);
        match evidence.evidence_type {
            0 => occurrence_0 += cov,
            1 => {
                let base = usize::from(SEQUENCE_NT6_TABLE[evidence.mismatch_base as usize]);
                if base < 6 {
                    occurrence_1[base] += cov;
                }
                difference += cov;
            }
            _ => {}
        }
        occurrence_2 += cov;
    }

    if occurrence_0 == 0 || difference <= 1 {
        return 0;
    }

    let mut allele_ids = [u64::MAX; 6];
    let mut statistics = 0usize;
    let mut homopolymer: Option<u32> = None;
    for base in 0..4 {
        if occurrence_1[base] >= 2 {
            let site = hap.list[start].site;
            let is_homopolymer = *homopolymer.get_or_insert_with(|| {
                if is_strict_homopolymer(site as i64, query) {
                    1
                } else {
                    0
                }
            });
            hap.snp_stat.push(SnpStats {
                identifier: hap.snp_stat.len() as u32,
                occ_0: (1 + occurrence_0) as u32,
                occ_1: occurrence_1[base] as u32,
                occ_2: (occurrence_2.wrapping_sub(1 + occurrence_0).wrapping_sub(occurrence_1[base]))
                    as u32,
                overlap_num: occurrence_2 as u32,
                homopolymer_num: 0,
                non_homopolymer_num: 0,
                score: -1,
                site,
                is_homopolymer: is_homopolymer as u8,
            });
            allele_ids[base] = hap.snp_stat.len() as u64 - 1;
            statistics += 1;
        }
    }
    if statistics == 0 {
        return 0;
    }

    let mut written = 0usize;
    for index in start..start + count {
        let mut evidence = hap.list[index];
        if evidence.evidence_type == 0 {
            evidence.overlap_site = hap.snp_stat.len() as u32 - 1;
        } else {
            let base = usize::from(SEQUENCE_NT6_TABLE[evidence.mismatch_base as usize]);
            if base >= 6 || allele_ids[base] == u64::MAX {
                continue;
            }
            evidence.coverage = evidence.overlap_site;
            evidence.overlap_site = allele_ids[base] as u32;
        }
        hap.list[output + written] = evidence;
        written += 1;
    }
    written
}

/// `partition_overlaps_advance`: the window loop that turns verified overlaps
/// into haplotype evidence and re-labels large indels.
pub fn partition_overlaps_advance(
    overlaps: &mut [OverlapRegion],
    reads: &ReadStore,
    query: &[u8],
    hap: &mut HaplotypeEvidenceAlloc,
    maximum_overlap_divergence: f64,
    informative_coverage: i64,
    single_haplotype_coverage: i64,
) {
    hap.resize_init();

    let mut pool = WindowPool::new(
        query.len() as i64,
        WINDOW_LENGTH,
        (1.0 / maximum_overlap_divergence) as i64,
    );
    let mut buffer = IntervalBuffer::with_capacity(overlaps.len());
    let mut flag = 0i32;
    while flag != INTERVAL_PAST_THE_END {
        let Some((window_start, window_end)) = pool.next_window() else {
            break;
        };
        buffer.length = 0;
        buffer.partial_length = 0;
        flag = get_available_interval(window_start, window_end, overlaps, &mut buffer);
        cluster_advance(
            query,
            window_start,
            window_end,
            overlaps,
            &buffer,
            reads,
            hap,
            SNP_THRESHOLD,
        );
    }

    hap.set_snp_matrix(hap.nn_snp, overlaps.len() as u64);

    let mut k = 1usize;
    let mut l = 0usize;
    let mut output = 0usize;
    while k <= hap.length {
        if k == hap.length || hap.list[k].site != hap.list[l].site {
            output += insert_snp_ee(hap, l, k - l, output, query);
            l = k;
        }
        k += 1;
    }
    hap.length = output;

    let mut mapped_overlaps_length = 0i64;
    generate_haplotypes_naive_hifi(
        hap,
        overlaps,
        0.04,
        1,
        0.0,
        u64::MAX,
        informative_coverage,
        single_haplotype_coverage,
        &mut mapped_overlaps_length,
    );

    lable_large_indels(overlaps, query.len() as i64, maximum_overlap_divergence);
}

/// `lable_large_indels`: marks overlaps with a long insertion/deletion or an
/// unaligned middle window, then re-labels the ones inside a read-level long
/// indel.
pub fn lable_large_indels(
    overlaps: &mut [OverlapRegion],
    read_length: i64,
    maximum_overlap_divergence: f64,
) {
    let mut is_delete = false;
    for overlap in overlaps.iter_mut() {
        if overlap.is_match != 1 || overlap.window_cigar.windows.len() < 3 {
            continue;
        }
        let window_count = overlap.window_cigar.windows.len();
        let mut marked = false;
        for window_index in 1..window_count - 1 {
            let window = overlap.window_cigar.windows[window_index];
            if window.y_end == -1 {
                overlap.is_match = 100;
                marked = true;
                break;
            }
            if window.cigar_length < 3 {
                continue;
            }
            for offset in 1..window.cigar_length - 1 {
                let (operation, length) = overlap.window_cigar.cell(window_index, offset);
                if length <= 5 {
                    continue;
                }
                if operation == 2 || operation == 3 {
                    overlap.is_match = 100;
                    marked = true;
                    break;
                }
            }
            if marked {
                break;
            }
        }
        if marked {
            is_delete = true;
        }
    }

    if is_delete {
        let mut pool = WindowPool::new(
            read_length,
            WINDOW_LENGTH,
            (1.0 / maximum_overlap_divergence) as i64,
        );
        let mut start_index = 0u64;
        let mut real_length = 0i64;
        let mut real_length_100 = 0i64;
        let mut to_recover = false;
        let mut flag = 0i32;
        while flag != -2 {
            let Some((window_start, window_end)) = pool.next_window() else {
                break;
            };
            flag = available_fully_covered_interval(
                window_start,
                window_end,
                overlaps,
                &mut start_index,
                &mut real_length,
                &mut real_length_100,
            );
            if real_length == 0 && real_length_100 > 0 {
                to_recover = true;
                break;
            }
        }
        if to_recover {
            for overlap in overlaps.iter_mut() {
                if overlap.is_match == 100 {
                    overlap.is_match = 1;
                }
            }
        }
    }

    for overlap in overlaps.iter_mut() {
        if overlap.is_match == 1 {
            overlap.without_large_indel = 1;
        }
        if overlap.is_match == 100 {
            overlap.is_match = 1;
            overlap.without_large_indel = 0;
        }
    }
}

//! Haplotype phasing of one corrected read, the counterpart of `rphase_hc`,
//! `hc_phase_robust_rr`, `extract_sub_cigar_hc` and `hpc_mask_ff` in
//! `Correct.cpp`.
//!
//! The query read is walked in sliding windows. Every overlap candidate that
//! intersects a window is aligned base by base against the query: the matching
//! cells of its stored window CIGAR contribute the query base, the mismatch
//! cells contribute the target base, and a per-base counter records how many
//! overlaps disagree with the query at that position. The counter is folded
//! into the site statistics every fixed number of windows, and the resulting
//! allele counts decide which overlaps support the second haplotype.
//!
//! Only the HiFi configuration is covered: the boundary trimming (`bd`) and the
//! homopolymer length (`hpc_len`) are carried through, but the indel pass
//! (`rphase_lidel`) is not ported yet and is rejected explicitly.

use crate::container::ksort::{radix_sort_by_key_32, radix_sort_u64};
use crate::correction::haplotype_evidence::{
    generate_haplotypes_naive_hifi, generate_haplotypes_sv, push_info, HaplotypeEvidenceAlloc,
};
use crate::correction::window::WindowList;
use crate::correction::window_chain::is_unaligned_window;
use crate::overlap::region::OverlapRegion;
use crate::reads::read_store::ReadStore;

/// `HPC_PL`: homopolymer flank used by the caller.
pub const HPC_FLANK: i64 = 12;
/// `HPC_RR`: homopolymer rounds.
pub const HPC_ROUNDS: i64 = 4;
/// `HPC_CC`: homopolymer cutoff.
pub const HPC_CUTOFF: i64 = 2;

/// `SEC_MODE`: mask of the `sec` field of `ul_ov_t`.
pub const SEC_MODE: u32 = 0x3fff_ffff;

/// `ul_ov_t`: one window of one overlap while it is phased.
///
/// The reference packs the strand flag and the boundary into a 32-bit
/// bit-field; here they are plain fields. `boundary` keeps only the low 30
/// bits, exactly like `ovlp_bd`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct UlOverlap {
    /// Target overlap id, `tn`, an index into the overlap list.
    pub overlap_identifier: u32,
    /// Index of the current window, `qn`.
    pub current_window: u32,
    /// Query offset, `qs`.
    pub current_x_offset: u32,
    /// Target offset, `ts`.
    pub current_y_offset: u32,
    /// Length of the target region spanned by the first pass, `te`.
    pub current_y_length: u32,
    /// Cigar offset, `qe`.
    pub current_cigar_offset: u32,
    /// Boundary trimming, `sec` masked with [`SEC_MODE`].
    pub boundary: u32,
}

/// Scratch buffers of [`rphase_hc`], the reference's `c_idx` and `idx`.
#[derive(Debug, Clone, Default)]
pub struct RphaseScratch {
    /// Candidate windows, `kv_ul_ov_t c_idx`.
    pub candidates: Vec<UlOverlap>,
    /// Sort keys; the low 32 bits address [`Self::candidates`].
    pub index: Vec<u64>,
    /// `b->clist.chainDP`: the ONT DP workspace.
    pub dp: crate::overlap::region::ChainData,
    /// `idx`/`buf` of `gen_rphase_dp`.
    pub dp_index: Vec<u64>,
    pub dp_result: Vec<u64>,
    /// `c_idx` reused as the variant vector of `rphase_lidel`.
    pub lidel_variants: Vec<crate::assembly::ul::UlOv>,
}

/// `hpc_mask_ff` with a NULL flag buffer, the only form the HiFi path calls.
///
/// Returns one when `position` sits in a homopolymer run of at least
/// `hpc_cutoff * round` bases for some round. The reference clears flag bytes
/// when a buffer is supplied; every call on this path passes NULL, so that
/// branch is omitted.
pub fn hpc_mask_ff(
    sequence: &[u8],
    position: i64,
    hpc_flank: i64,
    hpc_rounds: i64,
    hpc_cutoff: i64,
) -> bool {
    let length = sequence.len() as i64;
    let start = if position >= hpc_flank {
        position - hpc_flank
    } else {
        0
    };
    let end = if position + hpc_flank <= length {
        position + hpc_flank
    } else {
        length
    };

    for round in 1..=hpc_rounds {
        let required = round * hpc_cutoff;

        // Including `position`.
        let mut k = position + round;
        while k < end && (k - round) >= start && sequence[k as usize] == sequence[(k - round) as usize]
        {
            k += 1;
        }
        let mut run_end = k.min(end);
        k = position - 1;
        while k >= start && (k + round) < end && sequence[k as usize] == sequence[(k + round) as usize]
        {
            k -= 1;
        }
        let mut run_start = (k + 1).max(start);
        if (run_end - run_start) > round && (run_end - run_start) >= required {
            return true;
        }

        // The first run does not include `position` on the left.
        k = position + round + 1;
        while k < end && (k - round) >= start && sequence[k as usize] == sequence[(k - round) as usize]
        {
            k += 1;
        }
        run_start = (position + 1).max(start);
        run_end = k.min(end);
        if (run_end - run_start) > round && (run_end - run_start) >= required {
            return true;
        }

        // Including `position` again, walking the left run first.
        k = position - round;
        while k >= start && (k + round) < end && sequence[k as usize] == sequence[(k + round) as usize]
        {
            k -= 1;
        }
        run_start = (k + 1).max(start);
        k = position + 1;
        while k < end && (k - round) >= start && sequence[k as usize] == sequence[(k - round) as usize]
        {
            k += 1;
        }
        run_end = k.min(end);
        if (run_end - run_start) > round && (run_end - run_start) >= required {
            return true;
        }

        // The last run stops before `position`.
        k = position - round - 1;
        while k >= start && (k + round) < end && sequence[k as usize] == sequence[(k + round) as usize]
        {
            k -= 1;
        }
        run_start = (k + 1).max(start);
        run_end = position.min(end);
        if (run_end - run_start) > round && (run_end - run_start) >= required {
            return true;
        }
    }

    false
}

/// `hpc_mask_ff_region` (`Correct.cpp` 10454): the interval form of the
/// homopolymer mask used by the ONT indel pass.
#[allow(clippy::too_many_arguments)]
pub fn hpc_mask_ff_region(
    sequence: &[u8],
    s0: i64,
    e0: i64,
    hpc_flank: i64,
    hpc_rounds: i64,
    hpc_cutoff: i64,
    hpc_rate: f64,
) -> bool {
    assert!(e0 >= s0);
    let sn = sequence.len() as i64;
    let s = if s0 >= hpc_flank { s0 - hpc_flank } else { 0 };
    let e = if e0 + hpc_flank <= sn { e0 + hpc_flank } else { sn };
    let base = if e0 > s0 { e0 - 1 } else { e0 };

    let overlap_ok = |zs: i64, ze: i64| -> bool {
        let os = zs.max(s0);
        let oe = ze.min(e0);
        let ovlp = if oe > os { oe - os } else { 0 };
        if e0 > s0 {
            ovlp > 0 && (ovlp as f64) >= ((e0 - s0) as f64 * hpc_rate)
        } else {
            zs <= s0 && ze >= e0
        }
    };

    for r in 1..=hpc_rounds {
        let rc = r * hpc_cutoff.max(hpc_rounds);

        // Including the right end.
        let mut k = base + r;
        while k < e && (k - r) >= s && sequence[k as usize] == sequence[(k - r) as usize] {
            k += 1;
        }
        let ze = k.min(e);
        let mut k = base - 1;
        while k >= s && (k + r) < e && sequence[k as usize] == sequence[(k + r) as usize] {
            k -= 1;
        }
        let zs = (k + 1).max(s);
        if (ze - zs) > r && (ze - zs) >= rc && overlap_ok(zs, ze) {
            return true;
        }

        // Excluding the point, right side.
        if e0 <= s0 {
            let mut k = e0 + r + 1;
            while k < e && (k - r) >= s && sequence[k as usize] == sequence[(k - r) as usize] {
                k += 1;
            }
            let zs = (e0 + 1).max(s);
            let ze = k.min(e);
            if (ze - zs) > r && (ze - zs) >= rc {
                return true;
            }
        }

        // Including the left end.
        let mut k = s0 - r;
        while k >= s && (k + r) < e && sequence[k as usize] == sequence[(k + r) as usize] {
            k -= 1;
        }
        let zs = (k + 1).max(s);
        let mut k = s0 + 1;
        while k < e && (k - r) >= s && sequence[k as usize] == sequence[(k - r) as usize] {
            k += 1;
        }
        let ze = k.min(e);
        if (ze - zs) > r && (ze - zs) >= rc && overlap_ok(zs, ze) {
            return true;
        }

        // Excluding the point, left side.
        if e0 <= s0 {
            let mut k = s0 - r - 1;
            while k >= s && (k + r) < e && sequence[k as usize] == sequence[(k + r) as usize] {
                k -= 1;
            }
            let zs = (k + 1).max(s);
            let ze = s0.min(e);
            if (ze - zs) > r && (ze - zs) >= rc {
                return true;
            }
        }
    }
    false
}

/// `extract_sub_cigar_hc`: walks the CIGAR of one window over `[s, e)`.
///
/// With `set_flag != 0` the pass counts mismatches per query base into
/// [`HaplotypeEvidenceAlloc::flag`] and records where the target walk ended;
/// with `set_flag == 0` it turns the flagged positions into evidence records,
/// recovering target bases from the read on demand.
///
/// `flag_base` is the reference's `f = hp->flag + os - window_start`: the
/// offset of this candidate's slice inside the window counter array.
#[allow(clippy::too_many_arguments)]
pub fn extract_sub_cigar_hc(
    overlap: &OverlapRegion,
    reads: &ReadStore,
    hap: &mut HaplotypeEvidenceAlloc,
    query: &[u8],
    s_input: i64,
    e_input: i64,
    cursor: &mut UlOverlap,
    set_flag: i64,
    flag_base: i64,
    _occurrence_threshold: u8,
    hpc_length: u64,
) -> i64 {
    let window_index = cursor.current_window as usize;
    let boundary = cursor.boundary as i64;
    let window = overlap.window_cigar.windows[window_index];

    let mut s = s_input;
    let mut e = e_input;
    let window_start = window.x_start as i64 + boundary;
    let window_end = window.x_end as i64 + 1 - boundary;
    if s < window_start {
        s = window_start;
    }
    if e > window_end {
        e = window_end;
    }
    if s >= e {
        return -1;
    }
    let os = s.max(window_start);
    let oe = e.min(window_end);
    if oe <= os {
        return -1;
    }

    let cigar_length = window.cigar_length as i64;
    if cigar_length == 0 {
        return -1;
    }

    let mut xk = cursor.current_x_offset as i64;
    let mut yk = cursor.current_y_offset as i64;
    let mut ck = cursor.current_cigar_offset as i64;
    let assertion_x = xk;
    let assertion_y = yk;
    let assertion_cigar = ck;
    if ck < 0 || ck > cigar_length {
        ck = 0;
        xk = window.x_start as i64;
        yk = window.y_start as i64;
    }

    let cell = |index: i64| -> u16 {
        overlap.window_cigar.cells[(window.cigar_start as i64 + index) as usize]
    };

    // Rewind to the last cell that still starts at or after `s`.
    while ck > 0 && xk > s {
        ck -= 1;
        let operation = (cell(ck) >> 14) as u8;
        let length = (cell(ck) & 0x3fff) as i64;
        if operation != 2 {
            xk -= length;
        }
        if operation != 3 {
            yk -= length;
        }
    }

    // The cursor at the interleaved position, the reference's `xk0`/`yk0`/
    // `ck0` after the `set_f` branch saves them again.
    let cursor_x = xk;
    let cursor_y = yk;
    let cursor_cigar = ck;

    // In the second pass `yk1` is the end of the target region recovered by
    // the first pass; the start is recomputed per mismatch below.
    let mut target_end: i64 = -1;
    let mut walked_target_end: i64 = -1;
    if set_flag != 0 {
        cursor.current_y_length = 0;
    } else {
        assert_eq!(
            assertion_x, cursor_x,
            "the cursor must match between the two passes"
        );
        assert_eq!(
            assertion_y, cursor_y,
            "the cursor must match between the two passes"
        );
        assert_eq!(
            assertion_cigar, cursor_cigar,
            "the cursor must match between the two passes"
        );
        target_end = cursor_y + cursor.current_y_length as i64;
        let _read_length = reads.read_length(overlap.y_identifier as usize);
    }

    let mut target: Option<Vec<u8>> = None;
    // The reference's `yk0` inside the second pass: the read position the
    // recovered buffer starts at, fixed by the first mismatch.
    let mut target_start: i64 = 0;

    while ck < cigar_length && xk < e {
        let window_x = xk;
        let packed = cell(ck);
        let operation = (packed >> 14) as u8;
        let length = (packed & 0x3fff) as i64;
        if operation != 2 {
            xk += length;
        }
        if operation != 3 {
            yk += length;
        }
        ck += 1;
        let window_x_end = xk;
        if operation != 0 && operation != 1 {
            continue;
        }

        let os = s.max(window_x);
        let oe = e.min(window_x_end);
        let overlap_length = if oe > os { oe - os } else { 0 };
        if overlap_length == 0 {
            continue;
        }

        if set_flag != 0 {
            if operation == 1 {
                for t in os..oe {
                    let index = (flag_base + (t - s)) as usize;
                    hap.flag[index] = if hap.flag[index] <= 126 {
                        hap.flag[index] + 1
                    } else {
                        127
                    };
                }
                walked_target_end = oe - xk + yk;
            }
        } else if operation == 0 {
            for t in os..oe {
                let index = (flag_base + (t - s)) as usize;
                if hap.flag[index] != 0 {
                    let is_homopolymer = hap.flag[index] == 3;
                    hap.add(crate::correction::haplotype_evidence::HaplotypeEvidence {
                        site: t as u32,
                        overlap_identifier: cursor.overlap_identifier,
                        overlap_site: (t - xk + yk) as u32,
                        evidence_type: (is_homopolymer as u8) << 1,
                        coverage: 1,
                        mismatch_base: query[t as usize],
                    });
                }
            }
        } else if operation == 1 {
            for t in os..oe {
                let index = (flag_base + (t - s)) as usize;
                if hap.flag[index] == 0 {
                    continue;
                }
                if target.is_none() {
                    target_start = (t - xk + yk - hpc_length as i64).max(0);
                    let length = target_end - target_start;
                    target = Some(reads.read_sub_region(
                        overlap.y_identifier as usize,
                        target_start,
                        length,
                        overlap.y_strand as u8,
                    ));
                }
                let mut is_homopolymer = hap.flag[index] == 3;
                if !is_homopolymer
                    && hpc_length != 0
                    && hpc_mask_ff(
                        target.as_ref().expect("the target must be recovered"),
                        t - xk + yk - target_start,
                        hpc_length as i64,
                        HPC_ROUNDS,
                        HPC_CUTOFF,
                    )
                {
                    is_homopolymer = true;
                }
                hap.add(crate::correction::haplotype_evidence::HaplotypeEvidence {
                    site: t as u32,
                    overlap_identifier: cursor.overlap_identifier,
                    overlap_site: (t - xk + yk) as u32,
                    evidence_type: (is_homopolymer as u8) << 1 | 1,
                    coverage: 1,
                    mismatch_base: target.as_ref().expect("the target must be recovered")
                        [(t - xk + yk - target_start) as usize],
                });
            }
        }
    }

    if set_flag != 0 {
        cursor.current_x_offset = cursor_x as u32;
        cursor.current_y_offset = cursor_y as u32;
        cursor.current_cigar_offset = cursor_cigar as u32;
        if walked_target_end != -1 {
            let mut end = walked_target_end + hpc_length as i64;
            let read_length = reads.read_length(overlap.y_identifier as usize) as i64;
            if end > read_length {
                end = read_length;
            }
            cursor.current_y_length = (end - cursor_y) as u32;
        } else {
            cursor.current_y_length = 0;
        }
    } else {
        cursor.current_x_offset = xk as u32;
        cursor.current_y_offset = yk as u32;
        cursor.current_cigar_offset = ck as u32;
        cursor.current_y_length = 0;
    }

    1
}

/// `hc_phase_robust_rr`: runs one pass over the candidates of a window.
///
/// Returns one when some candidate ends at or before `window_end`, which tells
/// the caller that the candidate list can be compacted.
#[allow(clippy::too_many_arguments)]
pub fn hc_phase_robust_rr(
    overlaps: &[OverlapRegion],
    reads: &ReadStore,
    hap: &mut HaplotypeEvidenceAlloc,
    query: &[u8],
    candidates: &mut [UlOverlap],
    candidate_ids: &[u64],
    window_start: i64,
    window_end: i64,
    set_flag: i64,
    occurrence_threshold: u8,
    hpc_length: u64,
) -> u64 {
    let mut round_reached = 0u64;
    for identifier in candidate_ids {
        let candidate_index = *identifier as u32 as usize;
        let mut cursor = candidates[candidate_index];
        let overlap = &overlaps[cursor.overlap_identifier as usize];
        let window = overlap.window_cigar.windows[cursor.current_window as usize];
        let boundary = cursor.boundary as i64;
        let q0 = window.x_start as i64 + boundary;
        let q1 = window.x_end as i64 + 1 - boundary;
        if q1 <= window_end {
            round_reached = 1;
        }
        let os = q0.max(window_start);
        let oe = q1.min(window_end);
        if oe > os {
            extract_sub_cigar_hc(
                overlap,
                reads,
                hap,
                query,
                os,
                oe,
                &mut cursor,
                set_flag,
                os - window_start,
                occurrence_threshold,
                hpc_length,
            );
        }
        candidates[candidate_index] = cursor;
    }
    round_reached
}

/// `rphase_hc`: phases one corrected read against its overlaps.
///
/// `boundary` trims both ends of every window (zero on the HiFi path),
/// `window_length` is the sliding window size (`WINDOW_MAX_SIZE`), and
/// `occurrence_threshold` is the minimum number of disagreeing overlaps for a
/// position to count as a site. `standard_base` selects the stand-alone test of
/// the second correction round; `lindel` enables the indel pass, which is not
/// ported.
#[allow(clippy::too_many_arguments)]
pub fn rphase_hc(
    overlaps: &mut [OverlapRegion],
    reads: &ReadStore,
    hap: &mut HaplotypeEvidenceAlloc,
    query: &[u8],
    scratch: &mut RphaseScratch,
    boundary: i64,
    window_length: i64,
    occurrence_threshold: u8,
    hpc_length: u64,
    standard_base: u64,
    lindel: u8,
    mapped_overlaps_length: &mut i64,
    informative_coverage: i64,
    single_haplotype_coverage: i64,
    query_id: u64,
    het_cov: i64,
    hom_cov: i64,
    n_hap: i64,
    quality_reads: Option<&crate::reads::read_store::ReadStore>,
    use_dp: bool,
) {
    let query_length = query.len() as i64;

    // Collect the candidate windows, keyed by their query start.
    scratch.candidates.clear();
    scratch.index.clear();
    for (overlap_index, overlap) in overlaps.iter().enumerate() {
        let window_count = overlap.window_cigar.windows.len();
        if window_count == 0 {
            continue;
        }
        for window_index in 0..window_count {
            let window = overlap.window_cigar.windows[window_index];
            if is_unaligned_window(&window) {
                continue;
            }
            let q0 = window.x_start as i64 + boundary;
            let q1 = window.x_end as i64 - boundary;
            if q1 >= q0 {
                let key = ((q0 as u64) << 32).wrapping_add(scratch.candidates.len() as u64);
                scratch.index.push(key);
                scratch.candidates.push(UlOverlap {
                    overlap_identifier: overlap_index as u32,
                    current_window: window_index as u32,
                    current_x_offset: window.x_start as u32,
                    current_y_offset: window.y_start as u32,
                    current_y_length: 0,
                    current_cigar_offset: 0,
                    boundary: (boundary as u32) & SEC_MODE,
                });
            }
        }
    }

    // Sort by query start, then re-key every group of equal starts by query
    // end so that the candidates of one start are visited left to right.
    let sorted_count = scratch.index.len();
    radix_sort_u64(&mut scratch.index);
    let mut group_start = 0usize;
    for k in 1..=sorted_count {
        if k == sorted_count || (scratch.index[k] >> 32) != (scratch.index[group_start] >> 32) {
            if k - group_start > 1 {
                for t in group_start..k {
                    let candidate_index = scratch.index[t] as u32 as usize;
                    let cursor = scratch.candidates[candidate_index];
                    let window = overlaps[cursor.overlap_identifier as usize]
                        .window_cigar
                        .windows[cursor.current_window as usize];
                    let q1 = window.x_end as i64 + 1 - boundary;
                    scratch.index[t] =
                        ((q1 as u64) << 32).wrapping_add(scratch.index[t] as u32 as u64);
                }
                radix_sort_u64(&mut scratch.index[group_start..k]);
            }
            group_start = k;
        }
    }

    hap.resize_init();

    let mut round_reached: i64 = 0;
    let mut window_start: i64 = 0;
    let mut window_end: i64 = window_length.min(query_length);
    let mut candidate_cursor = 0usize;
    while window_start < query_length {
        // Drop the candidates that no longer reach the current window.
        if round_reached != 0 {
            let mut destination = sorted_count;
            for m in sorted_count..scratch.index.len() {
                let candidate_index = scratch.index[m] as u32 as usize;
                let cursor = scratch.candidates[candidate_index];
                let window = overlaps[cursor.overlap_identifier as usize]
                    .window_cigar
                    .windows[cursor.current_window as usize];
                let q0 = window.x_start as i64 + boundary;
                let q1 = window.x_end as i64 + 1 - boundary;
                let os = q0.max(window_start);
                let oe = q1.min(window_end);
                if oe > os {
                    scratch.index[destination] = scratch.index[m];
                    destination += 1;
                }
            }
            scratch.index.truncate(destination);
        }

        // Add the candidates that enter the window.
        while candidate_cursor < sorted_count {
            let key = scratch.index[candidate_cursor];
            let candidate_index = key as u32 as usize;
            let cursor = scratch.candidates[candidate_index];
            let window = overlaps[cursor.overlap_identifier as usize]
                .window_cigar
                .windows[cursor.current_window as usize];
            let q0 = window.x_start as i64 + boundary;
            if q0 >= window_end {
                break;
            }
            let q1 = window.x_end as i64 + 1 - boundary;
            let os = q0.max(window_start);
            let oe = q1.min(window_end);
            if oe > os {
                scratch.index.push(key as u32 as u64);
            }
            candidate_cursor += 1;
        }

        let evidence_before = hap.length;
        round_reached = hc_phase_robust_rr(
            overlaps,
            reads,
            hap,
            query,
            &mut scratch.candidates,
            &scratch.index[sorted_count..],
            window_start,
            window_end,
            1,
            occurrence_threshold,
            hpc_length,
        ) as i64;

        // Keep the positions that enough overlaps disagree with.
        let window_span = (window_end - window_start) as u64;
        let mut information_seen = 0u64;
        let mut last_informative_end = 0u64;
        let mut first_informative_start = u64::MAX;
        for offset in 0..window_span {
            let index = offset as usize;
            if hap.flag[index] == 0 {
                continue;
            }
            if hap.flag[index] > occurrence_threshold {
                information_seen = 1;
                hap.nn_snp += 1;
                hap.flag[index] = 1;
                if hpc_length != 0
                    && hpc_mask_ff(
                        query,
                        (offset + window_start as u64) as i64,
                        hpc_length as i64,
                        HPC_ROUNDS,
                        HPC_CUTOFF,
                    )
                {
                    hap.flag[index] = 3;
                }
                last_informative_end = offset + 1;
                if first_informative_start == u64::MAX {
                    first_informative_start = offset;
                }
            } else {
                hap.flag[index] = 0;
            }
        }

        if information_seen != 0 {
            hc_phase_robust_rr(
                overlaps,
                reads,
                hap,
                query,
                &mut scratch.candidates,
                &scratch.index[sorted_count..],
                window_start,
                window_end,
                0,
                occurrence_threshold,
                hpc_length,
            );
            if hap.length > evidence_before {
                radix_sort_by_key_32(&mut hap.list[evidence_before..hap.length], &|record| {
                    record.site
                });
            }
        }

        if last_informative_end > first_informative_start {
            for index in first_informative_start as usize..last_informative_end as usize {
                hap.flag[index] = 0;
            }
        }

        window_start += window_length;
        window_end += window_length;
        window_end = window_end.min(query_length);
    }

    // Fold the evidence of every site into allele statistics.
    hap.set_snp_matrix(hap.nn_snp, overlaps.len() as u64);
    let evidence_length = hap.length;
    // `z = (std_bs || t8) ? ol->list : NULL`.
    let statistics_source = if standard_base != 0 || quality_reads.is_some() {
        Some(&overlaps[..])
    } else {
        None
    };
    let mut out = 0usize;
    let mut group_start = 0usize;
    for k in 1..=evidence_length {
        if k == evidence_length || hap.list[k].site != hap.list[group_start].site {
            out += push_info(
                hap,
                group_start,
                k - group_start,
                out,
                statistics_source,
                quality_reads,
                5,
            );
            group_start = k;
        }
    }
    hap.length = out;

    // `if (dp) gen_rphase_dp(...)` runs before the naive promotion on ONT.
    if use_dp {
        crate::correction::rphase_dp::gen_rphase_dp(
            hap,
            overlaps,
            reads,
            query_id,
            if standard_base != 0 { 0.05 } else { 0.0 },
            if standard_base != 0 { 2 } else { u64::MAX },
            &mut scratch.dp,
            &mut scratch.dp_index,
            &mut scratch.dp_result,
            informative_coverage,
            single_haplotype_coverage,
            het_cov,
            hom_cov,
            n_hap,
        );
    }
    generate_haplotypes_naive_hifi(
        hap,
        overlaps,
        0.04,
        if standard_base != 0 { 0 } else { 1 },
        if standard_base != 0 { 0.05 } else { 0.0 },
        if standard_base != 0 { 2 } else { u64::MAX },
        informative_coverage,
        single_haplotype_coverage,
        mapped_overlaps_length,
    );

    if lindel != 0 {
        let found = crate::correction::rphase_lidel::rphase_lidel(
            overlaps,
            reads,
            hap,
            query,
            &mut scratch.lidel_variants,
            &mut scratch.dp_index,
            &mut scratch.dp_result,
            boundary,
            hpc_length,
            het_cov,
            hom_cov,
            n_hap,
        );
        if found != 0 {
            generate_haplotypes_sv(
                hap,
                overlaps,
                informative_coverage,
                single_haplotype_coverage,
                mapped_overlaps_length,
            );
        }
    }
}

/// The window of `cursor` inside `overlap`, used by the tests and by callers
/// that need the boundaries without borrowing the cigar.
pub fn cursor_window(overlap: &OverlapRegion, cursor: &UlOverlap) -> WindowList {
    overlap.window_cigar.windows[cursor.current_window as usize]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn homopolymer_masking_matches_the_reference_examples() {
        let homopolymer = b"AAAAAAAAAAAA";
        assert!(hpc_mask_ff(homopolymer, 5, 12, HPC_ROUNDS, HPC_CUTOFF));
        // The mask also fires on short-period repeats: at round two the
        // dinucleotide unit satisfies the `round * cutoff` length.
        let alternating = b"ACACACACACAC";
        assert!(hpc_mask_ff(alternating, 5, 12, HPC_ROUNDS, HPC_CUTOFF));
        // A two-base run is long enough for round one only.
        assert!(hpc_mask_ff(b"ACGTAAACGT", 5, 12, HPC_ROUNDS, HPC_CUTOFF));
        // A repeat that does not contain the position is invisible to the
        // position-centred scan.
        assert!(!hpc_mask_ff(b"TTATGCAGAAAA", 5, 12, HPC_ROUNDS, HPC_CUTOFF));
    }

    #[test]
    fn boundary_masking_keeps_the_low_30_bits() {
        let cursor = UlOverlap {
            boundary: 0xffff_ffff,
            ..UlOverlap::default()
        };
        assert_eq!(cursor.boundary & SEC_MODE, 0x3fff_ffff);
    }
}

#[cfg(test)]
mod probe_goldens {
    use super::*;
    use crate::correction::haplotype_evidence::WINDOW_MAXIMUM_SIZE;
    use crate::correction::window::WindowCigar;
    use crate::overlap::region::OverlapRegion;
    use crate::reads::read_store::ReadStore;

    struct CellSpec {
        operation: u8,
        length: u32,
    }

    struct WindowSpec {
        x_start: i32,
        x_end: i32,
        y_start: i32,
        y_end: i32,
        cells: Vec<CellSpec>,
    }

    struct OverlapSpec {
        y_identifier: u32,
        strand: u32,
        alignment_length: u32,
        windows: Vec<WindowSpec>,
    }

    /// `build_overlap`: an overlap with real window cigars, exactly as the C++
    /// probe constructs them.
    fn build_overlaps(specs: &[OverlapSpec]) -> Vec<OverlapRegion> {
        specs
            .iter()
            .map(|spec| {
                let mut window_cigar = WindowCigar::default();
                for window in &spec.windows {
                    let index = window_cigar.append_window(
                        window.x_start,
                        window.x_end,
                        window.y_start,
                        window.y_end,
                        0,
                        0,
                        0,
                        0,
                    );
                    let start = window_cigar.cells.len();
                    for cell in &window.cells {
                        window_cigar.push_cell(cell.operation, cell.length);
                    }
                    window_cigar.assign_cells(index, start);
                }
                let first = spec.windows.first().unwrap();
                let last = spec.windows.last().unwrap();
                OverlapRegion {
                    x_identifier: 0,
                    x_position_start: first.x_start as u32,
                    x_position_end: last.x_end as u32,
                    y_identifier: spec.y_identifier,
                    y_position_start: first.y_start as u32,
                    y_position_end: last.y_end as u32,
                    y_strand: spec.strand,
                    alignment_length: spec.alignment_length,
                    is_match: 1,
                    without_large_indel: 1,
                    strong: 0,
                    window_cigar,
                    ..OverlapRegion::default()
                }
            })
            .collect()
    }

    fn repeat(unit: &[u8], count: usize) -> Vec<u8> {
        let mut sequence = Vec::with_capacity(unit.len() * count);
        for _ in 0..count {
            sequence.extend_from_slice(unit);
        }
        sequence
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

    fn mutate(sequence: &mut [u8], position: usize, base: u8) {
        sequence[position] = base;
    }

    fn run_scenario(
        name: &str,
        read_sequences: &[Vec<u8>],
        query: &[u8],
        specs: &[OverlapSpec],
    ) -> Vec<String> {
        let mut reads = ReadStore::default();
        for (index, sequence) in read_sequences.iter().enumerate() {
            reads.push_read(&format!("read{index}"), sequence);
        }
        let mut overlaps = build_overlaps(specs);
        let mut hap = HaplotypeEvidenceAlloc::new();
        let mut scratch = RphaseScratch::default();
        let mut mapped: i64 = 0;

        rphase_hc(
            &mut overlaps,
            &reads,
            &mut hap,
            query,
            &mut scratch,
            0,
            WINDOW_MAXIMUM_SIZE as i64,
            1,
            0,
            0,
            0,
            &mut mapped,
            3,
            3,
            0,
            0,
            0,
            2,
            None,
            false,
        );

        let mut lines = Vec::new();
        lines.push(format!(
            "RPH name={} wmax={} length={} nn_snp={} statn={} overlap={} core={} snpsrt={}",
            name,
            WINDOW_MAXIMUM_SIZE,
            hap.length,
            hap.nn_snp,
            hap.snp_stat.len(),
            hap.overlap,
            hap.core_snp,
            hap.snp_srt.len()
        ));
        for (index, stats) in hap.snp_stat.iter().enumerate() {
            lines.push(format!(
                " SNP i={} id={} on={} occ0={} occ1={} occ2={} score={} site={} ishp={}",
                index,
                stats.identifier,
                stats.overlap_num,
                stats.occ_0,
                stats.occ_1,
                stats.occ_2,
                stats.score,
                stats.site,
                stats.is_homopolymer
            ));
        }
        for (index, evidence) in hap.list[..hap.length].iter().enumerate() {
            lines.push(format!(
                " HLE i={} site={} oid={} osite={} type={} cov={} base={}",
                index,
                evidence.site,
                evidence.overlap_identifier,
                evidence.overlap_site,
                evidence.evidence_type,
                evidence.coverage,
                evidence.mismatch_base
            ));
        }
        for (index, value) in hap.snp_srt.iter().enumerate() {
            lines.push(format!(" SSR i={} v={}", index, value));
        }
        lines.push(format!(" OLR mapped={} n={}", mapped, overlaps.len()));
        for (index, overlap) in overlaps.iter().enumerate() {
            lines.push(format!(
                " OLV i={} y={} strand={} aim={} strong={} wli={}",
                index,
                overlap.y_identifier,
                overlap.y_strand,
                overlap.is_match,
                overlap.strong,
                overlap.without_large_indel
            ));
        }
        lines
    }

    fn generate_lines() -> Vec<String> {
        let mut lines = Vec::new();

        // Scenario A.
        {
            let mut reads = vec![repeat(b"A", 60); 9];
            for read in reads.iter_mut().take(6).skip(3) {
                mutate(read, 20, b'C');
            }
            mutate(&mut reads[4], 40, b'G');
            mutate(&mut reads[5], 40, b'G');
            for read in reads.iter_mut().take(9).skip(6) {
                mutate(read, 40, b'T');
            }
            reads[8] = reverse_complement(&reads[8]);

            let mut specs = Vec::new();
            for i in 0..9u32 {
                let cells = if i < 3 {
                    vec![CellSpec { operation: 0, length: 60 }]
                } else if i == 3 {
                    vec![
                        CellSpec { operation: 0, length: 20 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 0, length: 39 },
                    ]
                } else if i <= 5 {
                    vec![
                        CellSpec { operation: 0, length: 20 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 0, length: 19 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 0, length: 19 },
                    ]
                } else {
                    vec![
                        CellSpec { operation: 0, length: 40 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 0, length: 19 },
                    ]
                };
                specs.push(OverlapSpec {
                    y_identifier: i,
                    strand: if i == 8 { 1 } else { 0 },
                    alignment_length: 60,
                    windows: vec![WindowSpec {
                        x_start: 0,
                        x_end: 59,
                        y_start: 0,
                        y_end: 59,
                        cells,
                    }],
                });
            }
            lines.extend(run_scenario(
                "strong_and_weak",
                &reads,
                &repeat(b"A", 60),
                &specs,
            ));
        }

        // Scenario B.
        {
            let mut reads = vec![repeat(b"A", 140); 7];
            for read in reads.iter_mut().skip(5) {
                mutate(read, 10, b'C');
                mutate(read, 50, b'G');
                mutate(read, 90, b'T');
                mutate(read, 130, b'G');
            }
            let mut specs = Vec::new();
            for i in 0..7u32 {
                let cells = if i < 5 {
                    vec![CellSpec { operation: 0, length: 140 }]
                } else {
                    vec![
                        CellSpec { operation: 0, length: 10 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 0, length: 39 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 0, length: 39 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 0, length: 39 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 0, length: 9 },
                    ]
                };
                specs.push(OverlapSpec {
                    y_identifier: i,
                    strand: 0,
                    alignment_length: 100,
                    windows: vec![WindowSpec {
                        x_start: 0,
                        x_end: 139,
                        y_start: 0,
                        y_end: 139,
                        cells,
                    }],
                });
            }
            lines.extend(run_scenario(
                "multi_check_weak",
                &reads,
                &repeat(b"A", 140),
                &specs,
            ));
        }

        // Scenario C.
        {
            let mut reads = vec![repeat(b"A", 500); 6];
            reads.push(repeat(b"A", 300));
            reads.push(repeat(b"A", 70));
            for read in reads.iter_mut().take(6).skip(3) {
                mutate(read, 100, b'C');
                mutate(read, 460, b'T');
            }
            let mut specs = Vec::new();
            for i in 0..8u32 {
                let mut alignment_length = 500;
                let window = if i < 3 {
                    WindowSpec {
                        x_start: 0,
                        x_end: 499,
                        y_start: 0,
                        y_end: 499,
                        cells: vec![CellSpec { operation: 0, length: 500 }],
                    }
                } else if i <= 5 {
                    WindowSpec {
                        x_start: 0,
                        x_end: 499,
                        y_start: 0,
                        y_end: 499,
                        cells: vec![
                            CellSpec { operation: 0, length: 100 },
                            CellSpec { operation: 1, length: 1 },
                            CellSpec { operation: 0, length: 359 },
                            CellSpec { operation: 1, length: 1 },
                            CellSpec { operation: 0, length: 39 },
                        ],
                    }
                } else if i == 6 {
                    alignment_length = 300;
                    WindowSpec {
                        x_start: 0,
                        x_end: 299,
                        y_start: 0,
                        y_end: 299,
                        cells: vec![CellSpec { operation: 0, length: 300 }],
                    }
                } else {
                    alignment_length = 70;
                    WindowSpec {
                        x_start: 430,
                        x_end: 499,
                        y_start: 0,
                        y_end: 69,
                        cells: vec![CellSpec { operation: 0, length: 70 }],
                    }
                };
                specs.push(OverlapSpec {
                    y_identifier: i,
                    strand: 0,
                    alignment_length,
                    windows: vec![window],
                });
            }
            lines.extend(run_scenario(
                "two_windows",
                &reads,
                &repeat(b"A", 500),
                &specs,
            ));
        }

        // Scenario D.
        {
            let mut indel = repeat(b"A", 28);
            indel.push(b'C');
            indel.extend(repeat(b"A", 16));
            indel.push(b'G');
            indel.extend(repeat(b"A", 14));
            let reads = vec![repeat(b"A", 60), repeat(b"A", 60), indel.clone(), indel];
            let mut specs = Vec::new();
            for i in 0..4u32 {
                let cells = if i < 2 {
                    vec![CellSpec { operation: 0, length: 60 }]
                } else {
                    vec![
                        CellSpec { operation: 0, length: 20 },
                        CellSpec { operation: 3, length: 2 },
                        CellSpec { operation: 0, length: 8 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 0, length: 9 },
                        CellSpec { operation: 2, length: 2 },
                        CellSpec { operation: 0, length: 5 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 0, length: 14 },
                    ]
                };
                specs.push(OverlapSpec {
                    y_identifier: i,
                    strand: 0,
                    alignment_length: 60,
                    windows: vec![WindowSpec {
                        x_start: 0,
                        x_end: 59,
                        y_start: 0,
                        y_end: 59,
                        cells,
                    }],
                });
            }
            lines.extend(run_scenario(
                "indel_shift",
                &reads,
                &repeat(b"A", 60),
                &specs,
            ));
        }

        // Scenario E.
        {
            let mut reads = vec![repeat(b"A", 60); 6];
            for read in reads.iter_mut().skip(4) {
                mutate(read, 20, b'C');
                mutate(read, 21, b'G');
                mutate(read, 40, b'T');
            }
            let mut specs = Vec::new();
            for i in 0..6u32 {
                let cells = if i < 4 {
                    vec![CellSpec { operation: 0, length: 60 }]
                } else {
                    vec![
                        CellSpec { operation: 0, length: 20 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 0, length: 18 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 0, length: 19 },
                    ]
                };
                specs.push(OverlapSpec {
                    y_identifier: i,
                    strand: 0,
                    alignment_length: 50,
                    windows: vec![WindowSpec {
                        x_start: 0,
                        x_end: 59,
                        y_start: 0,
                        y_end: 59,
                        cells,
                    }],
                });
            }
            lines.extend(run_scenario(
                "adjacent_sites",
                &reads,
                &repeat(b"A", 60),
                &specs,
            ));
        }

        // Scenario F.
        {
            let mut reads = vec![repeat(b"A", 80); 5];
            mutate(&mut reads[1], 20, b'C');
            for read in reads.iter_mut().skip(2) {
                mutate(read, 20, b'C');
                mutate(read, 60, b'T');
            }
            let mut specs = Vec::new();
            for i in 0..5u32 {
                let cells = if i == 0 {
                    vec![CellSpec { operation: 0, length: 80 }]
                } else if i == 1 {
                    vec![
                        CellSpec { operation: 0, length: 20 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 0, length: 59 },
                    ]
                } else {
                    vec![
                        CellSpec { operation: 0, length: 20 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 0, length: 39 },
                        CellSpec { operation: 1, length: 1 },
                        CellSpec { operation: 0, length: 19 },
                    ]
                };
                specs.push(OverlapSpec {
                    y_identifier: i,
                    strand: 0,
                    alignment_length: 50,
                    windows: vec![WindowSpec {
                        x_start: 0,
                        x_end: 79,
                        y_start: 0,
                        y_end: 79,
                        cells,
                    }],
                });
            }
            lines.extend(run_scenario(
                "strong_reset",
                &reads,
                &repeat(b"A", 80),
                &specs,
            ));
        }

        lines
    }

    #[test]
    #[ignore = "requires reference/ec_rphase_hc_probe.out from the C++ probe"]
    fn rphase_hc_matches_the_reference() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("reference")
            .join("ec_rphase_hc_probe.out");
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

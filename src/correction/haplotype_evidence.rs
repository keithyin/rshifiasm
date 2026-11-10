//! Haplotype evidence and per-site allele statistics, the counterpart of
//! `haplotype_evdience`, `SnpStats`, `haplotype_evdience_alloc`, `push_info`
//! and `generate_haplotypes_naive_HiFi` in `Correct.h` / `Correct.cpp`.
//!
//! The evidence list records one observation per query base: whether the
//! overlapping read agrees with the query at that base (`type & 1 == 0`) or
//! differs (`type & 1 == 1`), together with the observed base and the target
//! position. `push_info` folds the observations of one query site into per-site
//! allele counts, and `generate_haplotypes_naive_hifi` decides which overlaps
//! support a second haplotype.
//!
//! Several reference behaviours look like bugs and are reproduced on purpose:
//!
//! * the allele counts wrap modulo 2^32 when `occ_2` is derived,
//! * the 32-base multi-check filter reads stale `s`/`t` statistics because the
//!   reference reuses the pointers left over from earlier loops,
//! * the duplicate test at the end of the multi-check can never score a site
//!   because `l` is advanced at the end of every iteration, and
//! * the scores of overlaps that are not promoted are cleared only when the
//!   overlap ends the pass with `is_match == 1`.

use crate::container::ksort::{radix_sort_by_key_32, radix_sort_u64};
use crate::overlap::region::OverlapRegion;
use crate::reads::read_store::ReadStore;
use crate::sequence::encoding::SEQUENCE_NT6_TABLE;

/// `WINDOW_MAX_SIZE`: the largest window the phasing pass handles.
///
/// `WINDOW` is 375 and `HA_MIN_OV_DIFF` is 0.02, so the reference constant is
/// `375 + (int)(1.0 / 0.02) + 3`.
pub const WINDOW_MAXIMUM_SIZE: usize = 375 + 50 + 3;

/// `haplotype_evdience`: one observation of one query base against one
/// overlapping read.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct HaplotypeEvidence {
    /// Position of the observation in the query read.
    pub site: u32,
    /// Index of the overlapping read in the overlap list, `overlapID`.
    pub overlap_identifier: u32,
    /// Position of the observation in the overlapping read, or the allele id
    /// assigned by [`push_info`] for a match.
    pub overlap_site: u32,
    /// Bit 0: mismatch; bit 1: homopolymer.
    pub evidence_type: u8,
    /// Number of observations collapsed into this one.
    pub coverage: u32,
    /// Observed base, as the ASCII character recovered from the read.
    pub mismatch_base: u8,
}

impl HaplotypeEvidence {
    /// `hh_tp`: zero for a match, one for a mismatch.
    pub fn is_mismatch(&self) -> bool {
        self.evidence_type & 1 == 1
    }

    /// `hh_hp`: whether the observation is marked as a homopolymer.
    pub fn is_homopolymer(&self) -> bool {
        (self.evidence_type >> 1) & 1 == 1
    }

    /// `hh_bq`: the target base quality packed into bits 2-3.
    pub fn base_quality(&self) -> u8 {
        (self.evidence_type >> 2) & 3
    }

    /// `hh_wq`: the target window quality packed into bits 4-5.
    pub fn window_quality(&self) -> u8 {
        (self.evidence_type >> 4) & 3
    }
}

/// `SnpStats`: the folded allele counts of one query site.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SnpStats {
    /// Identifier assigned when the record was pushed, `id`.
    pub identifier: u32,
    /// Number of overlaps that carry an informative allele.
    pub overlap_num: u32,
    /// Observations that agree with the query.
    pub occ_0: u32,
    /// Observations that carry this allele.
    pub occ_1: u32,
    /// Remaining observations, `occ_2 - occ_0 - occ_1` with 32-bit wrap.
    pub occ_2: u32,
    pub homopolymer_num: u32,
    pub non_homopolymer_num: u32,
    /// One when the allele supports a second haplotype, minus one when reset.
    pub score: i32,
    /// Position of the site in the query read.
    pub site: u32,
    /// Whether the folded observations came from a homopolymer.
    pub is_homopolymer: u8,
}

/// `haplotype_evdience_alloc`: the evidence list and the statistics it is
/// folded into. The matrix and dynamic-programming buffers of the reference
/// are not used by the HiFi path and are left out.
pub struct HaplotypeEvidenceAlloc {
    /// Evidence records; [`Self::length`] entries are live.
    pub list: Vec<HaplotypeEvidence>,
    pub sub_list_start: u32,
    pub sub_list_length: u32,
    /// `length`: number of live evidence records.
    pub length: usize,
    /// Capacity the reference tracks for `list`.
    pub size: u32,
    /// `flag`: per-window count of mismatch observations.
    pub flag: [u8; WINDOW_MAXIMUM_SIZE],
    pub core_snp: u32,
    /// Number of overlaps the statistics were built from.
    pub overlap: u32,
    pub result_stat: SnpStats,
    /// Folded per-site allele statistics.
    pub snp_stat: Vec<SnpStats>,
    /// Number of positions that passed the window threshold.
    pub nn_snp: u32,
    /// Scratch keys of the haplotype pass.
    pub snp_srt: Vec<u64>,
}

impl Default for HaplotypeEvidenceAlloc {
    fn default() -> Self {
        Self::new()
    }
}

impl HaplotypeEvidenceAlloc {
    /// `InitHaplotypeEvdience`.
    pub fn new() -> Self {
        Self {
            list: Vec::new(),
            sub_list_start: 0,
            sub_list_length: 0,
            length: 0,
            size: 100,
            flag: [0; WINDOW_MAXIMUM_SIZE],
            core_snp: 0,
            overlap: 0,
            result_stat: SnpStats::default(),
            snp_stat: Vec::new(),
            nn_snp: 0,
            snp_srt: Vec::new(),
        }
    }

    /// `ResizeInitHaplotypeEvdience`.
    ///
    /// Note that the reference does not clear `snp_srt` here.
    pub fn resize_init(&mut self) {
        self.nn_snp = 0;
        self.snp_stat.clear();
        self.length = 0;
        self.sub_list_start = 0;
        self.sub_list_length = 0;
        self.flag.fill(0);
    }

    /// `addHaplotypeEvdience`: appends one record at the logical end.
    pub fn add(&mut self, evidence: HaplotypeEvidence) {
        let index = self.length;
        if index == self.list.len() {
            self.list.push(evidence);
        } else {
            self.list[index] = evidence;
        }
        self.length += 1;
        if self.length as u32 > self.size {
            self.size = (self.length as u32).next_power_of_two().max(1);
        }
    }

    /// `SetSnpMatrix` with `set_matrix == 0`: keeps the capacity of
    /// `snp_stat`, empties it and records the overlap count.
    pub fn set_snp_matrix(&mut self, nn_snp: u32, overlap_count: u64) {
        self.snp_stat.reserve(nn_snp as usize);
        self.snp_stat.clear();
        self.overlap = overlap_count as u32;
        self.core_snp = 0;
    }
}

/// `is_st_bs(s, rr, mm)`: the stand-alone test of the second correction round.
///
/// It is always false in the HiFi path because `mm` is `(uint64_t)-1` there.
pub(crate) fn is_standalone(stats: &SnpStats, rate: f64, maximum: u64) -> bool {
    maximum != u64::MAX
        && u64::from(stats.overlap_num) + maximum >= u64::from(stats.occ_0)
        && f64::from(stats.occ_0) * rate + f64::from(stats.overlap_num) >= f64::from(stats.occ_0)
}

/// `push_info`: folds the evidence of one query site into allele statistics.
///
/// `a_start`/`a_length` select the input records inside [`Self::list`];
/// `out_start` is where the compacted records are written. The output never
/// overtakes the input, so the compaction is safe in place, exactly like the
/// reference.
///
/// Returns the number of records written.
/// `get_wqual` (`Correct.cpp` 9616): the quality of one target base and the
/// minimum mean quality of the `window`-base sliding window around it.
pub fn window_quality(
    reads: &ReadStore,
    read_id: u32,
    position: u64,
    reverse: bool,
    window: u64,
) -> (u64, u64) {
    let length = reads.read_length(read_id as usize);
    let start = if position + 1 >= window {
        position + 1 - window
    } else {
        0
    };
    let end = (position + window).min(length);
    let levels = reads.read_quality(read_id as usize, start as i64, end as i64, reverse);
    if levels.is_empty() {
        return (0, 0);
    }
    let base_quality = u64::from(levels[(position - start) as usize]);
    let window = window as usize;
    let mut count = 0usize;
    let mut sum = 0u64;
    while count < window && count < levels.len() {
        sum += u64::from(levels[count]);
        count += 1;
    }
    let mut minimum = sum;
    for index in count..levels.len() {
        sum -= u64::from(levels[index - window]);
        sum += u64::from(levels[index]);
        if sum < minimum {
            minimum = sum;
        }
    }
    let mean = (minimum as f64 / count as f64).round() as u64;
    (base_quality, mean)
}

pub fn push_info(
    hap: &mut HaplotypeEvidenceAlloc,
    a_start: usize,
    a_length: usize,
    out_start: usize,
    overlaps: Option<&[OverlapRegion]>,
    quality_reads: Option<&ReadStore>,
    quality_window: u64,
) -> usize {
    let zone_end = a_start + a_length;
    radix_sort_by_key_32(&mut hap.list[a_start..zone_end], &|record| {
        record.overlap_identifier
    });

    let mut occ_0: u64 = 0;
    let mut occ_1 = [0u64; 6];
    let mut occ_2: u64 = 0;
    let mut difference: u64 = 0;
    let mut reverse_number: u64 = 0;
    let mut is_homopolymer: u8 = 0;

    // Collapse the records that share an overlap identifier, keeping the first
    // one of every run and accumulating its counts.
    let mut kept = 0usize;
    let mut block_start = 0usize;
    for k in 1..=a_length {
        let closes = k == a_length
            || hap.list[a_start + k].overlap_identifier
                != hap.list[a_start + block_start].overlap_identifier;
        if !closes {
            continue;
        }
        let record = hap.list[a_start + block_start];
        hap.list[a_start + kept] = record;
        if !record.is_mismatch() {
            occ_0 += u64::from(record.coverage);
            if let Some(overlaps) = overlaps {
                if overlaps[record.overlap_identifier as usize].y_strand == 0 {
                    reverse_number += u64::from(record.coverage);
                }
            }
            if is_homopolymer == 0 && record.is_homopolymer() {
                is_homopolymer = 1;
            }
        } else {
            occ_1[SEQUENCE_NT6_TABLE[record.mismatch_base as usize] as usize] +=
                u64::from(record.coverage);
            difference += u64::from(record.coverage);
        }
        occ_2 += u64::from(record.coverage);
        kept += 1;
        block_start = k;
    }
    let a_length = kept;

    // At least one read must agree with the query, and a single mismatch is a
    // sequencing error rather than an allele.
    if occ_0 == 0 || difference <= 1 {
        return 0;
    }
    if overlaps.is_some() && reverse_number == occ_0 {
        return 0;
    }
    if overlaps.is_none() {
        reverse_number = occ_2;
    } else {
        reverse_number += 1;
    }

    let site = hap.list[a_start].site;
    let mut allele_identifier = [u64::MAX; 6];
    let mut pushed = 0usize;
    for base in 0..4 {
        if occ_1[base] >= 2 {
            let identifier = hap.snp_stat.len() as u32;
            let occ_0_field = (1 + occ_0) as u32;
            let occ_1_field = occ_1[base] as u32;
            let occ_2_field = occ_2
                .wrapping_sub(u64::from(occ_0_field))
                .wrapping_sub(u64::from(occ_1_field)) as u32;
            hap.snp_stat.push(SnpStats {
                identifier,
                overlap_num: reverse_number as u32,
                occ_0: occ_0_field,
                occ_1: occ_1_field,
                occ_2: occ_2_field,
                homopolymer_num: 0,
                non_homopolymer_num: 0,
                score: -1,
                site,
                is_homopolymer,
            });
            allele_identifier[base] = u64::from(identifier);
            pushed += 1;
        } else {
            allele_identifier[base] = u64::MAX;
        }
    }
    allele_identifier[4] = u64::MAX;
    allele_identifier[5] = u64::MAX;
    if pushed == 0 {
        return 0;
    }

    // Keep only the observations that can be mapped to a pushed allele. A
    // match points at the last allele of the site; a mismatch keeps its target
    // position in `coverage` (the reference reuses the field).
    let mut out = 0usize;
    for index in 0..a_length {
        let mut record = hap.list[a_start + index];
        let target_position = record.overlap_site;
        if !record.is_mismatch() {
            record.overlap_site = hap.snp_stat.len() as u32 - 1;
        } else if allele_identifier[SEQUENCE_NT6_TABLE[record.mismatch_base as usize] as usize]
            != u64::MAX
        {
            record.coverage = target_position;
            record.overlap_site = allele_identifier
                [SEQUENCE_NT6_TABLE[record.mismatch_base as usize] as usize]
                as u32;
        } else {
            continue;
        }
        // `v8 != NULL`: pack the base and window qualities of the target into
        // the evidence type (`hh_bq` bits 2-3, `hh_wq` bits 4-5).
        if let Some(reads) = quality_reads {
            let overlaps = overlaps.expect("quality folding needs the overlap list");
            let target = &overlaps[record.overlap_identifier as usize];
            let (base_quality, mean_quality) = window_quality(
                reads,
                target.y_identifier,
                u64::from(target_position),
                target.y_strand != 0,
                quality_window,
            );
            record.evidence_type |= ((base_quality & 3) << 2) as u8;
            record.evidence_type |= ((mean_quality & 3) << 4) as u8;
        }
        hap.list[out_start + out] = record;
        out += 1;
    }
    out
}

/// `generate_haplotypes_naive_HiFi`: splits the overlaps of one read into the
/// haplotype they support.
///
/// `up` is the fraction of the alignment two informative alleles must cover,
/// `multi_check` enables the weak-allele pass, and `standard_rate` /
/// `standard_maximum` are the `is_st_bs` parameters (zero and `u64::MAX` in the
/// HiFi path). `mapped_overlaps_length` is decremented for every overlap that
/// turns out to support the second haplotype.
#[allow(clippy::too_many_arguments)]
pub fn generate_haplotypes_naive_hifi(
    hap: &mut HaplotypeEvidenceAlloc,
    overlaps: &mut [OverlapRegion],
    up: f64,
    multi_check: u64,
    standard_rate: f64,
    standard_maximum: u64,
    informative_coverage: i64,
    single_haplotype_coverage: i64,
    mapped_overlaps_length: &mut i64,
) {
    if hap.length == 0 {
        return;
    }

    // Drop the sites that neighbour another site, then fold the evidence into
    // the surviving sites. `snp_stat` and `list` are both ordered by site.
    let snp_count = hap.snp_stat.len();
    let mut k = 1usize;
    let mut l = 0usize;
    let mut i = 0usize;
    let mut kept_statistics = 0usize;
    let mut kept_evidence = 0usize;
    while k <= snp_count {
        let closes = k == snp_count || hap.snp_stat[k].site != hap.snp_stat[l].site;
        if closes {
            if l > 0 && hap.snp_stat[l].site == hap.snp_stat[l - 1].site.wrapping_add(1) {
                l = k;
                k += 1;
                continue;
            }
            if k < snp_count
                && hap.snp_stat[l].site.wrapping_add(1) == hap.snp_stat[k].site
            {
                l = k;
                k += 1;
                continue;
            }

            let site = hap.snp_stat[l].site;
            while i < hap.length && hap.list[i].site != site {
                i += 1;
            }
            assert!(
                i < hap.length && hap.list[i].site == site,
                "every kept site must have evidence"
            );
            let offset = (l - kept_statistics) as u32;
            while i < hap.length && hap.list[i].site == site {
                let mut record = hap.list[i];
                record.overlap_site = record.overlap_site.wrapping_sub(offset);
                hap.list[kept_evidence] = record;
                kept_evidence += 1;
                i += 1;
            }
            while l < k {
                hap.snp_stat[kept_statistics] = hap.snp_stat[l];
                kept_statistics += 1;
                l += 1;
            }
        }
        k += 1;
    }
    hap.snp_stat.truncate(kept_statistics);
    hap.length = kept_evidence;

    if hap.snp_stat.is_empty() || hap.length == 0 {
        return;
    }

    // Rank the overlaps by the number of strong alleles they carry. The rank is
    // stored inverted so that the largest count sorts first.
    hap.snp_srt.clear();
    radix_sort_by_key_32(&mut hap.list[..hap.length], &|record| {
        record.overlap_identifier
    });
    let length = hap.length;
    let mut k = 1usize;
    let mut l = 0usize;
    while k <= length {
        let closes =
            k == length || hap.list[k].overlap_identifier != hap.list[l].overlap_identifier;
        if closes {
            let mut count = 0u64;
            for index in l..k {
                let record = hap.list[index];
                if !record.is_mismatch() {
                    continue;
                }
                let stats = hap.snp_stat[record.overlap_site as usize];
                assert_eq!(stats.site, record.site, "evidence and statistics must agree");
                if stats.occ_0 < 2 || stats.occ_1 < 2 {
                    continue;
                }
                if is_standalone(&stats, standard_rate, standard_maximum) {
                    continue;
                }
                if i64::from(stats.occ_0) >= single_haplotype_coverage
                    && i64::from(stats.occ_1) >= informative_coverage
                {
                    count += 1;
                }
            }
            if count > 0 {
                let key = (u64::from(u32::MAX).wrapping_sub(count) << 32)
                    .wrapping_add(l as u64);
                hap.snp_srt.push(key);
            }
            l = k;
        }
        k += 1;
    }

    // The stale `s` and `t` pointers of the reference are modelled by indices
    // into `snp_stat`; they deliberately outlive the loops that assigned them.
    let mut stale_s: Option<usize> = None;
    let mut stale_t: Option<usize> = None;

    if !hap.snp_srt.is_empty() {
        radix_sort_u64(&mut hap.snp_srt);

        for key_index in 0..hap.snp_srt.len() {
            let l = hap.snp_srt[key_index] as u32 as usize;
            let mut count = 0u64;
            let mut index = l;
            while index < length
                && hap.list[index].overlap_identifier == hap.list[l].overlap_identifier
            {
                if hap.list[index].is_mismatch() {
                    stale_s = Some(hap.list[index].overlap_site as usize);
                    let stats = hap.snp_stat[hap.list[index].overlap_site as usize];
                    if stats.occ_0 >= 2
                        && stats.occ_1 >= 2
                        && !is_standalone(&stats, standard_rate, standard_maximum)
                        && i64::from(stats.occ_0) >= single_haplotype_coverage
                        && i64::from(stats.occ_1) >= informative_coverage
                    {
                        count += 1;
                    }
                }
                index += 1;
            }
            if count == 0 {
                continue;
            }

            let overlap_index = hap.list[l].overlap_identifier as usize;
            if overlaps[overlap_index].is_match == 1 {
                overlaps[overlap_index].is_match = 2;
            }
            let mut index = l;
            while index < length
                && hap.list[index].overlap_identifier == hap.list[l].overlap_identifier
            {
                let record = hap.list[index];
                if record.is_mismatch() {
                    stale_s = Some(record.overlap_site as usize);
                    hap.snp_stat[record.overlap_site as usize].score = 1;
                } else {
                    // A match of a promoted overlap no longer counts towards
                    // the allele totals of its site.
                    let first = record.overlap_site as usize;
                    stale_s = Some(first);
                    let site = hap.snp_stat[first].site;
                    let mut z = first as i64;
                    while z >= 0 {
                        stale_t = Some(z as usize);
                        if hap.snp_stat[z as usize].site != site {
                            break;
                        }
                        let stats = &mut hap.snp_stat[z as usize];
                        stats.occ_0 = stats.occ_0.wrapping_sub(record.coverage);
                        assert!(stats.occ_0 >= 1, "occ_0 must stay positive");
                        if standard_maximum != u64::MAX
                            && overlaps[overlap_index].y_strand == 0
                        {
                            stats.overlap_num =
                                stats.overlap_num.wrapping_sub(record.coverage);
                            assert!(stats.overlap_num >= 1, "overlap_num must stay positive");
                        }
                        z -= 1;
                    }
                }
                index += 1;
            }
        }

        // Second pass: an overlap is promoted when one of its mismatches now
        // carries a scored allele.
        for key_index in 0..hap.snp_srt.len() {
            let l = hap.snp_srt[key_index] as u32 as usize;
            let mut count = 0u64;
            let mut index = l;
            while index < length
                && hap.list[index].overlap_identifier == hap.list[l].overlap_identifier
            {
                if hap.list[index].is_mismatch() {
                    stale_s = Some(hap.list[index].overlap_site as usize);
                    let stats = hap.snp_stat[hap.list[index].overlap_site as usize];
                    if stats.occ_0 >= 2
                        && stats.occ_1 >= 2
                        && !is_standalone(&stats, standard_rate, standard_maximum)
                        && stats.score == 1
                    {
                        count += 1;
                    }
                }
                index += 1;
            }
            let overlap_index = hap.list[l].overlap_identifier as usize;
            if overlaps[overlap_index].is_match == 1 && count > 0 {
                overlaps[overlap_index].is_match = 2;
            }
        }

        // Clear the scores of every overlap that was not promoted.
        let mut k = 1usize;
        let mut l = 0usize;
        while k <= length {
            let closes =
                k == length || hap.list[k].overlap_identifier != hap.list[l].overlap_identifier;
            if closes {
                let overlap_index = hap.list[l].overlap_identifier as usize;
                if overlaps[overlap_index].is_match == 1 {
                    for index in l..k {
                        if hap.list[index].is_mismatch() {
                            hap.snp_stat[hap.list[index].overlap_site as usize].score = -1;
                        }
                    }
                }
                l = k;
            }
            k += 1;
        }
    }

    // Weak-allele pass: an overlap that disagrees with the query at several
    // well-separated weak sites supports the second haplotype as well.
    if multi_check != 0 {
        hap.snp_srt.clear();
        let mut k = 1usize;
        let mut l = 0usize;
        while k <= length {
            let closes =
                k == length || hap.list[k].overlap_identifier != hap.list[l].overlap_identifier;
            if closes {
                let overlap_index = hap.list[l].overlap_identifier as usize;
                if overlaps[overlap_index].is_match == 2 {
                    l = k;
                    k += 1;
                    continue;
                }

                let base = hap.snp_srt.len();
                let mut count = 0u64;
                for index in l..k {
                    let record = hap.list[index];
                    if !record.is_mismatch() {
                        continue;
                    }
                    stale_s = Some(record.overlap_site as usize);
                    let stats = hap.snp_stat[record.overlap_site as usize];
                    if stats.occ_0 < 2 || stats.occ_1 < 2 {
                        continue;
                    }
                    if is_standalone(&stats, standard_rate, standard_maximum) {
                        continue;
                    }
                    if i64::from(stats.occ_0) >= single_haplotype_coverage
                        && i64::from(stats.occ_1) >= informative_coverage
                    {
                        continue;
                    }
                    if stats.score == 1 {
                        continue;
                    }
                    count += 1;
                    hap.snp_srt.push(u64::from(record.overlap_site));
                }
                // The just-pushed block is scratch space in the reference: it
                // is taken out of the logical length again.
                let mut block: Vec<u64> = hap.snp_srt[base..].to_vec();
                hap.snp_srt.truncate(base);

                let alignment_length = overlaps[overlap_index].alignment_length;
                if count as f64 >= f64::from(alignment_length) * up {
                    radix_sort_u64(&mut block);
                    let mut written = 0usize;
                    for index in 0..count as usize {
                        let current = block[index] as usize;
                        if index > 0 {
                            stale_s = Some(block[index - 1] as usize);
                        }
                        if index + 1 < count as usize {
                            stale_t = Some(block[index + 1] as usize);
                        }
                        if let Some(previous) = stale_s {
                            if hap.snp_stat[previous].site.wrapping_add(32)
                                > hap.snp_stat[current].site
                            {
                                continue;
                            }
                        }
                        if let Some(next) = stale_t {
                            if hap.snp_stat[current].site.wrapping_add(32)
                                > hap.snp_stat[next].site
                            {
                                continue;
                            }
                        }
                        block[written] = block[index];
                        written += 1;
                    }
                    if written >= 2 {
                        hap.snp_srt.extend_from_slice(&block[..written]);
                    }
                }
                l = k;
            }
            k += 1;
        }

        if !hap.snp_srt.is_empty() {
            radix_sort_u64(&mut hap.snp_srt);
            // The reference advances `l` at the end of every iteration, so the
            // duplicate test cannot fire; it is reproduced as written.
            let mut k = 1usize;
            let mut l = 0usize;
            while k <= hap.snp_srt.len() {
                if k == hap.snp_srt.len() || hap.snp_srt[k] != hap.snp_srt[l] {
                    if k - l >= 2 {
                        hap.snp_stat[hap.snp_srt[l] as usize].score = 1;
                    }
                }
                l = k;
                k += 1;
            }
        }
    }

    // Final classification: every overlap that supports the second haplotype
    // is marked strong and removed from the mapped length.
    let mut k = 1usize;
    let mut l = 0usize;
    while k <= length {
        let closes =
            k == length || hap.list[k].overlap_identifier != hap.list[l].overlap_identifier;
        if closes {
            let overlap_index = hap.list[l].overlap_identifier as usize;
            if overlaps[overlap_index].is_match == 2 {
                overlaps[overlap_index].strong = 1;
                *mapped_overlaps_length -= overlap_span(&overlaps[overlap_index]);
            } else if overlaps[overlap_index].is_match == 1 {
                for index in l..k {
                    let record = hap.list[index];
                    // `hh_tp` is always zero or one, so both arms are taken.
                    let _ = record.evidence_type & 1;
                    let stats = hap.snp_stat[record.overlap_site as usize];
                    if stats.score == 1
                        && !(stats.occ_0 < 2 || stats.occ_1 < 2)
                        && !is_standalone(&stats, standard_rate, standard_maximum)
                    {
                        overlaps[overlap_index].strong = 1;
                        if record.is_mismatch() {
                            overlaps[overlap_index].is_match = 2;
                            *mapped_overlaps_length -= overlap_span(&overlaps[overlap_index]);
                            break;
                        }
                    }
                }
            }
            l = k;
        }
        k += 1;
    }
}

/// `generate_haplotypes_sv` (`Correct.cpp` 9115): the structural-variant
/// promotion used after the ONT indel pass. Unlike the HiFi variant it has no
/// stand-alone test and runs a second promotion round over the alleles that
/// the first round scored.
pub fn generate_haplotypes_sv(
    hap: &mut HaplotypeEvidenceAlloc,
    overlaps: &mut [OverlapRegion],
    informative_coverage: i64,
    single_haplotype_coverage: i64,
    mapped_overlaps_length: &mut i64,
) {
    if hap.snp_stat.is_empty() || hap.length == 0 {
        return;
    }
    let length = hap.length;
    hap.snp_srt.clear();
    radix_sort_by_key_32(&mut hap.list[..length], &|record| {
        record.overlap_identifier
    });

    let real_allele = |stats: &SnpStats| -> bool {
        stats.occ_0 >= 2
            && stats.occ_1 >= 2
            && i64::from(stats.occ_0) >= single_haplotype_coverage
            && i64::from(stats.occ_1) >= informative_coverage
    };

    let mut k = 1usize;
    let mut l = 0usize;
    while k <= length {
        if k == length || hap.list[k].overlap_identifier != hap.list[l].overlap_identifier {
            let mut count = 0u64;
            for index in l..k {
                let record = hap.list[index];
                if !record.is_mismatch() {
                    continue;
                }
                let stats = &hap.snp_stat[record.overlap_site as usize];
                assert_eq!(stats.site, record.site);
                if real_allele(stats) {
                    count += 1;
                }
            }
            if count > 0 {
                let key = (u64::from(u32::MAX).wrapping_sub(count) << 32)
                    .wrapping_add(l as u64);
                hap.snp_srt.push(key);
            }
            l = k;
        }
        k += 1;
    }

    if !hap.snp_srt.is_empty() {
        radix_sort_u64(&mut hap.snp_srt);

        // First promotion: score every real allele of the ranked overlaps.
        for key_index in 0..hap.snp_srt.len() {
            let l = hap.snp_srt[key_index] as u32 as usize;
            let mut count = 0u64;
            let mut index = l;
            while index < length
                && hap.list[index].overlap_identifier == hap.list[l].overlap_identifier
            {
                if hap.list[index].is_mismatch()
                    && real_allele(&hap.snp_stat[hap.list[index].overlap_site as usize])
                {
                    count += 1;
                }
                index += 1;
            }
            if count == 0 {
                continue;
            }
            let overlap_index = hap.list[l].overlap_identifier as usize;
            if overlaps[overlap_index].is_match == 1 {
                overlaps[overlap_index].is_match = 2;
            }
            let mut index = l;
            while index < length
                && hap.list[index].overlap_identifier == hap.list[l].overlap_identifier
            {
                let record = hap.list[index];
                if record.is_mismatch() {
                    hap.snp_stat[record.overlap_site as usize].score = 1;
                } else {
                    let site = record.overlap_site as usize;
                    hap.snp_stat[site].occ_0 -= 1;
                    assert!(hap.snp_stat[site].occ_0 >= 1);
                }
                index += 1;
            }
        }

        // Second promotion: an overlap still marked as a match but carrying a
        // scored allele joins the second haplotype.
        for key_index in 0..hap.snp_srt.len() {
            let l = hap.snp_srt[key_index] as u32 as usize;
            let mut count = 0u64;
            let mut index = l;
            while index < length
                && hap.list[index].overlap_identifier == hap.list[l].overlap_identifier
            {
                let record = hap.list[index];
                if record.is_mismatch() {
                    let stats = &hap.snp_stat[record.overlap_site as usize];
                    if stats.occ_0 >= 2 && stats.occ_1 >= 2 && stats.score == 1 {
                        count += 1;
                    }
                }
                index += 1;
            }
            let overlap_index = hap.list[l].overlap_identifier as usize;
            if overlaps[overlap_index].is_match == 1 && count > 0 {
                overlaps[overlap_index].is_match = 2;
            }
        }

        // Reset the scores of the overlaps that stayed in the first haplotype.
        let mut k = 1usize;
        let mut l = 0usize;
        while k <= length {
            if k == length || hap.list[k].overlap_identifier != hap.list[l].overlap_identifier {
                let overlap_index = hap.list[l].overlap_identifier as usize;
                if overlaps[overlap_index].is_match == 1 {
                    for index in l..k {
                        if hap.list[index].is_mismatch() {
                            hap.snp_stat[hap.list[index].overlap_site as usize].score = -1;
                        }
                    }
                }
                l = k;
            }
            k += 1;
        }
    }

    // Final marking: promoted overlaps become strong; a match overlap with a
    // scored real allele is promoted on the spot.
    let mut k = 1usize;
    let mut l = 0usize;
    while k <= length {
        if k == length || hap.list[k].overlap_identifier != hap.list[l].overlap_identifier {
            let overlap_index = hap.list[l].overlap_identifier as usize;
            if overlaps[overlap_index].is_match == 2 {
                overlaps[overlap_index].strong = 1;
                *mapped_overlaps_length -= overlap_span(&overlaps[overlap_index]);
            } else if overlaps[overlap_index].is_match == 1 {
                for index in l..k {
                    let record = hap.list[index];
                    let stats = &hap.snp_stat[record.overlap_site as usize];
                    if stats.score == 1 && !(stats.occ_0 < 2 || stats.occ_1 < 2) {
                        overlaps[overlap_index].strong = 1;
                        if record.is_mismatch() {
                            overlaps[overlap_index].is_match = 2;
                            *mapped_overlaps_length -= overlap_span(&overlaps[overlap_index]);
                            break;
                        }
                    }
                }
            }
            l = k;
        }
        k += 1;
    }
}

/// `x_pos_e + 1 - x_pos_s` with the 32-bit wrap of the reference.
fn overlap_span(overlap: &OverlapRegion) -> i64 {
    i64::from(
        overlap
            .x_position_end
            .wrapping_add(1)
            .wrapping_sub(overlap.x_position_start),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn window_quality_matches_the_sliding_minimum() {
        let mut reads = ReadStore::default();
        // Phred 40 ("I") is level 3, "!" is level 0.
        reads.push_read_with_quality("read", b"ACGTACGT", Some(b"I!IIIIII"));
        // position 1, window 3: levels [3,0,3,3], base 0, min window sum 6.
        assert_eq!(window_quality(&reads, 0, 1, false, 3), (0, 2));
        // A uniform track returns the level itself.
        reads.push_read_with_quality("flat", b"ACGTACGT", Some(b"IIIIIIII"));
        assert_eq!(window_quality(&reads, 1, 4, false, 5), (3, 3));
    }

    #[test]
    fn the_window_capacity_matches_the_reference() {
        assert_eq!(WINDOW_MAXIMUM_SIZE, 428);
        let hap = HaplotypeEvidenceAlloc::new();
        assert_eq!(hap.flag.len(), WINDOW_MAXIMUM_SIZE);
        assert_eq!(hap.list.len(), 0);
        assert_eq!(hap.snp_srt.len(), 0);
    }

    #[test]
    fn standalone_testing_is_disabled_for_the_hifi_path() {
        let stats = SnpStats {
            overlap_num: 100,
            occ_0: 10,
            ..SnpStats::default()
        };
        assert!(!is_standalone(&stats, 0.0, u64::MAX));
        assert!(is_standalone(&stats, 0.5, 0));
    }

    #[test]
    fn resize_keeps_the_sorted_keys() {
        let mut hap = HaplotypeEvidenceAlloc::new();
        hap.snp_srt.push(7);
        hap.resize_init();
        assert_eq!(hap.length, 0);
        assert!(hap.snp_stat.is_empty());
        assert_eq!(hap.snp_srt, vec![7], "the reference does not clear snp_srt");
    }
}

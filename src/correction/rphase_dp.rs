//! The ONT dynamic-programming phasing of `rphase_hc`, the counterparts of
//! `gen_rphase_dp`, `gen_rphase_dp0_single_path`, `comput_sc_rphase`,
//! `is_hpc_vec`, `get_hq_value` and `fill_incom` in `Correct.cpp`.
//!
//! The HiFi path never enters here (`dp == NULL`); ONT passes the chaining
//! `Chain_Data` and the query quality levels, so the allele scores are decided
//! by a longest-path DP over compatible sites instead of the naive vote.

use crate::container::ksort::{radix_sort_by_key_32, radix_sort_u64};
use crate::correction::haplotype_evidence::{
    is_standalone, HaplotypeEvidence, HaplotypeEvidenceAlloc, SnpStats,
};
use crate::overlap::region::{ChainData, OverlapRegion};
use crate::reads::read_store::ReadStore;
use crate::sequence::encoding::SEQUENCE_NT6_TABLE;

/// `sc_wn`: the quality window of `get_wqual`.
const QUALITY_WINDOW: u64 = 5;

/// `get_wqual` over an already-retrieved level array (the `va` variant).
fn window_quality_levels(levels: &[u8], position: u64, window: u64) -> (u64, u64) {
    let length = levels.len() as u64;
    let start = if position + 1 >= window {
        position + 1 - window
    } else {
        0
    };
    let end = (position + window).min(length);
    if start >= end {
        return (0, 0);
    }
    let base = u64::from(levels[(position - start) as usize]);
    let window = window as usize;
    let mut count = 0usize;
    let mut sum = 0u64;
    while count < window && (start as usize + count) < end as usize {
        sum += u64::from(levels[start as usize + count]);
        count += 1;
    }
    let mut minimum = sum;
    for index in count..(end - start) as usize {
        sum -= u64::from(levels[start as usize + index - window]);
        sum += u64::from(levels[start as usize + index]);
        if sum < minimum {
            minimum = sum;
        }
    }
    (base, (minimum as f64 / count as f64).round() as u64)
}

/// `comput_sc_rphase` (`Correct.cpp` 9244): whether two alleles share both a
/// matching and a disagreeing overlap, `i64::MIN` when they are incompatible.
fn comput_sc_rphase(
    stats: &[SnpStats],
    id: usize,
    jd: usize,
    evidence: &[HaplotypeEvidence],
) -> i64 {
    let ai = &stats[id];
    let aj = &stats[jd];
    if ai.site == aj.site {
        return i64::MIN;
    }
    let iz = &evidence[ai.non_homopolymer_num as usize..ai.homopolymer_num as usize];
    let jz = &evidence[aj.non_homopolymer_num as usize..aj.homopolymer_num as usize];
    let mut nn = [0i64; 2];
    let mut jk = 0usize;
    for ik in 0..iz.len() {
        while jk < jz.len() && jz[jk].overlap_identifier < iz[ik].overlap_identifier {
            jk += 1;
        }
        if jk < jz.len() && jz[jk].overlap_identifier == iz[ik].overlap_identifier {
            let mut fi = 2u8;
            if !iz[ik].is_mismatch() {
                fi = 0;
            } else if iz[ik].overlap_site == id as u32 {
                fi = 1;
            }
            let mut fj = 2u8;
            if !jz[jk].is_mismatch() {
                fj = 0;
            } else if jz[jk].overlap_site == jd as u32 {
                fj = 1;
            }
            if fi == 2
                && fj == 2
                && iz[ik].overlap_site != u32::MAX
                && jz[jk].overlap_site != u32::MAX
            {
                fi = 0;
                fj = 0;
            }
            if fi == 2 || fj == 2 || fi != fj {
                return i64::MIN;
            }
            nn[fi as usize] += 1;
        }
    }
    if nn[0] > 0 && nn[1] > 0 {
        1
    } else {
        i64::MIN
    }
}

/// `is_hpc_vec` (`Correct.cpp` 9398): whether removing the homopolymer
/// observations drops the allele below the real-allele thresholds.
fn is_hpc_vec(
    stat: &mut SnpStats,
    id: usize,
    evidence: &[HaplotypeEvidence],
    single_haplotype_coverage: i64,
    informative_coverage: i64,
) -> bool {
    let n0 = stat.occ_0;
    let n1 = stat.occ_1;
    let slice = &evidence[stat.non_homopolymer_num as usize..stat.homopolymer_num as usize];
    for record in slice {
        if !record.is_homopolymer() {
            continue;
        }
        if !record.is_mismatch() {
            stat.occ_0 -= record.coverage;
        } else if record.overlap_site == id as u32 {
            stat.occ_1 -= record.coverage;
        }
    }
    let failed = stat.occ_0 < 2
        || stat.occ_1 < 2
        || i64::from(stat.occ_0) < single_haplotype_coverage
        || i64::from(stat.occ_1) < informative_coverage;
    stat.occ_0 = n0;
    stat.occ_1 = n1;
    failed
}

/// `get_hq_value` (`Correct.cpp` 9415): the low/high-quality counts of the two
/// alleles, with the query base folded into the matching allele.
fn get_hq_value(
    stat: &SnpStats,
    id: usize,
    evidence: &[HaplotypeEvidence],
    query_levels: &[u8],
) -> (i64, i64, i64, i64) {
    let slice = &evidence[stat.non_homopolymer_num as usize..stat.homopolymer_num as usize];
    let mut occ0 = [0i64; 2];
    let mut occ1 = [0i64; 2];
    for record in slice {
        if record.is_mismatch() && record.overlap_site != id as u32 {
            continue;
        }
        let high = usize::from(record.base_quality() > 1);
        if !record.is_mismatch() {
            occ0[high] += 1;
        } else {
            occ1[high] += 1;
        }
    }
    if (stat.site as usize) < query_levels.len() {
        if query_levels[stat.site as usize] > 1 {
            occ0[1] += 1;
        } else {
            occ0[0] += 1;
        }
    }
    (occ0[0], occ0[1], occ1[0], occ1[1])
}

/// `fill_incom` (`Correct.cpp` 9597): records the sites an overlap missed
/// between its previous and current observation.
fn fill_incom(res: &mut Vec<u64>, oid: u64, site: u64, idx: &[u64]) {
    if res[oid as usize] == u64::MAX {
        res[oid as usize] = site;
        return;
    }
    let previous = res[oid as usize];
    let mut k = idx.len();
    while k > 0 && idx[k - 1] != previous {
        k -= 1;
    }
    for &value in &idx[k..] {
        res.push((oid << 32) | value);
    }
    res[oid as usize] = site;
}

/// `gen_rphase_dp0_single_path` (`Correct.cpp` 9429): the longest compatible
/// path over the alleles, scored with the quality-aware thresholds.
#[allow(clippy::too_many_arguments)]
fn dp0_single_path(
    stats: &mut [SnpStats],
    evidence: &[HaplotypeEvidence],
    dp: &mut ChainData,
    idx: &mut Vec<u64>,
    res: &mut Vec<u64>,
    het_cov: i64,
    hom_cov: i64,
    n_hap: i64,
    cut_rate: f64,
    cut_bd: u64,
    query_levels: &[u8],
    single_haplotype_coverage: i64,
    informative_coverage: i64,
) {
    let an = stats.len();
    if an == 0 {
        return;
    }
    dp.resize(an);
    idx.clear();
    res.clear();
    let mut cc = if het_cov > 0 { het_cov } else { hom_cov / n_hap } as f64;
    cc *= cut_rate;
    let cc = (cc as u64).max(cut_bd);

    let mut plus = 0i64;
    for i in 0..an {
        let mut max_f = 1i64;
        let mut max_j = -1i64;
        for j in (0..i).rev() {
            let sc = comput_sc_rphase(stats, i, j, evidence);
            if sc == i64::MIN {
                continue;
            }
            let sc = sc + i64::from(dp.score[j]);
            if sc > max_f {
                max_f = sc;
                max_j = j as i64;
            }
        }
        dp.score[i] = max_f as i32;
        dp.predecessor[i] = max_j;
        dp.occurrence[i] = 0;
        if max_f < plus {
            plus = max_f;
        }
    }

    for i in 0..an {
        let f = i64::from(dp.score[i]) - plus;
        let m = (u64::from(u32::MAX).wrapping_sub(f as u64) << 32) | i as u64;
        idx.push(m);
    }
    radix_sort_u64(idx);

    res.resize(an, 0);
    let mut rn = 0usize;
    for i in 0..an {
        let rn0 = rn;
        let mut k = idx[i] as u32 as i64;
        while k >= 0 && dp.occurrence[k as usize] == 0 {
            res[rn] = k as u64;
            rn += 1;
            dp.occurrence[k as usize] = 1;
            k = dp.predecessor[k as usize];
        }
        if rn0 == rn {
            continue;
        }
        let m = (u64::from(u32::MAX).wrapping_sub((rn - rn0) as u64) << 32) | rn0 as u64;
        idx.push(m);
    }

    radix_sort_u64(&mut idx[an..]);
    let ch_n = idx.len();
    for k in an..ch_n {
        let rn0 = idx[k] as u32 as usize;
        let rn = (u32::MAX - (idx[k] >> 32) as u32) as usize;
        let mut krn = rn;
        let mut i = 1usize;
        while i < rn
            && stats[res[rn0 + i] as usize]
                .site
                .wrapping_add(8)
                >= stats[res[rn0 + i - 1] as usize].site
        {
            i += 1;
        }
        if i >= rn {
            krn = 1;
        }

        for i in 0..rn {
            let site_index = res[rn0 + i] as usize;
            let (b0l, b0h, b1l, b1h) =
                get_hq_value(&stats[site_index], site_index, evidence, query_levels);
            let score = if krn > 1 {
                i64::from(
                    b0h > b0l
                        && b0h >= single_haplotype_coverage
                        && b1h > b1l
                        && b1h >= informative_coverage,
                ) * 2
                    - 1
            } else {
                let ok0 = b0h > b0l
                    && (b0h as f64) > (b0h + b0l) as f64 * 0.7
                    && b0h >= single_haplotype_coverage
                    && b0h >= cc as i64;
                let ok1 = b1h > b1l
                    && (b1h as f64) > (b1h + b1l) as f64 * 0.7
                    && b1h >= informative_coverage
                    && b1h >= cc as i64;
                i64::from(ok0 && ok1) * 2 - 1
            };
            stats[site_index].score = score as i32;
        }

        // A run of consecutive sites inherits the worst score of the run.
        let mut i = 0usize;
        while i < rn {
            let mut plus = stats[res[rn0 + i] as usize].score;
            let mut j = i + 1;
            while j < rn
                && stats[res[rn0 + i] as usize]
                    .site
                    .wrapping_sub(stats[res[rn0 + j] as usize].site)
                    == (j - i) as u32
            {
                if stats[res[rn0 + j] as usize].score == -1 {
                    plus = -1;
                }
                j += 1;
            }
            if plus == -1 {
                for index in i..j {
                    stats[res[rn0 + index] as usize].score = -1;
                }
            }
            i = j;
        }
    }
}

/// `gen_rphase_dp` (`Correct.cpp` 9648): the quality-aware allele filtering and
/// DP scoring of the ONT path.
#[allow(clippy::too_many_arguments)]
pub fn gen_rphase_dp(
    hap: &mut HaplotypeEvidenceAlloc,
    overlaps: &[OverlapRegion],
    reads: &ReadStore,
    query_id: u64,
    standard_rate: f64,
    standard_maximum: u64,
    dp: &mut ChainData,
    idx: &mut Vec<u64>,
    res: &mut Vec<u64>,
    informative_coverage: i64,
    single_haplotype_coverage: i64,
    het_cov: i64,
    hom_cov: i64,
    n_hap: i64,
) {
    if hap.length == 0 {
        return;
    }
    let query_levels = reads.read_quality(
        query_id as usize,
        0,
        reads.read_length(query_id as usize) as i64,
        false,
    );

    let real = |stat: &SnpStats| -> bool {
        !(stat.occ_0 < 2 || stat.occ_1 < 2)
            && !is_standalone(stat, standard_rate, standard_maximum)
            && i64::from(stat.occ_0) >= single_haplotype_coverage
            && i64::from(stat.occ_1) >= informative_coverage
    };

    idx.clear();
    res.clear();
    res.resize(overlaps.len(), u64::MAX);

    // First pass: keep the sites whose alleles carry enough high-quality
    // observations, folding the evidence into the surviving statistics.
    let snp_count = hap.snp_stat.len();
    let mut m_snp_stat = 0usize;
    let mut m_list = 0usize;
    let mut i = 0usize;
    let mut k = 1usize;
    let mut l = 0usize;
    while k <= snp_count {
        if k == snp_count || hap.snp_stat[k].site != hap.snp_stat[l].site {
            while i < hap.length && hap.list[i].site != hap.snp_stat[l].site {
                i += 1;
            }
            assert!(i < hap.length && hap.list[i].site == hap.snp_stat[l].site);
            let mut kept = 0usize;
            for o in l..k {
                if real(&hap.snp_stat[o]) {
                    kept += 1;
                }
            }
            if kept > 0 {
                let mut hq = [0u32; 5];
                let mut hp = [u32::MAX; 4];
                let i_start = i;
                let mut c0 = 0u64;
                let mut rev_n = 0u64;
                let mut i0 = i;
                while i0 < hap.length && hap.list[i0].site == hap.snp_stat[l].site {
                    hap.list[i0].site = u32::MAX;
                    let record = hap.list[i0];
                    if !record.is_mismatch() {
                        c0 += 1;
                        if overlaps[record.overlap_identifier as usize].y_strand == 0 {
                            rev_n += 1;
                        }
                        if record.base_quality() > 1 || record.window_quality() > 1 {
                            hq[4] += 1;
                        }
                    } else {
                        let base = SEQUENCE_NT6_TABLE[record.mismatch_base as usize] as usize;
                        hp[base] = record.overlap_site;
                        if record.base_quality() > 1 || record.window_quality() > 1 {
                            hq[base] += 1;
                        }
                    }
                    i0 += 1;
                }
                assert!(i0 > i_start);
                for base in 0..4 {
                    if hp[base] != u32::MAX && hq[base] < 2 {
                        hap.snp_stat[hp[base] as usize].occ_1 = 0;
                    }
                }
                let (tqual, wqual) =
                    window_quality_levels(&query_levels, hap.snp_stat[l].site as u64, QUALITY_WINDOW);
                if tqual > 1 || wqual > 1 {
                    hq[4] += 1;
                }

                if hq[4] >= 2 && (hq[0] >= 2 || hq[1] >= 2 || hq[2] >= 2 || hq[3] >= 2) {
                    let m_snp_stat0 = m_snp_stat;
                    for o in l..k {
                        let stat = hap.snp_stat[o];
                        if !real(&stat) {
                            continue;
                        }
                        hap.snp_stat[m_snp_stat] = stat;
                        let mut c1 = 0u64;
                        for ii in i_start..i0 {
                            let record = hap.list[ii];
                            if (record.overlap_site == o as u32 && record.site == u32::MAX)
                                || !record.is_mismatch()
                            {
                                hap.list[ii].overlap_site = m_snp_stat as u32;
                                hap.list[ii].site = hap.snp_stat[o].site;
                                if record.is_mismatch() {
                                    c1 += 1;
                                }
                            }
                        }
                        hap.snp_stat[m_snp_stat].occ_0 = (c0 + 1) as u32;
                        hap.snp_stat[m_snp_stat].occ_1 = c1 as u32;
                        hap.snp_stat[m_snp_stat].overlap_num = (rev_n + 1) as u32;
                        hap.snp_stat[m_snp_stat].non_homopolymer_num = m_list as u32;
                        hap.snp_stat[m_snp_stat].score = -1;
                        m_snp_stat += 1;
                    }
                    let m1 = m_list;
                    for ii in i_start..i0 {
                        if hap.list[ii].site == u32::MAX {
                            hap.list[ii].overlap_site = u32::MAX;
                            hap.list[ii].site = hap.snp_stat[l].site;
                        }
                        fill_incom(
                            res,
                            hap.list[ii].overlap_identifier as u64,
                            hap.list[ii].site as u64,
                            idx,
                        );
                        hap.list[m_list] = hap.list[ii];
                        m_list += 1;
                    }
                    radix_sort_by_key_32(&mut hap.list[m1..m_list], &|record| {
                        record.overlap_identifier
                    });
                    for ii in m_snp_stat0..m_snp_stat {
                        hap.snp_stat[ii].homopolymer_num = m_list as u32;
                    }
                    i = i0;
                    idx.push(hap.snp_stat[l].site as u64);
                }
            }
            l = k;
        }
        k += 1;
    }
    hap.snp_stat.truncate(m_snp_stat);
    hap.length = m_list;
    if hap.snp_stat.is_empty() || hap.length == 0 {
        return;
    }

    // Fill the observations an overlap missed between two kept sites.
    if res.len() > overlaps.len() {
        for k in overlaps.len()..res.len() {
            let value = res[k];
            hap.add(HaplotypeEvidence {
                site: value as u32,
                overlap_identifier: (value >> 32) as u32,
                overlap_site: u32::MAX,
                evidence_type: 0,
                coverage: 1,
                mismatch_base: 0,
            });
        }
        radix_sort_by_key_32(&mut hap.list[..hap.length], &|record| record.site);

        let snp_count = hap.snp_stat.len();
        let mut i = 0usize;
        let mut k = 1usize;
        let mut l = 0usize;
        while k <= snp_count {
            if k == snp_count || hap.snp_stat[k].site != hap.snp_stat[l].site {
                while i < hap.length && hap.list[i].site != hap.snp_stat[l].site {
                    i += 1;
                }
                assert!(i < hap.length && hap.list[i].site == hap.snp_stat[l].site);
                let mut m1 = u32::MAX;
                let mut mc = 0u8;
                let i_start = i;
                let mut c0 = 0u64;
                let mut rev_n = 0u64;
                let mut i0 = i;
                while i0 < hap.length && hap.list[i0].site == hap.snp_stat[l].site {
                    let record = hap.list[i0];
                    if !record.is_mismatch() {
                        c0 += 1;
                        if overlaps[record.overlap_identifier as usize].y_strand == 0 {
                            rev_n += 1;
                        }
                        if record.overlap_site != u32::MAX {
                            m1 = record.overlap_site;
                            mc = record.mismatch_base;
                        }
                    }
                    i0 += 1;
                }
                assert!(i0 > i_start);
                assert!(m1 != u32::MAX);
                assert!(mc != 0);

                let mut sorted = true;
                for ii in i_start..i0 {
                    if !hap.list[ii].is_mismatch() {
                        hap.list[ii].overlap_site = m1;
                        hap.list[ii].mismatch_base = mc;
                    }
                    if ii > i_start
                        && hap.list[ii].overlap_identifier <= hap.list[ii - 1].overlap_identifier
                    {
                        sorted = false;
                    }
                }
                if !sorted {
                    radix_sort_by_key_32(&mut hap.list[i_start..i0], &|record| {
                        record.overlap_identifier
                    });
                }

                for o in l..k {
                    let mut c1 = 0u64;
                    for ii in i_start..i0 {
                        if hap.list[ii].overlap_site == o as u32 && hap.list[ii].is_mismatch() {
                            c1 += 1;
                        }
                    }
                    hap.snp_stat[o].occ_0 = (c0 + 1) as u32;
                    hap.snp_stat[o].occ_1 = c1 as u32;
                    hap.snp_stat[o].overlap_num = (rev_n + 1) as u32;
                    hap.snp_stat[o].non_homopolymer_num = i_start as u32;
                    hap.snp_stat[o].homopolymer_num = i0 as u32;
                }
                i = i0;
            }
            l = k;
            k += 1;
        }
    }

    dp0_single_path(
        &mut hap.snp_stat,
        &hap.list[..hap.length],
        dp,
        idx,
        res,
        het_cov,
        hom_cov,
        n_hap,
        0.7,
        6,
        &query_levels,
        single_haplotype_coverage,
        informative_coverage,
    );

    // Third pass: drop the alleles the DP rejected and compact the evidence.
    let snp_count = hap.snp_stat.len();
    let mut m_snp_stat = 0usize;
    let mut m_list = 0usize;
    let mut i = 0usize;
    let mut k = 1usize;
    let mut l = 0usize;
    while k <= snp_count {
        if k == snp_count || hap.snp_stat[k].site != hap.snp_stat[l].site {
            while i < hap.length && hap.list[i].site != hap.snp_stat[l].site {
                i += 1;
            }
            assert!(i < hap.length && hap.list[i].site == hap.snp_stat[l].site);
            let mut kept = 0usize;
            for o in l..k {
                if hap.snp_stat[o].score != -1 {
                    kept += 1;
                }
            }
            if kept > 0 {
                let i_start = i;
                let mut c0 = 0u64;
                let mut rev_n = 0u64;
                let mut i0 = i;
                while i0 < hap.length && hap.list[i0].site == hap.snp_stat[l].site {
                    hap.list[i0].site = u32::MAX;
                    let record = hap.list[i0];
                    if !record.is_mismatch() {
                        c0 += 1;
                        if overlaps[record.overlap_identifier as usize].y_strand == 0 {
                            rev_n += 1;
                        }
                    }
                    i0 += 1;
                }
                assert!(i0 > i_start);
                let m_snp_stat0 = m_snp_stat;
                for o in l..k {
                    let stat = hap.snp_stat[o];
                    if stat.score == -1 {
                        continue;
                    }
                    hap.snp_stat[m_snp_stat] = stat;
                    let mut c1 = 0u64;
                    for ii in i_start..i0 {
                        let record = hap.list[ii];
                        if (record.overlap_site == o as u32 && record.site == u32::MAX)
                            || !record.is_mismatch()
                        {
                            hap.list[ii].overlap_site = m_snp_stat as u32;
                            hap.list[ii].site = hap.snp_stat[o].site;
                            if record.is_mismatch() {
                                c1 += 1;
                            }
                        }
                    }
                    hap.snp_stat[m_snp_stat].occ_0 = (c0 + 1) as u32;
                    hap.snp_stat[m_snp_stat].occ_1 = c1 as u32;
                    hap.snp_stat[m_snp_stat].overlap_num = (rev_n + 1) as u32;
                    hap.snp_stat[m_snp_stat].non_homopolymer_num = m_list as u32;
                    hap.snp_stat[m_snp_stat].score = -1;
                    m_snp_stat += 1;
                }
                for ii in i_start..i0 {
                    if hap.list[ii].site == u32::MAX {
                        continue;
                    }
                    hap.list[m_list] = hap.list[ii];
                    m_list += 1;
                }
                for ii in m_snp_stat0..m_snp_stat {
                    hap.snp_stat[ii].homopolymer_num = m_list as u32;
                }
                i = i0;
            }
            l = k;
        }
        k += 1;
    }
    hap.snp_stat.truncate(m_snp_stat);
    hap.length = m_list;
}

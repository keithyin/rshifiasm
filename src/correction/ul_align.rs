//! The ultra-long (`--ul`) read-to-unitig alignment pass, the counterpart of
//! `ul_map_lchain`, `ul_rid_lalign_adv` and `ul_local_aln` (`anchor.cpp` /
//! `Correct.cpp`).
//!
//! A UL read is chained against the unitig minimizer index with the qgen
//! chaining machinery, then every chain is refined against the unitig
//! sequence. The first pass keeps the chains whose refined windows cover the
//! overlap; the second pass splits each overlap into the HiFi reads of the
//! unitig and turns those sub-overlaps into the blocks of `all_ul_t`.

use crate::correction::hc_align::{align_ul_ed_post_extz, BitExtz, HcAlignParameters};
use crate::correction::target::TargetSource;
use crate::correction::window::WindowCigar;
use crate::counting::minimizer::sketch;
use crate::counting::position_index::PositionIndex;
use crate::overlap::lchain_parameters::set_lchain_parameters;
use crate::overlap::qgen_chain::{lchain_qgen_mcopy_fast, QgenChainParameters};
use crate::overlap::qgen_seeds::collect_qgen_seeds;
use crate::overlap::region::{CandidatesList, OverlapRegion};

/// `WINDOW` of the UL chaining (`Hash_Table.h` 9).
pub const UL_WINDOW: i64 = 375;
/// `THRESHOLD_MAX_SIZE` (`Hash_Table.h` 24).
const THRESHOLD_MAXIMUM_SIZE: f64 = 31.0;
/// `G_CHAIN_INDEL` (`inter.h` 27).
pub(crate) const CHAIN_INDEL: i64 = 128;
/// `OVERLAP_THRESHOLD_NOSI_FILTER` (`Hash_Table.h` 20).
const OVERLAP_THRESHOLD: f64 = 0.7;

/// Parameters of one UL alignment pass, the slice of `mg_idxopt_t` and
/// `asm_opt` the pass consumes.
#[derive(Debug, Clone, Copy)]
pub struct UlAlignParameters {
    /// `opt->k`, the unitig index k-mer length.
    pub kmer_length: u32,
    /// `opt->w`.
    pub minimizer_window: u32,
    /// `opt->diff_ec_ul`.
    pub error_rate: f64,
    /// `opt->bw_thres`.
    pub band_width_threshold: f64,
    /// `opt->max_n_chain`.
    pub maximum_chain_count: u32,
    /// `opt->is_HPC`.
    pub homopolymer_compression: bool,
    /// `asm_opt.mz_sample_dist`.
    pub sample_distance: u32,
    /// `asm_opt.mz_rewin`.
    pub rewindow_size: u32,
}

impl UlAlignParameters {
    /// `MIN(THRESHOLD_MAX_SIZE / diff_ec_ul, WINDOW)`.
    pub fn window_length(&self) -> i64 {
        (THRESHOLD_MAXIMUM_SIZE / self.error_rate).min(UL_WINDOW as f64) as i64
    }
}

/// The result of [`ul_map_lchain`]: the chains plus the anchor list their
/// refinement consumes.
#[derive(Debug, Clone, Default)]
pub struct UlMapResult {
    pub overlaps: Vec<OverlapRegion>,
    pub candidates: CandidatesList,
}

/// `ul_map_lchain` (`anchor.cpp` 2287) with `apend_be = 1`, `gen_off = 0`,
/// `is_accurate = 1` and the UL copy-number settings of
/// `worker_for_ul_scall_alignment`.
pub fn ul_map_lchain(
    query: &[u8],
    index: &PositionIndex,
    targets: &TargetSource<'_>,
    parameters: &UlAlignParameters,
) -> UlMapResult {
    let minimizers = sketch(
        query,
        parameters.minimizer_window as i32,
        parameters.kmer_length as i32,
        u32::MAX,
        parameters.homopolymer_compression,
        None,
        parameters.sample_distance as i32,
        parameters.rewindow_size as i32,
    )
    .minimizers;

    let target_lengths: Vec<u64> = (0..target_count(targets))
        .map(|identifier| targets.length(identifier) as u64)
        .collect();
    // `high_occ = 2`, `low_occ = NULL`.
    let seeds = collect_qgen_seeds(&minimizers, index, &target_lengths, 2, 0);

    let lchain = set_lchain_parameters(true, parameters.kmer_length);
    let chain_parameters = QgenChainParameters {
        maximum_chain_count: parameters.maximum_chain_count,
        maximum_skip: lchain.maximum_skip,
        maximum_iterations: lchain.maximum_iterations,
        maximum_distance: lchain.maximum_distance,
        gap_penalty: lchain.gap_penalty,
        skip_penalty: lchain.skip_penalty,
        band_width_rate: parameters.band_width_threshold,
        quick_check: lchain.quick_check != 0,
        append_both_ends: true,
        generate_offset: false,
        copy_number: 1,
        copy_rate: 0.2,
        chain_cutoff: 2,
        copy_kmer_cutoff: 3,
        coverage_window: 0,
    };
    let result = lchain_qgen_mcopy_fast(seeds, u32::MAX, query.len() as u64, &target_lengths, &chain_parameters);

    let mut candidates = CandidatesList {
        list: result.candidates,
        length: result.candidate_length as i64,
        size: result.candidate_length as i64,
        chain_data: Default::default(),
    };
    candidates.length = result.candidate_length as i64;
    UlMapResult {
        overlaps: result.overlaps,
        candidates,
    }
}

fn target_count(targets: &TargetSource<'_>) -> usize {
    match targets {
        TargetSource::Reads(reads) => reads.total_reads as usize,
        TargetSource::Unitigs { unitigs, .. } => unitigs.len(),
    }
}

/// `win_cluster_fliter` (`Correct.cpp` 24774): whether some window run of at
/// least `min_ovlp` query bases overlaps the aligned windows by `o_rate`.
pub fn win_cluster_filter(
    windows: &crate::correction::window::WindowCigar,
    min_ovlp: u64,
    o_rate: f64,
) -> bool {
    let list = &windows.windows;
    let count = list.len();
    let mut index = 0usize;
    while index < count && list[index].y_end == -1 {
        index += 1;
    }
    if index >= count {
        return false;
    }
    let mut window_start = list[index].x_start as u64;

    let mut aligned = 0u64;
    let mut skip = index;
    let mut window_end = window_start;
    while index < count {
        window_end = list[index].x_end as u64 + 1;
        if list[index].y_end != -1 {
            aligned += (list[index].x_end - list[index].x_start) as u64;
        }
        if window_end >= window_start + min_ovlp {
            break;
        }
        index += 1;
    }
    if aligned as f64 >= (window_end.saturating_sub(window_start) as f64) * o_rate
        && window_end.saturating_sub(window_start) >= min_ovlp
    {
        return true;
    }

    index += 1;
    while index < count {
        let window_end = list[index].x_end as u64 + 1;
        if list[index].y_end != -1 {
            aligned += (list[index].x_end - list[index].x_start) as u64;
        }
        while skip < count && (list[skip].x_start as u64) + min_ovlp < window_end {
            if list[skip].y_end != -1 {
                aligned -= (list[skip].x_end - list[skip].x_start) as u64;
            }
            skip += 1;
        }
        window_start = list[skip].x_start as u64;
        if aligned as f64 >= (window_end.saturating_sub(window_start) as f64) * o_rate
            && window_end.saturating_sub(window_start) >= min_ovlp
        {
            return true;
        }
        if skip > 0 && window_end.saturating_sub(window_start) < min_ovlp {
            skip -= 1;
            if list[skip].y_end != -1 {
                aligned += (list[skip].x_end - list[skip].x_start) as u64;
            }
            window_start = list[skip].x_start as u64;
            if aligned as f64 >= (window_end.saturating_sub(window_start) as f64) * o_rate
                && window_end.saturating_sub(window_start) >= min_ovlp
            {
                return true;
            }
        }
        index += 1;
    }
    false
}

/// The subset of `ul_ovec_buf_t` the alignment pass keeps between reads.
#[derive(Default)]
pub struct UlAlignScratch {
    pub extz: BitExtz,
    /// `ovlp_read`: the assembled target sequence reused by the second pass.
    pub target_sequence: Vec<u8>,
    /// `gen_aux_ovlp`: the scratch overlap the second pass writes.
    pub auxiliary: OverlapRegion,
    /// Reused window/cigar buffer of `ul_local_aln`.
    pub windows: WindowCigar,
}

/// `ul_local_aln` (`Correct.cpp` 24826) in its first (no sub-overlap) mode.
///
/// Returns whether the chain survives. The second mode, which splits the
/// overlap into sub-overlaps and links them into `ul_ov_t` chains, is not
/// ported yet.
pub fn ul_local_aln_first(
    region: &mut OverlapRegion,
    scratch: &mut UlAlignScratch,
    targets: &TargetSource<'_>,
    query: &[u8],
    parameters: &UlAlignParameters,
    minimum_overlap: u64,
) -> bool {
    let window_length = parameters.window_length();
    align_ul_ed_post_extz(
        region,
        &mut scratch.extz,
        targets,
        query,
        &HcAlignParameters {
            error_rate: parameters.error_rate,
            window_length,
            overlap_cutoff: -1.0,
            force_alignment: true,
        },
    );

    let overlap_length =
        i64::from(region.x_position_end) + 1 - i64::from(region.x_position_start);
    let aligned_length = i64::from(region.alignment_length);
    if aligned_length as f64 >= overlap_length as f64 * OVERLAP_THRESHOLD {
        return true;
    }
    if aligned_length <= minimum_overlap as i64 {
        return false;
    }
    win_cluster_filter(
        &region.window_cigar,
        minimum_overlap,
        parameters.error_rate / 2.0,
    )
}

/// `ul_rid_lalign_adv` (`Correct.cpp` 25083) in its first pass (`aln == NULL`).
///
/// Every overlap is refined and the survivors are compacted to the front,
/// with the same in-place swap order as the reference. Returns the new logical
/// length.
pub fn ul_rid_lalign_adv_first(
    overlaps: &mut [OverlapRegion],
    length: usize,
    scratch: &mut UlAlignScratch,
    targets: &TargetSource<'_>,
    query: &[u8],
    parameters: &UlAlignParameters,
    minimum_overlap: u64,
) -> usize {
    let mut kept = 0usize;
    for index in 0..length {
        overlaps[index].shared_seed = overlaps[index].non_homopolymer_errors as i32;
        if !ul_local_aln_first(
            &mut overlaps[index],
            scratch,
            targets,
            query,
            parameters,
            minimum_overlap,
        ) {
            continue;
        }
        if kept != index {
            overlaps.swap(kept, index);
        }
        overlaps[kept].is_match = 1;
        kept += 1;
    }
    kept
}

/// `get_chain_x_by_y` (`Correct.cpp` 22520): the query coordinate of target
/// position `q`, interpolated linearly between the two fake-cigar gap sites
/// that bracket it.
pub fn get_chain_x_by_y(region: &OverlapRegion, query: i64) -> i64 {
    let cigar = &region.fake_cigar;
    let last = cigar.len() - 1;
    let mut x = cigar.gap_position(last);
    let mut offset = cigar.gap_shift(last);
    let mut y = x - i64::from(region.x_position_start)
        + i64::from(region.y_position_start)
        + offset;
    if y == query {
        return x;
    }

    let mut left_x = -1i64;
    let mut left_y = -1i64;
    let mut index = 0usize;
    while index < cigar.len() {
        x = cigar.gap_position(index);
        offset = cigar.gap_shift(index);
        y = x - i64::from(region.x_position_start)
            + i64::from(region.y_position_start)
            + offset;
        if query < y {
            left_x = x;
            left_y = y;
            break;
        }
        index += 1;
    }
    debug_assert!(index != 0 && index != cigar.len());

    x = cigar.gap_position(index - 1);
    offset = cigar.gap_shift(index - 1);
    y = x - i64::from(region.x_position_start)
        + i64::from(region.y_position_start)
        + offset;
    let mut result =
        (((query - y) as f64) / ((left_y - y) as f64)) * ((left_x - x) as f64) + x as f64;
    if result < f64::from(region.x_position_start) {
        result = f64::from(region.x_position_start);
    }
    if result > f64::from(region.x_position_end) {
        result = f64::from(region.x_position_end);
    }
    result as i64
}

/// `extract_subov` (`Correct.cpp` 22896): clips one unitig read interval to
/// the overlap and checks that enough aligned windows cover it.
///
/// `window_hint` is the reference's `in_k` cursor and is updated in place.
pub fn extract_subov(
    target_start: i64,
    target_end: i64,
    region: &OverlapRegion,
    overlap_rate: f64,
    window_hint: Option<&mut i64>,
) -> Option<crate::assembly::ul::UlOv> {
    let mut target_start = target_start;
    let mut target_end = target_end;
    if target_start < i64::from(region.y_position_start) {
        target_start = i64::from(region.y_position_start);
    }
    if target_end > i64::from(region.y_position_end) + 1 {
        target_end = i64::from(region.y_position_end) + 1;
    }
    if target_start >= target_end {
        return None;
    }
    let query_start = get_chain_x_by_y(region, target_start);
    let query_end = get_chain_x_by_y(region, target_end - 1) + 1;
    debug_assert!(query_start < query_end);

    let mut result = crate::assembly::ul::UlOv {
        qs: query_start as u32,
        qe: query_end as u32,
        ts: target_start as u32,
        te: target_end as u32,
        rev: region.y_strand != 0,
        ..crate::assembly::ul::UlOv::default()
    };

    if overlap_rate >= 0.0 {
        let windows = &region.window_cigar.windows;
        if windows.is_empty() {
            return None;
        }
        let mut index = window_hint.as_deref().copied().unwrap_or(0);
        if index < 0 {
            index = 0;
        }
        if index >= windows.len() as i64 {
            index = windows.len() as i64 - 1;
        }
        while index < windows.len() as i64
            && query_start > i64::from(windows[index as usize].x_end)
        {
            index += 1;
        }
        if index < 0 {
            index = 0;
        }
        if index >= windows.len() as i64 {
            index = windows.len() as i64 - 1;
        }
        while index >= 0 && query_start < i64::from(windows[index as usize].x_start) {
            index -= 1;
        }
        if index < 0 {
            index = 0;
        }
        if let Some(hint) = window_hint {
            *hint = index;
        }

        let mut aligned = 0i64;
        while index < windows.len() as i64 {
            let window = windows[index as usize];
            if i64::from(window.x_start) >= query_end {
                break;
            }
            if window.y_end == -1
                || crate::correction::window_chain::is_unaligned_window(&window)
            {
                index += 1;
                continue;
            }
            let start = query_start.max(i64::from(window.x_start));
            let end = query_end.min(i64::from(window.x_end) + 1);
            if end > start {
                aligned += end - start;
            }
            index += 1;
        }
        if (query_end - query_start) as f64 * overlap_rate <= aligned as f64 && aligned > 0 {
            return Some(result);
        }
        return None;
    }
    Some(result)
}

/// `gen_sub_ov_adv` (`Correct.cpp` 22940): maps the unitig interval of an
/// overlap onto the reads and contained reads of the unitig.
#[allow(clippy::too_many_arguments)]
pub fn gen_sub_ov_adv(
    region: &OverlapRegion,
    overlap_rate: f64,
    unitig: &crate::assembly::asg::MaUtg,
    contain: &crate::assembly::ul_index::UlContain,
    contained: bool,
    read_lengths: &[u64],
    output: &mut Vec<crate::assembly::ul::UlOv>,
) -> usize {
    let before = output.len();
    let reverse = region.y_strand != 0;
    let unitig_length = i64::from(unitig.length);
    let (target_start, target_end) = if !reverse {
        (
            i64::from(region.y_position_start),
            i64::from(region.y_position_end) + 1,
        )
    } else {
        (
            unitig_length - (i64::from(region.y_position_end) + 1),
            unitig_length - i64::from(region.y_position_start),
        )
    };
    let mut window_hint = if !reverse {
        0i64
    } else {
        (region.window_cigar.windows.len() as i64 - 1).max(0)
    };

    if !contained {
        let mut position = 0i64;
        for (index, &packed) in unitig.reads.iter().enumerate() {
            let read_identifier = (packed >> 33) as usize;
            let start = position;
            let end = position + read_lengths[read_identifier] as i64;
            position += packed as u32 as i64;
            if end <= target_start {
                continue;
            }
            if start >= target_end {
                break;
            }
            let (clip_start, clip_end) = if reverse {
                (unitig_length - end, unitig_length - start)
            } else {
                (start, end)
            };
            if let Some(mut overlap) = extract_subov(
                clip_start,
                clip_end,
                region,
                overlap_rate,
                Some(&mut window_hint),
            ) {
                overlap.rev = reverse;
                overlap.el = false;
                overlap.tn = read_identifier as u32;
                overlap.qn = index as u32;
                output.push(overlap);
            }
        }
    } else {
        let count = contain.idx[region.y_identifier as usize] as u32 as usize;
        let begin = (contain.idx[region.y_identifier as usize] >> 32) as usize;
        for index in 0..count {
            let entry = contain.rids[begin + index];
            let read_identifier = (entry.x >> 1) as usize;
            let start = i64::from(entry.s);
            let end = i64::from(entry.e);
            if end <= target_start {
                continue;
            }
            if start >= target_end {
                break;
            }
            let (clip_start, clip_end) = if reverse {
                (unitig_length - end, unitig_length - start)
            } else {
                (start, end)
            };
            if let Some(mut overlap) = extract_subov(
                clip_start,
                clip_end,
                region,
                overlap_rate,
                Some(&mut window_hint),
            ) {
                overlap.rev = reverse;
                overlap.el = true;
                overlap.tn = read_identifier as u32;
                overlap.qn = index as u32;
                output.push(overlap);
            }
        }
    }
    output.len() - before
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::asg::MaUtg;
    use crate::correction::window::{WindowCigar, WindowList};
    use crate::reads::ReadStore;

    #[test]
    fn cluster_filter_requires_a_covered_run() {
        let mut windows = WindowCigar::default();
        for (start, end) in [(0i32, 99i32), (100, 199)] {
            windows.windows.push(WindowList {
                x_start: start,
                x_end: end,
                y_start: 0,
                y_end: 99,
                ..WindowList::default()
            });
        }
        assert!(win_cluster_filter(&windows, 50, 0.5));
        let empty = WindowCigar::default();
        assert!(!win_cluster_filter(&empty, 50, 0.5));
    }

    #[test]
    fn unitig_reads_are_materialised_from_the_read_store() {
        let unitig = MaUtg {
            length: 4,
            reads: vec![(0u64 << 33) | 4],
            ..MaUtg::default()
        };
        let mut reads = ReadStore::default();
        reads.read_lengths = vec![4];
        reads.sequences = vec![crate::sequence::encoding::PackedSequence::from_ascii(b"ACGT")];
        reads.total_reads = 1;
        reads.total_bases = 4;
        let store = crate::assembly::ul_index::materialize_unitig_reads(&[unitig], &reads);
        assert_eq!(store.read_lengths, vec![4]);
        assert_eq!(store.read_sequence(0), b"ACGT");
    }
}

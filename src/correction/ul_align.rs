//! The ultra-long (`--ul`) read-to-unitig alignment pass, the counterpart of
//! `ul_map_lchain`, `ul_rid_lalign_adv` and `ul_local_aln` (`anchor.cpp` /
//! `Correct.cpp`).
//!
//! A UL read is chained against the unitig minimizer index with the qgen
//! chaining machinery, then every chain is refined against the unitig
//! sequence. The first pass keeps the chains whose refined windows cover the
//! overlap; the second pass splits each overlap into the HiFi reads of the
//! unitig and turns those sub-overlaps into the blocks of `all_ul_t`.

use crate::alignment::edit_distance::TraceCigar;
use crate::assembly::asg::MaUtg;
use crate::assembly::ul::{UlOv, MAX_SEC_ERR};
use crate::correction::cal_exz::{
    adjust_ext_offset_fixed_t, cal_exact_simi_exz, cal_exz_extension_0, cal_exz_extension_1,
    cal_exz_global, cal_exz_infi_simi_adv, cal_exz_semi, scale_ed_thre, ExzAdvState, MODE_GLOBAL,
    MODE_FORWARD_EXTENSION, MODE_BACKWARD_EXTENSION,
};
use crate::correction::hc_align::{align_ul_ed_post_extz, BitExtz, HcAlignParameters};
use crate::correction::target::TargetSource;
use crate::correction::window::{WindowCigar, WindowList};
use crate::correction::window_chain::is_unaligned_window;
use crate::counting::minimizer::sketch;
use crate::counting::position_index::PositionIndex;
use crate::overlap::chaining::chain_refine;
use crate::overlap::lchain_parameters::set_lchain_parameters;
use crate::overlap::qgen_chain::{lchain_qgen_mcopy_fast, QgenChainParameters};
use crate::overlap::qgen_seeds::collect_qgen_seeds;
use crate::overlap::region::{CandidatesList, KmerHit, OverlapRegion};

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

    // the scall worker passes `is_accurate = 0`
    let lchain = set_lchain_parameters(false, parameters.kmer_length);
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
        // the scall worker passes `gen_off = 1`
        generate_offset: true,
        copy_number: 1,
        copy_rate: 0.2,
        chain_cutoff: 2,
        copy_kmer_cutoff: 3,
        // the non-fast `lchain_qgen_mcopy` passes `ocv_w = 1`: per-base
        // coverage windows for the type-3 threshold filter
        coverage_window: 1,
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

// ---------------------------------------------------------------------------
// The second-pass alignment cluster (`Correct.cpp` 15838-25134). The UL scall
// driver runs `ul_rid_lalign_adv` twice: the first pass only filters the
// chains, the second splits every surviving overlap into the reads and
// contained reads of its unitig and re-aligns each piece from the stored
// windows (`gen_raln`), or from scratch through the chain brackets
// (`cigar_gen_by_chain_adv_local`).
// ---------------------------------------------------------------------------

/// `MAX_CNS_L` (`Levenshtein_distance.h` 751).
pub(crate) const MAX_CNS_LENGTH: i64 = 3072;
/// `MAX_CNS_E` (`Levenshtein_distance.h` 752).
pub(crate) const MAX_CNS_ERROR: i64 = 1023;
/// `FORCE_CNS_L` (`Levenshtein_distance.h` 753).
pub(crate) const FORCED_CNS_LENGTH: i64 = 256;
/// `r_simi_w` (`Correct.cpp` 22).
const R_SIMI_W: f64 = 0.05;
/// The `alnw` sentinel of `readID`: a hit that no longer belongs to a chain.
pub(crate) const ALNW_SENTINEL: u32 = 0x7fff_ffff;

/// `gen_err_unaligned` (`Correct.cpp` 17978): the error charged to a span no
/// alignment covers. The `>` comparison of the second branch happens in
/// `double`, so the result converts through `double` before truncation.
pub(crate) fn gen_err_unaligned(query_length: i64, target_length: i64) -> i64 {
    const FORCE_SIN_L: i64 = 512;
    if query_length <= FORCE_SIN_L {
        return query_length.max(target_length);
    }
    let plain = query_length.min(target_length) as f64;
    let scaled = query_length as f64 * 0.51 + 1.0;
    (if plain > scaled { plain } else { scaled }) as i64
}

/// `set_bit_extz_t` (`Correct.h` 1375): a read-only view of one stored window
/// alignment. The reference points `cigar.a` into the shared cell buffer; the
/// port snapshots the window's cell slice instead, which every consumer only
/// reads.
pub(crate) struct WindowExtz<'a> {
    /// `ts`
    pub text_start: i32,
    /// `te`
    pub text_end: i32,
    /// `ps`
    pub pattern_start: i32,
    /// `pe`
    pub pattern_end: i32,
    /// `err`
    pub error: i32,
    /// `thre`
    pub threshold: i32,
    /// The window's packed cigar cells, `cigar.a[0..cigar.n)`.
    pub cells: &'a [u16],
}

impl<'a> WindowExtz<'a> {
    pub fn new(window: &'a WindowList, cells: &'a [u16]) -> Self {
        let start = window.cigar_start as usize;
        let end = start + window.cigar_length as usize;
        Self {
            text_start: window.x_start,
            text_end: window.x_end,
            pattern_start: window.y_start,
            pattern_end: window.y_end,
            error: i32::from(window.error),
            threshold: i32::from(window.error),
            cells: &cells[start..end],
        }
    }

    /// The standalone-`bit_extz_t` view of a freshly aligned `BitExtz`.
    pub fn of_extz(extz: &'a BitExtz) -> Self {
        Self {
            text_start: extz.result.text_start,
            text_end: extz.result.text_end,
            pattern_start: extz.result.pattern_start,
            pattern_end: extz.result.pattern_end,
            error: extz.result.error,
            threshold: extz.result.threshold,
            cells: &extz.cigar.packed,
        }
    }

    /// `is_align(ez)`
    pub fn is_aligned(&self) -> bool {
        self.error <= self.threshold
    }

    /// `pop_trace` (`Levenshtein_distance.h` 533).
    pub fn pop_trace(&self, mut index: usize) -> (usize, u16, u32) {
        let code = self.cells[index] >> 14;
        let mut length = u32::from(self.cells[index] & 0x3fff);
        index += 1;
        while index < self.cells.len() && code == self.cells[index] >> 14 {
            length += u32::from(self.cells[index] & 0x3fff);
            index += 1;
        }
        (index, code, length)
    }

    /// `pop_trace_back` (`Levenshtein_distance.h` 542).
    pub fn pop_trace_back(&self, mut index: i64) -> (i64, u16, u32) {
        let code = self.cells[index as usize] >> 14;
        let mut length = u32::from(self.cells[index as usize] & 0x3fff);
        index -= 1;
        while index >= 0 && code == self.cells[index as usize] >> 14 {
            length += u32::from(self.cells[index as usize] & 0x3fff);
            index -= 1;
        }
        (index, code, length)
    }
}

/// `rtrace_t` (`Correct.h` 1111): the trace index of one sub-overlap
/// alignment produced by `gen_raln`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Rtrace {
    /// `[c_qs, c_qe)` on the query and `[c_ts, c_te)` on the target.
    pub c_qs: i32,
    pub c_qe: i32,
    pub c_ts: i32,
    pub c_te: i32,
    /// `[c_wsid, c_weid]` window range and `[c_wsii, c_weii)` cigar range.
    pub c_wsid: i32,
    pub c_weid: i32,
    pub c_wsii: i32,
    pub c_weii: i32,
    pub rev: i32,
    pub sfx_e: i32,
    pub pfx_e: i32,
    pub mid_e: i32,
    pub oid: i32,
}

/// `rtrace_iter` (`Correct.h` 1386): the backwards error iterator of
/// `get_rid_backward_cigar_err`. A fresh iterator carries `k == INT32_MAX`,
/// which makes the first call initialise the rest.
#[derive(Clone, Copy, Debug)]
pub(crate) struct RtraceIter {
    pub k: i64,
    pub q: [i64; 2],
    pub t: [i64; 2],
    pub cq: [i64; 2],
    pub ct: [i64; 2],
    pub ci: [i64; 2],
    pub werr: i64,
    pub werr0: i64,
    pub cerr: i64,
    pub qoff: i64,
    pub f: i64,
    pub toff: i64,
    pub coff: i64,
    pub cur_qoff: i64,
}

impl Default for RtraceIter {
    fn default() -> Self {
        Self {
            k: i32::MAX as i64,
            q: [0; 2],
            t: [0; 2],
            cq: [0; 2],
            ct: [0; 2],
            ci: [0; 2],
            werr: 0,
            werr0: 0,
            cerr: 0,
            qoff: 0,
            f: 0,
            toff: 0,
            coff: 0,
            cur_qoff: 0,
        }
    }
}

/// `update_trace_idx` (`Correct.cpp` 23375): widens the trace index to cover
/// one more window.
#[allow(clippy::too_many_arguments)]
pub(crate) fn update_trace_idx(
    trace: &mut Rtrace,
    window_id: i64,
    window_cigar_start: i64,
    window_cigar_end: i64,
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
) {
    if query_start < i64::from(trace.c_qs) || target_start < i64::from(trace.c_ts) {
        trace.c_qs = query_start as i32;
        trace.c_ts = target_start as i32;
        trace.c_wsid = window_id as i32;
        trace.c_wsii = window_cigar_start as i32;
    }
    if query_end > i64::from(trace.c_qe) || target_end > i64::from(trace.c_te) {
        trace.c_qe = query_end as i32;
        trace.c_te = target_end as i32;
        trace.c_weid = window_id as i32;
        trace.c_weii = window_cigar_end as i32;
    }
}

/// `win_boundary_offset` (`Correct.cpp` 23073): maps one target coordinate
/// onto the query inside, before, or after window `wi`.
pub(crate) fn win_boundary_offset(
    windows: &[WindowList],
    window_index: usize,
    target_offset: i64,
    query_length: i64,
    query_offset: &mut i64,
    target_out: &mut i64,
) {
    let mut query = [-1i64; 2];
    let mut target = [-1i64; 2];
    let window = &windows[window_index];
    if target_offset >= i64::from(window.y_start) && target_offset <= i64::from(window.y_end) {
        query[0] = i64::from(window.x_start);
        query[1] = i64::from(window.x_end);
        target[0] = i64::from(window.y_start);
        target[1] = i64::from(window.y_end);
    } else if target_offset < i64::from(window.y_start) {
        if window_index > 0 {
            let previous = &windows[window_index - 1];
            query[0] = i64::from(previous.x_end) + 1;
            query[1] = i64::from(window.x_start) - 1;
            target[0] = i64::from(previous.y_end) + 1;
            target[1] = i64::from(window.y_start) - 1;
        } else {
            let query_distance = i64::from(window.x_start);
            let target_distance = i64::from(window.y_start) - target_offset;
            let distance = query_distance.min(target_distance);
            *query_offset = i64::from(window.x_start) - distance;
            *target_out = i64::from(window.y_start) - distance;
        }
    } else if target_offset > i64::from(window.y_end) {
        if window_index + 1 < windows.len() {
            let next = &windows[window_index + 1];
            query[0] = i64::from(window.x_end) + 1;
            query[1] = i64::from(next.x_start) - 1;
            target[0] = i64::from(window.y_end) + 1;
            target[1] = i64::from(next.y_start) - 1;
        } else {
            let query_distance = query_length - 1 - i64::from(window.x_end);
            let target_distance = target_offset - i64::from(window.y_end);
            let distance = query_distance.min(target_distance);
            *query_offset = i64::from(window.x_end) + distance;
            *target_out = i64::from(window.y_end) + distance;
        }
    }

    if query[0] >= 0 && query[1] >= 0 && target[0] >= 0 && target[1] >= 0 {
        *target_out = target_offset;
        *query_offset = query[0]
            + crate::assembly::ul::get_offset_adjust(
                (target_offset - target[0]) as u32,
                (target[1] - target[0]) as u32,
                (query[1] - query[0]) as u32,
            ) as i64;
    }
}

/// `hc_aln_exz_simi_adv` (`Correct.cpp` 23110): the threshold-ladder driver of
/// the UL second pass, running the no-traceback kernels unless
/// `generate_trace` asks for a cigar.
///
/// The reference derives the target identity from `uref`/`hpc_g`/`rref`; the
/// UL call sites always pass a unitig graph, so `targets` carries it and
/// `identifier`/`strand` select the unitig. `region` would only be consulted
/// for a semi-global mapping and no UL call site can reach that branch.
#[allow(clippy::too_many_arguments)]
pub(crate) fn hc_aln_exz_simi_adv(
    targets: &TargetSource,
    query: &[u8],
    identifier: usize,
    strand: u8,
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    query_minimum: i64,
    query_maximum: i64,
    target_minimum: i64,
    target_maximum: i64,
    mode: i64,
    extz: &mut BitExtz,
    error_rate: f64,
    maximum_length: i64,
    maximum_error: i64,
    forced_length: i64,
    mut estimate_error: i64,
    region: Option<&OverlapRegion>,
    generate_trace: bool,
) -> bool {
    extz.clear_align();
    extz.result.threshold = 0;
    let query_length = query_end - query_start;
    if estimate_error < 0 {
        estimate_error = (query_length as f64 * error_rate) as i64;
    }
    if query_length <= 0 {
        return false;
    }

    if query_length <= 16
        && cal_exact_simi_exz(
            targets,
            extz,
            query,
            identifier,
            strand,
            query_start,
            query_end,
            target_start,
            target_end,
            &mut ExzAdvState::default(),
            mode,
            region,
        )
    {
        return true;
    }

    if query_length <= maximum_length && (estimate_error >> 1) <= maximum_error {
        let mut state = ExzAdvState::default();
        let mut threshold = scale_ed_thre(estimate_error as u32, maximum_error as u32) as i64;
        if cal_exz_infi_simi_adv(
            targets,
            extz,
            query,
            identifier,
            strand,
            query_start,
            query_end,
            target_start,
            target_end,
            &mut state,
            threshold,
            query_minimum,
            query_maximum,
            target_minimum,
            target_maximum,
            mode,
            region,
        ) {
            return true;
        }

        let mut previous = threshold;
        threshold = (query_length as f64 * error_rate) as i64;
        threshold = scale_ed_thre(threshold as u32, maximum_error as u32) as i64;
        if threshold > previous
            && cal_exz_infi_simi_adv(
                targets,
                extz,
                query,
                identifier,
                strand,
                query_start,
                query_end,
                target_start,
                target_end,
                &mut state,
                threshold,
                query_minimum,
                query_maximum,
                target_minimum,
                target_maximum,
                mode,
                region,
            )
        {
            return true;
        }

        previous = threshold;
        threshold <<= 1;
        threshold = scale_ed_thre(threshold as u32, maximum_error as u32) as i64;
        if threshold > previous
            && cal_exz_infi_simi_adv(
                targets,
                extz,
                query,
                identifier,
                strand,
                query_start,
                query_end,
                target_start,
                target_end,
                &mut state,
                threshold,
                query_minimum,
                query_maximum,
                target_minimum,
                target_maximum,
                mode,
                region,
            )
        {
            return true;
        }

        previous = threshold;
        threshold = (query_length as f64 * 0.51) as i64;
        threshold = scale_ed_thre(threshold as u32, maximum_error as u32) as i64;
        if threshold > previous
            && cal_exz_infi_simi_adv(
                targets,
                extz,
                query,
                identifier,
                strand,
                query_start,
                query_end,
                target_start,
                target_end,
                &mut state,
                threshold,
                query_minimum,
                query_maximum,
                target_minimum,
                target_maximum,
                mode,
                region,
            )
        {
            return true;
        }

        if query_length <= forced_length {
            let threshold = maximum_error;
            if cal_exz_infi_simi_adv(
                targets,
                extz,
                query,
                identifier,
                strand,
                query_start,
                query_end,
                target_start,
                target_end,
                &mut state,
                threshold,
                query_minimum,
                query_maximum,
                target_minimum,
                target_maximum,
                mode,
                region,
            ) {
                return true;
            }
        }
    }

    false
}

/// The outputs of [`hc_aln_exz_by_exist_cigar_with_p`] besides the updated
/// `ts`/`te`: the cigar index and coordinate brackets the scan matched.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct ExistCigarBounds {
    pub error_delta: i64,
    pub cigar_start: i64,
    pub cigar_end: i64,
    pub pattern_start: i64,
    pub pattern_end: i64,
    pub text_start: i64,
    pub text_end: i64,
}

/// `hc_aln_exz_by_exist_cigar_with_p` (`Correct.cpp` 23239): re-derives the
/// error of `[ps, pe)` on the target from the stored cigar of a window view,
/// rewriting `*ts`/`*te` to the query coordinates that match.
#[allow(clippy::too_many_arguments)]
pub(crate) fn hc_aln_exz_by_exist_cigar_with_p(
    extz: &WindowExtz<'_>,
    mode: i64,
    pattern_start: i64,
    pattern_end: i64,
    text_bounds: &mut [i64; 2],
) -> ExistCigarBounds {
    assert!(i64::from(extz.pattern_start) <= pattern_start);
    assert!(i64::from(extz.pattern_end) + 1 >= pattern_end);
    assert!(mode > 0);
    let pattern_start0 = text_bounds[0];
    let pattern_end0 = text_bounds[1];
    let mut result = ExistCigarBounds::default();
    result.error_delta = 0;
    if mode == 1 {
        text_bounds[1] = -1;
    } else if mode == 2 {
        text_bounds[0] = -1;
    } else {
        text_bounds[0] = -1;
        text_bounds[1] = -1;
    }
    result.cigar_start = -1;
    result.cigar_end = -1;
    result.pattern_start = -1;
    result.pattern_end = -1;
    result.text_start = -1;
    result.text_end = -1;

    let mut pattern_index = i64::from(extz.pattern_start);
    let mut text_index = i64::from(extz.text_start);
    let mut error = [0i64; 2];
    let mut total_error = 0i64;
    let mut start_set = false;
    let mut end_set = false;
    let mut cell_index = 0usize;
    while cell_index < extz.cells.len() {
        let (next, code, length) = extz.pop_trace(cell_index);
        cell_index = next;
        let run_start = pattern_index;
        let run_pattern_start = pattern_index;
        let run_text_start = text_index;
        let length = i64::from(length);
        if code <= 1 {
            pattern_index += length;
            text_index += length;
        } else if code == 2 {
            pattern_index += length;
        } else if code == 3 {
            text_index += length;
        }
        let run_end = pattern_index;

        if run_end < pattern_start {
            if code != 0 {
                total_error += length;
            }
            continue;
        }
        if run_start > pattern_end {
            if code != 0 {
                total_error += length;
            }
            break;
        }

        if !start_set {
            if pattern_start >= run_start && pattern_start < run_end {
                if code <= 1 {
                    let text_start1 = text_index - (run_end - pattern_start);
                    if text_start1 == pattern_start0 || mode != 1 {
                        if code == 1 {
                            error[0] = total_error + (pattern_start - run_start);
                        } else {
                            error[0] = total_error;
                        }
                        start_set = true;
                        text_bounds[0] = text_start1;
                        result.cigar_start = cell_index as i64 - 1;
                        result.pattern_start = run_pattern_start;
                        result.text_start = run_text_start;
                    }
                } else if code == 2 {
                    let text_start1 = text_index;
                    if text_start1 == pattern_start0 || mode != 1 {
                        error[0] = total_error + (pattern_start - run_start);
                        start_set = true;
                        text_bounds[0] = text_start1;
                        result.cigar_start = cell_index as i64 - 1;
                        result.pattern_start = run_pattern_start;
                        result.text_start = run_text_start;
                    }
                } else if code == 3 {
                    let text_start1 = text_index;
                    if text_start1 == pattern_start0 || mode != 1 {
                        error[0] = total_error + length;
                        start_set = true;
                        text_bounds[0] = text_start1;
                        result.cigar_start = cell_index as i64 - 1;
                        result.pattern_start = run_pattern_start;
                        result.text_start = run_text_start;
                    }
                }
            }

            if code == 3 && pattern_start >= run_start && pattern_start <= run_end {
                let text_start1 = text_index - length;
                if text_start1 == pattern_start0 {
                    error[0] = total_error;
                    start_set = true;
                    text_bounds[0] = text_start1;
                    result.cigar_start = cell_index as i64 - 1;
                    result.pattern_start = run_pattern_start;
                    result.text_start = run_text_start;
                }
            }
        }

        if !end_set {
            if pattern_end > run_start && pattern_end <= run_end {
                if code <= 1 {
                    let text_end1 = text_index - (run_end - pattern_end);
                    if text_end1 == pattern_end0 || mode != 2 {
                        if code == 1 {
                            error[1] = total_error + (pattern_end - run_start);
                        } else {
                            error[1] = total_error;
                        }
                        end_set = true;
                        text_bounds[1] = text_end1;
                        result.cigar_end = cell_index as i64;
                        result.pattern_end = pattern_index;
                        result.text_end = text_index;
                    }
                } else if code == 2 {
                    let text_end1 = text_index;
                    if text_end1 == pattern_end0 || mode != 2 {
                        error[1] = total_error + (pattern_end - run_start);
                        end_set = true;
                        text_bounds[1] = text_end1;
                        result.cigar_end = cell_index as i64;
                        result.pattern_end = pattern_index;
                        result.text_end = text_index;
                    }
                } else if code == 3 {
                    let text_end1 = text_index - length;
                    if text_end1 == pattern_end0 || mode != 2 {
                        error[1] = total_error;
                        end_set = true;
                        text_bounds[1] = text_end1;
                        result.cigar_end = cell_index as i64;
                        result.pattern_end = pattern_index;
                        result.text_end = text_index;
                    }
                }
            }
            if code == 3 && pattern_end >= run_start && pattern_end <= run_end {
                let text_end1 = text_index;
                if text_end1 == pattern_end0 {
                    error[1] = total_error + length;
                    end_set = true;
                    text_bounds[1] = text_end1;
                    result.cigar_end = cell_index as i64;
                    result.pattern_end = pattern_index;
                    result.text_end = text_index;
                }
            }
            if end_set {
                break;
            }
        }

        if code != 0 {
            total_error += length;
        }
    }

    assert!(
        (mode == 1 && text_bounds[0] == pattern_start0)
            || (mode == 2 && text_bounds[1] == pattern_end0)
            || mode == 3
    );
    result.error_delta = error[1] - error[0];
    result
}

/// `extract_sub_werr` (`Correct.cpp` 23811): walks the cigar of `extz`
/// backwards from the iterator state and reports the error between
/// `spec_toff` and `tar_te`, updating the iterator so the next query of the
/// same overlap resumes from the last cell.
///
/// `extz` is either a stored window snapshot or the freshly aligned scratch
/// `bit_extz_t`; both expose the same read-only surface.
#[allow(clippy::too_many_arguments)]
pub(crate) fn extract_sub_werr(
    extz: &WindowExtz<'_>,
    cigar_start: i64,
    cigar_end: i64,
    cell_pattern_start: i64,
    cell_pattern_end: i64,
    cell_text_start: i64,
    cell_text_end: i64,
    target_pattern_start: i64,
    target_pattern_end: i64,
    target_text_start: i64,
    target_text_end: i64,
    spec_toff: i64,
    iterator: &mut RtraceIter,
) -> i64 {
    if iterator.toff >= target_pattern_end && iterator.qoff >= target_text_end {
        iterator.coff = cigar_end - 1;
        iterator.toff = cell_pattern_end;
        iterator.qoff = cell_text_end;
        iterator.cerr = 0;
    }
    if !extz.is_aligned() {
        return i32::MAX as i64;
    }
    if extz.cells.is_empty() {
        return (((target_text_end - spec_toff) as f64) / ((target_text_end - target_text_start) as f64)
            * f64::from(extz.error)) as i64;
    }

    let mut cell_index = iterator.coff;
    let mut pattern_index = iterator.toff;
    let mut text_index = iterator.qoff;
    let mut error = iterator.cerr;
    let mut result = i32::MAX as i64;
    let mut start_set = false;
    let mut end_set = false;
    while cell_index >= cigar_start {
        let cell_index0 = cell_index;
        let pattern_index0 = pattern_index;
        let text_index0 = text_index;
        let error0 = error;
        let (previous, code, length) = extz.pop_trace_back(cell_index);
        cell_index = previous;
        let length = i64::from(length);
        let run_text_end = text_index;
        if code <= 1 {
            pattern_index -= length;
            text_index -= length;
        } else if code == 2 {
            pattern_index -= length;
        } else if code == 3 {
            text_index -= length;
        }
        let run_text_start = text_index;
        let mut flag = [false; 2];
        let mut text = [run_text_start, run_text_end];

        if run_text_end < target_text_start {
            if code != 0 {
                error += length;
            }
            break;
        }
        if run_text_start > target_text_end {
            if code != 0 {
                error += length;
            }
            continue;
        }

        if !end_set {
            if target_text_end > run_text_start && target_text_end <= run_text_end {
                assert!(code != 2);
                if code <= 1 {
                    let pattern_end1 = pattern_index + (target_text_end - run_text_start);
                    if pattern_end1 == target_pattern_end {
                        if code == 1 {
                            error = target_text_end - run_text_start;
                        } else {
                            error = 0;
                        }
                        end_set = true;
                    }
                } else if code == 3 {
                    let pattern_end1 = pattern_index;
                    if pattern_end1 == target_pattern_end {
                        error = target_text_end - run_text_start;
                        end_set = true;
                    }
                }
            }

            if code == 2 && target_text_end >= run_text_start && target_text_end <= run_text_end {
                if target_pattern_end >= pattern_index && target_pattern_end < pattern_index + length {
                    error = target_pattern_end - pattern_index;
                    end_set = true;
                }
            }
            if end_set {
                flag[0] = true;
                text[1] = target_text_end;
            }
        }

        if !start_set {
            if target_text_start >= run_text_start && target_text_start < run_text_end {
                assert!(code != 2);
                if code <= 1 {
                    let pattern_start1 = pattern_index + (target_text_start - run_text_start);
                    if pattern_start1 == target_pattern_start {
                        if code == 1 {
                            error += run_text_end - target_text_start;
                        }
                        start_set = true;
                    }
                } else if code == 3 {
                    let pattern_start1 = pattern_index;
                    if pattern_start1 == target_pattern_start {
                        error += run_text_end - target_text_start;
                        start_set = true;
                    }
                }
            }

            if code == 2 && target_text_start >= run_text_start && target_text_start <= run_text_end {
                if target_pattern_start >= pattern_index && target_pattern_start < pattern_index + length {
                    error += pattern_index + length - target_pattern_start;
                    start_set = true;
                }
            }
            if start_set {
                flag[1] = true;
                text[0] = target_text_start;
            }
        }

        if spec_toff > text[1] {
            break;
        }

        if (spec_toff >= text[0] && spec_toff < text[1]) || (spec_toff == text[0] && spec_toff == text[1]) {
            result = error + i64::from(!flag[1] && !flag[0] && code != 0) * length;
            if code == 1 || code == 3 {
                result -= spec_toff - text[0];
            }
            iterator.coff = cell_index0;
            iterator.toff = pattern_index0;
            iterator.qoff = text_index0;
            iterator.cerr = error0;
        }
        if flag[1] {
            break;
        }
        if flag[0] {
            continue;
        }
        if code != 0 {
            error += length;
        }
    }
    result
}

/// `gen_clip_win_err` (`Correct.cpp` 23776): fills `extz` with the alignment
/// of one clipped window piece, keeping the previous cigar when the cached
/// boundaries match.
#[allow(clippy::too_many_arguments)]
pub(crate) fn gen_clip_win_err(
    region: &OverlapRegion,
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    error_rate: f64,
    estimated_error: i64,
    targets: &TargetSource,
    query: &[u8],
    extz: &mut BitExtz,
) {
    let query_length = query_end - query_start;
    let target_length = target_end - target_start;
    extz.result.text_start = query_start as i32;
    extz.result.text_end = (query_end - 1) as i32;
    extz.result.pattern_start = target_start as i32;
    extz.result.pattern_end = (target_end - 1) as i32;
    extz.cigar.packed.clear();
    extz.result.threshold = i32::MAX;
    extz.result.error = i32::MAX;
    if query_length == 0 && target_length == 0 {
        extz.result.error = 0;
        return;
    }
    if query_length == 0 && target_length != 0 {
        extz.result.error = target_length as i32;
        return;
    }
    if query_length != 0 && target_length == 0 {
        extz.result.error = query_length as i32;
        return;
    }
    let estimated_error = if estimated_error < 0 { -1 } else { estimated_error };
    if hc_aln_exz_simi_adv(
        targets,
        query,
        region.y_identifier as usize,
        region.y_strand as u8,
        query_start,
        query_end,
        target_start,
        target_end,
        query_start,
        query_end,
        target_start,
        target_end,
        0,
        extz,
        error_rate,
        MAX_CNS_LENGTH,
        MAX_CNS_ERROR,
        FORCED_CNS_LENGTH,
        estimated_error,
        None,
        true,
    ) {
        return;
    }
    extz.result.text_start = query_start as i32;
    extz.result.text_end = (query_end - 1) as i32;
    extz.result.pattern_start = target_start as i32;
    extz.result.pattern_end = (target_end - 1) as i32;
    extz.cigar.packed.clear();
    extz.result.threshold = i32::MAX;
    extz.result.error = gen_err_unaligned(query_length, target_length) as i32;
}

/// `get_rid_backward_cigar_err` (`Correct.cpp` 24127): the backwards error
/// iterator. Each call reports the error between `qs` and the previous query
/// point, caching the walked windows in `iterator`.
///
/// `trace`/`regions` stand in for the reference's NULL-defaulting: a missing
/// `tc` resolves through `trace[aln.qn]` and a missing region through
/// `regions[tc.oid]`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn get_rid_backward_cigar_err(
    iterator: &mut RtraceIter,
    aln: &UlOv,
    trace: Option<&[Rtrace]>,
    stored_trace: Option<&Rtrace>,
    targets: &TargetSource,
    query: &[u8],
    regions: Option<&[OverlapRegion]>,
    region: Option<&OverlapRegion>,
    extz: &mut BitExtz,
    error_rate: f64,
    qs: i64,
) -> i64 {
    let default_trace;
    let tc: &Rtrace = match stored_trace {
        Some(value) => value,
        None => {
            default_trace = trace.expect("either a trace table or a trace entry")[aln.qn as usize];
            &default_trace
        }
    };
    let default_region;
    let region: &OverlapRegion = match region {
        Some(value) => value,
        None => {
            default_region = regions.expect("either a region list or a region")
                [tc.oid as usize]
                .clone();
            &default_region
        }
    };

    if iterator.k < i64::from(tc.c_wsid) && qs > i64::from(aln.qs) {
        iterator.k = i32::MAX as i64;
    }

    if iterator.k == i32::MAX as i64 {
        iterator.k = i64::from(tc.c_weid);
        let query_end = i64::from(tc.c_qe).max(i64::from(aln.qe));
        let target_end = i64::from(tc.c_te).max(i64::from(aln.te));
        iterator.q[1] = query_end;
        iterator.t[1] = target_end;
        iterator.qoff = query_end;
        iterator.toff = target_end;
        iterator.cur_qoff = query_end;

        iterator.werr = 0;
        iterator.werr0 = 0;
        iterator.cerr = 0;
        iterator.f = 0;
        iterator.coff = i32::MAX as i64;
        extz.clear_align();
        extz.result.pattern_start = i32::MAX;
        extz.result.pattern_end = i32::MAX;
        extz.result.text_start = i32::MAX;
        extz.result.text_end = i32::MAX;
    }

    let windows = &region.window_cigar.windows;
    let cells = &region.window_cigar.cells;
    assert!(qs <= iterator.qoff);
    if qs == iterator.cur_qoff {
        return iterator.werr0;
    }

    while iterator.k >= i64::from(tc.c_wsid) {
        let window = &windows[iterator.k as usize];
        iterator.q[0] = i64::from(window.x_end) + 1;
        iterator.t[0] = i64::from(window.y_end) + 1;

        if iterator.f == 0 {
            if iterator.qoff >= iterator.q[0] {
                if iterator.q[0] <= iterator.q[1] && iterator.t[0] <= iterator.t[1] {
                    if i64::from(extz.result.pattern_start) != iterator.t[0]
                        || i64::from(extz.result.pattern_end) != iterator.t[1]
                        || i64::from(extz.result.text_start) != iterator.q[0]
                        || i64::from(extz.result.text_end) == iterator.q[1]
                    {
                        gen_clip_win_err(
                            region,
                            iterator.q[0],
                            iterator.q[1],
                            iterator.t[0],
                            iterator.t[1],
                            error_rate,
                            if iterator.k == i64::from(tc.c_weid) {
                                i64::from(tc.sfx_e)
                            } else {
                                -1
                            },
                            targets,
                            query,
                            extz,
                        );
                    }
                    let mut sub_error = i32::MAX as i64;
                    if qs >= iterator.q[0] && qs <= iterator.q[1] {
                        let view = WindowExtz::of_extz(extz);
                        sub_error = extract_sub_werr(
                            &view,
                            0,
                            extz.cigar.packed.len() as i64,
                            i64::from(extz.result.pattern_start),
                            i64::from(extz.result.pattern_end) + 1,
                            i64::from(extz.result.text_start),
                            i64::from(extz.result.text_end) + 1,
                            i64::from(extz.result.pattern_start),
                            i64::from(extz.result.pattern_end) + 1,
                            i64::from(extz.result.text_start),
                            i64::from(extz.result.text_end) + 1,
                            qs,
                            iterator,
                        );
                        if sub_error != i32::MAX as i64 {
                            if iterator.cur_qoff != qs || iterator.werr + sub_error > iterator.werr0 {
                                iterator.werr0 = iterator.werr + sub_error;
                            }
                            iterator.cur_qoff = qs;
                            if qs > iterator.q[0] && qs <= iterator.q[1] {
                                return iterator.werr0;
                            }
                        } else {
                            iterator.qoff = iterator.q[1];
                            iterator.toff = iterator.t[1];
                            iterator.cur_qoff = qs;
                            return iterator.werr0;
                        }
                    }

                    if qs <= iterator.q[0] {
                        iterator.qoff = iterator.q[0];
                        iterator.toff = iterator.t[0];
                        iterator.werr += i64::from(extz.result.error);
                        iterator.f = 1;
                        if iterator.cur_qoff != iterator.q[0] || iterator.werr > iterator.werr0 {
                            iterator.werr0 = iterator.werr;
                        }
                        iterator.cur_qoff = iterator.q[0];
                    }
                } else {
                    iterator.f = 1;
                }
            } else {
                iterator.f = 1;
            }
        }

        let window = &windows[iterator.k as usize];
        iterator.cq[0] = i64::from(window.x_start);
        iterator.q[0] = i64::from(window.x_start);
        iterator.cq[1] = i64::from(window.x_end) + 1;
        iterator.q[1] = i64::from(window.x_end) + 1;
        iterator.ct[0] = i64::from(window.y_start);
        iterator.t[0] = i64::from(window.y_start);
        iterator.ct[1] = i64::from(window.y_end) + 1;
        iterator.t[1] = i64::from(window.y_end) + 1;
        iterator.ci[0] = 0;
        iterator.ci[1] = i64::from(window.cigar_length);
        let query_window_length = iterator.q[1] - iterator.q[0];
        let target_window_length = iterator.t[1] - iterator.t[0];
        if is_unaligned_window(window) {
            if iterator.qoff >= iterator.q[0] && iterator.f == 1 {
                if qs == iterator.q[1] {
                    iterator.cur_qoff = qs;
                    return iterator.werr0;
                }
                if qs > iterator.q[0] && qs < iterator.q[1] {
                    let sub_error = (((iterator.q[1] - qs) as f64)
                        / ((iterator.q[1] - iterator.q[0]) as f64))
                        * gen_err_unaligned(query_window_length, target_window_length) as f64;
                    return iterator.werr + sub_error as i64;
                }
                iterator.qoff = iterator.q[0];
                iterator.toff = iterator.t[0];
                iterator.werr += gen_err_unaligned(query_window_length, target_window_length);
                iterator.f = 0;
                if iterator.cur_qoff != iterator.q[0] || iterator.werr > iterator.werr0 {
                    iterator.werr0 = iterator.werr;
                }
                iterator.cur_qoff = iterator.q[0];
            }
        } else {
            if iterator.k == i64::from(tc.c_wsid) {
                iterator.cq[0] = i64::from(tc.c_qs);
                iterator.ct[0] = i64::from(tc.c_ts);
                iterator.ci[0] = i64::from(tc.c_wsii);
                iterator.q[0] = i64::from(aln.qs).max(i64::from(tc.c_qs));
                iterator.t[0] = i64::from(aln.ts).max(i64::from(tc.c_ts));
            }
            if iterator.k == i64::from(tc.c_weid) {
                iterator.cq[1] = i64::from(tc.c_qe);
                iterator.ct[1] = i64::from(tc.c_te);
                iterator.ci[1] = i64::from(tc.c_weii);
                iterator.q[1] = i64::from(aln.qe).min(i64::from(tc.c_qe));
                iterator.t[1] = i64::from(aln.te).min(i64::from(tc.c_te));
            }

            if iterator.qoff >= iterator.q[0] && iterator.f == 1 {
                let mut sub_error = i32::MAX as i64;
                if qs >= iterator.q[0] && qs <= iterator.q[1] {
                    let view = WindowExtz::new(window, cells);
                    sub_error = extract_sub_werr(
                        &view,
                        iterator.ci[0],
                        iterator.ci[1],
                        iterator.ct[0],
                        iterator.ct[1],
                        iterator.cq[0],
                        iterator.cq[1],
                        iterator.t[0],
                        iterator.t[1],
                        iterator.q[0],
                        iterator.q[1],
                        qs,
                        iterator,
                    );
                    if sub_error != i32::MAX as i64 {
                        if iterator.cur_qoff != qs || iterator.werr + sub_error > iterator.werr0 {
                            iterator.werr0 = iterator.werr + sub_error;
                        }
                        iterator.cur_qoff = qs;
                        if qs > iterator.q[0] && qs <= iterator.q[1] {
                            return iterator.werr0;
                        }
                    } else {
                        iterator.qoff = iterator.q[1];
                        iterator.toff = iterator.t[1];
                        iterator.cur_qoff = qs;
                        return iterator.werr0;
                    }
                }
                if qs <= iterator.q[0] {
                    if qs == iterator.q[0] {
                        iterator.werr += sub_error;
                    } else {
                        let view = WindowExtz::new(window, cells);
                        iterator.werr += extract_sub_werr(
                            &view,
                            iterator.ci[0],
                            iterator.ci[1],
                            iterator.ct[0],
                            iterator.ct[1],
                            iterator.cq[0],
                            iterator.cq[1],
                            iterator.t[0],
                            iterator.t[1],
                            iterator.q[0],
                            iterator.q[1],
                            iterator.q[0],
                            iterator,
                        );
                    }
                    iterator.qoff = iterator.q[0];
                    iterator.toff = iterator.t[0];
                    iterator.f = 0;
                    if iterator.cur_qoff != iterator.q[0] || iterator.werr > iterator.werr0 {
                        iterator.werr0 = iterator.werr;
                    }
                    iterator.cur_qoff = iterator.q[0];
                }
            }
        }
        iterator.q[1] = iterator.q[0];
        iterator.t[1] = iterator.t[0];
        iterator.k -= 1;
    }

    iterator.q[0] = i64::from(aln.qs);
    iterator.t[0] = i64::from(aln.ts);
    if iterator.qoff >= iterator.q[0] && iterator.f == 0 {
        if iterator.q[0] <= iterator.q[1] && iterator.t[0] <= iterator.t[1] {
            if i64::from(extz.result.pattern_start) != iterator.t[0]
                || i64::from(extz.result.pattern_end) != iterator.t[1]
                || i64::from(extz.result.text_start) != iterator.q[0]
                || i64::from(extz.result.text_end) == iterator.q[1]
            {
                gen_clip_win_err(
                    region,
                    iterator.q[0],
                    iterator.q[1],
                    iterator.t[0],
                    iterator.t[1],
                    error_rate,
                    i64::from(tc.pfx_e),
                    targets,
                    query,
                    extz,
                );
            }
            let mut sub_error = i32::MAX as i64;
            if qs >= iterator.q[0] && qs <= iterator.q[1] {
                let view = WindowExtz::of_extz(extz);
                sub_error = extract_sub_werr(
                    &view,
                    0,
                    extz.cigar.packed.len() as i64,
                    i64::from(extz.result.pattern_start),
                    i64::from(extz.result.pattern_end) + 1,
                    i64::from(extz.result.text_start),
                    i64::from(extz.result.text_end) + 1,
                    i64::from(extz.result.pattern_start),
                    i64::from(extz.result.pattern_end) + 1,
                    i64::from(extz.result.text_start),
                    i64::from(extz.result.text_end) + 1,
                    qs,
                    iterator,
                );
                if sub_error != i32::MAX as i64 {
                    if iterator.cur_qoff != qs || iterator.werr + sub_error > iterator.werr0 {
                        iterator.werr0 = iterator.werr + sub_error;
                    }
                    iterator.cur_qoff = qs;
                    if qs > iterator.q[0] && qs <= iterator.q[1] {
                        return iterator.werr0;
                    }
                } else {
                    iterator.qoff = iterator.q[1];
                    iterator.toff = iterator.t[1];
                    iterator.cur_qoff = qs;
                    return iterator.werr0;
                }
            }

            if qs <= iterator.q[0] {
                iterator.qoff = iterator.q[0];
                iterator.toff = iterator.t[0];
                iterator.werr += i64::from(extz.result.error);
                iterator.f = 1;
                if iterator.cur_qoff != iterator.q[0] || iterator.werr > iterator.werr0 {
                    iterator.werr0 = iterator.werr;
                }
                iterator.cur_qoff = iterator.q[0];
            }
        }
    }
    iterator.werr
}

/// `update_ul_ov_t_coor` (`Correct.cpp` 23099): widens the result bracket.
/// Negative bounds wrap when assigned, exactly like the reference's `u32`
/// fields.
fn update_ul_ov_t_coor(result: &mut UlOv, query_begin: i64, query_end: i64, target_begin: i64, target_end: i64) {
    if i64::from(result.qs) > query_begin {
        result.qs = query_begin as u32;
    }
    if i64::from(result.qe) < query_end {
        result.qe = query_end as u32;
    }
    if i64::from(result.ts) > target_begin {
        result.ts = target_begin as u32;
    }
    if i64::from(result.te) < target_end {
        result.te = target_end as u32;
    }
}

/// Array form of [`adjust_ext_offset_fixed_t`]: the bounds are captured
/// before the mutation, matching the reference's argument evaluation.
#[allow(clippy::too_many_arguments)]
fn adjust_ext_offset_fixed_t_arr(
    query: &mut [i64; 2],
    target: &mut [i64; 2],
    query_minimum: i64,
    query_maximum: i64,
    target_minimum: i64,
    target_maximum: i64,
    threshold: i64,
    mode: i64,
) {
    let (mut query_start, mut query_end, mut target_start, mut target_end) =
        (query[0], query[1], target[0], target[1]);
    adjust_ext_offset_fixed_t(
        &mut query_start,
        &mut query_end,
        &mut target_start,
        &mut target_end,
        query_minimum,
        query_maximum,
        target_minimum,
        target_maximum,
        threshold,
        mode,
    );
    query[0] = query_start;
    query[1] = query_end;
    target[0] = target_start;
    target[1] = target_end;
}

/// `gen_raln` (`Correct.cpp` 24412): re-aligns `[ts, te)` of one overlap
/// window by window, reusing the stored cigars of fully covered windows and
/// running the similarity driver over the gaps and partial windows.
///
/// The reference's `rid` parameter feeds commented-out debug prints only.
#[allow(clippy::too_many_arguments)]
pub(crate) fn gen_raln(
    targets: &TargetSource,
    query_sequence: &[u8],
    region: &OverlapRegion,
    extz: &mut BitExtz,
    window_start_id: i64,
    window_end_id: i64,
    ts: i64,
    te: i64,
    query_length: i64,
    strand: u8,
    error_rate: f64,
    trace: &mut Rtrace,
    result: &mut UlOv,
) {
    let windows = &region.window_cigar.windows;
    *result = UlOv {
        qs: u32::MAX,
        ts: u32::MAX,
        ..UlOv::default()
    };
    *trace = Rtrace {
        c_qs: i32::MAX,
        c_ts: i32::MAX,
        c_qe: -1,
        c_te: -1,
        ..Rtrace::default()
    };

    let mut query = [-1i64; 2];
    let mut target = [ts, te];
    let mut total_error = 0i64;
    for k in window_start_id..window_end_id {
        query[1] = i64::from(windows[k as usize].x_start);
        target[1] = i64::from(windows[k as usize].y_start);
        let mut mode;
        let mut current_error;
        let mut is_aln = true;
        // the gap before the window
        if target[0] < target[1] && query[0] < query[1] {
            if query[0] < 0 {
                mode = 2;
                let (query_bound, target_bound_start, target_bound_end) =
                    (query[1], target[0], target[1]);
                adjust_ext_offset_fixed_t_arr(
                    &mut query,
                    &mut target,
                    0,
                    query_bound,
                    target_bound_start,
                    target_bound_end,
                    0,
                    mode,
                );
            } else {
                mode = 0;
            }
            let gap_query_length = query[1] - query[0];
            let gap_target_length = target[1] - target[0];
            if gap_query_length != 0 || gap_target_length != 0 {
                if gap_query_length == 0 && gap_target_length > 0 {
                    current_error = gap_target_length;
                } else if gap_target_length == 0 && gap_query_length > 0 {
                    current_error = gap_query_length;
                } else {
                    if hc_aln_exz_simi_adv(
                        targets,
                        query_sequence,
                        region.y_identifier as usize,
                        strand,
                        query[0],
                        query[1],
                        target[0],
                        target[1],
                        0,
                        query_length,
                        ts,
                        te,
                        mode,
                        extz,
                        error_rate,
                        MAX_CNS_LENGTH,
                        MAX_CNS_ERROR,
                        FORCED_CNS_LENGTH,
                        -1,
                        None,
                        false,
                    ) {
                        current_error = i64::from(extz.result.error);
                        query[0] = i64::from(extz.result.text_start);
                        query[1] = i64::from(extz.result.text_end) + 1;
                        target[0] = i64::from(extz.result.pattern_start);
                        target[1] = i64::from(extz.result.pattern_end) + 1;
                    } else {
                        current_error = gen_err_unaligned(gap_query_length, gap_target_length);
                        is_aln = false;
                    }
                }
                total_error += current_error;
                update_ul_ov_t_coor(result, query[0], query[1], target[0], target[1]);
                if k == window_start_id {
                    trace.pfx_e = (current_error * i64::from(is_aln)) as i32;
                }
            }
        }

        // inside the window
        query[0] = i64::from(windows[k as usize].x_start);
        target[0] = i64::from(windows[k as usize].y_start);
        query[1] = i64::from(windows[k as usize].x_end) + 1;
        target[1] = i64::from(windows[k as usize].y_end) + 1;
        let query_window_length = query[1] - query[0];
        let target_window_length = target[1] - target[0];
        current_error = 0;
        is_aln = true;
        if ts <= target[0] && te >= target[1] {
            if !is_unaligned_window(&windows[k as usize]) {
                current_error = i64::from(windows[k as usize].error);
            } else {
                current_error = gen_err_unaligned(query_window_length, target_window_length);
            }
            update_trace_idx(
                trace,
                k,
                0,
                i64::from(windows[k as usize].cigar_length),
                query[0],
                query[1],
                target[0],
                target[1],
            );
        } else {
            if ts > target[0] {
                target[0] = ts;
                query[0] = -1;
            }
            if te < target[1] {
                target[1] = te;
                query[1] = -1;
            }

            if query[0] != -1 && query[1] != -1 {
                mode = 0;
            } else if query[0] != -1 {
                mode = 1;
            } else if query[1] != -1 {
                mode = 2;
            } else {
                mode = 3;
            }

            if !is_unaligned_window(&windows[k as usize]) {
                let view = WindowExtz::new(&windows[k as usize], &region.window_cigar.cells);
                let bounds = hc_aln_exz_by_exist_cigar_with_p(
                    &view,
                    mode,
                    target[0],
                    target[1],
                    &mut query,
                );
                update_trace_idx(
                    trace,
                    k,
                    bounds.cigar_start,
                    bounds.cigar_end,
                    query[0],
                    query[1],
                    bounds.pattern_start,
                    bounds.pattern_end,
                );
                current_error = bounds.error_delta;
            } else {
                assert!(mode == 1 || mode == 2);
                let (query_bound_start, query_bound, target_bound_start, target_bound_end) =
                    (query[0], query[1], target[0], target[1]);
                adjust_ext_offset_fixed_t_arr(
                    &mut query,
                    &mut target,
                    query_bound_start,
                    query_bound,
                    target_bound_start,
                    target_bound_end,
                    0,
                    mode,
                );
                let clipped_query_length = query[1] - query[0];
                let clipped_target_length = target[1] - target[0];
                if clipped_query_length != 0 || clipped_target_length != 0 {
                    if clipped_query_length == 0 && clipped_target_length > 0 {
                        current_error = clipped_target_length;
                    } else if clipped_target_length == 0 && clipped_query_length > 0 {
                        current_error = clipped_query_length;
                    } else {
                        if hc_aln_exz_simi_adv(
                            targets,
                            query_sequence,
                            region.y_identifier as usize,
                            strand,
                            query[0],
                            query[1],
                            target[0],
                            target[1],
                            0,
                            query_length,
                            ts,
                            te,
                            mode,
                            extz,
                            error_rate,
                            MAX_CNS_LENGTH,
                            MAX_CNS_ERROR,
                            FORCED_CNS_LENGTH,
                            -1,
                            None,
                            false,
                        ) {
                            current_error = i64::from(extz.result.error);
                            query[0] = i64::from(extz.result.text_start);
                            query[1] = i64::from(extz.result.text_end) + 1;
                            target[0] = i64::from(extz.result.pattern_start);
                            target[1] = i64::from(extz.result.pattern_end) + 1;
                        } else {
                            current_error = gen_err_unaligned(clipped_query_length, clipped_target_length);
                            is_aln = false;
                        }
                    }
                    trace.sfx_e = (current_error * i64::from(is_aln)) as i32;
                }
            }
        }
        total_error += current_error;
        update_ul_ov_t_coor(result, query[0], query[1], target[0], target[1]);

        query[0] = query[1];
        target[0] = target[1];
    }

    // the forward extension past the last window
    if query[0] >= 0 && query[0] < query_length {
        query[1] = query_length;
        target[1] = te;
        let mut current_error = 0i64;
        let mut is_aln = true;
        let mode = 1;
        let (query_bound_start, query_bound, target_bound_start, target_bound_end) =
            (query[0], query[1], target[0], target[1]);
        adjust_ext_offset_fixed_t_arr(
            &mut query,
            &mut target,
            query_bound_start,
            query_bound,
            target_bound_start,
            target_bound_end,
            0,
            mode,
        );
        let tail_query_length = query[1] - query[0];
        let tail_target_length = target[1] - target[0];
        if tail_query_length != 0 || tail_target_length != 0 {
            if tail_query_length == 0 && tail_target_length > 0 {
                current_error = tail_target_length;
            } else if tail_target_length == 0 && tail_query_length > 0 {
                current_error = tail_query_length;
            } else if hc_aln_exz_simi_adv(
                targets,
                query_sequence,
                region.y_identifier as usize,
                strand,
                query[0],
                query[1],
                target[0],
                target[1],
                0,
                query_length,
                ts,
                te,
                mode,
                extz,
                error_rate,
                MAX_CNS_LENGTH,
                MAX_CNS_ERROR,
                FORCED_CNS_LENGTH,
                -1,
                None,
                false,
            ) {
                current_error = i64::from(extz.result.error);
                query[0] = i64::from(extz.result.text_start);
                query[1] = i64::from(extz.result.text_end) + 1;
                target[0] = i64::from(extz.result.pattern_start);
                target[1] = i64::from(extz.result.pattern_end) + 1;
            } else {
                current_error = gen_err_unaligned(tail_query_length, tail_target_length);
                is_aln = false;
            }
            update_ul_ov_t_coor(result, query[0], query[1], target[0], target[1]);
            trace.sfx_e = (current_error * i64::from(is_aln)) as i32;
        }
        total_error += current_error;
    }
    trace.mid_e = (total_error
        - i64::from(trace.pfx_e.abs())
        - i64::from(trace.sfx_e.abs())) as i32;
    result.sec = if total_error <= i64::from(MAX_SEC_ERR) {
        total_error as u32
    } else {
        MAX_SEC_ERR
    };
}

/// `extract_subov_cigar` (`Correct.cpp` 24637): finds the aligned windows
/// covering `[ts0, te0)` of the target and re-aligns that interval through
/// [`gen_raln`], accepting it when its error stays inside the similarity
/// threshold.
///
/// Returns whether the re-alignment succeeded; `result` and `trace` are only
/// meaningful on success. `window_hint` is the reference's `in_k` cursor.
#[allow(clippy::too_many_arguments)]
pub(crate) fn extract_subov_cigar(
    targets: &TargetSource,
    query: &[u8],
    region: &OverlapRegion,
    extz: &mut BitExtz,
    target_start0: i64,
    target_end0: i64,
    overlap_rate: f64,
    window_hint: Option<&mut i64>,
    query_length: i64,
    error_rate: f64,
    trace: &mut Rtrace,
    result: &mut UlOv,
) -> bool {
    let windows = &region.window_cigar.windows;
    let window_count = windows.len() as i64;
    if window_count <= 0 {
        return false;
    }
    let mut k = window_hint.as_deref().copied().unwrap_or(0);
    if k < 0 {
        k = 0;
    }
    if k >= window_count {
        k = window_count - 1;
    }
    while k < window_count && target_start0 > i64::from(windows[k as usize].y_end) {
        k += 1;
    }
    if k < 0 {
        k = 0;
    }
    if k >= window_count {
        k = window_count - 1;
    }
    while k >= 0 && target_start0 < i64::from(windows[k as usize].y_start) {
        k -= 1;
    }
    if k < 0 {
        k = 0;
    }
    if let Some(hint) = window_hint {
        *hint = k;
    }

    let mut first = i32::MAX as i64;
    let mut last = -1i64;
    let mut aligned_total = 0i64;
    while k < window_count {
        let window = &windows[k as usize];
        let window_target_start = i64::from(window.y_start);
        let window_target_end = i64::from(window.y_end) + 1;
        if window_target_start >= target_end0 {
            break;
        }
        if window.y_end == -1 || is_unaligned_window(window) {
            k += 1;
            continue;
        }
        let overlap_start = target_start0.max(window_target_start);
        let overlap_end = target_end0.min(window_target_end);
        let overlap = if overlap_end > overlap_start {
            overlap_end - overlap_start
        } else {
            0
        };
        if overlap == 0 {
            k += 1;
            continue;
        }
        if k < first {
            first = k;
        }
        if k > last {
            last = k;
        }
        aligned_total += overlap;
        k += 1;
    }
    if aligned_total == 0 || first == i32::MAX as i64 || last < 0 {
        return false;
    }
    if first == last && is_unaligned_window(&windows[first as usize]) {
        return false;
    }

    let mut query_bounds = [0i64; 2];
    let mut target_bounds = [0i64; 2];
    win_boundary_offset(
        windows,
        first as usize,
        target_start0,
        query_length,
        &mut query_bounds[0],
        &mut target_bounds[0],
    );
    win_boundary_offset(
        windows,
        last as usize,
        target_end0 - 1,
        query_length,
        &mut query_bounds[1],
        &mut target_bounds[1],
    );
    query_bounds[1] += 1;
    target_bounds[1] += 1;
    if (aligned_total as f64) < ((target_bounds[1] - target_bounds[0]) as f64) * overlap_rate {
        return false;
    }

    gen_raln(
        targets,
        query,
        region,
        extz,
        first,
        last + 1,
        target_start0,
        target_end0,
        query_length,
        region.y_strand as u8,
        error_rate,
        trace,
        result,
    );
    assert!(
        i64::from(result.ts) >= target_start0
            && i64::from(result.te) <= target_end0
            && i64::from(result.qs) >= 0
            && i64::from(result.qe) <= query_length
    );
    let similarity_threshold = error_rate + R_SIMI_W.min(error_rate / 2.0);
    (result.sec as f64) <= ((result.qe - result.qs) as f64) * similarity_threshold
        && (result.sec as f64) <= ((result.te - result.ts) as f64) * similarity_threshold
}

/// `gen_sub_ov_adv_cigar` (`Correct.cpp` 24695): maps every read (or
/// contained read) of the unitig onto the overlap and re-aligns it through
/// [`extract_subov_cigar`], collecting the sub-overlaps and their trace
/// indices.
///
/// The reference's `sid` parameter feeds commented-out debug prints only.
#[allow(clippy::too_many_arguments)]
pub(crate) fn gen_sub_ov_adv_cigar(
    targets: &TargetSource,
    unitig: &MaUtg,
    region: &OverlapRegion,
    query: &[u8],
    extz: &mut BitExtz,
    query_length: i64,
    overlap_rate: f64,
    error_rate: f64,
    contained: bool,
    contain: &crate::assembly::ul_index::UlContain,
    read_lengths: &[u64],
    overlap_id: i64,
    trace: &mut Vec<Rtrace>,
    result: &mut Vec<UlOv>,
) -> usize {
    let before = result.len();
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

    let mut entry_trace = Rtrace::default();
    let mut entry = UlOv::default();
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
            if extract_subov_cigar(
                targets,
                query,
                region,
                extz,
                clip_start,
                clip_end,
                overlap_rate,
                Some(&mut window_hint),
                query_length,
                error_rate,
                &mut entry_trace,
                &mut entry,
            ) {
                entry_trace.oid = overlap_id as i32;
                entry.ts -= clip_start as u32;
                entry.te -= clip_start as u32;
                entry.rev = u32::from(region.y_strand != 0) != ((packed >> 32) & 1) as u32;
                if entry.rev {
                    let read_length = read_lengths[read_identifier] as i64;
                    let flipped_start = i64::from(entry.ts);
                    let flipped_end = i64::from(entry.te);
                    entry.ts = (read_length - flipped_end) as u32;
                    entry.te = (read_length - flipped_start) as u32;
                }
                entry.el = true;
                entry.tn = read_identifier as u32;
                entry.qn = trace.len() as u32;
                result.push(entry);
                trace.push(entry_trace);
            }
        }
    } else {
        let count = contain.idx[region.y_identifier as usize] as u32 as usize;
        let begin = (contain.idx[region.y_identifier as usize] >> 32) as usize;
        for index in 0..count {
            let contained_entry = contain.rids[begin + index];
            let read_identifier = (contained_entry.x >> 1) as usize;
            let start = i64::from(contained_entry.s);
            let end = i64::from(contained_entry.e);
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
            if extract_subov_cigar(
                targets,
                query,
                region,
                extz,
                clip_start,
                clip_end,
                overlap_rate,
                Some(&mut window_hint),
                query_length,
                error_rate,
                &mut entry_trace,
                &mut entry,
            ) {
                entry_trace.oid = overlap_id as i32;
                entry.ts -= clip_start as u32;
                entry.te -= clip_start as u32;
                entry.rev = u32::from(region.y_strand != 0) != contained_entry.x & 1;
                if entry.rev {
                    let read_length = read_lengths[read_identifier] as i64;
                    let flipped_start = i64::from(entry.ts);
                    let flipped_end = i64::from(entry.te);
                    entry.ts = (read_length - flipped_end) as u32;
                    entry.te = (read_length - flipped_start) as u32;
                }
                entry.el = false;
                entry.tn = read_identifier as u32;
                entry.qn = trace.len() as u32;
                result.push(entry);
                trace.push(entry_trace);
            }
        }
    }
    result.len() - before
}

/// `is_alnw` (`Hash_Table.h` 221): the hit was refined away.
pub(crate) fn is_alnw(hit: &KmerHit) -> bool {
    hit.read_identifier == ALNW_SENTINEL
}

/// `is_pri_aln` (`Hash_Table.h` 222): the hit terminates a primary-alignment
/// block.
pub(crate) fn is_pri_aln(hit: &KmerHit) -> bool {
    hit.read_identifier == ALNW_SENTINEL || hit.count >= hit.read_identifier
}

/// `return_t_chain` (`Correct.cpp` 23050): refines the chain bracket of one
/// overlap in place, marking the hits past the refined chain with the `alnw`
/// sentinel. Returns the refined chain length.
pub(crate) fn return_t_chain(region: &mut OverlapRegion, candidates: &mut CandidatesList) -> i64 {
    let start = region.shared_seed as usize;
    let mut end = start;
    let first_identifier = candidates.list[start].read_identifier;
    while end < candidates.length as usize
        && candidates.list[end].read_identifier == first_identifier
        && candidates.list[end].read_identifier != ALNW_SENTINEL
    {
        end += 1;
    }
    let span = (end - start) as i64;
    let refined = chain_refine(
        &mut candidates.list[start..end],
        span,
        &mut candidates.chain_data,
        50,
        5000,
        512,
        16,
    );
    for hit in &mut candidates.list[start + refined as usize..end] {
        hit.read_identifier = ALNW_SENTINEL;
    }
    refined
}

/// `adjust_base_coordinates` (`Correct.cpp` 16547): snaps one candidate
/// interval onto the chain bracket that covers it, leaving `ts`/`te` open
/// (`(uint32_t)-1`) so the block pass treats it as semi-global.
pub(crate) fn adjust_base_coordinates(
    hits: &[KmerHit],
    hit_count: i64,
    entry: &mut UlOv,
    query_length: i64,
    mut cursor: i64,
) -> i64 {
    entry.ts = u32::MAX;
    entry.te = u32::MAX;
    entry.sec = 3;
    entry.qn = u32::MAX;
    entry.tn = u32::MAX;
    if hit_count == 0 {
        return cursor;
    }
    let start = i64::from(entry.qs);
    let end = i64::from(entry.qe);

    let mut index = cursor;
    while index >= 0 && i64::from(hits[index as usize].self_offset) >= start {
        index -= 1;
    }
    if index < 0 {
        index = 0;
    }
    while index < hit_count && i64::from(hits[index as usize].self_offset) < start {
        index += 1;
    }
    if index >= 0 && (index >= hit_count || i64::from(hits[index as usize].self_offset) > start) {
        index -= 1;
    }
    let begin = index;
    let mut index = if index >= 0 { index } else { 0 };
    while index < hit_count && i64::from(hits[index as usize].self_offset) < end {
        index += 1;
    }
    let end_index = index;
    cursor = index;

    if begin >= 0 {
        entry.qs = hits[begin as usize].self_offset;
        entry.ts = hits[begin as usize].offset;
        entry.qn = begin as u32;
    } else {
        entry.qs = 0;
        entry.qn = u32::MAX;
    }

    if end_index < hit_count {
        entry.qe = hits[end_index as usize].self_offset;
        entry.te = hits[end_index as usize].offset;
        entry.tn = end_index as u32;
    } else {
        entry.qe = query_length as u32;
        entry.tn = hit_count as u32;
    }

    cursor
}

/// `fusion_chain_ovlp` (`Correct.cpp` 16616): snaps every candidate interval
/// onto the chain bracket and merges the overlapping ones.
#[allow(clippy::too_many_arguments)]
pub(crate) fn fusion_chain_ovlp(
    hits: &[KmerHit],
    hit_count: i64,
    entries: &mut [UlOv],
    mut count: i64,
    query_length: i64,
) -> i64 {
    let mut sort_ok = true;
    let mut cursor = 0i64;
    for i in 0..count as usize {
        entries[i].sec = 6;
        cursor = adjust_base_coordinates(hits, hit_count, &mut entries[i], query_length, cursor);
        if i > 0 && entries[i].qs < entries[i - 1].qe {
            sort_ok = false;
        }
    }
    if count <= 1 {
        return count;
    }

    if !sort_ok {
        // `radix_sort_uov_srt_qs` is a stable 32-bit radix sort on `qs`.
        crate::container::ksort::radix_sort_by_key_32(entries, &|entry: &UlOv| entry.qs);
    }
    let mut write = 1usize;
    for read in 1..count as usize {
        let previous = entries[write - 1];
        let overlap_start = previous.qs.max(entries[read].qs);
        let overlap_end = previous.qe.min(entries[read].qe);
        let overlap = i64::from(overlap_end) - i64::from(overlap_start);
        if overlap >= 0 {
            let source = entries[read];
            let previous = &mut entries[write - 1];
            previous.qe = previous.qe.max(source.qe);
            previous.te = previous.te.max(source.te);
            previous.tn = previous.tn.max(source.tn);
        } else {
            entries[write] = entries[read];
            write += 1;
        }
    }
    count = write as i64;
    count
}

/// `ovlp_base_aln_all` (`Correct.cpp` 16652): the postprocess fallback of
/// [`ovlp_base_aln`], aligning every hit pair of the failed block.
///
/// The reference's `rid` parameter feeds commented-out debug prints only.
#[allow(clippy::too_many_arguments)]
pub(crate) fn ovlp_base_aln_all(
    region: &OverlapRegion,
    hits: &[KmerHit],
    hit_count: i64,
    block_start: i64,
    block_end: i64,
    targets: &TargetSource,
    query: &[u8],
    bracket: &super::ug_lalign::ChainBracket,
    query_total_length: i64,
    target_total_length: i64,
    window_length: i64,
    extz: &mut BitExtz,
    output: &mut WindowCigar,
    error_rate: f64,
) {
    for i in block_start + 1..=block_end {
        let previous = i - 1;
        let mut query_start;
        let mut query_end;
        let mut target_start;
        let mut target_end;
        let mode;

        if previous >= 0 {
            query_start = i64::from(hits[previous as usize].self_offset);
            target_start = i64::from(hits[previous as usize].offset);
        } else {
            query_start = i64::from(bracket.qs);
            target_start = -1;
        }
        if i < hit_count {
            query_end = i64::from(hits[i as usize].self_offset);
            target_end = i64::from(hits[i as usize].offset);
        } else {
            query_end = i64::from(bracket.qe);
            target_end = -1;
        }

        if target_start != -1 && target_end != -1 {
            mode = MODE_GLOBAL;
        } else if target_start != -1 {
            mode = MODE_FORWARD_EXTENSION;
        } else if target_end != -1 {
            mode = MODE_BACKWARD_EXTENSION;
        } else {
            mode = 3;
        }
        assert!(mode != 3);
        if mode == MODE_FORWARD_EXTENSION || mode == MODE_BACKWARD_EXTENSION {
            crate::correction::cal_exz::adjust_ext_offset(
                &mut query_start,
                &mut query_end,
                &mut target_start,
                &mut target_end,
                query_total_length,
                target_total_length,
                0,
                mode,
            );
        }
        // `at cns chain, the base alignment fails; there is no anchor between
        // soff and eoff`
        let aligned = if block_end - block_start <= 1
            && ((query_end - query_start) >> 1) < MAX_CNS_ERROR
        {
            false
        } else {
            crate::correction::aln_adv::hc_aln_exz_adv(
                region,
                targets,
                extz,
                query,
                query_start,
                query_end,
                target_start,
                target_end,
                mode,
                window_length,
                query_total_length,
                error_rate,
                crate::correction::cal_exz::MAXIMUM_SINGLE_LENGTH,
                crate::correction::cal_exz::MAXIMUM_SINGLE_ERROR,
                crate::correction::cal_exz::FORCED_SINGLE_LENGTH,
                -1,
                output,
            )
        };
        if !aligned {
            crate::correction::cal_exz::push_unmapped_window(
                output,
                query_start as i32,
                query_end - 1,
                target_start,
                target_end - 1,
                mode,
            );
        }
    }
}

/// `ovlp_base_aln` (`Correct.cpp` 16704): the block pass of `rechain_aln`.
/// Blocks break at primary-alignment hits; two facing aligned hits whose
/// strands flank the target span are replaced by the stored windows instead of
/// a fresh alignment, and failed blocks fall back to [`ovlp_base_aln_all`].
///
/// The reference's `rid` parameter feeds commented-out debug prints only.
#[allow(clippy::too_many_arguments)]
pub(crate) fn ovlp_base_aln(
    region: &OverlapRegion,
    hits: &[KmerHit],
    hit_count: i64,
    bracket: &super::ug_lalign::ChainBracket,
    window_length: i64,
    targets: &TargetSource,
    query: &[u8],
    extz: &mut BitExtz,
    output: &mut WindowCigar,
    error_rate: f64,
    query_total_length: i64,
    target_total_length: i64,
) {
    let first_block: i64 = if bracket.qn == u32::MAX {
        -1
    } else {
        i64::from(bracket.qn)
    };
    let last_block = i64::from(bracket.tn);
    assert!(last_block >= first_block + 1);

    let mut previous = first_block;
    let mut index = first_block + 1;
    while index <= last_block {
        if index == last_block || is_pri_aln(&hits[index as usize]) {
            let mut query_start;
            let mut query_end;
            let mut target_start;
            let mut target_end;
            let mode;
            let mut done;

            if previous >= 0 {
                query_start = i64::from(hits[previous as usize].self_offset);
                target_start = i64::from(hits[previous as usize].offset);
            } else {
                query_start = i64::from(bracket.qs);
                target_start = -1;
            }
            if index < hit_count {
                query_end = i64::from(hits[index as usize].self_offset);
                target_end = i64::from(hits[index as usize].offset);
            } else {
                query_end = i64::from(bracket.qe);
                target_end = -1;
            }

            if target_start != -1 && target_end != -1 {
                mode = MODE_GLOBAL;
            } else if target_start != -1 {
                mode = MODE_FORWARD_EXTENSION;
            } else if target_end != -1 {
                mode = MODE_BACKWARD_EXTENSION;
            } else {
                mode = 3;
            }

            if mode == MODE_GLOBAL
                && is_alnw(&hits[previous as usize])
                && is_alnw(&hits[index as usize])
                && hits[previous as usize].strand_value() == 0
                && hits[index as usize].strand_value() == 1
            {
                done = true;
                crate::correction::cal_exz::push_replace_window_adv(
                    region,
                    window_length,
                    output,
                    query_start as i32,
                    query_end - 1,
                    target_start,
                    target_end - 1,
                    mode,
                );
            } else if mode != 3 {
                if mode == MODE_FORWARD_EXTENSION || mode == MODE_BACKWARD_EXTENSION {
                    crate::correction::cal_exz::adjust_ext_offset(
                        &mut query_start,
                        &mut query_end,
                        &mut target_start,
                        &mut target_end,
                        query_total_length,
                        target_total_length,
                        0,
                        mode,
                    );
                }
                done = crate::correction::aln_adv::hc_aln_exz_adv(
                    region,
                    targets,
                    extz,
                    query,
                    query_start,
                    query_end,
                    target_start,
                    target_end,
                    mode,
                    window_length,
                    query_total_length,
                    error_rate,
                    MAX_CNS_LENGTH,
                    MAX_CNS_ERROR,
                    FORCED_CNS_LENGTH,
                    -1,
                    output,
                );
            } else {
                done = false;
            }

            if !done {
                ovlp_base_aln_all(
                    region,
                    hits,
                    hit_count,
                    previous,
                    index,
                    targets,
                    query,
                    bracket,
                    query_total_length,
                    target_total_length,
                    window_length,
                    extz,
                    output,
                    error_rate,
                );
            }
            previous = index;
        }
        index += 1;
    }
}

/// `cigar_gen_by_chain_adv_local` (`Correct.cpp` 17598): re-aligns one
/// overlap from its chain brackets. The brackets come from `entries` (already
/// snapped by [`fusion_chain_ovlp`]); unaligned windows are re-chained once,
/// and the assembled windows replace the overlap's own list.
///
/// The reference's `rid` parameter feeds commented-out debug prints only.
#[allow(clippy::too_many_arguments)]
pub(crate) fn cigar_gen_by_chain_adv_local(
    region: &mut OverlapRegion,
    candidates: &mut CandidatesList,
    entries: &mut [UlOv],
    mut count: i64,
    window_length: u64,
    targets: &TargetSource,
    query: &[u8],
    extz: &mut BitExtz,
    auxiliary: &mut OverlapRegion,
    error_rate: f64,
    query_total_length: i64,
    kmer_hit: i64,
) {
    if count <= 0 {
        return;
    }
    let chain_index = region.shared_seed as usize;
    let target_total_length = targets.length(region.y_identifier as usize);
    let mut chain_end = chain_index;
    let first_identifier = candidates.list[chain_index].read_identifier;
    while chain_end < candidates.length as usize
        && candidates.list[chain_end].read_identifier == first_identifier
    {
        chain_end += 1;
    }
    let hit_count = (chain_end - chain_index) as i64;
    let hits = candidates.list[chain_index..chain_end].to_vec();

    count = fusion_chain_ovlp(
        &hits,
        hit_count,
        entries,
        count,
        query_total_length,
    );

    auxiliary.window_cigar.windows.clear();
    auxiliary.window_cigar.cells.clear();
    auxiliary.y_identifier = region.y_identifier;
    auxiliary.y_strand = region.y_strand;
    auxiliary.x_position_start = region.x_position_start;
    auxiliary.x_position_end = region.x_position_end;
    auxiliary.y_position_start = region.y_position_start;
    auxiliary.y_position_end = region.y_position_end;

    for entry in entries.iter().take(count as usize) {
        let bracket = super::ug_lalign::ChainBracket {
            qn: entry.qn,
            tn: entry.tn,
            qs: entry.qs,
            qe: entry.qe,
            ts: entry.ts,
            te: entry.te,
        };
        super::ug_lalign::ovlp_base_direct(
            region,
            &hits,
            hit_count,
            &bracket,
            window_length as i64,
            targets,
            query,
            extz,
            &mut auxiliary.window_cigar,
            error_rate,
            query_total_length,
            target_total_length,
        );
    }

    let initial_window_count = auxiliary.window_cigar.windows.len();
    for i in 0..initial_window_count {
        if !is_unaligned_window(&auxiliary.window_cigar.windows[i]) {
            continue;
        }
        // will overwrite `hits`; does not matter
        super::ug_lalign::rechain_alignment_ul(
            region,
            candidates,
            &mut auxiliary.window_cigar,
            i,
            window_length as i64,
            targets,
            query,
            extz,
            error_rate,
            query_total_length,
            target_total_length,
            kmer_hit,
        );
    }
    if auxiliary.window_cigar.windows.len() > initial_window_count {
        let mut write = 0usize;
        for read in 0..auxiliary.window_cigar.windows.len() {
            if read < initial_window_count
                && is_unaligned_window(&auxiliary.window_cigar.windows[read])
            {
                continue;
            }
            auxiliary.window_cigar.windows[write] = auxiliary.window_cigar.windows[read];
            write += 1;
        }
        auxiliary.window_cigar.windows.truncate(write);
        // `radix_sort_window_list_xs_srt` is a stable 32-bit radix sort on
        // the signed `x_start`.
        crate::container::ksort::radix_sort_by_key_32(
            &mut auxiliary.window_cigar.windows,
            &|window: &WindowList| window.x_start as u32,
        );
    }

    // `update_overlap_region` (`Correct.cpp` 17249)
    region.update_from_windows(
        &auxiliary.window_cigar,
        query_total_length,
        target_total_length,
    );
    debug_assert!(region.x_position_end >= region.x_position_start);
    debug_assert!(region.y_position_end >= region.y_position_start);
}

/// `dedup_ul_ov_t` (`Correct.cpp` 25050): drops the contained-read sub-overlaps
/// that nearly duplicate an earlier one of the same target, keeping the one
/// with the smaller error.
pub(crate) fn dedup_ul_ov_t(entries: &mut Vec<UlOv>) {
    let count = entries.len();
    const RATIO: f64 = 0.95;
    for entry in entries.iter_mut() {
        entry.tn = entry.tn.wrapping_shl(1) | u32::from(entry.rev);
    }
    // `radix_sort_ul_ov_srt_tn1` is a stable 32-bit radix sort on `tn`.
    crate::container::ksort::radix_sort_by_key_32(entries, &|entry: &UlOv| entry.tn);

    let mut group_start = 0usize;
    let mut write = 0usize;
    let mut k = 1usize;
    while k <= count {
        if k == count || entries[k].tn != entries[group_start].tn {
            for z in group_start..k {
                let mut r = write as i64 - 1;
                while r >= 0
                    && entries[r as usize].tn == (entries[z].tn >> 1)
                    && entries[r as usize].rev == entries[z].rev
                {
                    let candidate = entries[z].clone();
                    let kept = entries[r as usize];
                    let query_overlap = (i64::from(candidate.qe.min(kept.qe))
                        - i64::from(candidate.qs.max(kept.qs)))
                    .max(0);
                    let target_overlap = (i64::from(candidate.te.min(kept.te))
                        - i64::from(candidate.ts.max(kept.ts)))
                    .max(0);
                    if query_overlap as f64 >= (i64::from(kept.qe - kept.qs) as f64) * RATIO
                        && query_overlap as f64
                            >= (i64::from(candidate.qe - candidate.qs) as f64) * RATIO
                        && target_overlap as f64 >= (i64::from(kept.te - kept.ts) as f64) * RATIO
                        && target_overlap as f64
                            >= (i64::from(candidate.te - candidate.ts) as f64) * RATIO
                    {
                        break;
                    }
                    r -= 1;
                }
                if r >= 0
                    && entries[r as usize].tn == (entries[z].tn >> 1)
                    && entries[r as usize].rev == entries[z].rev
                {
                    if entries[z].sec < entries[r as usize].sec {
                        entries[r as usize] = entries[z];
                        entries[r as usize].tn >>= 1;
                    }
                    continue;
                }
                entries[write] = entries[z];
                entries[write].tn >>= 1;
                write += 1;
            }
            group_start = k;
        }
        k += 1;
    }
    entries.truncate(write);
}

/// The `ul_local_aln` second-pass outputs: the query-side sub-overlaps
/// (`aln`), the contained-read sub-overlaps (`cln`) and their trace entries.
#[derive(Default)]
pub struct UlLocalAlignment {
    pub aligned: Vec<UlOv>,
    pub contained: Vec<UlOv>,
    pub trace: Vec<Rtrace>,
}

/// `ul_local_aln` (`Correct.cpp` 24826): refines one overlap. Without
/// `auxiliary` this is the first pass — align and filter. With `auxiliary`
/// the overlap is split into per-read sub-overlaps, re-aligned through the
/// chain brackets, and the per-read intervals are collected.
///
/// The reference's `rid`/`oid` parameters feed the trace bookkeeping and
/// commented-out debug prints.
#[allow(clippy::too_many_arguments)]
pub(crate) fn ul_local_aln(
    region: &mut OverlapRegion,
    candidates: &mut CandidatesList,
    extz: &mut BitExtz,
    targets: &TargetSource,
    unitigs: &[MaUtg],
    contain: &crate::assembly::ul_index::UlContain,
    read_lengths: &[u64],
    query: &[u8],
    error_rate: f64,
    window_length: u64,
    minimum_overlap: u64,
    overlap_rate: f64,
    output: &mut UlLocalAlignment,
    query_total_length: i64,
    overlap_id: usize,
    kmer_hit: i64,
    auxiliary: Option<&mut OverlapRegion>,
) -> bool {
    let unitig = &unitigs[region.y_identifier as usize];
    let Some(auxiliary) = auxiliary else {
        align_ul_ed_post_extz(
            region,
            extz,
            targets,
            query,
            &HcAlignParameters {
                error_rate,
                window_length: window_length as i64,
                overlap_cutoff: -1.0,
                force_alignment: true,
            },
        );
        let overlap_length =
            i64::from(region.x_position_end) + 1 - i64::from(region.x_position_start);
        let aligned_length = i64::from(region.alignment_length);
        if aligned_length as f64 >= overlap_length as f64 * overlap_rate {
            return true;
        }
        if (aligned_length as u64) <= minimum_overlap {
            return false;
        }
        return win_cluster_filter(&region.window_cigar, minimum_overlap, overlap_rate / 2.0);
    };

    let overlap_length = i64::from(region.x_position_end) + 1 - i64::from(region.x_position_start);
    let aligned_length = i64::from(region.alignment_length);
    if (aligned_length as u64) <= minimum_overlap {
        return false;
    }
    let output_start = output.aligned.len();

    if (overlap_length as f64) * overlap_rate <= aligned_length as f64 {
        // `[ts, te)` -> whole interval of the unitig adjusted by the reverse
        // strand
        let entry = UlOv {
            ts: region.y_position_start,
            te: region.y_position_end + 1,
            qs: region.x_position_start,
            qe: region.x_position_end + 1,
            rev: region.y_strand != 0,
            el: false,
            tn: u32::MAX,
            qn: u32::MAX,
            ..UlOv::default()
        };
        output.aligned.push(entry);
    } else {
        gen_sub_ov_adv(
            region,
            overlap_rate,
            unitig,
            contain,
            false,
            read_lengths,
            &mut output.aligned,
        );
        gen_sub_ov_adv(
            region,
            overlap_rate,
            unitig,
            contain,
            true,
            read_lengths,
            &mut output.aligned,
        );
    }

    if output.aligned.len() <= output_start {
        return false;
    }
    // `radix_sort_ul_ov_srt_qs1` is a stable 32-bit radix sort on `qs`.
    crate::container::ksort::radix_sort_by_key_32(
        &mut output.aligned[output_start..],
        &|entry: &UlOv| entry.qs,
    );
    // merge the overlapping candidates
    let mut write = output_start;
    let mut merged: Option<usize> = None;
    for read in output_start..output.aligned.len() {
        match merged {
            Some(previous) => {
                if i64::from(output.aligned[read].qs) < i64::from(output.aligned[previous].qe)
                    && i64::from(output.aligned[read].qe) > i64::from(output.aligned[previous].qe)
                {
                    output.aligned[previous].qe = output.aligned[read].qe;
                }
            }
            None => {
                output.aligned[write] = output.aligned[read];
                merged = Some(write);
                write += 1;
            }
        }
    }
    output.aligned.truncate(write);

    return_t_chain(region, candidates);
    let entry_count = (output.aligned.len() - output_start) as i64;
    cigar_gen_by_chain_adv_local(
        region,
        candidates,
        &mut output.aligned[output_start..],
        entry_count,
        window_length,
        targets,
        query,
        extz,
        auxiliary,
        error_rate,
        query_total_length,
        kmer_hit,
    );

    output.aligned.truncate(output_start);
    gen_sub_ov_adv_cigar(
        targets,
        unitig,
        region,
        query,
        extz,
        query_total_length,
        overlap_rate,
        error_rate,
        false,
        contain,
        read_lengths,
        overlap_id as i64,
        &mut output.trace,
        &mut output.aligned,
    );
    gen_sub_ov_adv_cigar(
        targets,
        unitig,
        region,
        query,
        extz,
        query_total_length,
        overlap_rate,
        error_rate,
        true,
        contain,
        read_lengths,
        overlap_id as i64,
        &mut output.trace,
        &mut output.contained,
    );
    true
}

/// `ul_rid_lalign_adv` (`Correct.cpp` 25083): the per-read driver. The first
/// call (`auxiliary == None`) filters the chain overlaps in place and returns
/// the compacted length; the second call splits every surviving overlap into
/// per-read sub-overlaps accumulated into `output`.
///
/// `minimum_overlap` is the reference's literal `1000`; `overlap_rate` is
/// `OVERLAP_THRESHOLD_NOSI_FILTER` (0.7). The reference's `sid` parameter
/// feeds commented-out debug prints only.
#[allow(clippy::too_many_arguments)]
pub(crate) fn ul_rid_lalign_adv(
    overlaps: &mut [OverlapRegion],
    length: usize,
    candidates: &mut CandidatesList,
    extz: &mut BitExtz,
    targets: &TargetSource,
    unitigs: &[MaUtg],
    contain: &crate::assembly::ul_index::UlContain,
    read_lengths: &[u64],
    query: &[u8],
    query_length: i64,
    error_rate: f64,
    window_length: i64,
    output: &mut UlLocalAlignment,
    kmer_length: i64,
    auxiliary: Option<&mut OverlapRegion>,
) -> usize {
    match auxiliary {
        None => {
            let mut kept = 0usize;
            for index in 0..length {
                overlaps[index].shared_seed = overlaps[index].non_homopolymer_errors as i32;
                if !ul_local_aln(
                    &mut overlaps[index],
                    candidates,
                    extz,
                    targets,
                    unitigs,
                    contain,
                    read_lengths,
                    query,
                    error_rate,
                    window_length as u64,
                    1000,
                    0.7,
                    output,
                    query_length,
                    index,
                    kmer_length,
                    None,
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
        Some(auxiliary) => {
            output.contained.clear();
            output.trace.clear();
            for index in 0..length {
                overlaps[index].shared_seed = overlaps[index].non_homopolymer_errors as i32;
                ul_local_aln(
                    &mut overlaps[index],
                    candidates,
                    extz,
                    targets,
                    unitigs,
                    contain,
                    read_lengths,
                    query,
                    error_rate,
                    window_length as u64,
                    1000,
                    0.7,
                    output,
                    query_length,
                    index,
                    kmer_length,
                    Some(auxiliary),
                );
            }
            if !output.contained.is_empty() {
                dedup_ul_ov_t(&mut output.contained);
            }
            if !output.contained.is_empty() {
                output
                    .aligned
                    .extend(output.contained.iter().copied());
            }
            length
        }
    }
}

/// `gen_aux_ovlp` (`inter.cpp` 7561): the scratch overlap the second pass
/// writes into, cloned from the first overlap of the list.
pub(crate) fn gen_aux_ovlp(overlaps: &[OverlapRegion], length: usize) -> OverlapRegion {
    if length > 0 {
        overlaps[0].clone()
    } else {
        OverlapRegion::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::asg::MaUtg;
    use crate::correction::window::{WindowCigar, WindowList};
    use crate::correction::window_chain::is_unaligned_window;
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

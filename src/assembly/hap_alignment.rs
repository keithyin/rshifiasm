//! The haplotype-alignment candidate scan of `purge_dups`, the counterpart of
//! `Purge_Dups.cpp` 611–3216 (the live `*_advance*` family only).
//!
//! Every unitig votes for the unitigs its reads overlap (the vote half stays
//! in [`crate::assembly::purge_dups`]); for each surviving candidate the scan
//! collects read-to-read arcs across the two unitigs, chains them into
//! haplotype candidates (`hap_chaining`), validates each chain against the
//! read graph (`calculate_pair_hap_similarity_advance`) and pushes the
//! surviving ones into the hap-overlap lists consumed by the purge graph.
//!
//! Float-precision notes, all verified against the C source:
//!
//! * `prefilter` compares in `float`: `min_count * Hap_rate` promotes both
//!   sides to `f32`.
//! * `determin_hap_alignment_boundary_single_side` compares
//!   `max_count > min_count * H_rate` in `f32` as well.
//! * `hap_chaining` keeps `band_width_penalty = 1 / band_width_threshold` in
//!   `f32`, and `(long long)(score * gap_rate * band_width_penalty)` converts
//!   the `i64` score to `f32` before each multiply.
//! * `calculate_pair_hap_similarity_advance` mixes `f64` totals with an
//!   `f32` rate; the product is computed in `f64`.
//! * `push_hap_can` never writes `index_beg` and `score`; the C comparator
//!   `cmp_hap_candidates` therefore ranks ties on uninitialized memory. The
//!   port keeps both at their default (0) and documents the divergence; the
//!   end-to-end gate catches any observable effect.

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::hap_cov::{ArcOffset, HapCov};
use crate::assembly::hit::{ma_hit_to_arc, CoverageCut, ReadToUnitig};
use crate::assembly::trans_cov::{compare_arc_offsets, get_x_offset, get_y_offset, quick_lis};
use crate::assembly::unitig::UnitigGraph;
use crate::correction::ec_driver::PafList;

/// `PLOID`.
pub const PLOID: u32 = 0;
/// `NON_PLOID`.
pub const NON_PLOID: u32 = 1;
/// `CUTOFF_THRES` of `get_pair_hap_similarity`.
const CUTOFF_THRES: u32 = 1000;
/// `CHAIN_UNMATCH` (`Process_Read.h`), the non-chain penalty of
/// `get_chain_score`; `CHAIN_MATCH` is 1.
const CHAIN_UNMATCH: f64 = 0.334;
/// `DUP_OVLP_RATE` of `hap_chaining`.
const DUP_OVLP_RATE: f64 = 0.75;
/// `CHAIN_FILTER_RATE` (`Purge_Dups.h`).
pub const CHAIN_FILTER_RATE: f64 = 0.7;
/// `TARGET_SGIFT` (sic) of `modify_target_interval`.
const TARGET_SHIFT: i64 = 3;
/// `REV_W`: the weight `init_pip_hh` gives a read with no reverse coverage.
pub const REV_W: i64 = 8;

/// `hap_candidates`: one chained haplotype-alignment candidate.
#[derive(Debug, Clone)]
pub struct HapCandidate {
    /// `weight`.
    pub weight: u64,
    /// `x_beg_pos` (read indices inside `x`).
    pub x_beg_pos: u32,
    /// `x_end_pos`.
    pub x_end_pos: u32,
    /// `y_beg_pos`.
    pub y_beg_pos: u32,
    /// `y_end_pos`.
    pub y_end_pos: u32,
    /// `index_beg` — the reference leaves this uninitialized in
    /// `push_hap_can`; see the module notes.
    pub index_beg: u32,
    /// `index_end`, also the overlap `type` (`Get_type`).
    pub index_end: u32,
    /// `score`.
    pub score: i64,
    /// `rev`, copied from the representative arc's `el`.
    pub rev: u8,
    /// `t`, the representative arc (the median-weight chain element).
    pub arc: AsgArc,
}

impl Default for HapCandidate {
    fn default() -> Self {
        Self {
            weight: 0,
            x_beg_pos: 0,
            x_end_pos: 0,
            y_beg_pos: 0,
            y_end_pos: 0,
            index_beg: 0,
            index_end: 0,
            score: 0,
            rev: 0,
            arc: AsgArc::default(),
        }
    }
}

/// Shared read-graph context of the scan; one struct keeps the C signatures
/// manageable without changing any semantics.
pub struct HapAlignmentContext<'a> {
    /// `read_g`.
    pub read_graph: &'a Asg,
    /// `reverse_sources`.
    pub reverse_sources: &'a [PafList],
    /// `ruIndex`.
    pub read_to_unitig: &'a ReadToUnitig,
    /// `coverage_cut`.
    pub coverage_cut: &'a [CoverageCut],
    /// `position_index`.
    pub position_index: &'a [u64],
    /// `max_hang`.
    pub maximum_hang: i32,
    /// `asm_opt.max_hang_rate`.
    pub maximum_hang_rate: f32,
    /// `min_ovlp`.
    pub minimum_overlap: i32,
    /// `asm_opt.purge_level_primary`.
    pub purge_level: u32,
}

/// `ma_hit2arc` + the shared target-side filters of the scan loops.
///
/// Returns `None` in exactly the cases the C loops `continue` on: a contained
/// overlap, a deleted target cut/read, or a target whose unitig is not
/// `y_uid`. On success the arc's `el` is *not* yet the strand byte — the
/// callers overwrite it through `get_xy_pos`.
fn scan_arc(
    ctx: &HapAlignmentContext,
    hit: &crate::overlap::paf::PafHit,
    y_uid: u32,
) -> Option<AsgArc> {
    let query_cut = &ctx.coverage_cut[hit.query_identifier as usize];
    let target_cut = &ctx.coverage_cut[hit.target_identifier as usize];
    if target_cut.deleted || ctx.read_graph.sequences[hit.target_identifier as usize].deleted {
        return None;
    }

    let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
    let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
    let (result, arc) = ma_hit_to_arc(
        hit,
        query_length,
        target_length,
        ctx.maximum_hang,
        ctx.maximum_hang_rate,
        ctx.minimum_overlap,
    );
    // A contained overlap cannot extend either unitig.
    if result < 0 {
        return None;
    }

    let target_read = arc.v >> 1;
    if ctx.read_graph.sequences[target_read as usize].deleted {
        return None;
    }
    let (hap_uid, is_unitig) = ctx.read_to_unitig.get(target_read);
    if is_unitig == 0 || hap_uid == u32::MAX || hap_uid != y_uid {
        return None;
    }
    Some(arc)
}

/// `prefilter` (`Purge_Dups.cpp:611`).
///
/// `dir == 0` aligns the two unitigs head-to-head, `dir == 1` head-to-tail;
/// either direction passing is enough for the caller.
pub fn prefilter(
    x_pos: u32,
    y_pos: u32,
    x_len: u32,
    y_len: u32,
    dir: u32,
    hap_rate: f32,
    seed_occ: u32,
) -> u32 {
    let x_left_len = x_pos;
    let x_right_beg = x_pos;
    let x_right_len = x_len - x_right_beg;
    let (y_left_beg, y_left_len, y_right_beg, y_right_len) = if dir == 0 {
        (0u32, y_pos, y_pos, y_len - y_pos)
    } else {
        (y_pos + 1, y_len - (y_pos + 1), 0, y_pos + 1)
    };

    let max_count = seed_occ;
    let min_count = x_left_len.min(y_left_len) + x_right_len.min(y_right_len);
    if min_count == 0 {
        return NON_PLOID;
    }
    // `float` comparison in the reference.
    if (max_count as f32) <= (min_count as f32) * hap_rate {
        return NON_PLOID;
    }
    PLOID
}

/// `get_xy_pos` (`Purge_Dups.cpp:641`), the `position_index` variant.
///
/// Returns `(packed, rev)`; the reference writes `rev` through the arc's
/// `el` pointer slot, so callers store it into `arc.exact_region`.
pub fn get_xy_pos(
    read_graph: &Asg,
    arc: &AsgArc,
    v_in_unitig: u32,
    w_in_unitig: u32,
    x_unitig_len: u32,
    y_unitig_len: u32,
    position_index: &[u64],
) -> (u64, u8) {
    let mut x_pos = u32::MAX;
    let mut y_pos = u32::MAX;
    let mut x_dir = 0u8;
    let mut y_dir = 0u8;

    if (arc.ul >> 32) as u32 == v_in_unitig {
        // The arc leaves the unitig read at its end.
        x_pos = ((position_index[(v_in_unitig >> 1) as usize] >> 32) as u32)
            .wrapping_add(read_graph.sequences[(v_in_unitig >> 1) as usize].length)
            .wrapping_sub(1);
        x_dir = 0;
    } else if (arc.ul >> 32) as u32 == (v_in_unitig ^ 1) {
        // The arc enters at its start.
        x_pos = (position_index[(v_in_unitig >> 1) as usize] >> 32) as u32;
        x_dir = 1;
    } else {
        eprintln!("ERROR-get_xy_pos");
    }

    if arc.v == w_in_unitig {
        y_pos = ((position_index[(w_in_unitig >> 1) as usize] >> 32) as u32)
            .wrapping_add(arc.overlap)
            .wrapping_sub(1);
        y_dir = 0;
    } else if arc.v == (w_in_unitig ^ 1) {
        y_pos = ((position_index[(w_in_unitig >> 1) as usize] >> 32) as u32)
            .wrapping_add(read_graph.sequences[(w_in_unitig >> 1) as usize].length)
            .wrapping_sub(arc.overlap);
        y_dir = 1;
    } else {
        eprintln!("ERROR-get_xy_pos");
    }

    let rev = x_dir ^ y_dir;
    if rev != 0 && y_unitig_len > y_pos {
        y_pos = y_unitig_len - y_pos - 1;
    } else if rev != 0 {
        y_pos = u32::MAX;
    }

    if x_pos >= x_unitig_len {
        x_pos = u32::MAX;
    }
    if y_pos >= y_unitig_len {
        y_pos = u32::MAX;
    }

    ((u64::from(x_pos) << 32) | u64::from(y_pos), rev)
}

/// `deduplicate_edge` (`Purge_Dups.cpp:710`): compacts the trailing run of
/// arcs that share the same query read, merging duplicates with the same
/// `el` strand and the same `Cal_Off` offset.
pub fn deduplicate_edge(u_buffer: &mut [ArcOffset]) {
    if u_buffer.is_empty() {
        return;
    }
    let mut i = u_buffer.len() as i64 - 1;
    let v = (u_buffer[i as usize].arc.ul >> 33) as u32;

    while i >= 0 {
        if ((u_buffer[i as usize].arc.ul >> 33) as u32) != v {
            break;
        }
        i -= 1;
    }
    let run_start = (i + 1) as usize;

    let mut m = run_start;
    for i in run_start..u_buffer.len() {
        if u_buffer[i].arc.deleted {
            continue;
        }
        let i_off = get_x_offset(u_buffer[i].offset) - get_y_offset(u_buffer[i].offset);
        for k in (i + 1)..u_buffer.len() {
            if u_buffer[k].arc.deleted {
                continue;
            }
            if u_buffer[i].arc.exact_region != u_buffer[k].arc.exact_region {
                continue;
            }
            let k_off = get_x_offset(u_buffer[k].offset) - get_y_offset(u_buffer[k].offset);
            if i_off != k_off {
                continue;
            }
            u_buffer[k].arc.deleted = true;
            u_buffer[i].weight += u_buffer[k].weight;
        }
        u_buffer.swap(m, i);
        m += 1;
    }
    // `swap`-compaction keeps the reference's `a[m++] = a[i]` permutation.
}

/// `cmp_hap_candidates` (`Purge_Dups.cpp:781`): weight descending, then
/// `index_beg` ascending (uninitialized in the reference, see the module
/// notes).
pub fn compare_hap_candidates(left: &HapCandidate, right: &HapCandidate) -> std::cmp::Ordering {
    right
        .weight
        .cmp(&left.weight)
        .then_with(|| left.index_beg.cmp(&right.index_beg))
}

/// `get_hap_overlapLen` (`Purge_Dups.cpp:791`): aligns the left ends, extends
/// the right ends by the smaller remainder, and returns the adjusted
/// coordinates plus the overlap length.
#[allow(clippy::too_many_arguments)]
pub fn get_hap_overlap_len(
    x_beg: i64,
    x_end: i64,
    x_len: i64,
    y_beg: i64,
    y_end: i64,
    y_len: i64,
) -> (i64, i64, i64, i64, i64) {
    let (mut x_beg, mut x_end, mut y_beg, mut y_end) = (x_beg, x_end, y_beg, y_end);
    if x_beg <= y_beg {
        y_beg -= x_beg;
        x_beg = 0;
    } else {
        x_beg -= y_beg;
        y_beg = 0;
    }

    let x_right_length = x_len - x_end - 1;
    let y_right_length = y_len - y_end - 1;

    if x_right_length <= y_right_length {
        x_end = x_len - 1;
        y_end += x_right_length;
    } else {
        x_end += y_right_length;
        y_end = y_len - 1;
    }

    (x_beg, x_end, y_beg, y_end, x_end - x_beg + 1)
}

/// `classify_hap_overlap` (`Purge_Dups.cpp:831`): returns the overlap type
/// (`X2Y`/`Y2X`/`XCY`/`YCX`) and the adjusted coordinates.
#[allow(clippy::too_many_arguments)]
pub fn classify_hap_overlap(
    x_beg: i64,
    x_end: i64,
    x_len: i64,
    y_beg: i64,
    y_end: i64,
    y_len: i64,
) -> (u32, i64, i64, i64, i64) {
    let (n_x_beg, n_x_end, n_y_beg, n_y_end, _) =
        get_hap_overlap_len(x_beg, x_end, x_len, y_beg, y_end, y_len);
    let kind = if n_x_beg == 0 && n_x_end == x_len - 1 {
        YCX
    } else if n_y_beg == 0 && n_y_end == y_len - 1 {
        XCY
    } else if n_y_beg == 0 && n_x_end == x_len - 1 {
        X2Y
    } else if n_x_beg == 0 && n_y_end == y_len - 1 {
        Y2X
    } else {
        XCY
    };
    (kind, n_x_beg, n_x_end, n_y_beg, n_y_end)
}

/// `X2Y`.
pub const X2Y: u32 = 0;
/// `Y2X`.
pub const Y2X: u32 = 1;
/// `XCY`.
pub const XCY: u32 = 2;
/// `YCX`.
pub const YCX: u32 = 3;

/// `get_pair_purge_coverage` (`Purge_Dups.cpp:870`): the mean per-read
/// coverage of the reads covered by the two aligned intervals. Each side
/// keeps the first run of reads whose overlap with the interval spans at
/// least 80% of the read; the score is `sum(cov) / sum(len)` over both runs.
#[allow(clippy::too_many_arguments)]
pub fn get_pair_purge_coverage(
    x_reads: &MaUtg,
    x_pos_beg: i64,
    x_pos_end: i64,
    y_reads: &MaUtg,
    mut y_pos_beg: i64,
    mut y_pos_end: i64,
    rev: u32,
    read_graph: &Asg,
    cov: &HapCov,
) -> u64 {
    if rev != 0 {
        // Mirror the y interval and restore the ascending order.
        let beg = i64::from(y_reads.length) - y_pos_beg - 1;
        y_pos_end = i64::from(y_reads.length) - y_pos_end - 1;
        y_pos_beg = y_pos_end;
        y_pos_end = beg;
    }

    let mut u_cov = 0u64;
    let mut u_len = 0u64;

    for (reads, i_beg, i_end) in [
        (x_reads, x_pos_beg, x_pos_end),
        (y_reads, y_pos_beg, y_pos_end),
    ] {
        let mut idx_beg = -1i64;
        let mut idx_end = -1i64;
        let mut offset = 0i64;
        for i in 0..reads.reads.len() {
            let r_id = (reads.reads[i] >> 33) as usize;
            let r_beg = offset;
            let r_end = offset + i64::from(read_graph.sequences[r_id].length) - 1;
            offset += reads.reads[i] as u32 as i64;

            let ovlp = r_end.min(i_end) - r_beg.max(i_beg) + 1;
            if ovlp <= 0 || (ovlp as f64) < (read_graph.sequences[r_id].length as f64) * 0.8 {
                if idx_beg != -1 && idx_end != -1 {
                    break;
                }
                continue;
            }
            if idx_beg == -1 {
                idx_beg = i as i64;
            }
            idx_end = i as i64;
        }
        if idx_beg != -1 && idx_end != -1 {
            for i in idx_beg..=idx_end {
                let r_id = (reads.reads[i as usize] >> 33) as usize;
                u_cov += u64::from(cov.cov[r_id]);
                u_len += u64::from(cov.read_lengths[r_id]);
            }
        }
    }

    if u_len == 0 {
        0
    } else {
        u_cov / u_len
    }
}

/// `get_pair_hap_similarity` (`Purge_Dups.cpp:986`): how many of the reads in
/// `read_ids` map to `target_uid`, with the reference's stale-`j` cutoff.
///
/// The stale loop variable matters: a read with no reverse hits skips the
/// inner loop without resetting `j`, so the cutoff decision reuses the
/// previous read's exit value. That behaviour is replicated verbatim.
pub fn get_pair_hap_similarity(
    read_ids: &[u64],
    target_uid: u32,
    ctx: &HapAlignmentContext,
) -> (f64, f64) {
    let mut min_count = 0u32;
    let mut max_count = 0u32;
    let mut cutoff = 0u32;
    // The C `j` lives across iterations.
    let mut j: usize = 0;

    for i in 0..read_ids.len() {
        if cutoff > CUTOFF_THRES && cutoff > (read_ids.len() as u32) >> 1 {
            max_count = 0;
            min_count = read_ids.len() as u32;
            break;
        }
        let qn = (read_ids[i] >> 33) as usize;
        let hits = &ctx.reverse_sources[qn].hits;
        if !hits.is_empty() {
            min_count += 1;
        }
        if hits.is_empty() {
            continue;
        }
        j = 0;
        while j < hits.len() {
            let hit = &hits[j];
            let mut tn = hit.target_identifier;
            if ctx.read_graph.sequences[tn as usize].deleted {
                let (mapped, is_unitig) = ctx.read_to_unitig.get(tn);
                tn = mapped;
                if tn == u32::MAX
                    || is_unitig == 1
                    || ctx.read_graph.sequences[tn as usize].deleted
                {
                    j += 1;
                    continue;
                }
            }

            let (u_id, is_unitig) = ctx.read_to_unitig.get(tn);
            if u_id != u32::MAX && is_unitig == 1 && u_id == target_uid {
                max_count += 1;
                break;
            }
            j += 1;
        }

        if j == hits.len() {
            // No match; a zero-hit read leaves `j` stale, exactly like C.
            cutoff += 1;
        } else {
            cutoff = 0;
        }
    }

    (f64::from(max_count), f64::from(min_count))
}

/// `get_pair_hap_similarity_by_base` (`Purge_Dups.cpp:945`): the same vote
/// restricted to the reads covered by `[xBegPos, xEndPos]`.
#[allow(clippy::too_many_arguments)]
pub fn get_pair_hap_similarity_by_base(
    x_reads: &MaUtg,
    target_uid: u32,
    x_beg_pos: i64,
    x_end_pos: i64,
    ctx: &HapAlignmentContext,
) -> (f64, f64) {
    let mut min_count = 0u32;
    let mut max_count = 0u32;

    let mut offset = 0i64;
    for &packed in &x_reads.reads {
        let qn = (packed >> 33) as usize;
        let r_beg = offset;
        let r_end = offset + i64::from(ctx.read_graph.sequences[qn].length) - 1;
        offset += packed as u32 as i64;

        let ovlp = r_end.min(x_end_pos) - r_beg.max(x_beg_pos) + 1;
        if ovlp <= 0 {
            continue;
        }

        let hits = &ctx.reverse_sources[qn].hits;
        if !hits.is_empty() {
            min_count += 1;
        }
        if hits.is_empty() {
            continue;
        }
        for hit in hits {
            let mut tn = hit.target_identifier;
            if ctx.read_graph.sequences[tn as usize].deleted {
                let (mapped, is_unitig) = ctx.read_to_unitig.get(tn);
                tn = mapped;
                if tn == u32::MAX
                    || is_unitig == 1
                    || ctx.read_graph.sequences[tn as usize].deleted
                {
                    continue;
                }
            }

            let (u_id, is_unitig) = ctx.read_to_unitig.get(tn);
            if u_id != u32::MAX && is_unitig == 1 && u_id == target_uid {
                max_count += 1;
                break;
            }
        }
    }

    (f64::from(max_count), f64::from(min_count))
}

/// `check_hap_match` (`Purge_Dups.cpp:1089`).
#[allow(clippy::too_many_arguments)]
fn check_hap_match(
    qn: u32,
    target_beg: u32,
    target_end: u32,
    target_id: u32,
    ctx: &HapAlignmentContext,
) -> (u32, u32) {
    let mut is_found = 0u32;
    let mut is_match = 0u32;
    let hits = &ctx.reverse_sources[qn as usize].hits;
    if !hits.is_empty() {
        is_found = 1;
    }
    for hit in hits {
        let mut tn = hit.target_identifier;
        if ctx.read_graph.sequences[tn as usize].deleted {
            let (mapped, is_unitig) = ctx.read_to_unitig.get(tn);
            tn = mapped;
            if tn == u32::MAX
                || is_unitig == 1
                || ctx.read_graph.sequences[tn as usize].deleted
            {
                continue;
            }
        }

        let (u_id, is_unitig) = ctx.read_to_unitig.get(tn);
        if u_id != u32::MAX && is_unitig == 1 && u_id == target_id {
            let offset = ctx.position_index[tn as usize] as u32;
            if offset >= target_beg && offset <= target_end {
                is_match = 1;
                break;
            }
        }
    }
    (is_found, is_match)
}

/// `determin_hap_alignment_boundary_single_side`
/// (`Purge_Dups.cpp:1156`): walks the reads of one side, finds the last
/// position where the match ratio exceeds `H_rate` (or the best local score),
/// then shrinks the boundary back to the first fully matching read.
///
/// Returns `(match_len, max_count, min_count)`; the reference writes all
/// three as zeros when `query_len == 0` or `min_count == 0`.
#[allow(clippy::too_many_arguments)]
pub fn determin_hap_alignment_boundary_single_side(
    read_ids: &[u64],
    target_beg: i64,
    target_end: i64,
    target_id: u32,
    e_match: u32,
    e_total: u32,
    dir: u32,
    hap_rate: f32,
    is_local: bool,
    ctx: &HapAlignmentContext,
) -> (u32, u32, u32) {
    if read_ids.is_empty() {
        return (0, 0, 0);
    }
    let query_len = read_ids.len() as i64;
    let mut min_count = i64::from(e_total);
    let mut max_count = i64::from(e_match);
    let mut match_len = 0i64;
    let mut max_id = 0i64;
    let (mut score, mut max_score) = (0i64, 0i64);

    let target_beg = target_beg as u32;
    let target_end = target_end as u32;

    if dir == 0 {
        for i in 0..query_len {
            let (is_found, is_match) =
                check_hap_match((read_ids[i as usize] >> 33) as u32, target_beg, target_end, target_id, ctx);
            min_count += i64::from(is_found);
            max_count += i64::from(is_match);
            if (max_count as f32) > (min_count as f32) * hap_rate {
                max_id = i;
            }
            if is_local && is_found != 0 {
                let r_len = i64::from(ctx.read_graph.sequences[(read_ids[i as usize] >> 33) as usize].length);
                score += if is_match != 0 { r_len } else { -r_len };
                if score >= max_score {
                    max_score = score;
                    max_id = i;
                }
            }
        }

        let mut i = max_id;
        while i >= 0 {
            let (is_found, is_match) =
                check_hap_match((read_ids[i as usize] >> 33) as u32, target_beg, target_end, target_id, ctx);
            if is_found == 1 && is_match == 1 {
                break;
            }
            min_count -= i64::from(is_found);
            max_count -= i64::from(is_match);
            i -= 1;
        }
        match_len = i + 1;
    } else {
        let mut i = query_len - 1;
        max_id = query_len - 1;
        while i >= 0 {
            let (is_found, is_match) =
                check_hap_match((read_ids[i as usize] >> 33) as u32, target_beg, target_end, target_id, ctx);
            min_count += i64::from(is_found);
            max_count += i64::from(is_match);
            if (max_count as f32) > (min_count as f32) * hap_rate {
                max_id = i;
            }
            if is_local && is_found != 0 {
                let r_len = i64::from(ctx.read_graph.sequences[(read_ids[i as usize] >> 33) as usize].length);
                score += if is_match != 0 { r_len } else { -r_len };
                if score >= max_score {
                    max_score = score;
                    max_id = i;
                }
            }
            i -= 1;
        }

        let mut i = max_id;
        while (i as usize) < read_ids.len() {
            let (is_found, is_match) =
                check_hap_match((read_ids[i as usize] >> 33) as u32, target_beg, target_end, target_id, ctx);
            if is_found == 1 && is_match == 1 {
                break;
            }
            min_count -= i64::from(is_found);
            max_count -= i64::from(is_match);
            i += 1;
        }
        match_len = query_len - i;
    }

    if min_count == 0 {
        return (0, 0, 0);
    }

    (match_len as u32, max_count as u32, min_count as u32)
}

/// `modify_target_interval` (`Purge_Dups.cpp:1253`).
fn modify_target_interval(beg: i64, end: i64, len: i64) -> (i64, i64) {
    let mut beg = beg - TARGET_SHIFT;
    let mut end = end + TARGET_SHIFT;
    if beg < 0 {
        beg = 0;
    }
    if end >= len {
        end = len - 1;
    }
    (beg, end)
}

/// `bi_direction_hap_alignment_extention` (`Purge_Dups.cpp:1266`).
#[allow(clippy::too_many_arguments)]
fn bi_direction_hap_alignment_extention(
    x_reads: &MaUtg,
    mut x_left_beg: u32,
    mut x_left_len: u32,
    mut x_right_beg: u32,
    mut x_right_len: u32,
    target_uid: u32,
    target_beg: u32,
    target_end: u32,
    hap_rate: f32,
    is_local: bool,
    ctx: &HapAlignmentContext,
    rev: u32,
) -> (i64, i64) {
    if rev != 0 {
        std::mem::swap(&mut x_left_beg, &mut x_right_beg);
        std::mem::swap(&mut x_left_len, &mut x_right_len);
    }

    let (mut n_match_len_left, mut x_max_left, mut x_min_left) = (0u32, 0u32, 0u32);
    let (mut n_match_len_right, mut x_max_right, mut x_min_right) = (0u32, 0u32, 0u32);

    let left = &x_reads.reads[x_left_beg as usize..(x_left_beg + x_left_len) as usize];
    (n_match_len_left, x_max_left, x_min_left) = determin_hap_alignment_boundary_single_side(
        left,
        i64::from(target_beg),
        i64::from(target_end),
        target_uid,
        x_max_left,
        x_min_left,
        1,
        hap_rate,
        is_local,
        ctx,
    );

    let right = &x_reads.reads[x_right_beg as usize..(x_right_beg + x_right_len) as usize];
    (n_match_len_right, x_max_right, x_min_right) = determin_hap_alignment_boundary_single_side(
        right,
        i64::from(target_beg),
        i64::from(target_end),
        target_uid,
        x_max_right,
        x_min_right,
        0,
        hap_rate,
        is_local,
        ctx,
    );

    if x_max_left >= x_max_right {
        (n_match_len_right, x_max_right, x_min_right) = determin_hap_alignment_boundary_single_side(
            right,
            i64::from(target_beg),
            i64::from(target_end),
            target_uid,
            x_max_left,
            x_min_left,
            0,
            hap_rate,
            is_local,
            ctx,
        );
    } else {
        (n_match_len_left, x_max_left, x_min_left) = determin_hap_alignment_boundary_single_side(
            left,
            i64::from(target_beg),
            i64::from(target_end),
            target_uid,
            x_max_right,
            x_min_right,
            1,
            hap_rate,
            is_local,
            ctx,
        );
    }

    let beg = i64::from(x_left_beg + x_left_len) - i64::from(n_match_len_left);
    let end = i64::from(x_right_beg + n_match_len_right) - 1;
    (beg, end)
}

/// `get_hap_alignment_boundary` (`Purge_Dups.cpp:1310`): the four-type
/// dispatch over the per-side match counts.
///
/// Returns `(x_beg, x_end, y_beg, y_end)`; negative entries mean the caller
/// should reject the chain.
#[allow(clippy::too_many_arguments)]
pub fn get_hap_alignment_boundary(
    x_reads: &MaUtg,
    y_reads: &MaUtg,
    kind: u32,
    x_left_match: u32,
    x_left_total: u32,
    y_left_match: u32,
    y_left_total: u32,
    x_right_match: u32,
    x_right_total: u32,
    y_right_match: u32,
    y_right_total: u32,
    x_left_beg: u32,
    x_left_len: u32,
    y_left_beg: u32,
    y_left_len: u32,
    x_right_beg: u32,
    x_right_len: u32,
    y_right_beg: u32,
    y_right_len: u32,
    x_uid: u32,
    y_uid: u32,
    hap_rate: f32,
    is_local: bool,
    ctx: &HapAlignmentContext,
    rev: u32,
) -> (i64, i64, i64, i64) {
    let mut x_interval_beg;
    let mut x_interval_end;
    let mut y_interval_beg;
    let mut y_interval_end;

    let unitig_count = |reads: &MaUtg| reads.reads.len() as i64;

    if kind == X2Y {
        // x side
        let x_max_count = x_right_match;
        let x_min_count = x_right_total;
        let (target_beg, target_end) = modify_target_interval(
            i64::from(y_left_beg),
            i64::from(y_left_beg) + i64::from(y_left_len) - 1,
            unitig_count(y_reads),
        );
        let (n_match_len, _, _) = determin_hap_alignment_boundary_single_side(
            &x_reads.reads[x_left_beg as usize..(x_left_beg + x_left_len) as usize],
            target_beg,
            target_end,
            y_uid,
            x_max_count,
            x_min_count,
            1,
            hap_rate,
            is_local,
            ctx,
        );
        x_interval_beg = i64::from(x_left_beg + x_left_len);
        x_interval_beg -= i64::from(n_match_len);
        x_interval_end = i64::from(x_right_beg + x_right_len) - 1;

        // y side
        let y_max_count = y_left_match;
        let y_min_count = y_left_total;
        let (target_beg, target_end) = modify_target_interval(
            i64::from(x_right_beg),
            i64::from(x_right_beg) + i64::from(x_right_len) - 1,
            unitig_count(x_reads),
        );
        let (n_match_len, _, _) = determin_hap_alignment_boundary_single_side(
            &y_reads.reads[y_right_beg as usize..(y_right_beg + y_right_len) as usize],
            target_beg,
            target_end,
            x_uid,
            y_max_count,
            y_min_count,
            rev,
            hap_rate,
            is_local,
            ctx,
        );
        if rev == 0 {
            y_interval_beg = i64::from(y_left_beg);
            y_interval_end = i64::from(y_right_beg + n_match_len) - 1;
        } else {
            y_interval_beg = i64::from(y_right_beg + y_right_len);
            y_interval_beg -= i64::from(n_match_len);
            y_interval_end = i64::from(y_left_beg + y_left_len) - 1;
        }
    } else if kind == Y2X {
        // x side
        let x_max_count = x_left_match;
        let x_min_count = x_left_total;
        let (target_beg, target_end) = modify_target_interval(
            i64::from(y_right_beg),
            i64::from(y_right_beg) + i64::from(y_right_len) - 1,
            unitig_count(y_reads),
        );
        let (n_match_len, _, _) = determin_hap_alignment_boundary_single_side(
            &x_reads.reads[x_right_beg as usize..(x_right_beg + x_right_len) as usize],
            target_beg,
            target_end,
            y_uid,
            x_max_count,
            x_min_count,
            0,
            hap_rate,
            is_local,
            ctx,
        );
        x_interval_beg = i64::from(x_left_beg);
        x_interval_end = i64::from(x_right_beg + n_match_len) - 1;

        // y side
        let y_max_count = y_right_match;
        let y_min_count = y_right_total;
        let (target_beg, target_end) = modify_target_interval(
            i64::from(x_left_beg),
            i64::from(x_left_beg) + i64::from(x_left_len) - 1,
            unitig_count(x_reads),
        );
        let (n_match_len, _, _) = determin_hap_alignment_boundary_single_side(
            &y_reads.reads[y_left_beg as usize..(y_left_beg + y_left_len) as usize],
            target_beg,
            target_end,
            x_uid,
            y_max_count,
            y_min_count,
            1 - rev,
            hap_rate,
            is_local,
            ctx,
        );
        if rev == 0 {
            y_interval_beg = i64::from(y_left_beg + y_left_len);
            y_interval_beg -= i64::from(n_match_len);
            y_interval_end = i64::from(y_right_beg + y_right_len) - 1;
        } else {
            y_interval_beg = i64::from(y_right_beg);
            y_interval_end = i64::from(y_left_beg + n_match_len) - 1;
        }
    } else if kind == XCY {
        (x_interval_beg, x_interval_end) = bi_direction_hap_alignment_extention(
            x_reads,
            x_left_beg,
            x_left_len,
            x_right_beg,
            x_right_len,
            y_uid,
            0,
            (y_reads.reads.len() - 1) as u32,
            hap_rate,
            is_local,
            ctx,
            0,
        );
        y_interval_beg = 0;
        y_interval_end = unitig_count(y_reads) - 1;
    } else if kind == YCX {
        x_interval_beg = 0;
        x_interval_end = unitig_count(x_reads) - 1;
        (y_interval_beg, y_interval_end) = bi_direction_hap_alignment_extention(
            y_reads,
            y_left_beg,
            y_left_len,
            y_right_beg,
            y_right_len,
            x_uid,
            0,
            (x_reads.reads.len() - 1) as u32,
            hap_rate,
            is_local,
            ctx,
            rev,
        );
    } else {
        panic!("classify_hap_overlap produced an unknown type");
    }

    (x_interval_beg, x_interval_end, y_interval_beg, y_interval_end)
}

/// `get_idx_by_base` (`Purge_Dups.cpp:2114`): the read-index range covered by
/// a base interval of a unitig.
fn get_idx_by_base(
    reads: &MaUtg,
    read_graph: &Asg,
    beg_base: i64,
    end_base: i64,
) -> (i64, i64) {
    let mut beg_idx = -1i64;
    let mut end_idx = -1i64;
    let mut offset = 0i64;
    for i in 0..reads.reads.len() {
        let r_id = (reads.reads[i] >> 33) as usize;
        let r_beg = offset;
        let r_end = offset + i64::from(read_graph.sequences[r_id].length) - 1;
        offset += reads.reads[i] as u32 as i64;
        if beg_base > r_end || r_beg > end_base {
            if beg_idx != -1 && end_idx != -1 {
                break;
            }
            continue;
        }
        if beg_idx == -1 {
            beg_idx = i as i64;
        }
        end_idx = i as i64;
    }
    (beg_idx, end_idx)
}

/// `get_base_boundary_chain` (`Purge_Dups.cpp:2135`): collects the arcs
/// between the two index ranges, merges equal offsets and returns whether the
/// longest increasing subsequence is non-empty.
///
/// The `u_buffer`/`tail_index`/`prev_index` scratch is reused across the two
/// rounds exactly like the reference.
#[allow(clippy::too_many_arguments)]
fn get_base_boundary_chain(
    ctx: &HapAlignmentContext,
    x_reads: &MaUtg,
    y_reads: &MaUtg,
    x_uid: u32,
    y_uid: u32,
    x_beg_index: i64,
    x_end_index: i64,
    y_beg_index: i64,
    y_end_index: i64,
    rev: u32,
    u_buffer: &mut Vec<ArcOffset>,
    tail_index: &mut Vec<i32>, prev_index: &mut Vec<i32>,
) -> bool {
    u_buffer.clear();
    for k in x_beg_index..=x_end_index {
        let k = k as usize;
        let x_read = (x_reads.reads[k] >> 33) as usize;
        for hit in &ctx.reverse_sources[x_read].hits {
            let Some(mut arc) = scan_arc(ctx, hit, y_uid) else {
                continue;
            };

            let v = (x_reads.reads[k] >> 32) as u32;
            let (v_unitig, is_unitig) = ctx.read_to_unitig.get(v >> 1);
            if is_unitig == 0 || v_unitig == u32::MAX || v_unitig != x_uid {
                continue;
            }
            if ctx.position_index[(v >> 1) as usize] as u32 != k as u32 {
                continue;
            }

            let w = (y_reads.reads[(ctx.position_index[(arc.v >> 1) as usize] as u32) as usize] >> 32) as u32;

            let v_dir = if (arc.ul >> 32) as u32 == v { 1u32 } else { 0u32 };
            let w_dir = if arc.v == w { 1u32 } else { 0u32 };
            if rev == 0 && v_dir != w_dir {
                continue;
            }
            if rev == 1 && v_dir == w_dir {
                continue;
            }

            // "may have bugs" in the reference: the offset window check.
            let offset = ctx.position_index[(arc.v >> 1) as usize] as u32 as i64;
            if offset < y_beg_index || offset > y_end_index {
                continue;
            }

            let (tmp, arc_rev) = get_xy_pos(
                ctx.read_graph,
                &arc,
                v,
                w,
                x_reads.length,
                y_reads.length,
                ctx.position_index,
            );
            if (tmp >> 32) == u32::MAX as u64 || (tmp as u32) == u32::MAX {
                continue;
            }
            arc.exact_region = arc_rev;

            u_buffer.push(ArcOffset {
                offset: tmp,
                arc,
                weight: 1,
            });
        }
    }
    if u_buffer.is_empty() {
        return false;
    }

    u_buffer.sort_by(compare_arc_offsets);

    let mut l = 0usize;
    let mut m = 0usize;
    let mut k = 1usize;
    while k <= u_buffer.len() {
        if k == u_buffer.len() || u_buffer[k].offset != u_buffer[l].offset {
            u_buffer[m] = u_buffer[l].clone();
            l += 1;
            while l < k {
                u_buffer[m].weight += u_buffer[l].weight;
                if u_buffer[l].arc.overlap > u_buffer[m].arc.overlap {
                    u_buffer[m].arc = u_buffer[l].arc.clone();
                }
                l += 1;
            }
            l = k;
            m += 1;
        }
        k += 1;
    }
    u_buffer.truncate(m);

    quick_lis(u_buffer, tail_index, prev_index);
    !tail_index.is_empty()
}

/// `get_base_boundary_advance` (`Purge_Dups.cpp:2233`): two rounds of
/// [`get_base_boundary_chain`], the first rounded through base coordinates.
#[allow(clippy::too_many_arguments)]
fn get_base_boundary_advance(
    ctx: &HapAlignmentContext,
    x_reads: &MaUtg,
    y_reads: &MaUtg,
    x_uid: u32,
    y_uid: u32,
    x_beg_index: i64,
    x_end_index: i64,
    y_beg_index: i64,
    y_end_index: i64,
    rev: u32,
    u_buffer: &mut Vec<ArcOffset>,
    tail_index: &mut Vec<i32>, prev_index: &mut Vec<i32>,
) -> (u32, u32, u32, u32) {
    let invalid = (u32::MAX, u32::MAX, u32::MAX, u32::MAX);
    if !get_base_boundary_chain(
        ctx,
        x_reads,
        y_reads,
        x_uid,
        y_uid,
        x_beg_index,
        x_end_index,
        y_beg_index,
        y_end_index,
        rev,
        u_buffer,
        tail_index,
        prev_index,
    ) {
        return invalid;
    }

    // Base coordinates of the LIS endpoints.
    let mut new_x_beg = get_x_offset(u_buffer[tail_index[0] as usize].offset);
    let mut new_y_beg = get_y_offset(u_buffer[tail_index[0] as usize].offset);
    let mut new_x_end = get_x_offset(u_buffer[tail_index[tail_index.len() - 1] as usize].offset);
    let mut new_y_end = get_y_offset(u_buffer[tail_index[tail_index.len() - 1] as usize].offset);

    if new_x_beg > new_x_end || new_y_beg > new_y_end {
        return invalid;
    }

    let (_, x_beg, x_end, y_beg, y_end) = classify_hap_overlap(
        new_x_beg,
        new_x_end,
        i64::from(x_reads.length),
        new_y_beg,
        new_y_end,
        i64::from(y_reads.length),
    );
    new_x_beg = x_beg;
    new_x_end = x_end;
    new_y_beg = y_beg;
    new_y_end = y_end;

    if rev != 0 {
        let len = i64::from(y_reads.length);
        new_y_beg = len - new_y_beg - 1;
        new_y_end = len - new_y_end - 1;
        std::mem::swap(&mut new_y_beg, &mut new_y_end);
    }

    // Index coordinates of the rounded interval.
    let (new_x_idx_beg, new_x_idx_end) = get_idx_by_base(x_reads, ctx.read_graph, new_x_beg, new_x_end);
    let (new_y_idx_beg, new_y_idx_end) = get_idx_by_base(y_reads, ctx.read_graph, new_y_beg, new_y_end);
    if new_x_idx_beg == -1 || new_x_idx_end == -1 || new_y_idx_beg == -1 || new_y_idx_end == -1 {
        return invalid;
    }

    if !get_base_boundary_chain(
        ctx,
        x_reads,
        y_reads,
        x_uid,
        y_uid,
        new_x_idx_beg,
        new_x_idx_end,
        new_y_idx_beg,
        new_y_idx_end,
        rev,
        u_buffer,
        tail_index,
        prev_index,
    ) {
        return invalid;
    }

    let x_beg = u_buffer[tail_index[0] as usize].offset >> 32;
    let y_beg = u_buffer[tail_index[0] as usize].offset as u32;
    let x_end = (u_buffer[tail_index[tail_index.len() - 1] as usize].offset >> 32) as u32;
    let y_end = u_buffer[tail_index[tail_index.len() - 1] as usize].offset as u32;
    (x_beg as u32, x_end, y_beg, y_end)
}

/// `get_base_boundary_advance_back` (`Purge_Dups.cpp:1983`): the single-round
/// variant used by `determine_hap_overlap_type_advance_back`.
#[allow(clippy::too_many_arguments)]
fn get_base_boundary_advance_back(
    ctx: &HapAlignmentContext,
    x_reads: &MaUtg,
    y_reads: &MaUtg,
    x_uid: u32,
    y_uid: u32,
    x_beg_index: i64,
    x_end_index: i64,
    y_beg_index: i64,
    y_end_index: i64,
    rev: u32,
    u_buffer: &mut Vec<ArcOffset>,
    tail_index: &mut Vec<i32>, prev_index: &mut Vec<i32>,
) -> (u32, u32, u32, u32) {
    let invalid = (u32::MAX, u32::MAX, u32::MAX, u32::MAX);
    u_buffer.clear();
    for k in x_beg_index..=x_end_index {
        let k = k as usize;
        let x_read = (x_reads.reads[k] >> 33) as usize;
        for hit in &ctx.reverse_sources[x_read].hits {
            let Some(mut arc) = scan_arc(ctx, hit, y_uid) else {
                continue;
            };

            let v = (x_reads.reads[k] >> 32) as u32;
            let (v_unitig, is_unitig) = ctx.read_to_unitig.get(v >> 1);
            if is_unitig == 0 || v_unitig == u32::MAX || v_unitig != x_uid {
                continue;
            }
            if ctx.position_index[(v >> 1) as usize] as u32 != k as u32 {
                continue;
            }

            let w = (y_reads.reads[(ctx.position_index[(arc.v >> 1) as usize] as u32) as usize] >> 32) as u32;

            let v_dir = if (arc.ul >> 32) as u32 == v { 1u32 } else { 0u32 };
            let w_dir = if arc.v == w { 1u32 } else { 0u32 };
            if rev == 0 && v_dir != w_dir {
                continue;
            }
            if rev == 1 && v_dir == w_dir {
                continue;
            }

            let offset = ctx.position_index[(arc.v >> 1) as usize] as u32 as i64;
            if offset < y_beg_index || offset > y_end_index {
                continue;
            }

            let (tmp, arc_rev) = get_xy_pos(
                ctx.read_graph,
                &arc,
                v,
                w,
                x_reads.length,
                y_reads.length,
                ctx.position_index,
            );
            if (tmp >> 32) == u32::MAX as u64 || (tmp as u32) == u32::MAX {
                continue;
            }
            arc.exact_region = arc_rev;

            u_buffer.push(ArcOffset {
                offset: tmp,
                arc,
                weight: 1,
            });
        }
    }
    if u_buffer.is_empty() {
        return invalid;
    }

    u_buffer.sort_by(compare_arc_offsets);

    let mut m = 1usize;
    for k in 1..u_buffer.len() {
        if u_buffer[m - 1].offset == u_buffer[k].offset {
            u_buffer[m - 1].weight += u_buffer[k].weight;
            if u_buffer[k].arc.overlap > u_buffer[m - 1].arc.overlap {
                u_buffer[m - 1].arc = u_buffer[k].arc.clone();
            }
            continue;
        }
        u_buffer[m] = u_buffer[k].clone();
        m += 1;
    }
    u_buffer.truncate(m);

    quick_lis(u_buffer, tail_index, prev_index);
    if tail_index.is_empty() {
        return invalid;
    }

    let x_beg = (u_buffer[tail_index[0] as usize].offset >> 32) as u32;
    let y_beg = u_buffer[tail_index[0] as usize].offset as u32;
    let last = tail_index[tail_index.len() - 1] as usize;
    let x_end = (u_buffer[last].offset >> 32) as u32;
    let y_end = u_buffer[last].offset as u32;
    (x_beg, x_end, y_beg, y_end)
}

/// `get_chain_score` (`Purge_Dups.cpp:2289`): the chain score against the
/// unchained reverse-covered reads of `x`.
#[allow(clippy::too_many_arguments)]
fn get_chain_score(
    x_reads: &MaUtg,
    u_buffer: &[ArcOffset],
    tail_index: &[i32],
    ctx: &HapAlignmentContext,
    x_beg_pos: i64,
    x_end_pos: i64,
) -> i64 {
    // `idx`: the sorted chain read indices inside `x`.
    let mut idx: Vec<i32> = Vec::with_capacity(tail_index.len());
    let mut k = 0usize;
    for &t in tail_index {
        let r_id = (u_buffer[t as usize].arc.ul >> 33) as usize;

        while k < x_reads.reads.len() {
            if r_id == (x_reads.reads[k] >> 33) as usize {
                break;
            }
            k += 1;
        }
        if k >= x_reads.reads.len() {
            k = 0;
            while k < x_reads.reads.len() {
                if r_id == (x_reads.reads[k] >> 33) as usize {
                    break;
                }
                k += 1;
            }
        }

        if k < x_reads.reads.len() {
            idx.push(k as i32);
        } else {
            eprintln!("\nERROR-get_chain_score");
        }
    }
    idx.sort_unstable();

    let mut inp_beg = -1i64;
    let mut inp_end = -2i64;
    let mut hap_beg = -1i64;
    let mut hap_end = -2i64;
    let (mut inp_match, mut hap_match) = (0i64, 0i64);
    let mut offset = 0i64;
    let mut k = 0usize;

    for i in 0..x_reads.reads.len() {
        let r_id = (x_reads.reads[i] >> 33) as usize;
        let r_beg = offset;
        let r_end = offset + i64::from(ctx.read_graph.sequences[r_id].length) - 1;
        offset += x_reads.reads[i] as u32 as i64;

        if !ctx.reverse_sources[r_id].hits.is_empty() {
            if r_beg <= hap_end {
                hap_end = hap_end.max(r_end);
            } else {
                let ovlp = hap_end.min(x_end_pos) - hap_beg.max(x_beg_pos) + 1;
                hap_match += ovlp.max(0);
                hap_beg = r_beg;
                hap_end = r_end;
            }
        }

        while k < idx.len() && i > idx[k] as usize {
            k += 1;
        }
        if k >= idx.len() {
            continue;
        }
        if i == idx[k] as usize {
            if r_beg <= inp_end {
                inp_end = inp_end.max(r_end);
            } else {
                let ovlp = inp_end.min(x_end_pos) - inp_beg.max(x_beg_pos) + 1;
                inp_match += ovlp.max(0);
                inp_beg = r_beg;
                inp_end = r_end;
            }
        }
    }

    let ovlp = inp_end.min(x_end_pos) - inp_beg.max(x_beg_pos) + 1;
    inp_match += ovlp.max(0);
    let ovlp = hap_end.min(x_end_pos) - hap_beg.max(x_beg_pos) + 1;
    hap_match += ovlp.max(0);

    ((inp_match as f64) * 1.0 - ((hap_match - inp_match) as f64) * CHAIN_UNMATCH) as i64
}

/// `determine_hap_overlap_type_advance` (`Purge_Dups.cpp:2384`).
///
/// On success returns the overlap type and fills the candidate's position,
/// weight, `index_beg`, `index_end` and `score` fields; the four position
/// outputs carry the base-coordinate interval from `classify_hap_overlap`.
#[allow(clippy::too_many_arguments)]
fn determine_hap_overlap_type_advance(
    hap_can: &mut HapCandidate,
    x_reads: &MaUtg,
    y_reads: &MaUtg,
    x_uid: u32,
    y_uid: u32,
    ctx: &HapAlignmentContext,
    u_buffer: &mut Vec<ArcOffset>,
    tail_index: &mut Vec<i32>, prev_index: &mut Vec<i32>,
) -> (u32, i64, i64, i64, i64) {
    let (x_pos_beg, x_pos_end, y_pos_beg, y_pos_end) = get_base_boundary_advance(
        ctx,
        x_reads,
        y_reads,
        x_uid,
        y_uid,
        i64::from(hap_can.x_beg_pos),
        i64::from(hap_can.x_end_pos),
        i64::from(hap_can.y_beg_pos),
        i64::from(hap_can.y_end_pos),
        u32::from(hap_can.rev),
        u_buffer,
        tail_index,
        prev_index,
    );
    if x_pos_beg == u32::MAX
        || y_pos_beg == u32::MAX
        || x_pos_end == u32::MAX
        || y_pos_end == u32::MAX
    {
        return (u32::MAX, 0, 0, 0, 0);
    }
    if x_pos_beg > x_pos_end || y_pos_beg > y_pos_end {
        return (u32::MAX, 0, 0, 0, 0);
    }

    let (kind, r_x_beg, r_x_end, r_y_beg, r_y_end) = classify_hap_overlap(
        i64::from(x_pos_beg),
        i64::from(x_pos_end),
        i64::from(x_reads.length),
        i64::from(y_pos_beg),
        i64::from(y_pos_end),
        i64::from(y_reads.length),
    );
    hap_can.index_end = kind;

    let first = tail_index[0] as usize;
    let last = tail_index[tail_index.len() - 1] as usize;
    let x_first = (u_buffer[first].arc.ul >> 33) as usize;
    let x_last = (u_buffer[last].arc.ul >> 33) as usize;
    let y_first = (u_buffer[first].arc.v >> 1) as usize;
    let y_last = (u_buffer[last].arc.v >> 1) as usize;
    hap_can.x_beg_pos = ctx.position_index[x_first] as u32;
    hap_can.x_end_pos = ctx.position_index[x_last] as u32;
    if hap_can.x_beg_pos > hap_can.x_end_pos {
        std::mem::swap(&mut hap_can.x_beg_pos, &mut hap_can.x_end_pos);
    }
    hap_can.y_beg_pos = ctx.position_index[y_first] as u32;
    hap_can.y_end_pos = ctx.position_index[y_last] as u32;
    if hap_can.y_beg_pos > hap_can.y_end_pos {
        std::mem::swap(&mut hap_can.y_beg_pos, &mut hap_can.y_end_pos);
    }

    let span = hap_can.x_end_pos as usize - hap_can.x_beg_pos as usize + 1;
    let (x_left_match, x_left_total) = get_pair_hap_similarity(
        &x_reads.reads[hap_can.x_beg_pos as usize..hap_can.x_beg_pos as usize + span],
        y_uid,
        ctx,
    );
    if x_left_match == 0.0 || x_left_total == 0.0 {
        return (u32::MAX, 0, 0, 0, 0);
    }
    hap_can.weight = x_left_match as u64;
    hap_can.index_beg = x_left_total as u32;
    hap_can.score = get_chain_score(x_reads, u_buffer, &tail_index, ctx, r_x_beg, r_x_end);
    if hap_can.score <= 0 {
        return (u32::MAX, 0, 0, 0, 0);
    }
    (kind, r_x_beg, r_x_end, r_y_beg, r_y_end)
}

/// `determine_hap_overlap_type_advance_back` (`Purge_Dups.cpp:2084`).
#[allow(clippy::too_many_arguments)]
fn determine_hap_overlap_type_advance_back(
    hap_can: &mut HapCandidate,
    x_reads: &MaUtg,
    y_reads: &MaUtg,
    x_uid: u32,
    y_uid: u32,
    ctx: &HapAlignmentContext,
    u_buffer: &mut Vec<ArcOffset>,
    tail_index: &mut Vec<i32>, prev_index: &mut Vec<i32>,
) -> (u32, i64, i64, i64, i64) {
    let (x_pos_beg, x_pos_end, y_pos_beg, y_pos_end) = get_base_boundary_advance_back(
        ctx,
        x_reads,
        y_reads,
        x_uid,
        y_uid,
        i64::from(hap_can.x_beg_pos),
        i64::from(hap_can.x_end_pos),
        i64::from(hap_can.y_beg_pos),
        i64::from(hap_can.y_end_pos),
        u32::from(hap_can.rev),
        u_buffer,
        tail_index,
        prev_index,
    );
    if x_pos_beg == u32::MAX
        || y_pos_beg == u32::MAX
        || x_pos_end == u32::MAX
        || y_pos_end == u32::MAX
    {
        return (u32::MAX, 0, 0, 0, 0);
    }
    if x_pos_beg > x_pos_end || y_pos_beg > y_pos_end {
        return (u32::MAX, 0, 0, 0, 0);
    }

    let (kind, r_x_beg, r_x_end, r_y_beg, r_y_end) = classify_hap_overlap(
        i64::from(x_pos_beg),
        i64::from(x_pos_end),
        i64::from(x_reads.length),
        i64::from(y_pos_beg),
        i64::from(y_pos_end),
        i64::from(y_reads.length),
    );
    (kind, r_x_beg, r_x_end, r_y_beg, r_y_end)
}

/// `calculate_pair_hap_similarity_advance` (`Purge_Dups.cpp:2479`).
///
/// Returns `Ok(true)` (`PLOID`) with the final coordinates and similarity, or
/// `Ok(false)` (`NON_PLOID`).
#[allow(clippy::too_many_arguments)]
pub fn calculate_pair_hap_similarity_advance(
    hap_can: &mut HapCandidate,
    x_uid: u32,
    y_uid: u32,
    x_reads: &MaUtg,
    y_reads: &MaUtg,
    ctx: &HapAlignmentContext,
    hap_rate: f32,
    is_local: bool,
    cov_threshold: u64,
    u_buffer: &mut Vec<ArcOffset>,
    tail_index: &mut Vec<i32>,
    prev_index: &mut Vec<i32>,
    cov: &HapCov,
) -> (bool, i64, i64, i64, i64, f32) {
    let x_len = x_reads.reads.len() as u32;
    let y_len = y_reads.reads.len() as u32;

    let anchor_arc = &hap_can.arc;
    let x_index = ctx.position_index[(anchor_arc.ul >> 33) as usize] as u32;
    let y_index = ctx.position_index[(anchor_arc.v >> 1) as usize] as u32;

    let (x_left_len, x_right_beg, x_right_len, y_left_len, y_right_beg, y_right_len) =
        if hap_can.rev == 0 {
            (
                x_index,
                x_index,
                x_len - x_index,
                y_index,
                y_index,
                y_len - y_index,
            )
        } else {
            (
                x_index,
                x_index,
                x_len - x_index,
                y_len - (y_index + 1),
                0,
                y_index + 1,
            )
        };
    let x_left_beg = 0u32;
    let y_left_beg = if hap_can.rev == 0 { 0u32 } else { y_index + 1 };

    let kind = hap_can.index_end;

    let (mut x_left_match, mut x_left_total, mut y_left_match, mut y_left_total) = (0.0, 0.0, 0.0, 0.0);
    let (mut x_right_match, mut x_right_total, mut y_right_match, mut y_right_total) = (0.0, 0.0, 0.0, 0.0);
    let mut max_count: u32;
    let mut min_count: u32;

    if kind == XCY {
        let (m, t) = get_pair_hap_similarity(
            &y_reads.reads[..y_len as usize],
            x_uid,
            ctx,
        );
        y_left_match = m;
        y_left_total = t;
        max_count = y_left_match as u32;
        min_count = y_left_total as u32;
    } else if kind == YCX {
        let (m, t) = get_pair_hap_similarity(
            &x_reads.reads[..x_len as usize],
            y_uid,
            ctx,
        );
        x_left_match = m;
        x_left_total = t;
        max_count = x_left_match as u32;
        min_count = x_left_total as u32;
    } else if kind == X2Y {
        let (m, t) = get_pair_hap_similarity(
            &y_reads.reads[y_left_beg as usize..(y_left_beg + y_left_len) as usize],
            x_uid,
            ctx,
        );
        y_left_match = m;
        y_left_total = t;
        let (m2, t2) = get_pair_hap_similarity(
            &x_reads.reads[x_right_beg as usize..(x_right_beg + x_right_len) as usize],
            y_uid,
            ctx,
        );
        x_right_match = m2;
        x_right_total = t2;
        max_count = (y_left_match + x_right_match) as u32;
        min_count = (y_left_total + x_right_total) as u32;
    } else if kind == Y2X {
        let (m, t) = get_pair_hap_similarity(
            &x_reads.reads[x_left_beg as usize..(x_left_beg + x_left_len) as usize],
            y_uid,
            ctx,
        );
        x_left_match = m;
        x_left_total = t;
        let (m2, t2) = get_pair_hap_similarity(
            &y_reads.reads[y_right_beg as usize..(y_right_beg + y_right_len) as usize],
            x_uid,
            ctx,
        );
        y_right_match = m2;
        y_right_total = t2;
        max_count = (x_left_match + y_right_match) as u32;
        min_count = (x_left_total + y_right_total) as u32;
    } else {
        panic!("calculate_pair_hap_similarity_advance: unknown type");
    }

    hap_can.weight = 0;
    hap_can.index_beg = 0;
    if min_count == 0 {
        return (false, 0, 0, 0, 0, 0.0);
    }
    // `float` comparison: `max_count > min_count*Hap_rate` in the reference
    // converts both sides to `float` (the totals are `double` there).
    if (max_count as f32) > (min_count as f32) * hap_rate || is_local {
        let (r_x_interval_beg, r_x_interval_end, r_y_interval_beg, r_y_interval_end) =
            get_hap_alignment_boundary(
                x_reads,
                y_reads,
                kind,
                x_left_match as u32,
                x_left_total as u32,
                y_left_match as u32,
                y_left_total as u32,
                x_right_match as u32,
                x_right_total as u32,
                y_right_match as u32,
                y_right_total as u32,
                x_left_beg,
                x_left_len,
                y_left_beg,
                y_left_len,
                x_right_beg,
                x_right_len,
                y_right_beg,
                y_right_len,
                x_uid,
                y_uid,
                hap_rate,
                is_local,
                ctx,
                u32::from(hap_can.rev),
            );

        if r_x_interval_beg < 0
            || r_x_interval_end < 0
            || r_y_interval_beg < 0
            || r_y_interval_end < 0
        {
            return (false, 0, 0, 0, 0, 0.0);
        }

        let span = (r_x_interval_end + 1 - r_x_interval_beg) as usize;
        let (match1, total1) = get_pair_hap_similarity(
            &x_reads.reads[r_x_interval_beg as usize..r_x_interval_beg as usize + span],
            y_uid,
            ctx,
        );
        if match1 == 0.0
            || total1 == 0.0
            || (!is_local && match1 <= total1 * f64::from(hap_rate))
        {
            return (false, 0, 0, 0, 0, 0.0);
        }

        hap_can.weight = match1 as u64;
        hap_can.index_beg = total1 as u32;
        hap_can.index_end = kind;
        hap_can.x_beg_pos = r_x_interval_beg as u32;
        hap_can.x_end_pos = r_x_interval_end as u32;
        hap_can.y_beg_pos = r_y_interval_beg as u32;
        hap_can.y_end_pos = r_y_interval_end as u32;
        let (kind2, r_x_beg, r_x_end, r_y_beg, r_y_end) = determine_hap_overlap_type_advance(
            hap_can,
            x_reads,
            y_reads,
            x_uid,
            y_uid,
            ctx,
            u_buffer,
            tail_index,
        prev_index,
        );
        if kind2 == XCY && y_reads.length > x_reads.length * 2 {
            return (false, 0, 0, 0, 0, 0.0);
        }
        if kind2 == YCX && x_reads.length > y_reads.length * 2 {
            return (false, 0, 0, 0, 0, 0.0);
        }
        if kind2 == u32::MAX {
            return (false, 0, 0, 0, 0, 0.0);
        }

        let ploid_coverage = get_pair_purge_coverage(
            x_reads,
            r_x_beg,
            r_x_end,
            y_reads,
            r_y_beg,
            r_y_end,
            u32::from(hap_can.rev),
            ctx.read_graph,
            cov,
        );
        if cov_threshold > 0 && ploid_coverage >= cov_threshold {
            return (false, 0, 0, 0, 0, 0.0);
        }

        let (match2, total2) = get_pair_hap_similarity_by_base(
            x_reads,
            y_uid,
            r_x_beg,
            r_x_end,
            ctx,
        );
        let sim = if total2 == 0.0 {
            0.0
        } else {
            match2 / total2
        };

        if match2 == 0.0 || total2 == 0.0 || (sim as f32) <= hap_rate {
            return (false, 0, 0, 0, 0, 0.0);
        }

        return (true, r_x_beg, r_x_end, r_y_beg, r_y_end, sim as f32);
    }
    (false, 0, 0, 0, 0, 0.0)
}

/// `get_max_index` (`Purge_Dups.cpp:2623`).
fn get_max_index(
    x: &[ArcOffset],
    scores: &[i32],
    flags: &[u8],
    x_read_len: i64,
    y_read_len: i64,
) -> i64 {
    let mut max_result = -1i64;
    let mut max_i = -1i64;
    let mut min_x_len = x_read_len * 2 + 2;

    for i in 0..x.len() {
        if flags[i] != 0 {
            continue;
        }
        let x_off = (x[i].offset >> 32) as i64;
        let y_off = (x[i].offset as u32) as i64;
        if i64::from(scores[i]) > max_result {
            max_result = i64::from(scores[i]);
            max_i = i as i64;
            let (_, _, _, _, span) =
                get_hap_overlap_len(x_off, x_off, x_read_len, y_off, y_off, y_read_len);
            min_x_len = span;
        } else if i64::from(scores[i]) == max_result {
            let (_, _, _, _, span) =
                get_hap_overlap_len(x_off, x_off, x_read_len, y_off, y_off, y_read_len);
            if span < min_x_len {
                max_result = i64::from(scores[i]);
                max_i = i as i64;
                min_x_len = span;
            }
        }
    }

    max_i
}

/// `push_hap_can` (`Purge_Dups.cpp:2673`).
fn push_hap_can(
    x: &[ArcOffset],
    u_can: &mut Vec<HapCandidate>,
    results: &[i32],
    x_read_len: i64,
    y_read_len: i64,
) {
    let chain_len = results.len();
    if chain_len == 0 {
        return;
    }
    let mut hap_can = HapCandidate::default();
    hap_can.rev = x[results[0] as usize].arc.exact_region;
    hap_can.x_end_pos = (x[results[0] as usize].offset >> 32) as u32;
    hap_can.x_beg_pos = hap_can.x_end_pos;
    hap_can.y_end_pos = x[results[0] as usize].offset as u32;
    hap_can.y_beg_pos = hap_can.y_end_pos;
    hap_can.weight = 0;

    for &r in results {
        let r = r as usize;
        hap_can.x_beg_pos = (x[r].offset >> 32) as u32;
        hap_can.y_beg_pos = x[r].offset as u32;
        hap_can.weight += x[r].weight;
    }

    let mut i = 0usize;
    let mut total_weight = 0u64;
    for &r in results {
        let r = r as usize;
        total_weight += x[r].weight;
        if total_weight >= hap_can.weight / 2 {
            break;
        }
        i += 1;
    }

    if i >= chain_len {
        i = chain_len - 1;
    }
    hap_can.arc = x[results[i] as usize].arc.clone();
    let (kind, _, _, _, _) = classify_hap_overlap(
        i64::from(hap_can.x_beg_pos),
        i64::from(hap_can.x_end_pos),
        x_read_len,
        i64::from(hap_can.y_beg_pos),
        i64::from(hap_can.y_end_pos),
        y_read_len,
    );
    hap_can.index_end = kind;
    u_can.push(hap_can);
}

/// `hap_chaining` (`Purge_Dups.cpp:2726`): the DP chain extraction with the
/// duplicate-overlap restart loop.
#[allow(clippy::too_many_arguments)]
pub fn hap_chaining(
    x: &[ArcOffset],
    scores: &mut Vec<i32>,
    prev_index: &mut Vec<i32>,
    beg_index: &mut Vec<i32>,
    flags: &mut Vec<u8>,
    band_width_threshold: f32,
    max_skip: i64,
    x_read_len: i64,
    y_read_len: i64,
    u_can: &mut Vec<HapCandidate>,
) {
    scores.clear();
    prev_index.clear();
    beg_index.clear();
    flags.clear();
    if x.is_empty() {
        return;
    }
    let n = x.len();
    scores.resize(n, 0);
    prev_index.resize(n, 0);
    beg_index.resize(n, 0);
    flags.resize(n, 0);

    let band_width_penalty = 1.0f32 / band_width_threshold;
    let mut max_result = -1i64;
    let mut max_i = -1i64;
    let mut min_x_len = x_read_len * 2 + 2;

    for i in 0..n {
        let mut n_max_skip = 0i64;
        let x_off = (x[i].offset >> 32) as i64;
        let y_off = (x[i].offset as u32) as i64;
        let mut max_j = -1i64;
        let mut max_score = x[i].weight as i64;
        let mut max_beg = i as i64;

        let mut j = i as i64 - 1;
        while j >= 0 {
            let jx = (x[j as usize].offset >> 32) as i64;
            let jy = (x[j as usize].offset as u32) as i64;
            let distance_x = x_off - jx;
            let distance_y = y_off - jy;
            // x has been sorted by x_off.
            if distance_x <= 0 || distance_y <= 0 {
                j -= 1;
                continue;
            }

            let bx = (x[beg_index[j as usize] as usize].offset >> 32) as i64;
            let by = (x[beg_index[j as usize] as usize].offset as u32) as i64;
            let total_distance_x = x_off - bx;
            let total_distance_y = y_off - by;

            let mut distance_gap = total_distance_x - total_distance_y;
            if distance_gap < 0 {
                distance_gap = -distance_gap;
            }
            // `float` comparison in the reference.
            if (distance_gap as f32) > band_width_threshold * (total_distance_x as f32) {
                j -= 1;
                continue;
            }

            let mut score = x[i].weight as i64;
            let gap_rate = (distance_gap as f32) / (total_distance_x as f32);
            score -= ((score as f32) * gap_rate * band_width_penalty) as i64;
            score += i64::from(scores[j as usize]);

            if score > max_score {
                max_score = score;
                max_j = j;
                max_beg = i64::from(beg_index[j as usize]);
                n_max_skip = 0;
            } else {
                n_max_skip += 1;
                if n_max_skip > max_skip {
                    break;
                }
            }
            j -= 1;
        }

        scores[i] = max_score as i32;
        prev_index[i] = max_j as i32;
        beg_index[i] = max_beg as i32;

        if i64::from(scores[i]) > max_result {
            max_result = i64::from(scores[i]);
            max_i = i as i64;
            let (_, _, _, _, span) =
                get_hap_overlap_len(x_off, x_off, x_read_len, y_off, y_off, y_read_len);
            min_x_len = span;
        } else if i64::from(scores[i]) == max_result {
            let (_, _, _, _, span) =
                get_hap_overlap_len(x_off, x_off, x_read_len, y_off, y_off, y_read_len);
            if span < min_x_len {
                max_result = i64::from(scores[i]);
                max_i = i as i64;
                min_x_len = span;
            }
        }
        flags[i] = 0;
    }

    while max_i != -1 {
        // `get_chain_details`: walks the prev chain into `beg_index` itself.
        let mut chain_len = 0usize;
        let mut dup = 0usize;
        let mut i = max_i;
        while i >= 0 {
            if flags[i as usize] == 1 {
                dup += 1;
            }
            beg_index[chain_len] = i as i32;
            i = i64::from(prev_index[i as usize]);
            chain_len += 1;
        }
        if chain_len == 0 {
            break;
        }
        let chain = &beg_index[..chain_len];
        if (dup as f64) > (chain_len as f64) * DUP_OVLP_RATE {
            for &b in chain {
                if flags[b as usize] == 1 {
                    continue;
                }
                flags[b as usize] = 2;
            }
        } else {
            push_hap_can(x, u_can, chain, x_read_len, y_read_len);
            for &b in chain {
                flags[b as usize] = 1;
            }
        }

        max_i = get_max_index(x, scores, flags, x_read_len, y_read_len);
    }
}

/// `get_candidate_hap_alignment` (`Purge_Dups.cpp:2847`): sort, merge equal
/// `(el, Off)` pairs, then chain the two strand groups separately.
#[allow(clippy::too_many_arguments)]
pub fn get_candidate_hap_alignment(
    u_can: &mut Vec<HapCandidate>,
    u_buffer: &mut Vec<ArcOffset>,
    scores: &mut Vec<i32>,
    prev_index: &mut Vec<i32>,
    beg_index: &mut Vec<i32>,
    flags: &mut Vec<u8>,
    band_width_threshold: f32,
    max_skip: i64,
    x_read_len: i64,
    y_read_len: i64,
) {
    u_can.clear();
    if u_buffer.is_empty() {
        return;
    }
    let mut break_point = u32::MAX;

    u_buffer.sort_by(compare_arc_offsets);

    let mut m = 1usize;
    for i in 1..u_buffer.len() {
        let mut is_merge = false;
        if u_buffer[m - 1].arc.exact_region == u_buffer[i].arc.exact_region
            && u_buffer[m - 1].offset == u_buffer[i].offset
        {
            is_merge = true;
        }
        if is_merge {
            u_buffer[m - 1].weight += u_buffer[i].weight;
            continue;
        }
        u_buffer[m] = u_buffer[i].clone();
        if u_buffer[m].arc.exact_region != u_buffer[m - 1].arc.exact_region {
            break_point = m as u32;
        }
        m += 1;
    }
    let n = m;
    u_buffer.truncate(n);
    if break_point > n as u32 {
        break_point = n as u32;
    }
    let break_point = break_point as usize;

    hap_chaining(
        &u_buffer[..break_point],
        scores,
        prev_index,
        beg_index,
        flags,
        band_width_threshold,
        max_skip,
        x_read_len,
        y_read_len,
        u_can,
    );
    hap_chaining(
        &u_buffer[break_point..],
        scores,
        prev_index,
        beg_index,
        flags,
        band_width_threshold,
        max_skip,
        x_read_len,
        y_read_len,
        u_can,
    );
}

/// `filter_secondary_chain` (`Purge_Dups.cpp:2898`).
pub fn filter_secondary_chain(max_score: i64, cur_score: i64, rate: f64) -> bool {
    if cur_score >= max_score {
        return true;
    }
    let diff = max_score - cur_score;
    let mut max_score = max_score;
    if max_score < 0 {
        max_score *= -1;
    }
    if (diff as f64) >= (max_score as f64) * (1.0 - rate) {
        return false;
    }
    true
}

/// `init_pip_hh` (`Purge_Dups.cpp:474`): the per-read reverse-coverage
/// weight; reads below the coverage floor weigh `REV_W - floor`.
pub fn init_pip_hh(
    read_graph: &Asg,
    reverse_sources: &[PafList],
    coverage_cut: &[CoverageCut],
    sc: i64,
) -> Vec<u8> {
    let mut c = vec![0u8; read_graph.sequences.len()];
    for (i, item) in c.iter_mut().enumerate() {
        if sc <= 0 {
            *item = 1;
            continue;
        }
        let r_base = i64::from(coverage_cut[i].end.wrapping_sub(coverage_cut[i].start));
        let mut c_base = 0i64;
        for hit in &reverse_sources[i].hits {
            c_base += i64::from(
                hit.query_end
                    .wrapping_sub(hit.query_start),
            );
        }
        c_base = if r_base != 0 { c_base / r_base } else { 0 };
        c_base /= sc;
        *item = 1;
        if c_base < REV_W {
            *item = (REV_W - c_base) as u8;
        }
    }
    c
}

// ---------------------------------------------------------------------------
// The per-unitig scan driver (`hap_alignment_advance_worker`, Purge_Dups.cpp
// 2964–3216).
// ---------------------------------------------------------------------------

/// The per-unitig scratch state; the C version is one `hap_alignment_struct`
/// per thread, the port runs single-threaded.
pub struct HapAlignmentScratch {
    /// `vote_counting`.
    pub votes: Vec<u64>,
    /// `visit`.
    pub visit: Vec<u8>,
    /// `u_vecs`: the sorted `(votes << 32) | unitig` candidate list.
    pub u_vecs: Vec<u64>,
    /// `u_buffer`.
    pub u_buffer: Vec<ArcOffset>,
    /// `u_buffer_tailIndex`, doubles as the LIS tail inside the candidate
    /// validation.
    pub tail_index: Vec<i32>,
    /// `u_buffer_prevIndex`.
    pub prev_index: Vec<i32>,
    /// `u_buffer_beg`.
    pub beg_index: Vec<i32>,
    /// `u_buffer_flag`.
    pub flags: Vec<u8>,
    /// `u_can`; the vector is reused across candidate batches like the
    /// reference's `kvec`.
    pub candidates: Vec<HapCandidate>,
}

impl HapAlignmentScratch {
    /// `init_hap_alignment_struct`.
    pub fn new(unitig_count: usize) -> Self {
        Self {
            votes: vec![0; unitig_count],
            visit: vec![0; unitig_count],
            u_vecs: Vec::new(),
            u_buffer: Vec::new(),
            tail_index: Vec::new(),
            prev_index: Vec::new(),
            beg_index: Vec::new(),
            flags: Vec::new(),
            candidates: Vec::new(),
        }
    }
}

/// The vote half plus the candidate scan of one unitig
/// (`hap_alignment_advance_worker`).
///
/// `hap_rate` is the reference's `Hap_rate` (`purge_simi_thres`),
/// `chain_rate` the fixed `0.1` of `init_hap_alignment_struct_pip`.
#[allow(clippy::too_many_arguments)]
pub fn hap_alignment_scan_worker(
    unitig_graph: &UnitigGraph,
    ctx: &HapAlignmentContext,
    x_uid: u32,
    hh: &[u8],
    hap_rate: f32,
    chain_rate: f32,
    cov_threshold: u64,
    cov: &HapCov,
    all_ovlp: &mut crate::assembly::purge_dups::HapOverlapsList,
    scratch: &mut HapAlignmentScratch,
    clean_visit: &dyn Fn(&mut [u8], u32, &[crate::overlap::paf::PafHit]),
) {
    let nsg = &unitig_graph.graph;
    if nsg.sequences[x_uid as usize].deleted
        || nsg.sequences[x_uid as usize].coverage_class == crate::assembly::hit::ALTER_LABEL
    {
        return;
    }
    let contig_count = nsg.sequences.len();

    scratch.votes.fill(0);
    scratch.visit.fill(0);

    // The vote half: every read of x_uid votes once per candidate unitig.
    let x_reads = &unitig_graph.unitigs[x_uid as usize];
    for &packed in &x_reads.reads {
        let x_reverse = &ctx.reverse_sources[(packed >> 33) as usize];
        for hit in &x_reverse.hits {
            let mut r_id = hit.target_identifier;
            if ctx.read_graph.sequences[r_id as usize].deleted {
                let (mapped, is_unitig) = ctx.read_to_unitig.get(r_id);
                r_id = mapped;
                if r_id == u32::MAX
                    || is_unitig == 1
                    || ctx.read_graph.sequences[r_id as usize].deleted
                {
                    continue;
                }
            }

            let (hap_uid, is_unitig) = ctx.read_to_unitig.get(r_id);
            if is_unitig == 0 || hap_uid == u32::MAX {
                continue;
            }
            if scratch.visit[hap_uid as usize] != 0 {
                continue;
            }
            scratch.visit[hap_uid as usize] = 1;
            if scratch.votes[hap_uid as usize] < u64::MAX {
                scratch.votes[hap_uid as usize] += 1;
            }
        }

        clean_visit(&mut scratch.visit, contig_count as u32, &x_reverse.hits);
    }

    // Collect the candidates and sort them descending (radix + reverse).
    scratch.u_vecs.clear();
    for i in 0..contig_count {
        if i as u32 == x_uid || scratch.votes[i] == 0 {
            continue;
        }
        scratch.u_vecs.push((scratch.votes[i] << 32) | i as u64);
    }
    if scratch.u_vecs.is_empty() {
        return;
    }
    scratch.u_vecs.sort_unstable();
    scratch.u_vecs.reverse();

    for &entry in &scratch.u_vecs {
        let y_uid = entry as u32;
        let seed_occ = (entry >> 32) as u32;
        let y_reads = &unitig_graph.unitigs[y_uid as usize];
        scratch.u_buffer.clear();

        for k in 0..x_reads.reads.len() {
            let x_read = (x_reads.reads[k] >> 33) as usize;
            let hhc = hh[x_read];
            for hit in &ctx.reverse_sources[x_read].hits {
                let Some(mut arc) = scan_arc(ctx, hit, y_uid) else {
                    continue;
                };

                let v = (x_reads.reads[k] >> 32) as u32;
                let (v_unitig, is_unitig) = ctx.read_to_unitig.get(v >> 1);
                if is_unitig == 0 || v_unitig == u32::MAX || v_unitig != x_uid {
                    continue;
                }
                if ctx.position_index[(v >> 1) as usize] as u32 != k as u32 {
                    continue;
                }

                let w =
                    (y_reads.reads[(ctx.position_index[(arc.v >> 1) as usize] as u32) as usize] >> 32) as u32;

                if ctx.purge_level <= 2 {
                    let drop = prefilter(
                        ctx.position_index[(v >> 1) as usize] as u32,
                        ctx.position_index[(arc.v >> 1) as usize] as u32,
                        x_reads.reads.len() as u32,
                        y_reads.reads.len() as u32,
                        0,
                        hap_rate,
                        seed_occ,
                    ) == NON_PLOID
                        && prefilter(
                            ctx.position_index[(v >> 1) as usize] as u32,
                            ctx.position_index[(arc.v >> 1) as usize] as u32,
                            x_reads.reads.len() as u32,
                            y_reads.reads.len() as u32,
                            1,
                            hap_rate,
                            seed_occ,
                        ) == NON_PLOID;
                    if drop {
                        continue;
                    }
                }

                let (packed_xy, arc_rev) = get_xy_pos(
                    ctx.read_graph,
                    &arc,
                    v,
                    w,
                    x_reads.length,
                    y_reads.length,
                    ctx.position_index,
                );
                if (packed_xy >> 32) == u32::MAX as u64 || (packed_xy as u32) == u32::MAX {
                    continue;
                }
                arc.exact_region = arc_rev;

                scratch.u_buffer.push(ArcOffset {
                    offset: packed_xy,
                    arc,
                    weight: u64::from(hhc),
                });
            }

            deduplicate_edge(&mut scratch.u_buffer);
        }

        if scratch.u_buffer.is_empty() {
            continue;
        }

        get_candidate_hap_alignment(
            &mut scratch.candidates,
            &mut scratch.u_buffer,
            &mut scratch.tail_index,
            &mut scratch.prev_index,
            &mut scratch.beg_index,
            &mut scratch.flags,
            chain_rate,
            50,
            i64::from(x_reads.length),
            i64::from(y_reads.length),
        );

        if scratch.candidates.is_empty() {
            continue;
        }

        scratch.candidates.sort_by(compare_hap_candidates);

        let m = all_ovlp.x[x_uid as usize].len();
        let mut max_weight = 0u64;
        let mut max_score = 0i64;
        for can_index in 0..scratch.candidates.len() {
            let can = scratch.candidates[can_index].clone();
            if (can.weight as f64) < (max_weight as f64) * 0.33 {
                continue;
            }
            let (ploid, r_x_beg, r_x_end, r_y_beg, r_y_end, sim) =
                calculate_pair_hap_similarity_advance(
                    &mut scratch.candidates[can_index],
                    x_uid,
                    y_uid,
                    x_reads,
                    y_reads,
                    ctx,
                    hap_rate,
                    ctx.purge_level > 2,
                    cov_threshold,
                    &mut scratch.u_buffer,
                    &mut scratch.tail_index,
                    &mut scratch.prev_index,
                    cov,
                );
            if !ploid {
                continue;
            }
            // The candidate is consumed in place; take the updated fields.
            let can = &scratch.candidates[can_index];
            if max_weight == 0 || max_score < can.score {
                max_score = can.score;
            }
            if max_weight < can.weight {
                max_weight = can.weight;
            }
            if !filter_secondary_chain(max_score, can.score, CHAIN_FILTER_RATE) {
                continue;
            }

            let mut hap_align = crate::assembly::purge_dups::HapOverlap::default();
            hap_align.reverse = can.rev;
            hap_align.kind = can.index_end as u8;
            hap_align.x_beg_id = can.x_beg_pos;
            hap_align.x_end_id = can.x_end_pos.wrapping_add(1);
            hap_align.y_beg_id = can.y_beg_pos;
            hap_align.y_end_id = can.y_end_pos.wrapping_add(1);
            hap_align.weight = can.weight as u32;
            hap_align.score = can.score;
            hap_align.x_beg_pos = r_x_beg as u32;
            hap_align.x_end_pos = (r_x_end + 1) as u32;
            if hap_align.reverse == 0 {
                hap_align.y_beg_pos = r_y_beg as u32;
                hap_align.y_end_pos = (r_y_end + 1) as u32;
            } else {
                hap_align.y_beg_pos = y_reads
                    .length
                    .wrapping_sub(r_y_end as u32)
                    .wrapping_sub(1);
                hap_align.y_end_pos = y_reads
                    .length
                    .wrapping_sub(r_y_beg as u32)
                    .wrapping_sub(1)
                    .wrapping_add(1);
            }
            hap_align.x_uid = x_uid;
            hap_align.y_uid = y_uid;
            hap_align.status = crate::assembly::purge_dups::SELF_EXIST;
            hap_align.similarity = sim;
            all_ovlp.x[hap_align.x_uid as usize].push(hap_align);
        }

        // Keep only the best chain with the maximal score, compacted to `m`.
        let list = &all_ovlp.x[x_uid as usize];
        let mut hap_align_x: Option<(usize, crate::assembly::purge_dups::HapOverlap)> = None;
        for k in m..list.len() {
            let entry = &list[k];
            if entry.score != max_score {
                continue;
            }
            match &hap_align_x {
                None => hap_align_x = Some((k, entry.clone())),
                Some((_, best)) => {
                    if entry.weight > best.weight {
                        hap_align_x = Some((k, entry.clone()));
                    } else if entry.weight == best.weight {
                        let span = entry.x_end_pos.wrapping_sub(entry.x_beg_pos);
                        let best_span = best.x_end_pos.wrapping_sub(best.x_beg_pos);
                        if span < best_span {
                            hap_align_x = Some((k, entry.clone()));
                        }
                    }
                }
            }
        }
        if let Some((_, best)) = hap_align_x {
            let list = &mut all_ovlp.x[x_uid as usize];
            if m < list.len() {
                list[m] = best;
            } else {
                list.push(best);
            }
            list.truncate(m + 1);
        }
    }
}

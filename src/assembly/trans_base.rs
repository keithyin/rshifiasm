//! `trans_base_infer` (`inter.cpp` 22574): base-level transitions inferred
//! by mapping every unitig against the unitig graph itself with a global
//! minimizer index.
//!
//! The index is built in two sketch passes (`gen_trans_base_count_comp`,
//! `inter.cpp` 22181): the first counts the `mz2` minimizers of every
//! non-homozygous unitig, the second fills the flat array. The array is then
//! sorted by hash, every entry's `x` is remapped to its sorted position, and
//! the occupancy quantile tightens `mini_cut`. The sorted array doubles as
//! the global lookup index of `minimizers_qgen_input`.

use crate::assembly::bubble_chain::BubbleType;
use crate::assembly::unitig::UnitigGraph;
use crate::command_line::options::AssemblyOptions;
use crate::container::ksort::radix_sort_by_key;
use crate::counting::minimizer::sketch;

/// One global minimizer, `ha_mzl_t`: `x` is the hash key before the remap
/// and the sorted position afterwards; `rid`/`pos`/`rev`/`span` locate the
/// k-mer on its unitig.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MzlEntry {
    pub x: u64,
    pub rid: u32,
    pub pos: u32,
    pub reverse: bool,
    pub span: u8,
}

/// The global minimizer index produced by the two count passes.
pub struct TransBaseIndex {
    /// `idx_a`: per-unitig minimizers, `x` remapped to the sorted position.
    pub entries: Vec<MzlEntry>,
    /// `idx_n`: prefix offsets, `entries[offsets[i]..offsets[i+1]]` belongs
    /// to unitig `i`.
    pub offsets: Vec<u64>,
    /// `srt_a`: the hash-sorted global index (original hash in `x`).
    pub sorted: Vec<MzlEntry>,
    /// `mini_cut` after the occupancy-quantile adjustment.
    pub mini_cut: u64,
}

/// `mz2_ha_sketch` through the shared `HA_SC_INIT` selection: the `mz2`
/// instantiation differs from `mz1` only in the field widths, and this call
/// path passes no filter table and no position index, so the `mz1` port
/// reproduces it exactly.
fn sketch_unitig(
    unitig_graph: &UnitigGraph,
    id: usize,
    options: &AssemblyOptions,
) -> Vec<MzlEntry> {
    let result = sketch(
        &unitig_graph.unitigs[id].sequence,
        options.transition_window_size,
        options.transition_k_mer_length,
        id as u32,
        false,
        None,
        options.minimizer_sample_distance,
        options.minimizer_rewindow_size,
    );
    result
        .minimizers
        .iter()
        .map(|m| MzlEntry {
            x: m.hash,
            rid: id as u32,
            pos: m.position,
            reverse: m.reverse,
            span: m.span as u8,
        })
        .collect()
}

/// Whether `worker_for_trans_ovlp_adv` skips the unitig: homozygous bubble
/// nodes and deleted sequences take no part in the base transitions.
fn skipped(unitig_graph: &UnitigGraph, bubble: &BubbleType, id: usize) -> bool {
    bubble.is_hom(id as u32) || unitig_graph.graph.sequences[id].deleted
}

/// The count/fill/sort/remap half of `gen_trans_base_count_comp`
/// (`inter.cpp` 22181): builds the global index and adjusts `mini_cut`.
pub fn build_trans_base_index(
    unitig_graph: &UnitigGraph,
    bubble: &BubbleType,
    options: &AssemblyOptions,
) -> TransBaseIndex {
    let n = unitig_graph.unitigs.len();
    let mut mini_cut = 512u64;

    // First pass: count the minimizers of every active unitig.
    let mut offsets = vec![0u64; n + 1];
    let mut total = 0u64;
    for i in 0..n {
        if !skipped(unitig_graph, bubble, i) {
            total += sketch_unitig(unitig_graph, i, options).len() as u64;
        }
        offsets[i + 1] = total;
    }

    // Second pass: fill the flat array in unitig order.
    let mut entries: Vec<MzlEntry> = Vec::with_capacity(total as usize);
    for i in 0..n {
        if skipped(unitig_graph, bubble, i) {
            continue;
        }
        entries.extend(sketch_unitig(unitig_graph, i, options));
    }
    assert_eq!(entries.len() as u64, total);

    // srt_a keeps the original index packed into rid/pos while sorting.
    let mut sorted: Vec<MzlEntry> = entries
        .iter()
        .enumerate()
        .map(|(i, e)| MzlEntry {
            x: e.x,
            rid: (i as u64 >> 32) as u32,
            pos: i as u32,
            reverse: e.reverse,
            span: e.span,
        })
        .collect();
    radix_sort_by_key(&mut sorted, &|e: &MzlEntry| e.x);

    // Remap every entry to its sorted position and collect group sizes.
    let mut cut: Vec<u64> = Vec::new();
    let mut l = 0usize;
    for k in 1..=sorted.len() {
        if k == sorted.len() || sorted[l].x != sorted[k].x {
            for i in l..k {
                let m = ((sorted[i].rid as u64) << 32) | sorted[i].pos as u64;
                assert_eq!(sorted[i].x, entries[m as usize].x);
                sorted[i] = entries[m as usize];
                entries[m as usize].x = i as u64;
            }
            cut.push((k - l) as u64);
            l = k;
        }
    }

    if !cut.is_empty() {
        radix_sort_by_key(&mut cut, &|&c: &u64| c);
        let m = (cut.len() as f64 * 0.0002) as usize;
        let mut cc = cut[cut.len() - 1] + 1;
        if m > 0 && m <= cut.len() {
            cc = cut[cut.len() - m] + 1;
        }
        if cc < mini_cut {
            mini_cut = cc;
        }
    }

    if let Some(directory) = crate::assembly::probe::dump_directory() {
        let mut out = format!("N {} {}\n", sorted.len(), mini_cut);
        for (i, e) in sorted.iter().enumerate() {
            out.push_str(&format!(
                "S {} {} {} {} {} {}\n",
                i, e.x, e.rid, e.pos, u8::from(e.reverse), e.span
            ));
        }
        for (i, &off) in offsets.iter().enumerate() {
            out.push_str(&format!("O {} {}\n", i, off));
        }
        for (i, e) in entries.iter().enumerate() {
            out.push_str(&format!(
                "I {} {} {} {} {} {}\n",
                i, e.x, e.rid, e.pos, u8::from(e.reverse), e.span
            ));
        }
        let _ = std::fs::write(directory.join("trans_base_idx.txt"), out);
    }

    TransBaseIndex {
        entries,
        offsets,
        sorted,
        mini_cut,
    }
}
// ---------------------------------------------------------------------------
// The `trans_base_infer` driver (`inter.cpp` 22574): maps every active unitig
// against the polished graph through the global minimizer index, refines the
// chains at base level and appends the `RC_3` transitions.
//
// The overlap pass (`worker_for_trans_ovlp_adv`, `inter.cpp` 10383) runs
// `gen_trans_adaptive_aln` (10297) per query: chain through
// `minimizers_qgen_input` / `lchain_qgen_mcopy_input`, deduplicate against
// the previous round's transitions, filter by the reliable ones, align in two
// error-rate slices (`split_ug_lalign`), and push the surviving overlaps and
// reliable transitions into the query's track list. The tail of
// `gen_trans_base_count_comp` (22181) then turns the track lists into
// `RC_3` transitions.

use crate::assembly::asg::Asg;
use crate::assembly::hap_cov::{CHAIN_MATCH, CHAIN_UNMATCH, RC_0, RC_1, RC_3, UTrans};
use crate::assembly::trans_clean::{clean_u_trans_t_idx_adv, u_trans_span};
use crate::assembly::ul::UlOv;
use crate::correction::hc_align::BitExtz;
use crate::correction::target::TargetSource;
use crate::correction::ug_lalign::split_ug_lalign;
use crate::overlap::lchain_parameters::set_lchain_parameters;
use crate::overlap::qgen_input::{
    InputChainParameters, SeedInput, lchain_qgen_mcopy_input, minimizers_qgen_input,
};
use crate::overlap::region::{CandidatesList, ChainData, OverlapRegion};
use crate::reads::read_store::ReadStore;

/// `trans_ovlp_connect` (`inter.cpp` 9537): whether a candidate transition's
/// overlap lengths are close enough to a graph arc between the two vertices.
pub(crate) fn trans_ovlp_connect(p: &UlOv, unitig_graph: &UnitigGraph) -> bool {
    let mut v = p.qn << 1;
    let mut w = (p.tn << 1) + u32::from(p.rev);

    for _pass in 0..2 {
        let mut arc_overlap = None;
        let arcs = unitig_graph.graph.arcs_of(v);
        for arc in arcs {
            if arc.deleted || arc.v != w {
                continue;
            }
            arc_overlap = Some(arc.overlap);
            break;
        }
        if let Some(dg) = arc_overlap {
            if interval_close(dg, p.qe - p.qs) || interval_close(dg, p.te - p.ts) {
                return false;
            }
        }
        v ^= 1;
        w ^= 1;
    }
    true
}

/// The `dif <= mm` check of `trans_ovlp_connect`, `mm = MAX(dg, da) * 0.06`
/// truncated to `uint32` and floored at `8`.
fn interval_close(dg: u32, da: u32) -> bool {
    let dif = if dg > da { dg - da } else { da - dg };
    let mut mm = ((dg.max(da) as f64) * 0.06) as u32;
    if mm < 8 {
        mm = 8;
    }
    dif <= mm
}

/// `filter_by_reliable_ovlp0` (`inter.cpp` 9578): whether `[s, e)` overlaps a
/// reliable (`RC_0`/`RC_1`) transition by more than `sec_rate`.
pub(crate) fn filter_by_reliable_ovlp0(
    transitions: &[UTrans],
    s: u32,
    e: u32,
    sec_rate: f64,
) -> bool {
    for a in transitions {
        if a.deleted {
            continue;
        }
        if a.flag == RC_0 || a.flag == RC_1 {
            let lower = a.query_start.max(s);
            let upper = a.query_end.min(e);
            let ovlp = if upper > lower { upper - lower } else { 0 };
            if ovlp == 0 {
                continue;
            }
            if ovlp as f64 > ((e - s) as f64 * sec_rate) {
                return true;
            }
        }
    }
    false
}

/// `filter_by_reliable_ovlp_adv` (`inter.cpp` 9799): deduplicates the chained
/// overlaps against the reliable transitions, drops the ones covered by
/// them, rotates the marked overlaps to the front and returns their count
/// (the `err_high` slice `ol_h`).
///
/// `overlaps` holds the reference's `ol->list` and `length` its `ol->length`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn filter_by_reliable_ovlp_adv(
    id: u32,
    transitions: &[UTrans],
    index: &[u64],
    scratch: &mut Vec<u64>,
    overlaps: &mut Vec<OverlapRegion>,
    length: &mut usize,
    target_lengths: &[u32],
    sec_rate: f64,
    avoid_dup_aln: bool,
    dedup_by_reliable_ovlp: bool,
) -> usize {
    let (off, n) = crate::assembly::trans_clean::u_trans_span(index, id);
    let a = &transitions[off..off + n];
    let mut marked = 0usize;

    if avoid_dup_aln {
        scratch.clear();
        for (k, t) in a.iter().enumerate() {
            if t.deleted {
                continue;
            }
            if t.flag == RC_0 || t.flag == RC_1 {
                let z = ((u64::from(t.target_identifier) << 1) | u64::from(t.reverse)) << 32;
                scratch.push(z);
            }
        }
        if !scratch.is_empty() {
            for (k, entry) in overlaps.iter().enumerate().take(*length) {
                let z = (((u64::from(entry.y_identifier) << 1)
                    | u64::from(entry.y_strand != 0))
                    << 32)
                    | k as u64
                    | 0x8000_0000;
                scratch.push(z);
            }

            crate::container::ksort::radix_sort_u64(scratch);
            let mut duplicates = 0usize;
            let mut l = 0usize;
            for k in 1..=scratch.len() {
                if k == scratch.len() || (scratch[l] >> 32) != (scratch[k] >> 32) {
                    if k - l > 1 && (scratch[l] & 0x8000_0000) == 0 {
                        for z in l..k {
                            if scratch[z] & 0x8000_0000 != 0 {
                                let slot = (scratch[z] & 0x7fff_ffff) as usize;
                                overlaps[slot].y_identifier = u32::MAX;
                                duplicates += 1;
                            }
                        }
                    }
                    l = k;
                }
            }

            if duplicates > 0 {
                let mut write = 0usize;
                for k in 0..*length {
                    if overlaps[k].y_identifier == u32::MAX {
                        continue;
                    }
                    if write != k {
                        overlaps.swap(write, k);
                    }
                    write += 1;
                }
                *length = write;
            }
        }
    }

    // The merged reliable query intervals of `id`.
    let mut reliable: Vec<u64> = Vec::new();
    if dedup_by_reliable_ovlp {
        for t in a {
            if t.deleted {
                continue;
            }
            if t.flag == RC_0 || t.flag == RC_1 {
                reliable.push((u64::from(t.query_start) << 32) | u64::from(t.query_end));
            }
        }
        if reliable.len() > 1 {
            crate::container::ksort::radix_sort_u64(&mut reliable);
            let mut z = 0usize;
            for k in 0..reliable.len() {
                let s = (reliable[k] >> 32) as u32;
                let e = reliable[k] as u32;
                if z > 0 && s <= (reliable[z - 1] as u32) {
                    let previous_end = reliable[z - 1] as u32;
                    if e > previous_end {
                        reliable[z - 1] = (reliable[z - 1] >> 32 << 32) | u64::from(e);
                    }
                } else {
                    reliable[z] = reliable[k];
                    z += 1;
                }
            }
            reliable.truncate(z);
        }
    }
    let spn = reliable.len();

    let mut write = 0usize;
    let mut rotated = 0usize;
    for k in 0..*length {
        let entry = &overlaps[k];
        let mut drop = false;
        if dedup_by_reliable_ovlp && entry.x_strand == 0 {
            let s = entry.x_position_start;
            let e = entry.x_position_end + 1;
            let mut covered = 0u64;
            let mut z = 0usize;
            while z < spn {
                let s0 = (reliable[z] >> 32) as u32;
                let e0 = reliable[z] as u32;
                let lower = s.max(s0);
                let upper = e.min(e0);
                let ov = if upper > lower { upper - lower } else { 0 };
                covered += u64::from(ov);
                if covered > 0 && covered as f64 > ((e - s) as f64 * sec_rate) {
                    break;
                }
                z += 1;
            }
            // Covered by the reliable transitions of `id`.
            if z < spn {
                drop = true;
            }
            if !drop {
                let target_length = target_lengths[entry.y_identifier as usize];
                let (ys, ye) = if entry.y_strand != 0 {
                    (
                        target_length - (entry.y_position_end + 1),
                        target_length - entry.y_position_start,
                    )
                } else {
                    (entry.y_position_start, entry.y_position_end + 1)
                };
                let (off, n) = crate::assembly::trans_clean::u_trans_span(index, entry.y_identifier);
                if filter_by_reliable_ovlp0(&transitions[off..off + n], ys, ye, sec_rate) {
                    drop = true;
                }
            }
        }
        if drop {
            continue;
        }

        if write != k {
            overlaps.swap(write, k);
        }
        if overlaps[write].x_strand != 0 {
            overlaps[write].x_strand = 0;
            if rotated != write {
                overlaps.swap(rotated, write);
            }
            rotated += 1;
        }
        write += 1;
    }
    *length = write;
    rotated
}

/// `cal_trans_ov_w` (`inter.cpp` 12562): the chain score of one overlap,
/// `min` of both sides.
pub(crate) fn cal_trans_ov_w(z: &UlOv) -> f64 {
    let uml = f64::from(z.sec);
    let mut ml = f64::from(z.qe - z.qs);
    ml -= uml;
    let forward = ml * CHAIN_MATCH - uml * CHAIN_UNMATCH;
    let mut ml = f64::from(z.te - z.ts);
    ml -= uml;
    let reverse = ml * CHAIN_MATCH - uml * CHAIN_UNMATCH;
    forward.min(reverse)
}

/// `test_het_aln` (`inter.cpp` 10037): whether the reliable transitions of
/// `rid` plus the chained overlaps tile at least 80% of the unitig.
pub(crate) fn test_het_aln(
    transitions: &[UTrans],
    scratch: &mut Vec<u64>,
    overlaps: &[OverlapRegion],
    length: usize,
    len: u32,
) -> bool {
    scratch.clear();
    for t in transitions {
        if t.deleted {
            continue;
        }
        scratch.push((u64::from(t.query_start) << 32) | u64::from(t.query_end));
    }
    for entry in overlaps.iter().take(length) {
        let s = entry.x_position_start;
        let e = entry.x_position_end + 1;
        scratch.push((u64::from(s) << 32) | u64::from(e));
    }

    if scratch.len() > 1 {
        crate::container::ksort::radix_sort_u64(scratch);
        let mut z = 0usize;
        for k in 0..scratch.len() {
            let s = (scratch[k] >> 32) as u32;
            let e = scratch[k] as u32;
            if z > 0 && s <= (scratch[z - 1] as u32) {
                let previous_end = scratch[z - 1] as u32;
                if e > previous_end {
                    scratch[z - 1] = (scratch[z - 1] >> 32 << 32) | u64::from(e);
                }
            } else {
                scratch[z] = scratch[k];
                z += 1;
            }
        }
        scratch.truncate(z);
    }
    let mut covered = 0u64;
    for &interval in scratch.iter() {
        covered += u64::from(interval as u32) - (interval >> 32) as u64;
    }

    covered > 0 && covered as f64 >= len as f64 * 0.8
}

/// `backward_dedup_ol` (`inter.cpp` 10237): drops the chained overlaps that
/// duplicate a track entry of the same query from an earlier round.
pub(crate) fn backward_dedup_ol(
    rid: u32,
    track: &mut [UlOv],
    track_length: usize,
    scratch: &mut Vec<u64>,
    overlaps: &mut Vec<OverlapRegion>,
    length: &mut usize,
) {
    scratch.clear();
    let mut run_start = track_length;
    while run_start > 0 && track[run_start - 1].qn == rid {
        run_start -= 1;
        let z = ((u64::from(track[run_start].tn) << 1) | u64::from(track[run_start].rev)) << 32;
        scratch.push(z);
    }
    if scratch.is_empty() {
        return; // no track overlap
    }

    for (k, entry) in overlaps.iter().enumerate().take(*length) {
        let z = (((u64::from(entry.y_identifier) << 1) | u64::from(entry.y_strand != 0)) << 32)
            | k as u64
            | 0x8000_0000;
        scratch.push(z);
    }

    crate::container::ksort::radix_sort_u64(scratch);
    let mut duplicates = 0usize;
    let mut l = 0usize;
    for k in 1..=scratch.len() {
        if k == scratch.len() || (scratch[l] >> 32) != (scratch[k] >> 32) {
            if k - l > 1 && (scratch[l] & 0x8000_0000) == 0 {
                for z in l..k {
                    if scratch[z] & 0x8000_0000 != 0 {
                        let slot = (scratch[z] & 0x7fff_ffff) as usize;
                        overlaps[slot].y_identifier = u32::MAX;
                        duplicates += 1;
                    }
                }
            }
            l = k;
        }
    }
    if duplicates > 0 {
        let mut write = 0usize;
        for k in 0..*length {
            if overlaps[k].y_identifier == u32::MAX {
                continue;
            }
            if write != k {
                overlaps.swap(write, k);
            }
            write += 1;
        }
        *length = write;
    }
}

/// `remove_trans_ovlp_connect` (`inter.cpp` 10284): filters the query's run
/// of track entries down to the ones a graph arc supports.
pub(crate) fn remove_trans_ovlp_connect(
    unitig_graph: &UnitigGraph,
    rid: u32,
    track: &mut Vec<UlOv>,
) {
    let mut run_start = track.len();
    while run_start > 0 && track[run_start - 1].qn == rid {
        run_start -= 1;
    }
    if run_start >= track.len() {
        return;
    }
    let mut write = run_start;
    for k in run_start..track.len() {
        if !trans_ovlp_connect(&track[k], unitig_graph) {
            continue;
        }
        track.swap(write, k);
        write += 1;
    }
    track.truncate(write);
}

/// `push_ul_ov_t` (`inter.cpp` 10144): appends the chained overlaps and the
/// reliable transitions the chains did not reproduce to the query's track
/// list.
///
/// `transitions` is `rid`'s slice of the filter table (`a`), empty on the
/// second round; `is_arc_filter` is the `test_het_aln` outcome of the first
/// round. Like the reference, nothing is appended when the overlap list is
/// empty.
#[allow(clippy::too_many_arguments)]
pub(crate) fn push_ul_ov_t(
    transitions: &[UTrans],
    scratch: &mut Vec<u64>,
    overlaps: &[OverlapRegion],
    length: usize,
    target_lengths: &[u32],
    unitig_graph: &UnitigGraph,
    is_arc_filter: bool,
    max_err: f64,
    track: &mut Vec<UlOv>,
) {
    scratch.clear();
    let mut reliable_transitions: Vec<usize> = Vec::new();
    if !transitions.is_empty() {
        for (k, t) in transitions.iter().enumerate() {
            if t.deleted {
                continue;
            }
            let z = ((u64::from(t.target_identifier) << 1) | u64::from(t.reverse)) << 32
                | k as u64
                | 0x8000_0000;
            scratch.push(z);
        }
        for (k, entry) in overlaps.iter().enumerate().take(length) {
            let z = ((u64::from(entry.y_identifier) << 1) | u64::from(entry.y_strand != 0)) << 32
                | k as u64;
            scratch.push(z);
        }

        crate::container::ksort::radix_sort_u64(scratch);
        let mut l = 0usize;
        for k in 1..=scratch.len() {
            if k == scratch.len() || (scratch[l] >> 32) != (scratch[k] >> 32) {
                // Only the groups with no chained overlap survive; the
                // transition-only entries lose their high bit.
                if scratch[l] & 0x8000_0000 != 0 {
                    for z in l..k {
                        reliable_transitions
                            .push((scratch[z] & 0x7fff_ffff) as usize);
                    }
                }
                l = k;
            }
        }
    }

    if length > 0 {
        for entry in overlaps.iter().take(length) {
            let target_length = target_lengths[entry.y_identifier as usize];
            let (ts, te) = if entry.y_strand != 0 {
                (
                    target_length - (entry.y_position_end + 1),
                    target_length - entry.y_position_start,
                )
            } else {
                (entry.y_position_start, entry.y_position_end + 1)
            };
            let pushed = UlOv {
                qn: entry.x_identifier,
                tn: entry.y_identifier,
                qs: entry.x_position_start,
                qe: entry.x_position_end + 1,
                ts,
                te,
                rev: entry.y_strand != 0,
                el: true,
                sec: entry.non_homopolymer_errors,
            };
            // The reference pushes first and pops on failure.
            if is_arc_filter && !trans_ovlp_connect(&pushed, unitig_graph) {
                continue;
            }
            track.push(pushed);
        }

        for k in reliable_transitions {
            let t = &transitions[k];
            let pushed = UlOv {
                qn: t.query_identifier,
                tn: t.target_identifier,
                qs: t.query_start,
                qe: t.query_end,
                ts: t.target_start,
                te: t.target_end,
                rev: t.reverse,
                el: false,
                sec: (((t.query_end - t.query_start) as f64) * max_err) as u32,
            };
            if is_arc_filter && !trans_ovlp_connect(&pushed, unitig_graph) {
                continue;
            }
            track.push(pushed);
        }
    }
}

/// `gen_trans_adaptive_aln` (`inter.cpp` 10297): one query of the overlap
/// pass — chain the unitig against the global index, drop the duplicates of
/// the earlier rounds, keep the reliable transitions, align the chains in two
/// error-rate slices and append the survivors to the query's track list.
///
/// The first round returns the `test_het_aln` outcome; the second round (run
/// with `first_round == false`) always returns `false` and filters the track
/// through the graph arcs.
#[allow(clippy::too_many_arguments)]
fn gen_trans_adaptive_aln(
    rid: u32,
    query: &[u8],
    query_length: u64,
    query_minimizers: &[MzlEntry],
    index: &TransBaseIndex,
    target_lengths: &[u32],
    unitig_graph: &UnitigGraph,
    transitions: &[UTrans],
    transition_index: &[u64],
    targets: &TargetSource,
    chain_parameters: &InputChainParameters,
    first_round: bool,
    error_low: f64,
    error_high: f64,
    high_occurrence: u32,
    sec_cutoff: f64,
    k_mer_length: i64,
    chain_cut: u64,
    track: &mut Vec<UlOv>,
    scratch: &mut Vec<u64>,
    extz: &mut BitExtz,
) -> bool {
    let (off, n) = u_trans_span(transition_index, rid);
    let known = &transitions[off..off + n];

    // `ug_map_lchain`: the query minimizers against the hash-sorted index.
    let seed = SeedInput {
        query_identifier: rid,
        query_minimizers,
        index: &index.sorted,
        minimizer_cutoff: index.mini_cut,
        chain_cutoff: chain_cut,
        known_transitions: known,
        target_lengths,
        high_occurrence: Some(high_occurrence),
        low_occurrence: None,
    };
    let anchors = minimizers_qgen_input(&seed);
    let (mut overlaps, candidates) = lchain_qgen_mcopy_input(
        anchors.candidates,
        rid,
        query_length,
        target_lengths,
        &anchors.group_starts,
        chain_parameters,
    );
    let mut candidates = CandidatesList {
        list: candidates,
        length: 0,
        size: 0,
        chain_data: ChainData::default(),
    };
    candidates.length = candidates.list.len() as i64;
    candidates.size = candidates.length;
    let mut length = overlaps.len();
    // `if(!fi) backward_dedup_ol(rid, bl, &(b->sp), &b->olist);`
    if !first_round {
        let track_length = track.len();
        backward_dedup_ol(
            rid,
            &mut track[..],
            track_length,
            scratch,
            &mut overlaps,
            &mut length,
        );
    }

    // `filter_by_reliable_ovlp_adv(..., sec_cutoff, 1, 1, &ol_h)`
    let mut head = filter_by_reliable_ovlp_adv(
        rid,
        transitions,
        transition_index,
        scratch,
        &mut overlaps,
        &mut length,
        target_lengths,
        sec_cutoff,
        true,
        true,
    );
    if !first_round {
        head = 0;
    }

    // The first slice aligns at `errh`, the second at `errl`; the auxiliary
    // overlap region is only generated between the two passes.
    head = split_ug_lalign(
        head,
        &mut overlaps,
        length,
        error_high,
        error_low,
        &mut candidates,
        targets,
        query,
        query_length as i64,
        extz,
        None,
        k_mer_length,
        chain_cut,
    );
    let mut auxiliary = OverlapRegion::default();
    split_ug_lalign(
        head,
        &mut overlaps,
        length,
        error_high,
        error_low,
        &mut candidates,
        targets,
        query,
        query_length as i64,
        extz,
        Some(&mut auxiliary),
        k_mer_length,
        chain_cut,
    );

    if first_round {
        let pass_aln = test_het_aln(known, scratch, &overlaps, length, query_length as u32);
        push_ul_ov_t(
            known,
            scratch,
            &overlaps,
            length,
            target_lengths,
            unitig_graph,
            pass_aln,
            error_high,
            track,
        );
        pass_aln
    } else {
        push_ul_ov_t(
            &[],
            scratch,
            &overlaps,
            length,
            target_lengths,
            unitig_graph,
            false,
            error_high,
            track,
        );
        remove_trans_ovlp_connect(unitig_graph, rid, track);
        false
    }
}

/// The tail of `gen_trans_base_count_comp` (`inter.cpp` 22234): turns the
/// per-query track lists into `RC_3` transitions.
///
/// With one thread the reference walks the single track list in query order,
/// so grouping the contiguous runs by `qn` reproduces its `srt_a` walk.
fn push_rc_3_transitions(track: &[UlOv], transitions: &mut Vec<UTrans>) {
    let mut start = 0usize;
    while start < track.len() {
        let rid = track[start].qn;
        let mut end = start + 1;
        while end < track.len() && track[end].qn == rid {
            end += 1;
        }
        for entry in &track[start..end] {
            if entry.qn == entry.tn {
                continue;
            }
            let weight = cal_trans_ov_w(entry);
            if weight <= 0.0 {
                continue;
            }
            transitions.push(UTrans {
                query_start: entry.qs,
                query_end: entry.qe,
                query_identifier: entry.qn,
                target_start: entry.ts,
                target_end: entry.te,
                target_identifier: entry.tn,
                // The reference leaves `occ` unset here; the following
                // `clean_u_trans_t_idx_filter_adv` resets it before any read.
                occurrence: 0,
                network_weight: weight,
                flag: RC_3,
                reverse: entry.rev,
                deleted: false,
            });
        }
        start = end;
    }
}

/// `trans_base_infer` (`inter.cpp` 22574): the base-level transition
/// inference of `--dual-scaf` (and, later, Hi-C).
///
/// The function rebuilds the transition table (`clean_u_trans_t_idx_adv`),
/// builds the global minimizer index over the polished graph, maps every
/// active unitig against it (`worker_for_trans_ovlp_adv`), and appends the
/// collected `RC_3` transitions to `transitions`.
#[allow(clippy::too_many_arguments)]
pub fn trans_base_infer(
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
    bubble: &BubbleType,
    transitions: &mut Vec<UTrans>,
    transition_index: &mut Vec<u64>,
    reads: &ReadStore,
    options: &mut AssemblyOptions,
) {
    // `init_aux_table(); ha_opt_update_cov(&asm_opt, asm_opt.hom_cov);`
    options.update_coverage(options.homozygous_coverage);

    // `clean_u_trans_t_idx_adv(res, p->ug, p->rg); p->filter = res;`
    clean_u_trans_t_idx_adv(transitions, transition_index, unitig_graph, read_graph);

    let index = build_trans_base_index(unitig_graph, bubble, options);

    let lchain = set_lchain_parameters(false, options.transition_k_mer_length as u32);
    let chain_parameters = InputChainParameters {
        maximum_chain_count: options.maximum_chain_count as u64,
        maximum_skip: lchain.maximum_skip,
        maximum_iterations: lchain.maximum_iterations,
        maximum_distance: lchain.maximum_distance,
        gap_penalty: lchain.gap_penalty,
        skip_penalty: lchain.skip_penalty,
        band_width_rate: 1.0 - f64::from(options.transition_base_rate),
        band_width_rate_sec: 1.0 - f64::from(options.secondary_transition_base_rate),
        append_both_ends: true,
        generate_offset: true,
    };
    let error_low = 1.0 - f64::from(options.transition_base_rate);
    let error_high = 1.0 - f64::from(options.secondary_transition_base_rate);
    let k_mer_length = i64::from(options.transition_k_mer_length);
    let high_occurrence = (options.polyploidy + 1) as u32;
    // `init_ug_trans_t(..., 0.85, ..., 512, 3, 1, ...)`
    let sec_cutoff = 0.85f64;
    let chain_cut = 3u64;

    let target_lengths: Vec<u32> = unitig_graph.unitigs.iter().map(|u| u.length).collect();
    let targets = TargetSource::Unitigs {
        unitigs: &unitig_graph.unitigs,
        reads,
    };

    let mut track: Vec<UlOv> = Vec::new();
    let mut scratch: Vec<u64> = Vec::new();
    let mut extz = BitExtz::new();

    for rid in 0..unitig_graph.unitigs.len() as u32 {
        if skipped(unitig_graph, bubble, rid as usize) {
            continue;
        }
        let start = index.offsets[rid as usize] as usize;
        let end = index.offsets[rid as usize + 1] as usize;
        let query = &unitig_graph.unitigs[rid as usize].sequence;
        let query_length = u64::from(unitig_graph.unitigs[rid as usize].length);

        // `if(!gen_trans_adaptive_aln(..., fi, ...)) gen_trans_adaptive_aln(..., NULL, ...)`
        let first = gen_trans_adaptive_aln(
            rid,
            query,
            query_length,
            &index.entries[start..end],
            &index,
            &target_lengths,
            unitig_graph,
            transitions,
            transition_index,
            &targets,
            &chain_parameters,
            true,
            error_low,
            error_high,
            high_occurrence,
            sec_cutoff,
            k_mer_length,
            chain_cut,
            &mut track,
            &mut scratch,
            &mut extz,
        );
        if !first {
            gen_trans_adaptive_aln(
                rid,
                query,
                query_length,
                &index.entries[start..end],
                &index,
                &target_lengths,
                unitig_graph,
                transitions,
                transition_index,
                &targets,
                &chain_parameters,
                false,
                error_high,
                error_high,
                high_occurrence,
                sec_cutoff,
                k_mer_length,
                chain_cut,
                &mut track,
                &mut scratch,
                &mut extz,
            );
        }
    }

    // The reference sizes `res` to `res->n + track entries`; the push loop
    // then appends the `RC_3` transitions in query order.
    let mut appended: Vec<UTrans> = Vec::with_capacity(track.len());
    push_rc_3_transitions(&track, &mut appended);
    transitions.reserve(track.len());
    transitions.extend_from_slice(&appended);

    if let Some(directory) = crate::assembly::probe::dump_directory() {
        let mut text = String::new();
        for (i, t) in transitions.iter().enumerate() {
            text.push_str(&format!(
                "T\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\n",
                i,
                t.query_identifier,
                t.query_start,
                t.query_end,
                t.target_identifier,
                t.target_start,
                t.target_end,
                u8::from(t.reverse),
                t.flag,
                t.occurrence,
                u8::from(t.deleted),
                t.network_weight
            ));
        }
        text.push_str(&format!("N\t{}\t{}\n", transitions.len(), track.len()));
        let _ = std::fs::write(directory.join("trans_base_trans.txt"), text);
    }
}

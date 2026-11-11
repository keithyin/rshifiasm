//! Ultra-long (`--ul`) index construction, the counterpart of `dedup_HiFis`
//! and the three indices it builds (`inter.cpp`).
//!
//! Ported here (line numbers in `inter.cpp`):
//! - `push_coverage_track` (20616): per-unitig coverage track.
//! - `ul_contain_gen` (20717): the unitig containment table `ul_contain`.
//! - `debug_append_inexact_edges` (20932): asymmetric-edge pruning.
//! - `append_inexact_edges` (20962): re-inserts inexact edges into the unitig
//!   graph.
//! - `update_gen_r_contain` (21049): both `is_src_cc` passes of the read
//!   containment index.
//! - `update_ug_uo_t` (21116): the `ou` field of every unitig-graph arc.
//! - `gen_r_contain` (21133): the read containment index `ucov_t`.
//! - `gen_cov_track` (21173): the unitig coverage track `ucov_t`.
//! - `dedup_HiFis` (21185): the whole HiFi-read deduplication and index build.
//! - `gen_r_ug_idx` with `ct != NULL` (21278): the read-to-unitig index; the
//!   `ct`-free variant stays in [`crate::assembly::ul`].
//!
//! The C++ code mutates `ma_hit_t::cc` inside `gen_r_contain`; `PafHit` has no
//! such field, so a parallel `Vec<Vec<u32>>` with the same shape as the overlap
//! lists carries the containment counter instead and is handed to
//! [`update_ug_uo_t`] afterwards.
//!
//! Deferred: the coverage-based coverage track window selection of
//! `push_coverage_track`'s commented-out block, `check_if_fully_contain`
//! (disabled in the reference, it returns 1 immediately), `hpc_g_gen`,
//! `gen_ul_idx`/`gen_ul_idx_t` glue and the `is_scaf` index variants.

use crate::assembly::asg::{Asg, AsgArc, MaUtg, OU_MASK};
use crate::assembly::graph_build::asg_arc_del_trans;
use crate::assembly::hit::{ma_hit_to_arc, MA_HT_QCONT, MA_HT_TCONT};
use crate::assembly::ul::{UtgRid, UtgRidDt};
use crate::assembly::unitig::{ma_ug_gen, UnitigGraph};
use crate::command_line::options::{AssemblyOptions, FLAG_NO_HPC};
use crate::container::ksort::{radix_sort_by_key_32, radix_sort_u64};
use crate::correction::ec_driver::PafList;
use crate::correction::target::unitig_sub_region;
use crate::counting::filter_table::KmerFilterTable;
use crate::counting::kmer_counter::{
    for_each_homopolymer_compressed_kmer, for_each_kmer, KmerCountTable, MAXIMUM_COUNT,
};
use crate::counting::minimizer::sketch;
use crate::counting::position_index::{PositionIndex, PositionRecord};
use crate::logging;
use crate::overlap::paf::PafHit;
use crate::reads::ReadStore;
use crate::sequence::encoding::PackedSequence;
use std::io::{Read, Write};

/// `utg_ct_t` (`Overlaps.h` 278): one contained read inside a unitig.
///
/// All three fields are `uint32_t` in the reference; the coordinates are
/// truncated when written and the truncation is reproduced here.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct UtgCt {
    /// `(read << 1) | orientation`.
    pub x: u32,
    /// Start on the unitig.
    pub s: u32,
    /// End on the unitig.
    pub e: u32,
}

/// `ul_contain` (`Overlaps.h` 295).
#[derive(Clone, Debug, Default)]
pub struct UlContain {
    /// `idx[t] = (first entry << 32) | entry count`. The reference allocates
    /// only `ug->u.n` entries.
    pub idx: Vec<u64>,
    pub rids: Vec<UtgCt>,
    /// `is_c`, one flag per read vertex.
    pub is_c: Vec<u8>,
}

/// `ucov_t` (`Overlaps.h` 281).
#[derive(Clone, Debug, Default)]
pub struct UCov {
    pub idx: Vec<u32>,
    /// Packed `(<coordinate> << 32) | <depth>` words, sorted by coordinate.
    pub interval: Vec<u64>,
}

/// `ul_idx_t` (`Overlaps.h` 350), reduced to the fields this path fills in.
#[derive(Clone, Debug, Default)]
pub struct UlIndex {
    pub unitig: UnitigGraph,
    pub contain: UlContain,
    pub coverage: UCov,
    pub read_contain: UCov,
    pub read_unitig: UtgRid,
}

/// The reference threads these values through every helper as `ug_opt_t` plus
/// the global `asm_opt`.
struct IndexContext<'a> {
    options: &'a AssemblyOptions,
    /// `R_INF` read lengths.
    read_lengths: &'a [u64],
    /// `uopt->sources`.
    sources: &'a [PafList],
    /// `uopt->min_ovlp`, an `int64_t` in the reference.
    minimum_overlap: i64,
    /// `uopt->max_hang`, an `int64_t` in the reference.
    maximum_hang: i64,
}

impl<'a> IndexContext<'a> {
    fn new(options: &'a AssemblyOptions, read_lengths: &'a [u64], sources: &'a [PafList]) -> Self {
        Self {
            options,
            read_lengths,
            sources,
            minimum_overlap: i64::from(options.minimum_overlap_length),
            maximum_hang: i64::from(options.maximum_hang_length),
        }
    }

    /// `ma_hit2arc` with the caller-chosen read lengths.
    fn hit_to_arc(&self, hit: &PafHit, query_length: i32, target_length: i32) -> (i32, AsgArc) {
        ma_hit_to_arc(
            hit,
            query_length,
            target_length,
            self.maximum_hang as i32,
            self.options.maximum_hang_rate,
            self.minimum_overlap as i32,
        )
    }

    /// The two `(Get_qe - Get_qs) < min_ovlp` rejections shared by every loop.
    fn overlap_spans_ok(&self, hit: &PafHit) -> bool {
        i64::from(hit.query_end.wrapping_sub(hit.query_start)) >= self.minimum_overlap
            && i64::from(hit.target_end.wrapping_sub(hit.target_start)) >= self.minimum_overlap
    }
}

/// `get_specific_overlap` (`Overlaps.cpp` 1075): the mirror entry of `qn -> tn`
/// inside `tn`'s overlap list.
fn get_specific_overlap(list: &PafList, qn: u32, tn: u32) -> Option<usize> {
    list.hits
        .iter()
        .position(|hit| hit.target_identifier == tn && hit.query_identifier == qn)
}

/// `push_coverage_track` (`inter.cpp` 20616): appends one unitig's coverage
/// track and turns the interval list into `(coordinate, depth)` pairs.
fn push_coverage_track(
    coverage: &mut UCov,
    contain: &UlContain,
    uid: usize,
    unitig: &MaUtg,
    ctx: &IndexContext,
) {
    coverage.idx[uid] = coverage.interval.len() as u32;

    let packed = contain.idx[uid];
    let entry_count = packed as u32 as usize;
    let entry_start = (packed >> 32) as usize;
    for z in 0..entry_count {
        let entry = contain.rids[entry_start + z];
        // `ct_a[z].s` is a `uint32_t`, so the shift happens in 32 bits.
        coverage.interval.push(u64::from(entry.s.wrapping_shl(1)));
        coverage
            .interval
            .push(u64::from((entry.e.wrapping_shl(1)) | 1));
    }

    let mut l: u64 = 0;
    for &read in &unitig.reads {
        coverage.interval.push(l << 1);
        let read_length = ctx.read_lengths[(read >> 33) as usize];
        coverage.interval.push(((l + read_length) << 1) | 1);
        l += u64::from(read as u32);
    }
    coverage.idx[uid + 1] = coverage.interval.len() as u32;

    let begin = coverage.idx[uid] as usize;
    let end = coverage.interval.len();
    radix_sort_u64(&mut coverage.interval[begin..end]);
    let mut depth: u64 = 0;
    for k in begin..end {
        if coverage.interval[k] & 1 != 0 {
            depth = depth.wrapping_sub(1);
        } else {
            depth += 1;
        }
        let coordinate = (coverage.interval[k] >> 1) << 32;
        coverage.interval[k] = coordinate.wrapping_add(depth);
    }
}

/// `ul_contain_gen` (`inter.cpp` 20717): the containment table of every
/// unitig, one entry per contained read plus a 5%-tolerant dedup pass.
fn ul_contain_gen(
    unitigs: &[MaUtg],
    rg: &Asg,
    ctx: &IndexContext,
    is_el: bool,
    is_del: bool,
) -> UlContain {
    let mut contain = UlContain {
        idx: vec![0u64; unitigs.len()],
        rids: Vec::new(),
        is_c: vec![0u8; rg.sequences.len()],
    };

    for (t, unitig) in unitigs.iter().enumerate() {
        let start = contain.rids.len();
        contain.idx[t] = (start as u64) << 32;

        let mut l: u64 = 0;
        for &read in &unitig.reads {
            let i = (read >> 33) as usize;
            let orientation = (read >> 32) & 1;
            for hit in &ctx.sources[i].hits {
                if is_el && !hit.exact_region {
                    continue;
                }
                if is_del && !hit.deleted {
                    continue;
                }
                let qn = hit.query_identifier as usize;
                let tn = hit.target_identifier as usize;
                if !rg.sequences[tn].deleted {
                    continue;
                }
                if !ctx.overlap_spans_ok(hit) {
                    continue;
                }
                let (r, _) = ctx.hit_to_arc(
                    hit,
                    rg.sequences[qn].length as i32,
                    rg.sequences[tn].length as i32,
                );
                if r != MA_HT_TCONT {
                    // `tn` is contained
                    continue;
                }
                contain.is_c[qn] = 1;
                let (qs, qe) = if orientation == 0 {
                    (u64::from(hit.query_start), u64::from(hit.query_end))
                } else {
                    let read_length = ctx.read_lengths[i];
                    (
                        read_length - u64::from(hit.query_end),
                        read_length - u64::from(hit.query_start),
                    )
                };
                let (qs, qe) = (qs + l, qe + l);
                contain.rids.push(UtgCt {
                    x: hit.target_identifier.wrapping_shl(1)
                        | u32::from(hit.reverse != (orientation == 1)),
                    s: qs as u32,
                    e: qe as u32,
                });
            }
            l += u64::from(read as u32);
        }

        // `radix_sort_utg_ct_t_x_srt`: the entries of this unitig by `x`.
        radix_sort_by_key_32(&mut contain.rids[start..], &|entry: &UtgCt| entry.x);

        // The 5%-tolerant overlap dedup, preserving the reference's backward
        // scan across earlier unitigs.
        let mut kept = start;
        let mut group = start;
        let mut k = start + 1;
        while k <= contain.rids.len() {
            if k == contain.rids.len() || contain.rids[k].x != contain.rids[group].x {
                for z in group..k {
                    let mut r = kept as i64 - 1;
                    while r >= 0 && contain.rids[r as usize].x == contain.rids[z].x {
                        let candidate = contain.rids[z];
                        let previous = contain.rids[r as usize];
                        let overlap = if candidate.e.min(previous.e) > candidate.s.max(previous.s) {
                            candidate.e.min(previous.e) - candidate.s.max(previous.s)
                        } else {
                            0
                        };
                        if overlap != 0 {
                            let o_z = u64::from(candidate.e.wrapping_sub(candidate.s));
                            let o_r = u64::from(previous.e.wrapping_sub(previous.s));
                            let o_o = o_z.min(o_r) as f64;
                            let overlap = f64::from(overlap);
                            if overlap <= o_o * 1.05 && overlap >= o_o * 0.95 {
                                break;
                            }
                        }
                        r -= 1;
                    }
                    if r >= 0 && contain.rids[r as usize].x == contain.rids[z].x {
                        continue;
                    }
                    contain.rids[kept] = contain.rids[z];
                    kept += 1;
                }
                group = k;
            }
            k += 1;
        }
        contain.rids.truncate(kept);

        // `radix_sort_utg_ct_t_s_srt`: this unitig's entries by `s`.
        radix_sort_by_key_32(&mut contain.rids[start..], &|entry: &UtgCt| entry.s);
        contain.idx[t] |= (contain.rids.len() - start) as u64;
    }

    contain
}

/// `gen_cov_track` (`inter.cpp` 21173).
fn gen_cov_track(unitigs: &[MaUtg], contain: &UlContain, ctx: &IndexContext) -> UCov {
    let mut coverage = UCov {
        idx: vec![0u32; unitigs.len() + 1],
        interval: Vec::new(),
    };
    let mut bases: u64 = 0;
    for (i, unitig) in unitigs.iter().enumerate() {
        bases += u64::from(unitig.length);
        push_coverage_track(&mut coverage, contain, i, unitig, ctx);
    }
    eprintln!("[M::gen_cov_track::] # bases: {}", bases);
    coverage
}

/// `gen_r_contain` (`inter.cpp` 21133): the read containment index.
///
/// Returns the index and the `cc` counter array threaded through in place of
/// `ma_hit_t::cc`.
#[allow(clippy::too_many_arguments)]
fn gen_r_contain(
    unitig: &mut UnitigGraph,
    rg: &Asg,
    ctx: &IndexContext,
    n_read: usize,
    is_el: bool,
    is_del: bool,
) -> UCov {
    let mut coverage = UCov {
        idx: vec![0u32; n_read + 1],
        interval: Vec::new(),
    };
    let mut cc: Vec<Vec<u32>> = ctx
        .sources
        .iter()
        .map(|list| vec![0u32; list.hits.len()])
        .collect();

    for i in 0..n_read {
        coverage.idx[i] = coverage.interval.len() as u32;
        if rg.sequences[i].deleted {
            continue;
        }
        for hit in &ctx.sources[i].hits {
            if is_el && !hit.exact_region {
                continue;
            }
            if is_del && !hit.deleted {
                continue;
            }
            if !ctx.overlap_spans_ok(hit) {
                continue;
            }
            let qn = hit.query_identifier as usize;
            let tn = hit.target_identifier as usize;
            if !rg.sequences[tn].deleted {
                continue;
            }
            let (r, _) = ctx.hit_to_arc(
                hit,
                ctx.read_lengths[qn] as i32,
                ctx.read_lengths[tn] as i32,
            );
            if r == MA_HT_TCONT {
                // `tn` is contained
                let start =
                    (u64::from(hit.query_start) << 1) << 32 | u64::from(hit.target_identifier);
                let end =
                    ((u64::from(hit.query_end) << 1) | 1) << 32 | u64::from(hit.target_identifier);
                coverage.interval.push(start);
                coverage.interval.push(end);
            }
        }
    }
    coverage.idx[n_read] = coverage.interval.len() as u32;

    // `kt_for(update_gen_r_contain, ..., is_src_cc = 1)` then `is_src_cc = 0`;
    // both passes must finish before the next one starts.
    for i in 0..n_read {
        read_contain_first_pass(&coverage, &mut cc, rg, ctx, i);
    }
    for i in 0..n_read {
        read_contain_second_pass(&mut coverage, &mut cc, ctx.sources, i);
    }

    update_ug_uo_t(&mut unitig.graph, &unitig.unitigs, rg, ctx, &cc);
    coverage
}

/// `update_gen_r_contain` with `is_src_cc == 1` (`inter.cpp` 21049).
fn read_contain_first_pass(
    coverage: &UCov,
    cc: &mut [Vec<u32>],
    rg: &Asg,
    ctx: &IndexContext,
    i: usize,
) {
    let begin = coverage.idx[i] as usize;
    let a_n = (coverage.idx[i + 1] - coverage.idx[i]) as usize;
    for z in 0..ctx.sources[i].hits.len() {
        cc[i][z] = 0;
        if a_n == 0 || rg.sequences[i].deleted {
            continue;
        }
        let hit = ctx.sources[i].hits[z];
        let qn = hit.query_identifier;
        let tn = hit.target_identifier;
        if qn > tn {
            continue;
        }
        // `aux.is_el = aux.is_del = 0` for both passes.
        if !ctx.overlap_spans_ok(&hit) {
            continue;
        }
        if rg.sequences[tn as usize].deleted {
            continue;
        }
        let b_begin = coverage.idx[tn as usize] as usize;
        let b_n = (coverage.idx[tn as usize + 1] - coverage.idx[tn as usize]) as usize;
        if b_n == 0 {
            continue;
        }
        let (r, _) = ctx.hit_to_arc(
            &hit,
            ctx.read_lengths[qn as usize] as i32,
            ctx.read_lengths[tn as usize] as i32,
        );
        if r < 0 {
            continue;
        }
        let qs = u64::from(hit.query_start);
        let qe = u64::from(hit.query_end);
        let mut k = 0usize;
        while k < a_n {
            let cs = coverage.interval[begin + k] >> 33;
            let ce = coverage.interval[begin + k + 1] >> 33;
            debug_assert!(rg.sequences[coverage.interval[begin + k] as u32 as usize].deleted);
            // 128 is the offset allowed for indels.
            if qs <= cs + 128 && qe + 128 >= ce {
                let target = coverage.interval[begin + k] as u32;
                let mut ti = 0usize;
                while ti < b_n {
                    if coverage.interval[b_begin + ti] as u32 == target {
                        let mut sum = u64::from(cc[i][z]) + (ce - cs);
                        if sum > 0x3fff_ffff {
                            sum = 0x3fff_ffff;
                        }
                        cc[i][z] = sum as u32;
                        break;
                    }
                    ti += 2;
                }
            }
            k += 2;
        }
    }
}

/// `update_gen_r_contain` with `is_src_cc == 0` (`inter.cpp` 21049): turns the
/// interval list into depth words and copies `cc` to the mirror entries.
fn read_contain_second_pass(
    coverage: &mut UCov,
    cc: &mut [Vec<u32>],
    sources: &[PafList],
    i: usize,
) {
    let begin = coverage.idx[i] as usize;
    let end = coverage.idx[i + 1] as usize;
    radix_sort_u64(&mut coverage.interval[begin..end]);
    let mut depth: u64 = 0;
    for k in begin..end {
        if (coverage.interval[k] >> 32) & 1 != 0 {
            depth = depth.wrapping_sub(1);
        } else {
            depth += 1;
        }
        let coordinate = (coverage.interval[k] >> 33) << 32;
        coverage.interval[k] = coordinate.wrapping_add(depth);
    }

    for z in 0..sources[i].hits.len() {
        let hit = sources[i].hits[z];
        let qn = hit.query_identifier;
        let tn = hit.target_identifier;
        if qn > tn {
            continue;
        }
        let value = cc[i][z];
        if value == 0 {
            continue;
        }
        if let Some(index) = get_specific_overlap(&sources[tn as usize], tn, qn) {
            cc[tn as usize][index] = value;
        }
    }
}

/// `update_ug_uo_t` (`inter.cpp` 21116): fills the `ou` field of every unitig
/// arc from the read containment counter.
fn update_ug_uo_t(
    graph: &mut Asg,
    unitigs: &[MaUtg],
    rg: &Asg,
    ctx: &IndexContext,
    cc: &[Vec<u32>],
) {
    for arc_index in 0..graph.arcs.len() {
        let uv = graph.arcs[arc_index].query_vertex();
        let uw = graph.arcs[arc_index].v;
        graph.arcs[arc_index].occurrence = 0;
        let source = &unitigs[(uv >> 1) as usize];
        let target = &unitigs[(uw >> 1) as usize];
        if source.circular || target.circular {
            continue;
        }
        let v = if uv & 1 != 0 {
            source.start ^ 1
        } else {
            source.end ^ 1
        };
        let w = if uw & 1 != 0 {
            target.end
        } else {
            target.start
        };
        let read = v >> 1;
        let neighbour = w >> 1;
        let list = &ctx.sources[read as usize];
        let mut found = false;
        for k in 0..list.hits.len() {
            let hit = list.hits[k];
            if hit.query_identifier == read && hit.target_identifier == neighbour {
                let (r, t) = ctx.hit_to_arc(
                    &hit,
                    rg.sequences[read as usize].length as i32,
                    rg.sequences[neighbour as usize].length as i32,
                );
                if r < 0 {
                    continue;
                }
                if t.query_vertex() != v || t.v != w {
                    continue;
                }
                let value = cc[read as usize][k];
                graph.arcs[arc_index].occurrence = if value > OU_MASK {
                    OU_MASK as u16
                } else {
                    value as u16
                };
                found = true;
                break;
            }
        }
        debug_assert!(found, "unitig arc without a matching read overlap");
    }
}

/// `get_ug_edge_src` (`inter.cpp` 3381): the read overlap behind one unitig
/// arc, or `None` when it is missing.
fn get_ug_edge_src<'a>(
    unitigs: &[MaUtg],
    ctx: &'a IndexContext,
    uv: u32,
    uw: u32,
) -> Option<&'a PafHit> {
    let source = &unitigs[(uv >> 1) as usize];
    let target = &unitigs[(uw >> 1) as usize];
    if source.circular || target.circular {
        return None;
    }
    let v = if uv & 1 != 0 {
        source.start ^ 1
    } else {
        source.end ^ 1
    };
    let w = if uw & 1 != 0 {
        target.end
    } else {
        target.start
    };
    let list = &ctx.sources[(v >> 1) as usize];
    for hit in &list.hits {
        if hit.query_identifier == (v >> 1) && hit.target_identifier == (w >> 1) {
            let (r, t) = ctx.hit_to_arc(
                hit,
                ctx.read_lengths[(v >> 1) as usize] as i32,
                ctx.read_lengths[(w >> 1) as usize] as i32,
            );
            if r < 0 {
                continue;
            }
            if t.query_vertex() != v || t.v != w {
                continue;
            }
            return Some(hit);
        }
    }
    None
}

/// `debug_append_inexact_edges` (`inter.cpp` 20932): drops the arcs whose
/// reverse complement is missing and reports disconnected edges.
fn debug_append_inexact_edges(unitig: &mut UnitigGraph, ctx: &IndexContext) {
    let mut asymmetric = 0u32;
    let mut disconnected = 0u32;
    for z in 0..unitig.graph.arcs.len() {
        if unitig.graph.arcs[z].deleted {
            continue;
        }
        let uv = unitig.graph.arcs[z].query_vertex();
        let uw = unitig.graph.arcs[z].v;
        if get_ug_edge_src(&unitig.unitigs, ctx, uv, uw).is_none() {
            disconnected += 1;
        }
        let v = uw ^ 1;
        let w = uv ^ 1;
        let (start, count) = unitig.graph.arc_range(v);
        let mut k = 0usize;
        while k < count as usize {
            let arc = unitig.graph.arcs[start + k];
            if !arc.deleted && arc.v == w {
                break;
            }
            k += 1;
        }
        if k == count as usize {
            unitig.graph.arcs[z].deleted = true;
            asymmetric += 1;
        }
    }

    if asymmetric != 0 || disconnected != 0 {
        unitig.graph.cleanup();
        eprintln!(
            "[M::debug_append_inexact_edges] # asymm edges: {}, # disconnect edges: {}",
            asymmetric, disconnected
        );
    }
}

/// One half of `append_inexact_edges` (`inter.cpp` 20962): walks the overlap
/// list of the read behind the unitig end `v` and re-inserts the inexact arcs.
#[allow(clippy::too_many_arguments)]
fn append_inexact_phase(
    unitigs: &[MaUtg],
    graph: &mut Asg,
    rg: &Asg,
    ctx: &IndexContext,
    idx: &[u32],
    v: u32,
    ut_v: u32,
    inserted: &mut u64,
) {
    let list = &ctx.sources[(v >> 1) as usize];
    for hit in &list.hits {
        if hit.exact_region {
            // we just need inexact edges
            continue;
        }
        let qn = hit.query_identifier as usize;
        let tn = hit.target_identifier;
        let tu = idx[tn as usize];
        let mut ut_w = u32::MAX;
        if tu == u32::MAX || graph.sequences[tu as usize].deleted {
            continue;
        }
        if !ctx.overlap_spans_ok(hit) {
            continue;
        }
        let (r, t) = ctx.hit_to_arc(
            hit,
            rg.sequences[qn].length as i32,
            rg.sequences[tn as usize].length as i32,
        );
        if r < 0 || t.query_vertex() != v {
            continue;
        }
        let target = &unitigs[tu as usize];
        if t.v == target.start {
            ut_w = tu << 1;
        }
        if t.v == target.end {
            ut_w = (tu << 1) + 1;
        }
        if ut_w == u32::MAX {
            continue;
        }
        let mut arc = t;
        arc.ul = (u64::from(ut_v) << 32) + u64::from(t.ul as u32);
        arc.v = ut_w;
        graph.push_arc(arc);
        *inserted += 1;
    }
}

/// `append_inexact_edges` (`inter.cpp` 20962).
fn append_inexact_edges(unitig: &mut UnitigGraph, rg: &Asg, ctx: &IndexContext) {
    let n_read = ctx.read_lengths.len();
    let mut idx = vec![u32::MAX; n_read];
    for (z, u) in unitig.unitigs.iter().enumerate() {
        if u.circular {
            continue;
        }
        idx[(u.start >> 1) as usize] = z as u32;
        idx[(u.end >> 1) as usize] = z as u32;
    }

    let mut inserted: u64 = 0;
    let UnitigGraph { unitigs, graph } = unitig;
    for z in 0..unitigs.len() {
        let u = &unitigs[z];
        if u.circular {
            continue;
        }
        append_inexact_phase(
            unitigs,
            graph,
            rg,
            ctx,
            &idx,
            u.end ^ 1,
            (z as u32) << 1,
            &mut inserted,
        );
        append_inexact_phase(
            unitigs,
            graph,
            rg,
            ctx,
            &idx,
            u.start ^ 1,
            ((z as u32) << 1) + 1,
            &mut inserted,
        );
    }

    if inserted != 0 {
        // `free(ug->g->idx); ug->g->idx = 0; ug->g->is_srt = 0;`
        graph.index.clear();
        graph.is_sorted = false;
        graph.cleanup();
    }

    debug_append_inexact_edges(unitig, ctx);
    eprintln!(
        "[M::append_inexact_edges] # inserted inexact edges: {}",
        inserted
    );
}

/// `gen_r_ug_idx` (`inter.cpp` 21278) with a containment table.
///
/// `keep_offset` mirrors the reference: the `--ul` index keeps the base offset
/// inside the unitig (`off = position`), the scaffolding index keeps the read
/// index (`off = k`).
pub fn gen_r_ug_idx_with_contain(
    unitigs: &[MaUtg],
    contain: &UlContain,
    read_count: usize,
    keep_offset: bool,
) -> UtgRid {
    let mut idx = vec![0u32; read_count + 1];
    for (i, unitig) in unitigs.iter().enumerate() {
        for &packed in &unitig.reads {
            idx[(packed >> 33) as usize] += 1;
        }
        let count = contain.idx[i] as u32 as usize;
        let start = (contain.idx[i] >> 32) as usize;
        for k in 0..count {
            idx[(contain.rids[start + k].x >> 1) as usize] += 1;
        }
    }
    let mut l = 0u32;
    for k in 0..=read_count {
        let count = idx[k];
        idx[k] = l;
        l += count;
    }

    let total = idx[read_count] as usize;
    let mut p = vec![UtgRidDt::default(); total];
    for (i, unitig) in unitigs.iter().enumerate() {
        let mut l = 0u32;
        for (k, &packed) in unitig.reads.iter().enumerate() {
            let rid = (packed >> 33) as usize;
            let base = idx[rid] as usize;
            let a_n = (idx[rid + 1] - idx[rid]) as usize;
            if a_n > 0 {
                let last = base + a_n - 1;
                let slot = if p[last].off as usize == a_n - 1 {
                    last
                } else {
                    let cursor = base + p[last].off as usize;
                    p[last].off += 1;
                    cursor
                };
                p[slot].u = ((i as u32) << 1) | (((packed >> 32) & 1) as u32);
                p[slot].pos = l;
                p[slot].off = if keep_offset { k as u32 } else { l };
            }
            l = l.wrapping_add(packed as u32);
        }

        let count = contain.idx[i] as u32 as usize;
        let start = (contain.idx[i] >> 32) as usize;
        for k in 0..count {
            let entry = contain.rids[start + k];
            let rid = (entry.x >> 1) as usize;
            let base = idx[rid] as usize;
            let a_n = (idx[rid + 1] - idx[rid]) as usize;
            if a_n == 0 {
                continue;
            }
            let last = base + a_n - 1;
            let slot = if p[last].off as usize == a_n - 1 {
                last
            } else {
                let cursor = base + p[last].off as usize;
                p[last].off += 1;
                cursor
            };
            p[slot].u = ((i as u32) << 1) | (entry.x & 1);
            p[slot].pos = entry.s;
            p[slot].off = if keep_offset { k as u32 } else { entry.s };
        }
    }

    UtgRid { idx, p }
}

/// `dedup_HiFis` (`inter.cpp` 21185): the HiFi contained-read removal, the
/// unitig graph and the three ultra-long indices.
pub fn dedup_hi_fis(
    options: &AssemblyOptions,
    reads: &ReadStore,
    sources: &[PafList],
    read_lengths: &[u64],
) -> UlIndex {
    // `dedup_HiFis(sl->uopt, 1, 0)`.
    let is_el = true;
    let is_del = false;
    let n_read = read_lengths.len();
    let ctx = IndexContext::new(options, read_lengths, sources);
    debug_assert_eq!(reads.read_lengths.len(), n_read);

    // `asg_init` plus one vertex per read.
    let mut rg = Asg::new();
    for i in 0..n_read {
        rg.set_sequence(i, reads.read_lengths[i] as u32, false);
    }
    rg.sequence_visited = vec![0u8; n_read * 2];

    // Pass one: mark the contained reads.
    for i in 0..n_read {
        if rg.sequences[i].deleted {
            continue;
        }
        for hit in &sources[i].hits {
            if is_el && !hit.exact_region {
                continue;
            }
            if is_del && !hit.deleted {
                continue;
            }
            let qn = hit.query_identifier as usize;
            let tn = hit.target_identifier as usize;
            if rg.sequences[qn].deleted || rg.sequences[tn].deleted {
                continue;
            }
            if !ctx.overlap_spans_ok(hit) {
                continue;
            }
            let (r, _) = ctx.hit_to_arc(
                hit,
                rg.sequences[qn].length as i32,
                rg.sequences[tn].length as i32,
            );
            if r == MA_HT_QCONT {
                rg.sequences[qn].deleted = true;
            } else if r == MA_HT_TCONT {
                rg.sequences[tn].deleted = true;
            }
            if rg.sequences[i].deleted {
                break;
            }
        }
    }

    // Pass two: the read graph.
    let mut cc_num: u64 = 0;
    for i in 0..n_read {
        if rg.sequences[i].deleted {
            cc_num += 1;
            continue;
        }
        for hit in &sources[i].hits {
            if is_el && !hit.exact_region {
                continue;
            }
            if is_del && !hit.deleted {
                continue;
            }
            let qn = hit.query_identifier as usize;
            let tn = hit.target_identifier as usize;
            if rg.sequences[qn].deleted || rg.sequences[tn].deleted {
                continue;
            }
            if !ctx.overlap_spans_ok(hit) {
                continue;
            }
            let (r, t) = ctx.hit_to_arc(
                hit,
                rg.sequences[qn].length as i32,
                rg.sequences[tn].length as i32,
            );
            if r >= 0 {
                rg.push_arc(t);
            }
        }
    }

    rg.cleanup();
    rg.symmetrize();
    asg_arc_del_trans(&mut rg, options.gap_fuzz);

    let mut unitig = ma_ug_gen(&mut rg);
    append_inexact_edges(&mut unitig, &rg, &ctx);

    let contain = ul_contain_gen(&unitig.unitigs, &rg, &ctx, is_el, is_del);
    let coverage = gen_cov_track(&unitig.unitigs, &contain, &ctx);
    let read_contain = gen_r_contain(&mut unitig, &rg, &ctx, n_read, is_el, is_del);
    let read_unitig = gen_r_ug_idx_with_contain(&unitig.unitigs, &contain, n_read, false);

    eprintln!(
        "[M::dedup_HiFis::] # unitigs: {}, # edges: {}, # cc_num: {}",
        unitig.unitigs.len(),
        unitig.graph.arcs.len(),
        cc_num
    );

    UlIndex {
        unitig,
        contain,
        coverage,
        read_contain,
        read_unitig,
    }
}

/// Materialises every unitig's sequence into a read store, the input of the
/// unitig-level minimizer index (`ha_ft_ul_gen`/`ha_pt_ul_gen` iterate
/// `ma_utg_v` through `retrieve_u_seq`).
pub fn materialize_unitig_reads(unitigs: &[MaUtg], reads: &ReadStore) -> ReadStore {
    let mut store = ReadStore::default();
    store.read_lengths = unitigs.iter().map(|unitig| unitig.length as u64).collect();
    store.read_names = (0..unitigs.len()).map(|id| format!("utg{id}")).collect();
    store.sequences = unitigs
        .iter()
        .enumerate()
        .map(|(id, unitig)| {
            unitig_sub_region(unitigs, reads, id, 0, unitig.length as i64, 0)
        })
        .map(|sequence| PackedSequence::from_ascii(&sequence))
        .collect();
    store.trio_labels = vec![0; unitigs.len()];
    store.qualities = vec![Vec::new(); unitigs.len()];
    store.total_reads = unitigs.len() as u64;
    store.total_bases = store.read_lengths.iter().sum();
    store
}

/// The unitig-level minimizer index of `uidx_l_build`: the high-occurrence
/// filter (`ha_ft_ul_gen`) and the position index (`ha_pt_ul_gen`).
pub struct UnitigMinimizerIndex {
    pub filter: KmerFilterTable,
    pub index: PositionIndex,
    /// The materialised unitig sequences, kept alive for the alignment pass.
    pub reads: ReadStore,
}

/// `ha_ft_ul_gen` (`htab.cpp` 1106): counts every k-mer of the unitigs,
/// drops the ones occurring at least `cutoff` times and builds the filter.
fn build_unitig_filter(
    options: &AssemblyOptions,
    reads: &ReadStore,
    cutoff: i64,
) -> KmerFilterTable {
    let use_homopolymer_compression = options.flags & FLAG_NO_HPC == 0;
    let mut counts = KmerCountTable::new(options.k_mer_length, options.bloom_filter_shift, 4);
    for read_id in 0..reads.total_reads as usize {
        let sequence = reads.read_sequence(read_id);
        if use_homopolymer_compression {
            for_each_homopolymer_compressed_kmer(&sequence, options.k_mer_length as u32, |hash| {
                counts.insert(hash);
            });
        } else {
            for_each_kmer(&sequence, options.k_mer_length as u32, |hash| {
                counts.insert(hash);
            });
        }
    }
    counts.set_total_bases(reads.total_bases);

    let mut boundary = cutoff;
    if boundary > MAXIMUM_COUNT as i64 - 1 {
        boundary = MAXIMUM_COUNT as i64 - 1;
    }
    let remaining = counts.shrink(boundary, MAXIMUM_COUNT as i64);
    eprintln!(
        "[M::ha_ct_shrink::{:.3}*{:.2}] ==> counted {} distinct minimizer k-mers",
        logging::realtime_seconds(),
        logging::cpu_usage(),
        remaining
    );
    let table = KmerFilterTable::from_counts(&counts, options.maximum_kmer_occurrence);
    eprintln!(
        "[M::ha_ft_ul_gen::{:.3}*{:.2}@{:.3}GB] ==> filtered out {} k-mers occurring {} or more times",
        logging::realtime_seconds(),
        logging::cpu_usage(),
        logging::peak_resident_memory_gigabytes(),
        table.len(),
        boundary
    );
    table
}

/// `ha_pt_ul_gen` (`htab.cpp` 1171): minimizer counting, shrink and the
/// two-pass position index over the unitig sequences.
fn build_unitig_position_index(
    options: &AssemblyOptions,
    reads: &ReadStore,
    filter: &KmerFilterTable,
) -> PositionIndex {
    let use_homopolymer_compression = options.flags & FLAG_NO_HPC == 0;
    let sketch_read = |read_id: usize| {
        let sequence = reads.read_sequence(read_id);
        sketch(
            &sequence,
            options.minimizer_window_size,
            options.k_mer_length,
            read_id as u32,
            use_homopolymer_compression,
            Some(filter),
            options.minimizer_sample_distance,
            options.minimizer_rewindow_size,
        )
        .minimizers
    };

    let mut counts = KmerCountTable::new(options.k_mer_length, 0, 0);
    let mut total_minimizers = 0u64;
    for read_id in 0..reads.total_reads as usize {
        let minimizers = sketch_read(read_id);
        for minimizer in &minimizers {
            counts.insert(minimizer.hash);
        }
        total_minimizers += minimizers.len() as u64;
    }
    counts.set_total_bases(reads.total_bases);
    if options.minimizer_window_size > 1 {
        eprintln!("[M::yak_count] collected {} minimizers", total_minimizers);
    }
    eprintln!(
        "[M::ha_pt_ul_gen::{:.3}*{:.2}] ==> counted {} distinct minimizer k-mers",
        logging::realtime_seconds(),
        logging::cpu_usage(),
        counts.distinct_count()
    );
    let remaining = counts.shrink(1, MAXIMUM_COUNT as i64 - 1);
    eprintln!(
        "[M::ha_ct_shrink::{:.3}*{:.2}] ==> counted {} distinct minimizer k-mers",
        logging::realtime_seconds(),
        logging::cpu_usage(),
        remaining
    );

    let mut index = PositionIndex::build(&counts, options.k_mer_length);
    for read_id in 0..reads.total_reads as usize {
        for minimizer in sketch_read(read_id) {
            index.insert(
                minimizer.hash,
                PositionRecord {
                    read_identifier: minimizer.read_identifier,
                    position: minimizer.position,
                    reverse: minimizer.reverse,
                    span: minimizer.span,
                },
            );
        }
    }
    if options.minimizer_window_size > 1 {
        eprintln!("[M::yak_count] collected {} minimizers", total_minimizers);
    }
    eprintln!(
        "[M::ha_pt_ul_gen::{:.3}*{:.2}] ==> indexed {} positions",
        logging::realtime_seconds(),
        logging::cpu_usage(),
        index.total_positions()
    );
    index
}

/// `uidx_l_build` (`inter.cpp` 515): the two unitig-level indices.
pub fn build_unitig_minimizer_index(
    options: &AssemblyOptions,
    unitigs: &[MaUtg],
    reads: &ReadStore,
    cutoff: i64,
) -> UnitigMinimizerIndex {
    let unitig_reads = materialize_unitig_reads(unitigs, reads);
    let filter = build_unitig_filter(options, &unitig_reads, cutoff);
    let index = build_unitig_position_index(options, &unitig_reads, &filter);
    UnitigMinimizerIndex {
        filter,
        index,
        reads: unitig_reads,
    }
}

/// `write_dbug` (`Overlaps.cpp` 14735): the unitig-graph fingerprint that
/// prefixes the caches which carry a `ma_ug_t`, so a stale cache is detected
/// when the graph changed.
pub fn write_unitig_fingerprint(
    unitigs: &[MaUtg],
    writer: &mut impl std::io::Write,
) -> std::io::Result<()> {
    writer.write_all(&(unitigs.len() as u64).to_le_bytes())?;
    for unitig in unitigs {
        writer.write_all(&unitig.length.to_le_bytes())?;
        writer.write_all(&(unitig.circular as u32).to_le_bytes())?;
        writer.write_all(&unitig.start.to_le_bytes())?;
        writer.write_all(&unitig.end.to_le_bytes())?;
        writer.write_all(&(unitig.reads.len() as u32).to_le_bytes())?;
        for read in &unitig.reads {
            writer.write_all(&read.to_le_bytes())?;
        }
    }
    Ok(())
}

/// `test_dbug` (`Overlaps.cpp` 14700): `true` when the fingerprint on disk
/// matches `unitigs`. On a mismatch the stream position is past the header,
/// so the caller abandons the whole cache.
pub fn test_unitig_fingerprint(
    unitigs: &[MaUtg],
    reader: &mut impl std::io::Read,
) -> std::io::Result<bool> {
    let mut buffer = [0u8; 8];
    reader.read_exact(&mut buffer)?;
    if u64::from_le_bytes(buffer) != unitigs.len() as u64 {
        return Ok(false);
    }
    for unitig in unitigs {
        let mut word = [0u8; 4];
        reader.read_exact(&mut word)?;
        if u32::from_le_bytes(word) != unitig.length {
            return Ok(false);
        }
        reader.read_exact(&mut word)?;
        if u32::from_le_bytes(word) != unitig.circular as u32 {
            return Ok(false);
        }
        reader.read_exact(&mut word)?;
        if u32::from_le_bytes(word) != unitig.start {
            return Ok(false);
        }
        reader.read_exact(&mut word)?;
        if u32::from_le_bytes(word) != unitig.end {
            return Ok(false);
        }
        reader.read_exact(&mut word)?;
        let read_count = u32::from_le_bytes(word) as usize;
        if read_count != unitig.reads.len() {
            return Ok(false);
        }
        let mut packed = vec![0u8; read_count * 8];
        reader.read_exact(&mut packed)?;
        for (index, value) in packed.chunks_exact(8).enumerate() {
            if u64::from_le_bytes(value.try_into().unwrap()) != unitig.reads[index] {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// `uidx_write` (`htab.cpp` 1558): the `%s.uidx.bin` cache, `ug`'s
/// fingerprint when the index belongs to a unitig graph.
pub fn uidx_write(
    index: &UnitigMinimizerIndex,
    file_name: &str,
    unitigs: Option<&[MaUtg]>,
) -> std::io::Result<()> {
    let path = format!("{}.uidx.bin", file_name);
    let file = std::fs::File::create(&path)?;
    let mut writer = std::io::BufWriter::new(file);
    if let Some(unitigs) = unitigs {
        write_unitig_fingerprint(unitigs, &mut writer)?;
    }
    writer.write_all(b"f")?;
    index.filter.write(&mut writer)?;
    writer.write_all(b"h")?;
    index.index.write(&mut writer)?;
    writer.flush()?;
    eprintln!("[M::uidx_write] Index has been written.");
    Ok(())
}

/// `uidx_load` (`htab.cpp` 1603): `None` when the cache is absent or stale.
///
/// The materialised unitig sequences are rebuilt from `unitigs`, which the
/// caller owns; the index itself never serialises them.
pub fn uidx_load(
    file_name: &str,
    unitigs: Option<&[MaUtg]>,
    reads: &ReadStore,
) -> std::io::Result<Option<UnitigMinimizerIndex>> {
    let path = format!("{}.uidx.bin", file_name);
    let Ok(file) = std::fs::File::open(&path) else {
        return Ok(None);
    };
    let mut reader = std::io::BufReader::new(file);
    if let Some(unitigs) = unitigs {
        if !test_unitig_fingerprint(unitigs, &mut reader)? {
            eprintln!("[M::uidx_load] Renew UL Index");
            return Ok(None);
        }
    }
    let mut mode = [0u8; 1];
    reader.read_exact(&mut mode)?;
    let filter = if mode[0] == b'f' {
        let filter = KmerFilterTable::read(&mut reader)?;
        reader.read_exact(&mut mode)?;
        filter
    } else {
        return Ok(None);
    };
    if mode[0] != b'h' {
        return Ok(None);
    }
    let index = match PositionIndex::read(&mut reader) {
        Ok(index) => index,
        Err(_) => return Ok(None),
    };
    eprintln!("[M::uidx_load] Index has been loaded.");
    Ok(Some(UnitigMinimizerIndex {
        filter,
        index,
        reads: materialize_unitig_reads(unitigs.unwrap_or(&[]), reads),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(qn: u32, qs: u32, qe: u32, tn: u32, ts: u32, te: u32, reverse: bool) -> PafHit {
        PafHit {
            query_identifier: qn,
            query_start: qs,
            query_end: qe,
            target_identifier: tn,
            target_start: ts,
            target_end: te,
            matching_length: 1,
            block_length: 0,
            reverse,
            without_large_indel: true,
            exact_region: true,
            exact_region_byte: 1,
            deleted: false,
        }
    }

    #[test]
    fn contained_read_is_dropped_and_all_indices_are_built() {
        let lengths = [1000u64, 1000, 300];
        let mut lists = vec![PafList::default(); 3];
        // A dovetail between read 0 and read 1, stored from both sides.
        lists[0].hits.push(hit(0, 800, 1000, 1, 0, 200, false));
        lists[1].hits.push(hit(1, 0, 200, 0, 800, 1000, false));
        // Read 2 (300 bp) sits inside read 0 at 200..500.
        lists[0].hits.push(hit(0, 200, 500, 2, 0, 300, false));
        lists[2].hits.push(hit(2, 0, 300, 0, 200, 500, false));

        let mut reads = ReadStore::default();
        reads.read_lengths = lengths.to_vec();
        let options = AssemblyOptions::default();

        let index = dedup_hi_fis(&options, &reads, &lists, &lengths);

        // One unitig of the two 1 kb reads; the contained read is gone.
        assert_eq!(index.unitig.unitigs.len(), 1);
        let unitig = &index.unitig.unitigs[0];
        assert_eq!(unitig.length, 1800);
        assert_eq!(unitig.reads.len(), 2);
        assert_eq!(unitig.reads[0] >> 33, 0);
        assert_eq!(unitig.reads[1] >> 33, 1);
        assert_eq!(index.unitig.graph.arcs.len(), 0);

        // The contained read as a `ul_contain` entry of the first unitig.
        assert_eq!(index.contain.idx, vec![1]);
        assert_eq!(index.contain.rids.len(), 1);
        assert_eq!(
            index.contain.rids[0],
            UtgCt {
                x: 2 << 1,
                s: 200,
                e: 500
            }
        );
        assert_eq!(index.contain.is_c[0], 1);

        // Coverage: the unitig reads plus the contained read's interval.
        assert_eq!(index.coverage.idx, vec![0, 6]);
        assert_eq!(
            index.coverage.interval,
            vec![
                (0u64 << 32) | 1,
                (200u64 << 32) | 2,
                (500u64 << 32) | 1,
                (800u64 << 32) | 2,
                (1000u64 << 32) | 1,
                (1800u64 << 32) | 0,
            ]
        );

        // Read containment: only read 0 covers the contained read 2.
        assert_eq!(index.read_contain.idx, vec![0, 2, 2, 2]);
        assert_eq!(
            index.read_contain.interval,
            vec![(200u64 << 32) | 1, (500u64 << 32) | 0]
        );

        // Read-to-unitig index, one entry per read.
        assert_eq!(index.read_unitig.idx, vec![0, 1, 2, 3]);
        let entries = &index.read_unitig.p;
        assert_eq!(entries[0].u, 0);
        assert_eq!(entries[0].pos, 0);
        assert_eq!(entries[1].pos, 800);
        assert_eq!(entries[2].u, 0);
        assert_eq!(entries[2].pos, 200);
    }

    #[test]
    fn inexact_edges_are_reinserted_between_unitig_ends() {
        let lengths = [1000u64, 1000, 1000, 1000];
        let mut lists = vec![PafList::default(); 4];

        fn link(lists: &mut [PafList], exact: bool, a: u32, b: u32) {
            let mut forward = hit(a, 800, 1000, b, 0, 200, false);
            forward.exact_region = exact;
            forward.exact_region_byte = u8::from(exact);
            let mut reverse = hit(b, 0, 200, a, 800, 1000, false);
            reverse.exact_region = exact;
            reverse.exact_region_byte = u8::from(exact);
            lists[a as usize].hits.push(forward);
            lists[b as usize].hits.push(reverse);
        }

        // A -> B exactly, B -> C and A -> D only inexactly.
        let mut a_to_b = hit(0, 800, 1000, 1, 0, 200, false);
        a_to_b.exact_region = true;
        lists[0].hits.push(a_to_b);
        lists[1].hits.push(hit(1, 0, 200, 0, 800, 1000, false));
        link(&mut lists, false, 1, 2);
        link(&mut lists, false, 0, 3);

        let mut reads = ReadStore::default();
        reads.read_lengths = lengths.to_vec();
        let options = AssemblyOptions::default();

        let index = dedup_hi_fis(&options, &reads, &lists, &lengths);

        // Only the exact A -> B overlap joins the read graph, so A+B form one
        // unitig and C and D stay singletons.
        assert_eq!(index.unitig.unitigs.len(), 3);
        assert_eq!(index.unitig.unitigs[0].reads.len(), 2);
        assert_eq!(index.unitig.unitigs[0].length, 1800);
        assert_eq!(index.unitig.unitigs[1].reads.len(), 1);
        assert_eq!(index.unitig.unitigs[2].reads.len(), 1);
        // The two inexact overlaps are re-inserted between the unitig ends.
        assert_eq!(index.unitig.graph.arcs.len(), 2);
        assert_eq!(index.unitig.graph.arcs[0].query_vertex(), 0);
        assert_eq!(index.unitig.graph.arcs[0].v, 2);
        assert_eq!(index.unitig.graph.arcs[1].query_vertex(), 3);
        assert_eq!(index.unitig.graph.arcs[1].v, 1);

        assert_eq!(index.read_unitig.idx, vec![0, 1, 2, 3, 4]);
        let entries = &index.read_unitig.p;
        assert_eq!((entries[0].u, entries[0].pos), (0, 0));
        assert_eq!((entries[1].u, entries[1].pos), (0, 800));
        assert_eq!((entries[2].u, entries[2].pos), (2, 0));
        assert_eq!((entries[3].u, entries[3].pos), (4, 0));
    }
}

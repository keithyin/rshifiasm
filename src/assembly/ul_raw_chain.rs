//! The raw-chain DP of the UL scall pass, the counterpart of `inter.cpp`
//! 7400-9218: `gl_rchain_lin_contain_adv` chains the per-unitig sub-overlaps
//! of one UL read into raw chains, threading through the contained reads of
//! each unitig (`connect_detect`), emitting the transitive alignments of the
//! skipped reads (`gen_trans_aln`), and flattening each chain
//! (`flat_contain_adv`). `select_clean_chain` deduplicates the chains and
//! renumbers their bodies, and `gen_rid_raw_chain` ties the pass to the scall
//! worker buffers.
//!
//! The C works on three aliased buffers (`ll->lo` the working copy, `ll->tk`
//! whose tail doubles as the chain-body output `ex`, and `ll->srt.a` whose
//! array doubles as the visibility flags); the port keeps them as separate
//! vectors, which is content-equivalent because every read in the C happens
//! after the corresponding write.
//!
//! Dead code skipped: `gl_rchain_lin` (7619), `gl_rchain_lin_contain` (8522),
//! `collapse_contain` (called from a disabled line), `prt_rid_raw_chain`/
//! `prt_all_chain` (debug printers).

use crate::assembly::asg::MaUtg;
use crate::assembly::ul::{infer_rovlp_reads, UlOv};
use crate::correction::ec_driver::PafList;
use crate::correction::rphase::SEC_MODE;
use crate::correction::ul_align::{
    gen_err_unaligned, get_rid_backward_cigar_err, Rtrace, RtraceIter,
};
use crate::correction::hc_align::BitExtz;
use crate::correction::target::TargetSource;
use crate::correction::window::WindowCigar;
use crate::assembly::hit::ma_hit_to_arc;
use crate::assembly::hit::{MA_HT_TCONT, MA_HT_QCONT};
use crate::overlap::region::ChainData;

/// `G_CHAIN_BW` (`inter.h` 7).
pub(crate) const G_CHAIN_BW: i64 = 16;
/// `G_CHAIN_INDEL` (`inter.h` 27), also in `ul_align` as `CHAIN_INDEL`.
pub(crate) const G_CHAIN_INDEL: u32 = 128;
/// `N_GCHAIN_RATE` (`inter.h` 29).
pub(crate) const N_GCHAIN_RATE: f64 = 0.04;
/// `UG_SKIP_N` (`inter.h` 18).
pub(crate) const UG_SKIP_N: i64 = 100;
/// `UG_ITER_N` (`inter.h` 19).
pub(crate) const UG_ITER_N: i64 = 5000;
/// `UG_DIS_N` (`inter.h` 20).
pub(crate) const UG_DIS_N: i64 = 50000;
/// `GC_OFFSET_POS` (`inter.cpp` 55).
const GC_OFFSET_POS: u64 = 8;
/// `rch_connect(x, i)` (`inter.cpp` 7764).
fn rch_connect(x: i32, i: i64) -> bool {
    ((x >> 2) as i64) == i && (x & 3) != 3
}

/// `ug_opt_t` as the raw-chain pass consumes it: the HiFi overlap lists plus
/// the two thresholds `gen_ug_opt_t` fills.
#[derive(Clone, Copy)]
pub struct UlGraphOptions<'a> {
    /// `uopt->sources`: the forward overlap lists of the HiFi reads.
    pub sources: &'a [PafList],
    /// `uopt->min_ovlp`.
    pub minimum_overlap: i32,
    /// `uopt->max_hang`.
    pub maximum_hang: i32,
    /// `asm_opt.max_hang_rate`, read by `ma_hit2arc`.
    pub maximum_hang_rate: f32,
}

/// The per-UL-read buffers of the scall worker, `glchain_t`.
#[derive(Default)]
pub struct GlChain {
    /// `lo`: the chain heads of the current read.
    pub lo: Vec<UlOv>,
    /// `tk`: every overlap of the current read; the raw chains replace the
    /// tail from `cha_idx` on.
    pub tk: Vec<UlOv>,
    /// `tc`: the trace entries of the current read's sub-overlaps.
    pub tc: Vec<Rtrace>,
    /// `srt.a`: the 64-bit scratch the DP shares with the visibility flags.
    pub srt: Vec<u64>,
}

/// `get_ecov_contain_adv` (`inter.cpp` 7576): checks the HiFi read graph for
/// an arc `v -> w` whose length corroborates the query distance `dq`, and
/// reports whether `w` is a contained read instead.
#[allow(clippy::too_many_arguments)]
pub(crate) fn get_ecov_contain_adv(
    options: &UlGraphOptions<'_>,
    reads: &crate::reads::read_store::ReadStore,
    v: u32,
    w: u32,
    band_width: i64,
    diff_ec_ul: f64,
    query_distance: i64,
    mut is_contain: Option<&mut bool>,
) -> bool {
    let mut distance = -1i64;
    let x = (v >> 1) as usize;
    if let Some(flag) = is_contain.as_deref_mut() {
        *flag = false;
    }
    let source = &options.sources[x];
    for hit in source.hits.iter().take(source.hits.len()) {
        let qn = hit.query_identifier;
        let tn = hit.target_identifier;
        if tn != w >> 1 {
            continue;
        }
        let (result, arc) = ma_hit_to_arc(
            hit,
            reads.read_length(qn as usize) as i32,
            reads.read_length(tn as usize) as i32,
            options.maximum_hang,
            options.maximum_hang_rate,
            options.minimum_overlap,
        );
        if result >= 0 {
            if (arc.ul >> 32) != u64::from(v) || u64::from(arc.v) != u64::from(w) {
                continue;
            }
            distance = i64::from(arc.overlap);
            if let Some(flag) = is_contain.as_deref_mut() {
                *flag = false;
            }
            break;
        } else if result == MA_HT_TCONT {
            // `tn` is contained in `qn`
            if hit.reverse == (((v ^ w) & 1) != 0) && tn == w >> 1 {
                distance = i64::from(hit.query_end.wrapping_sub(hit.query_start));
                let target_span = i64::from(hit.target_end.wrapping_sub(hit.target_start));
                if distance < target_span {
                    distance = target_span;
                }
                if let Some(flag) = is_contain.as_deref_mut() {
                    *flag = true;
                }
                break;
            }
        }
        let _ = MA_HT_QCONT;
    }

    if distance < 0 {
        return false;
    }
    let difference = if query_distance > distance {
        query_distance - distance
    } else {
        distance - query_distance
    };
    let mut allowed = query_distance.max(distance);
    allowed = (allowed as f64 * diff_ec_ul) as i64;
    if allowed < band_width {
        allowed = band_width;
    }
    difference <= allowed
}

/// `max_ovlp_src` (`inter.cpp` 1505): the longest outgoing arc of vertex `v`.
pub(crate) fn max_ovlp_src(
    options: &UlGraphOptions<'_>,
    reads: &crate::reads::read_store::ReadStore,
    v: u32,
) -> i64 {
    let mut longest = 0i64;
    let x = (v >> 1) as usize;
    let source = &options.sources[x];
    for hit in source.hits.iter().take(source.hits.len()) {
        let (result, arc) = ma_hit_to_arc(
            hit,
            reads.read_length(hit.query_identifier as usize) as i32,
            reads.read_length(hit.target_identifier as usize) as i32,
            options.maximum_hang,
            options.maximum_hang_rate,
            options.minimum_overlap,
        );
        if result < 0 {
            continue;
        }
        if (arc.ul >> 32) != u64::from(v) {
            continue;
        }
        if longest < i64::from(arc.overlap) {
            longest = i64::from(arc.overlap);
        }
    }
    longest
}

/// `max_ovlp_src_contain` (`inter.cpp` 7859): like [`max_ovlp_src`], but a
/// contained target counts with its flank length in the vertex direction.
pub(crate) fn max_ovlp_src_contain(
    options: &UlGraphOptions<'_>,
    reads: &crate::reads::read_store::ReadStore,
    v: u32,
) -> i64 {
    let mut longest: i64 = 0;
    let x = (v >> 1) as usize;
    let source = &options.sources[x];
    for hit in source.hits.iter().take(source.hits.len()) {
        let qn = hit.query_identifier;
        let tn = hit.target_identifier;
        let (result, arc) = ma_hit_to_arc(
            hit,
            reads.read_length(qn as usize) as i32,
            reads.read_length(tn as usize) as i32,
            options.maximum_hang,
            options.maximum_hang_rate,
            options.minimum_overlap,
        );
        if result >= 0 {
            if (arc.ul >> 32) != u64::from(v) {
                continue;
            }
            if longest < i64::from(arc.overlap) {
                longest = i64::from(arc.overlap);
            }
        } else if result == MA_HT_TCONT {
            // `tn` is contained in `qn`
            let flank = if v & 1 != 0 {
                i64::from(hit.query_end)
            } else {
                reads.read_length(qn as usize) as i64 - i64::from(hit.query_start)
            };
            if longest < flank {
                longest = flank;
            }
        }
    }
    longest
}

/// `find_ul_ov_max` (`inter.cpp` 2490): the last entry whose `qe` is below
/// `x`.
pub(crate) fn find_ul_ov_max(count: i32, entries: &[UlOv], x: u32) -> i32 {
    let n = count;
    if n == 0 {
        return -1;
    }
    if entries[(n - 1) as usize].qe < x {
        return n - 1;
    }
    if entries[0].qe >= x {
        return -1;
    }
    let mut start = 0i32;
    let mut end = n;
    while end > start {
        let middle = start + (end - start) / 2;
        if entries[middle as usize].qe >= x {
            end = middle;
        } else {
            start = middle + 1;
        }
    }
    debug_assert!(start == end);
    start
}

/// `connect_detect` (`inter.cpp` 7766): whether `a[ai]` may follow `a[aj]` in
/// a chain, caching the reachability test in `rch` and threading through the
/// contained reads between them. Returns the DP score or `INT32_MIN`.
///
/// The reference's `ridx` parameter is unused on this path (`infer_rovlp`
/// reads lengths from the unitig graph).
#[allow(clippy::too_many_arguments)]
pub(crate) fn connect_detect(
    entries: &mut [UlOv],
    count: i64,
    ai: i64,
    aj: i64,
    rch: &mut [i32],
    options: &UlGraphOptions<'_>,
    reads: &crate::reads::read_store::ReadStore,
    unitigs: &[MaUtg],
    band_width: i64,
    diff_ec_ul: f64,
    predecessor: &[i64],
    score: &[i32],
    trace_iterator: &mut RtraceIter,
    trace: &[Rtrace],
    targets: &TargetSource,
    query: &[u8],
    regions: &[crate::overlap::region::OverlapRegion],
    extz: &mut BitExtz,
    error_rate: f64,
    transition_score: i64,
    is_fc: &mut u32,
) -> i64 {
    let li = entries[ai as usize].clone();
    let lj = entries[aj as usize].clone();
    let li_v = (li.tn << 1) | u32::from(li.rev);
    let lj_v = (lj.tn << 1) | u32::from(lj.rev);
    *is_fc = 0;
    if li_v == lj_v {
        return i32::MIN as i64;
    }
    // even this pair has a overlap, its length will be very small; just ignore
    if lj.qe <= li.qs {
        return i32::MIN as i64;
    }

    let mut is_c;
    if (rch[aj as usize] >> 2) as i64 != ai {
        let query_overlap = infer_rovlp_reads(&li, &lj, &reads.read_lengths);
        rch[aj as usize] = ((ai as i32) << 2) | 3;
        let mut contained = false;
        if get_ecov_contain_adv(
            options,
            reads,
            li_v ^ 1,
            lj_v ^ 1,
            band_width,
            diff_ec_ul,
            query_overlap,
            Some(&mut contained),
        ) {
            rch[aj as usize] = ((ai as i32) << 2) | (contained as i32);
            is_c = contained as i64;
            if is_c != 0 {
                entries[aj as usize].el = false;
            }
        } else {
            return i32::MIN as i64;
        }
    }
    if !rch_connect(rch[aj as usize], ai) {
        return i32::MIN as i64;
    }
    is_c = i64::from(rch[aj as usize] & 1);
    let mut ak = aj;
    let mut afk = aj;
    let mut lowest_qs = i32::MAX;
    if is_c != 0 {
        lowest_qs = entries[aj as usize].qs as i32;
        ak = predecessor[aj as usize];
        while ak >= 0 {
            if (rch[ak as usize] >> 2) as i64 != ai {
                rch[ak as usize] = ((ak as i32) << 2) | 3;
                let lk = &entries[ak as usize];
                let lk_v = (lk.tn << 1) | u32::from(lk.rev);
                if li_v == lk_v {
                    break;
                }
                if lk.qe <= li.qs {
                    break;
                }
                let query_overlap = infer_rovlp_reads(&li, &lk, &reads.read_lengths);
                let mut contained = false;
                if !get_ecov_contain_adv(
                    options,
                    reads,
                    li_v ^ 1,
                    lk_v ^ 1,
                    band_width,
                    diff_ec_ul,
                    query_overlap,
                    Some(&mut contained),
                ) {
                    break;
                }
                rch[ak as usize] = ((ai as i32) << 2) | (contained as i32);
            }
            if !rch_connect(rch[ak as usize], ai) {
                break;
            }
            if lowest_qs > entries[ak as usize].qs as i32 {
                lowest_qs = entries[ak as usize].qs as i32;
            }
            is_c = i64::from(rch[ak as usize] & 1);
            if is_c == 0 {
                break;
            }
            ak = predecessor[ak as usize];
        }

        afk = ak;
        if ak >= 0 && !rch_connect(rch[ak as usize], ai) {
            return i32::MIN as i64; // go to a disconnected node
        }
    }

    // `afk >= 0`: reach to one non-contained read, `lqs` should be `<= li->qs`
    // `afk < 0`: all reads are contained, `lqs` is `>= li->qs`
    let mut score_value;
    if afk >= 0 || lowest_qs as i64 <= li.qs as i64 {
        let lj = &entries[aj as usize];
        let error = get_rid_backward_cigar_err(
            trace_iterator,
            &li,
            Some(trace),
            None,
            targets,
            query,
            Some(regions),
            None,
            extz,
            error_rate,
            i64::from(lj.qe),
        );
        score_value = i64::from(score[aj as usize]) + (li.qe as i64 - lj.qe as i64)
            - error * transition_score;
        *is_fc = 0;
        if afk < 0 {
            let mut alternative = li.qe as i64 - li.qs as i64;
            alternative -= (li.sec as i64) * transition_score;
            if score_value < alternative {
                score_value = alternative;
            }
        }
    } else {
        let mut value = li.qe as i64 - li.qs as i64;
        value -= (li.sec as i64) * transition_score;
        score_value = value;
        *is_fc = 1;
    }
    score_value
}

/// `quick_rechain_sc` (`inter.cpp` 7889): the score of one direct link for
/// [`quick_rechain`].
#[allow(clippy::too_many_arguments)]
pub(crate) fn quick_rechain_sc(
    entries: &[UlOv],
    regions: &[crate::overlap::region::OverlapRegion],
    options: &UlGraphOptions<'_>,
    reads: &crate::reads::read_store::ReadStore,
    unitigs: &[MaUtg],
    band_width: i64,
    diff_ec_ul: f64,
    trace: &[Rtrace],
    targets: &TargetSource,
    query: &[u8],
    trace_iterator: &mut RtraceIter,
    error_rate: f64,
    extz: &mut BitExtz,
    transition_score: i64,
    score: &[i32],
    ai: i64,
    aj: i64,
) -> i64 {
    let li = &entries[ai as usize];
    let lj = &entries[aj as usize];
    let li_v = (li.tn << 1) | u32::from(li.rev);
    let lj_v = (lj.tn << 1) | u32::from(lj.rev);
    if li_v == lj_v {
        return i32::MIN as i64;
    }
    if lj.qe <= li.qs {
        return i32::MIN as i64;
    }
    if lj.qs > li.qs {
        return i32::MIN as i64;
    }
    let query_overlap = infer_rovlp_reads(&li, &lj, &reads.read_lengths);
    let mut contained = false;
    if get_ecov_contain_adv(
        options,
        reads,
        li_v ^ 1,
        lj_v ^ 1,
        band_width,
        diff_ec_ul,
        query_overlap,
        Some(&mut contained),
    ) && !contained
    {
        let error = get_rid_backward_cigar_err(
            trace_iterator,
            &li,
            Some(trace),
            None,
            targets,
            query,
            Some(regions),
            None,
            extz,
            error_rate,
            i64::from(lj.qe),
        );
        return i64::from(score[aj as usize]) + (li.qe as i64 - lj.qe as i64)
            - error * transition_score;
    }
    i32::MIN as i64
}

/// `if_qchain_cnn` (`inter.cpp` 8476): whether the suffix entry `li` may
/// connect to the prefix entry `lj`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn if_qchain_cnn(
    options: &UlGraphOptions<'_>,
    reads: &crate::reads::read_store::ReadStore,
    unitigs: &[MaUtg],
    band_width: i64,
    diff_ec_ul: f64,
    li: &UlOv,
    lj: &UlOv,
    is_c: &mut i64,
) -> bool {
    let li_v = (li.tn << 1) | u32::from(li.rev);
    let lj_v = (lj.tn << 1) | u32::from(lj.rev);
    if li_v == lj_v || lj.qe <= li.qs {
        return false;
    }
    let query_overlap = infer_rovlp_reads(&li, &lj, &reads.read_lengths);
    let mut contained = false;
    if get_ecov_contain_adv(
        options,
        reads,
        li_v ^ 1,
        lj_v ^ 1,
        band_width,
        diff_ec_ul,
        query_overlap,
        Some(&mut contained),
    ) {
        *is_c = contained as i64;
        return true;
    }
    false
}

/// `quick_rechain_propagate_transitive_reduction` (`inter.cpp` 7915).
#[allow(clippy::too_many_arguments)]
pub(crate) fn quick_rechain_propagate_transitive_reduction(
    options: &UlGraphOptions<'_>,
    reads: &crate::reads::read_store::ReadStore,
    unitigs: &[MaUtg],
    band_width: i64,
    diff_ec_ul: f64,
    entries: &[UlOv],
    entry_count: i32,
    ai: i64,
    rch: &mut [i32],
    score: &mut [i32],
    predecessor: &[i64],
    occurrence: &mut [i32],
    scratch: &mut [i64],
    best_score: &mut i64,
    best_index: &mut i64,
    best_count: &mut i64,
) {
    if *best_index < 0 {
        return;
    }
    let best_index0 = *best_index;
    for j in best_index0 + 1..entry_count as i64 {
        scratch[j as usize] = best_index0 - 1;
        if predecessor[j as usize] < 0 {
            continue;
        }
        if (rch[j as usize] >> 2) as i64 != ai {
            rch[j as usize] = ((ai as i32) << 2) | 3;
            let mut is_c = 0i64;
            if entries[j as usize].el
                && if_qchain_cnn(
                    options,
                    reads,
                    unitigs,
                    band_width,
                    diff_ec_ul,
                    &entries[ai as usize],
                    &entries[j as usize],
                    &mut is_c,
                )
            {
                rch[j as usize] = ((ai as i32) << 2) | (is_c as i32);
            }
        }
        if !rch_connect(rch[j as usize], ai) {
            continue;
        }
        let mut k = predecessor[j as usize];
        while k >= 0 && k > best_index0 {
            if scratch[k as usize] == best_index0 {
                k = best_index0;
            } else {
                k = best_index0 - 1;
            }
            break;
        }
        if k != best_index0 {
            continue;
        }
        scratch[j as usize] = best_index0; // a[j] could reach mm_idx0
        if predecessor[j as usize] != k {
            let mut is_c = 0i64;
            if !if_qchain_cnn(
                options,
                reads,
                unitigs,
                band_width,
                diff_ec_ul,
                &entries[j as usize],
                &entries[k as usize],
                &mut is_c,
            ) {
                continue;
            }
        }
        let count = occurrence[j as usize] as i64 + 1;
        if count >= *best_count {
            // must >=
            *best_count = count;
            *best_index = j;
        }
    }
    let _ = score;
}

/// `quick_rechain` (`inter.cpp` 7983): the lightweight re-chain of the
/// contained prefix, returning the traced chain length in `scratch`.
///
/// The scratch arrays alias `Chain_Data` in the reference; the port slices
/// them from the caller's `ChainData`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn quick_rechain(
    regions: &[crate::overlap::region::OverlapRegion],
    entries: &mut [UlOv],
    entry_count: i64,
    trace: &[Rtrace],
    options: &UlGraphOptions<'_>,
    reads: &crate::reads::read_store::ReadStore,
    unitigs: &[MaUtg],
    band_width: i64,
    diff_ec_ul: f64,
    query_length: i64,
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
    extz: &mut BitExtz,
    transition_score: i64,
    targets: &TargetSource,
    query: &[u8],
    error_rate: f64,
    data: &mut ChainData,
) -> i64 {
    let score = &mut data.score;
    let predecessor = &mut data.predecessor;
    let occurrence = &mut data.occurrence;
    let self_length = &mut data.self_length;
    let scratch = &mut data.scratch;
    for value in scratch[..entry_count as usize].iter_mut() {
        *value = 0;
    }

    let mut st = 0i64;
    let mut max_ii = -1i64;
    for i in 0..entry_count {
        let li = entries[i as usize].clone();
        let li_v = (li.tn << 1) | u32::from(li.rev);
        self_length[i as usize] = i32::MAX;
        let mut current_score = li.qe as i64 - li.qs as i64;
        current_score -= (li.sec as i64) * transition_score;
        let mut best_score = i32::MIN as i64;
        let mut best_index = -1i64;
        let mut best_count = 0i64;
        let mut n_skip = 0i64;
        let mut end_j;
        let mut trace_iterator = RtraceIter::default();
        trace_iterator.k = i32::MAX as i64;

        if li.el {
            let maximum_overlap = max_ovlp_src_contain(options, reads, li_v ^ 1);
            let mut x = ((li.qs as i64 + maximum_overlap) as f64 * diff_ec_ul) as i64;
            if x < band_width {
                x = band_width;
            }
            x += li.qs as i64 + maximum_overlap;
            if x > query_length + 1 {
                x = query_length + 1;
            }
            let mut x = find_ul_ov_max(i as i32, entries, (x as u32).wrapping_add(G_CHAIN_INDEL)) as i64;
            if (x - st) > max_iter {
                st = x - max_iter;
            }
            let search_limit = x;
            let mut j = search_limit;
            while j >= st {
                // collect potential destination vertices
                let lj = &entries[j as usize];
                self_length[j as usize] = ((i as i32) << 2) | 3; // set it to unavailable
                if !lj.el {
                    j -= 1;
                    continue;
                }
                if lj.qe <= li.qs {
                    break;
                }
                if lj.qs > li.qs {
                    // no contain
                    j -= 1;
                    continue;
                }
                let sc = quick_rechain_sc(
                    entries,
                    regions,
                    options,
                    reads,
                    unitigs,
                    band_width,
                    diff_ec_ul,
                    trace,
                    targets,
                    query,
                    &mut trace_iterator,
                    error_rate,
                    extz,
                    transition_score,
                    score,
                    i,
                    j,
                );
                if sc == i32::MIN as i64 {
                    j -= 1;
                    continue;
                }
                self_length[j as usize] = (i as i32) << 2; // set it to available
                let count = occurrence[j as usize] as i64 + 1;
                if sc > best_score || (sc == best_score && count > best_count) {
                    best_score = sc;
                    best_index = j;
                    best_count = count;
                    if n_skip > 0 {
                        n_skip -= 1;
                    }
                } else if scratch[j as usize] == i {
                    n_skip += 1;
                    if n_skip > max_skip {
                        break;
                    }
                }
                if predecessor[j as usize] >= 0 {
                    scratch[predecessor[j as usize] as usize] = i;
                }
                j -= 1;
            }

            end_j = j;
            if max_ii < 0 || entries[i as usize].qe as i64 > entries[max_ii as usize].qe as i64 + max_dis {
                // too long
                let mut max = i32::MIN as i64;
                let mut max_n = i32::MIN as i64;
                max_ii = -1;
                let mut j = i - 1;
                while j >= st && (entries[i as usize].qe as i64) <= max_dis + entries[j as usize].qe as i64 {
                    if max < score[j as usize] as i64
                        || (max == score[j as usize] as i64 && max_n < occurrence[j as usize] as i64)
                    {
                        max = score[j as usize] as i64;
                        max_n = occurrence[j as usize] as i64;
                        max_ii = j;
                    }
                    j -= 1;
                }
            }

            if max_ii >= 0 && max_ii < end_j {
                // just have a try with a[i] <-> a[max_ii]
                let lj = entries[max_ii as usize].clone();
                if lj.el && lj.qe > li.qs && lj.qs <= li.qs {
                    let sc = quick_rechain_sc(
                        entries,
                        regions,
                        options,
                        reads,
                        unitigs,
                        band_width,
                        diff_ec_ul,
                        trace,
                        targets,
                        query,
                        &mut trace_iterator,
                        error_rate,
                        extz,
                        transition_score,
                        score,
                        i,
                        max_ii,
                    );
                    // as max_ii < end_j, get_rid_backward_cigar_err still works
                    if sc != i32::MIN as i64 {
                        let count = occurrence[max_ii as usize] as i64 + 1;
                        if sc > best_score || (sc == best_score && count > best_count) {
                            best_score = sc;
                            best_index = max_ii;
                            best_count = count;
                        }
                    }
                }
            }
            let sc = current_score;
            let sn = 1i64; // for i itself
            if sc > best_score || (sc == best_score && sn > best_count) {
                best_score = sc;
                best_index = -1;
                best_count = sn;
            }

            if best_score < 0 {
                best_score = current_score;
                best_index = -1;
                best_count = 1;
            }

            if best_index >= 0 {
                // no need to update mm_sc and is_fc
                let (mut best_score_ref, mut best_index_ref, mut best_count_ref) =
                    (best_score, best_index, best_count);
                let pred_copy = predecessor.clone();
                quick_rechain_propagate_transitive_reduction(
                    options,
                    reads,
                    unitigs,
                    band_width,
                    diff_ec_ul,
                    entries,
                    (search_limit + 1) as i32,
                    i,
                    self_length,
                    score,
                    &pred_copy,
                    occurrence,
                    scratch,
                    &mut best_score_ref,
                    &mut best_index_ref,
                    &mut best_count_ref,
                );
                best_score = best_score_ref;
                best_index = best_index_ref;
                best_count = best_count_ref;
            }
        }
        score[i as usize] = best_score as i32;
        predecessor[i as usize] = best_index;
        occurrence[i as usize] = best_count as i32;
        if max_ii < 0
            || ((entries[i as usize].qe as i64) <= max_dis + entries[max_ii as usize].qe as i64
                && (score[max_ii as usize] as i64) < best_score)
        {
            max_ii = i;
        }
    }

    let mut chain_length = 0i64;
    let mut i = entry_count - 1;
    while i >= 0 {
        scratch[chain_length as usize] = i;
        chain_length += 1;
        i = predecessor[i as usize];
    }
    chain_length
}

/// `push_trans_aln` (`inter.cpp` 8097): collects the contained reads between
/// two chain entries that corroborate both flanks, marking them as used.
///
/// `ai` is smaller than `aj`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn push_trans_aln(
    options: &UlGraphOptions<'_>,
    reads: &crate::reads::read_store::ReadStore,
    unitigs: &[MaUtg],
    band_width: i64,
    diff_ec_ul: f64,
    entries: &[UlOv],
    mut ai: i64,
    aj: i64,
    vis: &mut [u64],
    occurrence: &mut [i32],
) -> i64 {
    let mut m = 0i64;
    for k in ai + 1..aj {
        if vis[k as usize] & 1 != 0 {
            continue;
        }
        if !entries[k as usize].el {
            continue;
        }
        let lj = &entries[ai as usize];
        let li = &entries[k as usize];
        let li_v = (li.tn << 1) | u32::from(li.rev);
        let lj_v = (lj.tn << 1) | u32::from(lj.rev);
        if lj.qe <= li.qs {
            continue;
        }
        if lj.qs > li.qs {
            continue;
        }
        if li_v == lj_v {
            continue;
        }
        let query_overlap = infer_rovlp_reads(&li, &lj, &reads.read_lengths);
        let mut contained = false;
        if !get_ecov_contain_adv(
            options,
            reads,
            li_v ^ 1,
            lj_v ^ 1,
            band_width,
            diff_ec_ul,
            query_overlap,
            Some(&mut contained),
        ) {
            continue;
        }
        if contained {
            continue;
        }

        let lj = &entries[k as usize];
        let li = &entries[aj as usize];
        let li_v = (li.tn << 1) | u32::from(li.rev);
        let lj_v = (lj.tn << 1) | u32::from(lj.rev);
        if lj.qe <= li.qs {
            continue;
        }
        if lj.qs > li.qs {
            continue;
        }
        if li_v == lj_v {
            continue;
        }
        let query_overlap = infer_rovlp_reads(&li, &lj, &reads.read_lengths);
        let mut contained = false;
        if !get_ecov_contain_adv(
            options,
            reads,
            li_v ^ 1,
            lj_v ^ 1,
            band_width,
            diff_ec_ul,
            query_overlap,
            Some(&mut contained),
        ) {
            continue;
        }
        if contained {
            continue;
        }

        occurrence[m as usize] = k as i32;
        m += 1;
        vis[k as usize] |= 1;
        ai = k;
    }
    m
}

/// `flat_contain` (`inter.cpp` 8374): the plain DP over one chain body.
#[allow(clippy::too_many_arguments)]
pub(crate) fn flat_contain(
    reads: &crate::reads::read_store::ReadStore,
    options: &UlGraphOptions<'_>,
    unitigs: &[MaUtg],
    band_width: i64,
    diff_ec_ul: f64,
    query_length: i64,
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
    entries: &mut [UlOv],
    entry_count: i64,
    scratch: &mut [i32],
    occurrence: &mut [i32],
    predecessor: &mut [i64],
    score: &mut [i64],
    bounds: &mut UlOv,
) -> i64 {
    if entry_count <= 0 {
        return 0;
    }
    // group equal `qe` runs and sort them by `qs`
    {
        let mut j = 0usize;
        for i in 1..=entry_count as usize {
            if i == entry_count as usize || entries[i].qe != entries[j].qe {
                if i - j > 1 {
                    crate::container::ksort::radix_sort_by_key_32(
                        &mut entries[j..i],
                        &|entry: &UlOv| entry.qs,
                    );
                }
                j = i;
            }
        }
    }
    for value in scratch[..entry_count as usize].iter_mut() {
        *value = 0;
    }

    let mut st = 0i64;
    let mut max_ii = -1i64;
    let mut total_score = i32::MIN as i64;
    let mut total_count = i32::MIN as i64;
    let mut total_index = -1i64;
    for i in 0..entry_count {
        let li = entries[i as usize].clone();
        let li_v = (li.tn << 1) | u32::from(li.rev);
        let maximum_overlap = max_ovlp_src(options, reads, li_v ^ 1);
        let mut x = ((li.qs as i64 + maximum_overlap) as f64 * diff_ec_ul) as i64;
        if x < band_width {
            x = band_width;
        }
        x += li.qs as i64 + maximum_overlap;
        if x > query_length + 1 {
            x = query_length + 1;
        }
        let x = find_ul_ov_max(i as i32, entries, (x as u32).wrapping_add(G_CHAIN_INDEL)) as i64;
        let current_score = if li.el { li.qe as i64 - li.qs as i64 } else { 0 };
        let mut best_score = current_score;
        let mut best_index = -1i64;
        let mut best_count = 1i64;
        let mut n_skip = 0i64;
        let end_j;
        if (x - st) > max_iter {
            st = x - max_iter;
        }
        let mut j = x;
        while j >= st {
            // collect potential destination vertices
            let lj = &entries[j as usize];
            if lj.qe <= li.qs {
                break;
            }
            let query_overlap = infer_rovlp_reads(&li, lj, &reads.read_lengths);
            let lj_v = (lj.tn << 1) | u32::from(lj.rev);
            let mut contained = false;
            if li_v != lj_v
                && get_ecov_contain_adv(
                    options,
                    reads,
                    li_v ^ 1,
                    lj_v ^ 1,
                    band_width,
                    diff_ec_ul,
                    query_overlap,
                    Some(&mut contained),
                )
                && !contained
            {
                let sc = score[j as usize] + current_score;
                let sn = occurrence[j as usize] as i64 + 1;
                if sc > best_score || (sc == best_score && sn > best_count) {
                    best_score = sc;
                    best_index = j;
                    best_count = sn;
                    if n_skip > 0 {
                        n_skip -= 1;
                    }
                } else if scratch[j as usize] as i64 == i {
                    n_skip += 1;
                    if n_skip > max_skip {
                        break;
                    }
                }
                if predecessor[j as usize] >= 0 {
                    scratch[predecessor[j as usize] as usize] = i as i32;
                }
            }
            j -= 1;
        }

        end_j = j;
        if max_ii < 0 || entries[i as usize].qe as i64 > entries[max_ii as usize].qe as i64 + max_dis {
            // too long
            let mut max = i32::MIN as i64;
            let mut max_n = i32::MIN as i64;
            max_ii = -1;
            let mut j = i - 1;
            while j >= st && (entries[i as usize].qe as i64) <= max_dis + entries[j as usize].qe as i64 {
                if max < score[j as usize] || (max == score[j as usize] && max_n < occurrence[j as usize] as i64)
                {
                    max = score[j as usize];
                    max_n = occurrence[j as usize] as i64;
                    max_ii = j;
                }
                j -= 1;
            }
        }

        if max_ii >= 0 && max_ii < end_j {
            // just have a try with a[i] <-> a[max_ii]
            let lj = entries[max_ii as usize].clone();
            if lj.qe > li.qs && lj.qs < li.qs {
                let query_overlap = infer_rovlp_reads(&li, &lj, &reads.read_lengths);
                let lj_v = (lj.tn << 1) | u32::from(lj.rev);
                let mut contained = false;
                if li_v != lj_v
                    && get_ecov_contain_adv(
                        options,
                        reads,
                        li_v ^ 1,
                        lj_v ^ 1,
                        band_width,
                        diff_ec_ul,
                        query_overlap,
                        Some(&mut contained),
                    )
                    && !contained
                {
                    let sc = score[max_ii as usize] + current_score;
                    let sn = occurrence[max_ii as usize] as i64 + 1;
                    if sc > best_score || (sc == best_score && sn > best_count) {
                        best_score = sc;
                        best_index = max_ii;
                        best_count = sn;
                    }
                }
            }
        }

        score[i as usize] = best_score;
        predecessor[i as usize] = best_index;
        occurrence[i as usize] = best_count as i32;
        if max_ii < 0
            || ((entries[i as usize].qe as i64) <= max_dis + entries[max_ii as usize].qe as i64
                && score[max_ii as usize] < best_score)
        {
            max_ii = i;
        }

        if total_score < best_score || (total_score == best_score && total_count < best_count) {
            total_score = best_score;
            total_count = best_count;
            total_index = i;
        }
    }

    let mut chain_length = 0i64;
    let mut i = total_index;
    while i >= 0 {
        scratch[chain_length as usize] = i as i32;
        chain_length += 1;
        i = predecessor[i as usize];
    }
    bounds.qs = entries[scratch[(chain_length - 1) as usize] as usize].qs;
    bounds.qe = entries[scratch[(chain_length - 1) as usize] as usize].qe;
    for i in 0..chain_length as usize {
        entries[i] = entries[scratch[(chain_length - 1 - i as i64) as usize] as usize];
        if entries[i].qs < bounds.qs {
            bounds.qs = entries[i].qs;
        }
        if entries[i].qe > bounds.qe {
            bounds.qe = entries[i].qe;
        }
    }
    chain_length
}

/// `propagate_transitive_reduction` (`inter.cpp` 8485): the reachability
/// propagation of the raw-chain DP, without the contained-only restriction of
/// [`quick_rechain_propagate_transitive_reduction`].
#[allow(clippy::too_many_arguments)]
pub(crate) fn propagate_transitive_reduction(
    options: &UlGraphOptions<'_>,
    reads: &crate::reads::read_store::ReadStore,
    unitigs: &[MaUtg],
    band_width: i64,
    diff_ec_ul: f64,
    entries: &[UlOv],
    entry_count: i32,
    ai: i64,
    rch: &mut [i32],
    predecessor: &[i64],
    occurrence: &mut [i32],
    scratch: &mut [i64],
    best_score: &mut i64,
    best_index: &mut i64,
    best_count: &mut i64,
) {
    if *best_index < 0 {
        return;
    }
    let best_index0 = *best_index;
    for j in best_index0 + 1..entry_count as i64 {
        scratch[j as usize] = best_index0 - 1;
        if predecessor[j as usize] < 0 {
            continue;
        }
        if (rch[j as usize] >> 2) as i64 != ai {
            rch[j as usize] = ((ai as i32) << 2) | 3;
            let mut is_c = 0i64;
            if if_qchain_cnn(
                options,
                reads,
                unitigs,
                band_width,
                diff_ec_ul,
                &entries[ai as usize],
                &entries[j as usize],
                &mut is_c,
            ) {
                rch[j as usize] = ((ai as i32) << 2) | (is_c as i32);
            }
        }
        if !rch_connect(rch[j as usize], ai) {
            continue;
        }
        let mut k = predecessor[j as usize];
        while k >= 0 && k > best_index0 {
            if scratch[k as usize] == best_index0 {
                k = best_index0;
            } else {
                k = best_index0 - 1;
            }
            break;
        }
        if k != best_index0 {
            continue;
        }
        scratch[j as usize] = best_index0; // a[j] could reach mm_idx0
        if predecessor[j as usize] != k {
            let mut is_c = 0i64;
            if !if_qchain_cnn(
                options,
                reads,
                unitigs,
                band_width,
                diff_ec_ul,
                &entries[j as usize],
                &entries[k as usize],
                &mut is_c,
            ) {
                continue;
            }
        }
        let count = occurrence[j as usize] as i64 + 1;
        if count >= *best_count {
            // must >=
            *best_count = count;
            *best_index = j;
        }
    }
}

/// `flat_contain_adv` (`inter.cpp` 8170): the DP over one chain body that
/// also threads the transitive alignments and retries the contained prefix
/// through [`quick_rechain`].
///
/// The reference reuses the caller's scratch arrays and the `id` array as the
/// `flag` store; the port takes fresh vectors, which is content-equivalent
/// because every field is written before it is read within one call. The
/// reference's `ol` parameter is unused on this path.
#[allow(clippy::too_many_arguments)]
pub(crate) fn flat_contain_adv(
    reads: &crate::reads::read_store::ReadStore,
    options: &UlGraphOptions<'_>,
    unitigs: &[MaUtg],
    band_width: i64,
    diff_ec_ul: f64,
    query_length: i64,
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
    extz: &mut BitExtz,
    transition_score: i64,
    targets: &TargetSource,
    query: &[u8],
    error_rate: f64,
    trace: &[Rtrace],
    entries: &mut [UlOv],
    entry_count: i64,
    regions: &[crate::overlap::region::OverlapRegion],
    original: &[UlOv],
    flag: &mut [u64],
    rechain_window: u64,
    vis: &mut [u64],
    bounds: &mut UlOv,
    chain_data: &mut ChainData,
) -> i64 {
    if entry_count <= 0 {
        return 0;
    }
    // group equal `qe` runs and sort them by `qs`
    {
        let mut j = 0usize;
        for i in 1..=entry_count as usize {
            if i == entry_count as usize || entries[i].qe != entries[j].qe {
                if i - j > 1 {
                    crate::container::ksort::radix_sort_by_key_32(
                        &mut entries[j..i],
                        &|entry: &UlOv| entry.qs,
                    );
                }
                j = i;
            }
        }
    }

    let mut low_priority_count = 0i64;
    let mut chain_length = 0i64;
    let mut quick_chain_length = 0i64;
    let mut total_score = i32::MIN as i64;
    let mut total_count = i32::MIN as i64;
    let mut total_index = -1i64;
    {
        let mut scratch_t: Vec<i32> = vec![0; entry_count as usize];
        let mut dp_f: Vec<i64> = vec![0; entry_count as usize];
        let mut dp_p: Vec<i64> = vec![-1; entry_count as usize];
        let mut dp_c_n: Vec<i32> = vec![0; entry_count as usize];
        let mut st = 0i64;
        let mut max_ii = -1i64;
        for i in 0..entry_count {
            let li = entries[i as usize].clone();
            let li_v = (li.tn << 1) | u32::from(li.rev);
            flag[i as usize] = u64::from(li.sec);
            entries[i as usize].sec = original[(li.sec >> 1) as usize].sec;

            let current_score = if li.el { li.qe as i64 - li.qs as i64 } else { 0 };
            let mut best_score = current_score;
            let mut best_count = 1i64;
            let mut best_index = -1i64;
            let mut n_skip = 0i64;
            if flag[i as usize] & 1 != 0 {
                // low priority
                low_priority_count += 1;
                dp_f[i as usize] = i32::MIN as i64;
                dp_p[i as usize] = best_index;
                dp_c_n[i as usize] = 0;
                continue;
            }
            let maximum_overlap = max_ovlp_src(options, reads, li_v ^ 1);
            let mut x = ((li.qs as i64 + maximum_overlap) as f64 * diff_ec_ul) as i64;
            if x < band_width {
                x = band_width;
            }
            x += li.qs as i64 + maximum_overlap;
            if x > query_length + 1 {
                x = query_length + 1;
            }
            let x = find_ul_ov_max(i as i32, entries, (x as u32).wrapping_add(G_CHAIN_INDEL)) as i64;
            if (x - st) > max_iter {
                st = x - max_iter;
            }
            let mut j = x;
            while j >= st {
                // collect potential destination vertices
                let lj = &entries[j as usize];
                let lj_v = (lj.tn << 1) | u32::from(lj.rev);
                if flag[j as usize] & 1 != 0 {
                    // low priority
                    j -= 1;
                    continue;
                }
                if lj.qe <= li.qs {
                    break;
                }
                if lj.qs > li.qs {
                    // no contain, must be > instead of >=
                    j -= 1;
                    continue;
                }
                let query_overlap = infer_rovlp_reads(&li, lj, &reads.read_lengths);
                let mut contained = false;
                if li_v != lj_v
                    && get_ecov_contain_adv(
                        options,
                        reads,
                        li_v ^ 1,
                        lj_v ^ 1,
                        band_width,
                        diff_ec_ul,
                        query_overlap,
                        Some(&mut contained),
                    )
                    && !contained
                {
                    let sc = dp_f[j as usize] + current_score;
                    let sn = dp_c_n[j as usize] as i64 + 1;
                    if sc > best_score || (sc == best_score && sn > best_count) {
                        best_score = sc;
                        best_index = j;
                        best_count = sn;
                        if n_skip > 0 {
                            n_skip -= 1;
                        }
                    } else if scratch_t[j as usize] as i64 == i {
                        n_skip += 1;
                        if n_skip > max_skip {
                            break;
                        }
                    }
                    if dp_p[j as usize] >= 0 {
                        scratch_t[dp_p[j as usize] as usize] = i as i32;
                    }
                }
                j -= 1;
            }

            let end_j = j;
            if max_ii < 0
                || entries[i as usize].qe as i64 > entries[max_ii as usize].qe as i64 + max_dis
            {
                // too long
                let mut max = i32::MIN as i64;
                let mut max_n = i32::MIN as i64;
                max_ii = -1;
                let mut j = i - 1;
                while j >= st
                    && (entries[i as usize].qe as i64) <= max_dis + entries[j as usize].qe as i64
                {
                    if max < dp_f[j as usize]
                        || (max == dp_f[j as usize] && max_n < dp_c_n[j as usize] as i64)
                    {
                        max = dp_f[j as usize];
                        max_n = dp_c_n[j as usize] as i64;
                        max_ii = j;
                    }
                    j -= 1;
                }
            }

            if max_ii >= 0 && max_ii < end_j {
                // just have a try with a[i] <-> a[max_ii]
                let lj = entries[max_ii as usize].clone();
                if lj.qe > li.qs && lj.qs <= li.qs && flag[max_ii as usize] & 1 == 0 {
                    let query_overlap = infer_rovlp_reads(&li, &lj, &reads.read_lengths);
                    let lj_v = (lj.tn << 1) | u32::from(lj.rev);
                    let mut contained = false;
                    if li_v != lj_v
                        && get_ecov_contain_adv(
                            options,
                            reads,
                            li_v ^ 1,
                            lj_v ^ 1,
                            band_width,
                            diff_ec_ul,
                            query_overlap,
                            Some(&mut contained),
                        )
                        && !contained
                    {
                        let sc = dp_f[max_ii as usize] + current_score;
                        let sn = dp_c_n[max_ii as usize] as i64 + 1;
                        if sc > best_score || (sc == best_score && sn > best_count) {
                            best_score = sc;
                            best_index = max_ii;
                            best_count = sn;
                        }
                    }
                }
            }

            dp_f[i as usize] = best_score;
            dp_p[i as usize] = best_index;
            dp_c_n[i as usize] = best_count as i32;
            if max_ii < 0
                || ((entries[i as usize].qe as i64) <= max_dis + entries[max_ii as usize].qe as i64
                    && dp_f[max_ii as usize] < best_score)
            {
                max_ii = i;
            }

            if total_score < best_score || (total_score == best_score && total_count < best_count) {
                total_score = best_score;
                total_count = best_count;
                total_index = i;
            }
        }

        // trace the best chain into `flag`'s high half
        let mut i = total_index;
        while i >= 0 {
            flag[chain_length as usize] |= (i as u64) << 32;
            i = dp_p[i as usize];
            chain_length += 1;
        }
    }

    let mut st = (flag[(chain_length - 1) as usize] >> 32) as usize;
    let qs = entries[st].qs as u64;
    if qs <= rechain_window && low_priority_count > 0 {
        let mut i = 0usize;
        let mut full_flanks = 0i64;
        while i < chain_length as usize {
            st = (flag[(chain_length - 1 - i as i64) as usize] >> 32) as usize;
            if entries[st].el {
                break;
            }
            if entries[st].ts == 0
                && entries[st].te as i64 == reads.read_length(entries[st].tn as usize) as i64
            {
                full_flanks += 1;
            }
            i += 1;
        }
        if (i as i64) < chain_length && entries[st].qs as u64 > qs && full_flanks > 0 {
            quick_chain_length = quick_rechain(
                regions,
                entries,
                st as i64 + 1,
                trace,
                options,
                reads,
                unitigs,
                band_width,
                diff_ec_ul,
                query_length,
                max_skip,
                max_iter,
                max_dis,
                extz,
                transition_score,
                targets,
                query,
                error_rate,
                chain_data,
            );
            let tail = &chain_data.scratch;
            if quick_chain_length > 0
                && entries[tail[(quick_chain_length - 1) as usize] as usize].qs as u64 <= qs
            {
                let mut cl = chain_length - i as i64 - 1;
                debug_assert!(st as i64 == tail[0]);
                for value in tail.iter().take(quick_chain_length as usize) {
                    let slot = flag[cl as usize];
                    flag[cl as usize] = (slot << 32) >> 32;
                    flag[cl as usize] |= (*value as u64) << 32;
                    cl += 1;
                }
                chain_length = cl;
            }
        }
    }

    let mut t0: Vec<i64> = Vec::with_capacity(chain_length as usize);
    for i in 0..chain_length as usize {
        t0.push((flag[(chain_length - 1 - i as i64) as usize] >> 32) as i64);
        let raw = flag[t0[i] as usize];
        vis[((raw as u32) >> 1) as usize] |= 1; // set as used
    }

    let mut occurrence: Vec<i32> = Vec::new();
    let mut i = 0usize;
    if quick_chain_length > 1 {
        while i < quick_chain_length as usize {
            if flag[t0[i] as usize] & 1 != 0 {
                flag[t0[i] as usize] -= 1;
            }
            let raw = flag[t0[i] as usize];
            vis[((raw as u32) >> 1) as usize] |= 1;
            occurrence.push(((raw as u32) >> 1) as i32);

            if i < quick_chain_length as usize - 1 {
                let ori_i = ((flag[t0[i] as usize] as u32) >> 1) as i64;
                let ori_j = ((flag[t0[i + 1] as usize] as u32) >> 1) as i64;
                let appended = push_trans_aln(
                    options,
                    reads,
                    unitigs,
                    band_width,
                    diff_ec_ul,
                    original,
                    ori_i,
                    ori_j,
                    vis,
                    &mut occurrence,
                );
                let _ = appended;
            }
            i += 1;
        }
    }

    while i < chain_length as usize {
        if flag[t0[i] as usize] & 1 != 0 {
            flag[t0[i] as usize] -= 1;
        }
        let raw = flag[t0[i] as usize];
        vis[((raw as u32) >> 1) as usize] |= 1;
        occurrence.push(((raw as u32) >> 1) as i32);
        i += 1;
    }
    let new_chain_length = occurrence.len() as i64;
    debug_assert!(new_chain_length >= chain_length);
    if new_chain_length <= 0 {
        return new_chain_length;
    }

    bounds.qs = original[occurrence[0] as usize].qs;
    bounds.qe = original[occurrence[0] as usize].qe;
    for i in 0..new_chain_length as usize {
        entries[i] = original[occurrence[i] as usize];
        if entries[i].qs < bounds.qs {
            bounds.qs = entries[i].qs;
        }
        if entries[i].qe > bounds.qe {
            bounds.qe = entries[i].qe;
        }
    }
    new_chain_length
}

/// `gen_trans_aln` (`inter.cpp` 8717): fills the gap between two chain
/// elements with the contained reads that corroborate both flanks, deletes
/// the new candidates contradicted by the chain elements, and flattens the
/// result with [`flat_contain_adv`].
///
/// The reference's `ol` parameter is unused on this path.
#[allow(clippy::too_many_arguments)]
pub(crate) fn gen_trans_aln(
    reads: &crate::reads::read_store::ReadStore,
    options: &UlGraphOptions<'_>,
    unitigs: &[MaUtg],
    band_width: i64,
    diff_ec_ul: f64,
    query_length: i64,
    chain: &[usize],
    entries: &[UlOv],
    entry_count: i64,
    chain_qs: u64,
    chain_qe: u64,
    vis: &mut [u64],
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
    extz: &mut BitExtz,
    transition_score: i64,
    targets: &TargetSource,
    query: &[u8],
    error_rate: f64,
    trace: &[Rtrace],
    rechain_window: u32,
    result: &mut Vec<UlOv>,
    regions: &[crate::overlap::region::OverlapRegion],
    chain_data: &mut ChainData,
    bounds: &mut UlOv,
) -> i64 {
    let idn = chain.len() as i64;
    if idn == 0 {
        return 0;
    }
    let mut nv = 0i64;
    let mut tipl = 0u32;
    const HIGH_MARK: u32 = 0x2000_0000;
    // the reference writes into `ex + n_v0`; keep this call's region separate
    let region_base = result.len();

    // the contained reads before the first chain element with full flanks
    let mut full_flanks = 0i64;
    let mut first = idn;
    for i in 0..idn {
        let entry = &entries[chain[i as usize]];
        if entry.el {
            first = i;
            break;
        }
        if entry.ts == 0 && entry.te as i64 == reads.read_length(entry.tn as usize) as i64 {
            full_flanks += 1;
        }
    }
    if (first as i64) < idn && full_flanks > 0 && entries[chain[first as usize]].qs as u64 > chain_qs && chain_qs <= rechain_window as u64 {
        tipl = chain_qs as u32; // there are contained chains beforehand
    }

    let mut previous = -1i64;
    for index in 0..idn {
        let k = chain[index as usize] as i64;
        if k - previous > 1 {
            let mut l = previous;
            if l < 0 {
                l = k - 1;
                while l >= 0 && entries[l as usize].qe as u64 > chain_qs {
                    l -= 1;
                }
            }

            for z in l + 1..k {
                let li = &entries[z as usize];
                let li_v = (li.tn << 1) | u32::from(li.rev);
                if !li.el || vis[z as usize] & 1 != 0 {
                    continue;
                }
                if !((i64::from(li.qs) + i64::from(tipl)) >= chain_qs as i64
                    && (li.qe as u64) < chain_qe)
                {
                    continue;
                }
                let mut update = true;
                let mut j = index - 1;
                while j >= 0 {
                    let lj = &entries[chain[j as usize]];
                    let lj_v = (lj.tn << 1) | u32::from(lj.rev);
                    if lj.qe <= li.qs {
                        break;
                    }
                    let query_overlap = infer_rovlp_reads(&li, &lj, &reads.read_lengths);
                    let mut contained = false;
                    if li_v == lj_v
                        || !get_ecov_contain_adv(
                            options,
                            reads,
                            li_v ^ 1,
                            lj_v ^ 1,
                            band_width,
                            diff_ec_ul,
                            query_overlap,
                            Some(&mut contained),
                        )
                    {
                        update = false;
                        break;
                    }
                    j -= 1;
                }
                if update {
                    let mut candidate = entries[z as usize];
                    candidate.sec = ((z as u32) << 1) + u32::from((li.qs as u64) < chain_qs);
                    result.push(candidate);
                    nv += 1;
                }
            }
        }
        previous = k;

        let li = &entries[k as usize];
        let li_v = (li.tn << 1) | u32::from(li.rev);
        let mut j = nv - 1;
        while j >= 0 {
            let lj = &result[j as usize];
            let lj_v = (lj.tn << 1) | u32::from(lj.rev);
            if lj.qe <= li.qs {
                break;
            }
            if lj.qn == u32::MAX {
                // deleted
                j -= 1;
                continue;
            }
            if lj.sec & HIGH_MARK != 0 {
                // could not be deleted
                j -= 1;
                continue;
            }
            let query_overlap = infer_rovlp_reads(&li, &lj, &reads.read_lengths);
            let mut contained = false;
            if li_v == lj_v
                || !get_ecov_contain_adv(
                    options,
                    reads,
                    li_v ^ 1,
                    lj_v ^ 1,
                    band_width,
                    diff_ec_ul,
                    query_overlap,
                    Some(&mut contained),
                )
            {
                result[j as usize].qn = u32::MAX;
                result[j as usize].tn = u32::MAX;
            }
            j -= 1;
        }
        let mut candidate = entries[k as usize];
        candidate.sec = (k as u32) << 1;
        candidate.sec |= HIGH_MARK;
        result.push(candidate);
        nv += 1;
    }

    // the tail past the last chain element
    let k = entry_count;
    let i = idn;
    if k - previous > 1 {
        let mut l = previous;
        if l < 0 {
            l = k - 1;
            while l >= 0 && entries[l as usize].qe as u64 > chain_qs {
                l -= 1;
            }
        }
        for z in l + 1..k {
            let li = &entries[z as usize];
            let li_v = (li.tn << 1) | u32::from(li.rev);
            if !li.el || vis[z as usize] & 1 != 0 {
                continue;
            }
            if !((i64::from(li.qs) + i64::from(tipl)) >= chain_qs as i64 && (li.qe as u64) < chain_qe) {
                continue;
            }
            let mut update = true;
            let mut j = i - 1;
            while j >= 0 {
                let lj = &entries[chain[j as usize]];
                let lj_v = (lj.tn << 1) | u32::from(lj.rev);
                if lj.qe <= li.qs {
                    break;
                }
                let query_overlap = infer_rovlp_reads(&li, &lj, &reads.read_lengths);
                let mut contained = false;
                if li_v == lj_v
                    || !get_ecov_contain_adv(
                        options,
                        reads,
                        li_v ^ 1,
                        lj_v ^ 1,
                        band_width,
                        diff_ec_ul,
                        query_overlap,
                        Some(&mut contained),
                    )
                {
                    update = false;
                    break;
                }
                j -= 1;
            }
            if update {
                let mut candidate = entries[z as usize];
                candidate.sec = ((z as u32) << 1) + u32::from((li.qs as u64) < chain_qs);
                result.push(candidate);
                nv += 1;
            }
        }
    }

    // drop the deleted candidates and clear the high mark
    let region_length = result.len() - region_base;
    let mut write = 0usize;
    for read in 0..region_length {
        let entry = result[region_base + read];
        if entry.qn == u32::MAX {
            continue;
        }
        let mut entry = entry;
        if entry.sec & HIGH_MARK != 0 {
            entry.sec -= HIGH_MARK;
        }
        result[region_base + write] = entry;
        write += 1;
    }
    result.truncate(region_base + write);
    nv = write as i64;

    let working = &mut result[region_base..];
    flat_contain_adv(
        reads,
        options,
        unitigs,
        band_width,
        diff_ec_ul,
        query_length,
        max_skip,
        max_iter,
        max_dis,
        extz,
        transition_score,
        targets,
        query,
        error_rate,
        trace,
        working,
        nv,
        regions,
        entries,
        &mut vec![0u64; entry_count as usize],
        rechain_window as u64,
        vis,
        bounds,
        chain_data,
    )
}

/// `gl_rchain_lin_contain_adv` (`inter.cpp` 8848): the raw-chain DP over the
/// sub-overlaps of one UL read. Returns the total chain-body length and
/// rewrites `entries` into the chain heads (`ts`/`te` indexing `chain_bodies`).
///
/// The reference's `ol` and `rid` parameters are unused on this path.
#[allow(clippy::too_many_arguments)]
pub(crate) fn gl_rchain_lin_contain_adv(
    entries: &mut Vec<UlOv>,
    chain_bodies: &mut Vec<UlOv>,
    trace: &[Rtrace],
    options: &UlGraphOptions<'_>,
    reads: &crate::reads::read_store::ReadStore,
    unitigs: &[MaUtg],
    band_width: i64,
    diff_ec_ul: f64,
    query_length: i64,
    max_skip: i64,
    max_iter: i64,
    max_dis: i64,
    chain_data: &mut ChainData,
    extz: &mut BitExtz,
    id_scratch: &mut Vec<u64>,
    transition_score: i64,
    targets: &TargetSource,
    query: &[u8],
    error_rate: f64,
    need_srt: bool,
    regions: &[crate::overlap::region::OverlapRegion],
) -> i64 {
    if entries.is_empty() {
        return 0;
    }
    let res_n = entries.len() as i64;
    chain_data.resize((res_n << 1) as usize);
    if id_scratch.len() < res_n as usize {
        id_scratch.resize(res_n as usize, 0);
    }

    if need_srt {
        // `radix_sort_ul_ov_srt_qe` is a stable 32-bit radix sort on `qe`.
        crate::container::ksort::radix_sort_by_key_32(entries, &|entry: &UlOv| entry.qe);
        {
            let mut j = 0usize;
            for i in 1..=res_n as usize {
                entries[i - 1].qs = u32::MAX - entries[i - 1].qs;
                if i == res_n as usize || entries[i].qe != entries[j].qe {
                    if i - j > 1 {
                        crate::container::ksort::radix_sort_by_key_32(
                            &mut entries[j..i],
                            &|entry: &UlOv| entry.qs,
                        );
                    }
                    j = i;
                }
            }
        }
    }

    let mut sort_keys_backup: Vec<u64> = Vec::new();
    let mut vis: Vec<u64> = vec![0; res_n as usize];
    let mut min_qs = u32::MAX;
    let mut max_qe = 0u32;
    let mut plus = 0i64;
    let mut best: UlOv = UlOv {
        qn: u32::MAX,
        ..UlOv::default()
    };
    {
        let mut scratch: Vec<i64> = vec![0; res_n as usize];
        let mut score: Vec<i32> = vec![0; res_n as usize];
        let mut predecessor: Vec<i64> = vec![-1; res_n as usize];
        let mut occurrence: Vec<i32> = vec![0; res_n as usize];
        let mut self_length: Vec<i32> = vec![0; res_n as usize]; // c_sc == rch
        let mut ssc: Vec<i32> = vec![0; res_n as usize];
        let mut st = 0i64;
        let mut max_ii = -1i64;
        for i in 0..res_n {
            entries[i as usize].qs = u32::MAX - entries[i as usize].qs;
            let li = entries[i as usize].clone();
            let li_v = (li.tn << 1) | u32::from(li.rev);
            self_length[i as usize] = i32::MAX;
            ssc[i as usize] = i32::MIN;
            let maximum_overlap = max_ovlp_src_contain(options, reads, li_v ^ 1);
            let mut x = ((li.qs as i64 + maximum_overlap) as f64 * diff_ec_ul) as i64;
            if x < band_width {
                x = band_width;
            }
            x += li.qs as i64 + maximum_overlap;
            if x > query_length + 1 {
                x = query_length + 1;
            }
            let search_limit =
                find_ul_ov_max(i as i32, entries, (x as u32).wrapping_add(G_CHAIN_INDEL)) as i64;
            let mut current_score = li.qe as i64 - li.qs as i64;
            current_score -= (li.sec as i64) * transition_score;
            let mut best_score = i32::MIN as i64;
            let mut best_index = -1i64;
            let mut best_count = 0i64;
            let mut full_contained = 0u32;
            let mut n_skip = 0i64;
            let mut trace_iterator = RtraceIter::default();
            trace_iterator.k = i32::MAX as i64;
            if (search_limit - st) > max_iter {
                st = search_limit - max_iter;
            }
            let mut j = search_limit;
            while j >= st {
                // collect potential destination vertices
                let lj = &entries[j as usize];
                if lj.qe <= li.qs {
                    break;
                }
                let mut cc = 0u32;
                let sc = connect_detect(
                    entries,
                    res_n,
                    i,
                    j,
                    self_length.as_mut_slice(),
                    options,
                    reads,
                    unitigs,
                    band_width,
                    diff_ec_ul,
                    &predecessor,
                    &score,
                    &mut trace_iterator,
                    trace,
                    targets,
                    query,
                    regions,
                    extz,
                    error_rate,
                    transition_score,
                    &mut cc,
                );
                if sc == i32::MIN as i64 {
                    j -= 1;
                    continue;
                }
                let sn = occurrence[j as usize] as i64 + 1;
                if sc > best_score || (sc == best_score && sn > best_count) {
                    best_score = sc;
                    best_index = j;
                    best_count = sn;
                    full_contained = cc;
                    if n_skip > 0 {
                        n_skip -= 1;
                    }
                } else if scratch[j as usize] == i {
                    n_skip += 1;
                    if n_skip > max_skip {
                        break;
                    }
                }
                if predecessor[j as usize] >= 0 {
                    scratch[predecessor[j as usize] as usize] = i;
                }
                j -= 1;
            }

            let end_j = j;
            if max_ii < 0 || entries[i as usize].qe as i64 > entries[max_ii as usize].qe as i64 + max_dis
            {
                // too long
                let mut max = i32::MIN as i64;
                let mut max_n = i32::MIN as i64;
                max_ii = -1;
                let mut j = i - 1;
                while j >= st
                    && (entries[i as usize].qe as i64) <= max_dis + entries[j as usize].qe as i64
                {
                    if max < score[j as usize] as i64
                        || (max == score[j as usize] as i64 && max_n < occurrence[j as usize] as i64)
                    {
                        max = score[j as usize] as i64;
                        max_n = occurrence[j as usize] as i64;
                        max_ii = j;
                    }
                    j -= 1;
                }
            }

            if max_ii >= 0 && max_ii < end_j {
                // just have a try with a[i] <-> a[max_ii]
                let lj = entries[max_ii as usize].clone();
                if lj.qe > li.qs {
                    // as max_ii < end_j, get_rid_backward_cigar_err still works
                    let mut cc = 0u32;
                    let sc = connect_detect(
                        entries,
                        res_n,
                        i,
                        max_ii,
                        self_length.as_mut_slice(),
                        options,
                        reads,
                        unitigs,
                        band_width,
                        diff_ec_ul,
                        &predecessor,
                        &score,
                        &mut trace_iterator,
                        trace,
                        targets,
                        query,
                        regions,
                        extz,
                        error_rate,
                        transition_score,
                        &mut cc,
                    );
                    if sc != i32::MIN as i64 {
                        let sn = occurrence[max_ii as usize] as i64 + 1;
                        if sc > best_score || (sc == best_score && sn > best_count) {
                            best_score = sc;
                            best_index = max_ii;
                            best_count = sn;
                            full_contained = cc;
                        }
                    }
                }
            }

            let sc = current_score;
            let sn = 1i64; // for i itself
            if sc > best_score || (sc == best_score && sn > best_count) {
                best_score = sc;
                best_index = -1;
                best_count = sn;
                full_contained = 0;
            }

            if best_score < 0 {
                best_score = current_score;
                best_index = -1;
                best_count = 1;
            }

            if best_index >= 0 {
                // no need to update mm_sc and is_fc
                propagate_transitive_reduction(
                    options,
                    reads,
                    unitigs,
                    band_width,
                    diff_ec_ul,
                    entries,
                    (search_limit + 1) as i32,
                    i,
                    self_length.as_mut_slice(),
                    &predecessor,
                    occurrence.as_mut_slice(),
                    scratch.as_mut_slice(),
                    &mut best_score,
                    &mut best_index,
                    &mut best_count,
                );
            }

            score[i as usize] = best_score as i32;
            predecessor[i as usize] = best_index;
            occurrence[i as usize] = best_count as i32;
            if best_index < 0 || (ssc[best_index as usize] < best_score as i32 && full_contained == 0)
            {
                ssc[i as usize] = best_score as i32;
            } else {
                ssc[i as usize] = ssc[best_index as usize];
            }

            if max_ii < 0
                || ((entries[i as usize].qe as i64) <= max_dis + entries[max_ii as usize].qe as i64
                    && (score[max_ii as usize] as i64) < best_score)
            {
                max_ii = i;
            }
            if (ssc[i as usize] as i64) < plus {
                plus = ssc[i as usize] as i64; // minimum negative
            }
            if std::env::var("RSHIFIASM_DP_DEBUG").is_ok() {
                eprintln!(
                    "[dp-debug] i={} f={} p={} c_n={} ssc={} tn={} qs={} qe={} ts={} te={} sec={}",
                    i, score[i as usize], predecessor[i as usize], occurrence[i as usize],
                    ssc[i as usize], entries[i as usize].tn, entries[i as usize].qs,
                    entries[i as usize].qe, entries[i as usize].ts, entries[i as usize].te,
                    entries[i as usize].sec
                );
            }
        }

        // make all f[] positive and sort by the packed key
        for i in 0..res_n as usize {
            ssc[i] -= plus as i32;
            scratch[i] = (((ssc[i] as u32 as u64) << 32) + ((i as u64) << 1)) as i64;
        }
        let mut sort_keys: Vec<u64> = scratch.iter().map(|&value| value as u64).collect();
        crate::container::ksort::radix_sort_u64(&mut sort_keys);
        sort_keys_backup = sort_keys.clone();
        plus = 0;

        let mut n_v = 0i64;
        let mut n_u = 0i64;
        min_qs = u32::MAX;
        max_qe = 0;
        for k in (0..res_n as usize).rev() {
            let n_v0 = n_v;
            let mut idn = 0usize;
            let mut i = (((sort_keys_backup[k] as u32) >> 1) as i64);
            while i >= 0 && (sort_keys_backup[i as usize] & 1) == 0 {
                if entries[i as usize].qs < min_qs {
                    min_qs = entries[i as usize].qs;
                }
                if entries[i as usize].qe > max_qe {
                    max_qe = entries[i as usize].qe;
                }
                id_scratch[idn] = i as u64;
                idn += 1;
                sort_keys_backup[i as usize] |= 1;
                i = predecessor[i as usize];
                n_v += 1;
            }
            if n_v0 == n_v {
                continue;
            }
            let sc = if i < 0 {
                (sort_keys_backup[k] >> 32) as i64
            } else {
                (sort_keys_backup[k] >> 32) as i64 - score[i as usize] as i64
            };

            if best.qn != u32::MAX && (sc as f64) < (chain_data.self_length[best.qn as usize] as f64) * 0.98
            {
                let ovlp = (min_qe_cross(max_qe, best.qe) as i64)
                    .saturating_sub(i64::from(max_qs_cross(min_qs, best.qs)))
                    .max(0);
                if ovlp >= (((max_qe - min_qs) as f64) * 0.95) as i64 {
                    n_v = n_v0;
                    continue;
                }
            }

            // reverse the chain order
            let rev_n = idn >> 1;
            for i in 0..rev_n {
                let st = id_scratch[i];
                id_scratch[i] = id_scratch[idn - i - 1];
                id_scratch[idn - i - 1] = st;
            }
            let chain: Vec<usize> = id_scratch[..idn].iter().map(|&value| value as usize).collect();
            let mut bounds = UlOv {
                qn: u32::MAX,
                ..UlOv::default()
            };
            let produced = gen_trans_aln(
                reads,
                options,
                unitigs,
                band_width,
                diff_ec_ul,
                query_length,
                &chain,
                entries,
                res_n,
                u64::from(min_qs),
                u64::from(max_qe),
                vis.as_mut_slice(),
                max_skip,
                max_iter,
                max_dis,
                extz,
                transition_score,
                targets,
                query,
                error_rate,
                trace,
                1000,
                chain_bodies,
                regions,
                chain_data,
                &mut bounds,
            );
            let produced = produced as i64;
            debug_assert!(produced > 0);
            n_v = n_v0 + produced;
            chain_data.occurrence[n_u as usize] = (n_v - n_v0) as i32;
            chain_data.self_length[n_u as usize] = sc as i32;
            n_u += 1;
            if sc < plus {
                plus = sc;
            }
            if best.qn == u32::MAX
                || sc > chain_data.self_length[best.qn as usize] as i64
                || (sc == chain_data.self_length[best.qn as usize] as i64
                    && (bounds.qe - bounds.qs) > (best.qe - best.qs))
            {
                best = bounds;
                best.qn = (n_u - 1) as u32;
            }
        }
    }

    if std::env::var("RSHIFIASM_DP_DEBUG").is_ok() {
        for index in 0.. {
            if index >= res_n as usize {
                break;
            }
            if chain_data.occurrence[index] == 0 {
                break;
            }
            eprintln!(
                "[chain-record] {} c_n={} c_sc={}",
                index, chain_data.occurrence[index], chain_data.self_length[index]
            );
        }
    }
    // rebuild `entries` as the chain heads; `c_n`/`c_sc` were written for
    // `n_u` chains, so stop at the first entry with a zero body count past
    // the first (the reference tracks `n_u` explicitly)
    let mut n_v = 0i64;
    let mut heads: Vec<UlOv> = Vec::new();
    let mut n_u = 0i64;
    while n_u < res_n {
        let body_count = chain_data.occurrence[n_u as usize];
        let body_score = chain_data.self_length[n_u as usize];
        if body_count == 0 {
            break;
        }
        let n_v0 = n_v;
        n_v += body_count as i64;
        let mut head = UlOv::default();
        head.qn = (body_score - plus as i32) as u32; // score
        head.ts = n_v0 as u32;
        head.te = n_v as u32;
        head.qs = chain_bodies[n_v0 as usize].qs;
        head.qe = chain_bodies[(n_v - 1) as usize].qe;
        heads.push(head);
        n_u += 1;
    }
    *entries = heads;
    // `radix_sort_ul_ov_srt_qn` is a stable 32-bit radix sort on `qn`.
    crate::container::ksort::radix_sort_by_key_32(entries, &|entry: &UlOv| entry.qn);
    n_v
}

/// `MIN(mrr.qe, mqe)` of the reference's chain-overlap check.
fn min_qe_cross(a: u32, b: u32) -> u32 {
    a.min(b)
}

/// `MAX(mrr.qs, mqs)` of the reference's chain-overlap check.
fn max_qs_cross(a: u32, b: u32) -> u32 {
    a.max(b)
}

/// `select_clean_chain` (`inter.cpp` 9055): deduplicates the chain heads by
/// score/overlap, renumbers the chain bodies, and rewrites their `sec`/`qn`.
pub(crate) fn select_clean_chain(
    heads: &mut Vec<UlOv>,
    bodies: &mut [UlOv],
    body_count: i64,
    ulid_local: i64,
    b64: &mut Vec<u64>,
) -> i64 {
    let idx_n0 = heads.len();
    // reverse the heads and mark them all
    for k in 0..idx_n0 / 2 {
        let kp = heads[k];
        heads[k] = heads[idx_n0 - k - 1];
        heads[idx_n0 - k - 1] = kp;
        heads[k].tn = 1;
        heads[idx_n0 - k - 1].tn = 1;
    }
    if idx_n0 & 1 == 1 {
        heads[idx_n0 / 2].tn = 1;
    }

    // filter too-close chains
    for k in 0..heads.len() {
        let m = heads[k].clone();
        let om = m.qe - m.qs;
        if m.tn == 0 {
            continue;
        }
        for i in (0..k).rev() {
            let p = &heads[i];
            let ovlp = (i64::from(m.qe.min(p.qe)) - i64::from(m.qs.max(p.qs))).max(0);
            if ovlp == 0 {
                continue;
            }
            let min_sc = p.qn.min(m.qn);
            let max_sc = p.qn.max(m.qn);
            let ok = (p.qe - p.qs).max(om);
            if (min_sc as f64) < (max_sc as f64) * 0.98 {
                break;
            }
            if ovlp > GC_OFFSET_POS as i64
                && (min_sc as f64) > (max_sc as f64) * 0.98
                && ovlp > ((ok as f64) * 0.8) as i64
            {
                heads[k].tn = 0;
                heads[i].tn = 0;
            }
        }

        for i in k + 1..heads.len() {
            let p = &heads[i];
            let ovlp = (i64::from(m.qe.min(p.qe)) - i64::from(m.qs.max(p.qs))).max(0);
            if ovlp == 0 {
                continue;
            }
            let min_sc = p.qn.min(m.qn);
            let max_sc = p.qn.max(m.qn);
            let ok = (p.qe - p.qs).max(om);
            if (min_sc as f64) < (max_sc as f64) * 0.98 {
                break;
            }
            if ovlp > GC_OFFSET_POS as i64
                && (min_sc as f64) > (max_sc as f64) * 0.98
                && ovlp > ((ok as f64) * 0.8) as i64
            {
                heads[k].tn = 0;
                heads[i].tn = 0;
            }
        }
    }

    heads.retain(|entry| entry.tn != 0);
    let mut idx_n = heads.len();

    // drop chains that overlap the kept ones by more than 10%
    let mut mchain_count = 0usize;
    for k in 0..idx_n {
        let m = heads[k].clone();
        let om = m.qe - m.qs;
        let mut i = 0usize;
        while i < mchain_count {
            let p = heads[i];
            let ovlp = (i64::from(m.qe.min(p.qe)) - i64::from(m.qs.max(p.qs))).max(0);
            if ovlp == 0 {
                i += 1;
                continue;
            }
            let ok = p.qe - p.qs;
            if ovlp > ((ok as f64) * 0.1) as i64 || ovlp > ((om as f64) * 0.1) as i64 {
                break;
            }
            i += 1;
        }
        if i < mchain_count {
            continue;
        }
        heads[mchain_count] = heads[k];
        mchain_count += 1;
    }
    heads.truncate(mchain_count);
    idx_n = mchain_count;

    // order the chains by their body start
    b64.clear();
    b64.reserve(idx_n);
    for k in 0..idx_n {
        let mut om = heads[k].ts as u64;
        om <<= 32;
        om |= k as u64;
        b64.push(om);
    }
    crate::container::ksort::radix_sort_u64(b64);

    let mut res_n = 0i64;
    for k in 0..idx_n {
        let m = &heads[(b64[k] & 0xffff_ffff) as usize];
        let mut ok = SEC_MODE;
        for z in m.ts..m.te {
            let mut entry = bodies[z as usize];
            entry.el = true;
            entry.tn |= 0x8000_0000;
            entry.sec = ok;
            entry.qn = if idx_n <= 1 {
                ulid_local as u32
            } else {
                res_n as u32
            };
            ok = res_n as u32;
            bodies[res_n as usize] = entry;
            res_n += 1;
        }
    }

    if idx_n > 1 {
        // `radix_sort_ul_ov_srt_qe` then per-`qe`-group `radix_sort_ul_ov_srt_qs`
        crate::container::ksort::radix_sort_by_key_32(
            &mut bodies[..res_n as usize],
            &|entry: &UlOv| entry.qe,
        );
        {
            let mut k = 0usize;
            for i in 1..=res_n as usize {
                if i == res_n as usize || bodies[i].qe != bodies[k].qe {
                    if i - k > 1 {
                        crate::container::ksort::radix_sort_by_key_32(
                            &mut bodies[k..i],
                            &|entry: &UlOv| entry.qs,
                        );
                    }
                    k = i;
                }
            }
        }
        b64.clear();
        b64.reserve(res_n as usize);
        for _ in 0..res_n {
            b64.push(0);
        }
        for i in 0..res_n as usize {
            b64[bodies[i].qn as usize] = i as u64;
        }
        for i in 0..res_n as usize {
            let target = b64[i] as usize;
            if bodies[target].sec != SEC_MODE {
                bodies[target].sec = b64[bodies[target].sec as usize] as u32;
            }
            bodies[target].qn = ulid_local as u32;
        }
    }
    res_n
}

/// `gen_rid_raw_chain` (`inter.cpp` 9179): the raw chains of one UL read.
///
/// `chain_index` (`cha_idx`) is the offset of the current read's overlaps in
/// `chain.tk`; the chains replace `chain.tk[chain_index..]` and the heads go
/// to `chain.lo`. Returns the number of chain bodies whose query span covers
/// at least 95% of the read (accumulated into `bb->num_read_base`).
#[allow(clippy::too_many_arguments)]
pub fn gen_rid_raw_chain(
    chain: &mut GlChain,
    chain_index: usize,
    chain_data: &mut ChainData,
    options: &UlGraphOptions<'_>,
    reads: &crate::reads::read_store::ReadStore,
    unitigs: &[MaUtg],
    diff_ec_ul: f64,
    query_length: i64,
    targets: &TargetSource,
    query: &[u8],
    extz: &mut BitExtz,
    ulid_local: i64,
    max_chain: i64,
    regions: &[crate::overlap::region::OverlapRegion],
) -> i64 {
    let transition_score = if diff_ec_ul > 0.0 {
        (1.0 / diff_ec_ul) as i64
    } else {
        0
    };
    chain.lo.clear();
    if chain.tk.len() <= chain_index {
        return 0;
    }

    let working: Vec<UlOv> = chain.tk[chain_index..].to_vec();
    let mut entries = working;
    let mut bodies: Vec<UlOv> = Vec::new();
    let mut b64: Vec<u64> = Vec::new();
    let mut body_count = gl_rchain_lin_contain_adv(
        &mut entries,
        &mut bodies,
        &chain.tc,
        options,
        reads,
        unitigs,
        G_CHAIN_BW,
        N_GCHAIN_RATE,
        query_length,
        max_chain.max(UG_SKIP_N),
        UG_ITER_N,
        UG_DIS_N,
        chain_data,
        extz,
        &mut chain.srt,
        transition_score,
        targets,
        query,
        diff_ec_ul,
        true,
        regions,
    );

    // the reference's `res` (bl->lo) IS the heads buffer; move the heads in
    std::mem::swap(&mut chain.lo, &mut entries);
    body_count = select_clean_chain(&mut chain.lo, &mut bodies, body_count, ulid_local, &mut b64);
    chain.tk.truncate(chain_index);
    chain.tk.extend_from_slice(&bodies[..body_count as usize]);

    let mut full = 0i64;
    if !chain.lo.is_empty()
        && (chain.lo[0].qe - chain.lo[0].qs) as i64 >= ((query_length as f64) * 0.95) as i64
    {
        full += 1;
    }
    full
}

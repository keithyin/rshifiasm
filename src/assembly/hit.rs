//! The overlap-list operations of the assembly stage, the counterpart of the
//! `ma_hit_*` family in `Overlaps.cpp`.
//!
//! Every overlap list is a per-read `Vec<PafHit>`; the reference's
//! `ma_hit_t_alloc` arrays are the outer slices. Hit coordinates are the
//! reference's `qns`/`qe`/`tn`/`ts`/`te` fields, and the deletion flag is
//! [`PafHit::deleted`].
//!
//! The `ma_hit_t` bit fields (`ml`, `bl`, `del`) are plain fields here. The
//! weak-overlap marker still uses bit 30 of `bl`, exactly like the reference.

use crate::assembly::asg::AsgArc;
use crate::container::ksort::introsort_by;
use crate::logging;
use crate::overlap::paf::PafHit;

/// `MA_HT_DOUBLE`: an overlap that is nearly twice as long as expected.
pub const MA_HT_DOUBLE: i32 = -1024;
/// `MA_HT_INT`: the overhang or internal fraction rejects the overlap.
pub const MA_HT_INT: i32 = -1;
/// `MA_HT_QCONT`: the query is contained in the target.
pub const MA_HT_QCONT: i32 = -2;
/// `MA_HT_TCONT`: the target is contained in the query.
pub const MA_HT_TCONT: i32 = -3;
/// `MA_HT_SHORT_OVLP`: the overlap is shorter than `min_ovlp`.
pub const MA_HT_SHORT_OVLP: i32 = -4;

/// `PRIMARY_LABLE`: the read belongs to the primary assembly.
pub const PRIMARY_LABEL: u8 = 0;
/// `ALTER_LABLE`: the read belongs to an alternate assembly.
pub const ALTER_LABEL: u8 = 1;
/// `HAP_LABLE`: the read belongs to a haplotype.
pub const HAP_LABEL: u8 = 2;
/// `FAKE_LABLE`: the read is ignored by the read-to-unitig maps.
pub const FAKE_LABEL: u8 = 4;

/// Bit 30 of `bl`, the weak-overlap marker of `clean_weak_ma_hit_t`.
pub const WEAK_OVERLAP_MARKER: u32 = 0x4000_0000;
/// `0x3fffffff`, the mask that strips the weak-overlap marker from `bl`.
pub const BLOCK_LENGTH_MASK: u32 = 0x3fff_ffff;

/// `ma_sub_t`: the high-coverage region of one read and its coverage class.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CoverageCut {
    /// `s`, the first covered base.
    pub start: u32,
    /// `e`, one past the last covered base.
    pub end: u32,
    /// `del`, whether every overlap of the read has been removed.
    pub deleted: bool,
    /// `c`, the coverage class (`PRIMARY_LABLE`, `ALTER_LABLE`, ...).
    pub class: u8,
}

/// `get_specific_overlap`: the index of the `qn -> tn` hit in `hits`.
pub fn get_specific_overlap(hits: &[PafHit], query: u32, target: u32) -> Option<usize> {
    hits.iter().position(|hit| {
        hit.target_identifier == target && hit.query_identifier == query
    })
}

/// `get_specific_overlap_with_del`: like [`get_specific_overlap`], but skips
/// deleted hits and hits that touch a read whose coverage was removed.
pub fn get_specific_overlap_with_deleted(
    sources: &[Vec<PafHit>],
    coverage_cut: &[CoverageCut],
    query: u32,
    target: u32,
) -> Option<usize> {
    if coverage_cut[query as usize].deleted || coverage_cut[target as usize].deleted {
        return None;
    }

    sources[query as usize].iter().position(|hit| {
        !hit.deleted
            && !coverage_cut[hit.query_identifier as usize].deleted
            && !coverage_cut[hit.target_identifier as usize].deleted
            && hit.target_identifier == target
            && hit.query_identifier == query
    })
}

/// `delete_single_edge`: deletes the mirror of `qn -> tn`.
pub fn delete_single_edge(
    sources: &mut [Vec<PafHit>],
    coverage_cut: &mut [CoverageCut],
    query: u32,
    target: u32,
) {
    if let Some(index) = get_specific_overlap_with_deleted(sources, coverage_cut, query, target) {
        sources[query as usize][index].deleted = true;
    }
}

/// `delete_all_edges`: deletes every overlap of read `query` and the mirrors
/// pointing at it, then marks the read itself deleted.
pub fn delete_all_edges(
    sources: &mut [Vec<PafHit>],
    coverage_cut: &mut [CoverageCut],
    query: u32,
) {
    let hits = sources[query as usize].len();
    for j in 0..hits {
        sources[query as usize][j].deleted = true;
        let target = sources[query as usize][j].target_identifier;
        let source = sources[query as usize][j].query_identifier;
        delete_single_edge(sources, coverage_cut, target, source);
    }
    coverage_cut[query as usize].deleted = true;
}

/// `set_reverse_overlap`: fills `destination` with the mirror of `source`.
///
/// The reference deliberately leaves the deletion flag of `destination`
/// untouched.
pub fn set_reverse_overlap(destination: &mut PafHit, source: &PafHit) {
    destination.query_identifier = source.target_identifier;
    destination.query_start = source.target_start;
    destination.query_end = source.target_end;

    destination.target_identifier = source.query_identifier;
    destination.target_start = source.query_start;
    destination.target_end = source.query_end;

    destination.reverse = source.reverse;
    destination.exact_region = source.exact_region;
    destination.matching_length = source.matching_length;
    destination.without_large_indel = source.without_large_indel;
    destination.block_length = source.block_length;
}

/// `normalize_ma_hit_t_single_side_advance`: makes both halves of every
/// overlap pair agree, adding a reverse hit when a pair is one-sided.
///
/// `recuse_el` and `cmk`/`chemical_coverage` mirror the reference arguments;
/// the default run passes `is_ont == 0` and a null `cmk`.
pub fn normalize_ma_hit_t_single_side_advance(
    sources: &mut [Vec<PafHit>],
    recuse_el: bool,
    cmk: Option<&[u8]>,
    chemical_coverage: u32,
) {
    let start_time = logging::realtime_seconds();
    let number_of_sources = sources.len();

    for i in 0..number_of_sources {
        let mut j = 0;
        while j < sources[i].len() {
            let hit = sources[i][j];
            let query = hit.query_identifier;
            let target = hit.target_identifier;

            let mirror = get_specific_overlap(&sources[target as usize], target, query);
            if let Some(index) = mirror {
                let mut deleted = hit.deleted || sources[target as usize][index].deleted;

                let query_length_0 = hit.query_end.wrapping_sub(hit.query_start);
                let mirror_hit = sources[target as usize][index];
                let query_length_1 =
                    mirror_hit.query_end.wrapping_sub(mirror_hit.query_start);

                if query_length_0 == query_length_1 {
                    // `query` and `target` differ, so this orders the pair
                    // consistently in both directions.
                    if query < target {
                        set_reverse_overlap(&mut sources[target as usize][index], &hit);
                    }
                } else if query_length_0 > query_length_1 {
                    set_reverse_overlap(&mut sources[target as usize][index], &hit);
                }

                if recuse_el && hit.exact_region && mirror_hit.exact_region {
                    deleted = false;
                }
                if let Some(chemical) = cmk {
                    if chemical[query as usize] as u32 <= chemical_coverage
                        || chemical[target as usize] as u32 <= chemical_coverage
                    {
                        deleted = true;
                    }
                }

                sources[i][j].deleted = deleted;
                sources[target as usize][index].deleted = deleted;
            } else {
                // The edge only exists in one direction; synthesize its
                // mirror and push it to the other read.
                let mut reverse_hit = PafHit::default();
                set_reverse_overlap(&mut reverse_hit, &hit);
                let mut deleted = true;
                if recuse_el && hit.exact_region && reverse_hit.exact_region {
                    deleted = false;
                }
                if let Some(chemical) = cmk {
                    if chemical[query as usize] as u32 <= chemical_coverage
                        || chemical[target as usize] as u32 <= chemical_coverage
                    {
                        deleted = true;
                    }
                }
                sources[i][j].deleted = deleted;
                reverse_hit.deleted = deleted;
                sources[target as usize].push(reverse_hit);
            }

            j += 1;
        }
    }

    eprintln!(
        "[M::normalize_ma_hit_t_single_side_advance] takes {:.2}s\n",
        logging::realtime_seconds() - start_time
    );
}

/// `check_weak_ma_hit`: whether a weak overlap has no strong support.
///
/// The reference returns 0 when a strong overlap that contains the weak one
/// has a mirror in the reverse list, and 1 otherwise; the boolean here is the
/// "is weak" answer, i.e. `check_weak_ma_hit(...) != 0`.
pub fn is_weak_overlap(
    aim_hits: &[PafHit],
    reverse_sources: &[Vec<PafHit>],
    weak_identifier: u32,
    weak_query_start: u32,
    weak_query_end: u32,
) -> bool {
    for hit in aim_hits {
        if !hit.deleted
            && hit.matching_length == 1
            && hit.query_start <= weak_query_start
            && hit.query_end >= weak_query_end
        {
            let strong_identifier = hit.target_identifier;
            if get_specific_overlap(
                &reverse_sources[strong_identifier as usize],
                strong_identifier,
                weak_identifier,
            )
            .is_some()
            {
                return false;
            }
        }
    }

    true
}

/// `reset_weak_ovlp`: whether `source -> destination` is backed by a
/// full-length target overlap.
fn reset_weak_overlap(
    sources: &[Vec<PafHit>],
    read_lengths: &[u64],
    source: u32,
    destination: u32,
) -> bool {
    for hit in &sources[source as usize] {
        if hit.block_length & WEAK_OVERLAP_MARKER != 0 {
            continue;
        }
        if hit.deleted {
            continue;
        }
        let target = hit.target_identifier;
        // The reference compares the 32-bit `te` against the 64-bit read
        // length, so the comparison is widened here too.
        if hit.target_start == 0 && hit.target_end as u64 == read_lengths[target as usize] {
            if let Some(index) =
                get_specific_overlap(&sources[target as usize], target, destination)
            {
                let mirror = sources[target as usize][index];
                if !mirror.deleted && mirror.block_length & WEAK_OVERLAP_MARKER == 0 {
                    return true;
                }
            }
        }
    }

    false
}

/// `clean_weak_ma_hit_t`: marks weakly supported overlaps and removes them.
///
/// A weak overlap (`ml == 0`, and low occurrence when `ou_threshold` is not
/// `(uint32_t)-1`) is marked with bit 30 of `bl` on both sides when the
/// reference's `check_weak_ma_hit` returns 0, i.e. when a strong overlap that
/// contains it has a mirror in the reverse list. Note that this is the
/// reference's literal condition (`!check_weak_ma_hit(...)`), even though the
/// function name suggests the opposite. The first pass deletes every marked
/// hit; the third pass clears the marker again.
pub fn clean_weak_ma_hit_t(
    sources: &mut [Vec<PafHit>],
    reverse_sources: &[Vec<PafHit>],
    read_lengths: &[u64],
    ou_threshold: u32,
) {
    let start_time = logging::realtime_seconds();
    let number_of_sources = sources.len();

    for i in 0..number_of_sources {
        let mut j = 0;
        while j < sources[i].len() {
            let hit = sources[i][j];
            let query = hit.query_identifier;
            let target = hit.target_identifier;

            if !hit.deleted {
                let occurrence = hit.block_length & BLOCK_LENGTH_MASK;
                if hit.matching_length == 0
                    && (ou_threshold == u32::MAX || occurrence < ou_threshold)
                {
                    // The reference marks the overlap when
                    // `!check_weak_ma_hit(...)`, i.e. when the check found
                    // support; `is_weak_overlap` is the check's negated
                    // answer, so the extra `!` keeps that literal behaviour.
                    if !is_weak_overlap(
                        &sources[query as usize],
                        reverse_sources,
                        target,
                        hit.query_start,
                        hit.query_end,
                    ) {
                        sources[i][j].block_length |= WEAK_OVERLAP_MARKER;
                        if let Some(index) =
                            get_specific_overlap(&sources[target as usize], target, query)
                        {
                            sources[target as usize][index].block_length |= WEAK_OVERLAP_MARKER;
                        }
                    }
                }
            }

            j += 1;
        }
    }

    // `update_weak_by_contain` with `ul_occ == 0`: delete every marked hit.
    for hits in sources.iter_mut() {
        for hit in hits.iter_mut() {
            if hit.block_length & WEAK_OVERLAP_MARKER != 0 {
                hit.deleted = true;
            }
        }
    }

    // `ul_occ == 1` is skipped when `ou_threshold == -1`, but the reference
    // still performs it for the ultra-long path, so it is implemented here.
    if ou_threshold != u32::MAX {
        for i in 0..number_of_sources {
            let mut k = 0;
            while k < sources[i].len() {
                let hit = sources[i][k];
                let query = hit.query_identifier;
                let target = hit.target_identifier;
                if query > target {
                    k += 1;
                    continue;
                }
                if hit.deleted && hit.block_length & WEAK_OVERLAP_MARKER != 0 {
                    if reset_weak_overlap(sources, read_lengths, query, target) {
                        if let Some(index) =
                            get_specific_overlap(&sources[target as usize], target, query)
                        {
                            sources[target as usize][index].deleted = false;
                            sources[i][k].deleted = false;
                        }
                    }
                }
                k += 1;
            }
        }
    }

    // `ul_occ == 2`: drop the marker bit again.
    for hits in sources.iter_mut() {
        for hit in hits.iter_mut() {
            if hit.block_length & WEAK_OVERLAP_MARKER != 0 {
                hit.block_length -= WEAK_OVERLAP_MARKER;
            }
        }
    }

    let _ = start_time;
}

/// `ma_hit2arc`: converts one stored overlap into a string-graph arc.
///
/// The returned integer is the reference's return value: non-negative is the
/// node length (also stored in the low half of the arc's `ul`), negative is
/// one of the `MA_HT_*` rejections. The arc is still filled in for negative
/// non-rejection returns, exactly like the reference.
pub fn ma_hit_to_arc(
    hit: &PafHit,
    query_length: i32,
    target_length: i32,
    maximum_hang: i32,
    integration_fraction: f32,
    minimum_overlap: i32,
) -> (i32, AsgArc) {
    let query_start = hit.query_start as i32;

    let target_left;
    let target_right;
    if hit.reverse {
        target_left = (target_length as u32).wrapping_sub(hit.target_end) as i32;
        target_right = hit.target_start as i32;
    } else {
        target_left = hit.target_start as i32;
        target_right = (target_length as u32).wrapping_sub(hit.target_end) as i32;
    }

    let query_end_integer = hit.query_end as i32;
    let query_remainder = query_length.wrapping_sub(query_end_integer);

    let extension_five = query_start.min(target_left);
    let extension_three = query_remainder.min(target_right);

    // The reference computes these terms in unsigned 32-bit arithmetic and
    // only then converts them to `float`.
    let query_span = hit.query_end.wrapping_sub(query_start as u32);
    let query_total = query_span
        .wrapping_add(extension_five as u32)
        .wrapping_add(extension_three as u32);
    let target_span = hit.target_end.wrapping_sub(hit.target_start);
    let target_total = target_span
        .wrapping_add(extension_five as u32)
        .wrapping_add(extension_three as u32);

    if extension_five > maximum_hang
        || extension_three > maximum_hang
        || (query_span as f32) < (query_total as f32) * integration_fraction
        || (target_span as f32) < (target_total as f32) * integration_fraction
    {
        return (MA_HT_INT, AsgArc::default());
    }

    if query_start <= target_left && query_remainder <= target_right {
        return (MA_HT_QCONT, AsgArc::default());
    }
    if query_start >= target_left && query_remainder >= target_right {
        return (MA_HT_TCONT, AsgArc::default());
    }

    let (query_end, target_end, node_length);
    if query_start > target_left {
        query_end = 0u32;
        target_end = hit.reverse as u32;
        node_length = query_start.wrapping_sub(target_left) as u32;
    } else {
        query_end = 1u32;
        target_end = !hit.reverse as u32;
        node_length = query_remainder.wrapping_sub(target_right) as u32;
    }

    if query_end_integer.wrapping_sub(query_start).wrapping_add(extension_five).wrapping_add(extension_three)
        < minimum_overlap
        || (hit.target_end as i32)
            .wrapping_sub(hit.target_start as i32)
            .wrapping_add(extension_five)
            .wrapping_add(extension_three)
            < minimum_overlap
    {
        return (MA_HT_SHORT_OVLP, AsgArc::default());
    }

    let query_vertex = query_end | (hit.query_identifier << 1);
    let target_vertex = target_end | (hit.target_identifier << 1);
    let arc = AsgArc {
        ul: ((query_vertex as u64) << 32) | node_length as u64,
        v: target_vertex,
        // `ol` is a 31-bit bit field in the reference.
        overlap: (query_length as u32).wrapping_sub(node_length) & 0x7fff_ffff,
        deleted: false,
        occurrence: 0,
        strong: hit.matching_length != 0,
        without_large_indel: hit.without_large_indel,
        exact_region: hit.exact_region as u8,
    };

    (node_length as i32, arc)
}

/// `ma_hit_sub`: the coverage region of every read.
///
/// With `minimum_depth <= 1` (the default HiFi run) every read keeps its whole
/// length. Otherwise the longest stretch covered by at least `minimum_depth`
/// overlaps is selected.
pub fn ma_hit_sub(
    minimum_depth: i32,
    sources: &[Vec<PafHit>],
    read_lengths: &[u64],
    minimum_overlap_length: i64,
) -> Vec<CoverageCut> {
    let start_time = logging::realtime_seconds();
    let read_count = sources.len();
    let mut coverage_cut = vec![CoverageCut::default(); read_count];

    for i in 0..read_count {
        if minimum_depth <= 1 {
            coverage_cut[i].start = 0;
            coverage_cut[i].end = read_lengths[i] as u32;
            coverage_cut[i].deleted = false;
            continue;
        }

        let mut events: Vec<u32> = Vec::with_capacity(sources[i].len() * 2);
        for hit in &sources[i] {
            if hit.deleted {
                continue;
            }
            let query_start = hit.query_start;
            let query_end = hit.query_end;
            events.push(query_start.wrapping_shl(1));
            events.push(query_end.wrapping_shl(1) | 1);
        }

        introsort_by(&mut events, &|left: &u32, right: &u32| left < right);

        let mut maximum = CoverageCut::default();
        let mut second_maximum = CoverageCut::default();
        let mut depth = 0i32;
        let mut region_start = 0i32;

        for &event in &events {
            let old_depth = depth;
            if event & 1 != 0 {
                depth -= 1;
            } else {
                depth += 1;
            }

            if old_depth < minimum_depth && depth >= minimum_depth {
                region_start = (event >> 1) as i32;
            } else if old_depth >= minimum_depth && depth < minimum_depth {
                let length = (event >> 1) as i32 - region_start;
                if length > (maximum.end - maximum.start) as i32 {
                    second_maximum = maximum;
                    maximum.start = region_start as u32;
                    maximum.end = event >> 1;
                } else if length > (second_maximum.end - second_maximum.start) as i32 {
                    second_maximum.start = region_start as u32;
                    second_maximum.end = event >> 1;
                }
            }
        }

        if maximum.end > maximum.start {
            coverage_cut[i].start = maximum.start;
            coverage_cut[i].end = maximum.end;
            coverage_cut[i].deleted = false;
        } else {
            coverage_cut[i].start = 0;
            coverage_cut[i].end = 0;
            coverage_cut[i].deleted = true;
        }
    }

    let _ = (start_time, minimum_overlap_length);
    coverage_cut
}

/// `ma_hit_cut`: clips every overlap to the two reads' covered regions.
pub fn ma_hit_cut(
    sources: &mut [Vec<PafHit>],
    read_lengths: &[u64],
    minimum_overlap_length: i64,
    coverage_cut: &mut [CoverageCut],
) {
    let start_time = logging::realtime_seconds();
    let read_count = sources.len();

    for i in 0..read_count {
        let mut remaining = 0i64;
        let mut j = 0;
        while j < sources[i].len() {
            let hit = sources[i][j];
            if hit.deleted {
                j += 1;
                continue;
            }

            let query_cut = coverage_cut[hit.query_identifier as usize];
            let target_cut = coverage_cut[hit.target_identifier as usize];
            if query_cut.deleted || target_cut.deleted {
                j += 1;
                continue;
            }

            // The reference assigns the unsigned arithmetic to `int`, so the
            // intermediate values are computed with wrapping 32-bit
            // operations and only compared as signed integers afterwards.
            let query_start_pre;
            let query_end_pre;
            let target_start_pre;
            let target_end_pre;
            if hit.reverse {
                query_start_pre = if hit.target_end < target_cut.end {
                    hit.query_start
                } else {
                    hit.query_start
                        .wrapping_add(hit.target_end.wrapping_sub(target_cut.end))
                };
                query_end_pre = if hit.target_start > target_cut.start {
                    hit.query_end
                } else {
                    hit.query_end
                        .wrapping_sub(target_cut.start.wrapping_sub(hit.target_start))
                };
                target_start_pre = if hit.query_end < query_cut.end {
                    hit.target_start
                } else {
                    hit.target_start
                        .wrapping_add(hit.query_end.wrapping_sub(query_cut.end))
                };
                target_end_pre = if hit.query_start > query_cut.start {
                    hit.target_end
                } else {
                    hit.target_end
                        .wrapping_sub(query_cut.start.wrapping_sub(hit.query_start))
                };
            } else {
                query_start_pre = if hit.target_start > target_cut.start {
                    hit.query_start
                } else {
                    hit.query_start
                        .wrapping_add(target_cut.start.wrapping_sub(hit.target_start))
                };
                query_end_pre = if hit.target_end < target_cut.end {
                    hit.query_end
                } else {
                    hit.query_end
                        .wrapping_sub(hit.target_end.wrapping_sub(target_cut.end))
                };
                target_start_pre = if hit.query_start > query_cut.start {
                    hit.target_start
                } else {
                    hit.target_start
                        .wrapping_add(query_cut.start.wrapping_sub(hit.query_start))
                };
                target_end_pre = if hit.query_end < query_cut.end {
                    hit.target_end
                } else {
                    hit.target_end
                        .wrapping_sub(hit.query_end.wrapping_sub(query_cut.end))
                };
            }

            let query_start = if query_start_pre > query_cut.start {
                query_start_pre.wrapping_sub(query_cut.start)
            } else {
                0
            };
            let query_end = if query_end_pre < query_cut.end {
                query_end_pre.wrapping_sub(query_cut.start)
            } else {
                query_cut.end.wrapping_sub(query_cut.start)
            };
            let target_start = if target_start_pre > target_cut.start {
                target_start_pre.wrapping_sub(target_cut.start)
            } else {
                0
            };
            let target_end = if target_end_pre < target_cut.end {
                target_end_pre.wrapping_sub(target_cut.start)
            } else {
                target_cut.end.wrapping_sub(target_cut.start)
            };

            let query_length = (query_end as i32).wrapping_sub(query_start as i32) as i64;
            let target_length = (target_end as i32).wrapping_sub(target_start as i32) as i64;
            if query_length >= minimum_overlap_length
                && target_length >= minimum_overlap_length
            {
                sources[i][j].query_start = query_start;
                sources[i][j].query_end = query_end;
                sources[i][j].target_start = target_start;
                sources[i][j].target_end = target_end;
                sources[i][j].deleted = false;
                remaining += 1;
            } else {
                sources[i][j].deleted = true;
                delete_single_edge(
                    sources,
                    coverage_cut,
                    hit.target_identifier,
                    hit.query_identifier,
                );
            }

            j += 1;
        }

        if remaining == 0 {
            coverage_cut[i].deleted = true;
        }
    }

    let _ = (start_time, read_lengths);
}

/// `ma_hit_flt`: re-checks every remaining overlap with [`ma_hit_to_arc`].
pub fn ma_hit_flt(
    sources: &mut [Vec<PafHit>],
    coverage_cut: &mut [CoverageCut],
    maximum_hang: i32,
    minimum_overlap: i32,
    integration_fraction: f32,
) {
    let start_time = logging::realtime_seconds();
    let read_count = sources.len();

    for i in 0..read_count {
        let mut remaining = 0i64;
        let mut j = 0;
        while j < sources[i].len() {
            let hit = sources[i][j];
            if hit.deleted {
                j += 1;
                continue;
            }

            let query_cut = coverage_cut[hit.query_identifier as usize];
            let target_cut = coverage_cut[hit.target_identifier as usize];
            if query_cut.deleted || target_cut.deleted {
                j += 1;
                continue;
            }

            let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
            let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
            let result = ma_hit_to_arc(
                &hit,
                query_length,
                target_length,
                maximum_hang,
                integration_fraction,
                minimum_overlap,
            )
            .0;

            if result >= 0 || result == MA_HT_QCONT || result == MA_HT_TCONT {
                sources[i][j].deleted = false;
                remaining += 1;
            } else {
                sources[i][j].deleted = true;
                delete_single_edge(
                    sources,
                    coverage_cut,
                    hit.target_identifier,
                    hit.query_identifier,
                );
            }

            j += 1;
        }

        if remaining == 0 {
            coverage_cut[i].deleted = true;
        }
    }

    eprintln!(
        "[M::ma_hit_flt] takes {:.2}s\n",
        logging::realtime_seconds() - start_time
    );
}

/// `ma_hit_contained_advance`: folds contained reads into their containers.
pub fn ma_hit_contained_advance(
    sources: &mut [Vec<PafHit>],
    coverage_cut: &mut [CoverageCut],
    read_to_unitig: &mut ReadToUnitig,
    maximum_hang: i32,
    minimum_overlap: i32,
    integration_fraction: f32,
) {
    let start_time = logging::realtime_seconds();
    let read_count = sources.len();

    for i in 0..read_count {
        if coverage_cut[i].deleted {
            continue;
        }

        let mut j = 0;
        while j < sources[i].len() {
            let hit = sources[i][j];
            let query_cut = coverage_cut[hit.query_identifier as usize];
            let target_cut = coverage_cut[hit.target_identifier as usize];
            if query_cut.deleted || target_cut.deleted {
                j += 1;
                continue;
            }
            if hit.deleted {
                j += 1;
                continue;
            }

            let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
            let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
            let result = ma_hit_to_arc(
                &hit,
                query_length,
                target_length,
                maximum_hang,
                integration_fraction,
                minimum_overlap,
            )
            .0;

            if result == MA_HT_QCONT {
                sources[i][j].deleted = true;
                delete_single_edge(
                    sources,
                    coverage_cut,
                    hit.target_identifier,
                    hit.query_identifier,
                );
                delete_all_edges(sources, coverage_cut, hit.query_identifier);
                read_to_unitig.set(hit.query_identifier, hit.target_identifier, false, None);
            } else if result == MA_HT_TCONT {
                sources[i][j].deleted = true;
                delete_single_edge(
                    sources,
                    coverage_cut,
                    hit.target_identifier,
                    hit.query_identifier,
                );
                delete_all_edges(sources, coverage_cut, hit.target_identifier);
                read_to_unitig.set(hit.target_identifier, hit.query_identifier, false, None);
            }

            j += 1;
        }
    }

    read_to_unitig.transform();

    for i in 0..read_count {
        let mut remaining = 0i64;
        for j in 0..sources[i].len() {
            let hit = sources[i][j];
            if hit.deleted {
                continue;
            }
            if !coverage_cut[hit.query_identifier as usize].deleted
                && !coverage_cut[hit.target_identifier as usize].deleted
            {
                sources[i][j].deleted = false;
                remaining += 1;
            } else {
                sources[i][j].deleted = true;
            }
        }

        if remaining == 0 {
            coverage_cut[i].deleted = true;
        }
    }

    eprintln!(
        "[M::ma_hit_contained_advance] takes {:.2} s\n",
        logging::realtime_seconds() - start_time
    );
}

/// `get_real_sources_length`: how many live arcs of `query` leave it.
pub fn get_real_sources_length(
    sources: &[Vec<PafHit>],
    coverage_cut: &[CoverageCut],
    maximum_hang: i32,
    minimum_overlap: i32,
    integration_fraction: f32,
    query: u32,
) -> u32 {
    let read = (query >> 1) as usize;
    if coverage_cut[read].deleted {
        return 0;
    }

    let mut occurrence = 0u32;
    for hit in &sources[read] {
        if hit.deleted {
            continue;
        }
        let query_cut = coverage_cut[hit.query_identifier as usize];
        let target_cut = coverage_cut[hit.target_identifier as usize];
        if query_cut.deleted || target_cut.deleted {
            continue;
        }

        let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
        let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
        let (result, arc) = ma_hit_to_arc(
            hit,
            query_length,
            target_length,
            maximum_hang,
            integration_fraction,
            minimum_overlap,
        );
        if result < 0 {
            continue;
        }
        if arc.query_vertex() == query {
            occurrence += 1;
        }
    }

    occurrence
}

/// `delete_all_edges_carefully`: removes edges only when the read keeps a
/// unique link in one direction.
pub fn delete_all_edges_carefully(
    sources: &mut [Vec<PafHit>],
    coverage_cut: &mut [CoverageCut],
    maximum_hang: i32,
    minimum_overlap: i32,
    integration_fraction: f32,
    query: u32,
) -> u32 {
    if coverage_cut[query as usize].deleted {
        return 0;
    }

    let mut keep_edge = 0u32;
    let mut direction_flags = [0u32; 2];
    let mut j = 0;
    while j < sources[query as usize].len() {
        let hit = sources[query as usize][j];
        if hit.deleted {
            j += 1;
            continue;
        }
        let query_cut = coverage_cut[hit.query_identifier as usize];
        let target_cut = coverage_cut[hit.target_identifier as usize];
        if query_cut.deleted || target_cut.deleted {
            sources[query as usize][j].deleted = true;
            delete_single_edge(
                sources,
                coverage_cut,
                hit.target_identifier,
                hit.query_identifier,
            );
            j += 1;
            continue;
        }

        let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
        let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
        let (result, arc) = ma_hit_to_arc(
            &hit,
            query_length,
            target_length,
            maximum_hang,
            integration_fraction,
            minimum_overlap,
        );
        if result < 0 {
            sources[query as usize][j].deleted = true;
            delete_single_edge(
                sources,
                coverage_cut,
                hit.target_identifier,
                hit.query_identifier,
            );
            j += 1;
            continue;
        }

        if get_real_sources_length(
            sources,
            coverage_cut,
            maximum_hang,
            minimum_overlap,
            integration_fraction,
            arc.v ^ 1,
        ) == 1
        {
            let direction = ((arc.query_vertex() ^ 1) & 1) as usize;
            direction_flags[direction] = 1;
            keep_edge += 1;
        }

        j += 1;
    }

    if keep_edge != 0 {
        let mut j = 0;
        while j < sources[query as usize].len() {
            let hit = sources[query as usize][j];
            if hit.deleted {
                j += 1;
                continue;
            }
            let query_cut = coverage_cut[hit.query_identifier as usize];
            let target_cut = coverage_cut[hit.target_identifier as usize];
            if query_cut.deleted || target_cut.deleted {
                j += 1;
                continue;
            }

            let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
            let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
            let (result, arc) = ma_hit_to_arc(
                &hit,
                query_length,
                target_length,
                maximum_hang,
                integration_fraction,
                minimum_overlap,
            );
            if result < 0 {
                j += 1;
                continue;
            }

            if direction_flags[(arc.query_vertex() & 1) as usize] == 1 {
                j += 1;
                continue;
            }

            if get_real_sources_length(
                sources,
                coverage_cut,
                maximum_hang,
                minimum_overlap,
                integration_fraction,
                arc.v ^ 1,
            ) == 1
            {
                j += 1;
                continue;
            }

            sources[query as usize][j].deleted = true;
            delete_single_edge(
                sources,
                coverage_cut,
                hit.target_identifier,
                hit.query_identifier,
            );
            j += 1;
        }
    }

    if keep_edge == 0 {
        coverage_cut[query as usize].deleted = true;
    }

    keep_edge
}

/// `get_coverage`: mean number of overlaps per covered base.
pub fn get_coverage(
    sources: &[Vec<PafHit>],
    coverage_cut: &[CoverageCut],
    read_count: usize,
) -> i64 {
    let mut read_bases = 0i64;
    let mut covered_bases = 0i64;
    for i in 0..read_count {
        read_bases += coverage_cut[i].end.wrapping_sub(coverage_cut[i].start) as i64;
        for hit in &sources[i] {
            covered_bases += hit.query_end.wrapping_sub(hit.query_start) as i64;
        }
    }

    covered_bases / read_bases
}

/// `R_to_U`: the read-to-unitig map used while contained reads are folded.
///
/// `init_R_to_U` zero-fills the array; [`renew_graph_init`] is what marks every
/// entry as unmapped, exactly like the reference's `memset`.
#[derive(Debug, Clone, Default)]
pub struct ReadToUnitig {
    index: Vec<u32>,
}

impl ReadToUnitig {
    /// `init_R_to_U`.
    pub fn new(length: usize) -> Self {
        Self {
            index: vec![0u32; length],
        }
    }

    /// Builds the map from the raw backing array, as `load_debug_graph` seeds
    /// it from a dump.
    pub fn from_entries(entries: Vec<u32>) -> Self {
        Self { index: entries }
    }

    /// The raw backing array (`ruIndex->index`).
    pub fn entries(&self) -> &[u32] {
        &self.index
    }

    /// A direct `ruIndex->index[read] = value` store. The reference writes
    /// into the array without growing it; the map grows here instead, which
    /// only happens when a read was never stored before.
    pub fn set_raw(&mut self, read: u32, value: u32) {
        if read as usize >= self.index.len() {
            self.index.resize(read as usize + 1, u32::MAX);
        }
        self.index[read as usize] = value;
    }

    /// `renew_graph_init`'s `memset(ruIndex->index, -1, ...)`.
    pub fn reset(&mut self) {
        self.index.fill(u32::MAX);
    }

    /// `set_R_to_U`.
    pub fn set(&mut self, read: u32, unitig: u32, is_unitig: bool, flag: Option<u8>) {
        if let Some(class) = flag {
            if class == FAKE_LABEL {
                return;
            }
        }

        if read as usize >= self.index.len() {
            self.index.resize(read as usize + 1, u32::MAX);
        }

        let mut value = unitig & 0x7fff_ffff;
        if is_unitig {
            value |= 1 << 31;
        }
        self.index[read as usize] = value;
    }

    /// `get_R_to_U`; both outputs are `u32::MAX` when the read is unmapped.
    pub fn get(&self, read: u32) -> (u32, u32) {
        if read as usize >= self.index.len() || self.index[read as usize] == u32::MAX {
            return (u32::MAX, u32::MAX);
        }
        let value = self.index[read as usize];
        (value & 0x7fff_ffff, value >> 31)
    }

    /// `ruIndex->len`: the number of stored entries.
    pub fn len(&self) -> usize {
        self.index.len()
    }

    /// Whether the map holds no entries.
    pub fn is_empty(&self) -> bool {
        self.index.is_empty()
    }

    /// The final pass of the `rescue_*` functions: every entry that still
    /// points at a unitig becomes unmapped (`index[v] = -1`).
    pub fn clear_unitig_entries(&mut self) {
        for value in self.index.iter_mut() {
            // `get_R_to_U` reports `is_Unitig == -1` for an unmapped entry and
            // `1` for a live one; only the latter is cleared.
            if *value != u32::MAX && (*value >> 31) == 1 {
                *value = u32::MAX;
            }
        }
    }

    /// `transfor_R_to_U`: makes every entry point directly at its unitig.
    pub fn transform(&mut self) {
        for i in 0..self.index.len() {
            let (unitig, is_unitig) = self.get(i as u32);
            if unitig == u32::MAX {
                continue;
            }
            if is_unitig == 1 {
                continue;
            }

            let mut read = unitig;
            loop {
                let (next, next_is_unitig) = self.get(read);
                if next == u32::MAX {
                    break;
                }
                if next_is_unitig == 1 {
                    break;
                }
                read = next;
            }

            self.set(i as u32, read, false, None);
        }
    }
}

/// `renew_graph_init`: clears every deletion flag and the read-to-unitig map.
///
/// The reference also destroys the previous graph and coverage array; the
/// owning pipeline drops them here.
pub fn renew_graph_init(
    sources: &mut [Vec<PafHit>],
    reverse_sources: &mut [Vec<PafHit>],
    read_to_unitig: &mut ReadToUnitig,
) {
    read_to_unitig.reset();

    for hits in sources.iter_mut() {
        for hit in hits.iter_mut() {
            hit.deleted = false;
        }
    }
    for hits in reverse_sources.iter_mut() {
        for hit in hits.iter_mut() {
            hit.deleted = false;
        }
    }
}

/// `remove_overlaps`: drops the hits at `indices` from one list.
fn remove_overlaps(hits: &mut Vec<PafHit>, indices: &[usize]) {
    let mut removed = vec![false; hits.len()];
    for &index in indices {
        removed[index] = true;
    }
    let mut kept = Vec::with_capacity(hits.len() - indices.len());
    for (hit, was_removed) in hits.iter().zip(removed.iter()) {
        if !was_removed {
            kept.push(*hit);
        }
    }
    *hits = kept;
}

/// `add_overlaps_from_different_sources`: the reverse of each selected hit,
/// ready to be appended to the rescued read's list.
///
/// The reference appends while reading the source lists; because every
/// selected hit comes from another read's list, collecting first is
/// equivalent and keeps the borrow checker happy.
fn collect_overlaps_from_different_sources(
    sources: &[Vec<PafHit>],
    read_lengths: &[u64],
    source_indices: &[u64],
) -> Vec<PafHit> {
    let mut additions = Vec::with_capacity(source_indices.len());
    for &packed in source_indices {
        let source_read = (packed >> 32) as usize;
        let source_index = packed as u32 as usize;
        let source = sources[source_read][source_index];

        additions.push(PafHit {
            query_identifier: source.target_identifier,
            query_start: source.target_start,
            query_end: source.target_end,
            target_identifier: source.query_identifier,
            target_start: source.query_start,
            target_end: source.query_end,
            reverse: source.reverse,
            matching_length: source.matching_length,
            block_length: read_lengths[source.query_identifier as usize] as u32,
            without_large_indel: source.without_large_indel,
            exact_region: source.exact_region,
            deleted: false,
        });
    }

    additions
}

/// `try_rescue_overlaps`: moves chains of reverse overlaps whose mirrors exist
/// on the forward side into the forward lists.
pub fn try_rescue_overlaps(
    sources: &mut [Vec<PafHit>],
    reverse_sources: &mut [Vec<PafHit>],
    read_lengths: &[u64],
    rescue_threshold: u64,
    is_deleted: bool,
) {
    let start_time = logging::realtime_seconds();
    let read_count = reverse_sources.len();
    let mut revised = 0i64;

    for i in 0..read_count {
        let mut event_indices: Vec<u64> = Vec::new();
        let mut hit_indices: Vec<usize> = Vec::new();

        for (j, hit) in reverse_sources[i].iter().enumerate() {
            if is_deleted && hit.deleted {
                continue;
            }
            let query = hit.query_identifier;
            let target = hit.target_identifier;
            if let Some(mirror_index) =
                get_specific_overlap(&sources[target as usize], target, query)
            {
                if !is_deleted || !sources[target as usize][mirror_index].deleted {
                    event_indices.push(((target as u64) << 32) | mirror_index as u64);
                    hit_indices.push(j);
                }
            }
        }

        if (hit_indices.len() as u64) < rescue_threshold {
            continue;
        }

        let mut events: Vec<u32> = Vec::with_capacity(hit_indices.len());
        for &index in &hit_indices {
            let hit = reverse_sources[i][index];
            events.push(hit.query_start.wrapping_shl(1));
            events.push(hit.query_end.wrapping_shl(1) | 1);
        }
        introsort_by(&mut events, &|left: &u32, right: &u32| left < right);

        let mut depth = 0i64;
        let mut region_start = 0i64;
        let mut maximum_depth = 0i64;
        let mut maximum_interval = (0u32, 0u32);

        for &event in &events {
            let old_depth = depth;
            if event & 1 != 0 {
                depth -= 1;
            } else {
                depth += 1;
            }

            if old_depth < depth {
                if depth >= maximum_depth {
                    region_start = (event >> 1) as i64;
                    maximum_depth = depth;
                }
            } else if old_depth > depth {
                if old_depth == maximum_depth {
                    maximum_interval.0 = region_start as u32;
                    maximum_interval.1 = event >> 1;
                }
            }
        }

        if (maximum_depth as u64) < rescue_threshold {
            continue;
        }

        let mut kept = 0usize;
        for j in 0..hit_indices.len() {
            let hit = reverse_sources[i][hit_indices[j]];
            if hit.query_start <= maximum_interval.0 && hit.query_end >= maximum_interval.1 {
                hit_indices[kept] = hit_indices[j];
                event_indices[kept] = event_indices[j];
                kept += 1;
            }
        }
        hit_indices.truncate(kept);
        event_indices.truncate(kept);

        remove_overlaps(&mut reverse_sources[i], &hit_indices);
        let additions =
            collect_overlaps_from_different_sources(sources, read_lengths, &event_indices);
        sources[i].extend(additions);
        revised += kept as i64;
    }

    eprintln!(
        "[M::try_rescue_overlaps] took {:.2}s, rescue edges #: {}\n",
        logging::realtime_seconds() - start_time,
        revised
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(
        query: u32,
        query_start: u32,
        query_end: u32,
        target: u32,
        target_start: u32,
        target_end: u32,
        reverse: bool,
    ) -> PafHit {
        PafHit {
            query_identifier: query,
            query_start,
            query_end,
            target_identifier: target,
            target_start,
            target_end,
            matching_length: 1,
            block_length: 100,
            reverse,
            without_large_indel: true,
            exact_region: false,
            deleted: false,
        }
    }

    #[test]
    fn reverse_overlap_swaps_both_ends() {
        let source = hit(3, 10, 60, 7, 20, 70, true);
        let mut destination = PafHit::default();
        set_reverse_overlap(&mut destination, &source);
        assert_eq!(destination.query_identifier, 7);
        assert_eq!(destination.query_start, 20);
        assert_eq!(destination.query_end, 70);
        assert_eq!(destination.target_identifier, 3);
        assert_eq!(destination.target_start, 10);
        assert_eq!(destination.target_end, 60);
        assert!(destination.reverse);
        assert_eq!(destination.block_length, 100);
    }

    #[test]
    fn normalization_adds_the_missing_mirror() {
        let mut sources = vec![
            vec![hit(0, 0, 100, 1, 0, 100, false)],
            Vec::new(),
        ];
        normalize_ma_hit_t_single_side_advance(&mut sources, false, None, 0);

        assert_eq!(sources[1].len(), 1);
        let mirror = sources[1][0];
        assert_eq!(mirror.query_identifier, 1);
        assert_eq!(mirror.query_start, 0);
        assert_eq!(mirror.query_end, 100);
        assert_eq!(mirror.target_identifier, 0);
        assert!(mirror.deleted);
        assert!(sources[0][0].deleted);
    }

    #[test]
    fn normalization_keeps_pairs_alive() {
        let mut sources = vec![
            vec![hit(0, 0, 100, 1, 0, 100, false)],
            vec![hit(1, 0, 100, 0, 0, 100, false)],
        ];
        // query 0 is shorter than the mirror's 100, so the mirror is
        // rewritten in the direction of the longer overlap; both stay alive.
        normalize_ma_hit_t_single_side_advance(&mut sources, false, None, 0);
        assert_eq!(sources[1].len(), 1);
        assert!(!sources[0][0].deleted);
        assert!(!sources[1][0].deleted);
    }

    #[test]
    fn coverage_cut_keeps_the_whole_read_without_a_depth_cutoff() {
        let sources = vec![Vec::new(), Vec::new()];
        let cuts = ma_hit_sub(0, &sources, &[101, 202], 50);
        assert_eq!(cuts[0].start, 0);
        assert_eq!(cuts[0].end, 101);
        assert_eq!(cuts[1].end, 202);
        assert!(!cuts[0].deleted);
    }

    #[test]
    fn coverage_cut_picks_the_deepest_region() {
        // Three overlaps cover [0, 20), two cover [50, 100).
        let mut list = Vec::new();
        for _ in 0..3 {
            list.push(hit(0, 0, 20, 1, 0, 20, false));
        }
        for _ in 0..2 {
            list.push(hit(0, 50, 100, 2, 0, 50, false));
        }
        let cuts = ma_hit_sub(3, &[list], &[100], 50);
        assert_eq!(cuts[0].start, 0);
        assert_eq!(cuts[0].end, 20);
        assert!(!cuts[0].deleted);
    }

    #[test]
    fn contained_query_returns_qcont() {
        // The query read spans its whole length and does not stick out of the
        // target read's flanking regions, so it is contained.
        let overlapping = hit(0, 0, 100, 1, 10, 90, false);
        let (result, _) = ma_hit_to_arc(&overlapping, 100, 100, 1000, 0.8, 50);
        assert_eq!(result, MA_HT_QCONT);
    }

    #[test]
    fn contained_target_returns_tcont() {
        // Mirror image: the target's aligned region is interior to the query.
        let overlapping = hit(0, 10, 90, 1, 0, 100, false);
        let (result, _) = ma_hit_to_arc(&overlapping, 100, 100, 1000, 0.8, 50);
        assert_eq!(result, MA_HT_TCONT);
    }

    #[test]
    fn a_normal_overlap_becomes_an_arc() {
        // The query starts 5 bases before the target read, so the
        // query-to-target arc leaves query vertex 0 with a 5-base node.
        let overlapping = hit(0, 5, 195, 1, 0, 190, false);
        let (result, arc) = ma_hit_to_arc(&overlapping, 200, 200, 1000, 0.8, 50);
        assert_eq!(result, 5);
        assert_eq!(arc.query_vertex(), 0);
        assert_eq!(arc.v, 2);
        assert_eq!(arc.node_length(), 5);
        assert_eq!(arc.overlap, 195);
    }

    #[test]
    fn a_short_overlap_is_rejected() {
        let overlapping = hit(0, 10, 20, 1, 0, 10, false);
        let (result, _) = ma_hit_to_arc(&overlapping, 20, 20, 1000, 0.8, 50);
        assert_eq!(result, MA_HT_SHORT_OVLP);
    }

    #[test]
    fn deletion_cascades_through_the_mirror() {
        let mut sources = vec![
            vec![hit(0, 0, 100, 1, 0, 100, false)],
            vec![hit(1, 0, 100, 0, 0, 100, false)],
        ];
        let mut cuts = vec![
            CoverageCut {
                start: 0,
                end: 100,
                deleted: false,
                class: 0,
            },
            CoverageCut {
                start: 0,
                end: 100,
                deleted: false,
                class: 0,
            },
        ];
        delete_all_edges(&mut sources, &mut cuts, 0);
        assert!(sources[0][0].deleted);
        assert!(sources[1][0].deleted);
        assert!(cuts[0].deleted);
    }

    #[test]
    fn read_to_unitig_follows_chains() {
        let mut map = ReadToUnitig::new(4);
        map.reset();
        map.set(1, 0, false, None);
        map.set(2, 1, false, None);
        map.transform();
        assert_eq!(map.get(2), (0, 0));
        assert_eq!(map.get(1), (0, 0));
        assert_eq!(map.get(3), (u32::MAX, u32::MAX));
    }

    #[test]
    fn rescue_moves_consistent_reverse_hits() {
        let mut sources = vec![
            vec![hit(0, 0, 100, 1, 0, 100, false)],
            vec![hit(1, 0, 100, 0, 0, 100, false)],
        ];
        let mut reverse = vec![
            vec![
                hit(0, 0, 100, 1, 0, 100, false),
                hit(0, 0, 100, 1, 0, 100, false),
                hit(0, 0, 100, 1, 0, 100, false),
                hit(0, 0, 100, 1, 0, 100, false),
            ],
            Vec::new(),
        ];
        try_rescue_overlaps(&mut sources, &mut reverse, &[100, 100], 4, false);
        assert!(reverse[0].is_empty());
        assert_eq!(sources[0].len(), 5);
    }
}

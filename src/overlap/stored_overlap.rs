//! The stored-overlap half of `worker_hap_ec`: `push_nec_re`,
//! `extract_max_exact`, `push_ne_ovlp` and `check_well_cal` from
//! `ecovlp.cpp`.
//!
//! After a read has been corrected, `wcns_gen` leaves the read's own
//! consensus cigar in the auxiliary overlap's single window. `push_nec_re`
//! copies those cells out for the two stored-overlap lists (`scc` for the
//! forward list, `scb` for the reverse one). `push_ne_ovlp` then turns every
//! overlap whose `is_match` equals the requested flag into a `ma_hit_t`
//! record, and — for the forward list only — replaces its coordinates with the
//! maximal exact sub-region the consensus cigar shares with the overlap
//! (`extract_max_exact`). `check_well_cal` finally decides whether the read is
//! fully corrected and whether its stored overlap coverage is abnormal.

use crate::container::ksort::radix_sort_u64;
use crate::correction::window::WindowCigar;
use crate::overlap::paf::PafHit;
use crate::overlap::region::OverlapRegion;
use crate::reads::read_store::ReadStore;

/// `push_nec_re`: copies the read's consensus cigar, held in the auxiliary
/// overlap's only window, into the per-read cell list.
///
/// The reference asserts the list has exactly one window, and so does this
/// port.
pub fn push_nec_re(region: &OverlapRegion) -> Vec<u16> {
    push_nec_re_from_window(&region.window_cigar)
}

/// [`push_nec_re`] over the plain window cigar `wcns_gen` fills, which is all
/// the reference keeps in the auxiliary overlap.
pub fn push_nec_re_from_window(auxiliary: &WindowCigar) -> Vec<u16> {
    assert_eq!(auxiliary.windows.len(), 1);
    let window = auxiliary.windows[0];
    let start = window.cigar_start as usize;
    let length = window.cigar_length as usize;
    auxiliary.cells[start..start + length].to_vec()
}

/// The maximal exact sub-region found by [`extract_max_exact_sub`],
/// `[query_start, query_end)` and `[target_start, target_end)`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ExactSubRegion {
    pub query_start: u64,
    pub query_end: u64,
    pub target_start: u64,
    pub target_end: u64,
}

/// Position inside the read's consensus cigar, `exk`/`eyk`/`eck`.
///
/// The cursor survives between calls: `extract_max_exact` walks the overlap's
/// windows left to right and each call continues where the previous one
/// stopped.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ExactCursor {
    pub x: i64,
    pub y: i64,
    pub cell: i64,
}

/// `pop_trace_bp_f` over a plain cell slice: coalesces the cell at `index`
/// with the following cells of the same operation and the same two bases.
fn pop_trace_bp_f(cells: &[u16], mut index: usize) -> (usize, u8, u32) {
    let operation = (cells[index] >> 14) as u8;
    let (mut base_query, mut base_target, mut length) = match operation {
        2 | 3 => (-1i16, ((cells[index] >> 12) & 3) as i16, u32::from(cells[index] & 0xfff)),
        1 => (
            ((cells[index] >> 10) & 3) as i16,
            ((cells[index] >> 12) & 3) as i16,
            u32::from(cells[index] & 0x3ff),
        ),
        _ => (-1, -1, u32::from(cells[index] & 0x3fff)),
    };
    index += 1;
    while index < cells.len() && (cells[index] >> 14) as u8 == operation {
        let (next_query, next_target, next_length) = match operation {
            2 | 3 => (
                -1i16,
                ((cells[index] >> 12) & 3) as i16,
                u32::from(cells[index] & 0xfff),
            ),
            1 => (
                ((cells[index] >> 10) & 3) as i16,
                ((cells[index] >> 12) & 3) as i16,
                u32::from(cells[index] & 0x3ff),
            ),
            _ => (-1, -1, u32::from(cells[index] & 0x3fff)),
        };
        if base_query != next_query || base_target != next_target {
            break;
        }
        length += next_length;
        base_query = next_query;
        base_target = next_target;
        index += 1;
    }
    (index, operation, length)
}

/// `extract_max_exact_sub`: the longest match run inside one match cell of
/// the read's consensus cigar that also lies inside `[x_start, x_end)`.
///
/// The search resumes from `cursor`, which is advanced to the end of the
/// scanned span. The returned region, when present, is the strictly longest
/// one seen so far in this call; the reference resets its best region at every
/// entry, and so does this port.
pub fn extract_max_exact_sub(
    cells: &[u16],
    x_start: i64,
    x_end: i64,
    y_start: i64,
    y_end: i64,
    cursor: &mut ExactCursor,
) -> Option<ExactSubRegion> {
    let cell_count = cells.len() as i64;
    let mut x = cursor.x;
    let mut y = cursor.y;
    let mut cell = cursor.cell;
    // The reference takes `ye0` but never reads it; the target coordinates
    // are derived from `ys0` and the match run alone.
    let _ = y_end;

    if cell < 0 || cell > cell_count {
        cell = 0;
        x = 0;
        y = 0;
    }

    // Walk back to the start of the span; the reference re-reads the cell at
    // the decremented index and subtracts its length.
    while cell > 0 && x >= x_start {
        cell -= 1;
        let packed = cells[cell as usize];
        let operation = (packed >> 14) as u8;
        let length = match operation {
            2 | 3 => u32::from(packed & 0xfff),
            1 => u32::from(packed & 0x3ff),
            _ => u32::from(packed & 0x3fff),
        };
        if operation != 2 {
            x -= i64::from(length);
        }
        if operation != 3 {
            y -= i64::from(length);
        }
    }

    let mut best = ExactSubRegion::default();
    while cell < cell_count && x < x_end {
        let window_x_start = x;
        let window_y_start = y;
        let (next_cell, operation, length) = pop_trace_bp_f(cells, cell as usize);
        if operation != 2 {
            x += i64::from(length);
        }
        if operation != 3 {
            y += i64::from(length);
        }
        let window_x_end = x;
        if operation == 0 {
            let overlap_start = x_start.max(window_x_start);
            let overlap_end = x_end.min(window_x_end);
            let overlap_length = if overlap_end > overlap_start {
                overlap_end - overlap_start
            } else {
                0
            };
            if overlap_length > 0
                && overlap_length > (best.query_end - best.query_start) as i64
            {
                best.query_start = (window_y_start + overlap_start - window_x_start) as u64;
                best.query_end = (window_y_start + overlap_end - window_x_start) as u64;
                best.target_start = (y_start + overlap_start - x_start) as u64;
                best.target_end = (y_start + overlap_end - x_start) as u64;
            }
        }
        cell = next_cell as i64;
    }

    cursor.x = x;
    cursor.y = y;
    cursor.cell = cell;

    if best.query_end > best.query_start {
        Some(best)
    } else {
        None
    }
}

/// `pop_trace` over one window of the overlap's window cigar: coalesces
/// consecutive cells of the same operation, using the plain 14-bit length.
fn pop_trace_window(cells: &[u16], start: usize, length: usize, offset: usize) -> (usize, u8, u32) {
    let packed = cells[start + offset];
    let operation = (packed >> 14) as u8;
    let mut run_length = u32::from(packed & 0x3fff);
    let mut index = offset + 1;
    while index < length && (cells[start + index] >> 14) as u8 == operation {
        run_length += u32::from(cells[start + index] & 0x3fff);
        index += 1;
    }
    (index, operation, run_length)
}

/// `extract_max_exact`: the longest exact sub-region shared by the overlap's
/// window cigar and the read's consensus cigar (`exact_cells`).
///
/// The two cigars cover the same read, so the query coordinates of the first
/// one are matched against the target coordinates of the second. Returns
/// `None` when no match run overlaps.
pub fn extract_max_exact(region: &OverlapRegion, exact_cells: &[u16]) -> Option<ExactSubRegion> {
    let mut best_x = [0u64; 2];
    let mut best_y = [0u64; 2];
    let mut cursor = ExactCursor::default();

    for window in &region.window_cigar.windows {
        let start = window.cigar_start as usize;
        let length = window.cigar_length as usize;
        let mut x = i64::from(window.x_start);
        let mut y = i64::from(window.y_start);
        let mut offset = 0usize;
        while offset < length {
            let window_x_start = x;
            let window_y_start = y;
            let (next, operation, run_length) = pop_trace_window(&region.window_cigar.cells, start, length, offset);
            if operation != 2 {
                x += i64::from(run_length);
            }
            if operation != 3 {
                y += i64::from(run_length);
            }
            let window_x_end = x;
            offset = next;
            if window_x_end <= window_x_start {
                continue;
            }
            if (window_x_end - window_x_start) <= (best_x[1] - best_x[0]) as i64 {
                continue;
            }
            if operation != 0 {
                continue;
            }
            let found = extract_max_exact_sub(
                exact_cells,
                window_x_start,
                window_x_end,
                window_y_start,
                y,
                &mut cursor,
            );
            if let Some(found) = found {
                if (found.query_end - found.query_start) > best_x[1] - best_x[0] {
                    best_x[0] = found.query_start;
                    best_x[1] = found.query_end;
                    best_y[0] = found.target_start;
                    best_y[1] = found.target_end;
                }
            }
        }
    }

    if best_x[1] > best_x[0] {
        Some(ExactSubRegion {
            query_start: best_x[0],
            query_end: best_x[1],
            target_start: best_y[0],
            target_end: best_y[1],
        })
    } else {
        None
    }
}

/// `push_ne_ovlp`: stores every overlap whose refinement outcome equals
/// `flag`, in list order.
///
/// `exact_cells` is the read's consensus cigar ([`push_nec_re`]) for the
/// forward list; the reverse list passes `None`, and then no exact sub-region
/// is searched for. The stored coordinates are half-open and, for the reverse
/// list, stay in the reverse-complement frame the chain used.
pub fn push_ne_ovlp(
    overlaps: &[OverlapRegion],
    flag: u8,
    reads: &ReadStore,
    exact_cells: Option<&[u16]>,
) -> Vec<PafHit> {
    let selected = overlaps.iter().filter(|region| region.is_match == flag).count();
    let mut hits = Vec::with_capacity(selected);

    for region in overlaps {
        if region.is_match != flag {
            continue;
        }

        let mut hit = PafHit {
            query_identifier: region.x_identifier,
            query_start: region.x_position_start,
            query_end: region.x_position_end + 1,
            target_identifier: region.y_identifier,
            target_start: region.y_position_start,
            target_end: region.y_position_end + 1,
            // `ml` is a one-bit field in the reference.
            matching_length: (region.strong as u32) & 1,
            block_length: reads.read_length(region.y_identifier as usize) as u32,
            reverse: region.y_strand == 1,
            without_large_indel: region.without_large_indel != 0,
            exact_region: false,
            deleted: false,
        };

        if let Some(cells) = exact_cells {
            if let Some(found) = extract_max_exact(region, cells) {
                hit.query_start = found.query_start as u32;
                hit.query_end = found.query_end as u32;
                hit.target_start = found.target_start as u32;
                hit.target_end = found.target_end as u32;
                hit.exact_region = true;
            }
        }

        hits.push(hit);
    }

    hits
}

/// `check_well_cal`: decides whether the read is fully corrected and whether
/// its stored overlap coverage looks abnormal.
///
/// `cells` is the read's consensus cigar. `indices` is the reference's `idx`
/// scratch; it is cleared and refilled. `length` is the read length and
/// `minimum_depth` the coverage the read must keep over every stretch,
/// `MIN_COVERAGE_THRESHOLD * 2` at the call site.
pub fn check_well_cal(
    cells: &[u16],
    indices: &mut Vec<u64>,
    fully_corrected: &mut bool,
    abnormal: &mut u8,
    length: i64,
    minimum_depth: i64,
    hits: &[PafHit],
) {
    *fully_corrected = true;
    *abnormal = 0;
    indices.clear();

    for hit in hits {
        let start = u64::from(hit.query_start);
        let end = u64::from(hit.query_end);
        indices.push(start << 1);
        indices.push((end << 1) | 1);
    }

    radix_sort_u64(indices);

    let mut depth: i64 = 0;
    let mut start_position: i64 = 0;
    for &value in indices.iter() {
        let old_depth = depth;
        if value & 1 != 0 {
            depth -= 1;
        } else {
            depth += 1;
        }

        let end_position = (value >> 1) as i64;
        if end_position > start_position {
            if old_depth < minimum_depth {
                *fully_corrected = false;
            }
            if old_depth == 0 {
                if start_position > 0 && end_position < length {
                    *abnormal = 1;
                } else if *abnormal == 0 {
                    *abnormal = 2;
                }
            }
        }

        start_position = end_position;
    }

    let end_position = length;
    let old_depth = depth;
    if end_position > start_position {
        if old_depth < minimum_depth {
            *fully_corrected = false;
        }
        if old_depth == 0 {
            if start_position > 0 && end_position < length {
                *abnormal = 1;
            } else if *abnormal == 0 {
                *abnormal = 2;
            }
        }
    }

    if *fully_corrected {
        // Every cell must be a match once the read is considered corrected.
        if cells.iter().any(|cell| (cell >> 14) != 0) {
            *fully_corrected = false;
        }
    }
}

/// `pop_trace_bp_f` keeping the base fields: a base-aware run starting at
/// `index`, coalescing the following cells of the same operation and bases.
///
/// Returns the first unconsumed index, the operation, the query and target
/// base codes (`u16::MAX` when the operation does not carry one) and the run
/// length. The reference swaps the two bases for a deletion cell before it
/// returns; the same swap is applied here.
pub fn pop_trace_bp_f_full(cells: &[u16], mut index: usize) -> (usize, u16, u16, u16, u32) {
    let operation = cells[index] >> 14;
    let mut base_query = u16::MAX;
    let mut base_target = u16::MAX;
    let mut length;
    match operation {
        2 | 3 => {
            base_target = (cells[index] >> 12) & 3;
            length = u32::from(cells[index] & 0xfff);
        }
        1 => {
            base_target = (cells[index] >> 12) & 3;
            base_query = (cells[index] >> 10) & 3;
            length = u32::from(cells[index] & 0x3ff);
        }
        _ => length = u32::from(cells[index] & 0x3fff),
    }
    index += 1;
    while index < cells.len() && cells[index] >> 14 == operation {
        let mut next_query = u16::MAX;
        let mut next_target = u16::MAX;
        let next_length = match operation {
            2 | 3 => {
                next_target = (cells[index] >> 12) & 3;
                u32::from(cells[index] & 0xfff)
            }
            1 => {
                next_target = (cells[index] >> 12) & 3;
                next_query = (cells[index] >> 10) & 3;
                u32::from(cells[index] & 0x3ff)
            }
            _ => u32::from(cells[index] & 0x3fff),
        };
        if base_query != next_query || base_target != next_target {
            break;
        }
        length += next_length;
        index += 1;
    }
    if operation == 3 {
        base_query = base_target;
        base_target = u16::MAX;
    }
    (index, operation, base_query, base_target, length)
}

/// `pop_trace_bp_rev_f`: the backward form of [`pop_trace_bp_f_full`].
///
/// The walk starts at `index` and consumes cells towards the front, stopping
/// at the first cell that differs in operation or bases. The returned index is
/// the first cell that was **not** consumed and may be `-1`.
pub fn pop_trace_bp_rev_f(cells: &[u16], index: i64) -> (i64, u16, u16, u16, u32) {
    let mut position = index;
    let operation = cells[position as usize] >> 14;
    let mut base_query = u16::MAX;
    let mut base_target = u16::MAX;
    let mut length;
    match operation {
        2 | 3 => {
            base_target = (cells[position as usize] >> 12) & 3;
            length = u32::from(cells[position as usize] & 0xfff);
        }
        1 => {
            base_target = (cells[position as usize] >> 12) & 3;
            base_query = (cells[position as usize] >> 10) & 3;
            length = u32::from(cells[position as usize] & 0x3ff);
        }
        _ => length = u32::from(cells[position as usize] & 0x3fff),
    }
    position -= 1;
    while position >= 0 && cells[position as usize] >> 14 == operation {
        let mut next_query = u16::MAX;
        let mut next_target = u16::MAX;
        let next_length = match operation {
            2 | 3 => {
                next_target = (cells[position as usize] >> 12) & 3;
                u32::from(cells[position as usize] & 0xfff)
            }
            1 => {
                next_target = (cells[position as usize] >> 12) & 3;
                next_query = (cells[position as usize] >> 10) & 3;
                u32::from(cells[position as usize] & 0x3ff)
            }
            _ => u32::from(cells[position as usize] & 0x3fff),
        };
        if base_query != next_query || base_target != next_target {
            break;
        }
        length += next_length;
        position -= 1;
    }
    if operation == 3 {
        base_query = base_target;
        base_target = u16::MAX;
    }
    (position, operation, base_query, base_target, length)
}

/// `adjust_exact_match`: the longest match run of a consensus cigar that also
/// lies inside a stored overlap.
///
/// The cigar maps the read's pre-correction coordinates (`x`) to its
/// post-correction coordinates (`y`). `[x_start, x_end)` and `[y_start,
/// y_end)` describe the stored overlap in those two frames; the returned
/// [`ExactSubRegion`] holds the refined span in the `y` frame
/// (`query_start`/`query_end`) and the matching span of the stored overlap
/// (`target_start`/`target_end`), both half-open.
///
/// The reference scans the cigar forwards or backwards depending on `reverse`
/// and keeps the **strictly** longest match run it sees; the returned region
/// is all-zero when nothing matches.
pub fn adjust_exact_match(
    cells: &[u16],
    x_start: i64,
    x_end: i64,
    y_start: i64,
    y_end: i64,
    reverse: bool,
) -> ExactSubRegion {
    let mut result = ExactSubRegion::default();
    if x_end <= x_start || y_end <= y_start {
        return result;
    }

    let cell_count = cells.len() as i64;
    let mut x = 0i64;
    let mut y = 0i64;
    let mut best = 0i64;

    if !reverse {
        let mut cell = 0i64;
        while cell < cell_count && x < x_end {
            let window_x_start = x;
            let window_y_start = y;
            let (next, operation, length) = pop_trace_bp_f(cells, cell as usize);
            if operation != 2 {
                x += i64::from(length);
            }
            if operation != 3 {
                y += i64::from(length);
            }
            let window_x_end = x;
            if operation == 0 {
                let overlap_start = x_start.max(window_x_start);
                let overlap_end = x_end.min(window_x_end);
                let overlap = if overlap_end > overlap_start {
                    overlap_end - overlap_start
                } else {
                    0
                };
                if overlap > 0 && overlap > best {
                    result.query_start =
                        (window_y_start + overlap_start - window_x_start) as u64;
                    result.query_end = (window_y_start + overlap_end - window_x_start) as u64;
                    result.target_start = (y_start + overlap_start - x_start) as u64;
                    result.target_end = (y_start + overlap_end - x_start) as u64;
                    best = overlap;
                }
            }
            cell = next as i64;
        }
    } else {
        let mut cell = cell_count - 1;
        while cell >= 0 && x < x_end {
            let window_x_start = x;
            let window_y_start = y;
            let (next, operation, _, _, length) = pop_trace_bp_rev_f(cells, cell);
            if operation != 2 {
                x += i64::from(length);
            }
            if operation != 3 {
                y += i64::from(length);
            }
            let window_x_end = x;
            if operation == 0 {
                let overlap_start = x_start.max(window_x_start);
                let overlap_end = x_end.min(window_x_end);
                let overlap = if overlap_end > overlap_start {
                    overlap_end - overlap_start
                } else {
                    0
                };
                if overlap > 0 && overlap > best {
                    result.query_start =
                        (window_y_start + overlap_start - window_x_start) as u64;
                    result.query_end = (window_y_start + overlap_end - window_x_start) as u64;
                    result.target_start = (y_start + overlap_start - x_start) as u64;
                    result.target_end = (y_start + overlap_end - x_start) as u64;
                    best = overlap;
                }
            }
            cell = next;
        }
    }

    result
}

/// `quick_exact_match`: re-anchors a stored overlap through the target read's
/// consensus cigar and checks that the corrected reads still agree exactly.
///
/// The hit is rewritten in place in the target's post-correction frame even
/// when the check fails, exactly like the reference: the coordinates are
/// adjusted by [`adjust_exact_match`], clamped to the two read lengths, slid
/// towards the read ends so that at least one side starts at zero, and finally
/// extended to the ends. Returns whether the two corrected sequences agree over
/// the resulting span.
pub fn quick_exact_match(
    hit: &mut PafHit,
    reads: &ReadStore,
    query: &[u8],
    consensus_cigars: &[Vec<u16>],
) -> bool {
    let target_identifier = hit.target_identifier as usize;
    let mut matched = false;

    let adjusted = adjust_exact_match(
        &consensus_cigars[target_identifier],
        i64::from(hit.target_start),
        i64::from(hit.target_end),
        i64::from(hit.query_start),
        i64::from(hit.query_end),
        hit.reverse,
    );
    if adjusted.query_end > adjusted.query_start {
        hit.target_start = adjusted.query_start as u32;
        hit.target_end = adjusted.query_end as u32;
        hit.query_start = adjusted.target_start as u32;
        hit.query_end = adjusted.target_end as u32;
        matched = true;
    }

    let query_length = query.len() as i64;
    let target_length = reads.read_length(target_identifier) as i64;
    let mut query_start = i64::from(hit.query_start);
    let mut query_end = i64::from(hit.query_end);
    let mut target_start = i64::from(hit.target_start);
    let mut target_end = i64::from(hit.target_end);

    if query_start >= query_length {
        query_start = query_length;
    }
    if query_end > query_length {
        query_end = query_length;
    }
    if query_end <= query_start {
        matched = false;
    }
    if target_start >= target_length {
        target_start = target_length;
    }
    if target_end > target_length {
        target_end = target_length;
    }
    if target_end <= target_start {
        matched = false;
    }
    if query_end - query_start != target_end - target_start {
        matched = false;
    }

    if query_start <= target_start {
        target_start -= query_start;
        query_start = 0;
    } else {
        query_start -= target_start;
        target_start = 0;
    }

    let query_remaining = query_length - query_end;
    let target_remaining = target_length - target_end;
    if query_remaining <= target_remaining {
        query_end = query_length;
        target_end += query_remaining;
    } else {
        target_end = target_length;
        query_end += target_remaining;
    }

    hit.query_start = query_start as u32;
    hit.query_end = query_end as u32;
    hit.target_start = target_start as u32;
    hit.target_end = target_end as u32;

    if matched
        && target_end - target_start == query_end - query_start
        && query_end > query_start
    {
        let target = reads.read_sub_region(
            target_identifier,
            target_start,
            target_end - target_start,
            u8::from(hit.reverse),
        );
        if query[query_start as usize..query_end as usize] == target[..] {
            return true;
        }
    }

    false
}

#[cfg(test)]
mod probe_goldens {
    use super::*;
    use crate::correction::dedup_chains::dedup_chains;

    const READ_LENGTHS: [usize; 6] = [100, 200, 150, 300, 250, 180];

    fn read_store() -> ReadStore {
        let mut reads = ReadStore::default();
        for (index, &length) in READ_LENGTHS.iter().enumerate() {
            reads.push_read(&format!("read{}", index), &vec![b'A'; length]);
        }
        reads
    }

    // ------------------------------------------------------------ encoding

    /// `push_trace`: the plain 14-bit length, no base fields.
    fn push_plain(cells: &mut Vec<u16>, operation: u16, mut length: u32) {
        let encoded = operation << 14;
        while length >= 0x3fff {
            cells.push(encoded + 0x3fff);
            length -= 0x3fff;
        }
        if length > 0 {
            cells.push(encoded + length as u16);
        }
    }

    /// `push_trace_bp_f`: full base encoding, merging with the last cell when
    /// its operation and bases match (`is_append == 1`).
    fn push_encoded(
        cells: &mut Vec<u16>,
        operation: u16,
        mut base_query: u16,
        mut base_target: u16,
        mut length: u32,
    ) {
        if operation == 3 {
            base_target = base_query;
            base_query = u16::MAX;
        }
        if let Some(&last) = cells.last() {
            let last_operation = last >> 14;
            // `pop_trac_bpc_f` leaves `rbq` untouched for D/I cells.
            let last_query = if last_operation == 1 {
                (last >> 10) & 3
            } else {
                base_query
            };
            let last_target = match last_operation {
                1..=3 => (last >> 12) & 3,
                _ => base_target,
            };
            let last_length = match last_operation {
                1 => u32::from(last & 0x3ff),
                2 | 3 => u32::from(last & 0xfff),
                _ => u32::from(last & 0x3fff),
            };
            if operation == last_operation
                && base_query == last_query
                && base_target == last_target
            {
                cells.pop();
                length += last_length;
            }
        }

        let (mask, encoded) = match operation {
            2 | 3 => (
                0xfffu32,
                (u32::from(operation) << 14) + (u32::from(base_target & 3) << 12),
            ),
            1 => (
                0x3ff,
                (u32::from(operation) << 14)
                    + (u32::from(base_target & 3) << 12)
                    + (u32::from(base_query & 3) << 10),
            ),
            _ => (0x3fff, u32::from(operation) << 14),
        };
        while length >= mask {
            cells.push((encoded + mask) as u16);
            length -= mask;
        }
        if length > 0 {
            cells.push((encoded + length) as u16);
        }
    }

    /// `pop_trace_bp_f` over the whole cell list, including its coalescing and
    /// the D/I base swap on return.
    fn print_decoded_cells(lines: &mut Vec<String>, tag: &str, cells: &[u16]) {
        let mut line = format!("{} n={} :", tag, cells.len());
        let mut index = 0usize;
        while index < cells.len() {
            let operation = cells[index] >> 14;
            let (mut base_query, mut base_target, mut length) = match operation {
                2 | 3 => (
                    u16::MAX,
                    (cells[index] >> 12) & 3,
                    u32::from(cells[index] & 0xfff),
                ),
                1 => (
                    (cells[index] >> 10) & 3,
                    (cells[index] >> 12) & 3,
                    u32::from(cells[index] & 0x3ff),
                ),
                _ => (u16::MAX, u16::MAX, u32::from(cells[index] & 0x3fff)),
            };
            index += 1;
            while index < cells.len() && (cells[index] >> 14) == operation {
                let (next_query, next_target, next_length) = match operation {
                    2 | 3 => (
                        u16::MAX,
                        (cells[index] >> 12) & 3,
                        u32::from(cells[index] & 0xfff),
                    ),
                    1 => (
                        (cells[index] >> 10) & 3,
                        (cells[index] >> 12) & 3,
                        u32::from(cells[index] & 0x3ff),
                    ),
                    _ => (u16::MAX, u16::MAX, u32::from(cells[index] & 0x3fff)),
                };
                if next_query != base_query || next_target != base_target {
                    break;
                }
                length += next_length;
                base_query = next_query;
                base_target = next_target;
                index += 1;
            }
            if operation == 3 {
                base_query = base_target;
                base_target = u16::MAX;
            }
            line.push_str(&format!(
                " {}{}(q={},t={})",
                length,
                b"EMDI"[operation as usize] as char,
                base_query,
                base_target
            ));
        }
        lines.push(line);
    }

    // --------------------------------------------------------- construction

    struct DedupSpec {
        y_identifier: u32,
        x_start: u32,
        x_end: u32,
        is_match: u8,
        non_homopolymer_errors: u32,
    }

    fn build_dedup_region(spec: &DedupSpec) -> OverlapRegion {
        OverlapRegion {
            y_identifier: spec.y_identifier,
            x_position_start: spec.x_start,
            x_position_end: spec.x_end,
            is_match: spec.is_match,
            non_homopolymer_errors: spec.non_homopolymer_errors,
            ..OverlapRegion::default()
        }
    }

    struct WindowSpec {
        x_start: i32,
        x_end: i32,
        y_start: i32,
        y_end: i32,
        cells: Vec<(u16, u32)>,
    }

    fn build_window_region(
        y_identifier: u32,
        strand: u32,
        is_match: u8,
        windows: &[WindowSpec],
    ) -> OverlapRegion {
        let mut region = OverlapRegion {
            x_identifier: 0,
            y_identifier,
            y_strand: strand,
            is_match,
            ..OverlapRegion::default()
        };
        for window in windows {
            let index = region.window_cigar.append_window(
                window.x_start,
                window.x_end,
                window.y_start,
                window.y_end,
                0,
                0,
                0,
                0,
            );
            let start = region.window_cigar.cells.len();
            for &(operation, length) in &window.cells {
                push_plain(&mut region.window_cigar.cells, operation, length);
            }
            region.window_cigar.assign_cells(index, start);
        }
        region
    }

    fn build_consensus(spec: &[(u16, u16, u16, u32)]) -> Vec<u16> {
        let mut cells = Vec::new();
        for &(operation, base_query, base_target, length) in spec {
            if operation == 0 {
                push_plain(&mut cells, 0, length);
            } else {
                push_encoded(&mut cells, operation, base_query, base_target, length);
            }
        }
        cells
    }

    fn build_consensus_region(cells: &[u16]) -> OverlapRegion {
        let mut region = OverlapRegion::default();
        let index = region.window_cigar.append_window(0, 0, 0, 0, 0, 0, 0, 0);
        let start = region.window_cigar.cells.len();
        region.window_cigar.cells.extend_from_slice(cells);
        region.window_cigar.assign_cells(index, start);
        region
    }

    fn build_push_region(spec: &DedupSpec) -> OverlapRegion {
        let span = spec.x_end + 1 - spec.x_start;
        let mut region = OverlapRegion {
            x_identifier: 0,
            y_identifier: spec.y_identifier,
            x_position_start: spec.x_start,
            x_position_end: spec.x_end,
            y_position_start: 40,
            y_position_end: 40 + (spec.x_end - spec.x_start),
            y_strand: spec.y_identifier % 2,
            is_match: spec.is_match,
            strong: 1,
            without_large_indel: 1,
            non_homopolymer_errors: spec.non_homopolymer_errors,
            ..OverlapRegion::default()
        };
        let index = region.window_cigar.append_window(
            spec.x_start as i32,
            spec.x_end as i32,
            40,
            region.y_position_end as i32,
            0,
            0,
            0,
            0,
        );
        let start = region.window_cigar.cells.len();
        let mut offset = 0;
        while offset < span {
            let length = (span - offset).min(10);
            push_plain(&mut region.window_cigar.cells, 0, length);
            offset += 10;
        }
        region.window_cigar.assign_cells(index, start);
        region
    }

    // ------------------------------------------------------------ printing

    fn push_dedup_line(lines: &mut Vec<String>, index: usize, state: &str, region: &OverlapRegion) {
        lines.push(format!(
            " DDP i={} {} y={} x=[{},{}] aim={} nhpe={}",
            index,
            state,
            region.y_identifier,
            region.x_position_start,
            region.x_position_end,
            region.is_match,
            region.non_homopolymer_errors
        ));
    }

    fn run_dedup_batch(lines: &mut Vec<String>, _name: &str, specs: &[DedupSpec]) {
        let mut overlaps: Vec<OverlapRegion> = specs.iter().map(build_dedup_region).collect();
        let original = overlaps.len();
        lines.push(format!("DDB len={} on={}", overlaps.len(), original));
        for (index, region) in overlaps.iter().enumerate() {
            push_dedup_line(lines, index, "KEEP", region);
        }

        dedup_chains(&mut overlaps);

        lines.push(format!("DDA len={} on={}", overlaps.len(), original));
        for (index, region) in overlaps.iter().enumerate() {
            push_dedup_line(lines, index, "KEEP", region);
        }
    }

    fn run_extract_batch(
        lines: &mut Vec<String>,
        name: &str,
        consensus: &[u16],
        windows: &[WindowSpec],
    ) {
        let region = build_window_region(1, 0, 1, windows);
        match extract_max_exact(&region, consensus) {
            Some(found) => lines.push(format!(
                "EXB name={} ret=1 eq=[{},{}) et=[{},{})",
                name, found.query_start, found.query_end, found.target_start, found.target_end
            )),
            None => lines.push(format!(
                "EXB name={} ret=0 eq=[4294967295,4294967295) et=[4294967295,4294967295)",
                name
            )),
        }
    }

    fn run_push_ne_batch(
        lines: &mut Vec<String>,
        name: &str,
        consensus: &[u16],
        specs: &[DedupSpec],
        flag: u8,
        with_exact: bool,
    ) {
        let regions: Vec<OverlapRegion> = specs.iter().map(build_push_region).collect();
        let reads = read_store();
        let hits = push_ne_ovlp(
            &regions,
            flag,
            &reads,
            if with_exact { Some(consensus) } else { None },
        );

        lines.push(format!(
            "PNB name={} flag={} ec={}",
            name,
            flag,
            u32::from(with_exact)
        ));
        lines.push(format!("  PAF n={}", hits.len()));
        for (index, hit) in hits.iter().enumerate() {
            let mut line = format!(
                " PAF i={} q={} qs={} qe={} t={} ts={} te={} rev={} bl={} ml={}",
                index,
                hit.query_identifier,
                hit.query_start,
                hit.query_end,
                hit.target_identifier,
                hit.target_start,
                hit.target_end,
                u32::from(hit.reverse),
                hit.block_length,
                hit.matching_length
            );
            // `el` is only written when the exact cigar is supplied; the
            // reference leaves it untouched for the reverse list.
            if with_exact {
                line.push_str(&format!(" el={}", u32::from(hit.exact_region)));
            }
            line.push_str(&format!(" wli={}", u32::from(hit.without_large_indel)));
            lines.push(line);
        }
    }

    fn run_check_batch(
        lines: &mut Vec<String>,
        name: &str,
        intervals: &[(u32, u32)],
        length: i64,
        minimum_depth: i64,
        cells: &[(u16, u32)],
    ) {
        let mut cells_vec = Vec::new();
        for &(operation, cell_length) in cells {
            push_plain(&mut cells_vec, operation, cell_length);
        }
        let hits: Vec<PafHit> = intervals
            .iter()
            .map(|&(start, end)| PafHit {
                query_identifier: 0,
                query_start: start,
                query_end: end,
                target_identifier: 1,
                target_start: 0,
                target_end: end,
                matching_length: 0,
                block_length: 0,
                reverse: false,
                without_large_indel: false,
                exact_region: false,
                deleted: false,
            })
            .collect();

        let mut indices = Vec::new();
        let mut fully_corrected = false;
        let mut abnormal = 0u8;
        check_well_cal(
            &cells_vec,
            &mut indices,
            &mut fully_corrected,
            &mut abnormal,
            length,
            minimum_depth,
            &hits,
        );
        lines.push(format!(
            "CWB name={} len={} min={} f_ec={} abn={} idxn={}",
            name,
            length,
            minimum_depth,
            u32::from(fully_corrected),
            abnormal,
            indices.len()
        ));
        for (index, &key) in indices.iter().enumerate() {
            lines.push(format!("  CWI i={} key={:016x}", index, key));
        }
    }

    fn spec(
        y_identifier: u32,
        x_start: u32,
        x_end: u32,
        is_match: u8,
        non_homopolymer_errors: u32,
    ) -> DedupSpec {
        DedupSpec {
            y_identifier,
            x_start,
            x_end,
            is_match,
            non_homopolymer_errors,
        }
    }

    fn generate_lines() -> Vec<String> {
        let mut lines = Vec::new();

        // ---- dedup_chains ----
        run_dedup_batch(
            &mut lines,
            "dedup_mixed",
            &[
                spec(5, 0, 99, 1, 1),
                spec(3, 10, 39, 1, 0),
                spec(5, 0, 119, 1, 3),
                spec(7, 5, 55, 1, 0),
                spec(3, 10, 69, 0, 2),
                spec(3, 0, 19, 1, 0),
            ],
        );
        run_dedup_batch(
            &mut lines,
            "dedup_ties",
            &[
                spec(1, 0, 29, 1, 2),
                spec(1, 30, 59, 1, 2),
                spec(1, 60, 89, 1, 2),
                spec(1, 10, 39, 0, 100),
            ],
        );
        run_dedup_batch(
            &mut lines,
            "dedup_equal",
            &[
                spec(2, 60, 89, 1, 0),
                spec(2, 0, 29, 1, 0),
                spec(2, 90, 119, 1, 0),
                spec(2, 30, 59, 1, 0),
            ],
        );
        let mut many = Vec::new();
        for index in 0..70u32 {
            many.push(spec(9, index * 2, index * 2 + 19, 1, 0));
        }
        run_dedup_batch(&mut lines, "dedup_many", &many);
        run_dedup_batch(&mut lines, "dedup_single", &[spec(4, 0, 49, 1, 0)]);
        run_dedup_batch(&mut lines, "dedup_empty", &[]);

        // ---- push_nec_re / extract_max_exact ----
        let consensus_spec: [(u16, u16, u16, u32); 12] = [
            (0, 0, 0, 10),
            (1, 2, 1, 1),
            (0, 0, 0, 5),
            (2, 0, 3, 2),
            (0, 0, 0, 20),
            (3, 0, 1, 3),
            (0, 0, 0, 7),
            (1, 1, 2, 1),
            (1, 3, 2, 1),
            (0, 0, 0, 100),
            (1, 0, 0, 1),
            (0, 0, 0, 40),
        ];
        let consensus = build_consensus(&consensus_spec);
        let consensus_region = build_consensus_region(&consensus);
        let consensus_cells = push_nec_re(&consensus_region);
        print_decoded_cells(&mut lines, "NEC", &consensus_cells);

        run_extract_batch(
            &mut lines,
            "exact_full",
            &consensus,
            &[WindowSpec {
                x_start: 0,
                x_end: 49,
                y_start: 0,
                y_end: 49,
                cells: vec![(0, 50)],
            }],
        );
        run_extract_batch(
            &mut lines,
            "exact_two_windows",
            &consensus,
            &[
                WindowSpec {
                    x_start: 0,
                    x_end: 9,
                    y_start: 0,
                    y_end: 9,
                    cells: vec![(0, 10)],
                },
                WindowSpec {
                    x_start: 20,
                    x_end: 69,
                    y_start: 20,
                    y_end: 69,
                    cells: vec![(0, 50)],
                },
            ],
        );
        run_extract_batch(
            &mut lines,
            "exact_with_indel",
            &consensus,
            &[WindowSpec {
                x_start: 0,
                x_end: 29,
                y_start: 0,
                y_end: 29,
                cells: vec![(0, 10), (1, 5), (0, 15)],
            }],
        );
        run_extract_batch(
            &mut lines,
            "exact_none",
            &consensus,
            &[WindowSpec {
                x_start: 0,
                x_end: 7,
                y_start: 0,
                y_end: 7,
                cells: vec![(2, 5), (3, 3)],
            }],
        );
        run_extract_batch(
            &mut lines,
            "exact_short_second",
            &consensus,
            &[
                WindowSpec {
                    x_start: 0,
                    x_end: 39,
                    y_start: 0,
                    y_end: 39,
                    cells: vec![(0, 40)],
                },
                WindowSpec {
                    x_start: 50,
                    x_end: 54,
                    y_start: 50,
                    y_end: 54,
                    cells: vec![(0, 5)],
                },
            ],
        );
        run_extract_batch(
            &mut lines,
            "exact_offset",
            &consensus,
            &[WindowSpec {
                x_start: 100,
                x_end: 149,
                y_start: 60,
                y_end: 109,
                cells: vec![(0, 50)],
            }],
        );

        // ---- push_ne_ovlp ----
        let paf_specs = [
            spec(1, 0, 49, 1, 3),
            spec(2, 100, 149, 1, 0),
            spec(3, 20, 69, 2, 1),
            spec(5, 0, 29, 1, 0),
        ];
        run_push_ne_batch(&mut lines, "push_forward", &consensus, &paf_specs, 1, true);
        run_push_ne_batch(&mut lines, "push_reverse", &consensus, &paf_specs, 2, false);

        // ---- check_well_cal ----
        run_check_batch(
            &mut lines,
            "well_full_double",
            &[(0, 100), (0, 100)],
            100,
            2,
            &[(0, 100)],
        );
        run_check_batch(
            &mut lines,
            "well_interior_hole",
            &[(0, 20), (60, 100)],
            100,
            0,
            &[(0, 100)],
        );
        run_check_batch(&mut lines, "well_left_gap", &[(20, 100)], 100, 0, &[(0, 100)]);
        run_check_batch(&mut lines, "well_right_gap", &[(0, 80)], 100, 0, &[(0, 100)]);
        run_check_batch(&mut lines, "well_low_depth", &[(0, 100)], 100, 2, &[(0, 100)]);
        run_check_batch(&mut lines, "well_empty", &[], 10, 2, &[(0, 100)]);
        run_check_batch(
            &mut lines,
            "well_non_match_cell",
            &[(0, 100), (0, 100)],
            100,
            2,
            &[(0, 30), (1, 5), (0, 65)],
        );
        run_check_batch(
            &mut lines,
            "well_short_read",
            &[(0, 5), (0, 5)],
            5,
            2,
            &[(0, 5)],
        );

        lines
    }

    #[test]
    #[ignore = "requires reference/ec_stored_overlap_probe.out from the C++ probe"]
    fn stored_overlap_helpers_match_the_reference() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("reference")
            .join("ec_stored_overlap_probe.out");
        let contents = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("cannot read {}: {}", path.display(), error));
        // The reference prints its whole buffer after `dedup_chains`, so lines
        // past `ol->length` show the stale entries left by the in-place swaps.
        // A `Vec` is truncated to its length, so those lines cannot be
        // reproduced and are filtered out here.
        let golden: Vec<&str> = contents
            .lines()
            .filter(|line| !line.contains("TAIL"))
            .collect();
        let generated = generate_lines();

        assert_eq!(
            generated.len(),
            golden.len(),
            "the probe printed {} replayable lines, the replay generated {}",
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

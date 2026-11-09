//! Gap reassignment, `reassign_gaps` and its helpers in `Correct.cpp`.
//!
//! The window alignment [`crate::correction::fast_cigar`] produces is not
//! optimal inside homopolymer and tandem-repeat contexts: the aligner may
//! place an indel at any equivalent position and report mismatch cells whose
//! bases really agree once the indel sits elsewhere. `reassign_gaps` slides
//! every gap left as long as the target and query bases keep matching, and
//! turns the mismatch cells it can rescue into exact-match cells. The overlap's
//! substitution count shrinks by the number rescued, which is what the
//! corrector uses later.
//!
//! Each window of the overlap is re-emitted into an auxiliary overlap:
//! unaligned windows are copied verbatim, aligned windows are rebuilt gap by
//! gap. When at least one window changed, the auxiliary window list replaces
//! the original one and the overlap coordinates are re-derived from it.
//!
//! The reference leaves `extra_begin`/`extra_end` of the rebuilt windows
//! uninitialized (`kv_pushp` hands out an untouched slot); this port writes
//! zero there, matching the fresh-arena value the reference sees in practice.

use crate::correction::window::{WindowCigar, CELL_LENGTH_MASK};
use crate::correction::window_chain::is_unaligned_window;
use crate::overlap::region::OverlapRegion;
use crate::correction::target::TargetSource as AlignmentTargets;
use crate::reads::read_store::ReadStore;

/// Cigar operation for a matching base, `EAC_M`.
const EXACT_MATCH: u8 = 0;
/// Cigar operation for a mismatching base, `MIS_M`.
const MISMATCH: u8 = 1;
/// Cigar operation for a target base the query does not have, `MOR_YP`.
const MORE_TARGET: u8 = 2;
/// Cigar operation for a query base the target does not have, `MOR_XT`.
const MORE_QUERY: u8 = 3;

/// Where the target sequence comes from, the `tstr`/`rref` pair of the
/// reference.
pub enum TargetSource<'a> {
    /// `tstr` was provided by the caller.
    Provided(&'a [u8]),
    /// `tstr` was NULL: recover the target read from the read store,
    /// honouring the overlap's target strand.
    Reads(&'a ReadStore),
}

/// Packs one cell the way `push_trace` stores it: operation in the top two
/// bits, length in the low fourteen.
fn pack(operation: u8, length: u32) -> u16 {
    (u16::from(operation) << 14) | (length as u16 & CELL_LENGTH_MASK)
}

/// `push_trace`: appends one run to a packed cell buffer, splitting it into
/// chunks of at most `0x3fff` bases.
fn push_trace(cells: &mut Vec<u16>, operation: u8, mut length: u32) {
    let packed_operation = u16::from(operation) << 14;
    while length >= u32::from(CELL_LENGTH_MASK) {
        cells.push(packed_operation + CELL_LENGTH_MASK);
        length -= u32::from(CELL_LENGTH_MASK);
    }
    if length > 0 {
        cells.push(packed_operation + length as u16);
    }
}

/// `push_wcigar`: points the window at a copy of `cells` appended to the
/// shared buffer.
fn push_window_cells(cigar: &mut WindowCigar, index: usize, cells: &[u16]) {
    cigar.windows[index].cigar_start = cigar.cells.len() as u32;
    cigar.windows[index].cigar_length = cells.len() as u32;
    cigar.cells.extend_from_slice(cells);
}

/// `append_cigar`: appends one run to the window at `index`, merging it with
/// the previous run when both carry the same operation.
fn append_cigar(cigar: &mut WindowCigar, index: usize, operation: u8, length: u32) {
    if length == 0 {
        return;
    }
    let mut previous_operation = u16::MAX;
    let mut previous_length = 0u32;
    if cigar.windows[index].cigar_length > 0 {
        let last = *cigar
            .cells
            .last()
            .expect("a window with cells must end at the shared buffer tail");
        previous_operation = last >> 14;
        previous_length = u32::from(last & CELL_LENGTH_MASK);
    }
    let mut merged_length = length;
    if previous_operation == u16::from(operation) {
        merged_length += previous_length;
        cigar.cells.pop();
        cigar.windows[index].cigar_length -= 1;
    }
    push_trace(&mut cigar.cells, operation, merged_length);
    cigar.windows[index].cigar_length =
        cigar.cells.len() as u32 - cigar.windows[index].cigar_start;
}

/// `adjust_gap`: slides one gap cell left to the nearest position the
/// sequences support, rescuing mismatch cells on the way.
///
/// `target_position` and `query_position` are the alignment cursors just past
/// the gap. The window's cells built so far are scanned backwards and
/// `scratch` collects the pieces of the rebuilt tail in reverse order. Returns
/// whether anything moved together with the number of mismatch cells that
/// turned out to be matches; the count is only meaningful when the first
/// element is one.
fn adjust_gap(
    cigar: &mut WindowCigar,
    index: usize,
    target_sequence: &[u8],
    query_sequence: &[u8],
    mut target_position: i64,
    mut query_position: i64,
    gap_operation: u16,
    scratch: &mut Vec<u16>,
) -> (u16, i64) {
    let mut recovered = 0i64;

    if cigar.windows[index].cigar_length == 0 {
        append_cigar(cigar, index, gap_operation as u8, 1);
        return (0, recovered);
    }
    if gap_operation != u16::from(MORE_TARGET) && gap_operation != u16::from(MORE_QUERY) {
        return (0, recovered);
    }
    if gap_operation == u16::from(MORE_TARGET) {
        // The gap consumes a target base, so its query neighbour sits one
        // position earlier.
        query_position -= 1;
    } else {
        target_position -= 1;
    }

    scratch.clear();
    let cell_start = cigar.windows[index].cigar_start as usize;
    let cell_count = i64::from(cigar.windows[index].cigar_length);
    let mut first_length = [0i64; 2];
    let mut force = 0u16;
    let mut cell_index = cell_count - 1;
    while cell_index >= 0 {
        let packed = cigar.cells[cell_start + cell_index as usize];
        let operation = packed >> 14;
        let length = i64::from(packed & CELL_LENGTH_MASK);
        if operation == u16::from(MORE_TARGET) || operation == u16::from(MORE_QUERY) {
            // Stop at the previous gap: the moved gap lands just after it.
            scratch.push(pack(gap_operation as u8, 1));
            scratch.push(pack(operation as u8, length as u32));
            break;
        } else if operation == u16::from(EXACT_MATCH) {
            let mut k = length - 1;
            while k >= 0 {
                if target_sequence[target_position as usize]
                    != query_sequence[query_position as usize]
                {
                    break;
                }
                k -= 1;
                target_position -= 1;
                query_position -= 1;
            }
            first_length[1] = length - k - 1;
            first_length[0] = k + 1;
            if first_length[1] > 0 {
                scratch.push(pack(EXACT_MATCH, first_length[1] as u32));
                force = 1;
            }
            if first_length[0] > 0 {
                scratch.push(pack(gap_operation as u8, 1));
            }
            if first_length[0] > 0 {
                scratch.push(pack(EXACT_MATCH, first_length[0] as u32));
            }
        } else {
            // A mismatch run; the reference does not trust the cigar here, so
            // every agreeing base becomes an exact match and moves the gap one
            // base left.
            let mut k = length - 1;
            let mut first_part = length;
            while k >= 0 {
                if target_sequence[target_position as usize]
                    == query_sequence[query_position as usize]
                {
                    let second_part = first_part - k - 1;
                    first_part = k;
                    if second_part > 0 {
                        scratch.push(pack(MISMATCH, second_part as u32));
                        force = 1;
                    }
                    scratch.push(pack(EXACT_MATCH, 1));
                    recovered += 1;
                }
                k -= 1;
                target_position -= 1;
                query_position -= 1;
            }
            if first_part > 0 {
                scratch.push(pack(MISMATCH, first_part as u32));
            }
            first_length[0] = 0;
        }
        if first_length[0] > 0 {
            break;
        }
        cell_index -= 1;
    }

    if force == 0 {
        // Nothing could be rescued: the gap stays where it was.
        append_cigar(cigar, index, gap_operation as u8, 1);
        return (0, recovered);
    }
    let cell_start = cigar.windows[index].cigar_start as usize;
    if cell_index >= 0 {
        cigar.windows[index].cigar_length = cell_index as u32;
        cigar.cells.truncate(cell_start + cell_index as usize);
    } else {
        cigar.windows[index].cigar_length = 0;
        cigar.cells.truncate(cell_start);
        append_cigar(cigar, index, gap_operation as u8, 1);
    }
    for offset in (0..scratch.len()).rev() {
        let cell = scratch[offset];
        append_cigar(
            cigar,
            index,
            (cell >> 14) as u8,
            u32::from(cell & CELL_LENGTH_MASK),
        );
    }
    (1, recovered)
}

/// `ajust_end_cigar`: turns leading and trailing mismatch cells into
/// query-only cells and shrinks the window on the target side.
fn adjust_end_cigar(cigar: &mut WindowCigar, index: usize) -> u16 {
    let count = cigar.windows[index].cigar_length as usize;
    if count == 0 {
        return 0;
    }
    let start = cigar.windows[index].cigar_start as usize;
    let mut changed = 0u16;
    for offset in 0..count {
        let cell = cigar.cells[start + offset];
        if cell >> 14 != u16::from(MISMATCH) {
            break;
        }
        let length = cell & CELL_LENGTH_MASK;
        cigar.cells[start + offset] = (u16::from(MORE_QUERY) << 14) | length;
        cigar.windows[index].y_start += i32::from(length);
        changed = 1;
    }
    for offset in (0..count).rev() {
        let cell = cigar.cells[start + offset];
        if cell >> 14 != u16::from(MISMATCH) {
            break;
        }
        let length = cell & CELL_LENGTH_MASK;
        cigar.cells[start + offset] = (u16::from(MORE_QUERY) << 14) | length;
        cigar.windows[index].y_end -= i32::from(length);
        changed = 1;
    }
    changed
}

/// `move_wins`: rebuilds one window of `region` into `auxiliary`.
///
/// Returns whether the window moved, together with the number of mismatch
/// cells rescued by [`adjust_gap`].
fn move_windows(
    region: &OverlapRegion,
    window_index: usize,
    auxiliary: &mut OverlapRegion,
    target: &TargetSource<'_>,
    query_sequence: &[u8],
    scratch: &mut Vec<u16>,
) -> (u16, i64) {
    let window = region.window_cigar.windows[window_index];
    if is_unaligned_window(&window) {
        // Unaligned windows are carried over untouched; their cells were never
        // written, so the copied cigar pointer is never dereferenced.
        auxiliary.window_cigar.windows.push(window);
        return (0, 0);
    }

    let aux_index = auxiliary.window_cigar.append_window(
        window.x_start,
        window.x_end,
        window.y_start,
        window.y_end,
        window.error,
        0,
        0,
        0,
    );
    auxiliary.window_cigar.windows[aux_index].cigar_start =
        auxiliary.window_cigar.cells.len() as u32;
    auxiliary.window_cigar.windows[aux_index].cigar_length = 0;

    let source_start = window.cigar_start as usize;
    let source_count = window.cigar_length as usize;
    if window.error == 0 {
        push_window_cells(
            &mut auxiliary.window_cigar,
            aux_index,
            &region.window_cigar.cells[source_start..source_start + source_count],
        );
        return (0, 0);
    }

    // Look for the first gap that follows a match; without one there is
    // nothing to relocate.
    let mut scan = 0usize;
    let mut saw_match = false;
    while scan < source_count {
        let operation = region.window_cigar.cells[source_start + scan] >> 14;
        if operation < u16::from(MORE_TARGET) {
            saw_match = true;
        } else if saw_match {
            break;
        }
        scan += 1;
    }
    if scan >= source_count {
        push_window_cells(
            &mut auxiliary.window_cigar,
            aux_index,
            &region.window_cigar.cells[source_start..source_start + source_count],
        );
        if adjust_end_cigar(&mut auxiliary.window_cigar, aux_index) != 0 {
            return (1, 0);
        }
        return (0, 0);
    }

    let target_sequence: Vec<u8> = match target {
        TargetSource::Provided(sequence) => sequence.to_vec(),
        TargetSource::Reads(reads) => {
            if region.y_strand != 0 {
                reads.read_sequence_reverse_complement(auxiliary.y_identifier as usize)
            } else {
                reads.read_sequence(auxiliary.y_identifier as usize)
            }
        }
    };

    let mut target_position = i64::from(window.y_start);
    let mut query_position = i64::from(window.x_start);
    let mut moved = 0u16;
    let mut recovered_total = 0i64;
    let mut matched = 0u64;
    for scan in 0..source_count {
        let packed = region.window_cigar.cells[source_start + scan];
        let operation = packed >> 14;
        let length = i64::from(packed & CELL_LENGTH_MASK);
        if operation < u16::from(MORE_TARGET) {
            append_cigar(
                &mut auxiliary.window_cigar,
                aux_index,
                operation as u8,
                length as u32,
            );
            target_position += length;
            query_position += length;
            matched = 1;
        } else if matched == 0 {
            // Leading gaps stay put.
            append_cigar(
                &mut auxiliary.window_cigar,
                aux_index,
                operation as u8,
                length as u32,
            );
            if operation == u16::from(MORE_TARGET) {
                target_position += length;
            } else {
                query_position += length;
            }
        } else {
            for _ in 0..length {
                let (window_moved, recovered) = adjust_gap(
                    &mut auxiliary.window_cigar,
                    aux_index,
                    &target_sequence,
                    query_sequence,
                    target_position,
                    query_position,
                    operation,
                    scratch,
                );
                if window_moved != 0 {
                    moved = 1;
                    let error = auxiliary.window_cigar.windows[aux_index].error;
                    auxiliary.window_cigar.windows[aux_index].error =
                        error.wrapping_sub(recovered as i16);
                    recovered_total += recovered;
                }
                if operation == u16::from(MORE_TARGET) {
                    target_position += 1;
                } else {
                    query_position += 1;
                }
            }
        }
    }

    if adjust_end_cigar(&mut auxiliary.window_cigar, aux_index) != 0 {
        moved = 1;
    }
    (moved, recovered_total)
}

/// `reassign_gaps`: rebuilds every window of `region` into `auxiliary`,
/// relocating gaps to their leftmost equivalent positions.
///
/// When any window changed, `region` adopts the rebuilt window list and the
/// overlap coordinates are recomputed from it. The rescued mismatch count is
/// subtracted from `non_homopolymer_errors`. `target_length` is the `tl`
/// argument of the reference and is consulted only when the target sequence is
/// provided directly.
pub fn reassign_gaps(
    region: &mut OverlapRegion,
    auxiliary: &mut OverlapRegion,
    query_sequence: &[u8],
    query_length: i64,
    target: TargetSource<'_>,
    target_length: i64,
    scratch: &mut Vec<u16>,
) {
    if region.non_homopolymer_errors == 0 {
        return;
    }
    auxiliary.window_cigar.clear();
    auxiliary.y_identifier = region.y_identifier;
    auxiliary.y_strand = region.y_strand;
    auxiliary.x_position_start = region.x_position_start;
    auxiliary.x_position_end = region.x_position_end;
    auxiliary.y_position_start = region.y_position_start;
    auxiliary.y_position_end = region.y_position_end;

    let window_count = region.window_cigar.windows.len();
    let mut moved = 0u16;
    let mut recovered = 0i64;
    for window_index in 0..window_count {
        let (window_moved, window_recovered) = move_windows(
            region,
            window_index,
            auxiliary,
            &target,
            query_sequence,
            scratch,
        );
        if window_moved != 0 {
            moved = 1;
        }
        recovered += window_recovered;
    }
    if moved != 0 {
        let target_total_length = match &target {
            TargetSource::Reads(reads) => {
                reads.read_length(region.y_identifier as usize) as i64
            }
            TargetSource::Provided(_) => target_length,
        };
        region.update_from_windows(&auxiliary.window_cigar, query_length, target_total_length);
    }
    region.non_homopolymer_errors =
        region.non_homopolymer_errors.wrapping_sub(recovered as u32);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn push_trace_splits_runs_at_the_length_field() {
        let mut cells = Vec::new();
        push_trace(&mut cells, EXACT_MATCH, 0x4000);
        assert_eq!(cells, vec![0x3fff, 1]);

        cells.clear();
        push_trace(&mut cells, MORE_QUERY, 0);
        assert!(cells.is_empty());

        cells.clear();
        push_trace(&mut cells, MORE_QUERY, 0x3fff);
        assert_eq!(cells, vec![(3 << 14) | 0x3fff]);
    }

    #[test]
    fn append_cigar_merges_equal_operations() {
        let mut cigar = WindowCigar::default();
        let index = cigar.append_window(0, 9, 0, 9, 0, 0, 0, 0);
        append_cigar(&mut cigar, index, EXACT_MATCH, 4);
        append_cigar(&mut cigar, index, EXACT_MATCH, 3);
        append_cigar(&mut cigar, index, MISMATCH, 1);
        append_cigar(&mut cigar, index, EXACT_MATCH, 2);

        assert_eq!(
            cigar.window_cells(index).collect::<Vec<_>>(),
            vec![(EXACT_MATCH, 7), (MISMATCH, 1), (EXACT_MATCH, 2)]
        );
        assert_eq!(cigar.windows[index].cigar_length, 3);
    }

    #[test]
    fn adjust_end_cigar_converts_mismatch_ends() {
        let mut cigar = WindowCigar::default();
        let index = cigar.append_window(0, 9, 0, 9, 0, 0, 0, 0);
        append_cigar(&mut cigar, index, MISMATCH, 2);
        append_cigar(&mut cigar, index, EXACT_MATCH, 3);
        append_cigar(&mut cigar, index, MISMATCH, 1);

        assert_eq!(adjust_end_cigar(&mut cigar, index), 1);
        assert_eq!(
            cigar.window_cells(index).collect::<Vec<_>>(),
            vec![(MORE_QUERY, 2), (EXACT_MATCH, 3), (MORE_QUERY, 1)]
        );
        assert_eq!(cigar.windows[index].y_start, 2);
        assert_eq!(cigar.windows[index].y_end, 8);
    }
}

#[cfg(test)]
mod probe_goldens {
    use super::*;
    use crate::correction::fast_cigar::generate_fast_cigar;
    use crate::correction::hc_align::BitExtz;
    use crate::overlap::region::{CandidatesList, KmerHit};

    struct Lcg(u64);

    impl Lcg {
        fn next_u32(&mut self) -> u32 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 33) as u32
        }

        fn sequence(&mut self, length: usize) -> Vec<u8> {
            const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
            (0..length)
                .map(|_| BASES[(self.next_u32() & 3) as usize])
                .collect()
        }
    }

    fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                other => *other,
            })
            .collect()
    }

    /// `print_windows`: the reference zeroes the extra boundaries of every
    /// aligned window before printing, and that zeroing persists.
    fn push_windows(lines: &mut Vec<String>, tag: &str, cigar: &mut WindowCigar) {
        lines.push(format!(
            "{} n={} cn={}",
            tag,
            cigar.windows.len(),
            cigar.cells.len()
        ));
        for index in 0..cigar.windows.len() {
            if cigar.windows[index].error != i16::MAX {
                cigar.windows[index].extra_begin = 0;
                cigar.windows[index].extra_end = 0;
            }
            let window = cigar.windows[index];
            let mut line = format!(
                " {}W {} x=[{},{}] y=[{},{}] eb={} ee={} err={} thre={} cidx={} clen={} :",
                tag,
                index,
                window.x_start,
                window.x_end,
                window.y_start,
                window.y_end,
                window.extra_begin,
                window.extra_end,
                window.error,
                window.error_threshold,
                window.cigar_start,
                window.cigar_length
            );
            for offset in 0..window.cigar_length {
                let cell = cigar.cells[(window.cigar_start + offset) as usize];
                line.push_str(&format!(
                    " {}{}",
                    cell & 0x3fff,
                    b"EMDI"[(cell >> 14) as usize] as char
                ));
            }
            lines.push(line);
        }
    }

    #[derive(Clone, Copy)]
    struct WindowSpec {
        x_start: i32,
        x_end: i32,
        y_start: i32,
        y_end: i32,
        error: i16,
        threshold: i16,
    }

    struct ReassignSpec {
        name: &'static str,
        non_homopolymer_errors: u32,
        strand: u32,
        window_length: i64,
        hits: Vec<(u32, u32)>,
    }

    /// The `run_reassign` scenarios: a real alignment driven by
    /// `generate_fast_cigar`, followed by gap reassignment.
    fn run_reassign(
        sequences: &[Vec<u8>],
        x_end: u32,
        y_end: u32,
        windows: &[WindowSpec],
        spec: &ReassignSpec,
        lines: &mut Vec<String>,
    ) {
        let mut reads = ReadStore::default();
        for (index, sequence) in sequences.iter().enumerate() {
            reads.push_read(&format!("read{}", index), sequence);
        }

        let mut region = OverlapRegion {
            x_identifier: 0,
            x_position_start: 0,
            x_position_end: x_end,
            y_identifier: 1,
            y_position_start: 0,
            y_position_end: y_end,
            y_strand: spec.strand,
            overlap_length: x_end + 1,
            non_homopolymer_errors: spec.non_homopolymer_errors,
            ..OverlapRegion::default()
        };
        region.fake_cigar.push(0, 0);
        region.fake_cigar.push(x_end, 0);
        for window in windows {
            region.window_cigar.append_window(
                window.x_start,
                window.x_end,
                window.y_start,
                window.y_end,
                window.error,
                0,
                0,
                window.threshold,
            );
        }

        let mut extz = BitExtz::new();
        let mut auxiliary = OverlapRegion::default();
        let mut candidates = CandidatesList::default();
        for &(self_offset, offset) in &spec.hits {
            candidates.push(KmerHit {
                read_identifier: 1,
                reverse_strand: false,
                offset,
                self_offset,
                count: 1,
            });
        }

        let query = &sequences[0];
        let query_length = query.len() as i64;
        let (_, error) = generate_fast_cigar(
            &mut region,
            &mut candidates,
            spec.window_length,
            &AlignmentTargets::Reads(&reads),
            query,
            &mut extz,
            &mut auxiliary,
            0.04,
            query_length,
            10,
        );
        let reported_error = error.unwrap_or(-1);

        let mut scratch = Vec::new();
        reassign_gaps(
            &mut region,
            &mut auxiliary,
            query,
            query_length,
            TargetSource::Reads(&reads),
            -1,
            &mut scratch,
        );

        lines.push(format!(
            "RA name={} nhpe={} x=[{},{}] y=[{},{}] strand={} re={}",
            spec.name,
            region.non_homopolymer_errors,
            region.x_position_start,
            region.x_position_end,
            region.y_position_start,
            region.y_position_end,
            region.y_strand,
            reported_error
        ));
        push_windows(lines, "RAW", &mut region.window_cigar);
        lines.push(format!(
            "RAA name={} aux_n={} aux_cn={}",
            spec.name,
            auxiliary.window_cigar.windows.len(),
            auxiliary.window_cigar.cells.len()
        ));
    }

    struct ManualWindowSpec {
        x_start: i32,
        x_end: i32,
        y_start: i32,
        y_end: i32,
        error: i16,
        threshold: i16,
        extra_begin: i16,
        extra_end: i16,
        cells: Vec<(u8, u32)>,
    }

    /// The hand-built cigar scenarios: no alignment driver, the window cigars
    /// are written down directly so every `adjust_gap` branch is exercised.
    fn run_manual(
        name: &'static str,
        sequences: &[Vec<u8>],
        x_end: u32,
        y_end: u32,
        strand: u32,
        non_homopolymer_errors: u32,
        windows: &[ManualWindowSpec],
        lines: &mut Vec<String>,
    ) {
        let mut reads = ReadStore::default();
        for (index, sequence) in sequences.iter().enumerate() {
            reads.push_read(&format!("read{}", index), sequence);
        }

        let mut region = OverlapRegion {
            x_identifier: 0,
            x_position_start: 0,
            x_position_end: x_end,
            y_identifier: 1,
            y_position_start: 0,
            y_position_end: y_end,
            y_strand: strand,
            non_homopolymer_errors,
            ..OverlapRegion::default()
        };
        region.fake_cigar.push(0, 0);
        region.fake_cigar.push(x_end, 0);
        for window in windows {
            let index = region.window_cigar.append_window(
                window.x_start,
                window.x_end,
                window.y_start,
                window.y_end,
                window.error,
                window.extra_begin,
                window.extra_end,
                window.threshold,
            );
            let start = region.window_cigar.cells.len();
            for &(operation, length) in &window.cells {
                push_trace(&mut region.window_cigar.cells, operation, length);
            }
            region.window_cigar.assign_cells(index, start);
        }

        let mut auxiliary = OverlapRegion::default();
        let mut scratch = Vec::new();
        let query = &sequences[0];
        let query_length = query.len() as i64;
        reassign_gaps(
            &mut region,
            &mut auxiliary,
            query,
            query_length,
            TargetSource::Reads(&reads),
            -1,
            &mut scratch,
        );

        lines.push(format!(
            "RA name={} nhpe={} x=[{},{}] y=[{},{}] strand={} re=0",
            name,
            region.non_homopolymer_errors,
            region.x_position_start,
            region.x_position_end,
            region.y_position_start,
            region.y_position_end,
            region.y_strand
        ));
        push_windows(lines, "RAW", &mut region.window_cigar);
        lines.push(format!(
            "RAA name={} aux_n={} aux_cn={}",
            name,
            auxiliary.window_cigar.windows.len(),
            auxiliary.window_cigar.cells.len()
        ));
    }

    fn generate_lines() -> Vec<String> {
        let mut lines = Vec::new();
        let mut generator = Lcg(987654321);
        let query = generator.sequence(3000);

        let exact = query.clone();
        let mut substituted = query.clone();
        let mut position = 0;
        while position < substituted.len() {
            let base = &mut substituted[position];
            *base = if *base == b'A' { b'C' } else { b'A' };
            position += 97;
        }
        let mut indel = query.clone();
        indel.splice(1000..1000, b"AC".iter().copied());
        indel.drain(1700..1703);
        let reversed = reverse_complement(&substituted);
        let mut fragmented = query.clone();
        for base in &mut fragmented[1200..2000] {
            *base = if *base == b'A' { b'C' } else { b'A' };
        }

        // Insertions inside low-complexity contexts: the alignment shifts the
        // gap and leaves mismatch cells whose bases actually agree once the
        // gap moves.
        let mut homopolymer_query = query.clone();
        for base in &mut homopolymer_query[995..1005] {
            *base = b'A';
        }
        let mut homopolymer_insert = homopolymer_query.clone();
        homopolymer_insert.insert(1000, b'A');
        let mut homopolymer_delete = homopolymer_query.clone();
        homopolymer_delete.remove(1000);
        let mut tandem_query = query.clone();
        for (offset, base) in tandem_query[992..1008].iter_mut().enumerate() {
            *base = b"ACGT"[offset & 3];
        }
        let mut tandem_insert = tandem_query.clone();
        tandem_insert.splice(1000..1000, b"ACGT".iter().copied());

        let mut windows = Vec::new();
        let mut start = 0i32;
        while start < 3000 {
            let end = if start + 399 < 2999 {
                start + 399
            } else {
                2999
            };
            windows.push(WindowSpec {
                x_start: start,
                x_end: end,
                y_start: start,
                y_end: end,
                error: 4,
                threshold: 0,
            });
            start += 400;
        }
        let mut windows_y3000 = windows.clone();
        for window in &mut windows_y3000 {
            window.y_end = window.x_end + 1;
        }
        let windows_y2999 = windows.clone();

        let mut hits = Vec::new();
        let mut x = 0u32;
        while x < 3000 {
            hits.push((x, x));
            x += 100;
        }

        let spec = |name: &'static str, non_homopolymer_errors: u32, strand: u32| ReassignSpec {
            name,
            non_homopolymer_errors,
            strand,
            window_length: 400,
            hits: hits.clone(),
        };

        run_reassign(
            &[query.clone(), exact.clone()],
            2999,
            2999,
            &windows,
            &spec("exact", 16, 0),
            &mut lines,
        );
        run_reassign(
            &[query.clone(), substituted],
            2999,
            2999,
            &windows,
            &spec("subs", 16, 0),
            &mut lines,
        );
        run_reassign(
            &[query.clone(), indel],
            2999,
            2999,
            &windows,
            &spec("indel", 16, 0),
            &mut lines,
        );
        run_reassign(
            &[query.clone(), reversed],
            2999,
            2999,
            &windows,
            &spec("reverse", 16, 1),
            &mut lines,
        );
        run_reassign(
            &[query.clone(), fragmented],
            2999,
            2999,
            &windows,
            &spec("fragmented", 16, 0),
            &mut lines,
        );
        run_reassign(
            &[query.clone(), exact],
            2999,
            2999,
            &windows,
            &spec("no_errors", 0, 0),
            &mut lines,
        );
        run_reassign(
            &[homopolymer_query.clone(), homopolymer_insert],
            2999,
            3000,
            &windows_y3000,
            &spec("homopolymer_insert", 16, 0),
            &mut lines,
        );
        run_reassign(
            &[homopolymer_query, homopolymer_delete],
            2999,
            2998,
            &windows_y2999,
            &spec("homopolymer_delete", 16, 0),
            &mut lines,
        );
        run_reassign(
            &[tandem_query, tandem_insert],
            2999,
            3003,
            &windows_y3000,
            &spec("tandem_insert", 16, 0),
            &mut lines,
        );

        // -------------------------------------------------- hand-built cigars
        let poly_a = vec![b'A'; 50];

        let mut break_sequence = vec![b'A'; 50];
        for base in &mut break_sequence[3..=8] {
            *base = b'C';
        }

        let mut insert_sequence = vec![b'A'; 50];
        insert_sequence[5] = b'C';
        insert_sequence[6] = b'C';

        let mut no_move_sequence = vec![b'A'; 50];
        no_move_sequence[3] = b'G';
        no_move_sequence[4] = b'A';
        no_move_sequence[5] = b'C';

        // 10E 5M 1D 3E against identical poly-A reads: every "mismatch" is a
        // true match, so adjust_gap recovers them and the D moves to the
        // window start.
        run_manual(
            "manual_match_recovery",
            &[poly_a.clone(), poly_a.clone()],
            49,
            49,
            0,
            16,
            &[ManualWindowSpec {
                x_start: 0,
                x_end: 17,
                y_start: 0,
                y_end: 18,
                error: 5,
                threshold: 5,
                extra_begin: 0,
                extra_end: 0,
                cells: vec![(0, 10), (1, 5), (2, 1), (0, 3)],
            }],
            &mut lines,
        );

        // 5E 3M 1D 2E: the mismatch cells agree (recovered), then the scan
        // stops on a true mismatch inside the leading E run (the `ci >= 0`
        // truncation).
        run_manual(
            "manual_mismatch_break",
            &[break_sequence.clone(), break_sequence.clone()],
            49,
            49,
            0,
            16,
            &[ManualWindowSpec {
                x_start: 0,
                x_end: 9,
                y_start: 0,
                y_end: 10,
                error: 5,
                threshold: 5,
                extra_begin: 0,
                extra_end: 0,
                cells: vec![(0, 5), (1, 3), (2, 1), (0, 2)],
            }],
            &mut lines,
        );

        // 4E 2M 1I 3E: an insertion cell moves the other way (`op0 == 3`).
        run_manual(
            "manual_insertion",
            &[insert_sequence.clone(), insert_sequence.clone()],
            49,
            49,
            0,
            16,
            &[ManualWindowSpec {
                x_start: 0,
                x_end: 9,
                y_start: 0,
                y_end: 8,
                error: 5,
                threshold: 5,
                extra_begin: 0,
                extra_end: 0,
                cells: vec![(0, 4), (1, 2), (3, 1), (0, 3)],
            }],
            &mut lines,
        );

        // No rescued match at all: adjust_gap appends the gap unchanged.
        run_manual(
            "manual_no_move",
            &[no_move_sequence.clone(), no_move_sequence.clone()],
            49,
            49,
            0,
            16,
            &[ManualWindowSpec {
                x_start: 0,
                x_end: 9,
                y_start: 0,
                y_end: 10,
                error: 4,
                threshold: 4,
                extra_begin: 0,
                extra_end: 0,
                cells: vec![(0, 3), (1, 2), (2, 1), (0, 2)],
            }],
            &mut lines,
        );

        // Unaligned window followed by an exact one: no move, both are copied.
        let poly_a_long = vec![b'A'; 250];
        run_manual(
            "manual_unaligned",
            &[poly_a_long.clone(), poly_a_long.clone()],
            199,
            199,
            0,
            16,
            &[
                ManualWindowSpec {
                    x_start: 0,
                    x_end: 99,
                    y_start: 0,
                    y_end: 99,
                    error: i16::MAX,
                    threshold: 1,
                    extra_begin: 0,
                    extra_end: -1,
                    cells: Vec::new(),
                },
                ManualWindowSpec {
                    x_start: 100,
                    x_end: 199,
                    y_start: 100,
                    y_end: 199,
                    error: 0,
                    threshold: 0,
                    extra_begin: 0,
                    extra_end: 0,
                    cells: vec![(0, 100)],
                },
            ],
            &mut lines,
        );

        // Unaligned window plus a window that moves: update_overlap_region runs
        // with the unaligned window first in the list.
        run_manual(
            "manual_unaligned_move",
            &[poly_a_long.clone(), poly_a_long],
            199,
            199,
            0,
            16,
            &[
                ManualWindowSpec {
                    x_start: 0,
                    x_end: 99,
                    y_start: 0,
                    y_end: 99,
                    error: i16::MAX,
                    threshold: 1,
                    extra_begin: 0,
                    extra_end: -1,
                    cells: Vec::new(),
                },
                ManualWindowSpec {
                    x_start: 100,
                    x_end: 117,
                    y_start: 100,
                    y_end: 118,
                    error: 5,
                    threshold: 5,
                    extra_begin: 0,
                    extra_end: 0,
                    cells: vec![(0, 10), (1, 5), (2, 1), (0, 3)],
                },
            ],
            &mut lines,
        );

        lines
    }

    #[test]
    #[ignore = "requires reference/ec_reassign_probe.out from the C++ probe"]
    fn gap_reassignment_matches_the_reference() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("reference")
            .join("ec_reassign_probe.out");
        let contents = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("cannot read {}: {}", path.display(), error));
        let golden: Vec<&str> = contents.lines().collect();
        let generated = generate_lines();

        assert_eq!(
            generated.len(),
            golden.len(),
            "the probe printed {} lines, the replay generated {}",
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

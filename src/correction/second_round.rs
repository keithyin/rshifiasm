//! Second-round boundary correction, the counterpart of `process_boundary`,
//! `merge_cigars` and `add_cigar_to_cigar` in `Correct.cpp`.
//!
//! After a window has been corrected, the reference re-aligns a boundary
//! region straddling the previous window against the same overlaps, then
//! merges that temporary cigar back into the accumulated cigar. Only the
//! middle of the region (25 bases away from either edge) is allowed to change.

use crate::alignment::banded_bpm::NO_ALIGNMENT_ERROR;
use crate::alignment::banded_path::{BandedPath, reserve_banded_bpm_path};
use crate::correction::cigar::{CigarRecord, OPERATION_DELETION, OPERATION_INSERTION};
use crate::correction::cigar_generator::generate_window_cigar;
use crate::correction::consensus::CorrectedRead;
use crate::correction::consensus_graph::{GraphConsensusContext, get_sequence_from_graph};
use crate::correction::dag_con::DagCon;
use crate::correction::interval::IntervalBuffer;
use crate::correction::poa::PoaGraph;
use crate::correction::recalculate::double_error_threshold;
use crate::correction::region::{
    THRESHOLD_MAXIMUM, WINDOW_LENGTH, adjust_threshold, determine_overlap_region, fill_subregion,
};
use crate::overlap::region::OverlapRegion;
use crate::correction::target::TargetSource;
use crate::reads::read_store::ReadStore;
use crate::sequence::encoding::SEQUENCE_NT6_TABLE;

/// `WINDOW_BOUNDARY`: width of the straddling region.
pub const WINDOW_BOUNDARY: i64 = 375;

/// `WINDOW_UNCORRECT_SINGLE_SIDE_BOUNDARY`: bases at each edge that the second
/// round is not allowed to change.
pub const WINDOW_UNCORRECT_SINGLE_SIDE_BOUNDARY: i64 = 25;

/// `MIN_COVERAGE_THRESHOLD` of the second round.
pub const MINIMUM_SECOND_ROUND_COVERAGE: i64 = 3;

/// `Round2_alignment`: the temporary corrected read and the two cigars used
/// by the boundary pass.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SecondRoundAlignment {
    /// Corrected bases of the boundary region.
    pub temporary_read: CorrectedRead,
    /// Temporary cigar of the boundary region.
    pub temporary_cigar: CigarRecord,
    /// Accumulated cigar of the whole read.
    pub cigar: CigarRecord,
    /// Length already represented by `cigar`, `obtained_cigar_length`.
    pub obtained_cigar_length: i64,
}

impl SecondRoundAlignment {
    /// `clear_Round2_alignment`.
    pub fn clear(&mut self) {
        self.temporary_read = CorrectedRead::new();
        self.temporary_cigar = CigarRecord::new();
        self.cigar = CigarRecord::new();
        self.obtained_cigar_length = 0;
    }
}

/// `add_cigar_to_cigar`: merges the trimmed temporary cigar into the
/// accumulated cigar, converting substitutions into packed base pairs.
#[allow(clippy::too_many_arguments)]
fn add_cigar_to_cigar(
    backbone_read: &[u8],
    second_round: &mut SecondRoundAlignment,
    backbone_start: i64,
    backbone_length: i64,
    new_start: i64,
) {
    let x_offset = backbone_start.max(0) as usize;
    let y_offset = new_start.max(0) as usize;

    // Everything between the last recorded position and this boundary is an
    // unchanged match.
    let leading = backbone_start - second_round.obtained_cigar_length;
    let leading_bases = vec![b'A'; leading.max(0) as usize];
    second_round
        .cigar
        .add(&leading_bases, leading.max(0) as u32, 0);
    second_round.obtained_cigar_length = backbone_start + backbone_length;

    let mut x_index = 0usize;
    let mut y_index = 0usize;
    let cells: Vec<(u32, u32)> = second_round.temporary_cigar.cells().collect();
    for (operation, length) in cells {
        match operation {
            0 => {
                let bases = vec![b'A'; length as usize];
                second_round.cigar.add(&bases, length, 0);
                x_index += length as usize;
                y_index += length as usize;
            }
            1 => {
                for _ in 0..length {
                    let x_base = backbone_read.get(x_offset + x_index).copied().unwrap_or(b'N');
                    let y_base = second_round
                        .temporary_read
                        .sequence
                        .get(y_offset + y_index)
                        .copied()
                        .unwrap_or(b'N');
                    let merged = (SEQUENCE_NT6_TABLE[y_base as usize] << 3)
                        | SEQUENCE_NT6_TABLE[x_base as usize];
                    second_round.cigar.add(&[merged], 1, 1);
                    x_index += 1;
                    y_index += 1;
                }
            }
            2 => {
                let start = y_offset + y_index;
                let bases: Vec<u8> = second_round.temporary_read.sequence
                    [start..(start + length as usize).min(second_round.temporary_read.sequence.len())]
                    .to_vec();
                second_round.cigar.add(&bases, length, OPERATION_INSERTION);
                y_index += length as usize;
            }
            3 => {
                let start = x_offset + x_index;
                let bases: Vec<u8> = backbone_read
                    [start..(start + length as usize).min(backbone_read.len())]
                    .to_vec();
                second_round.cigar.add(&bases, length, OPERATION_DELETION);
                x_index += length as usize;
            }
            other => panic!("unknown cigar operation {other}"),
        }
    }
}

/// `merge_cigars`: trims the temporary cigar to the safe middle of the
/// boundary region and merges it into the accumulated cigar.
fn merge_cigars(
    backbone_read: &[u8],
    second_round: &mut SecondRoundAlignment,
    total_start_base: i64,
    total_end_base: i64,
    total_window_start: i64,
) -> bool {
    let temporary = &mut second_round.temporary_cigar;
    if temporary.cell_count() == 1 && temporary.cell_operation(0) == 0 {
        return true;
    }

    let start_base = total_start_base - total_window_start;
    let end_base = total_end_base - total_window_start;
    let mut x_index = 0i64;
    let mut y_index = 0i64;
    let mut x_start = -1i64;
    let mut x_end = -1i64;
    let mut y_start = -1i64;
    let mut start_cigar = -1i64;
    let mut end_cigar = -1i64;

    let cells = temporary.cells().collect::<Vec<_>>();
    for (index, (operation, length)) in cells.iter().enumerate() {
        match operation {
            0 | 1 => {
                for _ in 0..*length {
                    if x_index >= start_base && x_start == -1 {
                        x_start = x_index;
                        y_start = y_index;
                        start_cigar = index as i64;
                    }
                    if x_index >= end_base && x_end == -1 {
                        x_end = x_index;
                        end_cigar = index as i64;
                        break;
                    }
                    x_index += 1;
                    y_index += 1;
                }
            }
            2 => y_index += *length as i64,
            3 => x_index += *length as i64,
            other => panic!("unknown cigar operation {other}"),
        }
    }

    if x_end == -1 || x_start == -1 {
        return false;
    }

    // Shorten the first and last cells so the merged cigar covers exactly
    // [start_base, end_base].
    x_index = 0;
    // The reference also advances a target index in this loop but never reads
    // it, so only the query index is tracked here.
    y_index = 0;
    let _ = y_index;
    for (index, (operation, length)) in cells.iter().enumerate() {
        let index = index as i64;
        if index == start_cigar {
            let trimmed = (*length as i64 - (x_start - x_index)) as u32;
            temporary.set_cell(index as usize, *operation, trimmed);
            if index == end_cigar {
                let final_length = (x_end - x_start + 1) as u32;
                temporary.set_cell(index as usize, *operation, final_length);
                break;
            }
        } else if index == end_cigar {
            let final_length = (x_end - x_index + 1) as u32;
            temporary.set_cell(index as usize, *operation, final_length);
            break;
        }

        match operation {
            0 | 1 => x_index += *length as i64,
            3 => x_index += *length as i64,
            _ => {}
        }
    }

    let retained = (end_cigar - start_cigar + 1) as usize;
    temporary.retain_cells(start_cigar as usize, retained);

    let total_x_start = total_window_start + x_start;
    let x_length = x_end - x_start + 1;
    let total_y_start = y_start;
    add_cigar_to_cigar(
        backbone_read,
        second_round,
        total_x_start,
        x_length,
        total_y_start,
    );
    true
}

/// One alignment attempt for a boundary overlap.
struct BoundaryAlignment {
    target_start: i64,
    extra_begin: i32,
    extra_end: i32,
    threshold: i64,
    path: Vec<u8>,
    start_site: i32,
    end_site: i32,
    error: u32,
    target_sequence: Vec<u8>,
}

fn align_boundary_overlap(
    reads: &ReadStore,
    region: &OverlapRegion,
    anchor_target_start: i64,
    query: &[u8],
    query_length: i64,
    threshold: i64,
) -> Option<BoundaryAlignment> {
    let target_read_length = reads.read_length(region.y_identifier as usize) as i64;
    let target_start = anchor_target_start - WINDOW_BOUNDARY / 2;
    if target_start < 0 {
        return None;
    }
    let target_window_length = query_length + (threshold << 1);

    let mapped = determine_overlap_region(
        threshold,
        target_start,
        target_read_length,
        target_window_length,
    )?;
    let target_sequence = fill_subregion(
        &TargetSource::Reads(reads),
        region.y_identifier as usize,
        mapped.start,
        mapped.length,
        region.y_strand as u8,
        mapped.extra_begin,
        mapped.extra_end,
    );

    let BandedPath {
        end_site,
        start_site,
        error,
        path,
    } = reserve_banded_bpm_path(
        &target_sequence,
        target_window_length as usize,
        query,
        query_length as usize,
        threshold as usize,
        -1,
        -1,
    )?;
    if error == NO_ALIGNMENT_ERROR {
        return None;
    }

    Some(BoundaryAlignment {
        target_start: mapped.start,
        extra_begin: mapped.extra_begin,
        extra_end: mapped.extra_end,
        threshold,
        path,
        start_site,
        end_site,
        error,
        target_sequence,
    })
}

/// `process_boundary`: re-aligns a region straddling the previous window
/// boundary and merges the result into the second-round cigar.
#[allow(clippy::too_many_arguments)]
pub fn process_boundary(
    backbone_read: &[u8],
    overlaps: &[OverlapRegion],
    buffer: &IntervalBuffer,
    reads: &ReadStore,
    uncorrected_window_start: i64,
    last_boundary_length: i64,
    maximum_overlap_divergence: f64,
    backbone_graph: &mut PoaGraph,
    dag: &mut DagCon,
    second_round: &mut SecondRoundAlignment,
) -> bool {
    if last_boundary_length == 0 {
        return false;
    }
    let total_length = backbone_read.len() as i64;
    let corrected_window_start = (last_boundary_length - WINDOW_BOUNDARY / 2).max(0);
    let corrected_window_end =
        (last_boundary_length + WINDOW_BOUNDARY / 2 - 1).min(total_length - 1);
    if corrected_window_end < corrected_window_start {
        return false;
    }

    backbone_graph.clear();
    dag.graph.clear();
    dag.node_queue.clear();

    let backbone = &backbone_read[corrected_window_start as usize..=corrected_window_end as usize];
    backbone_graph.add_unmatched_sequence(backbone);

    let query_length = corrected_window_end - corrected_window_start + 1;
    let mut matched_coverage = 0i64;
    let mut total_error = 0i64;

    for identifier in buffer.full_windows() {
        let region = &overlaps[*identifier as usize];
        let aligned_origin =
            (region.x_position_start as i64 / WINDOW_LENGTH) * WINDOW_LENGTH;
        let window_id = (uncorrected_window_start - aligned_origin) / WINDOW_LENGTH;
        if window_id < 0 || window_id as usize >= region.window_cigar.windows.len() {
            continue;
        }
        let window = region.window_cigar.windows[window_id as usize];
        if window.y_end == -1 || window.x_start as i64 != uncorrected_window_start {
            continue;
        }

        let mut threshold = (query_length as f64 * maximum_overlap_divergence) as i64;
        threshold = adjust_threshold(threshold, query_length);
        let mut alignment = align_boundary_overlap(
            reads,
            region,
            window.y_start as i64,
            backbone,
            query_length,
            threshold,
        );

        if alignment.is_none() {
            // The reference retries with a doubled budget.
            let doubled = double_error_threshold(threshold, query_length)
                .min(THRESHOLD_MAXIMUM);
            alignment = align_boundary_overlap(
                reads,
                region,
                window.y_start as i64,
                backbone,
                query_length,
                doubled,
            );
        }

        let Some(mut alignment) = alignment else { continue };
        total_error += alignment.error as i64;
        matched_coverage += 1;

        // Rebuild the target subsequence from the traceback coordinates and
        // fold it into the boundary backbone graph.
        let mut start_site = alignment.start_site - alignment.extra_begin;
        let mut end_site = alignment.end_site - alignment.extra_begin;
        let target_length = end_site - start_site + 1;
        let target_start = alignment.target_start + start_site as i64;
        let target = reads.read_sub_region(
            region.y_identifier as usize,
            target_start,
            target_length as i64,
            region.y_strand as u8,
        );

        let mut temporary_cigar = crate::correction::window::WindowCigar::default();
        temporary_cigar.append_window(
            corrected_window_start as i32,
            corrected_window_end as i32,
            alignment.target_start as i32,
            (alignment.target_start + alignment.end_site as i64) as i32,
            alignment.error as i16,
            alignment.extra_begin as i16,
            alignment.extra_end as i16,
            alignment.threshold as i16,
        );
        generate_window_cigar(
            &mut alignment.path,
            0,
            &mut temporary_cigar,
            &mut start_site,
            &mut end_site,
            &mut alignment.error,
            backbone,
            &alignment.target_sequence,
        );
        let cells: Vec<(u32, u32)> = temporary_cigar
            .window_cells(0)
            .map(|(operation, length)| (operation as u32, length))
            .collect();
        backbone_graph.add_matched_sequence(0, &target, &cells);
    }

    if matched_coverage < MINIMUM_SECOND_ROUND_COVERAGE || total_error == 0 {
        return false;
    }

    second_round.temporary_cigar = CigarRecord::new();
    second_round.temporary_read = CorrectedRead::new();
    let start_base = corrected_window_start + WINDOW_UNCORRECT_SINGLE_SIDE_BOUNDARY;
    let end_base = corrected_window_end - WINDOW_UNCORRECT_SINGLE_SIDE_BOUNDARY;
    if end_base <= start_base {
        return false;
    }

    get_sequence_from_graph(
        backbone_graph,
        dag,
        &mut second_round.temporary_read,
        &GraphConsensusContext {
            query_sequence: backbone_read,
            query_window: backbone,
            window_start: corrected_window_start,
        },
    );
    // `get_sequence_from_graph` writes into `temporary_read.cigar`, which plays
    // the role of `second_round->tmp_cigar`.
    second_round.temporary_cigar = second_round.temporary_read.cigar.clone();

    merge_cigars(
        backbone_read,
        second_round,
        start_base,
        end_base,
        corrected_window_start,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_empty_boundary_is_not_processed() {
        let reads = ReadStore::default();
        let mut graph = PoaGraph::new();
        let mut dag = DagCon::default();
        let mut second_round = SecondRoundAlignment::default();
        let buffer = IntervalBuffer::with_capacity(1);

        assert!(!process_boundary(
            b"ACGT",
            &[],
            &buffer,
            &reads,
            0,
            0,
            0.04,
            &mut graph,
            &mut dag,
            &mut second_round,
        ));
    }

    /// Golden values printed by the real `merge_cigars` for a 50 base window
    /// with the safe middle `[10, 40)`.
    #[test]
    fn merge_cigars_matches_the_reference() {
        struct Case {
            name: &'static str,
            cells: Vec<(u32, u32)>,
            obtained: i64,
            first: (u32, u32),
            result: Vec<(u32, u32)>,
        }
        let cases = [
            Case {
                name: "trim_both",
                cells: vec![(0, 10), (1, 30), (0, 10)],
                obtained: 40,
                first: (1, 30),
                result: vec![(0, 10), (1, 30)],
            },
            Case {
                name: "pure_match",
                cells: vec![(0, 50)],
                obtained: 0,
                first: (0, 50),
                result: vec![],
            },
            Case {
                name: "with_insertion",
                cells: vec![(0, 10), (2, 2), (0, 38)],
                obtained: 40,
                first: (0, 30),
                result: vec![(0, 40)],
            },
        ];

        for case in cases {
            let backbone = vec![b'A'; 100];
            let mut second_round = SecondRoundAlignment::default();
            let mut first_cell = true;
            for (operation, length) in &case.cells {
                let bases = vec![0u8; *length as usize];
                second_round.temporary_cigar.add(&bases, *length, *operation);
                first_cell = false;
            }
            let _ = first_cell;
            second_round.temporary_read.sequence = vec![b'C'; 50];
            second_round.temporary_read.cigar = CigarRecord::new();

            let merged = merge_cigars(&backbone, &mut second_round, 10, 39, 0);
            assert!(merged, "{} must merge", case.name);
            assert_eq!(
                second_round.obtained_cigar_length, case.obtained,
                "obtained length of {}",
                case.name
            );
            assert_eq!(second_round.temporary_cigar.cell_count(), 1);
            assert_eq!(
                (
                    second_round.temporary_cigar.cell_operation(0),
                    second_round.temporary_cigar.cell_length(0)
                ),
                case.first,
                "first cell of {}",
                case.name
            );
            let result: Vec<(u32, u32)> = second_round.cigar.cells().collect();
            assert_eq!(result, case.result, "result cigar of {}", case.name);
        }
    }
}

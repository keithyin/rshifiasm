//! Turning a traceback path into packed window cigar cells, the counterpart of
//! `move_gap_greedy` and `generate_cigar` in `Correct.cpp`.
//!
//! A path is a sequence of operations stored most-recent-first. Before it is
//! run-length encoded, leading and trailing mismatches are trimmed (they would
//! otherwise extend the alignment outside the window) and adjacent gaps are
//! shifted greedily so that the same mismatch is never charged twice.

use crate::correction::window::WindowCigar;

/// `move_gap_greedy`: slides a gap towards the end of the path while the bases
/// it passes still agree, so that the reported cigar has fewer differences.
///
/// `query` and `target` are the two windows being compared; `query_index` and
/// `target_index` are the positions reached by the path so far. The path is
/// modified in place and `error` is decremented when a mismatch turns into a
/// match.
fn move_gap_greedy(
    path: &mut [u8],
    path_index: usize,
    query: &[u8],
    mut query_index: i64,
    target: &[u8],
    mut target_index: i64,
    error: &mut u32,
) {
    let operation = path[path_index];
    if operation < 2 {
        return;
    }

    let mut path_index = path_index + 1;
    if operation == 3 {
        target_index -= 1;
    } else {
        query_index -= 1;
    }

    while path_index < path.len() && query_index >= 0 && target_index >= 0 {
        let next_operation = path[path_index];
        let query_base = query.get(query_index as usize).copied().unwrap_or(b'N');
        let target_base = target.get(target_index as usize).copied().unwrap_or(b'N');

        if next_operation == 2
            || next_operation == 3
            || (next_operation == 0 && query_base != target_base)
        {
            break;
        }

        if next_operation == 1 && query_base == target_base {
            path[path_index - 1] = 0;
            *error = error.wrapping_sub(1);
        } else {
            path[path_index - 1] = next_operation;
        }
        path[path_index] = operation;

        path_index += 1;
        query_index -= 1;
        target_index -= 1;
    }
}

/// `generate_cigar`: writes the packed cells of one window into `window_cigar`.
///
/// `start` and `end` are the target coordinates of the alignment; both are
/// adjusted by leading/trailing mismatch trimming. `error` may also shrink
/// when gaps are moved greedily.
pub fn generate_window_cigar(
    path: &mut [u8],
    window_index: usize,
    window_cigar: &mut WindowCigar,
    start: &mut i32,
    end: &mut i32,
    error: &mut u32,
    query: &[u8],
    target: &[u8],
) {
    let cells_start = window_cigar.cells.len();

    if *error == 0 {
        let window = window_cigar.windows[window_index];
        window_cigar.push_cell(0, (window.x_end + 1 - window.x_start) as u32);
        window_cigar.assign_cells(window_index, cells_start);
        return;
    }

    // Leading mismatches (at the start of the path) become insertions and
    // shorten the alignment end; trailing ones shorten the alignment start.
    let path_length = path.len();
    let mut last_trailing_mismatch: i32 = -1;
    for index in 0..path_length {
        if path[index] == 1 {
            path[index] = 3;
            *end -= 1;
            last_trailing_mismatch = index as i32;
        } else {
            break;
        }
    }
    for index in (0..path_length).rev() {
        if path[index] == 1 {
            path[index] = 3;
            *start += 1;
        } else {
            break;
        }
    }

    let target_offset = (*start).max(0) as usize;
    let mut query_index: i64 = 0;
    let mut target_index: i64 = 0;
    for index in (0..path_length).rev() {
        if index as i32 <= last_trailing_mismatch {
            break;
        }
        match path[index] {
            0 | 1 => {
                query_index += 1;
                target_index += 1;
            }
            2 => {
                move_gap_greedy(
                    path,
                    index,
                    query,
                    query_index,
                    &target[target_offset..],
                    target_index,
                    error,
                );
                target_index += 1;
            }
            3 => {
                move_gap_greedy(
                    path,
                    index,
                    query,
                    query_index,
                    &target[target_offset..],
                    target_index,
                    error,
                );
                query_index += 1;
            }
            other => panic!("unknown traceback operation {other}"),
        }
    }

    // Run-length encode the path from the end towards the front.
    let mut previous_operation: u8 = 5;
    let mut previous_length: u32 = 0;
    for index in (0..path_length).rev() {
        if previous_operation != path[index] {
            if previous_length != 0 {
                window_cigar.push_cell(previous_operation, previous_length);
            }
            previous_operation = path[index];
            previous_length = 1;
        } else {
            previous_length += 1;
        }
    }
    if previous_length != 0 {
        window_cigar.push_cell(previous_operation, previous_length);
    }

    window_cigar.assign_cells(window_index, cells_start);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cigar_with_window(start: i32, end: i32) -> WindowCigar {
        let mut cigar = WindowCigar::default();
        cigar.append_window(start, end, start, end, 0, 0, 0, 2);
        cigar
    }

    #[test]
    fn an_exact_match_becomes_one_cell() {
        let mut cigar = cigar_with_window(0, 9);
        let mut path = [0u8; 10];
        let mut start = 2;
        let mut end = 11;
        let mut error = 0;

        generate_window_cigar(
            &mut path,
            0,
            &mut cigar,
            &mut start,
            &mut end,
            &mut error,
            b"GTACGTACGG",
            b"ACGTACGTACGG",
        );

        assert_eq!(cigar.window_cells(0).collect::<Vec<_>>(), vec![(0, 10)]);
        assert_eq!(start, 2);
        assert_eq!(end, 11);
    }

    #[test]
    fn a_mismatch_path_is_run_length_encoded() {
        let mut cigar = cigar_with_window(0, 9);
        // The real traceback of a single mismatch at the second base.
        let mut path = [0u8; 10];
        path[1] = 1;
        let mut start = 2;
        let mut end = 11;
        let mut error = 1;

        generate_window_cigar(
            &mut path,
            0,
            &mut cigar,
            &mut start,
            &mut end,
            &mut error,
            b"GTACGTACAG",
            b"ACGTACGTACGG",
        );

        let cells: Vec<(u8, u32)> = cigar.window_cells(0).collect();
        assert_eq!(cells.first(), Some(&(0, 8)));
        assert_eq!(cells.last(), Some(&(0, 1)));
        assert!(cells.contains(&(1, 1)), "one mismatch cell: {cells:?}");
    }

    #[test]
    fn leading_mismatches_are_trimmed_into_insertions() {
        let mut cigar = cigar_with_window(0, 9);
        let mut path = [0u8; 10];
        path[0] = 1;
        path[1] = 1;
        let mut start = 2;
        let mut end = 11;
        let mut error = 2;

        generate_window_cigar(
            &mut path,
            0,
            &mut cigar,
            &mut start,
            &mut end,
            &mut error,
            b"GTACGTACAA",
            b"ACGTACGTACGG",
        );

        assert_eq!(end, 9, "two leading mismatches shorten the end");
        assert_eq!(path[0], 3);
        assert_eq!(path[1], 3);
    }

    /// Golden values printed by the real `generate_cigar`.
    #[test]
    fn matches_the_reference_cigar_generation() {
        struct Case {
            name: &'static str,
            path: &'static str,
            query: &'static str,
            target: &'static str,
            start: i32,
            end: i32,
            error: u32,
            cells: Vec<(u8, u32)>,
            resulting_path: &'static str,
        }

        let cases = vec![
            Case {
                name: "exact",
                path: "0000000000",
                query: "GTACGTACGG",
                target: "ACGTACGTACGG",
                start: 2,
                end: 11,
                error: 0,
                cells: vec![(0, 10)],
                resulting_path: "0000000000",
            },
            Case {
                name: "one_mismatch",
                path: "0100000000",
                query: "GTACGTACAG",
                target: "ACGTACGTACGG",
                start: 2,
                end: 11,
                error: 1,
                cells: vec![(0, 8), (1, 1), (0, 1)],
                resulting_path: "0100000000",
            },
            Case {
                name: "two_leading",
                path: "1100000000",
                query: "GTACGTACAA",
                target: "ACGTACGTACGG",
                start: 2,
                end: 11,
                error: 2,
                cells: vec![(0, 8), (3, 2)],
                resulting_path: "3300000000",
            },
            Case {
                name: "middle_mismatch",
                path: "0000000001000000000000",
                query: "CTATCCATAATTAGACCCCCTG",
                target: "GCTATCCATAATTGGACCCCCTGCGATT",
                start: 1,
                end: 22,
                error: 1,
                cells: vec![(0, 12), (1, 1), (0, 9)],
                resulting_path: "0000000001000000000000",
            },
            Case {
                name: "indel",
                path: "300000000",
                query: "TGTCGCCGA",
                target: "TAATGATGTCGCCGTAAAT",
                start: 6,
                end: 13,
                error: 1,
                cells: vec![(0, 8), (3, 1)],
                resulting_path: "300000000",
            },
        ];

        for case in cases {
            let mut path: Vec<u8> = case.path.bytes().map(|code| code - b'0').collect();
            let window_end = path.len() as i32 - 1;
            let mut cigar = cigar_with_window(0, window_end);
            let mut start = case.start;
            let mut end = case.end;
            let mut error = case.error;

            generate_window_cigar(
                &mut path,
                0,
                &mut cigar,
                &mut start,
                &mut end,
                &mut error,
                case.query.as_bytes(),
                case.target.as_bytes(),
            );

            let cells: Vec<(u8, u32)> = cigar.window_cells(0).collect();
            assert_eq!(cells, case.cells, "cells of {}", case.name);
            assert_eq!(start, case.start, "start of {}", case.name);
            // C trims leading mismatches by shortening the alignment end.
            let expected_end = if case.name == "two_leading" { 9 } else { case.end };
            assert_eq!(end, expected_end, "end of {}", case.name);
            let resulting_path: String =
                path.iter().map(|code| char::from(b'0' + *code)).collect();
            assert_eq!(resulting_path, case.resulting_path, "path of {}", case.name);
        }
    }
}

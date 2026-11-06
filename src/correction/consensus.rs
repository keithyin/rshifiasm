//! Building blocks of `generate_consensus`: the corrected read buffer, the
//! strict homopolymer test and the selection of usable overlaps.

use crate::correction::cigar::{CigarRecord, OPERATION_MATCH};
use crate::correction::interval::{
    INTERVAL_HAS_OVERLAP, INTERVAL_NO_OVERLAP, INTERVAL_PAST_THE_END, IntervalBuffer,
    overlap_length,
};
use crate::overlap::region::OverlapRegion;

/// The corrected sequence and its cigar, the parts of `Correct_dumy` and
/// `Cigar_record` that `generate_consensus` writes.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CorrectedRead {
    /// Corrected bases produced so far, `corrected_read`.
    pub sequence: Vec<u8>,
    /// Number of bases changed by correction, `corrected_base`.
    pub corrected_base: u64,
    /// Length of the corrected read at the previous window boundary,
    /// `last_boundary_length`.
    pub last_boundary_length: u64,
    /// Cigar describing how the corrected read differs from the original.
    pub cigar: CigarRecord,
}

impl CorrectedRead {
    /// A fresh buffer with an empty cigar.
    pub fn new() -> Self {
        Self {
            sequence: Vec::new(),
            corrected_base: 0,
            last_boundary_length: 0,
            cigar: CigarRecord::new(),
        }
    }

    /// `add_base_to_correct_read_directly`.
    pub fn add_base_directly(&mut self, base: u8) {
        self.sequence.push(base);
    }

    /// `add_segment_to_correct_read`.
    pub fn add_segment(&mut self, segment: &[u8]) {
        self.sequence.extend_from_slice(segment);
    }

    /// Appends bases and records the corresponding cigar operation, the
    /// combination of `add_segment_to_correct_read` and `add_cigar_record`.
    pub fn add_operation(&mut self, bases: &[u8], operation: u32) {
        self.add_segment(bases);
        self.cigar.add(bases, bases.len() as u32, operation);
    }

    /// Appends one base as a match, `add_cigar_record(..., 0)`.
    pub fn add_match(&mut self, base: u8) {
        self.add_base_directly(base);
        self.cigar.add(&[base], 1, OPERATION_MATCH);
    }
}

/// `if_is_homopolymer_strict`: whether `site` sits in a run of at least three
/// identical bases, looking at most three bases in either direction.
///
/// A lower correction threshold is used at such sites, so getting this test
/// wrong changes which bases get corrected.
pub fn is_strict_homopolymer(site: i64, read: &[u8]) -> bool {
    const THRESHOLD: i64 = 3;
    if site < 0 || site as usize >= read.len() {
        return false;
    }

    let begin = (site - THRESHOLD).max(0);
    let end = (site + THRESHOLD).min(read.len() as i64 - 1);
    let site_base = read[site as usize];

    let mut forward_character = 0u8;
    let mut forward_length = 0i64;
    for index in site + 1..=end {
        let base = read[index as usize];
        if forward_character == 0 {
            forward_character = base;
            forward_length = 1;
        } else if base != forward_character {
            break;
        } else {
            forward_length += 1;
        }
    }

    let mut backward_character = 0u8;
    let mut backward_length = 0i64;
    for index in (begin..site).rev() {
        let base = read[index as usize];
        if backward_character == 0 {
            backward_character = base;
            backward_length = 1;
        } else if base != backward_character {
            break;
        } else {
            backward_length += 1;
        }
    }

    if forward_character == site_base {
        forward_length += 1;
    } else if backward_character == site_base {
        backward_length += 1;
    }

    if forward_length >= THRESHOLD || backward_length >= THRESHOLD {
        return true;
    }

    site_base == forward_character
        && backward_character == forward_character
        && forward_length + backward_length >= THRESHOLD
}

/// `get_available_interval`: collects the overlaps that cover the window and
/// are marked as real matches (`is_match == 1`).
///
/// Unlike `get_interval` this does not split by coverage length: every usable
/// overlap is appended to the front of the identifier array. The return value
/// distinguishes "no overlap touches the window" from "overlaps touch it but
/// none is usable", which the caller treats the same way in practice.
pub fn get_available_interval(
    window_start: i64,
    window_end: i64,
    overlaps: &[OverlapRegion],
    buffer: &mut IntervalBuffer,
) -> i32 {
    if window_start == 0 {
        buffer.start_index = 0;
    }

    let mut index = buffer.start_index as usize;
    while index < overlaps.len() {
        if window_end < overlaps[index].x_position_start as i64 {
            buffer.start_index = 0;
            buffer.length = 0;
            buffer.partial_length = 0;
            return INTERVAL_NO_OVERLAP;
        }
        buffer.start_index = index as u64;
        break;
    }

    if index >= overlaps.len() {
        buffer.start_index = overlaps.len() as u64;
        buffer.length = 0;
        buffer.partial_length = 0;
        return INTERVAL_PAST_THE_END;
    }

    buffer.length = 0;
    buffer.partial_length = 0;
    let mut touched = 0i64;
    let mut found = false;

    while index < overlaps.len() {
        let overlap = &overlaps[index];
        let covered = overlap_length(
            window_start,
            window_end,
            overlap.x_position_start as i64,
            overlap.x_position_end as i64,
        );

        if covered > 0 {
            touched += 1;
            if overlap.is_match == 1 {
                buffer.overlap_identifiers[buffer.length as usize] = index as u64;
                buffer.length += 1;
            }
            if !found {
                found = true;
                buffer.start_index = index as u64;
            }
        }

        if overlap.x_position_start as i64 > window_end {
            break;
        }
        index += 1;
    }

    if touched == 0 {
        INTERVAL_NO_OVERLAP
    } else {
        INTERVAL_HAS_OVERLAP
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_long_run_is_a_homopolymer() {
        assert!(is_strict_homopolymer(3, b"AAAACGT"));
        assert!(is_strict_homopolymer(2, b"ACCCGT"));
        assert!(is_strict_homopolymer(0, b"CCCCGT"));
    }

    #[test]
    fn a_short_or_broken_run_is_not_a_homopolymer() {
        assert!(!is_strict_homopolymer(2, b"AACCG"));
        assert!(!is_strict_homopolymer(3, b"ACGTACGT"));
        assert!(!is_strict_homopolymer(0, b"ACGT"));
        assert!(!is_strict_homopolymer(-1, b"AAAA"));
        assert!(!is_strict_homopolymer(4, b"AAAA"));
    }

    #[test]
    fn corrected_read_tracks_sequence_and_cigar() {
        let mut corrected = CorrectedRead::new();
        corrected.add_operation(b"ACGT", OPERATION_MATCH);
        corrected.add_match(b'A');

        assert_eq!(corrected.sequence, b"ACGTA");
        assert_eq!(corrected.cigar.new_read_length, 5);
        assert_eq!(corrected.cigar.cell_count(), 1);
        corrected.last_boundary_length = corrected.sequence.len() as u64;
        assert_eq!(corrected.last_boundary_length, 5);
    }

    #[test]
    fn available_interval_only_keeps_real_matches() {
        let make = |start: u32, is_match: u8| OverlapRegion {
            x_position_start: start,
            x_position_end: start + 99,
            is_match,
            ..OverlapRegion::default()
        };
        let overlaps = vec![make(0, 1), make(0, 0), make(0, 3), make(50, 1)];
        let mut buffer = IntervalBuffer::with_capacity(8);

        let flag = get_available_interval(0, 99, &overlaps, &mut buffer);
        assert_eq!(flag, INTERVAL_HAS_OVERLAP);
        assert_eq!(buffer.full_windows(), &[0, 3]);
        assert_eq!(buffer.partial_length, 0);
    }

    #[test]
    fn available_interval_reports_when_no_overlap_touches_the_window() {
        let overlaps = vec![OverlapRegion {
            x_position_start: 200,
            x_position_end: 299,
            is_match: 1,
            ..OverlapRegion::default()
        }];
        let mut buffer = IntervalBuffer::with_capacity(4);
        assert_eq!(
            get_available_interval(0, 99, &overlaps, &mut buffer),
            INTERVAL_NO_OVERLAP
        );
        assert_eq!(buffer.length, 0);
    }

    /// Golden values printed by the real `if_is_homopolymer_strict`.
    #[test]
    fn homopolymer_test_matches_the_reference() {
        let cases = [
            (3i64, "AAAACGT", true),
            (2, "ACCCGT", true),
            (0, "CCCCGT", true),
            (2, "AACCG", false),
            (3, "ACGTACGT", false),
            (0, "ACGT", false),
            (1, "AAACGT", true),
            (2, "AAACGT", true),
            (3, "AAACGT", true),
            (0, "AAAA", true),
            (3, "AAAA", true),
        ];
        for (site, read, expected) in cases {
            assert_eq!(
                is_strict_homopolymer(site, read.as_bytes()),
                expected,
                "site {site} of {read}"
            );
        }
        // Out-of-range positions are rejected.
        assert!(!is_strict_homopolymer(-1, b"AAAA"));
        assert!(!is_strict_homopolymer(4, b"AAAA"));
    }

    /// Golden values printed by the real `get_available_interval`.
    #[test]
    fn available_interval_matches_the_reference() {
        let make = |start: u32, is_match: u8| OverlapRegion {
            x_position_start: start,
            x_position_end: start + 99,
            is_match,
            ..OverlapRegion::default()
        };

        let overlaps = vec![make(0, 1), make(0, 0), make(0, 3), make(50, 1)];
        let mut buffer = IntervalBuffer::with_capacity(16);
        assert_eq!(
            get_available_interval(0, 99, &overlaps, &mut buffer),
            INTERVAL_HAS_OVERLAP
        );
        assert_eq!(buffer.full_windows(), &[0, 3]);
        assert_eq!(buffer.start_index, 0);

        // Overlaps touch the window but none is a real match: the flag still
        // reports `1`, only the identifier list stays empty.
        let overlaps = vec![make(0, 0), make(100, 0)];
        let mut buffer = IntervalBuffer::with_capacity(16);
        assert_eq!(
            get_available_interval(0, 99, &overlaps, &mut buffer),
            INTERVAL_HAS_OVERLAP
        );
        assert_eq!(buffer.length, 0);
        assert_eq!(buffer.start_index, 0);
    }
}

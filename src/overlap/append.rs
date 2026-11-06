//! Converting a chained overlap into the stored overlap record, the
//! counterpart of `append_inexact_overlap_region_alloc` in `Hash_Table.cpp`.

use crate::overlap::region::{OverlapRegion, chain_length};

/// `append_inexact_overlap_region_alloc`.
///
/// The chained interval is first widened to the read ends (the same
/// computation `get_chainLen` performs), mirrored when the chain runs on the
/// reverse strand, and finally appended unless the previous entry describes
/// the same target read with a better score.
///
/// `add_begin_and_end` is the reference `add_beg_end`: when set, the overlap
/// is padded with the first and last base positions so that later stages see
/// the full span.
pub fn append_inexact_overlap(
    overlap_list: &mut Vec<OverlapRegion>,
    chained: &OverlapRegion,
    x_read_length: i64,
    y_read_length: i64,
    add_begin_and_end: bool,
) -> bool {
    // The reference drops the new overlap when the previous one describes the
    // same target read with an equal or better score.
    if let Some(previous) = overlap_list.last() {
        if previous.y_identifier == chained.y_identifier
            && (previous.shared_seed > chained.shared_seed
                || (previous.shared_seed == chained.shared_seed
                    && previous.overlap_length <= chained.overlap_length))
        {
            return false;
        }
        if previous.y_identifier == chained.y_identifier {
            overlap_list.pop();
        }
    }

    let mut region = chained.clone();

    // Align the two starts, then extend both ends to the read boundaries.
    if region.x_position_start <= region.y_position_start {
        region.y_position_start -= region.x_position_start;
        region.x_position_start = 0;
    } else {
        region.x_position_start -= region.y_position_start;
        region.y_position_start = 0;
    }

    let x_right_length = x_read_length - region.x_position_end as i64 - 1;
    let y_right_length = y_read_length - region.y_position_end as i64 - 1;
    if x_right_length <= y_right_length {
        region.y_position_end = (region.y_position_end as i64 + x_right_length) as u32;
        region.x_position_end = (x_read_length - 1) as u32;
    } else {
        region.x_position_end = (region.x_position_end as i64 + y_right_length) as u32;
        region.y_position_end = (y_read_length - 1) as u32;
    }

    let mut fake_cigar = crate::overlap::region::FakeCigar::default();
    if region.x_strand == 1 {
        // Mirror the coordinates onto the reverse complement. Note that the
        // start takes the mirrored end and vice versa, exactly as the
        // reference does.
        let mirrored_start = x_read_length - region.x_position_end as i64 - 1;
        let mirrored_end = x_read_length - region.x_position_start as i64 - 1;
        let mirrored_y_start = y_read_length - region.y_position_end as i64 - 1;
        let mirrored_y_end = y_read_length - region.y_position_start as i64 - 1;
        region.x_position_start = mirrored_start as u32;
        region.x_position_end = mirrored_end as u32;
        region.x_strand = 0;
        region.y_position_start = mirrored_y_start as u32;
        region.y_position_end = mirrored_y_end as u32;
        region.y_strand = 1;

        if add_begin_and_end {
            fake_cigar.push(region.x_position_start, 0);
        }

        // The chained gap sites were produced on the reverse strand, so they
        // are mirrored and duplicates removed.
        let mut previous_shift = 0x0fff_ffff_ffff_ffffi64;
        for index in 0..chained.fake_cigar.len() {
            let shift = chained.fake_cigar.gap_shift(index);
            if shift != previous_shift {
                previous_shift = shift;
                let mirrored_site =
                    x_read_length - chained.fake_cigar.gap_position(index) - 1;
                fake_cigar.push(mirrored_site as u32, previous_shift as i32);
            }
        }

        if add_begin_and_end
            && fake_cigar.gap_position(fake_cigar.len() - 1) != region.x_position_end as i64
        {
            let last_shift = fake_cigar.gap_shift(fake_cigar.len() - 1);
            fake_cigar.push(region.x_position_end, last_shift as i32);
        }
    } else {
        if add_begin_and_end {
            fake_cigar.push(region.x_position_start, 0);
        }

        let distance_self = region.x_position_end as i64 - region.x_position_s() as i64;
        let distance_other = region.y_position_end as i64 - region.y_position_s() as i64;
        let initial_shift = distance_other - distance_self;

        let mut previous_shift = 0x0fff_ffff_ffff_ffffi64;
        for index in (0..chained.fake_cigar.len()).rev() {
            let shift = chained.fake_cigar.gap_shift(index);
            if shift != previous_shift {
                previous_shift = shift;
                fake_cigar.push(
                    chained.fake_cigar.gap_position(index) as u32,
                    (initial_shift - previous_shift) as i32,
                );
            }
        }

        if add_begin_and_end
            && fake_cigar.gap_position(fake_cigar.len() - 1) != region.x_position_end as i64
        {
            let last_shift = fake_cigar.gap_shift(fake_cigar.len() - 1);
            fake_cigar.push(region.x_position_end, last_shift as i32);
        }
    }
    region.fake_cigar = fake_cigar;

    region.alignment_length = 0;
    region.is_match = 0;
    region.non_homopolymer_errors = 0;
    region.strong = 0;

    overlap_list.push(region);
    true
}

#[allow(dead_code)]
fn unused_chain_length_helper() -> i64 {
    chain_length(0, 0, 1, 0, 0, 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chained(strand: u32, x_start: u32, x_end: u32, y_start: u32, y_end: u32, seed: i32) -> OverlapRegion {
        OverlapRegion {
            x_identifier: 0,
            x_position_start: x_start,
            x_position_end: x_end,
            x_strand: strand,
            y_identifier: 3,
            y_position_start: y_start,
            y_position_end: y_end,
            y_strand: 0,
            shared_seed: seed,
            overlap_length: 100,
            fake_cigar: {
                let mut cigar = crate::overlap::region::FakeCigar::default();
                cigar.push(x_start, 0);
                cigar
            },
            ..OverlapRegion::default()
        }
    }

    #[test]
    fn forward_overlap_extends_to_the_read_ends() {
        let mut list: Vec<OverlapRegion> = Vec::new();
        let appended = append_inexact_overlap(
            &mut list,
            &chained(0, 100, 500, 100, 500, 10),
            1000,
            1000,
            false,
        );

        assert!(appended);
        assert_eq!(list.len(), 1);
        let region = &list[0];
        assert_eq!(region.x_position_start, 0);
        assert_eq!(region.x_position_end, 999);
        assert_eq!(region.y_position_start, 0);
        assert_eq!(region.y_position_end, 999);
        assert_eq!(region.shared_seed, 10);
        assert_eq!(region.x_identifier, 0);
        assert_eq!(region.y_identifier, 3);
    }

    #[test]
    fn reverse_overlap_mirrors_the_query_coordinates() {
        let mut list: Vec<OverlapRegion> = Vec::new();
        // A reverse-strand chain: 300..700 on an 800 base query against
        // 100..700 on a 1000 base target.
        let mut source = chained(1, 300, 700, 100, 700, 5);
        source.fake_cigar.push(700, 2);

        assert!(append_inexact_overlap(&mut list, &source, 800, 1000, false));
        let region = &list[0];
        assert_eq!(region.x_strand, 0);
        assert_eq!(region.y_strand, 1);
        // After widening: query 200..799, target 0..799. Mirrored:
        // query start 800-799-1 = 0, end 800-200-1 = 599;
        // target start 1000-799-1 = 200, end 1000-0-1 = 999.
        assert_eq!(region.x_position_start, 0);
        assert_eq!(region.x_position_end, 599);
        assert_eq!(region.y_position_start, 200);
        assert_eq!(region.y_position_end, 999);
    }

    #[test]
    fn repeated_targets_keep_the_better_chain() {
        let mut list: Vec<OverlapRegion> = Vec::new();
        assert!(append_inexact_overlap(
            &mut list,
            &chained(0, 100, 500, 100, 500, 5),
            1000,
            1000,
            false
        ));
        // Same target, worse score: dropped.
        assert!(!append_inexact_overlap(
            &mut list,
            &chained(0, 100, 500, 100, 500, 4),
            1000,
            1000,
            false
        ));
        // Same target, better score: replaces.
        assert!(append_inexact_overlap(
            &mut list,
            &chained(0, 100, 500, 100, 500, 6),
            1000,
            1000,
            false
        ));
        assert_eq!(list.len(), 1);
        assert_eq!(list[0].shared_seed, 6);
    }

    /// Golden values printed by the real `append_inexact_overlap_region_alloc`
    /// object file for the same four inputs.
    #[test]
    fn matches_the_reference_implementation() {
        struct Case {
            name: &'static str,
            strand: u32,
            add_begin_and_end: bool,
            x_length: i64,
            y_length: i64,
            start: (u32, u32),
            end: (u32, u32),
            expected_x: (u32, u32),
            expected_y: (u32, u32),
            expected_strands: (u32, u32),
            expected_cigar: Vec<(i64, i64)>,
        }

        let cases = vec![
            Case {
                name: "forward_be0",
                strand: 0,
                add_begin_and_end: false,
                x_length: 1000,
                y_length: 1000,
                start: (100, 100),
                end: (500, 500),
                expected_x: (0, 999),
                expected_y: (0, 999),
                expected_strands: (0, 0),
                expected_cigar: vec![(500, -2), (100, 0)],
            },
            Case {
                name: "forward_be1",
                strand: 0,
                add_begin_and_end: true,
                x_length: 1000,
                y_length: 1000,
                start: (100, 100),
                end: (500, 500),
                expected_x: (0, 999),
                expected_y: (0, 999),
                expected_strands: (0, 0),
                expected_cigar: vec![(0, 0), (500, -2), (100, 0), (999, 0)],
            },
            Case {
                name: "reverse_be0",
                strand: 1,
                add_begin_and_end: false,
                x_length: 800,
                y_length: 1000,
                start: (300, 100),
                end: (700, 700),
                expected_x: (0, 599),
                expected_y: (200, 999),
                expected_strands: (0, 1),
                expected_cigar: vec![(499, 0), (99, 2)],
            },
            Case {
                name: "reverse_be1",
                strand: 1,
                add_begin_and_end: true,
                x_length: 800,
                y_length: 1000,
                start: (300, 100),
                end: (700, 700),
                expected_x: (0, 599),
                expected_y: (200, 999),
                expected_strands: (0, 1),
                expected_cigar: vec![(0, 0), (499, 0), (99, 2), (599, 2)],
            },
        ];

        for case in cases {
            let mut list: Vec<OverlapRegion> = Vec::new();
            let source = OverlapRegion {
                x_identifier: 0,
                y_identifier: 3,
                x_position_start: case.start.0,
                x_position_end: case.end.0,
                x_strand: case.strand,
                y_position_start: case.start.1,
                y_position_end: case.end.1,
                y_strand: 0,
                shared_seed: 5,
                overlap_length: 100,
                ..OverlapRegion::default()
            };
            let mut source = source;
            source.fake_cigar.push(source.x_position_start, 0);
            source.fake_cigar.push(source.x_position_end, 2);

            let appended = append_inexact_overlap(
                &mut list,
                &source,
                case.x_length,
                case.y_length,
                case.add_begin_and_end,
            );

            assert!(appended, "{} must be appended", case.name);
            let region = &list[0];
            assert_eq!(
                (region.x_position_start, region.x_position_end),
                case.expected_x,
                "x of {}",
                case.name
            );
            assert_eq!(
                (region.y_position_start, region.y_position_end),
                case.expected_y,
                "y of {}",
                case.name
            );
            assert_eq!(
                (region.x_strand, region.y_strand),
                case.expected_strands,
                "strands of {}",
                case.name
            );
            let entries: Vec<(i64, i64)> = (0..region.fake_cigar.len())
                .map(|index| {
                    (
                        region.fake_cigar.gap_position(index),
                        region.fake_cigar.gap_shift(index),
                    )
                })
                .collect();
            assert_eq!(entries, case.expected_cigar, "cigar of {}", case.name);
        }
    }

    #[test]
    fn begin_and_end_padding_extends_the_fake_cigar() {
        let mut list: Vec<OverlapRegion> = Vec::new();
        assert!(append_inexact_overlap(
            &mut list,
            &chained(0, 100, 500, 100, 500, 5),
            1000,
            1000,
            true
        ));
        let cigar = &list[0].fake_cigar;
        assert!(cigar.len() >= 2);
        assert_eq!(cigar.gap_position(0), 0);
        let last = cigar.len() - 1;
        assert_eq!(cigar.gap_position(last), 999);
    }
}

//! Stored overlap list and its PAF rendering, the counterpart of
//! `ma_hit_t`, `push_overlaps` and `Output_PAF`.
//!
//! The C structure packs a few flags into bit fields; here they are plain
//! fields because the value ranges are identical and the code reads better.

use std::io::{self, Write};

use crate::overlap::region::OverlapRegion;
use crate::reads::read_store::ReadStore;

/// One stored overlap, `ma_hit_t`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PafHit {
    /// Read the overlap starts on.
    pub query_identifier: u32,
    pub query_start: u32,
    pub query_end: u32,
    /// Read the overlap points at.
    pub target_identifier: u32,
    pub target_start: u32,
    pub target_end: u32,
    /// `ml`: whether the overlap survived the refinement, one bit in C.
    pub matching_length: u32,
    /// `bl`: length of the target read's block.
    pub block_length: u32,
    /// `rev`: whether the target is on the opposite strand.
    pub reverse: bool,
    /// `no_l_indel`: no gap larger than the configured threshold.
    pub without_large_indel: bool,
    /// `el`: the hit stores the maximal exact sub-region of the overlap
    /// rather than the whole overlap. `push_ne_ovlp` sets it after
    /// `extract_max_exact` finds one, and rewrites the four coordinates above
    /// to that sub-region, with the ends exclusive.
    pub exact_region: bool,
    /// The reference stores `shared_seed` itself in the one-byte `el` field;
    /// the pair-mode caches can therefore carry values other than one.
    pub exact_region_byte: u8,
    /// `del`: the assembly stage's deletion flag. PAF output ignores it.
    pub deleted: bool,
}

/// `push_overlaps`: copies the regions flagged `is_match == flag` into `paf`.
///
/// `reverse_coordinates` is the reference `if_reverse`: when set, every
/// coordinate is mirrored so that the list describes the reverse complement.
pub fn push_overlaps(
    paf: &mut Vec<PafHit>,
    overlap_list: &[OverlapRegion],
    flag: u8,
    reads: &ReadStore,
    reverse_coordinates: bool,
) {
    paf.clear();
    for region in overlap_list {
        if region.is_match != flag {
            continue;
        }

        let x_length = reads.read_length(region.x_identifier as usize);
        let y_length = reads.read_length(region.y_identifier as usize);

        let (query_start, query_end, target_start, target_end) = if reverse_coordinates {
            (
                (x_length - region.x_position_start as u64 - 1) as u32,
                (x_length - region.x_position_end as u64 - 1) as u32,
                (y_length - region.y_position_start as u64 - 1) as u32,
                (y_length - region.y_position_end as u64 - 1) as u32,
            )
        } else {
            (
                region.x_position_start,
                region.x_position_end,
                region.y_position_start,
                region.y_position_end,
            )
        };

        paf.push(PafHit {
            query_identifier: region.x_identifier,
            query_start,
            query_end,
            target_identifier: region.y_identifier,
            target_start,
            target_end,
            // The reference stores `strong` in the one-bit `ml` field.
            matching_length: region.strong as u32,
            block_length: y_length as u32,
            // Every stored overlap has the query on the forward strand, so the
            // query strand flag is the target strand flag.
            reverse: region.y_strand == 1,
            without_large_indel: region.without_large_indel != 0,
            // `push_overlaps` leaves `el` untouched in the reference; the
            // whole-overlap coordinate convention it stores is never consumed
            // as an exact region, so the flag stays clear.
            exact_region: false,
            exact_region_byte: 0,
            deleted: false,
        });
    }
}

/// `Output_PAF`: writes every stored overlap in PAF format.
///
/// The field order, the tabs and the trailing `255` map quality are part of
/// the byte-for-byte contract with the reference. Each element of `paf` is one
/// read's overlap list, in any type that exposes it as a `PafHit` slice, such
/// as `Vec<PafHit>` or the correction driver's `PafList`.
pub fn write_paf<W: Write, T: AsRef<[PafHit]>>(
    writer: &mut W,
    reads: &ReadStore,
    paf: &[T],
) -> io::Result<()> {
    for hits in paf {
        for hit in hits.as_ref() {
            let query_name = reads.read_name(hit.query_identifier as usize);
            let target_name = reads.read_name(hit.target_identifier as usize);
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255",
                query_name,
                reads.read_length(hit.query_identifier as usize),
                hit.query_start,
                hit.query_end,
                if hit.reverse { '-' } else { '+' },
                target_name,
                reads.read_length(hit.target_identifier as usize),
                hit.target_start,
                hit.target_end,
                hit.matching_length,
                hit.block_length
            )?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn store() -> ReadStore {
        let mut reads = ReadStore::default();
        reads.push_read("read_a", b"ACGTACGTAC");
        reads.push_read("read_b", b"TTTTGGGGCC");
        reads
    }

    fn region(start: u32, end: u32) -> OverlapRegion {
        OverlapRegion {
            x_identifier: 0,
            x_position_start: start,
            x_position_end: end,
            y_identifier: 1,
            y_position_start: start + 1,
            y_position_end: end + 1,
            y_strand: 0,
            strong: 1,
            is_match: 1,
            ..OverlapRegion::default()
        }
    }

    #[test]
    fn only_matching_regions_are_pushed() {
        let reads = store();
        let mut paf: Vec<PafHit> = Vec::new();
        let mut unmatched = region(1, 5);
        unmatched.is_match = 0;

        push_overlaps(&mut paf, &[region(1, 5), unmatched], 1, &reads, false);
        assert_eq!(paf.len(), 1);
        assert_eq!(paf[0].query_start, 1);
        assert_eq!(paf[0].query_end, 5);
        assert_eq!(paf[0].target_start, 2);
        assert_eq!(paf[0].block_length, 10);
        assert_eq!(paf[0].matching_length, 1);
    }

    #[test]
    fn reverse_coordinates_are_mirrored() {
        let reads = store();
        let mut paf: Vec<PafHit> = Vec::new();
        push_overlaps(&mut paf, &[region(1, 5)], 1, &reads, true);
        // Query length 10: 10-1-1 = 8 and 10-5-1 = 4.
        assert_eq!(paf[0].query_start, 8);
        assert_eq!(paf[0].query_end, 4);
    }

    #[test]
    fn paf_lines_follow_the_reference_layout() {
        let reads = store();
        let mut paf: Vec<PafHit> = Vec::new();
        push_overlaps(&mut paf, &[region(1, 5)], 1, &reads, false);

        let mut output: Vec<u8> = Vec::new();
        write_paf(&mut output, &reads, &[paf]).unwrap();
        let text = String::from_utf8(output).unwrap();
        assert_eq!(text, "read_a\t10\t1\t5\t+\tread_b\t10\t2\t6\t1\t10\t255\n");
    }
}

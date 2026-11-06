//! Overlap bookkeeping types, the counterpart of the structures in
//! `Hash_Table.h` that the chaining code fills in.

use crate::container::ksort::radix_sort_by_key_32;

/// `overlap_region_sort_y_id`: ascending order of `y_id`.
///
/// The reference instantiates `KRADIX_SORT_INIT` over the four-byte `y_id`,
/// so the port sorts an index permutation with the same 32-bit radix sort,
/// which yields the exact same permutation as sorting the records themselves.
pub fn sort_by_y_identifier(overlaps: &mut [OverlapRegion]) {
    let mut order: Vec<u32> = (0..overlaps.len() as u32).collect();
    radix_sort_by_key_32(&mut order, &|&index| overlaps[index as usize].y_identifier);
    let sorted: Vec<OverlapRegion> = order
        .iter()
        .map(|&index| overlaps[index as usize].clone())
        .collect();
    overlaps.clone_from_slice(&sorted);
}

/// One shared seed between two reads, `k_mer_hit`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct KmerHit {
    /// Read the seed points at.
    pub read_identifier: u32,
    /// Orientation of the hit.
    pub reverse_strand: bool,
    /// Position on the target read, the chaining sort key.
    pub offset: u32,
    /// Position on the query read.
    pub self_offset: u32,
    /// How many minimizers support the seed.
    pub count: u32,
}

impl KmerHit {
    /// Strand as the `0`/`1` value the reference stores in a bit field.
    pub fn strand_value(&self) -> u32 {
        u32::from(self.reverse_strand)
    }

    /// Low byte of `cnt`: the length of the seed run in bases.
    pub fn span(&self) -> u32 {
        self.count & 0xff
    }

    /// High bits of `cnt`: how many seeds support the run.
    pub fn coverage(&self) -> u32 {
        self.count >> 8
    }
}

/// A gap in the fake CIGAR, stored the way `Fake_Cigar` packs it: the site in
/// the high 32 bits and `(|shift| << 1 | is_negative)` in the low 32 bits.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FakeCigar {
    pub entries: Vec<u64>,
}

impl FakeCigar {
    /// `clear_fake_cigar`
    pub fn clear(&mut self) {
        self.entries.clear();
    }

    /// `add_fake_cigar`
    pub fn push(&mut self, gap_site: u32, gap_shift: i32) {
        let mut packed = (gap_site as u64) << 32;
        let shift_field: u32 = if gap_shift < 0 {
            ((gap_shift.unsigned_abs()) << 1) | 1
        } else {
            (gap_shift as u32) << 1
        };
        packed |= shift_field as u64;
        self.entries.push(packed);
    }

    /// `get_fake_gap_pos`
    pub fn gap_position(&self, index: usize) -> i64 {
        (self.entries[index] >> 32) as i64
    }

    /// `get_fake_gap_shift`
    pub fn gap_shift(&self, index: usize) -> i64 {
        let raw = self.entries[index] as u32;
        let magnitude = (raw >> 1) as i64;
        if raw & 1 == 1 {
            -magnitude
        } else {
            magnitude
        }
    }

    /// `y_start_offset`: the gap shift in effect at query position `x_start`.
    ///
    /// The gap structure is a step function over the query: between two gap
    /// sites the offset stays constant, so the offset of `x_start` is the shift
    /// of the last site at or before it. The reference treats the exact final
    /// site as a special case and aborts when `x_start` falls outside the
    /// recorded sites, which can only happen for a corrupt chain.
    pub fn y_start_offset(&self, x_start: i64) -> Option<i64> {
        let last = self.entries.len().checked_sub(1)?;
        if x_start == self.gap_position(last) {
            return Some(self.gap_shift(last));
        }

        let mut index = 0usize;
        while index < self.entries.len() {
            if x_start < self.gap_position(index) {
                break;
            }
            index += 1;
        }

        // `index == 0` means the position precedes every site, and
        // `index == len` means it follows the final one.
        if index == 0 || index == self.entries.len() {
            return None;
        }
        Some(self.gap_shift(index - 1))
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// The parts of `overlap_region` that chaining produces.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct OverlapRegion {
    /// Query read, `x_id`.
    pub x_identifier: u32,
    /// Query start and end.
    pub x_position_start: u32,
    pub x_position_end: u32,
    /// Strand of the query, `x_pos_strand`.
    pub x_strand: u32,
    /// Target read, `y_id`.
    pub y_identifier: u32,
    /// Target start and end.
    pub y_position_start: u32,
    pub y_position_end: u32,
    /// Strand of the target, `y_pos_strand`.
    pub y_strand: u32,
    /// Number of shared seeds supporting the chain.
    pub shared_seed: i32,
    /// Length of the overlap the chain spans.
    pub overlap_length: u32,
    /// Alignment length accumulated by the refinement, zero after chaining.
    pub alignment_length: u32,
    /// Whether every window passed the refinement.
    pub is_match: u8,
    /// No gap larger than the configured threshold.
    pub without_large_indel: u8,
    /// Refinement outcome, `strong`.
    pub strong: i8,
    /// Substitution count, `non_homopolymer_errors`.
    pub non_homopolymer_errors: u32,
    /// Gap structure of the chain.
    pub fake_cigar: FakeCigar,
    /// Windows verified against this overlap, `overlap_region::w_list`.
    pub window_cigar: crate::correction::window::WindowCigar,
    /// Cigars of the window boundaries, `overlap_region::boundary_cigars`.
    pub boundary_cigars: crate::correction::window::WindowCigar,
}


impl OverlapRegion {
    /// Read-only view used by the append code before mirroring.
    pub fn x_position_s(&self) -> u32 {
        self.x_position_start
    }

    /// Read-only view used by the append code before mirroring.
    pub fn x_position_e(&self) -> u32 {
        self.x_position_end
    }

    /// Read-only view used by the append code before mirroring.
    pub fn y_position_s(&self) -> u32 {
        self.y_position_start
    }

    /// Read-only view used by the append code before mirroring.
    pub fn y_position_e(&self) -> u32 {
        self.y_position_end
    }

    /// `update_overlap_region`: replaces the window list with `source` and
    /// re-derives the overlap's start and end coordinates from it.
    ///
    /// The coordinates are first slid so that one of the two reads starts at
    /// zero, then extended to the far read end along the same diagonal. The
    /// reference does this in `uint32_t` arithmetic, so a negative remaining
    /// length wraps; the wrapping operations mirror that.
    pub fn update_from_windows(
        &mut self,
        source: &crate::correction::window::WindowCigar,
        query_length: i64,
        target_length: i64,
    ) {
        self.window_cigar = source.clone();

        if let (Some(first), Some(last)) =
            (source.windows.first(), source.windows.last())
        {
            self.x_position_start = first.x_start as u32;
            self.x_position_end = last.x_end as u32;
            self.y_position_start = first.y_start as u32;
            self.y_position_end = last.y_end as u32;
        }

        if self.x_position_start <= self.y_position_start {
            self.y_position_start = self.y_position_start.wrapping_sub(self.x_position_start);
            self.x_position_start = 0;
        } else {
            self.x_position_start = self.x_position_start.wrapping_sub(self.y_position_start);
            self.y_position_start = 0;
        }

        let query_remaining = query_length - i64::from(self.x_position_end) - 1;
        let target_remaining = target_length - i64::from(self.y_position_end) - 1;
        if query_remaining <= target_remaining {
            self.x_position_end = (query_length - 1) as u32;
            self.y_position_end = (i64::from(self.y_position_end) + query_remaining) as u32;
        } else {
            self.y_position_end = (target_length - 1) as u32;
            self.x_position_end = (i64::from(self.x_position_end) + target_remaining) as u32;
        }
    }
}

/// DP arrays for chaining, `Chain_Data`.
#[derive(Debug, Clone, Default)]
pub struct ChainData {
    pub score: Vec<i32>,
    pub predecessor: Vec<i64>,
    pub indels: Vec<i32>,
    pub self_length: Vec<i32>,
    pub occurrence: Vec<i32>,
    /// Scratch marks used to skip redundant predecessors.
    pub scratch: Vec<i64>,
}

impl ChainData {
    /// `resize_Chain_Data`: makes room for `size` anchors.
    pub fn resize(&mut self, size: usize) {
        let needed = size + 1;
        if needed > self.score.len() {
            self.score.resize(needed, 0);
            self.predecessor.resize(needed, -1);
            self.indels.resize(needed, 0);
            self.self_length.resize(needed, 0);
            self.occurrence.resize(needed, 0);
            self.scratch.resize(needed, -1);
        }
    }
}

/// `Candidates_list`: the anchor buffer of one read plus the chaining
/// workspace reused across every overlap of that read.
///
/// `length` is the reference's fill pointer, not the buffer length: the
/// window-rechaining code temporarily rewinds it to the start of the region it
/// is about to rebuild and lets the algorithm overwrite the old anchors in
/// place.
#[derive(Debug, Clone, Default)]
pub struct CandidatesList {
    pub list: Vec<KmerHit>,
    pub length: i64,
    pub size: i64,
    pub chain_data: ChainData,
}

impl CandidatesList {
    /// `clear_Candidates_list`.
    pub fn clear(&mut self) {
        self.length = 0;
    }

    /// Appends `hit` at `length`, growing the buffer like `kv_pushp_cl`.
    pub fn push(&mut self, hit: KmerHit) {
        let index = self.length as usize;
        if self.list.len() <= index {
            self.list.resize(index + 1, KmerHit::default());
        }
        self.list[index] = hit;
        self.length += 1;
        if self.size < self.length {
            self.size = self.length;
        }
    }

    /// Makes sure `count` entries starting at `length` exist, like the
    /// `REALLOC` in `gen_single_khit`.
    pub fn reserve_at_length(&mut self, count: i64) {
        let needed = self.length + count;
        if needed > self.list.len() as i64 {
            self.list.resize(needed as usize, KmerHit::default());
        }
        if self.size < needed {
            self.size = needed;
        }
    }
}

/// `get_chainLen`: length of the overlap described by the two intervals.
///
/// The reference also tracks the resulting `y` interval, but the return value
/// only depends on the `x` side, so only that side is materialised here.
pub fn chain_length(
    x_begin: i64,
    x_end: i64,
    x_length: i64,
    y_begin: i64,
    y_end: i64,
    y_length: i64,
) -> i64 {
    // Align the two starts, then extend to whichever read ends first.
    let x_begin = if x_begin <= y_begin {
        0
    } else {
        x_begin - y_begin
    };

    let x_right_length = x_length - x_end - 1;
    let y_right_length = y_length - y_end - 1;

    let x_end = if x_right_length <= y_right_length {
        x_length - 1
    } else {
        x_end + y_right_length
    };

    x_end - x_begin + 1
}

/// `normal_w`: the anchor score normalised by the seed support.
pub fn normal_weight(value: i64, count: i64) -> i64 {
    if value >= count {
        value / count
    } else {
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normal_weight_never_drops_below_one() {
        assert_eq!(normal_weight(51, 1), 51);
        assert_eq!(normal_weight(51, 3), 17);
        assert_eq!(normal_weight(2, 5), 1);
    }

    #[test]
    fn chain_length_covers_both_reads() {
        // The interval is extended to the read ends on both sides.
        assert_eq!(chain_length(100, 500, 1000, 100, 500, 1000), 1000);
        assert_eq!(chain_length(0, 0, 100, 0, 0, 100), 100);
    }

    /// Golden values printed by the real `y_start_offset` for the gap sites
    /// `(0, 0), (100, 5), (200, -3), (300, 8)`. The two inputs outside that
    /// range make the reference abort, so they map to `None` here.
    #[test]
    fn y_start_offset_matches_the_reference_step_function() {
        let mut cigar = FakeCigar::default();
        for (site, shift) in [(0u32, 0i32), (100, 5), (200, -3), (300, 8)] {
            cigar.push(site, shift);
        }

        assert_eq!(cigar.y_start_offset(50), Some(0));
        assert_eq!(cigar.y_start_offset(100), Some(5), "a gap site takes its own shift");
        assert_eq!(cigar.y_start_offset(150), Some(5), "between sites the shift holds");
        assert_eq!(cigar.y_start_offset(300), Some(8), "the final site is a special case");
        assert_eq!(cigar.y_start_offset(350), None, "past the final site");
        assert_eq!(cigar.y_start_offset(-1), None, "before the first site");
        assert_eq!(FakeCigar::default().y_start_offset(0), None, "empty cigar");
    }

    #[test]
    fn fake_cigar_packs_site_and_signed_shift() {
        let mut cigar = FakeCigar::default();
        cigar.push(100, 5);
        assert_eq!(cigar.entries[0], (100u64 << 32) | 10);

        cigar.push(200, -3);
        assert_eq!(cigar.entries[1], (200u64 << 32) | 7);
    }
}

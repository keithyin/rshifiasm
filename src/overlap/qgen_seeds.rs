//! Candidate seeds of the new correction chain, `minimizers_qgen0` in
//! `anchor.cpp`.
//!
//! Unlike the older [`crate::overlap::candidates::collect_seeds`], this pass
//! sorts anchors by the query position, mirrors target positions on the
//! reverse strand, keeps the minimizer span in the low eight bits of the
//! weight, and derives occurrence bounds from the caller supplied
//! `high_occurrence` / `low_occurrence` instead of the homozygous coverage.

use crate::container::ksort::radix_sort_by_key;
use crate::counting::minimizer::Minimizer;
use crate::counting::position_index::PositionIndex;
use crate::overlap::region::KmerHit;

/// One anchor before it is converted into a `k_mer_hit`, `anchor1_t`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct QgenAnchor {
    /// `read_identifier << 33 | reverse << 32 | self_offset`.
    sort_key: u64,
    /// Position on the target, mirrored for the reverse strand.
    other_offset: u32,
    /// `occurrence << 8 | span`.
    packed_count: u32,
}

impl QgenAnchor {
    fn read_identifier(&self) -> u32 {
        (self.sort_key >> 33) as u32
    }

    fn reverse_strand(&self) -> bool {
        (self.sort_key >> 32) & 1 == 1
    }

    fn self_offset(&self) -> u32 {
        self.sort_key as u32
    }

    fn occurrence(&self) -> u32 {
        self.packed_count >> 8
    }

    fn span(&self) -> u8 {
        (self.packed_count & 0xff) as u8
    }
}

/// `high_occ`/`low_occ` weighting of `minimizers_qgen0`.
///
/// `max_occurrence` and `min_occurrence` are clamped to two by the reference,
/// and the packed span is preserved in the low eight bits of the result.
pub fn qgen_seed_weight(occurrence: u32, span: u8, max_occurrence: u32, min_occurrence: u32) -> u32 {
    let maximum = max_occurrence.max(2);
    let minimum = min_occurrence.max(2);
    let doubled_maximum = maximum << 1;

    let weight = if occurrence < maximum && occurrence > minimum {
        1u32
    } else if occurrence <= minimum {
        2u32
    } else {
        let scaled = 1 + (occurrence + doubled_maximum - 1) / doubled_maximum;
        (scaled as f64).powf(1.1) as u32
    };

    let capped = weight.min(0x00ff_ffff);
    (capped << 8) | span as u32
}

/// `minimizers_qgen0`: turns the minimizers of one read into the sorted seed
/// list consumed by `lchain_qgen_mcopy_fast`.
///
/// `read_lengths` gives the length of every target read, which is needed to
/// undo the reverse-strand mirroring of target positions.
pub fn collect_qgen_seeds(
    minimizers: &[Minimizer],
    index: &PositionIndex,
    read_lengths: &[u64],
    max_occurrence: u32,
    min_occurrence: u32,
) -> Vec<KmerHit> {
    let mut anchors: Vec<QgenAnchor> = Vec::new();

    for minimizer in minimizers {
        let positions = index.get(minimizer.hash);
        let occurrence = (positions.len() as u32).min(0x00ff_ffff);
        let span = minimizer.span.min(255) as u8;
        let packed_count = (occurrence << 8) | span as u32;

        for position in positions {
            let reverse = minimizer.reverse != position.reverse;
            let other_offset = if reverse {
                // `((uint32_t)-1)-1-(pos+1-span)`: the mirrored position of
                // the last base of the k-mer on the opposite strand.
                u32::MAX
                    .wrapping_sub(1)
                    .wrapping_sub(position.position.wrapping_add(1).wrapping_sub(position.span))
            } else {
                position.position
            };

            anchors.push(QgenAnchor {
                sort_key: ((position.read_identifier as u64) << 33)
                    | ((reverse as u64) << 32)
                    | minimizer.position as u64,
                other_offset,
                packed_count,
            });
        }
    }

    radix_sort_by_key(&mut anchors, &|anchor: &QgenAnchor| anchor.sort_key);

    // Inside each equal-key run sort by the mirrored target offset.
    let mut run_start = 0usize;
    for index in 1..=anchors.len() {
        if index == anchors.len() || anchors[index].sort_key != anchors[run_start].sort_key {
            if index - run_start > 1 {
                radix_sort_by_key(&mut anchors[run_start..index], &|anchor: &QgenAnchor| {
                    anchor.other_offset as u64
                });
            }
            run_start = index;
        }
    }

    anchors
        .iter()
        .map(|anchor| {
            let target_length = read_lengths
                .get(anchor.read_identifier() as usize)
                .copied()
                .unwrap_or(0) as u32;
            let offset = if anchor.reverse_strand() {
                let mirrored = u32::MAX - anchor.other_offset;
                target_length.wrapping_sub(mirrored)
            } else {
                anchor.other_offset
            };

            KmerHit {
                read_identifier: anchor.read_identifier(),
                reverse_strand: anchor.reverse_strand(),
                offset,
                self_offset: anchor.self_offset(),
                count: qgen_seed_weight(
                    anchor.occurrence(),
                    anchor.span(),
                    max_occurrence,
                    min_occurrence,
                ),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::counting::position_index::PositionRecord;

    #[test]
    fn weight_preserves_the_span_and_caps_the_weight() {
        assert_eq!(qgen_seed_weight(30, 51, 60, 12), (1 << 8) | 51);
        assert_eq!(qgen_seed_weight(5, 51, 60, 12), (2 << 8) | 51);
        assert_eq!(qgen_seed_weight(200, 51, 60, 12), (3 << 8) | 51);
        assert_eq!(qgen_seed_weight(1, 0, 2, 2), (2 << 8) | 0);
    }

    #[test]
    fn anchors_are_sorted_by_read_strand_and_query_position() {
        let minimizers = [
            Minimizer {
                hash: 42,
                read_identifier: 0,
                position: 100,
                reverse: false,
                span: 51,
            },
            Minimizer {
                hash: 43,
                read_identifier: 0,
                position: 50,
                reverse: true,
                span: 51,
            },
        ];

        let mut counts = crate::counting::KmerCountTable::new(51, 0, 0);
        for _ in 0..3 {
            counts.insert(42);
            counts.insert(43);
        }
        let mut index = PositionIndex::build(&counts, 51);
        index.insert(
            42,
            PositionRecord {
                read_identifier: 9,
                position: 300,
                reverse: false,
                span: 51,
            },
        );
        index.insert(
            43,
            PositionRecord {
                read_identifier: 9,
                position: 400,
                reverse: false,
                span: 51,
            },
        );

        let read_lengths = vec![1000u64; 10];
        let seeds = collect_qgen_seeds(&minimizers, &index, &read_lengths, 60, 12);

        assert_eq!(seeds.len(), 2);
        // The forward-strand anchor sorts before the reverse-strand one
        // because the strand bit sits above the query position in `srt`.
        assert_eq!(seeds[0].self_offset, 100);
        assert_eq!(seeds[1].self_offset, 50);
        assert!(seeds.iter().all(|seed| seed.read_identifier == 9));
        assert_eq!(seeds[0].reverse_strand, false, "both strands equal -> forward");
        assert_eq!(seeds[1].reverse_strand, true, "strands differ -> reverse");
    }
}

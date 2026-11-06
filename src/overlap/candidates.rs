//! Candidate seeds for one read, the counterpart of `ha_get_new_candidates`
//! in `anchor.cpp`.
//!
//! Every minimizer of the query is looked up in the position index; each hit
//! becomes an anchor that records where it sits on both reads. Anchors are
//! sorted by read, strand and target position so that the chaining code sees
//! the seeds of one candidate read in target order, and the seed support is
//! then turned into a weight that rewards seeds whose k-mer occurs a plausible
//! number of times.

use crate::container::ksort::radix_sort_by_key;
use crate::counting::minimizer::Minimizer;
use crate::counting::position_index::PositionIndex;
use crate::overlap::region::KmerHit;

/// `HA_KMER_GOOD_RATIO`: seeds below this fraction of the coverage are
/// considered too rare, above `2 - ratio` too repetitive.
const KMER_GOOD_RATIO: f64 = 0.333;

/// One anchor before it is turned into a `k_mer_hit`, `anchor1_t`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Anchor {
    /// `read_identifier << 33 | reverse << 32 | other_offset`.
    sort_key: u64,
    /// Position on the query read.
    self_offset: u32,
    /// Position on the target read.
    other_offset: u32,
    /// Number of positions this k-mer has in the index.
    count: u32,
}

/// Builds the seed list of one query read, mirroring `ha_get_new_candidates`.
///
/// `minimizers` are the query minimizers (already filtered by the
/// high-occurrence table), `index` the position index and `read_length` the
/// query length.
pub fn collect_seeds(
    minimizers: &[Minimizer],
    index: &PositionIndex,
    read_length: u32,
    homozygous_coverage: i32,
) -> Vec<KmerHit> {
    let mut anchors: Vec<Anchor> = Vec::new();

    for minimizer in minimizers {
        let positions = index.get(minimizer.hash);
        let position_count = positions.len() as u32;

        for position in positions {
            let reverse = minimizer.reverse != position.reverse;
            let self_offset = if reverse {
                // On the opposite strand the anchor is mirrored so that both
                // reads run in the same direction.
                read_length - 1 - (minimizer.position + 1 - minimizer.span)
            } else {
                minimizer.position
            };

            anchors.push(Anchor {
                sort_key: ((position.read_identifier as u64) << 33)
                    | ((reverse as u64) << 32)
                    | position.position as u64,
                self_offset,
                other_offset: position.position,
                count: position_count,
            });
        }
    }

    // Sort by read, strand and target position, then by query position inside
    // each equal-key run, so the chaining loop sees a stable anchor order.
    radix_sort_by_key(&mut anchors, &|anchor: &Anchor| anchor.sort_key);
    let mut run_start = 0usize;
    for index in 1..=anchors.len() {
        if index == anchors.len() || anchors[index].sort_key != anchors[run_start].sort_key {
            if index - run_start > 1 {
                radix_sort_by_key(&mut anchors[run_start..index], &|anchor: &Anchor| {
                    anchor.self_offset as u64
                });
            }
            run_start = index;
        }
    }

    anchors
        .iter()
        .map(|anchor| KmerHit {
            read_identifier: (anchor.sort_key >> 33) as u32,
            reverse_strand: (anchor.sort_key >> 32) & 1 == 1,
            offset: anchor.other_offset,
            self_offset: anchor.self_offset,
            count: seed_weight(anchor.count, homozygous_coverage),
        })
        .collect()
}

/// The `cnt` weighting applied when anchors become `k_mer_hit` entries.
fn seed_weight(position_count: u32, homozygous_coverage: i32) -> u32 {
    let low_occurrence = (homozygous_coverage as f64 * KMER_GOOD_RATIO) as u32;
    let high_occurrence = (homozygous_coverage as f64 * (2.0 - KMER_GOOD_RATIO)) as u32;

    if position_count > low_occurrence && position_count < high_occurrence {
        1
    } else if position_count <= low_occurrence {
        2
    } else {
        let doubled_high = (high_occurrence << 1).max(1);
        let scaled = 1 + (position_count + doubled_high - 1) / doubled_high;
        (scaled as f64).powf(1.1) as u32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seed_weight_rewards_rare_and_penalises_repetitive_seeds() {
        // With coverage 36: low = 11, high = 60.
        assert_eq!(seed_weight(5, 36), 2, "rare seeds are boosted");
        assert_eq!(seed_weight(30, 36), 1, "typical seeds keep weight one");
        // 100 positions: 1 + ceil(100 / 120) = 2, then 2^1.1 = 2.14 -> 2.
        assert_eq!(seed_weight(100, 36), 2);
        // 200 positions: 1 + ceil(200 / 120) = 3, then 3^1.1 = 3.34 -> 3.
        assert_eq!(seed_weight(200, 36), 3);
    }

    #[test]
    fn seeds_mirror_the_query_position_on_the_opposite_strand() {
        let minimizer = Minimizer {
            hash: 42,
            read_identifier: 7,
            position: 100,
            reverse: false,
            span: 51,
        };
        let mut counts = crate::counting::KmerCountTable::new(51, 0, 0);
        // Two occurrences keep the k-mer, and give it two positions.
        counts.insert(42);
        counts.insert(42);
        let mut index = crate::counting::position_index::PositionIndex::build(&counts, 51);
        index.insert(
            42,
            crate::counting::position_index::PositionRecord {
                read_identifier: 9,
                position: 300,
                reverse: false,
                span: 51,
            },
        );
        index.insert(
            42,
            crate::counting::position_index::PositionRecord {
                read_identifier: 9,
                position: 400,
                reverse: false,
                span: 51,
            },
        );

        let seeds = collect_seeds(&[minimizer], &index, 1000, 36);
        assert_eq!(seeds.len(), 2);
        assert!(seeds.iter().all(|seed| seed.read_identifier == 9));
        assert!(seeds.iter().all(|seed| !seed.reverse_strand));
        assert_eq!(seeds[0].self_offset, 100);
        assert_eq!(seeds[0].offset, 300);
        assert_eq!(seeds[1].offset, 400);

        let reverse_minimizer = Minimizer {
            reverse: true,
            ..minimizer
        };
        let seeds = collect_seeds(&[reverse_minimizer], &index, 1000, 36);
        assert!(seeds.iter().all(|seed| seed.reverse_strand));
        // 1000 - 1 - (100 + 1 - 51) = 949.
        assert_eq!(seeds[0].self_offset, 949);
    }
}

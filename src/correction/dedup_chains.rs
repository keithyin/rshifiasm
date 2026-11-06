//! `dedup_chains` (`ecovlp.cpp:2984`): keeps one overlap per target read.
//!
//! After the whole-overlap refinement, the chain list of a read may hold
//! several overlaps to the same target. The reference sorts the list by
//! `y_id`, then, for every group of two or more, keeps the entry with the
//! smallest `is_match`, the largest score `(x_end + 1 - x_start) -
//! 12 * non_homopolymer_errors`, and the longest query span, in that order of
//! priority. The survivors are compacted to the front with the reference's
//! in-place swap and the list is truncated.
//!
//! The groups are formed after `overlap_region_sort_y_id`, an MSD radix sort
//! over the four-byte `y_id`; the port sorts an index permutation with the
//! same 32-bit radix sort, which yields the exact same permutation as sorting
//! the records themselves.

use crate::overlap::region::{sort_by_y_identifier, OverlapRegion};

/// `dedup_chains`: compacts the list to one overlap per target read.
///
/// The list is sorted in place and truncated to the number of surviving
/// groups, exactly like the reference's `ol->length = m`.
pub fn dedup_chains(overlaps: &mut Vec<OverlapRegion>) {
    if overlaps.len() <= 1 {
        return;
    }

    sort_by_y_identifier(overlaps);

    let length = overlaps.len();
    let mut kept = 0usize;
    let mut group_start = 0usize;
    for group_end in 1..=length {
        if group_end == length
            || overlaps[group_end].y_identifier != overlaps[group_start].y_identifier
        {
            // The reference defaults to the group's first entry.
            let mut chosen = group_start;
            if group_end - group_start > 1 {
                let mut best_score: i64 = i64::from(i32::MIN);
                let mut best_match: u8 = 3;
                for index in group_start..group_end {
                    let region = &overlaps[index];
                    let plus = region.x_position_end.wrapping_add(1).wrapping_sub(region.x_position_start);
                    let minus = region.non_homopolymer_errors.wrapping_mul(12);
                    let score = i64::from(plus) - i64::from(minus);

                    let mut select = false;
                    if region.is_match < best_match {
                        select = true;
                    } else if region.is_match == best_match {
                        let chosen_span = overlaps[chosen]
                            .x_position_end
                            .wrapping_add(1)
                            .wrapping_sub(overlaps[chosen].x_position_start);
                        if score > best_score
                            || (score == best_score && plus > chosen_span)
                        {
                            select = true;
                        }
                    }
                    if select {
                        best_score = score;
                        chosen = index;
                        best_match = region.is_match;
                    }
                }
            }

            if chosen != kept {
                overlaps.swap(chosen, kept);
            }
            kept += 1;
            group_start = group_end;
        }
    }

    overlaps.truncate(kept);
}

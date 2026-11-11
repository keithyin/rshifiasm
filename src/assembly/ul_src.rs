//! `ul_load` post-processing: `update_ovlp_src` / `update_ovlp_src_bl`
//! (`inter.cpp` 19916 / 19976) label the reliable UL primary-chain linkages
//! onto the HiFi source overlaps (`ma_hit_t::bl`).
//!
//! For every HiFi read the UL read index (`ridx`) lists the UL blocks that
//! contain it. `determine_connective_backtrack` walks each such block's
//! primary chain (`pchain`/`pidx`) backwards and, for consecutive *full-length*
//! blocks, finds the matching HiFi source overlap and increments its `bl`.
//! `update_ovlp_src_bl` then symmetrises `bl` between the two directions.

use crate::assembly::hit::ma_hit_to_arc;
use crate::assembly::ul_store::{AllUl, UlReadVec};
use crate::correction::ec_driver::PafList;
use crate::reads::read_store::ReadStore;

/// The `ug_opt_t` slice the source labelling consumes.
#[derive(Debug, Clone, Copy)]
pub struct UlSourceOptions {
    pub maximum_hang: i32,
    pub maximum_hang_rate: f32,
    pub minimum_overlap: i32,
}

/// `update_ovlp_src` (`inter.cpp` 19916): zero every `bl`, then label the
/// reliable UL chains back onto the HiFi source overlaps.
pub fn update_ovlp_src(
    store: &AllUl,
    sources: &mut [PafList],
    reads: &ReadStore,
    options: &UlSourceOptions,
) {
    for read in sources.iter_mut() {
        for hit in &mut read.hits {
            hit.block_length = 0;
        }
    }
    for read_id in 0..reads.total_reads as usize {
        let start = store.ridx.idx[read_id] as usize;
        let end = store.ridx.idx[read_id + 1] as usize;
        for occurrence in &store.ridx.occ[start..end] {
            let ul_id = (occurrence >> 32) as usize;
            let block_index = *occurrence as u32;
            determine_connective_backtrack(
                &store.reads[ul_id],
                block_index,
                read_id as u32,
                sources,
                reads,
                options,
            );
        }
    }
}

/// `determine_connective_backtrack` (`inter.cpp` 19874): one reliable UL block
/// labels its primary-chain neighbours' HiFi source overlaps.
fn determine_connective_backtrack(
    ul_read: &UlReadVec,
    block_index: u32,
    read_id: u32,
    sources: &mut [PafList],
    reads: &ReadStore,
    options: &UlSourceOptions,
) {
    let li = &ul_read.bb[block_index as usize];
    debug_assert!(!li.base && li.hid == read_id && li.el);
    if block_index == 0 {
        return;
    }
    if li.pchain == 0 {
        return;
    }
    if li.pidx == u32::MAX {
        return;
    }
    let li_v = (((li.hid) << 1) | u32::from(li.rev)) ^ 1;
    if (li.te - li.ts) < reads.read_length(li.hid as usize) as u32 {
        return;
    }
    let li_hid = li_v >> 1;
    let source = &mut sources[li_hid as usize];

    let mut lk_index = li.pidx;
    loop {
        let lk = &ul_read.bb[lk_index as usize];
        let lk_v = (((lk.hid) << 1) | u32::from(lk.rev)) ^ 1;
        debug_assert!(!lk.base && lk.pchain != 0 && lk.el);
        if (lk.te - lk.ts) >= reads.read_length(lk.hid as usize) as u32 {
            let lk_hid = lk_v >> 1;
            let mut hit_index = 0usize;
            let hit_count = source.hits.len();
            while hit_index < hit_count {
                let hit = &source.hits[hit_index];
                if hit.query_identifier == li_hid && hit.target_identifier == lk_hid {
                    let (result, arc) = ma_hit_to_arc(
                        hit,
                        reads.read_length(li_hid as usize) as i32,
                        reads.read_length(lk_hid as usize) as i32,
                        options.maximum_hang,
                        options.maximum_hang_rate,
                        options.minimum_overlap,
                    );
                    if result >= 0
                        && (arc.ul >> 32) as u32 == li_v
                        && arc.v == lk_v
                    {
                        break;
                    }
                }
                hit_index += 1;
            }
            if hit_index < hit_count {
                source.hits[hit_index].block_length += 1;
            }
        }
        lk_index = if lk.pidx == u32::MAX {
            break;
        } else {
            lk.pidx
        };
    }
}

/// `update_ovlp_src_bl` (`inter.cpp` 19976): make `bl` symmetric, so both
/// directions of one HiFi overlap carry the same label (the sum).
pub fn update_ovlp_src_bl(sources: &mut [PafList], reads: &ReadStore) {
    for read_id in 0..reads.total_reads as usize {
        let hit_count = sources[read_id].hits.len();
        for hit_index in 0..hit_count {
            let (query, target, block_length) = {
                let hit = &sources[read_id].hits[hit_index];
                (hit.query_identifier, hit.target_identifier, hit.block_length)
            };
            if query > target {
                continue;
            }
            // `get_specific_overlap(&src[tn], tn, qn)`: the reverse hit.
            let reverse = sources[target as usize]
                .hits
                .iter()
                .position(|hit| hit.query_identifier == target && hit.target_identifier == query);
            let Some(reverse_index) = reverse else {
                continue;
            };
            let sum = block_length + sources[target as usize].hits[reverse_index].block_length;
            sources[read_id].hits[hit_index].block_length = sum;
            sources[target as usize].hits[reverse_index].block_length = sum;
        }
    }
}

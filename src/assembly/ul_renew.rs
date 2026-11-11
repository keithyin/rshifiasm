//! The `--ul` renewal cluster (`inter.cpp` 19557–23020, `gfa_ut.cpp`
//! 17598–17700): the second and later ultra-long passes. The read graph is
//! turned into a polished unitig graph (`gen_polished_ug`), the UL reads are
//! re-aligned against it (`ul_realignment`/`gen_UL_reovlps`/`rescall`), the
//! resulting block chains are filtered (`filter_ul_ug`) and labelled onto the
//! unitig arcs (`update_ug_arch_ul_mul`), and finally the read graph is
//! rebuilt from the UL-supported paths (`ul_realignment_gfa` → `renew_ng`).
//!
//! This module starts with the self-contained pieces; the rescall pipeline
//! and the integer/UL-vs-UL rebuilding live in the sibling modules.

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::hit::CoverageCut;
use crate::assembly::unitig::{ma_ug_gen, polish_unitig_advance_scaf, polish_unitig_scaf, UnitigGraph};
use crate::assembly::ul::UcBlock;
use crate::assembly::ul_store::{AllUl, UlReadVec};
use crate::command_line::options::AssemblyOptions;
use crate::correction::ec_driver::PafList;
use crate::reads::ReadStore;

/// `ugl_cover_check` (`inter.cpp` 19600): whether `[is, ie)` covers one whole
/// read of the unitig.
pub fn ugl_cover_check(is: u64, ie: u64, unitig: &MaUtg, reads: &ReadStore) -> bool {
    if is == 0 && ie == u64::from(unitig.length) {
        return true;
    }
    let mut l: u64 = 0;
    for &element in &unitig.reads {
        let us = l;
        let ue = l + reads.read_length((element >> 33) as usize);
        if is <= us && ie >= ue {
            return true;
        }
        if us >= ie {
            break;
        }
        l += (element & 0xFFFF_FFFF) as u64;
    }
    false
}

/// `clear_all_ul_t` (`inter.cpp` 22004): empties every block chain and the
/// read index while keeping the read lengths.
pub fn clear_all_ul_t(store: &mut AllUl) {
    for read in store.reads.iter_mut() {
        let rlen = read.rlen;
        *read = UlReadVec {
            rlen,
            ..UlReadVec::default()
        };
    }
    store.ridx = Default::default();
}

/// `gen_polished_ug` (`Overlaps.cpp` 10429): builds the unitig graph of `sg`
/// and polishes every unitig's consensus, then recomputes the arc lengths.
pub fn gen_polished_ug(
    options: &AssemblyOptions,
    reads: &ReadStore,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    sg: &mut Asg,
) -> UnitigGraph {
    let mut unitig_graph = ma_ug_gen(sg);
    let edges: Vec<AsgArc> = Vec::new();
    let mut new_edges: Vec<AsgArc> = Vec::new();
    let maximum_hang = options.maximum_hang_length;
    let minimum_overlap = options.minimum_overlap_length;
    let max_hang_rate = options.maximum_hang_rate;

    for unitig in unitig_graph.unitigs.iter_mut() {
        if unitig.reads.is_empty() {
            continue;
        }
        polish_unitig_scaf(
            unitig,
            sg,
            sources,
            coverage_cut,
            &edges,
            Some(&mut new_edges),
            maximum_hang,
            minimum_overlap,
            max_hang_rate,
        );
        polish_unitig_advance_scaf(
            unitig,
            sg,
            sources,
            coverage_cut,
            &edges,
            Some(&mut new_edges),
            reads,
            maximum_hang,
            minimum_overlap,
            max_hang_rate,
        );
    }
    for (i, unitig) in unitig_graph.unitigs.iter().enumerate() {
        unitig_graph.graph.sequences[i].length = unitig.length;
    }

    let vertex_count = unitig_graph.graph.sequences.len() * 2;
    for v in 0..vertex_count as u32 {
        if unitig_graph.graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        let (start, count) = unitig_graph.graph.arc_range(v);
        for i in 0..count {
            let arc = unitig_graph.graph.arcs[start + i];
            if arc.deleted {
                continue;
            }
            let source = (arc.ul >> 33) as usize;
            let node_length = unitig_graph.graph.sequences[source]
                .length
                .wrapping_sub(arc.overlap);
            unitig_graph.graph.arcs[start + i].ul = ((arc.ul >> 32) << 32) | u64::from(node_length);
        }
    }
    unitig_graph
}

/// `filter_short_ulalignments` (`inter.cpp` 19634) for one read: drops the
/// chain links whose blocks do not cover a whole unitig.
fn filter_short_ulalignments(unitig_graph: &UnitigGraph, read: &mut UlReadVec, reads: &ReadStore) {
    let a_n = read.bb.len();
    if a_n == 0 {
        return;
    }
    for k in (0..a_n).rev() {
        let p = read.bb[k];
        if p.base || !p.el || p.pchain == 0 {
            continue;
        }
        let unitig = &unitig_graph.unitigs[p.hid as usize];
        if p.pidx == u32::MAX {
            if !ugl_cover_check(u64::from(p.ts), u64::from(p.te), unitig, reads) {
                read.bb[k].pchain = 0;
                if p.aidx != u32::MAX {
                    read.bb[p.aidx as usize].pidx = u32::MAX;
                    read.bb[p.aidx as usize].pdis = u32::MAX;
                    read.bb[k].aidx = u32::MAX;
                }
            }
            continue;
        }
        if ugl_cover_check(u64::from(p.ts), u64::from(p.te), unitig, reads) {
            continue;
        }
        let mut z = p.pidx;
        while z != u32::MAX {
            if ugl_cover_check(
                u64::from(read.bb[z as usize].ts),
                u64::from(read.bb[z as usize].te),
                &unitig_graph.unitigs[read.bb[z as usize].hid as usize],
                reads,
            ) {
                break;
            }
            z = read.bb[z as usize].pidx;
        }
        let lz = z;
        let fz = p.aidx;
        let mut l: u64 = 0;
        if fz != u32::MAX {
            l = u64::from(read.bb[fz as usize].pdis);
        }
        let mut z = k as u32;
        while z != lz {
            let bz = read.bb[z as usize].pidx;
            l += u64::from(read.bb[z as usize].pdis);
            read.bb[z as usize].pidx = u32::MAX;
            read.bb[z as usize].pdis = u32::MAX;
            read.bb[z as usize].aidx = u32::MAX;
            read.bb[z as usize].pchain = 0;
            z = bz;
        }
        if fz != u32::MAX && lz != u32::MAX {
            read.bb[fz as usize].pdis = l as u32;
            read.bb[fz as usize].pidx = lz;
            read.bb[lz as usize].aidx = fz;
        } else if fz != u32::MAX {
            read.bb[fz as usize].pdis = u32::MAX;
            read.bb[fz as usize].pidx = u32::MAX;
        } else if lz != u32::MAX {
            read.bb[lz as usize].aidx = u32::MAX;
        }
    }

    for k in (0..a_n).rev() {
        let p = read.bb[k];
        if p.base || !p.el || p.pchain == 0 {
            continue;
        }
        if p.pidx != u32::MAX {
            let parent = read.bb[p.pidx as usize];
            if parent.aidx != k as u32 && parent.aidx == u32::MAX {
                read.bb[k].pidx = u32::MAX;
            } else {
                debug_assert_eq!(parent.aidx, k as u32);
                debug_assert!(parent.pchain != 0);
            }
        }
        if p.aidx != u32::MAX {
            let child = read.bb[p.aidx as usize];
            if child.pidx != k as u32 && child.pidx == u32::MAX {
                read.bb[k].aidx = u32::MAX;
            } else {
                debug_assert_eq!(child.pidx, k as u32);
                debug_assert!(child.pchain != 0);
            }
        }
    }
}

/// `filter_ul_ug` (`inter.cpp` 19764): the per-read pass above over every UL
/// read.
pub fn filter_ul_ug(store: &mut AllUl, unitig_graph: &UnitigGraph, reads: &ReadStore) {
    for read in store.reads.iter_mut() {
        filter_short_ulalignments(unitig_graph, read, reads);
    }
}

/// `update_ug_arch_ul` (`inter.cpp` 19613) for one arc: counts the UL chains
/// that pass through this unitig arc.
fn update_ug_arch_ul(
    unitig_graph: &mut UnitigGraph,
    store: &AllUl,
    arc_index: usize,
    ou: &mut u16,
) {
    let arc = unitig_graph.graph.arcs[arc_index];
    *ou = 0;
    let v = (arc.ul >> 32) as u32;
    let w = arc.v;
    let read_id = (v >> 1) as usize;
    let start = store.ridx.idx[read_id] as usize;
    let end = store.ridx.idx[read_id + 1] as usize;
    for &a in &store.ridx.occ[start..end] {
        let ul_read = &store.reads[(a >> 32) as usize];
        let p = &ul_read.bb[a as u32 as usize];
        if p.base || !p.el || p.pchain == 0 {
            continue;
        }
        let uv = (p.hid << 1) | u32::from(p.rev);
        if uv == v && p.aidx != u32::MAX {
            let n = &ul_read.bb[p.aidx as usize];
            debug_assert!(!n.base && n.el && n.pchain != 0 && n.pidx == a as u32);
            let uw = (n.hid << 1) | u32::from(n.rev);
            if uw == w {
                *ou = ou.saturating_add(1);
            }
        }
        if (uv ^ 1) == v && p.pidx != u32::MAX {
            let n = &ul_read.bb[p.pidx as usize];
            debug_assert!(!n.base && n.el && n.pchain != 0 && n.aidx == a as u32);
            let uw = ((n.hid << 1) | u32::from(n.rev)) ^ 1;
            if uw == w {
                *ou = ou.saturating_add(1);
            }
        }
    }
}

/// `update_ug_arch_ul_mul` (`inter.cpp` 19645): the whole-array version.
pub fn update_ug_arch_ul_mul(unitig_graph: &mut UnitigGraph, store: &AllUl) {
    for arc_index in 0..unitig_graph.graph.arcs.len() {
        let mut ou = 0u16;
        update_ug_arch_ul(unitig_graph, store, arc_index, &mut ou);
        unitig_graph.graph.arcs[arc_index].occurrence = ou;
    }
}

/// A convenience for the older call sites that still pass a bare read list.
pub fn filter_ul_ug_reads(store: &mut AllUl, unitig_graph: &UnitigGraph, reads: &ReadStore) {
    filter_ul_ug(store, unitig_graph, reads);
}

// The block type is re-exported for the callers that build chains by hand.
#[allow(unused_imports)]
use UcBlock as _UcBlock;

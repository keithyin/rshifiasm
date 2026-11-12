//! Unitig construction, polishing, printing and tip extension, the counterpart
//! of the `ma_ug_t` family in `Overlaps.cpp` and of `gen_ext_tip`/`ug_ext_gfa`
//! in `gfa_ut.cpp`.
//!
//! The reference keeps the unitig list (`ma_ug_v`) and the unitig string graph
//! (`asg_t`) in one `ma_ug_t` structure. [`UnitigGraph`] mirrors that with a
//! `Vec<MaUtg>` plus an [`Asg`]. Every helper that the reference hands a
//! pointer to (`R_to_U`, `ma_sub_t`, `ma_hit_t_alloc`) takes the Rust slice or
//! structure that owns the same data.
//!
//! One deliberate simplification: the reference's `edge` and `newE` parameters
//! of `ma_ug_seq` and of the polishing helpers are always the same
//! `new_rtg_edges` vector on the paths ported here, so the Rust functions take
//! one `edges: &mut Vec<AsgArc>` that plays both roles.

use std::io::{self, Write};

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::bubble::{BubbleBuffer, BubbleInfo, BubbleLabel};
use crate::assembly::bubble_identify::{
    asg_bub_pop1_label, asg_bub_pop1_primary_trio, asg_pop_bubble_primary_trio,
    get_s_bub_pop_max_dist_advance, get_unitig,
};
use crate::assembly::gfa_clean::{
    get_real_length, get_real_length_one, LOOP, MUL_INPUT, MUL_OUTPUT,
};
use crate::assembly::hap_cov::reset_trans_chain;
use crate::assembly::hit::{ma_hit_to_arc, CoverageCut, ReadToUnitig, ALTER_LABEL, PRIMARY_LABEL};
use crate::assembly::purge_dups::ALTER_COVERAGE_THRESHOLD;
use crate::assembly::rescue::asg_arc_del_simple_circle_untig;
use crate::assembly::trio_clean::{
    asg_seq_drop, check_different_haps, drop_semi_circle, DROP, PLOID,
};
use crate::correction::ec_driver::PafList;
use crate::overlap::paf::PafHit;
use crate::reads::ReadStore;

/// `AMBIGU`: no trio assignment.
pub const AMBIGU: u8 = 0;
/// `FATHER`: the read was assigned to the father.
pub const FATHER: u8 = 1;
/// `MOTHER`: the read was assigned to the mother.
pub const MOTHER: u8 = 2;

/// The sentinel the reference writes into `ma_utg_t::a` for a removed element.
const DELETED_ELEMENT: u64 = u64::MAX;

/// `comp_tab`: the complement table used when a read joins a unitig on the
/// reverse strand. Bytes above 127 map to themselves here, but the caller
/// turns them into `N` before the lookup, exactly like the reference.
const COMPLEMENT_TABLE: [u8; 128] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
    25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
    48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, b'@', b'T', b'V', b'G',
    b'H', b'E', b'F', b'C', b'D', b'I', b'J', b'M', b'L', b'K', b'N', b'O', b'P', b'Q', b'Y',
    b'S', b'A', b'A', b'B', b'W', b'X', b'R', b'Z', b'[', b'\\', b']', b'^', b'_', b'`', b't',
    b'v', b'g', b'h', b'e', b'f', b'c', b'd', b'i', b'j', b'm', b'l', b'k', b'n', b'o', b'p',
    b'q', b'y', b's', b'a', b'a', b'b', b'w', b'x', b'r', b'z', b'{', b'|', b'}', b'~', 127,
];

/// The complement of one read byte, `comp_tab[c]`.
fn complement_byte(base: u8) -> u8 {
    if base >= 128 {
        b'N'
    } else {
        COMPLEMENT_TABLE[base as usize]
    }
}

/// The unitig list and unitig graph of the reference's `ma_ug_t`.
#[derive(Debug, Clone, Default)]
/// The unitig graph, `ma_ug_t`.
pub struct UnitigGraph {
    /// `ma_ug_t::u`.
    pub unitigs: Vec<MaUtg>,
    /// `ma_ug_t::g`.
    pub graph: Asg,
}

/// `ma_ug_gen`: groups the read graph into linear (or circular) unitigs.
pub fn ma_ug_gen(graph: &mut Asg) -> UnitigGraph {
    ma_ug_gen_inner(graph, None)
}

/// `ma_ug_gen_primary`: [`ma_ug_gen`] with the coverage-class filter.
///
/// `flag` is `PRIMARY_LABLE` or `ALTER_LABLE`. Every `asg_seq_t::c` of the
/// default HiFi pipeline is `PRIMARY_LABLE`, so the filter rejects nothing.
pub fn ma_ug_gen_primary(graph: &mut Asg, flag: u8) -> UnitigGraph {
    ma_ug_gen_inner(graph, Some(flag))
}

fn ma_ug_gen_inner(graph: &mut Asg, coverage_flag: Option<u8>) -> UnitigGraph {
    graph.cleanup();

    let vertex_count = graph.sequences.len() * 2;
    let mut mark = vec![0i32; vertex_count];
    let mut unitigs: Vec<MaUtg> = Vec::new();
    let mut queue: Vec<u64> = Vec::new();

    for v in 0..vertex_count as u32 {
        if graph.sequences[(v >> 1) as usize].deleted || mark[v as usize] != 0 {
            continue;
        }
        if graph.arc_count(v) == 0 && graph.arc_count(v ^ 1) != 0 {
            continue;
        }
        if let Some(flag) = coverage_flag {
            let class = graph.sequences[(v >> 1) as usize].coverage_class;
            if flag == PRIMARY_LABEL && class == ALTER_LABEL {
                continue;
            }
            if flag == ALTER_LABEL && class != ALTER_LABEL {
                continue;
            }
        }

        mark[v as usize] = 1;
        queue.clear();
        let mut start = v;
        let mut end = v ^ 1;
        let mut length: u32 = 0;

        // Forward walk: the only suffix of `w` is `x`, and the only prefix of
        // `x` is `w`.
        let mut w = v;
        loop {
            if graph.arc_count(w) != 1 {
                break;
            }
            let x = graph.arcs_of(w)[0].v;
            if graph.arc_count(x ^ 1) != 1 {
                break;
            }
            mark[x as usize] = 1;
            mark[(w ^ 1) as usize] = 1;
            let edge_length = graph.arcs_of(w)[0].node_length();
            queue.push((u64::from(w) << 32) | u64::from(edge_length));
            end = x ^ 1;
            length = length.wrapping_add(edge_length);
            w = x;
            if x == v {
                break;
            }
        }

        let circular;
        if start != (end ^ 1) || queue.is_empty() {
            // Linear unitig: the last element contributes its sequence length.
            let sequence_length = graph.sequences[(end >> 1) as usize].length;
            queue.push((u64::from(end ^ 1) << 32) | u64::from(sequence_length));
            length = length.wrapping_add(sequence_length);
            circular = false;
        } else {
            // Circular unitig: no backward walk, no end vertices.
            start = u32::MAX;
            end = u32::MAX;
            circular = true;
        }

        if !circular {
            // Backward walk.
            let mut x = v;
            loop {
                if graph.arc_count(x ^ 1) != 1 {
                    break;
                }
                let w = graph.arcs_of(x ^ 1)[0].v ^ 1;
                if graph.arc_count(w) != 1 {
                    break;
                }
                mark[x as usize] = 1;
                mark[(w ^ 1) as usize] = 1;
                let edge_length = graph.arcs_of(w)[0].node_length();
                // The reference unshifts onto the front of the deque, so the
                // backward chain ends up in front of the forward chain.
                queue.insert(0, (u64::from(w) << 32) | u64::from(edge_length));
                start = w;
                length = length.wrapping_add(edge_length);
                x = w;
            }
        }

        if start != u32::MAX {
            mark[start as usize] = 1;
            mark[end as usize] = 1;
        }

        unitigs.push(MaUtg {
            length: length & 0x7fff_ffff,
            circular,
            start,
            end,
            reads: queue.iter().copied().collect(),
            sequence: Vec::new(),
        });
    }

    // Add the arcs between unitigs, reusing `mark` for the endpoint map.
    mark.fill(-1);
    for (i, unitig) in unitigs.iter().enumerate() {
        if unitig.circular {
            continue;
        }
        mark[unitig.start as usize] = ((i as i32) << 1) | 0;
        mark[unitig.end as usize] = ((i as i32) << 1) | 1;
    }

    let mut unitig_graph = Asg::new();
    for arc in &graph.arcs {
        if arc.deleted {
            continue;
        }
        // `p->ul>>32^1` binds as `(p->ul>>32)^1` in the reference.
        let query_mark = mark[(arc.query_vertex() ^ 1) as usize];
        let target_mark = mark[arc.v as usize];
        if query_mark >= 0 && target_mark >= 0 {
            let u = (query_mark as u32) ^ 1;
            let mut node_length = unitigs[(u >> 1) as usize]
                .length
                .wrapping_sub(arc.overlap) as i32;
            if node_length < 0 {
                node_length = 1;
            }
            unitig_graph.push_arc(AsgArc {
                ul: (u64::from(u) << 32) | u64::from(node_length as u32),
                v: target_mark as u32,
                overlap: arc.overlap,
                deleted: false,
                occurrence: 0,
                strong: false,
                without_large_indel: false,
                exact_region: arc.exact_region,
            });
        }
    }
    for (i, unitig) in unitigs.iter().enumerate() {
        unitig_graph.set_sequence(i, unitig.length, false);
    }
    unitig_graph.cleanup();

    UnitigGraph {
        unitigs,
        graph: unitig_graph,
    }
}

/// `PHASE_SEF`: the multiplier that turns the homozygous coverage into the
/// minimum occurrence of a phased unitig.
pub const PHASE_SEF: f64 = 2.0;
/// `PHASE_SEP`: the lower bound of the phase occurrence threshold.
pub const PHASE_SEP: u64 = 64;
/// `PHASE_SEP_RATE`: the ratio that decides whether a chain is split by the
/// trio labels.
pub const PHASE_SEP_RATE: f64 = 0.04;

/// `ma_ug_gen_phase`: [`ma_ug_gen`] with trio-label aware unitig splitting.
///
/// A chain whose reads carry a clear father or mother majority is split into
/// phase blocks at the label switches; the remaining chains stay whole.
pub fn ma_ug_gen_phase(
    graph: &mut Asg,
    min_occ: u32,
    cutoff: f64,
    trio_labels: &[u8],
) -> UnitigGraph {
    graph.cleanup();

    let normalize = |label: u8| -> u8 {
        if label == FATHER || label == MOTHER {
            label
        } else {
            AMBIGU
        }
    };
    let label_of = |packed: u64| -> u8 { normalize(trio_labels[(packed >> 33) as usize]) };

    let vertex_count = graph.sequences.len() * 2;
    let mut mark = vec![0i32; vertex_count];
    let mut unitigs: Vec<MaUtg> = Vec::new();
    let mut queue: Vec<u64> = Vec::new();
    let mut uidx: Vec<u64> = Vec::new();

    for v in 0..vertex_count as u32 {
        if graph.sequences[(v >> 1) as usize].deleted || mark[v as usize] != 0 {
            continue;
        }
        if graph.arc_count(v) == 0 && graph.arc_count(v ^ 1) != 0 {
            continue;
        }

        mark[v as usize] = 1;
        queue.clear();
        let mut start = v;
        let mut end = v ^ 1;
        let mut length: u32 = 0;
        let mut father_count: u32 = 0;
        let mut mother_count: u32 = 0;

        // Forward walk, exactly as in `ma_ug_gen`.
        let mut w = v;
        loop {
            if graph.arc_count(w) != 1 {
                break;
            }
            let x = graph.arcs_of(w)[0].v;
            if graph.arc_count(x ^ 1) != 1 {
                break;
            }
            mark[x as usize] = 1;
            mark[(w ^ 1) as usize] = 1;
            let edge_length = graph.arcs_of(w)[0].node_length();
            queue.push((u64::from(w) << 32) | u64::from(edge_length));
            match trio_labels[(w >> 1) as usize] {
                FATHER => father_count = father_count.wrapping_add(1),
                MOTHER => mother_count = mother_count.wrapping_add(1),
                _ => {}
            }
            end = x ^ 1;
            length = length.wrapping_add(edge_length);
            w = x;
            if x == v {
                break;
            }
        }

        let circular;
        if start != (end ^ 1) || queue.is_empty() {
            let sequence_length = graph.sequences[(end >> 1) as usize].length;
            queue.push((u64::from(end ^ 1) << 32) | u64::from(sequence_length));
            match trio_labels[(end >> 1) as usize] {
                FATHER => father_count = father_count.wrapping_add(1),
                MOTHER => mother_count = mother_count.wrapping_add(1),
                _ => {}
            }
            length = length.wrapping_add(sequence_length);
            circular = false;
        } else {
            start = u32::MAX;
            end = u32::MAX;
            circular = true;
        }

        if !circular {
            let mut x = v;
            loop {
                if graph.arc_count(x ^ 1) != 1 {
                    break;
                }
                let w = graph.arcs_of(x ^ 1)[0].v ^ 1;
                if graph.arc_count(w) != 1 {
                    break;
                }
                mark[x as usize] = 1;
                mark[(w ^ 1) as usize] = 1;
                let edge_length = graph.arcs_of(w)[0].node_length();
                queue.insert(0, (u64::from(w) << 32) | u64::from(edge_length));
                match trio_labels[(w >> 1) as usize] {
                    FATHER => father_count = father_count.wrapping_add(1),
                    MOTHER => mother_count = mother_count.wrapping_add(1),
                    _ => {}
                }
                start = w;
                length = length.wrapping_add(edge_length);
                x = w;
            }
        }

        if start != u32::MAX {
            mark[start as usize] = 1;
            mark[end as usize] = 1;
        }

        let common = father_count.min(mother_count);
        let is_phased = common > min_occ
            && f64::from(common) > f64::from(father_count.wrapping_add(mother_count)) * cutoff;
        if !is_phased || queue.is_empty() {
            unitigs.push(MaUtg {
                length: length & 0x7fff_ffff,
                circular,
                start,
                end,
                reads: queue.clone(),
                sequence: Vec::new(),
            });
            continue;
        }

        // Split the traversal into phase blocks at the label switches.
        uidx.clear();
        let queue_size = queue.len();
        let current = label_of(queue[0]);
        let mut previous = current;
        let first_father = father_count;
        let first_mother = mother_count;
        father_count = 0;
        mother_count = 0;
        if current == FATHER {
            father_count = 1;
        }
        if current == MOTHER {
            mother_count = 1;
        }
        let mut segment_start = 0usize;
        for k in 1..=queue_size {
            let mut split = false;
            let mut next = AMBIGU;
            if k == queue_size {
                split = true;
            } else {
                next = label_of(queue[k]);
                if next != AMBIGU && previous != AMBIGU && next != previous {
                    split = true;
                }
            }
            if split {
                let mut packed = (segment_start as u64) << 1;
                packed |= u64::from(father_count.max(mother_count)) << 32;
                if mother_count != 0 {
                    packed |= 1;
                }
                uidx.push(packed);
                father_count = 0;
                mother_count = 0;
                segment_start = k;
            }
            if next != AMBIGU {
                previous = next;
            }
            if next == FATHER {
                father_count = father_count.wrapping_add(1);
            }
            if next == MOTHER {
                mother_count = mother_count.wrapping_add(1);
            }
        }

        let mut father: u32 = 0;
        let mut mother: u32 = 0;
        let mut piece_father: u32 = 0;
        let mut piece_mother: u32 = 0;
        let mut emitted_reads: u32 = 0;
        let first = uidx[0];
        if first & 1 != 0 {
            mother = mother.wrapping_add((first >> 32) as u32);
        } else {
            father = father.wrapping_add((first >> 32) as u32);
        }
        let mut block_start = 0usize;
        for k in 1..=uidx.len() {
            let mut split = false;
            if k == uidx.len() {
                split = true;
            } else {
                let entry = uidx[k];
                if entry & 1 != 0 {
                    mother = mother.wrapping_add((entry >> 32) as u32);
                } else {
                    father = father.wrapping_add((entry >> 32) as u32);
                }
                let common = father.min(mother);
                if common > min_occ
                    && f64::from(common) > f64::from(father.wrapping_add(mother)) * cutoff
                {
                    split = true;
                }
            }
            if !split {
                continue;
            }

            let sz = ((uidx[block_start] as u32) >> 1) as usize;
            let ez = if k < uidx.len() {
                ((uidx[k] as u32) >> 1) as usize
            } else {
                queue_size
            };
            assert!(ez > sz);
            emitted_reads += (ez - sz) as u32;

            if circular && sz == 0 && ez == queue_size {
                let reads = queue.clone();
                for &packed in &reads {
                    let label = label_of(packed);
                    if label == FATHER {
                        piece_father = piece_father.wrapping_add(1);
                    }
                    if label == MOTHER {
                        piece_mother = piece_mother.wrapping_add(1);
                    }
                }
                unitigs.push(MaUtg {
                    length: length & 0x7fff_ffff,
                    circular: true,
                    start,
                    end,
                    reads,
                    sequence: Vec::new(),
                });
            } else {
                let mut piece_length: u32 = 0;
                let mut reads: Vec<u64> = Vec::with_capacity(ez - sz);
                for i in sz..ez - 1 {
                    let packed = queue[i];
                    reads.push(packed);
                    piece_length = piece_length.wrapping_add(packed as u32);
                    let label = label_of(packed);
                    if label == FATHER {
                        piece_father = piece_father.wrapping_add(1);
                    }
                    if label == MOTHER {
                        piece_mother = piece_mother.wrapping_add(1);
                    }
                }
                // The last element carries the sequence length of its vertex.
                let vertex = (queue[ez - 1] >> 32) as u32;
                let last = (u64::from(vertex) << 32)
                    | u64::from(graph.sequences[(vertex >> 1) as usize].length);
                reads.push(last);
                piece_length = piece_length.wrapping_add(last as u32);
                let label = label_of(last);
                if label == FATHER {
                    piece_father = piece_father.wrapping_add(1);
                }
                if label == MOTHER {
                    piece_mother = piece_mother.wrapping_add(1);
                }
                unitigs.push(MaUtg {
                    length: piece_length & 0x7fff_ffff,
                    circular: false,
                    start: (queue[sz] >> 32) as u32,
                    end: ((queue[ez - 1] >> 32) as u32) ^ 1,
                    reads,
                    sequence: Vec::new(),
                });
            }

            father = 0;
            mother = 0;
            block_start = k;
            if k < uidx.len() {
                let entry = uidx[k];
                if entry & 1 != 0 {
                    mother = mother.wrapping_add((entry >> 32) as u32);
                } else {
                    father = father.wrapping_add((entry >> 32) as u32);
                }
            }
        }
        assert_eq!(emitted_reads as usize, queue_size);
        assert_eq!(piece_father, first_father);
        assert_eq!(piece_mother, first_mother);
    }

    // Add the arcs between unitigs.
    mark.fill(-1);
    for (i, unitig) in unitigs.iter().enumerate() {
        if unitig.circular {
            continue;
        }
        mark[unitig.start as usize] = ((i as i32) << 1) | 0;
        mark[unitig.end as usize] = ((i as i32) << 1) | 1;
    }

    let mut unitig_graph = Asg::new();
    for arc in &graph.arcs {
        if arc.deleted {
            continue;
        }
        let query_mark = mark[(arc.query_vertex() ^ 1) as usize];
        let target_mark = mark[arc.v as usize];
        if query_mark >= 0 && target_mark >= 0 {
            let u = (query_mark as u32) ^ 1;
            let mut node_length = unitigs[(u >> 1) as usize]
                .length
                .wrapping_sub(arc.overlap) as i32;
            if node_length < 0 {
                node_length = 1;
            }
            unitig_graph.push_arc(AsgArc {
                ul: (u64::from(u) << 32) | u64::from(node_length as u32),
                v: target_mark as u32,
                overlap: arc.overlap,
                deleted: false,
                occurrence: 0,
                strong: false,
                without_large_indel: false,
                exact_region: arc.exact_region,
            });
        }
    }
    for (i, unitig) in unitigs.iter().enumerate() {
        unitig_graph.set_sequence(i, unitig.length, false);
    }
    unitig_graph.cleanup();

    UnitigGraph {
        unitigs,
        graph: unitig_graph,
    }
}
///
/// The scan includes hits whose `del` flag is set, exactly like the reference.
pub fn get_edge_from_source(
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    maximum_hang: i32,
    max_hang_rate: f32,
    minimum_overlap: i32,
    query: u32,
    target: u32,
) -> Option<AsgArc> {
    let query_read = (query >> 1) as usize;
    if query_read >= sources.len() {
        return None;
    }
    for hit in &sources[query_read].hits {
        let query_cut = &coverage_cut[hit.query_identifier as usize];
        let target_cut = &coverage_cut[hit.target_identifier as usize];
        let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
        let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
        let (result, arc) = ma_hit_to_arc(
            hit,
            query_length,
            target_length,
            maximum_hang,
            max_hang_rate,
            minimum_overlap,
        );
        if result < 0 {
            continue;
        }
        if arc.query_vertex() != query {
            continue;
        }
        if arc.v != target {
            continue;
        }
        return Some(arc);
    }
    None
}

/// `get_overlapLen`: the exact and inexact overlap bases stored for one read.
fn get_overlap_length(read: u32, sources: &[PafList]) -> (u32, u32) {
    let mut exact_length: u32 = 0;
    let mut inexact_length: u32 = 0;
    for hit in &sources[read as usize].hits {
        let length = hit.query_end.wrapping_sub(hit.query_start);
        if hit.exact_region {
            exact_length = exact_length.wrapping_add(length);
        } else {
            inexact_length = inexact_length.wrapping_add(length);
        }
    }
    (exact_length, inexact_length)
}

/// `get_specific_edge`: the live `query -> target` arc, from the read graph
/// first, then from the overlap list, then from the polished-edge buffer.
fn get_specific_edge(
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    edges: &[AsgArc],
    read_graph: Option<&Asg>,
    max_hang_rate: f32,
    maximum_hang: i32,
    minimum_overlap: i32,
    query: u32,
    target: u32,
) -> AsgArc {
    let mut result = AsgArc {
        ul: u64::MAX,
        v: u32::MAX,
        ..AsgArc::default()
    };

    if let Some(graph) = read_graph {
        for arc in graph.arcs_of(query) {
            if arc.deleted {
                continue;
            }
            if arc.v == target {
                result = *arc;
                break;
            }
        }
    }

    if result.ul == u64::MAX {
        if let Some(arc) = get_edge_from_source(
            sources,
            coverage_cut,
            maximum_hang,
            max_hang_rate,
            minimum_overlap,
            query,
            target,
        ) {
            result = arc;
        } else {
            result.ul = u64::MAX;
        }
    }

    if result.ul == u64::MAX {
        for arc in edges {
            if arc.deleted {
                continue;
            }
            if arc.query_vertex() == query && arc.v == target {
                result = *arc;
                break;
            }
        }
        // The reference prints `ERROR` when the buffer has no such arc and
        // leaves `ul` at -1; no golden run reaches that branch.
    }

    result
}

/// `detect_exact_ovec`: the first element after `dest_index` whose edge from
/// `src` is an exact overlap.
///
/// The reference takes the read graph but never reads it; the parameter is
/// kept so the call sites read the same as the original.
fn detect_exact_ovec(
    unitig: &MaUtg,
    _read_graph: &Asg,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    edges: &[AsgArc],
    max_hang_rate: f32,
    maximum_hang: i32,
    minimum_overlap: i32,
    src: u32,
    dest_index: usize,
) -> u32 {
    for i in dest_index..unitig.reads.len() {
        let dest = (unitig.reads[i] >> 32) as u32;

        let edge = match get_edge_from_source(
            sources,
            coverage_cut,
            maximum_hang,
            max_hang_rate,
            minimum_overlap,
            src,
            dest,
        ) {
            Some(arc) => Some(arc),
            None => edges.iter().find_map(|arc| {
                if !arc.deleted && arc.query_vertex() == src && arc.v == dest {
                    Some(*arc)
                } else {
                    None
                }
            }),
        };

        let edge = match edge {
            Some(edge) => edge,
            None => return u32::MAX,
        };
        if edge.exact_region != 1 {
            continue;
        }
        return i as u32;
    }

    u32::MAX
}

/// `get_consensus_rate`: how many of the reads before `cur_i` agree with the
/// sequence of the exact overlap that leaves `cur_i`.
fn get_consensus_rate(
    unitig: &MaUtg,
    cur_i: usize,
    next_i: usize,
    read_graph: &Asg,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    edges: &[AsgArc],
    reads: &ReadStore,
    max_hang_rate: f32,
    maximum_hang: i32,
    minimum_overlap: i32,
) -> Option<(i32, i32)> {
    if cur_i < 1 {
        return None;
    }
    if unitig.reads[cur_i] == DELETED_ELEMENT || unitig.reads[next_i] == DELETED_ELEMENT {
        return Some((0, 0));
    }

    let v = (unitig.reads[cur_i] >> 32) as u32;
    let mut w = (unitig.reads[next_i] >> 32) as u32;

    let edge = match get_edge_from_source(
        sources,
        coverage_cut,
        maximum_hang,
        max_hang_rate,
        minimum_overlap,
        v,
        w,
    ) {
        Some(arc) => Some(arc),
        None => edges.iter().find_map(|arc| {
            if !arc.deleted && arc.query_vertex() == v && arc.v == w {
                Some(*arc)
            } else {
                None
            }
        }),
    };
    let edge = edge?;

    let v_begin: u32 = 0;
    let v_end = edge.node_length().wrapping_sub(1);
    let mut recovered_reference: Option<Vec<u8>> = None;
    let mut matches: i32 = 0;
    let mut total: i32 = 0;

    let mut j = cur_i;
    loop {
        if j == 0 {
            break;
        }
        j -= 1;
        if unitig.reads[j] == DELETED_ELEMENT {
            continue;
        }

        w = (unitig.reads[j] >> 32) as u32;
        let edge = match get_edge_from_source(
            sources,
            coverage_cut,
            maximum_hang,
            max_hang_rate,
            minimum_overlap,
            w,
            v,
        ) {
            Some(arc) => Some(arc),
            None => edges.iter().find_map(|arc| {
                if !arc.deleted && arc.query_vertex() == w && arc.v == v {
                    Some(*arc)
                } else {
                    None
                }
            }),
        };
        let edge = match edge {
            Some(edge) => edge,
            None => break,
        };

        let w_begin = edge.node_length();
        let candidate = (w_begin.wrapping_add(v_end)) as i32;
        let limit = read_graph.sequences[(w >> 1) as usize]
            .length
            .wrapping_sub(1) as i32;
        let w_end = candidate.min(limit) as u32;

        // The overlap must cover the same span as the one that left `cur_i`.
        if w_end.wrapping_sub(w_begin) != v_end.wrapping_sub(v_begin) {
            break;
        }

        total += 1;

        if edge.exact_region == 1 {
            matches += 1;
        } else {
            let reference = recovered_reference.get_or_insert_with(|| {
                reads.read_sub_region(
                    (v >> 1) as usize,
                    i64::from(v_begin),
                    i64::from(v_end.wrapping_sub(v_begin).wrapping_add(1)),
                    (v & 1) as u8,
                )
            });
            let query = reads.read_sub_region(
                (w >> 1) as usize,
                i64::from(w_begin),
                i64::from(w_end.wrapping_sub(w_begin).wrapping_add(1)),
                (w & 1) as u8,
            );
            let length = query.len();
            if reference.len() >= length && reference[..length] == query[..] {
                matches += 1;
            }
        }
    }

    Some((matches, total))
}

/// `reduce_ma_utg_t_scaf`: compacts the deleted elements and recomputes every
/// link length from the read graph, the overlap lists or the polished edges.
fn reduce_ma_utg_t_scaf(
    unitig: &mut MaUtg,
    read_graph: &Asg,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    edges: &[AsgArc],
    mut new_edges: Option<&mut Vec<AsgArc>>,
    maximum_hang: i32,
    minimum_overlap: i32,
    max_hang_rate: f32,
) {
    // Drop the elements marked -1 by the polishing passes.
    unitig.reads.retain(|&element| element != DELETED_ELEMENT);
    let element_count = unitig.reads.len();

    let mut total_length: u32 = 0;
    let mut i = 0usize;
    while i + 1 < element_count {
        let v = (unitig.reads[i] >> 32) as u32;
        let w = (unitig.reads[i + 1] >> 32) as u32;
        let link_length = link_length_between(
            v,
            w,
            read_graph,
            sources,
            coverage_cut,
            edges,
            new_edges.as_deref_mut(),
            maximum_hang,
            minimum_overlap,
            max_hang_rate,
        );
        unitig.reads[i] = (u64::from(v) << 32) | u64::from(link_length);
        total_length = total_length.wrapping_add(link_length);
        i += 1;
    }

    if element_count > 0 && i < element_count {
        if unitig.circular {
            let v = (unitig.reads[i] >> 32) as u32;
            let w = (unitig.reads[0] >> 32) as u32;
            let link_length = link_length_between(
                v,
                w,
                read_graph,
                sources,
                coverage_cut,
                edges,
                new_edges.as_deref_mut(),
                maximum_hang,
                minimum_overlap,
                max_hang_rate,
            );
            unitig.reads[i] = (u64::from(v) << 32) | u64::from(link_length);
            total_length = total_length.wrapping_add(link_length);
        } else {
            let v = (unitig.reads[i] >> 32) as u32;
            let sequence_length = read_graph.sequences[(v >> 1) as usize].length;
            unitig.reads[i] = (u64::from(v) << 32) | u64::from(sequence_length);
            total_length = total_length.wrapping_add(sequence_length);
        }
    }

    unitig.length = total_length & 0x7fff_ffff;
    if !unitig.circular {
        unitig.start = (unitig.reads[0] >> 32) as u32;
        unitig.end = ((unitig.reads[element_count - 1] >> 32) as u32) ^ 1;
    }
}

/// The link length between two consecutive unitig elements, the body of the
/// per-pair loop of `reduce_ma_utg_t_scaf`.
#[allow(clippy::too_many_arguments)]
fn link_length_between(
    v: u32,
    w: u32,
    read_graph: &Asg,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    edges: &[AsgArc],
    mut new_edges: Option<&mut Vec<AsgArc>>,
    maximum_hang: i32,
    minimum_overlap: i32,
    max_hang_rate: f32,
) -> u32 {
    // The reference starts from the full read length here, then always
    // overwrites it because `IS_SCAF_READ` is false for corrected HiFi reads.
    for arc in read_graph.arcs_of(v) {
        if arc.deleted {
            continue;
        }
        if arc.v == w {
            return arc.node_length();
        }
    }

    if let Some(index) = edges
        .iter()
        .position(|arc| !arc.deleted && arc.query_vertex() == v && arc.v == w)
    {
        let forward = edges[index];
        let mut length = forward.node_length();
        if let Some(collection) = new_edges.as_deref_mut() {
            collection.push(forward);
            // The reference then searches the reverse complement in the same
            // buffer. When it misses, it pushes one element past the end
            // (undefined); a default arc keeps the behaviour deterministic and
            // the length stays the forward one.
            match edges.iter().find_map(|arc| {
                if !arc.deleted && arc.query_vertex() == (w ^ 1) && arc.v == (v ^ 1) {
                    Some(*arc)
                } else {
                    None
                }
            }) {
                Some(reverse) => {
                    length = reverse.node_length();
                    collection.push(reverse);
                }
                None => collection.push(AsgArc::default()),
            }
        }
        return length;
    }

    // The reference pushes the recovered edges into `newE` only when the
    // caller gave one; the output drivers pass `NULL`.
    match get_edge_from_source(
        sources,
        coverage_cut,
        maximum_hang,
        max_hang_rate,
        minimum_overlap,
        v,
        w,
    ) {
        Some(forward) => {
            let length = forward.node_length();
            if let Some(collection) = new_edges.as_deref_mut() {
                collection.push(forward);
                // `t_f` keeps the forward arc when the reverse lookup fails.
                let reverse = get_edge_from_source(
                    sources,
                    coverage_cut,
                    maximum_hang,
                    max_hang_rate,
                    minimum_overlap,
                    w ^ 1,
                    v ^ 1,
                )
                .unwrap_or(forward);
                collection.push(reverse);
            }
            length
        }
        None => {
            // The reference reads the length of an uninitialised `t_f` here and
            // prints `####ERROR1`. No golden run reaches it.
            eprintln!(
                "####ERROR1: v>>1: {}, v&1: {}, w>>1: {}, w&1: {}",
                v >> 1,
                v & 1,
                w >> 1,
                w & 1
            );
            u32::MAX
        }
    }
}

/// `polish_unitig_scaf`: removes elements whose links are poorly supported.
pub(crate) fn polish_unitig_scaf(
    unitig: &mut MaUtg,
    read_graph: &Asg,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    edges: &[AsgArc],
    new_edges: Option<&mut Vec<AsgArc>>,
    maximum_hang: i32,
    minimum_overlap: i32,
    max_hang_rate: f32,
) -> u32 {
    if unitig.reads.len() < 3 {
        return 0;
    }

    let element_count = unitig.reads.len();
    let mut skip: u32 = 0;

    // The reference loop is `for (z = 1, l = 0; z <= in->n; z++)` with the
    // body guarded by `z == in->n || IS_SCAF_READ(R_INF, (in->a[z])>>33)`.
    // `IS_SCAF_READ` tests `read_sperate[id] == NULL`, and every read of the
    // corrected corpus owns a packed sequence, so the guard is false until
    // `z == in->n` and the body runs exactly once, with `l == 0`.
    {
        let z = element_count;
        let window_start = 0usize;
        if z - window_start >= 3 {
            let mut pre = (unitig.reads[window_start] >> 32) as u32;
            let mut pre_i = window_start;
            for i in (window_start + 1)..(z - 1) {
                if unitig.reads[i] == DELETED_ELEMENT {
                    continue;
                }
                let v = (unitig.reads[i] >> 32) as u32;
                let mut afte = (unitig.reads[i + 1] >> 32) as u32;

                let pre_edge = get_specific_edge(
                    sources,
                    coverage_cut,
                    edges,
                    if pre_i == i - 1 { Some(read_graph) } else { None },
                    max_hang_rate,
                    maximum_hang,
                    minimum_overlap,
                    v ^ 1,
                    pre ^ 1,
                );
                let after_edge = get_specific_edge(
                    sources,
                    coverage_cut,
                    edges,
                    Some(read_graph),
                    max_hang_rate,
                    maximum_hang,
                    minimum_overlap,
                    v,
                    afte,
                );

                if pre_edge.exact_region == 1 && after_edge.exact_region == 1 {
                    pre = v;
                    pre_i = i;
                    continue;
                }

                let afte_i = detect_exact_ovec(
                    unitig,
                    read_graph,
                    sources,
                    coverage_cut,
                    edges,
                    max_hang_rate,
                    maximum_hang,
                    minimum_overlap,
                    pre,
                    i + 1,
                );
                if afte_i == u32::MAX {
                    pre = v;
                    pre_i = i;
                    continue;
                }
                afte = (unitig.reads[afte_i as usize] >> 32) as u32;

                let mut minimum_inexact: u32 = u32::MAX;
                let mut maximum_exact: u32 = 0;
                for k in i..(afte_i as usize) {
                    let (exact_length, inexact_length) =
                        get_overlap_length((unitig.reads[k] >> 33) as u32, sources);
                    if inexact_length < minimum_inexact {
                        minimum_inexact = inexact_length;
                        maximum_exact = exact_length;
                    }
                }

                let (exact_length, inexact_length) = get_overlap_length(pre >> 1, sources);
                if inexact_length > minimum_inexact
                    || (inexact_length == minimum_inexact && exact_length <= maximum_exact)
                {
                    pre = v;
                    pre_i = i;
                    continue;
                }

                let (exact_length, inexact_length) = get_overlap_length(afte >> 1, sources);
                if inexact_length > minimum_inexact
                    || (inexact_length == minimum_inexact && exact_length <= maximum_exact)
                {
                    pre = v;
                    pre_i = i;
                    continue;
                }

                for k in i..(afte_i as usize) {
                    unitig.reads[k] = DELETED_ELEMENT;
                    skip += 1;
                }
            }
        }
    }

    if skip == 0 {
        return 0;
    }
    reduce_ma_utg_t_scaf(
        unitig,
        read_graph,
        sources,
        coverage_cut,
        edges,
        new_edges,
        maximum_hang,
        minimum_overlap,
        max_hang_rate,
    );
    1
}

/// `polish_unitig_advance_scaf`: drops the prefix elements that the following
/// reads do not support.
pub(crate) fn polish_unitig_advance_scaf(
    unitig: &mut MaUtg,
    read_graph: &Asg,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    edges: &[AsgArc],
    new_edges: Option<&mut Vec<AsgArc>>,
    reads: &ReadStore,
    maximum_hang: i32,
    minimum_overlap: i32,
    max_hang_rate: f32,
) -> u32 {
    if unitig.reads.len() < 3 {
        return 0;
    }

    let element_count = unitig.reads.len();
    let mut skip: u32 = 0;

    // Same guard as `polish_unitig_scaf`: `IS_SCAF_READ` is false for every
    // read of the corrected corpus, so the body runs exactly once, at
    // `z == in_n` with `l == 0`.
    {
        let z = element_count;
        let window_start = 0usize;
        if z - window_start >= 3 {
            let mut i = window_start + 1;
            while i < z - 1 {
                let (mut match_count, mut total_count) = match get_consensus_rate(
                    unitig,
                    i,
                    i + 1,
                    read_graph,
                    sources,
                    coverage_cut,
                    edges,
                    reads,
                    max_hang_rate,
                    maximum_hang,
                    minimum_overlap,
                ) {
                    None => {
                        i += 1;
                        continue;
                    }
                    Some(counts) => counts,
                };
                if total_count == 0 {
                    // `match_rate` is zero, which never clears the support test.
                }

                let mut match_rate = if total_count == 0 {
                    0.0f64
                } else {
                    f64::from(match_count) / f64::from(total_count)
                };
                // Most reads support element `i`, so it is right.
                if f64::from(match_count) >= f64::from(total_count) * 0.5
                    && total_count > 0
                    && match_count > 0
                {
                    i += 1;
                    continue;
                }

                let mut max_i = i;
                let mut match_max = match_count;
                let mut match_rate_max = match_rate;

                let mut k = i;
                while k > 0 {
                    k -= 1;
                    if unitig.reads[k] == DELETED_ELEMENT {
                        continue;
                    }
                    let (previous_match, previous_total) = match get_consensus_rate(
                        unitig,
                        k,
                        i + 1,
                        read_graph,
                        sources,
                        coverage_cut,
                        edges,
                        reads,
                        max_hang_rate,
                        maximum_hang,
                        minimum_overlap,
                    ) {
                        None => break,
                        Some(counts) => counts,
                    };
                    match_count = previous_match;
                    total_count = previous_total;

                    if total_count == 0 {
                        break;
                    }
                    match_rate = f64::from(match_count) / f64::from(total_count);
                    if match_rate > match_rate_max
                        || (match_rate == match_rate_max && match_count > match_max)
                    {
                        max_i = k;
                        match_max = match_count;
                        match_rate_max = match_rate;
                    }
                }

                for k in (max_i + 1)..=i {
                    if unitig.reads[k] == DELETED_ELEMENT {
                        continue;
                    }
                    unitig.reads[k] = DELETED_ELEMENT;
                    skip += 1;
                }

                i += 1;
            }
        }
    }

    if skip == 0 {
        return 1;
    }
    reduce_ma_utg_t_scaf(
        unitig,
        read_graph,
        sources,
        coverage_cut,
        edges,
        new_edges,
        maximum_hang,
        minimum_overlap,
        max_hang_rate,
    );
    1
}

/// `ma_ug_seq`: polishes every unitig and materialises its sequence.
///
/// `edges` is the reference's `edge` (the caller's `new_rtg_edges`), which the
/// polishing passes only read. `new_edges` is the reference's `E`: when it is
/// given, the arcs recovered from the overlap lists during polishing are
/// collected there and appended to the read graph at the end. The output
/// drivers pass `NULL` for `E`, so their read graphs stay unchanged.
#[allow(clippy::too_many_arguments)]
pub fn ma_ug_seq(
    unitig_graph: &mut UnitigGraph,
    read_graph: &mut Asg,
    coverage_cut: &[CoverageCut],
    sources: &[PafList],
    edges: &[AsgArc],
    new_edges: Option<&mut Vec<AsgArc>>,
    maximum_hang: i32,
    minimum_overlap: i32,
    is_polish: bool,
    reads: &ReadStore,
    max_hang_rate: f32,
) {
    let mut new_edges = new_edges;
    let UnitigGraph { unitigs, graph } = unitig_graph;

    for (i, unitig) in unitigs.iter_mut().enumerate() {
        if is_polish {
            polish_unitig_scaf(
                unitig,
                read_graph,
                sources,
                coverage_cut,
                edges,
                new_edges.as_deref_mut(),
                maximum_hang,
                minimum_overlap,
                max_hang_rate,
            );
            polish_unitig_advance_scaf(
                unitig,
                read_graph,
                sources,
                coverage_cut,
                edges,
                new_edges.as_deref_mut(),
                reads,
                maximum_hang,
                minimum_overlap,
                max_hang_rate,
            );
        }

        graph.sequences[i].length = unitig.length;

        let mut sequence = vec![b'N'; unitig.length as usize];
        let mut l: u32 = 0;
        for &element in &unitig.reads {
            let read_id = (element >> 33) as u32;
            let orientation = ((element >> 32) & 1) as u8;
            let start = l;
            let element_length = element as u32;
            l = l.wrapping_add(element_length);

            if element_length == 0 {
                continue;
            }

            let read_sequence = recover_read(reads, read_id);
            let cut = &coverage_cut[read_id as usize];
            let read_slice = &read_sequence[cut.start as usize..cut.end as usize];
            let read_length = cut.end.wrapping_sub(cut.start) as usize;

            if orientation == 0 {
                for k in 0..element_length as usize {
                    sequence[start as usize + k] = read_slice[k];
                }
            } else {
                for k in 0..element_length as usize {
                    sequence[start as usize + k] = complement_byte(read_slice[read_length - 1 - k]);
                }
            }
        }
        unitig.sequence = sequence;
    }

    // The unitig node length is the drawn length minus the overlap, and the
    // low half of every inter-unitig arc becomes that node length.
    let vertex_count = graph.sequences.len() * 2;
    for v in 0..vertex_count as u32 {
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        let (start, count) = graph.arc_range(v);
        for index in start..start + count {
            let arc = &mut graph.arcs[index];
            if arc.deleted {
                continue;
            }
            let target_unitig = (arc.ul >> 33) as usize;
            let node_length = graph.sequences[target_unitig]
                .length
                .wrapping_sub(arc.overlap);
            arc.ul = ((arc.ul >> 32) << 32) | u64::from(node_length);
        }
    }

    if let Some(new_edges) = new_edges {
        if !new_edges.is_empty() {
            for arc in new_edges.iter() {
                read_graph.push_arc(*arc);
            }
            read_graph.index.clear();
            read_graph.is_sorted = false;
            read_graph.cleanup();
        }
    }
}

/// `recover_UC_Read` for a stored read. The reference falls back to
/// `recover_fake_read` for the synthetic reads of `asg_t::F_seq`; no stage of
/// the HiFi pipeline creates those.
fn recover_read(reads: &ReadStore, read_id: u32) -> Vec<u8> {
    assert!(
        u64::from(read_id) < reads.total_reads,
        "fake unitig reads are not supported yet (read {})",
        read_id
    );
    reads.read_sequence(read_id as usize)
}

/// `get_ug_coverage`: the read coverage of one unitig, as a percentage.
fn get_ug_coverage(
    unitig: &MaUtg,
    read_graph: &Asg,
    coverage_cut: &[CoverageCut],
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    read_flags: &mut [u8],
) -> u32 {
    let (read_bases, covered_bases) = get_ug_coverage_bases(
        unitig,
        read_graph,
        coverage_cut,
        sources,
        read_to_unitig,
        read_flags,
    );
    if read_bases == 0 {
        0
    } else {
        (covered_bases / read_bases) as u32
    }
}

/// The `rR`/`rC` out-parameters of `get_ug_coverage`
/// (`Overlaps.cpp` 10753): total read bases and covered bases.
pub(crate) fn get_ug_coverage_bases(
    unitig: &MaUtg,
    read_graph: &Asg,
    coverage_cut: &[CoverageCut],
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    read_flags: &mut [u8],
) -> (u64, u64) {
    let mut read_bases: i64 = 0;
    let mut covered_bases: i64 = 0;

    for &element in &unitig.reads {
        if element == DELETED_ELEMENT {
            continue;
        }
        read_flags[(element >> 33) as usize] = 1;
    }

    for &element in &unitig.reads {
        if element == DELETED_ELEMENT {
            continue;
        }
        let read_id = (element >> 33) as u32;
        let cut = &coverage_cut[read_id as usize];
        read_bases += cut.end.wrapping_sub(cut.start) as i64;

        for hit in &sources[read_id as usize].hits {
            if !hit.exact_region {
                continue;
            }
            let mut target = hit.target_identifier;
            if read_graph.sequences[target as usize].deleted {
                let (mapped, is_unitig) = read_to_unitig.get(target);
                target = mapped;
                if target == u32::MAX
                    || is_unitig == 1
                    || read_graph.sequences[target as usize].deleted
                {
                    continue;
                }
            }
            if read_graph.sequences[target as usize].deleted {
                continue;
            }
            if read_flags[target as usize] != 1 {
                continue;
            }
            covered_bases += hit.query_end.wrapping_sub(hit.query_start) as i64;
        }
    }

    for &element in &unitig.reads {
        if element == DELETED_ELEMENT {
            continue;
        }
        read_flags[(element >> 33) as usize] = 0;
    }

    (read_bases as u64, covered_bases as u64)
}

/// `cal_circle_ov`: the overlap that closes a circular unitig.
fn cal_circle_ov(unitig_graph: &UnitigGraph, unitig_id: u32, reverse: u32, graph: &Asg) -> u32 {
    let v = (unitig_id << 1) + u32::from(reverse != 0);
    let unitig = &unitig_graph.unitigs[unitig_id as usize];
    let first = unitig.reads[0];
    let last = unitig.reads[unitig.reads.len() - 1];

    let vx = if v & 1 != 0 {
        ((first >> 32) as u32) ^ 1
    } else {
        (last >> 32) as u32
    };
    let wx = if v & 1 != 0 {
        ((last >> 32) as u32) ^ 1
    } else {
        (first >> 32) as u32
    };

    for arc in graph.arcs_of(vx) {
        if arc.v == wx {
            return arc.overlap;
        }
    }
    0
}

/// `ma_ug_print2`: the GFA `S`/`A`/`L` lines of one unitig graph.
#[allow(clippy::too_many_arguments)]
pub fn ma_ug_print2(
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
    coverage_cut: Option<&[CoverageCut]>,
    sources: Option<&[PafList]>,
    read_to_unitig: Option<&ReadToUnitig>,
    print_sequence: bool,
    prefix: &str,
    output: &mut dyn Write,
    reads: &ReadStore,
    trio_labels: &[u8],
) -> io::Result<()> {
    let can_compute_coverage =
        coverage_cut.is_some() && sources.is_some() && read_to_unitig.is_some();
    let mut read_flags = vec![0u8; read_graph.sequences.len()];

    for (i, unitig) in unitig_graph.unitigs.iter().enumerate() {
        let name = format!(
            "{}{:06}{}",
            prefix,
            i + 1,
            if unitig.circular { 'c' } else { 'l' }
        );

        let coverage = if can_compute_coverage {
            get_ug_coverage(
                unitig,
                read_graph,
                coverage_cut.unwrap(),
                sources.unwrap(),
                read_to_unitig.unwrap(),
                &mut read_flags,
            )
        } else {
            0
        };

        if print_sequence {
            output.write_all(b"S\t")?;
            output.write_all(name.as_bytes())?;
            output.write_all(b"\t")?;
            output.write_all(&unitig.sequence)?;
            write!(
                output,
                "\tLN:i:{}\trd:i:{}\n",
                unitig.length, coverage
            )?;
        } else {
            write!(
                output,
                "S\t{}\t*\tLN:i:{}\trd:i:{}\n",
                name, unitig.length, coverage
            )?;
        }

        let mut l: u32 = 0;
        for &element in &unitig.reads {
            if element != DELETED_ELEMENT {
                let x = (element >> 33) as u32;
                let strand = if (element >> 32) & 1 == 0 { '+' } else { '-' };
                if u64::from(x) < reads.total_reads {
                    let cut = coverage_cut.unwrap()[x as usize];
                    let label = trio_labels[x as usize];
                    let hg_index = if label != FATHER && label != MOTHER {
                        AMBIGU
                    } else {
                        label
                    };
                    let hg = b"apmaaa"[hg_index as usize] as char;
                    write!(
                        output,
                        "A\t{}\t{}\t{}\t{}\t{}\t{}\tid:i:{}\tHG:A:{}\n",
                        name,
                        l as i32,
                        strand,
                        reads.read_name(x as usize),
                        cut.start,
                        cut.end,
                        x,
                        hg
                    )?;
                } else {
                    let cut = coverage_cut.unwrap()[x as usize];
                    write!(
                        output,
                        "A\t{}\t{}\t{}\tFAKE\t{}\t{}\tid:i:{}\tHG:A:*\n",
                        name, l as i32, strand, cut.start, cut.end, x
                    )?;
                }
            } else {
                write!(output, "A\t{}\t{}\t*\t*\t*\t*\tid:i:*\tHG:A:*\n", name, l as i32)?;
            }
            l = l.wrapping_add(element as u32);
        }
    }

    for (i, unitig) in unitig_graph.unitigs.iter().enumerate() {
        if unitig.circular {
            for reverse in 0..2u32 {
                let overlap = cal_circle_ov(unitig_graph, i as u32, reverse, read_graph);
                let l1 = if unitig.length >= overlap {
                    unitig.length - overlap
                } else {
                    0
                };
                let strand = if reverse == 0 { '+' } else { '-' };
                write!(
                    output,
                    "L\t{}{:06}c\t{}\t{}{:06}c\t{}\t{}M\tL1:i:{}\tL2:i:0\n",
                    prefix,
                    i + 1,
                    strand,
                    prefix,
                    i + 1,
                    strand,
                    overlap as i32,
                    l1
                )?;
            }
        } else {
            for u in [i as u32 * 2, i as u32 * 2 + 1] {
                let (start, count) = unitig_graph.graph.arc_range(u);
                for index in start..start + count {
                    let arc = unitig_graph.graph.arcs[index];
                    if arc.deleted {
                        continue;
                    }
                    let v = arc.v;
                    let target_unitig = &unitig_graph.unitigs[(v >> 1) as usize];
                    write!(
                        output,
                        "L\t{}{:06}{}\t{}\t{}{:06}{}\t{}\t{}M\tL1:i:{}\tL2:i:0\n",
                        prefix,
                        i + 1,
                        if unitig.circular { 'c' } else { 'l' },
                        if u & 1 == 0 { '+' } else { '-' },
                        prefix,
                        (v >> 1) + 1,
                        if target_unitig.circular { 'c' } else { 'l' },
                        if v & 1 == 0 { '+' } else { '-' },
                        arc.overlap as i32,
                        arc.node_length()
                    )?;
                }
            }
        }
    }

    Ok(())
}

/// `ma_ug_print`: [`ma_ug_print2`] with sequences.
#[allow(clippy::too_many_arguments)]
pub fn ma_ug_print(
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
    coverage_cut: Option<&[CoverageCut]>,
    sources: Option<&[PafList]>,
    read_to_unitig: Option<&ReadToUnitig>,
    prefix: &str,
    output: &mut dyn Write,
    reads: &ReadStore,
    trio_labels: &[u8],
) -> io::Result<()> {
    ma_ug_print2(
        unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        read_to_unitig,
        true,
        prefix,
        output,
        reads,
        trio_labels,
    )
}

/// `ma_ug_print_simple`: [`ma_ug_print2`] without sequences.
#[allow(clippy::too_many_arguments)]
pub fn ma_ug_print_simple(
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
    coverage_cut: Option<&[CoverageCut]>,
    sources: Option<&[PafList]>,
    read_to_unitig: Option<&ReadToUnitig>,
    prefix: &str,
    output: &mut dyn Write,
    reads: &ReadStore,
    trio_labels: &[u8],
) -> io::Result<()> {
    ma_ug_print2(
        unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        read_to_unitig,
        false,
        prefix,
        output,
        reads,
        trio_labels,
    )
}

/// `inter_interval`: the intersection of two inclusive intervals.
fn inter_interval(a_start: i32, a_end: i32, b_start: i32, b_end: i32) -> Option<(i32, i32)> {
    if a_start > b_end || b_start > a_end {
        return None;
    }
    Some((a_start.max(b_start), a_end.min(b_end)))
}

/// `if_exact_match`: whether the two sub-ranges hold the same bases.
#[allow(clippy::too_many_arguments)]
fn if_exact_match(
    x: &[u8],
    _x_length: usize,
    y: &[u8],
    _y_length: usize,
    x_begin: i32,
    x_end: i32,
    y_begin: i32,
    y_end: i32,
) -> bool {
    let overlap_length = x_end - x_begin + 1;
    if y_end - y_begin + 1 != overlap_length {
        return false;
    }
    let x_start = x_begin as usize;
    let y_start = y_begin as usize;
    let length = overlap_length as usize;
    x[x_start..x_start + length] == y[y_start..y_start + length]
}

/// Writes one low-quality interval of [`print_rough_inconsistent_sites`].
#[allow(clippy::too_many_arguments)]
fn write_bed_interval(
    output: &mut dyn Write,
    prefix: &str,
    unitig_id: usize,
    circular: bool,
    v: u32,
    v_begin: i32,
    v_end: i32,
    mut interval_begin: i32,
    mut interval_end: i32,
    total: u32,
    inexact: u32,
    coordinate_begin: u32,
    hits: &[PafHit],
    source_reads: &ReadStore,
) -> io::Result<()> {
    // The hit scan uses the interval in the read's forward coordinates, which
    // is what the reference keeps before it flips the printed coordinates.
    let forward_begin = interval_begin + v_begin;
    let forward_end = interval_end + v_begin;

    if v & 1 == 1 {
        let span = v_end - v_begin + 1;
        interval_begin = span - interval_begin - 1;
        interval_end = span - interval_end - 1;
        std::mem::swap(&mut interval_begin, &mut interval_end);
    }
    interval_end += 1;
    let rate = if total == 0 {
        100
    } else {
        inexact.wrapping_mul(100) / total
    };

    write!(
        output,
        "{}{:06}{}\t{}\t{}\t{}\t{}",
        prefix,
        unitig_id,
        if circular { 'c' } else { 'l' },
        (interval_begin as u32).wrapping_add(coordinate_begin),
        (interval_end as u32).wrapping_add(coordinate_begin),
        rate,
        source_reads.read_name((v >> 1) as usize),
    )?;
    for hit in hits {
        let intersects = inter_interval(
            forward_begin,
            forward_end,
            hit.query_start as i32,
            hit.query_end.wrapping_sub(1) as i32,
        );
        if intersects.is_none() {
            continue;
        }
        write!(
            output,
            ",{}",
            source_reads.read_name(hit.target_identifier as usize)
        )?;
    }
    writeln!(output)?;
    Ok(())
}

/// `print_rough_inconsistent_sites`: finds the parts of one unitig element
/// where the reads that cover it disagree and writes them as BED records.
#[allow(clippy::too_many_arguments)]
fn print_rough_inconsistent_sites(
    unitig: &MaUtg,
    cur_i: usize,
    next_i: usize,
    read_graph: &Asg,
    source_reads: &ReadStore,
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    edges: &[AsgArc],
    maximum_hang: i32,
    minimum_overlap: i32,
    coordinate_begin: u32,
    rate_threshold: u32,
    prefix: &str,
    unitig_id: usize,
    output: &mut dyn Write,
    max_hang_rate: f32,
) -> io::Result<()> {
    let element_count = unitig.reads.len();
    // The last element of a linear unitig has no successor.
    let next = if cur_i == element_count - 1 && next_i == element_count {
        if unitig.circular {
            0
        } else {
            usize::MAX
        }
    } else {
        next_i
    };

    let v = (unitig.reads[cur_i] >> 32) as u32;
    let read_id = (v >> 1) as usize;
    let read_length = source_reads.read_lengths[read_id] as i32;

    let (v_begin, v_end) = if next != usize::MAX {
        let w = (unitig.reads[next] >> 32) as u32;
        let edge = get_specific_edge(
            sources,
            coverage_cut,
            edges,
            Some(read_graph),
            max_hang_rate,
            maximum_hang,
            minimum_overlap,
            v,
            w,
        );
        let mut v_begin: i32 = 0;
        let mut v_end = edge.node_length().wrapping_sub(1) as i32;
        if v & 1 == 1 {
            v_begin = read_length - v_begin - 1;
            v_end = read_length - v_end - 1;
            std::mem::swap(&mut v_begin, &mut v_end);
        }
        (v_begin, v_end)
    } else {
        (0, read_length - 1)
    };

    let position_count = (v_end - v_begin + 1) as usize;
    let mut exact_count = vec![0u32; position_count];
    let mut total_count = vec![0u32; position_count];

    let reference = source_reads.read_sub_region(
        read_id,
        i64::from(v_begin),
        i64::from(v_end - v_begin + 1),
        0,
    );

    let hits = &sources[read_id].hits;
    for hit in hits {
        let Some((v_sub_begin, v_sub_end)) = inter_interval(
            v_begin,
            v_end,
            hit.query_start as i32,
            hit.query_end.wrapping_sub(1) as i32,
        ) else {
            continue;
        };

        if hit.exact_region {
            for j in v_sub_begin..=v_sub_end {
                exact_count[(j - v_begin) as usize] += 1;
                total_count[(j - v_begin) as usize] += 1;
            }
            continue;
        }

        let target_id = hit.target_identifier as usize;
        let target_length = source_reads.read_lengths[target_id] as i32;
        let mut w_begin = hit.target_start as i32;
        let mut w_end = hit.target_end.wrapping_sub(1) as i32;
        if hit.reverse {
            w_begin = target_length - w_begin - 1;
            w_end = target_length - w_end - 1;
            std::mem::swap(&mut w_begin, &mut w_end);
        }

        let mut w_sub_begin = w_begin + (v_sub_begin - hit.query_start as i32);
        if w_sub_begin >= target_length {
            w_sub_begin = target_length - 1;
        }
        let mut w_sub_end = w_end - (hit.query_end as i32 - 1 - v_sub_end);
        if w_sub_end < 0 {
            w_sub_end = 0;
        }

        if w_sub_begin > w_sub_end || (v_sub_end - v_sub_begin) != (w_sub_end - w_sub_begin) {
            for j in v_sub_begin..=v_sub_end {
                total_count[(j - v_begin) as usize] += 1;
            }
            continue;
        }

        let query = source_reads.read_sub_region(
            target_id,
            i64::from(w_sub_begin),
            i64::from(w_sub_end - w_sub_begin + 1),
            u8::from(hit.reverse),
        );
        if if_exact_match(
            &reference,
            reference.len(),
            &query,
            query.len(),
            v_sub_begin - v_begin,
            v_sub_end - v_begin,
            0,
            query.len() as i32 - 1,
        ) {
            for j in v_sub_begin..=v_sub_end {
                exact_count[(j - v_begin) as usize] += 1;
                total_count[(j - v_begin) as usize] += 1;
            }
        } else {
            for j in v_sub_begin..=v_sub_end {
                total_count[(j - v_begin) as usize] += 1;
            }
        }
    }

    let mut interval_begin: i32 = -1;
    let mut interval_end: i32 = -1;
    let mut inexact: u32 = 0;
    let mut total: u32 = 0;

    for i in 0..position_count {
        let rate = if total_count[i] == 0 {
            100
        } else {
            total_count[i].wrapping_sub(exact_count[i]).wrapping_mul(100) / total_count[i]
        };
        if rate >= rate_threshold {
            if interval_begin == -1 && interval_end == -1 {
                interval_begin = i as i32;
                interval_end = i as i32;
            } else {
                interval_end += 1;
            }
            total = total.wrapping_add(total_count[i]);
            inexact = inexact.wrapping_add(total_count[i].wrapping_sub(exact_count[i]));
        } else {
            if interval_begin != -1 && interval_end != -1 {
                write_bed_interval(
                    output,
                    prefix,
                    unitig_id,
                    unitig.circular,
                    v,
                    v_begin,
                    v_end,
                    interval_begin,
                    interval_end,
                    total,
                    inexact,
                    coordinate_begin,
                    hits,
                    source_reads,
                )?;
            }
            interval_begin = -1;
            interval_end = -1;
            inexact = 0;
            total = 0;
        }
    }

    if interval_begin != -1 && interval_end != -1 {
        write_bed_interval(
            output,
            prefix,
            unitig_id,
            unitig.circular,
            v,
            v_begin,
            v_end,
            interval_begin,
            interval_end,
            total,
            inexact,
            coordinate_begin,
            hits,
            source_reads,
        )?;
    }

    Ok(())
}

/// `ma_ug_print_bed`: the inconsistent regions of every unitig, as BED.
#[allow(clippy::too_many_arguments)]
pub fn ma_ug_print_bed(
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
    coverage_cut: &[CoverageCut],
    sources: &[PafList],
    edges: &[AsgArc],
    maximum_hang: i32,
    minimum_overlap: i32,
    rate_threshold: u32,
    prefix: &str,
    output: &mut dyn Write,
    reads: &ReadStore,
    max_hang_rate: f32,
) -> io::Result<()> {
    for (i, unitig) in unitig_graph.unitigs.iter().enumerate() {
        if unitig.reads.len() < 2 {
            continue;
        }
        let mut l: u32 = 0;
        for j in 0..unitig.reads.len() {
            let start = l;
            let element_length = unitig.reads[j] as u32;
            l = l.wrapping_add(element_length);
            print_rough_inconsistent_sites(
                unitig,
                j,
                j + 1,
                read_graph,
                reads,
                sources,
                coverage_cut,
                edges,
                maximum_hang,
                minimum_overlap,
                start,
                rate_threshold,
                prefix,
                i + 1,
                output,
                max_hang_rate,
            )?;
        }
    }
    Ok(())
}

/// `gen_spec_edge`: the stored overlap that becomes exactly `v -> w`.
pub(crate) fn gen_spec_edge(
    read_graph: &Asg,
    sources: &[PafList],
    maximum_hang: i32,
    max_hang_rate: f32,
    minimum_overlap: i32,
    v: u32,
    w: u32,
) -> Option<AsgArc> {
    for hit in &sources[(v >> 1) as usize].hits {
        let target = hit.target_identifier;
        if target != (w >> 1) {
            continue;
        }
        let query = hit.query_identifier;
        let (result, mut arc) = ma_hit_to_arc(
            hit,
            read_graph.sequences[query as usize].length as i32,
            read_graph.sequences[target as usize].length as i32,
            maximum_hang,
            max_hang_rate,
            minimum_overlap,
        );
        if result < 0 {
            continue;
        }
        if arc.query_vertex() != v || arc.v != w {
            continue;
        }
        arc.occurrence = 0;
        return Some(arc);
    }
    None
}

/// `gen_ext_tip`: one step of a tip extension.
///
/// Returns the reference's `rr` (`1` extended, `0` nothing to extend, `-1`
/// rejected), the target vertex and the updated path length. The reference
/// only resets `pmax.ol` on the second scan, so a stale `pmax.v` decides the
/// `(pmax.v ^ 1) == v0` test; that behaviour is preserved here.
#[allow(clippy::too_many_arguments)]
fn gen_ext_tip(
    mut v: u32,
    read_graph: &Asg,
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    final_flag: &mut [u8],
    min_overlap: u32,
    maximum_hang: i32,
    max_hang_rate: f32,
    minimum_overlap: i32,
    path_length: u64,
) -> (i32, u64, u64) {
    let mut best = AsgArc {
        overlap: 0,
        ..AsgArc::default()
    };

    for hit in &sources[(v >> 1) as usize].hits {
        let query = hit.query_identifier;
        let target = hit.target_identifier;
        let (result, arc) = ma_hit_to_arc(
            hit,
            read_graph.sequences[query as usize].length as i32,
            read_graph.sequences[target as usize].length as i32,
            maximum_hang,
            max_hang_rate,
            minimum_overlap,
        );
        if result < 0 {
            continue;
        }
        if arc.query_vertex() != v {
            continue;
        }
        if arc.overlap < min_overlap {
            continue;
        }
        if arc.overlap > best.overlap {
            best = arc;
        }
    }
    if best.overlap == 0 {
        return (0, 0, path_length);
    }

    let target_read = best.v >> 1;
    if !read_graph.sequences[target_read as usize].deleted || final_flag[target_read as usize] != 0
    {
        return (-1, 0, path_length);
    }
    let (contained, is_unitig) = read_to_unitig.get(target_read);
    if is_unitig == 1 {
        return (-1, 0, path_length);
    }
    if contained != u32::MAX
        && (!read_graph.sequences[contained as usize].deleted || final_flag[contained as usize] != 0)
    {
        return (-1, 0, path_length);
    }

    let mut path_length = path_length;
    if path_length >= u64::from(read_graph.sequences[(v >> 1) as usize].length) {
        path_length -= u64::from(read_graph.sequences[(v >> 1) as usize].length);
    } else {
        path_length = 0;
    }
    path_length += u64::from(best.node_length())
        + u64::from(read_graph.sequences[(best.v >> 1) as usize].length);

    let v0 = v;
    let result_vertex = u64::from(best.v);
    v = best.v ^ 1;
    best.overlap = 0;

    for hit in &sources[(v >> 1) as usize].hits {
        let query = hit.query_identifier;
        let target = hit.target_identifier;
        let (result, arc) = ma_hit_to_arc(
            hit,
            read_graph.sequences[query as usize].length as i32,
            read_graph.sequences[target as usize].length as i32,
            maximum_hang,
            max_hang_rate,
            minimum_overlap,
        );
        if result < 0 {
            continue;
        }
        if arc.query_vertex() != v {
            continue;
        }
        if arc.overlap < min_overlap {
            continue;
        }
        if arc.overlap > best.overlap {
            best = arc;
        }
    }

    if (best.v ^ 1) == v0 {
        final_flag[(result_vertex >> 1) as usize] = 1;
        return (1, result_vertex, path_length);
    }
    (-1, u64::MAX, path_length)
}

/// `ug_ext_gfa`: extends tips by following read overlaps backwards, then puts
/// the extension chains back into the graph.
#[allow(clippy::too_many_arguments)]
pub fn ug_ext_gfa(
    graph: &mut Asg,
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    tips_length: u32,
    maximum_hang: i32,
    minimum_overlap: i32,
    max_hang_rate: f32,
    maximum_length: u32,
) -> u32 {
    let vertex_count = graph.sequences.len() * 2;
    let tip_count = tips_length + 1;
    let mut final_flag = vec![0u8; graph.sequences.len()];
    let mut extension: Vec<u64> = Vec::new();
    let mut groups: Vec<u64> = Vec::new();

    for v in 0..vertex_count as u32 {
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if graph.arcs_of(v).iter().any(|arc| !arc.deleted) {
            continue;
        }

        let mut v0 = u64::from(v);
        let mut rr: i32 = 0;
        let group_start = extension.len();
        let mut path_length = u64::from(graph.sequences[(v >> 1) as usize].length);
        let mut z: u32 = 0;
        while z < tip_count || path_length < u64::from(maximum_length) {
            let (result, target, updated_length) = gen_ext_tip(
                v0 as u32,
                graph,
                sources,
                read_to_unitig,
                &mut final_flag,
                2000,
                maximum_hang,
                max_hang_rate,
                minimum_overlap,
                path_length,
            );
            rr = result;
            path_length = updated_length;
            if rr <= 0 {
                break;
            }
            extension.push((v0 << 32) | target);
            v0 = target;
            z += 1;
        }

        if rr < 0 || (z >= tip_count && path_length >= u64::from(maximum_length)) {
            for &packed in &extension[group_start..] {
                final_flag[(packed as u32 >> 1) as usize] = 0;
            }
            extension.truncate(group_start);
        }
        if extension.len() > group_start {
            groups.push(((group_start as u64) << 32) | extension.len() as u64);
        }
    }

    if !groups.is_empty() {
        for &group in &groups {
            let group_start = (group >> 32) as usize;
            let group_end = group as u32 as usize;
            for index in group_start..group_end {
                let v = (extension[index] >> 32) as u32;
                let w = extension[index] as u32;
                assert!(
                    !graph.sequences[(v >> 1) as usize].deleted
                        && graph.sequences[(w >> 1) as usize].deleted
                );
                graph.sequences[(w >> 1) as usize].deleted = false;

                let forward = gen_spec_edge(
                    graph,
                    sources,
                    maximum_hang,
                    max_hang_rate,
                    minimum_overlap,
                    v,
                    w,
                )
                .expect("ug_ext_gfa: the forward overlap must exist");
                graph.push_arc(forward);

                let reverse = gen_spec_edge(
                    graph,
                    sources,
                    maximum_hang,
                    max_hang_rate,
                    minimum_overlap,
                    w ^ 1,
                    v ^ 1,
                )
                .expect("ug_ext_gfa: the reverse overlap must exist");
                graph.push_arc(reverse);
            }
        }

        graph.index.clear();
        graph.is_sorted = false;
        graph.cleanup();
    }

    groups.len() as u32
}

/// `output_unitig_graph`: builds the unitigs, materialises them and writes the
/// `r_utg` GFA files.
#[allow(clippy::too_many_arguments)]
pub fn output_unitig_graph(
    read_graph: &mut Asg,
    coverage_cut: &[CoverageCut],
    output_file_name: &str,
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    maximum_hang: i32,
    minimum_overlap: i32,
    reads: &ReadStore,
    max_hang_rate: f32,
    bed_inconsistency_rate: u32,
) -> io::Result<()> {
    let edges: Vec<AsgArc> = Vec::new();
    let mut unitig_graph = ma_ug_gen(read_graph);
    ma_ug_seq(
        &mut unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        &edges,
        None,
        maximum_hang,
        minimum_overlap,
        true,
        reads,
        max_hang_rate,
    );

    let gfa_path = format!("{}.r_utg.gfa", output_file_name);
    let mut output = io::BufWriter::new(std::fs::File::create(&gfa_path)?);
    ma_ug_print(
        &unitig_graph,
        read_graph,
        Some(coverage_cut),
        Some(sources),
        Some(read_to_unitig),
        "utg",
        &mut output,
        reads,
        &reads.trio_labels,
    )?;
    output.flush()?;

    let noseq_path = format!("{}.r_utg.noseq.gfa", output_file_name);
    let mut output = io::BufWriter::new(std::fs::File::create(&noseq_path)?);
    ma_ug_print_simple(
        &unitig_graph,
        read_graph,
        Some(coverage_cut),
        Some(sources),
        Some(read_to_unitig),
        "utg",
        &mut output,
        reads,
        &reads.trio_labels,
    )?;
    output.flush()?;

    if bed_inconsistency_rate != 0 {
        let bed_path = format!("{}.r_utg.lowQ.bed", output_file_name);
        let mut output = io::BufWriter::new(std::fs::File::create(&bed_path)?);
        ma_ug_print_bed(
            &unitig_graph,
            read_graph,
            coverage_cut,
            sources,
            &edges,
            maximum_hang,
            minimum_overlap,
            bed_inconsistency_rate,
            "utg",
            &mut output,
            reads,
            max_hang_rate,
        )?;
        output.flush()?;
    }

    Ok(())
}

/// `output_contig_graph_primary_pre`: writes the primary-unitig GFA
/// (`p_utg.gfa`, `.noseq.gfa`, `.lowQ.bed`).
///
/// `bubble_dist` is `asm_opt.small_pop_bubble_size`; the reference pops small
/// unitig bubbles first when it is positive.
#[allow(clippy::too_many_arguments)]
pub fn output_contig_graph_primary_pre(
    read_graph: &mut Asg,
    coverage_cut: &[CoverageCut],
    output_file_name: &str,
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    bubble_dist: u64,
    maximum_hang: i32,
    minimum_overlap: i32,
    reads: &ReadStore,
    max_hang_rate: f32,
    bed_inconsistency_rate: u32,
) -> io::Result<()> {
    let edges: Vec<AsgArc> = Vec::new();
    let mut unitig_graph = ma_ug_gen_primary(read_graph, PRIMARY_LABEL);

    // Every live unitig is primary here; `EvaluateLen` is the unitig's
    // `start` field, which the bubble routines read as its read count.
    for v in 0..unitig_graph.graph.sequences.len() {
        if unitig_graph.graph.sequences[v].deleted {
            continue;
        }
        unitig_graph.graph.sequences[v].coverage_class = PRIMARY_LABEL;
        unitig_graph.unitigs[v].start = unitig_graph.unitigs[v].reads.len() as u32;
    }

    if bubble_dist > 0 {
        // `asg_pop_bubble_primary_trio(ug, &bubble_dist, (uint32_t)-1, DROP,
        // NULL, NULL, 0, NULL)`; the reference then drops the alternate
        // unitigs and rebuilds the unitig graph.
        asg_pop_bubble_primary_trio(
            &mut unitig_graph,
            Some(bubble_dist),
            u32::MAX,
            DROP,
            None,
            &reads.trio_labels,
            false,
            None,
        );
        delete_useless_nodes(&mut unitig_graph);
        renew_utg(&mut unitig_graph, read_graph, &edges);
    }

    ma_ug_seq(
        &mut unitig_graph,
        read_graph,
        coverage_cut,
        sources,
        &edges,
        None,
        maximum_hang,
        minimum_overlap,
        true,
        reads,
        max_hang_rate,
    );

    eprintln!("Writing processed unitig GFA to disk... ");

    let gfa_path = format!("{}.p_utg.gfa", output_file_name);
    let mut output = io::BufWriter::new(std::fs::File::create(&gfa_path)?);
    ma_ug_print(
        &unitig_graph,
        read_graph,
        Some(coverage_cut),
        Some(sources),
        Some(read_to_unitig),
        "utg",
        &mut output,
        reads,
        &reads.trio_labels,
    )?;
    output.flush()?;

    let noseq_path = format!("{}.p_utg.noseq.gfa", output_file_name);
    let mut output = io::BufWriter::new(std::fs::File::create(&noseq_path)?);
    ma_ug_print_simple(
        &unitig_graph,
        read_graph,
        Some(coverage_cut),
        Some(sources),
        Some(read_to_unitig),
        "utg",
        &mut output,
        reads,
        &reads.trio_labels,
    )?;
    output.flush()?;

    if bed_inconsistency_rate != 0 {
        let bed_path = format!("{}.p_utg.lowQ.bed", output_file_name);
        let mut output = io::BufWriter::new(std::fs::File::create(&bed_path)?);
        ma_ug_print_bed(
            &unitig_graph,
            read_graph,
            coverage_cut,
            sources,
            &edges,
            maximum_hang,
            minimum_overlap,
            bed_inconsistency_rate,
            "utg",
            &mut output,
            reads,
            max_hang_rate,
        )?;
        output.flush()?;
    }

    Ok(())
}

/// `hic_clean`: removes the alternate bubbles of a Hi-C-phased unitig graph
/// from the read graph.
///
/// The reference builds the primary unitig graph, finds every bubble whose
/// sides are made of at most two unitigs, and pops the small branches: a
/// branch is accepted when its interior holds at most three reads, at most
/// two unitigs, and fewer than a tenth of the enclosing bubble's reads. Every
/// read of an alternate unitig is then deleted from the read graph.
pub fn hic_clean(read_graph: &mut Asg) {
    let mut unitig_graph = ma_ug_gen_primary(read_graph, PRIMARY_LABEL);
    let vertex_count = unitig_graph.graph.sequences.len() * 2;
    let mut buffer = BubbleBuffer {
        info: vec![BubbleInfo::default(); vertex_count],
        ..BubbleBuffer::default()
    };
    let bubble_distance = get_s_bub_pop_max_dist_advance(&unitig_graph.graph, &mut buffer);
    let mut bubble_flag = vec![0u8; vertex_count];

    // Every live unitig becomes primary; `EvaluateLen` is the unitig's
    // `start` field, which is set to its read count.
    for v in 0..unitig_graph.graph.sequences.len() {
        if unitig_graph.graph.sequences[v].deleted {
            continue;
        }
        unitig_graph.graph.sequences[v].coverage_class = PRIMARY_LABEL;
        unitig_graph.unitigs[v].start = unitig_graph.unitigs[v].reads.len() as u32;
    }

    // First pass: locate the bubbles and flag their interior unitigs.
    for v in 0..vertex_count as u32 {
        if unitig_graph.graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if unitig_graph.graph.arc_count(v) < 2 {
            continue;
        }
        if bubble_flag[v as usize] != 0 {
            continue;
        }
        if asg_bub_pop1_primary_trio(
            &mut unitig_graph.graph,
            None,
            v,
            bubble_distance,
            &mut buffer,
            u32::MAX,
            u32::MAX,
            false,
            None,
            &[],
            false,
            false,
        ) != 0
        {
            let end = buffer.sources[0];
            for &node in &buffer.visited {
                if node == v || node == end {
                    continue;
                }
                bubble_flag[node as usize] = 1;
                bubble_flag[(node ^ 1) as usize] = 1;
            }
            bubble_flag[v as usize] = 2;
            bubble_flag[(end ^ 1) as usize] = 3;
        }
    }

    // Second pass: pop the small alternate branch of every bubble found.
    const BUBBLE_RATE: f64 = 0.1;
    let mut branch_unitigs: Vec<u32> = Vec::new();
    for v in 0..vertex_count as u32 {
        if bubble_flag[v as usize] != 2 {
            continue;
        }
        if asg_bub_pop1_primary_trio(
            &mut unitig_graph.graph,
            None,
            v,
            bubble_distance,
            &mut buffer,
            u32::MAX,
            u32::MAX,
            false,
            None,
            &[],
            false,
            false,
        ) == 0
        {
            continue;
        }
        branch_unitigs.clear();
        let mut primary_reads: u64 = 0;
        for &node in &buffer.visited {
            if node == v || node == buffer.sources[0] {
                continue;
            }
            primary_reads += unitig_graph.unitigs[(node >> 1) as usize].reads.len() as u64;
            branch_unitigs.push(node >> 1);
        }

        for i in 0..branch_unitigs.len() {
            let unitig = branch_unitigs[i];
            for k in 0..2u32 {
                let u = (unitig << 1) + k;
                if unitig_graph.graph.arc_count(u) < 2 {
                    continue;
                }
                if asg_bub_pop1_primary_trio(
                    &mut unitig_graph.graph,
                    None,
                    u,
                    bubble_distance,
                    &mut buffer,
                    u32::MAX,
                    u32::MAX,
                    false,
                    None,
                    &[],
                    false,
                    false,
                ) == 0
                {
                    continue;
                }
                let mut alternative_reads: u64 = 0;
                let mut alternative_unitigs: u64 = 0;
                for &node in &buffer.visited {
                    if node == u || node == buffer.sources[0] {
                        continue;
                    }
                    alternative_reads +=
                        unitig_graph.unitigs[(node >> 1) as usize].reads.len() as u64;
                    alternative_unitigs += 1;
                }
                if alternative_reads as f64 >= primary_reads as f64 * BUBBLE_RATE {
                    continue;
                }
                if alternative_reads > 3 {
                    continue;
                }
                if alternative_unitigs > 2 {
                    continue;
                }
                // The reference ignores the return value here; the pop is
                // recorded in the graph through `asg_bub_backtrack_primary`.
                asg_bub_pop1_primary_trio(
                    &mut unitig_graph.graph,
                    None,
                    u,
                    bubble_distance,
                    &mut buffer,
                    u32::MAX,
                    u32::MAX,
                    true,
                    None,
                    &[],
                    false,
                    false,
                );
            }
        }
    }

    // Delete every read of every alternate unitig.
    for v in 0..unitig_graph.graph.sequences.len() {
        if unitig_graph.graph.sequences[v].deleted {
            continue;
        }
        if unitig_graph.graph.sequences[v].coverage_class != ALTER_LABEL {
            continue;
        }
        let unitig = &unitig_graph.unitigs[v];
        if unitig.reads.is_empty() {
            continue;
        }
        for &entry in &unitig.reads {
            read_graph.delete_sequence((entry >> 33) as u32);
        }
    }

    read_graph.cleanup();
}

/// `copy_read_graph`: a deep copy of the read graph.
///
/// The reference allocates a fresh graph, copies the sequences and arcs,
/// drops the per-vertex index and runs `asg_cleanup`, so any deleted arc or
/// deleted sequence is gone from the copy. The sequence-visit flags are
/// copied when the source has them.
pub fn copy_read_graph(source: &Asg) -> Asg {
    let mut destination = Asg {
        arcs: source.arcs.clone(),
        // `asg_init` leaves `idx` null, so `asg_cleanup` rebuilds it.
        index: Vec::new(),
        is_sorted: source.is_sorted,
        is_symmetric: source.is_symmetric,
        sequences: source.sequences.clone(),
        read_sequence_count: source.read_sequence_count,
        sequence_visited: source.sequence_visited.clone(),
    };
    destination.cleanup();
    destination
}

/// `copy_untig_graph`: a deep copy of a unitig graph.
///
/// The unitig sequences (`ma_utg_t::s`) are not copied; the reference sets
/// them to null. The arc index is copied before `asg_cleanup`, matching the
/// reference, which only rebuilds it when removing deleted arcs invalidated
/// it.
pub fn copy_untig_graph(source: &UnitigGraph) -> UnitigGraph {
    let mut graph = Asg {
        arcs: source.graph.arcs.clone(),
        index: source.graph.index.clone(),
        is_sorted: source.graph.is_sorted,
        is_symmetric: source.graph.is_symmetric,
        sequences: source.graph.sequences.clone(),
        read_sequence_count: source.graph.read_sequence_count,
        sequence_visited: source.graph.sequence_visited.clone(),
    };
    graph.cleanup();

    let unitigs = source
        .unitigs
        .iter()
        .map(|unitig| MaUtg {
            length: unitig.length,
            circular: unitig.circular,
            start: unitig.start,
            end: unitig.end,
            reads: unitig.reads.clone(),
            sequence: Vec::new(),
        })
        .collect();

    UnitigGraph { unitigs, graph }
}

/// `fill_unitig`: rewrites a unitig's packed read list with the overlap
/// lengths between consecutive reads and returns the total length.
///
/// Every element of `buffer` is a `<vertex>:32 | <edge length>:32` word. For
/// each pair of neighbours the reference looks for the live arc `v -> w` in
/// the read graph and, failing that, in the extra `edges` list the caller
/// built while rescuing overlaps. The last element takes the arc back to the
/// first read when the unitig is circular, and the read's own length
/// otherwise. `total_length` is `uint32_t` in the reference, so it wraps.
pub fn fill_unitig(
    buffer: &mut [u64],
    read_graph: &Asg,
    edges: &[AsgArc],
    is_circle: bool,
) -> u64 {
    let mut total_length: u32 = 0;
    let last = buffer.len().saturating_sub(1);

    for i in 0..last {
        let v = (buffer[i] >> 32) as u32;
        let w = (buffer[i + 1] >> 32) as u32;
        let (length, found) = find_step_length(read_graph, edges, v, w);
        if !found {
            eprintln!(
                "####ERROR1-fill: i: {}, v>>1: {}, v&1: {}, w>>1: {}, w&1: {}",
                i,
                v >> 1,
                v & 1,
                w >> 1,
                w & 1
            );
        }
        buffer[i] = (u64::from(v) << 32) | u64::from(length);
        total_length = total_length.wrapping_add(length);
    }

    if last < buffer.len() {
        let v = (buffer[last] >> 32) as u32;
        let length = if is_circle {
            let w = (buffer[0] >> 32) as u32;
            let (length, found) = find_step_length(read_graph, edges, v, w);
            if !found {
                eprintln!(
                    "####ERROR2-fill: i: {}, v>>1: {}, v&1: {}, w>>1: {}, w&1: {}",
                    last,
                    v >> 1,
                    v & 1,
                    w >> 1,
                    w & 1
                );
            }
            length
        } else {
            read_graph.sequences[(v >> 1) as usize].length
        };
        buffer[last] = (u64::from(v) << 32) | u64::from(length);
        total_length = total_length.wrapping_add(length);
    }

    u64::from(total_length)
}

/// The `v -> w` step length used by [`fill_unitig`].
fn find_step_length(read_graph: &Asg, edges: &[AsgArc], v: u32, w: u32) -> (u32, bool) {
    for arc in read_graph.arcs_of(v) {
        if arc.deleted {
            continue;
        }
        if arc.v == w {
            return (arc.node_length(), true);
        }
    }
    for arc in edges {
        if arc.deleted {
            continue;
        }
        if arc.query_vertex() == v && arc.v == w {
            return (arc.node_length(), true);
        }
    }
    (0, false)
}

/// `merge_unitig_content`: replaces a high-level unitig's read list with the
/// concatenation of its child unitigs' read lists.
///
/// A child visited on the reverse strand (`ori == 1`) contributes its reads
/// backwards, with the direction bit of every entry flipped. The child read
/// arrays are freed (cleared) in the old unitig graph. The new read list then
/// goes through [`fill_unitig`], which fills in the overlap lengths, and the
/// unitig's length and endpoints are refreshed. The `circ` flag is set when
/// the only child was itself circular.
pub fn merge_unitig_content(
    collection: &mut MaUtg,
    unitigs: &mut [MaUtg],
    read_graph: &Asg,
    edges: &[AsgArc],
) {
    if collection.reads.is_empty() {
        return;
    }

    let mut index = 0usize;
    for &entry in &collection.reads {
        index += unitigs[(entry >> 33) as usize].reads.len();
    }

    let mut buffer: Vec<u64> = vec![0; index];
    let mut position = 0usize;
    for &entry in &collection.reads {
        let child = (entry >> 33) as usize;
        let orientation = (entry >> 32) & 1;
        let query = &mut unitigs[child];
        if orientation == 1 {
            let count = query.reads.len();
            for j in 0..count {
                buffer[position + count - j - 1] = query.reads[j] ^ 0x1_0000_0000;
            }
        } else {
            let count = query.reads.len();
            buffer[position..position + count].copy_from_slice(&query.reads);
        }
        position += query.reads.len();
        query.reads.clear();
    }

    if index == 0 {
        return;
    }

    let child_is_circular =
        collection.reads.len() == 1 && unitigs[(collection.reads[0] >> 33) as usize].circular;
    let total_length = fill_unitig(&mut buffer, read_graph, edges, child_is_circular);
    // `important. must be here`: the flag is raised before the endpoints are
    // computed, so a freshly circular unitig gets the UINT32_MAX endpoints.
    if child_is_circular {
        collection.circular = true;
    }

    collection.reads = buffer;
    collection.length = total_length as u32 & 0x7fff_ffff;
    if !collection.circular {
        collection.start = (collection.reads[0] >> 32) as u32;
        collection.end = ((collection.reads[collection.reads.len() - 1] >> 32) as u32) ^ 1;
    } else {
        collection.start = u32::MAX;
        collection.end = u32::MAX;
    }
}

/// `renew_utg`: rebuilds the unitig graph from its own string graph and
/// re-merges the reads of every surviving unitig.
///
/// The new unitigs all start out as primary. A unitig whose read list is
/// empty keeps its own entry; every other one gets its children's reads
/// concatenated by [`merge_unitig_content`] and its length refreshed.
pub fn renew_utg(unitig_graph: &mut UnitigGraph, read_graph: &Asg, edges: &[AsgArc]) {
    let mut high_level = ma_ug_gen(&mut unitig_graph.graph);
    let old_unitigs = &mut unitig_graph.unitigs;

    for i in 0..high_level.unitigs.len() {
        high_level.graph.sequences[i].coverage_class = PRIMARY_LABEL;
        if high_level.unitigs[i].reads.is_empty() {
            continue;
        }
        merge_unitig_content(&mut high_level.unitigs[i], old_unitigs, read_graph, edges);
        high_level.graph.sequences[i].length = high_level.unitigs[i].length;
    }

    *unitig_graph = high_level;
}

/// `delete_useless_nodes`: removes every unitig whose coverage class is
/// `ALTER_LABLE` from the unitig graph.
///
/// The reference marks the unitig's vertex deleted with `asg_seq_del` and
/// drops its read list, then cleans the graph up.
pub fn delete_useless_nodes(unitig_graph: &mut UnitigGraph) {
    let vertex_count = unitig_graph.graph.sequences.len();
    for v in 0..vertex_count {
        if unitig_graph.graph.sequences[v].deleted {
            continue;
        }
        if unitig_graph.graph.sequences[v].coverage_class == ALTER_LABEL {
            unitig_graph.graph.delete_sequence(v as u32);
            unitig_graph.unitigs[v].reads.clear();
        }
    }
    unitig_graph.graph.cleanup();
}

/// `discard_small_ctg`: removes every unitig whose read count is at most
/// `maximum_read_count`.
///
/// `renew_before`/`renew_after` are the reference's `is_renew0`/`is_renew1`;
/// the second renewal runs only when at least one unitig was removed.
pub fn discard_small_ctg(
    unitig_graph: &mut UnitigGraph,
    read_graph: &Asg,
    edges: &[AsgArc],
    maximum_read_count: i32,
    renew_before: bool,
    renew_after: bool,
) {
    if renew_before {
        renew_utg(unitig_graph, read_graph, edges);
    }

    let mut removed = 0u32;
    for i in 0..unitig_graph.unitigs.len() {
        let read_count = unitig_graph.unitigs[i].reads.len() as u32;
        // `(int32_t)u->n > max_rg_cut`: the reference casts the unsigned
        // read count to a signed integer.
        if (read_count as i32) > maximum_read_count {
            continue;
        }
        unitig_graph.graph.delete_sequence(i as u32);
        // `if (u->m != 0) { u->m = u->n = 0; free(u->a); ... }`: a unitig has
        // a nonzero capacity exactly when it has reads, so clearing the
        // whole record is the equivalent of the reference's field resets.
        unitig_graph.unitigs[i] = MaUtg::default();
        removed += 1;
    }

    if removed != 0 && renew_after {
        renew_utg(unitig_graph, read_graph, edges);
    }
}

/// `append_utg`: moves every live unitig of `atg` behind the unitigs of `ptg`.
///
/// The reference reallocates the unitig array up front, then transfers each
/// unitig's read list and materialised sequence by pointer. `is_r_het` stands
/// in for `t_ch->ir_het`, which aliases `hap_cov_t::is_r_het`; when given,
/// the reads of every appended unitig are reset to `N_HET`.
pub fn append_utg(ptg: &mut UnitigGraph, atg: &mut UnitigGraph, mut is_r_het: Option<&mut [u8]>) {
    let sequence_count = atg.graph.sequences.len();
    let mut node_count = 0u32;
    for v in 0..sequence_count {
        if atg.graph.sequences[v].deleted || atg.unitigs[v].reads.is_empty() {
            continue;
        }
        node_count += 1;
    }
    if node_count == 0 {
        return;
    }

    for v in 0..sequence_count {
        if atg.graph.sequences[v].deleted || atg.unitigs[v].reads.is_empty() {
            continue;
        }
        if let Some(flags) = is_r_het.as_deref_mut() {
            reset_trans_chain(flags, &atg.unitigs[v]);
        }

        let source = &mut atg.unitigs[v];
        let appended = MaUtg {
            length: source.length,
            circular: source.circular,
            start: source.start,
            end: source.end,
            reads: std::mem::take(&mut source.reads),
            sequence: std::mem::take(&mut source.sequence),
        };

        let index = ptg.unitigs.len();
        ptg.graph.set_sequence(index, appended.length, false);
        ptg.unitigs.push(appended);
    }

    ptg.graph.index.clear();
    ptg.graph.cleanup();
}

/// `recover_utg_by_coverage`: re-adds alternate unitigs whose reads are
/// covered densely enough by alternate-class reads.
///
/// `minimum`/`maximum` are `asm_opt.recover_atg_cov_min`/`_max`; the default
/// `-1024` range is resolved by the caller before this function runs. The
/// function returns early when either bound is `-1` or the interval is
/// inverted, exactly like the reference.
#[allow(clippy::too_many_arguments)]
pub fn recover_utg_by_coverage(
    ptg: &mut UnitigGraph,
    read_graph: &mut Asg,
    coverage_cut: &[CoverageCut],
    sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    minimum: i32,
    maximum: i32,
    is_r_het: Option<&mut [u8]>,
) {
    if minimum == -1 || maximum == -1 || minimum > maximum {
        return;
    }

    let mut atg = ma_ug_gen_primary(read_graph, ALTER_LABEL);
    let mut keep_atg = 0u32;

    for v in 0..atg.graph.sequences.len() {
        if atg.graph.sequences[v].deleted {
            continue;
        }
        if atg.unitigs[v].reads.is_empty() {
            continue;
        }

        let mut total_c_bases: i64 = 0;
        let mut total_c_bases_alter: i64 = 0;
        let mut available_reads = 0u32;

        for read_index in 0..atg.unitigs[v].reads.len() {
            let read_id = (atg.unitigs[v].reads[read_index] >> 33) as usize;
            let mut c_bases_primary: i64 = 0;
            let mut c_bases_alter: i64 = 0;
            // `R_bases = coverage_cut[rId].e - coverage_cut[rId].s` is an
            // unsigned 32-bit subtraction widened to long long.
            let r_bases = i64::from(
                coverage_cut[read_id].end.wrapping_sub(coverage_cut[read_id].start),
            );

            for hit in &sources[read_id].hits {
                if !hit.exact_region {
                    continue;
                }
                let mut target = hit.target_identifier;
                if read_graph.sequences[target as usize].deleted {
                    let (mapped, is_unitig) = read_to_unitig.get(target);
                    if mapped == u32::MAX
                        || is_unitig == 1
                        || read_graph.sequences[mapped as usize].deleted
                    {
                        continue;
                    }
                    target = mapped;
                }
                if read_graph.sequences[target as usize].deleted {
                    continue;
                }
                let covered = i64::from(hit.query_end) - i64::from(hit.query_start);
                if read_graph.sequences[target as usize].coverage_class == ALTER_LABEL {
                    c_bases_alter += covered;
                } else {
                    c_bases_primary += covered;
                }
            }

            let c_bases = c_bases_primary + c_bases_alter;
            total_c_bases += c_bases;
            // `C_bases = C_bases / R_bases` on long long operands.
            let c_bases = c_bases / r_bases;
            if c_bases >= i64::from(minimum) && c_bases <= i64::from(maximum) {
                if c_bases_alter as f64
                    >= ((c_bases_primary + c_bases_alter) as f64) * ALTER_COVERAGE_THRESHOLD
                {
                    available_reads += 1;
                }
                total_c_bases_alter += c_bases_alter;
            }
        }

        let read_count = atg.unitigs[v].reads.len() as u32;
        let sparse_alternate_reads = (available_reads as f64) < (read_count as f64) * 0.8
            && (total_c_bases_alter as f64) < (total_c_bases as f64) * 0.8;
        if sparse_alternate_reads || available_reads == 0 {
            atg.graph.delete_sequence(v as u32);
            atg.unitigs[v] = MaUtg::default();
        } else {
            keep_atg += 1;
        }
    }

    if keep_atg > 0 {
        atg.graph.cleanup();
        atg.graph.symmetrize();
        append_utg(ptg, &mut atg, is_r_het);

        for sequence in &mut read_graph.sequences {
            sequence.coverage_class = ALTER_LABEL;
        }

        let graph = &ptg.graph;
        for v in 0..graph.sequences.len() {
            if graph.sequences[v].deleted {
                continue;
            }
            if ptg.unitigs[v].reads.is_empty() {
                continue;
            }
            let class = graph.sequences[v].coverage_class;
            for &packed in &ptg.unitigs[v].reads {
                read_graph.sequences[(packed >> 33) as usize].coverage_class = class;
            }
        }

        for v in 0..read_graph.sequences.len() {
            if read_graph.sequences[v].coverage_class == ALTER_LABEL {
                asg_seq_drop(read_graph, v as u32);
            }
        }
    }

    // `ma_ug_destroy(atg)`: `atg` is dropped here.
}

/// `collect_ma_utg_ts`: appends the reads of unitig `w` to `result`.
///
/// When `w` points at the reverse orientation the reads are appended in
/// reverse order with the orientation bit of every packed word flipped. The
/// first call (`v == u32::MAX`) seeds the summary fields from `w`; only the
/// read list survives to the caller.
pub fn collect_ma_utg_ts(unitigs: &[MaUtg], v: u32, w: u32, result: &mut MaUtg) {
    let destination_reverse = (w & 1) == 1;
    let unitig = &unitigs[(w >> 1) as usize];

    if v == u32::MAX {
        result.start = unitig.start;
        result.circular = unitig.circular;
        result.end = unitig.end;
        result.length = unitig.length;
        result.reads.clear();
    }

    if destination_reverse {
        for i in (0..unitig.reads.len()).rev() {
            result.reads.push(unitig.reads[i] ^ 0x1_0000_0000);
        }
    } else {
        result.reads.extend_from_slice(&unitig.reads);
    }
}

/// `merge_ug_nodes`: merges the unitig chain `array` into one unitig.
///
/// `array` starts as the unitig chain collected by [`get_unitig`] and is
/// extended here with one `<new unitig>:32 | <overlap>:32,
/// <overlap>:32 | <target vertex>:32` pair per edge of the merged node. The
/// function walks `array` while it grows, exactly like the reference.
///
/// The reference only clears a chain member when its `m` (allocation
/// capacity) is nonzero; every unitig produced by [`ma_ug_gen`] or copied by
/// [`copy_untig_graph`] has at least one read when `m != 0`, so an empty read
/// list is used for that test here. The materialised sequence (`s`) is kept,
/// as in the reference.
pub fn merge_ug_nodes(
    unitigs: &mut [MaUtg],
    graph: &mut Asg,
    read_graph: &Asg,
    array: &mut Vec<u64>,
) {
    if array.is_empty() {
        return;
    }

    let mut result = MaUtg::default();
    let mut v = u32::MAX;
    for index in 0..array.len() {
        let w = array[index] as u32;
        collect_ma_utg_ts(unitigs, v, w, &mut result);
        v = w;
    }
    if result.reads.is_empty() {
        return;
    }

    // Recompute each step length from the read graph, then the merged length.
    result.length = 0;
    let mut i = 0usize;
    while i + 1 < result.reads.len() {
        let v = (result.reads[i] >> 32) as u32;
        let w = (result.reads[i + 1] >> 32) as u32;
        let mut length = 0u32;
        for arc in read_graph.arcs_of(v) {
            if arc.deleted {
                continue;
            }
            if arc.v == w {
                length = arc.node_length();
                break;
            }
        }
        result.reads[i] = (u64::from(v) << 32) | u64::from(length);
        result.length = result.length.wrapping_add(length);
        i += 1;
    }
    if i < result.reads.len() {
        let v = (result.reads[i] >> 32) as u32;
        let length = read_graph.sequences[(v >> 1) as usize].length;
        result.reads[i] = (u64::from(v) << 32) | u64::from(length);
        result.length = result.length.wrapping_add(length);
    }

    result.circular = false;
    result.start = (result.reads[0] >> 32) as u32;
    result.end = ((result.reads[result.reads.len() - 1] >> 32) as u32) ^ 1;

    let beg_uid = array[0] as u32;
    let end_uid = array[array.len() - 1] as u32;
    let real_length = array.len();
    let new_uid = (array[0] >> 1) as u32;

    // Edges entering the new node: walk the reverse arcs of `beg_uid`.
    let v = beg_uid ^ 1;
    let arcs: Vec<AsgArc> = graph.arcs_of(v).to_vec();
    for arc in arcs {
        if arc.deleted {
            continue;
        }
        let mut key = u64::from(new_uid.wrapping_shl(1) | 1) << 32;
        let w = arc.v ^ 1;
        let reverse = graph
            .arcs_of(w)
            .iter()
            .find(|candidate| !candidate.deleted && candidate.v == (v ^ 1));
        // The reference indexes `aw[i]` even when the search failed; a zero
        // overlap keeps the port memory-safe without changing valid inputs.
        let overlap = reverse.map_or(0, |candidate| candidate.overlap);
        key |= u64::from(overlap);
        array.push(key);
        array.push((u64::from(arc.overlap) << 32) | u64::from(arc.v));
    }

    // Edges leaving the new node: walk the arcs of `end_uid`.
    let v = end_uid;
    let arcs: Vec<AsgArc> = graph.arcs_of(v).to_vec();
    for arc in arcs {
        if arc.deleted {
            continue;
        }
        let mut key = u64::from(new_uid.wrapping_shl(1)) << 32;
        let w = arc.v ^ 1;
        let reverse = graph
            .arcs_of(w)
            .iter()
            .find(|candidate| !candidate.deleted && candidate.v == (v ^ 1));
        let overlap = reverse.map_or(0, |candidate| candidate.overlap);
        key |= u64::from(overlap);
        array.push(key);
        array.push((u64::from(arc.overlap) << 32) | u64::from(arc.v));
    }

    // Remove every unitig of the chain and let the first one host the merge.
    for index in 0..real_length {
        let w = array[index] as u32;
        let unitig = &mut unitigs[(w >> 1) as usize];
        if !unitig.reads.is_empty() {
            unitig.start = 0;
            unitig.end = 0;
            unitig.length = 0;
            unitig.circular = false;
            unitig.reads.clear();
        }
        graph.delete_sequence(w >> 1);
    }

    unitigs[(beg_uid >> 1) as usize] = result;
    let sequence_index = (beg_uid >> 1) as usize;
    graph.sequences[sequence_index].deleted = false;
    graph.sequences[sequence_index].length = unitigs[sequence_index].length;

    let mut index = real_length;
    while index + 1 < array.len() {
        let v = (array[index] >> 32) as u32;
        let mut w = array[index + 1] as u32;
        if (w >> 1) == (beg_uid >> 1) || (w >> 1) == (end_uid >> 1) {
            w = v;
        }

        let overlap = (array[index + 1] >> 32) as u32;
        graph.append_edge_to_sorted(
            v,
            unitigs[(v >> 1) as usize].length,
            w,
            overlap,
            false,
            0,
            false,
        );
        let reverse_overlap = array[index] as u32;
        graph.append_edge_to_sorted(
            w ^ 1,
            unitigs[(w >> 1) as usize].length,
            v ^ 1,
            reverse_overlap,
            false,
            0,
            false,
        );
        index += 2;
    }
}

/// `unroll_simple_case_advance`: merges the short unitig chains that connect
/// two simple neighbouring unitigs.
///
/// Every merge goes through [`merge_ug_nodes`]; the loop repeats until a full
/// pass makes no change.
#[allow(clippy::too_many_arguments)]
pub fn unroll_simple_case_advance(
    unitig_graph: &mut UnitigGraph,
    read_graph: &Asg,
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    mut labels: Option<&mut BubbleLabel>,
    is_r_het: &[u8],
    duplicate_length_threshold: f64,
    purge_similarity_threshold: f32,
) {
    let UnitigGraph { unitigs, graph } = unitig_graph;
    let vertex_count = graph.sequences.len() as u32 * 2;
    let mut n_reduce = 1u32;
    let mut unitig_vertices: Vec<u64> = Vec::new();
    let mut buffer_0: Vec<u32> = Vec::new();
    let mut buffer_1: Vec<u32> = Vec::new();

    if let Some(label) = labels.as_deref_mut() {
        let bubble_distance = get_s_bub_pop_max_dist_advance(graph, &mut label.buffers[0]);
        label.reset(graph.sequences.len(), bubble_distance, 0);
    }

    while n_reduce > 0 {
        n_reduce = 0;
        // Break nearly circular paths; the reference first cleans up, so the
        // vertex count stays what it was at entry.
        n_reduce += asg_arc_del_simple_circle_untig(graph, 100, false);

        for v in 0..vertex_count {
            if graph.sequences[(v >> 1) as usize].deleted {
                continue;
            }
            let v_left = v;
            if graph.arc_count(v_left) < 1 {
                continue;
            }
            if get_real_length(graph, v_left) != 1 {
                continue;
            }

            let mut v_right = u32::MAX;
            let mut node_length = 0i64;
            let mut reference_base = 0i64;
            let mut max_stop_node_length = 0i64;
            let mut max_stop_base_length = 0i64;
            let return_flag = get_unitig(
                graph,
                None,
                v_left ^ 1,
                &mut v_right,
                &mut node_length,
                &mut reference_base,
                &mut max_stop_node_length,
                &mut max_stop_base_length,
                1,
                None,
            );
            if return_flag == LOOP {
                continue;
            }
            if return_flag != MUL_INPUT {
                continue;
            }
            if graph.arc_count(v_right) < 1 {
                continue;
            }
            if get_real_length(graph, v_right) != 1 {
                continue;
            }

            let w_left = get_real_length_one(graph, v_left)
                .expect("a single live arc always reports its target");
            let w_right = get_real_length_one(graph, v_right)
                .expect("a single live arc always reports its target");
            if (v_left >> 1) == (w_left >> 1) || (v_left >> 1) == (w_right >> 1) {
                continue;
            }
            if (v_right >> 1) == (w_left >> 1) || (v_right >> 1) == (w_right >> 1) {
                continue;
            }

            if w_left != w_right {
                let mut convex = u32::MAX;
                let mut duplicate_base = 0i64;
                let return_flag = get_unitig(
                    graph,
                    None,
                    w_left,
                    &mut convex,
                    &mut node_length,
                    &mut duplicate_base,
                    &mut max_stop_node_length,
                    &mut max_stop_base_length,
                    1,
                    None,
                );
                if return_flag == LOOP {
                    continue;
                }
                if return_flag != MUL_OUTPUT {
                    continue;
                }
                if convex != (w_right ^ 1) {
                    continue;
                }

                let mut beg = u32::MAX;
                let mut end = u32::MAX;
                let mut live_count = 0u32;
                for arc in graph.arcs_of(w_left ^ 1) {
                    if arc.deleted {
                        continue;
                    }
                    live_count += 1;
                    if arc.v == (v_left ^ 1) {
                        continue;
                    }
                    beg = arc.v;
                }
                if live_count != 2 {
                    continue;
                }

                live_count = 0;
                for arc in graph.arcs_of(w_right ^ 1) {
                    if arc.deleted {
                        continue;
                    }
                    live_count += 1;
                    if arc.v == (v_right ^ 1) {
                        continue;
                    }
                    end = arc.v;
                }
                if live_count != 2 {
                    continue;
                }

                if get_real_length(graph, beg ^ 1) != 1 {
                    continue;
                }
                if get_real_length(graph, end ^ 1) != 1 {
                    continue;
                }
                if (beg >> 1) == (end >> 1) {
                    continue;
                }

                unitig_vertices.clear();
                unitig_vertices.push(u64::from(beg ^ 1));

                buffer_0.clear();
                get_unitig(
                    graph,
                    None,
                    w_left,
                    &mut convex,
                    &mut node_length,
                    &mut duplicate_base,
                    &mut max_stop_node_length,
                    &mut max_stop_base_length,
                    1,
                    Some(&mut buffer_0),
                );
                unitig_vertices.extend(buffer_0.iter().map(|&value| u64::from(value)));

                buffer_0.clear();
                get_unitig(
                    graph,
                    None,
                    v_right ^ 1,
                    &mut convex,
                    &mut node_length,
                    &mut duplicate_base,
                    &mut max_stop_node_length,
                    &mut max_stop_base_length,
                    1,
                    Some(&mut buffer_0),
                );
                unitig_vertices.extend(buffer_0.iter().map(|&value| u64::from(value)));

                buffer_0.clear();
                get_unitig(
                    graph,
                    None,
                    w_left,
                    &mut convex,
                    &mut node_length,
                    &mut duplicate_base,
                    &mut max_stop_node_length,
                    &mut max_stop_base_length,
                    1,
                    Some(&mut buffer_0),
                );
                unitig_vertices.extend(buffer_0.iter().map(|&value| u64::from(value)));

                unitig_vertices.push(u64::from(end));
                merge_ug_nodes(unitigs, graph, read_graph, &mut unitig_vertices);
                n_reduce += 1;
            } else {
                if get_real_length(graph, w_left ^ 1) != 2 {
                    continue;
                }
                let mut convex = u32::MAX;
                let mut duplicate_base = 0i64;
                let return_flag = get_unitig(
                    graph,
                    None,
                    w_left,
                    &mut convex,
                    &mut node_length,
                    &mut duplicate_base,
                    &mut max_stop_node_length,
                    &mut max_stop_base_length,
                    1,
                    None,
                );
                if return_flag == LOOP {
                    continue;
                }
                if return_flag != MUL_OUTPUT {
                    continue;
                }
                if get_real_length(graph, convex) != 2 {
                    continue;
                }
                if duplicate_base as f64 >= reference_base as f64 * duplicate_length_threshold {
                    continue;
                }

                let mut beg = u32::MAX;
                let mut end = u32::MAX;
                let mut live_count = 0u32;
                for arc in graph.arcs_of(convex) {
                    if arc.deleted {
                        continue;
                    }
                    if live_count == 0 {
                        beg = arc.v;
                    }
                    if live_count == 1 {
                        end = arc.v;
                    }
                    live_count += 1;
                }
                if live_count != 2 {
                    continue;
                }
                if (beg >> 1) == (end >> 1) {
                    continue;
                }
                if get_real_length(graph, beg ^ 1) != 1 {
                    continue;
                }
                if get_real_length(graph, end ^ 1) != 1 {
                    continue;
                }

                if let Some(label) = labels.as_deref_mut() {
                    let bubble_distance = label.bubble_distance;
                    if asg_bub_pop1_label(graph, convex, bubble_distance, &mut label.buffers[0]) != 0
                    {
                        continue;
                    }
                }

                if check_different_haps(
                    graph,
                    Some(&unitigs[..]),
                    read_graph,
                    beg,
                    end,
                    reverse_sources,
                    &mut buffer_0,
                    &mut buffer_1,
                    read_to_unitig,
                    is_r_het,
                    2,
                    1,
                    purge_similarity_threshold,
                ) == PLOID
                {
                    continue;
                }

                let w_term = convex ^ 1;
                unitig_vertices.clear();
                unitig_vertices.push(u64::from(beg ^ 1));

                buffer_0.clear();
                get_unitig(
                    graph,
                    None,
                    w_term,
                    &mut convex,
                    &mut node_length,
                    &mut duplicate_base,
                    &mut max_stop_node_length,
                    &mut max_stop_base_length,
                    1,
                    Some(&mut buffer_0),
                );
                unitig_vertices.extend(buffer_0.iter().map(|&value| u64::from(value)));

                buffer_0.clear();
                get_unitig(
                    graph,
                    None,
                    v_left ^ 1,
                    &mut convex,
                    &mut node_length,
                    &mut duplicate_base,
                    &mut max_stop_node_length,
                    &mut max_stop_base_length,
                    1,
                    Some(&mut buffer_0),
                );
                unitig_vertices.extend(buffer_0.iter().map(|&value| u64::from(value)));

                buffer_0.clear();
                get_unitig(
                    graph,
                    None,
                    w_left,
                    &mut convex,
                    &mut node_length,
                    &mut duplicate_base,
                    &mut max_stop_node_length,
                    &mut max_stop_base_length,
                    1,
                    Some(&mut buffer_0),
                );
                unitig_vertices.extend(buffer_0.iter().map(|&value| u64::from(value)));

                unitig_vertices.push(u64::from(end));
                merge_ug_nodes(unitigs, graph, read_graph, &mut unitig_vertices);
                n_reduce += 1;
            }
        }
    }
}

/// `adjust_utg_advance`: unrolls simple cases, drops semi-circles and
/// normalises the unitig graph.
#[allow(clippy::too_many_arguments)]
pub fn adjust_utg_advance(
    read_graph: &Asg,
    unitig_graph: &mut UnitigGraph,
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    labels: Option<&mut BubbleLabel>,
    is_r_het: &[u8],
    purge_similarity_threshold: f32,
) {
    unroll_simple_case_advance(
        unitig_graph,
        read_graph,
        reverse_sources,
        read_to_unitig,
        labels,
        is_r_het,
        2.5,
        purge_similarity_threshold,
    );
    drop_semi_circle(
        unitig_graph,
        read_graph,
        reverse_sources,
        read_to_unitig,
        is_r_het,
        purge_similarity_threshold,
    );
    unitig_graph.graph.cleanup();
    unitig_graph.graph.symmetrize();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::asg::AsgSeq;

    fn chain_graph(lengths: &[u32], overlaps: &[u32]) -> Asg {
        let mut graph = Asg::new();
        for (i, &length) in lengths.iter().enumerate() {
            graph.set_sequence(i, length, false);
            graph.sequences[i].coverage_class = PRIMARY_LABEL;
        }
        // A real string graph is symmetric: `0 -> 2` always comes with
        // `3 -> 1`, and the unitig walk relies on that.
        for (i, &overlap) in overlaps.iter().enumerate() {
            let query = (i as u32) << 1;
            let target = ((i + 1) as u32) << 1;
            let node_length = lengths[i] - overlap;
            graph.push_arc(AsgArc {
                ul: (u64::from(query) << 32) | u64::from(node_length),
                v: target,
                overlap,
                ..AsgArc::default()
            });
            graph.push_arc(AsgArc {
                ul: (u64::from(target ^ 1) << 32) | u64::from(node_length),
                v: query ^ 1,
                overlap,
                ..AsgArc::default()
            });
        }
        graph.cleanup();
        graph
    }

    #[test]
    fn a_linear_chain_becomes_one_unitig() {
        let mut graph = chain_graph(&[100, 200, 300], &[10, 20]);
        let unitigs = ma_ug_gen(&mut graph);
        assert_eq!(unitigs.unitigs.len(), 1);
        let unitig = &unitigs.unitigs[0];
        assert!(!unitig.circular);
        assert_eq!(unitig.reads.len(), 3);
        // 90 + 180 + 300
        assert_eq!(unitig.length, 570);
        assert_eq!(unitig.start, 0);
        assert_eq!(unitig.end, 5);
        assert_eq!(unitig_graph_arc_count(&unitigs), 0);
    }

    fn unitig_graph_arc_count(unitigs: &UnitigGraph) -> usize {
        unitigs.graph.arcs.len()
    }

    #[test]
    fn a_fork_leaves_three_singleton_unitigs() {
        let mut graph = Asg::new();
        for i in 0..3 {
            graph.set_sequence(i, 100, false);
        }
        // 0 -> 2 and 0 -> 4: vertex 0 has two outgoing arcs.
        for target in [2u32, 4] {
            graph.push_arc(AsgArc {
                ul: 50,
                v: target,
                overlap: 10,
                ..AsgArc::default()
            });
        }
        graph.cleanup();
        let unitigs = ma_ug_gen(&mut graph);
        assert_eq!(unitigs.unitigs.len(), 3);
        assert_eq!(unitigs.graph.arcs.len(), 2);
    }

    #[test]
    fn the_coverage_class_filter_skips_alternate_reads() {
        // Two unconnected reads, the first one marked as an alternate.
        let mut graph = Asg::new();
        graph.set_sequence(0, 100, false);
        graph.set_sequence(1, 100, false);
        graph.sequences[0].coverage_class = ALTER_LABEL;
        graph.sequences[1].coverage_class = PRIMARY_LABEL;
        graph.cleanup();

        let unitigs = ma_ug_gen_primary(&mut graph, PRIMARY_LABEL);
        assert_eq!(unitigs.unitigs.len(), 1);
        assert_eq!(unitigs.unitigs[0].reads.len(), 1);
        assert_eq!((unitigs.unitigs[0].reads[0] >> 32) as u32, 2);
    }

    #[test]
    fn the_complement_table_matches_the_reference() {
        assert_eq!(complement_byte(b'A'), b'T');
        assert_eq!(complement_byte(b'C'), b'G');
        assert_eq!(complement_byte(b'G'), b'C');
        assert_eq!(complement_byte(b'T'), b'A');
        assert_eq!(complement_byte(b'N'), b'N');
        assert_eq!(complement_byte(b'a'), b't');
        assert_eq!(complement_byte(200), b'N');
        assert_eq!(complement_byte(91), 91);
    }

    #[test]
    fn the_sequence_of_a_forward_unitig_is_concatenated() {
        use crate::reads::ReadStore;

        let mut store = ReadStore::default();
        store.read_names = vec!["r1".to_string(), "r2".to_string()];
        store.read_lengths = vec![8, 8];
        store.sequences = vec![
            crate::sequence::encoding::PackedSequence::from_ascii(b"ACGTACGT"),
            crate::sequence::encoding::PackedSequence::from_ascii(b"TTTTGGGG"),
        ];
        store.trio_labels = vec![AMBIGU, AMBIGU];
        store.total_reads = 2;
        store.total_bases = 16;

        let coverage_cut = vec![
            CoverageCut {
                start: 0,
                end: 8,
                deleted: false,
                class: PRIMARY_LABEL,
            };
            2
        ];

        let sources = vec![PafList::default(); 2];

        let mut graph = Asg::new();
        graph.set_sequence(0, 6, false);
        graph.set_sequence(1, 6, false);
        graph.push_arc(AsgArc {
            ul: (0u64 << 32) | 6,
            v: 2,
            overlap: 2,
            ..AsgArc::default()
        });
        graph.push_arc(AsgArc {
            ul: (3u64 << 32) | 6,
            v: 1,
            overlap: 2,
            ..AsgArc::default()
        });
        graph.cleanup();

        let mut unitig_graph = ma_ug_gen(&mut graph);
        let edges = Vec::new();
        ma_ug_seq(
            &mut unitig_graph,
            &mut graph,
            &coverage_cut,
            &sources,
            &edges,
            None,
            1000,
            50,
            false,
            &store,
            0.8,
        );
        assert_eq!(unitig_graph.unitigs[0].sequence, b"ACGTACTTTTGG");
    }

    #[test]
    fn unitig_sequences_are_reverse_complemented_on_the_reverse_strand() {
        use crate::reads::ReadStore;

        let mut store = ReadStore::default();
        store.read_names = vec!["r1".to_string(), "r2".to_string()];
        store.read_lengths = vec![4, 4];
        store.sequences = vec![
            crate::sequence::encoding::PackedSequence::from_ascii(b"ACGT"),
            crate::sequence::encoding::PackedSequence::from_ascii(b"AAAA"),
        ];
        store.trio_labels = vec![AMBIGU, AMBIGU];
        store.total_reads = 2;
        store.total_bases = 8;

        let coverage_cut = vec![
            CoverageCut {
                start: 0,
                end: 4,
                deleted: false,
                class: PRIMARY_LABEL,
            };
            2
        ];
        let sources = vec![PafList::default(); 2];

        let mut graph = Asg::new();
        graph.set_sequence(0, 4, false);
        graph.set_sequence(1, 4, false);
        // 0 -> (1^1) = 3: the second read enters on the reverse strand. The
        // symmetric counterpart is 2 -> 1.
        graph.push_arc(AsgArc {
            ul: (0u64 << 32) | 4,
            v: 3,
            overlap: 0,
            ..AsgArc::default()
        });
        graph.push_arc(AsgArc {
            ul: (2u64 << 32) | 4,
            v: 1,
            overlap: 0,
            ..AsgArc::default()
        });
        graph.cleanup();

        let mut unitig_graph = ma_ug_gen(&mut graph);
        let edges = Vec::new();
        ma_ug_seq(
            &mut unitig_graph,
            &mut graph,
            &coverage_cut,
            &sources,
            &edges,
            None,
            1000,
            50,
            false,
            &store,
            0.8,
        );
        assert_eq!(unitig_graph.unitigs[0].reads.len(), 2);
        assert_eq!(unitig_graph.unitigs[0].sequence, b"ACGTTTTT");
    }

    #[test]
    fn seq_entries_carry_the_unitig_length() {
        let mut graph = chain_graph(&[100, 200], &[10]);
        let unitigs = ma_ug_gen(&mut graph);
        // 90 for the first read after the overlap, 200 for the last one.
        assert_eq!(
            unitigs.graph.sequences[0],
            AsgSeq {
                length: 290,
                deleted: false,
                coverage_class: 0,
            }
        );
    }
}

//! The hybrid unitig graph over the `usg_t` vertex set (`gfa_ut.cpp`
//! 12524): `ma_ug_hybrid_gen` groups the hybrid graph's vertices into
//! unitigs through the two-direction arc walk, then re-derives the arc set
//! from the vertex marks.

use crate::assembly::asg::{Asg, AsgArc, AsgSeq, MaUtg};
use crate::assembly::ul2ul::Usg;

/// `ma_ug_hybrid_gen` (`gfa_ut.cpp` 12524).
pub fn ma_ug_hybrid_gen(g: &Usg) -> UnitigGraphParts {
    let n_vtx = (g.a.len() << 1) as u32;
    let mut mark = vec![0i32; n_vtx as usize];
    let mut unitigs: Vec<MaUtg> = Vec::new();
    let mut v = 0u32;
    while v < n_vtx {
        if g.a[(v >> 1) as usize].deleted || mark[v as usize] != 0 {
            v += 1;
            continue;
        }
        if g.arc_count_of(v) == 0 && g.arc_count_of(v ^ 1) != 0 {
            v += 1;
            continue;
        }
        mark[v as usize] = 1;
        let mut queue: Vec<u64> = Vec::new();
        let mut start = v;
        let mut end = v ^ 1;
        let mut len = 0u32;
        let mut w = v;
        loop {
            if g.arc_count_of(w) != 1 {
                break;
            }
            let x = g.arcs_of(w)[0].v;
            if g.arc_count_of(x ^ 1) != 1 {
                break;
            }
            let l = g.arcs_of(w)[0].ol;
            queue.push(((u64::from(w)) << 32) | u64::from(l));
            end = x ^ 1;
            len += l;
            w = x;
            if x == v {
                break;
            }
        }
        let mut circ = false;
        if start != end ^ 1 || queue.is_empty() {
            let l = g.a[(end >> 1) as usize].len;
            queue.push(((u64::from(end ^ 1)) << 32) | u64::from(l));
            len += l;
        } else {
            start = u32::MAX;
            end = u32::MAX;
            circ = true;
        }
        if !circ {
            let mut x = v;
            loop {
                if g.arc_count_of(x ^ 1) != 1 {
                    break;
                }
                let w = g.arcs_of(x ^ 1)[0].v ^ 1;
                if g.arc_count_of(w) != 1 {
                    break;
                }
                mark[x as usize] = 1;
                mark[(w ^ 1) as usize] = 1;
                let l = g.arcs_of(w)[0].ol;
                queue.insert(0, ((u64::from(w)) << 32) | u64::from(l));
                start = w;
                len += l;
                x = w;
            }
        }
        if start != u32::MAX {
            mark[start as usize] = 1;
            mark[end as usize] = 1;
        }
        unitigs.push(MaUtg {
            length: len,
            circular: circ,
            start,
            end,
            reads: queue,
            sequence: Vec::new(),
        });
        v += 1;
    }
    for slot in mark.iter_mut() {
        *slot = -1;
    }
    for (i, unitig) in unitigs.iter().enumerate() {
        if unitig.circular {
            continue;
        }
        mark[unitig.start as usize] = ((i as i32) << 1) | 0;
        mark[unitig.end as usize] = ((i as i32) << 1) | 1;
    }
    let mut graph = Asg::new();
    for v in 0..n_vtx {
        for arc in g.arcs_of(v) {
            if arc.deleted {
                continue;
            }
            let from = mark[((arc.ul >> 32) ^ 1) as usize];
            let to = mark[arc.v as usize];
            if from >= 0 && to >= 0 {
                let u = (from ^ 1) as usize;
                let mut l = unitigs[u >> 1].length as i32 - arc.ol as i32;
                if l < 0 {
                    l = 1;
                }
                graph.push_arc(AsgArc {
                    ul: ((u as u64) << 32) | l as u64,
                    v: to as u32,
                    overlap: arc.ol & 0x7fff_ffff,
                    deleted: false,
                    ..AsgArc::default()
                });
            }
        }
    }
    for unitig in &unitigs {
        graph.sequences.push(AsgSeq {
            length: unitig.length,
            deleted: false,
            coverage_class: 0,
        });
    }
    graph.cleanup();
    UnitigGraphParts {
        unitigs,
        graph,
        mark,
    }
}

/// The three outputs of [`ma_ug_hybrid_gen`] the callers consume.
pub struct UnitigGraphParts {
    pub unitigs: Vec<MaUtg>,
    pub graph: Asg,
    /// The `mark` array: `(unitig << 1) | side`, `-1` when unmapped.
    pub mark: Vec<i32>,
}

impl Usg {
    /// The live out-degree (`usg_arc_n`).
    pub fn arc_count_of(&self, v: u32) -> usize {
        self.arcs_of(v).iter().filter(|a| !a.deleted).count()
    }
}

/// `merge_unitig_content` (`Overlaps.cpp` 17814): concatenates the member
/// unitigs' read lists into `collection` (reversing per orientation), then
/// `fill_unitig` re-derives the coverage lengths. The Rust `MaUtg.reads`
/// field already carries the packed `(vertex << 32) | offset` words.
pub fn merge_unitig_content(
    collection: &mut MaUtg,
    unitigs: &[MaUtg],
    _read_graph: &Asg,
) {
    if collection.reads.is_empty() {
        return;
    }
    let mut index = 0usize;
    for &packed in &collection.reads {
        let uid = (packed >> 33) as usize;
        index += unitigs[uid].reads.len();
    }
    let mut buffer: Vec<u64> = Vec::with_capacity(index);
    for &packed in &collection.reads {
        let uid = (packed >> 33) as usize;
        let ori = ((packed >> 32) & 1) as usize;
        let query = &unitigs[uid];
        if ori == 1 {
            for entry in query.reads.iter().rev() {
                buffer.push(entry ^ 0x1_0000_0000);
            }
        } else {
            buffer.extend_from_slice(&query.reads);
        }
    }
    let circ = collection.reads.len() == 1 && unitigs[(collection.reads[0] >> 33) as usize].circular;
    collection.reads = buffer;
    let total_len: u64 = collection.reads.iter().map(|b| (b & 0xffff_ffff) as u64).sum();
    collection.length = total_len as u32;
    if !circ {
        collection.start = (collection.reads[0] >> 32) as u32;
        collection.end = ((collection.reads[collection.reads.len() - 1] >> 32) as u32) ^ 1;
    } else {
        collection.start = u32::MAX;
        collection.end = u32::MAX;
        collection.circular = true;
    }
}

/// `merge_hybrid_utg_content` (`gfa_ut.cpp` 16000 region): the hybrid
/// variant; the member ids remap through `ng.a[id].mm` first, and missing
/// inter-unitig arcs land in `edge` for `fill_unitig`.
pub fn merge_hybrid_utg_content(
    collection: &mut MaUtg,
    raw: &[MaUtg],
    rg: &Asg,
    ng: &crate::assembly::ul2ul::Usg,
    edge: &mut Vec<AsgArc>,
) {
    if collection.reads.is_empty() {
        return;
    }
    let mut index = 0usize;
    for packed in collection.reads.iter_mut() {
        let word = *packed;
        let id = (word >> 33) as usize;
        let ori = (word >> 32) & 1;
        let z = ((u64::from(ng.a[id].mm) << 1) | u64::from(ori)) << 32;
        *packed = z | (word & 0xffff_ffff);
        index += ng.a[id].occ as usize;
    }
    let mut buffer: Vec<u64> = Vec::with_capacity(index);
    for i in 0..collection.reads.len() {
        let packed = collection.reads[i];
        let uid = (packed >> 33) as usize;
        let ori = ((packed >> 32) & 1) as usize;
        let q = &raw[uid];
        if ori == 1 {
            for entry in q.reads.iter().rev() {
                buffer.push(entry ^ 0x1_0000_0000);
            }
        } else {
            buffer.extend_from_slice(&q.reads);
        }
        if i > 0 {
            let prev = collection.reads[i - 1];
            let uv = (prev >> 32) as u32;
            let uw = (packed >> 32) as u32;
            if crate::assembly::gfa_clean::get_real_length_one(rg, uv).is_none()
                || !rg.arcs_of(uv).iter().any(|a| !a.deleted && a.v == uw)
            {
                let mut bv = endpoint(raw, uv);
                if uv & 1 != 0 {
                    bv ^= 1;
                }
                let mut bw = endpoint(raw, uw);
                if uw & 1 != 0 {
                    bw ^= 1;
                }
                edge.push(AsgArc {
                    ul: (u64::from(bv) << 32) | u64::from(rg.sequences[(bv >> 1) as usize].length),
                    v: bw,
                    ..AsgArc::default()
                });
                let uv2 = uv ^ 1;
                let uw2 = uw ^ 1;
                let mut bv2 = endpoint(raw, uv2);
                if uv2 & 1 != 0 {
                    bv2 ^= 1;
                }
                let mut bw2 = endpoint(raw, uw2);
                if uw2 & 1 != 0 {
                    bw2 ^= 1;
                }
                edge.push(AsgArc {
                    ul: (u64::from(bv2) << 32) | u64::from(rg.sequences[(bv2 >> 1) as usize].length),
                    v: bw2,
                    ..AsgArc::default()
                });
            }
        }
    }
    if buffer.is_empty() {
        return;
    }
    let circ =
        collection.reads.len() == 1 && raw[(collection.reads[0] >> 33) as usize].circular;
    collection.reads = buffer;
    let total: u64 = collection.reads.iter().map(|b| (b & 0xffff_ffff) as u64).sum();
    collection.length = total as u32;
    if !circ {
        collection.start = (collection.reads[0] >> 32) as u32;
        collection.end = ((collection.reads[collection.reads.len() - 1] >> 32) as u32) ^ 1;
        collection.circular = false;
    } else {
        collection.start = u32::MAX;
        collection.end = u32::MAX;
        collection.circular = true;
    }
}

fn endpoint(unitigs: &[MaUtg], v: u32) -> u32 {
    let unitig = &unitigs[(v >> 1) as usize];
    if v & 1 != 0 {
        (unitig.reads[0] >> 32) as u32
    } else {
        (unitig.reads[unitig.reads.len() - 1] >> 32) as u32
    }
}

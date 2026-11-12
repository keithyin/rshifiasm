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

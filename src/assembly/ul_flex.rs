//! The `flex_asg_t` machinery of `gfa_ut.cpp` 474–2260: the ultra-long
//! containment cutter used by `ul_clean_gfa` when `is_ou` is set.
//!
//! `flex_asg_t` overlays the read graph with a small list of extra arcs
//! (re-inserted bridges) so the containment walk sees both the graph and the
//! pending edges. [`iter_contain_g`] walks a contained-read run, checks
//! whether the run is a tip (via [`detect_tip2`]) and deletes the run when it
//! is. [`asg_arc_cut_contain`] drives it over every contained vertex.
//!
//! Ported: `init_flex_asg_t` (2957), `get_flex_arcs` (476), `flex_arcs0`
//! (474), `iter_flex_asg` (1837), `push_flex_asg_t` (2005),
//! `detect_tip2` (1846), `trans_check` (1947), `append_trans_check` (2085),
//! `iter_contain_g` (2094), `flex_asg_t_cleanup` (2194),
//! `asg_arc_cut_contain` (2212). `append_notrans_e` (2019) is an empty
//! function in the reference and is not ported.

use crate::assembly::asg::{Asg, AsgArc};
use crate::assembly::hit::{ma_hit_to_arc, ReadToUnitig};
use crate::overlap::paf::PafHit;
use crate::reads::ReadStore;

/// `OU_MASK` (`Overlaps.h`).
const OU_MASK: u32 = 0x3fff;

/// `flex_asg_t` (`Overlaps.h` 238): the graph plus its extra arc list.
pub struct FlexAsg {
    pub minimum_overlap: i64,
    pub maximum_hang: i64,
    pub maximum_hang_rate: i64,
    pub need_srt: bool,
    pub gap_fuzz: i64,
    /// `idx[v]`: the head of vertex `v`'s extra-arc list, `u32::MAX` for none.
    pub idx: Vec<u32>,
    /// `pi[i]`: the next entry of extra arc `i`.
    pub pi: Vec<u32>,
    /// `a`: the extra arcs.
    pub a: Vec<AsgArc>,
}

/// `flex_asg_e_retrive_t` (`Overlaps.h` 247): the cursor of `iter_flex_asg`.
#[derive(Clone, Copy, Default)]
pub struct FlexCursor {
    pub i: [u32; 2],
}

/// `init_flex_asg_t` (`gfa_ut.cpp` 2957).
pub fn init_flex_asg_t(
    minimum_overlap: i64,
    maximum_hang: i64,
    maximum_hang_rate: i64,
    gap_fuzz: i64,
    vertex_count: usize,
) -> FlexAsg {
    FlexAsg {
        minimum_overlap,
        maximum_hang,
        maximum_hang_rate,
        need_srt: false,
        gap_fuzz,
        idx: vec![u32::MAX; vertex_count],
        pi: Vec::new(),
        a: Vec::new(),
    }
}

/// `get_flex_arcs` (`gfa_ut.cpp` 476): the number of live out-arcs of `v`
/// (graph arcs first, then extra arcs). The first one is written to `first`
/// when it is `Some`; extra arcs carry the `0x80000000` tag.
pub fn get_flex_arcs(
    graph: &Asg,
    fg: &FlexAsg,
    v: u32,
    mut first: Option<&mut u32>,
) -> u32 {
    let (start, count) = graph.arc_range(v);
    let mut kv = 0u32;
    for i in 0..count {
        if graph.arcs[start + i].deleted {
            continue;
        }
        if let Some(slot) = first.as_deref_mut() {
            if kv < 1 {
                *slot = (start + i) as u32;
            }
        }
        kv += 1;
    }
    let mut i = fg.idx[v as usize];
    while i != u32::MAX {
        if !fg.a[i as usize].deleted {
            if let Some(slot) = first.as_deref_mut() {
                if kv < 1 {
                    *slot = i | 0x8000_0000;
                }
            }
            kv += 1;
        }
        i = fg.pi[i as usize];
    }
    kv
}

/// `flex_arcs0` (`gfa_ut.cpp` 474): resolves an index produced by
/// [`get_flex_arcs`] to the arc itself.
fn flex_arcs0(graph: &Asg, fg: &FlexAsg, id: u32) -> AsgArc {
    if id & 0x8000_0000 == 0 {
        graph.arcs[id as usize]
    } else {
        fg.a[(id & 0x7fff_ffff) as usize]
    }
}

/// `iter_flex_asg` (`gfa_ut.cpp` 1837): the next arc of `v`, walking the
/// graph arcs first and the extra list afterwards.
fn iter_flex_asg(graph: &Asg, fg: &FlexAsg, cursor: &mut FlexCursor, v: u32) -> Option<AsgArc> {
    let (start, count) = graph.arc_range(v);
    if cursor.i[0] < count as u32 {
        let arc = graph.arcs[start + cursor.i[0] as usize];
        cursor.i[0] += 1;
        return Some(arc);
    }
    if cursor.i[0] != u32::MAX {
        cursor.i[0] = u32::MAX;
        cursor.i[1] = fg.idx[v as usize];
    }
    if cursor.i[1] != u32::MAX {
        let entry = cursor.i[1];
        let arc = fg.a[entry as usize];
        cursor.i[1] = fg.pi[entry as usize];
        return Some(arc);
    }
    None
}

/// `push_flex_asg_t` (`gfa_ut.cpp` 2005): puts `z` into the first deleted
/// graph arc slot at its source, or appends it to the extra list.
fn push_flex_asg(graph: &mut Asg, fg: &mut FlexAsg, z: AsgArc) {
    let v = (z.ul >> 32) as u32;
    let (start, count) = graph.arc_range(v);
    for k in 0..count {
        let slot = start + k;
        if graph.arcs[slot].deleted {
            graph.arcs[slot] = z;
            if !fg.need_srt {
                if k > 0 && graph.arcs[slot].ul < graph.arcs[slot - 1].ul {
                    fg.need_srt = true;
                }
                if k + 1 < count && graph.arcs[slot].ul > graph.arcs[slot + 1].ul {
                    fg.need_srt = true;
                }
            }
            return;
        }
    }
    let k = fg.a.len();
    fg.need_srt = true;
    fg.a.push(z);
    fg.pi.push(fg.idx[v as usize]);
    fg.idx[v as usize] = k as u32;
}

/// `detect_tip2` (`gfa_ut.cpp` 1846): whether the contained run rooted at
/// `id` looks like a tip. `true` corresponds to the reference's `1`.
pub fn detect_tip2(
    graph: &mut Asg,
    fg: &FlexAsg,
    id: u32,
    stack: &mut Vec<u64>,
    ou_rat: f32,
) -> bool {
    for direction in 0..2u32 {
        let v = (id << 1) | direction;
        let mut cursor = FlexCursor {
            i: [0, u32::MAX],
        };
        while let Some(sv) = iter_flex_asg(graph, fg, &mut cursor, v) {
            if sv.deleted {
                continue;
            }
            let w = sv.v ^ 1;
            if graph.sequence_visited[w as usize] & 128 != 0 {
                continue;
            }
            if graph.sequence_visited[(w >> 1) as usize] & 2 != 0 {
                continue;
            }

            let mut inner = FlexCursor {
                i: [0, u32::MAX],
            };
            let mut kw = 0u32;
            let mut ou_max = 0u32;
            let mut mm_ou = 0u32;
            while let Some(sw) = iter_flex_asg(graph, fg, &mut inner, w) {
                if sw.deleted {
                    continue;
                }
                if ou_max < u32::from(sw.occurrence) {
                    ou_max = u32::from(sw.occurrence);
                }
                if graph.sequence_visited[(sw.v >> 1) as usize] & 2 != 0 {
                    if mm_ou < u32::from(sw.occurrence) {
                        mm_ou = u32::from(sw.occurrence);
                    }
                    continue;
                }
                kw += 1;
            }
            if kw < 1 {
                return true;
            }
            if ou_rat >= 0.0 && (mm_ou as f32) > (ou_max as f32) * ou_rat {
                return true;
            }
            graph.sequence_visited[w as usize] |= 128;
            stack.push(u64::from(w));
        }
    }
    false
}

/// `trans_check` (`gfa_ut.cpp` 1947): `true` when `z` is not transitive.
fn trans_check(graph: &Asg, fg: &FlexAsg, z: &AsgArc, b: &mut Vec<u64>) -> bool {
    let v = (z.ul >> 32) as u32;
    {
        let mut rv = FlexCursor {
            i: [0, u32::MAX],
        };
        while let Some(sv) = iter_flex_asg(graph, fg, &mut rv, v) {
            if sv.deleted {
                continue;
            }
            if sv.v == z.v {
                return false;
            }
            let w = sv.v;
            let mut rw = FlexCursor {
                i: [0, u32::MAX],
            };
            while let Some(sw) = iter_flex_asg(graph, fg, &mut rw, w) {
                if sw.deleted {
                    continue;
                }
                if sw.v == z.v {
                    return false;
                }
            }
        }
    }

    let bn = b.len();
    let mut v = v;
    let mut vl: u64 = 0;
    let l = u64::from(z.node_length()) + fg.gap_fuzz as u64;
    b.push(u64::from(v));
    while b.len() > bn {
        let packed = b.pop().unwrap();
        v = packed as u32;
        vl = packed >> 32;
        let mut rv = FlexCursor {
            i: [0, u32::MAX],
        };
        while let Some(sv) = iter_flex_asg(graph, fg, &mut rv, v) {
            if sv.deleted {
                continue;
            }
            let d = vl + u64::from(sv.node_length());
            if d > l {
                continue;
            }
            if sv.v == z.v {
                b.truncate(bn);
                return false;
            }
            b.push((d << 32) | u64::from(sv.v));
        }
    }
    b.truncate(bn);
    true
}

/// One `get_specific_overlap` + `ma_hit2arc` pair with `ou` taken from `bl`.
fn specific_overlap_arc(
    fg: &FlexAsg,
    reads: &ReadStore,
    list: &[PafHit],
    qn: u32,
    tn: u32,
) -> Option<AsgArc> {
    let index = list
        .iter()
        .position(|hit| hit.query_identifier == qn && hit.target_identifier == tn)?;
    if list[index].deleted {
        return None;
    }
    let hit = &list[index];
    let (result, mut arc) = ma_hit_to_arc(
        hit,
        reads.read_length(qn as usize) as i32,
        reads.read_length(tn as usize) as i32,
        fg.maximum_hang as i32,
        fg.maximum_hang_rate as f32,
        fg.minimum_overlap as i32,
    );
    if result < 0 {
        return None;
    }
    arc.occurrence = list[index].block_length.min(OU_MASK) as u16;
    Some(arc)
}

/// `append_trans_check` (`gfa_ut.cpp` 2085).
fn append_trans_check(
    fg: &FlexAsg,
    sources: &[Vec<PafHit>],
    reads: &ReadStore,
    a: &[u64],
) -> bool {
    for (i, &v) in a.iter().enumerate() {
        let v = v as u32;
        for &w_raw in &a[i + 1..] {
            let w = (w_raw as u32) ^ 1;
            if specific_overlap_arc(fg, reads, &sources[(v >> 1) as usize], v, w).is_none() {
                return false;
            }
            if specific_overlap_arc(fg, reads, &sources[(w >> 1) as usize], w ^ 1, v ^ 1).is_none() {
                return false;
            }
        }
    }
    true
}

/// `iter_contain_g` (`gfa_ut.cpp` 2094): walks the contained run hanging off
/// `v0`, deleting it when [`detect_tip2`] flags it. Returns the number of
/// deleted reads.
pub fn iter_contain_g(
    graph: &mut Asg,
    fg: &FlexAsg,
    sources: &[Vec<PafHit>],
    reads: &ReadStore,
    read_to_unitig: &ReadToUnitig,
    v0: u32,
    b: &mut Vec<u64>,
    st: &mut Vec<u64>,
    ou_rat: f32,
    only_trans_nn: bool,
) -> u32 {
    let mut w: u32 = u32::MAX;
    let mut cnt = 0u32;
    let mut ulen;
    b.clear();
    st.clear();
    st.push(u64::from(v0));
    while !st.is_empty() {
        let mut v = st.pop().unwrap() as u32;
        let mut kv: u32;
        let mut kw: u32;
        ulen = 0u32;
        if graph.sequence_visited[(v >> 1) as usize] != 0 {
            continue;
        }
        while graph.sequence_visited[(v >> 1) as usize] == 0
            && is_contain_r(read_to_unitig, v >> 1)
        {
            // push a unitig
            let mut first = 0u32;
            kv = get_flex_arcs(graph, fg, v, Some(&mut first));
            ulen += 1;
            b.push(u64::from(v));
            if kv == 1 {
                let s = flex_arcs0(graph, fg, first);
                w = s.v;
                kw = get_flex_arcs(graph, fg, w ^ 1, None);
                if kw == 1 {
                    v = w;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        if ulen == 0
            || graph.sequence_visited[(v >> 1) as usize] != 0
            || !is_contain_r(read_to_unitig, v >> 1)
        {
            for i in b.len() - ulen as usize..b.len() {
                graph.sequence_visited[(b[i] as u32 >> 1) as usize] = 1;
                b[i] |= 0x1_0000_0000;
            }
            continue;
        }

        for i in 0..b.len() {
            if b[i] & 0x1_0000_0000 != 0 {
                continue;
            }
            graph.sequence_visited[((b[i] as u32) >> 1) as usize] = 2;
        }
        let st_n = st.len();
        let mut found = false;
        for i in 0..b.len() {
            if b[i] & 0x1_0000_0000 != 0 {
                continue;
            }
            if detect_tip2(graph, fg, (b[i] as u32) >> 1, st, ou_rat) {
                found = true;
                break;
            }
        }
        let mut is_purge = !found;

        // the newly pushed neighbours: strip the 128 marker
        let mut m = st_n;
        for i in st_n..st.len() {
            let id = st[i] as usize;
            if graph.sequence_visited[id] & 128 != 0 {
                graph.sequence_visited[id] -= 128;
            }
            st[m] = st[i];
            m += 1;
        }
        st.truncate(m);

        if is_purge && only_trans_nn {
            is_purge = append_trans_check(fg, sources, reads, &st[st_n..]);
        }

        if is_purge {
            for i in 0..b.len() {
                let x = ((b[i] as u32) >> 1) as usize;
                if graph.sequence_visited[x] & 2 != 0 {
                    graph.delete_sequence(x as u32);
                    cnt += 1;
                }
                graph.sequence_visited[x] = 0;
            }
            st.truncate(st_n);
            return cnt;
        }

        for i in 0..b.len() {
            if b[i] & 0x1_0000_0000 != 0 {
                continue;
            }
            graph.sequence_visited[((b[i] as u32) >> 1) as usize] = 1;
        }
        st.truncate(st_n);

        let mut rr = FlexCursor {
            i: [0, u32::MAX],
        };
        loop {
            let Some(s) = iter_flex_asg(graph, fg, &mut rr, v) else {
                break;
            };
            if s.deleted {
                continue;
            }
            if graph.sequence_visited[(s.v >> 1) as usize] != 0 {
                continue;
            }
            if !is_contain_r(read_to_unitig, s.v >> 1) {
                continue;
            }
            st.push(u64::from(s.v));
        }
    }
    for i in 0..b.len() {
        graph.sequence_visited[((b[i] as u32) >> 1) as usize] = 0;
    }
    cnt
}

/// `flex_asg_t_cleanup` (`gfa_ut.cpp` 2194): folds the surviving extra arcs
/// back into the graph.
pub fn flex_asg_t_cleanup(graph: &mut Asg, fg: &mut FlexAsg) {
    if !fg.a.is_empty() {
        for arc in &fg.a {
            if arc.deleted {
                continue;
            }
            graph.push_arc(*arc);
        }
        graph.index.clear();
        graph.is_sorted = false;
    } else if fg.need_srt {
        graph.is_sorted = false;
    }
    graph.cleanup();
}

/// `asg_arc_cut_contain` (`gfa_ut.cpp` 2212): deletes the contained runs that
/// look like tips.
#[allow(clippy::too_many_arguments)]
pub fn asg_arc_cut_contain(
    graph: &mut Asg,
    fg: &mut FlexAsg,
    sources: &[Vec<PafHit>],
    reads: &ReadStore,
    read_to_unitig: &ReadToUnitig,
    buffer: &mut Vec<u64>,
    buffer0: &mut Vec<u64>,
    ou_rat: f32,
    only_trans_nn: bool,
) -> u32 {
    let mut w: u32 = u32::MAX;
    let mut cnt = 0u32;
    let vertex_count = graph.sequences.len() * 2;
    buffer.clear();
    buffer0.clear();
    fg.idx.iter_mut().for_each(|slot| *slot = u32::MAX);
    fg.need_srt = false;
    for slot in graph.sequence_visited.iter_mut() {
        *slot = 0;
    }

    for v in 0..vertex_count as u32 {
        if graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if !is_contain_r(read_to_unitig, v >> 1) {
            continue;
        }
        let mut first = 0u32;
        if get_flex_arcs(graph, fg, v ^ 1, Some(&mut first)) == 1 {
            let s = flex_arcs0(graph, fg, first);
            if get_flex_arcs(graph, fg, s.v ^ 1, None) == 1 {
                continue;
            }
        }
        cnt += iter_contain_g(
            graph,
            fg,
            sources,
            reads,
            read_to_unitig,
            v,
            buffer,
            buffer0,
            ou_rat,
            only_trans_nn,
        );
    }

    if cnt > 0 {
        flex_asg_t_cleanup(graph, fg);
    }
    cnt
}

/// `is_contain_r` (`gfa_ut.h` 6): the read maps to a unitig *read* slot.
fn is_contain_r(read_to_unitig: &ReadToUnitig, read: u32) -> bool {
    if read as usize >= read_to_unitig.len() {
        return false;
    }
    let (unitig, is_unitig) = read_to_unitig.get(read);
    unitig != u32::MAX && is_unitig == 0
}

/// The unused `w` binding in [`asg_arc_cut_contain`] mirrors the reference,
/// which also never reads it.
#[allow(dead_code)]
fn _unused(w: &mut u32) {
    let _ = w;
    *w = u32::MAX;
}

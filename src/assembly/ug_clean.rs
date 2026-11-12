//! `hic_clean_adv` (`Overlaps.cpp` 14403) and the `ug_clean_t` cluster of
//! `deep_graph_clean` (`gfa_ut.cpp` 17697–18202): the advanced bubble-based
//! graph polishing that runs between the renewal rounds and inside every
//! `ul_realignment_gfa` call.
//!
//! The bubble primitives (`asg_bub_pop1_primary_trio`, `get_unitig`,
//! `get_bubbles`, `IF_HOM`) and the bubble-chain builder (`gen_bubble_chain`)
//! live in their own modules; this module only wires them the way the
//! reference does.

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::hit::{ALTER_LABEL, PRIMARY_LABEL};
use crate::assembly::bubble::BubbleBuffer;
use crate::assembly::bubble_chain::{get_bubbles, BubbleType};
use crate::assembly::ul_renew_gfa::if_hom;
use crate::assembly::bubble_identify::asg_bub_pop1_primary_trio;
use crate::assembly::gfa_clean::{get_real_length, get_real_length_one};
use crate::assembly::hit::{CoverageCut, ReadToUnitig};
use crate::assembly::trio_clean::{DOUBLE_CHECK_THRES, TRIO_THRES};
use crate::assembly::unitig::{
    ma_ug_gen, ma_ug_gen_primary, AMBIGU, FATHER, MOTHER,
};
use crate::assembly::ul_renew_gfa::gen_bubble_chain;
use crate::command_line::options::AssemblyOptions;
use crate::correction::ec_driver::PafList;
use crate::container::ksort::radix_sort_u64;

/// One `asg_bub_pop1_primary_trio` probe with the reference's constant
/// `(-1, -1)` vertex flags and no coverage/trio payloads.
macro_rules! probe_bubble {
    ($graph:expr, $v:expr, $tlen:expr, $buf:expr, $pop:expr) => {
        asg_bub_pop1_primary_trio(
            $graph,
            None,
            $v,
            $tlen,
            $buf,
            u32::MAX,
            u32::MAX,
            $pop,
            None,
            &[],
            false,
            false,
        )
    };
}

/// `update_ug_ou` (`Overlaps.cpp` 8694): copies every unitig arc's `ou` from
/// the string-graph arc that realizes it.
pub fn update_ug_ou(unitig_graph: &mut crate::assembly::unitig::UnitigGraph, sg: &Asg) {
    for k in 0..unitig_graph.graph.arcs.len() {
        let ue = unitig_graph.graph.arcs[k];
        let uv = (ue.ul >> 32) as u32;
        let uw = ue.v;
        let ru = &unitig_graph.unitigs[(uv >> 1) as usize];
        let rv = if uv & 1 != 0 {
            ru.start ^ 1
        } else {
            ru.end ^ 1
        };
        let rwu = &unitig_graph.unitigs[(uw >> 1) as usize];
        let rw = if uw & 1 != 0 {
            rwu.end
        } else {
            rwu.start
        };
        let mut found = 0u32;
        let mut ou = 0u16;
        for arc in sg.arcs_of(rv) {
            if arc.v == rw {
                found = 1;
                ou = arc.occurrence;
                break;
            }
        }
        assert!(found == 1);
        unitig_graph.graph.arcs[k].occurrence = ou;
    }
}

/// `get_ug_tip_trio_infor` (`gfa_ut.cpp` 1799) over the unitig graph.
fn get_ug_tip_trio_infor(
    unitig_graph: &crate::assembly::unitig::UnitigGraph,
    beg_node: u32,
    trio_labels: &[u8],
) -> u8 {
    let mut v = beg_node;
    let mut node_len = 0u32;
    let mut father_occ = 0u32;
    let mut mother_occ = 0u32;
    loop {
        let kv = get_real_length(&unitig_graph.graph, v);
        let unitig = &unitig_graph.unitigs[(v >> 1) as usize];
        node_len += unitig.reads.len() as u32;
        for &packed in &unitig.reads {
            let s = (packed >> 33) as usize;
            match trio_labels[s] {
                x if x == FATHER => father_occ += 1,
                x if x == MOTHER => mother_occ += 1,
                _ => {}
            }
        }
        if kv != 1 {
            break;
        }
        let w = get_real_length_one(&unitig_graph.graph, v).expect("kv == 1");
        if get_real_length(&unitig_graph.graph, w ^ 1) != 1 {
            break;
        }
        v = w;
        if v == beg_node {
            break;
        }
    }
    let u_len = node_len;
    let e_len = father_occ + mother_occ;
    if e_len == 0 {
        return AMBIGU;
    }
    if father_occ >= mother_occ {
        if (father_occ as f32) > TRIO_THRES * e_len as f32
            && (father_occ as f64) >= DOUBLE_CHECK_THRES * u_len as f64
        {
            return FATHER;
        }
    } else if (mother_occ as f32) > TRIO_THRES * e_len as f32
        && (mother_occ as f64) >= DOUBLE_CHECK_THRES * u_len as f64
    {
        return MOTHER;
    }
    AMBIGU
}

/// The `ug_clean_t` workspace of `deep_graph_clean`.
pub struct UgClean<'a> {
    /// `idx`: the best-het target of every vertex (`(u32)-1` when none).
    pub idx: Vec<u32>,
    /// `bid`: the bubble membership flags of every vertex.
    pub bid: Vec<u32>,
    /// `tlen`: the bubble-pop distance threshold.
    pub tlen: u64,
    /// The `buf_t` bubble workspace.
    pub buffer: BubbleBuffer,
    /// `ug`: the unitig graph under cleaning (a separate `ma_ug_gen`).
    pub unitig_graph: crate::assembly::unitig::UnitigGraph,
    /// `sg`: the read graph the unitig reads are deleted from.
    pub sg: &'a mut Asg,
    pub is_ou: bool,
    pub is_trio: bool,
    pub max_ext: i64,
    pub len_rat: f64,
    pub ou_rat: f64,
    pub min_ou: i64,
    pub trio_labels: Vec<u8>,
}

impl<'a> UgClean<'a> {
    /// `init_ug_clean_t` (`gfa_ut.cpp` 17697).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        sg: &'a mut Asg,
        trio_labels: Vec<u8>,
        is_ou: bool,
        is_trio: bool,
        max_ext: i64,
        len_rat: f64,
        ou_rat: f64,
        min_ou: i64,
    ) -> Self {
        let mut unitig_graph = ma_ug_gen(sg);
        let n_seq = unitig_graph.graph.sequences.len();
        for seq in unitig_graph.graph.sequences.iter_mut() {
            if seq.deleted {
                continue;
            }
            seq.coverage_class = PRIMARY_LABEL;
        }
        if is_ou {
            update_ug_ou(&mut unitig_graph, sg);
        }
        UgClean {
            idx: vec![u32::MAX; n_seq << 1],
            bid: vec![u32::MAX; n_seq << 1],
            tlen: 0,
            buffer: BubbleBuffer {
                info: vec![Default::default(); n_seq << 1],
                ..BubbleBuffer::default()
            },
            unitig_graph,
            sg,
            is_ou,
            is_trio,
            max_ext,
            len_rat,
            ou_rat,
            min_ou,
            trio_labels,
        }
    }
}

/// `update_ug_clean_t` (`gfa_ut.cpp` 17726): recompute `tlen`, `bid` and the
/// beg/sink index of every bubble.
pub fn update_ug_clean(sl: &mut UgClean<'_>) {
    let n_vtx = (sl.unitig_graph.graph.sequences.len() << 1) as u32;
    sl.tlen =
        crate::assembly::ul_renew_gfa::get_bub_pop_max_dist_advance(&mut sl.unitig_graph.graph);
    for k in 0..sl.unitig_graph.graph.sequences.len() {
        sl.idx[k << 1] = u32::MAX;
        sl.idx[(k << 1) + 1] = u32::MAX;
        sl.bid[k << 1] = u32::MAX;
        sl.bid[(k << 1) + 1] = u32::MAX;
        if sl.unitig_graph.graph.sequences[k].deleted {
            continue;
        }
        sl.unitig_graph.graph.sequences[k].coverage_class = PRIMARY_LABEL;
    }

    for v in 0..n_vtx {
        if sl.unitig_graph.graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if sl.unitig_graph.graph.arc_count(v) < 2 {
            continue;
        }
        if sl.bid[v as usize] != u32::MAX {
            continue;
        }
        if probe_bubble!(&mut sl.unitig_graph.graph, v, sl.tlen, &mut sl.buffer, false) != 0 {
            let sink = sl.buffer.sources[0];
            for &node in &sl.buffer.visited {
                if node == v || node == sink {
                    continue;
                }
                sl.bid[node as usize] = 1;
                sl.bid[(node ^ 1) as usize] = 1;
            }
            sl.bid[v as usize] = 2;
            sl.bid[(sink ^ 1) as usize] = 3;
        }
    }

    for v in 0..n_vtx {
        if sl.bid[v as usize] != 2 {
            continue;
        }
        if probe_bubble!(&mut sl.unitig_graph.graph, v, sl.tlen, &mut sl.buffer, false) != 0 {
            let sink = sl.buffer.sources[0];
            for &node in &sl.buffer.visited {
                if node == v || node == sink {
                    continue;
                }
                sl.idx[node as usize] = v;
                sl.idx[(node ^ 1) as usize] = sink ^ 1;
            }
            sl.idx[v as usize] = v;
            sl.idx[(sink ^ 1) as usize] = sink ^ 1;
        }
    }
    let bid = sl.bid.clone();
    sl.bid.copy_from_slice(&bid);
}

/// `get_best_het_node` (`gfa_ut.cpp` 17732): the best heterozygous successor
/// of `v`, or `(u32)-1`.
#[allow(clippy::too_many_arguments)]
fn get_best_het_node(sl: &UgClean<'_>, v: u32) -> u32 {
    let graph = &sl.unitig_graph.graph;
    if graph.sequences[(v >> 1) as usize].deleted {
        return u32::MAX;
    }
    let (start, count) = graph.arc_range(v);
    let arcs = &graph.arcs[start..start + count];
    let mut w = u32::MAX;
    let mut kv = 0u32;
    for arc in arcs {
        if arc.deleted {
            continue;
        }
        // Could not connect to the beg/sink node.
        if (arc.v >> 1) == (sl.bid[v as usize] >> 1)
            || (arc.v >> 1) == (sl.bid[(v ^ 1) as usize] >> 1)
        {
            return u32::MAX;
        }
        kv += 1;
        w = arc.v;
    }
    if kv < 2 {
        return w;
    }
    let mut trio_f = u32::MAX;
    let mut ntrio_f = u32::MAX;
    if sl.is_trio {
        let f = get_ug_tip_trio_infor(&sl.unitig_graph, v ^ 1, &sl.trio_labels);
        trio_f = f as u32;
        ntrio_f = if trio_f == FATHER as u32 {
            MOTHER as u32
        } else if trio_f == MOTHER as u32 {
            FATHER as u32
        } else {
            u32::MAX
        };
    }

    let mut ol_max = 0u32;
    let mut ol_k = usize::MAX;
    for (k, arc) in arcs.iter().enumerate() {
        if arc.deleted {
            continue;
        }
        if ol_max < arc.overlap {
            ol_max = arc.overlap;
            ol_k = k;
        }
    }
    if ol_k == usize::MAX {
        return u32::MAX;
    }
    for (k, arc) in arcs.iter().enumerate() {
        if arc.deleted || k == ol_k {
            continue;
        }
        if f64::from(arc.overlap) > f64::from(arcs[ol_k].overlap) * sl.len_rat {
            return u32::MAX;
        }
        if sl.is_ou
            && (u64::from(arc.occurrence) > sl.min_ou as u64)
            && (f64::from(arc.occurrence) > f64::from(arcs[ol_k].occurrence) * sl.ou_rat)
        {
            return u32::MAX;
        }
    }
    if sl.is_trio
        && get_ug_tip_trio_infor(&sl.unitig_graph, arcs[ol_k].v, &sl.trio_labels) as u32 == ntrio_f
    {
        return u32::MAX;
    }
    arcs[ol_k].v
}

#[inline]
fn is_best_arc(sl: &UgClean<'_>, v: u32, w: u32) -> bool {
    sl.idx[v as usize] == w && sl.idx[(w ^ 1) as usize] == (v ^ 1)
}

/// `cal_bub_best` (`gfa_ut.cpp` 17806) for unitig `i`.
fn cal_bub_best(sl: &mut UgClean<'_>, i: usize) {
    sl.idx[i << 1] = u32::MAX;
    sl.idx[(i << 1) + 1] = u32::MAX;
    if sl.unitig_graph.graph.sequences[i].deleted {
        return;
    }
    if sl.bid[i << 1] == u32::MAX || sl.bid[(i << 1) + 1] == u32::MAX {
        return; // not within a bubble
    }

    let mut v = (i << 1) as u32;
    if sl.bid[v as usize] != v {
        let w = get_best_het_node(sl, v);
        if w != u32::MAX && (v ^ 1) == get_best_het_node(sl, w ^ 1) {
            sl.idx[v as usize] = w;
        }
    }

    v = (i << 1) as u32 + 1;
    if sl.bid[v as usize] != v {
        let w = get_best_het_node(sl, v);
        if w != u32::MAX && (v ^ 1) == get_best_het_node(sl, w ^ 1) {
            sl.idx[v as usize] = w;
        }
    }
}

/// `usg_topocut_aux_unambi1` (`gfa_ut.cpp` 17826): the single live successor
/// of `v`, or `(u32)-1`.
fn usg_topocut_aux_unambi1(graph: &Asg, v: u32) -> u32 {
    let (start, count) = graph.arc_range(v);
    let arcs = &graph.arcs[start..start + count];
    let mut k = count;
    let mut kv = 0u32;
    for (i, arc) in arcs.iter().enumerate() {
        if !arc.deleted {
            kv += 1;
            k = i;
        }
    }
    if kv != 1 {
        return u32::MAX;
    }
    arcs[k].v
}

/// `usg_topocut_aux_del` (`gfa_ut.cpp` 17837): walks the unambiguous tip and
/// records every visited unitig vertex in `b`.
fn usg_topocut_aux_del(sl: &UgClean<'_>, mut v: u32, b: &mut Vec<u64>) -> i64 {
    let mut n_ext = 0i64;
    while n_ext < sl.max_ext && v != u32::MAX {
        if usg_topocut_aux_unambi1(&sl.unitig_graph.graph, v ^ 1) == u32::MAX {
            break;
        }
        b.push(u64::from(v));
        n_ext += sl.unitig_graph.unitigs[(v >> 1) as usize].reads.len() as i64;
        v = usg_topocut_aux_unambi1(&sl.unitig_graph.graph, v);
    }
    n_ext
}

/// `asg_arc_del_by_ug` (`gfa_ut.cpp` 17855): mark the string-graph arc that
/// realizes the unitig arc `uv -> uw` deleted.
fn asg_arc_del_by_ug(unitig_graph: &crate::assembly::unitig::UnitigGraph, sg: &mut Asg, uv: u32, uw: u32) {
    let ru = &unitig_graph.unitigs[(uv >> 1) as usize];
    let rv = if uv & 1 != 0 { ru.start ^ 1 } else { ru.end ^ 1 };
    let rwu = &unitig_graph.unitigs[(uw >> 1) as usize];
    let rw = if uw & 1 != 0 { rwu.end } else { rwu.start };
    let (start, count) = sg.arc_range(rv);
    for i in 0..count {
        if sg.arcs[start + i].v == rw {
            sg.arcs[start + i].deleted = true;
            return;
        }
    }
    panic!("asg_arc_del_by_ug: the unitig arc has no string-graph realization");
}

/// `cal_utg_occ` (`gfa_ut.cpp` 17870): the read count of the unitary chain
/// that starts at `beg_node`.
fn cal_utg_occ(unitig_graph: &crate::assembly::unitig::UnitigGraph, beg_node: u32) -> u32 {
    let mut v = beg_node;
    let mut occ = 0u32;
    loop {
        let kv = get_real_length(&unitig_graph.graph, v);
        occ += unitig_graph.unitigs[(v >> 1) as usize].reads.len() as u32;
        if kv != 1 {
            break;
        }
        let w = get_real_length_one(&unitig_graph.graph, v).expect("kv == 1");
        if get_real_length(&unitig_graph.graph, w ^ 1) != 1 {
            break;
        }
        v = w;
        if v == beg_node {
            break;
        }
    }
    occ
}

/// `cal_bub_best_by_len` (`gfa_ut.cpp` 17891): cut the longer alternate arc of
/// every best-arc-marked bubble whose tip does not extend beyond `max_ext`.
#[allow(clippy::too_many_arguments)]
pub fn cal_bub_best_by_len(
    sl: &mut UgClean<'_>,
    input: Option<&mut Vec<u64>>,
    max_ext: u32,
    is_trio: bool,
    is_ou: bool,
    len_rat: f64,
    ou_rat: f64,
    min_ou: u32,
    min_node: u32,
) {
    let mut owned: Vec<u64> = Vec::new();
    let b: &mut Vec<u64> = match input {
        Some(buffer) => buffer,
        None => &mut owned,
    };
    b.clear();
    let mut cnt = 0usize;

    let n_vtx = (sl.unitig_graph.graph.sequences.len() << 1) as u32;
    for v in 0..n_vtx {
        if sl.unitig_graph.graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        let (start, count) = sl.unitig_graph.graph.arc_range(v);
        if count < 2 {
            continue;
        }
        let mut kv = 0u32;
        let mut bb = 0u32;
        for i in 0..count {
            let arc = sl.unitig_graph.graph.arcs[start + i];
            if arc.deleted {
                continue;
            }
            if is_best_arc(sl, v, arc.v) {
                bb += 1;
            }
            kv += 1;
        }
        if kv < 2 || bb == 0 || kv <= bb {
            continue; // it is impossible that kv <= bb
        }

        for i in 0..count {
            let arc = sl.unitig_graph.graph.arcs[start + i];
            if arc.deleted {
                continue;
            }
            if is_best_arc(sl, v, arc.v) {
                continue;
            }
            b.push((u64::from(arc.overlap) << 32) | u64::from((start + i) as u32));
        }
    }

    radix_sort_u64(b);
    for k in 0..b.len() {
        let arc_index = b[k] as u32 as usize;
        if sl.unitig_graph.graph.arcs[arc_index].deleted {
            continue;
        }

        let v = (sl.unitig_graph.graph.arcs[arc_index].ul >> 32) as u32;
        let w = sl.unitig_graph.graph.arcs[arc_index].v ^ 1;
        if sl.unitig_graph.graph.sequences[(v >> 1) as usize].deleted
            || sl.unitig_graph.graph.sequences[(w >> 1) as usize].deleted
        {
            continue;
        }
        let (v_start, nv) = sl.unitig_graph.graph.arc_range(v);
        let (w_start, nw) = sl.unitig_graph.graph.arc_range(w);
        if nv <= 1 && nw <= 1 {
            continue;
        }

        let mut trio_f = u32::MAX;
        let mut ntrio_f = u32::MAX;
        if is_trio {
            if live_arc_count(&sl.unitig_graph.graph, v) <= 1
                && live_arc_count(&sl.unitig_graph.graph, w) <= 1
            {
                continue; // speedup
            }
            let f = get_ug_tip_trio_infor(&sl.unitig_graph, v ^ 1, &sl.trio_labels);
            trio_f = f as u32;
            ntrio_f = if trio_f == FATHER as u32 {
                MOTHER as u32
            } else if trio_f == MOTHER as u32 {
                FATHER as u32
            } else {
                u32::MAX
            };
        }

        // The paired reverse arc of `v -> w` at `w`.
        let ve_pair = sl.unitig_graph.graph.arcs[arc_index];
        let mut we_index = usize::MAX;
        for i in 0..nw {
            if sl.unitig_graph.graph.arcs[w_start + i].v == v ^ 1 {
                we_index = w_start + i;
                break;
            }
        }
        // The reference dereferences `we` unconditionally.
        assert!(we_index != usize::MAX);
        let mm_ol = ve_pair.overlap.min(sl.unitig_graph.graph.arcs[we_index].overlap);
        let mm_ou = ve_pair
            .occurrence
            .min(sl.unitig_graph.graph.arcs[we_index].occurrence);
        let mut lnid = 0u32;

        let mut kv = 0u32;
        let mut ol_max = 0u32;
        let mut ou_max = 0u16;
        for i in 0..nv {
            let arc = sl.unitig_graph.graph.arcs[v_start + i];
            if arc.deleted {
                continue;
            }
            kv += 1;
            if is_trio
                && get_ug_tip_trio_infor(&sl.unitig_graph, arc.v, &sl.trio_labels) as u32 == ntrio_f
            {
                continue;
            }
            if ol_max < arc.overlap {
                ol_max = arc.overlap;
            }
            if ou_max < arc.occurrence {
                ou_max = arc.occurrence;
            }
            if lnid == 0 && is_best_arc(sl, v, arc.v) {
                if cal_utg_occ(&sl.unitig_graph, v ^ 1) >= min_node
                    || cal_utg_occ(&sl.unitig_graph, arc.v) >= min_node
                {
                    lnid = 1;
                }
            }
        }
        if kv < 1 {
            continue;
        }
        if kv >= 2 {
            if f64::from(mm_ol) > f64::from(ol_max) * len_rat {
                continue;
            }
            if is_ou && (u64::from(mm_ou) > min_ou as u64) && (f64::from(mm_ou) > f64::from(ou_max) * ou_rat) {
                continue;
            }
        }

        let mut kw = 0u32;
        let mut ol_max = 0u32;
        let mut ou_max = 0u16;
        for i in 0..nw {
            let arc = sl.unitig_graph.graph.arcs[w_start + i];
            if arc.deleted {
                continue;
            }
            kw += 1;
            if is_trio
                && get_ug_tip_trio_infor(&sl.unitig_graph, arc.v, &sl.trio_labels) as u32 == ntrio_f
            {
                continue;
            }
            if ol_max < arc.overlap {
                ol_max = arc.overlap;
            }
            if ou_max < arc.occurrence {
                ou_max = arc.occurrence;
            }
            if lnid == 0 && is_best_arc(sl, w, arc.v) {
                if cal_utg_occ(&sl.unitig_graph, w ^ 1) >= min_node
                    || cal_utg_occ(&sl.unitig_graph, arc.v) >= min_node
                {
                    lnid = 1;
                }
            }
        }
        if kw < 1 {
            continue;
        }
        if kw >= 2 {
            if f64::from(mm_ol) > f64::from(ol_max) * len_rat {
                continue;
            }
            if is_ou && (u64::from(mm_ou) > min_ou as u64) && (f64::from(mm_ou) > f64::from(ou_max) * ou_rat) {
                continue;
            }
        }

        if kv <= 1 && kw <= 1 {
            continue;
        }
        if lnid == 0 {
            continue;
        }

        let mut to_del = false;
        let mut del_v = u32::MAX;
        let mut del_w = u32::MAX;
        if kv > 1 && kw > 1 {
            to_del = true;
        } else if kw == 1 {
            let tip = usg_topocut_aux_del(sl, w ^ 1, &mut Vec::new());
            if tip < i64::from(max_ext) {
                to_del = true;
                del_w = w ^ 1;
            }
        } else if kv == 1 {
            let tip = usg_topocut_aux_del(sl, v ^ 1, &mut Vec::new());
            if tip < i64::from(max_ext) {
                to_del = true;
                del_v = v ^ 1;
            }
        }

        if to_del {
            let bn = b.len();
            if del_v != u32::MAX {
                usg_topocut_aux_del(sl, del_v, b);
            }
            if del_w != u32::MAX {
                usg_topocut_aux_del(sl, del_w, b);
            }
            for i in bn..b.len() {
                let unitig_id = (b[i] >> 1) as usize;
                if sl.unitig_graph.unitigs[unitig_id].reads.is_empty() {
                    continue;
                }
                for &packed in &sl.unitig_graph.unitigs[unitig_id].reads {
                    sl.sg.delete_sequence((packed >> 33) as u32);
                }
                sl.unitig_graph.graph.delete_sequence(unitig_id as u32);
                sl.unitig_graph.unitigs[unitig_id].reads.clear();
            }
            b.truncate(bn);

            sl.unitig_graph.graph.arcs[arc_index].deleted = true;
            sl.unitig_graph.graph.arcs[we_index].deleted = true;
            let ve = sl.unitig_graph.graph.arcs[arc_index];
            let we = sl.unitig_graph.graph.arcs[we_index];
            asg_arc_del_by_ug(&sl.unitig_graph, sl.sg, (ve.ul >> 32) as u32, ve.v);
            asg_arc_del_by_ug(&sl.unitig_graph, sl.sg, (we.ul >> 32) as u32, we.v);
            cnt += 1;
        }
    }

    if cnt > 0 {
        sl.unitig_graph.graph.cleanup();
        sl.sg.cleanup();
    }
}

/// `get_arcs(g, v, NULL, 0)` — the live out-degree.
fn live_arc_count(graph: &Asg, v: u32) -> u32 {
    get_real_length(graph, v)
}

/// `hic_clean_adv` (`Overlaps.cpp` 14403): the advanced bubble polishing of
/// the string graph.
#[allow(clippy::too_many_arguments)]
/// `hic_clean_adv` (`Overlaps.cpp` 14403): the advanced bubble polishing of
/// the string graph.
#[allow(clippy::too_many_arguments)]
pub fn hic_clean_adv(
    options: &mut AssemblyOptions,
    sg: &mut Asg,
    sources: &[PafList],
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    coverage_cut: &mut [CoverageCut],
    bubble_labels: Option<&mut crate::assembly::bubble::BubbleLabel>,
) {
    let mut unitig_graph = ma_ug_gen_primary(sg, PRIMARY_LABEL);
    let n_vtx = unitig_graph.graph.sequences.len() << 1;
    let bub_rate = 0.1;
    let avoid_het = options.polyploidy > 2;
    // `gen_bubble_chain(sg, ug, uopt, &bf, ...)`: the het flags come back at
    // read granularity and are immediately replaced by the unitig-level
    // buffer (`REALLOC(bf, n_vtx); memset(bf, 0, ...)`).
    let (bub, _bf_het) = gen_bubble_chain(
        sg,
        &mut unitig_graph,
        sources,
        reverse_sources,
        read_to_unitig,
        options,
        coverage_cut,
        bubble_labels,
        avoid_het,
    );
    let mut bf = vec![0u8; n_vtx];
    let mut b = BubbleBuffer {
        info: vec![Default::default(); n_vtx],
        ..BubbleBuffer::default()
    };
    let mut buf: Vec<u64> = Vec::new();
    let mut n0 = 0u32;
    let mut n1 = 0u32;

    for seq in unitig_graph.graph.sequences.iter_mut() {
        if seq.deleted {
            continue;
        }
        seq.coverage_class = PRIMARY_LABEL;
    }

    let empty_chains: Vec<MaUtg> = Vec::new();
    let chains = bub.bubble_unitigs.as_deref().unwrap_or(&empty_chains);
    for chain in chains.iter() {
        if chain.reads.is_empty() {
            continue;
        }
        for &packed in &chain.reads {
            let bubble_id = (packed >> 33) as u64;
            let (beg, sink, ba, _path) = get_bubbles(&bub, bubble_id);
            let bn = ba.len();
            let mut t_len: u64 = 0;
            let mut vocc: u64 = 0;
            for &m in ba {
                bf[m as usize] = 1;
                bf[(m ^ 1) as usize] = 1;
                t_len += unitig_graph.unitigs[(m >> 1) as usize].length as u64;
                if if_hom(m >> 1, &bub) {
                    continue;
                }
                if unitig_graph.graph.sequences[(m >> 1) as usize].deleted {
                    continue;
                }
                vocc += unitig_graph.unitigs[(m >> 1) as usize].reads.len() as u64;
            }
            if beg != u32::MAX {
                t_len += unitig_graph.unitigs[(beg >> 1) as usize].length as u64;
            }
            if sink != u32::MAX {
                t_len += unitig_graph.unitigs[(sink >> 1) as usize].length as u64;
            }

            if vocc != 0 {
                // Collect every bubble vertex (both strands) whose branch is
                // fully inside this bubble chain and small enough.
                buf.clear();
                for &m in ba {
                    let mut v = m;
                    if unitig_graph.graph.sequences[(v >> 1) as usize].deleted {
                        continue;
                    }
                    for _direction in 0..2u32 {
                        if unitig_graph.graph.arc_count(v) < 2 {
                            v ^= 1;
                            continue;
                        }
                        if get_real_length(&unitig_graph.graph, v) < 2 {
                            v ^= 1;
                            continue;
                        }
                        if probe_bubble!(&mut unitig_graph.graph, v, t_len, &mut b, false) != 0 {
                            let b_sink = b.sources[0];
                            let mut socc = 0u64;
                            let mut z = 0usize;
                            while z < b.visited.len() {
                                let node = b.visited[z];
                                if node == v || node == b_sink {
                                    z += 1;
                                    continue;
                                }
                                socc += unitig_graph.unitigs[(node >> 1) as usize].reads.len()
                                    as u64;
                                z += 1;
                                if bf[node as usize] == 0 && bf[(node ^ 1) as usize] == 0 {
                                    break;
                                }
                            }
                            if z < b.visited.len() {
                                v ^= 1;
                                continue;
                            }
                            buf.push((socc << 32) | u64::from(v));
                        }
                        v ^= 1;
                    }
                }

                radix_sort_u64(&mut buf);
                let mut pocc = 0u64;
                for &entry in buf.iter() {
                    let v = entry as u32;
                    if unitig_graph.graph.sequences[(v >> 1) as usize].deleted {
                        continue;
                    }
                    if unitig_graph.graph.arc_count(v) < 2 {
                        continue;
                    }
                    if get_real_length(&unitig_graph.graph, v) < 2 {
                        continue;
                    }
                    if probe_bubble!(&mut unitig_graph.graph, v, t_len, &mut b, false) == 0 {
                        continue;
                    }
                    let b_sink = b.sources[0];
                    let mut socc = 0u64;
                    for &node in &b.visited {
                        if node == v || node == b_sink {
                            continue;
                        }
                        socc += unitig_graph.unitigs[(node >> 1) as usize].reads.len() as u64;
                    }
                    if (pocc + socc) as f64 >= vocc as f64 * bub_rate {
                        continue;
                    }
                    pocc += socc;
                    probe_bubble!(&mut unitig_graph.graph, v, t_len, &mut b, true);
                    for &node in &b.visited {
                        if node == v || node == b_sink {
                            continue;
                        }
                        if unitig_graph.graph.sequences[(node >> 1) as usize].deleted {
                            continue;
                        }
                        if unitig_graph.graph.sequences[(node >> 1) as usize].coverage_class
                            != ALTER_LABEL
                        {
                            continue;
                        }
                        let unitig = &unitig_graph.unitigs[(node >> 1) as usize];
                        if unitig.reads.is_empty() {
                            continue;
                        }
                        let reads: Vec<u64> = unitig.reads.clone();
                        for &read_packed in &reads {
                            sg.delete_sequence((read_packed >> 33) as u32);
                        }
                        unitig_graph.graph.delete_sequence(node >> 1);
                        unitig_graph.unitigs[(node >> 1) as usize].reads.clear();
                    }
                    n0 += 1;
                }
            }
            for &m in ba {
                bf[m as usize] = 0;
                bf[(m ^ 1) as usize] = 0;
            }
        }
    }

    // Every remaining alternate unitig loses its reads.
    for v in 0..unitig_graph.graph.sequences.len() {
        if unitig_graph.graph.sequences[v].deleted {
            continue;
        }
        if unitig_graph.graph.sequences[v].coverage_class != ALTER_LABEL {
            continue;
        }
        if unitig_graph.unitigs[v].reads.is_empty() {
            continue;
        }
        let reads: Vec<u64> = unitig_graph.unitigs[v].reads.clone();
        for &read_packed in &reads {
            sg.delete_sequence((read_packed >> 33) as u32);
        }
        unitig_graph.graph.delete_sequence(v as u32);
        unitig_graph.unitigs[v].reads.clear();
    }
    unitig_graph.graph.cleanup();

    // Type-1 pops: small bubbles whose both branches stay within short tips.
    let t_len =
        crate::assembly::ul_renew_gfa::get_bub_pop_max_dist_advance(&mut unitig_graph.graph);
    buf.clear();
    for v in 0..n_vtx as u32 {
        if unitig_graph.graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if unitig_graph.graph.arc_count(v) < 2 {
            continue;
        }
        if get_real_length(&unitig_graph.graph, v) < 2 {
            continue;
        }
        if probe_bubble!(&mut unitig_graph.graph, v, t_len, &mut b, false) != 0 {
            let b_sink = b.sources[0];
            let mut socc = 0u64;
            for &node in &b.visited {
                if node == v || node == b_sink {
                    continue;
                }
                socc += unitig_graph.unitigs[(node >> 1) as usize].reads.len() as u64;
            }
            if socc <= 16 {
                buf.push((socc << 32) | u64::from(v));
            }
        }
    }

    radix_sort_u64(&mut buf);
    for &entry in buf.iter() {
        let v = entry as u32;
        if unitig_graph.graph.sequences[(v >> 1) as usize].deleted {
            continue;
        }
        if unitig_graph.graph.arc_count(v) < 2 {
            continue;
        }
        if get_real_length(&unitig_graph.graph, v) < 2 {
            continue;
        }
        if probe_bubble!(&mut unitig_graph.graph, v, t_len, &mut b, false) == 0 {
            continue;
        }
        let b_sink = b.sources[0];
        let w = b_sink ^ 1;
        let mut socc = 0u64;
        for &node in &b.visited {
            if node == v || node == b_sink {
                continue;
            }
            socc += unitig_graph.unitigs[(node >> 1) as usize].reads.len() as u64;
        }
        if socc > 16 {
            continue;
        }

        b.visited.clear();
        let mut convex = 0u32;
        let mut node_len = 0i64;
        let mut base_len = 0i64;
        let mut max_stop_node_len = 0i64;
        let mut max_stop_base_len = 0i64;
        crate::assembly::bubble_identify::get_unitig(
            &unitig_graph.graph,
            None,
            v ^ 1,
            &mut convex,
            &mut node_len,
            &mut base_len,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            1,
            Some(&mut b.visited),
        );
        let mut vocc = 0u64;
        let mut stop = false;
        for &node in &b.visited {
            if if_hom(node >> 1, &bub) {
                stop = true;
                break;
            }
            vocc += unitig_graph.unitigs[(node >> 1) as usize].reads.len() as u64;
        }
        if stop {
            continue;
        }
        if socc as f64 >= vocc as f64 * bub_rate {
            continue;
        }

        b.visited.clear();
        crate::assembly::bubble_identify::get_unitig(
            &unitig_graph.graph,
            None,
            w ^ 1,
            &mut convex,
            &mut node_len,
            &mut base_len,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            1,
            Some(&mut b.visited),
        );
        let mut vocc = 0u64;
        let mut stop = false;
        for &node in &b.visited {
            if if_hom(node >> 1, &bub) {
                stop = true;
                break;
            }
            vocc += unitig_graph.unitigs[(node >> 1) as usize].reads.len() as u64;
        }
        if stop {
            continue;
        }
        if socc as f64 >= vocc as f64 * bub_rate {
            continue;
        }

        probe_bubble!(&mut unitig_graph.graph, v, t_len, &mut b, true);
        for &node in &b.visited {
            if node == v || node == b_sink {
                continue;
            }
            if unitig_graph.graph.sequences[(node >> 1) as usize].deleted {
                continue;
            }
            if unitig_graph.graph.sequences[(node >> 1) as usize].coverage_class != ALTER_LABEL {
                continue;
            }
            let unitig = &unitig_graph.unitigs[(node >> 1) as usize];
            if unitig.reads.is_empty() {
                continue;
            }
            let reads: Vec<u64> = unitig.reads.clone();
            for &read_packed in &reads {
                sg.delete_sequence((read_packed >> 33) as u32);
            }
            unitig_graph.graph.delete_sequence(node >> 1);
            unitig_graph.unitigs[(node >> 1) as usize].reads.clear();
        }
        n1 += 1;
    }

    // Every remaining alternate unitig loses its reads.
    for v in 0..unitig_graph.graph.sequences.len() {
        if unitig_graph.graph.sequences[v].deleted {
            continue;
        }
        if unitig_graph.graph.sequences[v].coverage_class != ALTER_LABEL {
            continue;
        }
        if unitig_graph.unitigs[v].reads.is_empty() {
            continue;
        }
        let reads: Vec<u64> = unitig_graph.unitigs[v].reads.clone();
        for &read_packed in &reads {
            sg.delete_sequence((read_packed >> 33) as u32);
        }
        unitig_graph.graph.delete_sequence(v as u32);
        unitig_graph.unitigs[v].reads.clear();
    }
    unitig_graph.graph.cleanup();

    sg.cleanup();
    let _ = (n0, n1);
}

/// `deep_graph_clean` (`gfa_ut.cpp` 18184).
#[allow(clippy::too_many_arguments)]
pub fn deep_graph_clean(
    options: &mut AssemblyOptions,
    sg: &mut Asg,
    trio_flags: &[u8],
    sources: &[PafList],
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    coverage_cut: &mut [CoverageCut],
    mut bubble_labels: Option<&mut crate::assembly::bubble::BubbleLabel>,
    is_ou: bool,
    is_trio: bool,
    max_ext: i64,
    min_ovlp_drop_ratio: f64,
    max_ovlp_drop_ratio: f64,
    ou_rat: f64,
    min_ou: i64,
    clean_round: i64,
    long_tip: u32,
) {
    hic_clean_adv(
        options,
        sg,
        sources,
        reverse_sources,
        read_to_unitig,
        coverage_cut,
        bubble_labels.as_deref_mut(),
    );
    let step = if clean_round == 1 {
        max_ovlp_drop_ratio
    } else {
        (max_ovlp_drop_ratio - min_ovlp_drop_ratio) / (clean_round - 1) as f64
    };
    let mut sl = UgClean::new(
        sg,
        trio_flags.to_vec(),
        is_ou,
        is_trio,
        max_ext,
        min_ovlp_drop_ratio,
        ou_rat,
        min_ou,
    );
    let mut b: Vec<u64> = Vec::new();
    for round in 0..clean_round {
        sl.len_rat = min_ovlp_drop_ratio + step * round as f64;
        if sl.len_rat > max_ovlp_drop_ratio {
            sl.len_rat = max_ovlp_drop_ratio;
        }
        update_ug_clean(&mut sl);
        for i in 0..sl.unitig_graph.graph.sequences.len() {
            cal_bub_best(&mut sl, i);
        }
        let (is_ou, len_rat, ou_rat) = (sl.is_ou, sl.len_rat, sl.ou_rat);
        cal_bub_best_by_len(
            &mut sl,
            Some(&mut b),
            max_ext as u32,
            is_trio,
            is_ou,
            len_rat,
            ou_rat,
            min_ou as u32,
            long_tip,
        );
    }
    hic_clean_adv(
        options,
        sg,
        sources,
        reverse_sources,
        read_to_unitig,
        coverage_cut,
        bubble_labels.as_deref_mut(),
    );
}

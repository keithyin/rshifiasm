//! Coverage inference and duplicate purging, the counterpart of
//! `Purge_Dups.cpp` and of `set_hom_global_coverage` in `Overlaps.cpp`.
//!
//! The cleaning stage calls `set_hom_global_coverage` once before the bubble
//! rescue round; it builds a unitig graph and runs `purge_dups` with
//! `just_coverage = 1`, which only infers the homozygous coverage from the
//! read overlap histogram and the per-unitig coverage bins, then exits through
//! `end_coverage`. The output stage (`output_bp_graph`) later runs the full
//! `purge_dups` on the copied primary graph.
//!
//! The driver, the haplotype-alignment scan, the hap-overlap normalization,
//! the maxcut partitioning, the purge-graph cleanup and the merge/link pass
//! are ported. The only remaining fallback is `if_ploid_sample`, reached only
//! when neither a k-mer nor a coverage peak is available.

use crate::assembly::asg::{Asg, AsgArc, MaUtg};
use crate::assembly::bubble_identify::get_unitig;
use crate::assembly::gfa_clean::{get_real_length, LOOP};
use crate::assembly::hap_alignment::{
    determine_hap_overlap_type_advance, filter_secondary_chain, get_xy_pos, HapAlignmentContext,
    HapCandidate,
};
use crate::assembly::hap_cov::{
    chain_origin_trans_uid_by_distance, ArcOffset, HapCov, N_HET, P_HET, RC_2,
};
use crate::assembly::hit::{
    get_specific_overlap, ma_hit_to_arc, CoverageCut, ReadToUnitig, ALTER_LABEL, PRIMARY_LABEL,
};
use crate::assembly::trans_cov::quick_lis;
use crate::assembly::trio_clean::asg_seq_drop;
use crate::assembly::unitig::{ma_ug_gen, UnitigGraph};
use crate::command_line::options::{AssemblyOptions, FLAG_PARTITION};
use crate::correction::ec_driver::PafList;
use crate::counting::peak_analysis::adjust_homozygous_peak;
use crate::container::ksort::radix_sort_by_key;
use crate::overlap::paf::PafHit;

/// `COV_COUNT`: the coverage histogram is one bin per supported depth.
pub const COVERAGE_COUNT: i64 = 1024;
/// `HOM_PEAK_RATE`: a homozygous peak implies this much read coverage.
pub const HOM_PEAK_RATE: f64 = 1.25;
/// `HET_PEAK_RATE`.
pub const HET_PEAK_RATE: f64 = HOM_PEAK_RATE * 2.0;
/// `ALTER_COV_THRES`: minimum alter fraction of a read's covered bases.
pub const ALTER_COVERAGE_THRESHOLD: f64 = 0.9;
/// `REAL_ALTER_THRES`: minimum alternate-to-primary base ratio.
pub const REAL_ALTER_THRESHOLD: f64 = 0.25;

/// The peaks reported by `get_read_peak`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadPeaks {
    /// `hom_peak`; `-1` when none was found.
    pub hom_peak: i64,
    /// `het_peak`; `-1` when none was found.
    pub het_peak: i64,
    /// `k_mer_only`: the homozygous k-mer peak.
    pub k_mer_only: i64,
    /// `coverage_only`: the homozygous coverage peak.
    pub coverage_only: i64,
}

impl Default for ReadPeaks {
    fn default() -> Self {
        Self {
            hom_peak: -1,
            het_peak: -1,
            k_mer_only: -1,
            coverage_only: -1,
        }
    }
}

/// `get_read_peak`: locates the homozygous and heterozygous peaks of the
/// read coverage histogram.
///
/// `topo_peak_cov` is the alternate peak reported by [`get_alter_peak`], or
/// `None` when the reference would pass a null pointer.
pub fn get_read_peak(
    read_graph: &Asg,
    coverage_histogram: &[i64],
    topo_peak_cov: Option<i64>,
    options: &AssemblyOptions,
) -> ReadPeaks {
    let cov_buf = coverage_histogram;
    let cov_buf_length = cov_buf.len() as i64;

    // `cov_buf[0]` is usually very large: find the first rising step, which
    // ends the error peak.
    let mut index = 1usize;
    while index < cov_buf.len() && cov_buf[index] <= cov_buf[index - 1] {
        index += 1;
    }
    // The reference reads `cov_buf[err_i + 1]` unconditionally. A histogram
    // without any rising step leaves that index one past the end, which is
    // undefined in C and never happens for a real coverage histogram.
    assert!(index < cov_buf.len(), "the coverage histogram has no peak");
    let err_i = index as i64 - 1;

    // The global highest peak.
    let mut max_i = err_i + 1;
    let mut max = cov_buf[max_i as usize];
    let mut i = max_i;
    while i < cov_buf_length {
        let value = cov_buf[i as usize];
        if value > max {
            max = value;
            max_i = i;
        }
        i += 1;
    }

    // Look for a smaller peak on the low end.
    let mut max2: i64 = -1;
    let mut max2_i: i64 = -1;
    let mut i = max_i - 1;
    while i > err_i {
        let value = cov_buf[i as usize];
        if value >= cov_buf[(i - 1) as usize] && value >= cov_buf[(i + 1) as usize] && value > max2
        {
            max2 = value;
            max2_i = i;
        }
        i -= 1;
    }
    if max2_i != -1 && max2_i > err_i && max2_i < max_i {
        let mut min = max;
        let mut j = max2_i + 1;
        while j < max_i {
            if cov_buf[j as usize] < min {
                min = cov_buf[j as usize];
            }
            j += 1;
        }
        // An insignificant second peak is discarded.
        if (max2 as f64) < (max as f64) * 0.05 || (min as f64) > (max2 as f64) * 0.95 {
            max2 = -1;
            max2_i = -1;
        }
    }
    if (max2 as f64) < (max as f64) * 0.0075 {
        // `max2` is not read again; only its index matters downstream.
        max2_i = -1;
    }

    // Look for a smaller peak on the high end. The last bin is ignored.
    let mut max3: i64 = -1;
    let mut max3_i: i64 = -1;
    let mut i = max_i + 1;
    while i < cov_buf_length - 1 {
        let value = cov_buf[i as usize];
        if value >= cov_buf[(i - 1) as usize] && value >= cov_buf[(i + 1) as usize] && value > max3
        {
            max3 = value;
            max3_i = i;
        }
        i += 1;
    }
    if max3 != -1 && max3_i > max_i {
        let mut min = max;
        let mut j = max_i + 1;
        while j < max3_i {
            if cov_buf[j as usize] < min {
                min = cov_buf[j as usize];
            }
            j += 1;
        }
        if (max3 as f64) < (max as f64) * 0.05
            || (min as f64) > (max3 as f64) * 0.95
            || max3_i > max_i * 3
        {
            max3 = -1;
            max3_i = -1;
        }
    }
    if (max3 as f64) < (max as f64) * 0.0075 {
        // `max3` is not read again; only its index matters downstream.
        max3_i = -1;
    }

    let mut peaks = ReadPeaks::default();

    let mut topo_peak_i: i64 = -1;
    if let Some(peak) = topo_peak_cov {
        // The reference dereferences the pointer only when the caller passed
        // one; [`get_alter_peak`] reports either `-1` or a valid bin.
        if (0..cov_buf_length).contains(&peak) {
            topo_peak_i = peak;
            let topo_peak = cov_buf[peak as usize];
            if (topo_peak as f64) <= (max as f64) * 0.05 {
                topo_peak_i = -1;
            }
        }
    }

    if options.purge_level_primary == 0 {
        peaks.hom_peak = max_i;
        return peaks;
    }

    let k_mer_het = i64::from(options.heterozygous_coverage);
    let k_mer_hom = i64::from(options.homozygous_coverage);
    let alter_peak = topo_peak_i;
    let mut coverage_het: i64;
    let mut coverage_hom: i64;
    if max3_i > 0 {
        coverage_het = max_i;
        coverage_hom = max3_i;
    } else {
        coverage_het = max2_i;
        coverage_hom = max_i;
    }

    if options.haploid_genome_size > 0 {
        let mut n_bs: i64 = 0;
        for sequence in &read_graph.sequences {
            n_bs += i64::from(sequence.length);
        }
        let m_peak_hom = n_bs / options.haploid_genome_size;
        if m_peak_hom > 0 {
            let mut p_ht: i32 = -1;
            coverage_hom = i64::from(adjust_homozygous_peak(
                m_peak_hom as i32,
                max_i as i32,
                max2_i as i32,
                max3_i as i32,
                &mut p_ht,
            ));
            coverage_het = i64::from(p_ht);
        }
    }

    if k_mer_het != -1 {
        peaks.het_peak = k_mer_het;
        peaks.hom_peak = k_mer_hom;
        return peaks;
    } else if coverage_het != -1 {
        peaks.het_peak = coverage_het;
        peaks.hom_peak = coverage_hom;
        return peaks;
    } else if (k_mer_hom as f64) > (coverage_hom as f64) * 1.5 {
        peaks.het_peak = coverage_hom;
        peaks.hom_peak = k_mer_hom;
        return peaks;
    } else if alter_peak != -1 {
        // On a heterozygous peak the coverage peak is more reliable; on a
        // homozygous one the k-mer peak is.
        if (coverage_hom as f64) >= (alter_peak as f64) * 0.8
            && (coverage_hom as f64) <= (alter_peak as f64) * 1.2
        {
            peaks.het_peak = coverage_hom;
            return peaks;
        } else if (k_mer_hom as f64) >= (alter_peak as f64) * 0.8 * 2.0
            && (k_mer_hom as f64) <= (alter_peak as f64) * 1.2 * 2.0
        {
            peaks.hom_peak = k_mer_hom;
            return peaks;
        }
    }

    peaks.k_mer_only = k_mer_hom;
    peaks.coverage_only = coverage_hom;
    peaks
}

/// `get_alter_peak`: the peak of the per-read coverage attributed to
/// alternate-unitig overlaps, or `-1`.
///
/// `position_index` mirrors `cov->pos_idx`: primary reads are marked with a
/// zero entry, every other read keeps `u64::MAX`. The reference resets the
/// whole array on entry and on exit.
pub fn get_alter_peak(
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
    read_to_unitig: &ReadToUnitig,
    position_index: &mut [u64],
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
) -> i64 {
    let unitig_read_graph = &unitig_graph.graph;
    let cov_buf_length = COVERAGE_COUNT;
    let mut cov_buf = vec![0i64; cov_buf_length as usize];
    let read_count = read_graph.sequences.len();
    let mut primary_bases: u64 = 0;
    let mut alter_bases: u64 = 0;

    position_index[..read_count].fill(u64::MAX);

    for v in 0..unitig_read_graph.sequences.len() {
        if unitig_read_graph.sequences[v].deleted {
            continue;
        }
        if unitig_read_graph.sequences[v].coverage_class == ALTER_LABEL {
            continue;
        }
        let unitig = &unitig_graph.unitigs[v];
        if unitig.reads.is_empty() {
            continue;
        }
        for entry in &unitig.reads {
            let query = (entry >> 33) as usize;
            position_index[query] = 0;
            let cut = &coverage_cut[query];
            let r_bases = cut.end.wrapping_sub(cut.start) as i64;
            primary_bases = primary_bases.wrapping_add(r_bases as u64);
        }
    }

    for query in 0..read_count {
        if position_index[query] == 0 {
            continue;
        }
        if read_graph.sequences[query].deleted {
            continue;
        }

        let mut c_bases_primary: i64 = 0;
        let mut c_bases_alter: i64 = 0;
        let cut = &coverage_cut[query];
        let r_bases = cut.end.wrapping_sub(cut.start) as i64;
        alter_bases = alter_bases.wrapping_add(r_bases as u64);

        for hit in &sources[query].hits {
            if !hit.exact_region {
                continue;
            }
            let mut target = hit.target_identifier;
            if read_graph.sequences[target as usize].deleted {
                // Resolve the deleted target to the read that contains it.
                let (resolved, is_unitig) = read_to_unitig.get(target);
                target = resolved;
                if target == u32::MAX || is_unitig == 1 {
                    continue;
                }
                if read_graph.sequences[target as usize].deleted {
                    continue;
                }
            }

            let covered = hit.query_end.wrapping_sub(hit.query_start) as i64;
            if position_index[target as usize] == 0 {
                c_bases_primary += covered;
            } else {
                c_bases_alter += covered;
            }
        }

        let c_bases = c_bases_primary + c_bases_alter;
        if (c_bases_alter as f64) < (c_bases as f64) * ALTER_COVERAGE_THRESHOLD {
            continue;
        }
        let c_bases = c_bases / r_bases;
        if c_bases < 0 || c_bases >= cov_buf_length {
            continue;
        }
        cov_buf[c_bases as usize] += 1;
    }

    let mut max_i: i64 = -1;
    let mut max: i64 = -1;
    for (bin, &count) in cov_buf.iter().enumerate() {
        if count > max {
            max = count;
            max_i = bin as i64;
        }
    }

    if (alter_bases as f64) < (primary_bases as f64) * REAL_ALTER_THRESHOLD {
        max_i = -1;
    }

    position_index[..read_count].fill(u64::MAX);
    max_i
}

/// `get_read_coverage_thres`: infers the homozygous read coverage.
///
/// Returns the threshold and the two peaks the reference reports separately
/// (`k_mer_only` and `coverage_only`) for its ploidy fallback.
pub fn get_read_coverage_thres(
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
    read_to_unitig: &ReadToUnitig,
    position_index: &mut [u64],
    sources: &[PafList],
    coverage_cut: &[CoverageCut],
    options: &AssemblyOptions,
) -> (i64, i64, i64) {
    let cov_buf_length = COVERAGE_COUNT;
    let mut cov_buf = vec![0i64; cov_buf_length as usize];
    let read_count = read_graph.sequences.len();

    for i in 0..read_count {
        let mut c_bases: i64 = 0;
        let r_bases = coverage_cut[i].end.wrapping_sub(coverage_cut[i].start) as i64;
        for hit in &sources[i].hits {
            if !hit.exact_region {
                continue;
            }
            c_bases += hit.query_end.wrapping_sub(hit.query_start) as i64;
        }
        let c_bases = c_bases / r_bases;
        if c_bases < 0 || c_bases >= cov_buf_length {
            continue;
        }
        cov_buf[c_bases as usize] += 1;
    }

    let alter_peak = get_alter_peak(
        unitig_graph,
        read_graph,
        read_to_unitig,
        position_index,
        sources,
        coverage_cut,
    );
    let topo_peak_cov = if alter_peak == -1 { None } else { Some(alter_peak) };
    let peaks = get_read_peak(read_graph, &cov_buf, topo_peak_cov, options);

    let mut threshold: i64 = -1;
    if peaks.hom_peak != -1 {
        threshold = (peaks.hom_peak as f64 * HOM_PEAK_RATE) as i64;
    } else if peaks.het_peak != -1 {
        threshold = (peaks.het_peak as f64 * HET_PEAK_RATE) as i64;
    }
    (threshold, peaks.k_mer_only, peaks.coverage_only)
}

/// `set_hom_global_coverage`: fills `asm_opt.hom_global_coverage` from the
/// read coverage before the bubble rescue round.
///
/// `maximum_hang` and `minimum_overlap` are the `purge_dups` arguments; the
/// coverage-only path never uses them, but they stay in the signature to keep
/// the call site self-describing.
///
/// # Panics
///
/// When neither a k-mer peak nor a coverage peak is available the reference
/// falls back to `if_ploid_sample`, which walks the hap-alignment overlap
/// pipeline (`Purge_Dups.cpp:4800`). That pipeline is not ported yet and the
/// supported pipelines always report a k-mer peak, so hitting the fallback
/// panics instead of silently using a wrong coverage.
pub fn set_hom_global_coverage(
    options: &mut AssemblyOptions,
    read_graph: &mut Asg,
    coverage_cut: &[CoverageCut],
    sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    _maximum_hang: i32,
    _minimum_overlap: i32,
) {
    if options.global_homozygous_coverage_is_set != 0 {
        return;
    }

    // `ma_ug_gen(sg)`: the unitig graph used by the purge loop and by the
    // alternate-peak scan.
    let mut unitig_graph = ma_ug_gen(read_graph);

    // `init_hap_cov_t`: of the state it allocates, the coverage-only exit only
    // ever uses `pos_idx`; `cov` and `is_r_het` feed the hap-alignment
    // machinery that `just_coverage` skips.
    let read_count = read_graph.sequences.len();
    let mut position_index = vec![u64::MAX; read_count];

    // The reference computes `hap_buf.cov_threshold` before filling the read
    // map, so `get_alter_peak` sees the read map as the rescues left it.
    let mut cov_threshold: i64 = if options.global_homozygous_coverage != -1 {
        if options.global_homozygous_coverage_is_set != 0 {
            (f64::from(options.global_homozygous_coverage) * HOM_PEAK_RATE) as i64
        } else {
            i64::from(options.global_homozygous_coverage)
        }
    } else {
        get_read_coverage_thres(
            &unitig_graph,
            read_graph,
            read_to_unitig,
            &mut position_index,
            sources,
            coverage_cut,
            options,
        )
        .0
    };

    // Map every read of every unitig to it; `purge_dups` also stores the
    // read's position inside the unitig, which `get_alter_peak` overwrites.
    for v in 0..unitig_graph.graph.sequences.len() {
        let unitig_id = v as u32;
        let mut offset: u32 = 0;
        let unitig = &unitig_graph.unitigs[v];
        for (i, entry) in unitig.reads.iter().enumerate() {
            let read = (entry >> 33) as u32;
            let class = read_graph.sequences[read as usize].coverage_class;
            read_to_unitig.set(read, unitig_id, true, Some(class));
            position_index[read as usize] = (u64::from(offset) << 32) | (i as u64);
            offset = offset.wrapping_add(*entry as u32);
        }
    }

    if cov_threshold < 0 {
        // `if_ploid_sample` decides between the coverage-only and k-mer-only
        // peaks from the purge overlap graph; see the function docs.
        panic!(
            "purge_dups needs if_ploid_sample (Purge_Dups.cpp:4800), which is not ported yet"
        );
    }
    if options.global_homozygous_coverage == -1 {
        options.global_homozygous_coverage = cov_threshold as i32;
    }
    if options.global_purge_coverage != -1 {
        cov_threshold = i64::from(options.global_purge_coverage);
    }

    eprintln!(
        "[M::purge_dups] homozygous read coverage threshold: {}",
        if options.global_homozygous_coverage_is_set != 0 {
            options.global_homozygous_coverage
        } else {
            (f64::from(options.global_homozygous_coverage) / HOM_PEAK_RATE) as i32
        }
    );
    eprintln!(
        "[M::purge_dups] purge duplication coverage threshold: {}",
        cov_threshold
    );

    // `if(just_coverage) goto end_coverage;`
    read_to_unitig.clear_unitig_entries();
    unitig_graph.graph.cleanup();
}

// ---------------------------------------------------------------------------
// The full `purge_dups`.
// ---------------------------------------------------------------------------

/// `SELF_EXIST`: the overlap was found by the scan of unitig `xUid`.
pub const SELF_EXIST: u8 = 0;
/// `REVE_EXIST`: the overlap is the reverse-side copy of another one.
pub const REVE_EXIST: u8 = 1;
/// `DELETE`: the overlap was rejected by a later filter.
pub const DELETE: u8 = 2;
/// `MIXED`.
pub const MIXED: u8 = 3;
/// `FLIP`: the overlap conflicts with the min-cut partition.
pub const FLIP: u8 = 4;

/// `X2Y`.
pub const X2Y: u8 = 0;
/// `Y2X`.
pub const Y2X: u8 = 1;
/// `XCY`.
pub const XCY: u8 = 2;
/// `YCX`.
pub const YCX: u8 = 3;

/// `REV_W`: the weight `init_pip_hh` gives a read with no reverse coverage.
pub const REV_W: i64 = 8;
/// `SET_TRIO`: the read was seen while the primary graph was partitioned.
pub const SET_TRIO: u8 = 8;

/// `hap_overlaps`: one haplotype-alignment chain between two unitigs.
#[derive(Debug, Clone, Copy, Default)]
pub struct HapOverlap {
    /// `rev`.
    pub reverse: u8,
    /// `type`.
    pub kind: u8,
    /// `status`.
    pub status: u8,
    /// `x_beg_pos`.
    pub x_beg_pos: u32,
    /// `x_end_pos`.
    pub x_end_pos: u32,
    /// `y_beg_pos`.
    pub y_beg_pos: u32,
    /// `y_end_pos`.
    pub y_end_pos: u32,
    /// `x_beg_id`.
    pub x_beg_id: u32,
    /// `x_end_id`.
    pub x_end_id: u32,
    /// `y_beg_id`.
    pub y_beg_id: u32,
    /// `y_end_id`.
    pub y_end_id: u32,
    /// `xUid`.
    pub x_uid: u32,
    /// `yUid`.
    pub y_uid: u32,
    /// `weight`.
    pub weight: u32,
    /// `score`.
    pub score: i64,
    /// `s`, the pair similarity.
    pub similarity: f32,
}

/// `hap_overlaps_list`: one overlap list per unitig.
#[derive(Debug, Clone, Default)]
pub struct HapOverlapsList {
    /// `x`.
    pub x: Vec<Vec<HapOverlap>>,
}

impl HapOverlapsList {
    /// `init_hap_overlaps_list`.
    pub fn new(num: usize) -> Self {
        Self { x: vec![Vec::new(); num] }
    }

    /// Whether every per-unitig list is empty.
    pub fn is_empty(&self) -> bool {
        self.x.iter().all(|list| list.is_empty())
    }
}

/// `p_node_t`: a run of reads of one unitig with the same heterozygosity flag.
#[derive(Debug, Clone, Copy, Default)]
pub struct PNode {
    /// `baseBeg`.
    pub base_beg: u32,
    /// `baseEnd`.
    pub base_end: u32,
    /// `nodeBeg`.
    pub node_beg: u32,
    /// `nodeEnd`.
    pub node_end: u32,
    /// `h_lev_idx`, unused by the reference.
    pub h_lev_idx: u32,
    /// `h_status`.
    pub h_status: u32,
    /// `c_ug_id`.
    pub c_ug_id: u32,
}

/// `p_g_in_t`: the slice of `pg_het_node` that belongs to one unitig.
#[derive(Debug, Clone, Copy, Default)]
pub struct PGIn {
    /// `beg`.
    pub beg: u32,
    /// `occ`.
    pub occ: u32,
}

/// `p_g_t`: the purge-level graph and its het/hom intervals.
#[derive(Debug, Clone, Default)]
pub struct PG {
    /// `pg_het_node`.
    pub pg_het_node: Vec<PNode>,
    /// `pg_h_lev`.
    pub pg_h_lev: Asg,
    /// `pg_h_lev_idx`.
    pub pg_h_lev_idx: Vec<PGIn>,
}

/// `init_p_g_t`: builds the unitig-level purge graph and the per-unitig
/// heterozygous/homozygous read runs.
fn init_p_g_t(unitig_graph: &UnitigGraph, read_graph: &Asg, is_r_het: &[u8]) -> PG {
    let nsg = &unitig_graph.graph;
    let mut pg = PG::default();

    for v in 0..nsg.sequences.len() {
        if nsg.sequences[v].deleted || nsg.sequences[v].coverage_class == ALTER_LABEL {
            // `asg_seq_set(pg->pg_h_lev, uId, 0, 1)` followed by the class
            // store. `asg_seq_set` writes the deletion flag and length; the
            // class byte is assigned separately in the reference.
            pg.pg_h_lev.set_sequence(v, 0, true);
            pg.pg_h_lev.sequences[v].coverage_class = ALTER_LABEL;
            continue;
        }
        pg.pg_h_lev.set_sequence(v, unitig_graph.unitigs[v].length, false);
        pg.pg_h_lev.sequences[v].coverage_class = PRIMARY_LABEL;
    }

    for v in 0..nsg.sequences.len() {
        if nsg.sequences[v].deleted || nsg.sequences[v].coverage_class == ALTER_LABEL {
            continue;
        }
        let unitig = &unitig_graph.unitigs[v];
        let run_start = pg.pg_het_node.len() as u32;
        let n = unitig.reads.len();
        let mut k = 1usize;
        let mut l = 0usize;
        let mut offset: u32 = 0;
        let mut l_pos: u32 = 0;
        while k <= n {
            let end = k == n
                || is_r_het[(unitig.reads[k] >> 33) as usize]
                    != is_r_het[(unitig.reads[l] >> 33) as usize];
            if end {
                pg.pg_het_node.push(PNode {
                    c_ug_id: v as u32,
                    h_status: u32::from(is_r_het[(unitig.reads[l] >> 33) as usize]),
                    base_beg: l_pos,
                    base_end: offset
                        .wrapping_add(read_graph.sequences[(unitig.reads[k - 1] >> 33) as usize].length)
                        .wrapping_sub(1),
                    node_beg: l as u32,
                    node_end: (k - 1) as u32,
                    h_lev_idx: 0,
                });
                l = k;
                l_pos = offset.wrapping_add(unitig.reads[k - 1] as u32);
            }
            offset = offset.wrapping_add(unitig.reads[k - 1] as u32);
            k += 1;
        }
        let occ = pg.pg_het_node.len() as u32 - run_start;
        pg.pg_h_lev_idx.push(PGIn { beg: run_start, occ });
    }

    pg
}

/// `clean_visit_flag`: clears the vote marks of the unitigs a read list votes
/// for.
pub(crate) fn clean_visit_flag(
    visit: &mut [u8],
    read_graph: &Asg,
    read_to_unitig: &ReadToUnitig,
    contig_count: usize,
    hits: &[PafHit],
) {
    if hits.len() * 2 > contig_count {
        visit.fill(0);
        return;
    }

    for hit in hits {
        let mut r_id = hit.target_identifier;
        if read_graph.sequences[r_id as usize].deleted {
            // `get_R_to_U` rewrites a deleted target to the read that
            // contains it.
            let (mapped, is_unitig) = read_to_unitig.get(r_id);
            r_id = mapped;
            if r_id == u32::MAX
                || is_unitig == 1
                || read_graph.sequences[r_id as usize].deleted
            {
                continue;
            }
        }

        let (hap_c_id, is_unitig) = read_to_unitig.get(r_id);
        if is_unitig == 0 || hap_c_id == u32::MAX {
            continue;
        }
        visit[hap_c_id as usize] = 0;
    }
}

/// `init_hap_alignment_struct_pip` + `kt_for(hap_alignment_advance_worker)`,
/// single-threaded: builds the per-read `hh` weights, then runs the vote and
/// candidate scan for every unitig into `all_ovlp`.
#[allow(clippy::too_many_arguments)]
fn hap_alignment_advance(
    unitig_graph: &UnitigGraph,
    read_graph: &Asg,
    coverage_cut: &[CoverageCut],
    sources: &[PafList],
    reverse_sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    position_index: &[u64],
    density: f32,
    max_hang: i32,
    min_ovlp: i32,
    cov_threshold: i64,
    cov: &HapCov,
    options: &AssemblyOptions,
    all_ovlp: &mut HapOverlapsList,
) -> Result<(), String> {
    let _ = sources;
    // `init_hap_alignment_struct_pip`'s `sc` derivation.
    let mut sc: i64 = -1;
    if options.global_homozygous_coverage_is_set != 0 {
        sc = (f64::from(options.global_homozygous_coverage) * 1.75) as i64;
    } else if options.global_homozygous_coverage > 0 {
        sc = ((((f64::from(options.global_homozygous_coverage)) / HOM_PEAK_RATE) as i32) as f64
            * 1.75) as i64;
    }
    if sc <= 0 {
        sc = -1;
    }
    let hh = crate::assembly::hap_alignment::init_pip_hh(read_graph, reverse_sources, coverage_cut, sc);

    let ctx = crate::assembly::hap_alignment::HapAlignmentContext {
        read_graph,
        reverse_sources,
        read_to_unitig,
        coverage_cut,
        position_index,
        maximum_hang: max_hang,
        maximum_hang_rate: options.maximum_hang_rate,
        minimum_overlap: min_ovlp,
        purge_level: options.purge_level_primary as u32,
    };
    let unitig_count = unitig_graph.graph.sequences.len();
    let mut scratch = crate::assembly::hap_alignment::HapAlignmentScratch::new(unitig_count);
    let cov_threshold_u64 = if cov_threshold < 0 { u64::MAX } else { cov_threshold as u64 };
    for v in 0..unitig_count {
        crate::assembly::hap_alignment::hap_alignment_scan_worker(
            unitig_graph,
            &ctx,
            v as u32,
            &hh,
            density,
            0.1,
            cov_threshold_u64,
            cov,
            all_ovlp,
            &mut scratch,
            &|visit, contig_count, hits| {
                clean_visit_flag(
                    visit,
                    read_graph,
                    read_to_unitig,
                    contig_count as usize,
                    hits,
                );
            },
        );
    }
    Ok(())
}

/// `filter_hap_overlaps_by_length`: drops chains shorter than `min_len`.
fn filter_hap_overlaps_by_length(all_ovlp: &mut HapOverlapsList, min_len: u32) {
    if min_len == 0 {
        return;
    }
    for list in &mut all_ovlp.x {
        list.retain(|overlap| overlap.x_end_id.wrapping_sub(overlap.x_beg_id) >= min_len);
    }
}

/// `get_specific_hap_overlap`: the index of the `qn -> tn` overlap in `list`.
pub(crate) fn get_specific_hap_overlap(list: &[HapOverlap], qn: u32, tn: u32) -> Option<usize> {
    list.iter().position(|overlap| overlap.x_uid == qn && overlap.y_uid == tn)
}

/// `set_reverse_hap_overlap`: copies the reverse-side view of an overlap,
/// mapping the type through `X2Y <-> Y2X` and `XCY <-> YCX`.
pub(crate) fn set_reverse_hap_overlap(
    dest: &mut HapOverlap,
    source: &HapOverlap,
    types: &[u8; 4],
) {
    dest.status = REVE_EXIST;
    dest.reverse = source.reverse;
    dest.kind = types[source.kind as usize];
    dest.weight = source.weight;
    dest.x_uid = source.y_uid;
    dest.y_uid = source.x_uid;
    dest.x_beg_pos = source.y_beg_pos;
    dest.x_end_pos = source.y_end_pos;
    dest.y_beg_pos = source.x_beg_pos;
    dest.y_end_pos = source.x_end_pos;
    dest.x_beg_id = source.y_beg_id;
    dest.x_end_id = source.y_end_id;
    dest.y_beg_id = source.x_beg_id;
    dest.y_end_id = source.x_end_id;
    dest.score = source.score;
}

/// The het/hom overlap lengths of `[qs, qe]` against the read runs of `uid`.
fn p_g_interval_lengths(pg: &PG, polyploidy: i32, uid: u32, qs: u32, qe: u32) -> (u64, u64) {
    let entry = pg.pg_h_lev_idx[uid as usize];
    let beg = entry.beg as usize;
    let mut het_len = 0u64;
    let mut hom_len = 0u64;
    for node in &pg.pg_het_node[beg..beg + entry.occ as usize] {
        let a_s = node.base_beg;
        let a_e = node.base_end;
        let ovlp = if qe.min(a_e) >= qs.max(a_s) {
            u64::from(qe.min(a_e) - qs.max(a_s) + 1)
        } else {
            0
        };
        if hom_len + het_len > 0 && ovlp == 0 {
            break;
        }
        if ovlp == 0 {
            continue;
        }
        if node.h_status == u32::from(N_HET)
            || (polyploidy <= 2 && node.h_status & u32::from(P_HET) != 0)
        {
            hom_len += ovlp;
        } else {
            het_len += ovlp;
        }
    }
    (het_len, hom_len)
}

/// `normalize_hap_overlaps_advance_by_p_g_t`: drops the overlaps that
/// entirely cover homozygous runs and mirrors every surviving overlap onto its
/// partner unitig's list.
///
/// The partner lists grow while this runs; the reference walks each list with
/// its current length, so mirrored entries of later unitigs are visited too.
fn normalize_hap_overlaps_advance_by_p_g_t(
    all_ovlp: &mut [Vec<HapOverlap>],
    back_all_ovlp: &mut [Vec<HapOverlap>],
    pg: &PG,
    polyploidy: i32,
    filter_rate: f64,
) {
    let types = [Y2X, X2Y, YCX, XCY];

    for list in all_ovlp.iter_mut() {
        for overlap in list.iter_mut() {
            let qn = overlap.x_uid;
            let (het_len, hom_len) = p_g_interval_lengths(
                pg,
                polyploidy,
                qn,
                overlap.x_beg_pos,
                overlap.x_end_pos.wrapping_sub(1),
            );
            if (het_len as f64) <= ((het_len + hom_len) as f64 * filter_rate) {
                overlap.status = DELETE;
                continue;
            }

            let tn = overlap.y_uid;
            let (het_len, hom_len) = p_g_interval_lengths(
                pg,
                polyploidy,
                tn,
                overlap.y_beg_pos,
                overlap.y_end_pos.wrapping_sub(1),
            );
            if (het_len as f64) <= ((het_len + hom_len) as f64 * filter_rate) {
                overlap.status = DELETE;
            }
        }
    }

    for list in all_ovlp.iter_mut() {
        list.retain(|overlap| overlap.status != DELETE);
    }

    let num = all_ovlp.len();
    for v in 0..num {
        let mut i = 0usize;
        while i < all_ovlp[v].len() {
            let x = all_ovlp[v][i];
            let qn = x.x_uid;
            let tn = x.y_uid;
            match get_specific_hap_overlap(&all_ovlp[tn as usize], tn, qn) {
                Some(index) => {
                    let y = all_ovlp[tn as usize][index];
                    if x.reverse == y.reverse && types[x.kind as usize] == y.kind {
                        i += 1;
                        continue;
                    }
                    if x.score >= y.score {
                        back_all_ovlp[tn as usize].push(y);
                        let mut reversed = y;
                        set_reverse_hap_overlap(&mut reversed, &x, &types);
                        all_ovlp[tn as usize][index] = reversed;
                    } else {
                        back_all_ovlp[qn as usize].push(x);
                        let mut reversed = x;
                        set_reverse_hap_overlap(&mut reversed, &y, &types);
                        all_ovlp[v][i] = reversed;
                    }
                }
                None => {
                    let mut reversed = HapOverlap::default();
                    set_reverse_hap_overlap(&mut reversed, &x, &types);
                    all_ovlp[tn as usize].push(reversed);
                }
            }
            i += 1;
        }
    }
}

/// The read-graph context of [`chain_origin_trans_uid_by_purge`].
struct PurgeAlignContext<'a> {
    unitigs: &'a [MaUtg],
    read_graph: &'a Asg,
    reverse_sources: &'a [PafList],
    read_to_unitig: &'a ReadToUnitig,
    coverage_cut: &'a [CoverageCut],
    max_hang: i32,
    max_hang_rate: f32,
    min_ovlp: i32,
    purge_level: u32,
}

/// `chain_origin_trans_uid_by_purge`: re-derives both directions of one hap
/// overlap and records the stronger side as an `RC_2` transition chain.
fn chain_origin_trans_uid_by_purge(
    x: &HapOverlap,
    ctx: &PurgeAlignContext,
    cov: &mut HapCov,
) {
    let unitigs = ctx.unitigs;
    let mut hap_for = HapCandidate {
        rev: x.reverse,
        x_beg_pos: x.x_beg_id,
        x_end_pos: x.x_end_id.wrapping_sub(1),
        y_beg_pos: x.y_beg_id,
        y_end_pos: x.y_end_id.wrapping_sub(1),
        ..HapCandidate::default()
    };
    let mut hap_rev = HapCandidate {
        rev: x.reverse,
        x_beg_pos: x.y_beg_id,
        x_end_pos: x.y_end_id.wrapping_sub(1),
        y_beg_pos: x.x_beg_id,
        y_end_pos: x.x_end_id.wrapping_sub(1),
        ..HapCandidate::default()
    };

    let align_ctx = HapAlignmentContext {
        read_graph: ctx.read_graph,
        reverse_sources: ctx.reverse_sources,
        read_to_unitig: ctx.read_to_unitig,
        coverage_cut: ctx.coverage_cut,
        position_index: &cov.pos_idx,
        maximum_hang: ctx.max_hang,
        maximum_hang_rate: ctx.max_hang_rate,
        minimum_overlap: ctx.min_ovlp,
        purge_level: ctx.purge_level,
    };

    let (r_x, ..) = determine_hap_overlap_type_advance(
        &mut hap_for,
        &unitigs[x.x_uid as usize],
        &unitigs[x.y_uid as usize],
        x.x_uid,
        x.y_uid,
        &align_ctx,
        &mut cov.u_buffer,
        &mut cov.tail_index,
        &mut cov.prev_index,
    );
    let (r_y, ..) = determine_hap_overlap_type_advance(
        &mut hap_rev,
        &unitigs[x.y_uid as usize],
        &unitigs[x.x_uid as usize],
        x.y_uid,
        x.x_uid,
        &align_ctx,
        &mut cov.u_buffer,
        &mut cov.tail_index,
        &mut cov.prev_index,
    );

    if r_x == u32::MAX && r_y == u32::MAX {
        eprintln!("ERROR-purge");
        return;
    }

    hap_for.x_beg_pos = x.x_beg_id;
    hap_for.x_end_pos = x.x_end_id.wrapping_sub(1);
    hap_for.y_beg_pos = x.y_beg_id;
    hap_for.y_end_pos = x.y_end_id.wrapping_sub(1);
    hap_rev.x_beg_pos = x.y_beg_id;
    hap_rev.x_end_pos = x.y_end_id.wrapping_sub(1);
    hap_rev.y_beg_pos = x.x_beg_id;
    hap_rev.y_end_pos = x.x_end_id.wrapping_sub(1);

    let (pri_uid, aux_uid, forward) = if r_x != u32::MAX && r_y == u32::MAX {
        (x.x_uid, x.y_uid, true)
    } else if r_x == u32::MAX && r_y != u32::MAX {
        (x.y_uid, x.x_uid, false)
    } else if hap_for.score >= hap_rev.score {
        (x.x_uid, x.y_uid, true)
    } else {
        (x.y_uid, x.x_uid, false)
    };
    let hap = if forward { &mut hap_for } else { &mut hap_rev };

    let (_, x_pos_beg, _, y_pos_beg, _) = determine_hap_overlap_type_advance(
        hap,
        &unitigs[pri_uid as usize],
        &unitigs[aux_uid as usize],
        pri_uid,
        aux_uid,
        &align_ctx,
        &mut cov.u_buffer,
        &mut cov.tail_index,
        &mut cov.prev_index,
    );

    let pri_len = unitigs[pri_uid as usize].length;
    let aux_len = unitigs[aux_uid as usize].length;
    let pri_vertex = pri_uid << 1;
    let aux_vertex = (aux_uid << 1) | u32::from(hap.rev);
    let score = hap.score as f64;
    let pri_a = [pri_vertex];
    let aux_a = [aux_vertex];

    if let Some(t_ch) = cov.t_ch.as_mut() {
        chain_origin_trans_uid_by_distance(
            t_ch,
            unitigs,
            &mut cov.pos_idx,
            &mut cov.u_buffer,
            &mut cov.tail_index,
            &pri_a,
            x_pos_beg as u32,
            Some(pri_len),
            &aux_a,
            y_pos_beg as u32,
            Some(aux_len),
            RC_2,
            score,
        );
    }
}

/// `collect_purge_trans_cov`: records the transition chain of every hap
/// overlap whose `xUid <= yUid`.
fn collect_purge_trans_cov(
    all_ovlp: &[Vec<HapOverlap>],
    ctx: &PurgeAlignContext,
    cov: &mut HapCov,
) {
    for list in all_ovlp {
        for x in list {
            if x.y_uid < x.x_uid {
                continue;
            }
            chain_origin_trans_uid_by_purge(x, ctx, cov);
        }
    }
}

/// `collect_trans_purge_cov`: adds the auxiliary chain's depth to the primary
/// chain's reads.
fn collect_trans_purge_cov(
    unitigs: &[MaUtg],
    cov: &mut HapCov,
    x: &HapOverlap,
    is_keep_x: bool,
) {
    if unitigs[x.x_uid as usize].reads.is_empty() || unitigs[x.y_uid as usize].reads.is_empty() {
        return;
    }
    let (pri_uid, pri_beg, pri_end, aux_uid, aux_beg, aux_end) = if is_keep_x {
        (x.x_uid, x.x_beg_id, x.x_end_id, x.y_uid, x.y_beg_id, x.y_end_id)
    } else {
        (x.y_uid, x.y_beg_id, x.y_end_id, x.x_uid, x.x_beg_id, x.x_end_id)
    };
    let pri = &unitigs[pri_uid as usize].reads[pri_beg as usize..pri_end as usize];
    let aux = &unitigs[aux_uid as usize].reads[aux_beg as usize..aux_end as usize];

    let mut u_cov = 0u64;
    let mut u_len = 0u64;
    for &packed in aux {
        u_cov = u_cov.wrapping_add(u64::from(cov.cov[(packed >> 33) as usize]));
    }
    for &packed in pri {
        u_len = u_len.wrapping_add(u64::from(cov.read_lengths[(packed >> 33) as usize]));
    }
    let u_cov = if u_len == 0 { 0 } else { u_cov / u_len };
    for &packed in pri {
        let read = (packed >> 33) as usize;
        let add = u_cov.wrapping_mul(u64::from(cov.read_lengths[read])) as u32;
        cov.cov[read] = cov.cov[read].wrapping_add(add);
    }
}

/// `collect_trans_purge_joint_cov`: swaps the depth of the two chains.
fn collect_trans_purge_joint_cov(unitigs: &[MaUtg], cov: &mut HapCov, x: &HapOverlap) {
    if unitigs[x.x_uid as usize].reads.is_empty() || unitigs[x.y_uid as usize].reads.is_empty() {
        return;
    }
    let a0 = &unitigs[x.x_uid as usize].reads[x.x_beg_id as usize..x.x_end_id as usize];
    let a1 = &unitigs[x.y_uid as usize].reads[x.y_beg_id as usize..x.y_end_id as usize];

    let mut u_cov = [0u64; 2];
    let mut u_len = [0u64; 2];
    for &packed in a0 {
        let read = (packed >> 33) as usize;
        u_cov[0] = u_cov[0].wrapping_add(u64::from(cov.cov[read]));
        u_len[0] = u_len[0].wrapping_add(u64::from(cov.read_lengths[read]));
    }
    for &packed in a1 {
        let read = (packed >> 33) as usize;
        u_cov[1] = u_cov[1].wrapping_add(u64::from(cov.cov[read]));
        u_len[1] = u_len[1].wrapping_add(u64::from(cov.read_lengths[read]));
    }
    let u_depth = [
        if u_len[0] == 0 { 0 } else { u_cov[1] / u_len[0] },
        if u_len[1] == 0 { 0 } else { u_cov[0] / u_len[1] },
    ];

    for (chain, depth) in [(a0, u_depth[0]), (a1, u_depth[1])] {
        for &packed in chain {
            let read = (packed >> 33) as usize;
            let add = depth.wrapping_mul(u64::from(cov.read_lengths[read])) as u32;
            cov.cov[read] = cov.cov[read].wrapping_add(add);
        }
    }
}

/// `remove_contained_haplotig`: every contained chain weaker than 95% of the
/// unitig's best chain marks its query unitig alternate and drops it from the
/// purge graph.
fn remove_contained_haplotig(
    all_ovlp: &mut [Vec<HapOverlap>],
    unitig_graph: &mut UnitigGraph,
    purge_g: &mut Asg,
    cov: &mut HapCov,
) {
    for list in all_ovlp.iter_mut() {
        if list.is_empty() {
            continue;
        }
        let mut max_score = i64::MIN;
        for overlap in list.iter() {
            if overlap.score > max_score {
                max_score = overlap.score;
            }
        }
        for overlap in list.iter_mut() {
            if overlap.kind != YCX {
                continue;
            }
            if !filter_secondary_chain(max_score, overlap.score, 0.95) {
                continue;
            }
            let x_uid = overlap.x_uid;
            unitig_graph.graph.sequences[x_uid as usize].coverage_class = ALTER_LABEL;
            purge_g.sequences[x_uid as usize].coverage_class = ALTER_LABEL;
            purge_g.sequences[x_uid as usize].deleted = true;
            overlap.status = DELETE;
            collect_trans_purge_cov(&unitig_graph.unitigs, cov, overlap, false);
        }
    }
}

/// `get_hap_arch`: the purge-graph arc of one hap overlap. The `score` replaces
/// the node length in the arc so the purge passes can rank the edges.
fn get_hap_arch(
    x: &HapOverlap,
    x_len: u32,
    y_len: u32,
    max_hang: i32,
    max_hang_rate: f32,
    min_ovlp: i32,
) -> Option<AsgArc> {
    let hit = PafHit {
        query_identifier: x.x_uid,
        query_start: x.x_beg_pos,
        query_end: x.x_end_pos,
        target_identifier: x.y_uid,
        target_start: x.y_beg_pos,
        target_end: x.y_end_pos,
        reverse: x.reverse != 0,
        ..PafHit::default()
    };
    let (r, mut arc) = ma_hit_to_arc(
        &hit,
        x_len as i32,
        y_len as i32,
        max_hang,
        max_hang_rate,
        min_ovlp,
    );
    if r < 0 {
        return None;
    }

    let mut score = x.score.unsigned_abs();
    if x.score >= 0 {
        score += 1u64 << 63;
    }
    arc.overlap = score as u32;
    arc.ul = ((arc.ul >> 32) << 32) | (score >> 32);
    Some(arc)
}

/// `decode_score`.
fn decode_score(h_bits: u32, l_bits: u32) -> i64 {
    let x = (u64::from(h_bits) << 32) | u64::from(l_bits);
    let magnitude = ((x << 1) >> 1) as i64;
    if x >> 63 == 0 {
        -magnitude
    } else {
        magnitude
    }
}

/// `encode_score`.
fn encode_score(value: i64) -> (u32, u32) {
    let mut score = value.unsigned_abs();
    if value >= 0 {
        score += 1u64 << 63;
    }
    ((score >> 32) as u32, score as u32)
}

/// `binfo_t`, the per-vertex state of a purge-graph bubble pop.
#[derive(Debug, Clone, Copy, Default)]
struct PurgeBinfo {
    parent: u32,
    distance: u32,
    remaining: u32,
    state: bool,
    score_h: u32,
    score_l: u32,
}

/// `buf_t`, the purge-graph pop workspace.
#[derive(Debug, Clone, Default)]
struct PurgePopBuffer {
    info: Vec<PurgeBinfo>,
    sources: Vec<u32>,
    visited: Vec<u32>,
    visited_edges: Vec<u32>,
}

/// `asg_bub_backtrack_primary` on a purge graph.
fn purge_bub_backtrack_primary(graph: &mut Asg, v0: u32, buffer: &PurgePopBuffer) {
    let sink = buffer.sources[0];
    let saved_class = graph.sequences[(sink >> 1) as usize].coverage_class;

    for &node in &buffer.visited {
        graph.sequences[(node >> 1) as usize].coverage_class = ALTER_LABEL;
    }

    let mut v = sink;
    loop {
        let u = buffer.info[v as usize].parent;
        graph.sequences[(v >> 1) as usize].coverage_class = crate::assembly::hit::HAP_LABEL;
        v = u;
        if v == v0 {
            break;
        }
    }
    graph.sequences[(sink >> 1) as usize].coverage_class = saved_class;

    for &edge in &buffer.visited_edges {
        let arc = graph.arcs[edge as usize];
        let query = (arc.ul >> 33) as u32;
        let target = arc.v >> 1;
        if graph.sequences[query as usize].coverage_class == ALTER_LABEL
            && graph.sequences[target as usize].coverage_class == ALTER_LABEL
        {
            continue;
        }
        graph.arcs[edge as usize].deleted = true;
        graph.delete_arc(arc.v ^ 1, (arc.ul >> 32) as u32 ^ 1, true);
    }

    let mut v = sink;
    loop {
        let u = buffer.info[v as usize].parent;
        graph.sequences[(v >> 1) as usize].deleted = false;
        graph.delete_arc(u, v, false);
        graph.delete_arc(v ^ 1, u ^ 1, false);
        v = u;
        if v == v0 {
            break;
        }
    }
}

/// `asg_bub_pop1_purge_graph`: the score-weighted bubble pop used by the
/// purge graph.
fn asg_bub_pop1_purge_graph(
    graph: &mut Asg,
    v0: u32,
    max_dist: u32,
    buffer: &mut PurgePopBuffer,
) -> u64 {
    if graph.sequences[(v0 >> 1) as usize].deleted
        || graph.sequences[(v0 >> 1) as usize].coverage_class == ALTER_LABEL
    {
        return 0;
    }
    if get_real_length(graph, v0) < 2 {
        return 0;
    }

    buffer.sources.clear();
    buffer.visited.clear();
    buffer.visited_edges.clear();
    {
        let info = &mut buffer.info[v0 as usize];
        info.score_h = 0;
        info.score_l = 0;
        info.distance = 0;
        info.remaining = 0;
    }
    buffer.sources.push(v0);

    let mut n_pending = 0u32;
    let mut n_tips = 0u32;
    let mut tip_end = u32::MAX;
    let mut success = false;

    'search: loop {
        let v = buffer.sources.pop().expect("the loop condition keeps it non-empty");
        let distance = buffer.info[v as usize].distance;
        let (start, count) = graph.arc_range(v);
        let total_score = decode_score(buffer.info[v as usize].score_h, buffer.info[v as usize].score_l);
        let mut early = false;

        let mut i = 0usize;
        while i < count {
            let arc = graph.arcs[start + i];
            i += 1;
            if arc.deleted {
                continue;
            }
            let w = arc.v;
            if (w >> 1) == (v0 >> 1) {
                early = true;
                break;
            }
            let arc_score = decode_score((arc.ul >> 32) as u32, arc.overlap);
            buffer.visited_edges.push((start + i - 1) as u32);
            if distance + 1 > max_dist {
                early = true;
                break;
            }

            let info = &mut buffer.info[w as usize];
            if !info.state {
                buffer.visited.push(w);
                info.parent = v;
                info.state = true;
                info.distance = distance + 1;
                let (h, l) = encode_score(total_score + arc_score);
                info.score_h = h;
                info.score_l = l;
                info.remaining = get_real_length(graph, w ^ 1);
                n_pending += 1;
            } else {
                let current = decode_score(info.score_h, info.score_l);
                if total_score + arc_score > current {
                    info.parent = v;
                    let (h, l) = encode_score(total_score + arc_score);
                    info.score_h = h;
                    info.score_l = l;
                }
                if distance + 1 < info.distance {
                    info.distance = distance + 1;
                }
            }

            info.remaining = info.remaining.wrapping_sub(1);
            if info.remaining == 0 {
                let outgoing = get_real_length(graph, w);
                if outgoing > 0 {
                    buffer.sources.push(w);
                } else {
                    if n_tips != 0 {
                        early = true;
                        break;
                    }
                    n_tips += 1;
                    tip_end = w;
                }
                n_pending = n_pending.wrapping_sub(1);
            }
        }

        if early {
            break 'search;
        }

        if n_tips == 1 {
            if tip_end != u32::MAX && n_pending == 0 && buffer.sources.is_empty() {
                buffer.sources.push(tip_end);
                success = true;
                break 'search;
            }
            break 'search;
        }

        if buffer.sources.is_empty() {
            break 'search;
        }
        if !(buffer.sources.len() > 1 || n_pending > 0) {
            success = true;
            break 'search;
        }
    }

    if success {
        purge_bub_backtrack_primary(graph, v0, buffer);
    }

    for &node in &buffer.visited {
        buffer.info[node as usize] = PurgeBinfo::default();
    }
    if success {
        1
    } else {
        0
    }
}

/// `asg_pop_bubble_purge_graph`.
fn asg_pop_bubble_purge_graph(graph: &mut Asg) -> u64 {
    if !graph.is_symmetric {
        graph.symmetrize();
    }
    let n_vtx = (graph.sequences.len() * 2) as u32;
    let mut buffer = PurgePopBuffer {
        info: vec![PurgeBinfo::default(); n_vtx as usize],
        ..PurgePopBuffer::default()
    };
    let mut n_pop = 0u64;
    for v in 0..n_vtx {
        let (_, count) = graph.arc_range(v);
        if count < 2
            || graph.sequences[(v >> 1) as usize].deleted
            || graph.sequences[(v >> 1) as usize].coverage_class == ALTER_LABEL
        {
            continue;
        }
        let live = graph.arcs_of(v).iter().filter(|arc| !arc.deleted).count();
        if live > 1 {
            n_pop += asg_bub_pop1_purge_graph(graph, v, graph.sequences.len() as u32, &mut buffer);
        }
    }
    if n_pop > 0 {
        graph.cleanup();
    }
    n_pop
}

/// `purge_g_arc_del_short_diploid_by_score`: drops the weak alternative arc of
/// every branching vertex.
fn purge_g_arc_del_short_diploid_by_score(graph: &mut Asg, drop_ratio: f32) -> u64 {
    let n_vtx = (graph.sequences.len() * 2) as u32;
    let mut scores: Vec<(u64, usize)> = Vec::new();
    for v in 0..n_vtx {
        if graph.sequences[(v >> 1) as usize].coverage_class == ALTER_LABEL
            || graph.sequences[(v >> 1) as usize].deleted
        {
            continue;
        }
        let (start, count) = graph.arc_range(v);
        if count < 2 {
            continue;
        }
        for i in 0..count {
            let arc = graph.arcs[start + i];
            let key = (u64::from((arc.ul >> 32) as u32) << 32) | u64::from(arc.overlap);
            scores.push((key, start + i));
        }
    }
    radix_sort_by_key(&mut scores, &|entry: &(u64, usize)| entry.0);

    let mut n_cut = 0u64;
    for &(_, eid) in &scores {
        let arc = graph.arcs[eid];
        let v = (arc.ul >> 32) as u32;
        let (start, count) = graph.arc_range(v);
        if count <= 1 {
            continue;
        }
        let mut kv = 0u32;
        let mut ovlp_max = 0i64;
        for i in 0..count {
            let other = graph.arcs[start + i];
            if other.deleted {
                continue;
            }
            let ovlp = decode_score((other.ul >> 32) as u32, other.overlap);
            if kv == 0 || ovlp_max < ovlp {
                ovlp_max = ovlp;
            }
            kv += 1;
        }
        if kv <= 1 {
            continue;
        }
        let ovlp = decode_score((arc.ul >> 32) as u32, arc.overlap);
        if ovlp >= 0 && ovlp_max >= 0 && (ovlp as f32) > (ovlp_max as f32) * drop_ratio {
            continue;
        }
        graph.arcs[eid].deleted = true;
        graph.delete_arc(arc.v ^ 1, ((arc.ul >> 32) as u32) ^ 1, true);
        n_cut += 1;
    }

    if n_cut > 0 {
        graph.cleanup();
        graph.symmetrize();
    }
    n_cut
}

/// `clean_purge_graph`.
fn clean_purge_graph(purge_g: &mut Asg, drop_ratio: f32, is_force_break: bool) {
    let mut operation = 1u64;
    while operation > 0 {
        operation = 0;
        operation += asg_pop_bubble_purge_graph(purge_g);
        operation += purge_g_arc_del_short_diploid_by_score(purge_g, drop_ratio);
    }
    if is_force_break {
        purge_g_arc_del_short_diploid_by_score(purge_g, 1.0);
    }
}

/// `fill_unitig`: measures every read-to-read step of a merged chain and
/// rewrites each buffer entry as `<vertex>:32 | <arc length>`.
fn fill_unitig(buffer: &mut [u64], read_graph: &Asg, edge: &[AsgArc]) -> u64 {
    let n = buffer.len();
    let mut total = 0u64;
    let mut i = 0usize;
    while i + 1 < n {
        let v = (buffer[i] >> 32) as u32;
        let w = (buffer[i + 1] >> 32) as u32;
        let mut length = 0u32;
        let mut found = false;
        for arc in read_graph.arcs_of(v) {
            if arc.deleted {
                continue;
            }
            if arc.v == w {
                length = arc.node_length();
                found = true;
                break;
            }
        }
        if !found {
            for arc in edge {
                if arc.deleted {
                    continue;
                }
                if (arc.ul >> 32) as u32 == v && arc.v == w {
                    length = arc.node_length();
                    found = true;
                    break;
                }
            }
            if !found {
                eprintln!("####ERROR1-fill: i: {}, v: {}, w: {}", i, v, w);
            }
        }
        buffer[i] = (u64::from(v) << 32) | u64::from(length);
        total += u64::from(length);
        i += 1;
    }

    if i < n {
        let v = (buffer[i] >> 32) as u32;
        let length = read_graph.sequences[(v >> 1) as usize].length;
        buffer[i] = (u64::from(v) << 32) | u64::from(length);
        total += u64::from(length);
    }
    total
}

/// `get_node_boundary_advance`: the forward and reverse arcs between two
/// consecutive unitigs of a merged chain.
#[allow(clippy::too_many_arguments)]
fn get_node_boundary_advance(
    read_graph: &Asg,
    reverse_sources: &[PafList],
    coverage_cut: &[CoverageCut],
    read_to_unitig: &ReadToUnitig,
    position_index: &[u64],
    max_hang: i32,
    max_hang_rate: f32,
    min_ovlp: i32,
    x_reads: &MaUtg,
    y_reads: &MaUtg,
    x_uid: u32,
    y_uid: u32,
    x_beg_index: i64,
    x_end_index: i64,
    y_beg_index: i64,
    y_end_index: i64,
    dir: u32,
    rev: u32,
    u_buffer: &mut Vec<ArcOffset>,
    tail_index: &mut Vec<i32>,
    prev_index: &mut Vec<i32>,
) -> Option<(AsgArc, AsgArc)> {
    u_buffer.clear();

    let mut k = x_beg_index;
    while k <= x_end_index {
        let x_read = (x_reads.reads[k as usize] >> 33) as usize;
        for hit in &reverse_sources[x_read].hits {
            let query_cut = &coverage_cut[hit.query_identifier as usize];
            let target_cut = &coverage_cut[hit.target_identifier as usize];
            if target_cut.deleted || read_graph.sequences[hit.target_identifier as usize].deleted {
                continue;
            }

            let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
            let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
            let (result, forward) = ma_hit_to_arc(
                hit,
                query_length,
                target_length,
                max_hang,
                max_hang_rate,
                min_ovlp,
            );
            if result < 0 {
                continue;
            }

            let hit_read = (forward.v >> 1) as usize;
            if read_graph.sequences[hit_read].deleted {
                continue;
            }
            let (unitig, is_unitig) = read_to_unitig.get(hit_read as u32);
            if is_unitig == 0 || unitig == u32::MAX || unitig != y_uid {
                continue;
            }

            let v = (x_reads.reads[k as usize] >> 32) as u32;
            let (unitig, is_unitig) = read_to_unitig.get(v >> 1);
            if is_unitig == 0 || unitig == u32::MAX || unitig != x_uid {
                continue;
            }
            // `(uint32_t)position_index[v>>1] != k`
            if position_index[(v >> 1) as usize] as u32 != k as u32 {
                continue;
            }

            let w = (y_reads.reads[position_index[hit_read] as u32 as usize] >> 32) as u32;
            let v_dir = u32::from((forward.ul >> 32) as u32 == v);
            let w_dir = u32::from(forward.v == w);
            if rev == 0 && v_dir != w_dir {
                continue;
            }
            if rev == 1 && v_dir == w_dir {
                continue;
            }
            if dir == v_dir {
                continue;
            }

            let offset = position_index[hit_read] as u32;
            if i64::from(offset) < y_beg_index || i64::from(offset) > y_end_index {
                continue;
            }

            let Some(index) =
                get_specific_overlap(&reverse_sources[hit.target_identifier as usize].hits, hit.target_identifier, hit.query_identifier)
            else {
                continue;
            };
            let reverse_hit = &reverse_sources[hit.target_identifier as usize].hits[index];
            let query_cut = &coverage_cut[reverse_hit.query_identifier as usize];
            let target_cut = &coverage_cut[reverse_hit.target_identifier as usize];
            if target_cut.deleted
                || read_graph.sequences[reverse_hit.target_identifier as usize].deleted
            {
                continue;
            }
            let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
            let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
            let (result, _backward) = ma_hit_to_arc(
                reverse_hit,
                query_length,
                target_length,
                max_hang,
                max_hang_rate,
                min_ovlp,
            );
            if result < 0 {
                continue;
            }

            let (tmp, reverse_flag) =
                get_xy_pos(read_graph, &forward, v, w, x_reads.length, y_reads.length, position_index);
            if (tmp >> 32) == u64::from(u32::MAX) || (tmp as u32) == u32::MAX {
                continue;
            }
            let mut arc = forward;
            arc.exact_region = reverse_flag;
            u_buffer.push(ArcOffset {
                arc,
                offset: tmp,
                weight: 1,
            });
        }
        k += 1;
    }

    if u_buffer.is_empty() {
        return None;
    }

    u_buffer.sort_by(crate::assembly::trans_cov::compare_arc_offsets);
    let mut m = 1usize;
    for k in 1..u_buffer.len() {
        if u_buffer[m - 1].offset == u_buffer[k].offset {
            u_buffer[m - 1].weight += u_buffer[k].weight;
            if u_buffer[k].arc.overlap > u_buffer[m - 1].arc.overlap {
                u_buffer[m - 1].arc = u_buffer[k].arc;
            }
            continue;
        }
        u_buffer[m] = u_buffer[k];
        m += 1;
    }
    u_buffer.truncate(m);

    quick_lis(u_buffer, tail_index, prev_index);
    if tail_index.is_empty() {
        return None;
    }

    let n_tail = tail_index.len();
    let mut found = None;
    for j in 0..n_tail {
        let k = if dir == 0 { j } else { n_tail - 1 - j };
        let entry = u_buffer[tail_index[k] as usize];
        let v = entry.arc.v >> 1;
        let w = (entry.arc.ul >> 33) as u32;
        let Some(index) = get_specific_overlap(&reverse_sources[v as usize].hits, v, w) else {
            continue;
        };
        let hit = &reverse_sources[v as usize].hits[index];
        let query_cut = &coverage_cut[hit.query_identifier as usize];
        let target_cut = &coverage_cut[hit.target_identifier as usize];
        if target_cut.deleted || read_graph.sequences[hit.target_identifier as usize].deleted {
            continue;
        }
        let query_length = query_cut.end.wrapping_sub(query_cut.start) as i32;
        let target_length = target_cut.end.wrapping_sub(target_cut.start) as i32;
        let (result, backward) = ma_hit_to_arc(
            hit,
            query_length,
            target_length,
            max_hang,
            max_hang_rate,
            min_ovlp,
        );
        if result < 0 {
            continue;
        }
        found = Some((entry.arc, backward));
        break;
    }
    found
}

/// `purge_merge`: concatenates the unitigs of one chain, collecting the
/// connecting arcs and rewriting the kept unitig's read list.
#[allow(clippy::too_many_arguments)]
fn purge_merge(
    purge_g: &mut Asg,
    unitig_graph: &mut UnitigGraph,
    all_ovlp: &[Vec<HapOverlap>],
    b_0: &[u32],
    read_graph: &Asg,
    coverage_cut: &[CoverageCut],
    reverse_sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    max_hang: i32,
    max_hang_rate: f32,
    min_ovlp: i32,
    edge: &mut Vec<AsgArc>,
    visit: &mut [u8],
    u_buffer: &mut Vec<ArcOffset>,
    tail_index: &mut Vec<i32>,
    prev_index: &mut Vec<i32>,
    cov: &mut HapCov,
) {
    let mut buffer: Vec<u64> = Vec::new();
    let mut i = 0usize;
    while i < b_0.len() {
        let mut cut_beg: u32 = 0;
        let mut cut_end: u32 = u32::MAX;
        buffer.clear();

        let v = b_0[i];
        let keep_uid = v >> 1;
        {
            let x_reads = &unitig_graph.unitigs[keep_uid as usize];
            if v & 1 != 0 {
                for k in (0..x_reads.reads.len()).rev() {
                    buffer.push(x_reads.reads[k] ^ 0x1_0000_0000);
                }
            } else {
                buffer.extend_from_slice(&x_reads.reads);
            }
            cut_beg = 0;
            cut_end = x_reads.reads.len().wrapping_sub(1) as u32;
        }
        i += 1;

        while i < b_0.len() {
            let v = b_0[i - 1];
            let w = b_0[i];
            let list = &all_ovlp[(v >> 1) as usize];
            let index = get_specific_hap_overlap(list, v >> 1, w >> 1)
                .expect("purge_merge needs an overlap for every chain step");
            let x = list[index];

            let x_reads_len = unitig_graph.unitigs[(v >> 1) as usize].reads.len() as u32;
            let y_reads_len = unitig_graph.unitigs[(w >> 1) as usize].reads.len() as u32;

            let mut beg_index = x.x_beg_id;
            if cut_beg > beg_index {
                beg_index = cut_beg;
            }
            let mut end_index = x.x_end_id.wrapping_sub(1);
            if cut_end < end_index {
                end_index = cut_end;
            }

            let (t_forward, t_backward) = get_node_boundary_advance(
                read_graph,
                reverse_sources,
                coverage_cut,
                read_to_unitig,
                &cov.pos_idx,
                max_hang,
                max_hang_rate,
                min_ovlp,
                &unitig_graph.unitigs[(v >> 1) as usize],
                &unitig_graph.unitigs[(w >> 1) as usize],
                v >> 1,
                w >> 1,
                i64::from(beg_index),
                i64::from(end_index),
                i64::from(x.y_beg_id),
                i64::from(x.y_end_id.wrapping_sub(1)),
                v & 1,
                u32::from(x.reverse),
                u_buffer,
                tail_index,
                prev_index,
            )
            .unwrap_or((AsgArc { deleted: true, ..AsgArc::default() }, AsgArc { deleted: true, ..AsgArc::default() }));

            if t_forward.deleted || t_backward.deleted {
                break;
            }
            edge.push(t_forward);
            edge.push(t_backward);

            let mut x_beg_index = 0u32;
            let mut x_end_index = x_reads_len - 1;
            let mut y_beg_index = 0u32;
            let mut y_end_index = y_reads_len - 1;

            if v & 1 == 0 {
                x_end_index = cov.pos_idx[(t_forward.ul >> 33) as usize] as u32;
                let drop = cut_end.wrapping_sub(x_end_index) as usize;
                let len = buffer.len().wrapping_sub(drop);
                buffer.truncate(len);
            } else {
                x_beg_index = cov.pos_idx[(t_forward.ul >> 33) as usize] as u32;
                let drop = x_beg_index.wrapping_sub(cut_beg) as usize;
                let len = buffer.len().wrapping_sub(drop);
                buffer.truncate(len);
            }

            if w & 1 == 1 {
                y_end_index = cov.pos_idx[(t_forward.v >> 1) as usize] as u32;
            } else {
                y_beg_index = cov.pos_idx[(t_forward.v >> 1) as usize] as u32;
            }
            cut_beg = y_beg_index;
            cut_end = y_end_index;

            {
                let y_reads = &unitig_graph.unitigs[(w >> 1) as usize];
                if w & 1 == 1 {
                    for k in (y_beg_index..=y_end_index).rev() {
                        buffer.push(y_reads.reads[k as usize] ^ 0x1_0000_0000);
                    }
                } else {
                    for k in y_beg_index..=y_end_index {
                        buffer.push(y_reads.reads[k as usize]);
                    }
                }
            }

            purge_g.sequences[(w >> 1) as usize].coverage_class = ALTER_LABEL;
            collect_trans_purge_joint_cov(&unitig_graph.unitigs, cov, &x);
            i += 1;
        }

        let total_len = fill_unitig(&mut buffer, read_graph, edge);

        let x_reads = &mut unitig_graph.unitigs[keep_uid as usize];
        x_reads.reads = buffer.clone();
        x_reads.length = total_len as u32;
        x_reads.circular = false;

        if let Some(&first_packed) = x_reads.reads.first() {
            let start = (first_packed >> 32) as u32;
            if x_reads.start != start {
                x_reads.start = start;
                let v = (keep_uid << 1) + 1;
                delete_vertex_arcs(&mut unitig_graph.graph, v);
            }
        }
        if let Some(&last_packed) = x_reads.reads.last() {
            let end = (last_packed >> 32) as u32 ^ 1;
            if x_reads.end != end {
                x_reads.end = end;
                let v = keep_uid << 1;
                delete_vertex_arcs(&mut unitig_graph.graph, v);
            }
        }
    }

    for &v in b_0 {
        visit[(v >> 1) as usize] = 1;
        // The reference indexes `seq[v>1]` here, not `seq[v>>1]`; both are the
        // same expression for every `v` except 0 and 1, where it picks seq[1]
        // and seq[0]. Replicated verbatim.
        let index = usize::from(v > 1);
        if purge_g.sequences[index].coverage_class != ALTER_LABEL {
            continue;
        }
        asg_seq_drop(purge_g, index as u32);
    }
}

/// Deletes every live outgoing arc of `v` and its reverse complement.
fn delete_vertex_arcs(graph: &mut Asg, v: u32) {
    let targets: Vec<u32> = graph
        .arcs_of(v)
        .iter()
        .filter(|arc| !arc.deleted)
        .map(|arc| arc.v)
        .collect();
    for target in targets {
        graph.delete_arc(v, target, true);
        graph.delete_arc(target ^ 1, v ^ 1, true);
    }
}

/// `link_unitigs`: walks every maximal chain of the purge graph and merges its
/// unitigs into the first one.
#[allow(clippy::too_many_arguments)]
fn link_unitigs(
    purge_g: &mut Asg,
    unitig_graph: &mut UnitigGraph,
    read_graph: &Asg,
    all_ovlp: &[Vec<HapOverlap>],
    coverage_cut: &[CoverageCut],
    reverse_sources: &[PafList],
    read_to_unitig: &ReadToUnitig,
    max_hang: i32,
    max_hang_rate: f32,
    min_ovlp: i32,
    edge: &mut Vec<AsgArc>,
    visit: &mut [u8],
    cov: &mut HapCov,
) {
    let n_vtx = (purge_g.sequences.len() * 2) as u32;
    visit.fill(0);

    let mut b_0: Vec<u32> = Vec::new();
    let mut u_buffer: Vec<ArcOffset> = Vec::new();
    let mut tail_index: Vec<i32> = Vec::new();
    let mut prev_index: Vec<i32> = Vec::new();

    for v in 0..n_vtx {
        if purge_g.sequences[(v >> 1) as usize].coverage_class == ALTER_LABEL
            || purge_g.sequences[(v >> 1) as usize].deleted
            || visit[(v >> 1) as usize] != 0
        {
            continue;
        }
        if get_real_length(purge_g, v) != 1 {
            continue;
        }
        if get_real_length(purge_g, v ^ 1) != 0 {
            continue;
        }

        b_0.clear();
        let mut end = u32::MAX;
        let mut node_len = 0i64;
        let mut base_len = 0i64;
        let mut max_stop_node_len = 0i64;
        let mut max_stop_base_len = 0i64;
        if get_unitig(
            purge_g,
            None,
            v,
            &mut end,
            &mut node_len,
            &mut base_len,
            &mut max_stop_node_len,
            &mut max_stop_base_len,
            1,
            Some(&mut b_0),
        ) == LOOP
        {
            continue;
        }

        purge_merge(
            purge_g,
            unitig_graph,
            all_ovlp,
            &b_0,
            read_graph,
            coverage_cut,
            reverse_sources,
            read_to_unitig,
            max_hang,
            max_hang_rate,
            min_ovlp,
            edge,
            visit,
            &mut u_buffer,
            &mut tail_index,
            &mut prev_index,
            cov,
        );
    }
}

/// `purge_dups`: the full duplicate-purging pass.
///
/// `just_coverage` selects the coverage-only exit used by
/// [`set_hom_global_coverage`]; the output stage calls it with
/// `just_coverage = false` to purge the copied primary graph.
///
/// `density` is the reference's `Hap_rate`; `max_hang` and `min_ovlp` are the
/// alignment bounds. `new_rtg_edges` collects the read-graph arcs the merge
/// pass creates, which the caller's `renew_utg` consumes.
///
/// # Errors
///
/// Returns an error when the ploidy fallback (`if_ploid_sample`) is needed,
/// which only happens when neither coverage peak exists.
#[allow(clippy::too_many_arguments)]
pub fn purge_dups(
    unitig_graph: &mut UnitigGraph,
    read_graph: &Asg,
    coverage_cut: &[CoverageCut],
    sources: &[PafList],
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
    new_rtg_edges: &mut Vec<AsgArc>,
    density: f32,
    purge_min_len: u32,
    max_hang: i32,
    min_ovlp: i32,
    drop_ratio: f32,
    just_contain: bool,
    just_coverage: bool,
    cov: &mut HapCov,
    collect_p_trans: bool,
    collect_p_trans_f: bool,
    options: &mut AssemblyOptions,
    trio_flags: &mut [u8],
) -> Result<(), String> {

    let read_count = read_graph.sequences.len();
    if cov.pos_idx.len() < read_count {
        cov.pos_idx.resize(read_count, u64::MAX);
    }
    cov.pos_idx[..read_count].fill(u64::MAX);

    let unitig_count = unitig_graph.graph.sequences.len();
    let mut all_ovlp = HapOverlapsList::new(unitig_count);
    let mut back_all_ovlp = HapOverlapsList::new(unitig_count);

    // `hap_buf.cov_threshold`. On this path `hom_global_coverage` was already
    // inferred by `set_hom_global_coverage`, so the fallback is only a
    // consistency check.
    let mut cov_threshold: i64 = if options.global_homozygous_coverage != -1 {
        if options.global_homozygous_coverage_is_set != 0 {
            (f64::from(options.global_homozygous_coverage) * HOM_PEAK_RATE) as i64
        } else {
            i64::from(options.global_homozygous_coverage)
        }
    } else {
        get_read_coverage_thres(
            unitig_graph,
            read_graph,
            read_to_unitig,
            &mut cov.pos_idx,
            sources,
            coverage_cut,
            options,
        )
        .0
    };

    if cov_threshold < 0 {
        return Err(
            "purge_dups needs if_ploid_sample (Purge_Dups.cpp:4800), which is not ported yet"
                .to_string(),
        );
    }

    // `init_hap_alignment_struct_pip` would build `hh` here; only the
    // candidate scan reads it.
    if options.global_homozygous_coverage == -1 {
        options.global_homozygous_coverage = cov_threshold as i32;
    }
    if options.global_purge_coverage != -1 {
        cov_threshold = i64::from(options.global_purge_coverage);
    }

    eprintln!(
        "[M::purge_dups] homozygous read coverage threshold: {}",
        if options.global_homozygous_coverage_is_set != 0 {
            options.global_homozygous_coverage
        } else {
            (f64::from(options.global_homozygous_coverage) / HOM_PEAK_RATE) as i32
        }
    );
    eprintln!(
        "[M::purge_dups] purge duplication coverage threshold: {}",
        cov_threshold
    );

    // Map every read of every unitig to it and store its position, even when
    // only the coverage is wanted: `end_coverage` then clears the unitig
    // entries again.
    for v in 0..unitig_graph.graph.sequences.len() {
        let unitig_id = v as u32;
        let mut offset: u32 = 0;
        let unitig = &unitig_graph.unitigs[v];
        for (i, entry) in unitig.reads.iter().enumerate() {
            let read = (entry >> 33) as u32;
            let class = read_graph.sequences[read as usize].coverage_class;
            read_to_unitig.set(read, unitig_id, true, Some(class));
            cov.pos_idx[read as usize] = (u64::from(offset) << 32) | i as u64;
            offset = offset.wrapping_add(*entry as u32);
        }
    }

    if !just_coverage {
        hap_alignment_advance(
            unitig_graph,
            read_graph,
            coverage_cut,
            sources,
            reverse_sources,
            read_to_unitig,
            &cov.pos_idx,
            density,
            max_hang,
            min_ovlp,
            cov_threshold,
            cov,
            options,
            &mut all_ovlp,
        )?;
        if std::env::var("RSHIFIASM_PURGE_DEBUG").is_ok() {
            let total: usize = all_ovlp.x.iter().map(|l| l.len()).sum();
            let nonempty = all_ovlp.x.iter().filter(|l| !l.is_empty()).count();
            eprintln!("[M::purge_dups_debug] after scan: total overlaps {} across {} unitig lists (pre-filter)", total, nonempty);
            for (_uid, list) in all_ovlp.x.iter().enumerate() {
                for o in list {
                    eprintln!("[M::purge_dups_debug]   ovlp x={} y={} xid=[{},{}] yid=[{},{}] xpos=[{},{}] ypos=[{},{}] rev={} kind={} w={} sc={}",
                        o.x_uid, o.y_uid, o.x_beg_id, o.x_end_id, o.y_beg_id, o.y_end_id,
                        o.x_beg_pos, o.x_end_pos, o.y_beg_pos, o.y_end_pos,
                        o.reverse, o.kind, o.weight, o.score);
                }
            }
        }
        filter_hap_overlaps_by_length(&mut all_ovlp, purge_min_len);
        if std::env::var("RSHIFIASM_PURGE_DEBUG").is_ok() {
            let total: usize = all_ovlp.x.iter().map(|l| l.len()).sum();
            eprintln!("[M::purge_dups_debug] after filter(min_len={}): total overlaps {}", purge_min_len, total);
        }

        let mut pg = init_p_g_t(unitig_graph, read_graph, &cov.is_r_het);
        normalize_hap_overlaps_advance_by_p_g_t(
            &mut all_ovlp.x,
            &mut back_all_ovlp.x,
            &pg,
            options.polyploidy,
            0.8,
        );

        let align_ctx = PurgeAlignContext {
            unitigs: &unitig_graph.unitigs,
            read_graph,
            reverse_sources,
            read_to_unitig,
            coverage_cut,
            max_hang,
            max_hang_rate: options.maximum_hang_rate,
            min_ovlp,
            purge_level: options.purge_level_primary as u32,
        };

        if collect_p_trans && !collect_p_trans_f {
            collect_purge_trans_cov(&all_ovlp.x, &align_ctx, cov);
        }

        if options.polyploidy <= 2 {
            let partition = options.flags & FLAG_PARTITION != 0 && cov.t_ch.is_some();
            crate::assembly::maxcut::mc_solve_purge(
                unitig_graph,
                &mut all_ovlp.x,
                &cov.is_r_het,
                &cov.read_lengths,
                0.8,
                options.polyploidy,
                options.perturbation_round_count,
                options.perturbation_fraction,
                options.random_seed,
                if partition { Some(trio_flags) } else { None },
            );
        }

        if collect_p_trans && collect_p_trans_f {
            collect_purge_trans_cov(&all_ovlp.x, &align_ctx, cov);
        }

        remove_contained_haplotig(&mut all_ovlp.x, unitig_graph, &mut pg.pg_h_lev, cov);

        if !just_contain {
            for u_id in 0..all_ovlp.x.len() {
                if pg.pg_h_lev.sequences[u_id].deleted
                    || pg.pg_h_lev.sequences[u_id].coverage_class == ALTER_LABEL
                {
                    continue;
                }
                for overlap in &all_ovlp.x[u_id] {
                    if overlap.status == DELETE {
                        continue;
                    }
                    let x_uid = overlap.x_uid as usize;
                    let y_uid = overlap.y_uid as usize;
                    if pg.pg_h_lev.sequences[x_uid].coverage_class == ALTER_LABEL
                        || pg.pg_h_lev.sequences[x_uid].deleted
                        || pg.pg_h_lev.sequences[y_uid].coverage_class == ALTER_LABEL
                        || pg.pg_h_lev.sequences[y_uid].deleted
                    {
                        continue;
                    }

                    let x_len = unitig_graph.unitigs[x_uid].length;
                    let y_len = unitig_graph.unitigs[y_uid].length;
                    let Some(arc) = get_hap_arch(
                        overlap,
                        x_len,
                        y_len,
                        max_hang,
                        options.maximum_hang_rate,
                        min_ovlp,
                    ) else {
                        continue;
                    };
                    pg.pg_h_lev.push_arc(arc);
                }
            }

            pg.pg_h_lev.cleanup();
            pg.pg_h_lev.symmetrize();
            clean_purge_graph(&mut pg.pg_h_lev, drop_ratio, true);

            let mut visit = vec![0u8; pg.pg_h_lev.sequences.len()];
            link_unitigs(
                &mut pg.pg_h_lev,
                unitig_graph,
                read_graph,
                &all_ovlp.x,
                coverage_cut,
                reverse_sources,
                read_to_unitig,
                max_hang,
                options.maximum_hang_rate,
                min_ovlp,
                new_rtg_edges,
                &mut visit,
                cov,
            );
        }

        for v in 0..unitig_count {
            if pg.pg_h_lev.sequences[v].coverage_class == ALTER_LABEL {
                unitig_graph.graph.sequences[v].coverage_class = ALTER_LABEL;
            }
        }
    }

    // `end_coverage`.
    read_to_unitig.clear_unitig_entries();
    unitig_graph.graph.cleanup();
    cov.pos_idx[..read_count].fill(u64::MAX);

    Ok(())
}

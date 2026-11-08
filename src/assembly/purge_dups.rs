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
//! The driver, the vote-counting half of `hap_alignment_advance_worker`, the
//! purge-graph bookkeeping and the partitioning side effects are ported. The
//! haplotype-alignment candidate scan and the non-empty purge-graph paths are
//! unreachable on the chr11-2M corpus (the copied primary graph holds a single
//! unitig) and return an explicit error instead of silently diverging.

use crate::assembly::asg::Asg;
use crate::assembly::hap_cov::{HapCov, N_HET, P_HET};
use crate::assembly::hit::{CoverageCut, ReadToUnitig, ALTER_LABEL, PRIMARY_LABEL};
use crate::assembly::unitig::{ma_ug_gen, UnitigGraph, FATHER};
use crate::command_line::options::{AssemblyOptions, FLAG_PARTITION};
use crate::correction::ec_driver::PafList;
use crate::counting::peak_analysis::adjust_homozygous_peak;
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

/// `normalize_hap_overlaps_advance_by_p_g_t`.
///
/// Every step of the reference walks the haplotype-overlap lists; with no
/// overlaps the function is a no-op.
fn normalize_hap_overlaps_advance_by_p_g_t(
    all_ovlp: &HapOverlapsList,
    back_all_ovlp: &HapOverlapsList,
    _pg: &PG,
) -> Result<(), String> {
    if all_ovlp.is_empty() && back_all_ovlp.is_empty() {
        return Ok(());
    }
    Err("purge_dups: hap-overlap normalization is not ported yet".to_string())
}

/// `collect_purge_trans_cov`.
///
/// The reference walks the hap-overlap lists and records their transitions in
/// `cov->t_ch->k_trans`; with no overlaps there is nothing to record.
fn collect_purge_trans_cov(
    all_ovlp: &HapOverlapsList,
    _position_index: &[u64],
) -> Result<(), String> {
    if all_ovlp.is_empty() {
        return Ok(());
    }
    Err("purge_dups: collect_purge_trans_cov is not ported yet".to_string())
}

/// `remove_contained_haplotig`.
///
/// With no hap overlaps no unitig is marked alternate.
fn remove_contained_haplotig(
    all_ovlp: &HapOverlapsList,
    _unitig_graph: &mut UnitigGraph,
    _pg: &mut PG,
) -> Result<(), String> {
    if all_ovlp.is_empty() {
        return Ok(());
    }
    Err("purge_dups: remove_contained_haplotig is not ported yet".to_string())
}

/// The partitioning side effect of `mc_solve` on an edgeless min-cut graph.
///
/// `init_mc_g_t` zeroes the solution vector, `update_mc_edges` builds no
/// edges, `mb_solve_core` returns immediately without a reference and
/// `mc_solve_core_adv` leaves every isolated component at `s = 0`. The only
/// observable effect left is `p_nodes`: every read of every unitig is marked
/// `SET_TRIO`, and the reads already flagged `P_HET` are also marked
/// `FATHER`.
fn mc_solve_edgeless(
    unitig_graph: &UnitigGraph,
    is_r_het: &[u8],
    trio_flags: &mut [u8],
    collect_trans: bool,
) {
    eprintln!("[M::mc_solve:: # edges: 0]");
    eprintln!("[M::mc_solve_core_adv::0.000] ==> Partition");

    if !collect_trans {
        return;
    }

    for unitig in &unitig_graph.unitigs {
        for &packed in &unitig.reads {
            let read = (packed >> 33) as usize;
            trio_flags[read] |= SET_TRIO;
            if is_r_het[read] == N_HET {
                continue;
            }
            // `set_p_flag` receives `s == 0` for every unitig here, so the
            // `s > 0 ? FATHER : MOTHER` branch is unreachable.
            if is_r_het[read] & P_HET != 0 {
                trio_flags[read] |= FATHER;
            }
        }
    }
}

/// `clean_purge_graph`.
fn clean_purge_graph(purge_g: &Asg, _drop_ratio: f32) -> Result<(), String> {
    if purge_g.arcs.iter().all(|arc| arc.deleted) {
        return Ok(());
    }
    Err("purge_dups: clean_purge_graph is not ported yet".to_string())
}

/// `link_unitigs`.
fn link_unitigs(purge_g: &Asg) -> Result<(), String> {
    if purge_g.arcs.iter().all(|arc| arc.deleted) {
        return Ok(());
    }
    Err("purge_dups: link_unitigs is not ported yet".to_string())
}

/// `purge_dups`: the full duplicate-purging pass.
///
/// `just_coverage` selects the coverage-only exit used by
/// [`set_hom_global_coverage`]; the output stage calls it with
/// `just_coverage = false` to purge the copied primary graph.
///
/// `density`, `max_hang` and `min_ovlp` are the reference's `Hap_rate`,
/// `max_hang` and `min_ovlp`; they are consumed by the haplotype-alignment
/// scan, which is not ported yet (see [`hap_alignment_advance_worker`]).
///
/// # Errors
///
/// Returns an error when the corpus reaches a code path that is not ported:
/// the ploidy fallback (`if_ploid_sample`), the candidate scan, a non-empty
/// hap-overlap list, or a purge graph with live arcs.
#[allow(clippy::too_many_arguments)]
pub fn purge_dups(
    unitig_graph: &mut UnitigGraph,
    read_graph: &Asg,
    coverage_cut: &[CoverageCut],
    sources: &[PafList],
    reverse_sources: &[PafList],
    read_to_unitig: &mut ReadToUnitig,
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
    let back_all_ovlp = HapOverlapsList::new(unitig_count);

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
            for (uid, list) in all_ovlp.x.iter().enumerate() {
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
        normalize_hap_overlaps_advance_by_p_g_t(&all_ovlp, &back_all_ovlp, &pg)?;

        if collect_p_trans && !collect_p_trans_f {
            collect_purge_trans_cov(&all_ovlp, &cov.pos_idx)?;
        }

        if options.polyploidy <= 2 {
            mc_solve_edgeless(
                unitig_graph,
                &cov.is_r_het,
                trio_flags,
                options.flags & FLAG_PARTITION != 0 && cov.t_ch.is_some(),
            );
        }

        if collect_p_trans && collect_p_trans_f {
            collect_purge_trans_cov(&all_ovlp, &cov.pos_idx)?;
        }

        remove_contained_haplotig(&all_ovlp, unitig_graph, &mut pg)?;

        if !just_contain {
            // `get_hap_arch` loop: no arcs because `all_ovlp` is empty.
            pg.pg_h_lev.cleanup();
            pg.pg_h_lev.symmetrize();
            clean_purge_graph(&pg.pg_h_lev, drop_ratio)?;
            link_unitigs(&pg.pg_h_lev)?;
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

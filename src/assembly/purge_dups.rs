//! Coverage inference for the purge stage, the counterpart of the peak
//! helpers in `Purge_Dups.cpp` and of `set_hom_global_coverage`.
//!
//! The cleaning stage calls `set_hom_global_coverage` once before the bubble
//! rescue round; it builds a unitig graph and runs `purge_dups` with
//! `just_coverage = 1`, which only infers the homozygous coverage from the
//! read overlap histogram and the per-unitig coverage bins, then exits through
//! `end_coverage`. The hap-alignment machinery of the full purge
//! (`hap_alignment_advance_worker`, `mc_solve`, `remove_contained_haplotig`,
//! ...) is never reached on that path and is not modelled here.

use crate::assembly::asg::Asg;
use crate::assembly::hit::{CoverageCut, ReadToUnitig, ALTER_LABEL};
use crate::assembly::unitig::{ma_ug_gen, UnitigGraph};
use crate::command_line::options::AssemblyOptions;
use crate::correction::ec_driver::PafList;
use crate::counting::peak_analysis::adjust_homozygous_peak;

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

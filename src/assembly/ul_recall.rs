//! The rescall/renewal drivers of the `--ul` second pass (`inter.cpp`
//! 20583-21600): `gen_UL_reovlps` rebuilds the UL overlaps against the
//! polished graph, `rescall_ul_pipeline` re-aligns every UL read with
//! `worker_for_ul_rescall_alignment`, and `ul_realignment` wraps the whole
//! renewal entry with the `.re.*` caches.

use crate::assembly::ul_binning::{build_binning_index, graph_ovlp_binning};
use crate::assembly::ul_glchain::{gl_chain, gl_chain_flter, ul_lalign, GdpChain};
use crate::assembly::ul_index::{
    build_unitig_minimizer_index, gen_ul_idx, uidx_load, uidx_write, UlIndex,
    UnitigMinimizerIndex,
};
use crate::assembly::ul_load::P_CHAIN_RATE;
use crate::assembly::ul_mask::{load_emask_table, write_emask_table, REA_ALIGN_CUTOFF};
use crate::assembly::ul_phase::{
    update_sketch_trace, EmaskTable, MaskUlOv, SecondaryWindows, UL_TOPN,
};
use crate::assembly::ul_renew::{
    clear_all_ul_t, filter_ul_ug, gen_polished_ug, update_ug_arch_ul_mul,
};
use crate::assembly::ul_rescall::{gen_shared_intervals, UlRescallContext};
use crate::assembly::ul_scall::ScallContext;
use crate::assembly::ul_store::{
    append_ul_t, gen_ul_vec_rid_t, load_all_ul_t, write_all_ul_t, AllUl, AppendParameters,
};
use crate::assembly::ul::UG_SKIP_GRAPH_N;
use crate::command_line::options::{AssemblyOptions, FLAG_NO_HPC};
use crate::correction::hc_align::BitExtz;
use crate::correction::haplotype_evidence::HaplotypeEvidenceAlloc;
use crate::correction::interval::IntervalBuffer;
use crate::correction::target::TargetSource;
use crate::correction::ul_align::{gen_aux_ovlp, ul_map_lchain_khit, UlAlignParameters};
use crate::logging;
use crate::overlap::paf::PafHit;
use crate::overlap::region::{ChainData, OverlapRegion};
use crate::overlap::region::CandidatesList;
use crate::reads::read_store::ReadStore;
use crate::sequence::fastx_reader::FastxReader;

/// `CHUNK_SIZE` (`inter.cpp` 59).
const CHUNK_SIZE: u64 = 1_000_000_000;

/// The `mg_idxopt_t` + `ug_opt_t` slices of the rescall pass.
#[derive(Debug, Clone, Copy)]
pub struct RescallParameters {
    pub kmer_length: u32,
    pub minimizer_window: u32,
    pub error_rate: f64,
    pub maximum_chain_count: u32,
    pub maximum_hang_rate: f32,
    pub maximum_hang: i64,
    pub minimum_overlap: i64,
    pub homopolymer_compression: bool,
    pub sample_distance: i32,
    pub rewindow_size: i32,
}

/// `init_mg_opt(..., 19, 10, cutoff, ...)` of every UL driver.
impl RescallParameters {
    pub fn new(options: &AssemblyOptions, cutoff: i64) -> Self {
        let _ = cutoff;
        Self {
            kmer_length: 19,
            minimizer_window: 10,
            error_rate: options.ultra_long_error_rate,
            maximum_chain_count: options.maximum_chain_count as u32,
            maximum_hang_rate: options.maximum_hang_rate,
            maximum_hang: i64::from(options.maximum_hang_length),
            minimum_overlap: i64::from(options.minimum_overlap_length),
            homopolymer_compression: options.flags & FLAG_NO_HPC == 0,
            sample_distance: options.minimizer_sample_distance,
            rewindow_size: options.minimizer_rewindow_size,
        }
    }
}

/// The per-read worker buffers (`ha_ovec_buf_t` + `glchain_t` + `gdpchain_t`).
#[derive(Default)]
pub struct RescallBuffers {
    pub overlaps: Vec<OverlapRegion>,
    pub candidates: CandidatesList,
    pub extz: BitExtz,
    pub dummy: IntervalBuffer,
    pub auxiliary: OverlapRegion,
    /// `b->r_buf.a` (`v_idx`).
    pub r_buf: Vec<u64>,
    /// `b->hap`.
    pub hap: HaplotypeEvidenceAlloc,
    /// `bl->lo`.
    pub lo: Vec<crate::assembly::ul::UlOv>,
    /// `bl->tk`.
    pub tk: Vec<crate::assembly::ul::UlOv>,
    /// `bl->srt.a`.
    pub ll_srt: Vec<u64>,
    /// `s->sps[tid]`.
    pub sps: Vec<u64>,
    /// `s->mk[tid]`.
    pub mk: MaskUlOv,
    /// `s->gdp[tid]`.
    pub gdp: GdpChain,
    /// `b->clist.chainDP`.
    pub chain_data: ChainData,
    /// `num_correct_base`.
    pub num_correct_base: u64,
    /// `num_read_base`.
    pub num_read_base: u64,
}

fn single_window_list(region: &OverlapRegion, error: u32) -> crate::correction::window::WindowList {
    crate::correction::window::WindowList {
        x_start: region.x_position_start as i32,
        x_end: region.x_position_end as i32 + 1,
        y_start: 0,
        y_end: 0,
        extra_begin: 0,
        extra_end: 0,
        error: 0,
        error_threshold: 0,
        cigar_start: 0,
        cigar_length: error,
    }
}

/// `worker_for_ul_rescall_alignment` (`inter.cpp` 9299) for one read.
#[allow(clippy::too_many_arguments)]
pub fn rescall_one_read(
    map_context: &ScallContext<'_>,
    context: &UlRescallContext<'_>,
    buffers: &mut RescallBuffers,
    store: &mut AllUl,
    read_id: usize,
    sequence: &[u8],
    parameters: &RescallParameters,
    emask: &EmaskTable,
    slots: &mut SecondaryWindows,
) {
    let diff_ec_ul = parameters.error_rate;
    let win_len = ((31.0 / diff_ec_ul) as i64).min(375);
    let targets = TargetSource::Unitigs {
        unitigs: &context.unitig_graph.unitigs,
        reads: context.reads,
    };
    let align_parameters = UlAlignParameters {
        kmer_length: parameters.kmer_length,
        minimizer_window: parameters.minimizer_window,
        error_rate: diff_ec_ul,
        band_width_threshold: diff_ec_ul,
        maximum_chain_count: parameters.maximum_chain_count,
        homopolymer_compression: parameters.homopolymer_compression,
        sample_distance: map_context.sample_distance as u32,
        rewindow_size: map_context.rewindow_size as u32,
    };

    // `ul_map_lchain(..., 2, 1)` — `mcopy_khit_cut = 1` for the rescall.
    let mapped = ul_map_lchain_khit(
        sequence,
        &map_context.unitig_index.index,
        &targets,
        &align_parameters,
        Some(&map_context.unitig_index.filter),
        1,
    );
    buffers.overlaps = mapped.overlaps;
    buffers.candidates = mapped.candidates;

    // first `ul_lalign` pass with `aln = NULL`
    ul_lalign(
        context,
        &mut buffers.overlaps,
        &mut buffers.candidates,
        &targets,
        sequence,
        &mut buffers.extz,
        &mut buffers.dummy,
        &mut buffers.r_buf,
        &mut buffers.auxiliary,
        diff_ec_ul,
        win_len,
        None,
        parameters.kmer_length as i64,
        &mut buffers.hap.snp_srt,
        &mut buffers.sps,
        emask,
        &mut buffers.mk,
        slots,
    );
    let ton = buffers.overlaps.len();
    buffers.auxiliary = gen_aux_ovlp(&buffers.overlaps, buffers.overlaps.len());
    let mut live = ton;

    let mut phase: u32 = 1;
    gl_chain_flter(
        context,
        &mut buffers.overlaps,
        &mut live,
        &mut Vec::new(),
        &mut buffers.ll_srt,
        &mut buffers.sps,
        &mut buffers.tk,
        &mut buffers.lo,
        diff_ec_ul,
        win_len,
        sequence.len() as i64,
        &mut phase,
    );

    let shared = if phase != 0 {
        gen_shared_intervals(
            &mut buffers.overlaps[..live],
            context,
            win_len,
            &mut buffers.lo,
        )
    } else {
        0
    };
    if phase != 0 && shared != 0 {
        crate::assembly::ul_glchain::filter_topN(
            &mut buffers.overlaps,
            &mut live,
            &mut buffers.lo,
            sequence.len() as u64,
            win_len as u64,
            UL_TOPN,
            &mut buffers.ll_srt,
        );
        // `copy_asg_arr(b0, hap->snp_srt); b1 = sps; b2 = r_buf` — the
        // reference threads these through `update_sketch_trace` and back.
        let mut b0 = std::mem::take(&mut buffers.hap.snp_srt);
        let mut b1 = std::mem::take(&mut buffers.sps);
        let mut b2 = std::mem::take(&mut buffers.r_buf);
        update_sketch_trace(
            &targets,
            &mut buffers.overlaps[..live],
            &mut buffers.lo,
            sequence.len() as i64,
            win_len,
            &mut b0,
            &mut b1,
            &mut b2,
            crate::assembly::ul_phase::max_lgap(sequence.len() as i64),
            diff_ec_ul,
        );
        buffers.hap.snp_srt = b0;
        buffers.sps = b1;
        buffers.r_buf = b2;

        let mut lo = std::mem::take(&mut buffers.lo);
        let live_now = live;
        ul_lalign(
            context,
            &mut buffers.overlaps[..live_now],
            &mut buffers.candidates,
            &targets,
            sequence,
            &mut buffers.extz,
            &mut buffers.dummy,
            &mut buffers.r_buf,
            &mut buffers.auxiliary,
            diff_ec_ul,
            win_len,
            Some(&mut lo),
            parameters.kmer_length as i64,
            &mut buffers.hap.snp_srt,
            &mut buffers.sps,
            emask,
            &mut buffers.mk,
            slots,
        );
        buffers.lo = lo;
        // recover the alignments the second pass dropped: `[ton, olist.length)`
        // keeps their pre-pass coordinates as a single unaligned window
        let grown = buffers.overlaps.len();
        if grown < ton {
            buffers.overlaps.resize(ton, OverlapRegion::default());
        }
        // `for (k = b->olist.length; k < ton; k++)`: repair the dropped
        // entries (still in the tail after the swaps) and restore the count
        let begin = live.min(ton);
        for region in buffers.overlaps[begin..ton].iter_mut() {
            let window = single_window_list(region, region.non_homopolymer_errors);
            region.window_cigar.windows.clear();
            region.window_cigar.windows.push(window);
            region.alignment_length = 0;
            region.overlap_length = region.x_position_end + 1 - region.x_position_start;
        }
        live = ton;
    } else {
        for region in buffers.overlaps.iter_mut() {
            let window = single_window_list(region, 0);
            region.window_cigar.windows.clear();
            region.window_cigar.windows.push(window);
            region.alignment_length = region.x_position_end + 1 - region.x_position_start;
            region.overlap_length = region.alignment_length;
            region.non_homopolymer_errors = 0;
        }
    }

    gl_chain(
        context,
        &mut store.reads[read_id],
        &mut buffers.overlaps[..live],
        &mut buffers.chain_data,
        &mut buffers.hap,
        &mut buffers.sps,
        &mut buffers.tk,
        &mut buffers.lo,
        &mut buffers.ll_srt,
        &mut buffers.gdp,
        diff_ec_ul,
        win_len,
        sequence.len() as i64,
    );
    buffers.num_read_base += 1;
}

/// `rescall_ul_pipeline` (`inter.cpp` 20018): registers, aligns and pushes
/// every read, then rebuilds the read index.
#[allow(clippy::too_many_arguments)]
pub fn rescall_ul_pipeline<'a>(
    options: &AssemblyOptions,
    map_context: &ScallContext<'_>,
    context: &UlRescallContext<'_>,
    parameters: &RescallParameters,
    store: &mut AllUl,
    mut ucr: Option<&mut std::fs::File>,
    emask: &EmaskTable,
) {
    let start_time = logging::realtime_seconds();
    let mut buffers = RescallBuffers::default();
    let mut slots: SecondaryWindows = Vec::new();

    let files: &[String] = options
        .ultra_long_read_files
        .as_ref()
        .map(|files| files.file_names.as_slice())
        .unwrap_or(&[]);
    let mut chunk: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    let mut chunk_bases = 0u64;
    for file in files {
        let mut reader = match FastxReader::open(file) {
            Ok(reader) => reader,
            Err(_) => return,
        };
        while let Some(record) = reader
            .next_record()
            .expect("the UL FASTQ reader must not fail")
        {
            if (record.sequence.len() as u64) < u64::from(parameters.kmer_length)
                || (record.sequence.len() as u64)
                    < u64::from(options.ultra_long_minimum_base_count)
            {
                continue;
            }
            chunk.push((record.name.as_bytes().to_vec(), record.sequence.clone()));
            chunk_bases += record.sequence.len() as u64;
            if chunk_bases >= CHUNK_SIZE {
                flush_rescall_chunk(
                    &chunk,
                    store,
                    &mut buffers,
                    map_context,
                    context,
                    parameters,
                    options,
                    ucr.as_deref_mut(),
                    emask,
                    &mut slots,
                );
                chunk.clear();
                chunk_bases = 0;
            }
        }
    }
    if !chunk.is_empty() {
        flush_rescall_chunk(
            &chunk,
            store,
            &mut buffers,
            map_context,
            context,
            parameters,
            options,
            ucr.as_deref_mut(),
            emask,
            &mut slots,
        );
    }

    eprintln!(
        "[M::rescall_ul_pipeline::{:.3}] ==> Qualification",
        logging::realtime_seconds() - start_time
    );
    eprintln!(
        "[M::rescall_ul_pipeline::] ==> # reads: {}, # bases: 0, # fully corrected reads: {}",
        store.reads.len(),
        buffers.num_correct_base
    );
    gen_ul_vec_rid_t(store, context.reads.total_reads as usize);
}

/// One chunk of the rescall pipeline: register, align, push, and dump the
/// `dd == 3` sequences to the ucr file.
#[allow(clippy::too_many_arguments)]
fn flush_rescall_chunk(
    chunk: &[(Vec<u8>, Vec<u8>)],
    store: &mut AllUl,
    buffers: &mut RescallBuffers,
    map_context: &ScallContext<'_>,
    context: &UlRescallContext<'_>,
    parameters: &RescallParameters,
    options: &AssemblyOptions,
    mut ucr: Option<&mut std::fs::File>,
    emask: &EmaskTable,
    slots: &mut SecondaryWindows,
) {
    for (name, _) in chunk {
        let append = AppendParameters {
            sources: context.sources,
            read_lengths: &context.reads.read_lengths,
            maximum_hang_rate: parameters.maximum_hang_rate,
            maximum_hang: parameters.maximum_hang,
            minimum_overlap: parameters.minimum_overlap,
        };
        append_ul_t(store, None, Some(name), None, &mut [], P_CHAIN_RATE, false, &append);
    }
    let chunk_len = chunk.len();
    for (local_id, (_, sequence)) in chunk.iter().enumerate() {
        let read_id = store.reads.len() - chunk_len + local_id;
        if read_id >= store.reads.len() {
            continue;
        }
        if slots.len() < store.reads.len() {
            slots.resize(store.reads.len(), Vec::new());
        }
        rescall_one_read(
            map_context,
            context,
            buffers,
            store,
            read_id,
            sequence,
            parameters,
            emask,
            slots,
        );
        if store.reads[read_id].dd == 3 {
            if let Some(file) = ucr.as_deref_mut() {
                write_compress_base_disk(file, read_id as u64, sequence);
            }
        }
    }
}

/// `write_compress_base_disk` (`Process_Read.cpp` 2166): the 2-bit encoded
/// read with its N-site list.
pub fn write_compress_base_disk(file: &mut std::fs::File, rid: u64, sequence: &[u8]) {
    use std::io::Write;
    let mut n_site: Vec<u32> = Vec::new();
    let mut packed = vec![0u8; sequence.len() / 4 + 1];
    crate::assembly::ul_store::encode_base_public(&mut packed, sequence, &mut n_site);
    let _ = file.write_all(&rid.to_le_bytes());
    let _ = file.write_all(&(sequence.len() as u32).to_le_bytes());
    let _ = file.write_all(&(n_site.len() as u64).to_le_bytes());
    for site in &n_site {
        let _ = file.write_all(&site.to_le_bytes());
    }
    let _ = file.write_all(&packed);
}

/// `gen_UL_reovlps` (`inter.cpp` 21537) with `is_emask = 1`.
#[allow(clippy::too_many_arguments)]
pub fn gen_ul_reovlps(
    options: &mut AssemblyOptions,
    map_context: &ScallContext<'_>,
    context: &UlRescallContext<'_>,
    gfa_name: &str,
    cutoff: i64,
    unitigs: &[crate::assembly::asg::MaUtg],
    read_graph: &crate::assembly::asg::Asg,
    store: &mut AllUl,
) {
    let parameters = RescallParameters::new(options, cutoff);
    let index = gen_ul_idx(
        options,
        unitigs,
        read_graph,
        context.sources,
        &context.reads.read_lengths,
    );
    let mut index_options = options.clone();
    index_options.k_mer_length = parameters.kmer_length as i32;
    index_options.minimizer_window_size = parameters.minimizer_window as i32;
    let unitig_index = if options.load_index_from_disk != 0 {
        match uidx_load(gfa_name, None, context.reads) {
            Ok(Some(index)) => Some(index),
            _ => None,
        }
    } else {
        None
    };
    let unitig_index = match unitig_index {
        Some(index) => index,
        None => {
            let built =
                build_unitig_minimizer_index(&index_options, unitigs, context.reads, i64::from(cutoff));
            let _ = uidx_write(&built, gfa_name, None);
            built
        }
    };
    let emask_loaded = if options.load_index_from_disk != 0 {
        load_emask_table(gfa_name, unitigs)
    } else {
        None
    };
    let emask = match emask_loaded {
        Some(table) => table,
        None => {
            let table = graph_ovlp_binning(context, options);
            write_emask_table(&table, gfa_name, unitigs);
            table
        }
    };
    let _ = &index;
    let _ = UG_SKIP_GRAPH_N;
    let ucr_path = format!("{}.uidx.ucr.bin", gfa_name);
    let mut ucr_file = std::fs::File::create(&ucr_path).ok();
    rescall_ul_pipeline(
        options,
        map_context,
        context,
        &parameters,
        store,
        ucr_file.as_mut(),
        &emask,
    );
}

/// `ul_realignment` (`inter.cpp` 22871): the renewal entry.
pub fn ul_realignment(
    options: &mut AssemblyOptions,
    map_context: &ScallContext<'_>,
    context: &UlRescallContext<'_>,
    source_lists: &[crate::correction::ec_driver::PafList],
    coverage_cut: &[crate::assembly::hit::CoverageCut],
    sg: &mut crate::assembly::asg::Asg,
    bin_file: &str,
    double_check_cache: u32,
) -> crate::assembly::unitig::UnitigGraph {
    eprintln!("[M::ul_realignment::] ==> starting UL");
    let gfa_name = format!("{}.{}", options.output_file_name, bin_file);
    options.update_coverage(options.homozygous_coverage);
    let cutoff = i64::from(REA_ALIGN_CUTOFF);
    let mut store = AllUl::new();
    let unitig_graph = gen_polished_ug(options, context.reads, source_lists, coverage_cut, sg);
    clear_all_ul_t(&mut store);
    let cached = load_all_ul_t(&gfa_name).ok().flatten();
    match cached {
        Some(cached_store) => {
            store = cached_store;
            let _ = double_check_cache;
        }
        None => {
            gen_ul_reovlps(
                options,
                map_context,
                context,
                &gfa_name,
                cutoff,
                &unitig_graph.unitigs,
                &unitig_graph.graph,
                &mut store,
            );
            let _ = write_all_ul_t(&store, &gfa_name);
        }
    }
    filter_ul_ug(&mut store, &unitig_graph, context.reads);
    gen_ul_vec_rid_t(&mut store, context.reads.total_reads as usize);
    let mut renewed = unitig_graph.clone();
    update_ug_arch_ul_mul(&mut renewed, &store);
    renewed
}

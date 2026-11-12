//! `ul_load` (`inter.cpp` 21693): the whole first-pass ultra-long
//! integration — build the deduplicated HiFi unitig graph and its minimizer
//! index, align every UL read against it (`scall_ul_pipeline`), cache the
//! result and label the UL-supported HiFi overlaps.
//!
//! The reference streams the UL files through a chunked `kt_pipeline`; the
//! port keeps the chunk structure (`CHUNK_SIZE` bases per chunk) but runs the
//! per-read steps sequentially, which is content-equivalent because reads are
//! independent and each chunk's chains are pushed in read order.

use crate::assembly::ul_index::{
    build_unitig_minimizer_index, dedup_hi_fis, uidx_load, uidx_write, UlIndex,
    UnitigMinimizerIndex,
};
use crate::assembly::ul_scall::{
    process_read, push_blocks, register_read, ScallBuffers, ScallContext, ScallParameters,
};
use crate::assembly::ul_src::{update_ovlp_src, update_ovlp_src_bl, UlSourceOptions};
use crate::assembly::ul_store::{gen_ul_vec_rid_t, load_all_ul_t, write_all_ul_t, AllUl};
use crate::command_line::options::AssemblyOptions;
use crate::overlap::paf::PafHit;
use crate::logging;
use crate::reads::ReadStore;
use crate::sequence::fastx_reader::FastxReader;

/// `CHUNK_SIZE` (`inter.cpp` 59): the bases-per-chunk bound of the pipeline.
const CHUNK_SIZE: u64 = 1_000_000_000;

/// `UL_COV_THRES` (`Overlaps.h` 38): how many UL chains must back an overlap
/// for it to survive the weak-overlap filter.
pub const UL_COV_THRES: u32 = 2;

/// `P_CHAIN_COV` (`inter.h` 14).
pub const P_CHAIN_RATE: f32 = 0.985;

/// `print_ovlp_src_bl_stat` (`inter.cpp` 19527): the two statistics lines the
/// reference prints after the source labelling.
pub fn print_ovlp_src_bl_stat(store: &AllUl, reads: &ReadStore, sources: &[Vec<PafHit>]) {
    let mut covered = 0u64;
    let mut chained = 0u64;
    for read_id in 0..reads.total_reads as usize {
        if store.ridx.idx[read_id + 1] - store.ridx.idx[read_id] == 0 {
            continue;
        }
        covered += 1;
        if sources[read_id].iter().any(|hit| hit.block_length != 0) {
            chained += 1;
        }
    }
    eprintln!(
        "[M::print_ovlp_src_bl_stat::] ==> # HiFi reads:{}, # covered HiFi reads:{}, # chained HiFi reads:{}",
        reads.total_reads, covered, chained
    );

    let mut status = [0u64; 4];
    for read in &store.reads {
        status[read.dd as usize] += 1;
    }
    eprintln!(
        "[M::print_ovlp_src_bl_stat::] ==> # passed UL reads:{}, # fully corrected UL reads:{}, # almost fully corrected UL reads:{}, # UL reads have primary chains:{}",
        status[0] + status[1] + status[2] + status[3],
        status[1],
        status[2],
        status[3]
    );
}

/// One chunk's step-1/2/3: register every read, align it, push the chains.
#[allow(clippy::too_many_arguments)]
fn flush_chunk(
    chunk: &[(Vec<u8>, Vec<u8>)],
    chunk_base_id: usize,
    store: &mut AllUl,
    buffers: &mut ScallBuffers,
    context: &ScallContext<'_>,
    parameters: &ScallParameters,
    sources: &[Vec<PafHit>],
    read_lengths: &[u64],
) {
    for (name, sequence) in chunk {
        register_read(
            store,
            name,
            sequence,
            buffers,
            parameters,
            sources,
            read_lengths,
        );
    }
    for (local_id, (_, sequence)) in chunk.iter().enumerate() {
        process_read(buffers, context, local_id, sequence, parameters);
    }
    push_blocks(
        store,
        buffers,
        parameters,
        sources,
        read_lengths,
        chunk_base_id,
    );
}

/// `scall_ul_pipeline` (`inter.cpp` 19992) minus the `init_all_ul_t`: streams
/// every UL file, registers/aligns/pushes each chunk and rebuilds the read
/// index. Returns `(read_count, base_count, fully_covered_reads)`.
#[allow(clippy::too_many_arguments)]
fn scall_ul_pipeline(
    options: &AssemblyOptions,
    reads: &ReadStore,
    sources: &[Vec<PafHit>],
    parameters: &ScallParameters,
    index: &UlIndex,
    unitig_index: &UnitigMinimizerIndex,
    store: &mut AllUl,
) {
    let start_time = logging::realtime_seconds();
    let read_lengths: Vec<u64> = reads.read_lengths.clone();
    let context = ScallContext {
        unitigs: &index.unitig.unitigs,
        contain: &index.contain,
        unitig_index,
        reads,
        sources,
        sample_distance: options.minimizer_sample_distance,
        rewindow_size: options.minimizer_rewindow_size,
    };
    let mut buffers = ScallBuffers::default();
    let mut total_pairs = 0u64;

    let files: &[String] = options
        .ultra_long_read_files
        .as_ref()
        .map(|files| files.file_names.as_slice())
        .unwrap_or(&[]);
    for file in files {
        // The reference pipelines one chunk at a time; with the identical
        // per-read work the chunking only decides when `push_uc_block_t`
        // runs, which does not change the accumulated store.
        let mut chunk: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
        let mut chunk_bases = 0u64;
        let mut reader = match FastxReader::open(file) {
            Ok(reader) => reader,
            Err(_) => return,
        };

        while let Some(record) = reader
            .next_record()
            .expect("the UL FASTQ reader must not fail")
        {
            // `if ((p->ks->seq.l < (uint64_t)p->opt->k) || (p->ks->seq.l <
            // asm_opt.ul_min_base)) continue;`: short reads never register.
            if (record.sequence.len() as u64) < u64::from(parameters.kmer_length)
                || (record.sequence.len() as u64)
                    < u64::from(options.ultra_long_minimum_base_count)
            {
                continue;
            }
            chunk.push((record.name.as_bytes().to_vec(), record.sequence.clone()));
            total_pairs += 1;
            chunk_bases += record.sequence.len() as u64;
            if chunk_bases >= CHUNK_SIZE {
                flush_chunk(
                    &chunk,
                    total_pairs as usize - chunk.len(),
                    store,
                    &mut buffers,
                    &context,
                    parameters,
                    sources,
                    &read_lengths,
                );
                chunk.clear();
                chunk_bases = 0;
            }
        }
        if !chunk.is_empty() {
            flush_chunk(
                &chunk,
                total_pairs as usize - chunk.len(),
                store,
                &mut buffers,
                &context,
                parameters,
                sources,
                &read_lengths,
            );
        }
    }

    eprintln!(
        "[M::scall_ul_pipeline::{:.3}] ==> Qualification",
        logging::realtime_seconds() - start_time
    );
    // `sl->total_base` is only accumulated by the recorrect pipeline, so the
    // scall print always says zero, exactly like the reference.
    eprintln!(
        "[M::scall_ul_pipeline::] ==> # reads: {}, # bases: 0",
        store.reads.len()
    );
    eprintln!(
        "[M::scall_ul_pipeline::] ==> # fully covered reads: {}",
        buffers.fully_covered_reads
    );
    gen_ul_vec_rid_t(store, reads.total_reads as usize);
}

/// `gen_UL_ovlps` (`inter.cpp` 21473): the dedup unitig graph, the cached
/// unitig minimizer index and the scall pass.
fn gen_ul_ovlps(
    options: &AssemblyOptions,
    reads: &ReadStore,
    sources: &[Vec<PafHit>],
    parameters: &ScallParameters,
    cutoff: i64,
    store: &mut AllUl,
) {
    let index = dedup_hi_fis(options, reads, sources, &reads.read_lengths);
    let mut index_options = options.clone();
    index_options.k_mer_length = parameters.kmer_length as i32;
    index_options.minimizer_window_size = parameters.minimizer_window as i32;
    let unitig_index = if options.load_index_from_disk != 0 {
        match uidx_load(&options.output_file_name, None, reads) {
            Ok(Some(index)) => Some(index),
            Ok(None) => None,
            Err(_) => None,
        }
    } else {
        None
    };
    let unitig_index = match unitig_index {
        Some(index) => index,
        None => {
            let built = build_unitig_minimizer_index(
                &index_options,
                &index.unitig.unitigs,
                reads,
                cutoff,
            );
            let _ = uidx_write(&built, &options.output_file_name, None);
            built
        }
    };
    scall_ul_pipeline(
        options,
        reads,
        sources,
        parameters,
        &index,
        &unitig_index,
        store,
    );
}

/// `ul_load` (`inter.cpp` 21693): `create_ul_info`'s driver. Takes and
/// returns the sources so the UL labels (`bl`) survive on the caller's lists.
pub fn ul_load(
    options: &mut AssemblyOptions,
    reads: &ReadStore,
    sources: &mut [Vec<PafHit>],
) -> AllUl {
    eprintln!("[M::ul_load::] ==> UL");
    // `ha_opt_update_cov(&asm_opt, asm_opt.hom_cov)`, then the shared
    // `init_mg_opt(..., 19, 10, cutoff, ...)` of every UL driver.
    options.update_coverage(options.homozygous_coverage);
    let cutoff = options.maximum_chain_count;
    let parameters = ScallParameters {
        kmer_length: 19,
        minimizer_window: 10,
        error_rate: options.ultra_long_error_rate,
        maximum_chain_count: options.maximum_chain_count as u32,
        maximum_hang_rate: options.maximum_hang_rate,
        maximum_hang: i64::from(options.maximum_hang_length),
        minimum_overlap: i64::from(options.minimum_overlap_length),
        p_chain_rate: P_CHAIN_RATE,
        homopolymer_compression: options.flags & crate::command_line::options::FLAG_NO_HPC == 0,
    };

    let mut store = AllUl::new();
    // `if(!load_all_ul_t(&UL_INF, asm_opt.output_file_name, &R_INF, NULL))`:
    // the cache is attempted regardless of `-i`, exactly like the reference.
    match load_all_ul_t(&options.output_file_name) {
        Ok(Some(cached)) => store = cached,
        _ => {
            gen_ul_ovlps(
                options,
                reads,
                sources,
                &parameters,
                i64::from(cutoff),
                &mut store,
            );
            let _ = write_all_ul_t(&store, &options.output_file_name);
        }
    }

    let source_options = UlSourceOptions {
        maximum_hang: options.maximum_hang_length,
        maximum_hang_rate: options.maximum_hang_rate,
        minimum_overlap: options.minimum_overlap_length,
    };
    update_ovlp_src(&store, sources, reads, &source_options);
    update_ovlp_src_bl(sources, reads);
    print_ovlp_src_bl_stat(&store, reads, sources);
    store
}

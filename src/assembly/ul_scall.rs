//! The UL scall pipeline (`ul_load` → `scall_ul_pipeline`): registers the UL
//! reads, runs the per-read worker sequence (`ul_map_lchain` →
//! `ul_rid_lalign_adv` twice with `gen_aux_ovlp` between them →
//! `gen_rid_raw_chain`) and pushes the raw chains into [`AllUl`] with
//! `push_uc_block_t` (`inter.cpp` 12886).
//!
//! The reference batches reads through a `kt_pipeline`; the port processes
//! them sequentially, which is content-equivalent because every read is
//! independent and the buffers accumulate in read order.

use crate::assembly::ul_index::{UlContain, UnitigMinimizerIndex};
use crate::assembly::ul_raw_chain::{gen_rid_raw_chain, GlChain, UlGraphOptions};
use crate::assembly::ul_store::{append_ul_t, AllUl, AppendParameters};
use crate::assembly::asg::MaUtg;
use crate::correction::hc_align::BitExtz;
use crate::correction::target::TargetSource;
use crate::correction::ul_align::{
    gen_aux_ovlp, ul_map_lchain, ul_rid_lalign_adv, UlAlignParameters, UlLocalAlignment,
};
use crate::overlap::paf::PafHit;
use crate::overlap::region::{CandidatesList, OverlapRegion};
use crate::reads::read_store::ReadStore;

/// The `mg_idxopt_t` and `ug_opt_t` slices the scall pass consumes.
#[derive(Debug, Clone, Copy)]
pub struct ScallParameters {
    /// `opt->k`.
    pub kmer_length: u32,
    /// `opt->w`.
    pub minimizer_window: u32,
    /// `opt->diff_ec_ul`.
    pub error_rate: f64,
    /// `opt->max_n_chain`.
    pub maximum_chain_count: u32,
    /// `asm_opt.max_hang_rate`.
    pub maximum_hang_rate: f32,
    /// `uopt->max_hang`.
    pub maximum_hang: i64,
    /// `uopt->min_ovlp`.
    pub minimum_overlap: i64,
    /// `P_CHAIN_COV` (`inter.h` 9).
    pub p_chain_rate: f32,
    /// `opt->is_HPC`.
    pub homopolymer_compression: bool,
}

impl ScallParameters {
    /// `MIN(THRESHOLD_MAX_SIZE / diff_ec_ul, WINDOW)` of the worker.
    pub fn window_length(&self) -> i64 {
        ((31.0 / self.error_rate).min(375.0)) as i64
    }
}

fn processed_reads_debug() -> usize {
    std::env::var("RSHIFIASM_SCALL_READ")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(0)
}

/// The unitig-side context of the pass: `s->uu` plus `s->uopt`'s overlap
/// lists and the HiFi read store.
pub struct ScallContext<'a> {
    /// `uu->ug->u.a`.
    pub unitigs: &'a [MaUtg],
    /// `uu->ct`.
    pub contain: &'a UlContain,
    /// the unitig minimizer position index, `s->ha_idx`.
    pub unitig_index: &'a UnitigMinimizerIndex,
    /// `R_INF`.
    pub reads: &'a ReadStore,
    /// `uopt->sources`.
    pub sources: &'a [Vec<PafHit>],
    /// `asm_opt.mz_sample_dist`.
    pub sample_distance: i32,
    /// `asm_opt.mz_rewin`.
    pub rewindow_size: i32,
}

/// The per-run buffers: the reference's `ha_ovec_buf_t` plus `glchain_t`.
#[derive(Default)]
pub struct ScallBuffers {
    /// `b->exz`.
    pub extz: BitExtz,
    /// `b->olist`: the current read's overlaps.
    pub overlaps: Vec<OverlapRegion>,
    /// `b->clist`: the current read's chain anchors.
    pub candidates: CandidatesList,
    /// `b->aux_o`: the scratch overlap of the second pass.
    pub auxiliary: OverlapRegion,
    /// The second-pass outputs of the current read; the chain buffers below
    /// accumulate their tails.
    pub local: UlLocalAlignment,
    /// `bl`: the accumulating chain buffers (`tk`/`tc`/`srt`; `lo` is reset
    /// per read by `gen_rid_raw_chain`).
    pub chain: GlChain,
    /// `s->seq` of the chunk, indexed by the local read id.
    pub sequences: Vec<Vec<u8>>,
    /// `bb->num_read_base`: the fully covered reads of this run.
    pub fully_covered_reads: u64,
}

fn append_parameters<'a>(
    parameters: &ScallParameters,
    sources: &'a [Vec<PafHit>],
    read_lengths: &'a [u64],
) -> AppendParameters<'a> {
    AppendParameters {
        sources,
        read_lengths,
        maximum_hang_rate: parameters.maximum_hang_rate,
        maximum_hang: parameters.maximum_hang,
        minimum_overlap: parameters.minimum_overlap,
    }
}

/// `scall_ul_pipeline` step 1 minus the chunking: registers one read by name.
pub fn register_read(
    store: &mut AllUl,
    name: &[u8],
    sequence: &[u8],
    buffers: &mut ScallBuffers,
    parameters: &ScallParameters,
    sources: &[Vec<PafHit>],
    read_lengths: &[u64],
) {
    let append = append_parameters(parameters, sources, read_lengths);
    append_ul_t(
        store,
        None,
        Some(name),
        None,
        &mut [],
        parameters.p_chain_rate,
        false,
        &append,
    );
    buffers.sequences.push(sequence.to_vec());
}

/// `worker_for_ul_scall_alignment` (`inter.cpp` 9220): aligns one UL read
/// against the unitig index and appends its raw chains to the buffers.
///
/// Returns whether some chain of the read covers at least 95% of it.
pub fn process_read(
    buffers: &mut ScallBuffers,
    context: &ScallContext<'_>,
    local_read_id: usize,
    sequence: &[u8],
    parameters: &ScallParameters,
) -> bool {
    let window_length = parameters.window_length();
    let align_parameters = UlAlignParameters {
        kmer_length: parameters.kmer_length,
        minimizer_window: parameters.minimizer_window,
        error_rate: parameters.error_rate,
        band_width_threshold: parameters.error_rate,
        maximum_chain_count: parameters.maximum_chain_count,
        homopolymer_compression: parameters.homopolymer_compression,
        sample_distance: context.sample_distance as u32,
        rewindow_size: context.rewindow_size as u32,
    };

    let targets = TargetSource::Unitigs {
        unitigs: context.unitigs,
        reads: context.reads,
    };

    // `ul_map_lchain(b->abl, (uint32_t)-1, ...)`
    let mapped = ul_map_lchain(
        sequence,
        &context.unitig_index.index,
        &targets,
        &align_parameters,
        Some(&context.unitig_index.filter),
    );
    buffers.overlaps = mapped.overlaps;
    buffers.candidates = mapped.candidates;

    if std::env::var("RSHIFIASM_SCALL_DEBUG").is_ok() {
        for (index, overlap) in buffers.overlaps.iter().take(6).enumerate() {
            eprintln!(
                "[scall-debug] map ov{} x=[{},{}] y_id={} y=[{},{}] strand={} fcigar={} aln_len={}",
                index,
                overlap.x_position_start,
                overlap.x_position_end,
                overlap.y_identifier,
                overlap.y_position_start,
                overlap.y_position_end,
                overlap.y_strand,
                overlap.fake_cigar.len(),
                overlap.alignment_length,
            );
        }
    }
    // first pass: keep the chains whose windows cover the overlap
    let overlap_count = buffers.overlaps.len();
    let kept = ul_rid_lalign_adv(
        &mut buffers.overlaps,
        overlap_count,
        &mut buffers.candidates,
        &mut buffers.extz,
        &targets,
        context.unitigs,
        context.contain,
        &context.reads.read_lengths,
        sequence,
        sequence.len() as i64,
        parameters.error_rate,
        window_length,
        &mut buffers.local,
        parameters.kmer_length as i64,
        None,
    );
    if std::env::var("RSHIFIASM_SCALL_DEBUG").is_ok() {
        eprintln!(
            "[scall-debug] after pass1 overlaps={} first={:?}",
            buffers.overlaps.len(),
            buffers
                .overlaps
                .first()
                .map(|o| (o.x_position_start, o.x_position_end, o.alignment_length, o.is_match)),
        );
    }
    if std::env::var("RSHIFIASM_SCALL_DEBUG").is_ok() {
        eprintln!(
            "[rs-pass1] kept={} of={} first={:?}",
            kept,
            buffers.overlaps.len(),
            buffers
                .overlaps
                .first()
                .map(|o| (o.x_position_start, o.x_position_end, o.y_identifier, o.alignment_length)),
        );
    }
    if std::env::var("RSHIFIASM_SCALL_DEBUG").is_ok()
        && local_read_id == processed_reads_debug()
        && std::env::var("RSHIFIASM_SCALL_P1").is_ok()
    {
        for (index, region) in buffers.overlaps.iter().enumerate().take(kept) {
            eprintln!(
                "[p1-reg] {} x=[{},{}] y_id={} y=[{},{}] strand={} aln={} nherr={} fc={}",
                index,
                region.x_position_start,
                region.x_position_end,
                region.y_identifier,
                region.y_position_start,
                region.y_position_end,
                region.y_strand,
                region.alignment_length,
                region.non_homopolymer_errors,
                region.fake_cigar.len()
            );
            for position in 0..region.fake_cigar.len() {
                eprintln!(
                    "[p1-fc] {} {} {} {}",
                    index,
                    position,
                    region.fake_cigar.gap_position(position),
                    region.fake_cigar.gap_shift(position)
                );
            }
            for (position, window) in region.window_cigar.windows.iter().enumerate() {
                eprintln!(
                    "[p1-win] {} {} x=[{},{}] y=[{},{}] eb={} ee={} err={} thre={} cidx={} clen={}",
                    index, position, window.x_start, window.x_end, window.y_start, window.y_end,
                    window.extra_begin, window.extra_end, window.error, window.error_threshold,
                    window.cigar_start, window.cigar_length
                );
            }
        }
    }
    if buffers.overlaps.is_empty() {
        return false;
    }

    // `aux_o = gen_aux_ovlp(&b->olist); cha_idx = bl->tk.n;`
    buffers.auxiliary = gen_aux_ovlp(&buffers.overlaps, buffers.overlaps.len());
    let chain_index = buffers.chain.tk.len();
    let aligned_start = buffers.local.aligned.len();

    // second pass: split every surviving overlap into per-read sub-overlaps
    let mut aux = std::mem::take(&mut buffers.auxiliary);
    ul_rid_lalign_adv(
        &mut buffers.overlaps,
        kept,
        &mut buffers.candidates,
        &mut buffers.extz,
        &targets,
        context.unitigs,
        context.contain,
        &context.reads.read_lengths,
        sequence,
        sequence.len() as i64,
        parameters.error_rate,
        window_length,
        &mut buffers.local,
        parameters.kmer_length as i64,
        Some(&mut aux),
    );
    buffers.auxiliary = aux;
    buffers
        .chain
        .tk
        .extend(buffers.local.aligned[aligned_start..].iter().copied());
    // the reference resets `trace->n = 0` at the start of every read's
    // second pass, so the trace buffer is per-read
    buffers.chain.tc.clear();
    buffers.chain.tc.extend(buffers.local.trace.iter().copied());

    if std::env::var("RSHIFIASM_SCALL_DEBUG").is_ok()
        && local_read_id == processed_reads_debug()
        && std::env::var("RSHIFIASM_SCALL_REGS").is_ok()
    {
        for (index, region) in buffers.overlaps.iter().enumerate() {
            eprintln!(
                "[rs-reg] {} x_id={} x=[{},{}] y_id={} y=[{},{}] strand={} aln={} nherr={} fc={}",
                index,
                region.x_identifier,
                region.x_position_start,
                region.x_position_end,
                region.y_identifier,
                region.y_position_start,
                region.y_position_end,
                region.y_strand,
                region.alignment_length,
                region.non_homopolymer_errors,
                region.fake_cigar.len()
            );
            for position in 0..region.fake_cigar.len() {
                eprintln!(
                    "[rs-fc] {} {} {} {}",
                    index,
                    position,
                    region.fake_cigar.gap_position(position),
                    region.fake_cigar.gap_shift(position)
                );
            }
            for (cell_index, cell) in region.window_cigar.cells.iter().enumerate() {
                eprintln!("[rs-cell] {} {} {}", index, cell_index, cell);
            }
            for (position, window) in region.window_cigar.windows.iter().enumerate() {
                eprintln!(
                    "[rs-win] {} {} x=[{},{}] y=[{},{}] eb={} ee={} err={} thre={} cidx={} clen={}",
                    index, position, window.x_start, window.x_end, window.y_start, window.y_end,
                    window.extra_begin, window.extra_end, window.error, window.error_threshold,
                    window.cigar_start, window.cigar_length
                );
            }
        }
    }
    if std::env::var("RSHIFIASM_SCALL_DEBUG").is_ok()
        && local_read_id == processed_reads_debug()
        && std::env::var("RSHIFIASM_SCALL_USEQ").is_ok()
    {
        for which in 3..6 {
            let region = &buffers.overlaps[which];
            let target_length = targets.length(region.y_identifier as usize);
            let mut span_s = i64::from(region.y_position_start) - 300;
            if span_s < 0 {
                span_s = 0;
            }
            let mut span_l =
                i64::from(region.y_position_end) + 1 - i64::from(region.y_position_start) + 600;
            if span_s + span_l > target_length {
                span_l = target_length - span_s;
            }
            let sequence = targets.sub_region(
                region.y_identifier as usize,
                span_s,
                span_l,
                region.y_strand as u8,
            );
            eprintln!(
                "[rs-useq] which={} y_id={} strand={} span=[{},{}] seq={}",
                which,
                region.y_identifier,
                region.y_strand,
                span_s,
                span_s + span_l,
                String::from_utf8_lossy(&sequence)
            );
        }
    }
    if std::env::var("RSHIFIASM_SCALL_DEBUG").is_ok() && local_read_id == processed_reads_debug() {
        for (index, entry) in buffers.chain.tk.iter().enumerate() {
            eprintln!(
                "[rs-tk] {} qn={} tn={} qs={} qe={} ts={} te={} sec={} el={} rev={}",
                index, entry.qn, entry.tn, entry.qs, entry.qe, entry.ts, entry.te,
                entry.sec, entry.el as u32, entry.rev as u32
            );
        }
    }
    // `gen_rid_raw_chain(&b->olist, bl, cha_idx, &(b->clist.chainDP), ...)`
    let graph_options = UlGraphOptions {
        sources: context.sources,
        minimum_overlap: parameters.minimum_overlap as i32,
        maximum_hang: parameters.maximum_hang as i32,
        maximum_hang_rate: parameters.maximum_hang_rate,
    };
    let fully_covered = gen_rid_raw_chain(
        &mut buffers.chain,
        chain_index,
        &mut buffers.candidates.chain_data,
        &graph_options,
        context.reads,
        context.unitigs,
        parameters.error_rate,
        sequence.len() as i64,
        &targets,
        sequence,
        &mut buffers.extz,
        local_read_id as i64,
        parameters.maximum_chain_count as i64,
        &buffers.overlaps,
    );
    buffers.fully_covered_reads += fully_covered as u64;
    if std::env::var("RSHIFIASM_SCALL_DEBUG").is_ok() {
        eprintln!(
            "[rs-rawchain] tk={} lo={} full={}",
            buffers.chain.tk.len(),
            buffers.chain.lo.len(),
            fully_covered
        );
        for (index, head) in buffers.chain.lo.iter().enumerate() {
            eprintln!(
                "[rs-head] {} qn={} tn={} qs={} qe={} ts={} te={}",
                index, head.qn, head.tn, head.qs, head.qe, head.ts, head.te
            );
        }
    }
    fully_covered != 0
}

/// `push_uc_block_t` (`inter.cpp` 12886) plus the step-3 fallback: pushes the
/// accumulated chains into `store` and registers reads without chains.
///
/// `read_ids` maps the local read index to its global UL id (`s->id + i`).
pub fn push_blocks(
    store: &mut AllUl,
    buffers: &mut ScallBuffers,
    parameters: &ScallParameters,
    sources: &[Vec<PafHit>],
    read_lengths: &[u64],
    read_id_base: usize,
) {
    let append = append_parameters(parameters, sources, read_lengths);
    let mut chains = std::mem::take(&mut buffers.chain.tk);
    let mut start = 0usize;
    for k in 1..=chains.len() {
        if k == chains.len() || chains[k].qn != chains[start].qn {
            let local_id = chains[start].qn as usize;
            let sequence = buffers
                .sequences
                .get(local_id)
                .map(|sequence| sequence.as_slice())
                .unwrap_or(&[]);
            append_ul_t(
                store,
                Some(read_id_base + local_id),
                None,
                Some(sequence),
                &mut chains[start..k],
                parameters.p_chain_rate,
                false,
                &append,
            );
            start = k;
        }
    }
    buffers.chain.tk = chains;

    // the reads without chains: `if((UL_INF.n <= rid) || (UL_INF.a[rid].rlen
    // != s->len[i])) append_ul_t(...)`
    for (local_id, sequence) in buffers.sequences.iter().enumerate() {
        let registered = store
            .reads
            .get(read_id_base + local_id)
            .map(|read| read.rlen == sequence.len() as u32)
            .unwrap_or(false);
        if !registered {
            append_ul_t(
                store,
                Some(read_id_base + local_id),
                None,
                Some(sequence),
                &mut [],
                parameters.p_chain_rate,
                false,
                &append,
            );
        }
    }
}

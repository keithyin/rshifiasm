//! Focused read-0 probe: dumps the second-pass `tk` entries in the same format
//! as `/tmp/ul_worker_probe.out` so they can be diffed against the reference.
//! Ignored by default; requires `/tmp/ulref` and `/tmp/ul.fq`.

use rshifiasm::assembly::hit::{
    normalize_ma_hit_t_single_side_advance, renew_graph_init, try_rescue_overlaps, ReadToUnitig,
};
use rshifiasm::assembly::read_cache;
use rshifiasm::assembly::ul_index::{build_unitig_minimizer_index, dedup_hi_fis};
use rshifiasm::assembly::ul_scall::{
    process_read, register_read, ScallBuffers, ScallContext, ScallParameters,
};
use rshifiasm::command_line::options::AssemblyOptions;
use rshifiasm::overlap::paf::PafHit;
use rshifiasm::sequence::fastx_reader::FastxReader;

#[test]
#[ignore]
fn dump_read0_tk() {
    let which: usize = std::env::var("PROBE_READ")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0);
    let mut options = AssemblyOptions::default();
    options.output_file_name = "/tmp/ulref".to_string();
    let (reads, state) = read_cache::load_all_data_from_disk(&mut options).unwrap();
    let lengths: Vec<u64> = reads.read_lengths.clone();
    let mut read_to_unitig = ReadToUnitig::new(reads.total_reads as usize);
    let (mut forward, forward_flags): (Vec<Vec<_>>, Vec<_>) = state
        .forward
        .iter()
        .map(|l| (l.hits.clone(), (l.is_fully_corrected, l.is_abnormal)))
        .unzip();
    let (mut reverse, _rf): (Vec<Vec<_>>, Vec<_>) = state
        .reverse
        .iter()
        .map(|l| (l.hits.clone(), (l.is_fully_corrected, l.is_abnormal)))
        .unzip();
    try_rescue_overlaps(&mut forward, &mut reverse, &reads.read_lengths, 4, false);
    renew_graph_init(&mut forward, &mut reverse, &mut read_to_unitig);
    normalize_ma_hit_t_single_side_advance(&mut forward, false, None, 0);
    normalize_ma_hit_t_single_side_advance(&mut reverse, false, None, 0);
    let sources: Vec<Vec<PafHit>> = forward;

    let index = dedup_hi_fis(&options, &reads, &sources, &lengths);
    let mut index_options = options.clone();
    index_options.k_mer_length = 19;
    index_options.minimizer_window_size = 10;
    let unitig_index =
        build_unitig_minimizer_index(&index_options, &index.unitig.unitigs, &reads, 180);

    let parameters = ScallParameters {
        kmer_length: 19,
        minimizer_window: 10,
        error_rate: options.ultra_long_error_rate as f64,
        maximum_chain_count: 180,
        maximum_hang_rate: options.maximum_hang_rate,
        maximum_hang: options.maximum_hang_length as i64,
        minimum_overlap: options.minimum_overlap_length as i64,
        p_chain_rate: 0.985,
        homopolymer_compression: options.flags & rshifiasm::command_line::options::FLAG_NO_HPC == 0,
    };

    let mut reader = FastxReader::open("/tmp/ul.fq").unwrap();
    let mut buffers = ScallBuffers::default();
    let mut store = rshifiasm::assembly::ul_store::AllUl::new();
    let mut sequences = Vec::new();
    while let Some(record) = reader.next_record().unwrap() {
        register_read(&mut store, record.name.as_bytes(), &record.sequence, &mut buffers,
            &parameters, &sources, &lengths);
        sequences.push(record.sequence);
    }
    drop(reader);

    let context = ScallContext {
        unitigs: &index.unitig.unitigs,
        contain: &index.contain,
        unitig_index: &unitig_index,
        reads: &reads,
        sources: &sources,
        sample_distance: options.minimizer_sample_distance,
        rewindow_size: options.minimizer_rewindow_size,
    };
    // process only the target read with fresh buffers so the tk dump is that
    // read's entries alone
    {
        process_read(&mut buffers, &context, which, &sequences[which], &parameters);
        eprintln!("RS_H heads={} bodies={}", buffers.chain.lo.len(), buffers.chain.tk.len());
        for (i, e) in buffers.chain.lo.iter().enumerate() {
            eprintln!("RS_HEAD {} qn={} tn={} qs={} qe={} ts={} te={} sec={} el={} rev={}",
                i, e.qn, e.tn, e.qs, e.qe, e.ts, e.te, e.sec, e.el as u32, e.rev as u32);
        }
        for (i, e) in buffers.chain.tk.iter().enumerate() {
            eprintln!("RS_BODY {} qn={} tn={} qs={} qe={} ts={} te={} sec={} el={} rev={}",
                i, e.qn, e.tn, e.qs, e.qe, e.ts, e.te, e.sec, e.el as u32, e.rev as u32);
        }
    }
}
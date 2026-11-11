//! Golden replay of the UL scall pass against the reference `--ul` run.
//!
//! Requires `/tmp/ulref` (the reference caches and unitig index), `/tmp/ul.fq`
//! (the synthetic UL reads) and `reference/ul.ul.ovlp.bin` (the reference's
//! `UL_INF` snapshot, written right after the scall pass). Ignored by
//! default.

use rshifiasm::assembly::hit::{
    normalize_ma_hit_t_single_side_advance, renew_graph_init, try_rescue_overlaps, ReadToUnitig,
};
use rshifiasm::assembly::read_cache;
use rshifiasm::assembly::ul_index::{build_unitig_minimizer_index, dedup_hi_fis};
use rshifiasm::assembly::ul_scall::{
    process_read, push_blocks, register_read, ScallBuffers, ScallContext, ScallParameters,
};
use rshifiasm::assembly::ul_store::load_all_ul_t;
use rshifiasm::command_line::options::AssemblyOptions;
use rshifiasm::correction::ec_driver::PafList;
use rshifiasm::sequence::fastx_reader::FastxReader;

fn source_lists(
    forward: Vec<Vec<rshifiasm::overlap::paf::PafHit>>,
    forward_flags: Vec<(bool, u8)>,
    reverse: Vec<Vec<rshifiasm::overlap::paf::PafHit>>,
    reverse_flags: Vec<(bool, u8)>,
) -> Vec<PafList> {
    let _ = reverse;
    let _ = reverse_flags;
    forward
        .into_iter()
        .zip(forward_flags)
        .map(|(hits, (is_fully_corrected, is_abnormal))| PafList {
            hits,
            is_fully_corrected,
            is_abnormal,
        })
        .collect()
}

#[test]
#[ignore]
fn scall_pass_reproduces_the_reference_ul_overlaps() {
    let mut options = AssemblyOptions::default();
    options.output_file_name = "/tmp/ulref".to_string();
    let (reads, state) = read_cache::load_all_data_from_disk(&mut options)
        .expect("the reference caches must load");
    let lengths: Vec<u64> = reads.read_lengths.clone();
    let mut read_to_unitig = ReadToUnitig::new(reads.total_reads as usize);
    let (mut forward, forward_flags): (Vec<Vec<_>>, Vec<_>) = state
        .forward
        .iter()
        .map(|list| (list.hits.clone(), (list.is_fully_corrected, list.is_abnormal)))
        .unzip();
    let (mut reverse, reverse_flags): (Vec<Vec<_>>, Vec<_>) = state
        .reverse
        .iter()
        .map(|list| (list.hits.clone(), (list.is_fully_corrected, list.is_abnormal)))
        .unzip();
    try_rescue_overlaps(&mut forward, &mut reverse, &reads.read_lengths, 4, false);
    renew_graph_init(&mut forward, &mut reverse, &mut read_to_unitig);
    normalize_ma_hit_t_single_side_advance(&mut forward, false, None, 0);
    normalize_ma_hit_t_single_side_advance(&mut reverse, false, None, 0);
    let sources = source_lists(forward, forward_flags, reverse, reverse_flags);

    let index = dedup_hi_fis(&options, &reads, &sources, &lengths);
    assert_eq!(index.unitig.unitigs.len(), 575);
    assert_eq!(index.unitig.graph.arcs.len(), 1910);

    let mut index_options = options.clone();
    index_options.k_mer_length = 19;
    index_options.minimizer_window_size = 10;
    let unitig_index =
        build_unitig_minimizer_index(&index_options, &index.unitig.unitigs, &reads, 180);
    assert_eq!(unitig_index.filter.len(), 128);

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

    let mut reader =
        FastxReader::open("/tmp/ul.fq").expect("the synthetic UL reads must be present");
    let mut buffers = ScallBuffers::default();
    let mut store = rshifiasm::assembly::ul_store::AllUl::new();
    let mut registered = 0usize;
    while let Some(record) = reader
        .next_record()
        .expect("the FASTQ reader must not fail")
    {
        register_read(
            &mut store,
            record.name.as_bytes(),
            &record.sequence,
            &mut buffers,
            &parameters,
            &sources,
            &lengths,
        );
        registered += 1;
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
    let mut reader =
        FastxReader::open("/tmp/ul.fq").expect("the synthetic UL reads must be present");
    let mut processed = 0usize;
    while let Some(record) = reader
        .next_record()
        .expect("the FASTQ reader must not fail")
    {
        process_read(
            &mut buffers,
            &context,
            processed,
            &record.sequence,
            &parameters,
        );
        processed += 1;
    }
    drop(reader);
    assert_eq!(registered, processed);

    push_blocks(
        &mut store,
        &mut buffers,
        &parameters,
        &sources,
        &lengths,
        0,
    );
    // `scall_ul_pipeline` builds the read index right after `push_uc_block_t`
    rshifiasm::assembly::ul_store::gen_ul_vec_rid_t(&mut store, reads.total_reads as usize);
    eprintln!(
        "reads={} fully_covered={} chain_entries={}",
        store.reads.len(),
        buffers.fully_covered_reads,
        buffers.chain.tk.len()
    );

    let reference = load_all_ul_t("reference/ul")
        .expect("the reference UL cache must load")
        .expect("the reference UL cache must exist");
    assert_eq!(reference.nid.len(), store.nid.len(), "UL read count");
    for (index, (expected, actual)) in reference.nid.iter().zip(store.nid.iter()).enumerate() {
        assert_eq!(expected, actual, "UL read {} name", index);
    }
    assert_eq!(
        reference.reads.len(),
        store.reads.len(),
        "UL read record count"
    );
    let mut matched_reads = 0usize;
    let mut matched_blocks = 0usize;
    let mut reference_blocks = 0usize;
    let mut actual_blocks = 0usize;
    let mut mismatched: Vec<usize> = Vec::new();
    for index in 0..reference.reads.len() {
        let expected = &reference.reads[index];
        let actual = &store.reads[index];
        let name = String::from_utf8_lossy(&store.nid[index]).to_string();
        reference_blocks += expected.bb.len();
        actual_blocks += actual.bb.len();
        let same = expected.dd == actual.dd
            && expected.bb == actual.bb
            && expected.r_base == actual.r_base
            && expected.rlen == actual.rlen;
        if same {
            matched_reads += 1;
        } else {
            mismatched.push(index);
        }
        for (expected_block, actual_block) in expected.bb.iter().zip(actual.bb.iter()) {
            if expected_block == actual_block {
                matched_blocks += 1;
            }
        }
        if index > 400 {
            continue;
        }
        let _ = name;
        assert_eq!(
            expected.rlen, actual.rlen,
            "UL read {} ({}) rlen",
            index, name
        );
        if std::env::var("RSHIFIASM_UL_BLOCKS").is_ok() && index < 2 {
            for block_index in 0..expected.bb.len().max(actual.bb.len()) {
                let reference_block = expected.bb.get(block_index);
                let actual_block = actual.bb.get(block_index);
                eprintln!(
                    "[blocks] read {} i={} ref={:?} act={:?}",
                    index,
                    block_index,
                    reference_block.map(|b| (b.pchain, b.base, b.el, b.hid, b.qs, b.qe, b.ts, b.te, b.rev)),
                    actual_block.map(|b| (b.pchain, b.base, b.el, b.hid, b.qs, b.qe, b.ts, b.te, b.rev)),
                );
            }
        }
        let _ = (expected.dd, actual.dd);
        let _ = name.as_bytes().len();
    }
    eprintln!(
        "[ul-scall-golden] matched_reads={}/{} matched_blocks={}/{} ref_blocks={} act_blocks={} first_mismatches={:?}",
        matched_reads,
        reference.reads.len(),
        matched_blocks,
        reference_blocks,
        reference_blocks,
        actual_blocks,
        &mismatched[..mismatched.len().min(20)],
    );
    assert_eq!(matched_reads, reference.reads.len(), "all UL reads must match");
    assert_eq!(reference.ridx.idx, store.ridx.idx, "UL read index CSR");
    assert_eq!(reference.ridx.occ, store.ridx.occ, "UL read index occ");
}

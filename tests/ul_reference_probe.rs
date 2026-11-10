//! Development probes against the reference `--ul` run in `/tmp/ulref`.
//! These tests are ignored because they need the large reference caches.

use rshifiasm::assembly::read_cache;
use rshifiasm::assembly::hit::{
    normalize_ma_hit_t_single_side_advance, renew_graph_init, try_rescue_overlaps, ReadToUnitig,
};
use rshifiasm::assembly::ul_index::{build_unitig_minimizer_index, dedup_hi_fis};
use rshifiasm::command_line::options::AssemblyOptions;
use rshifiasm::correction::ec_driver::PafList;

#[test]
#[ignore]
fn dedup_index_matches_the_reference_counts() {
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
    let (mut reverse, _): (Vec<Vec<_>>, Vec<_>) = state
        .reverse
        .iter()
        .map(|list| (list.hits.clone(), (list.is_fully_corrected, list.is_abnormal)))
        .unzip();
    try_rescue_overlaps(&mut forward, &mut reverse, &reads.read_lengths, 4, false);
    renew_graph_init(&mut forward, &mut reverse, &mut read_to_unitig);
    normalize_ma_hit_t_single_side_advance(&mut forward, false, None, 0);
    normalize_ma_hit_t_single_side_advance(&mut reverse, false, None, 0);
    let sources: Vec<PafList> = forward
        .into_iter()
        .zip(forward_flags)
        .map(|(hits, (is_fully_corrected, is_abnormal))| PafList {
            hits,
            is_fully_corrected,
            is_abnormal,
        })
        .collect();
    let index = dedup_hi_fis(&options, &reads, &sources, &lengths);
    let total_bases: u64 = index
        .unitig
        .unitigs
        .iter()
        .map(|unitig| unitig.length as u64)
        .sum();
    eprintln!(
        "unitigs={} edges={} bases={}",
        index.unitig.unitigs.len(),
        index.unitig.graph.arcs.len(),
        total_bases
    );
    assert_eq!(index.unitig.unitigs.len(), 575);
    assert_eq!(index.unitig.graph.arcs.len(), 1910);
    assert_eq!(total_bases, 13_197_295);

    let mut index_options = options.clone();
    index_options.k_mer_length = 19;
    index_options.minimizer_window_size = 10;
    let unitig_index =
        build_unitig_minimizer_index(&index_options, &index.unitig.unitigs, &reads, 180);
    assert_eq!(unitig_index.filter.len(), 128);
    assert_eq!(unitig_index.index.total_positions(), 1_716_900);
}

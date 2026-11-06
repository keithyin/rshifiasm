//! Counting of minimizers and construction of the position index.
//!
//! This is the `ha_pt_gen` pipeline of `htab.cpp`: a first pass counts the
//! selected minimizers, a histogram of those counts gives the coverage, and a
//! second pass stores every occurrence in the position index. Unlike the
//! high-occurrence pass it runs without a bloom filter (`HAF_COUNT_EXACT`).

use crate::command_line::options::{AssemblyOptions, FLAG_FAST, FLAG_NO_HPC};
use crate::counting::filter_table::KmerFilterTable;
use crate::counting::kmer_counter::{KmerCountTable, MAXIMUM_COUNT};
use crate::counting::minimizer::sketch;
use crate::counting::peak_analysis::analyze_count;
use crate::counting::position_index::{PositionIndex, PositionRecord};
use crate::logging;
use crate::reads::ReadStore;

/// Result of the minimizer counting pipeline.
pub struct MinimizerCounts {
    /// Counts of every selected minimizer, after singleton removal.
    pub counts: KmerCountTable,
    /// Position index, empty in fast mode.
    pub index: Option<PositionIndex>,
    /// Total number of minimizers selected over all reads, `n_mz`.
    pub total_minimizers: u64,
    /// Inferred homozygous coverage (`peak_hom`).
    pub homozygous_coverage: i32,
    /// Inferred heterozygous coverage (`peak_het`).
    pub heterozygous_coverage: i32,
}

/// `ha_pt_gen`: counts minimizers, infers the coverage and builds the index.
///
/// `filter_table` is `None` when the run disabled the high-occurrence filter
/// (`HA_F_NO_KMER_FLT`); the reference then derives the shrink cutoff from the
/// coverage instead of using `YAK_MAX_COUNT - 1`.
pub fn count_minimizers(
    options: &AssemblyOptions,
    reads: &ReadStore,
    filter_table: Option<&KmerFilterTable>,
) -> MinimizerCounts {
    let use_homopolymer_compression = options.flags & FLAG_NO_HPC == 0;
    let mut counts = KmerCountTable::new(options.k_mer_length, 0, 0);
    let mut total_minimizers = 0u64;

    for read_id in 0..reads.total_reads as usize {
        let sequence = reads.read_sequence(read_id);
        let sketch_result = sketch_read(options, &sequence, read_id, use_homopolymer_compression, filter_table);
        for minimizer in &sketch_result {
            counts.insert(minimizer.hash);
        }
        total_minimizers += sketch_result.len() as u64;
    }
    counts.set_total_bases(reads.total_bases);

    if options.minimizer_window_size > 1 {
        eprintln!("[M::yak_count] collected {} minimizers", total_minimizers);
    }
    eprintln!(
        "[M::ha_pt_gen::{:.3}*{:.2}] ==> counted {} distinct minimizer k-mers",
        logging::realtime_seconds(),
        logging::cpu_usage(),
        counts.distinct_count()
    );

    let histogram = counts.histogram();
    eprintln!(
        "[M::ha_pt_gen] count[{}] = {} (for sanity check)",
        MAXIMUM_COUNT, histogram[MAXIMUM_COUNT as usize]
    );

    let mut heterozygous_coverage = -1i32;
    // `asm_opt->hg_size > 0 ? ct->bs / asm_opt->hg_size : -1`.
    let expected_homozygous_peak = if options.haploid_genome_size > 0 {
        (counts.total_bases() / options.haploid_genome_size as u64) as i64
    } else {
        -1
    };
    let homozygous_coverage = analyze_count(
        &histogram,
        options.minimum_histogram_kmer_count,
        expected_homozygous_peak,
        &mut heterozygous_coverage,
    );
    if homozygous_coverage > 0 {
        eprintln!(
            "[M::ha_pt_gen] peak_hom: {}; peak_het: {}",
            homozygous_coverage, heterozygous_coverage
        );
    }

    // Without the high-occurrence filter the cutoff is the coverage times the
    // high-occurrence factor; with it, only the counts that saturate the
    // counter need to go.
    let cutoff = if filter_table.is_none() {
        let cutoff = (f64::from(homozygous_coverage) * options.high_occurrence_factor) as i64;
        cutoff.min(MAXIMUM_COUNT as i64 - 1)
    } else {
        MAXIMUM_COUNT as i64 - 1
    };
    // Singletons never help to anchor an overlap.
    let remaining = counts.shrink(2, cutoff);
    eprintln!(
        "[M::ha_ct_shrink::{:.3}*{:.2}] ==> counted {} distinct minimizer k-mers",
        logging::realtime_seconds(),
        logging::cpu_usage(),
        remaining
    );
    // `tot_cnt`: occurrences surviving the shrink, which the index must hold.
    let surviving_occurrences: u64 = (2..=cutoff)
        .map(|count| histogram[count as usize] as u64 * count as u64)
        .sum();

    let index = if options.flags & FLAG_FAST == 0 {
        eprintln!("[M::ha_pt_gen::] counting in normal mode");
        let mut index = PositionIndex::build(&counts, options.k_mer_length);
        for read_id in 0..reads.total_reads as usize {
            let sequence = reads.read_sequence(read_id);
            let minimizers =
                sketch_read(options, &sequence, read_id, use_homopolymer_compression, filter_table);
            for minimizer in minimizers {
                index.insert(
                    minimizer.hash,
                    PositionRecord {
                        read_identifier: minimizer.read_identifier,
                        position: minimizer.position,
                        reverse: minimizer.reverse,
                        span: minimizer.span,
                    },
                );
            }
        }
        if options.minimizer_window_size > 1 {
            eprintln!("[M::yak_count] collected {} minimizers", total_minimizers);
        }
        debug_assert_eq!(
            surviving_occurrences,
            index.total_positions(),
            "the second pass must index exactly the occurrences that survived the shrink"
        );
        eprintln!(
            "[M::ha_pt_gen::{:.3}*{:.2}] ==> indexed {} positions, counted {} distinct minimizer k-mers",
            logging::realtime_seconds(),
            logging::cpu_usage(),
            index.total_positions(),
            index.distinct_count()
        );
        Some(index)
    } else {
        eprintln!("[M::ha_pt_gen::] counting in fast mode");
        None
    };

    MinimizerCounts {
        counts,
        index,
        total_minimizers,
        homozygous_coverage,
        heterozygous_coverage,
    }
}

/// Runs the sketch with the same parameters the counting pipeline uses.
fn sketch_read(
    options: &AssemblyOptions,
    sequence: &[u8],
    read_id: usize,
    use_homopolymer_compression: bool,
    filter_table: Option<&KmerFilterTable>,
) -> Vec<crate::counting::minimizer::Minimizer> {
    sketch(
        sequence,
        options.minimizer_window_size,
        options.k_mer_length,
        read_id as u32,
        use_homopolymer_compression,
        filter_table,
        options.minimizer_sample_distance,
        options.minimizer_rewindow_size,
    )
    .minimizers
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temporary_fasta(name: &str, content: &[u8]) -> String {
        let path = std::env::temp_dir().join(name);
        let mut file = std::fs::File::create(&path).unwrap();
        file.write_all(content).unwrap();
        path.to_string_lossy().to_string()
    }

    #[test]
    fn builds_an_index_over_the_selected_minimizers() {
        // Reads must be longer than w + k - 1 to yield any minimizer at the
        // default k = 51, w = 51.
        let sequence = b"ACGTTGCAACGTAGCTAGCTAGGCTAACCGTACGATCGTAGCTAGCTA".repeat(6);
        let record = format!(
            ">r1\n{}\n>r2\n{}\n",
            String::from_utf8_lossy(&sequence),
            String::from_utf8_lossy(&sequence)
        );
        let path = write_temporary_fasta("rshifiasm_position_index.fa", record.as_bytes());

        let mut options = AssemblyOptions::default();
        options.input_file_names = vec![path];
        options.input_file_count = 1;
        options.read_length_cutoff = -1;

        let reads = ReadStore::load_from_files(&options).unwrap();
        let empty_counts = KmerCountTable::new(51, 0, 0);
        let filter_table = KmerFilterTable::from_counts(&empty_counts, 2000);

        let result = count_minimizers(&options, &reads, Some(&filter_table));
        assert!(result.total_minimizers > 0);

        // Every counted k-mer occurs at least twice, so the index must hold
        // exactly as many occurrences as there were minimizers.
        let index = result.index.expect("fast mode is off by default");
        assert_eq!(index.total_positions(), result.total_minimizers);
        assert_eq!(index.capacity(), result.total_minimizers);
    }
}

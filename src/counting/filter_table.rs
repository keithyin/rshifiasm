//! The high-occurrence filter table, the counterpart of `gen_hh` and
//! `ha_ft_cnt` in `htab.cpp`.
//!
//! K-mers that occur too often are repetitive and must not create overlaps.
//! The table stores their counts, saturating the ones above the configured
//! limit so that the lookup simply reports "filtered".

use crate::command_line::options::{AssemblyOptions, FLAG_NO_HPC};
use crate::container::khash::{identity_hash, U64HashTable};
use crate::counting::kmer_counter::{
    for_each_homopolymer_compressed_kmer, for_each_kmer, KmerCountTable, MAXIMUM_COUNT,
};
use crate::counting::peak_analysis::analyze_count;
use crate::logging;
use crate::reads::ReadStore;

/// Sentinel stored for counts above the configured limit (`INT16_MAX`).
pub const SATURATED_COUNT: u64 = i16::MAX as u64;

/// Value returned by [`KmerFilterTable::count`] for saturated entries
/// (`INT32_MAX`), which the sketch treats as "filtered".
pub const FILTERED_COUNT: i32 = i32::MAX;

/// Filter table mapping a k-mer hash to its (possibly saturated) count.
pub struct KmerFilterTable {
    table: U64HashTable,
}

impl KmerFilterTable {
    /// `gen_hh`: builds the table from a counting table.
    ///
    /// `maximum_count` is `--max-kocc`, clamped exactly like the reference.
    pub fn from_counts(counts: &KmerCountTable, maximum_count: i64) -> Self {
        let mut limit = maximum_count;
        if limit > MAXIMUM_COUNT as i64 - 1 {
            limit = MAXIMUM_COUNT as i64 - 1;
        }
        if limit > i16::MAX as i64 - 1 {
            limit = i16::MAX as i64 - 1;
        }

        let mut table = U64HashTable::new(identity_hash);
        table.resize((counts.distinct_count() as u32).saturating_mul(2));

        for (kmer_hash, count) in counts.iter_kmer_counts() {
            let (bucket, absent) = table.put(kmer_hash);
            if absent {
                let stored = if count > limit as u64 {
                    SATURATED_COUNT
                } else {
                    count
                };
                table.set_value_at(bucket, stored);
            }
        }

        Self { table }
    }

    /// Number of stored k-mers, `kh_size`.
    pub fn len(&self) -> usize {
        self.table.len()
    }

    pub fn is_empty(&self) -> bool {
        self.table.is_empty()
    }

    /// `ha_ft_cnt`: the count for `kmer_hash`, zero when absent and
    /// [`FILTERED_COUNT`] when saturated.
    pub fn count(&self, kmer_hash: u64) -> i32 {
        match self.table.get_value(kmer_hash) {
            None => 0,
            Some(value) if value == SATURATED_COUNT => FILTERED_COUNT,
            Some(value) => value as i32,
        }
    }
}

/// `ha_ft_gen`: counts every k-mer of the input reads, infers the homozygous
/// coverage and builds the high-occurrence filter table.
///
/// Returns the table together with the inferred `peak_hom`; a negative value
/// means that the coverage could not be inferred from the histogram.
pub fn generate_high_occurrence_filter(
    options: &AssemblyOptions,
    reads: &ReadStore,
) -> (KmerFilterTable, i32) {
    let use_homopolymer_compression = options.flags & FLAG_NO_HPC == 0;
    let mut counts = KmerCountTable::new(
        options.k_mer_length,
        options.bloom_filter_shift,
        4,
    );

    // `HAF_COUNT_ALL` forces the window to one, i.e. every k-mer is counted
    // rather than only minimizers.
    for read_id in 0..reads.total_reads as usize {
        let sequence = reads.read_sequence(read_id);
        if use_homopolymer_compression {
            for_each_homopolymer_compressed_kmer(
                &sequence,
                options.k_mer_length as u32,
                |kmer_hash| {
                    counts.insert(kmer_hash);
                },
            );
        } else {
            for_each_kmer(&sequence, options.k_mer_length as u32, |kmer_hash| {
                counts.insert(kmer_hash);
            });
        }
    }
    counts.set_total_bases(reads.total_bases);

    let histogram = counts.histogram();
    let expected_homozygous_peak = if options.haploid_genome_size > 0 {
        (counts.total_bases() / options.haploid_genome_size as u64) as i64
    } else {
        -1
    };
    let mut heterozygous_peak = -1i32;
    let homozygous_peak = analyze_count(
        &histogram,
        options.minimum_histogram_kmer_count,
        expected_homozygous_peak,
        &mut heterozygous_peak,
    );
    if homozygous_peak > 0 {
        eprintln!(
            "[M::ha_ft_gen] peak_hom: {}; peak_het: {}",
            homozygous_peak, heterozygous_peak
        );
    }

    let mut cutoff = (homozygous_peak as f64 * options.high_occurrence_factor) as i64;
    if cutoff > MAXIMUM_COUNT as i64 - 1 {
        cutoff = MAXIMUM_COUNT as i64 - 1;
    }

    let remaining = counts.shrink(cutoff, MAXIMUM_COUNT as i64);
    eprintln!(
        "[M::ha_ct_shrink::{:.3}*{:.2}] ==> counted {} distinct minimizer k-mers",
        logging::realtime_seconds(),
        logging::cpu_usage(),
        remaining
    );

    let table = KmerFilterTable::from_counts(&counts, options.maximum_kmer_occurrence);
    eprintln!(
        "[M::ha_ft_gen::{:.3}*{:.2}@{:.3}GB] ==> filtered out {} k-mers occurring {} or more times",
        logging::realtime_seconds(),
        logging::cpu_usage(),
        logging::peak_resident_memory_gigabytes(),
        table.len(),
        cutoff
    );

    (table, homozygous_peak)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stores_counts_and_saturates_above_the_limit() {
        let mut counts = KmerCountTable::new(3, 0, 0);
        for _ in 0..3 {
            counts.insert(111);
        }
        for _ in 0..30 {
            counts.insert(222);
        }

        let table = KmerFilterTable::from_counts(&counts, 10);
        assert_eq!(table.len(), 2);
        assert_eq!(table.count(111), 3);
        assert_eq!(table.count(222), FILTERED_COUNT);
        assert_eq!(table.count(333), 0);
    }

    #[test]
    fn maximum_count_is_clamped_like_the_reference() {
        let mut counts = KmerCountTable::new(3, 0, 0);
        counts.insert(7);
        // A limit above the saturating sentinel is clamped, so no entry can be
        // marked as filtered through it.
        let table = KmerFilterTable::from_counts(&counts, 1_000_000);
        assert_eq!(table.count(7), 1);
    }
}

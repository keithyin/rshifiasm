//! K-mer counting, the counterpart of the `yak_count` machinery in `htab.cpp`.
//!
//! The table is split into `1 << partition_bits` independent partitions, each a
//! `khashl` set whose key packs the occurrence count in the low bits. A blocked
//! bloom filter per partition decides whether a new key may be created, which
//! is what keeps singletons out of the table.

use crate::container::khash::{hash_without_counter_bits, U64HashTable};
use crate::counting::bloom_filter::BlockedBloomFilter;
use crate::counting::hash::hash_long_kmer;
use crate::sequence::encoding::SEQUENCE_NT6_TABLE;

/// `YAK_COUNTER_BITS`.
pub const COUNTER_BITS: u32 = 12;

/// Number of histogram buckets, `YAK_N_COUNTS`.
pub const COUNT_BUCKET_COUNT: usize = 1 << COUNTER_BITS;

/// Maximum stored count, `YAK_MAX_COUNT`.
pub const MAXIMUM_COUNT: u64 = (1 << COUNTER_BITS) - 1;

/// Number of partitions, `1 << pre`.
pub const PARTITION_COUNT: usize = 1 << COUNTER_BITS;

/// Hash used by the counting set: ignores the packed count.
fn count_table_hash(key: u64) -> u32 {
    hash_without_counter_bits(key, COUNTER_BITS)
}

/// A counting set: equality and hashing both ignore the packed counter, which
/// is what makes `key + 1` an in place counter increment.
fn new_count_table() -> U64HashTable {
    U64HashTable::with_ignored_low_bits(count_table_hash, COUNTER_BITS)
}

/// One partition: a counting set plus the optional bloom filter in front of it.
struct CountPartition {
    table: U64HashTable,
    bloom_filter: Option<BlockedBloomFilter>,
}

impl CountPartition {
    fn new(bloom_filter: Option<BlockedBloomFilter>) -> Self {
        Self {
            table: new_count_table(),
            bloom_filter,
        }
    }

    /// `ha_ct_insert_list` for a single k-mer, `create_new` semantics included.
    fn insert(&mut self, kmer_hash: u64) -> bool {
        let condensed = kmer_hash >> COUNTER_BITS;
        let counter_bit = u64::from(self.bloom_filter.is_some());

        let may_create = match &mut self.bloom_filter {
            None => true,
            Some(filter) => filter.insert(condensed) == filter.hash_count(),
        };

        if !may_create {
            return false;
        }

        let key = (condensed << COUNTER_BITS) | counter_bit;
        let (bucket, absent) = self.table.put(key);
        if self.table.key_at(bucket) & MAXIMUM_COUNT < MAXIMUM_COUNT {
            let incremented = self.table.key_at(bucket) + 1;
            self.table.set_key_at(bucket, incremented);
        }
        absent
    }

}

/// The whole counting table, `ha_ct_t`.
pub struct KmerCountTable {
    /// K-mer length the table was built for.
    kmer_length: i32,
    partitions: Vec<CountPartition>,
    /// Total bases seen, `ha_ct_t::bs`.
    total_bases: u64,
}

impl KmerCountTable {
    /// `ha_ct_init`.
    ///
    /// `bloom_shift` of zero (or not larger than the partition bits) disables
    /// the bloom filter, exactly like the reference implementation.
    pub fn new(kmer_length: i32, bloom_shift: i32, bloom_hash_count: i32) -> Self {
        let bloom_enabled = bloom_hash_count > 0 && bloom_shift as u32 > COUNTER_BITS;
        let partitions = (0..PARTITION_COUNT)
            .map(|_| {
                let filter = if bloom_enabled {
                    BlockedBloomFilter::new(
                        bloom_shift as u32 - COUNTER_BITS,
                        bloom_hash_count as u32,
                    )
                } else {
                    None
                };
                CountPartition::new(filter)
            })
            .collect();

        Self {
            kmer_length,
            partitions,
            total_bases: 0,
        }
    }

    pub fn kmer_length(&self) -> i32 {
        self.kmer_length
    }

    /// Adds one k-mer hash, returning `true` when the key was newly created.
    pub fn insert(&mut self, kmer_hash: u64) -> bool {
        let partition = (kmer_hash & (PARTITION_COUNT as u64 - 1)) as usize;
        self.partitions[partition].insert(kmer_hash)
    }

    /// Iterates one partition in bucket order, which is the order
    /// `worker_pt_gen` walks when it lays out the position index.
    pub fn partition_entries(&self, partition_index: usize) -> impl Iterator<Item = (u64, u64)> + '_ {
        self.partitions[partition_index].table.iter()
    }

    /// Number of distinct k-mers currently stored.
    pub fn distinct_count(&self) -> u64 {
        self.partitions
            .iter()
            .map(|partition| partition.table.len() as u64)
            .sum()
    }

    /// `ha_ct_hist`: counts how many k-mers occur exactly `index` times.
    pub fn histogram(&self) -> Vec<i64> {
        let mut counts = vec![0i64; COUNT_BUCKET_COUNT];
        for partition in &self.partitions {
            for (key, _) in partition.table.iter() {
                counts[(key & MAXIMUM_COUNT) as usize] += 1;
            }
        }
        counts
    }

    /// `ha_ct_shrink`: keeps only the k-mers whose count is within
    /// `[minimum, maximum]` and returns the number of survivors.
    pub fn shrink(&mut self, minimum: i64, maximum: i64) -> u64 {
        let mut survivors = 0u64;
        for partition in &mut self.partitions {
            let mut rebuilt = new_count_table();
            for (key, _) in partition.table.iter() {
                let count = (key & MAXIMUM_COUNT) as i64;
                if count >= minimum && count <= maximum {
                    let (bucket, _) = rebuilt.put(key);
                    rebuilt.set_value_at(bucket, 0);
                }
            }
            partition.table = rebuilt;
            survivors += partition.table.len() as u64;
        }
        survivors
    }

    /// Rebuilds every k-mer hash from its partition index and packed count.
    ///
    /// This is what `gen_hh` does while walking the partitions:
    /// `key >> pre << COUNTER_BITS | partition`.
    pub fn iter_kmers(&self) -> impl Iterator<Item = u64> + '_ {
        self.partitions.iter().enumerate().flat_map(|(index, partition)| {
            partition
                .table
                .iter()
                .map(move |(key, _)| (key >> COUNTER_BITS << COUNTER_BITS) | index as u64)
        })
    }

    /// Same as [`Self::iter_kmers`] but yields the occurrence count as well.
    pub fn iter_kmer_counts(&self) -> impl Iterator<Item = (u64, u64)> + '_ {
        self.partitions.iter().enumerate().flat_map(|(index, partition)| {
            partition.table.iter().map(move |(key, _)| {
                (
                    (key >> COUNTER_BITS << COUNTER_BITS) | index as u64,
                    key & MAXIMUM_COUNT,
                )
            })
        })
    }

    pub fn set_total_bases(&mut self, total_bases: u64) {
        self.total_bases = total_bases;
    }

    pub fn total_bases(&self) -> u64 {
        self.total_bases
    }

}

/// `sf_count_seq_buf`: enumerates every k-mer of a sequence without HPC.
pub fn for_each_kmer(sequence: &[u8], kmer_length: u32, mut visit: impl FnMut(u64)) {
    let mask = if kmer_length >= 64 {
        u64::MAX
    } else {
        (1u64 << kmer_length) - 1
    };
    let shift = kmer_length - 1;

    let mut encoded = [0u64; 4];
    let mut valid_length = 0u32;

    for base in sequence {
        let code = SEQUENCE_NT6_TABLE[*base as usize];
        if code < 4 {
            let code = code as u64;
            encoded[0] = (encoded[0] << 1 | (code & 1)) & mask;
            encoded[1] = (encoded[1] << 1 | (code >> 1)) & mask;
            encoded[2] = encoded[2] >> 1 | (1 - (code & 1)) << shift;
            encoded[3] = encoded[3] >> 1 | (1 - (code >> 1)) << shift;
            valid_length += 1;
            if valid_length >= kmer_length {
                visit(hash_long_kmer(&encoded));
            }
        } else {
            valid_length = 0;
            encoded = [0u64; 4];
        }
    }
}

/// `sf_count_seq_buf_HPC`: enumerates k-mers over the homopolymer-compressed
/// sequence, which is what the default HiFi configuration uses.
pub fn for_each_homopolymer_compressed_kmer(
    sequence: &[u8],
    kmer_length: u32,
    mut visit: impl FnMut(u64),
) {
    let mask = if kmer_length >= 64 {
        u64::MAX
    } else {
        (1u64 << kmer_length) - 1
    };
    let shift = kmer_length - 1;

    let mut encoded = [0u64; 4];
    let mut valid_length = 0u32;
    let mut last_base: i64 = -1;

    for base in sequence {
        let code = SEQUENCE_NT6_TABLE[*base as usize];
        if code >= 4 {
            valid_length = 0;
            last_base = -1;
            encoded = [0u64; 4];
            continue;
        }

        if i64::from(code) == last_base {
            continue;
        }

        let code = code as u64;
        encoded[0] = (encoded[0] << 1 | (code & 1)) & mask;
        encoded[1] = (encoded[1] << 1 | (code >> 1)) & mask;
        encoded[2] = encoded[2] >> 1 | (1 - (code & 1)) << shift;
        encoded[3] = encoded[3] >> 1 | (1 - (code >> 1)) << shift;
        valid_length += 1;
        if valid_length >= kmer_length {
            visit(hash_long_kmer(&encoded));
        }
        last_base = code as i64;
    }
}

/// Number of k-mers a sequence contributes, used for memory estimation.
pub fn count_kmers(sequence: &[u8], kmer_length: u32, use_hpc: bool) -> u64 {
    let mut total = 0u64;
    if use_hpc {
        for_each_homopolymer_compressed_kmer(sequence, kmer_length, |_| total += 1);
    } else {
        for_each_kmer(sequence, kmer_length, |_| total += 1);
    }
    total
}

/// Hash of a single k-mer window, exposed for unit tests.
pub fn kmer_hash_of(sequence: &[u8]) -> u64 {
    let mut hash = 0u64;
    for_each_kmer(sequence, sequence.len() as u32, |value| hash = value);
    hash
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn homopolymer_compressed_kmer_counts_collapse_runs() {
        // "AAAA" compresses to a single A, so no 2-mer is produced.
        let mut produced = 0;
        for_each_homopolymer_compressed_kmer(b"AAAA", 2, |_| produced += 1);
        assert_eq!(produced, 0);

        // "ACGT" compresses to ACGT, so two 3-mers are produced.
        let mut produced = 0;
        for_each_homopolymer_compressed_kmer(b"ACGT", 3, |_| produced += 1);
        assert_eq!(produced, 2);
    }

    #[test]
    fn ambiguous_bases_restart_the_window() {
        let mut produced = 0;
        for_each_kmer(b"ACGTNACGT", 3, |_| produced += 1);
        assert_eq!(produced, 4);
    }

    #[test]
    fn counting_saturates_at_the_maximum_count() {
        let mut table = KmerCountTable::new(3, 0, 0);
        // Without a bloom filter every occurrence is counted, up to the cap.
        for _ in 0..(MAXIMUM_COUNT + 10) {
            table.insert(12345);
        }
        let histogram = table.histogram();
        assert_eq!(histogram[MAXIMUM_COUNT as usize], 1);
        assert_eq!(table.distinct_count(), 1);
    }

    #[test]
    fn bloom_filter_delays_the_first_creation() {
        let mut table = KmerCountTable::new(3, 37, 4);
        // The first occurrence only sets bloom bits, so nothing is stored yet.
        table.insert(987_654);
        assert_eq!(table.distinct_count(), 0);

        table.insert(987_654);
        assert_eq!(table.distinct_count(), 1);
    }

    #[test]
    fn shrink_keeps_only_the_requested_count_range() {
        let mut table = KmerCountTable::new(3, 0, 0);
        for occurrence in 0..10 {
            for _ in 0..=occurrence {
                table.insert(1000 + occurrence);
            }
        }

        // Without a bloom filter the first insertion already counts as one,
        // so the k-mers seen 4..9 times survive the [5, 4095] window.
        let survivors = table.shrink(5, MAXIMUM_COUNT as i64);
        assert_eq!(survivors, 6);
        assert_eq!(table.distinct_count(), 6);
    }

    #[test]
    fn reverse_complementary_kmers_hash_identically() {
        // ACGTACG and CGTACGT are reverse complements of each other, so the
        // canonical strand choice must collapse them onto the same hash.
        assert_eq!(kmer_hash_of(b"ACGTACG"), kmer_hash_of(b"CGTACGT"));
        assert_ne!(kmer_hash_of(b"ACGTACG"), kmer_hash_of(b"ACGTACC"));
    }
}

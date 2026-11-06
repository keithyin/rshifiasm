//! Position index, the counterpart of `ha_pt_t` in `htab.cpp`.
//!
//! For every minimizer hash the index stores the list of places where it
//! occurs: read identifier, position, strand and span. The list of a hash is
//! filled in sketch order, which is what later stages rely on when they walk
//! the occurrences of a k-mer.

use crate::container::khash::{hash_without_counter_bits, U64HashTable};
use crate::counting::kmer_counter::{KmerCountTable, COUNTER_BITS, MAXIMUM_COUNT, PARTITION_COUNT};

/// One occurrence of a minimizer, `ha_idxpos_t`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PositionRecord {
    /// Read the minimizer was found in.
    pub read_identifier: u32,
    /// Position of the last base of the k-mer.
    pub position: u32,
    /// Strand of the k-mer.
    pub reverse: bool,
    /// Number of bases spanned by the k-mer.
    pub span: u32,
}

/// Hash used inside a partition: the low counter bits are ignored.
fn index_hash(key: u64) -> u32 {
    hash_without_counter_bits(key, COUNTER_BITS)
}

fn new_index_table() -> U64HashTable {
    U64HashTable::with_ignored_low_bits(index_hash, COUNTER_BITS)
}

/// One partition of the index: a hash from k-mer to its occurrence list plus
/// the occurrence storage itself.
struct IndexPartition {
    table: U64HashTable,
    occurrences: Vec<PositionRecord>,
}

/// The whole position index.
pub struct PositionIndex {
    kmer_length: i32,
    partitions: Vec<IndexPartition>,
    total_positions: u64,
}

impl PositionIndex {
    /// `ha_pt_gen`: builds an empty occurrence list for every counted k-mer.
    ///
    /// The start offset of each list is the running sum of the counts of all
    /// preceding k-mers in the same partition, exactly like `worker_pt_gen`.
    pub fn build(counts: &KmerCountTable, kmer_length: i32) -> Self {
        let mut partitions: Vec<IndexPartition> = (0..PARTITION_COUNT)
            .map(|_| IndexPartition {
                table: new_index_table(),
                occurrences: Vec::new(),
            })
            .collect();

        for (partition_index, partition) in partitions.iter_mut().enumerate() {
            let mut next_start = 0u64;
            for (key, _) in counts.partition_entries(partition_index) {
                let count = key & MAXIMUM_COUNT;
                let (bucket, _) = partition.table.put(key & !MAXIMUM_COUNT);
                partition.table.set_value_at(bucket, next_start);
                next_start += count;
            }
            partition.occurrences = vec![PositionRecord::default(); next_start as usize];
        }

        Self {
            kmer_length,
            partitions,
            total_positions: 0,
        }
    }

    pub fn kmer_length(&self) -> i32 {
        self.kmer_length
    }

    /// Number of stored occurrences, `ha_pt_t::tot_pos`.
    pub fn total_positions(&self) -> u64 {
        self.total_positions
    }

    /// Total number of occurrences over every hash, which equals the number of
    /// minimizers that were selected when the counts were collected.
    pub fn capacity(&self) -> u64 {
        self.partitions
            .iter()
            .map(|partition| partition.occurrences.len() as u64)
            .sum()
    }

    /// `ha_pt_insert_list` for a single minimizer.
    ///
    /// Returns `false` when the hash is not part of the index, which happens
    /// when the second pass selects a minimizer the first pass did not count.
    pub fn insert(&mut self, kmer_hash: u64, record: PositionRecord) -> bool {
        let partition_index = (kmer_hash & (PARTITION_COUNT as u64 - 1)) as usize;
        let partition = &mut self.partitions[partition_index];
        let bucket = match partition.table.get(kmer_hash & !MAXIMUM_COUNT) {
            Some(bucket) => bucket,
            None => return false,
        };

        let offset = partition.table.value_at(bucket);
        let count = partition.table.key_at(bucket) & MAXIMUM_COUNT;
        if count >= MAXIMUM_COUNT {
            return false;
        }

        partition.occurrences[(offset + count) as usize] = record;
        let incremented = partition.table.key_at(bucket) + 1;
        partition.table.set_key_at(bucket, incremented);
        self.total_positions += 1;
        true
    }

    /// `ha_pt_get`: the occurrence list of `kmer_hash`.
    pub fn get(&self, kmer_hash: u64) -> &[PositionRecord] {
        let partition_index = (kmer_hash & (PARTITION_COUNT as u64 - 1)) as usize;
        let partition = &self.partitions[partition_index];
        match partition.table.get(kmer_hash & !MAXIMUM_COUNT) {
            None => &[],
            Some(bucket) => {
                let offset = partition.table.value_at(bucket) as usize;
                let count = (partition.table.key_at(bucket) & MAXIMUM_COUNT) as usize;
                &partition.occurrences[offset..offset + count]
            }
        }
    }

    /// Number of distinct k-mers in the index, `ha_pt_t::tot`.
    pub fn distinct_count(&self) -> u64 {
        self.partitions
            .iter()
            .map(|partition| partition.table.len() as u64)
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn occurrences_are_stored_in_insertion_order() {
        let mut counts = KmerCountTable::new(3, 0, 0);
        for _ in 0..3 {
            counts.insert(500);
        }
        for _ in 0..2 {
            counts.insert(1500);
        }

        let mut index = PositionIndex::build(&counts, 3);
        assert_eq!(index.capacity(), 5);

        assert!(index.insert(
            500,
            PositionRecord {
                read_identifier: 7,
                position: 42,
                reverse: false,
                span: 3,
            }
        ));
        assert!(index.insert(
            500,
            PositionRecord {
                read_identifier: 8,
                position: 43,
                reverse: true,
                span: 3,
            }
        ));

        let occurrences = index.get(500);
        assert_eq!(occurrences.len(), 2);
        assert_eq!(occurrences[0].read_identifier, 7);
        assert_eq!(occurrences[1].read_identifier, 8);
        assert_eq!(index.get(1500).len(), 0);
        assert_eq!(index.total_positions(), 2);
    }

    #[test]
    fn unknown_hashes_are_ignored() {
        let counts = KmerCountTable::new(3, 0, 0);
        let mut index = PositionIndex::build(&counts, 3);
        assert!(!index.insert(12345, PositionRecord::default()));
        assert!(index.get(12345).is_empty());
    }
}

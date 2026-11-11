//! Position index, the counterpart of `ha_pt_t` in `htab.cpp`.
//!
//! For every minimizer hash the index stores the list of places where it
//! occurs: read identifier, position, strand and span. The list of a hash is
//! filled in sketch order, which is what later stages rely on when they walk
//! the occurrences of a k-mer.

use crate::container::khash::{hash_without_counter_bits, U64HashTable};
use crate::counting::filter_table::flag_size;
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
        Self::build_selected(counts, kmer_length, None)
    }

    /// `ha_pt_gen` restricted to the k-mers the `--fast` selection kept.
    pub fn build_selected(
        counts: &KmerCountTable,
        kmer_length: i32,
        keep: Option<&std::collections::HashSet<u64>>,
    ) -> Self {
        let mut partitions: Vec<IndexPartition> = (0..PARTITION_COUNT)
            .map(|_| IndexPartition {
                table: new_index_table(),
                occurrences: Vec::new(),
            })
            .collect();

        for (partition_index, partition) in partitions.iter_mut().enumerate() {
            let mut next_start = 0u64;
            for (key, _) in counts.partition_entries(partition_index) {
                let hash = (key >> COUNTER_BITS << COUNTER_BITS) | partition_index as u64;
                if let Some(keep) = keep {
                    if !keep.contains(&hash) {
                        continue;
                    }
                }
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

    /// The partition count, `1 << ha_pt_t::pre`.
    fn partition_count(&self) -> u32 {
        self.partitions.len() as u32
    }

    /// `uidx_write`'s position-index half (`htab.cpp` 1573): the header and
    /// then, per partition, the `khashl` map dump plus the occurrence list.
    pub fn write(&self, writer: &mut impl std::io::Write) -> std::io::Result<()> {
        writer.write_all(&self.kmer_length.to_le_bytes())?;
        writer.write_all(&self.partition_count().to_le_bytes())?;
        writer.write_all(&self.distinct_count().to_le_bytes())?;
        writer.write_all(&self.total_positions.to_le_bytes())?;
        for partition in &self.partitions {
            write_partition(partition, writer)?;
        }
        Ok(())
    }

    /// `uidx_load`'s position-index half.
    pub fn read(reader: &mut impl std::io::Read) -> std::io::Result<Self> {
        let kmer_length = read_i32_le(reader)?;
        let partition_count = read_i32_le(reader)?;
        let _tot = read_u64_le(reader)?;
        let total_positions = read_u64_le(reader)?;
        let mut partitions = Vec::with_capacity(partition_count as usize);
        for _ in 0..partition_count {
            partitions.push(read_partition(reader)?);
        }
        Ok(Self {
            kmer_length,
            partitions,
            total_positions,
        })
    }
}

fn read_i32_le(reader: &mut impl std::io::Read) -> std::io::Result<i32> {
    let mut buffer = [0u8; 4];
    reader.read_exact(&mut buffer)?;
    Ok(i32::from_le_bytes(buffer))
}

fn read_u64_le(reader: &mut impl std::io::Read) -> std::io::Result<u64> {
    let mut buffer = [0u8; 8];
    reader.read_exact(&mut buffer)?;
    Ok(u64::from_le_bytes(buffer))
}

/// One partition's `yak_pt_save` dump plus `g->n` and `g->al`.
fn write_partition(
    partition: &IndexPartition,
    writer: &mut impl std::io::Write,
) -> std::io::Result<()> {
    let bucket_count = partition.table.bucket_count();
    writer.write_all(&bucket_count.to_le_bytes())?;
    writer.write_all(&partition.table.bits().to_le_bytes())?;
    writer.write_all(&(partition.table.len() as u32).to_le_bytes())?;
    let allocated = bucket_count != 0;
    writer.write_all(&[u8::from(allocated)])?;
    if allocated {
        for word in partition.table.used_words() {
            writer.write_all(&word.to_le_bytes())?;
        }
        writer.write_all(&[1])?;
        for bucket in 0..bucket_count {
            writer.write_all(&partition.table.key_at(bucket).to_le_bytes())?;
            writer.write_all(&partition.table.value_at(bucket).to_le_bytes())?;
        }
    } else {
        writer.write_all(&[0])?;
    }
    // `g->n`: the occurrence array was allocated to the summed counts of the
    // counting pass, which is exactly the vector length here.
    writer.write_all(&(partition.occurrences.len() as u64).to_le_bytes())?;
    for record in &partition.occurrences {
        // `ha_idxposl_t`: `rid:31 | rev:1 | pos:32` packed into one `u64`,
        // then the one-byte span and the seven zero padding bytes of the
        // reference's `calloc`.
        let packed = u64::from(record.read_identifier & 0x7fff_ffff)
            | (u64::from(u8::from(record.reverse)) << 31)
            | (u64::from(record.position) << 32);
        writer.write_all(&packed.to_le_bytes())?;
        writer.write_all(&[record.span as u8])?;
        writer.write_all(&[0u8; 7])?;
    }
    Ok(())
}

fn read_partition(reader: &mut impl std::io::Read) -> std::io::Result<IndexPartition> {
    let bucket_count = read_i32_le(reader)? as u32;
    let bits = read_i32_le(reader)? as u32;
    let _count = read_i32_le(reader)?;
    let mut used_words = vec![0u32; flag_size(bucket_count)];
    let mut flag = [0u8; 1];
    reader.read_exact(&mut flag)?;
    if flag[0] != 0 {
        for word in used_words.iter_mut() {
            *word = read_i32_le(reader)? as u32;
        }
    } else {
        used_words.clear();
    }
    reader.read_exact(&mut flag)?;
    let mut keys = Vec::new();
    let mut values = Vec::new();
    if flag[0] != 0 {
        keys.reserve(bucket_count as usize);
        values.reserve(bucket_count as usize);
        for _ in 0..bucket_count {
            keys.push(read_u64_le(reader)?);
            values.push(read_u64_le(reader)?);
        }
    }
    let table = U64HashTable::from_raw_parts(
        bits,
        used_words,
        keys,
        values,
        index_hash,
        COUNTER_BITS,
    );
    let length = read_u64_le(reader)? as usize;
    let mut occurrences = Vec::with_capacity(length);
    let mut record = [0u8; 16];
    for _ in 0..length {
        reader.read_exact(&mut record)?;
        let packed = u64::from_le_bytes(record[..8].try_into().unwrap());
        occurrences.push(PositionRecord {
            read_identifier: (packed & 0x7fff_ffff) as u32,
            reverse: (packed >> 31) & 1 != 0,
            position: (packed >> 32) as u32,
            span: u32::from(record[8]),
        });
    }
    Ok(IndexPartition {
        table,
        occurrences,
    })
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

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

    /// `yak_ft_save` (`khashl.h` `__KHASHL_IMPL_S_L` on the packed
    /// `{key, int16}` buckets): the byte-exact on-disk form.
    pub fn write(&self, writer: &mut impl std::io::Write) -> std::io::Result<()> {
        write_hash_table(&self.table, writer, |value| {
            let value = i16::try_from(value).unwrap_or(i16::MAX);
            value.to_le_bytes().to_vec()
        })
    }

    /// `yak_ft_load`: the inverse of [`Self::write`].
    pub fn read(reader: &mut impl std::io::Read) -> std::io::Result<Self> {
        let table = read_hash_table(reader, |bytes| {
            i64::from(i16::from_le_bytes([bytes[8], bytes[9]])) as u64
        })?;
        Ok(Self { table })
    }
}

/// `__kh_fsize(m)`: words of the used bitmap.
pub(crate) fn flag_size(bucket_count: u32) -> usize {
    if bucket_count < 32 {
        1
    } else {
        (bucket_count >> 5) as usize
    }
}

/// The `khashl` map dump both `uidx_write` tables share: counts, the used
/// bitmap and then the packed `{key, value}` buckets in bucket order.
fn write_hash_table(
    table: &U64HashTable,
    writer: &mut impl std::io::Write,
    encode_value: impl Fn(u64) -> Vec<u8>,
) -> std::io::Result<()> {
    let bucket_count = table.bucket_count();
    writer.write_all(&bucket_count.to_le_bytes())?;
    writer.write_all(&table.bits().to_le_bytes())?;
    writer.write_all(&(table.len() as u32).to_le_bytes())?;
    // A never-allocated table writes neither array, exactly like the
    // reference's NULL `used`/`keys` pointers.
    let allocated = bucket_count != 0;
    writer.write_all(&[u8::from(allocated)])?;
    if allocated {
        for word in table.used_words() {
            writer.write_all(&word.to_le_bytes())?;
        }
        writer.write_all(&[1])?;
        for bucket in 0..bucket_count {
            let mut bucket_bytes = table.key_at(bucket).to_le_bytes().to_vec();
            bucket_bytes.extend_from_slice(&encode_value(table.value_at(bucket)));
            writer.write_all(&bucket_bytes)?;
        }
    } else {
        writer.write_all(&[0])?;
    }
    Ok(())
}

fn read_u32_le(reader: &mut impl std::io::Read) -> std::io::Result<u32> {
    let mut buffer = [0u8; 4];
    reader.read_exact(&mut buffer)?;
    Ok(u32::from_le_bytes(buffer))
}

/// The inverse of [`write_hash_table`]; `decode_value` reads the value part of
/// one packed bucket.
fn read_hash_table(
    reader: &mut impl std::io::Read,
    decode_value: impl Fn(&[u8]) -> u64,
) -> std::io::Result<U64HashTable> {
    let bucket_count = read_u32_le(reader)?;
    let bits = read_u32_le(reader)?;
    let _count = read_u32_le(reader)?;
    let mut used_words = vec![0u32; flag_size(bucket_count)];
    let mut flag = [0u8; 1];
    reader.read_exact(&mut flag)?;
    if flag[0] != 0 {
        for word in used_words.iter_mut() {
            *word = read_u32_le(reader)?;
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
        let mut bucket = [0u8; 16];
        for _ in 0..bucket_count {
            reader.read_exact(&mut bucket[..10])?;
            keys.push(u64::from_le_bytes(bucket[..8].try_into().unwrap()));
            values.push(decode_value(&bucket));
        }
    }
    Ok(U64HashTable::from_raw_parts(
        bits,
        used_words,
        keys,
        values,
        identity_hash,
        0,
    ))
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

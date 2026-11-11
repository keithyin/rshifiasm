//! A faithful port of `khashl.h`, the hash table used by the original
//! implementation for k-mer counts, position indexes and high-occurrence
//! filter tables.
//!
//! The layout matters: the bucket index of a key is derived with the same
//! Fibonacci hashing (`hash * 2654435769 >> (32 - bits)`), probing is linear,
//! the table grows at a 0.75 load factor, and iteration walks buckets from
//! zero upwards. Reproducing all of that keeps iteration order identical to
//! the C code, which several output orders depend on.

/// Load factor limit, `(n_buckets >> 1) + (n_buckets >> 2)`.
fn load_limit(bucket_count: u32) -> u32 {
    (bucket_count >> 1) + (bucket_count >> 2)
}

/// `__kh_h2b`: Fibonacci hashing of the 32-bit hash into `bits` buckets.
pub fn bucket_for_hash(hash: u32, bits: u32) -> u32 {
    hash.wrapping_mul(2_654_435_769) >> (32 - bits)
}

/// Used-bucket bitmap helpers (`__kh_used` and friends).
struct UsedBitmap {
    words: Vec<u32>,
}

impl UsedBitmap {
    fn new(bucket_count: u32) -> Self {
        let word_count = if bucket_count < 32 {
            1
        } else {
            (bucket_count >> 5) as usize
        };
        Self {
            words: vec![0u32; word_count],
        }
    }

    fn is_used(&self, bucket: u32) -> bool {
        (self.words[(bucket >> 5) as usize] >> (bucket & 0x1f)) & 1 != 0
    }

    fn set_used(&mut self, bucket: u32) {
        self.words[(bucket >> 5) as usize] |= 1u32 << (bucket & 0x1f);
    }

    fn set_unused(&mut self, bucket: u32) {
        self.words[(bucket >> 5) as usize] &= !(1u32 << (bucket & 0x1f));
    }
}

/// Open addressing map from `u64` keys to `u64` values, the concrete shape
/// used by the counting pipeline.
pub struct U64HashTable {
    /// log2 of the bucket count; zero means "not allocated yet".
    bits: u32,
    count: u32,
    used: UsedBitmap,
    keys: Vec<u64>,
    values: Vec<u64>,
    hash_function: fn(u64) -> u32,
    /// Lowest bits ignored by the equality predicate, mirroring the
    /// `yak_ct_eq` macro where the counter is packed into the key.
    ignored_low_bits: u32,
}

impl U64HashTable {
    /// `prefix_init`
    pub fn new(hash_function: fn(u64) -> u32) -> Self {
        Self::with_ignored_low_bits(hash_function, 0)
    }

    /// Builds a table whose equality predicate ignores the low `bits` bits of
    /// the key, as the counting table does for its packed counter.
    pub fn with_ignored_low_bits(hash_function: fn(u64) -> u32, bits: u32) -> Self {
        Self {
            bits: 0,
            count: 0,
            used: UsedBitmap::new(0),
            keys: Vec::new(),
            values: Vec::new(),
            hash_function,
            ignored_low_bits: bits,
        }
    }

    fn keys_equal(&self, left: u64, right: u64) -> bool {
        if self.ignored_low_bits == 0 {
            left == right
        } else {
            (left >> self.ignored_low_bits) == (right >> self.ignored_low_bits)
        }
    }

    /// `kh_size`
    pub fn len(&self) -> usize {
        self.count as usize
    }

    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Number of allocated buckets, `1 << bits`.
    pub fn bucket_count(&self) -> u32 {
        if self.keys.is_empty() {
            0
        } else {
            1u32 << self.bits
        }
    }

    /// `prefix_getp`: bucket of `key`, or `None`.
    pub fn get(&self, key: u64) -> Option<u32> {
        if self.keys.is_empty() {
            return None;
        }

        let bucket_count = 1u32 << self.bits;
        let mask = bucket_count - 1;
        let mut bucket = bucket_for_hash((self.hash_function)(key), self.bits);
        let first_probe = bucket;

        while self.used.is_used(bucket) && !self.keys_equal(self.keys[bucket as usize], key) {
            bucket = (bucket + 1) & mask;
            if bucket == first_probe {
                return None;
            }
        }

        if self.used.is_used(bucket) {
            Some(bucket)
        } else {
            None
        }
    }

    pub fn get_value(&self, key: u64) -> Option<u64> {
        self.get(key).map(|bucket| self.values[bucket as usize])
    }

    /// `prefix_putp`: returns the bucket for `key` and whether it is new.
    pub fn put(&mut self, key: u64) -> (u32, bool) {
        if self.keys.is_empty() || self.count + 1 > load_limit(1u32 << self.bits) {
            let requested = if self.keys.is_empty() {
                4
            } else {
                (1u32 << self.bits) + 1
            };
            self.resize(requested);
        }

        let bucket_count = 1u32 << self.bits;
        let mask = bucket_count - 1;
        let mut bucket = bucket_for_hash((self.hash_function)(key), self.bits);
        let first_probe = bucket;
        let mut absent = true;

        while self.used.is_used(bucket) {
            if self.keys_equal(self.keys[bucket as usize], key) {
                absent = false;
                break;
            }
            bucket = (bucket + 1) & mask;
            if bucket == first_probe {
                break;
            }
        }

        if absent {
            self.used.set_used(bucket);
            self.keys[bucket as usize] = key;
            self.values[bucket as usize] = 0;
            self.count += 1;
        }

        (bucket, absent)
    }

    pub fn key_at(&self, bucket: u32) -> u64 {
        self.keys[bucket as usize]
    }

    pub fn value_at(&self, bucket: u32) -> u64 {
        self.values[bucket as usize]
    }

    pub fn set_value_at(&mut self, bucket: u32, value: u64) {
        self.values[bucket as usize] = value;
    }

    /// Overwrites the key of an existing bucket.
    ///
    /// The counting table packs the occurrence count into the low bits of the
    /// key, so incrementing a counter means rewriting the key. This is safe
    /// because the hash and the equality predicate of that table both ignore
    /// those low bits.
    pub fn set_key_at(&mut self, bucket: u32, key: u64) {
        self.keys[bucket as usize] = key;
    }

    /// `prefix_del`
    pub fn remove(&mut self, bucket: u32) {
        if self.keys.is_empty() || !self.used.is_used(bucket) {
            return;
        }
        self.used.set_unused(bucket);
        self.count -= 1;
    }

    /// Whether `bucket` currently holds an entry; used when a bucket index
    /// itself participates in the algorithm.
    pub fn bucket_is_used(&self, bucket: u32) -> bool {
        !self.keys.is_empty() && bucket < self.bucket_count() && self.used.is_used(bucket)
    }

    /// `log2` of the bucket count, the `bits` field of the reference table.
    pub fn bits(&self) -> u32 {
        self.bits
    }

    /// Rebuilds a table from its serialized parts (`prefix##_load`).
    ///
    /// `used_words` holds the packed used-bucket bitmap, `keys` and `values`
    /// the per-bucket entries in bucket order. An empty `keys` slice
    /// reproduces a never-allocated table.
    pub fn from_raw_parts(
        bits: u32,
        used_words: Vec<u32>,
        keys: Vec<u64>,
        values: Vec<u64>,
        hash_function: fn(u64) -> u32,
        ignored_low_bits: u32,
    ) -> Self {
        let count = used_words
            .iter()
            .map(|word| word.count_ones())
            .sum::<u32>();
        Self {
            bits,
            count,
            used: UsedBitmap { words: used_words },
            keys,
            values,
            hash_function,
            ignored_low_bits,
        }
    }

    /// The used-bucket bitmap in word order, `h->used[0..__kh_fsize(n)]`.
    pub fn used_words(&self) -> &[u32] {
        &self.used.words
    }

    /// Iterates over live entries in ascending bucket order, matching the
    /// `for (k = 0; k < kh_end(h); ++k)` loops of the original code.
    pub fn iter(&self) -> impl Iterator<Item = (u64, u64)> + '_ {
        (0..self.bucket_count())
            .filter(move |bucket| self.used.is_used(*bucket))
            .map(move |bucket| (self.keys[bucket as usize], self.values[bucket as usize]))
    }

    /// `prefix_resize`
    pub fn resize(&mut self, requested_bucket_count: u32) {
        // `kroundup32` of zero stays zero and is then floored to four buckets.
        let requested_bucket_count = requested_bucket_count.max(1);
        let mut bits = 0u32;
        let mut remaining = requested_bucket_count;
        while {
            remaining >>= 1;
            remaining != 0
        } {
            bits += 1;
        }
        if requested_bucket_count & (requested_bucket_count - 1) != 0 {
            bits += 1;
        }
        if bits < 2 {
            bits = 2;
        }
        let new_bucket_count = 1u32 << bits;

        if self.count > load_limit(new_bucket_count) {
            return;
        }

        let mut new_used = UsedBitmap::new(new_bucket_count);
        let old_bucket_count = self.bucket_count();

        if old_bucket_count < new_bucket_count {
            self.keys.resize(new_bucket_count as usize, 0);
            self.values.resize(new_bucket_count as usize, 0);
        }

        let new_mask = new_bucket_count - 1;
        for old_bucket in 0..old_bucket_count {
            if !self.used.is_used(old_bucket) {
                continue;
            }

            let mut key = self.keys[old_bucket as usize];
            let mut value = self.values[old_bucket as usize];
            self.used.set_unused(old_bucket);

            loop {
                let mut probe = bucket_for_hash((self.hash_function)(key), bits);
                while new_used.is_used(probe) {
                    probe = (probe + 1) & new_mask;
                }
                new_used.set_used(probe);

                if probe < old_bucket_count && self.used.is_used(probe) {
                    std::mem::swap(&mut self.keys[probe as usize], &mut key);
                    std::mem::swap(&mut self.values[probe as usize], &mut value);
                    self.used.set_unused(probe);
                } else {
                    self.keys[probe as usize] = key;
                    self.values[probe as usize] = value;
                    break;
                }
            }
        }

        if old_bucket_count > new_bucket_count {
            self.keys.resize(new_bucket_count as usize, 0);
            self.values.resize(new_bucket_count as usize, 0);
        }

        self.bits = bits;
        self.used = new_used;
    }
}

/// Identity hash for keys that already are hashes (`kh_hash_dummy`).
pub fn identity_hash(key: u64) -> u32 {
    key as u32
}

/// Hash that ignores the counter packed in the low bits (`yak_ct_hash`).
pub fn hash_without_counter_bits(key: u64, counter_bits: u32) -> u32 {
    (key >> counter_bits) as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fibonacci_hashing_matches_the_reference_formula() {
        assert_eq!(bucket_for_hash(0, 4), 0);
        assert_eq!(bucket_for_hash(1, 4), (1u32.wrapping_mul(2_654_435_769)) >> 28);
        assert_eq!(bucket_for_hash(0xffff_ffff, 2), 1);
    }

    #[test]
    fn insertion_and_lookup_are_consistent_across_resizes() {
        let mut table = U64HashTable::new(identity_hash);
        for key in 0..5000u64 {
            let (bucket, inserted) = table.put(key * 7);
            assert!(inserted);
            table.set_value_at(bucket, key);
        }

        assert_eq!(table.len(), 5000);
        for key in 0..5000u64 {
            let bucket = table.get(key * 7).expect("key must be present");
            assert_eq!(table.value_at(bucket), key);
        }
        assert!(table.get(999_999).is_none());
    }

    #[test]
    fn duplicate_insertion_reports_the_existing_bucket() {
        let mut table = U64HashTable::new(identity_hash);
        let (first, inserted) = table.put(42);
        assert!(inserted);
        let (second, inserted_again) = table.put(42);
        assert!(!inserted_again);
        assert_eq!(first, second);
        assert_eq!(table.len(), 1);
    }

    #[test]
    fn iteration_is_stable_and_visits_every_entry_once() {
        let mut table = U64HashTable::new(identity_hash);
        for key in [10u64, 3, 77, 1] {
            table.put(key);
        }

        let first_pass: Vec<u64> = table.iter().map(|(key, _)| key).collect();
        let second_pass: Vec<u64> = table.iter().map(|(key, _)| key).collect();
        assert_eq!(first_pass, second_pass);

        let mut sorted = first_pass.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![1, 3, 10, 77]);
    }

    #[test]
    fn bucket_layout_matches_the_reference_khashl() {
        // Golden layout produced by the original `khashl.h` for twenty keys
        // `0, 3, ..., 57` inserted with an identity hash:
        //   n_buckets=32 order: 0@0 18@3 39@4 36@7 15@8 57@9 54@11 33@12
        //   12@13 51@16 9@17 30@18 27@21 6@22 48@23 45@25 24@26 3@27
        //   42@30 21@31
        let mut table = U64HashTable::new(identity_hash);
        for key in 0..20u64 {
            table.put(key * 3);
        }

        assert_eq!(table.bucket_count(), 32);
        assert_eq!(table.len(), 20);

        let layout: Vec<String> = (0..table.bucket_count())
            .filter(|bucket| table.bucket_is_used(*bucket))
            .map(|bucket| format!("{}@{}", table.key_at(bucket), bucket))
            .collect();
        assert_eq!(
            layout.join(" "),
            "0@0 18@3 39@4 36@7 15@8 57@9 54@11 33@12 12@13 51@16 9@17 30@18 27@21 6@22 48@23 45@25 24@26 3@27 42@30 21@31"
        );
    }

    #[test]
    fn removal_keeps_the_remaining_keys_reachable() {
        let mut table = U64HashTable::new(identity_hash);
        for key in 0..100u64 {
            table.put(key);
        }
        for key in 0..50u64 {
            let bucket = table.get(key).unwrap();
            table.remove(bucket);
        }
        assert_eq!(table.len(), 50);
        for key in 0..50u64 {
            assert!(table.get(key).is_none());
        }
        for key in 50..100u64 {
            assert!(table.get(key).is_some());
        }
    }
}

//! The 64-bit mixing function used to hash k-mers, `yak_hash64_64`.

/// Thomas Wang's 64-bit integer hash, as implemented in `htab.h`.
///
/// Everything downstream (counting, the bloom filter, the filter table) is
/// keyed by this value, so the exact sequence of shifts and additions matters.
pub fn hash64_64(key: u64) -> u64 {
    let mut key = key;
    key = (!key).wrapping_add(key << 21);
    key ^= key >> 24;
    key = key.wrapping_add(key << 3).wrapping_add(key << 8);
    key ^= key >> 14;
    key = key.wrapping_add(key << 2).wrapping_add(key << 4);
    key ^= key >> 28;
    key.wrapping_add(key << 31)
}

/// `yak_hash_long`: hashes the canonical orientation of a two-word k-mer.
///
/// `kmer[0..2]` is the forward encoding, `kmer[2..4]` the reverse complement
/// encoding. The strand is chosen with the same comparison the original uses.
pub fn hash_long_kmer(kmer: &[u64; 4]) -> u64 {
    let strand = if kmer[1] < kmer[3] { 0 } else { 1 };
    hash64_64(kmer[strand << 1]).wrapping_add(hash64_64(kmer[(strand << 1) | 1]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hash64_64_matches_the_reference_values() {
        // Values printed by the C implementation for the same inputs.
        assert_eq!(hash64_64(0), 8_633_297_058_295_171_728);
        assert_eq!(hash64_64(1), 6_614_235_796_240_398_542);
        assert_eq!(hash64_64(2), 13_228_483_051_453_548_148);
    }

    #[test]
    fn hash_long_kmer_chooses_the_canonical_strand() {
        let forward = [1u64, 2, 3, 4];
        let reversed = [4u64, 3, 1, 2];
        assert_eq!(hash_long_kmer(&forward), 1_395_974_773_984_395_074);
        assert_eq!(hash_long_kmer(&forward), hash_long_kmer(&reversed));
    }
}

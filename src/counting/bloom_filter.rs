//! Blocked bloom filter, the counterpart of `yak_bf_t` in `htab.cpp`.
//!
//! The filter is not an optimisation detail: the counting code only inserts a
//! k-mer once every probe bit is already set, which shifts the recorded counts
//! of rare k-mers. Reproducing the bit layout is therefore required for the
//! histogram (and the inferred coverage) to match the reference.

/// log2 of the block size in bits, `YAK_BLK_SHIFT` (64 bytes per block).
const BLOCK_SHIFT: u32 = 9;

/// Bits per block.
const BLOCK_BITS: u32 = 1 << BLOCK_SHIFT;

/// Mask selecting a bit inside one block.
const BLOCK_MASK: u32 = BLOCK_BITS - 1;

/// `yak_bf_t`: a blocked bloom filter with `hash_count` probes per insertion.
pub struct BlockedBloomFilter {
    /// log2 of the filter size in bits.
    shift: u32,
    /// Number of probe bits, four by default.
    hash_count: u32,
    /// Bit array, allocated on first use like the lazily touched partitions of
    /// the reference implementation.
    bits: Option<Vec<u8>>,
}

impl BlockedBloomFilter {
    /// `yak_bf_init`; returns `None` when the size is out of range.
    pub fn new(shift: u32, hash_count: u32) -> Option<Self> {
        if shift + BLOCK_SHIFT > 64 || shift < BLOCK_SHIFT {
            return None;
        }
        Some(Self {
            shift,
            hash_count,
            bits: None,
        })
    }

    /// Number of probe bits per insertion.
    pub fn hash_count(&self) -> u32 {
        self.hash_count
    }

    /// `yak_bf_insert`: sets the probe bits and returns how many of them were
    /// already set.
    pub fn insert(&mut self, hash: u64) -> u32 {
        let block_index_shift = self.shift - BLOCK_SHIFT;
        let block_index = (hash & ((1u64 << block_index_shift) - 1)) as usize;
        let first_probe = ((hash >> block_index_shift) & BLOCK_MASK as u64) as u32;
        let mut step = ((hash >> self.shift) & BLOCK_MASK as u64) as u32;

        // Avoid repeatedly probing the same few bits.
        if step & 31 == 0 {
            step = (step + 1) & BLOCK_MASK;
        }

        let byte_count = 1usize << (self.shift - 3);
        let bits = self
            .bits
            .get_or_insert_with(|| vec![0u8; byte_count]);

        let block = &mut bits[block_index * 64..block_index * 64 + 64];

        let mut already_set = 0;
        let mut probe = first_probe;
        for _ in 0..self.hash_count {
            let byte = &mut block[(probe >> 3) as usize];
            let mask = 1u8 << (probe & 7);
            if *byte & mask != 0 {
                already_set += 1;
            }
            *byte |= mask;
            probe = (probe + step) & BLOCK_MASK;
        }

        already_set
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repeated_insertion_reports_every_bit_as_set() {
        let mut filter = BlockedBloomFilter::new(25, 4).unwrap();
        let first = filter.insert(0x1234_5678);
        assert!(first <= 4);

        let second = filter.insert(0x1234_5678);
        assert_eq!(second, 4, "a repeated hash must report every bit as set");
    }

    #[test]
    fn out_of_range_sizes_are_rejected() {
        assert!(BlockedBloomFilter::new(8, 4).is_none());
        assert!(BlockedBloomFilter::new(60, 4).is_none());
        assert!(BlockedBloomFilter::new(25, 4).is_some());
    }
}

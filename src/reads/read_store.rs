//! In-memory read storage, the counterpart of `All_reads` (`Process_Read.h`).
//!
//! The original structure keeps every array in a hand managed C buffer with
//! parallel `index`/`size` counters. Here the same information is held in
//! `Vec`s, which keeps the indexing rules (`read_length`, `name_index`) while
//! removing the manual growth bookkeeping.

use crate::command_line::options::AssemblyOptions;
use crate::sequence::encoding::PackedSequence;
use crate::sequence::fastx_reader::FastxReader;

/// Label assigned to a read before trio binning (`AMBIGU`).
pub const TRIO_LABEL_UNKNOWN: u8 = 0;

/// `sc_bn`: the number of bits per ONT quality value.
pub const QUALITY_BITS: u32 = 2;

/// `sc_tb`: the Phred score of every 2-bit quality level.
pub const QUALITY_THRESHOLDS: [u8; 8] = [10, 20, 30, 40, 50, 60, 70, 80];

/// `ha_compress_qual`: packs `(quality - offset)` into 2-bit levels, four
/// bases per byte with the first base in the high bits and the tail padded on
/// the right, exactly like the reference.
pub fn compress_quality(quality: &[u8], offset: u8) -> Vec<u8> {
    let mut packed = Vec::with_capacity(quality.len().div_ceil(4));
    let mut byte = 0u8;
    let mut count = 0u32;
    for &value in quality {
        byte = (byte << QUALITY_BITS) | quality_level(value, offset);
        count += 1;
        if count == 4 {
            packed.push(byte);
            byte = 0;
            count = 0;
        }
    }
    if count > 0 {
        packed.push(byte << (8 - count * QUALITY_BITS));
    }
    packed
}

fn quality_level(value: u8, offset: u8) -> u8 {
    let score = value.wrapping_sub(offset);
    let levels = 1u8 << QUALITY_BITS;
    let mut level = 0u8;
    while level < levels && QUALITY_THRESHOLDS[level as usize] < score {
        level += 1;
    }
    level.min(levels - 1)
}

/// `ha_compress_qual_bit`: packs already-quantised 2-bit levels, four per
/// byte with the first position in the high bits and the tail padded on the
/// right.
pub fn pack_quality_levels(levels: &[u8]) -> Vec<u8> {
    let mut packed = Vec::with_capacity(levels.len().div_ceil(4));
    let mut byte = 0u8;
    let mut count = 0u32;
    for &level in levels {
        byte = (byte << QUALITY_BITS) | (level & 3);
        count += 1;
        if count == 4 {
            packed.push(byte);
            byte = 0;
            count = 0;
        }
    }
    if count > 0 {
        packed.push(byte << (8 - count * QUALITY_BITS));
    }
    packed
}

/// `retrive_bqual`: the 2-bit levels of `[start, end)`, reversed when
/// `reverse` is set. Positions past the packed length decode as zero.
pub fn retrieve_quality(packed: &[u8], start: i64, end: i64, reverse: bool) -> Vec<u8> {
    let mut levels = Vec::with_capacity((end - start).max(0) as usize);
    for position in start..end {
        if position < 0 {
            levels.push(0);
            continue;
        }
        let byte = position as usize / 4;
        let shift = 6 - (position as usize % 4) * 2;
        levels.push(packed.get(byte).map(|value| (value >> shift) & 3).unwrap_or(0));
    }
    if reverse {
        levels.reverse();
    }
    levels
}

/// A corpus of reads, ordered exactly as they appear in the input files.
#[derive(Debug, Clone, Default)]
pub struct ReadStore {
    /// Length of every kept read, after adapter trimming.
    pub read_lengths: Vec<u64>,
    /// Names of every kept read, without the leading `@`/`>`.
    pub read_names: Vec<String>,
    /// Packed sequences, filled in the second pass of the counting pipeline.
    pub sequences: Vec<PackedSequence>,
    /// Trio labels, `TRIO_LABEL_UNKNOWN` until trio binning runs.
    pub trio_labels: Vec<u8>,
    /// Packed 2-bit ONT qualities (`rsc`), one entry per read; empty when the
    /// run is not in `is_sc` mode.
    pub qualities: Vec<Vec<u8>>,
    /// Number of kept reads.
    pub total_reads: u64,
    /// Sum of `read_lengths`.
    pub total_bases: u64,
}

impl ReadStore {
    /// Reads every input file and stores the kept reads in input order,
    /// mirroring the filtering performed by `sf_worker_count` in `htab.cpp`.
    pub fn load_from_files(options: &AssemblyOptions) -> std::io::Result<Self> {
        let mut store = Self::default();

        for file_name in &options.input_file_names {
            let mut reader = FastxReader::open(file_name)?;
            while let Some(record) = reader.next_record()? {
                let adapter_length = options.adapter_length.max(0) as usize;
                if record.sequence.len() < adapter_length * 2 {
                    continue;
                }

                let trimmed_length = record.sequence.len() - adapter_length * 2;
                if trimmed_length == 0 {
                    continue;
                }
                if (trimmed_length as i64) < options.read_length_cutoff {
                    continue;
                }

                let trimmed_sequence = &record.sequence[adapter_length..adapter_length + trimmed_length];
                let quality = if options.single_cell_mode != 0 {
                    // `flt_quals(qual + adapterLen, l, 33, sc_cut)`: the
                    // mean-quality filter of the ONT preset.
                    let trimmed_quality =
                        &record.quality[adapter_length..adapter_length + trimmed_length];
                    if options.base_quality_cutoff > 0
                        && !quality_passes(
                            trimmed_quality,
                            33,
                            options.base_quality_cutoff,
                        )
                    {
                        continue;
                    }
                    Some(trimmed_quality)
                } else {
                    None
                };
                store.push_read_with_quality(&record.name, trimmed_sequence, quality);
            }
        }

        Ok(store)
    }

    /// Appends one read, the equivalent of `ha_insert_read_len` followed by
    /// the packed sequence insertion.
    pub fn push_read(&mut self, name: &str, sequence: &[u8]) {
        self.push_read_with_quality(name, sequence, None);
    }

    /// [`Self::push_read`] with an optional FASTQ quality string (already
    /// trimmed), packed into the `rsc` 2-bit levels.
    pub fn push_read_with_quality(&mut self, name: &str, sequence: &[u8], quality: Option<&[u8]>) {
        self.read_lengths.push(sequence.len() as u64);
        self.read_names.push(name.to_string());
        self.sequences.push(PackedSequence::from_ascii(sequence));
        self.trio_labels.push(TRIO_LABEL_UNKNOWN);
        self.qualities
            .push(quality.map(|values| compress_quality(values, 33)).unwrap_or_default());
        self.total_reads += 1;
        self.total_bases += sequence.len() as u64;
    }

    /// `Get_READ_LENGTH`.
    pub fn read_length(&self, read_id: usize) -> u64 {
        self.read_lengths[read_id]
    }

    /// `Get_NAME`.
    pub fn read_name(&self, read_id: usize) -> &str {
        &self.read_names[read_id]
    }

    /// `recover_UC_Read`: the read sequence in plain ASCII.
    pub fn read_sequence(&self, read_id: usize) -> Vec<u8> {
        self.sequences[read_id].unpack()
    }

    /// `recover_UC_Read_RC`: the reverse complement of the read sequence.
    pub fn read_sequence_reverse_complement(&self, read_id: usize) -> Vec<u8> {
        self.sequences[read_id].unpack_reverse_complement()
    }

    /// `recover_UC_Read_sub_region`: a slice of a read, in the requested
    /// orientation.
    pub fn read_sub_region(
        &self,
        read_id: usize,
        start_position: i64,
        length: i64,
        strand: u8,
    ) -> Vec<u8> {
        self.sequences[read_id].sub_region(start_position, length, strand == 1)
    }

    /// `retrive_bqual`: the 2-bit quality levels of `[start, end)` in the
    /// requested orientation.
    pub fn read_quality(&self, read_id: usize, start: i64, end: i64, reverse: bool) -> Vec<u8> {
        retrieve_quality(&self.qualities[read_id], start, end, reverse)
    }
}

/// `flt_quals`: whether the summed quality (minus `offset`) reaches the mean
/// `cutoff` before the end of the read.
pub fn quality_passes(quality: &[u8], offset: u8, cutoff: i64) -> bool {
    let minimum = quality.len() as i64 * cutoff;
    let mut total = 0i64;
    for &value in quality {
        if total >= minimum {
            break;
        }
        total += i64::from(value.wrapping_sub(offset));
    }
    total >= minimum
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
    fn keeps_reads_in_input_order() {
        let path = write_temporary_fasta(
            "rshifiasm_read_store_test.fa",
            b">first\nACGT\n>second\nACGTACGT\n",
        );

        let mut options = AssemblyOptions::default();
        options.input_file_names = vec![path];
        options.input_file_count = 1;
        options.read_length_cutoff = -1;

        let store = ReadStore::load_from_files(&options).unwrap();
        assert_eq!(store.total_reads, 2);
        assert_eq!(store.total_bases, 12);
        assert_eq!(store.read_name(0), "first");
        assert_eq!(store.read_length(1), 8);
        assert_eq!(store.read_sequence(1), b"ACGTACGT");
    }

    #[test]
    fn trims_adapters_from_both_ends() {
        let path = write_temporary_fasta("rshifiasm_adapter_test.fa", b">read\nXXACGTXX\n");

        let mut options = AssemblyOptions::default();
        options.input_file_names = vec![path];
        options.input_file_count = 1;
        options.adapter_length = 2;
        options.read_length_cutoff = -1;

        let store = ReadStore::load_from_files(&options).unwrap();
        assert_eq!(store.total_reads, 1);
        assert_eq!(store.read_length(0), 4);
        assert_eq!(store.read_sequence(0), b"ACGT");
    }

    #[test]
    fn quality_levels_pack_and_expand_like_the_reference() {
        // Phred 20 ("5") is level 1, Phred 40 ("I") is level 3.
        let packed = compress_quality(b"5I5!!", 33);
        assert_eq!(packed, vec![0b01_11_01_00, 0b00_00_00_00]);
        assert_eq!(retrieve_quality(&packed, 0, 5, false), vec![1, 3, 1, 0, 0]);
        assert_eq!(retrieve_quality(&packed, 0, 5, true), vec![0, 0, 1, 3, 1]);
        // Values past the end decode to zero.
        assert_eq!(retrieve_quality(&compress_quality(b"!!!!", 33), 0, 6, false), vec![0; 6]);
    }

    #[test]
    fn the_mean_quality_filter_follows_flt_quals() {
        assert!(quality_passes(b"IIII", 33, 10));
        assert!(!quality_passes(b"!!!!", 33, 10));
        // The scan stops as soon as the threshold is reached, so a low tail
        // does not matter.
        assert!(quality_passes(b"IIII!!!!", 33, 5));
    }
}

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
                store.push_read(&record.name, trimmed_sequence);
            }
        }

        Ok(store)
    }

    /// Appends one read, the equivalent of `ha_insert_read_len` followed by
    /// the packed sequence insertion.
    pub fn push_read(&mut self, name: &str, sequence: &[u8]) {
        self.read_lengths.push(sequence.len() as u64);
        self.read_names.push(name.to_string());
        self.sequences.push(PackedSequence::from_ascii(sequence));
        self.trio_labels.push(TRIO_LABEL_UNKNOWN);
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
}

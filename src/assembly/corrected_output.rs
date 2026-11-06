//! Optional output artefacts, the counterparts of `Output_corrected_reads`,
//! `Output_corrected_fastq` and `Output_PAF` in `Assembly.cpp`.

use std::fs::File;
use std::io::{self, BufWriter, Write};

use crate::command_line::options::{
    AssemblyOptions, FLAG_WRITE_ERROR_CORRECTED_READS, FLAG_WRITE_PAF,
};
use crate::correction::read_writeback::{write_corrected_fastq, write_corrected_reads};
use crate::overlap::paf::{PafHit, write_paf};
use crate::reads::ReadStore;

/// `Output_corrected_reads`: writes `<prefix>.ec.fa`.
pub fn output_corrected_reads(prefix: &str, reads: &ReadStore) -> io::Result<()> {
    let path = format!("{prefix}.ec.fa");
    let mut writer = BufWriter::new(File::create(&path)?);
    write_corrected_reads(&mut writer, reads)?;
    writer.flush()?;
    Ok(())
}

/// `Output_corrected_fastq`: writes `<prefix>.ec.fq`.
pub fn output_corrected_fastq(prefix: &str, reads: &ReadStore) -> io::Result<()> {
    let path = format!("{prefix}.ec.fq");
    let mut writer = BufWriter::new(File::create(&path)?);
    write_corrected_fastq(&mut writer, reads)?;
    writer.flush()?;
    Ok(())
}

/// `Output_PAF`: writes `<prefix>.ovlp.paf`.
pub fn output_paf<T: AsRef<[PafHit]>>(
    prefix: &str,
    reads: &ReadStore,
    paf: &[T],
) -> io::Result<()> {
    eprintln!("Writing PAF to disk ...... ");
    let path = format!("{prefix}.ovlp.paf");
    let mut writer = BufWriter::new(File::create(&path)?);
    write_paf(&mut writer, reads, paf)?;
    writer.flush()?;
    eprintln!("PAF has been written.");
    Ok(())
}

/// Writes every artefact requested by the command-line flags.
///
/// Returns the first I/O error, if any; the caller prints it with the same
/// text as the reference and exits with a non-zero status.
pub fn write_requested_outputs<T: AsRef<[PafHit]>>(
    options: &AssemblyOptions,
    reads: &ReadStore,
    paf: &[T],
    single_cell_quality: bool,
) -> io::Result<()> {
    if options.flags & FLAG_WRITE_ERROR_CORRECTED_READS != 0 {
        if single_cell_quality {
            output_corrected_fastq(&options.output_file_name, reads)?;
        } else {
            output_corrected_reads(&options.output_file_name, reads)?;
        }
    }
    if options.flags & FLAG_WRITE_PAF != 0 {
        output_paf(&options.output_file_name, reads, paf)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn store() -> ReadStore {
        let mut reads = ReadStore::default();
        reads.push_read("read_a", b"ACGT");
        reads
    }

    #[test]
    fn requested_outputs_are_written_to_the_prefix() {
        let directory = std::env::temp_dir().join("rshifiasm_output_test");
        std::fs::create_dir_all(&directory).unwrap();
        let prefix = directory.join("sample").to_string_lossy().to_string();

        let mut options = AssemblyOptions::default();
        options.output_file_name = prefix.clone();
        options.flags |= FLAG_WRITE_ERROR_CORRECTED_READS | FLAG_WRITE_PAF;

        let reads = store();
        let paf = vec![vec![PafHit {
            query_identifier: 0,
            query_start: 0,
            query_end: 4,
            target_identifier: 0,
            target_start: 0,
            target_end: 4,
            matching_length: 1,
            block_length: 4,
            reverse: false,
            without_large_indel: false,
            exact_region: false,
            deleted: false,
        }]];
        write_requested_outputs(&options, &reads, &paf, false).unwrap();

        assert_eq!(std::fs::read_to_string(format!("{prefix}.ec.fa")).unwrap(), ">read_a\nACGT\n");
        assert_eq!(
            std::fs::read_to_string(format!("{prefix}.ovlp.paf")).unwrap(),
            "read_a\t4\t0\t4\t+\tread_a\t4\t0\t4\t1\t4\t255\n"
        );
        std::fs::remove_file(format!("{prefix}.ec.fa")).unwrap();
        std::fs::remove_file(format!("{prefix}.ovlp.paf")).unwrap();
    }
}

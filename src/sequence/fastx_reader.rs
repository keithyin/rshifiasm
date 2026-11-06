//! Streaming FASTA/FASTQ reader with transparent gzip support.
//!
//! This is the counterpart of `kseq.h` as used by the original implementation:
//! the header is split at the first whitespace, sequences and qualities may
//! span several lines, and the reader never loads a whole file into memory.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};

/// One FASTA/FASTQ record, the equivalent of `kseq_t`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SequenceRecord {
    /// Text before the first whitespace of the header.
    pub name: String,
    /// Text after the first whitespace of the header.
    pub comment: String,
    /// Sequence bases in upper case as they appear in the file.
    pub sequence: Vec<u8>,
    /// Quality characters, empty for FASTA input.
    pub quality: Vec<u8>,
}

/// Streaming reader over one FASTA/FASTQ file, compressed or not.
pub struct FastxReader {
    reader: Box<dyn BufRead>,
    pending_line: Option<String>,
    line_buffer: String,
    finished: bool,
}

impl FastxReader {
    /// Opens `file_name`, detecting the gzip magic bytes if present.
    pub fn open(file_name: &str) -> io::Result<Self> {
        let mut raw_file = File::open(file_name)?;
        let mut magic = [0u8; 2];
        let read_count = raw_file.read(&mut magic)?;
        // The magic bytes have been consumed by the probe above; rewind so that
        // both the decompressor and the plain reader see the whole file.
        raw_file.seek(SeekFrom::Start(0))?;

        let reader: Box<dyn BufRead> = if read_count == 2 && magic == [0x1f, 0x8b] {
            let decoded = flate2::read::MultiGzDecoder::new(raw_file);
            Box::new(BufReader::with_capacity(1 << 20, decoded))
        } else {
            Box::new(BufReader::with_capacity(1 << 20, raw_file))
        };

        Ok(Self {
            reader,
            pending_line: None,
            line_buffer: String::new(),
            finished: false,
        })
    }

    fn read_line(&mut self) -> io::Result<Option<String>> {
        if let Some(line) = self.pending_line.take() {
            return Ok(Some(line));
        }
        self.line_buffer.clear();
        let bytes_read = self.reader.read_line(&mut self.line_buffer)?;
        if bytes_read == 0 {
            return Ok(None);
        }
        Ok(Some(std::mem::take(&mut self.line_buffer)))
    }

    fn push_back(&mut self, line: String) {
        self.pending_line = Some(line);
    }

    /// Reads the next record, or `None` at end of file.
    pub fn next_record(&mut self) -> io::Result<Option<SequenceRecord>> {
        if self.finished {
            return Ok(None);
        }

        let header_line = loop {
            match self.read_line()? {
                None => {
                    self.finished = true;
                    return Ok(None);
                }
                Some(line) if line.trim().is_empty() => continue,
                Some(line) => break line,
            }
        };

        let header = header_line.trim_end();
        let (name, comment) = split_header(&header[1..]);
        let is_fastq = header.starts_with('@');

        let mut sequence = Vec::new();
        loop {
            match self.read_line()? {
                None => break,
                Some(line) => {
                    let trimmed = line.trim_end();
                    if trimmed.starts_with('>') || trimmed.starts_with('@') {
                        self.push_back(line);
                        break;
                    }
                    if is_fastq && trimmed.starts_with('+') {
                        break;
                    }
                    sequence.extend(trimmed.bytes());
                }
            }
        }

        let mut quality = Vec::new();
        if is_fastq {
            loop {
                if quality.len() >= sequence.len() {
                    break;
                }
                match self.read_line()? {
                    None => break,
                    Some(line) => {
                        let trimmed = line.trim_end();
                        quality.extend(trimmed.bytes());
                    }
                }
            }
            if quality.len() != sequence.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "sequence and quality length mismatch for read {}: {} vs {}",
                        name,
                        sequence.len(),
                        quality.len()
                    ),
                ));
            }
        }

        Ok(Some(SequenceRecord {
            name,
            comment,
            sequence,
            quality,
        }))
    }
}

/// Splits `name comment` at the first whitespace, like `kseq.h`.
fn split_header(header: &str) -> (String, String) {
    let trimmed = header.trim();
    match trimmed.find(char::is_whitespace) {
        Some(position) => (
            trimmed[..position].to_string(),
            trimmed[position..].trim_start().to_string(),
        ),
        None => (trimmed.to_string(), String::new()),
    }
}

/// Reads only the first record of `file_name`, used by the option checker the
/// same way `get_queries` probes the input files.
pub fn read_first_record(file_name: &str) -> Result<SequenceRecord, String> {
    let mut reader = FastxReader::open(file_name).map_err(|error| error.to_string())?;
    match reader.next_record() {
        Ok(Some(record)) => Ok(record),
        Ok(None) => Ok(SequenceRecord::default()),
        Err(error) => Err(error.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temporary_file(name: &str, content: &[u8]) -> String {
        let path = std::env::temp_dir().join(name);
        let mut file = File::create(&path).unwrap();
        file.write_all(content).unwrap();
        path.to_string_lossy().to_string()
    }

    #[test]
    fn reads_multiline_fasta_records() {
        let path = write_temporary_file(
            "rshifiasm_fasta_test.fa",
            b">read1 comment here\nACGT\nACGT\n>read2\nTTTT\n",
        );
        let mut reader = FastxReader::open(&path).unwrap();

        let first = reader.next_record().unwrap().unwrap();
        assert_eq!(first.name, "read1");
        assert_eq!(first.comment, "comment here");
        assert_eq!(first.sequence, b"ACGTACGT");
        assert!(first.quality.is_empty());

        let second = reader.next_record().unwrap().unwrap();
        assert_eq!(second.name, "read2");
        assert_eq!(second.sequence, b"TTTT");

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn reads_fastq_records() {
        let path = write_temporary_file(
            "rshifiasm_fastq_test.fq",
            b"@read1\nACGT\n+\nIIII\n@read2\nTT\n+\n##\n",
        );
        let mut reader = FastxReader::open(&path).unwrap();

        let first = reader.next_record().unwrap().unwrap();
        assert_eq!(first.quality, b"IIII");

        let second = reader.next_record().unwrap().unwrap();
        assert_eq!(second.sequence, b"TT");
        assert_eq!(second.quality, b"##");
    }

    #[test]
    fn rejects_fastq_records_with_mismatched_quality() {
        let path = write_temporary_file("rshifiasm_broken_test.fq", b"@read1\nACGT\n+\nII\n");
        let mut reader = FastxReader::open(&path).unwrap();
        assert!(reader.next_record().is_err());
    }
}

//! Binary stage dumps in the layout of the reference probe.
//!
//! The reference replay (`/tmp/clean_probe.cpp`) calls `write_debug_graph`
//! after every major cleaning step and writes
//! `<tag>.all.debug.{source,reverse,coverage_cut,ruIndex,asg_t}.bin`. Setting
//! `RSHIFIASM_PIPELINE_DUMP_DIR` makes the Rust pipeline write the same files
//! at the same boundaries, so two runs can be compared field by field.
//!
//! The layouts mirror the reference structs byte for byte; the loaders in
//! `tests/*_golden.rs` are the ground truth.

use std::io;
use std::path::{Path, PathBuf};

use crate::assembly::asg::{Asg, AsgArc, AsgSeq};
use crate::assembly::hit::{CoverageCut, ReadToUnitig};
use crate::correction::ec_driver::PafList;
use crate::overlap::paf::PafHit;
use std::fs;

/// The directory requested by `RSHIFIASM_PIPELINE_DUMP_DIR`, if any.
pub fn dump_directory() -> Option<PathBuf> {
    std::env::var_os("RSHIFIASM_PIPELINE_DUMP_DIR").map(PathBuf::from)
}

/// One per-read overlap list as the two representations the pipeline holds.
enum SourceView<'a> {
    Lists(&'a [PafList]),
    Raw(&'a [Vec<PafHit>], &'a [(bool, u8)]),
}

impl SourceView<'_> {
    fn write(&self, path: &Path) -> io::Result<()> {
        let mut out = Vec::new();
        match self {
            SourceView::Lists(lists) => {
                out.extend_from_slice(&(lists.len() as u64).to_le_bytes());
                for list in lists.iter() {
                    out.push(list.is_fully_corrected as u8);
                    out.push(list.is_abnormal);
                    write_hits(&mut out, &list.hits);
                }
            }
            SourceView::Raw(lists, flags) => {
                out.extend_from_slice(&(lists.len() as u64).to_le_bytes());
                for (list, &(is_fully_corrected, is_abnormal)) in lists.iter().zip(flags.iter()) {
                    out.push(is_fully_corrected as u8);
                    out.push(is_abnormal);
                    write_hits(&mut out, list);
                }
            }
        }
        std::fs::write(path, out)
    }
}

/// `write_ma`/`ma_hit_t`: one 32-byte record per overlap.
fn write_hits(out: &mut Vec<u8>, hits: &[PafHit]) {
    out.extend_from_slice(&(hits.len() as u32).to_le_bytes());
    for hit in hits {
        let qns = ((hit.query_identifier as u64) << 32) | u64::from(hit.query_start);
        out.extend_from_slice(&qns.to_le_bytes());
        out.extend_from_slice(&hit.query_end.to_le_bytes());
        out.extend_from_slice(&hit.target_identifier.to_le_bytes());
        out.extend_from_slice(&hit.target_start.to_le_bytes());
        out.extend_from_slice(&hit.target_end.to_le_bytes());
        // `ml` in bit 30 of the field word, `rev` in bit 31.
        let field = ((hit.reverse as u32) << 31) | ((hit.matching_length & 1) << 30);
        out.extend_from_slice(&field.to_le_bytes());
        // `bl` in the low 31 bits, `del` in bit 31.
        let block = ((hit.deleted as u32) << 31) | (hit.block_length & 0x7fff_ffff);
        out.extend_from_slice(&block.to_le_bytes());
        out.push(hit.exact_region as u8);
        out.push(hit.without_large_indel as u8);
        out.extend_from_slice(&[0u8; 6]);
    }
}

/// `write_coverage_cut`.
fn write_coverage_cut(path: &Path, cuts: &[CoverageCut]) -> io::Result<()> {
    let mut out = Vec::with_capacity(8 + cuts.len() * 12);
    out.extend_from_slice(&(cuts.len() as u64).to_le_bytes());
    for cut in cuts {
        let start_word = ((cut.deleted as u32) << 31) | (cut.start & 0x7fff_ffff);
        out.extend_from_slice(&start_word.to_le_bytes());
        out.extend_from_slice(&cut.end.to_le_bytes());
        out.push(cut.class);
        out.extend_from_slice(&[0u8; 3]);
    }
    std::fs::write(path, out)
}

/// `write_ruIndex`.
///
/// The stored coverage thresholds are not part of this pipeline stage; they
/// are written as `(0, u32::MAX)` and only compared by the stages that set
/// them.
fn write_read_to_unitig(path: &Path, map: &ReadToUnitig, trio_labels: &[u8]) -> io::Result<()> {
    let mut out = Vec::with_capacity(12 + map.len() * 5);
    out.extend_from_slice(&(map.len() as u64).to_le_bytes());
    for &value in map.entries() {
        out.extend_from_slice(&value.to_le_bytes());
    }
    let mut labels = vec![0u8; map.len()];
    let available = trio_labels.len().min(map.len());
    labels[..available].copy_from_slice(&trio_labels[..available]);
    out.extend_from_slice(&labels);
    out.extend_from_slice(&0u32.to_le_bytes());
    out.extend_from_slice(&u32::MAX.to_le_bytes());
    std::fs::write(path, out)
}

/// `write_asg_t`.
fn write_asg(path: &Path, graph: &Asg) -> io::Result<()> {
    let sequence_count = graph.sequences.len();
    let vertex_count = sequence_count * 2;
    let mut out = Vec::new();
    out.extend_from_slice(&(graph.arcs.len() as u32).to_le_bytes());
    // `m_arc`: the reference writes its capacity; the loader ignores it.
    out.extend_from_slice(&(graph.arcs.len() as u32).to_le_bytes());
    out.extend_from_slice(&(graph.is_sorted as u32).to_le_bytes());
    out.extend_from_slice(&(sequence_count as u32).to_le_bytes());
    out.extend_from_slice(&(sequence_count as u32).to_le_bytes());
    out.extend_from_slice(&(graph.is_symmetric as u32).to_le_bytes());
    out.extend_from_slice(&graph.read_sequence_count.to_le_bytes());

    let mut visited = vec![0u8; vertex_count];
    let available = graph.sequence_visited.len().min(vertex_count);
    visited[..available].copy_from_slice(&graph.sequence_visited[..available]);
    out.extend_from_slice(&visited);

    assert!(
        graph.index.len() >= vertex_count,
        "the arc index must cover every vertex"
    );
    for &entry in &graph.index[..vertex_count] {
        out.extend_from_slice(&entry.to_le_bytes());
    }

    for arc in &graph.arcs {
        out.extend_from_slice(&arc.ul.to_le_bytes());
        out.extend_from_slice(&arc.v.to_le_bytes());
        let overlap_word = ((arc.deleted as u32) << 31) | (arc.overlap & 0x7fff_ffff);
        out.extend_from_slice(&overlap_word.to_le_bytes());
        let field_word = (arc.occurrence & 0x3fff)
            | ((arc.strong as u16) << 14)
            | ((arc.without_large_indel as u16) << 15);
        out.extend_from_slice(&field_word.to_le_bytes());
        out.push(arc.exact_region);
        out.extend_from_slice(&[0u8; 5]);
    }

    for sequence in &graph.sequences {
        let length_word = ((sequence.deleted as u32) << 31) | (sequence.length & 0x7fff_ffff);
        out.extend_from_slice(&length_word.to_le_bytes());
        out.push(sequence.coverage_class);
        out.extend_from_slice(&[0u8; 3]);
    }

    std::fs::write(path, out)
}

fn write_stage(
    directory: &Path,
    tag: &str,
    graph: Option<&Asg>,
    sources: SourceView<'_>,
    reverse_sources: SourceView<'_>,
    coverage_cut: Option<&[CoverageCut]>,
    read_to_unitig: &ReadToUnitig,
    trio_labels: &[u8],
) -> io::Result<()> {
    sources.write(&directory.join(format!("{}.all.debug.source.bin", tag)))?;
    reverse_sources.write(&directory.join(format!("{}.all.debug.reverse.bin", tag)))?;
    if let Some(cuts) = coverage_cut {
        write_coverage_cut(&directory.join(format!("{}.all.debug.coverage_cut.bin", tag)), cuts)?;
    }
    write_read_to_unitig(
        &directory.join(format!("{}.all.debug.ruIndex.bin", tag)),
        read_to_unitig,
        trio_labels,
    )?;
    if let Some(graph) = graph {
        write_asg(&directory.join(format!("{}.all.debug.asg_t.bin", tag)), graph)?;
    }
    Ok(())
}

/// Writes one stage after the overlap lists have been joined into
/// [`PafList`]s.
#[allow(clippy::too_many_arguments)]
pub fn dump_stage(
    directory: &Path,
    tag: &str,
    graph: Option<&Asg>,
    sources: &[PafList],
    reverse_sources: &[PafList],
    coverage_cut: Option<&[CoverageCut]>,
    read_to_unitig: &ReadToUnitig,
    trio_labels: &[u8],
) -> io::Result<()> {
    write_stage(
        directory,
        tag,
        graph,
        SourceView::Lists(sources),
        SourceView::Lists(reverse_sources),
        coverage_cut,
        read_to_unitig,
        trio_labels,
    )
}

/// Writes one stage while the pipeline still holds the raw hit vectors.
#[allow(clippy::too_many_arguments)]
pub fn dump_stage_raw(
    directory: &Path,
    tag: &str,
    graph: Option<&Asg>,
    sources: &[Vec<PafHit>],
    source_flags: &[(bool, u8)],
    reverse_sources: &[Vec<PafHit>],
    reverse_flags: &[(bool, u8)],
    coverage_cut: Option<&[CoverageCut]>,
    read_to_unitig: &ReadToUnitig,
    trio_labels: &[u8],
) -> io::Result<()> {
    write_stage(
        directory,
        tag,
        graph,
        SourceView::Raw(sources, source_flags),
        SourceView::Raw(reverse_sources, reverse_flags),
        coverage_cut,
        read_to_unitig,
        trio_labels,
    )
}

/// Ensures the dump directory exists; returns it when dumping is enabled.
pub fn prepare_directory() -> io::Result<Option<PathBuf>> {
    match dump_directory() {
        None => Ok(None),
        Some(directory) => {
            std::fs::create_dir_all(&directory)?;
            Ok(Some(directory))
        }
    }
}

/// Writes a stage only when dumping is enabled; errors are reported once to
/// stderr instead of aborting the assembly.
pub fn try_dump_stage(
    directory: Option<&Path>,
    tag: &str,
    graph: Option<&Asg>,
    sources: &[PafList],
    reverse_sources: &[PafList],
    coverage_cut: Option<&[CoverageCut]>,
    read_to_unitig: &ReadToUnitig,
    trio_labels: &[u8],
) {
    let Some(directory) = directory else {
        return;
    };
    if let Err(error) = dump_stage(
        directory,
        tag,
        graph,
        sources,
        reverse_sources,
        coverage_cut,
        read_to_unitig,
        trio_labels,
    ) {
        eprintln!("[WARN] failed to dump {}: {}", tag, error);
    }
}

/// Like [`try_dump_stage`], for the raw overlap representation.
#[allow(clippy::too_many_arguments)]
pub fn try_dump_stage_raw(
    directory: Option<&Path>,
    tag: &str,
    graph: Option<&Asg>,
    sources: &[Vec<PafHit>],
    source_flags: &[(bool, u8)],
    reverse_sources: &[Vec<PafHit>],
    reverse_flags: &[(bool, u8)],
    coverage_cut: Option<&[CoverageCut]>,
    read_to_unitig: &ReadToUnitig,
    trio_labels: &[u8],
) {
    let Some(directory) = directory else {
        return;
    };
    if let Err(error) = dump_stage_raw(
        directory,
        tag,
        graph,
        sources,
        source_flags,
        reverse_sources,
        reverse_flags,
        coverage_cut,
        read_to_unitig,
        trio_labels,
    ) {
        eprintln!("[WARN] failed to dump {}: {}", tag, error);
    }
}


// ---------------------------------------------------------------------------
// Readers for the dump format, shared with the stage golden test. They mirror
// `load_asg_t`, `load_ma_hit_ts`, `load_coverage_cut` and `load_ruIndex`.
// ---------------------------------------------------------------------------

fn read_u32(data: &[u8], position: &mut usize) -> u32 {
    let value = u32::from_le_bytes(data[*position..*position + 4].try_into().unwrap());
    *position += 4;
    value
}

fn read_u16(data: &[u8], position: &mut usize) -> u16 {
    let value = u16::from_le_bytes(data[*position..*position + 2].try_into().unwrap());
    *position += 2;
    value
}

fn read_u64(data: &[u8], position: &mut usize) -> u64 {
    let value = u64::from_le_bytes(data[*position..*position + 8].try_into().unwrap());
    *position += 8;
    value
}

fn read_u8(data: &[u8], position: &mut usize) -> u8 {
    let value = data[*position];
    *position += 1;
    value
}

/// `load_asg_t`.
pub fn load_asg(path: &Path) -> Asg {
    let data = fs::read(path).expect("the asg dump must be readable");
    let mut position = 0usize;
    let count_arc = read_u32(&data, &mut position) as usize;
    let _capacity_arc = read_u32(&data, &mut position);
    let is_sorted = read_u32(&data, &mut position) != 0;
    let count_sequence = read_u32(&data, &mut position) as usize;
    let _capacity_sequence = read_u32(&data, &mut position);
    let is_symmetric = read_u32(&data, &mut position) != 0;
    let read_sequence_count = read_u32(&data, &mut position);

    let mut sequence_visited = vec![0u8; count_sequence * 2];
    sequence_visited.copy_from_slice(&data[position..position + count_sequence * 2]);
    position += count_sequence * 2;

    let mut index = vec![0u64; count_sequence * 2];
    for entry in index.iter_mut() {
        *entry = read_u64(&data, &mut position);
    }

    let mut arcs = Vec::with_capacity(count_arc);
    for _ in 0..count_arc {
        let ul = read_u64(&data, &mut position);
        let v = read_u32(&data, &mut position);
        let overlap_word = read_u32(&data, &mut position);
        let field_word = read_u16(&data, &mut position);
        let exact_region = read_u8(&data, &mut position);
        position += 5;
        arcs.push(AsgArc {
            ul,
            v,
            overlap: overlap_word & 0x7fff_ffff,
            deleted: overlap_word >> 31 != 0,
            occurrence: field_word & 0x3fff,
            strong: (field_word >> 14) & 1 != 0,
            without_large_indel: (field_word >> 15) & 1 != 0,
            exact_region,
        });
    }

    let mut sequences = Vec::with_capacity(count_sequence);
    for _ in 0..count_sequence {
        let length_word = read_u32(&data, &mut position);
        let coverage_class = read_u8(&data, &mut position);
        position += 3;
        sequences.push(AsgSeq {
            length: length_word & 0x7fff_ffff,
            deleted: length_word >> 31 != 0,
            coverage_class,
        });
    }

    assert_eq!(position, data.len(), "the asg dump has trailing bytes");
    Asg {
        arcs,
        index,
        is_sorted,
        is_symmetric,
        sequences,
        read_sequence_count,
        sequence_visited,
    }
}

/// `load_ma_hit_ts`.
pub fn load_sources(path: &Path) -> Vec<PafList> {
    let data = fs::read(path).expect("the overlap dump must be readable");
    let mut position = 0usize;
    let read_count = read_u64(&data, &mut position) as usize;
    let mut lists = Vec::with_capacity(read_count);

    for _ in 0..read_count {
        let is_fully_corrected = read_u8(&data, &mut position) != 0;
        let is_abnormal = read_u8(&data, &mut position);
        let length = read_u32(&data, &mut position) as usize;
        let mut hits = Vec::with_capacity(length);
        for _ in 0..length {
            let qns = read_u64(&data, &mut position);
            let query_end = read_u32(&data, &mut position);
            let target_identifier = read_u32(&data, &mut position);
            let target_start = read_u32(&data, &mut position);
            let target_end = read_u32(&data, &mut position);
            let field_word = read_u32(&data, &mut position);
            let block_word = read_u32(&data, &mut position);
            let exact_region = read_u8(&data, &mut position);
            let without_large_indel = read_u8(&data, &mut position) != 0;
            position += 6;

            hits.push(PafHit {
                query_identifier: (qns >> 32) as u32,
                query_start: qns as u32,
                query_end,
                target_identifier,
                target_start,
                target_end,
                matching_length: (field_word >> 30) & 1,
                block_length: block_word & 0x7fff_ffff,
                reverse: field_word >> 31 != 0,
                without_large_indel,
                exact_region: exact_region != 0,
                exact_region_byte: exact_region,
                deleted: block_word >> 31 != 0,
            });
        }
        lists.push(PafList {
            hits,
            is_fully_corrected,
            is_abnormal,
        });
    }

    assert_eq!(position, data.len(), "the overlap dump has trailing bytes");
    lists
}

/// `load_coverage_cut`.
pub fn load_coverage_cut(path: &Path) -> Vec<CoverageCut> {
    let data = fs::read(path).expect("the coverage dump must be readable");
    let mut position = 0usize;
    let read_count = read_u64(&data, &mut position) as usize;
    let mut cuts = Vec::with_capacity(read_count);
    for _ in 0..read_count {
        let start_word = read_u32(&data, &mut position);
        let end = read_u32(&data, &mut position);
        let coverage_class = read_u8(&data, &mut position);
        position += 3;
        cuts.push(CoverageCut {
            start: start_word & 0x7fff_ffff,
            end,
            deleted: start_word >> 31 != 0,
            class: coverage_class,
        });
    }
    assert_eq!(position, data.len(), "the coverage dump has trailing bytes");
    cuts
}

/// `load_ruIndex`; the trio labels and the stored coverage thresholds travel
/// in the same file.
pub fn load_read_to_unitig(path: &Path) -> (ReadToUnitig, Vec<u8>, u32, u32) {
    let data = fs::read(path).expect("the read map dump must be readable");
    let mut position = 0usize;
    let length = read_u64(&data, &mut position) as usize;

    let mut values = Vec::with_capacity(length);
    for _ in 0..length {
        values.push(read_u32(&data, &mut position));
    }
    let mut trio_labels = vec![0u8; length];
    trio_labels.copy_from_slice(&data[position..position + length]);
    position += length;
    let hom_global_coverage_set = read_u32(&data, &mut position);
    let hom_global_coverage = read_u32(&data, &mut position);
    assert_eq!(position, data.len(), "the read map dump has trailing bytes");

    let mut read_to_unitig = ReadToUnitig::new(length);
    for (read, &value) in values.iter().enumerate() {
        read_to_unitig.set(read as u32, value & 0x7fff_ffff, value >> 31 != 0, None);
    }
    (
        read_to_unitig,
        trio_labels,
        hom_global_coverage_set,
        hom_global_coverage,
    )
}

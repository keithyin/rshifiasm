//! `gen_hc_fast_cigar0`: the per-overlap driver of the new HiFi correction
//! path.
//!
//! Chaining hands the corrector a list of candidate seeds grouped by target
//! read; `region.shared_seed` is the index of the first seed of the group the
//! overlap was built from. The driver aligns the whole overlap block by block,
//! retries every window that came out unaligned, drops the stale unaligned
//! windows and returns the remaining error count.
//!
//! The retry pass is where the dataflow matters: `region` is still the
//! original overlap, so its window list holds the tiling with no cached
//! cigars, and [`rechain_alignment_hc`] resolves those windows by backtracing
//! the target read. The first pass writes into `auxiliary`, not into `region`.

use crate::container::ksort::radix_sort_by_key_32;
use crate::correction::hc_align::BitExtz;
use crate::correction::ovlp_base::hc_ovlp_base_direct;
use crate::correction::window::WindowList;
use crate::correction::window_chain::{is_unaligned_window, rechain_alignment_hc};
use crate::overlap::region::{CandidatesList, OverlapRegion};
use crate::correction::target::TargetSource;

/// `gen_hc_fast_cigar0`: aligns one overlap and retries its failed windows.
///
/// Returns the summed error of the final window list, or `None` when the
/// candidate list holds no seed for the overlap (the reference returns `0`
/// without touching its output).
///
/// `read_identifier` of the reference is debug-only and is not carried over.
#[allow(clippy::too_many_arguments)]
pub fn generate_fast_cigar_direct(
    region: &mut OverlapRegion,
    candidates: &mut CandidatesList,
    window_length: i64,
    targets: &TargetSource,
    query: &[u8],
    extz: &mut BitExtz,
    auxiliary: &mut OverlapRegion,
    error_rate: f64,
    query_total_length: i64,
    h_kmer_hit: i64,
) -> Option<i64> {
    let chain_index = i64::from(region.shared_seed);
    if chain_index < 0 || chain_index >= candidates.length {
        return None;
    }
    let chain_start = chain_index as usize;
    let first_read_identifier = candidates.list[chain_start].read_identifier;
    let mut scan = chain_start;
    while scan < candidates.length as usize
        && candidates.list[scan].read_identifier == first_read_identifier
    {
        scan += 1;
    }
    let chain_count = scan as i64 - chain_index;
    if chain_count <= 0 {
        return None;
    }

    let target_total_length = targets.length(region.y_identifier as usize);

    auxiliary.window_cigar.clear();
    auxiliary.y_identifier = region.y_identifier;
    auxiliary.y_strand = region.y_strand;
    auxiliary.x_position_start = region.x_position_start;
    auxiliary.x_position_end = region.x_position_end;
    auxiliary.y_position_start = region.y_position_start;
    auxiliary.y_position_end = region.y_position_end;

    hc_ovlp_base_direct(
        region,
        &candidates.list[chain_start..chain_start + chain_count as usize],
        chain_count,
        window_length,
        targets,
        query,
        extz,
        &mut auxiliary.window_cigar,
        error_rate,
        query_total_length,
        target_total_length,
        -1,
    );

    let initial_window_count = auxiliary.window_cigar.windows.len();
    for index in 0..initial_window_count {
        if !is_unaligned_window(&auxiliary.window_cigar.windows[index]) {
            continue;
        }
        // The retry may overwrite the candidate entries the first pass used;
        // `rechain_alignment_hc` restores the list length on return.
        rechain_alignment_hc(
            region,
            candidates,
            &mut auxiliary.window_cigar,
            index,
            window_length,
            targets,
            query,
            extz,
            error_rate,
            query_total_length,
            target_total_length,
            h_kmer_hit,
        );
    }
    // The retry appended windows at the end. Drop the initial windows that are
    // still unaligned, keep everything else, and restore the query order.
    if auxiliary.window_cigar.windows.len() > initial_window_count {
        let mut write = 0usize;
        for read in 0..auxiliary.window_cigar.windows.len() {
            if read < initial_window_count
                && is_unaligned_window(&auxiliary.window_cigar.windows[read])
            {
                continue;
            }
            auxiliary.window_cigar.windows[write] = auxiliary.window_cigar.windows[read];
            write += 1;
        }
        auxiliary.window_cigar.windows.truncate(write);
        // `radix_sort_window_list_xs_srt` is a stable 32-bit radix sort on the
        // signed `x_start`; window starts are never negative.
        radix_sort_by_key_32(&mut auxiliary.window_cigar.windows, &|window: &WindowList| {
            window.x_start as u32
        });
    }

    region.update_from_windows(
        &auxiliary.window_cigar,
        query_total_length,
        target_total_length,
    );

    let mut total_error = 0i64;
    for window in &region.window_cigar.windows {
        if is_unaligned_window(window) {
            let query_extent = i64::from(window.x_end) + 1 - i64::from(window.x_start);
            let target_extent = i64::from(window.y_end) + 1 - i64::from(window.y_start);
            total_error += query_extent.max(target_extent);
        } else {
            total_error += i64::from(window.error);
        }
    }
    Some(total_error)
}

/// Sentinel `readID` that `return_t_chain` writes into every hit its refined
/// chain drops.
pub const DROPPED_HIT: u32 = 0x7fff_ffff;

/// `return_t_chain`: refines the first candidate group of the overlap.
///
/// The group is the run of hits starting at `region.shared_seed` that shares
/// the first hit's target read. `chain_refine` reorders the hits it keeps to
/// the front of the run and the dropped ones are marked [`DROPPED_HIT`], which
/// keeps them out of the group that `gen_hc_fast_cigar0` scans next.
///
/// Returns the number of hits in the refined chain.
pub fn return_target_chain(region: &OverlapRegion, candidates: &mut CandidatesList) -> i64 {
    let start = i64::from(region.shared_seed);
    if start < 0 || start >= candidates.length {
        return 0;
    }
    let start = start as usize;
    let read_identifier = candidates.list[start].read_identifier;
    let mut scan = start;
    while scan < candidates.length as usize
        && candidates.list[scan].read_identifier == read_identifier
        && candidates.list[scan].read_identifier != DROPPED_HIT
    {
        scan += 1;
    }
    let chain_size = scan as i64 - start as i64;
    let refined = crate::overlap::chaining::chain_refine(
        &mut candidates.list[start..scan],
        chain_size,
        &mut candidates.chain_data,
        50,
        5000,
        512,
        16,
    );
    for hit in &mut candidates.list[start + refined as usize..scan] {
        hit.read_identifier = DROPPED_HIT;
    }
    refined
}

/// `gen_hc_fast_cigar`: refines the candidate chain, then aligns the overlap.
///
/// The reference ignores the alignment's return value and always reports
/// success. The second element is the `re` output, `None` when the alignment
/// returned early without touching it.
#[allow(clippy::too_many_arguments)]
pub fn generate_fast_cigar(
    region: &mut OverlapRegion,
    candidates: &mut CandidatesList,
    window_length: i64,
    targets: &TargetSource,
    query: &[u8],
    extz: &mut BitExtz,
    auxiliary: &mut OverlapRegion,
    error_rate: f64,
    query_total_length: i64,
    h_kmer_hit: i64,
) -> (i64, Option<i64>) {
    return_target_chain(region, candidates);
    let total_error = generate_fast_cigar_direct(
        region,
        candidates,
        window_length,
        targets,
        query,
        extz,
        auxiliary,
        error_rate,
        query_total_length,
        h_kmer_hit,
    );
    (1, total_error)
}

#[cfg(test)]
mod probe_goldens {
    use super::*;
    use crate::reads::read_store::ReadStore;
    use crate::correction::window::WindowCigar;
    use crate::overlap::region::KmerHit;

    struct Lcg(u64);

    impl Lcg {
        fn next_u32(&mut self) -> u32 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 33) as u32
        }

        fn sequence(&mut self, length: usize) -> Vec<u8> {
            const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
            (0..length)
                .map(|_| BASES[(self.next_u32() & 3) as usize])
                .collect()
        }
    }

    fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                other => *other,
            })
            .collect()
    }

    /// `print_windows`: the reference zeroes the extra boundaries of every
    /// aligned window before printing, and that zeroing persists.
    fn push_windows(lines: &mut Vec<String>, tag: &str, cigar: &mut WindowCigar) {
        lines.push(format!(
            "{} n={} cn={}",
            tag,
            cigar.windows.len(),
            cigar.cells.len()
        ));
        for index in 0..cigar.windows.len() {
            if cigar.windows[index].error != i16::MAX {
                cigar.windows[index].extra_begin = 0;
                cigar.windows[index].extra_end = 0;
            }
            let window = cigar.windows[index];
            let mut line = format!(
                " {}W {} x=[{},{}] y=[{},{}] eb={} ee={} err={} thre={} cidx={} clen={} :",
                tag,
                index,
                window.x_start,
                window.x_end,
                window.y_start,
                window.y_end,
                window.extra_begin,
                window.extra_end,
                window.error,
                window.error_threshold,
                window.cigar_start,
                window.cigar_length
            );
            for offset in 0..window.cigar_length {
                let cell = cigar.cells[(window.cigar_start + offset) as usize];
                line.push_str(&format!(
                    " {}{}",
                    cell & 0x3fff,
                    b"EMDI"[(cell >> 14) as usize] as char
                ));
            }
            lines.push(line);
        }
    }

    #[derive(Clone, Copy)]
    struct WindowSpec {
        x_start: i32,
        x_end: i32,
        y_start: i32,
        y_end: i32,
        error: i16,
        threshold: i16,
    }

    struct ReadScenario {
        name: &'static str,
        reads: Vec<Vec<u8>>,
        x_start: u32,
        x_end: u32,
        y_start: u32,
        y_end: u32,
        y_strand: u32,
        non_homopolymer_errors: u32,
        gaps: Vec<(u32, i32)>,
        windows: Vec<WindowSpec>,
        hits: Vec<(u32, u32)>,
        window_length: i64,
    }

    fn build_base(scenario: &ReadScenario, non_homopolymer_errors: u32) -> OverlapRegion {
        let mut region = OverlapRegion {
            x_identifier: 0,
            x_position_start: scenario.x_start,
            x_position_end: scenario.x_end,
            y_identifier: 1,
            y_position_start: scenario.y_start,
            y_position_end: scenario.y_end,
            y_strand: scenario.y_strand,
            overlap_length: scenario.x_end + 1 - scenario.x_start,
            non_homopolymer_errors,
            ..OverlapRegion::default()
        };
        for &(site, shift) in &scenario.gaps {
            region.fake_cigar.push(site, shift);
        }
        for window in &scenario.windows {
            region.window_cigar.append_window(
                window.x_start,
                window.x_end,
                window.y_start,
                window.y_end,
                window.error,
                0,
                0,
                window.threshold,
            );
        }
        region
    }

    #[allow(clippy::too_many_arguments)]
    fn run_fast_cigar(
        scenario: &ReadScenario,
        suffix: &str,
        non_homopolymer_errors: u32,
        chain_prefix: i64,
        h_kmer_hit: i64,
        lines: &mut Vec<String>,
    ) {
        let mut reads = ReadStore::default();
        for (index, sequence) in scenario.reads.iter().enumerate() {
            reads.push_read(&format!("read{}", index), sequence);
        }

        let mut base = build_base(scenario, non_homopolymer_errors);
        base.shared_seed = chain_prefix as i32;

        let mut extz = BitExtz::new();
        let mut auxiliary = OverlapRegion::default();

        let mut candidates = CandidatesList::default();
        for index in 0..chain_prefix {
            candidates.push(KmerHit {
                read_identifier: 2,
                reverse_strand: false,
                offset: (3000 + index) as u32,
                self_offset: (3000 + index) as u32,
                count: 1,
            });
        }
        for &(self_offset, offset) in &scenario.hits {
            candidates.push(KmerHit {
                read_identifier: 1,
                reverse_strand: false,
                offset,
                self_offset,
                count: 1,
            });
        }

        let query = &scenario.reads[0];
        let query_total_length = query.len() as i64;
        let result = generate_fast_cigar_direct(
            &mut base,
            &mut candidates,
            scenario.window_length,
            &TargetSource::Reads(&reads),
            query,
            &mut extz,
            &mut auxiliary,
            0.04,
            query_total_length,
            h_kmer_hit,
        );

        lines.push(format!(
            "FC name={}_{} ret={} re={} ql={} wl={} hk={} nhpe={} ch_idx={} length={}",
            scenario.name,
            suffix,
            i32::from(result.is_some()),
            result.unwrap_or(-1),
            query_total_length,
            scenario.window_length,
            h_kmer_hit,
            non_homopolymer_errors,
            chain_prefix,
            candidates.length
        ));
        lines.push(format!(
            "FCZ x=[{},{}] y=[{},{}] strand={} n={} cn={}",
            base.x_position_start,
            base.x_position_end,
            base.y_position_start,
            base.y_position_end,
            base.y_strand,
            base.window_cigar.windows.len(),
            base.window_cigar.cells.len()
        ));
        lines.push(format!(
            "FCA x=[{},{}] y=[{},{}] strand={}",
            auxiliary.x_position_start,
            auxiliary.x_position_end,
            auxiliary.y_position_start,
            auxiliary.y_position_end,
            auxiliary.y_strand
        ));
        push_windows(lines, "FCW", &mut auxiliary.window_cigar);
        push_windows(lines, "FCZ", &mut base.window_cigar);
    }

    // ------------------------------------------------ return_t_chain replay

    struct HitSpec {
        rid: u32,
        x: u32,
        y: u32,
        count: u32,
    }

    fn hit(rid: u32, x: u32, y: u32) -> HitSpec {
        HitSpec {
            rid,
            x,
            y,
            count: 1,
        }
    }

    fn push_hits(lines: &mut Vec<String>, tag: &str, hits: &[KmerHit]) {
        for (index, entry) in hits.iter().enumerate() {
            lines.push(format!(
                "{} i={} x={} y={} cnt={} rid={} st={}",
                tag,
                index,
                entry.self_offset,
                entry.offset,
                entry.count,
                entry.read_identifier,
                entry.strand_value()
            ));
        }
    }

    fn run_chain(name: &str, shared_seed: i64, hits: &[HitSpec], lines: &mut Vec<String>) {
        let region = OverlapRegion {
            shared_seed: shared_seed as i32,
            ..OverlapRegion::default()
        };
        let mut candidates = CandidatesList::default();
        for entry in hits {
            candidates.push(KmerHit {
                read_identifier: entry.rid,
                reverse_strand: false,
                offset: entry.y,
                self_offset: entry.x,
                count: entry.count,
            });
        }
        let refined = return_target_chain(&region, &mut candidates);
        lines.push(format!(
            "RT name={} seed={} length={} ret={}",
            name, shared_seed, candidates.length, refined
        ));
        push_hits(
            lines,
            "RTH",
            &candidates.list[..candidates.length as usize],
        );
    }

    struct WrapperSpec {
        name: &'static str,
        non_homopolymer_errors: u32,
        window_length: i64,
        hits: Vec<HitSpec>,
    }

    #[allow(clippy::too_many_arguments)]
    fn run_wrapper(
        sequences: &[Vec<u8>],
        x_start: u32,
        x_end: u32,
        y_start: u32,
        y_end: u32,
        gaps: &[(u32, i32)],
        windows: &[WindowSpec],
        spec: &WrapperSpec,
        lines: &mut Vec<String>,
    ) {
        let scenario = ReadScenario {
            name: spec.name,
            reads: sequences.to_vec(),
            x_start,
            x_end,
            y_start,
            y_end,
            y_strand: 0,
            non_homopolymer_errors: spec.non_homopolymer_errors,
            gaps: gaps.to_vec(),
            windows: windows.to_vec(),
            hits: Vec::new(),
            window_length: spec.window_length,
        };
        let mut reads = ReadStore::default();
        for (index, sequence) in scenario.reads.iter().enumerate() {
            reads.push_read(&format!("read{}", index), sequence);
        }
        let mut base = build_base(&scenario, spec.non_homopolymer_errors);
        base.shared_seed = 0;

        let mut extz = BitExtz::new();
        let mut auxiliary = OverlapRegion::default();
        let mut candidates = CandidatesList::default();
        for entry in &spec.hits {
            candidates.push(KmerHit {
                read_identifier: entry.rid,
                reverse_strand: false,
                offset: entry.y,
                self_offset: entry.x,
                count: entry.count,
            });
        }

        let query = &scenario.reads[0];
        let query_total_length = query.len() as i64;
        let (ret, error) = generate_fast_cigar(
            &mut base,
            &mut candidates,
            spec.window_length,
            &TargetSource::Reads(&reads),
            query,
            &mut extz,
            &mut auxiliary,
            0.04,
            query_total_length,
            10,
        );

        lines.push(format!(
            "GC name={} ret={} re={} ql={} wl={} nhpe={} length={}",
            spec.name,
            ret,
            error.unwrap_or(-1),
            query_total_length,
            spec.window_length,
            spec.non_homopolymer_errors,
            candidates.length
        ));
        lines.push(format!(
            "GCZ x=[{},{}] y=[{},{}] strand={} n={} cn={}",
            base.x_position_start,
            base.x_position_end,
            base.y_position_start,
            base.y_position_end,
            base.y_strand,
            base.window_cigar.windows.len(),
            base.window_cigar.cells.len()
        ));
        lines.push(format!(
            "GCA x=[{},{}] y=[{},{}] strand={}",
            auxiliary.x_position_start,
            auxiliary.x_position_end,
            auxiliary.y_position_start,
            auxiliary.y_position_end,
            auxiliary.y_strand
        ));
        push_hits(
            lines,
            "GCH",
            &candidates.list[..candidates.length as usize],
        );
        push_windows(lines, "GCW", &mut auxiliary.window_cigar);
        push_windows(lines, "GCU", &mut base.window_cigar);
    }

    fn read_scenarios() -> Vec<ReadScenario> {
        let mut generator = Lcg(987654321);
        let query = generator.sequence(2000);
        let prefix = generator.sequence(600);
        let suffix = generator.sequence(400);
        let tile_query = generator.sequence(1400);

        let exact_target: Vec<u8> = [prefix.clone(), query.clone(), suffix.clone()].concat();
        let mut substituted = exact_target.clone();
        let mut position = 0usize;
        while position < query.len() {
            let base = &mut substituted[prefix.len() + position];
            *base = if *base == b'A' { b'C' } else { b'A' };
            position += 97;
        }
        let reversed = reverse_complement(&exact_target);

        let mut exact_tile = Vec::new();
        for k in 0..4 {
            let start = k * 350;
            exact_tile.push(WindowSpec {
                x_start: start,
                x_end: start + 349,
                y_start: start,
                y_end: start + 349,
                error: 0,
                threshold: 0,
            });
        }
        let subs_windows = vec![
            WindowSpec {
                x_start: 200,
                x_end: 599,
                y_start: 800,
                y_end: 1199,
                error: 4,
                threshold: 0,
            },
            WindowSpec {
                x_start: 600,
                x_end: 999,
                y_start: 1200,
                y_end: 1599,
                error: 4,
                threshold: 0,
            },
            WindowSpec {
                x_start: 1000,
                x_end: 1399,
                y_start: 1600,
                y_end: 1999,
                error: 4,
                threshold: 0,
            },
            WindowSpec {
                x_start: 1400,
                x_end: 1799,
                y_start: 2000,
                y_end: 2399,
                error: 5,
                threshold: 0,
            },
        ];
        let reverse_windows = vec![
            WindowSpec {
                x_start: 200,
                x_end: 599,
                y_start: 800,
                y_end: 1199,
                error: 0,
                threshold: 0,
            },
            WindowSpec {
                x_start: 600,
                x_end: 999,
                y_start: 1200,
                y_end: 1599,
                error: 0,
                threshold: 0,
            },
            WindowSpec {
                x_start: 1000,
                x_end: 1399,
                y_start: 1600,
                y_end: 1999,
                error: 0,
                threshold: 0,
            },
            WindowSpec {
                x_start: 1400,
                x_end: 1799,
                y_start: 2000,
                y_end: 2399,
                error: 0,
                threshold: 0,
            },
        ];

        let long_query = generator.sequence(20000);
        let long_exact = long_query.clone();
        let mut long_substituted = long_query.clone();
        for position in (0..long_substituted.len()).step_by(97) {
            let base = &mut long_substituted[position];
            *base = if *base == b'A' { b'C' } else { b'A' };
        }
        let mut long_corrupted = long_query.clone();
        for base in &mut long_corrupted[5000..8500] {
            *base = if *base == b'A' { b'C' } else { b'A' };
        }
        let mut long_windows = Vec::new();
        let mut start = 0i32;
        while start < 20000 {
            let end = if start + 399 < 19999 {
                start + 399
            } else {
                19999
            };
            long_windows.push(WindowSpec {
                x_start: start,
                x_end: end,
                y_start: start,
                y_end: end,
                error: 4,
                threshold: 0,
            });
            start += 400;
        }

        vec![
            ReadScenario {
                name: "W1_exact",
                reads: vec![tile_query.clone(), tile_query],
                x_start: 0,
                x_end: 1399,
                y_start: 0,
                y_end: 1399,
                y_strand: 0,
                non_homopolymer_errors: 0,
                gaps: vec![(0, 0), (1399, 0)],
                windows: exact_tile,
                hits: vec![(0, 0), (350, 350), (700, 700), (1050, 1050), (1399, 1399)],
                window_length: 350,
            },
            ReadScenario {
                name: "W2_subs",
                reads: vec![query.clone(), substituted],
                x_start: 200,
                x_end: 1799,
                y_start: 800,
                y_end: 2399,
                y_strand: 0,
                non_homopolymer_errors: 16,
                gaps: vec![(200, 0), (1799, 0)],
                windows: subs_windows,
                hits: vec![(200, 800), (600, 1200), (1000, 1600), (1400, 2000)],
                window_length: 400,
            },
            ReadScenario {
                name: "W3_reverse",
                reads: vec![query, reversed],
                x_start: 200,
                x_end: 1799,
                y_start: 800,
                y_end: 2399,
                y_strand: 1,
                non_homopolymer_errors: 0,
                gaps: vec![(200, 0), (1799, 0)],
                windows: reverse_windows,
                hits: vec![(200, 800), (600, 1200), (1000, 1600), (1400, 2000)],
                window_length: 400,
            },
            ReadScenario {
                name: "L_mixed_exact",
                reads: vec![long_query.clone(), long_exact],
                x_start: 0,
                x_end: 19999,
                y_start: 0,
                y_end: 19999,
                y_strand: 0,
                non_homopolymer_errors: 1,
                gaps: vec![(0, 0), (19999, 0)],
                windows: long_windows.clone(),
                hits: vec![(1000, 1000), (5000, 4000)],
                window_length: 400,
            },
            ReadScenario {
                name: "L_mixed_subs",
                reads: vec![long_query.clone(), long_substituted],
                x_start: 0,
                x_end: 19999,
                y_start: 0,
                y_end: 19999,
                y_strand: 0,
                non_homopolymer_errors: 16,
                gaps: vec![(0, 0), (19999, 0)],
                windows: long_windows.clone(),
                hits: vec![(1000, 1000), (5000, 4000)],
                window_length: 400,
            },
            ReadScenario {
                name: "X_mixed_corrupt",
                reads: vec![long_query, long_corrupted],
                x_start: 0,
                x_end: 19999,
                y_start: 0,
                y_end: 19999,
                y_strand: 0,
                non_homopolymer_errors: 1,
                gaps: vec![(0, 0), (19999, 0)],
                windows: long_windows,
                hits: vec![(1000, 1000), (12000, 12000)],
                window_length: 400,
            },
        ]
    }

    fn generate_lines() -> Vec<String> {
        let mut lines = Vec::new();
        for scenario in read_scenarios() {
            let variants: [(&str, u32, i64, i64); 4] = [
                ("plain", scenario.non_homopolymer_errors, 0, 10),
                ("prefix", scenario.non_homopolymer_errors, 3, 10),
                ("hk1000", scenario.non_homopolymer_errors, 0, 1000),
                (
                    "nhpe1",
                    if scenario.non_homopolymer_errors != 0 {
                        scenario.non_homopolymer_errors
                    } else {
                        1
                    },
                    0,
                    10,
                ),
            ];
            for (suffix, non_homopolymer_errors, chain_prefix, h_kmer_hit) in variants {
                run_fast_cigar(
                    &scenario,
                    suffix,
                    non_homopolymer_errors,
                    chain_prefix,
                    h_kmer_hit,
                    &mut lines,
                );
            }
        }
        lines
    }

    #[test]
    #[ignore = "requires reference/ec_fast_cigar_probe.out from the C++ probe"]
    fn fast_cigar_matches_the_reference() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("reference")
            .join("ec_fast_cigar_probe.out");
        let contents = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("cannot read {}: {}", path.display(), error));
        let golden: Vec<&str> = contents.lines().collect();
        let generated = generate_lines();

        assert_eq!(
            generated.len(),
            golden.len(),
            "the probe printed {} lines, the replay generated {}",
            golden.len(),
            generated.len()
        );
        for (index, (expected, actual)) in golden.iter().zip(generated.iter()).enumerate() {
            assert_eq!(
                expected, actual,
                "line {} differs\nexpected: {}\nactual:   {}",
                index + 1,
                expected,
                actual
            );
        }
    }

    #[test]
    #[ignore = "requires reference/ec_return_chain_probe.out from the C++ probe"]
    fn return_chain_matches_the_reference() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("reference")
            .join("ec_return_chain_probe.out");
        let contents = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("cannot read {}: {}", path.display(), error));
        let golden: Vec<&str> = contents.lines().collect();
        let generated = generate_return_chain_lines();

        assert_eq!(
            generated.len(),
            golden.len(),
            "the probe printed {} lines, the replay generated {}",
            golden.len(),
            generated.len()
        );
        for (index, (expected, actual)) in golden.iter().zip(generated.iter()).enumerate() {
            assert_eq!(
                expected, actual,
                "line {} differs\nexpected: {}\nactual:   {}",
                index + 1,
                expected,
                actual
            );
        }
    }

    fn generate_return_chain_lines() -> Vec<String> {
        let mut lines = Vec::new();

        let mut hits = Vec::new();
        for index in 0..30u32 {
            hits.push(hit(1, 100 + index * 10, 100 + index * 10));
        }
        run_chain("colinear30", 0, &hits, &mut lines);

        let mut hits = Vec::new();
        for index in 0..10u32 {
            hits.push(hit(1, 100 + index * 10, 100 + index * 10));
        }
        for index in 0..10u32 {
            hits.push(hit(1, 200 + index * 200, 200 + index * 250));
        }
        for index in 0..5u32 {
            hits.push(hit(1, 2300 + index * 10, 2800 + index * 10));
        }
        run_chain("diagonal_gap", 0, &hits, &mut lines);

        let mut hits = Vec::new();
        for index in 0..10u32 {
            hits.push(hit(1, 100 + index * 10, 100 + index * 10));
        }
        for index in 0..6u32 {
            hits.push(hit(1, 800 + index * 10, 900 + index * 10));
        }
        run_chain("max_distance_break", 0, &hits, &mut lines);

        let mut hits = Vec::new();
        for index in 0..20u32 {
            hits.push(hit(1, 100 + index * 10, 100 + index * 10));
        }
        for index in 0..5u32 {
            let y = (290 - (index as i32) * 200) as u32;
            hits.push(hit(1, 310 + index * 10, y));
        }
        run_chain("tail_drop", 0, &hits, &mut lines);

        let mut hits = Vec::new();
        for index in 0..10u32 {
            hits.push(hit(1, 100 + index * 10, 100 + index * 10));
        }
        for index in 0..5u32 {
            hits.push(hit(2, 100 + index * 10, 100 + index * 10));
        }
        run_chain("two_groups", 0, &hits, &mut lines);

        let mut hits = Vec::new();
        for index in 0..5u32 {
            hits.push(hit(2, 50 + index * 10, 50 + index * 10));
        }
        for index in 0..12u32 {
            hits.push(hit(1, 100 + index * 10, 100 + index * 10));
        }
        for index in 0..4u32 {
            hits.push(hit(3, 300 + index * 10, 300 + index * 10));
        }
        run_chain("group_middle", 5, &hits, &mut lines);

        run_chain("single", 0, &[hit(1, 500, 500)], &mut lines);

        let hits = vec![
            hit(1, 100, 100),
            hit(1, 110, 126),
            hit(1, 120, 137),
            hit(1, 130, 147),
            hit(1, 140, 157),
        ];
        run_chain("long_gap_boundary", 0, &hits, &mut lines);

        let mut hits = Vec::new();
        for index in 0..5u32 {
            hits.push(hit(1, 100 + index * 10, 100 + index * 10));
        }
        for index in 0..3u32 {
            hits.push(hit(0x7fff_ffff, 200 + index * 10, 200 + index * 10));
        }
        for index in 0..2u32 {
            hits.push(hit(1, 300 + index * 10, 300 + index * 10));
        }
        run_chain("sentinel_stop", 0, &hits, &mut lines);

        let mut hits = Vec::new();
        for index in 0..40u32 {
            let mut x = 100 + index * 10;
            let mut y = 100 + index * 10;
            if index % 7 == 3 {
                y += 200;
            }
            if index % 11 == 5 {
                x += 15;
            }
            hits.push(hit(1, x, y));
        }
        run_chain("noisy", 0, &hits, &mut lines);

        // ------------------------------------------------- fast cigar wrapper
        let mut generator = Lcg(987654321);
        let long_query = generator.sequence(20000);
        let long_exact = long_query.clone();

        let mut long_windows = Vec::new();
        let mut start = 0i32;
        while start < 20000 {
            let end = if start + 399 < 19999 {
                start + 399
            } else {
                19999
            };
            long_windows.push(WindowSpec {
                x_start: start,
                x_end: end,
                y_start: start,
                y_end: end,
                error: 4,
                threshold: 0,
            });
            start += 400;
        }

        let mut hits = Vec::new();
        for index in 10..=30u32 {
            let x = index * 100;
            if index == 11 {
                hits.push(hit(1, x + 50, 50));
            }
            if index == 21 {
                hits.push(hit(1, x + 50, 2800));
            }
            hits.push(hit(1, x, x));
        }
        let mut windows = Vec::new();
        for k in 0..10i32 {
            let start = k * 300;
            windows.push(WindowSpec {
                x_start: start,
                x_end: start + 299,
                y_start: start,
                y_end: start + 299,
                error: 4,
                threshold: 0,
            });
        }
        run_wrapper(
            &[long_query.clone(), long_exact.clone()],
            0,
            2999,
            0,
            2999,
            &[(0, 0), (2999, 0)],
            &windows,
            &WrapperSpec {
                name: "refined_decoys",
                non_homopolymer_errors: 1,
                window_length: 300,
                hits,
            },
            &mut lines,
        );

        let mut hits = Vec::new();
        for index in 1..=19u32 {
            hits.push(hit(1, index * 1000, index * 1000));
        }
        hits.push(hit(1, 19100, 100));
        hits.push(hit(1, 19200, 200));
        hits.push(hit(1, 19300, 300));
        run_wrapper(
            &[long_query, long_exact],
            0,
            19999,
            0,
            19999,
            &[(0, 0), (19999, 0)],
            &long_windows,
            &WrapperSpec {
                name: "tail_dropped",
                non_homopolymer_errors: 1,
                window_length: 400,
                hits,
            },
            &mut lines,
        );

        lines
    }
}

//! `hc_ovlp_base_direct`: the block-wise alignment of one overlap in the new
//! correction path.
//!
//! A chain is a list of seeds, each a pair of matching positions. Consecutive
//! seeds bracket a block of the overlap; the block between seed `l` and seed
//! `i` is aligned with [`hc_aln_exz_adv_hc`], and blocks that fail become
//! unaligned windows so that a later pass can retry them from scratch.
//!
//! Two details of the reference shape the port:
//!
//! * A block whose first seed is missing (`l < 0`) is an extension from the
//!   left read end (`mode 2`), and the block after the last seed (`i == ch_n`)
//!   is an extension from the right read end (`mode 1`); those two are widened
//!   to the read ends by [`adjust_ext_offset`] before alignment.
//! * When the caller has no preferred mode (`pre_mode < 0`) and the overlap
//!   already carries a zero-error window list that tiles it exactly, the whole
//!   overlap is accepted as a single match without running any DP, but only if
//!   the widened diagonal still starts and ends exactly on the recorded
//!   windows.

use crate::correction::aln_adv::hc_aln_exz_adv_hc;
use crate::correction::cal_exz::{
    adjust_ext_offset, push_alignment_window, push_unmapped_window, set_exact_exz,
    FORCED_SINGLE_LENGTH, MAXIMUM_SINGLE_ERROR, MAXIMUM_SINGLE_LENGTH, MODE_BACKWARD_EXTENSION,
    MODE_FORWARD_EXTENSION, MODE_GLOBAL, MODE_SEMI_GLOBAL,
};
use crate::correction::hc_align::BitExtz;
use crate::correction::window::WindowCigar;
use crate::overlap::region::{KmerHit, OverlapRegion};
use crate::correction::target::TargetSource;

/// `hc_ovlp_base_direct`: aligns every block of a chain, appending one window
/// per block to `output`.
///
/// `pre_mode` selects the blocks to align: negative aligns every block
/// (including both end extensions), `0` skips both extensions, `1` skips the
/// left one, `2` skips the right one. Any other value aligns every block.
/// `query_total_length` and `target_total_length` are the lengths of the two
/// reads, used to extend the end blocks and to clamp them.
///
/// `read_identifier` of the reference is debug-only and is not carried over.
#[allow(clippy::too_many_arguments)]
pub fn hc_ovlp_base_direct(
    region: &OverlapRegion,
    hits: &[KmerHit],
    hit_count: i64,
    window_length: i64,
    targets: &TargetSource,
    query: &[u8],
    extz: &mut BitExtz,
    output: &mut WindowCigar,
    error_rate: f64,
    query_total_length: i64,
    target_total_length: i64,
    pre_mode: i64,
) {
    // A zero-error window list that tiles the overlap exactly needs no DP: the
    // aligned region is widened to the read ends along the same diagonal and
    // accepted when the widening did not move either query end. The reference
    // compares successive `y_end` values rather than `y_start`, so a window
    // whose two spans differ is still accepted as exact.
    if pre_mode < 0 && region.non_homopolymer_errors == 0 {
        let windows = &region.window_cigar.windows;
        let window_count = windows.len();
        if window_count > 0 {
            let mut index = 1;
            while index < window_count {
                let previous = windows[index - 1];
                let current = windows[index];
                let collinear = current.error == 0
                    && previous.error == 0
                    && i64::from(current.x_start) == i64::from(previous.x_end) + 1
                    && i64::from(current.y_end)
                        == i64::from(previous.y_end)
                            + (i64::from(current.x_end) - i64::from(previous.x_end));
                if !collinear {
                    break;
                }
                index += 1;
            }

            if index >= window_count {
                let first = windows[0];
                let last = windows[window_count - 1];
                let mut query_start = i64::from(first.x_start);
                let mut query_end = i64::from(last.x_end);
                let mut target_start = i64::from(first.y_end)
                    - (i64::from(first.x_end) - i64::from(first.x_start));
                let mut target_end = i64::from(last.y_end);

                if query_start <= target_start {
                    target_start -= query_start;
                    query_start = 0;
                } else {
                    query_start -= target_start;
                    target_start = 0;
                }

                let query_remaining = query_total_length - query_end - 1;
                let target_remaining = target_total_length - target_end - 1;
                if query_remaining <= target_remaining {
                    query_end = query_total_length - 1;
                    target_end += query_remaining;
                } else {
                    target_end = target_total_length - 1;
                    query_end += target_remaining;
                }

                if query_start == i64::from(first.x_start) && query_end == i64::from(last.x_end) {
                    set_exact_exz(extz, query_start, query_end + 1, target_start, target_end + 1);
                    push_alignment_window(output, extz);
                    return;
                }
            }
        }
    }

    let mut first_block = 0i64;
    let mut last_block = hit_count;
    if pre_mode == 0 {
        first_block = 1;
        last_block = hit_count - 1;
    } else if pre_mode == 1 {
        first_block = 1;
    } else if pre_mode == 2 {
        last_block = hit_count - 1;
    }

    let mut previous = first_block - 1;
    let mut index = first_block;
    while index <= last_block {
        let mut query_start;
        let mut query_end;
        let mut target_start;
        let mut target_end;
        let mode;

        if previous >= 0 {
            let hit = &hits[previous as usize];
            query_start = i64::from(hit.self_offset);
            target_start = i64::from(hit.offset);
        } else {
            query_start = 0;
            target_start = -1;
        }

        if index < hit_count {
            let hit = &hits[index as usize];
            query_end = i64::from(hit.self_offset);
            target_end = i64::from(hit.offset);
        } else {
            query_end = query_total_length;
            target_end = -1;
        }

        if target_start != -1 && target_end != -1 {
            mode = MODE_GLOBAL;
        } else if target_start != -1 {
            mode = MODE_FORWARD_EXTENSION;
        } else if target_end != -1 {
            mode = MODE_BACKWARD_EXTENSION;
        } else {
            // No seed in range: the caller pushes a single unmapped window.
            mode = MODE_SEMI_GLOBAL;
        }

        if mode == MODE_FORWARD_EXTENSION || mode == MODE_BACKWARD_EXTENSION {
            adjust_ext_offset(
                &mut query_start,
                &mut query_end,
                &mut target_start,
                &mut target_end,
                query_total_length,
                target_total_length,
                0,
                mode,
            );
        }

        let aligned = hc_aln_exz_adv_hc(
            region,
            targets,
            extz,
            query,
            query_start,
            query_end,
            target_start,
            target_end,
            mode,
            window_length,
            query_total_length,
            error_rate,
            MAXIMUM_SINGLE_LENGTH,
            MAXIMUM_SINGLE_ERROR,
            FORCED_SINGLE_LENGTH,
            -1,
            output,
        );

        if !aligned {
            // `push_unmap_alnw` takes `[q0, q1]` on the query and `[t0, t1]` on
            // the target, so an exclusive block end loses one base; an unset
            // target end (mode 3) is recorded as `-1`/`-2`.
            push_unmapped_window(
                output,
                query_start as i32,
                query_end - 1,
                target_start,
                target_end - 1,
                mode,
            );
        }

        previous = index;
        index += 1;
    }
}

/// Replay of the reference golden probe.
///
/// `/tmp/ec_ovlpbase_probe.cpp` links against the unmodified hifiasm objects
/// and runs ten scenarios through `hc_ovlp_base_direct`: exact window tilings
/// with and without the shortcut, chains that walk each `pre_mode` selection,
/// a block too long to align, a target shorter than the recorded windows, and
/// a reverse-strand overlap. The replay rebuilds the same reads, windows and
/// seeds and compares every line, including the alignment left over from an
/// earlier call when a driver returns without touching `extz`.
#[cfg(test)]
mod probe_goldens {
    use super::*;
    use crate::reads::read_store::ReadStore;
    use crate::overlap::region::{FakeCigar, KmerHit, OverlapRegion};
    use crate::correction::target::TargetSource;

    #[derive(Clone)]
    struct WindowSpec {
        x_start: i32,
        x_end: i32,
        y_start: i32,
        y_end: i32,
        error: i16,
    }

    struct Scenario {
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
        pre_modes: Vec<i64>,
    }

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

    fn exz_suffix(extz: &BitExtz) -> String {
        let result = &extz.result;
        let mut out = format!(
            " err={} thre={} ps={} pe={} ts={} te={} pl={} tl={} cn={} C",
            result.error,
            result.threshold,
            result.pattern_start,
            result.pattern_end,
            result.text_start,
            result.text_end,
            result.pattern_length,
            result.text_length,
            extz.cigar.packed.len()
        );
        for cell in &extz.cigar.packed {
            out.push_str(&format!(
                " {}{}",
                cell & 0x3fff,
                b"EMDI"[(cell >> 14) as usize] as char
            ));
        }
        out
    }

    fn append_window_lines(lines: &mut Vec<String>, output: &WindowCigar) {
        lines.push(format!(
            " WN n={} cn={}",
            output.windows.len(),
            output.cells.len()
        ));
        for (index, window) in output.windows.iter().enumerate() {
            let mut line = format!(
                "  W {} x=[{},{}] y=[{},{}] eb={} ee={} err={} thre={} cidx={} clen={} :",
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
                let cell = output.cells[(window.cigar_start + offset) as usize];
                line.push_str(&format!(
                    " {}{}",
                    cell & 0x3fff,
                    b"EMDI"[(cell >> 14) as usize] as char
                ));
            }
            lines.push(line);
        }
    }

    fn tiling_windows(offset: i32, length: i32, error: i16) -> Vec<WindowSpec> {
        (0..4)
            .map(|index| {
                let start = offset + index * length;
                let end = start + length - 1;
                WindowSpec {
                    x_start: start,
                    x_end: end,
                    y_start: start,
                    y_end: end,
                    error,
                }
            })
            .collect()
    }

    fn run_scenario(scenario: &Scenario, extz: &mut BitExtz) -> Vec<String> {
        let mut reads = ReadStore::default();
        for (index, sequence) in scenario.reads.iter().enumerate() {
            reads.push_read(&format!("read{}", index), sequence);
        }

        let mut region = OverlapRegion {
            x_identifier: 0,
            x_position_start: scenario.x_start,
            x_position_end: scenario.x_end,
            y_identifier: 1,
            y_position_start: scenario.y_start,
            y_position_end: scenario.y_end,
            y_strand: scenario.y_strand,
            overlap_length: scenario.x_end + 1 - scenario.x_start,
            non_homopolymer_errors: scenario.non_homopolymer_errors,
            ..OverlapRegion::default()
        };
        region.fake_cigar = FakeCigar::default();
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
                0,
            );
        }

        let hits: Vec<KmerHit> = scenario
            .hits
            .iter()
            .map(|&(self_offset, offset)| KmerHit {
                read_identifier: 0,
                reverse_strand: false,
                offset,
                self_offset,
                count: 0,
            })
            .collect();

        let query = &scenario.reads[0];
        let query_total_length = query.len() as i64;
        let target_total_length = scenario.reads[1].len() as i64;
        let mut lines = Vec::new();

        for (index, &pre_mode) in scenario.pre_modes.iter().enumerate() {
            let mut output = WindowCigar::default();
            hc_ovlp_base_direct(
                &region,
                &hits,
                hits.len() as i64,
                400,
                &TargetSource::Reads(&reads),
                query,
                extz,
                &mut output,
                0.04,
                query_total_length,
                target_total_length,
                pre_mode,
            );
            lines.push(format!(
                "OVB {} i={} pre={} ch_n={} nhpe={} ql={} tl={}{}",
                scenario.name,
                index,
                pre_mode,
                hits.len(),
                scenario.non_homopolymer_errors,
                query_total_length,
                target_total_length,
                exz_suffix(extz)
            ));
            append_window_lines(&mut lines, &output);
        }

        lines
    }

    fn scenarios() -> Vec<Scenario> {
        let mut generator = Lcg(987654321);
        let query = generator.sequence(2000);
        let prefix = generator.sequence(600);
        let suffix = generator.sequence(400);
        let tile_query = generator.sequence(1400);
        let long_query = generator.sequence(20000);
        let long_target = generator.sequence(30000);

        let mut exact_target = prefix.clone();
        exact_target.extend_from_slice(&query);
        exact_target.extend_from_slice(&suffix);

        let mut substituted = exact_target.clone();
        for index in (0..query.len()).step_by(97) {
            let position = prefix.len() + index;
            substituted[position] = if substituted[position] == b'A' {
                b'C'
            } else {
                b'A'
            };
        }

        let reversed = reverse_complement(&exact_target);
        let short_target = tile_query[..1000].to_vec();

        let exact_tile = tiling_windows(0, 350, 0);
        let shifted_tile = tiling_windows(100, 325, 0);
        let mut error_tile = tiling_windows(0, 350, 0);
        error_tile[2].error = 1;

        let subs_windows = vec![
            WindowSpec { x_start: 200, x_end: 599, y_start: 800, y_end: 1199, error: 4 },
            WindowSpec { x_start: 600, x_end: 999, y_start: 1200, y_end: 1599, error: 4 },
            WindowSpec { x_start: 1000, x_end: 1399, y_start: 1600, y_end: 1999, error: 4 },
            WindowSpec { x_start: 1400, x_end: 1799, y_start: 2000, y_end: 2399, error: 5 },
        ];
        let reverse_windows = vec![
            WindowSpec { x_start: 200, x_end: 599, y_start: 800, y_end: 1199, error: 0 },
            WindowSpec { x_start: 600, x_end: 999, y_start: 1200, y_end: 1599, error: 0 },
            WindowSpec { x_start: 1000, x_end: 1399, y_start: 1600, y_end: 1999, error: 0 },
            WindowSpec { x_start: 1400, x_end: 1799, y_start: 2000, y_end: 2399, error: 0 },
        ];

        vec![
            Scenario {
                name: "A_exact_tiling",
                reads: vec![tile_query.clone(), tile_query.clone()],
                x_start: 0, x_end: 1399, y_start: 0, y_end: 1399, y_strand: 0,
                non_homopolymer_errors: 0,
                gaps: vec![(0, 0), (1399, 0)],
                windows: exact_tile.clone(),
                hits: vec![],
                pre_modes: vec![-1, 0, 1, 2],
            },
            Scenario {
                name: "B_tiling_shift",
                reads: vec![tile_query.clone(), tile_query.clone()],
                x_start: 0, x_end: 1399, y_start: 0, y_end: 1399, y_strand: 0,
                non_homopolymer_errors: 0,
                gaps: vec![(0, 0), (1399, 0)],
                windows: shifted_tile,
                hits: vec![(100, 100), (425, 425), (750, 750), (1075, 1075)],
                pre_modes: vec![-1, 0, 1, 2, 7],
            },
            Scenario {
                name: "C_tiling_error",
                reads: vec![tile_query.clone(), tile_query.clone()],
                x_start: 0, x_end: 1399, y_start: 0, y_end: 1399, y_strand: 0,
                non_homopolymer_errors: 0,
                gaps: vec![(0, 0), (1399, 0)],
                windows: error_tile,
                hits: vec![],
                pre_modes: vec![-1],
            },
            Scenario {
                name: "D_tiling_nonhpe",
                reads: vec![tile_query.clone(), tile_query.clone()],
                x_start: 0, x_end: 1399, y_start: 0, y_end: 1399, y_strand: 0,
                non_homopolymer_errors: 7,
                gaps: vec![(0, 0), (1399, 0)],
                windows: exact_tile.clone(),
                hits: vec![(0, 0), (350, 350), (700, 700), (1050, 1050)],
                pre_modes: vec![-1],
            },
            Scenario {
                name: "E_subs",
                reads: vec![query.clone(), substituted],
                x_start: 200, x_end: 1799, y_start: 800, y_end: 2399, y_strand: 0,
                non_homopolymer_errors: 16,
                gaps: vec![(200, 0), (1799, 0)],
                windows: subs_windows,
                hits: vec![(200, 800), (600, 1200), (1000, 1600), (1400, 2000)],
                pre_modes: vec![-1, 0, 1, 2],
            },
            Scenario {
                name: "F_long_fail",
                reads: vec![long_query, long_target],
                x_start: 0, x_end: 19999, y_start: 0, y_end: 19999, y_strand: 0,
                non_homopolymer_errors: 0,
                gaps: vec![],
                windows: vec![],
                hits: vec![(1000, 1000)],
                pre_modes: vec![-1],
            },
            Scenario {
                name: "G_target_clamp",
                reads: vec![tile_query.clone(), short_target],
                x_start: 0, x_end: 1399, y_start: 0, y_end: 999, y_strand: 0,
                non_homopolymer_errors: 0,
                gaps: vec![(0, 0), (1399, 0)],
                windows: exact_tile.clone(),
                hits: vec![(0, 0), (700, 700)],
                pre_modes: vec![-1],
            },
            Scenario {
                name: "H_no_iteration",
                reads: vec![tile_query.clone(), tile_query.clone()],
                x_start: 0, x_end: 1399, y_start: 0, y_end: 1399, y_strand: 0,
                non_homopolymer_errors: 0,
                gaps: vec![(0, 0), (1399, 0)],
                windows: vec![],
                hits: vec![],
                pre_modes: vec![0, 1, 2],
            },
            Scenario {
                name: "I_reverse",
                reads: vec![query.clone(), reversed],
                x_start: 200, x_end: 1799, y_start: 800, y_end: 2399, y_strand: 1,
                non_homopolymer_errors: 0,
                gaps: vec![(200, 0), (1799, 0)],
                windows: reverse_windows,
                hits: vec![(200, 800), (600, 1200), (1000, 1600), (1400, 2000)],
                pre_modes: vec![-1, 0],
            },
            Scenario {
                name: "J_both_unset",
                reads: vec![tile_query.clone(), tile_query],
                x_start: 0, x_end: 1399, y_start: 0, y_end: 1399, y_strand: 0,
                non_homopolymer_errors: 0,
                gaps: vec![(0, 0), (1399, 0)],
                windows: vec![],
                hits: vec![],
                pre_modes: vec![-1, 3, 9],
            },
        ]
    }

    #[test]
    #[ignore = "requires reference/ec_ovlpbase_probe.out from the C++ probe"]
    fn overlap_base_direct_matches_the_reference() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("reference")
            .join("ec_ovlpbase_probe.out");
        let contents = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("cannot read {}: {}", path.display(), error));
        let golden: Vec<&str> = contents.lines().collect();
        let mut generated = Vec::new();
        let mut extz = BitExtz::new();
        for scenario in scenarios() {
            generated.extend(run_scenario(&scenario, &mut extz));
        }

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
                index + 1, expected, actual
            );
        }
    }
}


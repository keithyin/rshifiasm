//! `hc_aln_exz_adv` / `hc_aln_exz_adv_hc`: the alignment driver that turns one
//! chain block of the new correction path into a window.
//!
//! The driver first tries to show the block is identical (via the exact
//! comparison or, in the `_hc` variant, via the window-based error estimate),
//! then walks a fixed ladder of thresholds – the estimate, the expected error,
//! twice that, half the block, and finally the full allowance for short blocks
//! – stopping at the first band that aligns. Every failed attempt leaves its
//! `pthre` behind, so later attempts with a narrower band are skipped inside
//! `cal_exz_infi_adv`.

use crate::correction::cal_exz::{
    cal_estimate_err, cal_estimate_err_hc, cal_exact_exz, cal_exz_infi_adv, push_alignment_window,
    scale_ed_thre, set_exact_exz, ExzAdvState, MODE_SEMI_GLOBAL,
};
use crate::correction::hc_align::BitExtz;
use crate::correction::window::WindowCigar;
use crate::overlap::region::OverlapRegion;
use crate::reads::read_store::ReadStore;

/// `hc_aln_exz_adv`: aligns one block, estimating the error from the windows
/// only when the block spans a whole window or more.
///
/// `estimate_error` is negative when the caller has no estimate yet; the
/// function then computes one. The alignment is appended to `output`, whose
/// previous window is merged with the new one when they are adjacent.
#[allow(clippy::too_many_arguments)]
pub fn hc_aln_exz_adv(
    region: &OverlapRegion,
    reads: &ReadStore,
    extz: &mut BitExtz,
    query: &[u8],
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    mode: i64,
    window_length: i64,
    query_total_length: i64,
    error_rate: f64,
    maximum_length: i64,
    maximum_error: i64,
    forced_length: i64,
    estimate_error: i64,
    output: &mut WindowCigar,
) -> bool {
    extz.clear_align();
    extz.result.threshold = 0;
    let mode = if target_start == -1 && target_end == -1 {
        MODE_SEMI_GLOBAL
    } else {
        mode
    };
    let query_length = query_end - query_start;
    if query_length == 0 && target_end - target_start == 0 {
        return true;
    }
    if query_length <= 0 || target_end - target_start <= 0 {
        return false;
    }

    let mut estimate_error = estimate_error;
    if estimate_error < 0 {
        estimate_error = if query_length > window_length {
            cal_estimate_err(region, window_length, query_start, query_end, error_rate)
        } else {
            (query_length as f64 * error_rate) as i64
        };
    }

    let mut state = ExzAdvState::default();
    if query_length <= 16
        && cal_exact_exz(
            region,
            reads,
            extz,
            query,
            query_start,
            query_end,
            target_start,
            target_end,
            &mut state,
            mode,
        )
    {
        push_alignment_window(output, extz);
        return true;
    }

    attempt_thresholds(
        region,
        reads,
        extz,
        query,
        query_start,
        query_end,
        target_start,
        target_end,
        mode,
        query_total_length,
        error_rate,
        maximum_length,
        maximum_error,
        forced_length,
        estimate_error,
        &mut state,
        output,
    )
}

/// `hc_aln_exz_adv_hc`: the window-based variant used by `gen_hc_r_alin`.
///
/// Unlike [`hc_aln_exz_adv`], the estimate always comes from
/// `cal_estimate_err_hc`, and an estimate of zero short-circuits to a
/// single-match window when the windows tile the target exactly.
#[allow(clippy::too_many_arguments)]
pub fn hc_aln_exz_adv_hc(
    region: &OverlapRegion,
    reads: &ReadStore,
    extz: &mut BitExtz,
    query: &[u8],
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    mode: i64,
    window_length: i64,
    query_total_length: i64,
    error_rate: f64,
    maximum_length: i64,
    maximum_error: i64,
    forced_length: i64,
    estimate_error: i64,
    output: &mut WindowCigar,
) -> bool {
    extz.clear_align();
    extz.result.threshold = 0;
    let mode = if target_start == -1 && target_end == -1 {
        MODE_SEMI_GLOBAL
    } else {
        mode
    };
    let query_length = query_end - query_start;
    if query_length == 0 && target_end - target_start == 0 {
        return true;
    }
    if query_length <= 0 || target_end - target_start <= 0 {
        return false;
    }

    let mut estimate_error = estimate_error;
    let mut full = false;
    if estimate_error < 0 {
        let (value, exact) = cal_estimate_err_hc(
            region,
            window_length,
            query_start,
            query_end,
            target_start,
            target_end,
            error_rate,
        );
        estimate_error = value;
        full = exact;
    }

    let mut state = ExzAdvState::default();
    if estimate_error == 0 {
        if full {
            set_exact_exz(extz, query_start, query_end, target_start, target_end);
            push_alignment_window(output, extz);
            return true;
        }
        if cal_exact_exz(
            region,
            reads,
            extz,
            query,
            query_start,
            query_end,
            target_start,
            target_end,
            &mut state,
            mode,
        ) {
            push_alignment_window(output, extz);
            return true;
        }
    }

    attempt_thresholds(
        region,
        reads,
        extz,
        query,
        query_start,
        query_end,
        target_start,
        target_end,
        mode,
        query_total_length,
        error_rate,
        maximum_length,
        maximum_error,
        forced_length,
        estimate_error,
        &mut state,
        output,
    )
}

/// The threshold ladder shared by both drivers.
///
/// `state` carries `pts`/`pte`/`pthre` across the attempts; the reference
/// initialises them to `-1` once, before the first attempt, and every
/// `cal_exz_infi_adv` call updates them.
#[allow(clippy::too_many_arguments)]
fn attempt_thresholds(
    region: &OverlapRegion,
    reads: &ReadStore,
    extz: &mut BitExtz,
    query: &[u8],
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    mode: i64,
    query_total_length: i64,
    error_rate: f64,
    maximum_length: i64,
    maximum_error: i64,
    forced_length: i64,
    estimate_error: i64,
    state: &mut ExzAdvState,
    output: &mut WindowCigar,
) -> bool {
    let query_length = query_end - query_start;
    if query_length > maximum_length || (estimate_error >> 1) > maximum_error {
        return false;
    }

    let mut threshold = scale_ed_thre(estimate_error as u32, maximum_error as u32) as i64;
    if attempt(
        region,
        reads,
        extz,
        query,
        query_start,
        query_end,
        target_start,
        target_end,
        mode,
        query_total_length,
        threshold,
        state,
        output,
    ) {
        return true;
    }

    let mut previous = threshold;
    threshold = (query_length as f64 * error_rate) as i64;
    threshold = scale_ed_thre(threshold as u32, maximum_error as u32) as i64;
    if threshold > previous
        && attempt(
            region,
            reads,
            extz,
            query,
            query_start,
            query_end,
            target_start,
            target_end,
            mode,
            query_total_length,
            threshold,
            state,
            output,
        )
    {
        return true;
    }

    previous = threshold;
    threshold <<= 1;
    threshold = scale_ed_thre(threshold as u32, maximum_error as u32) as i64;
    if threshold > previous
        && attempt(
            region,
            reads,
            extz,
            query,
            query_start,
            query_end,
            target_start,
            target_end,
            mode,
            query_total_length,
            threshold,
            state,
            output,
        )
    {
        return true;
    }

    previous = threshold;
    threshold = (query_length as f64 * 0.51) as i64;
    threshold = scale_ed_thre(threshold as u32, maximum_error as u32) as i64;
    if threshold > previous
        && attempt(
            region,
            reads,
            extz,
            query,
            query_start,
            query_end,
            target_start,
            target_end,
            mode,
            query_total_length,
            threshold,
            state,
            output,
        )
    {
        return true;
    }

    if query_length <= forced_length {
        threshold = maximum_error;
        if attempt(
            region,
            reads,
            extz,
            query,
            query_start,
            query_end,
            target_start,
            target_end,
            mode,
            query_total_length,
            threshold,
            state,
            output,
        ) {
            return true;
        }
    }

    false
}

/// One `cal_exz_infi_adv` attempt; on success the alignment is pushed.
#[allow(clippy::too_many_arguments)]
fn attempt(
    region: &OverlapRegion,
    reads: &ReadStore,
    extz: &mut BitExtz,
    query: &[u8],
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    mode: i64,
    query_total_length: i64,
    threshold: i64,
    state: &mut ExzAdvState,
    output: &mut WindowCigar,
) -> bool {
    let aligned = cal_exz_infi_adv(
        region,
        reads,
        extz,
        query,
        query_start,
        query_end,
        target_start,
        target_end,
        state,
        threshold,
        query_total_length,
        mode,
    );
    if aligned {
        push_alignment_window(output, extz);
    }
    aligned
}

/// Replay of the reference golden probe.
///
/// `/tmp/ec_alnadv_probe.cpp` links against the unmodified hifiasm objects and
/// runs both drivers over the same four scenarios as the `cal_exz` probe, with
/// sixteen parameter combinations each. The output is stored under
/// `reference/`; the test rebuilds the reads, windows and call order and
/// compares every line, including the alignment left over from an earlier call
/// when a driver returns without touching the cigar.
#[cfg(test)]
mod probe_goldens {
    use super::*;
    use crate::overlap::region::{FakeCigar, OverlapRegion};
    use crate::reads::read_store::ReadStore;

    struct WindowSpec {
        x_start: i32,
        x_end: i32,
        y_start: i32,
        y_end: i32,
        error: i16,
    }

    struct DriverSpec {
        query_start: i64,
        query_end: i64,
        target_start: i64,
        target_end: i64,
        mode: i64,
        maximum_length: i64,
        maximum_error: i64,
        forced_length: i64,
        estimate_error: i64,
        use_hc: bool,
    }

    struct Scenario {
        name: &'static str,
        reads: Vec<Vec<u8>>,
        x_start: u32,
        x_end: u32,
        y_start: u32,
        y_end: u32,
        y_strand: u32,
        gaps: Vec<(u32, i32)>,
        windows: Vec<WindowSpec>,
        specs: Vec<DriverSpec>,
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

    fn run_scenario(scenario: &Scenario) -> Vec<String> {
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

        let query = &scenario.reads[0];
        let query_total_length = query.len() as i64;
        let mut extz = BitExtz::new();
        let mut lines = Vec::new();

        for (index, spec) in scenario.specs.iter().enumerate() {
            let mut output = WindowCigar::default();
            let aligned = if spec.use_hc {
                hc_aln_exz_adv_hc(
                    &region,
                    &reads,
                    &mut extz,
                    query,
                    spec.query_start,
                    spec.query_end,
                    spec.target_start,
                    spec.target_end,
                    spec.mode,
                    400,
                    query_total_length,
                    0.04,
                    spec.maximum_length,
                    spec.maximum_error,
                    spec.forced_length,
                    spec.estimate_error,
                    &mut output,
                )
            } else {
                hc_aln_exz_adv(
                    &region,
                    &reads,
                    &mut extz,
                    query,
                    spec.query_start,
                    spec.query_end,
                    spec.target_start,
                    spec.target_end,
                    spec.mode,
                    400,
                    query_total_length,
                    0.04,
                    spec.maximum_length,
                    spec.maximum_error,
                    spec.forced_length,
                    spec.estimate_error,
                    &mut output,
                )
            };
            lines.push(format!(
                "HCA {} {} i={} qs={} qe={} ts={} te={} mode={} maxl={} maxe={} force={} est={} -> ret={}{}",
                scenario.name,
                if spec.use_hc { "hc" } else { "pl" },
                index,
                spec.query_start,
                spec.query_end,
                spec.target_start,
                spec.target_end,
                spec.mode,
                spec.maximum_length,
                spec.maximum_error,
                spec.forced_length,
                spec.estimate_error,
                i32::from(aligned),
                exz_suffix(&extz)
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

        let mut exact_target = prefix.clone();
        exact_target.extend_from_slice(&query);
        exact_target.extend_from_slice(&suffix);

        let mut substituted = exact_target.clone();
        for index in (0..query.len()).step_by(97) {
            let position = prefix.len() + index;
            substituted[position] = if substituted[position] == b'A' { b'C' } else { b'A' };
        }

        let mut inserted = prefix.clone();
        inserted.extend_from_slice(&query[..700]);
        inserted.extend_from_slice(b"ACGTAC");
        inserted.extend_from_slice(&query[700..]);
        inserted.extend_from_slice(&suffix);

        let reversed = reverse_complement(&exact_target);

        let exact_windows = |errors: [i16; 4]| {
            vec![
                WindowSpec { x_start: 200, x_end: 599, y_start: 800, y_end: 1199, error: errors[0] },
                WindowSpec { x_start: 600, x_end: 999, y_start: 1200, y_end: 1599, error: errors[1] },
                WindowSpec { x_start: 1000, x_end: 1399, y_start: 1600, y_end: 1999, error: errors[2] },
                WindowSpec { x_start: 1400, x_end: 1799, y_start: 2000, y_end: 2399, error: errors[3] },
            ]
        };

        let maximum_length = 10000;
        let maximum_error = 2047;
        let forced_length = 512;
        let specs = vec![
            DriverSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, mode: 0, maximum_length, maximum_error, forced_length, estimate_error: -1, use_hc: true },
            DriverSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, mode: 0, maximum_length, maximum_error, forced_length, estimate_error: 0, use_hc: true },
            DriverSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, mode: 0, maximum_length: 500, maximum_error, forced_length, estimate_error: -1, use_hc: true },
            DriverSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, mode: 0, maximum_length, maximum_error: 30, forced_length, estimate_error: -1, use_hc: true },
            DriverSpec { query_start: 200, query_end: 1000, target_start: 800, target_end: 1600, mode: 1, maximum_length, maximum_error, forced_length, estimate_error: -1, use_hc: true },
            DriverSpec { query_start: 1000, query_end: 1800, target_start: 1600, target_end: 2400, mode: 2, maximum_length, maximum_error, forced_length, estimate_error: -1, use_hc: true },
            DriverSpec { query_start: 200, query_end: 1800, target_start: -1, target_end: -1, mode: 0, maximum_length, maximum_error, forced_length, estimate_error: -1, use_hc: true },
            DriverSpec { query_start: 600, query_end: 610, target_start: 1200, target_end: 1210, mode: 0, maximum_length, maximum_error, forced_length, estimate_error: -1, use_hc: true },
            DriverSpec { query_start: 600, query_end: 613, target_start: 1200, target_end: 1213, mode: 0, maximum_length, maximum_error, forced_length, estimate_error: -1, use_hc: true },
            DriverSpec { query_start: 700, query_end: 1300, target_start: 1306, target_end: 1906, mode: 0, maximum_length, maximum_error, forced_length, estimate_error: -1, use_hc: true },
            DriverSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, mode: 0, maximum_length, maximum_error, forced_length, estimate_error: 100, use_hc: true },
            DriverSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, mode: 0, maximum_length, maximum_error, forced_length, estimate_error: 5, use_hc: false },
            DriverSpec { query_start: 200, query_end: 1800, target_start: 800, target_end: 2400, mode: 0, maximum_length, maximum_error, forced_length, estimate_error: 0, use_hc: false },
            DriverSpec { query_start: 600, query_end: 610, target_start: 1200, target_end: 1210, mode: 0, maximum_length, maximum_error, forced_length, estimate_error: -1, use_hc: false },
            DriverSpec { query_start: 600, query_end: 613, target_start: 1200, target_end: 1213, mode: 0, maximum_length, maximum_error, forced_length, estimate_error: -1, use_hc: false },
            DriverSpec { query_start: 200, query_end: 1800, target_start: -1, target_end: -1, mode: 3, maximum_length, maximum_error, forced_length, estimate_error: -1, use_hc: false },
        ];

        vec![
            Scenario {
                name: "A_exact",
                reads: vec![query.clone(), exact_target.clone()],
                x_start: 200,
                x_end: 1799,
                y_start: 800,
                y_end: 2399,
                y_strand: 0,
                gaps: vec![(200, 0), (1799, 0)],
                windows: exact_windows([0, 0, 0, 0]),
                specs: clone_specs(&specs),
            },
            Scenario {
                name: "B_subs",
                reads: vec![query.clone(), substituted],
                x_start: 200,
                x_end: 1799,
                y_start: 800,
                y_end: 2399,
                y_strand: 0,
                gaps: vec![(200, 0), (1799, 0)],
                windows: exact_windows([4, 4, 4, 5]),
                specs: clone_specs(&specs),
            },
            Scenario {
                name: "C_ins6",
                reads: vec![query.clone(), inserted],
                x_start: 200,
                x_end: 1799,
                y_start: 800,
                y_end: 2405,
                y_strand: 0,
                gaps: vec![(200, 0), (700, 6), (1799, 6)],
                windows: vec![
                    WindowSpec { x_start: 200, x_end: 699, y_start: 800, y_end: 1299, error: 1 },
                    WindowSpec { x_start: 700, x_end: 1099, y_start: 1306, y_end: 1705, error: 1 },
                    WindowSpec { x_start: 1100, x_end: 1399, y_start: 1706, y_end: 2005, error: 1 },
                    WindowSpec { x_start: 1400, x_end: 1799, y_start: 2006, y_end: 2405, error: 1 },
                ],
                specs: clone_specs(&specs),
            },
            Scenario {
                name: "D_reverse",
                reads: vec![query.clone(), reversed],
                x_start: 200,
                x_end: 1799,
                y_start: 800,
                y_end: 2399,
                y_strand: 1,
                gaps: vec![(200, 0), (1799, 0)],
                windows: exact_windows([0, 0, 0, 0]),
                specs: clone_specs(&specs),
            },
        ]
    }

    fn clone_specs(specs: &[DriverSpec]) -> Vec<DriverSpec> {
        specs
            .iter()
            .map(|spec| DriverSpec {
                query_start: spec.query_start,
                query_end: spec.query_end,
                target_start: spec.target_start,
                target_end: spec.target_end,
                mode: spec.mode,
                maximum_length: spec.maximum_length,
                maximum_error: spec.maximum_error,
                forced_length: spec.forced_length,
                estimate_error: spec.estimate_error,
                use_hc: spec.use_hc,
            })
            .collect()
    }

    #[test]
    #[ignore = "requires reference/ec_alnadv_probe.out from the C++ probe"]
    fn alignment_driver_matches_the_reference() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("reference")
            .join("ec_alnadv_probe.out");
        let contents = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("cannot read {}: {}", path.display(), error));
        let golden: Vec<&str> = contents.lines().collect();
        let mut generated = Vec::new();
        for scenario in scenarios() {
            generated.extend(run_scenario(&scenario));
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

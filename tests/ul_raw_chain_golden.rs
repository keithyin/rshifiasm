//! Field-exact replay of `gen_rid_raw_chain` for read 0 of the reference UL
//! run. The DP input (192 raw sub-overlaps), the trace table (234 entries)
//! and the C++ output (one chain head with score 44055 spanning the read and
//! 22 bodies) were dumped by `/tmp/ul_worker_probe.cpp` linked against the
//! unmodified hifiasm objects.
//!
//! The `sec`/`el`/`tn` field conventions of the fixture follow the C++ dump:
//! `sec` is the DP's per-entry error (`SEC_MODE` marks none), `el` is the
//! contained flag, and the replay mirrors what `ul_rid_lalign_adv` produced.

use rshifiasm::assembly::ul::{UlOv, MAX_SEC_ERR};
use rshifiasm::assembly::ul_raw_chain::{
    gen_rid_raw_chain, GlChain, UlGraphOptions,
};
use rshifiasm::overlap::paf::PafHit;
use rshifiasm::correction::hc_align::BitExtz;
use rshifiasm::correction::target::TargetSource;
use rshifiasm::correction::ul_align::Rtrace;
use rshifiasm::overlap::region::ChainData;

mod input {
    include!("ul_raw_chain_input.rs");
}

mod regions {
    include!("ul_raw_chain_regions.rs");
}

fn build_regions() -> Vec<rshifiasm::overlap::region::OverlapRegion> {
    regions::REGIONS
        .iter()
        .map(|fixture| {
            let mut region = rshifiasm::overlap::region::OverlapRegion {
                x_identifier: fixture.x_id,
                x_position_start: fixture.x_s,
                x_position_end: fixture.x_e,
                y_identifier: fixture.y_id,
                y_position_start: fixture.y_s,
                y_position_end: fixture.y_e,
                y_strand: fixture.strand,
                alignment_length: fixture.align_length,
                non_homopolymer_errors: fixture.non_homopolymer_errors,
                ..rshifiasm::overlap::region::OverlapRegion::default()
            };
            for &(pos, shift) in fixture.fake_cigar {
                region.fake_cigar.push(pos, shift);
            }
            for &(x_s, x_e, y_s, y_e, eb, ee, err, thre, cidx, clen) in fixture.windows {
                region.window_cigar.windows.push(rshifiasm::correction::window::WindowList {
                    x_start: x_s,
                    x_end: x_e,
                    y_start: y_s,
                    y_end: y_e,
                    extra_begin: eb,
                    extra_end: ee,
                    error: err,
                    error_threshold: thre,
                    cigar_start: cidx,
                    cigar_length: clen,
                });
            }
            region.window_cigar.cells = fixture.cells.to_vec();
            region
        })
        .collect()
}

fn build_entries() -> Vec<UlOv> {
    input::DP_INPUT
        .iter()
        .map(|&(qn, tn, qs, qe, ts, te, sec, el, rev)| UlOv {
            qn,
            qs,
            qe,
            tn,
            ts,
            te,
            sec: if sec == 0x3fffffff as u32 || sec >= MAX_SEC_ERR {
                sec
            } else {
                sec
            },
            el: el != 0,
            rev: rev != 0,
        })
        .collect()
}

fn build_trace() -> Vec<Rtrace> {
    input::TRACE_INPUT
        .iter()
        .map(|&(oid, c_qs, c_qe, c_ts, c_te, c_wsid, c_weid, c_wsii, c_weii, pfx, sfx, mid)| {
            Rtrace {
                oid,
                c_qs,
                c_qe,
                c_ts,
                c_te,
                c_wsid,
                c_weid,
                c_wsii,
                c_weii,
                pfx_e: pfx,
                sfx_e: sfx,
                mid_e: mid,
                ..Rtrace::default()
            }
        })
        .collect()
}

fn parse_golden() -> (Vec<(u32, u32, u32, u32, u32, u32, u32, u32, u32, u32)>, String) {
    let text = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/read0_golden.txt"),
    )
    .expect("tests/read0_golden.txt must exist");
    let mut heads = Vec::new();
    let mut bodies = Vec::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        let get = |key: &str| -> u32 {
            parts
                .iter()
                .find(|part| part.starts_with(key) && part.len() > key.len())
                .map(|part| part[key.len()..].parse().unwrap())
                .unwrap()
        };
        let row = (
            get("qn="),
            get("tn="),
            get("qs="),
            get("qe="),
            get("ts="),
            get("te="),
            get("sec="),
            get("el="),
            get("rev="),
            0,
        );
        if parts[0] == "HEAD" {
            heads.push(row);
        } else {
            bodies.push(row);
        }
    }
    (bodies, format!("heads={heads:?}"))
}

#[test]
#[ignore]
fn raw_chain_replays_the_reference_read_zero() {
    // The `sources`/read lengths feed `get_ecov_contain_adv`; the fixture's
    // entries reference HiFi reads 1810-1990, so the real reference caches
    // are required.
    let mut options = rshifiasm::command_line::options::AssemblyOptions::default();
    options.output_file_name = "/tmp/ulref".to_string();
    let (reads, state) =
        rshifiasm::assembly::read_cache::load_all_data_from_disk(&mut options)
            .expect("the reference caches must load");
    let sources: Vec<Vec<PafHit>> = state
        .forward
        .iter()
        .map(|list| list.hits.clone())
        .collect();
    let graph_options = UlGraphOptions {
        sources: &sources,
        minimum_overlap: options.minimum_overlap_length,
        maximum_hang: options.maximum_hang_length,
        maximum_hang_rate: options.maximum_hang_rate,
    };

    let entries = build_entries();
    let trace = build_trace();
    let mut chain = GlChain {
        lo: Vec::new(),
        tk: entries,
        tc: trace,
        srt: Vec::new(),
    };
    let mut chain_data = ChainData::default();
    let mut extz = BitExtz::new();
    let regions = build_regions();
    let targets = TargetSource::Reads(&reads);
    let _ = targets;
    let unitigs: Vec<rshifiasm::assembly::asg::MaUtg> = Vec::new();
    let _ = &unitigs;

    let query_length = 59620;
    let fully_covered = gen_rid_raw_chain(
        &mut chain,
        0,
        &mut chain_data,
        &graph_options,
        &reads,
        &unitigs,
        0.2,
        query_length,
        &TargetSource::Reads(&reads),
        &[],
        &mut extz,
        0,
        180,
        &regions,
    );

    let (golden_bodies, _) = parse_golden();
    assert_eq!(chain.lo.len(), 1, "one chain head");
    assert_eq!(chain.lo[0].qn, input::HEAD_QN, "chain head score");
    assert_eq!(chain.lo[0].qs, 0, "chain head query start");
    assert_eq!(chain.lo[0].qe, query_length as u32, "chain head query end");
    assert_eq!(chain.lo[0].te, input::HEAD_TE, "chain head body count");
    assert_eq!(fully_covered, 1, "the chain covers the whole read");
    assert_eq!(
        chain.tk.len(),
        golden_bodies.len(),
        "chain body count (expected the C++ bodies)"
    );
    for (index, (expected, actual)) in golden_bodies.iter().zip(chain.tk.iter()).enumerate() {
        let actual_row = (
            actual.qn,
            actual.tn,
            actual.qs,
            actual.qe,
            actual.ts,
            actual.te,
            actual.sec,
            actual.el as u32,
            actual.rev as u32,
            0,
        );
        assert_eq!(&actual_row, expected, "chain body {}", index);
    }
}

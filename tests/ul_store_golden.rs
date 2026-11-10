//! Layout check for the `--ul` overlap cache.
//!
//! The reference run with `--ul` writes `<prefix>.ul.ovlp.bin` via
//! `write_all_ul_t`, and `load_all_ul_t` reads it back. The test loads the
//! reference dump, writes it with the Rust serializer and loads it again; a
//! mismatch means the record layout no longer matches `Process_Read.h`, and a
//! parse failure means an earlier field is misaligned.
//!
//! The dump lives under `reference/` (git-ignored); the test reports a skip
//! when it is absent.

use std::fs;
use std::path::PathBuf;

use rshifiasm::assembly::ul_store::{load_all_ul_t, write_all_ul_t};

fn reference_stem() -> PathBuf {
    std::env::var_os("RSHIFIASM_UL_REFERENCE")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("reference/ul")
        })
}

#[test]
fn the_ul_cache_round_trips_through_the_reference_layout() {
    let stem = reference_stem();
    let Some(reference) = load_all_ul_t(stem.to_str().unwrap()).unwrap() else {
        eprintln!("skip: {} is absent", stem.display());
        return;
    };
    assert!(!reference.reads.is_empty(), "empty UL cache");

    let output_directory = std::env::temp_dir().join(format!("rshifiasm_ul_store_{}", std::process::id()));
    fs::create_dir_all(&output_directory).unwrap();
    let stem_out = output_directory.join("roundtrip");
    let stem_out = stem_out.to_str().unwrap();
    write_all_ul_t(&reference, stem_out).unwrap();
    let reloaded = load_all_ul_t(stem_out).unwrap().expect("writer produced no file");
    assert_eq!(reference, reloaded);
    let _ = fs::remove_dir_all(&output_directory);
}

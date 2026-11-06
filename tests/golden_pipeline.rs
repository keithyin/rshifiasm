//! End to end consistency tests.
//!
//! The tests in this file require the reference C++ outputs produced by
//! `scripts/reference_run.sh` and are therefore ignored by default. Run them
//! explicitly after generating the reference artefacts:
//!
//! ```text
//! scripts/reference_run.sh 1
//! cargo test --release --test golden_pipeline -- --ignored
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

fn repository_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).to_path_buf()
}

fn reference_directory() -> PathBuf {
    repository_root().join("reference")
}

fn candidate_prefix(profile: &str) -> PathBuf {
    reference_directory().join(format!("candidate_{profile}"))
}

fn rust_binary() -> PathBuf {
    repository_root().join("target").join("release").join("rshifiasm")
}

/// Runs the Rust binary with the arguments of the named reference profile.
fn run_rust_assembly(profile: &str, extra_arguments: &[&str]) {
    let prefix = candidate_prefix(profile);
    let mut command = Command::new(rust_binary());
    command.args(["-t1", "-o", prefix.to_str().unwrap()]);
    command.args(extra_arguments);
    command.arg("/home/keith/Downloads/chr11-2M.fa.gz");
    let status = command
        .status()
        .expect("the Rust binary must be built with --release first");
    assert!(status.success(), "the assembly run must succeed");
}

fn assert_identical(reference_profile: &str, candidate_profile: &str, suffix: &str) {
    let reference_file =
        PathBuf::from(format!("{}{}", reference_directory().join(reference_profile).display(), suffix));
    let candidate_file =
        PathBuf::from(format!("{}{}", candidate_prefix(candidate_profile).display(), suffix));

    let reference = std::fs::read(&reference_file)
        .unwrap_or_else(|_| panic!("missing reference file {}", reference_file.display()));
    let candidate = std::fs::read(&candidate_file)
        .unwrap_or_else(|_| panic!("missing candidate file {}", candidate_file.display()));

    if reference != candidate {
        let mismatch = reference
            .iter()
            .zip(candidate.iter())
            .position(|(left, right)| left != right)
            .unwrap_or(reference.len().min(candidate.len()));
        panic!(
            "{} differs at byte {} (reference {} bytes, candidate {} bytes)",
            suffix,
            mismatch,
            reference.len(),
            candidate.len()
        );
    }
}

#[test]
#[ignore = "requires the golden reference outputs and a release build"]
fn assembly_output_is_byte_identical_to_the_reference() {
    run_rust_assembly("default", &[]);
    for suffix in [
        ".bp.p_utg.gfa",
        ".bp.p_ctg.gfa",
        ".bp.hap1.p_ctg.gfa",
        ".bp.hap2.p_ctg.gfa",
        ".bp.p_ctg.noseq.gfa",
        ".bp.p_ctg.lowQ.bed",
    ] {
        assert_identical("reference", "default", suffix);
    }
}

#[test]
#[ignore = "requires the golden reference outputs and a release build"]
fn corrected_reads_are_byte_identical_to_the_reference() {
    run_rust_assembly("write_paf_ec", &["--write-paf", "--write-ec"]);
    assert_identical("write_paf_ec", "write_paf_ec", ".ec.fa");
}

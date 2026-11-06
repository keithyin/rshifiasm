//! Replays `clean_weak_ma_hit_t` from the reference d00 state and compares the
//! result with the reference d01 dump.
//!
//! This isolates the weak-overlap cleaner from the rest of the pipeline: the
//! inputs are the reference's own dumps, so any difference comes from the Rust
//! implementation of `renew_graph_init`, `normalize_ma_hit_t_single_side_advance`
//! or `clean_weak_ma_hit_t`.

use std::path::{Path, PathBuf};

use rshifiasm::assembly::hit::{
    clean_weak_ma_hit_t, normalize_ma_hit_t_single_side_advance, renew_graph_init, ReadToUnitig,
};
use rshifiasm::assembly::probe::load_sources;

fn reference_directory() -> PathBuf {
    std::env::var_os("RSHIFIASM_CLEANPROBE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("/tmp/cleanprobe"))
}

fn file(directory: &Path, tag: &str, suffix: &str) -> PathBuf {
    directory.join(format!("{}.all.debug.{}.bin", tag, suffix))
}

fn hits(path: &Path) -> Vec<Vec<rshifiasm::overlap::paf::PafHit>> {
    load_sources(path).into_iter().map(|list| list.hits).collect()
}

#[test]
fn the_weak_cleaner_reproduces_the_reference_after_normalization() {
    let reference = reference_directory();
    let forward_path = file(&reference, "d00_after_try_rescue", "source");
    let reverse_path = file(&reference, "d00_after_try_rescue", "reverse");
    let golden_path = file(&reference, "d01_after_weak_clean", "source");
    if !forward_path.exists() || !reverse_path.exists() || !golden_path.exists() {
        eprintln!("skipping: the reference dumps are missing from {}", reference.display());
        return;
    }

    let mut forward = hits(&forward_path);
    let mut reverse = hits(&reverse_path);
    let read_count = forward.len();
    let mut read_to_unitig = ReadToUnitig::new(read_count);

    renew_graph_init(&mut forward, &mut reverse, &mut read_to_unitig);
    normalize_ma_hit_t_single_side_advance(&mut forward, false, None, 0);
    normalize_ma_hit_t_single_side_advance(&mut reverse, false, None, 0);
    clean_weak_ma_hit_t(&mut forward, &reverse, &vec![0u64; read_count], u32::MAX);

    let golden = load_sources(&golden_path);
    assert_eq!(forward.len(), golden.len(), "read count");
    let mut differing_reads = 0usize;
    let mut first_report = 0usize;
    for (read, (produced, golden)) in forward.iter().zip(&golden).enumerate() {
        if produced.len() != golden.hits.len() {
            differing_reads += 1;
            if first_report < 10 {
                eprintln!(
                    "read {}: {} hits, reference has {}",
                    read,
                    produced.len(),
                    golden.hits.len()
                );
                first_report += 1;
            }
            continue;
        }
        let mut read_differs = false;
        for (index, (left, right)) in produced.iter().zip(&golden.hits).enumerate() {
            if left != right {
                read_differs = true;
                if first_report < 10 {
                    eprintln!("read {} hit {}:\n  rust: {:?}\n  ref:  {:?}", read, index, left, right);
                    first_report += 1;
                }
            }
        }
        if read_differs {
            differing_reads += 1;
        }
    }

    assert_eq!(differing_reads, 0, "the weak cleaner changed {} reads", differing_reads);
}

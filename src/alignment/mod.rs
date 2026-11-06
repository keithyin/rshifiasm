//! Sequence alignment kernels.
//!
//! `edit_distance` holds the bit-parallel banded edit distance used by the
//! overlap and correction paths, `ksw2` the extension alignment port.

pub mod banded_bpm;
pub mod banded_path;
pub mod edit_distance;
pub mod ksw2;

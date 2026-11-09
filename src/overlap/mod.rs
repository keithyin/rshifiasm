//! Overlap detection: candidate seeds, chaining, overlap records and the
//! per-read driver.

pub mod append;
pub mod candidates;
pub mod chain_score;
pub mod chaining;
pub mod engine;
pub mod final_overlap;
pub mod lchain_parameters;
pub mod paf;
pub mod qgen_chain;
pub mod qgen_input;
pub mod qgen_seeds;
pub mod region;
pub mod stored_overlap;

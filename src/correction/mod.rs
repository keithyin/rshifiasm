//! Read correction: the correction cigar, the window cigars, the sliding window
//! bookkeeping and the partial order alignment graph today; the alignment
//! passes on top of them next.

pub mod advance;
pub mod aln_adv;
pub mod boundary;
pub mod boundary_cigar;
pub mod cal_exz;
pub mod cigar;
pub mod coverage;
pub mod consensus;
pub mod consensus_gfa;
pub mod consensus_graph;
pub mod correct;
pub mod dag_con;
pub mod dedup_chains;
pub mod ec_driver;
pub mod cigar_generator;
pub mod fast_cigar;
pub mod gap_reassign;
pub mod hc_align;
pub mod hc_r_alin;
pub mod haplotype_evidence;
pub mod interval;
pub mod non_trim;
pub mod ovlp_base;
pub mod pipeline;
pub mod poa;
pub mod recalculate;
pub mod read_writeback;
pub mod region;
pub mod rphase;
pub mod second_round;
pub mod verify;
pub mod window;
pub mod window_chain;
pub mod window_consensus;

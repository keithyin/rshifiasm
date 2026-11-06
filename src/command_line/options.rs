//! Assembly options, field for field identical to `hifiasm_opt_t`
//! (`CommandLines.h`) together with its default initialisation
//! (`init_opt` in `CommandLines.cpp`).
//!
//! Field names are spelled out instead of using the original abbreviations so
//! that the translation reads naturally, while the doc comments keep the
//! original identifier available for cross checking against the C++ code.

/// Version string reported by `--version` and at the end of every run.
pub const VERSION: &str = "0.25.0-r726";

/// Default output prefix when `-o` is not given (`DEFAULT_OUTPUT`).
pub const DEFAULT_OUTPUT_PREFIX: &str = "hifiasm.asm";

/// Minimum sequence divergence allowed for an overlap (`HA_MIN_OV_DIFF`).
pub const MINIMUM_OVERLAP_DIVERGENCE: f64 = 0.02;

/// Lower bound of the adaptive chain count (`MIN_N_CHAIN`).
pub const MINIMUM_CHAIN_COUNT: i32 = 100;

// Option bit flags, mirroring the `HA_F_*` constants.
pub const FLAG_NO_HPC: u32 = 0x1;
pub const FLAG_NO_KMER_FILTER: u32 = 0x2;
pub const FLAG_VERBOSE_GFA: u32 = 0x4;
pub const FLAG_WRITE_ERROR_CORRECTED_READS: u32 = 0x8;
pub const FLAG_WRITE_PAF: u32 = 0x10;
pub const FLAG_SKIP_TRIOBIN: u32 = 0x20;
pub const FLAG_PURGE_CONTAINMENT: u32 = 0x40;
pub const FLAG_PURGE_JOIN: u32 = 0x80;
pub const FLAG_BAN_POST_JOIN: u32 = 0x100;
pub const FLAG_BAN_ASSEMBLY: u32 = 0x200;
pub const FLAG_HIGH_HETEROZYGOSITY: u32 = 0x400;
pub const FLAG_PARTITION: u32 = 0x800;
pub const FLAG_FAST: u32 = 0x1000;
pub const FLAG_UNSKEW: u32 = 0x2000;

/// Trio related read labels (`AMBIGU`, `FATHER`, `MOTHER`, ...).
pub const TRIO_UNKNOWN: u8 = 0;
pub const TRIO_FATHER: u8 = 1;
pub const TRIO_MOTHER: u8 = 2;

/// A comma separated list of file names, the equivalent of the C `enzyme`
/// structure. It is reused for read lists, Hi-C read pairs and restriction
/// enzyme motifs, exactly as in the original code.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FileList {
    /// Original names, without the separating commas.
    pub file_names: Vec<String>,
    /// Length of every name, kept for parity with `enzyme::l`.
    pub name_lengths: Vec<usize>,
}

impl FileList {
    /// Splits a comma separated argument, the equivalent of `get_hic_enzymes`.
    ///
    /// When `validate_as_sequence` is set, the argument must consist of
    /// A/C/G/T/N characters only; otherwise an empty list is returned, which
    /// the option checker later reports as a malformed argument.
    pub fn from_comma_separated(argument: &str, validate_as_sequence: bool) -> Self {
        if argument.is_empty() {
            return Self::default();
        }

        if validate_as_sequence {
            // The separators themselves are skipped by the original loop, so
            // only the characters inside each entry are validated.
            let has_invalid_character = argument.chars().any(|character| {
                character != ','
                    && !matches!(
                        character,
                        'A' | 'C' | 'G' | 'T' | 'a' | 'c' | 'g' | 't' | 'N' | 'n'
                    )
            });
            if has_invalid_character {
                return Self::default();
            }
        }

        let mut file_list = Self::default();
        for name in argument.split(',') {
            file_list.name_lengths.push(name.len());
            file_list.file_names.push(name.to_string());
        }
        file_list
    }

    /// Number of entries, the equivalent of `enzyme::n`.
    pub fn len(&self) -> usize {
        self.file_names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.file_names.is_empty()
    }
}

/// All tunables of a run; the equivalent of `hifiasm_opt_t`.
#[derive(Debug, Clone)]
pub struct AssemblyOptions {
    pub flags: u32,
    /// `num_reads`: number of input read files.
    pub input_file_count: i32,
    /// `read_file_names`: input read files in command line order.
    pub input_file_names: Vec<String>,
    /// `output_file_name`: prefix shared by every output file.
    pub output_file_name: String,
    /// `required_read_name`: `-b`, restrains the assembly to a read list.
    pub required_read_name: Option<String>,
    /// `fn_bin_yak`: `-1`/`-2`, paternal and maternal yak k-mer dumps.
    pub trio_yak_file_names: [Option<String>; 2],
    /// `fn_bin_list`: `-3`/`-4`, paternal and maternal read name lists.
    pub trio_list_file_names: [Option<String>; 2],
    /// `fn_bin_poy`: `-5`, list of trio binned reads.
    pub trio_polyploid_file_name: Option<String>,
    /// `extract_list`: `--ex-list`, dump reads selected in a given round.
    pub extract_list_file_name: Option<String>,
    /// `hic_reads`: `--h1`/`--h2`, the two Hi-C read files.
    pub hic_read_files: [Option<FileList>; 2],
    /// `hic_enzymes`: `--enzyme`, restriction enzyme motifs.
    pub hic_enzyme_motifs: Option<FileList>,
    /// `ar`: `--ul`, ultra long read files.
    pub ultra_long_read_files: Option<FileList>,
    /// `sec_in`: `--sec-in`, secondary input files.
    pub secondary_input_files: Option<FileList>,
    /// `extract_iter`: correction round used by `--ex-list`.
    pub extract_iteration: i32,
    /// `thread_num`
    pub thread_count: i32,
    /// `k_mer_length`
    pub k_mer_length: i32,
    /// `hic_mer_length`
    pub hic_k_mer_length: i32,
    /// `ul_mer_length`
    pub ultra_long_k_mer_length: i32,
    /// `trans_mer_length`
    pub transition_k_mer_length: i32,
    /// `bub_mer_length`
    pub bubble_k_mer_length: i32,
    /// `mz_win`: minimizer window size.
    pub minimizer_window_size: i32,
    /// `ul_mz_win`
    pub ultra_long_minimizer_window_size: i32,
    /// `trans_win`
    pub transition_window_size: i32,
    /// `mz_rewin`
    pub minimizer_rewindow_size: i32,
    /// `ul_mz_rewin`
    pub ultra_long_minimizer_rewindow_size: i32,
    /// `mz_sample_dist`
    pub minimizer_sample_distance: i32,
    /// `bf_shift`: number of bits of the bloom filter, 0 disables it.
    pub bloom_filter_shift: i32,
    /// `max_kmer_cnt`
    pub maximum_kmer_occurrence: i64,
    /// `high_factor`: k-mers above this multiple of the coverage are dropped.
    pub high_occurrence_factor: f64,
    /// `max_ov_diff_ec`
    pub maximum_overlap_divergence_error_correction: f64,
    /// `max_ov_diff_final`
    pub maximum_overlap_divergence_final: f64,
    /// `hom_cov`
    pub homozygous_coverage: i32,
    /// `het_cov`
    pub heterozygous_coverage: i32,
    /// `b_low_cov`
    pub breakpoint_low_coverage: i32,
    /// `b_high_cov`
    pub breakpoint_high_coverage: i32,
    /// `m_rate`
    pub exact_match_rate: f64,
    /// `max_n_chain`
    pub maximum_chain_count: i32,
    /// `min_hist_kmer_cnt`
    pub minimum_histogram_kmer_count: i32,
    /// `load_index_from_disk`
    pub load_index_from_disk: i32,
    /// `write_index_to_disk`
    pub write_index_to_disk: i32,
    /// `number_of_round`
    pub correction_round_count: i32,
    /// `number_of_pround`
    pub partial_correction_round_count: i32,
    /// `adapterLen`
    pub adapter_length: i32,
    /// `clean_round`
    pub cleaning_round_count: i32,
    /// `roundID`
    pub round_identifier: i32,
    /// `max_hang_Len`
    pub maximum_hang_length: i32,
    /// `gap_fuzz`
    pub gap_fuzz: i32,
    /// `min_overlap_Len`
    pub minimum_overlap_length: i32,
    /// `min_overlap_coverage`
    pub minimum_overlap_coverage: i32,
    /// `max_short_tip`
    pub maximum_short_tip_read_count: i32,
    /// `max_short_ul_tip`
    pub maximum_short_ultra_long_tip_read_count: i32,
    /// `max_contig_tip`
    pub maximum_contig_tip_read_count: i32,
    /// `min_cnt`
    pub minimum_binned_kmer_count: i32,
    /// `mid_cnt`
    pub middle_binned_kmer_count: i32,
    /// `purge_level_primary`
    pub purge_level_primary: i32,
    /// `purge_level_trio`
    pub purge_level_trio: i32,
    /// `purge_overlap_len`
    pub purge_overlap_length: i32,
    /// `recover_atg_cov_min`
    pub recover_alternate_contig_coverage_minimum: i64,
    /// `recover_atg_cov_max`
    pub recover_alternate_contig_coverage_maximum: i64,
    /// `hom_global_coverage`
    pub global_homozygous_coverage: i32,
    /// `hom_global_coverage_set`
    pub global_homozygous_coverage_is_set: i32,
    /// `pur_global_coverage`
    pub global_purge_coverage: i32,
    /// `bed_inconsist_rate`
    pub bed_inconsistency_rate: i32,
    /// `hic_inconsist_rate`
    pub hic_inconsistency_rate: i32,
    /// `max_hang_rate`
    pub maximum_hang_rate: f32,
    /// `min_drop_rate`
    pub minimum_drop_rate: f32,
    /// `max_drop_rate`
    pub maximum_drop_rate: f32,
    /// `purge_simi_rate_l2`
    pub purge_similarity_rate_level_two: f32,
    /// `purge_simi_rate_l3`
    pub purge_similarity_rate_level_three: f32,
    /// `purge_simi_thres`
    pub purge_similarity_threshold: f32,
    /// `trans_base_rate`
    pub transition_base_rate: f32,
    /// `trans_base_rate_sec`
    pub secondary_transition_base_rate: f32,
    /// `min_path_drop_rate`
    pub minimum_path_drop_rate: f32,
    /// `max_path_drop_rate`
    pub maximum_path_drop_rate: f32,
    /// `small_pop_bubble_size`
    pub small_bubble_population_size: i64,
    /// `large_pop_bubble_size`
    pub large_bubble_population_size: i64,
    /// `num_bases`
    pub total_bases: i64,
    /// `num_corrected_bases`
    pub corrected_bases: i64,
    /// `num_recorrected_bases`
    pub recorrected_bases: i64,
    /// `mem_buf`
    pub buffer_memory: i64,
    /// `coverage`
    pub coverage: i64,
    /// `hap_occ`
    pub haplotype_occurrence: i32,
    /// `polyploidy`
    pub polyploidy: i32,
    /// `trio_flag_occ_thres`
    pub trio_flag_occurrence_threshold: i32,
    /// `seed`
    pub random_seed: u64,
    /// `n_perturb`
    pub perturbation_round_count: i32,
    /// `f_perturb`
    pub perturbation_fraction: f64,
    /// `n_weight`
    pub reweighting_round_count: i32,
    /// `is_alt`
    pub alternate_assembly_only: u32,
    /// `misjoin_len`
    pub misjoin_length: u64,
    /// `scffold`
    pub scaffolding_length: u64,
    /// `dp_min_len`
    pub dynamic_programming_minimum_length: i32,
    /// `dp_e`
    pub dynamic_programming_error_rate: f32,
    /// `hg_size`
    pub haploid_genome_size: i64,
    /// `kpt_rate`
    pub kmer_peak_rate: f64,
    /// `infor_cov`
    pub informative_coverage: i64,
    /// `s_hap_cov`
    pub single_haplotype_coverage: i64,
    /// `trio_cov_het_ovlp`
    pub trio_heterozygous_overlap_coverage: i64,
    /// `ul_error_rate`
    pub ultra_long_error_rate: f64,
    /// `ul_error_rate_low`
    pub ultra_long_error_rate_low: f64,
    /// `ul_error_rate_hpc`
    pub ultra_long_error_rate_hpc: f64,
    /// `ul_ec_round`
    pub ultra_long_correction_round_count: i32,
    /// `ul_mod`
    pub ultra_long_mode: i32,
    /// `is_dbg_het_cnt`
    pub debug_heterozygosity_count: u8,
    /// `is_low_het_ul`
    pub low_heterozygosity_ultra_long: u8,
    /// `is_base_trans`
    pub base_level_transition: u8,
    /// `is_read_trans`
    pub read_level_transition: u8,
    /// `is_topo_trans`
    pub topology_level_transition: u8,
    /// `is_bub_trans`
    pub bubble_level_transition: u8,
    /// `bin_only`
    pub binning_only: u8,
    /// `ul_clean_round`
    pub ultra_long_cleaning_round_count: i32,
    /// `prt_dbg_gfa`
    pub print_debug_gfa: i32,
    /// `integer_correct_round`
    pub integer_correction_round: i32,
    /// `dbg_ovec_cal`
    pub debug_overlap_calculation: u8,
    /// `hifi_pst_join`
    pub hifi_post_join: u8,
    /// `ul_pst_join`
    pub ultra_long_post_join: u8,
    /// `ul_min_base`
    pub ultra_long_minimum_base_count: u32,
    /// `self_scaf`
    pub self_scaffolding: u8,
    /// `self_scaf_min`
    pub self_scaffolding_minimum_length: u64,
    /// `self_scaf_reliable_min`
    pub self_scaffolding_reliable_minimum_length: u64,
    /// `self_scaf_gap_max`
    pub self_scaffolding_maximum_gap_length: i64,
    /// `somatic_cov`
    pub somatic_coverage: i64,
    /// `telo_motif`
    pub telomere_motif: Option<String>,
    /// `telo_pen`
    pub telomere_penalty: i64,
    /// `telo_drop`
    pub telomere_maximum_drop: i64,
    /// `telo_mic_sc`
    pub telomere_minimum_score: i64,
    /// `is_ont`
    pub ont_mode: u64,
    /// `is_sc`
    pub single_cell_mode: u64,
    /// `chemical_cov`
    pub chemical_coverage: u64,
    /// `chemical_flank`
    pub chemical_flank_length: u64,
    /// `rl_cut`
    pub read_length_cutoff: i64,
    /// `sc_cut`
    pub base_quality_cutoff: i64,
}

impl Default for AssemblyOptions {
    /// `init_opt`: every default value the original binary starts from.
    fn default() -> Self {
        Self {
            flags: FLAG_PARTITION,
            input_file_count: 0,
            input_file_names: Vec::new(),
            output_file_name: DEFAULT_OUTPUT_PREFIX.to_string(),
            required_read_name: None,
            trio_yak_file_names: [None, None],
            trio_list_file_names: [None, None],
            trio_polyploid_file_name: None,
            extract_list_file_name: None,
            hic_read_files: [None, None],
            hic_enzyme_motifs: None,
            ultra_long_read_files: None,
            secondary_input_files: None,
            extract_iteration: 0,
            thread_count: 1,
            k_mer_length: 51,
            hic_k_mer_length: 31,
            ultra_long_k_mer_length: 19,
            transition_k_mer_length: 31,
            bubble_k_mer_length: 1_000_000,
            minimizer_window_size: 51,
            ultra_long_minimizer_window_size: 19,
            transition_window_size: 31,
            minimizer_rewindow_size: 1000,
            ultra_long_minimizer_rewindow_size: 360,
            minimizer_sample_distance: 500,
            bloom_filter_shift: 37,
            maximum_kmer_occurrence: 2000,
            high_occurrence_factor: 5.0,
            maximum_overlap_divergence_error_correction: 0.04,
            maximum_overlap_divergence_final: 0.03,
            homozygous_coverage: 20,
            heterozygous_coverage: -1024,
            breakpoint_low_coverage: 0,
            breakpoint_high_coverage: -1,
            exact_match_rate: 0.75,
            maximum_chain_count: MINIMUM_CHAIN_COUNT,
            minimum_histogram_kmer_count: 5,
            load_index_from_disk: 1,
            write_index_to_disk: 1,
            correction_round_count: 3,
            partial_correction_round_count: 0,
            adapter_length: 0,
            cleaning_round_count: 4,
            round_identifier: 0,
            maximum_hang_length: 1000,
            gap_fuzz: 1000,
            minimum_overlap_length: 50,
            minimum_overlap_coverage: 0,
            maximum_short_tip_read_count: 3,
            maximum_short_ultra_long_tip_read_count: 6,
            maximum_contig_tip_read_count: 3,
            minimum_binned_kmer_count: 2,
            middle_binned_kmer_count: 5,
            purge_level_primary: 3,
            purge_level_trio: 0,
            purge_overlap_length: 1,
            recover_alternate_contig_coverage_minimum: -1024,
            recover_alternate_contig_coverage_maximum: i32::MAX as i64,
            global_homozygous_coverage: -1,
            global_homozygous_coverage_is_set: 0,
            global_purge_coverage: -1,
            bed_inconsistency_rate: 70,
            hic_inconsistency_rate: 30,
            maximum_hang_rate: 0.8,
            minimum_drop_rate: 0.2,
            maximum_drop_rate: 0.8,
            purge_similarity_rate_level_two: 0.75,
            purge_similarity_rate_level_three: 0.55,
            purge_similarity_threshold: 0.0,
            transition_base_rate: 0.93,
            secondary_transition_base_rate: 0.5,
            minimum_path_drop_rate: 0.2,
            maximum_path_drop_rate: 0.6,
            small_bubble_population_size: 0,
            large_bubble_population_size: 10_000_000,
            total_bases: 0,
            corrected_bases: 0,
            recorrected_bases: 0,
            buffer_memory: 0,
            coverage: -1,
            haplotype_occurrence: 1,
            polyploidy: 2,
            trio_flag_occurrence_threshold: 60,
            random_seed: 11,
            perturbation_round_count: 10_000,
            perturbation_fraction: 0.1,
            reweighting_round_count: 3,
            alternate_assembly_only: 0,
            misjoin_length: 500_000,
            scaffolding_length: 0,
            dynamic_programming_minimum_length: 2000,
            dynamic_programming_error_rate: 0.0025,
            haploid_genome_size: -1,
            kmer_peak_rate: -1.0,
            informative_coverage: 3,
            single_haplotype_coverage: 3,
            trio_heterozygous_overlap_coverage: -1,
            ultra_long_error_rate: 0.2,
            ultra_long_error_rate_low: 0.1,
            ultra_long_error_rate_hpc: 0.2,
            ultra_long_correction_round_count: 3,
            ultra_long_mode: 0,
            debug_heterozygosity_count: 0,
            low_heterozygosity_ultra_long: 0,
            base_level_transition: 1,
            read_level_transition: 1,
            topology_level_transition: 1,
            bubble_level_transition: 1,
            binning_only: 0,
            ultra_long_cleaning_round_count: 1,
            print_debug_gfa: 0,
            integer_correction_round: 0,
            debug_overlap_calculation: 0,
            hifi_post_join: 1,
            ultra_long_post_join: 1,
            ultra_long_minimum_base_count: 0,
            self_scaffolding: 0,
            self_scaffolding_minimum_length: 250_000,
            self_scaffolding_reliable_minimum_length: 5_000_000,
            self_scaffolding_maximum_gap_length: 3_000_000,
            somatic_coverage: -1,
            telomere_motif: None,
            telomere_penalty: 1,
            telomere_maximum_drop: 2000,
            telomere_minimum_score: 500,
            ont_mode: 0,
            single_cell_mode: 0,
            chemical_coverage: 1,
            chemical_flank_length: 256,
            read_length_cutoff: 1000,
            base_quality_cutoff: 10,
        }
    }
}

impl AssemblyOptions {
    /// `ha_opt_reset_to_round`: prepares the counters of a correction round.
    pub fn reset_to_round(&mut self, round: i32) {
        self.total_bases = 0;
        self.corrected_bases = 0;
        self.recorrected_bases = 0;
        self.buffer_memory = 0;
        self.round_identifier = round;
    }

    /// `ha_opt_update_cov`: keeps the larger of the two chain count limits.
    pub fn update_coverage(&mut self, homozygous_coverage: i32) {
        let required_chain_count =
            (homozygous_coverage as f64 * self.high_occurrence_factor + 0.499) as i32;
        self.homozygous_coverage = homozygous_coverage;
        if self.maximum_chain_count < required_chain_count {
            self.maximum_chain_count = required_chain_count;
        }
        eprintln!(
            "[M::ha_opt_update_cov] updated max_n_chain to {}",
            self.maximum_chain_count
        );
    }

    /// `ha_opt_update_cov_min`: same as [`Self::update_coverage`] but the
    /// result is floored by `minimum_chain_count`.
    pub fn update_coverage_with_minimum(&mut self, homozygous_coverage: i32, minimum_chain_count: i32) {
        let required_chain_count =
            (homozygous_coverage as f64 * self.high_occurrence_factor + 0.499) as i32;
        self.homozygous_coverage = homozygous_coverage;
        self.maximum_chain_count = required_chain_count;
        if self.maximum_chain_count < minimum_chain_count {
            self.maximum_chain_count = minimum_chain_count;
        }
        eprintln!(
            "[M::ha_opt_update_cov_min] updated max_n_chain to {}",
            self.maximum_chain_count
        );
    }

    /// `ha_opt_triobin`: whether trio binning was requested.
    pub fn trio_binning_requested(&self) -> bool {
        let yak_ready = self.trio_yak_file_names[0].is_some() && self.trio_yak_file_names[1].is_some();
        let list_ready =
            self.trio_list_file_names[0].is_some() && self.trio_list_file_names[1].is_some();
        yak_ready || list_ready
    }

    /// `ha_opt_hic`: whether Hi-C reads were supplied.
    pub fn hic_partition_requested(&self) -> bool {
        self.hic_read_files[0].is_some() && self.hic_read_files[1].is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_match_the_cpp_initialisation() {
        let options = AssemblyOptions::default();
        assert_eq!(options.flags, FLAG_PARTITION);
        assert_eq!(options.coverage, -1);
        assert_eq!(options.thread_count, 1);
        assert_eq!(options.k_mer_length, 51);
        assert_eq!(options.homozygous_coverage, 20);
        assert_eq!(options.heterozygous_coverage, -1024);
        assert_eq!(options.correction_round_count, 3);
        assert_eq!(options.large_bubble_population_size, 10_000_000);
        assert_eq!(options.purge_level_primary, 3);
        assert_eq!(options.hifi_post_join, 1);
        assert_eq!(options.read_length_cutoff, 1000);
    }

    #[test]
    fn comma_separated_file_list_is_split() {
        let list = FileList::from_comma_separated("a.fq,b.fq", false);
        assert_eq!(list.file_names, vec!["a.fq", "b.fq"]);
        assert_eq!(list.name_lengths, vec![4, 4]);
    }

    #[test]
    fn enzyme_motif_validation_rejects_non_acgt() {
        assert!(FileList::from_comma_separated("A^GATCT", true).is_empty());
        assert_eq!(FileList::from_comma_separated("GATC,GANTC", true).len(), 2);
    }

    #[test]
    fn update_coverage_keeps_the_larger_chain_limit() {
        let mut options = AssemblyOptions::default();
        options.update_coverage(40);
        assert_eq!(options.maximum_chain_count, 200);
        options.update_coverage(10);
        assert_eq!(options.maximum_chain_count, 200);
    }

    #[test]
    fn update_coverage_with_minimum_respects_the_floor() {
        let mut options = AssemblyOptions::default();
        options.update_coverage_with_minimum(10, MINIMUM_CHAIN_COUNT);
        assert_eq!(options.homozygous_coverage, 10);
        assert_eq!(options.maximum_chain_count, MINIMUM_CHAIN_COUNT);
    }
}

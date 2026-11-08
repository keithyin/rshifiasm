//! Command line parsing.
//!
//! The original implementation uses `ketopt.h`, a small `getopt_long`
//! replacement with its own quirks (argument permutation, abbreviated long
//! options, `key=value` form). The parser below reproduces those quirks so that
//! every accepted invocation behaves identically.

use super::options::{AssemblyOptions, FileList, MINIMUM_OVERLAP_DIVERGENCE, VERSION};
use crate::utilities::printf_format::format_general;

/// Whether a long option takes an argument (`ko_*_argument`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArgumentRequirement {
    None,
    Required,
    Optional,
}

/// One entry of the long option table.
#[derive(Debug, Clone, Copy)]
pub struct LongOption {
    pub name: &'static str,
    pub requirement: ArgumentRequirement,
    pub value: i32,
}

pub const UNKNOWN_OPTION: i32 = '?' as i32;
pub const MISSING_ARGUMENT: i32 = ':' as i32;
pub const END_OF_OPTIONS: i32 = -1;

/// Result of a single [`Ketopt::next`] call.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsedOption {
    /// Short option character, or the `value` of the matched long option.
    pub value: i32,
    /// Option argument, when the option takes one.
    pub argument: Option<String>,
    /// Index of the argument currently being consumed; the original code
    /// refers to it as `opt.i` when it reports errors.
    pub argument_index: usize,
}

/// Parser state, the equivalent of `ketopt_t`.
#[derive(Debug, Clone)]
pub struct Ketopt {
    /// Index of the first non-option argument (`ind`).
    pub ind: usize,
    /// Last option returned (`opt`).
    pub opt: i32,
    /// Position inside the current short option cluster.
    pos: i32,
    /// Cursor over the argument vector.
    i: usize,
    /// Number of permuted non-option arguments seen so far.
    n_arguments: usize,
}

impl Default for Ketopt {
    fn default() -> Self {
        Self {
            ind: 1,
            opt: 0,
            pos: 0,
            i: 1,
            n_arguments: 0,
        }
    }
}

fn is_option_argument(argument: &str) -> bool {
    argument.len() >= 2 && argument.starts_with('-')
}

/// Moves `arguments[position]` left over `shift` elements, like
/// `ketopt_permute`.
fn permute(arguments: &mut [String], position: usize, shift: usize) {
    for step in 0..shift {
        arguments.swap(position - step, position - step - 1);
    }
}

impl Ketopt {
    /// Parses the next option, mirroring `ketopt()`.
    ///
    /// Returns `None` once every option has been consumed; `self.ind` then
    /// points at the first positional argument.
    pub fn next(
        &mut self,
        arguments: &mut [String],
        option_string: &str,
        long_options: &[LongOption],
    ) -> Option<ParsedOption> {
        let argument_count = arguments.len();

        while self.i < argument_count && !is_option_argument(&arguments[self.i]) {
            self.i += 1;
            self.n_arguments += 1;
        }

        self.opt = 0;
        let start_index = self.i;

        if self.i >= argument_count || !is_option_argument(&arguments[self.i]) {
            self.ind = self.i - self.n_arguments;
            return None;
        }

        let current = arguments[self.i].clone();
        let mut value: i32;
        let mut option_argument: Option<String> = None;

        if current.starts_with("--") {
            // Long options always consume the whole argument, which the
            // epilogue below detects through the negative position.
            self.pos = -1;
            if current.len() == 2 {
                permute(arguments, self.i, self.n_arguments);
                self.i += 1;
                self.ind = self.i - self.n_arguments;
                return None;
            }

            let name_end = current[2..]
                .find('=')
                .map(|offset| offset + 2)
                .unwrap_or(current.len());
            let requested_name = &current[2..name_end];

            let mut exact_match: Option<&LongOption> = None;
            let mut exact_match_count = 0;
            let mut partial_match: Option<&LongOption> = None;
            let mut partial_match_count = 0;
            for candidate in long_options {
                if candidate.name.starts_with(requested_name) {
                    if candidate.name.len() == requested_name.len() {
                        exact_match_count += 1;
                        exact_match = Some(candidate);
                    } else {
                        partial_match_count += 1;
                        partial_match = Some(candidate);
                    }
                }
            }

            value = UNKNOWN_OPTION;
            if exact_match_count > 1 || (exact_match_count == 0 && partial_match_count > 1) {
                return Some(ParsedOption {
                    value,
                    argument: None,
                    argument_index: self.i,
                });
            }

            let matched = if exact_match_count == 1 {
                exact_match
            } else if partial_match_count == 1 {
                partial_match
            } else {
                None
            };

            if let Some(option) = matched {
                value = option.value;
                self.opt = value;
                if current.as_bytes().get(name_end) == Some(&b'=') {
                    option_argument = Some(current[name_end + 1..].to_string());
                }
                if option.requirement == ArgumentRequirement::Required
                    && current.len() == name_end
                {
                    if self.i + 1 < argument_count {
                        self.i += 1;
                        option_argument = Some(arguments[self.i].clone());
                    } else {
                        value = MISSING_ARGUMENT;
                    }
                }
            }
        } else {
            if self.pos == 0 {
                self.pos = 1;
            }
            let current_bytes = current.as_bytes();
            let byte_position = self.pos as usize;
            value = current_bytes[byte_position] as i32;
            self.opt = value;
            self.pos += 1;

            match option_string.find(current_bytes[byte_position] as char) {
                None => value = UNKNOWN_OPTION,
                Some(matched_at) => {
                    let takes_argument = option_string.as_bytes().get(matched_at + 1) == Some(&b':');
                    if takes_argument {
                        if current.as_bytes().get(self.pos as usize).is_none() {
                            if self.i + 1 < argument_count {
                                self.i += 1;
                                option_argument = Some(arguments[self.i].clone());
                            } else {
                                value = MISSING_ARGUMENT;
                            }
                        } else {
                            option_argument = Some(current[self.pos as usize..].to_string());
                        }
                        self.pos = -1;
                    }
                }
            }
        }

        if self.pos < 0 || current.as_bytes().get(self.pos as usize).is_none() {
            self.i += 1;
            self.pos = 0;
            if self.n_arguments > 0 {
                for position in start_index..self.i {
                    permute(arguments, position, self.n_arguments);
                }
            }
        }

        self.ind = self.i - self.n_arguments;
        Some(ParsedOption {
            value,
            argument: option_argument,
            argument_index: self.i.saturating_sub(1),
        })
    }
}

/// `atoi`, including its "stop at the first non digit" behaviour.
pub fn parse_integer(argument: &str) -> i32 {
    parse_signed_prefix(argument) as i32
}

/// `atoll`, including its "stop at the first non digit" behaviour.
pub fn parse_long_integer(argument: &str) -> i64 {
    parse_signed_prefix(argument)
}

fn parse_signed_prefix(argument: &str) -> i64 {
    let trimmed = argument.trim_start();
    let mut characters = trimmed.char_indices().peekable();

    let mut sign: i64 = 1;
    if let Some((_, first)) = characters.peek().copied() {
        if first == '+' || first == '-' {
            if first == '-' {
                sign = -1;
            }
            characters.next();
        }
    }

    let mut value: i64 = 0;
    for (_, character) in characters {
        match character.to_digit(10) {
            Some(digit) => value = value.saturating_mul(10).saturating_add(digit as i64),
            None => break,
        }
    }

    sign * value
}

/// `atof`, including its "stop at the first invalid character" behaviour.
pub fn parse_double(argument: &str) -> f64 {
    let trimmed = argument.trim_start();
    let mut end = 0;
    let bytes = trimmed.as_bytes();
    let mut seen_digit = false;
    let mut seen_dot = false;
    let mut seen_exponent = false;

    while end < bytes.len() {
        let character = bytes[end];
        let valid = match character {
            b'+' | b'-' => end == 0 || bytes[end - 1] == b'e' || bytes[end - 1] == b'E',
            b'.' => !seen_dot && !seen_exponent,
            b'e' | b'E' => seen_digit && !seen_exponent,
            b'0'..=b'9' => true,
            _ => false,
        };
        if !valid {
            break;
        }
        match character {
            b'.' => seen_dot = true,
            b'e' | b'E' => seen_exponent = true,
            b'0'..=b'9' => seen_digit = true,
            _ => {}
        }
        end += 1;
    }

    trimmed[..end].parse().unwrap_or(0.0)
}

/// `inter_gsize`: parses a genome size with an optional `k`/`m`/`g` suffix.
pub fn parse_genome_size(argument: &str) -> i64 {
    if argument.len() <= 1 {
        return -2;
    }

    let last_character = argument.as_bytes()[argument.len() - 1] as char;
    let multiplier = match last_character {
        'k' | 'K' => 1_000.0,
        'm' | 'M' => 1_000_000.0,
        'g' | 'G' => 1_000_000_000.0,
        _ => return -2,
    };

    let numeric_part = &argument[..argument.len() - 1];
    (parse_double(numeric_part) * multiplier) as i64
}

/// The long option table, in the same order as `long_options` in
/// `CommandLines.cpp`.
pub const LONG_OPTIONS: &[LongOption] = &[
    LongOption { name: "version", requirement: ArgumentRequirement::None, value: 300 },
    LongOption { name: "dbg-gfa", requirement: ArgumentRequirement::None, value: 301 },
    LongOption { name: "write-paf", requirement: ArgumentRequirement::None, value: 302 },
    LongOption { name: "write-ec", requirement: ArgumentRequirement::None, value: 303 },
    LongOption { name: "skip-triobin", requirement: ArgumentRequirement::None, value: 304 },
    LongOption { name: "max-od-ec", requirement: ArgumentRequirement::Required, value: 305 },
    LongOption { name: "max-od-final", requirement: ArgumentRequirement::Required, value: 306 },
    LongOption { name: "ex-list", requirement: ArgumentRequirement::Required, value: 307 },
    LongOption { name: "ex-iter", requirement: ArgumentRequirement::Required, value: 308 },
    LongOption { name: "hom-cov", requirement: ArgumentRequirement::Required, value: 309 },
    LongOption { name: "pri-range", requirement: ArgumentRequirement::Required, value: 310 },
    LongOption { name: "lowQ", requirement: ArgumentRequirement::Required, value: 312 },
    LongOption { name: "min-hist-cnt", requirement: ArgumentRequirement::Required, value: 313 },
    LongOption { name: "h1", requirement: ArgumentRequirement::Required, value: 314 },
    LongOption { name: "h2", requirement: ArgumentRequirement::Required, value: 315 },
    LongOption { name: "enzyme", requirement: ArgumentRequirement::Required, value: 316 },
    LongOption { name: "b-cov", requirement: ArgumentRequirement::Required, value: 317 },
    LongOption { name: "h-cov", requirement: ArgumentRequirement::Required, value: 318 },
    LongOption { name: "m-rate", requirement: ArgumentRequirement::Required, value: 319 },
    LongOption { name: "primary", requirement: ArgumentRequirement::None, value: 320 },
    LongOption { name: "t-occ", requirement: ArgumentRequirement::Required, value: 321 },
    LongOption { name: "seed", requirement: ArgumentRequirement::Required, value: 322 },
    LongOption { name: "n-perturb", requirement: ArgumentRequirement::Required, value: 323 },
    LongOption { name: "f-perturb", requirement: ArgumentRequirement::Required, value: 324 },
    LongOption { name: "n-hap", requirement: ArgumentRequirement::Required, value: 325 },
    LongOption { name: "n-weight", requirement: ArgumentRequirement::Required, value: 326 },
    LongOption { name: "l-msjoin", requirement: ArgumentRequirement::Required, value: 327 },
    LongOption { name: "purge-max", requirement: ArgumentRequirement::Required, value: 328 },
    LongOption { name: "fast", requirement: ArgumentRequirement::None, value: 329 },
    LongOption { name: "dp-er", requirement: ArgumentRequirement::Required, value: 330 },
    LongOption { name: "max-kocc", requirement: ArgumentRequirement::Required, value: 331 },
    LongOption { name: "hg-size", requirement: ArgumentRequirement::Required, value: 332 },
    LongOption { name: "ul", requirement: ArgumentRequirement::Required, value: 333 },
    LongOption { name: "unskew", requirement: ArgumentRequirement::None, value: 334 },
    LongOption { name: "kpt-rate", requirement: ArgumentRequirement::Required, value: 335 },
    LongOption { name: "ul-rate", requirement: ArgumentRequirement::Required, value: 336 },
    LongOption { name: "dbg-het-cnt", requirement: ArgumentRequirement::None, value: 337 },
    LongOption { name: "ul-tip", requirement: ArgumentRequirement::Required, value: 338 },
    LongOption { name: "low-het", requirement: ArgumentRequirement::None, value: 339 },
    LongOption { name: "s-base", requirement: ArgumentRequirement::Required, value: 340 },
    LongOption { name: "bin-only", requirement: ArgumentRequirement::None, value: 341 },
    LongOption { name: "ul-round", requirement: ArgumentRequirement::Required, value: 342 },
    LongOption { name: "prt-raw", requirement: ArgumentRequirement::None, value: 343 },
    LongOption { name: "integer-correct", requirement: ArgumentRequirement::Required, value: 344 },
    LongOption { name: "dbg-ovec", requirement: ArgumentRequirement::None, value: 345 },
    LongOption { name: "path-max", requirement: ArgumentRequirement::Required, value: 346 },
    LongOption { name: "path-min", requirement: ArgumentRequirement::Required, value: 347 },
    LongOption { name: "trio-dual", requirement: ArgumentRequirement::None, value: 348 },
    LongOption { name: "ul-cut", requirement: ArgumentRequirement::Required, value: 349 },
    LongOption { name: "dual-scaf", requirement: ArgumentRequirement::None, value: 350 },
    LongOption { name: "scaf-gap", requirement: ArgumentRequirement::Required, value: 351 },
    LongOption { name: "sec-in", requirement: ArgumentRequirement::Required, value: 352 },
    LongOption { name: "somatic-cov", requirement: ArgumentRequirement::Required, value: 353 },
    LongOption { name: "telo-m", requirement: ArgumentRequirement::Required, value: 354 },
    LongOption { name: "telo-p", requirement: ArgumentRequirement::Required, value: 355 },
    LongOption { name: "telo-d", requirement: ArgumentRequirement::Required, value: 356 },
    LongOption { name: "telo-s", requirement: ArgumentRequirement::Required, value: 357 },
    LongOption { name: "ctg-n", requirement: ArgumentRequirement::Required, value: 358 },
    LongOption { name: "ont", requirement: ArgumentRequirement::None, value: 359 },
    LongOption { name: "chem-c", requirement: ArgumentRequirement::Required, value: 361 },
    LongOption { name: "chem-f", requirement: ArgumentRequirement::Required, value: 362 },
    LongOption { name: "ul-m", requirement: ArgumentRequirement::Required, value: 363 },
    LongOption { name: "rl-cut", requirement: ArgumentRequirement::Required, value: 364 },
    LongOption { name: "sc-cut", requirement: ArgumentRequirement::Required, value: 365 },
];

/// Short options understood by the parser, mirroring the `ostr` argument of
/// `ketopt` in `CommandLine_process`.
const SHORT_OPTIONS: &str = "hvt:o:k:w:m:n:r:a:b:z:x:y:p:c:d:M:P:if:D:FN:1:2:3:4:5:l:s:O:eu:";

/// Parses every option into `options`.
///
/// The error messages match the original text so that a diff of the standard
/// error stream stays meaningful.
pub fn parse_arguments(arguments: &[String], options: &mut AssemblyOptions) -> bool {
    let mut ordered_arguments: Vec<String> = arguments.to_vec();
    let mut state = Ketopt::default();

    loop {
        let parsed = match state.next(&mut ordered_arguments, SHORT_OPTIONS, LONG_OPTIONS) {
            Some(parsed) => parsed,
            None => break,
        };

        let argument = parsed.argument.clone().unwrap_or_default();

        match parsed.value {
            c if c == 'h' as i32 => {
                print_help(options);
                return false;
            }
            c if c == 'v' as i32 || c == 300 => {
                println!("{}", VERSION);
                return false;
            }
            c if c == 'f' as i32 => options.bloom_filter_shift = parse_integer(&argument),
            c if c == 't' as i32 => options.thread_count = parse_integer(&argument),
            c if c == 'o' as i32 => options.output_file_name = argument,
            c if c == 'r' as i32 => options.correction_round_count = parse_integer(&argument),
            c if c == 'k' as i32 => options.k_mer_length = parse_integer(&argument),
            c if c == 'i' as i32 => options.load_index_from_disk = 0,
            c if c == 'w' as i32 => options.minimizer_window_size = parse_integer(&argument),
            c if c == 'D' as i32 => options.high_occurrence_factor = parse_double(&argument),
            c if c == 'F' as i32 => options.flags |= super::options::FLAG_NO_KMER_FILTER,
            c if c == 'N' as i32 => options.maximum_chain_count = parse_integer(&argument),
            c if c == 'a' as i32 => options.cleaning_round_count = parse_integer(&argument),
            c if c == 'z' as i32 => options.adapter_length = parse_integer(&argument),
            c if c == 'b' as i32 => options.required_read_name = Some(argument),
            c if c == 'c' as i32 => options.minimum_binned_kmer_count = parse_integer(&argument),
            c if c == 'd' as i32 => options.middle_binned_kmer_count = parse_integer(&argument),
            c if c == '1' as i32 || c == 'P' as i32 => {
                options.trio_yak_file_names[0] = Some(argument)
            }
            c if c == '2' as i32 || c == 'M' as i32 => {
                options.trio_yak_file_names[1] = Some(argument)
            }
            c if c == '3' as i32 => options.trio_list_file_names[0] = Some(argument),
            c if c == '4' as i32 => options.trio_list_file_names[1] = Some(argument),
            c if c == '5' as i32 => options.trio_polyploid_file_name = Some(argument),
            c if c == 'x' as i32 => options.maximum_drop_rate = parse_double(&argument) as f32,
            c if c == 'y' as i32 => options.minimum_drop_rate = parse_double(&argument) as f32,
            c if c == 'p' as i32 => {
                options.small_bubble_population_size = parse_long_integer(&argument)
            }
            c if c == 'm' as i32 => {
                options.large_bubble_population_size = parse_long_integer(&argument)
            }
            c if c == 'n' as i32 => {
                options.maximum_short_tip_read_count = parse_long_integer(&argument) as i32
            }
            c if c == 'e' as i32 => options.flags |= super::options::FLAG_BAN_ASSEMBLY,
            c if c == 'u' as i32 => {
                let enabled = parse_long_integer(&argument) != 0;
                let value = u8::from(enabled);
                options.hifi_post_join = value;
                options.ultra_long_post_join = value;
            }
            c if c == 301 => options.flags |= super::options::FLAG_VERBOSE_GFA,
            c if c == 302 => options.flags |= super::options::FLAG_WRITE_PAF,
            c if c == 303 => options.flags |= super::options::FLAG_WRITE_ERROR_CORRECTED_READS,
            c if c == 304 => options.flags |= super::options::FLAG_SKIP_TRIOBIN,
            c if c == 305 => {
                options.maximum_overlap_divergence_error_correction = parse_double(&argument)
            }
            c if c == 306 => {
                options.maximum_overlap_divergence_final = parse_double(&argument)
            }
            c if c == 307 => options.extract_list_file_name = Some(argument),
            c if c == 308 => options.extract_iteration = parse_integer(&argument),
            c if c == 309 => {
                options.global_homozygous_coverage = parse_integer(&argument);
                options.global_homozygous_coverage_is_set = 1;
            }
            c if c == 310 => parse_primary_range(&argument, options),
            c if c == 312 => options.bed_inconsistency_rate = parse_integer(&argument),
            c if c == 313 => options.minimum_histogram_kmer_count = parse_integer(&argument),
            c if c == 314 => {
                options.hic_read_files[0] = Some(FileList::from_comma_separated(&argument, false))
            }
            c if c == 315 => {
                options.hic_read_files[1] = Some(FileList::from_comma_separated(&argument, false))
            }
            c if c == 316 => {
                options.hic_enzyme_motifs = Some(FileList::from_comma_separated(&argument, true))
            }
            c if c == 317 => options.breakpoint_low_coverage = parse_integer(&argument),
            c if c == 318 => options.breakpoint_high_coverage = parse_integer(&argument),
            c if c == 319 => options.exact_match_rate = parse_double(&argument),
            c if c == 320 => {
                options.flags -= super::options::FLAG_PARTITION;
                options.alternate_assembly_only = 1;
            }
            c if c == 321 => options.trio_flag_occurrence_threshold = parse_integer(&argument),
            c if c == 322 => options.random_seed = parse_long_integer(&argument) as u64,
            c if c == 323 => options.perturbation_round_count = parse_integer(&argument),
            c if c == 324 => options.perturbation_fraction = parse_double(&argument),
            c if c == 325 => options.polyploidy = parse_integer(&argument),
            c if c == 326 => options.reweighting_round_count = parse_integer(&argument),
            c if c == 327 => options.misjoin_length = parse_long_integer(&argument) as u64,
            c if c == 328 => options.global_purge_coverage = parse_integer(&argument),
            c if c == 329 => options.flags |= super::options::FLAG_FAST,
            c if c == 330 => options.dynamic_programming_error_rate = parse_double(&argument) as f32,
            c if c == 331 => options.maximum_kmer_occurrence = parse_long_integer(&argument),
            c if c == 332 => options.haploid_genome_size = parse_genome_size(&argument),
            c if c == 333 => {
                options.ultra_long_read_files = Some(FileList::from_comma_separated(&argument, false))
            }
            c if c == 334 => options.flags |= super::options::FLAG_UNSKEW,
            c if c == 335 => options.kmer_peak_rate = parse_double(&argument),
            c if c == 336 => options.ultra_long_error_rate = parse_double(&argument),
            c if c == 337 => options.debug_heterozygosity_count = 1,
            c if c == 338 => {
                options.maximum_short_ultra_long_tip_read_count =
                    parse_long_integer(&argument) as i32
            }
            c if c == 339 => options.low_heterozygosity_ultra_long = 1,
            c if c == 340 => {
                options.secondary_transition_base_rate = parse_double(&argument) as f32;
                if options.secondary_transition_base_rate < 0.0 {
                    options.base_level_transition = 0;
                }
            }
            c if c == 341 => options.binning_only = 1,
            c if c == 342 => {
                options.ultra_long_cleaning_round_count = parse_long_integer(&argument) as i32
            }
            c if c == 343 => options.print_debug_gfa = 1,
            c if c == 344 => options.integer_correction_round = parse_long_integer(&argument) as i32,
            c if c == 345 => options.debug_overlap_calculation = 1,
            c if c == 346 => options.maximum_path_drop_rate = parse_double(&argument) as f32,
            c if c == 347 => options.minimum_path_drop_rate = parse_double(&argument) as f32,
            c if c == 348 => options.trio_heterozygous_overlap_coverage = 1,
            c if c == 349 => {
                options.ultra_long_minimum_base_count = parse_long_integer(&argument) as u32
            }
            c if c == 350 => options.self_scaffolding = 1,
            c if c == 351 => options.self_scaffolding_maximum_gap_length = parse_long_integer(&argument),
            c if c == 352 => {
                options.secondary_input_files =
                    Some(FileList::from_comma_separated(&argument, false))
            }
            c if c == 353 => options.somatic_coverage = parse_long_integer(&argument),
            c if c == 354 => options.telomere_motif = Some(argument),
            c if c == 355 => options.telomere_penalty = parse_long_integer(&argument),
            c if c == 356 => options.telomere_maximum_drop = parse_long_integer(&argument),
            c if c == 357 => options.telomere_minimum_score = parse_long_integer(&argument),
            c if c == 358 => {
                options.maximum_contig_tip_read_count = parse_long_integer(&argument) as i32
            }
            c if c == 359 => {
                options.ont_mode = 1;
                options.maximum_overlap_divergence_error_correction = 0.07;
                options.single_cell_mode = 1;
            }
            c if c == 361 => options.chemical_coverage = parse_long_integer(&argument) as u64,
            c if c == 362 => options.chemical_flank_length = parse_long_integer(&argument) as u64,
            c if c == 364 => options.read_length_cutoff = parse_long_integer(&argument),
            c if c == 365 => options.base_quality_cutoff = parse_long_integer(&argument),
            c if c == 'l' as i32 => {
                let level = parse_integer(&argument);
                options.purge_level_primary = level;
                options.purge_level_trio = level;
            }
            c if c == 's' as i32 => {
                let rate = parse_double(&argument) as f32;
                options.purge_similarity_rate_level_two = rate;
                options.purge_similarity_rate_level_three = rate;
            }
            c if c == 'O' as i32 => options.purge_overlap_length = parse_long_integer(&argument) as i32,
            c if c == MISSING_ARGUMENT => {
                eprintln!(
                    "[ERROR] missing option argument in \"{}\"",
                    ordered_arguments[parsed.argument_index]
                );
                return false;
            }
            c if c == UNKNOWN_OPTION => {
                eprintln!(
                    "[ERROR] unknown option in \"{}\"",
                    ordered_arguments[parsed.argument_index]
                );
                return false;
            }
            _ => {}
        }
    }

    options.purge_similarity_threshold = if options.purge_level_primary > 2 {
        options.purge_similarity_rate_level_three
    } else {
        options.purge_similarity_rate_level_two
    };

    if state.ind >= ordered_arguments.len() {
        print_help(options);
        return false;
    }

    if !collect_input_files(&ordered_arguments[state.ind..], options) {
        return false;
    }

    let post_join_enabled = if options.ultra_long_read_files.is_some() {
        options.ultra_long_post_join
    } else {
        options.hifi_post_join
    };
    if post_join_enabled != 0 {
        options.flags &= !super::options::FLAG_BAN_POST_JOIN;
    } else {
        options.flags |= super::options::FLAG_BAN_POST_JOIN;
    }

    if options.ont_mode == 0 {
        options.read_length_cutoff = -1;
        options.base_quality_cutoff = 1;
    }

    check_options(options) && check_supported(options)
}

/// `--pri-range`: `INT1[,INT2]`, where `-1` disables the range filter.
fn parse_primary_range(argument: &str, options: &mut AssemblyOptions) {
    let mut parts = argument.split(',');
    let minimum = parts
        .next()
        .map(parse_long_integer)
        .unwrap_or(0);
    let maximum = parts.next().map(parse_long_integer);
    options.recover_alternate_contig_coverage_minimum = minimum;
    if let Some(maximum) = maximum {
        options.recover_alternate_contig_coverage_maximum = maximum;
    }
    if options.recover_alternate_contig_coverage_minimum == -1
        || options.recover_alternate_contig_coverage_maximum == -1
    {
        options.recover_alternate_contig_coverage_minimum = -1;
        options.recover_alternate_contig_coverage_maximum = -1;
    }
}

/// `get_queries`: records the positional arguments as input read files.
fn collect_input_files(remaining_arguments: &[String], options: &mut AssemblyOptions) -> bool {
    options.input_file_count = remaining_arguments.len() as i32;
    options.input_file_names = remaining_arguments.to_vec();

    for file_name in &options.input_file_names {
        match crate::sequence::fastx_reader::read_first_record(file_name) {
            Ok(record) => {
                if options.single_cell_mode != 0 && record.quality.is_empty() {
                    eprintln!(
                        "[ERROR] {} is in fasta format rather than fastq format",
                        file_name
                    );
                    options.single_cell_mode = 0;
                    return false;
                }
            }
            Err(message) => {
                eprintln!("[ERROR] Cannot find the input read file: {}", file_name);
                eprintln!("[ERROR] {}", message);
                return false;
            }
        }
    }

    true
}

/// `check_option`: validates the combination of options before the run starts.
fn check_options(options: &AssemblyOptions) -> bool {
    if options.input_file_count == 0 {
        eprintln!("[ERROR] missing input: please specify a read file");
        return false;
    }

    if options.thread_count < 1 {
        eprintln!("[ERROR] the number of threads must be > 0 (-t)");
        return false;
    }

    if options.correction_round_count < 1 {
        eprintln!("[ERROR] the number of rounds for correction must be > 0 (-r)");
        return false;
    }

    if options.cleaning_round_count < 1 {
        eprintln!("[ERROR] the number of rounds for assembly cleaning must be > 0 (-a)");
        return false;
    }

    if options.adapter_length < 0 {
        eprintln!("[ERROR] the length of removed adapters must be >= 0 (-z)");
        return false;
    }

    if options.k_mer_length >= 64 {
        eprintln!("[ERROR] the length of k_mer must be < 64 (-k)");
        return false;
    }

    if options.maximum_drop_rate < 0.0 || options.maximum_drop_rate >= 1.0 {
        eprintln!("[ERROR] max overlap drop ratio must be [0.0, 1.0) (-x)");
        return false;
    }

    if options.minimum_drop_rate < 0.0 || options.minimum_drop_rate >= 1.0 {
        eprintln!("[ERROR] min overlap drop ratio must be [0.0, 1.0) (-y)");
        return false;
    }

    if options.maximum_drop_rate <= options.minimum_drop_rate {
        eprintln!(
            "[ERROR] min overlap drop ratio must be less than max overlap drop ratio (-x/-y)"
        );
        return false;
    }

    if options.small_bubble_population_size < 0 {
        eprintln!("[ERROR] the size of popped small bubbles must be >= 0 (-p)");
        return false;
    }

    if options.large_bubble_population_size < 0 {
        eprintln!("[ERROR] the size of popped large bubbles must be >= 0 (-m)");
        return false;
    }

    if options.maximum_hang_length < 0 {
        eprintln!("[ERROR] max_hang_Len must be >= 0");
        return false;
    }

    if options.maximum_hang_rate < 0.0 {
        eprintln!("[ERROR] max_hang_rate must be >= 0");
        return false;
    }

    if options.gap_fuzz < 0 {
        eprintln!("[ERROR] gap_fuzz must be >= 0");
        return false;
    }

    if options.minimum_overlap_length < 0 {
        eprintln!("[ERROR] min_overlap_Len must be >= 0");
        return false;
    }

    if options.minimum_overlap_coverage < 0 {
        eprintln!("[ERROR] min_overlap_coverage must be >= 0");
        return false;
    }

    if options.maximum_overlap_divergence_error_correction
        < options.maximum_overlap_divergence_final
    {
        eprintln!("[ERROR] max_ov_diff_ec shouldn't be smaller than max_ov_diff_final");
        return false;
    }

    if options.maximum_overlap_divergence_error_correction < MINIMUM_OVERLAP_DIVERGENCE {
        eprintln!(
            "[ERROR] max_ov_diff_ec shouldn't be smaller than {}",
            MINIMUM_OVERLAP_DIVERGENCE
        );
        return false;
    }

    if options.maximum_short_tip_read_count < 0 {
        eprintln!("[ERROR] the length of removal tips must be >= 0 (-n)");
        return false;
    }

    if options.purge_level_primary < 0 || options.purge_level_primary > 3 {
        eprintln!("[ERROR] the level of purge-dup should be [0, 3] (-l)");
        return false;
    }

    if options.trio_binning_requested()
        && (options.purge_level_trio < 0 || options.purge_level_trio > 1)
    {
        eprintln!("[ERROR] the level of purge-dup for trio should be [0, 1] (-l)");
        return false;
    }

    if options.global_homozygous_coverage < 0 && options.global_homozygous_coverage != -1 {
        eprintln!("[ERROR] homozygous read coverage should be >= 0 (--hom-cov)");
        return false;
    }

    if options.global_purge_coverage < 0 && options.global_purge_coverage != -1 {
        eprintln!("[ERROR] purge duplication coverage threshold should be >= 0 (--purge-max)");
        return false;
    }

    if options.bed_inconsistency_rate < 0 || options.bed_inconsistency_rate > 100 {
        eprintln!("[ERROR] inconsistency rate should be [0, 100] (--lowQ)");
        return false;
    }

    if !check_optional_file(&options.trio_yak_file_names[0], "YAK1")
        || !check_optional_file(&options.trio_yak_file_names[1], "YAK2")
        || !check_optional_file(&options.trio_list_file_names[0], "LIST1")
        || !check_optional_file(&options.trio_list_file_names[1], "LIST2")
        || !check_optional_file(&options.required_read_name, "b")
    {
        return false;
    }

    for (index, option_name) in ["HIC1", "HIC2"].iter().enumerate() {
        if let Some(file_list) = &options.hic_read_files[index] {
            for file_name in &file_list.file_names {
                if !std::path::Path::new(file_name).exists() {
                    eprintln!("[ERROR] {} does not exist (-{})", file_name, option_name);
                    return false;
                }
            }
        }
    }

    if options.hic_read_files[0].is_some() && options.hic_read_files[1].is_none() {
        eprintln!("[ERROR] lack r2 of HiC reads (--h2)");
        return false;
    }

    if options.hic_read_files[1].is_some() && options.hic_read_files[0].is_none() {
        eprintln!("[ERROR] lack r1 of HiC reads (--h1)");
        return false;
    }

    if let (Some(first), Some(second)) = (&options.hic_read_files[0], &options.hic_read_files[1]) {
        if first.len() != second.len() {
            eprintln!("[ERROR] wrong r1 and r2 of HiC reads (--h1 && --h2)");
            return false;
        }
    }

    if let Some(motifs) = &options.hic_enzyme_motifs {
        if motifs.is_empty() {
            eprintln!("[ERROR] wrong HiC enzymes (--enzyme)");
            return false;
        }
    }

    if let Some(first) = &options.hic_read_files[0] {
        if first.is_empty() {
            eprintln!("[ERROR] wrong r1 of HiC reads (--h1)");
            return false;
        }
    }

    if let Some(second) = &options.hic_read_files[1] {
        if second.is_empty() {
            eprintln!("[ERROR] wrong r2 of HiC reads (--h2)");
            return false;
        }
    }

    if let Some(ultra_long_files) = &options.ultra_long_read_files {
        for file_name in &ultra_long_files.file_names {
            match crate::sequence::fastx_reader::read_first_record(file_name) {
                Ok(record) => {
                    if options.ultra_long_mode != 0 && record.quality.is_empty() {
                        eprintln!(
                            "[ERROR] {} is in fasta format rather than fastq format (--ul)",
                            file_name
                        );
                        eprintln!("[ERROR] set --ul-m 0 for fasta files");
                        return false;
                    }
                }
                Err(_) => {
                    eprintln!(
                        "[ERROR] Cannot find the input file: {} (--ul)",
                        file_name
                    );
                    return false;
                }
            }
        }
        if ultra_long_files.is_empty() {
            eprintln!("[ERROR] wrong UL reads (--ul)");
            return false;
        }
    }

    if options.breakpoint_low_coverage < 0 {
        eprintln!("[ERROR] must >= 0 (--b-cov)");
        return false;
    }

    if options.breakpoint_high_coverage != -1 && options.breakpoint_high_coverage < 0 {
        eprintln!("[ERROR] must >= 0 (--h-cov)");
        return false;
    }

    if options.exact_match_rate < 0.0 {
        eprintln!("[ERROR] must >= 0 (--m-rate)");
        return false;
    }

    if options.breakpoint_high_coverage != -1
        && options.breakpoint_high_coverage <= options.breakpoint_low_coverage
    {
        eprintln!("[ERROR] [--h-cov] must >= [--b-cov]");
        return false;
    }

    if options.purge_similarity_threshold < 0.0 {
        eprintln!("[ERROR] [-s] must >= 0");
        return false;
    }

    if options.maximum_kmer_occurrence < 0 {
        eprintln!("[ERROR] [--max-kocc] must >= 0");
        return false;
    }

    if options.haploid_genome_size < -1 {
        eprintln!("[ERROR] [--hg-size] wrong genome size");
        return false;
    }

    if options.ultra_long_mode != 0 && options.ultra_long_mode != 1 {
        eprintln!("[ERROR] must be 0 or 1 (--ul-m)");
        return false;
    }

    if let Some(motif) = &options.telomere_motif {
        if motif.len() > 32 {
            eprintln!("[ERROR] [--telo-m] must be no longer than 32");
            return false;
        }
        if !motif
            .chars()
            .all(|character| matches!(character, 'A' | 'C' | 'G' | 'T' | 'a' | 'c' | 'g' | 't'))
        {
            eprintln!("[ERROR] [--telo-m] must be A/C/G/T");
            return false;
        }
    }

    true
}

/// Rejects options that would enter pipeline branches this port has not
/// implemented yet. The reference would run those branches, so staying silent
/// would make such a run diverge without any warning. Parameter-only options
/// of the missing modes (`--ul-rate`, `--path-max`, `--chem-c`, `--kpt-rate`,
/// ...) are deliberately not rejected: without their mode switch the
/// reference ignores them too. `--low-het`, `--unskew`, `--skip-triobin`,
/// `--bin-only` and `--ul-m` are dead in the reference itself and stay
/// equally dead here.
fn check_supported(options: &AssemblyOptions) -> bool {
    if options.ultra_long_read_files.is_some() {
        eprintln!("[ERROR] the ultra-long integration (--ul) is not implemented yet");
        return false;
    }
    if options.flags & super::options::FLAG_VERBOSE_GFA != 0 {
        eprintln!(
            "[ERROR] --dbg-gfa (the .pt_flt index checkpoint and verbose GFA dump) is not implemented yet"
        );
        return false;
    }
    if options.required_read_name.is_some() {
        eprintln!("[ERROR] -b (single-read correction debugging) is not implemented yet");
        return false;
    }
    if options.ont_mode != 0 {
        eprintln!("[ERROR] the ONT correction branch is not implemented yet");
        return false;
    }
    if options.flags & super::options::FLAG_FAST != 0 {
        eprintln!("[ERROR] the fast-mode correction index is not implemented yet");
        return false;
    }
    true
}

fn check_optional_file(file_name: &Option<String>, option_name: &str) -> bool {
    match file_name {
        None => true,
        Some(name) => {
            if std::path::Path::new(name).exists() {
                true
            } else {
                eprintln!("[ERROR] {} does not exist (-{})", name, option_name);
                false
            }
        }
    }
}

/// `Print_H`: the help text, reproduced verbatim.
pub fn print_help(options: &AssemblyOptions) {
    eprintln!("Usage: hifiasm [options] <in_1.fq> <in_2.fq> <...>");
    eprintln!("Options:");
    eprintln!("  Input/Output:");
    eprintln!(
        "    -o STR       prefix of output files [{}]",
        options.output_file_name
    );
    eprintln!("    -t INT       number of threads [{}]", options.thread_count);
    eprintln!("    -h           show help information");
    eprintln!("    --version    show version number");
    eprintln!("  Preset options:");
    eprintln!("    --ont        assemble Oxford Nanopore reads");
    eprintln!("  Overlap/Error correction:");
    eprintln!(
        "    -k INT       k-mer length (must be <64) [{}]",
        options.k_mer_length
    );
    eprintln!(
        "    -w INT       minimizer window size [{}]",
        options.minimizer_window_size
    );
    eprintln!(
        "    -f INT       number of bits for bloom filter; 0 to disable [{}]",
        options.bloom_filter_shift
    );
    eprintln!(
        "    -D FLOAT     drop k-mers occurring >FLOAT*coverage times [{:.1}]",
        options.high_occurrence_factor
    );
    eprintln!(
        "    -N INT       consider up to max(-D*coverage,-N) overlaps for each oriented read [{}]",
        options.maximum_chain_count
    );
    eprintln!(
        "    -r INT       round of correction [{}]",
        options.correction_round_count
    );
    eprintln!(
        "    -z INT       length of adapters that should be removed [{}]",
        options.adapter_length
    );
    eprintln!("    --max-kocc   INT");
    eprintln!(
        "                 employ k-mers occurring <INT times to rescue repetitive overlaps [{}]",
        options.maximum_kmer_occurrence
    );
    eprintln!("    --hg-size    INT(k, m or g)");
    eprintln!("                 estimated haploid genome size used for inferring read coverage [auto]");
    eprintln!("  Assembly:");
    eprintln!(
        "    -a INT       round of assembly cleaning [{}]",
        options.cleaning_round_count
    );
    eprintln!(
        "    -m INT       pop bubbles of <INT in size in contig graphs [{}]",
        options.large_bubble_population_size
    );
    eprintln!(
        "    -p INT       pop bubbles of <INT in size in unitig graphs [{}]",
        options.small_bubble_population_size
    );
    eprintln!(
        "    -n INT       remove tip unitigs composed of <=INT reads [{}]",
        options.maximum_short_tip_read_count
    );
    eprintln!(
        "    -x FLOAT     max overlap drop ratio [{}]",
        format_general(options.maximum_drop_rate as f64, 2)
    );
    eprintln!(
        "    -y FLOAT     min overlap drop ratio [{}]",
        format_general(options.minimum_drop_rate as f64, 2)
    );
    eprintln!("    -i           ignore saved read correction and overlaps");
    eprintln!("    -u           post-join step for contigs which may improve N50; 0 to disable; 1 to enable");
    eprintln!(
        "                 [{}] and [{}] in default for the UL+HiFi assembly and the HiFi assembly, respectively",
        options.ultra_long_post_join, options.hifi_post_join
    );
    eprintln!("    --hom-cov    INT");
    eprintln!("                 homozygous read coverage [auto]");
    eprintln!("    --lowQ       INT");
    eprintln!(
        "                 output contig regions with >=INT% inconsistency in BED format; 0 to disable [{}]",
        options.bed_inconsistency_rate
    );
    eprintln!("    --b-cov      INT");
    eprintln!(
        "                 break contigs at positions with <INT-fold coverage; work with '--m-rate'; 0 to disable [{}]",
        options.breakpoint_low_coverage
    );
    eprintln!("    --h-cov      INT");
    eprintln!(
        "                 break contigs at positions with >INT-fold coverage; work with '--m-rate'; -1 to disable [{}]",
        options.breakpoint_high_coverage
    );
    eprintln!("    --m-rate     FLOAT");
    eprintln!("                 break contigs at positions with <=FLOAT*coverage exact overlaps;");
    eprintln!(
        "                 only work with '--b-cov' or '--h-cov'[{:.2}]",
        options.exact_match_rate
    );
    eprintln!("    --primary    output a primary assembly and an alternate assembly");
    eprintln!("    --ctg-n      INT");
    eprintln!(
        "                 remove tip contigs composed of <=INT reads [{}]",
        options.maximum_contig_tip_read_count
    );
    eprintln!("  Trio-partition:");
    eprintln!("    -1 FILE      hap1/paternal k-mer dump generated by \"yak count\" []");
    eprintln!("    -2 FILE      hap2/maternal k-mer dump generated by \"yak count\" []");
    eprintln!("    -3 FILE      list of hap1/paternal read names []");
    eprintln!("    -4 FILE      list of hap2/maternal read names []");
    eprintln!(
        "    -c INT       lower bound of the binned k-mer's frequency [{}]",
        options.minimum_binned_kmer_count
    );
    eprintln!(
        "    -d INT       upper bound of the binned k-mer's frequency [{}]",
        options.middle_binned_kmer_count
    );
    eprintln!("    --t-occ      INT");
    eprintln!("                 forcedly remove unitigs with >INT unexpected haplotype-specific reads;");
    eprintln!(
        "                 ignore graph topology; [{}]",
        options.trio_flag_occurrence_threshold
    );
    eprintln!("    --trio-dual  utilize homology information to correct trio phasing errors");
    eprintln!("  Purge-dups:");
    eprintln!("    -l INT       purge level. 0: no purging; 1: light; 2/3: aggressive [0 for trio; 3 for unzip]");
    eprintln!(
        "    -s FLOAT     similarity threshold for duplicate haplotigs in read-level [{} for -l1/-l2, {} for -l3]",
        options.purge_similarity_rate_level_two, options.purge_similarity_rate_level_three
    );
    eprintln!(
        "    -O INT       min number of overlapped reads for duplicate haplotigs [{}]",
        options.purge_overlap_length
    );
    eprintln!("    --purge-max  INT");
    eprintln!("                 coverage upper bound of Purge-dups [auto]");
    eprintln!("    --n-hap      INT");
    eprintln!("                 number of haplotypes [{}]", options.polyploidy);
    eprintln!("  Hi-C-partition:");
    eprintln!("    --h1 FILEs   file names of Hi-C R1  [r1_1.fq,r1_2.fq,...]");
    eprintln!("    --h2 FILEs   file names of Hi-C R2  [r2_1.fq,r2_2.fq,...]");
    eprintln!("    --seed INT   RNG seed [{}]", options.random_seed);
    eprintln!("    --s-base     FLOAT");
    eprintln!("                 similarity threshold for homology detection in base-level;");
    eprintln!(
        "                 -1 to disable [{}]; -s for read-level (see <Purge-dups>)",
        format_general(options.secondary_transition_base_rate as f64, 3)
    );
    eprintln!("    --n-weight   INT");
    eprintln!(
        "                 rounds of reweighting Hi-C links [{}]",
        options.reweighting_round_count
    );
    eprintln!("    --n-perturb  INT");
    eprintln!(
        "                 rounds of perturbation [{}]",
        options.perturbation_round_count
    );
    eprintln!("    --f-perturb  FLOAT");
    eprintln!(
        "                 fraction to flip for perturbation [{}]",
        format_general(options.perturbation_fraction, 3)
    );
    eprintln!("    --l-msjoin   INT");
    eprintln!(
        "                 detect misjoined unitigs of >=INT in size; 0 to disable [{}]",
        options.misjoin_length
    );
    eprintln!("  Ultra-Long-integration:");
    eprintln!("    --ul FILEs   file names of Ultra-Long reads [r1.fq,r2.fq,...]");
    eprintln!("    --ul-rate    FLOAT");
    eprintln!(
        "                 error rate of Ultra-Long reads [{}]",
        format_general(options.ultra_long_error_rate, 3)
    );
    eprintln!("    --ul-tip     INT");
    eprintln!(
        "                 remove tip unitigs composed of <=INT reads for the UL assembly [{}]",
        options.maximum_short_ultra_long_tip_read_count
    );
    eprintln!("    --path-max   FLOAT");
    eprintln!(
        "                 max path drop ratio [{}]; higher number may make the assembly cleaner",
        format_general(options.maximum_path_drop_rate as f64, 2)
    );
    eprintln!("                 but may lead to more misassemblies");
    eprintln!("    --path-min   FLOAT");
    eprintln!(
        "                 min path drop ratio [{}]; higher number may make the assembly cleaner",
        format_general(options.minimum_path_drop_rate as f64, 2)
    );
    eprintln!("                 but may lead to more misassemblies");
    eprintln!("    --ul-cut     INT");
    eprintln!(
        "                 filter out <INT UL reads during the UL assembly [{}]",
        options.ultra_long_minimum_base_count
    );
    eprintln!("  Dual-Scaffolding:");
    eprintln!("    --dual-scaf  output scaffolding");
    eprintln!("    --scaf-gap   INT");
    eprintln!(
        "                 max gap size for scaffolding [{}]",
        options.self_scaffolding_maximum_gap_length
    );
    eprintln!("  Telomere-identification:");
    eprintln!("    --telo-m     STR");
    eprintln!(
        "                 telomere motif at 5'-end; CCCTAA for human [{}]",
        options.telomere_motif.clone().unwrap_or_else(|| "NULL".to_string())
    );
    eprintln!("    --telo-p     INT");
    eprintln!("                 non-telomeric penalty [{}]", options.telomere_penalty);
    eprintln!("    --telo-d     INT");
    eprintln!("                 max drop [{}]", options.telomere_maximum_drop);
    eprintln!("    --telo-s     INT");
    eprintln!("                 min score for telomere reads [{}]", options.telomere_minimum_score);
    eprintln!("  ONT simplex assembly (beta):");
    eprintln!("    --ont        assemble ONT simplex reads in fastq format");
    eprintln!("    --chem-c     INT");
    eprintln!(
        "                 detect chimeric reads with <=INT other reads support [{}]",
        options.chemical_coverage
    );
    eprintln!("    --chem-f     INT");
    eprintln!(
        "                 length of flanking regions for chimeric read detection [{}]",
        options.chemical_flank_length
    );
    eprintln!("    --rl-cut     INT");
    eprintln!(
        "                 filter out ONT simplex reads shorter than <INT> for assembly [{}]",
        options.read_length_cutoff
    );
    eprintln!("    --sc-cut     INT");
    eprintln!(
        "                 filter out ONT simplex reads with a mean base quality score below <INT> [{}]",
        options.base_quality_cutoff
    );
    eprintln!("Example: ./hifiasm -o NA12878.asm -t 32 NA12878.fq.gz");
    eprintln!(
        "See `https://hifiasm.readthedocs.io/en/latest/' or `man ./hifiasm.1' for complete documentation."
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temporary_input_file(test_name: &str) -> String {
        let path = std::env::temp_dir().join(format!("rshifiasm_parser_{}.fa", test_name));
        std::fs::write(&path, b">read1\nACGTACGT\n").unwrap();
        path.to_string_lossy().to_string()
    }

    fn parse(arguments: &[&str]) -> (AssemblyOptions, bool) {
        let mut option_list = AssemblyOptions::default();
        // The first element is the program name, exactly as in `main`.
        let mut owned: Vec<String> = vec!["rshifiasm".to_string()];
        owned.extend(arguments.iter().map(|value| value.to_string()));
        let parsed = parse_arguments(&owned, &mut option_list);
        (option_list, parsed)
    }

    #[test]
    fn short_options_with_separate_argument_are_recognised() {
        let input = temporary_input_file("short_separate");
        let (options, ok) = parse(&["-t", "8", "-o", "out", &input]);
        assert!(ok);
        assert_eq!(options.thread_count, 8);
        assert_eq!(options.output_file_name, "out");
        assert_eq!(options.input_file_count, 1);
    }

    #[test]
    fn short_option_cluster_with_attached_argument_is_recognised() {
        let input = temporary_input_file("short_cluster");
        let (options, ok) = parse(&["-t4", "-k31", &input]);
        assert!(ok);
        assert_eq!(options.thread_count, 4);
        assert_eq!(options.k_mer_length, 31);
    }

    #[test]
    fn long_option_accepts_equals_form() {
        let input = temporary_input_file("long_equals");
        let (options, ok) = parse(&["--hom-cov=36", &input]);
        assert!(ok);
        assert_eq!(options.global_homozygous_coverage, 36);
        assert_eq!(options.global_homozygous_coverage_is_set, 1);
    }

    #[test]
    fn aborted_long_option_is_accepted_when_unambiguous() {
        let input = temporary_input_file("long_abbreviated");
        let (options, ok) = parse(&["--hom-c=36", &input]);
        assert!(ok);
        assert_eq!(options.global_homozygous_coverage, 36);
    }

    #[test]
    fn atoi_style_parsing_stops_at_the_first_invalid_character() {
        assert_eq!(parse_integer("42abc"), 42);
        assert_eq!(parse_integer("-7"), -7);
        assert_eq!(parse_integer("abc"), 0);
    }

    #[test]
    fn atof_style_parsing_handles_exponents_and_signs() {
        assert_eq!(parse_double("0.04x"), 0.04);
        assert_eq!(parse_double("-1.5e2"), -150.0);
        assert_eq!(parse_double("nonsense"), 0.0);
    }

    #[test]
    fn genome_size_suffixes_are_expanded() {
        assert_eq!(parse_genome_size("3g"), 3_000_000_000);
        assert_eq!(parse_genome_size("2.5m"), 2_500_000);
        assert_eq!(parse_genome_size("100"), -2);
    }
}

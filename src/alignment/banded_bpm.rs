//! The banded bit-parallel edit distance used by window verification,
//! `Reserve_Banded_BPM` in `Levenshtein_distance.h`.
//!
//! Myers' algorithm runs over a band of `2 * error_threshold + 1` pattern
//! positions. The pattern is the (longer) target window and the text is the
//! (shorter) query window; the result is the position in the pattern where the
//! best alignment ends, together with its edit distance.

/// Value of `return_err` when no alignment within the threshold exists, the
/// reference `(unsigned int)-1`.
pub const NO_ALIGNMENT_ERROR: u32 = u32::MAX;

/// Result of `Reserve_Banded_BPM`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BandedAlignment {
    /// End position in the pattern, `-1` when the alignment failed.
    pub end_site: i32,
    /// Edit distance of the reported alignment, [`NO_ALIGNMENT_ERROR`] on failure.
    pub error: u32,
}

impl BandedAlignment {
    /// The failed alignment.
    pub fn failure() -> Self {
        Self {
            end_site: -1,
            error: NO_ALIGNMENT_ERROR,
        }
    }
}

/// The four slots that are shifted as the band slides.
///
/// The reference keeps a 256 entry equality table: the initial band only
/// survives for A/C/G/T (the table is wiped afterwards), only those four are
/// shifted every step, but a pattern byte of any other value (notably `N`)
/// still receives a band bit when the band slides and can match the same byte
/// in the text. Keeping the full table reproduces that behaviour.
const SHIFTED_BASES: [usize; 4] = [b'A' as usize, b'C' as usize, b'G' as usize, b'T' as usize];

/// Advances one cell of the Myers matrix and returns the lowest bit of `D0`.
/// A zero bit means the edit distance grew by one.
fn advance_cell(
    equality_masks: &[u64; 256],
    text_base: u8,
    positive_vertical: &mut u64,
    negative_vertical: &mut u64,
) -> u64 {
    let equality = equality_masks[text_base as usize];
    let horizontal = equality | *negative_vertical;
    let d0 = ((positive_vertical.wrapping_add(horizontal & *positive_vertical))
        ^ *positive_vertical)
        | horizontal;
    let horizontal_negative = *positive_vertical & d0;
    let horizontal_positive = *negative_vertical | !(*positive_vertical | d0);
    let shifted = d0 >> 1;
    *negative_vertical = shifted & horizontal_positive;
    *positive_vertical = horizontal_negative | !(shifted | horizontal_positive);
    d0 & 1
}

/// `Reserve_Banded_BPM(pattern, p_length, text, t_length, errthold, &err)`.
///
/// `pattern` and `text` are read over exactly `pattern_length` and
/// `text_length` bytes. `error_threshold` is at most 31 in the reference, so
/// the band always fits into one 64-bit word.
pub fn reserve_banded_bpm(
    pattern: &[u8],
    pattern_length: usize,
    text: &[u8],
    text_length: usize,
    error_threshold: usize,
) -> BandedAlignment {
    if text_length == 0 || pattern_length == 0 {
        return BandedAlignment::failure();
    }

    // Build the band: bit i of the mask of base b is set when pattern[i] == b.
    let band_length = (error_threshold << 1) + 1;
    let mut equality_masks = [0u64; 256];
    let mut band_bit = 1u64;
    for index in 0..band_length {
        if let Some(base) = pattern.get(index) {
            // The reference builds the band into a 256 slot table but then
            // wipes everything except A/C/G/T, so an `N` in the first band
            // positions never matches. Later band positions are different:
            // they are added after the wipe and keep their bit.
            if SHIFTED_BASES.contains(&(*base as usize)) {
                equality_masks[*base as usize] |= band_bit;
            }
        }
        band_bit <<= 1;
    }

    // The highest bit of the band is reintroduced at every step as the band
    // slides one position to the right.
    let band_top_bit = 1u64 << (error_threshold << 1);
    let mut positive_vertical = 0u64;
    let mut negative_vertical = 0u64;
    let mut error: i32 = 0;
    let threshold = error_threshold as i32;
    let last_high = (error_threshold << 1) as i32;
    let mut band_position = error_threshold << 1;

    let mut text_index = 0usize;
    while text_index < text_length - 1 {
        let text_base = text.get(text_index).copied().unwrap_or(b'N');
        let lowest_bit = advance_cell(
            &equality_masks,
            text_base,
            &mut positive_vertical,
            &mut negative_vertical,
        );
        if lowest_bit == 0 {
            error += 1;
            if error - last_high > threshold {
                return BandedAlignment::failure();
            }
        }

        for slot in SHIFTED_BASES {
            equality_masks[slot] >>= 1;
        }
        text_index += 1;
        band_position += 1;
        if let Some(base) = pattern.get(band_position) {
            equality_masks[*base as usize] |= band_top_bit;
        }
    }

    // The last text base is handled outside the loop in the reference and does
    // not slide the band any further.
    let text_base = text.get(text_index).copied().unwrap_or(b'N');
    let lowest_bit = advance_cell(
        &equality_masks,
        text_base,
        &mut positive_vertical,
        &mut negative_vertical,
    );
    if lowest_bit == 0 {
        error += 1;
        if error - last_high > threshold {
            return BandedAlignment::failure();
        }
    }

    let site = text_length as i32 - 1;
    let mut return_site = -1i32;
    let mut return_error = NO_ALIGNMENT_ERROR;

    // The alignment may end at any of the `p_length - t_length` positions past
    // the text end; walking VP/VN gives the edit distance at each of them.
    let available_positions = pattern_length as i32 - text_length as i32;
    if error <= threshold && (error as u32) <= return_error {
        return_error = error as u32;
        return_site = site;
    }

    let mut ungap_error = NO_ALIGNMENT_ERROR;
    let mut position = 0i32;
    while position < available_positions {
        error += ((positive_vertical.checked_shr(position as u32).unwrap_or(0)) & 1) as i32;
        error -= ((negative_vertical.checked_shr(position as u32).unwrap_or(0)) & 1) as i32;
        position += 1;

        if error <= threshold && (error as u32) <= return_error {
            return_error = error as u32;
            return_site = site + position;
        }
        if position == threshold {
            ungap_error = error as u32;
        }
    }

    // Prefer the ungapped end when it has the same error as the best one.
    if ungap_error <= threshold as u32 && ungap_error == return_error {
        return_site = site + threshold;
    }

    BandedAlignment {
        end_site: return_site,
        error: return_error,
    }
}

/// Result of `alignment_extension`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExtensionAlignment {
    /// Edit distance of the reported alignment, [`NO_ALIGNMENT_ERROR`] on
    /// failure.
    pub error: u32,
    /// End position in the pattern, `-1` when the alignment failed.
    pub pattern_end: i32,
    /// End position in the text, `-1` when the alignment failed.
    pub text_end: i32,
    /// Number of text bases the alignment covers, `0` on failure.
    pub aligned_text_length: usize,
}

impl ExtensionAlignment {
    fn failure() -> Self {
        Self {
            error: NO_ALIGNMENT_ERROR,
            pattern_end: -1,
            text_end: -1,
            aligned_text_length: 0,
        }
    }
}

/// `get_error`: the best error of the current column when the alignment is
/// allowed to stop anywhere inside the band.
///
/// `VP`/`VN` describe how the edit distance changes as the end position moves
/// one pattern base to the right, so walking `2 * error_threshold` bits gives
/// the error of every ungapped end. Ties keep the earliest end, except that the
/// end at exactly `error_threshold` wins when its error equals the best one.
fn best_band_error(
    text_length: i32,
    error_threshold: i32,
    initial_error: i32,
    positive_vertical: u64,
    negative_vertical: u64,
) -> (u32, i32) {
    let site = text_length - 1;
    let available_positions = error_threshold << 1;
    let mut return_error = NO_ALIGNMENT_ERROR;
    let mut return_site = -1i32;

    let mut error = initial_error;
    if error <= error_threshold && (error as u32) <= return_error {
        return_error = error as u32;
        return_site = site;
    }

    let mut ungap_error = NO_ALIGNMENT_ERROR;
    let mut position = 0i32;
    while position < available_positions {
        error += ((positive_vertical >> position) & 1) as i32;
        error -= ((negative_vertical >> position) & 1) as i32;
        position += 1;

        if error <= error_threshold && (error as u32) <= return_error {
            return_error = error as u32;
            return_site = site + position;
        }
        if position == error_threshold {
            ungap_error = error as u32;
        }
    }

    if ungap_error <= error_threshold as u32 && ungap_error == return_error {
        return_site = site + error_threshold;
    }

    (return_error, return_site)
}

/// `Reserve_Banded_BPM_Extension`: like [`reserve_banded_bpm`], but reports
/// both ends of the alignment instead of only the pattern end.
///
/// Returns `(error, pattern_end, text_end)`; the two ends are `-1` when no
/// alignment within the threshold exists.
pub fn reserve_banded_bpm_extension(
    pattern: &[u8],
    pattern_length: usize,
    text: &[u8],
    text_length: usize,
    error_threshold: usize,
) -> (u32, i32, i32) {
    if text_length == 0 || pattern_length == 0 {
        return (NO_ALIGNMENT_ERROR, -1, -1);
    }

    let band_length = (error_threshold << 1) + 1;
    let mut equality_masks = [0u64; 256];
    let mut band_bit = 1u64;
    for index in 0..band_length {
        if let Some(base) = pattern.get(index) {
            if SHIFTED_BASES.contains(&(*base as usize)) {
                equality_masks[*base as usize] |= band_bit;
            }
        }
        band_bit <<= 1;
    }

    let band_top_bit = 1u64 << (error_threshold << 1);
    let mut positive_vertical = 0u64;
    let mut negative_vertical = 0u64;
    let mut error: i32 = 0;
    let threshold = error_threshold as i32;
    let last_high = (error_threshold << 1) as i32;
    let mut band_position = error_threshold << 1;

    let mut return_error = NO_ALIGNMENT_ERROR;
    let mut return_pattern_end = -1i32;
    let mut return_text_end = -1i32;

    let mut text_index = 0usize;
    while text_index + 1 < text_length {
        let text_base = text.get(text_index).copied().unwrap_or(b'N');
        let lowest_bit = advance_cell(
            &equality_masks,
            text_base,
            &mut positive_vertical,
            &mut negative_vertical,
        );
        if lowest_bit == 0 {
            error += 1;
            if error - last_high > threshold {
                return (return_error, return_pattern_end, return_text_end);
            }
        }

        let (line_error, return_site) = best_band_error(
            text_index as i32 + 1,
            threshold,
            error,
            positive_vertical,
            negative_vertical,
        );
        if line_error != NO_ALIGNMENT_ERROR {
            return_text_end = text_index as i32;
            return_pattern_end = return_site;
            return_error = line_error;
        }

        for slot in SHIFTED_BASES {
            equality_masks[slot] >>= 1;
        }
        text_index += 1;
        band_position += 1;
        if let Some(base) = pattern.get(band_position) {
            equality_masks[*base as usize] |= band_top_bit;
        }
    }

    let text_base = text.get(text_index).copied().unwrap_or(b'N');
    let lowest_bit = advance_cell(
        &equality_masks,
        text_base,
        &mut positive_vertical,
        &mut negative_vertical,
    );
    if lowest_bit == 0 {
        error += 1;
        if error - last_high > threshold {
            return (return_error, return_pattern_end, return_text_end);
        }
    }

    let (line_error, return_site) = best_band_error(
        text_index as i32 + 1,
        threshold,
        error,
        positive_vertical,
        negative_vertical,
    );
    if line_error != NO_ALIGNMENT_ERROR {
        return_text_end = text_index as i32;
        return_pattern_end = return_site;
        return_error = line_error;
    }

    (return_error, return_pattern_end, return_text_end)
}

/// `alignment_extension`: aligns from one end of the two windows.
///
/// `direction == 0` walks forwards from the start of both windows; any other
/// value walks backwards from their ends, which the reference implements by
/// reversing both strings around the call. The reported ends are translated
/// back into forward coordinates in that case.
pub fn alignment_extension(
    pattern: &[u8],
    pattern_length: usize,
    text: &[u8],
    text_length: usize,
    error_threshold: usize,
    direction: i32,
) -> ExtensionAlignment {
    if direction == 0 {
        let (error, pattern_end, text_end) = reserve_banded_bpm_extension(
            pattern,
            pattern_length,
            text,
            text_length,
            error_threshold,
        );
        if pattern_end != -1 && text_end != -1 {
            return ExtensionAlignment {
                error,
                pattern_end,
                text_end,
                aligned_text_length: text_end as usize + 1,
            };
        }
        return ExtensionAlignment::failure();
    }

    let reversed_pattern: Vec<u8> = pattern[..pattern_length].iter().rev().copied().collect();
    let reversed_text: Vec<u8> = text[..text_length].iter().rev().copied().collect();
    let (error, pattern_end, text_end) = reserve_banded_bpm_extension(
        &reversed_pattern,
        pattern_length,
        &reversed_text,
        text_length,
        error_threshold,
    );
    if pattern_end == -1 || text_end == -1 {
        return ExtensionAlignment::failure();
    }

    ExtensionAlignment {
        error,
        pattern_end: pattern_length as i32 - pattern_end,
        text_end: text_length as i32 - text_end,
        aligned_text_length: text_end as usize + 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn advance_cell_produces_the_myers_bit_vectors() {
        let mut masks = [0u64; 256];
        masks[b'A' as usize] = 1; // the first pattern base is A
        let mut positive = 0u64;
        let mut negative = 0u64;

        let lowest = advance_cell(&masks, b'A', &mut positive, &mut negative);
        assert_eq!(lowest, 1, "a match does not increment the error");

        let lowest = advance_cell(&masks, b'C', &mut positive, &mut negative);
        assert_eq!(lowest, 0, "a mismatch increments the error");
    }

    /// Golden values printed by the real `alignment_extension`, both
    /// directions. The second direction walks the windows backwards, which is
    /// how `verify_sub_window` extends a window to the left.
    #[test]
    fn matches_the_reference_extension_alignment() {
        let cases: &[(&str, &str, usize, usize, usize, i32, u32, i32, i32, usize)] = &[
            ("ACGTACGTACGGAA", "GTACGTACGG", 14usize, 10usize, 2usize, 0i32, 0u32, 11i32, 9i32, 10usize),
            ("ACGTACGTACGGAA", "GTACGTACGG", 14, 10, 2, 1, 0u32, 3, 1, 10),
            ("ACGTACGTACGGAA", "GTACGTACAG", 14, 10, 2, 0, 1u32, 11, 9, 10),
            ("ACGTACGTACGGAA", "GTACGTACAG", 14, 10, 2, 1, 1u32, 3, 1, 10),
            ("ACGTACGTACGGAA", "TTTTTTTTTT", 14, 10, 2, 0, 2u32, 4, 2, 3),
            ("ACGTACGTACGGAA", "TTTTTTTTTT", 14, 10, 2, 1, 2u32, 8, 8, 3),
            ("ACGTACGTACGGAA", "GTANGTACGG", 14, 10, 2, 0, 1u32, 11, 9, 10),
            ("ACGTACGTACGGAA", "GTANGTACGG", 14, 10, 2, 1, 1u32, 3, 1, 10),
            ("TTACGTACGTAAAA", "ACGTACGTAA", 14, 10, 2, 0, 0u32, 11, 9, 10),
            ("TTACGTACGTAAAA", "ACGTACGTAA", 14, 10, 2, 1, 0u32, 3, 1, 10),
            ("AACCGGTTAACCGGTT", "CCGGTTAACCGGTT", 16, 14, 4, 0, 0u32, 15, 13, 14),
            ("AACCGGTTAACCGGTT", "CCGGTTAACCGGTT", 16, 14, 4, 1, 0u32, 3, 1, 14),
            ("ACGTACGTACGTACGT", "ACGTACGTACGTACGT", 16, 16, 0, 0, 0u32, 15, 15, 16),
            ("ACGTACGTACGTACGT", "ACGTACGTACGTACGT", 16, 16, 0, 1, 0u32, 1, 1, 16),
            ("NNNNACGTACGTACGT", "ACGTACGTACGT", 16, 12, 3, 0, 0u32, 15, 11, 12),
            ("NNNNACGTACGTACGT", "ACGTACGTACGT", 16, 12, 3, 1, 0u32, 5, 1, 12),
        ];

        for (pattern_text, text_text, pattern_length, text_length, threshold, direction, error, pattern_end, text_end, aligned) in
            cases
        {
            let pattern = padded(pattern_text, *pattern_length);
            let text = padded(text_text, *text_length);
            let result = alignment_extension(
                &pattern,
                *pattern_length,
                &text,
                *text_length,
                *threshold,
                *direction,
            );
            let label = format!("{pattern_text}/{text_text} dir {direction}");
            assert_eq!(result.error, *error, "error of {label}");
            assert_eq!(result.pattern_end, *pattern_end, "pattern end of {label}");
            assert_eq!(result.text_end, *text_end, "text end of {label}");
            assert_eq!(
                result.aligned_text_length, *aligned,
                "aligned length of {label}"
            );
        }

        // Even a hopeless pair reports the best one-base alignment when the
        // band allows a single error: the reference returns `pe=1` forwards and
        // `te=10` backwards, the latter being one past the text end.
        let pattern = padded("AAAAAAAAAAAAAAAA", 16);
        let text = padded("TTTTTTTTTT", 10);
        let forward = alignment_extension(&pattern, 16, &text, 10, 1, 0);
        assert_eq!(
            (forward.error, forward.pattern_end, forward.text_end, forward.aligned_text_length),
            (1, 1, 0, 1)
        );
        let backward = alignment_extension(&pattern, 16, &text, 10, 1, 1);
        assert_eq!(
            (backward.error, backward.pattern_end, backward.text_end, backward.aligned_text_length),
            (1, 15, 10, 1)
        );
    }

    /// Builds the buffers exactly as the C probe does: the string bytes are
    /// copied and the remaining positions stay at `A`.
    fn padded(text: &str, length: usize) -> Vec<u8> {
        let mut buffer = vec![b'A'; length + 8];
        buffer[..text.len()].copy_from_slice(text.as_bytes());
        buffer
    }

    /// Golden values printed by the real `Reserve_Banded_BPM`.
    #[test]
    fn matches_the_reference_banded_alignment() {
        let cases = [
            // name, pattern, pattern length, text, text length, threshold,
            // expected end site, expected error
            ("identical", "ACGTACGTACGG", 14usize, "GTACGTACGG", 10usize, 2usize, 11i32, 0u32),
            ("one_mismatch", "ACGTACGTACGG", 14, "GTACGTACAG", 10, 2, 11, 1),
            (
                "three_mismatches",
                "ACGTACGTACGG",
                14,
                "TTTTTTTTTT",
                10,
                2,
                -1,
                NO_ALIGNMENT_ERROR,
            ),
            ("with_N", "ACGTACGTACGG", 14, "GTANGTACGG", 10, 2, 11, 1),
            ("threshold_zero", "ACGTACGTACGG", 10, "ACGTACGTAC", 10, 0, 9, 0),
            ("shifted", "TTACGTACGTAA", 14, "ACGTACGTAA", 10, 2, 11, 0),
            ("N_both", "ACNTACGTACGG", 14, "NTACGTACGG", 10, 2, 11, 1),
            ("N_in_text_only", "ACGTACGTACGG", 14, "GTANGTACGG", 10, 2, 11, 1),
            ("N_in_pattern_only", "ACNTACGTACGG", 14, "GTACGTACGG", 10, 2, 11, 1),
            ("longer_band", "AACCGGTTAACCGGTTAA", 24, "CCGGTTAACCGGTT", 20, 4, 21, 0),
            ("slide_N_both", "ACGTANGTACGG", 14, "GTANGTACGG", 10, 2, 11, 1),
            ("slide_N_pattern", "ACGTANGTACGG", 14, "GTACGTACGG", 10, 2, 11, 1),
            ("slide_N_text", "ACGTACGTACGG", 14, "GTANGTACGG", 10, 2, 11, 1),
            ("slide_N_late", "ACGTACGTANGG", 14, "GTACGTACGG", 10, 2, 11, 1),
        ];

        for (name, pattern_text, pattern_length, text_text, text_length, threshold, site, error) in
            cases
        {
            let pattern = padded(pattern_text, pattern_length);
            let text = padded(text_text, text_length);
            let result = reserve_banded_bpm(
                &pattern,
                pattern_length,
                &text,
                text_length,
                threshold,
            );
            assert_eq!(result.end_site, site, "end site of {name}");
            assert_eq!(result.error, error, "error of {name}");
        }
    }

    /// Forty deterministic random cases produced by the reference
    /// `Reserve_Banded_BPM`; the strings include `N` and the thresholds cover
    /// 1..8, so this exercises the band sliding and the ungapped-end tie-break.
    #[test]
    fn matches_the_reference_on_random_sequences() {
        let cases: &[(&str, &str, usize, usize, usize, i32, u32)] = &[
            ("NGGCTNGNNGCCCCGACCAGCNCG", "NGACTNGNNGCCCCGACCAGCN", 24usize, 22usize, 1usize, -1i32, NO_ALIGNMENT_ERROR),
            ("ACCNGAGNTTNNAAACNANCATTCTGNCATA", "GTCGTNGNNAAGNNTGCCGTGAN", 31usize, 23usize, 4usize, -1i32, NO_ALIGNMENT_ERROR),
            ("TCTCNGGCGANCNNTATNATCANCATGAACGGCTC", "CNGGCGANCNNTATNATCANCAT", 35usize, 23usize, 6usize, 25i32, 6u32),
            ("GTGGGNCGANTCCATNCGTACCCTGNNTACGTATTTNCTACGGCNANA", "TGAACNNGTNCTNCCANATTNANTTGATTAGN", 48usize, 32usize, 8usize, -1i32, NO_ALIGNMENT_ERROR),
            ("GNTGGGCNNTTGTTCCGTTGAANGAA", "NGCNNTTGTTCCGTTGAANG", 26usize, 20usize, 3usize, -1i32, NO_ALIGNMENT_ERROR),
            ("TTTTCCTCCTACNCATGNNNGANAGGTAACTGTTCGNTTNNAANAANTNT", "GGNTATAATGCNCGCNANNAAAACACAGACGGTTATCGTTAG", 50usize, 42usize, 4usize, -1i32, NO_ALIGNMENT_ERROR),
            ("GNNGCANCCGCAGNCNACGTGGGCTGATNTGTN", "NNGCANCCGCAGNCNACGTGGGC", 33usize, 23usize, 5usize, 23i32, 5u32),
            ("TNNACNGNTNNCAGCCTAGNNAAGAT", "GGCGGTCCGA", 26usize, 10usize, 8usize, 17i32, 5u32),
            ("AGNTNCANTAGCCGTGCCATNTGANACTATAGAACNCCTAANNC", "GCCGTGCCATNTGANACTATATAACNCCTA", 44usize, 30usize, 7usize, 39i32, 4u32),
            ("NNGTCCCAANGTGNACNCGGCAGTCNCCNACTGNTANGCNT", "NCNCGCAGTACNCANCGCTTANGACCGGTCGTTGCGGCA", 41usize, 39usize, 1usize, -1i32, NO_ALIGNMENT_ERROR),
            ("NANCNCACCTCTGCNTGTGGTATGCANACCCAAA", "TCTGCNTGTGGTATGCANAC", 34usize, 20usize, 7usize, 28i32, 2u32),
            ("ATCTCCGTCGCCTNTCCGGCTNTNAGCATTACN", "NNAACCTACTCCANNAGTNNAACTACATCGG", 33usize, 31usize, 1usize, -1i32, NO_ALIGNMENT_ERROR),
            ("GGTTCCATTCGANTGGTTTGACCGGCCANGTCACNNGTCTGNTCN", "CCATTAGANTGGTTTGACCGGCCANGTCA", 45usize, 29usize, 8usize, 32i32, 3u32),
            ("CAAAACTAANNTANNAC", "GTGCNNGTTGG", 17usize, 11usize, 3usize, -1i32, NO_ALIGNMENT_ERROR),
            ("NGTTTNCGNCAAGNCANNTTTAAGAAGACGTNCGCANANTATTNTNTCAAC", "CGNCAAGNCANNTTTAAGAAGACGTNCGCANANTATTNTNTCA", 51usize, 43usize, 4usize, -1i32, NO_ALIGNMENT_ERROR),
            ("CCCAATGCTGCGNGTCCGNTANTG", "CAGCNCTGAACACTAC", 24usize, 16usize, 4usize, -1i32, NO_ALIGNMENT_ERROR),
            ("CGATGGTGTTNNAGCANGANTTACCTCCTCGTGANATGAGGGGACGAGAT", "GTTNNAGCANGANTTACCTCCTCGTGANATGAGG", 50usize, 34usize, 8usize, 40i32, 5u32),
            ("TGTNGAGGNNGANTCGCATGGNNC", "AAACANTCNCGCGACN", 24usize, 16usize, 4usize, -1i32, NO_ALIGNMENT_ERROR),
            ("CGTNGNTTAAGCGGGTCCGANTNACCANGACAAATTTNNTTNCCNNGAGGTAANACGGAT", "GNTTAAGCGGGTCCGANTNACCANGACAAATTTNNTTNCCNNGAGG", 60usize, 46usize, 7usize, -1i32, NO_ALIGNMENT_ERROR),
            ("NTTTCCACGAACTTCCANANTTTCTCTGTATTCTCAGGTGTCNCGTCCGGT", "ANCGGCGNGGGCTTNNNGCTGATGGANNATCGAAA", 51usize, 35usize, 8usize, -1i32, NO_ALIGNMENT_ERROR),
            ("TNCGNTTACATTGNAANTCG", "NCGNTTACGTTGNA", 20usize, 14usize, 3usize, -1i32, NO_ALIGNMENT_ERROR),
            ("TTNTTAATGGNNAGGCGTACTCGCNATNACTCTAT", "CGAATTAGCCCGTCNCACGCAGTCG", 35usize, 25usize, 5usize, -1i32, NO_ALIGNMENT_ERROR),
            ("TTCCTTGATTTATTATNGTTGCTNNGNNCNTNACCGAAANGCTNCN", "GATTTATTATNGTTGCTNNGNNCNTNACCGAAANGCTNCN", 46usize, 40usize, 3usize, 45i32, 0u32),
            ("TNAACNAAAGGCNCGTNTCCNCTAAG", "TTTTCATNAGGATAGG", 26usize, 16usize, 5usize, -1i32, NO_ALIGNMENT_ERROR),
            ("GCTTNGGTTTGAGNGTGGNTTAGGTTGGGNTGCNATNNCGTTACGGGNGAA", "NGGTTTGAGNGTGGNTTAGGTTNGGNTGCNATNNCGTTACGGG", 51usize, 43usize, 4usize, -1i32, NO_ALIGNMENT_ERROR),
            ("GTACNANACTNNGACNNTATNGTGACNTAT", "CNCTTGNTNNNCCTTANNGG", 30usize, 20usize, 5usize, -1i32, NO_ALIGNMENT_ERROR),
            ("GATGNAAATTAGACCATNTGNNTAANTTCGANANTGCGTTGGNTTANNT", "GATGNAAATTAGACCATNTGNNTAANTTCGANANTGCGTTGGNTTAN", 49usize, 47usize, 1usize, -1i32, NO_ALIGNMENT_ERROR),
            ("CACCGTCNCCCTNCGAGTNCCTCATTGCAGNTTNACGA", "CTCCNNANNGTAACNNGTANNTGGTNNATCGTCG", 38usize, 34usize, 2usize, -1i32, NO_ALIGNMENT_ERROR),
            ("NATNGAACCGANGAATNNTNTANCTTTATCNCNTNGTATTNC", "TAACCGANGAATNNTNTANCTTTATCNC", 42usize, 28usize, 7usize, 31i32, 7u32),
            ("TANCACTTACTCTCGCACTGANGT", "CGTANNACGC", 24usize, 10usize, 7usize, 15i32, 5u32),
            ("TAACNCCGAACNGTCAATANAATTGCNCACTACAGGACCAGAATC", "AACNCCGAACNGTCAATANAATTGCNCACTACAGGACCAGAAT", 45usize, 43usize, 1usize, -1i32, NO_ALIGNMENT_ERROR),
            ("ANNATTATNTTNGACANNGTANGNCCTGANANGNCA", "NAGTGCCCTNGCAAGGANTNCAANNTTCTN", 36usize, 30usize, 3usize, -1i32, NO_ALIGNMENT_ERROR),
            ("TAGNNGGCTGGTACNGACNGCTGANCATT", "ACTGACNGCTGANCATT", 29usize, 17usize, 6usize, 28i32, 1u32),
            ("TGGTTACNCNCNGTNNCNTCTCCNAANGA", "CNCTCCNCGTGTGTC", 29usize, 15usize, 7usize, -1i32, NO_ALIGNMENT_ERROR),
            ("ACCANAGNNGTTGNTACTGGCA", "GTTGNTACTG", 22usize, 10usize, 6usize, 18i32, 1u32),
            ("GNNANATCTGTNCNCGTAGCGTNANGN", "CTANTNACATNGTACGGCGGACC", 27usize, 23usize, 2usize, -1i32, NO_ALIGNMENT_ERROR),
            ("GCNNTACCTCNGCNTGCNGACCANGAATTNGTAAAGTT", "ACCTCNGCNTGCNGACGANGAA", 38usize, 22usize, 8usize, 26i32, 5u32),
            ("CGATTTCNAGGGAAGAGAANTTCANACNGANGTNGATAATTC", "ANATAGNTGANGAGTNGNNCNATGCGAGAGGA", 42usize, 32usize, 5usize, -1i32, NO_ALIGNMENT_ERROR),
            ("AAGAGNNNAAGGACTCTAGNTGTNCCNCCCNNTNGGCGTNTCTGCTTCG", "AAGAGNNNAAGGACTCTAGNTGTNCCNCCCNNTNGGCGTNTCTGCTT", 49usize, 47usize, 1usize, -1i32, NO_ALIGNMENT_ERROR),
            ("CAGAGNNNNANGAGNGGTGNTTANGNNCNGAGGNTNA", "CTTGTGNATNAGCNGCGACCTCC", 37usize, 23usize, 7usize, -1i32, NO_ALIGNMENT_ERROR),
        ];

        for (pattern_text, text_text, pattern_length, text_length, threshold, site, error) in cases
        {
            let pattern = padded(pattern_text, *pattern_length);
            let text = padded(text_text, *text_length);
            let result = reserve_banded_bpm(
                &pattern,
                *pattern_length,
                &text,
                *text_length,
                *threshold,
            );
            assert_eq!(result.end_site, *site, "end site of {text_text}");
            assert_eq!(result.error, *error, "error of {text_text}");
        }
    }
}

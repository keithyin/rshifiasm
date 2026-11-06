//! Bit-parallel banded edit distance, the counterpart of the `_64_`, `_128_`
//! and `_infi_` families in `Levenshtein_distance.h`.
//!
//! The original implements the same Myers bit-vector recurrence once per word
//! width through macros: the 64-bit kernels are written by hand, the 128-bit
//! ones come from a macro instantiation, and the `infi` ones operate on
//! `nword` words at run time. All three variants perform the same arithmetic
//! on the meaningful bits of the band, so the width-generic kernels here stand
//! in for every width the reference dispatches to. Bits above the band are
//! sentinel bits that no read touches, which is why one implementation gives
//! the same result as the specialised 64-bit and 128-bit kernels.

use crate::sequence::encoding::SEQUENCE_NT4_TABLE;

/// Error threshold and result of a bit-vector alignment, the part of
/// `bit_extz_t` this module needs.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EditDistanceResult {
    /// Allowed number of errors.
    pub threshold: i32,
    /// Best number of errors found, `INT32_MAX` when nothing qualifies.
    pub error: i32,
    /// Pattern start in the reference, `-1` when unset.
    pub pattern_start: i32,
    /// Pattern end in the reference, `-1` when unset.
    pub pattern_end: i32,
    /// Text start in the reference, `0` when unset.
    pub text_start: i32,
    /// Text end in the reference, `-1` when unset.
    pub text_end: i32,
    /// Pattern length the alignment was run with.
    pub pattern_length: i32,
    /// Text length the alignment was run with.
    pub text_length: i32,
}

impl EditDistanceResult {
    /// `init_base_ed`
    fn reset(&mut self, threshold: i32, pattern_length: i32, text_length: i32) {
        self.threshold = threshold;
        self.error = i32::MAX;
        self.pattern_length = pattern_length;
        self.text_length = text_length;
    }
}

/// The three bit vectors `ed_core_64` produces and that callers need.
pub(crate) struct CoreVectors {
    /// `D0`, whose lowest bit tells whether the error count advances.
    pub diagonal_zero: u64,
    /// `HN`, the vertical negative vector before the update.
    pub horizontal_negative: u64,
    /// `HP`, the vertical positive vector before the update.
    pub horizontal_positive: u64,
}

/// `ed_core_64`: one Myers step over a single text base.
fn edit_core_64(
    peq: &[u64; 5],
    vertical_positive: &mut u64,
    vertical_negative: &mut u64,
    text_base: u8,
) -> CoreVectors {
    let code = SEQUENCE_NT4_TABLE[text_base as usize];
    let x = peq[code as usize] | *vertical_negative;
    let diagonal_zero =
        (vertical_positive.wrapping_add(x & *vertical_positive) ^ *vertical_positive) | x;
    let horizontal_negative = *vertical_positive & diagonal_zero;
    let horizontal_positive = *vertical_negative | !(*vertical_positive | diagonal_zero);
    let shifted = diagonal_zero >> 1;
    *vertical_negative = shifted & horizontal_positive;
    *vertical_positive = horizontal_negative | !(shifted | horizontal_positive);
    CoreVectors {
        diagonal_zero,
        horizontal_negative,
        horizontal_positive,
    }
}

/// Builds the initial pattern bit vectors (`ed_init_core` plus the loop around
/// it, which differs slightly per variant).
fn build_pattern_vectors(pattern: &[u8], first_length: i32, first_bit: i32) -> [u64; 5] {
    let mut peq = [0u64; 5];
    let mut mask = 1u64 << first_bit;
    for index in 0..first_length {
        let code = SEQUENCE_NT4_TABLE[pattern[index as usize] as usize];
        peq[code as usize] |= mask;
        mask <<= 1;
    }
    peq[4] = 0;
    peq
}

/// `ed_band_cal_global_64_w`: global alignment of `pattern` against `text`.
pub fn global_64(
    pattern: &[u8],
    pattern_length: i32,
    text: &[u8],
    text_length: i32,
    threshold: i32,
    result: &mut EditDistanceResult,
) {
    result.reset(threshold, pattern_length, text_length);
    result.pattern_start = 0;
    result.text_start = 0;
    let pn = pattern_length;
    let tn = text_length;
    if pn > tn + threshold || tn > pn + threshold {
        return;
    }

    let last_text_index = tn - 1;
    let cutoff = threshold + (threshold << 1);
    let first_length = (threshold + 1).min(pn);
    let mut peq = build_pattern_vectors(pattern, first_length, threshold);

    let mut error = threshold;
    let mut vertical_negative = (1u64 << threshold) - 1;
    let mut vertical_positive = (1u64 << ((threshold << 1) + 1)) - 1;
    vertical_positive ^= vertical_negative;

    let mut text_index = 0i32;
    let mut pattern_index = threshold;
    let new_bit = 1u64 << (threshold << 1);

    while text_index < last_text_index {
        let d0 = edit_core_64(
            &peq,
            &mut vertical_positive,
            &mut vertical_negative,
            text[text_index as usize],
        );
        if d0.diagonal_zero & 1 == 0 {
            error += 1;
            if error > cutoff {
                return;
            }
        }

        peq[0] >>= 1;
        peq[1] >>= 1;
        peq[2] >>= 1;
        peq[3] >>= 1;
        text_index += 1;
        pattern_index += 1;
        if pattern_index < pn {
            let code = SEQUENCE_NT4_TABLE[pattern[pattern_index as usize] as usize];
            if code < 4 {
                peq[code as usize] |= new_bit;
            }
        }
    }

    let d0 = edit_core_64(
        &peq,
        &mut vertical_positive,
        &mut vertical_negative,
        text[last_text_index as usize],
    );
    if d0.diagonal_zero & 1 == 0 {
        error += 1;
        if error > cutoff {
            return;
        }
    }

    let mut site = tn - 1 - threshold;
    let final_site = pn - 1;
    while site < final_site {
        error += (vertical_positive & 1) as i32;
        vertical_positive >>= 1;
        error -= (vertical_negative & 1) as i32;
        vertical_negative >>= 1;
        site += 1;
    }

    if site == final_site && error <= threshold {
        result.error = error;
        result.pattern_end = pn - 1;
        result.text_end = tn - 1;
    }
}

/// `ed_band_cal_semi_64_w`: semi-global alignment, the pattern must be matched
/// end to end while the text may be entered and left freely.
pub fn semi_global_64(
    pattern: &[u8],
    pattern_length: i32,
    text: &[u8],
    text_length: i32,
    threshold: i32,
    result: &mut EditDistanceResult,
) {
    result.reset(threshold, pattern_length, text_length);
    result.pattern_start = -1;
    result.pattern_end = -1;
    result.text_start = 0;
    result.text_end = text_length - 1;

    let pn = pattern_length;
    let tn = text_length;
    let last_high = threshold << 1;
    let cutoff = threshold + last_high;
    if pn > tn + cutoff || tn > pn + cutoff {
        return;
    }

    let first_length = last_high + 1;
    let first_length = if first_length <= pn { first_length } else { pn };
    let mut peq = build_pattern_vectors(pattern, first_length, 0);

    let last_text_index = tn - 1;
    let mut error = 0i32;
    let mut vertical_positive = 0u64;
    let mut vertical_negative = 0u64;
    let mut pattern_index = last_high;
    let new_bit = 1u64 << last_high;

    let mut text_index = 0i32;
    while text_index < last_text_index {
        let d0 = edit_core_64(
            &peq,
            &mut vertical_positive,
            &mut vertical_negative,
            text[text_index as usize],
        );
        if d0.diagonal_zero & 1 == 0 {
            error += 1;
            if error > cutoff {
                return;
            }
        }

        peq[0] >>= 1;
        peq[1] >>= 1;
        peq[2] >>= 1;
        peq[3] >>= 1;
        text_index += 1;
        pattern_index += 1;
        if pattern_index < pn {
            let code = SEQUENCE_NT4_TABLE[pattern[pattern_index as usize] as usize];
            if code < 4 {
                peq[code as usize] |= new_bit;
            }
        }
    }

    let d0 = edit_core_64(
        &peq,
        &mut vertical_positive,
        &mut vertical_negative,
        text[last_text_index as usize],
    );
    if d0.diagonal_zero & 1 == 0 {
        error += 1;
        if error > cutoff {
            return;
        }
    }

    let site = tn - 1;
    let available = pn - tn;
    let mut ungapped_error = i32::MAX;
    if error <= threshold && error <= result.error {
        result.error = error;
        result.pattern_end = site;
    }

    let mut index = 0i32;
    while index < available {
        error += ((vertical_positive >> index) & 1) as i32;
        error -= ((vertical_negative >> index) & 1) as i32;
        index += 1;
        if error <= threshold && error <= result.error {
            result.error = error;
            result.pattern_end = site + index;
        }
        if index == threshold {
            ungapped_error = error;
        }
    }

    if ungapped_error <= threshold && ungapped_error == result.error {
        result.pattern_end = site + threshold;
    }
}

/// `ed_band_cal_semi_64_w_absent_diag`: semi-global alignment that starts on
/// the anti-diagonal `abs_diagonal` of the pattern.
///
/// This is the variant used by the overlap code, where a chain already knows
/// how many bases of the pattern were skipped.
pub fn semi_global_64_absent_diagonal(
    pattern: &[u8],
    pattern_length: i32,
    text: &[u8],
    text_length: i32,
    threshold: i32,
    absolute_diagonal: i32,
    result: &mut EditDistanceResult,
) {
    result.reset(threshold, pattern_length, text_length);
    result.pattern_start = -1;
    result.pattern_end = -1;
    result.text_start = 0;
    result.text_end = text_length - 1;

    let pn = pattern_length;
    let tn = text_length;
    let mut error = absolute_diagonal;
    let cutoff = threshold + (threshold << 1);
    if pn > tn + cutoff || tn > pn + cutoff {
        return;
    }

    // `w_64_set_bit_lsub(ez->VN, abs_diag)` sets every bit below `abs_diag`.
    let mut vertical_positive = 0u64;
    let mut vertical_negative = (1u64 << absolute_diagonal) - 1;
    let mut peq = [0u64; 5];

    let mut first_length = (threshold << 1) + 1 - absolute_diagonal;
    if first_length > pn {
        first_length = pn;
    }
    let mut mask = 1u64 << absolute_diagonal;
    for offset in 0..first_length {
        let code = SEQUENCE_NT4_TABLE[pattern[offset as usize] as usize];
        peq[code as usize] |= mask;
        mask <<= 1;
    }
    let mut pattern_index = (threshold << 1) - absolute_diagonal;
    peq[4] = 0;

    let last_text_index = tn - 1;
    let last_high = threshold << 1;
    let new_bit = 1u64 << last_high;

    let mut text_index = 0i32;
    while text_index < last_text_index {
        let d0 = edit_core_64(
            &peq,
            &mut vertical_positive,
            &mut vertical_negative,
            text[text_index as usize],
        );
        if d0.diagonal_zero & 1 == 0 {
            error += 1;
            if error > cutoff {
                return;
            }
        }

        peq[0] >>= 1;
        peq[1] >>= 1;
        peq[2] >>= 1;
        peq[3] >>= 1;
        text_index += 1;
        pattern_index += 1;
        if pattern_index < pn {
            let code = SEQUENCE_NT4_TABLE[pattern[pattern_index as usize] as usize];
            if code < 4 {
                peq[code as usize] |= new_bit;
            }
        }
    }

    let d0 = edit_core_64(
        &peq,
        &mut vertical_positive,
        &mut vertical_negative,
        text[last_text_index as usize],
    );
    if d0.diagonal_zero & 1 == 0 {
        error += 1;
        if error > cutoff {
            return;
        }
    }

    let mut site = tn - 1 - absolute_diagonal;
    let available = pn - tn + absolute_diagonal;
    let mut ungapped_error = i32::MAX;

    let mut index = 0i32;
    while site < 0 && index < available {
        error += ((vertical_positive >> index) & 1) as i32;
        error -= ((vertical_negative >> index) & 1) as i32;
        index += 1;
        site += 1;
    }
    if error <= threshold && error <= result.error {
        result.error = error;
        result.pattern_end = site;
    }

    site -= index;
    while index < available {
        error += ((vertical_positive >> index) & 1) as i32;
        error -= ((vertical_negative >> index) & 1) as i32;
        index += 1;
        if error <= threshold && error <= result.error {
            result.error = error;
            result.pattern_end = site + index;
        }
        if index == threshold {
            ungapped_error = error;
        }
    }

    if ungapped_error <= threshold && ungapped_error == result.error {
        result.pattern_end = site + threshold;
    }
}

/// `ed_band_cal_extension_64_0_w`: extend the alignment forward from the
/// start of both sequences.
pub fn extension_forward_64(
    pattern: &[u8],
    pattern_length: i32,
    text: &[u8],
    text_length: i32,
    threshold: i32,
    result: &mut EditDistanceResult,
) {
    result.reset(threshold, pattern_length, text_length);
    result.pattern_start = 0;
    result.text_start = 0;
    result.pattern_end = -1;
    result.text_end = -1;

    let (pattern_length, text_length) = clamp_extensions(pattern_length, text_length, threshold);
    let pattern_end = pattern_length - 1;
    let last_text_index = text_length - 1;
    let cutoff = threshold + (threshold << 1);

    let first_length = (threshold + 1).min(pattern_length);
    let mut peq = build_pattern_vectors(pattern, first_length, threshold);

    let mut error = threshold;
    let mut vertical_negative = (1u64 << threshold) - 1;
    let mut vertical_positive = (1u64 << ((threshold << 1) + 1)) - 1;
    vertical_positive ^= vertical_negative;

    let mut temporary_error = i32::MAX;
    let mut text_index = 0i32;
    let mut pattern_index = threshold;
    let new_bit = 1u64 << (threshold << 1);

    while text_index < last_text_index {
        let vectors = edit_core_64(
            &peq,
            &mut vertical_positive,
            &mut vertical_negative,
            text[text_index as usize],
        );
        if vectors.diagonal_zero & 1 == 0 {
            error += 1;
            if error > cutoff {
                return;
            }
        }

        let mut offset = text_index - threshold;
        let mut shift = text_index + threshold - pattern_end;
        if shift >= 0 {
            if temporary_error == i32::MAX {
                temporary_error = error;
                let mut position = 0i32;
                while offset < pattern_end {
                    temporary_error += ((vertical_positive >> position) & 1) as i32;
                    temporary_error -= ((vertical_negative >> position) & 1) as i32;
                    offset += 1;
                    position += 1;
                }
            } else {
                shift = (threshold << 1) - shift;
                if shift >= 0 {
                    temporary_error += ((vectors.horizontal_positive >> shift) & 1) as i32;
                    temporary_error -= ((vectors.horizontal_negative >> shift) & 1) as i32;
                }
            }
            if temporary_error <= threshold && temporary_error < result.error {
                result.error = temporary_error;
                result.pattern_end = pattern_end;
                result.text_end = text_index;
            }
        }

        peq[0] >>= 1;
        peq[1] >>= 1;
        peq[2] >>= 1;
        peq[3] >>= 1;
        text_index += 1;
        pattern_index += 1;
        if pattern_index < pattern_length {
            let code = SEQUENCE_NT4_TABLE[pattern[pattern_index as usize] as usize];
            if code < 4 {
                peq[code as usize] |= new_bit;
            }
        }
    }

    let vectors = edit_core_64(
        &peq,
        &mut vertical_positive,
        &mut vertical_negative,
        text[last_text_index as usize],
    );
    if vectors.diagonal_zero & 1 == 0 {
        error += 1;
        if error > cutoff {
            return;
        }
    }

    let mut site = text_length - 1 - threshold;
    let final_site = pattern_length - 1;
    while site < final_site {
        error += (vertical_positive & 1) as i32;
        vertical_positive >>= 1;
        error -= (vertical_negative & 1) as i32;
        vertical_negative >>= 1;
        site += 1;
        if error <= threshold && error < result.error {
            result.error = error;
            result.pattern_end = site;
            result.text_end = text_length - 1;
        }
    }
    if error <= threshold && error < result.error {
        result.error = error;
        result.pattern_end = site;
        result.text_end = text_length - 1;
    }
}

/// `ed_band_cal_extension_64_1_w`: extend the alignment backwards from the end
/// of both sequences.
pub fn extension_backward_64(
    pattern: &[u8],
    pattern_length: i32,
    text: &[u8],
    text_length: i32,
    threshold: i32,
    result: &mut EditDistanceResult,
) {
    result.reset(threshold, pattern_length, text_length);
    result.pattern_start = i32::MAX;
    result.text_start = i32::MAX;
    result.pattern_end = pattern_length - 1;
    result.text_end = text_length - 1;

    // The pivot indices are taken before the length clamp, exactly like the
    // original, which stores them in `ez` and then clamps local variables.
    let pattern_pivot = result.pattern_end;
    let text_pivot = result.text_end;

    let (pattern_length, text_length) = clamp_extensions(pattern_length, text_length, threshold);
    let pattern_end = pattern_length - 1;
    let last_text_index = text_length - 1;
    let cutoff = threshold + (threshold << 1);

    let first_length = (threshold + 1).min(pattern_length);
    let mut peq = [0u64; 5];
    let mut mask = 1u64 << threshold;
    for offset in 0..first_length {
        let code = SEQUENCE_NT4_TABLE[pattern[(pattern_pivot - offset) as usize] as usize];
        peq[code as usize] |= mask;
        mask <<= 1;
    }
    peq[4] = 0;

    let mut error = threshold;
    let mut vertical_negative = (1u64 << threshold) - 1;
    let mut vertical_positive = (1u64 << ((threshold << 1) + 1)) - 1;
    vertical_positive ^= vertical_negative;

    let mut temporary_error = i32::MAX;
    let mut text_index = 0i32;
    let mut pattern_index = threshold;
    let new_bit = 1u64 << (threshold << 1);

    while text_index < last_text_index {
        let vectors = edit_core_64(
            &peq,
            &mut vertical_positive,
            &mut vertical_negative,
            text[(text_pivot - text_index) as usize],
        );
        if vectors.diagonal_zero & 1 == 0 {
            error += 1;
            if error > cutoff {
                return;
            }
        }

        let mut offset = text_index - threshold;
        let mut shift = text_index + threshold - pattern_end;
        if shift >= 0 {
            if temporary_error == i32::MAX {
                temporary_error = error;
                let mut position = 0i32;
                while offset < pattern_end {
                    temporary_error += ((vertical_positive >> position) & 1) as i32;
                    temporary_error -= ((vertical_negative >> position) & 1) as i32;
                    offset += 1;
                    position += 1;
                }
            } else {
                shift = (threshold << 1) - shift;
                if shift >= 0 {
                    temporary_error += ((vectors.horizontal_positive >> shift) & 1) as i32;
                    temporary_error -= ((vectors.horizontal_negative >> shift) & 1) as i32;
                }
            }
            if temporary_error <= threshold && temporary_error < result.error {
                result.error = temporary_error;
                result.pattern_start = pattern_pivot - pattern_end;
                result.text_start = text_pivot - text_index;
            }
        }

        peq[0] >>= 1;
        peq[1] >>= 1;
        peq[2] >>= 1;
        peq[3] >>= 1;
        text_index += 1;
        pattern_index += 1;
        if pattern_index < pattern_length {
            let code =
                SEQUENCE_NT4_TABLE[pattern[(pattern_pivot - pattern_index) as usize] as usize];
            if code < 4 {
                peq[code as usize] |= new_bit;
            }
        }
    }

    let vectors = edit_core_64(
        &peq,
        &mut vertical_positive,
        &mut vertical_negative,
        text[(text_pivot - text_index) as usize],
    );
    if vectors.diagonal_zero & 1 == 0 {
        error += 1;
        if error > cutoff {
            return;
        }
    }

    let mut site = text_length - 1 - threshold;
    let final_site = pattern_length - 1;
    while site < final_site {
        error += (vertical_positive & 1) as i32;
        vertical_positive >>= 1;
        error -= (vertical_negative & 1) as i32;
        vertical_negative >>= 1;
        site += 1;
        if error <= threshold && error < result.error {
            result.error = error;
            result.pattern_start = pattern_pivot - site;
            result.text_start = text_pivot + 1 - text_length;
        }
    }
    if error <= threshold && error < result.error {
        result.error = error;
        result.pattern_start = pattern_pivot - site;
        result.text_start = text_pivot + 1 - text_length;
    }
}

/// Shortens the longer sequence so that the two lengths differ by at most
/// `threshold`, the guard both extension variants apply.
fn clamp_extensions(pattern_length: i32, text_length: i32, threshold: i32) -> (i32, i32) {
    let mut pattern_length = pattern_length;
    let mut text_length = text_length;
    if pattern_length > text_length + threshold {
        pattern_length = text_length + threshold;
    } else if text_length > pattern_length + threshold {
        text_length = pattern_length + threshold;
    }
    (pattern_length, text_length)
}

/// One traceback operation: `code` is 0 for a match, 1 for a substitution,
/// 2 when the pattern has an extra base and 3 when the text does, and `length`
/// is the run length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TraceOperation {
    pub code: u16,
    pub length: u32,
}

/// Packed traceback, the equivalent of `asg16_v` as used for the CIGAR.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TraceCigar {
    /// Runs packed as `code << 14 | length`, exactly like `push_trace`.
    pub packed: Vec<u16>,
}

impl TraceCigar {
    /// `push_trace`: appends a run, splitting it when it does not fit in the
    /// fourteen bits available for the length.
    pub fn push(&mut self, code: u16, mut length: u32) {
        const MAXIMUM_RUN: u32 = 0x3fff;
        let shifted = code << 14;
        while length >= MAXIMUM_RUN {
            self.packed.push(shifted + MAXIMUM_RUN as u16);
            length -= MAXIMUM_RUN;
        }
        if length > 0 {
            self.packed.push(shifted + length as u16);
        }
    }

    /// Expands the packed runs into `(code, length)` pairs, merging the
    /// consecutive runs of the same operation the way `pop_trace` does.
    pub fn operations(&self) -> Vec<TraceOperation> {
        let mut operations: Vec<TraceOperation> = Vec::new();
        for entry in &self.packed {
            let code = entry >> 14;
            let length = (*entry & 0x3fff) as u32;
            match operations.last_mut() {
                Some(last) if last.code == code => last.length += length,
                _ => operations.push(TraceOperation { code, length }),
            }
        }
        operations
    }

    pub fn is_empty(&self) -> bool {
        self.packed.is_empty()
    }
}

/// `gen_trace`: walks the stored bit vectors backwards and rebuilds the CIGAR.
///
/// `trim` is `thre` for the global and extension variants and the absolute
/// diagonal for the semi-global variant. Returns the updated pattern start.
#[allow(clippy::too_many_arguments)]
fn generate_trace(
    path: &[u64],
    error: i32,
    threshold: i32,
    pattern_end: i32,
    text_end: i32,
    text_start: i32,
    pattern_start: i32,
    pattern_length: i32,
    trim: i32,
    reverse: bool,
) -> (TraceCigar, i32) {
    let mut cigar = TraceCigar::default();
    if error > threshold {
        return (cigar, pattern_start);
    }

    let text_length = text_end + 1 - text_start;
    let extended_pattern_length = text_length + (threshold << 1);
    let band_width = (threshold << 1) + 1;
    let column_size = path.len() as i32 / text_length;
    let vector_size = column_size / 5;

    let mut pattern_offset = pattern_end;
    let mut diagonal = band_width - (extended_pattern_length - pattern_end - trim);

    let mut text_index = text_length;
    let lowest_diagonal = band_width - 1;
    let mut direction: i32;
    let mut previous_direction: i32 = -1;
    let mut run_length: i32 = 0;
    let mut current_error = error;

    while text_index > 0 && current_error > 0 {
        let column = text_index - 1;
        let base = column as usize * column_size as usize;
        let word_of = |bit: i32| base + (bit >> 6) as usize;
        let shift_of = |bit: i32| (bit & 63) as usize;
        let vector = |index: usize, bit: i32| path[word_of(bit) + index * vector_size as usize];

        let diagonal_zero = vector(0, diagonal);
        let horizontal_positive = vector(3, diagonal);
        let horizontal_negative = vector(4, diagonal);

        let diagonal_shift = shift_of(diagonal);
        let diagonal_error = current_error - (((!(diagonal_zero >> diagonal_shift)) & 1) as i32);
        direction = 0;
        let mut minimum = diagonal_error;

        if diagonal != lowest_diagonal {
            let horizontal_error = current_error
                + (((horizontal_negative >> diagonal_shift) & 1) as i32)
                - (((horizontal_positive >> diagonal_shift) & 1) as i32);
            if horizontal_error + 1 == current_error && horizontal_error <= minimum {
                minimum = horizontal_error;
                direction = 3;
            }
        }
        if diagonal != 0 {
            let vertical_bit = diagonal - 1;
            let vertical_shift = shift_of(vertical_bit);
            let vertical_positive = vector(1, vertical_bit);
            let vertical_negative = vector(2, vertical_bit);
            let vertical_error = current_error + (((vertical_negative >> vertical_shift) & 1) as i32)
                - (((vertical_positive >> vertical_shift) & 1) as i32);
            if vertical_error + 1 == current_error && vertical_error <= minimum {
                minimum = vertical_error;
                direction = 2;
            }
        }

        if direction == 0 {
            if diagonal_error != current_error {
                direction = 1;
            }
            let _ = direction;
            text_index -= 1;
            pattern_offset -= 1;
        } else if direction == 2 {
            diagonal -= 1;
            pattern_offset -= 1;
        } else {
            text_index -= 1;
            diagonal += 1;
        }

        if direction == previous_direction {
            run_length += 1;
        } else {
            if run_length > 0 {
                cigar.push(previous_direction as u16, run_length as u32);
            }
            previous_direction = direction;
            run_length = 1;
        }
        current_error = minimum;
    }

    if text_index > 0 {
        direction = 0;
        pattern_offset -= text_index;
        if direction == previous_direction {
            run_length += text_index;
        } else {
            if run_length > 0 {
                cigar.push(previous_direction as u16, run_length as u32);
            }
            previous_direction = direction;
            run_length = text_index;
        }
    }

    pattern_offset += 1;
    let mut updated_pattern_start = pattern_start;
    if pattern_start < 0 || pattern_start >= pattern_length {
        updated_pattern_start = pattern_offset;
    } else if pattern_offset > pattern_start {
        direction = 2;
        let extra = pattern_offset - pattern_start;
        if direction == previous_direction {
            run_length += extra;
        } else {
            if run_length > 0 {
                cigar.push(previous_direction as u16, run_length as u32);
            }
            previous_direction = direction;
            run_length = extra;
        }
    }
    if run_length > 0 {
        cigar.push(previous_direction as u16, run_length as u32);
    }

    if reverse {
        cigar.packed.reverse();
    }

    (cigar, updated_pattern_start)
}

/// `ed_band_cal_global_64_w_trace`: global alignment that also produces a
/// CIGAR.
///
/// The caller may pass a result produced by an earlier alignment; when it
/// already has an acceptable error the traceback reuses it, exactly like the
/// reference.
pub fn global_64_trace(
    pattern: &[u8],
    pattern_length: i32,
    text: &[u8],
    text_length: i32,
    threshold: i32,
    result: &mut EditDistanceResult,
    cigar: &mut TraceCigar,
) {
    cigar.packed.clear();

    if result.error > threshold {
        result.reset(threshold, pattern_length, text_length);
        result.pattern_start = 0;
        result.text_start = 0;
    } else if result.error == 0 {
        cigar.push(0, (result.text_end + 1 - result.text_start) as u32);
        return;
    }

    let pn = pattern_length;
    let tn = text_length;
    if pn > tn + threshold || tn > pn + threshold {
        return;
    }

    let last_text_index = tn - 1;
    let cutoff = threshold + (threshold << 1);
    let first_length = (threshold + 1).min(pn);
    let mut peq = build_pattern_vectors(pattern, first_length, threshold);

    let mut error = threshold;
    let mut vertical_negative = (1u64 << threshold) - 1;
    let mut vertical_positive = (1u64 << ((threshold << 1) + 1)) - 1;
    vertical_positive ^= vertical_negative;

    let mut path: Vec<u64> = Vec::with_capacity((tn * 5) as usize);
    let mut text_index = 0i32;
    let mut pattern_index = threshold;
    let new_bit = 1u64 << (threshold << 1);

    while text_index < last_text_index {
        let vectors = edit_core_64(
            &peq,
            &mut vertical_positive,
            &mut vertical_negative,
            text[text_index as usize],
        );
        if vectors.diagonal_zero & 1 == 0 {
            error += 1;
            if error > cutoff {
                return;
            }
        }

        peq[0] >>= 1;
        peq[1] >>= 1;
        peq[2] >>= 1;
        peq[3] >>= 1;
        text_index += 1;
        pattern_index += 1;
        if pattern_index < pn {
            let code = SEQUENCE_NT4_TABLE[pattern[pattern_index as usize] as usize];
            if code < 4 {
                peq[code as usize] |= new_bit;
            }
        }

        path.push(vectors.diagonal_zero);
        path.push(vertical_positive);
        path.push(vertical_negative);
        path.push(vectors.horizontal_positive);
        path.push(vectors.horizontal_negative);
    }

    let vectors = edit_core_64(
        &peq,
        &mut vertical_positive,
        &mut vertical_negative,
        text[last_text_index as usize],
    );
    if vectors.diagonal_zero & 1 == 0 {
        error += 1;
        if error > cutoff {
            return;
        }
    }
    path.push(vectors.diagonal_zero);
    path.push(vertical_positive);
    path.push(vertical_negative);
    path.push(vectors.horizontal_positive);
    path.push(vectors.horizontal_negative);

    let mut site = tn - 1 - threshold;
    if result.error > threshold {
        let final_site = pn - 1;
        while site < final_site {
            error += (vertical_positive & 1) as i32;
            vertical_positive >>= 1;
            error -= (vertical_negative & 1) as i32;
            vertical_negative >>= 1;
            site += 1;
        }
        if site == final_site && error <= threshold {
            result.error = error;
            result.pattern_end = pn - 1;
            result.text_end = tn - 1;
        }
    }

    let (generated, pattern_start) = generate_trace(
        &path,
        result.error,
        threshold,
        result.pattern_end,
        result.text_end,
        result.text_start,
        result.pattern_start,
        result.pattern_length,
        threshold,
        true,
    );
    *cigar = generated;
    result.pattern_start = pattern_start;
}

/// `ed_band_cal_semi_64_w_absent_diag_trace`: the semi-global variant that
/// also produces a CIGAR.
pub fn semi_global_64_absent_diagonal_trace(
    pattern: &[u8],
    pattern_length: i32,
    text: &[u8],
    text_length: i32,
    threshold: i32,
    absolute_diagonal: i32,
    result: &mut EditDistanceResult,
    cigar: &mut TraceCigar,
) {
    cigar.packed.clear();

    if result.error > threshold {
        result.reset(threshold, pattern_length, text_length);
        result.pattern_start = -1;
        result.pattern_end = -1;
        result.text_start = 0;
        result.text_end = text_length - 1;
    } else if result.error == 0 {
        cigar.push(0, (result.text_end + 1 - result.text_start) as u32);
        result.pattern_start = result.pattern_end - (result.text_end - result.text_start);
        return;
    }

    let pn = pattern_length;
    let tn = text_length;
    let mut error = absolute_diagonal;
    let cutoff = threshold + (threshold << 1);
    if pn > tn + cutoff || tn > pn + cutoff {
        return;
    }

    let mut vertical_positive = 0u64;
    let mut vertical_negative = (1u64 << absolute_diagonal) - 1;
    let mut peq = [0u64; 5];

    let mut first_length = (threshold << 1) + 1 - absolute_diagonal;
    if first_length > pn {
        first_length = pn;
    }
    let mut mask = 1u64 << absolute_diagonal;
    for offset in 0..first_length {
        let code = SEQUENCE_NT4_TABLE[pattern[offset as usize] as usize];
        peq[code as usize] |= mask;
        mask <<= 1;
    }
    let mut pattern_index = (threshold << 1) - absolute_diagonal;
    peq[4] = 0;

    let last_text_index = tn - 1;
    let new_bit = 1u64 << (threshold << 1);

    let mut path: Vec<u64> = Vec::with_capacity((tn * 5) as usize);
    let mut text_index = 0i32;

    while text_index < last_text_index {
        let vectors = edit_core_64(
            &peq,
            &mut vertical_positive,
            &mut vertical_negative,
            text[text_index as usize],
        );
        if vectors.diagonal_zero & 1 == 0 {
            error += 1;
            if error > cutoff {
                return;
            }
        }

        peq[0] >>= 1;
        peq[1] >>= 1;
        peq[2] >>= 1;
        peq[3] >>= 1;
        text_index += 1;
        pattern_index += 1;
        if pattern_index < pn {
            let code = SEQUENCE_NT4_TABLE[pattern[pattern_index as usize] as usize];
            if code < 4 {
                peq[code as usize] |= new_bit;
            }
        }

        path.push(vectors.diagonal_zero);
        path.push(vertical_positive);
        path.push(vertical_negative);
        path.push(vectors.horizontal_positive);
        path.push(vectors.horizontal_negative);
    }

    let vectors = edit_core_64(
        &peq,
        &mut vertical_positive,
        &mut vertical_negative,
        text[last_text_index as usize],
    );
    if vectors.diagonal_zero & 1 == 0 {
        error += 1;
        if error > cutoff {
            return;
        }
    }
    path.push(vectors.diagonal_zero);
    path.push(vertical_positive);
    path.push(vertical_negative);
    path.push(vectors.horizontal_positive);
    path.push(vectors.horizontal_negative);

    let mut site = tn - 1 - absolute_diagonal;
    let available = pn - tn + absolute_diagonal;
    let mut ungapped_error = i32::MAX;

    if result.error > threshold {
        let mut index = 0i32;
        while site < 0 && index < available {
            error += ((vertical_positive >> index) & 1) as i32;
            error -= ((vertical_negative >> index) & 1) as i32;
            index += 1;
            site += 1;
        }
        if error <= threshold && error <= result.error {
            result.error = error;
            result.pattern_end = site;
        }

        site -= index;
        while index < available {
            error += ((vertical_positive >> index) & 1) as i32;
            error -= ((vertical_negative >> index) & 1) as i32;
            index += 1;
            if error <= threshold && error <= result.error {
                result.error = error;
                result.pattern_end = site + index;
            }
            if index == threshold {
                ungapped_error = error;
            }
        }

        if ungapped_error <= threshold && ungapped_error == result.error {
            result.pattern_end = site + threshold;
        }
    }

    let (generated, pattern_start) = generate_trace(
        &path,
        result.error,
        threshold,
        result.pattern_end,
        result.text_end,
        result.text_start,
        result.pattern_start,
        result.pattern_length,
        absolute_diagonal,
        true,
    );
    *cigar = generated;
    result.pattern_start = pattern_start;
}

/// Number of 64-bit words a band of `threshold` errors spans, the `nword` the
/// reference computes before it dispatches to a kernel of matching width.
pub fn band_word_count(threshold: i32) -> usize {
    let band = (threshold << 1) + 1;
    (band >> 6) as usize + usize::from((band & 63) != 0)
}

/// `w_infi_set_bit_lsub`: clears the vector and sets its lowest `length` bits.
fn set_low_bits(vector: &mut [u64], length: i32) {
    vector.fill(0);
    let length = length as usize;
    let full_words = (length >> 6).min(vector.len());
    vector[..full_words].fill(u64::MAX);
    let remainder = length & 63;
    if remainder != 0 && full_words < vector.len() {
        vector[full_words] = (1u64 << remainder) - 1;
    }
}

/// Reads one bit of a multi-word vector.
fn vector_bit(vector: &[u64], bit: i32) -> i32 {
    ((vector[(bit >> 6) as usize] >> (bit & 63)) & 1) as i32
}

/// The wide recurrence state, the `w_sig` vectors of `bit_extz_t` each
/// `word_count` words wide: the five pattern vectors (the fifth stays zero,
/// matching `Peq[4]`) and the two scratch and four mask vectors.
struct BandedVectors {
    word_count: usize,
    pattern: [Vec<u64>; 5],
    vertical_positive: Vec<u64>,
    vertical_negative: Vec<u64>,
    diagonal_zero: Vec<u64>,
    horizontal_negative: Vec<u64>,
    horizontal_positive: Vec<u64>,
}

impl BandedVectors {
    fn new(word_count: usize) -> Self {
        let zero = || vec![0u64; word_count];
        Self {
            word_count,
            pattern: std::array::from_fn(|_| zero()),
            vertical_positive: zero(),
            vertical_negative: zero(),
            diagonal_zero: zero(),
            horizontal_negative: zero(),
            horizontal_positive: zero(),
        }
    }

    fn clear_pattern(&mut self) {
        for vector in &mut self.pattern {
            vector.fill(0);
        }
    }

    /// `ed_init_core`: sets `pattern[offset]` at bit `first_bit + offset`.
    fn fill_pattern(&mut self, sequence: &[u8], length: i32, first_bit: i32) {
        self.clear_pattern();
        for offset in 0..length {
            let code = SEQUENCE_NT4_TABLE[sequence[offset as usize] as usize] as usize;
            if code < 4 {
                let bit = first_bit + offset;
                self.pattern[code][(bit >> 6) as usize] |= 1u64 << (bit & 63);
            }
        }
    }

    /// The mirrored `ed_init_core` of the backward extension, which walks the
    /// pattern down from `pivot`.
    fn fill_pattern_backward(&mut self, sequence: &[u8], pivot: i32, length: i32, first_bit: i32) {
        self.clear_pattern();
        for offset in 0..length {
            let code =
                SEQUENCE_NT4_TABLE[sequence[(pivot - offset) as usize] as usize] as usize;
            if code < 4 {
                let bit = first_bit + offset;
                self.pattern[code][(bit >> 6) as usize] |= 1u64 << (bit & 63);
            }
        }
    }

    /// `ed_infi_core`: one Myers step over a single text base.
    fn core(&mut self, text_base: u8) {
        let code = SEQUENCE_NT4_TABLE[text_base as usize] as usize;
        let mut carry = 0u64;
        for word in 0..self.word_count {
            let x = self.pattern[code][word] | self.vertical_negative[word];
            let mut sum = x & self.vertical_positive[word];
            sum = sum.wrapping_add(carry);
            carry = u64::from(sum < carry);
            sum = sum.wrapping_add(self.vertical_positive[word]);
            carry |= u64::from(sum < self.vertical_positive[word]);
            let diagonal = (sum ^ self.vertical_positive[word]) | x;
            self.diagonal_zero[word] = diagonal;
            let horizontal_negative = self.vertical_positive[word] & diagonal;
            self.horizontal_negative[word] = horizontal_negative;
            self.horizontal_positive[word] =
                self.vertical_negative[word] | !(self.vertical_positive[word] | diagonal);
        }
        carry = 0;
        for word in (0..self.word_count).rev() {
            let x = (self.diagonal_zero[word] >> 1) | carry;
            carry = self.diagonal_zero[word] << 63;
            self.vertical_negative[word] = x & self.horizontal_positive[word];
            self.vertical_positive[word] =
                self.horizontal_negative[word] | !(x | self.horizontal_positive[word]);
        }
    }

    /// `ed_infi_post_Peq`: shifts the four pattern vectors right by one bit,
    /// carrying bits down from the higher words.
    fn shift_pattern(&mut self) {
        let last = self.word_count - 1;
        for code in 0..4 {
            for word in 0..last {
                self.pattern[code][word] =
                    (self.pattern[code][word] >> 1) | (self.pattern[code][word + 1] << 63);
            }
            self.pattern[code][last] >>= 1;
        }
    }

    /// Appends the column to the path in the order `gen_trace` expects.
    fn push_path(&self, path: &mut Vec<u64>) {
        path.extend_from_slice(&self.diagonal_zero);
        path.extend_from_slice(&self.vertical_positive);
        path.extend_from_slice(&self.vertical_negative);
        path.extend_from_slice(&self.horizontal_positive);
        path.extend_from_slice(&self.horizontal_negative);
    }
}

/// `VN`/`VP` seed of the global and extension variants: the lowest `thre` bits
/// and the lowest `2*thre+1` bits exclusive-ored together.
fn initialise_seeded_masks(vectors: &mut BandedVectors, threshold: i32) {
    set_low_bits(&mut vectors.vertical_negative, threshold);
    set_low_bits(&mut vectors.vertical_positive, (threshold << 1) + 1);
    for word in 0..vectors.word_count {
        vectors.vertical_positive[word] ^= vectors.vertical_negative[word];
    }
}

/// Width-generic global alignment with traceback, the counterpart of
/// `ed_band_cal_global_infi_w_trace`. It stands in for the 64-bit and 128-bit
/// kernels as well, which perform the same arithmetic on the band bits.
pub fn global_trace(
    pattern: &[u8],
    pattern_length: i32,
    text: &[u8],
    text_length: i32,
    threshold: i32,
    result: &mut EditDistanceResult,
    cigar: &mut TraceCigar,
) {
    cigar.packed.clear();

    if result.error > threshold {
        result.reset(threshold, pattern_length, text_length);
        result.pattern_start = 0;
        result.text_start = 0;
    } else if result.error == 0 {
        cigar.push(0, (result.text_end + 1 - result.text_start) as u32);
        return;
    }

    let pn = pattern_length;
    let tn = text_length;
    if pn > tn + threshold || tn > pn + threshold {
        return;
    }

    let mut vectors = BandedVectors::new(band_word_count(threshold));
    let first_length = (threshold + 1).min(pn);
    vectors.fill_pattern(pattern, first_length, threshold);
    let mut error = threshold;
    initialise_seeded_masks(&mut vectors, threshold);

    let last_text_index = tn - 1;
    let cutoff = threshold + (threshold << 1);
    let new_bit = threshold << 1;
    let mut pattern_index = threshold;
    let mut text_index = 0i32;
    let mut path: Vec<u64> = Vec::with_capacity((vectors.word_count * 5 * tn as usize) as usize);

    while text_index < last_text_index {
        vectors.core(text[text_index as usize]);
        if vectors.diagonal_zero[0] & 1 == 0 {
            error += 1;
            if error > cutoff {
                return;
            }
        }

        vectors.shift_pattern();
        text_index += 1;
        pattern_index += 1;
        if pattern_index < pn {
            let code = SEQUENCE_NT4_TABLE[pattern[pattern_index as usize] as usize] as usize;
            if code < 4 {
                vectors.pattern[code][(new_bit >> 6) as usize] |= 1u64 << (new_bit & 63);
            }
        }
        vectors.push_path(&mut path);
    }

    vectors.core(text[last_text_index as usize]);
    if vectors.diagonal_zero[0] & 1 == 0 {
        error += 1;
        if error > cutoff {
            return;
        }
    }
    vectors.push_path(&mut path);

    let mut site = tn - 1 - threshold;
    if result.error > threshold {
        let final_site = pn - 1;
        let mut index = 0i32;
        while site < final_site {
            error += vector_bit(&vectors.vertical_positive, index);
            error -= vector_bit(&vectors.vertical_negative, index);
            index += 1;
            site += 1;
        }
        if site == final_site && error <= threshold {
            result.error = error;
            result.pattern_end = pn - 1;
            result.text_end = tn - 1;
        }
    }

    let (generated, pattern_start) = generate_trace(
        &path,
        result.error,
        threshold,
        result.pattern_end,
        result.text_end,
        result.text_start,
        result.pattern_start,
        result.pattern_length,
        threshold,
        true,
    );
    *cigar = generated;
    result.pattern_start = pattern_start;
}

/// Width-generic semi-global alignment on an absent anti-diagonal with
/// traceback, the counterpart of `ed_band_cal_semi_infi_w_absent_diag_trace`.
pub fn semi_global_absent_diagonal_trace(
    pattern: &[u8],
    pattern_length: i32,
    text: &[u8],
    text_length: i32,
    threshold: i32,
    absolute_diagonal: i32,
    result: &mut EditDistanceResult,
    cigar: &mut TraceCigar,
) {
    cigar.packed.clear();

    if result.error > threshold {
        result.reset(threshold, pattern_length, text_length);
        result.pattern_start = -1;
        result.pattern_end = -1;
        result.text_start = 0;
        result.text_end = text_length - 1;
    } else if result.error == 0 {
        cigar.push(0, (result.text_end + 1 - result.text_start) as u32);
        result.pattern_start = result.pattern_end - (result.text_end - result.text_start);
        return;
    }

    let pn = pattern_length;
    let tn = text_length;
    let cutoff = threshold + (threshold << 1);
    if pn > tn + cutoff || tn > pn + cutoff {
        return;
    }

    let mut vectors = BandedVectors::new(band_word_count(threshold));
    let first_length = ((threshold << 1) + 1 - absolute_diagonal).min(pn);
    vectors.fill_pattern(pattern, first_length, absolute_diagonal);
    set_low_bits(&mut vectors.vertical_negative, absolute_diagonal);

    let mut error = absolute_diagonal;
    let last_text_index = tn - 1;
    let new_bit = threshold << 1;
    let mut pattern_index = (threshold << 1) - absolute_diagonal;
    let mut text_index = 0i32;
    let mut path: Vec<u64> = Vec::with_capacity((vectors.word_count * 5 * tn as usize) as usize);

    while text_index < last_text_index {
        vectors.core(text[text_index as usize]);
        if vectors.diagonal_zero[0] & 1 == 0 {
            error += 1;
            if error > cutoff {
                return;
            }
        }

        vectors.shift_pattern();
        text_index += 1;
        pattern_index += 1;
        if pattern_index < pn {
            let code = SEQUENCE_NT4_TABLE[pattern[pattern_index as usize] as usize] as usize;
            if code < 4 {
                vectors.pattern[code][(new_bit >> 6) as usize] |= 1u64 << (new_bit & 63);
            }
        }
        vectors.push_path(&mut path);
    }

    vectors.core(text[last_text_index as usize]);
    if vectors.diagonal_zero[0] & 1 == 0 {
        error += 1;
        if error > cutoff {
            return;
        }
    }
    vectors.push_path(&mut path);

    let mut site = tn - 1 - absolute_diagonal;
    let available = pn - tn + absolute_diagonal;
    let mut ungapped_error = i32::MAX;
    if result.error > threshold {
        let mut index = 0i32;
        while site < 0 && index < available {
            error += vector_bit(&vectors.vertical_positive, index);
            error -= vector_bit(&vectors.vertical_negative, index);
            index += 1;
            site += 1;
        }
        if error <= threshold && error <= result.error {
            result.error = error;
            result.pattern_end = site;
        }

        site -= index;
        while index < available {
            error += vector_bit(&vectors.vertical_positive, index);
            error -= vector_bit(&vectors.vertical_negative, index);
            index += 1;
            if error <= threshold && error <= result.error {
                result.error = error;
                result.pattern_end = site + index;
            }
            if index == threshold {
                ungapped_error = error;
            }
        }

        if ungapped_error <= threshold && ungapped_error == result.error {
            result.pattern_end = site + threshold;
        }
    }

    let (generated, pattern_start) = generate_trace(
        &path,
        result.error,
        threshold,
        result.pattern_end,
        result.text_end,
        result.text_start,
        result.pattern_start,
        result.pattern_length,
        absolute_diagonal,
        true,
    );
    *cigar = generated;
    result.pattern_start = pattern_start;
}

/// Width-generic forward extension with traceback, the counterpart of
/// `ed_band_cal_extension_infi_0_w_trace`.
pub fn extension_forward_trace(
    pattern: &[u8],
    pattern_length: i32,
    text: &[u8],
    text_length: i32,
    threshold: i32,
    result: &mut EditDistanceResult,
    cigar: &mut TraceCigar,
) {
    cigar.packed.clear();

    let mut done = true;
    let mut pn = pattern_length;
    let mut tn = text_length;
    if result.error > threshold {
        result.reset(threshold, pattern_length, text_length);
        result.pattern_start = 0;
        result.text_start = 0;
        result.pattern_end = -1;
        result.text_end = -1;
        if pn > tn + threshold {
            pn = tn + threshold;
        } else if tn > pn + threshold {
            tn = pn + threshold;
        }
        done = false;
    } else if result.error == 0 {
        cigar.push(0, (result.text_end + 1 - result.text_start) as u32);
        return;
    } else {
        pn = result.pattern_end + 1 - result.pattern_start;
        tn = result.text_end + 1 - result.text_start;
    }

    let mut vectors = BandedVectors::new(band_word_count(threshold));
    let first_length = (threshold + 1).min(pn);
    vectors.fill_pattern(pattern, first_length, threshold);
    let mut error = threshold;
    initialise_seeded_masks(&mut vectors, threshold);

    let last_text_index = tn - 1;
    let cutoff = threshold + (threshold << 1);
    let new_bit = threshold << 1;
    let pattern_end = pn - 1;
    let mut temporary_error = i32::MAX;
    let mut pattern_index = threshold;
    let mut text_index = 0i32;
    let mut path: Vec<u64> = Vec::with_capacity((vectors.word_count * 5 * tn as usize) as usize);

    while text_index < last_text_index {
        vectors.core(text[text_index as usize]);
        if vectors.diagonal_zero[0] & 1 == 0 {
            error += 1;
            if error > cutoff {
                return;
            }
        }

        if !done {
            let mut pattern_offset = text_index - threshold;
            let mut shift = text_index + threshold - pattern_end;
            if shift >= 0 {
                if temporary_error == i32::MAX {
                    temporary_error = error;
                    let mut position = 0i32;
                    while pattern_offset < pattern_end {
                        temporary_error += vector_bit(&vectors.vertical_positive, position);
                        temporary_error -= vector_bit(&vectors.vertical_negative, position);
                        pattern_offset += 1;
                        position += 1;
                    }
                } else {
                    shift = (threshold << 1) - shift;
                    if shift >= 0 {
                        temporary_error += vector_bit(&vectors.horizontal_positive, shift);
                        temporary_error -= vector_bit(&vectors.horizontal_negative, shift);
                    }
                }
                if temporary_error <= threshold && temporary_error < result.error {
                    result.error = temporary_error;
                    result.pattern_end = pattern_end;
                    result.text_end = text_index;
                }
            }
        }

        vectors.shift_pattern();
        text_index += 1;
        pattern_index += 1;
        if pattern_index < pn {
            let code = SEQUENCE_NT4_TABLE[pattern[pattern_index as usize] as usize] as usize;
            if code < 4 {
                vectors.pattern[code][(new_bit >> 6) as usize] |= 1u64 << (new_bit & 63);
            }
        }
        vectors.push_path(&mut path);
    }

    vectors.core(text[last_text_index as usize]);
    if vectors.diagonal_zero[0] & 1 == 0 {
        error += 1;
        if error > cutoff {
            return;
        }
    }
    vectors.push_path(&mut path);

    let mut site = tn - 1 - threshold;
    if !done {
        let mut index = 0i32;
        let final_site = pn - 1;
        while site < final_site {
            error += vector_bit(&vectors.vertical_positive, index);
            error -= vector_bit(&vectors.vertical_negative, index);
            index += 1;
            site += 1;
            if error <= threshold && error < result.error {
                result.error = error;
                result.pattern_end = site;
                result.text_end = tn - 1;
            }
        }
        if error <= threshold && error < result.error {
            result.error = error;
            result.pattern_end = site;
            result.text_end = tn - 1;
        }
    }

    let traced_text_length = result.text_end + 1 - result.text_start;
    if traced_text_length != tn {
        let column_words = path.len() / tn as usize;
        // The reference multiplies a `size_t` by a possibly negative extent;
        // the result is only ever used while the error already exceeds the
        // threshold, so wrapping keeps the arithmetic faithful.
        path.truncate(column_words.wrapping_mul(traced_text_length as usize));
    }

    let (generated, pattern_start) = generate_trace(
        &path,
        result.error,
        threshold,
        result.pattern_end,
        result.text_end,
        result.text_start,
        result.pattern_start,
        result.pattern_length,
        threshold,
        true,
    );
    *cigar = generated;
    result.pattern_start = pattern_start;
}

/// Width-generic backward extension with traceback, the counterpart of
/// `ed_band_cal_extension_infi_1_w_trace`.
pub fn extension_backward_trace(
    pattern: &[u8],
    pattern_length: i32,
    text: &[u8],
    text_length: i32,
    threshold: i32,
    result: &mut EditDistanceResult,
    cigar: &mut TraceCigar,
) {
    cigar.packed.clear();

    let mut done = true;
    let mut pn = pattern_length;
    let mut tn = text_length;
    if result.error > threshold {
        result.reset(threshold, pattern_length, text_length);
        result.pattern_start = i32::MAX;
        result.text_start = i32::MAX;
        result.pattern_end = pattern_length - 1;
        result.text_end = text_length - 1;
        if pn > tn + threshold {
            pn = tn + threshold;
        } else if tn > pn + threshold {
            tn = pn + threshold;
        }
        done = false;
    } else if result.error == 0 {
        cigar.push(0, (result.text_end + 1 - result.text_start) as u32);
        return;
    } else {
        pn = result.pattern_end + 1 - result.pattern_start;
        tn = result.text_end + 1 - result.text_start;
    }

    // The pivot coordinates are the stored ends, taken after the reset above.
    let pattern_pivot = result.pattern_end;
    let text_pivot = result.text_end;

    let mut vectors = BandedVectors::new(band_word_count(threshold));
    let first_length = (threshold + 1).min(pn);
    vectors.fill_pattern_backward(pattern, pattern_pivot, first_length, threshold);
    let mut error = threshold;
    initialise_seeded_masks(&mut vectors, threshold);

    let last_text_index = tn - 1;
    let cutoff = threshold + (threshold << 1);
    let new_bit = threshold << 1;
    let pattern_end = pn - 1;
    let mut temporary_error = i32::MAX;
    let mut pattern_index = threshold;
    let mut text_index = 0i32;
    let mut path: Vec<u64> = Vec::with_capacity((vectors.word_count * 5 * tn as usize) as usize);

    while text_index < last_text_index {
        vectors.core(text[(text_pivot - text_index) as usize]);
        if vectors.diagonal_zero[0] & 1 == 0 {
            error += 1;
            if error > cutoff {
                return;
            }
        }

        if !done {
            let mut pattern_offset = text_index - threshold;
            let mut shift = text_index + threshold - pattern_end;
            if shift >= 0 {
                if temporary_error == i32::MAX {
                    temporary_error = error;
                    let mut position = 0i32;
                    while pattern_offset < pattern_end {
                        temporary_error += vector_bit(&vectors.vertical_positive, position);
                        temporary_error -= vector_bit(&vectors.vertical_negative, position);
                        pattern_offset += 1;
                        position += 1;
                    }
                } else {
                    shift = (threshold << 1) - shift;
                    if shift >= 0 {
                        temporary_error += vector_bit(&vectors.horizontal_positive, shift);
                        temporary_error -= vector_bit(&vectors.horizontal_negative, shift);
                    }
                }
                if temporary_error <= threshold && temporary_error < result.error {
                    result.error = temporary_error;
                    result.pattern_start = pattern_pivot - pattern_end;
                    result.text_start = text_pivot - text_index;
                }
            }
        }

        vectors.shift_pattern();
        text_index += 1;
        pattern_index += 1;
        if pattern_index < pn {
            let code =
                SEQUENCE_NT4_TABLE[pattern[(pattern_pivot - pattern_index) as usize] as usize]
                    as usize;
            if code < 4 {
                vectors.pattern[code][(new_bit >> 6) as usize] |= 1u64 << (new_bit & 63);
            }
        }
        vectors.push_path(&mut path);
    }

    vectors.core(text[(text_pivot - text_index) as usize]);
    if vectors.diagonal_zero[0] & 1 == 0 {
        error += 1;
        if error > cutoff {
            return;
        }
    }
    vectors.push_path(&mut path);

    let mut site = tn - 1 - threshold;
    if !done {
        let mut index = 0i32;
        let final_site = pn - 1;
        while site < final_site {
            error += vector_bit(&vectors.vertical_positive, index);
            error -= vector_bit(&vectors.vertical_negative, index);
            index += 1;
            site += 1;
            if error <= threshold && error < result.error {
                result.error = error;
                result.pattern_start = pattern_pivot - site;
                result.text_start = text_pivot + 1 - tn;
            }
        }
        if error <= threshold && error < result.error {
            result.error = error;
            result.pattern_start = pattern_pivot - site;
            result.text_start = text_pivot + 1 - tn;
        }
    }

    let traced_text_length = result.text_end + 1 - result.text_start;
    if traced_text_length != tn {
        let column_words = path.len() / tn as usize;
        // See the forward extension: the reference's `size_t` product wraps
        // when the extent is negative, and the result is not used then.
        path.truncate(column_words.wrapping_mul(traced_text_length as usize));
    }

    // The traceback runs on mirrored pattern coordinates, exactly like the
    // reference.
    let saved_start = result.pattern_start;
    result.pattern_start = pattern_pivot - result.pattern_end;
    result.pattern_end = pattern_pivot - saved_start;

    let (generated, pattern_start) = generate_trace(
        &path,
        result.error,
        threshold,
        result.pattern_end,
        result.text_end,
        result.text_start,
        result.pattern_start,
        result.pattern_length,
        threshold,
        false,
    );
    *cigar = generated;
    result.pattern_start = pattern_start;

    let saved_start = result.pattern_start;
    result.pattern_start = pattern_pivot - result.pattern_end;
    result.pattern_end = pattern_pivot - saved_start;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_sequences_have_zero_errors() {
        let pattern = b"ACGTACGTAC";
        let text = b"ACGTACGTAC";
        let mut result = EditDistanceResult::default();
        global_64(pattern, 10, text, 10, 2, &mut result);
        assert_eq!(result.error, 0);
    }

    #[test]
    fn substitutions_are_counted() {
        let pattern = b"ACGTACGTAC";
        let text = b"ACGTTCGTAG";
        let mut result = EditDistanceResult::default();
        global_64(pattern, 10, text, 10, 2, &mut result);
        assert_eq!(result.error, 2);
    }

    #[test]
    fn too_many_errors_exceed_the_threshold() {
        let pattern = b"AAAAAAAAAA";
        let text = b"GGGGGGGGGG";
        let mut result = EditDistanceResult::default();
        global_64(pattern, 10, text, 10, 2, &mut result);
        assert_eq!(result.error, i32::MAX);
    }

    #[test]
    fn sequences_far_outside_the_band_are_rejected() {
        // `text_length > pattern_length + 3 * threshold`, so the reference
        // gives up immediately and leaves the error untouched.
        let pattern = b"ACGTACGTAC";
        let text = b"TTTTACGTACGTACTTTT";
        let mut result = EditDistanceResult::default();
        semi_global_64(pattern, 10, text, 18, 2, &mut result);
        assert_eq!(result.error, i32::MAX);
        assert_eq!(result.pattern_end, -1);
    }

    /// Golden values printed by the reference `Levenshtein_distance.h`.
    #[test]
    fn matches_the_reference_64_bit_implementation() {
        let pattern = b"ACGTACGTACGT";
        let one_mismatch = b"ACGTACGTACGA";

        let mut result = EditDistanceResult::default();
        semi_global_64(pattern, 12, pattern, 12, 2, &mut result);
        assert_eq!((result.error, result.pattern_end), (0, 11));

        semi_global_64(pattern, 12, one_mismatch, 12, 2, &mut result);
        assert_eq!((result.error, result.pattern_end), (1, 11));

        semi_global_64(pattern, 12, b"ACGTACGTACGTAC", 14, 2, &mut result);
        assert_eq!((result.error, result.pattern_end), (2, 13));

        semi_global_64_absent_diagonal(pattern, 12, pattern, 12, 2, 0, &mut result);
        assert_eq!((result.error, result.pattern_end), (0, 11));

        semi_global_64_absent_diagonal(pattern, 12, one_mismatch, 12, 2, 0, &mut result);
        assert_eq!((result.error, result.pattern_end), (1, 11));

        semi_global_64_absent_diagonal(pattern, 12, pattern, 12, 2, 1, &mut result);
        assert_eq!((result.error, result.pattern_end), (0, 11));

        semi_global_64_absent_diagonal(pattern, 12, pattern, 12, 2, 2, &mut result);
        assert_eq!((result.error, result.pattern_end), (0, 11));

        global_64(pattern, 12, pattern, 12, 2, &mut result);
        assert_eq!((result.error, result.pattern_end, result.text_end), (0, 11, 11));

        global_64(pattern, 12, one_mismatch, 12, 2, &mut result);
        assert_eq!((result.error, result.pattern_end, result.text_end), (1, 11, 11));
    }

    /// Golden values printed by the reference `Levenshtein_distance.h` for the
    /// two extension variants.
    #[test]
    fn extensions_match_the_reference_implementation() {
        let pattern = b"ACGTACGTACGT";
        let mismatch = b"ACGTACGTACGA";
        let longer = b"ACGTACGTACGTAC";
        let mut result = EditDistanceResult::default();

        extension_forward_64(pattern, 12, pattern, 12, 2, &mut result);
        assert_eq!((result.error, result.pattern_end, result.text_end), (0, 11, 11));

        extension_forward_64(pattern, 12, mismatch, 12, 2, &mut result);
        assert_eq!((result.error, result.pattern_end, result.text_end), (1, 11, 10));

        extension_forward_64(pattern, 12, longer, 14, 2, &mut result);
        assert_eq!((result.error, result.pattern_end, result.text_end), (0, 11, 11));

        extension_forward_64(longer, 14, pattern, 12, 2, &mut result);
        assert_eq!((result.error, result.pattern_end, result.text_end), (0, 11, 11));

        extension_backward_64(pattern, 12, pattern, 12, 2, &mut result);
        assert_eq!((result.error, result.pattern_start, result.text_start), (0, 0, 0));

        extension_backward_64(pattern, 12, mismatch, 12, 2, &mut result);
        assert_eq!((result.error, result.pattern_start, result.text_start), (1, 0, 0));

        extension_backward_64(pattern, 12, longer, 14, 2, &mut result);
        assert_eq!((result.error, result.pattern_start, result.text_start), (2, 0, 4));

        extension_backward_64(longer, 14, pattern, 12, 2, &mut result);
        assert_eq!((result.error, result.pattern_start, result.text_start), (2, 0, 0));
    }

    fn cigar_text(cigar: &TraceCigar) -> String {
        let letters = ['M', 'S', 'D', 'I'];
        cigar
            .operations()
            .iter()
            .map(|operation| format!("{}:{}", operation.length, letters[operation.code as usize]))
            .collect::<Vec<String>>()
            .join(" ")
    }

    /// Golden values printed by the reference `Levenshtein_distance.h` after
    /// the non-trace call followed by the trace call, which is how the overlap
    /// code drives it.
    #[test]
    fn global_trace_matches_the_reference_implementation() {
        let pattern = b"ACGTACGTACGT";
        let cases: Vec<(&[u8], i32, i32, i32, &str)> = vec![
            (pattern, 12, 2, 0, "12:M"),
            (b"ACGTACGTACGA", 12, 2, 1, "11:M 1:S"),
            (b"ACGTACGTACGTAC", 14, 2, 2, "12:M 2:I"),
            (b"ACGTAACGTACGT", 13, 3, 1, "5:M 1:I 7:M"),
        ];

        for (text, text_length, threshold, expected_error, expected_cigar) in cases {
            let mut result = EditDistanceResult::default();
            global_64(pattern, 12, text, text_length, threshold, &mut result);
            let mut cigar = TraceCigar::default();
            global_64_trace(pattern, 12, text, text_length, threshold, &mut result, &mut cigar);

            assert_eq!(result.error, expected_error, "error of {:?}", text);
            assert_eq!(cigar_text(&cigar), expected_cigar, "cigar of {:?}", text);
            assert_eq!(result.pattern_start, 0);
        }
    }

    /// Golden values for the semi-global trace variant, again obtained by the
    /// reference after the non-trace call.
    #[test]
    fn semi_global_trace_matches_the_reference_implementation() {
        let pattern = b"ACGTACGTACGT";
        let cases: Vec<(&[u8], i32, &[u8], i32, i32, i32, i32, &str, i32, i32)> = vec![
            (pattern, 12, pattern, 12, 2, 0, 0, "12:M", 0, 11),
            (pattern, 12, b"ACGTACGTACGA", 12, 2, 0, 1, "11:M 1:S", 0, 11),
            (pattern, 12, pattern, 12, 2, 1, 0, "12:M", 0, 11),
            (b"ACGTACGTACGTAC", 14, pattern, 12, 2, 0, 0, "12:M", 0, 11),
        ];

        for (
            case_pattern,
            pattern_length,
            text,
            text_length,
            threshold,
            diagonal,
            expected_error,
            expected_cigar,
            expected_start,
            expected_end,
        ) in cases
        {
            let mut result = EditDistanceResult::default();
            semi_global_64_absent_diagonal(
                case_pattern,
                pattern_length,
                text,
                text_length,
                threshold,
                diagonal,
                &mut result,
            );
            let mut cigar = TraceCigar::default();
            semi_global_64_absent_diagonal_trace(
                case_pattern,
                pattern_length,
                text,
                text_length,
                threshold,
                diagonal,
                &mut result,
                &mut cigar,
            );

            assert_eq!(result.error, expected_error, "error of {:?}", text);
            assert_eq!(cigar_text(&cigar), expected_cigar, "cigar of {:?}", text);
            assert_eq!(result.pattern_start, expected_start);
            assert_eq!(result.pattern_end, expected_end);
        }
    }

    #[test]
    fn absent_diagonal_initialises_the_lower_bits() {
        // A diagonal offset of one leaves the first pattern base unaligned, so
        // the reported end slides with the diagonal.
        let pattern = b"ACGTACGTACGT";
        let mut result = EditDistanceResult::default();
        semi_global_64_absent_diagonal(pattern, 12, pattern, 12, 2, 1, &mut result);
        assert_eq!(result.error, 0);
        assert_eq!(result.pattern_end, 11);
    }
}

/// Replay of the reference golden probe.
///
/// `/tmp/ed_trace_probe.cpp` links against the unmodified hifiasm objects and
/// excercises the four trace kernels the correction code dispatches to, at
/// thresholds that span one, two, three and many words. The generated output
/// is stored under `reference/`; the test rebuilds the same sequences and
/// cases and compares every result field and every packed CIGAR cell.
#[cfg(test)]
mod trace_goldens {
    use super::*;

    struct ProbeCase {
        name: &'static str,
        pattern: Vec<u8>,
        text: Vec<u8>,
    }

    struct GoldenTrace {
        name: String,
        mode: i32,
        error: i32,
        pattern_start: i32,
        pattern_end: i32,
        text_start: i32,
        text_end: i32,
        pattern_length: i32,
        text_length: i32,
        cigar: Vec<(u16, u32)>,
    }

    struct Lcg(u64);

    impl Lcg {
        fn next_u32(&mut self) -> u32 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 33) as u32
        }

        fn sequence(&mut self, length: usize) -> Vec<u8> {
            const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
            (0..length)
                .map(|_| BASES[(self.next_u32() & 3) as usize])
                .collect()
        }
    }

    fn substitute(sequence: &[u8], step: usize) -> Vec<u8> {
        let mut result = sequence.to_vec();
        for index in (0..result.len()).step_by(step) {
            result[index] = if result[index] == b'A' { b'C' } else { b'A' };
        }
        result
    }

    fn probe_cases() -> Vec<ProbeCase> {
        let mut state = Lcg(88172645463325252);
        let core = state.sequence(500);
        let prefix = state.sequence(120);
        let suffix = state.sequence(90);

        let mut insertion = core.clone();
        insertion.splice(200..200, b"ACGTACGTAC".iter().copied());

        let mut deletion = core.clone();
        deletion.drain(180..192);

        let mut with_prefix = prefix.clone();
        with_prefix.extend_from_slice(&core);

        let mut tail = prefix;
        tail.extend_from_slice(&core[150..]);

        let mut longer = core.clone();
        longer.extend_from_slice(&suffix);
        let mut shorter = core.clone();
        shorter.extend_from_slice(&suffix);

        vec![
            ProbeCase { name: "exact", pattern: core.clone(), text: core.clone() },
            ProbeCase { name: "subs50", pattern: core.clone(), text: substitute(&core, 50) },
            ProbeCase { name: "subs7", pattern: core.clone(), text: substitute(&core, 7) },
            ProbeCase { name: "ins10", pattern: core.clone(), text: insertion },
            ProbeCase { name: "del12", pattern: core.clone(), text: deletion },
            ProbeCase {
                name: "shift_plus",
                pattern: core.clone(),
                text: with_prefix.clone(),
            },
            ProbeCase { name: "shift_minus", pattern: with_prefix, text: core.clone() },
            ProbeCase { name: "longer", pattern: core.clone(), text: longer },
            ProbeCase { name: "shorter", pattern: shorter, text: core.clone() },
            ProbeCase { name: "tail_only", pattern: core, text: tail },
        ]
    }

    fn parse_golden(contents: &str) -> Vec<GoldenTrace> {
        let mut records = Vec::new();
        let mut lines = contents.lines();
        while let Some(line) = lines.next() {
            if !line.starts_with("TRACE ") {
                continue;
            }
            let mut tokens = line.split_whitespace().skip(1);
            let name = tokens.next().expect("case name").to_string();
            let mut record = GoldenTrace {
                name,
                mode: 0,
                error: 0,
                pattern_start: 0,
                pattern_end: 0,
                text_start: 0,
                text_end: 0,
                pattern_length: 0,
                text_length: 0,
                cigar: Vec::new(),
            };
            let mut cigar_count = 0usize;
            for token in tokens {
                let (key, value) = token.split_once('=').expect("key=value");
                match key {
                    "mode" => record.mode = value.parse().expect("mode"),
                    "err" => record.error = value.parse().expect("err"),
                    "ps" => record.pattern_start = value.parse().expect("ps"),
                    "pe" => record.pattern_end = value.parse().expect("pe"),
                    "ts" => record.text_start = value.parse().expect("ts"),
                    "te" => record.text_end = value.parse().expect("te"),
                    "pl" => record.pattern_length = value.parse().expect("pl"),
                    "tl" => record.text_length = value.parse().expect("tl"),
                    "cigar_n" => cigar_count = value.parse().expect("cigar_n"),
                    // `path_n` is a `size_t` and can be a wrapped value in
                    // failed extensions; it is not part of the comparison.
                    _ => {}
                }
            }
            let cigar_line = lines.next().expect("cigar line");
            let mut cigar = Vec::new();
            for token in cigar_line.split_whitespace().skip(1) {
                let (length, letter) = token.split_at(token.len() - 1);
                let code = match letter {
                    "E" => 0,
                    "M" => 1,
                    "D" => 2,
                    "I" => 3,
                    other => panic!("unknown cigar letter {}", other),
                };
                cigar.push((code, length.parse::<u32>().expect("run length")));
            }
            assert_eq!(cigar.len(), cigar_count, "cigar of {}", record.name);
            record.cigar = cigar;
            records.push(record);
        }
        records
    }

    #[test]
    #[ignore = "requires reference/ed_trace_probe.out from the C++ probe"]
    fn trace_kernels_match_the_reference_dispatch() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("reference")
            .join("ed_trace_probe.out");
        let contents = std::fs::read_to_string(&path)
            .unwrap_or_else(|_| panic!("missing golden file {}", path.display()));
        let records = parse_golden(&contents);
        let cases = probe_cases();

        let thresholds: [i32; 17] = [
            1, 3, 16, 31, 32, 63, 64, 65, 100, 128, 255, 300, 400, 511, 512, 1000, 2047,
        ];
        let aux_values: [i32; 4] = [0, 1, 5, 15];

        let mut record_index = 0usize;
        for case in &cases {
            for &threshold in &thresholds {
                for mode in 0..4i32 {
                    let mut aux = 0i32;
                    if mode == 3 {
                        aux = aux_values[(((threshold + case.name.len() as i32) >> 2) & 3) as usize];
                        if aux >= threshold {
                            aux = threshold / 2;
                        }
                    }
                    let label = format!("{}|m{}|t{}|a{}", case.name, mode, threshold, aux);
                    let pattern_length = case.pattern.len() as i32;
                    let text_length = case.text.len() as i32;
                    if mode == 3 && pattern_length - text_length + aux > (threshold << 1) {
                        // Unreachable through `cal_exz_infi_adv`: it derives the
                        // semi-global window through `update_semi_coord`, and
                        // `init_waln` bounds `pn - tn + aux_beg` by `2*thre`.
                        // Past that the C kernel indexes `VP`/`VN` beyond
                        // `nword`, so the probe does not record such a call.
                        continue;
                    }

                    let mut result = EditDistanceResult::default();
                    // `clear_align`, the state the caller hands over.
                    result.error = i32::MAX;
                    let mut cigar = TraceCigar::default();
                    match mode {
                        0 => global_trace(
                            &case.pattern, pattern_length, &case.text, text_length,
                            threshold, &mut result, &mut cigar,
                        ),
                        1 => extension_forward_trace(
                            &case.pattern, pattern_length, &case.text, text_length,
                            threshold, &mut result, &mut cigar,
                        ),
                        2 => extension_backward_trace(
                            &case.pattern, pattern_length, &case.text, text_length,
                            threshold, &mut result, &mut cigar,
                        ),
                        _ => semi_global_absent_diagonal_trace(
                            &case.pattern, pattern_length, &case.text, text_length,
                            threshold, aux, &mut result, &mut cigar,
                        ),
                    }

                    let record = &records[record_index];
                    record_index += 1;
                    assert_eq!(record.name, label, "golden order");
                    assert_eq!(record.mode, mode, "mode of {}", label);
                    assert_eq!(record.error, result.error, "error of {}", label);
                    assert_eq!(record.pattern_start, result.pattern_start, "ps of {}", label);
                    assert_eq!(record.pattern_end, result.pattern_end, "pe of {}", label);
                    assert_eq!(record.text_start, result.text_start, "ts of {}", label);
                    assert_eq!(record.text_end, result.text_end, "te of {}", label);
                    assert_eq!(record.pattern_length, result.pattern_length, "pl of {}", label);
                    assert_eq!(record.text_length, result.text_length, "tl of {}", label);
                    let packed: Vec<(u16, u32)> = cigar
                        .packed
                        .iter()
                        .map(|cell| (cell >> 14, (*cell & 0x3fff) as u32))
                        .collect();
                    assert_eq!(record.cigar, packed, "cigar of {}", label);
                }
            }
        }
        assert_eq!(record_index, records.len(), "unconsumed golden records");
    }
}

//! Coverage inference from the k-mer histogram, the counterpart of `hist.cpp`.
//!
//! The routine returns the homozygous peak position (which becomes the read
//! coverage used for every downstream threshold) and reports the heterozygous
//! peak through an out parameter. The log lines it prints are part of the
//! regression gate, so they reproduce the original wording exactly.

/// Height of the printed histogram, `hist_max`.
const HISTOGRAM_HEIGHT: i64 = 100;

/// `ha_hist_line`: one bar of the ASCII histogram.
fn print_histogram_line(count: i32, bar_length: i64, exceeded: bool, kmer_count: i64) {
    let bar = "*".repeat(bar_length.max(0) as usize);
    let exceed_marker = if exceeded { ">" } else { "" };
    if count >= 0 {
        eprintln!(
            "[M::ha_hist_line] {:5}: {}{} {}",
            count, bar, exceed_marker, kmer_count
        );
    } else {
        eprintln!(
            "[M::ha_hist_line] {:>5}: {}{} {}",
            "rest", bar, exceed_marker, kmer_count
        );
    }
}

/// `adj_m_peak_hom`: reconciles the measured peaks with a user supplied
/// expected homozygous coverage.
pub(crate) fn adjust_homozygous_peak(
    expected_peak: i32,
    main_peak_index: i32,
    left_peak_index: i32,
    right_peak_index: i32,
    heterozygous_peak: &mut i32,
) -> i32 {
    let candidates = [left_peak_index, main_peak_index, right_peak_index];
    let mut chosen_index: i64 = -1;
    let mut smallest_distance: i64 = -1;

    for (index, candidate) in candidates.iter().enumerate() {
        if *candidate <= 0 {
            continue;
        }
        let distance = if *candidate >= expected_peak {
            (*candidate - expected_peak) as i64
        } else {
            (expected_peak - *candidate) as i64
        };
        if smallest_distance == -1
            || smallest_distance > distance
            || (smallest_distance == distance && index == 1)
        {
            smallest_distance = distance;
            chosen_index = index as i64;
        }
    }

    if chosen_index < 0 {
        return expected_peak;
    }

    let chosen = candidates[chosen_index as usize];
    if chosen < expected_peak {
        let distance = (expected_peak - chosen) as i64;
        if distance as f64 >= chosen as f64 * 0.51 {
            *heterozygous_peak = chosen;
            return expected_peak;
        }
    }

    if chosen_index >= 1 {
        for index in (0..chosen_index).rev() {
            if candidates[index as usize] <= 0 {
                continue;
            }
            *heterozygous_peak = candidates[index as usize];
            break;
        }
    }

    chosen
}

/// `ha_analyze_count`: finds the coverage peaks and prints the histogram.
///
/// `counts[occurrence]` holds the number of k-mers observed that many times.
/// Returns the homozygous peak, or `-1` when the coverage is too low.
pub fn analyze_count(
    counts: &[i64],
    start_count: i32,
    expected_homozygous_peak: i64,
    heterozygous_peak: &mut i32,
) -> i32 {
    let count_length = counts.len();
    *heterozygous_peak = -1;

    let start = if counts[1] > 0 { 1 } else { 2 };

    // Find the low point coming from the left.
    let mut low_index = if start > start_count { start } else { start_count };
    let mut index = low_index + 1;
    while index < count_length as i32 {
        if counts[index as usize] > counts[(index - 1) as usize] {
            break;
        }
        index += 1;
    }
    low_index = index - 1;
    eprintln!(
        "[M::ha_analyze_count] lowest: count[{}] = {}",
        low_index, counts[low_index as usize]
    );
    if low_index == count_length as i32 - 1 {
        return -1;
    }

    // Find the highest peak.
    let mut main_peak_index = low_index + 1;
    let mut main_peak_value = counts[main_peak_index as usize];
    for candidate in low_index + 1..count_length as i32 {
        if counts[candidate as usize] > main_peak_value {
            main_peak_value = counts[candidate as usize];
            main_peak_index = candidate;
        }
    }
    eprintln!(
        "[M::ha_analyze_count] highest: count[{}] = {}",
        main_peak_index, counts[main_peak_index as usize]
    );

    // Print the histogram.
    let mut index = start;
    while index < count_length as i32 {
        let scaled = (HISTOGRAM_HEIGHT as f64 * counts[index as usize] as f64
            / counts[main_peak_index as usize] as f64
            + 0.499) as i64;
        let mut bar_length = scaled;
        let mut exceeded = false;
        if bar_length > HISTOGRAM_HEIGHT {
            exceeded = true;
            bar_length = HISTOGRAM_HEIGHT;
        }
        if index > main_peak_index && bar_length == 0 {
            break;
        }
        print_histogram_line(index, bar_length, exceeded, counts[index as usize]);
        index += 1;
    }
    {
        let mut rest = 0i64;
        while index < count_length as i32 {
            rest += counts[index as usize];
            index += 1;
        }
        let mut bar_length = (HISTOGRAM_HEIGHT as f64 * rest as f64
            / counts[main_peak_index as usize] as f64
            + 0.499) as i64;
        let mut exceeded = false;
        if bar_length > HISTOGRAM_HEIGHT {
            exceeded = true;
            bar_length = HISTOGRAM_HEIGHT;
        }
        print_histogram_line(-1, bar_length, exceeded, rest);
    }

    // Look for a smaller peak on the low end.
    let mut left_peak_value = -1i64;
    let mut left_peak_index = -1i32;
    let mut index = main_peak_index - 1;
    while index > low_index {
        if counts[index as usize] >= counts[(index - 1) as usize]
            && counts[index as usize] >= counts[(index + 1) as usize]
            && counts[index as usize] > left_peak_value
        {
            left_peak_value = counts[index as usize];
            left_peak_index = index;
        }
        index -= 1;
    }
    if left_peak_index > low_index && left_peak_index < main_peak_index {
        let mut minimum = main_peak_value;
        let mut probe = left_peak_index + 1;
        while probe < main_peak_index {
            if counts[probe as usize] < minimum {
                minimum = counts[probe as usize];
            }
            probe += 1;
        }
        if (left_peak_value as f64) < main_peak_value as f64 * 0.05
            || minimum as f64 > left_peak_value as f64 * 0.95
        {
            left_peak_value = -1;
            left_peak_index = -1;
        }
    }
    if left_peak_value > 0 {
        eprintln!(
            "[M::ha_analyze_count] left: count[{}] = {}",
            left_peak_index, counts[left_peak_index as usize]
        );
    } else {
        eprintln!("[M::ha_analyze_count] left: none");
    }

    // Look for a smaller peak on the high end.
    let mut right_peak_value = -1i64;
    let mut right_peak_index = -1i32;
    let mut index = main_peak_index + 1;
    while index < count_length as i32 - 1 {
        if counts[index as usize] >= counts[(index - 1) as usize]
            && counts[index as usize] >= counts[(index + 1) as usize]
            && counts[index as usize] > right_peak_value
        {
            right_peak_value = counts[index as usize];
            right_peak_index = index;
        }
        index += 1;
    }
    if right_peak_index > main_peak_index {
        let mut minimum = main_peak_value;
        let mut probe = main_peak_index + 1;
        while probe < right_peak_index {
            if counts[probe as usize] < minimum {
                minimum = counts[probe as usize];
            }
            probe += 1;
        }
        if (right_peak_value as f64) < main_peak_value as f64 * 0.05
            || minimum as f64 > right_peak_value as f64 * 0.95
            || right_peak_index as f64 > main_peak_index as f64 * 2.5
        {
            right_peak_value = -1;
            right_peak_index = -1;
        }
    }
    if right_peak_value > 0 {
        eprintln!(
            "[M::ha_analyze_count] right: count[{}] = {}",
            right_peak_index, counts[right_peak_index as usize]
        );
    } else {
        eprintln!("[M::ha_analyze_count] right: none");
    }

    if expected_homozygous_peak > 0 {
        return adjust_homozygous_peak(
            expected_homozygous_peak as i32,
            main_peak_index,
            left_peak_index,
            right_peak_index,
            heterozygous_peak,
        );
    }

    if right_peak_index > 0 {
        *heterozygous_peak = main_peak_index;
        right_peak_index
    } else {
        if left_peak_index > 0 {
            *heterozygous_peak = left_peak_index;
        }
        main_peak_index
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn synthetic_histogram() -> Vec<i64> {
        // A valley at 5, a homozygous peak at 20 and a heterozygous peak at 10.
        let mut counts = vec![0i64; 64];
        for occurrence in 1..64usize {
            counts[occurrence] = 1;
        }
        counts[10] = 50;
        counts[20] = 100;
        counts
    }

    #[test]
    fn finds_the_highest_peak() {
        let counts = synthetic_histogram();
        let mut heterozygous_peak = -1;
        let homozygous_peak = analyze_count(&counts, 5, -1, &mut heterozygous_peak);
        assert_eq!(homozygous_peak, 20);
        assert_eq!(heterozygous_peak, 10);
    }

    #[test]
    fn reports_low_coverage_when_the_histogram_never_rises() {
        let mut counts = vec![0i64; 64];
        for occurrence in 1..64usize {
            counts[occurrence] = 1;
        }
        let mut heterozygous_peak = -1;
        assert_eq!(analyze_count(&counts, 5, -1, &mut heterozygous_peak), -1);
    }

    #[test]
    fn expected_peak_overrides_the_measurement() {
        let counts = synthetic_histogram();
        let mut heterozygous_peak = -1;
        let homozygous_peak = analyze_count(&counts, 5, 10, &mut heterozygous_peak);
        assert_eq!(homozygous_peak, 10);
    }
}

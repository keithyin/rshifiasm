//! Ports of the two sorting routines from `ksort.h`.
//!
//! `ks_introsort` (quicksort with a comb-sort fallback) is reproduced because
//! several call sites rely on the exact arrangement produced for equal keys,
//! and `radix_sort` because the position index is sorted with it.

/// Slices of at most this size are finished with an insertion sort.
const INSERTION_SORT_THRESHOLD: usize = 16;

/// `RS_MIN_SIZE`: buckets below this size are insertion sorted.
const RADIX_INSERTION_THRESHOLD: usize = 64;

/// `ks_introsort`: sorts `values` according to `less`.
pub fn introsort_by<T, F>(values: &mut [T], less: &F)
where
    F: Fn(&T, &T) -> bool,
{
    let length = values.len();
    if length < 1 {
        return;
    }
    if length == 2 {
        if less(&values[1], &values[0]) {
            values.swap(0, 1);
        }
        return;
    }

    // The depth limit is 2*ceil(log2(n)); once exhausted a comb sort finishes
    // the range, exactly like the original.
    let mut depth = 2u32;
    while (1u64 << depth) < length as u64 {
        depth += 1;
    }
    depth <<= 1;

    let mut stack: Vec<(usize, usize, u32)> = Vec::new();
    let mut start = 0usize;
    let mut end = length - 1;

    loop {
        if start < end {
            depth -= 1;
            if depth == 0 {
                comb_sort(&mut values[start..=end], less);
                end = start;
                continue;
            }

            let middle = start + ((end - start) >> 1) + 1;
            let mut pivot_index = middle;
            if less(&values[middle], &values[start]) {
                if less(&values[middle], &values[end]) {
                    pivot_index = end;
                }
            } else {
                pivot_index = if less(&values[end], &values[start]) {
                    start
                } else {
                    end
                };
            }

            if pivot_index != end {
                values.swap(pivot_index, end);
            }
            // The pivot stays at `end` for the whole partition, so it is
            // addressed by index and never copied.
            let mut low = start;
            let mut high = end;
            loop {
                loop {
                    low += 1;
                    if !less(&values[low], &values[end]) {
                        break;
                    }
                }
                loop {
                    high -= 1;
                    if low > high || !less(&values[end], &values[high]) {
                        break;
                    }
                }
                if high <= low {
                    break;
                }
                values.swap(low, high);
            }
            values.swap(low, end);

            if low - start > end - low {
                if low - start > INSERTION_SORT_THRESHOLD {
                    stack.push((start, low - 1, depth));
                }
                start = if end - low > INSERTION_SORT_THRESHOLD {
                    low + 1
                } else {
                    end
                };
            } else {
                if end - low > INSERTION_SORT_THRESHOLD {
                    stack.push((low + 1, end, depth));
                }
                end = if low - start > INSERTION_SORT_THRESHOLD {
                    low - 1
                } else {
                    start
                };
            }
        } else {
            match stack.pop() {
                None => {
                    insertion_sort(values, less);
                    return;
                }
                Some((saved_start, saved_end, saved_depth)) => {
                    start = saved_start;
                    end = saved_end;
                    depth = saved_depth;
                }
            }
        }
    }
}

fn insertion_sort<T, F>(values: &mut [T], less: &F)
where
    F: Fn(&T, &T) -> bool,
{
    for index in 1..values.len() {
        let mut position = index;
        while position > 0 && less(&values[position], &values[position - 1]) {
            values.swap(position, position - 1);
            position -= 1;
        }
    }
}

/// `ks_combsort`: used when the quicksort depth limit is exhausted.
fn comb_sort<T, F>(values: &mut [T], less: &F)
where
    F: Fn(&T, &T) -> bool,
{
    const SHRINK_FACTOR: f64 = 1.2473309501039786540366528676643;

    let length = values.len();
    let mut gap = length;
    loop {
        if gap > 2 {
            gap = (gap as f64 / SHRINK_FACTOR) as usize;
            if gap == 9 || gap == 10 {
                gap = 11;
            }
        }

        let mut swapped = false;
        if gap < length {
            for index in 0..length - gap {
                if less(&values[index + gap], &values[index]) {
                    values.swap(index, index + gap);
                    swapped = true;
                }
            }
        }

        if !swapped && gap <= 2 {
            break;
        }
    }

    if gap != 1 {
        insertion_sort(values, less);
    }
}

/// `radix_sort_ha64`: MSD radix sort, eight bits per pass, in place.
pub fn radix_sort_u64(values: &mut [u64]) {
    radix_sort_by_key(values, &|value: &u64| *value);
}

/// MSD radix sort over an arbitrary type, keyed by a 64-bit `key`.
pub fn radix_sort_by_key<T, F>(values: &mut [T], key: &F)
where
    T: Copy,
    F: Fn(&T) -> u64,
{
    radix_sort_by_key_bits(values, key, 64);
}

/// MSD radix sort over an arbitrary type, keyed by a 32-bit `key`.
///
/// `radix_sort_k_mer_hit_self` and `radix_sort_k_mer_hit_off` instantiate
/// `KRADIX_SORT_INIT` with a four-byte key, so the first pass starts at bit 24
/// instead of bit 56.
pub fn radix_sort_by_key_32<T, F>(values: &mut [T], key: &F)
where
    T: Copy,
    F: Fn(&T) -> u32,
{
    radix_sort_by_key_bits(values, &|value: &T| u64::from(key(value)), 32);
}

fn radix_sort_by_key_bits<T, F>(values: &mut [T], key: &F, total_bits: u32)
where
    T: Copy,
    F: Fn(&T) -> u64,
{
    if values.len() <= RADIX_INSERTION_THRESHOLD {
        radix_insertion_sort(values, key);
        return;
    }
    radix_sort_recursive(values, key, 8, total_bits - 8);
}

fn radix_insertion_sort<T, F>(values: &mut [T], key: &F)
where
    T: Copy,
    F: Fn(&T) -> u64,
{
    for index in 1..values.len() {
        let current = values[index];
        if key(&current) < key(&values[index - 1]) {
            let mut position = index;
            while position > 0 && key(&current) < key(&values[position - 1]) {
                values[position] = values[position - 1];
                position -= 1;
            }
            values[position] = current;
        }
    }
}

fn radix_sort_recursive<T, F>(values: &mut [T], key: &F, bits_per_pass: u32, shift: u32)
where
    T: Copy,
    F: Fn(&T) -> u64,
{
    let bucket_total = 1usize << bits_per_pass;
    let mask = (bucket_total - 1) as u64;

    let mut bucket_starts = vec![0usize; bucket_total];
    let mut bucket_ends = vec![0usize; bucket_total];

    for value in values.iter() {
        let bucket = ((key(value) >> shift) & mask) as usize;
        bucket_ends[bucket] += 1;
    }
    let mut running = 0usize;
    for bucket in 0..bucket_total {
        bucket_starts[bucket] = running;
        running += bucket_ends[bucket];
        bucket_ends[bucket] = running;
    }

    // Reorder in place with the same cycle-following permutation the C code
    // performs, so ties keep the original relative arrangement.
    let mut cursor = bucket_starts.clone();
    for bucket in 0..bucket_total {
        while cursor[bucket] != bucket_ends[bucket] {
            let position = cursor[bucket];
            let mut carried = values[position];
            let mut target_bucket = ((key(&carried) >> shift) & mask) as usize;
            while target_bucket != bucket {
                let target_position = cursor[target_bucket];
                let displaced = values[target_position];
                values[target_position] = carried;
                cursor[target_bucket] += 1;
                carried = displaced;
                target_bucket = ((key(&carried) >> shift) & mask) as usize;
            }
            values[position] = carried;
            cursor[bucket] += 1;
        }
    }

    if shift == 0 {
        return;
    }

    let next_shift = shift.saturating_sub(bits_per_pass);
    for bucket in 0..bucket_total {
        let start = bucket_starts[bucket];
        let end = bucket_ends[bucket];
        if end - start > RADIX_INSERTION_THRESHOLD {
            radix_sort_recursive(&mut values[start..end], key, bits_per_pass, next_shift);
        } else if end - start > 1 {
            radix_insertion_sort(&mut values[start..end], key);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reference_sort(values: &[u64]) -> Vec<u64> {
        let mut sorted = values.to_vec();
        sorted.sort_unstable();
        sorted
    }

    #[test]
    fn introsort_sorts_small_and_large_arrays() {
        let small = [5u64, 1, 4, 2, 3];
        let mut sorted_small = small;
        introsort_by(&mut sorted_small, &|a: &u64, b: &u64| a < b);
        assert_eq!(sorted_small.to_vec(), reference_sort(&small));

        let mut large: Vec<u64> = (0..10_000u64).map(|index| (index * 7919) % 10_003).collect();
        introsort_by(&mut large, &|a: &u64, b: &u64| a < b);
        assert_eq!(large, reference_sort(&(0..10_000u64).map(|index| (index * 7919) % 10_003).collect::<Vec<_>>()));
    }

    #[test]
    fn radix_sort_matches_a_comparison_sort() {
        let values: Vec<u64> = (0..5_000u64)
            .map(|index| (index.wrapping_mul(2_654_435_761)) & 0xffff_ffff)
            .collect();
        let mut sorted = values.clone();
        radix_sort_u64(&mut sorted);
        assert_eq!(sorted, reference_sort(&values));
    }

    #[test]
    fn radix_sort_handles_short_input() {
        let mut values = vec![3u64, 1, 2];
        radix_sort_u64(&mut values);
        assert_eq!(values, vec![1, 2, 3]);
    }

    #[test]
    fn radix_sort_is_deterministic_for_equal_keys() {
        let mut first = vec![(7u64, 0u64), (7, 1), (3, 2), (7, 3)];
        let mut second = first.clone();
        radix_sort_by_key(&mut first, &|pair: &(u64, u64)| pair.0);
        radix_sort_by_key(&mut second, &|pair: &(u64, u64)| pair.0);
        assert_eq!(first, second);
        assert_eq!(first[0].1, 2);
    }

    /// Golden records from a C++ probe of `ks_introsort` over `(i32, u32,
    /// u32)` tuples ordered by `(second << 32) | third`, the overlap list
    /// layout. Recorded on 2026-11-07 when the recursion-depth off-by-one in
    /// `introsort_klib_by` was fixed; the permutation among equal keys is the
    /// point of the test.
    const KLIB_PRESORT_INPUT: &str = "\
0 0 9181
1 0 9785
2 0 8761
3 0 10831
4 0 10586
5 0 10831
6 0 10831
7 0 7891
8 0 8575
9 0 10831
10 0 10831
11 0 9275
12 0 8942
13 0 9467
14 0 8387
15 0 10831
16 0 10831
17 0 10831
18 0 10831
19 0 10831
20 0 10831
21 0 9603
22 0 10831
23 0 10629
24 0 10831
26 20 10831
27 1032 10831
28 1215 10831
29 1458 10831
30 1584 10831
31 1675 10831
32 2054 10831
33 3039 10831
34 3045 10831
35 3076 10831
36 3655 10831
37 4061 10831
38 4127 10831
39 4253 10831
40 4349 10831
41 4784 10831
42 4942 10831
43 5233 10831
44 5335 10831
45 5459 10831
46 5837 10831
47 6042 10831
48 6515 10831
49 7045 10831
50 7098 10831
51 7459 10831
52 7515 10831
53 7553 10831
54 8453 10831
55 8624 10831
56 8784 10831
57 8954 10831
58 9139 10831
59 9638 10831
60 9681 10831
61 9827 10831";

    const KLIB_SORT_EXPECTED: &str = "\
7 0 7891
14 0 8387
8 0 8575
2 0 8761
12 0 8942
0 0 9181
11 0 9275
13 0 9467
21 0 9603
1 0 9785
4 0 10586
23 0 10629
10 0 10831
9 0 10831
6 0 10831
15 0 10831
16 0 10831
17 0 10831
18 0 10831
19 0 10831
20 0 10831
5 0 10831
22 0 10831
3 0 10831
24 0 10831
26 20 10831
27 1032 10831
28 1215 10831
29 1458 10831
30 1584 10831
31 1675 10831
32 2054 10831
33 3039 10831
34 3045 10831
35 3076 10831
36 3655 10831
37 4061 10831
38 4127 10831
39 4253 10831
40 4349 10831
41 4784 10831
42 4942 10831
43 5233 10831
44 5335 10831
45 5459 10831
46 5837 10831
47 6042 10831
48 6515 10831
49 7045 10831
50 7098 10831
51 7459 10831
52 7515 10831
53 7553 10831
54 8453 10831
55 8624 10831
56 8784 10831
57 8954 10831
58 9139 10831
59 9638 10831
60 9681 10831
61 9827 10831";

    #[test]
    fn introsort_klib_matches_the_c_reference_permutation() {
        let mut values: Vec<(i32, u32, u32)> = KLIB_PRESORT_INPUT
            .lines()
            .map(|line| {
                let fields: Vec<&str> = line.split_whitespace().collect();
                (
                    fields[0].parse().unwrap(),
                    fields[1].parse().unwrap(),
                    fields[2].parse().unwrap(),
                )
            })
            .collect();
        introsort_klib_by(&mut values, &|a: &(i32, u32, u32), b: &(i32, u32, u32)| {
            ((a.1 as u64) << 32 | a.2 as u64) < ((b.1 as u64) << 32 | b.2 as u64)
        });
        let actual: String = values
            .iter()
            .map(|value| format!("{} {} {}", value.0, value.1, value.2))
            .collect::<Vec<_>>()
            .join("\n");
        assert_eq!(actual, KLIB_SORT_EXPECTED, "the klib transcription diverges");
    }
}

/// A literal transcription of klib's `ks_introsort` (the `KSORT_INIT` macro in
/// the reference's `ksort.h`), element for element.
///
/// The generic [`introsort_by`] implements the same algorithm, but the two
/// diverge in the permutation they leave among *equal* elements, and the
/// reference's output order for the overlap lists depends on it. Used for the
/// two `ks_introsort_or_*` call sites of `lchain_qgen_mcopy_fast`.
pub fn introsort_klib_by<T, F>(values: &mut [T], less: &F)
where
    F: Fn(&T, &T) -> bool,
{
    const INSERTION_THRESHOLD: usize = 16;

    let length = values.len();
    if length < 1 {
        return;
    }
    if length == 2 {
        if less(&values[1], &values[0]) {
            values.swap(0, 1);
        }
        return;
    }

    let mut depth = 2u32;
    while (1u64 << depth) < length as u64 {
        depth += 1;
    }
    depth <<= 1;

    let mut stack: Vec<(usize, usize, u32)> = Vec::new();
    let mut left = 0usize;
    let mut right = length - 1;

    loop {
        if left < right {
            // `if (--d == 0)`: the depth is consumed on entry, so the comb
            // fallback fires one level earlier than a post-decrement check.
            depth -= 1;
            if depth == 0 {
                comb_sort_klib(&mut values[left..=right], less);
                right = left;
                continue;
            }

            // Median of three over `left`, `middle`, `right`; the median
            // element is swapped to `right` and becomes the pivot value.
            let middle = left + ((right - left) >> 1) + 1;
            let mut pivot_index = middle;
            if less(&values[middle], &values[left]) {
                if less(&values[middle], &values[right]) {
                    pivot_index = right;
                }
            } else {
                pivot_index = if less(&values[right], &values[left]) {
                    left
                } else {
                    right
                };
            }
            // The C copies the pivot value (`rp = *k`) before swapping it to
            // `right`; T is not `Copy`, so the pivot is addressed by index —
            // it sits at `right` for the whole partition, never moving.
            if pivot_index != right {
                values.swap(pivot_index, right);
            }

            let mut i = left;
            let mut j = right;
            loop {
                i += 1;
                while less(&values[i], &values[right]) {
                    i += 1;
                }
                j -= 1;
                while i <= j && less(&values[right], &values[j]) {
                    j -= 1;
                }
                if j <= i {
                    break;
                }
                values.swap(i, j);
            }
            values.swap(i, right);

            if i - left > right - i {
                if i - left > INSERTION_THRESHOLD {
                    stack.push((left, i - 1, depth));
                }
                left = if right - i > INSERTION_THRESHOLD { i + 1 } else { right };
            } else {
                if right - i > INSERTION_THRESHOLD {
                    stack.push((i + 1, right, depth));
                }
                right = if i - left > INSERTION_THRESHOLD { i - 1 } else { left };
            }
        } else {
            match stack.pop() {
                None => {
                    insertion_sort(values, less);
                    return;
                }
                Some((l, r, d)) => {
                    left = l;
                    right = r;
                    depth = d;
                }
            }
        }
    }
}

/// `ks_combsort`: the comb sort of the klib introsort's depth-limit fallback,
/// including its final insertion pass.
fn comb_sort_klib<T, F>(values: &mut [T], less: &F)
where
    F: Fn(&T, &T) -> bool,
{
    const SHRINK_FACTOR: f64 = 1.2473309501039786540366528676643;

    let length = values.len();
    let mut gap = length;
    loop {
        if gap > 2 {
            gap = (gap as f64 / SHRINK_FACTOR) as usize;
            if gap == 9 || gap == 10 {
                gap = 11;
            }
        }

        let mut swapped = false;
        if gap < length {
            for index in 0..length - gap {
                if less(&values[index + gap], &values[index]) {
                    values.swap(index, index + gap);
                    swapped = true;
                }
            }
        }

        if !swapped && gap <= 2 {
            break;
        }
    }

    if gap != 1 {
        insertion_sort(values, less);
    }
}

#[doc(hidden)]
pub fn comb_sort_for_test<T, F>(values: &mut [T], less: &F)
where
    F: Fn(&T, &T) -> bool,
{
    comb_sort_klib(values, less);
}

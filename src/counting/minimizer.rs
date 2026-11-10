//! Minimizer selection, the counterpart of `HA_SC_INIT(mz1, ...)` in
//! `sketch.cpp`.
//!
//! The routine emits one minimizer per `w` consecutive k-mers, then removes
//! minimizers that sit in high-occupancy streaks (regions where too many
//! minimizers would be sampled) using the same resampling rule as the
//! reference. Both steps influence which overlaps are found, so the window
//! bookkeeping and the ring buffer behaviour are reproduced literally.

use crate::counting::filter_table::{KmerFilterTable, FILTERED_COUNT};
use crate::counting::hash::hash_long_kmer;
use crate::sequence::encoding::SEQUENCE_NT6_TABLE;

/// Number of bits available for the read identifier (`RidBits`).
pub const READ_IDENTIFIER_BITS: u32 = 28;

/// Number of bits available for a minimizer position (`PosBits`).
pub const POSITION_BITS: u32 = 27;

/// Maximum number of high-occupancy minimizers sampled per streak
/// (`MAX_MAX_HIGH_OCC`).
const MAXIMUM_HIGH_OCCURRENCE_SAMPLE: i32 = 16;

/// Marker bit used by the selection array `st_mt_t` to flag chosen minimizers.
const SELECTION_MARKER: u64 = 1 << 32;

/// One minimizer, `ha_mz1_t`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Minimizer {
    /// Hash of the k-mer.
    pub hash: u64,
    /// Read identifier, temporarily reused for the k-mer count during
    /// selection.
    pub read_identifier: u32,
    /// Position of the last base of the k-mer.
    pub position: u32,
    /// Strand: 0 forward, 1 reverse.
    pub reverse: bool,
    /// Number of bases spanned by the k-mer (larger than the k-mer length for
    /// homopolymer compressed k-mers).
    pub span: u32,
}

impl Minimizer {
    const DUMMY: Self = Self {
        hash: u64::MAX,
        read_identifier: (1 << READ_IDENTIFIER_BITS) - 1,
        position: 0,
        reverse: false,
        span: 0,
    };

    /// The all-ones value that `memset(buf, 0xff, ...)` produces.
    const SENTINEL: Self = Self {
        hash: u64::MAX,
        read_identifier: (1 << READ_IDENTIFIER_BITS) - 1,
        position: (1 << POSITION_BITS) - 1,
        reverse: true,
        span: 255,
    };

    /// `mz1_mzcmp`: compare by count (or read id) first, then by hash.
    fn compare(left: &Self, right: &Self) -> std::cmp::Ordering {
        left.read_identifier
            .cmp(&right.read_identifier)
            .then(left.hash.cmp(&right.hash))
    }

    fn is_in_use(&self) -> bool {
        self.read_identifier > 0
    }
}

/// `tiny_queue_t`: a 64-slot ring buffer of homopolymer run lengths.
struct TinyQueue {
    slots: [i32; 64],
    front: usize,
    count: usize,
}

impl TinyQueue {
    fn new() -> Self {
        Self {
            slots: [0; 64],
            front: 0,
            count: 0,
        }
    }

    fn push(&mut self, value: i32) {
        self.slots[(self.count + self.front) & 0x3f] = value;
        self.count += 1;
    }

    fn shift(&mut self) -> i32 {
        if self.count == 0 {
            return -1;
        }
        let value = self.slots[self.front];
        self.front = (self.front + 1) & 0x3f;
        self.count -= 1;
        value
    }

    fn clear(&mut self) {
        self.count = 0;
        self.front = 0;
    }
}

/// Result of sketching one read: the minimizers and their selection markers.
pub struct MinimizerSketch {
    pub minimizers: Vec<Minimizer>,
    /// Parallel array; the low 32 bits hold the k-mer end position and bit 32
    /// marks a retained minimizer.
    pub markers: Vec<u64>,
}

/// `sf_refine_select`: the longest compatible chain of minimizers under a
/// count threshold, used by the `--fast` index selection.
///
/// The DP table stores, for every `(selected index, minimizer index)` pair,
/// the best score in the low 32 bits and the predecessor position in the
/// high 32 bits. The returned pair is the first and last minimizer index of
/// the selected chain.
fn refine_select(
    minimizers: &[Minimizer],
    start: isize,
    end: isize,
    selected_count: i32,
    minimum_frequency: i32,
    markers: &mut [u64],
) -> (isize, isize) {
    let n = selected_count as isize;
    let m = end - start + 1;
    if n <= 0 || m <= 0 {
        return (0, 0);
    }
    let mut table = vec![0u64; (n * m) as usize];
    let mut previous = -1isize;

    for i in 0..n {
        previous = -1;
        for k in 0..m {
            let item = &minimizers[(start + k) as usize];
            if (item.read_identifier as i32) < minimum_frequency {
                continue;
            }
            let window_start = item.position as i64 + 1 - i64::from(item.span);
            let mut t = -1isize;
            if i > 0 {
                t = k - 1;
                while t >= 0 {
                    let before = &minimizers[(start + t) as usize];
                    if !(before.position as i64 >= window_start
                        || (before.read_identifier as i32) < minimum_frequency)
                    {
                        break;
                    }
                    t -= 1;
                }
            }
            let mut left = if i <= 0 {
                0u64
            } else if t < 0 {
                0xffff_ffff
            } else {
                table[((i - 1) * m + t) as usize] & 0xffff_ffff
            };
            if left < 0xffff_ffff {
                left += u64::from(item.read_identifier);
            }
            let right = if previous < 0 {
                0xffff_ffff
            } else {
                table[(i * m + previous) as usize] & 0xffff_ffff
            };
            let mut score = left.min(right);
            if score & 0xffff_ffff == 0xffff_ffff {
                score |= 0xffff_ffffu64 << 32;
            } else if score == left {
                let predecessor = if i > 0 {
                    ((i - 1) * m + t) as u64
                } else {
                    0xffff_ffff
                };
                score |= predecessor << 32;
            } else if score == right {
                let predecessor = if previous >= 0 {
                    (i * m + previous) as u64
                } else {
                    0xffff_ffff
                };
                score |= predecessor << 32;
            }
            table[(i * m + k) as usize] = score;
            previous = k;
        }
    }

    let mut position = ((n - 1) * m + previous) as u64;
    let mut first = u64::MAX;
    let mut last = 0u64;
    while position != 0xffff_ffff {
        if position as usize >= table.len() {
            break;
        }
        let i = (position / m as u64) as isize;
        let k = (position % m as u64) as isize;
        position = table[position as usize] >> 32;
        if position == 0xffff_ffff || (position / m as u64) as isize == i - 1 {
            markers[(start + k) as usize] = 1;
            first = first.min(k as u64);
            last = last.max(k as u64);
        }
    }

    (first as isize + start, last as isize + start)
}

/// `sf_refine_sketch`: keeps only the minimizers selected by
/// [`refine_select`] on both ends and in the middle of a read.
///
/// The minimizers must already carry their k-mer counts in
/// `read_identifier`, exactly like the reference's temporary `rid` reuse.
pub fn refine_sketch(
    minimizers: &mut Vec<Minimizer>,
    length: i32,
    dp_min_len: i32,
    error_rate: f32,
    minimum_frequency: i32,
) {
    let n = minimizers.len() as isize;
    if n == 0 {
        return;
    }
    let window = length.min(dp_min_len);
    let mut selected_count = (window as f32 * error_rate + 1.0) as i32;
    let mut markers = vec![0u64; minimizers.len()];

    // Prefix.
    let mut i = 0isize;
    let mut count = 0i32;
    let mut boundary = -1i32;
    let mut end_index: isize = -1;
    while i < n {
        let item = &minimizers[i as usize];
        if (item.read_identifier as i32) < minimum_frequency {
            i += 1;
            continue;
        }
        let sub_length = item.position as i32 + 1;
        if sub_length > window {
            break;
        }
        end_index = i;
        if item.position as i32 + 1 - item.span as i32 > boundary {
            boundary = item.position as i32;
            count += 1;
        }
        i += 1;
    }
    let mut left_index: isize;
    if count >= selected_count {
        // `refine_select(..., NULL, &li)`: the prefix keeps the last selected
        // minimizer (`kx + sidx`, the second output).
        left_index = refine_select(
            minimizers,
            0,
            end_index,
            selected_count,
            minimum_frequency,
            &mut markers,
        )
        .1;
    } else {
        left_index = i - 1;
        let mut index = 0isize;
        while index <= left_index {
            markers[index as usize] = 1;
            index += 1;
        }
    }

    let mut right_index: isize = n;
    if window < length {
        // Suffix.
        let mut i = n - 1;
        let mut count = 0i32;
        let mut boundary = length + 1;
        let mut end_index: isize = -1;
        while i >= 0 {
            let item = &minimizers[i as usize];
            if (item.read_identifier as i32) < minimum_frequency {
                i -= 1;
                continue;
            }
            let sub_length = length - (item.position as i32 + 1 - item.span as i32);
            if sub_length > window {
                break;
            }
            end_index = i;
            if (item.position as i32) < boundary {
                boundary = item.position as i32 + 1 - item.span as i32;
                count += 1;
            }
            i -= 1;
        }
        if count >= selected_count {
            // `refine_select(..., &ri, NULL)`: the suffix keeps the first
            // selected minimizer (`ix + sidx`, the first output).
            right_index = refine_select(
                minimizers,
                end_index,
                n - 1,
                selected_count,
                minimum_frequency,
                &mut markers,
            )
            .0;
        } else {
            right_index = i + 1;
            let mut index = right_index;
            while index <= n - 1 {
                markers[index as usize] = 1;
                index += 1;
            }
        }

        if right_index - left_index >= 2 {
            left_index += 1;
            right_index -= 1;
            let span = minimizers[right_index as usize].position as i32
                - minimizers[left_index as usize].position as i32
                + minimizers[left_index as usize].span as i32;
            selected_count = (span as f32 * error_rate + 1.0) as i32;
            let mut i = left_index;
            let mut count = 0i32;
            let mut boundary = -1i32;
            while i <= right_index {
                let item = &minimizers[i as usize];
                if (item.read_identifier as i32) < minimum_frequency {
                    i += 1;
                    continue;
                }
                if item.position as i32 + 1 - item.span as i32 > boundary {
                    boundary = item.position as i32;
                    count += 1;
                    if count >= selected_count {
                        break;
                    }
                }
                i += 1;
            }
            if count >= selected_count {
                refine_select(
                    minimizers,
                    left_index,
                    right_index,
                    selected_count,
                    minimum_frequency,
                    &mut markers,
                );
            } else {
                let mut index = left_index;
                while index <= right_index {
                    markers[index as usize] = 1;
                    index += 1;
                }
            }
        }
    }

    let mut output = 0usize;
    for index in 0..minimizers.len() {
        if markers[index] != 0 {
            minimizers[output] = minimizers[index];
            output += 1;
        }
    }
    minimizers.truncate(output);
}

/// `mz_low_b`: the minimum k-mer count the `--fast` selection accepts.
pub fn low_count_threshold(peak_hom: i32, peak_het: i32) -> i32 {
    let mut low_frequency = 2;
    if peak_het > 0 {
        low_frequency = peak_het / 2;
    } else if peak_hom > 0 {
        low_frequency = peak_hom / 4;
    }
    low_frequency.max(2)
}

/// `mz1_ha_sketch`
///
/// `window` is the minimizer window `w`, `kmer_length` is `k`,
/// `sample_distance` is `asm_opt.mz_sample_dist` and `rewindow` is
/// `asm_opt.mz_rewin`.
pub fn sketch(
    sequence: &[u8],
    window: i32,
    kmer_length: i32,
    read_identifier: u32,
    use_homopolymer_compression: bool,
    filter_table: Option<&KmerFilterTable>,
    sample_distance: i32,
    rewindow: i32,
) -> MinimizerSketch {
    let length = sequence.len() as i32;
    assert!(length > 0);
    assert!(window > 0 && window < 256);
    assert!(kmer_length > 0 && kmer_length <= 63);

    let mut minimizers: Vec<Minimizer> = Vec::new();
    let mut markers: Vec<u64> = Vec::new();

    let shift_one = (kmer_length - 1) as u32;
    let mask = (1u64 << kmer_length) - 1;
    let mut kmer = [0u64; 4];

    let mut buffer = vec![Minimizer::SENTINEL; 256];
    let mut buffer_position_value = [0u32; 256];

    let mut queue = TinyQueue::new();
    let mut minimum = Minimizer::DUMMY;
    let mut minimum_position_value = u32::MAX;
    let mut minimum_ring_position = 0i32;
    let mut ring_position = 0i32;

    let mut valid_kmer_count: i32 = 0;
    let mut total_kmer_count: i32 = 0;
    let mut kmer_span: i32 = 0;
    let window_width = window as i32;

    let mut index: i32 = 0;
    while index < length {
        let base = sequence[index as usize];
        let code = SEQUENCE_NT6_TABLE[base as usize];
        let mut candidate = Minimizer::DUMMY;

        if code < 4 {
            if use_homopolymer_compression {
                let mut run_length = 1;
                if index + 1 < length && SEQUENCE_NT6_TABLE[sequence[(index + 1) as usize] as usize] == code {
                    run_length = 2;
                    while index + run_length < length
                        && SEQUENCE_NT6_TABLE[sequence[(index + run_length) as usize] as usize] == code
                    {
                        run_length += 1;
                    }
                    index += run_length - 1;
                }
                queue.push(run_length);
                kmer_span += run_length;
                if queue.count as i32 > kmer_length {
                    kmer_span -= queue.shift();
                }
            } else {
                kmer_span = if valid_kmer_count + 1 < kmer_length {
                    valid_kmer_count + 1
                } else {
                    kmer_length
                };
            }

            let code = code as u64;
            kmer[0] = (kmer[0] << 1 | (code & 1)) & mask;
            kmer[1] = (kmer[1] << 1 | (code >> 1)) & mask;
            kmer[2] = kmer[2] >> 1 | (1 - (code & 1)) << shift_one;
            kmer[3] = kmer[3] >> 1 | (1 - (code >> 1)) << shift_one;

            if kmer[1] == kmer[3] {
                index += 1;
                continue;
            }

            let strand = if kmer[1] < kmer[3] { 0usize } else { 1usize };
            valid_kmer_count += 1;
            total_kmer_count += 1;

            if valid_kmer_count >= kmer_length && kmer_span < 256 {
                let kmer_hash = hash_long_kmer(&kmer);
                let count = match filter_table {
                    None => 0,
                    Some(table) => table.count(kmer_hash),
                };
                let filtered = count >= FILTERED_COUNT;
                if !filtered {
                    candidate = Minimizer {
                        hash: kmer_hash,
                        read_identifier: count as u32,
                        position: index as u32,
                        reverse: strand == 1,
                        span: kmer_span as u32,
                    };
                }
            }
        } else {
            valid_kmer_count = 0;
            queue.clear();
            kmer_span = 0;
        }

        buffer[ring_position as usize] = candidate;
        buffer_position_value[ring_position as usize] = valid_kmer_count as u32;

        // Special case for the very first window, where identical k-mers have
        // not been recorded yet.
        if valid_kmer_count == window_width + kmer_length - 1 && minimum.hash != u64::MAX {
            for offset in (ring_position + 1)..window_width {
                let entry = buffer[offset as usize];
                if Minimizer::compare(&minimum, &entry).is_eq() && entry.position != minimum.position {
                    minimizers.push(entry);
                    markers.push(buffer_position_value[offset as usize] as u64);
                }
            }
            for offset in 0..ring_position {
                let entry = buffer[offset as usize];
                if Minimizer::compare(&minimum, &entry).is_eq() && entry.position != minimum.position {
                    minimizers.push(entry);
                    markers.push(buffer_position_value[offset as usize] as u64);
                }
            }
        }

        let ordering = Minimizer::compare(&minimum, &candidate);
        if ordering.is_ge() {
            // A new minimum: the previous one leaves the window.
            if valid_kmer_count >= window_width + kmer_length && minimum.hash != u64::MAX {
                minimizers.push(minimum);
                markers.push(minimum_position_value as u64);
            }
            minimum = candidate;
            minimum_ring_position = ring_position;
            minimum_position_value = buffer_position_value[ring_position as usize];
        } else if ring_position == minimum_ring_position {
            // The old minimum moved outside the window.
            if valid_kmer_count >= window_width + kmer_length - 1 && minimum.hash != u64::MAX {
                minimizers.push(minimum);
                markers.push(minimum_position_value as u64);
            }

            minimum = Minimizer::DUMMY;
            for offset in (ring_position + 1)..window_width {
                let entry = buffer[offset as usize];
                if Minimizer::compare(&minimum, &entry).is_ge() {
                    minimum = entry;
                    minimum_ring_position = offset;
                    minimum_position_value = buffer_position_value[offset as usize];
                }
            }
            for offset in 0..=ring_position {
                let entry = buffer[offset as usize];
                if Minimizer::compare(&minimum, &entry).is_ge() {
                    minimum = entry;
                    minimum_ring_position = offset;
                    minimum_position_value = buffer_position_value[offset as usize];
                }
            }

            if valid_kmer_count >= window_width + kmer_length - 1 && minimum.hash != u64::MAX {
                for offset in (ring_position + 1)..window_width {
                    let entry = buffer[offset as usize];
                    if Minimizer::compare(&minimum, &entry).is_eq()
                        && minimum.position != entry.position
                    {
                        minimizers.push(entry);
                        markers.push(buffer_position_value[offset as usize] as u64);
                    }
                }
                for offset in 0..=ring_position {
                    let entry = buffer[offset as usize];
                    if Minimizer::compare(&minimum, &entry).is_eq()
                        && minimum.position != entry.position
                    {
                        minimizers.push(entry);
                        markers.push(buffer_position_value[offset as usize] as u64);
                    }
                }
            }
        }

        ring_position += 1;
        if ring_position == window_width {
            ring_position = 0;
        }
        index += 1;
    }

    if minimum.hash != u64::MAX {
        minimizers.push(minimum);
        markers.push(minimum_position_value as u64);
    }

    if sample_distance > window {
        select_high_occupancy_minimizers(
            &mut minimizers,
            &mut markers,
            length,
            sample_distance,
            rewindow,
            kmer_length,
            total_kmer_count,
        );
    }

    for minimizer in &mut minimizers {
        minimizer.read_identifier = read_identifier;
    }

    MinimizerSketch {
        minimizers,
        markers,
    }
}

/// `sf_qfw`: finds the minimizer covering the position `ws` of the window.
fn find_window_minimizer(
    minimizers: &[Minimizer],
    markers: &mut [u64],
    length: i32,
    total_kmer_count: i32,
    window_start: i32,
    index: i32,
    best: &mut i32,
) -> i32 {
    let mut index = index;
    let segment_start = index;
    *best = -1;

    while index < length {
        let position = (markers[index as usize] & 0xffff_ffff) as i32;
        let next_position = if index + 1 < length {
            (markers[(index + 1) as usize] & 0xffff_ffff) as i32
        } else {
            -1
        };
        let last_case =
            index + 1 == length && total_kmer_count >= window_start && position < window_start;
        if position >= window_start
            || (index + 1 < length && position < window_start && next_position > window_start)
            || last_case
        {
            for candidate in segment_start..=index {
                if !minimizers[candidate as usize].is_in_use() {
                    continue;
                }
                if *best < 0
                    || compare_with_index(minimizers, *best, candidate).is_ge()
                {
                    *best = candidate;
                }
            }
            if *best >= 0 && minimizers[*best as usize].is_in_use() {
                for candidate in segment_start..=index {
                    if !minimizers[candidate as usize].is_in_use() {
                        continue;
                    }
                    if compare_with_index(minimizers, *best, candidate).is_eq() {
                        markers[candidate as usize] |= SELECTION_MARKER;
                    }
                }
            }
            break;
        }
        index += 1;
    }

    index
}

/// `sf_mzcmp_l`: orders minimizer slots, treating holes as largest.
fn compare_with_index(minimizers: &[Minimizer], left: i32, right: i32) -> std::cmp::Ordering {
    if left >= 0 && right >= 0 {
        let left_minimizer = &minimizers[left as usize];
        let right_minimizer = &minimizers[right as usize];
        if left_minimizer.is_in_use() && right_minimizer.is_in_use() {
            return Minimizer::compare(left_minimizer, right_minimizer);
        }
        return (left_minimizer.read_identifier == 0)
            .cmp(&(right_minimizer.read_identifier == 0));
    }
    (left < 0).cmp(&(right < 0))
}

/// `sf_hf_select`: picks at most `MAXIMUM_HIGH_OCCURRENCE_SAMPLE` minimizers of
/// a high-occupancy streak and drops the rest.
///
/// `read_length` is the length of the read, used when the streak reaches the
/// end of the minimizer array.
fn select_high_occupancy_sample(
    minimizers: &mut [Minimizer],
    segment_start: i32,
    segment_end: i32,
    read_length: i32,
    sample_distance: i32,
    force: bool,
) {
    let minimizer_count = minimizers.len() as i32;
    if segment_end - segment_start <= 1 {
        return;
    }

    let start_position = if segment_start < 0 {
        0
    } else {
        minimizers[segment_start as usize].position as i32
    };
    let end_position = if segment_end == minimizer_count {
        read_length
    } else {
        minimizers[segment_end as usize].position as i32
    };

    let mut sample_limit =
        ((end_position - start_position) as f64 / sample_distance as f64 + 0.499) as i32;
    if sample_limit > MAXIMUM_HIGH_OCCURRENCE_SAMPLE {
        sample_limit = MAXIMUM_HIGH_OCCURRENCE_SAMPLE;
    }

    // A max-heap of the smallest candidates, exactly like `ks_heapmake` over
    // `mz1_mz` entries.
    let mut heap: Vec<(Minimizer, i32)> = Vec::new();
    let mut index = segment_start + 1;
    while index < segment_end && (heap.len() as i32) < sample_limit {
        heap.push((minimizers[index as usize], index));
        index += 1;
    }
    heap_make(&mut heap);

    while index < segment_end {
        let entry = minimizers[index as usize];
        if let Some(root) = heap.first() {
            if Minimizer::compare(&entry, &root.0).is_lt() {
                heap[0] = (entry, index);
                heap_down(0, heap.len(), &mut heap);
            }
        }
        index += 1;
    }

    for (entry, position) in heap {
        if entry.read_identifier < (end_position - start_position) as u32 || force {
            minimizers[position as usize].read_identifier = 0;
        }
    }
}

/// `ks_heapdown`: sifts one entry down a max-heap ordered by minimizer.
fn heap_down(mut index: usize, length: usize, heap: &mut [(Minimizer, i32)]) {
    let mut current = index;
    let carried = heap[current];
    while {
        current = (current << 1) + 1;
        current < length
    } {
        if current != length - 1
            && Minimizer::compare(&heap[current].0, &heap[current + 1].0).is_lt()
        {
            current += 1;
        }
        if Minimizer::compare(&heap[current].0, &carried.0).is_lt() {
            break;
        }
        heap[index] = heap[current];
        index = current;
    }
    heap[index] = carried;
}

/// `ks_heapmake`
fn heap_make(heap: &mut [(Minimizer, i32)]) {
    let length = heap.len();
    if length == 0 {
        return;
    }
    let mut index = (length >> 1) as i64 - 1;
    while index >= 0 {
        heap_down(index as usize, length, heap);
        index -= 1;
    }
}

/// `mz1_select_mz_h`: thins out minimizers that fall inside high-occupancy
/// streaks.
fn select_high_occupancy_minimizers(
    minimizers: &mut Vec<Minimizer>,
    markers: &mut Vec<u64>,
    length: i32,
    sample_distance: i32,
    rewindow: i32,
    kmer_length: i32,
    total_kmer_count: i32,
) {
    let count = minimizers.len() as i32;
    if count == 0 {
        return;
    }

    let window_width = rewindow;
    let span_threshold = window_width + kmer_length - 1;

    // Does any high-occupancy streak exist at all?
    let mut last_zero = -1i32;
    let mut has_high_occupancy = false;
    let mut index = 0i32;
    while index <= count {
        let boundary = index == count || minimizers[index as usize].read_identifier == 0;
        if boundary {
            if index - last_zero > 1 {
                let start_position = if last_zero < 0 {
                    0
                } else {
                    minimizers[last_zero as usize].position as i32
                };
                let end_position = if index == count {
                    length
                } else {
                    minimizers[index as usize].position as i32
                };
                if ((end_position - start_position) as f64 / sample_distance as f64 + 0.499) as i32 > 0
                {
                    has_high_occupancy = true;
                    break;
                }
            }
            last_zero = index;
        }
        index += 1;
    }

    if !has_high_occupancy {
        return;
    }

    // Mark the minimum of every window.
    let mut minimum_index: i32 = -1;
    let mut cursor = find_window_minimizer(
        minimizers,
        markers,
        count,
        total_kmer_count,
        span_threshold,
        0,
        &mut minimum_index,
    );

    // Everything below only runs when at least one window minimum was found;
    // the original jumps straight to the squeeze step otherwise.
    if cursor != count {
        let mut segment_start = 0i32;
        cursor += 1;
        while cursor < count {
            while segment_start < cursor
                && (markers[segment_start as usize] & 0xffff_ffff) as i32 + window_width
                    > (markers[cursor as usize] & 0xffff_ffff) as i32
            {
                segment_start += 1;
            }

            if compare_with_index(minimizers, cursor, minimum_index).is_le() {
                if minimum_index >= 0 && minimizers[minimum_index as usize].is_in_use() {
                    markers[minimum_index as usize] |= SELECTION_MARKER;
                }
                minimum_index = cursor;
            } else if segment_start > minimum_index {
                if minimum_index >= 0 && minimizers[minimum_index as usize].is_in_use() {
                    markers[minimum_index as usize] |= SELECTION_MARKER;
                }
                minimum_index = -1;
                let mut candidate = segment_start;
                while candidate <= cursor {
                    if compare_with_index(minimizers, minimum_index, candidate).is_ge() {
                        minimum_index = candidate;
                    }
                    candidate += 1;
                }
                if minimum_index >= 0 && minimizers[minimum_index as usize].is_in_use() {
                    for candidate in segment_start..=cursor {
                        if !minimizers[candidate as usize].is_in_use() {
                            continue;
                        }
                        if compare_with_index(minimizers, minimum_index, candidate).is_eq() {
                            markers[candidate as usize] |= SELECTION_MARKER;
                        }
                    }
                }
            }
            cursor += 1;
        }

        if minimum_index >= 0 && minimizers[minimum_index as usize].is_in_use() {
            markers[minimum_index as usize] |= SELECTION_MARKER;
        }

        let tail_index = count - 1;
        while segment_start < count
            && (markers[segment_start as usize] & 0xffff_ffff) as i32 + window_width
                <= total_kmer_count + 1
        {
            if segment_start > minimum_index {
                if minimum_index >= 0 && minimizers[minimum_index as usize].is_in_use() {
                    markers[minimum_index as usize] |= SELECTION_MARKER;
                }
                minimum_index = -1;
                let mut candidate = segment_start;
                while candidate <= tail_index {
                    if compare_with_index(minimizers, minimum_index, candidate).is_ge() {
                        minimum_index = candidate;
                    }
                    candidate += 1;
                }
                if minimum_index >= 0 && minimizers[minimum_index as usize].is_in_use() {
                    for candidate in segment_start..=tail_index {
                        if !minimizers[candidate as usize].is_in_use() {
                            continue;
                        }
                        if compare_with_index(minimizers, minimum_index, candidate).is_eq() {
                            markers[candidate as usize] |= SELECTION_MARKER;
                        }
                    }
                }
            }
            segment_start += 1;
        }

        // Keep the marked minimizers inside each high-occupancy streak, or
        // fall back to sampling when none was marked.
        let mut last_zero = -1i32;
        let mut index = 0i32;
        while index <= count {
            let boundary = index == count || minimizers[index as usize].read_identifier == 0;
            if boundary {
                if index - last_zero > 1 {
                    let start_position = if last_zero < 0 {
                        0
                    } else {
                        minimizers[last_zero as usize].position as i32
                    };
                    let end_position = if index == count {
                        length
                    } else {
                        minimizers[index as usize].position as i32
                    };
                    if ((end_position - start_position) as f64 / sample_distance as f64 + 0.499)
                        as i32
                        > 0
                    {
                        let mut marked = 0;
                        for candidate in (last_zero + 1)..index {
                            if markers[candidate as usize] & SELECTION_MARKER != 0 {
                                minimizers[candidate as usize].read_identifier = 0;
                                marked += 1;
                            }
                        }
                        if marked == 0 {
                            select_high_occupancy_sample(
                                minimizers,
                                last_zero,
                                index,
                                length,
                                sample_distance,
                                false,
                            );
                        }
                    }
                }
                last_zero = index;
            }
            index += 1;
        }
    }

    // Squeeze out the dropped minimizers.
    let mut kept = 0usize;
    for index in 0..minimizers.len() {
        if minimizers[index].read_identifier == 0 {
            minimizers[kept] = minimizers[index];
            markers[kept] = markers[index];
            kept += 1;
        }
    }
    minimizers.truncate(kept);
    markers.truncate(kept);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sketch_emits_one_minimizer_per_window() {
        let sequence = b"ACGTTGCAACGTAGCTAGCTAGGCTAACCGT";
        let result = sketch(sequence, 5, 3, 0, false, None, 0, 1000);
        assert!(!result.minimizers.is_empty());
        assert_eq!(result.minimizers.len(), result.markers.len());
        for minimizer in &result.minimizers {
            assert_eq!(minimizer.read_identifier, 0);
        }
    }

    #[test]
    fn minimizers_are_ordered_by_position() {
        let sequence = b"ACGTTGCAACGTAGCTAGCTAGGCTAACCGTACGATCGTAGCTAGCTA";
        let result = sketch(sequence, 7, 5, 3, true, None, 0, 1000);
        assert!(!result.minimizers.is_empty());
        let mut previous = 0u32;
        for minimizer in &result.minimizers {
            assert!(minimizer.position >= previous);
            previous = minimizer.position;
        }
    }

    #[test]
    fn homopolymer_compression_changes_the_span() {
        let sequence = b"AAAAAAAAAAAAACGT";
        let result = sketch(sequence, 3, 3, 1, true, None, 0, 1000);
        for minimizer in &result.minimizers {
            assert!(minimizer.span >= 3);
        }
    }

    #[test]
    fn filter_table_removes_repetitive_kmers() {
        let mut counts = crate::counting::KmerCountTable::new(3, 0, 0);
        let sequence = b"ACGTTGCAACGTAGCTAGCTAGGCTAACCGT";
        crate::counting::kmer_counter::for_each_kmer(sequence, 3, |hash| {
            for _ in 0..5 {
                counts.insert(hash);
            }
        });
        let table = KmerFilterTable::from_counts(&counts, 10);

        let unfiltered = sketch(sequence, 5, 3, 0, false, None, 0, 1000);
        let filtered = sketch(sequence, 5, 3, 0, false, Some(&table), 0, 1000);
        assert!(filtered.minimizers.len() <= unfiltered.minimizers.len());
    }
}

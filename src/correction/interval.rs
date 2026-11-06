//! Sliding window bookkeeping for correction, the counterpart of
//! `Window_Pool`, `get_interval` and the overlap-selection part of
//! `Correct_dumy` in `Correct.cpp`.
//!
//! Correction walks a read from left to right. For every window it collects
//! the overlaps that reach into it, in two groups: those whose covered part has
//! exactly the window length (kept at the front of the identifier array) and
//! those that are only partially inside the window (kept growing from the
//! back). The two groups are separated because the first can be verified as a
//! batch while the second must be verified one by one.

use crate::overlap::region::OverlapRegion;

/// Return value of [`get_interval`] when the window reached no overlap.
pub const INTERVAL_NO_OVERLAP: i32 = 0;

/// Return value of [`get_interval`] when the window found overlaps.
pub const INTERVAL_HAS_OVERLAP: i32 = 1;

/// Return value of [`get_interval`] when the window starts after every overlap,
/// which also terminates the correction loop.
pub const INTERVAL_PAST_THE_END: i32 = -2;

/// `OVERLAP`: number of bases shared by two inclusive intervals.
///
/// Unlike a plain set intersection the reference does not clamp at zero, so a
/// negative value means the intervals are disjoint.
pub fn overlap_length(a_start: i64, a_end: i64, b_start: i64, b_end: i64) -> i64 {
    a_end.min(b_end) - a_start.max(b_start) + 1
}

/// Yields the windows of a read, `Window_Pool`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WindowPool {
    pub read_length: i64,
    pub window_length: i64,
    /// Number of windows the read is cut into, `ceil(read_length / window_length)`.
    pub window_count: i64,
    pub window_start: i64,
    pub window_end: i64,
    /// Length of the tail that the last window keeps as context.
    pub tail_length: i64,
    /// Set once the last window has been handed out.
    pub terminal: i32,
}

impl WindowPool {
    /// `init_Window_Pool`.
    pub fn new(read_length: i64, window_length: i64, tail_length: i64) -> Self {
        let mut window_end = window_length - 1;
        if window_end >= read_length {
            window_end = read_length - 1;
        }
        Self {
            read_length,
            window_length,
            window_count: (read_length + window_length - 1) / window_length,
            window_start: 0,
            window_end,
            tail_length,
            terminal: 0,
        }
    }

    /// `get_Window`: hands out the current window and advances.
    ///
    /// The window that contains the last base of the read is handed out once;
    /// the following call reports the end of the read.
    pub fn next_window(&mut self) -> Option<(i64, i64)> {
        let window_start = self.window_start;
        let window_end = self.window_end;

        if self.window_end == self.read_length - 1 {
            if self.terminal == 1 || self.read_length == 0 {
                return None;
            }
            self.terminal = 1;
        }

        self.window_start += self.window_length;
        self.window_end += self.window_length;
        if self.window_end >= self.read_length {
            self.window_end = self.read_length - 1;
        }

        Some((window_start, window_end))
    }
}

/// The per-read working set that survives between windows, the overlap
/// selection part of `Correct_dumy`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct IntervalBuffer {
    /// Overlap identifiers: full windows at the front, partial ones at the back.
    /// The length of this vector is the reference `size`.
    pub overlap_identifiers: Vec<u64>,
    /// Number of overlaps that cover the whole window.
    pub length: u64,
    /// Number of overlaps that only partly cover the window.
    pub partial_length: u64,
    /// Index the next scan starts from.
    pub start_index: u64,
    /// The corrected sequence built so far.
    pub corrected_read: Vec<u8>,
    /// Number of corrected bases, `corrected_base`.
    pub corrected_base: u64,
}

impl IntervalBuffer {
    /// A buffer with room for `size` overlap identifiers.
    pub fn with_capacity(size: usize) -> Self {
        Self {
            overlap_identifiers: vec![0; size],
            ..Self::default()
        }
    }

    /// `clear_Correct_dumy` for the fields this module owns.
    pub fn clear(&mut self) {
        self.length = 0;
        self.partial_length = 0;
        self.start_index = 0;
        self.corrected_read.clear();
        self.corrected_base = 0;
    }

    /// The overlaps that cover the whole window, in list order.
    pub fn full_windows(&self) -> &[u64] {
        &self.overlap_identifiers[..self.length as usize]
    }

    /// The overlaps that only partly cover the window, in list order.
    ///
    /// They are pushed backwards from the end of the array, so the first one
    /// added sits at `size - 1` and `verify_window` walks `size - i` for `i`
    /// from one, which restores the order of the overlap list.
    pub fn partial_windows(&self) -> Vec<u64> {
        let size = self.overlap_identifiers.len();
        (0..self.partial_length as usize)
            .map(|offset| self.overlap_identifiers[size - 1 - offset])
            .collect()
    }
}

/// `get_interval`: collects the overlaps that reach into `[window_start, window_end]`.
///
/// `block_length` is the window length; an overlap whose covered part equals a
/// whole window goes into the full group, every other one into the partial
/// group.
pub fn get_interval(
    window_start: i64,
    window_end: i64,
    overlaps: &[OverlapRegion],
    buffer: &mut IntervalBuffer,
    block_length: i64,
) -> i32 {
    if window_start == 0 {
        buffer.start_index = 0;
    }

    // Skip the overlaps that end before this window; the scan index is kept so
    // the next window resumes roughly where this one stopped.
    let mut index = buffer.start_index as usize;
    while index < overlaps.len() {
        if window_end < overlaps[index].x_position_start as i64 {
            buffer.start_index = 0;
            buffer.length = 0;
            buffer.partial_length = 0;
            return INTERVAL_NO_OVERLAP;
        }
        buffer.start_index = index as u64;
        break;
    }

    if index >= overlaps.len() {
        buffer.start_index = overlaps.len() as u64;
        buffer.length = 0;
        buffer.partial_length = 0;
        return INTERVAL_PAST_THE_END;
    }

    buffer.length = 0;
    buffer.partial_length = 0;
    let mut found = false;
    let size = buffer.overlap_identifiers.len() as u64;

    while index < overlaps.len() {
        let overlap = &overlaps[index];
        let covered = overlap_length(
            window_start,
            window_end,
            overlap.x_position_start as i64,
            overlap.x_position_end as i64,
        );

        if covered > 0 {
            if covered == block_length && window_end - window_start + 1 == block_length {
                buffer.overlap_identifiers[buffer.length as usize] = index as u64;
                buffer.length += 1;
            } else {
                buffer.partial_length += 1;
                buffer.overlap_identifiers[(size - buffer.partial_length) as usize] = index as u64;
            }
            if !found {
                found = true;
                buffer.start_index = index as u64;
            }
        }

        if overlap.x_position_start as i64 > window_end {
            break;
        }
        index += 1;
    }

    if buffer.length + buffer.partial_length == 0 {
        INTERVAL_NO_OVERLAP
    } else {
        INTERVAL_HAS_OVERLAP
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region(start: u32, end: u32) -> OverlapRegion {
        OverlapRegion {
            x_position_start: start,
            x_position_end: end,
            ..OverlapRegion::default()
        }
    }

    #[test]
    fn overlap_length_matches_the_reference_macro() {
        assert_eq!(overlap_length(0, 99, 0, 99), 100);
        assert_eq!(overlap_length(0, 99, 50, 149), 50);
        assert_eq!(overlap_length(0, 99, 100, 199), 0);
        assert_eq!(overlap_length(0, 99, 150, 199), -50, "disjoint stays negative");
    }

    #[test]
    fn windows_cover_the_read_including_the_last_one() {
        let mut pool = WindowPool::new(250, 100, 0);
        assert_eq!(pool.window_count, 3);

        assert_eq!(pool.next_window(), Some((0, 99)));
        assert_eq!(pool.next_window(), Some((100, 199)));
        assert_eq!(pool.next_window(), Some((200, 249)), "clamped to the read");
        assert_eq!(pool.next_window(), None, "the last window is handed out once");
    }

    #[test]
    fn a_read_shorter_than_the_window_yields_one_window() {
        let mut pool = WindowPool::new(30, 100, 0);
        assert_eq!(pool.window_count, 1);
        assert_eq!(pool.next_window(), Some((0, 29)));
        assert_eq!(pool.next_window(), None);
    }

    #[test]
    fn an_empty_read_yields_no_window() {
        let mut pool = WindowPool::new(0, 100, 0);
        assert_eq!(pool.next_window(), None);
    }

    #[test]
    fn a_whole_window_overlap_goes_to_the_full_group() {
        let overlaps = vec![region(0, 99), region(50, 149)];
        let mut buffer = IntervalBuffer::with_capacity(8);

        let flag = get_interval(0, 99, &overlaps, &mut buffer, 100);
        assert_eq!(flag, INTERVAL_HAS_OVERLAP);
        assert_eq!(buffer.full_windows(), &[0]);
        assert_eq!(buffer.partial_windows(), vec![1]);
        assert_eq!(buffer.length, 1);
        assert_eq!(buffer.partial_length, 1);
    }

    #[test]
    fn an_overlap_that_starts_after_the_window_ends_the_scan() {
        let overlaps = vec![region(200, 299)];
        let mut buffer = IntervalBuffer::with_capacity(4);

        assert_eq!(
            get_interval(0, 99, &overlaps, &mut buffer, 100),
            INTERVAL_NO_OVERLAP
        );
        assert_eq!(buffer.start_index, 0, "the next window restarts the scan");
        assert_eq!(buffer.length, 0);
    }

    #[test]
    fn a_window_after_every_overlap_terminates_the_loop() {
        let overlaps = vec![region(0, 49)];
        let mut buffer = IntervalBuffer::with_capacity(4);
        buffer.start_index = 1;

        assert_eq!(
            get_interval(100, 199, &overlaps, &mut buffer, 100),
            INTERVAL_PAST_THE_END
        );
        assert_eq!(buffer.start_index, 1);
    }

    #[test]
    fn partial_overlaps_keep_their_list_order() {
        // Three overlaps of different lengths inside one window.
        let overlaps = vec![region(10, 39), region(20, 29), region(0, 99)];
        let mut buffer = IntervalBuffer::with_capacity(8);

        get_interval(0, 99, &overlaps, &mut buffer, 100);
        assert_eq!(buffer.full_windows(), &[2]);
        // Pushed from the back but consumed from the far end, so the order of
        // the overlap list is preserved.
        assert_eq!(buffer.partial_windows(), vec![0, 1]);
    }
}

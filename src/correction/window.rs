//! Window cigars, the counterpart of `window_list`, `window_list_alloc` and
//! their helpers in `Correct.cpp` and `Correct.h`.
//!
//! An overlap is refined window by window. Each `WindowList` records the
//! window boundaries plus the outcome of its alignment, and points into one
//! shared cell buffer: cells of all windows are appended in order and a window
//! addresses its own run through `cigar_start` and `cigar_length`.
//!
//! A cell packs the operation in the top two bits and its length in the low
//! fourteen bits, exactly like the reference.

/// Mask of the length field inside a packed cell.
pub const CELL_LENGTH_MASK: u16 = 0x3fff;

/// One window of an overlap, `window_list`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct WindowList {
    /// Query window boundaries, inclusive.
    pub x_start: i32,
    pub x_end: i32,
    /// Target window boundaries, inclusive.
    pub y_start: i32,
    pub y_end: i32,
    /// Bases outside the window that were still aligned.
    pub extra_begin: i16,
    pub extra_end: i16,
    /// Edit distance measured inside the window.
    pub error: i16,
    /// Edit distance the window was allowed to have.
    pub error_threshold: i16,
    /// Index of the first cell of this window.
    pub cigar_start: u32,
    /// Number of cells of this window.
    pub cigar_length: u32,
}

/// The windows of one overlap plus the shared cell buffer, `window_list_alloc`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct WindowCigar {
    pub windows: Vec<WindowList>,
    /// Packed cells, `operation << 14 | length`.
    pub cells: Vec<u16>,
}

impl WindowCigar {
    /// `clear_window_list_alloc`.
    pub fn clear(&mut self) {
        self.windows.clear();
        self.cells.clear();
    }

    /// `append_window_list`: appends a window with an empty cigar.
    #[allow(clippy::too_many_arguments)]
    pub fn append_window(
        &mut self,
        x_start: i32,
        x_end: i32,
        y_start: i32,
        y_end: i32,
        error: i16,
        extra_begin: i16,
        extra_end: i16,
        error_threshold: i16,
    ) -> usize {
        self.windows.push(WindowList {
            x_start,
            x_end,
            y_start,
            y_end,
            extra_begin,
            extra_end,
            error,
            error_threshold,
            cigar_start: 0,
            cigar_length: 0,
        });
        self.windows.len() - 1
    }

    /// `push_cigar_cell`: appends one cell to the shared buffer.
    pub fn push_cell(&mut self, operation: u8, length: u32) {
        let packed = ((operation as u16) << 14) | (length as u16 & CELL_LENGTH_MASK);
        self.cells.push(packed);
    }

    /// Attaches the cells pushed since `start` to the window at `index`.
    ///
    /// The reference records `cidx` and `clen` by hand after every alignment;
    /// this helper makes the intent explicit and keeps the two fields in sync.
    pub fn assign_cells(&mut self, index: usize, start: usize) {
        let length = self.cells.len() - start;
        self.windows[index].cigar_start = start as u32;
        self.windows[index].cigar_length = length as u32;
    }

    /// `get_cigar_cell`: the cell at `offset` inside the window at `index`.
    pub fn cell(&self, index: usize, offset: u32) -> (u8, u32) {
        let packed = self.cells[(self.windows[index].cigar_start + offset) as usize];
        ((packed >> 14) as u8, (packed & CELL_LENGTH_MASK) as u32)
    }

    /// The cells of the window at `index`, in order.
    pub fn window_cells(&self, index: usize) -> impl Iterator<Item = (u8, u32)> + '_ {
        let window = self.windows[index];
        (0..window.cigar_length).map(move |offset| {
            let packed = self.cells[(window.cigar_start + offset) as usize];
            ((packed >> 14) as u8, (packed & CELL_LENGTH_MASK) as u32)
        })
    }

    /// `scan_cigar`: the error inside the first (or last) `scan_length` query
    /// bases of a window.
    ///
    /// `direction == 0` walks the cells forwards, anything else walks them
    /// backwards. Returns whether the requested length was reached, together
    /// with the error counted up to that point.
    pub fn scan_error(&self, index: usize, scan_length: i64, direction: i32) -> (bool, i64) {
        let window = self.windows[index];
        if window.cigar_length == 1 {
            let (operation, _) = self.cell(index, 0);
            if operation == 0 {
                return (true, 0);
            }
        }

        let count = window.cigar_length;
        let mut query_position = 0i64;
        let mut error = 0i64;

        for step in 0..count {
            let offset = if direction == 0 { step } else { count - 1 - step };
            let (operation, length) = self.cell(index, offset);
            let length = length as i64;
            match operation {
                0 => {
                    query_position += length;
                    if query_position >= scan_length {
                        return (true, error);
                    }
                }
                1 => {
                    for _ in 0..length {
                        query_position += 1;
                        error += 1;
                        if query_position >= scan_length {
                            return (true, error);
                        }
                    }
                }
                // The target has bases the query does not: no query base is
                // consumed, only the error grows.
                2 => {
                    error += length;
                }
                3 => {
                    for _ in 0..length {
                        query_position += 1;
                        error += 1;
                        if query_position >= scan_length {
                            return (true, error);
                        }
                    }
                }
                _ => {}
            }
        }

        (false, error)
    }

    /// `scan_cigar_interval`: the error between two query positions of a window,
    /// both inclusive.
    pub fn scan_error_interval(
        &self,
        index: usize,
        scan_start: i64,
        scan_end: i64,
    ) -> (bool, i64) {
        let window = self.windows[index];
        if window.cigar_length == 1 {
            let (operation, _) = self.cell(index, 0);
            if operation == 0 {
                return (true, 0);
            }
        }

        let count = window.cigar_length;
        let mut query_position = 0i64;
        let mut error = 0i64;

        for offset in 0..count {
            let (operation, length) = self.cell(index, offset);
            let length = length as i64;
            match operation {
                0 => {
                    for _ in 0..length {
                        if query_position == scan_start {
                            error = 0;
                        }
                        query_position += 1;
                        if query_position == scan_end + 1 {
                            return (true, error);
                        }
                    }
                }
                1 => {
                    for _ in 0..length {
                        if query_position == scan_start {
                            error = 0;
                        }
                        query_position += 1;
                        error += 1;
                        if query_position == scan_end + 1 {
                            return (true, error);
                        }
                    }
                }
                2 => {
                    error += length;
                }
                3 => {
                    for _ in 0..length {
                        if query_position == scan_start {
                            error = 0;
                        }
                        query_position += 1;
                        error += 1;
                        if query_position == scan_end + 1 {
                            return (true, error);
                        }
                    }
                }
                _ => {}
            }
        }

        (false, error)
    }

    /// `resize_window_list_alloc`: prepares `count` windows and empties the
    /// shared cell buffer.
    ///
    /// The reference marks every entry with `error = -1`, the "not computed"
    /// value, before the caller fills in the ones it keeps.
    pub fn reset_windows(&mut self, count: usize) {
        self.windows.clear();
        self.windows.resize(
            count,
            WindowList {
                error: -1,
                ..WindowList::default()
            },
        );
        self.cells.clear();
    }
}

/// `get_num_wins`: how many blocks of `block_length` cover `[start, end)`.
///
/// The first block starts at the multiple of `block_length` at or below
/// `start`, so the count is measured from that aligned origin.
pub fn number_of_windows(start: i64, end: i64, block_length: i64) -> i64 {
    let aligned_start = (start / block_length) * block_length;
    let covered = end - aligned_start;
    let mut windows = covered / block_length;
    if covered % block_length > 0 {
        windows += 1;
    }
    windows
}

/// `get_win_id_by_s`: window id containing `window_start`, plus its end.
pub fn window_id_by_start(
    overlap_start: i32,
    overlap_end: i32,
    window_start: i64,
    block_length: i64,
) -> (i64, i64) {
    let aligned_start = (overlap_start as i64 / block_length) * block_length;
    let identifier = (window_start - aligned_start) / block_length;
    let mut window_end = aligned_start + (identifier + 1) * block_length - 1;
    if window_end > overlap_end as i64 {
        window_end = overlap_end as i64;
    }
    (identifier, window_end)
}

/// `get_win_id_by_e`: window id containing `window_end`, plus its start.
pub fn window_id_by_end(
    overlap_start: i32,
    window_end: i64,
    block_length: i64,
) -> (i64, i64) {
    let aligned_start = (overlap_start as i64 / block_length) * block_length;
    let identifier = (window_end - aligned_start) / block_length;
    let mut window_start = aligned_start + identifier * block_length;
    if window_start < overlap_start as i64 {
        window_start = overlap_start as i64;
    }
    (identifier, window_start)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cells_pack_operation_and_length() {
        let mut cigar = WindowCigar::default();
        let index = cigar.append_window(0, 99, 0, 99, 0, 0, 0, 10);

        cigar.push_cell(0, 50);
        cigar.push_cell(1, 2);
        cigar.push_cell(2, 1);
        cigar.assign_cells(index, 0);

        assert_eq!(cigar.cells.len(), 3);
        assert_eq!(cigar.cells[0], 50);
        assert_eq!(cigar.cells[1], (1 << 14) | 2);
        assert_eq!(cigar.windows[index].cigar_length, 3);
        assert_eq!(cigar.cell(index, 1), (1, 2));

        let cells: Vec<(u8, u32)> = cigar.window_cells(index).collect();
        assert_eq!(cells, vec![(0, 50), (1, 2), (2, 1)]);
    }

    #[test]
    fn window_cells_are_independent_of_the_shared_buffer() {
        let mut cigar = WindowCigar::default();
        let first = cigar.append_window(0, 49, 0, 49, 0, 0, 0, 5);
        cigar.push_cell(0, 50);
        cigar.assign_cells(first, 0);

        let second_start = cigar.cells.len();
        let second = cigar.append_window(50, 99, 50, 99, 0, 0, 0, 5);
        cigar.push_cell(0, 40);
        cigar.push_cell(1, 10);
        cigar.assign_cells(second, second_start);

        assert_eq!(cigar.window_cells(first).collect::<Vec<_>>(), vec![(0, 50)]);
        assert_eq!(
            cigar.window_cells(second).collect::<Vec<_>>(),
            vec![(0, 40), (1, 10)]
        );
    }

    #[test]
    fn window_count_is_measured_from_the_aligned_origin() {
        // 0..1000 in blocks of 100: ten windows.
        assert_eq!(number_of_windows(0, 1000, 100), 10);
        // 1000..1500: the origin is 1000, five windows.
        assert_eq!(number_of_windows(1000, 1500, 100), 5);
        // 1050..1100 starts inside block 1000..1099, one window.
        assert_eq!(number_of_windows(1050, 1100, 100), 1);
        // 1050..1200 spans 1050..1199, two windows.
        assert_eq!(number_of_windows(1050, 1200, 100), 2);
        // An unaligned start still counts from the aligned origin.
        assert_eq!(number_of_windows(1234, 1500, 100), 3);
    }

    #[test]
    fn window_boundaries_clamp_to_the_overlap() {
        let (identifier, end) = window_id_by_start(1050, 1200, 1100, 100);
        assert_eq!(identifier, 1);
        assert_eq!(end, 1199);

        let (identifier, end) = window_id_by_start(1050, 1150, 1100, 100);
        assert_eq!(identifier, 1);
        assert_eq!(end, 1150, "clamped to the overlap end");

        let (identifier, start) = window_id_by_end(1050, 1199, 100);
        assert_eq!(identifier, 1);
        assert_eq!(start, 1100);

        let (identifier, start) = window_id_by_end(1050, 1099, 100);
        assert_eq!(identifier, 0);
        assert_eq!(start, 1050, "clamped to the overlap start");
    }
}

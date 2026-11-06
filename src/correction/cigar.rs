//! The correction cigar, the counterpart of `Cigar_record` and its `add_*`
//! helpers in `Correct.cpp`.
//!
//! Every cell packs the operation in the low two bits and its length in the
//! upper thirty bits. Mismatches, insertions and deletions additionally record
//! the bases they introduce in `lost_bases`, which is what turns a cigar back
//! into a corrected sequence.

/// Operation codes, matching the reference numbering.
pub const OPERATION_MATCH: u32 = 0;
pub const OPERATION_MISMATCH: u32 = 1;
pub const OPERATION_INSERTION: u32 = 2;
pub const OPERATION_DELETION: u32 = 3;

/// One `Cigar_record`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CigarRecord {
    /// Packed cells: `length << 2 | operation`.
    pub record: Vec<u32>,
    /// Operation of the last cell, `-1` before anything is added.
    pub current_operation: i32,
    /// Length of the sequence the cigar produces.
    pub new_read_length: u32,
    /// Bases of every mismatch, insertion and deletion cell, in order.
    pub lost_bases: Vec<u8>,
}

impl Default for CigarRecord {
    /// A default record is a fresh one, so that `derive(Default)` on the
    /// structs that contain a cigar does not accidentally start with
    /// `current_operation == 0` (which would mean "already a match cell").
    fn default() -> Self {
        Self::new()
    }
}

impl CigarRecord {
    /// A fresh record, with `current_operation` at its `-1` sentinel.
    pub fn new() -> Self {
        Self {
            record: Vec::new(),
            current_operation: -1,
            new_read_length: 0,
            lost_bases: Vec::new(),
        }
    }

    /// `clear_Cigar_record`.
    pub fn clear(&mut self) {
        self.record.clear();
        self.current_operation = -1;
        self.new_read_length = 0;
        self.lost_bases.clear();
    }

    /// Number of cells.
    pub fn cell_count(&self) -> usize {
        self.record.len()
    }

    /// Length of the cell at `index`.
    pub fn cell_length(&self, index: usize) -> u32 {
        self.record[index] >> 2
    }

    /// Operation of the cell at `index`.
    pub fn cell_operation(&self, index: usize) -> u32 {
        self.record[index] & 3
    }

    /// `add_new_cell_to_cigar_record`.
    fn push_cell(&mut self, length: u32, operation: u32) {
        self.record.push((length << 2) | operation);
    }

    /// `add_existing_cell_to_cigar_record`.
    fn extend_last_cell(&mut self, length: u32, operation: u32) {
        let last = self.record.len() - 1;
        let previous_length = self.record[last] >> 2;
        self.record[last] = ((previous_length + length) << 2) | operation;
    }

    /// Appends `length` bases of `bases` to the cigar.
    ///
    /// `operation` selects the branch of `add_cigar_record`: matches only
    /// accumulate length, the other three also store the bases and, except for
    /// deletions, grow the corrected read.
    pub fn add(&mut self, bases: &[u8], length: u32, operation: u32) {
        let same_operation = self.current_operation >= 0
            && self.current_operation as u32 == operation;

        if operation == OPERATION_MATCH {
            if same_operation {
                self.extend_last_cell(length, operation);
            } else {
                self.push_cell(length, operation);
            }
            self.new_read_length += length;
        } else {
            // Mismatch, insertion and deletion all keep the bases.
            if same_operation {
                self.extend_last_cell(length, operation);
            } else {
                self.push_cell(length, operation);
            }
            self.lost_bases
                .extend_from_slice(&bases[..length as usize]);
            if operation != OPERATION_DELETION {
                self.new_read_length += length;
            }
        }

        self.current_operation = operation as i32;
    }

    /// Rewrites one packed cell, used when the second round trims the first
    /// and last cells of a temporary cigar.
    pub fn set_cell(&mut self, index: usize, operation: u32, length: u32) {
        self.record[index] = (length << 2) | operation;
    }

    /// `memmove` of the reference: keeps only cells `[start, end]` and drops
    /// everything before `start`.
    pub fn retain_cells(&mut self, start: usize, count: usize) {
        self.record.drain(0..start);
        self.record.truncate(count);
    }
}

/// Iterates the cigar cells as `(operation, length)` pairs.
pub struct CigarCellIterator<'a> {
    record: &'a CigarRecord,
    index: usize,
}

impl<'a> Iterator for CigarCellIterator<'a> {
    type Item = (u32, u32);

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.record.record.len() {
            return None;
        }
        let operation = self.record.cell_operation(self.index);
        let length = self.record.cell_length(self.index);
        self.index += 1;
        Some((operation, length))
    }
}

impl CigarRecord {
    /// Iterates over the packed cells in insertion order.
    pub fn cells(&self) -> CigarCellIterator<'_> {
        CigarCellIterator {
            record: self,
            index: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adjacent_matches_merge_into_one_cell() {
        let mut cigar = CigarRecord::new();
        cigar.add(b"AAAA", 4, OPERATION_MATCH);
        cigar.add(b"CC", 2, OPERATION_MATCH);

        assert_eq!(cigar.cell_count(), 1);
        assert_eq!(cigar.cell_length(0), 6);
        assert_eq!(cigar.cell_operation(0), OPERATION_MATCH);
        assert_eq!(cigar.new_read_length, 6);
        assert!(cigar.lost_bases.is_empty());
    }

    #[test]
    fn mismatches_store_their_bases_and_grow_the_read() {
        let mut cigar = CigarRecord::new();
        cigar.add(b"AAAA", 4, OPERATION_MATCH);
        cigar.add(b"G", 1, OPERATION_MISMATCH);
        cigar.add(b"T", 1, OPERATION_MISMATCH);

        assert_eq!(cigar.cell_count(), 2);
        assert_eq!(cigar.cell_length(1), 2);
        assert_eq!(cigar.cell_operation(1), OPERATION_MISMATCH);
        assert_eq!(cigar.lost_bases, b"GT");
        assert_eq!(cigar.new_read_length, 6);
    }

    #[test]
    fn deletions_store_bases_without_growing_the_read() {
        let mut cigar = CigarRecord::new();
        cigar.add(b"AAAA", 4, OPERATION_MATCH);
        cigar.add(b"GG", 2, OPERATION_DELETION);

        assert_eq!(cigar.cell_operation(1), OPERATION_DELETION);
        assert_eq!(cigar.lost_bases, b"GG");
        assert_eq!(cigar.new_read_length, 4, "a deletion removes bases");

        cigar.add(b"C", 1, OPERATION_INSERTION);
        assert_eq!(cigar.new_read_length, 5, "an insertion adds a base");
        assert_eq!(cigar.lost_bases, b"GGC");
    }

    #[test]
    fn clear_resets_every_field() {
        let mut cigar = CigarRecord::new();
        cigar.add(b"ACGT", 4, OPERATION_MATCH);
        cigar.clear();

        assert_eq!(cigar, CigarRecord::new());
        assert_eq!(cigar.current_operation, -1);
        assert_eq!(cigar.cells().count(), 0);
    }

    #[test]
    fn cells_iterate_in_order() {
        let mut cigar = CigarRecord::new();
        cigar.add(b"AAAA", 4, OPERATION_MATCH);
        cigar.add(b"G", 1, OPERATION_MISMATCH);
        cigar.add(b"CC", 2, OPERATION_MATCH);

        let cells: Vec<(u32, u32)> = cigar.cells().collect();
        assert_eq!(
            cells,
            vec![
                (OPERATION_MATCH, 4),
                (OPERATION_MISMATCH, 1),
                (OPERATION_MATCH, 2)
            ]
        );
    }
}

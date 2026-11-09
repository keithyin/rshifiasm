//! The target-sequence source of the shared alignment stack, the counterpart
//! of the reference's `uref`/`rref` pair.
//!
//! Every alignment driver of the correction stack (`hc_aln_exz_adv`,
//! `aln_wlst_adv_exz`, `gen_extend_err_exz`, ...) extracts the target bases
//! through one of two sources: a read of [`ReadStore`]
//! (`recover_UC_Read_sub_region`) or a unitig of a unitig graph
//! (`return_str_seq_exz` → `retrieve_u_seq`). The unitig case never reads a
//! stored unitig string: `retrieve_u_seq` (`Process_Read.cpp` 1727) walks the
//! unitig's read list, copies every overlapping read segment in its stored
//! orientation, and reverse-complements the assembled span as a whole when the
//! request targets the reverse strand.

use crate::assembly::asg::MaUtg;
use crate::reads::read_store::ReadStore;

/// `const ul_idx_t *uref` versus `All_reads *rref`: which sequence store the
/// alignment stack reads its target bases from.
#[derive(Debug, Clone, Copy)]
pub enum TargetSource<'a> {
    /// `rref`: reads of the correction input.
    Reads(&'a ReadStore),
    /// `uref`: unitigs of a unitig graph, assembled from `reads` on demand.
    Unitigs {
        unitigs: &'a [MaUtg],
        reads: &'a ReadStore,
    },
}

impl<'a> TargetSource<'a> {
    /// `Get_READ_LENGTH` / `uref->ug->u.a[id].len`.
    pub fn length(&self, identifier: usize) -> i64 {
        match self {
            TargetSource::Reads(reads) => reads.read_length(identifier) as i64,
            TargetSource::Unitigs { unitigs, .. } => unitigs[identifier].length as i64,
        }
    }

    /// `recover_UC_Read_sub_region` / `return_str_seq_exz` (through
    /// [`unitig_sub_region`]): `length` bases of the target starting at
    /// `start`, reverse-complemented as a whole when `strand == 1`.
    pub fn sub_region(&self, identifier: usize, start: i64, length: i64, strand: u8) -> Vec<u8> {
        match self {
            TargetSource::Reads(reads) => {
                reads.read_sub_region(identifier, start, length, strand)
            }
            TargetSource::Unitigs { unitigs, reads } => {
                unitig_sub_region(unitigs, reads, identifier, start, length, strand)
            }
        }
    }
}

/// `retrieve_u_seq` (`Process_Read.cpp` 1727) with the `i_s` output buffer:
/// assembles `[start, start+length)` of unitig `identifier` from its reads.
///
/// The reference writes into a caller-owned buffer and leaves positions no
/// read covers untouched; here they decode to `A`, matching how the port
/// already renders `recover_UC_Read_sub_region`'s out-of-bounds reads.
pub(crate) fn unitig_sub_region(
    unitigs: &[MaUtg],
    reads: &ReadStore,
    identifier: usize,
    start: i64,
    length: i64,
    strand: u8,
) -> Vec<u8> {
    let unitig = &unitigs[identifier];
    // `if(u->m == 0 || u->n == 0) return;` — a read-less unitig leaves the
    // reference's buffer stale; the polished graphs always carry reads.
    let total_length = unitig.length as i64;
    let requested = if length < 0 { total_length } else { length };
    let mut region = vec![b'A'; requested.max(0) as usize];
    if requested <= 0 {
        return region;
    }
    if unitig.reads.is_empty() {
        return region;
    }

    let mut span_start = start;
    let mut span_end = start + requested;
    if strand == 1 {
        let sep = total_length - span_start;
        let ssp = total_length - span_end;
        span_start = ssp;
        span_end = sep;
    }

    // Walk the read segments, copying the pieces that intersect the span.
    let mut segment_begin = 0i64;
    let mut written = 0usize;
    for &packed in &unitig.reads {
        let read_identifier = (packed >> 33) as usize;
        let orientation = ((packed >> 32) & 1) as u8;
        let segment_length = (packed as u32) as i64;
        if segment_length == 0 {
            continue;
        }
        let segment_end = segment_begin + segment_length;
        if segment_end <= span_start {
            segment_begin = segment_end;
            continue;
        }
        if segment_begin >= span_end {
            break;
        }
        let take_start = segment_begin.max(span_start);
        let take_end = segment_end.min(span_end);
        segment_begin = segment_end;
        let piece_length = (take_end - take_start) as usize;
        if piece_length == 0 {
            continue;
        }
        let piece = reads.read_sub_region(
            read_identifier,
            take_start - (segment_end - segment_length),
            take_end - take_start,
            orientation,
        );
        region[written..written + piece_length]
            .copy_from_slice(&piece[..piece_length]);
        written += piece_length;
    }

    if strand == 1 {
        // The reference complements exactly the `e - s` extracted bases, in
        // place. `written` may fall short only when the span runs past the
        // unitig; the padded tail complements to its own 'A's.
        let complemented = requested as usize;
        for index in 0..complemented / 2 {
            let mirror = complemented - 1 - index;
            let first = complement_base(region[index]);
            let last = complement_base(region[mirror]);
            region[mirror] = first;
            region[index] = last;
        }
        if complemented & 1 == 1 {
            let middle = complemented / 2;
            region[middle] = complement_base(region[middle]);
        }
    }

    region
}

/// `RC_CHAR`: `rc_Table[seq_nt6_table[base]]` — ACGT complement, and every
/// other letter (including `N`) maps to `N`.
fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

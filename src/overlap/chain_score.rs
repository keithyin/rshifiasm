//! Score of a chaining link in the new correction chain, the counterpart of
//! `cal_bw` and `comput_sc_ch_ec` in `Hash_Table.cpp`.
//!
//! The packed `KmerHit::count` carries the supporting-occurrence weight in
//! its upper bits and the minimizer span in its low eight bits, matching
//! `minimizers_qgen0`.

use crate::overlap::region::KmerHit;

/// The reference `INT32_MIN`, used to reject a link entirely.
pub const INVALID_SCORE: i32 = i32::MIN;

/// `cal_bw`: allowable diagonal gap given how much of both reads remains.
pub fn chain_band_width(
    earlier: &KmerHit,
    later: &KmerHit,
    band_width_rate: f64,
    query_length: i64,
    target_length: i64,
) -> i32 {
    let mut query_start = earlier.self_offset as i64;
    let mut query_end = later.self_offset as i64 + 1;
    let target_start = earlier.offset as i64;
    let target_end = later.offset as i64 + 1;
    let query_right = query_length - query_end;
    let target_right = target_length - target_end;

    if query_start <= target_start {
        query_start = 0;
    } else {
        query_start -= target_start;
    }

    if query_right <= target_right {
        query_end = query_length;
    } else {
        query_end += target_right;
    }

    // C converts the double product to `int32_t`, truncating toward zero.
    ((query_end - query_start) as f64 * band_width_rate) as i32
}

/// `comput_sc_ch_ec`: score of extending a chain from `earlier` to `later`.
///
/// Returns [`INVALID_SCORE`] when the two seeds cannot be chained; otherwise
/// the positive span score minus the gap and skip penalties.
pub fn chain_link_score(
    earlier: &KmerHit,
    later: &KmerHit,
    band_width_rate: f64,
    gap_penalty: f64,
    skip_penalty: f64,
    query_length: i64,
    target_length: i64,
) -> i32 {
    chain_link_score_impl(
        earlier, later, band_width_rate, gap_penalty, skip_penalty, query_length, target_length,
        true,
    )
}

/// `comput_sc_ch`: the non-error-correction variant, which always keeps the
/// smaller of the linear and adaptive gap penalties.
pub fn chain_link_score_ch(
    earlier: &KmerHit,
    later: &KmerHit,
    band_width_rate: f64,
    gap_penalty: f64,
    skip_penalty: f64,
    query_length: i64,
    target_length: i64,
) -> i32 {
    chain_link_score_impl(
        earlier, later, band_width_rate, gap_penalty, skip_penalty, query_length, target_length,
        false,
    )
}

fn chain_link_score_impl(
    earlier: &KmerHit,
    later: &KmerHit,
    band_width_rate: f64,
    gap_penalty: f64,
    skip_penalty: f64,
    query_length: i64,
    target_length: i64,
    error_correction: bool,
) -> i32 {
    let query_distance = later.self_offset as i64 - earlier.self_offset as i64;
    if query_distance <= 0 {
        return INVALID_SCORE;
    }
    let target_distance = later.offset as i64 - earlier.offset as i64;
    if target_distance <= 0 {
        return INVALID_SCORE;
    }

    let diagonal_gap = (target_distance - query_distance).abs();
    if diagonal_gap > 16
        && diagonal_gap
            > chain_band_width(earlier, later, band_width_rate, query_length, target_length)
                as i64
    {
        return INVALID_SCORE;
    }

    let shared_length = target_distance.min(query_distance);
    let span = (later.count & 0xff) as i64;
    let occurrence_weight = (later.count >> 8) as i64;
    let mut score = if span < shared_length { span } else { shared_length };
    // `normal_w(x, y)`: integer division when `x >= y`, otherwise one.
    score = if score >= occurrence_weight.max(1) {
        score / occurrence_weight.max(1)
    } else {
        1
    };

    if diagonal_gap != 0 || (shared_length > span && shared_length > 0) {
        let shared_as_double = shared_length as f64;
        let mut linear_penalty = gap_penalty * diagonal_gap as f64;
        let adaptive_penalty =
            score as f64 * ((diagonal_gap as f64 / shared_as_double) / band_width_rate);
        // A short gap keeps the smaller of the two penalties, a long gap the
        // larger; this is inverted relative to the intuitive choice but is
        // exactly what the reference does. The non-ec variant always keeps
        // the smaller one.
        if !error_correction || diagonal_gap < 4 {
            // `lin_pen > a_pen ? a_pen : lin_pen` is the smaller value.
            linear_penalty = linear_penalty.min(adaptive_penalty);
        } else {
            // `lin_pen < a_pen ? a_pen : lin_pen` is the larger value.
            linear_penalty = linear_penalty.max(adaptive_penalty);
        }
        linear_penalty += skip_penalty * shared_as_double;
        score -= (linear_penalty as i32) as i64;
    }

    score as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(self_offset: u32, offset: u32, weight: u32, span: u8) -> KmerHit {
        KmerHit {
            read_identifier: 1,
            reverse_strand: false,
            self_offset,
            offset,
            count: (weight << 8) | span as u32,
        }
    }

    #[test]
    fn backwards_or_equal_positions_are_rejected() {
        let first = hit(100, 100, 2, 51);
        assert_eq!(
            chain_link_score(&first, &hit(100, 200, 2, 51), 0.02, 0.3, 0.0003, 1000, 1000),
            INVALID_SCORE
        );
        assert_eq!(
            chain_link_score(&first, &hit(200, 100, 2, 51), 0.02, 0.3, 0.0003, 1000, 1000),
            INVALID_SCORE
        );
    }

    /// Golden values printed by the real `cal_bw` / `comput_sc_ch_ec`.
    #[test]
    fn matches_the_reference_chain_scores() {
        const GAP: f64 = 0.30024778842926025;
        const SKIP: f64 = 0.00030024780269026496;
        const RATE: f64 = 0.02;

        struct Case {
            name: &'static str,
            query_start: u32,
            target_start: u32,
            query_end: u32,
            target_end: u32,
            weight: u32,
            span: u8,
            query_length: i64,
            band: i32,
            score: i32,
        }
        let cases = [
            Case { name: "collinear", query_start: 100, target_start: 100, query_end: 150, target_end: 150, weight: 2, span: 51, query_length: 1000, band: 20, score: 25 },
            Case { name: "gap2", query_start: 100, target_start: 100, query_end: 150, target_end: 152, weight: 2, span: 51, query_length: 1000, band: 19, score: 25 },
            Case { name: "gap5", query_start: 100, target_start: 100, query_end: 150, target_end: 155, weight: 2, span: 51, query_length: 1000, band: 19, score: -100 },
            Case { name: "gap20", query_start: 100, target_start: 100, query_end: 150, target_end: 170, weight: 2, span: 51, query_length: 1000, band: 19, score: INVALID_SCORE },
            Case { name: "gap100", query_start: 100, target_start: 100, query_end: 150, target_end: 250, weight: 2, span: 51, query_length: 1000, band: 18, score: INVALID_SCORE },
            Case { name: "near_end", query_start: 900, target_start: 900, query_end: 950, target_end: 950, weight: 2, span: 51, query_length: 1000, band: 20, score: 25 },
            Case { name: "heavy_weight", query_start: 100, target_start: 100, query_end: 150, target_end: 152, weight: 9, span: 51, query_length: 1000, band: 19, score: 5 },
            Case { name: "short_span", query_start: 100, target_start: 100, query_end: 150, target_end: 152, weight: 2, span: 10, query_length: 1000, band: 19, score: 5 },
        ];

        for case in cases {
            let earlier = hit(case.query_start, case.target_start, case.weight, case.span);
            let later = hit(case.query_end, case.target_end, case.weight, case.span);
            let band = chain_band_width(&earlier, &later, RATE, case.query_length, case.query_length);
            let score = chain_link_score(
                &earlier, &later, RATE, GAP, SKIP, case.query_length, case.query_length,
            );
            assert_eq!(band, case.band, "band of {}", case.name);
            assert_eq!(score, case.score, "score of {}", case.name);
        }
    }
}

//! Parameters of the new chaining pass, `set_lchain_dp_op` in `anchor.cpp`.
//!
//! The two penalty values are computed with single-precision `expf` and a
//! single-precision `0.0005f` literal. Both details are observable in the
//! chain scores, so the Rust port reproduces them exactly instead of using
//! `f64::exp`.

/// Parameters returned by `set_lchain_dp_op`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LchainParameters {
    pub maximum_skip: i64,
    pub maximum_iterations: i64,
    pub maximum_distance: i64,
    pub quick_check: i64,
    pub gap_penalty: f64,
    pub skip_penalty: f64,
}

/// `set_lchain_dp_op(is_accurate, mz_k, ...)`.
///
/// `accurate` selects the HiFi (`0.01`) or ONT (`0.1`) divergence; the
/// returned gap and skip penalties are the products of single-precision
/// constants with a single-precision exponential.
pub fn set_lchain_parameters(accurate: bool, kmer_length: u32) -> LchainParameters {
    let quick_check = i64::from(accurate);
    let divergence: f64 = if accurate { 0.01 } else { 0.1 };

    // `expf` takes a float argument and returns a float, which is then
    // widened to the double variable `tmp` in the reference.
    let exponent = (-divergence * kmer_length as f64) as f32;
    let decay = exponent.exp() as f64;

    let gap_penalty = 0.5f64 * decay;
    // `0.0005f` is a float literal widened to double, not `0.0005`.
    let skip_penalty = (0.0005f32 as f64) * decay;

    LchainParameters {
        maximum_skip: 25,
        maximum_iterations: 5000,
        maximum_distance: 5000,
        quick_check,
        gap_penalty,
        skip_penalty,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Golden values printed by the real `set_lchain_dp_op`.
    #[test]
    fn matches_the_reference_parameters() {
        let cases = [
            (false, 19u32, 0.074784308671951294f64, 7.4784312224015353e-05f64),
            (false, 27, 0.033602755516767502, 3.3602757112812746e-05),
            (false, 35, 0.015098691917955875, 1.509869263510526e-05),
            (false, 43, 0.0067842784337699413, 6.7842787560058759e-06),
            (false, 51, 0.003048373619094491, 3.0483737638844685e-06),
            (true, 19, 0.41347956657409668, 0.00041347958621332226),
            (true, 27, 0.38168975710868835, 0.00038168977523797901),
            (true, 35, 0.35234403610229492, 0.00035234405283773862),
            (true, 43, 0.32525455951690674, 0.00032525457496566934),
            (true, 51, 0.30024778842926025, 0.00030024780269026496),
        ];

        for (accurate, kmer_length, gap, skip) in cases {
            let parameters = set_lchain_parameters(accurate, kmer_length);
            assert_eq!(parameters.maximum_skip, 25);
            assert_eq!(parameters.maximum_iterations, 5000);
            assert_eq!(parameters.maximum_distance, 5000);
            assert_eq!(parameters.quick_check, i64::from(accurate));
            assert!(
                (parameters.gap_penalty - gap).abs() <= f64::EPSILON * gap.abs().max(1e-300),
                "gap penalty of accurate={accurate} k={kmer_length}: {} vs {gap}",
                parameters.gap_penalty
            );
            assert!(
                (parameters.skip_penalty - skip).abs() <= f64::EPSILON * skip.abs().max(1e-300),
                "skip penalty of accurate={accurate} k={kmer_length}: {} vs {skip}",
                parameters.skip_penalty
            );
        }
    }
}

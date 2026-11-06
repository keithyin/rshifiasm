//! Formatting helpers that reproduce the C `printf` conversions still used by
//! the original code, most notably `%g` which has no direct Rust equivalent.

/// C's `%g` conversion: `significant_digits` significant digits, scientific
/// notation only for very small or very large exponents, and trailing zeros
/// removed.
///
/// This is needed because the help text and several log lines embed values
/// formatted with `%.2g` / `%.3g`; using `{:.2}` instead would print `0.80`
/// where the reference prints `0.8`.
pub fn format_general(value: f64, significant_digits: usize) -> String {
    let precision = significant_digits.max(1);

    if value.is_nan() {
        return "nan".to_string();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.to_string();
    }
    if value == 0.0 {
        return if value.is_sign_negative() { "-0".to_string() } else { "0".to_string() };
    }

    // `%e` with `precision - 1` fractional digits tells us the decimal
    // exponent the way the C library computes it.
    let scientific = format!("{:.*e}", precision - 1, value);
    let exponent: i32 = scientific
        .rsplit('e')
        .next()
        .and_then(|text| text.parse().ok())
        .unwrap_or(0);

    if exponent < -4 || exponent >= precision as i32 {
        return trim_trailing_zeros(&normalise_exponent(&scientific));
    }

    let fraction_digits = (precision as i32 - 1 - exponent).max(0) as usize;
    let fixed = format!("{:.*}", fraction_digits, value);
    trim_trailing_zeros(&fixed)
}

/// Rewrites the exponent produced by Rust into the at-least-two-digits form
/// used by C: `1e-5` becomes `1e-05`.
fn normalise_exponent(scientific: &str) -> String {
    match scientific.find('e') {
        None => scientific.to_string(),
        Some(position) => {
            let mantissa = &scientific[..position];
            let exponent: i32 = scientific[position + 1..].parse().unwrap_or(0);
            format!("{}e{}{:02}", mantissa, if exponent < 0 { '-' } else { '+' }, exponent.abs())
        }
    }
}

/// Removes the decimal point and the trailing zeros of a fixed notation
/// number, as `%g` does.
fn trim_trailing_zeros(text: &str) -> String {
    let (mantissa, exponent) = match text.find('e') {
        Some(position) => (&text[..position], &text[position..]),
        None => (text, ""),
    };

    if !mantissa.contains('.') {
        return format!("{}{}", mantissa, exponent);
    }

    let trimmed = mantissa.trim_end_matches('0').trim_end_matches('.');
    format!("{}{}", trimmed, exponent)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_the_c_printf_g_conversion_for_default_values() {
        assert_eq!(format_general(0.8, 2), "0.8");
        assert_eq!(format_general(0.2, 2), "0.2");
        assert_eq!(format_general(0.6, 2), "0.6");
        assert_eq!(format_general(0.5, 3), "0.5");
        assert_eq!(format_general(0.1, 3), "0.1");
        assert_eq!(format_general(1.0, 3), "1");
    }

    #[test]
    fn switches_to_scientific_notation_outside_the_fixed_range() {
        assert_eq!(format_general(1e-5, 3), "1e-05");
        assert_eq!(format_general(123456.0, 3), "1.23e+05");
        assert_eq!(format_general(-1.5e7, 2), "-1.5e+07");
    }

    #[test]
    fn keeps_significant_digits_without_trailing_zeros() {
        assert_eq!(format_general(0.1234, 3), "0.123");
        assert_eq!(format_general(2.50, 3), "2.5");
        assert_eq!(format_general(0.0, 3), "0");
    }
}

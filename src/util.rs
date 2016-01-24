//! Some programming utilities

/// Returns a float rounded upto a certain number of
/// decimal digits
pub fn round_upto_digits(float: f64, decimal_digits: i32) -> f64 {
    let d = 10_f64.powi(decimal_digits);
    (float * d).round() / d
}

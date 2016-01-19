//! Some programming utilities

/// Returns an interger lesser than or equal to
/// the given number
pub fn int(x: f64) -> i64 {
    if x < 0_f64 {
        return x as i64 - 1;
    }
    x as i64
}

/// Returns a float rounded upto a certain number of
/// decimal digits
pub fn RoundUptoDigits(float: f64, decimal_digits: i32) -> f64 {
    let d = 10_f64.powi(decimal_digits);
    (float * d).round() / d
}

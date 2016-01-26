//! Some programming utilities

/// Returns a float rounded upto a certain number of
/// decimal digits
pub fn round_upto_digits(float: f64, decimal_digits: i32) -> f64 {
    let d = 10_f64.powi(decimal_digits);
    (float * d).round() / d
}

/**
Evaluates polynomials using Horner's method

# Arguments

* `$x`: The independent variable
* `$($a),*`: Sequence of coefficients. The first term in this
sequence of arguments should be the constant term, and
followed by the terms for `$x` in ascending powers of `$x`
**/
macro_rules! Horner_eval {
    ($x:expr, $($a:expr),*) => {
        {
            let mut y = 0_f64;
            let mut u = 1.0;

            $(
                y += u * $a;
                u *= $x;
            )*

            y
        }
    }
}

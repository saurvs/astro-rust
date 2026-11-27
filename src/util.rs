/*
Copyright (c) 2015, 2016 Saurav Sachidanand

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

//! Numerical utilities: rounding and polynomial evaluation (Meeus, chapters 1 and 2).
//!
//! # Example
//!
//! ```
//! use astro::util;
//!
//! // 2 + 3x - x² at x = 4, coefficients in ascending powers of x
//! assert_eq!(util::horner(4.0, &[2.0, 3.0, -1.0]), -2.0);
//! assert_eq!(astro::Horner_eval!(4.0, 2.0, 3.0, -1.0), -2.0);
//!
//! // Round only the final result (Meeus, chapter 2)
//! assert_eq!(util::round_upto_digits(1.4 + 1.4, 0), 3.0);
//! ```

/// Rounds a number to a given number of decimal places.
///
/// Rounds to the nearest multiple of 10^(−`decimal_digits`), with ties
/// rounded away from zero. As Meeus stresses (chapter 2), round only the
/// final result of a calculation, never its inputs or intermediate values.
///
/// The rounding is done on the binary value of `float`, so a number such as
/// 1.005, which is stored as slightly less than 1.005, rounds down.
///
/// # Arguments
///
/// * `float`: The number to round
/// * `decimal_digits`: Number of digits to keep after the decimal point
///
/// # Returns
///
/// The rounded number. Values that already have no digits beyond the
/// requested precision (including infinities, NaN, and any value when
/// `decimal_digits` is larger than an `f64` can resolve) are returned
/// unchanged.
#[inline]
pub fn round_upto_digits(float: f64, decimal_digits: u32) -> f64 {
    // 10^309 overflows to infinity, which is caught below.
    let mut d = 1.0;
    for _ in 0..decimal_digits.min(309) {
        d *= 10.0;
    }

    let scaled = float * d;
    // Every f64 of magnitude 2^52 or more is an integer: nothing to round.
    if !scaled.is_finite() || scaled.abs() >= 4_503_599_627_370_496.0 {
        return float;
    }
    scaled.round() / d
}

/// Evaluates a polynomial with Horner's method (Meeus, chapter 1).
///
/// Computes `c0 + x (c1 + x (c2 + ...))`, which needs no powers of `x` and
/// one multiplication and one addition per coefficient. Use this when the
/// coefficients are only known at run time; for literal coefficients the
/// [`Horner_eval!`](crate::Horner_eval) macro does the same without a slice.
///
/// # Arguments
///
/// * `x`: Value of the variable
/// * `coefficients`: Coefficients in ascending powers of `x`,
///   `[c0, c1, c2, ...]`
///
/// # Returns
///
/// The value of the polynomial, or `0.0` if `coefficients` is empty.
#[inline]
pub fn horner(x: f64, coefficients: &[f64]) -> f64 {
    coefficients.iter().rev().fold(0.0, |acc, &c| acc * x + c)
}

/// Evaluates a polynomial with Horner's method (Meeus, chapter 1).
///
/// `Horner_eval!(x, c0, c1, c2, ...)` expands to
/// `c0 + x * (c1 + x * (c2 + ...))`, where `x` is evaluated only once.
///
/// # Arguments
///
/// * `$x`: Value of the variable (`f32` or `f64`)
/// * `$c`: The constant term
/// * `$($a),*`: Coefficients of `x`, `x²`, ..., in ascending powers of `x`
#[macro_export]
macro_rules! Horner_eval
{
    (@nest $x:ident; $c:expr_2021) => { $c };
    (@nest $x:ident; $c:expr_2021, $($a:expr_2021),+) => {
        $c + $x * $crate::Horner_eval!(@nest $x; $($a),+)
    };
    ($x:expr_2021, $c:expr_2021 $(,)?) => { $c };
    ($x:expr_2021, $c:expr_2021, $($a:expr_2021),+ $(,)?) => {
        {
            let x = $x;
            $crate::Horner_eval!(@nest x; $c, $($a),+)
        }
    };
}

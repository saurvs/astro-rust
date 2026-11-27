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

//! Interpolation, extrema and zeros of tabulated functions (Meeus, chapter 3).
//!
//! Most functions here work on a table of values `y` of a function, given
//! for *equally spaced* values of the argument `x` (for example a daily
//! ephemeris). Positions between tabular entries are measured by the
//! *interpolating factor* `n`: the distance from the central tabular
//! argument, in units of the tabular interval, positive towards later
//! entries. If the central entry is at `x₀` and the interval is `h`, then
//! `x = x₀ + n·h`.
//!
//! * [`three_values`] and [`five_values`] interpolate a value.
//! * [`extremum_three_values`] and [`extremum_five_values`] find a maximum or
//!   minimum, and [`zero_three_values`] and [`zero_five_values`] find where
//!   the function is zero.
//! * [`four_values_halfway`] interpolates exactly half-way between two
//!   entries.
//! * [`lagrange`] and [`lagrange_polynomial`] work with unequally spaced
//!   arguments.
//! * [`unwrap_angles`] removes the jumps of a tabulated angle or time at
//!   360° or 24h, which must be done before interpolating it.
//!
//! Use three values when the second differences of the table are nearly
//! constant, and five values otherwise. Pick the central entry closest to the
//! point of interest, so that `|n| <= 0.5`; the formulas remain usable for
//! `|n| <= 1`, but results far outside that range (including extrema and
//! zeros found there) are unreliable.
//!
//! # Example
//!
//! The distance of Mars to the Earth (AU) on 1992 November 7, 8 and 9 at
//! 0h TD, interpolated to November 8 at 4h21m (Meeus, example 3.a):
//!
//! ```
//! use astro::interpol;
//!
//! let n = 4.35 / 24.0; // 4h21m after the central entry, in days
//! let dist = interpol::three_values(0.884226, 0.877366, 0.870531, n);
//! assert!((dist - 0.876125).abs() < 5e-7);
//!
//! // The least distance Sun-Mars in May 1992 from three values 4 days
//! // apart (Meeus, example 3.b)
//! let min = interpol::extremum_three_values(1.3814294, 1.3812213, 1.3812453).unwrap();
//! assert!((min.y - 1.3812030).abs() < 5e-8);
//! assert!((min.n - 0.39660).abs() < 5e-6); // 1.5864 days after May 16.0
//! ```

use crate::iteration;

/// Tolerance on the interpolating factor used when an extremum or a zero
/// has to be found by iteration. The factor is of order one, so this is
/// far below any tabular precision.
const N_TOLERANCE: f64 = 1e-12;

/// Step limit for those iterations. Newton's method converges in a handful
/// of steps when the formulas are used in their range of validity.
const N_MAX_ITERATIONS: u32 = 100;

/// An extremum (maximum or minimum) of a tabulated function.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Extremum {
    /// Interpolating factor of the extremum: its distance from the central
    /// tabular argument, in units of the tabular interval, positive towards
    /// later entries.
    pub n: f64,
    /// Value of the function at the extremum, in the units of the tabulated
    /// values.
    pub y: f64,
}

/// Interpolates a value of a function from three consecutive tabular values
/// (Meeus, eq. 3.3).
///
/// The function is approximated by the parabola through the three values.
///
/// # Arguments
///
/// * `y1`, `y2`, `y3`: Three consecutive values of the function, for
///   equally spaced arguments; `y2` is the central value
/// * `n`: Interpolating factor, measured from the argument of `y2` in units
///   of the tabular interval, positive towards `y3`. Keep `|n| <= 1`
///
/// # Returns
///
/// The interpolated value, in the units of the `y` values.
#[inline]
pub fn three_values(y1: f64, y2: f64, y3: f64, n: f64) -> f64 {
    let a = y2 - y1;
    let b = y3 - y2;
    let c = b - a;

    y2 + n * (a + b + n * c) / 2.0
}

/// Finds the extremum of a function from three consecutive tabular values
/// (Meeus, eqs. 3.4 and 3.5).
///
/// The extremum is that of the parabola through the three values.
///
/// # Arguments
///
/// * `y1`, `y2`, `y3`: Three consecutive values of the function, for
///   equally spaced arguments; `y2` is the central value
///
/// # Returns
///
/// The [`Extremum`], with its interpolating factor measured from the
/// argument of `y2`, or `None` if the three values lie on a straight line
/// (zero second difference) and there is no extremum. Check that `|n| <= 1`;
/// otherwise the result is an extrapolation and should not be trusted.
pub fn extremum_three_values(y1: f64, y2: f64, y3: f64) -> Option<Extremum> {
    let a = y2 - y1;
    let b = y3 - y2;
    let c = b - a;
    if c == 0.0 {
        return None;
    }
    let sum = a + b;
    let extremum = Extremum {
        n: -sum / (2.0 * c),
        y: y2 - sum * sum / (8.0 * c),
    };
    (extremum.n.is_finite() && extremum.y.is_finite()).then_some(extremum)
}

/// Finds where a function is zero, from three consecutive tabular values
/// (Meeus, eqs. 3.6 and 3.7).
///
/// Meeus finds the zero of the interpolating parabola by iterating eq. 3.6
/// or, better when the curvature is strong, by applying the Newton
/// correction of eq. 3.7, both started at `n = 0`. Since the parabola is
/// a quadratic in `n`, this function solves it directly instead, which
/// gives the same zero (the one nearest the central value) to full
/// precision and without iteration.
///
/// # Arguments
///
/// * `y1`, `y2`, `y3`: Three consecutive values of the function, for
///   equally spaced arguments; `y2` is the central value
///
/// # Returns
///
/// The interpolating factor of the zero nearest the central value, measured
/// from the argument of `y2` in units of the tabular interval, or `None` if
/// the interpolating parabola does not cross zero or a value is not finite.
/// If all three values are zero, returns `Some(0.0)`. Check that `|n| <= 1`; otherwise the result is
/// an extrapolation and should not be trusted.
pub fn zero_three_values(y1: f64, y2: f64, y3: f64) -> Option<f64> {
    if !(y1.is_finite() && y2.is_finite() && y3.is_finite()) {
        return None;
    }
    let a = y2 - y1;
    let b = y3 - y2;
    let c = b - a;

    // Twice the interpolated value is  c n² + (a + b) n + 2 y2.
    // Scale the coefficients to avoid overflow in the discriminant.
    let scale = c.abs().max((a + b).abs()).max((2.0 * y2).abs());
    if !scale.is_finite() {
        return None;
    }
    if y2 == 0.0 {
        return Some(0.0);
    }
    let qa = c / scale;
    let qb = (a + b) / scale;
    let qc = 2.0 * y2 / scale;

    if qa == 0.0 {
        // Straight line
        return if qb == 0.0 { None } else { Some(-qc / qb) };
    }

    let discriminant = qb * qb - 4.0 * qa * qc;
    if discriminant < 0.0 {
        return None;
    }
    // Numerically stable roots q/qa and qc/q; qc != 0 so q != 0.
    let q = -0.5 * (qb + qb.signum() * discriminant.sqrt());
    let n1 = q / qa;
    let n2 = qc / q;
    Some(if n1.abs() < n2.abs() { n1 } else { n2 })
}

/// Interpolates a value of a function from five consecutive tabular values
/// (Meeus, eq. 3.8).
///
/// The function is approximated by the polynomial of degree four through
/// the five values. Use this instead of [`three_values`] when the third
/// differences of the table cannot be neglected.
///
/// # Arguments
///
/// * `y1`, `y2`, `y3`, `y4`, `y5`: Five consecutive values of the function,
///   for equally spaced arguments; `y3` is the central value
/// * `n`: Interpolating factor, measured from the argument of `y3` in units
///   of the tabular interval, positive towards `y4`. Keep `|n| <= 1`
///
/// # Returns
///
/// The interpolated value, in the units of the `y` values.
pub fn five_values(y1: f64, y2: f64, y3: f64, y4: f64, y5: f64, n: f64) -> f64 {
    Quartic::new(y1, y2, y3, y4, y5).value(n)
}

/// The interpolating polynomial of five tabular values, written as
/// `y = y3 + q n + p n² + nn n³ + m n⁴` (Meeus' M, N, P and Q of eq. 3.11).
#[derive(Clone, Copy)]
struct Quartic {
    y3: f64,
    q: f64,
    p: f64,
    nn: f64,
    m: f64,
}

impl Quartic {
    fn new(y1: f64, y2: f64, y3: f64, y4: f64, y5: f64) -> Quartic {
        // Differences as in Meeus' table: A..D, E..G, H, J, K
        let a = y2 - y1;
        let b = y3 - y2;
        let c = y4 - y3;
        let d = y5 - y4;
        let e = b - a;
        let f = c - b;
        let g = d - c;
        let h = f - e;
        let j = g - f;
        let k = j - h;

        let m = k / 24.0;
        let nn = (h + j) / 12.0;
        Quartic {
            y3,
            q: (b + c) / 2.0 - nn,
            p: f / 2.0 - m,
            nn,
            m,
        }
    }

    fn value(&self, n: f64) -> f64 {
        Horner_eval!(n, self.y3, self.q, self.p, self.nn, self.m)
    }

    fn slope(&self, n: f64) -> f64 {
        Horner_eval!(n, self.q, 2.0 * self.p, 3.0 * self.nn, 4.0 * self.m)
    }

    fn curvature(&self, n: f64) -> f64 {
        Horner_eval!(n, 2.0 * self.p, 6.0 * self.nn, 12.0 * self.m)
    }
}

/// Finds an extremum of a function from five consecutive tabular values
/// (Meeus, eq. 3.9).
///
/// The extremum is that of the polynomial of degree four through the five
/// values. Meeus solves eq. 3.9 by simple iteration from `n = 0`; that
/// iteration can oscillate without converging (it does so for his own
/// example of sin x tabulated at 29°…33°), so this function solves the same
/// equation by Newton's method, also started from `n = 0`, which converges
/// to the extremum nearest the central value.
///
/// # Arguments
///
/// * `y1`, `y2`, `y3`, `y4`, `y5`: Five consecutive values of the function,
///   for equally spaced arguments; `y3` is the central value
///
/// # Returns
///
/// The [`Extremum`], with its interpolating factor measured from the
/// argument of `y3`, or `None` if no extremum was found (for example if the
/// values lie on a straight line). Check that `|n| <= 1`; otherwise the
/// result is an extrapolation and should not be trusted.
pub fn extremum_five_values(y1: f64, y2: f64, y3: f64, y4: f64, y5: f64) -> Option<Extremum> {
    let poly = Quartic::new(y1, y2, y3, y4, y5);
    let n = iteration::newton(
        |n| poly.slope(n),
        |n| poly.curvature(n),
        0.0,
        N_TOLERANCE,
        N_MAX_ITERATIONS,
    )
    .ok()?
    .x;
    // A stationary point with zero curvature is not an extremum.
    if poly.curvature(n) == 0.0 {
        return None;
    }
    let y = poly.value(n);
    y.is_finite().then_some(Extremum { n, y })
}

/// Finds where a function is zero, from five consecutive tabular values
/// (Meeus, eqs. 3.10 and 3.11).
///
/// The zero is that of the polynomial of degree four through the five
/// values. It is found with the Newton correction of eq. 3.11 started from
/// `n = 0`, which Meeus recommends over the simple iteration of eq. 3.10
/// because it also converges when the curvature is strong; where both
/// converge they give the same zero.
///
/// # Arguments
///
/// * `y1`, `y2`, `y3`, `y4`, `y5`: Five consecutive values of the function,
///   for equally spaced arguments; `y3` is the central value
///
/// # Returns
///
/// The interpolating factor of the zero, measured from the argument of `y3`
/// in units of the tabular interval, or `None` if the iteration did not
/// converge (for example because the polynomial has no zero near the
/// central value). Check that `|n| <= 1`; otherwise the result is an
/// extrapolation and should not be trusted.
pub fn zero_five_values(y1: f64, y2: f64, y3: f64, y4: f64, y5: f64) -> Option<f64> {
    let poly = Quartic::new(y1, y2, y3, y4, y5);
    iteration::newton(
        |n| poly.value(n),
        |n| poly.slope(n),
        0.0,
        N_TOLERANCE,
        N_MAX_ITERATIONS,
    )
    .ok()
    .map(|solution| solution.x)
}

/// Interpolates the value of a function exactly half-way between two
/// tabular values, from four consecutive tabular values (Meeus, eq. 3.12).
///
/// The result is the value at the midpoint of the arguments of `y2` and
/// `y3` of the cubic through the four values; it is accurate when the
/// fourth differences of the table are negligible.
///
/// # Arguments
///
/// * `y1`, `y2`, `y3`, `y4`: Four consecutive values of the function, for
///   equally spaced arguments
///
/// # Returns
///
/// The interpolated value, in the units of the `y` values.
#[inline]
pub fn four_values_halfway(y1: f64, y2: f64, y3: f64, y4: f64) -> f64 {
    (9.0 * (y2 + y3) - y1 - y4) / 16.0
}

/// Checks that `xs` and `ys` describe a valid set of points for Lagrange
/// interpolation.
fn lagrange_points_valid(xs: &[f64], ys: &[f64]) -> bool {
    if xs.is_empty() || xs.len() != ys.len() || xs.iter().any(|x| !x.is_finite()) {
        return false;
    }
    xs.iter()
        .enumerate()
        .all(|(i, xi)| xs[i + 1..].iter().all(|xj| xj != xi))
}

/// Interpolates a value of a function from points with arbitrary, unequally
/// spaced arguments, using Lagrange's formula (Meeus, eq. 3.13).
///
/// The result is the value at `x` of the unique polynomial of degree
/// `xs.len() - 1` through all the points. The points may be given in any
/// order. Accuracy is best for `x` inside the range of the `xs`;
/// extrapolating far outside it amplifies the rounding errors of the `ys`
/// enormously.
///
/// # Arguments
///
/// * `xs`: Arguments of the points; they must all be different and finite
/// * `ys`: Values of the function at those arguments, same length as `xs`
/// * `x`: Argument at which to interpolate, in the units of `xs`
///
/// # Returns
///
/// The interpolated value, or `None` if there are no points, if `xs` and
/// `ys` have different lengths, or if two arguments are equal or one is not
/// finite.
pub fn lagrange(xs: &[f64], ys: &[f64], x: f64) -> Option<f64> {
    if !lagrange_points_valid(xs, ys) {
        return None;
    }

    let mut sum = 0.0;
    for (i, (&xi, &yi)) in xs.iter().zip(ys).enumerate() {
        let mut factor = 1.0;
        for (j, &xj) in xs.iter().enumerate() {
            if j != i {
                factor *= (x - xj) / (xi - xj);
            }
        }
        sum += factor * yi;
    }
    Some(sum)
}

/// Computes the coefficients of the polynomial through a set of points,
/// using Lagrange's formula (Meeus, eq. 3.13 and example 3.g).
///
/// # Arguments
///
/// * `xs`: Arguments of the points; they must all be different and finite
/// * `ys`: Values of the function at those arguments, same length as `xs`
///
/// # Returns
///
/// The `xs.len()` coefficients `[c0, c1, c2, ...]` of the polynomial
/// `c0 + c1 x + c2 x² + ...`, in ascending powers of `x` (the order used by
/// [`crate::util::horner`] and [`Horner_eval!`](crate::Horner_eval)), or
/// `None` for the same invalid inputs as [`lagrange`].
pub fn lagrange_polynomial(xs: &[f64], ys: &[f64]) -> Option<Vec<f64>> {
    if !lagrange_points_valid(xs, ys) {
        return None;
    }

    let len = xs.len();
    let mut coefficients = vec![0.0; len];
    let mut basis = Vec::with_capacity(len);
    for (i, (&xi, &yi)) in xs.iter().zip(ys).enumerate() {
        // basis = product over j != i of (x - x_j), in ascending powers
        basis.clear();
        basis.push(1.0);
        let mut denominator = 1.0;
        for (j, &xj) in xs.iter().enumerate() {
            if j == i {
                continue;
            }
            basis.push(0.0);
            for k in (1..basis.len()).rev() {
                basis[k] = basis[k - 1] - xj * basis[k];
            }
            basis[0] *= -xj;
            denominator *= xi - xj;
        }
        let weight = yi / denominator;
        for (c, b) in coefficients.iter_mut().zip(&basis) {
            *c += weight * b;
        }
    }
    Some(coefficients)
}

/// Removes the jumps of a tabulated angle or time where it wraps around,
/// so that it can be interpolated (Meeus, chapter 3, "Important remarks").
///
/// A right ascension tabulated as 23h56m, then 0h01m, or a longitude going
/// from 11.48° to 358.25°, must be made continuous (0h01m → 24h01m, or
/// 358.25° → −1.75°) before any interpolation. This adds to each value,
/// after the first, the whole number of periods that brings it within half
/// a period of the value before it.
///
/// The interpolated result can be brought back into range afterwards with
/// `rem_euclid(period)` or the helpers in [`crate::angle`].
///
/// # Arguments
///
/// * `values`: Consecutive tabulated values, modified in place. The first
///   value is left unchanged.
/// * `period`: The period of the values, in their units: `2π` for radians,
///   `360.0` for degrees, `24.0` for hours. If it is not finite and
///   positive, `values` is left unchanged.
pub fn unwrap_angles(values: &mut [f64], period: f64) {
    if !period.is_finite() || period <= 0.0 {
        return;
    }
    for i in 1..values.len() {
        let turns = ((values[i] - values[i - 1]) / period).round();
        values[i] -= turns * period;
    }
}

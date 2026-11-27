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
//! Elliptic orbits
//!
//! Kepler's equation and the quantities of an unperturbed elliptic orbit
//! (Meeus, *Astronomical Algorithms*, 2nd ed., chapters 30, 33 and 39).
//!
//! Angles are in radians, distances in AU and times in Julian Ephemeris
//! Days (TT) unless a function says otherwise. Mean motions are in radians
//! per day.
//!
//! Several solvers for Kepler's equation `E = M + e sin(E)` are provided:
//!
//! * [`ecc_anom`] and [`ecc_anom_checked`]: a bracketed Newton solver that
//!   converges to full `f64` precision for every `0 ≤ e < 1` and every `M`.
//!   Use one of these unless you want a specific method from the book.
//! * [`ecc_anom_simple_iter`]: the plain fixed-point iteration (first method).
//! * [`ecc_anom_newton`]: Newton's correction formula (second method), with
//!   the optional limits on the size of the correction that Meeus describes.
//! * [`ecc_anom_binary_search`]: Sinnott's binary search (third method).
//! * [`ecc_anom_low_ecc`]: a closed-form approximation for small `e`
//!   (fourth method).
//!
//! # Example
//!
//! Solve Kepler's equation for e = 0.1 and M = 5° (Meeus, example 30.a),
//! then get the true anomaly and the radius vector of a body with a
//! semimajor axis of 1 AU.
//!
//! ```
//! use astro::orbit::elliptic;
//!
//! let (e, a) = (0.1, 1.0);
//! let E = elliptic::ecc_anom(5_f64.to_radians(), e, 1e-12);
//! assert!((E.to_degrees() - 5.554589).abs() < 1e-6);
//!
//! let v = elliptic::true_anom(E, e);
//! let r = elliptic::rad_vec_frm_ecc_anom(E, a, e);
//! assert!(v > E && r < a);
//! ```

use crate::consts;
use crate::orbit;
use std::f64::consts::{PI, TAU};

#[inline]
fn valid_ecc(ecc: f64) -> bool {
    (0.0..1.0).contains(&ecc)
}

/// Sign function with `sgn(0) = 0`, as BASIC's `SGN`
#[inline]
fn sgn(x: f64) -> f64 {
    if x > 0.0 {
        1.0
    } else if x < 0.0 {
        -1.0
    } else {
        0.0
    }
}

/// Splits `angl` into `(r, k)` with `angl = r + 2πk`, `r` in `[-π, π]` and
/// `k` a whole number
#[inline]
fn split_revolutions(angl: f64) -> (f64, f64) {
    let k = (angl / TAU).round();
    (angl - k * TAU, k)
}

/// Computes the true anomaly of a body in an elliptic orbit from its
/// eccentric anomaly
///
/// Uses Meeus, eq. 30.1, `tan(v/2) = sqrt((1 + e)/(1 - e)) tan(E/2)`.
///
/// # Returns
///
/// * `true_anom`: True anomaly of the body *| in radians*, in the range
///   `(-π, π]`. It is `NaN` if `ecc > 1`.
///
/// # Arguments
///
/// * `ecc_anom`: Eccentric anomaly of the body *| in radians*
/// * `ecc`     : Eccentricity of the orbit, `0 ≤ ecc < 1`
#[inline]
pub fn true_anom(ecc_anom: f64, ecc: f64) -> f64 {
    2.0 * ((1.0 + ecc).sqrt() * (ecc_anom / 2.0).tan()).atan2((1.0 - ecc).sqrt())
}

/// Computes the eccentric anomaly of a body in an elliptic orbit from its
/// true anomaly
///
/// Uses Meeus, eq. 39.1, `tan(E/2) = sqrt((1 - e)/(1 + e)) tan(v/2)`, the
/// inverse of [`true_anom`].
///
/// # Returns
///
/// * `ecc_anom`: Eccentric anomaly of the body *| in radians*, in the range
///   `(-π, π]`. It is `NaN` if `ecc > 1`.
///
/// # Arguments
///
/// * `true_anom`: True anomaly of the body *| in radians*
/// * `ecc`      : Eccentricity of the orbit, `0 ≤ ecc < 1`
#[inline]
pub fn ecc_anom_frm_true_anom(true_anom: f64, ecc: f64) -> f64 {
    2.0 * ((1.0 - ecc).sqrt() * (true_anom / 2.0).tan()).atan2((1.0 + ecc).sqrt())
}

/// Computes the mean anomaly of a body in an elliptic orbit from its
/// eccentric anomaly
///
/// Uses Kepler's equation in the direct form `M = E - e sin(E)` (Meeus,
/// eq. 39.2).
///
/// # Returns
///
/// * `mean_anom`: Mean anomaly of the body *| in radians*, in the same
///   revolution as `ecc_anom`
///
/// # Arguments
///
/// * `ecc_anom`: Eccentric anomaly of the body *| in radians*
/// * `ecc`     : Eccentricity of the orbit
#[inline]
pub fn mn_anom_frm_ecc_anom(ecc_anom: f64, ecc: f64) -> f64 {
    ecc_anom - ecc * ecc_anom.sin()
}

/// Computes the radius vector of a body in an elliptic orbit from its
/// eccentric anomaly
///
/// Uses Meeus, eq. 30.2, `r = a (1 - e cos(E))`.
///
/// # Returns
///
/// * `rad_vec`: Radius vector of the body *| in AU*
///
/// # Arguments
///
/// * `ecc_anom`: Eccentric anomaly of the body *| in radians*
/// * `a`       : Semimajor axis of the orbit *| in AU*
/// * `ecc`     : Eccentricity of the orbit
#[inline]
pub fn rad_vec_frm_ecc_anom(ecc_anom: f64, a: f64, ecc: f64) -> f64 {
    a * (1.0 - ecc * ecc_anom.cos())
}

/// Computes the radius vector of a body in an elliptic orbit from its
/// true anomaly
///
/// Uses Meeus, eq. 30.3, `r = a (1 - e²) / (1 + e cos(v))`. For an orbit
/// given by its perihelion distance instead of its semimajor axis, use
/// [`orbit::rad_vec_frm_true_anom`], which also works for parabolic and
/// hyperbolic orbits.
///
/// # Returns
///
/// * `rad_vec`: Radius vector of the body *| in AU*
///
/// # Arguments
///
/// * `true_anom`: True anomaly of the body *| in radians*
/// * `a`        : Semimajor axis of the orbit *| in AU*
/// * `ecc`      : Eccentricity of the orbit
#[inline]
pub fn rad_vec_frm_true_anom(true_anom: f64, a: f64, ecc: f64) -> f64 {
    a * (1.0 - ecc * ecc) / (1.0 + ecc * true_anom.cos())
}

/// Computes the eccentric anomaly of a body in an elliptic orbit
///
/// Solves Kepler's equation `E = M + e sin(E)` (Meeus, eq. 30.5) with a
/// Newton iteration that is kept inside an interval known to contain the
/// root, falling back to bisection whenever a Newton step would leave it.
/// The iteration always converges, for any `0 ≤ e < 1` (including `e` very
/// close to 1) and any mean anomaly, positive or negative, of any size.
///
/// The solution is always computed to (nearly) full `f64` precision, so it
/// is at least as accurate as `accuracy`. (Earlier releases used the plain
/// fixed-point iteration, whose stopping rule could return a value much
/// less accurate than `accuracy` when `e` was large, and which could loop
/// forever for `accuracy = 0`.)
///
/// See [`ecc_anom_checked`] for a version that reports invalid input with
/// `None`.
///
/// # Returns
///
/// * `ecc_anom`: Eccentric anomaly of the body *| in radians*. It lies in
///   the same revolution as `mean_anom` (that is, `|E - M| ≤ e`), so that
///   `mean_anom` = 370° gives `E` close to 370°, not 10°. It is `NaN` if
///   `ecc` is not in `[0, 1)` or if an argument is not finite.
///
/// # Arguments
///
/// * `mean_anom`: Mean anomaly of the body *| in radians*
/// * `ecc`      : Eccentricity of the orbit, `0 ≤ ecc < 1`
/// * `accuracy` : Desired accuracy for the eccentric anomaly *| in radians*,
///   e.g. `1e-6`. Kept for compatibility: the result is always more
///   accurate than this.
pub fn ecc_anom(mean_anom: f64, ecc: f64, accuracy: f64) -> f64 {
    let _ = accuracy;
    ecc_anom_checked(mean_anom, ecc).unwrap_or(f64::NAN)
}

/// Computes the eccentric anomaly of a body in an elliptic orbit, checking
/// the input
///
/// Same solver as [`ecc_anom`] (bracketed Newton iteration to full `f64`
/// precision, valid for any `0 ≤ e < 1` and any `M`).
///
/// # Returns
///
/// * `Some(ecc_anom)`: Eccentric anomaly of the body *| in radians*, in the
///   same revolution as `mean_anom` (`|E - M| ≤ e`)
/// * `None` if `ecc` is not in `[0, 1)` or if an argument is not finite
///
/// # Arguments
///
/// * `mean_anom`: Mean anomaly of the body *| in radians*
/// * `ecc`      : Eccentricity of the orbit, `0 ≤ ecc < 1`
///
/// # Example
///
/// A very eccentric orbit close to perihelion (Meeus, table 30.B):
///
/// ```
/// use astro::orbit::elliptic;
///
/// let E = elliptic::ecc_anom_checked(7_f64.to_radians(), 0.999).unwrap();
/// assert!((E.to_degrees() - 52.2702615).abs() < 1e-7);
///
/// assert_eq!(elliptic::ecc_anom_checked(1.0, 1.0), None);
/// ```
pub fn ecc_anom_checked(mean_anom: f64, ecc: f64) -> Option<f64> {
    if !mean_anom.is_finite() || !valid_ecc(ecc) {
        return None;
    }
    if ecc == 0.0 {
        return Some(mean_anom);
    }

    // Kepler's equation is odd and 2π-periodic in the sense that
    // E(M + 2πk) = E(M) + 2πk and E(-M) = -E(M), so solve for M in [0, π].
    let (r, k) = split_revolutions(mean_anom);
    let (sign, m) = if r < 0.0 { (-1.0, -r) } else { (1.0, r) };
    if m == 0.0 {
        return Some(k * TAU);
    }

    // f(E) = E - e sin(E) - m is increasing, and f(m) ≤ 0 ≤ f(min(m + e, π)),
    // so the root lies in [lo, hi].
    let mut lo = m;
    let mut hi = (m + ecc).min(PI);
    // Starting value suggested by Danby, kept inside the bracket.
    let mut E = (m + 0.85 * ecc).clamp(lo, hi);

    for _ in 0..200 {
        let f = E - ecc * E.sin() - m;
        if f == 0.0 {
            break;
        }
        if f < 0.0 {
            lo = E;
        } else {
            hi = E;
        }

        let newton = E - f / (1.0 - ecc * E.cos());
        let next = if newton > lo && newton < hi {
            newton
        } else {
            0.5 * (lo + hi)
        };

        let step = (next - E).abs();
        E = next;
        if step <= 2.0 * f64::EPSILON * E.abs() || hi - lo <= 2.0 * f64::EPSILON * hi {
            break;
        }
    }

    Some(sign * E + k * TAU)
}

/// Computes one step of the fixed-point iteration for Kepler's equation
///
/// Returns `M + e sin(E)`, the right-hand side of Meeus, eq. 30.5, which is
/// a better approximation of the eccentric anomaly than `ecc_anom`.
///
/// # Arguments
///
/// * `mean_anom`: Mean anomaly *| in radians*
/// * `ecc`      : Eccentricity of the orbit
/// * `ecc_anom` : Current approximation of the eccentric anomaly
///   *| in radians*
#[inline]
pub fn kepler_simple_step(mean_anom: f64, ecc: f64, ecc_anom: f64) -> f64 {
    mean_anom + ecc * ecc_anom.sin()
}

/// Solves Kepler's equation by simple fixed-point iteration
///
/// This is the first method of Meeus, chapter 30: starting from `E = M`,
/// `E` is replaced by `M + e sin(E)` (eq. 30.5) until two successive values
/// differ by at most `accuracy`.
///
/// The iteration converges for every `0 ≤ e < 1`, but very slowly when `e`
/// is close to 1 and `M` is close to 0 or π, and the difference between
/// the last two values then underestimates the true error by a factor of up
/// to `1 / (1 - e)`. Prefer [`ecc_anom_checked`] unless you specifically
/// want this method.
///
/// # Returns
///
/// * `Some(ecc_anom)`: Eccentric anomaly *| in radians*, in the same
///   revolution as `mean_anom`
/// * `None` if `ecc` is not in `[0, 1)`, if an argument is not finite, or if
///   the iteration has not met `accuracy` after `max_iter` steps
///
/// # Arguments
///
/// * `mean_anom`: Mean anomaly *| in radians*
/// * `ecc`      : Eccentricity of the orbit, `0 ≤ ecc < 1`
/// * `accuracy` : Stop when two successive values differ by at most this
///   *| in radians*
/// * `max_iter` : Maximum number of iterations
pub fn ecc_anom_simple_iter(mean_anom: f64, ecc: f64, accuracy: f64, max_iter: u32) -> Option<f64> {
    if !mean_anom.is_finite() || !valid_ecc(ecc) {
        return None;
    }

    let mut E = mean_anom;
    for _ in 0..max_iter {
        let next = kepler_simple_step(mean_anom, ecc, E);
        let diff = (next - E).abs();
        E = next;
        if diff <= accuracy {
            return Some(E);
        }
    }

    None
}

/// Computes the Newton correction for Kepler's equation
///
/// Returns `(M + e sin(E) - E) / (1 - e cos(E))`, the fraction of Meeus,
/// eq. 30.7, which is to be added to the current approximation `ecc_anom`.
///
/// # Arguments
///
/// * `mean_anom`: Mean anomaly *| in radians*
/// * `ecc`      : Eccentricity of the orbit
/// * `ecc_anom` : Current approximation of the eccentric anomaly
///   *| in radians*
#[inline]
pub fn kepler_newton_corr(mean_anom: f64, ecc: f64, ecc_anom: f64) -> f64 {
    (mean_anom + ecc * ecc_anom.sin() - ecc_anom) / (1.0 - ecc * ecc_anom.cos())
}

/// How [`ecc_anom_newton`] limits the size of each Newton correction
///
/// For a large eccentricity and a starting value `E = M`, the raw Newton
/// correction can be far larger than 1 radian, making the iteration jump
/// around before it settles. Since `|E - M| = |e sin(E)| < 1`, Meeus
/// (chapter 30) describes two ways of keeping each correction small.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NewtonCorrLimit {
    /// Use the raw correction of eq. 30.7
    Unlimited,
    /// Replace the correction `w` by `asin(sin(w))`, which is always in
    /// `[-π/2, π/2]` (a trick due to K. Leingärtner)
    ArcsinOfSin,
    /// Replace a correction larger than the given limit (in radians) by the
    /// limit, keeping its sign (a trick due to J. M. Steele, who suggests a
    /// limit of 0.5). The limit must be positive.
    Clamp(f64),
}

/// Solves Kepler's equation with Newton's correction formula
///
/// This is the second method of Meeus, chapter 30 (eq. 30.7). Starting from
/// `E = M`, the correction `(M + e sin(E) - E) / (1 - e cos(E))`, optionally
/// limited as set by `limit`, is added to `E` until it is at most `accuracy`
/// in absolute value.
///
/// Without a limit, the iteration converges very quickly for moderate `e`
/// but may wander or fail for `e` close to 1 and small `M`; the limits fix
/// that. Prefer [`ecc_anom_checked`] unless you want this method.
///
/// # Returns
///
/// * `Some(ecc_anom)`: Eccentric anomaly *| in radians*, in the same
///   revolution as `mean_anom`
/// * `None` if `ecc` is not in `[0, 1)`, if an argument is not finite, if
///   a `Clamp` limit is not positive, or if the correction has not dropped
///   to `accuracy` after `max_iter` steps
///
/// # Arguments
///
/// * `mean_anom`: Mean anomaly *| in radians*
/// * `ecc`      : Eccentricity of the orbit, `0 ≤ ecc < 1`
/// * `limit`    : How to limit each correction
/// * `accuracy` : Stop when a correction is at most this *| in radians*
/// * `max_iter` : Maximum number of iterations
pub fn ecc_anom_newton(
    mean_anom: f64,
    ecc: f64,
    limit: NewtonCorrLimit,
    accuracy: f64,
    max_iter: u32,
) -> Option<f64> {
    if !mean_anom.is_finite() || !valid_ecc(ecc) {
        return None;
    }
    if let NewtonCorrLimit::Clamp(c) = limit {
        if c.is_nan() || c <= 0.0 {
            return None;
        }
    }

    let mut E = mean_anom;
    for _ in 0..max_iter {
        let w = kepler_newton_corr(mean_anom, ecc, E);
        let corr = match limit {
            NewtonCorrLimit::Unlimited => w,
            NewtonCorrLimit::ArcsinOfSin => w.sin().asin(),
            NewtonCorrLimit::Clamp(c) => {
                if w.abs() > c {
                    c.copysign(w)
                } else {
                    w
                }
            }
        };
        E += corr;
        if !E.is_finite() {
            return None;
        }
        if corr.abs() <= accuracy {
            return Some(E);
        }
    }

    None
}

/// Solves Kepler's equation by binary search
///
/// This is the third method of Meeus, chapter 30, due to R. W. Sinnott:
/// the mean anomaly is reduced to `[0, π]` using the symmetries of the
/// equation, and `E` is found by halving an interval of `[0, π]` 53 times
/// (one step per bit of an `f64` mantissa). It cannot fail for any
/// `0 ≤ e < 1`, and gives `E` to about `1e-15` radians.
///
/// # Returns
///
/// * `Some(ecc_anom)`: Eccentric anomaly *| in radians*, in the same
///   revolution as `mean_anom`
/// * `None` if `ecc` is not in `[0, 1)` or `mean_anom` is not finite
///
/// # Arguments
///
/// * `mean_anom`: Mean anomaly *| in radians*
/// * `ecc`      : Eccentricity of the orbit, `0 ≤ ecc < 1`
pub fn ecc_anom_binary_search(mean_anom: f64, ecc: f64) -> Option<f64> {
    if !mean_anom.is_finite() || !valid_ecc(ecc) {
        return None;
    }

    let (r, k) = split_revolutions(mean_anom);
    let (sign, m) = if r < 0.0 { (-1.0, -r) } else { (1.0, r) };

    let mut E = PI / 2.0;
    let mut d = PI / 4.0;
    for _ in 0..53 {
        let m1 = E - ecc * E.sin();
        E += d * sgn(m - m1);
        d /= 2.0;
    }

    Some(sign * E + k * TAU)
}

/// Computes an approximate eccentric anomaly for a small eccentricity
///
/// This is the fourth method of Meeus, chapter 30 (eq. 30.8),
/// `tan(E) = sin(M) / (cos(M) - e)`, solved in the correct quadrant. It
/// needs no iteration but is only an approximation: the error is about
/// `0.2″` for the Earth's orbit (`e = 0.0167`), `0.03°` for `e = 0.15` and
/// grows quickly with `e`.
///
/// # Returns
///
/// * `ecc_anom`: Approximate eccentric anomaly *| in radians*, in the same
///   revolution as `mean_anom`
///
/// # Arguments
///
/// * `mean_anom`: Mean anomaly *| in radians*
/// * `ecc`      : Eccentricity of the orbit, small (say below 0.1)
pub fn ecc_anom_low_ecc(mean_anom: f64, ecc: f64) -> f64 {
    let E = mean_anom.sin().atan2(mean_anom.cos() - ecc);

    E + TAU * ((mean_anom - E) / TAU).round()
}

/// Computes the equation of the center of an elliptic orbit
///
/// The equation of the center is the difference `v - M` between the true
/// and the mean anomaly. This series in powers of the eccentricity, up to
/// `e⁵` (Meeus, chapter 33), avoids solving Kepler's equation but is only
/// accurate for small eccentricities: the error is below `0.0003″` for
/// `e = 0.03`, `0.45″` for `e = 0.1` and about `5″` for `e = 0.15`.
///
/// # Returns
///
/// * `eq_of_center`: `v - M` *| in radians*
///
/// # Arguments
///
/// * `mean_anom`: Mean anomaly *| in radians*
/// * `ecc`      : Eccentricity of the orbit (small)
pub fn equation_of_center(mean_anom: f64, ecc: f64) -> f64 {
    let e = ecc;
    let e2 = e * e;
    let e3 = e2 * e;
    let e4 = e3 * e;
    let e5 = e4 * e;
    let M = mean_anom;

    (2.0 * e - e3 / 4.0 + 5.0 * e5 / 96.0) * M.sin()
        + (5.0 * e2 / 4.0 - 11.0 * e4 / 24.0) * (2.0 * M).sin()
        + (13.0 * e3 / 12.0 - 43.0 * e5 / 64.0) * (3.0 * M).sin()
        + 103.0 * e4 / 96.0 * (4.0 * M).sin()
        + 1097.0 * e5 / 960.0 * (5.0 * M).sin()
}

/// Computes the radius vector of a body in an elliptic orbit from its mean
/// anomaly, by a series
///
/// This series in powers of the eccentricity, up to `e⁵` (Meeus,
/// chapter 33), avoids solving Kepler's equation but is only accurate for
/// small eccentricities.
///
/// # Returns
///
/// * `rad_vec`: Radius vector of the body *| in the unit of `a`*
///
/// # Arguments
///
/// * `mean_anom`: Mean anomaly *| in radians*
/// * `a`        : Semimajor axis of the orbit *| in AU*
/// * `ecc`      : Eccentricity of the orbit (small)
pub fn rad_vec_series(mean_anom: f64, a: f64, ecc: f64) -> f64 {
    let e = ecc;
    let e2 = e * e;
    let e3 = e2 * e;
    let e4 = e3 * e;
    let e5 = e4 * e;
    let M = mean_anom;

    a * (1.0 + e2 / 2.0
        - (e - 3.0 * e3 / 8.0 + 5.0 * e5 / 192.0) * M.cos()
        - (e2 / 2.0 - e4 / 3.0) * (2.0 * M).cos()
        - (3.0 * e3 / 8.0 - 45.0 * e5 / 128.0) * (3.0 * M).cos()
        - e4 / 3.0 * (4.0 * M).cos()
        - 125.0 * e5 / 384.0 * (5.0 * M).cos())
}

/// Computes the velocity of a body in an elliptic orbit
///
/// Uses the vis-viva equation in the form given by Meeus, chapter 33,
/// `V = 42.1219 sqrt(1/r - 1/(2a))`, neglecting the mass of the body.
///
/// # Returns
///
/// * `velocity`: Instantaneous velocity of the body
///   *| in kilometers per second*
///
/// # Arguments
///
/// * `r`: Radius vector of the body *| in AU*
/// * `a`: Semimajor axis of orbit *| in AU*
#[inline]
pub fn vel(r: f64, a: f64) -> f64 {
    42.1219 * (1.0 / r - 0.5 / a).sqrt()
}

/// Computes the velocity of a body at perihelion in an elliptic orbit
///
/// Meeus, chapter 33.
///
/// # Returns
///
/// * `velocity`: Velocity of the body at perihelion
///   *| in kilometers per second*
///
/// # Arguments
///
/// * `a`: Semimajor axis of orbit *| in AU*
/// * `e`: Eccentricity of orbit, `0 ≤ e < 1`
#[inline]
pub fn perih_vel(a: f64, e: f64) -> f64 {
    29.7847 * ((1.0 + e) / ((1.0 - e) * a)).sqrt()
}

/// Computes the velocity of a body at aphelion in an elliptic orbit
///
/// Meeus, chapter 33.
///
/// # Returns
///
/// * `velocity`: Velocity of the body at aphelion
///   *| in kilometers per second*
///
/// # Arguments
///
/// * `a`: Semimajor axis of orbit *| in AU*
/// * `e`: Eccentricity of orbit, `0 ≤ e < 1`
#[inline]
pub fn aph_vel(a: f64, e: f64) -> f64 {
    29.7847 * ((1.0 - e) / ((1.0 + e) * a)).sqrt()
}

/// Computes the approximate length of an ellipse using the Ramanujan
/// method
///
/// # Returns
///
/// * `approx_length`: An approximate value for the length of
///   the ellipse (same unit as that of `a`
///   and `b`), using a formula given by
///   [Ramanujan](https://en.wikipedia.org/wiki/Srinivasa_Ramanujan)
///   in 1914 (Meeus, chapter 33).
///
/// The result is too small by:
///
/// * 0% for a = b
/// * 0.4155% for e = 1
///
/// # Arguments
///
/// * `a`: Semimajor axis of the ellipse (same unit as that of `b`)
/// * `b`: Semiminor axis of the ellipse (same unit as that of `a`)
#[inline]
pub fn length_ramanujan(a: f64, b: f64) -> f64 {
    PI * (3.0 * (a + b) - ((a + 3.0 * b) * (3.0 * a + b)).sqrt())
}

/// Computes the approximate length of an ellipse
///
/// Uses a combination of the arithmetic, geometric and harmonic means of
/// the semi-axes (Meeus, chapter 33).
///
/// # Returns
///
/// * `approx_length`: An approximate value for the length of the ellipse
///   (same unit as that of `a` and `b`)
///
/// The error in `approx_length` is:
///
/// * less than 0.001% if e < 0.88
/// * less than 0.01% if e < 0.95
/// * 1% if e = 0.9997
/// * 3% if e = 1
///
/// # Arguments
///
/// * `a`: Semimajor axis of the ellipse (same unit as that of `b`)
/// * `b`: Semiminor axis of the ellipse (same unit as that of `a`)
pub fn length(a: f64, b: f64) -> f64 {
    let A = (a + b) / 2.0;
    let G = (a * b).sqrt();
    let H = (2.0 * a * b) / (a + b);

    PI * (21.0 * A - 2.0 * G - 3.0 * H) / 8.0
}

/// Computes the length of an ellipse accurately, by a series
///
/// Sums the series in powers of `m = (a - b)/(a + b)` given by Meeus,
/// chapter 33 (method 4, also known as the Gauss–Kummer series),
/// `L = π (a + b) [1 + (1/2)² m² + (1/(2·4))² m⁴ + (1·3/(2·4·6))² m⁶ + ...]`,
/// until the terms no longer change the sum.
///
/// The result is accurate to about `1e-15` (relative) for eccentricities
/// up to 0.99, and to better than `1e-11` for any eccentricity (the series
/// converges slowly for a very flat ellipse). A degenerate ellipse with
/// `b = 0` has length `4a`.
///
/// (Meeus also gives a series in powers of the eccentricity; it is not
/// provided because it converges more slowly than this one.)
///
/// # Returns
///
/// * `length`: Length of the ellipse (same unit as that of `a` and `b`)
///
/// # Arguments
///
/// * `a`: Semimajor axis of the ellipse, `a ≥ 0` (same unit as that of `b`)
/// * `b`: Semiminor axis of the ellipse, `b ≥ 0` (same unit as that of `a`)
///
/// # Example
///
/// The orbit of comet Halley (Meeus, example 33.d):
///
/// ```
/// use astro::orbit::elliptic;
///
/// let (a, e) = (17.9400782, 0.96727426);
/// let b: f64 = a * (1.0_f64 - e * e).sqrt();
/// let L = elliptic::length_series(a, b);
/// assert!((L - 77.07).abs() < 0.005);
/// ```
pub fn length_series(a: f64, b: f64) -> f64 {
    let s = a + b;
    if a == 0.0 || b == 0.0 {
        return 4.0 * s;
    }

    let m = (a - b) / s;
    let h = m * m;

    let mut sum = 1.0;
    let mut coeff = 1.0;
    let mut h_pow = 1.0;
    for n in 1..=100_000 {
        let n = n as f64;
        let ratio = (2.0 * n - 3.0) / (2.0 * n);
        coeff *= ratio * ratio;
        h_pow *= h;
        let term = coeff * h_pow;
        sum += term;
        if term.is_nan() || term <= 1e-17 * sum {
            break;
        }
    }

    PI * s * sum
}

/// Computes the semimajor axis of an elliptic orbit
///
/// Uses `a = q / (1 - e)` (Meeus, eq. 33.6).
///
/// # Returns
///
/// * `a`: Semimajor axis *| in the unit of `perih`* (AU)
///
/// # Arguments
///
/// * `perih`: Perihelion distance of the orbit *| in AU*
/// * `ecc`  : Eccentricity of the orbit, `0 ≤ ecc < 1`
#[inline(always)]
pub fn semimaj_axis(perih: f64, ecc: f64) -> f64 {
    perih / (1.0 - ecc)
}

/// Computes the perihelion distance of an elliptic orbit
///
/// Uses `q = a (1 - e)` (Meeus, chapter 30).
///
/// # Arguments
///
/// * `a`  : Semimajor axis of the orbit *| in AU*
/// * `ecc`: Eccentricity of the orbit
#[inline]
pub fn perih_dist(a: f64, ecc: f64) -> f64 {
    a * (1.0 - ecc)
}

/// Computes the aphelion distance of an elliptic orbit
///
/// Uses `Q = a (1 + e)` (Meeus, chapter 30).
///
/// # Arguments
///
/// * `a`  : Semimajor axis of the orbit *| in AU*
/// * `ecc`: Eccentricity of the orbit
#[inline]
pub fn aph_dist(a: f64, ecc: f64) -> f64 {
    a * (1.0 + ecc)
}

/// Computes the mean motion of an elliptic orbit
///
/// Uses `n = k / a^(3/2)` (Meeus, eq. 33.6), where `k` is the Gaussian
/// gravitational constant; the mass of the body is neglected.
///
/// # Returns
///
/// * `mean_motion`: Mean motion of the orbit
///   *| in radians per day of dynamical time*
///
/// # Arguments
///
/// * `semimaj_ax`: Semimajor axis of the orbit *| in AU*
#[inline(always)]
pub fn mn_motion(semimaj_ax: f64) -> f64 {
    consts::GAUSS_GRAV / semimaj_ax.powf(1.5)
}

/// Computes the time of passage of a body through a node of an elliptic
/// orbit, and its radius vector at that time
///
/// Meeus, chapter 39 (eqs. 39.1 to 39.4): the true anomaly at the node is
/// `-ω` (ascending) or `π - ω` (descending); it is converted to the
/// eccentric and then the mean anomaly, which gives the time.
///
/// The passage returned is the one closest to the time of perihelion `T`,
/// i.e. within half an orbital period of it (`|time_of_pass - T| ≤ π/n`).
/// Other passages follow at multiples of the period `2π/n`.
///
/// The node is on the plane of reference of the elements: for elements
/// referred to the equinox J2000.0, it is the node on the ecliptic of
/// J2000.0, not on the ecliptic of the date.
///
/// # Returns
///
/// `(time_of_pass, rad_vec)`
///
/// * `time_of_pass`: Time of passage through the node, in Julian
///   (Ephemeris) day
/// * `rad_vec`     : Radius vector of the body at the time of passage
///   *| in AU*
///
/// # Arguments
///
/// * `w`   : Argument of the perihelion *| in radians*
/// * `n`   : Mean motion of the orbit *| in radians per day*
/// * `a`   : Semimajor axis of the orbit *| in AU*
/// * `e`   : Eccentricity of the orbit, `0 ≤ e < 1`
/// * `T`   : Time of passage in perihelion, in Julian (Ephemeris) day
/// * `node`: `Ascend` or `Descend` node
///
/// # Example
///
/// Comet Halley at its ascending node in 1985 (Meeus, example 39.a):
///
/// ```
/// use astro::orbit::{self, elliptic};
///
/// let T = 2446470.95891; // 1986 Feb. 9.45891 TD
/// let (t, r) = elliptic::passage_through_node(
///     111.84644_f64.to_radians(),
///     0.01297082_f64.to_radians(),
///     17.9400782,
///     0.96727426,
///     T,
///     &orbit::Node::Ascend,
/// );
/// assert!((t - (T - 92.2998)).abs() < 1e-4);
/// assert!((r - 1.8045).abs() < 1e-4);
/// ```
#[inline]
pub fn passage_through_node(
    w: f64,
    n: f64,
    a: f64,
    e: f64,
    T: f64,
    node: &orbit::Node,
) -> (f64, f64) {
    match *node {
        orbit::Node::Ascend => pass_through_node(-w, n, a, e, T),
        orbit::Node::Descend => pass_through_node(PI - w, n, a, e, T),
    }
}

fn pass_through_node(v: f64, n: f64, a: f64, e: f64, T: f64) -> (f64, f64) {
    let E = ecc_anom_frm_true_anom(v, e);
    let M = mn_anom_frm_ecc_anom(E, e);

    (T + M / n, rad_vec_frm_ecc_anom(E, a, e))
}

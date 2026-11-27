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

//! Corrections for nutation
//!
//! Nutation in longitude (Δψ) and in obliquity (Δε) from the IAU 1980
//! theory, as given in Meeus, *Astronomical Algorithms*, 2nd ed.,
//! chapter 22 (terms of Table 22.A, i.e. all terms of 0.0003″ and more),
//! a low accuracy variant, and the effect of nutation on equatorial
//! coordinates (chapter 23).
//!
//! The mean and true obliquity of the ecliptic are in
//! [`ecliptic`](crate::ecliptic).
//!
//! # Example
//!
//! Meeus, example 22.a: 1987 April 10, 0h TD.
//!
//! ```
//! use astro::nutation;
//!
//! let (nut_in_long, nut_in_oblq) = nutation::nutation(2446895.5);
//! let arcsec = |x: f64| x.to_degrees() * 3600.0;
//! assert!((arcsec(nut_in_long) - -3.788).abs() < 0.0005);
//! assert!((arcsec(nut_in_oblq) - 9.443).abs() < 0.0005);
//! ```

use crate::coords;
use crate::time;
use std::f64::consts::{PI, TAU};

/// One arcsecond in radians
const ASEC: f64 = PI / 648000.0;

#[inline]
fn norm_two_pi(a: f64) -> f64 {
    let r = a.rem_euclid(TAU);
    if r >= TAU { 0.0 } else { r }
}

/// The five fundamental arguments of the IAU 1980 nutation theory
///
/// All are angles in radians, normalised to `[0, 2π)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NutationArgs {
    /// D: mean elongation of the Moon from the Sun
    pub mn_elong_moon: f64,
    /// M: mean anomaly of the Sun (Earth)
    pub mn_anom_sun: f64,
    /// M′: mean anomaly of the Moon
    pub mn_anom_moon: f64,
    /// F: the Moon's argument of latitude
    pub arg_lat_moon: f64,
    /// Ω: longitude of the ascending node of the Moon's mean orbit on the
    /// ecliptic, measured from the mean equinox of the date
    pub long_asc_node_moon: f64,
}

/// Computes the fundamental arguments D, M, M′, F and Ω of the nutation
///
/// Meeus, chapter 22 (cubic expressions in T, Julian centuries from
/// J2000.0). These differ slightly from the ones of Chapront's lunar
/// theory.
///
/// # Arguments
///
/// * `JD`: Julian Ephemeris Day
pub fn nutation_args(JD: f64) -> NutationArgs {
    let t = time::julian_cent(JD);

    let deg = |x: f64| norm_two_pi(x.rem_euclid(360.0).to_radians());

    NutationArgs {
        mn_elong_moon: deg(297.85036 + t * (445267.11148 - t * (0.0019142 - t / 189474.0))),
        mn_anom_sun: deg(357.52772 + t * (35999.05034 - t * (0.0001603 + t / 300000.0))),
        mn_anom_moon: deg(134.96298 + t * (477198.867398 + t * (0.0086972 + t / 56250.0))),
        arg_lat_moon: deg(93.27191 + t * (483202.017538 - t * (0.0036825 - t / 327270.0))),
        long_asc_node_moon: deg(125.04452 - t * (1934.136261 - t * (0.0020708 + t / 450000.0))),
    }
}

/// One periodic term of Table 22.A: the multiples of D, M, M′, F, Ω, then
/// the sine coefficient of Δψ and its rate per century, then the cosine
/// coefficient of Δε and its rate per century. Coefficients are in units
/// of 0.0001″, the rates in units of 0.00001″ per Julian century.
struct NutationTerm(i8, i8, i8, i8, i8, i32, i32, i32, i32);

#[rustfmt::skip]
const NUTATION_TERMS: [NutationTerm; 63] = [
    NutationTerm( 0,  0,  0,  0,  1, -171996, -1742, 92025,  89),
    NutationTerm(-2,  0,  0,  2,  2,  -13187,   -16,  5736, -31),
    NutationTerm( 0,  0,  0,  2,  2,   -2274,    -2,   977,  -5),
    NutationTerm( 0,  0,  0,  0,  2,    2062,     2,  -895,   5),
    NutationTerm( 0,  1,  0,  0,  0,    1426,   -34,    54,  -1),
    NutationTerm( 0,  0,  1,  0,  0,     712,     1,    -7,   0),
    NutationTerm(-2,  1,  0,  2,  2,    -517,    12,   224,  -6),
    NutationTerm( 0,  0,  0,  2,  1,    -386,    -4,   200,   0),
    NutationTerm( 0,  0,  1,  2,  2,    -301,     0,   129,  -1),
    NutationTerm(-2, -1,  0,  2,  2,     217,    -5,   -95,   3),
    NutationTerm(-2,  0,  1,  0,  0,    -158,     0,     0,   0),
    NutationTerm(-2,  0,  0,  2,  1,     129,     1,   -70,   0),
    NutationTerm( 0,  0, -1,  2,  2,     123,     0,   -53,   0),
    NutationTerm( 2,  0,  0,  0,  0,      63,     0,     0,   0),
    NutationTerm( 0,  0,  1,  0,  1,      63,     1,   -33,   0),
    NutationTerm( 2,  0, -1,  2,  2,     -59,     0,    26,   0),
    NutationTerm( 0,  0, -1,  0,  1,     -58,    -1,    32,   0),
    NutationTerm( 0,  0,  1,  2,  1,     -51,     0,    27,   0),
    NutationTerm(-2,  0,  2,  0,  0,      48,     0,     0,   0),
    NutationTerm( 0,  0, -2,  2,  1,      46,     0,   -24,   0),
    NutationTerm( 2,  0,  0,  2,  2,     -38,     0,    16,   0),
    NutationTerm( 0,  0,  2,  2,  2,     -31,     0,    13,   0),
    NutationTerm( 0,  0,  2,  0,  0,      29,     0,     0,   0),
    NutationTerm(-2,  0,  1,  2,  2,      29,     0,   -12,   0),
    NutationTerm( 0,  0,  0,  2,  0,      26,     0,     0,   0),
    NutationTerm(-2,  0,  0,  2,  0,     -22,     0,     0,   0),
    NutationTerm( 0,  0, -1,  2,  1,      21,     0,   -10,   0),
    NutationTerm( 0,  2,  0,  0,  0,      17,    -1,     0,   0),
    NutationTerm( 2,  0, -1,  0,  1,      16,     0,    -8,   0),
    NutationTerm(-2,  2,  0,  2,  2,     -16,     1,     7,   0),
    NutationTerm( 0,  1,  0,  0,  1,     -15,     0,     9,   0),
    NutationTerm(-2,  0,  1,  0,  1,     -13,     0,     7,   0),
    NutationTerm( 0, -1,  0,  0,  1,     -12,     0,     6,   0),
    NutationTerm( 0,  0,  2, -2,  0,      11,     0,     0,   0),
    NutationTerm( 2,  0, -1,  2,  1,     -10,     0,     5,   0),
    NutationTerm( 2,  0,  1,  2,  2,      -8,     0,     3,   0),
    NutationTerm( 0,  1,  0,  2,  2,       7,     0,    -3,   0),
    NutationTerm(-2,  1,  1,  0,  0,      -7,     0,     0,   0),
    NutationTerm( 0, -1,  0,  2,  2,      -7,     0,     3,   0),
    NutationTerm( 2,  0,  0,  2,  1,      -7,     0,     3,   0),
    NutationTerm( 2,  0,  1,  0,  0,       6,     0,     0,   0),
    NutationTerm(-2,  0,  2,  2,  2,       6,     0,    -3,   0),
    NutationTerm(-2,  0,  1,  2,  1,       6,     0,    -3,   0),
    NutationTerm( 2,  0, -2,  0,  1,      -6,     0,     3,   0),
    NutationTerm( 2,  0,  0,  0,  1,      -6,     0,     3,   0),
    NutationTerm( 0, -1,  1,  0,  0,       5,     0,     0,   0),
    NutationTerm(-2, -1,  0,  2,  1,      -5,     0,     3,   0),
    NutationTerm(-2,  0,  0,  0,  1,      -5,     0,     3,   0),
    NutationTerm( 0,  0,  2,  2,  1,      -5,     0,     3,   0),
    NutationTerm(-2,  0,  2,  0,  1,       4,     0,     0,   0),
    NutationTerm(-2,  1,  0,  2,  1,       4,     0,     0,   0),
    NutationTerm( 0,  0,  1, -2,  0,       4,     0,     0,   0),
    NutationTerm(-1,  0,  1,  0,  0,      -4,     0,     0,   0),
    NutationTerm(-2,  1,  0,  0,  0,      -4,     0,     0,   0),
    NutationTerm( 1,  0,  0,  0,  0,      -4,     0,     0,   0),
    NutationTerm( 0,  0,  1,  2,  0,       3,     0,     0,   0),
    NutationTerm( 0,  0, -2,  2,  2,      -3,     0,     0,   0),
    NutationTerm(-1, -1,  1,  0,  0,      -3,     0,     0,   0),
    NutationTerm( 0,  1,  1,  0,  0,      -3,     0,     0,   0),
    NutationTerm( 0, -1,  1,  2,  2,      -3,     0,     0,   0),
    NutationTerm( 2, -1, -1,  2,  2,      -3,     0,     0,   0),
    NutationTerm( 0,  0,  3,  2,  2,      -3,     0,     0,   0),
    NutationTerm( 2, -1,  0,  2,  2,      -3,     0,     0,   0),
];

/// Computes nutation in ecliptic longitude and obliquity
///
/// IAU 1980 theory of nutation, with the 63 terms of Meeus' Table 22.A
/// (terms smaller than 0.0003″ are neglected). Valid for several
/// centuries around J2000.0.
///
/// # Returns
///
/// `(nut_in_long, nut_in_oblq)`
///
/// * `nut_in_long`: Nutation in ecliptic longitude, Δψ *| in radians*
/// * `nut_in_oblq`: Nutation in obliquity of the ecliptic, Δε
///   *| in radians*
///
/// # Arguments
///
/// `JD`: Julian Ephemeris Day
pub fn nutation(JD: f64) -> (f64, f64) {
    let t = time::julian_cent(JD);
    let a = nutation_args(JD);

    let mut nut_in_long = 0.0;
    let mut nut_in_oblq = 0.0;

    // Sum the smallest terms first
    for x in NUTATION_TERMS.iter().rev() {
        let arg = (x.0 as f64) * a.mn_elong_moon
            + (x.1 as f64) * a.mn_anom_sun
            + (x.2 as f64) * a.mn_anom_moon
            + (x.3 as f64) * a.arg_lat_moon
            + (x.4 as f64) * a.long_asc_node_moon;

        let (sin_arg, cos_arg) = arg.sin_cos();

        nut_in_long += ((x.5 as f64) + t * (x.6 as f64) / 10.0) * sin_arg;
        nut_in_oblq += ((x.7 as f64) + t * (x.8 as f64) / 10.0) * cos_arg;
    }

    (nut_in_long * 0.0001 * ASEC, nut_in_oblq * 0.0001 * ASEC)
}

/// Computes nutation in ecliptic longitude and obliquity, with low
/// accuracy
///
/// Meeus, chapter 22: the four largest terms, in the mean longitudes of
/// the Sun (L) and the Moon (L′) and the longitude of the Moon's node Ω.
/// The error is below about 0.5″ in Δψ and 0.1″ in Δε.
///
/// (In the 2nd edition the first formula is printed with a doubled minus
/// sign in front of the `sin 2L` term; the term is negative.)
///
/// # Returns
///
/// `(nut_in_long, nut_in_oblq)` *| in radians*
///
/// # Arguments
///
/// `JD`: Julian Ephemeris Day
pub fn nutation_low_accuracy(JD: f64) -> (f64, f64) {
    let t = time::julian_cent(JD);

    let om = (125.04452 - 1934.136261 * t).to_radians();
    let L = (280.4665 + 36000.7698 * t).to_radians();
    let L1 = (218.3165 + 481267.8813 * t).to_radians();

    let nut_in_long = -17.20 * om.sin() - 1.32 * (2.0 * L).sin() - 0.23 * (2.0 * L1).sin()
        + 0.21 * (2.0 * om).sin();
    let nut_in_oblq = 9.20 * om.cos() + 0.57 * (2.0 * L).cos() + 0.10 * (2.0 * L1).cos()
        - 0.09 * (2.0 * om).cos();

    (nut_in_long * ASEC, nut_in_oblq * ASEC)
}

/// Computes nutation in equatorial coordinates
///
/// First-order corrections of Meeus, eq. 23.1, to add to a right
/// ascension and declination referred to the mean equator and equinox of
/// the date. They contain `tan(dec)` and are not valid close to the
/// celestial poles; use [`nutate_eq_coords`] there.
///
/// # Returns
///
/// `(nut_in_asc, nut_in_dec)`
///
/// * `nut_in_asc`: Nutation in right ascension *| in radians*
/// * `nut_in_dec`: Nutation in declination *| in radians*
///
/// # Arguments
///
/// * `eq_point`   : Equatorial point uncorrected for nutation *| in radians*
/// * `nut_in_long`: Nutation in longitude *| in radians*
/// * `nut_in_oblq`: Nutation in obliquity *| in radians*
/// * `tru_oblq`   : Obliquity of the ecliptic *| in radians* (mean or
///   true: the difference is negligible here)
pub fn nutation_in_eq_coords(
    eq_point: &coords::EqPoint,
    nut_in_long: f64,
    nut_in_oblq: f64,
    tru_oblq: f64,
) -> (f64, f64) {
    let (asc, dec) = (eq_point.asc, eq_point.dec);
    let (sin_a, cos_a) = asc.sin_cos();
    let (sin_e, cos_e) = tru_oblq.sin_cos();
    let tan_d = dec.tan();

    let nut_asc = nut_in_long * (cos_e + sin_e * sin_a * tan_d) - cos_a * tan_d * nut_in_oblq;

    let nut_dec = sin_e * cos_a * nut_in_long + sin_a * nut_in_oblq;

    (nut_asc, nut_dec)
}

/// Computes equatorial coordinates corrected for nutation, rigorously
///
/// Rotates a position referred to the mean equator and equinox of the
/// date into the true equator and equinox of the date: the position is
/// turned into ecliptical coordinates with the mean obliquity, Δψ is
/// added to the longitude, and the result is turned back into equatorial
/// coordinates with the true obliquity `mn_oblq + nut_in_oblq`. This is
/// the approach Meeus recommends (chapter 23) for stars close to the
/// celestial poles; elsewhere it agrees with [`nutation_in_eq_coords`] to
/// first order.
///
/// # Returns
///
/// `(asc, dec)`: right ascension in `[0, 2π)` and declination referred to
/// the true equator and equinox *| in radians*
///
/// # Arguments
///
/// * `asc`        : Right ascension, mean equinox of the date *| in radians*
/// * `dec`        : Declination, mean equator of the date *| in radians*
/// * `nut_in_long`: Nutation in longitude Δψ *| in radians*
/// * `nut_in_oblq`: Nutation in obliquity Δε *| in radians*
/// * `mn_oblq`    : Mean obliquity of the ecliptic *| in radians*
pub fn nutate_eq_coords(
    asc: f64,
    dec: f64,
    nut_in_long: f64,
    nut_in_oblq: f64,
    mn_oblq: f64,
) -> (f64, f64) {
    let (sin_a, cos_a) = asc.sin_cos();
    let (sin_d, cos_d) = dec.sin_cos();

    // Unit vector, mean equatorial frame
    let (x, y, z) = (cos_d * cos_a, cos_d * sin_a, sin_d);

    // To the ecliptic frame (rotation by the mean obliquity)
    let (sin_e0, cos_e0) = mn_oblq.sin_cos();
    let (x1, y1, z1) = (x, y * cos_e0 + z * sin_e0, -y * sin_e0 + z * cos_e0);

    // Add the nutation in longitude
    let (sin_p, cos_p) = nut_in_long.sin_cos();
    let (x2, y2) = (x1 * cos_p - y1 * sin_p, x1 * sin_p + y1 * cos_p);

    // Back to the equatorial frame with the true obliquity
    let (sin_e, cos_e) = (mn_oblq + nut_in_oblq).sin_cos();
    let (X, Y, Z) = (x2, y2 * cos_e - z1 * sin_e, y2 * sin_e + z1 * cos_e);

    (norm_two_pi(Y.atan2(X)), Z.atan2(X.hypot(Y)))
}

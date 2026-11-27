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

//! Corrections for aberration
//!
//! Annual aberration of the stars (Meeus, *Astronomical Algorithms*,
//! 2nd ed., chapter 23), by two methods:
//!
//! * the classical formulae 23.2 and 23.3, from the Sun's longitude and
//!   the eccentricity and perihelion of the Earth's orbit, in the frame of
//!   the date ([`annual_aberr_in_eq_coords`], [`annual_aberr_in_ecl_coords`],
//!   and their FK4 variants without the E-terms);
//! * the Ron–Vondrák expression for the velocity of the Earth relative to
//!   the barycentre of the solar system (Table 23.A), which gives the
//!   aberration in the fixed J2000.0 frame ([`stell_aberr_in_eq_coords`],
//!   [`earth_bary_vel`]).
//!
//! Also the aberration of the Sun ([`sol_aberr`], chapter 25).
//!
//! # Example
//!
//! Meeus, example 23.b: θ Persei on 2028 November 13.19 TD, position
//! referred to J2000.0 and corrected for proper motion.
//!
//! ```
//! use astro::{aberr, coords};
//!
//! let star = coords::EqPoint {
//!     asc: 41.0540613_f64.to_radians(),
//!     dec: 49.2277489_f64.to_radians(),
//! };
//! let (d_asc, d_dec) = aberr::stell_aberr_in_eq_coords(&star, 2462088.69);
//! assert!((d_asc - 0.000145252).abs() < 1e-9);
//! assert!((d_dec - 0.000032723).abs() < 1e-9);
//! ```

use crate::coords;
use crate::time;
use std::f64::consts::PI;

/// One arcsecond in radians
const ASEC: f64 = PI / 648000.0;

/// The constant of aberration κ = 20.49552″ *| in radians*
pub const CONST_OF_ABERR: f64 = 20.49552 * ASEC;

/// The speed of light in AU per day, in the units of [`earth_bary_vel`]
/// (Meeus, chapter 23: 17 314 463 350 × 10⁻⁸ AU/day)
pub(crate) const LIGHT_SPEED_AU_PER_DAY: f64 = 173.14463350;

/// Computes solar aberration in ecliptic longitude
///
/// Meeus, eq. 25.10: the aberration of the Sun, −20.4898″ / R. Add it
/// to the Sun's geometric longitude to get its apparent longitude (with
/// the nutation in longitude).
///
/// # Returns
///
/// * `abrr`: Solar aberration in ecliptic
///   longitude *| in radians* (negative)
///
/// # Arguments
///
/// * `R`: Sun-Earth distance *| in AU*
#[inline(always)]
pub fn sol_aberr(R: f64) -> f64 {
    -20.4898 * ASEC / R
}

/// Computes the eccentricity of the Earth's orbit
///
/// Meeus, chapter 23 (also eq. 25.4).
///
/// # Arguments
///
/// * `JD`: Julian Ephemeris Day
#[inline]
pub fn earth_orb_ecc(JD: f64) -> f64 {
    crate::sun::earth_orbit_eccentricity(JD)
}

/// Computes the longitude of the perihelion of the Earth's orbit
///
/// Meeus, chapter 23, referred to the mean equinox of the date.
///
/// # Returns
///
/// * `long`: Longitude of the perihelion *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `JD`: Julian Ephemeris Day
#[inline]
pub fn earth_perih_long(JD: f64) -> f64 {
    let t = time::julian_cent(JD);

    (102.93735 + t * (1.71946 + 0.00046 * t))
        .rem_euclid(360.0)
        .to_radians()
}

/// The velocity of the Earth in its unperturbed elliptical orbit, in units
/// of the speed of light, as `[x, y, z]` in the ecliptic frame of the date
/// (x towards the mean equinox, z towards the north ecliptic pole).
///
/// This is the velocity behind eqs. 23.2 and 23.3: its magnitude is the
/// constant of aberration κ, and the terms in `e` are the E-terms.
/// `sun_tru_long` is the Sun's true geometric longitude.
pub(crate) fn earth_vel_classical(sun_tru_long: f64, JD: f64) -> [f64; 3] {
    let k = CONST_OF_ABERR;
    let e = earth_orb_ecc(JD);
    let (sin_p, cos_p) = earth_perih_long(JD).sin_cos();
    let (sin_s, cos_s) = sun_tru_long.sin_cos();

    [k * (sin_s - e * sin_p), k * (e * cos_p - cos_s), 0.0]
}

/// Displaces the direction `(lon, lat)` by the aberration for an observer
/// moving with velocity `v` (in units of the speed of light, same frame):
/// the apparent direction is that of `u + v`, `u` being the unit vector of
/// the true direction.
///
/// To first order this is the same as eqs. 23.2–23.4, but it involves no
/// `1/cos(lat)` and so holds at the poles too. Returns the longitude in
/// `[0, 2π)` and the latitude.
pub(crate) fn aberrate(lon: f64, lat: f64, v: [f64; 3]) -> (f64, f64) {
    let (sin_b, cos_b) = lat.sin_cos();
    let (sin_l, cos_l) = lon.sin_cos();
    let x = cos_b * cos_l + v[0];
    let y = cos_b * sin_l + v[1];
    let z = sin_b + v[2];

    let r = y.atan2(x).rem_euclid(2.0 * PI);
    (if r >= 2.0 * PI { 0.0 } else { r }, z.atan2(x.hypot(y)))
}

/// The factors of κ in eq. 23.3 for a given longitude `l` (the Sun's or
/// the perihelion's): (for Δα, for Δδ).
fn aberr_eq_factors(asc: f64, dec: f64, l: f64, oblq: f64) -> (f64, f64) {
    let (sin_a, cos_a) = asc.sin_cos();
    let (sin_d, cos_d) = dec.sin_cos();
    let (sin_l, cos_l) = l.sin_cos();
    let (sin_e, cos_e) = oblq.sin_cos();

    let f_asc = (cos_a * cos_l * cos_e + sin_a * sin_l) / cos_d;
    let f_dec = cos_l * (sin_e * cos_d - cos_e * sin_a * sin_d) + cos_a * sin_d * sin_l;

    (f_asc, f_dec)
}

/// Computes the annual aberration of a star in equatorial coordinates
///
/// Classical formulae of Meeus, eq. 23.3, including the terms in the
/// eccentricity of the Earth's orbit (the E-terms). Use this for FK5
/// positions (and catalogues based on the FK5); for FK4 positions use
/// [`annual_aberr_in_eq_coords_fk4`]. The corrections are to be added to
/// the right ascension and declination referred to the equinox of the
/// date.
///
/// The right ascension correction contains `1/cos(dec)` and is not
/// finite at the celestial poles.
///
/// # Returns
///
/// `(abrr_in_asc, abrr_in_dec)` *| in radians*
///
/// # Arguments
///
/// * `asc`         : Right ascension of the star *| in radians*
/// * `dec`         : Declination of the star *| in radians*
/// * `sun_tru_long`: True (geometric) longitude of the Sun, mean equinox
///   of the date *| in radians*. An accuracy of 0.01° is enough.
/// * `oblq`        : Obliquity of the ecliptic *| in radians*
/// * `JD`          : Julian Ephemeris Day (for the eccentricity and the
///   perihelion of the Earth's orbit)
pub fn annual_aberr_in_eq_coords(
    asc: f64,
    dec: f64,
    sun_tru_long: f64,
    oblq: f64,
    JD: f64,
) -> (f64, f64) {
    let k = CONST_OF_ABERR;
    let e = earth_orb_ecc(JD);

    let (sun_a, sun_d) = aberr_eq_factors(asc, dec, sun_tru_long, oblq);
    let (per_a, per_d) = aberr_eq_factors(asc, dec, earth_perih_long(JD), oblq);

    (-k * sun_a + e * k * per_a, -k * sun_d + e * k * per_d)
}

/// Computes the annual aberration of a star in equatorial coordinates,
/// for positions referred to the FK4 system
///
/// Same as [`annual_aberr_in_eq_coords`] without the terms in the
/// eccentricity of the Earth's orbit (the E-terms), which are already
/// included in FK4 mean places (Meeus, chapters 21 and 23).
///
/// # Returns
///
/// `(abrr_in_asc, abrr_in_dec)` *| in radians*
///
/// # Arguments
///
/// * `asc`         : Right ascension of the star *| in radians*
/// * `dec`         : Declination of the star *| in radians*
/// * `sun_tru_long`: True (geometric) longitude of the Sun *| in radians*
/// * `oblq`        : Obliquity of the ecliptic *| in radians*
pub fn annual_aberr_in_eq_coords_fk4(
    asc: f64,
    dec: f64,
    sun_tru_long: f64,
    oblq: f64,
) -> (f64, f64) {
    let (sun_a, sun_d) = aberr_eq_factors(asc, dec, sun_tru_long, oblq);

    (-CONST_OF_ABERR * sun_a, -CONST_OF_ABERR * sun_d)
}

/// Computes the annual aberration of a star in ecliptical coordinates
///
/// Meeus, eq. 23.2, including the E-terms (use
/// [`annual_aberr_in_ecl_coords_fk4`] for FK4 positions). The latitude
/// enters through `1/cos(lat)`; the longitude correction is not finite at
/// the poles of the ecliptic.
///
/// # Returns
///
/// `(abrr_in_long, abrr_in_lat)` *| in radians*
///
/// # Arguments
///
/// * `ecl_long`    : Ecliptic longitude of the star *| in radians*
/// * `ecl_lat`     : Ecliptic latitude of the star *| in radians*
/// * `sun_tru_long`: True (geometric) longitude of the Sun *| in radians*
/// * `JD`          : Julian Ephemeris Day
pub fn annual_aberr_in_ecl_coords(
    ecl_long: f64,
    ecl_lat: f64,
    sun_tru_long: f64,
    JD: f64,
) -> (f64, f64) {
    let k = CONST_OF_ABERR;
    let e = earth_orb_ecc(JD);
    let per = earth_perih_long(JD);

    let (sin_s, cos_s) = (sun_tru_long - ecl_long).sin_cos();
    let (sin_p, cos_p) = (per - ecl_long).sin_cos();

    (
        (-k * cos_s + e * k * cos_p) / ecl_lat.cos(),
        -k * ecl_lat.sin() * (sin_s - e * sin_p),
    )
}

/// Computes the annual aberration of a star in ecliptical coordinates,
/// for positions referred to the FK4 system
///
/// Same as [`annual_aberr_in_ecl_coords`] without the E-terms.
///
/// # Returns
///
/// `(abrr_in_long, abrr_in_lat)` *| in radians*
///
/// # Arguments
///
/// * `ecl_long`    : Ecliptic longitude of the star *| in radians*
/// * `ecl_lat`     : Ecliptic latitude of the star *| in radians*
/// * `sun_tru_long`: True (geometric) longitude of the Sun *| in radians*
pub fn annual_aberr_in_ecl_coords_fk4(
    ecl_long: f64,
    ecl_lat: f64,
    sun_tru_long: f64,
) -> (f64, f64) {
    let k = CONST_OF_ABERR;
    let (sin_s, cos_s) = (sun_tru_long - ecl_long).sin_cos();

    (-k * cos_s / ecl_lat.cos(), -k * ecl_lat.sin() * sin_s)
}

/// Computes the velocity of the Earth relative to the barycentre of the
/// solar system
///
/// Ron–Vondrák expansion (Meeus, chapter 23, Table 23.A, 36 terms in the
/// mean longitudes of Venus to Neptune, of the Earth and of the Moon).
/// The components are in the equatorial frame of J2000.0 (FK5), not of
/// the date.
///
/// # Returns
///
/// `(x, y, z)`: velocity components *| in AU per day*
///
/// # Arguments
///
/// * `JD`: Julian Ephemeris Day
pub fn earth_bary_vel(JD: f64) -> (f64, f64, f64) {
    let t = time::julian_cent(JD);

    let l2 = 3.1761467 + 1021.3285546 * t;
    let l3 = 1.7534703 + 628.3075849 * t;
    let l4 = 6.2034809 + 334.0612431 * t;
    let l5 = 0.5995465 + 52.9690965 * t;
    let l6 = 0.8740168 + 21.3299095 * t;
    let l7 = 5.4812939 + 7.4781599 * t;
    let l8 = 5.3118863 + 3.8133036 * t;
    let l1 = 3.8103444 + 8399.6847337 * t;
    let d = 5.1984667 + 7771.3771486 * t;
    let m1 = 2.3555559 + 8328.6914289 * t;
    let f = 1.6279052 + 8433.4661601 * t;

    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;

    let mut A;
    let mut sinA;
    let mut cosA;

    // ROW 1
    A = l3;
    sinA = A.sin();
    cosA = A.cos();
    x += (-1719914.0 - 2.0 * t) * sinA - 25.0 * cosA;
    y += (25.0 - 13.0 * t) * sinA + (1578089.0 + 156.0 * t) * cosA;
    z += (10.0 + 32.0 * t) * sinA + (684185.0 - 358.0 * t) * cosA;

    // ROW 2
    A = 2.0 * l3;
    sinA = A.sin();
    cosA = A.cos();
    x += (6434.0 + 141.0 * t) * sinA + (28007.0 - 107.0 * t) * cosA;
    y += (25697.0 - 95.0 * t) * sinA + (-5904.0 - 130.0 * t) * cosA;
    z += (11141.0 - 48.0 * t) * sinA + (-2559.0 - 55.0 * t) * cosA;

    // ROW 3
    A = l5;
    sinA = A.sin();
    cosA = A.cos();
    x += 715.0 * sinA;
    y += 6.0 * sinA - 657.0 * cosA;
    z += -15.0 * sinA - 282.0 * cosA;

    // ROW 4
    A = l1;
    sinA = A.sin();
    cosA = A.cos();
    x += 715.0 * sinA;
    y += -656.0 * cosA;
    z += -285.0 * cosA;

    // ROW 5
    A = 3.0 * l3;
    sinA = A.sin();
    cosA = A.cos();
    x += (486.0 - 5.0 * t) * sinA + (-236.0 - 4.0 * t) * cosA;
    y += (-216.0 - 4.0 * t) * sinA + (-446.0 + 5.0 * t) * cosA;
    z += -94.0 * sinA - 193.0 * cosA;

    // ROW 6
    A = l6;
    sinA = A.sin();
    cosA = A.cos();
    x += 159.0 * sinA;
    y += 2.0 * sinA - 147.0 * cosA;
    z += -6.0 * sinA - 61.0 * cosA;

    // ROW 7
    A = f;
    cosA = A.cos();
    y += 26.0 * cosA;
    z += -59.0 * cosA;

    // ROW 8
    A = l1 + m1;
    sinA = A.sin();
    cosA = A.cos();
    x += 39.0 * sinA;
    y += -36.0 * cosA;
    z += -16.0 * cosA;

    // ROW 9
    A = 2.0 * l5;
    sinA = A.sin();
    cosA = A.cos();
    x += 33.0 * sinA - 10.0 * cosA;
    y += -9.0 * sinA - 30.0 * cosA;
    z += -5.0 * sinA - 13.0 * cosA;

    // ROW 10
    A = 2.0 * l3 - l5;
    sinA = A.sin();
    cosA = A.cos();
    x += 31.0 * sinA + cosA;
    y += sinA - 28.0 * cosA;
    z += -12.0 * cosA;

    // ROW 11
    A = 3.0 * l3 - 8.0 * l4 + 3.0 * l5;
    sinA = A.sin();
    cosA = A.cos();
    x += 8.0 * sinA - 28.0 * cosA;
    y += 25.0 * sinA + 8.0 * cosA;
    z += 11.0 * sinA + 3.0 * cosA;

    // ROW 12
    A = 5.0 * l3 - 8.0 * l4 + 3.0 * l5;
    sinA = A.sin();
    cosA = A.cos();
    x += 8.0 * sinA - 28.0 * cosA;
    y += -25.0 * sinA - 8.0 * cosA;
    z += -11.0 * sinA - 3.0 * cosA;

    // ROW 13
    A = 2.0 * l2 - l3;
    sinA = A.sin();
    cosA = A.cos();
    x += 21.0 * sinA;
    y += -19.0 * cosA;
    z += -8.0 * cosA;

    // ROW 14
    A = l2;
    sinA = A.sin();
    cosA = A.cos();
    x += -19.0 * sinA;
    y += 17.0 * cosA;
    z += 8.0 * cosA;

    // ROW 15
    A = l7;
    sinA = A.sin();
    cosA = A.cos();
    x += 17.0 * sinA;
    y += -16.0 * cosA;
    z += -7.0 * cosA;

    // ROW 16
    A = l3 - 2.0 * l5;
    sinA = A.sin();
    cosA = A.cos();
    x += 16.0 * sinA;
    y += 15.0 * cosA;
    z += sinA + 7.0 * cosA;

    // ROW 17
    A = l8;
    sinA = A.sin();
    cosA = A.cos();
    x += 16.0 * sinA;
    y += sinA - 15.0 * cosA;
    z += -3.0 * sinA - 6.0 * cosA;

    // ROW 18
    A = l3 + l5;
    sinA = A.sin();
    cosA = A.cos();
    x += 11.0 * sinA - cosA;
    y += -sinA - 10.0 * cosA;
    z += -sinA - 5.0 * cosA;

    // ROW 19
    A = 2.0 * (l2 - l3);
    sinA = A.sin();
    cosA = A.cos();
    x += -11.0 * cosA;
    y += -10.0 * sinA;
    z += -4.0 * sinA;

    // ROW 20
    A = l3 - l5;
    sinA = A.sin();
    cosA = A.cos();
    x += -11.0 * sinA - 2.0 * cosA;
    y += -2.0 * sinA + 9.0 * cosA;
    z += -sinA + 4.0 * cosA;

    // ROW 21
    A = 4.0 * l3;
    sinA = A.sin();
    cosA = A.cos();
    x += -7.0 * sinA - 8.0 * cosA;
    y += -8.0 * sinA + 6.0 * cosA;
    z += 3.0 * (cosA - sinA);

    // ROW 22
    A = 3.0 * l3 - 2.0 * l5;
    sinA = A.sin();
    cosA = A.cos();
    x += -10.0 * sinA;
    y += 9.0 * cosA;
    z += 4.0 * cosA;

    // ROW 23
    A = l2 - 2.0 * l3;
    sinA = A.sin();
    cosA = A.cos();
    x += -9.0 * sinA;
    y += -9.0 * cosA;
    z += -4.0 * cosA;

    // ROW 24
    A = 2.0 * l2 - 3.0 * l3;
    sinA = A.sin();
    cosA = A.cos();
    x += -9.0 * sinA;
    y += -8.0 * cosA;
    z += -4.0 * cosA;

    // ROW 25
    A = 2.0 * l6;
    sinA = A.sin();
    cosA = A.cos();
    x += -9.0 * cosA;
    y += -8.0 * sinA;
    z += -3.0 * sinA;

    // ROW 26
    A = 2.0 * l2 - 4.0 * l3;
    sinA = A.sin();
    cosA = A.cos();
    x += -9.0 * cosA;
    y += 8.0 * sinA;
    z += 3.0 * sinA;

    // ROW 27
    A = 3.0 * l3 - 2.0 * l4;
    sinA = A.sin();
    cosA = A.cos();
    x += 8.0 * sinA;
    y += -8.0 * cosA;
    z += -3.0 * cosA;

    // ROW 28
    A = l1 + 2.0 * d - m1;
    sinA = A.sin();
    cosA = A.cos();
    x += 8.0 * sinA;
    y += -7.0 * cosA;
    z += -3.0 * cosA;

    // ROW 29
    A = 8.0 * l2 - 12.0 * l3;
    sinA = A.sin();
    cosA = A.cos();
    x += -4.0 * sinA - 7.0 * cosA;
    y += -6.0 * sinA + 4.0 * cosA;
    z += -3.0 * sinA + 2.0 * cosA;

    // ROW 30
    A = 8.0 * l2 - 14.0 * l3;
    sinA = A.sin();
    cosA = A.cos();
    x += -4.0 * sinA - 7.0 * cosA;
    y += 6.0 * sinA - 4.0 * cosA;
    z += 3.0 * sinA - 2.0 * cosA;

    // ROW 31
    A = 2.0 * l4;
    sinA = A.sin();
    cosA = A.cos();
    x += -6.0 * sinA - 5.0 * cosA;
    y += -4.0 * sinA + 5.0 * cosA;
    z += 2.0 * (cosA - sinA);

    // ROW 32
    A = 3.0 * l2 - 4.0 * l3;
    sinA = A.sin();
    cosA = A.cos();
    x += -(sinA + cosA);
    y += -2.0 * sinA - 7.0 * cosA;
    z += sinA - 4.0 * cosA;

    // ROW 33
    A = 2.0 * (l3 - l5);
    sinA = A.sin();
    cosA = A.cos();
    x += 4.0 * sinA - 6.0 * cosA;
    y += -5.0 * sinA - 4.0 * cosA;
    z += -2.0 * (sinA + cosA);

    // ROW 34
    A = 3.0 * (l2 - l3);
    sinA = A.sin();
    cosA = A.cos();
    x += -7.0 * cosA;
    y += -6.0 * sinA;
    z += -3.0 * sinA;

    // ROW 35
    A = 2.0 * (l3 - l4);
    sinA = A.sin();
    cosA = A.cos();
    x += 5.0 * (sinA - cosA);
    y += -4.0 * sinA - 5.0 * cosA;
    z += -2.0 * (sinA + cosA);

    // ROW 36
    A = l1 - 2.0 * d;
    sinA = A.sin();
    cosA = A.cos();
    x += 5.0 * sinA;
    y += -5.0 * cosA;
    z += -2.0 * cosA;

    // Table 23.A is in units of 1e-8 AU per day
    (x * 1e-8, y * 1e-8, z * 1e-8)
}

/// Computes stellar aberration in equatorial coordinates
///
/// Ron–Vondrák method (Meeus, eq. 23.4), from the velocity of the Earth
/// relative to the barycentre of the solar system ([`earth_bary_vel`]).
/// More accurate than [`annual_aberr_in_eq_coords`].
///
/// The star's position and the corrections are referred to the fixed
/// equator and equinox of **J2000.0** (FK5), not to the equinox of the
/// date. So apply them to the J2000.0 position (after proper motion) and
/// *before* precession and nutation: FK5 position → proper motion →
/// aberration → precession → nutation.
///
/// # Returns
///
/// `(abrr_in_asc, abrr_in_dec)`
///
/// * `abrr_in_asc`: Aberration in right ascension *| in radians*
/// * `abrr_in_dec`: Aberration in declination *| in radians*
///
/// # Arguments
///
/// * `stell_eq_point`: Equatorial coordinates of the star, J2000.0
///   *| in radians*
/// * `JD`            : Julian Ephemeris Day
pub fn stell_aberr_in_eq_coords(stell_eq_point: &coords::EqPoint, JD: f64) -> (f64, f64) {
    let (x, y, z) = earth_bary_vel(JD);

    let c = LIGHT_SPEED_AU_PER_DAY;

    let (sin_a, cos_a) = stell_eq_point.asc.sin_cos();
    let (sin_d, cos_d) = stell_eq_point.dec.sin_cos();

    let delta_asc = (y * cos_a - x * sin_a) / (c * cos_d);
    let delta_dec = -((x * cos_a + y * sin_a) * sin_d - z * cos_d) / c;

    (delta_asc, delta_dec)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Eqs. 23.2 are the first-order expansion of `aberrate` with the
    // velocity of `earth_vel_classical`.
    #[test]
    fn vector_form_matches_eq_23_2() {
        let JD = 2462088.69;
        let sun = 231.328_f64.to_radians(); // Meeus, example 23.a
        for (l, b) in [(0.3, 0.2), (2.0, -1.0), (4.0, 1.3), (5.5, 0.0)] {
            let (dl, db) = annual_aberr_in_ecl_coords(l, b, sun, JD);
            let (l1, b1) = aberrate(l, b, earth_vel_classical(sun, JD));
            let d_l = (l1 - l + PI).rem_euclid(2.0 * PI) - PI;
            // second-order terms: κ² ≈ 1e-8 rad
            assert!((d_l - dl).abs() * b.cos() < 2e-8, "{l} {b}");
            assert!((b1 - b - db).abs() < 2e-8, "{l} {b}");
        }
    }

    #[test]
    fn aberrate_at_the_pole() {
        let v = [1e-4, -2e-4, 0.5e-4];
        let (l, b) = aberrate(0.7, PI / 2.0, v);
        assert!((0.0..2.0 * PI).contains(&l));
        // displaced towards the tangential part of v, by about its length
        let tangential = v[0].hypot(v[1]);
        assert!((PI / 2.0 - b - tangential.atan2(1.0 + v[2])).abs() < 1e-15);
        assert!((l - v[1].atan2(v[0]).rem_euclid(2.0 * PI)).abs() < 1e-12);
    }
}

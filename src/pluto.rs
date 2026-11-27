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

//! Pluto
//!
//! Heliocentric and geocentric positions of Pluto from the periodic-term
//! representation fitted to a numerical integration (Meeus, *Astronomical
//! Algorithms*, 2nd ed., chapter 37), plus Pluto's mean orbital elements
//! near 2000 and its visual magnitude (chapter 41).
//!
//! The periodic terms are only meaningful for the years **1885 to 2099**.
//! [`heliocent_pos`] does not check its argument and extrapolates the series
//! outside that span, where the result quickly becomes meaningless;
//! [`heliocent_pos_checked`] and [`geocent_astrometric_pos`] return `None`
//! instead. Use [`in_valid_range`] to test a date.
//!
//! # Example
//!
//! Pluto on 1992 October 13.0 TD (Meeus, example 37.a):
//!
//! ```
//! use astro::pluto;
//!
//! let jde = 2448908.5;
//! let (l, b, r) = pluto::heliocent_pos_checked(jde).unwrap();
//! assert!((l.to_degrees() - 232.74071).abs() < 1e-5);
//! assert!((b.to_degrees() - 14.58782).abs() < 1e-5);
//! assert!((r - 29.711111).abs() < 1e-6);
//!
//! // Geocentric rectangular equatorial coordinates of the Sun (J2000.0)
//! // for the same instant, e.g. from `sun` (Meeus, chapter 26).
//! let pos = pluto::geocent_astrometric_pos(jde, -0.9373959, -0.3131679, -0.1357792).unwrap();
//! assert!((pos.asc.to_degrees() - 232.93231).abs() < 1e-5);
//! assert!((pos.dec.to_degrees() - -4.45802).abs() < 1e-5);
//!
//! // Outside 1885-2099 the series is not valid.
//! assert!(pluto::heliocent_pos_checked(2299160.5).is_none());
//! ```

use crate::angle;
use crate::time;
use std::f64::consts::TAU;

/// Computes the geocentric equatorial semidiameter of Pluto
///
/// # Returns
///
/// * `semidiameter`: Geocentric equatorial semidiameter
///   *| in radians per AU*
///
/// # Arguments
///
/// * `pluto_earth_dist`: Pluto-Earth distance *| in AU*
#[deprecated(since = "3.0.0", note = "misspelt; use `pluto::semidiameter`")]
#[inline]
pub fn semdiameter(pluto_earth_dist: f64) -> f64 {
    semidiameter(pluto_earth_dist)
}

/// Computes the apparent semidiameter of Pluto
///
/// Uses a semidiameter of 2.07″ at a distance of 1 AU (Meeus, chapter 55).
///
/// # Returns
///
/// * `semidiameter`: Apparent semidiameter of Pluto *| in radians*
///
/// # Arguments
///
/// * `pluto_earth_dist`: Pluto-Earth distance *| in AU*
#[inline]
pub fn semidiameter(pluto_earth_dist: f64) -> f64 {
    angle::deg_frm_dms(0, 0, 2.07).to_radians() / pluto_earth_dist
}

/// Computes the apparent visual magnitude of Pluto
///
/// Uses the expression adopted by the *Astronomical Almanac* in 1984
/// (Meeus, chapter 41): `-1.00 + 5 log10(r Δ)`. The phase angle of Pluto
/// is always small and does not enter the formula.
///
/// # Returns
///
/// * `app_mag`: Apparent visual magnitude of Pluto
///
/// # Arguments
///
/// * `delta`: Pluto-Earth distance *| in AU*
/// * `r`    : Pluto-Sun distance *| in AU*
#[inline]
pub fn apprnt_mag_84(delta: f64, r: f64) -> f64 {
    (r * delta).log10() * 5.0 - 1.0
}

/// First Julian Ephemeris Day for which the Pluto series is valid
/// (1885 January 1.0 TD)
pub const VALID_FROM_JDE: f64 = 2409542.5;

/// Julian Ephemeris Day at which the Pluto series stops being valid
/// (2100 January 1.0 TD, i.e. the end of 2099)
pub const VALID_UNTIL_JDE: f64 = 2488069.5;

/// Checks whether a date lies in the span 1885–2099 covered by the Pluto
/// series
///
/// Returns `true` for `VALID_FROM_JDE <= JD < VALID_UNTIL_JDE`, and `false`
/// otherwise (including for NaN).
///
/// # Arguments
///
/// * `JD`: Julian Ephemeris Day
#[inline]
pub fn in_valid_range(JD: f64) -> bool {
    (VALID_FROM_JDE..VALID_UNTIL_JDE).contains(&JD)
}

/// One row of the periodic terms for Pluto (Meeus, table 37.A).
///
/// The argument is `j*J + s*S + p*P`. The longitude and latitude
/// coefficients are in units of 1e-6 degree and the radius vector
/// coefficients in units of 1e-7 AU, as printed in the book.
struct Term {
    j: i8,
    s: i8,
    p: i8,
    lon: (i32, i32),
    lat: (i32, i32),
    rad: (i32, i32),
}

const fn t(j: i8, s: i8, p: i8, lon: (i32, i32), lat: (i32, i32), rad: (i32, i32)) -> Term {
    Term {
        j,
        s,
        p,
        lon,
        lat,
        rad,
    }
}

#[rustfmt::skip]
const TERMS: [Term; 43] = [
    t(0,  0,  1, (-19799805, 19850055), (-5452852, -14974862), (66865439, 68951812)),
    t(0,  0,  2, (   897144, -4954829), ( 3527812,   1672790), (-11827535, -332538)),
    t(0,  0,  3, (   611149,  1211027), (-1050748,    327647), ( 1593179, -1438890)),
    t(0,  0,  4, (  -341243,  -189585), (  178690,   -292153), (  -18444,   483220)),
    t(0,  0,  5, (   129287,   -34992), (   18650,    100340), (  -65977,   -85431)),
    t(0,  0,  6, (   -38164,    30893), (  -30697,    -25823), (   31174,    -6032)),
    t(0,  1, -1, (    20442,    -9987), (    4878,     11248), (   -5794,    22161)),
    t(0,  1,  0, (    -4063,    -5071), (     226,       -64), (    4601,     4032)),
    t(0,  1,  1, (    -6016,    -3336), (    2030,      -836), (   -1729,      234)),
    t(0,  1,  2, (    -3956,     3039), (      69,      -604), (    -415,      702)),
    t(0,  1,  3, (     -667,     3572), (    -247,      -567), (     239,      723)),
    t(0,  2, -2, (     1276,      501), (     -57,         1), (      67,      -67)),
    t(0,  2, -1, (     1152,     -917), (    -122,       175), (    1034,     -451)),
    t(0,  2,  0, (      630,    -1277), (     -49,      -164), (    -129,      504)),
    t(1, -1,  0, (     2571,     -459), (    -197,       199), (     480,     -231)),
    t(1, -1,  1, (      899,    -1449), (     -25,       217), (       2,     -441)),
    t(1,  0, -3, (    -1016,     1043), (     589,      -248), (   -3359,      265)),
    t(1,  0, -2, (    -2343,    -1012), (    -269,       711), (    7856,    -7832)),
    t(1,  0, -1, (     7042,      788), (     185,       193), (      36,    45763)),
    t(1,  0,  0, (     1199,     -338), (     315,       807), (    8663,     8547)),
    t(1,  0,  1, (      418,      -67), (    -130,       -43), (    -809,     -769)),
    t(1,  0,  2, (      120,     -274), (       5,         3), (     263,     -144)),
    t(1,  0,  3, (      -60,     -159), (       2,        17), (    -126,       32)),
    t(1,  0,  4, (      -82,      -29), (       2,         5), (     -35,      -16)),
    t(1,  1, -3, (      -36,      -29), (       2,         3), (     -19,       -4)),
    t(1,  1, -2, (      -40,        7), (       3,         1), (     -15,        8)),
    t(1,  1, -1, (      -14,       22), (       2,        -1), (      -4,       12)),
    t(1,  1,  0, (        4,       13), (       1,        -1), (       5,        6)),
    t(1,  1,  1, (        5,        2), (       0,        -1), (       3,        1)),
    t(1,  1,  3, (       -1,        0), (       0,         0), (       6,       -2)),
    t(2,  0, -6, (        2,        0), (       0,        -2), (       2,        2)),
    t(2,  0, -5, (       -4,        5), (       2,         2), (      -2,       -2)),
    t(2,  0, -4, (        4,       -7), (      -7,         0), (      14,       13)),
    t(2,  0, -3, (       14,       24), (      10,        -8), (     -63,       13)),
    t(2,  0, -2, (      -49,      -34), (      -3,        20), (     136,     -236)),
    t(2,  0, -1, (      163,      -48), (       6,         5), (     273,     1065)),
    t(2,  0,  0, (        9,      -24), (      14,        17), (     251,      149)),
    t(2,  0,  1, (       -4,        1), (      -2,         0), (     -25,       -9)),
    t(2,  0,  2, (       -3,        1), (       0,         0), (       9,       -2)),
    t(2,  0,  3, (        1,        3), (       0,         0), (      -8,        7)),
    t(3,  0, -2, (       -3,       -1), (       0,         1), (       2,      -10)),
    t(3,  0, -1, (        5,       -3), (       0,         0), (      19,       35)),
    t(3,  0,  0, (        0,        0), (       1,         0), (      10,        3)),
];

/// Sums of the periodic terms for Julian centuries `T` from J2000.0:
/// `(longitude in degrees, latitude in degrees, radius vector in AU)`.
fn periodic_sums(T: f64) -> (f64, f64, f64) {
    // Mean longitudes of Jupiter, Saturn and Pluto, in degrees.
    let J = 34.35 + 3034.9057 * T;
    let S = 50.08 + 1222.1138 * T;
    let P = 238.96 + 144.96 * T;

    let (mut lon, mut lat, mut rad) = (0.0, 0.0, 0.0);
    for term in TERMS.iter() {
        let alpha =
            (f64::from(term.j) * J + f64::from(term.s) * S + f64::from(term.p) * P).to_radians();
        let (sin, cos) = alpha.sin_cos();
        lon += f64::from(term.lon.0) * sin + f64::from(term.lon.1) * cos;
        lat += f64::from(term.lat.0) * sin + f64::from(term.lat.1) * cos;
        rad += f64::from(term.rad.0) * sin + f64::from(term.rad.1) * cos;
    }

    (lon * 1e-6, lat * 1e-6, rad * 1e-7)
}

/// Computes the heliocentric coordinates of Pluto, referred to the
/// ecliptic and standard equinox of J2000.0
///
/// The position is heliocentric (not barycentric) and geometric (not
/// corrected for light-time). Meeus, chapter 37.
///
/// # Returns
///
/// `(long, lat, rad_vec)`
///
/// * `long`   : Heliocentric longitude of Pluto *| in radians*, in `[0, 2π)`
/// * `lat`    : Heliocentric latitude of Pluto *| in radians*
/// * `rad_vec`: Heliocentric radius vector of Pluto *| in AU*
///
/// Within 1885–2099, and compared with the numerical integration the series
/// was fitted to, the error is below 0.07″ in longitude, 0.02″ in latitude
/// and 0.000006 AU in radius vector.
///
/// **Valid range:** the years 1885 to 2099 only. This function does not
/// check the date: outside that span it silently extrapolates the series and
/// the result is meaningless. Use [`heliocent_pos_checked`] to get `None`
/// instead.
///
/// # Arguments
///
/// * `JD`: Julian Ephemeris Day
pub fn heliocent_pos(JD: f64) -> (f64, f64, f64) {
    let T = time::julian_cent(JD);
    let (sum_lon, sum_lat, sum_rad) = periodic_sums(T);

    let long = 238.958116 + 144.96 * T + sum_lon;
    let lat = -3.908239 + sum_lat;
    let rad_vec = 40.7241346 + sum_rad;

    (norm_2pi(long.to_radians()), lat.to_radians(), rad_vec)
}

/// Computes the heliocentric coordinates of Pluto, referred to the
/// ecliptic and standard equinox of J2000.0, checking the date
///
/// Same as [`heliocent_pos`], but returns `None` when `JD` is outside the
/// valid span 1885–2099 (see [`in_valid_range`]).
///
/// # Arguments
///
/// * `JD`: Julian Ephemeris Day
pub fn heliocent_pos_checked(JD: f64) -> Option<(f64, f64, f64)> {
    if in_valid_range(JD) {
        Some(heliocent_pos(JD))
    } else {
        None
    }
}

/// Mean obliquity of the ecliptic at J2000.0, 23°26'21.448" (IAU 1976),
/// in degrees. Its sine and cosine are 0.397777156 and 0.917482062, the
/// values used in Meeus, chapter 37.
const OBLQ_J2000_DEG: f64 = 23.0 + 26.0 / 60.0 + 21.448 / 3600.0;

/// Converts heliocentric ecliptic coordinates referred to J2000.0 into
/// heliocentric rectangular equatorial coordinates referred to J2000.0
///
/// Meeus, chapter 37, equation 37.1, using the mean obliquity of J2000.0.
///
/// # Returns
///
/// `(x, y, z)` *| in AU*: the x-axis points to the J2000.0 equinox, the
/// z-axis to the J2000.0 north celestial pole.
///
/// # Arguments
///
/// * `long`   : Heliocentric ecliptic longitude (J2000.0) *| in radians*
/// * `lat`    : Heliocentric ecliptic latitude (J2000.0) *| in radians*
/// * `rad_vec`: Heliocentric radius vector *| in AU*
pub fn heliocent_rect_eq_coords(long: f64, lat: f64, rad_vec: f64) -> (f64, f64, f64) {
    let (sin_l, cos_l) = long.sin_cos();
    let (sin_b, cos_b) = lat.sin_cos();
    let (sin_e, cos_e) = OBLQ_J2000_DEG.to_radians().sin_cos();

    let x = rad_vec * cos_l * cos_b;
    let y = rad_vec * (sin_l * cos_b * cos_e - sin_b * sin_e);
    let z = rad_vec * (sin_l * cos_b * sin_e + sin_b * cos_e);

    (x, y, z)
}

/// Geocentric astrometric position of Pluto, referred to J2000.0
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeocentPos {
    /// Astrometric right ascension (J2000.0) *| in radians*, in `[0, 2π)`
    pub asc: f64,
    /// Astrometric declination (J2000.0) *| in radians*
    pub dec: f64,
    /// Distance from the Earth at the given instant to Pluto at the
    /// instant the light left it (light-time corrected) *| in AU*
    pub dist: f64,
    /// True (geometric) Earth-Pluto distance at the given instant *| in AU*
    pub true_dist: f64,
    /// Light-time from Pluto to the Earth *| in days*
    pub light_time: f64,
}

/// Light-time in days for a distance in AU (Meeus, eq. 33.3)
#[inline]
fn light_time(dist: f64) -> f64 {
    0.0057755183 * dist
}

/// Geocentric rectangular equatorial coordinates (J2000.0) of Pluto at
/// `jde`, given the Sun's geocentric coordinates.
fn geocent_rect(jde: f64, sun_x: f64, sun_y: f64, sun_z: f64) -> (f64, f64, f64) {
    let (l, b, r) = heliocent_pos(jde);
    let (x, y, z) = heliocent_rect_eq_coords(l, b, r);
    (sun_x + x, sun_y + y, sun_z + z)
}

/// Computes the geocentric astrometric right ascension and declination of
/// Pluto, referred to the standard equinox of J2000.0
///
/// The heliocentric position from [`heliocent_pos`] is converted to
/// rectangular equatorial coordinates (eq. 37.1) and combined with the
/// Sun's geocentric rectangular coordinates. The effect of light-time is
/// taken into account by iterating on the instant `JD - τ` at which the
/// light left Pluto (Meeus, chapters 33 and 37). The result is
/// *astrometric*: it is not corrected for aberration or nutation, which
/// makes it directly comparable with star catalogue positions for J2000.0.
///
/// # Returns
///
/// `Some(GeocentPos)`, or `None` if `JD` is outside the valid span
/// 1885–2099 of the Pluto series (see [`in_valid_range`]).
///
/// # Arguments
///
/// * `JD`   : Julian Ephemeris Day
/// * `sun_x`, `sun_y`, `sun_z`: Geocentric rectangular equatorial
///   coordinates of the Sun at `JD`, referred to J2000.0 *| in AU*
///   (Meeus, chapter 26). The x-axis points to the J2000.0 equinox and
///   the z-axis to the J2000.0 north celestial pole.
pub fn geocent_astrometric_pos(JD: f64, sun_x: f64, sun_y: f64, sun_z: f64) -> Option<GeocentPos> {
    if !in_valid_range(JD) || !(sun_x.is_finite() && sun_y.is_finite() && sun_z.is_finite()) {
        return None;
    }

    let dist_of = |(x, y, z): (f64, f64, f64)| (x * x + y * y + z * z).sqrt();

    let mut rect = geocent_rect(JD, sun_x, sun_y, sun_z);
    let true_dist = dist_of(rect);
    let mut tau = light_time(true_dist);

    // Pluto moves so slowly that two or three iterations reach full
    // double precision; the limit only guards against pathological input.
    for _ in 0..10 {
        rect = geocent_rect(JD - tau, sun_x, sun_y, sun_z);
        let new_tau = light_time(dist_of(rect));
        let converged = (new_tau - tau).abs() < 1e-12;
        tau = new_tau;
        if converged {
            break;
        }
    }

    let (xi, eta, zeta) = rect;
    Some(GeocentPos {
        asc: norm_2pi(eta.atan2(xi)),
        dec: zeta.atan2(xi.hypot(eta)),
        dist: dist_of(rect),
        true_dist,
        light_time: tau,
    })
}

/// Returns the mean orbital elements of Pluto near 2000 AD
///
/// The angular elements are referred to the ecliptic and equinox of
/// J2000.0 (Meeus, chapter 37). They describe a mean orbit only and are not
/// suitable for an accurate position; use [`heliocent_pos`] for that.
///
/// # Returns
///
/// `(a, e, i, omega, w)`
///
/// * `a`    : Semimajor axis of the orbit *| in AU*
/// * `e`    : Eccentricity of the orbit
/// * `i`    : Inclination of the plane of the orbit with the
///   plane of the ecliptic *| in radians*
/// * `omega`: Longitude of the ascending node *| in radians*
/// * `w`    : Argument of the perihelion *| in radians*
pub fn mn_orb_elements_2000AD() -> (f64, f64, f64, f64, f64) {
    (
        39.543,                   // a
        0.249,                    // e
        17.140_f64.to_radians(),  // i
        110.307_f64.to_radians(), // omega
        113.768_f64.to_radians(), // w
    )
}

/// Reduces an angle in radians to `[0, 2π)`.
#[inline]
fn norm_2pi(angle: f64) -> f64 {
    let r = angle.rem_euclid(TAU);
    if r >= TAU { 0.0 } else { r }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Meeus, example 37.a: the sums of the periodic terms.
    #[test]
    fn periodic_sums_example_37a() {
        let T = time::julian_cent(2448908.5);
        assert!((T - -0.0721834360).abs() < 1e-10);

        let (lon, lat, rad) = periodic_sums(T);
        assert!((lon - 4.246306).abs() < 1e-6, "{lon}");
        assert!((lat - 18.496056).abs() < 1e-6, "{lat}");
        assert!((rad - -11.0130236).abs() < 1e-7, "{rad}");
    }

    #[test]
    fn obliquity_matches_book_constants() {
        let (sin_e, cos_e) = OBLQ_J2000_DEG.to_radians().sin_cos();
        assert!((sin_e - 0.397777156).abs() < 5e-10, "{sin_e}");
        assert!((cos_e - 0.917482062).abs() < 5e-10, "{cos_e}");
    }

    #[test]
    fn norm_2pi_range() {
        assert_eq!(norm_2pi(0.0), 0.0);
        assert_eq!(norm_2pi(-1e-300), 0.0);
        assert!((norm_2pi(-0.5) - (TAU - 0.5)).abs() < 1e-15);
        assert!((norm_2pi(TAU + 0.5) - 0.5).abs() < 1e-15);
    }
}

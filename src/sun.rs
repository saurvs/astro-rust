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

//! The Sun
//!
//! Positions of the Sun and related quantities, following Jean Meeus,
//! *Astronomical Algorithms*, 2nd edition:
//!
//! - chapter 25, solar coordinates: a low accuracy method that assumes a
//!   purely elliptical motion of the Earth (about 0.01°,
//!   [`low_accuracy_position`]) and a high accuracy method based on the
//!   complete VSOP87 theory with the corrections for the FK5 frame, nutation
//!   and aberration ([`apparent_position`]);
//! - chapter 26, geocentric rectangular equatorial coordinates referred to the
//!   mean equinox of the date, of J2000.0, of B1950.0 or of any other epoch;
//! - chapter 28, the equation of time ([`equation_of_time`]);
//! - chapter 29, the ephemeris for physical observations of the Sun
//!   ([`physical_ephemeris`]) and the start of Carrington's synodic rotations.
//!
//! The equinoxes and solstices (chapter 27) are in [`crate::season`].
//!
//! Unless stated otherwise, times are Julian Ephemeris Days (dynamical time,
//! TT); convert a Julian Day in UT with [`time::julian_ephemeris_day`].
//! Angles are in radians.
//!
//! # Example
//!
//! ```
//! use astro::sun;
//!
//! // 1992 October 13.0 TD (Meeus, examples 25.a and 25.b)
//! let pos = sun::apparent_position(2448908.5);
//!
//! // apparent right ascension 13h13m30.75s, declination -7°47'01.7"
//! let ra_hours = pos.right_ascension.to_degrees() / 15.0;
//! let dec_deg = pos.declination.to_degrees();
//! assert!((ra_hours - 13.225208).abs() < 1e-5);
//! assert!((dec_deg + 7.783817).abs() < 1e-5);
//! ```

use crate::angle;
use crate::coords;
use crate::ecliptic;
use crate::nutation;
use crate::planet;
use crate::precess;
use crate::time;
use std;

/// Computes the Sun's equatorial semidiameter
///
/// # Arguments
///
/// * `sun_earth_dist`: Sun-Earth distance *| in AU*
#[inline]
pub fn semidiameter(sun_earth_dist: f64) -> f64 {
    angle::deg_frm_dms(0, 0, 959.63) / sun_earth_dist
}

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

const TAU: f64 = std::f64::consts::TAU;
const PI: f64 = std::f64::consts::PI;

/// Julian Day of the epoch J2000.0 (2000 January 1.5 TT).
const J2000: f64 = 2451545.0;

/// Reduces an angle in radians to `[0, 2π)`.
#[inline]
pub(crate) fn norm_two_pi(a: f64) -> f64 {
    let r = a.rem_euclid(TAU);
    // rem_euclid can round up to exactly 2π for tiny negative inputs
    if r >= TAU { 0.0 } else { r }
}

/// Reduces an angle in radians to `[-π, π)`.
#[inline]
pub(crate) fn norm_pi(a: f64) -> f64 {
    norm_two_pi(a + PI) - PI
}

/// Reduces an angle in degrees to `[0, 360)` and converts it to radians.
#[inline]
fn deg_to_norm_rad(deg: f64) -> f64 {
    norm_two_pi(deg.rem_euclid(360.0).to_radians())
}

/// Arcseconds to radians.
#[inline]
fn arcsec(sec: f64) -> f64 {
    (sec / 3600.0).to_radians()
}

// ---------------------------------------------------------------------------
// Chapter 25, low accuracy
// ---------------------------------------------------------------------------

/// Computes the geometric mean longitude of the Sun, referred to the mean
/// equinox of the date (Meeus, eq. 25.2)
///
/// This is the low-accuracy polynomial in Julian centuries. A slightly more
/// accurate value, from VSOP87, is given by [`mean_longitude_vsop87`].
///
/// # Returns
///
/// * `L0`: Geometric mean longitude *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
pub fn geometric_mean_longitude(jde: f64) -> f64 {
    let t = time::julian_cent(jde);
    deg_to_norm_rad(Horner_eval!(t, 280.46646, 36000.76983, 0.0003032))
}

/// Computes the mean anomaly of the Sun, which equals the mean anomaly of the
/// Earth (Meeus, eq. 25.3)
///
/// # Returns
///
/// * `M`: Mean anomaly *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
pub fn mean_anomaly(jde: f64) -> f64 {
    let t = time::julian_cent(jde);
    deg_to_norm_rad(Horner_eval!(t, 357.52911, 35999.05029, -0.0001537))
}

/// Computes the eccentricity of the Earth's orbit (Meeus, eq. 25.4)
///
/// # Returns
///
/// * `e`: Eccentricity (dimensionless, about 0.0167 at present)
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
pub fn earth_orbit_eccentricity(jde: f64) -> f64 {
    let t = time::julian_cent(jde);
    Horner_eval!(t, 0.016708634, -0.000042037, -0.0000001267)
}

/// Computes the Sun's equation of the center, the difference between its
/// true and its mean anomaly for an unperturbed elliptical orbit (Meeus,
/// chapter 25)
///
/// # Returns
///
/// * `C`: Equation of the center *| in radians* (within ±0.034 rad, about
///   ±1.92°)
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
pub fn equation_of_center(jde: f64) -> f64 {
    let t = time::julian_cent(jde);
    let m = mean_anomaly(jde);

    (Horner_eval!(t, 1.914602, -0.004817, -0.000014) * m.sin()
        + (0.019993 - 0.000101 * t) * (2.0 * m).sin()
        + 0.000289 * (3.0 * m).sin())
    .to_radians()
}

/// Longitude of the ascending node of the Moon's mean orbit, the low accuracy
/// expression used by the corrections of chapter 25 (degrees, unreduced).
#[inline]
fn low_accuracy_omega(t: f64) -> f64 {
    (125.04 - 1934.136 * t).to_radians()
}

/// The Sun's position by the low accuracy method of Meeus, chapter 25
///
/// All longitudes are referred to the mean equinox of the date, except
/// `apparent_longitude` which is referred to the true equinox of the date.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LowAccuracyPosition {
    /// Geometric mean longitude `L0` *| in radians*, in `[0, 2π)`
    pub mean_longitude: f64,
    /// Mean anomaly `M` *| in radians*, in `[0, 2π)`
    pub mean_anomaly: f64,
    /// Eccentricity `e` of the Earth's orbit
    pub eccentricity: f64,
    /// Equation of the center `C` *| in radians*
    pub equation_of_center: f64,
    /// True geometric longitude `☉ = L0 + C` *| in radians*, in `[0, 2π)`
    pub true_longitude: f64,
    /// True anomaly `v = M + C` *| in radians*, in `[0, 2π)`
    pub true_anomaly: f64,
    /// Sun-Earth distance `R` *| in AU* (eq. 25.5)
    pub radius_vector: f64,
    /// Apparent longitude `λ`, corrected for nutation and aberration with
    /// the low accuracy expression of chapter 25 *| in radians*, in `[0, 2π)`
    pub apparent_longitude: f64,
}

/// Computes the Sun's position by the low accuracy method (Meeus,
/// chapter 25)
///
/// The Earth is assumed to move on an unperturbed ellipse; the accuracy is
/// about 0.01°. The Sun's latitude is taken as zero (it never exceeds
/// 1.2″ with respect to the ecliptic of the date).
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
///
/// # Example
///
/// ```
/// use astro::sun;
///
/// // Meeus, example 25.a: 1992 October 13.0 TD
/// let pos = sun::low_accuracy_position(2448908.5);
/// assert!((pos.true_longitude.to_degrees() - 199.90988).abs() < 1e-5);
/// assert!((pos.radius_vector - 0.99766).abs() < 1e-5);
/// ```
pub fn low_accuracy_position(jde: f64) -> LowAccuracyPosition {
    let t = time::julian_cent(jde);

    let mean_longitude = geometric_mean_longitude(jde);
    let mean_anomaly = mean_anomaly(jde);
    let eccentricity = earth_orbit_eccentricity(jde);
    let equation_of_center = equation_of_center(jde);

    let true_longitude = norm_two_pi(mean_longitude + equation_of_center);
    let true_anomaly = norm_two_pi(mean_anomaly + equation_of_center);

    // eq. 25.5
    let radius_vector = 1.000001018 * (1.0 - eccentricity * eccentricity)
        / (1.0 + eccentricity * true_anomaly.cos());

    let omega = low_accuracy_omega(t);
    let apparent_longitude =
        norm_two_pi(true_longitude - (0.00569 + 0.00478 * omega.sin()).to_radians());

    LowAccuracyPosition {
        mean_longitude,
        mean_anomaly,
        eccentricity,
        equation_of_center,
        true_longitude,
        true_anomaly,
        radius_vector,
        apparent_longitude,
    }
}

/// Computes the Sun's apparent right ascension and declination by the low
/// accuracy method (Meeus, eqs. 25.6, 25.7 and 25.8)
///
/// Uses the low accuracy apparent longitude of [`low_accuracy_position`],
/// a zero latitude, and the mean obliquity of the IAU (eq. 22.2) corrected
/// by `+0.00256° cos Ω`. The accuracy is about 0.01°.
///
/// # Returns
///
/// * `eq_point`: Apparent equatorial point of the Sun *| in radians*; right
///   ascension in `[0, 2π)`, declination in `[-π/2, π/2]`
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
pub fn low_accuracy_apparent_eq_point(jde: f64) -> coords::EqPoint {
    let t = time::julian_cent(jde);
    let lambda = low_accuracy_position(jde).apparent_longitude;
    let oblq = ecliptic::mn_oblq_IAU(jde) + (0.00256 * low_accuracy_omega(t).cos()).to_radians();

    coords::EqPoint {
        asc: norm_two_pi((oblq.cos() * lambda.sin()).atan2(lambda.cos())),
        dec: (oblq.sin() * lambda.sin()).asin(),
    }
}

/// Converts a solar longitude referred to the mean equinox of the date to
/// the standard equinox of J2000.0, with the low accuracy expression of
/// Meeus, chapter 25 (`☉2000 = ☉ − 0.01397° (year − 2000)`)
///
/// Sufficient (about 0.01°) between the years 1900 and 2100; for better
/// accuracy use [`geocent_ecl_pos_j2000`].
///
/// # Returns
///
/// * `long_2000`: Longitude referred to the equinox of J2000.0
///   *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `long`: Longitude referred to the mean equinox of the date
///   *| in radians*
/// * `year`: The date, as a year with decimals (e.g. `1992.78`)
pub fn low_accuracy_long_to_j2000(long: f64, year: f64) -> f64 {
    norm_two_pi(long - (0.01397 * (year - 2000.0)).to_radians())
}

// ---------------------------------------------------------------------------
// Chapter 25, higher accuracy
// ---------------------------------------------------------------------------

/// Computes the Sun's geocentric geometric ecliptic position from the
/// complete VSOP87 theory, referred to the mean ecliptic and equinox of
/// the date in the dynamical (VSOP) frame
///
/// The result is *not* converted to the FK5 frame (see
/// [`ecl_coords_to_FK5`] or [`geocent_ecl_pos_fk5`]) and not corrected for
/// nutation or aberration (see [`apparent_position`]). The Earth's
/// heliocentric VSOP87 coordinates are turned into geocentric solar ones
/// by adding 180° to the longitude and changing the sign of the latitude
/// (Meeus, chapter 25).
///
/// # Returns
///
/// `(sun_ecl_point, sun_earth_dist)`
///
/// * `sun_ecl_point` : Ecliptic point of the Sun *| in radians*; longitude
///   in `[0, 2π)`, latitude in `[-π/2, π/2]` (in practice within ±1.2″)
/// * `sun_earth_dist`: Sun-Earth distance *| in AU*
///
/// # Arguments
///
/// * `JD`: Julian Ephemeris Day. VSOP87 is accurate to about 1″ over the
///   years −2000 to +6000.
pub fn geocent_ecl_pos(JD: f64) -> (coords::EclPoint, f64) {
    let (L, B, R) = planet::heliocent_coords(&planet::Planet::Earth, JD);

    // Reduce B to [-π, π) first: planet::heliocent_coords has returned the
    // latitude reduced to [0, 2π), so a small negative B came back as ~2π.
    let ecl_point = coords::EclPoint {
        long: norm_two_pi(L + PI),
        lat: -norm_pi(B),
    };

    (ecl_point, R)
}

/// Converts the Sun's geocentric ecliptic coordinates from the dynamical
/// (VSOP) frame to the FK5 system (Meeus, eq. 25.9)
///
/// The longitude is shifted by −0.09033″ and the latitude by
/// `+0.03916″ (cos λ′ − sin λ′)`, where `λ′ = ☉ − 1.397° T − 0.00031° T²`.
/// These corrections are only significant with the complete VSOP87.
///
/// # Returns
///
/// `(ecl_long_FK5, ecl_lat_FK5)`
///
/// * `ecl_long_FK5`: Ecliptic longitude of the Sun *| in radians*,
///   converted to the FK5 system, in `[0, 2π)`
/// * `ecl_lat_FK5` : Ecliptic latitude of the Sun *| in radians*,
///   converted to the FK5 system
///
/// # Arguments
///
/// * `JD`      : Julian Ephemeris Day
/// * `ecl_long`: Ecliptic longitude of the Sun on `JD`
///   *| in radians*, referred to the mean equinox
///   of the date (VSOP frame)
/// * `ecl_lat` : Ecliptic latitude of the Sun `JD`
///   *| in radians*, referred to the mean equinox
///   of the date (VSOP frame)
pub fn ecl_coords_to_FK5(JD: f64, ecl_long: f64, ecl_lat: f64) -> (f64, f64) {
    let ecl_long_FK5 = norm_two_pi(ecl_long - arcsec(0.09033));

    let JC = time::julian_cent(JD);
    let lambda1 = ecl_long - JC * (1.397 + JC * 0.00031).to_radians();

    let ecl_lat_FK5 = ecl_lat + arcsec(0.03916) * (lambda1.cos() - lambda1.sin());

    (ecl_long_FK5, ecl_lat_FK5)
}

/// Computes the Sun's geocentric geometric ecliptic position from the
/// complete VSOP87 theory, referred to the mean equinox of the date and
/// converted to the FK5 system
///
/// This is [`geocent_ecl_pos`] followed by [`ecl_coords_to_FK5`]: the
/// "true geometric" longitude ☉ and latitude β of Meeus, chapter 25,
/// without nutation and aberration.
///
/// # Returns
///
/// `(sun_ecl_point, sun_earth_dist)`
///
/// * `sun_ecl_point` : Ecliptic point *| in radians*; longitude in `[0, 2π)`
/// * `sun_earth_dist`: Sun-Earth distance *| in AU*
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
pub fn geocent_ecl_pos_fk5(jde: f64) -> (coords::EclPoint, f64) {
    let (p, r) = geocent_ecl_pos(jde);
    let (long, lat) = ecl_coords_to_FK5(jde, p.long, p.lat);

    (coords::EclPoint { long, lat }, r)
}

/// Periodic terms of the daily variation of the Sun's geocentric longitude
/// (Meeus, chapter 25): `(A, B, C, n)` for `A τⁿ sin(B + C τ)` in
/// arcseconds per day, B and C in degrees, τ in Julian millennia.
const DAILY_VARIATION_TERMS: [(f64, f64, f64, i32); 21] = [
    (118.568, 87.5287, 359993.7286, 0),
    (2.476, 85.0561, 719987.4571, 0),
    (1.376, 27.8502, 4452671.1152, 0),
    (0.119, 73.1375, 450368.8564, 0),
    (0.114, 337.2264, 329644.6718, 0),
    (0.086, 222.5400, 659289.3436, 0),
    (0.078, 162.8136, 9224659.7915, 0),
    (0.054, 82.5823, 1079981.1857, 0),
    (0.052, 171.5189, 225184.4282, 0),
    (0.034, 30.3214, 4092677.3866, 0),
    (0.033, 119.8105, 337181.4711, 0),
    (0.023, 247.5418, 299295.6151, 0),
    (0.023, 325.1526, 315559.5560, 0),
    (0.021, 155.1241, 675553.2846, 0),
    (7.311, 333.4515, 359993.7286, 1),
    (0.305, 330.9814, 719987.4571, 1),
    (0.010, 328.5170, 1079981.1857, 1),
    (0.309, 241.4518, 359993.7286, 2),
    (0.021, 205.0482, 719987.4571, 2),
    (0.004, 297.8610, 4452671.1152, 2),
    (0.010, 154.7066, 359993.7286, 3),
];

fn daily_variation_arcsec(jde: f64, constant: f64) -> f64 {
    let tau = time::julian_mill(jde);

    constant
        + DAILY_VARIATION_TERMS
            .iter()
            .map(|&(a, b, c, n)| {
                a * tau.powi(n) * (b + c * tau).rem_euclid(360.0).to_radians().sin()
            })
            .sum::<f64>()
}

/// Computes the daily variation of the Sun's geocentric longitude in a
/// fixed reference frame (the equinox of J2000.0) (Meeus, chapter 25)
///
/// Only the main periodic terms are kept: the error is below 0.1″ per day.
/// This is the quantity needed by [`aberration_in_longitude`].
///
/// # Returns
///
/// * `Δλ`: Daily variation of the longitude *| in radians per day*
///   (about 3548″ per day)
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
pub fn daily_variation_in_longitude(jde: f64) -> f64 {
    arcsec(daily_variation_arcsec(jde, 3548.193))
}

/// Computes the daily variation of the Sun's geocentric longitude referred
/// to the mean equinox of the date (Meeus, chapter 25)
///
/// Same as [`daily_variation_in_longitude`] plus the general precession
/// (the constant term 3548.193″ becomes 3548.330″).
///
/// # Returns
///
/// * `Δλ`: Daily variation of the longitude *| in radians per day*
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
pub fn daily_variation_in_longitude_of_date(jde: f64) -> f64 {
    arcsec(daily_variation_arcsec(jde, 3548.330))
}

/// Computes the aberration in the Sun's longitude from the daily variation
/// of the longitude (Meeus, eq. 25.11)
///
/// The correction is `−0.005775518 R Δλ`, where 0.005775518 day is the
/// light-time for unit distance. Unlike the simpler `−20.4898″ / R`
/// ([`crate::aberr::sol_aberr`], eq. 25.10), it accounts for the perturbed
/// motion of the Earth; the error is below 0.001″.
///
/// # Returns
///
/// * `aberr`: Aberration correction to add to the geometric longitude
///   *| in radians* (about −20.5″)
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
/// * `sun_earth_dist`: Sun-Earth distance *| in AU*
pub fn aberration_in_longitude(jde: f64, sun_earth_dist: f64) -> f64 {
    -0.005775518 * sun_earth_dist * daily_variation_in_longitude(jde)
}

/// The Sun's apparent position and the quantities used to compute it
/// (Meeus, chapter 25, higher accuracy)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ApparentPosition {
    /// True geometric longitude `☉`, FK5, mean equinox of the date
    /// *| in radians*, in `[0, 2π)`
    pub geometric_longitude: f64,
    /// Latitude `β`, FK5 *| in radians*. The apparent and the geometric
    /// latitudes of the Sun are equal to the precision of the theory.
    pub latitude: f64,
    /// Sun-Earth distance `R` *| in AU*
    pub radius_vector: f64,
    /// Nutation in longitude `Δψ` *| in radians*
    pub nutation_in_longitude: f64,
    /// Nutation in obliquity `Δε` *| in radians*
    pub nutation_in_obliquity: f64,
    /// Aberration in longitude *| in radians* (eq. 25.11, about −20.5″)
    pub aberration: f64,
    /// Mean obliquity of the ecliptic `ε0` *| in radians* (Laskar,
    /// eq. 22.3)
    pub mean_obliquity: f64,
    /// True obliquity of the ecliptic `ε = ε0 + Δε` *| in radians*
    pub true_obliquity: f64,
    /// Apparent longitude `λ = ☉ + Δψ + aberration`, true equinox of the
    /// date *| in radians*, in `[0, 2π)`
    pub apparent_longitude: f64,
    /// Apparent right ascension `α` *| in radians*, in `[0, 2π)`
    pub right_ascension: f64,
    /// Apparent declination `δ` *| in radians*, in `[-π/2, π/2]`
    pub declination: f64,
}

/// Computes the Sun's apparent geocentric position with high accuracy
/// (Meeus, chapter 25)
///
/// Steps: the complete VSOP87 theory for the Earth gives the geometric
/// position of the Sun referred to the mean equinox of the date; it is
/// converted to the FK5 system (eq. 25.9), corrected for nutation (IAU 1980
/// theory, chapter 22) and for aberration (eq. 25.11), and converted to
/// right ascension and declination with the true obliquity of the ecliptic
/// (eqs. 13.3 and 13.4).
///
/// The mean obliquity comes from Laskar's formula (eq. 22.3), which agrees
/// with the IAU formula (eq. 22.2) used in the book's examples to better
/// than 0.001″ for 1900–2100 and is more accurate far from J2000.0.
///
/// The accuracy is about 0.01″ near the present, limited by the nutation
/// series; VSOP87 itself is good to about 1″ over −2000 to +6000.
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
///
/// # Example
///
/// ```
/// use astro::sun;
///
/// // Meeus, chapter 25: 1992 October 13.0 TD; the values the book obtains
/// // with the complete VSOP87 theory.
/// let pos = sun::apparent_position(2448908.5);
///
/// let lambda = pos.apparent_longitude.to_degrees();
/// assert!((lambda - (199.0 + 54.0 / 60.0 + 21.56 / 3600.0)).abs() < 0.02 / 3600.0);
/// assert!((pos.latitude.to_degrees() * 3600.0 - 0.72).abs() < 0.01);
/// assert!((pos.radius_vector - 0.99760853).abs() < 2e-8);
/// ```
pub fn apparent_position(jde: f64) -> ApparentPosition {
    let (geo, radius_vector) = geocent_ecl_pos_fk5(jde);
    let (nutation_in_longitude, nutation_in_obliquity) = nutation::nutation(jde);
    let aberration = aberration_in_longitude(jde, radius_vector);

    let mean_obliquity = ecliptic::mn_oblq_laskar(jde);
    let true_obliquity = mean_obliquity + nutation_in_obliquity;

    let apparent_longitude = norm_two_pi(geo.long + nutation_in_longitude + aberration);

    let right_ascension = norm_two_pi(coords::asc_frm_ecl(
        apparent_longitude,
        geo.lat,
        true_obliquity,
    ));
    let declination = coords::dec_frm_ecl(apparent_longitude, geo.lat, true_obliquity);

    ApparentPosition {
        geometric_longitude: geo.long,
        latitude: geo.lat,
        radius_vector,
        nutation_in_longitude,
        nutation_in_obliquity,
        aberration,
        mean_obliquity,
        true_obliquity,
        apparent_longitude,
        right_ascension,
        declination,
    }
}

/// Computes the Sun's apparent right ascension and declination with high
/// accuracy
///
/// Shorthand for the `right_ascension` and `declination` of
/// [`apparent_position`].
///
/// # Returns
///
/// * `eq_point`: Apparent equatorial point *| in radians*; right ascension
///   in `[0, 2π)`
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
pub fn apparent_eq_point(jde: f64) -> coords::EqPoint {
    let pos = apparent_position(jde);

    coords::EqPoint {
        asc: pos.right_ascension,
        dec: pos.declination,
    }
}

// ---------------------------------------------------------------------------
// Chapter 26, rectangular coordinates
// ---------------------------------------------------------------------------

/// Rectangular coordinates *| in AU*
///
/// For the equatorial frames of chapter 26 the origin is the centre of the
/// Earth, the x-axis points to the vernal equinox (0° longitude), the
/// y-axis lies in the plane of the equator towards 90° longitude and the
/// z-axis points to the north celestial pole. For the ecliptic frame the
/// x-y plane is the ecliptic and z points to its north pole.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RectCoords {
    /// X coordinate *| in AU*
    pub x: f64,
    /// Y coordinate *| in AU*
    pub y: f64,
    /// Z coordinate *| in AU*
    pub z: f64,
}

/// Computes the Sun's geocentric rectangular equatorial coordinates from its
/// ecliptic coordinates (Meeus, eq. 26.1)
///
/// # Returns
///
/// `(x, y z)`
///
/// * `x`: The X coordinate *| in AU*
/// * `y`: The Y coordinate *| in AU*
/// * `z`: The Z coordinate *| in AU*
///
/// * The positive x-axis is directed towards the Earth's vernal equinox
///   (0 degrees longitude)
/// * The positive y-axis lies in the plane of the Earth's equator and is
///   directed towards 90 degrees longitude
/// * The positive z-axis is directed towards the Earth's northern
///   celestial pole
///
/// # Arguments
///
/// * `sun_geo_long`: The Sun's geometric longitude *| in radians*,
///   converted to the FK5 system, *without* corrections for nutation and
///   aberration
/// * `sun_geo_lat` : The Sun's geometric latitude *| in radians*,
///   converted to the FK5 system, *without* corrections for nutation and
///   aberration
/// * `sun_rad_vec` : The Sun's geometric radius vector *| in AU*
/// * `mn_oblq`     : Mean obliquity of the ecliptic *| in radians*
pub fn geocent_rect_coords(
    sun_geo_long: f64,
    sun_geo_lat: f64,
    sun_rad_vec: f64,
    mn_oblq: f64,
) -> (f64, f64, f64) {
    let x = sun_rad_vec * sun_geo_lat.cos() * sun_geo_long.cos();

    let y = sun_rad_vec
        * (sun_geo_lat.cos() * sun_geo_long.sin() * mn_oblq.cos()
            - sun_geo_lat.sin() * mn_oblq.sin());

    let z = sun_rad_vec
        * (sun_geo_lat.cos() * sun_geo_long.sin() * mn_oblq.sin()
            + sun_geo_lat.sin() * mn_oblq.cos());

    (x, y, z)
}

/// Computes the Sun's geocentric rectangular equatorial coordinates,
/// referred to the mean equator and equinox of the date (Meeus, chapter 26)
///
/// Uses the geometric FK5 position of [`geocent_ecl_pos_fk5`] (no nutation,
/// no aberration) and the mean obliquity of the date (Laskar, eq. 22.3) in
/// eq. 26.1.
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
pub fn rect_coords_mean_equinox_of_date(jde: f64) -> RectCoords {
    let (p, r) = geocent_ecl_pos_fk5(jde);
    let (x, y, z) = geocent_rect_coords(p.long, p.lat, r, ecliptic::mn_oblq_laskar(jde));

    RectCoords { x, y, z }
}

/// Computes the Sun's geocentric geometric ecliptic position referred to the
/// ecliptic and dynamical equinox of J2000.0 (VSOP87 frame)
///
/// The VSOP87 position of the date ([`geocent_ecl_pos`]) is precessed to
/// J2000.0 with the rigorous ecliptic precession of Meeus, eq. 21.5. This
/// is equivalent to the J2000 series of Meeus, chapter 26 (VSOP87 uses the
/// same IAU 1976 general precession). To obtain the longitude in the FK5
/// frame, subtract 0.09033″.
///
/// # Returns
///
/// `(sun_ecl_point, sun_earth_dist)`
///
/// * `sun_ecl_point` : Ecliptic point *| in radians*; longitude in `[0, 2π)`
/// * `sun_earth_dist`: Sun-Earth distance *| in AU*
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
pub fn geocent_ecl_pos_j2000(jde: f64) -> (coords::EclPoint, f64) {
    let (p, r) = geocent_ecl_pos(jde);
    let (long, lat) = precess::precess_ecl_coords(p.long, p.lat, jde, J2000);

    (
        coords::EclPoint {
            long: norm_two_pi(long),
            lat,
        },
        r,
    )
}

/// Computes the Sun's geocentric rectangular *ecliptic* coordinates
/// referred to the ecliptic and dynamical equinox of J2000.0 (Meeus,
/// eq. 26.2)
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
pub fn ecl_rect_coords_j2000(jde: f64) -> RectCoords {
    let (p, r) = geocent_ecl_pos_j2000(jde);

    RectCoords {
        x: r * p.lat.cos() * p.long.cos(),
        y: r * p.lat.cos() * p.long.sin(),
        z: r * p.lat.sin(),
    }
}

fn mat_mul(m: &[[f64; 3]; 3], c: RectCoords) -> RectCoords {
    RectCoords {
        x: m[0][0] * c.x + m[0][1] * c.y + m[0][2] * c.z,
        y: m[1][0] * c.x + m[1][1] * c.y + m[1][2] * c.z,
        z: m[2][0] * c.x + m[2][1] * c.y + m[2][2] * c.z,
    }
}

/// Rotation from the VSOP87 ecliptic frame of J2000.0 to the FK5 equatorial
/// frame of J2000.0 (Meeus, eq. 26.3)
const VSOP_TO_FK5_J2000: [[f64; 3]; 3] = [
    [1.0, 0.000000440360, -0.000000190919],
    [-0.000000479966, 0.917482137087, -0.397776982902],
    [0.0, 0.397776982902, 0.917482137087],
];

/// Rotation from the VSOP87 ecliptic frame of J2000.0 to the FK5 equatorial
/// frame of B1950.0 (Meeus, chapter 26)
const VSOP_TO_FK5_B1950: [[f64; 3]; 3] = [
    [0.999925702634, 0.012189716217, 0.000011134016],
    [-0.011179418036, 0.917413998946, -0.397777041885],
    [-0.004859003787, 0.397747363646, 0.917482111428],
];

/// Transforms rectangular coordinates from the VSOP87 ecliptic dynamical
/// frame of J2000.0 to the equatorial FK5 frame of J2000.0 (Meeus,
/// eq. 26.3)
///
/// # Arguments
///
/// * `ecl`: Rectangular ecliptic coordinates, VSOP87 frame of J2000.0
pub fn vsop_to_fk5_j2000(ecl: RectCoords) -> RectCoords {
    mat_mul(&VSOP_TO_FK5_J2000, ecl)
}

/// Transforms rectangular coordinates from the VSOP87 ecliptic dynamical
/// frame of J2000.0 to the equatorial FK5 frame of B1950.0 (Meeus,
/// chapter 26)
///
/// The result refers to the mean equator and equinox of B1950.0 in the FK5
/// system, not in the FK4 system (which has the "equinox error").
///
/// # Arguments
///
/// * `ecl`: Rectangular ecliptic coordinates, VSOP87 frame of J2000.0
pub fn vsop_to_fk5_b1950(ecl: RectCoords) -> RectCoords {
    mat_mul(&VSOP_TO_FK5_B1950, ecl)
}

/// Precesses rectangular equatorial FK5 coordinates from the mean equinox
/// of J2000.0 to the mean equinox of another epoch (Meeus, chapter 26,
/// with the angles ζ, z, θ of eq. 21.3)
///
/// # Arguments
///
/// * `eq`        : Rectangular equatorial coordinates, FK5, J2000.0
/// * `jd_equinox`: Julian (Ephemeris) Day of the epoch of the new equinox
pub fn precess_rect_coords_frm_j2000(eq: RectCoords, jd_equinox: f64) -> RectCoords {
    let t = time::julian_cent(jd_equinox);

    // Meeus, eq. 21.3 with the starting epoch at J2000.0 (T = 0)
    let zeta = arcsec(t * (2306.2181 + t * (0.30188 + t * 0.017998)));
    let z = arcsec(t * (2306.2181 + t * (1.09468 + t * 0.018203)));
    let theta = arcsec(t * (2004.3109 - t * (0.42665 + t * 0.041833)));

    let (sz, cz) = zeta.sin_cos();
    let (s_z, c_z) = z.sin_cos();
    let (st, ct) = theta.sin_cos();

    // Rows are (Xx, Yx, Zx), (Xy, Yy, Zy), (Xz, Yz, Zz) in Meeus' notation.
    let m = [
        [
            cz * c_z * ct - sz * s_z,
            -cz * s_z - sz * c_z * ct,
            -c_z * st,
        ],
        [
            sz * c_z + cz * s_z * ct,
            cz * c_z - sz * s_z * ct,
            -s_z * st,
        ],
        [cz * st, -sz * st, ct],
    ];

    mat_mul(&m, eq)
}

/// Computes the Sun's geocentric rectangular equatorial coordinates,
/// referred to the mean equator and equinox of J2000.0 in the FK5 system
/// (Meeus, chapter 26)
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
///
/// # Example
///
/// ```
/// use astro::sun;
///
/// // Meeus, example 26.b: 1992 October 13.0 TD; the book's values from the
/// // complete VSOP87 theory.
/// let c = sun::rect_coords_j2000(2448908.5);
/// assert!((c.x + 0.93739707).abs() < 2e-8);
/// assert!((c.y + 0.31316725).abs() < 2e-8);
/// assert!((c.z + 0.13577842).abs() < 2e-8);
/// ```
pub fn rect_coords_j2000(jde: f64) -> RectCoords {
    vsop_to_fk5_j2000(ecl_rect_coords_j2000(jde))
}

/// Computes the Sun's geocentric rectangular equatorial coordinates,
/// referred to the mean equator and equinox of B1950.0 in the FK5 system
/// (Meeus, chapter 26)
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
pub fn rect_coords_b1950(jde: f64) -> RectCoords {
    vsop_to_fk5_b1950(ecl_rect_coords_j2000(jde))
}

/// Computes the Sun's geocentric rectangular equatorial coordinates,
/// referred to the mean equator and equinox of an arbitrary epoch in the FK5
/// system (Meeus, chapter 26)
///
/// The coordinates are those of the instant `jde`, but referred to the
/// equinox of `jd_equinox` (for example J2044.0 = JD 2467616.0).
///
/// # Arguments
///
/// * `jde`       : Julian Ephemeris Day of the position
/// * `jd_equinox`: Julian (Ephemeris) Day of the epoch of the equinox
pub fn rect_coords_equinox(jde: f64, jd_equinox: f64) -> RectCoords {
    precess_rect_coords_frm_j2000(rect_coords_j2000(jde), jd_equinox)
}

// ---------------------------------------------------------------------------
// Chapter 28, equation of time
// ---------------------------------------------------------------------------

/// Computes the Sun's mean longitude from VSOP87, referred to the mean
/// equinox of the date (Meeus, eq. 28.2)
///
/// A more accurate version of [`geometric_mean_longitude`].
///
/// # Returns
///
/// * `L0`: Mean longitude *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
pub fn mean_longitude_vsop87(jde: f64) -> f64 {
    let tau = time::julian_mill(jde);

    deg_to_norm_rad(Horner_eval!(
        tau,
        280.4664567,
        360007.6982779,
        0.03032028,
        1.0 / 49931.0,
        -1.0 / 15300.0,
        -1.0 / 2000000.0
    ))
}

/// Computes the equation of time with high accuracy (Meeus, eq. 28.1)
///
/// The equation of time is apparent solar time minus mean solar time: the
/// hour angle of the true Sun minus that of the mean Sun. When it is
/// positive the true Sun crosses the meridian before the mean Sun. (Some
/// almanacs and older books use the opposite sign.)
///
/// Uses the apparent right ascension, nutation and true obliquity of
/// [`apparent_position`] in [`crate::planet::earth::equation_of_time`].
///
/// # Returns
///
/// * `E`: Equation of time *| in radians* of hour angle, in `[-π, π)`
///   (in practice within ±0.0873 rad = ±20 minutes of time). Multiply
///   `E.to_degrees()` by 4 to get minutes of time.
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
///
/// # Example
///
/// ```
/// use astro::sun;
///
/// // Meeus, example 28.a: 1992 October 13.0 TD, E = +13m42.6s
/// let minutes = sun::equation_of_time(2448908.5).to_degrees() * 4.0;
/// assert!((minutes - 13.7094).abs() < 0.001);
/// ```
pub fn equation_of_time(jde: f64) -> f64 {
    let pos = apparent_position(jde);

    planet::earth::equation_of_time(
        jde,
        pos.right_ascension,
        pos.nutation_in_longitude,
        pos.true_obliquity,
    )
}

/// Computes the equation of time with Smart's formula (Meeus, eq. 28.3)
///
/// A trigonometric series in the Sun's mean longitude and mean anomaly, the
/// eccentricity of the Earth's orbit and the obliquity of the ecliptic
/// (chapter 25 low accuracy quantities and IAU mean obliquity, eq. 22.2).
/// It neglects the perturbations and the nutation; the error is of the
/// order of a few seconds of time. Same sign convention as
/// [`equation_of_time`].
///
/// # Returns
///
/// * `E`: Equation of time *| in radians* of hour angle
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day
pub fn equation_of_time_smart(jde: f64) -> f64 {
    let oblq = ecliptic::mn_oblq_IAU(jde);
    let l0 = geometric_mean_longitude(jde);
    let e = earth_orbit_eccentricity(jde);
    let m = mean_anomaly(jde);

    let y = (oblq / 2.0).tan().powi(2);

    y * (2.0 * l0).sin() - 2.0 * e * m.sin() + 4.0 * e * y * m.sin() * (2.0 * l0).cos()
        - 0.5 * y * y * (4.0 * l0).sin()
        - 1.25 * e * e * (2.0 * m).sin()
}

// ---------------------------------------------------------------------------
// Chapter 29, ephemeris for physical observations
// ---------------------------------------------------------------------------

/// Carrington's longitude θ of the prime heliographic meridian (rotation
/// with a sidereal period of 25.38 days, fixed by convention), in radians
/// reduced to `[0, 2π)`.
fn carrington_theta(jd: f64) -> f64 {
    deg_to_norm_rad((jd - 2398220.0) * 360.0 / 25.38)
}

/// Longitude K of the ascending node of the solar equator on the ecliptic,
/// in radians.
fn solar_equator_node(jd: f64) -> f64 {
    (73.6667 + 1.3958333 * (jd - 2396758.0) / 36525.0).to_radians()
}

/// Inclination I of the solar equator on the ecliptic (7.25°, by
/// convention constant), in radians.
fn solar_equator_inclination() -> f64 {
    7.25_f64.to_radians()
}

/// Return quantites used in the ephemeris for physical observations of
/// the Sun, from given solar coordinates (Meeus, chapter 29)
///
/// Carrington's elements are used: sidereal rotation period 25.38 days,
/// inclination of the solar equator 7.25°. See [`physical_ephemeris`] to
/// compute everything from the time only.
///
/// # Returns
///
/// `(P, B0, L0)`
///
/// * `P` : Position angle of the northern extremity of the axis of
///   rotation, measured eastwards from the North point of the
///   solar disk *| in radians*, in `(-π/2, π/2)` (positive when the axis
///   is tilted to the east; within ±26.3° in practice)
/// * `B0`: Heliographic latitude of the center of the solar
///   disk *| in radians* (within ±7.25°)
/// * `L0`: Heliographic longitude of the center of the solar
///   disk *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `JD`      : Julian Ephemeris Day
/// * `app_long`: Apparent longitude of the Sun *| in radians*,
///   including the effect of abberation and *not* that
///   of nutation
/// * `app_long_with_nut`: Apparent longitude of the Sun *| in radians*,
///   including the effect of abberation *and* that
///   of nutation
/// * `oblq_eclip`: True obliquity of the ecliptic *| in radians*
pub fn ephemeris(
    JD: f64,
    app_long: f64,
    app_long_with_nut: f64,
    oblq_eclip: f64,
) -> (f64, f64, f64) {
    let theta = carrington_theta(JD);
    let I = solar_equator_inclination();
    let K = solar_equator_node(JD);

    let z = app_long - K;
    let (sin_z, cos_z) = z.sin_cos();

    // atan gives angles in (-π/2, π/2), the range the book asks for
    let x = (-app_long_with_nut.cos() * oblq_eclip.tan()).atan();
    let y = (-cos_z * I.tan()).atan();

    let B_0 = (sin_z * I.sin()).asin();
    let eta = (-sin_z * I.cos()).atan2(-cos_z);
    let L_0 = norm_two_pi(eta - theta);

    let P = x + y;

    (P, B_0, L_0)
}

/// Quantities of the ephemeris for physical observations of the Sun
/// (Meeus, chapter 29)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhysicalEphemeris {
    /// `P`: position angle of the northern extremity of the rotation axis,
    /// measured eastwards from the north point of the disk *| in radians*,
    /// in `(-π/2, π/2)`; positive when the axis is tilted to the east
    pub axis_position_angle: f64,
    /// `B0`: heliographic latitude of the centre of the disk *| in radians*;
    /// positive when the Sun's north pole is tilted towards the Earth
    pub center_helio_lat: f64,
    /// `L0`: heliographic (Carrington) longitude of the centre of the disk
    /// *| in radians*, in `[0, 2π)`
    pub center_helio_long: f64,
}

/// Computes the ephemeris for physical observations of the Sun (Meeus,
/// chapter 29)
///
/// The apparent longitude (with aberration, without nutation) comes from
/// [`apparent_position`] and is passed to [`ephemeris`] together with the
/// nutation in longitude and the true obliquity.
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day. For an instant in UT add ΔT first.
///
/// # Example
///
/// ```
/// use astro::sun;
///
/// // Meeus, example 29.a: 1992 October 13 at 0h UT, ΔT = 59 s
/// let eph = sun::physical_ephemeris(2448908.50068);
/// assert!((eph.axis_position_angle.to_degrees() - 26.27).abs() < 0.005);
/// assert!((eph.center_helio_lat.to_degrees() - 5.99).abs() < 0.005);
/// assert!((eph.center_helio_long.to_degrees() - 238.63).abs() < 0.005);
/// ```
pub fn physical_ephemeris(jde: f64) -> PhysicalEphemeris {
    let pos = apparent_position(jde);
    let app_long = pos.geometric_longitude + pos.aberration;
    let app_long_with_nut = app_long + pos.nutation_in_longitude;

    let (p, b0, l0) = ephemeris(jde, app_long, app_long_with_nut, pos.true_obliquity);

    PhysicalEphemeris {
        axis_position_angle: p,
        center_helio_lat: b0,
        center_helio_long: l0,
    }
}

/// Mean synodic period of the Sun's rotation, in days (Meeus, eq. 29.1).
const SYNODIC_PERIOD: f64 = 27.2752316;

/// Computes an approximate time for the beginning of a solar synodic
/// rotation (Meeus, eqs. 29.1 and 29.2)
///
/// A synodic rotation begins when the heliographic longitude `L0` of the
/// centre of the disk passes through 0. Carrington's rotation No. 1 began
/// on 1853 November 9.
///
/// # Returns
///
/// * `JDE`: Julian *Ephemeris* Day (dynamical time) of the beginning of the
///   rotation. Subtract ΔT to obtain UT. Between the years 1850 and 2100
///   the error is below 0.002 day; see [`carrington_rot_start`] for an
///   exact value.
///
/// # Arguments
///
/// * `C`: Carrington's synodic rotation number
pub fn synodic_rot(C: i64) -> f64 {
    let M = (281.96 + 26.882476 * (C as f64)).to_radians();

    2398140.227 + SYNODIC_PERIOD * (C as f64) + 0.1454 * M.sin()
        - 0.0085 * (2.0 * M).sin()
        - 0.0141 * (2.0 * M).cos()
}

/// Computes the time of the beginning of a solar synodic rotation exactly,
/// by solving `L0 = 0` (Meeus, chapter 29)
///
/// Starts from [`synodic_rot`] and refines it with Newton iterations on the
/// heliographic longitude `L0` of [`physical_ephemeris`] (which decreases by
/// about 13.2° per day).
///
/// # Returns
///
/// * `JDE`: Julian Ephemeris Day of the beginning of the rotation
///   (subtract ΔT to obtain UT)
///
/// # Arguments
///
/// * `C`: Carrington's synodic rotation number
pub fn carrington_rot_start(C: i64) -> f64 {
    let rate = TAU / SYNODIC_PERIOD; // radians per day at which L0 decreases
    let mut jde = synodic_rot(C);

    for _ in 0..20 {
        let l0 = norm_pi(physical_ephemeris(jde).center_helio_long);
        let step = l0 / rate;
        jde += step;
        if step.abs() < 1e-7 {
            break;
        }
    }

    jde
}

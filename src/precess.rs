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

//! Corrections for precession
//!
//! Reduction of mean positions from one equinox to another, following
//! Meeus, *Astronomical Algorithms*, 2nd ed., chapters 21 (precession of
//! equatorial and ecliptical coordinates) and 24 (reduction of the
//! orbital elements i, ω, Ω).
//!
//! Two precession models are available:
//!
//! * the IAU 1976 (Lieske) expressions, for positions referred to the
//!   **FK5** system and the J2000.0 standard epoch: [`precess_eq_coords`],
//!   [`eq_precess_angles`], [`precess_ecl_coords`], [`ecl_precess_angles`],
//!   [`precess_ecl_elements`];
//! * Newcomb's expressions, for old positions referred to the **FK4**
//!   system and Besselian epochs: [`precess_eq_coords_newcomb`],
//!   [`eq_precess_angles_newcomb`]. When turning FK4 mean places into
//!   apparent places, also drop the E-terms of the aberration
//!   ([`aberr::annual_aberr_in_eq_coords_fk4`](crate::aberr::annual_aberr_in_eq_coords_fk4))
//!   and add [`fk4_equinox_correction`] to the final right ascension.
//!
//! These polynomials are only valid within a few thousand years of
//! J2000.0 (Meeus warns that around the year 32700 they put Polaris at a
//! declination of −87°).
//!
//! All the functions here deal with *mean* places: proper motion must be
//! added separately (see [`star::eq_coords_frm_motion`](crate::star::eq_coords_frm_motion)
//! and [`star::apparent_place`](crate::star::apparent_place)).
//!
//! # Example
//!
//! Meeus, example 21.b: θ Persei, J2000.0 → 2028 November 13.19 TD
//! (proper motion already applied).
//!
//! ```
//! use astro::precess;
//!
//! let (asc, dec) = precess::precess_eq_coords(
//!     41.054063_f64.to_radians(),
//!     49.227750_f64.to_radians(),
//!     precess::J2000,
//!     2462088.69,
//! );
//! assert!((asc.to_degrees() - 41.547214).abs() < 1e-6);
//! assert!((dec.to_degrees() - 49.348483).abs() < 1e-6);
//! ```

use crate::time;
use std::f64::consts::{PI, TAU};

/// One arcsecond in radians
const ASEC: f64 = PI / 648000.0;

#[inline]
fn norm_two_pi(a: f64) -> f64 {
    let r = a.rem_euclid(TAU);
    if r >= TAU { 0.0 } else { r }
}

// ---------------------------------------------------------------------------
// Epochs
// ---------------------------------------------------------------------------

/// Julian Ephemeris Day of the standard epoch J2000.0 (2000 January 1.5 TD)
pub const J2000: f64 = 2451545.0;

/// Julian Ephemeris Day of the Besselian epoch B1900.0, as given by Meeus
/// (chapter 21)
pub const B1900: f64 = 2415020.3135;

/// Julian Ephemeris Day of the Besselian epoch B1950.0, as given by Meeus
/// (chapter 21)
pub const B1950: f64 = 2433282.4235;

/// Length of the Julian year, in days
pub const JULIAN_YEAR: f64 = 365.25;

/// Length of the Besselian (tropical) year around 1900, in days
pub const BESSELIAN_YEAR: f64 = 365.242198781;

/// Computes the Julian Ephemeris Day of a Julian epoch
///
/// A Julian epoch such as J2050.0 is J2000.0 plus a whole number of
/// Julian years of 365.25 days.
///
/// # Arguments
///
/// * `epoch`: Julian epoch, as a decimal year (e.g. `2050.0` for J2050.0)
#[inline]
pub fn jd_frm_julian_epoch(epoch: f64) -> f64 {
    J2000 + (epoch - 2000.0) * JULIAN_YEAR
}

/// Computes the Julian epoch (decimal year) of a Julian Ephemeris Day
///
/// Inverse of [`jd_frm_julian_epoch`].
///
/// # Arguments
///
/// * `JD`: Julian Ephemeris Day
#[inline]
pub fn julian_epoch_frm_jd(JD: f64) -> f64 {
    2000.0 + (JD - J2000) / JULIAN_YEAR
}

/// Computes the Julian Ephemeris Day of a Besselian epoch
///
/// Besselian years were used for star catalogues before 1984; each one
/// begins when the Sun's mean longitude (aberration included, from the
/// mean equinox of the date) reaches 280°. This uses B1900.0 = JDE 2415020.31352 and a year of
/// 365.242198781 days (the usual IAU convention), which reproduces Meeus'
/// rounded values of B1900.0 and B1950.0.
///
/// # Arguments
///
/// * `epoch`: Besselian epoch, as a decimal year (e.g. `1950.0`)
#[inline]
pub fn jd_frm_besselian_epoch(epoch: f64) -> f64 {
    2415020.31352 + (epoch - 1900.0) * BESSELIAN_YEAR
}

/// Computes the Besselian epoch (decimal year) of a Julian Ephemeris Day
///
/// Inverse of [`jd_frm_besselian_epoch`].
///
/// # Arguments
///
/// * `JD`: Julian Ephemeris Day
#[inline]
pub fn besselian_epoch_frm_jd(JD: f64) -> f64 {
    1900.0 + (JD - 2415020.31352) / BESSELIAN_YEAR
}

// ---------------------------------------------------------------------------
// Low accuracy
// ---------------------------------------------------------------------------

/// Computes the annual precession in equatorial coordinates (low accuracy)
///
/// Meeus, eq. 21.1. This is the rate of change of the right ascension
/// and declination of a fixed point caused by precession, per Julian
/// year. For a star, add its annual proper motion to these rates before
/// multiplying by the number of years.
///
/// Only suitable when the two epochs are not far apart (a few decades)
/// and the point is not close to a celestial pole (the right ascension
/// rate contains `tan(dec)`).
///
/// # Returns
///
/// `(ann_precess_asc, ann_precess_dec)`
///
/// * `ann_precess_asc`: Annual precession in right ascension *| in radians
///   per year*
/// * `ann_precess_dec`: Annual precession in declination *| in radians per
///   year*
///
/// # Arguments
///
/// * `asc`: Right ascension *| in radians*
/// * `dec`: Declination *| in radians*
/// * `JD` : Julian (Ephemeris) day at which the rates `m` and `n` are
///   evaluated; they change very slowly, so any epoch between the two
///   being reduced will do
pub fn annual_precess(asc: f64, dec: f64, JD: f64) -> (f64, f64) {
    let JC = time::julian_cent(JD);

    // m and n are in seconds of time; 1 s = 15″
    let m = (3.07496 + 0.00186 * JC) * 15.0 * ASEC;
    let n = (1.33621 - 0.00057 * JC) * 15.0 * ASEC;

    (m + n * asc.sin() * dec.tan(), n * asc.cos())
}

/// Computes equatorial coordinates reduced to a different epoch, with the
/// low accuracy method
///
/// Applies the annual precession of [`annual_precess`] (evaluated at the
/// mean of the two epochs and at the initial position) linearly over the
/// interval. Proper motion is not included. See [`annual_precess`] for
/// the limitations; use [`precess_eq_coords`] otherwise.
///
/// # Returns
///
/// `(new_asc, new_dec)` *| in radians*, `new_asc` in `[0, 2π)`
///
/// # Arguments
///
/// * `asc`: Right ascension at the old epoch *| in radians*
/// * `dec`: Declination at the old epoch *| in radians*
/// * `JD1`: Julian (Ephemeris) day of the old epoch
/// * `JD2`: Julian (Ephemeris) day of the new epoch
pub fn precess_eq_coords_low_accuracy(asc: f64, dec: f64, JD1: f64, JD2: f64) -> (f64, f64) {
    let years = (JD2 - JD1) / JULIAN_YEAR;
    let (d_asc, d_dec) = annual_precess(asc, dec, 0.5 * (JD1 + JD2));

    (norm_two_pi(asc + d_asc * years), dec + d_dec * years)
}

// ---------------------------------------------------------------------------
// Rigorous method, equatorial coordinates
// ---------------------------------------------------------------------------

/// The three angles ζ, z and θ that rotate equatorial coordinates from one
/// mean equinox to another
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EqPrecessAngles {
    /// ζ *| in radians*
    pub zeta: f64,
    /// z *| in radians*
    pub z: f64,
    /// θ *| in radians*
    pub theta: f64,
}

/// Computes the IAU 1976 precession angles ζ, z, θ between two epochs
///
/// Meeus, eq. 21.2 (Lieske's expressions, FK5 system). For a starting
/// epoch of J2000.0 these reduce to eq. 21.3.
///
/// # Arguments
///
/// * `JD1`: Julian (Ephemeris) day of the starting epoch
/// * `JD2`: Julian (Ephemeris) day of the final epoch
pub fn eq_precess_angles(JD1: f64, JD2: f64) -> EqPrecessAngles {
    let T = time::julian_cent(JD1);
    let t = (JD2 - JD1) / 36525.0;

    let c = 2306.2181 + T * (1.39656 - 0.000139 * T);

    let zeta = t * (c + t * ((0.30188 - 0.000344 * T) + 0.017998 * t));
    let z = t * (c + t * ((1.09468 + 0.000066 * T) + 0.018203 * t));
    let theta = t
        * ((2004.3109 - T * (0.85330 + 0.000217 * T))
            - t * ((0.42665 + 0.000217 * T) + 0.041833 * t));

    EqPrecessAngles {
        zeta: zeta * ASEC,
        z: z * ASEC,
        theta: theta * ASEC,
    }
}

/// Computes Newcomb's precession angles ζ, z, θ between two epochs
///
/// Meeus, chapter 21, "The old precessional elements". These are the
/// expressions on which the FK4 system and the B1950.0 catalogues are
/// based; time is counted in tropical centuries from B1900.0.
///
/// # Arguments
///
/// * `JD1`: Julian (Ephemeris) day of the starting epoch
/// * `JD2`: Julian (Ephemeris) day of the final epoch
pub fn eq_precess_angles_newcomb(JD1: f64, JD2: f64) -> EqPrecessAngles {
    let T = (JD1 - 2415020.3135) / 36524.2199;
    let t = (JD2 - JD1) / 36524.2199;

    let zeta = t * ((2304.250 + 1.396 * T) + t * (0.302 + 0.018 * t));
    let z = zeta + t * t * (0.791 + 0.001 * t);
    let theta = t * ((2004.682 - 0.853 * T) - t * (0.426 + 0.042 * t));

    EqPrecessAngles {
        zeta: zeta * ASEC,
        z: z * ASEC,
        theta: theta * ASEC,
    }
}

/// Applies the rotation ζ, z, θ (Meeus, eq. 21.4). The declination is
/// taken from atan2 so that it stays accurate near the poles.
fn rotate_eq_coords(asc: f64, dec: f64, angles: &EqPrecessAngles) -> (f64, f64) {
    let (sin_dec, cos_dec) = dec.sin_cos();
    let (sin_h, cos_h) = (asc + angles.zeta).sin_cos();
    let (sin_th, cos_th) = angles.theta.sin_cos();

    let A = cos_dec * sin_h;
    let B = cos_th * cos_dec * cos_h - sin_th * sin_dec;
    let C = sin_th * cos_dec * cos_h + cos_th * sin_dec;

    (norm_two_pi(A.atan2(B) + angles.z), C.atan2(A.hypot(B)))
}

/// Computes equatorial coordinates reduced to a different epoch (FK5)
///
/// Rigorous method of Meeus, eqs. 21.2 and 21.4, with the IAU 1976
/// precession. Use it for positions referred to the FK5 system (or
/// catalogues based on it, J2000.0). For FK4 / B1950.0 positions use
/// [`precess_eq_coords_newcomb`].
///
/// Proper motion is not included: apply it before precessing.
///
/// # Returns
///
/// `(new_asc, new_dec)`
///
/// * `new_asc`: Right ascension in the new epoch *| in radians*, in `[0, 2π)`
/// * `new_dec`: Declination in the new epoch *| in radians*
///
/// # Arguments
///
/// * `old_asc`: Right ascension in the old epoch *| in radians*,
///   referred to the FK5 system
/// * `old_dec`: Declination in the old epoch *| in radians*,
///   referred to the FK5 system
/// * `JD1`    : Julian (Ephemeris) day corresponding to the old epoch
/// * `JD2`    : Julian (Ephemeris) day corresponding to the new epoch
///
/// Valid for a few thousand years around J2000.0.
pub fn precess_eq_coords(old_asc: f64, old_dec: f64, JD1: f64, JD2: f64) -> (f64, f64) {
    rotate_eq_coords(old_asc, old_dec, &eq_precess_angles(JD1, JD2))
}

/// Computes equatorial coordinates referred to the FK4 system reduced to
/// a different epoch, with Newcomb's precession
///
/// Meeus, chapter 21, "The old precessional elements", with the rotation
/// of eq. 21.4. Use it for old (B1950.0, FK4) catalogue positions. The
/// result is still in the FK4 system: to compare it with FK5 positions,
/// add [`fk4_equinox_correction`] to the right ascension (and do not use
/// the E-terms of aberration when computing apparent places).
///
/// # Returns
///
/// `(new_asc, new_dec)` *| in radians*, `new_asc` in `[0, 2π)`
///
/// # Arguments
///
/// * `old_asc`: Right ascension in the old epoch *| in radians*, FK4
/// * `old_dec`: Declination in the old epoch *| in radians*, FK4
/// * `JD1`    : Julian (Ephemeris) day of the old epoch
/// * `JD2`    : Julian (Ephemeris) day of the new epoch
pub fn precess_eq_coords_newcomb(old_asc: f64, old_dec: f64, JD1: f64, JD2: f64) -> (f64, f64) {
    rotate_eq_coords(old_asc, old_dec, &eq_precess_angles_newcomb(JD1, JD2))
}

/// Computes equatorial coordinates referred to the FK4 system reduced to
/// a different epoch
///
/// Despite its name, this function uses Newcomb's (FK4) precession, not
/// the FK5 one. It is the same as [`precess_eq_coords_newcomb`]. For FK5
/// positions use [`precess_eq_coords`].
///
/// Up to astro 2.0.0 this function used the precession angles in
/// degrees as if they were radians, giving wrong results.
///
/// # Returns
///
/// `(new_asc, new_dec)` *| in radians*
///
/// # Arguments
///
/// * `old_asc`: Right ascension in the old epoch *| in radians*, FK4
/// * `old_dec`: Declination in the old epoch *| in radians*, FK4
/// * `JD1`    : Julian (Ephemeris) day corresponding to the old epoch
/// * `JD2`    : Julian (Ephemeris) day corresponding to the new epoch
#[deprecated(note = "this uses Newcomb's FK4 precession despite its name; use \
            `precess_eq_coords_newcomb` (FK4) or `precess_eq_coords` (FK5)")]
pub fn precess_eq_coords_FK5(old_asc: f64, old_dec: f64, JD1: f64, JD2: f64) -> (f64, f64) {
    precess_eq_coords_newcomb(old_asc, old_dec, JD1, JD2)
}

/// Computes the equinox correction of the FK4 system
///
/// The zero point of right ascension in the FK4 is slightly off. Add this
/// correction to a right ascension computed from FK4 data (with
/// [`precess_eq_coords_newcomb`]) to put it on the FK5 equinox.
/// Meeus, chapter 21: 0.0775 s + 0.0850 s per Julian century from J2000.0.
///
/// # Returns
///
/// * `corr`: Correction to add to the right ascension *| in radians*
///
/// # Arguments
///
/// * `JD`: Julian (Ephemeris) day
pub fn fk4_equinox_correction(JD: f64) -> f64 {
    (0.0775 + 0.0850 * time::julian_cent(JD)) * 15.0 * ASEC
}

// ---------------------------------------------------------------------------
// Ecliptical coordinates
// ---------------------------------------------------------------------------

/// The three angles η, Π and p used to reduce ecliptical coordinates and
/// orbital elements from one mean equinox to another
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EclPrecessAngles {
    /// η: inclination of the final ecliptic on the starting one
    /// *| in radians*
    pub eta: f64,
    /// Π: longitude, on the starting ecliptic, of the ascending node of the
    /// final ecliptic *| in radians*
    pub pi: f64,
    /// p: general precession in longitude *| in radians*
    pub p: f64,
}

/// Computes the angles η, Π and p for the precession of ecliptical
/// coordinates between two epochs
///
/// Meeus, eq. 21.5 (eq. 21.6 when the starting epoch is J2000.0).
///
/// # Arguments
///
/// * `JD1`: Julian (Ephemeris) day of the starting epoch
/// * `JD2`: Julian (Ephemeris) day of the final epoch
pub fn ecl_precess_angles(JD1: f64, JD2: f64) -> EclPrecessAngles {
    let T = time::julian_cent(JD1);
    let t = (JD2 - JD1) / 36525.0;

    let eta = t
        * ((47.0029 - T * (0.06603 - 0.000598 * T))
            + t * ((-0.03302 + 0.000598 * T) + 0.000060 * t));

    let pi_asec = T * (3289.4789 + 0.60622 * T) - t * ((869.8089 + 0.50491 * T) - 0.03536 * t);

    let p = t
        * ((5029.0966 + T * (2.22226 - 0.000042 * T))
            + t * ((1.11113 - 0.000042 * T) - 0.000006 * t));

    EclPrecessAngles {
        eta: eta * ASEC,
        pi: 174.876384_f64.to_radians() + pi_asec * ASEC,
        p: p * ASEC,
    }
}

/// Computes ecliptic coordinates reduced to a different epoch
///
/// Rigorous method of Meeus, eqs. 21.5 and 21.7 (IAU 1976 precession).
/// The latitude is taken from atan2 so it stays accurate near the poles
/// of the ecliptic.
///
/// # Returns
///
/// `(new_long, new_lat)`
///
/// * `new_long`: Ecliptic longitude in the new epoch *| in radians*,
///   in `[0, 2π)`
/// * `new_lat` : Ecliptic latitude in the new epoch *| in radians*
///
/// # Arguments
///
/// * `old_long`: Ecliptic longitude in the old epoch *| in radians*
/// * `old_lat` : Ecliptic latitude in the old epoch *| in radians*
/// * `JD_old`  : Julian (Ephemeris) day corresponding to the old epoch
/// * `JD_new`  : Julian (Ephemeris) day corresponding to the new epoch
pub fn precess_ecl_coords(old_long: f64, old_lat: f64, JD_old: f64, JD_new: f64) -> (f64, f64) {
    let EclPrecessAngles { eta, pi, p } = ecl_precess_angles(JD_old, JD_new);

    let (sin_eta, cos_eta) = eta.sin_cos();
    let (sin_lat, cos_lat) = old_lat.sin_cos();
    let (sin_d, cos_d) = (pi - old_long).sin_cos();

    let A = cos_eta * cos_lat * sin_d - sin_eta * sin_lat;
    let B = cos_lat * cos_d;
    let C = cos_eta * sin_lat + sin_eta * cos_lat * sin_d;

    (norm_two_pi(p + pi - A.atan2(B)), C.atan2(A.hypot(B)))
}

/// Computes the proper motion of a star in ecliptical coordinates from its
/// proper motion in equatorial coordinates
///
/// Meeus, chapter 21 (formulae following example 21.c). Both rates must
/// be in the same unit of time; the result uses that unit too.
///
/// # Returns
///
/// `(pmotion_long, pmotion_lat)`
///
/// * `pmotion_long`: Rate of change of the ecliptic longitude, dλ/dt
///   *| in radians per unit time* (not multiplied by `cos(lat)`)
/// * `pmotion_lat` : Rate of change of the ecliptic latitude, dβ/dt
///   *| in radians per unit time*
///
/// # Arguments
///
/// * `asc`        : Right ascension of the star *| in radians*
/// * `dec`        : Declination of the star *| in radians*
/// * `pmotion_asc`: Proper motion in right ascension, dα/dt *| in radians
///   per unit time* (an angle of right ascension, not multiplied by
///   `cos(dec)`; 1 second of time per year is 15″ per year)
/// * `pmotion_dec`: Proper motion in declination, dδ/dt *| in radians per
///   unit time*
/// * `oblq_eclip` : Obliquity of the ecliptic *| in radians*
///
/// The longitude rate is undefined at the poles of the ecliptic, where
/// the result is not finite.
pub fn proper_motion_in_ecl_coords(
    asc: f64,
    dec: f64,
    pmotion_asc: f64,
    pmotion_dec: f64,
    oblq_eclip: f64,
) -> (f64, f64) {
    let (sin_e, cos_e) = oblq_eclip.sin_cos();
    let (sin_a, cos_a) = asc.sin_cos();
    let (sin_d, cos_d) = dec.sin_cos();

    let sin_lat = sin_d * cos_e - cos_d * sin_e * sin_a;
    let cos_lat = (1.0 - sin_lat * sin_lat).max(0.0).sqrt();

    let k = cos_e * cos_d + sin_e * sin_d * sin_a;

    let pmotion_long =
        (pmotion_dec * sin_e * cos_a + pmotion_asc * cos_d * k) / (cos_lat * cos_lat);
    let pmotion_lat = (pmotion_dec * k - pmotion_asc * sin_e * cos_a * cos_d) / cos_lat;

    (pmotion_long, pmotion_lat)
}

// ---------------------------------------------------------------------------
// Orbital elements (chapter 24)
// ---------------------------------------------------------------------------

/// The three orbital elements that depend on the reference ecliptic and
/// equinox
///
/// The semimajor axis, eccentricity, time of perihelion and so on do not
/// change when an orbit is referred to another equinox.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EclElements {
    /// Inclination *| in radians*, in `[0, π]` (above `π/2` for a
    /// retrograde orbit)
    pub inc: f64,
    /// Argument of perihelion *| in radians*
    pub arg_perih: f64,
    /// Longitude of the ascending node *| in radians*
    pub long_ascend_node: f64,
}

/// Rotates i, ω, Ω to a new ecliptic (Meeus, eqs. 24.1 to 24.3).
///
/// `pi` is the longitude, on the old ecliptic, of the ascending node of
/// the new ecliptic, `psi` the same point measured on the new ecliptic,
/// and `eta` (given by its sine and cosine) the angle between the
/// ecliptics.
fn rotate_ecl_elements(
    el: &EclElements,
    sin_eta: f64,
    cos_eta: f64,
    pi: f64,
    psi: f64,
) -> EclElements {
    let (sin_i0, cos_i0) = el.inc.sin_cos();
    let (sin_w, cos_w) = (el.long_ascend_node - pi).sin_cos();

    // sin i sin(Ω - ψ), sin i cos(Ω - ψ) and cos i
    let A = sin_i0 * sin_w;
    let B = cos_eta * sin_i0 * cos_w - sin_eta * cos_i0;
    let cos_i = cos_i0 * cos_eta + sin_i0 * sin_eta * cos_w;
    let sin_i = A.hypot(B);

    if sin_i == 0.0 {
        // The orbit lies in the new ecliptic: the node is undefined. Keep
        // the directions of the node and of the perihelion fixed in space;
        // only the origin of longitudes moves.
        return EclElements {
            inc: if cos_i < 0.0 { PI } else { 0.0 },
            arg_perih: norm_two_pi(el.arg_perih),
            long_ascend_node: norm_two_pi(el.long_ascend_node + psi - pi),
        };
    }

    // sin i sin Δω and sin i cos Δω
    let dw_sin = -sin_eta * sin_w;
    let dw_cos = sin_i0 * cos_eta - cos_i0 * sin_eta * cos_w;

    EclElements {
        inc: sin_i.atan2(cos_i),
        arg_perih: norm_two_pi(el.arg_perih + dw_sin.atan2(dw_cos)),
        long_ascend_node: norm_two_pi(psi + A.atan2(B)),
    }
}

/// Computes orbital elements reduced to a different equinox
///
/// Meeus, chapter 24, eqs. 24.1 to 24.3, with the IAU 1976 precession
/// angles η, Π, p of [`ecl_precess_angles`]. The orbit itself is not
/// changed: the elements remain valid for their original epoch, only the
/// reference ecliptic and equinox change.
///
/// Quadrants are resolved with atan2 throughout, so retrograde orbits
/// (inclination above 90°) and orbits of zero inclination are handled;
/// for a zero initial inclination the new node lies at ψ + 180° and the
/// longitude of perihelion is preserved.
///
/// # Arguments
///
/// * `elements`: Inclination, argument of perihelion and longitude of the
///   ascending node referred to the old equinox *| in radians*
/// * `JD1`     : Julian (Ephemeris) day of the old equinox
/// * `JD2`     : Julian (Ephemeris) day of the new equinox
///
/// Angles of the result are in `[0, 2π)` (inclination in `[0, π]`).
pub fn precess_ecl_elements(elements: &EclElements, JD1: f64, JD2: f64) -> EclElements {
    let EclPrecessAngles { eta, pi, p } = ecl_precess_angles(JD1, JD2);
    let (sin_eta, cos_eta) = eta.sin_cos();

    rotate_ecl_elements(elements, sin_eta, cos_eta, pi, pi + p)
}

/// Computes orbital elements reduced to a different equinox
///
/// Tuple version of [`precess_ecl_elements`]; see there for details.
///
/// Up to astro 2.0.0 this function took the new inclination from its
/// sine only, which is wrong for retrograde orbits, and the change in the
/// argument of perihelion from its sine only.
///
/// # Returns
///
/// `(new_inc, new_arg_perih, new_long_ascend_node)`
///
/// * `new_inc`             : Inclination in the new
///   equinox *| in radians*, in `[0, π]`
/// * `new_arg_perih`       : Argument of perihelion in the
///   new equinox *| in radians*, in `[0, 2π)`
/// * `new_long_ascend_node`: Longitude of ascending node in
///   the new equinox *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `old_inc`             : Inclination in the old equinox *| in radians*
/// * `old_arg_perih`       : Argument of perihelion in the old
///   equinox *| in radians*
/// * `old_long_ascend_node`: Longitude of ascending node in the old
///   equinox *| in radians*
/// * `JD1`                 : Julian (Ephemeris) day corresponding to the old
///   equinox
/// * `JD2`                 : Julian (Ephemeris) day corresponding to the new
///   equinox
pub fn precess_orb_elements(
    old_inc: f64,
    old_arg_perih: f64,
    old_long_ascend_node: f64,
    JD1: f64,
    JD2: f64,
) -> (f64, f64, f64) {
    let el = precess_ecl_elements(
        &EclElements {
            inc: old_inc,
            arg_perih: old_arg_perih,
            long_ascend_node: old_long_ascend_node,
        },
        JD1,
        JD2,
    );

    (el.inc, el.arg_perih, el.long_ascend_node)
}

/// Computes orbital elements referred to the equinox B1950.0 reduced to
/// the equinox J2000.0
///
/// Meeus, eq. 24.4: the general method of [`precess_ecl_elements`] with
/// the constants for B1950.0 → J2000.0 already evaluated. The elements
/// must be referred to the FK5 system; for elements in the old FK4
/// system use [`ecl_elements_fk4_b1950_to_fk5_j2000`].
///
/// # Arguments
///
/// * `elements`: i, ω, Ω referred to the mean ecliptic and equinox of
///   B1950.0 *| in radians*
pub fn ecl_elements_b1950_to_j2000(elements: &EclElements) -> EclElements {
    const SIN_ETA: f64 = 0.0001139788;
    const COS_ETA: f64 = 0.9999999935;

    rotate_ecl_elements(
        elements,
        SIN_ETA,
        COS_ETA,
        174.298782_f64.to_radians(),
        174.997194_f64.to_radians(),
    )
}

/// Computes orbital elements referred to the FK4 system and the equinox
/// B1950.0, converted to the FK5 system and the equinox J2000.0
///
/// Algorithm due to D. K. Yeomans, given in Meeus, chapter 24. It
/// accounts for both the precession and the change of reference system.
///
/// # Arguments
///
/// * `elements`: i, ω, Ω referred to the FK4 system, B1950.0
///   *| in radians*
pub fn ecl_elements_fk4_b1950_to_fk5_j2000(elements: &EclElements) -> EclElements {
    let L1 = 4.50001688_f64.to_radians(); // L'
    let L = 5.19856209_f64.to_radians();
    let J = 0.00651966_f64.to_radians();

    let (sin_i0, cos_i0) = elements.inc.sin_cos();
    let (sin_j, cos_j) = J.sin_cos();
    let (sin_w, cos_w) = (L + elements.long_ascend_node).sin_cos();

    // sin i sin(ω - ω0), sin i cos(ω - ω0)
    let dw_sin = sin_j * sin_w;
    let dw_cos = sin_i0 * cos_j + cos_i0 * sin_j * cos_w;

    let cos_i = cos_i0 * cos_j - sin_i0 * sin_j * cos_w;

    // sin i sin(L' + Ω), sin i cos(L' + Ω)
    let n_sin = sin_i0 * sin_w;
    let n_cos = cos_i0 * sin_j + sin_i0 * cos_j * cos_w;

    let sin_i = n_sin.hypot(n_cos);

    EclElements {
        inc: sin_i.atan2(cos_i),
        arg_perih: norm_two_pi(elements.arg_perih + dw_sin.atan2(dw_cos)),
        long_ascend_node: norm_two_pi(n_sin.atan2(n_cos) - L1),
    }
}

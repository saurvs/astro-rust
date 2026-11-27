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

//! Jupiter: ephemeris for physical observations
//!
//! Implements Meeus, *Astronomical Algorithms*, 2nd ed., chapter 43. For a
//! given instant, [`physical_ephemeris`] gives
//!
//! * the jovicentric declinations of the Earth (`D_E`) and of the Sun
//!   (`D_S`), positive when Jupiter's north pole is tilted towards them,
//! * the position angle `P` of Jupiter's rotation axis, measured from the
//!   north point of the disk towards the east,
//! * the longitudes `ω1` and `ω2` of the central meridian in rotation
//!   Systems I (the equatorial zone, 877.90°/day) and II (higher latitudes,
//!   870.27°/day), for the geometric disk and for the illuminated disk
//!   (the two differ by the *correction for phase*, at most 0.61°).
//!
//! [`physical_ephemeris_low_accuracy`] is the book's shorter method, good to
//! 0.1–0.2° in `ω1` and `ω2`, which needs no planetary theory. The older
//! [`ephemeris`] lets the caller supply the obliquity and nutation.
//!
//! Jovigraphic longitudes increase towards the west, so `ω1` and `ω2`
//! increase with time. All angles are in radians; longitudes and position
//! angles are reduced to `[0, 2π)`.
//!
//! # Example
//!
//! Meeus' example 43.a, 1992 December 16 at 0h UT (JDE 2448972.50068):
//!
//! ```
//! use astro::planet::jupiter;
//!
//! let eph = jupiter::physical_ephemeris(2448972.50068);
//!
//! assert!((eph.earth_dec.to_degrees() - -2.48).abs() < 0.005);
//! assert!((eph.sun_dec.to_degrees() - -2.20).abs() < 0.005);
//! assert!((eph.axis_pos_angle.to_degrees() - 24.80).abs() < 0.005);
//! // central meridians of the geometric disk, Systems I and II
//! assert!((eph.central_meridian_i.to_degrees() - 267.63).abs() < 0.005);
//! assert!((eph.central_meridian_ii.to_degrees() - 72.31).abs() < 0.005);
//! // and of the illuminated disk
//! assert!((eph.illum_central_meridian_i.to_degrees() - 268.06).abs() < 0.005);
//! assert!((eph.illum_central_meridian_ii.to_degrees() - 72.74).abs() < 0.005);
//! ```

pub mod moon;

use crate::angle;
use crate::coords;
use crate::nutation;
use crate::planet;

/// Computes Jupiter's apparent equatorial semidiameter
///
/// `98.44″ / Δ` (Meeus, *Astronomical Algorithms*, 2nd ed., chapter 55).
///
/// # Returns
///
/// * `eq_semidia`: Equatorial semidiameter *| in radians*
///
/// # Arguments
///
/// * `jup_earth_dist`: Jupiter-Earth distance *| in AU*
#[inline(always)]
pub fn eq_semidiameter(jup_earth_dist: f64) -> f64 {
    angle::deg_frm_dms(0, 0, 98.44).to_radians() / jup_earth_dist
}

/// Computes Jupiter's polar semidiameter
///
/// `92.06″ / Δ`, the true polar semidiameter at distance Δ (Meeus,
/// chapter 55). Because the Jovicentric latitude of the Earth never
/// exceeds about 4°, this is also the apparent polar semidiameter to good
/// accuracy; [`apprnt_pol_semidiameter`] gives the rigorous value.
///
/// # Returns
///
/// * `pol_semidia`: Polar semidiameter *| in radians*
///
/// # Arguments
///
/// * `jup_earth_dist`: Jupiter-Earth distance *| in AU*
#[inline(always)]
pub fn pol_semidiameter(jup_earth_dist: f64) -> f64 {
    angle::deg_frm_dms(0, 0, 92.06).to_radians() / jup_earth_dist
}

/// Computes Jupiter's apparent polar semidiameter, allowing for the tilt
/// of its axis towards the Earth
///
/// `s_p = s_e sqrt(1 - k cos²D_E)` with `k = 1 - (b/a)²`, the same
/// procedure Meeus gives for Saturn (chapter 55). It is between
/// [`pol_semidiameter`] (axis in the sky plane) and [`eq_semidiameter`]
/// (pole facing the Earth).
///
/// # Returns
///
/// * Apparent polar semidiameter *| in radians*
///
/// # Arguments
///
/// * `jup_earth_dist`: Jupiter-Earth distance *| in AU*
/// * `earth_lat`     : Jovicentric declination of the Earth `D_E`
///   *| in radians* (see `Ephemeris::De`)
pub fn apprnt_pol_semidiameter(jup_earth_dist: f64, earth_lat: f64) -> f64 {
    let a = eq_semidiameter(1.0);
    let b = pol_semidiameter(1.0);
    let k = 1.0 - (b / a).powi(2);

    (a / jup_earth_dist) * (1.0 - k * earth_lat.cos().powi(2)).sqrt()
}

/// The constant of annual aberration (about 20.49″), in degrees, as used
/// by Meeus for the aberration of the planet in chapters 42, 43 and 45.
const ABERRATION: f64 = 0.005_693;

/// Holds Jupiter's ephemeris values for physical observations
///
/// Returned by [`ephemeris`]; [`PhysicalEphemeris`] holds the same values
/// and a few more. All angles are in radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ephemeris {
    /// Jupiter-centric declination of the Earth, `[-π/2, π/2]`
    pub De: f64,
    /// Jupiter-centric declination of the Sun, `[-π/2, π/2]`
    pub Ds: f64,
    /// Geocentric position angle of Jupiter's northern
    /// rotation pole, or also called, position angle
    /// of the axis, `[0, 2π)`
    pub P: f64,
    /// Longitude of the central meridian of the illuminated disk for
    /// Rotational System I, `[0, 2π)`
    pub w1: f64,
    /// Longitude of the central meridian of the illuminated disk for
    /// Rotational System II, `[0, 2π)`
    pub w2: f64,
}

/// Quantities describing the appearance of Jupiter's disk at one instant
///
/// Returned by [`physical_ephemeris`]. All angles are in radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhysicalEphemeris {
    /// Jovicentric declination of the Earth, `D_E`, in `[-π/2, π/2]`
    /// (it never exceeds about ±3.4°)
    pub earth_dec: f64,
    /// Jovicentric declination of the Sun, `D_S`, in `[-π/2, π/2]`
    /// (it never exceeds about ±3.12°)
    pub sun_dec: f64,
    /// Position angle `P` of Jupiter's northern rotation pole, measured
    /// from the north towards the east, in `[0, 2π)`
    pub axis_pos_angle: f64,
    /// Longitude `ω1` of the central meridian of the geometric disk in
    /// System I, in `[0, 2π)`
    pub central_meridian_i: f64,
    /// Longitude `ω2` of the central meridian of the geometric disk in
    /// System II, in `[0, 2π)`
    pub central_meridian_ii: f64,
    /// Correction for phase `C`: the amount to add to a central meridian
    /// longitude of the geometric disk to get that of the illuminated
    /// disk. Its sign is that of `sin(l - l0)` (`l`, `l0` the heliocentric
    /// longitudes of Jupiter and the Earth); `|C| ≤ 0.61°`.
    pub phase_corr: f64,
    /// `ω1 + C`: System I longitude of the central meridian of the
    /// illuminated disk, in `[0, 2π)`
    pub illum_central_meridian_i: f64,
    /// `ω2 + C`: System II longitude of the central meridian of the
    /// illuminated disk, in `[0, 2π)`
    pub illum_central_meridian_ii: f64,
    /// Distance from the Earth to Jupiter, `Δ` *| in AU*
    pub earth_dist: f64,
    /// Distance from the Sun to Jupiter, `r` *| in AU*
    pub sun_dist: f64,
}

/// Results of the shorter, lower-accuracy method of Meeus, ch. 43
///
/// Returned by [`physical_ephemeris_low_accuracy`]. All angles are in
/// radians. There is no position angle of the axis in this method.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LowAccuracyEphemeris {
    /// Jovicentric declination of the Earth, `D_E`
    pub earth_dec: f64,
    /// Jovicentric declination of the Sun, `D_S`
    pub sun_dec: f64,
    /// Longitude `ω1` of the central meridian of the geometric disk in
    /// System I, in `[0, 2π)`
    pub central_meridian_i: f64,
    /// Longitude `ω2` of the central meridian of the geometric disk in
    /// System II, in `[0, 2π)`
    pub central_meridian_ii: f64,
    /// Correction for phase `C`, `±57.3° sin²(ψ/2)` with the sign opposite
    /// to that of `sin K` (`K` being about `l0 - l`)
    pub phase_corr: f64,
    /// `ω1 + C`, in `[0, 2π)`
    pub illum_central_meridian_i: f64,
    /// `ω2 + C`, in `[0, 2π)`
    pub illum_central_meridian_ii: f64,
    /// Phase angle `ψ` (the angle Earth–Jupiter–Sun), signed: it has the
    /// sign of `sin K`. Always within ±12°.
    pub phase_angle: f64,
    /// Heliocentric longitude `λ` of Jupiter, referred to the equinox of
    /// J2000.0, in `[0, 2π)`
    pub helio_long: f64,
    /// Distance from the Earth to Jupiter, `Δ` *| in AU*
    pub earth_dist: f64,
    /// Distance from the Sun to Jupiter, `r` *| in AU*
    pub sun_dist: f64,
    /// Distance from the Sun to the Earth, `R` *| in AU*
    pub earth_sun_dist: f64,
}

/// Computes the ephemeris for physical observations of Jupiter
///
/// The more accurate method of Meeus, ch. 43. The Earth and Jupiter come
/// from the crate's VSOP87 series; light-time is allowed for by
/// correcting Jupiter's heliocentric longitude and the rotation angles,
/// and the aberration and nutation are included in `P`. The nutation and
/// the mean obliquity are computed internally with [`nutation::nutation`]
/// and [`crate::ecliptic::mn_oblq_IAU`]. The aberration of the Sun as seen
/// from Jupiter is neglected in `D_S` (error below 0.5").
///
/// The result agrees with Meeus' example 43.a to the printed precision
/// (0.01°). The rotation elements are those adopted at the end of the
/// 19th century; System I and II longitudes are by definition exact, but
/// the pole has been refined since, and its expressions are linear in
/// time, so `P`, `D_E` and `D_S` are best within a few centuries of 2000
/// (by the year -1000 `D_S` is off by about half a degree).
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day (TT)
pub fn physical_ephemeris(jde: f64) -> PhysicalEphemeris {
    let (nut_in_long, nut_in_oblq) = nutation::nutation(jde);
    let s = compute(
        jde,
        planet::heliocent_coords(&planet::Planet::Earth, jde),
        planet::heliocent_coords(&planet::Planet::Jupiter, jde),
        crate::ecliptic::mn_oblq_IAU(jde),
        nut_in_long,
        nut_in_oblq,
    );

    PhysicalEphemeris {
        earth_dec: s.de,
        sun_dec: s.ds,
        axis_pos_angle: s.p,
        central_meridian_i: s.omega1,
        central_meridian_ii: s.omega2,
        phase_corr: s.phase_corr,
        illum_central_meridian_i: s.illum_omega1,
        illum_central_meridian_ii: s.illum_omega2,
        earth_dist: s.delta,
        sun_dist: s.r,
    }
}

/// Computes quantities used in the ephemeris for physical observations
/// of Jupiter
///
/// The same algorithm as [`physical_ephemeris`] (Meeus, ch. 43, more
/// accurate method), with the obliquity and nutation supplied by the
/// caller. The central meridian longitudes returned are those of the
/// illuminated disk (correction for phase included).
///
/// # Returns
///
/// * `ephemeris`: Jupiter's ephemeris. *All angles are in radians*
///
/// # Arguments
///
/// * `JD`         : Julian Ephemeris Day (TT)
/// * `mn_oblq`    : Mean obliquity of the ecliptic on `JD` *| in radians*
/// * `nut_in_long`: Nutation in ecliptic longitude on `JD` *| in radians*
/// * `nut_in_oblq`: Nutation in obliquity of the ecliptic on `JD` *| in radians*
pub fn ephemeris(JD: f64, mn_oblq: f64, nut_in_long: f64, nut_in_oblq: f64) -> Ephemeris {
    let s = compute(
        JD,
        planet::heliocent_coords(&planet::Planet::Earth, JD),
        planet::heliocent_coords(&planet::Planet::Jupiter, JD),
        mn_oblq,
        nut_in_long,
        nut_in_oblq,
    );

    Ephemeris {
        De: s.de,
        Ds: s.ds,
        P: s.p,
        w1: s.illum_omega1,
        w2: s.illum_omega2,
    }
}

/// Computes the ephemeris for physical observations of Jupiter with the
/// shorter, lower-accuracy method
///
/// Meeus, ch. 43, "lower accuracy". The positions of the Earth and Jupiter
/// come from a few periodic terms, so no planetary theory is needed. The
/// central meridian longitudes can be in error by 0.1–0.2° in the second
/// half of the 20th century (about 0.3° a century further from 2000), and
/// the declinations by a few hundredths of a degree. Light-time is
/// included (1/173 day per AU).
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day (TT)
pub fn physical_ephemeris_low_accuracy(jde: f64) -> LowAccuracyEphemeris {
    let s = compute_low_accuracy(jde);
    LowAccuracyEphemeris {
        earth_dec: s.de.to_radians(),
        sun_dec: s.ds.to_radians(),
        central_meridian_i: s.omega1.to_radians(),
        central_meridian_ii: s.omega2.to_radians(),
        phase_corr: s.phase_corr.to_radians(),
        illum_central_meridian_i: s.illum_omega1.to_radians(),
        illum_central_meridian_ii: s.illum_omega2.to_radians(),
        phase_angle: s.psi.to_radians(),
        helio_long: angle::limit_to_360(s.lambda).to_radians(),
        earth_dist: s.delta,
        sun_dist: s.r,
        earth_sun_dist: s.R,
    }
}

/// Every intermediate quantity of the accurate method of Meeus, ch. 43,
/// named after the step that produces it. Angles in radians unless noted,
/// distances in AU. Only some of them make up the public results; the
/// others are kept so that the unit tests can check every step.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(not(test), allow(dead_code))]
struct Steps {
    // step 1: days since 1950 Jan 1.0 TD and the pole, mean equinox of date
    d: f64,
    asc0: f64,
    dec0: f64,
    // step 2, in degrees, reduced to [0, 360)
    w1_deg: f64,
    w2_deg: f64,
    // step 5: geometric distance
    delta_geomet: f64,
    // steps 6-7: longitude and distance corrected for light-time
    l: f64,
    r: f64,
    delta: f64,
    // step 9: the Sun seen from Jupiter
    sun_asc: f64,
    sun_dec: f64,
    // step 10
    ds: f64,
    // step 11
    asc: f64,
    dec: f64,
    zeta: f64,
    // step 12
    de: f64,
    // steps 13-14: central meridians and the correction for phase
    omega1: f64,
    omega2: f64,
    phase_corr: f64,
    illum_omega1: f64,
    illum_omega2: f64,
    // step 16: corrections for aberration
    asc_ab_corr: f64,
    dec_ab_corr: f64,
    // step 17: apparent planet and pole
    asc_app: f64,
    dec_app: f64,
    asc0_app: f64,
    dec0_app: f64,
    // step 18
    p: f64,
}

/// Runs the accurate method of Meeus, ch. 43.
///
/// * `earth`: heliocentric `(l0, b0, R)` of the Earth at `jde`
/// * `jupiter`: heliocentric `(l, b, r)` of Jupiter at `jde` (not corrected
///   for light-time; that is done here)
fn compute(
    jde: f64,
    earth: (f64, f64, f64),
    jupiter: (f64, f64, f64),
    mn_oblq: f64,
    nut_in_long: f64,
    nut_in_oblq: f64,
) -> Steps {
    // Step 1
    let d = jde - 2_433_282.5;
    let T1 = d / 36525.0;
    let asc0 = (268.00 + 0.1061 * T1).to_radians();
    let dec0 = (64.50 - 0.0164 * T1).to_radians();

    // Step 2
    let w1_deg = angle::limit_to_360(17.710 + 877.900_035_39 * d);
    let w2_deg = angle::limit_to_360(16.838 + 870.270_035_39 * d);

    // Steps 3-5
    let (l0, _, R) = earth;
    let (l, b, r) = jupiter;
    let (_, _, _, delta_geomet) = geocentric(earth, (l, b, r));

    // Steps 6-7: allow for light-time in the longitude only
    let l = l - (0.012_990 * delta_geomet / (r * r)).to_radians();
    let (x, y, z, delta) = geocentric(earth, (l, b, r));

    // Step 9
    let sun_asc = (mn_oblq.cos() * l.sin() - mn_oblq.sin() * b.tan()).atan2(l.cos());
    let sun_dec = (mn_oblq.cos() * b.sin() + mn_oblq.sin() * b.cos() * l.sin())
        .clamp(-1.0, 1.0)
        .asin();

    // Step 10
    let ds = planetocentric_dec(asc0, dec0, sun_asc, sun_dec);

    // Step 11, as for Mars
    let u = y * mn_oblq.cos() - z * mn_oblq.sin();
    let v = y * mn_oblq.sin() + z * mn_oblq.cos();
    let asc = u.atan2(x);
    let dec = v.atan2(x.hypot(u));
    let zeta = (dec0.sin() * dec.cos() * (asc0 - asc).cos() - dec.sin() * dec0.cos())
        .atan2(dec.cos() * (asc0 - asc).sin());

    // Step 12
    let de = planetocentric_dec(asc0, dec0, asc, dec);

    // Step 13: the last terms are the rotation during the light-time
    let zeta_deg = zeta.to_degrees();
    let omega1_deg = w1_deg - zeta_deg - 5.070_33 * delta;
    let omega2_deg = w2_deg - zeta_deg - 5.026_26 * delta;

    // Step 14: the correction for phase. 57.2958 × (...) in degrees is
    // (...) in radians.
    let phase_corr = ((2.0 * r * delta + R * R - r * r - delta * delta) / (4.0 * r * delta))
        .copysign((l - l0).sin());
    let omega1 = angle::limit_to_two_PI(omega1_deg.to_radians());
    let omega2 = angle::limit_to_two_PI(omega2_deg.to_radians());
    let illum_omega1 = angle::limit_to_two_PI(omega1 + phase_corr);
    let illum_omega2 = angle::limit_to_two_PI(omega2 + phase_corr);

    // Step 15
    let tru_oblq = mn_oblq + nut_in_oblq;

    // Step 16: Jupiter's aberration; both corrections use the uncorrected
    // α and δ
    let k_ab = ABERRATION.to_radians();
    let asc_ab_corr =
        k_ab * (asc.cos() * l0.cos() * tru_oblq.cos() + asc.sin() * l0.sin()) / dec.cos();
    let dec_ab_corr = k_ab
        * (l0.cos() * tru_oblq.cos() * (tru_oblq.tan() * dec.cos() - asc.sin() * dec.sin())
            + asc.cos() * dec.sin() * l0.sin());
    let (asc_ab, dec_ab) = (asc + asc_ab_corr, dec + dec_ab_corr);

    // Step 17: nutation, eq. 23.1
    let (d_asc, d_dec) = nutation::nutation_in_eq_coords(
        &coords::EqPoint {
            asc: asc_ab,
            dec: dec_ab,
        },
        nut_in_long,
        nut_in_oblq,
        tru_oblq,
    );
    let (asc_app, dec_app) = (asc_ab + d_asc, dec_ab + d_dec);
    let (d_asc0, d_dec0) = nutation::nutation_in_eq_coords(
        &coords::EqPoint {
            asc: asc0,
            dec: dec0,
        },
        nut_in_long,
        nut_in_oblq,
        tru_oblq,
    );
    let (asc0_app, dec0_app) = (asc0 + d_asc0, dec0 + d_dec0);

    // Step 18, eq. 42.4
    let p = angle::limit_to_two_PI((dec0_app.cos() * (asc0_app - asc_app).sin()).atan2(
        dec0_app.sin() * dec_app.cos()
            - dec0_app.cos() * dec_app.sin() * (asc0_app - asc_app).cos(),
    ));

    Steps {
        d,
        asc0,
        dec0,
        w1_deg,
        w2_deg,
        delta_geomet,
        l,
        r,
        delta,
        sun_asc,
        sun_dec,
        ds,
        asc,
        dec,
        zeta,
        de,
        omega1,
        omega2,
        phase_corr,
        illum_omega1,
        illum_omega2,
        asc_ab_corr,
        dec_ab_corr,
        asc_app,
        dec_app,
        asc0_app,
        dec0_app,
        p,
    }
}

/// Intermediate quantities of the lower-accuracy method, in degrees and AU,
/// as printed in Meeus' example 43.b.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(not(test), allow(dead_code))]
struct LowAccuracySteps {
    d: f64,
    V: f64,
    M: f64,
    N: f64,
    J: f64,
    A: f64,
    B: f64,
    K: f64,
    R: f64,
    r: f64,
    delta: f64,
    psi: f64,
    d_light: f64,
    omega1: f64,
    omega2: f64,
    phase_corr: f64,
    illum_omega1: f64,
    illum_omega2: f64,
    lambda: f64,
    ds: f64,
    de: f64,
}

/// Runs the lower-accuracy method of Meeus, ch. 43. All angles in degrees.
fn compute_low_accuracy(jde: f64) -> LowAccuracySteps {
    let sin = |x: f64| x.to_radians().sin();
    let cos = |x: f64| x.to_radians().cos();

    // days since J2000.0
    let d = jde - 2_451_545.0;
    // argument of the long-period term in Jupiter's motion
    let V = 172.74 + 0.001_115_88 * d;
    // mean anomalies of the Earth and Jupiter
    let M = 357.529 + 0.985_600_3 * d;
    let N = 20.020 + 0.083_085_3 * d + 0.329 * sin(V);
    // difference of the mean heliocentric longitudes, Earth minus Jupiter
    let J = 66.115 + 0.902_517_9 * d - 0.329 * sin(V);
    // equations of the centre
    let A = 1.915 * sin(M) + 0.020 * sin(2.0 * M);
    let B = 5.555 * sin(N) + 0.168 * sin(2.0 * N);
    let K = J + A - B;
    // radius vectors and distance
    let R = 1.000_14 - 0.016_71 * cos(M) - 0.000_14 * cos(2.0 * M);
    let r = 5.208_72 - 0.252_08 * cos(N) - 0.006_11 * cos(2.0 * N);
    let delta = (r * r + R * R - 2.0 * r * R * cos(K)).sqrt();
    // phase angle, signed like sin K
    let psi = (R / delta * sin(K)).clamp(-1.0, 1.0).asin().to_degrees();

    // central meridians of the geometric disk; Δ/173 is the light-time
    let d_light = d - delta / 173.0;
    let omega1 = angle::limit_to_360(210.98 + 877.816_908_8 * d_light + psi - B);
    let omega2 = angle::limit_to_360(187.23 + 870.186_908_8 * d_light + psi - B);

    // correction for phase, with the sign opposite to that of sin K
    let phase_corr = (57.3 * sin(psi / 2.0).powi(2)).copysign(-sin(K));
    let illum_omega1 = angle::limit_to_360(omega1 + phase_corr);
    let illum_omega2 = angle::limit_to_360(omega2 + phase_corr);

    // heliocentric longitude of Jupiter, equinox of J2000.0
    let lambda = 34.35 + 0.083_091 * d + 0.329 * sin(V) + B;
    // the amplitudes are tilts: equator on orbit (3.12°), equator on
    // ecliptic (2.22°) and orbit on ecliptic (1.30°)
    let ds = 3.12 * sin(lambda + 42.8);
    let de = ds
        - 2.22 * sin(psi) * cos(lambda + 22.0)
        - 1.30 * (r - delta) / delta * sin(lambda - 100.5);

    LowAccuracySteps {
        d,
        V,
        M,
        N,
        J,
        A,
        B,
        K,
        R,
        r,
        delta,
        psi,
        d_light,
        omega1,
        omega2,
        phase_corr,
        illum_omega1,
        illum_omega2,
        lambda,
        ds,
        de,
    }
}

/// Planetocentric declination of the direction `(asc, dec)` for a planet
/// whose north pole is at `(asc0, dec0)` (Meeus, ch. 43, steps 10 and 12).
fn planetocentric_dec(asc0: f64, dec0: f64, asc: f64, dec: f64) -> f64 {
    (-dec0.sin() * dec.sin() - dec0.cos() * dec.cos() * (asc0 - asc).cos())
        .clamp(-1.0, 1.0)
        .asin()
}

/// Geocentric ecliptic rectangular coordinates and distance of a planet
/// (Meeus, eq. 42.2 and 42.3), from heliocentric `(l, b, r)` of the Earth
/// and of the planet.
fn geocentric(earth: (f64, f64, f64), planet: (f64, f64, f64)) -> (f64, f64, f64, f64) {
    let (l0, b0, R) = earth;
    let (l, b, r) = planet;
    let x = r * b.cos() * l.cos() - R * b0.cos() * l0.cos();
    let y = r * b.cos() * l.sin() - R * b0.cos() * l0.sin();
    let z = r * b.sin() - R * b0.sin();
    (x, y, z, (x * x + y * y + z * z).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[track_caller]
    fn close(actual_rad: f64, expected_deg: f64, tol_deg: f64) {
        close_deg(actual_rad.to_degrees(), expected_deg, tol_deg);
    }

    #[track_caller]
    fn close_deg(actual: f64, expected: f64, tol: f64) {
        let diff = (actual - expected + 180.0).rem_euclid(360.0) - 180.0;
        assert!(
            diff.abs() <= tol,
            "expected {expected}° ± {tol}, got {actual}° (off by {diff:e})"
        );
    }

    #[track_caller]
    fn close_val(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() <= tol,
            "expected {expected} ± {tol}, got {actual}"
        );
    }

    const JDE: f64 = 2_448_972.500_68;

    /// Meeus' example 43.a from the book's own inputs: the Earth and
    /// Jupiter from the full VSOP87 (steps 3-4), its obliquity (step 8) and
    /// nutation (step 15).
    fn book_steps() -> Steps {
        compute(
            JDE,
            (
                84.285_703_f64.to_radians(),
                0.000_197_f64.to_radians(),
                0.984_123_16,
            ),
            (
                181.882_168_f64.to_radians(),
                1.290_464_f64.to_radians(),
                5.446_423_20,
            ),
            (23.0 + 26.0 / 60.0 + 24.745 / 3600.0_f64).to_radians(),
            (16.86 / 3600.0_f64).to_radians(),
            (-1.79 / 3600.0_f64).to_radians(),
        )
    }

    // Meeus, example 43.a, steps 1 to 12
    #[test]
    fn example_43a_steps_1_to_12() {
        let s = book_steps();
        close_val(s.d, 15_690.000_68, 1e-6);
        close(s.asc0, 268.045_58, 5e-6);
        close(s.dec0, 64.492_96, 5e-6);
        close_deg(s.w1_deg, 309.8622, 5e-5);
        close_deg(s.w2_deg, 114.2851, 5e-5);
        close_val(s.delta_geomet, 5.661_164_5, 1e-7);
        close(s.l, 181.879_689, 1e-6);
        close_val(s.r, 5.446_423_20, 1e-9);
        close_val(s.delta, 5.661_123_9, 1e-7);
        close(s.sun_asc, 182.237_749, 2e-6);
        close(s.sun_dec, 0.436_472, 2e-6);
        close(s.ds, -2.20, 0.005);
        close(s.asc, 191.340_327, 2e-6);
        close(s.dec, -3.524_749, 2e-6);
        close(s.zeta, 13.5238, 5e-5);
        close(s.de, -2.48, 0.005);
    }

    // Meeus, example 43.a, steps 13 to 18
    #[test]
    fn example_43a_steps_13_to_18() {
        let s = book_steps();
        close(s.omega1, 267.63, 0.005);
        close(s.omega2, 72.31, 0.005);
        close(s.phase_corr, 0.43, 0.005);
        close(s.illum_omega1, 268.06, 0.005);
        close(s.illum_omega2, 72.74, 0.005);
        close(s.asc_ab_corr, -0.001_627, 5e-7);
        close(s.dec_ab_corr, 0.000_560, 5e-7);
        close(s.asc_app, 191.343_05, 5e-6);
        close(s.dec_app, -3.525_92, 5e-6);
        close(s.asc0_app, 268.045_94, 5e-6);
        close(s.dec0_app, 64.493_39, 5e-6);
        close(s.p, 24.80, 0.005);
    }

    // The δ aberration correction of step 16 has the term sin α sin δ.
    // The old code had sin α cos α there and gave +0.000466° instead of
    // the book's +0.000560°.
    #[test]
    fn example_43a_aberration_in_declination() {
        let s = book_steps();
        close(s.dec_ab_corr, 0.000_560, 5e-7);
    }

    // Meeus, example 43.b
    #[test]
    fn example_43b_steps() {
        let s = compute_low_accuracy(JDE);
        close_val(s.d, -2_572.499_32, 1e-6);
        close_deg(s.V, 169.87, 0.005);
        close_deg(s.M, 342.073, 5e-4);
        close_deg(s.N, -193.659, 5e-4);
        close_deg(s.J, 264.330, 5e-4);
        close_val(s.A, -0.601, 5e-4);
        close_val(s.B, 1.235, 5e-4);
        // the book adds J, A and B rounded to 0.001°; unrounded K = 262.4946°
        close_deg(s.K, 262.494, 1e-3);
        close_val(s.R, 0.984_13, 5e-6);
        close_val(s.r, 5.448_24, 5e-6);
        // from the rounded K the book gets 5.661507; unrounded, 5.661497
        close_val(s.delta, 5.661_51, 2e-5);
        close_val(s.psi.to_radians().sin(), -0.172_34, 5e-6);
        close_val(s.psi, -9.924, 5e-4);
        close_val(s.d_light, -2_572.532_05, 5e-6);
        close_deg(s.omega1, 267.69, 0.005);
        close_deg(s.omega2, 72.36, 0.005);
        close_val(s.phase_corr, 0.43, 0.005);
        close_deg(s.lambda, -178.11, 0.005);
        close_val(s.ds, -2.194, 5e-4);
        close_val(s.de, -2.50, 0.005);
    }
}

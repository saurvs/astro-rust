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

//! Mars: ephemeris for physical observations
//!
//! Implements Meeus, *Astronomical Algorithms*, 2nd ed., chapter 42. For a
//! given instant, [`physical_ephemeris`] describes how the disk of Mars is
//! oriented and lit as seen from the Earth:
//!
//! * the planetocentric declinations of the Earth (`D_E`) and of the Sun
//!   (`D_S`): how far Mars' north pole is tilted towards (positive) or away
//!   from (negative) each of them,
//! * the position angle `P` of Mars' rotation axis on the sky,
//! * the areographic longitude `ω` of the central meridian,
//! * the size `q` and position angle `Q` of the greatest defect of
//!   illumination, the illuminated fraction `k`, and the apparent diameter.
//!
//! Positions of the Earth and Mars come from the crate's VSOP87 series
//! ([`planet::heliocent_coords`]), corrected for light-time; nutation comes
//! from [`nutation::nutation`] and the mean obliquity from
//! [`ecliptic::mn_oblq_IAU`]. [`ephemeris`] is the older entry point that
//! lets the caller supply the obliquity, the nutation and the pole.
//!
//! Position angles are measured from the north point of the disk towards
//! the east (counter-clockwise on the sky), in `[0, 2π)`.
//!
//! # Example
//!
//! Meeus' example 42.a, 1992 November 9 at 0h UT (JDE 2448935.500683):
//!
//! ```
//! use astro::planet::mars;
//!
//! let eph = mars::physical_ephemeris(2448935.500683);
//!
//! assert!((eph.earth_dec.to_degrees() - 12.44).abs() < 0.005);
//! assert!((eph.sun_dec.to_degrees() - -2.76).abs() < 0.005);
//! assert!((eph.axis_pos_angle.to_degrees() - 347.64).abs() < 0.005);
//! assert!((eph.central_meridian.to_degrees() - 111.55).abs() < 0.005);
//! // apparent diameter 10.75", greatest defect of illumination 1.06"
//! assert!((eph.apparent_diameter.to_degrees() * 3600.0 - 10.75).abs() < 0.005);
//! assert!((eph.defect_of_illum.to_degrees() * 3600.0 - 1.06).abs() < 0.005);
//! ```

use core::f64::consts::PI;

use crate::angle;
use crate::coords;
use crate::ecliptic;
use crate::nutation;
use crate::planet;
use crate::time;

/// Light-time for a distance of one AU, in days (Meeus, eq. 33.3).
const LIGHT_TIME_PER_AU: f64 = 0.005_775_518_3;

/// Apparent diameter of Mars at a distance of one AU, in arcseconds.
const DIAMETER_AT_1_AU: f64 = 9.36;

/// The constant of annual aberration (about 20.49″), in degrees, as used
/// by Meeus for the aberration of the planet in chapters 42, 43 and 45.
const ABERRATION: f64 = 0.005_693;

/// Returns the equatorial coordinates of Mars's north pole for the epoch
/// J1950.0
///
/// The values are the IAU 1982 ones quoted by Meeus (ch. 42), referred to
/// the mean equator and equinox of B1950.0.
///
/// # Returns
///
/// * `eq_coords`: Equatorial coordinates of Mars's north pole for the
///   epoch J1950.0 *| in radians*
#[inline(always)]
pub fn north_pol_eq_coords_J1950() -> coords::EqPoint {
    coords::EqPoint {
        asc: 317.342_f64.to_radians(),
        dec: 52.7110_f64.to_radians(),
    }
}

/// Returns the equatorial coordinates of Mars's north pole for the epoch
/// J2000.0
///
/// The values are the IAU 1982 ones quoted by Meeus (ch. 42), referred to
/// the mean equator and equinox of J2000.0.
///
/// # Returns
///
/// * `eq_coords`: Equatorial coordinates of Mars's north pole for the
///   epoch J2000.0 *| in radians*
#[inline(always)]
pub fn north_pol_eq_coords_J2000() -> coords::EqPoint {
    coords::EqPoint {
        asc: 317.681_f64.to_radians(),
        dec: 52.8860_f64.to_radians(),
    }
}

/// Computes the ecliptic coordinates of Mars's north pole, referred to
/// the mean equinox of the date
///
/// Meeus, eq. 42.1. The expressions include the precession of the
/// rotation axes of both the Earth and Mars.
///
/// # Returns
///
/// * `ecl_coords`: Ecliptic longitude `λ0`, in `[0, 2π)`, and latitude
///   `β0` of Mars's north pole *| in radians*
///
/// # Arguments
///
/// * `JC`: Julian (Ephemeris) centuries since J2000.0, see
///   [`time::julian_cent`]
#[inline(always)]
pub fn north_pol_ecl_coords(JC: f64) -> coords::EclPoint {
    coords::EclPoint {
        long: angle::limit_to_360(352.9065 + 1.17330 * JC).to_radians(),
        lat: (63.28180 - 0.00394 * JC).to_radians(),
    }
}

/// Holds Mars's ephemeris values for physical observations
///
/// Returned by [`ephemeris`]; [`PhysicalEphemeris`] holds the same values
/// and a few more. All angles are in radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ephemeris {
    /// Mars-centric declination of the Earth, `[-π/2, π/2]`
    pub De: f64,
    /// Mars-centric declination of the Sun, `[-π/2, π/2]`
    pub Ds: f64,
    /// Geocentric position angle of Mars's northern rotation pole,
    /// or also called the position angle of the axis, `[0, 2π)`
    pub P: f64,
    /// Angular amount of the greatest defect of illumination
    pub q: f64,
    /// Longitude of the central meridian, as seen from the Earth,
    /// `[0, 2π)`
    pub w: f64,
    /// Apparent diameter of Mars
    pub d: f64,
}

/// Quantities describing the appearance of Mars' disk at one instant
///
/// Returned by [`physical_ephemeris`]. All angles are in radians;
/// position angles are measured from the north point of the disk towards
/// the east.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhysicalEphemeris {
    /// Planetocentric declination of the Earth, `D_E`, in `[-π/2, π/2]`.
    /// Positive when Mars' north pole is tilted towards the Earth.
    pub earth_dec: f64,
    /// Planetocentric declination of the Sun, `D_S`, in `[-π/2, π/2]`.
    /// Positive when Mars' north pole is lit.
    pub sun_dec: f64,
    /// Position angle `P` of Mars' northern rotation pole, in `[0, 2π)`
    pub axis_pos_angle: f64,
    /// Areographic longitude `ω` of the central meridian as seen from the
    /// Earth, in `[0, 2π)`. Areographic longitudes increase towards the
    /// west, so `ω` increases with time as Mars rotates.
    pub central_meridian: f64,
    /// Position angle `Q` of the greatest defect of illumination, in
    /// `[0, 2π)`. It is opposite (180° away from) the position angle of the
    /// midpoint of the bright limb.
    pub defect_pos_angle: f64,
    /// Angular width `q` of the greatest defect of illumination,
    /// `(1 - k) d`
    pub defect_of_illum: f64,
    /// Apparent (equatorial) diameter `d` of the disk
    pub apparent_diameter: f64,
    /// Illuminated fraction `k` of the disk, in `[0, 1]` (Meeus, eq. 41.2)
    pub illum_frac: f64,
    /// Distance from the Earth to Mars, `Δ`, at the time the light left
    /// Mars *| in AU*
    pub earth_dist: f64,
    /// Distance from the Sun to Mars, `r`, at the time the light left
    /// Mars *| in AU*
    pub sun_dist: f64,
    /// Light-time from Mars to the Earth, `τ` *| in days*
    pub light_time: f64,
}

/// Computes the ephemeris for physical observations of Mars
///
/// Meeus, ch. 42. The Earth and Mars come from the crate's VSOP87 series,
/// Mars being taken at the time `jde - τ` its light left it (the
/// light-time `τ` is found by iteration); the aberration of the Sun as
/// seen from Mars is included in `D_S`, and the nutation and the
/// aberration of Mars in `P`. The nutation and mean obliquity of the date
/// are computed internally with [`nutation::nutation`] and
/// [`ecliptic::mn_oblq_IAU`].
///
/// The result agrees with Meeus' example 42.a to the printed precision
/// (0.01°, and 0.01" for `q` and `d`). The accuracy is limited by the
/// pole and rotation elements (IAU 1982) rather than by the positions,
/// so the result is best within a few centuries of 2000.
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day (TT)
pub fn physical_ephemeris(jde: f64) -> PhysicalEphemeris {
    let pole = north_pol_ecl_coords(time::julian_cent(jde));
    let (nut_in_long, nut_in_oblq) = nutation::nutation(jde);
    let s = compute(
        jde,
        (pole.long, pole.lat),
        planet::heliocent_coords(&planet::Planet::Earth, jde),
        |t| planet::heliocent_coords(&planet::Planet::Mars, t),
        ecliptic::mn_oblq_IAU(jde),
        nut_in_long,
        nut_in_oblq,
    );

    PhysicalEphemeris {
        earth_dec: s.de,
        sun_dec: s.ds,
        axis_pos_angle: s.p,
        central_meridian: s.omega,
        defect_pos_angle: s.defect_pos_angle,
        defect_of_illum: s.q,
        apparent_diameter: s.d,
        illum_frac: s.k,
        earth_dist: s.delta,
        sun_dist: s.r,
        light_time: s.tau,
    }
}

/// Computes quantities used in the ephemeris for physical observations of
/// Mars
///
/// The same algorithm as [`physical_ephemeris`] (Meeus, ch. 42), with the
/// pole, obliquity and nutation supplied by the caller. Positions of the
/// Earth and Mars come from the crate's VSOP87 series.
///
/// # Returns
///
/// * `ephemeris`: Mars's ephemeris. *All angles are in radians*
///
/// # Arguments
///
/// * `JD`                   : Julian Ephemeris Day (TT)
/// * `north_pole_ecl_coords`: Ecliptic coordinates of Mars's north pole
///   on `JD`, see [`north_pol_ecl_coords`] *| in radians*
/// * `mn_oblq`              : Mean obliquity of the ecliptic on `JD`
///   *| in radians*
/// * `nut_in_long`          : Nutation in ecliptic longitude on `JD`
///   *| in radians*
/// * `nut_in_oblq`          : Nutation in obliquity of the ecliptic on
///   `JD` *| in radians*
pub fn ephemeris(
    JD: f64,
    north_pole_ecl_coords: &coords::EclPoint,
    mn_oblq: f64,
    nut_in_long: f64,
    nut_in_oblq: f64,
) -> Ephemeris {
    let s = compute(
        JD,
        (north_pole_ecl_coords.long, north_pole_ecl_coords.lat),
        planet::heliocent_coords(&planet::Planet::Earth, JD),
        |t| planet::heliocent_coords(&planet::Planet::Mars, t),
        mn_oblq,
        nut_in_long,
        nut_in_oblq,
    );

    Ephemeris {
        De: s.de,
        Ds: s.ds,
        P: s.p,
        q: s.q,
        w: s.omega,
        d: s.d,
    }
}

/// Every intermediate quantity of Meeus' procedure in ch. 42, named after
/// the step that produces it. Angles in radians, distances in AU, times in
/// days. Only some of them make up the public results; the others are
/// kept so that the unit tests can check every step against the book.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(not(test), allow(dead_code))]
struct Steps {
    // steps 3-4: Mars at t - τ, geocentric rectangular coordinates
    l: f64,
    b: f64,
    r: f64,
    delta: f64,
    tau: f64,
    // step 5: geocentric ecliptic coordinates, mean equinox of date
    lambda: f64,
    beta: f64,
    // step 6
    de: f64,
    // step 7: Mars' node and the heliocentric position corrected for the
    // aberration of the Sun as seen from Mars
    node: f64,
    l1: f64,
    b1: f64,
    // step 8
    ds: f64,
    // step 9: rotation angle of Mars' prime meridian
    w_angle: f64,
    // step 10: pole, mean equator and equinox of date
    asc0: f64,
    dec0: f64,
    // step 11
    asc: f64,
    dec: f64,
    zeta: f64,
    // step 12
    omega: f64,
    // step 14: aberration-corrected geocentric position
    lambda_ab: f64,
    beta_ab: f64,
    // step 15: longitudes corrected for nutation, true obliquity
    lambda0_app: f64,
    lambda_app: f64,
    tru_oblq: f64,
    // step 16: apparent pole and planet, true equator and equinox
    asc0_app: f64,
    dec0_app: f64,
    asc_app: f64,
    dec_app: f64,
    // step 17
    p: f64,
    // step 18: the Sun, and the position angles of the bright limb and of
    // the greatest defect of illumination
    sun_asc: f64,
    sun_dec: f64,
    bright_limb_pos_angle: f64,
    defect_pos_angle: f64,
    // step 19
    k: f64,
    d: f64,
    q: f64,
}

/// Runs Meeus' ch. 42 procedure.
///
/// * `pole`: ecliptic `(λ0, β0)` of Mars' north pole, mean equinox of date
/// * `earth`: heliocentric `(l0, b0, R)` of the Earth at `jde`
/// * `mars_at`: heliocentric `(l, b, r)` of Mars at a given JDE
fn compute(
    jde: f64,
    pole: (f64, f64),
    earth: (f64, f64, f64),
    mars_at: impl Fn(f64) -> (f64, f64, f64),
    mn_oblq: f64,
    nut_in_long: f64,
    nut_in_oblq: f64,
) -> Steps {
    let (lambda0, beta0) = pole;
    let (l0, b0, R) = earth;

    // Steps 3-4
    let pos = light_time_corrected(jde, earth, mars_at);
    let (l, b, r) = (pos.l, pos.b, pos.r);
    let (x, y, z) = (pos.x, pos.y, pos.z);
    let (delta, tau) = (pos.delta, pos.tau);

    // Step 5
    let lambda = y.atan2(x);
    let beta = z.atan2(x.hypot(y));

    // Step 6
    let de = planetocentric_dec(beta0, lambda0, beta, lambda);

    // Steps 7-8
    let T = time::julian_cent(jde);
    let node = (49.5581 + 0.7721 * T).to_radians();
    let l1 = l - (0.00697 / r).to_radians();
    let b1 = b - (0.000225 * (l - node).cos() / r).to_radians();
    let ds = planetocentric_dec(beta0, lambda0, b1, l1);

    // Step 9
    let w_angle = (11.504 + 350.892_000_25 * (jde - tau - 2_433_282.5))
        .rem_euclid(360.0)
        .to_radians();

    // Step 10
    let asc0 = coords::asc_frm_ecl(lambda0, beta0, mn_oblq);
    let dec0 = coords::dec_frm_ecl(lambda0, beta0, mn_oblq);

    // Step 11
    let u = y * mn_oblq.cos() - z * mn_oblq.sin();
    let v = y * mn_oblq.sin() + z * mn_oblq.cos();
    let asc = u.atan2(x);
    let dec = v.atan2(x.hypot(u));
    let zeta = (dec0.sin() * dec.cos() * (asc0 - asc).cos() - dec.sin() * dec0.cos())
        .atan2(dec.cos() * (asc0 - asc).sin());

    // Step 12
    let omega = angle::limit_to_two_PI(w_angle - zeta);

    // Step 14: both corrections use the uncorrected λ and β
    let k_ab = ABERRATION.to_radians();
    let lambda_ab = lambda + k_ab * (l0 - lambda).cos() / beta.cos();
    let beta_ab = beta + k_ab * (l0 - lambda).sin() * beta.sin();

    // Steps 15-16
    let tru_oblq = mn_oblq + nut_in_oblq;
    let (lambda0_app, lambda_app) = (lambda0 + nut_in_long, lambda_ab + nut_in_long);
    let asc0_app = coords::asc_frm_ecl(lambda0_app, beta0, tru_oblq);
    let dec0_app = coords::dec_frm_ecl(lambda0_app, beta0, tru_oblq);
    let asc_app = coords::asc_frm_ecl(lambda_app, beta_ab, tru_oblq);
    let dec_app = coords::dec_frm_ecl(lambda_app, beta_ab, tru_oblq);

    // Step 17, eq. 42.4
    let p = angle::limit_to_two_PI(pos_angle(asc0_app, dec0_app, asc_app, dec_app));

    // Step 18: the Sun seen from the Earth is opposite the Earth seen from
    // the Sun. It is referred to the same frame (mean equinox of date) as
    // the α, δ of Mars from step 11; the position angle of the bright limb
    // follows from eq. 48.5, and the greatest defect is opposite it.
    let sun_long = l0 + PI;
    let sun_asc = coords::asc_frm_ecl(sun_long, -b0, mn_oblq);
    let sun_dec = coords::dec_frm_ecl(sun_long, -b0, mn_oblq);
    let bright_limb_pos_angle = angle::limit_to_two_PI(pos_angle(sun_asc, sun_dec, asc, dec));
    let defect_pos_angle = angle::limit_to_two_PI(bright_limb_pos_angle + PI);

    // Step 19, with the illuminated fraction from eq. 41.2
    let k = (((r + delta).powi(2) - R * R) / (4.0 * r * delta)).clamp(0.0, 1.0);
    let d = (DIAMETER_AT_1_AU / 3600.0).to_radians() / delta;
    let q = (1.0 - k) * d;

    Steps {
        l,
        b,
        r,
        delta,
        tau,
        lambda,
        beta,
        de,
        node,
        l1,
        b1,
        ds,
        w_angle,
        asc0,
        dec0,
        asc,
        dec,
        zeta,
        omega,
        lambda_ab,
        beta_ab,
        lambda0_app,
        lambda_app,
        tru_oblq,
        asc0_app,
        dec0_app,
        asc_app,
        dec_app,
        p,
        sun_asc,
        sun_dec,
        bright_limb_pos_angle,
        defect_pos_angle,
        k,
        d,
        q,
    }
}

/// Planetocentric declination of a body seen in the direction
/// `(lambda, beta)` from a planet whose north pole is at
/// `(lambda0, beta0)` (Meeus, ch. 42, steps 6 and 8).
fn planetocentric_dec(beta0: f64, lambda0: f64, beta: f64, lambda: f64) -> f64 {
    (-beta0.sin() * beta.sin() - beta0.cos() * beta.cos() * (lambda0 - lambda).cos())
        .clamp(-1.0, 1.0)
        .asin()
}

/// Position angle of the direction `(asc0, dec0)` as seen at `(asc, dec)`,
/// measured from the north towards the east (Meeus, eq. 42.4 and 48.5).
/// The result is in `(-π, π]`.
fn pos_angle(asc0: f64, dec0: f64, asc: f64, dec: f64) -> f64 {
    (dec0.cos() * (asc0 - asc).sin())
        .atan2(dec0.sin() * dec.cos() - dec0.cos() * dec.sin() * (asc0 - asc).cos())
}

/// A planet's heliocentric and geocentric positions, corrected for
/// light-time.
#[derive(Debug, Clone, Copy)]
struct LightTimePos {
    l: f64,
    b: f64,
    r: f64,
    x: f64,
    y: f64,
    z: f64,
    delta: f64,
    tau: f64,
}

/// Iterates Meeus' steps 3 and 4: the planet is taken at `jde - τ`, where
/// `τ` is the light-time for the planet's distance from the Earth at
/// `jde`. The iteration stops when `τ` changes by less than 1e-9 day
/// (0.1 ms).
fn light_time_corrected(
    jde: f64,
    earth: (f64, f64, f64),
    planet_at: impl Fn(f64) -> (f64, f64, f64),
) -> LightTimePos {
    let (l0, b0, R) = earth;
    let at = |tau: f64| {
        let (l, b, r) = planet_at(jde - tau);
        // eq. 42.2 and 42.3
        let x = r * b.cos() * l.cos() - R * b0.cos() * l0.cos();
        let y = r * b.cos() * l.sin() - R * b0.cos() * l0.sin();
        let z = r * b.sin() - R * b0.sin();
        let delta = (x * x + y * y + z * z).sqrt();
        LightTimePos {
            l,
            b,
            r,
            x,
            y,
            z,
            delta,
            tau,
        }
    };

    let mut pos = at(0.0);
    for _ in 0..10 {
        let tau = LIGHT_TIME_PER_AU * pos.delta;
        if (tau - pos.tau).abs() < 1e-9 {
            break;
        }
        pos = at(tau);
    }
    pos
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dms(d: f64, m: f64, s: f64) -> f64 {
        d.signum() * (d.abs() + m / 60.0 + s / 3600.0)
    }

    #[track_caller]
    fn close(actual_rad: f64, expected_deg: f64, tol_deg: f64) {
        let a = actual_rad.to_degrees();
        let diff = (a - expected_deg + 180.0).rem_euclid(360.0) - 180.0;
        assert!(
            diff.abs() <= tol_deg,
            "expected {expected_deg}° ± {tol_deg}, got {a}° (off by {diff:e})"
        );
    }

    #[track_caller]
    fn close_val(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() <= tol,
            "expected {expected} ± {tol}, got {actual}"
        );
    }

    /// Three-point interpolation (Meeus, eq. 3.3) in a table of Mars'
    /// heliocentric coordinates at 1992 Nov 8.0, 9.0 and 10.0 TD, as
    /// printed in Meeus' example 42.a (from the full VSOP87).
    fn mars_book(jde: f64) -> (f64, f64, f64) {
        let table = [
            (dms(77.0, 57.0, 48.45), dms(0.0, 52.0, 54.74), 1.540_379_7),
            (dms(78.0, 28.0, 24.28), dms(0.0, 53.0, 46.72), 1.541_658_5),
            (dms(78.0, 58.0, 57.09), dms(0.0, 54.0, 38.36), 1.542_934_7),
        ];
        let n = jde - 2_448_935.5;
        let interp = |y1: f64, y2: f64, y3: f64| {
            let (a, b) = (y2 - y1, y3 - y2);
            y2 + n / 2.0 * (a + b + n * (b - a))
        };
        let (t0, t1, t2) = (table[0], table[1], table[2]);
        (
            interp(t0.0, t1.0, t2.0).to_radians(),
            interp(t0.1, t1.1, t2.1).to_radians(),
            interp(t0.2, t1.2, t2.2),
        )
    }

    const JDE: f64 = 2_448_935.500_683;

    /// Meeus' example 42.a, computed from the book's own inputs (Earth and
    /// Mars from the full VSOP87, its nutation and obliquity), so every
    /// intermediate value can be compared at the precision printed.
    fn book_steps() -> Steps {
        let pole = north_pol_ecl_coords(time::julian_cent(JDE));
        let earth = (
            46.843_861_f64.to_radians(),
            (-0.000_167_f64).to_radians(),
            0.990_413_01,
        );
        compute(
            JDE,
            (pole.long, pole.lat),
            earth,
            mars_book,
            dms(23.0, 26.0, 24.793).to_radians(),
            (15.42_f64 / 3600.0).to_radians(),
            (-1.00_f64 / 3600.0).to_radians(),
        )
    }

    // Meeus, example 42.a, step 1
    #[test]
    fn example_42a_pole() {
        let T = time::julian_cent(JDE);
        close_val(T, -0.071_444_197_6, 1e-10);
        let pole = north_pol_ecl_coords(T);
        close(pole.long, 352.822_67, 5e-6);
        close(pole.lat, 63.282_08, 5e-6);
    }

    // Meeus, example 42.a, the interpolated table values of step 3
    #[test]
    fn example_42a_book_table() {
        let (l, b, r) = mars_book(JDE);
        close(l, 78.473_759, 2e-6);
        close(b, 0.896_321, 2e-6);
        close_val(r, 1.541_659_4, 1e-7);
    }

    // Meeus, example 42.a, steps 3 to 12.
    //
    // Tolerances: the book's inputs are rounded to 1e-6° and 1e-7 AU, which
    // leaves a few 1e-6° of play in the geocentric angles. The book also
    // stopped the light-time iteration one step early (τ = 0.005028 d from
    // the first Δ; the converged value is 0.0050275 d), which moves W by
    // 0.00018°.
    #[test]
    fn example_42a_steps_3_to_12() {
        let s = book_steps();
        close_val(s.tau, 0.005_028, 1e-6);
        close(s.l, 78.471_197, 2e-6);
        close(s.b, 0.896_249, 2e-6);
        close_val(s.r, 1.541_652_9, 1e-7);
        close_val(s.delta, 0.870_480_1, 1e-7);
        close(s.lambda, 115.117_321, 5e-6);
        close(s.beta, 1.587_619, 5e-6);
        close(s.de, 12.44, 0.005);
        close(s.node, 49.5029, 5e-5);
        close(s.l1, 78.466_676, 2e-6);
        close(s.b1, 0.896_121, 2e-6);
        close(s.ds, -2.76, 0.005);
        close(s.w_angle, 2.4593, 2.5e-4);
        close(s.asc0, 317.632_606, 5e-6);
        close(s.dec0, 52.860_916, 5e-6);
        close(s.asc, 117.377_075, 5e-6);
        close(s.dec, 22.672_176, 5e-6);
        close(s.zeta, 250.9052, 5e-5);
        close(s.omega, 111.55, 0.005);
    }

    // Meeus, example 42.a, steps 14 to 19
    #[test]
    fn example_42a_steps_14_to_19() {
        let s = book_steps();
        close(s.lambda_ab, 115.119_429, 5e-6);
        close(s.beta_ab, 1.587_472, 5e-6);
        // the book added Δψ to λ0 rounded to five decimals
        close(s.lambda0_app, 352.826_95, 1e-5);
        close(s.lambda_app, 115.123_712, 5e-6);
        close(s.tru_oblq, 23.439_942, 5e-7);
        close(s.asc0_app, 317.635_29, 1e-5);
        close(s.dec0_app, 52.862_36, 1e-5);
        close(s.asc_app, 117.383_80, 1e-5);
        close(s.dec_app, 22.670_62, 1e-5);
        close(s.p, 347.64, 0.005);
        close(s.sun_asc, 224.378, 5e-4);
        close(s.sun_dec, -16.869, 5e-4);
        close(s.bright_limb_pos_angle, 99.91, 0.005);
        close(s.defect_pos_angle, 279.91, 0.005);
        close_val(s.k, 0.9012, 5e-5);
        close(s.q, 1.06 / 3600.0, 0.005 / 3600.0);
        close(s.d, 10.75 / 3600.0, 0.005 / 3600.0);
    }

    // The aberration corrections of step 14 must both be computed from the
    // uncorrected longitude. The old code used the corrected λ for β.
    #[test]
    fn aberration_uses_uncorrected_longitude() {
        let s = book_steps();
        let k = ABERRATION.to_radians();
        let l0 = 46.843_861_f64.to_radians();
        let expected = s.beta + k * (l0 - s.lambda).sin() * s.beta.sin();
        assert_eq!(s.beta_ab, expected);
    }

    // The Earth's heliocentric latitude of step 2 is small and negative;
    // `planet::heliocent_coords` returns it signed, so no re-wrapping is
    // needed before the steps that print `b0`.
    #[test]
    fn earth_heliocent_latitude_is_signed() {
        let (_, b, _) = planet::heliocent_coords(&planet::Planet::Earth, JDE);
        close(b, -0.000_167, 1e-5);
        assert!(b < 0.0);
    }
}

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

//! The ring system of Saturn
//!
//! Implements Meeus, *Astronomical Algorithms*, 2nd ed., chapter 45. For a
//! given instant, [`ring_elements`] describes how the ring (treated as one
//! flat system lying in Saturn's equatorial plane) appears from the Earth:
//!
//! * `B`, the Saturnicentric latitude of the Earth above the ring plane
//!   (positive when the northern face of the ring is visible),
//! * `B′`, the same for the Sun (positive when the northern face is lit),
//! * `ΔU`, the difference of the Saturnicentric longitudes of the Sun and
//!   the Earth in the ring plane (used for Saturn's magnitude, ch. 41),
//! * `P`, the position angle of the ring's northern semiminor axis, which
//!   is also that of Saturn's north pole,
//! * `a` and `b`, the major and minor axes of the outer edge of the outer
//!   ring. Multiply both by [`inn_edge_outer_ring`] and friends to get the
//!   other edges.
//!
//! [`elements`] is the older entry point that lets the caller supply the
//! nutation in longitude and the true obliquity.
//!
//! # Example
//!
//! Meeus' example 45.a, 1992 December 16 at 0h UT (JDE 2448972.50068):
//!
//! ```
//! use astro::planet::saturn::ring;
//!
//! let ring = ring::ring_elements(2448972.50068);
//!
//! assert!((ring.earth_lat.to_degrees() - 16.442).abs() < 0.0005);
//! assert!((ring.sun_lat.to_degrees() - 14.679).abs() < 0.0005);
//! assert!((ring.long_diff.to_degrees() - 4.198).abs() < 0.0005);
//! assert!((ring.axis_pos_angle.to_degrees() - 6.741).abs() < 0.0005);
//! // the axes of the outer edge, 35.87" and 10.15"
//! assert!((ring.major_axis.to_degrees() * 3600.0 - 35.87).abs() < 0.005);
//! assert!((ring.minor_axis.to_degrees() * 3600.0 - 10.15).abs() < 0.005);
//!
//! // the inner edge of the inner ring (ring B)
//! let (a, b) = ring::inn_edge_inner_ring(ring.major_axis, ring.minor_axis);
//! assert!((a.to_degrees() * 3600.0 - 23.85).abs() < 0.01);
//! # let _ = b;
//! ```

use core::f64::consts::{PI, TAU};

use crate::angle;
use crate::coords;
use crate::ecliptic;
use crate::nutation;
use crate::planet;
use crate::time;

/// Light-time for a distance of one AU, in days (Meeus, eq. 33.3).
const LIGHT_TIME_PER_AU: f64 = 0.005_775_518_3;

/// Major axis of the outer edge of the outer ring at one AU, in
/// arcseconds.
const MAJOR_AXIS_AT_1_AU: f64 = 375.35;

/// The constant of annual aberration (about 20.49″), in degrees, as used
/// by Meeus for the aberration of the planet in chapters 42, 43 and 45.
const ABERRATION: f64 = 0.005_693;

/// Computes the inclination of the plane of Saturn's ring
///
/// Meeus, eq. 45.1, from G. Dourneau's elements. The inclination is
/// referred to the ecliptic and mean equinox of the date.
///
/// # Returns
///
/// * `inc`: Inclination of the ring plane to the ecliptic *| in radians*
///
/// # Arguments
///
/// * `JC`: Julian (Ephemeris) centuries since J2000.0, see
///   [`time::julian_cent`]
#[inline]
pub fn inc(JC: f64) -> f64 {
    (28.075216 - JC * (0.012998 + JC * 0.000004)).to_radians()
}

/// Computes the longitude of the ascending node of the plane of Saturn's
/// ring
///
/// Meeus, eq. 45.1, from G. Dourneau's elements. The node is on the
/// ecliptic and referred to the mean equinox of the date.
///
/// # Returns
///
/// * `ascend_node`: Longitude of the ascending node of the ring plane
///   *| in radians*
///
/// # Arguments
///
/// * `JC`: Julian (Ephemeris) centuries since J2000.0, see
///   [`time::julian_cent`]
#[inline]
pub fn ascend_node(JC: f64) -> f64 {
    (169.50847 + JC * (1.394681 + JC * 0.000412)).to_radians()
}

/// Holds the elements for the ring system of Saturn
///
/// Returned by [`elements`]; [`RingElements`] holds the same values and a
/// few more. All angles are in radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Elements {
    /// Saturnicentric latitude of the Earth, referred to the plane
    /// of the ring, `[-π/2, π/2]`
    pub B: f64,
    /// Saturnicentric latitude of the Sun, referred to the plane of
    /// the ring, `[-π/2, π/2]`
    pub B1: f64,
    /// Geocentric position angle of the northern semiminor axis of
    /// the apparent ellipse of the ring (and of Saturn's north pole),
    /// `[0, 2π)`
    pub P: f64,
    /// Difference between Saturnicentric longitudes of the Sun and
    /// the Earth, measured in the plane of the ring, `[0, π]`
    pub deltaU: f64,
    /// Major axis of the outer edge of the outer ring
    pub a: f64,
    /// Minor axis of the outer edge of the outer ring
    pub b: f64,
}

/// The appearance of Saturn's ring at one instant
///
/// Returned by [`ring_elements`]. All angles are in radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RingElements {
    /// Saturnicentric latitude `B` of the Earth referred to the ring
    /// plane, positive towards the north, in `[-π/2, π/2]`. When positive,
    /// the northern face of the ring is the visible one.
    pub earth_lat: f64,
    /// Saturnicentric latitude `B′` of the Sun referred to the ring plane,
    /// in `[-π/2, π/2]`. When positive, the northern face is the lit one.
    pub sun_lat: f64,
    /// Difference `ΔU` between the Saturnicentric longitudes of the Sun
    /// and the Earth, measured in the ring plane, in `[0, π]` (in practice
    /// it never exceeds about 7°)
    pub long_diff: f64,
    /// Position angle `P` of the northern semiminor axis of the ring's
    /// apparent ellipse, measured from the north towards the east, in
    /// `[0, 2π)`. Saturn's north pole has the same position angle.
    pub axis_pos_angle: f64,
    /// Major axis `a` of the outer edge of the outer ring (the full width,
    /// not the semi-axis)
    pub major_axis: f64,
    /// Minor axis `b = a sin|B|` of the outer edge of the outer ring
    pub minor_axis: f64,
    /// Distance from the Earth to Saturn, `Δ`, at the time the light left
    /// Saturn *| in AU*
    pub earth_dist: f64,
    /// Distance from the Sun to Saturn, `r`, at the time the light left
    /// Saturn *| in AU*
    pub sun_dist: f64,
    /// Light-time from Saturn to the Earth, `τ` *| in days*
    pub light_time: f64,
}

/// Computes the elements of Saturn's ring
///
/// Meeus, ch. 45. The Earth and Saturn come from the crate's VSOP87
/// series, Saturn being taken at the time its light left it (light-time
/// found by iteration). The aberration of the Sun as seen from Saturn is
/// included in `B′`, and Saturn's aberration and the nutation in `P`; the
/// nutation and the obliquity are computed internally with
/// [`nutation::nutation`] and [`ecliptic::mn_oblq_IAU`].
///
/// The result agrees with Meeus' example 45.a to the printed precision
/// (0.001°, 0.01"). The ring plane's elements are given to about ±0.004°,
/// and the method is best within a few centuries of 2000.
///
/// The book asks for the Earth's position in the FK5 system; the FK5
/// correction (under 0.1") is not applied, as in the book's example.
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day (TT)
pub fn ring_elements(jde: f64) -> RingElements {
    let (nut_in_long, nut_in_oblq) = nutation::nutation(jde);
    let s = compute(
        jde,
        planet::heliocent_coords(&planet::Planet::Earth, jde),
        |t| planet::heliocent_coords(&planet::Planet::Saturn, t),
        nut_in_long,
        ecliptic::mn_oblq_IAU(jde) + nut_in_oblq,
    );

    RingElements {
        earth_lat: s.B,
        sun_lat: s.B1,
        long_diff: s.delta_U,
        axis_pos_angle: s.p,
        major_axis: s.a,
        minor_axis: s.b_axis,
        earth_dist: s.delta,
        sun_dist: s.r,
        light_time: s.tau,
    }
}

/// Computes the elements for the ring system of Saturn
///
/// The same algorithm as [`ring_elements`] (Meeus, ch. 45), with the
/// nutation in longitude and the true obliquity supplied by the caller.
///
/// # Returns
///
/// `elements`: Elements for the ring system. *All angles are in radians*
///
/// # Arguments
///
/// * `JD`         : Julian Ephemeris Day (TT)
/// * `nut_in_long`: Nutation in longitude on `JD` *| in radians*
/// * `true_oblq`  : True obliquity of the ecliptic on `JD` *| in radians*
pub fn elements(JD: f64, nut_in_long: f64, true_oblq: f64) -> Elements {
    let s = compute(
        JD,
        planet::heliocent_coords(&planet::Planet::Earth, JD),
        |t| planet::heliocent_coords(&planet::Planet::Saturn, t),
        nut_in_long,
        true_oblq,
    );

    Elements {
        B: s.B,
        B1: s.B1,
        P: s.p,
        deltaU: s.delta_U,
        a: s.a,
        b: s.b_axis,
    }
}

/// Factor for the inner edge of the outer ring (ring A).
const INN_EDGE_OUTER_RING: f64 = 0.8801;
/// Factor for the outer edge of the inner ring (ring B).
const OUT_EDGE_INNER_RING: f64 = 0.8599;
/// Factor for the inner edge of the inner ring (ring B).
const INN_EDGE_INNER_RING: f64 = 0.6650;
/// Factor for the inner edge of the dusky ring (ring C).
const INN_EDGE_DUSK_RING: f64 = 0.5486;

/// Computes the axes of the inner edge of the outer ring (ring A)
///
/// # Returns
///
/// `(a, b)`: Major and minor axes of the inner edge of the outer ring, in
/// the units of the arguments
///
/// # Arguments
///
/// * `a`: Major axis of the outer edge of the outer ring, see
///   [`RingElements::major_axis`]
/// * `b`: Minor axis of the outer edge of the outer ring, see
///   [`RingElements::minor_axis`]
#[inline(always)]
pub fn inn_edge_outer_ring(a: f64, b: f64) -> (f64, f64) {
    (a * INN_EDGE_OUTER_RING, b * INN_EDGE_OUTER_RING)
}

/// Computes the axes of the outer edge of the inner ring (ring B)
///
/// # Returns
///
/// `(a, b)`: Major and minor axes of the outer edge of the inner ring, in
/// the units of the arguments
///
/// # Arguments
///
/// * `a`: Major axis of the outer edge of the outer ring, see
///   [`RingElements::major_axis`]
/// * `b`: Minor axis of the outer edge of the outer ring, see
///   [`RingElements::minor_axis`]
#[inline(always)]
pub fn out_edge_inner_ring(a: f64, b: f64) -> (f64, f64) {
    (a * OUT_EDGE_INNER_RING, b * OUT_EDGE_INNER_RING)
}

/// Computes the axes of the outer edge of the inner ring (ring B)
///
/// Misspelt name kept for compatibility; use [`out_edge_inner_ring`].
#[deprecated(since = "3.0.0", note = "misspelt; use `out_edge_inner_ring`")]
#[inline(always)]
pub fn out_edge_inner_ing(a: f64, b: f64) -> (f64, f64) {
    out_edge_inner_ring(a, b)
}

/// Computes the axes of the inner edge of the inner ring (ring B)
///
/// # Returns
///
/// `(a, b)`: Major and minor axes of the inner edge of the inner ring, in
/// the units of the arguments
///
/// # Arguments
///
/// * `a`: Major axis of the outer edge of the outer ring, see
///   [`RingElements::major_axis`]
/// * `b`: Minor axis of the outer edge of the outer ring, see
///   [`RingElements::minor_axis`]
#[inline(always)]
pub fn inn_edge_inner_ring(a: f64, b: f64) -> (f64, f64) {
    (a * INN_EDGE_INNER_RING, b * INN_EDGE_INNER_RING)
}

/// Computes the axes of the inner edge of the dusky ring (ring C)
///
/// # Returns
///
/// `(a, b)`: Major and minor axes of the inner edge of the dusky ring, in
/// the units of the arguments
///
/// # Arguments
///
/// * `a`: Major axis of the outer edge of the outer ring, see
///   [`RingElements::major_axis`]
/// * `b`: Minor axis of the outer edge of the outer ring, see
///   [`RingElements::minor_axis`]
#[inline(always)]
pub fn inn_edge_dusk_ring(a: f64, b: f64) -> (f64, f64) {
    (a * INN_EDGE_DUSK_RING, b * INN_EDGE_DUSK_RING)
}

/// Every intermediate quantity of Meeus' procedure in ch. 45, named after
/// the step that produces it. Angles in radians, distances in AU, times in
/// days. Only some of them make up the public results; the others are
/// kept so that the unit tests can check every step against the book.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(not(test), allow(dead_code))]
struct Steps {
    // step 1
    inc: f64,
    node: f64,
    // steps 3-4
    l: f64,
    b: f64,
    r: f64,
    delta: f64,
    tau: f64,
    // step 5
    lambda: f64,
    beta: f64,
    // step 6
    B: f64,
    a: f64,
    b_axis: f64,
    // step 7
    N: f64,
    l1: f64,
    b1: f64,
    // step 8
    B1: f64,
    // step 9
    U1: f64,
    U2: f64,
    delta_U: f64,
    // step 11
    lambda0: f64,
    beta0: f64,
    // step 12
    lambda_ab: f64,
    beta_ab: f64,
    // step 13
    lambda0_app: f64,
    lambda_app: f64,
    // step 14
    asc0: f64,
    dec0: f64,
    asc: f64,
    dec: f64,
    // step 15
    p: f64,
}

/// Runs Meeus' ch. 45 procedure.
///
/// * `earth`: heliocentric `(l0, b0, R)` of the Earth at `jde`
/// * `saturn_at`: heliocentric `(l, b, r)` of Saturn at a given JDE
fn compute(
    jde: f64,
    earth: (f64, f64, f64),
    saturn_at: impl Fn(f64) -> (f64, f64, f64),
    nut_in_long: f64,
    tru_oblq: f64,
) -> Steps {
    // Step 1
    let T = time::julian_cent(jde);
    let inc = inc(T);
    let node = ascend_node(T);

    // Steps 3-4
    let (l0, b0, R) = earth;
    let geocentric = |tau: f64| {
        let (l, b, r) = saturn_at(jde - tau);
        let x = r * b.cos() * l.cos() - R * b0.cos() * l0.cos();
        let y = r * b.cos() * l.sin() - R * b0.cos() * l0.sin();
        let z = r * b.sin() - R * b0.sin();
        (l, b, r, x, y, z, (x * x + y * y + z * z).sqrt())
    };
    let mut tau = 0.0;
    let mut pos = geocentric(tau);
    for _ in 0..10 {
        let new_tau = LIGHT_TIME_PER_AU * pos.6;
        if (new_tau - tau).abs() < 1e-9 {
            break;
        }
        tau = new_tau;
        pos = geocentric(tau);
    }
    let (l, b, r, x, y, z, delta) = pos;

    // Step 5
    let lambda = y.atan2(x);
    let beta = z.atan2(x.hypot(y));

    // Step 6
    let B = ring_plane_lat(inc, node, lambda, beta);
    let a = (MAJOR_AXIS_AT_1_AU / 3600.0).to_radians() / delta;
    let b_axis = a * B.abs().sin();

    // Steps 7-8: the Sun seen from Saturn, corrected for aberration
    let N = (113.6655 + 0.8771 * T).to_radians();
    let l1 = l - (0.01759 / r).to_radians();
    let b1 = b - (0.000764 * (l - N).cos() / r).to_radians();
    let B1 = ring_plane_lat(inc, node, l1, b1);

    // Step 9
    let U1 = ring_plane_long(inc, node, l1, b1);
    let U2 = ring_plane_long(inc, node, lambda, beta);
    let delta_U = long_diff(U1, U2);

    // Step 11: the north pole of the ring plane
    let lambda0 = node - PI / 2.0;
    let beta0 = PI / 2.0 - inc;

    // Step 12: Saturn's aberration; both corrections use the uncorrected λ
    let k_ab = ABERRATION.to_radians();
    let lambda_ab = lambda + k_ab * (l0 - lambda).cos() / beta.cos();
    let beta_ab = beta + k_ab * (l0 - lambda).sin() * beta.sin();

    // Step 13
    let lambda0_app = lambda0 + nut_in_long;
    let lambda_app = lambda_ab + nut_in_long;

    // Step 14
    let asc0 = coords::asc_frm_ecl(lambda0_app, beta0, tru_oblq);
    let dec0 = coords::dec_frm_ecl(lambda0_app, beta0, tru_oblq);
    let asc = coords::asc_frm_ecl(lambda_app, beta_ab, tru_oblq);
    let dec = coords::dec_frm_ecl(lambda_app, beta_ab, tru_oblq);

    // Step 15
    let p = angle::limit_to_two_PI(
        (dec0.cos() * (asc0 - asc).sin())
            .atan2(dec0.sin() * dec.cos() - dec0.cos() * dec.sin() * (asc0 - asc).cos()),
    );

    Steps {
        inc,
        node,
        l,
        b,
        r,
        delta,
        tau,
        lambda,
        beta,
        B,
        a,
        b_axis,
        N,
        l1,
        b1,
        B1,
        U1,
        U2,
        delta_U,
        lambda0,
        beta0,
        lambda_ab,
        beta_ab,
        lambda0_app,
        lambda_app,
        asc0,
        dec0,
        asc,
        dec,
        p,
    }
}

/// Latitude above the ring plane of the direction with ecliptic
/// coordinates `(long, lat)` (Meeus, ch. 45, steps 6 and 8).
fn ring_plane_lat(inc: f64, node: f64, long: f64, lat: f64) -> f64 {
    (inc.sin() * lat.cos() * (long - node).sin() - inc.cos() * lat.sin())
        .clamp(-1.0, 1.0)
        .asin()
}

/// Longitude in the ring plane, from the ring's ascending node on the
/// ecliptic, of the direction `(long, lat)` (Meeus, ch. 45, step 9).
fn ring_plane_long(inc: f64, node: f64, long: f64, lat: f64) -> f64 {
    (inc.sin() * lat.sin() + inc.cos() * lat.cos() * (long - node).sin())
        .atan2(lat.cos() * (long - node).cos())
}

/// `|U1 - U2|` reduced to `[0, π]`: the angle between two longitudes.
///
/// The old code took `|U1 - U2|` of the raw atan2 values, which is close to
/// 2π whenever `U1` and `U2` straddle the ±π branch cut.
fn long_diff(U1: f64, U2: f64) -> f64 {
    ((U1 - U2 + PI).rem_euclid(TAU) - PI).abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dms(d: f64, m: f64, s: f64) -> f64 {
        let v = d.abs() + m / 60.0 + s / 3600.0;
        if d < 0.0 || d.is_sign_negative() {
            -v
        } else {
            v
        }
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

    const JDE: f64 = 2_448_972.500_68;

    /// Three-point interpolation (Meeus, eq. 3.3) in the table of Saturn's
    /// heliocentric coordinates at 1992 Dec 15.0, 16.0 and 17.0 TD printed
    /// in Meeus' example 45.a (from the full VSOP87).
    fn saturn_book(jde: f64) -> (f64, f64, f64) {
        let table = [
            (dms(319.0, 9.0, 44.23), dms(-1.0, 4.0, 26.52), 9.868_084_6),
            (dms(319.0, 11.0, 36.61), dms(-1.0, 4.0, 30.92), 9.867_869_0),
            (dms(319.0, 13.0, 28.99), dms(-1.0, 4.0, 35.31), 9.867_653_4),
        ];
        let n = jde - 2_448_972.5;
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

    /// Meeus' example 45.a from the book's own inputs: the Earth and
    /// Saturn from the full VSOP87 (steps 2-3) and its nutation and true
    /// obliquity (step 10).
    fn book_steps() -> Steps {
        compute(
            JDE,
            (
                84.285_703_f64.to_radians(),
                0.000_197_f64.to_radians(),
                0.984_123_16,
            ),
            saturn_book,
            (16.86 / 3600.0_f64).to_radians(),
            dms(23.0, 26.0, 22.96).to_radians(),
        )
    }

    // Meeus, example 45.a, the interpolated table values of step 3
    #[test]
    fn example_45a_book_table() {
        let (l, b, r) = saturn_book(2_448_972.500_68 - 0.051_98);
        close(l, 319.191_900, 2e-6);
        close(b, -1.075_192, 2e-6);
        close_val(r, 9.867_880_1, 1e-7);
    }

    // Meeus, example 45.a, steps 1 to 9.
    //
    // The geocentric angles carry ~2e-6° of play from the rounding of the
    // book's printed inputs (1e-7 AU at 10.5 AU is 5e-7°, 0.01" is 3e-6°).
    #[test]
    fn example_45a_steps_1_to_9() {
        let s = book_steps();
        close_val(time::julian_cent(JDE), -0.070_431_193, 1e-9);
        close(s.inc, 28.076_131, 5e-7);
        close(s.node, 169.410_243, 5e-7);
        close_val(s.tau, 0.060_44, 5e-6);
        close(s.l, 319.191_636, 2e-6);
        close(s.b, -1.075_183, 2e-6);
        close_val(s.r, 9.867_881_9, 1e-7);
        close_val(s.delta, 10.464_605_9, 1e-7);
        close(s.lambda, 314.777_850, 2e-6);
        close(s.beta, -1.013_885, 2e-6);
        close(s.B, 16.442, 5e-4);
        close(s.a, 35.87 / 3600.0, 0.005 / 3600.0);
        close(s.b_axis, 10.15 / 3600.0, 0.005 / 3600.0);
        close(s.N, 113.6037, 5e-5);
        close(s.l1, 319.189_853, 2e-6);
        close(s.b1, -1.075_113, 2e-6);
        close(s.B1, 14.679, 5e-4);
        close(s.U1, 153.2645, 5e-5);
        close(s.U2, 149.0663, 5e-5);
        close(s.delta_U, 4.198, 5e-4);
    }

    // Meeus, example 45.a, steps 11 to 15
    #[test]
    fn example_45a_steps_11_to_15() {
        let s = book_steps();
        close(s.lambda0, 79.410_243, 5e-7);
        close(s.beta0, 61.923_869, 5e-7);
        close(s.lambda_ab, 314.774_228, 2e-6);
        close(s.beta_ab, -1.013_963, 2e-6);
        close(s.lambda0_app, 79.414_926, 2e-6);
        close(s.lambda_app, 314.778_911, 2e-6);
        close(s.asc0, 40.363_65, 5e-6);
        close(s.dec0, 83.484_86, 5e-6);
        close(s.asc, 317.554_21, 5e-6);
        close(s.dec, -17.370_56, 5e-6);
        close(s.p, 6.741, 5e-4);
    }

    // ΔU must be the angle between U1 and U2 even when they straddle the
    // ±180° branch cut of atan2; the old code returned 359° for the first
    // case.
    #[test]
    fn long_diff_across_branch_cut() {
        let deg = |x: f64| x.to_radians();
        close(long_diff(deg(179.5), deg(-179.5)), 1.0, 1e-12);
        close(long_diff(deg(-179.5), deg(179.5)), 1.0, 1e-12);
        close(long_diff(deg(153.2645), deg(149.0663)), 4.1982, 1e-12);
        close(long_diff(deg(149.0663), deg(153.2645)), 4.1982, 1e-12);
        close(long_diff(0.0, 0.0), 0.0, 0.0);
        let u1 = ring_plane_long(28_f64.to_radians(), 0.0, deg(179.5), 0.0);
        let u2 = ring_plane_long(28_f64.to_radians(), 0.0, deg(-179.5), 0.0);
        assert!(u1 > 0.0 && u2 < 0.0);
        assert!(long_diff(u1, u2).to_degrees() < 1.2);
    }

    // Saturn's heliocentric latitude is negative in example 45.a;
    // `planet::heliocent_coords` returns it signed (no re-wrapping needed).
    #[test]
    fn saturn_heliocent_latitude_is_signed() {
        let (_, b, _) = planet::heliocent_coords(&planet::Planet::Saturn, JDE);
        assert!(b < 0.0 && b.to_degrees() > -1.1);
    }
}

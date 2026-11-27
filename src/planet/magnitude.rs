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

//! Phase angle and illuminated fraction of a planet's disk
//!
//! Additional algorithms of Meeus, *Astronomical Algorithms*, 2nd ed.,
//! chapter 41, complementing the functions in [`crate::planet`]
//! ([`illum_frac_frm_phase_angl`](super::illum_frac_frm_phase_angl),
//! [`illum_frac_frm_dist`](super::illum_frac_frm_dist),
//! [`phase_angl`](super::phase_angl),
//! [`apprnt_mag_muller`](super::apprnt_mag_muller),
//! [`apprnt_mag_84`](super::apprnt_mag_84)), in [`crate::planet::saturn`]
//! (Saturn's magnitude with the ring) and [`crate::pluto::apprnt_mag_84`]:
//!
//! * the phase angle from heliocentric coordinates (eq. 41.3) or from
//!   geocentric rectangular coordinates (eq. 41.4), and
//! * a low-accuracy illuminated fraction of Venus that needs only the date.
//!
//! # Example
//!
//! Venus on 1992 December 20.0 TD (Meeus, examples 41.a and 41.c):
//!
//! ```
//! use astro::planet::{self, magnitude, Planet};
//!
//! // Heliocentric coordinates of Venus (light-time corrected) and of the
//! // Earth, and the Earth-Venus distance, from Meeus' example 33.a.
//! let i = magnitude::phase_angl_frm_heliocent(
//!     26.10588_f64.to_radians(),
//!     (-2.62102_f64).to_radians(),
//!     0.724604,
//!     88.35704_f64.to_radians(),
//!     0.983824,
//!     0.910947,
//! );
//! let k = planet::illum_frac_frm_phase_angl(i);
//! assert!((k - 0.647).abs() < 5e-4);
//!
//! let m = planet::apprnt_mag_muller(&Planet::Venus, i, 0.910947, 0.724604).unwrap();
//! assert!((m - -3.8).abs() < 0.05);
//!
//! // The low-accuracy method needs only the date (example 41.b).
//! let approx = magnitude::venus_illum_frac_approx(2448976.5);
//! assert!((approx.illum_frac - 0.640).abs() < 5e-4);
//! ```

use crate::time;

/// Computes a planet's phase angle from heliocentric coordinates
///
/// `cos i = (R − R0 cos B cos(L − L0)) / Δ` (Meeus, eq. 41.3), for a
/// position obtained by the "first method" of Meeus' chapter 33. The
/// Earth's heliocentric latitude (always below 1.2″) is neglected.
///
/// # Returns
///
/// * `phase_angl`: Phase angle `i` (Sun-planet-Earth) *| in radians*, in
///   `[0, π]`
///
/// # Arguments
///
/// * `planet_long`      : Heliocentric longitude `L` of the planet
///   *| in radians* (light-time corrected, as for the position)
/// * `planet_lat`       : Heliocentric latitude `B` of the planet
///   *| in radians*
/// * `planet_rad_vec`   : Heliocentric radius vector `R` of the planet
///   *| in AU*
/// * `earth_long`       : Heliocentric longitude `L0` of the Earth
///   *| in radians*
/// * `earth_rad_vec`    : Heliocentric radius vector `R0` of the Earth
///   *| in AU*
/// * `planet_earth_dist`: Planet-Earth distance `Δ` *| in AU*
pub fn phase_angl_frm_heliocent(
    planet_long: f64,
    planet_lat: f64,
    planet_rad_vec: f64,
    earth_long: f64,
    earth_rad_vec: f64,
    planet_earth_dist: f64,
) -> f64 {
    let cos_i = (planet_rad_vec
        - earth_rad_vec * planet_lat.cos() * (planet_long - earth_long).cos())
        / planet_earth_dist;

    cos_i.clamp(-1.0, 1.0).acos()
}

/// Computes a planet's phase angle from its geocentric rectangular
/// ecliptic coordinates
///
/// `cos i = (x cos B cos L + y cos B sin L + z sin B) / Δ` with
/// `Δ = √(x² + y² + z²)` (Meeus, eq. 41.4).
///
/// # Returns
///
/// * `phase_angl`: Phase angle `i` (Sun-planet-Earth) *| in radians*, in
///   `[0, π]`
///
/// # Arguments
///
/// * `x`, `y`, `z`: Geocentric rectangular ecliptic coordinates of the
///   planet (Meeus, eq. 33.1) *| in AU*
/// * `planet_long`: Heliocentric longitude `L` of the planet *| in radians*
/// * `planet_lat` : Heliocentric latitude `B` of the planet *| in radians*
pub fn phase_angl_frm_geocent_rect(
    x: f64,
    y: f64,
    z: f64,
    planet_long: f64,
    planet_lat: f64,
) -> f64 {
    let dist = (x * x + y * y + z * z).sqrt();
    let (sin_l, cos_l) = planet_long.sin_cos();
    let (sin_b, cos_b) = planet_lat.sin_cos();
    let cos_i = (x * cos_b * cos_l + y * cos_b * sin_l + z * sin_b) / dist;

    cos_i.clamp(-1.0, 1.0).acos()
}

/// Low-accuracy illuminated fraction, distance and elongation of Venus
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VenusIllumApprox {
    /// Illuminated fraction `k` of the disk, in `[0, 1]`
    pub illum_frac: f64,
    /// Approximate Earth-Venus distance `Δ` *| in AU*
    pub dist: f64,
    /// Approximate elongation `ψ` of Venus from the Sun *| in radians*, in
    /// `[0, π]`
    pub elongation: f64,
}

/// Angles of the Venus approximation, in degrees: `(V, M, N, W)`.
fn venus_angles(T: f64) -> (f64, f64, f64, f64) {
    let V = 261.51 + 22518.443 * T;
    let M = 177.53 + 35999.050 * T;
    let N = 50.42 + 58517.811 * T;
    let W = V + 1.91 * M.to_radians().sin() + 0.78 * N.to_radians().sin();
    (V, M, N, W)
}

/// Computes an approximate illuminated fraction of the disk of Venus
///
/// A quick method of Meeus, chapter 41, that needs only the date: Venus and
/// the Earth are put on circular orbits (radii 0.72333 and 1 AU) and the
/// difference of their mean longitudes is corrected by the two largest
/// equation-of-centre terms. It also yields approximate values of the
/// distance and the elongation. For an accurate result compute the phase
/// angle from positions and use
/// [`illum_frac_frm_phase_angl`](super::illum_frac_frm_phase_angl); in
/// Meeus' example 41.b this method gives 0.640 where the accurate value is
/// 0.647.
///
/// # Returns
///
/// See [`VenusIllumApprox`].
///
/// # Arguments
///
/// * `JD`: Julian Ephemeris Day
pub fn venus_illum_frac_approx(JD: f64) -> VenusIllumApprox {
    let T = time::julian_cent(JD);
    let (_, _, _, W) = venus_angles(T);

    let dist = (1.52321 + 1.44666 * W.to_radians().cos()).sqrt();
    let a = 0.72333 + dist;
    // The rounded constants can push k a few 1e-7 above 1 near superior
    // conjunction; clamp to the physical range.
    let illum_frac = ((a * a - 1.0) / (2.89332 * dist)).clamp(0.0, 1.0);
    let cos_psi = (dist * dist + 0.4768) / (2.0 * dist);

    VenusIllumApprox {
        illum_frac,
        dist,
        elongation: cos_psi.clamp(-1.0, 1.0).acos(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planet::{self, Planet};

    fn assert_close(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() <= tol,
            "expected {expected} ± {tol}, got {actual}"
        );
    }

    // Meeus, example 41.a: cos i = 0.29312 from eq. 41.3, with L, B, R of
    // Venus (B) and L0, R0 of the Earth (A) of example 33.a.
    #[test]
    fn phase_angl_frm_heliocent_example_41a() {
        let i = phase_angl_frm_heliocent(
            26.10588_f64.to_radians(),
            (-2.62102_f64).to_radians(),
            0.724604,
            88.35704_f64.to_radians(),
            0.983824,
            0.910947,
        );
        assert_close(i.cos(), 0.29312, 5e-6);
        assert_close(planet::illum_frac_frm_phase_angl(i), 0.647, 5e-4);
        assert_close(i.to_degrees(), 72.96, 0.005); // example 41.c
    }

    // Meeus, example 41.a: the same from eq. 41.4 with x, y, z (C).
    #[test]
    fn phase_angl_frm_geocent_rect_example_41a() {
        let i = phase_angl_frm_geocent_rect(
            0.621794,
            -0.664905,
            -0.033138,
            26.10588_f64.to_radians(),
            (-2.62102_f64).to_radians(),
        );
        assert_close(i.cos(), 0.29312, 5e-6);
        assert_close(planet::illum_frac_frm_phase_angl(i), 0.647, 5e-4);
    }

    // The three ways of computing the phase angle agree for a consistent
    // geometry, built here from arbitrary heliocentric positions.
    #[test]
    fn phase_angl_methods_agree() {
        let cases = [
            (26.1_f64, -2.6_f64, 0.7246_f64, 88.36_f64, 0.9838_f64),
            (200.0, 1.5, 1.52, 10.0, 1.0167),
            (350.0, 0.5, 5.2, 5.0, 0.9833), // across 0°
            (100.0, 6.0, 0.39, 280.0, 1.0), // inferior planet beyond the Sun
        ];
        for &(l, b, r, l0, r0) in &cases {
            let (l, b, l0) = (l.to_radians(), b.to_radians(), l0.to_radians());
            let x = r * b.cos() * l.cos() - r0 * l0.cos();
            let y = r * b.cos() * l.sin() - r0 * l0.sin();
            let z = r * b.sin();
            let delta = (x * x + y * y + z * z).sqrt();

            let i1 = phase_angl_frm_heliocent(l, b, r, l0, r0, delta);
            let i2 = phase_angl_frm_geocent_rect(x, y, z, l, b);
            let i3 = planet::phase_angl(r, delta, r0);
            assert_close(i1, i3, 1e-12);
            assert_close(i2, i3, 1e-12);
            assert!((0.0..=std::f64::consts::PI).contains(&i1));
        }
    }

    // Meeus, example 41.b.
    #[test]
    fn venus_illum_frac_approx_example_41b() {
        let T = time::julian_cent(2448976.5);
        assert_close(T, -0.070321697, 1e-9);

        let (V, M, N, W) = venus_angles(T);
        assert_close(V, -1322.025, 5e-4);
        assert_close(M, -2353.984, 5e-4);
        assert_close(N, -4064.652, 5e-4);
        assert_close(W - V, 0.462 - 0.755, 1e-3);
        assert_close((W - V).abs(), (0.462_f64 - 0.755).abs(), 1e-3);
        assert_close(crate::angle::limit_to_360(W), 117.682, 5e-4);

        let v = venus_illum_frac_approx(2448976.5);
        // The book prints Δ = 0.922575; with W = 117.682° exactly the
        // formula gives 0.92258, so the printed value carries the rounding
        // of the intermediate angles.
        assert_close(v.dist, 0.922575, 1e-5);
        assert_close(v.illum_frac, 0.640, 5e-4);
        // Accurate elongation on that date (example 33.a): about 44.7°.
        assert_close(v.elongation.to_degrees(), 44.7, 1.0);
    }

    // Consistency of the approximation with its own geometry: k from the
    // phase angle of the triangle Sun (0) - Venus (0.72333) - Earth (1).
    #[test]
    fn venus_illum_frac_approx_geometry() {
        let mut jd = 2415020.0; // 1900
        while jd < 2488070.0 {
            let v = venus_illum_frac_approx(jd);
            assert!(v.dist > 0.27 && v.dist < 1.73, "{jd}: {}", v.dist);
            assert!(
                (0.0..=1.0).contains(&v.illum_frac),
                "{jd}: {}",
                v.illum_frac
            );
            // Greatest elongation of Venus is about 47°.
            assert!(v.elongation.to_degrees() < 47.5, "{jd}: {}", v.elongation);
            let i = planet::phase_angl(0.72333, v.dist, 1.0);
            assert_close(v.illum_frac, planet::illum_frac_frm_phase_angl(i), 1e-4);
            jd += 37.7;
        }
    }

    // Example 41.c through the generic magnitude function (phase angle in
    // radians).
    #[test]
    fn venus_magnitude_example_41c() {
        let i = 72.96_f64.to_radians();
        let m = planet::apprnt_mag_muller(&Planet::Venus, i, 0.910947, 0.724604).unwrap();
        assert_close(m, -3.8, 0.05);
    }
}

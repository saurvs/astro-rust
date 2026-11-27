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

//! Miscellaneous routines
//!
//! The parallactic angle and related angles between great circles on the
//! sky (Meeus, *Astronomical Algorithms*, 2nd ed., chapter 14):
//!
//! * the parallactic angle: [`parllc_angl`], [`parllc_angl_on_hz`];
//! * the ecliptic and the horizon: [`eclip_hz_points`] here, and
//!   [`ecliptic::eclip_points_on_hz`](crate::ecliptic::eclip_points_on_hz),
//!   [`ecliptic::angl_betwn_eclip_and_hz`](crate::ecliptic::angl_betwn_eclip_and_hz);
//! * the ecliptic and the equator: [`angl_betwn_north_celes_and_eclip_poles`],
//!   [`angl_betwn_eclip_and_east_west`];
//! * the diurnal path and the horizon:
//!   [`planet::earth::angl_betwn_diurnal_path_and_hz`](crate::planet::earth::angl_betwn_diurnal_path_and_hz).
//!
//! All angles are in radians. Hour angles are measured westwards from the
//! meridian, latitudes are positive north.
//!
//! # Example
//!
//! ```
//! use astro::misc;
//!
//! let lat = 51_f64.to_radians();
//! let dec = 20_f64.to_radians();
//!
//! // On the meridian, north is straight up from the zenith's direction.
//! assert_eq!(misc::parllc_angl(lat, 0.0, dec), 0.0);
//!
//! // Two hours before and after the meridian the angle has opposite signs.
//! let before = misc::parllc_angl(lat, -30_f64.to_radians(), dec);
//! let after = misc::parllc_angl(lat, 30_f64.to_radians(), dec);
//! assert!(before < 0.0 && after > 0.0);
//! assert!((before + after).abs() < 1e-15);
//! ```

use crate::coords;

/// Computes the parallactic angle of a celestial body
///
/// The parallactic angle `q` is the angle, at the body, between the
/// direction of the north celestial pole and the direction of the zenith
/// (Meeus, eq. 14.1). It is the position angle of the zenith measured from
/// the north through the east, so it is negative before the body crosses
/// the meridian, zero on the meridian south of the zenith, and positive
/// after. Use it, for example, to orient the image of the Sun or the Moon
/// relative to the vertical.
///
/// # Returns
///
/// * `parallac_angle`: Parallactic angle, in `(−π, π]`
///   *| in radians*. It is undefined at the zenith (where `0` is
///   returned) and jumps from `−π/2` to `+π/2` when a body passes exactly
///   through the zenith.
///
/// # Arguments
///
/// * `observer_lat`: The observer's geographical latitude *| in radians*
/// * `hour_angle`  : Local hour angle *| in radians*
/// * `dec`         : Declination of the celestial body *| in radians*
#[inline]
pub fn parllc_angl(observer_lat: f64, hour_angle: f64, dec: f64) -> f64 {
    // Eq. 14.1 with numerator and denominator multiplied by cos φ ≥ 0,
    // which keeps it finite at the poles.
    let (sin_h, cos_h) = hour_angle.sin_cos();
    let (sin_d, cos_d) = dec.sin_cos();
    let (sin_p, cos_p) = observer_lat.sin_cos();

    (sin_h * cos_p).atan2(sin_p * cos_d - cos_p * sin_d * cos_h)
}

/// Computes the parallactic angle of a celestial body on the
/// horizon
///
/// When the body rises or sets, `cos q = sin φ / cos δ` (Meeus,
/// chapter 14), and the hour angle is not needed.
///
/// # Returns
///
/// * `parallac_angle`: Magnitude of the parallactic angle at rising or
///   setting, in `[0, π]` *| in radians*. The angle is positive at setting
///   and negative at rising, so negate the result for a rising body. The
///   result is `NaN` if the body never reaches the horizon (circumpolar or
///   never rising, `|sin φ| > cos δ`).
///
/// # Arguments
///
/// * `observer_lat`: The observer's geographical latitude *| in radians*
/// * `dec`         : Declination of the celestial body *| in radians*
#[inline]
pub fn parllc_angl_on_hz(observer_lat: f64, dec: f64) -> f64 {
    (observer_lat.sin() / dec.cos()).acos()
}

/// Computes the angle between the directions of the north celestial pole
/// and the north ecliptic pole, at a point on the sky
///
/// Meeus, chapter 14 ("Ecliptic and Equator"). The angle is the position
/// angle of the north pole of the ecliptic measured from the north
/// celestial pole through the east: `tan q = cos λ tan ε /
/// (sin β sin λ tan ε − cos β)`, with the quadrant chosen so that
/// `q = 0` at the summer solstice point (λ = 90°, β = 0).
///
/// # Returns
///
/// * `angl`: The angle, in `(−π, π]` *| in radians*
///
/// # Arguments
///
/// * `ecl_long`: Ecliptic longitude of the point *| in radians*
/// * `ecl_lat`: Ecliptic latitude of the point *| in radians*
/// * `oblq_eclip`: Obliquity of the ecliptic *| in radians*
pub fn angl_betwn_north_celes_and_eclip_poles(ecl_long: f64, ecl_lat: f64, oblq_eclip: f64) -> f64 {
    // Book formula with numerator and denominator multiplied by
    // −cos ε (cos ε > 0), which fixes the quadrant and avoids tan ε.
    let (sin_l, cos_l) = ecl_long.sin_cos();
    let (sin_b, cos_b) = ecl_lat.sin_cos();
    let (sin_e, cos_e) = oblq_eclip.sin_cos();

    let q = (-cos_l * sin_e).atan2(cos_b * cos_e - sin_b * sin_l * sin_e);
    if q == -std::f64::consts::PI {
        std::f64::consts::PI
    } else {
        q
    }
}

/// Computes the angle between the ecliptic and the east-west direction at
/// a point of the ecliptic
///
/// Meeus, chapter 14: `tan q₀ = −cos λ tan ε`, the value of
/// [`angl_betwn_north_celes_and_eclip_poles`] for β = 0. Along the
/// ecliptic, the direction of increasing longitude points `−q₀` north of
/// due east: near the vernal equinox (`q₀ = −ε`) the ecliptic climbs
/// northwards towards the east, near the autumnal equinox (`q₀ = +ε`) it
/// descends. Useful, for example, to draw the path of the Moon through
/// the Earth's shadow in a lunar eclipse diagram.
///
/// # Returns
///
/// * `angl`: The angle `q₀`, in `[−ε, ε]` *| in radians*
///
/// # Arguments
///
/// * `ecl_long`: Ecliptic longitude of the point *| in radians*
/// * `oblq_eclip`: Obliquity of the ecliptic *| in radians*
pub fn angl_betwn_eclip_and_east_west(ecl_long: f64, oblq_eclip: f64) -> f64 {
    let (sin_e, cos_e) = oblq_eclip.sin_cos();
    (-ecl_long.cos() * sin_e).atan2(cos_e)
}

/// The two points of the ecliptic that are on the horizon
///
/// See [`eclip_hz_points`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EclipHzPoints {
    /// Ecliptic longitude of the point on the eastern horizon (the rising
    /// point, or ascendant), in `[0, 2π)` *| in radians*
    pub east_long: f64,
    /// Ecliptic longitude of the point on the western horizon (the setting
    /// point, or descendant), in `[0, 2π)` *| in radians*
    pub west_long: f64,
}

/// Computes the ecliptic longitudes of the two points of the ecliptic on
/// the horizon, and tells which is which
///
/// Meeus, eq. 14.2. Same as
/// [`ecliptic::eclip_points_on_hz`](crate::ecliptic::eclip_points_on_hz),
/// with named results. The two points are 180° apart. The eastern and
/// western points are identified correctly at all latitudes, including
/// inside the polar circles, where the usual `atan2` form of eq. 14.2
/// sometimes returns the eastern point in place of the western one. The
/// points are undefined at the instants when the ecliptic coincides with
/// the horizon (possible only where `|φ| = 90° − ε`) and east and west are
/// undefined at the poles.
///
/// # Arguments
///
/// * `oblq_eclip`  : Obliquity of the ecliptic *| in radians*
/// * `observer_lat`: The observer's geographical latitude *| in radians*
/// * `loc_sidreal` : Local sidereal time *| in radians*
///
/// # Example
///
/// Meeus, example 14.a:
///
/// ```
/// use astro::misc;
///
/// let p = misc::eclip_hz_points(
///     23.44_f64.to_radians(),
///     51_f64.to_radians(),
///     75_f64.to_radians(),
/// );
/// assert!((p.east_long.to_degrees() - 169.36).abs() < 0.01); // 169°21′
/// assert!((p.west_long.to_degrees() - 349.36).abs() < 0.01); // 349°21′
/// ```
pub fn eclip_hz_points(oblq_eclip: f64, observer_lat: f64, loc_sidreal: f64) -> EclipHzPoints {
    let west_long = descendant(oblq_eclip, observer_lat, loc_sidreal);
    EclipHzPoints {
        east_long: coords::wrap_two_pi(west_long + std::f64::consts::PI),
        west_long,
    }
}

/// Ecliptic longitude of the point of the ecliptic on the western horizon,
/// in `[0, 2π)`.
///
/// Eq. 14.2 gives `tan λ`, whose two solutions are the two points. With the
/// ecliptic point `(cos λ, sin λ cos ε, sin λ sin ε)` and the zenith in
/// equatorial coordinates, the horizon condition reads
/// `a cos λ + b sin λ = 0` with `a = cos θ cos φ` and
/// `b = cos ε sin θ cos φ + sin ε sin φ` (the book's form multiplied by
/// cos φ ≥ 0, which keeps it finite at the poles). The solution
/// `(cos λ, sin λ) ∝ (b, −a)` lies on the western side exactly when
/// `k = cos ε cos φ + sin ε sin φ sin θ > 0`. That holds everywhere outside
/// the polar circles, but not always inside them, so the sign of `k` picks
/// the solution.
pub(crate) fn descendant(oblq_eclip: f64, observer_lat: f64, loc_sidreal: f64) -> f64 {
    let (sin_t, cos_t) = loc_sidreal.sin_cos();
    let (sin_e, cos_e) = oblq_eclip.sin_cos();
    let (sin_p, cos_p) = observer_lat.sin_cos();

    let a = cos_t * cos_p;
    let b = cos_e * sin_t * cos_p + sin_e * sin_p;
    let k = cos_e * cos_p + sin_e * sin_p * sin_t;
    let sign = if k < 0.0 { -1.0 } else { 1.0 };

    coords::wrap_two_pi((-a * sign).atan2(b * sign))
}

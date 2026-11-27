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

//! Transform between coordinate systems
//!
//! Conversions between the equatorial, ecliptic, local horizontal and
//! galactic systems (Meeus, *Astronomical Algorithms*, 2nd ed., chapter 13),
//! and the hour angle of a body from the sidereal time.
//!
//! # Conventions
//!
//! All angles are in **radians**.
//!
//! * **Geographic longitude is positive WEST** of Greenwich and negative
//!   east, as in Meeus. With that convention the local hour angle is
//!   `H = θ₀ − L − α`. If your longitudes are positive east (the IAU and
//!   ISO 6709 convention), use [`hr_angl_frm_observer_east_long`] or convert
//!   with [`west_long_frm_east_long`].
//! * **Azimuth is measured from the SOUTH, increasing towards the west**
//!   (S = 0, W = π/2, N = π, E = 3π/2), as in Meeus. Navigators and most
//!   other software measure it from the north, increasing towards the east;
//!   the two differ by exactly π. Convert with [`az_north_frm_south`] and
//!   [`az_south_frm_north`].
//! * The hour angle is measured westwards from the observer's meridian
//!   (positive after meridian transit).
//!
//! Results are normalised as follows: right ascension and ecliptic or
//! galactic longitude to `[0, 2π)`; hour angle and azimuth to `(−π, π]`;
//! declination, altitude and latitudes to `[−π/2, π/2]`.
//!
//! The transformations are rigid rotations of the celestial sphere, so they
//! are computed from Cartesian components with `atan2`, which keeps them
//! accurate near the poles of either system (where the textbook forms with
//! `tan` and `asin` lose precision or divide by zero). Where a longitude-like
//! angle is undefined (a point exactly at a pole), `0` is returned.
//!
//! # Example
//!
//! Venus seen from the U.S. Naval Observatory, Washington (Meeus,
//! example 13.b):
//!
//! ```
//! use astro::{angle, coords};
//!
//! // Apparent sidereal time at Greenwich, apparent position of Venus.
//! let greenw_sidr = angle::deg_frm_hms(8, 34, 56.853).to_radians();
//! let asc = angle::deg_frm_hms(23, 9, 16.641).to_radians();
//! let dec = angle::deg_frm_dms(-6, 43, 11.61).to_radians();
//!
//! // Longitude positive WEST: Washington is at +77°03′56″.
//! let washington = coords::GeographPoint {
//!     long: angle::deg_frm_hms(5, 8, 15.7).to_radians(),
//!     lat: angle::deg_frm_dms(38, 55, 17.0).to_radians(),
//! };
//!
//! let hour_angle = coords::hr_angl_frm_observer_long(greenw_sidr, washington.long, asc);
//! let hz = coords::hz_frm_hr_angl_and_dec(hour_angle, dec, washington.lat);
//!
//! // Azimuth from the south (towards the west), as in Meeus.
//! assert!((hz.az.to_degrees() - 68.0337).abs() < 1e-4);
//! assert!((hz.alt.to_degrees() - 15.1249).abs() < 1e-4);
//!
//! // The same direction on a compass (from the north, towards the east).
//! let compass = coords::az_north_frm_south(hz.az).to_degrees();
//! assert!((compass - 248.0337).abs() < 1e-4);
//! ```

use crate::angle;
use std::f64::consts::{PI, TAU};

/// Reduces an angle to `[0, 2π)`.
#[inline]
pub(crate) fn wrap_two_pi(angl: f64) -> f64 {
    if (0.0..TAU).contains(&angl) {
        return angl;
    }
    let a = angl.rem_euclid(TAU);
    // rem_euclid can round a tiny negative angle up to exactly 2π.
    if a >= TAU { 0.0 } else { a }
}

/// Reduces an angle to `(−π, π]`.
#[inline]
pub(crate) fn wrap_pi(angl: f64) -> f64 {
    if angl > -PI && angl <= PI {
        return angl;
    }
    let a = angl.rem_euclid(TAU);
    if a > PI { a - TAU } else { a }
}

/// `atan2` with the result in `(−π, π]` (`atan2(−0, x < 0)` is `−π`).
#[inline]
fn atan2_pi(y: f64, x: f64) -> f64 {
    let a = y.atan2(x);
    if a == -PI { PI } else { a }
}

/// Mean obliquity of the ecliptic at the standard epoch J2000.0,
/// 23°26′21.448″ *| in radians*
///
/// Use it with [`ecl_long_frm_eq`] and friends for coordinates referred to
/// the mean equator and equinox of J2000.0 (Meeus, chapter 13).
pub const OBLQ_J2000: f64 = 23.439_291_111_111_11 * PI / 180.0;

/// Mean obliquity of the ecliptic at the standard epoch B1950.0,
/// 23.445 7889° *| in radians*
///
/// Use it for coordinates referred to the mean equator and equinox of
/// B1950.0 (Meeus, chapter 13).
pub const OBLQ_B1950: f64 = 23.445_788_9 * PI / 180.0;

/// A point on the geographical surface of the Earth
///
/// The longitude is **positive west** of Greenwich, negative east (Meeus'
/// convention). The latitude is the geographic (geodetic) latitude,
/// positive north.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeographPoint {
    /// Geographical longitude, positive west *| in radians*
    pub long: f64,
    /// Geographical latitude, positive north *| in radians*
    pub lat: f64,
}

impl GeographPoint {
    /// Returns the angle between two points as seen from the centre of a
    /// spherical Earth *| in radians*
    ///
    /// Multiply by a mean Earth radius to get a distance; see
    /// [`planet::earth::approx_geodesic_dist`](crate::planet::earth::approx_geodesic_dist).
    pub fn anglr_sepr(&self, other_point: &GeographPoint) -> f64 {
        angle::anglr_sepr(self.long, self.lat, other_point.long, other_point.lat)
    }
}

/// A point in the equatorial coordinate system
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EqPoint {
    /// Right ascension *| in radians*
    pub asc: f64,
    /// Declination *| in radians*
    pub dec: f64,
}

impl EqPoint {
    /// Returns the angular separation between two equatorial points
    /// *| in radians*
    pub fn anglr_sepr(&self, other_point: &EqPoint) -> f64 {
        angle::anglr_sepr(self.asc, self.dec, other_point.asc, other_point.dec)
    }

    /// Converts the point to ecliptic coordinates
    ///
    /// Same as [`ecl_long_frm_eq`] and [`ecl_lat_frm_eq`] together.
    ///
    /// # Arguments
    ///
    /// * `oblq_eclip`: Obliquity of the ecliptic for the equinox the point
    ///   is referred to (true obliquity for apparent places) *| in radians*
    pub fn to_ecl(&self, oblq_eclip: f64) -> EclPoint {
        let (x, y, z) = rotate_x(self.asc, self.dec, oblq_eclip);
        EclPoint {
            long: wrap_two_pi(y.atan2(x)),
            lat: z.atan2(x.hypot(y)),
        }
    }
}

/// A point in the ecliptic coordinate system
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EclPoint {
    /// Ecliptic longitude *| in radians*
    pub long: f64,
    /// Ecliptic latitude *| in radians*
    pub lat: f64,
}

impl EclPoint {
    /// Returns the angular separation between two ecliptic points
    /// *| in radians*
    pub fn anglr_sepr(&self, other_point: &EclPoint) -> f64 {
        angle::anglr_sepr(self.long, self.lat, other_point.long, other_point.lat)
    }

    /// Converts the point to equatorial coordinates
    ///
    /// Same as [`asc_frm_ecl`] and [`dec_frm_ecl`] together.
    ///
    /// # Arguments
    ///
    /// * `oblq_eclip`: Obliquity of the ecliptic for the equinox the point
    ///   is referred to (true obliquity for apparent places) *| in radians*
    pub fn to_eq(&self, oblq_eclip: f64) -> EqPoint {
        let (x, y, z) = rotate_x(self.long, self.lat, -oblq_eclip);
        EqPoint {
            asc: wrap_two_pi(y.atan2(x)),
            dec: z.atan2(x.hypot(y)),
        }
    }
}

/// A direction in the local horizontal coordinate system
///
/// The azimuth is measured from the **south**, increasing towards the west
/// (Meeus' convention); see [`az_north_frm_south`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HzPoint {
    /// Azimuth, from the south towards the west, in `(−π, π]`
    /// *| in radians*
    pub az: f64,
    /// Altitude above the horizon, in `[−π/2, π/2]` *| in radians*
    pub alt: f64,
}

/// A direction given by its local hour angle and declination
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HrAnglDec {
    /// Local hour angle, westwards from the meridian, in `(−π, π]`
    /// *| in radians*
    pub hour_angle: f64,
    /// Declination, in `[−π/2, π/2]` *| in radians*
    pub dec: f64,
}

/// Cartesian components of the unit vector (`long`, `lat`) after a rotation
/// of the frame by `angl` about the x axis (the axis through longitude 0).
///
/// With `angl = ε` this turns equatorial into ecliptic components, with
/// `angl = −ε` the reverse.
#[inline]
fn rotate_x(long: f64, lat: f64, angl: f64) -> (f64, f64, f64) {
    let (sin_long, cos_long) = long.sin_cos();
    let (sin_lat, cos_lat) = lat.sin_cos();
    let (sin_a, cos_a) = angl.sin_cos();

    (
        cos_lat * cos_long,
        cos_lat * sin_long * cos_a + sin_lat * sin_a,
        sin_lat * cos_a - cos_lat * sin_long * sin_a,
    )
}

/// Converts an east-positive longitude to the west-positive convention used
/// throughout this crate *| in radians*
///
/// The two conventions differ only in sign. The result is not reduced to
/// any range.
#[inline]
pub fn west_long_frm_east_long(east_long: f64) -> f64 {
    -east_long
}

/// Converts a west-positive longitude (this crate's convention) to the
/// east-positive one *| in radians*
///
/// The two conventions differ only in sign. The result is not reduced to
/// any range.
#[inline]
pub fn east_long_frm_west_long(west_long: f64) -> f64 {
    -west_long
}

/// Computes the local sidereal time from the sidereal time at Greenwich
///
/// `θ = θ₀ − L`, with the longitude `L` positive **west** (Meeus,
/// chapter 13).
///
/// # Returns
///
/// * `local_sidreal`: Local sidereal time, in `[0, 2π)` *| in radians*
///
/// # Arguments
///
/// * `green_sidreal`: Sidereal time at Greenwich *| in radians*
/// * `observer_long`: Observer's geographical longitude, **positive west**
///   *| in radians*
#[inline]
pub fn loc_sidr_frm_green_sidr(green_sidreal: f64, observer_long: f64) -> f64 {
    wrap_two_pi(green_sidreal - observer_long)
}

/// Computes the hour angle from geographical longitude and Greenwich
/// sidereal time
///
/// `H = θ₀ − L − α` (Meeus, chapter 13). **The longitude `L` must be
/// positive WEST** of Greenwich and negative east, like
/// [`GeographPoint::long`] everywhere in this crate. Passing an
/// east-positive longitude here gives a wrong hour angle; for east-positive
/// longitudes (`H = θ₀ + L_east − α`) use [`hr_angl_frm_observer_east_long`].
///
/// If the right ascension is apparent (affected by nutation), the sidereal
/// time must be the apparent one as well.
///
/// # Returns
///
/// * `hour_angle`: Local hour angle, measured westwards from the meridian,
///   in `(−π, π]` *| in radians*
///
/// # Arguments
///
/// * `green_sidreal`: Sidereal time at Greenwich *| in radians*
/// * `observer_long`: Observer's geographical longitude, **positive west**
///   *| in radians*
/// * `asc`: Right ascension *| in radians*
///
/// # Example
///
/// ```
/// use astro::{angle, coords};
///
/// // Washington, D.C. is at 77°03′56″ WEST, i.e. a positive longitude here.
/// let washington_long = angle::deg_frm_hms(5, 8, 15.7).to_radians();
/// let greenw_sidr = angle::deg_frm_hms(8, 34, 56.853).to_radians();
/// let asc = angle::deg_frm_hms(23, 9, 16.641).to_radians();
///
/// let h = coords::hr_angl_frm_observer_long(greenw_sidr, washington_long, asc);
/// assert!((h.to_degrees() - 64.352133).abs() < 1e-5);
///
/// // Same result with an east-positive longitude.
/// let h_east = coords::hr_angl_frm_observer_east_long(greenw_sidr, -washington_long, asc);
/// assert!((h - h_east).abs() < 1e-12);
/// ```
#[inline]
pub fn hr_angl_frm_observer_long(green_sidreal: f64, observer_long: f64, asc: f64) -> f64 {
    wrap_pi(green_sidreal - observer_long - asc)
}

/// Computes the hour angle from an **east-positive** geographical longitude
/// and the Greenwich sidereal time
///
/// `H = θ₀ + L_east − α`. This is [`hr_angl_frm_observer_long`] for callers
/// whose longitudes are positive east of Greenwich (the IAU and ISO 6709
/// convention); the two functions agree when `observer_east_long` is
/// `−observer_long`.
///
/// # Returns
///
/// * `hour_angle`: Local hour angle, measured westwards from the meridian,
///   in `(−π, π]` *| in radians*
///
/// # Arguments
///
/// * `green_sidreal`: Sidereal time at Greenwich *| in radians*
/// * `observer_east_long`: Observer's geographical longitude, **positive
///   east** *| in radians*
/// * `asc`: Right ascension *| in radians*
#[inline]
pub fn hr_angl_frm_observer_east_long(
    green_sidreal: f64,
    observer_east_long: f64,
    asc: f64,
) -> f64 {
    wrap_pi(green_sidreal + observer_east_long - asc)
}

/// Computes the hour angle from local sidereal time and right
/// ascension
///
/// `H = θ − α` (Meeus, chapter 13).
///
/// # Returns
///
/// * `hour_angle`: Local hour angle, measured westwards from the meridian,
///   in `(−π, π]` *| in radians*
///
/// # Arguments
///
/// * `local_sidreal`: Local sidereal time *| in radians*
/// * `asc`: Right ascension *| in radians*
#[inline]
pub fn hr_angl_frm_loc_sidr(local_sidreal: f64, asc: f64) -> f64 {
    wrap_pi(local_sidreal - asc)
}

/// Computes the ecliptic longitude from equatorial coordinates
///
/// Meeus, eq. 13.1.
///
/// # Returns
///
/// * `ecl_long`: Ecliptic longitude, in `[0, 2π)` *| in radians*
///
/// # Arguments
///
/// * `asc`: Right ascension *| in radians*
/// * `dec`: Declination *| in radians*
/// * `oblq_eclip`: If `asc` and `dec` are corrected for
///   nutation, then *true* obliquity. If not, then
///   *mean* obliquity. *| in radians*
///
/// # Example
///
/// Pollux, J2000.0 (Meeus, example 13.a):
///
/// ```
/// use astro::{angle, coords};
///
/// let asc = angle::deg_frm_hms(7, 45, 18.946).to_radians();
/// let dec = angle::deg_frm_dms(28, 1, 34.26).to_radians();
///
/// let long = coords::ecl_long_frm_eq(asc, dec, coords::OBLQ_J2000);
/// let lat = coords::ecl_lat_frm_eq(asc, dec, coords::OBLQ_J2000);
/// assert!((long.to_degrees() - 113.215630).abs() < 1e-6);
/// assert!((lat.to_degrees() - 6.684170).abs() < 1e-6);
/// ```
pub fn ecl_long_frm_eq(asc: f64, dec: f64, oblq_eclip: f64) -> f64 {
    let (x, y, _) = rotate_x(asc, dec, oblq_eclip);
    wrap_two_pi(y.atan2(x))
}

/// Computes the ecliptic latitude from equatorial coordinates
///
/// Meeus, eq. 13.2.
///
/// # Returns
///
/// * `ecl_lat`: Ecliptic latitude, in `[−π/2, π/2]` *| in radians*
///
/// # Arguments
///
/// * `asc`: Right ascension *| in radians*
/// * `dec`: Declination *| in radians*
/// * `oblq_eclip`: If `asc` and `dec` are corrected for
///   nutation, then *true* obliquity. If not, then
///   *mean* obliquity. *| in radians*
pub fn ecl_lat_frm_eq(asc: f64, dec: f64, oblq_eclip: f64) -> f64 {
    let (x, y, z) = rotate_x(asc, dec, oblq_eclip);
    z.atan2(x.hypot(y))
}

/// Computes ecliptic coordinates from equatorial coordinates
///
/// # Returns
///
/// `(ecl_long, ecl_lat)`
///
/// * `ecl_long`: Ecliptic longitude, in `[0, 2π)` *| in radians*
/// * `ecl_lat`: Ecliptic latitude *| in radians*
///
/// # Arguments
///
/// * `$asc`: Right ascension *| in radians*
/// * `$dec`: Declination *| in radians*
/// * `$oblq_eclip`: If `$asc` and `$dec` are corrected for
///   nutation, then *true* obliquity. If not, then
///   *mean* obliquity. *| in radians*
///
/// See also [`EqPoint::to_ecl`](crate::coords::EqPoint::to_ecl).
#[macro_export]
macro_rules! ecl_frm_eq {
    ($asc: expr_2021, $dec: expr_2021, $oblq_eclip: expr_2021) => {{
        (
            $crate::coords::ecl_long_frm_eq($asc, $dec, $oblq_eclip),
            $crate::coords::ecl_lat_frm_eq($asc, $dec, $oblq_eclip),
        )
    }};
}

/// Computes the right ascension from ecliptic coordinates
///
/// Meeus, eq. 13.3.
///
/// # Returns
///
/// * `asc`: Right ascension, in `[0, 2π)` *| in radians*
///
/// # Arguments
///
/// * `ecl_long`: Ecliptic longitude *| in radians*
/// * `ecl_lat`: Ecliptic latitude *| in radians*
/// * `oblq_eclip`: If `ecl_long` and `ecl_lat` are corrected
///   for nutation, then *true* obliquity. If not, then
///   *mean* obliquity. *| in radians*
pub fn asc_frm_ecl(ecl_long: f64, ecl_lat: f64, oblq_eclip: f64) -> f64 {
    let (x, y, _) = rotate_x(ecl_long, ecl_lat, -oblq_eclip);
    wrap_two_pi(y.atan2(x))
}

/// Computes the declination from ecliptic coordinates
///
/// Meeus, eq. 13.4.
///
/// # Returns
///
/// * `dec`: Declination, in `[−π/2, π/2]` *| in radians*
///
/// # Arguments
///
/// * `ecl_long`: Ecliptic longitude *| in radians*
/// * `ecl_lat`: Ecliptic latitude *| in radians*
/// * `oblq_eclip`: If `ecl_long` and `ecl_lat` are corrected
///   for nutation, then *true* obliquity. If not, then
///   *mean* obliquity. *| in radians*
pub fn dec_frm_ecl(ecl_long: f64, ecl_lat: f64, oblq_eclip: f64) -> f64 {
    let (x, y, z) = rotate_x(ecl_long, ecl_lat, -oblq_eclip);
    z.atan2(x.hypot(y))
}

/// Computes equatorial coordinates from ecliptic coordinates
///
/// # Returns
///
/// `(asc, dec)`
///
/// * `asc`: Right ascension, in `[0, 2π)` *| in radians*
/// * `dec`: Declination *| in radians*
///
/// # Arguments
///
/// * `$ecl_long`: Ecliptic longitude *| in radians*
/// * `$ecl_lat`: Ecliptic latitude *| in radians*
/// * `$oblq_eclip`: If `$ecl_long` and `$ecl_lat` are corrected for
///   nutation, then *true* obliquity. If not, then
///   *mean* obliquity. *| in radians*
///
/// See also [`EclPoint::to_eq`](crate::coords::EclPoint::to_eq).
#[macro_export]
macro_rules! eq_frm_ecl {
    ($ecl_long: expr_2021, $ecl_lat: expr_2021, $oblq_eclip: expr_2021) => {{
        (
            $crate::coords::asc_frm_ecl($ecl_long, $ecl_lat, $oblq_eclip),
            $crate::coords::dec_frm_ecl($ecl_long, $ecl_lat, $oblq_eclip),
        )
    }};
}

/// Cartesian components of a direction in the horizontal frame
/// (x towards the south, y towards the west, z towards the zenith) from its
/// hour angle and declination.
#[inline]
fn hz_vec_frm_eq(hour_angle: f64, dec: f64, observer_lat: f64) -> (f64, f64, f64) {
    let (sin_h, cos_h) = hour_angle.sin_cos();
    let (sin_d, cos_d) = dec.sin_cos();
    let (sin_p, cos_p) = observer_lat.sin_cos();

    (
        cos_h * cos_d * sin_p - sin_d * cos_p,
        sin_h * cos_d,
        sin_p * sin_d + cos_p * cos_d * cos_h,
    )
}

/// Cartesian components of a direction in the hour-angle frame (x towards
/// the meridian on the equator, y towards the west point, z towards the
/// north celestial pole) from its azimuth (from the south) and altitude.
#[inline]
fn eq_vec_frm_hz(az: f64, alt: f64, observer_lat: f64) -> (f64, f64, f64) {
    let (sin_a, cos_a) = az.sin_cos();
    let (sin_h, cos_h) = alt.sin_cos();
    let (sin_p, cos_p) = observer_lat.sin_cos();

    (
        cos_a * cos_h * sin_p + sin_h * cos_p,
        sin_a * cos_h,
        sin_p * sin_h - cos_p * cos_h * cos_a,
    )
}

/// Computes the azimuth from equatorial coordinates
///
/// Meeus, eq. 13.5. The azimuth is measured from the **south**, increasing
/// towards the west; use [`az_north_frm_south`] for a compass bearing.
///
/// # Returns
///
/// * `az`: Azimuth, from the south towards the west, in `(−π, π]`
///   (negative values are east of the meridian) *| in radians*. `0` is
///   returned at the zenith and nadir, where the azimuth is undefined.
///
/// # Arguments
///
/// * `hour_angle`: Local hour angle *| in radians*
/// * `dec`: Declination *| in radians*
/// * `observer_lat`: Observer's geographical latitude *| in radians*
pub fn az_frm_eq(hour_angle: f64, dec: f64, observer_lat: f64) -> f64 {
    let (x, y, _) = hz_vec_frm_eq(hour_angle, dec, observer_lat);
    atan2_pi(y, x)
}

/// Computes the altitude from equatorial coordinates
///
/// Meeus, eq. 13.6. No correction is made for atmospheric refraction
/// (see [`atmos`](crate::atmos)), parallax or the dip of the horizon.
///
/// # Returns
///
/// * `alt`: Altitude, in `[−π/2, π/2]` *| in radians*
///
/// # Arguments
///
/// * `hour_angle`: Local hour angle *| in radians*
/// * `dec`: Declination *| in radians*
/// * `observer_lat`: Observer's geographical latitude *| in radians*
pub fn alt_frm_eq(hour_angle: f64, dec: f64, observer_lat: f64) -> f64 {
    let (x, y, z) = hz_vec_frm_eq(hour_angle, dec, observer_lat);
    z.atan2(x.hypot(y))
}

/// Computes the azimuth and altitude from the hour angle and declination
///
/// Same as [`az_frm_eq`] and [`alt_frm_eq`], computed together (cheaper
/// when both are needed, e.g. in a loop).
///
/// # Returns
///
/// The azimuth (from the south towards the west, in `(−π, π]`) and the
/// altitude *| in radians*
///
/// # Arguments
///
/// * `hour_angle`: Local hour angle *| in radians*
/// * `dec`: Declination *| in radians*
/// * `observer_lat`: Observer's geographical latitude *| in radians*
pub fn hz_frm_hr_angl_and_dec(hour_angle: f64, dec: f64, observer_lat: f64) -> HzPoint {
    let (x, y, z) = hz_vec_frm_eq(hour_angle, dec, observer_lat);
    HzPoint {
        az: atan2_pi(y, x),
        alt: z.atan2(x.hypot(y)),
    }
}

/// Computes local horizontal coordinates from equatorial coordinates
///
/// # Returns
///
/// `(az, alt)`
///
/// * `az`: Azimuth, from the south towards the west *| in radians*
/// * `alt`: Altitude *| in radians*
///
/// # Arguments
///
/// * `$hour_angle`: Hour angle *| in radians*
/// * `$dec`: Declination *| in radians*
/// * `$observer_lat`: Observer's geographical latitude *| in radians*
///
/// See also
/// [`hz_frm_hr_angl_and_dec`](crate::coords::hz_frm_hr_angl_and_dec).
#[macro_export]
macro_rules! loc_hz_frm_eq {
    ($hour_angle: expr_2021, $dec: expr_2021, $observer_lat: expr_2021) => {{
        (
            $crate::coords::az_frm_eq($hour_angle, $dec, $observer_lat),
            $crate::coords::alt_frm_eq($hour_angle, $dec, $observer_lat),
        )
    }};
}

/// Computes the hour angle from local horizontal coordinates
///
/// Meeus, chapter 13 (inverse of eq. 13.5). **The azimuth is measured from
/// the south**, increasing towards the west; convert a compass azimuth
/// (from the north) with [`az_south_frm_north`] first.
///
/// # Returns
///
/// * `hour_angle`: Local hour angle, in `(−π, π]` *| in radians*. `0` is
///   returned for a direction exactly at a celestial pole.
///
/// # Arguments
///
/// * `az`: Azimuth, from the south towards the west *| in radians*
/// * `alt`: Altitude *| in radians*
/// * `observer_lat`: Observer's geographical latitude *| in radians*
pub fn hr_angl_frm_hz(az: f64, alt: f64, observer_lat: f64) -> f64 {
    let (x, y, _) = eq_vec_frm_hz(az, alt, observer_lat);
    atan2_pi(y, x)
}

/// Computes the declination from local horizontal coordinates
///
/// Meeus, chapter 13 (inverse of eq. 13.6):
/// `sin δ = sin φ sin h − cos φ cos h cos A`. **The azimuth is measured from
/// the south**, increasing towards the west.
///
/// Earlier versions of this crate used `cos A · cos A` in place of
/// `cos h · cos A` here and returned wrong declinations (GitHub issue #19).
///
/// # Returns
///
/// * `dec`: Declination, in `[−π/2, π/2]` *| in radians*
///
/// # Arguments
///
/// * `az`: Azimuth, from the south towards the west *| in radians*
/// * `alt`: Altitude *| in radians*
/// * `observer_lat`: Observer's geographical latitude *| in radians*
pub fn dec_frm_hz(az: f64, alt: f64, observer_lat: f64) -> f64 {
    let (x, y, z) = eq_vec_frm_hz(az, alt, observer_lat);
    z.atan2(x.hypot(y))
}

/// Computes the hour angle and declination from local horizontal
/// coordinates
///
/// Same as [`hr_angl_frm_hz`] and [`dec_frm_hz`], computed together and
/// sharing the trigonometric work, for use in tight loops. **The azimuth is
/// measured from the south**, increasing towards the west; for a compass
/// azimuth (from the north, as in ERFA/SOFA `ae2hd`) convert it with
/// [`az_south_frm_north`] first.
///
/// # Returns
///
/// The local hour angle (in `(−π, π]`) and the declination *| in radians*
///
/// # Arguments
///
/// * `az`: Azimuth, from the south towards the west *| in radians*
/// * `alt`: Altitude *| in radians*
/// * `observer_lat`: Observer's geographical latitude *| in radians*
///
/// # Example
///
/// ```
/// use astro::coords;
///
/// let lat = 40_f64.to_radians();
/// // Due east (compass azimuth 90°), on the horizon.
/// let az = coords::az_south_frm_north(90_f64.to_radians());
/// let p = coords::hr_angl_and_dec_frm_hz(az, 0.0, lat);
///
/// // The east point of the horizon is on the celestial equator, six hours
/// // before the meridian.
/// assert!((p.hour_angle.to_degrees() + 90.0).abs() < 1e-12);
/// assert!(p.dec.abs() < 1e-12);
/// ```
pub fn hr_angl_and_dec_frm_hz(az: f64, alt: f64, observer_lat: f64) -> HrAnglDec {
    let (x, y, z) = eq_vec_frm_hz(az, alt, observer_lat);
    HrAnglDec {
        hour_angle: atan2_pi(y, x),
        dec: z.atan2(x.hypot(y)),
    }
}

/// Converts an azimuth measured from the south (Meeus' convention) to one
/// measured from the north (the compass convention)
///
/// Both conventions increase in the same sense (south → west → north for
/// Meeus, north → east → south for the compass, i.e. clockwise seen from
/// above); they differ by π.
///
/// # Returns
///
/// * `az_north`: Azimuth from the north, increasing towards the east, in
///   `[0, 2π)` *| in radians*
///
/// # Arguments
///
/// * `az_south`: Azimuth from the south, increasing towards the west
///   *| in radians*
#[inline]
pub fn az_north_frm_south(az_south: f64) -> f64 {
    wrap_two_pi(az_south + PI)
}

/// Converts an azimuth measured from the north (the compass convention) to
/// one measured from the south (Meeus' convention)
///
/// # Returns
///
/// * `az_south`: Azimuth from the south, increasing towards the west, in
///   `(−π, π]` (the range [`az_frm_eq`] uses) *| in radians*
///
/// # Arguments
///
/// * `az_north`: Azimuth from the north, increasing towards the east
///   *| in radians*
#[inline]
pub fn az_south_frm_north(az_north: f64) -> f64 {
    wrap_pi(az_north + PI)
}

// The galactic system of the IAU (1959), defined in the B1950.0 equatorial
// frame: north galactic pole at α = 192.25°, δ = +27.4°; the galactic
// longitude of the ascending node of the galactic equator on the B1950.0
// equator is 33°, which gives the constants 303° (= 33° + 270°) and
// 123° (= 33° + 90°) below. Meeus, eqs. 13.7 and 13.8.
const GAL_POLE_ASC_B1950: f64 = 192.25 * PI / 180.0;
const GAL_POLE_DEC_B1950: f64 = 27.4 * PI / 180.0;
const GAL_LONG_OFFSET: f64 = 303.0 * PI / 180.0;
const GAL_NCP_LONG: f64 = 123.0 * PI / 180.0;
const GAL_ASC_OFFSET: f64 = 12.25 * PI / 180.0;

/// Components (numerator of tan x, denominator of tan x, sin b) of eqs.
/// 13.7 and 13.8, each multiplied by cos δ so that they stay finite at the
/// celestial poles.
#[inline]
fn gal_vec_frm_eq(asc: f64, dec: f64) -> (f64, f64, f64) {
    let (sin_da, cos_da) = (GAL_POLE_ASC_B1950 - asc).sin_cos();
    let (sin_d, cos_d) = dec.sin_cos();
    let (sin_p, cos_p) = GAL_POLE_DEC_B1950.sin_cos();

    (
        sin_da * cos_d,
        cos_da * cos_d * sin_p - sin_d * cos_p,
        sin_d * sin_p + cos_d * cos_p * cos_da,
    )
}

/// Same as [`gal_vec_frm_eq`] for the inverse transformation.
#[inline]
fn eq_vec_frm_gal(gal_long: f64, gal_lat: f64) -> (f64, f64, f64) {
    let (sin_dl, cos_dl) = (gal_long - GAL_NCP_LONG).sin_cos();
    let (sin_b, cos_b) = gal_lat.sin_cos();
    let (sin_p, cos_p) = GAL_POLE_DEC_B1950.sin_cos();

    (
        sin_dl * cos_b,
        cos_dl * cos_b * sin_p - sin_b * cos_p,
        sin_b * sin_p + cos_b * cos_p * cos_dl,
    )
}

/// Computes the galactic longitude from equatorial coordinates
///
/// Meeus, eq. 13.7, for the IAU 1959 galactic system.
///
/// **The equatorial coordinates must be referred to the standard equinox of
/// B1950.0**, in which the galactic pole was defined (α = 12ʰ49ᵐ,
/// δ = +27.4°). Precess J2000.0 coordinates to B1950.0 first (see
/// [`precess`](crate::precess)).
///
/// # Returns
///
/// * `gal_long`: Galactic longitude, in `[0, 2π)` *| in radians*
///
/// # Arguments
///
/// * `asc`: Right ascension (B1950.0) *| in radians*
/// * `dec`: Declination (B1950.0) *| in radians*
///
/// # Example
///
/// Nova Serpentis 1978 (Meeus, chapter 13):
///
/// ```
/// use astro::{angle, coords};
///
/// let asc = angle::deg_frm_hms(17, 48, 59.74).to_radians();
/// let dec = angle::deg_frm_dms(-14, 43, 8.2).to_radians();
///
/// let (l, b) = astro::gal_frm_eq!(asc, dec);
/// assert!((l.to_degrees() - 12.9593).abs() < 1e-4);
/// assert!((b.to_degrees() - 6.0463).abs() < 1e-4);
/// ```
pub fn gal_long_frm_eq(asc: f64, dec: f64) -> f64 {
    let (x, y, _) = gal_vec_frm_eq(asc, dec);
    wrap_two_pi(GAL_LONG_OFFSET - x.atan2(y))
}

/// Computes the galactic latitude from equatorial coordinates
///
/// Meeus, eq. 13.8, for the IAU 1959 galactic system.
///
/// **The equatorial coordinates must be referred to the standard equinox of
/// B1950.0.**
///
/// # Returns
///
/// * `gal_lat`: Galactic latitude, in `[−π/2, π/2]` *| in radians*
///
/// # Arguments
///
/// * `asc`: Right ascension (B1950.0) *| in radians*
/// * `dec`: Declination (B1950.0) *| in radians*
pub fn gal_lat_frm_eq(asc: f64, dec: f64) -> f64 {
    let (x, y, z) = gal_vec_frm_eq(asc, dec);
    z.atan2(x.hypot(y))
}

/// Computes galactic coordinates from equatorial coordinates
///
/// # Returns
///
/// `(gal_long, gal_lat)`
///
/// * `gal_long`: Galactic longitude, in `[0, 2π)` *| in radians*
/// * `gal_lat`: Galactic latitude *| in radians*
///
/// # Arguments
///
/// * `$asc`: Right ascension *| in radians*
/// * `$dec`: Declination *| in radians*
///
/// The equatorial coordinates passed are assumed to be referred to the
/// standard equinox of B1950.0.
#[macro_export]
macro_rules! gal_frm_eq {
    ($asc: expr_2021, $dec: expr_2021) => {{
        (
            $crate::coords::gal_long_frm_eq($asc, $dec),
            $crate::coords::gal_lat_frm_eq($asc, $dec),
        )
    }};
}

/// Computes the right ascension from galactic coordinates
///
/// Meeus, chapter 13 (inverse of eq. 13.7), for the IAU 1959 galactic
/// system.
///
/// # Returns
///
/// * `asc`: Right ascension, referred to the standard equinox of B1950.0,
///   in `[0, 2π)` *| in radians*
///
/// # Arguments
///
/// * `gal_long`: Galactic longitude *| in radians*
/// * `gal_lat`: Galactic latitude *| in radians*
pub fn asc_frm_gal(gal_long: f64, gal_lat: f64) -> f64 {
    let (x, y, _) = eq_vec_frm_gal(gal_long, gal_lat);
    wrap_two_pi(GAL_ASC_OFFSET + x.atan2(y))
}

/// Computes the declination from galactic coordinates
///
/// Meeus, chapter 13 (inverse of eq. 13.8), for the IAU 1959 galactic
/// system.
///
/// # Returns
///
/// * `dec`: Declination, referred to the standard equinox of B1950.0, in
///   `[−π/2, π/2]` *| in radians*
///
/// # Arguments
///
/// * `gal_long`: Galactic longitude *| in radians*
/// * `gal_lat`: Galactic latitude *| in radians*
pub fn dec_frm_gal(gal_long: f64, gal_lat: f64) -> f64 {
    let (x, y, z) = eq_vec_frm_gal(gal_long, gal_lat);
    z.atan2(x.hypot(y))
}

/// Computes equatorial coordinates from galactic coordinates
///
/// # Returns
///
/// `(asc, dec)`
///
/// * `asc`: Right ascension, in `[0, 2π)` *| in radians*
/// * `dec`: Declination *| in radians*
///
/// The equatorial coordinates returned are referred to the standard
/// equinox of B1950.0.
///
/// # Arguments
///
/// * `$gal_long`: Galactic longitude *| in radians*
/// * `$gal_lat`: Galactic latitude *| in radians*
#[macro_export]
macro_rules! eq_frm_gal {
    ($gal_long: expr_2021, $gal_lat: expr_2021) => {{
        (
            $crate::coords::asc_frm_gal($gal_long, $gal_lat),
            $crate::coords::dec_frm_gal($gal_long, $gal_lat),
        )
    }};
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wrap_two_pi_range() {
        assert_eq!(wrap_two_pi(0.0), 0.0);
        assert_eq!(wrap_two_pi(-1e-300), 0.0);
        assert_eq!(wrap_two_pi(TAU), 0.0);
        assert!((wrap_two_pi(-PI / 2.0) - 1.5 * PI).abs() < 1e-15);
        assert!((wrap_two_pi(7.0 * TAU + 1.0) - 1.0).abs() < 1e-13);
        for i in -1000..1000 {
            let a = wrap_two_pi(i as f64 * 0.37);
            assert!((0.0..TAU).contains(&a));
        }
    }

    #[test]
    fn wrap_pi_range() {
        assert_eq!(wrap_pi(PI), PI);
        assert_eq!(wrap_pi(-PI), PI);
        assert_eq!(wrap_pi(-1e-300), -1e-300);
        assert!((wrap_pi(1.5 * PI) + 0.5 * PI).abs() < 1e-15);
        for i in -1000..1000 {
            let a = wrap_pi(i as f64 * 0.37);
            assert!(a > -PI && a <= PI);
        }
    }

    #[test]
    fn gal_constants_consistent() {
        // 303° − 123° = 180°, and α_pole − 12.25° = 180°.
        assert!((GAL_LONG_OFFSET - GAL_NCP_LONG - PI).abs() < 1e-15);
        assert!((GAL_POLE_ASC_B1950 - GAL_ASC_OFFSET - PI).abs() < 1e-15);
    }
}

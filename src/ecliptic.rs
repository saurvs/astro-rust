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

//! The Earth's ecliptic

use crate::angle;
use crate::time;
use std::f64::consts::PI;

/// Computes the mean obliquity of the ecliptic using
/// J. Laskar's formula
///
/// Meeus, eq. 22.3: a polynomial of degree 10 in U, the time in units of
/// 10 000 Julian years from J2000.0.
///
/// # Returns
///
/// * `mn_oblq`: Mean obliquity of the ecliptic *| in radians*
///
/// The accuracy of `mn_oblq` is estimated to be 0.01 arcseconds
/// for 1000 years before and after 2000 AD, and a few arcseconds
/// for 10000 years before and after 2000 AD.
///
/// The formula is only valid for |U| < 1, that is between the years
/// −8000 and +12000 roughly. Outside that range the polynomial diverges
/// quickly and the result is meaningless (it reaches 90° near
/// U = +2.834), although it is still returned. Use
/// [`mn_oblq_laskar_checked`] to get `None` outside the valid range.
///
/// # Arguments
///
/// * `JD`: Julian (Ephemeris) day
pub fn mn_oblq_laskar(JD: f64) -> f64 {
    let u = time::julian_cent(JD) / 100.0;

    Horner_eval!(
        u,
        angle::deg_frm_dms(23, 26, 21.448),
        -angle::deg_frm_dms(0, 0, 4680.93),
        -angle::deg_frm_dms(0, 0, 1.55),
        angle::deg_frm_dms(0, 0, 1999.25),
        -angle::deg_frm_dms(0, 0, 51.38),
        -angle::deg_frm_dms(0, 0, 249.67),
        -angle::deg_frm_dms(0, 0, 39.05),
        angle::deg_frm_dms(0, 0, 7.12),
        angle::deg_frm_dms(0, 0, 27.87),
        angle::deg_frm_dms(0, 0, 5.79),
        angle::deg_frm_dms(0, 0, 2.45)
    )
    .to_radians()
}

/// Computes the mean obliquity of the ecliptic using J. Laskar's
/// formula, checking the time span
///
/// Same as [`mn_oblq_laskar`], but returns `None` when `JD` is more than
/// 10 000 Julian years from J2000.0 (|U| > 1), where the formula is not
/// valid, or when `JD` is not finite.
///
/// # Arguments
///
/// * `JD`: Julian (Ephemeris) day
pub fn mn_oblq_laskar_checked(JD: f64) -> Option<f64> {
    let u = time::julian_cent(JD) / 100.0;

    if u.abs() <= 1.0 {
        Some(mn_oblq_laskar(JD))
    } else {
        None
    }
}

/// Computes the mean obliquity of the ecliptic using
/// the IAU formula
///
/// Meeus, eq. 22.2 (IAU 1980), a cubic in T, the time in Julian
/// centuries from J2000.0. This is the obliquity that goes with the
/// IAU 1980 nutation of [`nutation::nutation`](crate::nutation::nutation)
/// and the IAU 1976 precession.
///
/// # Returns
///
/// * `mn_oblq`: Mean obliquity of the ecliptic *| in radians*
///
/// The error in `mn_oblq` reaches 1 arcsecond over a period of
/// 2000 years from 2000 AD, and about 10 arcseconds over a period of
/// 4000 years from 2000 AD. For longer spans use [`mn_oblq_laskar`].
///
/// # Arguments
///
/// * `JD`: Julian (Ephemeris) day
///
/// # Example
///
/// Meeus, example 22.a: 1987 April 10, 0h TD.
///
/// ```
/// use astro::ecliptic;
///
/// let eps0 = ecliptic::mn_oblq_IAU(2446895.5).to_degrees();
/// assert!((eps0 - (23.0 + 26.0 / 60.0 + 27.407 / 3600.0)).abs() < 0.0005 / 3600.0);
/// ```
pub fn mn_oblq_IAU(JD: f64) -> f64 {
    // The IAU 1980 polynomial (Meeus, eq. 22.2) is in Julian centuries,
    // not in units of ten thousand years as Laskar's is.
    let t = time::julian_cent(JD);

    Horner_eval!(
        t,
        angle::deg_frm_dms(23, 26, 21.448),
        -angle::deg_frm_dms(0, 0, 46.815),
        -angle::deg_frm_dms(0, 0, 0.00059),
        angle::deg_frm_dms(0, 0, 0.001813)
    )
    .to_radians()
}

/// Computes the true obliquity of the ecliptic
///
/// The true obliquity is the angle between the ecliptic and the true
/// (instantaneous) equator: ε = ε0 + Δε, the mean obliquity of
/// [`mn_oblq_IAU`] plus the nutation in obliquity of
/// [`nutation::nutation`](crate::nutation::nutation) (Meeus, chapter 22).
///
/// # Returns
///
/// * `tru_oblq`: True obliquity of the ecliptic *| in radians*
///
/// # Arguments
///
/// * `JD`: Julian Ephemeris Day
pub fn tru_oblq(JD: f64) -> f64 {
    mn_oblq_IAU(JD) + crate::nutation::nutation(JD).1
}

/// Computes the longitudes of the two ecliptic points on
/// a horizon on Earth
///
/// Meeus, eq. 14.2. See also
/// [`misc::eclip_hz_points`](crate::misc::eclip_hz_points), which returns
/// the same values with named fields.
///
/// # Returns
///
/// `(long_point_1, long_point_2)`, 180° apart, each in `[0, 2π)`
///
/// * `long_point_1`: Longitude of the ecliptic point on the **western**
///   horizon (the setting point, or descendant) *| in radians*
/// * `long_point_2`: Longitude of the ecliptic point on the **eastern**
///   horizon (the rising point, or ascendant) *| in radians*
///
/// # Arguments
///
/// * `oblq_eclip`  : Obliquity of the ecliptic *| in radians*
/// * `observer_lat`: The observer's geographical latitude *| in radians*
/// * `loc_sidreal` : Local sidereal time *| in radians*
pub fn eclip_points_on_hz(oblq_eclip: f64, observer_lat: f64, loc_sidreal: f64) -> (f64, f64) {
    let west = crate::misc::descendant(oblq_eclip, observer_lat, loc_sidreal);

    (west, crate::coords::wrap_two_pi(west + PI))
}

/// Computes the angle between the ecliptic and a horizon
/// on Earth
///
/// Meeus, eq. 14.3: `cos I = cos ε sin φ − sin ε cos φ sin θ`. `I` is the
/// angle between the north pole of the ecliptic and the zenith, so it is
/// usually acute in the northern hemisphere and obtuse in the southern one
/// (where `π − I` is the acute angle between the two circles). Over one
/// sidereal day it swings between `90° − φ − ε` (at θ = 270°) and
/// `90° − φ + ε` (at θ = 90°). It is not the angle the Sun's daily path
/// makes with the horizon; for that see
/// [`planet::earth::angl_betwn_diurnal_path_and_hz`](crate::planet::earth::angl_betwn_diurnal_path_and_hz).
///
/// # Returns
///
/// * `angle`: Angle between the ecliptic and the horizon, in `[0, π]`
///   *| in radians*
///
/// # Arguments
///
/// * `oblq_eclip`  : Obliquity of the ecliptic *| in radians*
/// * `observer_lat`: The observer's geographical latitude *| in radians*
/// * `loc_sidreal` : Local sidereal time *| in radians*
pub fn angl_betwn_eclip_and_hz(oblq_eclip: f64, observer_lat: f64, loc_sidreal: f64) -> f64 {
    (oblq_eclip.cos() * observer_lat.sin()
        - oblq_eclip.sin() * observer_lat.cos() * loc_sidreal.sin())
    .clamp(-1.0, 1.0)
    .acos()
}

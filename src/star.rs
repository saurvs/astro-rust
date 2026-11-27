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

//! Stars

/// Computes the combined magnitude of two stars
///
/// # Arguments
///
/// * `m1`: Magnitude of star 1
/// * `m2`: Magnitude of star 2
#[inline]
pub fn combined_mag(m1: f64, m2: f64) -> f64 {
    m2 - 2.5 * (brightness_ratio(m1, m2) + 1.0)
}

/// Computes the combined magnitude of two or more stars
///
/// # Arguments
///
/// * `m`: Array of magnitudes of stars
pub fn combined_mag_of_many(m: &[f64]) -> f64 {
    let mut sum = 0.0;

    for i in m.iter() {
        sum += 10_f64.powf(-0.4 * i);
    }

    -2.5 * sum.log10()
}

/// Computes the brightness ratio of two stars
///
/// # Arguments
///
/// * `m1`: Magnitude of star 1
/// * `m2`: Magnitude of star 2
#[inline]
pub fn brightness_ratio(m1: f64, m2: f64) -> f64 {
    10.0_f64.powf(0.4 * (m2 - m1))
}

/// Computes the difference in magnitude of two stars
///
/// # Arguments
///
/// * `br`: Brightness ratio of two stars
#[inline]
pub fn mag_diff(br: f64) -> f64 {
    2.5 * br.log10()
}

/// Computes the absolute magnitude of a star from its parallax
///
/// # Arguments
///
/// * `par`: Parallax of the star
/// * `am`: Apparent magnitude of the star
#[inline]
pub fn abs_mag_frm_parallax(par: f64, am: f64) -> f64 {
    am + 5.0 + 5.0 * (par.to_degrees() * 3600.0).log10()
}

/// Computes the absolute magnitude of a star from its distance from earth
///
/// # Arguments
///
/// * `d`: The star's to earth *(parsecs)*
/// * `am`: Apparent magnitude of the star
#[inline]
pub fn abs_mag_frm_dist(d: f64, am: f64) -> f64 {
    am + 5.0 - 5.0 * d.log10()
}

/// Computes the angle, at a star, between the direction of the north
/// celestial pole and the direction of the north pole of the ecliptic
///
/// Meeus, chapter 14 ("Ecliptic and equator"). The angle `q` is the
/// position angle of the north ecliptic pole as seen from the star:
/// measured from the direction of the north celestial pole, positive
/// towards the east. For a star on the ecliptic (latitude 0) it reduces to
/// `tan(q0) = -cos(eclip_long) tan(oblq_eclip)`, the angle between the
/// ecliptic and the east-west direction at that point.
///
/// Up to astro 2.0.0 this function returned `q + π` (wrong quadrant).
///
/// # Returns
///
/// * `angle`: The angle `q` *| in radians*, in `(-π, π]`
///
/// # Arguments
///
/// * `eclip_long`: The star's ecliptical longitude *| in radians*
/// * `eclip_lat`: The star's ecliptical latitude *| in radians*
/// * `oblq_eclip`: Obliquity of the ecliptic *| in radians*
#[inline]
pub fn angl_between_north_celes_and_eclip_pole(
    eclip_long: f64,
    eclip_lat: f64,
    oblq_eclip: f64,
) -> f64 {
    let tan_e = oblq_eclip.tan();

    (-eclip_long.cos() * tan_e).atan2(eclip_lat.cos() - eclip_lat.sin() * eclip_long.sin() * tan_e)
}

/// Computes the equatorial coordinates of a star at a different time,
/// from its motion in space
///
/// Meeus, chapter 21 ("Motion in space"). Instead of assuming a uniform
/// proper motion, this combines the proper motion with the star's
/// distance and radial velocity into a straight-line motion in space
/// relative to the Sun. Over many millennia this is more accurate than
/// a linear proper motion, provided the distance and radial velocity are
/// well known. Precession is not included: the result stays in the
/// reference frame of the initial position.
///
/// Useful conversions (Meeus): a radial velocity in km/s divided by
/// 977 792 gives parsecs per year; a distance in light-years times
/// 0.306 60 gives parsecs; a parallax of `p` arcseconds means a distance
/// of `1/p` parsecs.
///
/// Up to astro 2.0.0 this function multiplied the proper motion in right
/// ascension by 15, as if it had been given in seconds of time.
///
/// # Returns
///
/// `(new_asc, new_dec)`
///
/// * `new_asc`: Right ascension at the different
///   time *| in radians*, in `[0, 2π)`
/// * `new_dec`: Declination at the different
///   time *| in radians*
///
/// # Arguments
///
/// * `asc0`: Right ascension of the star initially *| in radians*
/// * `dec0`: Declination of the star initially *| in radians*
/// * `r`: Distance of the star *| in parsecs*
/// * `delta_r`: Radial velocity of the star *| in parsecs per Julian
///   year*, positive when receding
/// * `proper_motion_asc`: Annual proper motion of the star in right
///   ascension, dα/dt *| in radians per year* (not multiplied by
///   `cos(dec0)`; one second of time per year is 15″ per year)
/// * `proper_motion_dec`: Annual proper motion of the star in declination
///   *| in radians per year*
/// * `t`: Julian years from the initial time; negative in the past
///   and positive in the future
pub fn eq_coords_frm_motion(
    asc0: f64,
    dec0: f64,
    r: f64,
    delta_r: f64,
    proper_motion_asc: f64,
    proper_motion_dec: f64,
    t: f64,
) -> (f64, f64) {
    let (sin_a, cos_a) = asc0.sin_cos();
    let (sin_d, cos_d) = dec0.sin_cos();

    let x = r * cos_d * cos_a;
    let y = r * cos_d * sin_a;
    let z = r * sin_d;

    let delta_asc = proper_motion_asc;
    let delta_dec = proper_motion_dec;

    // (x/r), (y/r), (z/r) written without the division so that r = 0 works
    let delta_x = cos_d * cos_a * delta_r - z * delta_dec * cos_a - y * delta_asc;
    let delta_y = cos_d * sin_a * delta_r - z * delta_dec * sin_a + x * delta_asc;
    let delta_z = sin_d * delta_r + r * delta_dec * cos_d;

    let x1 = x + t * delta_x;
    let y1 = y + t * delta_y;
    let z1 = z + t * delta_z;

    let asc = y1.atan2(x1).rem_euclid(std::f64::consts::TAU);
    let dec = z1.atan2(x1.hypot(y1));

    (asc, dec)
}

/// Computes the proper motion of a star in ecliptical coordinates
///
/// Despite its name, this function converts a proper motion given in
/// equatorial coordinates into ecliptical coordinates. See
/// [`proper_motion_in_ecl_coords`](crate::precess::proper_motion_in_ecl_coords),
/// which does not need the ecliptic latitude as an extra argument.
///
/// # Returns
///
/// `(pmotion_long, pmotion_lat)`: dλ/dt and dβ/dt *| in radians per unit
/// time*
///
/// # Arguments
///
/// * `asc`: Right ascension of the star *| in radians*
/// * `dec`: Declination of the star *| in radians*
/// * `pmotion_asc`: Proper motion in right ascension, dα/dt *| in radians
///   per unit time*
/// * `pmotion_dec`: Proper motion in declination *| in radians per unit
///   time*
/// * `ecl_lat`: Ecliptic latitude of the star *| in radians*
/// * `oblq_eclip`: Obliquity of the ecliptic *| in radians*
#[deprecated(note = "converts to ecliptical, not equatorial, coordinates; use \
            `precess::proper_motion_in_ecl_coords`")]
pub fn proper_motion_in_eq_coords(
    asc: f64,
    dec: f64,
    pmotion_asc: f64,
    pmotion_dec: f64,
    ecl_lat: f64,
    oblq_eclip: f64,
) -> (f64, f64) {
    let ecl_lat_cos = ecl_lat.cos();

    let k = oblq_eclip.cos() * dec.cos() + oblq_eclip.sin() * dec.sin() * asc.sin();

    let pmotion_long = (pmotion_dec * oblq_eclip.sin() * asc.cos() + pmotion_asc * dec.cos() * k)
        / (ecl_lat_cos * ecl_lat_cos);

    let pmotion_lat =
        (pmotion_dec * k - pmotion_asc * oblq_eclip.sin() * asc.cos() * dec.cos()) / ecl_lat_cos;

    (pmotion_long, pmotion_lat)
}

// Apparent place of a star (Meeus, chapter 23), kept in its own file.
mod place;
pub use place::{AberrMethod, ApparentPlace, MeanPlace, apparent_place};

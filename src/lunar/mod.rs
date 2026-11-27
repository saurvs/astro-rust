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

//! The Earth's moon

pub mod apsis;
pub mod declination;
pub mod illumination;
pub mod libration;
pub mod node;
pub mod phase;
pub mod position;
mod sequence;

pub use illumination::Illumination;
pub use libration::{
    Librations, PhysicalEphemeris, SelenographicSun, TopocentLibrations, physical_ephemeris,
};
pub use node::time_of_passage_through_nodes;
pub use phase::{Phase, time_of_phase};
pub use position::{
    ApparentPosition, GeocentricPosition, apparent_position, geocentric_position,
    mean_ascending_node, mean_perigee, true_ascending_node,
};

use crate::angle;
use crate::coords;
use crate::time;

/// Computes the equatorial horizontal parallax of the Moon from its
/// distance, using `sin π = 6378.14 km / Δ` (Meeus, chapter 47)
///
/// # Returns
///
/// * `eq_hz_parllx`: Equatorial horizontal parallax of the
///   Moon *| in radians*. `NaN` if the distance is smaller than the
///   Earth's equatorial radius.
///
/// # Arguments
///
/// * `earth_moon_dist`: Distance between the centres of the Earth and the
///   Moon *| in kilometers*
#[inline]
pub fn eq_hz_parllx(earth_moon_dist: f64) -> f64 {
    (6378.14 / earth_moon_dist).asin()
}

/// Computes the geocentric semidiameter of the Moon
///
/// Rigorous form `sin s = k sin π`, with `sin π = 6378.14 / Δ` and the
/// ratio `k = 0.272481` of the Moon's mean radius to the Earth's equatorial
/// radius (Meeus, *Astronomical Algorithms*, 2nd ed., chapter 55). This is
/// the semidiameter seen from the center of the Earth; see
/// [`topocent_semidiameter`] for an observer on the surface.
///
/// # Returns
///
/// * Geocentric semidiameter *| in radians*
///
/// # Arguments
///
/// * `earth_moon_dist`: Distance between the centers of the Earth and the
///   Moon *| in kilometers*
///
/// # Example
///
/// ```
/// // At the mean distance of 384 400 km the Moon's semidiameter is
/// // about 932.55″ (Meeus' approximation: 358 473 400″ km / Δ).
/// let s = astro::lunar::semidiameter(384_400.0);
/// assert!((s.to_degrees() * 3600.0 - 358_473_400.0 / 384_400.0).abs() < 0.0005);
/// ```
#[inline]
pub fn semidiameter(earth_moon_dist: f64) -> f64 {
    // Bug fix: returned sin s instead of s (about 0.003″ too small) and was
    // documented as "radians per kilometer".
    (0.272481 * eq_hz_parllx(earth_moon_dist).sin()).asin()
}

/// Computes the topocentric semidiameter of the Moon (rigorous)
///
/// An observer on the Earth's surface is nearer to the Moon than the
/// Earth's center (unless the Moon is on the horizon), so the Moon looks
/// slightly larger. With `q = Δ′ / Δ` the ratio of topocentric to
/// geocentric distance, `sin s′ = sin s / q` (Meeus, chapter 55), where `q`
/// is found from the geocentric hour angle and declination as in Meeus'
/// eq. 40.6 and 40.7.
///
/// # Returns
///
/// * Topocentric semidiameter *| in radians*
///
/// # Arguments
///
/// * `geocent_semidia`   : Geocentric semidiameter of the Moon
///   *| in radians* (eg: from [`semidiameter`])
/// * `eq_hz_parllx`      : Equatorial horizontal parallax of the Moon
///   *| in radians* (eg: from [`eq_hz_parllx`])
/// * `geocent_dec`       : Geocentric declination of the Moon
///   *| in radians*
/// * `geocent_hour_angle`: Geocentric hour angle of the Moon, measured
///   westwards from the observer's meridian *| in radians*
/// * `rho_sin_phi`       : `ρ sin φ′` of the observer, in Earth equatorial
///   radii (eg: from `planet::earth::rho_sin_cos_phi`)
/// * `rho_cos_phi`       : `ρ cos φ′` of the observer, in Earth equatorial
///   radii
pub fn topocent_semidiameter(
    geocent_semidia: f64,
    eq_hz_parllx: f64,
    geocent_dec: f64,
    geocent_hour_angle: f64,
    rho_sin_phi: f64,
    rho_cos_phi: f64,
) -> f64 {
    let sin_pi = eq_hz_parllx.sin();
    let a = geocent_dec.cos() * geocent_hour_angle.sin();
    let b = geocent_dec.cos() * geocent_hour_angle.cos() - rho_cos_phi * sin_pi;
    let c = geocent_dec.sin() - rho_sin_phi * sin_pi;
    let q = (a * a + b * b + c * c).sqrt();

    (geocent_semidia.sin() / q).asin()
}

/// Computes the topocentric semidiameter of the Moon from its altitude
/// (approximate)
///
/// `s′ ≈ s (1 + sin h sin π)` (Meeus, chapter 55). The increase is zero
/// with the Moon on the horizon and 14″ to 18″ at the zenith. The error
/// is a small fraction of an arcsecond; use [`topocent_semidiameter`] for
/// the rigorous value.
///
/// # Returns
///
/// * Topocentric semidiameter *| in radians*
///
/// # Arguments
///
/// * `geocent_semidia`: Geocentric semidiameter of the Moon *| in radians*
/// * `eq_hz_parllx`   : Equatorial horizontal parallax of the Moon
///   *| in radians*
/// * `alt`            : Altitude of the Moon above the observer's horizon
///   *| in radians*
#[inline]
pub fn topocent_semidiameter_frm_alt(geocent_semidia: f64, eq_hz_parllx: f64, alt: f64) -> f64 {
    geocent_semidia * (1.0 + alt.sin() * eq_hz_parllx.sin())
}

/// Computes the inclination of the mean lunar equator with the
/// ecliptic
///
/// This is the IAU value `I = 1°32′32.7″` used for the librations (Meeus,
/// chapter 53).
///
/// # Returns
///
/// * `inc`: Inclination of the mean lunar equator with the
///   ecliptic *| in radians*
#[inline]
pub fn inc_of_mn_lunar_eq() -> f64 {
    angle::deg_frm_dms(1, 32, 32.7).to_radians()
}

/// Computes the optical librations of the Moon in longitude and
/// latitude (Meeus, eq. 53.1)
///
/// See [`libration::librations`] for a function that returns the optical,
/// physical and total librations together, and [`libration`] for the sign
/// conventions.
///
/// # Returns
///
/// `(optical_libr_in_long, optical_libr_in_lat)`
///
/// * `optical_libr_in_long`: The optical libration in
///   longitude `l′` *| in radians*, in `(-π, π]`
/// * `optical_libr_in_lat` : The optical libration in
///   latitude `b′` *| in radians*
///
/// # Arguments
///
/// * `JD`                 : Julian Ephemeris Day (TD)
/// * `mn_ecl_long_moon`   : Geocentric ecliptic longitude of the Moon
///   referred to the mean equinox of the date *| in radians*, i.e.
///   *without* the nutation in longitude (`λ − Δψ`), as returned by
///   [`geocent_ecl_pos`]
/// * `apprnt_ecl_lat_moon`: Geocentric ecliptic latitude of the Moon
///   *| in radians*
pub fn optical_libr(JD: f64, mn_ecl_long_moon: f64, apprnt_ecl_lat_moon: f64) -> (f64, f64) {
    let l = libration::librations(JD, mn_ecl_long_moon, apprnt_ecl_lat_moon);

    (l.optical_long, l.optical_lat)
}

/// Computes the physical librations of the Moon in longitude and
/// latitude (Meeus, eq. 53.2)
///
/// # Returns
///
/// `(physical_libr_in_long, physical_libr_in_lat)`
///
/// * `physical_libr_in_long`: The physical libration in
///   longitude `l″` *| in radians*
/// * `physical_libr_in_lat`: The physical libration in
///   latitude `b″` *| in radians*
///
/// Both are small: at most about 0.045° in longitude and 0.065° in
/// latitude.
///
/// # Arguments
///
/// * `JD`                 : Julian Ephemeris Day (TD)
/// * `mn_ecl_long_moon`   : Geocentric ecliptic longitude of the Moon
///   *without* the nutation in longitude (`λ − Δψ`) *| in radians*
/// * `apprnt_ecl_lat_moon`: Geocentric ecliptic latitude of the Moon
///   *| in radians*
/// * `optical_lib_lat`    : The optical libration in latitude `b′`
///   *| in radians*, from [`optical_libr`]
pub fn physical_libr(
    JD: f64,
    mn_ecl_long_moon: f64,
    apprnt_ecl_lat_moon: f64,
    optical_lib_lat: f64,
) -> (f64, f64) {
    let T = time::julian_cent(JD);
    let args = position::fundamental_args_frm_cent(T);
    let node = position::mn_ascend_node_frm_cent(T);
    let rst = libration::rho_sigma_tau(T, &args, node);
    let opt = libration::optical(node, args.arg_of_lat, mn_ecl_long_moon, apprnt_ecl_lat_moon);

    libration::physical(&rst, opt.A, optical_lib_lat)
}

/// Computes the total librations of the Moon in longitude and
/// latitude, the sums of the optical and physical librations (Meeus,
/// chapter 53)
///
/// # Returns
///
/// `(total_libr_in_long, total_libr_in_lat)`
///
/// * `total_libr_in_long`: The total libration in longitude `l`
///   *| in radians*, in `(-π, π]`
/// * `total_libr_in_lat` : The total libration in latitude `b`
///   *| in radians*
///
/// # Arguments
///
/// * `JD`                 : Julian Ephemeris Day (TD)
/// * `mn_ecl_long_moon`   : Geocentric ecliptic longitude of the Moon
///   *without* the nutation in longitude (`λ − Δψ`) *| in radians*
/// * `apprnt_ecl_lat_moon`: Geocentric ecliptic latitude of the Moon
///   *| in radians*
pub fn total_libr(JD: f64, mn_ecl_long_moon: f64, apprnt_ecl_lat_moon: f64) -> (f64, f64) {
    let l = libration::librations(JD, mn_ecl_long_moon, apprnt_ecl_lat_moon);

    (l.total_long, l.total_lat)
}

/// Computes the position angle of the axis of rotation of the Moon,
/// including the effect of the physical libration (Meeus, chapter 53)
///
/// See also [`libration::pos_angl_of_axis`], which needs only the Moon's
/// position.
///
/// # Returns
///
/// * `pos_angl_of_axis_of_rot`: The position angle `P` of the axis
///   of rotation of the Moon *| in radians*, measured eastward from the
///   north point of the disk, in `[0, 2π)` (it is always within about 25°
///   of north)
///
/// # Arguments
///
/// * `JD`                 : Julian Ephemeris Day (TD)
/// * `mn_ascent_node_long`: Longitude of the mean ascending node of the
///   Moon *| in radians* (see [`mean_ascending_node`])
/// * `total_lib_lat`      : Total libration of the Moon in
///   latitude `b` *| in radians*
/// * `nut_in_long`        : Nutation in longitude `Δψ` *| in radians*
/// * `true_oblq`          : True obliquity of the ecliptic `ε`
///   *| in radians*
/// * `apprnt_moon_asc`    : Apparent geocentric right ascension of the
///   Moon *| in radians*
pub fn pos_angl_of_axis_of_rot(
    JD: f64,
    mn_ascent_node_long: f64,
    total_lib_lat: f64,
    nut_in_long: f64,
    true_oblq: f64,
    apprnt_moon_asc: f64,
) -> f64 {
    let T = time::julian_cent(JD);
    let args = position::fundamental_args_frm_cent(T);
    // ρ and σ do not depend on the node; τ, which does, is not used here.
    let rst = libration::rho_sigma_tau(T, &args, mn_ascent_node_long);

    libration::pos_angl(
        mn_ascent_node_long,
        nut_in_long,
        &rst,
        total_lib_lat,
        apprnt_moon_asc,
        true_oblq,
    )
}

/// Computes the corrections that turn the geocentric librations and
/// position angle of the axis into topocentric ones (Meeus, chapter 53,
/// "differential corrections")
///
/// Add the returned corrections to the geocentric values. See also
/// [`libration::topocent_librations`] and
/// [`libration::topocent_librations_by_diff_corrections`], which work from
/// a Julian Ephemeris Day.
///
/// # Returns
///
/// `(topocent_libr_in_long, topocent_libr_in_lat, topocent_libr_in_P)`
///
/// * `topocent_libr_in_long`: Correction `Δl` to the libration in
///   longitude *| in radians*
/// * `topocent_libr_in_lat` : Correction `Δb` to the libration in
///   latitude *| in radians*
/// * `topocent_libr_in_P`   : Correction `ΔP` to the position angle of the
///   axis of rotation *| in radians*
///
/// # Arguments
///
/// * `observer_lat`          : Geographic latitude of the observer
///   *| in radians*
/// * `geocent_dec_moon`      : Geocentric declination of the Moon
///   *| in radians*
/// * `local_hour_angl`       : Local hour angle of the Moon, from the local
///   sidereal time and the *geocentric* right ascension *| in radians*
/// * `geocent_hz_parllx_moon`: Geocentric equatorial horizontal parallax of
///   the Moon *| in radians*
/// * `pos_angl_axis_of_rot`  : Geocentric position angle of the axis of
///   rotation of the Moon *| in radians*
/// * `total_lib_lat`         : Geocentric total libration of the Moon in
///   latitude *| in radians*
pub fn topocent_libr_by_diff_corrections(
    observer_lat: f64,
    geocent_dec_moon: f64,
    local_hour_angl: f64,
    geocent_hz_parllx_moon: f64,
    pos_angl_axis_of_rot: f64,
    total_lib_lat: f64,
) -> (f64, f64, f64) {
    libration::topocent_corrections(
        observer_lat,
        geocent_dec_moon,
        local_hour_angl,
        geocent_hz_parllx_moon,
        pos_angl_axis_of_rot,
        total_lib_lat,
    )
}

/// Computes the geocentric ecliptic position of the Moon,
/// referred to the mean equinox of the date (Meeus, chapter 47)
///
/// See [`position::geocentric_position`], which also returns the parallax,
/// and [`position::apparent_position`] for the apparent position.
///
/// # Returns
///
/// `(moon_ecl_point, rad_vec)`
///
/// * `moon_ecl_point`: Ecliptic point of the Moon *| in radians*: the
///   longitude, in `[0, 2π)`, does *not* include nutation; the latitude
///   is in `[-π/2, π/2]`
/// * `rad_vec`: Distance between the centres of the Earth and the Moon
///   *| in kilometers*
///
/// This function uses the principal terms of the Chapront ELP-2000/82
/// lunar theory, so the accuracy of:
///
/// * `moon_ecl_point.long`: is about `10` arcseconds
/// * `moon_ecl_point.lat`: is about `4` arcseconds
///
/// # Arguments
///
/// * `JD`: Julian Ephemeris Day (TD)
///
/// # Example
///
/// ```
/// // Meeus, example 47.a: 1992 April 12, 0h TD
/// let (moon, dist) = astro::lunar::geocent_ecl_pos(2448724.5);
/// assert!((moon.long.to_degrees() - 133.162655).abs() < 1e-6);
/// assert!((moon.lat.to_degrees() - -3.229126).abs() < 1e-6);
/// assert!((dist - 368409.7).abs() < 0.1);
/// ```
pub fn geocent_ecl_pos(JD: f64) -> (coords::EclPoint, f64) {
    let p = position::geocentric_position(JD);

    (
        coords::EclPoint {
            long: p.long,
            lat: p.lat,
        },
        p.dist,
    )
}

/// Computes the longitude of the mean ascending node of the Moon
/// (Meeus, eq. 47.7)
///
/// Note that this function takes Julian *centuries*; see
/// [`mean_ascending_node`] for the same quantity from a Julian Ephemeris
/// Day.
///
/// # Returns
///
/// * `long_mn_ascend_node`: Longitude of the mean ascending
///   node *| in radians*, in `[0, 2π)`, measured from the mean equinox of
///   the date
///
/// # Arguments
///
/// * `JC`: Julian centuries (TD) since J2000.0, as given by
///   [`time::julian_cent`]
pub fn mn_ascend_node(JC: f64) -> f64 {
    position::mn_ascend_node_frm_cent(JC)
}

/// Computes the longitude of the true ascending node of the Moon
/// (Meeus, chapter 47)
///
/// Note that this function takes Julian *centuries*; see
/// [`true_ascending_node`] for the same quantity from a Julian Ephemeris
/// Day.
///
/// # Returns
///
/// * `long_true_ascend_node`: Longitude of the true ascending
///   node *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `JC`: Julian centuries (TD) since J2000.0, as given by
///   [`time::julian_cent`]
pub fn true_ascend_node(JC: f64) -> f64 {
    position::true_ascend_node_frm_cent(JC)
}

/// Computes the longitude of the mean perigee of the Moon
/// (Meeus, chapter 47)
///
/// Note that this function takes Julian *centuries*; see
/// [`mean_perigee`] for the same quantity from a Julian Ephemeris Day.
///
/// # Returns
///
/// * `long_mn_perigee`: Longitude of mean perigee
///   *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `JC`: Julian centuries (TD) since J2000.0, as given by
///   [`time::julian_cent`]
pub fn mn_perigee(JC: f64) -> f64 {
    position::mn_perigee_frm_cent(JC)
}

/// Computes the position angle of the Moon's bright limb
/// (Meeus, eq. 48.5)
///
/// See [`illumination::bright_limb_pos_angl`].
///
/// # Returns
///
/// * `pos_angl_of_bright_limb` The position angle of the midpoint
///   of the illuminated limb of the Moon, measured eastward from the north
///   point of the disk *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `sun_eq_point`: Apparent equatorial coordinates of the Sun
///   *| in radians*
/// * `moon_eq_point`: Apparent equatorial coordinates of the Moon
///   *| in radians*
pub fn bright_limb(sun_eq_point: coords::EqPoint, moon_eq_point: coords::EqPoint) -> f64 {
    illumination::bright_limb_pos_angl(
        sun_eq_point.asc,
        sun_eq_point.dec,
        moon_eq_point.asc,
        moon_eq_point.dec,
    )
}

/// Computes the illuminated fraction of the lunar disk, using equatorial
/// coordinates (Meeus, eqs. 48.1 to 48.3)
///
/// # Returns
///
/// The illuminated fraction `k`, in `[0, 1]`.
///
/// # Arguments
///
/// * `sun_eq_point`   : Geocentric equatorial coordinates of the Sun
///   *| in radians*
/// * `moon_eq_point`  : Geocentric equatorial coordinates of the Moon
///   *| in radians*
/// * `earth_moon_dist`: Earth-Moon distance (in any unit, but same as
///   that of `earth_sun_dist`)
/// * `earth_sun_dist` : Earth-Sun distance (in any unit, but same as
///   that of `earth_moon_dist`)
pub fn illum_frac_frm_eq_coords(
    sun_eq_point: &coords::EqPoint,
    moon_eq_point: &coords::EqPoint,
    earth_moon_dist: f64,
    earth_sun_dist: f64,
) -> f64 {
    illuminated_frac(
        illumination::geocent_elong_frm_eq(
            sun_eq_point.asc,
            sun_eq_point.dec,
            moon_eq_point.asc,
            moon_eq_point.dec,
        ),
        earth_moon_dist,
        earth_sun_dist,
    )
}

/// Computes the illuminated fraction of the lunar disk, using ecliptic
/// coordinates (Meeus, eqs. 48.1 to 48.3)
///
/// The Sun's ecliptic latitude is neglected.
///
/// # Returns
///
/// The illuminated fraction `k`, in `[0, 1]`.
///
/// # Arguments
///
/// * `moon_long`      : Geocentric ecliptic longitude of the Moon
///   *| in radians*
/// * `moon_lat`       : Geocentric ecliptic latitude of the Moon
///   *| in radians*
/// * `sun_long`       : Geocentric ecliptic longitude of the Sun
///   *| in radians*
/// * `earth_moon_dist`: Earth-Moon distance (in any unit, but same as
///   that of `earth_sun_dist`)
/// * `earth_sun_dist` : Earth-Sun distance (in any unit, but same as
///   that of `earth_moon_dist`)
pub fn illum_frac_frm_ecl_coords(
    moon_long: f64,
    moon_lat: f64,
    sun_long: f64,
    earth_moon_dist: f64,
    earth_sun_dist: f64,
) -> f64 {
    illuminated_frac(
        illumination::geocent_elong_frm_ecl(moon_long, moon_lat, sun_long),
        earth_moon_dist,
        earth_sun_dist,
    )
}

#[inline]
fn illuminated_frac(moon_geocent_elong: f64, earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    let i = (earth_sun_dist * moon_geocent_elong.sin())
        .atan2(earth_moon_dist - earth_sun_dist * moon_geocent_elong.cos());

    (1.0 + i.cos()) / 2.0
}

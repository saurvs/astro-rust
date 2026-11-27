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
pub mod position;

pub use illumination::Illumination;
pub use libration::{
    Librations, PhysicalEphemeris, SelenographicSun, TopocentLibrations, physical_ephemeris,
};
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

/// Computes the equatorial semidiameter of the Moon
///
/// # Returns
///
/// * `eq_semidiameter`: Geocentric equatorial semidiameter
///   *| in radians per kilometers*
///
/// # Arguments
///
/// * `earth_moon_dist`: Earth-Moon distance *| in kilometers*
#[inline]
pub fn semidiameter(earth_moon_dist: f64) -> f64 {
    0.272481 * eq_hz_parllx(earth_moon_dist).sin()
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

#[inline]
fn E(JC: f64) -> f64 {
    1.0 - JC * (0.002_516 + JC * 0.000_007_4)
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
/// Both are always smaller than about 0.04°.
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

/// Computes the times of passage of the Moon through the ascending and
/// descending nodes, close to a given date
///
/// # Returns
///
/// `(time_of_ascend_node, time_of_descend_node)`
///
/// * `time_of_ascend_node` : Time of passage through the ascending node
/// * `time_of_descend_node`: Time of passage through the descending node
///
/// # Arguments
///
/// `date`: The Date
#[inline]
pub fn time_of_passage_through_nodes(date: &time::Date) -> (f64, f64) {
    let k = 13.4223 * (time::decimal_year(date) - 2000.05);
    let T = k / 1342.23;
    let k1 = (k as i32) as f64;
    let k2 = k1 + 0.5;

    (
        time_of_passage_through_node(k1, T),
        time_of_passage_through_node(k2, T),
    )
}

fn time_of_passage_through_node(k: f64, T: f64) -> f64 {
    let D = (183.638 + 331.73735682 * k + T * T * (0.0014852 + T * (0.00000209 - T * 0.00000001)))
        .to_radians();
    let D_times_2 = 2.0 * D;
    let M = (17.4006 + 26.8203725 * k + T * T * (0.0001186 + T * 0.00000006)).to_radians();
    let M1 =
        (38.3776 + 355.52747313 * k + T * T * (0.0123499 + T * (0.000014627 - T * 0.000000069)))
            .to_radians();
    let sigma = (123.9767 - 1.44098956 * k
        + T * T * (0.0020608 + T * (0.00000214 - T * 0.000000016)))
        .to_radians();
    let P = sigma + (272.75 - T * 2.3).to_radians();
    let V = (299.75 + T * (132.85 - T * 0.009173)).to_radians();

    2451565.1619 + 27.212220817 * k + T * T * (0.0002762 + T * (0.000000021 - T * 0.000000000088))
        - 0.4721 * M1.sin()
        - 0.1649 * D_times_2.sin()
        - 0.0868 * (D_times_2 - M1).sin()
        + 0.0084 * (D_times_2 + M1).sin()
        - 0.0083 * (D_times_2 - M).sin()
        - 0.0039 * (D_times_2 - M - M1).sin()
        + 0.0034 * (2.0 * M1).sin()
        - 0.0031 * (2.0 * (D - M1)).sin()
        + 0.0030 * (D_times_2 + M).sin()
        + 0.0028 * (M - M1).sin()
        + 0.0026 * M.sin()
        + 0.0025 * (2.0 * D_times_2).sin()
        + 0.0024 * D.sin()
        + 0.0022 * (M + M1).sin()
        + 0.0017 * sigma.sin()
        + 0.0014 * (2.0 * D_times_2 - M1).sin()
        + 0.0005 * (D_times_2 + M - M1).sin()
        + 0.0004 * (D_times_2 - M + M1).sin()
        + 0.0003 * ((2.0 * (M - D)).sin() + (2.0 * D_times_2 - M).sin() + V.sin() + P.sin())
}

/// Represents a phase of the Moon
pub enum Phase {
    /// New Moon
    New,
    /// First Quarter
    First,
    /// Full Moon
    Full,
    /// Last Quarter
    Last,
}

/// Computes the Julian day corresponding to one of the four phases
/// of the Moon
///
/// # Returns
///
/// * `JD`: Julian day corresponding to the exact time of the phase
///   that is closest to `date`
///
/// *Meeus* says the mean error in `JD` for phases between 1980 AD
/// and mid-2020 AD is 3.8 seconds.
///
/// # Arguments
///
/// * `date` : Date of interest, close to the phase
/// * `phase`: The [Phase](./enum.Phase.html)
pub fn time_of_phase(date: &time::Date, phase: &Phase) -> f64 {
    let mut K = 12.3685 * (time::decimal_year(date) - 2000.0);
    K = (K as i64) as f64;

    let k = match *phase {
        Phase::New => K,
        Phase::First => K + 0.25,
        Phase::Full => K + 0.5,
        Phase::Last => K + 0.75,
    };

    let T = k / 1236.85;

    let mut JD = 2451550.09766
        + k * 29.530588861
        + T * Horner_eval!(T, 0.0, 0.00015437, -0.00000015, 0.00000000073);

    let E = E(T);

    let M =
        (2.5534 + k * 29.1053567 + T * Horner_eval!(T, 0.0, -0.0000014, -0.00000011)).to_radians();

    let M1 = (201.5643
        + k * 385.81693528
        + T * Horner_eval!(T, 0.0, 0.0107582, 0.00001238, -0.000000058))
    .to_radians();

    let F = (160.7108
        + k * 390.67050284
        + T * Horner_eval!(T, 0.0, -0.0016118, -0.00000227, 0.000000011))
    .to_radians();

    let omega =
        (124.7746 - k * 1.56375588 + T * Horner_eval!(T, 0.0, 0.0020672, 0.00000215)).to_radians();

    let A1 = (299.77 + 0.107408 * k - 0.009173 * T * T).to_radians();
    let A2 = (251.88 + 0.016321 * k).to_radians();
    let A3 = (251.83 + 26.651886 * k).to_radians();
    let A4 = (349.42 + 36.412478 * k).to_radians();
    let A5 = (84.66 + 18.206239 * k).to_radians();
    let A6 = (141.74 + 53.303771 * k).to_radians();
    let A7 = (207.14 + 2.453732 * k).to_radians();
    let A8 = (154.84 + 7.306860 * k).to_radians();
    let A9 = (34.52 + 27.261239 * k).to_radians();
    let A10 = (207.19 + 0.121824 * k).to_radians();
    let A11 = (291.34 + 1.844379 * k).to_radians();
    let A12 = (161.72 + 24.198154 * k).to_radians();
    let A13 = (239.56 + 25.513099 * k).to_radians();
    let A14 = (331.55 + 3.592518 * k).to_radians();

    let is_quarter = match phase {
        &Phase::Last | &Phase::First => true,
        &_ => false,
    };

    if is_quarter {
        let W = 0.00306 - 0.00038 * E * M.cos() + 0.00026 * M1.cos()
            - 0.00002 * ((M1 - M).cos() - (M1 + M).cos() - (2.0 * F).cos());

        JD += match *phase {
            Phase::Last => -W,
            Phase::First => W,
            _ => 0.0,
        };

        let corrections = [
            [-0.62801, M1],
            [0.17172 * E, M],
            [-0.01183 * E, M1 + M],
            [0.00862, 2.0 * M1],
            [0.00804, 2.0 * F],
            [0.00454 * E, M1 - M],
            [0.00204 * E * E, 2.0 * M],
            [-0.0018, M1 - 2.0 * F],
            [-0.0007, M1 + 2.0 * F],
            [-0.0004, 3.0 * M1],
            [-0.00034, 2.0 * M1 - M],
            [0.00032 * E, M + 2.0 * F],
            [0.00032 * E, M - 2.0 * F],
            [-0.00028 * E * E, M1 + 2.0 * M],
            [0.00027 * E, 2.0 * M1 + M],
            [-0.00017, omega],
            [-0.00005, M1 - M - 2.0 * F],
            [0.00004, 2.0 * (M1 + F)],
            [-0.00004, M1 + M + 2.0 * F],
            [0.00004, M1 - 2.0 * M],
            [0.00003, M1 + M - 2.0 * F],
            [0.00003, 3.0 * M],
            [0.00002, 2.0 * (M1 - F)],
            [0.00002, M1 - M + 2.0 * F],
            [-0.00002, 3.0 * M1 + M],
        ];

        for &x in corrections.iter() {
            JD += x[0] * x[1].sin();
        }
    } else {
        let is_new = matches!(*phase, Phase::New);

        let sine_arguments = [
            M1,
            M,
            2.0 * M1,
            2.0 * F,
            M1 - M,
            M1 + M,
            2.0 * M,
            M1 - 2.0 * F,
            M1 + 2.0 * F,
            2.0 * M1 + M,
            3.0 * M1,
            M + 2.0 * F,
            M - 2.0 * F,
            2.0 * M1 - M,
            omega,
            M1 + 2.0 * M,
            2.0 * (M1 - F),
            3.0 * M,
            M1 + M - 2.0 * F,
            2.0 * (M1 + F),
            M1 + M + 2.0 * F,
            M1 - M + 2.0 * F,
            M1 - M - 2.0 * F,
            3.0 * M1 + M,
            4.0 * M1,
        ];

        let multipliers = if is_new {
            [
                -0.4072,
                0.17241 * E,
                0.01608,
                0.01039,
                0.00739 * E,
                -0.00514 * E,
                0.00208 * E * E,
                -0.00111,
                -0.00057,
                0.00056 * E,
                -0.00042,
                0.00042 * E,
                0.00038 * E,
                -0.00024 * E,
                -0.00017,
                -0.00007,
                0.00004,
                0.00004,
                0.00003,
                0.00003,
                -0.00003,
                0.00003,
                -0.00002,
                -0.00002,
                0.00002,
            ]
        } else {
            [
                -0.40614,
                0.17302 * E,
                0.01614,
                0.01043,
                0.00734 * E,
                -0.00515 * E,
                0.00209 * E * E,
                -0.00111,
                -0.00057,
                0.00056 * E,
                -0.00042,
                0.00042 * E,
                0.00038 * E,
                -0.00024 * E,
                -0.00017,
                -0.00007,
                0.00004,
                0.00004,
                0.00003,
                0.00003,
                -0.00003,
                0.00003,
                -0.00002,
                -0.00002,
                0.00002,
            ]
        };

        for mx in multipliers.iter().zip(sine_arguments.iter()) {
            JD += mx.0 * (mx.1).sin();
        }
    }

    let additional_corrections = [
        [0.000_325, A1],
        [0.000_165, A2],
        [0.000_164, A3],
        [0.000_126, A4],
        [0.000_110, A5],
        [0.000_062, A6],
        [0.000_06, A7],
        [0.000_056, A8],
        [0.000_047, A9],
        [0.000_042, A10],
        [0.000_04, A11],
        [0.000_037, A12],
        [0.000_035, A13],
        [0.000_023, A14],
    ];
    for &x in additional_corrections.iter() {
        JD += x[0] * x[1].sin();
    }

    JD
}

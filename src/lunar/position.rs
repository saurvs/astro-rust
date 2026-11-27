//! Geocentric position of the Moon (Meeus, *Astronomical Algorithms*,
//! 2nd ed., chapter 47).
//!
//! The position is computed from the principal periodic terms of the
//! Chapront ELP-2000/82 lunar theory, with the improved polynomial mean
//! arguments of Chapront, Chapront-Touzé and Francou (1998). Meeus quotes an
//! accuracy of about 10″ in longitude and 4″ in latitude.
//!
//! The entry points are [`geocentric_position`] (geometric, mean equinox of
//! date) and [`apparent_position`] (with nutation, plus right ascension and
//! declination). The longitudes of the lunar node and perigee are given by
//! [`mean_ascending_node`], [`true_ascending_node`] and [`mean_perigee`].
//!
//! # Example
//!
//! ```
//! use astro::lunar::position;
//!
//! // Meeus, example 47.a: 1992 April 12, 0h TD
//! let moon = position::geocentric_position(2448724.5);
//! assert!((moon.long.to_degrees() - 133.162655).abs() < 1e-6);
//! assert!((moon.lat.to_degrees() - -3.229126).abs() < 1e-6);
//! assert!((moon.dist - 368409.7).abs() < 0.1);
//!
//! let app = position::apparent_position(2448724.5);
//! assert!((app.asc.to_degrees() - 134.688470).abs() < 1e-5);
//! assert!((app.dec.to_degrees() - 13.768368).abs() < 1e-5);
//! ```

use std::f64::consts::{PI, TAU};

use crate::coords;
use crate::ecliptic;
use crate::nutation;
use crate::time;

/// Equatorial radius of the Earth used by Meeus for the lunar parallax
/// (eq. 47 text), in kilometres.
const EARTH_EQ_RAD_KM: f64 = 6378.14;

/// Constant part of the Earth–Moon distance in the ELP-2000/82 truncation,
/// in kilometres.
const MEAN_DIST_KM: f64 = 385000.56;

/// Fundamental arguments of the lunar theory at an instant.
///
/// All angles are in radians, normalised to `[0, 2π)`. They are referred to
/// the mean equinox of the date.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FundamentalArgs {
    /// The Moon's mean longitude `L′` (eq. 47.1). It already contains the
    /// constant part (−0.70″) of the effect of light-time.
    pub mn_long: f64,
    /// Mean elongation of the Moon from the Sun, `D` (eq. 47.2)
    pub mn_elong: f64,
    /// Mean anomaly of the Sun (of the Earth), `M` (eq. 47.3)
    pub sun_mn_anom: f64,
    /// Mean anomaly of the Moon, `M′` (eq. 47.4)
    pub mn_anom: f64,
    /// The Moon's argument of latitude `F`, its mean distance from the
    /// ascending node (eq. 47.5)
    pub arg_of_lat: f64,
    /// Argument `A1` of the additive terms due to Venus
    pub a1: f64,
    /// Argument `A2` of the additive term due to Jupiter
    pub a2: f64,
    /// Argument `A3` of an additive term in latitude
    pub a3: f64,
    /// Factor `E` (eq. 47.6) for the decreasing eccentricity of the Earth's
    /// orbit. Dimensionless, close to 1.
    pub ecc_fac: f64,
}

/// Sums of the periodic terms of the lunar theory at an instant, including
/// the additive terms due to Venus, Jupiter and the flattening of the Earth.
///
/// Meeus tabulates these sums (`Σl`, `Σb`, `Σr`) in units of 10⁻⁶ degree and
/// 10⁻³ km; here they are already converted to radians and kilometres.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PeriodicSums {
    /// Periodic part of the ecliptic longitude, `Σl` *| in radians*
    pub long: f64,
    /// Ecliptic latitude, `Σb` *| in radians*
    pub lat: f64,
    /// Periodic part of the Earth–Moon distance, `Σr` *| in kilometres*
    pub dist: f64,
}

/// Geometric geocentric position of the Moon, referred to the mean ecliptic
/// and equinox of the date.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeocentricPosition {
    /// Ecliptic longitude `λ` *| in radians*, in `[0, 2π)`, without
    /// nutation
    pub long: f64,
    /// Ecliptic latitude `β` *| in radians*, in `[-π/2, π/2]`
    pub lat: f64,
    /// Distance between the centres of the Earth and the Moon, `Δ`
    /// *| in kilometres*
    pub dist: f64,
    /// Equatorial horizontal parallax `π` *| in radians*
    pub eq_hz_parllx: f64,
}

/// Apparent geocentric position of the Moon, referred to the true equator
/// and equinox of the date.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ApparentPosition {
    /// Apparent ecliptic longitude `λ + Δψ` *| in radians*, in `[0, 2π)`
    pub long: f64,
    /// Ecliptic latitude `β` *| in radians*, in `[-π/2, π/2]`
    pub lat: f64,
    /// Distance between the centres of the Earth and the Moon
    /// *| in kilometres*
    pub dist: f64,
    /// Equatorial horizontal parallax *| in radians*
    pub eq_hz_parllx: f64,
    /// Apparent right ascension *| in radians*, in `[0, 2π)`
    pub asc: f64,
    /// Apparent declination *| in radians*, in `[-π/2, π/2]`
    pub dec: f64,
    /// Nutation in longitude `Δψ` that was applied *| in radians*
    pub nut_long: f64,
    /// True obliquity of the ecliptic `ε = ε0 + Δε` used for the
    /// conversion to equatorial coordinates *| in radians*
    pub true_oblq: f64,
}

/// Periodic terms for the longitude and distance (Meeus, table 47.A).
///
/// Each row is: multiples of `D`, `M`, `M′`, `F`; coefficient of the sine
/// for `Σl` (10⁻⁶ degree); coefficient of the cosine for `Σr` (10⁻³ km).
const TERMS_LR: [(i8, i8, i8, i8, i32, i32); 60] = [
    (0, 0, 1, 0, 6288774, -20905355),
    (2, 0, -1, 0, 1274027, -3699111),
    (2, 0, 0, 0, 658314, -2955968),
    (0, 0, 2, 0, 213618, -569925),
    (0, 1, 0, 0, -185116, 48888),
    (0, 0, 0, 2, -114332, -3149),
    (2, 0, -2, 0, 58793, 246158),
    (2, -1, -1, 0, 57066, -152138),
    (2, 0, 1, 0, 53322, -170733),
    (2, -1, 0, 0, 45758, -204586),
    (0, 1, -1, 0, -40923, -129620),
    (1, 0, 0, 0, -34720, 108743),
    (0, 1, 1, 0, -30383, 104755),
    (2, 0, 0, -2, 15327, 10321),
    (0, 0, 1, 2, -12528, 0),
    (0, 0, 1, -2, 10980, 79661),
    (4, 0, -1, 0, 10675, -34782),
    (0, 0, 3, 0, 10034, -23210),
    (4, 0, -2, 0, 8548, -21636),
    (2, 1, -1, 0, -7888, 24208),
    (2, 1, 0, 0, -6766, 30824),
    (1, 0, -1, 0, -5163, -8379),
    (1, 1, 0, 0, 4987, -16675),
    (2, -1, 1, 0, 4036, -12831),
    (2, 0, 2, 0, 3994, -10445),
    (4, 0, 0, 0, 3861, -11650),
    (2, 0, -3, 0, 3665, 14403),
    (0, 1, -2, 0, -2689, -7003),
    (2, 0, -1, 2, -2602, 0),
    (2, -1, -2, 0, 2390, 10056),
    (1, 0, 1, 0, -2348, 6322),
    (2, -2, 0, 0, 2236, -9884),
    (0, 1, 2, 0, -2120, 5751),
    (0, 2, 0, 0, -2069, 0),
    (2, -2, -1, 0, 2048, -4950),
    (2, 0, 1, -2, -1773, 4130),
    (2, 0, 0, 2, -1595, 0),
    (4, -1, -1, 0, 1215, -3958),
    (0, 0, 2, 2, -1110, 0),
    (3, 0, -1, 0, -892, 3258),
    (2, 1, 1, 0, -810, 2616),
    (4, -1, -2, 0, 759, -1897),
    (0, 2, -1, 0, -713, -2117),
    (2, 2, -1, 0, -700, 2354),
    (2, 1, -2, 0, 691, 0),
    (2, -1, 0, -2, 596, 0),
    (4, 0, 1, 0, 549, -1423),
    (0, 0, 4, 0, 537, -1117),
    (4, -1, 0, 0, 520, -1571),
    (1, 0, -2, 0, -487, -1739),
    (2, 1, 0, -2, -399, 0),
    (0, 0, 2, -2, -381, -4421),
    (1, 1, 1, 0, 351, 0),
    (3, 0, -2, 0, -340, 0),
    (4, 0, -3, 0, 330, 0),
    (2, -1, 2, 0, 327, 0),
    (0, 2, 1, 0, -323, 1165),
    (1, 1, -1, 0, 299, 0),
    (2, 0, 3, 0, 294, 0),
    (2, 0, -1, -2, 0, 8752),
];

/// Periodic terms for the latitude (Meeus, table 47.B).
///
/// Each row is: multiples of `D`, `M`, `M′`, `F`; coefficient of the sine
/// for `Σb` (10⁻⁶ degree).
const TERMS_B: [(i8, i8, i8, i8, i32); 60] = [
    (0, 0, 0, 1, 5128122),
    (0, 0, 1, 1, 280602),
    (0, 0, 1, -1, 277693),
    (2, 0, 0, -1, 173237),
    (2, 0, -1, 1, 55413),
    (2, 0, -1, -1, 46271),
    (2, 0, 0, 1, 32573),
    (0, 0, 2, 1, 17198),
    (2, 0, 1, -1, 9266),
    (0, 0, 2, -1, 8822),
    (2, -1, 0, -1, 8216),
    (2, 0, -2, -1, 4324),
    (2, 0, 1, 1, 4200),
    (2, 1, 0, -1, -3359),
    (2, -1, -1, 1, 2463),
    (2, -1, 0, 1, 2211),
    (2, -1, -1, -1, 2065),
    (0, 1, -1, -1, -1870),
    (4, 0, -1, -1, 1828),
    (0, 1, 0, 1, -1794),
    (0, 0, 0, 3, -1749),
    (0, 1, -1, 1, -1565),
    (1, 0, 0, 1, -1491),
    (0, 1, 1, 1, -1475),
    (0, 1, 1, -1, -1410),
    (0, 1, 0, -1, -1344),
    (1, 0, 0, -1, -1335),
    (0, 0, 3, 1, 1107),
    (4, 0, 0, -1, 1021),
    (4, 0, -1, 1, 833),
    (0, 0, 1, -3, 777),
    (4, 0, -2, 1, 671),
    (2, 0, 0, -3, 607),
    (2, 0, 2, -1, 596),
    (2, -1, 1, -1, 491),
    (2, 0, -2, 1, -451),
    (0, 0, 3, -1, 439),
    (2, 0, 2, 1, 422),
    (2, 0, -3, -1, 421),
    (2, 1, -1, 1, -366),
    (2, 1, 0, 1, -351),
    (4, 0, 0, 1, 331),
    (2, -1, 1, 1, 315),
    (2, -2, 0, -1, 302),
    (0, 0, 1, 3, -283),
    (2, 1, 1, -1, -229),
    (1, 1, 0, -1, 223),
    (1, 1, 0, 1, 223),
    (0, 1, -2, -1, -220),
    (2, 1, -1, -1, -220),
    (1, 0, 1, 1, -185),
    (2, -1, -2, -1, 181),
    (0, 1, 2, 1, -177),
    (4, 0, -2, -1, 176),
    (4, -1, -1, -1, 166),
    (1, 0, 1, -1, -164),
    (4, 0, 1, -1, 132),
    (1, 0, -1, -1, -119),
    (4, -1, 0, -1, 115),
    (2, -2, 0, 1, 107),
];

/// Reduces an angle in radians to `[0, 2π)`.
#[inline]
pub(crate) fn norm_2pi(x: f64) -> f64 {
    let y = x.rem_euclid(TAU);
    // rem_euclid can round up to exactly 2π for tiny negative inputs
    if y >= TAU { 0.0 } else { y }
}

/// Reduces an angle in radians to `(-π, π]`.
#[inline]
pub(crate) fn norm_pi(x: f64) -> f64 {
    let y = norm_2pi(x);
    if y > PI { y - TAU } else { y }
}

/// Reduces an angle in degrees to `[0, 360)` and converts it to radians.
#[inline]
fn deg_to_rad_2pi(deg: f64) -> f64 {
    norm_2pi(deg.rem_euclid(360.0).to_radians())
}

/// `L′`, eq. 47.1, in degrees (not reduced).
#[inline]
fn mn_long_deg(T: f64) -> f64 {
    Horner_eval!(
        T,
        218.3164477,
        481267.88123421,
        -0.0015786,
        1.0 / 538841.0,
        -1.0 / 65194000.0
    )
}

/// `D`, eq. 47.2, in degrees (not reduced).
#[inline]
fn mn_elong_deg(T: f64) -> f64 {
    Horner_eval!(
        T,
        297.8501921,
        445267.1114034,
        -0.0018819,
        1.0 / 545868.0,
        -1.0 / 113065000.0
    )
}

/// `M`, eq. 47.3, in degrees (not reduced).
#[inline]
fn sun_mn_anom_deg(T: f64) -> f64 {
    Horner_eval!(T, 357.5291092, 35999.0502909, -0.0001536, 1.0 / 24490000.0)
}

/// `M′`, eq. 47.4, in degrees (not reduced).
#[inline]
fn mn_anom_deg(T: f64) -> f64 {
    Horner_eval!(
        T,
        134.9633964,
        477198.8675055,
        0.0087414,
        1.0 / 69699.0,
        -1.0 / 14712000.0
    )
}

/// `F`, eq. 47.5, in degrees (not reduced).
#[inline]
fn arg_of_lat_deg(T: f64) -> f64 {
    Horner_eval!(
        T,
        93.272095,
        483202.0175233,
        -0.0036539,
        -1.0 / 3526000.0,
        1.0 / 863310000.0
    )
}

/// `E`, eq. 47.6.
#[inline]
pub(crate) fn ecc_fac(T: f64) -> f64 {
    1.0 - T * (0.002516 + T * 0.0000074)
}

/// The fundamental arguments for `T` Julian centuries (TD) from J2000.0.
pub(crate) fn fundamental_args_frm_cent(T: f64) -> FundamentalArgs {
    FundamentalArgs {
        mn_long: deg_to_rad_2pi(mn_long_deg(T)),
        mn_elong: deg_to_rad_2pi(mn_elong_deg(T)),
        sun_mn_anom: deg_to_rad_2pi(sun_mn_anom_deg(T)),
        mn_anom: deg_to_rad_2pi(mn_anom_deg(T)),
        arg_of_lat: deg_to_rad_2pi(arg_of_lat_deg(T)),
        a1: deg_to_rad_2pi(119.75 + 131.849 * T),
        a2: deg_to_rad_2pi(53.09 + 479264.29 * T),
        a3: deg_to_rad_2pi(313.45 + 481266.484 * T),
        ecc_fac: ecc_fac(T),
    }
}

/// Computes the fundamental arguments of the lunar theory (Meeus, eqs. 47.1
/// to 47.6, and the arguments `A1`, `A2`, `A3`).
///
/// # Arguments
///
/// * `JDE`: Julian Ephemeris Day (TD)
pub fn fundamental_args(JDE: f64) -> FundamentalArgs {
    fundamental_args_frm_cent(time::julian_cent(JDE))
}

/// Power of `E` that multiplies a term whose argument contains `m` times `M`.
#[inline]
fn ecc_power(m: i8, E: f64) -> f64 {
    match m.abs() {
        1 => E,
        2 => E * E,
        _ => 1.0,
    }
}

/// The periodic sums for the given fundamental arguments.
pub(crate) fn periodic_sums_frm_args(a: &FundamentalArgs) -> PeriodicSums {
    let (D, M, M1, F, E) = (
        a.mn_elong,
        a.sun_mn_anom,
        a.mn_anom,
        a.arg_of_lat,
        a.ecc_fac,
    );
    let L1 = a.mn_long;

    let arg = |d: i8, m: i8, m1: i8, f: i8| {
        f64::from(d) * D + f64::from(m) * M + f64::from(m1) * M1 + f64::from(f) * F
    };

    let mut sum_l = 0.0;
    let mut sum_r = 0.0;
    for &(d, m, m1, f, cl, cr) in TERMS_LR.iter() {
        let x = arg(d, m, m1, f);
        let e = ecc_power(m, E);
        sum_l += f64::from(cl) * e * x.sin();
        sum_r += f64::from(cr) * e * x.cos();
    }

    let mut sum_b = 0.0;
    for &(d, m, m1, f, cb) in TERMS_B.iter() {
        sum_b += f64::from(cb) * ecc_power(m, E) * arg(d, m, m1, f).sin();
    }

    // Additive terms: A1 (Venus), A2 (Jupiter), and the terms in L′ due to
    // the flattening of the Earth.
    sum_l += 3958.0 * a.a1.sin() + 1962.0 * (L1 - F).sin() + 318.0 * a.a2.sin();
    sum_b += -2235.0 * L1.sin()
        + 382.0 * a.a3.sin()
        + 175.0 * (a.a1 - F).sin()
        + 175.0 * (a.a1 + F).sin()
        + 127.0 * (L1 - M1).sin()
        - 115.0 * (L1 + M1).sin();

    PeriodicSums {
        long: (sum_l * 1e-6).to_radians(),
        lat: (sum_b * 1e-6).to_radians(),
        dist: sum_r * 1e-3,
    }
}

/// Computes the sums of the periodic terms of the Moon's longitude,
/// latitude and distance (Meeus, tables 47.A and 47.B, with the additive
/// terms).
///
/// Useful mainly to check intermediate values against the book; use
/// [`geocentric_position`] to get the position.
///
/// # Arguments
///
/// * `JDE`: Julian Ephemeris Day (TD)
pub fn periodic_sums(JDE: f64) -> PeriodicSums {
    periodic_sums_frm_args(&fundamental_args(JDE))
}

/// Computes the equatorial horizontal parallax of the Moon from its
/// distance, `sin π = 6378.14 / Δ` (Meeus, chapter 47).
///
/// Returns `NaN` if `dist` is smaller than the Earth's equatorial radius.
#[inline]
pub(crate) fn eq_hz_parllx_frm_dist(dist: f64) -> f64 {
    (EARTH_EQ_RAD_KM / dist).asin()
}

/// Computes the geometric geocentric position of the Moon, referred to the
/// mean ecliptic and equinox of the date (Meeus, chapter 47).
///
/// The longitude does *not* include nutation; add the nutation in longitude
/// (or use [`apparent_position`]) to get the apparent longitude. The constant
/// part of the light-time correction is already included.
///
/// Accuracy: about 10″ in longitude and 4″ in latitude according to Meeus
/// (abridged ELP-2000/82). The polynomial mean arguments are those of
/// Chapront et al. (1998); the accuracy slowly degrades for dates far from
/// the present (many centuries).
///
/// # Arguments
///
/// * `JDE`: Julian Ephemeris Day (TD)
pub fn geocentric_position(JDE: f64) -> GeocentricPosition {
    let args = fundamental_args(JDE);
    let sums = periodic_sums_frm_args(&args);
    let dist = MEAN_DIST_KM + sums.dist;

    GeocentricPosition {
        long: norm_2pi(args.mn_long + sums.long),
        lat: sums.lat,
        dist,
        eq_hz_parllx: eq_hz_parllx_frm_dist(dist),
    }
}

/// Computes the apparent geocentric position of the Moon, in ecliptic and
/// equatorial coordinates (Meeus, chapter 47 and example 47.a).
///
/// The apparent longitude is the geometric longitude of
/// [`geocentric_position`] plus the nutation in longitude (IAU 1980, from
/// [`nutation::nutation`]). Right ascension and declination are referred to
/// the true equator of the date, using the true obliquity
/// `ε0 + Δε` with `ε0` from the IAU formula ([`ecliptic::mn_oblq_IAU`]).
/// The Moon's aberration is included through the light-time constant of
/// `L′`.
///
/// # Arguments
///
/// * `JDE`: Julian Ephemeris Day (TD)
pub fn apparent_position(JDE: f64) -> ApparentPosition {
    let geo = geocentric_position(JDE);
    let (nut_long, nut_oblq) = nutation::nutation(JDE);
    let true_oblq = ecliptic::mn_oblq_IAU(JDE) + nut_oblq;
    let long = norm_2pi(geo.long + nut_long);

    ApparentPosition {
        long,
        lat: geo.lat,
        dist: geo.dist,
        eq_hz_parllx: geo.eq_hz_parllx,
        asc: norm_2pi(coords::asc_frm_ecl(long, geo.lat, true_oblq)),
        dec: coords::dec_frm_ecl(long, geo.lat, true_oblq),
        nut_long,
        true_oblq,
    }
}

/// `Ω`, eq. 47.7, in radians in `[0, 2π)`, for `T` Julian centuries.
pub(crate) fn mn_ascend_node_frm_cent(T: f64) -> f64 {
    deg_to_rad_2pi(Horner_eval!(
        T,
        125.0445479,
        -1934.1362891,
        0.0020754,
        1.0 / 467441.0,
        -1.0 / 60616000.0
    ))
}

/// True node, in radians in `[0, 2π)`, for `T` Julian centuries.
pub(crate) fn true_ascend_node_frm_cent(T: f64) -> f64 {
    let a = fundamental_args_frm_cent(T);
    let (D, M, M1, F) = (a.mn_elong, a.sun_mn_anom, a.mn_anom, a.arg_of_lat);

    let corr = -1.4979 * (2.0 * (D - F)).sin() - 0.1500 * M.sin() - 0.1226 * (2.0 * D).sin()
        + 0.1176 * (2.0 * F).sin()
        - 0.0801 * (2.0 * (M1 - F)).sin();

    norm_2pi(mn_ascend_node_frm_cent(T) + corr.to_radians())
}

/// `Π`, mean perigee, in radians in `[0, 2π)`, for `T` Julian centuries.
pub(crate) fn mn_perigee_frm_cent(T: f64) -> f64 {
    deg_to_rad_2pi(Horner_eval!(
        T,
        83.3532465,
        4069.0137287,
        -0.01032,
        -1.0 / 80053.0,
        1.0 / 18999000.0
    ))
}

/// Computes the longitude of the mean ascending node of the lunar orbit,
/// `Ω` (Meeus, eq. 47.7).
///
/// The longitude is measured from the mean equinox of the date.
///
/// # Returns
///
/// The longitude *| in radians*, in `[0, 2π)`.
///
/// # Arguments
///
/// * `JDE`: Julian Ephemeris Day (TD)
pub fn mean_ascending_node(JDE: f64) -> f64 {
    mn_ascend_node_frm_cent(time::julian_cent(JDE))
}

/// Computes the longitude of the true ascending node of the lunar orbit,
/// the node of the instantaneous (osculating) orbit.
///
/// This is the mean node of [`mean_ascending_node`] plus the five main
/// periodic terms given by Meeus (chapter 47, from Chapront's tables). The
/// difference from the mean node never exceeds about 2°.
///
/// # Returns
///
/// The longitude *| in radians*, in `[0, 2π)`, measured from the mean
/// equinox of the date.
///
/// # Arguments
///
/// * `JDE`: Julian Ephemeris Day (TD)
pub fn true_ascending_node(JDE: f64) -> f64 {
    true_ascend_node_frm_cent(time::julian_cent(JDE))
}

/// Computes the longitude of the mean perigee of the lunar orbit, `Π`
/// (Meeus, chapter 47).
///
/// # Returns
///
/// The longitude *| in radians*, in `[0, 2π)`, measured from the mean
/// equinox of the date.
///
/// # Arguments
///
/// * `JDE`: Julian Ephemeris Day (TD)
pub fn mean_perigee(JDE: f64) -> f64 {
    mn_perigee_frm_cent(time::julian_cent(JDE))
}

/// Apparent geocentric position of the Sun used by the lunar routines.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ApparentSun {
    /// Apparent ecliptic longitude, radians in `[0, 2π)`
    pub long: f64,
    /// Ecliptic latitude, radians in `(-π, π]`
    pub lat: f64,
    /// Earth–Sun distance, km
    pub dist: f64,
    /// Apparent right ascension, radians in `[0, 2π)`
    pub asc: f64,
    /// Apparent declination, radians
    pub dec: f64,
}

/// Astronomical unit in kilometres (IAU 2012).
pub(crate) const AU_KM: f64 = 149_597_870.7;

/// Computes the apparent position of the Sun from VSOP87 (via
/// [`crate::sun::geocent_ecl_pos`]), converted to the FK5 frame, with the
/// given nutation in longitude and aberration (Meeus, chapter 25).
pub(crate) fn apparent_sun(JDE: f64, nut_long: f64, true_oblq: f64) -> ApparentSun {
    let (geo, R) = crate::sun::geocent_ecl_pos(JDE);
    // sun::geocent_ecl_pos reports the latitude in [0, 2π); bring it back
    // to a signed angle (it is only a fraction of an arcsecond).
    let (long_fk5, lat_fk5) = crate::sun::ecl_coords_to_FK5(JDE, geo.long, norm_pi(geo.lat));
    let long = norm_2pi(long_fk5 + nut_long + crate::aberr::sol_aberr(R));

    ApparentSun {
        long,
        lat: lat_fk5,
        dist: R * AU_KM,
        asc: norm_2pi(coords::asc_frm_ecl(long, lat_fk5, true_oblq)),
        dec: coords::dec_frm_ecl(long, lat_fk5, true_oblq),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn norm_helpers() {
        assert_eq!(norm_2pi(-1e-300), 0.0);
        assert!(norm_2pi(-1e-17) < TAU);
        assert!((norm_2pi(-0.5) - (TAU - 0.5)).abs() < 1e-15);
        assert!((norm_pi(TAU - 0.5) + 0.5).abs() < 1e-15);
        assert_eq!(norm_pi(PI), PI);
        assert!((norm_pi(-PI) - PI).abs() < 1e-15);
    }

    #[test]
    fn ecc_power_by_multiple_of_M() {
        let e = 0.99;
        assert_eq!(ecc_power(0, e), 1.0);
        assert_eq!(ecc_power(1, e), e);
        assert_eq!(ecc_power(-1, e), e);
        assert_eq!(ecc_power(2, e), e * e);
        assert_eq!(ecc_power(-2, e), e * e);
    }

    // Meeus, example 47.a: the Sun's apparent position at 1992 April 12,
    // 0h TD, as quoted in examples 48.a and 53.b.
    #[test]
    fn apparent_sun_1992_apr_12() {
        let jde = 2448724.5;
        let (nut_long, nut_oblq) = nutation::nutation(jde);
        let eps = ecliptic::mn_oblq_IAU(jde) + nut_oblq;
        let s = apparent_sun(jde, nut_long, eps);
        assert!((s.long.to_degrees() - 22.33978).abs() < 2e-5);
        assert!((s.asc.to_degrees() - 20.6579).abs() < 1e-4);
        assert!((s.dec.to_degrees() - 8.6964).abs() < 1e-4);
        assert!((s.dist / AU_KM - 1.0024977).abs() < 1e-7);
    }
}

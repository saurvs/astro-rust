//! Illuminated fraction of the Moon's disk and position angle of its bright
//! limb (Meeus, *Astronomical Algorithms*, 2nd ed., chapter 48).
//!
//! The illuminated fraction `k` follows from the phase angle `i` (the
//! Sun–Moon–Earth angle) as `k = (1 + cos i) / 2`. The phase angle can be
//! computed accurately from the geocentric elongation of the Moon and the
//! distances of the Sun and the Moon ([`phase_angle`]), approximately as the
//! supplement of the elongation ([`phase_angle_frm_elong_approx`]), or
//! without any position at all from the mean arguments of the lunar theory
//! ([`phase_angle_approx`]).
//!
//! [`illumination`] does everything for a given instant.
//!
//! # Example
//!
//! ```
//! use astro::lunar::illumination;
//!
//! // Meeus, example 48.a: 1992 April 12, 0h TD
//! let ill = illumination::illumination(2448724.5);
//! assert!((ill.illum_frac - 0.6786).abs() < 1e-4);
//! assert!((ill.bright_limb.to_degrees() - 285.0).abs() < 0.05);
//! ```

use std::f64::consts::PI;

use super::position::{self, norm_2pi, norm_pi};

/// Illumination of the Moon's disk at an instant, as seen from the centre of
/// the Earth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Illumination {
    /// Geocentric elongation `ψ` of the Moon from the Sun *| in radians*,
    /// in `[0, π]`
    pub geocent_elong: f64,
    /// Phase angle `i`, the selenocentric elongation of the Earth from the
    /// Sun *| in radians*, in `[0, π]`
    pub phase_angle: f64,
    /// Illuminated fraction `k` of the disk, in `[0, 1]`
    pub illum_frac: f64,
    /// Position angle `χ` of the midpoint of the bright limb, measured
    /// eastward from the north point of the disk *| in radians*, in
    /// `[0, 2π)`
    pub bright_limb: f64,
}

/// Computes the geocentric elongation `ψ` of the Moon from the Sun, from
/// equatorial coordinates (Meeus, eq. 48.2, first form).
///
/// # Returns
///
/// The elongation *| in radians*, in `[0, π]`.
///
/// # Arguments
///
/// * `sun_asc`, `sun_dec`  : Right ascension and declination of the Sun
///   *| in radians*
/// * `moon_asc`, `moon_dec`: Right ascension and declination of the Moon
///   *| in radians*
pub fn geocent_elong_frm_eq(sun_asc: f64, sun_dec: f64, moon_asc: f64, moon_dec: f64) -> f64 {
    let c = sun_dec.sin() * moon_dec.sin()
        + sun_dec.cos() * moon_dec.cos() * (sun_asc - moon_asc).cos();
    c.clamp(-1.0, 1.0).acos()
}

/// Computes the geocentric elongation `ψ` of the Moon from the Sun, from
/// ecliptic coordinates (Meeus, eq. 48.2, second form). The Sun's latitude
/// is neglected.
///
/// # Returns
///
/// The elongation *| in radians*, in `[0, π]`.
///
/// # Arguments
///
/// * `moon_long`: Ecliptic longitude of the Moon *| in radians*
/// * `moon_lat` : Ecliptic latitude of the Moon *| in radians*
/// * `sun_long` : Ecliptic longitude of the Sun *| in radians*
pub fn geocent_elong_frm_ecl(moon_long: f64, moon_lat: f64, sun_long: f64) -> f64 {
    (moon_lat.cos() * (moon_long - sun_long).cos())
        .clamp(-1.0, 1.0)
        .acos()
}

/// Computes the phase angle `i` of the Moon from its geocentric elongation
/// (Meeus, eq. 48.3).
///
/// # Returns
///
/// The phase angle *| in radians*, in `[0, π]`.
///
/// # Arguments
///
/// * `geocent_elong`: Geocentric elongation `ψ` of the Moon *| in radians*,
///   in `[0, π]`
/// * `moon_dist`    : Earth–Moon distance `Δ`
/// * `sun_dist`     : Earth–Sun distance `R`, in the same unit as
///   `moon_dist`
pub fn phase_angle(geocent_elong: f64, moon_dist: f64, sun_dist: f64) -> f64 {
    (sun_dist * geocent_elong.sin()).atan2(moon_dist - sun_dist * geocent_elong.cos())
}

/// Computes the phase angle of the Moon approximately as `i = π − ψ`, i.e.
/// with `cos i = −cos ψ` (Meeus, chapter 48).
///
/// The resulting error in the illuminated fraction never exceeds 0.0014.
///
/// # Returns
///
/// The phase angle *| in radians*, in `[0, π]`.
///
/// # Arguments
///
/// * `geocent_elong`: Geocentric elongation `ψ` of the Moon *| in radians*,
///   in `[0, π]`
#[inline]
pub fn phase_angle_frm_elong_approx(geocent_elong: f64) -> f64 {
    PI - geocent_elong
}

/// Computes an approximate phase angle of the Moon from the mean arguments
/// `D`, `M` and `M′` alone (Meeus, eq. 48.4). The Moon's latitude is
/// neglected; no position of the Sun or the Moon is needed.
///
/// # Returns
///
/// The phase angle *| in radians*, in `[0, π]`. (The formula gives a signed
/// angle; its absolute value is returned, which has the same cosine and so
/// the same illuminated fraction.)
///
/// # Arguments
///
/// * `JDE`: Julian Ephemeris Day (TD)
pub fn phase_angle_approx(JDE: f64) -> f64 {
    let a = position::fundamental_args(JDE);
    let (D, M, M1) = (a.mn_elong, a.sun_mn_anom, a.mn_anom);

    let corr = -6.289 * M1.sin() + 2.100 * M.sin()
        - 1.274 * (2.0 * D - M1).sin()
        - 0.658 * (2.0 * D).sin()
        - 0.214 * (2.0 * M1).sin()
        - 0.110 * D.sin();

    norm_pi(PI - D + corr.to_radians()).abs()
}

/// Computes the illuminated fraction `k` of the Moon's disk from the phase
/// angle (Meeus, eq. 48.1).
///
/// `k` is both the ratio of the illuminated area to the area of the whole
/// disk and the illuminated fraction of the diameter perpendicular to the
/// line of cusps.
///
/// # Returns
///
/// The illuminated fraction, in `[0, 1]`.
///
/// # Arguments
///
/// * `phase_angle`: Phase angle `i` *| in radians*
#[inline]
pub fn illum_frac_frm_phase_angle(phase_angle: f64) -> f64 {
    (1.0 + phase_angle.cos()) / 2.0
}

/// Computes the position angle `χ` of the midpoint of the Moon's bright
/// limb (Meeus, eq. 48.5).
///
/// The angle is measured from the north point of the disk towards the east
/// (not from the Moon's axis of rotation, nor from the zenith; subtract the
/// parallactic angle to refer it to the zenith). It is near 270° around
/// First Quarter and near 90° after Full Moon. The cusps are at `χ ± 90°`.
/// The same formula applies to a planet.
///
/// # Returns
///
/// The position angle *| in radians*, in `[0, 2π)`.
///
/// # Arguments
///
/// * `sun_asc`, `sun_dec`  : Apparent right ascension and declination of
///   the Sun *| in radians*
/// * `moon_asc`, `moon_dec`: Apparent right ascension and declination of
///   the Moon *| in radians*
pub fn bright_limb_pos_angl(sun_asc: f64, sun_dec: f64, moon_asc: f64, moon_dec: f64) -> f64 {
    let d_asc = sun_asc - moon_asc;
    let n = sun_dec.cos() * d_asc.sin();
    let d = sun_dec.sin() * moon_dec.cos() - sun_dec.cos() * moon_dec.sin() * d_asc.cos();

    norm_2pi(n.atan2(d))
}

/// Computes the illumination of the Moon's disk for a given instant (Meeus,
/// chapter 48, accurate method).
///
/// The Moon's apparent position comes from
/// [`position::apparent_position`] and the Sun's from VSOP87 (with
/// nutation and aberration). The phase angle is computed with eq. 48.3 from
/// the elongation of eq. 48.2.
///
/// # Arguments
///
/// * `JDE`: Julian Ephemeris Day (TD)
pub fn illumination(JDE: f64) -> Illumination {
    let moon = position::apparent_position(JDE);
    let sun = position::apparent_sun(JDE, moon.nut_long, moon.true_oblq);

    let geocent_elong = geocent_elong_frm_eq(sun.asc, sun.dec, moon.asc, moon.dec);
    let i = phase_angle(geocent_elong, moon.dist, sun.dist);

    Illumination {
        geocent_elong,
        phase_angle: i,
        illum_frac: illum_frac_frm_phase_angle(i),
        bright_limb: bright_limb_pos_angl(sun.asc, sun.dec, moon.asc, moon.dec),
    }
}

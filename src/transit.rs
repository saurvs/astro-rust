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

//! Time of rise, transit and set for a celestial body
//!
//! This module implements the method of Meeus, *Astronomical Algorithms*,
//! 2nd ed., chapter 15. Given the apparent right ascension and declination
//! of a body at 0h Dynamical Time on three consecutive days `D - 1`, `D` and
//! `D + 1`, and the apparent sidereal time at Greenwich at 0h UT on day `D`,
//! it finds the instants (in UT, as fractions of day `D`) at which the body
//! rises, crosses the local meridian (upper culmination) and sets.
//!
//! A rise or a set is the moment the *geometric* altitude of the body's
//! centre equals a "standard altitude" `h0` that allows for refraction at
//! the horizon (and for the Sun and the Moon, for their semidiameter and the
//! Moon's parallax); see [`standard_altitude`]. Twilights are obtained the
//! same way with `h0` = −6°, −12° or −18° for the Sun.
//!
//! Conventions: angles in radians; geographic longitude positive **west** of
//! Greenwich (negative east), as in Meeus; the result is a fraction of the UT
//! day `D` (multiply by 24 for hours).
//!
//! If you would rather not compute the inputs yourself, see
//! [`observer::rise_transit_set`](crate::observer::rise_transit_set), which
//! does everything from a date and a place.
//!
//! # Example
//!
//! Meeus' example 15.a, the rising of Venus at Boston on 1988 March 20:
//!
//! ```
//! use astro::{coords, transit};
//!
//! let deg = |d: f64| d.to_radians();
//! let boston = coords::GeographPoint { long: deg(71.0833), lat: deg(42.3333) };
//! let venus = [
//!     coords::EqPoint { asc: deg(40.68021), dec: deg(18.04761) }, // March 19, 0h TD
//!     coords::EqPoint { asc: deg(41.73129), dec: deg(18.44092) }, // March 20
//!     coords::EqPoint { asc: deg(42.78204), dec: deg(18.82742) }, // March 21
//! ];
//! let h0 = transit::standard_altitude(&transit::TransitBody::StarOrPlanet, 0.0);
//! let m = transit::day_fraction(
//!     &transit::TransitType::Rise, h0, &boston,
//!     &venus[0], &venus[1], &venus[2],
//!     deg(177.74208), // apparent sidereal time at Greenwich, 0h UT
//!     56.0,           // ΔT in seconds
//! )
//! .unwrap();
//! assert!((m - 0.51766).abs() < 1e-5); // 12h25m UT
//! ```

use crate::coords;
use crate::interpol;
use std::f64::consts::{PI, TAU};

const DEG: f64 = PI / 180.0;

/// Rate of the Greenwich sidereal time, in radians per mean solar day
/// (Meeus' 360.985647° per day).
const SIDEREAL_RATE: f64 = 360.985647 * DEG;

/// Stop iterating when the correction to `m` is smaller than this (in days,
/// about 0.01 s).
const M_TOLERANCE: f64 = 1e-7;

/// Maximum number of correction steps.
const MAX_STEPS: usize = 20;

/// The standard altitude of a star or a planet at rise and set, −0°34′, in
/// radians: the mean refraction at the horizon (Meeus, ch. 15).
pub const STD_ALT_STAR_OR_PLANET: f64 = -34.0 / 60.0 * DEG;

/// The standard altitude of the Sun at sunrise and sunset, −0°50′, in
/// radians: refraction at the horizon plus the Sun's semidiameter, so that
/// the times refer to the upper limb (Meeus, ch. 15).
pub const STD_ALT_SUN: f64 = -50.0 / 60.0 * DEG;

/// A mean standard altitude of the Moon, +0.125°, in radians, for when no
/// great accuracy is needed. Prefer [`standard_altitude`] with the Moon's
/// actual equatorial horizontal parallax.
pub const STD_ALT_MOON_MEAN: f64 = 0.125 * DEG;

/// Represents a celestial body in transit
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransitBody {
    /// A star or a planet
    StarOrPlanet,
    /// The Sun
    Sun,
    /// The Moon
    Moon,
}

/// Represents a transit type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransitType {
    /// Rise
    Rise,
    /// Transit across the local meridian (upper culmination)
    Transit,
    /// Set
    Set,
}

/// Why a rise or a set could not be computed
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoEvent {
    /// The body stays above the standard altitude for the whole day
    /// (it is circumpolar, e.g. the midnight Sun).
    AlwaysAbove,
    /// The body stays below the standard altitude for the whole day (e.g.
    /// the polar night).
    AlwaysBelow,
    /// The iteration did not settle on a time. This happens when the body
    /// only just reaches the standard altitude, so that its motion during the
    /// day decides whether it rises or sets at all; it is also returned for
    /// non-finite inputs.
    NoConvergence,
}

impl std::fmt::Display for NoEvent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            NoEvent::AlwaysAbove => "the body stays above the altitude all day",
            NoEvent::AlwaysBelow => "the body stays below the altitude all day",
            NoEvent::NoConvergence => "the time of the event could not be determined",
        })
    }
}

impl std::error::Error for NoEvent {}

/// Computes the standard altitude `h0` of a body at rise and set
///
/// This is the geometric altitude of the body's centre at the moment of its
/// apparent rise or set (Meeus, ch. 15):
///
/// * stars and planets: −0°34′ ([`STD_ALT_STAR_OR_PLANET`]), the mean
///   refraction at the horizon;
/// * the Sun: −0°50′ ([`STD_ALT_SUN`]), refraction plus the semidiameter, so
///   that the times refer to the Sun's upper limb;
/// * the Moon: `0.7275 π − 0°34′`, where `π` is the Moon's equatorial
///   horizontal parallax; this accounts for the Moon's parallax and
///   semidiameter, which vary with its distance.
///
/// # Returns
///
/// * `h0`: Standard altitude *| in radians*
///
/// # Arguments
///
/// * `transit_body`: The `TransitBody`
/// * `moon_eq_hz_parallax`: Equatorial horizontal parallax of the Moon
///   *| in radians*. Only used for `TransitBody::Moon`.
pub fn standard_altitude(transit_body: &TransitBody, moon_eq_hz_parallax: f64) -> f64 {
    match *transit_body {
        TransitBody::StarOrPlanet => STD_ALT_STAR_OR_PLANET,
        TransitBody::Sun => STD_ALT_SUN,
        TransitBody::Moon => 0.7275 * moon_eq_hz_parallax + STD_ALT_STAR_OR_PLANET,
    }
}

/// Computes the hour angle at which a body reaches a given altitude
///
/// This is Meeus' equation 15.1, `cos H0 = (sin h0 − sin φ sin δ) / (cos φ
/// cos δ)`, with the declination held fixed. `H0` is half the time, as an
/// angle, that the body spends above the altitude `h0` (the semidiurnal arc).
///
/// # Returns
///
/// * `Ok(H0)`: The hour angle *| in radians*, in `[0, π]`. The body is at
///   altitude `h0` at the hour angles `−H0` (rising) and `+H0` (setting).
/// * `Err(NoEvent::AlwaysAbove)` or `Err(NoEvent::AlwaysBelow)` if the body
///   never crosses `h0` (the right-hand side of 15.1 is outside `[−1, 1]`).
/// * `Err(NoEvent::NoConvergence)` for non-finite inputs.
///
/// # Arguments
///
/// * `alt`: The altitude `h0` *| in radians*
/// * `observer_lat`: Observer's geographical latitude *| in radians*
/// * `dec`: Declination of the body *| in radians*
pub fn hour_angle_at_altitude(alt: f64, observer_lat: f64, dec: f64) -> Result<f64, NoEvent> {
    let num = alt.sin() - observer_lat.sin() * dec.sin();
    let den = observer_lat.cos() * dec.cos();

    if !(num.is_finite() && den.is_finite()) {
        return Err(NoEvent::NoConvergence);
    }
    // den >= 0 for latitudes and declinations in [-90°, 90°].
    let den = den.abs();
    if num > den {
        return Err(NoEvent::AlwaysBelow);
    }
    if num < -den {
        return Err(NoEvent::AlwaysAbove);
    }
    if den == 0.0 {
        // num == 0 too: the body stays exactly at h0 (at a pole).
        return Err(NoEvent::NoConvergence);
    }

    Ok((num / den).clamp(-1.0, 1.0).acos())
}

/// Computes the time of rise, transit or set of a body, as a fraction of
/// the UT day
///
/// Meeus, ch. 15. The approximate times of equation 15.2 are corrected
/// iteratively, interpolating the body's coordinates between the three given
/// days, until the correction is below 0.01 s.
///
/// # Returns
///
/// * `Ok(m)`: The time of the event on day `D`, in UT, as a fraction of a day
///   (`24 m` is the time in hours). Normally `m` is in `[0, 1)`. For an event
///   within a few minutes of 0h UT the correction may push `m` just outside
///   that range; a negative value, or one of 1 or more, then means that the
///   event nearest to 0h actually falls on the previous or next UT day. Use
///   `m.rem_euclid(1.0)` if you want a time on day `D` regardless.
/// * `Err(NoEvent)`: For a rise or a set, if the body stays above
///   (`AlwaysAbove`) or below (`AlwaysBelow`) the altitude `h0` all day, as
///   judged from its declination on day `D` (Meeus' note 2), or if the
///   iteration does not converge (`NoConvergence`). A transit is always
///   computed, even for a circumpolar body, unless the inputs are not finite.
///
/// # Arguments
///
/// * `transit_type`: A `TransitType`
/// * `h0`: Standard altitude of the body *| in radians*; see
///   [`standard_altitude`]. For twilights use −6° (civil), −12° (nautical) or
///   −18° (astronomical) with the Sun. Ignored for a transit.
/// * `geograph_point`: Geographic point of the observer *| in radians*,
///   with the longitude measured **positive west** of Greenwich
/// * `eq_point1`: Apparent equatorial point of the body at 0h TD on day
///   `D - 1` *| in radians*
/// * `eq_point2`: Same, on day `D`
/// * `eq_point3`: Same, on day `D + 1`
/// * `apprnt_greenwich_sidr`: Apparent sidereal time at Greenwich at 0h UT
///   on day `D` *| in radians*
/// * `delta_t`: ΔT = TD − UT for day `D` *| in seconds*
///
/// The right ascensions may be given in any range; a step across 0h between
/// consecutive days is handled. The three-point interpolation is accurate
/// for the Sun, the planets and the stars; for the Moon, which moves fast,
/// the result is good to about a minute.
pub fn day_fraction(
    transit_type: &TransitType,
    h0: f64,
    geograph_point: &coords::GeographPoint,
    eq_point1: &coords::EqPoint,
    eq_point2: &coords::EqPoint,
    eq_point3: &coords::EqPoint,
    apprnt_greenwich_sidr: f64,
    delta_t: f64,
) -> Result<f64, NoEvent> {
    let path = Path::new(eq_point1, eq_point2, eq_point3);
    let lat = geograph_point.lat;
    let long = geograph_point.long;

    let m0 = (path.asc[1] + long - apprnt_greenwich_sidr) / TAU;
    let m_start = match *transit_type {
        TransitType::Transit => m0,
        TransitType::Rise => m0 - hour_angle_at_altitude(h0, lat, eq_point2.dec)? / TAU,
        TransitType::Set => m0 + hour_angle_at_altitude(h0, lat, eq_point2.dec)? / TAU,
    };
    if !m_start.is_finite() || !delta_t.is_finite() {
        return Err(NoEvent::NoConvergence);
    }
    // Meeus: bring the approximate times into [0, 1).
    let m_start = m_start.rem_euclid(1.0);

    let mut m = m_start;
    for _ in 0..MAX_STEPS {
        let step = correction(
            transit_type,
            m,
            h0,
            lat,
            long,
            &path,
            apprnt_greenwich_sidr,
            delta_t,
        );
        m += step.delta_m;
        // A correction of a quarter of a day means the linearisation broke
        // down (the body only grazes h0).
        if !m.is_finite() || (m - m_start).abs() > 0.25 {
            return Err(NoEvent::NoConvergence);
        }
        if step.delta_m.abs() < M_TOLERANCE {
            return Ok(m);
        }
    }

    Err(NoEvent::NoConvergence)
}

/// Computes the time of transit for a celestial body
///
/// # Returns
///
/// * `(hour, min, sec)`: Time of the event on the day of interest, in UT.
///   If the body does not rise or set on that day (it is circumpolar or
///   stays below the horizon), `(0, 0, NaN)` is returned.
///
/// # Arguments
///
/// * `transit_type`  : A `TransitType`
/// * `transit_body`  : The `TransitBody`
/// * `geograph_point`: Geographic point of the observer *| in radians*,
///   longitude positive west
///
/// Let `JD` be the Julian day of interest,
///
/// * `eq_point1`: Apparent equatorial point of the transit body at 0h TD on
///   `JD - 1` *| in radians*
/// * `eq_point2`: Apparent equatorial point of the transit body at 0h TD on
///   `JD` *| in radians*
/// * `eq_point3`: Apparent equatorial point of the transit body at 0h TD on
///   `JD + 1` *| in radians*
/// * `apprnt_greenwhich_sidr`: Apparent sidereal time at Greenwich at 0h UT
///   on `JD` *| in radians*
/// * `delta_t`: ΔT for `JD` *| in seconds*
/// * `moon_eq_hz_parallax`: Equatorial horizontal parallax of the Moon on `JD`
///   *| in radians*. *Pass a meaningful value here only when*
///   `TransitBody::Moon` *is passed for* `transit_body`.
#[deprecated(
    since = "3.0.0",
    note = "use `transit::day_fraction`, which reports bodies that do not rise or set \
            instead of returning NaN, or `observer::rise_transit_set`"
)]
pub fn time(
    transit_type: &TransitType,
    transit_body: &TransitBody,
    geograph_point: &coords::GeographPoint,
    eq_point1: &coords::EqPoint,
    eq_point2: &coords::EqPoint,
    eq_point3: &coords::EqPoint,
    apprnt_greenwhich_sidr: f64,
    delta_t: f64,
    moon_eq_hz_parallax: f64,
) -> (i64, i64, f64) {
    let h0 = standard_altitude(transit_body, moon_eq_hz_parallax);

    match day_fraction(
        transit_type,
        h0,
        geograph_point,
        eq_point1,
        eq_point2,
        eq_point3,
        apprnt_greenwhich_sidr,
        delta_t,
    ) {
        Ok(m) => {
            let hours = 24.0 * m.rem_euclid(1.0);
            let hour = hours.floor();
            let minutes = (hours - hour) * 60.0;
            let minute = minutes.floor();
            (hour as i64, minute as i64, (minutes - minute) * 60.0)
        }
        Err(_) => (0, 0, f64::NAN),
    }
}

/// The body's coordinates on the three days, with the right ascensions made
/// continuous across 0h.
struct Path {
    asc: [f64; 3],
    dec: [f64; 3],
}

impl Path {
    fn new(p1: &coords::EqPoint, p2: &coords::EqPoint, p3: &coords::EqPoint) -> Path {
        let a2 = p2.asc;
        Path {
            asc: [a2 + wrap_pi(p1.asc - a2), a2, a2 + wrap_pi(p3.asc - a2)],
            dec: [p1.dec, p2.dec, p3.dec],
        }
    }
}

/// One correction step of Meeus' method, with its intermediate quantities
/// (only `delta_m` is used outside the unit tests).
#[derive(Debug)]
#[cfg_attr(not(test), allow(dead_code))]
struct Step {
    /// Sidereal time at Greenwich at the trial instant *| in radians*
    theta0: f64,
    /// Interpolating factor
    n: f64,
    /// Interpolated right ascension *| in radians*
    asc: f64,
    /// Interpolated declination *| in radians*
    dec: f64,
    /// Local hour angle, in `[−π, π)` *| in radians*
    hour_angle: f64,
    /// Altitude (not computed for a transit) *| in radians*
    alt: f64,
    /// Correction to `m`, in days
    delta_m: f64,
}

fn correction(
    transit_type: &TransitType,
    m: f64,
    h0: f64,
    lat: f64,
    long: f64,
    path: &Path,
    theta0_0h: f64,
    delta_t: f64,
) -> Step {
    let theta0 = theta0_0h + SIDEREAL_RATE * m;
    let n = m + delta_t / 86400.0;
    let asc = interpol::three_values(path.asc[0], path.asc[1], path.asc[2], n);
    let dec = interpol::three_values(path.dec[0], path.dec[1], path.dec[2], n);
    let hour_angle = wrap_pi(coords::hr_angl_frm_observer_long(theta0, long, asc));

    let (alt, delta_m) = match *transit_type {
        TransitType::Transit => (f64::NAN, -hour_angle / TAU),
        TransitType::Rise | TransitType::Set => {
            let alt = coords::alt_frm_eq(hour_angle, dec, lat);
            let dm = (alt - h0) / (TAU * dec.cos() * lat.cos() * hour_angle.sin());
            (alt, dm)
        }
    };

    Step {
        theta0,
        n,
        asc,
        dec,
        hour_angle,
        alt,
        delta_m,
    }
}

/// Reduces an angle to `[−π, π)`.
#[inline]
fn wrap_pi(a: f64) -> f64 {
    (a + PI).rem_euclid(TAU) - PI
}

#[cfg(test)]
mod tests {
    use super::*;

    fn deg(d: f64) -> f64 {
        d.to_radians()
    }

    fn venus() -> (coords::GeographPoint, Path) {
        let boston = coords::GeographPoint {
            long: deg(71.0833),
            lat: deg(42.3333),
        };
        let p = |a, d| coords::EqPoint {
            asc: deg(a),
            dec: deg(d),
        };
        let path = Path::new(
            &p(40.68021, 18.04761),
            &p(41.73129, 18.44092),
            &p(42.78204, 18.82742),
        );
        (boston, path)
    }

    fn close(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() <= tol,
            "expected {expected} ± {tol}, got {actual}"
        );
    }

    fn close_deg_mod_360(actual_rad: f64, expected_deg: f64, tol: f64) {
        let d = (actual_rad.to_degrees() - expected_deg).rem_euclid(360.0);
        assert!(
            d.min(360.0 - d) <= tol,
            "expected {expected_deg}° ± {tol}, got {}°",
            actual_rad.to_degrees()
        );
    }

    // Meeus, example 15.a: the intermediate values of the first correction.
    #[test]
    fn example_15_a_intermediate_values() {
        let (boston, path) = venus();
        let theta0 = deg(177.74208);
        // The book rounds h0 to -0.5667°; use the same value to reproduce its
        // digits (the exact -0°34' changes cos H0 by 1e-6).
        let h0 = deg(-0.5667);

        let h_0 = hour_angle_at_altitude(h0, boston.lat, deg(18.44092)).unwrap();
        close(h_0.cos(), -0.3178735, 1e-7);
        close(h_0.to_degrees(), 108.5344, 1e-4);

        // (type, m, θ0, n, α, δ, H, h, Δm)
        let rows = [
            (
                TransitType::Rise,
                0.51817,
                4.79401,
                0.51882,
                42.27648,
                18.64229,
                -108.56577,
                -0.44393,
                -0.00051,
            ),
            (
                TransitType::Transit,
                0.81965,
                113.62397,
                0.82030,
                42.59324,
                f64::NAN,
                -0.05257,
                f64::NAN,
                0.00015,
            ),
            (
                TransitType::Set,
                0.12113,
                221.46827,
                0.12178,
                41.85927,
                18.48835,
                108.52570,
                -0.52711,
                0.00017,
            ),
        ];
        for (tt, m, th, n, a, d, hh, h, dm) in rows {
            let s = correction(&tt, m, h0, boston.lat, boston.long, &path, theta0, 56.0);
            close_deg_mod_360(s.theta0, th, 1e-5);
            close(s.n, n, 1e-5);
            close(s.asc.to_degrees(), a, 1e-5);
            if !d.is_nan() {
                close(s.dec.to_degrees(), d, 1e-5);
                close(s.alt.to_degrees(), h, 1e-5);
            }
            close(s.hour_angle.to_degrees(), hh, 1e-5);
            close(s.delta_m, dm, 1e-5);
        }
    }

    #[test]
    fn right_ascension_is_unwrapped_across_0h() {
        let p = |a: f64| coords::EqPoint {
            asc: deg(a),
            dec: 0.0,
        };
        let path = Path::new(&p(359.0), &p(0.0), &p(1.0));
        close(path.asc[0], deg(-1.0), 1e-12);
        close(path.asc[2], deg(1.0), 1e-12);
        let path = Path::new(&p(358.0), &p(359.0), &p(0.0));
        close(path.asc[2], deg(360.0), 1e-12);
    }

    #[test]
    fn wrap_pi_range() {
        for a in [-7.0, -PI, -1.0, 0.0, 1.0, PI, 7.0, 100.0] {
            let w = wrap_pi(a);
            assert!((-PI..PI).contains(&w));
            close(
                (a - w).rem_euclid(TAU).min(TAU - (a - w).rem_euclid(TAU)),
                0.0,
                1e-12,
            );
        }
    }
}

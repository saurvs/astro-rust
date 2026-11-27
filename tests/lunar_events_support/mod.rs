//! Helpers for the tests of Meeus chapters 49–52 (lunar phases, apsides,
//! node passages and extreme declinations).
//!
//! The results of those chapters are checked against an independent
//! reference built from the chapter 47 lunar series and the VSOP87 Sun of
//! this crate, by locating the events numerically.

#![allow(dead_code)]

use astro::*;

/// Gregorian calendar date to Julian day.
pub fn jd(year: i16, month: time::Month, day: f64) -> f64 {
    time::julian_day(&time::Date {
        year,
        month,
        decimal_day: day,
        cal_type: time::CalType::Gregorian,
    })
}

/// Moon's geocentric ecliptic longitude and latitude (radians, mean
/// equinox of date) and distance (km) from chapter 47.
pub fn moon(jde: f64) -> (f64, f64, f64) {
    let (p, r) = lunar::geocent_ecl_pos(jde);
    (p.long, p.lat, r)
}

/// Sun's apparent longitude minus nutation (radians): VSOP87, converted to
/// FK5 and corrected for aberration (Meeus, ch. 25). Nutation cancels in
/// the Moon-minus-Sun difference used for the phases.
pub fn sun_long_without_nutation(jde: f64) -> f64 {
    let (p, r) = sun::geocent_ecl_pos(jde);
    p.long - (0.09033_f64 / 3600.0).to_radians() - (20.4898_f64 / 3600.0).to_radians() / r
}

/// Moon's apparent geocentric declination (radians).
pub fn moon_apparent_dec(jde: f64) -> f64 {
    let (l, b, _) = moon(jde);
    let (nut_long, nut_oblq) = nutation::nutation(jde);
    let oblq = ecliptic::mn_oblq_laskar(jde) + nut_oblq;
    coords::dec_frm_ecl(l + nut_long, b, oblq)
}

/// Reduces an angle in radians to (-π, π].
pub fn wrap_pi(x: f64) -> f64 {
    let tau = std::f64::consts::TAU;
    let y = x.rem_euclid(tau);
    if y > std::f64::consts::PI { y - tau } else { y }
}

/// Golden-section search for the minimum of a unimodal `f` on [a, b].
pub fn golden_min(mut a: f64, mut b: f64, f: impl Fn(f64) -> f64) -> f64 {
    let g = 0.5 * (5f64.sqrt() - 1.0);
    let mut c = b - g * (b - a);
    let mut d = a + g * (b - a);
    while b - a > 1e-6 {
        if f(c) < f(d) {
            b = d;
        } else {
            a = c;
        }
        c = b - g * (b - a);
        d = a + g * (b - a);
    }
    0.5 * (a + b)
}

/// Bisection for a root of `f` on [a, b], where `f(a)` and `f(b)` have
/// opposite signs.
pub fn bisect(mut a: f64, mut b: f64, f: impl Fn(f64) -> f64) -> f64 {
    let fa = f(a);
    for _ in 0..50 {
        let m = 0.5 * (a + b);
        if (f(m) > 0.0) == (fa > 0.0) {
            a = m;
        } else {
            b = m;
        }
    }
    0.5 * (a + b)
}

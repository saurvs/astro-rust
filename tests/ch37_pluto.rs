//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 37: Pluto.

#![allow(non_snake_case)]

mod common;

use astro::*;
use common::{assert_close, assert_close_deg};

const JDE_37A: f64 = 2448908.5; // 1992 October 13.0 TD
const SUN_37A: (f64, f64, f64) = (-0.9373959, -0.3131679, -0.1357792); // example 26.b

fn jd(year: i16, month: time::Month, day: f64) -> f64 {
    time::julian_day(&time::Date {
        year,
        month,
        decimal_day: day,
        cal_type: time::CalType::Gregorian,
    })
}

// Meeus, example 37.a (1): geometric heliocentric coordinates.
#[test]
fn heliocent_pos_example_37a() {
    let (l, b, r) = pluto::heliocent_pos(JDE_37A);
    assert_close_deg(l.to_degrees(), 232.74071, 6e-6);
    assert_close(b.to_degrees(), 14.58782, 6e-6);
    assert_close(r, 29.711111, 6e-7);
}

// Meeus, example 37.a (2): the position repeated for the light-time
// corrected instant October 12.82368.
#[test]
fn heliocent_pos_example_37a_retarded() {
    let (l, b, r) = pluto::heliocent_pos(JDE_37A - 0.17632);
    assert_close_deg(l.to_degrees(), 232.73949, 6e-6);
    assert_close(b.to_degrees(), 14.58801, 6e-6);
    assert_close(r, 29.711094, 6e-7);
}

// Meeus, example 37.a (2): heliocentric rectangular equatorial coordinates
// from formula 37.1, for both instants of the example. For the first
// instant the book evidently used l and b rounded to 0.00001°, so those
// rounded values are the input here.
#[test]
fn heliocent_rect_eq_coords_example_37a() {
    let (x, y, z) = pluto::heliocent_rect_eq_coords(
        232.74071_f64.to_radians(),
        14.58782_f64.to_radians(),
        29.711111,
    );
    assert_close(x, -17.4079141, 2e-7);
    assert_close(y, -23.9730804, 2e-7);
    assert_close(z, -2.2374228, 2e-7);

    let (l, b, r) = pluto::heliocent_pos(JDE_37A - 0.17632);
    let (x, y, z) = pluto::heliocent_rect_eq_coords(l, b, r);
    assert_close(x, -17.4083780, 2e-7);
    assert_close(y, -23.9727452, 2e-7);
    assert_close(z, -2.2371797, 2e-7);
}

// Meeus, example 37.a (2): geocentric astrometric J2000.0 position.
#[test]
fn geocent_astrometric_pos_example_37a() {
    let (X, Y, Z) = SUN_37A;
    let pos = pluto::geocent_astrometric_pos(JDE_37A, X, Y, Z).unwrap();

    assert_close(pos.true_dist, 30.528746, 6e-7);
    assert_close(pos.dist, 30.528739, 6e-7);
    assert_close(pos.light_time, 0.17632, 6e-6);

    assert_close_deg(pos.asc.to_degrees(), 232.93231, 6e-6);
    assert_close(pos.dec.to_degrees(), -4.45802, 6e-6);

    // 15h31m43.8s, -4°27'29"
    assert_close_deg(
        pos.asc.to_degrees(),
        common::hms(15, 31, 43.8),
        common::arcsec(0.75),
    );
    assert_close(
        pos.dec.to_degrees(),
        common::dms(true, 4, 27, 29.0),
        common::arcsec(0.5),
    );
}

#[test]
fn geocent_astrometric_pos_light_time_consistency() {
    let (X, Y, Z) = SUN_37A;
    let pos = pluto::geocent_astrometric_pos(JDE_37A, X, Y, Z).unwrap();
    // tau is the light-time for the distance at which Pluto is seen.
    assert_close(pos.light_time, 0.0057755183 * pos.dist, 1e-12);
    assert!((0.0..std::f64::consts::TAU).contains(&pos.asc));
}

// With the Sun at the geocentre, the "geocentric" true distance is Pluto's
// heliocentric radius vector and the direction is the heliocentric one,
// converted to the equator of J2000.0.
#[test]
fn geocent_astrometric_pos_sun_at_origin() {
    let jde = 2455000.5;
    let pos = pluto::geocent_astrometric_pos(jde, 0.0, 0.0, 0.0).unwrap();
    let (_, _, r) = pluto::heliocent_pos(jde);
    assert_close(pos.true_dist, r, 1e-12);

    let (l2, b2, r2) = pluto::heliocent_pos(jde - pos.light_time);
    let eps = 0.397777156_f64.asin();
    let asc = coords::asc_frm_ecl(l2, b2, eps);
    let dec = coords::dec_frm_ecl(l2, b2, eps);
    common::assert_close_rad(pos.asc, asc, 1e-8);
    assert_close(pos.dec, dec, 1e-8);
    assert_close(pos.dist, r2, 1e-12);
    // The light-time is that of Pluto's distance to the Sun (about 5.5 hours
    // per 40 AU).
    assert!(
        pos.light_time > 0.15 && pos.light_time < 0.3,
        "{}",
        pos.light_time
    );
}

#[test]
fn heliocent_rect_eq_coords_round_trip() {
    let eps = 0.397777156_f64.asin();
    for &(l, b, r) in &[(0.3, 0.2, 30.0), (4.0, -0.3, 45.0), (2.5, 0.0, 1.0)] {
        let (x, y, z) = pluto::heliocent_rect_eq_coords(l, b, r);
        assert_close((x * x + y * y + z * z).sqrt(), r, 1e-12 * r);
        let asc = y.atan2(x);
        let dec = z.atan2(x.hypot(y));
        common::assert_close_rad(asc, coords::asc_frm_ecl(l, b, eps), 1e-9);
        assert_close(dec, coords::dec_frm_ecl(l, b, eps), 1e-9);
    }
}

#[test]
fn valid_range_constants_match_calendar() {
    assert_eq!(pluto::VALID_FROM_JDE, jd(1885, time::Month::Jan, 1.0));
    assert_eq!(pluto::VALID_UNTIL_JDE, jd(2100, time::Month::Jan, 1.0));
}

#[test]
fn valid_range_boundaries() {
    let start = pluto::VALID_FROM_JDE;
    let end = pluto::VALID_UNTIL_JDE;

    assert!(pluto::in_valid_range(start));
    assert!(!pluto::in_valid_range(start - 1e-6));
    assert!(pluto::in_valid_range(end - 1e-3));
    assert!(!pluto::in_valid_range(end));
    assert!(!pluto::in_valid_range(f64::NAN));
    assert!(!pluto::in_valid_range(f64::INFINITY));

    assert!(pluto::heliocent_pos_checked(start).is_some());
    assert!(pluto::heliocent_pos_checked(start - 0.5).is_none());
    assert!(pluto::heliocent_pos_checked(end - 0.5).is_some());
    assert!(pluto::heliocent_pos_checked(end).is_none());
    // Julian calendar date before 1582 and a negative year.
    assert!(pluto::heliocent_pos_checked(2299160.5).is_none());
    assert!(pluto::heliocent_pos_checked(1000000.5).is_none());
    assert!(pluto::heliocent_pos_checked(f64::NAN).is_none());

    let (X, Y, Z) = SUN_37A;
    assert!(pluto::geocent_astrometric_pos(start - 1.0, X, Y, Z).is_none());
    assert!(pluto::geocent_astrometric_pos(end, X, Y, Z).is_none());
    assert!(pluto::geocent_astrometric_pos(JDE_37A, f64::NAN, Y, Z).is_none());
}

#[test]
fn heliocent_pos_checked_matches_unchecked() {
    let jde = jd(1950, time::Month::June, 15.0);
    assert_eq!(
        pluto::heliocent_pos_checked(jde),
        Some(pluto::heliocent_pos(jde))
    );
}

// The longitude is reduced to [0, 2π) over the whole valid span; before the
// fix it exceeded 2π late in the 21st century.
#[test]
fn heliocent_long_normalised() {
    let (l, _, _) = pluto::heliocent_pos(jd(2099, time::Month::Dec, 1.0));
    assert!((0.0..std::f64::consts::TAU).contains(&l), "{l}");

    let mut jde = pluto::VALID_FROM_JDE;
    while jde < pluto::VALID_UNTIL_JDE {
        let (l, b, r) = pluto::heliocent_pos(jde);
        assert!((0.0..std::f64::consts::TAU).contains(&l), "{jde}: {l}");
        assert!(b.abs() < 18_f64.to_radians(), "{jde}: {b}");
        assert!(r > 29.0 && r < 50.0, "{jde}: {r}");
        jde += 97.3;
    }
}

// Pluto moves at most about 0.04° per day as seen from the Sun; consecutive
// days must therefore agree closely, including across the 0/360° wrap.
#[test]
fn heliocent_pos_continuity() {
    let mut jde = pluto::VALID_FROM_JDE;
    let (mut l0, mut b0, mut r0) = pluto::heliocent_pos(jde);
    while jde < pluto::VALID_UNTIL_JDE - 5.0 {
        jde += 5.0;
        let (l, b, r) = pluto::heliocent_pos(jde);
        assert_close_deg(l.to_degrees(), l0.to_degrees(), 0.2);
        assert_close(b, b0, 0.2_f64.to_radians());
        assert_close(r, r0, 0.05);
        (l0, b0, r0) = (l, b, r);
    }
}

// Pluto's perihelion passage of 1989 September (r ≈ 29.66 AU): the radius
// vector has its minimum then.
#[test]
fn heliocent_perihelion_1989() {
    let r = |jde: f64| pluto::heliocent_pos(jde).2;
    let mut best = (0.0, f64::MAX);
    let mut jde = jd(1985, time::Month::Jan, 1.0);
    while jde < jd(1995, time::Month::Jan, 1.0) {
        if r(jde) < best.1 {
            best = (jde, r(jde));
        }
        jde += 10.0;
    }
    assert!(best.0 > jd(1989, time::Month::Jan, 1.0) && best.0 < jd(1990, time::Month::Jan, 1.0));
    assert_close(best.1, 29.66, 0.02);
}

// Meeus, chapter 37: mean orbital elements near 2000.
#[test]
fn mn_orb_elements_2000AD() {
    let (a, e, i, omega, w) = pluto::mn_orb_elements_2000AD();
    assert_close(a, 39.543, 1e-12);
    assert_close(e, 0.2490, 1e-12);
    assert_close(i.to_degrees(), 17.140, 1e-12);
    assert_close(omega.to_degrees(), 110.307, 1e-12);
    assert_close(w.to_degrees(), 113.768, 1e-12);
}

// Meeus, chapter 41: Pluto's magnitude -1.00 + 5 log(r Δ).
#[test]
fn apprnt_mag_84() {
    assert_close(pluto::apprnt_mag_84(1.0, 1.0), -1.0, 1e-12);
    assert_close(pluto::apprnt_mag_84(10.0, 10.0), 9.0, 1e-12);
    // Example 37.a distances: about magnitude 13.8.
    let m = pluto::apprnt_mag_84(30.528739, 29.711094);
    assert_close(m, -1.0 + 5.0 * (30.528739_f64 * 29.711094).log10(), 1e-12);
    assert_close(m, 13.79, 0.01);
}

#[test]
#[allow(deprecated)]
fn semidiameter() {
    let s = pluto::semidiameter(1.0);
    assert_close(s.to_degrees() * 3600.0, 2.07, 1e-12);
    assert_close(pluto::semidiameter(30.0), s / 30.0, 1e-18);
    assert_eq!(pluto::semdiameter(30.0), pluto::semidiameter(30.0));
}

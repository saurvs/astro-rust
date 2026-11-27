//! Meeus, chapter 12: sidereal time at Greenwich.

mod common;

use astro::time;
use std::f64::consts::TAU;

/// Hours, minutes, seconds of time to radians.
fn hms_rad(h: f64, m: f64, s: f64) -> f64 {
    ((h + m / 60.0 + s / 3600.0) * 15.0).to_radians()
}

/// Seconds of time to radians.
fn sec_rad(s: f64) -> f64 {
    (s / 240.0).to_radians()
}

const JD_12A: f64 = 2446895.5; // 1987 April 10, 0h UT
const JD_12B: f64 = 2446896.30625; // 1987 April 10, 19h21m UT

// Meeus, example 12.a: mean sidereal time at 0h UT
#[test]
fn example_12a_mean() {
    common::assert_close(time::julian_cent(JD_12A), -0.127296372348, 1e-12);
    let expected = hms_rad(13.0, 10.0, 46.3668);
    common::assert_close_rad(time::mn_sidr(JD_12A), expected, sec_rad(1e-4));
    common::assert_close_rad(time::mean_sidereal_time_0h(JD_12A), expected, sec_rad(1e-4));
}

// Meeus, example 12.a: apparent sidereal time, from Δψ = -3.788″ and
// ε = 23°26′36.85″ (example 22.a)
#[test]
fn example_12a_apparent() {
    let nut_in_long = (-3.788_f64 / 3600.0).to_radians();
    let true_oblq = common::dms(false, 23, 26, 36.85).to_radians();

    let eq_eq = time::equation_of_equinoxes(nut_in_long, true_oblq);
    common::assert_close(eq_eq.to_degrees() * 240.0, -0.2317, 5e-5);

    let mean = time::mn_sidr(JD_12A);
    common::assert_close_rad(
        time::apprnt_sidr(mean, nut_in_long, true_oblq),
        hms_rad(13.0, 10.0, 46.1351),
        sec_rad(1e-4),
    );

    // with the crate's own nutation and obliquity
    let expected = hms_rad(13.0, 10.0, 46.1351);
    common::assert_close_rad(
        time::apparent_sidereal_time(JD_12A),
        expected,
        sec_rad(1e-4),
    );
    common::assert_close_rad(astro::apprnt_sidr!(JD_12A), expected, sec_rad(1e-4));
}

// Meeus, example 12.b: mean sidereal time at 19h21m UT
#[test]
fn example_12b_any_instant() {
    common::assert_close(time::julian_cent(JD_12B), -0.12727430, 1e-8);
    let theta = time::mn_sidr(JD_12B);
    // Meeus prints 128.7378734°; exact arithmetic gives 128.73787328°, and
    // the f64 representation of the JD alone is worth about 1e-7°.
    common::assert_close_rad(theta, 128.7378734_f64.to_radians(), 3e-7_f64.to_radians());
    common::assert_close_rad(theta, hms_rad(8.0, 34.0, 57.0896), sec_rad(1e-4));

    // First method: sidereal time at 0h plus the elapsed UT times 1.00273790935
    let ut_seconds = 19.0 * 3600.0 + 21.0 * 60.0;
    common::assert_close(ut_seconds * time::SIDEREAL_RATE, 69850.7228, 1e-4);
    let theta0 = time::mean_sidereal_time_0h(JD_12B);
    common::assert_close_rad(theta0, hms_rad(13.0, 10.0, 46.3668), sec_rad(1e-4));
    common::assert_close_rad(
        theta0 + sec_rad(ut_seconds * time::SIDEREAL_RATE),
        theta,
        sec_rad(1e-4),
    );
}

#[test]
fn sidereal_rate_matches_eq_12_4() {
    // Meeus rounds the ratio to 11 decimals
    common::assert_close(time::SIDEREAL_RATE, 360.98564736629 / 360.0, 1e-12);
}

// Eq. 12.3 at 0h and eq. 12.4 agree, and 12.3 plus the sidereal rate
// reproduces 12.4 during the day. The rounded constants of the two
// expressions drift apart by about 1e-12 degree per day, which reaches
// 0.001 s of time a thousand years from 2000.
#[test]
fn zero_hour_formula_agrees_with_any_instant_formula() {
    let mut jd0: f64 = 2_000_000.5;
    while jd0 < 2_800_000.0 {
        let tol = if (jd0 - 2451545.0).abs() < 40_000.0 {
            sec_rad(1e-4)
        } else {
            1e-7
        };
        let theta0 = time::mean_sidereal_time_0h(jd0);
        common::assert_close_rad(time::mn_sidr(jd0), theta0, tol);
        for frac in [0.1, 0.5, 0.9] {
            let jd = jd0 + frac;
            // any instant of the day gives the 0h value
            assert_eq!(time::mean_sidereal_time_0h(jd), theta0);
            let theta = theta0 + (frac * 360.0 * time::SIDEREAL_RATE).to_radians();
            common::assert_close_rad(time::mn_sidr(jd), theta, tol);
        }
        jd0 += 3_652.0;
    }
}

#[test]
fn results_are_normalised() {
    let mut jd = -1.0e6;
    while jd < 5.0e6 {
        for v in [
            time::mn_sidr(jd),
            time::mean_sidereal_time_0h(jd),
            time::apparent_sidereal_time(jd),
        ] {
            assert!((0.0..TAU).contains(&v), "JD {jd}: {v}");
        }
        jd += 9876.54321;
    }
    // apparent sidereal time wraps around 0
    assert!((0.0..TAU).contains(&time::apprnt_sidr(1e-7, -1e-5, 0.4)));
    common::assert_close(
        time::apprnt_sidr(1e-7, -1e-5, 0.0),
        TAU - (1e-5 - 1e-7),
        1e-12,
    );
    common::assert_close(time::apprnt_sidr(TAU - 1e-7, 1e-5, 0.0), 1e-5 - 1e-7, 1e-12);
}

#[test]
fn local_sidereal_time_longitude_positive_west() {
    let theta0 = time::mn_sidr(JD_12B);
    // Washington, 77°03′56″ W (Meeus, example 13.b)
    let washington = common::dms(false, 77, 3, 56.0).to_radians();
    common::assert_close_rad(
        time::local_sidereal_time(theta0, washington),
        theta0 - washington,
        1e-12,
    );
    // an eastern observer (negative longitude) is ahead of Greenwich
    let east = (-30.0_f64).to_radians();
    common::assert_close_rad(
        time::local_sidereal_time(theta0, east),
        theta0 + 30_f64.to_radians(),
        1e-12,
    );
    // wrap-around and normalisation
    common::assert_close(time::local_sidereal_time(0.1, 0.2), TAU - 0.1, 1e-12);
    common::assert_close(time::local_sidereal_time(6.2, -0.2), 6.4 - TAU, 1e-12);
    for l in [-3.0, -1.0, 0.0, 1.0, 3.0] {
        let v = time::local_sidereal_time(theta0, l);
        assert!((0.0..TAU).contains(&v));
    }
    assert_eq!(time::local_sidereal_time(theta0, 0.0), theta0);
}

// The equation of the equinoxes never exceeds about 1.2 s of time.
#[test]
fn equation_of_equinoxes_is_small() {
    let mut jd = 2_400_000.5;
    while jd < 2_500_000.0 {
        let d = time::apparent_sidereal_time(jd) - time::mn_sidr(jd);
        let d = (d + TAU / 2.0).rem_euclid(TAU) - TAU / 2.0;
        assert!(
            d.abs() < sec_rad(1.2),
            "JD {jd}: {} s",
            d.to_degrees() * 240.0
        );
        jd += 333.3;
    }
}

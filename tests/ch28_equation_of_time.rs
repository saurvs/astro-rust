//! Meeus, chapter 28: equation of time.

mod common;

use astro::*;
use common::*;

use std::f64::consts::PI;

/// 1992 October 13.0 TD, the date of examples 28.a and 28.b.
const JDE_1992_OCT_13: f64 = 2448908.5;

/// Radians of hour angle to minutes of time.
fn minutes(e: f64) -> f64 {
    e.to_degrees() * 4.0
}

// Meeus, example 28.a: L0, then eq. 28.1 with the book's α, Δψ and ε
#[test]
fn example_28a_with_book_inputs() {
    assert_close(time::julian_mill(JDE_1992_OCT_13), -0.0072183436, 1e-10);
    assert_close_deg(
        sun::mean_longitude_vsop87(JDE_1992_OCT_13).to_degrees(),
        201.807193,
        5e-7,
    );

    let e = planet::earth::equation_of_time(
        JDE_1992_OCT_13,
        198.378178_f64.to_radians(),
        arcsec(15.908).to_radians(),
        23.4401443_f64.to_radians(),
    );
    assert_close(e.to_degrees(), 3.427351, 1e-6);
    assert_close(minutes(e), 13.70940, 1e-5);
}

// Meeus, example 28.a with the complete VSOP87: α is 0.2" smaller than the
// book's (Appendix III) value, so E is larger by 0.014 s.
#[test]
fn example_28a_complete_vsop87() {
    let e = sun::equation_of_time(JDE_1992_OCT_13);
    assert_close(minutes(e), 13.70940, 0.0005);
    // 13m42.6s
    assert_close(minutes(e) * 60.0, 13.0 * 60.0 + 42.6, 0.05);
}

// Meeus, example 28.b: Smart's formula
#[test]
fn example_28b_smart() {
    let e = sun::equation_of_time_smart(JDE_1992_OCT_13);
    assert_close(e, 0.059825572, 1e-8);
    assert_close(e.to_degrees(), 3.427753, 1e-6);
    // 13m42.7s
    assert_close(minutes(e) * 60.0, 13.0 * 60.0 + 42.7, 0.05);

    // y = tan²(ε/2) = 0.0430381 with the IAU mean obliquity
    let y = (ecliptic::mn_oblq_IAU(JDE_1992_OCT_13) / 2.0).tan().powi(2);
    assert_close(y, 0.0430381, 5e-8);
}

// Bug: earth::equation_of_time subtracted α and L0 each reduced to one
// turn, so just after the March equinox (α a little above 0°, L0 a little
// below 360°) it returned about +1433 minutes instead of about −7 minutes.
#[test]
fn wraps_near_march_equinox() {
    // 2000 March 21.0 TD
    let jde = 2451624.5;
    let pos = sun::apparent_position(jde);
    assert!(pos.asc.to_degrees() < 1.0);

    let e = planet::earth::equation_of_time(jde, pos.asc, pos.nut_long, pos.true_oblq);
    assert_close(minutes(e), -7.23, 0.01);
    assert_eq!(e, sun::equation_of_time(jde));

    // a contrived case: the result is always reduced to [-π, π)
    let e = planet::earth::equation_of_time(2451545.0, 359.9_f64.to_radians(), 0.0, 0.4);
    assert!((-PI..PI).contains(&e));
}

// Bug: the τ² coefficient of the mean longitude (eq. 28.2) was 0.030328
// instead of 0.03032028; it matters far from J2000.0 (0.25" at τ = −3).
#[test]
fn mean_longitude_tau2_coefficient() {
    let jde = 2451545.0 - 3.0 * 365250.0; // year −1000
    let tau: f64 = -3.0;
    let l0 = 280.4664567 + 360007.6982779 * tau + 0.03032028 * tau.powi(2) + tau.powi(3) / 49931.0
        - tau.powi(4) / 15300.0
        - tau.powi(5) / 2000000.0;

    // with α = Δψ = 0, E = L0 − 0.0057183°
    let e = planet::earth::equation_of_time(jde, 0.0, 0.0, 0.4);
    assert_close_rad(e, (l0 - 0.0057183).to_radians(), 1e-9);
    assert_close_rad(sun::mean_longitude_vsop87(jde), l0.to_radians(), 1e-9);
}

// Through the year 2000: |E| stays below 17 minutes, with a minimum of
// about −14.2 min in mid-February, a maximum of about +16.4 min early in
// November, zeros near April 15, June 13, September 1 and December 25, and
// no jumps. Smart's formula stays within a few seconds of the full method.
#[test]
fn annual_curve_2000() {
    let start = 2451544.5; // 2000 January 1.0
    let values: Vec<f64> = (0..367)
        .map(|d| minutes(sun::equation_of_time(start + d as f64)))
        .collect();

    let (imin, emin) = values
        .iter()
        .enumerate()
        .fold((0, f64::MAX), |a, (i, &v)| if v < a.1 { (i, v) } else { a });
    let (imax, emax) = values
        .iter()
        .enumerate()
        .fold((0, f64::MIN), |a, (i, &v)| if v > a.1 { (i, v) } else { a });

    assert!((-14.5..-14.0).contains(&emin), "minimum {emin}");
    assert!((40..45).contains(&imin), "minimum on day {imin}"); // ~Feb 11
    assert!((16.2..16.6).contains(&emax), "maximum {emax}");
    assert!((305..310).contains(&imax), "maximum on day {imax}"); // ~Nov 3

    let mut zeros = Vec::new();
    for d in 1..values.len() {
        assert!((values[d] - values[d - 1]).abs() < 0.6, "jump on day {d}");
        if values[d - 1].signum() != values[d].signum() {
            zeros.push(d);
        }
    }
    // day of year of the sign changes: Apr 15, Jun 13, Sep 1, Dec 25
    let expected = [105, 164, 244, 359];
    assert_eq!(zeros.len(), 4, "zeros at {zeros:?}");
    for (z, e) in zeros.iter().zip(expected) {
        assert!((*z as i32 - e).abs() <= 2, "zero on day {z}, expected ~{e}");
    }

    for (d, full) in values.iter().enumerate() {
        let smart = minutes(sun::equation_of_time_smart(start + d as f64));
        assert_close(smart * 60.0, full * 60.0, 4.0);
    }
}

// The equation of time stays below 20 minutes in absolute value over many
// centuries (Meeus, chapter 28), and Smart's formula follows the full
// method.
#[test]
fn bounded_over_centuries() {
    for year in [-2000, -1000, 0, 1000, 1600, 2100, 3000, 4000, 5000] {
        let jd0 = 2451545.0 + (year as f64 - 2000.0) * 365.25;
        for i in 0..52 {
            let jde = jd0 + 7.0 * i as f64;
            let full = minutes(sun::equation_of_time(jde));
            let smart = minutes(sun::equation_of_time_smart(jde));
            assert!(full.abs() < 20.0, "{year}: {full} min");
            assert!(smart.abs() < 20.0, "{year}: {smart} min");
            assert_close(smart, full, 0.25);
        }
    }
}

// Meeus, Astronomical Algorithms, 2nd ed., chapter 35: near-parabolic motion

#![allow(non_snake_case)]

mod common;

use astro::consts;
use astro::orbit::{elliptic, near_parabolic, parabolic};
use common::{assert_close, assert_close_deg, assert_close_rad};
use std::f64::consts::{PI, TAU};

fn landgraf(t: f64, e: f64, q: f64) -> Result<(f64, f64), &'static str> {
    near_parabolic::true_anom_and_rad_vec(t, 0.0, e, q, 1e-10)
}

/// Reference solution for an elliptic orbit, through Kepler's equation
fn kepler_elliptic(t: f64, e: f64, q: f64) -> (f64, f64) {
    let a = elliptic::semimaj_axis(q, e);
    let M = elliptic::mn_motion(a) * t;
    let E = elliptic::ecc_anom_checked(M, e).unwrap();
    (
        elliptic::true_anom(E, e).rem_euclid(TAU),
        elliptic::rad_vec_frm_ecc_anom(E, a, e),
    )
}

/// Reference solution for a hyperbolic orbit, through the hyperbolic form
/// of Kepler's equation, e sinh(H) - H = M, solved by bisection
fn kepler_hyperbolic(t: f64, e: f64, q: f64) -> (f64, f64) {
    let a = q / (e - 1.0);
    let M = consts::GAUSS_GRAV / a.powf(1.5) * t;
    let (mut lo, mut hi) = (-50.0_f64, 50.0_f64);
    for _ in 0..300 {
        let mid = 0.5 * (lo + hi);
        if e * mid.sinh() - mid < M {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let H = 0.5 * (lo + hi);
    let v = 2.0 * (((e + 1.0) / (e - 1.0)).sqrt() * (H / 2.0).tanh()).atan();
    (v.rem_euclid(TAU), a * (e * H.cosh() - 1.0))
}

#[test]
fn landgraf_exercises() {
    // Meeus, chapter 35: the table of cases following Landgraf's program
    let cases = [
        (0.921326, 1.0, 138.4783, 102.74426, 2.364192),
        (0.1, 0.987, 254.9, 164.50029, 4.063777),
        (0.123456, 0.99997, -30.47, 221.91190, 0.965053),
        (3.363943, 1.05731, 1237.1, 109.40598, 10.668551),
        (0.5871018, 0.9672746, 20.0, 52.85331, 0.729116),
        (0.5871018, 0.9672746, 0.0, 0.0, 0.5871018),
    ];
    for (q, e, t, v, r) in cases {
        // with the book's accuracy, and with a tighter one
        for acc in [1e-9, 1e-12] {
            let (tv, tr) = near_parabolic::true_anom_and_rad_vec(t, 0.0, e, q, acc).unwrap();
            assert_close(tv.to_degrees(), v, 1e-5);
            assert_close(tr, r, 1e-6);
        }
        // t and T only enter through t - T
        let T = 2451545.0;
        let (tv, tr) = near_parabolic::true_anom_and_rad_vec(T + t, T, e, q, 1e-12).unwrap();
        assert_close(tv.to_degrees(), v, 1e-5);
        assert_close(tr, r, 1e-6);
    }
}

#[test]
fn landgraf_convergence_limits() {
    // Meeus, chapter 35: the table of true anomalies and of the cases that
    // do not converge. The book rounds v to 1° (0.01° for the last block).
    // The limits depend a little on the accuracy asked for; use the book's.
    let v_deg = |q: f64, e: f64, t: f64| {
        near_parabolic::true_anom_and_rad_vec(t, 0.0, e, q, 1e-9).map(|(v, _)| v.to_degrees())
    };

    for (t, v) in [(10.0, 126.0), (20.0, 142.0)] {
        assert_close(v_deg(0.1, 0.9, t).unwrap(), v, 0.5);
    }
    assert_eq!(v_deg(0.1, 0.9, 30.0), Err(near_parabolic::NO_CONVERGENCE));

    let cases = [
        (10.0, 123.0),
        (20.0, 137.0),
        (30.0, 143.0),
        (60.0, 152.0),
        (100.0, 157.0),
        (200.0, 163.0),
        (400.0, 167.0),
    ];
    for (t, v) in cases {
        assert_close(v_deg(0.1, 0.987, t).unwrap(), v, 0.5);
    }
    assert_eq!(
        v_deg(0.1, 0.987, 500.0),
        Err(near_parabolic::NO_CONVERGENCE)
    );

    for (t, v) in [
        (100.0, 156.0),
        (200.0, 161.0),
        (500.0, 166.0),
        (1000.0, 169.0),
        (5000.0, 174.0),
    ] {
        assert_close(v_deg(0.1, 0.999, t).unwrap(), v, 0.5);
    }

    assert_close(v_deg(1.0, 0.99999, 1e5).unwrap(), 172.5, 0.05);
    for (t, v) in [(1e7, 178.41), (1.4e7, 178.58), (1.7e7, 178.68)] {
        assert_close(v_deg(1.0, 0.99999, t).unwrap(), v, 0.005);
    }
    // 17 million days after the perihelion the body is 7220 AU from the Sun
    assert_close(landgraf(1.7e7, 0.99999, 1.0).unwrap().1, 7220.0, 0.5);
    assert_eq!(
        v_deg(1.0, 0.99999, 1.8e7),
        Err(near_parabolic::NO_CONVERGENCE)
    );
}

#[test]
fn divergence_threshold_matches_landgraf() {
    // Bug fix: the size above which a term of the series counts as
    // divergent was 1000 instead of Landgraf's 10000, so these cases
    // returned "no convergence" although the method works (checked against
    // Kepler's equation).
    for (q, e, t) in [(0.1, 0.999, 10000.0), (0.5, 0.999, 1e5), (2.0, 0.999, -1e6)] {
        let (v, r) = landgraf(t, e, q).unwrap();
        let (kv, kr) = kepler_elliptic(t, e, q);
        assert_close_rad(v, kv, 1e-9);
        assert_close(r, kr, 1e-9 * kr);
    }
    let (v, r) = landgraf(10000.0, 0.999, 0.1).unwrap();
    assert_close(v, 3.06229984031017, 1e-9);
    assert_close(r, 48.2980165523957, 1e-8);
}

#[test]
fn agrees_with_kepler_for_elliptic_orbits() {
    for &(q, e, tmax) in &[
        (0.5871018, 0.9672746, 40.0),
        (0.1, 0.987, 400.0),
        (0.123456, 0.99997, 3000.0),
        (1.0, 0.9999, 1e5),
    ] {
        for i in -20..=20 {
            let t = tmax * i as f64 / 20.0;
            let (v, r) = landgraf(t, e, q).unwrap();
            let (kv, kr) = kepler_elliptic(t, e, q);
            assert_close_rad(v, kv, 1e-9);
            assert_close(r, kr, 1e-9 * kr);
        }
    }
}

#[test]
fn agrees_with_kepler_for_hyperbolic_orbits() {
    for &(q, e, tmax) in &[
        (3.363943, 1.05731, 2000.0),
        (0.5, 1.001, 3000.0),
        (1.2, 1.1, 300.0),
    ] {
        for i in -20..=20 {
            let t = tmax * i as f64 / 20.0;
            let (v, r) = landgraf(t, e, q).unwrap();
            let (kv, kr) = kepler_hyperbolic(t, e, q);
            assert_close_rad(v, kv, 1e-9);
            assert_close(r, kr, 1e-9 * kr);
            // the true anomaly stays inside the asymptotes
            let v_signed = if v > PI { v - TAU } else { v };
            assert!(v_signed.abs() < (-1.0 / e).acos());
        }
    }
    // Meeus's hyperbolic case, through the hyperbolic Kepler equation
    let (v, r) = kepler_hyperbolic(1237.1, 1.05731, 3.363943);
    assert_close(v.to_degrees(), 109.40598, 1e-5);
    assert_close(r, 10.668551, 1e-6);
}

#[test]
fn parabolic_case_is_barker() {
    let q = 1.487469;
    for t in [-1e5, -365.25, -1.0, 1e-3, 112.5642, 1e4] {
        let (v, r) = landgraf(t, 1.0, q).unwrap();
        let (pv, pr) = parabolic::true_anom_and_rad_vec(t, 0.0, q);
        assert_close_rad(v, pv, 1e-14);
        assert_close(r, pr, 1e-14 * pr);
    }
}

#[test]
fn before_and_after_perihelion() {
    // v(-t) = 2π - v(t) (the result is in [0, 2π)) and r(-t) = r(t)
    for &(q, e) in &[(0.5871018, 0.9672746), (1.0, 1.0), (3.363943, 1.05731)] {
        for t in [0.5, 5.0, 20.0, 100.0] {
            let (v1, r1) = landgraf(t, e, q).unwrap();
            let (v2, r2) = landgraf(-t, e, q).unwrap();
            assert!((0.0..PI).contains(&v1));
            assert!((PI..TAU).contains(&v2));
            assert_close(v2, TAU - v1, 1e-12);
            assert_close(r2, r1, 1e-12 * r1);
        }
    }
    // at the perihelion
    assert_eq!(landgraf(0.0, 0.97, 0.6), Ok((0.0, 0.6)));
    // continuity across the perihelion
    let (v, r) = landgraf(1e-6, 0.97, 0.6).unwrap();
    assert!(v < 1e-6 && (r - 0.6).abs() < 1e-12);
    let (v, _) = landgraf(-1e-6, 0.97, 0.6).unwrap();
    assert_close_deg(v.to_degrees(), 360.0, 1e-4);
}

#[test]
fn invalid_input_and_termination() {
    let f = near_parabolic::true_anom_and_rad_vec;
    let invalid = Err(near_parabolic::INVALID_INPUT);
    assert_eq!(f(10.0, 0.0, 1.0, 0.0, 1e-9), invalid);
    assert_eq!(f(10.0, 0.0, 1.0, -1.0, 1e-9), invalid);
    assert_eq!(f(10.0, 0.0, -0.5, 1.0, 1e-9), invalid);
    assert_eq!(f(f64::NAN, 0.0, 0.99, 1.0, 1e-9), invalid);
    assert_eq!(f(10.0, 0.0, f64::INFINITY, 1.0, 1e-9), invalid);
    assert_eq!(f(10.0, 0.0, 0.99, 1.0, f64::NAN), invalid);
    // an unreachable accuracy ends (with a result or an error) instead of
    // looping: the inner Newton iteration used to have no limit
    for t in [1.0, 10.0, 100.0, -50.0] {
        let _ = f(t, 0.0, 0.99, 1.0, 0.0);
    }
    assert_eq!(
        f(30.0, 0.0, 0.98, 0.1, -1.0),
        Err(near_parabolic::NO_CONVERGENCE)
    );
    // the error messages are unchanged from earlier releases
    assert_eq!(
        near_parabolic::NO_CONVERGENCE,
        "No convergence at orbit::near_parabolic::true_anom_and_rad_vec()"
    );
}

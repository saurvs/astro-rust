//! Meeus, *Astronomical Algorithms*, 2nd ed., ch. 43: ephemeris for
//! physical observations of Jupiter.
//!
//! The intermediate steps of examples 43.a and 43.b, computed from the
//! book's own inputs, are checked by the unit tests in
//! `src/planet/jupiter/mod.rs`. These tests use the public API, with the
//! crate's VSOP87 positions.

#![allow(non_snake_case)]

mod common;

use astro::*;
use common::{arcsec, assert_close, assert_close_deg};
use std::f64::consts::{PI, TAU};

/// JDE of 1992 December 16, 0h UT, with ΔT = 59 s (examples 43.a, 43.b).
const JDE_43A: f64 = 2_448_972.500_68;

fn in_two_pi(a: f64) -> bool {
    (0.0..TAU).contains(&a)
}

/// Signed difference of two angles in degrees, in [-180°, 180°).
fn ang_diff_deg(a: f64, b: f64) -> f64 {
    ((a - b + PI).rem_euclid(TAU) - PI).to_degrees()
}

/// Every 7.3 days for the given span of Julian years from 2000.
fn dates(from_year: f64, to_year: f64) -> impl Iterator<Item = f64> {
    let start = 2_451_545.0 + 365.25 * from_year;
    let n = (365.25 * (to_year - from_year) / 7.3) as usize;
    (0..n).map(move |i| start + 7.3 * i as f64)
}

// Meeus, example 43.a
#[test]
fn example_43a_physical_ephemeris() {
    let e = planet::jupiter::physical_ephemeris(JDE_43A);

    assert_close_deg(e.earth_dec.to_degrees(), -2.48, 0.005);
    assert_close_deg(e.sun_dec.to_degrees(), -2.20, 0.005);
    assert_close_deg(e.axis_pos_angle.to_degrees(), 24.80, 0.005);
    assert_close_deg(e.central_meridian_i.to_degrees(), 267.63, 0.005);
    assert_close_deg(e.central_meridian_ii.to_degrees(), 72.31, 0.005);
    assert_close(e.phase_corr.to_degrees(), 0.43, 0.005);
    // the book adds the rounded 0.43° to the rounded ω, hence 0.01°
    assert_close_deg(e.illum_central_meridian_i.to_degrees(), 268.06, 0.01);
    assert_close_deg(e.illum_central_meridian_ii.to_degrees(), 72.74, 0.01);
    // step 7
    assert_close(e.earth_dist, 5.661_123_9, 2e-7);
    assert_close(e.sun_dist, 5.446_423_2, 2e-7);
}

// Meeus, example 43.a, through the older API with the book's obliquity
// (step 8) and nutation (step 15).
//
// Before the fix the old code applied the light-time twice (it both took
// Jupiter at t - τ and subtracted the light-time correction from l), which
// put ω about 0.0026° too far: System I came out as 268.0658° (book 268.06,
// geometric 267.63 + 0.43; correct value 268.0632°).
#[test]
fn example_43a_old_api() {
    let e = planet::jupiter::ephemeris(
        JDE_43A,
        23.440_206_9_f64.to_radians(),
        arcsec(16.86).to_radians(),
        arcsec(-1.79).to_radians(),
    );

    assert_close_deg(e.De.to_degrees(), -2.48, 0.005);
    assert_close_deg(e.Ds.to_degrees(), -2.20, 0.005);
    assert_close(e.P.to_degrees(), 24.80, 0.005);
    assert_close(e.w1.to_degrees(), 268.0632, 0.0005);
    assert_close(e.w2.to_degrees(), 72.7355, 0.0005);
}

// Meeus, example 43.b (lower accuracy)
#[test]
fn example_43b_low_accuracy() {
    let e = planet::jupiter::physical_ephemeris_low_accuracy(JDE_43A);

    assert_close_deg(e.central_meridian_i.to_degrees(), 267.69, 0.005);
    assert_close_deg(e.central_meridian_ii.to_degrees(), 72.36, 0.005);
    assert_close(e.phase_corr.to_degrees(), 0.43, 0.005);
    assert_close_deg(e.illum_central_meridian_i.to_degrees(), 268.12, 0.01);
    assert_close_deg(e.illum_central_meridian_ii.to_degrees(), 72.79, 0.01);
    assert_close(e.phase_angle.to_degrees(), -9.924, 5e-4);
    assert_close_deg(e.helio_long.to_degrees(), -178.11, 0.005);
    assert!(in_two_pi(e.helio_long));
    assert_close(e.sun_dec.to_degrees(), -2.194, 5e-4);
    assert_close(e.earth_dec.to_degrees(), -2.50, 0.005);
    assert_close(e.earth_sun_dist, 0.984_13, 5e-6);
    assert_close(e.sun_dist, 5.448_24, 5e-6);
    // the book rounds K before computing Δ (5.66151); unrounded, 5.661497
    assert_close(e.earth_dist, 5.661_51, 2e-5);
    // and the book's own comparison with the accurate method
    let acc = planet::jupiter::physical_ephemeris(JDE_43A);
    assert!(ang_diff_deg(e.central_meridian_i, acc.central_meridian_i).abs() < 0.07);
    assert!(ang_diff_deg(e.central_meridian_ii, acc.central_meridian_ii).abs() < 0.07);
}

// The two accurate entry points share the algorithm: given the same
// obliquity and nutation they must give identical results.
#[test]
fn old_and_new_api_agree() {
    for jde in [JDE_43A, 2_451_545.0, 2_460_000.25, 2_430_000.75] {
        let (nut_long, nut_oblq) = nutation::nutation(jde);
        let old = planet::jupiter::ephemeris(jde, ecliptic::mn_oblq_IAU(jde), nut_long, nut_oblq);
        let new = planet::jupiter::physical_ephemeris(jde);
        assert_eq!(old.De, new.earth_dec);
        assert_eq!(old.Ds, new.sun_dec);
        assert_eq!(old.P, new.axis_pos_angle);
        assert_eq!(old.w1, new.illum_central_meridian_i);
        assert_eq!(old.w2, new.illum_central_meridian_ii);
    }
}

// Ranges and the limits quoted in the chapter: |D_S| ≤ 3.12°, |D_E| ≤ 3.4°
// (we allow a few hundredths), |C| ≤ 0.61°. Before the fix the old API
// could return ω outside [0, 360°) after adding C.
#[test]
fn ranges_and_limits() {
    let mut max = [0.0_f64; 3];
    for jde in dates(-100.0, 100.0) {
        let e = planet::jupiter::physical_ephemeris(jde);
        for a in [
            e.axis_pos_angle,
            e.central_meridian_i,
            e.central_meridian_ii,
            e.illum_central_meridian_i,
            e.illum_central_meridian_ii,
        ] {
            assert!(in_two_pi(a), "angle {a} out of range at {jde}");
        }
        max[0] = max[0].max(e.sun_dec.abs().to_degrees());
        max[1] = max[1].max(e.earth_dec.abs().to_degrees());
        max[2] = max[2].max(e.phase_corr.abs().to_degrees());

        let lo = planet::jupiter::physical_ephemeris_low_accuracy(jde);
        for a in [
            lo.central_meridian_i,
            lo.central_meridian_ii,
            lo.illum_central_meridian_i,
            lo.illum_central_meridian_ii,
            lo.helio_long,
        ] {
            assert!(in_two_pi(a), "low-accuracy angle {a} out of range at {jde}");
        }
        assert!(lo.phase_angle.abs().to_degrees() < 12.0);

        let old = planet::jupiter::ephemeris(jde, ecliptic::mn_oblq_IAU(jde), 0.0, 0.0);
        assert!(in_two_pi(old.P) && in_two_pi(old.w1) && in_two_pi(old.w2));
    }
    assert!(max[0] > 3.10 && max[0] < 3.13, "max |D_S| = {}", max[0]);
    assert!(max[1] > 3.3 && max[1] < 3.45, "max |D_E| = {}", max[1]);
    assert!(max[2] > 0.60 && max[2] < 0.611, "max |C| = {}", max[2]);
}

// The correction for phase is sin²(i/2) in radians, i being the phase
// angle, with the sign of sin(l - l0).
#[test]
fn phase_correction_sign_and_size() {
    for jde in dates(-20.0, 20.0) {
        let e = planet::jupiter::physical_ephemeris(jde);
        let (l0, _, R) = planet::heliocent_coords(&planet::Planet::Earth, jde);
        let (l, _, _) = planet::heliocent_coords(&planet::Planet::Jupiter, jde);
        let (r, delta) = (e.sun_dist, e.earth_dist);
        let cos_i = (r * r + delta * delta - R * R) / (2.0 * r * delta);
        assert_close(e.phase_corr.abs(), (1.0 - cos_i) / 2.0, 1e-12);
        if (l - l0).sin().abs() > 1e-3 {
            assert_eq!(e.phase_corr.signum(), (l - l0).sin().signum(), "at {jde}");
        }
        assert_close(
            ang_diff_deg(e.illum_central_meridian_i, e.central_meridian_i),
            e.phase_corr.to_degrees(),
            1e-9,
        );
    }
}

// Systems I and II differ only through their rotation rates and epochs:
// ω1 - ω2 = 0.872° + 7.63° d - 0.04407° Δ, with d days from 1950 Jan 1.0.
#[test]
fn system_i_minus_system_ii() {
    for jde in dates(-50.0, 50.0).step_by(5) {
        let e = planet::jupiter::physical_ephemeris(jde);
        let d = jde - 2_433_282.5;
        let expected = 0.872 + 7.63 * d - 0.044_07 * e.earth_dist;
        let got = ang_diff_deg(e.central_meridian_i, e.central_meridian_ii);
        assert_close_deg(got, expected, 1e-6);
    }
}

// The central meridians advance by about 877.9°/day (System I) and
// 870.3°/day (System II), i.e. 36.58° and 36.26° per hour, less a small
// amount from the changing geometry.
#[test]
fn central_meridian_rates() {
    for jde in dates(-10.0, 10.0).step_by(20) {
        let a = planet::jupiter::physical_ephemeris(jde);
        let b = planet::jupiter::physical_ephemeris(jde + 1.0 / 24.0);
        let r1 = (b.central_meridian_i - a.central_meridian_i)
            .rem_euclid(TAU)
            .to_degrees();
        let r2 = (b.central_meridian_ii - a.central_meridian_ii)
            .rem_euclid(TAU)
            .to_degrees();
        assert!(
            (36.55..36.60).contains(&r1),
            "System I rate {r1}°/h at {jde}"
        );
        assert!(
            (36.23..36.28).contains(&r2),
            "System II rate {r2}°/h at {jde}"
        );
    }
}

// The lower-accuracy method agrees with the accurate one to 0.1-0.2° in ω
// near the present, as the book says, and to about 0.33° a century away;
// the declinations agree to a few hundredths of a degree.
#[test]
fn low_accuracy_agrees_with_accurate_method() {
    let check = |from: f64, to: f64, tol_omega: f64| {
        for jde in dates(from, to) {
            let a = planet::jupiter::physical_ephemeris(jde);
            let b = planet::jupiter::physical_ephemeris_low_accuracy(jde);
            let d1 = ang_diff_deg(a.central_meridian_i, b.central_meridian_i);
            let d2 = ang_diff_deg(a.central_meridian_ii, b.central_meridian_ii);
            assert!(d1.abs() < tol_omega && d2.abs() < tol_omega, "ω at {jde}");
            assert!(
                (a.sun_dec - b.sun_dec).abs().to_degrees() < 0.012,
                "D_S at {jde}"
            );
            assert!(
                (a.earth_dec - b.earth_dec).abs().to_degrees() < 0.05,
                "D_E at {jde}"
            );
            assert!((a.phase_corr - b.phase_corr).abs().to_degrees() < 0.003);
            assert!((a.earth_dist - b.earth_dist).abs() < 0.008, "Δ at {jde}");
        }
    };
    check(-50.0, 0.0, 0.22);
    check(-100.0, 100.0, 0.35);
}

// D_E, D_S and P change smoothly from day to day.
#[test]
fn smooth_day_to_day() {
    let mut prev = planet::jupiter::physical_ephemeris(2_448_000.5);
    for i in 1..900 {
        let jde = 2_448_000.5 + i as f64;
        let e = planet::jupiter::physical_ephemeris(jde);
        assert!((e.earth_dec - prev.earth_dec).abs().to_degrees() < 0.05);
        assert!((e.sun_dec - prev.sun_dec).abs().to_degrees() < 0.01);
        assert!(ang_diff_deg(e.axis_pos_angle, prev.axis_pos_angle).abs() < 0.1);
        prev = e;
    }
}

// Far from 2000 (before the Gregorian reform, negative years, and the
// fourth millennium) the results remain finite and in range. The pole of
// step 1 is linear in time, so by -1000 it has drifted enough for |D_S| to
// reach 3.6° (physically it stays below 3.12°); hence the loose bounds.
#[test]
fn distant_dates() {
    // -1000 Jan 1, 1000 Jan 1, 1582 Oct 4 (Julian), 3000 Jan 1
    for jde in [1_355_807.5, 2_086_307.5, 2_299_159.5, 2_816_787.5] {
        let e = planet::jupiter::physical_ephemeris(jde);
        assert!(in_two_pi(e.axis_pos_angle) && in_two_pi(e.illum_central_meridian_i));
        assert!(e.sun_dec.abs().to_degrees() < 4.0 && e.earth_dec.abs().to_degrees() < 4.5);
        let lo = planet::jupiter::physical_ephemeris_low_accuracy(jde);
        assert!(in_two_pi(lo.central_meridian_i) && in_two_pi(lo.central_meridian_ii));
        assert!(lo.earth_dec.is_finite() && lo.sun_dec.is_finite());
    }
}

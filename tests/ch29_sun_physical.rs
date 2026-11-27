//! Meeus, chapter 29: ephemeris for physical observations of the Sun.

mod common;

use astro::*;
use common::*;

use std::f64::consts::{PI, TAU};

fn wrap_pi(a: f64) -> f64 {
    (a + PI).rem_euclid(TAU) - PI
}

/// Meeus, example 29.a: 1992 October 13 at 0h UT with ΔT = 59 s.
const JDE_29A: f64 = 2448908.50068;

// Meeus, example 29.a, with the book's apparent longitudes and obliquity
#[test]
fn example_29a_with_book_inputs() {
    // λ = L + 180° − 20.4898"/R with L = 19.908045°, R = 0.997608
    let ab = aberr::sol_aberr(0.997608).to_degrees();
    assert_close(ab, -0.005705, 5e-7);
    let lambda = 19.908045 + 180.0 + ab;
    assert_close(lambda, 199.902340, 5e-7);
    let lambda_nut = lambda + 0.004419;
    assert_close(lambda_nut, 199.906759, 5e-7);

    let (p, b0, l0) = sun::ephemeris(
        JDE_29A,
        lambda.to_radians(),
        lambda_nut.to_radians(),
        23.440144_f64.to_radians(),
    );

    // P = x + y = 22.1790° + 4.0945°, sin B0 = 0.104324, L0 = −121.3683°
    assert_close(p.to_degrees(), 22.1790 + 4.0945, 1e-4);
    assert_close(b0.sin(), 0.104324, 1e-6);
    assert_close_deg(l0.to_degrees(), 360.0 - 121.3683, 1e-4);
    assert!((0.0..TAU).contains(&l0));

    // rounded as in the book: 26.27°, 5.99°, 238.63°
    assert_close(p.to_degrees(), 26.27, 0.005);
    assert_close(b0.to_degrees(), 5.99, 0.005);
    assert_close(l0.to_degrees(), 238.63, 0.005);
}

// Meeus, example 29.a from the time only, with the complete VSOP87
#[test]
fn example_29a_physical_ephemeris() {
    let eph = sun::physical_ephemeris(JDE_29A);
    assert_close(eph.axis_position_angle.to_degrees(), 26.27, 0.005);
    assert_close(eph.center_helio_lat.to_degrees(), 5.99, 0.005);
    assert_close(eph.center_helio_long.to_degrees(), 238.63, 0.005);
}

// Meeus, example 29.b: the beginning of synodic rotation No. 1699
#[test]
fn example_29b_synodic_rot() {
    let jde = sun::synodic_rot(1699);
    assert_close(jde, 2444480.8455 - 0.1225, 1e-4);

    // ΔT = 51 s in 1980: JD 2444480.7224, i.e. 1980 August 29.22 UT
    let jd = jde - 51.0 / 86400.0;
    assert_close(jd, 2444480.7224, 1e-4);
    let (y, m, d) = time::date_frm_julian_day(jd).unwrap();
    assert_eq!((y, m), (1980, 8));
    assert_close(d, 29.22, 0.005);
}

// The exact start of a rotation (L0 = 0) is within 0.002 day of eqs. 29.1
// and 29.2 for every rotation between 1853 and 2100 (Meeus, chapter 29).
#[test]
fn carrington_rot_start_vs_approximation() {
    for c in (1..=3300).step_by(7) {
        let exact = sun::carrington_rot_start(c);
        let approx = sun::synodic_rot(c);
        assert!(
            (exact - approx).abs() < 0.002,
            "rotation {c}: {exact} vs {approx}"
        );

        let l0 = sun::physical_ephemeris(exact).center_helio_long;
        assert!(wrap_pi(l0).abs() < 1e-7, "rotation {c}: L0 = {l0}");
    }
}

// Rotation No. 1 began on 1853 November 9; consecutive rotations are one
// synodic period apart.
#[test]
fn carrington_rotation_1() {
    let jde = sun::carrington_rot_start(1);
    let (y, m, d) = time::date_frm_julian_day(jde).unwrap();
    assert_eq!((y, m, d.floor() as i64), (1853, 11, 9));

    for c in [1, 1000, 2000, 2300] {
        let period = sun::carrington_rot_start(c + 1) - sun::carrington_rot_start(c);
        assert_close(period, 27.2753, 0.2);
    }
}

// L0 decreases by about 13.2° per day; just before the start of a rotation
// it is close to 0°, and just after close to 360°.
#[test]
fn l0_decreases() {
    let start = sun::carrington_rot_start(2000);
    let before = sun::physical_ephemeris(start - 0.01).center_helio_long;
    let after = sun::physical_ephemeris(start + 0.01).center_helio_long;
    assert_close(before.to_degrees(), 0.132, 0.01);
    assert_close(after.to_degrees(), 360.0 - 0.132, 0.01);

    for i in 0..100 {
        let jde = 2451545.0 + 3.3 * i as f64;
        let a = sun::physical_ephemeris(jde).center_helio_long;
        let b = sun::physical_ephemeris(jde + 1.0).center_helio_long;
        assert!((0.0..TAU).contains(&a));
        assert_close(wrap_pi(a - b).to_degrees(), 13.2, 0.1);
    }
}

// Meeus, chapter 29: over a year P stays within ±26.4°, reaching about
// −26.3° near April 7 and +26.3° near October 11, and is zero near January
// 5 and July 7; B0 stays within ±7.26°, reaching −7.25° near March 6 and
// +7.25° near September 8, and is zero near June 6 and December 7.
#[test]
fn annual_extremes_2000() {
    let start = 2451544.5; // 2000 January 1.0
    let eph: Vec<_> = (0..366)
        .map(|d| sun::physical_ephemeris(start + d as f64))
        .collect();
    let p: Vec<f64> = eph
        .iter()
        .map(|e| e.axis_position_angle.to_degrees())
        .collect();
    let b: Vec<f64> = eph
        .iter()
        .map(|e| e.center_helio_lat.to_degrees())
        .collect();

    let argmin = |v: &[f64]| (0..v.len()).fold(0, |i, j| if v[j] < v[i] { j } else { i });
    let argmax = |v: &[f64]| (0..v.len()).fold(0, |i, j| if v[j] > v[i] { j } else { i });

    // day of year (0 = January 1)
    let near = |day: usize, expected: usize, tol: usize| {
        assert!(
            day.abs_diff(expected) <= tol,
            "day {day}, expected ~{expected}"
        );
    };

    assert!(p.iter().all(|x| x.abs() < 26.4));
    assert_close(p[argmin(&p)], -26.3, 0.1);
    near(argmin(&p), 97, 3); // April 7 (leap year)
    assert_close(p[argmax(&p)], 26.3, 0.1);
    near(argmax(&p), 284, 3); // October 11

    assert!(b.iter().all(|x| x.abs() < 7.26));
    assert_close(b[argmin(&b)], -7.25, 0.01);
    near(argmin(&b), 65, 3); // March 6
    assert_close(b[argmax(&b)], 7.25, 0.01);
    near(argmax(&b), 251, 3); // September 8

    let zeros = |v: &[f64]| -> Vec<usize> {
        (1..v.len())
            .filter(|&i| v[i - 1].signum() != v[i].signum())
            .collect()
    };
    let pz = zeros(&p);
    assert_eq!(pz.len(), 2, "{pz:?}");
    near(pz[0], 4, 3); // January 5
    near(pz[1], 188, 3); // July 7
    let bz = zeros(&b);
    assert_eq!(bz.len(), 2, "{bz:?}");
    near(bz[0], 157, 3); // June 6
    near(bz[1], 341, 3); // December 7
}

// The physical ephemeris is defined for any date; the angles stay in their
// ranges far from the present too: |P| <= ε + I, since |x| <= ε and |y| <= I.
#[test]
fn ranges_far_from_present() {
    for jde in [1355807.5, 2086307.5, 2299159.5, 2816787.5] {
        for d in 0..40 {
            let t = jde + 9.0 * d as f64;
            let e = sun::physical_ephemeris(t);
            let eps = sun::apparent_position(t).true_oblq;
            assert!(e.axis_position_angle.abs() <= eps + 7.25_f64.to_radians());
            assert!(e.center_helio_lat.abs() < 7.26_f64.to_radians());
            assert!((0.0..TAU).contains(&e.center_helio_long));
        }
    }
}

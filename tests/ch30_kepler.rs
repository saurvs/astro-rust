// Meeus, Astronomical Algorithms, 2nd ed., chapter 30: the equation of Kepler

#![allow(non_snake_case)]

mod common;

use astro::orbit::elliptic::{self, NewtonCorrLimit};
use common::assert_close;
use std::f64::consts::{PI, TAU};

fn rad(deg: f64) -> f64 {
    deg.to_radians()
}

/// Residual of Kepler's equation
fn residual(E: f64, M: f64, e: f64) -> f64 {
    E - e * E.sin() - M
}

/// A reference solution by plain bisection, independent of the crate
fn reference_ecc_anom(M: f64, e: f64) -> f64 {
    let k = (M / TAU).round();
    let m = M - k * TAU;
    let (sign, m) = if m < 0.0 { (-1.0, -m) } else { (1.0, m) };
    let (mut lo, mut hi) = (0.0_f64, PI);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid - e * mid.sin() < m {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    sign * 0.5 * (lo + hi) + k * TAU
}

// --- Example 30.a: first method (simple iteration) ---------------------------

#[test]
fn example_30a_simple_iteration_steps() {
    // Meeus, example 30.a: e = 0.1, M = 5°; the successive values in degrees
    let (M, e) = (rad(5.0), 0.1);
    let expected = [
        5.499366, 5.549093, 5.554042, 5.554535, 5.554584, 5.554589, 5.554589,
    ];
    let mut E = M;
    for x in expected {
        E = elliptic::kepler_simple_step(M, e, E);
        assert_close(E.to_degrees(), x, 1e-6);
    }
}

#[test]
fn example_30a_simple_iteration_solver() {
    // Meeus, example 30.a, to an accuracy of 0.000001°
    let E = elliptic::ecc_anom_simple_iter(rad(5.0), 0.1, rad(1e-6), 100).unwrap();
    assert_close(E.to_degrees(), 5.554589, 1e-6);
    // Table 30.A: 6 iterations suffice for e = 0.1, M = 5° (plus one to see
    // that the value no longer changes)
    assert!(elliptic::ecc_anom_simple_iter(rad(5.0), 0.1, rad(1e-6), 7).is_some());
}

#[test]
fn simple_iteration_is_slow_for_high_eccentricity() {
    // Meeus, chapter 30: e = 0.99, M = 2°. After 50 iterations the result is
    // 32.345452°, still more than 0.01° from the solution 32.361007°.
    let (M, e) = (rad(2.0), 0.99);
    let first = [2.000000, 3.979598, 5.936635, 7.866758, 9.763644, 11.619294];
    let mut E = M;
    assert_close(E.to_degrees(), first[0], 1e-6);
    for x in &first[1..] {
        E = elliptic::kepler_simple_step(M, e, E);
        assert_close(E.to_degrees(), *x, 1e-6);
    }
    let mut E = M;
    for _ in 0..50 {
        E = elliptic::kepler_simple_step(M, e, E);
    }
    assert_close(E.to_degrees(), 32.345452, 1e-6);
    assert!((E.to_degrees() - 32.361007).abs() > 0.01);

    assert_eq!(elliptic::ecc_anom_simple_iter(M, e, rad(1e-6), 50), None);
    let E = elliptic::ecc_anom_simple_iter(M, e, 1e-15, 10_000).unwrap();
    assert_close(E.to_degrees(), 32.361007, 1e-6);
}

#[test]
fn simple_iteration_oscillates_beyond_m0() {
    // Meeus, chapter 30: e = 0.75, M = 70°, between M0 and 180°, so the
    // successive values alternate above and below the solution
    let (M, e) = (rad(70.0), 0.75);
    let expected = [
        110.380316, 110.281870, 110.307524, 110.300850, 110.302587, 110.302135,
    ];
    let mut E = M;
    for x in expected {
        E = elliptic::kepler_simple_step(M, e, E);
        assert_close(E.to_degrees(), x, 1e-6);
    }
    // The "valley": for M0 = π/2 - e the solution is exactly E = 90°
    let E = elliptic::ecc_anom_checked(PI / 2.0 - e, e).unwrap();
    assert_close(E, PI / 2.0, 1e-15);
}

// --- Example 30.b: second method (Newton) ------------------------------------

#[test]
fn example_30b_newton_steps() {
    // Meeus, example 30.b: corrections and new values, in degrees
    let (M, e) = (rad(5.0), 0.1);
    let steps = [
        (0.554616193, 5.554616193),
        (-0.000026939, 5.554589254),
        (-0.000000001, 5.554589253),
    ];
    let mut E = M;
    for (corr, new) in steps {
        let c = elliptic::kepler_newton_corr(M, e, E);
        assert_close(c.to_degrees(), corr, 1e-9);
        E += c;
        assert_close(E.to_degrees(), new, 1e-9);
    }
    let E = elliptic::ecc_anom_newton(M, e, NewtonCorrLimit::Unlimited, rad(1e-9), 3).unwrap();
    assert_close(E.to_degrees(), 5.554589253, 1e-9);
}

#[test]
fn table_30b_newton_high_eccentricity() {
    // Meeus, table 30.B: Newton's formula starting from E = M
    let newton_seq = |M: f64, e: f64, n: usize| {
        let mut E = M;
        (0..n)
            .map(|_| {
                E += elliptic::kepler_newton_corr(M, e, E);
                E.to_degrees()
            })
            .collect::<Vec<_>>()
    };

    // e = 0.99, M = 2°: converges after 8 iterations
    let seq = newton_seq(rad(2.0), 0.99, 10);
    let book = [
        188.700250865,
        90.0043959725,
        58.7251974236,
        41.762008288,
        34.1821261793,
        32.4485414136,
        32.361223124,
        32.3610074734,
        32.3610074722,
        32.3610074722,
    ];
    for (x, b) in seq.iter().zip(book) {
        // the book was computed with 12 significant digits
        assert_close(*x, b, 1e-8);
    }

    // e = 0.999, M = 6°: the first iterations jump wildly. The book's
    // values were computed with 12 digits and the sequence is chaotic, so
    // only the first values can be compared.
    let seq = newton_seq(rad(6.0), 0.999, 3);
    assert_close(seq[0], 930.362114752, 1e-7);
    assert_close(seq[1], 418.384869795, 1e-6);
    assert_close(seq[2], -345.064633754, 1e-5);
    let E =
        elliptic::ecc_anom_newton(rad(6.0), 0.999, NewtonCorrLimit::Unlimited, 1e-14, 100).unwrap();
    assert_close(E.to_degrees(), 49.5696248539, 1e-9);

    // e = 0.999, M = 7°
    let seq = newton_seq(rad(7.0), 0.999, 3);
    assert_close(seq[0], 832.86912333, 1e-7);
    assert_close(seq[1], 275.954959759, 1e-5);
    assert_close(seq[2], -87.610596019, 1e-4);
    let E = elliptic::ecc_anom_checked(rad(7.0), 0.999).unwrap();
    assert_close(E.to_degrees(), 52.2702615, 1e-7);
}

#[test]
fn newton_with_arcsin_of_sin() {
    // Meeus, chapter 30 (Leingärtner's trick): e = 0.99, M = 0.2 rad
    let (M, e) = (0.2, 0.99);
    assert_close(elliptic::kepler_newton_corr(M, e, M), 6.614719035698, 1e-11);

    let steps = [
        (0.331533728518, 0.531533728518),
        (1.161431415069, 1.692965143587),
        (-0.455401365518, 1.237563778069),
        (-0.150884433942, 1.086679344127),
        (-0.019368331549, 1.067311012578),
        (-0.000313565645, 1.066997446933),
        (-0.000000081651, 1.066997365282),
    ];
    let mut E = M;
    for (corr, new) in steps {
        let c = elliptic::kepler_newton_corr(M, e, E).sin().asin();
        assert_close(c, corr, 1e-11);
        E += c;
        assert_close(E, new, 1e-11);
    }

    // The library applies the same limit: one step with an infinite
    // tolerance stops after the first correction
    let E1 =
        elliptic::ecc_anom_newton(M, e, NewtonCorrLimit::ArcsinOfSin, f64::INFINITY, 1).unwrap();
    assert_close(E1, 0.531533728518, 1e-11);
    let E = elliptic::ecc_anom_newton(M, e, NewtonCorrLimit::ArcsinOfSin, 1e-14, 20).unwrap();
    assert_close(E, 1.066997365282, 1e-12);
}

#[test]
fn newton_with_clamped_correction() {
    // Meeus, chapter 30 (Steele's trick, limit 0.5): e = 0.99, M = 0.2 rad
    let (M, e) = (0.2, 0.99);
    let steps = [
        (6.614719035698, 0.5, 0.7),
        (0.567429870979, 0.5, 1.2),
        (-0.120513681086, -0.120513681086, 1.079486318914),
        (-0.012361504682, -0.012361504682, 1.067124814232),
        (-0.000127435465, -0.000127435465, 1.066997378767),
        (-0.000000013485, -0.000000013485, 1.066997365282),
    ];
    let mut E = M;
    for (raw, used, new) in steps {
        let w = elliptic::kepler_newton_corr(M, e, E);
        assert_close(w, raw, 1e-11);
        E += if w.abs() > 0.5 {
            0.5_f64.copysign(w)
        } else {
            w
        };
        assert_close(E, used + (E - used), 1e-15);
        assert_close(E, new, 1e-11);
    }

    let lim = NewtonCorrLimit::Clamp(0.5);
    let E1 = elliptic::ecc_anom_newton(M, e, lim, f64::INFINITY, 1).unwrap();
    assert_close(E1, 0.7, 1e-15);
    let E = elliptic::ecc_anom_newton(M, e, lim, 1e-14, 20).unwrap();
    assert_close(E, 1.066997365282, 1e-12);
    // Steele notes that a limit of 1 works too
    let E = elliptic::ecc_anom_newton(M, e, NewtonCorrLimit::Clamp(1.0), 1e-14, 30).unwrap();
    assert_close(E, 1.066997365282, 1e-12);
}

#[test]
fn newton_rejects_bad_input() {
    let n = |lim| elliptic::ecc_anom_newton(1.0, 0.5, lim, 1e-12, 50);
    assert_eq!(n(NewtonCorrLimit::Clamp(0.0)), None);
    assert_eq!(n(NewtonCorrLimit::Clamp(-1.0)), None);
    assert_eq!(n(NewtonCorrLimit::Clamp(f64::NAN)), None);
    assert!(n(NewtonCorrLimit::Clamp(0.5)).is_some());
    assert_eq!(
        elliptic::ecc_anom_newton(1.0, 1.0, NewtonCorrLimit::Unlimited, 1e-12, 50),
        None
    );
    assert_eq!(
        elliptic::ecc_anom_newton(f64::NAN, 0.5, NewtonCorrLimit::Unlimited, 1e-12, 50),
        None
    );
    // not converged within max_iter
    assert_eq!(
        elliptic::ecc_anom_newton(rad(7.0), 0.999, NewtonCorrLimit::Unlimited, 1e-12, 10),
        None
    );
}

// --- Third method: binary search ---------------------------------------------

#[test]
fn binary_search_matches_book_values() {
    let bs = |M: f64, e: f64| elliptic::ecc_anom_binary_search(M, e).unwrap();
    assert_close(bs(rad(5.0), 0.1).to_degrees(), 5.554589253, 1e-9);
    assert_close(bs(rad(2.0), 0.99).to_degrees(), 32.3610074722, 1e-9);
    assert_close(bs(rad(6.0), 0.999).to_degrees(), 49.5696248539, 1e-9);
    assert_close(bs(rad(7.0), 0.999).to_degrees(), 52.2702615, 1e-7);
    assert_close(bs(0.2, 0.99), 1.066997365282, 1e-12);
    // Table 30.A
    assert_close(bs(rad(5.0), 0.5).to_degrees(), 9.950063, 1e-6);
    assert_close(bs(rad(1.0), 0.99).to_degrees(), 24.725822, 1e-6);
    assert_close(bs(rad(33.0), 0.99).to_degrees(), 89.722155, 1e-6);
}

#[test]
fn binary_search_quadrants_and_revolutions() {
    for &e in &[0.0, 0.3, 0.9, 0.999999] {
        for i in -40..=40 {
            let M = i as f64 * 0.37;
            let E = elliptic::ecc_anom_binary_search(M, e).unwrap();
            assert_close(E, reference_ecc_anom(M, e), 1e-14 * M.abs().max(1.0));
        }
    }
    assert_eq!(elliptic::ecc_anom_binary_search(1.0, -0.1), None);
    assert_eq!(elliptic::ecc_anom_binary_search(1.0, 1.0), None);
    assert_eq!(elliptic::ecc_anom_binary_search(f64::INFINITY, 0.5), None);
}

// --- Fourth method: approximate formula for small e --------------------------

#[test]
fn low_eccentricity_formula() {
    // Meeus, chapter 30, formula 30.8 for the data of example 30.a
    let (M, e) = (rad(5.0), 0.1);
    assert_close(M.sin() / (M.cos() - e), 0.09725090, 1e-8);
    let E = elliptic::ecc_anom_low_ecc(M, e);
    assert_close(E.to_degrees(), 5.554599, 1e-6);
    // error 0.035″ here, but about 35″ for M = 82°
    let exact = elliptic::ecc_anom_checked(M, e).unwrap();
    assert_close((E - exact).to_degrees() * 3600.0, 0.035, 0.002);
    let M = rad(82.0);
    let err = (elliptic::ecc_anom_low_ecc(M, e) - elliptic::ecc_anom_checked(M, e).unwrap())
        .to_degrees()
        .abs()
        * 3600.0;
    assert_close(err, 35.0, 0.5);
}

#[test]
fn low_eccentricity_greatest_errors() {
    // Meeus, chapter 30: greatest error of formula 30.8
    let max_err = |e: f64| {
        (0..=3600)
            .map(|i| {
                let M = rad(i as f64 * 0.05);
                (elliptic::ecc_anom_low_ecc(M, e) - elliptic::ecc_anom_checked(M, e).unwrap())
                    .abs()
                    .to_degrees()
            })
            .fold(0.0, f64::max)
    };
    assert_close(max_err(0.15), 0.0327, 0.0001);
    assert_close(max_err(0.20), 0.0783, 0.0001);
    assert_close(max_err(0.25), 0.1552, 0.0001);
    assert_close(max_err(0.50), 1.42, 0.005);
    // The book gives 24.7° for e = 0.99, but the error of formula 30.8
    // reaches 43.1° (near M = 2.7°); already at M = 1° it is 35.8°. The
    // other values above are reproduced exactly, so 24.7° seems to be a
    // misprint or a coarser search.
    assert_close(max_err(0.99), 43.11, 0.01);
    // Earth's orbit: less than 0.2″
    assert!(max_err(0.0167) * 3600.0 < 0.2);
}

#[test]
fn low_eccentricity_quadrants() {
    // the approximation stays in the revolution of M, in every quadrant
    for i in -50..=50 {
        let M = i as f64 * 0.3;
        let E = elliptic::ecc_anom_low_ecc(M, 0.05);
        assert!((E - reference_ecc_anom(M, 0.05)).abs() < 1e-3, "M = {M}");
    }
}

// --- The robust solver -------------------------------------------------------

#[test]
fn ecc_anom_example_30a() {
    // Meeus, example 30.a, via the original function
    let E = elliptic::ecc_anom(rad(5.0), 0.1, 1e-7);
    assert_close(E.to_degrees(), 5.554589253, 1e-9);
}

#[test]
fn table_30a_values() {
    // Meeus, table 30.A
    let cases = [
        (0.1, 5.0, 5.554589),
        (0.2, 5.0, 6.246908),
        (0.3, 5.0, 7.134960),
        (0.4, 5.0, 8.313903),
        (0.5, 5.0, 9.950063),
        (0.6, 5.0, 12.356653),
        (0.7, 5.0, 16.167990),
        (0.8, 5.0, 22.656579),
        (0.9, 5.0, 33.344447),
        (0.99, 5.0, 45.361023),
        (0.99, 1.0, 24.725822),
        (0.99, 33.0, 89.722155),
    ];
    for (e, M, E) in cases {
        let x = elliptic::ecc_anom_checked(rad(M), e).unwrap();
        assert_close(x.to_degrees(), E, 1e-6);
        assert_close(elliptic::ecc_anom(rad(M), e, 1e-6).to_degrees(), E, 1e-6);
    }
}

#[test]
fn ecc_anom_honours_accuracy_for_high_eccentricity() {
    // Bug fix: the old fixed-point solver stopped when two successive values
    // differed by less than `accuracy`, and returned 0.1707883 here, wrong
    // by 6.3e-5 rad although 1e-6 was requested.
    let (M, e) = (0.001, 0.999);
    let E = elliptic::ecc_anom(M, e, 1e-6);
    assert_close(E, 0.17085095632357872, 1e-14);
    assert!(residual(E, M, e).abs() < 1e-16);
}

#[test]
fn ecc_anom_terminates_with_zero_accuracy() {
    // Bug fix: with accuracy = 0 the old fixed-point solver could cycle
    // between two floating-point values forever (this case hung).
    let (M, e) = (2.438552371830389, 0.9202082482294778);
    let E = elliptic::ecc_anom(M, e, 0.0);
    assert_close(E, reference_ecc_anom(M, e), 1e-14);
}

#[test]
fn ecc_anom_invalid_input() {
    assert!(elliptic::ecc_anom(1.0, 1.0, 1e-9).is_nan());
    assert!(elliptic::ecc_anom(1.0, 1.5, 1e-9).is_nan());
    assert!(elliptic::ecc_anom(1.0, -0.1, 1e-9).is_nan());
    assert!(elliptic::ecc_anom(f64::NAN, 0.5, 1e-9).is_nan());
    assert_eq!(elliptic::ecc_anom_checked(1.0, f64::NAN), None);
    assert_eq!(elliptic::ecc_anom_checked(f64::INFINITY, 0.5), None);
    assert_eq!(elliptic::ecc_anom_checked(1.0, 1.0), None);
    assert_eq!(elliptic::ecc_anom_simple_iter(1.0, 1.0, 1e-9, 100), None);
    assert_eq!(elliptic::ecc_anom_checked(1.25, 0.0), Some(1.25));
}

#[test]
fn ecc_anom_converges_everywhere() {
    let eccs = [
        0.0,
        1e-12,
        0.0167,
        0.1,
        0.5,
        0.9,
        0.99,
        0.999,
        0.999999,
        1.0 - 1e-12,
    ];
    for &e in &eccs {
        for i in -2000..=2000 {
            // covers several revolutions, both signs, and points near 0 and π
            let M = i as f64 * 0.00873 + if i % 7 == 0 { 1e-9 } else { 0.0 };
            let E = elliptic::ecc_anom_checked(M, e).unwrap();
            let scale = M.abs().max(1.0);
            assert!(
                residual(E, M, e).abs() <= 8.0 * f64::EPSILON * scale,
                "e = {e}, M = {M}: residual {}",
                residual(E, M, e)
            );
            assert!((E - M).abs() <= e + 1e-15 * scale);
            assert_close(E, reference_ecc_anom(M, e), 1e-9 * scale);
        }
    }
}

#[test]
fn ecc_anom_near_parabolic_small_mean_anomaly() {
    // e close to 1 and M close to 0: E ≈ (6M)^(1/3)
    for &(M, e) in &[(1e-6, 0.999999), (1e-9, 1.0 - 1e-12), (1e-12, 0.9999)] {
        let E = elliptic::ecc_anom_checked(M, e).unwrap();
        assert_close(E, reference_ecc_anom(M, e), 1e-13);
        assert!(residual(E, M, e).abs() <= 4.0 * f64::EPSILON * E);
    }
}

#[test]
fn ecc_anom_symmetries() {
    for &e in &[0.2, 0.7, 0.97] {
        for i in 0..100 {
            let M = i as f64 * 0.0631;
            let E = elliptic::ecc_anom_checked(M, e).unwrap();
            // odd function of M
            assert_close(elliptic::ecc_anom_checked(-M, e).unwrap(), -E, 1e-15);
            // one more revolution in M gives one more revolution in E
            assert_close(
                elliptic::ecc_anom_checked(M + TAU, e).unwrap(),
                E + TAU,
                1e-13,
            );
            assert_close(
                elliptic::ecc_anom_checked(M - 3.0 * TAU, e).unwrap(),
                E - 3.0 * TAU,
                1e-13,
            );
            // E(2π - M) = 2π - E(M)
            assert_close(
                elliptic::ecc_anom_checked(TAU - M, e).unwrap(),
                TAU - E,
                1e-13,
            );
        }
        // fixed points at perihelion and aphelion
        assert_eq!(elliptic::ecc_anom_checked(0.0, e), Some(0.0));
        assert_close(elliptic::ecc_anom_checked(PI, e).unwrap(), PI, 1e-15);
    }
}

#[test]
fn ecc_anom_is_monotonic() {
    for &e in &[0.5, 0.999] {
        let mut prev = f64::NEG_INFINITY;
        for i in -1000..=1000 {
            let E = elliptic::ecc_anom_checked(i as f64 * 0.01, e).unwrap();
            assert!(E > prev);
            prev = E;
        }
    }
}

#[test]
fn all_methods_agree() {
    for &e in &[0.05, 0.3, 0.6, 0.9, 0.98] {
        for i in -30..=30 {
            let M = i as f64 * 0.2;
            let robust = elliptic::ecc_anom_checked(M, e).unwrap();
            let bs = elliptic::ecc_anom_binary_search(M, e).unwrap();
            // near M = π the fixed-point iteration oscillates and rounding
            // keeps it from settling below ~1e-14, so ask for less
            let it = elliptic::ecc_anom_simple_iter(M, e, 1e-12, 100_000).unwrap();
            let nw =
                elliptic::ecc_anom_newton(M, e, NewtonCorrLimit::ArcsinOfSin, 1e-14, 200).unwrap();
            let cl =
                elliptic::ecc_anom_newton(M, e, NewtonCorrLimit::Clamp(0.5), 1e-14, 200).unwrap();
            for x in [bs, nw, cl] {
                assert_close(x, robust, 1e-13);
            }
            // the simple iteration's error can exceed its stopping step
            assert_close(it, robust, 1e-10);
        }
    }
}

// --- True anomaly and radius vector (eqs. 30.1 to 30.4) ----------------------

#[test]
fn true_anomaly_and_radius_vector() {
    let (a, e) = (2.5, 0.6);
    let q = elliptic::perih_dist(a, e);
    assert_close(q, 1.0, 1e-15);
    assert_close(elliptic::aph_dist(a, e), 4.0, 1e-15);
    assert_close(elliptic::semimaj_axis(q, e), a, 1e-15);
    assert_close(q + elliptic::aph_dist(a, e), 2.0 * a, 1e-15);

    for i in 1..180 {
        let M = rad(i as f64);
        let E = elliptic::ecc_anom_checked(M, e).unwrap();
        let v = elliptic::true_anom(E, e);
        // between perihelion and aphelion, 0 < M < E < v < 180°
        assert!(M < E && E < v && v < PI);
        // three expressions for the radius vector agree
        let r1 = elliptic::rad_vec_frm_ecc_anom(E, a, e);
        let r2 = elliptic::rad_vec_frm_true_anom(v, a, e);
        let r3 = astro::orbit::rad_vec_frm_true_anom(v, q, e);
        assert_close(r2, r1, 1e-13);
        assert_close(r3, r1, 1e-13);
        // after aphelion the true anomaly is negative and symmetric
        assert_close(elliptic::true_anom(-E, e), -v, 1e-15);
        assert_close(elliptic::true_anom(TAU - E, e), -v, 1e-12);
    }
    // perihelion and aphelion
    assert_eq!(elliptic::true_anom(0.0, e), 0.0);
    assert_close(elliptic::true_anom(PI, e), PI, 1e-15);
    assert_close(elliptic::rad_vec_frm_ecc_anom(0.0, a, e), q, 1e-15);
    assert_close(elliptic::rad_vec_frm_ecc_anom(PI, a, e), 4.0, 1e-15);
}

#[test]
fn ecc_anom_for_large_mean_anomalies() {
    // many revolutions away from the epoch, in both directions
    for &M in &[1e3, -1e3, 123456.789, -987654.321, 1e8] {
        for &e in &[0.1, 0.9, 0.9999] {
            let E = elliptic::ecc_anom_checked(M, e).unwrap();
            // the reduction of M to [-π, π] costs about one ulp of M
            assert!(residual(E, M, e).abs() <= 4.0 * f64::EPSILON * M.abs());
            assert!((E - M).abs() <= e + 1e-6);
            let bs = elliptic::ecc_anom_binary_search(M, e).unwrap();
            assert_close(bs, E, 8.0 * f64::EPSILON * M.abs());
        }
    }
}

/// `E - e sin(E)` evaluated without cancellation, for checking solutions
/// very close to the perihelion of a nearly parabolic orbit: written as
/// `(1 - e) sin(E) + (E - sin(E))`, the last part from its Taylor series.
fn accurate_mean_anom(E: f64, e: f64) -> f64 {
    assert!(E.abs() < 0.5);
    let x2 = E * E;
    let (mut term, mut sum) = (E * x2 / 6.0, 0.0);
    for n in (3..60).step_by(2) {
        sum += term;
        let n = n as f64;
        term *= -x2 / ((n + 1.0) * (n + 2.0));
    }
    (1.0 - e) * E.sin() + sum
}

#[test]
fn ecc_anom_full_precision_close_to_parabolic() {
    // Bug fix: f(E) = E - e sin(E) - M was evaluated directly, which loses
    // most of its digits when E is small and e close to 1 (E and e sin(E)
    // nearly cancel). For e = 1 - 1e-10 and M = 1e-14 the eccentric anomaly
    // was wrong by 1e-7 of its value although the docs promise full f64
    // precision for any e < 1.
    for &e in &[0.9999, 1.0 - 1e-6, 1.0 - 1e-8, 1.0 - 1e-10, 1.0 - 1e-12] {
        for &M in &[1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3] {
            let E = elliptic::ecc_anom_checked(M, e).unwrap();
            let rel = (accurate_mean_anom(E, e) - M) / M;
            assert!(
                rel.abs() < 1e-14,
                "e = {e}, M = {M}: relative residual {rel:e}"
            );
            let E_neg = elliptic::ecc_anom_checked(-M, e).unwrap();
            assert_eq!(E_neg, -E);
        }
    }
}

#[test]
fn ecc_anom_huge_mean_anomaly_does_not_panic() {
    // Bug fix: the mean anomaly was reduced with M - round(M / 2π) 2π,
    // which for |M| above about 1e16 left a remainder outside [-π, π]; the
    // bracket of the root was then empty and `f64::clamp` panicked (for
    // example M = 1e16 or 1e18 with e = 0.5).
    for &M in &[
        1e15,
        1e16,
        1e18,
        1e20,
        -1e16,
        -1e18,
        1e300,
        f64::MAX,
        -f64::MAX,
    ] {
        for &e in &[0.0, 0.5, 0.999] {
            let E = elliptic::ecc_anom_checked(M, e).unwrap();
            assert!(E.is_finite());
            // one ulp of M is far larger than the e sin(E) term
            assert!((E - M).abs() <= e + 2.0 * f64::EPSILON * M.abs(), "M = {M}");
            let bs = elliptic::ecc_anom_binary_search(M, e).unwrap();
            assert!(
                (bs - M).abs() <= e + 2.0 * f64::EPSILON * M.abs(),
                "M = {M}"
            );
        }
    }
    // the reduction is exact: 2π-shifts of a moderately large M still give
    // the same eccentric anomaly to within the precision of M itself
    let (M, e) = (1e12 + 0.3, 0.7);
    let E = elliptic::ecc_anom_checked(M, e).unwrap();
    assert!(residual(E, M, e).abs() <= 2.0 * f64::EPSILON * M);
}

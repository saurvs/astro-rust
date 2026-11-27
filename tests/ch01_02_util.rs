//! Meeus, chapters 1 (hints and tips) and 2 (about accuracy): polynomial
//! evaluation and rounding helpers in `astro::util`.

mod common;

use astro::Horner_eval;
use astro::util;
use common::assert_close;

// --- Horner's method (chapter 1) -------------------------------------------

#[test]
fn horner_matches_explicit_polynomial() {
    // y = A + Bx + Cx² + Dx³ + Ex⁴ written out term by term
    let c = [1.5, -2.25, 0.75, 3.0, -0.5];
    for &x in &[-3.0_f64, -1.0, -0.3, 0.0, 0.7, 2.0, 10.0] {
        let explicit = c[0] + c[1] * x + c[2] * x * x + c[3] * x.powi(3) + c[4] * x.powi(4);
        assert_close(
            util::horner(x, &c),
            explicit,
            1e-12 * explicit.abs().max(1.0),
        );
        assert_close(
            Horner_eval!(x, 1.5, -2.25, 0.75, 3.0, -0.5),
            explicit,
            1e-12 * explicit.abs().max(1.0),
        );
    }
}

#[test]
fn horner_edge_cases() {
    assert_eq!(util::horner(3.0, &[]), 0.0);
    assert_eq!(util::horner(3.0, &[7.0]), 7.0);
    assert_eq!(util::horner(0.0, &[7.0, 1.0, 2.0]), 7.0);
    assert_eq!(Horner_eval!(3.0, 7.0), 7.0);
    assert_eq!(Horner_eval!(3.0, 7.0,), 7.0);
    assert_eq!(Horner_eval!(2.0, 1.0, 2.0, 3.0,), 17.0);
    // Works for f32 too
    let y: f32 = Horner_eval!(2.0_f32, 1.0, 2.0, 3.0);
    assert_eq!(y, 17.0);
}

#[test]
fn horner_eval_is_true_horner_scheme() {
    // The old macro accumulated explicit powers of x, so x³ overflowed to
    // infinity and inf * 0 gave NaN even though the polynomial is just x.
    let big = 1e200;
    assert_eq!(Horner_eval!(big, 0.0, 1.0, 0.0), big);
    assert_eq!(util::horner(big, &[0.0, 1.0, 0.0]), big);
}

#[test]
fn horner_eval_evaluates_its_argument_once() {
    let mut calls = 0;
    let y = Horner_eval!(
        {
            calls += 1;
            2.0
        },
        1.0,
        1.0,
        1.0,
        1.0
    );
    assert_eq!(y, 15.0);
    assert_eq!(calls, 1);
}

#[test]
fn horner_uranus_eccentricity_out_of_range() {
    // Meeus, eq. 1.1: e = 0.04638122 − 0.000027293 T + 0.0000000789 T².
    // Far outside its range of validity, at T = −3307.9, it gives e = 1.
    let e = util::horner(-3307.9, &[0.04638122, -0.000027293, 0.0000000789]);
    assert_close(e, 1.0, 5e-5);
    assert_close(
        util::horner(0.0, &[0.04638122, -0.000027293, 0.0000000789]),
        0.04638122,
        0.0,
    );
}

// --- Rounding (chapter 2) ----------------------------------------------------

#[test]
fn round_to_nearest() {
    // Meeus, chapter 2, "Rounding the final result"
    assert_eq!(util::round_upto_digits(15.88, 1), 15.9);
    assert_eq!(util::round_upto_digits(15.88, 0), 16.0);
    // Round the final result, not the inputs: 1.4 + 1.4 rounds to 3, not 2
    assert_eq!(util::round_upto_digits(1.4 + 1.4, 0), 3.0);
    assert_eq!(
        util::round_upto_digits(1.4, 0) + util::round_upto_digits(1.4, 0),
        2.0
    );
}

#[test]
fn round_jupiter_magnitude() {
    // Meeus, chapter 2: m = −8.93 + 5 log(r Δ) for Jupiter on 1992 May 14
    let m = -8.93 + 5.0 * (5.417149_f64 * 5.125382).log10();
    assert_close(m, -1.712514898, 5e-10);
    assert_eq!(util::round_upto_digits(m, 2), -1.71);
}

#[test]
fn round_neptune_altitude() {
    // Meeus, chapter 2: altitude at transit h = 90° − φ + δ for Neptune
    // (δ = −20°24′) seen from φ = +50°23′. Rounding at the end gives 19°;
    // rounding φ and δ first would give the wrong 20°.
    let phi = 50.0 + 23.0 / 60.0;
    let delta = -(20.0 + 24.0 / 60.0);
    let h = 90.0 - phi + delta;
    assert_close(h, 19.0 + 13.0 / 60.0, 1e-12);
    assert_eq!(util::round_upto_digits(h, 0), 19.0);
    let early = 90.0 - util::round_upto_digits(phi, 0) + util::round_upto_digits(delta, 0);
    assert_eq!(early, 20.0);
}

#[test]
fn round_negative_and_ties() {
    assert_eq!(util::round_upto_digits(-2.345678, 3), -2.346);
    assert_eq!(util::round_upto_digits(-2.345278, 3), -2.345);
    assert_eq!(util::round_upto_digits(2.5, 0), 3.0);
    assert_eq!(util::round_upto_digits(-2.5, 0), -3.0);
    assert_eq!(util::round_upto_digits(0.0, 4), 0.0);
    assert_eq!(util::round_upto_digits(123.0, 0), 123.0);
}

#[test]
fn round_non_finite_and_huge() {
    assert!(util::round_upto_digits(f64::NAN, 2).is_nan());
    assert_eq!(util::round_upto_digits(f64::INFINITY, 2), f64::INFINITY);
    assert_eq!(
        util::round_upto_digits(f64::NEG_INFINITY, 2),
        f64::NEG_INFINITY
    );
    // Already integral at the requested precision: returned unchanged
    assert_eq!(util::round_upto_digits(1e300, 5), 1e300);
    assert_eq!(util::round_upto_digits(f64::MAX, 3), f64::MAX);
}

#[test]
fn round_with_absurd_digit_count_does_not_panic() {
    // Used to overflow the loop bound `decimal_digits + 1` (a panic in debug
    // builds) and to loop four billion times.
    assert_eq!(util::round_upto_digits(1.5, u32::MAX), 1.5);
    assert_eq!(util::round_upto_digits(0.123456789, 400), 0.123456789);
    assert_eq!(util::round_upto_digits(0.123456789, 309), 0.123456789);
}

#[test]
fn round_matches_previous_behaviour_for_ordinary_digits() {
    // The power of ten is still built by exact repeated multiplication, so
    // ordinary roundings are bit-for-bit unchanged.
    for digits in 0..=15 {
        let d = 10f64.powi(digits as i32);
        for &x in &[
            0.876125,
            -1.712514898,
            13.369480613,
            54.0 * 60.0 + 13.369,
            1e-7,
        ] {
            assert_eq!(util::round_upto_digits(x, digits), (x * d).round() / d);
        }
    }
}

// --- Accuracy of the arithmetic (chapter 2) ---------------------------------

#[test]
fn squaring_test_shows_binary_representation_error() {
    // Meeus, chapter 2: squaring 1.0000001 27 times should give 674530.4707
    // to ten significant figures. Binary double precision gives 674530.4755,
    // the value Meeus lists for the HP-Integral and QuickBASIC 4.5: 1.0000001
    // has no exact binary representation, and squaring 27 times multiplies
    // that tiny input error by 2^27.
    let mut x = 1.0000001_f64;
    for _ in 0..27 {
        x *= x;
    }
    assert_close(x, 674530.4755, 5e-5);
    assert!((x - 674530.4707).abs() > 1e-3);
}

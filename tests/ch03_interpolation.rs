//! Meeus, chapter 3 (interpolation): `astro::interpol`.

mod common;

use astro::interpol::{self, Extremum};
use astro::util;
use common::{assert_close, dms, hms};
use std::f64::consts::TAU;

/// Arcseconds of an angle given in degrees, arcminutes and arcseconds.
fn arcsec(neg: bool, d: u32, m: u32, s: f64) -> f64 {
    dms(neg, d, m, s) * 3600.0
}

// --- Three tabular values ---------------------------------------------------------

#[test]
fn example_3a_three_values() {
    // Meeus, example 3.a: distance Mars–Earth (AU) on 1992 November 7, 8, 9
    // at 0h TD, interpolated to November 8 at 4h21m TD
    let n = 4.35 / 24.0;
    assert_close(n, 0.18125, 1e-15);
    let y = interpol::three_values(0.884226, 0.877366, 0.870531, n);
    assert_close(y, 0.876125, 5e-7);
}

#[test]
fn three_values_reproduces_table_and_exact_parabola() {
    let (y1, y2, y3) = (0.884226, 0.877366, 0.870531);
    assert_close(interpol::three_values(y1, y2, y3, -1.0), y1, 1e-15);
    assert_close(interpol::three_values(y1, y2, y3, 0.0), y2, 0.0);
    assert_close(interpol::three_values(y1, y2, y3, 1.0), y3, 1e-15);

    // Exact for any parabola: y = 3 + 2x − 3x² (Meeus, example 3.d)
    let p = |x: f64| 3.0 + 2.0 * x - 3.0 * x * x;
    for &n in &[-1.0, -0.7, -0.25, 0.0, 0.3, 0.5, 1.0] {
        assert_close(
            interpol::three_values(p(-1.0), p(0.0), p(1.0), n),
            p(n),
            1e-14,
        );
    }
}

#[test]
fn three_values_reversed_table() {
    // Reading the table backwards is the same as negating n
    let (y1, y2, y3) = (0.884226, 0.877366, 0.870531);
    for &n in &[-0.8, -0.3, 0.1, 0.45] {
        assert_close(
            interpol::three_values(y3, y2, y1, -n),
            interpol::three_values(y1, y2, y3, n),
            1e-15,
        );
    }
}

#[test]
fn example_3b_extremum_three_values() {
    // Meeus, example 3.b: distance Sun–Mars on 1992 May 12, 16 and 20 at 0h TD
    let ext = interpol::extremum_three_values(1.3814294, 1.3812213, 1.3812453).unwrap();
    assert_close(ext.y, 1.3812030, 5e-8);
    assert_close(ext.n, 0.39660, 5e-6);
    // Perihelion 4 × n days after May 16.0: May 17 at 14h TD
    let days = 4.0 * ext.n;
    assert_close(days, 1.58640, 2e-5);
    assert_eq!(util::round_upto_digits((days - 1.0) * 24.0, 0), 14.0);
    // The extremum is the vertex of the parabola: same value from eq. 3.3
    assert_close(
        interpol::three_values(1.3814294, 1.3812213, 1.3812453, ext.n),
        ext.y,
        1e-15,
    );
}

#[test]
fn extremum_three_values_invariants() {
    // Vertex of y = 3 + 2x − 3x² is at x = 1/3, y = 10/3
    let ext = interpol::extremum_three_values(-2.0, 3.0, 2.0).unwrap();
    assert_close(ext.n, 1.0 / 3.0, 1e-15);
    assert_close(ext.y, 10.0 / 3.0, 1e-15);
    // Symmetric table: extremum at the central value
    assert_eq!(
        interpol::extremum_three_values(1.0, 3.0, 1.0),
        Some(Extremum { n: 0.0, y: 3.0 })
    );
    // Reversing the table mirrors the extremum
    let rev = interpol::extremum_three_values(2.0, 3.0, -2.0).unwrap();
    assert_close(rev.n, -ext.n, 1e-15);
    assert_close(rev.y, ext.y, 1e-15);
    // Straight line or constant: no extremum
    assert_eq!(interpol::extremum_three_values(1.0, 2.0, 3.0), None);
    assert_eq!(interpol::extremum_three_values(5.0, 5.0, 5.0), None);
}

#[test]
fn extremum_far_outside_table_is_absurd() {
    // Meeus, chapter 3, important remark 3: the three central values of
    // sin 30°, 31°, 32° give an absurd maximum 1.22827 at 31° + 95.35°
    let s = [0.5000000000, 0.5150380749, 0.5299192642];
    let ext = interpol::extremum_three_values(s[0], s[1], s[2]).unwrap();
    assert_close(ext.y, 1.22827, 5e-6);
    assert_close(ext.n, 95.35, 5e-3);
}

#[test]
fn example_3c_zero_three_values() {
    // Meeus, example 3.c: declination of Mercury on 1973 February 26, 27
    // and 28 at 0h TD, in arcseconds
    let y1 = arcsec(true, 0, 28, 13.4);
    let y2 = arcsec(false, 0, 6, 46.3);
    let y3 = arcsec(false, 0, 38, 23.2);
    assert_close(y1, -1693.4, 1e-9);
    assert_close(y2, 406.3, 1e-9);
    assert_close(y3, 2303.2, 1e-9);

    let n0 = interpol::zero_three_values(y1, y2, y3).unwrap();
    assert_close(n0, -0.20127, 5e-6);
    // 1973 February 26.79873 = February 26 at 19h10m TD
    assert_close(27.0 + n0, 26.79873, 5e-6);
    assert_eq!(
        util::round_upto_digits((27.0 + n0 - 26.0) * 24.0 * 60.0, 0),
        19.0 * 60.0 + 10.0
    );
    assert_close(interpol::three_values(y1, y2, y3, n0), 0.0, 1e-9);
}

#[test]
fn example_3d_zero_three_values_strong_curvature() {
    // Meeus, example 3.d: y = 3 + 2x − 3x² at x = −1, 0, +1; eq. 3.7 gives
    // −0.720759220056 to 12 decimals in six steps
    let n0 = interpol::zero_three_values(-2.0, 3.0, 2.0).unwrap();
    assert_close(n0, -0.720759220056, 5e-13);
    assert_close(n0, (1.0 - 10f64.sqrt()) / 3.0, 1e-15);
}

#[test]
fn zero_three_values_edge_cases() {
    // No real zero
    assert_eq!(interpol::zero_three_values(5.0, 2.0, 1.0), None);
    assert_eq!(interpol::zero_three_values(1.0, 2.0, 4.0), None);
    // Constant non-zero and constant zero
    assert_eq!(interpol::zero_three_values(1.0, 1.0, 1.0), None);
    assert_eq!(interpol::zero_three_values(0.0, 0.0, 0.0), Some(0.0));
    // Zero exactly at the central value
    assert_eq!(interpol::zero_three_values(-1.0, 0.0, 4.0), Some(0.0));
    // Straight line y = 2 + 3n
    assert_close(
        interpol::zero_three_values(-1.0, 2.0, 5.0).unwrap(),
        -2.0 / 3.0,
        1e-15,
    );
    // Double root (parabola touching zero) at n = 0.5: y = (n − 0.5)²
    let p = |n: f64| (n - 0.5) * (n - 0.5);
    assert_close(
        interpol::zero_three_values(p(-1.0), p(0.0), p(1.0)).unwrap(),
        0.5,
        1e-15,
    );
    // Huge values do not overflow
    let n = interpol::zero_three_values(-1e300, 1e300, 3e300).unwrap();
    assert_close(n, -0.5, 1e-15);
    // Non-finite input
    assert_eq!(interpol::zero_three_values(f64::NAN, 1.0, 2.0), None);
    assert_eq!(interpol::zero_three_values(f64::INFINITY, 1.0, 2.0), None);
}

#[test]
fn zero_three_values_picks_root_nearest_centre() {
    // y = (n − 0.2)(n − 5): zeros at 0.2 and 5
    let p = |n: f64| (n - 0.2) * (n - 5.0);
    let n0 = interpol::zero_three_values(p(-1.0), p(0.0), p(1.0)).unwrap();
    assert_close(n0, 0.2, 1e-15);
    // y = (n + 0.3)(n + 4)
    let p = |n: f64| (n + 0.3) * (n + 4.0);
    let n0 = interpol::zero_three_values(p(-1.0), p(0.0), p(1.0)).unwrap();
    assert_close(n0, -0.3, 1e-15);
}

#[test]
fn zero_three_values_symmetries() {
    let (y1, y2, y3) = (-1693.4, 406.3, 2303.2);
    let n0 = interpol::zero_three_values(y1, y2, y3).unwrap();
    // Negating or scaling the values does not move the zero
    assert_close(
        interpol::zero_three_values(-y1, -y2, -y3).unwrap(),
        n0,
        1e-15,
    );
    assert_close(
        interpol::zero_three_values(y1 / 3600.0, y2 / 3600.0, y3 / 3600.0).unwrap(),
        n0,
        1e-15,
    );
    // Reversing the table mirrors it
    assert_close(interpol::zero_three_values(y3, y2, y1).unwrap(), -n0, 1e-15);
}

// --- Five tabular values ----------------------------------------------------------

#[test]
fn example_3e_five_values() {
    // Meeus, example 3.e: equatorial horizontal parallax of the Moon on
    // 1992 February 27.0 to 29.0 every 12 hours, in arcseconds; value on
    // February 28 at 3h20m TD
    let y = [
        arcsec(false, 0, 54, 36.125),
        arcsec(false, 0, 54, 24.606),
        arcsec(false, 0, 54, 15.486),
        arcsec(false, 0, 54, 8.694),
        arcsec(false, 0, 54, 4.133),
    ];
    let n = (3.0 + 20.0 / 60.0) / 12.0;
    assert_close(n, 0.2777778, 5e-8);
    let p = interpol::five_values(y[0], y[1], y[2], y[3], y[4], n);
    assert_close(p, arcsec(false, 0, 54, 13.369), 5e-4);
    // Neglecting the third differences (three values only) changes the
    // result in the third decimal of the arcsecond
    let p3 = interpol::three_values(y[1], y[2], y[3], n);
    assert!((p3 - p).abs() > 0.002 && (p3 - p).abs() < 0.1, "{}", p3 - p);
}

#[test]
fn five_values_reproduces_table_and_exact_quartic() {
    let q = |x: f64| 1.0 - 0.5 * x + 0.25 * x * x - 0.125 * x.powi(3) + 0.0625 * x.powi(4);
    let y = [q(-2.0), q(-1.0), q(0.0), q(1.0), q(2.0)];
    for (i, &n) in [-2.0, -1.0, 0.0, 1.0, 2.0].iter().enumerate() {
        assert_close(
            interpol::five_values(y[0], y[1], y[2], y[3], y[4], n),
            y[i],
            1e-14,
        );
    }
    for &n in &[-1.0, -0.6, -0.1, 0.25, 0.5, 0.9] {
        assert_close(
            interpol::five_values(y[0], y[1], y[2], y[3], y[4], n),
            q(n),
            1e-14,
        );
    }
    // For a parabola, five and three values agree
    let p = |x: f64| 3.0 + 2.0 * x - 3.0 * x * x;
    for &n in &[-0.9, -0.4, 0.2, 0.7] {
        assert_close(
            interpol::five_values(p(-2.0), p(-1.0), p(0.0), p(1.0), p(2.0), n),
            interpol::three_values(p(-1.0), p(0.0), p(1.0), n),
            1e-14,
        );
    }
}

#[test]
fn extremum_five_values_exact_quartic() {
    // y = 5 − (n − 0.3)² + 0.05 (n − 0.3)⁴ has a maximum 5 at n = 0.3
    let q = |n: f64| 5.0 - (n - 0.3).powi(2) + 0.05 * (n - 0.3).powi(4);
    let ext = interpol::extremum_five_values(q(-2.0), q(-1.0), q(0.0), q(1.0), q(2.0)).unwrap();
    assert_close(ext.n, 0.3, 1e-13);
    assert_close(ext.y, 5.0, 1e-13);
}

#[test]
fn extremum_five_values_matches_three_values_for_parabola() {
    // Sun–Mars distances of example 3.b continued as an exact parabola
    let (y1, y2, y3) = (1.3814294, 1.3812213, 1.3812453);
    let p = |n: f64| interpol::three_values(y1, y2, y3, n);
    let e3 = interpol::extremum_three_values(y1, y2, y3).unwrap();
    let e5 = interpol::extremum_five_values(p(-2.0), p(-1.0), p(0.0), p(1.0), p(2.0)).unwrap();
    assert_close(e5.n, e3.n, 1e-9);
    assert_close(e5.y, e3.y, 1e-15);
}

#[test]
fn extremum_five_values_sine_example() {
    // Meeus, chapter 3, important remark 3: sin 29° ... sin 33° (ten
    // decimals) and eq. 3.9. Meeus prints n = +57.30 (maximum at 88.30°)
    // and a maximum of 0.99348. The root of eq. 3.9 for these values is in
    // fact n = 57.267, where the interpolating polynomial reaches 0.99330;
    // 0.99348 exceeds the largest value that polynomial takes anywhere. The
    // simple iteration of eq. 3.9 from n = 0 oscillates between about 8.9
    // and 94.2 for these data, so the book's figures are not reproducible.
    // The point of the remark stands: the result is far from 90° and 1.
    let s = [
        0.4848096202,
        0.5000000000,
        0.5150380749,
        0.5299192642,
        0.5446390350,
    ];
    let ext = interpol::extremum_five_values(s[0], s[1], s[2], s[3], s[4]).unwrap();
    assert_close(ext.n, 57.2669, 5e-4);
    assert_close(ext.y, 0.99330, 5e-6);
    // It is a stationary point of the interpolating polynomial
    let f = |n: f64| interpol::five_values(s[0], s[1], s[2], s[3], s[4], n);
    assert!(f(ext.n) >= f(ext.n - 0.01) && f(ext.n) >= f(ext.n + 0.01));
}

#[test]
fn extremum_five_values_none() {
    // Straight line: no extremum
    assert_eq!(
        interpol::extremum_five_values(1.0, 2.0, 3.0, 4.0, 5.0),
        None
    );
    // Constant
    assert_eq!(
        interpol::extremum_five_values(2.0, 2.0, 2.0, 2.0, 2.0),
        None
    );
    // y = n³ has a stationary inflection at 0, not an extremum
    assert_eq!(
        interpol::extremum_five_values(-8.0, -1.0, 0.0, 1.0, 8.0),
        None
    );
}

#[test]
fn exercise_zero_five_values() {
    // Meeus, chapter 3, exercise: heliocentric latitude of Mercury on 1988
    // January 25.0 to 29.0 TD, in arcseconds; it was zero for n = −0.361413,
    // on January 26.638587 = January 26 at 15h20m TD
    let y = [
        arcsec(true, 1, 11, 21.23),
        arcsec(true, 0, 28, 12.31),
        arcsec(false, 0, 16, 7.02),
        arcsec(false, 1, 1, 0.13),
        arcsec(false, 1, 45, 46.33),
    ];
    let n0 = interpol::zero_five_values(y[0], y[1], y[2], y[3], y[4]).unwrap();
    assert_close(n0, -0.361413, 5e-7);
    assert_close(27.0 + n0, 26.638587, 5e-7);
    assert_eq!(
        util::round_upto_digits((27.0 + n0 - 26.0) * 24.0 * 60.0, 0),
        15.0 * 60.0 + 20.0
    );
    assert_close(
        interpol::five_values(y[0], y[1], y[2], y[3], y[4], n0),
        0.0,
        1e-9,
    );

    // From the three central values (eq. 3.6) instead: −0.362166, about
    // 1.1 minutes later
    let n3 = interpol::zero_three_values(y[1], y[2], y[3]).unwrap();
    assert_close(n3, -0.362166, 5e-7);
    assert_close((n0 - n3) * 1440.0, 1.1, 0.05);
}

#[test]
fn zero_five_values_edge_cases() {
    // Exact quartic with a zero at n = 0.4: y = (n − 0.4)(n + 3)(n² + 1)
    let q = |n: f64| (n - 0.4) * (n + 3.0) * (n * n + 1.0);
    let n0 = interpol::zero_five_values(q(-2.0), q(-1.0), q(0.0), q(1.0), q(2.0)).unwrap();
    assert_close(n0, 0.4, 1e-14);
    // No zero: y = 1 + n²
    assert_eq!(interpol::zero_five_values(5.0, 2.0, 1.0, 2.0, 5.0), None);
    // Constant zero and constant non-zero
    assert_eq!(
        interpol::zero_five_values(0.0, 0.0, 0.0, 0.0, 0.0),
        Some(0.0)
    );
    assert_eq!(interpol::zero_five_values(3.0, 3.0, 3.0, 3.0, 3.0), None);
    // Zero exactly at the centre
    assert_eq!(
        interpol::zero_five_values(-2.0, -1.0, 0.0, 1.0, 2.0),
        Some(0.0)
    );
    // NaN input
    assert_eq!(
        interpol::zero_five_values(f64::NAN, 1.0, 2.0, 3.0, 4.0),
        None
    );
}

#[test]
fn zero_five_values_agrees_with_three_values_for_parabola() {
    let p = |n: f64| 3.0 + 2.0 * n - 3.0 * n * n;
    let n5 = interpol::zero_five_values(p(-2.0), p(-1.0), p(0.0), p(1.0), p(2.0)).unwrap();
    let n3 = interpol::zero_three_values(p(-1.0), p(0.0), p(1.0)).unwrap();
    assert_close(n5, n3, 1e-14);
}

// --- Interpolation to halves --------------------------------------------------------

#[test]
fn example_3f_halfway() {
    // Meeus, example 3.f: apparent right ascension of the Moon on 1994 March
    // 25 at 8h, 10h, 12h, 14h TD, seconds after 10h; value at 11h TD
    let y = [1128.732, 1402.835, 1677.247, 1951.983];
    let ra = interpol::four_values_halfway(y[0], y[1], y[2], y[3]);
    assert_close(ra, 1540.001, 5e-4);
    // α = 10h25m40.001s
    assert_close(
        10.0 * 3600.0 + ra,
        10.0 * 3600.0 + 25.0 * 60.0 + 40.001,
        5e-4,
    );
    assert_close(
        hms(10, 25, 40.001),
        (10.0 * 3600.0 + ra) / 240.0,
        5e-4 / 240.0,
    );
}

#[test]
fn halfway_exact_for_cubic_and_agrees_with_lagrange() {
    let c = |x: f64| 2.0 - x + 0.5 * x * x + 0.1 * x.powi(3);
    let xs = [-1.0, 0.0, 1.0, 2.0];
    let ys = xs.map(c);
    let mid = interpol::four_values_halfway(ys[0], ys[1], ys[2], ys[3]);
    assert_close(mid, c(0.5), 1e-14);
    assert_close(mid, interpol::lagrange(&xs, &ys, 0.5).unwrap(), 1e-14);
}

// --- Lagrange interpolation ----------------------------------------------------------

#[test]
fn lagrange_sine_exercise() {
    // Meeus, chapter 3: six unequally spaced, unordered values of sin x
    let xs = [29.43, 30.97, 27.69, 28.11, 31.58, 33.05];
    let ys = [
        0.4913598528,
        0.5145891926,
        0.4646875083,
        0.4711658342,
        0.5236885653,
        0.5453707057,
    ];
    assert_close(interpol::lagrange(&xs, &ys, 30.0).unwrap(), 0.5, 5e-11);
    // Far outside the data Meeus quotes +0.0000482 at 0° and +1.00007 at
    // 90°. Evaluating the formula in exact rational arithmetic on the
    // printed data gives +0.0000512248 and +0.9999648078 instead: the
    // Lagrange coefficients reach about 10⁶ there, so the last-digit
    // rounding of the inputs dominates and the book's figures are not
    // reproducible. In double precision the binary rounding of the decimal
    // abscissae is amplified the same way, by up to a few units in 10⁻⁹.
    assert_close(
        interpol::lagrange(&xs, &ys, 0.0).unwrap(),
        0.0000512248,
        5e-9,
    );
    assert_close(
        interpol::lagrange(&xs, &ys, 90.0).unwrap(),
        0.9999648078,
        5e-9,
    );
    // The order of the points does not matter
    let mut pts: Vec<(f64, f64)> = xs.iter().copied().zip(ys.iter().copied()).collect();
    pts.sort_by(|a, b| a.0.total_cmp(&b.0));
    let (sx, sy): (Vec<f64>, Vec<f64>) = pts.into_iter().unzip();
    for &x in &[28.0, 30.0, 32.5] {
        assert_close(
            interpol::lagrange(&sx, &sy, x).unwrap(),
            interpol::lagrange(&xs, &ys, x).unwrap(),
            1e-15,
        );
    }
}

#[test]
fn example_3g_lagrange_polynomial() {
    // Meeus, example 3.g: the cubic through (1, −6), (3, 6), (4, 9), (6, 15)
    // is (x³ − 13x² + 69x − 87)/5
    let xs = [1.0, 3.0, 4.0, 6.0];
    let ys = [-6.0, 6.0, 9.0, 15.0];
    let c = interpol::lagrange_polynomial(&xs, &ys).unwrap();
    let expected = [-87.0 / 5.0, 69.0 / 5.0, -13.0 / 5.0, 1.0 / 5.0];
    assert_eq!(c.len(), 4);
    for (got, want) in c.iter().zip(expected) {
        assert_close(*got, want, 1e-13);
    }
    // The polynomial passes through the points and agrees with lagrange()
    for (&x, &y) in xs.iter().zip(&ys) {
        assert_close(util::horner(x, &c), y, 1e-12);
    }
    for &x in &[-2.0, 0.0, 2.5, 5.0, 10.0] {
        assert_close(
            util::horner(x, &c),
            interpol::lagrange(&xs, &ys, x).unwrap(),
            1e-11,
        );
    }
}

#[test]
fn lagrange_small_and_invalid_inputs() {
    // One point: a constant; two points: a straight line
    assert_eq!(interpol::lagrange(&[2.0], &[7.0], 100.0), Some(7.0));
    assert_eq!(
        interpol::lagrange_polynomial(&[2.0], &[7.0]),
        Some(vec![7.0])
    );
    assert_close(
        interpol::lagrange(&[0.0, 2.0], &[1.0, 5.0], 1.5).unwrap(),
        4.0,
        1e-15,
    );
    let line = interpol::lagrange_polynomial(&[0.0, 2.0], &[1.0, 5.0]).unwrap();
    assert_close(line[0], 1.0, 1e-15);
    assert_close(line[1], 2.0, 1e-15);
    // Invalid: empty, mismatched lengths, repeated or non-finite abscissae
    assert_eq!(interpol::lagrange(&[], &[], 1.0), None);
    assert_eq!(interpol::lagrange(&[1.0, 2.0], &[1.0], 1.0), None);
    assert_eq!(
        interpol::lagrange(&[1.0, 2.0, 1.0], &[1.0, 2.0, 3.0], 1.5),
        None
    );
    assert_eq!(interpol::lagrange(&[1.0, f64::NAN], &[1.0, 2.0], 1.5), None);
    assert_eq!(interpol::lagrange_polynomial(&[], &[]), None);
    assert_eq!(
        interpol::lagrange_polynomial(&[3.0, 3.0], &[1.0, 2.0]),
        None
    );
    assert_eq!(interpol::lagrange_polynomial(&[1.0], &[1.0, 2.0]), None);
}

#[test]
fn lagrange_matches_three_and_five_values_for_equal_spacing() {
    let y = [36.125, 24.606, 15.486, 8.694, 4.133];
    let xs5 = [-2.0, -1.0, 0.0, 1.0, 2.0];
    for &n in &[-0.5, 0.2777778, 0.9] {
        assert_close(
            interpol::lagrange(&xs5, &y, n).unwrap(),
            interpol::five_values(y[0], y[1], y[2], y[3], y[4], n),
            1e-12,
        );
        assert_close(
            interpol::lagrange(&xs5[1..4], &y[1..4], n).unwrap(),
            interpol::three_values(y[1], y[2], y[3], n),
            1e-12,
        );
    }
}

// --- Angles and times that wrap around -------------------------------------------------

#[test]
fn unwrap_central_meridian_longitudes() {
    // Meeus, chapter 3, important remark 2: longitude of the Sun's central
    // meridian on 1992 June 14–17 at 0h UT; 358.25° must be taken as −1.75°
    let mut l = [37.96, 24.72, 11.48, 358.25];
    interpol::unwrap_angles(&mut l, 360.0);
    assert_eq!(&l[..3], &[37.96, 24.72, 11.48]);
    assert_close(l[3], -1.75, 1e-12);
    // The daily change is then about −13.24° throughout
    for w in l.windows(2) {
        assert_close(w[1] - w[0], -13.24, 0.02);
    }
}

#[test]
fn unwrap_right_ascension_hours() {
    // Meeus, chapter 3, important remark 2: right ascension of Mercury on
    // 1992 April 5.0, 6.0 and 7.0 TD; 0h01m00.71s must be read as 24h01m.
    let h = |hh: f64, mm: f64, ss: f64| hh + mm / 60.0 + ss / 3600.0;
    let mut ra = [
        h(23.0, 51.0, 56.04),
        h(23.0, 56.0, 28.49),
        h(0.0, 1.0, 0.71),
    ];
    let naive = interpol::three_values(ra[0], ra[1], ra[2], 0.2743);
    interpol::unwrap_angles(&mut ra, 24.0);
    assert_close(ra[2], h(24.0, 1.0, 0.71), 1e-12);
    // At April 6.2743 the right ascension lies between the 6th and 7th
    let ra_t = interpol::three_values(ra[0], ra[1], ra[2], 0.2743).rem_euclid(24.0);
    assert!(ra_t > ra[1] && ra_t < 24.0, "{ra_t}");
    // Without unwrapping the result is hours off
    assert!((naive - ra_t).abs() > 1.0);
}

#[test]
fn unwrap_angles_radians_and_edge_cases() {
    // Increasing through 2π
    let mut a = [6.0, 6.2, 0.1, 0.3];
    interpol::unwrap_angles(&mut a, TAU);
    assert_close(a[2], 0.1 + TAU, 1e-15);
    assert_close(a[3], 0.3 + TAU, 1e-15);
    // Several whole turns between entries are removed
    let mut b = [10.0, 10.0 + 3.0 * 360.0 + 5.0, 20.0 - 720.0];
    interpol::unwrap_angles(&mut b, 360.0);
    assert_close(b[1], 15.0, 1e-12);
    assert_close(b[2], 20.0, 1e-12);
    // Already continuous data are unchanged
    let mut c = [1.0, 2.0, 3.0];
    interpol::unwrap_angles(&mut c, 360.0);
    assert_eq!(c, [1.0, 2.0, 3.0]);
    // Empty and single-value slices, and invalid periods, are no-ops
    let mut empty: [f64; 0] = [];
    interpol::unwrap_angles(&mut empty, 360.0);
    let mut one = [400.0];
    interpol::unwrap_angles(&mut one, 360.0);
    assert_eq!(one, [400.0]);
    let mut d = [10.0, 350.0];
    for period in [0.0, -360.0, f64::NAN, f64::INFINITY] {
        interpol::unwrap_angles(&mut d, period);
        assert_eq!(d, [10.0, 350.0]);
    }
}

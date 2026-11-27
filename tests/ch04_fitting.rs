//! Meeus, chapter 4 (curve fitting): `astro::fitting`.

mod common;

use astro::fitting::{self, LinearFit, QuadraticFit};
use common::assert_close;

/// Meeus, table 4.A: comet Wild 2 (1978); x = log r, y = m − 5 log Δ.
const WILD2_X: [f64; 19] = [
    0.2982, 0.2969, 0.2918, 0.2905, 0.2707, 0.2574, 0.2485, 0.2287, 0.2238, 0.2156, 0.1992, 0.1948,
    0.1931, 0.1889, 0.1781, 0.1772, 0.1770, 0.1755, 0.1746,
];
const WILD2_Y: [f64; 19] = [
    10.92, 11.01, 10.99, 10.78, 10.87, 10.80, 10.75, 10.14, 10.21, 9.97, 9.69, 9.57, 9.66, 9.63,
    9.65, 9.44, 9.44, 9.32, 9.20,
];

/// Meeus, table 4.B: sunspot maxima 1761–1989; x = rise time (months),
/// y = height of the maximum.
const RISE_X: [f64; 22] = [
    73., 38., 35., 42., 78., 68., 74., 42., 52., 54., 39., 61., 42., 49., 50., 62., 44., 39., 43.,
    54., 44., 37.,
];
const RISE_Y: [f64; 22] = [
    90.4, 125.3, 161.8, 143.4, 52.5, 50.8, 71.5, 152.8, 131.3, 98.5, 144.8, 78.1, 89.5, 63.9,
    112.1, 82.0, 119.8, 161.2, 208.4, 111.6, 167.1, 162.1,
];

/// Meeus, table 4.C: yearly mean sunspot number (x) and annual rainfall at
/// Uccle in mm (y), 1901–1996.
const SUNSPOTS: [f64; 96] = [
    2.7, 5.0, 24.4, 42.0, 63.5, 53.8, 62.0, 48.5, 43.9, 18.6, 5.7, 3.6, 1.4, 9.6, 47.4, 57.1,
    103.9, 80.6, 63.6, 37.6, 26.1, 14.2, 5.8, 16.7, 44.3, 63.9, 69.0, 77.8, 64.9, 35.7, 21.2, 11.1,
    5.7, 8.7, 36.1, 79.7, 114.4, 109.6, 88.8, 67.8, 47.5, 30.6, 16.3, 9.6, 33.2, 92.6, 151.6,
    136.3, 134.7, 83.9, 69.4, 31.5, 13.9, 4.4, 38.0, 141.7, 190.2, 184.8, 159.0, 112.3, 53.9, 37.5,
    27.9, 10.2, 15.1, 47.0, 93.8, 105.9, 105.5, 104.5, 66.6, 68.9, 38.0, 34.5, 15.5, 12.6, 27.5,
    92.5, 155.4, 154.6, 140.5, 115.9, 66.6, 45.9, 17.9, 13.4, 29.2, 100.2, 157.6, 142.6, 145.7,
    94.3, 54.6, 29.9, 17.5, 8.6,
];
const RAINFALL: [f64; 96] = [
    700., 762., 854., 663., 912., 821., 622., 678., 842., 990., 741., 941., 801., 877., 910.,
    1054., 851., 848., 980., 760., 417., 938., 917., 849., 1075., 896., 837., 882., 688., 953.,
    858., 858., 738., 707., 916., 763., 900., 711., 928., 837., 744., 841., 738., 766., 745., 861.,
    640., 792., 521., 951., 878., 926., 557., 741., 616., 795., 801., 834., 560., 962., 903., 862.,
    713., 785., 1073., 1054., 707., 776., 776., 727., 691., 710., 690., 1039., 734., 541., 855.,
    767., 839., 913., 1016., 800., 689., 931., 758., 946., 908., 1005., 639., 759., 794., 916.,
    857., 894., 763., 745.,
];

fn sum(v: impl IntoIterator<Item = f64>) -> f64 {
    v.into_iter().sum()
}

// --- Linear regression and correlation ------------------------------------------------

#[test]
fn example_4a_comet_wild_2() {
    let (x, y) = (&WILD2_X, &WILD2_Y);
    // The sums printed by Meeus (checks the data)
    assert_close(sum(x.iter().copied()), 4.2805, 5e-5);
    assert_close(sum(y.iter().copied()), 192.0400, 5e-5);
    assert_close(sum(x.iter().map(|v| v * v)), 1.0031, 5e-5);
    assert_close(sum(x.iter().zip(y).map(|(a, b)| a * b)), 43.7943, 5e-5);

    // y = 13.67 x + 7.03, i.e. m = 7.03 + 5 log Δ + 13.67 log r
    let fit = fitting::linear_regression(x, y).unwrap();
    assert_close(fit.slope, 13.67, 5e-3);
    assert_close(fit.intercept, 7.03, 5e-3);
    assert_close(fit.value(0.25), 13.67 * 0.25 + 7.03, 5e-3);
}

#[test]
fn example_4b_sunspot_rise_time() {
    let (x, y) = (&RISE_X, &RISE_Y);
    assert_close(sum(x.iter().copied()), 1120.0, 1e-9);
    assert_close(sum(y.iter().copied()), 2578.9, 1e-9);
    assert_close(sum(x.iter().map(|v| v * v)), 60608.0, 1e-9);
    assert_close(sum(y.iter().map(|v| v * v)), 340225.91, 1e-6);
    assert_close(sum(x.iter().zip(y).map(|(a, b)| a * b)), 122337.1, 1e-6);

    // Eq. 4.4: y = 244.18 − 2.49 x, and r = −0.767
    let fit = fitting::linear_regression(x, y).unwrap();
    assert_close(fit.intercept, 244.18, 5e-3);
    assert_close(fit.slope, -2.49, 5e-3);
    assert_close(fitting::correlation(x, y).unwrap(), -0.767, 5e-4);
}

#[test]
fn exercise_sunspots_and_rainfall_uncorrelated() {
    // Meeus, chapter 4 exercise: r = −0.054 for table 4.C
    assert_close(sum(SUNSPOTS.iter().copied()), 5865.7, 1e-9);
    assert_close(sum(RAINFALL.iter().copied()), 78119.0, 1e-9);
    let r = fitting::correlation(&SUNSPOTS, &RAINFALL).unwrap();
    assert_close(r, -0.054, 5e-4);
}

#[test]
fn correlation_zero_for_symmetric_parabola() {
    // Meeus, figure 6: seven points exactly on y = 2 − 2x − x², yet r = 0
    let x = [-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0];
    let y = [-6.0, -1.0, 2.0, 3.0, 2.0, -1.0, -6.0];
    for (&a, &b) in x.iter().zip(&y) {
        assert_eq!(2.0 - 2.0 * a - a * a, b);
    }
    assert_close(fitting::correlation(&x, &y).unwrap(), 0.0, 1e-15);
    // ...while the quadratic fit recovers the parabola exactly
    let q = fitting::quadratic_fit(&x, &y).unwrap();
    assert_close(q.a, -1.0, 1e-13);
    assert_close(q.b, -2.0, 1e-13);
    assert_close(q.c, 2.0, 1e-13);
}

#[test]
fn perfect_linear_relations() {
    let x = [1.0, 2.0, 3.0, 4.0, 5.0];
    let up = x.map(|v| 3.0 * v - 2.0);
    let down = x.map(|v| -0.5 * v + 7.0);
    assert_close(fitting::correlation(&x, &up).unwrap(), 1.0, 2e-16);
    assert_close(fitting::correlation(&x, &down).unwrap(), -1.0, 2e-16);
    // Never outside [-1, 1], even with rounding
    for scale in [1e-3, 0.1, 7.0, 1e6] {
        let xs = x.map(|v| v * scale + 0.1);
        let r = fitting::correlation(&xs, &xs.map(|v| 3.0 * v + 0.7)).unwrap();
        assert!(r <= 1.0 && (r - 1.0).abs() < 1e-15);
        let r = fitting::correlation(&xs, &xs.map(|v| -0.3 * v)).unwrap();
        assert!(r >= -1.0 && (r + 1.0).abs() < 1e-15);
    }
    let fit = fitting::linear_regression(&x, &up).unwrap();
    assert_close(fit.slope, 3.0, 1e-14);
    assert_close(fit.intercept, -2.0, 1e-14);
    assert_eq!(
        fitting::linear_regression(&[0.0, 2.0], &[1.0, 5.0]),
        Some(LinearFit {
            slope: 2.0,
            intercept: 1.0
        })
    );
}

#[test]
fn correlation_invariants() {
    let (x, y) = (&RISE_X, &RISE_Y);
    let r = fitting::correlation(x, y).unwrap();
    // Symmetric in x and y
    assert_close(fitting::correlation(y, x).unwrap(), r, 1e-15);
    // Unchanged by a change of units or origin; sign flips with the scale
    let y2: Vec<f64> = y.iter().map(|v| 1000.0 + 0.01 * v).collect();
    assert_close(fitting::correlation(x, &y2).unwrap(), r, 1e-12);
    let y3: Vec<f64> = y.iter().map(|v| -4.0 * v).collect();
    assert_close(fitting::correlation(x, &y3).unwrap(), -r, 1e-15);
    // The slopes of y on x and of x on y multiply to r²
    let a = fitting::linear_regression(x, y).unwrap().slope;
    let b = fitting::linear_regression(y, x).unwrap().slope;
    assert_close(a * b, r * r, 1e-14);
    // The fitted line passes through the centroid
    let fit = fitting::linear_regression(x, y).unwrap();
    let n = x.len() as f64;
    assert_close(
        fit.value(sum(x.iter().copied()) / n),
        sum(y.iter().copied()) / n,
        1e-12,
    );
}

#[test]
fn regression_with_large_abscissae() {
    // Julian-day-like abscissae. Evaluated literally, the raw sums of eq. 4.2
    // give the fitted value at 2451547.0 as 3.50030 instead of 3.5 (the
    // intercept is the difference of two numbers near 10¹⁹); sums about the
    // mean are good to 10⁻⁹.
    let x: Vec<f64> = (0..10).map(|i| 2_451_545.0 + i as f64 * 0.5).collect();
    let y: Vec<f64> = x.iter().map(|v| 0.25 * (v - 2_451_545.0) + 3.0).collect();
    let fit = fitting::linear_regression(&x, &y).unwrap();
    assert_close(fit.slope, 0.25, 1e-12);
    assert_close(fit.value(2_451_547.0), 3.5, 1e-9);
    assert_close(fitting::correlation(&x, &y).unwrap(), 1.0, 1e-14);
}

#[test]
fn degenerate_data() {
    // Mismatched lengths, too few points, too few distinct values
    assert_eq!(fitting::linear_regression(&[1.0, 2.0], &[1.0]), None);
    assert_eq!(fitting::linear_regression(&[], &[]), None);
    assert_eq!(fitting::linear_regression(&[1.0], &[1.0]), None);
    // Repeated x that does not average back exactly (0.1 + 0.1 + 0.1)
    assert_eq!(
        fitting::linear_regression(&[0.1, 0.1, 0.1], &[1.0, 2.0, 3.0]),
        None
    );
    assert_eq!(
        fitting::correlation(&[0.1, 0.1, 0.1], &[1.0, 2.0, 3.0]),
        None
    );
    assert_eq!(
        fitting::correlation(&[1.0, 2.0, 3.0], &[5.0, 5.0, 5.0]),
        None
    );
    assert_eq!(fitting::correlation(&[1.0, 2.0], &[1.0]), None);
    assert_eq!(
        fitting::quadratic_fit(&[1.0, 2.0, 1.0, 2.0], &[1.0, 2.0, 3.0, 4.0]),
        None
    );
    assert_eq!(fitting::quadratic_fit(&[1.0, 2.0, 3.0], &[1.0, 2.0]), None);
    // NaN data give None rather than NaN coefficients
    assert_eq!(
        fitting::linear_regression(&[1.0, f64::NAN], &[1.0, 2.0]),
        None
    );
    assert_eq!(
        fitting::quadratic_fit(&[1.0, 2.0, 3.0], &[1.0, f64::NAN, 2.0]),
        None
    );
}

// --- Quadratic and general least-squares fits ------------------------------------------

#[test]
fn quadratic_fit_matches_three_function_fit() {
    // Meeus: with f₀ = x², f₁ = x, f₂ = 1, eq. 4.7 reduces to eqs. 4.5–4.6
    let q = fitting::quadratic_fit(&RISE_X, &RISE_Y).unwrap();
    let g = fitting::fit_three_functions(&RISE_X, &RISE_Y, |x| x * x, |x| x, |_| 1.0).unwrap();
    assert_close(q.a, g.a, 1e-9 * g.a.abs());
    assert_close(q.b, g.b, 1e-9 * g.b.abs());
    assert_close(q.c, g.c, 1e-9 * g.c.abs());
    // The parabola fits at least as well as the straight line
    let line = fitting::linear_regression(&RISE_X, &RISE_Y).unwrap();
    let ss = |f: &dyn Fn(f64) -> f64| -> f64 {
        RISE_X
            .iter()
            .zip(&RISE_Y)
            .map(|(&x, &y)| (y - f(x)).powi(2))
            .sum()
    };
    assert!(ss(&|x| q.value(x)) <= ss(&|x| line.value(x)));
}

#[test]
fn quadratic_fit_least_squares_property() {
    // Perturbing any coefficient of the fit increases the sum of squares
    let q = fitting::quadratic_fit(&SUNSPOTS, &RAINFALL).unwrap();
    let ss = |p: QuadraticFit| -> f64 {
        SUNSPOTS
            .iter()
            .zip(&RAINFALL)
            .map(|(&x, &y)| (y - p.value(x)).powi(2))
            .sum()
    };
    let best = ss(q);
    for (da, db, dc) in [(1e-5, 0.0, 0.0), (0.0, 1e-3, 0.0), (0.0, 0.0, 0.1)] {
        for s in [-1.0, 1.0] {
            let p = QuadraticFit {
                a: q.a + s * da,
                b: q.b + s * db,
                c: q.c + s * dc,
            };
            assert!(ss(p) > best);
        }
    }
}

#[test]
fn quadratic_fit_exact_parabola_in_years() {
    // Abscissae around 1950, where raw power sums are badly conditioned
    let x: Vec<f64> = (1901..=1996).map(f64::from).collect();
    let p = |t: f64| 0.002 * (t - 1950.0).powi(2) - 0.3 * (t - 1950.0) + 12.0;
    let y: Vec<f64> = x.iter().map(|&t| p(t)).collect();
    let q = fitting::quadratic_fit(&x, &y).unwrap();
    for &t in &[1900.0, 1950.0, 1975.5, 2000.0] {
        assert_close(q.value(t), p(t), 1e-8);
    }
    assert_close(q.a, 0.002, 1e-13);
}

#[test]
fn example_4c_three_sines() {
    // Meeus, example 4.c: y = a sin x + b sin 2x + c sin 3x; the data were
    // computed from a = 1.2, b = −0.77, c = 0.39 and rounded to 4 decimals
    let x = [
        3., 20., 34., 50., 75., 88., 111., 129., 143., 160., 183., 200., 218., 230., 248., 269.,
        290., 303., 320., 344.,
    ];
    let y = [
        0.0433, 0.2532, 0.3386, 0.3560, 0.4983, 0.7577, 1.4585, 1.8628, 1.8264, 1.2431, -0.2043,
        -1.2431, -1.8422, -1.8726, -1.4889, -0.8372, -0.4377, -0.3640, -0.3508, -0.2126,
    ];
    let fit = fitting::fit_three_functions(
        &x,
        &y,
        |d: f64| d.to_radians().sin(),
        |d: f64| (2.0 * d).to_radians().sin(),
        |d: f64| (3.0 * d).to_radians().sin(),
    )
    .unwrap();
    assert_close(fit.a, 1.2, 1e-4);
    assert_close(fit.b, -0.77, 1e-4);
    assert_close(fit.c, 0.39, 1e-4);
}

#[test]
fn three_function_fit_degenerate() {
    let x = [1.0, 2.0, 3.0, 4.0];
    let y = [2.0, 4.0, 6.0, 8.0];
    // f₁ = 2 f₀: linearly dependent functions
    assert_eq!(
        fitting::fit_three_functions(&x, &y, |t| t, |t| 2.0 * t, |_| 1.0),
        None
    );
    // Too few points
    assert_eq!(
        fitting::fit_three_functions(&x[..2], &y[..2], |t| t * t, |t| t, |_| 1.0),
        None
    );
    // Mismatched lengths
    assert_eq!(
        fitting::fit_three_functions(&x, &y[..3], |t| t * t, |t| t, |_| 1.0),
        None
    );
    // A function that vanishes everywhere
    assert_eq!(
        fitting::fit_three_functions(&x, &y, |_| 0.0, |t| t, |_| 1.0),
        None
    );
}

#[test]
fn example_4d_square_root() {
    // Meeus, example 4.d: y = a √x through six points; Σ y√x = 15.2437,
    // Σ x = 15, a = 1.016
    let x = [0.0_f64, 1.0, 2.0, 3.0, 4.0, 5.0];
    let y = [0.0, 1.2, 1.4, 1.7, 2.1, 2.2];
    assert_close(
        sum(x.iter().zip(&y).map(|(a, b)| a.sqrt() * b)),
        15.2437,
        5e-5,
    );
    let a = fitting::fit_one_function(&x, &y, f64::sqrt).unwrap();
    assert_close(a, 15.2437 / 15.0, 5e-6);
    assert_close(a, 1.016, 5e-4);
}

#[test]
fn one_function_fit_edge_cases() {
    // Exact proportionality
    let x = [1.0, 2.0, 3.0];
    assert_close(
        fitting::fit_one_function(&x, &[3.0, 6.0, 9.0], |t| t).unwrap(),
        3.0,
        1e-15,
    );
    assert_eq!(fitting::fit_one_function(&x, &[1.0, 2.0], |t| t), None);
    assert_eq!(fitting::fit_one_function(&[], &[], |t| t), None);
    assert_eq!(
        fitting::fit_one_function(&x, &[1.0, 2.0, 3.0], |_| 0.0),
        None
    );
}

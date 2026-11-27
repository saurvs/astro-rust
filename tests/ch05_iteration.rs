//! Meeus, chapter 5 (iteration): fixed-point iteration, Newton's method and
//! bisection in `astro::iteration`, plus the iterations of chapter 3.

mod common;

use astro::iteration::{self, IterationError, Solution};
use common::assert_close;

/// The estimate after exactly `steps` steps of a fixed-point iteration that
/// has not converged yet.
fn fixed_point_iterate(g: impl FnMut(f64) -> f64, x0: f64, steps: u32) -> f64 {
    match iteration::fixed_point(g, x0, 0.0, steps) {
        Err(IterationError::NoConvergence { last, iterations }) => {
            assert_eq!(iterations, steps);
            last
        }
        other => panic!("expected no convergence after {steps} steps, got {other:?}"),
    }
}

/// The estimate after exactly `steps` Newton steps that have not converged.
fn newton_iterate(
    f: impl FnMut(f64) -> f64,
    df: impl FnMut(f64) -> f64,
    x0: f64,
    steps: u32,
) -> f64 {
    match iteration::newton(f, df, x0, 0.0, steps) {
        Err(IterationError::NoConvergence { last, iterations }) => {
            assert_eq!(iterations, steps);
            last
        }
        other => panic!("expected no convergence after {steps} steps, got {other:?}"),
    }
}

// --- Fixed-point iteration ----------------------------------------------------

#[test]
fn example_5a_square_root() {
    // Meeus, example 5.a: √159 by n → (n + 159/n)/2, starting from 12
    let g = |n: f64| (n + 159.0 / n) / 2.0;
    assert_close(fixed_point_iterate(g, 12.0, 1), 12.625, 0.0);
    assert_close(fixed_point_iterate(g, 12.0, 2), 12.6095297, 5e-8);
    assert_close(fixed_point_iterate(g, 12.0, 3), 12.60952021, 5e-9);
    let root = iteration::fixed_point(g, 12.0, 1e-10, 20).unwrap();
    assert_close(root.x, 12.60952021, 5e-9);
    assert_close(root.x, 159f64.sqrt(), 1e-14);
}

#[test]
fn example_5b_fifth_degree() {
    // Meeus, example 5.b: x⁵ + 17x − 8 = 0 as x = (8 − x⁵)/17
    let g = |x: f64| (8.0 - x.powi(5)) / 17.0;
    assert_close(fixed_point_iterate(g, 0.0, 1), 0.470588235, 5e-10);
    assert_close(fixed_point_iterate(g, 0.0, 2), 0.469230684, 5e-10);
    let root = iteration::fixed_point(g, 0.0, 1e-10, 50).unwrap();
    assert_close(root.x, 0.469249878, 5e-10);
    assert_close(root.x.powi(5) + 17.0 * root.x - 8.0, 0.0, 1e-9);

    // The other form suggested by Meeus, x = 8/(x⁴ + 17), gives the same root
    let other = iteration::fixed_point(|x| 8.0 / (x.powi(4) + 17.0), 0.0, 1e-12, 50).unwrap();
    assert_close(other.x, root.x, 1e-10);
}

#[test]
fn example_5c_diverges() {
    // Meeus, example 5.c: x⁵ + 3x − 8 = 0 as x = (8 − x⁵)/3 diverges
    let g = |x: f64| (8.0 - x.powi(5)) / 3.0;
    assert_close(fixed_point_iterate(g, 0.0, 1), 2.6667, 5e-5);
    assert_close(fixed_point_iterate(g, 0.0, 2), -42.2826, 5e-5);
    assert_close(fixed_point_iterate(g, 0.0, 3), 45049099.0, 0.5);
    assert_close(fixed_point_iterate(g, 0.0, 4) / 1e37, -6.18, 5e-3);
    // Then x⁵ overflows: the iteration stops instead of producing garbage
    assert!(matches!(
        iteration::fixed_point(g, 0.0, 1e-10, 100),
        Err(IterationError::NonFinite { .. })
    ));
}

#[test]
fn example_5d_rewritten_form_converges() {
    // Meeus, example 5.d: x = (8 − 3x)^(1/5), root reached after 14 steps
    let g = |x: f64| (8.0 - 3.0 * x).powf(0.2);
    let root = iteration::fixed_point(g, 0.0, 1e-9, 100).unwrap();
    assert_close(root.x, 1.321785627, 5e-10);
    assert_eq!(root.iterations, 14);
}

#[test]
fn oscillating_iteration_is_reported() {
    // Meeus, chapter 5: x = 8/(x⁴ + 3) jumps between two values forever
    let g = |x: f64| 8.0 / (x.powi(4) + 3.0);
    match iteration::fixed_point(g, 0.0, 1e-10, 100) {
        Err(IterationError::NoConvergence { last, iterations }) => {
            assert_eq!(iterations, 100);
            assert_close(last, 0.149436927, 5e-10);
        }
        other => panic!("expected no convergence, got {other:?}"),
    }
    assert_close(fixed_point_iterate(g, 0.0, 101), 2.666223459, 5e-10);
}

#[test]
fn sine_cosine_equations() {
    // Meeus, chapter 5: sin φ = 3 cos φ. Iterating sin φ = 3 cos φ directly
    // asks for arcsin 3, an impossibility...
    let bad = iteration::fixed_point(|p: f64| (3.0 * p.cos()).asin(), 0.0, 1e-12, 100);
    assert!(matches!(
        bad,
        Err(IterationError::NonFinite { iterations: 1 })
    ));
    // ...while cos φ = (sin φ)/3 converges to 71.565051°
    let good = iteration::fixed_point(|p: f64| (p.sin() / 3.0).acos(), 0.0, 1e-12, 200).unwrap();
    assert_close(good.x.to_degrees(), 71.565051, 5e-7);
    assert_close(good.x.tan(), 3.0, 1e-10);

    // sin φ = cos 2φ from 29°: φ = arcsin(cos 2φ) diverges, but
    // 2φ = arccos(sin φ) converges to 30°
    let start = 29f64.to_radians();
    let bad = iteration::fixed_point(|p: f64| (2.0 * p).cos().asin(), start, 1e-12, 200);
    assert!(bad.is_err());
    let good = iteration::fixed_point(|p: f64| p.sin().acos() / 2.0, start, 1e-13, 200).unwrap();
    assert_close(good.x.to_degrees(), 30.0, 1e-9);
}

// --- Newton's method ------------------------------------------------------------

#[test]
fn newton_bad_starting_value_oscillates() {
    // Meeus, chapter 5: x⁵ − 3x − 8 = 0 from x = 0. The function has a
    // maximum near x = −0.88 and the estimates jump about without settling.
    let f = |x: f64| x.powi(5) - 3.0 * x - 8.0;
    let df = |x: f64| 5.0 * x.powi(4) - 3.0;
    let expected = [
        -2.666666667,
        -2.126929222,
        -1.672392941,
        -1.227532073,
        -0.376965299,
        -2.749036974,
        -2.194266642,
        -1.731201846,
        -1.293218529,
        -0.588844800,
        -3.216865068,
        -2.572967056,
        -2.049930312,
        -1.603831481,
        -1.145086796,
    ];
    for (k, &x) in expected.iter().enumerate() {
        assert_close(newton_iterate(f, df, 0.0, k as u32 + 1), x, 5e-9);
    }
}

#[test]
fn newton_good_starting_value_converges() {
    // Meeus, chapter 5: the same equation from x = 1 converges to
    // 1.670579156 (9 decimals) after 11 steps
    let f = |x: f64| x.powi(5) - 3.0 * x - 8.0;
    let df = |x: f64| 5.0 * x.powi(4) - 3.0;
    assert_close(newton_iterate(f, df, 1.0, 1), 6.0, 0.0);
    assert_close(newton_iterate(f, df, 1.0, 2), 4.803458391, 5e-10);
    assert_close(newton_iterate(f, df, 1.0, 10), 1.670579511, 5e-10);
    assert_close(newton_iterate(f, df, 1.0, 11), 1.670579156, 5e-10);
    let root = iteration::newton(f, df, 1.0, 1e-9, 50).unwrap();
    assert_close(root.x, 1.670579156, 5e-10);
    assert_eq!(root.iterations, 12); // the 12th step confirms the 11th
}

#[test]
fn smaller_than_test_stops_oscillation() {
    // Meeus, chapter 5, "Test on smaller than": s³ + 3s − W = 0 solved by
    // s → (2s³ + W) / (3(s² + 1)) from s = 0. On the HP-85, W = 1.5 made the
    // last digit oscillate for ever. Here even a zero tolerance terminates,
    // because estimates a few ulps apart count as equal.
    let g = |w: f64| move |s: f64| (2.0 * s.powi(3) + w) / (3.0 * (s * s + 1.0));

    let s = iteration::fixed_point(g(0.9), 0.0, 0.0, 100).unwrap();
    assert_close(fixed_point_iterate(g(0.9), 0.0, 1), 0.3, 1e-15);
    assert_close(fixed_point_iterate(g(0.9), 0.0, 2), 0.291743119266, 5e-13);
    // Meeus prints 0.291724443548 from the 12-digit HP-85; the exact root is
    // 0.29172444354709, so the book is one unit off in the last digit.
    assert_close(s.x, 0.291724443548, 1e-12);
    assert_close(s.x, 0.29172444354708566, 1e-16);

    let s = iteration::fixed_point(g(1.5), 0.0, 0.0, 100).unwrap();
    assert_close(fixed_point_iterate(g(1.5), 0.0, 2), 0.466666666667, 5e-13);
    assert_close(fixed_point_iterate(g(1.5), 0.0, 3), 0.466220600162, 5e-13);
    assert_close(s.x, 0.466220523910, 1e-12);
    assert_close(s.x, 0.46622052391077343, 1e-16);
    assert!(s.iterations < 10);

    // The iteration is Newton's method for s³ + 3s − W
    let newton = iteration::newton(
        |s: f64| s.powi(3) + 3.0 * s - 1.5,
        |s: f64| 3.0 * s * s + 3.0,
        0.0,
        0.0,
        100,
    )
    .unwrap();
    assert_close(newton.x, s.x, 1e-15);
}

#[test]
fn newton_zero_derivative_and_no_real_root() {
    // f = x² + 1 has no real root, and f' = 0 at the start
    let f = |x: f64| x * x + 1.0;
    let df = |x: f64| 2.0 * x;
    assert_eq!(
        iteration::newton(f, df, 0.0, 1e-12, 50),
        Err(IterationError::ZeroDerivative {
            x: 0.0,
            iterations: 0
        })
    );
    // From elsewhere it wanders without converging
    assert!(iteration::newton(f, df, 0.5, 1e-12, 200).is_err());
}

#[test]
fn newton_exact_root_at_start() {
    let sol = iteration::newton(|x| x - 2.0, |_| 1.0, 2.0, 1e-12, 10).unwrap();
    assert_eq!(
        sol,
        Solution {
            x: 2.0,
            iterations: 1
        }
    );
}

// --- Chapter 3 iterations, run through the generic solvers ------------------------

#[test]
fn example_3c_zero_by_eq_3_6() {
    // Meeus, example 3.c: n₀ = −2y₂ / (a + b + c n₀) for the declination of
    // Mercury, in arcseconds: y = −1693.4, +406.3, +2303.2
    let (y2, a, b) = (406.3, 2099.7, 1896.9);
    let c = b - a;
    let g = |n: f64| -2.0 * y2 / (a + b + c * n);
    assert_close(fixed_point_iterate(g, 0.0, 1), -0.20332, 5e-6);
    assert_close(fixed_point_iterate(g, 0.0, 2), -0.20125, 5e-6);
    assert_close(fixed_point_iterate(g, 0.0, 3), -0.20127, 5e-6);
    let n0 = iteration::fixed_point(g, 0.0, 1e-12, 100).unwrap().x;
    let direct = astro::interpol::zero_three_values(-1693.4, 406.3, 2303.2).unwrap();
    assert_close(n0, direct, 1e-12);
}

#[test]
fn example_3d_eq_3_6_versus_eq_3_7() {
    // Meeus, example 3.d: y = 3 + 2x − 3x² tabulated at x = −1, 0, +1
    let (y1, y2, y3) = (-2.0, 3.0, 2.0);
    let (a, b) = (y2 - y1, y3 - y2);
    let c = b - a;

    // Eq. 3.6 converges slowly (Meeus prints these truncated, not rounded)
    let g = |n: f64| -2.0 * y2 / (a + b + c * n);
    let slow = [-1.5, -0.461538, -0.886363, -0.643902, -0.763027, -0.699450];
    for (k, &n) in slow.iter().enumerate() {
        assert_close(fixed_point_iterate(g, 0.0, k as u32 + 1), n, 1e-6);
    }
    // ...needing 24 steps before the sixth decimal is right for good
    let rounded = |k: u32| astro::util::round_upto_digits(fixed_point_iterate(g, 0.0, k), 6);
    assert_eq!(rounded(23), -0.72076);
    for k in 24..=40 {
        assert_eq!(rounded(k), -0.720759);
    }

    // Eq. 3.7 is Newton's method on 2y₂ + n(a + b + c n) and is fast
    let f = |n: f64| 2.0 * y2 + n * (a + b + c * n);
    let df = |n: f64| a + b + 2.0 * c * n;
    let fast = [
        -1.5,
        -0.886363636364,
        -0.732001693959,
        -0.720818540935,
        -0.720759221726,
        -0.720759220056,
    ];
    // Meeus' values come from a 12-digit machine; two of them differ from
    // the double-precision iterates by one unit in the 12th decimal.
    for (k, &n) in fast.iter().enumerate() {
        assert_close(newton_iterate(f, df, 0.0, k as u32 + 1), n, 1e-12);
    }

    let n0 = iteration::newton(f, df, 0.0, 1e-13, 50).unwrap().x;
    assert_close(
        n0,
        astro::interpol::zero_three_values(y1, y2, y3).unwrap(),
        1e-15,
    );
    assert_close(
        iteration::fixed_point(g, 0.0, 1e-14, 500).unwrap().x,
        n0,
        1e-12,
    );
}

// --- Bisection --------------------------------------------------------------------

#[test]
fn binary_search_example() {
    // Meeus, chapter 5, "The binary search": x(x⁴ + 17) − 8 on [0, 1]
    let f = |x: f64| x * (x.powi(4) + 17.0) - 8.0;
    assert_eq!(f(0.0), -8.0);
    assert_eq!(f(1.0), 10.0);
    assert_close(f(0.5), 0.53125, 0.0);

    // After 32 halvings the bracket is 2^-32 wide, so its midpoint is within
    // 2^-33 of the root: nine exact decimals, as Meeus says.
    let tolerance = 2f64.powi(-33);
    let root = iteration::bisection(f, 0.0, 1.0, tolerance).unwrap();
    assert_eq!(root.iterations, 32);
    assert_close(root.x, 0.4692498784547387, tolerance);
    assert_close(root.x, 0.469249878, 1e-9);

    // Bracket ends in either order give the same result
    let reversed = iteration::bisection(f, 1.0, 0.0, tolerance).unwrap();
    assert_eq!(reversed, root);
}

#[test]
fn bisection_full_precision() {
    let f = |x: f64| x * (x.powi(4) + 17.0) - 8.0;
    let root = iteration::bisection(f, 0.0, 1.0, 0.0).unwrap();
    let newton = iteration::newton(f, |x| 5.0 * x.powi(4) + 17.0, 0.5, 0.0, 50).unwrap();
    assert!((root.x - newton.x).abs() <= 2.0 * f64::EPSILON);
    assert!(root.iterations <= 60);
}

#[test]
fn bisection_root_at_bracket_end_or_midpoint() {
    let f = |x: f64| x - 1.0;
    assert_eq!(
        iteration::bisection(f, 1.0, 3.0, 1e-9).unwrap(),
        Solution {
            x: 1.0,
            iterations: 0
        }
    );
    assert_eq!(
        iteration::bisection(f, -3.0, 1.0, 1e-9).unwrap(),
        Solution {
            x: 1.0,
            iterations: 0
        }
    );
    assert_eq!(
        iteration::bisection(f, 0.0, 2.0, 1e-9).unwrap(),
        Solution {
            x: 1.0,
            iterations: 1
        }
    );
}

#[test]
fn bisection_needs_a_sign_change() {
    let f = |x: f64| x * x + 1.0;
    assert_eq!(
        iteration::bisection(f, -1.0, 1.0, 1e-9),
        Err(IterationError::NoSignChange)
    );
    assert_eq!(
        iteration::bisection(|_| f64::NAN, -1.0, 1.0, 1e-9),
        Err(IterationError::NoSignChange)
    );
}

#[test]
fn bisection_converges_to_a_discontinuity() {
    // Meeus' footnote: tan 86° > 0 and tan 93° < 0, but tan has no zero
    // there; bisection homes in on the pole at 90° instead.
    let sol = iteration::bisection(|x: f64| x.to_radians().tan(), 86.0, 93.0, 1e-9).unwrap();
    assert_close(sol.x, 90.0, 1e-9);
    assert!(sol.x.to_radians().tan().abs() > 1e6);
}

#[test]
fn bisection_nan_inside_bracket() {
    let f = |x: f64| {
        if x > 0.2 && x < 0.8 {
            f64::NAN
        } else {
            x - 0.5
        }
    };
    assert!(matches!(
        iteration::bisection(f, 0.0, 1.0, 1e-9),
        Err(IterationError::NonFinite { iterations: 1 })
    ));
}

#[test]
fn bisection_huge_bracket() {
    let sol = iteration::bisection(|x| x - 3.0, -f64::MAX, f64::MAX, 1e-9).unwrap();
    assert_close(sol.x, 3.0, 1e-9);
}

// --- Invalid input and error reporting ------------------------------------------

#[test]
fn invalid_inputs() {
    let g = |x: f64| x / 2.0;
    for tol in [f64::NAN, -1.0] {
        assert_eq!(
            iteration::fixed_point(g, 1.0, tol, 10),
            Err(IterationError::InvalidInput)
        );
        assert_eq!(
            iteration::newton(g, |_| 0.5, 1.0, tol, 10),
            Err(IterationError::InvalidInput)
        );
        assert_eq!(
            iteration::bisection(g, -1.0, 1.0, tol),
            Err(IterationError::InvalidInput)
        );
    }
    for x0 in [f64::NAN, f64::INFINITY] {
        assert_eq!(
            iteration::fixed_point(g, x0, 1e-9, 10),
            Err(IterationError::InvalidInput)
        );
        assert_eq!(
            iteration::newton(g, |_| 0.5, x0, 1e-9, 10),
            Err(IterationError::InvalidInput)
        );
        assert_eq!(
            iteration::bisection(g, x0, 1.0, 1e-9),
            Err(IterationError::InvalidInput)
        );
    }
}

#[test]
fn zero_step_limit() {
    assert_eq!(
        iteration::fixed_point(|x| x, 1.0, 1e-9, 0),
        Err(IterationError::NoConvergence {
            last: 1.0,
            iterations: 0
        })
    );
    assert_eq!(
        iteration::newton(|x| x, |_| 1.0, 1.0, 1e-9, 0),
        Err(IterationError::NoConvergence {
            last: 1.0,
            iterations: 0
        })
    );
}

#[test]
fn non_finite_function_values() {
    assert!(matches!(
        iteration::fixed_point(|_| f64::NAN, 1.0, 1e-9, 10),
        Err(IterationError::NonFinite { iterations: 1 })
    ));
    assert!(matches!(
        iteration::newton(|x| x, |_| f64::INFINITY, 1.0, 1e-9, 10),
        Err(IterationError::NonFinite { iterations: 1 })
    ));
    assert!(matches!(
        iteration::newton(|_| f64::NAN, |_| 1.0, 1.0, 1e-9, 10),
        Err(IterationError::NonFinite { iterations: 1 })
    ));
}

#[test]
fn errors_display_and_box() {
    let err = iteration::fixed_point(|x| 8.0 / (x.powi(4) + 3.0), 0.0, 1e-10, 10).unwrap_err();
    assert!(err.to_string().contains("did not converge in 10 steps"));
    let boxed: Box<dyn std::error::Error> = Box::new(err);
    assert!(!boxed.to_string().is_empty());
    for e in [
        IterationError::NonFinite { iterations: 3 },
        IterationError::ZeroDerivative {
            x: 1.0,
            iterations: 2,
        },
        IterationError::NoSignChange,
        IterationError::InvalidInput,
    ] {
        assert!(!e.to_string().is_empty());
    }
}

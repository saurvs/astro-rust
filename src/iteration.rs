//! Iterative solution of equations: fixed-point iteration, Newton's method and bisection (Meeus, chapter 5).
//!
//! Each solver repeats a step until two successive estimates agree to within a
//! tolerance chosen by the caller, and gives up with an
//! [`IterationError`] if it cannot get there. Two safeguards from Meeus,
//! *Astronomical Algorithms*, 2nd ed., chapters 1 and 5, are built in:
//!
//! * a step limit, so that a diverging or oscillating iteration stops;
//! * a "smaller than" convergence test instead of an equality test. Two
//!   estimates count as equal when they differ by at most the tolerance *or*
//!   by at most a few units in the last place, so a tolerance that is too
//!   small for `f64` never makes an iteration oscillate forever in its
//!   last digit.
//!
//! [`bisection`] (Meeus' "binary search") cannot diverge: it only needs a
//! bracket on which the function changes sign.
//!
//! # Example
//!
//! The real root of x⁵ + 17x − 8 = 0 (Meeus, example 5.b), first by
//! rewriting the equation as x = (8 − x⁵) / 17, then by Newton's method and
//! by bisection:
//!
//! ```
//! use astro::iteration;
//!
//! let fixed = iteration::fixed_point(|x| (8.0 - x.powi(5)) / 17.0, 0.0, 1e-12, 100).unwrap();
//! assert!((fixed.x - 0.469249878).abs() < 1e-9);
//!
//! let f = |x: f64| x.powi(5) + 17.0 * x - 8.0;
//! let df = |x: f64| 5.0 * x.powi(4) + 17.0;
//! let newton = iteration::newton(f, df, 0.0, 1e-12, 100).unwrap();
//! assert!((newton.x - fixed.x).abs() < 1e-12);
//!
//! let bisect = iteration::bisection(f, 0.0, 1.0, 1e-12).unwrap();
//! assert!((bisect.x - fixed.x).abs() < 1e-12);
//! ```

use std::fmt;

/// A converged solution of an iteration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Solution {
    /// The solution: the last estimate, which agrees with the previous one
    /// to within the requested tolerance.
    pub x: f64,
    /// Number of iteration steps performed (evaluations of the iteration
    /// formula, or bisection halvings).
    pub iterations: u32,
}

/// Reasons why an iteration failed.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum IterationError {
    /// The estimates did not settle within the allowed number of steps. The
    /// iteration may be diverging slowly, oscillating, or merely converging
    /// too slowly for the step limit.
    NoConvergence {
        /// The last estimate computed.
        last: f64,
        /// Number of steps performed.
        iterations: u32,
    },
    /// An estimate or a function value became infinite or NaN: the
    /// iteration diverged or left the domain of the function.
    NonFinite {
        /// The step at which the non-finite value appeared.
        iterations: u32,
    },
    /// Newton's method hit a point where the derivative is exactly zero, so
    /// the next estimate is undefined.
    ZeroDerivative {
        /// The estimate at which the derivative vanished.
        x: f64,
        /// Number of steps performed before that point.
        iterations: u32,
    },
    /// Bisection was given a bracket on whose ends the function has the same
    /// sign (or is NaN), so a root is not known to lie inside it.
    NoSignChange,
    /// An argument was invalid: a NaN or negative tolerance, or a starting
    /// value or bracket end that is not finite.
    InvalidInput,
}

impl fmt::Display for IterationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            IterationError::NoConvergence { last, iterations } => write!(
                f,
                "iteration did not converge in {iterations} steps (last estimate {last})"
            ),
            IterationError::NonFinite { iterations } => {
                write!(
                    f,
                    "iteration produced a non-finite value at step {iterations}"
                )
            }
            IterationError::ZeroDerivative { x, iterations } => {
                write!(f, "zero derivative at x = {x} after {iterations} steps")
            }
            IterationError::NoSignChange => {
                write!(f, "function does not change sign over the bracket")
            }
            IterationError::InvalidInput => write!(f, "invalid input to iteration"),
        }
    }
}

impl std::error::Error for IterationError {}

/// Largest difference between two successive estimates that is accepted as
/// convergence: the caller's tolerance, but never less than four units in
/// the last place of `x`.
#[inline]
fn converged(x_new: f64, x_old: f64, tolerance: f64) -> bool {
    (x_new - x_old).abs() <= tolerance.max(4.0 * f64::EPSILON * x_new.abs())
}

#[inline]
fn valid_tolerance(tolerance: f64) -> bool {
    tolerance >= 0.0 && !tolerance.is_nan()
}

/// Solves `x = g(x)` by fixed-point iteration (Meeus, chapter 5).
///
/// Starting from `x0`, computes x₁ = g(x₀), x₂ = g(x₁), … until two
/// successive estimates differ by at most `tolerance` (or by at most a few
/// units in the last place, whichever is larger).
///
/// The iteration converges only if |g′(x)| < 1 near the solution, so the
/// way an equation is rewritten in the form x = g(x) matters: Meeus shows
/// that x = (8 − x⁵)/3 diverges for the root of x⁵ + 3x − 8 = 0, while
/// x = (8 − 3x)^(1/5) converges.
///
/// # Arguments
///
/// * `g`: The iteration function
/// * `x0`: Starting estimate (must be finite)
/// * `tolerance`: Absolute convergence tolerance, `>= 0`, in the units of `x`
/// * `max_iterations`: Maximum number of evaluations of `g`
///
/// # Returns
///
/// The converged [`Solution`], or
///
/// * [`IterationError::NoConvergence`] after `max_iterations` steps,
/// * [`IterationError::NonFinite`] if `g` returns an infinite or NaN value,
/// * [`IterationError::InvalidInput`] for a non-finite `x0` or a negative or
///   NaN `tolerance`.
pub fn fixed_point<G>(
    mut g: G,
    x0: f64,
    tolerance: f64,
    max_iterations: u32,
) -> Result<Solution, IterationError>
where
    G: FnMut(f64) -> f64,
{
    if !x0.is_finite() || !valid_tolerance(tolerance) {
        return Err(IterationError::InvalidInput);
    }

    let mut x = x0;
    for k in 1..=max_iterations {
        let x_new = g(x);
        if !x_new.is_finite() {
            return Err(IterationError::NonFinite { iterations: k });
        }
        if converged(x_new, x, tolerance) {
            return Ok(Solution {
                x: x_new,
                iterations: k,
            });
        }
        x = x_new;
    }

    Err(IterationError::NoConvergence {
        last: x,
        iterations: max_iterations,
    })
}

/// Solves `f(x) = 0` by Newton's method (Meeus, chapter 5).
///
/// Starting from `x0`, computes xₙ₊₁ = xₙ − f(xₙ)/f′(xₙ) until two
/// successive estimates differ by at most `tolerance` (or by at most a few
/// units in the last place, whichever is larger).
///
/// Convergence is fast near a simple root, but the starting value matters:
/// near a maximum or minimum of `f` the tangent is almost horizontal and the
/// estimates can jump about without settling (Meeus' example with
/// x⁵ − 3x − 8 = 0 started from x = 0).
///
/// # Arguments
///
/// * `f`: The function whose root is wanted
/// * `df`: Its derivative
/// * `x0`: Starting estimate (must be finite)
/// * `tolerance`: Absolute convergence tolerance, `>= 0`, in the units of `x`
/// * `max_iterations`: Maximum number of Newton steps
///
/// # Returns
///
/// The converged [`Solution`], or
///
/// * [`IterationError::NoConvergence`] after `max_iterations` steps,
/// * [`IterationError::ZeroDerivative`] if `df` is exactly zero at an
///   estimate where `f` is not,
/// * [`IterationError::NonFinite`] if `f`, `df` or an estimate becomes
///   infinite or NaN,
/// * [`IterationError::InvalidInput`] for a non-finite `x0` or a negative or
///   NaN `tolerance`.
pub fn newton<F, D>(
    mut f: F,
    mut df: D,
    x0: f64,
    tolerance: f64,
    max_iterations: u32,
) -> Result<Solution, IterationError>
where
    F: FnMut(f64) -> f64,
    D: FnMut(f64) -> f64,
{
    if !x0.is_finite() || !valid_tolerance(tolerance) {
        return Err(IterationError::InvalidInput);
    }

    let mut x = x0;
    for k in 1..=max_iterations {
        let fx = f(x);
        if !fx.is_finite() {
            return Err(IterationError::NonFinite { iterations: k });
        }
        if fx == 0.0 {
            return Ok(Solution { x, iterations: k });
        }
        let dfx = df(x);
        if !dfx.is_finite() {
            return Err(IterationError::NonFinite { iterations: k });
        }
        if dfx == 0.0 {
            return Err(IterationError::ZeroDerivative {
                x,
                iterations: k - 1,
            });
        }
        let x_new = x - fx / dfx;
        if !x_new.is_finite() {
            return Err(IterationError::NonFinite { iterations: k });
        }
        if converged(x_new, x, tolerance) {
            return Ok(Solution {
                x: x_new,
                iterations: k,
            });
        }
        x = x_new;
    }

    Err(IterationError::NoConvergence {
        last: x,
        iterations: max_iterations,
    })
}

/// Upper bound on bisection steps. Halving a finite interval of `f64`s
/// exhausts the representable values in it after at most about 2100 steps,
/// so this is never reached in practice; it only makes termination obvious.
const MAX_BISECTION_STEPS: u32 = 4096;

/// Finds a root of `f` in the bracket `[a, b]` by bisection, Meeus' "binary
/// search" (chapter 5).
///
/// The bracket is halved repeatedly, keeping the half on whose ends `f`
/// still has opposite signs, until its half-width is at most `tolerance`.
/// With `tolerance = 0.0` the search continues until the two ends are
/// adjacent floating-point numbers. The method cannot diverge: after `n`
/// steps the root is known to within `|b − a| / 2ⁿ⁺¹`.
///
/// `f` must be continuous on the bracket. If it is not, a sign change may
/// come from a discontinuity instead of a root (for example tan x between 86°
/// and 93°), and bisection then converges to the discontinuity.
///
/// # Arguments
///
/// * `f`: The function whose root is wanted
/// * `a`, `b`: Ends of the bracket, in either order (must be finite)
/// * `tolerance`: Absolute tolerance on the root, `>= 0`, in the units of `x`
///
/// # Returns
///
/// The [`Solution`]: the midpoint of the final bracket (or a point where `f`
/// is exactly zero), and the number of halvings. Errors:
///
/// * [`IterationError::NoSignChange`] if `f(a)` and `f(b)` have the same
///   sign or either is NaN,
/// * [`IterationError::NonFinite`] if `f` returns NaN inside the bracket,
/// * [`IterationError::InvalidInput`] for a non-finite bracket end or a
///   negative or NaN `tolerance`.
pub fn bisection<F>(mut f: F, a: f64, b: f64, tolerance: f64) -> Result<Solution, IterationError>
where
    F: FnMut(f64) -> f64,
{
    if !a.is_finite() || !b.is_finite() || !valid_tolerance(tolerance) {
        return Err(IterationError::InvalidInput);
    }

    let fa = f(a);
    if fa == 0.0 {
        return Ok(Solution {
            x: a,
            iterations: 0,
        });
    }
    let fb = f(b);
    if fb == 0.0 {
        return Ok(Solution {
            x: b,
            iterations: 0,
        });
    }
    if fa.is_nan() || fb.is_nan() || (fa < 0.0) == (fb < 0.0) {
        return Err(IterationError::NoSignChange);
    }

    let (mut lo, mut hi, lo_negative) = if a < b {
        (a, b, fa < 0.0)
    } else {
        (b, a, fb < 0.0)
    };

    for k in 1..=MAX_BISECTION_STEPS {
        // Halving each end separately cannot overflow.
        let mid = 0.5 * lo + 0.5 * hi;
        if mid <= lo || mid >= hi {
            // The ends are adjacent floating-point numbers.
            return Ok(Solution {
                x: mid,
                iterations: k - 1,
            });
        }
        let fm = f(mid);
        if fm.is_nan() {
            return Err(IterationError::NonFinite { iterations: k });
        }
        if fm == 0.0 {
            return Ok(Solution {
                x: mid,
                iterations: k,
            });
        }
        if (fm < 0.0) == lo_negative {
            lo = mid;
        } else {
            hi = mid;
        }
        if 0.5 * (hi - lo) <= tolerance {
            return Ok(Solution {
                x: 0.5 * lo + 0.5 * hi,
                iterations: k,
            });
        }
    }

    Err(IterationError::NoConvergence {
        last: 0.5 * lo + 0.5 * hi,
        iterations: MAX_BISECTION_STEPS,
    })
}

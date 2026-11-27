//! Curve fitting: linear regression, correlation and least-squares fits (Meeus, chapter 4).
//!
//! All fits minimise the sum of the squared vertical deviations
//! `Σ (yᵢ − f(xᵢ))²`; the `xs` are taken as exact and the `ys` as measured.
//! Data are given as two slices of equal length, `xs` and `ys`.
//!
//! * [`linear_regression`]: straight line `y = a x + b` (eq. 4.2)
//! * [`correlation`]: linear correlation coefficient (eq. 4.3)
//! * [`quadratic_fit`]: parabola `y = a x² + b x + c` (eqs. 4.5 and 4.6)
//! * [`fit_three_functions`]: `y = a f₀(x) + b f₁(x) + c f₂(x)` for any three
//!   chosen functions (eq. 4.7)
//! * [`fit_one_function`]: `y = a f(x)` (eq. 4.8)
//!
//! The results are the same as those of Meeus' formulas, but the sums for the
//! straight line, the correlation and the parabola are taken about the mean of
//! the data. Mathematically this changes nothing; numerically it avoids the
//! loss of precision that the raw sums `N Σx² − (Σx)²` suffer when the `xs`
//! are large compared with their spread (years, Julian days).
//!
//! Every function returns `None` when the fit is not defined: slices of
//! different lengths, too few points, or too few distinct `xs`.
//!
//! # Example
//!
//! The rise time of the solar cycle (months, x) against the height of the
//! following sunspot maximum (y), for the 22 maxima from 1761 to 1989
//! (Meeus, example 4.b):
//!
//! ```
//! use astro::fitting;
//!
//! let x = [73., 38., 35., 42., 78., 68., 74., 42., 52., 54., 39., 61., 42., 49., 50., 62., 44.,
//!          39., 43., 54., 44., 37.];
//! let y = [90.4, 125.3, 161.8, 143.4, 52.5, 50.8, 71.5, 152.8, 131.3, 98.5, 144.8, 78.1, 89.5,
//!          63.9, 112.1, 82.0, 119.8, 161.2, 208.4, 111.6, 167.1, 162.1];
//!
//! let line = fitting::linear_regression(&x, &y).unwrap();
//! assert!((line.slope - -2.49).abs() < 0.005);
//! assert!((line.intercept - 244.18).abs() < 0.005);
//!
//! let r = fitting::correlation(&x, &y).unwrap();
//! assert!((r - -0.767).abs() < 0.0005);
//! ```

/// A straight line `y = slope · x + intercept` fitted to data.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearFit {
    /// Slope `a` of the line, in units of y per unit of x.
    pub slope: f64,
    /// Value `b` of the line at `x = 0`, in units of y.
    pub intercept: f64,
}

impl LinearFit {
    /// Evaluates the fitted line at `x`.
    #[inline]
    pub fn value(&self, x: f64) -> f64 {
        self.slope * x + self.intercept
    }
}

/// A parabola `y = a x² + b x + c` fitted to data.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuadraticFit {
    /// Coefficient of `x²`.
    pub a: f64,
    /// Coefficient of `x`.
    pub b: f64,
    /// Constant term.
    pub c: f64,
}

impl QuadraticFit {
    /// Evaluates the fitted parabola at `x`.
    #[inline]
    pub fn value(&self, x: f64) -> f64 {
        (self.a * x + self.b) * x + self.c
    }
}

/// The coefficients of a fit `y = a f₀(x) + b f₁(x) + c f₂(x)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThreeFunctionFit {
    /// Coefficient of `f₀(x)`.
    pub a: f64,
    /// Coefficient of `f₁(x)`.
    pub b: f64,
    /// Coefficient of `f₂(x)`.
    pub c: f64,
}

/// Returns `true` if `values` holds at least `count` (at most 3) different
/// values.
fn has_distinct(values: &[f64], count: usize) -> bool {
    let mut found: [f64; 3] = [0.0; 3];
    let count = count.min(found.len());
    let mut n = 0;
    for &v in values {
        if !found[..n].contains(&v) {
            found[n] = v;
            n += 1;
            if n >= count {
                return true;
            }
        }
    }
    false
}

#[inline]
fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

/// Sums of squares and products of deviations from the means:
/// `(x̄, ȳ, Σ(x−x̄)², Σ(y−ȳ)², Σ(x−x̄)(y−ȳ))`.
fn centred_sums(xs: &[f64], ys: &[f64]) -> (f64, f64, f64, f64, f64) {
    let mx = mean(xs);
    let my = mean(ys);
    let (mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0);
    for (&x, &y) in xs.iter().zip(ys) {
        let dx = x - mx;
        let dy = y - my;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    (mx, my, sxx, syy, sxy)
}

/// Fits a straight line `y = a x + b` to data by least squares (linear
/// regression; Meeus, eqs. 4.1 and 4.2).
///
/// # Arguments
///
/// * `xs`: Values of the independent variable (taken as exact)
/// * `ys`: Measured values, same length as `xs`
///
/// # Returns
///
/// The fitted [`LinearFit`], or `None` if the slices have different
/// lengths, if there are fewer than two distinct `xs`, or if the result is
/// not finite.
pub fn linear_regression(xs: &[f64], ys: &[f64]) -> Option<LinearFit> {
    if xs.len() != ys.len() || !has_distinct(xs, 2) {
        return None;
    }
    let (mx, my, sxx, _, sxy) = centred_sums(xs, ys);
    let slope = sxy / sxx;
    let fit = LinearFit {
        slope,
        intercept: my - slope * mx,
    };
    (fit.slope.is_finite() && fit.intercept.is_finite()).then_some(fit)
}

/// Computes the linear correlation coefficient of data (Meeus, eq. 4.3).
///
/// The coefficient measures how well the points follow a straight line:
/// `+1` or `−1` for points exactly on a line rising or falling with `x`,
/// near `0` when there is no *linear* relation. A value near zero does not
/// rule out a strong non-linear relation: points on a parabola symmetric
/// about the middle of the data give exactly zero.
///
/// # Arguments
///
/// * `xs`: Values of the first variable
/// * `ys`: Values of the second variable, same length as `xs`
///
/// # Returns
///
/// The correlation coefficient, in `[-1, 1]`, or `None` if the slices have
/// different lengths, if the `xs` or the `ys` take fewer than two distinct
/// values, or if the sums of squares overflow or underflow (deviations from
/// the mean beyond about 10^±150).
pub fn correlation(xs: &[f64], ys: &[f64]) -> Option<f64> {
    if xs.len() != ys.len() || !has_distinct(xs, 2) || !has_distinct(ys, 2) {
        return None;
    }
    let (_, _, sxx, syy, sxy) = centred_sums(xs, ys);
    // One square root rounds less than two; fall back to two if the product
    // overflows or underflows.
    let product = sxx * syy;
    let denominator = if product.is_finite() && product > 0.0 {
        product.sqrt()
    } else {
        sxx.sqrt() * syy.sqrt()
    };
    let r = sxy / denominator;
    r.is_finite().then(|| r.clamp(-1.0, 1.0))
}

/// Solves the 3×3 normal equations of a least-squares fit with three
/// functions, given Meeus' sums (eq. 4.7):
/// `M = Σf₀², P = Σf₀f₁, Q = Σf₀f₂, R = Σf₁², S = Σf₁f₂, T = Σf₂²`,
/// `U = Σyf₀, V = Σyf₁, W = Σyf₂`.
fn solve_normal_equations(
    m: f64,
    p: f64,
    q: f64,
    r: f64,
    s: f64,
    t: f64,
    u: f64,
    v: f64,
    w: f64,
) -> Option<ThreeFunctionFit> {
    let d = m * r * t + 2.0 * p * q * s - m * s * s - r * q * q - t * p * p;
    // D is a Gram determinant, so 0 <= D <= M·R·T. A value that is zero to
    // within rounding means the functions are linearly dependent on the data.
    if !d.is_finite() || d <= 64.0 * f64::EPSILON * m * r * t {
        return None;
    }
    let fit = ThreeFunctionFit {
        a: (u * (r * t - s * s) + v * (q * s - p * t) + w * (p * s - q * r)) / d,
        b: (u * (s * q - p * t) + v * (m * t - q * q) + w * (p * q - m * s)) / d,
        c: (u * (p * s - r * q) + v * (p * q - m * s) + w * (m * r - p * p)) / d,
    };
    (fit.a.is_finite() && fit.b.is_finite() && fit.c.is_finite()).then_some(fit)
}

/// Fits a parabola `y = a x² + b x + c` to data by least squares (Meeus,
/// eqs. 4.5 and 4.6).
///
/// # Arguments
///
/// * `xs`: Values of the independent variable (taken as exact)
/// * `ys`: Measured values, same length as `xs`
///
/// # Returns
///
/// The fitted [`QuadraticFit`], or `None` if the slices have different
/// lengths, if there are fewer than three distinct `xs`, or if the result
/// is not finite.
pub fn quadratic_fit(xs: &[f64], ys: &[f64]) -> Option<QuadraticFit> {
    if xs.len() != ys.len() || !has_distinct(xs, 3) {
        return None;
    }

    // Fit y = a u² + b' u + c' with u = x − x̄, then expand back to x.
    let mx = mean(xs);
    let (mut s1, mut s2, mut s3, mut s4) = (0.0, 0.0, 0.0, 0.0);
    let (mut sy, mut suy, mut su2y) = (0.0, 0.0, 0.0);
    for (&x, &y) in xs.iter().zip(ys) {
        let u = x - mx;
        let u2 = u * u;
        s1 += u;
        s2 += u2;
        s3 += u2 * u;
        s4 += u2 * u2;
        sy += y;
        suy += u * y;
        su2y += u2 * y;
    }
    let n = xs.len() as f64;
    // With f₀ = u², f₁ = u, f₂ = 1, eq. 4.7 reduces to eqs. 4.5 and 4.6.
    let centred = solve_normal_equations(s4, s3, s2, s2, s1, n, su2y, suy, sy)?;

    let fit = QuadraticFit {
        a: centred.a,
        b: centred.b - 2.0 * centred.a * mx,
        c: centred.c - centred.b * mx + centred.a * mx * mx,
    };
    (fit.a.is_finite() && fit.b.is_finite() && fit.c.is_finite()).then_some(fit)
}

/// Fits `y = a f₀(x) + b f₁(x) + c f₂(x)` to data by least squares, for
/// three chosen functions of `x` (multiple linear regression; Meeus,
/// eq. 4.7).
///
/// For example `f₀ = sin x`, `f₁ = sin 2x`, `f₂ = sin 3x`, or
/// `f₀ = x²`, `f₁ = x`, `f₂ = 1` for a parabola (but prefer
/// [`quadratic_fit`] for that, which is better conditioned).
///
/// # Arguments
///
/// * `xs`: Values of the independent variable (taken as exact)
/// * `ys`: Measured values, same length as `xs`
/// * `f0`, `f1`, `f2`: The three functions
///
/// # Returns
///
/// The [`ThreeFunctionFit`] coefficients, or `None` if the slices have
/// different lengths, if there are fewer than three points, or if the three
/// functions are linearly dependent on the given `xs` (to within rounding),
/// so that the coefficients are not determined.
pub fn fit_three_functions<F0, F1, F2>(
    xs: &[f64],
    ys: &[f64],
    f0: F0,
    f1: F1,
    f2: F2,
) -> Option<ThreeFunctionFit>
where
    F0: Fn(f64) -> f64,
    F1: Fn(f64) -> f64,
    F2: Fn(f64) -> f64,
{
    if xs.len() != ys.len() || xs.len() < 3 {
        return None;
    }
    let (mut m, mut p, mut q, mut r, mut s, mut t) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut u, mut v, mut w) = (0.0, 0.0, 0.0);
    for (&x, &y) in xs.iter().zip(ys) {
        let (g0, g1, g2) = (f0(x), f1(x), f2(x));
        m += g0 * g0;
        p += g0 * g1;
        q += g0 * g2;
        r += g1 * g1;
        s += g1 * g2;
        t += g2 * g2;
        u += y * g0;
        v += y * g1;
        w += y * g2;
    }
    solve_normal_equations(m, p, q, r, s, t, u, v, w)
}

/// Fits `y = a f(x)` to data by least squares, for one chosen function of
/// `x` (Meeus, eq. 4.8).
///
/// # Arguments
///
/// * `xs`: Values of the independent variable (taken as exact)
/// * `ys`: Measured values, same length as `xs`
/// * `f`: The function
///
/// # Returns
///
/// The coefficient `a = Σ y f(x) / Σ f(x)²`, or `None` if the slices have
/// different lengths, or if `f` is zero at every `x` (including when there
/// are no points).
pub fn fit_one_function<F>(xs: &[f64], ys: &[f64], f: F) -> Option<f64>
where
    F: Fn(f64) -> f64,
{
    if xs.len() != ys.len() {
        return None;
    }
    let (mut syf, mut sff) = (0.0, 0.0);
    for (&x, &y) in xs.iter().zip(ys) {
        let fx = f(x);
        syf += y * fx;
        sff += fx * fx;
    }
    if sff == 0.0 {
        return None;
    }
    let a = syf / sff;
    a.is_finite().then_some(a)
}

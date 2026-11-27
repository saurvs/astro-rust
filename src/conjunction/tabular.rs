//! Quantities tabulated at equally spaced instants (crate-private helpers
//! for the conjunction and geometry modules).
//!
//! A quantity is given at 1 (constant), 3 or 5 instants and interpolated
//! with Meeus' formulas 3.3 and 3.8 through [`crate::interpol`]. Instants
//! are expressed as the interpolating factor `n`: 0 at the central tabulated
//! instant, ±1 one tabular interval later/earlier.

use crate::interpol;

/// Spacing of the search grid, in tabular intervals. A power of two, so the
/// tabulated instants are grid points exactly.
const GRID_STEP: f64 = 1.0 / 16.0;

/// A quantity tabulated at 1, 3 or 5 equally spaced instants
#[derive(Debug, Clone, Copy)]
pub(crate) struct Tabular {
    y: [f64; 5],
    len: usize,
}

impl Tabular {
    /// `None` unless `values` has 1, 3 or 5 elements.
    pub(crate) fn new(values: &[f64]) -> Option<Self> {
        let len = values.len();
        if !matches!(len, 1 | 3 | 5) {
            return None;
        }
        let mut y = [0.0; 5];
        y[..len].copy_from_slice(values);
        Some(Self { y, len })
    }

    /// Interpolated value at `n` (Meeus, eq. 3.3 or 3.8).
    pub(crate) fn value(&self, n: f64) -> f64 {
        let y = &self.y;
        match self.len {
            3 => interpol::three_values(y[0], y[1], y[2], n),
            5 => interpol::five_values(y[0], y[1], y[2], y[3], y[4], n),
            _ => y[0],
        }
    }

    /// Derivative with respect to `n` of the interpolating polynomial, per
    /// tabular interval.
    pub(crate) fn deriv(&self, n: f64) -> f64 {
        let y = &self.y;
        match self.len {
            3 => {
                // y = y2 + n/2 (a + b + n c)
                let a = y[1] - y[0];
                let b = y[2] - y[1];
                0.5 * (a + b) + n * (b - a)
            }
            5 => {
                // Second form of eq. 3.8:
                // y = y3 + n q + n² p + n³ r + n⁴ m
                let (a, b, c, d) = (y[1] - y[0], y[2] - y[1], y[3] - y[2], y[4] - y[3]);
                let (e, f, g) = (b - a, c - b, d - c);
                let (h, j) = (f - e, g - f);
                let k = j - h;
                let m = k / 24.0;
                let r = (h + j) / 12.0;
                let p = 0.5 * f - m;
                let q = 0.5 * (b + c) - r;
                q + n * (2.0 * p + n * (3.0 * r + n * 4.0 * m))
            }
            _ => 0.0,
        }
    }
}

/// The number of tabulated instants (3 or 5) shared by several series of
/// the given lengths. Every length must be either 1 (a constant, such as a
/// star) or that common number, and at least one series must vary.
pub(crate) fn common_len(lens: &[usize]) -> Option<usize> {
    let mut common = None;
    for &len in lens {
        match len {
            1 => {}
            3 | 5 => match common {
                None => common = Some(len),
                Some(c) if c == len => {}
                Some(_) => return None,
            },
            _ => return None,
        }
    }
    common
}

/// Half-width of the tabulated span, in tabular intervals.
pub(crate) fn half_span(len: usize) -> f64 {
    ((len - 1) / 2) as f64
}

/// The `i`-th tabulated element of `s`; a single element is constant.
pub(crate) fn nth<P>(s: &[P], i: usize) -> &P {
    if s.len() == 1 { &s[0] } else { &s[i] }
}

/// Makes a sequence of angles continuous (consecutive values differ by
/// less than π) by adding multiples of 2π; the first is reduced to
/// `[-π, π)`.
pub(crate) fn unwrap(angles: &mut [f64]) {
    let mut prev: Option<f64> = None;
    for a in angles.iter_mut() {
        *a = match prev {
            None => crate::angle::limit_to_pm_pi(*a),
            // Whole turns only, so values needing no change keep all digits
            Some(p) => *a + ((p - *a) / crate::angle::TWO_PI).round() * crate::angle::TWO_PI,
        };
        prev = Some(*a);
    }
}

/// Grid points from `-half` to `half` inclusive.
fn grid(half: f64) -> impl Iterator<Item = f64> {
    let steps = (2.0 * half / GRID_STEP).round() as i64;
    (0..=steps).map(move |i| -half + i as f64 * GRID_STEP)
}

/// Refines a root of `f` between `a` and `b`, where `f(a)` has sign
/// `fa < 0` or `fa > 0` and `f(b)` the other sign (or zero), by bisection
/// down to adjacent floating point numbers.
fn bisect<F: Fn(f64) -> f64>(f: &F, mut a: f64, fa: f64, mut b: f64) -> f64 {
    let a_negative = fa < 0.0;
    for _ in 0..200 {
        let m = 0.5 * (a + b);
        if m == a || m == b {
            break;
        }
        let fm = f(m);
        if fm == 0.0 {
            return m;
        }
        if (fm < 0.0) == a_negative {
            a = m;
        } else {
            b = m;
        }
    }
    0.5 * (a + b)
}

/// The root of `f` in `[-half, half]` that is nearest to `n = 0`, if any.
///
/// Sign changes are located on a grid of step 1/16 and refined by
/// bisection. A root where `f` touches zero without changing sign is found
/// only if it falls on a grid point.
pub(crate) fn root_nearest_centre<F: Fn(f64) -> f64>(f: F, half: f64) -> Option<f64> {
    let mut best: Option<f64> = None;
    let mut consider = |x: f64| {
        if x.is_finite() && best.is_none_or(|b| x.abs() < b.abs()) {
            best = Some(x);
        }
    };

    let mut prev: Option<(f64, f64)> = None;
    for x in grid(half) {
        let fx = f(x);
        if fx == 0.0 {
            consider(x);
        } else if let Some((x0, f0)) = prev {
            if f0.is_finite() && fx.is_finite() && f0 != 0.0 && (f0 < 0.0) != (fx < 0.0) {
                consider(bisect(&f, x0, f0, x));
            }
        }
        prev = Some((x, fx));
    }

    best
}

/// The local minima in `[-half, half]` of a function whose derivative is
/// `df`: the places where `df` goes from negative to zero or positive.
pub(crate) fn minima<F: Fn(f64) -> f64>(df: F, half: f64) -> Vec<f64> {
    let mut found = Vec::new();
    let mut prev: Option<(f64, f64)> = None;
    for x in grid(half) {
        let d = df(x);
        if let Some((x0, d0)) = prev {
            if d0 < 0.0 && d >= 0.0 {
                let root = if d == 0.0 { x } else { bisect(&df, x0, d0, x) };
                if root.is_finite() {
                    found.push(root);
                }
            }
        }
        prev = Some((x, d));
    }
    found
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derivatives_match_finite_differences() {
        let t3 = Tabular::new(&[1.0, 2.5, 0.5]).unwrap();
        let t5 = Tabular::new(&[36.125, 24.606, 15.486, 8.694, 4.133]).unwrap();
        let h = 1e-6;
        for i in -20..=20 {
            let n = i as f64 * 0.1;
            for t in [t3, t5] {
                let fd = (t.value(n + h) - t.value(n - h)) / (2.0 * h);
                assert!((t.deriv(n) - fd).abs() < 1e-6, "{n}");
            }
        }
        let c = Tabular::new(&[7.0]).unwrap();
        assert_eq!((c.value(0.3), c.deriv(0.3)), (7.0, 0.0));
        assert!(Tabular::new(&[1.0, 2.0]).is_none());
    }

    #[test]
    fn interpolation_reproduces_tabulated_values() {
        let v = [0.3, -1.2, 4.0, 2.2, 9.1];
        let t = Tabular::new(&v).unwrap();
        for (i, &y) in v.iter().enumerate() {
            assert!((t.value(i as f64 - 2.0) - y).abs() < 1e-12);
        }
    }

    // Meeus, example 3.d: the zero of y = 3 + 2x − 3x² nearest the centre
    // is −0.720759220056, where formula 3.6 converges slowly.
    #[test]
    fn root_of_strongly_curved_function() {
        let t = Tabular::new(&[-2.0, 3.0, 2.0]).unwrap();
        let n = root_nearest_centre(|n| t.value(n), 1.0).unwrap();
        assert!((n - -0.720759220056).abs() < 1e-12, "{n}");
    }

    #[test]
    fn common_lengths() {
        assert_eq!(common_len(&[3, 3]), Some(3));
        assert_eq!(common_len(&[1, 5, 5]), Some(5));
        assert_eq!(common_len(&[5, 1]), Some(5));
        assert_eq!(common_len(&[1, 1]), None);
        assert_eq!(common_len(&[3, 5]), None);
        assert_eq!(common_len(&[4, 4]), None);
        assert_eq!(common_len(&[0, 3]), None);
    }

    #[test]
    fn unwrap_across_zero_hours() {
        let mut a = [6.2, 0.05, 0.2];
        unwrap(&mut a);
        assert!((a[0] - (6.2 - crate::angle::TWO_PI)).abs() < 1e-15);
        assert_eq!(&a[1..], &[0.05, 0.2]);
        // Continuing past +π
        let mut a = [3.0, -3.1, 2.9 + 4.0 * crate::angle::TWO_PI];
        unwrap(&mut a);
        assert_eq!(a[0], 3.0);
        assert!((a[1] - (crate::angle::TWO_PI - 3.1)).abs() < 1e-15);
        assert!((a[2] - 2.9).abs() < 1e-13);
    }

    #[test]
    fn minima_and_roots() {
        // (n − 0.3)² has its minimum at 0.3; sin(πn) has roots at −1, 0, 1
        let m = minima(|n| 2.0 * (n - 0.3), 1.0);
        assert_eq!(m.len(), 1);
        assert!((m[0] - 0.3).abs() < 1e-15);
        let r = root_nearest_centre(|n| (std::f64::consts::PI * (n - 0.01)).sin(), 2.0).unwrap();
        assert!((r - 0.01).abs() < 1e-15);
        assert!(root_nearest_centre(|n| n * n + 1.0, 1.0).is_none());
        assert!(minima(|n| n.exp(), 1.0).is_empty());
    }
}

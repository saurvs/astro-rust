/*
Copyright (c) 2015, 2016 Saurav Sachidanand

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
//! Parabolic orbits
//!
//! Motion of a body (usually a comet) in an unperturbed parabolic orbit
//! around the Sun (Meeus, *Astronomical Algorithms*, 2nd ed., chapters 34
//! and 39). Times are Julian Ephemeris Days (TT), angles are in radians and
//! distances in AU.
//!
//! The position follows from Barker's equation `s³ + 3s = W`, where
//! `s = tan(v/2)` and `W` is proportional to the time since perihelion
//! ([`barker_w`]). [`solve_barker`] solves it exactly; the three methods of
//! Meeus are also available ([`solve_barker_iter`],
//! [`solve_barker_bauschinger`] and [`solve_barker_cubic`]). For the
//! geocentric position of a comet in a parabolic orbit, use
//! [`orbit::Motion::Parabolic`] with [`orbit::geocent_eq_pos`].
//!
//! # Example
//!
//! Comet Stonehouse (C/1998 H1) on 1998 August 5.0 TD (Meeus,
//! example 34.a):
//!
//! ```
//! use astro::orbit::parabolic;
//!
//! let T = 2450917.9358; // time of perihelion, 1998 April 14.4358 TD
//! let t = 2451030.5;    // 1998 August 5.0 TD
//! let q = 1.487469;     // perihelion distance, AU
//!
//! let (v, r) = parabolic::true_anom_and_rad_vec(t, T, q);
//! assert!((v.to_degrees() - 66.78862).abs() < 1e-5);
//! assert!((r - 2.133911).abs() < 1e-6);
//! ```

use crate::consts;
use crate::orbit;
use std::f64::consts::PI;

/// `3k/√2 = 0.03649116245...`, where `k` is the Gaussian gravitational
/// constant (Meeus, eq. 34.1). Computed from `k` rather than rounded, so
/// that parabolic, near-parabolic and elliptic motions agree exactly.
const W_COEFF: f64 = 3.0 * consts::GAUSS_GRAV / std::f64::consts::SQRT_2;

/// Computes the parameter `W` of Barker's equation for a parabolic orbit
///
/// Uses `W = 3k / √2 × (t - T) / q^(3/2)`, where `k` is the Gaussian
/// gravitational constant (Meeus, eq. 34.1).
///
/// # Returns
///
/// * `W`: The right-hand side of Barker's equation `s³ + 3s = W`; negative
///   before the perihelion
///
/// # Arguments
///
/// * `t`: Time of interest, in Julian (Ephemeris) day
/// * `T`: Time of passage through the perihelion, in Julian (Ephemeris) day
/// * `q`: Perihelion distance *| in AU*
#[inline]
pub fn barker_w(t: f64, T: f64, q: f64) -> f64 {
    W_COEFF * (t - T) / (q * q.sqrt())
}

/// Solves Barker's equation `s³ + 3s = W` for `s = tan(v/2)`
///
/// The cubic has a single real root for every `W`. This function computes
/// it as `s = 2 sinh(asinh(W/2)/3)`, which is the closed form of Meeus,
/// eq. 34.6 written so that it is accurate for any `W`, positive or
/// negative, small or large.
///
/// # Returns
///
/// * `s`: `tan(v/2)`, where `v` is the true anomaly; it has the sign of `W`
///
/// # Arguments
///
/// * `W`: Parameter of Barker's equation, from [`barker_w`]
#[inline]
pub fn solve_barker(W: f64) -> f64 {
    2.0 * ((W / 2.0).asinh() / 3.0).sinh()
}

/// Solves Barker's equation by iteration
///
/// This is the first method of Meeus, chapter 34: starting from `s = 0`,
/// `s` is replaced by `(2s³ + W) / (3(s² + 1))` (eq. 34.4, which is Newton's
/// method for the cubic) until two successive values differ by at most
/// `accuracy`. The iteration converges for any `W`.
///
/// # Returns
///
/// * `Some(s)`: `tan(v/2)`, where `v` is the true anomaly
/// * `None` if `W` is not finite, or if `accuracy` has not been reached
///   after `max_iter` iterations
///
/// # Arguments
///
/// * `W`       : Parameter of Barker's equation, from [`barker_w`]
/// * `accuracy`: Stop when two successive values differ by at most this
/// * `max_iter`: Maximum number of iterations
pub fn solve_barker_iter(W: f64, accuracy: f64, max_iter: u32) -> Option<f64> {
    if !W.is_finite() {
        return None;
    }

    let mut s = 0.0;
    for _ in 0..max_iter {
        let next = (2.0 * s * s * s + W) / (3.0 * (s * s + 1.0));
        let diff = (next - s).abs();
        s = next;
        if diff <= accuracy {
            return Some(s);
        }
    }

    None
}

/// Solves Barker's equation with Bauschinger's trigonometric formulae
///
/// This is the second method of Meeus, chapter 34 (eq. 34.5):
/// `tan(β) = 2/W`, `tan(γ) = ∛tan(β/2)`, `s = 2/tan(2γ)`, computed with
/// `|W|` and the sign of `W` restored at the end, and with `s = 0` for
/// `W = 0`.
///
/// # Returns
///
/// * `s`: `tan(v/2)`, where `v` is the true anomaly
///
/// # Arguments
///
/// * `W`: Parameter of Barker's equation, from [`barker_w`]
pub fn solve_barker_bauschinger(W: f64) -> f64 {
    if W == 0.0 {
        return 0.0;
    }

    let beta = (2.0 / W.abs()).atan();
    let gamma = (beta / 2.0).tan().cbrt().atan();
    let s = 2.0 / (2.0 * gamma).tan();

    s.copysign(W)
}

/// Solves Barker's equation with Meeus's algebraic formulae
///
/// This is the third method of Meeus, chapter 34 (eq. 34.6):
/// `G = W/2`, `Y = ∛(G + √(G² + 1))`, `s = Y - 1/Y`. It is evaluated for
/// `|W|`, and the sign restored at the end, to avoid the cancellation in
/// `G + √(G² + 1)` when `W` is large and negative. [`solve_barker`] gives
/// the same root with better accuracy for very small `|W|`.
///
/// # Returns
///
/// * `s`: `tan(v/2)`, where `v` is the true anomaly
///
/// # Arguments
///
/// * `W`: Parameter of Barker's equation, from [`barker_w`]
pub fn solve_barker_cubic(W: f64) -> f64 {
    let G = W.abs() / 2.0;
    let Y = (G + (G * G + 1.0).sqrt()).cbrt();

    (Y - 1.0 / Y).copysign(W)
}

/// Computes the true anomaly and radius vector of a body in a parabolic
/// orbit at a given time
///
/// Meeus, eqs. 34.1 to 34.3: Barker's equation is solved with
/// [`solve_barker`], then `v = 2 atan(s)` and `r = q (1 + s²)`.
///
/// # Returns
///
/// `(true_anom, rad_vec)`
///
/// * `true_anom`: True anomaly of the body at time `t` *| in radians*, in
///   `(-π, π)`: negative before the perihelion, positive after it
/// * `rad_vec`  : Radius vector of the body at time `t` *| in AU*
///
/// # Arguments
///
/// * `t`: Time of interest, in Julian (Ephemeris) day
/// * `T`: Time of passage through the perihelion, in Julian (Ephemeris) day
/// * `q`: Perihelion distance *| in AU*
pub fn true_anom_and_rad_vec(t: f64, T: f64, q: f64) -> (f64, f64) {
    let s = solve_barker(barker_w(t, T, q));
    let v = 2.0 * s.atan();
    let r = q * (1.0 + s * s);

    (v, r)
}

/// Computes the time of passage of a body through a node of a parabolic
/// orbit and its radius vector at that time
///
/// Meeus, chapter 39: with `v = -ω` (ascending node) or `v = π - ω`
/// (descending node) and `s = tan(v/2)`, the time is
/// `T + 27.403895 (s³ + 3s) q^(3/2)` days (the constant is `√2/(3k)`, `k`
/// being the Gaussian gravitational constant), and `r = q (1 + s²)`.
///
/// The node is on the plane of reference of the elements (e.g. the
/// ecliptic of J2000.0 for elements referred to that equinox).
///
/// If the node is at `v = ±π` (`ω = π` for the ascending node, `ω = 0`
/// for the descending one) the body never reaches it; the results are then
/// very large or infinite.
///
/// # Returns
///
/// `(time_of_pass, rad_vec)`
///
/// * `time_of_pass`: Time of passage through the node, in Julian
///   (Ephemeris) day
/// * `rad_vec`     : Radius vector of the body at the time of passage
///   *| in AU*
///
/// # Arguments
///
/// * `w`   : Argument of the perihelion *| in radians*
/// * `q`   : Perihelion distance *| in AU*
/// * `T`   : Time of passage in perihelion, in Julian (Ephemeris) day
/// * `node`: `Ascend` or `Descend` node
///
/// # Example
///
/// Comet Helin-Roman at its descending node (Meeus, example 39.b):
///
/// ```
/// use astro::orbit::{self, parabolic};
///
/// let T = 2447758.791; // 1989 Aug. 20.2910 TD
/// let (t, r) = parabolic::passage_through_node(
///     154.9103_f64.to_radians(),
///     1.324502,
///     T,
///     &orbit::Node::Descend,
/// );
/// assert!((t - (T + 28.3454)).abs() < 1e-4);
/// assert!((r - 1.3901).abs() < 1e-4);
/// ```
#[inline]
pub fn passage_through_node(w: f64, q: f64, T: f64, node: &orbit::Node) -> (f64, f64) {
    match *node {
        orbit::Node::Ascend => pass_through_node(-w, q, T),
        orbit::Node::Descend => pass_through_node(PI - w, q, T),
    }
}

fn pass_through_node(v: f64, q: f64, T: f64) -> (f64, f64) {
    let s = (v / 2.0).tan();
    // 1 / W_COEFF = √2 / (3k) = 27.403895...
    let T_node = T + q.powf(1.5) * (s * (s * s + 3.0)) / W_COEFF;
    let rad_vec = q * (1.0 + s * s);

    (T_node, rad_vec)
}

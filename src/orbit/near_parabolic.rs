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
//! Near-parabolic orbits
//!
//! Motion of a body in an orbit with an eccentricity close to 1, either
//! slightly elliptic or slightly hyperbolic, by the method of W. Landgraf
//! (Meeus, *Astronomical Algorithms*, 2nd ed., chapter 35). Times are
//! Julian Ephemeris Days (TT), angles are in radians and distances in AU.
//!
//! The method solves, by iteration, a series in `s = tan(v/2)` that
//! generalises Barker's equation of the parabolic orbit. The series
//! converges only while the body is not too far from the perihelion: the
//! farther `e` is from 1, the shorter the time span over which it works
//! (for `q = 0.1` AU, about 20 days for `e = 0.9`, 400 days for
//! `e = 0.987`, more than 5000 days for `e = 0.999`). Outside that span the
//! functions return an error; for an elliptic orbit, use Kepler's equation
//! instead (see [`orbit::elliptic`]).
//!
//! # Example
//!
//! The first case of Meeus's table in chapter 35, a parabolic orbit:
//!
//! ```
//! use astro::orbit::near_parabolic;
//!
//! let (v, r) = near_parabolic::true_anom_and_rad_vec(138.4783, 0.0, 1.0, 0.921326, 1e-10)
//!     .unwrap();
//! assert!((v.to_degrees() - 102.74426).abs() < 1e-5);
//! assert!((r - 2.364192).abs() < 1e-6);
//! ```

use crate::consts;
use crate::orbit::{self, parabolic};

/// The message of the error returned when Landgraf's method does not
/// converge
pub const NO_CONVERGENCE: &str = "No convergence at orbit::near_parabolic::true_anom_and_rad_vec()";

/// The message of the error returned for invalid arguments
pub const INVALID_INPUT: &str = "Invalid input at orbit::near_parabolic::true_anom_and_rad_vec()";

/// Largest number of terms of the series, and of outer iterations
const MAX_TERMS: u32 = 50;

/// A term of the series larger than this means that it diverges
const MAX_TERM: f64 = 10_000.0;

/// Largest number of Newton steps in the inner iteration
const MAX_NEWTON: u32 = 100;

/// Computes the true anomaly and radius vector of a body in a near-parabolic
/// orbit at a given time
///
/// Implements Landgraf's method (Meeus, chapter 35, eq. 35.1 and the
/// program that follows it). It is valid for elliptic orbits with `e`
/// slightly below 1, for hyperbolic orbits with `e` slightly above 1, and
/// for the parabola (`e = 1`, where it reduces to Barker's equation).
///
/// # Returns
///
/// * `Ok((true_anom, rad_vec))`
///   * `true_anom`: True anomaly of the body at time `t` *| in radians*, in
///     `[0, 2π)` (an instant before the perihelion gives a value above π)
///   * `rad_vec`  : Radius vector of the body at time `t` *| in AU*
/// * `Err(NO_CONVERGENCE)` if the series of eq. 35.1 does not converge
///   (the body is too far from the perihelion for its eccentricity), or the
///   iteration does not reach `accuracy`
/// * `Err(INVALID_INPUT)` if `q` is not positive, `ecc` is negative, or an
///   argument is not finite
///
/// # Arguments
///
/// * `t`       : Time of interest, in Julian (Ephemeris) day
/// * `T`       : Time of passage through the perihelion, in Julian
///   (Ephemeris) day
/// * `ecc`     : Eccentricity of the near-parabolic orbit
/// * `q`       : Perihelion distance *| in AU*
/// * `accuracy`: Desired accuracy for `s = tan(v/2)` (and the size of the
///   last term kept in the series), e.g. `1e-9` as in the book. A smaller
///   value needs more terms of the series and so slightly narrows the time
///   span over which the method converges; below about `1e-15` it may
///   prevent convergence altogether.
pub fn true_anom_and_rad_vec(
    t: f64,
    T: f64,
    ecc: f64,
    q: f64,
    accuracy: f64,
) -> Result<(f64, f64), &'static str> {
    let finite = t.is_finite() && T.is_finite() && ecc.is_finite() && q.is_finite();
    if !finite || q <= 0.0 || ecc < 0.0 || accuracy.is_nan() {
        return Err(INVALID_INPUT);
    }

    let days_frm_perih = t - T;
    if days_frm_perih == 0.0 {
        return Ok((0.0, q));
    }

    let s = solve_s(days_frm_perih, ecc, q, accuracy)?;

    let v = orbit::limit_to_two_pi(2.0 * s.atan());
    let r = q * (1.0 + ecc) / (1.0 + ecc * v.cos());

    Ok((v, r))
}

/// Solves eq. 35.1 for `s = tan(v/2)`
fn solve_s(days_frm_perih: f64, ecc: f64, q: f64, accuracy: f64) -> Result<f64, &'static str> {
    let k = consts::GAUSS_GRAV;

    let q1 = k * ((1.0 + ecc) / q).sqrt() / (2.0 * q);
    let g = (1.0 - ecc) / (1.0 + ecc);
    let q2 = q1 * days_frm_perih;

    // Starting value: the root for a parabola, for which s³/3 + s = q2.
    // (Landgraf's program evaluates the same root with Bauschinger's
    // trigonometric formulae, which lose some precision.)
    let mut s = parabolic::solve_barker(3.0 * q2);
    if ecc == 1.0 {
        return Ok(s);
    }

    let mut outer = 0;
    loop {
        let s0 = s;

        // the right-hand side of s³/3 + s = q3, summed term by term
        let y = s * s;
        let mut g1 = -y * s;
        let mut q3 = q2 + g * s * y * 2.0 / 3.0;
        let mut z = 1.0;
        loop {
            z += 1.0;
            g1 = -g1 * g * y;
            let z1 = (z - (z + 1.0) * g) / (2.0 * z + 1.0);
            let f = z1 * g1;
            q3 += f;
            if z > MAX_TERMS as f64 || f.is_nan() || f.abs() > MAX_TERM {
                return Err(NO_CONVERGENCE);
            }
            if f.abs() <= accuracy {
                break;
            }
        }

        outer += 1;
        if outer > MAX_TERMS {
            return Err(NO_CONVERGENCE);
        }

        // Newton's method for s³/3 + s = q3
        let mut newton = 0;
        loop {
            let s1 = s;
            s = (s * s * s * 2.0 / 3.0 + q3) / (s * s + 1.0);
            if (s - s1).abs() <= accuracy {
                break;
            }
            newton += 1;
            if newton > MAX_NEWTON {
                return Err(NO_CONVERGENCE);
            }
        }

        if (s - s0).abs() <= accuracy {
            return Ok(s);
        }
    }
}

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

//! Diameters of asteroids
//!
//! True diameter of an asteroid from its absolute magnitude and albedo,
//! and its apparent diameter seen from the Earth (Meeus, *Astronomical
//! Algorithms*, 2nd ed., chapter 55). Because asteroids are often
//! irregular, these give only a rough size.
//!
//! ```
//! use astro::asteroid;
//!
//! // An asteroid with H = 3.34 and albedo 0.09 (about the size of Ceres)
//! let d = asteroid::diameter(3.34, 0.09); // meters
//! assert!((d / 1000.0 - 944.0).abs() < 1.0);
//!
//! // Seen from 1.6 AU it spans about 0.8″
//! let app = asteroid::apparent_diameter(d / 1000.0, 1.6); // radians
//! assert!((app.to_degrees() * 3600.0 - 0.813).abs() < 0.001);
//! ```

use crate::consts;

/// Computes the diameter of an asteroid
///
/// `log10 d = 3.12 - H/5 - 0.5 log10 A`, with `d` in kilometers (Meeus,
/// chapter 55). Typical albedos: 0.04 for many dark asteroids, 0.10 for
/// Ceres, 0.38 for Vesta. A non-positive albedo gives NaN or +∞.
///
/// # Returns
///
/// * `diameter`: Diameter of the asteroid *| in meters*
///
/// # Arguments
///
/// * `abs_mag`: Absolute magnitude `H` of the asteroid
/// * `albedo` : Albedo `A` (reflective power) of the asteroid, between 0
///   and 1
#[inline]
pub fn diameter(abs_mag: f64, albedo: f64) -> f64 {
    // Bug fix: the albedo term used the natural-log coefficient 0.217147
    // with a base-10 logarithm; with log10 the coefficient is 0.5. For an
    // albedo of 0.04 the diameter came out 2.5 times too small.
    1000.0 * 10_f64.powf(3.12 - abs_mag / 5.0 - 0.5 * albedo.log10())
}

/// Computes the apparent diameter of an asteroid
///
/// The angle subtended by a body of diameter `d` at distance `Δ`, which
/// is `0.0013788″ × d / Δ` for `d` in km and `Δ` in AU (Meeus, chapter 55).
///
/// # Returns
///
/// * `apparent_diameter`: Apparent diameter of the asteroid
///   *| in radians*
///
/// # Arguments
///
/// * `true_diameter`      : True diameter of the asteroid *| in kilometers*
///   (note: [`diameter`] returns meters)
/// * `asteroid_earth_dist`: Asteroid-Earth distance *| in AU*
#[inline]
pub fn apparent_diameter(true_diameter: f64, asteroid_earth_dist: f64) -> f64 {
    // Bug fix: this returned 1.3788 d / Δ, which is the angle in
    // milli-arcseconds, documented as "meters". It now returns radians.
    true_diameter * 1000.0 / (asteroid_earth_dist * consts::ASTRONOMICAL_UNIT)
}

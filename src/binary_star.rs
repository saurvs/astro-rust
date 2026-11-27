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

//! Binary stars
//!
//! Apparent position of the companion of a visual binary star from the
//! orbital elements of its true orbit, and the eccentricity of its
//! apparent orbit (Meeus, *Astronomical Algorithms*, 2nd ed., chapter 57).
//!
//! The elements follow the usual conventions for visual binaries: the
//! period in mean solar years, the epoch of periastron as a year with
//! decimals, the semimajor axis as an angle, the inclination of the true
//! orbit to the plane of the sky (0 to 90° for direct apparent motion, 90°
//! to 180° for retrograde), the position angle of the ascending node, and
//! the longitude of periastron measured in the plane of the true orbit
//! from the ascending node in the direction of motion. Position angles are
//! counted from north through east (0 north, π/2 east).
//!
//! ```
//! use astro::binary_star::{self, OrbitalElements};
//!
//! // η Coronae Borealis (Meeus, example 57.a), semimajor axis in arcseconds
//! let eta_crb = OrbitalElements {
//!     period: 41.623,
//!     time_of_periastron: 1934.008,
//!     eccentricity: 0.2763,
//!     semimajor_axis: 0.907,
//!     inclination: 59.025_f64.to_radians(),
//!     asc_node_position_angle: 23.717_f64.to_radians(),
//!     periastron_longitude: 219.907_f64.to_radians(),
//! };
//!
//! let pos = binary_star::apparent_position(&eta_crb, 1980.0);
//! assert!((pos.position_angle.to_degrees() - 318.4).abs() < 0.05);
//! assert!((pos.separation - 0.411).abs() < 0.0005); // arcseconds
//!
//! let e_app = eta_crb.apparent_eccentricity();
//! assert!((e_app - 0.860).abs() < 0.0005);
//! ```

use crate::angle;
use crate::orbit;

/// Computes the mean annual motion of the companion star
///
/// `n = 2π / P` (Meeus, chapter 57).
///
/// # Returns
///
/// * Mean annual motion *| in radians per year*
///
/// # Arguments
///
/// * `P`: Period of revolution of the binary star *| in mean solar years*
#[inline]
pub fn mn_ann_motion_of_compan(P: f64) -> f64 {
    angle::TWO_PI / P
}

/// Computes the mean anomaly of the companion star
///
/// `M = n (t - T)`, not reduced to `[0, 2π)` (Meeus, chapter 57).
///
/// # Returns
///
/// * Mean anomaly *| in radians*
///
/// # Arguments
///
/// * `n`: Mean annual motion of the companion star *| in radians per
///   year*
/// * `t`: Current time, as a year with decimals (eg: 1980.0)
/// * `T`: Time of periastron passage, as a year with decimals
///   (eg: 1945.62)
#[inline]
pub fn mn_anom_of_compan(n: f64, t: f64, T: f64) -> f64 {
    n * (t - T)
}

/// Computes the radius vector of the companion star in its true orbit
///
/// `r = a (1 - e cos E)` (Meeus, chapter 57).
///
/// # Returns
///
/// * Radius vector, in the same unit as `a`
///
/// # Arguments
///
/// * `a`       : Semimajor axis of the true orbit (an angle, eg: in
///   arcseconds or radians)
/// * `e`       : Eccentricity of the true orbit
/// * `ecc_anom`: Eccentric anomaly of the companion star *| in radians*
#[inline]
pub fn rad_vec(a: f64, e: f64, ecc_anom: f64) -> f64 {
    a * (1.0 - e * ecc_anom.cos())
}

/// Computes the true anomaly of the companion star
///
/// `tan(v/2) = sqrt((1 + e) / (1 - e)) tan(E/2)` (Meeus, chapter 57).
///
/// # Returns
///
/// * True anomaly *| in radians*, in `(-π, π]`
///
/// # Arguments
///
/// * `e`       : Eccentricity of the true orbit
/// * `ecc_anom`: Eccentric anomaly of the companion star *| in radians*
#[inline]
pub fn true_anom(e: f64, ecc_anom: f64) -> f64 {
    2.0 * (((1.0 + e) / (1.0 - e)).sqrt() * (ecc_anom / 2.0).tan()).atan()
}

/// Computes the apparent position angle of the companion star
///
/// Uses `tan(θ - Ω) = sin(v + ω) cos i / cos(v + ω)` with the quadrant
/// taken from the signs of numerator and denominator (Meeus, eq. 57.1).
///
/// # Returns
///
/// * Position angle `θ`, from north through east *| in radians*, in
///   `[0, 2π)`. Between 0 and π the companion "follows" the primary in the
///   diurnal motion; between π and 2π it "precedes" it.
///
/// # Arguments
///
/// * `asc_node_coords`: Position angle of the ascending node *| in
///   radians*
/// * `true_anom`      : True anomaly of the companion star *| in radians*
/// * `w`              : Longitude of periastron *| in radians*
/// * `i`              : Inclination of the true orbit to the plane of the
///   sky (the plane at right angles to the line of sight) *| in radians*
pub fn apprnt_coords_angl(asc_node_coords: f64, true_anom: f64, w: f64, i: f64) -> f64 {
    let x = ((true_anom + w).sin() * i.cos()).atan2((true_anom + w).cos());

    angle::limit_to_two_PI(x + asc_node_coords)
}

/// Computes the apparent angular separation of the two components
///
/// Uses Greaney's form `ρ = r sqrt(sin²(v + ω) cos²i + cos²(v + ω))`, which
/// avoids the division by `cos(θ - Ω)` (Meeus, chapter 57).
///
/// # Returns
///
/// * Angular separation `ρ`, in the same unit as `rad_vec`
///
/// # Arguments
///
/// * `rad_vec`  : Radius vector of the companion star (an angle, eg: in
///   arcseconds or radians)
/// * `true_anom`: True anomaly of the companion star *| in radians*
/// * `w`        : Longitude of periastron *| in radians*
/// * `i`        : Inclination of the true orbit to the plane of the sky
///   *| in radians*
pub fn anglr_sepr(rad_vec: f64, true_anom: f64, w: f64, i: f64) -> f64 {
    rad_vec * (((true_anom + w).sin() * i.cos()).powi(2) + (true_anom + w).cos().powi(2)).sqrt()
}

/// Computes the eccentricity of the apparent orbit of a binary star
///
/// The projection of the true orbit on the sky is an ellipse whose
/// eccentricity `e'` depends only on `e`, `ω` and `i` (Meeus, chapter 57,
/// after Meeus 1979). It may be larger or smaller than `e`; it equals `e`
/// for a face-on orbit (`i = 0`) and 1 for an edge-on orbit (`i = π/2`).
///
/// # Arguments
///
/// * `e`: Eccentricity of the true orbit
/// * `w`: Longitude of periastron *| in radians*
/// * `i`: Inclination of the true orbit to the plane of the sky
///   *| in radians*
pub fn ecc_of_apprnt_orb(e: f64, w: f64, i: f64) -> f64 {
    let i_cos = i.cos();
    let e2 = e * e;

    let a = (1.0 - e2 * w.cos().powi(2)) * i_cos * i_cos;
    let b = e2 * w.sin() * w.cos() * i_cos;
    // Bug fix: C uses sin²ω (it used cos²ω, which made e' = 0.8578 instead
    // of 0.860 in Meeus' example 57.b).
    let c = 1.0 - e2 * w.sin().powi(2);
    let sqrt_d = ((a - c) * (a - c) + 4.0 * b * b).sqrt();

    ((2.0 * sqrt_d) / (a + c + sqrt_d)).sqrt()
}

/// Orbital elements of the true orbit of a visual binary star
///
/// See the [module documentation](self) for the conventions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrbitalElements {
    /// Period of revolution `P` *| in mean solar years*
    pub period: f64,
    /// Time of periastron passage `T`, as a year with decimals
    /// (eg: 1934.008)
    pub time_of_periastron: f64,
    /// Eccentricity `e` of the true orbit, in `[0, 1)`
    pub eccentricity: f64,
    /// Semimajor axis `a`, as an angle. Any unit may be used (the
    /// separation is returned in the same unit); radians keep with the
    /// rest of the crate, arcseconds with catalogues.
    pub semimajor_axis: f64,
    /// Inclination `i` of the true orbit to the plane of the sky
    /// *| in radians*, in `[0, π]` (more than π/2 for retrograde motion)
    pub inclination: f64,
    /// Position angle `Ω` of the ascending node *| in radians*
    pub asc_node_position_angle: f64,
    /// Longitude of periastron `ω`, measured in the plane of the true
    /// orbit from the ascending node in the direction of motion
    /// *| in radians*
    pub periastron_longitude: f64,
}

/// Apparent position of the companion relative to the primary star
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ApparentPosition {
    /// Position angle `θ`, from north through east *| in radians*, in
    /// `[0, 2π)`
    pub position_angle: f64,
    /// Angular separation `ρ`, in the unit of the semimajor axis
    pub separation: f64,
}

impl OrbitalElements {
    /// Computes the apparent position of the companion at time `t`, a year
    /// with decimals. See [`apparent_position`].
    pub fn apparent_position(&self, t: f64) -> ApparentPosition {
        apparent_position(self, t)
    }

    /// Computes the eccentricity of the apparent orbit. See
    /// [`ecc_of_apprnt_orb`].
    pub fn apparent_eccentricity(&self) -> f64 {
        ecc_of_apprnt_orb(
            self.eccentricity,
            self.periastron_longitude,
            self.inclination,
        )
    }
}

/// Computes the apparent position angle and separation of a binary star
/// at a given time
///
/// Chains the steps of Meeus' chapter 57: mean motion and mean anomaly,
/// Kepler's equation (with [`orbit::elliptic::ecc_anom_checked`], to full
/// double precision for any eccentricity below 1), radius vector, true
/// anomaly, then position angle (eq. 57.1) and separation. The position
/// angle refers to the same equinox as the element `Ω`.
///
/// For an eccentricity outside `[0, 1)` or non-finite input the result
/// contains NaN.
///
/// # Arguments
///
/// * `elements`: Orbital elements of the true orbit
/// * `t`       : Time, as a year with decimals (eg: 1980.0)
pub fn apparent_position(elements: &OrbitalElements, t: f64) -> ApparentPosition {
    let e = elements.eccentricity;
    let n = mn_ann_motion_of_compan(elements.period);
    let M = mn_anom_of_compan(n, t, elements.time_of_periastron);
    // The orbit module's solver converges for any 0 <= e < 1 and returns
    // None for an eccentricity outside that range or non-finite input.
    let E = orbit::elliptic::ecc_anom_checked(M, e).unwrap_or(f64::NAN);
    let r = rad_vec(elements.semimajor_axis, e, E);
    let v = true_anom(e, E);
    let w = elements.periastron_longitude;
    let i = elements.inclination;

    ApparentPosition {
        position_angle: apprnt_coords_angl(elements.asc_node_position_angle, v, w, i),
        separation: anglr_sepr(r, v, w, i),
    }
}

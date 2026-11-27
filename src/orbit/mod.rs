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
//! Elliptic, parabolic and near-parabolic orbits
//!
//! Positions of minor planets and comets from their orbital elements
//! (Meeus, *Astronomical Algorithms*, 2nd ed., chapters 30, 33 to 35 and 39).
//!
//! The submodules solve for the true anomaly and the radius vector in each
//! kind of orbit: [`elliptic`] (Kepler's equation, `e < 1`), [`parabolic`]
//! (Barker's equation, `e = 1`) and [`near_parabolic`] (Landgraf's method,
//! `e` close to 1 on either side). This module turns those into
//! heliocentric and geocentric positions:
//!
//! * [`Elements`] and [`Motion`] describe an orbit;
//! * [`geocent_eq_pos`] gives the geocentric right ascension and declination
//!   from the Sun's geocentric rectangular equatorial coordinates (Meeus's
//!   second method, eqs. 33.7 to 33.12), with the correction for light-time;
//! * [`geocent_ecl_pos`] gives the geocentric ecliptic longitude and latitude
//!   from the Earth's heliocentric ecliptic coordinates (as in Meeus's first
//!   method, eqs. 33.1 to 33.5, with the body's heliocentric position taken
//!   from its elements), also with the correction for light-time;
//! * lower-level pieces of both methods are public too
//!   ([`eq_orbit_consts`], [`gauss_vectors`], [`heliocent_ecl_rect`], ...).
//!
//! Conventions: angles in radians, distances in AU, times in Julian
//! Ephemeris Days (TT). The angular elements (`incl`, `arg_perih`,
//! `long_asc_node`), the Sun's or the Earth's coordinates and the obliquity
//! must all refer to the same equinox (for example J2000.0); the results
//! then refer to that equinox too.
//!
//! # Example
//!
//! The astrometric position of comet Encke on 1990 October 6.0 TD, using
//! the Sun's rectangular coordinates for the equinox J2000.0 (Meeus,
//! example 33.b):
//!
//! ```
//! use astro::orbit::{self, Elements, Motion};
//!
//! let encke = Elements {
//!     motion: Motion::elliptic_frm_perih_time(
//!         2.2091404,    // semimajor axis, AU
//!         0.8502196,    // eccentricity
//!         2448192.5 + 0.54502, // time of perihelion, 1990 Oct. 28.54502 TD
//!     ),
//!     incl: 11.94524_f64.to_radians(),
//!     arg_perih: 186.23352_f64.to_radians(),
//!     long_asc_node: 334.75006_f64.to_radians(),
//! };
//!
//! let sun = (-0.9756732, -0.2003254, -0.0868566);
//! let pos = orbit::geocent_eq_pos(&encke, 2448170.5, sun, orbit::OBLQ_J2000).unwrap();
//!
//! assert!((pos.asc.to_degrees() - 158.558965).abs() < 1e-5);
//! assert!((pos.dec.to_degrees() - 19.158496).abs() < 1e-5);
//! assert!((pos.earth_dist - 0.8242811).abs() < 1e-6);
//! assert!((pos.elong.to_degrees() - 40.51).abs() < 0.01);
//! ```

pub mod elliptic;
pub mod near_parabolic;
pub mod parabolic;

use std::f64::consts::{PI, TAU};

/// Reduces an angle to `[0, 2π)`
#[inline]
pub(crate) fn limit_to_two_pi(angl: f64) -> f64 {
    let r = angl.rem_euclid(TAU);
    if r >= TAU { 0.0 } else { r }
}

/// Represents an orbital node
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Node {
    /// Ascending node
    Ascend,
    /// Descending node
    Descend,
}

/// Mean obliquity of the ecliptic at the epoch J2000.0, 23°26′21.448″
/// *| in radians*
///
/// Use it with orbital elements and solar coordinates referred to the
/// standard equinox J2000.0 (Meeus, chapter 33).
pub const OBLQ_J2000: f64 = (23.0 + 26.0 / 60.0 + 21.448 / 3600.0) * (PI / 180.0);

/// Light-time for a distance of one AU *| in days* (Meeus, eq. 33.3)
const LIGHT_TIME_PER_AU: f64 = 0.0057755183;

/// Computes the light-time for a given distance
///
/// Uses Meeus, eq. 33.3: the light-time is 0.0057755183 days per AU.
///
/// # Returns
///
/// * `light_time`: Time taken by light to cover `dist` *| in days*
///
/// # Arguments
///
/// * `dist`: Distance *| in AU*
#[inline]
pub fn light_time(dist: f64) -> f64 {
    LIGHT_TIME_PER_AU * dist
}

/// Computes the radius vector of a body in any conic orbit from its true
/// anomaly
///
/// Uses `r = q (1 + e) / (1 + e cos(v))` (Meeus, eq. 30.4), which holds for
/// elliptic, parabolic and hyperbolic orbits.
///
/// # Returns
///
/// * `rad_vec`: Radius vector of the body *| in AU*. It is negative or
///   infinite for a true anomaly the body can never reach (beyond the
///   asymptotes of a parabola or a hyperbola).
///
/// # Arguments
///
/// * `true_anom` : True anomaly of the body *| in radians*
/// * `perih_dist`: Perihelion distance *| in AU*
/// * `ecc`       : Eccentricity of the orbit
#[inline]
pub fn rad_vec_frm_true_anom(true_anom: f64, perih_dist: f64, ecc: f64) -> f64 {
    perih_dist * (1.0 + ecc) / (1.0 + ecc * true_anom.cos())
}

/// Auxiliary constants `a, b, c, A, B, C` of an orbit, for computing
/// heliocentric rectangular equatorial coordinates
///
/// These are Meeus's quantities of eqs. 33.7 and 33.8. They depend only on
/// the inclination, the longitude of the ascending node and the obliquity
/// of the ecliptic, so they can be computed once for a whole ephemeris.
/// Do not confuse `a` with the semimajor axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EqOrbitConsts {
    /// Amplitude for the x coordinate (dimensionless, positive)
    pub a: f64,
    /// Amplitude for the y coordinate (dimensionless, positive)
    pub b: f64,
    /// Amplitude for the z coordinate (dimensionless, positive)
    pub c: f64,
    /// Phase angle for the x coordinate *| in radians*, in `[0, 2π)`
    pub A: f64,
    /// Phase angle for the y coordinate *| in radians*, in `[0, 2π)`
    pub B: f64,
    /// Phase angle for the z coordinate *| in radians*, in `[0, 2π)`
    pub C: f64,
}

/// Computes the auxiliary constants of an orbit used for heliocentric
/// rectangular equatorial coordinates
///
/// Meeus, eqs. 33.7 and 33.8.
///
/// # Arguments
///
/// * `long_asc_node`: Longitude of the ascending node *| in radians*
/// * `incl`         : Inclination of the orbit *| in radians*
/// * `oblq`         : Obliquity of the ecliptic for the equinox of the
///   elements *| in radians* (e.g. [`OBLQ_J2000`])
pub fn eq_orbit_consts(long_asc_node: f64, incl: f64, oblq: f64) -> EqOrbitConsts {
    let (sin_node, cos_node) = long_asc_node.sin_cos();
    let (sin_i, cos_i) = incl.sin_cos();
    let (sin_e, cos_e) = oblq.sin_cos();

    let F = cos_node;
    let G = sin_node * cos_e;
    let H = sin_node * sin_e;
    let P = -sin_node * cos_i;
    let Q = cos_node * cos_i * cos_e - sin_i * sin_e;
    let R = cos_node * cos_i * sin_e + sin_i * cos_e;

    EqOrbitConsts {
        a: F.hypot(P),
        b: G.hypot(Q),
        c: H.hypot(R),
        A: limit_to_two_pi(F.atan2(P)),
        B: limit_to_two_pi(G.atan2(Q)),
        C: limit_to_two_pi(H.atan2(R)),
    }
}

impl EqOrbitConsts {
    /// Computes heliocentric rectangular equatorial coordinates
    ///
    /// Meeus, eq. 33.9. The x-axis points to the vernal equinox, the
    /// z-axis to the north celestial pole (of the equinox of the elements).
    ///
    /// # Returns
    ///
    /// `(x, y, z)` *| in AU*
    ///
    /// # Arguments
    ///
    /// * `arg_perih`: Argument of the perihelion *| in radians*
    /// * `true_anom`: True anomaly of the body *| in radians*
    /// * `rad_vec`  : Radius vector of the body *| in AU*
    pub fn heliocent_rect(&self, arg_perih: f64, true_anom: f64, rad_vec: f64) -> (f64, f64, f64) {
        let u = arg_perih + true_anom;

        (
            rad_vec * self.a * (self.A + u).sin(),
            rad_vec * self.b * (self.B + u).sin(),
            rad_vec * self.c * (self.C + u).sin(),
        )
    }
}

/// Orientation vectors `P` and `Q` of an orbit in equatorial coordinates
///
/// `P` is the unit vector from the Sun towards the perihelion and `Q` the
/// unit vector in the orbital plane 90° ahead of it, in the direction of
/// motion, both in rectangular equatorial coordinates (Meeus, chapter 33,
/// the alternative to eqs. 33.7 to 33.9).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GaussVectors {
    /// x component of `P`
    pub Px: f64,
    /// y component of `P`
    pub Py: f64,
    /// z component of `P`
    pub Pz: f64,
    /// x component of `Q`
    pub Qx: f64,
    /// y component of `Q`
    pub Qy: f64,
    /// z component of `Q`
    pub Qz: f64,
}

/// Computes the orientation vectors `P` and `Q` of an orbit
///
/// Meeus, chapter 33 (the alternative to eqs. 33.7 and 33.8, given just
/// before eq. 33.11).
///
/// # Arguments
///
/// * `arg_perih`    : Argument of the perihelion *| in radians*
/// * `long_asc_node`: Longitude of the ascending node *| in radians*
/// * `incl`         : Inclination of the orbit *| in radians*
/// * `oblq`         : Obliquity of the ecliptic for the equinox of the
///   elements *| in radians*
pub fn gauss_vectors(arg_perih: f64, long_asc_node: f64, incl: f64, oblq: f64) -> GaussVectors {
    let (sin_w, cos_w) = arg_perih.sin_cos();
    let (sin_n, cos_n) = long_asc_node.sin_cos();
    let (sin_i, cos_i) = incl.sin_cos();
    let (sin_e, cos_e) = oblq.sin_cos();

    let px = cos_w * cos_n - sin_w * sin_n * cos_i;
    let p_tmp = cos_w * sin_n + sin_w * cos_n * cos_i;
    let qx = -sin_w * cos_n - cos_w * sin_n * cos_i;
    let q_tmp = cos_w * cos_n * cos_i - sin_w * sin_n;

    GaussVectors {
        Px: px,
        Py: cos_e * p_tmp - sin_e * sin_w * sin_i,
        Pz: sin_e * p_tmp + cos_e * sin_w * sin_i,
        Qx: qx,
        Qy: cos_e * q_tmp - sin_e * cos_w * sin_i,
        Qz: sin_e * q_tmp + cos_e * cos_w * sin_i,
    }
}

impl GaussVectors {
    /// Computes heliocentric rectangular equatorial coordinates
    ///
    /// `r (P cos(v) + Q sin(v))` (Meeus, chapter 33).
    ///
    /// # Returns
    ///
    /// `(x, y, z)` *| in AU*
    ///
    /// # Arguments
    ///
    /// * `true_anom`: True anomaly of the body *| in radians*
    /// * `rad_vec`  : Radius vector of the body *| in AU*
    pub fn heliocent_rect(&self, true_anom: f64, rad_vec: f64) -> (f64, f64, f64) {
        let (sin_v, cos_v) = true_anom.sin_cos();

        (
            rad_vec * (self.Px * cos_v + self.Qx * sin_v),
            rad_vec * (self.Py * cos_v + self.Qy * sin_v),
            rad_vec * (self.Pz * cos_v + self.Qz * sin_v),
        )
    }
}

/// Computes heliocentric rectangular ecliptic coordinates of a body
///
/// Meeus, chapter 33 (section on heliocentric ecliptical coordinates). The
/// x-axis points to the vernal equinox and the z-axis to the north pole of
/// the ecliptic of the elements.
///
/// # Returns
///
/// `(x, y, z)` *| in AU*
///
/// # Arguments
///
/// * `arg_perih`    : Argument of the perihelion *| in radians*
/// * `long_asc_node`: Longitude of the ascending node *| in radians*
/// * `incl`         : Inclination of the orbit *| in radians*
/// * `true_anom`    : True anomaly of the body *| in radians*
/// * `rad_vec`      : Radius vector of the body *| in AU*
pub fn heliocent_ecl_rect(
    arg_perih: f64,
    long_asc_node: f64,
    incl: f64,
    true_anom: f64,
    rad_vec: f64,
) -> (f64, f64, f64) {
    let (sin_u, cos_u) = (arg_perih + true_anom).sin_cos();
    let (sin_n, cos_n) = long_asc_node.sin_cos();
    let (sin_i, cos_i) = incl.sin_cos();

    (
        rad_vec * (cos_n * cos_u - sin_n * sin_u * cos_i),
        rad_vec * (sin_n * cos_u + cos_n * sin_u * cos_i),
        rad_vec * sin_i * sin_u,
    )
}

/// Computes the heliocentric ecliptic longitude and latitude of a body
///
/// Meeus, chapter 33 (section on heliocentric ecliptical coordinates).
///
/// # Returns
///
/// `(long, lat)`
///
/// * `long`: Heliocentric ecliptic longitude *| in radians*, in `[0, 2π)`
/// * `lat` : Heliocentric ecliptic latitude *| in radians*, in
///   `[-π/2, π/2]`
///
/// # Arguments
///
/// * `arg_perih`    : Argument of the perihelion *| in radians*
/// * `long_asc_node`: Longitude of the ascending node *| in radians*
/// * `incl`         : Inclination of the orbit *| in radians*
/// * `true_anom`    : True anomaly of the body *| in radians*
pub fn heliocent_ecl_coords(
    arg_perih: f64,
    long_asc_node: f64,
    incl: f64,
    true_anom: f64,
) -> (f64, f64) {
    let (x, y, z) = heliocent_ecl_rect(arg_perih, long_asc_node, incl, true_anom, 1.0);

    (limit_to_two_pi(y.atan2(x)), z.atan2(x.hypot(y)))
}

/// Computes geocentric equatorial coordinates from heliocentric
/// rectangular equatorial coordinates
///
/// Meeus, eq. 33.10: the geocentric vector is the sum of the body's
/// heliocentric vector and the Sun's geocentric vector.
///
/// # Returns
///
/// `(asc, dec, earth_dist)`
///
/// * `asc`       : Geocentric right ascension *| in radians*, in `[0, 2π)`
/// * `dec`       : Geocentric declination *| in radians*
/// * `earth_dist`: Distance from the Earth *| in AU*
///
/// # Arguments
///
/// * `body_rect`: Heliocentric rectangular equatorial coordinates
///   `(x, y, z)` of the body *| in AU*
/// * `sun_rect` : Geocentric rectangular equatorial coordinates `(X, Y, Z)`
///   of the Sun, for the same equinox *| in AU*
pub fn geocent_eq_frm_rect(
    body_rect: (f64, f64, f64),
    sun_rect: (f64, f64, f64),
) -> (f64, f64, f64) {
    let xi = sun_rect.0 + body_rect.0;
    let eta = sun_rect.1 + body_rect.1;
    let zeta = sun_rect.2 + body_rect.2;

    let dist = (xi * xi + eta * eta + zeta * zeta).sqrt();

    (
        limit_to_two_pi(eta.atan2(xi)),
        zeta.atan2(xi.hypot(eta)),
        dist,
    )
}

/// Computes the elongation of a body from the Sun from three distances
///
/// Uses `cos(ψ) = (R² + Δ² - r²) / (2RΔ)` (Meeus, eq. 33.11).
///
/// # Returns
///
/// * `elong`: Elongation *| in radians*, in `[0, π]`
///
/// # Arguments
///
/// * `sun_earth_dist` : Distance of the Sun from the Earth, `R` *| in AU*
/// * `body_earth_dist`: Distance of the body from the Earth, `Δ` *| in AU*
/// * `body_sun_dist`  : Distance of the body from the Sun, `r` *| in AU*
pub fn elong_frm_dist(sun_earth_dist: f64, body_earth_dist: f64, body_sun_dist: f64) -> f64 {
    let (R, D, r) = (sun_earth_dist, body_earth_dist, body_sun_dist);

    ((R * R + D * D - r * r) / (2.0 * R * D))
        .clamp(-1.0, 1.0)
        .acos()
}

/// Computes the phase angle of a body (the angle Sun–body–Earth) from three
/// distances
///
/// Uses `cos(β) = (r² + Δ² - R²) / (2rΔ)` (Meeus, eq. 33.12).
///
/// # Returns
///
/// * `phase_angl`: Phase angle *| in radians*, in `[0, π]`
///
/// # Arguments
///
/// * `body_sun_dist`  : Distance of the body from the Sun, `r` *| in AU*
/// * `body_earth_dist`: Distance of the body from the Earth, `Δ` *| in AU*
/// * `sun_earth_dist` : Distance of the Sun from the Earth, `R` *| in AU*
pub fn phase_angl_frm_dist(body_sun_dist: f64, body_earth_dist: f64, sun_earth_dist: f64) -> f64 {
    let (r, D, R) = (body_sun_dist, body_earth_dist, sun_earth_dist);

    ((r * r + D * D - R * R) / (2.0 * r * D))
        .clamp(-1.0, 1.0)
        .acos()
}

/// Computes the elongation of a body from the Sun from ecliptic coordinates
///
/// Uses `cos(ψ) = cos(β) cos(λ - λ0)` (Meeus, eq. 33.5), which neglects the
/// Sun's latitude (never more than about 1.2″).
///
/// # Returns
///
/// * `elong`: Elongation *| in radians*, in `[0, π]`
///
/// # Arguments
///
/// * `long`    : Geocentric ecliptic longitude of the body *| in radians*
/// * `lat`     : Geocentric ecliptic latitude of the body *| in radians*
/// * `sun_long`: Geocentric ecliptic longitude of the Sun *| in radians*
pub fn elong_frm_ecl(long: f64, lat: f64, sun_long: f64) -> f64 {
    (lat.cos() * (long - sun_long).cos())
        .clamp(-1.0, 1.0)
        .acos()
}

/// How a body moves along its orbit
///
/// Together with the orientation of the orbit in [`Elements`], this
/// determines the position of the body at any time (as an unperturbed
/// two-body orbit). Times are Julian Ephemeris Days (TT).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Motion {
    /// Elliptic orbit (`0 ≤ e < 1`), solved with Kepler's equation
    Elliptic {
        /// Semimajor axis *| in AU*
        semimaj_axis: f64,
        /// Eccentricity, `0 ≤ ecc < 1`
        ecc: f64,
        /// Mean motion *| in radians per day*
        mn_motion: f64,
        /// Mean anomaly at `epoch` *| in radians*
        mn_anom: f64,
        /// Epoch of `mn_anom`, in Julian Ephemeris Day
        epoch: f64,
    },
    /// Parabolic orbit (`e = 1`), solved with Barker's equation
    Parabolic {
        /// Perihelion distance *| in AU*
        perih_dist: f64,
        /// Time of passage through the perihelion, in Julian Ephemeris Day
        perih_time: f64,
    },
    /// Orbit with an eccentricity close to 1 (elliptic or hyperbolic),
    /// solved with Landgraf's method (see [`near_parabolic`])
    NearParabolic {
        /// Perihelion distance *| in AU*
        perih_dist: f64,
        /// Eccentricity, close to 1
        ecc: f64,
        /// Time of passage through the perihelion, in Julian Ephemeris Day
        perih_time: f64,
    },
}

impl Motion {
    /// Describes an elliptic motion from the semimajor axis, the
    /// eccentricity and the time of perihelion
    ///
    /// The mean motion is computed from the semimajor axis with
    /// [`elliptic::mn_motion`] (Meeus, eq. 33.6).
    ///
    /// # Arguments
    ///
    /// * `semimaj_axis`: Semimajor axis *| in AU*
    /// * `ecc`         : Eccentricity, `0 ≤ ecc < 1`
    /// * `perih_time`  : Time of passage through the perihelion, in Julian
    ///   Ephemeris Day
    pub fn elliptic_frm_perih_time(semimaj_axis: f64, ecc: f64, perih_time: f64) -> Motion {
        Motion::Elliptic {
            semimaj_axis,
            ecc,
            mn_motion: elliptic::mn_motion(semimaj_axis),
            mn_anom: 0.0,
            epoch: perih_time,
        }
    }

    /// Computes the true anomaly and the radius vector at a given time
    ///
    /// # Returns
    ///
    /// * `Some((true_anom, rad_vec))`: True anomaly *| in radians* and
    ///   radius vector *| in AU*
    /// * `None` if the elements are invalid (non-positive distance,
    ///   eccentricity out of range for the kind of motion, non-finite
    ///   values) or, for [`Motion::NearParabolic`], if Landgraf's method
    ///   does not converge
    ///
    /// # Arguments
    ///
    /// * `JD`: Julian Ephemeris Day
    pub fn true_anom_and_rad_vec(&self, JD: f64) -> Option<(f64, f64)> {
        match *self {
            Motion::Elliptic {
                semimaj_axis,
                ecc,
                mn_motion,
                mn_anom,
                epoch,
            } => {
                if semimaj_axis <= 0.0 || !semimaj_axis.is_finite() {
                    return None;
                }
                let M = mn_anom + mn_motion * (JD - epoch);
                let E = elliptic::ecc_anom_checked(M, ecc)?;
                Some((
                    elliptic::true_anom(E, ecc),
                    elliptic::rad_vec_frm_ecc_anom(E, semimaj_axis, ecc),
                ))
            }
            Motion::Parabolic {
                perih_dist,
                perih_time,
            } => {
                if perih_dist <= 0.0 || !perih_dist.is_finite() {
                    return None;
                }
                let (v, r) = parabolic::true_anom_and_rad_vec(JD, perih_time, perih_dist);
                (v.is_finite() && r.is_finite()).then_some((v, r))
            }
            Motion::NearParabolic {
                perih_dist,
                ecc,
                perih_time,
            } => {
                let (v, r) =
                    near_parabolic::true_anom_and_rad_vec(JD, perih_time, ecc, perih_dist, 1e-10)
                        .ok()?;
                (v.is_finite() && r.is_finite()).then_some((v, r))
            }
        }
    }
}

/// Orbital elements of a body moving around the Sun
///
/// The angular elements must refer to a standard equinox (for example
/// J2000.0). They are taken as constant: the orbit is unperturbed, which is
/// a good approximation for weeks or months around the epoch of osculating
/// elements, but not over long periods (Meeus, chapter 33).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Elements {
    /// Shape of the orbit and motion along it
    pub motion: Motion,
    /// Inclination *| in radians*, in `[0, π]`; above `π/2` the motion is
    /// retrograde
    pub incl: f64,
    /// Argument of the perihelion *| in radians*
    pub arg_perih: f64,
    /// Longitude of the ascending node *| in radians*
    pub long_asc_node: f64,
}

impl Elements {
    /// Computes the heliocentric rectangular equatorial coordinates of the
    /// body at a given time
    ///
    /// Meeus, eqs. 33.7 to 33.9.
    ///
    /// # Returns
    ///
    /// * `Some((x, y, z))` *| in AU*, referred to the equinox of the
    ///   elements
    /// * `None` if [`Motion::true_anom_and_rad_vec`] fails
    ///
    /// # Arguments
    ///
    /// * `JD`  : Julian Ephemeris Day
    /// * `oblq`: Obliquity of the ecliptic for the equinox of the elements
    ///   *| in radians*
    pub fn heliocent_eq_rect(&self, JD: f64, oblq: f64) -> Option<(f64, f64, f64)> {
        let (v, r) = self.motion.true_anom_and_rad_vec(JD)?;

        Some(
            eq_orbit_consts(self.long_asc_node, self.incl, oblq).heliocent_rect(
                self.arg_perih,
                v,
                r,
            ),
        )
    }

    /// Computes the heliocentric rectangular ecliptic coordinates of the
    /// body at a given time
    ///
    /// Meeus, chapter 33 (section on heliocentric ecliptical coordinates).
    ///
    /// # Returns
    ///
    /// * `Some((x, y, z))` *| in AU*, referred to the ecliptic and equinox
    ///   of the elements
    /// * `None` if [`Motion::true_anom_and_rad_vec`] fails
    ///
    /// # Arguments
    ///
    /// * `JD`: Julian Ephemeris Day
    pub fn heliocent_ecl_rect(&self, JD: f64) -> Option<(f64, f64, f64)> {
        let (v, r) = self.motion.true_anom_and_rad_vec(JD)?;

        Some(heliocent_ecl_rect(
            self.arg_perih,
            self.long_asc_node,
            self.incl,
            v,
            r,
        ))
    }
}

/// Geocentric equatorial position of a body computed from its orbital
/// elements
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeocentEqPos {
    /// Right ascension *| in radians*, in `[0, 2π)`
    pub asc: f64,
    /// Declination *| in radians*, in `[-π/2, π/2]`
    pub dec: f64,
    /// Distance from the Earth, `Δ`, at the time the light left the body
    /// *| in AU*
    pub earth_dist: f64,
    /// Distance from the Sun, `r`, at the time the light left the body
    /// *| in AU*
    pub sun_dist: f64,
    /// Light-time, `τ` *| in days*
    pub light_time: f64,
    /// Elongation from the Sun, `ψ` *| in radians*, in `[0, π]`
    pub elong: f64,
    /// Phase angle (angle Sun–body–Earth), `β` *| in radians*, in `[0, π]`
    pub phase_angl: f64,
}

/// Geocentric ecliptic position of a body computed from its orbital
/// elements
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeocentEclPos {
    /// Ecliptic longitude *| in radians*, in `[0, 2π)`
    pub long: f64,
    /// Ecliptic latitude *| in radians*, in `[-π/2, π/2]`
    pub lat: f64,
    /// Distance from the Earth, `Δ`, at the time the light left the body
    /// *| in AU*
    pub earth_dist: f64,
    /// Distance from the Sun, `r`, at the time the light left the body
    /// *| in AU*
    pub sun_dist: f64,
    /// Light-time, `τ` *| in days*
    pub light_time: f64,
    /// Elongation from the Sun, `ψ` *| in radians*, in `[0, π]`
    pub elong: f64,
    /// Phase angle (angle Sun–body–Earth), `β` *| in radians*, in `[0, π]`
    pub phase_angl: f64,
}

/// A rectangular vector `(x, y, z)`
type Vec3 = (f64, f64, f64);

#[inline]
fn dot(a: Vec3, b: Vec3) -> f64 {
    a.0 * b.0 + a.1 * b.1 + a.2 * b.2
}

#[inline]
fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

/// Angle between two vectors, in `[0, π]`, accurate for small and large
/// angles alike
fn angle_between(a: Vec3, b: Vec3) -> f64 {
    let cross = (
        a.1 * b.2 - a.2 * b.1,
        a.2 * b.0 - a.0 * b.2,
        a.0 * b.1 - a.1 * b.0,
    );

    norm(cross).atan2(dot(a, b))
}

/// Iterates on the light-time: `geocent(t)` returns the body's heliocentric
/// and geocentric vectors at time `t`
fn with_light_time<F>(JD: f64, geocent: F) -> Option<(Vec3, Vec3, f64)>
where
    F: Fn(f64) -> Option<(Vec3, Vec3)>,
{
    let mut tau = 0.0;
    let (mut helio, mut geo) = geocent(JD)?;

    for _ in 0..20 {
        let new_tau = light_time(norm(geo));
        if (new_tau - tau).abs() <= 1e-12 {
            break;
        }
        tau = new_tau;
        (helio, geo) = geocent(JD - tau)?;
    }

    Some((helio, geo, tau))
}

/// Computes the geocentric equatorial position of a body from its orbital
/// elements and the Sun's geocentric rectangular equatorial coordinates
///
/// This is Meeus's second method of chapter 33 (eqs. 33.7 to 33.12), with
/// the iteration for the light-time: the body's position is computed for
/// the time `JD - τ`, the Sun's position is kept at `JD`.
///
/// The result is the *astrometric* position: corrected for light-time, but
/// not for aberration or nutation, and referred to the equinox of the
/// elements. It can be compared directly with star catalogue positions for
/// that equinox.
///
/// # Returns
///
/// * `Some(GeocentEqPos)`
/// * `None` if [`Motion::true_anom_and_rad_vec`] fails
///
/// # Arguments
///
/// * `elements`: Orbital elements of the body
/// * `JD`      : Julian Ephemeris Day
/// * `sun_rect`: Geocentric rectangular equatorial coordinates `(X, Y, Z)`
///   of the Sun at `JD`, for the equinox of the elements *| in AU*
///   (see Meeus, chapter 26)
/// * `oblq`    : Obliquity of the ecliptic for the equinox of the elements
///   *| in radians* (e.g. [`OBLQ_J2000`])
pub fn geocent_eq_pos(
    elements: &Elements,
    JD: f64,
    sun_rect: (f64, f64, f64),
    oblq: f64,
) -> Option<GeocentEqPos> {
    let consts = eq_orbit_consts(elements.long_asc_node, elements.incl, oblq);

    let (helio, geo, tau) = with_light_time(JD, |t| {
        let (v, r) = elements.motion.true_anom_and_rad_vec(t)?;
        let h = consts.heliocent_rect(elements.arg_perih, v, r);
        Some((h, (h.0 + sun_rect.0, h.1 + sun_rect.1, h.2 + sun_rect.2)))
    })?;

    Some(GeocentEqPos {
        asc: limit_to_two_pi(geo.1.atan2(geo.0)),
        dec: geo.2.atan2(geo.0.hypot(geo.1)),
        earth_dist: norm(geo),
        sun_dist: norm(helio),
        light_time: tau,
        elong: angle_between(geo, sun_rect),
        phase_angl: angle_between(helio, geo),
    })
}

/// Computes the geocentric ecliptic position of a body from its orbital
/// elements and the Earth's heliocentric ecliptic coordinates
///
/// The body's heliocentric ecliptic position is computed from its elements
/// and combined with the Earth's as in Meeus's first method of chapter 33
/// (eqs. 33.1 to 33.4), with the iteration for the light-time: the body's
/// position is computed for the time `JD - τ`, the Earth's position is kept
/// at `JD`.
///
/// The result is corrected for light-time only, and is referred to the
/// ecliptic and equinox of the elements. (Passing the Earth's coordinates
/// for `JD - τ` instead also accounts for the aberration, as Meeus points
/// out in example 33.a.)
///
/// # Returns
///
/// * `Some(GeocentEclPos)`
/// * `None` if [`Motion::true_anom_and_rad_vec`] fails
///
/// # Arguments
///
/// * `elements`      : Orbital elements of the body
/// * `JD`            : Julian Ephemeris Day
/// * `earth_heliocent`: Heliocentric ecliptic coordinates `(L, B, R)` of
///   the Earth at `JD`, for the ecliptic and equinox of the elements:
///   longitude and latitude *| in radians*, radius vector *| in AU*
pub fn geocent_ecl_pos(
    elements: &Elements,
    JD: f64,
    earth_heliocent: (f64, f64, f64),
) -> Option<GeocentEclPos> {
    let (L0, B0, R0) = earth_heliocent;
    let earth = (
        R0 * B0.cos() * L0.cos(),
        R0 * B0.cos() * L0.sin(),
        R0 * B0.sin(),
    );

    let (helio, geo, tau) = with_light_time(JD, |t| {
        let h = elements.heliocent_ecl_rect(t)?;
        Some((h, (h.0 - earth.0, h.1 - earth.1, h.2 - earth.2)))
    })?;

    Some(GeocentEclPos {
        long: limit_to_two_pi(geo.1.atan2(geo.0)),
        lat: geo.2.atan2(geo.0.hypot(geo.1)),
        earth_dist: norm(geo),
        sun_dist: norm(helio),
        light_time: tau,
        elong: angle_between(geo, (-earth.0, -earth.1, -earth.2)),
        phase_angl: angle_between(helio, geo),
    })
}

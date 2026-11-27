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

//! The Earth
//!
//! The figure of the Earth as an ellipsoid of revolution (Meeus,
//! *Astronomical Algorithms*, 2nd ed., chapter 11): geocentric position of
//! an observer, geocentric latitude, radii of the parallel and of the
//! meridian, rotational velocity, and the distance between two points on
//! the surface.
//!
//! # Ellipsoids and units
//!
//! The free functions of this module (`eq_rad`, `geodesic_dist`,
//! `rad_of_parll_lat`, ...) use the **WGS 84** ellipsoid and return
//! **kilometres**, as they always have. The [`Ellipsoid`] type lets you
//! choose the ellipsoid and works in **metres**. Meeus' worked examples use
//! the **IAU 1976** ellipsoid (a = 6378.140 km, f = 1/298.257), available
//! as [`Ellipsoid::IAU_1976`]; WGS 84 differs from it by 3 m in the
//! equatorial radius, which changes the book's results by a few metres
//! (e.g. a radius of the parallel of 4746.999 km instead of 4747.001 km at
//! latitude 42°).
//!
//! Geographic longitudes are positive **west** of Greenwich (see
//! [`coords::GeographPoint`]); latitudes are geographic (geodetic).
//!
//! # Example
//!
//! Distance between the Paris and the U.S. Naval observatories (Meeus,
//! example 11.c):
//!
//! ```
//! use astro::{angle, coords, planet::earth};
//!
//! // Longitudes positive WEST: Paris is east of Greenwich.
//! let paris = coords::GeographPoint {
//!     long: angle::deg_frm_dms(-2, 20, 14.0).to_radians(),
//!     lat: angle::deg_frm_dms(48, 50, 11.0).to_radians(),
//! };
//! let washington = coords::GeographPoint {
//!     long: angle::deg_frm_dms(77, 3, 56.0).to_radians(),
//!     lat: angle::deg_frm_dms(38, 55, 17.0).to_radians(),
//! };
//!
//! // Andoyer's method on the IAU 1976 ellipsoid, in metres.
//! let s = earth::Ellipsoid::IAU_1976.geodesic_dist(&paris, &washington);
//! assert!((s / 1000.0 - 6181.63).abs() < 0.005);
//!
//! // Same on WGS 84, in kilometres.
//! let s = earth::geodesic_dist(&paris, &washington);
//! assert!((s - 6181.63).abs() < 0.005);
//!
//! // On a sphere of radius 6371 km.
//! let s = earth::approx_geodesic_dist(&paris, &washington);
//! assert!((s - 6166.0).abs() < 0.5);
//! ```

use crate::angle;
use crate::consts;
use crate::coords;
use crate::time;
use std::f64::consts::PI;

/// A reference ellipsoid of revolution for the figure of the Earth
///
/// Defined by its equatorial radius `a` and flattening `f = (a − b) / a`,
/// where `b` is the polar radius. Lengths are in **metres**.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipsoid {
    /// Equatorial radius (semi-major axis) *| in metres*
    pub eq_rad: f64,
    /// Flattening `f = (a − b) / a`
    pub flat: f64,
}

impl Ellipsoid {
    /// The ellipsoid adopted by the IAU in 1976: a = 6 378 140 m,
    /// f = 1/298.257. Meeus' examples in chapter 11 (and elsewhere) use it.
    pub const IAU_1976: Ellipsoid = Ellipsoid {
        eq_rad: 6_378_140.0,
        flat: 1.0 / 298.257,
    };

    /// The World Geodetic System 1984 ellipsoid (the GPS datum):
    /// a = 6 378 137 m, f = 1/298.257 223 563.
    pub const WGS_84: Ellipsoid = Ellipsoid {
        eq_rad: consts::wgs84::EQUATORIAL_RADIUS,
        flat: consts::wgs84::FLATTENING,
    };

    /// The World Geodetic System 1972 ellipsoid: a = 6 378 135 m,
    /// f = 1/298.26.
    pub const WGS_72: Ellipsoid = Ellipsoid {
        eq_rad: consts::wgs72::EQUATORIAL_RADIUS,
        flat: consts::wgs72::FLATTENING,
    };

    /// Returns the polar radius `b = a (1 − f)` *| in metres*
    #[inline]
    pub fn pol_rad(&self) -> f64 {
        self.eq_rad * (1.0 - self.flat)
    }

    /// Returns the eccentricity of the meridian, `e = √(2f − f²)`
    #[inline]
    pub fn ecc(&self) -> f64 {
        (self.flat * (2.0 - self.flat)).sqrt()
    }

    /// Computes `ρ sin φ′` and `ρ cos φ′` for an observer
    ///
    /// `ρ` is the observer's distance from the centre of the Earth, in units
    /// of the equatorial radius, and `φ′` the observer's geocentric
    /// latitude. These are the quantities needed for diurnal parallax,
    /// eclipses and occultations (Meeus, chapter 11).
    ///
    /// # Returns
    ///
    /// `(rho_sin_phi, rho_cos_phi)`: `ρ sin φ′` has the sign of the
    /// latitude; `ρ cos φ′` is never negative.
    ///
    /// # Arguments
    ///
    /// * `geograph_lat`: Observer's geographic latitude *| in radians*
    /// * `height`: Observer's height above the ellipsoid (sea level)
    ///   *| in metres*
    pub fn rho_sin_cos_phi(&self, geograph_lat: f64, height: f64) -> (f64, f64) {
        let b_over_a = 1.0 - self.flat;
        let (sin_lat, cos_lat) = geograph_lat.sin_cos();
        // Reduced (parametric) latitude u, with tan u = (b/a) tan φ.
        let u = (b_over_a * sin_lat).atan2(cos_lat);
        let x = height / self.eq_rad;

        (b_over_a * u.sin() + x * sin_lat, u.cos() + x * cos_lat)
    }

    /// Computes the geocentric latitude of a point at sea level
    ///
    /// `tan φ′ = (b²/a²) tan φ` (Meeus, chapter 11). For an observer above
    /// sea level, use `atan2` of the two values returned by
    /// [`rho_sin_cos_phi`](Self::rho_sin_cos_phi).
    ///
    /// # Returns
    ///
    /// * `geocent_lat`: Geocentric latitude, in `[−π/2, π/2]`
    ///   *| in radians*
    ///
    /// # Arguments
    ///
    /// * `geograph_lat`: Geographic latitude *| in radians*
    pub fn geocent_lat(&self, geograph_lat: f64) -> f64 {
        let b_over_a = 1.0 - self.flat;
        let (sin_lat, cos_lat) = geograph_lat.sin_cos();
        (b_over_a * b_over_a * sin_lat).atan2(cos_lat)
    }

    /// Computes the distance from the centre of the Earth to a point at sea
    /// level, in units of the equatorial radius
    ///
    /// This is the exact value; [`rho`] is Meeus' short series for it.
    ///
    /// # Arguments
    ///
    /// * `geograph_lat`: Geographic latitude *| in radians*
    pub fn rho(&self, geograph_lat: f64) -> f64 {
        let (rho_sin_phi, rho_cos_phi) = self.rho_sin_cos_phi(geograph_lat, 0.0);
        rho_sin_phi.hypot(rho_cos_phi)
    }

    /// Computes the geographic and geocentric latitudes at which their
    /// difference is greatest
    ///
    /// That happens where `tan φ₀ = a/b`, `tan φ₀′ = b/a` (so
    /// `φ₀ + φ₀′ = 90°`), i.e. at a reduced latitude of 45° (Meeus,
    /// chapter 11). For the IAU 1976 ellipsoid φ₀ = 45°05′46.36″ and
    /// φ₀′ = 44°54′13.64″.
    ///
    /// # Returns
    ///
    /// `(geograph_lat, geocent_lat)` *| in radians*
    pub fn lat_of_max_geocent_diff(&self) -> (f64, f64) {
        let b_over_a = 1.0 - self.flat;
        (1.0_f64.atan2(b_over_a), b_over_a.atan())
    }

    /// Computes the radius of the parallel of a latitude
    ///
    /// `R_p = a cos φ / √(1 − e² sin² φ)` (Meeus, chapter 11), the distance
    /// of a point at sea level from the Earth's axis.
    ///
    /// # Returns
    ///
    /// * `rad`: Radius of the parallel *| in metres*
    ///
    /// # Arguments
    ///
    /// * `geograph_lat`: Geographic latitude *| in radians*
    pub fn rad_of_parll_lat(&self, geograph_lat: f64) -> f64 {
        let e2 = self.flat * (2.0 - self.flat);
        let (sin_lat, cos_lat) = geograph_lat.sin_cos();
        self.eq_rad * cos_lat / (1.0 - e2 * sin_lat * sin_lat).sqrt()
    }

    /// Computes the length of one degree of longitude at a latitude
    ///
    /// `(π/180) R_p` (Meeus, chapter 11).
    ///
    /// # Returns
    ///
    /// * `len`: Length of 1° of longitude along the parallel *| in metres*
    ///
    /// # Arguments
    ///
    /// * `geograph_lat`: Geographic latitude *| in radians*
    #[inline]
    pub fn len_of_deg_long(&self, geograph_lat: f64) -> f64 {
        self.rad_of_parll_lat(geograph_lat).to_radians()
    }

    /// Computes the linear velocity of a point at sea level due to the
    /// rotation of the Earth
    ///
    /// `ω R_p`, with `ω` from [`rot_angular_velocity`] (Meeus, chapter 11).
    ///
    /// # Returns
    ///
    /// * `lin_vel`: Linear velocity *| in metres per second*
    ///
    /// # Arguments
    ///
    /// * `geograph_lat`: Geographic latitude *| in radians*
    #[inline]
    pub fn linear_velocity_at_lat(&self, geograph_lat: f64) -> f64 {
        rot_angular_velocity() * self.rad_of_parll_lat(geograph_lat)
    }

    /// Computes the radius of curvature of the meridian at a latitude
    ///
    /// `R_m = a (1 − e²) / (1 − e² sin² φ)^(3/2)` (Meeus, chapter 11). It
    /// is smallest at the equator, `a (1 − e²)`, and largest at the poles,
    /// `a / √(1 − e²)`.
    ///
    /// # Returns
    ///
    /// * `rad`: Radius of curvature of the meridian *| in metres*
    ///
    /// # Arguments
    ///
    /// * `geograph_lat`: Geographic latitude *| in radians*
    pub fn rad_curv_of_meridian(&self, geograph_lat: f64) -> f64 {
        let e2 = self.flat * (2.0 - self.flat);
        let sin_lat = geograph_lat.sin();
        self.eq_rad * (1.0 - e2) / (1.0 - e2 * sin_lat * sin_lat).powf(1.5)
    }

    /// Computes the length of one degree of latitude at a latitude
    ///
    /// `(π/180) R_m` (Meeus, chapter 11).
    ///
    /// # Returns
    ///
    /// * `len`: Length of 1° of latitude along the meridian *| in metres*
    ///
    /// # Arguments
    ///
    /// * `geograph_lat`: Geographic latitude *| in radians*
    #[inline]
    pub fn len_of_deg_lat(&self, geograph_lat: f64) -> f64 {
        self.rad_curv_of_meridian(geograph_lat).to_radians()
    }

    /// Computes the geodesic distance between two points at sea level,
    /// with H. Andoyer's method
    ///
    /// Meeus, chapter 11. The relative error is of the order of the square
    /// of the flattening: a few tens of metres at most for points that are
    /// not nearly antipodal (Meeus quotes about 50 m for example 11.c).
    ///
    /// The method breaks down near antipodal points: within about a degree
    /// of the antipode the error can reach tens of kilometres. Identical
    /// points give `0`; exactly antipodal points (to within rounding) give
    /// half the length of a meridian, the true shortest distance between
    /// them on an oblate ellipsoid.
    ///
    /// # Returns
    ///
    /// * `dist`: Distance along the surface *| in metres*
    ///
    /// # Arguments
    ///
    /// * `p1`, `p2`: The two points. Only the difference of their longitudes
    ///   matters, so either longitude convention works if both points use
    ///   the same one.
    pub fn geodesic_dist(&self, p1: &coords::GeographPoint, p2: &coords::GeographPoint) -> f64 {
        andoyer(self, p1, p2).dist
    }

    /// Half the length of a meridian (pole to pole), from the series in the
    /// third flattening `n = (a − b) / (a + b)`.
    fn half_meridian(&self) -> f64 {
        let n = self.flat / (2.0 - self.flat);
        let n2 = n * n;
        PI * (self.eq_rad + self.pol_rad()) / 2.0 * (1.0 + n2 / 4.0 + n2 * n2 / 64.0)
    }
}

/// Intermediate quantities of Andoyer's method (Meeus, example 11.c).
/// Only `dist` is used outside the unit tests.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(not(test), allow(dead_code))]
struct Andoyer {
    s: f64,
    c: f64,
    omega: f64,
    r: f64,
    d: f64,
    dist: f64,
}

fn andoyer(ell: &Ellipsoid, p1: &coords::GeographPoint, p2: &coords::GeographPoint) -> Andoyer {
    let (sin_f, cos_f) = ((p1.lat + p2.lat) / 2.0).sin_cos();
    let (sin_g, cos_g) = ((p1.lat - p2.lat) / 2.0).sin_cos();
    let (sin_l, cos_l) = ((p1.long - p2.long) / 2.0).sin_cos();

    let s = (sin_g * cos_l).powi(2) + (cos_f * sin_l).powi(2);
    let c = (cos_g * cos_l).powi(2) + (sin_f * sin_l).powi(2);

    if s == 0.0 {
        // The same point.
        return Andoyer {
            s,
            c,
            omega: 0.0,
            r: 0.0,
            d: 0.0,
            dist: 0.0,
        };
    }
    if c < 1e-20 {
        // Antipodal points, to within rounding (C = 0 at the antipode; 1e-20
        // is about 0.6 mm from it on the Earth). Andoyer's correction terms
        // become 0 · ∞ there and their limit depends on the direction of
        // approach, so return the exact answer instead.
        let dist = ell.half_meridian();
        return Andoyer {
            s,
            c,
            omega: PI / 2.0,
            r: 0.0,
            d: PI * ell.eq_rad,
            dist,
        };
    }

    let omega = (s / c).sqrt().atan();
    let r = (s * c).sqrt() / omega;
    let d = 2.0 * omega * ell.eq_rad;
    let h1 = (3.0 * r - 1.0) / (2.0 * c);
    let h2 = (3.0 * r + 1.0) / (2.0 * s);
    let f = ell.flat;

    let dist = d * (1.0 + f * h1 * (sin_f * cos_g).powi(2) - f * h2 * (cos_f * sin_g).powi(2));

    Andoyer {
        s,
        c,
        omega,
        r,
        d,
        dist,
    }
}

/// Returns the flattening factor of the Earth (WGS 84)
///
/// Reference: World Geodetic System 1984. Meeus uses the IAU 1976 value
/// 1/298.257; see [`Ellipsoid`].
#[inline(always)]
pub fn flat_fac() -> f64 {
    Ellipsoid::WGS_84.flat
}

/// Returns the equatorial radius of the Earth (WGS 84) *| in kilometres*
///
/// Reference: World Geodetic System 1984. Meeus uses the IAU 1976 value
/// 6378.140 km; see [`Ellipsoid`].
#[inline(always)]
pub fn eq_rad() -> f64 {
    Ellipsoid::WGS_84.eq_rad / 1000.0
}

/// Returns the polar radius of the Earth (WGS 84) *| in kilometres*
///
/// `b = a (1 − f)`, from [`eq_rad`] and [`flat_fac`].
#[inline]
pub fn pol_rad() -> f64 {
    eq_rad() * (1.0 - flat_fac())
}

/// Returns the eccentricity of the Earth's meridian (WGS 84)
///
/// `e = √(2f − f²)`.
#[inline]
pub fn ecc_of_meridian() -> f64 {
    Ellipsoid::WGS_84.ecc()
}

/// Computes a low accuracy geodesic distance between two points
/// on the Earth's surface *| in kilometres*
///
/// Treats the Earth as a sphere of radius 6371 km (Meeus, eqs. 11.1 and
/// 11.2). The error can reach several tens of kilometres over long
/// distances (16 km in Meeus' example 11.c); use [`geodesic_dist`] for
/// better accuracy.
///
/// # Arguments
///
/// * `p1`: `GeographPoint` 1
/// * `p2`: `GeographPoint` 2
#[inline]
pub fn approx_geodesic_dist(p1: &coords::GeographPoint, p2: &coords::GeographPoint) -> f64 {
    6371.0 * p1.anglr_sepr(p2)
}

/// Computes a high accuracy geodesic distance between two points on the
/// Earth's surface (WGS 84) *| in kilometres*
///
/// Andoyer's method (Meeus, chapter 11) on the WGS 84 ellipsoid; see
/// [`Ellipsoid::geodesic_dist`] for the accuracy, the limiting cases and
/// other ellipsoids.
///
/// # Arguments
///
/// * `p1`: `GeographPoint` 1
/// * `p2`: `GeographPoint` 2
pub fn geodesic_dist(p1: &coords::GeographPoint, p2: &coords::GeographPoint) -> f64 {
    Ellipsoid::WGS_84.geodesic_dist(p1, p2) / 1000.0
}

/// Computes two quantities that are used elsewhere in the library
/// (WGS 84)
///
/// # Returns
///
/// `Rho` here denotes the distance from the Earth's center to the
/// observer, in units of the equatorial radius, and `Phi` here denotes the
/// geocentric latitude of the observer (Meeus, chapter 11).
///
/// `(rho_sin_phi, rho_cos_phi)`
///
/// * `rho_sin_phi`: Rho sin phi, positive in the northern hemisphere
/// * `rho_cos_phi`: Rho cos phi, never negative
///
/// See [`Ellipsoid::rho_sin_cos_phi`] for other ellipsoids.
///
/// # Arguments
///
/// * `geograph_lat`: Observer's geographic latitude *| in radians*
/// * `height`      : Observer's height above sea level *| in metres*
pub fn rho_sin_cos_phi(geograph_lat: f64, height: f64) -> (f64, f64) {
    Ellipsoid::WGS_84.rho_sin_cos_phi(geograph_lat, height)
}

/// Computes the distance from the Earth's center to a point on the
/// ellipsoid
///
/// Meeus' short series for the IAU 1976 ellipsoid, accurate to about
/// 1e-7; [`Ellipsoid::rho`] gives the exact value for any ellipsoid.
///
/// # Returns
///
/// * `rho`: Distance from the Earth's center to the point on the
///   ellipsoid (*fraction of the equatorial radius*)
///
/// # Arguments
///
/// * `geograph_lat`: Geographic latitude of a point on the
///   ellipsoid *| in radians*
pub fn rho(geograph_lat: f64) -> f64 {
    0.9983271 + 0.0016764 * (2.0 * geograph_lat).cos() - 0.0000035 * (4.0 * geograph_lat).cos()
}

/// Returns the rotational angular velocity of the Earth
/// *| in radian per second*
///
/// 7.292 114 992 × 10⁻⁵ rad/s, with respect to the stars (not to the moving
/// equinox), at the epoch 1996.5 (Meeus, chapter 11). It decreases slowly
/// as the Earth's rotation slows down.
#[inline(always)]
pub fn rot_angular_velocity() -> f64 {
    0.00007292114992
}

/// Computes the radius of the parallel of a latitude (WGS 84)
///
/// # Returns
///
/// * `rad`: Radius of the parallel of the latitude
///   *| in kilometres*
///
/// See [`Ellipsoid::rad_of_parll_lat`] for the formula and other
/// ellipsoids.
///
/// # Arguments
///
/// * `geograph_lat`: Geographic latitude of a point
///   on the ellipsoid *| in radians*
pub fn rad_of_parll_lat(geograph_lat: f64) -> f64 {
    Ellipsoid::WGS_84.rad_of_parll_lat(geograph_lat) / 1000.0
}

/// Computes the linear velocity of a point at a latitude (WGS 84)
///
/// # Returns
///
/// * `lin_vel`: Linear velocity at the latitude due to the rotation of
///   the Earth *| in kilometres per second*
///
/// # Arguments
///
/// * `geograph_lat`: Geographic latitude of a point on
///   the ellipsoid *| in radians*
#[inline(always)]
pub fn linear_velocity_at_lat(geograph_lat: f64) -> f64 {
    rot_angular_velocity() * rad_of_parll_lat(geograph_lat)
}

/// Computes the radius of curvature of the Earth's meridian
/// at a latitude (WGS 84)
///
/// # Returns
///
/// * `rad`: Radius of curvature of the Earth's meridian at the
///   latitude *| in kilometres*
///
/// See [`Ellipsoid::rad_curv_of_meridian`] for the formula and other
/// ellipsoids.
///
/// # Arguments
///
/// * `geograph_lat`: Geographic latitude of a point on the
///   ellipsoid *| in radians*
pub fn rad_curv_of_meridian(lat: f64) -> f64 {
    Ellipsoid::WGS_84.rad_curv_of_meridian(lat) / 1000.0
}

/// Computes the difference between the geographic latitude and
/// geocentric latitude
///
/// Meeus' series `φ − φ′ = 692.73″ sin 2φ − 1.16″ sin 4φ` for the IAU 1976
/// ellipsoid, accurate to about 0.01″; [`Ellipsoid::geocent_lat`] gives
/// the exact geocentric latitude. Earlier versions of this crate returned
/// the value in degrees despite documenting radians.
///
/// # Returns
///
/// * `diff`: Geographic latitude minus geocentric latitude
///   *| in radians*
///
/// # Arguments
///
/// * `geograph_lat`: Geographic latitude *| in radians*
pub fn geograph_geocent_lat_diff(geograph_lat: f64) -> f64 {
    (angle::deg_frm_dms(0, 0, 692.73) * (2.0 * geograph_lat).sin()
        - angle::deg_frm_dms(0, 0, 1.16) * (4.0 * geograph_lat).sin())
    .to_radians()
}

/// Computes the equation of time *| in radians*
///
/// # Arguments
///
/// * `JD`      : Julian (Ephemeris) day
/// * `sun_asc` : Right ascension of the Sun *| in radians*
/// * `nut_log` : Nutation correction for longitude *| in radians*
/// * `tru_oblq`: True obliquity of the ecliptic *| in radians*
pub fn equation_of_time(JD: f64, sun_asc: f64, nut_long: f64, tru_oblq: f64) -> f64 {
    let t = time::julian_mill(JD);
    let L = angle::limit_to_360(
        280.4664567
            + t * (360007.6982779
                + t * (0.030328 + t * (1.0 / 49931.0 - t * (1.0 / 15300.0 + t / 2000000.0)))),
    );

    (L - 0.0057183 - sun_asc.to_degrees() + nut_long.to_degrees() * tru_oblq.cos()).to_radians()
}

/// Computes the angle between diurnal path and the horizon
///
/// The acute angle `J` at which a body's daily path (its parallel of
/// declination) crosses the horizon when it rises or sets (Meeus,
/// chapter 14): `B = tan δ tan φ`, `C = √(1 − B²)`,
/// `tan J = C cos δ / |tan φ|`. The declination is taken as constant and
/// refraction is ignored. On the equator `J = 90°`; for a body on the
/// celestial equator `J = 90° − |φ|`. The result is the same for northern
/// and southern observers at the same `|φ|`.
///
/// Earlier versions of this crate returned `π − J` for southern latitudes.
///
/// # Returns
///
/// * `angl`: Angle between the diurnal path of a celestial
///   body and the horizon, in `[0, π/2]` *| in radians*. `NaN` if the body
///   never rises or never sets (`|tan δ tan φ| > 1`).
///
/// # Arguments
///
/// * `dec`         : Declination of the celestial body
///   *| in radians*
/// * `observer_lat`: Observer's geographic latitude
///   *| in radians*
pub fn angl_betwn_diurnal_path_and_hz(dec: f64, observer_lat: f64) -> f64 {
    // tan J = C cos δ / |tan φ|, multiplied through by |cos φ| (so the
    // equator needs no special case): C cos δ |cos φ| / |sin φ|.
    let B = dec.tan() * observer_lat.tan();
    let C = (1.0 - B * B).sqrt();
    let (sin_p, cos_p) = observer_lat.sin_cos();

    (C * dec.cos() * cos_p.abs()).atan2(sin_p.abs())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn point(long_deg: f64, lat_deg: f64) -> coords::GeographPoint {
        coords::GeographPoint {
            long: long_deg.to_radians(),
            lat: lat_deg.to_radians(),
        }
    }

    // Meeus, example 11.c: the intermediate quantities of Andoyer's method
    // on the IAU 1976 ellipsoid.
    #[test]
    fn andoyer_intermediates_example_11_c() {
        let paris = point(
            -angle::deg_frm_dms(2, 20, 14.0),
            angle::deg_frm_dms(48, 50, 11.0),
        );
        let washington = point(
            angle::deg_frm_dms(77, 3, 56.0),
            angle::deg_frm_dms(38, 55, 17.0),
        );
        let a = andoyer(&Ellipsoid::IAU_1976, &paris, &washington);

        assert!((a.s - 0.21642696).abs() < 5e-9, "{}", a.s);
        assert!((a.c - 0.78357304).abs() < 5e-9, "{}", a.c);
        assert!((a.omega.to_degrees() - 27.724274).abs() < 5e-7);
        assert!((a.omega - 0.48387987).abs() < 5e-9);
        assert!((a.r - 0.8510555).abs() < 5e-8);
        assert!((a.d / 1000.0 - 6172.507).abs() < 5e-4);
        assert!((a.dist / 1000.0 - 6181.63).abs() < 5e-3);
    }

    #[test]
    fn half_meridian_length() {
        // WGS 84 meridian quadrant: 10 001 965.729 m.
        let q = Ellipsoid::WGS_84.half_meridian() / 2.0;
        assert!((q - 10_001_965.729).abs() < 1e-3, "{q}");
    }
}

//! Ephemeris for physical observations of the Moon (Meeus, *Astronomical
//! Algorithms*, 2nd ed., chapter 53): librations, position angle of the
//! axis, topocentric librations, the selenographic position of the Sun and
//! sunrise and sunset at a point on the Moon.
//!
//! # Conventions
//!
//! * Selenographic longitudes are measured from the mean central meridian of
//!   the disk and are positive towards Mare Crisium (that is, towards the
//!   west of the sky as seen from the Earth). Selenographic latitudes are
//!   positive towards the north (towards Mare Serenitatis).
//! * A positive libration in longitude means that the mean centre of the
//!   disk is displaced to the east on the sky, so that more of the western
//!   limb (the Mare Crisium side) is visible; a positive libration in
//!   latitude exposes more of the northern limb.
//! * Librations are returned in `(-π, π]`. Position angles are measured
//!   eastward from the north point of the disk, in `[0, 2π)`.
//! * The mean lunar equator is inclined by `I = 1°32′32.7″` (IAU) to the
//!   ecliptic; the physical libration uses the series of D. H. Eckhardt
//!   (1981) as abridged by Meeus.
//!
//! # Example
//!
//! ```
//! use astro::lunar::libration;
//!
//! // Meeus, examples 53.a and 53.b: 1992 April 12, 0h TD
//! let eph = libration::physical_ephemeris(2448724.5);
//! assert!((eph.libr.total_long.to_degrees() - -1.23).abs() < 0.005);
//! assert!((eph.libr.total_lat.to_degrees() - 4.20).abs() < 0.005);
//! assert!((eph.pos_angl_axis.to_degrees() - 15.08).abs() < 0.005);
//! assert!((eph.sun.colong.to_degrees() - 22.11).abs() < 0.005);
//! ```

use std::f64::consts::{FRAC_PI_2, PI};

use super::inc_of_mn_lunar_eq;
use super::position::{self, FundamentalArgs, norm_2pi, norm_pi};
use crate::coords;
use crate::parallax;
use crate::time;

/// Mean daily motion of the Sun's selenographic colongitude, equal to that
/// of the Moon's mean elongation `D` *| in degrees per day*.
const COLONG_RATE_DEG_PER_DAY: f64 = 12.190749;

/// Optical, physical and total librations of the Moon.
///
/// All values are *| in radians*, in `(-π, π]`. Each total is the sum of the
/// optical and the physical part.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Librations {
    /// Optical libration in longitude `l′`
    pub optical_long: f64,
    /// Optical libration in latitude `b′`
    pub optical_lat: f64,
    /// Physical libration in longitude `l″`
    pub physical_long: f64,
    /// Physical libration in latitude `b″`
    pub physical_lat: f64,
    /// Total libration in longitude `l = l′ + l″`, the selenographic
    /// longitude of the Earth
    pub total_long: f64,
    /// Total libration in latitude `b = b′ + b″`, the selenographic
    /// latitude of the Earth
    pub total_lat: f64,
}

/// Selenographic position of the Sun: the point of the lunar surface that
/// has the Sun in its zenith (Meeus, chapter 53).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SelenographicSun {
    /// Selenographic longitude `l0` of the subsolar point *| in radians*,
    /// in `(-π, π]`
    pub long: f64,
    /// Selenographic latitude `b0` of the subsolar point *| in radians*
    pub lat: f64,
    /// Selenographic colongitude of the Sun, `c0 = 90° − l0`
    /// *| in radians*, in `[0, 2π)`. It increases with time by about
    /// 12.19° per day; it is near 0°, 90°, 180° and 270° at First
    /// Quarter, Full Moon, Last Quarter and New Moon.
    pub colong: f64,
}

/// Physical ephemeris of the Moon for an instant, as seen from the centre
/// of the Earth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhysicalEphemeris {
    /// Optical, physical and total librations
    pub libr: Librations,
    /// Position angle `P` of the Moon's axis of rotation, eastward from the
    /// north point of the disk *| in radians*, in `[0, 2π)` (always within
    /// about 25° of north)
    pub pos_angl_axis: f64,
    /// Selenographic position of the Sun
    pub sun: SelenographicSun,
}

/// Librations and position angle of the axis as seen by an observer on the
/// surface of the Earth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TopocentLibrations {
    /// Topocentric (total) libration in longitude *| in radians*, in
    /// `(-π, π]`
    pub long: f64,
    /// Topocentric (total) libration in latitude *| in radians*
    pub lat: f64,
    /// Topocentric position angle of the axis *| in radians*, in `[0, 2π)`
    pub pos_angl_axis: f64,
}

/// Optical libration and the auxiliary angles `W` and `A` of eq. 53.1.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Optical {
    pub W: f64,
    pub A: f64,
    pub long: f64,
    pub lat: f64,
}

/// Eckhardt's quantities `ρ`, `σ`, `τ`, in radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct RhoSigmaTau {
    pub rho: f64,
    pub sigma: f64,
    pub tau: f64,
}

/// Eq. 53.1 for a body (the Earth, or the Sun) at geocentric ecliptic
/// longitude `long` (without nutation) and latitude `lat` as seen from the
/// Moon's centre in the reverse direction.
pub(crate) fn optical(node: f64, F: f64, long: f64, lat: f64) -> Optical {
    let I = inc_of_mn_lunar_eq();
    let W = long - node;
    let (sin_W, cos_W) = W.sin_cos();
    let (sin_b, cos_b) = lat.sin_cos();
    let (sin_I, cos_I) = I.sin_cos();

    let A = (sin_W * cos_b * cos_I - sin_b * sin_I).atan2(cos_W * cos_b);
    let b1 = (-sin_W * cos_b * sin_I - sin_b * cos_I)
        .clamp(-1.0, 1.0)
        .asin();

    Optical {
        W: norm_2pi(W),
        A: norm_2pi(A),
        long: norm_pi(A - F),
        lat: b1,
    }
}

/// `K1` and `K2` (degrees) of the physical libration, for `T` centuries.
#[inline]
fn K1_K2(T: f64) -> (f64, f64) {
    (119.75 + 131.849 * T, 72.56 + 20.186 * T)
}

/// Eckhardt's `ρ`, `σ`, `τ` (Meeus, chapter 53), in radians.
pub(crate) fn rho_sigma_tau(T: f64, a: &FundamentalArgs, node: f64) -> RhoSigmaTau {
    let (D, M, M1, F, E) = (
        a.mn_elong,
        a.sun_mn_anom,
        a.mn_anom,
        a.arg_of_lat,
        a.ecc_fac,
    );
    let (K1, K2) = K1_K2(T);
    let (K1, K2) = (K1.to_radians(), K2.to_radians());
    let x2F = 2.0 * F;
    let x2D = 2.0 * D;

    let rho = -0.02752 * M1.cos() - 0.02245 * F.sin() + 0.00684 * (M1 - x2F).cos()
        - 0.00293 * x2F.cos()
        - 0.00085 * (x2F - x2D).cos()
        - 0.00054 * (M1 - x2D).cos()
        - 0.00020 * (M1 + F).sin()
        - 0.00020 * (M1 + x2F).cos()
        - 0.00020 * (M1 - F).cos()
        + 0.00014 * (M1 + x2F - x2D).cos();

    let sigma = -0.02816 * M1.sin() + 0.02244 * F.cos()
        - 0.00682 * (M1 - x2F).sin()
        - 0.00279 * x2F.sin()
        - 0.00083 * (x2F - x2D).sin()
        + 0.00069 * (M1 - x2D).sin()
        + 0.00040 * (M1 + F).cos()
        - 0.00025 * (2.0 * M1).sin()
        - 0.00023 * (M1 + x2F).sin()
        + 0.00020 * (M1 - F).cos()
        + 0.00019 * (M1 - F).sin()
        + 0.00013 * (M1 + x2F - x2D).sin()
        - 0.00010 * (M1 - 3.0 * F).cos();

    let tau = 0.02520 * E * M.sin() + 0.00473 * (2.0 * M1 - x2F).sin() - 0.00467 * M1.sin()
        + 0.00396 * K1.sin()
        + 0.00276 * (2.0 * M1 - x2D).sin()
        + 0.00196 * node.sin()
        - 0.00183 * (M1 - F).cos()
        + 0.00115 * (M1 - x2D).sin()
        - 0.00096 * (M1 - D).sin()
        + 0.00046 * (x2F - x2D).sin()
        - 0.00039 * (M1 - F).sin()
        - 0.00032 * (M1 - M - D).sin()
        + 0.00027 * (2.0 * M1 - M - x2D).sin()
        + 0.00023 * K2.sin()
        - 0.00014 * x2D.sin()
        + 0.00014 * (2.0 * M1 - x2F).cos()
        - 0.00012 * (M1 - x2F).sin()
        - 0.00012 * (2.0 * M1).sin()
        + 0.00011 * (2.0 * M1 - 2.0 * M - x2D).sin();

    RhoSigmaTau {
        rho: rho.to_radians(),
        sigma: sigma.to_radians(),
        tau: tau.to_radians(),
    }
}

/// Eq. 53.2: physical librations `(l″, b″)` in radians.
#[inline]
pub(crate) fn physical(rst: &RhoSigmaTau, A: f64, optical_lat: f64) -> (f64, f64) {
    let (sin_A, cos_A) = A.sin_cos();
    (
        -rst.tau + (rst.rho * cos_A + rst.sigma * sin_A) * optical_lat.tan(),
        rst.sigma * cos_A - rst.rho * sin_A,
    )
}

/// Auxiliary quantities `V`, `X`, `Y` of the position angle of the axis.
#[inline]
fn pos_angl_aux(node: f64, nut_long: f64, rst: &RhoSigmaTau, true_oblq: f64) -> (f64, f64, f64) {
    let I = inc_of_mn_lunar_eq();
    let V = node + nut_long + rst.sigma / I.sin();
    let (sin_Ir, cos_Ir) = (I + rst.rho).sin_cos();
    let (sin_e, cos_e) = true_oblq.sin_cos();
    let X = sin_Ir * V.sin();
    let Y = sin_Ir * V.cos() * cos_e - cos_Ir * sin_e;
    (V, X, Y)
}

/// Position angle of the axis (Meeus, chapter 53), in `[0, 2π)`.
pub(crate) fn pos_angl(
    node: f64,
    nut_long: f64,
    rst: &RhoSigmaTau,
    total_lat: f64,
    app_asc: f64,
    true_oblq: f64,
) -> f64 {
    let (_, X, Y) = pos_angl_aux(node, nut_long, rst, true_oblq);
    let w = X.atan2(Y);
    let sin_P = (X.hypot(Y) * (app_asc - w).cos() / total_lat.cos()).clamp(-1.0, 1.0);
    // asin gives (-π/2, π/2]; the book wants the angle in the first or
    // fourth quadrant, i.e. in [0°, 90°] ∪ [270°, 360°).
    norm_2pi(sin_P.asin())
}

/// Everything that depends only on the instant: `T`, the fundamental
/// arguments, the mean node and `ρ, σ, τ`.
struct Epoch {
    args: FundamentalArgs,
    node: f64,
    rst: RhoSigmaTau,
}

impl Epoch {
    fn new(JDE: f64) -> Epoch {
        let T = time::julian_cent(JDE);
        let args = position::fundamental_args_frm_cent(T);
        let node = position::mn_ascend_node_frm_cent(T);
        let rst = rho_sigma_tau(T, &args, node);
        Epoch { args, node, rst }
    }

    fn librations(&self, long: f64, lat: f64) -> Librations {
        let opt = optical(self.node, self.args.arg_of_lat, long, lat);
        let (l2, b2) = physical(&self.rst, opt.A, opt.lat);
        Librations {
            optical_long: opt.long,
            optical_lat: opt.lat,
            physical_long: l2,
            physical_lat: b2,
            total_long: norm_pi(opt.long + l2),
            total_lat: opt.lat + b2,
        }
    }
}

/// Computes the optical, physical and total librations of the Moon (Meeus,
/// eqs. 53.1 and 53.2).
///
/// # Arguments
///
/// * `JDE`      : Julian Ephemeris Day (TD)
/// * `moon_long`: Geocentric ecliptic longitude of the Moon *| in radians*,
///   referred to the mean equinox of the date, i.e. the apparent longitude
///   *minus* the nutation in longitude (`λ − Δψ`). This is the `long` of
///   [`position::geocentric_position`].
/// * `moon_lat` : Geocentric ecliptic latitude of the Moon *| in radians*
pub fn librations(JDE: f64, moon_long: f64, moon_lat: f64) -> Librations {
    Epoch::new(JDE).librations(moon_long, moon_lat)
}

/// Computes the position angle of the Moon's axis of rotation, including
/// the effect of the physical libration (Meeus, chapter 53).
///
/// # Returns
///
/// The position angle `P` *| in radians*, in `[0, 2π)`, measured eastward
/// from the north point of the disk.
///
/// # Arguments
///
/// * `JDE`      : Julian Ephemeris Day (TD)
/// * `moon_long`: Ecliptic longitude of the Moon without nutation
///   (`λ − Δψ`) *| in radians*
/// * `moon_lat` : Ecliptic latitude of the Moon *| in radians*
/// * `moon_asc` : Apparent right ascension of the Moon *| in radians*
/// * `nut_long` : Nutation in longitude `Δψ` *| in radians*
/// * `true_oblq`: True obliquity of the ecliptic `ε` *| in radians*
pub fn pos_angl_of_axis(
    JDE: f64,
    moon_long: f64,
    moon_lat: f64,
    moon_asc: f64,
    nut_long: f64,
    true_oblq: f64,
) -> f64 {
    let ep = Epoch::new(JDE);
    let libr = ep.librations(moon_long, moon_lat);
    pos_angl(
        ep.node,
        nut_long,
        &ep.rst,
        libr.total_lat,
        moon_asc,
        true_oblq,
    )
}

/// Heliocentric ecliptic coordinates `(λH, βH)` of the Moon, in the
/// approximation of Meeus (chapter 53). `sun_long` and `moon_long` are
/// apparent longitudes; distances in the same unit.
pub(crate) fn moon_helio_ecl(
    sun_long: f64,
    moon_long: f64,
    moon_lat: f64,
    moon_dist: f64,
    sun_dist: f64,
) -> (f64, f64) {
    let ratio = moon_dist / sun_dist;
    (
        norm_2pi(sun_long + PI + ratio * moon_lat.cos() * (sun_long - moon_long).sin()),
        ratio * moon_lat,
    )
}

impl Epoch {
    fn selenographic_sun(
        &self,
        sun_long: f64,
        moon_long: f64,
        moon_lat: f64,
        moon_dist: f64,
        sun_dist: f64,
        nut_long: f64,
    ) -> SelenographicSun {
        let (lH, bH) = moon_helio_ecl(sun_long, moon_long, moon_lat, moon_dist, sun_dist);
        let libr = self.librations(lH - nut_long, bH);
        SelenographicSun {
            long: libr.total_long,
            lat: libr.total_lat,
            colong: norm_2pi(FRAC_PI_2 - libr.total_long),
        }
    }
}

/// Computes the selenographic position of the Sun (Meeus, chapter 53).
///
/// The subsolar point is found with the libration formulae, replacing the
/// geocentric position of the Moon by its (approximate) heliocentric
/// position.
///
/// # Arguments
///
/// * `JDE`      : Julian Ephemeris Day (TD)
/// * `sun_long` : Apparent geocentric ecliptic longitude of the Sun `λ0`
///   *| in radians*
/// * `moon_long`: Apparent geocentric ecliptic longitude of the Moon `λ`
///   (with nutation) *| in radians*
/// * `moon_lat` : Geocentric ecliptic latitude of the Moon `β`
///   *| in radians*
/// * `moon_dist`: Earth–Moon distance `Δ`
/// * `sun_dist` : Earth–Sun distance `R`, in the same unit as `moon_dist`
/// * `nut_long` : Nutation in longitude `Δψ` *| in radians*
pub fn selenographic_sun(
    JDE: f64,
    sun_long: f64,
    moon_long: f64,
    moon_lat: f64,
    moon_dist: f64,
    sun_dist: f64,
    nut_long: f64,
) -> SelenographicSun {
    Epoch::new(JDE).selenographic_sun(sun_long, moon_long, moon_lat, moon_dist, sun_dist, nut_long)
}

/// Computes the physical ephemeris of the Moon for an instant: librations,
/// position angle of the axis and selenographic position of the Sun, all
/// geocentric (Meeus, chapter 53).
///
/// The Moon's position comes from [`position::apparent_position`]
/// (chapter 47) and the Sun's from VSOP87 with nutation and aberration.
///
/// # Arguments
///
/// * `JDE`: Julian Ephemeris Day (TD)
pub fn physical_ephemeris(JDE: f64) -> PhysicalEphemeris {
    let moon = position::apparent_position(JDE);
    let sun = position::apparent_sun(JDE, moon.nut_long, moon.true_oblq);
    let ep = Epoch::new(JDE);

    let libr = ep.librations(moon.long - moon.nut_long, moon.lat);
    let pos_angl_axis = pos_angl(
        ep.node,
        moon.nut_long,
        &ep.rst,
        libr.total_lat,
        moon.asc,
        moon.true_oblq,
    );
    let sun = ep.selenographic_sun(
        sun.long,
        moon.long,
        moon.lat,
        moon.dist,
        sun.dist,
        moon.nut_long,
    );

    PhysicalEphemeris {
        libr,
        pos_angl_axis,
        sun,
    }
}

/// Computes the altitude of the Sun above the local horizon of a point on
/// the Moon (Meeus, eq. 53.3).
///
/// The Sun rises at a point of selenographic longitude `η` roughly when the
/// colongitude `c0` is `−η` (mod 360°), culminates near `c0 = 90° − η` and
/// sets near `c0 = 180° − η`.
///
/// # Returns
///
/// The altitude of the centre of the Sun *| in radians*, in
/// `[-π/2, π/2]`.
///
/// # Arguments
///
/// * `sun` : Selenographic position of the Sun, from [`selenographic_sun`]
///   or [`physical_ephemeris`]
/// * `long`: Selenographic longitude `η` of the point *| in radians*,
///   positive towards Mare Crisium
/// * `lat` : Selenographic latitude `θ` of the point *| in radians*
pub fn sun_alt_at_lunar_point(sun: &SelenographicSun, long: f64, lat: f64) -> f64 {
    (sun.lat.sin() * lat.sin() + sun.lat.cos() * lat.cos() * (sun.colong + long).sin())
        .clamp(-1.0, 1.0)
        .asin()
}

/// Sunrise or sunset at a point on the Moon, iterating eq. 53.4.
fn lunar_sun_event(JDE: f64, long: f64, lat: f64, rising: bool) -> Option<f64> {
    if !(JDE.is_finite() && long.is_finite() && lat.is_finite()) {
        return None;
    }
    let cos_lat = lat.cos();
    if cos_lat < 1e-6 {
        return None;
    }

    // Start from the instant at which the colongitude reaches −η (sunrise)
    // or 180° − η (sunset), taking the nearest such instant.
    let target = if rising { -long } else { PI - long };
    let c0 = physical_ephemeris(JDE).sun.colong;
    let mut t = JDE + norm_pi(target - c0).to_degrees() / COLONG_RATE_DEG_PER_DAY;

    for _ in 0..30 {
        let sun = physical_ephemeris(t).sun;
        let h = sun_alt_at_lunar_point(&sun, long, lat).to_degrees();
        // Eq. 53.4 (Meeus uses 12.19075 °/day here)
        let step = h / (12.19075 * cos_lat);
        let dt = if rising { -step } else { step };
        if !dt.is_finite() {
            return None;
        }
        t += dt;
        // Far from the start the event is not the one asked for (near the
        // poles the Sun may stay up or down for months).
        if (t - JDE).abs() > 20.0 {
            return None;
        }
        if dt.abs() < 1e-6 {
            return Some(t);
        }
    }
    None
}

/// Computes the instant of sunrise at a point on the Moon nearest to a
/// given instant (Meeus, eqs. 53.3 and 53.4).
///
/// The first estimate is the instant at which the Sun's colongitude equals
/// `−η` (mod 360°); it is then refined by Newton-like iteration on the
/// altitude of the Sun. The result refers to the centre of the solar disk.
///
/// # Returns
///
/// The Julian Ephemeris Day of sunrise, or `None` if the point is at (or
/// very near) a pole, the Sun does not cross the horizon there near the
/// given date, or the iteration does not converge.
///
/// # Arguments
///
/// * `JDE` : Julian Ephemeris Day (TD) near the wanted sunrise. The
///   sunrise found is normally the one within half a lunation of `JDE`;
///   results more than 20 days away are rejected.
/// * `long`: Selenographic longitude `η` of the point *| in radians*,
///   positive towards Mare Crisium
/// * `lat` : Selenographic latitude `θ` of the point *| in radians*
pub fn lunar_sunrise(JDE: f64, long: f64, lat: f64) -> Option<f64> {
    lunar_sun_event(JDE, long, lat, true)
}

/// Computes the instant of sunset at a point on the Moon nearest to a
/// given instant (Meeus, eqs. 53.3 and 53.4).
///
/// See [`lunar_sunrise`]; the first estimate is the instant at which the
/// Sun's colongitude equals `180° − η`.
///
/// # Returns
///
/// The Julian Ephemeris Day of sunset, or `None` (see [`lunar_sunrise`]).
///
/// # Arguments
///
/// * `JDE` : Julian Ephemeris Day (TD) near the wanted sunset
/// * `long`: Selenographic longitude `η` of the point *| in radians*,
///   positive towards Mare Crisium
/// * `lat` : Selenographic latitude `θ` of the point *| in radians*
pub fn lunar_sunset(JDE: f64, long: f64, lat: f64) -> Option<f64> {
    lunar_sun_event(JDE, long, lat, false)
}

/// Differential corrections for the topocentric librations and position
/// angle (Meeus, chapter 53, method b), as `(Δl, Δb, ΔP)` in radians.
pub(crate) fn topocent_corrections(
    observer_lat: f64,
    geocent_dec: f64,
    loc_hr_angl: f64,
    eq_hz_parllx: f64,
    pos_angl_axis: f64,
    total_lat: f64,
) -> (f64, f64, f64) {
    let (sin_phi, cos_phi) = observer_lat.sin_cos();
    let (sin_d, cos_d) = geocent_dec.sin_cos();
    let (sin_H, cos_H) = loc_hr_angl.sin_cos();

    let Q = (cos_phi * sin_H).atan2(cos_d * sin_phi - sin_d * cos_phi * cos_H);
    let z = (sin_d * sin_phi + cos_d * cos_phi * cos_H)
        .clamp(-1.0, 1.0)
        .acos();
    let pi1 = eq_hz_parllx * (z.sin() + 0.0084 * (2.0 * z).sin());

    let dl = -pi1 * (Q - pos_angl_axis).sin() / total_lat.cos();
    let db = pi1 * (Q - pos_angl_axis).cos();
    let dP = dl * (total_lat + db).sin() - pi1 * Q.sin() * geocent_dec.tan();

    (dl, db, dP)
}

/// Greenwich apparent sidereal time for the UT instant `JDE − ΔT`.
fn apprnt_sidr_ut(JDE: f64, delta_t: f64, nut_long: f64, true_oblq: f64) -> f64 {
    let jd_ut = JDE - delta_t / 86400.0;
    time::apprnt_sidr(time::mn_sidr(jd_ut), nut_long, true_oblq)
}

/// Computes the topocentric librations and position angle of the axis by
/// direct calculation (Meeus, chapter 53, method a).
///
/// The geocentric right ascension and declination of the Moon are replaced
/// by topocentric ones (Meeus, eqs. 40.2 and 40.3, via
/// [`parallax::topocent_eq_coords`]), converted back to ecliptic
/// coordinates, and used in the libration and position-angle formulae.
/// The differences from the geocentric values can reach about 1°.
///
/// # Arguments
///
/// * `JDE`     : Julian Ephemeris Day (TD)
/// * `delta_t` : ΔT = TD − UT *| in seconds*, used to get the sidereal time
/// * `observer`: Geographic position of the observer *| in radians*, with
///   the longitude measured positive **westward** from Greenwich
/// * `height`  : Observer's height above sea level *| in metres*
pub fn topocent_librations(
    JDE: f64,
    delta_t: f64,
    observer: &coords::GeographPoint,
    height: f64,
) -> TopocentLibrations {
    let moon = position::apparent_position(JDE);
    let sidr = apprnt_sidr_ut(JDE, delta_t, moon.nut_long, moon.true_oblq);

    let geo_eq = coords::EqPoint {
        asc: moon.asc,
        dec: moon.dec,
    };
    let topo = parallax::topocent_eq_coords(&geo_eq, moon.eq_hz_parllx, observer, height, sidr);
    let long = coords::ecl_long_frm_eq(topo.asc, topo.dec, moon.true_oblq);
    let lat = coords::ecl_lat_frm_eq(topo.asc, topo.dec, moon.true_oblq);

    let ep = Epoch::new(JDE);
    let libr = ep.librations(long - moon.nut_long, lat);
    let pos_angl_axis = pos_angl(
        ep.node,
        moon.nut_long,
        &ep.rst,
        libr.total_lat,
        topo.asc,
        moon.true_oblq,
    );

    TopocentLibrations {
        long: libr.total_long,
        lat: libr.total_lat,
        pos_angl_axis,
    }
}

/// Computes the topocentric librations and position angle of the axis by
/// applying differential corrections to the geocentric values (Meeus,
/// chapter 53, method b).
///
/// The agreement with [`topocent_librations`] is about 0.005° in the
/// librations and 0.01° in the position angle.
///
/// # Arguments
///
/// * `JDE`     : Julian Ephemeris Day (TD)
/// * `delta_t` : ΔT = TD − UT *| in seconds*, used to get the sidereal time
/// * `observer`: Geographic position of the observer *| in radians*, with
///   the longitude measured positive **westward** from Greenwich (the
///   height does not enter this method)
pub fn topocent_librations_by_diff_corrections(
    JDE: f64,
    delta_t: f64,
    observer: &coords::GeographPoint,
) -> TopocentLibrations {
    let moon = position::apparent_position(JDE);
    let sidr = apprnt_sidr_ut(JDE, delta_t, moon.nut_long, moon.true_oblq);
    let H = coords::hr_angl_frm_observer_long(sidr, observer.long, moon.asc);

    let ep = Epoch::new(JDE);
    let libr = ep.librations(moon.long - moon.nut_long, moon.lat);
    let P = pos_angl(
        ep.node,
        moon.nut_long,
        &ep.rst,
        libr.total_lat,
        moon.asc,
        moon.true_oblq,
    );
    let (dl, db, dP) = topocent_corrections(
        observer.lat,
        moon.dec,
        H,
        moon.eq_hz_parllx,
        P,
        libr.total_lat,
    );

    TopocentLibrations {
        long: norm_pi(libr.total_long + dl),
        lat: libr.total_lat + db,
        pos_angl_axis: norm_2pi(P + dP),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn deg(x: f64) -> f64 {
        x.to_radians()
    }

    #[track_caller]
    fn close_deg(actual_rad: f64, expected_deg: f64, tol_deg: f64) {
        let d = (actual_rad.to_degrees() - expected_deg + 180.0).rem_euclid(360.0) - 180.0;
        assert!(
            d.abs() <= tol_deg,
            "expected {expected_deg}° ± {tol_deg}°, got {}°",
            actual_rad.to_degrees()
        );
    }

    const JDE: f64 = 2448724.5;

    // Meeus, example 53.a: intermediate values for 1992 April 12, 0h TD,
    // with the book's inputs λ − Δψ = 133.162655°, β = −3.229126°.
    #[test]
    fn example_53a_intermediate() {
        let T = time::julian_cent(JDE);
        let args = position::fundamental_args_frm_cent(T);
        let node = position::mn_ascend_node_frm_cent(T);
        close_deg(node, 274.400656, 1e-6);

        let (K1, K2) = K1_K2(T);
        assert!((K1 - 109.57).abs() < 0.005);
        assert!((K2 - 71.00).abs() < 0.005);

        let opt = optical(node, args.arg_of_lat, deg(133.162655), deg(-3.229126));
        close_deg(opt.W, 218.761999, 1e-6);
        close_deg(opt.A, 218.683932, 1e-6);
        close_deg(opt.long, -1.206, 5e-4);
        close_deg(opt.lat, 4.194, 5e-4);

        let rst = rho_sigma_tau(T, &args, node);
        close_deg(rst.rho, -0.01042, 5e-6);
        close_deg(rst.sigma, -0.01574, 5e-6);
        close_deg(rst.tau, 0.02673, 5e-6);

        let (l2, b2) = physical(&rst, opt.A, opt.lat);
        close_deg(l2, -0.025, 5e-4);
        close_deg(b2, 0.006, 5e-4);

        let nut_long = deg(0.004610);
        let eps = deg(23.440636);
        let (V, X, Y) = pos_angl_aux(node, nut_long, &rst, eps);
        // The book prints V = 273.820507°; recomputing σ from the book's own
        // D, M′, F gives 273.820473°. The 3.4e-5° gap is 9e-7° in σ, below
        // the precision of σ's coefficients, and does not affect X, Y, ω or
        // P at the printed precision. Likewise the book's I + ρ = 1.53200°
        // is the sum of the rounded I and ρ (unrounded: 1.531993°).
        close_deg(V, 273.820507, 5e-5);
        close_deg(inc_of_mn_lunar_eq() + rst.rho, 1.53200, 1e-5);
        assert!((X - -0.026676).abs() < 5e-7);
        assert!((Y - -0.396022).abs() < 5e-7);
        close_deg(X.atan2(Y), 183.8536, 5e-5);

        let P = pos_angl(node, nut_long, &rst, opt.lat + b2, deg(134.688470), eps);
        close_deg(P, 15.08, 5e-3);
    }

    // Meeus, example 53.b: the selenographic position of the Sun at
    // 1992 April 12, 0h TD, with the book's inputs.
    #[test]
    fn example_53b_intermediate() {
        let (lH, bH) = moon_helio_ecl(
            deg(22.33978),
            deg(133.167265),
            deg(-3.229126),
            368406.0,
            149971500.0,
        );
        close_deg(lH, 202.208438, 5e-6);
        close_deg(bH, -0.007932, 5e-7);

        let T = time::julian_cent(JDE);
        let args = position::fundamental_args_frm_cent(T);
        let node = position::mn_ascend_node_frm_cent(T);
        let opt = optical(node, args.arg_of_lat, lH - deg(0.004610), bH);
        close_deg(opt.W, 287.803172, 5e-6);
        close_deg(opt.A, 287.809283, 5e-6);
        close_deg(opt.long, 67.920, 5e-4);
        close_deg(opt.lat, 1.476, 5e-4);

        let rst = rho_sigma_tau(T, &args, node);
        let (l2, b2) = physical(&rst, opt.A, opt.lat);
        close_deg(l2, -0.026, 5e-4);
        close_deg(b2, -0.015, 5e-4);

        let sun = selenographic_sun(
            JDE,
            deg(22.33978),
            deg(133.167265),
            deg(-3.229126),
            368406.0,
            149971500.0,
            deg(0.004610),
        );
        close_deg(sun.long, 67.89, 5e-3);
        close_deg(sun.lat, 1.46, 5e-3);
        close_deg(sun.colong, 22.11, 5e-3);
    }

    #[test]
    fn diff_corrections_vanish_at_zenith() {
        // Moon in the zenith: H = 0, δ = φ
        let (dl, db, dP) = topocent_corrections(deg(30.0), deg(30.0), 0.0, deg(0.95), 0.3, 0.05);
        assert!(dl.abs() < 1e-12 && db.abs() < 1e-12 && dP.abs() < 1e-12);
    }

    #[test]
    fn lunar_sun_event_rejects_poles_and_nan() {
        assert_eq!(lunar_sunrise(JDE, 0.0, FRAC_PI_2), None);
        assert_eq!(lunar_sunset(f64::NAN, 0.0, 0.0), None);
    }
}

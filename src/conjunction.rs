//! Conjunctions and close approaches of two bodies (Meeus, chapters 17
//! and 18)
//!
//! The functions take the positions of each body tabulated at **three or
//! five equally spaced instants** (from an ephemeris), or a **single
//! position** for a body that may be taken as fixed, such as a star (Meeus
//! asks that a star's coordinates be entered only once). The tabulated
//! coordinates are interpolated with Meeus' formulas 3.3 and 3.8 — the
//! coordinates themselves, never derived quantities such as distances, as
//! chapter 17 insists — and the event is searched for within the tabulated
//! span.
//!
//! Instants are returned as the **interpolating factor `n`**, in units of
//! the tabular interval, measured from the central tabulated instant
//! (positive towards later instants). With 3 positions `n` is in
//! `[-1, 1]`, with 5 in `[-2, 2]`. The instant is `t_central + n * interval`.
//!
//! Positions may be equatorial ([`EqPoint`](crate::coords::EqPoint)),
//! ecliptic ([`EclPoint`](crate::coords::EclPoint)) or `(longitude,
//! latitude)` tuples in radians; see
//! [`SphericalPoint`]. Both bodies must be
//! referred to the same equinox: for a planet's apparent place, use the
//! star's apparent place too (chapter 23); for astrometric J2000 planet
//! positions, the star's J2000 position corrected for proper motion only.
//!
//! A conjunction (equal right ascension, or equal longitude) and the least
//! angular separation are different events and generally happen at
//! different instants.
//!
//! # Example
//!
//! The conjunction in right ascension of Mercury and Venus on 1991 August 7
//! (Meeus, example 18.a), from apparent positions at 0h TD on August 5 to 9:
//!
//! ```
//! use astro::{angle, conjunction};
//!
//! let p = |h, m, s, d, dm, ds| {
//!     (angle::deg_frm_hms(h, m, s).to_radians(), angle::deg_frm_dms(d, dm, ds).to_radians())
//! };
//! let mercury = [
//!     p(10, 24, 30.125, 6, 26, 32.05),
//!     p(10, 25, 0.342, 6, 10, 57.72),
//!     p(10, 25, 12.515, 5, 57, 33.08),
//!     p(10, 25, 6.235, 5, 46, 27.07),
//!     p(10, 24, 41.185, 5, 37, 48.45),
//! ];
//! let venus = [
//!     p(10, 27, 27.175, 4, 4, 41.83),
//!     p(10, 26, 32.410, 3, 55, 54.66),
//!     p(10, 25, 29.042, 3, 48, 3.51),
//!     p(10, 24, 17.191, 3, 41, 10.25),
//!     p(10, 22, 57.024, 3, 35, 16.61),
//! ];
//!
//! let c = conjunction::conjunction(&mercury, &venus).unwrap();
//!
//! // August 7.0 is the central instant and the interval is one day.
//! let day = 7.0 + c.n;
//! assert!((day - 7.23797).abs() < 1e-5); // August 7 at 5h42m TD
//! // Mercury was 2.13940° (2°08′) north of Venus.
//! assert!((c.lat_diff.to_degrees() - 2.13940).abs() < 1e-5);
//! ```

pub(crate) mod tabular;

use crate::angle::{self, SphericalPoint};
use tabular::Tabular;

/// A conjunction in right ascension (or longitude) of two bodies
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Conjunction {
    /// Instant of the conjunction, as the interpolating factor: tabular
    /// intervals from the central tabulated instant (positive later)
    pub n: f64,
    /// Latitude of body 1 minus latitude of body 2 (declination difference
    /// for equatorial coordinates) at the conjunction; positive if body 1
    /// is north of body 2 *| in radians*
    pub lat_diff: f64,
}

/// The least angular separation of two bodies
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LeastSeparation {
    /// Instant of the least separation, as the interpolating factor:
    /// tabular intervals from the central tabulated instant (positive later)
    pub n: f64,
    /// Least angular separation *| in radians*
    pub sepr: f64,
}

/// The two bodies' coordinates as interpolable series: the difference in
/// longitude (made continuous) and the two latitudes.
struct Pair {
    d_lon: Tabular,
    lat1: Tabular,
    lat2: Tabular,
    half: f64,
}

impl Pair {
    fn new<P: SphericalPoint>(body1: &[P], body2: &[P]) -> Option<Self> {
        let len = tabular::common_len(&[body1.len(), body2.len()])?;
        let mut d_lon = [0.0; 5];
        let mut lat1 = [0.0; 5];
        let mut lat2 = [0.0; 5];
        for i in 0..len {
            let p1 = tabular::nth(body1, i);
            let p2 = tabular::nth(body2, i);
            d_lon[i] = p1.lon() - p2.lon();
            lat1[i] = p1.lat();
            lat2[i] = p2.lat();
        }
        tabular::unwrap(&mut d_lon[..len]);

        Some(Self {
            d_lon: Tabular::new(&d_lon[..len])?,
            lat1: Tabular::new(&lat1[..len])?,
            lat2: Tabular::new(&lat2[..len])?,
            half: tabular::half_span(len),
        })
    }

    /// Exact angular separation of the interpolated positions at `n`
    fn sepr(&self, n: f64) -> f64 {
        angle::anglr_sepr(
            self.d_lon.value(n),
            self.lat1.value(n),
            0.0,
            self.lat2.value(n),
        )
    }

    /// Half the derivative with respect to `n` of the squared chord between
    /// the interpolated positions. Zero at a minimum of the separation.
    fn chord2_deriv(&self, n: f64) -> f64 {
        let (dl, dl_d) = (self.d_lon.value(n), self.d_lon.deriv(n));
        let (b1, b1_d) = (self.lat1.value(n), self.lat1.deriv(n));
        let (b2, b2_d) = (self.lat2.value(n), self.lat2.deriv(n));
        let (sin_b1, cos_b1) = b1.sin_cos();
        let (sin_b2, cos_b2) = b2.sin_cos();
        let (sin_dl, cos_dl) = dl.sin_cos();

        // Unit vectors p1 = (cos b1 cos dl, cos b1 sin dl, sin b1) and
        // p2 = (cos b2, 0, sin b2); w = p1 − p2, written without
        // cancellation for close positions.
        let s_half_diff = (0.5 * (b1 - b2)).sin();
        let (s_mean, c_mean) = (0.5 * (b1 + b2)).sin_cos();
        let s_half_dl = (0.5 * dl).sin();
        let w = [
            -2.0 * s_mean * s_half_diff - 2.0 * cos_b1 * s_half_dl * s_half_dl,
            cos_b1 * sin_dl,
            2.0 * c_mean * s_half_diff,
        ];
        let dw = [
            -sin_b1 * cos_dl * b1_d - cos_b1 * sin_dl * dl_d + sin_b2 * b2_d,
            -sin_b1 * sin_dl * b1_d + cos_b1 * cos_dl * dl_d,
            cos_b1 * b1_d - cos_b2 * b2_d,
        ];

        w[0] * dw[0] + w[1] * dw[1] + w[2] * dw[2]
    }
}

/// Computes the conjunction in right ascension (or longitude) of two bodies
///
/// Meeus, chapter 18. The difference in right ascension is interpolated
/// from the tabulated positions and the instant at which it vanishes is
/// found; the difference in declination is then interpolated to that
/// instant. With ecliptic coordinates this gives the conjunction in
/// longitude and the difference in latitude.
///
/// Right ascensions may be given in any range: a body crossing 0h between
/// two tabulated instants is handled.
///
/// # Returns
///
/// `Some(Conjunction)`, or `None` if there is no conjunction within the
/// tabulated span (or the input is invalid, see below). If there are
/// several, the one nearest to the central instant is returned.
///
/// # Arguments
///
/// * `body1`: Positions of body 1 at 3 or 5 equally spaced instants, or a
///   single fixed position *| in radians*
/// * `body2`: Positions of body 2, at the same instants, or a single fixed
///   position *| in radians*
///
/// The lengths must be 1, 3 or 5, and the ones that are not 1 must be
/// equal; at least one body must have 3 or 5 positions.
///
/// # Example
///
/// The conjunction of the minor planet Vesta with the star β Librae in
/// 1996 (Meeus, ch. 18, exercise). J2000 positions of Vesta at 0h TD on
/// February 7, 12, 17, 22 and 27; the star's J2000 position corrected for
/// proper motion to 1996.13.
///
/// ```
/// use astro::{angle, conjunction};
///
/// let p = |h, m, s, d, dm, ds| {
///     (angle::deg_frm_hms(h, m, s).to_radians(), angle::deg_frm_dms(d, dm, ds).to_radians())
/// };
/// let vesta = [
///     p(15, 3, 51.937, -8, 57, 34.51),
///     p(15, 9, 57.327, -9, 9, 3.88),
///     p(15, 15, 37.898, -9, 17, 37.94),
///     p(15, 20, 50.632, -9, 23, 16.25),
///     p(15, 25, 32.695, -9, 26, 1.01),
/// ];
/// let beta_lib = p(15, 17, 0.446, -9, 22, 58.47);
///
/// let c = conjunction::conjunction(&vesta, &[beta_lib]).unwrap();
///
/// // February 17.0 is the central instant; the interval is 5 days.
/// let day = 17.0 + 5.0 * c.n;
/// assert!((day - (18.0 + (6.0 + 37.0 / 60.0) / 24.0)).abs() < 0.5 / 1440.0); // Feb 18, 6h37m
/// // Vesta passed 0°03′38″ north of the star.
/// assert!((c.lat_diff.to_degrees() * 3600.0 - 218.0).abs() < 0.5);
/// ```
pub fn conjunction<P: SphericalPoint>(body1: &[P], body2: &[P]) -> Option<Conjunction> {
    let pair = Pair::new(body1, body2)?;
    let n = tabular::root_nearest_centre(|n| pair.d_lon.value(n), pair.half)?;

    Some(Conjunction {
        n,
        lat_diff: pair.lat1.value(n) - pair.lat2.value(n),
    })
}

/// Computes the least angular separation of two moving bodies
///
/// Meeus, chapter 17. Following the book, the coordinates (the two
/// latitudes and the difference in longitude) are interpolated, not the
/// separation, which is not well represented by a polynomial near a close
/// approach. The separation of the interpolated positions is then
/// minimised exactly (the book does this by trial; here the minimum is
/// located on a fine grid and refined to full precision), so the result
/// is valid for close and wide approaches alike.
///
/// # Returns
///
/// `Some(LeastSeparation)` for the smallest local minimum of the
/// separation within the tabulated span, or `None` if the separation has
/// no minimum there — for example if the bodies are still approaching at
/// the last tabulated instant — or if the input is invalid (see
/// [`conjunction`]).
///
/// # Arguments
///
/// * `body1`: Positions of body 1 at 3 or 5 equally spaced instants, or a
///   single fixed position *| in radians*
/// * `body2`: Positions of body 2, at the same instants, or a single fixed
///   position *| in radians*
///
/// # Example
///
/// Mercury and Saturn in September 1978 (Meeus, ch. 17, exercise):
/// positions at 0h TD on September 12 to 16.
///
/// ```
/// use astro::{angle, conjunction};
///
/// let p = |h, m, s, d, dm, ds| {
///     (angle::deg_frm_hms(h, m, s).to_radians(), angle::deg_frm_dms(d, dm, ds).to_radians())
/// };
/// let mercury = [
///     p(10, 23, 17.65, 11, 31, 46.3),
///     p(10, 29, 44.27, 11, 2, 5.9),
///     p(10, 36, 19.63, 10, 29, 51.7),
///     p(10, 43, 1.75, 9, 55, 16.7),
///     p(10, 49, 48.85, 9, 18, 34.7),
/// ];
/// let saturn = [
///     p(10, 33, 1.23, 10, 42, 53.5),
///     p(10, 33, 29.64, 10, 40, 13.2),
///     p(10, 33, 57.97, 10, 37, 33.4),
///     p(10, 34, 26.22, 10, 34, 53.9),
///     p(10, 34, 54.39, 10, 32, 14.9),
/// ];
///
/// let m = conjunction::least_angular_separation(&mercury, &saturn).unwrap();
///
/// // 0°03′44″ on September 13.6295 (13 at 15h06.5m TD)
/// assert!((m.sepr.to_degrees() * 3600.0 - 224.0).abs() < 0.5);
/// assert!((14.0 + m.n - 13.629498).abs() < 1e-5);
/// ```
pub fn least_angular_separation<P: SphericalPoint>(
    body1: &[P],
    body2: &[P],
) -> Option<LeastSeparation> {
    let pair = Pair::new(body1, body2)?;

    tabular::minima(|n| pair.chord2_deriv(n), pair.half)
        .into_iter()
        .map(|n| LeastSeparation {
            n,
            sepr: pair.sepr(n),
        })
        .min_by(|a, b| a.sepr.total_cmp(&b.sepr))
}

/// Danjon's rectangular coordinates `(u, v)` of body 2 relative to body 1
/// (Meeus, ch. 17), in radians rather than arcseconds.
fn danjon_uv(a1: f64, d1: f64, a2: f64, d2: f64) -> (f64, f64) {
    let d_a = angle::limit_to_pm_pi(a2 - a1);
    let d_d = d2 - d1;
    let (sin_d1, cos_d1) = d1.sin_cos();
    let tan_da = d_a.tan();
    let tan_half_da = (0.5 * d_a).tan();

    let k = 1.0 / (1.0 + sin_d1 * sin_d1 * tan_da * tan_half_da);
    let u = -k * (1.0 - d1.tan() * d_d.sin()) * cos_d1 * tan_da;
    let v = k * (d_d.sin() + sin_d1 * cos_d1 * tan_da * tan_half_da);

    (u, v)
}

/// Computes the least angular separation of two moving bodies by Danjon's
/// method of rectangular coordinates
///
/// Meeus, chapter 17, after Danjon. The position of body 2 relative to
/// body 1 is expressed in rectangular coordinates `u`, `v` at each
/// tabulated instant; these are interpolated, and the instant where
/// `u² + v²` is least is found by the iteration
/// `Δn = −(u u′ + v v′) / (u′² + v′²)` starting from `n = 0`. The least
/// separation is `√(u² + v²)`.
///
/// This reproduces the book's method; it needs no search, but `u`, `v` are
/// a projection valid only for small separations (up to a few degrees), and
/// the tangent factors fail when the bodies are 90° apart in right
/// ascension or near a pole. [`least_angular_separation`] is exact and
/// should normally be preferred. Meeus applies the method to 3 positions;
/// 5 are accepted too.
///
/// # Returns
///
/// `Some(LeastSeparation)`, or `None` if the iteration does not converge
/// within the tabulated span or the input is invalid (see [`conjunction`]).
///
/// # Arguments
///
/// * `body1`: Positions of body 1 at 3 or 5 equally spaced instants, or a
///   single fixed position *| in radians*
/// * `body2`: Positions of body 2, at the same instants, or a single fixed
///   position *| in radians*
pub fn least_angular_separation_danjon<P: SphericalPoint>(
    body1: &[P],
    body2: &[P],
) -> Option<LeastSeparation> {
    let len = tabular::common_len(&[body1.len(), body2.len()])?;
    let mut u = [0.0; 5];
    let mut v = [0.0; 5];
    for i in 0..len {
        let p1 = tabular::nth(body1, i);
        let p2 = tabular::nth(body2, i);
        (u[i], v[i]) = danjon_uv(p1.lon(), p1.lat(), p2.lon(), p2.lat());
    }
    let u = Tabular::new(&u[..len])?;
    let v = Tabular::new(&v[..len])?;
    let half = tabular::half_span(len);

    let mut n = 0.0;
    for _ in 0..100 {
        let (uu, vv) = (u.value(n), v.value(n));
        let (du, dv) = (u.deriv(n), v.deriv(n));
        let denom = du * du + dv * dv;
        // Also rejects NaN
        if denom.is_nan() || denom <= 0.0 {
            return None;
        }
        let dn = -(uu * du + vv * dv) / denom;
        n += dn;
        if !n.is_finite() || n.abs() > half + 0.5 {
            return None;
        }
        if dn.abs() < 1e-12 {
            if n.abs() > half {
                return None;
            }
            return Some(LeastSeparation {
                n,
                sepr: u.value(n).hypot(v.value(n)),
            });
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rad_hms(h: f64, m: f64, s: f64) -> f64 {
        (15.0 * (h + m / 60.0 + s / 3600.0)).to_radians()
    }

    fn rad_dms(d: f64, m: f64, s: f64) -> f64 {
        (d + m / 60.0 + s / 3600.0).to_radians()
    }

    const ARCSEC: f64 = 206_264.806_2;

    // Meeus, ch. 17: Mercury and Saturn on 1978 September 13, 14, 15
    fn mercury_saturn_uv() -> [(f64, f64); 3] {
        let mercury = [
            (rad_hms(10.0, 29.0, 44.27), rad_dms(11.0, 2.0, 5.9)),
            (rad_hms(10.0, 36.0, 19.63), rad_dms(10.0, 29.0, 51.7)),
            (rad_hms(10.0, 43.0, 1.75), rad_dms(9.0, 55.0, 16.7)),
        ];
        let saturn = [
            (rad_hms(10.0, 33.0, 29.64), rad_dms(10.0, 40.0, 13.2)),
            (rad_hms(10.0, 33.0, 57.97), rad_dms(10.0, 37.0, 33.4)),
            (rad_hms(10.0, 34.0, 26.22), rad_dms(10.0, 34.0, 53.9)),
        ];
        let mut uv = [(0.0, 0.0); 3];
        for i in 0..3 {
            let (u, v) = danjon_uv(mercury[i].0, mercury[i].1, saturn[i].0, saturn[i].1);
            uv[i] = (u * ARCSEC, v * ARCSEC);
        }
        uv
    }

    // Meeus, ch. 17: u, v in arcseconds and the first two Newton steps
    #[test]
    fn danjon_intermediate_values() {
        let uv = mercury_saturn_uv();
        let expected = [(-3322.44, -1307.48), (2088.54, 463.66), (7605.36, 2401.71)];
        for (got, exp) in uv.iter().zip(expected.iter()) {
            assert!((got.0 - exp.0).abs() < 0.01, "u {got:?} {exp:?}");
            assert!((got.1 - exp.1).abs() < 0.01, "v {got:?} {exp:?}");
        }

        let u = Tabular::new(&[uv[0].0, uv[1].0, uv[2].0]).unwrap();
        let v = Tabular::new(&[uv[0].1, uv[1].1, uv[2].1]).unwrap();
        let step = |n: f64| {
            let (du, dv) = (u.deriv(n), v.deriv(n));
            -(u.value(n) * du + v.value(n) * dv) / (du * du + dv * dv)
        };

        // n = 0: u′ = 5463.90, v′ = 1854.595, Δn = −0.368582
        assert!((u.deriv(0.0) - 5463.90).abs() < 0.01);
        assert!((v.deriv(0.0) - 1854.595).abs() < 0.01);
        let n1 = step(0.0);
        assert!((n1 - -0.368582).abs() < 1e-6, "{n1}");

        // u = 81.83, v = −208.57, u′ = 5424.89, v′ = 1793.07, Δn = −0.002142
        assert!((u.value(n1) - 81.83).abs() < 0.01);
        assert!((v.value(n1) - -208.57).abs() < 0.01);
        assert!((u.deriv(n1) - 5424.89).abs() < 0.01);
        assert!((v.deriv(n1) - 1793.07).abs() < 0.01);
        let n2 = n1 + step(n1);
        assert!((n2 - -0.370724).abs() < 1e-6, "{n2}");

        // Δn = −0.000003, n = −0.370727
        let n3 = n2 + step(n2);
        assert!((n3 - -0.370727).abs() < 1e-6, "{n3}");
    }
}

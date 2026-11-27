//! Spherical geometry of several bodies: bodies in a straight line and the
//! smallest circle containing three bodies (Meeus, chapters 19 and 20)
//!
//! A "straight line" on the celestial sphere is an arc of a great circle.
//! Positions may be equatorial ([`EqPoint`](crate::coords::EqPoint)),
//! ecliptic ([`EclPoint`](crate::coords::EclPoint)) or `(longitude,
//! latitude)` tuples in radians (see
//! [`SphericalPoint`]); all the points of one
//! call must be in the same system and, if stars are mixed with planets,
//! referred to the same equinox.
//!
//! # Example
//!
//! The three stars of Orion's Belt are not quite in a straight line
//! (Meeus, ch. 19): the angle at ε Ori is 172°29′, and ε Ori lies 324″
//! south of the great circle through δ and ζ Ori.
//!
//! ```
//! use astro::{angle, geometry};
//! use astro::angle::{Dms, Hms};
//!
//! let star = |asc: Hms, dec: Dms| (asc.to_deg().to_radians(), dec.to_deg().to_radians());
//! let delta = star(Hms::new(false, 5, 32, 0.40), Dms::new(true, 0, 17, 56.9));
//! let eps = star(Hms::new(false, 5, 36, 12.81), Dms::new(true, 1, 12, 7.0));
//! let zeta = star(Hms::new(false, 5, 40, 45.52), Dms::new(true, 1, 56, 33.3));
//!
//! let angle = geometry::vertex_angle(&delta, &eps, &zeta);
//! assert!((angle.to_degrees() - 172.4831).abs() < 1e-4);
//!
//! let omega = geometry::distance_from_great_circle(&delta, &zeta, &eps).unwrap();
//! assert!((omega.to_degrees() - -0.089876).abs() < 1e-6); // negative: south
//! ```

use crate::angle::{self, SphericalPoint};
use crate::conjunction::tabular::{self, Tabular};

type Vec3 = [f64; 3];

/// Unit vector of a point given by longitude and latitude
fn unit(lon: f64, lat: f64) -> Vec3 {
    let (sin_b, cos_b) = lat.sin_cos();
    let (sin_l, cos_l) = lon.sin_cos();
    [cos_b * cos_l, cos_b * sin_l, sin_b]
}

fn unit_of<P: SphericalPoint>(p: &P) -> Vec3 {
    unit(p.lon(), p.lat())
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: Vec3) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}

/// Computes Meeus' criterion for three bodies to be in a straight line
///
/// Meeus, ch. 19, eq. 19.1: the expression
/// `tan δ₁ sin(α₂ − α₃) + tan δ₂ sin(α₃ − α₁) + tan δ₃ sin(α₁ − α₂)`,
/// which is zero when the three bodies lie on one great circle. Its sign
/// tells on which side of the great circle through bodies 1 and 2 the
/// third body is; it changes sign when the bodies pass through alignment.
///
/// The expression equals `−[p₁, p₂, p₃] / (cos δ₁ cos δ₂ cos δ₃)`, where
/// `[p₁, p₂, p₃]` is the triple product of the unit vectors of the bodies;
/// it is undefined if a body is at a pole.
///
/// # Arguments
///
/// * `p1`, `p2`, `p3`: Positions of the three bodies *| in radians*
pub fn straight_line_criterion<P: SphericalPoint>(p1: &P, p2: &P, p3: &P) -> f64 {
    let (a1, a2, a3) = (p1.lon(), p2.lon(), p3.lon());

    p1.lat().tan() * (a2 - a3).sin()
        + p2.lat().tan() * (a3 - a1).sin()
        + p3.lat().tan() * (a1 - a2).sin()
}

/// Computes the instant when three bodies are in a straight line
///
/// Meeus, ch. 19. Each body is given by 3 or 5 positions at equally spaced
/// instants, or by a single fixed position (a star, given once). The
/// coordinates are interpolated (formulas 3.3, 3.8) and the instant at which
/// the three bodies lie on one great circle is found. (Meeus interpolates
/// the values of eq. 19.1 instead; the result is the same to the accuracy
/// of the tabulation. Here the vanishing triple product of the interpolated
/// positions is used, which also works near the poles.)
///
/// The result is for the observer for whom the positions are given: for
/// geocentric positions of the Moon, use topocentric ones (chapter 40)
/// for an observer on the Earth's surface.
///
/// # Returns
///
/// The interpolating factor `n` of the instant (tabular intervals from the
/// central tabulated instant, positive later), or `None` if the bodies are
/// not in a straight line within the tabulated span or the input is
/// invalid (lengths other than 1, 3 or 5; varying lengths that differ; all
/// three fixed). If they are aligned more than once, the instant nearest
/// the central one is returned.
///
/// # Arguments
///
/// * `body1`, `body2`, `body3`: Positions of the three bodies *| in radians*
///
/// # Example
///
/// Mars in line with Castor and Pollux in 1994 (Meeus, example 19.a): the
/// stars are fixed, Mars is given at 0h TD on September 29 to October 3.
///
/// ```
/// use astro::geometry;
///
/// let p = |a: f64, d: f64| (a.to_radians(), d.to_radians());
/// let castor = p(113.56833, 31.89756);
/// let pollux = p(116.25042, 28.03681);
/// let mars = [
///     p(118.98067, 21.68417),
///     p(119.59396, 21.58983),
///     p(120.20413, 21.49394),
///     p(120.81108, 21.39653),
///     p(121.41475, 21.29761),
/// ];
///
/// let n = geometry::time_in_straight_line(&[castor], &[pollux], &mars).unwrap();
/// // October 1.0 is the central instant: October 1.2233, at 5h TD
/// assert!((1.0 + n - 1.2233).abs() < 1e-4);
/// ```
pub fn time_in_straight_line<P: SphericalPoint>(
    body1: &[P],
    body2: &[P],
    body3: &[P],
) -> Option<f64> {
    let len = tabular::common_len(&[body1.len(), body2.len(), body3.len()])?;

    // Longitudes relative to body 2 (the triple product does not depend on
    // the origin of longitudes), made continuous across 0h.
    let mut dl1 = [0.0; 5];
    let mut dl3 = [0.0; 5];
    let mut b = [[0.0; 5]; 3];
    for i in 0..len {
        let (q1, q2, q3) = (
            tabular::nth(body1, i),
            tabular::nth(body2, i),
            tabular::nth(body3, i),
        );
        dl1[i] = q1.lon() - q2.lon();
        dl3[i] = q3.lon() - q2.lon();
        b[0][i] = q1.lat();
        b[1][i] = q2.lat();
        b[2][i] = q3.lat();
    }
    tabular::unwrap(&mut dl1[..len]);
    tabular::unwrap(&mut dl3[..len]);

    let dl1 = Tabular::new(&dl1[..len])?;
    let dl3 = Tabular::new(&dl3[..len])?;
    let b1 = Tabular::new(&b[0][..len])?;
    let b2 = Tabular::new(&b[1][..len])?;
    let b3 = Tabular::new(&b[2][..len])?;

    let triple = |n: f64| {
        let u1 = unit(dl1.value(n), b1.value(n));
        let u2 = unit(0.0, b2.value(n));
        let u3 = unit(dl3.value(n), b3.value(n));
        dot(u1, cross(u2, u3))
    };

    tabular::root_nearest_centre(triple, tabular::half_span(len))
}

/// Computes the angle at the middle one of three bodies
///
/// Meeus, ch. 19: the angle `S₁–S₂–S₃` that the arc `S₂S₁` makes with the
/// arc `S₂S₃`, which Meeus computes as the sum `C₁ + C₂` of two angles
/// measured from the meridian through `S₂`. It is swept from the direction
/// of `S₁` to that of `S₃` as seen from `S₂`, in the sense of increasing
/// position angle (from the north through the east), i.e. it is the
/// difference of the position angles of `S₃` and `S₁` relative to `S₂`.
///
/// When the bodies are taken in increasing order of right ascension
/// (`α₁ < α₂ < α₃`, over a small region of sky), this is the angle on the
/// north side of `S₂`: 180° when the three are in a straight line, less
/// than 180° when `S₂` is south of the great circle through `S₁` and `S₃`.
///
/// # Returns
///
/// * `angle`: The angle, in `[0, 2π)`. Meaningless if `S₁` or `S₃`
///   coincides with `S₂`. *| in radians*
///
/// # Arguments
///
/// * `p1`, `p2`, `p3`: Positions of `S₁`, `S₂` (the vertex) and `S₃`
///   *| in radians*
pub fn vertex_angle<P: SphericalPoint>(p1: &P, p2: &P, p3: &P) -> f64 {
    let pa1 = angle::relative_position_angle(p1.lon(), p1.lat(), p2.lon(), p2.lat());
    let pa3 = angle::relative_position_angle(p3.lon(), p3.lat(), p2.lon(), p2.lat());

    angle::limit_to_two_PI(pa3 - pa1)
}

/// Computes the angular distance of a body from the great circle through
/// two others
///
/// Meeus, ch. 19: the distance `ω` of `S₀` from the great circle through
/// `S₁` and `S₂` (extended beyond them as needed). Computed with the
/// normal to the plane of the great circle (the vector form of Meeus'
/// formula, also given by Piessens), as `atan2`, so it is accurate for all
/// distances and positions, including near the poles where Meeus' `tan α₀`
/// and `tan δ₀ / cos α₀` are unbounded.
///
/// # Returns
///
/// `Some(ω)`, in `[-π/2, π/2]` *| in radians*. Meeus' value, taken
/// between 0° and 90°, is `ω.abs()`. The sign tells the side: `ω` is
/// positive when `S₀` is on the side of the pole `p₁ × p₂` of the great
/// circle — the north side if `S₁` and `S₂` are near the equator with `S₂`
/// east of `S₁`.
///
/// `None` if `S₁` and `S₂` coincide or are diametrically opposite, so that
/// they do not define a great circle.
///
/// # Arguments
///
/// * `p1`, `p2`: Positions of the two bodies defining the great circle
///   *| in radians*
/// * `p0`: Position of the third body *| in radians*
pub fn distance_from_great_circle<P: SphericalPoint>(p1: &P, p2: &P, p0: &P) -> Option<f64> {
    let normal = cross(unit_of(p1), unit_of(p2));
    if norm(normal) == 0.0 || !norm(normal).is_finite() {
        return None;
    }
    let u0 = unit_of(p0);

    Some(dot(normal, u0).atan2(norm(cross(normal, u0))))
}

/// Computes the angle between two great circles
///
/// Meeus, ch. 19, after Piessens: the angle `ψ` between the great circle
/// through `A₁` and `A₂` and the great circle through `B₁` and `B₂`,
/// which is the angle between the normals `a₁ × a₂` and `b₁ × b₂` to their
/// planes. The circles are oriented (from `A₁` to `A₂`, from `B₁` to
/// `B₂`): reversing one of them gives `π − ψ`.
///
/// For three bodies `S₁`, `S₂`, `S₃`, the angle between the circles
/// `S₁S₂` and `S₂S₃` is the supplement of [`vertex_angle`] when the latter
/// is below 180°.
///
/// # Returns
///
/// `Some(ψ)`, in `[0, π]` *| in radians*, or `None` if either pair of
/// points does not define a great circle (coincident or antipodal points).
///
/// # Arguments
///
/// * `a1`, `a2`: Two points of the first great circle *| in radians*
/// * `b1`, `b2`: Two points of the second great circle *| in radians*
pub fn angle_between_great_circles<P: SphericalPoint>(
    a1: &P,
    a2: &P,
    b1: &P,
    b2: &P,
) -> Option<f64> {
    let na = cross(unit_of(a1), unit_of(a2));
    let nb = cross(unit_of(b1), unit_of(b2));
    if norm(na) == 0.0 || norm(nb) == 0.0 || !(norm(na) * norm(nb)).is_finite() {
        return None;
    }

    Some(norm(cross(na, nb)).atan2(dot(na, nb)))
}

/// Which of Meeus' two cases the smallest circle is
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SmallestCircleKind {
    /// Type I: the longest side of the triangle is a diameter of the
    /// circle; the third body is inside the circle.
    TypeI,
    /// Type II: the circle passes through all three bodies.
    TypeII,
}

/// The smallest circle containing three bodies
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SmallestCircle {
    /// Angular diameter of the circle *| in radians*
    pub diameter: f64,
    /// Whether it is a type I or type II circle
    pub kind: SmallestCircleKind,
}

/// Computes the smallest circle containing three bodies
///
/// Meeus, ch. 20. The sides `a ≥ b, c` of the triangle are the angular
/// separations of the bodies. If `a² ≥ b² + c²` (an obtuse or right
/// triangle) the circle is of type I and its diameter is `a`; otherwise it
/// is of type II, the circle through the three bodies, of diameter
/// `Δ = 2abc / √((a+b+c)(a+b−c)(b+c−a)(a+c−b))` (eq. 20.1).
///
/// As in the book, the triangle is treated as plane, which is adequate for
/// bodies within a few degrees of each other (Meeus suggests closer than
/// about 6°): for three planets about 4° apart the type II diameter differs
/// from that of the exact spherical circle by under 1″. Coincident bodies
/// are allowed (a type I circle).
///
/// To find the least diameter over a period, compute the diameter from
/// interpolated positions at many instants; the diameter itself is not a
/// smooth function of time and must not be interpolated (Meeus, ch. 20).
///
/// # Arguments
///
/// * `p1`, `p2`, `p3`: Positions of the three bodies *| in radians*
///
/// # Example
///
/// Mercury, Jupiter and Saturn on 1981 September 11 at 0h TD (Meeus,
/// example 20.a): a type II grouping of diameter 4.26364° (4°16′).
///
/// ```
/// use astro::{angle, geometry};
///
/// let p = |h, m, s, d, dm, ds| {
///     (angle::deg_frm_hms(h, m, s).to_radians(), angle::deg_frm_dms(d, dm, ds).to_radians())
/// };
/// let mercury = p(12, 41, 8.63, -5, 37, 54.2);
/// let jupiter = p(12, 52, 5.21, -4, 22, 26.2);
/// let saturn = p(12, 39, 28.11, -1, 50, 3.7);
///
/// let c = geometry::smallest_circle(&mercury, &jupiter, &saturn);
/// assert_eq!(c.kind, geometry::SmallestCircleKind::TypeII);
/// assert!((c.diameter.to_degrees() - 4.26364).abs() < 5e-6);
/// ```
pub fn smallest_circle<P: SphericalPoint>(p1: &P, p2: &P, p3: &P) -> SmallestCircle {
    let sepr = |p: &P, q: &P| angle::anglr_sepr(p.lon(), p.lat(), q.lon(), q.lat());
    let mut sides = [sepr(p2, p3), sepr(p1, p3), sepr(p1, p2)];
    sides.sort_by(|x, y| y.total_cmp(x));
    let [a, b, c] = sides;

    if a * a >= b * b + c * c {
        SmallestCircle {
            diameter: a,
            kind: SmallestCircleKind::TypeI,
        }
    } else {
        // Here a < b + c and b, c > 0, so every factor is positive.
        let d = 2.0 * a * b * c / ((a + b + c) * (a + b - c) * (b + c - a) * (a + c - b)).sqrt();
        SmallestCircle {
            diameter: d,
            kind: SmallestCircleKind::TypeII,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The identity relating eq. 19.1 to the triple product
    #[test]
    fn criterion_is_scaled_triple_product() {
        let pts = [
            (0.3, 0.2),
            (1.1, -0.4),
            (2.5, 0.9),
            (-1.0, 1.2),
            (4.0, -1.3),
        ];
        for &q1 in &pts {
            for &q2 in &pts {
                for &q3 in &pts {
                    let crit = straight_line_criterion(&q1, &q2, &q3);
                    let triple = dot(unit(q1.0, q1.1), cross(unit(q2.0, q2.1), unit(q3.0, q3.1)));
                    let scale = q1.1.cos() * q2.1.cos() * q3.1.cos();
                    assert!(
                        (crit + triple / scale).abs() < 1e-13,
                        "{q1:?} {q2:?} {q3:?}"
                    );
                }
            }
        }
    }
}

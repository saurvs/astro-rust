//! Solar and lunar eclipses
//!
//! Prediction of the main circumstances of eclipses of the Sun and the
//! Moon from the mean lunar phases, following Meeus, *Astronomical
//! Algorithms*, 2nd ed., chapter 54: the instant of greatest eclipse, the
//! least distance `γ` between the shadow axis and the center of the Earth
//! (solar) or of the Moon (lunar), the radius `u` of the Moon's umbral cone,
//! the type of eclipse, magnitudes and, for lunar eclipses, the
//! semidurations of the phases.
//!
//! Local circumstances of solar eclipses (which differ from place to place)
//! are not computed; these need Besselian elements.
//!
//! # Accuracy
//!
//! This is a fast, low-precision method. For the solar eclipses of
//! 1951–2050 Meeus finds a mean error of 0.36 minute (greatest 1.1 minute)
//! in the time of greatest eclipse, and a mean error of 0.00065 (greatest
//! 0.0024, about 15 km) in `γ`. In limiting cases (an eclipse barely
//! happening, or `γ` very close to the central limit) the classification
//! may be wrong and an accurate computation is needed. Lunar eclipse
//! results are adequate for historical work.
//!
//! All times are Julian Ephemeris Days (Dynamical Time, TD); subtract
//! ΔT (`time::delta_t`) to get Universal Time.
//!
//! # Example
//!
//! ```
//! use astro::eclipse::{self, SolarEclipseKind};
//!
//! // First solar eclipse after 2017 August 1, 0h TD
//! let e = eclipse::next_solar_eclipse(2457966.5).unwrap();
//!
//! // The total eclipse of 2017 August 21, greatest at about 18h27m TD
//! assert_eq!(e.kind, SolarEclipseKind::Total);
//! assert!(e.central);
//! assert!((e.jde - 2457987.2685).abs() < 0.001);
//! assert!((e.gamma - 0.4367).abs() < 0.003); // positive: northern track
//! ```

use crate::lunar::phase;

/// The node of the Moon's orbit near which an eclipse takes place: the
/// ascending node when the Moon's argument of latitude `F` is near 0°, the
/// descending node when it is near 180°
///
/// This is the same type as [`lunar::node::Node`](crate::lunar::node::Node).
pub use crate::lunar::node::Node;

/// Type of a solar eclipse
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolarEclipseKind {
    /// Only the Moon's penumbra touches the Earth
    Partial,
    /// The Moon's apparent disk is smaller than the Sun's at greatest
    /// eclipse: a ring of sunlight remains
    Annular,
    /// The Moon's apparent disk covers the Sun entirely
    Total,
    /// Annular-total (hybrid) eclipse: total along part of the track and
    /// annular along the rest. Only central eclipses are classified so.
    Hybrid,
}

/// Main circumstances of a solar eclipse
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolarEclipse {
    /// Lunation number `k` of the New Moon (an integer; 0 is the New Moon
    /// of 2000 January 6)
    pub k: f64,
    /// Instant of greatest eclipse for the Earth as a whole
    /// *| in Julian Ephemeris Days (TD)*
    pub jde: f64,
    /// Least distance from the axis of the Moon's shadow to the center of
    /// the Earth, *| in units of the Earth's equatorial radius*. Positive
    /// when the axis passes north of the Earth's center.
    pub gamma: f64,
    /// Radius of the Moon's umbral cone in the fundamental plane (the plane
    /// through the Earth's center perpendicular to the shadow axis)
    /// *| in Earth equatorial radii*. Negative when the umbral cone's
    /// vertex lies beyond the fundamental plane (a total eclipse).
    pub u: f64,
    /// Type of the eclipse
    pub kind: SolarEclipseKind,
    /// `true` if the shadow axis meets the Earth's surface, so that there
    /// is a central line (`|γ| < 0.9972`). Non-central total or annular
    /// eclipses (the umbra grazing a polar region) have `central == false`.
    pub central: bool,
    /// Node of the Moon's orbit near which the eclipse occurs
    pub node: Node,
    /// Greatest magnitude, for partial eclipses only (`None` for total,
    /// annular and hybrid eclipses). Reached at the point of the Earth's
    /// surface nearest to the shadow axis (Meeus, eq. 54.2).
    pub magnitude: Option<f64>,
}

impl SolarEclipse {
    /// Radius of the Moon's penumbral cone in the fundamental plane,
    /// `u + 0.5461` *| in Earth equatorial radii*
    pub fn penumbra_radius(&self) -> f64 {
        self.u + 0.5461
    }
}

/// Type of a lunar eclipse
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LunarEclipseKind {
    /// The Moon enters only the Earth's penumbra
    Penumbral,
    /// Part of the Moon enters the Earth's umbra
    Partial,
    /// The whole Moon enters the Earth's umbra
    Total,
}

/// Rule used to enlarge the Earth's shadow for the effect of its
/// atmosphere in lunar eclipse computations
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ShadowRule {
    /// The rule used since 1951 by the French *Connaissance des Temps*
    /// (after Danjon), Meeus' default. It better matches observations.
    #[default]
    French,
    /// The traditional rule enlarging the geometric shadow radii by 1/50.
    /// Compared with the French rule it gives magnitudes about 0.005
    /// larger for umbral and 0.026 larger for penumbral eclipses.
    Traditional,
}

/// Main circumstances of a lunar eclipse
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LunarEclipse {
    /// Lunation number `k` of the Full Moon (an integer plus 0.5)
    pub k: f64,
    /// Instant of greatest eclipse *| in Julian Ephemeris Days (TD)*
    pub jde: f64,
    /// Least distance from the center of the Moon to the axis of the
    /// Earth's shadow *| in Earth equatorial radii*. Positive when the
    /// Moon's center passes north of the axis.
    pub gamma: f64,
    /// The quantity `u` of the chapter (the Moon's umbral cone radius at
    /// the Earth), used to correct the shadow radii for the varying
    /// distances
    pub u: f64,
    /// Type of the eclipse
    pub kind: LunarEclipseKind,
    /// Node of the Moon's orbit near which the eclipse occurs
    pub node: Node,
    /// Radius of the Earth's penumbra at the Moon's distance
    /// *| in Earth equatorial radii*
    pub penumbra_radius: f64,
    /// Radius of the Earth's umbra at the Moon's distance
    /// *| in Earth equatorial radii*
    pub umbra_radius: f64,
    /// Penumbral magnitude (always positive, since otherwise there is no
    /// eclipse)
    pub penumbral_mag: f64,
    /// Umbral magnitude; zero or negative for a penumbral eclipse, 1 or more
    /// for a total eclipse
    pub umbral_mag: f64,
    /// Half the time between first and last contact with the penumbra
    /// *| in days*
    pub penumbral_semidur: f64,
    /// Half the time between first and last contact with the umbra
    /// *| in days*, or `None` for a penumbral eclipse
    pub partial_semidur: Option<f64>,
    /// Half the duration of totality *| in days*, or `None` unless the
    /// eclipse is total
    pub total_semidur: Option<f64>,
}

/// Instants of the contacts of a lunar eclipse *| in Julian Ephemeris Days
/// (TD)*
///
/// The Moon's entry into and exit from the penumbra are not visible to the
/// eye; P1 and P4 are computed quantities only.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LunarEclipseContacts {
    /// First contact with the penumbra (P1)
    pub penumbra_begin: f64,
    /// First contact with the umbra (U1), if the eclipse is umbral
    pub umbra_begin: Option<f64>,
    /// Beginning of totality (U2), if the eclipse is total
    pub total_begin: Option<f64>,
    /// Greatest eclipse
    pub greatest: f64,
    /// End of totality (U3), if the eclipse is total
    pub total_end: Option<f64>,
    /// Last contact with the umbra (U4), if the eclipse is umbral
    pub umbra_end: Option<f64>,
    /// Last contact with the penumbra (P4)
    pub penumbra_end: f64,
}

impl LunarEclipse {
    /// Computes the instants of the contacts, as the instant of greatest
    /// eclipse plus or minus the semidurations
    pub fn contacts(&self) -> LunarEclipseContacts {
        LunarEclipseContacts {
            penumbra_begin: self.jde - self.penumbral_semidur,
            umbra_begin: self.partial_semidur.map(|d| self.jde - d),
            total_begin: self.total_semidur.map(|d| self.jde - d),
            greatest: self.jde,
            total_end: self.total_semidur.map(|d| self.jde + d),
            umbra_end: self.partial_semidur.map(|d| self.jde + d),
            penumbra_end: self.jde + self.penumbral_semidur,
        }
    }
}

/// Computes the solar eclipse at a given New Moon, if there is one
///
/// # Returns
///
/// `None` if there is no eclipse at that New Moon.
///
/// # Arguments
///
/// * `lunation`: Lunation number `k` of the New Moon; 0 is the New Moon of
///   2000 January 6, negative numbers are earlier. About
///   `(year - 2000) × 12.3685` (Meeus, eq. 49.2).
pub fn solar_eclipse(lunation: i64) -> Option<SolarEclipse> {
    solar_eclipse_at(lunation as f64)
}

/// Computes the lunar eclipse at a given Full Moon, if there is one, using
/// the French shadow rule
///
/// # Returns
///
/// `None` if there is no eclipse (not even penumbral) at that Full Moon.
///
/// # Arguments
///
/// * `lunation`: Number `n` of the New Moon preceding the Full Moon; the
///   Full Moon has `k = n + 0.5` (so `k = -328.5` is `lunation = -329`).
pub fn lunar_eclipse(lunation: i64) -> Option<LunarEclipse> {
    lunar_eclipse_with_rule(lunation, ShadowRule::French)
}

/// Computes the lunar eclipse at a given Full Moon, if there is one, using
/// a chosen rule for the enlargement of the Earth's shadow
///
/// See [`lunar_eclipse`] for the arguments.
pub fn lunar_eclipse_with_rule(lunation: i64, rule: ShadowRule) -> Option<LunarEclipse> {
    lunar_eclipse_at(lunation as f64 + 0.5, rule)
}

/// Finds the first solar eclipse whose greatest phase comes after a given
/// instant
///
/// # Returns
///
/// `None` only for a non-finite or absurdly large `jde`.
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day to search from
pub fn next_solar_eclipse(jde: f64) -> Option<SolarEclipse> {
    search(jde, 0.0, true, solar_eclipse_at, |e| e.jde)
}

/// Finds the last solar eclipse whose greatest phase comes before a given
/// instant
///
/// See [`next_solar_eclipse`].
pub fn prev_solar_eclipse(jde: f64) -> Option<SolarEclipse> {
    search(jde, 0.0, false, solar_eclipse_at, |e| e.jde)
}

/// Finds the first lunar eclipse (penumbral, partial or total; French
/// shadow rule) whose greatest phase comes after a given instant
///
/// # Returns
///
/// `None` only for a non-finite or absurdly large `jde`.
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day to search from
pub fn next_lunar_eclipse(jde: f64) -> Option<LunarEclipse> {
    search(
        jde,
        0.5,
        true,
        |k| lunar_eclipse_at(k, ShadowRule::French),
        |e| e.jde,
    )
}

/// Finds the last lunar eclipse (French shadow rule) whose greatest phase
/// comes before a given instant
///
/// See [`next_lunar_eclipse`].
pub fn prev_lunar_eclipse(jde: f64) -> Option<LunarEclipse> {
    search(
        jde,
        0.5,
        false,
        |k| lunar_eclipse_at(k, ShadowRule::French),
        |e| e.jde,
    )
}

/// Steps through the lunations `k = n + frac` from just before `jde`
/// (forwards) or just after it (backwards) and returns the first eclipse
/// strictly after (or before) `jde`. Consecutive eclipses of either kind
/// are never more than six lunations apart, so a few steps suffice.
fn search<T>(
    jde: f64,
    frac: f64,
    forward: bool,
    eclipse: impl Fn(f64) -> Option<T>,
    time: impl Fn(&T) -> f64,
) -> Option<T> {
    if !jde.is_finite() {
        return None;
    }
    let n = phase::k_from_jd(jde) - frac;
    let (start, step) = if forward {
        (n.floor() - 1.0, 1.0)
    } else {
        (n.ceil() + 1.0, -1.0)
    };

    (0..24)
        .map(|i| start + step * f64::from(i) + frac)
        .filter_map(eclipse)
        .find(|e| {
            if forward {
                time(e) > jde
            } else {
                time(e) < jde
            }
        })
}

/// Quantities common to solar and lunar eclipses
struct Core {
    jde: f64,
    #[cfg(test)]
    p: f64,
    #[cfg(test)]
    q: f64,
    gamma: f64,
    u: f64,
    /// Moon's mean anomaly (radians), needed for lunar semidurations
    mp: f64,
    node: Node,
}

/// Computes the time of greatest eclipse, γ and u at lunation `k`, or
/// `None` if the Moon is too far from a node for any eclipse
/// (`|sin F| > 0.36`). Meeus, chapter 54.
fn core(k: f64, lunar: bool) -> Option<Core> {
    // Mean phase and mean arguments, as for the phases of the Moon
    // (eqs. 49.1 to 49.7 and 47.6).
    let mean = phase::mean_args(k);
    let (t, m, mp, f, om, e) = (mean.T, mean.M, mean.M1, mean.F, mean.omega, mean.E);

    if !mean.jde.is_finite() || f.sin().abs() > 0.36 {
        return None;
    }

    let f1 = f - 0.02665_f64.to_radians() * om.sin();
    let a1 = (299.77 + 0.107408 * k - 0.009173 * t * t)
        .rem_euclid(360.0)
        .to_radians();

    let (c1, c2) = if lunar {
        (-0.4065, 0.1727)
    } else {
        (-0.4075, 0.1721)
    };

    // Correction from the mean to the true instant of greatest eclipse
    // (eq. 54.1), in days.
    let dt = c1 * mp.sin() + c2 * e * m.sin() + 0.0161 * (2.0 * mp).sin()
        - 0.0097 * (2.0 * f1).sin()
        + 0.0073 * e * (mp - m).sin()
        - 0.0050 * e * (mp + m).sin()
        - 0.0023 * (mp - 2.0 * f1).sin()
        + 0.0021 * e * (2.0 * m).sin()
        + 0.0012 * (mp + 2.0 * f1).sin()
        + 0.0006 * e * (2.0 * mp + m).sin()
        - 0.0004 * (3.0 * mp).sin()
        - 0.0003 * e * (m + 2.0 * f1).sin()
        + 0.0003 * a1.sin()
        - 0.0002 * e * (m - 2.0 * f1).sin()
        - 0.0002 * e * (2.0 * mp - m).sin()
        - 0.0002 * om.sin();

    let p = 0.2070 * e * m.sin() + 0.0024 * e * (2.0 * m).sin() - 0.0392 * mp.sin()
        + 0.0116 * (2.0 * mp).sin()
        - 0.0073 * e * (mp + m).sin()
        + 0.0067 * e * (mp - m).sin()
        + 0.0118 * (2.0 * f1).sin();

    let q = 5.2207 - 0.0048 * e * m.cos() + 0.0020 * e * (2.0 * m).cos()
        - 0.3299 * mp.cos()
        - 0.0060 * e * (mp + m).cos()
        + 0.0041 * e * (mp - m).cos();

    let w = f1.cos().abs();
    let gamma = (p * f1.cos() + q * f1.sin()) * (1.0 - 0.0048 * w);

    let u = 0.0059 + 0.0046 * e * m.cos() - 0.0182 * mp.cos() + 0.0004 * (2.0 * mp).cos()
        - 0.0005 * (m + mp).cos();

    let node = if f.cos() > 0.0 {
        Node::Ascending
    } else {
        Node::Descending
    };

    Some(Core {
        jde: mean.jde + dt,
        #[cfg(test)]
        p,
        #[cfg(test)]
        q,
        gamma,
        u,
        mp,
        node,
    })
}

fn solar_eclipse_at(k: f64) -> Option<SolarEclipse> {
    let c = core(k, false)?;
    let g = c.gamma.abs();
    let u = c.u;

    // Neither the umbra nor the penumbra reaches the Earth.
    if g > 1.5433 + u {
        return None;
    }

    let central = g < 0.9972;
    let kind = if central {
        if u < 0.0 {
            SolarEclipseKind::Total
        } else if u > 0.0047 {
            SolarEclipseKind::Annular
        } else if u < 0.00464 * (1.0 - c.gamma * c.gamma).sqrt() {
            SolarEclipseKind::Hybrid
        } else {
            SolarEclipseKind::Annular
        }
    } else if g < 0.9972 + u.abs() {
        // The edge of the umbral cone touches a polar region while its axis
        // misses the Earth.
        if u < 0.0 {
            SolarEclipseKind::Total
        } else {
            SolarEclipseKind::Annular
        }
    } else {
        SolarEclipseKind::Partial
    };

    let magnitude = match kind {
        SolarEclipseKind::Partial => Some((1.5433 + u - g) / (0.5461 + 2.0 * u)),
        _ => None,
    };

    Some(SolarEclipse {
        k,
        jde: c.jde,
        gamma: c.gamma,
        u,
        kind,
        central,
        node: c.node,
        magnitude,
    })
}

/// Shadow constants (Meeus, chapter 54) for a given rule:
/// penumbra radius, umbra radius, penumbral magnitude, umbral magnitude
/// and total-phase constants.
struct ShadowConsts {
    penumbra: f64,
    umbra: f64,
    pen_mag: f64,
    umb_mag: f64,
    total: f64,
}

fn shadow_consts(rule: ShadowRule) -> ShadowConsts {
    match rule {
        ShadowRule::French => ShadowConsts {
            penumbra: 1.2848,
            umbra: 0.7403,
            pen_mag: 1.5573,
            umb_mag: 1.0128,
            total: 0.4678,
        },
        ShadowRule::Traditional => ShadowConsts {
            penumbra: 1.2985,
            umbra: 0.7432,
            pen_mag: 1.5710,
            umb_mag: 1.0157,
            total: 0.4707,
        },
    }
}

fn lunar_eclipse_at(k: f64, rule: ShadowRule) -> Option<LunarEclipse> {
    let c = core(k, true)?;
    let sc = shadow_consts(rule);
    let g = c.gamma.abs();
    let u = c.u;

    let penumbral_mag = (sc.pen_mag + u - g) / 0.5450;
    if penumbral_mag.is_nan() || penumbral_mag <= 0.0 {
        return None;
    }
    let umbral_mag = (sc.umb_mag - u - g) / 0.5450;

    // Semidurations: (60 / n) sqrt(x² - γ²) minutes, here in days.
    let n = 0.5458 + 0.0400 * c.mp.cos();
    let semidur = |x: f64| {
        if x > g {
            Some((x * x - c.gamma * c.gamma).sqrt() / n / 24.0)
        } else {
            None
        }
    };
    let p = sc.umb_mag - u;
    let t = sc.total - u;
    let h = sc.pen_mag + u;

    let partial_semidur = semidur(p);
    let total_semidur = semidur(t);
    let kind = if total_semidur.is_some() {
        LunarEclipseKind::Total
    } else if partial_semidur.is_some() {
        LunarEclipseKind::Partial
    } else {
        LunarEclipseKind::Penumbral
    };

    Some(LunarEclipse {
        k,
        jde: c.jde,
        gamma: c.gamma,
        u,
        kind,
        node: c.node,
        penumbra_radius: sc.penumbra + u,
        umbra_radius: sc.umbra - u,
        penumbral_mag,
        umbral_mag,
        penumbral_semidur: semidur(h).unwrap_or(0.0),
        partial_semidur,
        total_semidur,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() <= tol,
            "expected {expected} ± {tol}, got {actual}"
        );
    }

    fn check_mean_phase(k: f64, jde: f64, m: f64, mp: f64, f: f64) {
        let p = phase::mean_args(k);
        close(p.jde, jde, 0.00005);
        close(p.M.to_degrees(), m, 0.00005);
        close(p.M1.to_degrees(), mp, 0.00005);
        close(p.F.to_degrees(), f, 0.00005);
    }

    // Meeus, example 54.a: solar eclipse of 1993 May 21
    #[test]
    fn example_54a_intermediate_values() {
        check_mean_phase(-82.0, 2449128.5894, 135.9142, 244.5757, 165.7296);
        let p = phase::mean_args(-82.0);
        close(p.omega.to_degrees(), 253.0026, 0.00005);
        let f1 = p.F - 0.02665_f64.to_radians() * p.omega.sin();
        close(f1.to_degrees(), 165.7550, 0.00005);
        let c = core(-82.0, false).unwrap();
        // The book adds the rounded 2449128.5894 and 0.5085; unrounded the
        // sum is 2449128.58937 + 0.50846 = 2449129.09784.
        close(c.jde, 2449129.0979, 0.0001);
        close(c.gamma, 1.1348, 0.00005);
        close(c.u, 0.0097, 0.00005);
    }

    // Meeus, example 54.b: solar eclipse of 2009 July 22
    #[test]
    fn example_54b_intermediate_values() {
        check_mean_phase(118.0, 2455034.7071, 196.9855, 7.9628, 179.8301);
        let c = core(118.0, false).unwrap();
        close(c.jde, 2455034.6088, 0.00005);
        close(c.gamma, 0.0695, 0.00005);
        close(c.u, -0.0157, 0.00005);
    }

    // Meeus, example 54.c: lunar eclipse of 1973 June
    #[test]
    fn example_54c_intermediate_values() {
        check_mean_phase(-328.5, 2441849.2992, 161.4437, 180.7018, 345.4505);
        let c = core(-328.5, true).unwrap();
        close(c.jde, 2441849.3687, 0.00005);
        close(c.gamma, -1.3249, 0.00005);
        close(c.u, 0.0197, 0.00005);
    }

    // Meeus, example 54.d: first lunar eclipse after 1997 July 1
    #[test]
    fn example_54d_intermediate_values() {
        // The two preceding Full Moons are too far from a node.
        close(phase::mean_args(-30.5).F.to_degrees(), 125.2605, 0.00005);
        close(phase::mean_args(-29.5).F.to_degrees(), 155.9310, 0.00005);
        assert!(core(-30.5, true).is_none());
        assert!(core(-29.5, true).is_none());

        check_mean_phase(-28.5, 2450708.4759, 253.0507, 5.7817, 186.6015);
        let c = core(-28.5, true).unwrap();
        close(c.jde, 2450708.2835, 0.00005);
        close(c.gamma, -0.3791, 0.00005);
        close(c.u, -0.0131, 0.00005);
        let n = 0.5458 + 0.0400 * c.mp.cos();
        close(n, 0.5856, 0.00005);
    }

    // P and Q of examples 54.a and 54.b
    #[test]
    fn p_and_q() {
        for &(k, p, q) in &[(-82.0, 0.1842, 5.3589), (118.0, -0.0573, 4.9016)] {
            let c = core(k, false).unwrap();
            close(c.p, p, 0.00005);
            close(c.q, q, 0.00005);
        }
    }

    #[test]
    fn non_finite_lunations_give_none() {
        assert!(core(f64::NAN, false).is_none());
        assert!(core(f64::INFINITY, true).is_none());
        assert!(solar_eclipse_at(1e300).is_none());
    }
}

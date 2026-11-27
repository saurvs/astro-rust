//! Passages of the Moon through the nodes of its orbit (Meeus,
//! *Astronomical Algorithms*, 2nd ed., chapter 51).
//!
//! The Moon passes through a node when its geocentric ecliptic latitude is
//! zero: northward at the ascending node, southward at the descending node.
//! Meeus gives a mean instant, a polynomial in the passage number, plus
//! periodic terms.
//!
//! # Numbering
//!
//! Passages are numbered by the integer part of Meeus' `k`: `k = n` is the
//! ascending node passage number `n` and `k = n + 0.5` the descending node
//! passage that follows it. Number 0 is the ascending node passage of
//! 2000 January 21; earlier passages have negative numbers.
//!
//! # Time scale and accuracy
//!
//! Instants are Julian Ephemeris Days (Dynamical Time, TD). Meeus compared
//! 551 passages per node, both for 1980–2020 and for the years 0–40, with
//! an accurate theory: the largest error was about 2.4 minutes and most
//! errors were below one minute.
//!
//! # Example
//!
//! ```
//! use astro::lunar::node::{self, Node};
//!
//! // Meeus, example 51.a: ascending node passage of 1987 May, k = -170.
//! let jde = node::node_passage(-170, Node::Ascending);
//! assert!((jde - 2446938.76803).abs() < 1e-5); // 1987 May 23, 6h26m TD
//!
//! // The same passage, found from a date (1987 May 15.0 = JD 2446930.5).
//! let p = node::nearest_node_passage(2446930.5, Node::Ascending);
//! assert_eq!(p.number, -170);
//! ```

use super::sequence::{self, SlotIter, Timed};
use crate::time;

/// Mean draconic month in days (the coefficient of `k` in the formula for
/// the instant).
const DRACONIC_MONTH: f64 = 27.212_220_817;

/// JDE of the mean ascending node passage with `k = 0` (2000 January 21).
const EPOCH: f64 = 2_451_565.161_9;

/// A node of the Moon's orbit on the ecliptic
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Node {
    /// Ascending node: the Moon crosses the ecliptic going north
    Ascending,
    /// Descending node: the Moon crosses the ecliptic going south
    Descending,
}

impl Node {
    /// Fractional part of Meeus' `k` for this node (0 or 0.5)
    pub fn fraction(self) -> f64 {
        match self {
            Node::Ascending => 0.0,
            Node::Descending => 0.5,
        }
    }
}

/// A passage of the Moon through one of the nodes of its orbit
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NodePassage {
    /// Which node
    pub node: Node,
    /// Passage number: the integer part of Meeus' `k` (see the
    /// [module docs](self))
    pub number: i64,
    /// Instant of the passage, as a Julian Ephemeris Day (TD). NaN if the
    /// passage was requested for a non-finite Julian day.
    pub jde: f64,
}

impl NodePassage {
    /// Meeus' `k` for this passage: `number + node.fraction()`
    pub fn k(&self) -> f64 {
        self.number as f64 + self.node.fraction()
    }
}

impl Timed for NodePassage {
    fn jde(&self) -> f64 {
        self.jde
    }
}

/// Approximate (non-integer) `k` for a year with decimals (Meeus,
/// eq. 51.1), e.g. `1987.37` for mid-May 1987
///
/// To get the passage nearest to a given instant, use
/// [`nearest_node_passage`] instead.
pub fn k_from_year(year: f64) -> f64 {
    (year - 2000.05) * 13.4223
}

/// Approximate (non-integer) `k` corresponding to a Julian day: the number
/// of mean draconic months since the mean ascending node passage of 2000
/// January 21
pub fn k_from_jd(jd: f64) -> f64 {
    (jd - EPOCH) / DRACONIC_MONTH
}

/// Instant of a passage of the Moon through a node, as a Julian Ephemeris
/// Day (TD)
///
/// # Arguments
///
/// * `number`: Passage number (see the [module docs](self))
/// * `node`  : The node
pub fn node_passage(number: i64, node: Node) -> f64 {
    passage_terms(number as f64 + node.fraction()).jde
}

/// The passage through the given node that is nearest in time to `jd`
///
/// # Arguments
///
/// * `jd`  : Julian day. The difference between UT and TD is ignored when
///   choosing the nearest passage.
/// * `node`: The node wanted
///
/// If `jd` is not finite, the returned `jde` is NaN.
pub fn nearest_node_passage(jd: f64, node: Node) -> NodePassage {
    if !jd.is_finite() {
        return NodePassage {
            node,
            number: 0,
            jde: f64::NAN,
        };
    }
    sequence::nearest(jd, k_from_jd(jd) - node.fraction(), |number| NodePassage {
        node,
        number,
        jde: node_passage(number, node),
    })
}

/// Iterator over all the passages of the Moon through both nodes, in
/// chronological order, whose instants lie in `[start_jde, end_jde)`
///
/// Both bounds are Julian Ephemeris Days. The iterator is empty if
/// `start_jde` is not finite or not less than `end_jde`; an infinite
/// `end_jde` gives an unbounded iterator. Ascending and descending
/// passages alternate.
pub fn node_passages_between(start_jde: f64, end_jde: f64) -> NodePassages {
    let first_slot = 2.0 * (k_from_jd(start_jde).floor() - 1.0);
    NodePassages(SlotIter::new(
        start_jde,
        end_jde,
        first_slot,
        passage_at_slot,
    ))
}

/// Iterator returned by [`node_passages_between`]
#[derive(Debug, Clone)]
pub struct NodePassages(SlotIter<NodePassage>);

impl Iterator for NodePassages {
    type Item = NodePassage;

    fn next(&mut self) -> Option<NodePassage> {
        self.0.next()
    }
}

impl std::iter::FusedIterator for NodePassages {}

fn passage_at_slot(slot: i64) -> NodePassage {
    let number = slot.div_euclid(2);
    let node = if slot.rem_euclid(2) == 0 {
        Node::Ascending
    } else {
        Node::Descending
    };
    NodePassage {
        node,
        number,
        jde: node_passage(number, node),
    }
}

/// Computes the times of passage of the Moon through the ascending and
/// descending nodes, close to a given date
///
/// # Returns
///
/// `(time_of_ascend_node, time_of_descend_node)`
///
/// * `time_of_ascend_node` : Instant of the passage through the ascending
///   node nearest to `date`, as a Julian Ephemeris Day (TD)
/// * `time_of_descend_node`: Instant of the passage through the descending
///   node nearest to `date`, as a Julian Ephemeris Day (TD)
///
/// The two passages are chosen independently, so the descending node
/// passage may be before or after the ascending one.
///
/// # Arguments
///
/// `date`: The date
///
/// Versions before 3.0 truncated Meeus' approximate `k` towards zero and
/// always returned the descending passage that follows the ascending one;
/// they also left out the factor E on the solar terms and evaluated `T`
/// with the unrounded `k`. See [`nearest_node_passage`] for a typed API.
pub fn time_of_passage_through_nodes(date: &time::Date) -> (f64, f64) {
    let jd = time::julian_day(date);
    (
        nearest_node_passage(jd, Node::Ascending).jde,
        nearest_node_passage(jd, Node::Descending).jde,
    )
}

/// Intermediate quantities of the calculation; the unit tests compare them
/// with Meeus' example 51.a.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(not(test), allow(dead_code))]
struct PassageTerms {
    T: f64,
    /// Moon's mean elongation, radians in [0, 2π)
    D: f64,
    /// Sun's mean anomaly, radians in [0, 2π)
    M: f64,
    /// Moon's mean anomaly, radians in [0, 2π)
    M1: f64,
    /// Longitude of the Moon's mean ascending node, radians in [0, 2π)
    omega: f64,
    /// Planetary argument V, radians in [0, 2π)
    V: f64,
    /// Argument P = Ω + 272.75° − 2.3° T, radians in [0, 2π)
    P: f64,
    E: f64,
    jde: f64,
}

fn passage_terms(k: f64) -> PassageTerms {
    use sequence::reduced_rad;

    // T must be computed from the admissible (integer or half-integer) k.
    let T = k / 1342.23;
    let TT = T * T;

    let D_deg =
        183.6380 + 331.737_356_82 * k + TT * (0.001_485_2 + T * (0.000_002_09 - T * 0.000_000_010));
    let M_deg = 17.4006 + 26.820_372_50 * k + TT * (0.000_118_6 + T * 0.000_000_06);
    let M1_deg =
        38.3776 + 355.527_473_13 * k + TT * (0.012_349_9 + T * (0.000_014_627 - T * 0.000_000_069));
    let omega_deg =
        123.9767 - 1.440_989_56 * k + TT * (0.002_060_8 + T * (0.000_002_14 - T * 0.000_000_016));
    let V_deg = 299.75 + T * (132.85 - 0.009_173 * T);
    let P_deg = omega_deg + 272.75 - 2.3 * T;

    let D = reduced_rad(D_deg);
    let M = reduced_rad(M_deg);
    let M1 = reduced_rad(M1_deg);
    let omega = reduced_rad(omega_deg);
    let V = reduced_rad(V_deg);
    let P = reduced_rad(P_deg);
    let E = sequence::eccentricity_factor(T);

    let D2 = 2.0 * D;
    let D4 = 4.0 * D;

    let mean = EPOCH
        + DRACONIC_MONTH * k
        + TT * (0.000_276_2 + T * (0.000_000_021 - T * 0.000_000_000_088));

    // Terms that contain the Sun's mean anomaly M are multiplied by E.
    let periodic = -0.4721 * M1.sin() - 0.1649 * D2.sin() - 0.0868 * (D2 - M1).sin()
        + 0.0084 * (D2 + M1).sin()
        - 0.0083 * E * (D2 - M).sin()
        - 0.0039 * E * (D2 - M - M1).sin()
        + 0.0034 * (2.0 * M1).sin()
        - 0.0031 * (D2 - 2.0 * M1).sin()
        + 0.0030 * E * (D2 + M).sin()
        + 0.0028 * E * (M - M1).sin()
        + 0.0026 * E * M.sin()
        + 0.0025 * D4.sin()
        + 0.0024 * D.sin()
        + 0.0022 * E * (M + M1).sin()
        + 0.0017 * omega.sin()
        + 0.0014 * (D4 - M1).sin()
        + 0.0005 * E * (D2 + M - M1).sin()
        + 0.0004 * E * (D2 - M + M1).sin()
        - 0.0003 * E * (D2 - 2.0 * M).sin()
        + 0.0003 * E * (D4 - M).sin()
        + 0.0003 * V.sin()
        + 0.0003 * P.sin();

    PassageTerms {
        T,
        D,
        M,
        M1,
        omega,
        V,
        P,
        E,
        jde: mean + periodic,
    }
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

    // Meeus, example 51.a: ascending node passage of 1987 May, k = -170
    #[test]
    fn example_51a_intermediate_values() {
        let t = passage_terms(-170.0);
        close(t.T, -0.126655, 0.0000005);
        close(t.D.to_degrees(), 308.28736, 0.000005);
        close(t.M.to_degrees(), 137.93728, 0.000005);
        close(t.M1.to_degrees(), 78.70737, 0.000005);
        // 8.94496: the book truncates to 8.9449.
        close(t.omega.to_degrees(), 8.9449, 0.0001);
        close(t.V.to_degrees(), 282.92, 0.005);
        close(t.P.to_degrees(), 281.99, 0.005);
        close(t.E, 1.000319, 0.0000005);
        close(t.jde, 2446938.76803, 0.000005);
    }
}

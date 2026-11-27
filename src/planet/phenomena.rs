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

//! Planetary phenomena: conjunctions, oppositions, greatest elongations and
//! stations (Meeus, *Astronomical Algorithms*, 2nd ed., chapter 36).
//!
//! The times are found without computing any planetary position. A *mean*
//! configuration, which assumes circular orbits and uniform motion, recurs
//! every synodic period; its instant is `JDE₀ = A + k·B`, where the integer
//! `k` counts synodic periods from the first such configuration after
//! 2000 January 1 (`k < 0` before it). A sum of periodic terms in the mean
//! anomaly `M` of the Earth at that instant (and, for Jupiter to Neptune, in a
//! few slowly varying angles) then gives the correction to the instant of the
//! *true* configuration.
//!
//! The phenomena covered, and the planets for which the method provides
//! terms, are:
//!
//! | Phenomenon | Planets |
//! |------------|---------|
//! | [`InferiorConjunction`](Phenomenon::InferiorConjunction), [`SuperiorConjunction`](Phenomenon::SuperiorConjunction) | Mercury, Venus |
//! | [`Opposition`](Phenomenon::Opposition), [`Conjunction`](Phenomenon::Conjunction) (with the Sun) | Mars, Jupiter, Saturn, Uranus, Neptune |
//! | [`GreatestElongation`](Phenomenon::GreatestElongation) (time and value) | Mercury, Venus |
//! | [`Station`](Phenomenon::Station) in longitude | Mercury, Venus, Mars, Jupiter, Saturn |
//!
//! Greatest elongations and stations have no table row of their own: they
//! are counted from the *mean* inferior conjunction (Mercury, Venus) or the
//! *mean* opposition (Mars, Jupiter, Saturn) with the same `k`.
//!
//! # Time scale and accuracy
//!
//! All instants are Julian Ephemeris Days (dynamical time, TD). Subtract
//! ΔT to get Universal Time.
//!
//! The periodic terms are empirical fits, meant for the years
//! −2000 to +4000; do not use them outside that span. Near the ends of that
//! span the errors grow, mostly for Jupiter and Saturn. Between A.D. 1800
//! and 2200 Meeus gives the maximum error as under 1 hour for Mercury and
//! Venus, 2 hours for Saturn, Uranus and Neptune, 3 hours for Mars and 4
//! hours for Jupiter.
//!
//! Comparing every occurrence from 1800 to 2200 with instants found from
//! the apparent geocentric VSOP87 positions of
//! [`planet::heliocent_coords`](super::heliocent_coords), the crate's tests
//! measure these largest differences, in hours:
//!
//! | Planet  | Conjunctions, oppositions        | Greatest elongations | Stations |
//! |---------|----------------------------------|----------------------|----------|
//! | Mercury | 0.6                              | 0.4                  | 0.5      |
//! | Venus   | 0.4 (inferior), 1.0 (superior)   | 2.0                  | 1.1      |
//! | Mars    | 3.7 (opposition), 1.7 (conj.)    |                      | 5.6      |
//! | Jupiter | 4.0 (opposition), 4.7 (conj.)    |                      | 4.6      |
//! | Saturn  | 1.5 (opposition), 1.9 (conj.)    |                      | 4.2      |
//! | Uranus  | 2.7 (opposition), 3.6 (conj.)    |                      |          |
//! | Neptune | 2.0 (opposition), 1.3 (conj.)    |                      |          |
//!
//! A greatest elongation is a very flat maximum and a station a very flat
//! extremum of the longitude, so their instants are less well defined than
//! those of conjunctions. The elongation *values* are angular distances
//! between the planet and the centre of the Sun (not differences of
//! longitude) and agree with VSOP87 to 0.01°. Stations are stations in
//! ecliptic longitude, not in right ascension.
//!
//! # Example
//!
//! Mercury's inferior conjunction nearest to 1993 October 1 (Meeus,
//! example 36.a):
//!
//! ```
//! use astro::planet::Planet;
//! use astro::planet::phenomena::{self, Phenomenon};
//!
//! let phenomenon = Phenomenon::InferiorConjunction;
//! let k = phenomena::nearest_k(&Planet::Mercury, phenomenon, 1993.75).unwrap();
//! assert_eq!(k, -20);
//!
//! let jde = phenomena::time_of_phenomenon(&Planet::Mercury, phenomenon, k).unwrap();
//! assert!((jde - 2449297.644).abs() < 0.001); // 1993 November 6, 3h TD
//! ```

use std::f64::consts::TAU;

use crate::planet::Planet;

/// Side of the Sun on which a greatest elongation of Mercury or Venus
/// occurs
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Elongation {
    /// Greatest eastern elongation: the planet is east of the Sun and is
    /// seen in the evening sky.
    East,
    /// Greatest western elongation: the planet is west of the Sun and is
    /// seen in the morning sky.
    West,
}

/// Which of the two stations in ecliptic longitude around an inferior
/// conjunction or an opposition
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Station {
    /// The station that precedes the inferior conjunction (Mercury, Venus)
    /// or the opposition (Mars, Jupiter, Saturn): the planet stops and then
    /// starts to move westward (retrograde) among the stars.
    First,
    /// The station that follows the inferior conjunction or the opposition:
    /// the planet resumes its direct (eastward) motion.
    Second,
}

/// A planetary configuration whose time can be computed with
/// [`time_of_phenomenon`]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Phenomenon {
    /// Inferior conjunction with the Sun (Mercury and Venus): the planet
    /// passes between the Earth and the Sun.
    InferiorConjunction,
    /// Superior conjunction with the Sun (Mercury and Venus): the planet is
    /// on the far side of the Sun.
    SuperiorConjunction,
    /// Opposition (Mars to Neptune): the geocentric longitudes of the planet
    /// and the Sun differ by 180°.
    Opposition,
    /// Conjunction with the Sun (Mars to Neptune): the geocentric longitudes
    /// of the planet and the Sun are equal.
    Conjunction,
    /// Greatest elongation, east or west, of Mercury or Venus: the maximum
    /// angular distance between the planet and the Sun.
    GreatestElongation(Elongation),
    /// Station in ecliptic longitude (Mercury to Saturn): the planet's
    /// geocentric longitude is momentarily stationary.
    Station(Station),
}

/// The mean configuration from which a phenomenon is computed
///
/// For conjunctions and oppositions this is the mean conjunction or
/// opposition itself. For greatest elongations and stations it is the mean
/// inferior conjunction (Mercury, Venus) or the mean opposition (Mars,
/// Jupiter, Saturn) with the same `k`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeanPhenomenon {
    /// Number of synodic periods since the first mean configuration of this
    /// kind after 2000 January 1 (negative before it)
    pub k: i64,
    /// Instant of the mean configuration, `A + k·B` *| Julian Ephemeris Day*
    pub jde: f64,
    /// Mean anomaly of the Earth at `jde`, `M₀ + k·M₁` *| in radians, in
    /// `[0, 2π)`*
    pub earth_mean_anomaly: f64,
    /// Time of `jde` from J2000.0 *| in Julian centuries*
    pub julian_centuries: f64,
}

/// Instant and value of a greatest elongation of Mercury or Venus
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GreatestElongation {
    /// Instant of the greatest elongation *| Julian Ephemeris Day*
    pub jde: f64,
    /// Angular distance between the planet and the centre of the Sun at
    /// that instant *| in radians*
    pub elongation: f64,
}

/// Returns `true` if the method of this module provides the phenomenon for
/// the planet (see the [module documentation](self) for the list)
///
/// # Arguments
///
/// * `planet`    : Any variant of [`Planet`]
/// * `phenomenon`: The [`Phenomenon`]
///
/// # Example
///
/// ```
/// use astro::planet::Planet;
/// use astro::planet::phenomena::{self, Phenomenon, Station};
///
/// assert!(phenomena::is_available(&Planet::Mars, Phenomenon::Opposition));
/// assert!(!phenomena::is_available(&Planet::Venus, Phenomenon::Opposition));
/// assert!(!phenomena::is_available(&Planet::Uranus, Phenomenon::Station(Station::First)));
/// ```
pub fn is_available(planet: &Planet, phenomenon: Phenomenon) -> bool {
    terms(planet, phenomenon).is_some()
}

/// Computes the approximate, real-valued `k` of a phenomenon for a given
/// date
///
/// This is `(365.2425·year + 1721060 − A) / B` (Meeus, eq. 36.1), where `A`
/// and `B` belong to the configuration the phenomenon is counted from (the
/// conjunction or opposition itself, or the inferior conjunction or
/// opposition for greatest elongations and stations). Round it to the
/// nearest integer for the configuration nearest to the date, or round up
/// for the first one after it.
///
/// # Returns
///
/// `None` if the phenomenon is not available for the planet (see
/// [`is_available`]).
///
/// # Arguments
///
/// * `planet`    : Any variant of [`Planet`]
/// * `phenomenon`: The [`Phenomenon`]
/// * `year`      : Year with decimals, e.g. `1993.75` for 1993 October 1
pub fn approx_k(planet: &Planet, phenomenon: Phenomenon, year: f64) -> Option<f64> {
    let terms = terms(planet, phenomenon)?;
    Some((jde_frm_year(year) - terms.mean.a) / terms.mean.b)
}

/// Computes the integer `k` of the occurrence of a phenomenon nearest to a
/// given date
///
/// Unlike rounding [`approx_k`], this takes into account the offset of the
/// phenomenon from its mean configuration (about ±70 days for the greatest
/// elongations of Venus, for example) and the periodic correction, so that
/// [`time_of_phenomenon`] with the returned `k` is the occurrence closest in
/// time to `year`.
///
/// # Returns
///
/// `None` if the phenomenon is not available for the planet (see
/// [`is_available`]), or if `year` is not finite.
///
/// # Arguments
///
/// * `planet`    : Any variant of [`Planet`]
/// * `phenomenon`: The [`Phenomenon`]
/// * `year`      : Year with decimals, e.g. `1993.75` for 1993 October 1
pub fn nearest_k(planet: &Planet, phenomenon: Phenomenon, year: f64) -> Option<i64> {
    let terms = terms(planet, phenomenon)?;
    if !year.is_finite() {
        return None;
    }
    let target = jde_frm_year(year);
    let offset = terms.correction.constant[0];
    let guess = ((target - terms.mean.a - offset) / terms.mean.b).round() as i64;

    let mut best = guess;
    let mut best_dist = f64::INFINITY;
    for k in guess.saturating_sub(1)..=guess.saturating_add(1) {
        let dist = (time_frm_terms(&terms, k) - target).abs();
        if dist < best_dist {
            best = k;
            best_dist = dist;
        }
    }
    Some(best)
}

/// Computes the mean configuration that a phenomenon is counted from
///
/// # Returns
///
/// `None` if the phenomenon is not available for the planet (see
/// [`is_available`]).
///
/// # Arguments
///
/// * `planet`    : Any variant of [`Planet`]
/// * `phenomenon`: The [`Phenomenon`]
/// * `k`         : Integer number of synodic periods since 2000 January 1
///
/// # Example
///
/// ```
/// use astro::planet::Planet;
/// use astro::planet::phenomena::{self, Phenomenon};
///
/// // Meeus, example 36.a
/// let mean = phenomena::mean_phenomenon(&Planet::Mercury, Phenomenon::InferiorConjunction, -20)
///     .unwrap();
/// assert!((mean.jde - 2449294.473).abs() < 0.001);
/// assert!((mean.earth_mean_anomaly.to_degrees() - 299.4092).abs() < 0.0001);
/// ```
pub fn mean_phenomenon(planet: &Planet, phenomenon: Phenomenon, k: i64) -> Option<MeanPhenomenon> {
    let terms = terms(planet, phenomenon)?;
    Some(mean_frm_terms(&terms, k))
}

/// Computes the instant of a phenomenon
///
/// The instant is the mean configuration `A + k·B` plus the sum of the
/// periodic terms for the phenomenon (Meeus, tables 36.B, 36.C and 36.D).
/// For greatest elongations, [`greatest_elongation`] also returns the value
/// of the elongation.
///
/// # Returns
///
/// * `jde`: Instant of the phenomenon *| Julian Ephemeris Day*
///
/// `None` if the phenomenon is not available for the planet (see
/// [`is_available`]).
///
/// # Arguments
///
/// * `planet`    : Any variant of [`Planet`]
/// * `phenomenon`: The [`Phenomenon`]
/// * `k`         : Integer number of synodic periods since 2000 January 1,
///   for example from [`nearest_k`]. For greatest elongations and stations,
///   the `k` of the inferior conjunction (Mercury, Venus) or the opposition
///   (Mars, Jupiter, Saturn) they belong to.
///
/// # Example
///
/// ```
/// use astro::planet::Planet;
/// use astro::planet::phenomena::{self, Phenomenon, Station};
///
/// // Mars' second station after the opposition of 1997 (Meeus, example 36.d)
/// let jde =
///     phenomena::time_of_phenomenon(&Planet::Mars, Phenomenon::Station(Station::Second), -2)
///         .unwrap();
/// assert!((jde - 2450566.255).abs() < 0.001); // 1997 April 27, 18h TD
/// ```
pub fn time_of_phenomenon(planet: &Planet, phenomenon: Phenomenon, k: i64) -> Option<f64> {
    let terms = terms(planet, phenomenon)?;
    Some(time_frm_terms(&terms, k))
}

/// Computes the instant of the occurrence of a phenomenon nearest to a
/// given date
///
/// Equivalent to [`time_of_phenomenon`] with the `k` from [`nearest_k`].
///
/// # Returns
///
/// * `jde`: Instant of the phenomenon *| Julian Ephemeris Day*
///
/// `None` if the phenomenon is not available for the planet (see
/// [`is_available`]), or if `year` is not finite.
///
/// # Arguments
///
/// * `planet`    : Any variant of [`Planet`]
/// * `phenomenon`: The [`Phenomenon`]
/// * `year`      : Year with decimals, e.g. `2125.0` for the beginning of
///   2125
pub fn nearest_time_of_phenomenon(
    planet: &Planet,
    phenomenon: Phenomenon,
    year: f64,
) -> Option<f64> {
    let k = nearest_k(planet, phenomenon, year)?;
    time_of_phenomenon(planet, phenomenon, k)
}

/// Computes the instant and the value of a greatest elongation of Mercury
/// or Venus
///
/// The elongation is the angular distance between the planet and the
/// centre of the Sun (definition used by the *Astronomical Almanac*), not
/// the difference of their ecliptic longitudes.
///
/// # Returns
///
/// `None` if `planet` is not Mercury or Venus.
///
/// # Arguments
///
/// * `planet`: [`Planet::Mercury`] or [`Planet::Venus`]
/// * `side`  : Eastern (evening) or western (morning) elongation
/// * `k`     : The `k` of the inferior conjunction the elongation belongs
///   to; the eastern elongation precedes that conjunction and the western
///   one follows it.
///
/// # Example
///
/// ```
/// use astro::planet::Planet;
/// use astro::planet::phenomena::{self, Elongation};
///
/// // Meeus, example 36.c
/// let e = phenomena::greatest_elongation(&Planet::Mercury, Elongation::West, -20).unwrap();
/// assert!((e.jde - 2449314.14).abs() < 0.01); // 1993 November 22, 15h TD
/// assert!((e.elongation.to_degrees() - 19.7506).abs() < 0.0001);
/// ```
pub fn greatest_elongation(
    planet: &Planet,
    side: Elongation,
    k: i64,
) -> Option<GreatestElongation> {
    let terms = terms(planet, Phenomenon::GreatestElongation(side))?;
    let value = match (planet, side) {
        (Planet::Mercury, Elongation::East) => &MERCURY_EAST_ELONGATION_VALUE,
        (Planet::Mercury, Elongation::West) => &MERCURY_WEST_ELONGATION_VALUE,
        (Planet::Venus, Elongation::East) => &VENUS_EAST_ELONGATION_VALUE,
        (Planet::Venus, Elongation::West) => &VENUS_WEST_ELONGATION_VALUE,
        _ => return None,
    };
    let mean = mean_frm_terms(&terms, k);
    let elongation = value.eval(mean.earth_mean_anomaly, mean.julian_centuries);

    Some(GreatestElongation {
        jde: mean.jde + correction(&terms, &mean),
        elongation: elongation.to_radians(),
    })
}

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

/// Julian Ephemeris Day of a year with decimals, as used in eq. 36.1
#[inline]
fn jde_frm_year(year: f64) -> f64 {
    365.2425 * year + 1721060.0
}

/// Polynomial `c₀ + c₁·T + c₂·T²` in Julian centuries `T`
type Poly = [f64; 3];

#[inline]
fn poly(p: &Poly, t: f64) -> f64 {
    p[0] + t * (p[1] + t * p[2])
}

/// Elements of a mean configuration (one line of Meeus' table 36.A)
struct MeanElements {
    /// JDE of the mean configuration for k = 0
    a: f64,
    /// Mean synodic period, in days
    b: f64,
    /// Mean anomaly of the Earth for k = 0, in degrees
    m0: f64,
    /// Increase of that mean anomaly per synodic period, in degrees
    m1: f64,
}

/// A sum of periodic terms in the mean anomaly `M` of the Earth:
/// `constant + Σⱼ (Sⱼ·sin jM + Cⱼ·cos jM)`, every coefficient a polynomial
/// in `T`.
struct Series {
    constant: Poly,
    /// `[Sⱼ, Cⱼ]` for j = 1, 2, 3, ...
    harmonics: &'static [[Poly; 2]],
}

impl Series {
    fn eval(&self, m: f64, t: f64) -> f64 {
        let mut sum = poly(&self.constant, t);
        for (j, [s, c]) in self.harmonics.iter().enumerate() {
            let x = (j + 1) as f64 * m;
            sum += poly(s, t) * x.sin() + poly(c, t) * x.cos();
        }
        sum
    }
}

/// A periodic term in an angle other than `M` (Jupiter to Neptune):
/// `S·sin θ + C·cos θ` with `θ = θ₀ + θ₁·T` degrees.
struct ExtraTerm {
    /// `[θ₀, θ₁]`, in degrees and degrees per Julian century
    angle: [f64; 2],
    sin: Poly,
    cos: Poly,
}

impl ExtraTerm {
    fn eval(&self, t: f64) -> f64 {
        let theta = (self.angle[0] + self.angle[1] * t).to_radians();
        poly(&self.sin, t) * theta.sin() + poly(&self.cos, t) * theta.cos()
    }
}

/// Everything needed to compute one phenomenon of one planet
struct Terms {
    mean: &'static MeanElements,
    correction: &'static Series,
    extra: &'static [ExtraTerm],
}

fn terms(planet: &Planet, phenomenon: Phenomenon) -> Option<Terms> {
    use self::Station::{First, Second};
    use Elongation::{East, West};
    use Phenomenon::{
        Conjunction, GreatestElongation, InferiorConjunction, Opposition, SuperiorConjunction,
    };

    let (mean, correction, extra): (&MeanElements, &Series, &[ExtraTerm]) =
        match (planet, phenomenon) {
            (Planet::Mercury, InferiorConjunction) => (&MERCURY_INF, &MERCURY_INF_CONJ, &[]),
            (Planet::Mercury, SuperiorConjunction) => (&MERCURY_SUP, &MERCURY_SUP_CONJ, &[]),
            (Planet::Mercury, GreatestElongation(East)) => {
                (&MERCURY_INF, &MERCURY_EAST_ELONGATION, &[])
            }
            (Planet::Mercury, GreatestElongation(West)) => {
                (&MERCURY_INF, &MERCURY_WEST_ELONGATION, &[])
            }
            (Planet::Mercury, Phenomenon::Station(First)) => {
                (&MERCURY_INF, &MERCURY_STATION_1, &[])
            }
            (Planet::Mercury, Phenomenon::Station(Second)) => {
                (&MERCURY_INF, &MERCURY_STATION_2, &[])
            }

            (Planet::Venus, InferiorConjunction) => (&VENUS_INF, &VENUS_INF_CONJ, &[]),
            (Planet::Venus, SuperiorConjunction) => (&VENUS_SUP, &VENUS_SUP_CONJ, &[]),
            (Planet::Venus, GreatestElongation(East)) => (&VENUS_INF, &VENUS_EAST_ELONGATION, &[]),
            (Planet::Venus, GreatestElongation(West)) => (&VENUS_INF, &VENUS_WEST_ELONGATION, &[]),
            (Planet::Venus, Phenomenon::Station(First)) => (&VENUS_INF, &VENUS_STATION_1, &[]),
            (Planet::Venus, Phenomenon::Station(Second)) => (&VENUS_INF, &VENUS_STATION_2, &[]),

            (Planet::Mars, Opposition) => (&MARS_OPP, &MARS_OPPOSITION, &[]),
            (Planet::Mars, Conjunction) => (&MARS_CONJ, &MARS_CONJUNCTION, &[]),
            (Planet::Mars, Phenomenon::Station(First)) => (&MARS_OPP, &MARS_STATION_1, &[]),
            (Planet::Mars, Phenomenon::Station(Second)) => (&MARS_OPP, &MARS_STATION_2, &[]),

            (Planet::Jupiter, Opposition) => (&JUPITER_OPP, &JUPITER_OPPOSITION, JUPITER_EXTRA),
            (Planet::Jupiter, Conjunction) => (&JUPITER_CONJ, &JUPITER_CONJUNCTION, JUPITER_EXTRA),
            (Planet::Jupiter, Phenomenon::Station(First)) => {
                (&JUPITER_OPP, &JUPITER_STATION_1, JUPITER_EXTRA)
            }
            (Planet::Jupiter, Phenomenon::Station(Second)) => {
                (&JUPITER_OPP, &JUPITER_STATION_2, JUPITER_EXTRA)
            }

            (Planet::Saturn, Opposition) => (&SATURN_OPP, &SATURN_OPPOSITION, SATURN_EXTRA),
            (Planet::Saturn, Conjunction) => (&SATURN_CONJ, &SATURN_CONJUNCTION, SATURN_EXTRA),
            (Planet::Saturn, Phenomenon::Station(First)) => {
                (&SATURN_OPP, &SATURN_STATION_1, SATURN_EXTRA)
            }
            (Planet::Saturn, Phenomenon::Station(Second)) => {
                (&SATURN_OPP, &SATURN_STATION_2, SATURN_EXTRA)
            }

            (Planet::Uranus, Opposition) => (&URANUS_OPP, &URANUS_OPPOSITION, URANUS_EXTRA),
            (Planet::Uranus, Conjunction) => (&URANUS_CONJ, &URANUS_CONJUNCTION, URANUS_EXTRA),

            (Planet::Neptune, Opposition) => (&NEPTUNE_OPP, &NEPTUNE_OPPOSITION, NEPTUNE_EXTRA),
            (Planet::Neptune, Conjunction) => (&NEPTUNE_CONJ, &NEPTUNE_CONJUNCTION, NEPTUNE_EXTRA),

            _ => return None,
        };

    Some(Terms {
        mean,
        correction,
        extra,
    })
}

fn mean_frm_terms(terms: &Terms, k: i64) -> MeanPhenomenon {
    let k_f = k as f64;
    let jde = terms.mean.a + k_f * terms.mean.b;
    // Reduce k·M₁ before adding it, so that large |k| loses no precision.
    let m = (terms.mean.m0 + (k_f * terms.mean.m1).rem_euclid(360.0)).rem_euclid(360.0);
    let earth_mean_anomaly = m.to_radians().rem_euclid(TAU);

    MeanPhenomenon {
        k,
        jde,
        earth_mean_anomaly,
        julian_centuries: (jde - 2451545.0) / 36525.0,
    }
}

/// Sum of the periodic terms, in days
fn correction(terms: &Terms, mean: &MeanPhenomenon) -> f64 {
    let t = mean.julian_centuries;
    let mut sum = terms.correction.eval(mean.earth_mean_anomaly, t);
    for term in terms.extra {
        sum += term.eval(t);
    }
    sum
}

fn time_frm_terms(terms: &Terms, k: i64) -> f64 {
    let mean = mean_frm_terms(terms, k);
    mean.jde + correction(terms, &mean)
}

// ---------------------------------------------------------------------------
// Table 36.A: mean configurations (A, B, M₀, M₁)
// ---------------------------------------------------------------------------

const MERCURY_INF: MeanElements = MeanElements {
    a: 2451612.023,
    b: 115.8774771,
    m0: 63.5867,
    m1: 114.2088742,
};
const MERCURY_SUP: MeanElements = MeanElements {
    a: 2451554.084,
    b: 115.8774771,
    m0: 6.4822,
    m1: 114.2088742,
};
const VENUS_INF: MeanElements = MeanElements {
    a: 2451996.706,
    b: 583.921361,
    m0: 82.7311,
    m1: 215.513058,
};
const VENUS_SUP: MeanElements = MeanElements {
    a: 2451704.746,
    b: 583.921361,
    m0: 154.9745,
    m1: 215.513058,
};
const MARS_OPP: MeanElements = MeanElements {
    a: 2452097.382,
    b: 779.936104,
    m0: 181.9573,
    m1: 48.705244,
};
const MARS_CONJ: MeanElements = MeanElements {
    a: 2451707.414,
    b: 779.936104,
    m0: 157.6047,
    m1: 48.705244,
};
const JUPITER_OPP: MeanElements = MeanElements {
    a: 2451870.628,
    b: 398.884046,
    m0: 318.4681,
    m1: 33.140229,
};
const JUPITER_CONJ: MeanElements = MeanElements {
    a: 2451671.186,
    b: 398.884046,
    m0: 121.8980,
    m1: 33.140229,
};
const SATURN_OPP: MeanElements = MeanElements {
    a: 2451870.170,
    b: 378.091904,
    m0: 318.0172,
    m1: 12.647487,
};
const SATURN_CONJ: MeanElements = MeanElements {
    a: 2451681.124,
    b: 378.091904,
    m0: 131.6934,
    m1: 12.647487,
};
const URANUS_OPP: MeanElements = MeanElements {
    a: 2451764.317,
    b: 369.656035,
    m0: 213.6884,
    m1: 4.333093,
};
const URANUS_CONJ: MeanElements = MeanElements {
    a: 2451579.489,
    b: 369.656035,
    m0: 31.5219,
    m1: 4.333093,
};
const NEPTUNE_OPP: MeanElements = MeanElements {
    a: 2451753.122,
    b: 367.486703,
    m0: 202.6544,
    m1: 2.194998,
};
const NEPTUNE_CONJ: MeanElements = MeanElements {
    a: 2451569.379,
    b: 367.486703,
    m0: 21.5569,
    m1: 2.194998,
};

// ---------------------------------------------------------------------------
// Angles other than M (Jupiter to Neptune). The same terms apply to every
// phenomenon of a given planet.
// ---------------------------------------------------------------------------

/// Jupiter: angle a
const JUPITER_EXTRA: &[ExtraTerm] = &[ExtraTerm {
    angle: [82.74, 40.76],
    sin: [0.0, 0.0144, -0.00008],
    cos: [0.3642, -0.0019, -0.00029],
}];

/// Saturn: angles a, b, c, d
const SATURN_EXTRA: &[ExtraTerm] = &[
    ExtraTerm {
        angle: [82.74, 40.76],
        sin: [0.0, -0.0337, 0.00018],
        cos: [-0.8510, 0.0044, 0.00068],
    },
    ExtraTerm {
        angle: [29.86, 1181.36],
        sin: [0.0, -0.0064, 0.00004],
        cos: [0.2397, -0.0012, -0.00008],
    },
    ExtraTerm {
        angle: [14.13, 590.68],
        sin: [0.0, -0.0010, 0.0],
        cos: [0.1245, 0.0006, 0.0],
    },
    ExtraTerm {
        angle: [220.02, 1262.87],
        sin: [0.0, 0.0024, -0.00003],
        cos: [0.0477, -0.0005, -0.00006],
    },
];

/// Uranus: angles e, f
const URANUS_EXTRA: &[ExtraTerm] = &[
    ExtraTerm {
        angle: [207.83, 8.51],
        sin: [0.0; 3],
        cos: [0.8850, 0.0, 0.0],
    },
    ExtraTerm {
        angle: [108.84, 419.96],
        sin: [0.0; 3],
        cos: [0.2153, 0.0, 0.0],
    },
];

/// Neptune: angles e, g
const NEPTUNE_EXTRA: &[ExtraTerm] = &[
    ExtraTerm {
        angle: [207.83, 8.51],
        sin: [0.0; 3],
        cos: [-0.5964, 0.0, 0.0],
    },
    ExtraTerm {
        angle: [276.74, 209.98],
        sin: [0.0; 3],
        cos: [0.0728, 0.0, 0.0],
    },
];

// ---------------------------------------------------------------------------
// Table 36.B: corrections (days) to the mean conjunctions and oppositions
// ---------------------------------------------------------------------------

const MERCURY_INF_CONJ: Series = Series {
    constant: [0.0545, 0.0002, 0.0],
    harmonics: &[
        [[-6.2008, 0.0074, 0.00003], [-3.2750, -0.0197, 0.00001]],
        [[0.4737, -0.0052, -0.00001], [0.8111, 0.0033, -0.00002]],
        [[0.0037, 0.0018, 0.0], [-0.1768, 0.0, 0.00001]],
        [[-0.0211, -0.0004, 0.0], [0.0326, -0.0003, 0.0]],
        [[0.0083, 0.0001, 0.0], [-0.0040, 0.0001, 0.0]],
    ],
};

const MERCURY_SUP_CONJ: Series = Series {
    constant: [-0.0548, -0.0002, 0.0],
    harmonics: &[
        [[7.3894, -0.0100, -0.00003], [3.2200, 0.0197, -0.00001]],
        [[0.8383, -0.0064, -0.00001], [0.9666, 0.0039, -0.00003]],
        [[0.0770, -0.0026, 0.0], [0.2758, 0.0002, -0.00002]],
        [[-0.0128, -0.0008, 0.0], [0.0734, -0.0004, -0.00001]],
        [[-0.0122, -0.0002, 0.0], [0.0173, -0.0002, 0.0]],
    ],
};

const VENUS_INF_CONJ: Series = Series {
    constant: [-0.0096, 0.0002, -0.00001],
    harmonics: &[
        [[2.0009, -0.0033, -0.00001], [0.5980, -0.0104, 0.00001]],
        [[0.0967, -0.0018, -0.00003], [0.0913, 0.0009, -0.00002]],
        [[0.0046, -0.0002, 0.0], [0.0079, 0.0001, 0.0]],
    ],
};

const VENUS_SUP_CONJ: Series = Series {
    constant: [0.0099, -0.0002, -0.00001],
    harmonics: &[
        [[4.1991, -0.0121, -0.00003], [-0.6095, 0.0102, -0.00002]],
        [[0.2500, -0.0028, -0.00003], [0.0063, 0.0025, -0.00002]],
        [[0.0232, -0.0005, -0.00001], [0.0031, 0.0004, 0.0]],
    ],
};

const MARS_OPPOSITION: Series = Series {
    constant: [-0.3088, 0.0, 0.00002],
    harmonics: &[
        [[-17.6965, 0.0363, 0.00005], [18.3131, 0.0467, -0.00006]],
        [[-0.2162, -0.0198, -0.00001], [-4.5028, -0.0019, 0.00007]],
        [[0.8987, 0.0058, -0.00002], [0.7666, -0.0050, -0.00003]],
        [[-0.3636, -0.0001, 0.00002], [0.0402, 0.0032, 0.0]],
        [[0.0737, -0.0008, 0.0], [-0.0980, -0.0011, 0.0]],
    ],
};

const MARS_CONJUNCTION: Series = Series {
    constant: [0.3102, -0.0001, 0.00001],
    harmonics: &[
        [[9.7273, -0.0156, 0.00001], [-18.3195, -0.0467, 0.00009]],
        [[-1.6488, -0.0133, 0.00001], [-2.6117, -0.0020, 0.00004]],
        [[-0.6827, -0.0026, 0.00001], [0.0281, 0.0035, 0.00001]],
        [[-0.0823, 0.0006, 0.00001], [0.1584, 0.0013, 0.0]],
        [[0.0270, 0.0005, 0.0], [0.0433, 0.0, 0.0]],
    ],
};

const JUPITER_OPPOSITION: Series = Series {
    constant: [-0.1029, 0.0, -0.00009],
    harmonics: &[
        [[-1.9658, -0.0056, 0.00007], [6.1537, 0.0210, -0.00006]],
        [[-0.2081, -0.0013, 0.0], [-0.1116, -0.0010, 0.0]],
        [[0.0074, 0.0001, 0.0], [-0.0097, -0.0001, 0.0]],
    ],
};

const JUPITER_CONJUNCTION: Series = Series {
    constant: [0.1027, 0.0002, -0.00009],
    harmonics: &[
        [[-2.2637, 0.0163, -0.00003], [-6.1540, -0.0210, 0.00008]],
        [[-0.2021, -0.0017, 0.00001], [0.1310, -0.0008, 0.0]],
        [[0.0086, 0.0, 0.0], [0.0087, 0.0002, 0.0]],
    ],
};

const SATURN_OPPOSITION: Series = Series {
    constant: [-0.0209, 0.0006, 0.00023],
    harmonics: &[
        [[4.5795, -0.0312, -0.00017], [1.1462, -0.0351, 0.00011]],
        [[0.0985, -0.0015, 0.0], [0.0733, -0.0031, 0.00001]],
        [[0.0025, -0.0001, 0.0], [0.0050, -0.0002, 0.0]],
    ],
};

const SATURN_CONJUNCTION: Series = Series {
    constant: [0.0172, -0.0006, 0.00023],
    harmonics: &[
        [[-8.5885, 0.0411, 0.00020], [-1.1470, 0.0352, -0.00011]],
        [[0.3331, -0.0034, -0.00001], [0.1145, -0.0045, 0.00002]],
        [[-0.0169, 0.0002, 0.0], [-0.0109, 0.0004, 0.0]],
    ],
};

const URANUS_OPPOSITION: Series = Series {
    constant: [0.0844, -0.0006, 0.0],
    harmonics: &[
        [[-0.1048, 0.0246, 0.0], [-5.1221, 0.0104, 0.00003]],
        [[-0.1428, 0.0005, 0.0], [-0.0148, -0.0013, 0.0]],
        [[0.0, 0.0, 0.0], [0.0055, 0.0, 0.0]],
    ],
};

const URANUS_CONJUNCTION: Series = Series {
    constant: [-0.0859, 0.0003, 0.0],
    harmonics: &[
        [[-3.8179, -0.0148, 0.00003], [5.1228, -0.0105, -0.00002]],
        [[-0.0803, 0.0011, 0.0], [-0.1905, -0.0006, 0.0]],
        [[0.0088, 0.0001, 0.0], [0.0, 0.0, 0.0]],
    ],
};

const NEPTUNE_OPPOSITION: Series = Series {
    constant: [-0.0140, 0.0, 0.00001],
    harmonics: &[
        [[-1.3486, 0.0010, 0.00001], [0.8597, 0.0037, 0.0]],
        [[-0.0082, -0.0002, 0.00001], [0.0037, -0.0003, 0.0]],
    ],
};

const NEPTUNE_CONJUNCTION: Series = Series {
    constant: [0.0168, 0.0, 0.0],
    harmonics: &[
        [[-2.5606, 0.0088, 0.00002], [-0.8611, -0.0037, 0.00002]],
        [[0.0118, -0.0004, 0.00001], [0.0307, -0.0003, 0.0]],
    ],
};

// ---------------------------------------------------------------------------
// Table 36.C: greatest elongations of Mercury and Venus. Corrections (days)
// to the mean inferior conjunction, and the elongation (degrees).
// ---------------------------------------------------------------------------

const MERCURY_EAST_ELONGATION: Series = Series {
    constant: [-21.6101, 0.0002, 0.0],
    harmonics: &[
        [[-1.9803, -0.0060, 0.00001], [1.4151, -0.0072, -0.00001]],
        [[0.5528, -0.0005, -0.00001], [0.2905, 0.0034, 0.00001]],
        [[-0.1121, -0.0001, 0.00001], [-0.0098, -0.0015, 0.0]],
        [[0.0192, 0.0, 0.0], [0.0111, 0.0004, 0.0]],
        [[-0.0061, 0.0, 0.0], [-0.0032, -0.0001, 0.0]],
    ],
};

const MERCURY_EAST_ELONGATION_VALUE: Series = Series {
    constant: [22.4697, 0.0, 0.0],
    harmonics: &[
        [[-4.2666, 0.0054, 0.00002], [-1.8537, -0.0137, 0.0]],
        [[0.3598, 0.0008, -0.00001], [-0.0680, 0.0026, 0.0]],
        [[-0.0524, -0.0003, 0.0], [0.0052, -0.0006, 0.0]],
        [[0.0107, 0.0001, 0.0], [-0.0013, 0.0001, 0.0]],
        [[-0.0021, 0.0, 0.0], [0.0003, 0.0, 0.0]],
    ],
};

const MERCURY_WEST_ELONGATION: Series = Series {
    constant: [21.6249, -0.0002, 0.0],
    harmonics: &[
        [[0.1306, 0.0065, 0.0], [-2.7661, -0.0011, 0.00001]],
        [[0.2438, -0.0024, -0.00001], [0.5767, 0.0023, 0.0]],
        [[0.1041, 0.0, 0.0], [-0.0184, 0.0007, 0.0]],
        [[-0.0051, -0.0001, 0.0], [0.0048, 0.0001, 0.0]],
        [[0.0026, 0.0, 0.0], [0.0037, 0.0, 0.0]],
    ],
};

const MERCURY_WEST_ELONGATION_VALUE: Series = Series {
    constant: [22.4143, -0.0001, 0.0],
    harmonics: &[
        [[4.3651, -0.0048, -0.00002], [2.3787, 0.0121, -0.00001]],
        [[0.2674, 0.0022, 0.0], [-0.3873, 0.0008, 0.00001]],
        [[-0.0369, -0.0001, 0.0], [0.0017, -0.0001, 0.0]],
        [[0.0059, 0.0, 0.0], [0.0061, 0.0001, 0.0]],
        [[0.0007, 0.0, 0.0], [-0.0011, 0.0, 0.0]],
    ],
};

const VENUS_EAST_ELONGATION: Series = Series {
    constant: [-70.7600, 0.0002, -0.00001],
    harmonics: &[
        [[1.0282, -0.0010, -0.00001], [0.2761, -0.0060, 0.0]],
        [[-0.0438, -0.0023, 0.00002], [0.1660, -0.0037, -0.00004]],
        [[0.0036, 0.0001, 0.0], [-0.0011, 0.0, 0.00001]],
    ],
};

const VENUS_EAST_ELONGATION_VALUE: Series = Series {
    constant: [46.3173, 0.0001, 0.0],
    harmonics: &[
        [[0.6916, -0.0024, 0.0], [0.6676, -0.0045, 0.0]],
        [[0.0309, -0.0002, 0.0], [0.0036, -0.0001, 0.0]],
    ],
};

const VENUS_WEST_ELONGATION: Series = Series {
    constant: [70.7462, 0.0, -0.00001],
    harmonics: &[
        [[1.1218, -0.0025, -0.00001], [0.4538, -0.0066, 0.0]],
        [[0.1320, 0.0020, -0.00003], [-0.0702, 0.0022, 0.00004]],
        [[0.0062, -0.0001, 0.0], [0.0015, 0.0, -0.00001]],
    ],
};

const VENUS_WEST_ELONGATION_VALUE: Series = Series {
    constant: [46.3245, 0.0, 0.0],
    harmonics: &[
        [[-0.5366, -0.0003, 0.00001], [0.3097, 0.0016, -0.00001]],
        [[-0.0163, 0.0, 0.0], [-0.0075, 0.0001, 0.0]],
    ],
};

// ---------------------------------------------------------------------------
// Table 36.D: stations in longitude. Corrections (days) to the mean
// inferior conjunction (Mercury, Venus) or the mean opposition (Mars,
// Jupiter, Saturn).
// ---------------------------------------------------------------------------

const MERCURY_STATION_1: Series = Series {
    constant: [-11.0761, 0.0003, 0.0],
    harmonics: &[
        [[-4.7321, 0.0023, 0.00002], [-1.3230, -0.0156, 0.0]],
        [[0.2270, -0.0046, 0.0], [0.7184, 0.0013, -0.00002]],
        [[0.0638, 0.0016, 0.0], [-0.1655, 0.0007, 0.0]],
        [[-0.0395, -0.0003, 0.0], [0.0247, -0.0006, 0.0]],
        [[0.0131, 0.0, 0.0], [0.0008, 0.0002, 0.0]],
    ],
};

const MERCURY_STATION_2: Series = Series {
    constant: [11.1343, -0.0001, 0.0],
    harmonics: &[
        [[-3.9137, 0.0073, 0.00002], [-3.3861, -0.0128, 0.00001]],
        [[0.5222, -0.0040, -0.00002], [0.5929, 0.0039, -0.00002]],
        [[-0.0593, 0.0018, 0.0], [-0.1733, -0.0007, 0.00001]],
        [[-0.0053, -0.0006, 0.0], [0.0476, -0.0001, 0.0]],
        [[0.0070, 0.0002, 0.0], [-0.0115, 0.0001, 0.0]],
    ],
};

const VENUS_STATION_1: Series = Series {
    constant: [-21.0672, 0.0002, -0.00001],
    harmonics: &[
        [[1.9396, -0.0029, -0.00001], [1.0727, -0.0102, 0.0]],
        [[0.0404, -0.0023, -0.00001], [0.1305, -0.0004, -0.00003]],
        [[-0.0007, -0.0002, 0.0], [0.0098, 0.0, 0.0]],
    ],
};

const VENUS_STATION_2: Series = Series {
    constant: [21.0623, 0.0, -0.00001],
    harmonics: &[
        [[1.9913, -0.0040, -0.00001], [-0.0407, -0.0077, 0.0]],
        [[0.1351, -0.0009, -0.00004], [0.0303, 0.0019, 0.0]],
        [[0.0089, -0.0002, 0.0], [0.0043, 0.0001, 0.0]],
    ],
};

const MARS_STATION_1: Series = Series {
    constant: [-37.0790, -0.0009, 0.00002],
    harmonics: &[
        [[-20.0651, 0.0228, 0.00004], [14.5205, 0.0504, -0.00001]],
        [[1.1737, -0.0169, 0.0], [-4.2550, -0.0075, 0.00008]],
        [[0.4897, 0.0074, -0.00001], [1.1151, -0.0021, -0.00005]],
        [[-0.3636, -0.0020, 0.00001], [-0.1769, 0.0028, 0.00002]],
        [[0.1437, -0.0004, 0.0], [-0.0383, -0.0016, 0.0]],
    ],
};

const MARS_STATION_2: Series = Series {
    constant: [36.7191, 0.0016, 0.00003],
    harmonics: &[
        [[-12.6163, 0.0417, -0.00001], [20.1218, 0.0379, -0.00006]],
        [[-1.6360, -0.0190, 0.0], [-3.9657, 0.0045, 0.00007]],
        [[1.1546, 0.0029, -0.00003], [0.2888, -0.0073, -0.00002]],
        [[-0.3128, 0.0017, 0.00002], [0.2513, 0.0026, -0.00002]],
        [[-0.0021, -0.0016, 0.0], [-0.1497, -0.0006, 0.0]],
    ],
};

const JUPITER_STATION_1: Series = Series {
    constant: [-60.3670, -0.0001, -0.00009],
    harmonics: &[
        [[-2.3144, -0.0124, 0.00007], [6.7439, 0.0166, -0.00006]],
        [[-0.2259, -0.0010, 0.0], [-0.1497, -0.0014, 0.0]],
        [[0.0105, 0.0001, 0.0], [-0.0098, 0.0, 0.0]],
    ],
};

const JUPITER_STATION_2: Series = Series {
    constant: [60.3023, 0.0002, -0.00009],
    harmonics: &[
        [[0.3506, -0.0034, 0.00004], [5.3635, 0.0247, -0.00007]],
        [[-0.1872, -0.0016, 0.0], [-0.0037, -0.0005, 0.0]],
        [[0.0012, 0.0001, 0.0], [-0.0096, -0.0001, 0.0]],
    ],
};

const SATURN_STATION_1: Series = Series {
    constant: [-68.8840, 0.0009, 0.00023],
    harmonics: &[
        [[5.5452, -0.0279, -0.00020], [3.0727, -0.0430, 0.00007]],
        [[0.1101, -0.0006, -0.00001], [0.1654, -0.0043, 0.00001]],
        [[0.0010, 0.0001, 0.0], [0.0095, -0.0003, 0.0]],
    ],
};

const SATURN_STATION_2: Series = Series {
    constant: [68.8720, -0.0007, 0.00023],
    harmonics: &[
        [[5.9399, -0.0400, -0.00015], [-0.7998, -0.0266, 0.00014]],
        [[0.1738, -0.0032, 0.0], [-0.0039, -0.0024, 0.00001]],
        [[0.0073, -0.0002, 0.0], [0.0020, -0.0002, 0.0]],
    ],
};

#[cfg(test)]
mod tests {
    use super::*;

    const ALL_MEAN: [&MeanElements; 14] = [
        &MERCURY_INF,
        &MERCURY_SUP,
        &VENUS_INF,
        &VENUS_SUP,
        &MARS_OPP,
        &MARS_CONJ,
        &JUPITER_OPP,
        &JUPITER_CONJ,
        &SATURN_OPP,
        &SATURN_CONJ,
        &URANUS_OPP,
        &URANUS_CONJ,
        &NEPTUNE_OPP,
        &NEPTUNE_CONJ,
    ];

    /// The columns of table 36.A are not independent: M₀ is the Earth's mean
    /// anomaly at A, and M₁ is its increase over B days. Checking this
    /// guards against transcription errors in any of the four columns.
    #[test]
    fn table_36a_is_self_consistent() {
        // Mean anomaly of the Earth (Meeus, eq. 47.3, linear part)
        let earth_m = |jde: f64| 357.5291092 + 0.985600281 * (jde - 2451545.0);
        for e in ALL_MEAN {
            let m0 = earth_m(e.a).rem_euclid(360.0);
            assert!((m0 - e.m0).abs() < 0.001, "M0 {} vs {}", m0, e.m0);
            let m1 = (0.985600281 * e.b).rem_euclid(360.0);
            assert!((m1 - e.m1).abs() < 1e-5, "M1 {} vs {}", m1, e.m1);
        }
        // The two configurations of a planet are half a synodic period apart.
        for pair in ALL_MEAN.chunks(2) {
            let (first, second) = (pair[0], pair[1]);
            assert_eq!(first.b, second.b);
            assert_eq!(first.m1, second.m1);
            let half = (first.a - second.a).rem_euclid(first.b);
            assert!((half - first.b / 2.0).abs() < 0.002, "{half}");
        }
    }

    // Meeus, example 36.b: the angles a, b, c, d for Saturn at T = +1.25627
    #[test]
    fn saturn_angles_example_36b() {
        let t = 1.25627;
        let expected = [133.95, 73.97, 36.18, 6.53];
        for (term, e) in SATURN_EXTRA.iter().zip(expected) {
            let angle = (term.angle[0] + term.angle[1] * t).rem_euclid(360.0);
            assert!((angle - e).abs() < 0.006, "{angle} vs {e}");
        }
    }

    #[test]
    fn mean_anomaly_is_reduced() {
        let terms = terms(&Planet::Mercury, Phenomenon::InferiorConjunction).unwrap();
        for k in [-100_000, -20, -1, 0, 1, 20, 100_000] {
            let m = mean_frm_terms(&terms, k).earth_mean_anomaly;
            assert!((0.0..TAU).contains(&m), "{m}");
        }
    }
}

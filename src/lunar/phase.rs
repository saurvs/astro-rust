//! Phases of the Moon (Meeus, *Astronomical Algorithms*, 2nd ed.,
//! chapter 49).
//!
//! New Moon, First Quarter, Full Moon and Last Quarter are the instants at
//! which the apparent geocentric longitude of the Moon is greater than that
//! of the Sun by 0°, 90°, 180° and 270°. Meeus' method starts from the
//! *mean* phase, a polynomial in the lunation number, and adds periodic
//! corrections derived from the ELP-2000/82 lunar theory and the VSOP87
//! solar theory.
//!
//! # Numbering
//!
//! Lunations are numbered as in Meeus: lunation 0 starts with the New Moon
//! of 2000 January 6, and lunations before it have negative numbers. Meeus'
//! parameter `k` is the lunation number plus [`Phase::fraction`], e.g.
//! `k = -282.75` is the First Quarter of lunation −283. The functions here
//! take the integer lunation and a [`Phase`], so that the "meaningless"
//! values of `k` the book warns about cannot be passed.
//!
//! # Time scale and accuracy
//!
//! All instants are Julian Ephemeris Days (Dynamical Time, TD). Subtract
//! ΔT to obtain Universal Time. Meeus reports that, for 1980–2020, the
//! results differ from the full theories by about 4 seconds on average and
//! by at most about 17 seconds. The method is intended for a few millennia
//! around the year 2000.
//!
//! # Example
//!
//! ```
//! use astro::lunar::phase::{self, Phase};
//!
//! // Meeus, example 49.a: the New Moon of 1977 February is lunation -283.
//! let jde = phase::true_phase(-283, Phase::New);
//! assert!((jde - 2443192.65118).abs() < 1e-5); // 1977 Feb 18, 3h37m TD
//!
//! // The Full Moon nearest to 2024 March 1.0 (JD 2460370.5).
//! let full = phase::nearest_phase(2460370.5, Phase::Full);
//! assert_eq!(full.phase, Phase::Full);
//! assert!((full.jde - 2460365.02).abs() < 0.01); // 2024 Feb 24, 12h32m TD
//!
//! // All phases during 2024: 13 New Moons and 12 Full Moons.
//! let phases: Vec<_> = phase::phases_between(2460310.5, 2460676.5).collect();
//! assert_eq!(phases.iter().filter(|p| p.phase == Phase::New).count(), 13);
//! assert_eq!(phases.iter().filter(|p| p.phase == Phase::Full).count(), 12);
//! ```

use super::sequence::{self, SlotIter, Timed};
use crate::time;

/// Mean synodic month in days (the coefficient of `k` in eq. 49.1).
const SYNODIC_MONTH: f64 = 29.530_588_861;

/// JDE of the mean New Moon with `k = 0` (2000 January 6).
const EPOCH: f64 = 2_451_550.097_66;

/// One of the four principal phases of the Moon
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Phase {
    /// New Moon
    New,
    /// First Quarter
    First,
    /// Full Moon
    Full,
    /// Last Quarter
    Last,
}

impl Phase {
    /// The four phases in the order in which they occur in a lunation
    pub const ALL: [Phase; 4] = [Phase::New, Phase::First, Phase::Full, Phase::Last];

    /// Fraction of a lunation after New Moon at which the phase falls
    /// (0, 0.25, 0.5 or 0.75): the fractional part of Meeus' `k`.
    pub fn fraction(self) -> f64 {
        match self {
            Phase::New => 0.0,
            Phase::First => 0.25,
            Phase::Full => 0.5,
            Phase::Last => 0.75,
        }
    }

    /// Excess of the Moon's apparent geocentric ecliptic longitude over the
    /// Sun's at this phase, in radians (0, π/2, π or 3π/2)
    pub fn longitude_excess(self) -> f64 {
        self.fraction() * std::f64::consts::TAU
    }
}

/// A principal phase of the Moon and its instant
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhaseEvent {
    /// Which phase
    pub phase: Phase,
    /// Lunation number (0 is the lunation that begins with the New Moon of
    /// 2000 January 6; negative before)
    pub lunation: i64,
    /// Instant of the phase, as a Julian Ephemeris Day (TD). NaN if the
    /// event was requested for a non-finite Julian day.
    pub jde: f64,
}

impl PhaseEvent {
    /// Meeus' `k` for this event: `lunation + phase.fraction()`
    pub fn k(&self) -> f64 {
        self.lunation as f64 + self.phase.fraction()
    }
}

impl Timed for PhaseEvent {
    fn jde(&self) -> f64 {
        self.jde
    }
}

/// Approximate (non-integer) `k` for a year with decimals (Meeus,
/// eq. 49.2), e.g. `1977.13` for mid-February 1977
///
/// Meeus picks the nearest admissible `k` (an integer for New Moon, an
/// integer plus 0.25, 0.5 or 0.75 for the other phases) by hand. To get
/// the phase nearest to a given instant, use [`nearest_phase`] instead.
pub fn k_from_year(year: f64) -> f64 {
    (year - 2000.0) * 12.3685
}

/// Approximate (non-integer) `k` corresponding to a Julian day: the number
/// of mean synodic months since the mean New Moon of 2000 January 6
pub fn k_from_jd(jd: f64) -> f64 {
    (jd - EPOCH) / SYNODIC_MONTH
}

/// Instant of a *mean* phase of the Moon (Meeus, eq. 49.1), as a Julian
/// Ephemeris Day (TD)
///
/// # Arguments
///
/// * `lunation`: Lunation number (see the [module docs](self))
/// * `phase`   : The phase
pub fn mean_phase(lunation: i64, phase: Phase) -> f64 {
    let k = lunation as f64 + phase.fraction();
    mean_jde(k, k / 1236.85)
}

/// Instant of a *true* (apparent) phase of the Moon, as a Julian Ephemeris
/// Day (TD)
///
/// This is the mean phase plus all the periodic corrections of chapter 49,
/// including the 14 planetary terms.
///
/// # Arguments
///
/// * `lunation`: Lunation number (see the [module docs](self))
/// * `phase`   : The phase
pub fn true_phase(lunation: i64, phase: Phase) -> f64 {
    terms(lunation as f64 + phase.fraction(), phase).jde()
}

/// The phase of the given kind that is nearest in time to `jd`
///
/// # Arguments
///
/// * `jd`   : Julian day. The difference between UT and TD is ignored when
///   choosing the nearest phase, so for a `jd` that is almost exactly
///   between two phases either one may be returned.
/// * `phase`: The phase wanted
///
/// If `jd` is not finite, the returned `jde` is NaN.
pub fn nearest_phase(jd: f64, phase: Phase) -> PhaseEvent {
    if !jd.is_finite() {
        return PhaseEvent {
            phase,
            lunation: 0,
            jde: f64::NAN,
        };
    }
    sequence::nearest(jd, k_from_jd(jd) - phase.fraction(), |lunation| {
        PhaseEvent {
            phase,
            lunation,
            jde: true_phase(lunation, phase),
        }
    })
}

/// Iterator over all the principal phases of the Moon, in chronological
/// order, whose instants lie in `[start_jde, end_jde)`
///
/// Both bounds are Julian Ephemeris Days. The iterator is empty if
/// `start_jde` is not finite or not less than `end_jde`. An infinite
/// `end_jde` gives an unbounded iterator.
///
/// # Example
///
/// ```
/// use astro::lunar::phase::{self, Phase};
///
/// // Phases between 2000 Jan 1.0 and Feb 1.0
/// let list: Vec<_> = phase::phases_between(2451544.5, 2451575.5).collect();
/// let kinds: Vec<_> = list.iter().map(|p| p.phase).collect();
/// assert_eq!(
///     kinds,
///     [Phase::New, Phase::First, Phase::Full, Phase::Last]
/// );
/// assert_eq!(list[0].lunation, 0); // 2000 Jan 6
/// ```
pub fn phases_between(start_jde: f64, end_jde: f64) -> Phases {
    // Four slots per lunation; start one lunation early to be safe.
    let first_slot = 4.0 * (k_from_jd(start_jde).floor() - 1.0);
    Phases(SlotIter::new(start_jde, end_jde, first_slot, phase_at_slot))
}

/// Iterator returned by [`phases_between`]
#[derive(Debug, Clone)]
pub struct Phases(SlotIter<PhaseEvent>);

impl Iterator for Phases {
    type Item = PhaseEvent;

    fn next(&mut self) -> Option<PhaseEvent> {
        self.0.next()
    }
}

impl std::iter::FusedIterator for Phases {}

fn phase_at_slot(slot: i64) -> PhaseEvent {
    let lunation = slot.div_euclid(4);
    let phase = Phase::ALL[slot.rem_euclid(4) as usize];
    PhaseEvent {
        phase,
        lunation,
        jde: true_phase(lunation, phase),
    }
}

/// Computes the Julian Ephemeris Day of the phase of the Moon nearest to
/// a date
///
/// # Returns
///
/// * `JDE`: Instant of the true phase, of the kind given by `phase`, that
///   is nearest to `date`, as a Julian Ephemeris Day (TD)
///
/// # Arguments
///
/// * `date` : Date of interest (see [`nearest_phase`] for how the nearest
///   phase is chosen)
/// * `phase`: The [`Phase`]
///
/// Versions before 3.0 truncated Meeus' approximate `k` towards zero, so
/// they did not always return the nearest phase: for dates after 2000 they
/// returned the phase of the lunation that began with the mean New Moon
/// before the date, for dates before 2000 the one of the lunation after it.
pub fn time_of_phase(date: &time::Date, phase: &Phase) -> f64 {
    nearest_phase(time::julian_day(date), *phase).jde
}

/// Mean phase polynomial (eq. 49.1)
fn mean_jde(k: f64, T: f64) -> f64 {
    EPOCH + SYNODIC_MONTH * k + T * T * (0.000_154_37 + T * (-0.000_000_150 + T * 0.000_000_000_73))
}

/// The parts that make up the instant of a true phase; kept separate so
/// that the unit tests can check each one against Meeus' examples.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(not(test), allow(dead_code))]
struct Terms {
    /// Mean phase (eq. 49.1)
    mean: f64,
    /// Eccentricity factor E (eq. 47.6)
    E: f64,
    /// Sun's mean anomaly, radians in [0, 2π)
    M: f64,
    /// Moon's mean anomaly, radians in [0, 2π)
    M1: f64,
    /// Moon's argument of latitude, radians in [0, 2π)
    F: f64,
    /// Longitude of the ascending node, radians in [0, 2π)
    omega: f64,
    /// Sum of the first group of periodic terms (days)
    periodic: f64,
    /// Quarter-phase correction: +W for First, −W for Last Quarter, 0 for
    /// New and Full Moon (days)
    quarter: f64,
    /// Sum of the 14 planetary terms (days)
    planetary: f64,
}

impl Terms {
    fn jde(&self) -> f64 {
        self.mean + self.periodic + self.quarter + self.planetary
    }
}

/// Periodic terms for New and Full Moon: (coefficient for New Moon,
/// coefficient for Full Moon, power of E, multiples of M, M′, F, Ω).
const NEW_FULL_TERMS: [(f64, f64, i32, [f64; 4]); 25] = [
    (-0.40720, -0.40614, 0, [0.0, 1.0, 0.0, 0.0]),
    (0.17241, 0.17302, 1, [1.0, 0.0, 0.0, 0.0]),
    (0.01608, 0.01614, 0, [0.0, 2.0, 0.0, 0.0]),
    (0.01039, 0.01043, 0, [0.0, 0.0, 2.0, 0.0]),
    (0.00739, 0.00734, 1, [-1.0, 1.0, 0.0, 0.0]),
    (-0.00514, -0.00515, 1, [1.0, 1.0, 0.0, 0.0]),
    (0.00208, 0.00209, 2, [2.0, 0.0, 0.0, 0.0]),
    (-0.00111, -0.00111, 0, [0.0, 1.0, -2.0, 0.0]),
    (-0.00057, -0.00057, 0, [0.0, 1.0, 2.0, 0.0]),
    (0.00056, 0.00056, 1, [1.0, 2.0, 0.0, 0.0]),
    (-0.00042, -0.00042, 0, [0.0, 3.0, 0.0, 0.0]),
    (0.00042, 0.00042, 1, [1.0, 0.0, 2.0, 0.0]),
    (0.00038, 0.00038, 1, [1.0, 0.0, -2.0, 0.0]),
    (-0.00024, -0.00024, 1, [-1.0, 2.0, 0.0, 0.0]),
    (-0.00017, -0.00017, 0, [0.0, 0.0, 0.0, 1.0]),
    (-0.00007, -0.00007, 0, [2.0, 1.0, 0.0, 0.0]),
    (0.00004, 0.00004, 0, [0.0, 2.0, -2.0, 0.0]),
    (0.00004, 0.00004, 0, [3.0, 0.0, 0.0, 0.0]),
    (0.00003, 0.00003, 0, [1.0, 1.0, -2.0, 0.0]),
    (0.00003, 0.00003, 0, [0.0, 2.0, 2.0, 0.0]),
    (-0.00003, -0.00003, 0, [1.0, 1.0, 2.0, 0.0]),
    (0.00003, 0.00003, 0, [-1.0, 1.0, 2.0, 0.0]),
    (-0.00002, -0.00002, 0, [-1.0, 1.0, -2.0, 0.0]),
    (-0.00002, -0.00002, 0, [1.0, 3.0, 0.0, 0.0]),
    (0.00002, 0.00002, 0, [0.0, 4.0, 0.0, 0.0]),
];

/// Periodic terms for First and Last Quarter: (coefficient, power of E,
/// multiples of M, M′, F, Ω).
const QUARTER_TERMS: [(f64, i32, [f64; 4]); 25] = [
    (-0.62801, 0, [0.0, 1.0, 0.0, 0.0]),
    (0.17172, 1, [1.0, 0.0, 0.0, 0.0]),
    (-0.01183, 1, [1.0, 1.0, 0.0, 0.0]),
    (0.00862, 0, [0.0, 2.0, 0.0, 0.0]),
    (0.00804, 0, [0.0, 0.0, 2.0, 0.0]),
    (0.00454, 1, [-1.0, 1.0, 0.0, 0.0]),
    (0.00204, 2, [2.0, 0.0, 0.0, 0.0]),
    (-0.00180, 0, [0.0, 1.0, -2.0, 0.0]),
    (-0.00070, 0, [0.0, 1.0, 2.0, 0.0]),
    (-0.00040, 0, [0.0, 3.0, 0.0, 0.0]),
    // Versions before 3.0 omitted the factor E on this term.
    (-0.00034, 1, [-1.0, 2.0, 0.0, 0.0]),
    (0.00032, 1, [1.0, 0.0, 2.0, 0.0]),
    (0.00032, 1, [1.0, 0.0, -2.0, 0.0]),
    (-0.00028, 2, [2.0, 1.0, 0.0, 0.0]),
    (0.00027, 1, [1.0, 2.0, 0.0, 0.0]),
    (-0.00017, 0, [0.0, 0.0, 0.0, 1.0]),
    (-0.00005, 0, [-1.0, 1.0, -2.0, 0.0]),
    (0.00004, 0, [0.0, 2.0, 2.0, 0.0]),
    (-0.00004, 0, [1.0, 1.0, 2.0, 0.0]),
    (0.00004, 0, [-2.0, 1.0, 0.0, 0.0]),
    (0.00003, 0, [1.0, 1.0, -2.0, 0.0]),
    (0.00003, 0, [3.0, 0.0, 0.0, 0.0]),
    (0.00002, 0, [0.0, 2.0, -2.0, 0.0]),
    (0.00002, 0, [-1.0, 1.0, 2.0, 0.0]),
    (-0.00002, 0, [1.0, 3.0, 0.0, 0.0]),
];

/// The 14 planetary terms A1..A14: (coefficient in days, constant and rate
/// per lunation of the argument in degrees). A1 also has a T² term.
const PLANETARY_TERMS: [(f64, f64, f64); 14] = [
    (0.000_325, 299.77, 0.107_408),
    (0.000_165, 251.88, 0.016_321),
    (0.000_164, 251.83, 26.651_886),
    (0.000_126, 349.42, 36.412_478),
    (0.000_110, 84.66, 18.206_239),
    (0.000_062, 141.74, 53.303_771),
    (0.000_060, 207.14, 2.453_732),
    (0.000_056, 154.84, 7.306_860),
    (0.000_047, 34.52, 27.261_239),
    (0.000_042, 207.19, 0.121_824),
    (0.000_040, 291.34, 1.844_379),
    (0.000_037, 161.72, 24.198_154),
    (0.000_035, 239.56, 25.513_099),
    (0.000_023, 331.55, 3.592_518),
];

/// Quarter-phase correction W, which is added for First Quarter and
/// subtracted for Last Quarter.
fn quarter_w(E: f64, M: f64, M1: f64, F: f64) -> f64 {
    // Releases before PR #12 lost the last three terms to a stray `;`.
    0.00306 - 0.00038 * E * M.cos() + 0.00026 * M1.cos() - 0.00002 * (M1 - M).cos()
        + 0.00002 * (M1 + M).cos()
        + 0.00002 * (2.0 * F).cos()
}

/// The mean phase with Meeus' `k` (an integer for New Moon, plus 0.25, 0.5
/// or 0.75 for the other phases) and the mean arguments at that instant
/// (eqs. 49.1 to 49.7 and 47.6). Chapter 54 (eclipses) starts from the
/// same quantities.
#[derive(Debug, Clone, Copy)]
pub(crate) struct MeanArgs {
    /// Julian centuries (TD) from J2000.0, `k / 1236.85`
    pub T: f64,
    /// JDE of the mean phase (eq. 49.1)
    pub jde: f64,
    /// Eccentricity factor E (eq. 47.6)
    pub E: f64,
    /// Sun's mean anomaly, radians in [0, 2π) (eq. 49.4)
    pub M: f64,
    /// Moon's mean anomaly, radians in [0, 2π) (eq. 49.5)
    pub M1: f64,
    /// Moon's argument of latitude, radians in [0, 2π) (eq. 49.6)
    pub F: f64,
    /// Longitude of the ascending node, radians in [0, 2π) (eq. 49.7)
    pub omega: f64,
}

pub(crate) fn mean_args(k: f64) -> MeanArgs {
    use sequence::reduced_rad;

    let T = k / 1236.85;

    MeanArgs {
        T,
        jde: mean_jde(k, T),
        E: sequence::eccentricity_factor(T),
        M: reduced_rad(2.5534 + 29.105_356_70 * k + T * T * (-0.000_001_4 - T * 0.000_000_11)),
        M1: reduced_rad(
            201.5643
                + 385.816_935_28 * k
                + T * T * (0.010_758_2 + T * (0.000_012_38 - T * 0.000_000_058)),
        ),
        F: reduced_rad(
            160.7108
                + 390.670_502_84 * k
                + T * T * (-0.001_611_8 + T * (-0.000_002_27 + T * 0.000_000_011)),
        ),
        omega: reduced_rad(124.7746 - 1.563_755_88 * k + T * T * (0.002_067_2 + T * 0.000_002_15)),
    }
}

fn terms(k: f64, phase: Phase) -> Terms {
    use sequence::reduced_rad;

    let MeanArgs {
        T,
        jde: mean,
        E,
        M,
        M1,
        F,
        omega,
    } = mean_args(k);

    let args = [M, M1, F, omega];
    let arg = |m: &[f64; 4]| m.iter().zip(args.iter()).map(|(a, b)| a * b).sum::<f64>();

    let (periodic, quarter) = match phase {
        Phase::New | Phase::Full => {
            let new = phase == Phase::New;
            let sum = NEW_FULL_TERMS
                .iter()
                .map(|(c_new, c_full, e_pow, m)| {
                    let c = if new { *c_new } else { *c_full };
                    c * E.powi(*e_pow) * arg(m).sin()
                })
                .sum();
            (sum, 0.0)
        }
        Phase::First | Phase::Last => {
            let sum = QUARTER_TERMS
                .iter()
                .map(|(c, e_pow, m)| c * E.powi(*e_pow) * arg(m).sin())
                .sum();
            let W = quarter_w(E, M, M1, F);
            (sum, if phase == Phase::First { W } else { -W })
        }
    };

    let planetary = PLANETARY_TERMS
        .iter()
        .enumerate()
        .map(|(i, (c, a0, rate))| {
            let mut a = a0 + rate * k;
            if i == 0 {
                a -= 0.009_173 * T * T;
            }
            c * reduced_rad(a).sin()
        })
        .sum();

    Terms {
        mean,
        E,
        M,
        M1,
        F,
        omega,
        periodic,
        quarter,
        planetary,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn deg(rad: f64) -> f64 {
        rad.to_degrees()
    }

    fn close(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() <= tol,
            "expected {expected} ± {tol}, got {actual}"
        );
    }

    // Meeus, example 49.a: New Moon of 1977 February, k = -283
    #[test]
    fn example_49a_intermediate_values() {
        let t = terms(-283.0, Phase::New);
        close(-283.0 / 1236.85, -0.22881, 0.000005);
        close(t.mean, 2443192.94102, 0.000005);
        close(t.E, 1.0005753, 0.00000005);
        close(deg(t.M), 45.7375, 0.00005);
        close(deg(t.M1), 95.3722, 0.00005);
        close(deg(t.F), 120.9584, 0.00005);
        close(deg(t.omega), 207.3176, 0.00005);
        close(t.periodic, -0.28916, 0.000005);
        close(t.quarter, 0.0, 0.0);
        close(t.planetary, -0.00068, 0.000005);
        close(t.jde(), 2443192.65118, 0.000005);
    }

    // Meeus, example 49.b: first Last Quarter of 2044, k = 544.75
    #[test]
    fn example_49b_intermediate_values() {
        let t = terms(544.75, Phase::Last);
        close(t.mean, 2467636.88597, 0.000005);
        close(t.periodic, -0.39153, 0.000005);
        close(t.quarter, -0.00251, 0.000005);
        close(t.planetary, -0.00007, 0.000005);
        close(t.jde(), 2467636.49186, 0.000005);
    }

    // W written out independently from the book; pins the fix of PR #12
    // (the last three cosine terms used to be dropped).
    #[test]
    fn quarter_correction_has_all_terms() {
        for &(E, M, M1, F) in &[
            (1.0, 0.0, 0.0, 0.0),
            (1.0005, 0.3, 1.7, 2.9),
            (0.9995, 4.0, 5.5, 0.1),
        ] {
            let expected = 0.00306 - 0.00038 * E * f64::cos(M) + 0.00026 * f64::cos(M1)
                - 0.00002 * f64::cos(M1 - M)
                + 0.00002 * f64::cos(M1 + M)
                + 0.00002 * f64::cos(2.0 * F);
            close(quarter_w(E, M, M1, F), expected, 1e-15);
        }
        // With every angle zero the three late terms add +0.00002.
        close(
            quarter_w(1.0, 0.0, 0.0, 0.0),
            0.00306 - 0.00038 + 0.00026 + 0.00002,
            1e-15,
        );
    }

    #[test]
    fn quarters_apply_w_with_opposite_signs() {
        for &k0 in &[-1000.0, -283.0, 0.0, 12.0, 544.0] {
            let first = terms(k0 + 0.25, Phase::First);
            let w = quarter_w(first.E, first.M, first.M1, first.F);
            close(first.quarter, w, 1e-15);
            let last = terms(k0 + 0.75, Phase::Last);
            let w = quarter_w(last.E, last.M, last.M1, last.F);
            close(last.quarter, -w, 1e-15);
            assert!(first.quarter > 0.0 && last.quarter < 0.0);
        }
    }

    // The quarter-phase term -0.00034 sin(2M' - M) carries a factor E in
    // the book; versions before 3.0 left it out.
    #[test]
    fn quarter_term_2m1_minus_m_has_factor_e() {
        assert_eq!(QUARTER_TERMS[10], (-0.00034, 1, [-1.0, 2.0, 0.0, 0.0]));
        // Far from 2000 E is clearly different from 1, and so is the sum.
        let t = terms(-3000.25, Phase::First);
        let arg = 2.0 * t.M1 - t.M;
        let without_e = t.periodic + 0.00034 * (t.E - 1.0) * arg.sin();
        assert!((t.E - 1.0).abs() > 0.005);
        assert!((t.periodic - without_e).abs() > 1e-7);
    }
}

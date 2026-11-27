//! Perigee and apogee of the Moon (Meeus, *Astronomical Algorithms*, 2nd
//! ed., chapter 50).
//!
//! The Moon is at perigee when its distance from the centre of the Earth is
//! a minimum and at apogee when it is a maximum. Meeus gives the mean
//! instants as a polynomial in the apsis number and corrects them with
//! periodic terms derived from ELP-2000/82. A second set of terms gives the
//! Moon's equatorial horizontal parallax at that instant, from which the
//! distance follows.
//!
//! # Numbering
//!
//! Apsides are numbered by the integer part of Meeus' `k`: `k = n` is
//! perigee number `n` and `k = n + 0.5` the apogee that follows it.
//! Perigee 0 is that of 1999 December 22; earlier apsides have negative
//! numbers. For example, Meeus' `k = -148.5` is the apogee with number
//! −149.
//!
//! # Time scale and accuracy
//!
//! Instants are Julian Ephemeris Days (Dynamical Time, TD). Meeus compared
//! 600 perigees and 600 apogees between 1977 and 2022 with ELP-2000/82:
//! the largest errors were 31 minutes (perigee) and 3 minutes (apogee) in
//! time, and 0.124″ (perigee) and 0.051″ (apogee) in parallax, i.e. about
//! 12 km and 6 km in distance.
//!
//! The perigee is much more strongly perturbed by the Sun than the apogee:
//! the true perigee can be up to about 45 hours from the mean one (the
//! apogee only about 13 hours), successive perigees are between about 24.7
//! and 28.5 days apart (apogees between 26.98 and 27.90 days), and perigee
//! distances range over roughly 356 370–370 350 km against 404 050–406 720
//! km for the apogee. The most extreme perigees and apogees occur around
//! the northern winter, when the Earth is nearest the Sun.
//!
//! # Example
//!
//! ```
//! use astro::lunar::apsis::{self, Apsis};
//!
//! // Meeus, example 50.a: the apogee of 1988 October (k = -148.5). The book
//! // prints 2447442.3537, but its own table gives 2447442.3543 (20h30m TD).
//! let apogee = apsis::true_apsis(-149, Apsis::Apogee);
//! assert!((apogee.jde - 2447442.3543).abs() < 1e-4); // 1988 Oct 7, 20h30m TD
//! let parallax_arcsec = apogee.parallax.to_degrees() * 3600.0;
//! assert!((parallax_arcsec - 3240.679).abs() < 1e-3);
//! assert!((apogee.distance - 405_977.0).abs() < 1.0); // km
//! ```

use super::sequence::{self, SlotIter, Timed};

/// Mean anomalistic month in days (the coefficient of `k` in eq. 50.1).
const ANOMALISTIC_MONTH: f64 = 27.554_549_89;

/// JDE of the mean perigee with `k = 0` (1999 December 22).
const EPOCH: f64 = 2_451_534.669_8;

/// Earth's equatorial radius in km, as used by Meeus to relate the Moon's
/// parallax to its distance (same value as `lunar::eq_hz_parllx`).
const EARTH_EQUATORIAL_RADIUS_KM: f64 = 6378.14;

/// One of the two apsides of the Moon's orbit
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Apsis {
    /// Least distance from the Earth
    Perigee,
    /// Greatest distance from the Earth
    Apogee,
}

impl Apsis {
    /// Fractional part of Meeus' `k` for this apsis (0 or 0.5)
    pub fn fraction(self) -> f64 {
        match self {
            Apsis::Perigee => 0.0,
            Apsis::Apogee => 0.5,
        }
    }
}

/// A passage of the Moon through perigee or apogee
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ApsisEvent {
    /// Perigee or apogee
    pub apsis: Apsis,
    /// Apsis number: the integer part of Meeus' `k` (see the
    /// [module docs](self))
    pub number: i64,
    /// Instant, as a Julian Ephemeris Day (TD). NaN if the event was
    /// requested for a non-finite Julian day.
    pub jde: f64,
    /// Moon's equatorial horizontal parallax at that instant *| in radians*
    pub parallax: f64,
    /// Distance between the centres of the Earth and the Moon at that
    /// instant *| in kilometres*, computed from the parallax with an
    /// equatorial radius of the Earth of 6378.14 km
    pub distance: f64,
}

impl ApsisEvent {
    /// Meeus' `k` for this event: `number + apsis.fraction()`
    pub fn k(&self) -> f64 {
        self.number as f64 + self.apsis.fraction()
    }
}

impl Timed for ApsisEvent {
    fn jde(&self) -> f64 {
        self.jde
    }
}

/// Approximate (non-integer) `k` for a year with decimals (Meeus,
/// eq. 50.2), e.g. `1988.75` for the beginning of October 1988
///
/// To get the apsis nearest to a given instant, use [`nearest_apsis`]
/// instead.
pub fn k_from_year(year: f64) -> f64 {
    (year - 1999.97) * 13.2555
}

/// Approximate (non-integer) `k` corresponding to a Julian day: the number
/// of mean anomalistic months since the mean perigee of 1999 December 22
pub fn k_from_jd(jd: f64) -> f64 {
    (jd - EPOCH) / ANOMALISTIC_MONTH
}

/// Instant of a *mean* perigee or apogee (Meeus, eq. 50.1), as a Julian
/// Ephemeris Day (TD)
///
/// # Arguments
///
/// * `number`: Apsis number (see the [module docs](self))
/// * `apsis` : Perigee or apogee
pub fn mean_apsis(number: i64, apsis: Apsis) -> f64 {
    let k = number as f64 + apsis.fraction();
    mean_jde(k, k / 1325.55)
}

/// Instant of a *true* perigee or apogee, with the Moon's parallax and
/// distance at that instant
///
/// # Arguments
///
/// * `number`: Apsis number (see the [module docs](self))
/// * `apsis` : Perigee or apogee
pub fn true_apsis(number: i64, apsis: Apsis) -> ApsisEvent {
    let t = apsis_terms(number as f64 + apsis.fraction(), apsis);
    let parallax = (t.parallax_arcsec / 3600.0).to_radians();
    ApsisEvent {
        apsis,
        number,
        jde: t.mean + t.periodic,
        parallax,
        distance: EARTH_EQUATORIAL_RADIUS_KM / parallax.sin(),
    }
}

/// The perigee or apogee nearest in time to `jd`
///
/// # Arguments
///
/// * `jd`   : Julian day. The difference between UT and TD is ignored when
///   choosing the nearest event.
/// * `apsis`: Perigee or apogee
///
/// If `jd` is not finite, the returned `jde`, `parallax` and `distance`
/// are NaN.
pub fn nearest_apsis(jd: f64, apsis: Apsis) -> ApsisEvent {
    if !jd.is_finite() {
        return ApsisEvent {
            apsis,
            number: 0,
            jde: f64::NAN,
            parallax: f64::NAN,
            distance: f64::NAN,
        };
    }
    sequence::nearest(jd, k_from_jd(jd) - apsis.fraction(), |number| {
        true_apsis(number, apsis)
    })
}

/// Iterator over all the perigees and apogees of the Moon, in chronological
/// order, whose instants lie in `[start_jde, end_jde)`
///
/// Both bounds are Julian Ephemeris Days. The iterator is empty if
/// `start_jde` is not finite or not less than `end_jde`; an infinite
/// `end_jde` gives an unbounded iterator. Perigees and apogees alternate.
pub fn apsides_between(start_jde: f64, end_jde: f64) -> Apsides {
    let first_slot = 2.0 * (k_from_jd(start_jde).floor() - 1.0);
    Apsides(SlotIter::new(start_jde, end_jde, first_slot, apsis_at_slot))
}

/// Iterator returned by [`apsides_between`]
#[derive(Debug, Clone)]
pub struct Apsides(SlotIter<ApsisEvent>);

impl Iterator for Apsides {
    type Item = ApsisEvent;

    fn next(&mut self) -> Option<ApsisEvent> {
        self.0.next()
    }
}

impl std::iter::FusedIterator for Apsides {}

fn apsis_at_slot(slot: i64) -> ApsisEvent {
    let apsis = if slot.rem_euclid(2) == 0 {
        Apsis::Perigee
    } else {
        Apsis::Apogee
    };
    true_apsis(slot.div_euclid(2), apsis)
}

/// Mean apsis polynomial (eq. 50.1)
fn mean_jde(k: f64, T: f64) -> f64 {
    EPOCH
        + ANOMALISTIC_MONTH * k
        + T * T * (-0.000_669_1 + T * (-0.000_001_098 + T * 0.000_000_005_2))
}

/// Intermediate quantities; the unit tests compare them with example 50.a.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(not(test), allow(dead_code))]
struct ApsisTerms {
    T: f64,
    mean: f64,
    /// Moon's mean elongation, radians in [0, 2π)
    D: f64,
    /// Sun's mean anomaly, radians in [0, 2π)
    M: f64,
    /// Moon's argument of latitude, radians in [0, 2π)
    F: f64,
    /// Sum of the periodic terms for the instant (days)
    periodic: f64,
    /// Equatorial horizontal parallax (arcseconds)
    parallax_arcsec: f64,
}

/// A periodic term: coefficient, its rate of change per Julian century
/// (T), and the multiples of D, M and F in the argument.
type Term = (f64, f64, [f64; 3]);

/// Periodic terms (sine, days) for the instant of perigee.
const PERIGEE_TIME: [Term; 60] = [
    (-1.6769, 0.0, [2.0, 0.0, 0.0]),
    (0.4589, 0.0, [4.0, 0.0, 0.0]),
    (-0.1856, 0.0, [6.0, 0.0, 0.0]),
    (0.0883, 0.0, [8.0, 0.0, 0.0]),
    (-0.0773, 0.00019, [2.0, -1.0, 0.0]),
    (0.0502, -0.00013, [0.0, 1.0, 0.0]),
    (-0.0460, 0.0, [10.0, 0.0, 0.0]),
    (0.0422, -0.00011, [4.0, -1.0, 0.0]),
    (-0.0256, 0.0, [6.0, -1.0, 0.0]),
    (0.0253, 0.0, [12.0, 0.0, 0.0]),
    (0.0237, 0.0, [1.0, 0.0, 0.0]),
    (0.0162, 0.0, [8.0, -1.0, 0.0]),
    (-0.0145, 0.0, [14.0, 0.0, 0.0]),
    (0.0129, 0.0, [0.0, 0.0, 2.0]),
    (-0.0112, 0.0, [3.0, 0.0, 0.0]),
    (-0.0104, 0.0, [10.0, -1.0, 0.0]),
    (0.0086, 0.0, [16.0, 0.0, 0.0]),
    (0.0069, 0.0, [12.0, -1.0, 0.0]),
    (0.0066, 0.0, [5.0, 0.0, 0.0]),
    (-0.0053, 0.0, [2.0, 0.0, 2.0]),
    (-0.0052, 0.0, [18.0, 0.0, 0.0]),
    (-0.0046, 0.0, [14.0, -1.0, 0.0]),
    (-0.0041, 0.0, [7.0, 0.0, 0.0]),
    (0.0040, 0.0, [2.0, 1.0, 0.0]),
    (0.0032, 0.0, [20.0, 0.0, 0.0]),
    (-0.0032, 0.0, [1.0, 1.0, 0.0]),
    (0.0031, 0.0, [16.0, -1.0, 0.0]),
    (-0.0029, 0.0, [4.0, 1.0, 0.0]),
    (0.0027, 0.0, [9.0, 0.0, 0.0]),
    (0.0027, 0.0, [4.0, 0.0, 2.0]),
    (-0.0027, 0.0, [2.0, -2.0, 0.0]),
    (0.0024, 0.0, [4.0, -2.0, 0.0]),
    (-0.0021, 0.0, [6.0, -2.0, 0.0]),
    (-0.0021, 0.0, [22.0, 0.0, 0.0]),
    (-0.0021, 0.0, [18.0, -1.0, 0.0]),
    (0.0019, 0.0, [6.0, 1.0, 0.0]),
    (-0.0018, 0.0, [11.0, 0.0, 0.0]),
    (-0.0014, 0.0, [8.0, 1.0, 0.0]),
    (-0.0014, 0.0, [4.0, 0.0, -2.0]),
    (-0.0014, 0.0, [6.0, 0.0, 2.0]),
    (0.0014, 0.0, [3.0, 1.0, 0.0]),
    (-0.0014, 0.0, [5.0, 1.0, 0.0]),
    (0.0013, 0.0, [13.0, 0.0, 0.0]),
    (0.0013, 0.0, [20.0, -1.0, 0.0]),
    (0.0011, 0.0, [3.0, 2.0, 0.0]),
    (-0.0011, 0.0, [4.0, -2.0, 2.0]),
    (-0.0010, 0.0, [1.0, 2.0, 0.0]),
    (-0.0009, 0.0, [22.0, -1.0, 0.0]),
    (-0.0008, 0.0, [0.0, 0.0, 4.0]),
    (0.0008, 0.0, [6.0, 0.0, -2.0]),
    (0.0008, 0.0, [2.0, 1.0, -2.0]),
    (0.0007, 0.0, [0.0, 2.0, 0.0]),
    (0.0007, 0.0, [0.0, -1.0, 2.0]),
    (0.0007, 0.0, [2.0, 0.0, 4.0]),
    (-0.0006, 0.0, [0.0, -2.0, 2.0]),
    (-0.0006, 0.0, [2.0, 2.0, -2.0]),
    (0.0006, 0.0, [24.0, 0.0, 0.0]),
    (0.0005, 0.0, [4.0, 0.0, -4.0]),
    (0.0005, 0.0, [2.0, 2.0, 0.0]),
    (-0.0004, 0.0, [1.0, -1.0, 0.0]),
];

/// Periodic terms (sine, days) for the instant of apogee.
const APOGEE_TIME: [Term; 32] = [
    (0.4392, 0.0, [2.0, 0.0, 0.0]),
    (0.0684, 0.0, [4.0, 0.0, 0.0]),
    (0.0456, -0.00011, [0.0, 1.0, 0.0]),
    (0.0426, -0.00011, [2.0, -1.0, 0.0]),
    (0.0212, 0.0, [0.0, 0.0, 2.0]),
    (-0.0189, 0.0, [1.0, 0.0, 0.0]),
    (0.0144, 0.0, [6.0, 0.0, 0.0]),
    (0.0113, 0.0, [4.0, -1.0, 0.0]),
    (0.0047, 0.0, [2.0, 0.0, 2.0]),
    (0.0036, 0.0, [1.0, 1.0, 0.0]),
    (0.0035, 0.0, [8.0, 0.0, 0.0]),
    (0.0034, 0.0, [6.0, -1.0, 0.0]),
    (-0.0034, 0.0, [2.0, 0.0, -2.0]),
    (0.0022, 0.0, [2.0, -2.0, 0.0]),
    (-0.0017, 0.0, [3.0, 0.0, 0.0]),
    (0.0013, 0.0, [4.0, 0.0, 2.0]),
    (0.0011, 0.0, [8.0, -1.0, 0.0]),
    (0.0010, 0.0, [4.0, -2.0, 0.0]),
    (0.0009, 0.0, [10.0, 0.0, 0.0]),
    (0.0007, 0.0, [3.0, 1.0, 0.0]),
    (0.0006, 0.0, [0.0, 2.0, 0.0]),
    (0.0005, 0.0, [2.0, 1.0, 0.0]),
    (0.0005, 0.0, [2.0, 2.0, 0.0]),
    (0.0004, 0.0, [6.0, 0.0, 2.0]),
    (0.0004, 0.0, [6.0, -2.0, 0.0]),
    (0.0004, 0.0, [10.0, -1.0, 0.0]),
    (-0.0004, 0.0, [5.0, 0.0, 0.0]),
    (-0.0004, 0.0, [4.0, 0.0, -2.0]),
    (0.0003, 0.0, [0.0, 1.0, 2.0]),
    (0.0003, 0.0, [12.0, 0.0, 0.0]),
    (0.0003, 0.0, [2.0, -1.0, 2.0]),
    (-0.0003, 0.0, [1.0, -1.0, 0.0]),
];

/// Constant term of the parallax at perigee (arcseconds).
const PERIGEE_PARALLAX_0: f64 = 3629.215;

/// Periodic terms (cosine, arcseconds) for the parallax at perigee.
const PERIGEE_PARALLAX: [Term; 46] = [
    (63.224, 0.0, [2.0, 0.0, 0.0]),
    (-6.990, 0.0, [4.0, 0.0, 0.0]),
    (2.834, -0.0071, [2.0, -1.0, 0.0]),
    (1.927, 0.0, [6.0, 0.0, 0.0]),
    (-1.263, 0.0, [1.0, 0.0, 0.0]),
    (-0.702, 0.0, [8.0, 0.0, 0.0]),
    (0.696, -0.0017, [0.0, 1.0, 0.0]),
    (-0.690, 0.0, [0.0, 0.0, 2.0]),
    (-0.629, 0.0016, [4.0, -1.0, 0.0]),
    (-0.392, 0.0, [2.0, 0.0, -2.0]),
    (0.297, 0.0, [10.0, 0.0, 0.0]),
    (0.260, 0.0, [6.0, -1.0, 0.0]),
    (0.201, 0.0, [3.0, 0.0, 0.0]),
    (-0.161, 0.0, [2.0, 1.0, 0.0]),
    (0.157, 0.0, [1.0, 1.0, 0.0]),
    (-0.138, 0.0, [12.0, 0.0, 0.0]),
    (-0.127, 0.0, [8.0, -1.0, 0.0]),
    (0.104, 0.0, [2.0, 0.0, 2.0]),
    (0.104, 0.0, [2.0, -2.0, 0.0]),
    (-0.079, 0.0, [5.0, 0.0, 0.0]),
    (0.068, 0.0, [14.0, 0.0, 0.0]),
    (0.067, 0.0, [10.0, -1.0, 0.0]),
    (0.054, 0.0, [4.0, 1.0, 0.0]),
    (-0.038, 0.0, [12.0, -1.0, 0.0]),
    (-0.038, 0.0, [4.0, -2.0, 0.0]),
    (0.037, 0.0, [7.0, 0.0, 0.0]),
    (-0.037, 0.0, [4.0, 0.0, 2.0]),
    (-0.035, 0.0, [16.0, 0.0, 0.0]),
    (-0.030, 0.0, [3.0, 1.0, 0.0]),
    (0.029, 0.0, [1.0, -1.0, 0.0]),
    (-0.025, 0.0, [6.0, 1.0, 0.0]),
    (0.023, 0.0, [0.0, 2.0, 0.0]),
    (0.023, 0.0, [14.0, -1.0, 0.0]),
    (-0.023, 0.0, [2.0, 2.0, 0.0]),
    (0.022, 0.0, [6.0, -2.0, 0.0]),
    (-0.021, 0.0, [2.0, -1.0, -2.0]),
    (-0.020, 0.0, [9.0, 0.0, 0.0]),
    (0.019, 0.0, [18.0, 0.0, 0.0]),
    (0.017, 0.0, [6.0, 0.0, 2.0]),
    (0.014, 0.0, [0.0, -1.0, 2.0]),
    (-0.014, 0.0, [16.0, -1.0, 0.0]),
    (0.013, 0.0, [4.0, 0.0, -2.0]),
    (0.012, 0.0, [8.0, 1.0, 0.0]),
    (0.011, 0.0, [11.0, 0.0, 0.0]),
    (0.010, 0.0, [5.0, 1.0, 0.0]),
    (-0.010, 0.0, [20.0, 0.0, 0.0]),
];

/// Constant term of the parallax at apogee (arcseconds).
const APOGEE_PARALLAX_0: f64 = 3245.251;

/// Periodic terms (cosine, arcseconds) for the parallax at apogee.
const APOGEE_PARALLAX: [Term; 17] = [
    (-9.147, 0.0, [2.0, 0.0, 0.0]),
    (-0.841, 0.0, [1.0, 0.0, 0.0]),
    (0.697, 0.0, [0.0, 0.0, 2.0]),
    (-0.656, 0.0016, [0.0, 1.0, 0.0]),
    (0.355, 0.0, [4.0, 0.0, 0.0]),
    (0.159, 0.0, [2.0, -1.0, 0.0]),
    (0.127, 0.0, [1.0, 1.0, 0.0]),
    (0.065, 0.0, [4.0, -1.0, 0.0]),
    (0.052, 0.0, [6.0, 0.0, 0.0]),
    (0.043, 0.0, [2.0, 1.0, 0.0]),
    (0.031, 0.0, [2.0, 0.0, 2.0]),
    (-0.023, 0.0, [2.0, 0.0, -2.0]),
    (0.022, 0.0, [2.0, -2.0, 0.0]),
    (0.019, 0.0, [2.0, 2.0, 0.0]),
    (-0.016, 0.0, [0.0, 2.0, 0.0]),
    (0.014, 0.0, [6.0, -1.0, 0.0]),
    (0.010, 0.0, [8.0, 0.0, 0.0]),
];

fn apsis_terms(k: f64, apsis: Apsis) -> ApsisTerms {
    use sequence::reduced_rad;

    let T = k / 1325.55;
    let TT = T * T;

    let D = reduced_rad(
        171.9179
            + 335.910_604_6 * k
            + TT * (-0.010_038_3 + T * (-0.000_011_56 + T * 0.000_000_055)),
    );
    let M = reduced_rad(347.3477 + 27.157_772_1 * k + TT * (-0.000_813_0 - T * 0.000_001_0));
    let F = reduced_rad(316.6109 + 364.528_791_1 * k + TT * (-0.012_505_3 - T * 0.000_014_8));

    let sum = |terms: &[Term], f: fn(f64) -> f64| -> f64 {
        terms
            .iter()
            .map(|(c, c_T, m)| (c + c_T * T) * f(m[0] * D + m[1] * M + m[2] * F))
            .sum()
    };

    let (periodic, parallax_arcsec) = match apsis {
        Apsis::Perigee => (
            sum(&PERIGEE_TIME, f64::sin),
            PERIGEE_PARALLAX_0 + sum(&PERIGEE_PARALLAX, f64::cos),
        ),
        Apsis::Apogee => (
            sum(&APOGEE_TIME, f64::sin),
            APOGEE_PARALLAX_0 + sum(&APOGEE_PARALLAX, f64::cos),
        ),
    };

    ApsisTerms {
        T,
        mean: mean_jde(k, T),
        D,
        M,
        F,
        periodic,
        parallax_arcsec,
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

    // Meeus, example 50.a: apogee of 1988 October, k = -148.5.
    //
    // Erratum: the book prints -0.4654 for the sum of the apogee terms of
    // Table 50.A and hence JDE 2447442.3537 (20h29m TD), but the 32 terms as
    // printed sum to -0.4648, giving JDE 2447442.3543 = 20h30m TD, which is
    // the "exact" value quoted in the same example. Every other value of the
    // example (T, mean JDE, D, M, F, parallax) is reproduced.
    #[test]
    fn example_50a_intermediate_values() {
        let t = apsis_terms(-148.5, Apsis::Apogee);
        close(t.T, -0.112029, 0.0000005);
        close(t.mean, 2447442.8191, 0.00005);
        close(t.D.to_degrees(), 329.1930, 0.00005);
        close(t.M.to_degrees(), 274.4185, 0.00005);
        close(t.F.to_degrees(), 184.0853, 0.00005);
        close(t.periodic, -0.4648, 0.00005);
        close(t.parallax_arcsec, 3240.679, 0.0005);
        close(t.mean + t.periodic, 2447442.3543, 0.00005);
    }
}

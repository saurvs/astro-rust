//! Maximum declinations of the Moon (Meeus, *Astronomical Algorithms*, 2nd
//! ed., chapter 52).
//!
//! The Moon's orbit is inclined by only about 5° to the ecliptic, so its
//! declination swings between a northern and a southern extreme, about two
//! weeks apart, in every tropical month (27.3 days). How far it swings
//! depends on where the nodes of the orbit are in their 18.6-year cycle:
//! the extremes range from roughly 18° to 28.6°. Meeus expresses both the
//! instants and the extreme values as a mean term plus periodic terms. The
//! values refer to the centre of the Moon as seen from the centre of the
//! Earth (geocentric).
//!
//! # Numbering
//!
//! Each tropical month has one northern and one southern extreme with the
//! same integer number `k`, the southern one coming first; `k = 0` is the
//! month of January 2000 (southern extreme on January 6, northern on
//! January 19). Earlier extremes have negative numbers. Unlike chapters
//! 49–51, `k` is always an integer here.
//!
//! # Time scale and accuracy
//!
//! Instants are Julian Ephemeris Days (Dynamical Time, TD). Meeus compared
//! 600 northern and 600 southern extremes between 1977 and 2022 with an
//! accurate theory: the errors were at most 10 minutes in time and 26″ in
//! declination, and usually much smaller. Meeus adds that the coefficients
//! assume the obliquity of the ecliptic of 2000.0, and that this alone
//! keeps the error in time below half an hour between −1000 and +5000.
//!
//! In practice the series is only that good near its 1977–2022 fitting
//! interval. Compared with the chapter 47 lunar theory (and, for one date,
//! with JPL Horizons) the errors grow to about 50 minutes in time and 3′ in
//! declination by 1910 and by 2085, while staying within Meeus' figures
//! for 1977–2022. For other epochs, use the result as a starting point and
//! refine it with an accurate lunar position.
//!
//! # Example
//!
//! ```
//! use astro::lunar::declination::{self, Extreme};
//!
//! // Meeus, example 52.a: greatest northern declination of 1988 December.
//! let max = declination::max_declination(-148, Extreme::North);
//! assert!((max.jde - 2447518.3347).abs() < 1e-4); // 1988 Dec 22, 20h02m TD
//! assert!((max.declination.to_degrees() - 28.1562).abs() < 1e-4);
//!
//! // Meeus, example 52.b: a southern extreme. Its declination is negative.
//! let min = declination::max_declination(659, Extreme::South);
//! assert!((min.jde - 2469553.0834).abs() < 1e-4); // 2049 Apr 21, 14h TD
//! assert!((min.declination.to_degrees() + 22.1384).abs() < 1e-4);
//! ```

use super::sequence::{self, SlotIter, Timed};

/// Mean tropical month in days (the coefficient of `k` for the instant).
const TROPICAL_MONTH: f64 = 27.321_582_247;

/// JDE of the mean northern extreme with `k = 0` (2000 January 19).
const EPOCH_NORTH: f64 = 2_451_562.589_7;

/// JDE of the mean southern extreme with `k = 0` (2000 January 5).
const EPOCH_SOUTH: f64 = 2_451_548.928_9;

/// Greatest northern or greatest southern declination
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Extreme {
    /// Greatest northern declination
    North,
    /// Greatest southern declination
    South,
}

/// An extreme declination of the Moon and its instant
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaxDeclination {
    /// Northern or southern extreme
    pub extreme: Extreme,
    /// Meeus' `k` (an integer; see the [module docs](self))
    pub number: i64,
    /// Instant, as a Julian Ephemeris Day (TD). NaN if the event was
    /// requested for a non-finite Julian day.
    pub jde: f64,
    /// Geocentric declination of the Moon at that instant *| in radians*:
    /// positive for a northern extreme, negative for a southern one (Meeus
    /// prints the absolute value)
    pub declination: f64,
}

impl Timed for MaxDeclination {
    fn jde(&self) -> f64 {
        self.jde
    }
}

/// Approximate (non-integer) `k` for a year with decimals (Meeus,
/// eq. 52.1), e.g. `1988.95` for mid-December 1988
///
/// Round it to an integer to obtain the extreme numbers of that month. To
/// get the extreme nearest to a given instant, use
/// [`nearest_max_declination`] instead.
pub fn k_from_year(year: f64) -> f64 {
    (year - 2000.03) * 13.3686
}

/// Approximate (non-integer) `k` corresponding to a Julian day, for the
/// given kind of extreme: the number of mean tropical months since its
/// mean instant for `k = 0`
pub fn k_from_jd(jd: f64, extreme: Extreme) -> f64 {
    let epoch = match extreme {
        Extreme::North => EPOCH_NORTH,
        Extreme::South => EPOCH_SOUTH,
    };
    (jd - epoch) / TROPICAL_MONTH
}

/// Instant and value of a greatest northern or southern declination of the
/// Moon
///
/// # Arguments
///
/// * `number` : Meeus' `k` (see the [module docs](self))
/// * `extreme`: Northern or southern extreme
pub fn max_declination(number: i64, extreme: Extreme) -> MaxDeclination {
    let t = declination_terms(number as f64, extreme);
    let declination = t.declination_deg.to_radians();
    MaxDeclination {
        extreme,
        number,
        jde: t.jde,
        declination: match extreme {
            Extreme::North => declination,
            Extreme::South => -declination,
        },
    }
}

/// The northern or southern extreme of declination nearest in time to `jd`
///
/// # Arguments
///
/// * `jd`     : Julian day. The difference between UT and TD is ignored
///   when choosing the nearest event.
/// * `extreme`: Northern or southern extreme
///
/// If `jd` is not finite, the returned `jde` and `declination` are NaN.
pub fn nearest_max_declination(jd: f64, extreme: Extreme) -> MaxDeclination {
    if !jd.is_finite() {
        return MaxDeclination {
            extreme,
            number: 0,
            jde: f64::NAN,
            declination: f64::NAN,
        };
    }
    sequence::nearest(jd, k_from_jd(jd, extreme), |number| {
        max_declination(number, extreme)
    })
}

/// Iterator over all the northern and southern extremes of the Moon's
/// declination, in chronological order, whose instants lie in
/// `[start_jde, end_jde)`
///
/// Both bounds are Julian Ephemeris Days. The iterator is empty if
/// `start_jde` is not finite or not less than `end_jde`; an infinite
/// `end_jde` gives an unbounded iterator. Northern and southern extremes
/// alternate; within one `number` the southern extreme comes first.
pub fn max_declinations_between(start_jde: f64, end_jde: f64) -> MaxDeclinations {
    let first_slot = 2.0 * (k_from_jd(start_jde, Extreme::South).floor() - 1.0);
    MaxDeclinations(SlotIter::new(
        start_jde,
        end_jde,
        first_slot,
        extreme_at_slot,
    ))
}

/// Iterator returned by [`max_declinations_between`]
#[derive(Debug, Clone)]
pub struct MaxDeclinations(SlotIter<MaxDeclination>);

impl Iterator for MaxDeclinations {
    type Item = MaxDeclination;

    fn next(&mut self) -> Option<MaxDeclination> {
        self.0.next()
    }
}

impl std::iter::FusedIterator for MaxDeclinations {}

fn extreme_at_slot(slot: i64) -> MaxDeclination {
    let extreme = if slot.rem_euclid(2) == 0 {
        Extreme::South
    } else {
        Extreme::North
    };
    max_declination(slot.div_euclid(2), extreme)
}

/// Intermediate quantities; the unit tests compare them with Meeus'
/// examples.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(not(test), allow(dead_code))]
struct DeclinationTerms {
    T: f64,
    /// Moon's mean elongation, radians in [0, 2π)
    D: f64,
    /// Sun's mean anomaly, radians in [0, 2π)
    M: f64,
    /// Moon's mean anomaly, radians in [0, 2π)
    M1: f64,
    /// Moon's argument of latitude, radians in [0, 2π)
    F: f64,
    E: f64,
    jde: f64,
    /// Absolute value of the extreme declination (degrees)
    declination_deg: f64,
}

/// A periodic term: coefficient for the northern extreme, coefficient for
/// the southern one, whether it is multiplied by E, whether it is a cosine
/// (otherwise a sine), and the multiples of D, M, M′ and F.
type Term = (f64, f64, bool, bool, [f64; 4]);

const SIN: bool = false;
const COS: bool = true;

/// Periodic terms (days) for the instant of the extreme declination.
const TIME_TERMS: [Term; 44] = [
    (0.8975, -0.8975, false, COS, [0.0, 0.0, 0.0, 1.0]),
    (-0.4726, -0.4726, false, SIN, [0.0, 0.0, 1.0, 0.0]),
    (-0.1030, -0.1030, false, SIN, [0.0, 0.0, 0.0, 2.0]),
    (-0.0976, -0.0976, false, SIN, [2.0, 0.0, -1.0, 0.0]),
    (-0.0462, 0.0541, false, COS, [0.0, 0.0, 1.0, -1.0]),
    (-0.0461, 0.0516, false, COS, [0.0, 0.0, 1.0, 1.0]),
    (-0.0438, -0.0438, false, SIN, [2.0, 0.0, 0.0, 0.0]),
    (0.0162, 0.0112, true, SIN, [0.0, 1.0, 0.0, 0.0]),
    (-0.0157, 0.0157, false, COS, [0.0, 0.0, 0.0, 3.0]),
    (0.0145, 0.0023, false, SIN, [0.0, 0.0, 1.0, 2.0]),
    (0.0136, -0.0136, false, COS, [2.0, 0.0, 0.0, -1.0]),
    (-0.0095, 0.0110, false, COS, [2.0, 0.0, -1.0, -1.0]),
    (-0.0091, 0.0091, false, COS, [2.0, 0.0, -1.0, 1.0]),
    (-0.0089, 0.0089, false, COS, [2.0, 0.0, 0.0, 1.0]),
    (0.0075, 0.0075, false, SIN, [0.0, 0.0, 2.0, 0.0]),
    (-0.0068, -0.0030, false, SIN, [0.0, 0.0, 1.0, -2.0]),
    (0.0061, -0.0061, false, COS, [0.0, 0.0, 2.0, -1.0]),
    (-0.0047, -0.0047, false, SIN, [0.0, 0.0, 1.0, 3.0]),
    (-0.0043, -0.0043, true, SIN, [2.0, -1.0, -1.0, 0.0]),
    (-0.0040, 0.0040, false, COS, [0.0, 0.0, 1.0, -2.0]),
    (-0.0037, -0.0037, false, SIN, [2.0, 0.0, -2.0, 0.0]),
    (0.0031, -0.0031, false, SIN, [0.0, 0.0, 0.0, 1.0]),
    (0.0030, 0.0030, false, SIN, [2.0, 0.0, 1.0, 0.0]),
    (-0.0029, 0.0029, false, COS, [0.0, 0.0, 1.0, 2.0]),
    (-0.0029, -0.0029, true, SIN, [2.0, -1.0, 0.0, 0.0]),
    (-0.0027, -0.0027, false, SIN, [0.0, 0.0, 1.0, 1.0]),
    (0.0024, 0.0024, true, SIN, [0.0, 1.0, -1.0, 0.0]),
    (-0.0021, -0.0021, false, SIN, [0.0, 0.0, 1.0, -3.0]),
    (0.0019, -0.0019, false, SIN, [0.0, 0.0, 2.0, 1.0]),
    (0.0018, -0.0006, false, COS, [2.0, 0.0, -2.0, -1.0]),
    (0.0018, -0.0018, false, SIN, [0.0, 0.0, 0.0, 3.0]),
    (0.0017, -0.0017, false, COS, [0.0, 0.0, 1.0, 3.0]),
    (0.0017, 0.0017, false, COS, [0.0, 0.0, 2.0, 0.0]),
    (-0.0014, 0.0014, false, COS, [2.0, 0.0, -1.0, 0.0]),
    (0.0013, -0.0013, false, COS, [2.0, 0.0, 1.0, 1.0]),
    (0.0013, -0.0013, false, COS, [0.0, 0.0, 1.0, 0.0]),
    (0.0012, 0.0012, false, SIN, [0.0, 0.0, 3.0, 1.0]),
    (0.0011, 0.0011, false, SIN, [2.0, 0.0, -1.0, 1.0]),
    (-0.0011, 0.0011, false, COS, [2.0, 0.0, -2.0, 0.0]),
    (0.0010, 0.0010, false, COS, [1.0, 0.0, 0.0, 1.0]),
    (0.0010, 0.0010, true, SIN, [0.0, 1.0, 1.0, 0.0]),
    (-0.0009, -0.0009, false, SIN, [2.0, 0.0, 0.0, -2.0]),
    (0.0007, -0.0007, false, COS, [0.0, 0.0, 2.0, 1.0]),
    (-0.0007, -0.0007, false, COS, [0.0, 0.0, 3.0, 1.0]),
];

/// Periodic terms (degrees) for the value of the extreme declination.
///
/// Two rows are not symmetric between north and south in the book:
/// sin(M′ − F) (−0.0124 / −0.0015) and sin(M′ + F) (+0.0006 / +0.0037).
/// They are used as printed; examples 52.a–c are reproduced with them.
const DECLINATION_TERMS: [Term; 37] = [
    (5.1093, -5.1093, false, SIN, [0.0, 0.0, 0.0, 1.0]),
    (0.2658, 0.2658, false, COS, [0.0, 0.0, 0.0, 2.0]),
    (0.1448, -0.1448, false, SIN, [2.0, 0.0, 0.0, -1.0]),
    (-0.0322, 0.0322, false, SIN, [0.0, 0.0, 0.0, 3.0]),
    (0.0133, 0.0133, false, COS, [2.0, 0.0, 0.0, -2.0]),
    (0.0125, 0.0125, false, COS, [2.0, 0.0, 0.0, 0.0]),
    (-0.0124, -0.0015, false, SIN, [0.0, 0.0, 1.0, -1.0]),
    (-0.0101, 0.0101, false, SIN, [0.0, 0.0, 1.0, 2.0]),
    (0.0097, -0.0097, false, COS, [0.0, 0.0, 0.0, 1.0]),
    (-0.0087, 0.0087, true, SIN, [2.0, 1.0, 0.0, -1.0]),
    (0.0074, 0.0074, false, SIN, [0.0, 0.0, 1.0, 3.0]),
    (0.0067, 0.0067, false, SIN, [1.0, 0.0, 0.0, 1.0]),
    (0.0063, -0.0063, false, SIN, [0.0, 0.0, 1.0, -2.0]),
    (0.0060, -0.0060, true, SIN, [2.0, -1.0, 0.0, -1.0]),
    (-0.0057, 0.0057, false, SIN, [2.0, 0.0, -1.0, -1.0]),
    (-0.0056, -0.0056, false, COS, [0.0, 0.0, 1.0, 1.0]),
    (0.0052, -0.0052, false, COS, [0.0, 0.0, 1.0, 2.0]),
    (0.0041, -0.0041, false, COS, [0.0, 0.0, 2.0, 1.0]),
    (-0.0040, -0.0040, false, COS, [0.0, 0.0, 1.0, -3.0]),
    (0.0038, -0.0038, false, COS, [0.0, 0.0, 2.0, -1.0]),
    (-0.0034, 0.0034, false, COS, [0.0, 0.0, 1.0, -2.0]),
    (-0.0029, -0.0029, false, SIN, [0.0, 0.0, 2.0, 0.0]),
    (0.0029, 0.0029, false, SIN, [0.0, 0.0, 3.0, 1.0]),
    (-0.0028, 0.0028, true, COS, [2.0, 1.0, 0.0, -1.0]),
    (-0.0028, -0.0028, false, COS, [0.0, 0.0, 1.0, -1.0]),
    (-0.0023, 0.0023, false, COS, [0.0, 0.0, 0.0, 3.0]),
    (-0.0021, 0.0021, false, SIN, [2.0, 0.0, 0.0, 1.0]),
    (0.0019, 0.0019, false, COS, [0.0, 0.0, 1.0, 3.0]),
    (0.0018, 0.0018, false, COS, [1.0, 0.0, 0.0, 1.0]),
    (0.0017, -0.0017, false, SIN, [0.0, 0.0, 2.0, -1.0]),
    (0.0015, 0.0015, false, COS, [0.0, 0.0, 3.0, 1.0]),
    (0.0014, 0.0014, false, COS, [2.0, 0.0, 2.0, 1.0]),
    (-0.0012, 0.0012, false, SIN, [2.0, 0.0, -2.0, -1.0]),
    (-0.0012, -0.0012, false, COS, [0.0, 0.0, 2.0, 0.0]),
    (-0.0010, 0.0010, false, COS, [0.0, 0.0, 1.0, 0.0]),
    (-0.0010, -0.0010, false, SIN, [0.0, 0.0, 0.0, 2.0]),
    (0.0006, 0.0037, false, SIN, [0.0, 0.0, 1.0, 1.0]),
];

fn declination_terms(k: f64, extreme: Extreme) -> DeclinationTerms {
    use sequence::reduced_rad;

    let T = k / 1336.86;
    let TT = T * T;
    let north = extreme == Extreme::North;

    // The polynomials differ between north and south only in the constant.
    let (D0, M0, M10, F0, jde0) = if north {
        (152.2029, 14.8591, 4.6881, 325.8867, EPOCH_NORTH)
    } else {
        (345.6676, 1.3951, 186.2100, 145.1633, EPOCH_SOUTH)
    };
    let D = reduced_rad(D0 + 333.070_554_6 * k + TT * (-0.000_421_4 + T * 0.000_000_11));
    let M = reduced_rad(M0 + 26.928_159_2 * k + TT * (-0.000_035_5 - T * 0.000_000_10));
    let M1 = reduced_rad(M10 + 356.956_279_4 * k + TT * (0.010_306_6 + T * 0.000_012_51));
    let F = reduced_rad(F0 + 1.446_780_7 * k + TT * (-0.002_069_0 - T * 0.000_002_15));
    let E = sequence::eccentricity_factor(T);

    let sum = |terms: &[Term]| -> f64 {
        terms
            .iter()
            .map(|&(c_north, c_south, with_e, cos, m)| {
                let c = if north { c_north } else { c_south };
                let arg = m[0] * D + m[1] * M + m[2] * M1 + m[3] * F;
                let trig = if cos { arg.cos() } else { arg.sin() };
                c * trig * if with_e { E } else { 1.0 }
            })
            .sum()
    };

    let jde =
        jde0 + TROPICAL_MONTH * k + TT * (0.000_119_804 - T * 0.000_000_141) + sum(&TIME_TERMS);
    let declination_deg = 23.6961 - 0.013_004 * T + sum(&DECLINATION_TERMS);

    DeclinationTerms {
        T,
        D,
        M,
        M1,
        F,
        E,
        jde,
        declination_deg,
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

    // Meeus, example 52.a: greatest northern declination, 1988 December
    #[test]
    fn example_52a_intermediate_values() {
        let t = declination_terms(-148.0, Extreme::North);
        close(t.T, -0.110707, 0.0000005);
        close(t.D.to_degrees(), 177.7608, 0.00005);
        close(t.M.to_degrees(), 349.4915, 0.00005);
        close(t.M1.to_degrees(), 95.1589, 0.00005);
        close(t.F.to_degrees(), 111.7631, 0.00005);
        close(t.E, 1.000278, 0.0000005);
        // 2447518.33465: the book rounds it up to .3347.
        close(t.jde, 2447518.3347, 0.0001);
        close(t.declination_deg, 28.1562, 0.00005);
    }

    // Meeus, example 52.b: greatest southern declination, k = 659
    #[test]
    fn example_52b_southern() {
        let t = declination_terms(659.0, Extreme::South);
        close(t.jde, 2469553.0834, 0.00005);
        close(t.declination_deg, 22.1384, 0.00005);
    }

    // Meeus, example 52.c: greatest northern declination, March of -4.
    // The printed terms give JDE 1719672.14125, 13 s before the book's
    // 1719672.1414. No single misread term explains the difference without
    // breaking examples 52.a and 52.b, so it is probably rounding in the
    // book's own computation this far (20 centuries) from the epoch.
    #[test]
    fn example_52c_ancient() {
        let t = declination_terms(-26788.0, Extreme::North);
        close(t.jde, 1719672.1414, 0.0002);
        close(t.declination_deg, 28.9739, 0.00005);
    }
}

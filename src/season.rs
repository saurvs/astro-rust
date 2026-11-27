//! Equinoxes and solstices (Meeus, *Astronomical Algorithms*, 2nd ed.,
//! chapter 27)
//!
//! An equinox or a solstice is the instant when the Sun's *apparent*
//! geocentric longitude (aberration and nutation included) is a multiple
//! of 90°. Three methods are provided:
//!
//! - [`mean_instant`]: the instant of the "mean" equinox or solstice, from
//!   polynomials in the year (Meeus, tables 27.A and 27.B);
//! - [`approx_instant`]: the mean instant corrected by 24 periodic terms
//!   (table 27.C), within about one minute of the true instant for
//!   1951–2050;
//! - [`instant`]: the true instant, found by iterating on the Sun's apparent
//!   longitude from the complete VSOP87 theory (Meeus, eq. 27.1), within
//!   about a second for the years −2000 to +6000.
//!
//! [`duration`] gives the length of the astronomical seasons.
//!
//! All instants are Julian Ephemeris Days (dynamical time, TT); subtract
//! ΔT to obtain UT. Years are astronomical years (1 BC = year 0, 2 BC =
//! year −1, ...).
//!
//! # Example
//!
//! ```
//! use astro::season::{self, Event};
//!
//! // Meeus, examples 27.a and 27.b: the June solstice of 1962
//! let jde = season::instant(1962, Event::JuneSolstice);
//!
//! // 1962 June 21 at 21h24m42s TD
//! let expected = 2437836.5 + (21.0 + 24.0 / 60.0 + 42.0 / 3600.0) / 24.0;
//! assert!((jde - expected).abs() * 86400.0 < 1.5);
//! ```

use crate::sun;
use crate::time;

/// An equinox or a solstice
///
/// The names refer to the month in which the event takes place, which does
/// not depend on the hemisphere.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Event {
    /// March (northward) equinox: apparent solar longitude 0°
    MarchEquinox,
    /// June solstice: apparent solar longitude 90°
    JuneSolstice,
    /// September (southward) equinox: apparent solar longitude 180°
    SeptemberEquinox,
    /// December solstice: apparent solar longitude 270°
    DecemberSolstice,
}

impl Event {
    /// The four events in the order in which they happen in a calendar year
    pub const ALL: [Event; 4] = [
        Event::MarchEquinox,
        Event::JuneSolstice,
        Event::SeptemberEquinox,
        Event::DecemberSolstice,
    ];

    fn index(self) -> usize {
        match self {
            Event::MarchEquinox => 0,
            Event::JuneSolstice => 1,
            Event::SeptemberEquinox => 2,
            Event::DecemberSolstice => 3,
        }
    }

    /// Returns the Sun's apparent geocentric longitude at the event
    /// *| in radians*: 0, π/2, π or 3π/2
    pub fn solar_longitude(self) -> f64 {
        self.index() as f64 * std::f64::consts::FRAC_PI_2
    }
}

/// An astronomical season, named as in the northern hemisphere
///
/// Each season runs from one equinox or solstice to the next; in the
/// southern hemisphere the names are swapped (spring ↔ autumn,
/// summer ↔ winter).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Season {
    /// From the March equinox to the June solstice
    Spring,
    /// From the June solstice to the September equinox
    Summer,
    /// From the September equinox to the December solstice
    Autumn,
    /// From the December solstice to the March equinox of the next year
    Winter,
}

impl Season {
    /// The four seasons, starting with the one that begins in March
    pub const ALL: [Season; 4] = [
        Season::Spring,
        Season::Summer,
        Season::Autumn,
        Season::Winter,
    ];

    /// Returns the event at which the season begins
    pub fn start(self) -> Event {
        match self {
            Season::Spring => Event::MarchEquinox,
            Season::Summer => Event::JuneSolstice,
            Season::Autumn => Event::SeptemberEquinox,
            Season::Winter => Event::DecemberSolstice,
        }
    }

    /// Returns the event at which the season ends (for [`Season::Winter`],
    /// the March equinox of the following year)
    pub fn end(self) -> Event {
        match self {
            Season::Spring => Event::JuneSolstice,
            Season::Summer => Event::SeptemberEquinox,
            Season::Autumn => Event::DecemberSolstice,
            Season::Winter => Event::MarchEquinox,
        }
    }
}

/// Polynomial coefficients of the mean instants, in `Y = year / 1000`, for
/// the years −1000 to +1000 (Meeus, table 27.A). One row per [`Event`].
const MEAN_BEFORE_1000: [[f64; 5]; 4] = [
    [1721139.29189, 365242.13740, 0.06134, 0.00111, -0.00071],
    [1721233.25401, 365241.72562, -0.05323, 0.00907, 0.00025],
    [1721325.70455, 365242.49558, -0.11677, -0.00297, 0.00074],
    [1721414.39987, 365242.88257, -0.00769, -0.00933, -0.00006],
];

/// Polynomial coefficients of the mean instants, in
/// `Y = (year − 2000) / 1000`, for the years +1000 to +3000 (Meeus,
/// table 27.B). One row per [`Event`].
const MEAN_AFTER_1000: [[f64; 5]; 4] = [
    [2451623.80984, 365242.37404, 0.05169, -0.00411, -0.00057],
    [2451716.56767, 365241.62603, 0.00325, 0.00888, -0.00030],
    [2451810.21715, 365242.01767, -0.11575, 0.00337, 0.00078],
    [2451900.05952, 365242.74049, -0.06223, -0.00823, 0.00032],
];

/// Periodic terms `A cos(B + C T)` (B, C in degrees, T in Julian centuries)
/// correcting the mean instants (Meeus, table 27.C).
const PERIODIC_TERMS: [[f64; 3]; 24] = [
    [485.0, 324.96, 1934.136],
    [203.0, 337.23, 32964.467],
    [199.0, 342.08, 20.186],
    [182.0, 27.85, 445267.112],
    [156.0, 73.14, 45036.886],
    [136.0, 171.52, 22518.443],
    [77.0, 222.54, 65928.934],
    [74.0, 296.72, 3034.906],
    [70.0, 243.58, 9037.513],
    [58.0, 119.81, 33718.147],
    [52.0, 297.17, 150.678],
    [50.0, 21.02, 2281.226],
    [45.0, 247.54, 29929.562],
    [44.0, 325.15, 31555.956],
    [29.0, 60.93, 4443.417],
    [18.0, 155.12, 67555.328],
    [17.0, 288.79, 4562.452],
    [16.0, 198.04, 62894.029],
    [14.0, 199.76, 31436.921],
    [12.0, 95.39, 14577.848],
    [12.0, 287.11, 31931.756],
    [12.0, 320.81, 34777.259],
    [9.0, 227.73, 1222.114],
    [8.0, 15.45, 16859.074],
];

/// Computes the instant of the mean equinox or solstice of a year (Meeus,
/// tables 27.A and 27.B)
///
/// Table 27.A is used for years before 1000 and table 27.B from 1000 on.
/// They are meant for the years −1000 to +3000; a few centuries beyond
/// that range the error is still small, but it keeps growing.
///
/// # Returns
///
/// * `JDE0`: Julian Ephemeris Day of the mean event
///
/// # Arguments
///
/// * `year`: Astronomical year (an integer: 1 BC is year 0)
/// * `event`: The equinox or solstice
pub fn mean_instant(year: i32, event: Event) -> f64 {
    let (c, y) = if year < 1000 {
        (&MEAN_BEFORE_1000[event.index()], year as f64 / 1000.0)
    } else {
        (
            &MEAN_AFTER_1000[event.index()],
            (year as f64 - 2000.0) / 1000.0,
        )
    };

    Horner_eval!(y, c[0], c[1], c[2], c[3], c[4])
}

/// Correction in days added to a mean instant: `0.00001 S / Δλ`, where `S`
/// is the sum of the periodic terms and `Δλ` accounts for the variable
/// speed of the Sun (Meeus, chapter 27).
fn periodic_correction(jde0: f64) -> f64 {
    let t = time::julian_cent(jde0);

    let w = (35999.373 * t - 2.47).to_radians();
    let delta_lambda = 1.0 + 0.0334 * w.cos() + 0.0007 * (2.0 * w).cos();

    let s: f64 = PERIODIC_TERMS
        .iter()
        .map(|&[a, b, c]| a * (b + c * t).rem_euclid(360.0).to_radians().cos())
        .sum();

    0.00001 * s / delta_lambda
}

/// Computes an approximate instant of an equinox or solstice (Meeus,
/// chapter 27)
///
/// The mean instant of [`mean_instant`] plus the periodic correction of
/// table 27.C. For 1951–2050 the error is below one minute (Meeus, table
/// 27.D); use [`instant`] for a more accurate value.
///
/// # Returns
///
/// * `JDE`: Julian Ephemeris Day of the event
///
/// # Arguments
///
/// * `year`: Astronomical year (an integer: 1 BC is year 0)
/// * `event`: The equinox or solstice
pub fn approx_instant(year: i32, event: Event) -> f64 {
    let jde0 = mean_instant(year, event);

    jde0 + periodic_correction(jde0)
}

/// Computes the instant of an equinox or solstice accurately (Meeus,
/// chapter 27, eq. 27.1)
///
/// Starting from [`approx_instant`], the Sun's apparent longitude `λ` is
/// computed with [`sun::apparent_position`] (complete VSOP87, FK5,
/// nutation and aberration) and the time is corrected by
/// `58 sin(k·90° − λ)` days until the correction is below 0.01 s.
///
/// The result is limited by VSOP87 (about 1″ in longitude, i.e. about 24 s
/// in time, at the ends of the range −2000 to +6000) and is within a
/// second or so near the present.
///
/// # Returns
///
/// * `JDE`: Julian Ephemeris Day of the event
///
/// # Arguments
///
/// * `year`: Astronomical year (an integer: 1 BC is year 0)
/// * `event`: The equinox or solstice
pub fn instant(year: i32, event: Event) -> f64 {
    let target = event.solar_longitude();
    let mut jde = approx_instant(year, event);

    for _ in 0..20 {
        let lambda = sun::apparent_position(jde).apparent_longitude;
        let correction = 58.0 * (target - lambda).sin();
        jde += correction;
        if correction.abs() < 1e-7 {
            break;
        }
    }

    jde
}

/// Computes the duration of an astronomical season (Meeus, chapter 27)
///
/// The difference between the accurate instants ([`instant`]) of the
/// events that begin and end the season. The winter of `year` begins at
/// the December solstice of `year` and ends at the March equinox of
/// `year + 1`.
///
/// # Returns
///
/// * `days`: Duration of the season in days (of dynamical time)
///
/// # Arguments
///
/// * `year`: Astronomical year in which the season begins
/// * `season`: The season (northern hemisphere names)
///
/// # Example
///
/// ```
/// use astro::season::{self, Season};
///
/// // In 2000 the northern summer lasted about 93.65 days and the
/// // winter about 88.99 days.
/// assert!((season::duration(2000, Season::Summer) - 93.65).abs() < 0.01);
/// assert!((season::duration(2000, Season::Winter) - 88.99).abs() < 0.01);
/// ```
pub fn duration(year: i32, season: Season) -> f64 {
    let end_year = if season == Season::Winter {
        year.saturating_add(1)
    } else {
        year
    };

    instant(end_year, season.end()) - instant(year, season.start())
}

#[cfg(test)]
mod tests {
    use super::*;

    // Meeus, example 27.a: intermediate values for the June solstice of 1962
    #[test]
    fn example_27a_intermediate_values() {
        let jde0 = mean_instant(1962, Event::JuneSolstice);
        assert!((jde0 - 2437837.38589).abs() < 1e-5);

        let t = time::julian_cent(jde0);
        assert!((t + 0.375294021).abs() < 1e-9);

        let w = (35999.373 * t - 2.47).to_radians();
        let delta_lambda = 1.0 + 0.0334 * w.cos() + 0.0007 * (2.0 * w).cos();
        assert!((delta_lambda - 0.9681).abs() < 1e-4);

        // S = +635, and the correction is 0.00001 S / Δλ = 0.00656 day
        let correction = periodic_correction(jde0);
        assert!((correction * delta_lambda / 0.00001 - 635.0).abs() < 0.5);
        assert!((jde0 + correction - 2437837.39245).abs() < 1e-5);
    }

    #[test]
    fn event_longitudes_and_indices() {
        for (k, e) in Event::ALL.iter().enumerate() {
            assert_eq!(e.index(), k);
            assert!((e.solar_longitude() - k as f64 * std::f64::consts::FRAC_PI_2).abs() < 1e-15);
        }
    }
}

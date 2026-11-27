/*
Copyright (c) 2016 Saurav Sachidanand

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

//! Astronomical, physical and timekeeping constants
//!
//! Most values follow Meeus, *Astronomical Algorithms*, 2nd ed., Appendix I.
//! Where a quantity now has an exact modern definition (the astronomical
//! unit of IAU 2012, the SI speed of light) that definition is used; the
//! difference from the book's rounded values is noted on each item.
//!
//! Conventions: lengths in metres, times in days (ephemeris days of
//! 86400 SI seconds), angles in radians, unless the item says otherwise.
//!
//! ```
//! use astro::consts;
//!
//! // Light-time across one astronomical unit, in seconds (book: 499.0048 s)
//! let secs = consts::LIGHT_TIME_PER_AU * 86400.0;
//! assert!((secs - 499.0048).abs() < 1e-4);
//!
//! // One parsec is about 3.2616 light-years
//! assert!((consts::PARSEC_IN_LIGHT_YEARS - 3.2616).abs() < 1e-4);
//! ```

use std::f64::consts::PI;

pub mod iau1976;
pub mod wgs72;
pub mod wgs84;

// ---------------------------------------------------------------------------
// Mathematical
// ---------------------------------------------------------------------------

/// Number of arcseconds in one radian, `648000 / π` (≈ 206264.806247)
pub const ARCSEC_PER_RADIAN: f64 = 648000.0 / PI;

// ---------------------------------------------------------------------------
// Distances
// ---------------------------------------------------------------------------

/// Astronomical unit *| in meters*
///
/// The exact value fixed by IAU 2012 resolution B2. Meeus' Appendix I
/// uses 149 597 870 km (the IAU 1976 value), which differs by 0.7 km.
pub const ASTRONOMICAL_UNIT: f64 = 149_597_870_700.0;

/// Time taken by light to travel one astronomical unit *| in days*
///
/// ≈ 0.0057755183 day = 499.0048 s, the factor used in light-time
/// corrections (`τ = LIGHT_TIME_PER_AU × Δ`, Δ in AU).
pub const LIGHT_TIME_PER_AU: f64 = ASTRONOMICAL_UNIT / SPEED_OF_LIGHT / 86400.0;

/// One parsec *| in astronomical units*
///
/// A parsec is the distance at which one AU subtends one arcsecond, so this
/// equals [`ARCSEC_PER_RADIAN`] (≈ 206264.8 AU).
pub const PARSEC_IN_AU: f64 = ARCSEC_PER_RADIAN;

/// One parsec *| in meters* (≈ 30.8568 × 10¹⁵ m)
pub const PARSEC: f64 = ASTRONOMICAL_UNIT * PARSEC_IN_AU;

/// One light-year, the distance light travels in a Julian year of 365.25
/// days *| in meters* (≈ 9.4607 × 10¹⁵ m)
pub const LIGHT_YEAR: f64 = SPEED_OF_LIGHT * JULIAN_YEAR * 86400.0;

/// One parsec *| in light-years* (≈ 3.2616)
pub const PARSEC_IN_LIGHT_YEARS: f64 = PARSEC / LIGHT_YEAR;

/// One light-year *| in astronomical units* (≈ 63241)
pub const LIGHT_YEAR_IN_AU: f64 = LIGHT_YEAR / ASTRONOMICAL_UNIT;

/// Mean distance between the centers of the Earth and the Moon
/// *| in meters*
pub const MEAN_EARTH_MOON_DIST: f64 = 384_400_000.0;

/// Diameter of the Sun *| in meters*
pub const SUN_DIAMETER: f64 = 1_392_000_000.0;

/// Diameter of the Moon *| in meters*
pub const MOON_DIAMETER: f64 = 3_476_000.0;

/// Ratio of the Moon's mean radius to the Earth's equatorial radius
///
/// The value 0.272481, adopted for eclipse work in the *Astronomical
/// Ephemeris* from 1963 (Meeus, chapter 55).
pub const MOON_EARTH_RADIUS_RATIO: f64 = 0.272481;

// ---------------------------------------------------------------------------
// Time
// ---------------------------------------------------------------------------

/// Julian Ephemeris Day of the standard epoch J2000.0
/// (2000 January 1.5 TT)
pub const J2000: f64 = 2451545.0;

/// Length of the Julian year *| in days*
pub const JULIAN_YEAR: f64 = 365.25;

/// Length of the Julian century *| in days*
pub const JULIAN_CENTURY: f64 = 36525.0;

/// Length of the mean sidereal day *| in mean solar days*
/// (23h 56m 04.0905s)
pub const SIDEREAL_DAY: f64 = 0.9972695663;

/// Length of the mean solar day *| in sidereal days*
pub const SOLAR_DAY_IN_SIDEREAL_DAYS: f64 = 1.00273790935;

/// Tropical year (equinox to equinox) at epoch 2000.0 *| in days*
pub const TROPICAL_YEAR: f64 = 365.242190;

/// Sidereal year (fixed star to fixed star) at epoch 2000.0 *| in days*
pub const SIDEREAL_YEAR: f64 = 365.256363;

/// Anomalistic year (perihelion to perihelion) at epoch 2000.0
/// *| in days*
pub const ANOMALISTIC_YEAR: f64 = 365.259636;

/// Tropical month of the Moon (equinox to equinox) at epoch 2000.0
/// *| in days*
pub const TROPICAL_MONTH: f64 = 27.321582;

/// Sidereal month of the Moon (fixed star to fixed star) at epoch 2000.0
/// *| in days*
pub const SIDEREAL_MONTH: f64 = 27.321662;

/// Anomalistic month of the Moon (perigee to perigee) at epoch 2000.0
/// *| in days*
pub const ANOMALISTIC_MONTH: f64 = 27.554550;

/// Draconic month of the Moon (node to node) at epoch 2000.0 *| in days*
pub const DRACONIC_MONTH: f64 = 27.212221;

/// Synodic month of the Moon (New Moon to New Moon) at epoch 2000.0
/// *| in days*
///
/// The mean-phase formula of Meeus' chapter 49 uses the more precise
/// 29.530588861 days.
pub const SYNODIC_MONTH: f64 = 29.530589;

// ---------------------------------------------------------------------------
// Varia
// ---------------------------------------------------------------------------

/// Mean equatorial horizontal parallax of the Sun (the solar parallax),
/// 8.79415″ *| in radians*
pub const SOLAR_PARALLAX: f64 = 8.79415 / ARCSEC_PER_RADIAN;

/// Constant of annual aberration, 20.4955″ *| in radians*
pub const ABERRATION_CONSTANT: f64 = 20.4955 / ARCSEC_PER_RADIAN;

/// Gaussian gravitational constant `k` *| in radians*
///
/// The daily mean motion, in radians per day, of a massless body in an
/// unperturbed orbit of 1 AU around the Sun.
pub const GAUSS_GRAV: f64 = 0.01720209895;

/// Gaussian gravitational constant `k` converted to degrees
/// *| in degrees* (≈ 0.9856076686)
pub const GAUSS_GRAV_DEG: f64 = GAUSS_GRAV * 180.0 / PI;

/// Speed of light in vacuum *| in meters per second* (exact, SI)
pub const SPEED_OF_LIGHT: f64 = 299792458.0;

/// Ratio of the mass of the Earth to the mass of the Moon
pub const EARTH_MOON_MASS_RATIO: f64 = 81.3007;

/// Ratio of the mass of the Sun to the mass of the Earth
pub const SUN_EARTH_MASS_RATIO: f64 = 332946.0;

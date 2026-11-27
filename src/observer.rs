//! High-level helpers for an observer on the Earth: where a body is in the sky, and when it rises, transits and sets.
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

//!
//! The lower-level modules of this crate each implement one step of a
//! calculation (a planet's heliocentric position, nutation, sidereal time,
//! parallax, ...). This module chains them together so that the two most
//! common questions can be answered in one call:
//!
//! * **Where is it?** [`topocentric_position`] gives the azimuth and altitude
//!   (elevation) of the Sun, the Moon, a planet or a star for an observer at a
//!   given instant, with and without atmospheric refraction.
//! * **When does it rise and set?** [`rise_transit_set`] gives the times of
//!   rising, transit and setting during a day, handles bodies that never rise
//!   or never set (polar day and night), and computes civil, nautical and
//!   astronomical twilight.
//!
//! # Conventions
//!
//! * Angles are in radians.
//! * Geographic longitude is measured **positive west** of Greenwich, as in
//!   Meeus. Most maps and GPS receivers give it positive *east*: use
//!   [`Observer::from_degrees_east`] to avoid sign mistakes.
//! * Times are Julian Days in UT (`jd_ut`) unless the argument is called
//!   `jde` (Julian Ephemeris Day, Dynamical Time). Conversions between the two
//!   use the crate's ΔT model, [`time::delta_t_from_decimal_year`].
//! * Azimuths are given both in Meeus' convention (measured westward from the
//!   south) and in the compass convention (measured eastward from the north).
//!
//! # Accuracy
//!
//! Positions use the full VSOP87 theory for the Sun and the planets (better
//! than 1″), Meeus' abridged ELP-2000/82 for the Moon (about 10″), and the
//! IAU 1980 nutation. Rise and set times follow Meeus' definitions (chapter
//! 15) and are limited by the unpredictable refraction near the horizon,
//! which can shift them by a minute or more; the computed values themselves
//! are good to a few seconds.
//!
//! # Examples
//!
//! Sunrise and sunset at Greenwich on the June solstice of 2024:
//!
//! ```
//! use astro::observer::{self, Body, Horizon, Observer};
//! use astro::time;
//!
//! let greenwich = Observer::from_degrees_east(51.4769, -0.0005, 46.0);
//! let jd = time::julian_day(&time::Date {
//!     year: 2024,
//!     month: time::Month::June,
//!     decimal_day: 21.0, // 0h UT
//!     cal_type: time::CalType::Gregorian,
//! });
//!
//! let sun = observer::rise_transit_set(Body::Sun, &greenwich, jd, Horizon::Standard);
//! let hours_ut = |t: f64| (t - jd) * 24.0;
//! let rise = hours_ut(sun.rise.unwrap()); // 3.72, i.e. 03:43 UT
//! let set = hours_ut(sun.set.unwrap()); // 20.36, i.e. 20:21 UT
//! assert!((rise - (3.0 + 43.0 / 60.0)).abs() < 2.0 / 60.0);
//! assert!((set - (20.0 + 21.0 / 60.0)).abs() < 2.0 / 60.0);
//! ```
//!
//! The azimuth and elevation of Mars seen from Palomar Observatory:
//!
//! ```
//! use astro::observer::{self, Body, Observer};
//!
//! // 33°21′22″ N, 116°51′47″ W, 1706 m
//! let palomar = Observer::from_degrees_east(33.356111, -116.863056, 1706.0);
//! let jd_ut = 2452879.63681; // 2003 August 28, 3h17m UT
//! let mars = observer::topocentric_position(Body::Mars, &palomar, jd_ut);
//!
//! println!(
//!     "azimuth {:.2}° (from north), elevation {:.2}°",
//!     mars.azimuth_from_north.to_degrees(),
//!     mars.apparent_altitude.to_degrees()
//! );
//! assert!(mars.altitude > 0.0); // Mars is up
//! ```

use crate::coords;
use crate::ecliptic;
use crate::lunar;
use crate::nutation;
use crate::parallax;
use crate::planet;
use crate::precess;
use crate::sun;
use crate::time;
use crate::transit;
use crate::{aberr, atmos};
use std::f64::consts::{PI, TAU};

const DEG: f64 = PI / 180.0;

/// The astronomical unit in kilometres (IAU 2012).
const AU_KM: f64 = 149_597_870.7;

/// Julian Ephemeris Day of the epoch J2000.0.
const J2000: f64 = 2451545.0;

/// Standard pressure (mbar) and temperature (°C) for refraction.
const STD_PRESSURE: f64 = 1010.0;
const STD_TEMPERATURE: f64 = 10.0;

/// Below this true altitude no refraction is applied.
const MIN_REFRACTION_ALT: f64 = -DEG;

/// An observer on the surface of the Earth
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Observer {
    /// Geographic (geodetic) latitude *| in radians*, positive north
    pub lat: f64,
    /// Geographic longitude *| in radians*, measured **positive west** of
    /// Greenwich and negative east (Meeus' convention)
    pub long: f64,
    /// Height above sea level *| in metres*. It is used for the parallax of
    /// the body; rise and set times are for a sea-level horizon.
    pub height: f64,
}

impl Observer {
    /// Creates an observer from a latitude and a longitude in radians
    ///
    /// # Arguments
    ///
    /// * `lat`: Geographic latitude *| in radians*, positive north
    /// * `long`: Geographic longitude *| in radians*, **positive west**
    /// * `height`: Height above sea level *| in metres*
    pub fn new(lat: f64, long: f64, height: f64) -> Observer {
        Observer { lat, long, height }
    }

    /// Creates an observer from a latitude and a longitude in degrees, with
    /// the longitude measured **positive east** of Greenwich as on maps and
    /// GPS receivers
    ///
    /// The longitude is converted to the crate's west-positive convention.
    ///
    /// # Arguments
    ///
    /// * `lat_deg`: Geographic latitude *| in degrees*, positive north
    /// * `east_long_deg`: Geographic longitude *| in degrees*, positive
    ///   **east**
    /// * `height`: Height above sea level *| in metres*
    ///
    /// # Example
    ///
    /// ```
    /// use astro::observer::Observer;
    ///
    /// let sydney = Observer::from_degrees_east(-33.8688, 151.2093, 0.0);
    /// assert!(sydney.long < 0.0); // east of Greenwich is negative here
    /// ```
    pub fn from_degrees_east(lat_deg: f64, east_long_deg: f64, height: f64) -> Observer {
        Observer {
            lat: lat_deg.to_radians(),
            long: -east_long_deg.to_radians(),
            height,
        }
    }

    /// Returns the observer's latitude and longitude as a
    /// [`coords::GeographPoint`] (longitude positive west)
    pub fn geograph_point(&self) -> coords::GeographPoint {
        coords::GeographPoint {
            long: self.long,
            lat: self.lat,
        }
    }
}

/// A celestial body whose position can be computed
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Body {
    /// The Sun
    Sun,
    /// The Moon
    Moon,
    /// Mercury
    Mercury,
    /// Venus
    Venus,
    /// Mars
    Mars,
    /// Jupiter
    Jupiter,
    /// Saturn
    Saturn,
    /// Uranus
    Uranus,
    /// Neptune
    Neptune,
    /// A star, or any other body far outside the solar system, given by
    /// its mean right ascension and declination for the equator and equinox
    /// of J2000.0 (as in modern catalogues), in radians
    ///
    /// Precession, nutation and annual aberration are applied. Proper motion
    /// and annual parallax are not: if the star's proper motion matters, add
    /// it to the coordinates for the date of interest first.
    Star {
        /// Right ascension, J2000.0 *| in radians*
        asc: f64,
        /// Declination, J2000.0 *| in radians*
        dec: f64,
    },
}

/// The apparent geocentric position of a body
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ApparentPosition {
    /// Apparent right ascension, referred to the true equator and equinox of
    /// the date *| in radians*, in `[0, 2π)`
    pub asc: f64,
    /// Apparent declination, referred to the true equator of the date
    /// *| in radians*
    pub dec: f64,
    /// Distance from the centre of the Earth *| in AU*. For a planet this is
    /// the light-time distance; for a star it is `f64::INFINITY`.
    pub distance: f64,
    /// Equatorial horizontal parallax *| in radians*; zero for a star
    pub eq_hz_parallax: f64,
}

/// The position of a body in the sky of an observer
///
/// All quantities are topocentric (as seen from the observer's place, with
/// the parallax taken into account) and apparent (corrected for light-time,
/// aberration and nutation).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TopocentricPosition {
    /// Azimuth measured westward from the **south** (Meeus' convention:
    /// south 0, west π/2, north π, east 3π/2) *| in radians*, in `[0, 2π)`
    pub azimuth_from_south: f64,
    /// Azimuth measured eastward from the **north** (the compass convention:
    /// north 0, east π/2, south π, west 3π/2) *| in radians*, in `[0, 2π)`
    pub azimuth_from_north: f64,
    /// Altitude (elevation) of the body's centre above the horizon, without
    /// atmospheric refraction ("airless" or "true" altitude) *| in radians*,
    /// in `[−π/2, π/2]`
    pub altitude: f64,
    /// Altitude of the body's centre as seen through the atmosphere, i.e.
    /// `altitude` plus the refraction for a pressure of 1010 mbar and a
    /// temperature of 10 °C *| in radians*. No refraction is added when the
    /// true altitude is below −1°. See [`refracted_altitude`] for other
    /// conditions.
    pub apparent_altitude: f64,
    /// Local hour angle, measured westward from the meridian *| in radians*,
    /// in `[−π, π)`
    pub hour_angle: f64,
    /// Topocentric apparent right ascension *| in radians*, in `[0, 2π)`
    pub asc: f64,
    /// Topocentric apparent declination *| in radians*
    pub dec: f64,
    /// Distance from the observer *| in AU*; `f64::INFINITY` for a star
    pub distance: f64,
}

/// The altitude that defines the rise and set of a body
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Horizon {
    /// The apparent rise and set of the body, using Meeus' standard
    /// altitudes for its centre: −0°34′ for stars and planets (refraction),
    /// −0°50′ for the Sun (refraction and semidiameter: the upper limb
    /// touches the horizon), and `0.7275 π − 0°34′` for the Moon, where `π`
    /// is its horizontal parallax. See [`transit::standard_altitude`].
    Standard,
    /// Civil twilight: the centre of the body (normally the Sun) at −6°
    CivilTwilight,
    /// Nautical twilight: the centre of the body at −12°
    NauticalTwilight,
    /// Astronomical twilight: the centre of the body at −18°
    AstronomicalTwilight,
    /// A custom geometric (airless, geocentric) altitude of the body's
    /// centre *| in radians*
    Altitude(f64),
}

/// Whether a body crossed the chosen altitude during the day
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HorizonStatus {
    /// The body rose, set, or both during the day.
    RisesOrSets,
    /// The body stayed above the altitude all day (polar day, midnight Sun,
    /// circumpolar star, or a twilight that lasts all night).
    AlwaysAbove,
    /// The body stayed below the altitude all day (polar night).
    AlwaysBelow,
}

/// The times of rising, transit and setting of a body during one day
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RiseTransitSet {
    /// Julian Day (UT) at which the body rose above the altitude, if it did
    /// during the day. For a twilight horizon this is dawn.
    pub rise: Option<f64>,
    /// Julian Day (UT) of the upper transit (culmination across the local
    /// meridian), if it happened during the day. Circumpolar bodies transit
    /// too.
    pub transit: Option<f64>,
    /// Julian Day (UT) at which the body set below the altitude, if it did
    /// during the day. For a twilight horizon this is dusk.
    pub set: Option<f64>,
    /// Whether the body crossed the altitude, or stayed above or below it
    /// all day
    pub status: HorizonStatus,
}

/// Computes the apparent geocentric position of a body
///
/// The position is referred to the true equator and equinox of the date,
/// and includes light-time, annual aberration and nutation:
///
/// * the Sun: VSOP87, conversion to FK5, nutation and aberration (Meeus,
///   ch. 25, "higher accuracy");
/// * the Moon: Meeus' ELP-2000/82 series plus nutation in longitude (ch. 47);
/// * a planet: VSOP87, with the Earth and the planet both taken at the
///   time the light left the planet, which accounts for light-time and
///   aberration together; then FK5 and nutation (ch. 33);
/// * a star: precession from J2000.0 (ch. 21), annual aberration including
///   the eccentricity terms (eq. 23.2, applied as a vector so that it holds
///   at any latitude) and nutation in longitude (ch. 23).
///
/// # Arguments
///
/// * `body`: The [`Body`]
/// * `jde`: Julian Ephemeris Day (Dynamical Time)
///
/// # Example
///
/// The Moon on 1992 April 12 at 0h TD (Meeus, example 47.a):
///
/// ```
/// use astro::observer::{self, Body};
///
/// let moon = observer::geocentric_apparent_position(Body::Moon, 2448724.5);
/// assert!((moon.asc.to_degrees() - 134.688470).abs() < 1e-5);
/// assert!((moon.dec.to_degrees() - 13.768368).abs() < 1e-5);
/// ```
pub fn geocentric_apparent_position(body: Body, jde: f64) -> ApparentPosition {
    geocentric(body, jde, &Nutation::at(jde))
}

/// Computes where a body appears in the sky of an observer at an instant
///
/// This chains the whole calculation: the apparent geocentric position of
/// the body ([`geocentric_apparent_position`]) at `jd_ut + ΔT`, the apparent
/// sidereal time at Greenwich (Meeus, ch. 12), the topocentric correction
/// for parallax (ch. 40), the local horizontal coordinates (ch. 13) and the
/// atmospheric refraction (ch. 16). ΔT comes from
/// [`time::delta_t_from_decimal_year`]; to use
/// your own value, call [`geocentric_apparent_position`] and
/// [`topocentric_from_geocentric`] directly.
///
/// # Arguments
///
/// * `body`: The [`Body`]
/// * `observer`: The [`Observer`] (longitude positive west)
/// * `jd_ut`: Julian Day (UT)
///
/// # Example
///
/// ```
/// use astro::observer::{self, Body, Observer};
///
/// // Sirius from Cape Town, 2024 January 1, 22h UT
/// let cape_town = Observer::from_degrees_east(-33.9249, 18.4241, 0.0);
/// let sirius = Body::Star {
///     asc: (6.0_f64 + 45.0 / 60.0 + 8.917 / 3600.0).to_radians() * 15.0,
///     dec: -(16.0_f64 + 42.0 / 60.0 + 58.02 / 3600.0).to_radians(),
/// };
/// let pos = observer::topocentric_position(sirius, &cape_town, 2460311.416667);
/// assert!(pos.altitude > 0.0);
/// println!("Sirius: azimuth {:.1}° from north", pos.azimuth_from_north.to_degrees());
/// ```
pub fn topocentric_position(body: Body, observer: &Observer, jd_ut: f64) -> TopocentricPosition {
    let jde = jd_ut + delta_t_secs(jd_ut) / 86400.0;
    let nut = Nutation::at(jde);
    let geo = geocentric(body, jde, &nut);
    topocentric_from_geocentric(&geo, observer, nut.apparent_sidereal(jd_ut))
}

/// Computes the topocentric position of a body from its apparent geocentric
/// position
///
/// Applies the parallax for the observer's place and height (Meeus, eqs.
/// 40.2 and 40.3), then transforms to horizontal coordinates (eqs. 13.5 and
/// 13.6) and adds the refraction for standard conditions (eq. 16.4).
///
/// # Arguments
///
/// * `position`: Apparent geocentric position, e.g. from
///   [`geocentric_apparent_position`]
/// * `observer`: The [`Observer`] (longitude positive west)
/// * `apprnt_greenwich_sidr`: Apparent sidereal time at Greenwich at the
///   same instant *| in radians*
pub fn topocentric_from_geocentric(
    position: &ApparentPosition,
    observer: &Observer,
    apprnt_greenwich_sidr: f64,
) -> TopocentricPosition {
    let geo_point = coords::EqPoint {
        asc: position.asc,
        dec: position.dec,
    };

    let (eq_point, dist_ratio) = if position.eq_hz_parallax == 0.0 {
        (geo_point, 1.0)
    } else {
        let topo = parallax::topocent_eq_coords(
            &geo_point,
            position.eq_hz_parallax,
            &observer.geograph_point(),
            observer.height,
            apprnt_greenwich_sidr,
        );
        (
            topo,
            topocentric_distance_ratio(position, observer, apprnt_greenwich_sidr),
        )
    };

    let asc = norm_2pi(eq_point.asc);
    let dec = eq_point.dec;
    let hour_angle = wrap_pi(coords::hr_angl_frm_observer_long(
        apprnt_greenwich_sidr,
        observer.long,
        asc,
    ));
    let azimuth_from_south = norm_2pi(coords::az_frm_eq(hour_angle, dec, observer.lat));
    let altitude = coords::alt_frm_eq(hour_angle, dec, observer.lat);

    TopocentricPosition {
        azimuth_from_south,
        azimuth_from_north: norm_2pi(azimuth_from_south + PI),
        altitude,
        apparent_altitude: refracted_altitude(altitude, STD_PRESSURE, STD_TEMPERATURE),
        hour_angle,
        asc,
        dec,
        distance: position.distance * dist_ratio,
    }
}

/// Computes the altitude at which a body appears through the atmosphere
///
/// Adds to the true (airless) altitude the refraction given by
/// Sæmundsson's formula (Meeus, eq. 16.4), scaled for pressure and
/// temperature. Refraction this close to the horizon is uncertain by several
/// arcminutes; below a true altitude of −1° the formula is meaningless and
/// the true altitude is returned unchanged.
///
/// # Returns
///
/// * `apparent_alt`: Apparent altitude *| in radians*
///
/// # Arguments
///
/// * `true_alt`: True (airless) altitude *| in radians*
/// * `pressure`: Atmospheric pressure at the observer *| in millibars*
///   (1010 is standard)
/// * `temperature`: Air temperature at the observer *| in degrees Celsius*
///   (10 is standard)
pub fn refracted_altitude(true_alt: f64, pressure: f64, temperature: f64) -> f64 {
    if !(MIN_REFRACTION_ALT..=PI / 2.0).contains(&true_alt) {
        return true_alt;
    }
    let refraction = atmos::refrac_frm_true_alt(true_alt)
        * atmos::refrac_pressr_temp_factor(pressure, temperature);

    true_alt + refraction.max(0.0)
}

/// Computes the times of rising, transit and setting of a body during one
/// day
///
/// The body's altitude is followed through the 24 hours starting at
/// `jd_start`, and every instant at which it crosses the altitude chosen by
/// `horizon` is found to better than a second. Unlike the classical method
/// of Meeus' chapter 15 (available as [`transit::day_fraction`]), which
/// assumes one rise and one set per day, this also gives correct answers for
/// the Moon (which skips a rise or a set about once a month), near the polar
/// circles, and at the poles. The definitions are Meeus': geocentric
/// altitude of the body's centre, with his standard altitudes.
///
/// # Arguments
///
/// * `body`: The [`Body`]
/// * `observer`: The [`Observer`] (longitude positive west). The height is
///   ignored: times are for a sea-level horizon.
/// * `jd_start`: Julian Day (UT) at which the search window starts. For the
///   events of a UT calendar date pass the Julian Day of 0h UT of that date
///   (a number ending in .5). For a local calendar date pass the Julian Day of
///   local midnight: `jd_0h_ut - utc_offset_hours / 24`.
/// * `horizon`: The altitude that defines rising and setting; see
///   [`Horizon`]
///
/// # Returns
///
/// A [`RiseTransitSet`] with the first rise, transit and set in the window
/// `[jd_start, jd_start + 1)`, each `None` if it did not happen then, and a
/// status telling whether the body stayed above or below the altitude all
/// day. Consecutive windows share no events, so iterating over days finds
/// every event exactly once. Non-finite inputs give no events (and an
/// unspecified status) rather than a panic.
///
/// # Example
///
/// The midnight Sun and the polar night at Tromsø (69.65° N):
///
/// ```
/// use astro::observer::{self, Body, Horizon, HorizonStatus, Observer};
///
/// let tromso = Observer::from_degrees_east(69.6492, 18.9553, 0.0);
/// let june_21 = 2460482.5; // 2024 June 21, 0h UT
/// let dec_21 = 2460665.5; // 2024 December 21, 0h UT
///
/// let day = observer::rise_transit_set(Body::Sun, &tromso, june_21, Horizon::Standard);
/// assert_eq!(day.status, HorizonStatus::AlwaysAbove);
/// assert_eq!((day.rise, day.set), (None, None));
/// assert!(day.transit.is_some());
///
/// let night = observer::rise_transit_set(Body::Sun, &tromso, dec_21, Horizon::Standard);
/// assert_eq!(night.status, HorizonStatus::AlwaysBelow);
///
/// // ... but there is civil twilight around noon.
/// let twilight = observer::rise_transit_set(Body::Sun, &tromso, dec_21, Horizon::CivilTwilight);
/// assert!(twilight.rise.is_some() && twilight.set.is_some());
/// ```
pub fn rise_transit_set(
    body: Body,
    observer: &Observer,
    jd_start: f64,
    horizon: Horizon,
) -> RiseTransitSet {
    let track = Track::new(body, jd_start);
    let eval = |t: f64| {
        let (geo, gast) = track.at(t);
        let hour_angle = wrap_pi(coords::hr_angl_frm_observer_long(
            gast,
            observer.long,
            geo.asc,
        ));
        let alt = coords::alt_frm_eq(hour_angle, geo.dec, observer.lat);
        let h0 = match horizon {
            Horizon::Standard => match body {
                Body::Sun => transit::STD_ALT_SUN,
                Body::Moon => {
                    transit::standard_altitude(&transit::TransitBody::Moon, geo.eq_hz_parallax)
                }
                _ => transit::STD_ALT_STAR_OR_PLANET,
            },
            Horizon::CivilTwilight => -6.0 * DEG,
            Horizon::NauticalTwilight => -12.0 * DEG,
            Horizon::AstronomicalTwilight => -18.0 * DEG,
            Horizon::Altitude(h) => h,
        };
        Sample {
            t,
            hour_angle,
            above: alt - h0,
        }
    };

    let end = jd_start + 1.0;
    let mut samples: Vec<Sample> = (0..=SAMPLES_PER_DAY)
        .map(|i| eval(jd_start + i as f64 / SAMPLES_PER_DAY as f64))
        .collect();

    // Upper and lower culminations: where the hour angle crosses 0 and π.
    // The altitude is (nearly) monotonic between them, so adding them to the
    // samples makes sure no brief rise above the horizon is missed.
    let mut transit = None;
    let mut culminations = Vec::new();
    for w in samples.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        if a.hour_angle < 0.0 && b.hour_angle >= 0.0 && b.hour_angle - a.hour_angle < PI {
            let t = find_root(|t| eval(t).hour_angle, a.t, a.hour_angle, b.t, b.hour_angle);
            if transit.is_none() && t < end {
                transit = Some(t);
            }
            culminations.push(t);
        }
        let (ga, gb) = (wrap_pi(a.hour_angle - PI), wrap_pi(b.hour_angle - PI));
        if ga < 0.0 && gb >= 0.0 && gb - ga < PI {
            culminations.push(find_root(
                |t| wrap_pi(eval(t).hour_angle - PI),
                a.t,
                ga,
                b.t,
                gb,
            ));
        }
    }
    samples.extend(culminations.into_iter().map(eval));
    samples.sort_by(|a, b| a.t.total_cmp(&b.t));

    let mut rise = None;
    let mut set = None;
    for w in samples.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        let rising = a.above < 0.0 && b.above >= 0.0;
        let setting = a.above >= 0.0 && b.above < 0.0;
        if (rising && rise.is_none()) || (setting && set.is_none()) {
            let t = find_root(|t| eval(t).above, a.t, a.above, b.t, b.above);
            if t < end {
                if rising {
                    rise = Some(t);
                } else {
                    set = Some(t);
                }
            }
        }
    }

    let status = if rise.is_some() || set.is_some() {
        HorizonStatus::RisesOrSets
    } else if samples[0].above >= 0.0 {
        HorizonStatus::AlwaysAbove
    } else {
        HorizonStatus::AlwaysBelow
    };

    RiseTransitSet {
        rise,
        transit,
        set,
        status,
    }
}

/// Number of altitude samples per day before refinement (hourly).
const SAMPLES_PER_DAY: usize = 24;

/// Root-finding tolerance, in days (about 1 ms).
const ROOT_TOLERANCE: f64 = 1e-8;

struct Sample {
    t: f64,
    hour_angle: f64,
    /// Altitude minus the horizon altitude; >= 0 means above.
    above: f64,
}

/// Nutation quantities needed for an apparent position.
struct Nutation {
    long: f64,
    mn_oblq: f64,
    tru_oblq: f64,
}

impl Nutation {
    fn at(jde: f64) -> Nutation {
        let (long, oblq) = nutation::nutation(jde);
        let mn_oblq = ecliptic::mn_oblq_IAU(jde);
        Nutation {
            long,
            mn_oblq,
            tru_oblq: mn_oblq + oblq,
        }
    }

    /// Apparent sidereal time at Greenwich (Meeus, ch. 12): the mean
    /// sidereal time plus the nutation in right ascension.
    fn apparent_sidereal(&self, jd_ut: f64) -> f64 {
        time::mn_sidr(jd_ut) + self.long * self.tru_oblq.cos()
    }
}

/// Where the body's geocentric position comes from during a rise/set search.
enum Track {
    /// Computed afresh at every instant (the Moon).
    Direct(Body),
    /// Interpolated from positions at the start, middle and end of the day
    /// (the Sun, the planets and the stars, whose apparent motion over a day
    /// is very smooth).
    Interpolated {
        mid: f64,
        pos: [ApparentPosition; 3],
    },
}

impl Track {
    fn new(body: Body, jd_start: f64) -> Track {
        let pos_at = |jd_ut: f64| {
            let jde = jd_ut + delta_t_secs(jd_ut) / 86400.0;
            geocentric_apparent_position(body, jde)
        };
        match body {
            Body::Moon => Track::Direct(body),
            _ => Track::Interpolated {
                mid: jd_start + 0.5,
                pos: [
                    pos_at(jd_start),
                    pos_at(jd_start + 0.5),
                    pos_at(jd_start + 1.0),
                ],
            },
        }
    }

    /// The apparent geocentric position and the apparent sidereal time at
    /// Greenwich at an instant (UT).
    fn at(&self, jd_ut: f64) -> (ApparentPosition, f64) {
        let jde = jd_ut + delta_t_secs(jd_ut) / 86400.0;
        let nut = Nutation::at(jde);
        let pos = match self {
            Track::Direct(body) => geocentric(*body, jde, &nut),
            Track::Interpolated { mid, pos } => {
                let n = (jd_ut - mid) / 0.5;
                let interp =
                    |y1: f64, y2: f64, y3: f64| crate::interpol::three_values(y1, y2, y3, n);
                let a2 = pos[1].asc;
                ApparentPosition {
                    asc: norm_2pi(interp(
                        a2 + wrap_pi(pos[0].asc - a2),
                        a2,
                        a2 + wrap_pi(pos[2].asc - a2),
                    )),
                    dec: interp(pos[0].dec, pos[1].dec, pos[2].dec),
                    distance: if pos[1].distance.is_finite() {
                        interp(pos[0].distance, pos[1].distance, pos[2].distance)
                    } else {
                        pos[1].distance
                    },
                    eq_hz_parallax: interp(
                        pos[0].eq_hz_parallax,
                        pos[1].eq_hz_parallax,
                        pos[2].eq_hz_parallax,
                    ),
                }
            }
        };
        (pos, nut.apparent_sidereal(jd_ut))
    }
}

fn geocentric(body: Body, jde: f64, nut: &Nutation) -> ApparentPosition {
    let planet = match body {
        Body::Sun => return sun_position(jde, nut),
        Body::Moon => return moon_position(jde, nut),
        Body::Star { asc, dec } => return star_position(asc, dec, jde, nut),
        Body::Mercury => planet::Planet::Mercury,
        Body::Venus => planet::Planet::Venus,
        Body::Mars => planet::Planet::Mars,
        Body::Jupiter => planet::Planet::Jupiter,
        Body::Saturn => planet::Planet::Saturn,
        Body::Uranus => planet::Planet::Uranus,
        Body::Neptune => planet::Planet::Neptune,
    };
    planet_position(&planet, jde, nut)
}

/// Meeus, ch. 25, "higher accuracy".
fn sun_position(jde: f64, nut: &Nutation) -> ApparentPosition {
    let (l, b, r) = planet::heliocent_coords(&planet::Planet::Earth, jde);
    let (long, lat) = sun::ecl_coords_to_FK5(jde, l + PI, wrap_pi(-b));
    let long = long + nut.long + aberr::sol_aberr(r);

    from_ecliptic(long, lat, nut.tru_oblq, r, parallax::eq_hz_parallax(r))
}

/// Meeus, ch. 47: the series give the position referred to the mean equinox
/// of the date; adding the nutation in longitude makes it apparent.
fn moon_position(jde: f64, nut: &Nutation) -> ApparentPosition {
    let (ecl, dist_km) = lunar::geocent_ecl_pos(jde);

    from_ecliptic(
        ecl.long + nut.long,
        wrap_pi(ecl.lat),
        nut.tru_oblq,
        dist_km / AU_KM,
        lunar::eq_hz_parllx(dist_km),
    )
}

/// Meeus, ch. 33, first method, with the Earth also taken at `t − τ` so
/// that light-time and aberration are corrected together.
fn planet_position(planet: &planet::Planet, jde: f64, nut: &Nutation) -> ApparentPosition {
    let earth = |t| planet::heliocent_coords(&planet::Planet::Earth, t);
    let (l0, b0, r0) = earth(jde);

    let mut tau = 0.0;
    let mut dist = 0.0;
    for _ in 0..10 {
        let (l, b, r) = planet::heliocent_coords(planet, jde - tau);
        let (_, _, d, new_tau) = planet::geocent_geomet_ecl_coords(l0, b0, r0, l, b, r);
        let done = (new_tau - tau).abs() < 1e-9;
        tau = new_tau;
        dist = d;
        if done {
            break;
        }
    }

    let (l0, b0, r0) = earth(jde - tau);
    let (l, b, r) = planet::heliocent_coords(planet, jde - tau);
    let (long, lat, _, _) = planet::geocent_geomet_ecl_coords(l0, b0, r0, l, b, r);
    let (long, lat) = planet::ecl_coords_to_FK5(jde, long, lat);

    from_ecliptic(
        long + nut.long,
        lat,
        nut.tru_oblq,
        dist,
        parallax::eq_hz_parallax(dist),
    )
}

/// Meeus, chs. 21 and 23.
fn star_position(asc: f64, dec: f64, jde: f64, nut: &Nutation) -> ApparentPosition {
    let (asc, dec) = precess::precess_eq_coords(asc, dec, J2000, jde);
    let long = coords::ecl_long_frm_eq(asc, dec, nut.mn_oblq);
    let lat = coords::ecl_lat_frm_eq(asc, dec, nut.mn_oblq);

    let (long, lat) = annual_aberration(long, lat, jde);

    from_ecliptic(long + nut.long, lat, nut.tru_oblq, f64::INFINITY, 0.0)
}

/// Annual aberration of a star in ecliptic coordinates of the date.
///
/// Meeus' eq. 23.2 gives the first-order displacement caused by the Earth's
/// orbital velocity (including the terms in the eccentricity of the orbit);
/// it is applied here in vector form (the direction of `u + v`, with `v`
/// the Earth's velocity in units of the speed of light), which is
/// equivalent to 23.2 to first order but has no singularity at the
/// ecliptic poles.
fn annual_aberration(long: f64, lat: f64, jde: f64) -> (f64, f64) {
    let (sun_ecl, _) = sun::geocent_ecl_pos(jde);

    aberr::aberrate(long, lat, aberr::earth_vel_classical(sun_ecl.long, jde))
}

fn from_ecliptic(
    long: f64,
    lat: f64,
    tru_oblq: f64,
    distance: f64,
    eq_hz_parallax: f64,
) -> ApparentPosition {
    ApparentPosition {
        asc: norm_2pi(coords::asc_frm_ecl(long, lat, tru_oblq)),
        dec: coords::dec_frm_ecl(long, lat, tru_oblq),
        distance,
        eq_hz_parallax,
    }
}

/// Ratio of the topocentric to the geocentric distance (Meeus, eqs. 40.6
/// and 40.7).
fn topocentric_distance_ratio(pos: &ApparentPosition, observer: &Observer, gast: f64) -> f64 {
    let (rho_sin, rho_cos) = planet::earth::rho_sin_cos_phi(observer.lat, observer.height);
    let h = coords::hr_angl_frm_observer_long(gast, observer.long, pos.asc);
    let sin_par = pos.eq_hz_parallax.sin();
    let a = pos.dec.cos() * h.sin();
    let b = pos.dec.cos() * h.cos() - rho_cos * sin_par;
    let c = pos.dec.sin() - rho_sin * sin_par;

    (a * a + b * b + c * c).sqrt()
}

/// ΔT in seconds for a Julian Day, from [`time::delta_t_from_decimal_year`]
/// (the same model as [`time::delta_t`], evaluated at the instant itself so
/// that it varies continuously).
fn delta_t_secs(jd: f64) -> f64 {
    let year = 2000.0 + (jd - 2451544.5) / 365.2425;
    if !year.is_finite() {
        return 0.0;
    }
    time::delta_t_from_decimal_year(year)
}

/// Finds a zero of `f` between `a` and `b`, where `f(a)` and `f(b)` are on
/// opposite sides of zero (`f >= 0` counts as positive), by the Illinois
/// variant of the regula falsi.
fn find_root(
    mut f: impl FnMut(f64) -> f64,
    mut a: f64,
    mut fa: f64,
    mut b: f64,
    mut fb: f64,
) -> f64 {
    let positive = |v: f64| v >= 0.0;
    let mut last_side = 0;
    for _ in 0..200 {
        if (b - a).abs() < ROOT_TOLERANCE {
            break;
        }
        let mut c = (a * fb - b * fa) / (fb - fa);
        if !(c > a.min(b) && c < a.max(b)) {
            c = 0.5 * (a + b);
        }
        let fc = f(c);
        if positive(fc) == positive(fa) {
            a = c;
            fa = fc;
            if last_side == -1 {
                fb *= 0.5;
            }
            last_side = -1;
        } else {
            b = c;
            fb = fc;
            if last_side == 1 {
                fa *= 0.5;
            }
            last_side = 1;
        }
    }
    0.5 * (a + b)
}

/// Reduces an angle to `[−π, π)`.
#[inline]
fn wrap_pi(a: f64) -> f64 {
    (a + PI).rem_euclid(TAU) - PI
}

/// Reduces an angle to `[0, 2π)`.
#[inline]
fn norm_2pi(a: f64) -> f64 {
    let r = a.rem_euclid(TAU);
    if r >= TAU { 0.0 } else { r }
}

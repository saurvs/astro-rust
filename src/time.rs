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

//! Time for astronomy: Julian Day, calendar dates, ΔT and sidereal time.
//!
//! This module follows Meeus, *Astronomical Algorithms*, 2nd ed.:
//!
//! - chapter 7: Julian Day ([`julian_day`], [`Date::from_julian_day`]), day of
//!   the week, day of the year, Modified Julian Day, leap years;
//! - chapter 10: Dynamical Time and Universal Time ([`delta_t`] and the
//!   approximations of the book);
//! - chapter 12: sidereal time at Greenwich ([`mn_sidr`], [`apprnt_sidr`],
//!   [`mean_sidereal_time_0h`], [`local_sidereal_time`]).
//!
//! # Conventions
//!
//! - A Julian Day (JD) is a count of days from -4712 January 1, 12h. When
//!   the instant is measured in Universal Time it is a JD (UT); when it is
//!   measured in Dynamical (Terrestrial) Time it is a Julian Ephemeris Day
//!   (JDE). [`julian_ephemeris_day`] converts from one to the other with ΔT.
//! - Years are numbered astronomically: the year before 1 is 0, the one
//!   before that is -1, and so on (so 585 B.C. is the year -584).
//! - The calendar of a [`Date`] is explicit ([`CalType`]). The Julian
//!   calendar is extended indefinitely into the past and the Gregorian
//!   calendar may be used proleptically. Functions that pick a calendar
//!   themselves use the historical switch: 1582 October 4 (Julian) was
//!   followed by 1582 October 15 (Gregorian).
//! - Longitudes on the Earth are measured **positive westward** from
//!   Greenwich, as in Meeus.
//!
//! # Example
//!
//! ```
//! use astro::time::{self, CalType, Date, Month};
//!
//! // Launch of Sputnik 1: 1957 October 4.81 (Meeus, example 7.a)
//! let date = Date {
//!     year: 1957,
//!     month: Month::Oct,
//!     decimal_day: 4.81,
//!     cal_type: CalType::Gregorian,
//! };
//! let jd = time::julian_day(&date);
//! assert!((jd - 2436116.31).abs() < 1e-9);
//!
//! // ...and back (Meeus, example 7.c)
//! let back = Date::from_julian_day(jd).unwrap();
//! assert_eq!((back.year, back.month, back.cal_type), (1957, Month::Oct, CalType::Gregorian));
//! assert!((back.decimal_day - 4.81).abs() < 1e-6);
//!
//! // Mean sidereal time at Greenwich, 1987 April 10 at 0h UT (Meeus, example 12.a)
//! let theta0 = time::mn_sidr(2446895.5).to_degrees();
//! assert!((theta0 / 15.0 - (13.0 + 10.0 / 60.0 + 46.3668 / 3600.0)).abs() < 1e-7);
//! ```

use std::f64::consts::TAU;

/// JD of the first day of the Gregorian calendar, 1582 October 15 at 0h.
///
/// Instants before it are converted to dates in the Julian calendar by
/// [`Date::from_julian_day`] and [`date_frm_julian_day`].
pub const GREGORIAN_REFORM_JD: f64 = 2299160.5;

/// Ratio of the rate of mean sidereal time to the rate of mean solar time
/// (Universal Time), 1.00273790935.
///
/// A UT interval multiplied by this ratio gives the corresponding interval
/// of mean sidereal time (Meeus, chapter 12).
pub const SIDEREAL_RATE: f64 = 1.002_737_909_35;

/// Represents a calendar type
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CalType {
    /// Gregorian calendar
    Gregorian,
    /// Julian calendar
    Julian,
}

/// Represents a month in the Gregorian and Julian calendars
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Month {
    /// January
    Jan = 1,
    /// February
    Feb = 2,
    /// March
    Mar = 3,
    /// April
    Apr = 4,
    /// May
    May = 5,
    /// June
    June = 6,
    /// July
    July = 7,
    /// August
    Aug = 8,
    /// September
    Sept = 9,
    /// October
    Oct = 10,
    /// November
    Nov = 11,
    /// December
    Dec = 12,
}

impl Month {
    /// Returns the month with the given number (1 = January, ...,
    /// 12 = December), or `None` if `number` is not in `1..=12`.
    pub fn from_number(number: u8) -> Option<Month> {
        Some(match number {
            1 => Month::Jan,
            2 => Month::Feb,
            3 => Month::Mar,
            4 => Month::Apr,
            5 => Month::May,
            6 => Month::June,
            7 => Month::July,
            8 => Month::Aug,
            9 => Month::Sept,
            10 => Month::Oct,
            11 => Month::Nov,
            12 => Month::Dec,
            _ => return None,
        })
    }

    /// Returns the number of the month (1 = January, ..., 12 = December).
    pub fn number(self) -> u8 {
        self as u8
    }
}

/// Represents a date with year, month, decimal day and calendar type
///
/// The fraction of `decimal_day` is the fraction of the day elapsed since
/// 0h, so 4.81 is the 4th of the month at 19h26m24s. Functions that take a
/// `Date` treat the day arithmetically: a `decimal_day` of 0.5 is the last
/// day of the previous month at 12h, and 32.0 in January is February 1.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Date {
    /// Year, numbered astronomically (1 B.C. is 0, 2 B.C. is -1)
    pub year: i16,
    /// Month
    pub month: Month,
    /// Day of the month, with the fraction of the day
    ///
    /// Normally in the range `[1.0, 32.0)`.
    pub decimal_day: f64,
    /// Calendar type
    pub cal_type: CalType,
}

/// Represents a day of a month with hours, minutes and seconds
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DayOfMonth {
    /// Day of month
    ///
    /// range: *1 - 31*
    pub day: u8,
    /// Hour of day
    ///
    /// range: *0 - 23*
    pub hr: u8,
    /// Minute of hour
    ///
    /// range: *0 - 59*
    pub min: u8,
    /// Second of minute
    ///
    /// range: *0.0 - 60.0*
    pub sec: f64,
    /// Time zone *| in decimal hours*, positive east of Greenwich
    ///
    /// The local time is `time_zone` hours ahead of UT. Example: Pacific
    /// Standard Time is -8.0, Central European Time is +1.0.
    pub time_zone: f64,
}

/// Represents a day of the week
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Weekday {
    /// Sunday
    Sunday,
    /// Monday
    Monday,
    /// Tuesday
    Tuesday,
    /// Wednesday
    Wednesday,
    /// Thursday
    Thursday,
    /// Friday
    Friday,
    /// Saturday
    Saturday,
}

impl Date {
    /// Computes the calendar date corresponding to a Julian Day.
    ///
    /// The date is given in the Julian calendar before 1582 October 15
    /// (JD < [`GREGORIAN_REFORM_JD`]) and in the Gregorian calendar from then
    /// on. Negative Julian Days (dates before -4712 January 1.5) are
    /// supported.
    ///
    /// Returns `None` if `jd` is not finite or if the year does not fit in
    /// an `i16`.
    ///
    /// Meeus, chapter 7, "Calculation of the Calendar Date from the JD".
    ///
    /// # Arguments
    ///
    /// * `jd`: Julian Day (any time scale)
    ///
    /// # Example
    ///
    /// ```
    /// use astro::time::{CalType, Date, Month};
    ///
    /// // Meeus, exercise after example 7.c
    /// let d = Date::from_julian_day(1507900.13).unwrap();
    /// assert_eq!((d.year, d.month, d.cal_type), (-584, Month::May, CalType::Julian));
    /// assert!((d.decimal_day - 28.63).abs() < 1e-6);
    /// ```
    pub fn from_julian_day(jd: f64) -> Option<Date> {
        date_frm_jd(jd, None)
    }

    /// Computes the date corresponding to a Julian Day in a given calendar.
    ///
    /// Unlike [`Date::from_julian_day`], the calendar is not switched in
    /// 1582: the Gregorian calendar is used proleptically before 1582, and
    /// the Julian calendar after it, if so requested.
    ///
    /// Returns `None` if `jd` is not finite or if the year does not fit in
    /// an `i16`.
    ///
    /// # Arguments
    ///
    /// * `jd`      : Julian Day (any time scale)
    /// * `cal_type`: Calendar of the result
    pub fn from_julian_day_in(jd: f64, cal_type: CalType) -> Option<Date> {
        date_frm_jd(jd, Some(cal_type))
    }
}

fn date_frm_jd(jd: f64, cal_type: Option<CalType>) -> Option<Date> {
    // Far beyond the i16 year range; also keeps the integer casts exact.
    if !jd.is_finite() || jd.abs() > 1.0e9 {
        return None;
    }

    let jd = jd + 0.5;
    let z = jd.floor();
    let f = jd - z;
    let z = z as i64;

    let cal_type = cal_type.unwrap_or(if z < 2299161 {
        CalType::Julian
    } else {
        CalType::Gregorian
    });

    // Meeus' algorithm is stated for non-negative Julian Days. Earlier days
    // are moved forward by whole calendar cycles (4 Julian years = 1461 days,
    // 400 Gregorian years = 146097 days), and the year is moved back again.
    let (cycle_days, cycle_years) = match cal_type {
        CalType::Julian => (1461, 4),
        CalType::Gregorian => (146_097, 400),
    };
    let (z, year_shift) = if z < 0 {
        let k = -z / cycle_days + 1;
        (z + k * cycle_days, k * cycle_years)
    } else {
        (z, 0)
    };

    let a = match cal_type {
        CalType::Julian => z,
        CalType::Gregorian => {
            let alpha = ((z as f64 - 1867216.25) / 36524.25).floor() as i64;
            z + 1 + alpha - alpha.div_euclid(4)
        }
    };

    let b = a + 1524;
    let c = ((b as f64 - 122.1) / 365.25).floor() as i64;
    let d = (365.25 * c as f64).floor() as i64;
    let e = ((b - d) as f64 / 30.6001).floor() as i64;

    let day = (b - d - (30.6001 * e as f64).floor() as i64) as f64 + f;
    let month = if e < 14 { e - 1 } else { e - 13 };
    let year = if month > 2 { c - 4716 } else { c - 4715 } - year_shift;

    Some(Date {
        year: i16::try_from(year).ok()?,
        month: Month::from_number(u8::try_from(month).ok()?)?,
        decimal_day: day,
        cal_type,
    })
}

/// Computes the day of the week of a date
///
/// The weekday is that of the calendar day `date.decimal_day.floor()` in the
/// calendar `date.cal_type`. The week cycle was not interrupted by the
/// Gregorian reform: Thursday 1582 October 4 (Julian) was followed by Friday
/// 1582 October 15 (Gregorian).
///
/// Meeus, chapter 7, "Day of the week".
///
/// # Arguments
///
/// * `date`: A `Date`
pub fn weekday_frm_date(date: &Date) -> Weekday {
    let date_0h = Date {
        decimal_day: date.decimal_day.floor(),
        ..*date
    };
    weekday_from_julian_day(julian_day(&date_0h))
}

/// Computes the day of the week for a Julian Day
///
/// Returns the weekday of the calendar day (from 0h to 24h) that contains
/// the instant `jd`, so the result is in the same time scale as `jd`
/// (usually UT). For a non-finite `jd` the result is meaningless.
///
/// Meeus, chapter 7, "Day of the week".
///
/// # Arguments
///
/// * `jd`: Julian Day
pub fn weekday_from_julian_day(jd: f64) -> Weekday {
    // Days start at JD n.5, so jd + 1.5 runs from one integer to the next
    // over a day; JD 0.5 + k·7 days (e.g. -4712 January 2) is a Tuesday.
    match (jd + 1.5).floor().rem_euclid(7.0) as u8 {
        0 => Weekday::Sunday,
        1 => Weekday::Monday,
        2 => Weekday::Tuesday,
        3 => Weekday::Wednesday,
        4 => Weekday::Thursday,
        5 => Weekday::Friday,
        _ => Weekday::Saturday,
    }
}

/// Computes decimal day for a `DayOfMonth`
///
/// The result is the day of the month plus the fraction of the day, in UT:
/// the time zone offset is removed, so the result may fall outside the
/// range of the month (for example 0.9 or 31.1), which [`julian_day`]
/// handles correctly.
///
/// # Arguments
///
/// * `day_of_month`: A `DayOfMonth` struct
///
/// # Example
///
/// ```
/// use astro::time::{self, DayOfMonth};
///
/// // 19h26m24s UT on the 4th (Meeus, example 7.a)
/// let d = DayOfMonth { day: 4, hr: 19, min: 26, sec: 24.0, time_zone: 0.0 };
/// assert!((time::decimal_day(&d) - 4.81).abs() < 1e-12);
/// ```
pub fn decimal_day(day: &DayOfMonth) -> f64 {
    (day.day as f64)
        + (day.hr as f64) / 24.0
        + (day.min as f64) / (60.0 * 24.0)
        + day.sec / (60.0 * 60.0 * 24.0)
        - day.time_zone / 24.0
}

/// Computes decimal year for a `Date`
///
/// The result is `year + t / L`, where `t` is the time in days elapsed since
/// January 1 at 0h of `date.year` and `L` is the length of that year in
/// days (365 or 366) in the date's calendar. January 1.0 thus gives exactly
/// `year`, and the whole year maps onto `[year, year + 1)`.
///
/// # Arguments
///
/// * `date`: A `Date` struct
pub fn decimal_year(date: &Date) -> f64 {
    let days = days_in_year(date.year, date.cal_type) as f64;
    (date.year as f64) + (day_of_year(date) - 1.0) / days
}

/// Checks if a year is a leap year
///
/// In the Julian calendar every year divisible by 4 is a leap year. The
/// Gregorian calendar follows the same rule except for the century years
/// not divisible by 400 (1700, 1800, 1900, 2100, ...), which are common
/// years. Negative (astronomically numbered) years follow the same rules,
/// so the years 0, -4, -8, ... are leap years.
///
/// # Arguments
///
/// * `year`    : Year
/// * `cal_type`: `CalType` enum
pub fn is_leap_year(year: i16, cal_type: &CalType) -> bool {
    match *cal_type {
        CalType::Julian => year % 4 == 0,
        CalType::Gregorian => {
            if year % 100 == 0 {
                year % 400 == 0
            } else {
                year % 4 == 0
            }
        }
    }
}

fn days_in_year(year: i16, cal_type: CalType) -> u16 {
    if is_leap_year(year, &cal_type) {
        366
    } else {
        365
    }
}

/// Computes the day of the year (ordinal date) of a `Date`
///
/// Returns `N + f`, where `N` is the number of the day in the year (1 on
/// January 1, up to 365 or 366 on December 31) and `f` is the fraction of
/// the day in `date.decimal_day`. Use `.floor()` for the integer day
/// number. The day of the month is used arithmetically, so days outside the
/// month extend linearly.
///
/// Meeus, chapter 7, "Day of the Year".
///
/// # Arguments
///
/// * `date`: A `Date`
///
/// # Example
///
/// ```
/// use astro::time::{self, CalType, Date, Month};
///
/// // Meeus, example 7.f
/// let date = Date { year: 1978, month: Month::Nov, decimal_day: 14.0, cal_type: CalType::Gregorian };
/// assert_eq!(time::day_of_year(&date), 318.0);
/// ```
pub fn day_of_year(date: &Date) -> f64 {
    let k = if is_leap_year(date.year, &date.cal_type) {
        1
    } else {
        2
    };
    let m = date.month as i32;
    let n = (275 * m) / 9 - k * ((m + 9) / 12) - 30;

    n as f64 + date.decimal_day
}

/// Computes the date from a year and a day of the year
///
/// This is the inverse of [`day_of_year`]: the integer part of
/// `day_of_year` is the number of the day in the year (1 = January 1), and
/// its fraction is carried into the `decimal_day` of the result.
///
/// Returns `None` if `day_of_year` is not finite or is outside
/// `[1, L + 1)`, where `L` is the number of days in the year (365 or 366).
///
/// Meeus, chapter 7, "Day of the Year" (algorithm by A. Pouplier).
///
/// # Arguments
///
/// * `year`       : Year
/// * `day_of_year`: Day of the year, with fraction
/// * `cal_type`   : Calendar of `year` and of the result
pub fn date_from_day_of_year(year: i16, day_of_year: f64, cal_type: CalType) -> Option<Date> {
    let days = days_in_year(year, cal_type) as f64;
    if !(1.0..days + 1.0).contains(&day_of_year) {
        return None;
    }

    let n = day_of_year.floor() as i32;
    let k = if days == 366.0 { 1 } else { 2 };

    // INT(9 (K + N) / 275 + 0.98), in exact integer arithmetic
    let m = if n < 32 {
        1
    } else {
        (900 * (k + n) + 26950) / 27500
    };
    let d = n - (275 * m) / 9 + k * ((m + 9) / 12) + 30;

    Some(Date {
        year,
        month: Month::from_number(m as u8)?,
        decimal_day: d as f64 + (day_of_year - n as f64),
        cal_type,
    })
}

/// Computes Julian century for a Julian day
///
/// Returns the number of Julian centuries (of 36525 days) elapsed since the
/// epoch J2000.0 (JD 2451545.0); negative before it.
///
/// # Arguments
///
/// * `JD`: Julian (Ephemeris) day
#[inline]
pub fn julian_cent(JD: f64) -> f64 {
    (JD - 2451545.0) / 36525.0
}

/// Computes Julian millennium for a Julian day
///
/// Returns the number of Julian millennia (of 365250 days) elapsed since
/// the epoch J2000.0 (JD 2451545.0); negative before it.
///
/// # Arguments
///
/// * `JD`: Julian (Ephemeris) day
#[inline]
pub fn julian_mill(JD: f64) -> f64 {
    (JD - 2451545.0) / 365250.0
}

/// Computes Julian day from a `Date`
///
/// Valid for all dates, in the Julian or the (proleptic) Gregorian
/// calendar as given by `date.cal_type`, including negative years and
/// dates that give a negative Julian Day. The calendar is not switched
/// automatically in 1582. The day of the month is used arithmetically, so a
/// `decimal_day` outside the month (as produced by [`decimal_day`] with a
/// time zone) still gives the right instant.
///
/// The time scale of the result is that of the date: a date in UT gives a
/// JD (UT), a date in Dynamical Time gives a JDE.
///
/// Meeus, chapter 7, eq. 7.1.
///
/// # Arguments
///
/// `date`: A `Date`
pub fn julian_day(date: &Date) -> f64 {
    let month = date.month as u8;
    let (y, m) = if month == 1 || month == 2 {
        ((date.year as f64) - 1.0, (month + 12) as f64)
    } else {
        (date.year as f64, month as f64)
    };

    let a = (y / 100.0).floor();
    let b = match date.cal_type {
        CalType::Gregorian => 2.0 - a + (a / 4.0).floor(),
        CalType::Julian => 0.0,
    };

    (365.25 * (y + 4716.0)).floor() + (30.6001 * (m + 1.0)).floor() + date.decimal_day + b - 1524.5
}

/// Computes the Julian Day of January 0.0 of a year
///
/// January 0.0 is the same instant as December 31.0 (0h) of the previous
/// year, so the day of the year of an instant is its JD minus this value.
///
/// Meeus, chapter 7 (JD₀).
///
/// # Arguments
///
/// * `year`    : Year
/// * `cal_type`: Calendar of `year`
pub fn julian_day_of_jan0(year: i16, cal_type: CalType) -> f64 {
    julian_day(&Date {
        year,
        month: Month::Jan,
        decimal_day: 0.0,
        cal_type,
    })
}

/// Computes the number of days from one date to another
///
/// Returns `JD(to) - JD(from)`, negative if `to` is earlier than `from`.
/// The two dates may be in different calendars.
///
/// Meeus, chapter 7, "Time interval in days".
///
/// # Arguments
///
/// * `from`: Start date
/// * `to`  : End date
pub fn days_between(from: &Date, to: &Date) -> f64 {
    julian_day(to) - julian_day(from)
}

/// Computes the Modified Julian Day for a Julian Day
///
/// MJD = JD - 2400000.5. Unlike the JD, the MJD starts at midnight; MJD 0.0
/// is 1858 November 17 at 0h.
///
/// Meeus, chapter 7.
///
/// # Arguments
///
/// * `jd`: Julian Day
#[inline]
pub fn modified_julian_day(jd: f64) -> f64 {
    jd - 2400000.5
}

/// Computes the Julian Day for a Modified Julian Day
///
/// This is the inverse of [`modified_julian_day`].
///
/// # Arguments
///
/// * `mjd`: Modified Julian Day
#[inline]
pub fn julian_day_from_modified(mjd: f64) -> f64 {
    mjd + 2400000.5
}

/// Computes the Julian Ephemeris day
///
/// JDE = JD + ΔT, with ΔT = TD - UT converted from seconds to days.
///
/// Meeus, chapter 10.
///
/// # Arguments
///
/// * `JD`     : Julian day (UT)
/// * `delta_t`: ΔT *| in seconds* (see [`delta_t`])
#[inline]
pub fn julian_ephemeris_day(JD: f64, delta_t: f64) -> f64 {
    delta_t / 86400.0 + JD
}

/// Computes the Julian Day (UT) for a Julian Ephemeris Day
///
/// JD = JDE - ΔT, with ΔT = TD - UT converted from seconds to days. This is
/// the inverse of [`julian_ephemeris_day`].
///
/// Meeus, chapter 10 (UT = TD - ΔT).
///
/// # Arguments
///
/// * `jde`    : Julian Ephemeris Day (Dynamical Time)
/// * `delta_t`: ΔT *| in seconds* (see [`delta_t`])
#[inline]
pub fn julian_day_from_ephemeris_day(jde: f64, delta_t: f64) -> f64 {
    jde - delta_t / 86400.0
}

/// Computes a year, month and decimal day equivalent to a given Julian day
///
/// The date is in the Julian calendar before 1582 October 15 (JD
/// 2299160.5) and in the Gregorian calendar from then on. Negative Julian
/// Days are supported. [`Date::from_julian_day`] gives the same result as a
/// [`Date`].
///
/// # Returns
///
/// `Ok((year, month, decimal_day))`
///
/// * `year`       : Year
/// * `month`      : Month (1 = January, ..., 12 = December)
/// * `decimal_day`: Decimal day
///
/// `Err` if `JD` is not finite, or if the year does not fit in an `i16`.
///
/// Meeus, chapter 7.
///
/// # Arguments
///
/// `JD`: Julian Day
pub fn date_frm_julian_day(JD: f64) -> Result<(i16, u8, f64), &'static str> {
    match Date::from_julian_day(JD) {
        Some(date) => Ok((date.year, date.month as u8, date.decimal_day)),
        None => Err("time::date_frm_julian_day(): JD is not finite or its year is out of range"),
    }
}

/// Normalises an angle in radians to `[0, 2π)`.
fn limit_to_tau(angle: f64) -> f64 {
    let a = angle.rem_euclid(TAU);
    // rem_euclid can round up to exactly 2π for tiny negative inputs
    if a >= TAU { 0.0 } else { a }
}

/// Normalises an angle in degrees to `[0, 360)`.
fn limit_to_360(angle: f64) -> f64 {
    let a = angle.rem_euclid(360.0);
    if a >= 360.0 { 0.0 } else { a }
}

/// Computes the equation of the equinoxes (nutation in right ascension)
///
/// Returns Δψ cos ε *| in radians*, the difference between apparent and
/// mean sidereal time. Divide by 15 after converting to arcseconds to get
/// seconds of time.
///
/// Meeus, chapter 12.
///
/// # Arguments
///
/// * `nut_in_long`: Nutation in longitude Δψ *| in radians*
/// * `true_oblq`  : True obliquity of the ecliptic ε *| in radians*
#[inline]
pub fn equation_of_equinoxes(nut_in_long: f64, true_oblq: f64) -> f64 {
    nut_in_long * true_oblq.cos()
}

/// Computes apparent sidereal time from the mean sidereal time
///
/// Adds the equation of the equinoxes (Δψ cos ε) to the mean sidereal time.
///
/// # Returns
///
/// * `apprnt_sidr`: Apparent sidereal time *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `mn_sidr`    : Mean sidereal time *| in radians*
/// * `nut_in_long`: Nutatation in longitude *| in radians*
/// * `true_oblq`  : True obliquity of the ecliptic *| in radians*
pub fn apprnt_sidr(mn_sidr: f64, nut_in_long: f64, true_oblq: f64) -> f64 {
    limit_to_tau(mn_sidr + equation_of_equinoxes(nut_in_long, true_oblq))
}

/// Computes apparent sidereal time at Greenwich for a Julian day
///
/// The nutation is from [`crate::nutation::nutation`] and the true obliquity
/// is J. Laskar's mean obliquity ([`crate::ecliptic::mn_oblq_laskar`]) plus
/// the nutation in obliquity. As in Meeus, the nutation is evaluated at the
/// UT instant: the difference with the TD instant is negligible here.
///
/// # Returns
///
/// * `apprnt_sidr`: Apparent sidereal time *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `jd`: Julian day (UT)
///
/// # Example
///
/// ```
/// use astro::time;
///
/// // 1987 April 10, 0h UT: 13h10m46.1351s (Meeus, example 12.a)
/// let hours = time::apparent_sidereal_time(2446895.5).to_degrees() / 15.0;
/// assert!((hours - (13.0 + 10.0 / 60.0 + 46.1351 / 3600.0)).abs() < 1e-4 / 3600.0);
/// ```
pub fn apparent_sidereal_time(jd: f64) -> f64 {
    let (nut_in_long, nut_in_oblq) = crate::nutation::nutation(jd);
    let true_oblq = crate::ecliptic::mn_oblq_laskar(jd) + nut_in_oblq;
    apprnt_sidr(mn_sidr(jd), nut_in_long, true_oblq)
}

/// Computes apparent sidereal time for a Julian day
///
/// This macro is equivalent to [`time::apparent_sidereal_time`]; it uses
/// J. Laskar's formula for the mean obliquity of the ecliptic.
///
/// [`time::apparent_sidereal_time`]: crate::time::apparent_sidereal_time
///
/// # Returns
///
/// * `apprnt_sidr`: Apparent sidereal time *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `$JD`: Julian day (UT)
#[macro_export]
macro_rules! apprnt_sidr {
    ($JD: expr_2021) => {{ $crate::time::apparent_sidereal_time($JD) }};
}

/// Computes mean sidereal time for a Julian day
///
/// This is the Greenwich hour angle of the mean vernal equinox, from the
/// expression adopted by the IAU in 1982, valid for any instant (Meeus,
/// eq. 12.4).
///
/// # Returns
///
/// * `mn_sidr`: Mean sidereal time *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `JD`: Julian day (UT)
pub fn mn_sidr(JD: f64) -> f64 {
    let JC = julian_cent(JD);

    limit_to_360(
        280.46061837
            + 360.98564736629 * (JD - 2451545.0)
            + JC * JC * (0.000387933 - JC / 38710000.0),
    )
    .to_radians()
}

/// Computes mean sidereal time at Greenwich at 0h UT
///
/// Returns Θ₀, the mean sidereal time at 0h UT of the UT calendar day that
/// contains `jd`, from Meeus' eq. 12.3 (which is only valid at 0h UT). The
/// mean sidereal time at an instant `UT` hours later is Θ₀ plus
/// [`SIDEREAL_RATE`] × UT, which agrees with [`mn_sidr`].
///
/// # Returns
///
/// * Mean sidereal time at 0h UT *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `jd`: Julian day (UT), any instant of the day
pub fn mean_sidereal_time_0h(jd: f64) -> f64 {
    let jd0 = (jd - 0.5).floor() + 0.5;
    let t = julian_cent(jd0);

    limit_to_360(100.46061837 + 36000.770053608 * t + t * t * (0.000387933 - t / 38710000.0))
        .to_radians()
}

/// Computes local sidereal time from the sidereal time at Greenwich
///
/// θ = θ₀ - L, where L is the observer's longitude, positive west of
/// Greenwich. Works with mean and apparent sidereal time alike.
///
/// # Returns
///
/// * Local sidereal time *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `greenwich_sidr`: Sidereal time at Greenwich *| in radians*
/// * `geog_long`     : Geographical longitude of the observer *| in
///   radians*, **positive west**, negative east of Greenwich
pub fn local_sidereal_time(greenwich_sidr: f64, geog_long: f64) -> f64 {
    limit_to_tau(greenwich_sidr - geog_long)
}

/// Evaluates a polynomial with coefficients in ascending powers of `x`.
fn poly(x: f64, coeffs: &[f64]) -> f64 {
    coeffs.iter().rev().fold(0.0, |acc, &c| acc * x + c)
}

/// Computes an approximate value of ΔT for a given year and month
///
/// ΔT = TD - UT, in seconds. This uses the polynomial expressions of
/// Espenak and Meeus (NASA, *Five Millennium Canon of Solar Eclipses*,
/// 2006), fitted to the values of Morrison & Stephenson (2004) and of the
/// Astronomical Almanac. They are more accurate, and cover a wider time
/// range, than the approximations given in the *Meeus* book (which are
/// available as [`delta_t_chapront_francou`], [`delta_t_1800_1997`],
/// [`delta_t_1800_1899`] and [`delta_t_1900_1997`]).
///
/// The polynomials are evaluated at the middle of the month, `y = year +
/// (month - 0.5) / 12`; see [`delta_t_from_decimal_year`].
///
/// # Arguments
///
/// * `year` : Year (astronomical numbering)
/// * `month`: Month *range: 1 - 12*
pub fn delta_t(year: i32, month: u8) -> f64 {
    delta_t_from_decimal_year(year as f64 + (month as f64 - 0.5) / 12.0)
}

/// Computes an approximate value of ΔT for a decimal year
///
/// ΔT = TD - UT, in seconds, from the Espenak & Meeus (NASA) polynomials
/// (see [`delta_t`]). The pieces join up to within 0.3 s.
///
/// # Validity
///
/// - The polynomials are meant for -1999 to +3000. Before -500 and after
///   +2150 they reduce to the long-term parabola `-20 + 32 u²`, with
///   `u = (y - 1820) / 100`, which is used for any year outside that range
///   too; such values are rough extrapolations.
/// - The uncertainty grows from below a second in the 20th century to
///   minutes in the Middle Ages and hours before 1000 B.C.
/// - After 2005 the values are predictions; by the 2020s they exceed the
///   observed ΔT by several seconds. Use observed values (IERS) for
///   precise recent work.
/// - The values assume a secular acceleration of the Moon of -26″/cy².
///
/// # Arguments
///
/// * `year`: Decimal year, e.g. 1990.5 for the middle of 1990
pub fn delta_t_from_decimal_year(year: f64) -> f64 {
    let y = year;
    let long_term = |y: f64| {
        let u = (y - 1820.0) / 100.0;
        -20.0 + 32.0 * u * u
    };

    if y < -500.0 {
        long_term(y)
    } else if y < 500.0 {
        poly(
            y / 100.0,
            &[
                10583.6,
                -1014.41,
                33.78311,
                -5.952053,
                -0.1798452,
                0.022174192,
                0.0090316521,
            ],
        )
    } else if y < 1600.0 {
        poly(
            (y - 1000.0) / 100.0,
            &[
                1574.2,
                -556.01,
                71.23472,
                0.319781,
                -0.8503463,
                -0.005050998,
                0.0083572073,
            ],
        )
    } else if y < 1700.0 {
        poly(y - 1600.0, &[120.0, -0.9808, -0.01532, 1.0 / 7129.0])
    } else if y < 1800.0 {
        poly(
            y - 1700.0,
            &[8.83, 0.1603, -0.0059285, 0.00013336, -1.0 / 1174000.0],
        )
    } else if y < 1860.0 {
        poly(
            y - 1800.0,
            &[
                13.72,
                -0.332447,
                0.0068612,
                0.0041116,
                -0.00037436,
                0.0000121272,
                -0.0000001699,
                0.000000000875,
            ],
        )
    } else if y < 1900.0 {
        poly(
            y - 1860.0,
            &[
                7.62,
                0.5737,
                -0.251754,
                0.01680668,
                -0.0004473624,
                1.0 / 233174.0,
            ],
        )
    } else if y < 1920.0 {
        poly(
            y - 1900.0,
            &[-2.79, 1.494119, -0.0598939, 0.0061966, -0.000197],
        )
    } else if y < 1941.0 {
        poly(y - 1920.0, &[21.20, 0.84493, -0.076100, 0.0020936])
    } else if y < 1961.0 {
        poly(y - 1950.0, &[29.07, 0.407, -1.0 / 233.0, 1.0 / 2547.0])
    } else if y < 1986.0 {
        poly(y - 1975.0, &[45.45, 1.067, -1.0 / 260.0, -1.0 / 718.0])
    } else if y < 2005.0 {
        poly(
            y - 2000.0,
            &[
                63.86,
                0.3345,
                -0.060374,
                0.0017275,
                0.000651814,
                0.00002373599,
            ],
        )
    } else if y < 2050.0 {
        poly(y - 2000.0, &[62.92, 0.32217, 0.005589])
    } else if y < 2150.0 {
        long_term(y) - 0.5628 * (2150.0 - y)
    } else {
        // also reached for NaN, which propagates
        long_term(y)
    }
}

/// Computes ΔT from the expressions of Chapront, Chapront-Touzé and Francou
///
/// ΔT = TD - UT, in seconds, from Meeus' eq. 10.1 (before 948) and eq.
/// 10.2 (948 to 1600, and after 2000). Between 2000 and 2100 the correction
/// `0.37 (year - 2100)` is added, as Meeus advises, to avoid a jump at 2000.
///
/// These are rough approximations: back to 4000 B.C. the error can reach
/// two hours. [`delta_t`] is more accurate and should normally be
/// preferred.
///
/// Returns `None` for 1600 < `year` < 2000, where the expressions are not
/// meant to be used (Meeus gives observed values there; use [`delta_t`]).
///
/// Meeus, chapter 10.
///
/// # Arguments
///
/// * `year`: Decimal year
pub fn delta_t_chapront_francou(year: f64) -> Option<f64> {
    let t = (year - 2000.0) / 100.0;
    if year < 948.0 {
        Some(2177.0 + 497.0 * t + 44.1 * t * t)
    } else if year <= 1600.0 || year >= 2000.0 {
        let mut dt = 102.0 + 102.0 * t + 25.3 * t * t;
        if (2000.0..=2100.0).contains(&year) {
            dt += 0.37 * (year - 2100.0);
        }
        Some(dt)
    } else {
        None
    }
}

/// Computes ΔT for 1800 to 1997 from Meeus' single polynomial
///
/// ΔT = TD - UT, in seconds. The polynomial in θ = (year - 1900) / 100
/// (Julian centuries from 1900.0) represents the observed values with a
/// maximum error of 2.3 s.
///
/// Returns `None` outside 1800.0 ≤ `year` < 1998.0: the expression diverges
/// quickly outside its range.
///
/// Meeus, chapter 10.
///
/// # Arguments
///
/// * `year`: Decimal year
pub fn delta_t_1800_1997(year: f64) -> Option<f64> {
    if !(1800.0..1998.0).contains(&year) {
        return None;
    }
    Some(poly(
        (year - 1900.0) / 100.0,
        &[
            -1.02, 91.02, 265.90, -839.16, -1545.20, 3603.62, 4385.98, -6993.23, -6090.04, 6298.12,
            4102.86, -2137.64, -1081.51,
        ],
    ))
}

/// Computes ΔT for 1800 to 1899 from Meeus' polynomial
///
/// ΔT = TD - UT, in seconds, with a maximum error of 0.9 s. The polynomial
/// is in θ = (year - 1900) / 100 (Julian centuries from 1900.0).
///
/// Returns `None` outside 1800.0 ≤ `year` < 1900.0: the expression diverges
/// quickly outside its range (it gives about 70000 s for 1945).
///
/// Meeus, chapter 10.
///
/// # Arguments
///
/// * `year`: Decimal year
pub fn delta_t_1800_1899(year: f64) -> Option<f64> {
    if !(1800.0..1900.0).contains(&year) {
        return None;
    }
    Some(poly(
        (year - 1900.0) / 100.0,
        &[
            -2.50, 228.95, 5218.61, 56282.84, 324011.78, 1061660.75, 2087298.89, 2513807.78,
            1818961.41, 727058.63, 123563.95,
        ],
    ))
}

/// Computes ΔT for 1900 to 1997 from Meeus' polynomial
///
/// ΔT = TD - UT, in seconds, with a maximum error of 0.9 s. The polynomial
/// is in θ = (year - 1900) / 100 (Julian centuries from 1900.0).
///
/// Returns `None` outside 1900.0 ≤ `year` < 1998.0: the expression diverges
/// quickly outside its range.
///
/// Meeus, chapter 10.
///
/// # Arguments
///
/// * `year`: Decimal year
pub fn delta_t_1900_1997(year: f64) -> Option<f64> {
    if !(1900.0..1998.0).contains(&year) {
        return None;
    }
    Some(poly(
        (year - 1900.0) / 100.0,
        &[
            -2.44, 87.24, 815.20, -2637.80, -18756.33, 124906.15, -303191.19, 372919.88,
            -232424.66, 58353.42,
        ],
    ))
}

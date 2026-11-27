//! Calendars: the date of Easter and the Jewish and Islamic calendars.
//!
//! This module follows Meeus, *Astronomical Algorithms*, 2nd ed.:
//!
//! - chapter 8: the date of Easter Sunday in the Gregorian
//!   ([`easter_gregorian`]) and Julian ([`easter_julian`]) calendars;
//! - chapter 9: the Jewish calendar (date of Pesach, Jewish New Year, type
//!   and length of a Jewish year) and the Islamic (Moslem) calendar
//!   ([`date_from_islamic`], [`islamic_from_date`]).
//!
//! Dates in the Julian and Gregorian calendars are [`time::Date`]s. Years
//! are numbered astronomically (the year before 1 is 0).
//!
//! # Example
//!
//! ```
//! use astro::calendar;
//! use astro::time::{CalType, Month};
//!
//! // Easter Sunday 2000 was on April 23 (Meeus, chapter 8)
//! let easter = calendar::easter_gregorian(2000).unwrap();
//! assert_eq!((easter.month, easter.decimal_day), (Month::Apr, 23.0));
//!
//! // Pesach (15 Nisan) in 1990 was on April 10 (Meeus, example 9.a)
//! let pesach = calendar::pesach(1990);
//! assert_eq!((pesach.month, pesach.decimal_day), (Month::Apr, 10.0));
//!
//! // 1 Muharram 1421 was 2000 April 6 (Meeus, example 9.b)
//! let date = calendar::date_from_islamic(&calendar::IslamicDate { year: 1421, month: 1, day: 1 })
//!     .unwrap();
//! assert_eq!((date.year, date.month, date.decimal_day), (2000, Month::Apr, 6.0));
//! assert_eq!(date.cal_type, CalType::Gregorian);
//! ```

use crate::time::{self, CalType, Date, Month};

/// JD of the start of the Islamic era: 1 Muharram A.H. 1, which is 622 July
/// 16 (Julian calendar) at 0h.
pub const ISLAMIC_EPOCH_JD: f64 = 1948439.5;

/// Builds the date `day` March of `year`, normalising days past the end of
/// March into the following months.
fn date_in_march(year: i16, day: i64, cal_type: CalType) -> Date {
    let march = Date {
        year,
        month: Month::Mar,
        decimal_day: day as f64,
        cal_type,
    };
    Date::from_julian_day_in(time::julian_day(&march), cal_type).unwrap_or(march)
}

/// Computes the date of Easter Sunday in the Gregorian calendar
///
/// Uses the method published by Spencer Jones (from Butcher's
/// *Ecclesiastical Calendar*, 1876). It needs no special cases, unlike
/// Gauss' rule, and works for any Gregorian year. Easter Sunday falls
/// between March 22 and April 25.
///
/// Returns `None` for years before 1583, the first full year of the
/// Gregorian calendar.
///
/// Meeus, chapter 8.
///
/// # Arguments
///
/// * `year`: Year in the Gregorian calendar
pub fn easter_gregorian(year: i16) -> Option<Date> {
    if year < 1583 {
        return None;
    }
    let x = year as i32;

    let a = x % 19;
    let (b, c) = (x / 100, x % 100);
    let (d, e) = (b / 4, b % 4);
    let f = (b + 8) / 25;
    let g = (b - f + 1) / 3;
    let h = (19 * a + b - d - g + 15) % 30;
    let (i, k) = (c / 4, c % 4);
    let l = (32 + 2 * e + 2 * i - h - k) % 7;
    let m = (a + 11 * h + 22 * l) / 451;
    let n = (h + l - 7 * m + 114) / 31;
    let p = (h + l - 7 * m + 114) % 31;

    Some(Date {
        year,
        month: if n == 3 { Month::Mar } else { Month::Apr },
        decimal_day: (p + 1) as f64,
        cal_type: CalType::Gregorian,
    })
}

/// Computes the date of Easter Sunday in the Julian calendar
///
/// The result is a date in the Julian calendar, between March 22 and April
/// 25 (Julian). The Julian Easter repeats every 532 years. The rule is
/// applied to any year, although it was only fixed in the 4th century.
///
/// The Easter of the Eastern churches is this date; convert it to the
/// Gregorian calendar with [`time::julian_day`] and
/// [`Date::from_julian_day_in`].
///
/// Meeus, chapter 8.
///
/// # Arguments
///
/// * `year`: Year in the Julian calendar
pub fn easter_julian(year: i16) -> Date {
    let x = year as i32;

    let a = x.rem_euclid(4);
    let b = x.rem_euclid(7);
    let c = x.rem_euclid(19);
    let d = (19 * c + 15) % 30;
    let e = (2 * a + 4 * b - d + 34) % 7;
    let f = (d + e + 114) / 31;
    let g = (d + e + 114) % 31;

    Date {
        year,
        month: if f == 3 { Month::Mar } else { Month::Apr },
        decimal_day: (g + 1) as f64,
        cal_type: CalType::Julian,
    }
}

/// Gauss' formula for Pesach; returns `(Q, j, r, D)`, where `D` is the day
/// of March (it may exceed 31).
fn pesach_terms(x: i64) -> (f64, i64, f64, i64) {
    let c = x.div_euclid(100);
    let s = if x >= 1583 {
        (3 * c - 5).div_euclid(4)
    } else {
        0
    };
    let a = (12 * x + 12).rem_euclid(19);
    let b = x.rem_euclid(4);

    let q = -1.904412361576 + 1.554241796621 * a as f64 + 0.25 * b as f64
        - 0.003177794022 * x as f64
        + s as f64;
    let int_q = q.floor();
    let j = (int_q as i64 + 3 * x + 5 * b + 2 - s).rem_euclid(7);
    let r = q - int_q;

    let d = int_q as i64
        + match j {
            2 | 4 | 6 => 23,
            1 if a > 6 && r >= 0.632870370 => 24,
            0 if a > 11 && r >= 0.897723765 => 23,
            _ => 22,
        };

    (q, j, r, d)
}

/// Computes the date of Pesach (15 Nisan) in a Julian or Gregorian year
///
/// Uses Gauss' formula as given by Meeus. The result is in the Julian
/// calendar for years before 1583 and in the Gregorian calendar from 1583
/// on. Pesach falls in the Jewish year `year + 3760`. Like every Jewish
/// day, 15 Nisan begins at sunset on the evening before the date returned.
///
/// The present Jewish calendar dates from the 4th century; for earlier
/// years the result is proleptic.
///
/// Meeus, chapter 9.
///
/// # Arguments
///
/// * `year`: Year in the Julian or Gregorian calendar (astronomical
///   numbering)
pub fn pesach(year: i16) -> Date {
    let cal_type = if year < 1583 {
        CalType::Julian
    } else {
        CalType::Gregorian
    };
    date_in_march(year, pesach_terms(year as i64).3, cal_type)
}

/// Returns the Julian Day of 1 Tishri of a Jewish year, if the civil years
/// involved fit in an `i16`.
fn jewish_new_year_jd(jewish_year: i16) -> Option<f64> {
    // 1 Tishri of year A falls 163 days after Pesach of the civil year
    // A - 3761 (which is 15 Nisan of the Jewish year A - 1).
    let civil = i16::try_from(jewish_year as i32 - 3761).ok()?;
    Some(time::julian_day(&pesach(civil)) + 163.0)
}

/// Computes the date of the Jewish New Year (1 Tishri)
///
/// Returns the date of 1 Tishri of the Jewish year `jewish_year`, which
/// falls in the civil year `jewish_year - 3761`. The date is in the Julian
/// calendar before 1582 October 15 and in the Gregorian calendar from then
/// on. It is in September or October for every Gregorian year up to about
/// 8100 (the Jewish year is slightly longer than the Gregorian one, so the
/// date drifts later, into November after that); in the Julian calendar it
/// can be as early as August 25 (in many years from the 1st to the 16th
/// century).
///
/// Returns `None` if the civil year does not fit in an `i16`.
///
/// Meeus, chapter 9.
///
/// # Arguments
///
/// * `jewish_year`: Year of the Jewish calendar (Anno Mundi)
pub fn jewish_new_year(jewish_year: i16) -> Option<Date> {
    Date::from_julian_day(jewish_new_year_jd(jewish_year)?)
}

/// Checks if a Jewish year is embolismic (a leap year of 13 months)
///
/// The years 3, 6, 8, 11, 14, 17 and 19 of each 19-year cycle are
/// embolismic.
///
/// Meeus, chapter 9.
///
/// # Arguments
///
/// * `jewish_year`: Year of the Jewish calendar (Anno Mundi)
pub fn is_jewish_leap_year(jewish_year: i16) -> bool {
    matches!(
        (jewish_year as i32).rem_euclid(19),
        0 | 3 | 6 | 8 | 11 | 14 | 17
    )
}

/// Represents the three lengths a Jewish year can have
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum JewishYearKind {
    /// Deficient year: 353 days (common) or 383 days (embolismic)
    Deficient,
    /// Regular year: 354 days (common) or 384 days (embolismic)
    Regular,
    /// Complete year: 355 days (common) or 385 days (embolismic)
    Complete,
}

/// Represents a year of the Jewish calendar
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JewishYear {
    /// Date of 1 Tishri, the first day of the year, in the Julian (before
    /// 1582 October 15) or Gregorian calendar
    pub new_year: Date,
    /// Number of days in the year: 353, 354, 355, 383, 384 or 385
    pub days: u16,
    /// Whether the year is embolismic (13 months instead of 12)
    pub embolismic: bool,
    /// Deficient, regular or complete
    pub kind: JewishYearKind,
}

/// Computes the beginning, length and type of a Jewish year
///
/// The length is the number of days from 1 Tishri of `jewish_year` to 1
/// Tishri of the next year.
///
/// Returns `None` if the civil years involved do not fit in an `i16`.
///
/// Meeus, chapter 9.
///
/// # Arguments
///
/// * `jewish_year`: Year of the Jewish calendar (Anno Mundi)
///
/// # Example
///
/// ```
/// use astro::calendar::{self, JewishYearKind};
/// use astro::time::Month;
///
/// // Meeus, example 9.a: the year 5751 began on 1990 September 20 and has 354 days
/// let y = calendar::jewish_year(5751).unwrap();
/// assert_eq!((y.new_year.year, y.new_year.month, y.new_year.decimal_day), (1990, Month::Sept, 20.0));
/// assert_eq!((y.days, y.embolismic, y.kind), (354, false, JewishYearKind::Regular));
/// ```
pub fn jewish_year(jewish_year: i16) -> Option<JewishYear> {
    let start = jewish_new_year_jd(jewish_year)?;
    let end = jewish_new_year_jd(jewish_year.checked_add(1)?)?;
    let days = (end - start) as u16;
    let kind = match days {
        353 | 383 => JewishYearKind::Deficient,
        354 | 384 => JewishYearKind::Regular,
        355 | 385 => JewishYearKind::Complete,
        _ => return None,
    };
    Some(JewishYear {
        new_year: Date::from_julian_day(start)?,
        days,
        embolismic: is_jewish_leap_year(jewish_year),
        kind,
    })
}

/// Represents a date in the Islamic (Moslem) calendar
///
/// The tabular (arithmetical) Islamic calendar of Meeus, chapter 9: a
/// purely lunar calendar of 12 months of alternately 30 and 29 days, the
/// last month having 30 days in leap years. Years are counted from the
/// Hegira (A.H.).
#[doc(alias = "moslem")]
#[doc(alias = "hijri")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct IslamicDate {
    /// Year (Anno Hegirae), from 1
    pub year: i16,
    /// Month, 1 (Muharram) to 12 (Dhu l-Hijja)
    pub month: u8,
    /// Day of the month, 1 to 29 or 30
    pub day: u8,
}

/// Checks if an Islamic year is a leap year of 355 days
///
/// The years 2, 5, 7, 10, 13, 16, 18, 21, 24, 26 and 29 of each 30-year
/// cycle are leap years; the others have 354 days.
///
/// Meeus, chapter 9.
///
/// # Arguments
///
/// * `year`: Year of the Islamic calendar (A.H.)
pub fn is_islamic_leap_year(year: i16) -> bool {
    islamic_year_days(year as i64) == 355
}

fn islamic_year_days(year: i64) -> i64 {
    if (11 * year.rem_euclid(30) + 3) % 30 > 18 {
        355
    } else {
        354
    }
}

fn islamic_month_days(year: i64, month: i64) -> i64 {
    if month % 2 == 1 || (month == 12 && islamic_year_days(year) == 355) {
        30
    } else {
        29
    }
}

/// Meeus' conversion of an Islamic date to the year `X` and the day number
/// `J` within that year of the Julian calendar.
fn islamic_to_julian_year_and_day(h: i64, m: i64, d: i64) -> (i64, i64) {
    let n = d + (29.5001 * (m - 1) as f64 + 0.99).floor() as i64;
    let q = h.div_euclid(30);
    let r = h.rem_euclid(30);
    let a = (11 * r + 3) / 30;
    let w = 404 * q + 354 * r + 208 + a;
    let q1 = w.div_euclid(1461);
    let q2 = w.rem_euclid(1461);
    let g = 621 + 4 * (7 * q + q1);
    let k = (q2 as f64 / 365.2422).floor() as i64;
    let e = (365.2422 * k as f64).floor() as i64;
    let mut j = q2 - e + n - 1;
    let mut x = g + k;

    if j > 366 && x.rem_euclid(4) == 0 {
        j -= 366;
        x += 1;
    }
    if j > 365 && x.rem_euclid(4) > 0 {
        j -= 365;
        x += 1;
    }
    (x, j)
}

/// Converts a date of the Islamic calendar to the Julian or Gregorian
/// calendar
///
/// The result is in the Julian calendar before 1582 October 15 and in the
/// Gregorian calendar from then on, at 0h of the civil day. The Islamic day
/// begins at sunset on the evening before.
///
/// Returns `None` if the date is not valid: `year` < 1, `month` not in
/// `1..=12`, or `day` not within the length of the month.
///
/// Meeus, chapter 9 (method of M. Francœur, modified by D. Savoie and J.
/// Meeus).
///
/// # Arguments
///
/// * `date`: Date in the Islamic calendar
pub fn date_from_islamic(date: &IslamicDate) -> Option<Date> {
    let (h, m, d) = (date.year as i64, date.month as i64, date.day as i64);
    if h < 1 || !(1..=12).contains(&m) || d < 1 || d > islamic_month_days(h, m) {
        return None;
    }

    let (x, j) = islamic_to_julian_year_and_day(h, m, d);
    // Julian Day number (noon) of day J of the Julian year X, then 0h
    let jd = (365.25 * (x - 1) as f64).floor() + 1721423.0 + j as f64 - 0.5;
    Date::from_julian_day(jd)
}

/// Meeus' conversion of a Julian-calendar date to the Islamic year `H` and
/// the day number `JJ` in that year, before the final adjustment.
fn julian_to_islamic_year_and_day(x: i64, m: i64, d: i64) -> (i64, i64) {
    let w = if x.rem_euclid(4) == 0 { 1 } else { 2 };
    let n = (275 * m) / 9 - w * ((m + 9) / 12) + d - 30;
    let a = x - 623;
    let b = a.div_euclid(4);
    let c = a.rem_euclid(4);
    let c1 = 365.2501 * c as f64;
    let mut c2 = c1.floor() as i64;
    if c1 - c2 as f64 > 0.5 {
        c2 += 1;
    }
    let d1 = 1461 * b + 170 + c2;
    let q = d1.div_euclid(10631);
    let r = d1.rem_euclid(10631);
    let j = r / 354;
    let k = r % 354;
    let o = (11 * j + 14) / 30;
    (30 * q + j + 1, k - o + n - 1)
}

/// Converts a Julian or Gregorian date to the Islamic calendar
///
/// The civil day is `date.decimal_day.floor()` in the calendar of `date`;
/// days outside the month are carried over arithmetically.
///
/// Returns `None` for dates before the start of the Islamic era (622 July
/// 16, Julian calendar), or if the Islamic year does not fit in an `i16`.
///
/// Meeus, chapter 9. Meeus adjusts the day number of the Islamic year
/// (his `JJ`) by at most one year, and maps a zero day number to day 355 of
/// the previous year. That is wrong for some dates: the day number can be
/// negative (early January), can exceed the length of two Islamic years
/// (late December), and after a common year of 354 days the zero day is
/// day 354. This function adjusts the year and the day number until the
/// day lies within the year, which agrees with the arithmetical Islamic
/// calendar on every day.
///
/// # Arguments
///
/// * `date`: Date in the Julian or Gregorian calendar
///
/// # Example
///
/// ```
/// use astro::calendar::{self, IslamicDate};
/// use astro::time::{CalType, Date, Month};
///
/// // Meeus, example 9.c: 1991 August 13 is 2 Safar 1412
/// let date = Date { year: 1991, month: Month::Aug, decimal_day: 13.0, cal_type: CalType::Gregorian };
/// assert_eq!(
///     calendar::islamic_from_date(&date),
///     Some(IslamicDate { year: 1412, month: 2, day: 2 })
/// );
/// ```
pub fn islamic_from_date(date: &Date) -> Option<IslamicDate> {
    let jd = time::julian_day(&Date {
        decimal_day: date.decimal_day.floor(),
        ..*date
    });
    if !jd.is_finite() || jd < ISLAMIC_EPOCH_JD {
        return None;
    }
    // Meeus first converts a Gregorian date to the Julian calendar.
    let julian = Date::from_julian_day_in(jd, CalType::Julian)?;
    let (mut h, mut jj) = julian_to_islamic_year_and_day(
        julian.year as i64,
        julian.month as i64,
        julian.decimal_day as i64,
    );

    // Bring the day number within the Islamic year (see the docs above).
    while jj > islamic_year_days(h) {
        jj -= islamic_year_days(h);
        h += 1;
    }
    while jj < 1 {
        h -= 1;
        jj += islamic_year_days(h);
    }

    let (month, day) = if jj == 355 {
        (12, 30)
    } else {
        let s = ((jj - 1) as f64 / 29.5).floor() as i64;
        (1 + s, (jj as f64 - 29.5 * s as f64).floor() as i64)
    };

    Some(IslamicDate {
        year: i16::try_from(h).ok()?,
        month: month as u8,
        day: day as u8,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    // Meeus, example 9.a: intermediate values
    #[test]
    fn pesach_1990_intermediate_values() {
        let (q, j, r, d) = pesach_terms(1990);
        assert!((q - 19.2599537042).abs() < 1e-9, "{q}");
        assert_eq!(j, 3);
        assert!((r - 0.2599537042).abs() < 1e-9);
        assert_eq!(d, 19 + 22);
    }

    // Meeus, example 9.b: intermediate values (N = 1, ..., J = 84, X = 2000)
    #[test]
    fn islamic_1421_intermediate_values() {
        assert_eq!(islamic_to_julian_year_and_day(1421, 1, 1), (2000, 84));
        let jdn = (365.25 * 1999.0_f64).floor() + 1721423.0 + 84.0;
        assert_eq!(jdn, 2451641.0);
    }

    // Meeus, example 9.c: intermediate values (H = 1411 and JJ = 386 before
    // the adjustment to H = 1412, JJ = 32)
    #[test]
    fn julian_1991_july_31_intermediate_values() {
        assert_eq!(julian_to_islamic_year_and_day(1991, 7, 31), (1411, 386));
        assert_eq!(islamic_year_days(1411), 354);
    }
}

//! Meeus, chapter 7: Julian Day, calendar dates, day of week and of year.

mod common;

use astro::time::{self, CalType, Date, DayOfMonth, Month, Weekday};

fn date(year: i16, month: Month, day: f64, cal_type: CalType) -> Date {
    Date {
        year,
        month,
        decimal_day: day,
        cal_type,
    }
}

fn greg(year: i16, month: Month, day: f64) -> Date {
    date(year, month, day, CalType::Gregorian)
}

fn jul(year: i16, month: Month, day: f64) -> Date {
    date(year, month, day, CalType::Julian)
}

fn days_in_month(year: i16, month: Month, cal_type: CalType) -> u8 {
    match month {
        Month::Feb if time::is_leap_year(year, &cal_type) => 29,
        Month::Feb => 28,
        Month::Apr | Month::June | Month::Sept | Month::Nov => 30,
        _ => 31,
    }
}

// Meeus, example 7.a
#[test]
fn example_7a_sputnik() {
    let d = DayOfMonth {
        day: 4,
        hr: 19,
        min: 26,
        sec: 24.0,
        time_zone: 0.0,
    };
    let dd = time::decimal_day(&d);
    common::assert_close(dd, 4.81, 1e-12);
    common::assert_close(
        time::julian_day(&greg(1957, Month::Oct, dd)),
        2436116.31,
        1e-8,
    );
}

// Meeus, example 7.b
#[test]
fn example_7b_year_333() {
    assert_eq!(time::julian_day(&jul(333, Month::Jan, 27.5)), 1842713.0);
}

// Meeus, chapter 7, list of test dates (p. 62), via the date → JD → date round trip
#[test]
fn list_of_test_dates_round_trip() {
    let cases = [
        (greg(2000, Month::Jan, 1.5), 2451545.0),
        (greg(1999, Month::Jan, 1.0), 2451179.5),
        (greg(1987, Month::Jan, 27.0), 2446822.5),
        (greg(1987, Month::June, 19.5), 2446966.0),
        (greg(1988, Month::Jan, 27.0), 2447187.5),
        (greg(1988, Month::June, 19.5), 2447332.0),
        (greg(1900, Month::Jan, 1.0), 2415020.5),
        (greg(1600, Month::Jan, 1.0), 2305447.5),
        (greg(1600, Month::Dec, 31.0), 2305812.5),
        (jul(837, Month::Apr, 10.3), 2026871.8),
        (jul(-123, Month::Dec, 31.0), 1676496.5),
        (jul(-122, Month::Jan, 1.0), 1676497.5),
        (jul(-1000, Month::July, 12.5), 1356001.0),
        (jul(-1000, Month::Feb, 29.0), 1355866.5),
        (jul(-1001, Month::Aug, 17.9), 1355671.4),
        (jul(-4712, Month::Jan, 1.5), 0.0),
    ];
    for (d, jd) in cases {
        common::assert_close(time::julian_day(&d), jd, 1e-8);
        let back = Date::from_julian_day(jd).unwrap();
        assert_eq!(
            (back.year, back.month, back.cal_type),
            (d.year, d.month, d.cal_type)
        );
        common::assert_close(back.decimal_day, d.decimal_day, 1e-6);
    }
}

// Meeus, example 7.c and the exercise that follows it
#[test]
fn example_7c_date_from_jd() {
    let cases = [
        (2436116.31, 1957, Month::Oct, 4.81, CalType::Gregorian),
        (1842713.0, 333, Month::Jan, 27.5, CalType::Julian),
        (1507900.13, -584, Month::May, 28.63, CalType::Julian),
    ];
    for (jd, y, m, d, c) in cases {
        let date = Date::from_julian_day(jd).unwrap();
        assert_eq!((date.year, date.month, date.cal_type), (y, m, c));
        common::assert_close(date.decimal_day, d, 1e-6);

        let (y2, m2, d2) = time::date_frm_julian_day(jd).unwrap();
        assert_eq!((y2, m2), (y, m as u8));
        common::assert_close(d2, d, 1e-6);
    }
}

// Meeus, example 7.d
#[test]
fn example_7d_halley_interval() {
    let from = greg(1910, Month::Apr, 20.0);
    let to = greg(1986, Month::Feb, 9.0);
    assert_eq!(time::julian_day(&from), 2418781.5);
    assert_eq!(time::julian_day(&to), 2446470.5);
    assert_eq!(time::days_between(&from, &to), 27689.0);
    assert_eq!(time::days_between(&to, &from), -27689.0);
}

// Meeus, exercise after example 7.d: 10000 days after 1991 July 11
#[test]
fn exercise_10000_days_later() {
    let jd = time::julian_day(&greg(1991, Month::July, 11.0)) + 10000.0;
    let d = Date::from_julian_day(jd).unwrap();
    assert_eq!(d, greg(2018, Month::Nov, 26.0));
}

// Meeus, example 7.e
#[test]
fn example_7e_weekday() {
    let d = greg(1954, Month::June, 30.0);
    assert_eq!(time::julian_day(&d), 2434923.5);
    assert_eq!(time::weekday_frm_date(&d), Weekday::Wednesday);
    // any instant of that day
    for frac in [0.0, 0.25, 0.5, 0.999] {
        let dd = greg(1954, Month::June, 30.0 + frac);
        assert_eq!(time::weekday_frm_date(&dd), Weekday::Wednesday);
        assert_eq!(
            time::weekday_from_julian_day(time::julian_day(&dd)),
            Weekday::Wednesday
        );
    }
}

// The week runs on across the Gregorian reform (Meeus, chapter 7).
#[test]
fn weekday_across_gregorian_reform() {
    assert_eq!(
        time::weekday_frm_date(&jul(1582, Month::Oct, 4.0)),
        Weekday::Thursday
    );
    assert_eq!(
        time::weekday_frm_date(&greg(1582, Month::Oct, 15.0)),
        Weekday::Friday
    );
    // The same instant named in both calendars has the same weekday.
    assert_eq!(
        time::weekday_frm_date(&jul(2000, Month::Jan, 1.0)),
        time::weekday_frm_date(&greg(2000, Month::Jan, 14.0))
    );
}

#[test]
fn weekday_cycles_through_negative_jd() {
    let order = [
        Weekday::Sunday,
        Weekday::Monday,
        Weekday::Tuesday,
        Weekday::Wednesday,
        Weekday::Thursday,
        Weekday::Friday,
        Weekday::Saturday,
    ];
    // JD 0 (-4712 January 1.5) was a Monday
    assert_eq!(time::weekday_from_julian_day(0.0), Weekday::Monday);
    assert_eq!(time::weekday_from_julian_day(-0.5), Weekday::Monday);
    assert_eq!(time::weekday_from_julian_day(-0.51), Weekday::Sunday);
    for k in -30_i32..30 {
        let jd = 0.0 + k as f64;
        let expected = order[(1 + k).rem_euclid(7) as usize];
        assert_eq!(time::weekday_from_julian_day(jd), expected, "JD {jd}");
    }
    // non-finite input does not panic
    let _ = time::weekday_from_julian_day(f64::NAN);
    let _ = time::weekday_from_julian_day(f64::INFINITY);
}

// Meeus, examples 7.f and 7.g
#[test]
fn examples_7f_7g_day_of_year() {
    assert_eq!(time::day_of_year(&greg(1978, Month::Nov, 14.0)), 318.0);
    assert_eq!(time::day_of_year(&greg(1988, Month::Apr, 22.0)), 113.0);
    // the fraction of the day is kept
    assert_eq!(time::day_of_year(&greg(1988, Month::Apr, 22.25)), 113.25);
    // limits
    assert_eq!(time::day_of_year(&greg(1978, Month::Jan, 1.0)), 1.0);
    assert_eq!(time::day_of_year(&greg(1978, Month::Dec, 31.0)), 365.0);
    assert_eq!(time::day_of_year(&greg(1988, Month::Dec, 31.0)), 366.0);
    assert_eq!(time::day_of_year(&greg(1900, Month::Dec, 31.0)), 365.0);
    assert_eq!(time::day_of_year(&jul(1900, Month::Dec, 31.0)), 366.0);
}

// Day of the year from the JD of January 0.0 and its inverse
// (Meeus, chapter 7, Pouplier's algorithm)
#[test]
fn day_of_year_round_trip() {
    for (year, cal_type) in [
        (1978, CalType::Gregorian),
        (1988, CalType::Gregorian),
        (1900, CalType::Gregorian),
        (2000, CalType::Gregorian),
        (1900, CalType::Julian),
        (-1, CalType::Julian),
        (0, CalType::Julian),
        (-100, CalType::Gregorian),
    ] {
        let jd0 = time::julian_day_of_jan0(year, cal_type);
        let days = if time::is_leap_year(year, &cal_type) {
            366
        } else {
            365
        };
        let mut count = 0;
        for month_number in 1..=12 {
            let month = Month::from_number(month_number).unwrap();
            for day in 1..=days_in_month(year, month, cal_type) {
                count += 1;
                let d = date(year, month, day as f64, cal_type);
                let n = time::day_of_year(&d);
                assert_eq!(n, count as f64);
                assert_eq!(time::julian_day(&d) - jd0, n);
                assert_eq!(time::date_from_day_of_year(year, n, cal_type), Some(d));
                let with_frac = time::date_from_day_of_year(year, n + 0.75, cal_type).unwrap();
                assert_eq!(with_frac, date(year, month, day as f64 + 0.75, cal_type));
            }
        }
        assert_eq!(count, days);
        assert_eq!(time::date_from_day_of_year(year, 0.999, cal_type), None);
        assert_eq!(
            time::date_from_day_of_year(year, days as f64 + 1.0, cal_type),
            None
        );
        assert_eq!(time::date_from_day_of_year(year, f64::NAN, cal_type), None);
    }
}

// Meeus, chapter 7: JD of January 0.0 of a Gregorian year
#[test]
fn julian_day_of_jan0_matches_book_formulas() {
    for year in 1583_i16..3000 {
        let y = (year - 1) as f64;
        let a = (y / 100.0).floor();
        let jd0 = (365.25 * y).floor() - a + (a / 4.0).floor() + 1721424.5;
        assert_eq!(time::julian_day_of_jan0(year, CalType::Gregorian), jd0);
        if (1901..=2099).contains(&year) {
            assert_eq!(jd0, 1721409.5 + (365.25 * y).floor());
        }
    }
    // 2000 January 0.0 = 1999 December 31.0
    assert_eq!(
        time::julian_day_of_jan0(2000, CalType::Gregorian),
        time::julian_day(&greg(1999, Month::Dec, 31.0))
    );
}

// Meeus, chapter 7: between 1900 March 1 and 2100 February 28, B = -13.
#[test]
fn gregorian_minus_julian_is_13_days_in_1900_2100() {
    for year in [1900_i16, 1950, 2000, 2099, 2100] {
        for m in 1..=12 {
            let month = Month::from_number(m).unwrap();
            if (year == 1900 && m < 3) || (year == 2100 && m > 2) {
                continue;
            }
            let g = time::julian_day(&greg(year, month, 10.0));
            let j = time::julian_day(&jul(year, month, 10.0));
            assert_eq!(g - j, -13.0, "{year}-{m}");
        }
    }
    let g = time::julian_day(&greg(1900, Month::Feb, 28.0));
    let j = time::julian_day(&jul(1900, Month::Feb, 28.0));
    assert_eq!(g - j, -12.0);
}

// Leap years (Meeus, chapter 7)
#[test]
fn leap_years() {
    for y in [900, 1236, 1700, 1800, 1900, 0, -4, -8, -1000] {
        assert!(time::is_leap_year(y, &CalType::Julian), "{y}");
    }
    for y in [750, 1429, -1, -3, -5, 1, 2023] {
        assert!(!time::is_leap_year(y, &CalType::Julian), "{y}");
    }
    for y in [1600, 2000, 2400, 1996, 0, -400, -4] {
        assert!(time::is_leap_year(y, &CalType::Gregorian), "{y}");
    }
    for y in [1700, 1800, 1900, 2100, -100, -1, 2023] {
        assert!(!time::is_leap_year(y, &CalType::Gregorian), "{y}");
    }
}

// The Gregorian reform: 1582 October 4 (Julian) is followed by October 15 (Gregorian).
#[test]
fn gregorian_reform_switch() {
    let oct4 = jul(1582, Month::Oct, 4.0);
    let oct15 = greg(1582, Month::Oct, 15.0);
    assert_eq!(time::julian_day(&oct4), 2299159.5);
    assert_eq!(time::julian_day(&oct15), time::GREGORIAN_REFORM_JD);
    assert_eq!(time::days_between(&oct4, &oct15), 1.0);

    assert_eq!(
        Date::from_julian_day(time::GREGORIAN_REFORM_JD),
        Some(oct15)
    );
    let before = Date::from_julian_day(time::GREGORIAN_REFORM_JD - 1e-6).unwrap();
    assert_eq!(
        (before.year, before.month, before.cal_type),
        (1582, Month::Oct, CalType::Julian)
    );
    common::assert_close(before.decimal_day, 5.0, 1e-5);
    assert!(before.decimal_day < 5.0);

    // Explicit calendars ignore the switch.
    let jd = time::julian_day(&greg(1000, Month::Mar, 1.0));
    assert_eq!(
        Date::from_julian_day_in(jd, CalType::Gregorian),
        Some(greg(1000, Month::Mar, 1.0))
    );
    assert_eq!(
        Date::from_julian_day_in(jd, CalType::Julian),
        Some(jul(1000, Month::Feb, 24.0))
    );
    let jd = time::julian_day(&greg(2024, Month::Jan, 1.0));
    assert_eq!(
        Date::from_julian_day_in(jd, CalType::Julian),
        Some(jul(2023, Month::Dec, 19.0))
    );
}

// JD 0 and negative Julian Days
#[test]
fn jd_zero_and_negative() {
    assert_eq!(
        Date::from_julian_day(0.0),
        Some(jul(-4712, Month::Jan, 1.5))
    );
    assert_eq!(
        Date::from_julian_day(-0.5),
        Some(jul(-4712, Month::Jan, 1.0))
    );
    assert_eq!(
        Date::from_julian_day(-1.0),
        Some(jul(-4713, Month::Dec, 31.5))
    );
    // -4713 is not a leap year, -4716 is
    assert_eq!(
        Date::from_julian_day(time::julian_day(&jul(-4713, Month::Mar, 1.0))),
        Some(jul(-4713, Month::Mar, 1.0))
    );
    assert_eq!(
        Date::from_julian_day(time::julian_day(&jul(-4716, Month::Feb, 29.0))),
        Some(jul(-4716, Month::Feb, 29.0))
    );
    // JD 0 in the proleptic Gregorian calendar: -4713 November 24.5
    assert_eq!(
        Date::from_julian_day_in(0.0, CalType::Gregorian),
        Some(greg(-4713, Month::Nov, 24.5))
    );
    assert_eq!(time::julian_day(&greg(-4713, Month::Nov, 24.5)), 0.0);

    let (y, m, d) = time::date_frm_julian_day(-1_000_000.0).unwrap();
    let back = time::julian_day(&jul(y, Month::from_number(m).unwrap(), d));
    assert_eq!(back, -1_000_000.0);
}

// Date → JD → date for every day over a wide span, and day-by-day continuity
#[test]
fn round_trip_every_day() {
    for cal_type in [None, Some(CalType::Julian), Some(CalType::Gregorian)] {
        let mut prev: Option<Date> = None;
        let mut jd = -2_000_000.5;
        while jd < 4_000_000.0 {
            let d = match cal_type {
                None => Date::from_julian_day(jd),
                Some(c) => Date::from_julian_day_in(jd, c),
            }
            .unwrap();
            assert_eq!(time::julian_day(&d), jd, "{d:?}");
            assert_eq!(d.decimal_day.fract(), 0.0);
            let dim = days_in_month(d.year, d.month, d.cal_type) as f64;
            assert!(d.decimal_day >= 1.0 && d.decimal_day <= dim, "{d:?}");
            if let Some(p) = prev {
                let consecutive = if d.decimal_day == 1.0 {
                    p.decimal_day == days_in_month(p.year, p.month, p.cal_type) as f64
                } else {
                    p.decimal_day + 1.0 == d.decimal_day && p.month == d.month
                };
                let reform = p.cal_type != d.cal_type;
                assert!(consecutive || reform, "{p:?} -> {d:?}");
            }
            prev = Some(d);
            jd += 1.0;
        }
    }
}

#[test]
fn round_trip_fractional_jd() {
    let mut jd = -1.0e6 + 0.123;
    while jd < 5.0e6 {
        let d = Date::from_julian_day(jd).unwrap();
        common::assert_close(time::julian_day(&d), jd, 1e-8);
        let g = Date::from_julian_day_in(jd, CalType::Gregorian).unwrap();
        common::assert_close(time::julian_day(&g), jd, 1e-8);
        jd += 1234.567;
    }
}

#[test]
fn from_julian_day_bad_input() {
    assert_eq!(Date::from_julian_day(f64::NAN), None);
    assert_eq!(Date::from_julian_day(f64::INFINITY), None);
    assert_eq!(Date::from_julian_day(f64::NEG_INFINITY), None);
    // years outside i16
    assert_eq!(Date::from_julian_day(1.0e8), None);
    assert_eq!(Date::from_julian_day(-1.0e8), None);
    assert_eq!(Date::from_julian_day(1.0e300), None);
    assert!(time::date_frm_julian_day(f64::NAN).is_err());
    assert!(time::date_frm_julian_day(-1.0e8).is_err());
    // the limits of i16 still work
    let max = time::julian_day(&greg(i16::MAX, Month::Dec, 31.0));
    assert_eq!(
        Date::from_julian_day(max),
        Some(greg(i16::MAX, Month::Dec, 31.0))
    );
    assert_eq!(Date::from_julian_day(max + 1.0), None);
    let min = time::julian_day(&jul(i16::MIN, Month::Jan, 1.0));
    assert_eq!(
        Date::from_julian_day(min),
        Some(jul(i16::MIN, Month::Jan, 1.0))
    );
    assert_eq!(Date::from_julian_day(min - 1.0), None);
}

// Modified Julian Day (Meeus, chapter 7)
#[test]
fn modified_julian_day() {
    let jd = time::julian_day(&greg(1858, Month::Nov, 17.0));
    assert_eq!(time::modified_julian_day(jd), 0.0);
    assert_eq!(time::modified_julian_day(2451545.0), 51544.5);
    assert_eq!(time::julian_day_from_modified(51544.5), 2451545.0);
    for mjd in [-1.0e6, -0.5, 0.0, 12345.678, 60000.0] {
        common::assert_close(
            time::modified_julian_day(time::julian_day_from_modified(mjd)),
            mjd,
            1e-9,
        );
    }
}

// decimal_day: issues #8 and #17 (minutes and seconds must be fractions of a
// day, and the time zone is subtracted to give UT).
#[test]
fn decimal_day_hours_minutes_seconds_and_zone() {
    let d = |day, hr, min, sec, time_zone| {
        time::decimal_day(&DayOfMonth {
            day,
            hr,
            min,
            sec,
            time_zone,
        })
    };
    assert_eq!(d(1, 0, 0, 0.0, 0.0), 1.0);
    assert_eq!(d(1, 12, 0, 0.0, 0.0), 1.5);
    common::assert_close(d(1, 0, 1, 0.0, 0.0), 1.0 + 1.0 / 1440.0, 1e-15);
    common::assert_close(d(1, 0, 0, 1.0, 0.0), 1.0 + 1.0 / 86400.0, 1e-15);
    common::assert_close(d(15, 23, 59, 59.999, 0.0), 16.0, 1e-7);
    // 12h local time at UTC+5:30 is 6h30m UT
    common::assert_close(d(1, 12, 0, 0.0, 5.5), 1.0 + 6.5 / 24.0, 1e-12);
    // 22h local time at UTC-8 is 6h UT on the next day
    common::assert_close(d(1, 22, 0, 0.0, -8.0), 2.25, 1e-12);
    // 1h local time at UTC+3 is 22h UT on the previous day: the decimal day
    // leaves the month, and julian_day still gives the right instant
    let dd = d(1, 1, 0, 0.0, 3.0);
    common::assert_close(dd, 22.0 / 24.0, 1e-12);
    common::assert_close(
        time::julian_day(&greg(2000, Month::Jan, dd)),
        time::julian_day(&greg(1999, Month::Dec, 31.0 + 22.0 / 24.0)),
        1e-9,
    );
}

#[test]
fn decimal_year() {
    common::assert_close(
        time::decimal_year(&greg(2001, Month::Jan, 1.0)),
        2001.0,
        1e-12,
    );
    // 2001 July 2.5 is 182.5 days after January 1.0 of a 365-day year
    common::assert_close(
        time::decimal_year(&greg(2001, Month::July, 2.5)),
        2001.5,
        1e-12,
    );
    // 2000 is a leap year: July 2.0 is 183 days after January 1.0
    common::assert_close(
        time::decimal_year(&greg(2000, Month::July, 2.0)),
        2000.5,
        1e-12,
    );
    // continuous and increasing across a year boundary
    let end = time::decimal_year(&greg(1999, Month::Dec, 31.9999));
    let start = time::decimal_year(&greg(2000, Month::Jan, 1.0));
    assert!(end < start && start - end < 1e-6);
    // Julian calendar, negative year
    common::assert_close(
        time::decimal_year(&jul(-100, Month::Jan, 1.0)),
        -100.0,
        1e-12,
    );
}

// Dates and Julian Ephemeris Days (Meeus, chapter 7: 1977 April 26.4)
#[test]
fn julian_ephemeris_day_round_trip() {
    let jd = time::julian_day(&greg(1977, Month::Apr, 26.4));
    common::assert_close(jd, 2443259.9, 1e-8);
    let jde = time::julian_ephemeris_day(jd, 48.0);
    common::assert_close(jde - jd, 48.0 / 86400.0, 1e-9);
    common::assert_close(time::julian_day_from_ephemeris_day(jde, 48.0), jd, 1e-9);
}

#[test]
fn julian_centuries_and_millennia() {
    assert_eq!(time::julian_cent(2451545.0), 0.0);
    assert_eq!(time::julian_cent(2451545.0 + 36525.0), 1.0);
    assert_eq!(time::julian_mill(2451545.0 - 365250.0), -1.0);
    // Meeus, example 12.a
    common::assert_close(time::julian_cent(2446895.5), -0.127296372348, 1e-12);
}

#[test]
fn month_numbers() {
    for n in 1..=12 {
        assert_eq!(Month::from_number(n).unwrap().number(), n);
    }
    assert_eq!(Month::from_number(0), None);
    assert_eq!(Month::from_number(13), None);
    assert!(Month::Jan < Month::Dec);
}

#[test]
fn julian_day_extreme_years_are_finite() {
    for y in [i16::MIN, i16::MAX] {
        for c in [CalType::Julian, CalType::Gregorian] {
            assert!(time::julian_day(&date(y, Month::Jan, 1.0, c)).is_finite());
        }
    }
}

// Review R1: for Julian Days close to zero the fraction of the day is so
// fine that adding it to the day of the month could round up to the next
// day, giving -4713 December 32.0. The instant is now carried into the
// next day instead, so the day stays within the month.
#[test]
fn from_julian_day_fraction_never_rounds_past_month_end() {
    let just_before = f64::from_bits((-0.5_f64).to_bits() + 1); // -0.5000000000000001
    assert!(just_before < -0.5);
    for cal_type in [CalType::Julian, CalType::Gregorian] {
        let d = Date::from_julian_day_in(just_before, cal_type).unwrap();
        let dim = days_in_month(d.year, d.month, d.cal_type) as f64;
        assert!(d.decimal_day < dim + 1.0, "{d:?}");
        common::assert_close(time::julian_day(&d), just_before, 1e-15);
    }
    let d = Date::from_julian_day(just_before).unwrap();
    assert_eq!(d, jul(-4712, Month::Jan, 1.0));

    // Every instant within a few ulps of a day boundary near JD 0
    let mut x = -40.5_f64;
    while x < 40.0 {
        for k in 1..4u64 {
            let jd = if x < 0.0 {
                f64::from_bits(x.to_bits() + k)
            } else {
                f64::from_bits(x.to_bits() - k)
            };
            let d = Date::from_julian_day(jd).unwrap();
            let dim = days_in_month(d.year, d.month, d.cal_type) as f64;
            assert!(
                d.decimal_day >= 1.0 && d.decimal_day < dim + 1.0,
                "{jd}: {d:?}"
            );
            // `decimal_day` (up to 32) resolves 7e-15 day, so instants a few
            // ulps of JD from midnight are rounded to midnight.
            common::assert_close(time::julian_day(&d), jd, 3e-14);
        }
        x += 1.0;
    }
}

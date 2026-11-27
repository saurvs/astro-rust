//! Meeus, chapter 8: date of Easter.

use astro::calendar;
use astro::time::{self, CalType, Date, Month, Weekday};

fn md(d: &Date) -> (Month, f64) {
    (d.month, d.decimal_day)
}

// Meeus, chapter 8: test years
#[test]
fn gregorian_easter_test_years() {
    for (year, month, day) in [
        (1991, Month::Mar, 31.0),
        (1992, Month::Apr, 19.0),
        (1993, Month::Apr, 11.0),
        (1954, Month::Apr, 18.0),
        (2000, Month::Apr, 23.0),
        (1818, Month::Mar, 22.0),
        // 1967: Easter on March 26 although the Full Moon was that day
        (1967, Month::Mar, 26.0),
    ] {
        let e = calendar::easter_gregorian(year).unwrap();
        assert_eq!(e.year, year);
        assert_eq!(e.cal_type, CalType::Gregorian);
        assert_eq!(md(&e), (month, day), "{year}");
    }
}

// Meeus, chapter 8: the extreme dates
#[test]
fn gregorian_easter_extremes() {
    for year in [1818, 2285] {
        assert_eq!(
            md(&calendar::easter_gregorian(year).unwrap()),
            (Month::Mar, 22.0)
        );
    }
    for year in [1886, 1943, 2038] {
        assert_eq!(
            md(&calendar::easter_gregorian(year).unwrap()),
            (Month::Apr, 25.0)
        );
    }
}

#[test]
fn gregorian_easter_before_1583_is_none() {
    assert!(calendar::easter_gregorian(1582).is_none());
    assert!(calendar::easter_gregorian(0).is_none());
    assert!(calendar::easter_gregorian(i16::MIN).is_none());
    assert!(calendar::easter_gregorian(1583).is_some());
    assert!(calendar::easter_gregorian(i16::MAX).is_some());
}

/// Independent algorithm (Lilius and Clavius, as given by Knuth, *The Art
/// of Computer Programming*, vol. 1, § 1.3.2, exercise 14).
fn knuth_easter(y: i32) -> (Month, f64) {
    let g = y % 19 + 1;
    let c = y / 100 + 1;
    let x = 3 * c / 4 - 12;
    let z = (8 * c + 5) / 25 - 5;
    let d = 5 * y / 4 - x - 10;
    let mut e = (11 * g + 20 + z - x).rem_euclid(30);
    if (e == 25 && g > 11) || e == 24 {
        e += 1;
    }
    let mut n = 44 - e;
    if n < 21 {
        n += 30;
    }
    n = n + 7 - (d + n) % 7;
    if n > 31 {
        (Month::Apr, (n - 31) as f64)
    } else {
        (Month::Mar, n as f64)
    }
}

// Every Gregorian Easter is a Sunday between March 22 and April 25, and
// agrees with Knuth's independent algorithm.
#[test]
fn gregorian_easter_all_years() {
    let mut april19 = 0;
    for year in 1583..=i16::MAX {
        let e = calendar::easter_gregorian(year).unwrap();
        assert_eq!(md(&e), knuth_easter(year as i32), "{year}");
        assert_eq!(time::weekday_frm_date(&e), Weekday::Sunday, "{year}");
        let n = time::day_of_year(&e)
            - time::day_of_year(&Date {
                decimal_day: 1.0,
                month: Month::Mar,
                ..e
            });
        assert!((21.0..=55.0).contains(&n), "{year}: {e:?}");
        if md(&e) == (Month::Apr, 19.0) {
            april19 += 1;
        }
    }
    // April 19 is the most frequent date in the long run (about 3.9 %).
    let share = april19 as f64 / (i16::MAX as f64 - 1582.0);
    assert!(share > 0.035 && share < 0.042, "{share}");
}

// Meeus, chapter 8: Julian Easter on April 12 in 179, 711 and 1243
#[test]
fn julian_easter_periodicity_examples() {
    for year in [179, 711, 1243] {
        let e = calendar::easter_julian(year);
        assert_eq!(e.year, year);
        assert_eq!(e.cal_type, CalType::Julian);
        assert_eq!(md(&e), (Month::Apr, 12.0), "{year}");
    }
}

#[test]
fn julian_easter_all_years() {
    for year in i16::MIN..=i16::MAX {
        let e = calendar::easter_julian(year);
        assert_eq!(time::weekday_frm_date(&e), Weekday::Sunday, "{year}");
        let ok = match e.month {
            Month::Mar => e.decimal_day >= 22.0 && e.decimal_day <= 31.0,
            Month::Apr => e.decimal_day >= 1.0 && e.decimal_day <= 25.0,
            _ => false,
        };
        assert!(ok, "{year}: {e:?}");
        if let Some(later) = year.checked_add(532) {
            assert_eq!(md(&calendar::easter_julian(later)), md(&e), "{year}");
        }
    }
}

// Eastern (Orthodox) Easter is the Julian Easter, shown in the Gregorian calendar.
#[test]
fn orthodox_easter_in_gregorian_calendar() {
    for (year, month, day) in [
        (2008, Month::Apr, 27.0),
        (2010, Month::Apr, 4.0),
        (2024, Month::May, 5.0),
    ] {
        let e = calendar::easter_julian(year);
        let g = Date::from_julian_day_in(time::julian_day(&e), CalType::Gregorian).unwrap();
        assert_eq!((g.year, g.month, g.decimal_day), (year, month, day));
    }
    // 2010: both churches celebrated on the same day
    let w = calendar::easter_gregorian(2010).unwrap();
    let e = calendar::easter_julian(2010);
    assert_eq!(time::julian_day(&w), time::julian_day(&e));
}

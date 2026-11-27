//! Regression tests for bugs fixed in the time module (group A2).

#![allow(non_snake_case)]

mod common;

use astro::time::{self, CalType, Date, Month, Weekday};

// ΔT: the NASA (Espenak & Meeus) polynomials must join up at the segment
// boundaries; wrong signs in the nested form broke several segments.
#[test]
fn delta_t_is_continuous_at_segment_boundaries() {
    for b in [
        -500.0, 500.0, 1600.0, 1700.0, 1800.0, 1860.0, 1900.0, 1920.0, 1941.0, 1961.0, 1986.0,
        2005.0, 2050.0, 2150.0,
    ] {
        let before = time::delta_t_from_decimal_year(b - 1e-9);
        let after = time::delta_t_from_decimal_year(b);
        common::assert_close(after, before, 0.3);
    }
}

#[test]
fn delta_t_segments_with_sign_errors() {
    // Values of the published polynomials (ascending-power form) inside the
    // segments whose nested form had wrong signs.
    let dt = time::delta_t_from_decimal_year;
    common::assert_close(dt(-100.0), 11637.55, 0.05);
    common::assert_close(dt(1750.0), 13.37, 0.005);
    common::assert_close(dt(1799.5), 13.968, 0.001);
    common::assert_close(dt(1859.5), 7.410, 0.001);
    common::assert_close(dt(1899.5), -3.387, 0.001);
    common::assert_close(dt(1919.5), 21.033, 0.001);
}

#[test]
fn weekday_respects_julian_calendar() {
    // 1582 October 4 (Julian) was a Thursday, followed by Friday 15 October (Gregorian).
    let d = Date {
        year: 1582,
        month: Month::Oct,
        decimal_day: 4.0,
        cal_type: CalType::Julian,
    };
    assert_eq!(time::weekday_frm_date(&d), Weekday::Thursday);
}

#[test]
fn weekday_does_not_panic_for_negative_jd() {
    // JD 0 = -4712 January 1.5 (Julian) was a Monday.
    let d = Date {
        year: -4712,
        month: Month::Jan,
        decimal_day: 1.5,
        cal_type: CalType::Julian,
    };
    assert_eq!(time::weekday_frm_date(&d), Weekday::Monday);
    let d = Date {
        year: -4713,
        month: Month::Dec,
        decimal_day: 25.0,
        cal_type: CalType::Julian,
    };
    assert_eq!(time::weekday_frm_date(&d), Weekday::Monday);
}

#[test]
fn decimal_year_starts_at_zero_on_january_first() {
    let jan1 = Date {
        year: 2001,
        month: Month::Jan,
        decimal_day: 1.0,
        cal_type: CalType::Gregorian,
    };
    common::assert_close(time::decimal_year(&jan1), 2001.0, 1e-12);
    let dec31 = Date {
        year: 2001,
        month: Month::Dec,
        decimal_day: 31.999,
        cal_type: CalType::Gregorian,
    };
    assert!(time::decimal_year(&dec31) < 2002.0);
}

#[test]
fn julian_day_does_not_overflow_for_extreme_years() {
    let d = Date {
        year: i16::MIN,
        month: Month::Jan,
        decimal_day: 1.0,
        cal_type: CalType::Julian,
    };
    assert!(time::julian_day(&d).is_finite());
}

#[test]
fn date_frm_julian_day_handles_negative_and_bad_input() {
    // -4713 December 31.5 (Julian) is JD -1.
    let (y, m, d) = time::date_frm_julian_day(-1.0).unwrap();
    assert_eq!((y, m), (-4713, 12));
    common::assert_close(d, 31.5, 1e-9);
    assert!(time::date_frm_julian_day(f64::NAN).is_err());
    assert!(time::date_frm_julian_day(1.0e8).is_err());
}

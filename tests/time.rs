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

#![allow(non_snake_case)]

#[macro_use]
extern crate astro;

use astro::*;

#[test]
fn sidreal_time() {

    let (h1, m1, s1) = angle::hms_frm_deg(time::mn_sidr(2446895.5).to_degrees());
    assert_eq!((h1, m1, util::round_upto_digits(s1, 4)), (13, 10, 46.3668));

    let (h2, m2, s2) = angle::hms_frm_deg(apprnt_sidr!(2446895.5).to_degrees());
    assert_eq!((h2, m2, util::round_upto_digits(s2, 4)), (13, 10, 46.1351));

}

#[test]
fn julian_day() {

    // Test taken from Meeus 2nd ed. on pages 61-62

    let mut date = time::Date {
        year        : 0,
        month       : time::Month::Jan,
        decimal_day : 0.0,
        cal_type    : time::CalType::Gregorian
    };

    struct TestDate(i16, time::Month, f64, f64, time::CalType);
    let gregorian_dates = [
        TestDate( 1957, time::Month::Oct, 4.81, 2436116.31, time::CalType::Gregorian),
        TestDate( 2000, time::Month::Jan,  1.5,  2451545.0,  time::CalType::Gregorian),
        TestDate( 1999, time::Month::Jan,  1.0,  2451179.5,  time::CalType::Gregorian),
        TestDate( 1987, time::Month::Jan,  27.0, 2446822.5,  time::CalType::Gregorian),
        TestDate( 1987, time::Month::June,  19.5, 2446966.0,  time::CalType::Gregorian),
        TestDate( 1988, time::Month::Jan,  27.0, 2447187.5,  time::CalType::Gregorian),
        TestDate( 1988, time::Month::June,  19.5, 2447332.0,  time::CalType::Gregorian),
        TestDate( 1900, time::Month::Jan,  1.0,  2415020.5,  time::CalType::Gregorian),
        TestDate( 1600, time::Month::Jan,  1.0,  2305447.5,  time::CalType::Gregorian),
        TestDate( 1600, time::Month::Dec, 31.0, 2305812.5,  time::CalType::Gregorian),
        TestDate( 837,  time::Month::Apr,  10.3, 2026871.8,  time::CalType::Julian),
        TestDate(-123,  time::Month::Dec, 31.0, 1676496.5,  time::CalType::Julian),
        TestDate(-122,  time::Month::Jan,  1.0,  1676497.5,  time::CalType::Julian),
        TestDate(-1000, time::Month::July,  12.5, 1356001.0,  time::CalType::Julian),
        TestDate(-1000, time::Month::Feb,  29.0, 1355866.5,  time::CalType::Julian),
        TestDate(-1001, time::Month::Aug,  17.9, 1355671.4,  time::CalType::Julian),
        TestDate(-4712, time::Month::Jan,  1.5,  0.0,        time::CalType::Julian)
    ];

    for date_fields in gregorian_dates.iter() {
        date.year = date_fields.0;
        date.month = date_fields.1;
        date.decimal_day = date_fields.2;
        date.cal_type = match date_fields.4 {
            time::CalType::Gregorian => time::CalType::Gregorian,
            time::CalType::Julian    => time::CalType::Julian,
        };
        assert_eq!(time::julian_day(&date), date_fields.3);
    }

}

#[test]
fn weekday_frm_date()
{
    let date = time::Date {
        year: 1954,
        month: time::Month::June,
        decimal_day: 30.0,
        cal_type: time::CalType::Gregorian
    };
    match time::weekday_frm_date(&date) {
        time::Weekday::Wednesday => {},
        _ => panic!("time::weekday_frm_date failed")
    }
}

#[test]
fn date_frm_julian_day() {

    // Test taken from Meeus 2nd ed. on page 64

    struct TestData(i16, u8, f64, f64);
    let test_data = [
        TestData( 1957, 10, 4.81,  2436116.31),
        TestData( 333,  1,  27.5,  1842713.0),
        TestData(-584,  5,  28.63, 1507900.13),
    ];

    for data in test_data.iter() {

        let (year, month, day) = time::date_frm_julian_day(data.3).unwrap();

        assert_eq!(year, data.0);
        assert_eq!(month, data.1);
        assert_eq!(util::round_upto_digits(day, 2), data.2);

    }

}

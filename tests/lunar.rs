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

extern crate astro;
use astro::*;

#[test]
fn geocent_ecl_pos() {

    let (moon_eq_point, rad_vec) = lunar::geocent_ecl_pos(2448724.5);

    assert_eq!(
        util::round_upto_digits(moon_eq_point.long.to_degrees(), 6),
        133.162655
    );
    assert_eq!(
        util::round_upto_digits(moon_eq_point.lat.to_degrees(), 6),
        -3.229126
    );
    assert_eq!(
        util::round_upto_digits(rad_vec, 1), 368409.7
    );

}

#[test]
#[allow(unused_variables)]
fn time_of_passage_through_nodes() {

    let date = time::Date {
        year        : 1987,
        month       : 5,
        decimal_day : 15.0,
        cal_type    : time::CalType::Gregorian
    };

    let (ascend_JD, desend_JD) = lunar::time_of_passage_through_nodes(&date);

    assert_eq!(
        util::round_upto_digits(ascend_JD, 5), 2446938.76803
    );

}

#[test]
fn liberations() {

    let day_of_month = time::DayOfMonth {
        day       : 12,
        hr        : 0,
        min       : 0,
        sec       : 0.0,
        time_zone : 0.0
    };
    let date = time::Date {
        year        : 1992,
        month       : 4,
        decimal_day : time::decimal_day(&day_of_month),
        cal_type    : time::CalType::Gregorian
    };

    let (opt_long, opt_lat) = lunar::optical_libr (
        time::julian_day(&date),
        133.162655_f64.to_radians(),
       -3.229126_f64.to_radians()
    );
    assert_eq!(
        util::round_upto_digits(opt_long.to_degrees(), 3),
        angle::limit_to_360(-1.206)
    );
    assert_eq!(
        util::round_upto_digits(opt_lat.to_degrees(), 3), 4.194
    );

    let (phy_long, phys_lat) = lunar::physical_libr (
        time::julian_day(&date),
        133.162655_f64.to_radians(),
       -3.229126_f64.to_radians(),
        opt_lat
    );
    assert_eq!(util::round_upto_digits(phy_long.to_degrees(), 3), -0.025);
    assert_eq!(util::round_upto_digits(phys_lat.to_degrees(), 3), 0.006);

}

#[test]
fn phases() {
    let date_last_quarter = time::Date {
        year        : 2044,
        month       : 1,
        decimal_day : 0.0,
        cal_type    : time::CalType::Gregorian
    };
    let JD_last_quarter = lunar::time_of_phase (
        &date_last_quarter, &lunar::Phase::Last
    );
    assert_eq!(util::round_upto_digits(JD_last_quarter, 5), 2467636.49186);

    let date_new_moon = time::Date {
        year        : 1977,
        month       : 2,
        decimal_day : 0.0,
        cal_type    : time::CalType::Gregorian};
    let JD_new_moon = lunar::time_of_phase(&date_new_moon, &lunar::Phase::New);
    assert_eq!(util::round_upto_digits(JD_new_moon, 5), 2443192.65118);
}

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
fn true_anom_and_rad_vec() {

    let t_date = time::Date {
        year        : 1998,
        month       : time::Month::Aug,
        decimal_day : 5.0,
        cal_type    : time::CalType::Gregorian
    };
    let t = time::julian_day(&t_date);

    let T_date = time::Date {
        year        : 1998,
        month       : time::Month::Apr,
        decimal_day : 14.4358,
        cal_type    : time::CalType::Gregorian
    };
    let T = time::julian_day(&T_date);

    let (tru_anom, rad_vec) = orbit::parabolic::true_anom_and_rad_vec(
        t, T, 1.487469
    );
    assert_eq!(util::round_upto_digits(tru_anom.to_degrees(), 5), 66.78862);
    assert_eq!(util::round_upto_digits(rad_vec, 6), 2.133911);

}

#[test]
fn passage_through_nodes() {

    let T = time::julian_day(&time::Date{
        year        : 1989,
        month       : time::Month::Aug,
        decimal_day : 20.291,
        cal_type    : time::CalType::Gregorian
    });
    let w = 154.9103_f64.to_radians();
    let q = 1.324502;

    let (ascen, r_a) = orbit::parabolic::passage_through_node(
        w, q, T, &orbit::Node::Ascend
    );
    assert_eq!(util::round_upto_digits((T - ascen), 2), 4354.65);
    assert_eq!(util::round_upto_digits(r_a, 2), 28.07);

    let (descend, r_b) = orbit::parabolic::passage_through_node(
        w, q, T, &orbit::Node::Descend
    );
    assert_eq!(util::round_upto_digits((T - descend), 4), -28.3454);
    assert_eq!(util::round_upto_digits(r_b, 4), 1.3901);

}

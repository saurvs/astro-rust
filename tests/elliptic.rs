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
fn ecc_anom() {

    let ecc_anom = orbit::elliptic::ecc_anom(5_f64.to_radians(), 0.1, 0.0000001);

    assert_eq!(util::round_upto_digits(ecc_anom.to_degrees(), 6), 5.554589);

}

#[test]
fn vel() {

    let a = 17.9400782;
    let e = 0.96727426;
    let (vel_p, vel_a) = (orbit::elliptic::perih_vel(a, e), orbit::elliptic::aph_vel(a, e));
    assert_eq!(util::round_upto_digits(vel_p, 2), 54.52);
    assert_eq!(util::round_upto_digits(vel_a, 2), 0.91);

    let V = orbit::elliptic::vel(1.0, a);
    assert_eq!(util::round_upto_digits(V, 2), 41.53);

}

#[test]
fn passage_through_nodes() {

    let T = time::julian_day (
        &time::Date{
            year        : 1986,
            month       : 2,
            decimal_day : 9.45891,
            cal_type    : time::CalType::Gregorian
        }
    );
    let w = 111.84644_f64.to_radians();
    let e = 0.96727426;
    let n = 0.01297082_f64.to_radians();
    let a = 17.9400782;

    let (ascen, r_a) = orbit::elliptic::passage_through_node (
        w, n, a, e, T, &orbit::Node::Ascend
    );
    assert_eq!(util::round_upto_digits((T - ascen), 4), 92.2998);
    assert_eq!(util::round_upto_digits(r_a, 4), 1.8045);

    let (descend, r_b) = orbit::elliptic::passage_through_node (
        w, n, a, e, T, &orbit::Node::Descend
    );
    assert_eq!(util::round_upto_digits((T - descend), 4), -28.9105);
    assert_eq!(util::round_upto_digits(r_b, 4), 0.8493);

}

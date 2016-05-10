/*
Copyright (c) 2015 Saurav Sachidanand

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
fn nutation() {
    let (nut_in_long, nut_in_oblq) = nutation::nutation(2446895.5);

    let (d1, m1, s1) = angle::dms_frm_deg(nut_in_long.to_degrees());
    assert_eq!((d1, m1, util::round_upto_digits(s1, 3)), (0, 0, -3.788));

    let (d2, m2, s2) = angle::dms_frm_deg(nut_in_oblq.to_degrees());
    assert_eq!((d2, m2, util::round_upto_digits(s2, 3)), (0, 0, 9.443));
}

#[test]
fn nutation_in_eq_coords() {
    let eq_point = coords::EqPoint{asc: 41.5555635_f64.to_radians(), dec: 49.3503415_f64.to_radians()};
    let (a, b) = nutation::nutation_in_eq_coords(
        &eq_point,
        angle::deg_frm_dms(0, 0, 14.861).to_radians(),
        angle::deg_frm_dms(0, 0, 2.705).to_radians(),
        23.436_f64.to_radians()
    );

    assert_eq!(util::round_upto_digits(a.to_degrees(), 7), 0.0044011);
    assert_eq!(util::round_upto_digits(b.to_degrees(), 7), 0.001727);
}

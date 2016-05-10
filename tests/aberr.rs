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
fn stell_aberr_in_eq_coords() {
    let d = time::Date{year: 2028, month: 11, decimal_day: 13.19, cal_type: time::CalType::Gregorian};
    let stell_eq_point = coords::EqPoint{
        asc: 41.0540613_f64.to_radians(),
        dec: 49.2277489_f64.to_radians()
    };
    let (a, b) = aberr::stell_aberr_in_eq_coords(&stell_eq_point, time::julian_day(&d));

    assert_eq!(util::round_upto_digits(a.to_degrees(), 7), 0.0083223);
    assert_eq!(util::round_upto_digits(b.to_degrees(), 7), 0.0018749);
}

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
fn annual_precess() {

    let d = time::Date {
        year        : 1978,
        month       : 1,
        decimal_day : 0.0,
        cal_type    : time::CalType::Gregorian};

    let (new_asc, new_dec) = precess::annual_precess(
        angle::deg_frm_hms(10, 8, 22.3).to_radians(),
        angle::deg_frm_dms(11, 58, 2.0).to_radians(),
        time::julian_day(&d)
    );

    let (a, b, c) = angle::hms_frm_deg(new_asc.to_degrees());
    assert_eq!((a, b), (0, 0));
    assert_eq!(util::round_upto_digits(c, 2), util::round_upto_digits(3.208, 2));

    let (d, e, f) = angle::dms_frm_deg(new_dec.to_degrees());
    assert_eq!((d, e), (0, 0));
    assert_eq!(util::round_upto_digits(f, 2), -17.71);

}

#[test]
fn precess_eq_coords() {

    let (new_asc, new_dec) = precess::precess_eq_coords(
        41.054063_f64.to_radians(),
        49.22775_f64.to_radians(),
        2451545.0,
        2462088.69
    );
    assert_eq!(
        (util::round_upto_digits(new_asc.to_degrees(), 6),
        util::round_upto_digits(new_dec.to_degrees(), 6)),
        (41.547214, 49.348483)
    );

}

#[test]
fn precess_orb_elements() {

    let (new_inc, new_arg_perih, new_long_ascend_node) =
        precess::precess_orb_elements (
            47.122_f64.to_radians(),
            151.4486_f64.to_radians(),
            45.7481_f64.to_radians(),
            2358042.5305,
            2433282.4235
        );

    assert_eq!(util::round_upto_digits(new_inc.to_degrees(), 4), 47.138);
    assert_eq!(util::round_upto_digits(new_arg_perih.to_degrees(), 4), 151.4782);
    assert_eq!(util::round_upto_digits(new_long_ascend_node.to_degrees(), 4), 48.6037);

}

#[test]
fn precess_ecl_coords() {

    let (new_asc, new_dec) = precess::precess_ecl_coords(
        149.48194_f64.to_radians(),
        1.76549_f64.to_radians(),
        2451545.0,
        1643074.5
    );
    
    assert_eq!((util::round_upto_digits(new_asc.to_degrees(), 3),
                util::round_upto_digits(new_dec.to_degrees(), 3)), (118.704, 1.615));
}

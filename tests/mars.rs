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
fn ephemeris() {
    let JD = 2448935.500638;
    let mut ephemeris = planet::mars::ephemeris(
        JD,

        &planet::mars::north_pol_ecl_coords(time::julian_cent(JD)),

        23.44022_f64.to_radians(),

        angle::deg_frm_dms(0, 0, 15.42).to_radians(),
        angle::deg_frm_dms(0, 0, -1.0).to_radians(),
    );

    ephemeris.De = util::round_upto_digits(ephemeris.De.to_degrees(), 2);
    ephemeris.Ds = util::round_upto_digits(ephemeris.Ds.to_degrees(), 2);
    ephemeris.P = util::round_upto_digits(angle::limit_to_360(ephemeris.P.to_degrees()), 2);
    ephemeris.w = util::round_upto_digits(ephemeris.w.to_degrees(), 1);

    assert_eq!(ephemeris.De, 12.44);
    assert_eq!(ephemeris.Ds, -2.76);
    assert_eq!(ephemeris.P, 347.64);
    assert_eq!(ephemeris.w, 111.5);

    let (h1, m1, s1) = angle::dms_frm_deg(ephemeris.d.to_degrees());
    assert_eq!((h1, m1), (0, 0));
    assert_eq!(util::round_upto_digits(s1, 2), 10.75);

    let (h2, m2, s2) = angle::dms_frm_deg(ephemeris.q.to_degrees());
    assert_eq!((h2, m2), (0, 0));
    assert_eq!(util::round_upto_digits(s2, 2), 1.06);
}

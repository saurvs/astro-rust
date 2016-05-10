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

extern crate astro;

use astro::*;

#[test]
fn eq_coords_frm_motion() {

    let (asc, dec) = star::eq_coords_frm_motion(
        101.286962_f64.to_radians(),
       -16.716108_f64.to_radians(),
        2.64,
       -0.000007773,
        (-0.03847 / 3600.0 as f64).to_radians(),
        (-1.20530 / 3600.0 as f64).to_radians(),
       -1000.0
    );

    let (h1, m1, s1) = angle::hms_frm_deg(asc.to_degrees());
    assert_eq!((h1, m1, util::round_upto_digits(s1, 2)), (6, 45, 47.16));

    let (d2, m2, s2) = angle::dms_frm_deg(dec.to_degrees());
    assert_eq!((d2, m2, util::round_upto_digits(s2, 1)), (-16, -22, -56.0));

}

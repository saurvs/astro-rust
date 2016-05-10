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
pub fn topocent_eq_coords() {

    let eq_point = coords::EqPoint{
        asc: 339.530208_f64.to_radians(),
        dec: -15.771083_f64.to_radians()
    };
    let geograph_point = coords::GeographPoint{
        long: angle::deg_frm_hms(7, 47, 27.0).to_radians(),
        lat: 33.356111_f64.to_radians(),
    };
    let topo_eq_point = parallax::topocent_eq_coords(
        &eq_point,
        angle::deg_frm_dms(0, 0, 23.592).to_radians(),
        &geograph_point,
        1706.0,
        angle::deg_frm_hms(1, 40, 45.0).to_radians()
    );

    let (h, m1, s1) = angle::hms_frm_deg(topo_eq_point.asc.to_degrees());
    assert_eq!((h, m1, util::round_upto_digits(s1, 2)), (22, 38, 8.54));

    let (d, m2, s2) = angle::dms_frm_deg(topo_eq_point.dec.to_degrees());
    assert_eq!((d, m2, util::round_upto_digits(s2, 1)), (-15, -46, -30.0));

}

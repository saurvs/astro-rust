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
fn ring_elements() {
    let mut elements = planet::saturn::ring::elements(
        2448972.50068,

        angle::deg_frm_dms(0, 0, 16.86).to_radians(),
        23.43971_f64.to_radians()
    );

    elements.B = util::round_upto_digits(elements.B.to_degrees(), 3);
    elements.B1 = util::round_upto_digits(elements.B1.to_degrees(), 3);
    elements.P = util::round_upto_digits(angle::limit_to_360(elements.P.to_degrees()), 3);
    elements.deltaU = util::round_upto_digits(elements.deltaU.to_degrees(), 3);

    assert_eq!(elements.B, 16.442);
    assert_eq!(elements.B1, 14.679);
    assert_eq!(elements.P, 6.741);
    assert_eq!(elements.deltaU, 4.198);

    let (d_a, m_a, s_a) = angle::dms_frm_deg(elements.a.to_degrees());
    assert_eq!((d_a, m_a), (0, 0));
    assert_eq!(util::round_upto_digits(s_a, 2), 35.87);

    let (d_b, m_b, s_b) = angle::dms_frm_deg(elements.b.to_degrees());
    assert_eq!((d_b, m_b), (0, 0));
    assert_eq!(util::round_upto_digits(s_b, 2), 10.15);
}

#[test]
#[allow(unused_variables)]
fn moons() {
    let data = [
        (3.102, -0.204, planet::saturn::moon::Moon::Mimas),
        (3.823, 0.318, planet::saturn::moon::Moon::Enceladus),
        (4.027, -1.061, planet::saturn::moon::Moon::Tethys),
        (-5.365, -1.148, planet::saturn::moon::Moon::Dione),
        (-1.122, -3.123, planet::saturn::moon::Moon::Rhea),
        (14.568, 4.738, planet::saturn::moon::Moon::Titan),
        (-18.001, -5.328, planet::saturn::moon::Moon::Hyperion),
        (-48.76/* Meeus gives -48.759 */, 4.137/* Meeus gives 4.136 */, planet::saturn::moon::Moon::Iapetus),
    ];

    for tuple in data.iter() {
        let (X, Y, Z) = planet::saturn::moon::apprnt_rect_coords(
            2451439.50074,
            &tuple.2
        );
        assert_eq!(util::round_upto_digits(X, 3), tuple.0);
        assert_eq!(util::round_upto_digits(Y, 3), tuple.1);
    }
}

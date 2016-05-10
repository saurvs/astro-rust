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

    let mut ephemeris = planet::jupiter::ephemeris(
        2448972.50068,

        23.4402069_f64.to_radians(),

        angle::deg_frm_dms(0, 0, 16.86).to_radians(),
        angle::deg_frm_dms(0, 0, -1.79).to_radians(),
    );

    ephemeris.De = util::round_upto_digits(ephemeris.De.to_degrees(), 2);
    ephemeris.Ds = util::round_upto_digits(ephemeris.Ds.to_degrees(), 2);
    ephemeris.P = util::round_upto_digits(angle::limit_to_360(ephemeris.P.to_degrees()), 2);
    ephemeris.w1 = util::round_upto_digits(ephemeris.w1.to_degrees(), 0);
    ephemeris.w2 = util::round_upto_digits(ephemeris.w2.to_degrees(), 2);

    assert_eq!(ephemeris.De, -2.48);
    assert_eq!(ephemeris.Ds, -2.20);
    assert_eq!(ephemeris.P, 24.80);
    assert_eq!(ephemeris.w1, 268.0);
    assert_eq!(ephemeris.w2, 72.74);
}

#[test]
fn moons() {

    let data = [
        (-3.44, 0.21, planet::jupiter::moon::Moon::Io),
        (7.44, 0.25, planet::jupiter::moon::Moon::Europa),
        (1.24, 0.65, planet::jupiter::moon::Moon::Ganymede),
        (7.08, 1.1, planet::jupiter::moon::Moon::Callisto),
    ];

    for tuple in data.iter() {
        let (X, Y) = planet::jupiter::moon::apprnt_rect_coords(
            2448972.50068,
            &tuple.2
        );

        assert_eq!(util::round_upto_digits(X, 2), tuple.0);
        assert_eq!(util::round_upto_digits(Y, 2), tuple.1);
    }

}

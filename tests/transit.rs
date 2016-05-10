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
#[allow(unused_variables)]
fn time() {

    let eq_point1 = coords::EqPoint{
        asc: 40.68021_f64.to_radians(),
        dec: 18.04761_f64.to_radians()
    };
    let eq_point2 = coords::EqPoint{
        asc: 41.73129_f64.to_radians(),
        dec: 18.44092_f64.to_radians()
    };
    let eq_point3 = coords::EqPoint{
        asc: 42.78204_f64.to_radians(),
        dec: 18.82742_f64.to_radians()
    };

    let geograph_point = coords::GeographPoint{
        long: 71.0833_f64.to_radians(),
        lat: 42.3333_f64.to_radians(),
    };

    let Theta0 = 177.74208_f64.to_radians();

    let deltaT = time::delta_t(1988, 3);

    let (h_rise, m_rise, s_rise) = transit::time(
        &transit::TransitType::Rise,
        &transit::TransitBody::StarOrPlanet,
        &geograph_point,
        &eq_point1,
        &eq_point2,
        &eq_point3,
        Theta0,
        deltaT,
        0.0
    );

    assert_eq!((h_rise, m_rise), (12, 25));

    let (h_transit, m_transit, s_transit) = transit::time(
        &transit::TransitType::Transit,
        &transit::TransitBody::StarOrPlanet,
        &geograph_point,
        &eq_point1,
        &eq_point2,
        &eq_point3,
        Theta0,
        deltaT,
        0.0
    );

    assert_eq!((h_transit, m_transit), (19, 40));

    let (h_set, m_set, s_set) = transit::time(
        &transit::TransitType::Set,
        &transit::TransitBody::StarOrPlanet,
        &geograph_point,
        &eq_point1,
        &eq_point2,
        &eq_point3,
        Theta0,
        deltaT,
        0.0
    );

    assert_eq!((h_set, m_set), (2, 54));

}

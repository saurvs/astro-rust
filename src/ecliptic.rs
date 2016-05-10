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

//! The Earth's ecliptic

use angle;
use std::f64::consts::PI;
use time;

/**
Computes the mean obliquity of the ecliptic using
J. Laskar's formula

# Returns

* `mn_oblq`: Mean obliquity of the ecliptic *| in radians*

The accuracy of `mn_oblq` is estimated to be 0.01 arcseconds
for 1000 years before and after 2000 AD, and a few arcseconds
for 10000 years before and after 2000 AD.

# Arguments

* `JD`: Julian (Ephemeris) day
**/
pub fn mn_oblq_laskar(JD: f64) -> (f64) {

    let u = time::julian_cent(JD) / 100.0;

    Horner_eval!(
         u,
         angle::deg_frm_dms(23, 26, 21.448),
        -angle::deg_frm_dms(0, 0, 4680.93),
        -angle::deg_frm_dms(0, 0, 1.55),
         angle::deg_frm_dms(0, 0, 1999.25),
        -angle::deg_frm_dms(0, 0, 51.38),
        -angle::deg_frm_dms(0, 0, 249.67),
        -angle::deg_frm_dms(0, 0, 39.05),
         angle::deg_frm_dms(0, 0, 7.12),
         angle::deg_frm_dms(0, 0, 27.87),
         angle::deg_frm_dms(0, 0, 5.79),
         angle::deg_frm_dms(0, 0, 2.45)
    ).to_radians()

}

/**
Computes the mean obliquity of the ecliptic using
the IAU formula

# Returns

* `mn_oblq`: Mean obliquity of the ecliptic *| in radians*

The error in `mn_oblq` reaches 1 arcsecond over a period of
2000 years from 2000 AD, and about 10 arcseconds over a period of
4000 years from 2000 AD.

# Arguments

* `JD`: Julian (Ephemeris) day
**/
pub fn mn_oblq_IAU(JD: f64) -> (f64) {

    let u = time::julian_cent(JD) / 100.0;

    Horner_eval!(
         u,
         angle::deg_frm_dms(23, 26, 21.448),
        -angle::deg_frm_dms(0, 0, 46.815),
        -angle::deg_frm_dms(0, 0, 0.00059),
         angle::deg_frm_dms(0, 0, 0.001813)
    ).to_radians()

}

/**
Computes the longitudes of the two ecliptic points on
a horizon on Earth

# Returns

`(long_point_1, long_point_2)`

* `long_point_1`: Longitude of ecliptic point 1 *| in radians*
* `long_point_2`: Longitude of ecliptic point 2 *| in radians*

# Arguments

* `oblq_eclip`  : Obliquity of the ecliptic *| in radians*
* `observer_lat`: The observer's geographical latitude *| in radians*
* `loc_sidreal` : Local sidereal time *| in radians*
**/
pub fn eclip_points_on_hz(oblq_eclip: f64, observer_lat: f64, loc_sidreal: f64) -> (f64, f64) {

    let p = (-loc_sidreal.cos()).atan2 (
        oblq_eclip.sin() * observer_lat.tan() +
        oblq_eclip.cos() * loc_sidreal.sin()
    );

    (p, p + PI)

}

/**
Computes the angle between the ecliptic and a horizon
on Earth

# Returns

* `angle`: Angle between the ecliptic and the horizon *| in radians*

# Arguments

* `oblq_eclip`  : Obliquity of the ecliptic *| in radians*
* `observer_lat`: The observer's geographical latitude *| in radians*
* `loc_sidreal` : Local sidereal time *| in radians*
**/
pub fn angl_betwn_eclip_and_hz(oblq_eclip: f64, observer_lat: f64, loc_sidreal: f64) -> f64 {

    (
        oblq_eclip.cos() * observer_lat.sin() -
        oblq_eclip.sin() * observer_lat.cos() * loc_sidreal.sin()
    ).acos()

}

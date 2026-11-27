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

//! The Sun

use crate::angle;
use crate::time;
use std;
use crate::planet;
use crate::coords;

/**
Computes the Sun's equatorial semidiameter

# Arguments

* `sun_earth_dist`: Sun-Earth distance *| in AU*
**/
#[inline]
pub fn semidiameter(sun_earth_dist: f64) -> f64 {

    angle::deg_frm_dms(0, 0, 959.63) / sun_earth_dist

}

/**
Computes the Sun's geocentric ecliptic position, referred to the mean
equinox of the date

# Returns

`(sun_ecl_point, sun_earth_dist)`

* `sun_ecl_point` : Ecliptic point of the Sun *| in radians*
* `sun_earth_dist`: Sun-Earth distance *| in AU*

# Arguments

* `JD`: Julian (Ephemeris) day
**/
pub fn geocent_ecl_pos(JD: f64) -> (coords::EclPoint, f64) {

    let (L, B, R) = planet::heliocent_coords(&planet::Planet::Earth, JD);

    let ecl_point = coords::EclPoint {
        long: angle::limit_to_two_PI(L + std::f64::consts::PI),
        lat:  angle::limit_to_two_PI(-B)
    };

    (ecl_point, R)

}

/**
Computes the Sun's geocentric ecliptic coordinates converted to the
FK5 system

# Returns

`(ecl_long_FK5, ecl_lat_FK5)`

* `ecl_long_FK5`: Ecliptic longitude of the Sun *| in radians*,
                  converted to the FK5 system
* `ecl_lat_FK5` : Ecliptic latitude of the Sun *| in radians*,
                  converted to the FK5 system

# Arguments

* `JD`      : Julian (Ephemeris) day
* `ecl_long`: Ecliptic longitude of the Sun on `JD`
              *| in radians*, referred to the mean equinox
              of the date
* `ecl_lat` : Ecliptic latitude of the Sun `JD`
              *| in radians*, referred to the mean equinox
              of the date
**/
pub fn ecl_coords_to_FK5(JD: f64, ecl_long: f64, ecl_lat: f64) -> (f64, f64) {

    let ecl_long_FK5 =
        ecl_long
      - angle::deg_frm_dms(0, 0, 0.09033).to_radians();

     let JC = time::julian_cent(JD);
     let lambda1 =
         ecl_long
       - JC * (1.397 + JC*0.00031).to_radians();

    let ecl_lat_FK5 =
        ecl_lat
      + angle::deg_frm_dms(0, 0, 0.03916).to_radians() * (
          lambda1.cos() - lambda1.sin()
        );

    (ecl_long_FK5, ecl_lat_FK5)

}

/**
Computes the Sun's geocentric rectangular coordinates, referred to
the mean equinox of the date

# Returns

`(x, y z)`

* `x`: The X coordinate *| in AU*
* `y`: The Y coordinate *| in AU*
* `z`: The Z coordinate *| in AU*

* The positive x-axis is directed towards the Earth's vernal equinox
(0 degrees longitude)
* The positive y-axis lies in the plane of the Earth's equator and is
directed towards 90 degrees longitude
* The positive z-axis is directed towards the Earth's northern
celestial pole

# Arguments

* `sun_geo_long`: The Sun's geometric longitude *| in radians*,
                      *without* corrections for nutation and abberation
* `sun_geo_lat` : The Sun's geometric latitude *| in radians*,
                     *without* corrections for nutation and abberation
* `sun_rad_vec` : The Sun's geometric radius vector *| in AU*
* `mn_oblq`     : Mean obliquity of the ecliptic
**/
pub fn geocent_rect_coords (

    sun_geo_long : f64,
    sun_geo_lat  : f64,
    sun_rad_vec  : f64,
    mn_oblq      : f64

) -> (f64, f64, f64) {

    let x = sun_rad_vec * sun_geo_lat.cos() * sun_geo_long.cos();

    let y = sun_rad_vec * (
        sun_geo_lat.cos() * sun_geo_long.sin() * mn_oblq.cos()
      - sun_geo_lat.sin() * mn_oblq.sin()
    );

    let z = sun_rad_vec * (
        sun_geo_lat.cos() * sun_geo_long.sin() * mn_oblq.sin()
      + sun_geo_lat.sin() * mn_oblq.cos()
    );

    (x, y, z)

}

/**
Return quantites used in the ephemeris for physical observations of
the Sun

# Returns

`(P, B0, L0)`

* `P` : Position angle of the northern extremity of the axis of
        rotation, measured eastwards from the North point of the
        solar disk *| in radians*
* `B0`: Heliographic latitude of the center of the solar
        disk *| in radians*
* `L0`: Heliographic longitude of the center of the solar
        disk *| in radians*

# Arguments

* `JD`      : Julian (Ephemeris) day
* `app_long`: Apparent longitude of the Sun *| in radians*,
                  including the effect of abberation and *not* that
                  of nutation
* `app_long_with_nut`: Apparent longitude of the Sun *| in radians*,
                  including the effect of abberation *and* that
                  of nutation
* `oblq_eclip`: True obliquity of the ecliptic *| in radians*
**/
pub fn ephemeris (

    JD                : f64,
    app_long          : f64,
    app_long_with_nut : f64,
    oblq_eclip        : f64

) -> (f64, f64, f64) {

    let theta = angle::limit_to_360((JD - 2398220.0) * 360.0/25.38).to_radians();
    let I = 7.25_f64.to_radians();
    let K = (73.6667 + 1.3958333*(JD - 2396758.0)/36525.0).to_radians();

    let z = app_long - K;
    let sin_z = z.sin();
    let cos_z = z.cos();

    let mut x = (-app_long_with_nut.cos() * oblq_eclip.tan()).atan();
    let mut y = (-cos_z * I.tan()).atan();
    x = magnitude_limited_to_less_than_PI(x);
    y = magnitude_limited_to_less_than_PI(y);

    let B_0 = (sin_z * I.sin()).asin();
    let nu = (-sin_z * I.cos()).atan2(-cos_z);
    let L_0 = angle::limit_to_two_PI(nu - theta);

    let P = x + y;

    (P, B_0, L_0)

}

#[inline]
fn magnitude_limited_to_less_than_PI(a: f64) -> f64 {

    let PI_INTO_THREE_BY_TWO = std::f64::consts::PI * 3.0/2.0;

    if a > PI_INTO_THREE_BY_TWO { a - angle::TWO_PI }
    else                        { a }

}

/**
Computes an approximate time for the beginning of a solar synodic
rotation

# Returns

* `JD`: The Julian Day corresponding to the approximate time for
        the beginning of a solar synodic rotation

Between the years 1850 and 2100, `JD` will be in less than 0.002 days
in error.

# Arguments

* `C`: Carrington's synodic rotation number
**/
pub fn synodic_rot(C: i64) -> f64 {

    let M = (281.96 + 26.882476*(C as f64)).to_radians();

    2398140.227 + 27.2752316*(C as f64)
  + 0.1454 * M.sin()
  - 0.0085 * (2.0 * M).sin()
  - 0.0141 * (2.0 * M).cos()

}

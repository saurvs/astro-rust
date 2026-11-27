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

//! Find topocentric coordinates

use crate::angle;
use crate::coords;
use crate::planet;

/**
Computes the equatorial horizontal parallax of a celestial body

# Returns

* `eq_hz_parllx`: Equatorial horizontal parallax of the celestial
                  body *| in radians*

# Arguments

* `dist_to_earth`: The celestial body's distance to the Earth *| in AU*
**/
#[inline]
pub fn eq_hz_parallax(dist_to_earth: f64) -> f64
{
    (angle::deg_frm_dms(0, 0, 8.794).to_radians().sin() / dist_to_earth).asin()
}

/**
Computes the topocentric equatorial coordinates of a celestial body

# Returns

* `topocent_eq_point`: Topocentric equatorial point of the
                      celestial body *| in radians*

# Arguments

* `eq_point`          : Equatorial point of the celestial body
                        *| in radians*
* `eq_hz_parllx`      : Equatorial horizontal parallax of the celestial
                        body *| in radians*
* `geograph_point`    : Geographic point of the observer *| in radians*
* `observer_ht`       : Height of the observer above sea level
                        *| in meters*
* `greenw_sidr`: Sidereal time at Greenwhich *| in radians*
**/
pub fn topocent_eq_coords(eq_point: &coords::EqPoint,
                          eq_hz_parllx: f64,
                          geograph_point: &coords::GeographPoint,
                          observer_ht: f64,
                          greenw_sidr: f64) -> coords::EqPoint
{
    let (rho_sin, rho_cos) = planet::earth::rho_sin_cos_phi (
        geograph_point.lat, observer_ht
    );

    let geocent_hr_angl = coords::hr_angl_frm_observer_long (
        greenw_sidr, geograph_point.long, eq_point.asc
    );

    let eq_hz_parllx_sin = eq_hz_parllx.sin();

    let del_asc = (-rho_cos * eq_hz_parllx_sin * geocent_hr_angl.sin()).atan2(
        eq_point.dec.cos()
      - rho_cos * eq_hz_parllx_sin * geocent_hr_angl.cos()
    );

    let dec_1 = (
        (eq_point.dec.sin() - rho_sin * eq_hz_parllx_sin) * del_asc.cos()
    ).atan2(
        eq_point.dec.cos()
      - rho_cos * eq_hz_parllx_sin * geocent_hr_angl.cos()
    );

    coords::EqPoint {
        asc: eq_point.asc + del_asc,
        dec: dec_1
    }
}

/**
Computes the topocentric ecliptic coordinates of a celestial body

# Returns

`(topocent_ecl_point, topocent_geocent_semidia)`

* `topocent_ecl_point`     : Topocentric ecliptic point of the
                            celestial body *| in radians*
* `topocent_geocent_semidia`: Topocentric semidiameter of the
                            celestial body *| in radians*

# Arguments

* `ecl_point`     : Ecliptic point of the celestial body *| in radians*
* `eq_hz_parllx`  : Equatorial horizontal parallax of the celestial
                    body *| in radians*
* `geograph_point`: Geographic point of the observer *| in radians*
* `observer_ht`   : Height of the observer above sea level *| in meters*
* `loc_sidr`      : Local sidereal time *| in radians*
* `eclip_oblq`    : Obliquity of the ecliptic *| in radians*
* `geocent_semdia` : Geocentric semidiameter of the celestial body *| in radians*
**/
pub fn topopcent_ecl_coords(ecl_point: &coords::EclPoint,
                            eq_hz_parllx: f64,
                            geograph_point: &coords::GeographPoint,
                            observer_ht: f64,
                            loc_sidr: f64,
                            eclip_oblq: f64,
                            geocent_semdia: f64) -> (coords::EclPoint, f64)
{
    let (rho_sin, rho_cos) = planet::earth::rho_sin_cos_phi (
        geograph_point.lat, observer_ht
    );

    let eq_hz_parllx_sin = eq_hz_parllx.sin();
    let loc_sidr_sin = loc_sidr.sin();

    let eclip_oblq_sin = eclip_oblq.sin();
    let eclip_oblq_cos = eclip_oblq.cos();

    let ecl_point_lat_cos = ecl_point.lat.cos();

    let N =
        ecl_point.long.cos() * ecl_point_lat_cos
      - rho_cos * eq_hz_parllx_sin * loc_sidr.cos();

    let ecl_long_1 = (
        ecl_point.long.sin() * ecl_point_lat_cos
      - eq_hz_parllx_sin * (
            rho_sin * eclip_oblq_sin +
            rho_cos * eclip_oblq_cos * loc_sidr_sin
        )
    ).atan2(N);

    let ecl_lat_1 = (
        ecl_long_1.cos() * (
            ecl_point.lat.sin()
          - eq_hz_parllx_sin * (
                rho_sin * eclip_oblq_cos -
                rho_cos * eclip_oblq_sin * loc_sidr_sin
            )
        )
    ).atan2(N);

    let geocent_semdia_1 = (
        ecl_long_1.cos() * ecl_lat_1.cos() * geocent_semdia.sin() / N
    ).asin();

    (
        coords::EclPoint {
            long: angle::limit_to_two_PI(ecl_long_1),
            lat: ecl_lat_1
        },
        geocent_semdia_1
    )
}

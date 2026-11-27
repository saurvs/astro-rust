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

//! Mars

use crate::angle;
use crate::planet;
use crate::time;
use crate::coords;

/**
Returns the equatorial coordinates of Mars's north pole for the epoch
J1950.0

# Returns

* `eq_coords`: Equatorial coordinates of Mars's north pole for the
               epoch J1950.0 *| in radians*
**/
#[inline(always)]
pub fn north_pol_eq_coords_J1950() -> coords::EqPoint {

    coords::EqPoint {
        asc: 317.342_f64.to_radians(),
        dec: 52.7110_f64.to_radians()
    }

}

/**
Returns the equatorial coordinates of Mars's north pole for the epoch
J2000.0

# Returns

* `eq_coords`: Equatorial coordinates of Mars's north pole for the
               epoch J2000.0 *| in radians*
**/
#[inline(always)]
pub fn north_pol_eq_coords_J2000() -> coords::EqPoint {

    coords::EqPoint {
        asc: 317.681_f64.to_radians(),
        dec: 52.8860_f64.to_radians()
    }

}

/**
Computes the ecliptic coordinates of Mars's north pole, referred to
the mean equinox of the date

# Returns

* `ecl_coords`: Ecliptic coordinates of Mars's north pole
                *| in radians*

# Arguments

* `JC`: Julian (Ephemeris) century
**/
#[inline(always)]
pub fn north_pol_ecl_coords(JC: f64) -> coords::EclPoint {

    coords::EclPoint {
        long: (352.9065 + 1.17330*JC).to_radians(),
        lat:  (63.28180 - 0.00394*JC).to_radians()
    }
    
}

/// Holds Mar's ephemeris values for physical observations
#[derive(Debug)]
pub struct Ephemeris {
    /// Mars-centric declination of the Earth
    pub De: f64,
    /// Mars-centric declination of the Sun
    pub Ds: f64,
    /// Geocentric position angle of Mars's northern rotation pole,
    /// or also called the position angle of the axis
    pub P : f64,
    /// Angular amount of the greatest defect of illumination
    pub q : f64,
    /// Longitude of the central meridian, as seen from the Earth
    pub w : f64,
    /// Apparent diameter of Mars
    pub d : f64,
}

/**
Computes quantites used in the ephemeris for physical observations of
Mars

# Returns

* `ephemeris`: Mar's ephemeris. *All angles are in radians*

# Arguments

* `JD`                   : Julian (Ephemeris) day
* `north_pole_ecl_coords`: Ecliptic coordinates of Mars's north pole
                           on `JD` *| in radians*
* `mn_oblq`              : Mean obliquity of the ecliptic on `JD`
                           *| in radians*
* `nut_in_long`          : Nutation in ecliptic longitude on `JD`
                           *| in radians*
* `nut_in_oblq`          : Nutation in obliquity of the ecliptic on
                           `JD` *| in radians*
**/
pub fn ephemeris (

    JD                    : f64,
    north_pole_ecl_coords : &coords::EclPoint,
    mn_oblq               : f64,
    nut_in_long           : f64,
    nut_in_oblq           : f64

) -> Ephemeris {

    let (mut lambda0, beta0) = (north_pole_ecl_coords.long, north_pole_ecl_coords.lat);

    let (l0, b0, R) = planet::heliocent_coords(&planet::Planet::Earth, JD);

    let (mut l, mut b, mut r) = (0.0, 0.0, 0.0);
    let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
    let mut mars_earth_dist = 0.0;
    let mut light_time = 0.0;

    let mut i: u8 = 1;
    while i <= 2 {

        let (new_l, new_b, new_r) = planet::heliocent_coords(&planet::Planet::Mars, JD - light_time);
        l = new_l; b = new_b; r = new_r;

        let (new_x, new_y, new_z) = planet::geocent_ecl_rect_coords(l0, b0, R, l, b, r);
        x = new_x; y = new_y; z = new_z;

        mars_earth_dist = planet::dist_frm_ecl_rect_coords(x, y, z);
        light_time = planet::light_time(mars_earth_dist);

        i += 1;

    }

    let (mut lambda, mut beta) = planet::ecl_coords_frm_ecl_rect_coords(x, y, z);

    let D_e = (
        -beta0.sin() * beta.sin() -
         beta0.cos() * beta.cos() * (lambda0 - lambda).cos()
    ).asin();

    let JC = time::julian_cent(JD);
    let N = (49.5581 + 0.7721*JC).to_radians();

    let l1 = l - (0.00697 / r).to_radians();
    let b1 = b - (0.000225 * (l - N).cos()/r).to_radians();
    let D_s = (
        -beta0.sin() * b1.sin()
       - beta0.cos() * b1.cos() * (lambda0 - l1).cos()
    ).asin();

    let W = angle::limit_to_360(
        11.504
      + 350.89200025 * (JD - light_time - 2433282.5)
    ).to_radians();

    let asc0 = coords::asc_frm_ecl(lambda0, beta0, mn_oblq);
    let dec0 = coords::dec_frm_ecl(lambda0, beta0, mn_oblq);

    let u = y*mn_oblq.cos() - z*mn_oblq.sin();
    let v = y*mn_oblq.sin() + z*mn_oblq.cos();
    let asc = u.atan2(x);
    let dec = v.atan2((x*x + u*u).sqrt());
    let zeta = (
        dec0.sin() * dec.cos() * (asc0 - asc).cos()
      - dec.sin() * dec0.cos()
    ).atan2(dec.cos() * (asc0 - asc).sin());
    let w = W - zeta;

    lambda += 0.005693_f64.to_radians() * (l0 - lambda).cos()/beta.cos();
    beta += 0.005693_f64.to_radians() * (l0 - lambda).sin() * beta.sin();

    lambda += nut_in_long;
    lambda0 += nut_in_long;
    let true_oblq = mn_oblq + nut_in_oblq;

    let asc01 = coords::asc_frm_ecl(lambda0, beta0, true_oblq);
    let dec01 = coords::dec_frm_ecl(lambda0, beta0, true_oblq);
    let asc1 = coords::asc_frm_ecl(lambda, beta, true_oblq);
    let dec1 = coords::dec_frm_ecl(lambda, beta, true_oblq);

    let P = (dec01.cos() * (asc01 - asc1).sin()).atan2 (
        dec01.sin() * dec1.cos()
      - dec01.cos() * dec1.sin() * (asc01 - asc1).cos()
    );

    let d =
        angle::deg_frm_dms(0, 0, 9.36).to_radians()
      / mars_earth_dist;
    let k = planet::illum_frac_frm_dist(r, mars_earth_dist, R);
    let q = (1.0 - k) * d;

    Ephemeris {
        De: D_e,
        Ds: D_s,
        P: P,
        q: q,
        w: w,
        d: d
    }

}

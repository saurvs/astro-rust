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

//! The ring system of Saturn

use crate::angle;
use crate::coords;
use crate::planet;
use crate::time;

#[inline]
pub fn inc(JC: f64) -> f64 {

    (
        28.075216
      - JC * (
            0.012998
          + JC * 0.000004
        )
    ).to_radians()

}

#[inline]
pub fn ascend_node(JC: f64) -> f64 {

    (
        169.50847
      + JC * (
            1.394681
          + JC * 0.000412
        )
    ).to_radians()

}

/// Holds the elements for the ring system of Saturn
#[derive(Debug)]
pub struct Elements {
    /// Saturnicentric latitude of the Earth, referred to the plane
    /// of the ring
    pub B : f64,
    /// Saturnicentric latitude of the Sun, referred to the plane of
    /// the ring
    pub B1 : f64,
    /// Geocentric position angle of the northern semiminor axis of
    /// the apparent ellipse of the ring of the axis
    pub P  : f64,
    /// Difference between Saturnicentric longitudes of the Sun and
    /// the Earth, measured in the plane of the ring
    pub deltaU : f64,
    /// Major axis of the outer edge of the outer ring
    pub a : f64,
    /// Minor axis of the outer edge of the outer ring
    pub b : f64
}

/**
Computes the elements for the ring system of Saturn

# Returns

`elements`: Elements for the ring system. *All angles are in radians*

# Arguments

* `JD`         : Julian (Ephemeris) day
* `nut_in_long`: Nutation in longitude on `JD` *| in radians*
* `true_oblq`  : True obliquity of the ecliptic on `JD` *| in radians*
**/
pub fn elements(JD: f64, nut_in_long: f64, true_oblq: f64) -> Elements {

    let (l0, b0, R) = planet::heliocent_coords(&planet::Planet::Earth, JD);

    let (mut l, mut b, mut r) = (0.0, 0.0, 0.0);
    let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
    let mut saturn_earth_dist = 0.0;
    let mut light_time = 0.0;

    let mut i: u8 = 1;
    while i <= 2 {

        let (new_l, new_b, new_r) =
            planet::heliocent_coords(&planet::Planet::Saturn, JD - light_time);
        l = new_l; b = new_b; r = new_r;

        let (new_x, new_y, new_z) =
            planet::geocent_ecl_rect_coords(l0, b0, R, l, b, r);
        x = new_x; y = new_y; z = new_z;

        saturn_earth_dist = planet::dist_frm_ecl_rect_coords(x, y, z);
        light_time = planet::light_time(saturn_earth_dist);

        i += 1;

    }

    let JC = time::julian_cent(JD);
    let inc = inc(JC);
    let ascend_node = ascend_node(JC);

    let (mut lambda, mut beta) = planet::ecl_coords_frm_ecl_rect_coords(x, y, z);
    let B = (  inc.sin() * beta.cos() * (lambda - ascend_node).sin()
             - inc.cos() * beta.sin()
            ).asin();
    let semi_maj = angle::deg_frm_dms(0, 0, 375.35).to_radians() / saturn_earth_dist;
    let semi_min = semi_maj * B.abs().sin();

    let N = (113.6655 + 0.8771*JC).to_radians();

    let l1 = l - (0.01759/r).to_radians();
    let b1 = b - (0.000764*(l - N).cos()/r).to_radians();

    let B1 = (
        inc.sin() * b1.cos() * (l1 - ascend_node).sin()
        - inc.cos() * b1.sin()
    ).asin();
    let U1 = (
        inc.sin() * b1.sin()
        + inc.cos() * b1.cos() * (l1 - ascend_node).sin()
    ).atan2(
        b1.cos() * (l1 - ascend_node).cos()
    );
    let U2 = (
        inc.sin() * beta.sin()
        + inc.cos() * beta.cos() * (lambda - ascend_node).sin()
    ).atan2(
        beta.cos() * (lambda - ascend_node).cos()
    );
    let deltaU = (U1 - U2).abs();

    let mut lambda0 = ascend_node - 90.0_f64.to_radians();
    let beta0 = 90.0_f64.to_radians() - inc;

    let q = 0.005693_f64.to_radians();
    lambda += q * (l0 - lambda).cos() / beta.cos();
    beta += q * (l0 - lambda).sin() * beta.sin();

    lambda0 += nut_in_long;
    lambda  += nut_in_long;

    let asc0 = coords::asc_frm_ecl(lambda0, beta0, true_oblq);
    let dec0 = coords::dec_frm_ecl(lambda0, beta0, true_oblq);
    let asc = coords::asc_frm_ecl(lambda, beta, true_oblq);
    let dec = coords::dec_frm_ecl(lambda, beta, true_oblq);

    let P = (dec0.cos() * (asc0 - asc).sin())
            .atan2(dec0.sin()*dec.cos() - dec0.cos()*dec.sin()*(asc0 - asc).cos());

    Elements {
        B      : B,
        B1     : B1,
        P      : P,
        deltaU : deltaU,
        a      : semi_maj,
        b      : semi_min
    }

}

#[inline(always)]
pub fn inn_edge_outer_ring(a: f64, b: f64) -> (f64, f64) {

    (a * 0.8801, b * 0.8801)

}

#[inline(always)]
pub fn out_edge_inner_ing(a: f64, b: f64) -> (f64, f64) {

    (a * 0.8599, b * 0.8599)

}

#[inline(always)]
pub fn inn_edge_inner_ring(a: f64, b: f64) -> (f64, f64) {

    (a * 0.665, b * 0.665)

}

#[inline(always)]
pub fn inn_edge_dusk_ring(a: f64, b: f64) -> (f64, f64) {

    (a * 0.5486, b * 0.5486)

}

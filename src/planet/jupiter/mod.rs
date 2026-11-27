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

//! Jupiter

pub mod moon;

use crate::angle;
use crate::nutation;
use crate::planet;
use crate::coords;

/**
Computes Jupiter's equatorial semidiameter

# Returns

* `eq_semidia`: Equatorial semidiameter *| in radians per AU*

# Arguments

* `jup_earth_dist`: Jupiter-Earth distance *| in AU*
**/
#[inline(always)]
pub fn eq_semidiameter(jup_earth_dist: f64) -> f64 {

    angle::deg_frm_dms(0, 0, 98.44).to_radians() / jup_earth_dist

}

/**
Computes Jupiter's polar semidiameter

# Returns

* `pol_semidia`: Polar semidiameter *| in radians per AU*

# Arguments

* `jup_earth_dist`: Jupiter-Earth distance *| in AU*
**/
#[inline(always)]
pub fn pol_semidiameter(jup_earth_dist: f64) -> f64 {

    angle::deg_frm_dms(0, 0, 92.06).to_radians() / jup_earth_dist

}

/// Holds Jupiter's ephemeris values for physical observations
#[derive(Debug)]
pub struct Ephemeris {
    /// Jupiter-centric declination of the Earth
    pub De : f64,
    /// Jupiter-centric declination of the Sun
    pub Ds : f64,
    /// Geocentric position angle of Jupiter's northern
    /// rotation pole, or also called, position angle
    /// of the axis
    pub P  : f64,
    /// Longitude of the central meridian for Rotational System I
    pub w1 : f64,
    /// Longitude of the central meridian for Rotational System II
    pub w2 : f64,
}

/**
Return quantites used in the ephemeris for physical observations
of Jupiter

# Returns

* `ephemeris`: Jupiter's ephemeris. *All angles are in radians*

# Arguments

* `JD`         : Julian (Ephemeris) day
* `mn_oblq`    : Mean obliquity of the ecliptic on `JD` *| in radians*
* `nut_in_long`: Nutation in ecliptic longitude on `JD` *| in radians*
* `nut_in_oblq`: Nutation in obliquity of the ecliptic on `JD` *| in radians*
**/
pub fn ephemeris (

    JD          : f64,
    mn_oblq     : f64,
    nut_in_long : f64,
    nut_in_oblq : f64

) -> Ephemeris {

    let d = JD - 2433282.5;
    let T1 = d / 36525.0;

    let asc0 = (268.0 + 0.1061*T1).to_radians();
    let dec0 = (64.5  - 0.0164*T1).to_radians();

    let W1 = angle::limit_to_360(17.710 + 877.90003539*d).to_radians();
    let W2 = angle::limit_to_360(16.838 + 870.27003539*d).to_radians();

    let (l0, b0, R) = planet::heliocent_coords(&planet::Planet::Earth, JD);

    let (mut l, mut b, mut r) = (0.0, 0.0, 0.0);
    let mut x;
    let mut y;
    let mut z;
    let mut jup_earth_dist = 0.0;
    let mut light_time = 0.0;

    let mut i: u8 = 1;
    while i <= 2 {

        let (new_l, new_b, new_r) = planet::heliocent_coords(&planet::Planet::Jupiter, JD - light_time);
        l = new_l; b = new_b; r = new_r;

        let (new_x, new_y, new_z) = planet::geocent_ecl_rect_coords(l0, b0, R, l, b, r);
        x = new_x; y = new_y; z = new_z;

        jup_earth_dist = planet::dist_frm_ecl_rect_coords(x, y, z);
        light_time = planet::light_time(jup_earth_dist);

        i += 1;

    }

    l -= 0.01299_f64.to_radians()*jup_earth_dist / (r*r);
    let (x, y, z) = planet::geocent_ecl_rect_coords(l0, b0, R, l, b, r);
    jup_earth_dist = planet::dist_frm_ecl_rect_coords(x, y, z);

    let asc_s = (mn_oblq.cos()*l.sin() - mn_oblq.sin()*b.tan()).atan2(l.cos());
    let dec_s = (mn_oblq.cos()*b.sin() + mn_oblq.sin()*b.cos()*l.sin()).asin();

    let D_s = (-dec0.sin()*dec_s.sin() - dec0.cos()*dec_s.cos()*(asc0 - asc_s).cos()).asin();

    let u = y*mn_oblq.cos() - z*mn_oblq.sin();
    let v = y*mn_oblq.sin() + z*mn_oblq.cos();
    let mut asc = u.atan2(x);
    let mut dec = v.atan2((x*x + u*u).sqrt());
    let zeta =
        (dec0.sin()*dec.cos()*(asc0 - asc).cos() - dec.sin()*dec0.cos())
        .atan2(dec.cos()*(asc0 - asc).sin());

    let D_e = (-dec0.sin()*dec.sin() - dec0.cos()*dec.cos()*(asc0 - asc).cos()).asin();

    let mut w1 = angle::limit_to_360(W1.to_degrees() - zeta.to_degrees() - 5.07033*jup_earth_dist);
    let mut w2 = angle::limit_to_360(W2.to_degrees() - zeta.to_degrees() - 5.02626*jup_earth_dist);

    let mut C =
        57.2958 * (2.0*r*jup_earth_dist + R*R - r*r - jup_earth_dist*jup_earth_dist)
      / (4.0 * r * jup_earth_dist);
    if (l - l0).sin() < 0.0 {
        C *= -1.0
    }
    w1 = (w1 + C).to_radians();
    w2 = (w2 + C).to_radians();

    let tru_oblq = mn_oblq + nut_in_oblq;

    let q = 0.005693_f64.to_radians();
    asc += q * (asc.cos()*l0.cos()*tru_oblq.cos() + asc.sin()*l0.sin()) / dec.cos();
    dec += q * (  l0.cos()*tru_oblq.cos()*(tru_oblq.tan()*dec.cos()
                - asc.sin()*asc.cos())
                + asc.cos()*dec.sin()*l0.sin());

    let (asc_nut, dec_nut) = nutation::nutation_in_eq_coords(
        &coords::EqPoint{asc: asc, dec: dec},
        nut_in_long,
        nut_in_oblq,
        tru_oblq
    );
    let asc1 = asc + asc_nut;
    let dec1 = dec + dec_nut;

    let (asc0_nut, dec0_nut) = nutation::nutation_in_eq_coords (
        &coords::EqPoint{asc: asc0, dec: dec0},
        nut_in_long,
        nut_in_oblq,
        tru_oblq
    );
    let asc01 = asc0 + asc0_nut;
    let dec01 = dec0 + dec0_nut;

    let P = (dec01.cos() * (asc01 - asc1).sin())
            .atan2(dec01.sin()*dec1.cos() - dec01.cos()*dec1.sin()*(asc01 - asc1).cos());

    Ephemeris {
        De: D_e,
        Ds: D_s,
        P : P,
        w1: w1,
        w2: w2
    }

}

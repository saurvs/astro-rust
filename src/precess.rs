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

//! Corrections for precession

use crate::angle;
use std;
use crate::time;

/**
Computes annual precession in equatorial coordinates towards a new
epoch

# Returns

`(ann_precess_asc, ann_precess_dec)`

* `ann_precess_asc`: Annual precession in right ascension towards
                     the new epoch *| in radians*
* `ann_precess_dec`: Annual precession in declination towards
                     the new epoch *| in radians*

In the case of a star, the precession returned here does not take
into account it's annual proper motion. The annual proper motion
in right ascension and declination must simply be added to the
corresponding annual precession returned here in order to reduce
the equatorial coordinates to the new epoch.

# Arguments

* `asc`: Right ascension for the old epoch *| in radians*
* `dec`: Declination for the old epoch *| in radians*. Shouldn't
         be too close to the celestial poles of the Earth.
* `JD` : Julian (Ephemeris) day corresponding to the new epoch.
         Shouldn't be more than a few hundred years away from
         the old epoch.
**/
pub fn annual_precess(asc: f64, dec: f64, JD: f64) -> (f64, f64)
{
    let JC = time::julian_cent(JD);

    let m = (
          angle::deg_frm_hms(0, 0, 3.07496)
        + angle::deg_frm_hms(0, 0, 0.00186) * JC
    ).to_radians();

    let n = (
          angle::deg_frm_hms(0, 0, 1.33621)
        - angle::deg_frm_hms(0, 0, 0.00057) * JC
    ).to_radians();

    (m + n*asc.sin()*dec.tan(), n * asc.cos())
}

/**
Computes equatorial coordinates reduced to a different epoch

# Returns

`(new_asc, new_dec)`

* `new_asc`: Right ascension in the new epoch *| in radians*
* `new_dec`: Declination in the new epoch *| in radians*

# Arguments

* `old_asc`: Right ascension in the old epoch *| in radians*,
             referred to the FK5 system
* `old_dec`: Declination in the old epoch *| in radians*,
             referred to the FK5 system
* `JD1`    : Julian (Ephemeris) day corresponding to the old epoch
* `JD2`    : Julian (Ephemeris) day corresponding to the new epoch
**/
pub fn precess_eq_coords(old_asc: f64,
                         old_dec: f64,
                         JD1: f64,
                         JD2: f64) -> (f64, f64)
{
    let T = time::julian_cent(JD1);
    let t = (JD2 - JD1) / 36525.0;

    let x = t * (angle::deg_frm_dms(0, 0, 2306.2181) +
                  T * (angle::deg_frm_dms(0, 0, 1.39656) -
                       T*angle::deg_frm_dms(0, 0, 0.000139)));
    let xi = (x + t*t*((angle::deg_frm_dms(0, 0, 0.30188) -
                       T*angle::deg_frm_dms(0, 0, 0.000344)) +
                      t*angle::deg_frm_dms(0, 0, 0.017998))).to_radians();

    let zeta = (x + t*t*((angle::deg_frm_dms(0, 0, 1.09468) -
                       T*angle::deg_frm_dms(0, 0, 0.000066)) +
                      t*angle::deg_frm_dms(0, 0, 0.018203))).to_radians();

    let y = T * angle::deg_frm_dms(0, 0, 0.000217);
    let theta = (t * (angle::deg_frm_dms(0, 0, 2004.3109) +
                   T * (angle::deg_frm_dms(0, 0, 0.8533) - y) -
                  t * ((angle::deg_frm_dms(0, 0, 0.42665) + y) +
                       t*angle::deg_frm_dms(0, 0, 0.041833)))).to_radians();

    let A = old_dec.cos() * (old_asc + xi).sin();

    let B =
          theta.cos() * old_dec.cos() * (old_asc + xi).cos()
        - theta.sin() * old_dec.sin();

    let C =
          theta.sin() * old_dec.cos() * (old_asc + xi).cos()
        + theta.cos() * old_dec.sin();

    (A.atan2(B) + zeta, C.asin())
}

/**
Computes equatorial coordinates, from coordinates referred to the
FK4 system, reduced to a different epoch

# Returns

`(new_asc, new_dec)`

* `new_asc`: Right ascension in the new epoch *| in radians*
* `new_dec`: Declination in the new epoch *| in radians*

# Arguments

* `old_asc`: Right ascension in the old epoch *| in radians*,
             referred to the FK4 system
* `old_dec`: Declination for in old epoch *| in radians*,
             referred to the FK4 system
* `JD1`    : Julian (Ephemeris) day corresponding to the old epoch
* `JD2`    : Julian (Ephemeris) day corresponding to the new epoch
**/
pub fn precess_eq_coords_FK5(old_asc: f64,
                             old_dec: f64,
                             JD1: f64,
                             JD2: f64) -> (f64, f64)
{
    let T = (JD1 - 2415020.3135) / 36524.2199;
    let t = (JD2 - JD1) / 36524.2199;

    let xi = t * (
          angle::deg_frm_dms(0, 0, 2304.25) + T*angle::deg_frm_dms(0, 0, 1.396)
        + t * (
            angle::deg_frm_dms(0, 0, 0.302) + t*angle::deg_frm_dms(0, 0, 0.018)
        )
    );

    let zeta =
        xi
      + t * t * (
          angle::deg_frm_dms(0, 0, 0.791) + t*angle::deg_frm_dms(0, 0, 0.001)
        );

    let theta = t * (
        angle::deg_frm_dms(0, 0, 2004.682)
      - T * angle::deg_frm_dms(0, 0, 0.853)
      - t * (
            angle::deg_frm_dms(0, 0, 0.426)
          + t * angle::deg_frm_dms(0, 0, 0.042)
        )
    );

    let A = old_dec.cos() * (old_asc + xi).sin();

    let B =
          theta.cos() * old_dec.cos() * (old_asc + xi).cos()
        - theta.sin() * old_dec.sin();

    let C =
          theta.sin() * old_dec.cos() * (old_asc + xi).cos()
        + theta.cos() * old_dec.sin();

    (A.atan2(B) + zeta, C.asin())
}

/**
Computes ecliptic coordinates reduced to a different epoch

# Returns

`(new_long, new_lat)`

* `new_long`: Ecliptic longitude in the new epoch *| in radians*
* `new_lat` : Ecliptic latitude in the new epoch *| in radians*

# Arguments

* `old_long`: Ecliptic longitude in the old epoch *| in radians*
* `old_lat` : Ecliptic latitude in the old epoch *| in radians*
* `JD_old`  : Julian (Ephemeris) day corresponding to the old epoch
* `JD_new`  : Julian (Ephemeris) day corresponding to the new epoch
**/
pub fn precess_ecl_coords(old_long: f64,
                          old_lat: f64,
                          JD_old: f64,
                          JD_new: f64) -> (f64, f64)
{
    let T = time::julian_cent(JD_old);
    let t = (JD_new - JD_old) / 36525.0;

    let (nu, Pi, rho) = angles_for_ecl_change(t, T);

    let A =
          nu.cos() * old_lat.cos() * (Pi - old_long).sin()
        - nu.sin() * old_lat.sin();

    let B = old_lat.cos() * (Pi - old_long).cos();

    let C =
          nu.cos() * old_lat.sin()
        + nu.sin() * old_lat.cos() * (Pi - old_long).sin();

    let new_long = rho + Pi - A.atan2(B);
    let new_lat = C.asin();

    (new_long, new_lat)
}

#[inline]
fn angles_for_ecl_change(t: f64, T: f64) -> (f64, f64, f64)
{
    let x = T * angle::deg_frm_dms(0, 0, 0.000598);
    let nu = (
        t * (
            angle::deg_frm_dms(0, 0, 47.0029)
          - T * (angle::deg_frm_dms(0, 0, 0.06603) - x)
          + t * (
                (angle::deg_frm_dms(0, 0, -0.03302) + x)
              + t * angle::deg_frm_dms(0, 0, 0.00006)
            )
        )
    ).to_radians();

    let Pi = (
        174.876384
      + T * (
            angle::deg_frm_dms(0, 0, 3289.4789)
          + T * angle::deg_frm_dms(0, 0, 0.60622)
        )
      - t * (
                (
                    angle::deg_frm_dms(0, 0, 869.8089)
                  + T * angle::deg_frm_dms(0, 0, 0.50491)
                )
              - t * angle::deg_frm_dms(0, 0, 0.03536)
        )
    ).to_radians();

    let y = T * angle::deg_frm_dms(0, 0, 0.000042);
    let rho = (
        t * (
            angle::deg_frm_dms(0, 0, 5029.0966)
          + T * (angle::deg_frm_dms(0, 0, 2.22226) - y)
          + t * (
                (angle::deg_frm_dms(0, 0, 1.11113) - y)
              - t * angle::deg_frm_dms(0, 0, 0.000006)
            )
        )
    ).to_radians();

    (nu, Pi, rho)
}

/**
Computes orbital elements reduced to a different equinox

# Returns

`(new_inc, new_arg_perih, new_long_ascend_node)`

* `new_inc`             : Inclination in the new
                          equinox *| in radians*
* `new_arg_perih`       : Argument of perihelion in the
                          new equinox *| in radians*
* `new_long_ascend_node`: Longitude of ascending node in
                          the new equinox *| in radians*

# Arguments

* `old_inc`             : Inclination in the old equinox *| in radians*
* `old_arg_perih`       : Argument of perihelion in the old
                          equinox *| in radians*
* `old_long_ascend_node`: Longitude of ascending node in the old
                          equinox *| in radians*
* `JD1`                 : Julian (Ephemeris) day corresponding to the old
                          equinox
* `JD2`                 : Julian (Ephemeris) day corresponding to the new
                          equinox
**/
pub fn precess_orb_elements(old_inc: f64,
                            old_arg_perih: f64,
                            old_long_ascend_node: f64,
                            JD1: f64,
                            JD2: f64) -> (f64, f64, f64)
{
    let T = time::julian_cent(JD1);
    let t = (JD2 - JD1) / 36525.0;

    let (nu, Pi, rho) = angles_for_ecl_change(t, T);

    let new_inc;
    let new_long_ascend_node;
    let phi = Pi + rho;

    if old_inc == 0.0 {
        new_inc = nu;
        new_long_ascend_node = phi + std::f64::consts::PI;
    } else {
        let A = old_inc.sin() * (old_long_ascend_node - Pi).sin();

        let B =
             -nu.sin() * old_inc.cos()
            + nu.cos() * old_inc.sin() * (old_long_ascend_node - Pi).cos();

        new_inc = (A*A + B*B).sqrt().asin();
        new_long_ascend_node = phi + A.atan2(B);
    }

    let delta_w =
        (-nu.sin() * (old_long_ascend_node - Pi).sin()/new_inc.sin()).asin();

    (new_inc, old_arg_perih + delta_w, new_long_ascend_node)
}

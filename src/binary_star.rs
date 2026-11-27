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

//! Binary stars

use crate::angle;

/**
Computes mean annual motion of companion star

* `P`: Period of revolution of binary star
       (*mean solar year*)
**/

#[inline]
pub fn mn_ann_motion_of_compan(P: f64) -> f64
{
    angle::TWO_PI / P
}

/**
Computes mean anomaly of companion star

# Arguments

* `n`: Mean annual motion of companion star
* `t`: Current time, given as year with
       decimals (eg: 1945.62)
* `T`: Time of periastron passage, given as
       a year with decimals (eg: 1945.62)
**/

#[inline]
pub fn mn_anom_of_compan(n: f64, t: f64, T: f64) -> f64
{
    n * (t - T)
}

/**
Computes radius vector of a binary star

# Arguments

* `a`       : Apparent semimajor axis
* `e`       : Eccentricity of true orbit
* `ecc_anom`: Eccentric anomaly of binary star
**/

#[inline]
pub fn rad_vec(a: f64, e: f64, ecc_anom: f64) -> f64
{
    a * (1.0 - e*ecc_anom.cos())
}

/**
Computes true anomaly of a binary star

# Arguments

* `e`       : Eccentricity of true orbit
* `ecc_anom`: Eccentric anomaly of binary star
**/

#[inline]
pub fn true_anom(e: f64, ecc_anom: f64) -> f64
{
    2.0 * (
        ((1.0 + e)/(1.0 - e)).sqrt() * (ecc_anom / 2.0).tan()
    ).atan()
}

/**
Computes apparent position angle of a binary star

# Arguments

* `asc_node_coords`: Position angle of ascending node
* `true_anom`   : True anomaly of binary star
* `w`           : Longitude of periastron
* `i`           : Inclination of true orbit to a
                  plane at right angles to the
                  line of sight *| in radians*
**/
pub fn apprnt_coords_angl(asc_node_coords: f64, true_anom: f64,
                          w: f64, i: f64) -> f64
{
    let x = (
        (true_anom + w).sin() * i.cos()
    ).atan2((true_anom + w).cos());

    angle::limit_to_two_PI(x + asc_node_coords)
}

/**
Computes angular separation of a binary star

# Arguments

* `rad_vec`  : Radius vector of a binary star
* `true_anom`: True anomaly of a binary star
* `w`        : Longitude of periastron
* `i`        : Inclination of true orbit to a
               plane at right angles to the
               line of sight *| in radians*
**/
pub fn anglr_sepr(rad_vec: f64, true_anom: f64, w: f64, i: f64) -> f64
{
    rad_vec * (
        ((true_anom + w).sin() * i.cos()).powi(2)
      + (true_anom + w).cos().powi(2)
    ).sqrt()
}

/**
Computes eccentricity of an apparent orbit

# Arguments

* `e`: Eccentricity of the true orbit
* `i`: Inclination of true orbit
       to plane at right angles to line
       of sight *| in radians*
* `w`: Longitude of periastron
**/
pub fn ecc_of_apprnt_orb(e: f64, w: f64, i: f64) -> f64
{
    let i_cos = i.cos();
    let e_w_cos = e * w.cos();
    let e_w_cos_sqr = e_w_cos * e_w_cos;

    let a = (1.0 - e_w_cos_sqr) * i_cos * i_cos;
    let b = e * w.sin() * e_w_cos * i_cos;
    let c = 1.0 - e_w_cos_sqr;
    let d = ((a - c)*(a - c) + 4.0*b*b).sqrt();

    ((2.0 * d) / (a + c + d)).sqrt()
}

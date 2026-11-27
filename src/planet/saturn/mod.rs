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

//! Saturn

use crate::angle;

pub mod moon;
pub mod ring;

/**
Computes Saturn's apparent magnitude using G. Muller's formula

# Returns

* `app_mag`: Apparent magnitude of Saturn *| in radians*

# Arguments

* `delta` : Saturn-Earth distance *| in AU*
* `r`     : Saturn-Sun distance *| in AU*
* `deltaU`: Difference between Saturnicentric longitudes of
            the Sun and the Earth, measured in the plane of
            Saturn's ring *| in radians*
* `B`     : Saturnicentric latitude of the Earth *| in radians*
**/
pub fn apprnt_mag_muller(delta: f64, r: f64, delU: f64, B: f64) -> f64 {

    - 8.68
    + 5.0*(r*delta).log10()
    + 0.044*delU.abs()
    - 2.6*B.abs().sin()
    + 1.25*B.sin().powi(2)

}

/**
Computes Saturn's apparent magnitude using the Astronomical
Almanac's formula adopted in 1984

# Returns

* `app_mag`: Apparent magnitude of Saturn *| in radians*

# Arguments

* `delta` : Saturn-Earth distance *| in AU*
* `r`     : Saturn-Sun distance *| in AU*
* `deltaU`: Difference between Saturnicentric longitudes of
            the Sun and the Earth, measured in the plane of
            Saturn's ring *| in radians*
* `B`     : Saturnicentric latitude of the Earth *| in radians*
**/
pub fn apprnt_mag_84(delta: f64, r: f64, delU: f64, B: f64) -> f64 {

    - 8.88
    + 5.0*(r*delta).log10()
    + 0.044*delU.abs()
    - 2.6*B.abs().sin()
    + 1.25*B.sin().powi(2)

}

#[inline(always)]
fn equatorial_unit_semidiameter() -> f64 {

    angle::deg_frm_dms(0, 0, 82.73).to_radians()

}

#[inline(always)]
fn polar_unit_semidiameter() -> f64 {

    angle::deg_frm_dms(0, 0, 73.82).to_radians()

}

/**
Computes Saturn's polar semidiameter

# Returns

* `pol_semidiameter`: Polar semidiameter *| in radians per AU*

# Arguments

* `saturn_earth_dist`: Saturn-Earth distance *| in AU*
* `earth_lat`        : Saturnicentric latitude of Earth *| in radians*
**/
pub fn pol_semidiameter(saturn_earth_dist: f64, earth_lat: f64) -> f64 {

    let a = equatorial_unit_semidiameter();
    let b = polar_unit_semidiameter();
    let k = 1.0 - (b/a).powi(2);

    (a / saturn_earth_dist) * (1.0 - k*earth_lat.cos().powi(2)).sqrt()

}

/**
Computes Saturn's equatorial semidiameter

# Returns

* `eq_semidiameter`: Equatorial semidiameter *| in radians per AU*

# Arguments

* `saturn_earth_dist`: Saturn-Earth distance *| in AU*
**/
#[inline(always)]
pub fn eq_semidiameter(saturn_earth_dist: f64) -> f64 {

    equatorial_unit_semidiameter() / saturn_earth_dist

}

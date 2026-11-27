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

//! Corrections for atmospheric refraction

use crate::angle;
use std::f64::consts::PI;

/**
Computes the refraction term for true altitudes greater than 15
degrees

# Returns

* `refrac_term`: The refraction term *| in radians*, that needs to be
                 subtracted from the apparent altitude to get the
                 true altitude

# Arguments

* `apprnt_alt`: Apparent altitude *| in radians*
**/
pub fn refrac_frm_apprnt_alt_15(apprnt_alt: f64) -> f64
{
    let x =
        angle::deg_frm_dms(0, 0, 0.0668).to_radians()
      * (PI - apprnt_alt).tan();

      angle::deg_frm_dms(0, 0, 58.294).to_radians() * (PI - apprnt_alt).tan()
    - x * x * x
}

/**
Computes the refraction term for apparent altitudes greater than 15
degrees

# Returns

* `refrac_term`: The refraction term *| in radians*, that needs to be
                 added to the true altitude to get the apparent
                 altitude

# Arguments

* `true_alt`: True altitude *| in radians*
**/
pub fn refrac_frm_true_alt_15(true_alt: f64) -> f64
{
    let x =
        angle::deg_frm_dms(0, 0, 0.0824).to_radians()
      * (PI - true_alt).tan();

      angle::deg_frm_dms(0, 0, 58.276).to_radians() * (PI - true_alt).tan()
    - x * x * x
}

/**
Computes the refraction term for true altitude

# Returns

* `refrac_term`: The refraction term *| in radians*, that needs to be
                 subtracted from the apparent altitude to get the
                 rue altitude

The accuracy of `refrac_term` is upto 0.07 arcminutes.

# Arguments

* `apprnt_alt`: Apparent altitude *| in radians*
**/
pub fn refrac_frm_apprnt_alt(apprnt_alt: f64) -> f64
{
    if apprnt_alt == PI { 0.0 }
    else {
        let apprnt_alt_deg = apprnt_alt.to_degrees();
        let a = apprnt_alt_deg + 7.31/(apprnt_alt_deg + 4.4);
        let R = 1.0 / a.to_radians().tan();

        (R / 60.0).to_radians()
    }
}

/**
Computes the refraction term for apparent altitude

This function is consistent with `refrac_frm_apprnt_alt()` to
within 4 arcseconds.

# Returns

* `refrac_term`: The refraction term *| in radians*, that needs to be
                 added to the true altitude to get the apparent
                 altitude

# Arguments

* `true_alt`: True altitude *| in radians*
**/
pub fn refrac_frm_true_alt(true_alt: f64) -> f64
{
    if true_alt == PI { 0.0 }
    else {
        let true_alt_deg = true_alt.to_degrees();
        let a = true_alt_deg + 10.3/(true_alt_deg + 5.11);
        let R = 1.02 / a.to_radians().tan();

        (R / 60.0).to_radians()
    }
}

/**
Computes the refraction term modifier for local pressure

# Returns

* `refrac_term_modifier`: The value that needs to be multiplied by
                          the refraction term to account for local
                          pressure

# Arguments

* `pressure`: Local pressure *| in millibars*
**/
#[inline(always)]
pub fn refrac_by_pressr(pressure: f64) -> f64
{
    pressure / 1010.0
}

/**
Computes the refraction term modifier for local temperature

# Returns

* `refrac_term_modifier`: The value that needs to be multiplied by
                          the refraction term to account for local
                          temperature

# Arguments

* `temp`: Local temperature *| in kelvins*
**/
#[inline(always)]
pub fn refrac_by_temp(temp: f64) -> f64
{
    283.0 / temp
}

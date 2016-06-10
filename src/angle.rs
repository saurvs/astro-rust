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

//! Angles for astronomy

use std;

pub const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

/**
Computes the angular separation between two angular points

# Arguments

* `p1a1`: Angle 1 of point 1 *| in radians*
* `p1a2`: Angle 2 of point 1 *| in radians*
* `p2a1`: Angle 1 of point 2 *| in radians*
* `p2a2`: Angle 2 of point 2 *| in radians*

Angle 1 may be right ascension or longitude.
Angle 2 may be declination or latitude.
**/
#[inline]
pub fn anglr_sepr(p1a1: f64, p1a2: f64, p2a1: f64, p2a2: f64) -> f64
{
    (
          p1a2.sin() * p2a2.sin()
        + p1a2.cos() * p2a2.cos() * (p1a1 - p2a1).cos()
    ).acos()
}

/**
Computes an angle in degrees with decimals, from an angle
expressed in degrees, arcminutes and arcseconds

# Returns

* `deg`: Angle in degrees with decimals

# Arguments

* `deg`: Degrees
* `min`: Arcminutes
* `sec`: Arcseconds
**/
#[inline]
pub fn deg_frm_dms(deg: i64, min: i64, sec: f64) -> f64
{
    let (M, S) =
        if deg < 0 { (-min.abs(), -sec.abs()) }
        else       { (min, sec) };

    (deg as f64) + (M as f64)/60.0 + S/3600.0
}

/**
Computes an angle expressed in degrees, arcminutes and
arcseconds, from an angle in degrees with decimals

# Returns

`(deg, min, sec)`

* `deg`: Degrees
* `min`: Arcminutes
* `sec`: Arcseconds

# Arguments

* `deg`: Angle in degrees with decimals
**/
#[inline]
pub fn dms_frm_deg(deg: f64) -> (i64, i64, f64)
{
    let degree = deg as i64;

    let minutes = (deg - (degree as f64)) * 60.0;
    let minute = minutes as i64;

    let seconds = (minutes - (minute as f64)) * 60.0;

    (degree, minute, seconds)
}

/**
Computes an angle in degrees with decimals, from an angle
expressed in hours, minutes and seconds

# Arguments

* `hours`: Hours
* `min`: Minutes
* `sec`: Seconds
**/
#[inline]
pub fn deg_frm_hms(hour: i64, min: i64, sec: f64) -> f64
{
    15.0 * ((hour as f64) + (min as f64)/60.0 + sec/3600.0)
}

/**
Computes an angle expressed in hours, minutes and
seconds, from an angle in degrees with decimals

# Returns

`(deg, min, sec)`

* `hour`: Hours
* `min`: Minutes
* `sec`: Seconds

# Arguments

* `deg`: Angle in degrees with decimals
**/
#[inline]
pub fn hms_frm_deg(deg: f64) -> (i64, i64, f64)
{
    let hours = deg / 15.0;
    let hour = hours as i64;

    let minutes = (hours - (hour as f64)) * 60.0;
    let minute = minutes as i64;

    let seconds = (minutes - (minute as f64)) * 60.0;

    (hour, minute, seconds)
}

/**
Computes the equivalent angle in [0, 360] degree range

# Arguments

* `angl`: Angle *| in degrees*
**/
#[inline]
pub fn limit_to_360(angl: f64) -> f64
{
    let n = (angl / 360.0) as i64;
    let limited_angl = angl - (360.0 * (n as f64));

    if limited_angl < 0.0 { limited_angl + 360.0 }
    else                  { limited_angl }
}

/**
Computes the equivalent angle in [0, 2π] radian range

# Arguments

* `angl`: Angle *| in radians*
**/
#[inline]
pub fn limit_to_two_PI(angl: f64) -> f64
{
    let n = (angl / TWO_PI) as i64;
    let limited_angl = angl - (TWO_PI * (n as f64));

    if limited_angl < 0.0 { limited_angl + TWO_PI }
    else                  { limited_angl }
}

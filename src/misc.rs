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

//! Miscellaneous routines

/**
Computes the parallactic angle of a celestial body

# Returns

* `parallac_angle`: Parallactic angle of the celestial
                    body *| in radians*

# Arguments

* `observer_lat`: The observer's geographical latitude *| in radians*
* `hour_angle`  : Local hour angle *| in radians*
* `dec`         : Declination of the celestial body *| in radians*
**/
#[inline]
pub fn parllc_angl(observer_lat: f64, hour_angle: f64, dec: f64) -> f64
{
    hour_angle.sin().atan2 (
        observer_lat.tan() * dec.cos()
      - hour_angle.cos() * dec.sin()
    )
}

/**
Computes the parallactic angle of a celestial body on the
horizon

# Returns

* `parallac_angle`: Parallactic angle of the celestial body
                    on the horizon *| in radians*

# Arguments

* `observer_lat`: The observer's geographical latitude *| in radians*
* `dec`         : Declination of the celestial body *| in radians*
**/
#[inline]
pub fn parllc_angl_on_hz(observer_lat: f64, dec: f64) -> f64
{
    (observer_lat.sin() / dec.cos()).acos()
}

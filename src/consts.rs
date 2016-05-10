/*
Copyright (c) 2015 Saurav Sachidanand

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

//! Some useful physical constants

/// Returns the Gaussian gravitational constant *| in radians*
#[inline(always)]
pub fn gauss_grav() -> f64 {

    0.01720209895

}

/// Returns the speed of light in vaccum *| in meters per second*
#[inline(always)]
pub fn speed_of_light() -> f64 {

    299_792_458.0

}

/// Returns the Earth-Moon mass ratio
#[inline(always)]
pub fn earth_moon_mass_ratio() -> f64 {

    81.3007

}

/// Returns the Sun-Earth mass ratio
#[inline(always)]
pub fn sun_earth_mass_ratio() -> f64 {

    332_946.0

}

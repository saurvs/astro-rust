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

//! Diameters of asteroids

/**
Computes the diameter of an asteroid

# Returns

* `diameter`: Diameter of the asteroid *| in meters*

# Arguments

* `abs_mag`: Absolute magnitude of the asteroid
* `albedo` : Reflective power of the asteroid
**/
#[inline]
pub fn diameter(abs_mag: f64, albedo: f64) -> f64 {

    1000.0 * 10_f64.powf(3.12 - abs_mag/5.0 - 0.217147*albedo.log10())

}

/**
Computes the apparent diameter of an asteroid

# Returns

* `apparent_diameter`: Apparent diameter of the asteroid *| in meters*

# Arguments

* `true_diameter`      : True diameter of the asteroid *| in kilometers*
* `asteroid_earth_dist`: Asteroid-Earth distance *| in AU*
**/
#[inline]
pub fn apparent_diameter(true_diameter: f64, asteroid_earth_dist: f64) -> f64 {

    1.3788 * true_diameter/asteroid_earth_dist
    
}

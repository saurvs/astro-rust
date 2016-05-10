/*
Copyright (c) 2016 Saurav Sachidanand

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

//! World Geodetic System 1972

/// Equatorial radius of the Earth (Semimajor axis) *| in meters*
pub const EQUATORIAL_RADIUS: f64 = 6378135.0;

/// Polar radius of the Earth (Semiminor axis) *| in meters*
pub const POLAR_RADIUS: f64 = EQUATORIAL_RADIUS * (1.0 - FLATTENING);

/// Flattening
pub const FLATTENING: f64 = 1.0 / 298.26;

/// Angular velocity of the Earth *| in radians / second*
pub const ANGULAR_VELOCITY: f64 = 7.292115147e-5;

/// Gravitational constant (including the mass of the Earth's
/// atmosphere) *| in meters^3 / second^2*
pub const GRAV_CONST: f64 = 3.986008e+14;

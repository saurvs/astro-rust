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

//! Corrections for nutation

use crate::angle;
use crate::time;
use crate::coords;

/**
Computes nutation in ecliptic longitude and obliquity

# Returns

`(nut_in_long, nut_in_oblq)`

* `nut_in_long`: Nutation in ecliptic longitude *| in radians*
* `nut_in_oblq`: Nutation in obliquity of the ecliptic *| in radians*

# Arguments

`JD`: Julian (Ephemeris) day
**/
pub fn nutation(JD: f64) -> (f64, f64)
{
    struct terms(i8, i8, i8, i8, i8, i32, i32, i32, i16);
    let terms_for_nutation = [
        terms( 0,  0,  0,  0,  1, -171996, -1742, 92025,  89),
        terms(-2,  0,  0,  2,  2,  -13187,   -16,  5736, -31),
        terms( 0,  0,  0,  2,  2,   -2274,    -2,   977,  -5),
        terms( 0,  0,  0,  0,  2,    2062,     2,  -895,   5),
        terms( 0,  1,  0,  0,  0,    1426,   -34,    54,  -1),
        terms( 0,  0,  1,  0,  0,     712,     1,    -7,   0),
        terms(-2,  1,  0,  2,  2,    -517,    12,   224,  -6),
        terms( 0,  0,  0,  2,  1,    -386,    -4,   200,   0),
        terms( 0,  0,  1,  2,  2,    -301,     0,   129,  -1),
        terms(-2, -1,  0,  2,  2,     217,    -5,   -95,   3),
        terms(-2,  0,  1,  0,  0,    -158,     0,     0,   0),
        terms(-2,  0,  0,  2,  1,     129,     1,   -70,   0),
        terms( 0,  0, -1,  2,  2,     123,     0,   -53,   0),
        terms( 2,  0,  0,  0,  0,      63,     0,     0,   0),
        terms( 0,  0,  1,  0,  1,      63,     1,   -33,   0),
        terms( 2,  0, -1,  2,  2,     -59,     0,    26,   0),
        terms( 0,  0, -1,  0,  1,     -58,    -1,    32,   0),
        terms( 0,  0,  1,  2,  1,     -51,     0,    27,   0),
        terms(-2,  0,  2,  0,  0,      48,     0,     0,   0),
        terms( 0,  0, -2,  2,  1,      46,     0,   -24,   0),
        terms( 2,  0,  0,  2,  2,     -38,     0,    16,   0),
        terms( 0,  0,  2,  2,  2,     -31,     0,    13,   0),
        terms( 0,  0,  2,  0,  0,      29,     0,     0,   0),
        terms(-2,  0,  1,  2,  2,      29,     0,   -12,   0),
        terms( 0,  0,  0,  2,  0,      26,     0,     0,   0),
        terms(-2,  0,  0,  2,  0,     -22,     0,     0,   0),
        terms( 0,  0, -1,  2,  1,      21,     0,   -10,   0),
        terms( 0,  2,  0,  0,  0,      17,    -1,     0,   0),
        terms( 2,  0, -1,  0,  1,      16,     0,    -8,   0),
        terms(-2,  2,  0,  2,  2,     -16,     1,     7,   0),
        terms( 0,  1,  0,  0,  1,     -15,     0,     9,   0),
        terms(-2,  0,  1,  0,  1,     -13,     0,     7,   0),
        terms( 0, -1,  0,  0,  1,     -12,     0,     6,   0),
        terms( 0,  0,  2, -2,  0,      11,     0,     0,   0),
        terms( 2,  0, -1,  2,  1,     -10,     0,     5,   0),
        terms( 2,  0,  1,  2,  2,      -8,     0,     3,   0),
        terms( 0,  1,  0,  2,  2,       7,     0,    -3,   0),
        terms(-2,  1,  1,  0,  0,      -7,     0,     0,   0),
        terms( 0, -1,  0,  2,  2,      -7,     0,     3,   0),
        terms( 2,  0,  0,  2,  1,      -7,     0,     3,   0),
        terms( 2,  0,  1,  0,  0,       6,     0,     0,   0),
        terms(-2,  0,  2,  2,  2,       6,     0,    -3,   0),
        terms(-2,  0,  1,  2,  1,       6,     0,    -3,   0),
        terms( 2,  0, -2,  0,  1,      -6,     0,     3,   0),
        terms( 2,  0,  0,  0,  1,      -6,     0,     3,   0),
        terms( 0, -1,  1,  0,  0,       5,     0,     0,   0),
        terms(-2, -1,  0,  2,  1,      -5,     0,     3,   0),
        terms(-2,  0,  0,  0,  1,      -5,     0,     3,   0),
        terms( 0,  0,  2,  2,  1,      -5,     0,     3,   0),
        terms(-2,  0,  2,  0,  1,       4,     0,     0,   0),
        terms(-2,  1,  0,  2,  1,       4,     0,     0,   0),
        terms( 0,  0,  1, -2,  0,       4,     0,     0,   0),
        terms(-1,  0,  1,  0,  0,      -4,     0,     0,   0),
        terms(-2,  1,  0,  0,  0,      -4,     0,     0,   0),
        terms( 1,  0,  0,  0,  0,      -4,     0,     0,   0),
        terms( 0,  0,  1,  2,  0,       3,     0,     0,   0),
        terms( 0,  0, -2,  2,  2,      -3,     0,     0,   0),
        terms(-1, -1,  1,  0,  0,      -3,     0,     0,   0),
        terms( 0,  1,  1,  0,  0,      -3,     0,     0,   0),
        terms( 0, -1,  1,  2,  2,      -3,     0,     0,   0),
        terms( 2, -1, -1,  2,  2,      -3,     0,     0,   0),
        terms( 0,  0,  3,  2,  2,      -3,     0,     0,   0),
        terms( 2, -1,  0,  2,  2,      -3,     0,     0,   0),
    ];

    let t = time::julian_cent(JD);

    let M1 = angle::limit_to_360(
        134.96298 + t*(477198.867398 + t*(0.0086972 + t/56250.0))
    ).to_radians();
    let M  = angle::limit_to_360(
        357.52772 + t*(35999.05034   - t*(0.0001603 + t/300000.0))
    ).to_radians();
    let D  = angle::limit_to_360(
        297.85036 + t*(445267.11148  - t*(0.0019142 - t/189474.0))
    ).to_radians();
    let F  = angle::limit_to_360(
        93.27191  + t*(483202.017538 - t*(0.0036825 - t/327270.0))
    ).to_radians();
    let om = angle::limit_to_360(
        125.04452 - t*(1934.136261   - t*(0.0020708 + t/450000.0))
    ).to_radians();

    let mut nut_in_long = 0.0;
    let mut nut_in_oblq = 0.0;

    let div = 0.0001/3600.0;

    for x in terms_for_nutation.iter() {
        let arg =
            (x.0 as f64) * D   +
            (x.1 as f64) * M   +
            (x.2 as f64) * M1  +
            (x.3 as f64) * F   +
            (x.4 as f64) * om;

        nut_in_long += ((x.5 as f64) + t*(x.6 as f64)/10.0) * arg.sin() * div;
        nut_in_oblq += ((x.7 as f64) + t*(x.8 as f64)/10.0) * arg.cos() * div;
    }

    (nut_in_long.to_radians(), nut_in_oblq.to_radians())
}

/**
Computes nutation in equatorial coordinates

# Returns

`(nut_in_asc, nut_in_dec)`

* `nut_in_asc`: Nutation in right ascension *| in radians*
* `nut_in_dec`: Nutation in declination *| in radians*

# Arguments

* `eq_point`   : Equatorial point uncorrected for nutation *| in radians*
* `nut_in_long`: Nutation in longitude *| in radians*
* `nut_in_oblq`: Nutation in obliquity *| in radians*
* `tru_oblq`   : True obliquity of the ecliptic *| in radians*

The declination of `eq_point` should not be close to either of the
two of the celestial poles, as the values of nutation returned here
are only from first-order corrections.
**/
pub fn nutation_in_eq_coords(eq_point: &coords::EqPoint, nut_in_long: f64,
                             nut_in_oblq : f64, tru_oblq: f64) -> (f64, f64)
{
    let (asc, dec) = (eq_point.asc, eq_point.dec);

    let nut_asc = nut_in_long * (
        tru_oblq.cos()
      + tru_oblq.sin() * asc.sin() * dec.tan()
    ) - asc.cos() * dec.tan() * nut_in_oblq;

    let nut_dec =
        tru_oblq.sin() * asc.cos() * nut_in_long
      + asc.sin() * nut_in_oblq;

    (nut_asc, nut_dec)
}

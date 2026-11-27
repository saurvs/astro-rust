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

//! Near-parabolic orbits

use crate::angle;
use crate::consts;

/**
Computes the true anomaly and radius vector of a body in a near-parabolic
orbit at a given time

# Returns

`(true_anom, rad_vec)`

* `true_anom`: True anomaly of the body at time `t` *| in radians*
* `rad_vec`  : Radius vector of the body at time `t` *| in AU*

# Arguments

* `t`       : Time of interest, in Julian (Ephemeris) day
* `T`       : Time of passage through the perihelion, in Julian (Ephemeris) day
* `ecc`     : Eccentricity of the near-parabolic orbit
* `q`       : Perihelion distance *| in AU*
* `accuracy`: Desired accuracy for the results. For example, passing `0.000001`
              gives that much accuracy in radians and AU.
**/
pub fn true_anom_and_rad_vec<'a> (

    t        : f64,
    T        : f64,
    ecc      : f64,
    q        : f64,
    accuracy : f64

) -> Result<(f64, f64), &'a str>  {

    let days_frm_perih = t - T;

    if days_frm_perih == 0.0 {
        return Ok((0.0, q))
    }

    let k = consts::GAUSS_GRAV;
    let d1 = 1000.0;

    let q1 = k * ((1.0 + ecc)/q).sqrt() / (2.0*q);
    let g = (1.0 - ecc) / (1.0 + ecc);

    let q2 = q1 * days_frm_perih;
    let mut s = 2.0 / (3.0 * q2.abs());
    s = 2.0 / (2.0 * (s.atan() / 2.0).tan().cbrt().atan()).tan();

    if days_frm_perih < 0.0 { s = -s; }
    if ecc != 1.0 {
        let mut l = 0.0;
        loop {
            let s0 = s;
            let mut z = 1.0;
            let y = s * s;
            let mut g1 = -y * s;
            let mut q3 = q2 + g*s*y*2.0/3.0;

            loop {
                z += 1.0;
                g1 = -g1 * g * y;
                let z1 = (z - (z + 1.0)*g) / (2.0*z + 1.0);
                let f = z1 * g1;
                q3 += f;
                if z > 50.0 || f.abs() > d1 {
                    return Err("No convergence at orbit::near_parabolic::true_anom_and_rad_vec()");
                }
                if f.abs() <= accuracy { break; }
            }

            l += 1.0;
            if l > 50.0 {
                return Err("No convergence at orbit::near_parabolic::true_anom_and_rad_vec()");
            }

            loop {
                let s1 = s;
                s = (s*s*s*2.0/3.0 + q3)/(s*s + 1.0);
                if (s - s1).abs() <= accuracy { break; }
            }

            if (s - s0).abs() <= accuracy { break; }
        }
    }

    let v = angle::limit_to_two_PI(2.0 * s.atan());
    let r = q * (1.0 + ecc) / (1.0 + ecc*v.cos());

    Ok((v, r))

}

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

//! Parabolic orbits

use std;
use crate::orbit;

/**
Computes the true anomaly and radius vector of a body in a parabolic
orbit at a given time

# Returns

`(tru_anom, rad_vec)`

* `true_anom`: True anomaly of the body at time `t` *| in radians*
* `rad_vec`  : Radius vector of the body at time `t` *| in AU*

# Arguments

* `t`: Time of interest, in Julian (Ephemeris) day
* `T`: Time of passage through the perihelion, in Julian (Ephemeris) day
* `q`: Perihelion distance *| in AU*
**/
pub fn true_anom_and_rad_vec(t: f64, T: f64, q: f64) -> (f64, f64) {

    let W = 0.03649116245 * (t - T) / q.powf(1.5);
    let G = W / 2.0;
    let Y = (G + (G*G + 1.0).sqrt()).cbrt();
    let s = Y - 1.0/Y;
    let v = 2.0 * s.atan();
    let r = q * (1.0 + s*s);

    (v, r)

}

/**
Computes the time of passage of a body through a node of a parabolic
orbit and it's radius vector at that time

# Returns

`(time_of_pass, rad_vec)`

* `time_of_pass`: Time of passage through the node, in Julian
                  (Ephemeris) day
* `rad_vec`     : Radius vector of the body at the time of passage
                  *| in AU*

# Arguments

* `w`   : Argument of the perihelion *| in radians*
* `q`   : Perihelion distance *| in AU*
* `T`   : Time of passage in perihelion, in Julian (Ephemeris) day
* `node`: `Ascend` or `Descend` node
**/

#[inline]
pub fn passage_through_node (

    w    : f64,
    q    : f64,
    T    : f64,
    node : &orbit::Node

) -> (f64, f64) {

    match *node {
        orbit::Node::Ascend  => pass_through_node(-w, q, T),
        orbit::Node::Descend => pass_through_node(std::f64::consts::PI - w, q, T)
    }

}

fn pass_through_node(v: f64, q: f64, T: f64) -> (f64, f64) {

    let s = (v / 2.0).tan();
    let T_node = T + q.powf(1.5) * (s * (s*s + 3.0)) *  27.403895;
    let rad_vec = q * (1.0 + s*s);

    (T_node, rad_vec)

}

//! Parabolic orbits

use std::*;
use orbit;

/**
Returns the true anomaly and radius vector of a body in a
parabolic orbit

# Returns

`(tru_anom, rad_vec)`

* `tru_anom`: True anomaly of the body *| in radians*
* `rad_vec` : Radius vector of the body *| in AU*

# Arguments

* `t`: The time of interest, in Julian (Ephemeris) day
* `T`: Time of passage in perihelion, in Julian (Ephemeris) day
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
Returns the time of passage of a body through a node,
along with it's radius vector at that time

# Returns

`(time_of_pass, rad_vec)`

* `time_of_pass`: Time of passage through the node,
                  in Julian (Ephemeris) day
* `rad_vec`     : Radius vector of the body at the
                  time of passage *| in AU*

# Arguments
E
* `w`   : Argument of the perihelion *| in radians*
* `q`   : Perihelion distance *| in AU*
* `T`   : Time of passage in perihelion, in Julian
          (Ephemeris) day
* `node`: `Ascend` or `Descend` node
**/
pub fn passage_through_node(w: f64, q: f64, T: f64, node: &orbit::Node) -> (f64, f64) {
    match node {
        &orbit::Node::Ascend  => pass_through_node(-w, q, T),
        &orbit::Node::Descend => pass_through_node(f64::consts::PI - w, q, T)
    }
}
fn pass_through_node(v: f64, q: f64, T: f64) -> (f64, f64) {
    let s = (v/2.0).tan();

    (
        T + 27.403895*(s*(s*s + 3.0))*q.powf(1.5),
        q*(1.0 + s*s)
    )
}

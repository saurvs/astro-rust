/**
Returns the true anomaly and radius vector of a body in a parabolic orbit

# Returns

```(tru_anom, rad_vec)```

* ```tru_anom```: True anomaly of the body (*radians*)
* ```rad_vec```: Radius vector of the body from the Sun (*AU*)

# Arguments

* ```t```: Current time, in Julian (Ephemeris) day
* ```T```: Time of passage in perihelion, in Julian (Ephemeris) day
* ```q```: Perihelion distance (*AU*)
**/
pub fn TruAnomAndRadVec(t: f64, T: f64, q: f64) -> (f64, f64) {
    let W = 0.03649116245 * (t - T) / q.powf(1.5);
    let G = W / 2.0;
    let Y = (G + (G*G + 1.0).sqrt()).powf(1.0/3.0);
    let s = Y - 1.0/Y;
    let v = 2.0 * s.atan();
    let r = q * (1.0 + s*s);

    (v, r)
}

pub fn TimeOfPassThroughAscendNode(w: f64, q: f64, T: f64) -> (f64, f64) {
    time_of_passage_through_node(-1.0 * w, q, T)
}

pub fn TimeOfPassThroughDescendNode(w: f64, q: f64, T: f64) -> (f64, f64) {
    time_of_passage_through_node(180_f64.to_radians() * w, q, T)
}

fn time_of_passage_through_node(v: f64, q: f64, T: f64) -> (f64, f64) {
    let s = (v / 2.0).tan();
    let ss = s * s;
    (T + 27.403895*(ss*(s + 3.0)) * q.powf(1.5),
     q*(1.0 + ss))
}

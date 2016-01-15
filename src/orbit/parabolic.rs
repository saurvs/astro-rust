fn ABC_abc_terms(obl_eclp: f64, long_asc_node: f64, inc: f64, perih_arg: f64) -> (f64, f64, f64, f64, f64, f64) {
    let sin_obl_eclp = obl_eclp.sin();
    let cos_obl_eclp = obl_eclp.cos();
    let cos_long_asc_node = long_asc_node.cos();
    let sin_long_asc_node = long_asc_node.sin();
    let cos_inc = inc.cos();
    let sin_inc = inc.sin();

    let f = cos_long_asc_node;
    let g = sin_long_asc_node*cos_obl_eclp;
    let h = sin_long_asc_node*sin_obl_eclp;
    let p = -1.0*sin_long_asc_node*sin_inc;
    let q = cos_long_asc_node*cos_inc*cos_obl_eclp - sin_inc*sin_obl_eclp;
    let r = cos_long_asc_node*cos_inc*sin_obl_eclp + sin_inc*cos_obl_eclp;

    let A = f.atan2(p);
    let B = g.atan2(q);
    let C = h.atan2(r);
    let a = (f*f + p*p).sqrt();
    let b = (g*g + q*q).sqrt();
    let c = (h*h + r*r).sqrt();

    (A, B, C, a, b, c)
}

/**
Returns the **true anomaly** and **radius vector** (*AU*) of a body in a parabolic orbit

# Returns

```(true_anomaly, radius_vector)```

* ```true_anomaly```: True anomaly of the body
* ```radius_vector```: Radius vector of the body from the Sun (*AU*)

# Arguments

* ```t```: Time
* ```time_passg_perih```: Time of passage in perihelion
* ```perih_dist```: Perihelion distance (*AU*)
**/
pub fn TruAnomAndRadVec(obl_eclp: f64, long_asc_node: f64, inc: f64, perih_arg: f64, time_passg_perih: f64, t: f64, perih_dist: f64) -> (f64, f64) {
    let W = 0.03649116245 * (t - time_passg_perih) / perih_dist.powf(1.5);

    let G = W / 2.0;
    let Y = G + (G*G + 1.0).sqrt();
    let s = Y - 1.0/Y;
    let v = 2.0 * s.atan();
    let r = perih_dist * (1.0 + s*s);

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

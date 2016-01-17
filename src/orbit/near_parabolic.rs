use consts;

/**
Returns the **true anomaly** and **radius vector** of a body in a near-parabolic orbit

# Returns

```(true_anomaly, radius_vector)```

* ```true_anomaly```: True anomaly of the body
* ```radius_vector```: Radius vector of the body from the Sun (*AU*)

# Arguments

* ```ecc```: Eccentricity of the orbit
* ```t```: Time of perihelion passage in days
* ```perih_dist```: Perihelion distance (*AU*)
**/
pub fn TruAnomAndRadVec(ecc: f64, t: f64, perih_dist: f64) -> (f64, f64) {

    let Q = consts::GaussGrav()*(1.0 + ecc).sqrt()/(2.0 * perih_dist.powf(1.5));
    let y = (1.0 - ecc)/(1.0 + ecc);
    let s: f64 = 0.2;

    let v = 2.0 * s.atan();
    let r = perih_dist * (1.0 + ecc)/(1.0 + ecc*v.cos());

    (v, r)

}

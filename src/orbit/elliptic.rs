//! Elliptic orbits

use std::*;
use orbit;

/**
Returns the **true anomaly** of a body in an elliptic orbit

# Returns

* ```true_anom```: True anomaly of the body *| in radians*

# Arguments

* ```ecc_anom```: Eccentric anomaly of the body *| in radians*
* ```ecc```: Eccentricity of the orbit
**/
pub fn true_anom(ecc_anom: f64, ecc: f64) -> f64 {
    2.0 * ((1.0 + ecc).sqrt() * (ecc_anom/2.0).tan())
          .atan2((1.0 - ecc).sqrt())
}

/**
Returns the **radius vector** of a body in an elliptic orbit from
it's **eccentric anomaly**

# Returns

* ```rad_vec```: Radius vector of the body *| in AU*

# Arguments

* ```ecc_anom```: Eccentric anomaly of the body *| in radians*
* ```a```: Semimajor axis of the orbit *| in AU*
* ```ecc```: Eccentricity of the orbit
**/
pub fn rad_vec_frm_ecc_anom(ecc_anom: f64, a: f64, ecc: f64) -> f64 {
    a*(1.0 - ecc*ecc_anom.cos())
}

/**
Returns the **radius vector** of a body in an elliptic orbit from
it's **true anomaly**

# Returns

* ```rad_vec```: Radius vector of the body *| in AU*

# Arguments

* ```true_anom```: True anomaly of the body *| in radians*
* ```a```: Semimajor axis of the orbit *| in AU*
* ```ecc```: Rccentricity of the orbit
**/
pub fn rad_vec_frm_true_anom(true_anom: f64, a: f64, ecc: f64) -> f64 {
    a*(1.0 - ecc*ecc) / (1.0 + ecc*true_anom.cos())
}

/**
Returns the **eccentric anomaly** of a body in an elliptic orbit

# Returns

* ```ecc_anom```: Eccentric anomaly of the body *| in radians*

# Arguments

* ```mean_anom```: Mean anomaly of the body *| in radians*
* ```ecc```: Eccentricity of the orbit
* ```accuracy```: Desired accuracy for the eccentric anomaly. Eg: 0.000001 radians
**/
pub fn ecc_anom(mean_anom: f64, ecc: f64, accuracy: f64) -> f64 {
    let mut prev_E = 0.0;
    let mut current_E = mean_anom;
    while (current_E - prev_E).abs() > accuracy {
        prev_E = current_E;
        current_E = mean_anom + ecc*current_E.sin();
    }

    current_E
}

/**
Returns the **velocity** of a body in an elliptic orbit

# Returns

* ```velocity```: Instantaneous velocity of the body
                    *| in kilometers per second*

# Arguments

* ```r```: Radius vector of the body *| in AU*
* ```a```: Semimajor axis of orbit *| in AU*
**/
pub fn velocity(r: f64, a:f64) -> f64 {
    42.1219 * (1.0/r - 1.0/(2.0 * a)).sqrt()
}

/**
Returns the **velocity** of a body at **perihelion**
in an elliptic orbit

# Returns
* ```velocity```: Velocity of the body at perihelion
                    *| in kilometers per second*

# Arguments

* ```a```: Semimajor axis of orbit *| in AU*
* ```e```: Eccentricity of orbit
**/
pub fn perih_velocity(a:f64, e:f64) -> f64 {
    29.7847 * ((1.0 + e) / ((1.0 - e) * a)).sqrt()
}

/**
Returns the **velocity** of a body at **aphelion**
in an elliptic orbit

# Returns

* ```velocity```: Velocity of the body at aphelion
                    *| in kilometers per second*

# Arguments

* ```a```: Semimajor axis of orbit *| in AU*
* ```e```: Eccentricity of orbit
**/
pub fn aph_velocity(a:f64, e:f64) -> f64 {
    29.7847 * ((1.0 - e) / ((1.0 + e) * a)).sqrt()
}

/**
Returns the approximate **length** of an ellipse using the **Ramanujan**
method

# Returns

* ```approx_length```: An approximate value for the length of
                              the ellipse (same unit as that of ```a```
                              and ```b```), using a formula given by
                              [Ramanujan](https://en.wikipedia.org/wiki/Srinivasa_Ramanujan)
                              in 1914.

The **error** in ```(approximate_length)``` is:

* 0% for a = b
* 0.4155% for e = 1

# Arguments

* ```a```: Semimajor axis of the ellipse (same unit as that of ```b```)
* ```b```: Semiminor axis of the ellipse (same unit as that of ```a```)
* ```e```: Eccentricity of the ellipse
**/
pub fn length_ramanujan(a: f64, b: f64, e: f64) -> f64 {
    f64::consts::PI * (3.0*(a + b) - ((a + 3.0*b)*(3.0*a + b)).sqrt())
}

/**
Returns the approximate **length** of an **ellipse**

# Returns

* ```approx_length```: An approximate value for the length of
                              the ellipse (same unit as that of ```a```
                              and ```b```)

The **error** in ```(approximate_length)``` is:

* less than 0.001% if e < 0.88
* less than 0.01% if e < 0.95
* 1% if e = 0.9997
* 3% if e = 1

# Arguments

* ```a```: Semimajor axis of the ellipse (same unit as that of ```b```)
* ```b```: Semiminor axis of the ellipse (same unit as that of ```a```)
* ```e```: Eccentricity of the ellipse
**/
pub fn length(a: f64, b: f64, e: f64) -> f64 {
    let A = (a + b)/2.0;
    let G = (a * b).sqrt();
    let H = (2.0 * a * b)/(a + b);
    f64::consts::PI * (21.0*A - 2.0*G - 3.0*H) / 8.0
}

/**
Returns the **semimajor axis** of an elliptic orbit

# Arguments

* ```perih```: Perihelion of the orbit
* ```ecc```: Eccentricity of the orbit
**/
pub fn semimajor_axis(perih: f64, ecc: f64) -> f64 {
    perih / (1.0 - ecc)
}

/**
Returns the **mean motion** of an elliptic orbit

# Returns

* ```mean_motion```: Mean motion of the orbit
                     *| in radians per days of dynamical time*

# Arguments

* ```semimaj_ax```: Semimajor axis of the orbit
**/
pub fn mn_motion(semimaj_ax: f64) -> f64 {
    0.01720209895 / (semimaj_ax.powf(1.5))
}

/**
Returns the **time of passage** of a body through a **node**,
along with it's **radius vector** at that time

# Returns

```(time_of_pass, rad_vec)```

* ```time_of_pass```: Time of passage through the node, in Julian (Ephemeris) day
* ```rad_vec```: Radius vector of the body *| in AU*

# Arguments

* ```w```: Argument of the perihelion *| in radians*
* ```n```: Mean motion of the orbit (*radians per day*)
* ```a```: Semimajor axis of the orbit *| in AU*
* ```e```: Eccentricity of the orbit
* ```T```: Time of passage in perihelion, in Julian (Ephemeris) day
* ```node```: ```Ascend``` or ```Descend``` node
**/
pub fn passage_through_node(w: f64, n: f64, a: f64, e: f64, T: f64, node: &orbit::Node) -> (f64, f64) {
    match node {
        &orbit::Node::Ascend  => pass_through_node(-w, n, a, e, T),
        &orbit::Node::Descend => pass_through_node(f64::consts::PI - w, n, a, e, T)
    }
}

fn pass_through_node(v: f64, n: f64, a: f64, e: f64, T: f64)  -> (f64, f64) {
    let E = 2.0 * ((1.0 - e).sqrt()*(v/2.0).tan()).atan2((1.0 + e).sqrt());
    let M = E - e*E.sin();

    (T + M/n, a*(1.0 - e*E.cos()))
}

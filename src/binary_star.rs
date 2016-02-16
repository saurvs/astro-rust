//! Binary stars

use angle;
use std;

/**
Returns mean annual motion of companion star

* `P`: Period of revolution of binary star
       (*mean solar year*)
**/
pub fn mn_ann_motion_of_compan(P: f64) -> f64 {
    2.0 * std::f64::consts::PI / P
}

/**
Returns mean anomaly of companion star

# Arguments

* `n`: Mean annual motion of companion star
* `t`: Current time, given as year with
       decimals (eg: 1945.62)
* `T`: Time of periastron passage, given as
       a year with decimals (eg: 1945.62)
**/
pub fn mn_anom_of_compan(n: f64, t: f64, T: f64) -> f64 {
    n * (t - T)
}

/**
Returns radius vector of a binary star

# Arguments

* `a`       : Apparent semimajor axis
* `e`       : Eccentricity of true orbit
* `ecc_anom`: Eccentric anomaly of binary star
**/
pub fn rad_vec(a: f64, e: f64, ecc_anom: f64) -> f64 {
    a * (1.0 - e*ecc_anom.cos())
}

/**
Returns true anomaly of a binary star

# Arguments

* `e`       : Eccentricity of true orbit
* `ecc_anom`: Eccentric anomaly of binary star
**/
pub fn true_anom(e: f64, ecc_anom: f64) -> f64 {
    2.0 * (((1.0 + e)/(1.0 - e)).sqrt() * (ecc_anom / 2.0).tan()).atan()
}

/**
Returns apparent position angle of a binary star

# Arguments

* `asc_node_pos`: Position angle of ascending node
* `true_anom`   : True anomaly of binary star
* `w`           : Longitude of periastron
* `i`           : Inclination of true orbit to a
                  plane at right angles to the
                  line of sight *| in radians*
**/
pub fn apprnt_pos_angl(asc_node_pos: f64, true_anom: f64, w: f64, i: f64) -> f64 {
    let x = ((true_anom + w).sin() * i.cos()).atan2((true_anom + w).cos());
    angle::limit_to_360(x.to_degrees() + asc_node_pos.to_degrees()).to_radians()
}

/**
Returns angular separation of a binary star

# Arguments

* `rad_vec`  : Radius vector of a binary star
* `true_anom`: True anomaly of a binary star
* `w`        : Longitude of periastron
* `i`        : Inclination of true orbit to a
               plane at right angles to the
               line of sight *| in radians*
**/
pub fn angular_sepr(rad_vec: f64, true_anom: f64, w: f64, i: f64) -> f64 {
    rad_vec *
    (
      ((true_anom + w).sin() * i.cos()).powi(2) +
      (true_anom + w).cos().powi(2)
    ).sqrt()
}

/**
Returns eccentricity of an apparent orbit

# Arguments

* `e`: Eccentricity of the true orbit
* `i`: Inclination of true orbit
       to plane at right angles to line
       of sight *| in radians*
* `w`: Longitude of periastron
**/
pub fn ecc_of_apprnt_orb(e: f64, w: f64, i: f64) -> f64 {
    let a = (1.0 - (e * w.cos()).powi(2)) * i.cos().powi(2);
    let b = e.powi(2) * w.sin() * w.cos() * i.cos();
    let c = 1.0 - (e * w.sin()).powi(2);
    let d = ((a - c).powi(2) + 4.0*b.powi(2)).sqrt();

    ((2.0 * d) / (a + c + d)).sqrt()
}

use angle;

/**
Returns the **mean annual motion** of the **companion star**

* ```period_rev```: Period of revolution of the binary star
                    (*mean solar year*)
**/
pub fn MnAnnMotionOfCompan(period_rev: f64) -> f64 {
    360_f64.to_radians() / period_rev
}

/**
Returns the **mean anomaly** of the **companion star**

# Arguments

* ```n```: Mean annual motion of the companion star
* ```t```: Current time, given as a year with
           decimals (eg: 1945.62)
* ```T```: Time of periastron passage, given as a
           year with decimals (eg: 1945.62)
**/
pub fn MnAnomOfCompan(n: f64, t: f64, T: f64) -> f64 {
    n * (t - T)
}

/**
Returns the **radius vector** of a binary star

# Arguments

* ```a```: Semimajor axis of the binary star
* ```e```: Eccentricity of the true orbit
* ```ecc_anom```: Eccentric anomaly of the binary star
**/
pub fn RadVec(a: f64, e: f64, ecc_anom: f64) -> f64 {
    a * (1.0 - e*ecc_anom.cos())
}

/**
Returns the **true anomaly** of a binary star

# Arguments

* ```ecc_true_orb```: Eccentricity of the true orbit
* ```ecc_anom```: Eccentric anomaly of binary star
**/
pub fn TruAnom(ecc_true_orb: f64, ecc_anom: f64) -> f64 {
    2.0 * (((1.0 + ecc_true_orb)/(1.0 - ecc_true_orb)).sqrt() * (ecc_anom / 2.0).tan()).atan()
}

/**
Returns the **apparent position angle** of a binary star

# Arguments

* ```asc_node_pos```: Position angle of the ascending node
* ```true_anom```: True anomaly of the binary star
* ```periastron_long```: Longitude of the periastron
* ```i```: Inclination of the plane of the true orbit to
           the plane at right angles to the line of sight
**/
pub fn ApprntPosAngl(asc_node_pos: f64, true_anom: f64, periastron_long: f64, i: f64) -> f64 {
    let x = ((true_anom + periastron_long).sin() * i.cos()).atan2((true_anom + periastron_long).cos());
    angle::LimitTo360(x.to_degrees() + asc_node_pos.to_degrees()).to_radians()
}

/**
Returns the **angular separation** of a binary star

# Arguments

* ```rad_vec```: Radius vector of the binary star
* ```true_anom```: True anomaly of the binary star
* ```periastron_long```: Longitude of the periastron
* ```i```: Inclination of the plane of the true orbit
           to the plane at right angles to the line
           of sight
**/
pub fn AnglSepr(rad_vec: f64, true_anom: f64, periastron_long: f64, i: f64) -> f64 {
    rad_vec *
    (
      ((true_anom + periastron_long).sin() * i.cos()).powi(2) +
      (true_anom + periastron_long).cos().powi(2)
    ).sqrt()
}

/**
Returns the **eccentricity** of the apparent orbit

# Arguments

* ```ecc_true_orb```: Eccentricity of the true orbit
* ```i```: Inclination of the plane of the true orbit
           to the plane at right angles to the line
           of sight
* ```periastron_long```: Longitude of the periastron
**/
pub fn EccOfApprntOrb(ecc_true_orb: f64, periastron_long: f64, i: f64) -> f64 {
    let a = (1.0 - (ecc_true_orb * periastron_long.cos()).powi(2)) * i.cos().powi(2);
    let b = ecc_true_orb.powi(2) * periastron_long.sin() * periastron_long.cos() * i.cos();
    let c = 1.0 - (ecc_true_orb * periastron_long.sin()).powi(2);
    let d = ((a - c).powi(2) + 4.0*b.powi(2)).sqrt();

    ((2.0 * d) / (a + c + d)).sqrt()
}

use angle;

/*

    NOTE: All angles passed as arguments, and those returned,
          are assumed to be radians, even if the comments
          describe them with degrees.

*/

/**
Returns the mean annual motion of the companion star

* ```period_rev```: Period of revolution expressed in mean solar years
**/

pub fn comp_mean_annual_motion(period_rev: f64) -> f64 {
    360_f64.to_radians() / period_rev
}

/**
Returns the mean anomaly of the companion star

* ```mean_annual_motion```: Mean annual motion of the companion star
* ```periastron_pass```: Time of periastron passage, given as a
                           year with decimals (eg: 1945.62)
* ```time```: Current time, given as a year with decimals (eg: 1945.62)
**/

pub fn comp_mean_anomaly(mean_annual_motion: f64, periastron_pass: f64, time: f64) -> f64 {
    mean_annual_motion * (time - periastron_pass)
}

/**
Returns the radius vector of the binary star

* ```semimajor_axis```: Semimajor axis of the binary star
* ```ecc_true_orb```: Eccentricity of true orbit
* ```ecc_anomaly```: Eccentric anomaly of the binary star
**/

pub fn radius_vec(semimajor_axis: f64, ecc_true_orb: f64, ecc_anomaly: f64) -> f64 {
    semimajor_axis * (1.0 - ecc_true_orb * ecc_anomaly.cos())
}

/**
Returns the true anomaly of the binary star

* ```ecc_true_orb```: Eccentricity of true orbit
* ```ecc_anomaly```: Eccentric anomaly of the binary star
**/

pub fn true_anomaly(ecc_true_orb: f64, ecc_anomaly: f64) -> f64 {
    2.0 * (((1.0 + ecc_true_orb) / (1.0 - ecc_true_orb)).sqrt() * (ecc_anomaly / 2.0).tan()).atan()
}

/**
Returns the apparent position angle of the binary star

* ```asc_node_pos```: Position angle of the ascending node
* ```true_anomaly```: True anomaly of the binary star
* ```periastron_long```: Longitude of periastron
* ```inc```: Inclination of the plane of true orbit to the plane at
             right angles to the line of sight
**/

pub fn app_pos(asc_node_pos: f64, true_anomaly: f64, periastron_long: f64, inc: f64) -> f64 {
    let x = ((true_anomaly + periastron_long).sin() * inc.cos()).atan2((true_anomaly + periastron_long).cos());
    angle::limited_to_360(x.to_degrees() + asc_node_pos.to_degrees()).to_radians()
}

/**
Returns the angular separation of the binary star

* ```radius_vec```: Radius vector of the binary star
* ```true_anomaly```: True anomaly of the binary star
* ```periastron_long```: Longitude of periastron
* ```inc```: Inclination of the plane of true orbit to the plane at
             right angles to the line of sight
**/

pub fn angular_sep(radius_vec: f64, true_anomaly: f64, periastron_long: f64, inc: f64) -> f64 {
    radius_vec *
    (
      ((true_anomaly + periastron_long).sin() * inc.cos()).powi(2) +
      (true_anomaly + periastron_long).cos().powi(2)
    ).sqrt()
}

/**
Returns the eccentricity of the apparent orbit

* ```ecc_true_orb```: Eccentricity of true orbit
* ```true_anomaly```: True anomaly of the binary star
* ```inc```: Inclination of the plane of true orbit to the plane at
             right angles to the line of sight
* ```periastron_long```: Longitude of periastron

**/

pub fn ecc_app_orb(ecc_true_orb: f64, true_anomaly: f64, periastron_long: f64, inc: f64) -> f64 {
    let a = (1.0 - (ecc_true_orb * periastron_long.cos()).powi(2)) * inc.cos().powi(2);
    let b = ecc_true_orb.powi(2) * periastron_long.sin() * periastron_long.cos() * inc.cos();
    let c = 1.0 - (ecc_true_orb * periastron_long.sin()).powi(2);
    let d = ((a - c).powi(2) + 4.0 * b.powi(2)).sqrt();

    ((2.0 * d) / (a + c + d)).sqrt()
}

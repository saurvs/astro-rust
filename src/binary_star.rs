use angle;

/**
Computes **mean annual motion** of companion star

* ```period_rev```: Period of revolution expressed in mean solar years
**/
pub fn MeanAnnualMotionOfCompanion(period_rev: f64) -> f64 {
    360_f64.to_radians() / period_rev
}

/**
Computes **mean anomaly** of companion star

# Arguments

* ```mean_annual_motion```: Mean annual motion of companion star
* ```periastron_pass```: Time of periastron passage, given as a
                           year with decimals (eg: 1945.62)
* ```time```: Current time, given as a year with decimals (eg: 1945.62)
**/
pub fn MeanAnomalyOfCompanion(mean_annual_motion: f64, periastron_pass: f64, time: f64) -> f64 {
    mean_annual_motion * (time - periastron_pass)
}

/**
Computes **radius vector** of binary star

# Arguments

* ```semimajor_axis```: Semimajor axis of binary star
* ```ecc_true_orb```: Eccentricity of true orbit
* ```ecc_anomaly```: Eccentric anomaly of binary star
**/
pub fn RadiusVector(semimajor_axis: f64, ecc_true_orb: f64, ecc_anomaly: f64) -> f64 {
    semimajor_axis * (1.0 - ecc_true_orb * ecc_anomaly.cos())
}

/**
Computes **true anomaly** of binary star

# Arguments

* ```ecc_true_orb```: Eccentricity of true orbit
* ```ecc_anomaly```: Eccentric anomaly of binary star
**/
pub fn TrueAnomaly(ecc_true_orb: f64, ecc_anomaly: f64) -> f64 {
    2.0 * (((1.0 + ecc_true_orb) / (1.0 - ecc_true_orb)).sqrt() * (ecc_anomaly / 2.0).tan()).atan()
}

/**
Computes **apparent position angle** of binary star

# Arguments

* ```asc_node_pos```: Position angle of ascending node
* ```true_anomaly```: True anomaly of binary star
* ```periastron_long```: Longitude of periastron
* ```inc```: Inclination of the plane of true orbit to the plane at
             right angles to the line of sight
**/
pub fn ApparentPositionAngle(asc_node_pos: f64, true_anomaly: f64, periastron_long: f64, inc: f64) -> f64 {
    let x = ((true_anomaly + periastron_long).sin() * inc.cos()).atan2((true_anomaly + periastron_long).cos());
    angle::LimitedTo360(x.to_degrees() + asc_node_pos.to_degrees()).to_radians()
}

/**
Computes **angular separation** of binary star

# Arguments

* ```radius_vec```: Radius vector of binary star
* ```true_anomaly```: True anomaly of binary star
* ```periastron_long```: Longitude of periastron
* ```inc```: Inclination of the plane of true orbit to the plane at
             right angles to the line of sight
**/
pub fn AngularSep(radius_vec: f64, true_anomaly: f64, periastron_long: f64, inc: f64) -> f64 {
    radius_vec *
    (
      ((true_anomaly + periastron_long).sin() * inc.cos()).powi(2) +
      (true_anomaly + periastron_long).cos().powi(2)
    ).sqrt()
}

/**
Computes **eccentricity** of the **apparent orbit**

# Arguments

* ```ecc_true_orb```: Eccentricity of true orbit
* ```true_anomaly```: True anomaly of binary star
* ```inc```: Inclination of the plane of true orbit to the plane at
             right angles to the line of sight
* ```periastron_long```: Longitude of periastron
**/
pub fn EccentricityOfApparentOrbit(ecc_true_orb: f64, true_anomaly: f64, periastron_long: f64, inc: f64) -> f64 {
    let a = (1.0 - (ecc_true_orb * periastron_long.cos()).powi(2)) * inc.cos().powi(2);
    let b = ecc_true_orb.powi(2) * periastron_long.sin() * periastron_long.cos() * inc.cos();
    let c = 1.0 - (ecc_true_orb * periastron_long.sin()).powi(2);
    let d = ((a - c).powi(2) + 4.0 * b.powi(2)).sqrt();

    ((2.0 * d) / (a + c + d)).sqrt()
}

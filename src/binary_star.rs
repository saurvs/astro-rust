use angle;

/*

    NOTE: All angles passed as arguments, and those returned,
          are assumed to be radians, even though the comments
          describe them with degrees.

*/

/*

    Returns the mean annual motion of the companion star
    -----------------------------------------------------------------
        period_rev: Period of revolution expressed in mean solar years

*/

pub fn comp_mean_annual_motion(period_rev: f64) -> f64 {
    360_f64.to_radians() / period_rev
}

/*

    Returns the mean anomaly of the companion star
    -----------------------------------------------------------------
        mean_annual_motion: Mean annual motion of the companion star
           periastron_pass: Time of the periastron passage, given as
                            a year with decimals (eg: 1945.62)
                      time: current time given as a year with decimals

*/

pub fn comp_mean_anomaly(mean_annual_motion: f64, periastron_pass: f64, time: f64) -> f64 {
    mean_annual_motion * (time - periastron_pass)
}

/*

    Returns the radius vector
    -----------------------------------------------------------------
                   semimajor_axis: The semimajor axis
                     ecc_true_orb: The eccentricity of the true orbit
                      ecc_anomaly: The eccentric anomaly

*/

pub fn radius_vec(semimajor_axis: f64, ecc_true_orb: f64, ecc_anomaly: f64) -> f64 {
    semimajor_axis * (1.0 - ecc_true_orb * ecc_anomaly.cos())
}

/*

    Returns the true anomaly
    -----------------------------------------------------------------
                     ecc_true_orb: The eccentricity of the true orbit
                      ecc_anomaly: The eccentric anomaly

*/

pub fn true_anomaly(ecc_true_orb: f64, ecc_anomaly: f64) -> f64 {
    2.0 * (((1.0 + ecc_true_orb) / (1.0 - ecc_true_orb)).sqrt() * (ecc_anomaly / 2.0).tan()).atan()
}

/*

    Returns the apparent position angle
    -----------------------------------------------------------------
       asc_node_pos: The position angle of the ascending node
       true_anomaly: The true anomaly
    periastron_long: The longitude of the periastron
                inc: The inclination of the plane of the true orbit to
                     the plane at right angles to the line of sight

*/

pub fn app_pos(asc_node_pos: f64, true_anomaly: f64, periastron_long: f64, inc: f64) -> f64 {
    let x = ((true_anomaly + periastron_long).sin() * inc.cos()).atan2((true_anomaly + periastron_long).cos());
    angle::limited_to_360(x.to_degrees() + asc_node_pos.to_degrees()).to_radians()
}

/*

    Returns the angular separation
    -----------------------------------------------------------------
         radius_vec: The radius vector
       true_anomaly: The true anomaly
    periastron_long: The longitude of the periastron
                inc: the inclination of the plane of the true orbit to
                     the plane at right angles to the line of sight

*/

pub fn angular_sep(radius_vec: f64, true_anomaly: f64, periastron_long: f64, inc: f64) -> f64 {
    radius_vec *
    (
      ((true_anomaly + periastron_long).sin() * inc.cos()).powi(2) +
      (true_anomaly + periastron_long).cos().powi(2)
    ).sqrt()
}

/*

    Returns the eccentricity of the apparent orbit
    -----------------------------------------------------------------
       ecc_true_orb: The eccentricity of the true orbit
                inc: The inclination of the plane of the true orbit
                     to the plane at right angles to the line of sight
     periastron_long: The longitude of the periastron

*/

pub fn ecc_app_orb(ecc_true_orb: f64, true_anomaly: f64, periastron_long: f64, inc: f64) -> f64 {
    let A = (1.0 - (ecc_true_orb * periastron_long.cos()).powi(2)) * inc.cos().powi(2);
    let B = ecc_true_orb.powi(2) * periastron_long.sin() * periastron_long.cos() * inc.cos();
    let C = 1.0 - (ecc_true_orb * periastron_long.sin()).powi(2);
    let D = ((A - C).powi(2) + 4.0 * B.powi(2)).sqrt();

    ((2.0 * D) / (A + C + D)).sqrt()
}

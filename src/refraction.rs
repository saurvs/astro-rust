use angle;

/*

    NOTE: All angles passed as arguments, and those returned,
          are assumed to be radians, even though the comments
          describe them with degrees.

    The refraction correction angle R needs to be found to account
    for the effects of atmospheric refraction.
    We will only use the horizontal coordinate system for convenience
    (altitude and azimuth).

    R, at altitudes very near 0 or 180 degrees, that is, very near to
    the horizon, fluctuates unpredictably by 0.3 degrees or more,
    and so no sensible estimate can be made for it

*/

/*

    The following two functions may be used to find R when the
    altitude is larger that 15 degrees.

    The first function refrc_from_app_alt() is used to find R given
    an apparent altitude, in which case the returned R needs to be
    subtracted from the apparent altitude.

    The second function refrc_from_true_alt() is used to find R given
    a true altitude, in which case the returned R needs to be added
    to the true altitude.

*/

pub fn refrc_from_app_alt(app_alt: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 58.294).to_radians() * (90f64.to_radians() - app_alt).tan() -
    angle::pure_degrees(0.0, 0.0, 0.0668).to_radians() * (90f64.to_radians() - app_alt).tan().powi(3)
}

pub fn refrc_from_true_alt(true_alt: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 58.276).to_radians() * (90f64.to_radians() - true_alt).tan() -
    angle::pure_degrees(0.0, 0.0, 0.0824).to_radians() * (90f64.to_radians() - true_alt).tan().powi(3)
}

use angle;

fn close_to_90_degrees(a: f64) -> bool {
    if (90.0 - a).abs() < 0.0001 { true }
    else                         { false }
}

fn pt_correction(prss: f64, temp: f64) -> f64 {
    (prss/1010.0) * (283.0 / (273.0+temp))
}

/*

    NOTE: All angles passed as arguments, and those returned,
          are assumed to be radians, even though the comments
          describe them with degrees.

*/

/*

    The refraction correction angle R needs to be found to account
    for the effects of atmospheric refraction.
    We will only use the horizontal coordinate system for convenience
    (that is, altitude and azimuth).

    R, at altitudes very near to 0 or to 180 degrees (that is, very
    near to the horizon), fluctuates unpredictably by 0.3 degrees or
    more, and so no sensible estimate can be made for R there.

*/

/*

    The following two functions can be used to find R when the
    altitude is larger that 15 degrees.

    The first function approx_refrc_from_app_alt() is used to find R,
    given an apparent altitude (app_alt), in which case the returned
    R needs to be subtracted from the apparent altitude to obtain the
    true altitude.

    The second function approx_refrc_from_true_alt() is used to find R,
    given a true altitude (true_alt), in which case the returned R needs
    to be added to the true altitude to obtain the apparent altitude.

*/

pub fn approx_refraction_from_app_alt(app_alt: f64) -> f64 {
    angle::PureDegrees(0.0, 0.0, 58.294).to_radians() * (90_f64.to_radians()-app_alt).tan() -
    angle::PureDegrees(0.0, 0.0, 0.0668).to_radians() * (90_f64.to_radians()-app_alt).tan().powi(3)
}

pub fn approx_refraction_from_true_alt(true_alt: f64) -> f64 {
    angle::PureDegrees(0.0, 0.0, 58.276).to_radians() * (90_f64.to_radians()-true_alt).tan() -
    angle::PureDegrees(0.0, 0.0, 0.0824).to_radians() * (90_f64.to_radians()-true_alt).tan().powi(3)
}

/*

    The next two functions can be used to find R for all values of
    altitude from 0 to 90 degrees.

    The first function refrc_from_app_alt() is used to find R, given
    an apparent altitude (app_alt), pressure in millibars (prss) and
    temperature in Celsius (temp). In this case the returned R
    needs to be subtracted from the apparent altitude to obtain the
    true altitude.

    The second function refrc_from_true_alt() is used to find R,
    given a true altitude (app_alt), pressure in millibars (prss)
    and temperature in Celsius (temp). In this case the returned R
    needs to be added from the true altitude to obtain the apparent
    altitude.

*/

pub fn refraction_from_app_alt(mut app_alt: f64, prss: f64, temp: f64) -> f64 {
    app_alt = app_alt.to_degrees();
    if close_to_90_degrees(app_alt) { return 0.0; }
    (1.0 / (60.0 * (app_alt + (7.31 / (app_alt + 4.4))).to_radians().tan())) * pt_correction(prss, temp)
}

pub fn refraction_from_true_alt(mut true_alt: f64, prss: f64, temp: f64) -> f64 {
    true_alt = true_alt.to_degrees();
    if close_to_90_degrees(true_alt) { return 0.0; }
    (1.02 / (60.0 * ((true_alt + (10.3 / (true_alt + 5.11)))).to_radians().tan())) *
    pt_correction(prss, temp)
}

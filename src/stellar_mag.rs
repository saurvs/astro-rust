use angle;

/*

    NOTE: All angles passed as arguments, and those returned,
          are assumed to be radians, even though the comments
          describe them with degrees.

*/

/*

    Returns the combined magnitude of two stars
    -----------------------------------------------------------------
        m1: Magnitude of star 1
        m2: Magnitude of star 2

*/

pub fn comb_mag(m1: f64, m2: f64) -> f64 {
    m2 - 2.5 * (bright_ratio(m1, m2) + 1.0)
}

/*

    Returns the combined magnitude of more than two stars
    -----------------------------------------------------------------
        m: An array of magnitudes of stars

*/

pub fn comb_mag_from_many(m: &[f64]) -> f64 {
    let mut sum = 0.0;
    for i in m.iter() {
        sum += 10_f64.powf(-0.4 * i);
    }
    -2.5 * sum.log(10.0)
}

/*

    Returns the brightness ratio of two stars
    -----------------------------------------------------------------
        m1: Magnitude of star 1
        m2: Magnitude of star 2

*/

pub fn bright_ratio(m1: f64, m2: f64) -> f64 {
    10.0_f64.powf(0.4 * (m2 - m1))
}

/*

    Returns the magnitude difference of two stars
    -----------------------------------------------------------------
        br: The brightness ratio of two stars

*/

pub fn mag_diff(br: f64) -> f64 {
    2.5 * br.log(10.0)
}

/*

    Returns the absolute magnitude of a star
    -----------------------------------------------------------------
        par: The parallax of a star
          m: The apparent magnitude of a star

*/

pub fn abs_mag_from_app_mag(mut par: f64, m: f64) -> f64 {
    par = par.to_degrees() * 3600.0;
    m + 5.0 + 5.0 * par.log(10.0)
}

/*

    Returns the absolute magnitude of a star
    -----------------------------------------------------------------
        d: The distance of a star in parsecs
        m: The apparent magnitude of a star

*/

pub fn abs_mag_from_dist(d: f64, m: f64) -> f64 {
    m + 5.0 - 5.0 * d.log(10.0)
}

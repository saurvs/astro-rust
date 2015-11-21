/*

    combined_mag(m1, m2)
    -----------------------------------------------------------------
    Returns the combined magnitude of two stars

        m1: Magnitude of star 1
        m2: Magnitude of star 2

*/

pub fn combined_mag(m1: f64, m2: f64) -> f64 {
    m2 - 2.5 * (brightness_ratio(m1, m2) + 1.0)
}

/*

    combined_mag_from_many(m)
    -----------------------------------------------------------------
    Returns the combined magnitude of more than two stars

        m: An array of magnitudes of stars

*/

pub fn combined_mag_from_many(m: &[f64]) -> f64 {
    let mut sum = 0.0;
    for i in m.iter() {
        sum += 10_f64.powf(-0.4 * i);
    }
    -2.5 * sum.log(10.0)
}

/*

    brightness_ratio(m1, m2)
    -----------------------------------------------------------------
    Returns the brightness ratio of two stars

        m1: Magnitude of star 1
        m2: Magnitude of star 2

*/

pub fn brightness_ratio(m1: f64, m2: f64) -> f64 {
    10.0_f64.powf(0.4 * (m2 - m1))
}

/*

    mag_diff(br)
    -----------------------------------------------------------------
    Returns the magnitude difference of two stars

        br: The brightness ratio of two stars

*/

pub fn mag_diff(br: f64) -> f64 {
    2.5 * br.log(10.0)
}

/*

    absolute_mag_from_apparent_mag(par, am)
    -----------------------------------------------------------------
    Returns the absolute magnitude of a star from its parallax

        par: The parallax of a star
         am: The apparent magnitude of a star

*/

pub fn absolute_mag_from_parallax(mut par: f64, am: f64) -> f64 {
    par = par.to_degrees() * 3600.0;
    am + 5.0 + 5.0 * par.log(10.0)
}

/*

    absolute_mag_from_dist(d, am)
    -----------------------------------------------------------------
    Returns the absolute magnitude of a star from its distance

         d: The distance of a star in parsecs
        am: The apparent magnitude of a star

*/

pub fn absolute_mag_from_dist(d: f64, am: f64) -> f64 {
    am + 5.0 - 5.0 * d.log(10.0)
}

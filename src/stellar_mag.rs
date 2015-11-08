use angle;

/*

    NOTE: All angles passed as arguments, and those returned,
          are assumed to be radians, even though the comments
          describe them with degrees.

*/

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

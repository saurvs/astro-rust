use angle;

fn unit_semidia_saturn_eq() -> f64 {
    angle::pure_degrees(0.0, 0.0, 82.73)
}

fn unit_semidia_saturn_pol() -> f64 {
    angle::pure_degrees(0.0, 0.0, 73.82)
}

fn unit_semidia_jupiter_pol() -> f64 {
    angle::pure_degrees(0.0, 0.0, 92.06)
}

/*

    The following eleven functions return the angular semidiameters
    at of various bodies in the solar system.
    -----------------------------------------------------------------
        del: The body's distance to Earth in AU

*/

pub fn semidia_sun(del: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 959.63) / del
}

pub fn semidia_mercury(del: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 3.36) / del
}

pub fn semidia_venus(del: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 8.41) / del
}

pub fn semidia_mars(del: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 4.68) / del
}

pub fn semidia_jupiter_eq(del: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 98.44) / del
}

pub fn semidia_jupiter_pol(del: f64) -> f64 {
    unit_semidia_jupiter_pol() / del
}

pub fn semidia_saturn_eq(del: f64) -> f64 {
    unit_semidia_saturn_eq() / del
}

pub fn semidia_saturn_pol(del: f64) -> f64 {
    unit_semidia_saturn_pol() / del
}

pub fn semidia_uranus(del: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 35.02) / del
}

pub fn semidia_neptune(del: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 33.5) / del
}

pub fn semidia_pluto(del: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 2.07) / del
}

/*

    Returns the apparent polar semidiameter of Saturn
    -----------------------------------------------------------------
            del: The body's distance to Earth in AU
        sat_lat: The Saturnicentric latitude of Earth

*/

pub fn app_pol_semidia_saturn(del: f64, sat_lat: f64) -> f64 {
    let a = unit_semidia_saturn_eq();
    let b = unit_semidia_saturn_pol();
    let k = 1.0 - (b / a).powi(2);
    (a / del) * (1.0 - k * sat_lat.cos().powi(2)).sqrt()
}

/*

    Returns the apparent polar semidiameter of Jupiter
    -----------------------------------------------------------------
            del: The body's distance to Earth in AU

*/

pub fn app_pol_semidia_jupiter(del: f64) -> f64 {
    unit_semidia_jupiter_pol() / del
}

/*

    Returns the diameter of an asteroid in kilometers
    -----------------------------------------------------------------
           abs_mag: The absolute magnitude of an asteroid
            albedo: The albedo or reflective power

*/

pub fn astroid_diameter(abs_mag: f64, albedo: f64) -> f64 {
    10.0_f64.powf(3.12 - (abs_mag / 5.0) - (0.217147 * albedo.log(10.0)))
}

/*

    Returns the apparent diameter of an asteroid in kilometers
    -----------------------------------------------------------------
        diameter: The true diameter of an asteroid in kilometers
             del: The asteroid's distance to Earth in AU

*/

pub fn astroid_app_diameter(diameter: f64, del: f64) -> f64 {
    0.0013788 * (diameter / del)
}

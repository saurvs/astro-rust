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

    semidia_body(distance_to_earth) -> (angular_semidiameter)

    distance_to_earth should be in AU
    -----------------------------------------------------------------
    The following eleven functions return the angular semidiameters
    of various bodies in the solar system.

*/

pub fn semidia_sun(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 959.63) / distance_to_earth
}

pub fn semidia_mercury(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 3.36) / distance_to_earth
}

pub fn semidia_venus(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 8.41) / distance_to_earth
}

pub fn semidia_mars(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 4.68) / distance_to_earth
}

pub fn semidia_jupiter_eq(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 98.44) / distance_to_earth
}

pub fn semidia_jupiter_pol(distance_to_earth: f64) -> f64 {
    unit_semidia_jupiter_pol() / distance_to_earth
}

pub fn semidia_saturn_eq(distance_to_earth: f64) -> f64 {
    unit_semidia_saturn_eq() / distance_to_earth
}

pub fn semidia_saturn_pol(distance_to_earth: f64) -> f64 {
    unit_semidia_saturn_pol() / distance_to_earth
}

pub fn semidia_uranus(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 35.02) / distance_to_earth
}

pub fn semidia_neptune(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 33.5) / distance_to_earth
}

pub fn semidia_pluto(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 2.07) / distance_to_earth
}

/*

    app_polar_semidia_saturn(distance_to_earth,
                            saturnicentric_latitude_of_Earth)
                            -> (apparent_polar_semidiameter_of_saturn)

    distance_to_earth should be in AU
    -----------------------------------------------------------------
    Returns the apparent polar semidiameter of Saturn

*/

pub fn app_polar_semidia_saturn(distance_to_earth: f64, sat_lat: f64) -> f64 {
    let a = unit_semidia_saturn_eq();
    let b = unit_semidia_saturn_pol();
    let k = 1.0 - (b / a).powi(2);
    (a / distance_to_earth) * (1.0 - k * sat_lat.cos().powi(2)).sqrt()
}

/*

    app_polar_semidia_jupiter(distance_to_earth)
                           -> (apparent_polar_semidiameter_of_jupiter)

    distance_to_earth should be in AU
    -----------------------------------------------------------------
    Returns the apparent polar semidiameter of Jupiter

*/

pub fn app_polar_semidia_jupiter(distance_to_earth: f64) -> f64 {
    unit_semidia_jupiter_pol() / distance_to_earth
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

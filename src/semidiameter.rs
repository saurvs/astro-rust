use angle;

fn unit_semidia_saturn_eq() -> f64 {
    angle::pure_degrees(0.0, 0.0, 82.73)
}

fn unit_semidia_saturn_pol() -> f64 {
    angle::pure_degrees(0.0, 0.0, 73.82)
}

/*

    NOTE: All angles passed as arguments, and those returned,
          are assumed to be radians, even if the comments
          describe them with degrees.

*/

/*

    Semidiameter: The apparent equatorial or polar radius of a
    celestial body when viewed as a disc from Earth, expressed
    as an angle.

    For Saturn and Jupiter, the both apparent equatorial and
    polar radii are available.

    For the Sun, Pluto, and the rest of the planets, only the
    apparent equatorial radius is available.

*/

/*

    semidia_body(distance_to_earth) -> (angular_semidiameter)
    -----------------------------------------------------------------
    The following eleven functions return the angular semidiameters
    of various bodies in the solar system.

    distance_to_earth: The distance of the body from the Earth in AU

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

pub fn equatorial_semidia_jupiter(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 98.44) / distance_to_earth
}

pub fn polar_semidia_jupiter(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 92.06) / distance_to_earth
}

pub fn equatorial_semidia_saturn(distance_to_earth: f64) -> f64 {
    unit_semidia_saturn_eq() / distance_to_earth
}

pub fn polar_semidia_saturn(distance_to_earth: f64) -> f64 {
    let a = unit_semidia_saturn_eq();
    let b = unit_semidia_saturn_pol();
    let k = 1.0 - (b / a).powi(2);
    (a / distance_to_earth) * (1.0 - k * sat_lat.cos().powi(2)).sqrt()
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

    Returns the diameter of an asteroid in kilometers.
    -----------------------------------------------------------------
    astroid_diameter(absolute_magnitude_of_asteroid, albedo)
                                            -> (astroid_diameter)

    abs_mag: The absolute magnitude of the asteroid
    albedo: The albedo or reflective power of the asteroid

*/

pub fn astroid_diameter(abs_mag: f64, albedo: f64) -> f64 {
    10.0_f64.powf(3.12 - (abs_mag / 5.0) - (0.217147 * albedo.log(10.0)))
}

/*

    Returns the apparent diameter of an asteroid in kilometers.
    -----------------------------------------------------------------
    app_astroid_diameter(true_diameter_of_asteroid, distance_to_earth)
                                    -> (apparent_diameter_of_asteroid)

    true_diameter: The true diameter of the asteroid in kilometers
    distance_to_earth: The asteroid's distance to Earth in AU

*/

pub fn app_astroid_diameter(true_diameter: f64, distance_to_earth: f64) -> f64 {
    0.0013788 * (true_diameter / distance_to_earth)
}

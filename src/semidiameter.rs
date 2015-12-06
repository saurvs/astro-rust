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

    For Saturn and Jupiter, both the apparent equatorial and
    polar radii are available.

    For the Sun, Pluto, and the rest of the planets, only the
    apparent equatorial radius is available.

*/

/**
Returns the Sun's equatorial angular semidiameter

# Arguments

* ```distance_to_earth```: The Sun's distance from the Earth (in AU)
**/
pub fn semidia_sun(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 959.63) / distance_to_earth
}

/**
Returns Mercury's equatorial angular semidiameter

# Arguments

* ```distance_to_earth```: Mercury's distance from the Earth (in AU)
**/
pub fn semidia_mercury(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 3.36) / distance_to_earth
}

/**
Returns Venus's equatorial angular semidiameter

# Arguments

* ```distance_to_earth```: Venus's  from  Earth (in AU)
**/
pub fn semidia_venus(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 8.41) / distance_to_earth
}

/**
Returns Mars's equatorial angular semidiameter

# Arguments

* ```distance_to_earth```: Mars's distance from Earth (in AU)
**/
pub fn semidia_mars(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 4.68) / distance_to_earth
}

/**
Returns Jupiter's equatorial angular semidiameter

# Arguments

* ```distance_to_earth```: Jupiter's distance from Earth (in AU)
**/
pub fn equatorial_semidia_jupiter(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 98.44) / distance_to_earth
}

/**
Returns Jupiter's polar angular semidiameter

# Arguments

* ```distance_to_earth```: Jupiter's distance from Earth (in AU)
**/
pub fn polar_semidia_jupiter(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 92.06) / distance_to_earth
}

/**
Returns Saturn's equatorial angular semidiameter

# Arguments

* ```distance_to_earth```: Saturn's distance from Earth (in AU)
**/
pub fn equatorial_semidia_saturn(distance_to_earth: f64) -> f64 {
    unit_semidia_saturn_eq() / distance_to_earth
}

/**
Returns Saturn's polar angular semidiameter

# Arguments

* ```distance_to_earth```: Saturn's distance from Earth (in AU)
* ```earth_lat```: Earth's Saturnicentric latitude (in radians)
**/
pub fn polar_semidia_saturn(distance_to_earth: f64, earth_lat: f64) -> f64 {
    let a = unit_semidia_saturn_eq();
    let b = unit_semidia_saturn_pol();
    let k = 1.0 - (b / a).powi(2);
    (a / distance_to_earth) * (1.0 - k * earth_lat.cos().powi(2)).sqrt()
}

/**
Returns Uranus's equatorial angular semidiameter

# Arguments

* ```distance_to_earth```: Uranus's distance from Earth (in AU)
**/
pub fn semidia_uranus(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 35.02) / distance_to_earth
}

/**
Returns Neptune's equatorial angular semidiameter

# Arguments

* ```distance_to_earth```: Neptune's distance from Earth (in AU)
**/
pub fn semidia_neptune(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 33.5) / distance_to_earth
}

/**
Returns Pluto's equatorial angular semidiameter

# Arguments

* ```distance_to_earth```: Pluto's distance from Earth (in AU)
**/
pub fn semidia_pluto(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 2.07) / distance_to_earth
}

/**
Returns the diameter of an asteroid (in kilometers)

# Arguments

* ```abs_mag```: Absolute magnitude of asteroid
* ```albedo```: Reflective power of asteroid
**/
pub fn astroid_diameter(abs_mag: f64, albedo: f64) -> f64 {
    10.0_f64.powf(3.12 - (abs_mag / 5.0) - (0.217147 * albedo.log(10.0)))
}

/**
Returns the apparent diameter of an asteroid (in kilometers)

# Arguments

* ```true_diameter```: True diameter of asteroid (in kilometers)
* ```distance_to_earth```: Asteroid's distance to Earth (in AU)
**/
pub fn app_astroid_diameter(true_diameter: f64, distance_to_earth: f64) -> f64 {
    0.0013788 * (true_diameter / distance_to_earth)
}

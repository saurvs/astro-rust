//! Diameters of asteroids

/**
Returns the **diameter** of an asteroid

# Returns

* ```diameter```: Diameter of the asteroid *| in meters*

# Arguments

* ```abs_mag```: Absolute magnitude of the asteroid
* ```albedo```: Reflective power of the asteroid
**/
pub fn Diameter(abs_mag: f64, albedo: f64) -> f64 {
    1000.0 * 10_f64.powf(3.12 - abs_mag/5.0 - 0.217147*albedo.log10())
}

/**
Returns the **apparent diameter** of an asteroid

# Returns

* ```apprnt_diameter```: Apparent diameter of the asteroid *| in meters*

# Arguments

* ```true_diameter```: True diameter of the asteroid *| in kilometers*
* ```asteroid_earth_dist```: Asteroid-Earth distance *| in AU*
**/
pub fn ApprntDiameter(true_diameter: f64, asteroid_earth_dist: f64) -> f64 {
    1.3788 * (true_diameter / asteroid_earth_dist)
}

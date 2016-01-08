/**
Returns the **diameter** of an asteroid

# Returns

* ```diameter```: Diameter of the asteroid *(meters)*

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

* ```apparent_diameter```: Apparent diameter of the asteroid *(meters)*

# Arguments

* ```true_diameter```: True diameter of the asteroid *(kilometers)*
* ```distance_to_earth```: Asteroid's distance to Earth *(AU)*
**/
pub fn ApparentDiameter(true_diameter: f64, distance_to_earth: f64) -> f64 {
    1.3788 * (true_diameter / distance_to_earth)
}

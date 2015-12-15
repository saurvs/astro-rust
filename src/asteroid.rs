/**
Computes **diameter** of asteroid *(meters)*

# Arguments

* ```abs_mag```: Absolute magnitude of asteroid
* ```albedo```: Reflective power of asteroid
**/
pub fn diameter(abs_mag: f64, albedo: f64) -> f64 {
    10.0_f64.powf(3.12 - (abs_mag / 5.0) - (0.217147 * albedo.log(10.0)))  / 1000.0
}

/**
Computes **apparent diameter** of asteroid *(meters)*

# Arguments

* ```true_diameter```: True diameter of asteroid (in kilometers)
* ```distance_to_earth```: Asteroid's distance to Earth (in AU)
**/
pub fn apparent_diameter(true_diameter: f64, distance_to_earth: f64) -> f64 {
    0.0013788 * (true_diameter / distance_to_earth) / 1000.0
}

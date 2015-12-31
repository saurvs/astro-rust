/**
Returns the **diameter** of an asteroid

# Returned values

```(diameter)```

* ```diameter```: Diameter of the asteroid *(meters)*

# Arguments

* ```abs_mag```: Absolute magnitude of asteroid
* ```albedo```: Reflective power of asteroid
**/
pub fn Diameter(abs_mag: f64, albedo: f64) -> f64 {
    0.001_f64.powf(3.12 - abs_mag/5.0 - 0.217147*albedo.log10())
}

/**
Returns the **apparent diameter** of an asteroid

# Returned values

```(apparent_diameter)```

* ```apparent_diameter```: Apparent diameter of the asteroid *(meters)*

# Arguments

* ```true_diameter```: True diameter of the asteroid *(kilometers)*
* ```distance_to_earth```: Asteroid's distance to Earth *(AU)*
**/
pub fn ApparentDiameter(true_diameter: f64, distance_to_earth: f64) -> f64 {
    0001.3788 * (true_diameter/distance_to_earth)
}

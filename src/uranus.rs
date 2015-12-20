use angle;

/**
Computes Uranus's **equatorial semidiameter**

# Arguments

* ```distance_to_earth```: Uranus's distance from Earth *(AU)*
**/
pub fn Semidiameter(distance_to_earth: f64) -> f64 {
    angle::PureDegrees(0.0, 0.0, 35.02) / distance_to_earth
}

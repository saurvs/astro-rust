use angle;

/**
Computes Neptune's **equatorial semidiameter**

# Arguments

* ```distance_to_earth```: Neptune's distance from Earth *(AU)*
**/
pub fn Semidiameter(distance_to_earth: f64) -> f64 {
    angle::PureDegrees(0.0, 0.0, 33.5) / distance_to_earth
}

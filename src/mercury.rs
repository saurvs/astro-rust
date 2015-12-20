use angle;

/**
Computes Mercury's **equatorial semidiameter**

# Arguments

* ```distance_to_earth```: Mercury's distance from the Earth *(AU)*
**/
pub fn Semidiameter(distance_to_earth: f64) -> f64 {
    angle::PureDegrees(0.0, 0.0, 3.36) / distance_to_earth
}

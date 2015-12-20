use angle;

/**
Computes Mars's **equatorial semidiameter**

# Arguments

* ```distance_to_earth```: Mars's distance from Earth *(AU)*
**/
pub fn Semidiameter(distance_to_earth: f64) -> f64 {
    angle::PureDegrees(0.0, 0.0, 4.68) / distance_to_earth
}

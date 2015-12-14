use angle;

/**
Returns Mars's **equatorial semidiameter**

# Arguments

* ```distance_to_earth```: Mars's distance from Earth *(AU)*
**/
pub fn semidiameter(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 4.68) / distance_to_earth
}

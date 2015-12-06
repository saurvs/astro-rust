use angle;

/**
Returns Uranus's equatorial angular semidiameter

# Arguments

* ```distance_to_earth```: Uranus's distance from Earth (in AU)
**/
pub fn semidiameter(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 35.02) / distance_to_earth
}

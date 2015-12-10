use angle;

/**
Returns Mercury's equatorial angular semidiameter

# Arguments

* ```distance_to_earth```: Mercury's distance from the Earth (in AU)
**/
pub fn semidiameter(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 3.36) / distance_to_earth
}
use angle;

/**
Returns Venus's equatorial angular semidiameter

# Arguments

* ```distance_to_earth```: Venus's distance from  Earth (in AU)
**/
pub fn semidiameter(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 8.41) / distance_to_earth
}
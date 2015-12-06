use angle;

/**
Returns Jupiter's equatorial angular semidiameter

# Arguments

* ```distance_to_earth```: Jupiter's distance from Earth (in AU)
**/
pub fn eq_semidiameter(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 98.44) / distance_to_earth
}

/**
Returns Jupiter's polar angular semidiameter

# Arguments

* ```distance_to_earth```: Jupiter's distance from Earth (in AU)
**/
pub fn pol_semidiameter(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 92.06) / distance_to_earth
}

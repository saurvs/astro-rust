use angle;

/**
Returns Jupiter's **equatorial semidiameter**

# Arguments

* ```DistanceToEarth```: Jupiter's distance from Earth *(AU)*
**/
pub fn equatorial_semidiameter(DistanceToEarth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 98.44) / DistanceToEarth
}

/**
Returns Jupiter's **polar semidiameter**

# Arguments

* ```DistanceToEarth```: Jupiter's distance from Earth *(AU)*
**/
pub fn polar_semidiameter(DistanceToEarth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 92.06) / DistanceToEarth
}

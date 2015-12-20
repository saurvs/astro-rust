use angle;

fn equatorial_unit_semidiameter() -> f64 {
    angle::PureDegrees(0.0, 0.0, 82.73)
}

fn polar_unit_semidiameter() -> f64 {
    angle::PureDegrees(0.0, 0.0, 73.82)
}

/**
Computes Saturn's **polar semidiameter**

# Arguments

* ```distance_to_earth```: Saturn's distance from Earth *(AU)*
* ```earth_lat```: Earth's Saturnicentric latitude *(radians)*
**/
pub fn PolarSemidiameter(distance_to_earth: f64, earth_lat: f64) -> f64 {
    let a = equatorial_unit_semidiameter();
    let b = polar_unit_semidiameter();
    let k = 1.0 - (b / a).powi(2);
    (a / distance_to_earth) * (1.0 - k * earth_lat.cos().powi(2)).sqrt()
}

/**
Computes Saturn's **equatorial semidiameter**

# Arguments

* ```distance_to_earth```: Saturn's distance from Earth *(AU)*
**/
pub fn EquatorialSemidiameter(distance_to_earth: f64) -> f64 {
    equatorial_unit_semidiameter() / distance_to_earth
}

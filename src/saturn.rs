use angle;

fn equatorial_unit_semidiameter() -> f64 {
    angle::PureDegrees(0.0, 0.0, 82.73)
}

fn polar_unit_semidiameter() -> f64 {
    angle::PureDegrees(0.0, 0.0, 73.82)
}

/**
Returns Saturn's **polar semidiameter**

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
Returns Saturn's **equatorial semidiameter**

# Arguments

* ```distance_to_earth```: Saturn's distance from Earth *(AU)*
**/
pub fn EquatorialSemidiameter(distance_to_earth: f64) -> f64 {
    equatorial_unit_semidiameter() / distance_to_earth
}

/**
Returns Saturn's **orbital elements**

# Returned values

```(L, a, e, i, omega, pi, M, w)```

* ```L```: The mean longitude of Saturn *(radians)*
* ```a```: The semimajor axis of Saturn's orbit *(AU)*
* ```e```: The eccentricity of Saturn's orbit
* ```i```: The inclination on the plane of the Earth's ecliptic *(radians)*
* ```omega```: The longitude of the ascending node *(radians)*
* ```pi```: The longitude of the perihelion *(radians)*
* ```M```: The mean anomaly of Saturn *(radians)*
* ```w```: The argument of the perihelion *(radians)*

# Arguments

* ```T```: Time in Julian century
**/
pub fn OrbitalElements(T: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let TT = T * T;
    let TTT = TT * T;

    let L = 50.077444 + 1223.5110686*T + 0.00051908*TT - 0.00000003*TTT;
    let a = 9.554909192 - 0.0000021390*T + 0.000000004*TT;
    let e = 0.05554814 - 0.000346641*T - 0.0000006436*TTT + 0.0000000034*TTT;
    let i = 2.488879 - 0.0037362*T - 0.00001519*TT + 0.000000087*TTT;
    let omega = 113.665503 + 0.877088*T - 0.00012176*TT - 0.000002249*TTT;
    let pi = 93.057237 + 1.9637613*T + 0.00083753*TT + 0.000004928*TTT;

    (angle::LimitedTo360(L).to_radians(),
     a, e,
     angle::LimitedTo360(i).to_radians(),
     angle::LimitedTo360(omega).to_radians(),
     angle::LimitedTo360(pi).to_radians(),
     angle::LimitedTo360(L - pi).to_radians(),
     angle::LimitedTo360(pi - omega).to_radians()
    )
}

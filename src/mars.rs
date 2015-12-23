use angle;

/**
Computes Mars's **equatorial semidiameter**

# Arguments

* ```distance_to_earth```: Mars's distance from Earth *(AU)*
**/
pub fn Semidiameter(distance_to_earth: f64) -> f64 {
    angle::PureDegrees(0.0, 0.0, 4.68) / distance_to_earth
}

/**
Returns Mars's **orbital elements**

# Returned values

```(L, a, e, i, omega, pi, M, w)```

* ```L```: The mean longitude of Mars *(radians)*
* ```a```: The semimajor axis of Mars's orbit *(AU)*
* ```e```: The eccentricity of Mars's orbit
* ```i```: The inclination on the plane of the Earth's ecliptic *(radians)*
* ```omega```: The longitude of the ascending node *(radians)*
* ```pi```: The longitude of the perihelion *(radians)*
* ```M```: The mean anomaly of Mars *(radians)*
* ```w```: The argument of the perihelion *(radians)*

# Arguments

* ```T```: Time in Julian century
**/
pub fn OrbitalElements(T: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let TT = T * T;
    let TTT = TT * T;

    let L = 355.433 + 19141.6964471*T + 0.00031052*TT + 0.000000016*TTT;
    let a = 1.523679342;
    let e = 0.09340065 + 0.000090484*T - 0.0000000806*TTT - 0.00000000025*TTT;
    let i = 1.849726 - 0.0006011*T + 0.00001276*TT - 0.000000007*TTT;
    let omega = 49.558093 + 0.7720959*T + 0.00001557*TT - 0.000002267*TTT;
    let pi = 336.060234 + 1.8410449*T + 0.00013477*TT + 0.000000536*TTT;

    (angle::LimitedTo360(L).to_radians(),
     a, e,
     angle::LimitedTo360(i).to_radians(),
     angle::LimitedTo360(omega).to_radians(),
     angle::LimitedTo360(pi).to_radians(),
     angle::LimitedTo360(L - pi).to_radians(),
     angle::LimitedTo360(pi - omega).to_radians()
    )
}

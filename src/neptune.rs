use angle;

/**
Computes Neptune's **equatorial semidiameter**

# Arguments

* ```distance_to_earth```: Neptune's distance from Earth *(AU)*
**/
pub fn Semidiameter(distance_to_earth: f64) -> f64 {
    angle::PureDegrees(0.0, 0.0, 33.5) / distance_to_earth
}

/**
Returns Neptune's **orbital elements**

# Returned values

```(L, a, e, i, omega, pi, M, w)```

* ```L```: The mean longitude of Neptune *(radians)*
* ```a```: The semimajor axis of Neptune's orbit *(AU)*
* ```e```: The eccentricity of Neptune's orbit
* ```i```: The inclination on the plane of the Earth's ecliptic *(radians)*
* ```omega```: The longitude of the ascending node *(radians)*
* ```pi```: The longitude of the perihelion *(radians)*
* ```M```: The mean anomaly of Neptune *(radians)*
* ```w```: The argument of the perihelion *(radians)*

# Arguments

* ```T```: Time in Julian century
**/
pub fn OrbitalElements(T: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let TT = T * T;
    let TTT = TT * T;

    let L = 304.348665 + 219.8833092*T + 0.00030882*TT + 0.000000018*TTT;
    let a = 30.110386869 - 0.0000001663*T + 0.00000000069*TT;
    let e = 0.00945575 + 0.000006033*T - 0.00000000005*TTT;
    let i = 1.769953 - 0.0093082*T - 0.00000708*TT + 0.000000027*TTT;
    let omega = 131.784057 + 1.1022039*T + 0.00025952*TT - 0.000000637*TTT;
    let pi = 48.120276 + 1.4262957*T + 0.00038434*TT + 0.00000002*TTT;

    (angle::LimitedTo360(L).to_radians(),
     a, e,
     angle::LimitedTo360(i).to_radians(),
     angle::LimitedTo360(omega).to_radians(),
     angle::LimitedTo360(pi).to_radians(),
     angle::LimitedTo360(L - pi).to_radians(),
     angle::LimitedTo360(pi - omega).to_radians()
    )
}

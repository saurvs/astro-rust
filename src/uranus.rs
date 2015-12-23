use angle;

/**
Computes Uranus's **equatorial semidiameter**

# Arguments

* ```distance_to_earth```: Uranus's distance from Earth *(AU)*
**/
pub fn Semidiameter(distance_to_earth: f64) -> f64 {
    angle::PureDegrees(0.0, 0.0, 35.02) / distance_to_earth
}

/**
Returns Uranus's **orbital elements**

# Returned values

```(L, a, e, i, omega, pi, M, w)```

* ```L```: The mean longitude of Uranus *(radians)*
* ```a```: The semimajor axis of Uranus's orbit *(AU)*
* ```e```: The eccentricity of Uranus's orbit
* ```i```: The inclination on the plane of the Earth's ecliptic *(radians)*
* ```omega```: The longitude of the ascending node *(radians)*
* ```pi```: The longitude of the perihelion *(radians)*
* ```M```: The mean anomaly of Uranus *(radians)*
* ```w```: The argument of the perihelion *(radians)*

# Arguments

* ```T```: Time in Julian century
**/
pub fn OrbitalElements(T: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let TT = T * T;
    let TTT = TT * T;

    let L = 314.055005 + 429.8640561*T + 0.0003039*TT - 0.000000026*TTT;
    let a = 19.218446062 - 0.0000000372*T + 0.00000000098*TT;
    let e = 0.04638122 - 0.000027293*T + 0.0000000789*TTT + 0.00000000024*TTT;
    let i = 0.773197 + 0.0007744*T + 0.00003749*TT - 0.000000092*TTT;
    let omega = 74.005957 + 0.5211278*T + 0.00133947*TT + 0.000018484*TTT;
    let pi = 173.005291 + 1.486379*T + 0.00021406*TT + 0.000000434*TTT;

    (angle::LimitedTo360(L).to_radians(),
     a, e,
     angle::LimitedTo360(i).to_radians(),
     angle::LimitedTo360(omega).to_radians(),
     angle::LimitedTo360(pi).to_radians(),
     angle::LimitedTo360(L - pi).to_radians(),
     angle::LimitedTo360(pi - omega).to_radians()
    )
}

use angle;

/**
Returns Mercury's **equatorial semidiameter**

# Arguments

* ```distance_to_earth```: Mercury's distance from the Earth *(AU)*
**/
pub fn Semidiameter(distance_to_earth: f64) -> f64 {
    angle::PureDegrees(0.0, 0.0, 3.36) / distance_to_earth
}

/**
Returns Mercury's **orbital elements**

# Returned values

```(L, a, e, i, omega, pi, M, w)```

* ```L```: The mean longitude of Mercury *(radians)*
* ```a```: The semimajor axis of Mercury's orbit *(AU)*
* ```e```: The eccentricity of Mercury's orbit
* ```i```: The inclination on the plane of the Earth's ecliptic *(radians)*
* ```omega```: The longitude of the ascending node *(radians)*
* ```pi```: The longitude of the perihelion *(radians)*
* ```M```: The mean anomaly of Mercury *(radians)*
* ```w```: The argument of the perihelion *(radians)*

# Arguments

* ```T```: Time in Julian century
**/
pub fn OrbitalElements(T: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let TT = T * T;
    let TTT = TT * T;

    let L = 252.250906 + 149474.0722491*T + 0.0003035*TT + 0.000000018*TTT;
    let a = 0.038709831;
    let e = 0.20563175 + 0.000020407*T - 0.0000000283*TT + 0.00000000018*TTT;
    let i = 7.004986 + 0.0018215*T - 0.0000181*TT + 0.000000056*TTT;
    let omega = 48.330893 + 1.1861883*T + 0.00017542*TT + 0.000000215*TTT;
    let pi = 77.456119 + 1.5564776*T + 0.00029544*TT + 0.000000009*TTT;

    (angle::LimitedTo360(L).to_radians(),
     a, e,
     angle::LimitedTo360(i).to_radians(),
     angle::LimitedTo360(omega).to_radians(),
     angle::LimitedTo360(pi).to_radians(),
     angle::LimitedTo360(L - pi).to_radians(),
     angle::LimitedTo360(pi - omega).to_radians()
    )
}

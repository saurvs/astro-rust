use angle;

/**
Computes Venus's **equatorial semidiameter**

# Arguments

* ```distance_to_earth```: Venus's distance from  Earth *(AU)*
**/
pub fn Semidiameter(distance_to_earth: f64) -> f64 {
    angle::PureDegrees(0.0, 0.0, 8.41) / distance_to_earth
}

/**
Returns Venus's **orbital elements**

# Returned values

```(L, a, e, i, omega, pi, M, w)```

* ```L```: The mean longitude of Venus *(radians)*
* ```a```: The semimajor axis of Venus's orbit *(AU)*
* ```e```: The eccentricity of Venus's orbit
* ```i```: The inclination on the plane of the Earth's ecliptic *(radians)*
* ```omega```: The longitude of the ascending node *(radians)*
* ```pi```: The longitude of the perihelion *(radians)*
* ```M```: The mean anomaly of Venus *(radians)*
* ```w```: The argument of the perihelion *(radians)*

# Arguments

* ```T```: Time in Julian century
**/
pub fn OrbitalElements(T: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let TT = T * T;
    let TTT = TT * T;

    let L = 181.979801 + 58519.2130302*T + 0.00031014*TT + 0.000000015*TTT;
    let a = 0.72332982;
    let e = 0.00677192 - 0.000047765*T + 0.0000000981*TTT + 0.00000000046*TTT;
    let i = 3.394662 + 0.0010037*T - 0.00000088*TT - 0.000000007*TTT;
    let omega = 76.67992 + 0.9011206*T + 0.00040618*TT - 0.000000093*TTT;
    let pi = 131.563703 + 1.4022288*T - 0.00107618*TT - 0.000005678*TTT;

    (angle::LimitedTo360(L).to_radians(),
     a, e,
     angle::LimitedTo360(i).to_radians(),
     angle::LimitedTo360(omega).to_radians(),
     angle::LimitedTo360(pi).to_radians(),
     angle::LimitedTo360(L - pi).to_radians(),
     angle::LimitedTo360(pi - omega).to_radians()
    )
}

pub mod moon;
pub mod ecliptic;

use coordinates;
use angle;
use time;

/**
Returns the Earth's **orbital elements**

# Returned values

```(L, a, e, pi, M)```

* ```L```: The mean longitude of the Earth *(radians)*
* ```a```: The semimajor axis of the Earth's orbit *(AU)*
* ```e```: The eccentricity of the Earth's orbit
* ```pi```: The longitude of the perihelion *(radians)*
* ```M```: The mean anomaly of the Earth *(radians)*

# Arguments

* ```T```: Time in Julian century
**/
pub fn OrbitalElements(T: f64) -> (f64, f64, f64, f64, f64) {
    let TT = T * T;
    let TTT = TT * T;

    let L = 100.466457 + 36000.7698278*T + 0.00030322*TT + 0.00000002*TTT;
    let a = 1.000001018;
    let e = 0.01670863 - 0.000042037*T - 0.0000001267*TTT + 0.00000000014*TTT;
    let pi = 102.937348 + 1.7195366*T + 0.00045688*TT - 0.000000018*TTT;

    (angle::LimitedTo360(L).to_radians(),
     a, e,
     angle::LimitedTo360(pi).to_radians(),
     angle::LimitedTo360(L - pi).to_radians()
    )
}

/**
Returns the **flattening factor** of the Earth

Reference: [World Geodetic System 1984](https://confluence.qps.nl/pages/viewpage.action?pageId=29855173)
**/
pub fn Flattening() -> f64 {
    1.0 / 298.257223563
}

/**
Returns the **equatorial radius** of the Earth *(meters)*

Reference: [World Geodetic System 1984](https://confluence.qps.nl/pages/viewpage.action?pageId=29855173)
**/
pub fn EquatorialRadius() -> f64 {
    PolarRadius() / (1.0 - Flattening())
}

/**

Computes the **polar radius** of the Earth *(meters)*

Calculated using [```flattening()```](./fn.flattening.html) and
[```eq_radius()```](./fn.eq_radius.html)
**/
pub fn PolarRadius() -> f64 {
    6378137.0
}

/**
Computes the **eccentricity** of the Earth's **meridian**

Calculated using [```flattening()```](./fn.flattening.html)
**/
pub fn EccentricityOfMeridian() -> f64 {
    ((2.0 - Flattening()) * Flattening()).sqrt()
}

/**
Computes **angular distance** between two points on Earth's
surface

# Arguments

* ```p1```: Point 1
* ```p2```: Point 2
**/
pub fn AngularDistance(p1: coordinates::SurfacePoint, p2: coordinates::SurfacePoint) -> f64 {
    (p1.lat.sin() * p2.lat.sin() +
     p1.lat.cos() * p2.lat.cos() * (p1.long - p2.long).cos()
    ).acos()
}

/**
Computes a **low accuracy geodesic** between two points on the Earth's
surface *(meters)*

Assumes that the Earth is a sphere.

# Arguments

* ```p1```: Point 1
* ```p2```: Point 2
**/
pub fn ApproxGeodesic(p1: coordinates::SurfacePoint, p2: coordinates::SurfacePoint) -> f64 {
    6371.0 * AngularDistance(p1, p2)
}

/**
Computes a **high accuracy geodesic** between two points on the Earth's
surface *(meters)*

# Arguments

* ```p1```: Point 1
* ```p2```: Point 2
**/
pub fn Geodesic(p1: coordinates::SurfacePoint, p2: coordinates::SurfacePoint) -> f64 {
    let f = (p1.lat + p2.lat) / 2.0;
    let g = (p1.lat - p2.lat) / 2.0;
    let lam = (p1.long - p2.long) / 2.0;
    let s = (g.sin() * lam.cos()).powi(2) +
            (f.cos() * lam.sin()).powi(2);
    let c = (g.cos() * lam.cos()).powi(2) +
            (f.sin() * lam.sin()).powi(2);
    let om = ((s / c).sqrt()).atan();
    let r = (s * c).sqrt() / om;
    let d = 2.0 * om * EquatorialRadius();
    let h1 = (3.0 * r - 1.0) / (2.0 * c);
    let h2 = (3.0 * r + 1.0) / (2.0 * s);

    d * (1.0 +
         Flattening() * h1 * (f.sin() * g.cos()).powi(2) -
         Flattening() * h2 * (f.cos() * g.sin()).powi(2))

}

/**
Computes two quantities that are used elsewhere in the library

```Rho``` here denotes the geocentric radius vector, and ```Phi```
here denotes the geocentric latitude, both of an observer on the
Earth's surface.

# Arguments

* ```height```: Observer's height above sea level (in meters)
* ```geograph_lat```: Observer's geographical latitude *(radians)*
**/
pub fn RhoSinAndCosPhi(height: f64, geograph_lat: f64) -> (f64, f64) {
    let u = (geograph_lat.tan() * PolarRadius() / EquatorialRadius()).atan();
    let x = height / EquatorialRadius();
    let rho_sin_phi = (u.sin() * PolarRadius() / EquatorialRadius()) +
                      (geograph_lat.sin() * x);
    let rho_cos_phi = u.cos() + (geograph_lat.cos() * x);

    (rho_sin_phi, rho_cos_phi)
}

/**
Computes **nutation correction** *(radians)* for ecliptical longitude
and obliquity

# Return variables

Computes the nutation correction *(radians)*, that needs to be added to
the ecliptical longitude and the obliquity of the ecliptic to adjust for
the Earth's nutation. Nutation does *not* affect ecliptical latitude.

```NutationCorrections() -> (nutation_in_longitude, nutation_in_obliquity)```

# Arguments

```julian_ephemeris_day```: Julian Ephemeris day
**/
pub fn NutationCorrections(julian_ephemeris_day: f64) -> (f64, f64) {

    struct terms(i8, i8, i8, i8, i8, f64, f64, f64, f64);
    let tuple_terms = [
        terms(0, 0, 0, 0, 1, -171996.0, -174.2, 92025.0, 8.9),
        terms(-2, 0, 0, 2, 2, -13187.0, -1.6, 5736.0, -3.1),
        terms(0, 0, 0, 2, 2, -2274.0, -0.2, 977.0, -0.5),
        terms(0, 0, 0, 0, 2, 2062.0, 0.2, -895.0, 0.5),
        terms(0, 1, 0, 0, 0, 1426.0, -3.4, 54.0, -0.1),
        terms(0, 0, 1, 0, 0, 712.0, 0.1, -7.0, 0.0),
        terms(-2, 1, 0, 2, 2, -517.0, 1.2, 224.0, -0.6),
        terms(0, 0, 0, 2, 1, -386.0, -0.4, 200.0, 0.0),
        terms(0, 0, 1, 2, 2, -301.0, 0.0, 129.0, -0.1),
        terms(-2, -1, 0, 2, 2, 217.0, -0.5, -95.0, 0.3),
        terms(-2, 0, 1, 0, 0, -158.0, 0.0, 0.0, 0.0),
        terms(-2, 0, 0, 2, 1, 129.0, 0.1, -70.0, 0.0),
        terms(0, 0, -1, 2, 2, 123.0, 0.0, -53.0, 0.0),
        terms(2, 0, 0, 0, 0, 63.0, 0.0, 0.0, 0.0),
        terms(0, 0, 1, 0, 1, 63.0, 0.1, -33.0, 0.0),
        terms(2, 0, -1, 2, 2, -59.0, 0.0, 26.0, 0.0),
        terms(0, 0, -1, 0, 1, -58.0, -0.1, 32.0, 0.0),
        terms(0, 0, 1, 2, 1, -51.0, 0.0, 27.0, 0.0),
        terms(-2, 0, 2, 0, 0, 48.0, 0.0, 0.0, 0.0),
        terms(0, 0, -2, 2, 1, 46.0, 0.0, -24.0, 0.0),
        terms(2, 0, 0, 2, 2, -38.0, 0.0, 16.0, 0.0),
        terms(0, 0, 2, 0, 0, 29.0, 0.0, 0.0, 0.0),
        terms(-2, 0, 1, 2, 2, 29.0, 0.0, -12.0, 0.0),
        terms(0, 0, 0, 2, 0, 26.0, 0.0, 0.0, 0.0),
        terms(-2, 0, 0, 2, 0, -22.0, 0.0, 0.0, 0.0),
        terms(0, 0, -1, 2, 1, 21.0, 0.0, -10.0, 0.0),
        terms(0, 2, 0, 0, 0, 17.0, -0.1, 0.0, 0.0),
        terms(2, 0, -1, 0, 1, 16.0, 0.0, -8.0, 0.0),
        terms(-2, 2, 0, 2, 2, -16.0, 0.1, 7.0, 0.0),
        terms(0, 1, 0, 0, 1, -15.0, 0.0, 9.0, 0.0),
        terms(-2, 0, 1, 0, 1, -13.0, 0.0, 7.0, 0.0),
        terms(0, -1, 0, 0, 1, -12.0, 0.0, 6.0, 0.0),
        terms(0, 0, 2, -2, 0, 11.0, 0.0, 0.0, 0.0),
        terms(2, 0, -1, 2, 1, -10.0, 0.0, 5.0, 0.0),
        terms(2, 0, 1, 2, 2, -8.0, 0.0, 3.0, 0.0),
        terms(0, 1, 0, 2, 2, 7.0, 0.0, -3.0, 0.0),
        terms(-2, 1, 1, 0, 0, -7.0, 0.0, 0.0, 0.0),
        terms(0, -1, 0, 2, 2, -7.0, 0.0, 3.0, 0.0),
        terms(2, 0, 0, 2, 1, -7.0, 0.0, 3.0, 0.0),
        terms(2, 0, 1, 0, 0, 6.0, 0.0, 0.0, 0.0),
        terms(-2, 0, 2, 2, 2, 6.0, 0.0, -3.0, 0.0),
        terms(-2, 0, 1, 2, 1, 6.0, 0.0, -3.0, 0.0),
        terms(2, 0, -2, 0, 1, -6.0, 0.0, 3.0, 0.0),
        terms(2, 0, 0, 0, 1, -6.0, 0.0, 3.0, 0.0),
        terms(0, -1, 1, 0, 0, 5.0, 0.0, 0.0, 0.0),
        terms(-2, -1, 0, 2, 1, -5.0, 0.0, 3.0, 0.0),
        terms(-2, 0, 0, 0, 1, -5.0, 0.0, 3.0, 0.0),
        terms(0, 0, 2, 2, 1, -5.0, 0.0, 3.0, 0.0),
        terms(-2, 0, 2, 0, 1, 4.0, 0.0, 0.0, 0.0),
        terms(-2, 1, 0, 2, 1, 4.0, 0.0, 0.0, 0.0),
        terms(0, 0, 1, -2, 0, 4.0, 0.0, 0.0, 0.0),
        terms(-1, 0, 1, 0, 0, -4.0, 0.0, 0.0, 0.0),
        terms(-2, 1, 0, 0, 0, -4.0, 0.0, 0.0, 0.0),
        terms(1, 0, 0, 0, 0, -4.0, 0.0, 0.0, 0.0),
        terms(0, 0, 1, 2, 0, 3.0, 0.0, 0.0, 0.0),
        terms(0, 0, -2, 2, 2, -3.0, 0.0, 0.0, 0.0),
        terms(-1, -1, 1, 0, 0, -3.0, 0.0, 0.0, 0.0),
        terms(0, 1, 1, 0, 0, -3.0, 0.0, 0.0, 0.0),
        terms(0, -1, 1, 2, 2, -3.0, 0.0, 0.0, 0.0),
        terms(2, -1, -1, 2, 2, -3.0, 0.0, 0.0, 0.0),
        terms(0, 0, 3, 2, 2, -3.0, 0.0, 0.0, 0.0),
        terms(2, -1, 0, 2, 2, -3.0, 0.0, 0.0, 0.0),
    ];

    let t = time::JulianCentury(julian_ephemeris_day);

    let M1 = angle::LimitedTo360((134.96298 + t*(477198.867398 + t*(0.0086972 + t/5620.0)))).to_radians();
    let M = angle::LimitedTo360((357.52772 + t*(35999.050340 - t*(0.0001603 + t/300000.0)))).to_radians();
    let D = angle::LimitedTo360((297.85036 + t*(445267.11148 - t*(0.0019142 + t/189474.0)))).to_radians();
    let F = angle::LimitedTo360((93.27191 + t*(483202.017538 - t*(-0.0036825 - t/327270.0)))).to_radians();
    let om = angle::LimitedTo360((125.04452 - t*(1934.136261 - t*(0.0020708 + t/450000.0)))).to_radians();

    let mut nut_in_long = 0.0;
    let mut nut_in_obl = 0.0;

    for x in tuple_terms.iter() {
        let arg = (x.0 as f64) * D +
                  (x.1 as f64) * M +
                  (x.2 as f64) * M1 +
                  (x.3 as f64) * F +
                  (x.4 as f64) * om;
        nut_in_long += ((x.5 as f64) + t*(x.6 as f64)) * arg.sin() * (0.0001 / 3600.0);
        nut_in_obl += ((x.7 as f64) + t*(x.8 as f64)) * arg.cos() * (0.0001 / 3600.0);
    }

    (nut_in_long.to_radians(), nut_in_obl.to_radians())
}

/**
Computes **equation of time** *(radians)*

# Arguments

* ```jed```: Julian Ephemeris day
* ```sun_asc```: Right ascension of the Sun *(radians)*
* ```nut_log```: Nutation correction for longitude *(radians)*
* ```tru_obl```: *True* obliquity of the ecliptic *(radians)*
**/
pub fn EquationOfTime(jed: f64, sun_asc: f64, nut_long: f64, tru_obl: f64) -> f64 {
    let t = time::JulianCentury(jed) / 10.0;
    let L = angle::LimitedTo360(
            280.4664567 +
            t * (360007.6982779 +
            t * (0.030328 +
            t * (1.0/49931.0 -
            t * (1.0/15300.0 +
            t * (1.0/2000000.0)
            ))))                 );

    (L - 0.0057183 -
     sun_asc.to_degrees() +
     nut_long.to_degrees()*tru_obl.cos()
    ).to_radians()
}

#[macro_export]
macro_rules! EquationOfTime {
    ($x: expr, $y: expr) => {{
            let (nut_long, nut_obl) = earth::NutationCorrection($x);
            let true_obl = earth::MeanObliquity($x) + nut_obl;
            earth::EquationOfTime($x, $y, nut_long, true_obl)
    }};
}

pub fn GeocentricCoordsOfPlanet(L: f64, B: f64, R: f64, L0: f64, B0: f64, R0: f64) -> (f64, f64, f64) {
    let x = R*B.cos()*L.cos() - R0*B0.cos()*L0.cos();
    let y = R*B.cos()*L.sin() - R0*B0.cos()*L0.sin();
    let z = R*B.sin() - R0*B0.sin();

    (y.atan2(x),
     z/(x*x + y*y).sqrt(),
     0.0057755183 * (x*x + y*y + z*z).sqrt())
}

pub fn what(obl_eclp: f64, long_asc_node: f64, inc: f64) {
    let sin_obl_eclp = obl_eclp.sin();
    let cos_obl_eclp = obl_eclp.cos();
    let cos_long_asc_node = long_asc_node.cos();
    let sin_long_asc_node = long_asc_node.sin();
    let cos_inc = inc.cos();
    let sin_inc = inc.sin();

    let f = cos_long_asc_node;
    let g = sin_long_asc_node*cos_obl_eclp;
    let h = sin_long_asc_node*sin_obl_eclp;
    let p = -1.0*sin_long_asc_node*sin_inc;
    let q = cos_long_asc_node*cos_inc*cos_obl_eclp - sin_inc*sin_obl_eclp;
    let r = cos_long_asc_node*cos_inc*sin_obl_eclp + sin_inc*cos_obl_eclp;

    let A = f.atan2(p);
    let B = g.atan2(q);
    let C = h.atan2(r);
    let a = (f*f + p*p).sqrt();
    let b = (g*g + q*q).sqrt();
    let c = (h*h + r*r).sqrt();

    let x = r * a * (A + perih_arg + v);
    let y = r * b * (B + perih_arg + v);
    let z = r * c * (C + perih_arg + v);

    let xi = X + x;
    let nu = Y + y;
    let et = Z + z;

    let mut asc = nu.atan2(xi);
    let dec = et.atan2((xi*xi + nu*nu).sqrt());

    if asc < 0.0 {
        asc += 360f64.to_radians();
    }
}

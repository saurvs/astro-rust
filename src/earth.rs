use coordinates;
use angle;
use time;

/**
Returns the flattening factor of the Earth

Reference: [World Geodetic System 1984](https://confluence.qps.nl/pages/viewpage.action?pageId=29855173)
**/
pub fn flattening() -> f64 {
    1.0 / 298.257223563
}

/**
Returns the equatorial radius of the Earth (in meters)

Reference: [World Geodetic System 1984](https://confluence.qps.nl/pages/viewpage.action?pageId=29855173)
**/
pub fn eq_radius() -> f64 {
    polar_radius() / (1.0 - flattening())
}

/**

Returns the polar radius of the Earth (in meters)

Calculated using [```flattening()```](./fn.flattening.html) and
[```eq_radius()```](./fn.eq_radius.html)
**/
pub fn polar_radius() -> f64 {
    6378137.0
}

/**
Returns the eccentricity of the Earth's meridian

Calculated using [```flattening()```](./fn.flattening.html)
**/
pub fn ecc_mer() -> f64 {
    ((2.0 - flattening()) * flattening()).sqrt()
}

/**
Returns the angular distance between two points on Earth's
surface

# Arguments

* ```p1```: Point 1
* ```p2```: Point 2
**/
pub fn angular_dist(p1: coordinates::surf_point, p2: coordinates::surf_point) -> f64 {
    (p1.lat.sin() * p2.lat.sin() +
     p1.lat.cos() * p2.lat.cos() * (p1.long - p2.long).cos()
    ).acos()
}

/**
Returns a low accuracy geodesic distance between two points on Earth's
surface (in meters)

Assumes that the Earth is a sphere.

# Arguments

* ```p1```: Point 1
* ```p2```: Point 2
**/

pub fn geodesic_approx_dist(p1: coordinates::surf_point, p2: coordinates::surf_point) -> f64 {
    6371.0 * angular_dist(p1, p2)
}

/**
Returns a high accuracy geodesic distance between two points on Earth's
surface (in meters)

# Arguments

* ```p1```: Point 1
* ```p2```: Point 2
**/

pub fn geodesic_dist(p1: coordinates::surf_point, p2: coordinates::surf_point) -> f64 {
    let f = (p1.lat + p2.lat) / 2.0;
    let g = (p1.lat - p2.lat) / 2.0;
    let lam = (p1.long - p2.long) / 2.0;
    let s = (g.sin() * lam.cos()).powi(2) +
            (f.cos() * lam.sin()).powi(2);
    let c = (g.cos() * lam.cos()).powi(2) +
            (f.sin() * lam.sin()).powi(2);
    let om = ((s / c).sqrt()).atan();
    let r = (s * c).sqrt() / om;
    let d = 2.0 * om * eq_radius();
    let h1 = (3.0 * r - 1.0) / (2.0 * c);
    let h2 = (3.0 * r + 1.0) / (2.0 * s);

    d * (1.0 +
         flattening() * h1 * (f.sin() * g.cos()).powi(2) -
         flattening() * h2 * (f.cos() * g.sin()).powi(2))

}

/**
Returns two quantities that are used elsewhere in the library

```rho``` here denotes the geocentric radius vector, and ```phi```
here denotes the geocentric latitude, both of an observer on the
Earth's surface.

# Arguments

* ```height```: Observer's height above sea level (in meters)
* ```geograph_lat```: Observer's geographical latitude (in radians)
**/

pub fn rho_sin_and_cos_phi(height: f64, geograph_lat: f64) -> (f64, f64) {
    let u = (geograph_lat.tan() * polar_radius() / eq_radius()).atan();
    let x = height / eq_radius();
    let rho_sin_phi = (u.sin() * polar_radius() / eq_radius()) +
                      (geograph_lat.sin() * x);
    let rho_cos_phi = u.cos() + (geograph_lat.cos() * x);

    (rho_sin_phi, rho_cos_phi)
}

/**
Returns the nutation correction (in radians) for longitude and obliquity

# Return variables

Returns the nutation correction (in radians), that needs to be added to
the longitude and the obliquity of the ecliptic to adjust for the
Earth's nutation. Nutation does not affect latitudes.

```nutation() -> (nutation_in_longitude, nutation_in_obliquity)```

# Arguments

```julian_ephemeris_day```: Time in Julian Ephemeris Days
**/

pub fn nutation(julian_ephemeris_day: f64) -> (f64, f64) {

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

    let t = time::julian_century(julian_ephemeris_day);

    let M1 = angle::limited_to_360((134.96298 + t*(477198.867398 + t*(0.0086972 + t/5620.0)))).to_radians();
    let M = angle::limited_to_360((357.52772 + t*(35999.050340 - t*(0.0001603 + t/300000.0)))).to_radians();
    let D = angle::limited_to_360((297.85036 + t*(445267.11148 - t*(0.0019142 + t/189474.0)))).to_radians();
    let F = angle::limited_to_360((93.27191 + t*(483202.017538 - t*(-0.0036825 - t/327270.0)))).to_radians();
    let om = angle::limited_to_360((125.04452 - t*(1934.136261 - t*(0.0020708 + t/450000.0)))).to_radians();

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
Returns the mean obliquity of the Earth's ecliptic (in radians)

The obliquity of the ecliptic, or the inclination of the Earth's axis
of rotation, is the angle between the equator and the ecliptic.

The *mean* obliquity isn't corrected for nutation.
To obtain the *true* obliquity, use [```nutation()```]
(./fn.nutation.html) to get the nutation correction for obliquity,
and add it to the *mean* obliquity.

# Arguments

```julian_ephemeris_day```: Time in Julian Ephemeris days
**/

pub fn mean_obliquity(jed: f64) -> (f64) {
    let u = time::julian_century(jed) / 100.0;

    (angle::pure_degrees(23.0, 26.0, 21.448) -
    u * (angle::pure_degrees(0.0, 0.0, 4680.93) +
    u * (angle::pure_degrees(0.0, 0.0, 1.55) +
    u * (angle::pure_degrees(0.0, 0.0, 1999.25) -
    u * (angle::pure_degrees(0.0, 0.0, 51.38) -
    u * (angle::pure_degrees(0.0, 0.0, 249.67) +
    u * (angle::pure_degrees(0.0, 0.0, 39.05) +
    u * (angle::pure_degrees(0.0, 0.0, 7.12) -
    u * (angle::pure_degrees(0.0, 0.0, 27.87) +
    u * (angle::pure_degrees(0.0, 0.0, 5.79) +
    u * angle::pure_degrees(23.0, 26.0, 2.45)
    )))))))))).to_radians()
}

/**
Returns the equation of time (in radians)

# Arguments

* ```jed```: Time in Julian Ephemeris days
* ```sun_asc```: Right ascension of the Sun (in radians)
* ```nut_log```: Nutation correction for longitude (in radians)
* ```tru_obl```: *True* obliquity of the ecliptic (in radians)
**/

pub fn equation_of_time(jed: f64, sun_asc: f64, nut_long: f64, tru_obl: f64) -> f64 {
    let t = time::julian_century(jed) / 10.0;
    let L = angle::limited_to_360(
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

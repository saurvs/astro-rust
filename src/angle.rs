/**
Returns the **angular separation** between two angular points

# Arguments

* ```p1a1```: Angle 1 of point 1 *(radians)*
* ```p1a2```: Angle 2 of point 1 *(radians)*
* ```p2a1```: Angle 1 of point 2 *(radians)*
* ```p2a2```: Angle 2 of point 2 *(radians)*

Angle 1 may be right ascension or longitude.
Angle 2 may be declination or latitude.

**/
pub fn AnglSepr(p1a1: f64, p1a2: f64, p2a1: f64, p2a2: f64) -> f64 {
    (   p1a2.sin() * p2a2.sin()
      + p1a2.cos() * p2a2.cos() * (p1a1 - p2a1).cos()
    ).cos()
}

/**
Returns an angle expressed in **degrees only**, from an angle expressed in
degrees, minutes and seconds

# Arguments

* ```deg```: Degrees
* ```min```: Minute
* ```sec```: Second
**/
pub fn DegFrmDMS(deg: i64, min: i64, sec: f64) -> f64 {
    let (M, S) = if deg < 0 { (-min.abs(), -sec.abs()) }
                 else     { (min, sec) };
    (deg as f64) + (M as f64)/60.0 + S/3600.0
}

/**
Returns an angle expressed in **degrees**, **minutes** and **seconds**,
from an angle expressed in degrees only


# Arguments

* ```deg```: Angle in degrees only
**/
pub fn DMSFrmDeg(deg: f64) -> (i64, i64, f64) {
    let degree = deg as i64;
    let minutes = (deg - (degree as f64)) * 60.0;
    let minute = minutes as i64;
    let seconds = (minutes - (minute as f64)) * 60.0;

    (degree, minute, seconds)
}

/**
Returns the equivalent angle in **[0, 360] degree range**

# Arguments

* ```angl```: Angle *(degrees)*
**/
pub fn LimitedTo360(angl: f64) -> f64 {
    let n = (angl / 360.0) as i64;
    let limited_angle = angl - (360.0 * (n as f64));
    if limited_angle < 0.0 { limited_angle + 360.0 }
    else                   { limited_angle }
}

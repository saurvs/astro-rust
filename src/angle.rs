use coordinates;

/**
Returns **angular separation** between two equatorial points

# Arguments

* ```p1```: Equatorial point 1 *(radians)*
* ```p2```: Equatorial point 2 *(radians)*
**/
pub fn angular_sep(p1: coordinates::equator_point, p2: coordinates::equator_point) -> f64 {
    (p1.dec.sin() * p2.dec.sin() +
     p1.dec.cos() * p2.dec.cos() * (p1.asc - p2.asc).cos()
    ).cos()
}

/**
Returns an angle expressed in **degrees only**

# Arguments

* ```d```: Degrees
* ```m```: Minute
* ```s```: Second
**/
pub fn pure_degrees(d: f64, mut m: f64, mut s: f64) -> f64 {
    if d < 0.0 {
        m = -1.0 * m.abs();
        s = -1.0 * s.abs();
    }
    d + (m / 60.0) + (s / 3600.0)
}

/**
Returns the equivalent angle in **[0, 360] degree range**

# Arguments

* ```angle```: Angle *(degrees)*
**/
pub fn limited_to_360(angle: f64) -> f64 {
    let n = (angle / 360.0) as i64;
    let mut limited_angle = angle - (360.0 * (n as f64));
    if limited_angle < 0.0 {
        limited_angle += 360.0;
    }
    limited_angle
}

/**
Checks wether an angle is "small"

Checks wether an angle in degrees is sufficiently "small", for an
an arbritrary definition of "small", chosen for adequate purposes
for this library.

# Arguments

* ```angle```: Angle (in degrees)
**/
fn small_angle(angle: f64) -> bool {
    if angle < 0.003 { return true; }
    false
}

use coordinates;

/**
Returns **angular separation** between two equatorial points

# Arguments

* ```p1```: Equatorial point 1 *(radians)*
* ```p2```: Equatorial point 2 *(radians)*
**/
pub fn AngularSep(p1: coordinates::EquatorialPoint, p2: coordinates::EquatorialPoint) -> f64 {
    (p1.dec.sin() * p2.dec.sin() +
     p1.dec.cos() * p2.dec.cos() * (p1.asc-p2.asc).cos()
    ).cos()
}

/**
Returns an angle expressed in **degrees only**

# Arguments

* ```d```: Degrees
* ```m```: Minute
* ```s```: Second
**/
pub fn PureDegrees(d: f64, m: f64, s: f64) -> f64 {
    let (M, S) = if d < 0.0 { (-m.abs(), -s.abs()) }
                 else       { (m, s) };
    d + M/60.0 + S/3600.0
}

/**
Returns the equivalent angle in **[0, 360] degree range**

# Arguments

* ```angle```: Angle *(degrees)*
**/
pub fn LimitedTo360(angle: f64) -> f64 {
    let n = (angle/360.0) as i64;
    let limited_angle = angle - (360.0 * (n as f64));
    if limited_angle < 0.0 { limited_angle + 360.0 }
    else                   { limited_angle }
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
    if angle < 0.003 { true }
    else             { false }
}

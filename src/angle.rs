use coordinates;

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
pub fn AngularSeparation(p1a1: f64, p1a2: f64, p2a1: f64, p2a2: f64) -> f64 {
    (   p1a2.sin() * p2a2.sin()
      + p1a2.cos() * p2a2.cos() * (p1a1 - p2a1).cos()
    ).cos()
}

/**
Returns an angle expressed in **degrees only**

# Arguments

* ```d```: Degrees
* ```m```: Minute
* ```s```: Second
**/
pub fn PureDegrees(d: i64, m: i64, s: f64) -> f64 {
    let (M, S) = if d < 0 { (-m.abs(), -s.abs()) }
                 else     { (m, s) };
    (d as f64) + (M as f64)/60.0 + S/3600.0
}

/**
Returns an angle expressed in **degrees**, **minutes** and **seconds**

# Arguments

* ```degrees```: Angle in degrees only
**/
pub fn DegreesMinutesSecondsFromPureDegrees(degrees: f64) -> (i64, i64, f64) {
    let degree = degrees as i64;
    let minutes = (degrees - (degree as f64)) * 60.0;
    let minute = minutes as i64;
    let seconds = (minutes - (minute as f64)) * 60.0;

    (degree, minute, seconds)
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

fn small_angle(angle: f64) -> bool {
    if angle < 0.003 { true }
    else             { false }
}

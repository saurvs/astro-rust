//! Angles for astronomy

/**
Returns the **angular separation** between two angular points

# Arguments

* `p1a1`: Angle 1 of point 1 *| in radians*
* `p1a2`: Angle 2 of point 1 *| in radians*
* `p2a1`: Angle 1 of point 2 *| in radians*
* `p2a2`: Angle 2 of point 2 *| in radians*

Angle 1 may be right ascension or longitude.
Angle 2 may be declination or latitude.

**/
pub fn angl_sepr(p1a1: f64, p1a2: f64, p2a1: f64, p2a2: f64) -> f64 {
    (   p1a2.sin() * p2a2.sin()
      + p1a2.cos() * p2a2.cos() * (p1a1 - p2a1).cos()
    ).acos()
}

/**
Returns an angle in **degrees with decimals**, from an angle
expressed in **degrees, arcminutes** and **arcseconds**

# Arguments

* `deg`: Degrees
* `min`: Arcminutes
* `sec`: Arcseconds
**/
pub fn deg_frm_dms(deg: i64, min: i64, sec: f64) -> f64 {
    let (M, S) = if deg < 0 { (-min.abs(), -sec.abs()) }
                 else     { (min, sec) };
    (deg as f64) + (M as f64)/60.0 + S/3600.0
}

/**
Returns an **angle** expressed in **degrees, arcminutes** and
**arcseconds**, from an angle in **degrees with decimals**

# Returns

`(deg, min, sec)`

* `deg`: Degrees
* `min`: Arcminutes
* `sec`: Arcseconds

# Arguments

* `deg`: Angle in degrees with decimals
**/
pub fn dms_frm_deg(deg: f64) -> (i64, i64, f64) {
    let degree = deg as i64;
    let minutes = (deg - (degree as f64)) * 60.0;
    let minute = minutes as i64;
    let seconds = (minutes - (minute as f64)) * 60.0;

    (degree, minute, seconds)
}

/**
Returns an **angle** expressed in **hours, minutes** and
**seconds**, from an angle in **degrees with decimals**

# Returns

`(deg, min, sec)`

* `hour`: Hours
* `min`: Minutes
* `sec`: Seconds

# Arguments

* `deg`: Angle in degrees with decimals
**/
pub fn hms_frm_deg(deg: f64) -> (i64, i64, f64) {
    let hours = deg / 15.0;
    let hour = hours as i64;

    let minutes = (hours - (hour as f64)) * 60.0;
    let minute = minutes as i64;

    let seconds = (minutes - (minute as f64)) * 60.0;

    (hour, minute, seconds)
}

/**
Returns an angle in **degrees with decimals**, from an angle
expressed in **hours, minutes** and **seconds**

# Arguments

* `hours`: Hours
* `min`: Minutes
* `sec`: Seconds
**/
pub fn deg_frm_hms(hour: i64, min: i64, sec: f64) -> f64 {
    15.0 * ((hour as f64) + (min as f64)/60.0 + sec/3600.0)
}

/**
Returns the equivalent angle in **[0, 360]** degree range

# Arguments

* `angl`: Angle *(degrees)*
**/
pub fn limit_to_360(angl: f64) -> f64 {
    let n = (angl / 360.0) as i64;
    let limited_angle = angl - (360.0 * (n as f64));
    if limited_angle < 0.0 { limited_angle + 360.0 }
    else                   { limited_angle }
}

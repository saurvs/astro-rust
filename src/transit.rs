//! Time of rise, transit and set for a celestial body

use angle;
use coords;
use interpol;
use std::*;

/// Represents a **celestial body** in transit
pub enum TransitBody {
    /// A star or a planet
    StarOrPlanet,
    /// The Sun
    Sun,
    /// The Moon
    Moon
}

/// Represents a **transit type**
pub enum TransitType {
    /// Rise
    Rise,
    /// Transit at zenith
    Transit,
    /// Set
    Set
}

/**
Returns the **time** of **transit** for a celestial body

# Returns

* ```(hour, minute, second)```: Time of the transit on the day of interest

# Arguments

* ```transit_type```: A ```TransitType```
* ```transit_body```: The ```TransitBody```
* ```geograph_point```: Geographic point of the observer *| in radians*

Let ```JD``` be the Julian (Ephemeris) day of interest,

* ```eq_point1```: Equatorial point of the transit body on ```JD - 1``` *| in radians*
* ```eq_point2```: Equatorial point of the transit body on ```JD``` *| in radians*
* ```eq_point3```: Equatorial point of the transit body on ```JD + 1``` *| in radians*
* ```app_green_sidr```: Apparent sidereal time at Greenwhich on ```JD``` *| in radians*
* ```deltaT```: Delta T on ```JD```
* ```moon_eq_hz_parallax```: Equatorial horizontal parallax of the Moon on ```JD```
                             *| in radians*. *Pass a meaningfull value here only when*
                             ```TransitBody::Moon``` *is passed for* ```transit_body```.

**/
pub fn Time(
    transit_type: &TransitType,
    transit_body: &TransitBody,
    geograph_point: &coords::GeographPoint,
    eq_point1: &coords::EqPoint,
    eq_point2: &coords::EqPoint,
    eq_point3: &coords::EqPoint,
    app_green_sidr: f64,
    deltaT: f64,
    moon_eq_hz_parallax: f64) -> (i64, i64, f64) {

    let h0 = match transit_body {
        &TransitBody::StarOrPlanet => -0.5667_f64.to_radians(),
        &TransitBody::Sun => -0.8333_f64.to_radians(),
        &TransitBody::Moon =>  0.7275*moon_eq_hz_parallax
                             - 0.5667_f64.to_radians(),
    };

    let mut H0 = ( (h0.sin() - geograph_point.lat.sin()*eq_point2.dec.sin())
                   / (geograph_point.lat.cos() * eq_point2.dec.cos())
                 ).acos();
    H0 = angle::LimitTo360(H0.to_degrees()).to_radians();

    let rad360 = 2.0 * f64::consts::PI;
    let mut m = m(&transit_type, H0, eq_point2.asc, geograph_point.long, app_green_sidr, rad360);
    let theta0 = app_green_sidr + m*360.985647_f64.to_radians();

    let d = m + deltaT/86400.0;

    let asc = interpol::ThreeValues(eq_point1.asc, eq_point2.asc, eq_point3.asc, d);

    let dec = match transit_type {
        &TransitType::Transit => 0.0,
        &TransitType::Rise    => interpol::ThreeValues(eq_point1.dec, eq_point2.dec, eq_point3.dec, d),
        &TransitType::Set     => interpol::ThreeValues(eq_point1.dec, eq_point2.dec, eq_point3.dec, d)
    };

    let mut H = coords::HrAnglFrmObserverLong(theta0, geograph_point.long, asc).to_degrees();
    H = angle::LimitTo360(H);
    if H > 180.0 { H = H - 360.0; }
    H = H.to_radians();

    let h = match transit_type {
        &TransitType::Transit => 0.0,
        &TransitType::Rise    => coords::AltFrmEqCoords(H, dec, geograph_point.lat),
        &TransitType::Set     => coords::AltFrmEqCoords(H, dec, geograph_point.lat)
    };

    m += match transit_type {
        &TransitType::Transit => -H/rad360,
        &TransitType::Rise    => (h - h0) / (rad360 * dec.cos() * geograph_point.lat.cos() * H.sin()),
        &TransitType::Set     => (h - h0) / (rad360 * dec.cos() * geograph_point.lat.cos() * H.sin())
    };

    let h = 24.0 * m;
    let hour = h as i64;
    let m = (h - (hour as f64)) * 60.0;
    let minute = m as i64;
    let second = (m - (minute as f64)) * 60.0;

    (hour, minute, second)
}

fn m(transit_type: &TransitType, H0: f64, asc: f64, L: f64, Theta0: f64, rad360: f64) -> f64 {
    let mut m = (asc + L - Theta0)/rad360;

    let p = H0/rad360;
    m += match transit_type {
        &TransitType::Transit => 0.0,
        &TransitType::Rise    => -p,
        &TransitType::Set     =>  p
    };

    if      m < 0.0 { m += 1.0 }
    else if m > 1.0 { m -= 1.0 }

    m
}

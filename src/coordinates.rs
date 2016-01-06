/// Represents a point on **Earth's surface**
pub struct SurfacePoint {
    /// Latitude
    pub lat: f64,
    /// Longitude
    pub long: f64,
}

/// Represents a point on the celestial sphere, using the **equatorial
/// coordinate system**
pub struct EquatorialPoint {
    /// Right ascension
    pub asc: f64,
    /// Declination
    pub dec: f64,
}

//-------------------------------------------------------------------

/**
Returns hour angle from sidereal time at Greenwhich

# Arguments

* ```green_sid```: Sidereal time at Greenwhich
* ```obv_long```: Observer's longitude *(radians)*
* ```right_asc```: Right ascension *(radians)*
**/
pub fn HourAngleFromGreenwhichSidereal(green_sid: f64, obv_long: f64, right_asc: f64) -> f64 {
    green_sid - obv_long - right_asc
}

/**
Returns hour angle from local sidereal time

# Arguments

* ```loc_sid```: Local sidereal time
* ```right_asc```: Right ascension *(radians)*
**/
pub fn HourAngleFromLocalSidereal(loc_sid: f64, right_asc: f64) -> f64 {
    loc_sid - right_asc
}

//-------------------------------------------------------------------
// Ecliptical coordinates to equatorial coordinates

/**
Returns **ecliptical longitude** *(radians)* from **equatorial coordinates**

# Arguments

* ```right_asc```: Right ascension *(radians)*
* ```dec```: Declination *(radians)*
* ```oblq_ecl```: If ```right_asc``` and ```dec``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. *(radians)*
**/
pub fn EclipticalLongitudeFromEquatorialCoords(right_asc: f64, dec: f64, oblq_ecl: f64,) -> f64 {
    ((right_asc.sin()*oblq_ecl.cos() + dec.tan()*oblq_ecl.sin())).atan2(right_asc.cos())
}

/**
Returns **ecliptical latitude** *(radians)* from **equatorial coordinates**

# Arguments

* ```right_asc```: Right ascension *(radians)*
* ```dec```: Declination *(radians)*
* ```oblq_ecl```: If ```right_asc``` and ```dec``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. *(radians)*
**/
pub fn EclipticalLatitudeFromEquatorialCoords(right_asc: f64, dec: f64, oblq_ecl: f64) -> f64 {
    (dec.sin()*oblq_ecl.cos() - dec.cos()*oblq_ecl.sin()*right_asc.sin()).asin()
}

/**
Returns **ecliptical coordinates** *(radians)* from **equatorial coordinates**

# Returns

```(ecliptical_longitude, ecliptical_latitude)```

* ```ecliptical_longitude```: Ecliptical longitude *(radians)*
* ```ecliptical_latitude```: Ecliptical latitude *(radians)*

# Arguments

* ```$x```: Right ascension *(radians)*
* ```$y```: Declination *(radians)*
* ```$z```: If ```$x``` and ```$y``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. *(radians)*
**/
#[macro_export]
macro_rules! EclipticalCoordsFromEquatorialCoords {
    ($x: expr, $y: expr, $z: expr) => {{
        (astro::earth::eclip_long_from_equatorial($x, $y, $z),
         astro::earth::eclip_lat_from_equatorial($x, $y, $z))
    }};
}

//-------------------------------------------------------------------
// Equatorial coordinates to ecliptical coordinates

/**
Returns **right ascension** *(radians)* from **ecliptical coordinates**

# Arguments

* ```ecl_long```: Ecliptical longitude *(radians)*
* ```ecl_lat```: Ecliptical latitude *(radians)*
* ```oblq_ecl```: If ```ecl_long``` and ```ecl_lat``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. *(radians)*
**/
pub fn RightAscensionFromEclipticalCoords(ecl_long: f64, ecl_lat: f64, oblq_ecl: f64) -> f64 {
    ((ecl_long.sin()*oblq_ecl.cos() - ecl_lat.tan()*oblq_ecl.sin())).atan2(ecl_long.cos())
}

/**
Returns **declination** *(radians)* from **ecliptical coordinates**

# Arguments

* ```ecl_long```: Ecliptical longitude *(radians)*
* ```ecl_lat```: Ecliptical latitude *(radians)*
* ```oblq_ecl```: If ```ecl_long``` and ```ecl_lat``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. *(radians)*
**/
pub fn DeclinationFromEclipticalCoords(ecl_long: f64, ecl_lat: f64, oblq_ecl: f64) -> f64 {
    (ecl_lat.sin()*oblq_ecl.cos() - ecl_lat.cos()*oblq_ecl.sin()*ecl_long.sin()).asin()
}

/**
Returns **equatorial coordinates** *(radians)* from **ecliptical coordinates**

# Returns

```(equatorial_longitude, equatorial_latitude)```

* ```equatorial_longitude```: Equatorial longitude *(radians)*
* ```equatorial_latitude```: Equatorial latitude *(radians)*

# Arguments

* ```$x```: Ecliptical longitude *(radians)*
* ```$y```: Ecliptical latitude *(radians)*
* ```$z```: If ```$x``` and ```$y``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. *(radians)*
**/
#[macro_export]
macro_rules! EquatorialCoordsFromEclipticalCoords {
    ($x: expr, $y: expr, $z: expr) => {{
        (astro::earth::right_ascen_from_eclip($x, $y, $z),
         astro::earth::declin_from_eclip($x, $y, $z))
    }};
}

//-------------------------------------------------------------------
// local horizontal coordinates from equatorial coordinates

pub fn AzimuthFromEquatorialCoords(hour_angle: f64, declination: f64, observer_lat: f64) -> f64 {
    hour_angle.sin().atan2(hour_angle.cos()*observer_lat.sin() - declination.tan()*observer_lat.cos())
}

pub fn AltitudeFromEquatorialCoords(hour_angle: f64, declination: f64, observer_lat: f64) -> f64 {
    (observer_lat.sin()*declination.sin() + observer_lat.cos()*declination.cos()*hour_angle.cos()).asin()
}

//-------------------------------------------------------------------
// horizontal coordinates to equatorial coordinates

pub fn HourAngleFromEquatorialCoords(azimuth: f64, altitude: f64, obv_lat: f64) -> f64 {
    azimuth.sin().atan2(azimuth.cos()*obv_lat.sin() + altitude.tan()*obv_lat.cos())
}

pub fn eq_dec_from_equatorialial(azimuth: f64, altitude: f64, obv_lat: f64) -> f64 {
    (obv_lat.sin()*altitude.sin() - obv_lat.cos()*azimuth.cos() * azimuth.cos()).asin()
}

//-------------------------------------------------------------------
// Equatorial coordinates to galactic coordinates

/**
Returns **galactic longitude** *(radians)* from **equatorial coordinates**

# Arguments

* ```right_asc```: Right ascension *(radians)*
* ```dec```: Declination *(radians)*
**/
pub fn GalacticLongitudeFromEquatorialCoords(right_asc: f64, dec: f64) -> f64 {
      5.28835
    - ((3.35539549 - right_asc).sin())
      .atan2((  (3.35539549 - right_asc).cos() * 0.4782202_f64.sin()
              - dec.tan()                      * 0.4782202_f64.cos()
            ))
}

/**
Returns **galactic latitude** *(radians)* from **equatorial coordinates**

# Arguments

* ```right_asc```: Right ascension *(radians)*
* ```dec```: Declination *(radians)*
**/
pub fn GalacticLatitudeFromEquatorialCoords(right_asc: f64, dec: f64) -> f64 {
      dec.sin() * 0.4782202_f64.sin()
    + dec.cos() * 0.4782202_f64.cos() * (3.35539549 - right_asc).asin()
}

/**
Returns **galactic coordinates** *(radians)* from **equatorial coordinates**

# Returns

```(galactic_longitude, galactic_latitude) = galactic_from_equatorial!()```
# Returns

```(galactic_longitude, galactic_latitude)```

* ```galactic_longitude```: Galactic longitude *(radians)*
* ```galactic_latitude```: Declination *(radians)*

# Arguments

* ```$x```: Right ascension *(radians)*
* ```$y```: Declination *(radians)*
**/
#[macro_export]
macro_rules! GalacticCoordsFromEquatorialCoords {
    ($x: expr, $y: expr) => {{
        (astro::earth::galac_long_from_equatorial($x, $y),
         astro::earth::galac_lat_from_equatorial($x, $y))
    }};
}

//-------------------------------------------------------------------
// Galactic coordinates to equatorial coordinates

/**
Returns **right ascension** *(radians)* from **galactic coordinates**

# Arguments

* ```gal_long```: Galactic longitude *(radians)*
* ```gal_lat```: Galactic latitude *(radians)*
**/
pub fn RightAscensionFromGalacticCoords(gal_long: f64, gal_lat: f64) -> f64 {
      0.21380283
    + ((gal_long - 2.14675).sin())
      .atan2((   (gal_long - 2.14675).cos() * 0.4782202_f64.sin()
               - gal_lat.tan()              * 0.4782202_f64.cos()
            ))
}

/**
Returns **declination** *(radians)* from **galactic coordinates**

# Arguments

* ```gal_long```: Galactic longitude *(radians)*
* ```gal_lat```: Galactic latitude *(radians)*
**/
pub fn DeclinationFromGalacticCoords(gal_long: f64, gal_lat: f64) -> f64 {
      gal_lat.sin() * 0.4782202_f64.sin()
    + gal_lat.cos() * 0.4782202_f64.cos() * (gal_long - 2.14675).asin()
}

/**
Returns **equatorial coordinates** from **galactic coordinates**

# Returns

```(right_ascension, declination)```

* ```right_ascension```: Right ascension *(radians)*
* ```declination```: Declination *(radians)*

# Arguments

* ```$x```: Galactic longitude *(radians)*
* ```$y```: Galactic latitude *(radians)*
**/
#[macro_export]
macro_rules! EquatorialCoordsFromGalacticCoords {
    ($x: expr, $y: expr) => {{
        (astro::earth::right_ascen_from_galac($x, $y),
         astro::earth::declin_from_galac($x, $y))
    }};
}

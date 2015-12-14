/// Represents a point on **Earth's surface**
pub struct surface_point {
    /// Latitude
    pub lat: f64,
    /// Longitude
    pub long: f64,
}

/// Represents a point on the celestial sphere, using the **equatorial
/// coordinate system**
pub struct equator_point {
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
pub fn hour_angle_from_greenwhich_sidereal(green_sid: f64, obv_long: f64, right_asc: f64) -> f64 {
    green_sid - obv_long - right_asc
}

/**
Returns hour angle from local sidereal time

# Arguments

* ```loc_sid```: Local sidereal time
* ```right_asc```: Right ascension *(radians)*
**/
pub fn hour_angle_from_local_sidereal(loc_sid: f64, right_asc: f64) -> f64 {
    loc_sid - right_asc
}

//-------------------------------------------------------------------

/**
Returns **ecliptical longitude** *(radians)* from **equatorial coordinates**

# Arguments

* ```right_asc```: Right ascension *(radians)*
* ```dec```: Declination *(radians)*
* ```oblq_ecl```: If ```right_asc``` and ```dec``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. *(radians)*
**/
pub fn eclip_long_from_equatorial(right_asc: f64, dec: f64, oblq_ecl: f64,) -> f64 {
    ((right_asc.sin() * oblq_ecl.cos() + dec.tan() * oblq_ecl.sin())).atan2(right_asc.cos())
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
pub fn eclip_lat_from_equatorial(right_asc: f64, dec: f64, oblq_ecl: f64) -> f64 {
    (dec.sin() * oblq_ecl.cos() - dec.cos() * oblq_ecl.sin() * right_asc.sin()).asin()
}

/**
Returns **ecliptical coordinates** *(radians)* from **equatorial coordinates**

# Return values

```(ecliptical_longitude, ecliptical_latitude) = ecliptical_from_equatorial!()```

# Arguments

* ```$x```: Right ascension *(radians)*
* ```$y```: Declination *(radians)*
* ```$z```: If ```$x``` and ```$y``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. *(radians)*
**/
#[macro_export]
macro_rules! ecliptical_from_equatorial {
    ($x: expr, $y: expr, $z: expr) => {{
        (earth::eclip_long_from_equatorial($x, $y, $z),
         earth::eclip_lat_from_equatorial($x, $y, $z))
    }};
}

//-------------------------------------------------------------------

/**
Returns **right ascension** *(radians)* from **ecliptical coordinates**

# Arguments

* ```ecl_long```: Ecliptical longitude *(radians)*
* ```ecl_lat```: Ecliptical latitude *(radians)*
* ```oblq_ecl```: If ```ecl_long``` and ```ecl_lat``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. *(radians)*
**/
pub fn right_ascen_from_eclip(ecl_long: f64, ecl_lat: f64, oblq_ecl: f64) -> f64 {
    ((ecl_long.sin() * oblq_ecl.cos() - ecl_lat.tan() * oblq_ecl.sin())).atan2(ecl_long.cos())
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
pub fn declin_from_eclip(ecl_long: f64, ecl_lat: f64, oblq_ecl: f64) -> f64 {
    (ecl_lat.sin() * oblq_ecl.cos() - ecl_lat.cos() * oblq_ecl.sin() * ecl_long.sin()).asin()
}

/**
Returns **equatorial coordinates** *(radians)* from **ecliptical coordinates**

# Return values

```(equatorial_longitude, equatorial_latitude) = equatorial_from_ecliptical!()```

# Arguments

* ```$x```: Ecliptical longitude *(radians)*
* ```$y```: Ecliptical latitude *(radians)*
* ```$z```: If ```$x``` and ```$y``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. *(radians)*
**/
#[macro_export]
macro_rules! equatorial_from_ecliptical {
    ($x: expr, $y: expr, $z: expr) => {{
        (earth::right_ascen_from_eclip($x, $y, $z),
         earth::declin_from_eclip($x, $y, $z))
    }};
}

//-------------------------------------------------------------------
// local horizontal coordinates from equatorial coordinates

pub fn azimuth_from_equatorial(hour_angle: f64, eq_dec: f64, obv_lat: f64) -> f64 {
    hour_angle.sin().atan2(hour_angle.cos() * obv_lat.sin() - eq_dec.tan() * obv_lat.cos())
}

pub fn altitude_from_equatorial(hour_angle: f64, eq_dec: f64, obv_lat: f64) -> f64 {
    (obv_lat.sin() * eq_dec.sin() + obv_lat.cos() * eq_dec.cos() * hour_angle.cos()).asin()
}

//-------------------------------------------------------------------
// horizontal coordinates to equatorial coordinates

pub fn hour_angle_from_equatorialial(azimuth: f64, altitude: f64, obv_lat: f64) -> f64 {
    azimuth.sin().atan2(azimuth.cos() * obv_lat.sin() + altitude.tan() * obv_lat.cos())
}

pub fn eq_dec_from_equatorialial(azimuth: f64, altitude: f64, obv_lat: f64) -> f64 {
    (obv_lat.sin() * altitude.sin() - obv_lat.cos() * azimuth.cos() * azimuth.cos()).asin()
}

//-------------------------------------------------------------------

/**
Returns **galactic longitude** *(radians)* from **equatorial coordinates**

# Arguments

* ```right_asc```: Right ascension *(radians)*
* ```dec```: Declination *(radians)*
**/
pub fn galac_long_from_equatorial(right_asc: f64, dec: f64) -> f64 {
    5.28835 - ((3.35539549 - right_asc).sin())
            .atan2(((3.35539549 - right_asc).cos() * (0.4782202_f64).sin() -
                     dec.tan() * (0.4782202_f64).cos()
                   ))
}

/**
Returns **galactic latitude** *(radians)* from **equatorial coordinates**

# Arguments

* ```right_asc```: Right ascension *(radians)*
* ```dec```: Declination *(radians)*
**/
pub fn galac_lat_from_equatorial(right_asc: f64, dec: f64) -> f64 {
    dec.sin() * (0.4782202_f64).sin() +
    dec.cos() * (0.4782202_f64).cos() * (3.35539549 - right_asc).asin()
}

/**
Returns **galactic coordinates** *(radians)* from **equatorial coordinates**

# Return values

```(galactic_longitude, galactic_latitude) = galactic_from_equatorial!()```

# Arguments

* ```$x```: Right ascension *(radians)*
* ```$y```: Declination *(radians)*
**/
#[macro_export]
macro_rules! galactic_from_equatorial {
    ($x: expr, $y: expr) => {{
        (earth::galac_long_from_equatorial($x, $y),
         earth::galac_lat_from_equatorial($x, $y))
    }};
}

//-------------------------------------------------------------------

/**
Returns **right ascension** *(radians)* from **galactic coordinates**

# Arguments

* ```gal_long```: Galactic longitude *(radians)*
* ```gal_lat```: Galactic latitude *(radians)*
**/
pub fn right_ascen_from_galac(gal_long: f64, gal_lat: f64) -> f64 {
    0.21380283 + ((gal_long - 2.14675).sin())
            .atan2(((gal_long - 2.14675).cos() * (0.4782202_f64).sin() -
                     gal_lat.tan() * (0.4782202_f64).cos()
                   ))
}

/**
Returns **declination** *(radians)* from **galactic coordinates**

# Arguments

* ```gal_long```: Galactic longitude *(radians)*
* ```gal_lat```: Galactic latitude *(radians)*
**/
pub fn declin_from_galac(gal_long: f64, gal_lat: f64) -> f64 {
    gal_lat.sin() * (0.4782202_f64).sin() +
    gal_lat.cos() * (0.4782202_f64).cos() * (gal_long - 2.14675).asin()
}

/**
Returns **equatorial coordinates** *(radians)* from **galactic coordinates**

# Return values

```(right_ascension, declination) = equatorial_from_galactic!()```

# Arguments

* ```$x```: Galactic longitude *(radians)*
* ```$y```: Galactic latitude *(radians)*
**/
#[macro_export]
macro_rules! equatorial_from_galactic {
    ($x: expr, $y: expr) => {{
        (earth::right_ascen_from_galac($x, $y),
         earth::declin_from_galac($x, $y))
    }};
}

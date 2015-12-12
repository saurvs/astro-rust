/// Represents a point on Earth's surface
pub struct surf_point {
    /// Latitude
    pub lat: f64,
    /// Longitude
    pub long: f64,
}

/// Represents a point on the celestial sphere, using the equatorial
/// coordinate system
pub struct equa_point {
    /// Right ascension
    pub asc: f64,
    /// Declination
    pub dec: f64,
}

/// Returns the obliquity of the ecliptic for the epoch J2000.0
pub fn oblq_2000() -> f64 {
    23.4392911_f64.to_radians()
}

/// Returns the obliquity of the ecliptic for the epoch J1950.0
pub fn oblq_1950() -> f64 {
    23.4457889_f64.to_radians()
}

/**
Returns hour angle from sidereal time at Greenwhich

# Arguments

* ```green_sid```: Sidereal time at Greenwhich
* ```obv_long```: Observer's longitude (in radians)
* ```right_asc```: Right ascension (in radians)
**/
pub fn hour_angle_from_greenwhich_sid(green_sid: f64, obv_long: f64, right_asc: f64) -> f64 {
    green_sid - obv_long - right_asc
}

/**
Returns hour angle from local sidereal time

# Arguments

* ```loc_sid```: Local sidereal time
* ```right_asc```: Right ascension (in radians)
**/
pub fn hour_angle_from_loc_sid(loc_sid: f64, right_asc: f64) -> f64 {
    loc_sid - right_asc
}

/**
Returns ecliptical longitude (in radians) from equatorial coordinates

# Arguments

* ```right_asc```: Right ascension (in radians)
* ```dec```: Declination (in radians)
* ```oblq_ecl```: If ```right_asc``` and ```dec``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. (in radians)
**/
pub fn eclip_long_from_equa(right_asc: f64, dec: f64, oblq_ecl: f64,) -> f64 {
    ((right_asc.sin() * oblq_ecl.cos() + dec.tan() * oblq_ecl.sin())).atan2(right_asc.cos())
}

/**
Returns ecliptical latitude (in radians) from equatorial coordinates

# Arguments

* ```right_asc```: Right ascension (in radians)
* ```dec```: Declination (in radians)
* ```oblq_ecl```: If ```right_asc``` and ```dec``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. (in radians)
**/
pub fn eclip_lat_from_equa(right_asc: f64, dec: f64, oblq_ecl: f64) -> f64 {
    (dec.sin() * oblq_ecl.cos() - dec.cos() * oblq_ecl.sin() * right_asc.sin()).asin()
}

/**
Returns right ascension (in radians) from ecliptical coordinates

# Arguments

* ```ecl_long```: Ecliptical longitude (in radians)
* ```ecl_lat```: Ecliptical latitude (in radians)
* ```oblq_ecl```: If ```ecl_long``` and ```ecl_lat``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. (in radians)
**/
pub fn right_asc_from_eclip(ecl_long: f64, ecl_lat: f64, oblq_ecl: f64) -> f64 {
    ((ecl_long.sin() * oblq_ecl.cos() - ecl_lat.tan() * oblq_ecl.sin())).atan2(ecl_long.cos())
}

/**
Returns declination (in radians) from ecliptical coordinates

# Arguments

* ```ecl_long```: Ecliptical longitude (in radians)
* ```ecl_lat```: Ecliptical latitude (in radians)
* ```oblq_ecl```: If ```ecl_long``` and ```ecl_lat``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. (in radians)
**/
pub fn dec_from_eclip(ecl_long: f64, ecl_lat: f64, oblq_ecl: f64) -> f64 {
    (ecl_lat.sin() * oblq_ecl.cos() - ecl_lat.cos() * oblq_ecl.sin() * ecl_long.sin()).asin()
}

//-------------------------------------------------------------------
// local horizontal coordinates from equatorial coordinates

pub fn azimuth(hour_angle: f64, eq_dec: f64, obv_lat: f64) -> f64 {
    hour_angle.sin().atan2(hour_angle.cos() * obv_lat.sin() - eq_dec.tan() * obv_lat.cos())
}

pub fn altitude(hour_angle: f64, eq_dec: f64, obv_lat: f64) -> f64 {
    (obv_lat.sin() * eq_dec.sin() + obv_lat.cos() * eq_dec.cos() * hour_angle.cos()).asin()
}

//-------------------------------------------------------------------
// horizontal coordinates to equatorial coordinates

pub fn hour_angle(azimuth: f64, altitude: f64, obv_lat: f64) -> f64 {
    azimuth.sin().atan2(azimuth.cos() * obv_lat.sin() + altitude.tan() * obv_lat.cos())
}

pub fn eq_dec(azimuth: f64, altitude: f64, obv_lat: f64) -> f64 {
    (obv_lat.sin() * altitude.sin() - obv_lat.cos() * azimuth.cos() * azimuth.cos()).asin()
}

/**
Returns galactic longitude (in radians) from equatorial coordinates

# Arguments

* ```right_asc```: Right ascension (in radians)
* ```dec```: Declination (in radians)
**/
pub fn gal_long_from_equa(right_asc: f64, dec: f64) -> f64 {
    5.28835 - ((3.35539549 - right_asc).sin())
            .atan2(((3.35539549 - right_asc).cos() * (0.4782202_f64).sin() -
                     dec.tan() * (0.4782202_f64).cos()
                   ))
}

/**
Returns galactic latitude (in radians) from equatorial coordinates

# Arguments

* ```right_asc```: Right ascension (in radians)
* ```dec```: Declination (in radians)
**/
pub fn gal_lat_from_equa(right_asc: f64, dec: f64) -> f64 {
    dec.sin() * (0.4782202_f64).sin() +
    dec.cos() * (0.4782202_f64).cos() * (3.35539549 - right_asc).asin()
}

/**
Returns right ascension (in radians) from galactic coordinates

# Arguments

* ```gal_long```: Galactic longitude (in radians)
* ```gal_lat```: Galactic latitude (in radians)
**/
pub fn right_asc_from_gal(gal_long: f64, gal_lat: f64) -> f64 {
    0.21380283 + ((gal_long - 2.14675).sin())
            .atan2(((gal_long - 2.14675).cos() * (0.4782202_f64).sin() -
                     gal_lat.tan() * (0.4782202_f64).cos()
                   ))
}

/**
Returns declination (in radians) from galactic coordinates

# Arguments

* ```gal_long```: Galactic longitude (in radians)
* ```gal_lat```: Galactic latitude (in radians)
**/
pub fn dec_from_gal(gal_long: f64, gal_lat: f64) -> f64 {
    gal_lat.sin() * (0.4782202_f64).sin() +
    gal_lat.cos() * (0.4782202_f64).cos() * (gal_long - 2.14675).asin()
}

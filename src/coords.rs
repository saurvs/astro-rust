//! Transform between coordinate systems

use angle;
use time;

/// Represents a point on the geographical surface of the Earth
pub struct GeographPoint {
    /// Geographical longitude
    pub long: f64,
    /// Geographical latitude
    pub lat: f64,
}

impl GeographPoint {
    pub fn angl_sepr(&self, other_point: &GeographPoint) -> f64 {
        angle::angl_sepr(self.long, self.lat,
                      other_point.long, other_point.lat)
    }
}

/// Represents a point in the equatorial coordinate system
pub struct EqPoint {
    /// Right ascension
    pub asc: f64,
    /// Declination
    pub dec: f64,
}

impl EqPoint {
    pub fn angl_sepr(&self, other_point: &EqPoint) -> f64 {
        angle::angl_sepr(self.asc, self.dec,
                      other_point.asc, other_point.dec)
    }
}

/// Represents a point in the ecliptic coordinate system
pub struct EclPoint {
    /// Ecliptic longitude
    pub long: f64,
    /// Ecliptic latitude
    pub lat: f64,
}

impl EclPoint {
    pub fn angl_sepr(&self, other_point: &EclPoint) -> f64 {
        angle::angl_sepr(self.long, self.lat,
                      other_point.long, other_point.lat)
    }
}

//-------------------------------------------------------------------

/**
Returns the **hour angle** from **geographical longitude** and **Greenwhich
sidereal time**

# Returns

* `hour_angle`: Hour angle *| in radians*

# Arguments

* `green_sidreal`: Sidereal time at Greenwhich *| in radians*
* `observer_long`: Observer's geographical longitude *| in radians*
* `asc`: Right ascension *| in radians*
**/
pub fn HrAnglFrmObserverLong(green_sidreal: f64, observer_long: f64, asc: f64) -> f64 {
    green_sidreal - observer_long - asc
}

/**
Returns the **hour angle** from **local sidereal time** and **right
ascension**

# Returns

* `hour_angle`: Hour angle *| in radians*

# Arguments

* `local_sidreal`: Local sidereal time *| in radians*
* `asc`: Right ascension *| in radians*
**/
pub fn HrAnglFrmLocSidr(local_sidreal: f64, asc: f64) -> f64 {
    local_sidreal - asc
}

//-------------------------------------------------------------------
// Ecliptic coordinates to equatorial coordinates

/**
Returns the **ecliptic longitude** from **equatorial** coordinates

# Returns

* `ecl_long`: Ecliptic longitude *| in radians*

# Arguments

* `asc`: Right ascension *| in radians*
* `dec`: Declination *| in radians*
* `oblq_eclip`: If `asc` and `dec` are corrected for
                    nutation, then *true* obliquity. If not, then
                    *mean* obliquity. *| in radians*
**/
pub fn EclLongFrmEq(asc: f64, dec: f64, oblq_eclip: f64,) -> f64 {
    (   asc.sin() * oblq_eclip.cos()
      + dec.tan() * oblq_eclip.sin()
    ).atan2(asc.cos())
}

/**
Returns the **ecliptic latitude** from **equatorial** coordinates

# Returns

* `ecl_lat`: Ecliptic latitude *| in radians*

# Arguments

* `asc`: Right ascension *| in radians*
* `dec`: Declination *| in radians*
* `oblq_eclip`: If `asc` and `dec` are corrected for
                    nutation, then *true* obliquity. If not, then
                    *mean* obliquity. *| in radians*
**/
pub fn EclLatFrmEq(asc: f64, dec: f64, oblq_eclip: f64) -> f64 {
    (   dec.sin() * oblq_eclip.cos()
      - dec.cos() * oblq_eclip.sin() * asc.sin()
    ).asin()
}

/**
Returns **ecliptic coordinates** from **equatorial** coordinates

# Returns

`(ecl_long, ecl_lat)`

* `ecl_long`: Ecliptic longitude *| in radians*
* `ecl_lat`: Ecliptic latitude *| in radians*

# Arguments

* `$asc`: Right ascension *| in radians*
* `$dec`: Declination *| in radians*
* `$oblq_eclip`: If `$asc` and `$dec` are corrected for
                     nutation, then *true* obliquity. If not, then
                     *mean* obliquity. *| in radians*
**/

macro_rules! EclFrmEq {
    ($asc: expr, $dec: expr, $oblq_eclip: expr) => {{
        (astro::coords::EclLongFrmEq($asc, $dec, $oblq_eclip),
         astro::coords::EclLatFrmEq($asc, $dec, $oblq_eclip))
    }};
}

//-------------------------------------------------------------------
// Equatorial coordinates to ecliptic coordinates

/**
Returns the **right ascension** from **ecliptic** coordinates

# Returns

* `asc`: Right ascension *| in radians*

# Arguments

* `ecl_long`: Ecliptic longitude *| in radians*
* `ecl_lat`: Ecliptic latitude *| in radians*
* `oblq_eclip`: If `ecl_long` and `ecl_lat` are corrected
                    for nutation, then *true* obliquity. If not, then
                    *mean* obliquity. *| in radians*
**/
pub fn AscFrmEcl(ecl_long: f64, ecl_lat: f64, oblq_eclip: f64) -> f64 {
    (   ecl_long.sin() * oblq_eclip.cos()
      - ecl_lat.tan()  * oblq_eclip.sin()
    ).atan2(ecl_long.cos())
}

/**
Returns the **declination** from **ecliptic** coordinates

# Returns

* `dec`: Declination *| in radians*

# Arguments

* `ecl_long`: Ecliptic longitude *| in radians*
* `ecl_lat`: Ecliptic latitude *| in radians*
* `oblq_eclip`: If `ecl_long` and `ecl_lat` are corrected
                    for nutation, then *true* obliquity. If not, then
                    *mean* obliquity. *| in radians*
**/
pub fn DecFrmEcl(ecl_long: f64, ecl_lat: f64, oblq_eclip: f64) -> f64 {
    (   ecl_lat.sin() * oblq_eclip.cos()
      + ecl_lat.cos() * oblq_eclip.sin() * ecl_long.sin()
    ).asin()
}

/**
Returns **equatorial** coordinates from **ecliptic** coordinates

# Returns

`(asc, dec)`

* `asc`: Right ascension *| in radians*
* `dec`: Declination *| in radians*

# Arguments

* `$ecl_long`: Ecliptic longitude *| in radians*
* `$ecl_lat`: Ecliptic latitude *| in radians*
* `$oblq_eclip`: If `$x` and `$y` are corrected for
                     nutation, then *true* obliquity. If not, then
                     *mean* obliquity. *| in radians*
**/

macro_rules! EqFrmEcl {
    ($ecl_long: expr, $y: expr, $oblq_eclip: expr) => {{
        (astro::coords::AscFrmEcl($ecl_long, $ecl_lat, $oblq_eclip),
         astro::coords::DecFrmEcl($ecl_long, $ecl_lat, $oblq_eclip))
    }};
}

//-------------------------------------------------------------------
// Local horizontal coordinates from equatorial coordinates

/**
Returns the **azimuth** from **equatorial** coordinates

# Returns

* `az`: Azimuth *| in radians*

# Arguments

* `hour_angle`: Hour angle *| in radians*
* `dec`: Declination *| in radians*
* `observer_lat`: Observer's geographical latitude *| in radians*
**/
pub fn AzFrmEqCoords(hour_angle: f64, dec: f64, observer_lat: f64) -> f64 {
    hour_angle.sin()
    .atan2(   hour_angle.cos()  * observer_lat.sin()
            - dec.tan() * observer_lat.cos()
          )
}

/**
Returns the **altitude** from **equatorial** coordinates

# Returns

* `alt`: Altitude *| in radians*

# Arguments

* `hour_angle`: Hour angle *| in radians*
* `dec`: Declination *| in radians*
* `observer_lat`: Observer's geographical latitude *| in radians*
**/
pub fn AltFrmEqCoords(hour_angle: f64, dec: f64, observer_lat: f64) -> f64 {
    (   observer_lat.sin() * dec.sin()
      + observer_lat.cos() * dec.cos() * hour_angle.cos()
    ).asin()
}

/**
Returns **local horizontal** coordinates from **equatorial** coordinates

# Returns

`(az, alt)`

* `az`: Azimuth *| in radians*
* `alt`: Altitude *| in radians*

# Arguments

* `$hour_angle`: Hour angle *| in radians*
* `$dec`: Declination *| in radians*
* `$observer_lat`: Observer's geographical latitude *| in radians*
**/

macro_rules! LocHzFrmEq {
    ($hour_angle: expr, $dec: expr, $observer_lat: expr) => {{
        (astro::coords::AzFrmEqCoords($hour_angle, $dec, $observer_lat),
         astro::coords::AltFrmEqCoords($hour_angle, $dec, $observer_lat))
    }};
}

//-------------------------------------------------------------------
// Equatorial coordinates from local horizontal coordinates

/**
Returns the **hour angle** from **local horizontal** coordinates

# Returns

* `hour_angle`: Hour angle *| in radians*

# Arguments

* `az`: Azimuth *| in radians*
* `alt`: Altitude *| in radians*
* `observer_lat`: Observer's geographical latitude *| in radians*
**/
pub fn HourAnglFrmHz(az: f64, alt: f64, observer_lat: f64) -> f64 {
    az.sin()
    .atan2(   az.cos()  * observer_lat.sin()
            + alt.tan() * observer_lat.cos()
          )
}

/**
Returns the **declination** from **local horizontal** coordinates

# Returns

* `dec`: Declination *| in radians*

# Arguments

* `az`: Azimuth *| in radians*
* `alt`: Altitude *| in radians*
* `observer_lat`: Observer's geographical latitude *| in radians*
**/
pub fn DecFrmHz(az: f64, alt: f64, observer_lat: f64) -> f64 {
    (   observer_lat.sin() * alt.sin()
      - observer_lat.cos() * az.cos() * az.cos()
    ).asin()
}

//-------------------------------------------------------------------
// Equatorial coordinates to galactic coordinates

/**
Returns the **galactic longitude** from **equatorial** coordinates

# Returns

* `gal_long`: Galactic longitude *| in radians*

# Arguments

* `asc`: Right ascension *| in radians*
* `dec`: Declination *| in radians*

The equatorial coordinates passed are assumed to be referred to the
standard equinox of B1950.0.
**/
pub fn GalLongFrmEq(asc: f64, dec: f64) -> f64 {
      303_f64.to_radians()
    - (192.25_f64.to_radians() - asc).sin()
      .atan2(   27.4_f64.to_radians().sin() * (192.25_f64.to_radians() - asc).cos()
              - 27.4_f64.to_radians().cos() * dec.tan()
            )
}

/**
Returns the **galactic latitude** from **equatorial** coordinates

# Returns

* `gal_lat`: Galactic latitude *| in radians*

# Arguments

* `asc`: Right ascension *| in radians*
* `dec`: Declination *| in radians*

The equatorial coordinates passed are assumed to be referred to the
standard equinox of B1950.0.
**/
pub fn GalLatFrmEq(asc: f64, dec: f64) -> f64 {
    (   dec.sin() * 27.4_f64.to_radians().sin()
      + dec.cos() * 27.4_f64.to_radians().cos() * (192.25_f64.to_radians() - asc).cos()
    ).asin()
}

/**
Returns **galactic** coordinates from **equatorial** coordinates

# Returns

`(gal_long, gal_lat)`

* `gal_long`: Galactic longitude *| in radians*
* `gal_lat`: Galactic latitude *| in radians*

# Arguments

* `$asc`: Right ascension *| in radians*
* `$dec`: Declination *| in radians*

The equatorial coordinates passed are assumed to be referred to the
standard equinox of B1950.0.
**/

macro_rules! GalFrmEq {
    ($asc: expr, $dec: expr) => {{
        (astro::coords::GalLongFrmEq($asc, $dec),
         astro::coords::GalLatFrmEq($asc, $dec))
    }};
}

//-------------------------------------------------------------------
// Galactic coordinates to equatorial coordinates

/**
Returns the **right ascension** from **galactic** coordinates

# Returns

* `asc`: Right ascension *| in radians*

The right ascension returned here is referred to the standard equinox
of  B1950.0.

# Arguments

* `gal_long`: Galactic longitude *| in radians*
* `gal_lat`: Galactic latitude *| in radians*
**/
pub fn AscFrmGal(gal_long: f64, gal_lat: f64) -> f64 {
      12.25_f64.to_radians()
    + (gal_long - 123_f64.to_radians()).sin()
      .atan2(   27.4_f64.to_radians().sin()*(gal_long - 123_f64.to_radians()).cos()
              - 27.4_f64.to_radians().cos()*gal_lat.tan()
            )
}

/**
Returns the **declination** from **galactic** coordinates

# Returns

* `dec`: Declination *| in radians*

The declination returned here is referred to the standard equinox
of  B1950.0.

# Arguments

* `gal_long`: Galactic longitude *| in radians*
* `gal_lat`: Galactic latitude *| in radians*
**/
pub fn DecFrmGal(gal_long: f64, gal_lat: f64) -> f64 {
    (   gal_lat.sin()*27.4_f64.to_radians().sin()
      + gal_lat.cos()*27.4_f64.to_radians().cos()*(gal_long - 123_f64.to_radians()).cos()
    ).asin()
}

/**
Returns **equatorial** coordinates from galactic** coordinates

# Returns

`(asc, dec)`

* `asc`: Right ascension *| in radians*
* `dec`: Declination *| in radians*

The equatorial coordinates returned here is referred to the standard
equinox of  B1950.0.

# Arguments

* `$gal_long`: Galactic longitude *| in radians*
* `$gal_lat`: Galactic latitude *| in radians*
**/

macro_rules! EqFrmGal {
    ($gal_long: expr, $gal_lat: expr) => {{
        (astro::coords::AscFrmGal($gal_long, $gal_lat),
         astro::coords::DecFrmGal($gal_long, $gal_lat))
    }};
}

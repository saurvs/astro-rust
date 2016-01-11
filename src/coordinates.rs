use angle;

/// Represents a point on the **Earth's surface**
pub struct GeographicalPoint {
    /// Right ascension
    pub long: f64,
    /// Declination
    pub lat: f64,
}

impl GeographicalPoint {
    pub fn AngularSeparation(&self, other_point: GeographicalPoint) -> f64 {
        angle::AngularSeparation(self.long, self.lat,
                                 other_point.long, other_point.lat)
    }
}

/// Represents a point on the **equatorial coordinate system**
pub struct EquatorialPoint {
    /// Right ascension
    pub right_ascen: f64,
    /// Declination
    pub declin: f64,
}

impl EquatorialPoint {
    pub fn AngularSeparation(&self, other_point: EquatorialPoint) -> f64 {
        angle::AngularSeparation(self.right_ascen, self.declin,
                                 other_point.right_ascen, other_point.declin)
    }
}

//-------------------------------------------------------------------

/**
Returns the **hour angle** from **local longitude** and **right ascension**

If right ascension is corrected for nutation, then sidereal time
too must be corrected for it.

# Returns

* ```hour_angle```: Hour angle *(radians)*

# Arguments

* ```green_sid```: Sidereal time at Greenwhich *(radians)*
* ```obv_long```: Observer's geographical longitude *(radians)*
* ```right_asc```: Right ascension *(radians)*
**/
pub fn HourAngleFromObserverLong(green_sid: f64, obv_long: f64, right_asc: f64) -> f64 {
    green_sid - obv_long - right_asc
}

/**
Returns the **hour angle** from **local sidereal time** and **right ascension**

If right ascension is corrected for nutation, then sidereal time
too must be corrected for it.

# Returns

* ```hour_angle```: Hour angle *(radians)*

# Arguments

* ```loc_sid```: Local sidereal time *(radians)*
* ```right_asc```: Right ascension *(radians)*
**/
pub fn HourAngleFromLocalSidereal(loc_sid: f64, right_asc: f64) -> f64 {
    loc_sid - right_asc
}

//-------------------------------------------------------------------
// Ecliptical coordinates to equatorial coordinates

/**
Returns the **ecliptical longitude** from **equatorial coordinates**

# Returns

* ```eclip_long```: Ecliptical longitude *(radians)*

# Arguments

* ```right_ascen```: Right ascension *(radians)*
* ```declin```: Declination *(radians)*
* ```oblq_eclip```: If ```right_asc``` and ```dec``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. *(radians)*
**/
pub fn EclipticalLongFromEquatorialCoords(right_ascen: f64, declin: f64, oblq_eclip: f64,) -> f64 {
    (   right_ascen.sin()  * oblq_eclip.cos()
      + declin.tan()       * oblq_eclip.sin()
    ).atan2(right_ascen.cos())
}

/**
Returns the **ecliptical latitude** from **equatorial coordinates**

# Returns

* ```ecliptical_latitude```: Ecliptical latitude *(radians)*

# Arguments

* ```right_ascen```: Right ascension *(radians)*
* ```declin```: Declination *(radians)*
* ```oblq_eclip```: If ```right_asc``` and ```dec``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. *(radians)*
**/
pub fn EclipticalLatFromEquatorialCoords(right_ascen: f64, declin: f64, oblq_eclip: f64) -> f64 {
    (   declin.sin() * oblq_eclip.cos()
      - declin.cos() * oblq_eclip.sin() * right_ascen.sin()
    ).asin()
}

/**
Returns **ecliptical coordinates** from **equatorial coordinates**

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
        (astro::coordinates::EclipticalLongFromEquatorialCoords($x, $y, $z),
         astro::coordinates::EclipticalLatFromEquatorialCoords($x, $y, $z))
    }};
}

//-------------------------------------------------------------------
// Equatorial coordinates to ecliptical coordinates

/**
Returns the **right ascension** from **ecliptical coordinates**

# Returns

* ```right_ascension```: Right ascension *(radians)*

# Arguments

* ```eclip_long```: Ecliptical longitude *(radians)*
* ```eclip_lat```: Ecliptical latitude *(radians)*
* ```oblq_eclip```: If ```ecl_long``` and ```ecl_lat``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. *(radians)*
**/
pub fn RightAscensionFromEclipticalCoords(eclip_long: f64, eclip_lat: f64, oblq_eclip: f64) -> f64 {
    (   eclip_long.sin() * oblq_eclip.cos()
      - eclip_lat.tan()  * oblq_eclip.sin()
    ).atan2(eclip_long.cos())
}

/**
Returns the **declination** from **ecliptical coordinates**

# Returns

* ```declination```: Declination *(radians)*

# Arguments

* ```eclip_long```: Ecliptical longitude *(radians)*
* ```eclip_lat```: Ecliptical latitude *(radians)*
* ```oblq_eclip```: If ```ecl_long``` and ```ecl_lat``` are corrected for
                  nutation, then *true* obliquity. If not, then
                  *mean* obliquity. *(radians)*
**/
pub fn DeclinationFromEclipticalCoords(eclip_long: f64, eclip_lat: f64, oblq_eclip: f64) -> f64 {
    (   eclip_lat.sin() * oblq_eclip.cos()
      + eclip_lat.cos() * oblq_eclip.sin() * eclip_long.sin()
    ).asin()
}

/**
Returns **equatorial coordinates** from **ecliptical coordinates**

# Returns

```(right_ascension, declination)```

* ```right_ascension```: Right ascension *(radians)*
* ```declination```: Declination *(radians)*

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
        (astro::coordinates::RightAscensionFromEclipticalCoords($x, $y, $z),
         astro::coordinates::DeclinationFromEclipticalCoords($x, $y, $z))
    }};
}

//-------------------------------------------------------------------
// Local horizontal coordinates from equatorial coordinates

/**
Returns the **azimuth**

# Returns

* ```azimuth```: Azimuth *(radians)*

# Arguments

* ```hour_angle```: Hour angle *(radians)*
* ```declin```: Declination *(radians)*
* ```observer_lat```: Observer's geographical latitude *(radians)*
**/
pub fn Azimuth(hour_angle: f64, declination: f64, observer_lat: f64) -> f64 {
    hour_angle.sin()
    .atan2(   hour_angle.cos()  * observer_lat.sin()
            - declination.tan() * observer_lat.cos()
          )
}

/**
Returns the **altitude**

# Returns

* ```altitude```: Altitude *(radians)*

# Arguments

* ```hour_angle```: Hour angle *(radians)*
* ```declin```: Declination *(radians)*
* ```observer_lat```: Observer's geographical latitude *(radians)*
**/
pub fn Altitude(hour_angle: f64, declination: f64, observer_lat: f64) -> f64 {
    (   observer_lat.sin() * declination.sin()
      + observer_lat.cos() * declination.cos() * hour_angle.cos()
    ).asin()
}

//-------------------------------------------------------------------
// Equatorial coordinates from local horizontal coordinates

/**
Returns the **hour angle** from **horizontal coordinates**

# Returns

* ```hour_angle```: Hour angle *(radians)*

# Arguments

* ```azimuth```: Azimuth *(radians)*
* ```altitude```: Altitude *(radians)*
* ```observer_lat```: Observer's geographical latitude *(radians)*
**/
pub fn HourAngleFromHorizontalCoords(azimuth: f64, altitude: f64, observer_lat: f64) -> f64 {
    azimuth.sin()
    .atan2(   azimuth.cos()  * observer_lat.sin()
            + altitude.tan() * observer_lat.cos()
          )
}

/**
Returns the **declination** from **horizontal coordinates**

# Returns

* ```declination```: Declination *(radians)*

# Arguments

* ```azimuth```: Azimuth *(radians)*
* ```altitude```: Altitude *(radians)*
* ```observer_lat```: Observer's geographical latitude *(radians)*
**/
pub fn DeclinationFromHorizontalCoords(azimuth: f64, altitude: f64, observer_lat: f64) -> f64 {
    (   observer_lat.sin() * altitude.sin()
      - observer_lat.cos() * azimuth.cos() * azimuth.cos()
    ).asin()
}

//-------------------------------------------------------------------
// Equatorial coordinates to galactic coordinates

/**
Returns the **galactic longitude** from **equatorial coordinates**

The equatorial coordinates must be referred to the standard equinox
of B1950.0

# Returns

* ```galactic_longitude```: Galactic longitude *(radians)*

# Arguments

* ```right_ascen```: Right ascension *(radians)*
* ```declin```: Declination *(radians)*
**/
pub fn GalacticLongFromEquatorialCoords(right_ascen: f64, declin: f64) -> f64 {
      303_f64.to_radians()
    - (192.25_f64.to_radians() - right_ascen).sin()
      .atan2(   27.4_f64.to_radians().sin() * (192.25_f64.to_radians() - right_ascen).cos()
              - 27.4_f64.to_radians().cos() * declin.tan()
            )
}

/**
Returns the **galactic latitude** from **equatorial coordinates**

The equatorial coordinates must be referred to the standard equinox
of B1950.0

# Returns

* ```galactic_latitude```: Galactic latitude *(radians)*

# Arguments

* ```right_ascen```: Right ascension *(radians)*
* ```declin```: Declination *(radians)*
**/
pub fn GalacticLatFromEquatorialCoords(right_ascen: f64, declin: f64) -> f64 {
    (   declin.sin() * 27.4_f64.to_radians().sin()
      + declin.cos() * 27.4_f64.to_radians().cos() * (192.25_f64.to_radians() - right_ascen).cos()
    ).asin()
}

/**
Returns **galactic coordinates** from **equatorial coordinates**

# Returns

```(galactic_longitude, galactic_latitude)```

* ```galactic_longitude```: Galactic longitude *(radians)*
* ```galactic_latitude```: Galactic latitude *(radians)*

# Arguments

* ```$x```: Right ascension *(radians)*
* ```$y```: Declination *(radians)*
**/
#[macro_export]
macro_rules! GalacticCoordsFromEquatorialCoords {
    ($x: expr, $y: expr) => {{
        (astro::coordinates::GalacticLongFromEquatorialCoords($x, $y),
         astro::coordinates::GalacticLatFromEquatorialCoords($x, $y))
    }};
}

//-------------------------------------------------------------------
// Galactic coordinates to equatorial coordinates

/**
Returns the **right ascension** from **galactic coordinates**

The right ascension is referred to the standard equinox
of B1950.0

# Returns

* ```right_ascension```: Right ascension *(radians)*

# Arguments

* ```gal_long```: Galactic longitude *(radians)*
* ```gal_lat```: Galactic latitude *(radians)*
**/
pub fn RightAscensionFromGalacticCoords(gal_long: f64, gal_lat: f64) -> f64 {
      12.25_f64.to_radians()
    + (gal_long - 123_f64.to_radians()).sin()
      .atan2(   27.4_f64.to_radians().sin() * (gal_long - 123_f64.to_radians()).cos()
              - 27.4_f64.to_radians().cos() * gal_lat.tan()
            )
}

/**
Returns the **declination** from **galactic coordinates**

The declination is referred to the standard equinox
of B1950.0

# Returns

* ```declination```: Declination *(radians)*

# Arguments

* ```gal_long```: Galactic longitude *(radians)*
* ```gal_lat```: Galactic latitude *(radians)*
**/
pub fn DeclinationFromGalacticCoords(gal_long: f64, gal_lat: f64) -> f64 {
    (   gal_lat.sin() * 27.4_f64.to_radians().sin()
      + gal_lat.cos() * 27.4_f64.to_radians().cos() * (gal_long - 123_f64.to_radians()).cos()
    ).asin()
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
        (astro::coordinates::RightAscensionFromGalacticCoords($x, $y),
         astro::coordinates::DeclinationFromGalacticCoords($x, $y))
    }};
}

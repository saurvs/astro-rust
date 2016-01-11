# astro-rust (work in progress)

[![License](https://img.shields.io/packagist/l/doctrine/orm.svg)](https://github.com/saurvs/astro-rust/blob/master/LICENSE.md) [![Build Status](https://travis-ci.org/saurvs/astro-rust.svg?branch=master)](https://travis-ci.org/saurvs/astro-rust)

Essential algorithms for fundamental astronomy.

* Find the positions of all the planets (and pluto), their moons (and the rings of Saturn), and the Sun.
* Calculate times of solar and lunar eclipses and passage through orbital nodes.
* Calculate detailed ephemeris for observations from the Earth.
* Correct for tiny effects in observation such as nutation,
   aberration, precession and atmospheric refraction.
* Convert points between different coordinate systems and epochs.
* And more related to asteroids, stars, binary stars, elliptic, parabolic and near-parabolic orbits, and satellites.

See full list of algorithms below.

Also, see [API Documentation](https://saurvs.github.io/astro-rust/) for this Cargo library.

## Usage

* Include the crate ```astro``` in your code
  ```rust
  extern crate astro;

  use astro::*;
  ```

* Find the Julian day (the most important step for almost everything)
  ```rust
  // time of the Apollo 11 moon landing

  let day_of_month = time::DayOfMonth{day: 20,
				 			          hr : 20, // UTC
                                      min: 18,
                                      sec: 4.0};

  let date = time::Date{year       : 1969,
                        month      : 7, // July
                        decimal_day: time::DecimalDay(day_of_month),
                        cal_type   : time::CalType::Gregorian};

  let julian_day = time::JulDay(&date);

  // to be super accurate, get the Julian Ephemeris day;
  // first find delta T, or get an observed value of delta T from 
  // the Astronomical Almanac
  
  let deltaT = time::ApproxDelT(date.year, date.month);

  let julian_ephm_day = time::JulEphmDay(&date, julian_day);

  ```

* Find the ecliptical *geocentric* coordinates of the Sun
  ```rust
  let (long, lat, rad_vec) = sun::EclGeocenCoords(julian_day);
  ```

* Also for the Moon
  ```rust
  let (long, lat, rad_vec) = planet::earth::moon::EclGeocenCoords(julian_day);
  ```

* Find the *heliocentric* coordinates of Saturn
  ```rust
  let (long, lat, rad_vec) = planet::HeliocenCoords(planet::Planet::Saturn, julian_day);
  ```

* Find the corrections for nutation (of the Earth) in ecliptical longitude and obliquity of the ecliptic
  ```rust
  let (nut_in_long, nut_in_oblq) = nutation::Corrections(julian_day);
  ```

* Find the geodesic distance between two locations on the Earth
  ```rust
	// geodesic distance between the Observatoire de Paris and
    // the US Naval Observatory at Washington DC

    let paris = coords::GeographPoint{long: angle::DegFrmDMS(-2, 20, 14.0).to_radians(),
                                      lat : angle::DegFrmDMS(48, 50, 11.0).to_radians()};

    let washington = coords::GeographPoint{long: angle::DegFrmDMS(77,  3, 56.0).to_radians(),
                                           lat : angle::DegFrmDMS(38, 55, 17.0).to_radians()};

	// angle::DegFrmDMS() converts degrees expressed in degrees,
	// minutes and seconds into degrees with decimals

    let distance = planet::earth::GeodesicDist(paris, washington); // in meters
  ```

* Convert equatorial coordinates to ecliptical coordinates
  ```rust
	// equatorial coordinates of the star Pollux

    let right_ascension = 116.328942_f64.to_radians();
    let declination = 28.026183_f64.to_radians();

    // mean obliquity of the ecliptic

    let oblq_eclip = 23.4392911_f64.to_radians();

    // you can also get oblq_eclip from ecliptic::MnOblq(julian_day)
    // for the same Julian day when the coordinates of the star
    // were observed

    // also make sure to type #[macro_use] before including the crate
    // to use macros

    // now, convert equatorial coordinates to ecliptical coordinates

    let (ecl_long, ecl_lat) = EclFrmEq!(right_ascension, declination, oblq_eclip);
  ```

* Convert equatorial coordinates to galactic coordinates
  ```rust
	// equatorial coordinates of the Nova Serpentis 1978

    let right_ascension = angle::DegFrmHMS(17, 48, 59.74).to_radians();
    let declination = angle::DegFrmDMS(-14, 43, 8.2).to_radians();

    // and to galactic coordinates

    let (gal_long, gal_lat) = GalFrmEq!(right_ascension, declination);
  ```

## Things you can find/do
* Heliocentric coordinates of Mercury, Venus, Earth, Mars, Jupiter, Saturn, Neptune, and Uranus (and Pluto).
* Geocentric coordinates of the Sun.
* Transformation of heliocentric coordinates to geocentric coordinates, and to topocentric coordinates.
* Correction in coordinates for nutation, precession, atmospheric refraction, and aberration.
* Transformation between the ecliptic, equatorial, horizontal, and galactic coordinates systems.
* Transformation of coordinates between different epochs and equinoxes.
* Distances between points on the Earth and obliquity of the Earth's ecliptic.
* Geocentric position, illuminated fraction and bright limb of the Moon.
* Semidiameters of Solar System bodies.
* Describing elliptic, parabolic, and near-parabolic orbits, and the motion of bodies in them.
* Orbital elements of all eight planets.
* Transformation of Gregorian and Julian dates to Julian day, and vice-versa.
* Delta T, Julian Ephemeris day and Julian century.
* Sidereal time.
* Equation of time.
* Stellar magnitudes.
* Describe binary stars.
* Asteroid diameters.

## References
* [Astronomical Algorithms, by Jean Meeus (2nd edition)](http://www.willbell.com/math/mc1.htm)
* [World Geodetic System 1984](https://confluence.qps.nl/pages/viewpage.action?pageId=29855173)
* [Five Millennium Canon of Solar Eclipses [Espenak and Meeus]](http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html)
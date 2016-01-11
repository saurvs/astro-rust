# astro-rust (work in progress)

[![License](https://img.shields.io/packagist/l/doctrine/orm.svg)](https://github.com/saurvs/astro-rust/blob/master/LICENSE.md) [![Build Status](https://travis-ci.org/saurvs/astro-rust.svg?branch=master)](https://travis-ci.org/saurvs/astro-rust)

Essential algorithms for fundamental astronomy.

* Find the positions of the planets (and pluto), their moons (and the rings of Saturn), and the Sun.
* Calculate times of solar and lunar eclipses and passage through orbital nodes.
* Calculate detailed ephemeris for observations from the Earth.
* Correct for tiny effects in observation such as nutation,
   aberration, precession and atmospheric refraction.
* Convert points between different coordinate systems and epochs.
* And more related to asteroids, stars, binary stars, elliptic, parabolic and near-parabolic orbits, satellites.

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

  let julian_day = time::JulDay(date);
  ```

* Find the ecliptical geocentric coordinates of the Sun
  ```rust
  let (long, lat, rad_vec) = sun::EclGeocenCoords(julian_day);
  ```

* And the Moon
  ```rust
  let (long, lat, rad_vec) = planet::earth::moon::EclGeocenCoords(julian_day);
  ```

* Find the heliocentric coordinates of Saturn
  ```rust
  let (long, lat, rad_vec) = planet::HeliocenCoords(planet::Planet::Mars, julian_day);
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
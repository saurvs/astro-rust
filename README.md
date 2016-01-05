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

* Add the dependency to your ```Cargo.toml```
   ```toml
   [dependencies]
   astro-rust = "1.0"
   ```

* Include the crate ```astro``` in your code
  ```rust
  extern crate astro;

  use astro::*;
  ```

* Find the Julian day (the most important step for almost everything)
  ```rust
  let day_of_month = time::DayOfMonth{ day: 17,
				 					   hour: 12,
                                       minute: 0, 
                                       second: 0.0 };

  let date = time::Date{ year: 2016,
                         month: 1,
                         day_of_month: day_of_month,
                         calendar_type: time::Gregorian };

  let JD = time::JulianDay(date);
  ```

* Find the equatorial coordinates of the Sun
* Find the topocentric coordinates (for viewing on the sky from Earth) of the Sun
* Similarly, for the Moon
* Find the heliocentric coordinates of Mars
  ```rust
  let JC = time::JulianCentury(JD);

  let (longitude, latitude, radius_vector) = planet::mars::HeliocentricCoords(JC);
  ```

* Find the topocentric coordinates (for viewing on the sky from Earth) of Mars
* Find the corrections for nutation
  ```rust
  let JED = time::JulianEphemerisDay(JD);

  let (nutation_in_longitude, nutation_in_obliquity) = nutation::Corrections(JC);
  ```

* Find the corrections for aberration
* Find the corrections for precession
* Find the corrections for atmospheric refraction
* Find the orbital elements of Jupiter
* Find the elements to describe the rings of Saturn
* Find the positions of the satellites of Saturn
* Find the time of solar and lunar eclipse close to a certain date

## Detailed list of things you can find/do
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
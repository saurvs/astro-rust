# astro-rust

[![License](https://img.shields.io/packagist/l/doctrine/orm.svg)](https://github.com/saurvs/astro-rust/blob/master/LICENSE.md) [![Build Status](https://travis-ci.org/saurvs/astro-rust.svg?branch=master)](https://travis-ci.org/saurvs/astro-rust)

```astro-rust``` is an MIT licensed library of algorithms focused on rigorous and accurate astronomical calculations.

Several such libraries have existed for a long time, being written in popular languages such as [C/C++](http://www.projectpluto.com/source.htm), [Python](https://pypi.python.org/pypi/astronomia/0.4.1) and [Java](http://celestjava.sourceforge.net/), which have comprehensive, well-tested coverage of various astronomical algorithms, albeit occasionally distributed with slightly restrictive licenses.

 This library however is written in the Rust programming language; a modern systems programming language which is fast, safe and expressive. It presents several improvements over other low level languages like C and C++, like excellent memory safety without a garbage collector, strong static typing, better concurrency support, better module system, and a [blazingly fast](http://benchmarksgame.alioth.debian.org/u64q/rust.html) runtime. This new language justifies a new astronomical library for use in the future.

 And, the MIT license adopted here is as liberal as open source licenses get, with the permission to do pretty much anything imaginable, as long as due credit is given to the original authors(s).

 Most of the algorithms implemented in this library are those described in the book *Astronomical Algorithms by Jean Meeus*, which includes things like planetary, solar and lunar positioning, corrections of precession, nutation, parallax, and aberration, times of conjunctions, elongations, and oppositions, physical ephemeris of Mars, Jupiter, Saturn and the Moon, finding position angles, illuminated fractions, and visual magnitudes. Even Pluto gets a chapter.

 However, the 2nd edition of the book was published in 1998, with only corrections for typos published since. And so, some of the algorithms (and physical constants) used in the book may differ from those used in this library, which were in turn adopted in recent times by NASA and the IAU.

For information related to the programming aspects of this library on the modules and functions available, see the [Rust API Documentation](https://saurvs.github.io/astro-rust/).

A high level list of algorithms is given at the bottom of this page.

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
                        decimal_day: time::DecimalDay(&day_of_month),
                        cal_type   : time::CalType::Gregorian};

  let julian_day = time::JulDay(&date);

  // to be super accurate, get the Julian Ephemeris day;
  // first find delta T, or get an observed value of delta T from
  // the Astronomical Almanac

  let deltaT = time::ApproxDelT(date.year, date.month);

  let julian_ephm_day = time::JulEphmDay(julian_day, deltaT);

  ```

* Find the ecliptic *geocentric* coordinates of the Sun
  ```rust
  let (long, lat, rad_vec) = sun::EclGeocenCoords(julian_day);

  // long    - ecliptic longitude (radians)
  // lat     - ecliptic latitude (radians)
  // rad_vec - distance between the Sun and the Earth (AU)
  ```

* Also for the Moon
  ```rust
  let (long, lat, rad_vec) = planet::earth::moon::EclGeocenCoords(julian_day);
  ```

* Find the *heliocentric* coordinates of Jupiter
  ```rust
  let (long, lat, rad_vec) = planet::HeliocenCoords(planet::Planet::Jupiter, julian_day);
  ```

* Find the corrections for nutation in ecliptic longitude and obliquity of the ecliptic
  ```rust
  let (nut_in_long, nut_in_oblq) = nutation::Nutation(julian_day);
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

    let distance = planet::earth::GeodesicDist(&paris, &washington); // in meters
  ```

* Convert equatorial coordinates to ecliptic coordinates
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

    // now, convert equatorial coordinates to ecliptic coordinates

    let (ecl_long, ecl_lat) = EclFrmEq!(right_ascension, declination, oblq_eclip);
  ```

* Convert equatorial coordinates to galactic coordinates
  ```rust
	// equatorial coordinates of the Nova Serpentis 1978

    let right_ascension = angle::DegFrmHMS(17, 48, 59.74).to_radians();
    let declination = angle::DegFrmDMS(-14, 43, 8.2).to_radians();

    // convert to galactic coordinates

    let (gal_long, gal_lat) = GalFrmEq!(right_ascension, declination);
  ```

## Algorithms

Algorithms implemented in this library allow you to calculate or perform the following:

The 8 Solar System Planets

* heliocentric coordinates (using the *full* VSOP87-D solution)
* orbital elements
* apparent magnitude
* equatorial semidiameter
* polar semidiameter for Saturn and Jupiter
* time of passage through the nodes
* illuminated fraction of the planetary disk

Sun

* ecliptic geocentric coordinates
* rectangular geocentric coordinates
* ephemeris for physical observations
* time of the begining of Carrington's synodic rotation

Moon

* ecliptic geocentric coordinates
* optical, physical and topocentric liberations
* times of passage through nodes
* illuminated fraction of the lunar disk
* equatorial horizontal parallax

Earth

* high accuracy geodesic
* eccentricity of the meridian
* angle between the diurnal path and the horizon

Pluto

* heliocentric coordinates
* orbital elements
* apparent magnitude
* equatorial semidiameter

Mars

* ecliptic coordinates of the North pole
* ephemeris for physical observations

Jupiter

* ephemeris for physical observations

Transform

* ecliptic coordinates to equatorial coordinates, and vice-versa
* equatorial coordinates to local horizontal coordinates, and vice-versa
* ecliptic coordinates to galactic coordinates, and vice-versa
* heliocentric coordinates to ecliptic/equatorial geocentric coordinates
* orbital elements to heliocentric coordinates
* ecliptic/equatorial coordinates from one epoch to another
  (due to precession)
* orbital elements from one equinox to another

Time

* Julian day from Gregorian and Julian dates, and vice-versa
* a good analytic approximation of delta T over several centuries
* mean and apparent sidereal time

Ecliptic

* mean obliquity (IAU and Laskar formulae)
* angle between the ecliptic and the horizon
* longitudes of the two ecliptic points on the horizon

Nutation in

* ecliptic longitude
* obliquity of the ecliptic
* right ascension and declination

Atmospheric refraction

* apparent altitude from true altitude, and vice-versa
* effect of local pressure and temperature

Stars

* combined magnitude of two or more stars
* absolute magnitude from parallax or distance from Earth
* brightness ratio from difference in magnitude, and vice-versa
* angle between a vector from a star to the Earth's north
  celestial pole and a vector from the same star to the north
  pole of the ecliptic
* equatorial coordinates at a different time due to proper
  motion radial velocity

Binary stars

* angular separation
* apparent position angle
* eccentric of the apparent orbit
* radius vector
* true anomaly
* mean annual motion of the companion star from period of revolution
* mean anomaly of the companion star

Asteroids

* true diameter from absolute magnitude and albedo
* apparent diameter from true diameter and distance to the
  Earth

## Contributing

The goal of this project is to build a polished, production-quality library with a well-tested implementation of algorithms useful in astronomy, with a well-designed Rust API and an excellent documentation of code and coverage.

Anyone interested to contribute in any way possible is encouraged to do so.

A good start would be to go through Meeus's book or documents published by the IAU, browse the API documentation, read through the code, and submit a pull request for a new algorithm or modification of an existing one. Refactored code and minor optimizations are also accepted.

One suggested addition is the recent IAU 2000/2006 precession-nutation model. This method improves upon the existing model implemented here so far by taking into
account the effect of mantle anelasticity, ocean tides, electromagnetic couplings
produced between the fluid outer core and the mantle as well as between the solid
inner core and fluid outer core

## References
* [Astronomical Algorithms, by Jean Meeus (2nd edition)](http://www.willbell.com/math/mc1.htm)
* [World Geodetic System 1984](https://confluence.qps.nl/pages/viewpage.action?pageId=29855173)
* [Five Millennium Canon of Solar Eclipses [Espenak and Meeus]](http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html)

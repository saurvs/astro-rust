# astro-rust

[![License](https://img.shields.io/packagist/l/doctrine/orm.svg)](https://github.com/saurvs/astro-rust/blob/master/LICENSE.md) [![Build Status](https://travis-ci.org/saurvs/astro-rust.svg?branch=master)](https://travis-ci.org/saurvs/astro-rust)

```astro-rust``` implements, in the rust programming language, essential algorithms used in astronomy for accurate calculations.

The main reference used for implementation is the comprehensive set of algorithms described in book *Astronomical Algorithms by Jean Meeus*, that includes methods for calculating the positions of the planets, the sun, the moon, converting between different celestial coordinates systems, finding corrections for nutation and atmospheric refraction, and so on.

Detailed information about functions and modules available in this library is in the [Rust API Documentation](https://saurvs.github.io/astro-rust/), and the full listing of algorithms available is given below.

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

Algorithms implemented in this library so far allow you to calculate or perform the following:

* **Time** Julian day from Gregorian and Julian dates, and vice-versa. Analytic approximation to delta T. Mean and apparent sidereal time.
* **Transform** coordinates between equatorial, ecliptic, topocentric, local horizontal, and galactic coordinate systems
* **Transform** heliocentric coordinates to ecliptic and equatorial geocentric coordinates
* **All 8 Planets** heliocentric coordinates (full VSOP87-D solution), orbital elements, apparent magnitudes, equatorial semidiameters, times of passage through nodes, illuminated fraction of the disk
* **Pluto** heliocentric coordinates, apparent magnitude, equatorial semidiameter
* **Earth** high-accuracy geodesic distances
* **Moon** ecliptic geocentric coordinates, optical, physical and topocentric librations, times of passage through nodes,
illuminated fraction of the disk
* **Sun** ecliptic geocentric coordinates, rectangular geocentric coordinates, ephemeris for physical observations, times of Carrington's synodic rotation
* **Mars** ecliptic coordinates of the North pole, ephemeris for physical observations
* **Stars** combined magnitude of two or more stars, absolute magnitude from parallax
* **Binary stars** radius vector, angular separation, eccentricity of the apparent orbit, apparent position angle, and mean annual motion and mean anomaly of the companion star
* **Asteroids** true and apparent diameters
* **Atmospheric refraction** apparent altitude from true altitude and vice-versa, effect of local pressure and temperature
* **Ecliptic** mean obliquity by IAU and Laskar
* **Nutation** nutation in ecliptic longitude, obliquity, and equatorial coordinates

## Contributing

The goal of this project is to build

* a polished, production-quality library with a well-tested implementation of algorithms used in fundamental astronomy,
* a well-designed Rust API,
* and an excellent documentation of code and coverage.

Anyone interested to contribute in any way possible is encouraged to do so.

A good start would be to go through Meeus's book and browse the API documentation, read the code, and submit a pull request if anything warrants an addition or a modification.

Even refactoring code and minor optimizations are accepted.

## References
* [Astronomical Algorithms, by Jean Meeus (2nd edition)](http://www.willbell.com/math/mc1.htm)
* [World Geodetic System 1984](https://confluence.qps.nl/pages/viewpage.action?pageId=29855173)
* [Five Millennium Canon of Solar Eclipses [Espenak and Meeus]](http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html)

# astro-rust [![](http://meritbadge.herokuapp.com/astro)](https://crates.io/crates/astro) [![](https://travis-ci.org/saurvs/astro-rust.svg?branch=master)](https://travis-ci.org/saurvs/astro-rust) [![](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/saurvs/astro-rust/blob/master/LICENSE.md)

**Contents**

* [About](#about)
* [Usage](#usage)
* [Contributing](#contributing)
* [References](#references)

[API Docs](https://saurvs.github.io/astro-rust/)

> This project is currently a work in progress.
> Some APIs may be updated occasionally.

## About

```astro-rust``` is an MIT licensed library of astronomical algorithms for the Rust programming language.

Implemented capabilities include planetary, solar, and lunar positioning, corrections for precession, nutation, parallax, and aberration, calculation of physical ephemeris of Mars, Jupiter, and the ring system of Saturn, finding times of rise, set and transit of celestial bodies and [much more](https://saurvs.github.io/astro-rust/).

## Usage

* Add the dependency ```astro``` in your ```Cargo.toml```
  ```toml
  [dependencies]
  astro = "1.0.4"
  ```

* Include the crate ```astro``` in your code
  ```rust
  extern crate astro;

  use astro::*;
  ```

* Specify your time of interest using the [Julian day](http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html)
  ```rust
  // for example, the time of the Apollo 11 moon landing

  let day_of_month = time::DayOfMonth{day      : 20,  // 1 - 31
				 			          hr       : 20,  // 0 - 24
                                      min      : 18,  // 0 - 60
                                      sec      : 4.0, // 0.0 - 60.0
                                      time_zone: 0.0  // at the Greenwhich meridian
                                     };

  let date = time::Date{year       : 1969, // +ve for AD, -ve for BC
                        month      : 7, // July
                        decimal_day: time::decimal_day(&day_of_month),
                        cal_type   : time::CalType::Gregorian};

  let julian_day = time::julian_day(&date);

  // for higher accuracy in specifying the time of interest,
  // find the Julian Ephemeris day; this slightly differs from
  // the Julian day by usually a few seconds, denoted as ΔT.
  // you can get a reported value of that from the Astronomical
  // Almanac, or calculate it using a built-in function

  let delta_t = time::approx_delta_t(date.year, date.month); // a good one, actually

  let julian_ephm_day = time::julian_emph_day(julian_day, delta_t);
  ```

* Find the position of the Sun and the Moon with respect to the Earth
  ```rust

  // geocentric ecliptic point and radius vector of the Sun
  let (sun_ecl_point, rad_vec_sun) = sun::geocen_ecl_pos(julian_day);

  // sun_ecl_point.long    - ecliptic longitude (radians)
  // sun_ecl_point.lat     - ecliptic latitude  (radians)
  // rad_vec_sun - distance between the Sun and the Earth (AU)

  // and similarly for the Moon
  let (moon_ecl_point, rad_vec_moon) = lunar::geocen_ecl_pos(julian_day);

  ```

* Find the position of a planet with respect to the Sun
  ```rust
  // the heliocentric point and radius vector of a planet, like Jupiter
  let (jup_long, jup_lat, rad_vec) = planet::heliocen_pos(&planet::Planet::Jupiter, julian_day);

  // or neptune
  let (nep_long, nep_lat, rad_vec) = planet::heliocen_pos(&planet::Planet::Neptune, julian_day);

  // all eight planets (and (the dwarf planet) Pluto) are supported

  let (plut_long, plut_lat, rad_vec) = pluto::heliocen_pos(julian_day);
  ```

* Find the geodesic distance between two locations on Earth
  ```rust
	// geodesic distance between the Observatoire de Paris and
    // the US Naval Observatory at Washington DC

    let paris = coords::GeographPoint{long: angle::deg_frm_dms(-2, 20, 14.0).to_radians(),
                                      lat : angle::deg_frm_dms(48, 50, 11.0).to_radians()};

    let washington = coords::GeographPoint{long: angle::deg_frm_dms(77,  3, 56.0).to_radians(),
                                           lat : angle::deg_frm_dms(38, 55, 17.0).to_radians()};

	// angle::deg_frm_dms() converts degrees expressed in degrees,
	// minutes and seconds into a fractional degree

    let distance = planet::earth::geodesic_dist(&paris, &washington); // in meters
  ```

* Convert equatorial coordinates to ecliptic coordinates
  ```rust
	// equatorial coordinates of the star Pollux

    let right_ascension = 116.328942_f64.to_radians();
    let declination = 28.026183_f64.to_radians();

    // mean obliquity of the ecliptic

    let oblq_eclip = 23.4392911_f64.to_radians();

    // you can also get oblq_eclip from ecliptic::mn_oblq_IAU(julian_day)
    // for the Julian day on which the coordinates of the star
    // were observed

    // also make sure to type #[macro_use] before including the crate
    // to use macros

    // now, convert equatorial coordinates to ecliptic coordinates

    let (ecl_long, ecl_lat) = ecl_frm_eq!(right_ascension, declination, oblq_eclip);
  ```

* Convert equatorial coordinates to galactic coordinates
  ```rust
	// equatorial coordinates of the Nova Serpentis 1978

    let right_ascension = angle::deg_frm_hms(17, 48, 59.74).to_radians();
    let declination = angle::deg_frm_dms(-14, 43, 8.2).to_radians();

    // convert to galactic coordinates

    let (gal_long, gal_lat) = gal_frm_eq!(right_ascension, declination);
  ```

* Correct for nutation in different coordinate systems
  ```rust
  // nutation in ecliptic longitude and obliquity of the ecliptic
  let (nut_in_long, nut_in_oblq) = nutation::nutation(julian_day);

  // nutation in equatorial coordinates
  let (nut_in_asc, nut_in_dec) = nutation::nutation_in_eq_coords(julian_day);
  ```

## Contributing

Anyone interested to contribute in any way possible is encouraged to do so. Not all the algorithms in Meeus's book have been implemented yet. Documentation and tests need to be written for them as well. Refactored code and minor optimizations for the existing code are also welcome.

The end goal (of this project) is to build a modern, well-tested, well-documented library of algorithms for future use in astronomy. And doing it all in Rust is very much the right choice for that.

A fun suggestion is the addition of the recent [IAU 2000/2006 precession-nutation model](http://62.161.69.131/iers/conv2010/conv2010_c5.html). This method improves upon the existing model implemented here *"by taking into account the effect of mantle anelasticity, ocean tides, electromagnetic couplings produced between the fluid outer core and the mantle as well as between the solid inner core and fluid outer core"*.

## References


The main reference used as the source of algorithms is the famous book *Astronomical Algorithms by Jean Meeus*, whose almost every chapter has been addressed here, with functions that are well-documented and tests that use example data from the book; in some cases, such as ΔT approximation and planetary heliocentric positioning, more accurate methods have been implemented.

* Most algorithms: [Astronomical Algorithms, 2nd edition (Meeus)](http://www.willbell.com/math/mc1.htm)
* Planetary heliocentric positioning: [VSOP87-D](http://cdsarc.u-strasbg.fr/viz-bin/qcat?VI/81/)
* Approximating ΔT: [Five Millennium Canon of Solar Eclipses (Espenak and Meeus)](http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html)
* Some physical constants: [World Geodetic System 1984](https://confluence.qps.nl/pages/viewpage.action?pageId=29855173)

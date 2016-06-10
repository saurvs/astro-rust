/*
Copyright (c) 2015, 2016 Saurav Sachidanand

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

//! Stars

/**
Computes the combined magnitude of two stars

# Arguments

* `m1`: Magnitude of star 1
* `m2`: Magnitude of star 2
**/
#[inline]
pub fn combined_mag(m1: f64, m2: f64) -> f64
{
    m2 - 2.5 * (brightness_ratio(m1, m2) + 1.0)
}

/**
Computes the combined magnitude of two or more stars

# Arguments

* `m`: Array of magnitudes of stars
**/
pub fn combined_mag_of_many(m: &[f64]) -> f64
{
    let mut sum = 0.0;

    for i in m.iter() {
        sum += 10_f64.powf(-0.4 * i);
    }

    -2.5 * sum.log10()
}

/**
Computes the brightness ratio of two stars

# Arguments

* `m1`: Magnitude of star 1
* `m2`: Magnitude of star 2
**/
#[inline]
pub fn brightness_ratio(m1: f64, m2: f64) -> f64
{
    10.0_f64.powf(0.4 * (m2 - m1))
}

/**
Computes the difference in magnitude of two stars

# Arguments

* `br`: Brightness ratio of two stars
**/
#[inline]
pub fn mag_diff(br: f64) -> f64
{
    2.5 * br.log10()
}

/**
Computes the absolute magnitude of a star from its parallax

# Arguments

* `par`: Parallax of the star
* `am`: Apparent magnitude of the star
**/
#[inline]
pub fn abs_mag_frm_parallax(par: f64, am: f64) -> f64
{
    am + 5.0 + 5.0*(par.to_degrees() * 3600.0).log10()
}

/**
Computes the absolute magnitude of a star from its distance from earth

# Arguments

* `d`: The star's to earth *(parsecs)*
* `am`: Apparent magnitude of the star
**/

#[inline]
pub fn abs_mag_frm_dist(d: f64, am: f64) -> f64
{
    am + 5.0 - 5.0*d.log10()
}

/**
Computes the angle between a vector from a star to the
north celestial pole of the Earth and a vector from the
same star to the north pole of the ecliptic

# Returns

* `angle`: The desired angle *| in radians*

# Arguments

* `eclip_long`: The star's ecliptical longitude *| in radians*
* `eclip_lat`: The star's ecliptical latitude *| in radians*
* `oblq_eclip`: Obliquity of the ecliptic *| in radians*
**/

#[inline]
pub fn angl_between_north_celes_and_eclip_pole(eclip_long: f64,
                                               eclip_lat: f64,
                                               oblq_eclip: f64) -> f64
{
    (eclip_long.cos() * oblq_eclip.tan()).atan2 (
        eclip_lat.sin() * eclip_long.sin() * oblq_eclip.tan()
      - eclip_lat.cos()
    )
}

/**
Computes the equatorial coordinates of a star at
at a different time from it's motion in space

This function Computes the equatorial coordinates
of a star at a different time by taking into account
it's proper motion, distance and radial velocity.

# Returns

`(new_asc, new_dec)`

* `new_asc`: Right ascension at the different
                 time *| in radians*
* `new_dec`: Declination at the different
                 time *| in radians*

# Arguments

* `asc0`: Right ascension of the star initially *| in radians*
* `dec0`: Declination of the star initially *| in radians*
* `r`: Distance of the star (*parsecs*)
* `delta_r`: Radial velocity of the star (*parsecs/second*)
* `proper_motion_asc`: Proper motion of the star in right ascension
                           *| in radians*
* `proper_motion_dec`: Proper motion of the star in declination
                           *| in radians*
* `t`: Decimal years from the inital time; negative in the past
          and positive in the future
**/
pub fn eq_coords_frm_motion(asc0: f64,
                            dec0: f64,
                            r: f64,
                            delta_r: f64,
                            proper_motion_asc: f64,
                            proper_motion_dec: f64,
                            t: f64) -> (f64, f64)
{
    let x = r * dec0.cos() * asc0.cos();
    let y = r * dec0.cos() * asc0.sin();
    let z = r * dec0.sin();

    let delta_asc = 3600.0 * proper_motion_asc.to_degrees()/13751.0;
    let delta_dec = 3600.0 * proper_motion_dec.to_degrees()/206265.0;

    let delta_x = (x / r)*delta_r - z*delta_dec*asc0.cos() - y*delta_asc;
    let delta_y = (y / r)*delta_r - z*delta_dec*asc0.sin() + x*delta_asc;
    let delta_z = (z / r)*delta_r + r*delta_dec*dec0.cos();

    let x1 = x + t*delta_x;
    let y1 = y + t*delta_y;
    let z1 = z + t*delta_z;

    let asc = y1.atan2(x1);
    let dec = z1.atan2((x1*x1 + y1*y1).sqrt());

    (asc, dec)
}

pub fn proper_motion_in_eq_coords(asc: f64,
                                  dec: f64,
                                  pmotion_asc: f64,
                                  pmotion_dec: f64,
                                  ecl_lat: f64,
                                  oblq_eclip: f64) -> (f64, f64)
{
    let ecl_lat_cos = ecl_lat.cos();

    let pmotion_long = (
        pmotion_dec * oblq_eclip.sin() * asc.cos()
      + pmotion_asc * dec.cos() * (
            oblq_eclip.cos() * dec.cos()
          + oblq_eclip.sin() * dec.sin() * asc.sin()
        )
    ) / (ecl_lat_cos * ecl_lat_cos);

    let pmotion_lat = (
        pmotion_dec * (
            oblq_eclip.cos() * dec.cos()
          + oblq_eclip.sin() * dec.sin() * asc.sin()
        )
      - pmotion_asc * oblq_eclip.sin() * asc.cos() * dec.cos()
    ) / ecl_lat_cos;

    (pmotion_long, pmotion_lat)
}

use angle;

/**
Returns the combined magnitude of two stars

# Arguments

* ```m1```: Magnitude of star 1
* ```m2```: Magnitude of star 2
**/
pub fn CombinedMag(m1: f64, m2: f64) -> f64 {
    m2 - 2.5*(BrightnessRatio(m1, m2) + 1.0)
}

/**
Returns the combined magnitude of two or more stars

# Arguments

* ```m```: Array of magnitudes of stars
**/
pub fn CombinedMagOfMany(m: &[f64]) -> f64 {
    let mut sum = 0.0;
    for i in m.iter() {
        sum += 10_f64.powf(-0.4 * i);
    }
    -2.5 * sum.log10()
}

/**
Returns the brightness ratio of two stars

# Arguments

* ```m1```: Magnitude of star 1
* ```m2```: Magnitude of star 2
**/
pub fn BrightnessRatio(m1: f64, m2: f64) -> f64 {
    10.0_f64.powf(0.4 * (m2-m1))
}

/**
Returns the difference in magnitude of two stars

# Arguments

* ```br```: Brightness ratio of two stars
**/
pub fn MagDiff(br: f64) -> f64 {
    2.5 * br.log10()
}

/**
Returns the absolute magnitude of a star from its parallax

# Arguments

* ```par```: Parallax of the star
* ```am```: Apparent magnitude of the star
**/
pub fn AbsMagFrmParallax(mut par: f64, am: f64) -> f64 {
    par = par.to_degrees() * 3600.0;
    am + 5.0 + 5.0*par.log10()
}

/**
Returns the absolute magnitude of a star from its distance from earth

# Arguments

* ```d```: The star's to earth *(parsecs)*
* ```am```: Apparent magnitude of the star
**/
pub fn AbsMagFrmDist(d: f64, am: f64) -> f64 {
    am + 5.0 - 5.0*d.log10()
}

/**
Returns the angle between a vector from a star to the
Earth's north celestial pole and a vector from the
same star to the north pole of the ecliptic

# Returns

* ```angle```: The desired angle *| in radians*

# Arguments

* ```eclip_long```: The star's ecliptical longitude *| in radians*
* ```eclip_lat```: The star's ecliptical latitude *| in radians*
* ```oblq_eclip```: Obliquity of the ecliptic *| in radians*
**/
pub fn AnglBetweenNorthCelesAndEclipticPole(eclip_long: f64, eclip_lat: f64, oblq_eclip: f64) -> f64 {
    (eclip_long.cos() * oblq_eclip.tan())
    .atan2(   eclip_lat.sin() * eclip_long.sin() * oblq_eclip.tan()
            - eclip_lat.cos())
}

/**
Returns the equatorial coordinates of a star at
at a different time from it's motion in space

This function returns the equatorial coordinates
of a star at a different time by taking into account
it's proper motion, distance and radial velocity.

# Returns

```(new_asc, new_dec)```

* ```new_asc```: Right ascension at the different
                 time *| in radians*
* ```new_dec```: Declination at the different
                 time *| in radians*

# Arguments

* ```asc0```: Right ascension of the star initially *| in radians*
* ```dec0```: Declination of the star initially *| in radians*
* ```r```: Distance of the star (*parsecs*)
* ```delta_r```: Radial velocity of the star (*parsecs/second*)
* ```proper_motion_asc```: Proper motion of the star in right ascension
                           *| in radians*
* ```proper_motion_dec```: Proper motion of the star in declination
                           *| in radians*
* ```t```: Decimal years from the inital time; negative in the past
          and positive in the future
**/
pub fn EqCoordsFrmMotion(asc0: f64, dec0: f64, r: f64, delta_r: f64,
                     proper_motion_asc: f64, proper_motion_dec: f64, t: f64) -> (f64, f64) {

    let x = r*dec0.cos()*asc0.cos();
    let y = r*dec0.cos()*asc0.sin();
    let z = r*dec0.sin();

    let delta_asc = angle::TimeSecFrmDeg(proper_motion_asc.to_degrees()) / 13751.0;
    let delta_dec = angle::ArcSecFrmDeg(proper_motion_dec.to_degrees()) / 206265.0;

    let delta_x = (x/r)*delta_r - z*delta_dec*asc0.cos() - y*delta_asc;
    let delta_y = (y/r)*delta_r - z*delta_dec*asc0.sin() + x*delta_asc;
    let delta_z = (z/r)*delta_r + r*delta_dec*dec0.cos();

    let x1 = x + t*delta_x;
    let y1 = y + t*delta_y;
    let z1 = z + t*delta_z;

    let asc = y1.atan2(x1);
    let dec = z1.atan2((x1*x1 + y1*y1).sqrt());

    (asc, dec)

}

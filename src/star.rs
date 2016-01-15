use time;

/**
Returns the **combined magnitude** of two stars

# Arguments

* ```m1```: Magnitude of star 1
* ```m2```: Magnitude of star 2
**/
pub fn CombinedMag(m1: f64, m2: f64) -> f64 {
    m2 - 2.5*(BrightnessRatio(m1, m2) + 1.0)
}

/**
Returns the **combined magnitude** of two or more stars

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
Returns the **brightness ratio** of two stars

# Arguments

* ```m1```: Magnitude of star 1
* ```m2```: Magnitude of star 2
**/
pub fn BrightnessRatio(m1: f64, m2: f64) -> f64 {
    10.0_f64.powf(0.4 * (m2-m1))
}

/**
Returns the **difference in magnitude** of two stars

# Arguments

* ```br```: Brightness ratio of two stars
**/
pub fn MagDiff(br: f64) -> f64 {
    2.5 * br.log10()
}

/**
Returns the **absolute magnitude** of a star from its parallax

# Arguments

* ```par```: Parallax of the star
* ```am```: Apparent magnitude of the star
**/
pub fn AbsMagFrmParallax(mut par: f64, am: f64) -> f64 {
    par = par.to_degrees() * 3600.0;
    am + 5.0 + 5.0*par.log10()
}

/**
Returns the **absolute magnitude** of a star from its distance from earth

# Arguments

* ```d```: The star's to earth *(parsecs)*
* ```am```: Apparent magnitude of the star
**/
pub fn AbsMagFrmDist(d: f64, am: f64) -> f64 {
    am + 5.0 - 5.0*d.log10()
}

/**
Returns the **angle** between the vector from a star to the
**Earth's northern celestial pole** and vector from the same star to the
**north pole** of the **ecliptic**

# Returns

* ```angle```: The desired angle *(radians)*

# Arguments

* ```eclip_long```: The star's ecliptical longitude *(radians)*
* ```eclip_lat```: The star's ecliptical latitude *(radians)*
* ```oblq_eclip```: Obliquity of the ecliptic *(radians)*
**/
pub fn AnglBetweenNorthCelesAndEclipticPole(eclip_long: f64, eclip_lat: f64, oblq_eclip: f64) -> f64 {
    (eclip_long.cos() * oblq_eclip.tan())
    .atan2(   eclip_lat.sin() * eclip_long.sin() * oblq_eclip.tan()
            - eclip_lat.cos())
}

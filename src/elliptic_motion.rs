use std::*;

/**
Computes instantaneous **velocity** *(meters per second)* of a body in an
unperturbed elliptic orbit

# Arguments

* ```dist_to_sun```: Body's distance to Sun *(AU)*
* ```semimaj_axis```: Semimajor axis of orbit *(AU)*
**/
pub fn Velocity(dist_to_sun: f64, semimaj_axis:f64) -> f64 {
    0.0421219 * (1.0/dist_to_sun - 1.0/(2.0 * semimaj_axis)).sqrt()
}

/**
Computes **velocity** *(meters per second)* of a body at **perihelion**
in an unperturbed elliptic orbit

# Arguments

* ```semimaj_axis```: Semimajor axis of orbit *(AU)*
* ```orb_eccen```: Eccentricity of orbit
**/
pub fn PerihelionVelocity(semimaj_axis:f64, orb_eccen:f64) -> f64 {
    0.0297847 * ((1.0 + orb_eccen) / ((1.0 - orb_eccen) * semimaj_axis)).sqrt()
}

/**
Computes **velocity** *(meters per second)* of a body at **aphelion**
in an unperturbed elliptic orbit

# Arguments

* ```semimaj_axis```: Semimajor axis of orbit *(AU)*
* ```orb_eccen```: Eccentricity of orbit
**/
pub fn AphelionVelocity(semimaj_axis:f64, orb_eccen:f64) -> f64 {
    0.0297847 * ((1.0 - orb_eccen) / ((1.0 + orb_eccen) * semimaj_axis)).sqrt()
}

/**
Computes the approximate **length** of an ellipse using the Ramanujan
method

# Return values

This function returns an approximate value for the length of an
ellipse, using a formula given by [Ramanujan](https://en.wikipedia.org/wiki/Srinivasa_Ramanujan) in 1914.

The **error** is:

* 0% for a = b
* 0.4155% for e = 1

# Arguments

* ```a```: Semimajor axis of the ellipse (same unit as ```b```)
* ```b```: Semiminor axis of the ellipse (same unit as ```a```)
* ```e```: Eccentricity of the ellipse
**/
pub fn LengthOfEllipse_Ramanujan(a: f64, b: f64, e: f64) -> f64 {
    f64::consts::PI * (3.0*(a + b) - ((a + 3.0*b)*(3.0*a + b)).sqrt())
}

/**
Computes the approximate **length** of an ellipse

# Return values

This function returns an approximate value for the length of an
ellipse.

The **error** is:

* less than 0.001% if e < 0.88
* less than 0.01% if e < 0.95
* 1% if e = 0.9997
* 3% if e = 1

# Arguments

* ```a```: Semimajor axis of the ellipse (same unit as ```b```)
* ```b```: Semiminor axis of the ellipse (same unit as ```a```)
* ```e```: Eccentricity of the ellipse
**/
pub fn LengthOfEllipse(a: f64, b: f64, e: f64) -> f64 {
    let A = (a + b)/2.0;
    let G = (a*b).sqrt();
    let H = (2.0 * a * b)/(a + b);
    f64::consts::PI * (21.0*A - 2.0*G - 3.0*H) / 8.0
}

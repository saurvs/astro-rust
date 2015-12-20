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

* ```a```: Semimajor axis of the ellipse (same unit as that of ```b```)
* ```b```: Semiminor axis of the ellipse (same unit as that of ```a```)
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

* ```a```: Semimajor axis of the ellipse (same unit as that of ```b```)
* ```b```: Semiminor axis of the ellipse (same unit as that of ```a```)
* ```e```: Eccentricity of the ellipse
**/
pub fn LengthOfEllipse(a: f64, b: f64, e: f64) -> f64 {
    let A = (a + b)/2.0;
    let G = (a*b).sqrt();
    let H = (2.0 * a * b)/(a + b);
    f64::consts::PI * (21.0*A - 2.0*G - 3.0*H) / 8.0
}

/**
Computes the **semimajor axis**

# Arguments

* ```perih```: Perihelion of the orbit
* ```ecc```: Eccentricity of the orbit
**/
pub fn SemimajorAxis(perih: f64, ecc: f64) -> f64 {
    perih / (1.0 - ecc)
}

/**
Computes the **mean motion** *(radians/day)*

# Arguments

* ```semimaj_ax```: Semimajor axis of the orbit
**/
pub fn MeanMotion(semimaj_ax: f64) -> f64 {
    0.01720209895 / (semimaj_ax.powf(1.5))
}

pub fn TimeOfPassageThroughAscendingNode(w: f64, n: f64, a: f64, e: f64, T: f64) -> (f64, f64) {
    time_of_passage_through_node(-1.0 * w, n, a, e, T)
}

pub fn TimeOfPassageThroughDescendingNode(w: f64, n: f64, a: f64, e: f64, T: f64) -> (f64, f64) {
    time_of_passage_through_node(180_f64.to_radians() * w, n, a, e, T)
}

fn time_of_passage_through_node(v: f64, n: f64, a: f64, e: f64, T: f64)  -> (f64, f64) {
    let E = 2.0 * ((1.0 - e).sqrt()*(v/2.0).tan()).atan2((1.0 + e).sqrt());
    let M = E - e*E.sin();

    (T + M/n,
     a*(1.0 - e*E.cos()))
}

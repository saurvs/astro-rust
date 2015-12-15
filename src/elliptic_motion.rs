/**
Returns instantaneous **velocity** *(meters per second)* of a body in an
unperturbed elliptic orbit

# Arguments

* ```dist_to_sun```: Body's distance to Sun *(AU)*
* ```semimaj_axis```: Semimajor axis of orbit *(AU)*
**/
pub fn velocity(dist_to_sun: f64, semimaj_axis:f64) -> f64 {
    0.0421219 * (1.0/dist_to_sun - 1.0/(2.0 * semimaj_axis)).sqrt()
}

/**
Returns **velocity** *(meters per second)* of a body at **perihelion**
in an unperturbed elliptic orbit

# Arguments

* ```semimaj_axis```: Semimajor axis of orbit *(AU)*
* ```orb_eccen```: Eccentricity of orbit
**/
pub fn perih_velocity(semimaj_axis:f64, orb_eccen:f64) -> f64 {
    0.0297847 * ((1.0 + orb_eccen) / ((1.0 - orb_eccen) * semimaj_axis)).sqrt()
}

/**
Returns **velocity** *(meters per second)* of a body at **aphelion**
in an unperturbed elliptic orbit

# Arguments

* ```semimaj_axis```: Semimajor axis of orbit *(AU)*
* ```orb_eccen```: Eccentricity of orbit
**/
pub fn aph_velocity(semimaj_axis:f64, orb_eccen:f64) -> f64 {
    0.0297847 * ((1.0 - orb_eccen) / ((1.0 + orb_eccen) * semimaj_axis)).sqrt()
}

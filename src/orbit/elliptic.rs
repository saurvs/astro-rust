use std::*;

/**
Returns the **true anomaly** of a body in an elliptic orbit

# Returned value

```(true_anomaly)```

* ```(true_anomaly)```: True anomaly of the body *(AU)*

# Arguments

* ```eccentric_anomaly```: Eccentric anomaly of the body *(radians)*
* ```eccentricity```: Eccentricity of the orbit
**/
pub fn TruAnom(eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    2.0*((1.0 + eccentric_anomaly).sqrt() * (eccentricity/2.0).tan()).atan2((1.0 - eccentric_anomaly).sqrt())
}

/**
Returns the **radius vector** of a body in an elliptic orbit, using
it's eccentric anomaly

# Returned value

```(radius_vector)```

* ```(radius_vector)```: Radius vector of the body *(AU)*

# Arguments

* ```eccentric_anomaly```: Eccentric anomaly of the body *(radians)*
* ```semimaj_ax```: Semimajor axis of the orbit *(AU)*
* ```eccentricity```: Eccentricity of the orbit
**/
pub fn RadVecFrmEccAnom(eccentric_anomaly: f64, semimaj_ax: f64, eccentricity: f64) -> f64 {
    semimaj_ax*(1.0 - eccentricity*eccentric_anomaly.cos())
}

/**
Returns the **radius vector** of a body in an elliptic orbit, using
it's true anomaly

# Returned value

```(radius_vector)```

* ```(radius_vector)```: Radius vector of the body *(AU)*

# Arguments

* ```true_anomaly```: True anomaly of the body *(radians)*
* ```semimaj_ax```: Semimajor axis of the orbit *(AU)*
* ```eccentricity```: Rccentricity of the orbit
**/
pub fn RadVecFrmTruAnom(true_anomaly: f64, semimaj_ax: f64, eccentricity: f64) -> f64 {
    semimaj_ax*(1.0 - eccentricity*eccentricity) / (1.0 + eccentricity*true_anomaly.cos())
}

/**
Returns the **eccentric anomaly** of a body in an elliptic orbit

# Returned value

```(eccentric_anomaly)```

* ```(eccentric_anomaly)```: Eccentric anomaly of the body *(radians)*

# Arguments

* ```mean_anomaly```: Mean anomaly of the body *(radians)*
* ```eccentricity```: Eccentricity of the orbit
* ```accuracy```: Desired accuracy for eccentric anomaly. *Eg: 0.000001 radians*
**/
pub fn EccAnom(mean_anomaly: f64, eccentricity: f64, accuracy: f64) -> f64 {
    let mut prev_E = 0.0;
    let mut current_E = mean_anomaly;
    while (current_E - prev_E).abs() > accuracy {
        prev_E = current_E;
        current_E = mean_anomaly + eccentricity*current_E.sin();
    }

    (current_E)
}

/**
Returns the instantaneous **velocity** of a body in an
elliptic orbit

# Returns

```(velocity)```

* ```(velocity)```: Instantaneous velocity of the body
                    *(meters per second)*

# Arguments

* ```dist_to_sun```: Body's distance to Sun *(AU)*
* ```semimaj_axis```: Semimajor axis of orbit *(AU)*
**/
pub fn Vel(dist_to_sun: f64, semimaj_axis:f64) -> f64 {
    0.0421219 * (1.0/dist_to_sun - 1.0/(2.0 * semimaj_axis)).sqrt()
}

/**
Returns the **velocity** of a body at **perihelion**
in an elliptic orbit

# Returns

```(velocity)```

* ```(velocity)```: Velocity of the body at perihelion
                    *(meters per second)*

# Arguments

* ```semimaj_axis```: Semimajor axis of orbit *(AU)*
* ```orb_eccen```: Eccentricity of orbit
**/
pub fn PerihVel(semimaj_axis:f64, orb_eccen:f64) -> f64 {
    0.0297847 * ((1.0 + orb_eccen) / ((1.0 - orb_eccen) * semimaj_axis)).sqrt()
}

/**
Returns the **velocity** of a body at **aphelion** in an elliptic orbit

# Returns

```(velocity)```

* ```(velocity)```: Velocity of the body at aphelion
                    *(meters per second)*

# Arguments

* ```semimaj_axis```: Semimajor axis of orbit *(AU)*
* ```orb_eccen```: Eccentricity of orbit
**/
pub fn AphVel(semimaj_axis:f64, orb_eccen:f64) -> f64 {
    0.0297847 * ((1.0 - orb_eccen) / ((1.0 + orb_eccen) * semimaj_axis)).sqrt()
}

/**
Returns the approximate **length** of an ellipse using the Ramanujan
method

# Returns

```(approximate_length)```

* ```(approximate_length)```: An approximate value for the length of
                              the ellipse (same unit as that of ```a```
                              and ```b```), using a formula given by
                              [Ramanujan](https://en.wikipedia.org/wiki/Srinivasa_Ramanujan)
                              in 1914.

The **error** in ```(approximate_length)``` is:

* 0% for a = b
* 0.4155% for e = 1

# Arguments

* ```a```: Semimajor axis of the ellipse (same unit as that of ```b```)
* ```b```: Semiminor axis of the ellipse (same unit as that of ```a```)
* ```e```: Eccentricity of the ellipse
**/
pub fn Length_Ramanujan(a: f64, b: f64, e: f64) -> f64 {
    f64::consts::PI * (3.0*(a + b) - ((a + 3.0*b)*(3.0*a + b)).sqrt())
}

/**
Returns the approximate **length** of an ellipse

# Returns

```(approximate_length)```

* ```(approximate_length)```: An approximate value for the length of
                              the ellipse (same unit as that of ```a```
                              and ```b```)

The **error** in ```(approximate_length)``` is:

* less than 0.001% if e < 0.88
* less than 0.01% if e < 0.95
* 1% if e = 0.9997
* 3% if e = 1

# Arguments

* ```a```: Semimajor axis of the ellipse (same unit as that of ```b```)
* ```b```: Semiminor axis of the ellipse (same unit as that of ```a```)
* ```e```: Eccentricity of the ellipse
**/
pub fn Length(a: f64, b: f64, e: f64) -> f64 {
    let A = (a + b)/2.0;
    let G = (a * b).sqrt();
    let H = (2.0 * a * b)/(a + b);
    f64::consts::PI * (21.0*A - 2.0*G - 3.0*H) / 8.0
}

/**
Returns the **semimajor axis** of an elliptic orbit

# Arguments

* ```perih```: Perihelion of the orbit
* ```ecc```: Eccentricity of the orbit
**/
pub fn SemimajAx(perih: f64, ecc: f64) -> f64 {
    perih / (1.0 - ecc)
}

/**
Returns the **mean motion** of an elliptic orbit

# Returns

```(mean_motion)```

* ```(mean_motion)```: Mean motion of the orbit *(radians/day)*

# Arguments

* ```semimaj_ax```: Semimajor axis of the orbit
**/
pub fn MnMotion(semimaj_ax: f64) -> f64 {
    0.01720209895 / (semimaj_ax.powf(1.5))
}

pub fn TimeOfPassageThroughAscendNode(w: f64, n: f64, a: f64, e: f64, T: f64) -> (f64, f64) {
    time_of_passage_through_node(-w, n, a, e, T)
}

pub fn TimeOfPassageThroughDescendNode(w: f64, n: f64, a: f64, e: f64, T: f64) -> (f64, f64) {
    time_of_passage_through_node(180_f64.to_radians()*w, n, a, e, T)
}

fn time_of_passage_through_node(v: f64, n: f64, a: f64, e: f64, T: f64)  -> (f64, f64) {
    let E = 2.0 * ((1.0 - e).sqrt()*(v/2.0).tan()).atan2((1.0 + e).sqrt());
    let M = E - e*E.sin();

    (T + M/n, a*(1.0 - e*E.cos()))
}

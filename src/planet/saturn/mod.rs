//! Saturn

use angle;
use coords;
use time;

pub mod moon;
pub mod ring;

/**
Returns Saturn's apparent magnitude using G. Muller's formula

# Returns

* `app_mag`: Apparent magnitude of Saturn *| in radians*

# Arguments

* `i`     : Phase angle of Saturn *| in radians*
* `delta` : Saturn-Earth distance *| in AU*
* `r`     : Saturn-Sun distance *| in AU*
* `deltaU`: Difference between Saturnicentric longitudes of
            the Sun and the Earth, measured in the plane of
            Saturn's ring *| in radians*
* `B`     : Saturnicentric latitude of the Earth *| in radians*
**/
pub fn apprnt_mag_muller(i: f64, delta: f64, r: f64, delU: f64, B: f64) -> f64 {
    - 8.68
    + 5.0*(r*delta).log10()
    + 0.044*delU.abs()
    - 2.6*B.abs().sin()
    + 1.25*B.sin().powi(2)
}

/**
Returns Saturn's apparent magnitude using the Astronomical
Almanac's formula adopted in 1984

# Returns

* `app_mag`: Apparent magnitude of Saturn *| in radians*

# Arguments

* `i`     : Phase angle of Saturn *| in radians*
* `delta` : Saturn-Earth distance *| in AU*
* `r`     : Saturn-Sun distance *| in AU*
* `deltaU`: Difference between Saturnicentric longitudes of
            the Sun and the Earth, measured in the plane of
            Saturn's ring *| in radians*
* `B`     : Saturnicentric latitude of the Earth *| in radians*
**/
pub fn apprnt_mag_84(i: f64, delta: f64, r: f64, delU: f64, B: f64) -> f64 {
    - 8.88
    + 5.0*(r*delta).log10()
    + 0.044*delU.abs()
    - 2.6*B.abs().sin()
    + 1.25*B.sin().powi(2)
}

fn equatorial_unit_semidiameter() -> f64 { angle::deg_frm_dms(0, 0, 82.73) }
fn polar_unit_semidiameter() -> f64 { angle::deg_frm_dms(0, 0, 73.82) }

/**
Returns Saturn's polar semidiameter

# Returns

* `pol_semidiameter`: Polar semidiameter *| in radians per AU*

# Arguments

* `saturn_earth_dist`: Saturn-Earth distance *| in AU*
* `earth_lat`        : Saturnicentric latitude of Earth *| in radians*
**/
pub fn pol_semidiameter(saturn_earth_dist: f64, earth_lat: f64) -> f64 {
    let a = equatorial_unit_semidiameter();
    let b = polar_unit_semidiameter();
    let k = 1.0 - (b/a).powi(2);
    (a / saturn_earth_dist) * (1.0 - k*earth_lat.cos().powi(2)).sqrt()
}

/**
Returns Saturn's equatorial semidiameter

# Returns

* `eq_semidiameter`: Equatorial semidiameter *| in radians per AU*

# Arguments

* `saturn_earth_dist`: Saturn-Earth distance *| in AU*
**/
pub fn eq_semidiameter(saturn_earth_dist: f64) -> f64 {
    equatorial_unit_semidiameter() / saturn_earth_dist
}

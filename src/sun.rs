//! The Sun

use angle;
use time;
use std;
use planet;

/**
Returns the **equatorial semidiameter** of the Sun

# Arguments

* ```sun_earth_dist```: Sun-Earth distance *| in AU*
**/
pub fn semdiameter(sun_earth_dist: f64) -> f64 {
    angle::DegFrmDMS(0, 0, 959.63) / sun_earth_dist
}

/**
Returns **geocentric ecliptic coordinates** of the Sun,
referred to the mean equinox of the date

# Returns

```(longitude, latitude, distance)```

* ```longitude```: Ecliptic longitude of the Sun *| in radians*
* ```latitude```: Ecliptic latitude of the Sun *| in radians*
* ```distance```: Distance between the Sun and the Earth *| in AU*

# Arguments

* ```JD```: Julian (Ephemeris) day
**/
pub fn geocen_ecl_pos(JD: f64) -> (f64, f64, f64) {
    let (L, B, R) = planet::heliocen_pos(&planet::Planet::Earth, JD);

    let L_sun = angle::LimitTo360((L + std::f64::consts::PI).to_degrees());
    let B_sun = angle::LimitTo360(-B.to_degrees());

    (L_sun.to_radians(), B_sun.to_radians(), R)
}

/**
Returns geocentric ecliptic coordinates of the Sun
converted to the **FK5** system

# Returns

```(ecl_long_FK5, ecl_lat_FK5)```

* ```ecl_long_FK5```: Ecliptic longitude of the Sun *| in radians*,
                      converted to the FK5 system
* ```ecl_lat_FK5```: Ecliptic latitude of the Sun *| in radians*,
                     converted to the FK5 system

# Arguments

* ```JD```: Julian (Ephemeris) day
* ```ecl_long```: Ecliptic longitude of the Sun on ```JD```
                  *| in radians*, referred to the mean equinox
                  of the date
* ```ecl_lat```: Ecliptic latitude of the Sun ```JD```
                 *| in radians*, referred to the mean equinox
                 of the date
**/
pub fn ecl_coords_to_FK5(JD: f64, ecl_long: f64, ecl_lat: f64) -> (f64, f64) {
    let JC = time::JulCent(JD);
    let lambda1 = ecl_long - JC*(1.397 + JC*0.00031).to_radians();

    (ecl_long - angle::DegFrmDMS(0, 0, 0.09033).to_radians(),
     ecl_lat  + angle::DegFrmDMS(0, 0, 0.03916).to_radians()*(lambda1.cos() - lambda1.sin()))
}

//#[macro_export]
macro_rules! ApprntGeocenEclPos {
    ($planet: expr, $JD: expr) => {{
        3
    }};
}

/**
Returns the geocentric **rectangular** coordinates of the Sun,
referred to the **mean equinox of the date**

* The positive x-axis is directed towards the Earth's vernal equinox
(0 degrees longitude)
* The positive y-axis lies in the plane of the Earth's equator and is
directed towards 90 degrees longitude
* The positive z-axis is directed towards the Earth's northern
celestial pole
* The unit for all three axes is AU

# Returns

```(x, y z)```

* ```x```: The X coordinate *| in AU*
* ```y```: The Y coordinate *| in AU*
* ```z```: The Z coordinate *| in AU*

# Arguments

* ```sun_geo_long```: The Sun's geometric longitude *| in radians*,
                      *without* corrections for nutation and abberation
* ```sun_geo_lat```: The Sun's geometric latitude *| in radians*,
                     *without* corrections for nutation and abberation
* ```sun_rad_vec```: The Sun's geometric radius vector *| in AU*
* ```mean_oblq```: Mean obliquity of the ecliptic
**/
pub fn geocen_rect_coords(sun_geo_long: f64, sun_geo_lat: f64, sun_rad_vec: f64, mean_oblq: f64) -> (f64, f64, f64) {
    let x = sun_rad_vec * sun_geo_lat.cos() * sun_geo_long.cos();
    let y = sun_rad_vec * (sun_geo_lat.cos()*sun_geo_long.sin()*mean_oblq.cos() - sun_geo_lat.sin()*mean_oblq.sin());
    let z = sun_rad_vec * (sun_geo_lat.cos()*sun_geo_long.sin()*mean_oblq.sin() + sun_geo_lat.sin()*mean_oblq.cos());
    (x, y, z)
}

/**
Return quantites used in the **ephemeris** for **physical observations**
of the Sun

# Returns

```(P, B0, L0)```

* ```P```: Position angle of the northern extremity of the axis of
           rotation, measured eastwards from the North point of the
           solar disk *| in radians*
* ```B0```: Heliographic latitude of the center of the solar
            disk *| in radians*
* ```L0```: Heliographic longitude of the center of the solar
            disk *| in radians*

# Arguments

* ```JD```: Julian (Ephemeris) day
* ```app_long```: Apparent longitude of the Sun *| in radians*,
                  including the effect of abberation and *not* that
                  of nutation
* ```app_long_with_nut```: Apparent longitude of the Sun *| in radians*,
                  including the effect of abberation *and* that
                  of nutation
* ```oblq_eclip```: True obliquity of the ecliptic *| in radians*
**/
pub fn ephemeris(JD: f64, app_long: f64, app_long_with_nut: f64, oblq_eclip: f64) -> (f64, f64, f64) {
    let theta = angle::LimitTo360((JD - 2398220.0) * (360.0/25.38)).to_radians();
    let I = 7.25_f64.to_radians();
    let K = (73.6667 + 1.3958333*((JD - 2396758.0) / 36525.0)).to_radians();

    let z = app_long - K;
    let sin_z = z.sin();
    let cos_z = z.cos();

    let mut x = (-app_long_with_nut.cos() * oblq_eclip.tan()).atan();
    let mut y = (-cos_z * I.tan()).atan();
    x = (magnitude_limited_to_less_than_90(x.to_degrees())).to_radians();
    y = (magnitude_limited_to_less_than_90(y.to_degrees())).to_radians();

    let B_0 = (sin_z * I.sin()).asin();
    let nu = (-sin_z * I.cos()).atan2(-cos_z);
    let L_0 = angle::LimitTo360((nu - theta).to_degrees()).to_radians();

    (x + y, B_0, L_0)
}

pub fn carring_synd_rot(C: i64) -> f64 {
    let M = (281.96 + 26.882476*(C as f64)).to_radians();

      2398140.227 + 27.2752316*(C as f64)
    + 0.1454 * M.sin()
    - 0.0085 * (2.0*M).sin()
    - 0.0141 * (2.0*M).cos()
}

fn magnitude_limited_to_less_than_90(a: f64) -> f64 {
    if a > 270.0 {a - 360.0}
    else         {a}
}

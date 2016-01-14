use angle;
use planet;
use std;

/**
Returns the **equatorial semidiameter** of the Sun

# Arguments

* ```distance_to_earth```: The Sun's distance from the Earth *(AU)*
**/
pub fn Semdia(distance_to_earth: f64) -> f64 {
    angle::DegFrmDMS(0, 0, 959.63) / distance_to_earth
}

/**
Returns the **ecliptic geocentric coordinates** of the Sun

# Returns

```(longitude, latitude, distance)```

* ```longitude```: Ecliptic longitude of the Sun *(radians)*
* ```latitude```: Ecliptic latitude of the Sun *(radians)*
* ```distance```: Distance between the Sun and the Earth *(kilometers)*

# Arguments

* ```JD```: Julian (Emphemeris) day
**/
pub fn EclGeocenCoords(JD: f64) -> (f64, f64, f64) {
    let (L, B, R) = planet::HeliocenCoords(planet::Planet::Earth, JD);

    let L_sun = angle::LimitTo360((L + std::f64::consts::PI).to_degrees());
    let B_sun = angle::LimitTo360(-1.0 * B.to_degrees());

    (L_sun.to_radians(), B_sun.to_radians(), R)
}

/**
Returns the **rectangular geocentric coordinates** of the Sun

* The positive x-axis is directed towards the Earth's vernal equinox
(0 degrees longitude)
* The positive y-axis lies in the plane of the Earth's equator and is
directed towards 90 degrees longitude
* The positive z-axis is directed towards the Earth's northern
celestial pole
* The unit for all three axes is AU

# Returns

```(x, y z)```

* ```x```: The X coordinate *(AU)*
* ```y```: The Y coordinate *(AU)*
* ```z```: The Z coordinate *(AU)*

# Arguments

* ```sun_geo_long```: The Sun's geometric longitude *(radians)*,
                      *without* corrections for nutation and abberation
* ```sun_geo_lat```: The Sun's geometric latitude *(radians)*,
                     *without* corrections for nutation and abberation
* ```sun_rad_vec```: The Sun's geometric radius vector *(AU)*
* ```mean_obl```: The *mean* obliquity of the Earth's ecliptic;
                  not *true* obliquity
**/
pub fn RectGeocenCoords(sun_geo_long: f64, sun_geo_lat: f64, sun_rad_vec: f64, mean_obl: f64) -> (f64, f64, f64) {
    let x = sun_rad_vec * sun_geo_lat.cos() * sun_geo_long.cos();
    let y = sun_rad_vec * (sun_geo_lat.cos()*sun_geo_long.sin()*mean_obl.cos() - sun_geo_lat.sin()*mean_obl.sin());
    let z = sun_rad_vec * (sun_geo_lat.cos()*sun_geo_long.sin()*mean_obl.sin() + sun_geo_lat.sin()*mean_obl.cos());
    (x, y, z)
}

/**
Return quantites used in the **emphemeris** for physical observations
of the Sun

# Returns

```(P, B0, L0)```

* ```P```: Position angle of the northern extremity of the axis of
           rotation, measured eastwards from the North point of the
           solar disk *(radians)*
* ```B0```: Heliographic latitude of the center of the solar
            disk *(radians)*
* ```L0```: Heliographic longitude of the center of the solar
            disk *(radians)*

# Arguments

* ```JD```: Julian (Emphemeris) day
* ```app_long```: Apparent longitude of the Sun *(radians)*,
                  including the effect of abberation and *not* that
                  of nutation
* ```app_long_with_nut```: Apparent longitude of the Sun *(radians)*,
                  including the effect of abberation *and* that
                  of nutation
* ```oblq_eclip```: *True* obliquity of the Earth's ecliptic *(radians)*,
                    i.e, *with* correction for nutation
**/
pub fn Ephm(JD: f64, app_long: f64, app_long_with_nut: f64, oblq_eclip: f64) -> (f64, f64, f64) {
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

pub fn CarringSyndRot(C: i64) -> f64 {
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

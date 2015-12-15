use angle;

/**
Computes the Sun's **equatorial semidiameter**

# Arguments

* ```distance_to_earth```: The Sun's distance from the Earth *(AU)*
**/
pub fn semidiameter(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 959.63) / distance_to_earth
}

/**
Computes **rectangular geocentric equatorial coordinates** of the Sun

# Return variables

* The positive x-axis is directed towards the Earth's vernal equinox
(0 degrees longitude)
* The positive y-axis lies in the plane of the Earth's equator and is directed
towards 90 degrees longitude
* The positive z-axis is directed towards the Earth's northern celestial pole
* The unit for all three axes is AU

```rect_geocen_coords() -> (x, y z)```

# Arguments

* ```sun_geo_long```: Sun's geometric longitude *(radians)*, *without* corrections for nutation and abberation
* ```sun_geo_lat```: Sun's geometric latitude *(radians)*, *without* corrections for nutation and abberation
* ```sun_rad_vec```: Sun's geometric radius vector *(AU)*
* ```mean_obl```: *Mean* obliquity of Earth's ecliptic; not *true* obliquity

**/

pub fn rect_geocen_coords(sun_geo_long: f64, sun_geo_lat: f64, sun_rad_vec: f64, mean_obl: f64) -> (f64, f64, f64) {
    let x = sun_rad_vec * sun_geo_lat.cos() * sun_geo_long.cos();
    let y = sun_rad_vec * (sun_geo_lat.cos()*sun_geo_long.sin()*mean_obl.cos() - sun_geo_lat.sin()*mean_obl.sin());
    let z = sun_rad_vec * (sun_geo_lat.cos()*sun_geo_long.sin()*mean_obl.sin() + sun_geo_lat.sin()*mean_obl.cos());
    (x, y, z)
}

pub fn disk_ephemeris(jde: f64, app_long: f64, app_long_with_nut: f64, obl_eclip: f64) -> (f64, f64, f64) {
    let theta = (jde - 2398220.0) * (360.0 / 25.38);
    let I: f64 = 7.25;
    let K = 73.6667 + 1.3958333*((jde - 2396758.0) / 36525.0);

    let z = app_long - K;
    let sin_z = z.sin();
    let cos_z = z.cos();

    let mut x = (-1.0 * app_long_with_nut.cos() * obl_eclip.tan()).atan();
    let mut y = (-1.0 * cos_z * I.tan()).atan();
    x = (magnitude_limited_to_less_than_90(x.to_degrees())).to_radians();
    y = (magnitude_limited_to_less_than_90(y.to_degrees())).to_radians();

    let P = x + y;
    let B_0 = (sin_z * I.sin()).asin();
    let nu = (-1.0 * sin_z * I.cos()).atan2(-1.0 * cos_z);
    let L_0 = angle::limited_to_360((nu - theta).to_degrees()).to_radians();

    (P, B_0, L_0)
}

fn magnitude_limited_to_less_than_90(a: f64) -> f64 {
    if a > 270.0 {a - 360.0}
    else {a}
}

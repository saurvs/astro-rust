use angle;

/**
Returns the Sun's equatorial angular semidiameter

# Arguments

* ```distance_to_earth```: The Sun's distance from the Earth (in AU)
**/
pub fn semidiameter(distance_to_earth: f64) -> f64 {
    angle::pure_degrees(0.0, 0.0, 959.63) / distance_to_earth
}

/**
Returns rectangular geocentric equatorial coordinates of the Sun

# Return variables

* The positive x-axis is directed towards Earth's vernal equinox
(0 degrees longitude)
* The positive y-axis lies in the plane of Earth's equator and is directed
towards 90 degrees longitude
* The positive z-axis is directed towards Earth's northern celestial pole
* The unit for all three axes is AU

```rect_coords() -> (x, y z)```

# Arguments

* ```sun_geo_long```: Sun's geometric longitude (in radians), without corrections for nutation and abberation
* ```sun_geo_lat```: Sun's geometric latitude (in radians), without corrections for nutation and abberation
* ```sun_rad_vec```: Sun's geometric radius vector (in AU)
* ```mean_obl```: *Mean* obliquity of Earth's ecliptic; not *true* obliquity

**/

pub fn rect_coords(sun_geo_long: f64, sun_geo_lat: f64, sun_rad_vec: f64, mean_obl: f64) -> (f64, f64, f64) {
    let x = sun_rad_vec * sun_geo_lat.cos() * sun_geo_long.cos();
    let y = sun_rad_vec * (sun_geo_lat.cos()*sun_geo_long.sin()*mean_obl.cos() - sun_geo_lat.sin()*mean_obl.sin());
    let z = sun_rad_vec * (sun_geo_lat.cos()*sun_geo_long.sin()*mean_obl.sin() + sun_geo_lat.sin()*mean_obl.cos());
    (x, y, z)
}

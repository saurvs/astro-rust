/*

    NOTE: All angles passed as arguments, and those returned,
          are assumed to be radians, even if the comments
          describe them with degrees.

*/

/*

    rect_coord(sun_geocentric_longitude, sun_geocentric_latitude,
               sun_radius_vector, mean_obliquity_of_ecliptic)
               -> (x, y z)
    -----------------------------------------------------------------
    Returns rectangular coordinates of the sun.

    The origin of these coordinates is the center of the Earth.

    The positive x-axis is directed towards the vernal equinox, that is,
    longitude 0 degrees.
    The positive y-axis lies in the plane of the equator and is directed
    towards longitude 90 degrees.
    The positive z-axis is directed towards the northern celestial pole.

*/

pub fn rect_coord(sun_geo_long: f64, sun_geo_lat: f64, sun_rad_vec: f64, mean_obl: f64) -> (f64, f64, f64) {
    let x = sun_rad_vec * sun_geo_lat.cos() * sun_geo_long.cos();
    let y = sun_rad_vec * (sun_geo_lat.cos()*sun_geo_long.sin()*mean_obl.cos() - sun_geo_lat.sin()*mean_obl.sin());
    let z = sun_rad_vec * (sun_geo_lat.cos()*sun_geo_long.sin()*mean_obl.sin() + sun_geo_lat.sin()*mean_obl.cos());
    (x, y, z)
}

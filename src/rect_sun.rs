/*

    Returns the geocentric rectangle coordinates of the sun.
    -----------------------------------------------------------------
        geo_long_sun: The geocentric longitude of the sun
         rad_vec_sun: The geocentric latitude of the sun
         rad_vec_sun: The geocentric radius vector of the sun
            mean_obl: The mean obliquity of the ecliptic

*/

pub fn rect_coord(geo_long_sun: f64, geo_lat_sun: f64, rad_vec_sun: f64, mean_obl: f64) -> (f64, f64, f64) {
    let x = rad_vec_sun * geo_lat_sun.cos() * geo_long_sun.cos();
    let y = rad_vec_sun * (geo_lat_sun.cos()*geo_long_sun.sin()*mean_obl.cos() - geo_lat_sun.sin()*mean_obl.sin());
    let z = rad_vec_sun * (geo_lat_sun.cos()*geo_long_sun.sin()*mean_obl.sin() + geo_lat_sun.sin()*mean_obl.cos());
    (x, y, z)
}

use coordinates;

/*

    NOTE: All angles passed as arguments, and those returned,
          are assumed to be radians, even if the comments
          describe them with degrees.

*/

/*
    This code uses revised values for flattening factor and
    equatorial radius of the earth.
    See: http://www.unoosa.org/pdf/icg/2012/template/WGS_84.pdf
    or https://confluence.qps.nl/pages/viewpage.action?pageId=29855173
*/

/*

    flattening() -> (flattening_factor_of_earth)
    -----------------------------------------------------------------
    Returns the flattening factor of the earth

*/

pub fn flattening() -> f64 {
    1.0 / 298.257223563
}

/*

    eq_radius() -> (equatorial_radius_of_earth)
    -----------------------------------------------------------------
    Returns the equatorial radius of the earth (in meters)

*/

pub fn eq_radius() -> f64 {
    polar_radius() / (1.0 - flattening())
}

/*

    polar_radius() -> (polar_radius_of_earth)
    -----------------------------------------------------------------
    Returns the polar radius of the earth (in meters)

*/

pub fn polar_radius() -> f64 {
    6378137.0
}

/*

    ecc() -> (eccentricity_of_earth_meridian)
    -----------------------------------------------------------------
    Returns the eccentricity of the earth's meridian

*/

pub fn ecc() -> f64 {
    let f = flattening();
    ((2.0 - f) * f).sqrt()
}

/*

    angular_dist(point_1, point_2) -> (angular_distance)
    -----------------------------------------------------------------
    Returns the angular distance between two points on Earth's
    surface

*/

pub fn angular_dist(p1: coordinates::surf_point, p2: coordinates::surf_point) -> f64 {
    (p1.lat.sin() * p2.lat.sin() +
     p1.lat.cos() * p2.lat.cos() * (p1.long - p2.long).cos()
    ).acos()
}

/*

    approx_dist(point_1, point_2) -> (distance)
    -----------------------------------------------------------------
    Returns a low accuracy distance between two points on earth's
    surface (in meters). Assumes the Earth is spherical.

*/

pub fn approx_dist(p1: coordinates::surf_point, p2: coordinates::surf_point) -> f64 {
    6371.0 * angular_dist(p1, p2)
}

/*

    dist(point_1, point_2) -> (distance)
    -----------------------------------------------------------------
    Returns a high accuracy distance between two points on earth's
    surface (in meters)

*/

pub fn dist(p1: coordinates::surf_point, p2: coordinates::surf_point) -> f64 {
    let f = (p1.lat + p2.lat) / 2.0;
    let g = (p1.lat - p2.lat) / 2.0;
    let lam = (p1.long - p2.long) / 2.0;
    let s = (g.sin() * lam.cos()).powi(2) +
            (f.cos() * lam.sin()).powi(2);
    let c = (g.cos() * lam.cos()).powi(2) +
            (f.sin() * lam.sin()).powi(2);
    let om = ((s / c).sqrt()).atan();
    let r = (s * c).sqrt() / om;
    let d = 2.0 * om * eq_radius();
    let h1 = (3.0 * r - 1.0) / (2.0 * c);
    let h2 = (3.0 * r + 1.0) / (2.0 * s);

    d * (1.0 +
         flattening() * h1 * (f.sin() * g.cos()).powi(2) -
         flattening() * h2 * (f.cos() * g.sin()).powi(2))

}

/*

    rho_sin_and_cos_phi(height, geograph_lat)
                                        -> (rho_sin_phi, rho_cos_phi)
    -----------------------------------------------------------------
    Returns two quantities that are used elsewhere in the library.

    'rho' here denotes the geocentric radius vector, and 'phi' here
    denotes the geocentric latitude, both of an observer on the
    Earth's surface.

    height: The observer's height above sea level (meters)
    geograph_lat: The observer's geographical latitude

*/

pub fn rho_sin_and_cos_phi(height: f64, geograph_lat: f64) -> (f64, f64) {
    let u = (geograph_lat.tan() * polar_radius() / eq_radius()).atan();
    let x = height / eq_radius();
    let rho_sin_phi = (u.sin() * polar_radius() / eq_radius()) +
                      (geograph_lat.sin() * x);
    let rho_cos_phi = u.cos() + (geograph_lat.cos() * x);

    (rho_sin_phi, rho_cos_phi)
}

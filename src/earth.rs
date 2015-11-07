// flattening factor and equatorial radius of the earth used here
// aren't those given in the book; they have been revised since 2004
// by the World Geodetic System.
// See: http://www.unoosa.org/pdf/icg/2012/template/WGS_84.pdf
// or https://confluence.qps.nl/pages/viewpage.action?pageId=29855173

// flattening factor of the earth
pub fn flattening() -> f64 {
    1.0 / 298.257223563
}

// equatorial radius of the earth (in meters)
pub fn eq_radius() -> f64 {
    6378137.0
}

// polar radius of the earth (in meters)
pub fn pol_radius() -> f64 {
    eq_radius() / (1.0 - flattening())
}

// eccentricity of the earth's meridian
pub fn ecc() -> f64 {
    let f = flattening();
    ((2.0 - f) * f).sqrt()
}

// struct for representing a point on the surface of the earth using
// latitude and longitude (in radians)
pub struct surf_point {
    pub lat: f64,
    pub long: f64,
}

// angular distance between two points on earth's surface (in radians)
pub fn angular_dist(p1: surf_point, p2: surf_point) -> f64 {
    (p1.lat.sin() * p2.lat.sin() +
     p1.lat.cos() * p2.lat.cos() * (p1.long - p2.long).cos()
    ).acos()
}

// low accuracy distance between two points on earth's surface (in meters)
// assumes earth is spherical
pub fn approx_dist(p1: surf_point, p2: surf_point) -> f64 {
    6371.0 * angular_dist(p1, p2)
}

// high accuracy distance between two points on earth's surface (in meters)
pub fn dist(p1: surf_point, p2: surf_point) -> f64 {
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

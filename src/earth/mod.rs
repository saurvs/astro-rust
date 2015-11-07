pub mod earth;

// flattening factor and equatorial radius of the earth
// aren't those given in the book; they have been revised
// since 2004 by the World Geodetic System.
// See: http://www.unoosa.org/pdf/icg/2012/template/WGS_84.pdf

// flattening factor of the earth
fn flattening() -> f64 {
    298.257223563
}

// equatorial radius of the earth (in kilometers)
fn eq_radius() -> f64 {
    6378.137
}

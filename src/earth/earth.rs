// flattening factor and equatorial radius of the earth
// aren't those given in the book; they have been revised
// since 2004 by the World Geodetic System.
// See: http://www.unoosa.org/pdf/icg/2012/template/WGS_84.pdf
// or https://confluence.qps.nl/pages/viewpage.action?pageId=29855173

// flattening factor of the earth
fn flattening() -> f64 {
    298.257223563
}

// equatorial radius of the earth (in kilometers)
fn eq_radius() -> f64 {
    6378.137
}

// polar radius of the earth (in kilometers)
fn pol_radius() -> f64 {
    eq_radius() / (1.0 - flattening())
}

// eccentricity the earth's meridian
fn ecc() -> f64 {
    let f = flattening();
    (2.0 * f - f * f).sqrt()
}

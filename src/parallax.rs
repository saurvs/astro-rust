use earth;

pub fn TopocentricEquatorialCoords(asc: f64, dec: f64, height: f64, geograph_lat: f64) -> (f64, f64) {
    let (rho_sin, rho_cos) = earth::RhoSinAndCosPhi(height, geograph_lat);
}

pub fn TopocentricEclipticalCoords(asc: f64, dec: f64, height: f64, geograph_lat: f64) -> (f64, f64) {
    let (rho_sin, rho_cos) = earth::RhoSinAndCosPhi(height, geograph_lat);
}

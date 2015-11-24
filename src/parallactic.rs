pub fn parallactic_angle(geograph_lat: f64, local_hour_angle: f64, dec_celes_body: f64) -> f64 {
    (local_hour_angle.sin() /
     geograph_lat.tan()*dec_celes_body.cos() - dec_celes_body.sin()*local_hour_angle.cos()).atan()

}

pub fn parallactic_angle_on_horizon(geograph_lat: f64, dec_celes_body: f64) -> f64 {
    (geograph_lat.sin() / dec_celes_body.cos()).acos()
}

pub fn longitudes_of_ecl_on_horzion(obl_ecl: f64, observ_lat: f64, loc_sid: f64) -> (f64, f64) {
    let lambda = (-1.0 * loc_sid.cos()).atan2(obl_ecl.sin()*observ_lat.tan() + obl_ecl.cos()*loc_sid.sin());
    (lambda, 180_f64.to_radians() + lambda)
}

pub fn angle_between_ecl_and_horizon(obl_ecl: f64, observ_lat: f64, loc_sid: f64) -> f64 {
    (obl_ecl.cos()*observ_lat.sin() - obl_ecl.sin()*observ_lat.cos()*loc_sid.sin()).acos()
}

pub fn angle_between_north_celes_pole_and_north_ecl_pole(ecl_long: f64, ecl_lat: f64, obl_ecl: f64) -> f64 {
    let x = obl_ecl.tan();
    (ecl_long.cos()*x) / (ecl_lat.sin()*ecl_long.sin()*x - ecl_lat.cos())
}

pub fn angle_between_diurnal_path_and_horizon(dec_celes_body: f64, observ_lat: f64) -> f64 {
    let x = observ_lat.tan();
    let B = dec_celes_body.tan() * x;
    let C = (1.0 - B*B).sqrt();
    (C * dec_celes_body.cos()).atan2(x)
}

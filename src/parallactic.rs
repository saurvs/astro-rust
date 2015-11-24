pub fn parallactic_angle(geograph_lat: f64, local_hour_angle: f64, dec_celes_body: f64) -> f64 {
    (local_hour_angle.sin() /
     geograph_lat.tan()*dec_celes_body.cos() - dec_celes_body.sin()*local_hour_angle.cos()).atan()

}

pub fn horizon_parallactic_angle(geograph_lat: f64, dec_celes_body: f64) -> f64 {
    (geograph_lat.sin() / dec_celes_body.cos()).acos()
}

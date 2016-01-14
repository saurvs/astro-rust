use angle;
use planet;

pub fn EqHzParllx(dist_to_earth: f64) -> f64 {
    (angle::DegFrmDMS(0, 0, 8.794).to_radians().sin() / dist_to_earth).asin()
}

pub fn TopocenEqCoords(asc: f64, dec: f64, observer_ht: f64, geograph_lat: f64,
                       eq_hz_parllx: f64, geocen_hr_angl: f64,
                       ) -> (f64, f64) {

    let (rho_sin, rho_cos) = planet::earth::RhoSinAndCosPhi(observer_ht, geograph_lat);

    let eq_hz_parllx_sin = eq_hz_parllx.sin();
    let del_asc = (-rho_cos * eq_hz_parllx_sin * geocen_hr_angl.sin())
                  .atan2(dec.cos() - rho_cos*eq_hz_parllx_sin*geocen_hr_angl.cos());println!("{:?}", del_asc.to_degrees());
    let dec_1 = ((dec.sin() - rho_sin*eq_hz_parllx_sin) * del_asc.cos())
                .atan2(dec.cos() - rho_cos*eq_hz_parllx_sin*geocen_hr_angl.cos());

    (asc + del_asc, dec_1)
}

pub fn TopocenEclCoords(ecl_long: f64, ecl_lat: f64, eq_hz_parllx: f64,
                        observer_ht: f64, geograph_lat: f64,
                        loc_sidr: f64, eclip_oblq: f64,
                        geocen_semdia: f64) -> (f64, f64, f64) {

    let (rho_sin, rho_cos) = planet::earth::RhoSinAndCosPhi(observer_ht, geograph_lat);

    let eq_hz_parllx_sin = eq_hz_parllx.sin();
    let N = ecl_long.cos()*ecl_lat.cos() - rho_cos*eq_hz_parllx_sin*loc_sidr.cos();
    let ecl_long_1 = (  ecl_long.sin()*ecl_lat.cos()
                      - eq_hz_parllx_sin*(rho_sin*eclip_oblq.sin() + rho_cos*eclip_oblq.cos()*loc_sidr.sin()))
                     .atan2(N);
    let ecl_lat_1 = (ecl_long_1.cos()*(ecl_lat.sin() - eq_hz_parllx_sin*(  rho_sin*eclip_oblq.cos()
                                                                         - rho_cos*eclip_oblq.sin()*loc_sidr.sin())))
                    .atan2(N);
    let geocen_semdia_1 = (ecl_long_1.cos()*ecl_lat_1.cos()*geocen_semdia.sin() / N).asin();

    (ecl_long_1, ecl_lat_1, geocen_semdia_1)
}

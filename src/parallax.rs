//! Find topocentric coordinates

use angle;
use coords;
use planet;

/**
Returns the **equatorial horizontal parallax** of a celestial body

# Returns

* `eq_hz_parllx`: Equatorial horizontal parallax of the celestial
                      body *| in radians*

# Arguments

* `dist_to_earth`: The celestial body's distance to the Earth *| in AU*
**/
pub fn eq_hz_parallax(dist_to_earth: f64) -> f64 {
    (angle::deg_frm_dms(0, 0, 8.794).to_radians().sin() / dist_to_earth).asin()
}

/**
Returns the **topocentric equatorial coordinates** of a celestial body

# Returns

* `topocen_eq_point`: Topocentric equatorial point of the
                     celestial body *| in radians*

# Arguments

* `eq_point`: Equatorial point of the celestial body *| in radians*
* `dist_to_earth`: The celestial body's distance to the Earth *| in AU*
* `geograph_point`: Geographic point of the observer *| in radians*
* `observer_ht`: Height of the observer above sea level *| in meters*
* `greenwhich_sidreal`: Sidereal time at Greenwhich *| in radians*
**/
pub fn topocen_eq_coords(
    eq_point: &coords::EqPoint,
    dist_to_earth: f64,
    geograph_point: &coords::GeographPoint,
    observer_ht: f64,
    greenwhich_sidreal: f64) -> coords::EqPoint {

    let (rho_sin, rho_cos) = planet::earth::rho_sin_cos_phi(geograph_point.lat, observer_ht);
    let geocen_hr_angl = coords::hr_angl_frm_observer_long(greenwhich_sidreal, geograph_point.long, eq_point.asc);

    let eq_hz_parllx_sin = eq_hz_parallax(dist_to_earth).sin();
    let del_asc = (-rho_cos*eq_hz_parllx_sin*geocen_hr_angl.sin())
                  .atan2(eq_point.dec.cos() - rho_cos*eq_hz_parllx_sin*geocen_hr_angl.cos());
    let dec_1 = ((eq_point.dec.sin() - rho_sin*eq_hz_parllx_sin) * del_asc.cos())
                .atan2(eq_point.dec.cos() - rho_cos*eq_hz_parllx_sin*geocen_hr_angl.cos());

    coords::EqPoint{
        asc: eq_point.asc + del_asc,
        dec: dec_1
    }
}

/**
Returns the **topocentric ecliptic coordinates** of a celestial body

# Returns

`(topocen_ecl_point, topocen_geocen_semidia)`

* `topocen_ecl_point`: Topocentric ecliptic point of the
                     celestial body *| in radians*
* `topocen_geocen_semidia`: Topocentric geocentric semidiameter of the celestial body *| in radians*

# Arguments

* `ecl_point`: Ecliptic point of the celestial body *| in radians*
* `eq_hz_parllx`: Equatorial horizontal parallax of the celestial body *| in radians*
* `geograph_point`: Geographic point of the observer *| in radians*
* `observer_ht`: Height of the observer above sea level *| in meters*
* `loc_sidr`: Local sidereal time *| in radians*
* `eclip_oblq`: Obliquity of the ecliptic *| in radians*
* `geocen_semdia`: Geocentric semidiameter of the celestial body *| in radians*
**/
pub fn topopcen_ecl_coords(
    ecl_point: &coords::EclPoint,
    eq_hz_parllx: f64,
    geograph_point: &coords::GeographPoint,
    observer_ht: f64,
    loc_sidr: f64,
    eclip_oblq: f64,
    geocen_semdia: f64) -> (coords::EclPoint, f64) {

    let (rho_sin, rho_cos) = planet::earth::rho_sin_cos_phi(geograph_point.lat, observer_ht);

    let eq_hz_parllx_sin = eq_hz_parllx.sin();
    let N = ecl_point.long.cos()*ecl_point.lat.cos() - rho_cos*eq_hz_parllx_sin*loc_sidr.cos();
    let ecl_long_1 = (  ecl_point.long.sin()*ecl_point.lat.cos()
                      - eq_hz_parllx_sin*(rho_sin*eclip_oblq.sin() + rho_cos*eclip_oblq.cos()*loc_sidr.sin()))
                     .atan2(N);
    let ecl_lat_1 = (ecl_long_1.cos()*(ecl_point.lat.sin() - eq_hz_parllx_sin*(  rho_sin*eclip_oblq.cos()
                                                                         - rho_cos*eclip_oblq.sin()*loc_sidr.sin())))
                    .atan2(N);
    let geocen_semdia_1 = (ecl_long_1.cos()*ecl_lat_1.cos()*geocen_semdia.sin() / N).asin();

    (coords::EclPoint{
         long: angle::limit_to_360(ecl_long_1.to_degrees()).to_radians(),
         lat: ecl_lat_1
     },
      geocen_semdia_1
    )
}

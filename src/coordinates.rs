// struct for representing a point on the surface of the earth using
// latitude and longitude (in radians)
pub struct surf_point {
    pub lat: f64,
    pub long: f64,
}

/*

    A struct for representing a point on the celestial sphere using
    in the equatorial coordinate system.
    -----------------------------------------------------------------
        asc: Right ascension
        dec: Declination

*/

pub struct celes_point {
    pub asc: f64,
    pub dec: f64,
}

// obliquity of the ecliptic for the 2000 epoch
pub fn oblq_ecl_2000() -> f64 {
    23.4392911_f64.to_radians()
}

// obliquity of the ecliptic for the 1950 epoch
pub fn oblq_ecl_1950() -> f64 {
    23.4457889_f64.to_radians()
}

// return the hour angle from Greenwhich sidereal time, observer longitude and right ascension
pub fn hour_angle_from_greenwhich_sidereal(green_sid: f64, obv_long: f64, right_asc: f64) -> f64 {
    green_sid - obv_long - right_asc
}

// return the hour angle from local sidereal time and right ascension
pub fn hour_angle_from_local_sidereal(local_sid: f64, right_asc: f64) -> f64 {
    local_sid - right_asc
}

//-------------------------------------------------------------------
// equatorial coordinates to ecliptical coordinates

pub fn ecl_long(right_asc: f64, dec: f64, oblq_ecl: f64,) -> f64 {
    ((right_asc.sin() * oblq_ecl.cos() + dec.tan() * oblq_ecl.sin())).atan2(right_asc.cos())
}

pub fn ecl_lat(right_asc: f64, dec: f64, oblq_ecl: f64) -> f64 {
    (dec.sin() * oblq_ecl.cos() - dec.cos() * oblq_ecl.sin() * right_asc.sin()).asin()
}

//-------------------------------------------------------------------
// ecliptic coordinates to equatorial coordinates

pub fn right_asc(ecl_long: f64, ecl_lat: f64, oblq_ecl: f64) -> f64 {
    ((ecl_long.sin() * oblq_ecl.cos() - ecl_lat.tan() * oblq_ecl.sin())).atan2(ecl_long.cos())
}

pub fn dec(ecl_long: f64, ecl_lat: f64, oblq_ecl: f64) -> f64 {
    (ecl_lat.sin() * oblq_ecl.cos() - ecl_lat.cos() * oblq_ecl.sin() * ecl_long.sin()).asin()
}

//-------------------------------------------------------------------
// local horizontal coordinates from equatorial coordinates

pub fn azimuth(hour_angle: f64, eq_dec: f64, obv_lat: f64) -> f64 {
    hour_angle.sin().atan2(hour_angle.cos() * obv_lat.sin() - eq_dec.tan() * obv_lat.cos())
}

pub fn altitude(hour_angle: f64, eq_dec: f64, obv_lat: f64) -> f64 {
    (obv_lat.sin() * eq_dec.sin() + obv_lat.cos() * eq_dec.cos() * hour_angle.cos()).asin()
}

//-------------------------------------------------------------------
// horizontal coordinates to equatorial coordinates

pub fn hour_angle(azimuth: f64, altitude: f64, obv_lat: f64) -> f64 {
    azimuth.sin().atan2(azimuth.cos() * obv_lat.sin() + altitude.tan() * obv_lat.cos())
}

pub fn eq_dec(azimuth: f64, altitude: f64, obv_lat: f64) -> f64 {
    (obv_lat.sin() * altitude.sin() - obv_lat.cos() * azimuth.cos() * azimuth.cos()).asin()
}

//-------------------------------------------------------------------
// equatorial coordinates to galactic coordinates

pub fn gal_long(right_asc: f64, dec: f64) -> f64 {

    let x = ((3.35539549 - right_asc).sin())
            .atan2(((3.35539549 - right_asc).cos() * (0.4782202_f64).sin() -
                     dec.tan() * (0.4782202_f64).cos()
                   ));

    5.28835 - x
}

pub fn gal_lat(right_asc: f64, dec: f64) -> f64 {
    dec.sin() * (0.4782202_f64).sin() +
    dec.cos() * (0.4782202_f64).cos() * (3.35539549 - right_asc).asin()
}

//-------------------------------------------------------------------
// galactic coordinates to equatorial coordinates

pub fn right_asc_from_gal(gal_long: f64, gal_lat: f64) -> f64 {

    let y = ((gal_long - 2.14675).sin())
            .atan2(((gal_long - 2.14675).cos() * (0.4782202_f64).sin() -
                     gal_lat.tan() * (0.4782202_f64).cos()
                   ));

    0.21380283 + y
}

pub fn dec_from_gal(gal_long: f64, gal_lat: f64) -> f64 {
    gal_lat.sin() * (0.4782202_f64).sin() +
    gal_lat.cos() * (0.4782202_f64).cos() * (gal_long - 2.14675).asin()
}

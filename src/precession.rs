use time;
use angle;

/**
Returns stellar **equatorial coordinates** for a **different equinox**

# Returns

```(new_right_ascension, new_declination)```

* ```new_right_ascension```: Right ascension for the new epoch *(radians)*
* ```new_declination```: Declination for the new epoch *(radians)*

# Arguments

* ```jd_1```: Julian day corresponding to the old epoch
* ```jd_2```: Julian day corresponding to the new epoch
* ```asc_old```: Right ascension for the old epoch *(radians)*
* ```dec_old```: Declination for the old epoch *(radians)*
**/
pub fn ChangeEpochForEquatorialCoords(jd_1: f64, jd_2: f64, asc_old: f64, dec_old: f64) -> (f64, f64) {
    let T = time::JulianCentury(jd_1);
    let t = (jd_2 - jd_1) / 36525.0;

    let x = t * (angle::PureDegrees(0.0, 0.0, 2306.2181) +
                  T * (angle::PureDegrees(0.0, 0.0, 1.39656) -
                       T*angle::PureDegrees(0.0, 0.0, 0.000139)));
    let a = (x + t*t*((angle::PureDegrees(0.0, 0.0, 0.30188) -
                       T*angle::PureDegrees(0.0, 0.0, 0.000344)) +
                      t*angle::PureDegrees(0.0, 0.0, 0.017998))).to_radians();

    let b = (x + t*t*((angle::PureDegrees(0.0, 0.0, 1.09468) -
                       T*angle::PureDegrees(0.0, 0.0, 0.000066)) +
                      t*angle::PureDegrees(0.0, 0.0, 0.018203))).to_radians();

    let y = T * angle::PureDegrees(0.0, 0.0, 0.000217);
    let c = (t * (angle::PureDegrees(0.0, 0.0, 2004.3109) +
                   T * (angle::PureDegrees(0.0, 0.0, 0.8533) - y) -
                  t * ((angle::PureDegrees(0.0, 0.0, 0.42665) + y) +
                       t*angle::PureDegrees(0.0, 0.0, 0.041833)))).to_radians();

    let sin_dec_old = dec_old.sin();
    let cos_dec_old = dec_old.cos();
    let cos_asc_old_plus_a = (asc_old + a).cos();
    let cos_c = c.cos();
    let sin_c = c.sin();

    let A = cos_dec_old * (asc_old + a).sin();
    let B = cos_c*cos_dec_old*cos_asc_old_plus_a - sin_c*sin_dec_old;
    let C = sin_c*cos_dec_old*cos_asc_old_plus_a + cos_c*sin_dec_old;

    (A.atan2(B) + b, C.asin())
}

/**
Returns stellar **ecliptical coordinates** for a **different equinox**

# Returns

```(new_longitude, new_latitude)```

* ```new_longitude```: Longitude for the new epoch *(radians)*
* ```new_latitude```: Latitude for the new epoch *(radians)*

# Arguments

* ```jd_1```: Julian day corresponding to the old epoch
* ```jd_2```: Julian day corresponding to the new epoch
* ```long_old```: Longitude for the old epoch *(radians)*
* ```lat_old```: Latitude for the old epoch *(radians)*
**/
pub fn ChangeEpochForEclipticalCoords(jd_1: f64, jd_2: f64, long_old: f64, lat_old: f64) -> (f64, f64) {
    let T = time::JulianCentury(jd_1);
    let t = (jd_2 - jd_1) / 36525.0;

    let x = T * angle::PureDegrees(0.0, 0.0, 0.000598);
    let a = (t * (angle::PureDegrees(0.0, 0.0, 47.0029) -
                   T * (angle::PureDegrees(0.0, 0.0, 0.06603) - x) +
                  t * ((angle::PureDegrees(0.0, 0.0, -0.03302) + x) +
                       t * angle::PureDegrees(0.0, 0.0, 0.00006)))).to_radians();

    let b = (174.876384 + T * (angle::PureDegrees(0.0, 0.0, 3289.4789) +
                              T * angle::PureDegrees(0.0, 0.0, 0.60622)) -
             t * ((angle::PureDegrees(0.0, 0.0, 869.8089) +
                   T * angle::PureDegrees(0.0, 0.0, 0.50491)) -
                  t * angle::PureDegrees(0.0, 0.0, 0.03536))).to_radians();

    let y = T * angle::PureDegrees(0.0, 0.0, 0.000042);
    let c = (t * (angle::PureDegrees(0.0, 0.0, 5029.0966) +
                  T * (angle::PureDegrees(0.0, 0.0, 2.22226) - y) +
                 t * ((angle::PureDegrees(0.0, 0.0, 1.11113) - y) -
                      t * angle::PureDegrees(0.0, 0.0, 0.000006)))).to_radians();

    let sin_lat_old = lat_old.sin();
    let cos_lat_old = lat_old.cos();
    let sin_b_minus_long_old = (b - long_old).sin();
    let cos_a = a.cos();
    let sin_a = a.sin();

    let A = cos_a*cos_lat_old*sin_b_minus_long_old - sin_a*sin_lat_old;
    let B = cos_lat_old*(b - long_old).cos();
    let C = cos_a*sin_lat_old + sin_a*cos_lat_old*(b - long_old).sin();

    (c + b - A.atan2(B), C.asin())
}

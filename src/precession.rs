use time;
use angle;

/**
Returns the stellar equatorial coordinates for a different equinox

# Return variables

Returns the equatorial coordinates (in radians) of a star for a
different equinox

```change_epoch_for_equa() -> (new_right_ascension, new_declination)```

# Arguments

* ```jd_1```: Julian day corresponding to the old epoch
* ```jd_2```: Julian day corresponding to the new epoch
* ```asc_old```: Right ascension for the old epoch (in radians)
* ```dec_old```: Declination for the old epoch (in radians)

**/

pub fn change_epoch_for_equa(jd_1: f64, jd_2: f64, asc_old: f64, dec_old: f64) -> (f64, f64) {
    let T = time::julian_century(jd_1);
    let t = (jd_2 - jd_1) / 36525.0;

    let x = t * (angle::pure_degrees(0.0, 0.0, 2306.2181) +
                  T * (angle::pure_degrees(0.0, 0.0, 1.39656) -
                       T*angle::pure_degrees(0.0, 0.0, 0.000139)));
    let a = (x + t*t*((angle::pure_degrees(0.0, 0.0, 0.30188) -
                       T*angle::pure_degrees(0.0, 0.0, 0.000344)) +
                      t*angle::pure_degrees(0.0, 0.0, 0.017998))).to_radians();

    let b = (x + t*t*((angle::pure_degrees(0.0, 0.0, 1.09468) -
                       T*angle::pure_degrees(0.0, 0.0, 0.000066)) +
                      t*angle::pure_degrees(0.0, 0.0, 0.018203))).to_radians();

    let y = T * angle::pure_degrees(0.0, 0.0, 0.000217);
    let c = (t * (angle::pure_degrees(0.0, 0.0, 2004.3109) +
                   T * (angle::pure_degrees(0.0, 0.0, 0.8533) - y) -
                  t * ((angle::pure_degrees(0.0, 0.0, 0.42665) + y) +
                       t*angle::pure_degrees(0.0, 0.0, 0.041833)))).to_radians();

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

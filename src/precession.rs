use time;
use angle;

/**
Returns the **annual precession** in equatorial coordinates at a certain epoch

# Returns

```(ann_precess_asc, ann_precess_dec)```

* ```ann_precess_asc```: Annual precession in right ascension for the new
                         epoch *(radians)*
* ```ann_precess_dec```: Annual precession in declination for the new
                         epoch *(radians)*

# Arguments

* ```asc```: Right ascension for the old epoch *(radians)*
* ```dec```: Declination for the old epoch *(radians)*; should not
             be too close to the poles
* ```JD```: Julian (Ephemeris) day corresponding to the new epoch;
            should not be more than a few hundred years away from
            the year 2000 AD.
**/
pub fn AnnualPrecess(asc: f64, dec: f64, JD: f64) -> (f64, f64) {
    let JC = time::JulCent(JD);
    let m = (angle::DegFrmHMS(0, 0, 3.07496) + JC*angle::DegFrmHMS(0, 0, 0.00186)).to_radians();
    let n = (angle::DegFrmHMS(0, 0, 1.33621) - JC*angle::DegFrmHMS(0, 0, 0.00057)).to_radians();

    (m + n*asc.sin()*dec.tan(), n*asc.cos())
}

/**
Returns **equatorial coordinates** for a **different epoch**

# Returns

```(new_asc, new_dec)```

* ```new_asc```: Right ascension for the new epoch *(radians)*
* ```new_dec```: Declination for the new epoch *(radians)*

# Arguments

* ```JD1```: Julian (Ephemeris) day corresponding to the old epoch
* ```JD2```: Julian (Ephemeris) day corresponding to the new epoch
* ```old_asc```: Right ascension for the old epoch *(radians)*,
                 referred to the FK5 system
* ```old_dec```: Declination for the old epoch *(radians)*,
                 referred to the FK5 system
**/
pub fn ChangeEpochEqCoords(old_asc: f64, old_dec: f64, JD1: f64, JD2: f64) -> (f64, f64) {
    let T = time::JulCent(JD1);
    let t = (JD2 - JD1) / 36525.0;

    let x = t * (angle::DegFrmDMS(0, 0, 2306.2181) +
                  T * (angle::DegFrmDMS(0, 0, 1.39656) -
                       T*angle::DegFrmDMS(0, 0, 0.000139)));
    let xi = (x + t*t*((angle::DegFrmDMS(0, 0, 0.30188) -
                       T*angle::DegFrmDMS(0, 0, 0.000344)) +
                      t*angle::DegFrmDMS(0, 0, 0.017998))).to_radians();

    let zeta = (x + t*t*((angle::DegFrmDMS(0, 0, 1.09468) -
                       T*angle::DegFrmDMS(0, 0, 0.000066)) +
                      t*angle::DegFrmDMS(0, 0, 0.018203))).to_radians();

    let y = T * angle::DegFrmDMS(0, 0, 0.000217);
    let theta = (t * (angle::DegFrmDMS(0, 0, 2004.3109) +
                   T * (angle::DegFrmDMS(0, 0, 0.8533) - y) -
                  t * ((angle::DegFrmDMS(0, 0, 0.42665) + y) +
                       t*angle::DegFrmDMS(0, 0, 0.041833)))).to_radians();

    let A = old_dec.cos() * (old_asc + xi).sin();
    let B = theta.cos()*old_dec.cos()*(old_asc + xi).cos() - theta.sin()*old_dec.sin();
    let C = theta.sin()*old_dec.cos()*(old_asc + xi).cos() + theta.cos()*old_dec.sin();

    (A.atan2(B) + zeta, C.asin())
}

/**
Returns **ecliptical coordinates** for a **different epoch**

# Returns

```(new_long, new_lat)```

* ```new_long```: Longitude for the new epoch *(radians)*
* ```new_lat```: Latitude for the new epoch *(radians)*

# Arguments

* ```JD1```: Julian (Ephemeris) day corresponding to the old epoch
* ```JD2```: Julian (Ephemeris) day corresponding to the new epoch
* ```old_long```: Longitude for the old epoch *(radians)*
* ```old_lat```: Latitude for the old epoch *(radians)*
**/
pub fn ChangeEpochEclCoords(old_long: f64, old_lat: f64, JD1: f64, JD2: f64) -> (f64, f64) {
    let T = time::JulCent(JD1);
    let t = (JD2 - JD1) / 36525.0;

    let x = T * angle::DegFrmDMS(0, 0, 0.000598);
    let a = (t * (angle::DegFrmDMS(0, 0, 47.0029) -
                   T * (angle::DegFrmDMS(0, 0, 0.06603) - x) +
                  t * ((angle::DegFrmDMS(0, 0, -0.03302) + x) +
                       t * angle::DegFrmDMS(0, 0, 0.00006)))).to_radians();

    let b = (174.876384 + T * (angle::DegFrmDMS(0, 0, 3289.4789) +
                              T * angle::DegFrmDMS(0, 0, 0.60622)) -
             t * ((angle::DegFrmDMS(0, 0, 869.8089) +
                   T * angle::DegFrmDMS(0, 0, 0.50491)) -
                  t * angle::DegFrmDMS(0, 0, 0.03536))).to_radians();

    let y = T * angle::DegFrmDMS(0, 0, 0.000042);
    let c = (t * (angle::DegFrmDMS(0, 0, 5029.0966) +
                  T * (angle::DegFrmDMS(0, 0, 2.22226) - y) +
                 t * ((angle::DegFrmDMS(0, 0, 1.11113) - y) -
                      t * angle::DegFrmDMS(0, 0, 0.000006)))).to_radians();

    let sin_old_lat = old_lat.sin();
    let cos_old_lat = old_lat.cos();
    let sin_b_minus_old_long = (b - old_long).sin();
    let cos_a = a.cos();
    let sin_a = a.sin();

    let A = cos_a*cos_old_lat*sin_b_minus_old_long - sin_a*sin_old_lat;
    let B = cos_old_lat*(b - old_long).cos();
    let C = cos_a*sin_old_lat + sin_a*cos_old_lat*(b - old_long).sin();

    (c + b - A.atan2(B), C.asin())
}

/**
Returns **equatorial coordinates**
for a **different epoch**

# Returns

```(new_asc, new_dec)```

* ```new_asc```: Right ascension for the new epoch *(radians)*
* ```new_dec```: Declination for the new epoch *(radians)*

# Arguments

* ```JD1```: Julian (Ephemeris) day corresponding to the old epoch
* ```JD2```: Julian (Ephemeris) day corresponding to the new epoch
* ```old_asc```: Right ascension for the old epoch *(radians)*,
                 referred to the FK4 system
* ```old_dec```: Declination for the old epoch *(radians)*,
                 referred to the FK4 system
**/
pub fn ChangeEpochEqCoords_FK4(old_asc: f64, old_dec: f64, JD1: f64, JD2: f64) -> (f64, f64) {
    let T = (JD1 - 2415020.3135) / 36524.2199;
    let t = (JD2 - JD1) / 36524.2199;

    let xi = t*(angle::DegFrmDMS(0, 0, 2304.25) + T*angle::DegFrmDMS(0, 0, 1.396)
                + t*(angle::DegFrmDMS(0, 0, 0.302) + t*angle::DegFrmDMS(0, 0, 0.018)));

    let zeta = xi + t*t*(angle::DegFrmDMS(0, 0, 0.791) + t*angle::DegFrmDMS(0, 0, 0.001));

    let theta = t*(angle::DegFrmDMS(0, 0, 2004.682) - T*angle::DegFrmDMS(0, 0, 0.853)
                   - t*(angle::DegFrmDMS(0, 0, 0.426) + t*angle::DegFrmDMS(0, 0, 0.042)));


    let A = old_dec.cos() * (old_asc + xi).sin();
    let B = theta.cos()*old_dec.cos()*(old_asc + xi).cos() - theta.sin()*old_dec.sin();
    let C = theta.sin()*old_dec.cos()*(old_asc + xi).cos() + theta.cos()*old_dec.sin();

    (A.atan2(B) + zeta, C.asin())
}

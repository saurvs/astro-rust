use angle;
use std;
use time;

/**
Returns the annual precession in equatorial coordinates
towards a new epoch

# Returns

```(ann_precess_asc, ann_precess_dec)```

* ```ann_precess_asc```: Annual precession in right ascension from the
                         new epoch *| in radians*
* ```ann_precess_dec```: Annual precession in declination from the
                         new epoch *| in radians*

In the case of a star, the precession returned here does not take
into account the annual proper motion of the star. The annual proper
motion in right ascension and declination must simply be added to the
corresponding annual precession returned in order to reduce the equatorial
coordinates to the new epoch.

# Arguments

* ```asc```: Right ascension for the old epoch *| in radians*
* ```dec```: Declination for the old epoch *| in radians*; should not
             be too close to the celestial poles
* ```JD```: Julian (Ephemeris) day corresponding to the new epoch;
            should not be more than a few hundred years away from
            the year 2000 AD.
**/
pub fn AnnualPrecess(asc: f64, dec: f64, JD: f64) -> (f64, f64) {

    let JC = time::JulCent(JD);

    let m = (     angle::DegFrmHMS(0, 0, 3.07496)
             + JC*angle::DegFrmHMS(0, 0, 0.00186)).to_radians();
    let n = (     angle::DegFrmHMS(0, 0, 1.33621)
             - JC*angle::DegFrmHMS(0, 0, 0.00057)).to_radians();

    (m + n*asc.sin()*dec.tan(), n*asc.cos())

}

/**
Returns equatorial coordinates reduced to a different epoch

# Returns

```(new_asc, new_dec)```

* ```new_asc```: Right ascension for the new epoch *| in radians*
* ```new_dec```: Declination for the new epoch *| in radians*

# Arguments

* ```old_asc```: Right ascension for the old epoch *| in radians*,
                 referred to the FK5 system
* ```old_dec```: Declination for the old epoch *| in radians*,
                 referred to the FK5 system
* ```JD1```: Julian (Ephemeris) day corresponding to the old epoch
* ```JD2```: Julian (Ephemeris) day corresponding to the new epoch
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
Returns equatorial coordinates, from coordinates referred to the
FK4 system, reduced to a different epoch

# Returns

```(new_asc, new_dec)```

* ```new_asc```: Right ascension for the new epoch *| in radians*
* ```new_dec```: Declination for the new epoch *| in radians*

# Arguments

* ```old_asc```: Right ascension for the old epoch *| in radians*,
                 referred to the FK4 system
* ```old_dec```: Declination for the old epoch *| in radians*,
                 referred to the FK4 system
* ```JD1```: Julian (Ephemeris) day corresponding to the old epoch
* ```JD2```: Julian (Ephemeris) day corresponding to the new epoch
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

/**
Returns ecliptical coordinates reduced to a different epoch

# Returns

```(new_long, new_lat)```

* ```new_long```: Longitude for the new epoch *| in radians*
* ```new_lat```: Latitude for the new epoch *| in radians*

# Arguments

* ```old_long```: Longitude for the old epoch *| in radians*
* ```old_lat```: Latitude for the old epoch *| in radians*
* ```JD1```: Julian (Ephemeris) day corresponding to the old epoch
* ```JD2```: Julian (Ephemeris) day corresponding to the new epoch
**/
pub fn ChangeEpochEclCoords(old_long: f64, old_lat: f64, JD1: f64, JD2: f64) -> (f64, f64) {
    let T = time::JulCent(JD1);
    let t = (JD2 - JD1) / 36525.0;

    let (nu, Pi, rho) = AnglesForEclChange(t, T);

    let A = nu.cos()*old_lat.cos()*(Pi - old_long).sin() - nu.sin()*old_lat.sin();
    let B = old_lat.cos()*(Pi - old_long).cos();
    let C = nu.cos()*old_lat.sin() + nu.sin()*old_lat.cos()*(Pi - old_long).sin();

    (rho + Pi - A.atan2(B), C.asin())
}

fn AnglesForEclChange(t: f64, T: f64) -> (f64, f64, f64) {
    let x = T * angle::DegFrmDMS(0, 0, 0.000598);
    let nu = (t * (angle::DegFrmDMS(0, 0, 47.0029) -
                   T * (angle::DegFrmDMS(0, 0, 0.06603) - x) +
                  t * ((angle::DegFrmDMS(0, 0, -0.03302) + x) +
                       t * angle::DegFrmDMS(0, 0, 0.00006)))).to_radians();

    let Pi = (174.876384 + T * (angle::DegFrmDMS(0, 0, 3289.4789) +
                              T * angle::DegFrmDMS(0, 0, 0.60622)) -
             t * ((angle::DegFrmDMS(0, 0, 869.8089) +
                   T * angle::DegFrmDMS(0, 0, 0.50491)) -
                  t * angle::DegFrmDMS(0, 0, 0.03536))).to_radians();

    let y = T * angle::DegFrmDMS(0, 0, 0.000042);
    let rho = (t * (angle::DegFrmDMS(0, 0, 5029.0966) +
                  T * (angle::DegFrmDMS(0, 0, 2.22226) - y) +
                 t * ((angle::DegFrmDMS(0, 0, 1.11113) - y) -
                      t * angle::DegFrmDMS(0, 0, 0.000006)))).to_radians();

    (nu, Pi, rho)
}

/**
Returns orbital elements reduced to a different equinox

# Returns

```(new_inc, new_arg_perih, new_long_ascend_node)```

* ```new_inc```: Inclination for the new equinox *| in radians*
* ```new_arg_perih```: Argument of perihelion for the new equinox *| in radians*
* ```new_long_ascend_node```: Longitude of ascending node for the new equinox *| in radians*

# Arguments

* ```old_inc```: Inclination for the old equinox *| in radians*
* ```old_arg_perih```: Argument of perihelion for the old equinox *| in radians*
* ```old_long_ascend_node```: Longitude of ascending node for the old equinox *| in radians*
* ```JD1```: Julian (Ephemeris) day corresponding to the old equinox
* ```JD2```: Julian (Ephemeris) day corresponding to the new equinox
**/
pub fn ChangeOrbElements(old_inc: f64, old_arg_perih: f64, old_long_ascend_node: f64,
                         JD1: f64, JD2: f64) -> (f64, f64, f64) {

    let T = time::JulCent(JD1);
    let t = (JD2 - JD1) / 36525.0;

    let (nu, Pi, rho) = AnglesForEclChange(t, T);

    let mut new_inc;
    let mut new_long_ascend_node;
    let phi = Pi + rho;

    if old_inc == 0.0 {
        new_inc = nu;
        new_long_ascend_node = phi + std::f64::consts::PI;
    }
    else {
        let A = old_inc.sin() * (old_long_ascend_node - Pi).sin();println!("A = {:?}", A);
        let B = - nu.sin() * old_inc.cos()
                + nu.cos() * old_inc.sin() * (old_long_ascend_node - Pi).cos();

        new_inc = (A*A + B*B).sqrt().asin();
        new_long_ascend_node = phi + A.atan2(B);
    }

    let delta_w = (-nu.sin()*(old_long_ascend_node - Pi).sin()/new_inc.sin()).asin();

    (new_inc, old_arg_perih + delta_w, new_long_ascend_node)

}

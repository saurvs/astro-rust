//! The Earth's moon

use angle;
use coords;
use time;

/**
Returns the **equatorial horizontal parallax** of the Moon

# Returns

* `eq_hz_parllx`: Equatorial horizontal parallax of the Moon *| in radians*

# Arguments

* `earth_moon_dist`: Earth-Moon distance *| in kilometers*
**/
pub fn eq_hz_parllx(earth_moon_dist: f64) -> f64 {
    (6378.14 / earth_moon_dist).asin()
}

/**
Returns the **equatorial semidiameter** of the Moon

# Returns

* `eq_semidiameter`: Geocentric equatorial semidiameter
                         *| in radians per kilometers*

# Arguments

* `earth_moon_dist`: Earth-Moon distance *| in kilometers*
**/
pub fn semidiameter(earth_moon_dist: f64) -> f64 {
    0.272481 * eq_hz_parllx(earth_moon_dist).sin()
}

/**
Returns the **inclination** of the **mean lunar equator** with the
ecliptic

# Returns

* `inc`: Inclination of the mean lunar equator with the
             ecliptic *| in radians*
**/
pub fn inc_of_mn_lunar_eq() -> f64 {
    angle::deg_frm_dms(1, 32, 32.7).to_radians()
}

fn A(mean_geocen_moon_long: f64, app_geocen_moon_lat: f64,
     longitude_of_mean_ascen_node: f64) -> f64 {
    let I = inc_of_mn_lunar_eq();
    let W = mean_geocen_moon_long - longitude_of_mean_ascen_node;

    (   W.sin() * app_geocen_moon_lat.cos() * I.cos()
      - app_geocen_moon_lat.sin() * I.sin()
    ).atan2(W.cos() * app_geocen_moon_lat.cos())
}

fn F(JC: f64) -> f64 {
    angle::limit_to_360(Horner_eval!(
             JC,
             93.272095,
             483202.0175233,
             -0.0036539,
             -1.0/3526000.0,
             1.0/863310000.0
    )).to_radians()
}

fn E(JC: f64) -> f64 {
    1.0 - JC*(0.002_516 + JC*0.000_0074)
}

fn DMM1(JC: f64) -> (f64, f64, f64) {
    let D = angle::limit_to_360(Horner_eval!(
             JC,
             297.8501921,
             445267.1114034,
             -0.0018819,
             1.0/545868.0,
             -1.0/113065000.0
    )).to_radians();

    let M = angle::limit_to_360(Horner_eval!(
             JC,
             357.5291092,
             35999.0502909,
             -0.0001536,
             1.0/24490000.0
    )).to_radians();

    let M1 = angle::limit_to_360(Horner_eval!(
             JC,
             134.9633964,
             477198.8675055,
             0.0087414,
             1.0/69699.0,
             -1.0/14712000.0
    )).to_radians();

    (D, M, M1)
}

fn rho_sig(D: f64, M: f64, M1: f64, F: f64) -> (f64, f64) {
    let x2F = 2.0 * F;
    let x2F2D = x2F - 2.0*D;
    (( - 0.02752 * M1.cos()
                - 0.02245 * F.sin()
                + 0.00684 * (M1 - x2F).cos()
                - 0.00293 * (x2F).cos()
                - 0.00085 * (x2F2D).cos()
                - 0.00054 * (M1 - 2.0*D).cos()
                - 0.0002  * (  (M1 + F).sin()
                             + (M1 + x2F).cos()
                             + (M1 - F).cos()
                            )
                + 0.00014 * (M1 + x2F2D).cos()
              ).to_radians(),
    ( - 0.02816 * M1.sin()
                + 0.02244 * F.cos()
                - 0.00682 * (M1 - x2F).sin()
                - 0.00279 * (x2F).sin()
                - 0.00083 * (x2F2D).sin()
                + 0.00069 * (M1 - 2.0*D).sin()
                + 0.0004  * (M1 + F).cos()
                - 0.00025 * (2.0 * M1).sin()
                - 0.00023 * (M1 + x2F).sin()
                + 0.0002  * (M1 - F).cos()
                + 0.00019 * (M1 - F).sin()
                + 0.00013 * (M1 + x2F2D).sin()
                - 0.0001  * (M1 - 3.0*F).cos()
              ).to_radians()
    )
}

/**
Returns the **optical librations** of the Moon in longitude and
latitude

# Returns

`(optical_libr_in_long, optical_libr_in_lat)`

* `optical_libr_in_long`: The optical libration in
longitude *| in radians*
* `optical_libr_in_lat`: The optical libration in
latitude *| in radians*

# Arguments

* `mn_ecl_long_moon`: The mean ecliptic longitude of the
Moon *| in radians*, i.e, *without* the correction for nutation
* `apprnt_ecl_lat_moon`: The apparent ecliptic latitude of the
Moon *| in radians*, i.e, *with* the correction for nutation
* `JD`: Julian (Ephemeris) day
**/
pub fn optical_libr(mn_ecl_long_moon: f64, apprnt_ecl_lat_moon: f64,
                         JD: f64) -> (f64, f64) {
    let JC = time::julian_cent(JD);
    let F = F(JC);
    let I = inc_of_mn_lunar_eq();

    let long_of_mn_ascend_node = mn_ascend_node(JC);
    let W = mn_ecl_long_moon - long_of_mn_ascend_node;

    let A = A(mn_ecl_long_moon, apprnt_ecl_lat_moon,
              long_of_mn_ascend_node);

    let b1 = ( - W.sin() * apprnt_ecl_lat_moon.cos() * I.sin()
               - apprnt_ecl_lat_moon.sin() * I.cos()
             ).asin();

    (angle::limit_to_360((A - F).to_degrees()).to_radians(),
     b1)
}

/**
Returns the **physical librations** of the Moon in longitude and
latitude

# Returns

`(physical_libr_in_long, physical_libr_in_lat)`

* `physical_libr_in_long`: The physical libration in
longitude *| in radians*
* `physical_libr_in_lat`: The physical libration in
latitude *| in radians*

# Arguments

* `mn_ecl_long_moon`: The mean ecliptic longitude of the
Moon *| in radians*, i.e, *without* the correction for nutation
* `apprnt_ecl_lat_moon`: The apparent ecliptic latitude of the
Moon *| in radians*, i.e, *with* the correction for nutation
* `optical_lib_lat`: The optical libration in latitude *| in radians*
* `JD`: Julian (Ephemeris) day
**/
pub fn physical_libr(mn_ecl_long_moon: f64, apprnt_ecl_lat_moon: f64,
                          optical_lib_lat: f64, JD: f64) -> (f64, f64) {
    let JC = time::julian_cent(JD);
    let K1 = (119.75 + 131.849*JC).to_radians();
    let K2 = (72.56 + 20.186*JC).to_radians();

    let long_of_mn_ascend_node = mn_ascend_node(JC);
    let (D, M, M1) = DMM1(JC);
    let F = F(JC);
    let E = E(JC);

    let x2F = 2.0 * F;
    let x2D = 2.0 * D;
    let x2F2D = x2F - x2D;

    let (rho, sig) = rho_sig(D, M, M1, F);
    let tau = (  0.0252  * E * M.sin()
               + 0.00473 * (2.0*M1 - x2F).sin()
               - 0.00467 * M1.sin()
               + 0.00396 * K1.sin()
               + 0.00276 * (2.0*M1 - x2D).sin()
               + 0.00196 * long_of_mn_ascend_node.sin()
               - 0.00183 * (M1 - F).cos()
               + 0.00115 * (M1 - x2D).sin()
               - 0.00096 * (M1 - D).sin()
               + 0.00046 * (x2F2D).sin()
               - 0.00039 * (M1 - F).sin()
               - 0.00032 * (M1 - M - D).sin()
               + 0.00027 * (2.0*M1 - M - x2D).sin()
               + 0.00023 * K2.sin()
               - 0.00014 * (  x2D.sin()
                            - (2.0*M1 - x2F).cos()
                           )
               - 0.00012 * (  (M1 - x2F).sin()
                            + (2.0*M1).sin()
                           )
               + 0.00011 * (2.0*(M1 - M - D)).sin()
              ).to_radians();

    let W = mn_ecl_long_moon - long_of_mn_ascend_node;
    let A = A(mn_ecl_long_moon, apprnt_ecl_lat_moon,
              long_of_mn_ascend_node);

    (-tau + optical_lib_lat.tan()*(rho*A.cos() + sig*A.sin()),
      sig*A.cos() - rho*A.sin())
}

/**
Returns the **total librations** of the Moon in longitude and
latitude

# Returns

`(total_libr_in_long, total_libr_in_lat)`

* `total_libr_in_long`: The total libration in
longitude *| in radians*
* `total_libr_in_lat`: The total libration in
latitude *| in radians*

# Arguments

* `mn_ecl_long_moon`: The mean ecliptic longitude of the
Moon *| in radians*, i.e, *without* the correction for nutation
* `apprnt_ecl_lat_moon`: The apparent ecliptic latitude of the
Moon *| in radians*, i.e, *with* the correction for nutation
* `JD`: Julian (Ephemeris) day
**/
pub fn total_libr(mn_ecl_long_moon: f64, apprnt_ecl_lat_moon: f64,
                         JD: f64) -> (f64, f64) {
    let (opt_long, opt_lat) = optical_libr(mn_ecl_long_moon,
                                                apprnt_ecl_lat_moon, JD);
    let (phys_long, phys_lat) = physical_libr(mn_ecl_long_moon,
                                                apprnt_ecl_lat_moon, JD, opt_lat);

    (opt_long + phys_long,
     opt_lat + phys_lat)
}

/**
Returns the **position angle** of the **axis of rotation** of the Moon

# Returns

* `pos_angl_of_axis_of_rot`: The position angle of the axis
                                            of rotation of the Moon *| in radians*

# Arguments

* `mn_ascen_node_long`: Longitude of the mean ascending node of the Moon *| in radians*
* `total_lib_lat`: Total libration of the Moon in latitude *| in radians*
* `nut_in_long`: Nutation correction for longitude *| in radians*
* `true_oblq`: True obliquity of the ecliptic *| in radians*
* `apprnt_moon_asc`: Apparent geocentric right ascension of the Moon *| in radians*
* `JD`: Julian (Ephemeris) day
**/
pub fn pos_angl_of_axis_of_rot(mn_ascen_node_long: f64, total_lib_lat: f64,
                                     nut_in_long: f64, true_oblq: f64,
                                     apprnt_moon_asc: f64, JD: f64) -> f64 {
    let JC = time::julian_cent(JD);
    let (D, M, M1) = DMM1(JC);
    let F = F(JC);
    let (rho, sig) = rho_sig(D, M, M1, F);

    let I = inc_of_mn_lunar_eq();
    let V = mn_ascen_node_long + nut_in_long + sig/I.sin();
    let X = (I + rho).sin() * V.sin();
    let Y =   (I + rho).sin()*V.cos()*true_oblq.cos()
            - (I + rho).cos()*true_oblq.sin();
    let w = X.atan2(Y);

    ( (((X*X + Y*Y).sqrt()*(apprnt_moon_asc - w).cos()) / total_lib_lat.cos()).asin() )
}

/**
Returns the **topocentric librations** of the Moon

# Returns

`(topocen_libr_in_long, topocen_libr_in_lat, topocen_libr_in_P)`

* `topocen_libr_in_long`: Topocentric libration in longitude *| in radians*
* `topocen_libr_in_lat`: Topocentric libration in latitude *| in radians*
* `topocen_libr_in_P`: Topocentric libration in position angle of the axis
                                            of rotation *| in radians*

# Arguments

* `observer_lat`: Latitude of the observer *| in radians*
* `geocen_dec_moon`: Geocentric decation of the Moon *| in radians*
* `local_hour_angle`: Local hour angle of the Moon *| in radians*
* `geocen_horizontal_parallax_moon`: Geocentric equatorial horizontal parallax of the Moon *| in radians*
* `pos_angle_axis_of_rot`: Position angle of the axis of rotation of the Moon *| in radians*
* `total_lib_lat`: Total libration of the Moon in latitude *| in radians*
**/
pub fn topocen_libr_by_diff_corrections(observer_lat: f64, geocen_dec_moon: f64,
                                                      local_hour_angle: f64, geocen_horizontal_parallax_moon: f64,
                                                      pos_angle_axis_of_rot: f64, total_lib_lat: f64) -> (f64, f64, f64) {
    let Q = (observer_lat.cos() * local_hour_angle.sin())
            .atan2(geocen_dec_moon.cos()*observer_lat.sin() - geocen_dec_moon.sin()*observer_lat.cos()*local_hour_angle.cos());
    let z = (geocen_dec_moon.sin()*observer_lat.sin() + geocen_dec_moon.cos()*observer_lat.cos()*local_hour_angle.cos()).acos();
    let pi1 = geocen_horizontal_parallax_moon * (z.sin() + 0.0084*(2.0 * z).sin());

    let delta_l = -pi1 * (Q - pos_angle_axis_of_rot).sin() / total_lib_lat.cos();
    let delta_b = -pi1 * (Q - pos_angle_axis_of_rot).cos();
    let delta_P = delta_l*(total_lib_lat + delta_b).sin() - pi1*Q.sin()*geocen_dec_moon.tan();

    (delta_l, delta_b, delta_P)
}

/**
Returns the **geocentric ecliptic** position of the Moon,
referred to the **mean equinox of the date**

# Returns

`(moon_ecl_point, rad_vec)`

* `moon_ecl_point`: Ecliptic point of the Moon *| in radians*
* `rad_vec`: Moon-Earth distance *| in kilometers*

This function uses the partial Chapront ELP-2000/82
lunar theory, so the accruacy of:

* `moon_ecl_point.long`: is `10` arcseconds
* `moon_ecl_point.lat`: is `4` arcseconds

# Arguments

* `JD`: Julian (Ephemeris) day
**/
pub fn geocen_ecl_pos(JD: f64) -> (coords::EclPoint, f64) {
    let JC = time::julian_cent(JD);
    let (D, M, M1) = DMM1(JC);
    let F = F(JC);
    let E = E(JC);
    let L1 = angle::limit_to_360(Horner_eval!(
             JC,
             218.3164477,
             481267.88123421,
             -0.0015786,
             1.0/538841.0,
             -1.0/65194000.0
    )).to_radians();

    let A1 = angle::limit_to_360(119.75 + 131.849*JC).to_radians();
    let A2 = angle::limit_to_360(53.09 + 479264.29*JC).to_radians();
    let A3 = angle::limit_to_360(313.45 + 481266.484*JC).to_radians();

    struct lrterms(i8, i8, i8, i8, i32, i32);
    let terms_for_lr = [
        lrterms(0, 0, 1, 0, 6288774, -20905355),
        lrterms(2, 0, -1, 0, 1274027, -3699111),
        lrterms(2, 0, 0, 0, 658314, -2955968),
        lrterms(0, 0, 2, 0, 213618, -569925),
        lrterms(0, 1, 0, 0, -185116, 48888),
        lrterms(0, 0, 0, 2, -114332, -3149),
        lrterms(2, 0, -2, 0, 58793, 246158),
        lrterms(2, -1, -1, 0, 57066, -152138),
        lrterms(2, 0, 1, 0, 53322, -170733),
        lrterms(2, -1, 0, 0, 45758, -204586),
        lrterms(0, 1, -1, 0, -40923, -129620),
        lrterms(1, 0, 0, 0, -34720, 108743),
        lrterms(0, 1, 1, 0, -30383, 104755),
        lrterms(2, 0, 0, -2, 15327, 10321),
        lrterms(0, 0, 1, 2, -12528, 0),
        lrterms(0, 0, 1, -2, 10980, 79661),
        lrterms(4, 0, -1, 0, 10675, -34782),
        lrterms(0, 0, 3, 0, 10034, -23210),
        lrterms(4, 0, -2, 0, 8548, -21636),
        lrterms(2, 1, -1, 0, -7888, 24208),
        lrterms(2, 1, 0, 0, -6766, 30824),
        lrterms(1, 0, -1, 0, -5163, -8379),
        lrterms(1, 1, 0, 0, 4987, -16675),
        lrterms(2, -1, 1, 0, 4036, -12831),
        lrterms(2, 0, 2, 0, 3994, -10445),
        lrterms(4, 0, 0, 0, 3861, -11650),
        lrterms(2, 0, -3, 0, 3665, 14403),
        lrterms(0, 1, -2, 0, -2689, -7003),
        lrterms(2, 0, -1, 2, -2602, 0),
        lrterms(2, -1, -2, 0, 2390, 10056),
        lrterms(1, 0, 1, 0, -2348, 6322),
        lrterms(2, -2, 0, 0, 2236, -9884),
        lrterms(0, 1, 2, 0, -2120, 5751),
        lrterms(0, 2, 0, 0, -2069, 0),
        lrterms(2, -2, -1, 0, 2048, -4950),
        lrterms(2, 0, 1, -2, -1773, 4130),
        lrterms(2, 0, 0, 2, -1595, 0),
        lrterms(4, -1, -1, 0, 1215, -3958),
        lrterms(0, 0, 2, 2, -1110, 0),
        lrterms(3, 0, -1, 0, -892, 3258),
        lrterms(2, 1, 1, 0, -810, 2616),
        lrterms(4, -1, -2, 0, 759, -1897),
        lrterms(0, 2, -1, 0, -713, -2117),
        lrterms(2, 2, -1, 0, -700, 2354),
        lrterms(2, 1, -2, 0, 691, 0),
        lrterms(2, -1, 0, -2, 596, 0),
        lrterms(4, 0, 1, 0, 549, -1423),
        lrterms(0, 0, 4, 0, 537, -1117),
        lrterms(4, -1, 0, 0, 520, -1571),
        lrterms(1, 0, -2, 0, -487, -1739),
        lrterms(2, 1, 0, -2, -399, 0),
        lrterms(0, 0, 2, -2, -381, -4421),
        lrterms(1, 1, 1, 0, 351, 0),
        lrterms(3, 0, -2, 0, -340, 0),
        lrterms(4, 0, -3, 0, 330, 0),
        lrterms(2, -1, 2, 0, 327, 0),
        lrterms(0, 2, 1, 0, -323, 1165),
        lrterms(1, 1, -1, 0, 299, 0),
        lrterms(2, 0, 3, 0, 294, 0),
        lrterms(2, 0, -1, -2, 0, 8752),
    ];
    struct bterms(i8, i8, i8, i8, i32);
    let terms_for_b = [
        bterms(0, 0, 0, 1, 5128122),
    	bterms(0, 0, 1, 1, 280602),
    	bterms(0, 0, 1, -1, 277693),
    	bterms(2, 0, 0, -1, 173237),
    	bterms(2, 0, -1, 1, 55413),
    	bterms(2, 0, -1, -1, 46271),
    	bterms(2, 0, 0, 1, 32573),
    	bterms(0, 0, 2, 1, 17198),
    	bterms(2, 0, 1, -1, 9266),
    	bterms(0, 0, 2, -1, 8822),
    	bterms(2, -1, 0, -1, 8216),
    	bterms(2, 0, -2, -1, 4324),
    	bterms(2, 0, 1, 1, 4200),
    	bterms(2, 1, 0, -1, -3359),
    	bterms(2, -1, -1, 1, 2463),
    	bterms(2, -1, 0, 1, 2211),
    	bterms(2, -1, -1, -1, 2065),
    	bterms(0, 1, -1, -1, -1870),
    	bterms(4, 0, -1, -1, 1828),
    	bterms(0, 1, 0, 1, -1794),
    	bterms(0, 0, 0, 3, -1749),
    	bterms(0, 1, -1, 1, -1565),
    	bterms(1, 0, 0, 1, -1491),
    	bterms(0, 1, 1, 1, -1475),
    	bterms(0, 1, 1, -1, -1410),
    	bterms(0, 1, 0, -1, -1344),
    	bterms(1, 0, 0, -1, -1335),
    	bterms(0, 0, 3, 1, 1107),
    	bterms(4, 0, 0, -1, 1021),
    	bterms(4, 0, -1, 1, 833),
    	bterms(0, 0, 1, -3, 777),
    	bterms(4, 0, -2, 1, 671),
    	bterms(2, 0, 0, -3, 607),
    	bterms(2, 0, 2, -1, 596),
    	bterms(2, -1, 1, -1, 491),
    	bterms(2, 0, -2, 1, -451),
    	bterms(0, 0, 3, -1, 439),
    	bterms(2, 0, 2, 1, 422),
    	bterms(2, 0, -3, -1, 421),
    	bterms(2, 1, -1, 1, -366),
    	bterms(2, 1, 0, 1, -351),
    	bterms(4, 0, 0, 1, 331),
    	bterms(2, -1, 1, 1, 315),
    	bterms(2, -2, 0, -1, 302),
    	bterms(0, 0, 1, 3, -283),
    	bterms(2, 1, 1, -1, -229),
    	bterms(1, 1, 0, -1, 223),
    	bterms(1, 1, 0, 1, 223),
    	bterms(0, 1, -2, -1, -220),
    	bterms(2, 1, -1, -1, -220),
    	bterms(1, 0, 1, 1, -185),
    	bterms(2, -1, -2, -1, 181),
    	bterms(0, 1, 2, 1, -177),
    	bterms(4, 0, -2, -1, 176),
    	bterms(4, -1, -1, -1, 166),
    	bterms(1, 0, 1, -1, -164),
    	bterms(4, 0, 1, -1, 132),
    	bterms(1, 0, -1, -1, -119),
    	bterms(4, -1, 0, -1, 115),
    	bterms(2, -2, 0, 1, 107),
    ];

    let mut l = 0.0;
    let mut r = 0.0;
    let mut b = 0.0;

    for x in terms_for_lr.iter() {
        let arg =   (x.0 as f64) * D
                  + (x.1 as f64) * M
                  + (x.2 as f64) * M1
                  + (x.3 as f64) * F;

        let t = if      (x.1).abs() == 1 { E }
                else if (x.1).abs() == 2 { E * E }
                else                     { 1.0 };

        l += (x.4 as f64) * t * arg.sin();
        r += (x.5 as f64) * t * arg.cos();
    }

    for x in terms_for_b.iter() {
        let arg =   (x.0 as f64) * D
                  + (x.1 as f64) * M
                  + (x.2 as f64) * M1
                  + (x.3 as f64) * F;
        let mut t = (x.4 as f64) * arg.sin();

        t *= if      (x.1).abs() == 1 { E }
             else if (x.1).abs() == 2 { E * E }
             else                     { 1.0 };

        b += t;
    }

    l +=   3958.0 * A1.sin()
         + 1962.0 * (L1 - F).sin()
         + 318.0  * A2.sin();

    b += - 2235.0 * L1.sin()
         + 382.0  * A3.sin()
         + 175.0  * (  (A1 - F).sin()
                     + (A1 + F).sin()
                    )
         + 127.0  * (L1 - M1).sin()
         - 115.0  * (L1 + M1).sin();

    l = l.to_radians();
    b = b.to_radians();

    let ecl_point = coords::EclPoint {
        long: L1 + l/1000000.0,
        lat: b/1000000.0
    };

    let rad_vec = 385000.56 + r/1000.0;

    (ecl_point, rad_vec)
}

/**
Returns the **longitude** of the **mean ascending node** of the Moon

# Returns

* `long_mn_ascend_node`: Longitude of the mean ascending
                                       node *| in radians*

# Arguments

* `JC`: Julian century
**/
pub fn mn_ascend_node(JC: f64) -> f64 {
    angle::limit_to_360(
        Horner_eval!(
            JC,
            125.0445479,
            -1934.1362891,
            0.0020754,
            1.0/467441.0,
            -1.0/60616000.0
        )
    ).to_radians()
}

/**
Returns the **longitude** of the **true ascending node** of the Moon

# Returns

* `long_true_ascend_node`: Longitude of the true ascending
                                       node *| in radians*

# Arguments

* `JC`: Julian century
**/
pub fn true_ascend_node(JC: f64) -> f64 {
    let (D, M, M1) = DMM1(JC);
    let F = F(JC);
    mn_ascend_node(JC) +
    (
        - 1.4979 * (2.0*(D - F)).sin()
        - 0.15   * M.sin()
        - 0.1226 * (2.0 * D).sin()
        + 0.1176 * (2.0 * F).sin()
        - 0.0801 * (2.0*(M1 - F)).sin()
    ).to_radians()
}

/**
Returns the **longitude** of the **mean perigee** of the Moon

# Returns

* `long_mn_perigee`: Longitude of mean perigee *| in radians*

# Arguments

* `JC`: Julian century
**/
pub fn mn_perigee(JC: f64) -> f64 {
    angle::limit_to_360(
        Horner_eval!(
            JC,
            83.3532465,
            4069.0137287,
            -0.01032,
            -1.0/80053.0,
            1.0/18999000.0
        )
    ).to_radians()
}

/**
Returns the **position angle** of the Moon's **bright limb**

# Returns

* `pos_angl_of_bright_limb` The position angle of the midpoint
                                       of the illuminated limb of the Moon
                                       *| in radians*

# Arguments

* `sun_eq_point`: Equatorial coordinate of the Sun *| in radians*
* `moon_eq_point`: Equatorial coordinate of the Moon *| in radians*
**/
pub fn bright_limb(sun_eq_point: coords::EqPoint,
                                 moon_eq_point: coords::EqPoint) -> f64 {
    let a = sun_eq_point.dec.cos();
    let n = a * (sun_eq_point.asc - moon_eq_point.asc).sin();
    let d =   sun_eq_point.dec.sin()  * moon_eq_point.dec.cos()
            - moon_eq_point.dec.sin() * (sun_eq_point.asc - moon_eq_point.asc).cos() * a;

    n.atan2(d)
}

/**
Returns the **illuminated fraction** of the lunar **disk**, using equatorial coordinates

# Arguments

* `sun_eq_point`: Equatorial coordinate of the Sun *| in radians*
* `moon_eq_point`: Equatorial coordinate of the Moon *| in radians*
* `earth_moon_dist`: Earth-Moon distance
                         (in any unit, but same as that of `earth_sun_dist`)
* `earth_sun_dist`: Earth-Sun distance
                        (in any unit, but same as that of `earth_moon_dist`)
**/
pub fn illum_frac_frm_eq_coords(sun_eq_point: &coords::EqPoint,
                                               moon_eq_point: &coords::EqPoint,
                                               earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    illuminated_frac(sun_eq_point.angl_sepr(&moon_eq_point).acos(),
                         earth_moon_dist, earth_sun_dist)
}

/**
Returns the **illuminated fraction** of the lunar **disk**, using eclipctic coordinates

# Arguments

* `moon_long`: Eclipctic longitude of the Moon *| in radians*
* `moon_lat`: Eclipctic latitude of the Moon *| in radians*
* `sun_long`: Eclipctic longitude of the Sun *| in radians*
* `earth_moon_dist`: Earth-Moon distance
                         (in any unit, but same as that of `earth_sun_dist`)
* `earth_sun_dist`: Earth-Sun distance
                        (in any unit, but same as that of `earth_moon_dist`)
**/
pub fn illum_frac_frm_ecl_coords(moon_long: f64, moon_lat: f64, sun_long: f64,
                                               earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    illuminated_frac((moon_lat.cos()*(moon_long - sun_long).cos()).acos(),
                         earth_moon_dist, earth_sun_dist)
}

fn illuminated_frac(moon_geocen_elong: f64, earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    let i = (earth_sun_dist * moon_geocen_elong.sin()).atan2(earth_moon_dist - earth_sun_dist*moon_geocen_elong.cos());
    (1.0 + i.cos()) / 2.0
}

/**
Returns the **times of passage** of the Moon through the **ascending**
and **descending nodes**, close to a given date

# Returns

`(time_of_ascend_node, time_of_descend_node)`

* `time_of_ascend_node`: Time of passage through the ascending node
* `time_of_descend_node`: Time of passage through the descending node

# Arguments

`date`: The Date
**/
pub fn time_of_passage_through_nodes(date: &time::Date) -> (f64, f64) {
    let k = (time::decimal_year(&date) - 2000.05)*13.4223;
    let T = k / 1342.23;
    let k1 = (k as i32) as f64;
    let k2 = (k1 as f64) + 0.5;

    (time_of_passage_through_node(k1, T), time_of_passage_through_node(k2, T))
}

fn time_of_passage_through_node(k: f64, T: f64) -> f64 {
    let D = (183.638 + 331.73735682*k + T*T*(0.0014852 + T*(0.00000209 - T*0.00000001))).to_radians();
    let D_times_2 = 2.0 * D;
    let M = (17.4006 + 26.8203725*k + T*T*(0.0001186 + T*0.00000006)).to_radians();
    let M1 = (38.3776 + 355.52747313*k + T*T*(0.0123499 + T*(0.000014627 - T*0.000000069))).to_radians();
    let sigma = (123.9767 - 1.44098956*k + T*T*(0.0020608 + T*(0.00000214 - T*0.000000016))).to_radians();
    let P = sigma + (272.75 - T*2.3).to_radians();
    let V = (299.75 + T*(132.85 - T*0.009173)).to_radians();

    (2451565.1619 + 27.212220817*k +
     T*T*(0.0002762 + T*(0.000000021 - T*0.000000000088)) -
     0.4721 * M1.sin() -
     0.1649 * D_times_2.sin() -
     0.0868 * (D_times_2 - M1).sin() +
     0.0084 * (D_times_2 + M1).sin() -
     0.0083 * (D_times_2 - M).sin() -
     0.0039 * (D_times_2 - M - M1).sin() +
     0.0034 * (2.0*M1).sin() -
     0.0031 * (2.0 * (D - M1)).sin() +
     0.003 * (D_times_2 + M).sin() +
     0.0028 * (M - M1).sin() +
     0.0026 * M.sin() +
     0.0025 * (2.0 * D_times_2).sin() +
     0.0024 * D.sin() +
     0.0022 * (M + M1).sin() +
     0.0017 * sigma.sin() +
     0.0014 * (2.0*D_times_2 - M1).sin() +
     0.0005 * (D_times_2 + M - M1).sin() +
     0.0004 * (D_times_2 - M + M1).sin() +
     0.0003 * (-1.0 * (2.0 * (D - M)).sin() +
               (2.0*D_times_2 - M).sin() + V.sin() + P.sin()))
}

/// Represents a phase of the Moon
pub enum Phase {
    /// New Moon
    New,
    /// First Quarter
    First,
    /// Full Moon
    Full,
    /// Last Quarter
    Last
}

/**
Returns the Julian day corresponding to one of the four
**phases** of the Moon

# Returns

* `JD`: Julian day corresponding to the exact time of the phase
        that is closest to `date`

*Meeus* says the mean error in `JD` for phases between 1980 AD
and mid-2020 AD is 3.8 seconds.

# Arguments

* `date`: Date of interest, close to the phase
* `phase`: The [Phase](./enum.Phase.html)
**/
pub fn time_of_phase(date: &time::Date, phase: &Phase) -> f64 {
    let K = (((time::decimal_year(&date) - 2000.0)*12.3685) as i64) as f64;
    let k = match phase {
        &Phase::New   => K,
        &Phase::First => K + 0.25,
        &Phase::Full  => K + 0.5,
        &Phase::Last  => K + 0.75,
    };
    let T = k/1236.85;

    let mut JD =
        2451550.09766
        + k*29.530588861
        + T*Horner_eval!(
            T, 0.0,
            0.000_154_37,
            -0.000_000_15,
            0.000_000_000_73
        );
    let E = E(T);
    let M = (
        2.5534
        + k*29.105_356_7
        + T*Horner_eval!(
            T, 0.0,
            -0.000_0014,
            -0.000_000_11
        )).to_radians();
    let M1 = (
        201.5643
        + k*385.816_935_28
        + T*Horner_eval!(
            T, 0.0,
            0.010_7582,
            0.000_012_38,
            -0.000_000_058
        )).to_radians();
    let F = (
        160.7108
        + k*390.670_502_84
        + T*Horner_eval!(
            T, 0.0,
            -0.001_6118,
            -0.000_002_27,
            0.000_000_011
        )).to_radians();
    let omega = (
        124.7746
        - k*1.563_755_88
        + T*Horner_eval!(
            T, 0.0,
            0.002_0672,
            0.000_002_15
        )).to_radians();

    let A1  = (299.77 +  0.107_408*k - 0.009_173*T*T).to_radians();
    let A2  = (251.88 +  0.016_321*k).to_radians();
    let A3  = (251.83 + 26.651_886*k).to_radians();
    let A4  = (349.42 + 36.412_478*k).to_radians();
    let A5  = ( 84.66 + 18.206_239*k).to_radians();
    let A6  = (141.74 + 53.303_771*k).to_radians();
    let A7  = (207.14 +  2.453_732*k).to_radians();
    let A8  = (154.84 +  7.306_860*k).to_radians();
    let A9  = ( 34.52 + 27.261_239*k).to_radians();
    let A10 = (207.19 +  0.121_824*k).to_radians();
    let A11 = (291.34 +  1.844_379*k).to_radians();
    let A12 = (161.72 + 24.198_154*k).to_radians();
    let A13 = (239.56 + 25.513_099*k).to_radians();
    let A14 = (331.55 +  3.592_518*k).to_radians();

    let is_quarter = match phase {
        &Phase::Last  => true,
        &Phase::First => true,
        &_  => false,
    };

    if is_quarter {
        let W =
            0.003_06
            - 0.000_38*E*M.cos()
            + 0.000_26*M1.cos();
            - 0.000_02*((M1 - M).cos() - (M1 + M).cos() - (2.0*F).cos());
        JD += match phase {
            &Phase::Last  => -W,
            &Phase::First => W,
            &_  => 0.0,
        };
        let corrections = [
            [-0.62801, M1],
            [0.17172*E, M],
            [-0.01183*E, M1 + M],
            [0.00862, 2.0*M1],
            [0.00804, 2.0*F],
            [0.00454*E, M1 - M],
            [0.00204*E*E, 2.0*M],
            [-0.0018, M1 - 2.0*F],
            [-0.0007, M1 + 2.0*F],
            [-0.0004, 3.0*M1],
            [-0.00034, 2.0*M1 - M],
            [0.00032*E, M + 2.0*F],
            [0.00032*E, M - 2.0*F],
            [-0.00028*E*E, M1 + 2.0*M],
            [0.00027*E, 2.0*M1 + M],
            [-0.00017, omega],
            [-0.00005, M1 - M - 2.0*F],
            [0.00004, 2.0*(M1 + F)],
            [-0.00004, M1 + M + 2.0*F],
            [0.00004, M1 - 2.0*M],
            [0.00003, M1 + M - 2.0*F],
            [0.00003, 3.0*M],
            [0.00002, 2.0*(M1 - F)],
            [0.00002, M1 - M + 2.0*F],
            [-0.00002, 3.0*M1 + M],
        ];
        for &x in corrections.iter() { JD += x[0]*x[1].sin(); }
    }
    else {
        let is_new = match phase {
            &Phase::New   => true,
            &_ => false,
        };

        let sine_arguments = [
            M1,
            M,
            2.0*M1,
            2.0*F,
            M1 - M,
            M1 + M,
            2.0*M,
            M1 - 2.0*F,
            M1 + 2.0*F,
            2.0*M1 + M,
            3.0*M1,
            M + 2.0*F,
            M - 2.0*F,
            2.0*M1 - M,
            omega,
            M1 + 2.0*M,
            2.0*(M1 - F),
            3.0*M,
            M1 + M - 2.0*F,
            2.0*(M1 + F),
            M1 + M + 2.0*F,
            M1 - M + 2.0*F,
            M1 - M - 2.0*F,
            3.0*M1 + M,
            4.0*M1
        ];

        let multipliers = if is_new {
            [-0.4072,
            0.17241*E,
            0.01608,
            0.01039,
            0.00739*E,
            -0.00514*E,
            0.00208*E*E,
            -0.00111,
            -0.00057,
            0.00056*E,
            -0.00042,
            0.00042*E,
            0.00038*E,
            -0.00024*E,
            -0.00017,
            -0.00007,
            0.00004,
            0.00004,
            0.00003,
            0.00003,
            -0.00003,
            0.00003,
            -0.00002,
            -0.00002,
            0.00002]
        }
        else {
            [-0.40614,
            0.17302*E,
            0.01614,
            0.01043,
            0.00734*E,
            -0.00515*E,
            0.00209*E*E,
            -0.00111,
            -0.00057,
            0.00056*E,
            -0.00042,
            0.00042*E,
            0.00038*E,
            -0.00024*E,
            -0.00017,
            -0.00007,
            0.00004,
            0.00004,
            0.00003,
            0.00003,
            -0.00003,
            0.00003,
            -0.00002,
            -0.00002,
            0.00002]
        };
        for mx in multipliers.iter().zip(sine_arguments.iter()) { JD += mx.0*(mx.1).sin(); }
    }

    let additional_corrections = [
        [0.000_325, A1],
        [0.000_165, A2],
        [0.000_164, A3],
        [0.000_126, A4],
        [0.000_110, A5],
        [0.000_062, A6],
        [0.000_06, A7],
        [0.000_056, A8],
        [0.000_047, A9],
        [0.000_042, A10],
        [0.000_04, A11],
        [0.000_037, A12],
        [0.000_035, A13],
        [0.000_023, A14],
    ];
    for &x in additional_corrections.iter() { JD += x[0]*x[1].sin(); }

    JD
}

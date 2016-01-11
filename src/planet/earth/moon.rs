use angle;
use coordinates;
use time;

/**
Returns the geocentric **equatorial horizontal parallax** of the Moon

# Arguments

* ```moon_earth_distance```: Distance between the Moon and Earth *(kilometers)*
**/
pub fn EquatorialHorizontalParallax(moon_earth_distance: f64) -> f64 {
    (6378.14 / moon_earth_distance).asin()
}

/**
Returns the Moon's geocentric **equatorial semidiameter**

# Returns

* ```equatorial_semidiameter```: Semidiameter *(radians per kilometers)*

# Arguments

* ```distance_to_earth```: The moon's distance to Earth *(kilometers)*
**/
pub fn GeocentricSemidiameter(distance_to_earth: f64) -> f64 {
    0.272481 * EquatorialHorizontalParallax(distance_to_earth).sin()
}

/**
Returns the **inclination** of the **mean lunar equator** with the
**ecliptic**

# Returns

```(inclination_of_mean_lunar_equator_with_ecliptic)```

* ```inclination_of_mean_lunar_equator_with_ecliptic```: Inclination
of the mean lunar equator with the ecliptic *(radians)*
**/
pub fn InclinationOfMeanLunarEquatorWithEcliptic() -> f64 {
    angle::PureDegrees(1, 32, 32.7).to_radians()
}

fn A(mean_geocen_moon_long: f64, app_geocen_moon_lat: f64,
     longitude_of_mean_ascen_node: f64) -> f64 {
    let I = InclinationOfMeanLunarEquatorWithEcliptic();
    let W = mean_geocen_moon_long - longitude_of_mean_ascen_node;

    (   W.sin() * app_geocen_moon_lat.cos() * I.cos()
      - app_geocen_moon_lat.sin() * I.sin()
    ).atan2(W.cos() * app_geocen_moon_lat.cos())
}

fn F(JC: f64) -> f64 {
    angle::LimitedTo360(93.272095 + JC*(483202.0175233 -
                                        JC*(0.0036539 +
                                        JC*(1.0/3526000.0 -
                                        JC/863310000.0)))).to_radians()
}

fn E(JC: f64) -> f64 {
    1.0 - JC*(0.002516 + JC*0.0000074)
}

fn DMM1(JC: f64) -> (f64, f64, f64) {
    (angle::LimitedTo360(297.8501921 + JC*(445267.1114034 -
                                                  JC*(0.0018819 -
                                                  JC*(1.0/545868.0 +
                                                  JC/113065000.0)))).to_radians(),
     angle::LimitedTo360(357.5291092 + JC*(35999.0502909 -
                                                  JC*(0.0001536 -
                                                  JC/24490000.0))).to_radians(),
     angle::LimitedTo360(134.9633964 + JC*(477198.8675055 +
                                                   JC*(0.0087414 +
                                                   JC*(1.0/69699.0 -
                                                   JC/14712000.0)))).to_radians())
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
Returns the **optical librations** of the Moon in **longitude** and
**latitude**

# Returns

```(optical_libration_in_longitude, optical_libration_in_longitude)```

* ```optical_libration_in_longitude```: The optical libration in
longitude *(radians)*
* ```optical_libration_in_longitude```: The optical libration in
latitude *(radians)*

# Arguments

* ```mean_geocen_moon_long```: The mean geocentric longitude of the
Moon *(radians)*, i.e, *without* the correction for nutation
* ```app_geocen_moon_lat```: The apparent geocentric longitude of the
Moon *(radians)*, i.e, *with* the correction for nutation
* ```JED```: Julian Emphemeris day
**/
pub fn OpticalLibrations(mean_geocen_moon_long: f64, app_geocen_moon_lat: f64,
                         JED: f64) -> (f64, f64) {
    let JC = time::JulianCentury(JED);
    let F = F(JC);
    let I = InclinationOfMeanLunarEquatorWithEcliptic();

    let longitude_of_mean_ascen_node = MeanAscendingNode(JC);
    let W = mean_geocen_moon_long - longitude_of_mean_ascen_node;

    let A = A(mean_geocen_moon_long, app_geocen_moon_lat,
              longitude_of_mean_ascen_node);

    let b1 = ( - W.sin() * app_geocen_moon_lat.cos() * I.sin()
               - app_geocen_moon_lat.sin() * I.cos()
             ).asin();

    (angle::LimitedTo360((A - F).to_degrees()).to_radians(),
     b1)
}

/**
Returns the **physical librations** of the Moon in **longitude** and
**latitude**

# Returns

```(physical_libration_in_longitude, physical_libration_in_longitude)```

* ```physical_libration_in_longitude```: The physical libration in
longitude *(radians)*
* ```physical_libration_in_longitude```: The physical libration in
latitude *(radians)*

# Arguments

* ```mean_geocen_moon_long```: The mean geocentric longitude of the
Moon *(radians)*, i.e, *without* the correction for nutation
* ```app_geocen_moon_lat```: The apparent geocentric longitude of the
Moon *(radians)*, i.e, *with* the correction for nutation
* ```optical_lib_lat```: The optical libration in latitude *(radians)*
* ```JED```: Julian Emphemeris day
**/
pub fn PhysicalLibrations(mean_geocen_moon_long: f64, app_geocen_moon_lat: f64,
                          JED: f64, optical_lib_lat: f64) -> (f64, f64) {
    let JC = time::JulianCentury(JED);
    let K1 = (119.75 + 131.849*JC).to_radians();
    let K2 = (72.56 + 20.186*JC).to_radians();

    let longitude_of_mean_ascen_node = MeanAscendingNode(JC);
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
               + 0.00196 * longitude_of_mean_ascen_node.sin()
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

    let W = mean_geocen_moon_long - longitude_of_mean_ascen_node;
    let A = A(mean_geocen_moon_long, app_geocen_moon_lat,
              longitude_of_mean_ascen_node);

    (-tau + optical_lib_lat.tan()*(rho*A.cos() + sig*A.sin()),
      sig*A.cos() - rho*A.sin())
}

/**
Returns the **total librations** of the Moon in **longitude** and
**latitude**

# Returns

```(total_libration_in_longitude, total_libration_in_longitude)```

* ```total_libration_in_longitude```: The total libration in
longitude *(radians)*
* ```total_libration_in_longitude```: The total libration in
latitude *(radians)*

# Arguments

* ```mean_geocen_moon_long```: The mean geocentric longitude of the
Moon *(radians)*, i.e, *without* the correction for nutation
* ```app_geocen_moon_lat```: The apparent geocentric longitude of the
Moon *(radians)*, i.e, *with* the correction for nutation
* ```optical_lib_lat```: The optical libration in latitude *(radians)*
* ```JED```: Julian Emphemeris day
**/
pub fn TotalLibrations(mean_geocen_moon_long: f64, app_geocen_moon_lat: f64,
                         JED: f64, optical_lib_lat: f64) -> (f64, f64) {
    let (opt_long, opt_lat) = OpticalLibrations(mean_geocen_moon_long,
                                                app_geocen_moon_lat, JED);
    let (phys_long, phys_lat) = PhysicalLibrations(mean_geocen_moon_long,
                                                app_geocen_moon_lat, JED, optical_lib_lat);

    (opt_long + phys_long,
     opt_lat + phys_lat)
}

/**
Returns the **position angle** of the **axis of rotation** of the Moon

# Returns

```position_angle_of_axis_of_rotation```

* ```position_angle_of_axis_of_rotation```: The position angle of the axis
                                            of rotation of the Moon *(radians)*

# Arguments

* ```mean_ascen_node_long```: Longitude of the mean ascending node of the Moon *(radians)*
* ```total_lib_lat```: Total libration of the Moon in latitude *(radians)*
* ```nut_in_long```: Nutation correction for longitude *(radians)*
* ```true_oblq_eclip```: True obliquity of the ecliptic *(radians)*
* ```app_moon_right_ascen```: Apparent geocentric right ascension of the Moon *(radians)*
* ```JED```: Julian Emphemeris day
**/
pub fn AxisOfRotation(mean_ascen_node_long: f64, total_lib_lat: f64,
                                     nut_in_long: f64, true_oblq_eclip: f64,
                                     app_moon_right_ascen: f64, JED: f64) -> f64 {
    let JC = time::JulianCentury(JED);
    let (D, M, M1) = DMM1(JC);
    let F = F(JC);
    let (rho, sig) = rho_sig(D, M, M1, F);

    let I = InclinationOfMeanLunarEquatorWithEcliptic();
    let V = mean_ascen_node_long + nut_in_long + sig/I.sin();
    let X = (I + rho).sin() * V.sin();
    let Y =   (I + rho).sin()*V.cos()*true_oblq_eclip.cos()
            - (I + rho).cos()*true_oblq_eclip.sin();
    let w = X.atan2(Y);

    ( (((X*X + Y*Y).sqrt()*(app_moon_right_ascen - w).cos()) / total_lib_lat.cos()).asin() )
}

/**
Returns the **topocentric librations** of the Moon

# Returns

```(topocentric_libration_in_long, topocentric_libration_in_lat, topocentric_libration_in_P)```

* ```topocentric_libration_in_long```: Topocentric libration in longitude *(radians)*
* ```topocentric_libration_in_lat```: Topocentric libration in latitude *(radians)*
* ```topocentric_libration_in_P```: Topocentric libration in position angle of the axis
                                            of rotation *(radians)*

# Arguments

* ```observer_lat```: Latitude of the observer *(radians)*
* ```geocen_declin_moon```: Geocentric declination of the Moon *(radians)*
* ```local_hour_angle```: Local hour angle of the Moon *(radians)*
* ```geocen_horizontal_parallax_moon```: Geocentric equatorial horizontal parallax of the Moon *(radians)*
* ```pos_angle_axis_of_rot```: Position angle of the axis of rotation of the Moon *(radians)*
* ```total_lib_lat```: Total libration of the Moon in latitude *(radians)*
**/
pub fn TopocentricLibrationsByDifferentialCorrections(observer_lat: f64, geocen_declin_moon: f64,
                                                      local_hour_angle: f64, geocen_horizontal_parallax_moon: f64,
                                                      pos_angle_axis_of_rot: f64, total_lib_lat: f64) -> (f64, f64, f64) {
    let Q = (observer_lat.cos() * local_hour_angle.sin())
            .atan2(geocen_declin_moon.cos()*observer_lat.sin() - geocen_declin_moon.sin()*observer_lat.cos()*local_hour_angle.cos());
    let z = (geocen_declin_moon.sin()*observer_lat.sin() + geocen_declin_moon.cos()*observer_lat.cos()*local_hour_angle.cos()).acos();
    let pi1 = geocen_horizontal_parallax_moon * (z.sin() + 0.0084*(2.0 * z).sin());

    let delta_l = -pi1 * (Q - pos_angle_axis_of_rot).sin() / total_lib_lat.cos();
    let delta_b = -pi1 * (Q - pos_angle_axis_of_rot).cos();
    let delta_P = delta_l*(total_lib_lat + delta_b).sin() - pi1*Q.sin()*geocen_declin_moon.tan();

    (delta_l, delta_b, delta_P)
}

/**
Returns the **geocentric ecliptical coordinates** of the moon

# Returns

```(longitude, latitude, distance)```

* ```longitude```: Ecliptical longitude of the Moon *(radians)*
* ```latitude```: Ecliptical latitude of the Moon *(radians)*
* ```distance```: Distance between the Moon and the Earth *(kilometers)*

# Arguments

* ```JED```: Julian Emphemeris day
**/
pub fn EclipticalGeocentricCoords(JED: f64) -> (f64, f64, f64) {
    let JC = time::JulianCentury(JED);
    let (D, M, M1) = DMM1(JC);
    let F = F(JC);
    let E = E(JC);
    let L1 = angle::LimitedTo360(218.3164477 + JC*(481267.88123421 -
                                                   JC*(0.0015786 -
                                                   JC*(1.0/538841.0 +
                                                   JC/65194000.0)))).to_radians();

    let A1 = angle::LimitedTo360(119.75 + 131.849*JC).to_radians();
    let A2 = angle::LimitedTo360(53.09 + 479264.29*JC).to_radians();
    let A3 = angle::LimitedTo360(313.45 + 481266.484*JC).to_radians();

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

    (L1 + l/1000000.0,
     b/1000000.0,
     385000.56 + r/1000.0)
}

/**
Returns the longitude of the **mean ascending node** of the Moon

# Returns

```(longitude_mean_ascending_node)```

* ```longitude_mean_ascending_node```: Longitude of the mean ascending
                                       node *(radians)*

# Arguments

* ```JC```: Julian century
**/
pub fn MeanAscendingNode(JC: f64) -> f64 {
    angle::LimitedTo360(125.0445479 - JC*(1934.1362891 -
                                          JC*(0.0020754 +
                                          JC*(1.0/467441.0 -
                                          JC/60616000.0)))).to_radians()
}

/**
Returns the longitude of the **true ascending node** of the Moon

# Returns

```(longitude_true_ascending_node)```

* ```longitude_true_ascending_node```: Longitude of the true ascending
                                       node *(radians)*

# Arguments

* ```JC```: Julian century
**/
pub fn TrueAscendingNode(JC: f64) -> f64 {
    let (D, M, M1) = DMM1(JC);
    let F = F(JC);
    MeanAscendingNode(JC) +
    (
        - 1.4979 * (2.0*(D - F)).sin()
        - 0.15   * M.sin()
        - 0.1226 * (2.0 * D).sin()
        + 0.1176 * (2.0 * F).sin()
        - 0.0801 * (2.0*(M1 - F)).sin()
    ).to_radians()
}

/**
Returns the longitude of the **mean perigee** of the Moon

# Returns

```(longitude_mean_perigee)```

* ```longitude_mean_perigee```: Longitude of mean perigee *(radians)*

# Arguments

* ```JC```: Julian century
**/
pub fn MeanPerigee(JC: f64) -> f64 {
    angle::LimitedTo360(83.3532465 + JC*(4069.0137287 -
                                         JC*(0.01032 +
                                         JC*(1.0/80053.0 -
                                         JC/18999000.0)))).to_radians()
}

/**
Returns the **position angle** of the **bright limb** of the Moon

# Returns

```position_angle_of_bright_limb```

* ```position_angle_of_bright_limb```: The position angle of the midpoint
                                       of the illuminated limb of the Moon
                                       *(radians)*

# Arguments

* ```sun_equa_point```: Equatorial coordinate of the Sun *(radians)*
* ```moon_equa_point```: Equatorial coordinate of the Moon *(radians)*
**/
pub fn BrightLimb(sun_equa_point: coordinates::EquatorialPoint,
                                 moon_equa_point: coordinates::EquatorialPoint) -> f64 {
    let a = sun_equa_point.declin.cos();
    let n = a * (sun_equa_point.right_ascen - moon_equa_point.right_ascen).sin();
    let d =   sun_equa_point.declin.sin()  * moon_equa_point.declin.cos()
            - moon_equa_point.declin.sin() * (sun_equa_point.right_ascen - moon_equa_point.right_ascen).cos() * a;
    n.atan2(d)
}

/**
Returns the **illuminated fraction** of the Moon, using **equatorial coordinates**

# Arguments

* ```sun_equa_point```: Equatorial coordinate of the Sun *(radians)*
* ```moon_equa_point```: Equatorial coordinate of the Moon *(radians)*
* ```earth_moon_dist```: Distance between the Earth and it's Moon
                         (in any unit, but same as that of ```earth_sun_dist```)
* ```earth_sun_dist```: Distance between the Earth and the Sun
                        (in any unit, but same as that of ```earth_moon_dist```)
**/
pub fn IlluminatedFractionFromEquatorCoords(sun_equa_point: coordinates::EquatorialPoint,
                                               moon_equa_point: coordinates::EquatorialPoint,
                                               earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    illuminated_fraction(sun_equa_point.AngularSeparation(moon_equa_point).acos(),
                         earth_moon_dist, earth_sun_dist)
}

/**
Return the **illuminated fraction** of the Moon, using **eclipctical coordinates**

# Arguments

* ```moon_long```: Eclipctical longitude of the Moon *(radians)*
* ```moon_lat```: Eclipctical latitude of the Moon *(radians)*
* ```sun_long```: Eclipctical longitude of the Sun *(radians)*
* ```earth_moon_dist```: Distance between the Earth and it's Moon
                         (in any unit, but same as that of ```earth_sun_dist```)
* ```earth_sun_dist```: Distance between the Earth and the Sun
                        (in any unit, but same as that of ```earth_moon_dist```)
**/
pub fn IlluminatedFractionFromEclipCoords(moon_long: f64, moon_lat: f64, sun_long: f64,
                                               earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    illuminated_fraction((moon_lat.cos()*(moon_long - sun_long).cos()).acos(),
                         earth_moon_dist, earth_sun_dist)
}

fn illuminated_fraction(moon_geocen_elong: f64, earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    let i = (earth_sun_dist * moon_geocen_elong.sin()).atan2(earth_moon_dist - earth_sun_dist*moon_geocen_elong.cos());
    (1.0 + i.cos()) / 2.0
}

/**
Returns the **times of passage** of the Moon through the **ascending**
and **descending nodes**, close to a given date

# Returns

```(time_of_ascending_node, time_of_descending_node)```

* ```time_of_ascending_node```: Time of passage through the ascending node
* ```time_of_descending_node```: Time of passage through the descending node

# Arguments

```date```: The Date
**/
pub fn TimesOfPassageThroughNodes(date: time::Date) -> (f64, f64) {
    let k = (time::DecimalYear(date) - 2000.05)*13.4223;
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

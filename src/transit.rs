use coordinates;

enum CelestialBody {
    StarOrPlanet,
    Sun,
    Moon
}

pub fn func(celestial_body: CelestialBody,
            apparent_sidreal_green: f64
            observer_long: f64, observer_lat: f64
            eq_pnt_1: Equa, eq_pnt_2: f64, eq_pnt_3: f64,
            moon_eq_hz_parallax: f64) {
    let standard_alt = match celestial_body {
        CelestialBody::StarOrPlanet => angle::DegFrmDMS(-0, 34, 0.0),
        CelestialBody::Sun => angle::DegFrmDMS(-0, 50, 0.0),
        CelestialBody::Moon =>   0.7275 * moon_eq_hz_parallax
                               - angle::DegFrmDMS(-0, 34, 0.0),
    }

    let mut H0 = ( (standard_alt.sin() - observer_lat.sin()*eq_pnt_2.declin.sin())
                   / (observer_lat.cos() * declin_2.cos())
                 ).acos();
    H0 = angle::LimitedTo360(H0.to_degrees()).to_radians();

    let rad360 = 360.0_f64.to_radians();

    let mut m_0 = right_ascen_2 + observer_long - apparent_sidreal_green;
    let mut m_1 = m_0 - H0 / rad360;
    let mut m_2 = m_0 + H0 / rad360;

    let r = 360.985647_f64.to_radians();
    let sid_0 = apparent_sidreal_green + r*m_0; // for transit
    let sid_1 = apparent_sidreal_green + r*m_1; // for rising
    let sid_2 = apparent_sidreal_green + r*m_2; // for setting

    let right_ascen = 0.0;
    let declin = 0.0;

    // for transit
    let local_hour_angle_0 = coords::HrAngFrmObserverLong(sid_0,
                                                                    observer_long, right_ascen);
    m_0 -= local_hour_angle_0 / rad360;

    // for rising
    let local_hour_angle_1 = coords::HrAngFrmObserverLong(sid_1,
                                                                    observer_long, right_ascen);
    let altitude_1 = coords::Az(local_hour_angle_1, declin, observer_lat);
    m_1 +=   (altitude_1 - standard_alt)
           / (rad360 * declin.cos() * observer_lat.cos() * local_hour_angle_1.sin());

    // for setting
    let local_hour_angle_2 = coords::HrAngFrmObserverLong(sid_2,
                                                                    observer_long, right_ascen);
    let altitude_2 = coords::Az(local_hour_angle_2, declin, observer_lat);
    m_2 +=   (altitude_2 - standard_alt)
           / (rad360 * declin.cos() * observer_lat.cos() * local_hour_angle_2.sin());

}

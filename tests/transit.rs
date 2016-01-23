extern crate astro;

use astro::*;

#[test]
fn time() {
    let deltaT = time::ApproxDelT(1988, 3);println!("deltaT = {:?}", deltaT);
    let eq_point1 = coords::EqPoint{
        asc: 40.68021_f64.to_radians(),
        dec: 18.04761_f64.to_radians()
    };
    let eq_point2 = coords::EqPoint{
        asc: 41.73129_f64.to_radians(),
        dec: 18.44092_f64.to_radians()
    };
    let eq_point3 = coords::EqPoint{
        asc: 42.78204_f64.to_radians(),
        dec: 18.82742_f64.to_radians()
    };
    let geograph_point = coords::GeographPoint{
        long: 71.0833_f64.to_radians(),
        lat: 42.3333_f64.to_radians(),
    };
    let Theta0 = 177.74208_f64.to_radians();

    let (h_rise, m_rise, s_rise) = transit::time(
        &transit::TransitType::Rise,
        &transit::TransitBody::StarOrPlanet,
        &geograph_point,
        &eq_point1,
        &eq_point2,
        &eq_point3,
        Theta0,
        deltaT,
        0.0
    );
    assert_eq!((h_rise, m_rise), (12, 25));

    let (h_transit, m_transit, s_transit) = transit::time(
        &transit::TransitType::Transit,
        &transit::TransitBody::StarOrPlanet,
        &geograph_point,
        &eq_point1,
        &eq_point2,
        &eq_point3,
        Theta0,
        deltaT,
        0.0
    );
    assert_eq!((h_transit, m_transit), (19, 40));

    let (h_set, m_set, s_set) = transit::time(
        &transit::TransitType::Set,
        &transit::TransitBody::StarOrPlanet,
        &geograph_point,
        &eq_point1,
        &eq_point2,
        &eq_point3,
        Theta0,
        deltaT,
        0.0
    );
    assert_eq!((h_set, m_set), (2, 54));
}

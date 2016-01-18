extern crate astro;

use astro::*;

#[test]
fn Time() {
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

    let m_rise = transit::Time(
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
    assert_eq!(util::RoundUptoDigits(m_rise, 5), 0.51766);

    let m_transit = transit::Time(
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
    assert_eq!(util::RoundUptoDigits(m_transit, 5), 0.81980);

    let m_set = transit::Time(
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
    assert_eq!(util::RoundUptoDigits(m_set, 5), 0.12130);
}

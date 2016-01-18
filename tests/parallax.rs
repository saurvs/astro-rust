extern crate astro;

use astro::*;

#[test]
pub fn TopocenEqCoords() {

    let eq_point = coords::EqPoint{
        asc: 339.530208_f64.to_radians(),
        dec: -15.771083_f64.to_radians()
    };
    let geograph_point = coords::GeographPoint{
        long: angle::DegFrmHMS(7, 47, 27.0).to_radians(),
        lat: 33.356111_f64.to_radians(),
    };
    let topo_eq_point = parallax::TopocenEqCoords(
        &eq_point,
        0.37276,
        &geograph_point,
        1706.0,
        angle::DegFrmHMS(1, 40, 45.0).to_radians()
    );

    let (h, m1, s1) = angle::HMSFrmDeg(topo_eq_point.asc.to_degrees());
    assert_eq!((h, m1, util::RoundUptoDigits(s1, 2)), (22, 38, 8.54));

    let (d, m2, s2) = angle::DMSFrmDeg(topo_eq_point.dec.to_degrees());
    assert_eq!((d, m2, util::RoundUptoDigits(s2, 1)), (-15, -46, -30.0));

}
extern crate astro;

use astro::*;

#[test]
pub fn topocen_eq_coords() {

    let eq_point = coords::EqPoint{
        asc: 339.530208_f64.to_radians(),
        dec: -15.771083_f64.to_radians()
    };
    let geograph_point = coords::GeographPoint{
        long: angle::deg_frm_hms(7, 47, 27.0).to_radians(),
        lat: 33.356111_f64.to_radians(),
    };
    let topo_eq_point = parallax::topocen_eq_coords(
        &eq_point,
        angle::deg_frm_dms(0, 0, 23.592).to_radians(),
        &geograph_point,
        1706.0,
        angle::deg_frm_hms(1, 40, 45.0).to_radians()
    );

    let (h, m1, s1) = angle::hms_frm_deg(topo_eq_point.asc.to_degrees());
    assert_eq!((h, m1, util::round_upto_digits(s1, 2)), (22, 38, 8.54));

    let (d, m2, s2) = angle::dms_frm_deg(topo_eq_point.dec.to_degrees());
    assert_eq!((d, m2, util::round_upto_digits(s2, 1)), (-15, -46, -30.0));

}

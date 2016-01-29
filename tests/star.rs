extern crate astro;

use astro::*;

#[test]
fn eq_coords_frm_motion() {
    let (asc, dec) = star::eq_coords_frm_motion(101.286962_f64.to_radians(),
                                             -16.716108_f64.to_radians(),
                                             2.64, -0.000007773,
                                             (-0.03847/3600.0 as f64).to_radians(),
                                             (-1.2053/3600.0 as f64).to_radians(),
                                             -1000.0);

    let (h1, m1, s1) = angle::hms_frm_deg(asc.to_degrees());
    assert_eq!((h1, m1, util::round_upto_digits(s1, 2)), (6, 45, 47.16));

    let (d2, m2, s2) = angle::dms_frm_deg(dec.to_degrees());
    assert_eq!((d2, m2, util::round_upto_digits(s2, 1)), (-16, -22, -56.0));
}

extern crate astro;

use astro::*;

#[test]
fn EqCoordsFrmMotion() {
    let (asc, dec) = star::EqCoordsFrmMotion(101.286962_f64.to_radians(),
                                             -16.716108_f64.to_radians(),
                                             2.64, -0.000007773,
                                             angle::DegFrmTimeSec(-0.03847).to_radians(),
                                             angle::DegFrmArcSec(-1.2053).to_radians(),
                                             -1000.0);
                                             
    let (h1, m1, s1) = angle::HMSFrmDeg(asc.to_degrees());
    assert_eq!((h1, m1, util::RoundUptoDigits(s1, 2)), (6, 45, 47.16));

    let (d2, m2, s2) = angle::DMSFrmDeg(dec.to_degrees());
    assert_eq!((d2, m2, util::RoundUptoDigits(s2, 1)), (-16, -22, -56.0));
}

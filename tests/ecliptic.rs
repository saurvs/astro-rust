extern crate astro;

use astro::*;

#[test]
fn MnOblq_Laskar() {
    let (d, m, s) = angle::DMSFrmDeg(ecliptic::MnOblq_Laskar(2446895.5).to_degrees());
    assert_eq!((d, m, util::RoundUptoDigits(s, 3)), (23, 26, 27.407));
}

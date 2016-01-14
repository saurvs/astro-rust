extern crate astro;

use astro::*;

#[test]
fn nutation() {
    let (nut_in_long, nut_in_oblq) = nutation::Corrections(2446895.5);

    let (d1, m1, s1) = angle::DMSFrmDeg(nut_in_long.to_degrees());
    assert_eq!((d1, m1, util::RoundUptoDigits(s1, 3)), (0, 0, -3.788));

    let (d2, m2, s2) = angle::DMSFrmDeg(nut_in_oblq.to_degrees());
    assert_eq!((d2, m2, util::RoundUptoDigits(s2, 3)), (0, 0, 9.443));
}

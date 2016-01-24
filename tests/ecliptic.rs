extern crate astro;

use astro::*;

#[test]
fn mn_oblq_Laskar() {
    let (d, m, s) = angle::dms_frm_deg(ecliptic::mn_oblq_laskar(2446895.5).to_degrees());
    assert_eq!((d, m, util::round_upto_digits(s, 3)), (23, 26, 27.407));
}

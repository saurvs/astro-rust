extern crate astro;

use astro::*;

#[test]
fn nutation() {
    let (nut_in_long, nut_in_oblq) = nutation::nutation(2446895.5);

    let (d1, m1, s1) = angle::dms_frm_deg(nut_in_long.to_degrees());
    assert_eq!((d1, m1, util::round_upto_digits(s1, 3)), (0, 0, -3.788));

    let (d2, m2, s2) = angle::dms_frm_deg(nut_in_oblq.to_degrees());
    assert_eq!((d2, m2, util::round_upto_digits(s2, 3)), (0, 0, 9.443));
}

#[test]
fn nutation_in_eq_coords() {
    let d = time::Date{year: 2028, month: 11, decimal_day: 13.19, cal_type: time::CalType::Gregorian};
    let (a, b) = nutation::nutation_in_eq_coords(
                        41.5555635_f64.to_radians(), 49.3503415_f64.to_radians(),
                        angle::deg_frm_dms(0, 0, 14.861).to_radians(), angle::deg_frm_dms(0, 0, 2.705).to_radians(),
                        23.436_f64.to_radians());

    assert_eq!(util::round_upto_digits(a.to_degrees(), 7), 0.0044011);
    assert_eq!(util::round_upto_digits(b.to_degrees(), 7), 0.001727);
}

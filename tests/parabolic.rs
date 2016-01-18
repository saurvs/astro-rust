extern crate astro;

use astro::*;

#[test]
fn TruAnomAndRadVec() {
    let t_date = time::Date{year: 1998, month: 8, decimal_day: 5.0, cal_type: time::CalType::Gregorian};
    let t = time::JulDay(&t_date);

    let T_date = time::Date{year: 1998, month: 4, decimal_day: 14.4358, cal_type: time::CalType::Gregorian};
    let T = time::JulDay(&T_date);

    let (tru_anom, rad_vec) = orbit::parabolic::TruAnomAndRadVec(t, T, 1.487469);
    assert_eq!(util::RoundUptoDigits(tru_anom.to_degrees(), 5), 66.78862);
    assert_eq!(util::RoundUptoDigits(rad_vec, 6), 2.133911);
}

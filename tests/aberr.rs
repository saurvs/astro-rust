extern crate astro;

use astro::*;

#[test]
fn StellAberrInEqCoords() {
    let d = time::Date{year: 2028, month: 11, decimal_day: 13.19, cal_type: time::CalType::Gregorian};
    let (a, b) = aberr::StellAberrInEqCoords(41.0540613_f64.to_radians(), 49.2277489_f64.to_radians(), time::JulDay(&d));

    assert_eq!(util::RoundUptoDigits(a.to_degrees(), 7), 0.0083223);
    assert_eq!(util::RoundUptoDigits(b.to_degrees(), 7), 0.0018749);
}

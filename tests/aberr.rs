extern crate astro;

use astro::*;

#[test]
fn stell_aberr_in_eq_coords() {
    let d = time::Date{year: 2028, month: 11, decimal_day: 13.19, cal_type: time::CalType::Gregorian};
    let stell_eq_point = coords::EqPoint{
        asc: 41.0540613_f64.to_radians(),
        dec: 49.2277489_f64.to_radians()
    };
    let (a, b) = aberr::stell_aberr_in_eq_coords(&stell_eq_point, time::JulDay(&d));

    assert_eq!(util::RoundUptoDigits(a.to_degrees(), 7), 0.0083223);
    assert_eq!(util::RoundUptoDigits(b.to_degrees(), 7), 0.0018749);
}

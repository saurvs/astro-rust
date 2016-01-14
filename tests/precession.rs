extern crate astro;

use astro::*;

#[test]
fn AnnualPrecess() {
    let d = time::Date{year: 1978, month: 1, decimal_day: 0.0, cal_type: time::CalType::Gregorian};
    
    let (new_asc, new_dec) = precession::AnnualPrecess(
        angle::DegFrmHMS(10, 8, 22.3).to_radians(),
        angle::DegFrmDMS(11, 58, 2.0).to_radians(),
        time::JulDay(&d)
    );

    let (a, b, c) = angle::HMSFrmDeg(new_asc.to_degrees());
    assert_eq!(util::RoundUptoDigits(c, 2), util::RoundUptoDigits(3.208, 2));

    let (d, e, f) = angle::DMSFrmDeg(new_dec.to_degrees());
    assert_eq!(util::RoundUptoDigits(f, 2), -17.71);
}

#[test]
fn ChangeEpochEqCoords() {
    let (new_asc, new_dec) = precession::ChangeEpochEqCoords(
        41.054063_f64.to_radians(),
        49.22775_f64.to_radians(),
        2451545.0,
        2462088.69
    );
    assert_eq!((util::RoundUptoDigits(new_asc.to_degrees(), 6),
                util::RoundUptoDigits(new_dec.to_degrees(), 6)), (41.547214, 49.348483));
}

#[test]
fn ChangeEpochEclCoords() {
    let (new_asc, new_dec) = precession::ChangeEpochEclCoords(
        149.48194_f64.to_radians(),
        1.76549_f64.to_radians(),
        2451545.0,
        1643074.5
    );
    assert_eq!((util::RoundUptoDigits(new_asc.to_degrees(), 3),
                util::RoundUptoDigits(new_dec.to_degrees(), 3)), (118.704, 1.615));
}

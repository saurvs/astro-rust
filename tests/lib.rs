#[macro_use]
extern crate astro;

use astro::*;

fn RoundUptoDigits(number: f64, frac_digits: i32) -> f64 {
    let d = 10_f64.powi(frac_digits);
    (number * d).round() / d
}

#[test]
fn SidrealTime() {
    let (h1, m1, s1) = angle::HMSFrmDeg(time::MnSidr(2446895.5).to_degrees());
    assert_eq!((h1, m1, RoundUptoDigits(s1, 4)), (13, 10, 46.3668));

    let (h2, m2, s2) = angle::HMSFrmDeg(AppSidr!(2446895.5).to_degrees());
    assert_eq!((h2, m2, RoundUptoDigits(s2, 4)), (13, 10, 46.1351));
}

#[test]
fn nutation() {
    let (nut_in_long, nut_in_oblq) = nutation::Corrections(2446895.5);

    let (d1, m1, s1) = angle::DMSFrmDeg(nut_in_long.to_degrees());
    assert_eq!((d1, m1, RoundUptoDigits(s1, 3)), (0, 0, -3.788));

    let (d2, m2, s2) = angle::DMSFrmDeg(nut_in_oblq.to_degrees());
    assert_eq!((d2, m2, RoundUptoDigits(s2, 3)), (0, 0, 9.443));
}

#[test]
fn MnOblq() {
    let (d, m, s) = angle::DMSFrmDeg(ecliptic::MnOblq(2446895.5).to_degrees());
    assert_eq!((d, m, RoundUptoDigits(s, 3)), (23, 26, 27.407));
}

#[test]
fn JulDay() {

    let mut date = time::Date{year: 0, month: 0,
                              decimal_day: 0.0,
                              cal_type: time::CalType::Gregorian};

    struct test_date(i32, u8, f64, f64, time::CalType);
    let gregorian_dates = [test_date(1957, 10, 4.81, 2436116.31, time::CalType::Gregorian),
                           test_date(2000, 1, 1.5, 2451545.0, time::CalType::Gregorian),
                           test_date(1999, 1, 1.0, 2451179.5, time::CalType::Gregorian),
                           test_date(1987, 1, 27.0, 2446822.5, time::CalType::Gregorian),
                           test_date(1987, 6, 19.5, 2446966.0, time::CalType::Gregorian),
                           test_date(1988, 1, 27.0, 2447187.5, time::CalType::Gregorian),
                           test_date(1988, 6, 19.5, 2447332.0, time::CalType::Gregorian),
                           test_date(1900, 1, 1.0, 2415020.5, time::CalType::Gregorian),
                           test_date(1600, 1, 1.0, 2305447.5, time::CalType::Gregorian),
                           test_date(1600, 12, 31.0, 2305812.5, time::CalType::Gregorian),
                           test_date(837, 4, 10.3, 2026871.8, time::CalType::Julian),
                           test_date(-123, 12, 31.0, 1676496.5, time::CalType::Julian),
                           test_date(-122, 1, 1.0, 1676497.5, time::CalType::Julian),
                           test_date(-1000, 7, 12.5, 1356001.0, time::CalType::Julian),
                           test_date(-1000, 2, 29.0, 1355866.5, time::CalType::Julian),
                           test_date(-1001, 8, 17.9, 1355671.4, time::CalType::Julian),
                           test_date(-4712, 1, 1.5, 0.0, time::CalType::Julian)];

    for date_fields in gregorian_dates.iter() {
        date.year = date_fields.0;
        date.month = date_fields.1;
        date.decimal_day = date_fields.2;
        date.cal_type = match date_fields.4 {
            time::CalType::Gregorian => time::CalType::Gregorian,
            time::CalType::Julian    => time::CalType::Julian,
        };
        assert_eq!(time::JulDay(&date), date_fields.3);
    }

}

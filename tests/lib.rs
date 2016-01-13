extern crate astro;

use astro::*;

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

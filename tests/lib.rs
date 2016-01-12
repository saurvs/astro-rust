extern crate astro;

use astro::*;

#[test]
fn JulDay() {

    let mut date = time::Date{year       : 1957,
                              month      : 10,
                              decimal_day: 4.81,
                              cal_type   : time::CalType::Gregorian};
    assert_eq!(time::JulDay(&date), 2436116.31);

    date.year = 333;
    date.month = 1;
    date.decimal_day = time::DecimalDay(&time::DayOfMonth{day: 27, hr : 12, min: 0, sec: 0.0});
    date.cal_type = time::CalType::Julian;
    assert_eq!(time::JulDay(&date), 1842713.0);

    date.cal_type = time::CalType::Gregorian;

    date.year = 2000; date.month = 1; date.decimal_day = 1.5;
    assert_eq!(time::JulDay(&date), 2451545.0);

    date.year = 1999; date.month = 1; date.decimal_day = 1.0;
    assert_eq!(time::JulDay(&date), 2451179.5);

    date.year = 1987; date.month = 1; date.decimal_day = 27.0;
    assert_eq!(time::JulDay(&date), 2446822.5);

    date.year = 1987; date.month = 6; date.decimal_day = 19.5;
    assert_eq!(time::JulDay(&date), 2446966.0);

    date.year = 1988; date.month = 1; date.decimal_day = 27.0;
    assert_eq!(time::JulDay(&date), 2447187.5);

    date.year = 1988; date.month = 6; date.decimal_day = 19.5;
    assert_eq!(time::JulDay(&date), 2447332.0);

    date.year = 1900; date.month = 1; date.decimal_day = 1.0;
    assert_eq!(time::JulDay(&date), 2415020.5);

    date.year = 1600; date.month = 1; date.decimal_day = 1.0;
    assert_eq!(time::JulDay(&date), 2305447.5);

    date.year = 1600; date.month = 12; date.decimal_day = 31.0;
    assert_eq!(time::JulDay(&date), 2305812.5);

    date.cal_type = time::CalType::Julian;

    date.year = 837; date.month = 4; date.decimal_day = 10.3;
    assert_eq!(time::JulDay(&date), 2026871.8);

    date.year = -123; date.month = 12; date.decimal_day = 31.0;
    assert_eq!(time::JulDay(&date), 1676496.5);

    date.year = -122; date.month = 1; date.decimal_day = 1.0;
    assert_eq!(time::JulDay(&date), 1676497.5);

    date.year = -1000; date.month = 7; date.decimal_day = 12.5;
    assert_eq!(time::JulDay(&date), 1356001.0);

    date.year = -1000; date.month = 2; date.decimal_day = 29.0;
    assert_eq!(time::JulDay(&date), 1355866.5);

    date.year = -1001; date.month = 8; date.decimal_day = 17.9;
    assert_eq!(time::JulDay(&date), 1355671.4);

    date.year = -4712; date.month = 1; date.decimal_day = 1.5;
    assert_eq!(time::JulDay(&date), 0.0);

}

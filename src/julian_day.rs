use time;
use util;

/*

    Returns the Julian day.
    -----------------------------------------------------------------
        date: A given date (with a decimal day)

*/

pub fn julian_day(mut date: time::date) -> f64 {

    if date.m == 1 || date.m == 2 {
        date.y = date.y - 1;
        date.m = date.m + 12;
    }

    let a = util::int((date.y as f64) / 100.0) as f64;
    let mut b;
    match date.t {
        time::calendar_type::gregorian => b = 2.0 - a + (util::int(a / 4.0) as f64),
        time::calendar_type::julian => b = 0.0,
    };

    (util::int(365.25 * ((date.y as f64) + 4716.0)) as f64) +
    (util::int(30.6001 * ((date.m as f64) + 1.0)) as f64) +
    (date.d as f64) +
    (b as f64) -
    1524.5

}

/*

    Returns a date (year, month, decimal day).
    -----------------------------------------------------------------
        jd: Julian Day (can't be negative)

*/

fn date_from_julian_day(mut jd: f64) -> (i16, i8, f64) {
    if jd < 0.0 {
        // panic
    }
    jd += 0.5;
    let Z = jd as i64;
    let F = jd - (Z as f64);
    let mut A;
    if Z < 2299161 {
        A = Z;
    }
    else {
        let alpha = util::int(((Z as f64) - 1867216.25)/36524.25);
        A = Z + 1 + alpha - util::int((alpha as f64)/4.0);
    }
    let B = A + 1524;
    let C = util::int(((B as f64) - 122.1)/365.25);
    let D = util::int(365.25 * (C as f64));
    let E = util::int(((B - D) as f64)/30.6001);
    let day = ((B - D) as f64) - (util::int(30.6001 * (E as f64)) as f64) + F;
    let month = if E < 14 {
                    E - 1
                }
                else if E == 14 || E == 15 { E - 13 }
                else {// panic
                0
                };
    let year = if month > 2 { C - 4716 }
               else if month == 1 || month == 2 { C - 4715 }
               else {// panic
               0
               };

    (year as i16, month as i8, day)
}

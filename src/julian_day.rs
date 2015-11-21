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

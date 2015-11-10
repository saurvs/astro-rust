use time;
use util;

/*

    Returns the Julian day.
    -----------------------------------------------------------------
        date: Gregorian date (with a decimal day)

*/

pub fn julian_day(mut date: time::greg_date) -> f64 {

    if date.m == 1 || date.m == 2 {
        date.y = date.y - 1;
        date.m = date.m + 12;
    }
    let a: i64 = util::int((date.y as f64) / 100.0);
    let b: i64 = 2 - a + util::int(a as f64 / 4.0);

    (util::int(365.25 * ((date.y as f64) + 4716.0)) as f64) +
    (util::int(30.6001 * ((date.m as f64) + 1.0)) as f64) +
    (date.d as f64) +
    (b as f64) -
    1524.5

}

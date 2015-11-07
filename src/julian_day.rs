use time;

// returns the julian day for a given Gregorian date (with a decimal day)
pub fn julian_day(mut date: time::greg_date) -> f64 {

    if date.m == 1 || date.m == 2 {
        date.y = date.y - 1;
        date.m = date.m + 12;
    }
    let a: i64 = int((date.y as f64) / 100.0);
    let b: i64 = 2 - a + int(a as f64 / 4.0);

    (int(365.25 * ((date.y as f64) + 4716.0)) as f64) +
    (int(30.6001 * ((date.m as f64) + 1.0)) as f64) +
    (date.d as f64) +
    (b as f64) -
    1524.5

}

// returns the largest integer less than or equal to (x)
fn int(x: f64) -> i64 {
    if x < 0_f64 {
        return -1 + x as i64;
    }
    x as i64
}

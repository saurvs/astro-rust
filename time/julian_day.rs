fn main() {
    let y = 1969;
    let m = 7;
    let d = decimal_day(20, 20, 18, 4.0);
    println!("{}", julian_day(y, m, d)); // day when the Apollo 11 Lunar Lander landed on the moon
}

// returns the julian day for a given (year, month, decimal day)
// currently assumes that given date is from the Gregorian calendar
fn julian_day(mut y: i32, mut m: i16, d: f64) -> f64 {

    if m == 1 || m == 2 {
        y = y - 1;
        m = m + 12;
    }
    let a: i64 = int((y as f64) / 100.0);
    let b: i64 = 2 - a + int(a as f64 / 4.0);

    (int(365.25 * ((y as f64) + 4716.0)) as f64) +
    (int(30.6001 * ((m as f64) + 1.0)) as f64) +
    (d as f64) +
    (b as f64) -
    1524.5

}

// returns the decimal day for a given (day, hour, minute, second)
// currently assumes hour, minute and second are in UTC
fn decimal_day(d: i32, h: i32, m: i32, s: f64) -> f64 {
    (d as f64) +
    (h as f64) / 24.0 +
    (m as f64) / 60.0 +
    s / 60.0
}

// returns the largest integer less than or equal to (x)
fn int(x: f64) -> i64 {
    if x < 0_f64 {
        return -1 + x as i64;
    }
    x as i64
}

pub mod julian_day;

// returns the decimal day for a given (day, hour, minute, second)
// currently assumes hour, minute and second are in UTC
pub fn decimal_day(d: i32, h: i32, m: i32, s: f64) -> f64 {
    (d as f64) +
    (h as f64) / 24.0 +
    (m as f64) / 60.0 +
    s / 60.0
}

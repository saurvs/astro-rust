// struct for representing a Gregorian date
// note that d here means decimal day
pub struct greg_date {
    pub y: i32,
    pub m: i16,
    pub d: f64,
}

// struct for representing a day in the usual sense
// using hours, minutes and seconds
pub struct usual_day {
    pub d: i16, // day of month
    pub h: i16, // hours of day
    pub m: i16, // minute of hour
    pub s: f64 // second of minute
}

// returns the decimal day for a given day expressed in the usual sense
// currently assumes hour, minute and second are in UTC
pub fn decimal_day(day: usual_day) -> f64 {
    (day.d as f64) +
    (day.h as f64) / 24.0 +
    (day.m as f64) / 60.0 +
    day.s / 60.0
}

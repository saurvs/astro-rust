/*

    An enum for representing different calendar types.
    -----------------------------------------------------------------
        gregorian: The Gregorian calendar
           julian: The Julian calendar

*/

pub enum calendar_type {
    gregorian,
    julian
}

/*

    A struct for representing a date.
    -----------------------------------------------------------------
        y: The year
        m: The month (1 - 12)
        d: The decimal day
        t: The calenday type

*/

pub struct date {
    pub y: i32,
    pub m: i8,
    pub d: f64,
    pub t: calendar_type,
}

/*

    A struct for representing a day in the usual sense.
    -----------------------------------------------------------------
        d: The day of the month (1 - 31)
        h: The hour of the day (0 - 24)
        m: The minute of the hour (0 - 60)
        s: The second of the minute (0.0 - 60.0)

*/

pub struct usual_day {
    pub d: i16,
    pub h: i16,
    pub m: i16,
    pub s: f64
}

/*

    Returns the decimal day for a given day expressed in the
    usual sense.
    -----------------------------------------------------------------
        usual_day: The an instance of the struct usual_day

*/

pub fn decimal_day(day: usual_day) -> f64 {
    (day.d as f64) +
    (day.h as f64) / 24.0 +
    (day.m as f64) / 60.0 +
    day.s / 60.0
}

/*

    Returns the time measured in Julian centuries from the
    Epoch J2000.
    -----------------------------------------------------------------
        jed: The Julian Emphemeris Day

*/

pub fn julian_centuries(jed: f64) -> f64 {
    (jed - 2451545.0) / 36525.0
}

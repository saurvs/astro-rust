/// Represents different calendar types
pub enum calendar_type {
    /// Gregorian calendar
    gregorian,
    /// Julian calendar
    julian
}

/// Represents a date
pub struct date {
    /// Year
    pub y: i32,
    /// Month
    ///
    /// range: *1 - 12*
    pub m: u8,
    /// Decimal day
    pub d: f64,
    /// Calenday type
    pub t: calendar_type,
}

/// Represents a day of month with hours, minutes and seconds
pub struct usual_day {
    /// Day of month
    ///
    /// range: *1 - 31*
    pub d: i16,
    /// Hour of day
    ///
    /// range: *0 - 60*
    pub h: u8,
    /// Minute of hour
    ///
    /// range: *0 - 60*
    pub m: u8,
    /// Second of minute
    ///
    /// range: *0.0 - 60.0*
    pub s: f64
}

/**
Returns the decimal day for a ```usual_day```

* ```usual_day```: A ```usual_day``` struct
**/
pub fn decimal_day(day: usual_day) -> f64 {
    (day.d as f64) +
    (day.h as f64) / 24.0 +
    (day.m as f64) / 60.0 +
    day.s / 60.0
}

/**
Returns time between the epoch J2000.0 and a given Julian Emphemeris Day, measured in Julian centuries

* ```jed```: Julian Emphemeris Day
**/
pub fn julian_century(jed: f64) -> f64 {
    (jed - 2451545.0) / 36525.0
}

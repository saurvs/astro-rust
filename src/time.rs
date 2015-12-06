use util;

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

/**
Returns a Julian day

# Arguments

```date```: A ```date``` struct
**/
pub fn julian_day(mut date: date) -> f64 {

    if date.m == 1 || date.m == 2 {
        date.y = date.y - 1;
        date.m = date.m + 12;
    }

    let a = util::int((date.y as f64) / 100.0) as f64;
    let mut b;
    match date.t {
        calendar_type::gregorian => b = 2.0 - a + (util::int(a / 4.0) as f64),
        calendar_type::julian => b = 0.0,
    };

    (util::int(365.25 * ((date.y as f64) + 4716.0)) as f64) +
    (util::int(30.6001 * ((date.m as f64) + 1.0)) as f64) +
    (date.d as f64) +
    (b as f64) -
    1524.5

}

/**
Returns a date equivalent to a given Julian day

# Return variables

```date_from_julian_day() -> (year, month, decimal_day)```

# Arguments

```jd```: Julian Day. **Can't be a negative value.**
**/
pub fn date_from_julian_day(mut jd: f64) -> (i16, i8, f64) {
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

/**
Returns mean sidereal time at any instant of Universal Time

Mean sidereal time is at the Greenwhich meridian.

# Arguments

* ```date```: A ```date``` struct
**/
pub fn mean_sidereal(date: date) -> f64 {
    let jd = julian_day(date);
    let t = julian_century(jd);

    280.46061837 +
    360.98564736629 * (jd - 2451545.0) +
    (t * t) * (0.000387933  - t / 38710000.0)
}

/**
Returns mean sidereal time at 0th hour of Universal Time

Mean sidereal time is at the Greenwhich meridian.

# Arguments

* ```date```: A ```date``` struct
**/
pub fn mean_sidereal_ut(date: date) -> f64 {
    let t = julian_century(julian_day(date));

    (100.46061837 +
    t * (36000.770053608 +
    t * (0.000387933 -
    t / 38710000.0))).to_radians()
}

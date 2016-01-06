use angle;
use util;
use nutation;
use planet;

/// Represents different **calendar types**
pub enum CalendarType {
    /// Gregorian calendar
    Gregorian,
    /// Julian calendar
    Julian
}

/// Represents a **date** with **year, month, decimal day** and **calendar type**
pub struct Date {
    /// Year
    pub year: i32,
    /// Month
    ///
    /// range: *1 - 12*
    pub month: u8,
    /// Decimal day
    pub decimal_day: f64,
    /// Calenday type
    pub calendar_type: CalendarType,
}

/// Represents a **day of month** with **hours, minutes and seconds**
pub struct DayOfMonth {
    /// Day of month
    ///
    /// range: *1 - 31*
    pub day: i16,
    /// Hour of day
    ///
    /// range: *0 - 60*
    pub hour: u8,
    /// Minute of hour
    ///
    /// range: *0 - 60*
    pub minute: u8,
    /// Second of minute
    ///
    /// range: *0.0 - 60.0*
    pub second: f64
}

/**
Returns the **decimal day** for a ```DayOfMonth```

* ```day_of_month```: A ```day_of_month``` struct
**/
pub fn DecimalDay(day: DayOfMonth) -> f64 {
      (day.day as f64)
    + (day.hour as f64) / 24.0
    + (day.minute as f64) / 60.0
    +  day.second / 60.0
}

/**
Returns the **decimal year** for a ```Date```

* ```date```: A ```date``` struct
**/
pub fn DecimalYear(date: Date) -> f64 {
    let mut y = 0;
    let mut days = 365.0;

    if date.month > 1 { y += 31; }
    if date.month > 2 {
        y += 28;
        if IsLeapYear(date.year, date.calendar_type) {
            y += 1;
            days += 1.0;
        }
    }
    if date.month > 3 { y += 31; }
    if date.month > 4 { y += 30; }
    if date.month > 5 { y += 31; }
    if date.month > 6 { y += 30; }
    if date.month > 7 { y += 31; }
    if date.month > 8 { y += 31; }
    if date.month > 9 { y += 30; }
    if date.month > 10 { y += 31; }
    if date.month > 11 { y += 30; }

    (date.year as f64) + ((y as f64) + date.decimal_day)/days
}

/**
Checks if a **year** is a **leap year**

# Arguments

* ```year```: Year
* ```calendar_type```: ```CalendarType``` enum
**/
pub fn IsLeapYear(year: i32, calendar_type: CalendarType) -> (bool) {
    match calendar_type {
        CalendarType::Julian => year % 4 == 0,
        CalendarType::Gregorian => {
            if year%100 == 0 { year%400 == 0 }
            else { year%4 == 0 }
        },
    };
    false
}

/**
Returns **Julian century**

# Arguments

* ```JED```: Julian Emphemeris day
**/
pub fn JulianCentury(JED: f64) -> f64 {
    (JED - 2451545.0) / 36525.0
}

/**
Returns **Julian millenium**

# Arguments

* ```JED```: Julian Emphemeris day
**/
pub fn JulianMillenium(JED: f64) -> f64 {
    (JED - 2451545.0) / 365250.0
}

/**
Returns **Julian day**

# Arguments

```date```: A ```date``` struct
**/
pub fn JulianDay(mut date: Date) -> f64 {

    if date.month == 1 || date.month == 2 {
        date.year = date.year - 1;
        date.month = date.month + 12;
    }

    let a = util::int((date.year as f64) / 100.0) as f64;
    let mut b;
    match date.calendar_type {
        CalendarType::Gregorian => b = 2.0 - a + (util::int(a/4.0) as f64),
        CalendarType::Julian => b = 0.0,
    };

      (util::int(365.25 * ((date.year as f64) + 4716.0)) as f64)
    + (util::int(30.6001 * ((date.month as f64) + 1.0)) as f64)
    + (date.decimal_day as f64)
    + (b as f64)
    - 1524.5

}

/**
Returns **Julian Emphemeris day**

# Arguments

```date```: A ```date``` struct
**/
pub fn JulianEmphemerisDay(mut date: Date) -> f64 {
    ApproximateDeltaT(date.year, date.month)/86400.0 + JulianDay(date)
}

#[macro_export]
macro_rules! JulianDay {
    ($a: expr, $b: expr, $c: expr, $d: expr, $e: expr, $f: expr, $g: expr) => {{
        let day = DayOfMonth{};
        let date = Date{};
        astro::time::JulianDay()
    }};
}

/**
Returns a ```Date``` **equivalent** to a given **Julian day**

# Returns

```(year, month, decimal_day)```

* ```year```: Year
* ```month```: Month
* ```decimal_day```: Decimal day

# Arguments

```jd```: Julian Day. **Can't be a negative value.**
**/
pub fn DateFromJulianDay(mut jd: f64) -> (i16, i8, f64) {
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
Returns **apparent sidereal time** at any instant of Universal Time

# Arguments

* ```JD```: Julian day
**/
pub fn ApparentSiderealTime(JD: f64) -> (i8, i8, f64) {
    let (hour, minute, seconds) = MeanSiderealTime(JD);

    let (nut_in_long, nut_in_oblq) = nutation::Corrections(JD);
    let eclip_oblq = planet::earth::ecliptic::MeanObliquity(JD);

    let seconds_correction =   nut_in_long.to_degrees()*3600.0
                             * (eclip_oblq + nut_in_oblq).cos()
                             / 15.0;

    (hour, minute, seconds + seconds_correction)
}

/**
Returns **mean sidereal time** at any instant of Universal Time

Mean sidereal time is at the Greenwhich meridian.

# Returns

```(hour, minute, seconds)```

* ```hour```: Hour *(range: 0 - 24)*
* ```minute```: Minute *(range: 0 - 60)*
* ```seconds```: Seconds *(range: 0.0 - 60.0)*

# Arguments

* ```JD```: Julian day
**/
pub fn MeanSiderealTime(JD: f64) -> (i8, i8, f64) {
    let JC = JulianCentury(JD);
    let angle = angle::LimitedTo360(  280.46061837
                                    + 360.98564736629 * (JD - 2451545.0)
                                    + JC*JC * (0.000387933 - JC/38710000.0)
                                   );

    HoursMinutesSecondsFromDegrees(angle)
}

pub fn HoursMinutesSecondsFromDegrees(angle: f64) -> (i8, i8, f64) {
    let hours = angle / 15.0;
    let hour = hours as i8;

    let minutes = (hours - (hour as f64)) * 60.0;
    let minute = minutes as i8;

    let seconds = (minutes - (minute as f64)) * 60.0;

    (hour, minute, seconds)
}

pub fn DegreesFromHoursMinutesSeconds(hour: i8, minute: i8, seconds: f64) -> f64 {
    (hour as f64)*15.0 + (minute as f64)/60.0 + seconds/3600.0
}

/**
Returns an approximate value of **ΔT** for a given year and month

This function approximates **ΔT** from polynomial expressions using a
method different from that given in the *Meeus* book. The method
used is given [here](http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html);
it covers a far wider time range.

# Arguments

* ```year```: Year
* ```month```: Month *(range: 1 - 12)*
**/
pub fn ApproximateDeltaT(year: i32, month: u8) -> f64 {
    let y = (year as f64) + ((month as f64) - 0.5)/12.0;

    if y < -500.0 {
        let u = (y-1820.0) / 100.0;
        return -20.0 + 32.0*u*u;
    }
    else if y < 500.0 {
        let u = y / 100.0;
        return 10583.6 -
               u * (1014.41 -
               u * (33.78311 -
               u * (5.952053 +
		       u * (0.1798452 -
               u * (0.022174192 -
                    u * 0.0090316521)
                   ))));
    }
    else if y < 1600.0 {
        let u = (y-1000.0) / 100.0;
        return 1574.2 -
               u * (556.01 -
               u * (71.23472 +
               u * (0.319781 -
		       u * (0.8503463 +
               u * (0.005050998 +
                    u * 0.0083572073)
                   ))));
    }
    else if y < 1700.0 {
        let u = y - 1600.0;
        return 120.0 -
               u * (0.9808 +
               u * (0.01532 -
                    u / 7129.0
                   ));
    }
    else if y < 1800.0 {
        let u = y - 1700.0;
        return 8.83 +
               u * (0.1603 -
               u * (0.0059285 -
               u * (0.00013336 -
                    u / 1174000.0
                   )));
    }
    else if y < 1860.0 {
        let u = y - 1800.0;
        return 13.72 -
               u * (0.332447 -
               u * (0.0068612 +
               u * (0.0041116 -
               u * (0.00037436 -
               u * (0.0000121272 -
               u * (0.0000001699 -
                    u * 0.000000000875
                   ))))));
    }
    else if y < 1900.0 {
        let u = y - 1860.0;
        return 7.62 +
               u * (0.5737 -
               u * (0.251754 -
               u * (0.01680668 -
               u * (0.0004473624 -
                    u / 233174.0
                   ))));
    }
    else if y < 1920.0 {
        let u = y - 1900.0;
        return -2.79 +
               u * (1.494119 -
               u * (0.0598939 -
               u * (0.0061966 +
                    u * 0.000197
                   )));
    }
    else if y < 1941.0 {
        let u = y - 1920.0;
        return 21.20 +
               u * (0.84493 -
               u * (0.076100 -
                    u * 0.0020936
                   ));
    }
    else if y < 1961.0 {
        let u = y - 1950.0;
        return 29.07 +
               u * (0.407 -
               u * ((1.0 / 233.0) -
                    u / 2547.0
                   ));
    }
    else if y < 1986.0 {
        let u = y - 1975.0;
        return 45.45 +
               u * (1.067 -
               u * ((1.0 / 260.0) +
                    u / 718.0)
                   );
    }
    else if y < 2005.0 {
        let u = y - 2000.0;
        return 63.86 +
               u * (0.3345 -
               u * (0.060374 -
               u * (0.0017275 +
               u * (0.000651814 +
                    u * 0.00002373599)
                   )));
    }
    else if y < 2050.0 {
        let u = y - 2000.0;
        return 62.92 +
               u * (0.32217 -
                    u * 0.005589);
    }
    else if y < 2150.0 {
        let u = (y-1820.0) / 100.0;
        return -20.0 + 32.0*u*u - 0.5628*(2150.0 - y);
    }
    else if y >= 2150.0 {
        let u = (y-1820.0) / 100.0;
        return -20.0 + 32.0*u*u;
    }

    0.0
}

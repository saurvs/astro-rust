use util;

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
    pub y: i32,
    /// Month
    ///
    /// range: *1 - 12*
    pub m: u8,
    /// Decimal day
    pub d: f64,
    /// Calenday type
    pub t: CalendarType,
}

/// Represents a **day of month** with **hours, minutes and seconds**
pub struct UsualDay {
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
Computes the **decimal day** for a ```UsualDay```

* ```usual_day```: A ```usual_day``` struct
**/
pub fn DecimalDay(day: UsualDay) -> f64 {
    (day.d as f64) +
    (day.h as f64) / 24.0 +
    (day.m as f64) / 60.0 +
    day.s / 60.0
}

/**
Computes the **decimal year** for a ```Date```

* ```date```: A ```date``` struct
**/
pub fn DecimalYear(date: Date) -> f64 {
    let mut y = 0;
    let mut days = 365.0;

    if date.m > 1 {
        y += 31;
    }
    if date.m > 2 {
        y += 28;
        if IsLeapYear(date.y, date.t) {
            y += 1;
            days += 1.0;
        }
    }
    if date.m > 3 {
        y += 31;
    }
    if date.m > 4 {
        y += 30;
    }
    if date.m > 5 {
        y += 31;
    }
    if date.m > 6 {
        y += 30;
    }
    if date.m > 7 {
        y += 31;
    }
    if date.m > 8 {
        y += 31;
    }
    if date.m > 9 {
        y += 30;
    }
    if date.m > 10 {
        y += 31;
    }
    if date.m > 11 {
        y += 30;
    }

    (date.y as f64) + ((y as f64) + date.d)/days
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
            if year%100 == 0 {
                year%400 == 0
            }
            else {
                year%4 == 0
            }
        },
    };
    false
}

/**
Computes a **Julian century**, the time between the epoch J2000.0 and a given Julian Emphemeris Day

* ```jed```: Julian Emphemeris day
**/
pub fn JulianCentury(jed: f64) -> f64 {
    (jed - 2451545.0) / 36525.0
}

/**
Computes a **Julian day**

# Arguments

```date```: A ```date``` struct
**/
pub fn JulianDay(mut date: Date) -> f64 {

    if date.m == 1 || date.m == 2 {
        date.y = date.y - 1;
        date.m = date.m + 12;
    }

    let a = util::int((date.y as f64) / 100.0) as f64;
    let mut b;
    match date.t {
        CalendarType::Gregorian => b = 2.0 - a + (util::int(a / 4.0) as f64),
        CalendarType::Julian => b = 0.0,
    };

    (util::int(365.25 * ((date.y as f64) + 4716.0)) as f64) +
    (util::int(30.6001 * ((date.m as f64) + 1.0)) as f64) +
    (date.d as f64) +
    (b as f64) -
    1524.5

}

/**
Returns a ```Date``` **equivalent** to a given **Julian day**

# Return variables

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
Computes **mean sidereal time** at any instant of Universal Time

Mean sidereal time is at the Greenwhich meridian.

# Arguments

* ```date```: A ```date``` struct
**/
pub fn MeanSiderealTime(date: Date) -> f64 {
    let jd = JulianDay(date);
    let t = JulianCentury(jd);

    280.46061837 +
    360.98564736629 * (jd - 2451545.0) +
    (t * t) * (0.000387933  - t / 38710000.0)
}

/**
Computes **mean sidereal time** at 0th hour of Universal Time

Mean sidereal time is at the Greenwhich meridian.

# Arguments

* ```date```: A ```date``` struct
**/
pub fn MeanSiderealTimeAt0thHour(date: Date) -> f64 {
    let t = JulianCentury(JulianDay(date));

    (100.46061837 +
    t * (36000.770053608 +
    t * (0.000387933 -
    t / 38710000.0))).to_radians()
}

/**
Computes an approximate value of **ΔT** for a given year and month

This function approximates **ΔT** from polynomial expressions using a
method different from that given in the *Meeus* book. The method
used is given [here](http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html);
it covers a far wider time range.

# Arguments

* ```year```: Year
* ```month```: Month *(range: 1 - 12)*
**/
pub fn ApproximateDeltaT(year: i32, month: i8) -> f64 {
    let y = (year as f64) + ((month as f64) - 0.5)/12.0;

    if y < -500.0 {
        let u = (y - 1820.0)/100.0;
        return -20.0 + 32.0*u*u;
    }
    else if y < 500.0 {
        let u = y/100.0;
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
        let u = (y - 1000.0)/100.0;
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
                    u * (1.0 / 7129.0)
                   ));
    }
    else if y < 1800.0 {
        let u = y - 1700.0;
        return 8.83 +
               u * (0.1603 -
               u * (0.0059285 -
               u * (0.00013336 -
                    u * (1.0 / 1174000.0)
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
                    u * (1.0 / 233174.0)
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
                    u * 0.0020936)
                   );
    }
    else if y < 1961.0 {
        let u = y - 1950.0;
        return 29.07 +
               u * (0.407 -
               u * ((1.0 / 233.0) -
                    u * (1.0 / 2547.0))
                   );
    }
    else if y < 1986.0 {
        let u = y - 1975.0;
        return 45.45 +
               u * (1.067 -
               u * ((1.0 / 260.0) +
                    u * (1.0 / 718.0))
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
        let u = (y - 1820.0) / 100.0;
        return -20.0 + 32.0*u*u - 0.5628*(2150.0 - y);
    }
    else if y >= 2150.0 {
        let u = (y - 1820.0) / 100.0;
        return -20.0 + 32.0*u*u;
    }

    0.0
}

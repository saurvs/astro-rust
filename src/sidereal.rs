use time;
use julian_day;
use nutation;

/*

    NOTE: All angles passed as arguments, and those returned,
          are assumed to be radians, even if the comments
          describe them with degrees.

*/

fn t(jd: f64) -> f64 {
    (jd - 2451545.0) / 36525.0
}

/*

    mean_sidereal_greenwhich_zero_ut(date)
        -> (mean_sidereal_at_greenwhich_at_zero_hour_universal_time)

    mean_sidereal_at_greenwhich_at_zero_hour_universal_time is
    returned in radians, which then needs to be converted to degrees,
    and then to hours, minutes and seconds.
    -----------------------------------------------------------------
    Returns the mean sidereal time at the meridian of Greenwhich at
    the 0th hour of Universal Time of the given date

*/

pub fn mean_sidereal_greenwhich_zero_ut(date: time::date) -> f64 {
    let jd = julian_day::julian_day(date);
    let t = t(jd);

    (100.46061837 +
    t * (36000.770053608 +
    t * (0.000387933 -
    t / 38710000.0))).to_radians()
}

/*

    mean_sidereal_greenwhich(date)
        -> (mean_sidereal_at_greenwhich_at_any_instant_universal_time)

    mean_sidereal_at_greenwhich_at_any_instant_universal_time is
    returned in radians, which then needs to be converted to degrees,
    and then to hours, minutes and seconds.
    -----------------------------------------------------------------
    Returns the mean sidereal time at the meridian of Greenwhich at
    any instant of Universal Time of a given date

*/

pub fn mean_sidereal_greenwhich(date: time::date) -> f64 {
    let jd = julian_day::julian_day(date);
    let t = t(jd);

    280.46061837 +
    360.98564736629 * (jd - 2451545.0) +
    (t * t) * (0.000387933  - t / 38710000.0)
}

/*

    app_sidereal_greenwhich(date)
        -> (app_sidereal_at_greenwhich_at_any_instant_universal_time)

    app_sidereal_at_greenwhich_at_any_instant_universal_time is
    returned in radians, which then needs to be converted to degrees,
    and then to hours, minutes and seconds.
    -----------------------------------------------------------------
    Returns the apparent sidereal time at the meridian of Greenwhich
    at any instant of Universal Time of a given date

*/

pub fn app_sidereal_greenwhich(date: time::date) -> f64 {
    mean_sidereal_greenwhich(date)
}

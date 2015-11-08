use time;
use julian_day;

fn t(jd: f64) -> f64 {
    (jd - 2451545.0) / 36525.0
}

// returns the sidereal time at the meridian of Greenwhich at the
// zero-th hour of Universal Time of a given Gregorian date
pub fn mean_sidereal_greenwhich_zero_ut(date: time::greg_date) -> f64 {
    let jd = julian_day::julian_day(date);
    let t = t(jd);

    (100.46061837 +
    t * (36000.770053608 + t * (0.000387933 - t / 38710000.0))).to_radians()
}

// returns the sidereal time at the meridian of Greenwhich at the
// any instant of Universal Time of a given Gregorian date
pub fn mean_sidereal_greenwhich(date: time::greg_date) -> f64 {
    let jd = julian_day::julian_day(date);
    let t = t(jd);

    280.46061837 +
    360.98564736629 * (jd - 2451545.0) +
    (t * t) * (0.000387933  - t / 38710000.0)
}

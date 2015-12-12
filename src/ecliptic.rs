use angle;
use time;

/**
Returns the mean obliquity of the Earth's ecliptic (in radians)

The obliquity of the ecliptic, or the inclination of the Earth's axis
of rotation, is the angle between the equator and the ecliptic.

The *mean* obliquity isn't corrected for nutation.
To obtain the *true* obliquity, use [```nutation()```]
(./fn.nutation.html) to get the nutation correction for obliquity,
and add it to the *mean* obliquity.

# Arguments

```julian_ephemeris_day```: Time in Julian Ephemeris days
**/
pub fn mean_obliquity(jed: f64) -> (f64) {
    let u = time::julian_century(jed) / 100.0;

    (angle::pure_degrees(23.0, 26.0, 21.448) -
    u * (angle::pure_degrees(0.0, 0.0, 4680.93) +
    u * (angle::pure_degrees(0.0, 0.0, 1.55) +
    u * (angle::pure_degrees(0.0, 0.0, 1999.25) -
    u * (angle::pure_degrees(0.0, 0.0, 51.38) -
    u * (angle::pure_degrees(0.0, 0.0, 249.67) +
    u * (angle::pure_degrees(0.0, 0.0, 39.05) +
    u * (angle::pure_degrees(0.0, 0.0, 7.12) -
    u * (angle::pure_degrees(0.0, 0.0, 27.87) +
    u * (angle::pure_degrees(0.0, 0.0, 5.79) +
    u * angle::pure_degrees(0.0, 0.0, 2.45)
    )))))))))).to_radians()
}

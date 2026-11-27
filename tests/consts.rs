//! Meeus, *Astronomical Algorithms*, 2nd ed., Appendix I: constants.

mod common;

use astro::consts;
use common::assert_close;

#[test]
fn distances_match_appendix_i() {
    // 1 AU = 149 597 870 km (IAU 1976) = 499.0048 light-seconds
    //      = 8.32 light-minutes = 0.0057755183 light-day.
    assert_close(consts::ASTRONOMICAL_UNIT / 1000.0, 149_597_870.0, 1.0);
    assert_close(consts::LIGHT_TIME_PER_AU * 86400.0, 499.0048, 0.00005);
    assert_close(consts::LIGHT_TIME_PER_AU * 1440.0, 8.32, 0.005);
    assert_close(consts::LIGHT_TIME_PER_AU, 0.0057755183, 5e-11);

    // 1 parsec = 30.8568e12 km = 3.2616 ly = 206264.8 AU.
    assert_close(consts::PARSEC / 1e15, 30.8568, 0.00005);
    assert_close(consts::PARSEC_IN_LIGHT_YEARS, 3.2616, 0.00005);
    assert_close(consts::PARSEC_IN_AU, 206264.8, 0.05);

    // 1 light-year = 9.4607e12 km = 0.30660 pc = 63241 AU.
    assert_close(consts::LIGHT_YEAR / 1e15, 9.4607, 0.00005);
    assert_close(1.0 / consts::PARSEC_IN_LIGHT_YEARS, 0.30660, 0.000005);
    assert_close(consts::LIGHT_YEAR_IN_AU, 63241.0, 0.5);

    assert_close(consts::MEAN_EARTH_MOON_DIST / 1000.0, 384400.0, 1e-9);
    assert_close(consts::SUN_DIAMETER / 1000.0, 1_392_000.0, 1e-9);
    assert_close(consts::MOON_DIAMETER / 1000.0, 3476.0, 1e-9);
}

#[test]
fn earth_figure_iau1976() {
    assert_close(consts::iau1976::EQUATORIAL_RADIUS / 1000.0, 6378.14, 1e-9);
    assert_close(consts::iau1976::POLAR_RADIUS / 1000.0, 6356.76, 0.005);
    assert_close(1.0 / consts::iau1976::FLATTENING, 298.257, 1e-9);
    // WGS 84 and IAU 1976 agree to a few metres.
    assert!((consts::wgs84::EQUATORIAL_RADIUS - consts::iau1976::EQUATORIAL_RADIUS).abs() < 5.0);
    assert!((consts::wgs84::POLAR_RADIUS - consts::iau1976::POLAR_RADIUS).abs() < 5.0);
}

#[test]
fn mathematical_constants() {
    assert_close(consts::ARCSEC_PER_RADIAN, 206264.806247, 5e-7);
    assert_close(
        1.0_f64.to_degrees() * 3600.0,
        consts::ARCSEC_PER_RADIAN,
        1e-6,
    );
}

#[test]
fn time_constants() {
    // 1 sidereal day = 23h 56m 04.0905s of mean solar time.
    assert_close(
        consts::SIDEREAL_DAY * 86400.0,
        23.0 * 3600.0 + 56.0 * 60.0 + 4.0905,
        0.00005,
    );
    assert_close(
        consts::SIDEREAL_DAY * consts::SOLAR_DAY_IN_SIDEREAL_DAYS,
        1.0,
        1e-10,
    );
    assert_close(consts::JULIAN_CENTURY, 100.0 * consts::JULIAN_YEAR, 0.0);
    assert_close(consts::J2000, 2451545.0, 0.0);

    // Ordering of the year and month lengths.
    const { assert!(consts::TROPICAL_YEAR < consts::SIDEREAL_YEAR) };
    const { assert!(consts::SIDEREAL_YEAR < consts::ANOMALISTIC_YEAR) };
    const { assert!(consts::DRACONIC_MONTH < consts::TROPICAL_MONTH) };
    const { assert!(consts::TROPICAL_MONTH < consts::SIDEREAL_MONTH) };
    const { assert!(consts::SIDEREAL_MONTH < consts::ANOMALISTIC_MONTH) };

    // Synodic month from the sidereal month and year: 1/S = 1/P - 1/Y.
    let syn = 1.0 / (1.0 / consts::SIDEREAL_MONTH - 1.0 / consts::SIDEREAL_YEAR);
    assert_close(syn, consts::SYNODIC_MONTH, 0.000005);
}

#[test]
fn varia() {
    assert_close(
        consts::SOLAR_PARALLAX * consts::ARCSEC_PER_RADIAN,
        8.79415,
        1e-12,
    );
    assert_close(
        consts::ABERRATION_CONSTANT * consts::ARCSEC_PER_RADIAN,
        20.4955,
        1e-12,
    );
    assert_close(consts::GAUSS_GRAV_DEG, 0.9856076686, 5e-11);
    assert_close(consts::SPEED_OF_LIGHT / 1000.0, 299792.458, 1e-9);
    assert_close(consts::EARTH_MOON_MASS_RATIO, 81.3007, 0.0);
    assert_close(consts::SUN_EARTH_MASS_RATIO, 332946.0, 0.0);
    assert_close(consts::MOON_EARTH_RADIUS_RATIO, 0.272481, 0.0);

    // The solar parallax is the angle subtended by the Earth's equatorial
    // radius at 1 AU.
    let p = (consts::iau1976::EQUATORIAL_RADIUS / consts::ASTRONOMICAL_UNIT).asin();
    assert_close(p * consts::ARCSEC_PER_RADIAN, 8.794, 0.0005);
}

//! Ephemeris for physical observations of the Moon (Meeus, chapter 53).

#![allow(non_snake_case)]

mod common;

use astro::lunar::{self, libration, position};
use astro::*;
use common::{assert_close, assert_close_deg};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

/// 1992 April 12, 0h TD
const JDE: f64 = 2448724.5;

fn rad(x: f64) -> f64 {
    x.to_radians()
}

fn deg(x: f64) -> f64 {
    x.to_degrees()
}

/// Signed difference of two angles in degrees, in (−180, 180].
fn diff_deg(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(360.0);
    if d > 180.0 { d - 360.0 } else { d }
}

// Inputs of example 53.a (from example 47.a)
const MOON_LONG_NO_NUT: f64 = 133.162655;
const MOON_LAT: f64 = -3.229126;
const NUT_LONG: f64 = 0.004610;
const TRUE_OBLQ: f64 = 23.440636;
const MOON_ASC: f64 = 134.688470;

// Meeus, example 53.a, with the old tuple-returning functions.
#[test]
fn example_53a_old_api() {
    let (l1, b1) = lunar::optical_libr(JDE, rad(MOON_LONG_NO_NUT), rad(MOON_LAT));
    assert_close(deg(l1), -1.206, 5e-4);
    assert_close(deg(b1), 4.194, 5e-4);

    let (l2, b2) = lunar::physical_libr(JDE, rad(MOON_LONG_NO_NUT), rad(MOON_LAT), b1);
    assert_close(deg(l2), -0.025, 5e-4);
    assert_close(deg(b2), 0.006, 5e-4);

    let (l, b) = lunar::total_libr(JDE, rad(MOON_LONG_NO_NUT), rad(MOON_LAT));
    assert_close(deg(l), -1.23, 5e-3);
    assert_close(deg(b), 4.20, 5e-3);
    assert_close(l, l1 + l2, 1e-15);
    assert_close(b, b1 + b2, 1e-15);

    let node = lunar::mn_ascend_node(time::julian_cent(JDE));
    assert_close_deg(deg(node), 274.400656, 1e-6);
    let P =
        lunar::pos_angl_of_axis_of_rot(JDE, node, b, rad(NUT_LONG), rad(TRUE_OBLQ), rad(MOON_ASC));
    assert_close(deg(P), 15.08, 5e-3);
}

// Meeus, example 53.a, with the struct-returning functions.
#[test]
fn example_53a_librations() {
    let l = libration::librations(JDE, rad(MOON_LONG_NO_NUT), rad(MOON_LAT));
    assert_close(deg(l.optical_long), -1.206, 5e-4);
    assert_close(deg(l.optical_lat), 4.194, 5e-4);
    assert_close(deg(l.physical_long), -0.025, 5e-4);
    assert_close(deg(l.physical_lat), 0.006, 5e-4);
    assert_close(deg(l.total_long), -1.23, 5e-3);
    assert_close(deg(l.total_lat), 4.20, 5e-3);

    let P = libration::pos_angl_of_axis(
        JDE,
        rad(MOON_LONG_NO_NUT),
        rad(MOON_LAT),
        rad(MOON_ASC),
        rad(NUT_LONG),
        rad(TRUE_OBLQ),
    );
    assert_close(deg(P), 15.08, 5e-3);
}

// Meeus, example 53.b, with the book's inputs: λ0 = 22.33978°,
// Δ = 368406 km, R = 149971500 km.
#[test]
fn example_53b_selenographic_sun() {
    let s = libration::selenographic_sun(
        JDE,
        rad(22.33978),
        rad(MOON_LONG_NO_NUT + NUT_LONG),
        rad(MOON_LAT),
        368406.0,
        149971500.0,
        rad(NUT_LONG),
    );
    assert_close(deg(s.long), 67.89, 5e-3);
    assert_close(deg(s.lat), 1.46, 5e-3);
    assert_close(deg(s.colong), 22.11, 5e-3);
    assert_close_deg(deg(s.colong), 90.0 - deg(s.long), 1e-12);

    // R given in AU and Δ in AU: only the ratio matters
    let t = libration::selenographic_sun(
        JDE,
        rad(22.33978),
        rad(MOON_LONG_NO_NUT + NUT_LONG),
        rad(MOON_LAT),
        368406.0 / 149597870.7,
        1.00249769,
        rad(NUT_LONG),
    );
    assert_close(deg(t.long), deg(s.long), 1e-6);
    assert_close(deg(t.lat), deg(s.lat), 1e-6);
}

// Meeus, examples 53.a and 53.b, everything computed from the instant.
#[test]
fn examples_53a_53b_from_jde() {
    let e = lunar::physical_ephemeris(JDE);
    assert_close(deg(e.libr.optical_long), -1.206, 5e-4);
    assert_close(deg(e.libr.optical_lat), 4.194, 5e-4);
    assert_close(deg(e.libr.physical_long), -0.025, 5e-4);
    assert_close(deg(e.libr.physical_lat), 0.006, 5e-4);
    assert_close(deg(e.libr.total_long), -1.23, 5e-3);
    assert_close(deg(e.libr.total_lat), 4.20, 5e-3);
    assert_close(deg(e.pos_angl_axis), 15.08, 5e-3);
    assert_close(deg(e.sun.long), 67.89, 5e-3);
    assert_close(deg(e.sun.lat), 1.46, 5e-3);
    assert_close(deg(e.sun.colong), 22.11, 5e-3);

    // same as feeding the computed positions to the lower-level functions
    let moon = lunar::apparent_position(JDE);
    let l = libration::librations(JDE, moon.long - moon.nut_long, moon.lat);
    assert_eq!(l, e.libr);
    let P = libration::pos_angl_of_axis(
        JDE,
        moon.long - moon.nut_long,
        moon.lat,
        moon.asc,
        moon.nut_long,
        moon.true_oblq,
    );
    assert_close(P, e.pos_angl_axis, 1e-15);
}

// Meeus, example 53.c: sunrise at the crater Copernicus (η = −20.0°,
// θ = +9.7°) in April 1992.
#[test]
fn example_53c_sunrise_copernicus() {
    let (eta, theta) = (rad(-20.0), rad(9.7));

    // with the rounded b0 = 1.46°, c0 = 22.11° of example 53.b
    let sun = libration::SelenographicSun {
        long: rad(67.89),
        lat: rad(1.46),
        colong: rad(22.11),
    };
    let h = libration::sun_alt_at_lunar_point(&sun, eta, theta);
    assert_close(deg(h), 2.3253, 1e-4);
    // correction of eq. 53.4 for sunrise
    assert_close(-deg(h) / (12.19075 * theta.cos()), -0.1935, 1e-4);

    // with the unrounded b0, c0
    let e = lunar::physical_ephemeris(JDE);
    assert_close(
        deg(libration::sun_alt_at_lunar_point(&e.sun, eta, theta)),
        2.3253,
        0.005,
    );

    // 1992 April 11.8065 TD
    let t = libration::lunar_sunrise(JDE, eta, theta).unwrap();
    assert_close(t, 2448724.5 - 0.1935, 1e-4);
    let at = lunar::physical_ephemeris(t);
    assert_close(deg(at.sun.lat), 1.46, 5e-3);
    assert_close(deg(at.sun.colong), 19.75, 5e-3);
    assert_close(
        libration::sun_alt_at_lunar_point(&at.sun, eta, theta),
        0.0,
        1e-6,
    );
}

fn sun_alt(jde: f64, eta: f64, theta: f64) -> f64 {
    libration::sun_alt_at_lunar_point(&lunar::physical_ephemeris(jde).sun, eta, theta)
}

// Sunrise and sunset for points all over the visible hemisphere: the Sun is
// on the horizon, rising or setting as appropriate, and the lunar day lasts
// about half a synodic month.
#[test]
fn sunrise_and_sunset() {
    for &(eta, theta) in &[
        (-20.0, 9.7),
        (45.0, -30.0),
        (-70.0, 60.0),
        (80.0, -75.0),
        (0.0, 0.0),
        (179.0, 10.0),
    ] {
        let (eta, theta) = (rad(eta), rad(theta));
        for k in 0..6 {
            let start = JDE + k as f64 * 5.3;
            let rise = libration::lunar_sunrise(start, eta, theta).unwrap();
            assert!((rise - start).abs() < 16.0);
            assert_close(sun_alt(rise, eta, theta), 0.0, 1e-6);
            assert!(sun_alt(rise - 0.05, eta, theta) < 0.0);
            assert!(sun_alt(rise + 0.05, eta, theta) > 0.0);

            let set = libration::lunar_sunset(rise + 14.77, eta, theta).unwrap();
            assert_close(sun_alt(set, eta, theta), 0.0, 1e-6);
            assert!(sun_alt(set - 0.05, eta, theta) > 0.0);
            assert!(sun_alt(set + 0.05, eta, theta) < 0.0);
            // half a synodic month, give or take the effect of b0 and the
            // variable speed of the Moon at high latitudes
            let day = set - rise;
            assert!(day > 13.0 && day < 16.5, "day length {day}");
        }
    }
}

#[test]
fn sunrise_without_answer() {
    // at the pole the method does not apply
    assert_eq!(libration::lunar_sunrise(JDE, 0.0, FRAC_PI_2), None);
    assert_eq!(libration::lunar_sunset(JDE, 0.0, -FRAC_PI_2), None);
    // near the north pole with b0 = +1.46° the Sun stays up
    assert_eq!(libration::lunar_sunset(JDE, 0.0, rad(89.9)), None);
    assert_eq!(libration::lunar_sunrise(f64::NAN, 0.0, 0.0), None);
    assert_eq!(libration::lunar_sunrise(JDE, f64::INFINITY, 0.0), None);
}

#[test]
fn sun_alt_limits() {
    let sun = libration::SelenographicSun {
        long: rad(30.0),
        lat: rad(1.2),
        colong: rad(60.0),
    };
    // the subsolar point has the Sun in the zenith
    let h = libration::sun_alt_at_lunar_point(&sun, sun.long, sun.lat);
    assert_close(deg(h), 90.0, 1e-6);
    // its antipode has it at the nadir
    let h = libration::sun_alt_at_lunar_point(&sun, sun.long + PI, -sun.lat);
    assert_close(deg(h), -90.0, 1e-6);
    // longitudes are taken modulo 360°
    let a = libration::sun_alt_at_lunar_point(&sun, rad(-20.0), rad(9.7));
    let b = libration::sun_alt_at_lunar_point(&sun, rad(340.0), rad(9.7));
    assert_close(a, b, 1e-14);
}

// Ranges and rates over 30 years, plus a few remote dates.
#[test]
fn ranges_and_rates() {
    let mut max_l = 0.0_f64;
    let mut max_b = 0.0_f64;
    let mut max_b0 = 0.0_f64;
    for k in 0..3000 {
        let jde = 2444000.0 + k as f64 * 3.65;
        let e = lunar::physical_ephemeris(jde);
        let l = e.libr;
        assert!(l.total_long > -PI && l.total_long <= PI);
        // Meeus says the physical libration never exceeds 0.04°, but his own
        // series give up to 0.044° in longitude and 0.061° in latitude: the
        // −0.02245 sin F term of ρ and the +0.02244 cos F term of σ add an
        // almost constant +0.0224° to b″.
        assert!(
            deg(l.physical_long).abs() < 0.045,
            "{}",
            deg(l.physical_long)
        );
        assert!(deg(l.physical_lat).abs() < 0.065, "{}", deg(l.physical_lat));
        assert_close(l.total_long, l.optical_long + l.physical_long, 1e-15);
        assert_close(l.total_lat, l.optical_lat + l.physical_lat, 1e-15);
        max_l = max_l.max(deg(l.total_long).abs());
        max_b = max_b.max(deg(l.total_lat).abs());

        let P = deg(e.pos_angl_axis);
        assert!(P <= 25.0 || (335.0..360.0).contains(&P), "P = {P}");

        assert!((0.0..TAU).contains(&e.sun.colong));
        assert!(e.sun.long > -PI && e.sun.long <= PI);
        max_b0 = max_b0.max(deg(e.sun.lat).abs());

        // colongitude ≈ D − 90°, within the optical libration
        let D = deg(position::fundamental_args(jde).mn_elong);
        assert!(diff_deg(deg(e.sun.colong), D - 90.0).abs() < 10.0);

        // and it increases by about 12.19° per day
        let e1 = lunar::physical_ephemeris(jde + 1.0);
        let rate = diff_deg(deg(e1.sun.colong), deg(e.sun.colong));
        assert!((rate - 12.19).abs() < 0.6, "rate {rate}");
    }
    assert!(max_l > 7.0 && max_l < 8.5, "{max_l}");
    assert!(max_b > 6.5 && max_b < 7.0, "{max_b}");
    // |b0| is at most I = 1.54° plus the physical libration
    assert!(max_b0 > 1.45 && max_b0 < 1.61, "{max_b0}");

    // dates before 1582 and negative years still give sane values
    for &jde in &[0.0, 1000000.5, 2000000.5, 2299160.5, 3000000.5] {
        let e = lunar::physical_ephemeris(jde);
        assert!(deg(e.libr.total_long).abs() < 9.0);
        assert!(deg(e.libr.total_lat).abs() < 7.5);
        assert!(e.pos_angl_axis.is_finite());
    }
}

// The libration functions take the longitude modulo 2π.
#[test]
fn longitude_wrap_around() {
    let a = libration::librations(JDE, rad(MOON_LONG_NO_NUT), rad(MOON_LAT));
    let b = libration::librations(JDE, rad(MOON_LONG_NO_NUT - 720.0), rad(MOON_LAT));
    assert_close(a.total_long, b.total_long, 1e-12);
    assert_close(a.total_lat, b.total_lat, 1e-12);
}

const DELTA_T: f64 = 58.7;

fn observers() -> Vec<coords::GeographPoint> {
    // longitude positive west
    [
        (45.0, 0.0),
        (-33.87, -151.21), // Sydney: southern, eastern
        (51.48, 0.0),
        (19.82, 155.47), // Mauna Kea: western
        (-70.0, 30.0),
        (0.0, -90.0),
        (64.0, -120.0),
    ]
    .iter()
    .map(|&(lat, long)| coords::GeographPoint {
        long: rad(long),
        lat: rad(lat),
    })
    .collect()
}

// The direct and the differential methods agree, and the topocentric
// values differ from the geocentric ones by at most about 1°.
#[test]
fn topocentric_methods_agree() {
    for obs in observers() {
        for k in 0..120 {
            let jde = JDE + k as f64 * 0.29;
            let a = libration::topocent_librations(jde, DELTA_T, &obs, 0.0);
            let b = libration::topocent_librations_by_diff_corrections(jde, DELTA_T, &obs);
            assert_close(deg(a.long), deg(b.long), 0.006);
            assert_close(deg(a.lat), deg(b.lat), 0.006);
            assert_close_deg(deg(a.pos_angl_axis), deg(b.pos_angl_axis), 0.012);

            let g = lunar::physical_ephemeris(jde);
            assert!(diff_deg(deg(a.long), deg(g.libr.total_long)).abs() < 1.05);
            assert!((deg(a.lat) - deg(g.libr.total_lat)).abs() < 1.05);
            assert!(a.long > -PI && a.long <= PI);
            assert!((0.0..TAU).contains(&a.pos_angl_axis));
        }
    }
}

// The differential-correction convenience function agrees with the old
// low-level function fed with the same quantities.
#[test]
fn topocentric_diff_corrections_old_api() {
    let obs = coords::GeographPoint {
        long: rad(-151.21),
        lat: rad(-33.87),
    };
    let moon = lunar::apparent_position(JDE);
    let g = lunar::physical_ephemeris(JDE);
    let sidr = time::apprnt_sidr(
        time::mn_sidr(JDE - DELTA_T / 86400.0),
        moon.nut_long,
        moon.true_oblq,
    );
    let H = coords::hr_angl_frm_observer_long(sidr, obs.long, moon.asc);
    let (dl, db, dP) = lunar::topocent_libr_by_diff_corrections(
        obs.lat,
        moon.dec,
        H,
        moon.eq_hz_parllx,
        g.pos_angl_axis,
        g.libr.total_lat,
    );
    let t = libration::topocent_librations_by_diff_corrections(JDE, DELTA_T, &obs);
    assert_close(t.long, g.libr.total_long + dl, 1e-12);
    assert_close(t.lat, g.libr.total_lat + db, 1e-12);
    assert_close_deg(deg(t.pos_angl_axis), deg(g.pos_angl_axis + dP), 1e-9);
}

// An observer at the sub-lunar point sees the Moon in the zenith: no
// parallax, so the topocentric librations equal the geocentric ones (up to
// the difference between geographic and geocentric latitude). This also
// checks the west-positive longitude convention.
#[test]
fn topocentric_at_sublunar_point() {
    for k in 0..20 {
        let jde = JDE + k as f64 * 1.9;
        let moon = lunar::apparent_position(jde);
        let sidr = time::apprnt_sidr(
            time::mn_sidr(jde - DELTA_T / 86400.0),
            moon.nut_long,
            moon.true_oblq,
        );
        // H = θ0 − L − α = 0
        let obs = coords::GeographPoint {
            long: sidr - moon.asc,
            lat: moon.dec,
        };
        let g = lunar::physical_ephemeris(jde);
        let a = libration::topocent_librations(jde, DELTA_T, &obs, 0.0);
        let b = libration::topocent_librations_by_diff_corrections(jde, DELTA_T, &obs);
        for t in [a, b] {
            assert_close(deg(t.long), deg(g.libr.total_long), 0.005);
            assert_close(deg(t.lat), deg(g.libr.total_lat), 0.005);
            assert_close_deg(deg(t.pos_angl_axis), deg(g.pos_angl_axis), 0.005);
        }
    }
}

// Observers north of the Moon's geocentric direction see more of its
// northern regions: the topocentric libration in latitude is larger;
// southern observers see the opposite.
#[test]
fn topocentric_latitude_direction() {
    let moon = lunar::apparent_position(JDE);
    let g = lunar::physical_ephemeris(JDE);
    let sidr = time::apprnt_sidr(
        time::mn_sidr(JDE - DELTA_T / 86400.0),
        moon.nut_long,
        moon.true_oblq,
    );
    // observers on the meridian of the Moon, at various latitudes
    for lat in [-60.0, -30.0, 0.0, 40.0, 70.0] {
        let obs = coords::GeographPoint {
            long: sidr - moon.asc,
            lat: rad(lat),
        };
        let t = libration::topocent_librations(JDE, DELTA_T, &obs, 0.0);
        let db = deg(t.lat - g.libr.total_lat);
        // δ = 13.77°: observers north of it gain, south of it lose
        if lat > 20.0 {
            assert!(db > 0.1, "lat {lat}: Δb = {db}");
        } else {
            assert!(db < -0.05, "lat {lat}: Δb = {db}");
        }
    }
}

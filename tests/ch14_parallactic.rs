// Tests for Meeus, *Astronomical Algorithms*, 2nd ed., chapter 14
// (the parallactic angle, and three other topics).

mod common;

use astro::{coords, ecliptic, misc, planet::earth};
use common::{assert_close, assert_close_deg, assert_close_rad};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

fn rad(deg: f64) -> f64 {
    deg.to_radians()
}

type V3 = [f64; 3];

fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Unit vector of (right ascension, declination).
fn unit(asc: f64, dec: f64) -> V3 {
    [dec.cos() * asc.cos(), dec.cos() * asc.sin(), dec.sin()]
}

/// Position angle of the direction `target` as seen from `s`, measured
/// from the north celestial pole through the east, computed in the plane
/// tangent to the sphere at `s`.
fn position_angle(s: V3, target: V3) -> f64 {
    let pole = [0.0, 0.0, 1.0];
    let east = cross(pole, s);
    let north = cross(s, east);
    dot(target, east).atan2(dot(target, north))
}

// ---------------------------------------------------------------------
// Parallactic angle
// ---------------------------------------------------------------------

// The parallactic angle is the position angle of the zenith (eq. 14.1);
// compare with a direct vector computation for northern and southern
// observers and all hour angles.
#[test]
fn parllc_angl_matches_position_angle_of_zenith() {
    let theta = rad(123.0); // local sidereal time, arbitrary
    for lat in (-17..=17).map(|i| i as f64 * 5.0 + 0.3) {
        let zenith = unit(theta, rad(lat));
        for h in (-18..=18).map(|i| i as f64 * 10.0 + 0.9) {
            for dec in (-17..=17).map(|i| i as f64 * 5.0 + 0.6) {
                let star = unit(theta - rad(h), rad(dec));
                let expected = position_angle(star, zenith);
                let q = misc::parllc_angl(rad(lat), rad(h), rad(dec));
                assert_close_rad(q, expected, 1e-10);
            }
        }
    }
}

#[test]
fn parllc_angl_signs_and_symmetry() {
    let lat = rad(45.0);
    // Negative before, zero at, positive after the meridian.
    assert!(misc::parllc_angl(lat, rad(-20.0), rad(10.0)) < 0.0);
    assert_close(misc::parllc_angl(lat, 0.0, rad(10.0)), 0.0, 0.0);
    assert!(misc::parllc_angl(lat, rad(20.0), rad(10.0)) > 0.0);
    // North of the zenith at upper culmination: north is down, q = 180°.
    assert_close_rad(misc::parllc_angl(lat, 0.0, rad(70.0)), PI, 1e-15);
    // Symmetric about the meridian.
    for h in (1..18).map(|i| rad(i as f64 * 10.0)) {
        for dec in (-8..=8).map(|i| rad(i as f64 * 10.0)) {
            let a = misc::parllc_angl(lat, h, dec);
            let b = misc::parllc_angl(lat, -h, dec);
            assert_close(a, -b, 1e-14);
        }
    }
    // Passing through the zenith, q jumps from −90° to +90°.
    let before = misc::parllc_angl(lat, -1e-9, lat);
    let after = misc::parllc_angl(lat, 1e-9, lat);
    assert_close(before.to_degrees(), -90.0, 1e-4);
    assert_close(after.to_degrees(), 90.0, 1e-4);
    // At the zenith it is undefined; 0 is returned, not NaN.
    assert_eq!(misc::parllc_angl(lat, 0.0, lat), 0.0);
    // At the north pole the formula stays finite.
    assert!(misc::parllc_angl(FRAC_PI_2, 1.0, 0.3).is_finite());
    // On the equator, for a body on the celestial equator, |q| = 90° except
    // at the zenith.
    assert_close(
        misc::parllc_angl(0.0, rad(45.0), 0.0).to_degrees(),
        90.0,
        1e-12,
    );
}

// When the body is on the horizon, cos q = sin φ / cos δ (Meeus, ch. 14).
#[test]
fn parllc_angl_on_hz_matches_general_formula() {
    for lat in (-13..=13).map(|i| i as f64 * 5.0 + 0.4) {
        for dec in (-5..=5).map(|i| i as f64 * 5.0 + 0.2) {
            let (lat, dec) = (rad(lat), rad(dec));
            if (lat.tan() * dec.tan()).abs() >= 1.0 {
                // Never rises or never sets.
                assert!(misc::parllc_angl_on_hz(lat, dec).is_nan());
                continue;
            }
            // Hour angle of setting for a geometric altitude of 0.
            let h0 = (-lat.tan() * dec.tan()).acos();
            let q_set = misc::parllc_angl(lat, h0, dec);
            let q_rise = misc::parllc_angl(lat, -h0, dec);
            let q_hz = misc::parllc_angl_on_hz(lat, dec);
            assert!((0.0..=PI).contains(&q_hz));
            assert_close(q_set, q_hz, 1e-12);
            assert_close(q_rise, -q_hz, 1e-12);
            // Check that h0 is really on the horizon.
            assert_close(coords::alt_frm_eq(h0, dec, lat), 0.0, 1e-12);
        }
    }
    // A circumpolar body never reaches the horizon.
    assert!(misc::parllc_angl_on_hz(rad(60.0), rad(50.0)).is_nan());
}

// ---------------------------------------------------------------------
// Ecliptic and horizon
// ---------------------------------------------------------------------

// Meeus, example 14.a: ε = 23.44°, φ = +51°, θ = 75°.
#[test]
fn example_14_a() {
    let (e, lat, theta) = (rad(23.44), rad(51.0), rad(75.0));

    let (west, east) = ecliptic::eclip_points_on_hz(e, lat, theta);
    // tan λ = −0.1879, λ = 169°21′ and 349°21′.
    assert_close(west.tan(), -0.1879, 5e-5);
    assert_close(east.tan(), -0.1879, 5e-5);
    assert_close(east.to_degrees(), 169.0 + 21.0 / 60.0, 1.0 / 60.0);
    assert_close(west.to_degrees(), 349.0 + 21.0 / 60.0, 1.0 / 60.0);

    let p = misc::eclip_hz_points(e, lat, theta);
    assert_eq!((p.west_long, p.east_long), (west, east));

    // I = 62°.
    let i = ecliptic::angl_betwn_eclip_and_hz(e, lat, theta);
    assert_close(i.to_degrees(), 62.0, 0.5);
}

// Meeus, chapter 14: for φ = +48°, ε = 23°26′, I ranges from 18°34′
// (θ = 270°) to 65°26′ (θ = 90°).
#[test]
fn angl_betwn_eclip_and_hz_extremes() {
    let (e, lat) = (rad(23.0 + 26.0 / 60.0), rad(48.0));
    let max = ecliptic::angl_betwn_eclip_and_hz(e, lat, rad(90.0));
    let min = ecliptic::angl_betwn_eclip_and_hz(e, lat, rad(270.0));
    assert_close(max.to_degrees(), 65.0 + 26.0 / 60.0, 1e-12);
    assert_close(min.to_degrees(), 18.0 + 34.0 / 60.0, 1e-12);
    for t in (0..360).map(|d| rad(d as f64)) {
        let i = ecliptic::angl_betwn_eclip_and_hz(e, lat, t);
        assert!(i >= min - 1e-12 && i <= max + 1e-12);
    }
}

// I is the angle between the pole of the ecliptic and the zenith.
#[test]
fn angl_betwn_eclip_and_hz_matches_vectors() {
    let e = rad(23.4392911);
    let pole = unit(rad(270.0), FRAC_PI_2 - e);
    for lat in (-18..=18).map(|i| i as f64 * 5.0) {
        for t in (0..36).map(|i| i as f64 * 10.0 + 4.0) {
            let zenith = unit(rad(t), rad(lat));
            let expected = dot(pole, zenith).clamp(-1.0, 1.0).acos();
            let i = ecliptic::angl_betwn_eclip_and_hz(e, rad(lat), rad(t));
            assert_close(i, expected, 1e-12);
            assert!((0.0..=PI).contains(&i));
        }
    }
}

// The first point returned is on the western horizon, the second on the
// eastern one; both are on the horizon. Checked by converting them to
// horizontal coordinates.
#[test]
fn eclip_points_on_hz_are_west_and_east() {
    for &e in &[rad(23.44), rad(23.4457889), rad(10.0)] {
        for lat in (-16..=16).map(|i| i as f64 * 5.0 + 0.5) {
            for t in (-18..54).map(|i| i as f64 * 10.0 + 2.5) {
                let (lat, t) = (rad(lat), rad(t));
                let (west, east) = ecliptic::eclip_points_on_hz(e, lat, t);
                assert!((0.0..TAU).contains(&west) && (0.0..TAU).contains(&east));
                assert_close_rad(east - west, PI, 1e-12);

                for (long, is_east) in [(west, false), (east, true)] {
                    let asc = coords::asc_frm_ecl(long, 0.0, e);
                    let dec = coords::dec_frm_ecl(long, 0.0, e);
                    let h = coords::hr_angl_frm_loc_sidr(t, asc);
                    let hz = coords::hz_frm_hr_angl_and_dec(h, dec, lat);
                    assert_close(hz.alt, 0.0, 1e-10);
                    // Azimuth from the south: the east is at negative values.
                    assert_eq!(hz.az < 0.0, is_east, "lat {lat} t {t}");
                }
            }
        }
    }
}

#[test]
fn eclip_points_on_hz_at_the_poles() {
    // At the north pole the horizon is the equator: the equinoxes.
    let (west, east) = ecliptic::eclip_points_on_hz(rad(23.44), FRAC_PI_2, rad(33.0));
    assert!(west.is_finite() && east.is_finite());
    assert_close_rad(west.min(east), 0.0, 1e-12);
    assert_close_rad(west.max(east), PI, 1e-12);
}

// ---------------------------------------------------------------------
// Ecliptic and equator
// ---------------------------------------------------------------------

// q is the position angle of the north ecliptic pole (Meeus, ch. 14).
#[test]
fn angl_betwn_poles_matches_position_angle() {
    let e = rad(23.4392911);
    let ecl_pole = unit(rad(270.0), FRAC_PI_2 - e);
    for long in (-18..54).map(|i| i as f64 * 10.0 + 1.5) {
        for lat in (-17..=17).map(|i| i as f64 * 5.0 + 0.5) {
            let (l, b) = (rad(long), rad(lat));
            let star = unit(coords::asc_frm_ecl(l, b, e), coords::dec_frm_ecl(l, b, e));
            let expected = position_angle(star, ecl_pole);
            let q = misc::angl_betwn_north_celes_and_eclip_poles(l, b, e);
            assert_close_rad(q, expected, 1e-10);
            assert!(q > -PI && q <= PI);
            // The book's tan q.
            let tan_q = l.cos() * e.tan() / (b.sin() * l.sin() * e.tan() - b.cos());
            assert_close(q.tan(), tan_q, 1e-9 * (1.0 + tan_q.abs()));
        }
    }
}

#[test]
fn angl_betwn_eclip_and_east_west_values() {
    let e = rad(23.44);
    assert_close(misc::angl_betwn_eclip_and_east_west(0.0, e), -e, 1e-15);
    assert_close(misc::angl_betwn_eclip_and_east_west(PI, e), e, 1e-15);
    assert_close(
        misc::angl_betwn_eclip_and_east_west(FRAC_PI_2, e),
        0.0,
        1e-15,
    );
    for long in (-36..72).map(|i| rad(i as f64 * 10.0 + 3.0)) {
        let q0 = misc::angl_betwn_eclip_and_east_west(long, e);
        assert!(q0.abs() <= e + 1e-15);
        // tan q₀ = −cos λ tan ε, the β = 0 case of the general angle.
        assert_close(q0.tan(), -long.cos() * e.tan(), 1e-14);
        assert_close(
            q0,
            misc::angl_betwn_north_celes_and_eclip_poles(long, 0.0, e),
            1e-14,
        );
    }
}

// ---------------------------------------------------------------------
// Diurnal path and horizon
// ---------------------------------------------------------------------

// Meeus, chapter 14: for the Sun at latitude 40° north or south, J is 50°
// at the equinoxes and 45°31′ at the solstices. Southern latitudes used
// to give 180° − J.
#[test]
fn diurnal_path_sun_at_40_degrees() {
    let e = rad(23.44);
    for &lat in &[rad(40.0), rad(-40.0)] {
        let j = earth::angl_betwn_diurnal_path_and_hz(0.0, lat);
        assert_close(j.to_degrees(), 50.0, 1e-12);
        for &dec in &[e, -e] {
            let j = earth::angl_betwn_diurnal_path_and_hz(dec, lat);
            assert_close_deg(j.to_degrees(), 45.0 + 31.0 / 60.0, 1.0 / 60.0);
        }
    }
}

#[test]
fn diurnal_path_limits_and_vectors() {
    // On the equator the paths are vertical.
    for dec in (-8..=8).map(|i| rad(i as f64 * 10.0)) {
        assert_close(
            earth::angl_betwn_diurnal_path_and_hz(dec, 0.0),
            FRAC_PI_2,
            1e-15,
        );
    }
    // δ = 0: J = 90° − |φ|.
    for lat in (-17..=17).map(|i| rad(i as f64 * 5.0)) {
        assert_close(
            earth::angl_betwn_diurnal_path_and_hz(0.0, lat),
            FRAC_PI_2 - lat.abs(),
            1e-15,
        );
    }
    // From the velocity vector of the body at rising: sin J = cos φ · C.
    for lat in (-16..=16).map(|i| rad(i as f64 * 5.0 + 0.5)) {
        for dec in (-6..=6).map(|i| rad(i as f64 * 5.0 + 0.3)) {
            if (lat.tan() * dec.tan()).abs() >= 1.0 {
                assert!(earth::angl_betwn_diurnal_path_and_hz(dec, lat).is_nan());
                continue;
            }
            let h0 = (-lat.tan() * dec.tan()).acos();
            // Direction of motion (dH > 0) in the hour-angle frame, and the
            // zenith in the same frame.
            let v = [-dec.cos() * h0.sin(), dec.cos() * h0.cos(), 0.0];
            let z = [lat.cos(), 0.0, lat.sin()];
            let norm = dot(v, v).sqrt();
            let expected = (dot(v, z).abs() / norm).asin();
            let j = earth::angl_betwn_diurnal_path_and_hz(dec, lat);
            assert_close(j, expected, 1e-12);
            assert!((0.0..=FRAC_PI_2).contains(&j));
        }
    }
    // Circumpolar: never crosses the horizon.
    assert!(earth::angl_betwn_diurnal_path_and_hz(rad(60.0), rad(50.0)).is_nan());
}

// Inside the polar circles the plain atan2 form of eq. 14.2,
// atan2(−cos θ, sin ε tan φ + cos ε sin θ), can give the eastern point
// rather than the western one; the crate's function checks which is which.
#[test]
fn eclip_points_on_hz_inside_polar_circle() {
    let (e, lat, t) = (rad(23.44), rad(-79.5), rad(32.5));
    let plain = (-t.cos())
        .atan2(e.sin() * lat.tan() + e.cos() * t.sin())
        .rem_euclid(TAU);

    let (west, east) = ecliptic::eclip_points_on_hz(e, lat, t);
    assert_close_rad(plain, east, 1e-12);

    let hz_of = |long: f64| {
        let asc = coords::asc_frm_ecl(long, 0.0, e);
        let dec = coords::dec_frm_ecl(long, 0.0, e);
        coords::hz_frm_hr_angl_and_dec(coords::hr_angl_frm_loc_sidr(t, asc), dec, lat)
    };
    assert!(hz_of(west).az > 0.0, "western point");
    assert!(hz_of(east).az < 0.0, "eastern point");
}

// Review R1: the documented range is (−π, π], but a negative zero hour
// angle (e.g. the result of rounding a tiny negative value) gave −π for
// a body on the meridian north of the zenith.
#[test]
fn parllc_angl_range_excludes_minus_pi() {
    let lat = rad(45.0);
    assert_eq!(misc::parllc_angl(lat, -0.0, rad(70.0)), PI);
    assert_eq!(misc::parllc_angl(-lat, -0.0, rad(-70.0)), 0.0);
    for &h in &[-0.0, 0.0, -PI, PI] {
        for dec in (-9..=9).map(|i| rad(i as f64 * 10.0)) {
            let q = misc::parllc_angl(lat, h, dec);
            assert!(q > -PI && q <= PI, "{h} {dec}: {q}");
        }
    }
}

// Meeus, Astronomical Algorithms, 2nd ed., chapter 24: reduction of
// ecliptical elements from one equinox to another one.

#![allow(non_snake_case)]

mod common;

use astro::precess::{self, EclElements};
use astro::time;
use common::{assert_close, assert_close_deg};
use std::f64::consts::{PI, TAU};

fn el(inc: f64, arg_perih: f64, node: f64) -> EclElements {
    EclElements {
        inc: inc.to_radians(),
        arg_perih: arg_perih.to_radians(),
        long_ascend_node: node.to_radians(),
    }
}

fn assert_elements(e: &EclElements, inc: f64, arg_perih: f64, node: f64, tol: f64) {
    assert_close(e.inc.to_degrees(), inc, tol);
    assert_close_deg(e.arg_perih.to_degrees(), arg_perih, tol);
    assert_close_deg(e.long_ascend_node.to_degrees(), node, tol);
}

fn in_range(e: &EclElements) -> bool {
    (0.0..=PI).contains(&e.inc)
        && (0.0..TAU).contains(&e.arg_perih)
        && (0.0..TAU).contains(&e.long_ascend_node)
}

/// Ecliptic (longitude, latitude) of the pole of an orbit and of its
/// perihelion
fn pole_and_perihelion(e: &EclElements) -> ((f64, f64), (f64, f64)) {
    let (si, ci) = e.inc.sin_cos();
    let (sn, cn) = e.long_ascend_node.sin_cos();
    let (sw, cw) = e.arg_perih.sin_cos();

    let pole = (e.long_ascend_node - PI / 2.0, PI / 2.0 - e.inc);

    let (x, y, z) = (cw * cn - sw * sn * ci, cw * sn + sw * cn * ci, sw * si);
    let perih = (y.atan2(x), z.atan2(x.hypot(y)));

    (pole, perih)
}

/// The elements of an orbit from the directions of its pole and perihelion
fn elements_frm_directions(pole: (f64, f64), perih: (f64, f64)) -> EclElements {
    let inc = PI / 2.0 - pole.1;
    let node = pole.0 + PI / 2.0;

    let (sl, cl) = perih.0.sin_cos();
    let (sb, cb) = perih.1.sin_cos();
    let p = (cb * cl, cb * sl, sb);

    let (si, ci) = inc.sin_cos();
    let (sn, cn) = node.sin_cos();
    let cos_w = p.0 * cn + p.1 * sn;
    let sin_w = -p.0 * sn * ci + p.1 * cn * ci + p.2 * si;

    EclElements {
        inc,
        arg_perih: sin_w.atan2(cos_w).rem_euclid(TAU),
        long_ascend_node: node.rem_euclid(TAU),
    }
}

// Meeus, example 24.a: comet Klinkenberg, from the mean equinox of
// B1744.0 to B1950.0
#[test]
fn example_24_a() {
    let JD0 = 2433282.4235 - 206.0 * 365.2421988;
    assert_close(JD0, 2358042.5305, 1e-4);
    let JD0 = 2358042.5305;
    let JD = precess::B1950;

    assert_close(time::julian_cent(JD0), -2.559958097, 1e-9);
    assert_close((JD - JD0) / 36525.0, 2.059956002, 1e-9);

    let a = precess::ecl_precess_angles(JD0, JD);
    assert_close(a.eta.to_degrees() * 3600.0, 97.0341, 0.0001);
    assert_close(a.eta.to_degrees(), 0.026954, 1e-6);
    assert_close(a.pi.to_degrees(), 172.041409, 1e-6);
    assert_close(a.p.to_degrees() * 3600.0, 10352.7137, 0.0001);
    assert_close(a.p.to_degrees(), 2.875754, 1e-6);
    let psi = a.pi + a.p;
    assert_close(psi.to_degrees(), 174.917163, 1e-6);

    let old = el(47.1220, 151.4486, 45.7481);
    let new = precess::precess_ecl_elements(&old, JD0, JD);

    // sin i = 0.73299372, i = 47.1380°, Ω - ψ = -126.313473°, Ω = 48.6037°
    assert_close(new.inc.sin(), 0.73299372, 1e-8);
    assert_close(new.inc.to_degrees(), 47.1380, 0.00005);
    assert_close_deg((new.long_ascend_node - psi).to_degrees(), -126.313473, 1e-6);
    assert_close(new.long_ascend_node.to_degrees(), 48.6037, 0.00005);

    // Δω = +0.0296°, ω = 151.4782°
    assert_close(
        (new.arg_perih - old.arg_perih).to_degrees(),
        0.0296,
        0.00005,
    );
    assert_close(new.arg_perih.to_degrees(), 151.4782, 0.00005);

    // The tuple version gives the same
    let (i, w, n) =
        precess::precess_orb_elements(old.inc, old.arg_perih, old.long_ascend_node, JD0, JD);
    assert_eq!((i, w, n), (new.inc, new.arg_perih, new.long_ascend_node));
}

// Meeus, example 24.b: comet Encke, B1950.0 -> J2000.0 with eq. 24.4
#[test]
fn example_24_b() {
    let old = el(11.93911, 186.24444, 334.04096);
    let new = precess::ecl_elements_b1950_to_j2000(&old);

    // W = 159.742178°, sin i = 0.2069767, x = 159.752866°, Δω = -0.01092°
    assert_close_deg(334.04096 - 174.298782, 159.742178, 1e-9);
    assert_close(new.inc.sin(), 0.2069767, 1e-7);
    assert_close_deg(
        new.long_ascend_node.to_degrees() - 174.997194,
        159.752866,
        1e-6,
    );
    assert_close(
        (new.arg_perih - old.arg_perih).to_degrees(),
        -0.01092,
        0.000005,
    );

    assert_elements(&new, 11.94524, 186.23352, 334.75006, 0.000005);

    // Ω increases by about 0.7° and Δω is near 0°, not 180°
    let d_node = (new.long_ascend_node - old.long_ascend_node).to_degrees();
    assert!((d_node - 0.7).abs() < 0.05);

    // The constants of eq. 24.4 are those of the general method for
    // B1950.0 -> J2000.0
    let general = precess::precess_ecl_elements(&old, precess::B1950, precess::J2000);
    assert_elements(
        &general,
        new.inc.to_degrees(),
        new.arg_perih.to_degrees(),
        new.long_ascend_node.to_degrees(),
        1e-6,
    );
    let a = precess::ecl_precess_angles(precess::B1950, precess::J2000);
    assert_close(a.eta.sin(), 0.0001139788, 1e-10);
    assert_close(a.eta.cos(), 0.9999999935, 1e-10);
    assert_close(a.pi.to_degrees(), 174.298782, 1e-6);
    assert_close((a.pi + a.p).to_degrees(), 174.997194, 1e-6);
}

// Meeus, example 24.c: comet Encke, FK4 B1950.0 -> FK5 J2000.0 (Yeomans)
#[test]
fn example_24_c() {
    let old = el(11.93911, 186.24444, 334.04096);
    let new = precess::ecl_elements_fk4_b1950_to_fk5_j2000(&old);
    assert_elements(&new, 11.94521, 186.23327, 334.75043, 0.000005);
    assert!(in_range(&new));
}

// The new elements must describe the same orbit: the pole and the
// perihelion of the orbit, precessed as points with eq. 21.7, give back
// the same elements.
#[test]
fn elements_follow_the_orbit_in_space() {
    let orbits = [
        el(47.122, 151.4486, 45.7481),
        el(11.93911, 186.24444, 334.04096),
        el(162.2627, 111.3325, 58.4201), // Halley: retrograde
        el(90.0, 10.0, 200.0),
        el(0.5, 300.0, 359.9),
        el(179.5, 45.0, 0.1),
    ];
    let epochs = [
        (precess::J2000, precess::B1950),
        (precess::J2000, precess::jd_frm_julian_epoch(-1000.0)),
        (
            precess::jd_frm_julian_epoch(1744.0),
            precess::jd_frm_julian_epoch(2300.0),
        ),
    ];

    for old in &orbits {
        for &(JD1, JD2) in &epochs {
            let new = precess::precess_ecl_elements(old, JD1, JD2);
            assert!(in_range(&new));

            let (pole, perih) = pole_and_perihelion(old);
            let pole = precess::precess_ecl_coords(pole.0, pole.1, JD1, JD2);
            let perih = precess::precess_ecl_coords(perih.0, perih.1, JD1, JD2);
            let expected = elements_frm_directions(pole, perih);

            assert_close(new.inc, expected.inc, 1e-10);
            assert_close_deg(
                new.long_ascend_node.to_degrees(),
                expected.long_ascend_node.to_degrees(),
                1e-8,
            );
            assert_close_deg(
                new.arg_perih.to_degrees(),
                expected.arg_perih.to_degrees(),
                1e-8,
            );
        }
    }
}

// Before the fix the inclination came from its sine alone, so retrograde
// orbits came out prograde (Halley: 17.7° instead of 162.2°).
#[test]
fn retrograde_orbit() {
    let halley = el(162.2627, 111.3325, 58.4201);
    let (i, _, _) = precess::precess_orb_elements(
        halley.inc,
        halley.arg_perih,
        halley.long_ascend_node,
        precess::J2000,
        precess::B1950,
    );
    assert!((i.to_degrees() - 162.26).abs() < 0.05, "{}", i.to_degrees());

    let fk5 = precess::ecl_elements_b1950_to_j2000(&halley);
    assert!(fk5.inc.to_degrees() > 162.0);
    let fk5 = precess::ecl_elements_fk4_b1950_to_fk5_j2000(&halley);
    assert!(fk5.inc.to_degrees() > 162.0);
}

#[test]
fn round_trip() {
    for old in [
        el(47.122, 151.4486, 45.7481),
        el(162.2627, 111.3325, 58.4201),
        el(3.0, 359.99, 0.01),
    ] {
        let JD = precess::jd_frm_julian_epoch(1744.0);
        let new = precess::precess_ecl_elements(&old, precess::J2000, JD);
        let back = precess::precess_ecl_elements(&new, JD, precess::J2000);
        assert_elements(
            &back,
            old.inc.to_degrees(),
            old.arg_perih.to_degrees(),
            old.long_ascend_node.to_degrees(),
            1e-6,
        );
    }
}

// Zero inclination: the new inclination is |η|, the node lies at ψ + 180°
// (ψ when η < 0, going back in time), and the perihelion keeps its
// direction in space.
#[test]
fn zero_inclination() {
    let old = el(0.0, 45.0, 123.0);
    for JD2 in [
        precess::jd_frm_julian_epoch(2300.0),
        precess::jd_frm_julian_epoch(1500.0),
    ] {
        let a = precess::ecl_precess_angles(precess::J2000, JD2);
        let new = precess::precess_ecl_elements(&old, precess::J2000, JD2);
        assert!(in_range(&new));
        assert_close(new.inc, a.eta.abs(), 1e-12);

        let psi = (a.pi + a.p).to_degrees();
        let node = if a.eta > 0.0 { psi + 180.0 } else { psi };
        assert_close_deg(new.long_ascend_node.to_degrees(), node, 1e-9);

        let (_, perih) = pole_and_perihelion(&old);
        let perih = precess::precess_ecl_coords(perih.0, perih.1, precess::J2000, JD2);
        let (_, new_perih) = pole_and_perihelion(&new);
        assert_close_deg(new_perih.0.to_degrees(), perih.0.to_degrees(), 1e-9);
        assert_close(new_perih.1, perih.1, 1e-12);
    }
}

#[test]
fn identity() {
    for old in [
        el(0.0, 45.0, 123.0),
        el(47.122, 151.4486, 45.7481),
        el(180.0, 10.0, 20.0),
    ] {
        let new = precess::precess_ecl_elements(&old, 2440000.5, 2440000.5);
        assert_elements(
            &new,
            old.inc.to_degrees(),
            old.arg_perih.to_degrees(),
            old.long_ascend_node.to_degrees(),
            1e-12,
        );
    }
}

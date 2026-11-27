//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 41: illuminated
//! fraction of the disk and magnitude of a planet. (The additions in
//! `planet::magnitude` are tested in that module.)

#![allow(non_snake_case)]

mod common;

use astro::planet::{self, Planet};
use astro::*;
use common::{assert_close, assert_close_deg};

// Example 33.a / 41.a: Venus on 1992 December 20.0 TD.
const R_VENUS: f64 = 0.724604;
const R_EARTH: f64 = 0.983824;
const DELTA_VENUS: f64 = 0.910947;

fn rad(deg: f64) -> f64 {
    deg.to_radians()
}

// Meeus, example 41.a: k = 0.647 by formula 41.2.
#[test]
fn illum_frac_frm_dist_example_41a() {
    let k = planet::illum_frac_frm_dist(R_VENUS, DELTA_VENUS, R_EARTH);
    assert_close(k, 0.647, 5e-4);
}

// Meeus, examples 41.a and 41.c: cos i = 0.29312, i = 72.96°. Before the
// fix `phase_angl` returned cos i instead of i.
#[test]
fn phase_angl_example_41a() {
    let i = planet::phase_angl(R_VENUS, DELTA_VENUS, R_EARTH);
    assert_close(i.cos(), 0.29312, 5e-6);
    assert_close(i.to_degrees(), 72.96, 0.005);
    assert_close(planet::illum_frac_frm_phase_angl(i), 0.647, 5e-4);
}

#[test]
fn phase_angl_limits_and_consistency() {
    // Planet in opposition: Sun, Earth, planet in line, i = 0.
    assert_close(planet::phase_angl(5.2, 4.2, 1.0), 0.0, 1e-7);
    // Inferior planet at inferior conjunction: i = 180°.
    assert_close(
        planet::phase_angl(0.72, 0.28, 1.0).to_degrees(),
        180.0,
        1e-5,
    );
    // Right angle: r² + Δ² = R².
    assert_close(planet::phase_angl(0.6, 0.8, 1.0).to_degrees(), 90.0, 1e-12);
    // Rounding slightly outside [-1, 1] must not give NaN.
    assert!(!planet::phase_angl(5.2, 4.2, 1.0 - 1e-17).is_nan());
    assert!(!planet::phase_angl(0.72, 0.28 + 1e-16, 1.0).is_nan());

    for &(r, delta, R) in &[
        (0.39, 0.8, 0.99),
        (1.52, 0.6, 1.01),
        (30.0, 29.5, 0.99),
        (0.72, 1.6, 1.0),
    ] {
        let i = planet::phase_angl(r, delta, R);
        assert_close(
            planet::illum_frac_frm_phase_angl(i),
            planet::illum_frac_frm_dist(r, delta, R),
            1e-12,
        );
    }
}

#[test]
fn illum_frac_frm_phase_angl_limits() {
    assert_close(planet::illum_frac_frm_phase_angl(0.0), 1.0, 1e-15);
    assert_close(planet::illum_frac_frm_phase_angl(rad(90.0)), 0.5, 1e-15);
    assert_close(planet::illum_frac_frm_phase_angl(rad(180.0)), 0.0, 1e-15);
    // Even in the phase angle.
    assert_close(
        planet::illum_frac_frm_phase_angl(rad(-40.0)),
        planet::illum_frac_frm_phase_angl(rad(40.0)),
        1e-15,
    );
}

// Meeus, example 41.c: Müller's formula gives −3.8 for Venus. Before the
// fix the phase angle (in radians) entered the degree-based formula and
// gave −4.9.
#[test]
fn apprnt_mag_muller_example_41c() {
    let i = planet::phase_angl(R_VENUS, DELTA_VENUS, R_EARTH);
    let m = planet::apprnt_mag_muller(&Planet::Venus, i, DELTA_VENUS, R_VENUS).unwrap();
    assert_close(m, -3.8, 0.05);
}

/// Müller's expressions, phase angle in degrees, written out independently.
fn muller(planet: &Planet, i: f64, delta: f64, r: f64) -> Option<f64> {
    let x = 5.0 * (r * delta).log10();
    match planet {
        Planet::Mercury => Some(1.16 + x + 0.02838 * (i - 50.0) + 0.0001023 * (i - 50.0).powi(2)),
        Planet::Venus => Some(-4.00 + x + 0.01322 * i + 0.0000004247 * i.powi(3)),
        Planet::Mars => Some(-1.30 + x + 0.01486 * i),
        Planet::Jupiter => Some(-8.93 + x),
        Planet::Uranus => Some(-6.85 + x),
        Planet::Neptune => Some(-7.05 + x),
        Planet::Earth | Planet::Saturn => None,
    }
}

/// The 1984 Astronomical Almanac expressions, phase angle in degrees.
fn aa84(planet: &Planet, i: f64, delta: f64, r: f64) -> Option<f64> {
    let x = 5.0 * (r * delta).log10();
    match planet {
        Planet::Mercury => {
            Some(-0.42 + x + 0.0380 * i - 0.000273 * i.powi(2) + 0.000002 * i.powi(3))
        }
        Planet::Venus => {
            Some(-4.40 + x + 0.0009 * i + 0.000239 * i.powi(2) - 0.00000065 * i.powi(3))
        }
        Planet::Mars => Some(-1.52 + x + 0.016 * i),
        Planet::Jupiter => Some(-9.40 + x + 0.005 * i),
        Planet::Uranus => Some(-7.19 + x),
        Planet::Neptune => Some(-6.87 + x),
        Planet::Earth | Planet::Saturn => None,
    }
}

const ALL: [Planet; 8] = [
    Planet::Mercury,
    Planet::Venus,
    Planet::Earth,
    Planet::Mars,
    Planet::Jupiter,
    Planet::Saturn,
    Planet::Uranus,
    Planet::Neptune,
];

// Every planet's formula, at several phase angles, against the expressions
// of Meeus chapter 41 (this also checks the Mercury coefficient 0.0001023,
// which was 0.000102 before the fix).
#[test]
fn apprnt_mag_formulas() {
    for planet in &ALL {
        for &i_deg in &[0.0, 12.5, 50.0, 72.96, 120.0, 150.0] {
            let (delta, r) = (0.91, 0.72);
            let muller_lib = planet::apprnt_mag_muller(planet, rad(i_deg), delta, r).ok();
            let aa84_lib = planet::apprnt_mag_84(planet, rad(i_deg), delta, r).ok();
            match (muller_lib, muller(planet, i_deg, delta, r)) {
                (Some(a), Some(b)) => assert_close(a, b, 1e-12),
                (None, None) => {}
                other => panic!("Müller mismatch: {other:?}"),
            }
            match (aa84_lib, aa84(planet, i_deg, delta, r)) {
                (Some(a), Some(b)) => assert_close(a, b, 1e-12),
                (None, None) => {}
                other => panic!("AA 1984 mismatch: {other:?}"),
            }
        }
    }
}

#[test]
fn apprnt_mag_errors_are_static() {
    let e: &'static str = planet::apprnt_mag_muller(&Planet::Earth, 0.0, 1.0, 1.0).unwrap_err();
    assert!(e.contains("Earth"));
    let e: &'static str = planet::apprnt_mag_84(&Planet::Saturn, 0.0, 1.0, 1.0).unwrap_err();
    assert!(e.contains("saturn::apprnt_mag_84"));
    assert!(planet::apprnt_mag_muller(&Planet::Saturn, 0.0, 1.0, 1.0).is_err());
    assert!(planet::apprnt_mag_84(&Planet::Earth, 0.0, 1.0, 1.0).is_err());
}

#[test]
fn apprnt_mag_distance_law() {
    // Ten times the product r Δ makes the planet five magnitudes fainter.
    for planet in &[Planet::Mars, Planet::Jupiter, Planet::Neptune] {
        let m1 = planet::apprnt_mag_84(planet, rad(10.0), 1.0, 2.0).unwrap();
        let m2 = planet::apprnt_mag_84(planet, rad(10.0), 10.0, 2.0).unwrap();
        assert_close(m2 - m1, 5.0, 1e-12);
    }
    // Jupiter at opposition at mean distances, Müller: about −2.2.
    let m = planet::apprnt_mag_muller(&Planet::Jupiter, 0.0, 4.2, 5.2).unwrap();
    assert_close(m, -2.24, 0.05);
}

// Meeus, example 41.d: Saturn on 1992 December 16.0 TD, Müller: +0.9.
// Before the fix ΔU (in radians) entered the degree-based formula, giving
// +0.76, i.e. +0.8.
#[test]
fn saturn_apprnt_mag_muller_example_41d() {
    let m = planet::saturn::apprnt_mag_muller(10.464606, 9.867882, rad(4.198), rad(16.442));
    assert_close(m, 0.9, 0.05);
    let expected = -8.68 + 5.0 * (9.867882_f64 * 10.464606).log10() + 0.044 * 4.198
        - 2.60 * rad(16.442).sin()
        + 1.25 * rad(16.442).sin().powi(2);
    assert_close(m, expected, 1e-12);
}

#[test]
fn saturn_apprnt_mag_84() {
    let (delta, r, du, b) = (10.464606, 9.867882, rad(4.198), rad(16.442));
    let m84 = planet::saturn::apprnt_mag_84(delta, r, du, b);
    let mm = planet::saturn::apprnt_mag_muller(delta, r, du, b);
    // The two formulae differ only by the constant (−8.88 against −8.68).
    assert_close(m84 - mm, -0.20, 1e-12);
    assert_close(m84, 0.7, 0.05);
}

#[test]
fn saturn_apprnt_mag_ring_symmetries() {
    let (delta, r) = (9.0, 9.5);
    let m = planet::saturn::apprnt_mag_muller(delta, r, rad(4.198), rad(16.442));
    // Only the magnitudes of ΔU and B matter.
    assert_close(
        planet::saturn::apprnt_mag_muller(delta, r, rad(-4.198), rad(16.442)),
        m,
        1e-12,
    );
    assert_close(
        planet::saturn::apprnt_mag_muller(delta, r, rad(4.198), rad(-16.442)),
        m,
        1e-12,
    );
    // ΔU given as a longitude difference across the 0/360° wrap.
    assert_close(
        planet::saturn::apprnt_mag_muller(delta, r, rad(355.802), rad(16.442)),
        m,
        1e-12,
    );
    assert_close(
        planet::saturn::apprnt_mag_84(delta, r, rad(-355.802), rad(-16.442)),
        m - 0.2,
        1e-12,
    );
    // Ring edge-on and Sun-Earth aligned: only the ball.
    assert_close(
        planet::saturn::apprnt_mag_84(delta, r, 0.0, 0.0),
        -8.88 + 5.0 * (delta * r).log10(),
        1e-12,
    );
    // An open ring makes Saturn brighter.
    assert!(
        planet::saturn::apprnt_mag_84(delta, r, 0.0, rad(26.0))
            < planet::saturn::apprnt_mag_84(delta, r, 0.0, 0.0)
    );
}

// Meeus, example 48.a (formula 48.5, which holds for planets too):
// tan χ = −0.90283 / +0.24266, χ = 285.0°. Before the fix the angle was
// returned as −75.0°.
#[test]
fn pos_angle_of_bright_limb_example_48a() {
    let sun = coords::EqPoint {
        asc: rad(20.6579),
        dec: rad(8.6964),
    };
    let moon = coords::EqPoint {
        asc: rad(134.6885),
        dec: rad(13.7684),
    };
    let chi = planet::pos_angle_of_bright_limb(sun, moon);
    assert_close(chi.tan(), -0.90283 / 0.24266, 5e-4);
    assert_close_deg(chi.to_degrees(), 285.0, 0.05);
    assert!((0.0..std::f64::consts::TAU).contains(&chi));
}

#[test]
fn pos_angle_of_bright_limb_directions() {
    let at = |asc: f64, dec: f64| coords::EqPoint {
        asc: rad(asc),
        dec: rad(dec),
    };
    // Sun to the east of the planet (larger RA): bright limb at χ = 90°.
    let chi = planet::pos_angle_of_bright_limb(at(100.0, 0.0), at(80.0, 0.0));
    assert_close_deg(chi.to_degrees(), 90.0, 1e-9);
    // Sun to the west: 270°.
    let chi = planet::pos_angle_of_bright_limb(at(60.0, 0.0), at(80.0, 0.0));
    assert_close_deg(chi.to_degrees(), 270.0, 1e-9);
    // Sun due north: 0°; due south: 180°.
    let chi = planet::pos_angle_of_bright_limb(at(80.0, 30.0), at(80.0, 10.0));
    assert_close_deg(chi.to_degrees(), 0.0, 1e-9);
    assert!((0.0..std::f64::consts::TAU).contains(&chi));
    let chi = planet::pos_angle_of_bright_limb(at(80.0, -30.0), at(80.0, -10.0));
    assert_close_deg(chi.to_degrees(), 180.0, 1e-9);
    // Across the 0/360° RA wrap: Sun at 5°, planet at 355° => Sun east.
    let chi = planet::pos_angle_of_bright_limb(at(5.0, 0.0), at(355.0, 0.0));
    assert_close_deg(chi.to_degrees(), 90.0, 1e-9);
}

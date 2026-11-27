//! Illuminated fraction of the Moon's disk (Meeus, chapter 48).

#![allow(non_snake_case)]

mod common;

use astro::lunar::{self, illumination, position};
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

// Inputs of Meeus, example 48.a
const MOON_ASC: f64 = 134.6885;
const MOON_DEC: f64 = 13.7684;
const MOON_DIST: f64 = 368410.0;
const SUN_ASC: f64 = 20.6579;
const SUN_DEC: f64 = 8.6964;
const SUN_DIST: f64 = 149971520.0;

// Meeus, example 48.a: accurate method, with intermediate values.
#[test]
fn example_48a_accurate() {
    let psi = illumination::geocent_elong_frm_eq(
        rad(SUN_ASC),
        rad(SUN_DEC),
        rad(MOON_ASC),
        rad(MOON_DEC),
    );
    assert_close(psi.cos(), -0.354991, 1e-6);
    assert_close(deg(psi), 110.7929, 1e-4);

    let i = illumination::phase_angle(psi, MOON_DIST, SUN_DIST);
    assert_close(i.tan(), 2.615404, 2e-5);
    assert_close(deg(i), 69.0756, 1e-4);

    let k = illumination::illum_frac_frm_phase_angle(i);
    assert_close(k, 0.6786, 5e-5);

    // the old API gives the same
    let sun = coords::EqPoint {
        asc: rad(SUN_ASC),
        dec: rad(SUN_DEC),
    };
    let moon = coords::EqPoint {
        asc: rad(MOON_ASC),
        dec: rad(MOON_DEC),
    };
    let k_old = lunar::illum_frac_frm_eq_coords(&sun, &moon, MOON_DIST, SUN_DIST);
    assert_close(k_old, k, 1e-15);
}

// Meeus, example 48.a: with cos i = −cos ψ, k = 0.6775.
#[test]
fn example_48a_supplement_of_elongation() {
    let psi = illumination::geocent_elong_frm_eq(
        rad(SUN_ASC),
        rad(SUN_DEC),
        rad(MOON_ASC),
        rad(MOON_DEC),
    );
    let i = illumination::phase_angle_frm_elong_approx(psi);
    assert_close(i.cos(), -psi.cos(), 1e-15);
    assert_close(illumination::illum_frac_frm_phase_angle(i), 0.6775, 5e-5);
}

// Meeus, example 48.a: approximate phase angle from D, M, M′ (eq. 48.4):
// i = 68.88°, k = 0.6802.
#[test]
fn example_48a_from_mean_arguments() {
    let i = illumination::phase_angle_approx(JDE);
    assert_close(deg(i), 68.88, 0.005);
    // The book's k = 0.6802 follows from the rounded i = 68.88°; the
    // unrounded i = 68.8834° gives k = 0.68013.
    assert_close(illumination::illum_frac_frm_phase_angle(i), 0.6802, 1e-4);
    assert_close(illumination::illum_frac_frm_phase_angle(i), 0.68013, 1e-5);
}

// Meeus, example 48.a: position angle of the bright limb.
#[test]
fn example_48a_bright_limb() {
    let chi = illumination::bright_limb_pos_angl(
        rad(SUN_ASC),
        rad(SUN_DEC),
        rad(MOON_ASC),
        rad(MOON_DEC),
    );
    // tan χ = −0.90283 / +0.24266
    let d_asc = rad(SUN_ASC - MOON_ASC);
    let num = rad(SUN_DEC).cos() * d_asc.sin();
    let den = rad(SUN_DEC).sin() * rad(MOON_DEC).cos()
        - rad(SUN_DEC).cos() * rad(MOON_DEC).sin() * d_asc.cos();
    assert_close(num, -0.90283, 5e-6);
    assert_close(den, 0.24266, 5e-6);
    assert_close_deg(deg(chi), 285.0, 0.05);
    assert!((0.0..TAU).contains(&chi));
}

// The ecliptic form of eq. 48.2, with the apparent longitudes of example
// 47.a and the Sun's longitude of example 53.b.
#[test]
fn example_48a_from_ecliptic_coords() {
    let k = lunar::illum_frac_frm_ecl_coords(
        rad(133.167265),
        rad(-3.229126),
        rad(22.33978),
        MOON_DIST,
        SUN_DIST,
    );
    assert_close(k, 0.6786, 5e-5);
    let psi = illumination::geocent_elong_frm_ecl(rad(133.167265), rad(-3.229126), rad(22.33978));
    assert_close(deg(psi), 110.7929, 2e-4);
}

// Meeus, example 48.a, everything computed from the instant.
#[test]
fn example_48a_from_jde() {
    let ill = illumination::illumination(JDE);
    assert_close(deg(ill.geocent_elong), 110.7929, 2e-4);
    assert_close(deg(ill.phase_angle), 69.0756, 2e-4);
    assert_close(ill.illum_frac, 0.6786, 5e-5);
    assert_close_deg(deg(ill.bright_limb), 285.0, 0.05);
}

#[test]
fn phase_angle_limits() {
    let (moon, sun) = (384400.0, 149597870.7);
    // full Moon: elongation 180°, phase angle 0, fully lit
    assert_close(illumination::phase_angle(PI, moon, sun), 0.0, 1e-15);
    assert_close(illumination::illum_frac_frm_phase_angle(0.0), 1.0, 0.0);
    // new Moon: elongation 0, phase angle 180°, dark
    assert_close(illumination::phase_angle(0.0, moon, sun), PI, 1e-15);
    assert_close(illumination::illum_frac_frm_phase_angle(PI), 0.0, 1e-16);
    // quarter: elongation 90° gives a phase angle slightly less than 90°
    let i = illumination::phase_angle(FRAC_PI_2, moon, sun);
    assert!(i < FRAC_PI_2);
    assert_close(deg(FRAC_PI_2 - i), deg((moon / sun).atan()), 1e-12);
    assert_close(
        illumination::illum_frac_frm_phase_angle(FRAC_PI_2),
        0.5,
        1e-16,
    );
}

// The accurate phase angle never differs from 180° − ψ by more than the
// angle subtended by the Earth–Moon distance at the Sun (≈ 0.16°), so the
// error of the cos i = −cos ψ shortcut in k stays below 0.0014.
#[test]
fn supplement_approximation_error_bound() {
    for k in 0..=1800 {
        let psi = rad(k as f64 / 10.0);
        for &moon in &[356400.0, 384400.0, 406700.0] {
            for &sun in &[1.471e8, 1.521e8] {
                let i = illumination::phase_angle(psi, moon, sun);
                let j = illumination::phase_angle_frm_elong_approx(psi);
                assert!((0.0..=PI).contains(&i));
                assert!(deg((i - j).abs()) < 0.16);
                let dk = illumination::illum_frac_frm_phase_angle(i)
                    - illumination::illum_frac_frm_phase_angle(j);
                assert!(dk.abs() < 0.0014);
            }
        }
    }
}

// Over many lunations: the three methods agree within their accuracy,
// k is in [0, 1], and the bright limb faces the Sun: east (χ near 270°)
// while the Moon waxes, west (near 90°) while it wanes.
#[test]
fn lunations() {
    let mut max_dk_approx = 0.0_f64;
    for n in 0..1500 {
        let jde = 2444000.0 + n as f64 * 1.37;
        let ill = illumination::illumination(jde);
        assert!((0.0..=1.0).contains(&ill.illum_frac));
        assert!((0.0..=PI).contains(&ill.phase_angle));
        assert!((0.0..TAU).contains(&ill.bright_limb));

        let i_approx = illumination::phase_angle_approx(jde);
        assert!((0.0..=PI).contains(&i_approx));
        let dk = illumination::illum_frac_frm_phase_angle(i_approx) - ill.illum_frac;
        max_dk_approx = max_dk_approx.max(dk.abs());

        let D = deg(position::fundamental_args(jde).mn_elong);
        let chi = deg(ill.bright_limb);
        if (15.0..165.0).contains(&D) {
            assert!(chi > 180.0, "waxing, D = {D}, χ = {chi}");
        } else if (195.0..345.0).contains(&D) {
            assert!(chi < 180.0, "waning, D = {D}, χ = {chi}");
        }
    }
    // eq. 48.4 neglects the latitude: k is still good to about 0.01
    assert!(max_dk_approx < 0.01, "{max_dk_approx}");
}

// Old ecliptic- and equatorial-coordinate functions agree with each other.
#[test]
fn ecliptic_and_equatorial_forms_agree() {
    for n in 0..300 {
        let jde = 2451545.0 + n as f64 * 2.9;
        let ill = illumination::illumination(jde);
        let moon = lunar::apparent_position(jde);
        let (sun, r) = sun::geocent_ecl_pos(jde);
        let k = lunar::illum_frac_frm_ecl_coords(
            moon.long,
            moon.lat,
            sun.long + moon.nut_long,
            moon.dist,
            r * 149597870.7,
        );
        // the ecliptic form neglects the Sun's latitude and here aberration
        assert_close(k, ill.illum_frac, 2e-4);
    }
}

// The bright-limb angle is unchanged when right ascensions are given
// modulo 2π, and the formula also works for a southern declination.
#[test]
fn bright_limb_wrap_around() {
    let a = illumination::bright_limb_pos_angl(rad(359.0), rad(-20.0), rad(1.0), rad(-20.0));
    let b = illumination::bright_limb_pos_angl(rad(-1.0), rad(-20.0), rad(361.0), rad(-20.0));
    assert_close_deg(deg(a), deg(b), 1e-9);
    // Sun 2° west of the Moon, same declination: limb towards the west
    // (position angle near 270°)
    assert_close_deg(deg(a), 270.0, 1.0);
    // and towards the east (90°) when the Sun is east of the Moon
    let c = illumination::bright_limb_pos_angl(rad(1.0), rad(-20.0), rad(359.0), rad(-20.0));
    assert_close_deg(deg(c), 90.0, 1.0);
    // Sun due north of the Moon: 0°; due south: 180°
    let n = illumination::bright_limb_pos_angl(rad(10.0), rad(30.0), rad(10.0), rad(-20.0));
    assert_close_deg(deg(n), 0.0, 1e-9);
    let s = illumination::bright_limb_pos_angl(rad(10.0), rad(-30.0), rad(10.0), rad(20.0));
    assert_close_deg(deg(s), 180.0, 1e-9);
}

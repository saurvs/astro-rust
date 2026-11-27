// Tests for Meeus, *Astronomical Algorithms*, 2nd ed., chapter 16
// (atmospheric refraction).

mod common;

use astro::atmos;
use common::assert_close;
use std::f64::consts::FRAC_PI_2;

fn rad(deg: f64) -> f64 {
    deg.to_radians()
}

fn to_arcmin(r: f64) -> f64 {
    r.to_degrees() * 60.0
}

fn to_arcsec(r: f64) -> f64 {
    r.to_degrees() * 3600.0
}

fn frm_arcmin(m: f64) -> f64 {
    (m / 60.0).to_radians()
}

/// Altitudes from `lo` to `hi` degrees in steps of 0.01°, in radians.
fn alts(lo: f64, hi: f64) -> impl Iterator<Item = f64> {
    let n = ((hi - lo) * 100.0).round() as i64;
    (0..=n).map(move |i| rad(lo + i as f64 / 100.0))
}

// Meeus, example 16.a: apparent flattening of the Sun's disk with the
// lower limb at an apparent altitude of 0°30′.
#[test]
fn example_16_a() {
    let lower_apprnt = frm_arcmin(30.0);
    let r = atmos::refrac_frm_apprnt_alt(lower_apprnt);
    assert_close(to_arcmin(r), 28.754, 5e-4);

    let lower_true = lower_apprnt - r;
    assert_close(to_arcmin(lower_true), 1.246, 5e-4);
    let upper_true = lower_true + frm_arcmin(32.0);
    assert_close(to_arcmin(upper_true), 33.246, 5e-4);
    assert_close(upper_true.to_degrees(), 0.5541, 5e-5);

    let r = atmos::refrac_frm_true_alt(upper_true);
    assert_close(to_arcmin(r), 24.618, 5e-4);
    let upper_apprnt = upper_true + r;
    assert_close(to_arcmin(upper_apprnt), 57.864, 5e-4);

    let vert_diam = upper_apprnt - lower_apprnt;
    assert_close(to_arcmin(vert_diam), 27.864, 5e-4);
    assert_close(vert_diam / frm_arcmin(32.0), 0.871, 5e-4);
}

// ---------------------------------------------------------------------
// Formulas for altitudes above 15°
// ---------------------------------------------------------------------

// Eqs. 16.1 and 16.2 at 45°, where tan z = 1. The old code used
// tan(180° − h) and the cube of the coefficient, and gave −58.29″ here.
#[test]
fn smart_formulas_at_45_degrees() {
    let r1 = atmos::refrac_frm_apprnt_alt_15(rad(45.0));
    assert_close(to_arcsec(r1), 58.294 - 0.0668, 1e-9);
    let r2 = atmos::refrac_frm_true_alt_15(rad(45.0));
    assert_close(to_arcsec(r2), 58.276 - 0.0824, 1e-9);
}

#[test]
fn smart_formulas_at_zenith_and_values() {
    assert_close(atmos::refrac_frm_apprnt_alt_15(FRAC_PI_2), 0.0, 1e-20);
    assert_close(atmos::refrac_frm_true_alt_15(FRAC_PI_2), 0.0, 1e-20);
    // At 30°: tan z = √3.
    let t = 3_f64.sqrt();
    assert_close(
        to_arcsec(atmos::refrac_frm_apprnt_alt_15(rad(30.0))),
        58.294 * t - 0.0668 * t * t * t,
        1e-9,
    );
    // Refraction decreases towards the zenith.
    let mut prev = f64::INFINITY;
    for h in alts(15.0, 90.0) {
        let r = atmos::refrac_frm_apprnt_alt_15(h);
        assert!(r <= prev && r >= 0.0);
        prev = r;
    }
}

// Eq. 16.2 was derived from 16.1, so they invert each other closely.
#[test]
fn smart_formulas_are_inverse() {
    for h0 in alts(15.0, 90.0) {
        let r1 = atmos::refrac_frm_apprnt_alt_15(h0);
        let r2 = atmos::refrac_frm_true_alt_15(h0 - r1);
        assert_close(to_arcsec(r2), to_arcsec(r1), 0.03);
    }
}

// Above 15° Bennett's formula is within its stated 0.07′ of Smart's, and
// the refined version within its stated 0.015′ (0.9″).
#[test]
fn bennett_agrees_with_smart_above_15_degrees() {
    for h0 in alts(15.0, 90.0) {
        let smart = to_arcmin(atmos::refrac_frm_apprnt_alt_15(h0));
        let bennett = to_arcmin(atmos::refrac_frm_apprnt_alt(h0));
        let refined = to_arcmin(atmos::refrac_frm_apprnt_alt_refined(h0));
        assert_close(bennett, smart, 0.07);
        assert_close(refined, smart, 0.015);
    }
}

// ---------------------------------------------------------------------
// Bennett and Sæmundsson
// ---------------------------------------------------------------------

#[test]
fn bennett_and_saemundsson_values() {
    // About 35′ on the horizon, about 1′ at 45°.
    assert_close(to_arcmin(atmos::refrac_frm_apprnt_alt(0.0)), 34.478, 1e-3);
    assert_close(
        to_arcmin(atmos::refrac_frm_apprnt_alt(rad(45.0))),
        0.995,
        1e-3,
    );
    assert_close(to_arcmin(atmos::refrac_frm_true_alt(0.0)), 28.982, 1e-3);
}

// Meeus, chapter 16: at the zenith eq. 16.3 gives −0.08″ and eq. 16.4
// −0.12″; the constants 0.0013515′ and 0.0019279′ cancel that. The
// refined Bennett formula gives −0.89″ there.
#[test]
fn zenith_values_and_corrections() {
    let rb = atmos::refrac_frm_apprnt_alt(FRAC_PI_2);
    assert_close(to_arcsec(rb), -0.08, 0.005);
    assert_close(to_arcsec(rb + atmos::BENNETT_ZENITH_CORR), 0.0, 1e-5);
    assert_close(to_arcmin(atmos::BENNETT_ZENITH_CORR), 0.0013515, 1e-15);

    let rs = atmos::refrac_frm_true_alt(FRAC_PI_2);
    assert_close(to_arcsec(rs), -0.116, 0.005);
    assert_close(to_arcsec(rs + atmos::SAEMUNDSSON_ZENITH_CORR), 0.0, 1e-5);
    assert_close(to_arcmin(atmos::SAEMUNDSSON_ZENITH_CORR), 0.0019279, 1e-15);

    let rr = atmos::refrac_frm_apprnt_alt_refined(FRAC_PI_2);
    assert_close(to_arcsec(rr), -0.89, 0.005);
}

// Sæmundsson's formula is consistent with Bennett's to within 4″
// (Meeus, chapter 16).
#[test]
fn saemundsson_inverts_bennett() {
    for h0 in alts(0.0, 90.0) {
        let rb = atmos::refrac_frm_apprnt_alt(h0);
        let rs = atmos::refrac_frm_true_alt(h0 - rb);
        assert_close(to_arcsec(rs), to_arcsec(rb), 4.0);
    }
}

#[test]
fn refined_bennett_close_to_plain() {
    for h0 in alts(0.0, 90.0) {
        let plain = to_arcmin(atmos::refrac_frm_apprnt_alt(h0));
        let refined = to_arcmin(atmos::refrac_frm_apprnt_alt_refined(h0));
        assert!((plain - refined).abs() <= 0.06 + 1e-12);
    }
}

// Down to the turning points of the formulas (−1.70° for Bennett, −1.90°
// for Sæmundsson) the refraction keeps growing as the altitude decreases.
#[test]
fn formulas_monotonic_down_to_turning_points() {
    let mut prev = -1.0;
    for h in alts(-1.69, 90.0).collect::<Vec<_>>().into_iter().rev() {
        let r = atmos::refrac_frm_apprnt_alt(h);
        assert!(r > prev, "Bennett at {}°", h.to_degrees());
        prev = r;
    }
    let mut prev = -1.0;
    for h in alts(-1.89, 90.0).collect::<Vec<_>>().into_iter().rev() {
        let r = atmos::refrac_frm_true_alt(h);
        assert!(r > prev, "Sæmundsson at {}°", h.to_degrees());
        prev = r;
    }
    // At their singularities the formulas give non-finite values, but do
    // not panic.
    let _ = atmos::refrac_frm_apprnt_alt(rad(-4.4));
    let _ = atmos::refrac_frm_true_alt(rad(-5.11));
}

// ---------------------------------------------------------------------
// Pressure and temperature
// ---------------------------------------------------------------------

#[test]
fn pressure_and_temperature_factors() {
    assert_close(atmos::refrac_by_pressr(1010.0), 1.0, 0.0);
    assert_close(atmos::refrac_by_temp(283.0), 1.0, 0.0);
    assert_close(atmos::refrac_pressr_temp_factor(1010.0, 10.0), 1.0, 1e-15);
    for &p in &[600.0, 950.0, 1010.0, 1040.0] {
        for &t in &[-30.0, 0.0, 10.0, 35.0] {
            assert_close(
                atmos::refrac_pressr_temp_factor(p, t),
                atmos::refrac_by_pressr(p) * atmos::refrac_by_temp(273.0 + t),
                1e-15,
            );
        }
    }
    // More refraction for higher pressure, less for higher temperature.
    assert!(atmos::refrac_pressr_temp_factor(1030.0, 10.0) > 1.0);
    assert!(atmos::refrac_pressr_temp_factor(1010.0, 30.0) < 1.0);
}

// ---------------------------------------------------------------------
// True <-> apparent altitude
// ---------------------------------------------------------------------

#[test]
fn alt_conversions_standard_atmosphere() {
    for h in alts(0.0, 90.0) {
        let h0 = atmos::apprnt_alt_frm_true_alt(h, 1010.0, 10.0);
        let expected = if h.to_degrees() < 90.0 {
            h + atmos::refrac_frm_true_alt(h) + atmos::SAEMUNDSSON_ZENITH_CORR
        } else {
            h
        };
        assert_close(h0, expected, 1e-15);

        let back = atmos::true_alt_frm_apprnt_alt(h0, 1010.0, 10.0);
        assert_close(to_arcsec(back), to_arcsec(h), 4.0);
    }
    // No refraction at the zenith.
    assert_close(
        atmos::apprnt_alt_frm_true_alt(FRAC_PI_2, 1010.0, 10.0),
        FRAC_PI_2,
        1e-12,
    );
    assert_close(
        atmos::true_alt_frm_apprnt_alt(FRAC_PI_2, 1010.0, 10.0),
        FRAC_PI_2,
        1e-12,
    );
}

#[test]
fn alt_conversions_scale_with_conditions() {
    let h = rad(5.0);
    let r_std = atmos::apprnt_alt_frm_true_alt(h, 1010.0, 10.0) - h;
    let r = atmos::apprnt_alt_frm_true_alt(h, 980.0, -5.0) - h;
    assert_close(
        r,
        r_std * atmos::refrac_pressr_temp_factor(980.0, -5.0),
        1e-15,
    );

    let h0 = rad(5.0);
    let r_std = h0 - atmos::true_alt_frm_apprnt_alt(h0, 1010.0, 10.0);
    let r = h0 - atmos::true_alt_frm_apprnt_alt(h0, 1025.0, 25.0);
    assert_close(
        r,
        r_std * atmos::refrac_pressr_temp_factor(1025.0, 25.0),
        1e-15,
    );
}

// The conversions are finite, continuous and increasing over the whole
// range of altitudes, including below the horizon.
#[test]
fn alt_conversions_defined_everywhere() {
    let mut prev_app = f64::NEG_INFINITY;
    let mut prev_true = f64::NEG_INFINITY;
    for h in alts(-90.0, 90.0) {
        let a = atmos::apprnt_alt_frm_true_alt(h, 1010.0, 10.0);
        let t = atmos::true_alt_frm_apprnt_alt(h, 1010.0, 10.0);
        assert!(a.is_finite() && t.is_finite());
        assert!(a > prev_app && t > prev_true, "at {}°", h.to_degrees());
        assert!(a >= h - 1e-15 && t <= h + 1e-15);
        // Continuity: steps of 0.01° in, at most 0.02° out.
        if prev_app.is_finite() {
            assert!(a - prev_app < rad(0.02) && t - prev_true < rad(0.02));
        }
        prev_app = a;
        prev_true = t;
    }
    // Below −1° the refraction is held at its value at −1°.
    let r_lim = atmos::apprnt_alt_frm_true_alt(rad(-1.0), 1010.0, 10.0) - rad(-1.0);
    let r = atmos::apprnt_alt_frm_true_alt(rad(-30.0), 1010.0, 10.0) - rad(-30.0);
    assert_close(r, r_lim, 1e-15);
    assert_close(to_arcmin(r_lim), 38.795 + 0.0019279, 1e-3);
}

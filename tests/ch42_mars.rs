//! Meeus, *Astronomical Algorithms*, 2nd ed., ch. 42: ephemeris for
//! physical observations of Mars.
//!
//! The intermediate steps of example 42.a, computed from the book's own
//! inputs, are checked by the unit tests in `src/planet/mars.rs`. These
//! tests use the public API, with the crate's VSOP87 positions.

#![allow(non_snake_case)]

mod common;

use astro::*;
use common::{arcsec, assert_close, assert_close_deg};
use std::f64::consts::TAU;

/// JDE of 1992 November 9, 0h UT, with ΔT = 59 s (example 42.a).
const JDE_42A: f64 = 2_448_935.500_683;

fn in_two_pi(a: f64) -> bool {
    (0.0..TAU).contains(&a)
}

/// Sun-Earth distance from VSOP87, for the phase-angle checks.
fn earth_sun_dist(jde: f64) -> f64 {
    planet::heliocent_coords(&planet::Planet::Earth, jde).2
}

/// Dates every 3.7 days from 1970 to 2024, covering many oppositions and
/// both Martian solstices and equinoxes.
fn sample_dates() -> impl Iterator<Item = f64> {
    (0..5400).map(|i| 2_440_000.5 + 3.7 * i as f64)
}

// Meeus, example 42.a
#[test]
fn example_42a_physical_ephemeris() {
    let e = planet::mars::physical_ephemeris(JDE_42A);

    assert_close_deg(e.earth_dec.to_degrees(), 12.44, 0.005);
    assert_close_deg(e.sun_dec.to_degrees(), -2.76, 0.005);
    assert_close_deg(e.axis_pos_angle.to_degrees(), 347.64, 0.005);
    assert_close_deg(e.central_meridian.to_degrees(), 111.55, 0.005);
    assert_close_deg(e.defect_pos_angle.to_degrees(), 279.91, 0.005);
    assert_close(e.illum_frac, 0.9012, 5e-5);
    assert_close(e.defect_of_illum.to_degrees(), arcsec(1.06), arcsec(0.005));
    assert_close(
        e.apparent_diameter.to_degrees(),
        arcsec(10.75),
        arcsec(0.005),
    );
    // Δ and τ of steps 3-4 (the book's τ = 0.005028 d is from the first
    // iteration; the converged value is 0.0050275 d)
    assert_close(e.earth_dist, 0.870_480_1, 2e-7);
    assert_close(e.sun_dist, 1.541_652_9, 2e-7);
    assert_close(e.light_time, 0.005_028, 1e-6);
}

// Meeus, example 42.a, through the older API with the book's obliquity and
// nutation (steps 10 and 13)
#[test]
fn example_42a_old_api() {
    let jc = time::julian_cent(JDE_42A);
    let e = planet::mars::ephemeris(
        JDE_42A,
        &planet::mars::north_pol_ecl_coords(jc),
        23.440_220_f64.to_radians(),
        arcsec(15.42).to_radians(),
        arcsec(-1.00).to_radians(),
    );

    assert_close_deg(e.De.to_degrees(), 12.44, 0.005);
    assert_close_deg(e.Ds.to_degrees(), -2.76, 0.005);
    assert_close_deg(e.w.to_degrees(), 111.55, 0.005);
    assert_close(e.q.to_degrees(), arcsec(1.06), arcsec(0.005));
    assert_close(e.d.to_degrees(), arcsec(10.75), arcsec(0.005));
    // Before the fix, P came out as -12.36° instead of 347.64°.
    assert!(in_two_pi(e.P), "P = {} rad", e.P);
    assert_close(e.P.to_degrees(), 347.64, 0.005);
}

// Meeus, example 42.a, step 1 (eq. 42.1)
#[test]
fn example_42a_north_pole() {
    let pole = planet::mars::north_pol_ecl_coords(time::julian_cent(JDE_42A));
    assert_close(pole.long.to_degrees(), 352.822_67, 5e-6);
    assert_close(pole.lat.to_degrees(), 63.282_08, 5e-6);

    // At J2000.0 the expressions reduce to their constant terms.
    let pole = planet::mars::north_pol_ecl_coords(0.0);
    assert_close(pole.long.to_degrees(), 352.9065, 1e-12);
    assert_close(pole.lat.to_degrees(), 63.2818, 1e-12);
}

// The ecliptic pole of eq. 42.1 must agree with the IAU J2000.0 equatorial
// pole quoted in the chapter, converted with the J2000.0 obliquity.
#[test]
fn north_pole_ecliptic_matches_equatorial_j2000() {
    let eq = planet::mars::north_pol_eq_coords_J2000();
    let ecl = planet::mars::north_pol_ecl_coords(0.0);
    let eps = ecliptic::mn_oblq_IAU(2_451_545.0);
    let asc = coords::asc_frm_ecl(ecl.long, ecl.lat, eps);
    let dec = coords::dec_frm_ecl(ecl.long, ecl.lat, eps);
    assert_close_deg(asc.to_degrees(), eq.asc.to_degrees(), 0.001);
    assert_close(dec.to_degrees(), eq.dec.to_degrees(), 0.001);

    let eq50 = planet::mars::north_pol_eq_coords_J1950();
    assert_close(eq50.asc.to_degrees(), 317.342, 1e-12);
    assert_close(eq50.dec.to_degrees(), 52.711, 1e-12);
}

// The two entry points share the algorithm: given the same obliquity and
// nutation they must give identical results.
#[test]
fn old_and_new_api_agree() {
    for jde in [2_448_935.500_683, 2_451_545.0, 2_460_000.25, 2_430_000.75] {
        let (nut_long, nut_oblq) = nutation::nutation(jde);
        let pole = planet::mars::north_pol_ecl_coords(time::julian_cent(jde));
        let old =
            planet::mars::ephemeris(jde, &pole, ecliptic::mn_oblq_IAU(jde), nut_long, nut_oblq);
        let new = planet::mars::physical_ephemeris(jde);
        assert_eq!(old.De, new.earth_dec);
        assert_eq!(old.Ds, new.sun_dec);
        assert_eq!(old.P, new.axis_pos_angle);
        assert_eq!(old.w, new.central_meridian);
        assert_eq!(old.q, new.defect_of_illum);
        assert_eq!(old.d, new.apparent_diameter);
    }
}

// Angles are reduced to their documented ranges. Before the fix the old
// API returned ω anywhere in (-180°, 540°) and P in (-180°, 180°].
#[test]
fn angles_in_documented_ranges() {
    let half_pi = std::f64::consts::FRAC_PI_2;
    for jde in sample_dates().chain((0..100).map(|i| JDE_42A + 0.02 * i as f64)) {
        let e = planet::mars::physical_ephemeris(jde);
        assert!(in_two_pi(e.axis_pos_angle), "P at {jde}");
        assert!(in_two_pi(e.central_meridian), "ω at {jde}");
        assert!(in_two_pi(e.defect_pos_angle), "Q at {jde}");
        assert!(e.earth_dec.abs() <= half_pi && e.sun_dec.abs() <= half_pi);

        let jc = time::julian_cent(jde);
        let old = planet::mars::ephemeris(
            jde,
            &planet::mars::north_pol_ecl_coords(jc),
            ecliptic::mn_oblq_IAU(jde),
            0.0,
            0.0,
        );
        assert!(in_two_pi(old.P), "old P at {jde}");
        assert!(in_two_pi(old.w), "old ω at {jde}");
    }
}

// Physical limits: D_S cannot exceed Mars' obliquity (25.19°); D_E can
// exceed it by up to about a degree because the Earth is not in Mars'
// orbital plane. The disk is never less than 84 % lit, its diameter
// stays between 3.5" and 25.2", and the defect never exceeds 16 % of it.
#[test]
fn physical_limits() {
    let mut max_ds: f64 = 0.0;
    for jde in sample_dates() {
        let e = planet::mars::physical_ephemeris(jde);
        max_ds = max_ds.max(e.sun_dec.abs().to_degrees());
        assert!(e.earth_dec.abs().to_degrees() < 26.5, "D_E at {jde}");
        assert!((0.83..=1.0).contains(&e.illum_frac), "k at {jde}");
        let d = e.apparent_diameter.to_degrees() * 3600.0;
        assert!((3.4..25.2).contains(&d), "d = {d}\" at {jde}");
        assert!(e.defect_of_illum >= 0.0 && e.defect_of_illum < 0.17 * e.apparent_diameter);
        assert!((0.37..2.7).contains(&e.earth_dist) && (1.38..1.67).contains(&e.sun_dist));
    }
    assert!(max_ds > 25.0 && max_ds < 25.25, "max |D_S| = {max_ds}");
}

// The Earth and the Sun are separated by the phase angle i as seen from
// Mars, so their planetocentric declinations can differ by at most i
// (plus the 0.007° aberration correction applied to D_S). Also checks the
// internal consistency of k, q, d and τ.
#[test]
fn consistency_relations() {
    for jde in sample_dates() {
        let e = planet::mars::physical_ephemeris(jde);
        let (r, delta) = (e.sun_dist, e.earth_dist);
        let R = earth_sun_dist(jde);
        let cos_i = (r * r + delta * delta - R * R) / (2.0 * r * delta);
        let i = cos_i.acos();
        assert!(
            (e.earth_dec - e.sun_dec).abs() <= i + 0.01_f64.to_radians(),
            "|D_E - D_S| > i at {jde}"
        );
        assert_close(e.illum_frac, (1.0 + cos_i) / 2.0, 1e-6);
        assert_close(
            e.defect_of_illum,
            (1.0 - e.illum_frac) * e.apparent_diameter,
            1e-15,
        );
        assert_close(
            e.apparent_diameter.to_degrees() * delta,
            arcsec(9.36),
            1e-12,
        );
        assert_close(e.light_time, 0.005_775_518_3 * delta, 1e-9);
    }
}

// Mars turns 350.892° a day; over one hour the central meridian advances by
// about 14.62°, less a small amount from the changing geometry.
#[test]
fn central_meridian_rate() {
    for jde in sample_dates().step_by(50) {
        let w1 = planet::mars::physical_ephemeris(jde).central_meridian;
        let w2 = planet::mars::physical_ephemeris(jde + 1.0 / 24.0).central_meridian;
        let step = (w2 - w1).rem_euclid(TAU).to_degrees();
        assert!((14.55..14.68).contains(&step), "Δω = {step}° at {jde}");
    }
}

// D_E, D_S and P change smoothly from day to day (no branch jumps).
#[test]
fn smooth_day_to_day() {
    let mut prev = planet::mars::physical_ephemeris(2_448_000.5);
    for i in 1..800 {
        let jde = 2_448_000.5 + i as f64;
        let e = planet::mars::physical_ephemeris(jde);
        assert!((e.earth_dec - prev.earth_dec).abs().to_degrees() < 0.6);
        assert!((e.sun_dec - prev.sun_dec).abs().to_degrees() < 0.3);
        let dp = (e.axis_pos_angle - prev.axis_pos_angle + std::f64::consts::PI).rem_euclid(TAU)
            - std::f64::consts::PI;
        assert!(dp.abs().to_degrees() < 1.0, "P jumps at {jde}");
        prev = e;
    }
}

// Mars' northern spring equinox (areocentric solar longitude 0°, start of
// Mars Year 37) was on 2022 December 26: D_S changes sign from negative to
// positive then.
#[test]
fn sun_dec_crosses_zero_at_northern_spring_equinox() {
    let dec_26 = 2_459_939.5; // 2022 Dec 26, 0h
    assert!(planet::mars::physical_ephemeris(dec_26 - 2.0).sun_dec < 0.0);
    assert!(planet::mars::physical_ephemeris(dec_26 + 2.0).sun_dec > 0.0);
}

// Far from 2000 (before the Gregorian reform, negative years, and the
// fourth millennium) the results remain finite and in range.
#[test]
fn distant_dates() {
    // -1000 Jan 1, 1000 Jan 1, 1582 Oct 4 (Julian), 3000 Jan 1
    for jde in [1_355_807.5, 2_086_307.5, 2_299_159.5, 2_816_787.5] {
        let e = planet::mars::physical_ephemeris(jde);
        assert!(in_two_pi(e.axis_pos_angle) && in_two_pi(e.central_meridian));
        assert!(in_two_pi(e.defect_pos_angle));
        assert!(e.earth_dec.abs().to_degrees() < 27.0 && e.sun_dec.abs().to_degrees() < 25.3);
        assert!((0.8..=1.0).contains(&e.illum_frac));
    }
}

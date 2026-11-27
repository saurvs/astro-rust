//! Meeus, *Astronomical Algorithms*, 2nd ed., ch. 45: the ring of Saturn.
//!
//! The intermediate steps of example 45.a, computed from the book's own
//! inputs, are checked by the unit tests in `src/planet/saturn/ring.rs`.
//! These tests use the public API, with the crate's VSOP87 positions.

#![allow(non_snake_case)]

mod common;

use astro::planet::saturn::ring;
use astro::*;
use common::{arcsec, assert_close, assert_close_deg};
use std::f64::consts::{PI, TAU};

/// JDE of 1992 December 16, 0h UT, with ΔT = 59 s (example 45.a).
const JDE_45A: f64 = 2_448_972.500_68;

fn in_two_pi(a: f64) -> bool {
    (0.0..TAU).contains(&a)
}

/// Every 5 days from 1960 to 2040.
fn sample_dates() -> impl Iterator<Item = f64> {
    (0..5844).map(|i| 2_436_934.5 + 5.0 * i as f64)
}

// Meeus, example 45.a
#[test]
fn example_45a_ring_elements() {
    let e = ring::ring_elements(JDE_45A);

    assert_close(e.earth_lat.to_degrees(), 16.442, 5e-4);
    assert_close(e.sun_lat.to_degrees(), 14.679, 5e-4);
    assert_close(e.long_diff.to_degrees(), 4.198, 5e-4);
    assert_close_deg(e.axis_pos_angle.to_degrees(), 6.741, 5e-4);
    assert_close(e.major_axis.to_degrees(), arcsec(35.87), arcsec(0.005));
    assert_close(e.minor_axis.to_degrees(), arcsec(10.15), arcsec(0.005));
    // steps 3-4
    assert_close(e.earth_dist, 10.464_605_9, 2e-7);
    assert_close(e.sun_dist, 9.867_881_9, 2e-7);
    assert_close(e.light_time, 0.060_44, 5e-6);
}

// Meeus, example 45.a, through the older API with the book's nutation and
// true obliquity (step 10)
#[test]
fn example_45a_old_api() {
    let e = ring::elements(
        JDE_45A,
        arcsec(16.86).to_radians(),
        (23.0 + 26.0 / 60.0 + 22.96 / 3600.0_f64).to_radians(),
    );

    assert_close(e.B.to_degrees(), 16.442, 5e-4);
    assert_close(e.B1.to_degrees(), 14.679, 5e-4);
    assert_close(e.P.to_degrees(), 6.741, 5e-4);
    assert_close(e.deltaU.to_degrees(), 4.198, 5e-4);
    assert_close(e.a.to_degrees(), arcsec(35.87), arcsec(0.005));
    assert_close(e.b.to_degrees(), arcsec(10.15), arcsec(0.005));
}

// Meeus, example 45.a, step 1 (eq. 45.1)
#[test]
fn example_45a_ring_plane() {
    let T = time::julian_cent(JDE_45A);
    assert_close(T, -0.070_431_193, 1e-9);
    assert_close(ring::inc(T).to_degrees(), 28.076_131, 5e-7);
    assert_close(ring::ascend_node(T).to_degrees(), 169.410_243, 5e-7);
    // at J2000.0 the expressions reduce to their constant terms
    assert_close(ring::inc(0.0).to_degrees(), 28.075_216, 1e-12);
    assert_close(ring::ascend_node(0.0).to_degrees(), 169.508_47, 1e-12);
}

// The axes of the other ring edges are fixed fractions of those of the
// outer edge (Meeus, ch. 45, step 6).
#[test]
fn ring_edges() {
    let e = ring::ring_elements(JDE_45A);
    let (a, b) = (e.major_axis, e.minor_axis);
    for (f, (x, y)) in [
        (0.8801, ring::inn_edge_outer_ring(a, b)),
        (0.8599, ring::out_edge_inner_ring(a, b)),
        (0.6650, ring::inn_edge_inner_ring(a, b)),
        (0.5486, ring::inn_edge_dusk_ring(a, b)),
    ] {
        assert_close(x, f * a, 1e-15);
        assert_close(y, f * b, 1e-15);
    }
    // example 45.a: 35.87" × 0.8801 = 31.57" etc.
    let (a1, _) = ring::inn_edge_outer_ring(a, b);
    assert_close(a1.to_degrees(), arcsec(31.57), arcsec(0.01));

    // The misspelt name still works and gives the same result.
    #[allow(deprecated)]
    let old = ring::out_edge_inner_ing(a, b);
    assert_eq!(old, ring::out_edge_inner_ring(a, b));
}

// The two entry points share the algorithm: given the same nutation and
// obliquity they must give identical results.
#[test]
fn old_and_new_api_agree() {
    for jde in [JDE_45A, 2_451_545.0, 2_460_000.25, 2_430_000.75] {
        let (nut_long, nut_oblq) = nutation::nutation(jde);
        let old = ring::elements(jde, nut_long, ecliptic::mn_oblq_IAU(jde) + nut_oblq);
        let new = ring::ring_elements(jde);
        assert_eq!(old.B, new.earth_lat);
        assert_eq!(old.B1, new.sun_lat);
        assert_eq!(old.P, new.axis_pos_angle);
        assert_eq!(old.deltaU, new.long_diff);
        assert_eq!(old.a, new.major_axis);
        assert_eq!(old.b, new.minor_axis);
    }
}

// ΔU is small (the book says at most about 7°). Before the fix it came out
// as 357.2° on 1995 April 13 (JD 2449820.5), when the Saturnicentric
// longitudes of the Sun and the Earth straddled the ±180° branch cut.
#[test]
fn long_diff_is_small() {
    let e = ring::elements(2_449_820.5, 0.0, 23.44_f64.to_radians());
    assert!(
        e.deltaU.to_degrees() < 7.0,
        "ΔU = {}°",
        e.deltaU.to_degrees()
    );

    let mut max: f64 = 0.0;
    for jde in sample_dates() {
        let e = ring::ring_elements(jde);
        assert!(e.long_diff >= 0.0);
        max = max.max(e.long_diff.to_degrees());
    }
    // the in-plane difference can slightly exceed the 6.5° maximum phase
    // angle, by a factor of up to 1/cos B
    assert!(max > 6.5 && max < 7.3, "max ΔU = {max}°");
}

// B, B′ and ΔU are the ring-plane latitudes and longitude difference of
// two directions (Saturn to the Earth and to the Sun) separated by the
// phase angle i, so |B - B′| ≤ i and, by the spherical law of cosines,
// cos ΔU = (cos i - sin B sin B′) / (cos B cos B′). The Sun's direction
// is corrected for aberration (≤ 0.002°), hence a small slack on i.
#[test]
fn bounded_by_phase_angle() {
    let slack = 0.003_f64.to_radians();
    let du_for = |i: f64, b1: f64, b2: f64| {
        ((i.cos() - b1.sin() * b2.sin()) / (b1.cos() * b2.cos()))
            .clamp(-1.0, 1.0)
            .acos()
    };
    for jde in sample_dates() {
        let e = ring::ring_elements(jde);
        let (r, delta) = (e.sun_dist, e.earth_dist);
        let R = planet::heliocent_coords(&planet::Planet::Earth, jde).2;
        let i = ((r * r + delta * delta - R * R) / (2.0 * r * delta)).acos();
        let (b1, b2) = (e.earth_lat, e.sun_lat);
        assert!((b1 - b2).abs() <= i + slack, "|B - B′| > i at {jde}");
        assert!(
            e.long_diff <= du_for(i + slack, b1, b2),
            "ΔU too large at {jde}"
        );
        assert!(
            e.long_diff >= du_for((i - slack).max(0.0), b1, b2),
            "ΔU too small at {jde}"
        );
    }
}

// Ranges and physical limits: |B′| cannot exceed the 26.7° tilt of the ring
// plane to Saturn's orbit; |B| can exceed it by about the phase angle.
// P stays within about ±7° of north, so it wraps around 0. Before the fix P
// came out negative half of the time.
#[test]
fn ranges_and_limits() {
    let half_pi = PI / 2.0;
    let mut max_sun_lat: f64 = 0.0;
    for jde in sample_dates() {
        let e = ring::ring_elements(jde);
        assert!(in_two_pi(e.axis_pos_angle), "P at {jde}");
        let p = (e.axis_pos_angle + PI).rem_euclid(TAU) - PI;
        assert!(
            p.abs().to_degrees() < 7.5,
            "P = {}° at {jde}",
            p.to_degrees()
        );
        assert!(e.earth_lat.abs() <= half_pi && e.sun_lat.abs() <= half_pi);
        assert!(e.earth_lat.abs().to_degrees() < 28.0);
        max_sun_lat = max_sun_lat.max(e.sun_lat.abs().to_degrees());

        let old = ring::elements(jde, 0.0, 23.44_f64.to_radians());
        assert!(in_two_pi(old.P), "old P at {jde}");
    }
    assert!(
        max_sun_lat > 26.6 && max_sun_lat < 26.8,
        "max |B′| = {max_sun_lat}"
    );
}

// b = a sin|B|, a·Δ = 375.35", and τ = 0.0057755183 Δ.
#[test]
fn consistency_relations() {
    for jde in sample_dates().step_by(7) {
        let e = ring::ring_elements(jde);
        assert_close(e.minor_axis, e.major_axis * e.earth_lat.abs().sin(), 1e-15);
        assert_close(
            e.major_axis.to_degrees() * e.earth_dist,
            arcsec(375.35),
            1e-12,
        );
        assert_close(e.light_time, 0.005_775_518_3 * e.earth_dist, 1e-9);
        assert!((7.9..11.2).contains(&e.earth_dist) && (9.0..10.1).contains(&e.sun_dist));
    }
}

// The Earth and the Sun crossed Saturn's ring plane in 1995-96 and 2009.
// Earth: 1995 May 22, 1995 Aug 11, 1996 Feb 11, 2009 Sep 4.
// Sun: 1995 Nov 19, 2009 Aug 11.
// B (or B′) must change sign within a day or two of these dates.
#[test]
fn ring_plane_crossings() {
    let changes_sign = |jde: f64, f: fn(&ring::RingElements) -> f64| {
        let before = f(&ring::ring_elements(jde - 2.0));
        let after = f(&ring::ring_elements(jde + 2.0));
        before.signum() != after.signum()
    };
    let earth = |e: &ring::RingElements| e.earth_lat;
    let sun = |e: &ring::RingElements| e.sun_lat;

    for jde in [2_449_859.5, 2_449_940.5, 2_450_124.5, 2_455_078.5] {
        assert!(changes_sign(jde, earth), "Earth crossing near {jde}");
        assert!(!changes_sign(jde, sun), "no Sun crossing near {jde}");
    }
    for jde in [2_450_040.5, 2_455_054.5] {
        assert!(changes_sign(jde, sun), "Sun crossing near {jde}");
    }
    // between the first two 1995 Earth crossings the Earth was south of the
    // ring plane while the Sun was still north of it, so the unlit face was
    // turned towards us
    let e = ring::ring_elements(2_449_900.5);
    assert!(e.earth_lat < 0.0 && e.sun_lat > 0.0);
    // and near the crossings the ring is seen edge-on
    let e = ring::ring_elements(2_449_859.5);
    assert!(e.minor_axis.to_degrees() * 3600.0 < 0.05);
}

// B, B′, ΔU and P change smoothly from day to day.
#[test]
fn smooth_day_to_day() {
    let mut prev = ring::ring_elements(2_449_000.5);
    for i in 1..1500 {
        let jde = 2_449_000.5 + i as f64;
        let e = ring::ring_elements(jde);
        assert!((e.earth_lat - prev.earth_lat).abs().to_degrees() < 0.1);
        assert!((e.sun_lat - prev.sun_lat).abs().to_degrees() < 0.02);
        assert!((e.long_diff - prev.long_diff).abs().to_degrees() < 0.15);
        let dp = (e.axis_pos_angle - prev.axis_pos_angle + PI).rem_euclid(TAU) - PI;
        assert!(dp.abs().to_degrees() < 0.05, "P jumps at {jde}");
        prev = e;
    }
}

// Far from 2000 (before the Gregorian reform, negative years, and the
// fourth millennium) the results remain finite and in range.
#[test]
fn distant_dates() {
    // -1000 Jan 1, 1000 Jan 1, 1582 Oct 4 (Julian), 3000 Jan 1
    for jde in [1_355_807.5, 2_086_307.5, 2_299_159.5, 2_816_787.5] {
        let e = ring::ring_elements(jde);
        assert!(in_two_pi(e.axis_pos_angle));
        assert!(e.earth_lat.abs().to_degrees() < 30.0 && e.sun_lat.abs().to_degrees() < 28.0);
        assert!(e.long_diff.to_degrees() < 7.5);
        assert!(e.major_axis > 0.0 && e.minor_axis >= 0.0);
    }
}

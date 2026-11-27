//! Position of the Moon (Meeus, chapter 47).

#![allow(non_snake_case)]

mod common;

use astro::lunar::{self, position};
use astro::*;
use common::{arcsec, assert_close, assert_close_deg, dms, hms};
use std::f64::consts::TAU;

/// 1992 April 12, 0h TD
const JDE: f64 = 2448724.5;

fn deg(x: f64) -> f64 {
    x.to_degrees()
}

// Meeus, example 47.a: T and the fundamental arguments.
#[test]
fn example_47a_fundamental_args() {
    assert_close(time::julian_cent(JDE), -0.077221081451, 1e-12);

    let a = position::fundamental_args(JDE);
    assert_close_deg(deg(a.mn_long), 134.290182, 1e-6);
    assert_close_deg(deg(a.mn_elong), 113.842304, 1e-6);
    assert_close_deg(deg(a.sun_mn_anom), 97.643514, 1e-6);
    assert_close_deg(deg(a.mn_anom), 5.150833, 1e-6);
    assert_close_deg(deg(a.arg_of_lat), 219.889721, 1e-6);
    assert_close_deg(deg(a.a1), 109.57, 0.005);
    assert_close_deg(deg(a.a2), 123.78, 0.005);
    assert_close_deg(deg(a.a3), 229.53, 0.005);
    assert_close(a.ecc_fac, 1.000194, 5e-7);
}

// Meeus, example 47.a: the sums of the periodic terms, including the
// additive terms. The book gives them as integers in units of 10⁻⁶ degree
// and 10⁻³ km.
#[test]
fn example_47a_periodic_sums() {
    let s = position::periodic_sums(JDE);
    assert_close(deg(s.long) * 1e6, -1127527.0, 1.0);
    assert_close(deg(s.lat) * 1e6, -3229126.0, 1.0);
    assert_close(s.dist * 1e3, -16590875.0, 1.0);
}

// Meeus, example 47.a: geometric position and parallax.
#[test]
fn example_47a_geocentric_position() {
    let p = lunar::geocentric_position(JDE);
    assert_close_deg(deg(p.long), 133.162655, 1e-6);
    assert_close(deg(p.lat), -3.229126, 1e-6);
    assert_close(p.dist, 368409.7, 0.05);
    assert_close(deg(p.eq_hz_parllx), 0.991990, 1e-6);
    assert_close(deg(p.eq_hz_parllx), dms(false, 0, 59, 31.2), arcsec(0.05));

    // the old tuple-returning API gives the same
    let (q, r) = lunar::geocent_ecl_pos(JDE);
    assert_eq!((q.long, q.lat, r), (p.long, p.lat, p.dist));
}

// Meeus, example 47.a: apparent position.
#[test]
fn example_47a_apparent_position() {
    let p = lunar::apparent_position(JDE);
    assert_close(deg(p.nut_long), 0.004610, 1e-6);
    assert_close(deg(p.nut_long) * 3600.0, 16.595, 0.01);
    assert_close_deg(deg(p.long), 133.167265, 2e-6);
    assert_close(deg(p.true_oblq), 23.440636, 2e-6);
    assert_close_deg(deg(p.asc), 134.688470, 5e-6);
    assert_close(deg(p.dec), 13.768368, 5e-6);
    assert_close_deg(deg(p.asc), hms(8, 58, 45.2), 15.0 * arcsec(0.05));
    assert_close(deg(p.dec), dms(false, 13, 46, 6.0), arcsec(0.5));
    assert_eq!(p.lat, lunar::geocentric_position(JDE).lat);
}

// Meeus, example 47.a also quotes the values of the complete ELP-2000/82
// theory; the abridged series must agree within its stated accuracy
// (10″ in longitude, 4″ in latitude).
#[test]
fn example_47a_agrees_with_full_theory() {
    let p = lunar::apparent_position(JDE);
    assert_close_deg(deg(p.long), dms(false, 133, 10, 0.0), arcsec(10.0));
    assert_close(deg(p.lat), dms(true, 3, 13, 45.0), arcsec(4.0));
    assert_close(p.dist, 368405.6, 10.0);
    assert_close_deg(deg(p.asc), hms(8, 58, 45.1), arcsec(15.0));
    assert_close(deg(p.dec), dms(false, 13, 46, 6.0), arcsec(4.0));
}

// sin π = 6378.14 / Δ (Meeus, chapter 47)
#[test]
fn eq_hz_parllx_formula() {
    assert_close(deg(lunar::eq_hz_parllx(368409.7)), 0.991990, 1e-6);
    // same quantity as the general formula of chapter 40 (sin π =
    // sin 8.794″ / Δ in AU), to the precision of the rounded constant
    // 8.794″ (6378.14 km / 1 AU = 8.79415″)
    let dist_au = 368409.7 / 149597870.7;
    assert_close(
        deg(lunar::eq_hz_parllx(368409.7)),
        deg(parallax::eq_hz_parallax(dist_au)),
        arcsec(0.07),
    );
    for k in 0..100 {
        let p = lunar::geocentric_position(JDE + k as f64 * 0.77);
        assert_eq!(p.eq_hz_parllx, lunar::eq_hz_parllx(p.dist));
    }
    assert!(lunar::eq_hz_parllx(1000.0).is_nan());
}

// Ranges and continuity over a long span, including dates before 1582 and
// negative years.
#[test]
fn ranges_and_continuity() {
    for &start in &[
        0.0,       // -4712 Jan 1.5
        1355807.5, // -1000
        2299159.5, // 1582 Oct 4/15
        2451545.0, // J2000
        2816787.5, // 3000
    ] {
        let mut prev = lunar::geocentric_position(start);
        for h in 1..(24 * 40) {
            let p = lunar::geocentric_position(start + h as f64 / 24.0);
            assert!((0.0..TAU).contains(&p.long));
            assert!(deg(p.lat).abs() < 5.35, "lat {}", deg(p.lat));
            assert!(p.dist > 355_000.0 && p.dist < 407_500.0, "dist {}", p.dist);
            // the Moon moves 0.45°–0.68° per hour in longitude
            let dl = (p.long - prev.long).rem_euclid(TAU).to_degrees();
            assert!(dl > 0.4 && dl < 0.75, "step {dl} at {start} + {h} h");
            assert!((p.lat - prev.lat).abs().to_degrees() < 0.1);
            assert!((p.dist - prev.dist).abs() < 300.0);
            prev = p;
        }
    }
}

// The Moon's longitude advances by 360° in a tropical month of 27.32158
// days on average.
#[test]
fn mean_motion() {
    let a0 = position::fundamental_args(2451545.0);
    assert_close_deg(deg(a0.mn_long), 218.3164477, 1e-7);
    assert_close_deg(deg(a0.mn_elong), 297.8501921, 1e-7);
    assert_close_deg(deg(a0.sun_mn_anom), 357.5291092, 1e-7);
    assert_close_deg(deg(a0.mn_anom), 134.9633964, 1e-7);
    assert_close_deg(deg(a0.arg_of_lat), 93.2720950, 1e-7);
    assert_close(a0.ecc_fac, 1.0, 0.0);

    let a1 = position::fundamental_args(2451545.0 + 27.321582);
    assert_close_deg(deg(a1.mn_long), 218.3164477, 0.001);
}

// Longitude of the mean ascending node, eq. 47.7.
#[test]
fn mean_ascending_node() {
    // value at J2000.0 and in example 53.a
    assert_close_deg(
        deg(lunar::mean_ascending_node(2451545.0)),
        125.0445479,
        1e-7,
    );
    assert_close_deg(deg(lunar::mean_ascending_node(JDE)), 274.400656, 1e-6);

    // the old function takes centuries and gives the same values
    for k in -50..50 {
        let jde = 2451545.0 + k as f64 * 3652.5;
        let T = time::julian_cent(jde);
        assert_eq!(lunar::mn_ascend_node(T), lunar::mean_ascending_node(jde));
        let n = lunar::mean_ascending_node(jde);
        assert!((0.0..TAU).contains(&n));
    }

    // regression of the nodes in 18.6 years (6798.38 days), westward
    let n0 = lunar::mean_ascending_node(2451545.0);
    let n1 = lunar::mean_ascending_node(2451545.0 + 6798.38);
    assert_close_deg(deg(n1), deg(n0), 0.01);
    let n2 = lunar::mean_ascending_node(2451545.0 + 10.0);
    assert_close_deg(deg(n0 - n2), 10.0 * 1934.1362891 / 36525.0, 1e-6);

    // Ω equals L′ − F up to a constant 0.000195° (same secular terms)
    for k in -20..20 {
        let jde = 2451545.0 + k as f64 * 3652.5;
        let a = position::fundamental_args(jde);
        let n = lunar::mean_ascending_node(jde);
        assert_close_deg(deg(n), deg(a.mn_long - a.arg_of_lat) + 0.0001952, 1e-6);
    }
}

// Longitude of the true ascending node.
#[test]
fn true_ascending_node() {
    // 1992 April 12: the five periodic terms, evaluated independently from
    // the book's D, M, M′, F of example 47.a, sum to −0.663128°.
    assert_close_deg(
        deg(lunar::true_ascending_node(JDE)),
        274.400656 - 0.663128,
        2e-6,
    );

    for k in 0..2000 {
        let jde = 2440000.0 + k as f64 * 3.7;
        let mean = lunar::mean_ascending_node(jde);
        let tru = lunar::true_ascending_node(jde);
        assert!((0.0..TAU).contains(&tru));
        // the sum of the amplitudes is 1.9682°
        assert_close_deg(deg(tru), deg(mean), 1.9683);
        assert_eq!(lunar::true_ascend_node(time::julian_cent(jde)), tru);
    }
}

// Longitude of the mean perigee.
#[test]
fn mean_perigee() {
    assert_close_deg(deg(lunar::mean_perigee(2451545.0)), 83.3532465, 1e-7);
    // one revolution of the line of apsides in 8.85 years (T = 0.0885),
    // apart from the quadratic term
    let T = 360.0 / 4069.0137287;
    let p0 = lunar::mean_perigee(2451545.0);
    let p1 = lunar::mean_perigee(2451545.0 + T * 36525.0);
    assert_close_deg(deg(p1), deg(p0) - 0.01032 * T * T, 1e-8);

    // Π equals L′ − M′ up to a constant 0.000195° (same secular terms)
    for k in -20..20 {
        let jde = 2451545.0 + k as f64 * 3652.5;
        let a = position::fundamental_args(jde);
        let p = lunar::mean_perigee(jde);
        assert!((0.0..TAU).contains(&p));
        assert_close_deg(deg(p), deg(a.mn_long - a.mn_anom) + 0.0001952, 1e-6);
        assert_eq!(lunar::mn_perigee(time::julian_cent(jde)), p);
    }
}

// The apparent position is the geometric one plus nutation, and α, δ are
// consistent with λ, β.
#[test]
fn apparent_position_consistency() {
    for k in 0..200 {
        let jde = 2415020.0 + k as f64 * 365.2422 + 0.37 * k as f64;
        let g = lunar::geocentric_position(jde);
        let a = lunar::apparent_position(jde);
        let (dpsi, _) = nutation::nutation(jde);
        assert_close_deg(deg(a.long), deg(g.long + dpsi), 1e-10);
        assert!((0.0..TAU).contains(&a.asc));
        let long = coords::ecl_long_frm_eq(a.asc, a.dec, a.true_oblq);
        let lat = coords::ecl_lat_frm_eq(a.asc, a.dec, a.true_oblq);
        assert_close_deg(deg(long), deg(a.long), 1e-9);
        assert_close(deg(lat), deg(a.lat), 1e-9);
        assert!(deg(a.dec).abs() < 23.45 + 5.35);
    }
}

//! Regression tests for bugs fixed in the pre-existing lunar API
//! (Meeus, chapters 47, 48 and 53). Every test here uses only functions that
//! existed before the fixes, and failed against the old code.

#![allow(non_snake_case)]

mod common;

use astro::*;
use common::{assert_close, assert_close_deg};
use std::f64::consts::TAU;

const JDE_1992_APR_12: f64 = 2448724.5;

fn eq(asc_deg: f64, dec_deg: f64) -> coords::EqPoint {
    coords::EqPoint {
        asc: asc_deg.to_radians(),
        dec: dec_deg.to_radians(),
    }
}

// Bug: `geocent_ecl_pos` returned L′ + Σl without reducing it, so the
// longitude could be slightly negative (or ≥ 2π) when L′ was near 0.
#[test]
fn geocent_ecl_pos_longitude_is_normalised() {
    // hourly over two months covers several passages of L′ through 0
    for h in 0..(24 * 60) {
        let jde = 2451545.0 + h as f64 / 24.0;
        let (p, _) = lunar::geocent_ecl_pos(jde);
        assert!(
            (0.0..TAU).contains(&p.long),
            "longitude {} at JDE {jde}",
            p.long
        );
    }
}

// Bug: `true_ascend_node` added the periodic terms to the mean node without
// reducing the sum, so it left [0, 2π) when the mean node was near 0.
#[test]
fn true_ascend_node_is_normalised() {
    // the mean node passes through 0 at about 1987.85 and 2006.46
    for d in 0..2000 {
        let jde = 2446700.0 + d as f64 * 0.5;
        let n = lunar::true_ascend_node(time::julian_cent(jde));
        assert!((0.0..TAU).contains(&n), "node {n} at JDE {jde}");
        let jde = 2453500.0 + d as f64 * 0.5;
        let n = lunar::true_ascend_node(time::julian_cent(jde));
        assert!((0.0..TAU).contains(&n), "node {n} at JDE {jde}");
    }
}

// Bug: the optical libration in longitude was reduced to [0, 2π), giving
// 358.794° instead of the book's −1.206° (Meeus, example 53.a).
#[test]
fn optical_libr_is_signed() {
    let (l, b) = lunar::optical_libr(
        JDE_1992_APR_12,
        133.162655_f64.to_radians(),
        -3.229126_f64.to_radians(),
    );
    assert_close(l.to_degrees(), -1.206, 5e-4);
    assert_close(b.to_degrees(), 4.194, 5e-4);
}

// Bug: `total_libr` passed its arguments to `optical_libr` and
// `physical_libr` in the wrong order (longitude as the Julian day, etc.).
// Meeus, example 53.a: l = −1.23°, b = +4.20°.
#[test]
fn total_libr_example_53a() {
    let (l, b) = lunar::total_libr(
        JDE_1992_APR_12,
        133.162655_f64.to_radians(),
        -3.229126_f64.to_radians(),
    );
    assert_close(l.to_degrees(), -1.23, 5e-3);
    assert_close(b.to_degrees(), 4.20, 5e-3);
}

// Bug: the position angle of the axis came straight from asin, i.e. in
// (−90°, 90°], while it must be given in [0°, 360°) (first or fourth
// quadrant).
#[test]
fn pos_angl_of_axis_of_rot_is_normalised() {
    let mut seen_fourth_quadrant = false;
    for k in 0..400 {
        let jde = 2448000.0 + k as f64 * 3.3;
        let (p, _) = lunar::geocent_ecl_pos(jde);
        let (nut_long, nut_oblq) = nutation::nutation(jde);
        let eps = ecliptic::mn_oblq_IAU(jde) + nut_oblq;
        let asc = coords::asc_frm_ecl(p.long + nut_long, p.lat, eps);
        let (_, b) = lunar::total_libr(jde, p.long, p.lat);
        let node = lunar::mn_ascend_node(time::julian_cent(jde));
        let P = lunar::pos_angl_of_axis_of_rot(jde, node, b, nut_long, eps, asc);
        assert!((0.0..TAU).contains(&P), "P = {P} at JDE {jde}");
        seen_fourth_quadrant |= P > std::f64::consts::PI;
    }
    assert!(seen_fourth_quadrant);
}

// Bug: the correction to the libration in latitude had the wrong sign
// (Meeus gives Δb = +π′ cos(Q − P)). For an observer in the northern
// hemisphere looking at the Moon on the meridian near the celestial
// equator (Q = 0), with the axis pointing north (P = 0), the observer sees
// further over the Moon's north limb: Δb must be positive.
#[test]
fn topocent_diff_corrections_sign_of_delta_b() {
    let par = 0.95_f64.to_radians();
    let (dl, db, _) =
        lunar::topocent_libr_by_diff_corrections(45_f64.to_radians(), 0.0, 0.0, par, 0.0, 0.0);
    let z = 45_f64.to_radians();
    let pi1 = par * (z.sin() + 0.0084 * (2.0 * z).sin());
    assert_close(db, pi1, 1e-15);
    assert_close(dl, 0.0, 1e-15);
}

// Bug: `illum_frac_frm_eq_coords` took the arc cosine of the angular
// separation, which already is an angle. Meeus, example 48.a: k = 0.6786.
#[test]
fn illum_frac_frm_eq_coords_example_48a() {
    let k = lunar::illum_frac_frm_eq_coords(
        &eq(20.6579, 8.6964),
        &eq(134.6885, 13.7684),
        368410.0,
        149971520.0,
    );
    assert_close(k, 0.6786, 5e-5);
}

// Bug: the position angle of the bright limb came straight from atan2,
// i.e. in (−180°, 180°], instead of [0°, 360°). Meeus, example 48.a:
// χ = 285.0°.
#[test]
fn bright_limb_example_48a() {
    let chi = lunar::bright_limb(eq(20.6579, 8.6964), eq(134.6885, 13.7684));
    assert!((0.0..TAU).contains(&chi));
    assert_close_deg(chi.to_degrees(), 285.0, 0.05);
}

//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 31: mean elements of
//! the planetary orbits (tables 31.A and 31.B).

#![allow(non_snake_case)]

mod common;

use astro::planet::{self, OrbElements, Planet};
use common::assert_close;

const PLANETS: [Planet; 8] = [
    Planet::Mercury,
    Planet::Venus,
    Planet::Earth,
    Planet::Mars,
    Planet::Jupiter,
    Planet::Saturn,
    Planet::Uranus,
    Planet::Neptune,
];

const J2000: f64 = 2451545.0;

fn jd_frm_cent(t: f64) -> f64 {
    J2000 + 36525.0 * t
}

/// Signed difference of two angles in degrees, in (-180, 180].
fn diff_deg(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(360.0);
    if d > 180.0 { d - 360.0 } else { d }
}

/// Absolute difference of two angles in radians, modulo 2π.
fn angle_diff_rad(a: f64, b: f64) -> f64 {
    let tau = std::f64::consts::TAU;
    let d = (a - b).rem_euclid(tau);
    d.min(tau - d)
}

fn assert_angles_in_range(el: &OrbElements) {
    let tau = std::f64::consts::TAU;
    for x in [
        el.mn_long,
        el.long_ascend_node,
        el.long_perih,
        el.mn_anom,
        el.arg_perih,
    ] {
        assert!((0.0..tau).contains(&x), "{el:?}");
    }
}

// Meeus, example 31.a: Mercury on 2065 June 24 at 0h TD.
#[test]
fn mean_elements_of_date_example_31a() {
    let jd = 2475460.5;
    assert_close(astro::time::julian_cent(jd), 0.654770704997, 1e-12);

    let el = planet::mean_orb_elements_of_date(&Planet::Mercury, jd);
    assert_close(el.mn_long.to_degrees(), 203.494701, 1e-6);
    assert_close(el.semimaj_axis, 0.387098310, 1e-9);
    assert_close(el.ecc, 0.20564510, 1e-8);
    assert_close(el.inc.to_degrees(), 7.006171, 1e-6);
    assert_close(el.long_ascend_node.to_degrees(), 49.107650, 1e-6);
    assert_close(el.long_perih.to_degrees(), 78.475382, 1e-6);
    assert_close(el.mn_anom.to_degrees(), 125.019319, 1e-6);
    assert_close(el.arg_perih.to_degrees(), 29.367732, 1e-6);
}

// The tuple API must give the same numbers. Before the fix the semimajor
// axis of Mercury was 0.0387 AU instead of 0.387 AU.
#[test]
fn orb_elements_tuple_example_31a() {
    let (L, a, e, i, omega, pi, M, w) = planet::orb_elements(&Planet::Mercury, 2475460.5);
    assert_close(L.to_degrees(), 203.494701, 1e-6);
    assert_close(a, 0.387098310, 1e-9);
    assert_close(e, 0.20564510, 1e-8);
    assert_close(i.to_degrees(), 7.006171, 1e-6);
    assert_close(omega.to_degrees(), 49.107650, 1e-6);
    assert_close(pi.to_degrees(), 78.475382, 1e-6);
    assert_close(M.to_degrees(), 125.019319, 1e-6);
    assert_close(w.to_degrees(), 29.367732, 1e-6);
}

#[test]
fn orb_elements_tuple_matches_struct() {
    for p in PLANETS {
        for t in [-20.0, -1.3, 0.0, 0.25, 7.5] {
            let jd = jd_frm_cent(t);
            let el = planet::mean_orb_elements_of_date(&p, jd);
            let tup = planet::orb_elements(&p, jd);
            assert_eq!(
                tup,
                (
                    el.mn_long,
                    el.semimaj_axis,
                    el.ecc,
                    el.inc,
                    el.long_ascend_node,
                    el.long_perih,
                    el.mn_anom,
                    el.arg_perih,
                ),
                "{p:?} T = {t}"
            );
        }
    }
}

// Semimajor axes and eccentricities of table 31.A (they do not depend on
// the reference frame). Evaluating away from J2000.0 checks every power of
// T: before the fix several eccentricities used T³ in place of T², one
// a3 sign was wrong and Mercury's semimajor axis was off by a factor 10.
#[test]
fn semimaj_axis_and_ecc_polynomials() {
    #[rustfmt::skip]
    let table: [(Planet, [f64; 4], [f64; 4]); 8] = [
        (Planet::Mercury, [0.387098310, 0.0, 0.0, 0.0], [0.20563175, 0.000020407, -0.0000000283, -0.00000000018]),
        (Planet::Venus, [0.723329820, 0.0, 0.0, 0.0], [0.00677192, -0.000047765, 0.0000000981, 0.00000000046]),
        (Planet::Earth, [1.000001018, 0.0, 0.0, 0.0], [0.01670863, -0.000042037, -0.0000001267, 0.00000000014]),
        (Planet::Mars, [1.523679342, 0.0, 0.0, 0.0], [0.09340065, 0.000090484, -0.0000000806, -0.00000000025]),
        (Planet::Jupiter, [5.202603209, 0.0000001913, 0.0, 0.0], [0.04849793, 0.000163225, -0.0000004714, -0.00000000201]),
        (Planet::Saturn, [9.554909192, -0.0000021390, 0.000000004, 0.0], [0.05554814, -0.000346641, -0.0000006436, 0.0000000034]),
        (Planet::Uranus, [19.218446062, -0.0000000372, 0.00000000098, 0.0], [0.04638122, -0.000027293, 0.0000000789, 0.00000000024]),
        (Planet::Neptune, [30.110386869, -0.0000001663, 0.00000000069, 0.0], [0.00945575, 0.000006033, 0.0, -0.00000000005]),
    ];
    let poly = |c: &[f64; 4], t: f64| c[0] + c[1] * t + c[2] * t * t + c[3] * t * t * t;

    for (p, a, e) in table {
        for t in [-10.0, -3.0, 0.0, 2.0, 10.0] {
            let jd = jd_frm_cent(t);
            for el in [
                planet::mean_orb_elements_of_date(&p, jd),
                planet::mean_orb_elements_j2000(&p, jd),
            ] {
                assert_close(el.semimaj_axis, poly(&a, t), 1e-12);
                assert_close(el.ecc, poly(&e, t), 1e-12);
            }
        }
    }
}

// At J2000.0 both tables give the same angles (except the Earth's node,
// which is undefined in the ecliptic of date).
#[test]
fn tables_agree_at_j2000() {
    for p in PLANETS {
        let a = planet::mean_orb_elements_of_date(&p, J2000);
        let b = planet::mean_orb_elements_j2000(&p, J2000);
        assert_eq!(a.mn_long, b.mn_long, "{p:?}");
        assert_eq!(a.inc, b.inc, "{p:?}");
        assert_eq!(a.long_perih, b.long_perih, "{p:?}");
        assert_eq!(a.mn_anom, b.mn_anom, "{p:?}");
        if p != Planet::Earth {
            assert_eq!(a.long_ascend_node, b.long_ascend_node, "{p:?}");
        }
    }
}

// The mean anomaly is measured from the perihelion and does not depend on
// the equinox: the difference L - ϖ must be the same in both tables. This
// checks the L and ϖ rows of tables 31.A and 31.B against each other (the
// residual comes from the rounding of the printed coefficients).
#[test]
fn mean_anomaly_is_frame_independent() {
    for p in PLANETS {
        for k in -20..=20 {
            let jd = jd_frm_cent(k as f64 * 0.5);
            let a = planet::mean_orb_elements_of_date(&p, jd);
            let b = planet::mean_orb_elements_j2000(&p, jd);
            let d = diff_deg(a.mn_anom.to_degrees(), b.mn_anom.to_degrees());
            assert!(d.abs() < 5e-6, "{p:?} T = {}: {d}", k as f64 * 0.5);
        }
    }
}

// Precessing the J2000.0 elements (table 31.B) to the ecliptic and equinox
// of date (Meeus, ch. 24) must give the elements of table 31.A. The small
// residuals come from the different precession models.
#[test]
fn j2000_elements_precess_to_elements_of_date() {
    for p in PLANETS {
        if p == Planet::Earth {
            continue;
        }
        for k in -20..=20 {
            let t = k as f64 * 0.5;
            let jd = jd_frm_cent(t);
            let date = planet::mean_orb_elements_of_date(&p, jd);
            let fixed = planet::mean_orb_elements_j2000(&p, jd);
            let (inc, arg_perih, node) = astro::precess::precess_orb_elements(
                fixed.inc,
                fixed.arg_perih,
                fixed.long_ascend_node,
                J2000,
                jd,
            );
            let (tol_inc, tol_node) = if t.abs() <= 1.0 {
                (1e-5, 3e-4)
            } else {
                (1e-4, 1e-2)
            };
            let d_inc = inc.to_degrees() - date.inc.to_degrees();
            let d_node = diff_deg(node.to_degrees(), date.long_ascend_node.to_degrees());
            let d_arg = diff_deg(arg_perih.to_degrees(), date.arg_perih.to_degrees());
            assert!(d_inc.abs() < tol_inc, "{p:?} T = {t}: Δi = {d_inc}");
            assert!(d_node.abs() < tol_node, "{p:?} T = {t}: ΔΩ = {d_node}");
            assert!(d_arg.abs() < tol_node, "{p:?} T = {t}: Δω = {d_arg}");
        }
    }
}

// Meeus, ch. 31: on the fixed ecliptic of J2000.0 the inclination of the
// Earth's orbit is counted negative before 2000, so that i and Ω vary
// continuously through the epoch.
#[test]
fn earth_j2000_inclination_changes_sign_at_j2000() {
    let before = planet::mean_orb_elements_j2000(&Planet::Earth, jd_frm_cent(-0.01));
    let at = planet::mean_orb_elements_j2000(&Planet::Earth, J2000);
    let after = planet::mean_orb_elements_j2000(&Planet::Earth, jd_frm_cent(0.01));

    assert!(before.inc < 0.0);
    assert_eq!(at.inc, 0.0);
    assert!(after.inc > 0.0);
    assert_close(after.inc.to_degrees(), 0.0130548 * 0.01, 1e-9);
    assert_close(at.long_ascend_node.to_degrees(), 174.873176, 1e-9);
    assert!(
        diff_deg(
            after.long_ascend_node.to_degrees(),
            before.long_ascend_node.to_degrees()
        )
        .abs()
            < 0.01
    );
}

// In the ecliptic of date the Earth's node is undefined: by convention it
// is 0 and the argument of the perihelion equals the longitude of the
// perihelion.
#[test]
fn earth_of_date_conventions() {
    for t in [-5.0, 0.0, 3.0] {
        let el = planet::mean_orb_elements_of_date(&Planet::Earth, jd_frm_cent(t));
        assert_eq!(el.inc, 0.0);
        assert_eq!(el.long_ascend_node, 0.0);
        assert_close(el.arg_perih, el.long_perih, 1e-12);
    }
}

#[test]
fn derived_elements_and_ranges() {
    let tau = std::f64::consts::TAU;
    for p in PLANETS {
        // From about 1000 B.C. to A.D. 5000, including negative years.
        for k in -30..=30 {
            let jd = jd_frm_cent(k as f64) + 0.3;
            for el in [
                planet::mean_orb_elements_of_date(&p, jd),
                planet::mean_orb_elements_j2000(&p, jd),
            ] {
                assert_angles_in_range(&el);
                let m = (el.mn_long - el.long_perih).rem_euclid(tau);
                let w = (el.long_perih - el.long_ascend_node).rem_euclid(tau);
                assert!(angle_diff_rad(el.mn_anom, m) < 1e-9);
                assert!(angle_diff_rad(el.arg_perih, w) < 1e-9);
                assert!(el.ecc > 0.0 && el.ecc < 0.25, "{p:?} {el:?}");
                assert!(el.inc.abs() < 8f64.to_radians(), "{p:?} {el:?}");
            }
        }
    }
}

// The mean longitude increases by one revolution per sidereal period: the
// a1 coefficient of L (degrees per Julian century) gives the mean motion.
#[test]
fn mean_longitude_advances_with_mean_motion() {
    // Sidereal periods in days, for a rough sanity check.
    let periods = [
        (Planet::Mercury, 87.969),
        (Planet::Venus, 224.701),
        (Planet::Earth, 365.256),
        (Planet::Mars, 686.980),
    ];
    for (p, period) in periods {
        let a = planet::mean_orb_elements_j2000(&p, J2000);
        let b = planet::mean_orb_elements_j2000(&p, J2000 + period);
        let d = diff_deg(b.mn_long.to_degrees(), a.mn_long.to_degrees());
        assert!(d.abs() < 0.01, "{p:?}: {d}");
    }
}

/// Solves Kepler's equation by Newton's method.
fn ecc_anom(m: f64, e: f64) -> f64 {
    let mut ea = m;
    for _ in 0..50 {
        ea -= (ea - e * ea.sin() - m) / (1.0 - e * ea.cos());
    }
    ea
}

// A Keplerian orbit with the mean elements of date tracks the complete
// VSOP87 position, apart from the periodic perturbations that the mean
// elements leave out (largest for the giant planets: the great inequality
// of Jupiter and Saturn). This checks every element of table 31.A against
// an independent theory, and heliocent_coords_frm_orb_elements as well.
#[test]
fn kepler_orbit_from_mean_elements_follows_vsop87() {
    // Maximum differences in longitude and latitude (degrees) over ±10
    // centuries; observed values are about 70% of these.
    let tolerances = [
        (Planet::Mercury, 0.015, 0.003),
        (Planet::Venus, 0.015, 0.003),
        (Planet::Earth, 0.015, 0.003),
        (Planet::Mars, 0.06, 0.003),
        (Planet::Jupiter, 0.6, 0.02),
        (Planet::Saturn, 1.6, 0.08),
        (Planet::Uranus, 1.3, 0.02),
        (Planet::Neptune, 0.9, 0.03),
    ];
    for (p, tol_long, tol_lat) in tolerances {
        for k in -100..=100 {
            let jd = jd_frm_cent(k as f64 * 0.1) + 0.37 * k as f64;
            let el = planet::mean_orb_elements_of_date(&p, jd);
            let ea = ecc_anom(el.mn_anom, el.ecc);
            let v = 2.0 * (((1.0 + el.ecc) / (1.0 - el.ecc)).sqrt() * (ea / 2.0).tan()).atan();
            let r = el.semimaj_axis * (1.0 - el.ecc * ea.cos());
            let (l, b) = planet::heliocent_coords_frm_orb_elements(
                el.inc,
                el.long_ascend_node,
                el.arg_perih,
                v,
                r,
            );
            let (L, B, R) = planet::heliocent_coords(&p, jd);
            let d_long = diff_deg(l.to_degrees(), L.to_degrees());
            let d_lat = b.to_degrees() - B.to_degrees();
            assert!(d_long.abs() < tol_long, "{p:?} JD {jd}: Δl = {d_long}");
            assert!(d_lat.abs() < tol_lat, "{p:?} JD {jd}: Δb = {d_lat}");
            assert!(((r - R) / el.semimaj_axis).abs() < 0.006, "{p:?} JD {jd}");
        }
    }
}

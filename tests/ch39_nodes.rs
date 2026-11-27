// Meeus, Astronomical Algorithms, 2nd ed., chapter 39: passages through the
// nodes

#![allow(non_snake_case)]

mod common;

use astro::orbit::{self, Elements, Motion, Node, elliptic, parabolic};
use astro::time;
use common::{assert_close, assert_close_deg};
use std::f64::consts::{PI, TAU};

fn rad(deg: f64) -> f64 {
    deg.to_radians()
}

fn jd(year: i16, month: time::Month, day: f64) -> f64 {
    time::julian_day(&time::Date {
        year,
        month,
        decimal_day: day,
        cal_type: time::CalType::Gregorian,
    })
}

// Meeus, example 39.a: comet Halley, elements by W. Landgraf (1950.0)
const HALLEY_W: f64 = 111.84644;
const HALLEY_E: f64 = 0.96727426;
const HALLEY_N: f64 = 0.01297082; // degrees per day
const HALLEY_A: f64 = 17.9400782;

fn halley_T() -> f64 {
    jd(1986, time::Month::Feb, 9.45891)
}

#[test]
fn example_39a_ascending_node() {
    let e = HALLEY_E;
    // v = 360° - ω or -ω
    let v = rad(360.0 - HALLEY_W);
    assert_close(v.to_degrees(), 248.15356, 1e-9);
    assert_close(
        ((1.0 - e) / (1.0 + e)).sqrt() * (v / 2.0).tan(),
        -0.1906646,
        1e-7,
    );
    let E = elliptic::ecc_anom_frm_true_anom(v, e);
    assert_close(E.to_degrees(), -21.5894332, 1e-7);
    let M = elliptic::mn_anom_frm_ecc_anom(E, e);
    assert_close(M.to_degrees(), -1.1972043, 1e-7);
    assert_close(M.to_degrees() / HALLEY_N, -92.2998, 1e-4);
    // -ω gives the same result as 360° - ω
    assert_close(
        elliptic::ecc_anom_frm_true_anom(rad(-HALLEY_W), e),
        E,
        1e-15,
    );

    let T = halley_T();
    let (t, r) =
        elliptic::passage_through_node(rad(HALLEY_W), rad(HALLEY_N), HALLEY_A, e, T, &Node::Ascend);
    assert_close(t - T, -92.2998, 1e-4);
    // 1985 November 9.16 TD
    assert_close(t, jd(1985, time::Month::Nov, 9.16), 0.005);
    assert_close(r, 1.8045, 1e-4);
}

#[test]
fn example_39a_descending_node() {
    let e = HALLEY_E;
    let v = rad(180.0 - HALLEY_W);
    assert_close(v.to_degrees(), 68.15356, 1e-9);
    let E = elliptic::ecc_anom_frm_true_anom(v, e);
    assert_close(E.to_degrees(), 9.9726067, 1e-7);
    let M = elliptic::mn_anom_frm_ecc_anom(E, e);
    assert_close(M.to_degrees(), 0.3749928, 1e-7);

    let T = halley_T();
    let (t, r) = elliptic::passage_through_node(
        rad(HALLEY_W),
        rad(HALLEY_N),
        HALLEY_A,
        e,
        T,
        &Node::Descend,
    );
    assert_close(t - T, 28.9105, 1e-4);
    // 1986 March 10.37 TD
    assert_close(t, jd(1986, time::Month::Mar, 10.37), 0.005);
    assert_close(r, 0.8493, 1e-4);
}

#[test]
fn example_39b_parabolic_orbit() {
    // Meeus, example 39.b: comet Helin-Roman (1989 IX)
    let T = jd(1989, time::Month::Aug, 20.291);
    let (w, q) = (rad(154.9103), 1.324502);

    // ascending node: v = -154.9103°, s = -4.4940577
    assert_close((-w / 2.0).tan(), -4.4940577, 1e-7);
    let (t, r) = parabolic::passage_through_node(w, q, T, &Node::Ascend);
    // The book prints T - 4354.66 days; the formula gives 4354.651 (the
    // date, 1977 September 17, agrees either way).
    assert_close(t - T, -4354.65, 0.005);
    assert_close(t, jd(1977, time::Month::Sept, 17.64), 0.01);
    assert_close(r, 28.07, 0.005);

    // descending node: v = +25.0897°, s = +0.2225161
    assert_close(((PI - w) / 2.0).tan(), 0.2225161, 1e-7);
    let (t, r) = parabolic::passage_through_node(w, q, T, &Node::Descend);
    assert_close(t - T, 28.3454, 1e-4);
    assert_close(t, jd(1989, time::Month::Sept, 17.636), 1e-3);
    assert_close(r, 1.3901, 1e-4);
}

#[test]
fn example_39c_venus() {
    // Meeus, example 39.c: Venus at the ascending node nearest to 1979.0
    let a = 0.723329820;
    let n = elliptic::mn_motion(a);
    assert_close(n.to_degrees(), 1.602137, 1e-6);

    // e and ω for T = -0.21 from the polynomials quoted in the example
    let Tc: f64 = -0.21;
    let e = 0.00677192 - 0.000047765 * Tc + 0.0000000981 * Tc * Tc;
    let w_deg = 54.883783 + 0.5011082 * Tc - 0.0014824 * Tc * Tc;
    assert_close(e, 0.00678195, 1e-8);
    assert_close(w_deg, 54.778485, 1e-6);
    // continue, as the book does, with the rounded values
    let (e, w_deg) = (0.00678195, 54.778485);

    let v = rad(-w_deg);
    let E = elliptic::ecc_anom_frm_true_anom(v, e);
    assert_close(E.to_degrees(), -54.461662, 1e-6);
    let M = elliptic::mn_anom_frm_ecc_anom(E, e);
    assert_close(M.to_degrees(), -54.145467, 1e-6);

    // T = 1978 December 31.204 (example 38.a): node on 1978 Nov. 27.408
    let T = jd(1978, time::Month::Dec, 31.204);
    let (t, _) = elliptic::passage_through_node(rad(w_deg), n, a, e, T, &Node::Ascend);
    assert_close(t - T, -33.7958, 1e-4);
    assert_close(t, jd(1978, time::Month::Nov, 27.408), 1e-3);
}

#[test]
fn true_and_eccentric_anomaly_round_trip() {
    for &e in &[0.0, 0.0167, 0.5, 0.967, 0.9999] {
        for i in -179..=180 {
            let v = rad(i as f64);
            let E = elliptic::ecc_anom_frm_true_anom(v, e);
            assert!(E > -PI && E <= PI);
            assert_close(elliptic::true_anom(E, e), v, 1e-12);
            let M = elliptic::mn_anom_frm_ecc_anom(E, e);
            assert_close(elliptic::ecc_anom_checked(M, e).unwrap(), E, 1e-12);
        }
    }
}

/// Heliocentric ecliptic latitude of a body at a given time
fn latitude(el: &Elements, t: f64) -> f64 {
    let (x, y, z) = el.heliocent_ecl_rect(t).unwrap();
    z.atan2(x.hypot(y))
}

fn longitude(el: &Elements, t: f64) -> f64 {
    let (x, y, _) = el.heliocent_ecl_rect(t).unwrap();
    y.atan2(x).rem_euclid(TAU)
}

#[test]
fn body_is_in_the_ecliptic_at_the_nodes_elliptic() {
    // At the computed times the body's heliocentric latitude is zero and
    // its longitude is Ω (ascending node) or Ω + 180° (descending node)
    let T = halley_T();
    let node = rad(58.14397);
    for w_deg in [HALLEY_W, 0.0, 45.0, 179.0, 181.0, 300.0, -60.0, 400.0] {
        let w = rad(w_deg);
        let n = elliptic::mn_motion(HALLEY_A);
        let el = Elements {
            motion: Motion::Elliptic {
                semimaj_axis: HALLEY_A,
                ecc: HALLEY_E,
                mn_motion: n,
                mn_anom: 0.0,
                epoch: T,
            },
            incl: rad(162.23932),
            arg_perih: w,
            long_asc_node: node,
        };
        let period = TAU / n;
        for (which, long, sign) in [(Node::Ascend, node, 1.0), (Node::Descend, node + PI, -1.0)] {
            let (t, r) = elliptic::passage_through_node(w, n, HALLEY_A, HALLEY_E, T, &which);
            assert!((t - T).abs() <= period / 2.0 + 1e-6, "ω = {w_deg}");
            assert_close(latitude(&el, t), 0.0, 1e-9);
            assert_close_deg(longitude(&el, t).to_degrees(), long.to_degrees(), 1e-7);
            let (_, rv) = el.motion.true_anom_and_rad_vec(t).unwrap();
            assert_close(r, rv, 1e-9 * r);
            // the latitude changes sign there: northward at the ascending
            // node, southward at the descending one
            let dt = 1e-3;
            assert!(sign * (latitude(&el, t + dt) - latitude(&el, t - dt)) > 0.0);
        }
    }
}

#[test]
fn body_is_in_the_ecliptic_at_the_nodes_parabolic() {
    let T = jd(1989, time::Month::Aug, 20.291);
    let q = 1.324502;
    let node = rad(83.3);
    for w_deg in [154.9103, 10.0, 90.0, 200.0, 350.0, -20.0] {
        let w = rad(w_deg);
        let el = Elements {
            motion: Motion::Parabolic {
                perih_dist: q,
                perih_time: T,
            },
            incl: rad(9.7),
            arg_perih: w,
            long_asc_node: node,
        };
        for (which, long) in [(Node::Ascend, node), (Node::Descend, node + PI)] {
            let (t, r) = parabolic::passage_through_node(w, q, T, &which);
            assert_close(latitude(&el, t), 0.0, 1e-9);
            assert_close_deg(longitude(&el, t).to_degrees(), long.to_degrees(), 1e-7);
            let (_, rv) = el.motion.true_anom_and_rad_vec(t).unwrap();
            assert_close(r, rv, 1e-9 * r);
        }
    }
}

#[test]
fn passages_are_periodic_in_the_argument_of_perihelion() {
    let T = halley_T();
    let n = rad(HALLEY_N);
    for which in [Node::Ascend, Node::Descend] {
        for w in [0.3, 2.0, 4.0, 6.0] {
            let a = elliptic::passage_through_node(w, n, HALLEY_A, HALLEY_E, T, &which);
            let b = elliptic::passage_through_node(w + TAU, n, HALLEY_A, HALLEY_E, T, &which);
            let c = elliptic::passage_through_node(w - TAU, n, HALLEY_A, HALLEY_E, T, &which);
            assert_close(b.0, a.0, 1e-6);
            assert_close(c.0, a.0, 1e-6);
            assert_close(b.1, a.1, 1e-9);
            let p = parabolic::passage_through_node(w, 1.3, T, &which);
            let p2 = parabolic::passage_through_node(w + TAU, 1.3, T, &which);
            assert_close(p2.0, p.0, 1e-6 * (p.0 - T).abs().max(1.0));
        }
    }
}

#[test]
fn circular_orbit_nodes() {
    // e = 0: the node passages are ω/n before and (π - ω)/n after the
    // perihelion, half a period apart
    let (w, n, a, T) = (rad(40.0), rad(1.0), 1.0, 2451545.0);
    let (ta, ra) = elliptic::passage_through_node(w, n, a, 0.0, T, &Node::Ascend);
    let (td, rd) = elliptic::passage_through_node(w, n, a, 0.0, T, &Node::Descend);
    assert_close(ta, T - 40.0, 1e-9);
    assert_close(td, T + 140.0, 1e-9);
    assert_close(td - ta, 180.0, 1e-9);
    assert_eq!((ra, rd), (1.0, 1.0));
}

#[test]
fn node_enum_traits() {
    let n = Node::Ascend;
    let m = n;
    assert_eq!(n, m);
    assert_ne!(Node::Ascend, Node::Descend);
    assert_eq!(format!("{:?}", Node::Descend), "Descend");
    let _ = orbit::Node::Ascend;
}

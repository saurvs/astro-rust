//! Meeus, chapter 51: passages of the Moon through the nodes.

#![allow(non_snake_case)]

mod common;
mod lunar_events_support;

use astro::lunar::node::{self, Node};
use astro::{lunar, time};
use common::assert_close;
use lunar_events_support as support;
use time::Month::*;

const DRACONIC_MONTH: f64 = 27.212220817;

fn date(year: i16, month: time::Month, day: f64) -> time::Date {
    time::Date {
        year,
        month,
        decimal_day: day,
        cal_type: time::CalType::Gregorian,
    }
}

// Meeus, example 51.a: passage through the ascending node in May 1987.
#[test]
fn example_51a_ascending_node_1987_may() {
    // Mid-May 1987 = 1987.37 gives k ~ -170.19, hence k = -170.
    assert_close(node::k_from_year(1987.37), -170.19, 0.005);
    let jde = node::node_passage(-170, Node::Ascending);
    assert_close(jde, 2446938.76803, 0.000005);
    // 1987 May 23.26803 TD = 6h26m TD
    assert_close(
        (jde - support::jd(1987, May, 23.0)) * 24.0,
        6.0 + 26.0 / 60.0,
        0.01,
    );

    let p = node::nearest_node_passage(support::jd(1987, May, 15.0), Node::Ascending);
    assert_eq!(p.number, -170);
    assert_eq!(p.k(), -170.0);
    assert_eq!(p.jde, jde);

    let (ascending, descending) = lunar::time_of_passage_through_nodes(&date(1987, May, 15.0));
    assert_close(ascending, 2446938.76803, 0.000005);
    // The nearest descending passage is the one before (k = -170.5).
    assert_eq!(descending, node::node_passage(-171, Node::Descending));
}

#[test]
fn numbering() {
    // k = 0 is the ascending node passage of 2000 January 21.
    let jde = node::node_passage(0, Node::Ascending);
    assert!((jde - support::jd(2000, Jan, 21.0)).abs() < 1.0);
    // Meeus' examples of admissible k values.
    for (number, n, k) in [
        (223, Node::Ascending, 223.0),
        (-147, Node::Ascending, -147.0),
        (223, Node::Descending, 223.5),
        (-47, Node::Descending, -46.5),
    ] {
        let p = node::nearest_node_passage(node::node_passage(number, n), n);
        assert_eq!((p.number, p.node, p.k()), (number, n, k));
        assert_close(node::k_from_jd(p.jde), k, 0.05);
    }
}

// Both elements of the tuple are the passages nearest to the date, each
// within half a draconic month of it.
#[test]
fn time_of_passage_through_nodes_returns_nearest_passages() {
    for (y, m, d) in [
        (1987, May, 15.0),
        (2000, Jan, 1.0),
        (1582, Oct, 15.0),
        (-500, Mar, 3.5),
        (2100, Dec, 31.0),
    ] {
        let jd = support::jd(y, m, d);
        let (a, b) = lunar::time_of_passage_through_nodes(&date(y, m, d));
        assert_eq!(a, node::nearest_node_passage(jd, Node::Ascending).jde);
        assert_eq!(b, node::nearest_node_passage(jd, Node::Descending).jde);
        assert!((a - jd).abs() < DRACONIC_MONTH / 2.0 + 0.5);
        assert!((b - jd).abs() < DRACONIC_MONTH / 2.0 + 0.5);
    }
}

#[test]
fn nearest_passage_is_never_farther_than_its_neighbours() {
    let start = support::jd(1985, Jan, 1.0);
    let mut jd = start;
    while jd < start + 400.0 {
        for n in [Node::Ascending, Node::Descending] {
            let p = node::nearest_node_passage(jd, n);
            let d = (p.jde - jd).abs();
            for k in [p.number - 1, p.number + 1] {
                assert!(d <= (node::node_passage(k, n) - jd).abs());
            }
        }
        jd += 0.31;
    }
    assert!(
        node::nearest_node_passage(f64::NAN, Node::Ascending)
            .jde
            .is_nan()
    );
}

#[test]
fn node_passages_between_alternate() {
    let start = support::jd(2023, Jan, 1.0);
    let list: Vec<_> = node::node_passages_between(start, start + 365.25).collect();
    assert!((26..=28).contains(&list.len()));
    for w in list.windows(2) {
        assert_ne!(w[0].node, w[1].node);
        let gap = w[1].jde - w[0].jde;
        assert!((12.0..15.2).contains(&gap), "{gap}");
    }
    let e = list[0];
    assert_eq!(
        node::node_passages_between(e.jde, e.jde + 1.0).collect::<Vec<_>>(),
        [e]
    );
    assert_eq!(node::node_passages_between(e.jde - 1.0, e.jde).count(), 0);
    assert_eq!(node::node_passages_between(f64::NAN, 1.0).count(), 0);
    assert_eq!(
        node::node_passages_between(0.0, f64::INFINITY)
            .take(5)
            .count(),
        5
    );
}

/// Instant at which the chapter 47 latitude of the Moon is zero, near
/// `jde`, and the sign of the latitude just after it.
fn ch47_node(jde: f64) -> (f64, f64) {
    let lat = |t: f64| support::moon(t).1;
    let t = support::bisect(jde - 0.3, jde + 0.3, lat);
    (t, lat(t + 0.1).signum())
}

// Independent check against the latitude from the chapter 47 series,
// 1980-2020. Meeus quotes largest errors of 142 s (ascending) and 132 s
// (descending) against an accurate theory.
#[test]
fn passages_agree_with_chapter_47_latitude() {
    for n in (-270..270).step_by(2) {
        for nd in [Node::Ascending, Node::Descending] {
            let jde = node::node_passage(n, nd);
            let (t, sign_after) = ch47_node(jde);
            let dt = (t - jde) * 86400.0;
            assert!(dt.abs() < 160.0, "{n} {nd:?}: {dt} s");
            let expected = if nd == Node::Ascending { 1.0 } else { -1.0 };
            assert_eq!(sign_after, expected);
        }
    }
}

// The terms of the series that contain the Sun's mean anomaly M must be
// multiplied by E. Versions before 3.0 did not; the effect grows away from
// 2000. Over the years 0-40 the passages with E agree clearly better with
// the chapter 47 latitude than the same series without E.
#[test]
fn factor_e_improves_agreement_far_from_2000() {
    let rad = |deg: f64| deg.rem_euclid(360.0).to_radians();
    let (mut with_e, mut without_e, mut count) = (0.0, 0.0, 0.0);
    for n in (-26900..-26360).step_by(3) {
        for nd in [Node::Ascending, Node::Descending] {
            let jde = node::node_passage(n, nd);
            let k = n as f64 + nd.fraction();
            let T = k / 1342.23;
            let D = rad(183.6380
                + 331.73735682 * k
                + T * T * (0.0014852 + T * (0.00000209 - T * 0.000000010)));
            let M = rad(17.4006 + 26.82037250 * k + T * T * (0.0001186 + T * 0.00000006));
            let M1 = rad(38.3776
                + 355.52747313 * k
                + T * T * (0.0123499 + T * (0.000014627 - T * 0.000000069)));
            let E = 1.0 - T * (0.002516 + T * 0.0000074);
            let (D2, D4) = (2.0 * D, 4.0 * D);
            let solar = -0.0083 * (D2 - M).sin() - 0.0039 * (D2 - M - M1).sin()
                + 0.0030 * (D2 + M).sin()
                + 0.0028 * (M - M1).sin()
                + 0.0026 * M.sin()
                + 0.0022 * (M + M1).sin()
                + 0.0005 * (D2 + M - M1).sin()
                + 0.0004 * (D2 - M + M1).sin()
                - 0.0003 * (D2 - 2.0 * M).sin()
                + 0.0003 * (D4 - M).sin();
            let no_e = jde - (E - 1.0) * solar;
            let (t, _) = ch47_node(jde);
            with_e += ((t - jde) * 86400.0).powi(2);
            without_e += ((t - no_e) * 86400.0).powi(2);
            count += 1.0;
        }
    }
    let (rms_e, rms_no_e) = ((with_e / count).sqrt(), (without_e / count).sqrt());
    assert!(rms_e < rms_no_e - 4.0, "{rms_e} s vs {rms_no_e} s");
    // Meeus quotes errors up to 144 s for the years 0-40.
    assert!(rms_e < 60.0);
}

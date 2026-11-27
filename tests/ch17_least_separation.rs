// Meeus, chapter 17: least angular separation of two moving bodies.

mod common;

use astro::{angle, conjunction, interpol};
use common::{assert_close, dms, hms};

/// Positions of a body at five instants, `(α, δ)` in radians
type Table = [(f64, f64); 5];

fn p(h: u32, m: u32, s: f64, neg: bool, d: u32, dm: u32, ds: f64) -> (f64, f64) {
    (hms(h, m, s).to_radians(), dms(neg, d, dm, ds).to_radians())
}

/// Mercury and Saturn at 0h TD on 1978 September 12 to 16 (Meeus, ch. 17)
fn mercury_saturn() -> (Table, Table) {
    let mercury = [
        p(10, 23, 17.65, false, 11, 31, 46.3),
        p(10, 29, 44.27, false, 11, 2, 5.9),
        p(10, 36, 19.63, false, 10, 29, 51.7),
        p(10, 43, 1.75, false, 9, 55, 16.7),
        p(10, 49, 48.85, false, 9, 18, 34.7),
    ];
    let saturn = [
        p(10, 33, 1.23, false, 10, 42, 53.5),
        p(10, 33, 29.64, false, 10, 40, 13.2),
        p(10, 33, 57.97, false, 10, 37, 33.4),
        p(10, 34, 26.22, false, 10, 34, 53.9),
        p(10, 34, 54.39, false, 10, 32, 14.9),
    ];
    (mercury, saturn)
}

/// Separation of the positions interpolated (eq. 3.8) at `n`, interpolating
/// δ₁, δ₂ and α₁ − α₂ as the book prescribes.
fn interpolated_sepr(b1: &Table, b2: &Table, n: f64) -> f64 {
    let f = |g: &dyn Fn(usize) -> f64| interpol::five_values(g(0), g(1), g(2), g(3), g(4), n);
    let d_a = f(&|i| b1[i].0 - b2[i].0);
    let d1 = f(&|i| b1[i].1);
    let d2 = f(&|i| b2[i].1);
    angle::anglr_sepr(d_a, d1, 0.0, d2)
}

// Meeus, ch. 17: the distances at the five tabulated instants
#[test]
fn mercury_saturn_tabulated_distances() {
    let (m, s) = mercury_saturn();
    let expected = [2.5211, 0.9917, 0.5943, 2.2145, 3.8710];
    for i in 0..5 {
        let d = angle::anglr_sepr(m[i].0, m[i].1, s[i].0, s[i].1);
        assert_close(d.to_degrees(), expected[i], 5e-5);
    }
}

// Meeus, ch. 17: distances from the interpolated coordinates
#[test]
fn mercury_saturn_interpolated_distances() {
    let (m, s) = mercury_saturn();
    let table = [
        (-0.50, 0.21437),
        (-0.45, 0.14057),
        (-0.40, 0.07790),
        (-0.35, 0.07028),
        (-0.30, 0.12815),
        (-0.38, 0.06408),
        (-0.37, 0.06229),
        (-0.36, 0.06448),
    ];
    for &(n, d) in &table {
        assert_close(interpolated_sepr(&m, &s, n).to_degrees(), d, 5e-6);
    }
}

// Meeus, ch. 17, exercise: least separation 0°03′44″ (0.06228°) on 1978
// September 13 at 15h06.5m TD, n = −0.370502 from five positions.
#[test]
fn mercury_saturn_least_separation_five_positions() {
    let (m, s) = mercury_saturn();
    let r = conjunction::least_angular_separation(&m, &s).unwrap();
    // The book's n comes from a parabola through distances rounded to
    // 1e-5°, which limits it to about ±1e-5.
    assert_close(r.n, -0.370502, 1e-5);
    assert_close(14.0 + r.n, 13.629498, 1e-5);
    assert_close(r.sepr.to_degrees(), 0.06228, 5e-6);
    assert_close(r.sepr.to_degrees() * 3600.0, 224.0, 0.5);

    // It is the minimum of the interpolated distance
    for i in -40..=40 {
        let n = r.n + i as f64 * 1e-3;
        assert!(interpolated_sepr(&m, &s, n) >= r.sepr - 1e-15);
    }
    // and its derivative vanishes there: the minimum is located to ~1e-9
    let h = 1e-6;
    let lo = interpolated_sepr(&m, &s, r.n - h);
    let hi = interpolated_sepr(&m, &s, r.n + h);
    assert!((lo - hi).abs() < 1e-14, "{lo} {hi}");
}

// Meeus, ch. 17: Danjon's method with three positions (September 13, 14,
// 15) gives n = −0.370727, u = +70.20″, v = −212.42″, distance 224″.
#[test]
fn mercury_saturn_danjon_three_positions() {
    let (m, s) = mercury_saturn();
    let r = conjunction::least_angular_separation_danjon(&m[1..4], &s[1..4]).unwrap();
    assert_close(r.n, -0.370727, 1e-6);
    assert_close(r.sepr.to_degrees() * 3600.0, 70.20_f64.hypot(-212.42), 0.01);
    assert_close(r.sepr.to_degrees() * 3600.0, 224.0, 0.5);
    // Differs from the five-position value by about 0.000225 day (19 s)
    let r5 = conjunction::least_angular_separation(&m, &s).unwrap();
    assert_close(r5.n - r.n, 0.000225, 1e-5);
}

#[test]
fn exact_and_danjon_methods_agree_for_close_approaches() {
    let (m, s) = mercury_saturn();
    let e = conjunction::least_angular_separation(&m, &s).unwrap();
    let d = conjunction::least_angular_separation_danjon(&m, &s).unwrap();
    assert_close(e.n, d.n, 1e-5);
    assert_close(
        e.sepr.to_degrees() * 3600.0,
        d.sepr.to_degrees() * 3600.0,
        0.01,
    );
    let e3 = conjunction::least_angular_separation(&m[1..4], &s[1..4]).unwrap();
    assert_close(e3.n, -0.3707, 1e-4);
    assert_close(e3.sepr.to_degrees() * 3600.0, 224.0, 0.5);
}

#[test]
fn least_separation_symmetries() {
    let (m, s) = mercury_saturn();
    let r = conjunction::least_angular_separation(&m, &s).unwrap();

    // Swapping the bodies
    let q = conjunction::least_angular_separation(&s, &m).unwrap();
    assert_close(q.n, r.n, 1e-12);
    assert_close(q.sepr, r.sepr, 1e-15);

    // Rotating all right ascensions, so that both bodies cross 0h/24h
    // within the table, and wrapping them to [0, 2π)
    let shift =
        |b: &Table| b.map(|(a, d)| (angle::limit_to_two_PI(a - (159.0_f64).to_radians()), d));
    let q = conjunction::least_angular_separation(&shift(&m), &shift(&s)).unwrap();
    assert_close(q.n, r.n, 1e-12);
    assert_close(q.sepr, r.sepr, 1e-14);
    let q = conjunction::least_angular_separation_danjon(&shift(&m), &shift(&s)).unwrap();
    assert_close(q.n, -0.3705, 1e-4);

    // Mirror image in the southern sky
    let south = |b: &Table| b.map(|(a, d)| (a, -d));
    let q = conjunction::least_angular_separation(&south(&m), &south(&s)).unwrap();
    assert_close(q.n, r.n, 1e-12);
    assert_close(q.sepr, r.sepr, 1e-15);

    // Time reversal: n changes sign
    let rev = |b: &Table| [b[4], b[3], b[2], b[1], b[0]];
    let q = conjunction::least_angular_separation(&rev(&m), &rev(&s)).unwrap();
    assert_close(q.n, -r.n, 1e-12);
}

// A body moving uniformly in right ascension along the equator past a fixed
// star 0.002 rad north of it: the closest approach is exactly when the RAs
// are equal, at a separation of exactly 0.002 rad.
#[test]
fn uniform_motion_past_a_star() {
    let star = (1.0, 0.002);
    for &n0 in &[-0.9, -0.3, 0.0, 0.25, 0.8] {
        let body: Vec<(f64, f64)> = (-1..=1)
            .map(|i| (1.0 + 0.01 * (i as f64 - n0), 0.0))
            .collect();
        let r = conjunction::least_angular_separation(&body, &[star]).unwrap();
        assert_close(r.n, n0, 1e-12);
        assert_close(r.sepr, 0.002, 1e-15);
        // Danjon's projected coordinates are not exactly linear in time
        // (u varies as tan Δα), so his method is only approximate here.
        let r = conjunction::least_angular_separation_danjon(&[star], &body).unwrap();
        assert_close(r.n, n0, 2e-5);
        assert_close(r.sepr, 0.002, 2e-9);
    }
}

// Bodies that pass exactly through each other
#[test]
fn exact_occultation_has_zero_separation() {
    let star = (0.5, -0.3);
    let body: Vec<(f64, f64)> = (-2..=2)
        .map(|i| {
            let t = i as f64 - 0.4;
            (0.5 + 0.003 * t, -0.3 + 0.002 * t)
        })
        .collect();
    let r = conjunction::least_angular_separation(&body, &[star]).unwrap();
    assert_close(r.n, 0.4, 1e-9);
    assert!(r.sepr < 1e-11, "{}", r.sepr);
}

// A wide approach (about 20°), where Danjon's flat projection is poor but
// the exact method is still right. Body 1 moves along the equator, body 2
// is fixed at declination 20°: the least separation is 20° exactly.
#[test]
fn wide_approach() {
    let star = (2.0, 20_f64.to_radians());
    let body: Vec<(f64, f64)> = (-2..=2)
        .map(|i| (2.0 + 0.05 * (i as f64 - 1.3), 0.0))
        .collect();
    let r = conjunction::least_angular_separation(&body, &[star]).unwrap();
    assert_close(r.n, 1.3, 1e-10);
    assert_close(r.sepr.to_degrees(), 20.0, 1e-12);
}

#[test]
fn least_separation_none_cases() {
    let (m, s) = mercury_saturn();
    // Still approaching at the end of the table (September 12–14 only)
    assert!(conjunction::least_angular_separation(&m[0..3], &s[0..3]).is_some());
    let approaching: Vec<(f64, f64)> = (0..5).map(|i| (1.0 + 0.01 * i as f64, 0.0)).collect();
    assert!(conjunction::least_angular_separation(&approaching, &[(1.2, 0.0)]).is_none());
    assert!(conjunction::least_angular_separation_danjon(&approaching, &[(1.2, 0.0)]).is_none());
    // Invalid table lengths
    assert!(conjunction::least_angular_separation(&m[0..4], &s[0..4]).is_none());
    assert!(conjunction::least_angular_separation(&m, &s[0..3]).is_none());
    assert!(conjunction::least_angular_separation(&m[0..1], &s[0..1]).is_none());
    let empty: [(f64, f64); 0] = [];
    assert!(conjunction::least_angular_separation(&empty, &empty).is_none());
    assert!(conjunction::least_angular_separation_danjon(&m[0..2], &s[0..2]).is_none());
    // Non-finite input does not panic
    let mut bad = m;
    bad[2].1 = f64::NAN;
    let _ = conjunction::least_angular_separation(&bad, &s);
    let _ = conjunction::least_angular_separation_danjon(&bad, &s);
}

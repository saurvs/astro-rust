// Meeus, chapter 19: bodies in a straight line.

mod common;

use astro::angle::{self, Dms};
use astro::{coords, geometry};
use common::{assert_close, dms, hms};
use std::f64::consts::PI;

/// `(α, δ)` in radians
type Point = (f64, f64);

fn deg(a: f64, d: f64) -> (f64, f64) {
    (a.to_radians(), d.to_radians())
}

fn p(h: u32, m: u32, s: f64, neg: bool, d: u32, dm: u32, ds: f64) -> (f64, f64) {
    (hms(h, m, s).to_radians(), dms(neg, d, dm, ds).to_radians())
}

/// Castor, Pollux and Mars at 0h TD on 1994 September 29 to October 3, in
/// degrees as given by Meeus (example 19.a)
fn castor_pollux_mars() -> (Point, Point, [Point; 5]) {
    let castor = deg(113.56833, 31.89756);
    let pollux = deg(116.25042, 28.03681);
    let mars = [
        deg(118.98067, 21.68417),
        deg(119.59396, 21.58983),
        deg(120.20413, 21.49394),
        deg(120.81108, 21.39653),
        deg(121.41475, 21.29761),
    ];
    (castor, pollux, mars)
}

/// Stars of Orion's Belt, J2000 (Meeus, ch. 19)
fn orion_belt() -> [(f64, f64); 3] {
    [
        p(5, 32, 0.40, true, 0, 17, 56.9),  // δ Ori
        p(5, 36, 12.81, true, 1, 12, 7.0),  // ε Ori
        p(5, 40, 45.52, true, 1, 56, 33.3), // ζ Ori
    ]
}

// Meeus, example 19.a: the values of eq. 19.1 at the five instants
#[test]
fn example_19_a_criterion_values() {
    let (castor, pollux, mars) = castor_pollux_mars();
    let expected = [0.0019767, 0.0010851, 0.0001976, -0.0006855, -0.0015641];
    for (m, &e) in mars.iter().zip(expected.iter()) {
        assert_close(
            geometry::straight_line_criterion(&castor, &pollux, m),
            e,
            5e-8,
        );
    }
}

// Meeus, example 19.a: in a straight line on 1994 October 1.2233, at 5h TD
#[test]
fn example_19_a_time() {
    let (castor, pollux, mars) = castor_pollux_mars();
    let n = geometry::time_in_straight_line(&[castor], &[pollux], &mars).unwrap();
    assert_close(1.0 + n, 1.2233, 5e-5);
    assert_close((n * 24.0).round(), 5.0, 0.0);

    // The order of the bodies does not matter
    for n2 in [
        geometry::time_in_straight_line(&mars, &[castor], &[pollux]),
        geometry::time_in_straight_line(&[pollux], &mars, &[castor]),
        geometry::time_in_straight_line(&[pollux], &[castor], &mars),
    ] {
        assert_close(n2.unwrap(), n, 1e-12);
    }

    // The fixed stars given at every instant: same result
    let n5 = geometry::time_in_straight_line(&[castor; 5], &[pollux; 5], &mars).unwrap();
    assert_close(n5, n, 1e-12);

    // From the three central positions only
    let n3 = geometry::time_in_straight_line(&[castor], &[pollux], &mars[1..4]).unwrap();
    assert_close(n3, n, 1e-4);

    // From the sexagesimal values of the book
    let castor = p(7, 34, 16.40, false, 31, 53, 51.2);
    let pollux = p(7, 45, 0.10, false, 28, 2, 12.5);
    let mars = [
        p(7, 55, 55.36, false, 21, 41, 3.0),
        p(7, 58, 22.55, false, 21, 35, 23.4),
        p(8, 0, 48.99, false, 21, 29, 38.2),
        p(8, 3, 14.66, false, 21, 23, 47.5),
        p(8, 5, 39.54, false, 21, 17, 51.4),
    ];
    let n = geometry::time_in_straight_line(&[castor], &[pollux], &mars).unwrap();
    assert_close(1.0 + n, 1.2233, 5e-5);
}

#[test]
fn time_in_straight_line_with_equatorial_points() {
    let (castor, pollux, mars) = castor_pollux_mars();
    let eq = |q: (f64, f64)| coords::EqPoint { asc: q.0, dec: q.1 };
    let mars: Vec<coords::EqPoint> = mars.iter().map(|&q| eq(q)).collect();
    let n = geometry::time_in_straight_line(&[eq(castor)], &[eq(pollux)], &mars).unwrap();
    assert_close(1.0 + n, 1.2233, 5e-5);
}

#[test]
fn time_in_straight_line_invariants() {
    let (castor, pollux, mars) = castor_pollux_mars();
    let n = geometry::time_in_straight_line(&[castor], &[pollux], &mars).unwrap();

    // Rotating the right ascensions so that Mars crosses 0h in the table
    let rot = |q: (f64, f64)| (angle::limit_to_two_PI(q.0 - 120.0_f64.to_radians()), q.1);
    let mars_r = mars.map(rot);
    assert!(mars_r[0].0 > 6.0 && mars_r[4].0 < 0.1);
    let nr = geometry::time_in_straight_line(&[rot(castor)], &[rot(pollux)], &mars_r).unwrap();
    assert_close(nr, n, 1e-11);

    // Mirror image in the southern sky
    let south = |q: (f64, f64)| (q.0, -q.1);
    let ns = geometry::time_in_straight_line(&[south(castor)], &[south(pollux)], &mars.map(south))
        .unwrap();
    assert_close(ns, n, 1e-12);

    // Time reversal
    let rev = [mars[4], mars[3], mars[2], mars[1], mars[0]];
    let nt = geometry::time_in_straight_line(&[castor], &[pollux], &rev).unwrap();
    assert_close(nt, -n, 1e-12);
}

// A body moving uniformly in declination along a meridian crosses the
// celestial equator (the great circle through two equatorial stars) when
// its declination is zero — also when the stars are far apart and when the
// motion passes over the pole region.
#[test]
fn crossing_the_equator() {
    let s1 = deg(10.0, 0.0);
    let s2 = deg(100.0, 0.0);
    for &n0 in &[-1.7, -0.5, 0.0, 0.3, 1.9] {
        let body: Vec<(f64, f64)> = (-2..=2)
            .map(|i| deg(250.0, 2.0 * (i as f64 - n0)))
            .collect();
        let n = geometry::time_in_straight_line(&[s1], &[s2], &body).unwrap();
        assert_close(n, n0, 1e-12);
    }
    // A body 1° from the north pole moving uniformly in right ascension
    // crosses the great circle through the pole along the 0h and 12h
    // meridians when its right ascension is 0h.
    let a = deg(0.0, 80.0);
    let b = deg(180.0, 85.0);
    let body: Vec<(f64, f64)> = (-1..=1)
        .map(|i| deg(20.0 * (i as f64 - 0.25), 89.0))
        .collect();
    let n = geometry::time_in_straight_line(&[a], &[b], &body).unwrap();
    assert_close(n, 0.25, 1e-12);
    let at = deg(20.0 * (n - 0.25), 89.0);
    let omega = geometry::distance_from_great_circle(&a, &b, &at).unwrap();
    assert!(omega.abs() < 1e-14, "{omega}");
}

#[test]
fn time_in_straight_line_none_cases() {
    let (castor, pollux, mars) = castor_pollux_mars();
    // Not yet aligned on September 29 – October 1 (the criterion stays positive)
    assert!(geometry::time_in_straight_line(&[castor], &[pollux], &mars[0..3]).is_none());
    // All fixed
    assert!(geometry::time_in_straight_line(&[castor], &[pollux], &mars[0..1]).is_none());
    // Inconsistent or invalid lengths
    assert!(geometry::time_in_straight_line(&mars[0..3], &[pollux], &mars).is_none());
    assert!(geometry::time_in_straight_line(&[castor], &[pollux], &mars[0..4]).is_none());
    let empty: [(f64, f64); 0] = [];
    assert!(geometry::time_in_straight_line(&empty, &[pollux], &mars).is_none());
}

// Meeus, ch. 19: Orion's Belt. C1 + C2 = 172.4831° = 172°29′; ε Ori is
// 0.089876° = 324″ (south) from the great circle δ–ζ; Piessens' angle
// between the circles δε and εζ is 7°31′, the supplement of 172°29′.
#[test]
fn orion_belt_geometry() {
    let [delta, eps, zeta] = orion_belt();

    // The book's 172.4831 is the sum of C1 and C2 each rounded to four
    // decimals, so it is good to 1e-4 only (each part matches to 5e-5 in
    // tests/angle.rs).
    let c = geometry::vertex_angle(&delta, &eps, &zeta);
    assert_close(c.to_degrees(), 172.4831, 1e-4);
    let r = Dms::from_deg_rounded(c.to_degrees(), 0);
    assert_eq!((r.deg, r.min), (172, 28));
    assert_close(c.to_degrees() * 60.0, 172.0 * 60.0 + 29.0, 0.5);

    let omega = geometry::distance_from_great_circle(&delta, &zeta, &eps).unwrap();
    assert_close(omega.abs().to_degrees(), 0.089876, 5e-7);
    assert_close(omega.abs().to_degrees() * 3600.0, 324.0, 0.5);
    assert!(omega < 0.0, "ε Ori is south of the line");

    let psi = geometry::angle_between_great_circles(&delta, &eps, &eps, &zeta).unwrap();
    assert_close(psi.to_degrees() * 60.0, 7.0 * 60.0 + 31.0, 0.5);
    assert_close(psi + c, PI, 1e-12);

    // Cross-check with Meeus' own formula (X, Y, Z, A, B, C, m, n) ...
    let xyz = |q: (f64, f64)| [q.1.cos() * q.0.cos(), q.1.cos() * q.0.sin(), q.1.sin()];
    let ([x1, y1, z1], [x2, y2, z2]) = (xyz(delta), xyz(zeta));
    let (a, b, cc) = (y1 * z2 - z1 * y2, z1 * x2 - x1 * z2, x1 * y2 - y1 * x2);
    let m = eps.0.tan();
    let n = eps.1.tan() / eps.0.cos();
    let sin_w =
        (a + b * m + cc * n) / ((a * a + b * b + cc * cc).sqrt() * (1.0 + m * m + n * n).sqrt());
    assert_close(sin_w.abs().asin(), omega.abs(), 1e-15);
    // ... and with Piessens' sin ω for the angle between Oε and the plane
    // through δ and ζ (direction numbers a2, b2, c2 and l3, m3, n3).
    let ([a1, b1, c1], [a2, b2, c2], [a3, b3, c3]) = (xyz(delta), xyz(eps), xyz(zeta));
    let (l3, m3, n3) = (b1 * c3 - b3 * c1, c1 * a3 - c3 * a1, a1 * b3 - a3 * b1);
    let sin_w = (a2 * l3 + b2 * m3 + c2 * n3)
        / ((a2 * a2 + b2 * b2 + c2 * c2).sqrt() * (l3 * l3 + m3 * m3 + n3 * n3).sqrt());
    assert_close(sin_w.asin(), omega, 1e-15);
    // Piessens' cos ψ with l1, m1, n1 (δ, ε) and l2, m2, n2 (ε, ζ)
    let (l1, m1, n1) = (b1 * c2 - b2 * c1, c1 * a2 - c2 * a1, a1 * b2 - a2 * b1);
    let (l2, m2, n2) = (b2 * c3 - b3 * c2, c2 * a3 - c3 * a2, a2 * b3 - a3 * b2);
    let cos_psi = (l1 * l2 + m1 * m2 + n1 * n2)
        / ((l1 * l1 + m1 * m1 + n1 * n1).sqrt() * (l2 * l2 + m2 * m2 + n2 * n2).sqrt());
    assert_close(cos_psi, psi.cos(), 1e-12);

    // The criterion of eq. 19.1 is small but not zero, and has the sign
    // of the triple product's opposite.
    let crit = geometry::straight_line_criterion(&delta, &eps, &zeta);
    assert!(crit.abs() < 1e-4 && crit != 0.0);
}

// Meeus, ch. 19, exercise: the line from α to β UMa, extended, misses
// Polaris by ω = 1°55′.
#[test]
fn pointers_miss_polaris() {
    let alpha_uma = p(11, 3, 43.666, false, 61, 45, 3.22);
    let beta_uma = p(11, 1, 50.482, false, 56, 22, 56.65);
    let alpha_umi = p(2, 31, 48.704, false, 89, 15, 50.72);
    let omega = geometry::distance_from_great_circle(&alpha_uma, &beta_uma, &alpha_umi).unwrap();
    assert_close(omega.abs().to_degrees() * 60.0, 60.0 + 55.0, 0.5);
    // Same with the circle's points swapped, with the opposite sign
    let omega2 = geometry::distance_from_great_circle(&beta_uma, &alpha_uma, &alpha_umi).unwrap();
    assert_close(omega2, -omega, 1e-15);
}

#[test]
fn geometry_on_a_great_circle() {
    // Three points on the great circle inclined 30° to the equator with
    // its node at RA 40°
    let on_circle = |t: f64| {
        let (i, node) = (30_f64.to_radians(), 40_f64.to_radians());
        let x = t.cos();
        let y = t.sin() * i.cos();
        let z = t.sin() * i.sin();
        (node + y.atan2(x), z.asin())
    };
    let (a, b, c) = (on_circle(0.2), on_circle(0.5), on_circle(1.3));
    assert!(geometry::straight_line_criterion(&a, &b, &c).abs() < 1e-15);
    assert!(
        geometry::distance_from_great_circle(&a, &c, &b)
            .unwrap()
            .abs()
            < 1e-15
    );
    assert_close(geometry::vertex_angle(&a, &b, &c), PI, 1e-14);
    let psi = geometry::angle_between_great_circles(&a, &b, &b, &c).unwrap();
    assert!(psi < 1e-14);
    // Reversing the second circle gives π − ψ
    let psi = geometry::angle_between_great_circles(&a, &b, &c, &b).unwrap();
    assert_close(psi, PI, 1e-14);
    // The angle between the equator and this circle is the inclination
    let psi = geometry::angle_between_great_circles(&deg(0.0, 0.0), &deg(90.0, 0.0), &a, &b);
    assert_close(psi.unwrap().to_degrees(), 30.0, 1e-12);
}

#[test]
fn distance_from_great_circle_values() {
    // From the equator: the declination, positive to the north
    let e1 = deg(10.0, 0.0);
    let e2 = deg(50.0, 0.0);
    for &d in &[-89.0, -30.0, -0.001, 0.0, 12.5, 90.0] {
        let w = geometry::distance_from_great_circle(&e1, &e2, &deg(200.0, d)).unwrap();
        assert_close(w.to_degrees(), d, 1e-12);
    }
    // Degenerate great circles
    assert!(geometry::distance_from_great_circle(&e1, &e1, &e2).is_none());
    assert!(geometry::distance_from_great_circle(&deg(0.0, 90.0), &deg(0.0, 90.0), &e2).is_none());
    assert!(geometry::angle_between_great_circles(&e1, &e1, &e1, &e2).is_none());
    // A great circle through the pole (a meridian) works
    let w =
        geometry::distance_from_great_circle(&deg(30.0, 90.0), &deg(30.0, 0.0), &deg(31.0, 0.0));
    assert_close(w.unwrap().to_degrees().abs(), 1.0, 1e-12);
}

#[test]
fn vertex_angle_orientation() {
    let [delta, eps, zeta] = orion_belt();
    let c = geometry::vertex_angle(&delta, &eps, &zeta);
    // Taking the bodies in the opposite order gives the other side's angle
    let c_rev = geometry::vertex_angle(&zeta, &eps, &delta);
    assert_close(c + c_rev, 2.0 * PI, 1e-12);
    // A right angle at the vertex: north and east of it
    let v = deg(80.0, -20.0);
    let a = geometry::vertex_angle(&deg(80.0, -19.0), &v, &deg(81.0, -20.0));
    assert_close(a.to_degrees(), 90.0, 0.2);
    assert!((0.0..2.0 * PI).contains(&a));
}

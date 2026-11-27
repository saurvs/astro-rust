// Meeus, chapter 20: smallest circle containing three celestial bodies.

mod common;

use astro::geometry::{self, SmallestCircleKind};
use astro::{angle, interpol};
use common::{assert_close, dms, hms};

fn p(h: u32, m: u32, s: f64, neg: bool, d: u32, dm: u32, ds: f64) -> (f64, f64) {
    (hms(h, m, s).to_radians(), dms(neg, d, dm, ds).to_radians())
}

fn deg(a: f64, d: f64) -> (f64, f64) {
    (a.to_radians(), d.to_radians())
}

fn sepr(a: (f64, f64), b: (f64, f64)) -> f64 {
    angle::anglr_sepr(a.0, a.1, b.0, b.1)
}

// Meeus, example 20.a: Mercury, Jupiter, Saturn on 1981 September 11, 0h TD
#[test]
fn example_20_a() {
    let mercury = p(12, 41, 8.63, true, 5, 37, 54.2);
    let jupiter = p(12, 52, 5.21, true, 4, 22, 26.2);
    let saturn = p(12, 39, 28.11, true, 1, 50, 3.7);

    // Sides
    let b = sepr(mercury, jupiter).to_degrees();
    let c = sepr(mercury, saturn).to_degrees();
    let a = sepr(jupiter, saturn).to_degrees();
    assert_close(b, 3.00152, 5e-6);
    assert_close(c, 3.82028, 5e-6);
    assert_close(a, 4.04599, 5e-6);
    assert_close(b.hypot(c), 4.85836, 5e-6);

    let r = geometry::smallest_circle(&mercury, &jupiter, &saturn);
    assert_eq!(r.kind, SmallestCircleKind::TypeII);
    assert_close(r.diameter.to_degrees(), 4.26364, 5e-6);
    assert_close(r.diameter.to_degrees() * 60.0, 4.0 * 60.0 + 16.0, 0.5);
}

// Meeus, ch. 20, exercise: Venus, Mars, Jupiter on 1991 June 20, 0h TD:
// type I, Δ = 2°19′.
#[test]
fn venus_mars_jupiter_type_i() {
    let venus = p(9, 5, 41.44, false, 18, 30, 30.0);
    let mars = p(9, 9, 29.00, false, 17, 43, 56.7);
    let jupiter = p(8, 59, 47.14, false, 17, 49, 36.8);
    let r = geometry::smallest_circle(&venus, &mars, &jupiter);
    assert_eq!(r.kind, SmallestCircleKind::TypeI);
    assert_close(r.diameter.to_degrees() * 60.0, 2.0 * 60.0 + 19.0, 0.5);
    // The diameter is the longest side
    assert_eq!(r.diameter, sepr(mars, jupiter));
}

#[test]
fn order_of_the_bodies_does_not_matter() {
    let a = p(12, 41, 8.63, true, 5, 37, 54.2);
    let b = p(12, 52, 5.21, true, 4, 22, 26.2);
    let c = p(12, 39, 28.11, true, 1, 50, 3.7);
    let r = geometry::smallest_circle(&a, &b, &c);
    for q in [
        geometry::smallest_circle(&a, &c, &b),
        geometry::smallest_circle(&b, &a, &c),
        geometry::smallest_circle(&b, &c, &a),
        geometry::smallest_circle(&c, &a, &b),
        geometry::smallest_circle(&c, &b, &a),
    ] {
        assert_eq!(q.kind, r.kind);
        assert_close(q.diameter, r.diameter, 1e-15);
    }
}

#[test]
fn special_triangles() {
    let e = 0.01;
    // Equilateral (to first order): Δ = 2a/√3
    let a = deg(100.0, 0.0);
    let b = deg(100.0 + e, 0.0);
    let c = deg(100.0 + e / 2.0, e * 3_f64.sqrt() / 2.0);
    let r = geometry::smallest_circle(&a, &b, &c);
    assert_eq!(r.kind, SmallestCircleKind::TypeII);
    let side = sepr(a, b);
    assert_close(r.diameter, 2.0 * side / 3_f64.sqrt(), 1e-8);

    // Three bodies on one great circle: type I, the outer two on a diameter
    let r = geometry::smallest_circle(&deg(10.0, 0.0), &deg(12.0, 0.0), &deg(13.0, 0.0));
    assert_eq!(r.kind, SmallestCircleKind::TypeI);
    assert_close(r.diameter.to_degrees(), 3.0, 1e-12);

    // Two coincident bodies, and all three coincident
    let r = geometry::smallest_circle(&deg(10.0, 5.0), &deg(10.0, 5.0), &deg(11.0, 5.0));
    assert_eq!(r.kind, SmallestCircleKind::TypeI);
    assert_close(r.diameter, sepr(deg(10.0, 5.0), deg(11.0, 5.0)), 1e-15);
    let r = geometry::smallest_circle(&deg(1.0, 2.0), &deg(1.0, 2.0), &deg(1.0, 2.0));
    assert_eq!(r.kind, SmallestCircleKind::TypeI);
    assert_eq!(r.diameter, 0.0);
}

// Across the boundary between the two types (a right angle at the third
// body) the diameter is continuous.
#[test]
fn continuity_between_types() {
    // B and C 2° apart on the equator; A on the circle with diameter BC
    // has a right angle; moving A slightly in or out switches the type.
    let b = deg(50.0, 0.0);
    let c = deg(52.0, 0.0);
    let mut kinds = Vec::new();
    let mut prev: Option<f64> = None;
    for i in -20..=20 {
        let r_off = 1.0 + i as f64 * 1e-4;
        let t = 1.0_f64;
        let a = deg(51.0 + r_off * t.cos(), r_off * t.sin());
        let r = geometry::smallest_circle(&a, &b, &c);
        kinds.push(r.kind);
        if let Some(d) = prev {
            assert!((r.diameter - d).abs() < 1e-5, "{i}");
        }
        prev = Some(r.diameter);
    }
    assert_eq!(kinds[0], SmallestCircleKind::TypeI);
    assert_eq!(kinds[40], SmallestCircleKind::TypeII);
}

/// Ecliptic positions of Mercury, Jupiter, Saturn, 1981 September 7 to 14
/// at 0h TD (Meeus, example 20.b): (longitude, latitude) in degrees
const MERCURY: [(f64, f64); 8] = [
    (186.045, -0.560),
    (187.482, -0.696),
    (188.897, -0.833),
    (190.290, -0.971),
    (191.661, -1.109),
    (193.008, -1.246),
    (194.332, -1.384),
    (195.631, -1.521),
];
const JUPITER: [(f64, f64); 8] = [
    (192.866, 1.117),
    (193.069, 1.116),
    (193.272, 1.114),
    (193.476, 1.113),
    (193.681, 1.112),
    (193.886, 1.110),
    (194.092, 1.109),
    (194.299, 1.108),
];
const SATURN: [(f64, f64); 8] = [
    (189.324, 2.226),
    (189.439, 2.225),
    (189.555, 2.224),
    (189.671, 2.223),
    (189.788, 2.222),
    (189.906, 2.221),
    (190.023, 2.220),
    (190.142, 2.219),
];

/// Position interpolated from the 20.b table at September `t`, with five
/// values centred on the nearest day for which five are available.
fn at(table: &[(f64, f64); 8], t: f64) -> (f64, f64) {
    let k = (t.round() as i64).clamp(9, 12);
    let i = (k - 7) as usize;
    let n = t - k as f64;
    let f = |g: fn(&(f64, f64)) -> f64| {
        let v: Vec<f64> = table[i - 2..=i + 2].iter().map(g).collect();
        interpol::five_values(v[0], v[1], v[2], v[3], v[4], n)
    };
    deg(f(|q| q.0), f(|q| q.1))
}

/// Samples the smallest circle at steps of 0.001 day from `t0` to `t1`;
/// returns the instants where the type changes and the least diameter.
fn scan(
    t0: f64,
    t1: f64,
    circle: impl Fn(f64) -> geometry::SmallestCircle,
) -> (Vec<(f64, SmallestCircleKind)>, (f64, f64)) {
    let mut changes = Vec::new();
    let mut prev = None;
    let mut least = (f64::INFINITY, 0.0);
    let steps = ((t1 - t0) / 0.001).round() as i64;
    for s in 0..=steps {
        let t = t0 + s as f64 * 0.001;
        let c = circle(t);
        if prev.is_some_and(|k| k != c.kind) {
            changes.push((t, c.kind));
        }
        prev = Some(c.kind);
        if c.diameter < least.0 {
            least = (c.diameter, t);
        }
    }
    (changes, least)
}

// Meeus, example 20.b: type I from September 7.00 to 8.81 (Δ from 7°01′ to
// about 5°00′), type II from 8.81 to 12.19 with a least Δ of 4°14′ on
// September 10.53, then type I again.
#[test]
fn example_20_b_grouping_over_time() {
    let circle =
        |t: f64| geometry::smallest_circle(&at(&MERCURY, t), &at(&JUPITER, t), &at(&SATURN, t));

    let start = circle(7.0);
    assert_eq!(start.kind, SmallestCircleKind::TypeI);
    assert_close(start.diameter.to_degrees() * 60.0, 7.0 * 60.0 + 1.0, 0.5);

    let (changes, (least, t_least)) = scan(7.0, 14.0, circle);
    assert_eq!(changes.len(), 2, "{changes:?}");
    assert_eq!(changes[0].1, SmallestCircleKind::TypeII);
    assert_close(changes[0].0, 8.81, 0.005);
    assert_eq!(changes[1].1, SmallestCircleKind::TypeI);
    assert_close(changes[1].0, 12.19, 0.005);
    assert_close(circle(8.81).diameter.to_degrees() * 60.0, 5.0 * 60.0, 1.0);

    assert_close(least.to_degrees() * 60.0, 4.0 * 60.0 + 14.0, 0.5);
    assert_close(t_least, 10.53, 0.005);
}

// Meeus, example 20.c (fictitious planets with uniform motion in
// longitude): type II only during March 15.87–15.91 and 17.93–18.05; the
// least Δ is 1°55′ on March 17.94.
#[test]
fn example_20_c_discontinuities() {
    let planet = |l0: f64, b: f64, rate: f64, t: f64| deg(l0 + rate * (t - 12.0), b);
    let circle = |t: f64| {
        geometry::smallest_circle(
            &planet(214.23, 0.29, 0.11, t),
            &planet(211.79, 0.48, 0.20, t),
            &planet(208.41, 0.75, 1.08, t),
        )
    };
    let (changes, (least, t_least)) = scan(12.0, 22.0, circle);
    let t: Vec<f64> = changes.iter().map(|c| c.0).collect();
    assert_eq!(t.len(), 4, "{changes:?}");
    assert_eq!(changes[0].1, SmallestCircleKind::TypeII);
    for (got, want) in t.iter().zip([15.87, 15.91, 17.93, 18.05]) {
        assert_close(*got, want, 0.005);
    }
    assert_close(least.to_degrees() * 60.0, 60.0 + 55.0, 0.5);
    assert_close(t_least, 17.94, 0.005);
    // About 5.8° on March 12 (figure)
    assert_close(circle(12.0).diameter.to_degrees(), 5.8, 0.1);
}

#[test]
fn smallest_circle_wraps_and_mirrors() {
    let a = p(23, 58, 0.0, true, 1, 0, 0.0);
    let b = p(0, 3, 0.0, true, 2, 30, 0.0);
    let c = p(0, 0, 30.0, false, 0, 40, 0.0);
    let r = geometry::smallest_circle(&a, &b, &c);
    // Shift the right ascensions by 12h: same result
    let sh = |q: (f64, f64)| (q.0 + std::f64::consts::PI, q.1);
    let q = geometry::smallest_circle(&sh(a), &sh(b), &sh(c));
    assert_eq!(q.kind, r.kind);
    assert_close(q.diameter, r.diameter, 1e-14);
    // Mirror in declination: same result
    let m = |q: (f64, f64)| (q.0, -q.1);
    let q = geometry::smallest_circle(&m(a), &m(b), &m(c));
    assert_eq!(q.kind, r.kind);
    assert_close(q.diameter, r.diameter, 1e-14);
    assert!(r.diameter.to_degrees() < 4.0);
}

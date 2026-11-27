//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 46: positions of the
//! satellites of Saturn.

#![allow(clippy::approx_constant)] // 0.318 is a book value, not 1/π

mod common;

use astro::planet;
use astro::planet::saturn::moon::{self, Moon};
use common::{assert_close, assert_close_deg};

/// Example 46.a: 1999 September 18 at 0h UT, ΔT = 64 s
const JDE_46: f64 = 2451439.50074;

// Meeus, example 46.a: geometry of Saturn
#[test]
fn saturn_geometry_example_46a() {
    let p = moon::positions(JDE_46);
    assert_close(p.light_time, 0.04942, 1e-5);
    assert_close(p.earth_dist, 8.557599, 2e-6);
    // The book's positions are in the FK5 frame, 0.09" (2.5e-5°) away from
    // the dynamical frame of VSOP87 used here.
    assert_close_deg(p.saturn_long.to_degrees(), 46.865071, 5e-5);
    assert_close(p.saturn_lat.to_degrees(), -2.539334, 5e-6);
    // Reducing the book's own λ0 = 46.865071° to B1950.0 with eqs. 21.5 and
    // 21.7, as the text prescribes, gives 46.170517°, not the 46.170287°
    // printed; the 0.83" difference moves X and Y by less than 1.5e-4.
    assert_close_deg(p.saturn_long_b1950.to_degrees(), 46.170287, 3e-4);
    assert_close(p.saturn_lat_b1950.to_degrees(), -2.544441, 5e-6);
}

// Meeus, example 46.a: final coordinates of the eight satellites
#[test]
fn satellites_example_46a() {
    let x = [
        3.102, 3.823, 4.027, -5.365, -1.122, 14.568, -18.001, -48.759,
    ];
    let y = [-0.204, 0.318, -1.061, -1.148, -3.123, 4.738, -5.328, 4.136];
    // The book prints three decimals. For Iapetus the book's run used the
    // opposite sign for the term 0.00650 sin(u5 + ψ) of λ'; the code keeps
    // the printed (and geometrically consistent) "+", with which X and Y
    // are −48.7604 and 4.1372, as an independent Go implementation
    // (soniakeys/meeus) also finds. See the comment in `iapetus()` in
    // src/planet/saturn/moon.rs.
    let tol = [5e-4, 5e-4, 5e-4, 5e-4, 5e-4, 5e-4, 5e-4, 2e-3];

    let p = moon::positions(JDE_46);
    for m in Moon::ALL {
        let k = m.number() as usize - 1;
        let s = p.moon(m);
        assert_close(s.x, x[k], tol[k]);
        assert_close(s.y, y[k], tol[k]);
        // the original entry point gives the same coordinates
        assert_eq!(moon::apprnt_rect_coords(JDE_46, &m), (s.x, s.y, s.z));
    }
    // pin down the Iapetus values
    let iapetus = p.moon(Moon::Iapetus);
    assert_close(iapetus.x, -48.7604, 1e-4);
    assert_close(iapetus.y, 4.1372, 1e-4);
}

// Meeus, example 46.a: λ, γ, Ω and r, through the public API. `positions`
// uses the unrounded light-time 0.0494243 d, while the book worked at
// JDE − 0.04942 = 2451439.45132, 0.37 s later; that moves λ of the fast inner
// satellites by up to 0.0017° (Mimas), hence the λ tolerance of 0.002°.
// Titan's λ and r and Iapetus' λ and r carry the discrepancies described in
// the source. The unit tests in the module check every value at the book's
// exact instant.
#[test]
fn orbital_quantities_example_46a() {
    let expected = [
        (320.0015, 1.5630, 47.7110, 3.1224),
        (300.4638, 0.0262, 123.2135, 3.9257),
        (347.9653, 1.0976, 51.0615, 4.8810),
        (102.6463, 0.0139, 128.2982, 6.2561),
        (50.9947, 0.3359, 118.2589, 8.7054),
        (270.5289, 0.3701, 8.4467, 19.9254),
        (91.9269, 1.0461, 21.5992, 24.1160),
        (103.1318, 15.5062, 22.3756, 58.9396),
    ];
    let lambda_tol = [2e-3, 2e-3, 2e-3, 2e-3, 2e-3, 4.5e-3, 2e-3, 3e-3];
    let r_tol = [5e-5, 5e-5, 5e-5, 5e-5, 5e-5, 1.5e-3, 5e-5, 1e-4];
    let p = moon::positions(JDE_46);
    for (k, (lambda, gamma, node, r)) in expected.into_iter().enumerate() {
        let s = p.moons[k];
        assert_close_deg(s.longitude.to_degrees(), lambda, lambda_tol[k]);
        assert_close(s.inclination.to_degrees(), gamma, 5e-5);
        assert_close_deg(s.node.to_degrees(), node, 5e-5);
        assert_close(s.radius, r, r_tol[k]);
    }
}

// Bug fixed: the light-time was hard-coded to the 0.04942 day of example
// 46.a, so Mimas (382°/day) was misplaced by up to 7° when Saturn was near
// conjunction. It is now iterated from the Earth–Saturn distance.
#[test]
fn light_time_follows_the_distance() {
    // Saturn near opposition (1999 Nov 6) and near conjunction (2000 May 10)
    for (jde, near) in [(2451488.5, true), (2451674.5, false)] {
        let p = moon::positions(jde);

        // independent light-time iteration with the planetary theory
        let (l0, b0, r0) = planet::heliocent_coords(&planet::Planet::Earth, jde);
        let mut tau = 0.0;
        let mut dist = 0.0;
        for _ in 0..5 {
            let (l, b, r) = planet::heliocent_coords(&planet::Planet::Saturn, jde - tau);
            let x = r * b.cos() * l.cos() - r0 * b0.cos() * l0.cos();
            let y = r * b.cos() * l.sin() - r0 * b0.cos() * l0.sin();
            let z = r * b.sin() - r0 * b0.sin();
            dist = (x * x + y * y + z * z).sqrt();
            tau = 0.0057755183 * dist;
        }
        assert_close(p.earth_dist, dist, 1e-7);
        assert_close(p.light_time, tau, 1e-9);
        // 8.21 AU (τ = 0.0474 d) and 10.17 AU (τ = 0.0587 d): the old fixed
        // value was wrong by 0.002 and 0.009 day, i.e. by 0.8° and 3.6° in
        // the longitude of Mimas.
        if near {
            assert!(p.earth_dist < 8.3);
        } else {
            assert!(p.earth_dist > 10.1);
        }
        assert!((p.light_time - 0.04942).abs() > 0.0019);
    }
}

#[test]
fn coordinates_are_consistent() {
    let mut jde = 2433282.5;
    while jde < 2469807.5 {
        let p = moon::positions(jde);
        for (k, s) in p.moons.iter().enumerate() {
            // the corrections change the length by a fraction of a percent
            let n = (s.x * s.x + s.y * s.y + s.z * s.z).sqrt();
            assert!((n / s.radius - 1.0).abs() < 0.005, "{k} {n} {}", s.radius);
            for a in [s.longitude, s.node] {
                assert!((0.0..std::f64::consts::TAU).contains(&a));
            }
            // Iapetus' γ ranges over 14.7°–16.2° in 1950–2050
            assert!(s.inclination >= 0.0 && s.inclination < 17_f64.to_radians());
        }
        // the satellites never overtake each other in distance
        for k in 1..8 {
            assert!(p.moons[k].radius > p.moons[k - 1].radius);
        }
        for a in [p.saturn_long, p.saturn_long_b1950] {
            assert!((0.0..std::f64::consts::TAU).contains(&a));
        }
        jde += 97.3;
    }
}

// Enceladus and Dione move almost exactly in Saturn's equatorial (ring)
// plane, which we see at the Saturnicentric latitude B of the Earth. A point
// of that plane has Y = Z tan B, so Y/Z of both satellites must give the B
// of chapter 45.
#[test]
fn equatorial_satellites_lie_in_the_ring_plane() {
    for jde in [2448972.50068, JDE_46, 2453000.5, 2455000.5] {
        let p = moon::positions(jde);
        let ring = planet::saturn::ring::elements(jde, 0.0, 23.4393_f64.to_radians());
        for m in [Moon::Enceladus, Moon::Dione] {
            let s = p.moon(m);
            if s.z.abs() < 0.5 * s.radius {
                continue; // near elongation Y/Z is ill-conditioned
            }
            let b = (s.y / s.z).atan().to_degrees();
            assert_close(b, ring.B.to_degrees(), 0.05);
        }
    }
}

#[test]
fn positions_are_continuous() {
    let h = 0.01;
    let start = 2451430.0;
    let mut prev2 = moon::positions(start - h);
    let mut prev = moon::positions(start);
    for n in 1..=500 {
        let cur = moon::positions(start + n as f64 * h);
        for k in 0..8 {
            for (a, b, c) in [
                (prev2.moons[k].x, prev.moons[k].x, cur.moons[k].x),
                (prev2.moons[k].y, prev.moons[k].y, cur.moons[k].y),
                (prev2.moons[k].z, prev.moons[k].z, cur.moons[k].z),
            ] {
                // |second difference| ≤ r ω² h², 0.014 for Mimas
                assert!(
                    (a - 2.0 * b + c).abs() < 0.02,
                    "step {n}, satellite {}",
                    k + 1
                );
            }
        }
        prev2 = prev;
        prev = cur;
    }
}

#[test]
fn far_dates_are_finite() {
    // 500 BC, 1000 AD, 3000 AD
    for jde in [1538432.5, 2086302.5, 2816787.5] {
        let p = moon::positions(jde);
        for s in p.moons.iter() {
            assert!(s.x.is_finite() && s.y.is_finite() && s.z.is_finite());
        }
    }
}

// Bad input must propagate as NaN, not panic or come out finite. Before
// the fix the inclination of satellites V-VIII came out as exactly 90°.
#[test]
fn nan_and_infinite_dates() {
    for jde in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let p = moon::positions(jde);
        assert!(p.light_time.is_nan() && p.earth_dist.is_nan());
        for m in Moon::ALL {
            let s = p.moon(m);
            assert!(s.x.is_nan() && s.y.is_nan() && s.z.is_nan(), "{m:?}");
            assert!(s.longitude.is_nan() && s.node.is_nan(), "{m:?}");
            // Tethys moves on a circle of fixed radius
            if m != Moon::Tethys {
                assert!(s.radius.is_nan(), "{m:?}");
            }
            // satellites I-IV have constant inclinations in this theory
            if m.number() >= 5 {
                assert!(s.inclination.is_nan(), "{m:?}: {}", s.inclination);
            }
        }
    }
}

#[test]
fn moon_numbers() {
    let n: Vec<u8> = Moon::ALL.iter().map(|m| m.number()).collect();
    assert_eq!(n, (1..=8).collect::<Vec<u8>>());
}

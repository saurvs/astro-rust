// Meeus, Astronomical Algorithms, 2nd ed., chapter 34: parabolic motion

#![allow(non_snake_case)]

mod common;

use astro::orbit::{self, Elements, Motion, elliptic, near_parabolic, parabolic};
use astro::{consts, time};
use common::assert_close;

fn jd(year: i16, month: time::Month, day: f64) -> f64 {
    time::julian_day(&time::Date {
        year,
        month,
        decimal_day: day,
        cal_type: time::CalType::Gregorian,
    })
}

// Meeus, example 34.a: comet Stonehouse (C/1998 H1)
fn stonehouse() -> (f64, f64, f64) {
    let T = jd(1998, time::Month::Apr, 14.4358);
    let t = jd(1998, time::Month::Aug, 5.0);
    (t, T, 1.487469)
}

#[test]
fn constants() {
    // 0.03649116245 = 3k/√2 (eq. 34.1); 54.807791 = 2√2/(3k) (eq. 34.5)
    let k = consts::GAUSS_GRAV;
    assert_close(3.0 * k / 2_f64.sqrt(), 0.03649116245, 1e-11);
    assert_close(2.0 * 2_f64.sqrt() / (3.0 * k), 54.807791, 1e-6);
    let q: f64 = 1.7;
    assert_close(
        2.0 / parabolic::barker_w(10.0, 0.0, q),
        54.807791 * q * q.sqrt() / 10.0,
        1e-6,
    );
}

#[test]
fn example_34a_w() {
    let (t, T, q) = stonehouse();
    assert_close(t - T, 112.5642, 1e-9);
    assert_close(parabolic::barker_w(t, T, q), 2.264206862, 1e-9);
}

#[test]
fn example_34a_iteration() {
    // Meeus, example 34.a: successive values of eq. 34.4, from s = 0
    let (t, T, q) = stonehouse();
    let W = parabolic::barker_w(t, T, q);
    let expected = [0.7547356, 0.6634364, 0.6592441, 0.6592360, 0.6592360];
    let mut s = 0.0;
    for x in expected {
        s = (2.0 * s * s * s + W) / (3.0 * (s * s + 1.0));
        assert_close(s, x, 1e-7);
    }
    // the library's iteration
    let s = parabolic::solve_barker_iter(W, 1e-7, 5).unwrap();
    assert_close(s, 0.6592360, 1e-7);
    let s = parabolic::solve_barker_iter(W, 1e-15, 20).unwrap();
    assert_close(s, parabolic::solve_barker(W), 1e-15);
    // not converged after 3 steps
    assert_eq!(parabolic::solve_barker_iter(W, 1e-7, 3), None);
    assert_eq!(parabolic::solve_barker_iter(f64::NAN, 1e-7, 30), None);
}

#[test]
fn example_34a_direct_formulae() {
    // Meeus, example 34.a with eq. 34.6: G = 1.132103431, Y = 1.382541577
    let (t, T, q) = stonehouse();
    let W = parabolic::barker_w(t, T, q);
    let G = W / 2.0;
    assert_close(G, 1.132103431, 1e-9);
    let Y = (G + (G * G + 1.0).sqrt()).cbrt();
    assert_close(Y, 1.382541577, 1e-9);
    assert_close(Y - 1.0 / Y, 0.6592360, 1e-7);

    assert_close(parabolic::solve_barker_cubic(W), 0.6592360, 1e-7);
    assert_close(parabolic::solve_barker_bauschinger(W), 0.6592360, 1e-7);
    assert_close(parabolic::solve_barker(W), 0.6592360, 1e-7);
}

#[test]
fn example_34a_true_anomaly_and_radius_vector() {
    // Meeus, example 34.a: v = +66.78862°, r = 2.133911 AU
    let (t, T, q) = stonehouse();
    let (v, r) = parabolic::true_anom_and_rad_vec(t, T, q);
    assert_close(v.to_degrees(), 66.78862, 1e-5);
    assert_close(r, 2.133911, 1e-6);
    // eq. 30.4 with e = 1 gives the same radius vector
    assert_close(orbit::rad_vec_frm_true_anom(v, q, 1.0), r, 1e-14);
}

#[test]
fn at_perihelion() {
    let (v, r) = parabolic::true_anom_and_rad_vec(2451000.0, 2451000.0, 0.7);
    assert_eq!(v, 0.0);
    assert_eq!(r, 0.7);
    assert_eq!(parabolic::solve_barker(0.0), 0.0);
    assert_eq!(parabolic::solve_barker_bauschinger(0.0), 0.0);
    assert_eq!(parabolic::solve_barker_cubic(0.0), 0.0);
    assert_eq!(parabolic::solve_barker_iter(0.0, 0.0, 2), Some(0.0));
}

#[test]
fn barker_residual_over_wide_range() {
    for k in -40..=40 {
        for sign in [-1.0, 1.0] {
            let W = sign * 10_f64.powf(k as f64 * 0.3);
            let s = parabolic::solve_barker(W);
            assert_eq!(s.signum(), W.signum());
            // relative residual of s³ + 3s = W
            let res = (s * s * s + 3.0 * s - W) / W;
            assert!(res.abs() < 1e-14, "W = {W}: residual {res}");
            // all methods agree
            assert_close(
                parabolic::solve_barker_cubic(W),
                s,
                1e-12 * s.abs().max(1e-3),
            );
            assert_close(
                parabolic::solve_barker_bauschinger(W),
                s,
                1e-9 * s.abs().max(1e-3),
            );
            let it = parabolic::solve_barker_iter(W, 1e-15 * s.abs(), 400).unwrap();
            assert_close(it, s, 1e-13 * s.abs());
        }
    }
}

#[test]
fn symmetric_before_and_after_perihelion() {
    // Bug fix: eq. 34.6 was evaluated directly with a negative W, where
    // G + √(G² + 1) suffers from cancellation. Before, r one million days
    // before the perihelion differed from r one million days after it by
    // 3.7e-8 (relative).
    let (T, q) = (2451000.0, 1.3);
    for dt in [1e-3, 1.0, 112.5642, 1e3, 1e5, 1e6, 1e7, 1e8] {
        let (v1, r1) = parabolic::true_anom_and_rad_vec(T + dt, T, q);
        let (v2, r2) = parabolic::true_anom_and_rad_vec(T - dt, T, q);
        assert_eq!(v2, -v1);
        assert_close(r2, r1, 4.0 * f64::EPSILON * r1);
        assert!(v1 > 0.0 && v1 < std::f64::consts::PI);
        assert!(r1 > q);
    }
    for W in [1e-9, 0.5, 2.0, 1e4, 1e9] {
        assert_eq!(
            parabolic::solve_barker_cubic(-W),
            -parabolic::solve_barker_cubic(W)
        );
        assert_eq!(
            parabolic::solve_barker_bauschinger(-W),
            -parabolic::solve_barker_bauschinger(W)
        );
    }
}

#[test]
fn true_anomaly_increases_with_time() {
    let (T, q) = (0.0, 0.5);
    let mut prev = -std::f64::consts::PI;
    for i in -1000..=1000 {
        let (v, _) = parabolic::true_anom_and_rad_vec(i as f64 * 3.7, T, q);
        assert!(v > prev);
        prev = v;
    }
}

#[test]
fn limit_of_elliptic_and_near_parabolic_orbits() {
    // An ellipse with the same q and e very close to 1 moves almost like
    // the parabola near the perihelion
    let q = 1.487469;
    let e = 1.0 - 1e-9;
    let a = elliptic::semimaj_axis(q, e);
    let ell = Motion::elliptic_frm_perih_time(a, e, 0.0);
    for dt in [-300.0, -50.0, 1.0, 112.5642, 400.0] {
        let (v, r) = parabolic::true_anom_and_rad_vec(dt, 0.0, q);
        let (ve, re) = ell.true_anom_and_rad_vec(dt).unwrap();
        assert_close(ve, v, 1e-7);
        assert_close(re, r, 1e-7 * r);
        // Landgraf's method with e = 1 is exactly the parabolic case
        let (vn, rn) = near_parabolic::true_anom_and_rad_vec(dt, 0.0, 1.0, q, 1e-12).unwrap();
        common::assert_close_rad(vn, v, 1e-12);
        assert_close(rn, r, 1e-12 * r);
    }
}

#[test]
fn geocentric_position_of_a_parabolic_comet() {
    // The orbit machinery of chapter 33 applies unchanged (Meeus, ch. 34)
    let (t, T, q) = stonehouse();
    let el = Elements {
        motion: Motion::Parabolic {
            perih_dist: q,
            perih_time: T,
        },
        incl: 104.6948_f64.to_radians(),
        arg_perih: 344.9632_f64.to_radians(),
        long_asc_node: 222.1103_f64.to_radians(),
    };
    let sun = (0.6, -0.75, -0.33);
    let p = orbit::geocent_eq_pos(&el, t, sun, orbit::OBLQ_J2000).unwrap();
    let (v, r) = parabolic::true_anom_and_rad_vec(t - p.light_time, T, q);
    let h = orbit::eq_orbit_consts(el.long_asc_node, el.incl, orbit::OBLQ_J2000).heliocent_rect(
        el.arg_perih,
        v,
        r,
    );
    let (asc, dec, dist) = orbit::geocent_eq_frm_rect(h, sun);
    assert_close(p.asc, asc, 1e-12);
    assert_close(p.dec, dec, 1e-12);
    assert_close(p.earth_dist, dist, 1e-12);
    assert_close(p.sun_dist, r, 1e-12);
    assert_close(p.light_time, orbit::light_time(dist), 1e-12);
}

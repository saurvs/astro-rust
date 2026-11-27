// Meeus, Astronomical Algorithms, 2nd ed., chapter 33: elliptic motion
// (bodies given by orbital elements), and the magnitudes of minor planets
// and comets

#![allow(non_snake_case)]

mod common;

use astro::orbit::{self, Elements, Motion, elliptic};
use astro::{asteroid, time};
use common::{assert_close, assert_close_deg};
use std::f64::consts::PI;

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

// Meeus, example 33.b: periodic comet Encke, elements for J2000.0
const ENCKE_A: f64 = 2.2091404;
const ENCKE_E: f64 = 0.8502196;
const ENCKE_I: f64 = 11.94524;
const ENCKE_NODE: f64 = 334.75006;
const ENCKE_W: f64 = 186.23352;
// The Sun's geocentric rectangular equatorial coordinates, J2000.0, on
// 1990 Oct. 6.0 TD
const SUN: (f64, f64, f64) = (-0.9756732, -0.2003254, -0.0868566);

fn encke_T() -> f64 {
    jd(1990, time::Month::Oct, 28.54502)
}

fn encke() -> Elements {
    Elements {
        motion: Motion::elliptic_frm_perih_time(ENCKE_A, ENCKE_E, encke_T()),
        incl: rad(ENCKE_I),
        arg_perih: rad(ENCKE_W),
        long_asc_node: rad(ENCKE_NODE),
    }
}

#[test]
fn example_33b_orbit_constants() {
    // Meeus, example 33.b, eqs. 33.7 and 33.8
    let c = orbit::eq_orbit_consts(rad(ENCKE_NODE), rad(ENCKE_I), orbit::OBLQ_J2000);
    assert_close(c.A.to_degrees(), 65.230615, 1e-6);
    assert_close(c.B.to_degrees(), 331.787680, 1e-6);
    assert_close(c.C.to_degrees(), 342.613052, 1e-6);
    assert_close(c.a, 0.99609485, 1e-8);
    assert_close(c.b, 0.82787174, 1e-8);
    assert_close(c.c, 0.56782342, 1e-8);
    // F² + G² + H² = 1 and P² + Q² + R² = 1, hence a² + b² + c² = 2
    assert_close(c.a * c.a + c.b * c.b + c.c * c.c, 2.0, 1e-15);
    // sin ε and cos ε for J2000.0
    assert_close(orbit::OBLQ_J2000.sin(), 0.397777156, 1e-9);
    assert_close(orbit::OBLQ_J2000.cos(), 0.917482062, 1e-9);
}

#[test]
fn example_33b_intermediate_values() {
    // Meeus, example 33.b
    let T = encke_T();
    let t = jd(1990, time::Month::Oct, 6.0);
    assert_eq!(t, 2448170.5);
    assert_close(t - T, -22.54502, 1e-9);

    let n = elliptic::mn_motion(ENCKE_A);
    assert_close(n.to_degrees(), 0.300171252, 1e-9);
    let M = n * (t - T);
    assert_close(M.to_degrees(), -6.767367, 1e-6);

    let e = ENCKE_E;
    let E = elliptic::ecc_anom_checked(M, e).unwrap();
    assert_close(E.to_degrees(), -34.026714, 1e-6);
    let v = elliptic::true_anom(E, e);
    assert_close(v.to_degrees(), -94.163310, 1e-6);
    let r = elliptic::rad_vec_frm_ecc_anom(E, ENCKE_A, e);
    assert_close(r, 0.6524867, 1e-7);

    let consts = orbit::eq_orbit_consts(rad(ENCKE_NODE), rad(ENCKE_I), orbit::OBLQ_J2000);
    let (x, y, z) = consts.heliocent_rect(rad(ENCKE_W), v, r);
    assert_close(x, 0.2508066, 1e-7);
    assert_close(y, 0.4849175, 1e-7);
    assert_close(z, 0.3573373, 1e-7);

    let (_, _, dist) = orbit::geocent_eq_frm_rect((x, y, z), SUN);
    assert_close(dist, 0.8243689, 1e-7);
    let tau = orbit::light_time(dist);
    assert_close(tau, 0.00476, 1e-5);

    // again for t - τ = 1990 Oct. 5.99524 (the book's second-pass values
    // correspond to the unrounded first-pass τ = 0.0047611)
    let t2 = t - tau;
    assert_close(t2, jd(1990, time::Month::Oct, 5.99524), 1e-5);
    let M = n * (t2 - T);
    assert_close(M.to_degrees(), -6.768796, 1e-6);
    let E = elliptic::ecc_anom_checked(M, e).unwrap();
    assert_close(E.to_degrees(), -34.031552, 1e-6);
    let v = elliptic::true_anom(E, e);
    assert_close(v.to_degrees(), -94.171933, 1e-6);
    let r = elliptic::rad_vec_frm_ecc_anom(E, ENCKE_A, e);
    assert_close(r, 0.6525755, 1e-7);
    let (x, y, z) = consts.heliocent_rect(rad(ENCKE_W), v, r);
    assert_close(x, 0.2509310, 1e-7);
    assert_close(y, 0.4849477, 1e-7);
    assert_close(z, 0.3573712, 1e-7);
    assert_close(SUN.0 + x, -0.7247422, 1e-7);
    assert_close(SUN.1 + y, 0.2846223, 1e-7);
    assert_close(SUN.2 + z, 0.2705146, 1e-7);

    let (asc, dec, dist) = orbit::geocent_eq_frm_rect((x, y, z), SUN);
    assert_close(dist, 0.8242811, 1e-7);
    // α and δ follow from the rounded ξ, η, ζ only to a few 1e-6°
    assert_close(asc.to_degrees(), 158.558965, 3e-6);
    assert_close(dec.to_degrees(), 19.158496, 3e-6);
}

#[test]
fn example_33b_geocentric_position() {
    // Meeus, example 33.b: α2000 = 10h34m14.2s, δ2000 = +19°09′31″, ψ = 40.51°
    let t = jd(1990, time::Month::Oct, 6.0);
    let p = orbit::geocent_eq_pos(&encke(), t, SUN, orbit::OBLQ_J2000).unwrap();

    // The Sun's coordinates are given to 1e-7 AU, i.e. about 7e-6° as seen
    // from the comet's distance, and the book stops the light-time
    // iteration after one pass: compare to 1e-5° (0.04″).
    assert_close(p.asc.to_degrees(), 158.558965, 1e-5);
    assert_close(p.asc.to_degrees(), common::hms(10, 34, 14.2), 0.05 / 240.0);
    assert_close(p.dec.to_degrees(), 19.158496, 1e-5);
    assert_close(
        p.dec.to_degrees(),
        common::dms(false, 19, 9, 31.0),
        0.5 / 3600.0,
    );
    assert_close(p.earth_dist, 0.8242811, 1e-7);
    assert_close(p.sun_dist, 0.6525755, 1e-7);
    assert_close(p.light_time, 0.00476, 1e-5);
    assert_close(p.light_time, orbit::light_time(p.earth_dist), 1e-11);
    assert_close(p.elong.to_degrees(), 40.51, 0.005);

    // eqs. 33.11 and 33.12 give the same angles from the distances
    let R = (SUN.0 * SUN.0 + SUN.1 * SUN.1 + SUN.2 * SUN.2).sqrt();
    let elong = orbit::elong_frm_dist(R, p.earth_dist, p.sun_dist);
    let phase = orbit::phase_angl_frm_dist(p.sun_dist, p.earth_dist, R);
    assert_close(elong, p.elong, 1e-12);
    assert_close(phase, p.phase_angl, 1e-12);
    // the three angles of the Sun–Earth–comet triangle add up to π
    let sun_angle = orbit::phase_angl_frm_dist(R, p.sun_dist, p.earth_dist);
    assert_close(p.elong + p.phase_angl + sun_angle, PI, 1e-12);
}

#[test]
fn gauss_vectors_agree_with_orbit_constants() {
    // Meeus, chapter 33: the P, Q vectors are an alternative to eq. 33.9
    let cases = [
        (rad(ENCKE_W), rad(ENCKE_NODE), rad(ENCKE_I)),
        (0.3, 4.0, 2.8), // retrograde
        (5.9, 0.0, 0.0), // in the ecliptic
        (1.0, 2.0, PI / 2.0),
    ];
    for (w, node, i) in cases {
        let pq = orbit::gauss_vectors(w, node, i, orbit::OBLQ_J2000);
        let c = orbit::eq_orbit_consts(node, i, orbit::OBLQ_J2000);
        let (P, Q) = ((pq.Px, pq.Py, pq.Pz), (pq.Qx, pq.Qy, pq.Qz));
        // orthonormal
        assert_close(P.0 * P.0 + P.1 * P.1 + P.2 * P.2, 1.0, 1e-15);
        assert_close(Q.0 * Q.0 + Q.1 * Q.1 + Q.2 * Q.2, 1.0, 1e-15);
        assert_close(P.0 * Q.0 + P.1 * Q.1 + P.2 * Q.2, 0.0, 1e-15);
        for k in 0..36 {
            let v = rad(10.0 * k as f64 - 170.0);
            let a = pq.heliocent_rect(v, 1.7);
            let b = c.heliocent_rect(w, v, 1.7);
            assert_close(a.0, b.0, 1e-14);
            assert_close(a.1, b.1, 1e-14);
            assert_close(a.2, b.2, 1e-14);
        }
    }
}

#[test]
fn ecliptic_and_equatorial_rectangular_coordinates() {
    // Rotating the heliocentric ecliptic coordinates about the x-axis by
    // the obliquity gives the equatorial ones (eq. 33.9).
    let eps = orbit::OBLQ_J2000;
    let el = encke();
    for k in 0..50 {
        let t = encke_T() + 20.0 * (k as f64 - 25.0);
        let (x, y, z) = el.heliocent_ecl_rect(t).unwrap();
        let eq = el.heliocent_eq_rect(t, eps).unwrap();
        assert_close(eq.0, x, 1e-14);
        assert_close(eq.1, y * eps.cos() - z * eps.sin(), 1e-14);
        assert_close(eq.2, y * eps.sin() + z * eps.cos(), 1e-14);

        // heliocentric longitude and latitude
        let (v, r) = el.motion.true_anom_and_rad_vec(t).unwrap();
        let (l, b) = orbit::heliocent_ecl_coords(el.arg_perih, el.long_asc_node, el.incl, v);
        assert_close(l, y.atan2(x).rem_euclid(std::f64::consts::TAU), 1e-14);
        assert_close(b.sin(), z / r, 1e-14);
        assert!((0.0..std::f64::consts::TAU).contains(&l));
    }
}

#[test]
fn heliocentric_ecliptic_coordinates_at_the_nodes() {
    // At the ascending node (u = ω + v = 0) the latitude is 0 and the
    // longitude is Ω; at u = 90° the latitude equals i (direct orbit).
    let (w, node, i) = (rad(40.0), rad(100.0), rad(30.0));
    let (l, b) = orbit::heliocent_ecl_coords(w, node, i, -w);
    assert_close_deg(l.to_degrees(), 100.0, 1e-12);
    assert_close(b, 0.0, 1e-15);
    let (_, b) = orbit::heliocent_ecl_coords(w, node, i, rad(90.0) - w);
    assert_close(b.to_degrees(), 30.0, 1e-12);
    // a retrograde orbit (i = 150°) reaches latitude 30° too
    let (_, b) = orbit::heliocent_ecl_coords(w, node, rad(150.0), rad(90.0) - w);
    assert_close(b.to_degrees(), 30.0, 1e-12);
    // descending node: longitude Ω + 180°
    let (l, b) = orbit::heliocent_ecl_coords(w, node, i, PI - w);
    assert_close_deg(l.to_degrees(), 280.0, 1e-12);
    assert_close(b, 0.0, 1e-15);
}

/// The Earth's heliocentric ecliptic coordinates (J2000.0) implied by the
/// Sun's geocentric equatorial coordinates of example 33.b
fn earth_frm_sun(sun: (f64, f64, f64)) -> (f64, f64, f64) {
    let eps = orbit::OBLQ_J2000;
    let (x, y, z) = (-sun.0, -sun.1, -sun.2);
    let (xe, ye, ze) = (
        x,
        y * eps.cos() + z * eps.sin(),
        -y * eps.sin() + z * eps.cos(),
    );
    let R = (xe * xe + ye * ye + ze * ze).sqrt();
    (ye.atan2(xe), ze.atan2(xe.hypot(ye)), R)
}

#[test]
fn ecliptic_method_agrees_with_equatorial_method() {
    // Meeus, example 33.b computed through ecliptic coordinates instead
    let t = jd(1990, time::Month::Oct, 6.0);
    let earth = earth_frm_sun(SUN);
    let p = orbit::geocent_ecl_pos(&encke(), t, earth).unwrap();

    // convert to equatorial coordinates (eqs. 13.3 and 13.4)
    let eps = orbit::OBLQ_J2000;
    let (l, b) = (p.long, p.lat);
    let asc = (l.sin() * eps.cos() - b.tan() * eps.sin()).atan2(l.cos());
    let dec = (b.sin() * eps.cos() + b.cos() * eps.sin() * l.sin()).asin();
    assert_close_deg(asc.to_degrees(), 158.558965, 1e-5);
    assert_close(dec.to_degrees(), 19.158496, 1e-5);
    assert_close(p.earth_dist, 0.8242811, 1e-7);
    assert_close(p.sun_dist, 0.6525755, 1e-7);
    assert_close(p.elong.to_degrees(), 40.51, 0.005);

    let q = orbit::geocent_eq_pos(&encke(), t, SUN, eps).unwrap();
    assert_close(p.light_time, q.light_time, 1e-14);
    assert_close(p.phase_angl, q.phase_angl, 1e-12);
    assert_close(p.elong, q.elong, 1e-12);

    // eq. 33.5 neglects the Sun's latitude (at most 1.2″ = 6e-6 rad)
    let sun_long = earth.0 + PI;
    assert_close(orbit::elong_frm_ecl(p.long, p.lat, sun_long), p.elong, 6e-6);
}

#[test]
fn light_time_iteration() {
    // A body at rest relative to the Sun, far away: the light-time is
    // Δ × 0.0057755183 days
    assert_close(orbit::light_time(1.0), 0.0057755183, 1e-16);
    assert_close(orbit::light_time(0.910947), 0.0052612, 1e-7); // example 33.a

    // geometric position at t - τ equals the returned one
    let t = jd(1990, time::Month::Oct, 6.0);
    let p = orbit::geocent_eq_pos(&encke(), t, SUN, orbit::OBLQ_J2000).unwrap();
    let h = encke()
        .heliocent_eq_rect(t - p.light_time, orbit::OBLQ_J2000)
        .unwrap();
    let (asc, dec, dist) = orbit::geocent_eq_frm_rect(h, SUN);
    assert_close(asc, p.asc, 1e-12);
    assert_close(dec, p.dec, 1e-12);
    assert_close(dist, p.earth_dist, 1e-12);
}

#[test]
fn motion_from_mean_anomaly_at_epoch() {
    // The same orbit given by a mean anomaly at another epoch
    let n = elliptic::mn_motion(ENCKE_A);
    let epoch = encke_T() - 100.0;
    let el = Elements {
        motion: Motion::Elliptic {
            semimaj_axis: ENCKE_A,
            ecc: ENCKE_E,
            mn_motion: n,
            mn_anom: -100.0 * n,
            epoch,
        },
        ..encke()
    };
    for k in 0..20 {
        let t = epoch + 37.0 * k as f64;
        let a = el.motion.true_anom_and_rad_vec(t).unwrap();
        let b = encke().motion.true_anom_and_rad_vec(t).unwrap();
        assert_close(a.0, b.0, 1e-12);
        assert_close(a.1, b.1, 1e-12);
    }
}

#[test]
fn motion_rejects_invalid_elements() {
    let bad = [
        Motion::elliptic_frm_perih_time(-1.0, 0.5, 0.0),
        Motion::elliptic_frm_perih_time(1.0, 1.0, 0.0),
        Motion::elliptic_frm_perih_time(f64::NAN, 0.5, 0.0),
        Motion::Parabolic {
            perih_dist: 0.0,
            perih_time: 0.0,
        },
        Motion::NearParabolic {
            perih_dist: 1.0,
            ecc: f64::NAN,
            perih_time: 0.0,
        },
    ];
    for m in bad {
        assert_eq!(m.true_anom_and_rad_vec(10.0), None, "{m:?}");
        let el = Elements {
            motion: m,
            ..encke()
        };
        assert_eq!(orbit::geocent_eq_pos(&el, 10.0, SUN, 0.4), None);
        assert_eq!(orbit::geocent_ecl_pos(&el, 10.0, (0.0, 0.0, 1.0)), None);
    }
}

#[test]
fn parabolic_and_near_parabolic_motion() {
    // Meeus, example 34.a and the first case of Landgraf's table (ch. 35)
    let T = jd(1998, time::Month::Apr, 14.4358);
    let t = jd(1998, time::Month::Aug, 5.0);
    let (v, r) = Motion::Parabolic {
        perih_dist: 1.487469,
        perih_time: T,
    }
    .true_anom_and_rad_vec(t)
    .unwrap();
    assert_close(v.to_degrees(), 66.78862, 1e-5);
    assert_close(r, 2.133911, 1e-6);

    let (v, r) = Motion::NearParabolic {
        perih_dist: 0.921326,
        ecc: 1.0,
        perih_time: 0.0,
    }
    .true_anom_and_rad_vec(138.4783)
    .unwrap();
    assert_close(v.to_degrees(), 102.74426, 1e-5);
    assert_close(r, 2.364192, 1e-6);

    // An elliptic orbit with e close to 1 gives the same result through
    // Kepler's equation and through Landgraf's method
    let (q, e) = (0.5871018, 0.9672746);
    let a = elliptic::semimaj_axis(q, e);
    let ell = Motion::elliptic_frm_perih_time(a, e, 0.0);
    let near = Motion::NearParabolic {
        perih_dist: q,
        ecc: e,
        perih_time: 0.0,
    };
    for t in [-30.0, -5.0, 0.0, 1.0, 20.0, 45.0] {
        let (v1, r1) = ell.true_anom_and_rad_vec(t).unwrap();
        let (v2, r2) = near.true_anom_and_rad_vec(t).unwrap();
        common::assert_close_rad(v1, v2, 1e-9);
        assert_close(r1, r2, 1e-9);
    }
}

#[test]
fn equation_of_center_accuracy() {
    // Meeus, chapter 33: greatest error of the series up to e⁵
    let max_err = |e: f64| {
        (1..3600)
            .map(|i| {
                let M = rad(i as f64 * 0.05);
                let E = elliptic::ecc_anom_checked(M, e).unwrap();
                let v = elliptic::true_anom(E, e);
                (elliptic::equation_of_center(M, e) - (v - M))
                    .abs()
                    .to_degrees()
                    * 3600.0
            })
            .fold(0.0, f64::max)
    };
    assert_close(max_err(0.03), 0.0003, 0.00005);
    assert_close(max_err(0.05), 0.007, 0.0005);
    assert_close(max_err(0.10), 0.45, 0.005);
    assert_close(max_err(0.15), 5.0, 0.5);
    assert_close(max_err(0.20), 29.0, 0.5);
    assert_close(max_err(0.25), 111.0, 0.5);
    assert_close(max_err(0.30), 331.0, 1.0);

    // odd in M, zero at perihelion and aphelion, and zero for e = 0
    let e = 0.0934;
    assert_close(
        elliptic::equation_of_center(-1.1, e),
        -elliptic::equation_of_center(1.1, e),
        1e-17,
    );
    assert_close(elliptic::equation_of_center(0.0, e), 0.0, 1e-17);
    assert_close(elliptic::equation_of_center(PI, e), 0.0, 1e-14);
    assert_eq!(elliptic::equation_of_center(1.0, 0.0), 0.0);
}

#[test]
fn radius_vector_series() {
    for &(e, tol) in &[(0.0167, 5e-11), (0.05, 2e-8), (0.0934, 1e-6)] {
        let a = 1.523679;
        for i in -36..=36 {
            let M = rad(10.0 * i as f64);
            let E = elliptic::ecc_anom_checked(M, e).unwrap();
            let r = elliptic::rad_vec_frm_ecc_anom(E, a, e);
            assert_close(elliptic::rad_vec_series(M, a, e), r, tol * a);
        }
    }
    // perihelion and aphelion distances to order e⁵
    let (a, e) = (1.0, 0.01);
    assert_close(elliptic::rad_vec_series(0.0, a, e), a * (1.0 - e), 1e-13);
    assert_close(elliptic::rad_vec_series(PI, a, e), a * (1.0 + e), 1e-13);
}

#[test]
fn example_33c_velocities() {
    // Meeus, example 33.c: comet Halley, 1986
    let (a, e) = (17.9400782, 0.96727426);
    assert_close(elliptic::perih_vel(a, e), 54.52, 0.005);
    assert_close(elliptic::aph_vel(a, e), 0.91, 0.005);
    assert_close(elliptic::vel(1.0, a), 41.53, 0.005);
    // the vis-viva equation gives the perihelion and aphelion velocities too
    let (q, Q) = (elliptic::perih_dist(a, e), elliptic::aph_dist(a, e));
    assert_close(elliptic::vel(q, a), elliptic::perih_vel(a, e), 1e-3);
    assert_close(elliptic::vel(Q, a), elliptic::aph_vel(a, e), 1e-4);
    // Kepler's second law: q Vp = Q Va
    assert_close(
        q * elliptic::perih_vel(a, e),
        Q * elliptic::aph_vel(a, e),
        1e-12,
    );
    // circular orbit of 1 AU: the Earth's mean orbital velocity
    assert_close(elliptic::vel(1.0, 1.0), 29.7847, 1e-4);
}

#[test]
fn example_33d_length_of_ellipse() {
    // Meeus, example 33.d: the orbit of comet Halley is 77.07 AU long, or
    // 11530 million km
    let (a, e) = (17.9400782_f64, 0.96727426_f64);
    let b = a * (1.0 - e * e).sqrt();
    let L = elliptic::length_series(a, b);
    assert_close(L, 77.07, 0.005);
    assert_close(L * 149.5978707, 11530.0, 0.5);
    // the approximate formulae
    assert_close(elliptic::length_ramanujan(a, b), L, L * 2e-4);
    assert_close(elliptic::length(a, b), L, L * 3e-4);
}

#[test]
fn length_of_ellipse_limits() {
    // circle
    assert_close(elliptic::length_series(2.0, 2.0), 4.0 * PI, 1e-15);
    assert_close(elliptic::length_ramanujan(2.0, 2.0), 4.0 * PI, 1e-14);
    assert_close(elliptic::length(2.0, 2.0), 4.0 * PI, 1e-14);
    // segment (e = 1): 4a; Ramanujan's formula is 0.4155% short and the
    // A-G-H formula about 3% long
    assert_eq!(elliptic::length_series(3.0, 0.0), 12.0);
    assert_close(
        elliptic::length_ramanujan(1.0, 0.0) / 4.0 - 1.0,
        -0.004155,
        1e-6,
    );
    assert_close(elliptic::length(1.0, 0.0) / 4.0 - 1.0, 0.0308, 1e-4);
    // nearly flat: the series converges slowly but stays accurate
    let L = elliptic::length_series(1.0, 1e-9);
    assert_close(L, 4.0, 1e-9);
    // symmetric in a and b
    assert_close(
        elliptic::length_series(1.0, 3.0),
        elliptic::length_series(3.0, 1.0),
        1e-14,
    );
    assert_eq!(elliptic::length_series(0.0, 0.0), 0.0);
}

#[test]
fn length_of_ellipse_series_and_error_bounds() {
    // Meeus, chapter 33: the bracket of the series in e is 0.99937 for
    // e = 0.05, 0.99750 for e = 0.10 and 2/π for e = 1, i.e. L / (2πa)
    let ratio = |e: f64| elliptic::length_series(1.0, (1.0 - e * e).sqrt()) / (2.0 * PI);
    assert_close(ratio(0.05), 0.99937, 5e-6);
    assert_close(ratio(0.10), 0.99750, 5e-6);
    assert_close(ratio(1.0), 2.0 / PI, 1e-15);

    // error bounds of the A-G-H formula: < 0.001% for e < 0.88,
    // < 0.01% for e < 0.95, 1% for e = 0.9997
    let rel_err = |e: f64| {
        let b = (1.0 - e * e).sqrt();
        (elliptic::length(1.0, b) / elliptic::length_series(1.0, b) - 1.0).abs()
    };
    for i in 0..88 {
        assert!(rel_err(i as f64 / 100.0) < 1e-5);
    }
    assert!(rel_err(0.95) < 1e-4);
    assert_close(rel_err(0.9997), 0.01, 0.001);

    // Ramanujan's error grows monotonically to 0.4155%
    let ram_err = |e: f64| {
        let b = (1.0 - e * e).sqrt();
        1.0 - elliptic::length_ramanujan(1.0, b) / elliptic::length_series(1.0, b)
    };
    assert!(ram_err(0.5) < 1e-6);
    assert!(ram_err(0.99) < 0.004155);
}

#[test]
fn semimajor_axis_and_mean_motion() {
    // eq. 33.6, with example 33.b
    let q = elliptic::perih_dist(ENCKE_A, ENCKE_E);
    assert_close(elliptic::semimaj_axis(q, ENCKE_E), ENCKE_A, 1e-14);
    assert_close(
        elliptic::mn_motion(ENCKE_A).to_degrees(),
        0.9856076686 / (ENCKE_A * ENCKE_A.sqrt()),
        1e-12,
    );
    // Earth: one revolution in a Gaussian year
    assert_close(2.0 * PI / elliptic::mn_motion(1.0), 365.2568983, 1e-6);
}

#[test]
fn minor_planet_magnitude() {
    // Meeus, eq. 33.14, with the H and G of (1) Ceres from the book
    let (H, G) = (3.34, 0.12);
    let m = asteroid::apparent_mag_HG(H, G, 2.77, 1.85, rad(15.0));
    // value computed independently from eq. 33.14
    assert_close(m, 7.754492466457433, 1e-12);
    // at zero phase angle Φ1 = Φ2 = 1
    assert_close(
        asteroid::apparent_mag_HG(H, G, 2.0, 1.0, 0.0),
        H + 5.0 * 2.0_f64.log10(),
        1e-14,
    );
    // H is the magnitude at r = Δ = 1 AU and zero phase
    assert_close(asteroid::apparent_mag_HG(H, G, 1.0, 1.0, 0.0), H, 1e-14);
    // fainter at larger phase angles; sign of the phase angle ignored
    let m30 = asteroid::apparent_mag_HG(H, G, 2.0, 1.5, rad(30.0));
    let m10 = asteroid::apparent_mag_HG(H, G, 2.0, 1.5, rad(10.0));
    assert!(m30 > m10);
    assert_eq!(asteroid::apparent_mag_HG(H, G, 2.0, 1.5, rad(-30.0)), m30);
    // invalid input
    assert!(asteroid::apparent_mag_HG(H, G, 0.0, 1.5, 0.1).is_nan());
    assert!(asteroid::apparent_mag_HG(H, G, 2.0, -1.0, 0.1).is_nan());
    assert!(asteroid::apparent_mag_HG(H, G, 2.0, 1.0, PI).is_nan());
}

#[test]
fn comet_magnitude() {
    // Meeus, eq. 33.13
    assert_close(asteroid::comet_total_mag(5.5, 10.0, 1.0, 1.0), 5.5, 1e-15);
    assert_close(asteroid::comet_total_mag(5.5, 10.0, 10.0, 1.0), 15.5, 1e-14);
    assert_close(asteroid::comet_total_mag(5.5, 10.0, 1.0, 10.0), 10.5, 1e-14);
    assert_close(
        asteroid::comet_total_mag(6.0, 8.0, 0.6525755, 0.8242811),
        6.0 + 5.0 * 0.8242811_f64.log10() + 8.0 * 0.6525755_f64.log10(),
        1e-14,
    );
}

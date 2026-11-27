//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 33: geocentric and
//! apparent positions of the planets (first method), and the formulae of
//! the second method that live in `planet`.

#![allow(non_snake_case)]

mod common;

use astro::coords::EclPoint;
use astro::planet::{self, Planet};
use common::{assert_close, assert_close_deg, dms, hms};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

const PLANETS: [Planet; 7] = [
    Planet::Mercury,
    Planet::Venus,
    Planet::Mars,
    Planet::Jupiter,
    Planet::Saturn,
    Planet::Uranus,
    Planet::Neptune,
];

const ARCSEC: f64 = 1.0 / 3600.0;

fn arcsec(rad: f64) -> f64 {
    rad.to_degrees() * 3600.0
}

// Meeus, example 33.a: Venus on 1992 December 20 at 0h TD.
//
// The book carries out the reduction with the abridged series of appendix
// III, so its intermediate angles differ from the complete VSOP87 by up to
// about 1". It also quotes the values of the complete theory for the final
// result: α = 21h04m41.454s, δ = -18°53'16.84", true distance 0.91084596.
#[test]
fn geocent_apprnt_pos_example_33a() {
    let pos = planet::geocent_apprnt_pos(&Planet::Venus, 2448976.5).unwrap();

    // True distance (complete VSOP87) and light-time.
    assert_close(pos.geomet_dist, 0.91084596, 1e-8);
    assert_close(pos.light_time, 0.0052612, 1e-7);
    assert_close(pos.dist, 0.910947, 1.5e-6);
    assert_close(pos.light_time, 0.0057755183 * pos.dist, 1e-12);

    // Corrected for light-time only: λ = 313.08102°, β = -2.08474°.
    assert_close_deg(pos.light_time_long.to_degrees(), 313.08102, ARCSEC);
    assert_close(pos.light_time_lat.to_degrees(), -2.08474, ARCSEC);

    // Apparent: λ = 313.08151°, β = -2.08487°.
    assert_close_deg(pos.apprnt_long.to_degrees(), 313.08151, ARCSEC);
    assert_close(pos.apprnt_lat.to_degrees(), -2.08487, ARCSEC);

    // The corrections themselves do not depend on the series used:
    // aberration Δλ = -14.868", Δβ = -0.531"; FK5 Δλ = -0.09027",
    // Δβ = +0.05535"; nutation Δψ = +16.749".
    let d_long = arcsec(pos.apprnt_long - pos.light_time_long);
    let d_lat = arcsec(pos.apprnt_lat - pos.light_time_lat);
    assert_close(d_long, -14.868 - 0.09027 + 16.749, 0.003);
    assert_close(d_lat, -0.531 + 0.05535, 0.003);

    // Final values of the complete VSOP87 theory.
    assert_close_deg(
        pos.apprnt_asc.to_degrees(),
        hms(21, 4, 41.454),
        0.005 / 240.0,
    );
    assert_close(
        pos.apprnt_dec.to_degrees(),
        dms(true, 18, 53, 16.84),
        0.05 * ARCSEC,
    );

    // The geometric position is ahead of the light-time corrected one by
    // the planet's motion during τ.
    assert!(pos.geomet_long > pos.light_time_long);
    assert!(pos.geomet_dist < pos.dist);

    let eq = pos.apprnt_eq_point();
    assert_eq!((eq.asc, eq.dec), (pos.apprnt_asc, pos.apprnt_dec));
    let ecl = pos.apprnt_ecl_point();
    assert_eq!((ecl.long, ecl.lat), (pos.apprnt_long, pos.apprnt_lat));
}

// Meeus, example 18.a: apparent positions of Mercury and Venus at 0h TD,
// 1991 August 5-9, "from an accurate ephemeris" (0.001s, 0.01").
#[test]
fn geocent_apprnt_pos_example_18a() {
    #[rustfmt::skip]
    let mercury = [
        ((10, 24, 30.125), (6, 26, 32.05)),
        ((10, 25, 0.342), (6, 10, 57.72)),
        ((10, 25, 12.515), (5, 57, 33.08)),
        ((10, 25, 6.235), (5, 46, 27.07)),
        ((10, 24, 41.185), (5, 37, 48.45)),
    ];
    #[rustfmt::skip]
    let venus = [
        ((10, 27, 27.175), (4, 4, 41.83)),
        ((10, 26, 32.410), (3, 55, 54.66)),
        ((10, 25, 29.042), (3, 48, 3.51)),
        ((10, 24, 17.191), (3, 41, 10.25)),
        ((10, 22, 57.024), (3, 35, 16.61)),
    ];
    for (p, table) in [(Planet::Mercury, mercury), (Planet::Venus, venus)] {
        for (k, ((h, m, s), (d, am, as_))) in table.into_iter().enumerate() {
            let jd = 2448473.5 + k as f64;
            let pos = planet::geocent_apprnt_pos(&p, jd).unwrap();
            assert_close_deg(pos.apprnt_asc.to_degrees(), hms(h, m, s), 0.003 / 240.0);
            assert_close(
                pos.apprnt_dec.to_degrees(),
                dms(false, d, am, as_),
                0.02 * ARCSEC,
            );
        }
    }
}

// Meeus, example 15.a: apparent positions of Venus at 0h TD, 1988 March
// 19-21, from an accurate ephemeris (0.01s, 0.1").
#[test]
fn geocent_apprnt_pos_example_15a() {
    let venus = [
        (2447239.5, (2, 42, 43.25), (18, 2, 51.4)),
        (2447240.5, (2, 46, 55.51), (18, 26, 27.3)),
        (2447241.5, (2, 51, 7.69), (18, 49, 38.7)),
    ];
    for (jd, (h, m, s), (d, am, as_)) in venus {
        let pos = planet::geocent_apprnt_pos(&Planet::Venus, jd).unwrap();
        assert_close_deg(pos.apprnt_asc.to_degrees(), hms(h, m, s), 0.01 / 240.0);
        assert_close(
            pos.apprnt_dec.to_degrees(),
            dms(false, d, am, as_),
            0.06 * ARCSEC,
        );
    }
}

#[test]
fn geocent_apprnt_pos_of_earth_is_none() {
    assert_eq!(planet::geocent_apprnt_pos(&Planet::Earth, 2451545.0), None);
}

// Meeus, ch. 33: instead of applying the aberration formula, the Earth can
// be taken at t - τ as well. Both methods agree to first order; the
// residual (a few 0.01") reflects the Keplerian approximation of the
// Earth's motion in eq. 23.2.
#[test]
fn aberration_equals_earth_at_retarded_time() {
    for p in PLANETS {
        for jd in [2433282.5, 2448976.5, 2457023.5] {
            let pos = planet::geocent_apprnt_pos(&p, jd).unwrap();
            let tau = pos.light_time;
            let (L0, B0, R0) = planet::heliocent_coords(&Planet::Earth, jd - tau);
            let (L, B, R) = planet::heliocent_coords(&p, jd - tau);
            let (lam, beta, _, _) = planet::geocent_geomet_ecl_coords(L0, B0, R0, L, B, R);
            let (lam, beta) = planet::ecl_coords_to_FK5(jd, lam, beta);
            let (nut_long, _) = astro::nutation::nutation(jd);

            let d_long = arcsec(lam + nut_long - pos.apprnt_long) * pos.apprnt_lat.cos();
            let d_lat = arcsec(beta - pos.apprnt_lat);
            assert!(d_long.abs() < 0.05, "{p:?} {jd}: {d_long}");
            assert!(d_lat.abs() < 0.05, "{p:?} {jd}: {d_lat}");
        }
    }
}

// geocent_apprnt_ecl_coords is the light-time-only part of the reduction.
#[test]
fn geocent_apprnt_ecl_coords_is_light_time_corrected() {
    for p in PLANETS {
        for jd in [2415020.5, 2448976.5, 2460676.5] {
            let pos = planet::geocent_apprnt_pos(&p, jd).unwrap();
            let (ecl, dist) = planet::geocent_apprnt_ecl_coords(&p, jd);
            assert_eq!(ecl.long, pos.light_time_long);
            assert_eq!(ecl.lat, pos.light_time_lat);
            assert_eq!(dist, pos.dist);
        }
    }

    // Meeus, example 33.a (abridged series in the book).
    let (ecl, dist) = planet::geocent_apprnt_ecl_coords(&Planet::Venus, 2448976.5);
    assert_close_deg(ecl.long.to_degrees(), 313.08102, ARCSEC);
    assert_close(ecl.lat.to_degrees(), -2.08474, ARCSEC);
    assert_close(dist, 0.910947, 1.5e-6);
}

// Meeus, example 33.a, using the book's own heliocentric coordinates (A).
#[test]
fn geocent_geomet_ecl_coords_example_33a() {
    let (lambda, beta, delta, tau) = planet::geocent_geomet_ecl_coords(
        88.35704_f64.to_radians(),
        0.00014_f64.to_radians(),
        0.983824,
        26.11428_f64.to_radians(),
        -2.62070_f64.to_radians(),
        0.724603,
    );
    assert!((0.0..TAU).contains(&lambda));
    assert_close(delta, 0.910845, 1e-6);
    assert_close(tau, 0.0052606, 1e-7);
    // tan λ = y/x and tan β = z/√(x²+y²) with the book's x, y, z.
    let (x, y, z): (f64, f64, f64) = (0.621746, -0.664810, -0.033134);
    // The book's x, y, z are rounded to 1e-6 AU, hence the tolerance.
    assert_close_deg(lambda.to_degrees(), y.atan2(x).to_degrees(), 1e-4);
    assert_close(beta.to_degrees(), z.atan2(x.hypot(y)).to_degrees(), 1e-4);
}

// The elongation agrees with the triangle Sun-Earth-planet
// (cos ψ = (R² + Δ² - r²) / 2RΔ, Meeus eq. 33.11); the small residual is
// due to light-time and aberration.
#[test]
fn elongation_matches_distances() {
    for p in PLANETS {
        for k in 0..40 {
            let jd = 2440000.5 + k as f64 * 97.3;
            let pos = planet::geocent_apprnt_pos(&p, jd).unwrap();
            let (_, _, R) = planet::heliocent_coords(&Planet::Earth, jd);
            let (_, _, r) = planet::heliocent_coords(&p, jd - pos.light_time);
            let delta = pos.dist;
            let cos_psi = (R * R + delta * delta - r * r) / (2.0 * R * delta);
            let psi = cos_psi.clamp(-1.0, 1.0).acos();
            assert!((0.0..=PI).contains(&pos.elongation));
            assert!(
                (pos.elongation - psi).to_degrees().abs() < 0.02,
                "{p:?} {jd}: {} vs {}",
                pos.elongation.to_degrees(),
                psi.to_degrees()
            );
        }
    }
    // The inner planets never get far from the Sun.
    for k in 0..200 {
        let jd = 2450000.5 + k as f64 * 11.0;
        let merc = planet::geocent_apprnt_pos(&Planet::Mercury, jd).unwrap();
        let ven = planet::geocent_apprnt_pos(&Planet::Venus, jd).unwrap();
        assert!(merc.elongation.to_degrees() < 28.5);
        assert!(ven.elongation.to_degrees() < 47.9);
    }
}

#[test]
fn elongation_function() {
    let sun = 10f64.to_radians();
    let p = EclPoint {
        long: 70f64.to_radians(),
        lat: 0.0,
    };
    assert_close(planet::elongation(&p, sun).to_degrees(), 60.0, 1e-12);
    // Symmetric in longitude, and works across 0°/360°.
    let p = EclPoint {
        long: 350f64.to_radians(),
        lat: 0.0,
    };
    assert_close(planet::elongation(&p, sun).to_degrees(), 20.0, 1e-9);
    // Opposition and a point at the pole of the ecliptic.
    let p = EclPoint {
        long: 190f64.to_radians(),
        lat: 0.0,
    };
    assert_close(planet::elongation(&p, sun).to_degrees(), 180.0, 1e-6);
    let p = EclPoint {
        long: 123f64.to_radians(),
        lat: FRAC_PI_2,
    };
    assert_close(planet::elongation(&p, sun).to_degrees(), 90.0, 1e-9);
}

// Ranges and continuity of the apparent positions, over two millennia
// either side of J2000 (dates before 1582 and negative years included).
#[test]
fn apprnt_pos_ranges() {
    for p in PLANETS {
        for k in -40..=40 {
            let jd = 2451545.0 + k as f64 * 18262.5 + 0.25;
            let pos = planet::geocent_apprnt_pos(&p, jd).unwrap();
            for long in [
                pos.geomet_long,
                pos.light_time_long,
                pos.apprnt_long,
                pos.apprnt_asc,
            ] {
                assert!((0.0..TAU).contains(&long), "{p:?} {jd}: {pos:?}");
            }
            for lat in [
                pos.geomet_lat,
                pos.light_time_lat,
                pos.apprnt_lat,
                pos.apprnt_dec,
            ] {
                assert!((-FRAC_PI_2..=FRAC_PI_2).contains(&lat), "{p:?} {jd}");
            }
            assert!(pos.light_time > 0.0 && pos.light_time < 0.3);
            // The aberration is at most about 20.5" / cos β, plus the
            // nutation (< 20") and the FK5 offset.
            let shift = arcsec(pos.apprnt_long - pos.light_time_long);
            let shift = (shift + 648000.0).rem_euclid(1296000.0) - 648000.0;
            assert!(shift.abs() < 45.0, "{p:?} {jd}: {shift}");
        }
    }
}

// Meeus, example 33.b (second method): periodic comet Encke on 1990
// October 6.0 TD, elements for the equinox J2000.0. The second call uses
// the anomaly and radius vector for t - τ and gives the astrometric
// position α = 158.558965°, δ = +19.158496°, Δ = 0.8242811.
#[test]
fn geocent_eq_coords_example_33b() {
    let i = 11.94524_f64.to_radians();
    let node = 334.75006_f64.to_radians();
    let w = 186.23352_f64.to_radians();
    let oblq = dms(false, 23, 26, 21.448).to_radians();
    let (X, Y, Z) = (-0.9756732, -0.2003254, -0.0868566);

    // Geometric position: Δ = 0.8243689, τ = 0.00476 day.
    let (_, _, tau) = planet::geocent_eq_coords(
        X,
        Y,
        Z,
        i,
        w,
        node,
        oblq,
        -94.163310_f64.to_radians(),
        0.6524867,
    );
    assert_close(tau / 0.0057755183, 0.8243689, 1e-7);
    assert_close(tau, 0.00476, 5e-6);

    let (asc, dec, tau) = planet::geocent_eq_coords(
        X,
        Y,
        Z,
        i,
        w,
        node,
        oblq,
        -94.171933_f64.to_radians(),
        0.6525755,
    );
    assert_close_deg(asc.to_degrees(), 158.558965, 2e-6);
    assert_close(dec.to_degrees(), 19.158496, 2e-6);
    assert_close(tau / 0.0057755183, 0.8242811, 1e-7);
}

// Meeus, example 33.b, through heliocent_coords_frm_orb_elements: the
// ecliptic coordinates rotated to the equator of J2000.0 must give the
// book's x = +0.2508066, y = +0.4849175, z = +0.3573373.
#[test]
fn heliocent_coords_frm_orb_elements_example_33b() {
    let i = 11.94524_f64.to_radians();
    let node = 334.75006_f64.to_radians();
    let w = 186.23352_f64.to_radians();
    let eps = dms(false, 23, 26, 21.448).to_radians();
    let (v, r) = (-94.163310_f64.to_radians(), 0.6524867);

    let (l, b) = planet::heliocent_coords_frm_orb_elements(i, node, w, v, r);
    assert!((0.0..TAU).contains(&l));
    let x = r * b.cos() * l.cos();
    let y = r * (b.cos() * l.sin() * eps.cos() - b.sin() * eps.sin());
    let z = r * (b.cos() * l.sin() * eps.sin() + b.sin() * eps.cos());
    assert_close(x, 0.2508066, 1e-7);
    assert_close(y, 0.4849175, 1e-7);
    assert_close(z, 0.3573373, 1e-7);
}

// Simple cases of the ecliptic formulae: a circular orbit in the ecliptic
// and a polar orbit.
#[test]
fn heliocent_coords_frm_orb_elements_limits() {
    let (l, b) = planet::heliocent_coords_frm_orb_elements(0.0, 1.0, 2.0, 3.5, 1.0);
    assert_close(l, (1.0f64 + 2.0 + 3.5).rem_euclid(TAU), 1e-12);
    assert_close(b, 0.0, 1e-12);

    // i = 90°, node 0, u = 30°: the body is 30° above the ecliptic, on
    // the node's meridian.
    let (l, b) =
        planet::heliocent_coords_frm_orb_elements(FRAC_PI_2, 0.0, 0.0, 30f64.to_radians(), 2.0);
    assert_close(l, 0.0, 1e-12);
    assert_close(b.to_degrees(), 30.0, 1e-12);

    // Retrograde orbit (i = 180°): motion in the opposite direction.
    let (l, b) = planet::heliocent_coords_frm_orb_elements(PI, 0.0, 0.0, 30f64.to_radians(), 1.0);
    assert_close(l.to_degrees(), 330.0, 1e-9);
    assert_close(b, 0.0, 1e-12);
}

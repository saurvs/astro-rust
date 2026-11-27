//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 32: heliocentric
//! positions of the planets from VSOP87, and the conversion to FK5.

#![allow(non_snake_case)]

mod common;

use astro::planet::{self, Planet};
use common::{assert_close, assert_close_deg, dms};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

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

fn arcsec(rad: f64) -> f64 {
    rad.to_degrees() * 3600.0
}

// Meeus, example 32.a: Venus on 1992 December 20 at 0h TD. The book uses the
// abridged series of appendix III; the complete VSOP87 series used here
// must agree to within the error of the abridged series (about 1").
#[test]
fn heliocent_coords_example_32a() {
    let (L, B, R) = planet::heliocent_coords(&Planet::Venus, 2448976.5);
    assert_close_deg(L.to_degrees(), 26.11428, 1.0 / 3600.0);
    assert_close(B.to_degrees(), -2.62070, 1.0 / 3600.0);
    assert_close(R, 0.724603, 2e-6);
}

// Regression: the latitude used to be reduced to [0, 2π), so a planet
// south of the ecliptic had a latitude close to 2π.
#[test]
fn heliocent_lat_is_signed() {
    let (_, B, _) = planet::heliocent_coords(&Planet::Venus, 2448976.5);
    assert!(B < 0.0, "{B}");
    assert!(B > -FRAC_PI_2);

    let mut seen_south = [false; 8];
    let mut seen_north = [false; 8];
    for (n, p) in PLANETS.iter().enumerate() {
        for k in 0..400 {
            // 1600 years of dates, half of them before J2000.
            let jd = 2451545.0 + (k as f64 - 200.0) * 1461.17;
            let (L, B, R) = planet::heliocent_coords(p, jd);
            assert!((0.0..TAU).contains(&L), "{p:?} {jd}: L = {L}");
            assert!((-FRAC_PI_2..=FRAC_PI_2).contains(&B), "{p:?} {jd}: B = {B}");
            // No planet strays more than 8° from the ecliptic of date.
            assert!(B.abs() < 8f64.to_radians(), "{p:?} {jd}: B = {B}");
            assert!(R > 0.3 && R < 31.0, "{p:?} {jd}: R = {R}");
            seen_south[n] |= B < 0.0;
            seen_north[n] |= B > 0.0;
        }
    }
    assert!(seen_south.iter().all(|&s| s) && seen_north.iter().all(|&s| s));
}

// Meeus, example 25.a: the book gives the Sun's position on 1992 October 13
// at 0h TD from the complete VSOP87 theory (FK5, mean equinox of date):
// longitude 199°54'26.18", latitude +0.72", radius vector 0.99760853 AU.
#[test]
fn earth_complete_vsop87_example_25a() {
    let jd = 2448908.5;
    let (L, B, R) = planet::heliocent_coords(&Planet::Earth, jd);
    let (sun_long, sun_lat) = planet::ecl_coords_to_FK5(jd, L + PI, -B);

    assert_close_deg(
        sun_long.to_degrees(),
        dms(false, 199, 54, 26.18),
        0.006 / 3600.0,
    );
    assert_close(arcsec(sun_lat), 0.72, 0.006);
    assert_close(R, 0.99760853, 1.5e-8);
}

// Meeus, example 32.b: the heliocentric longitude of Saturn on 1999 July 26
// at 0h TD is 39.9723901° (FK5, mean equinox of date), which the book says
// is the value of the complete VSOP87 theory. The polynomial the book
// starts from is accurate to 0.02", hence the tolerance.
#[test]
fn saturn_complete_vsop87_example_32b() {
    let jd = 2451385.5;
    let (L, B, _) = planet::heliocent_coords(&Planet::Saturn, jd);
    let (L_fk5, _) = planet::ecl_coords_to_FK5(jd, L, B);
    assert_close_deg(L_fk5.to_degrees(), 39.9723901, 0.02 / 3600.0);
}

// Meeus, example 33.a: the FK5 corrections (eq. 32.3) to the geocentric
// position of Venus are Δλ = -0.09027" and Δβ = +0.05535".
#[test]
fn ecl_coords_to_FK5_example_33a() {
    let jd = 2448976.5;
    let (long, lat) = (313.07689_f64.to_radians(), -2.08489_f64.to_radians());
    let (long_fk5, lat_fk5) = planet::ecl_coords_to_FK5(jd, long, lat);
    assert_close(arcsec(long_fk5 - long), -0.09027, 0.00002);
    assert_close(arcsec(lat_fk5 - lat), 0.05535, 0.00002);
    assert_close_deg(long_fk5.to_degrees(), 313.07686, 1e-5);
    assert_close(lat_fk5.to_degrees(), -2.08487, 1e-5);
}

// Meeus, example 25.b: for the Sun (β ≈ 0) the longitude correction is the
// constant -0.09033" and Δβ = -0.023".
#[test]
fn ecl_coords_to_FK5_example_25b() {
    let jd = 2448908.5;
    let (long, lat) = (199.907372_f64.to_radians(), 0.644_f64 / 3600.0);
    let (long_fk5, lat_fk5) = planet::ecl_coords_to_FK5(jd, long, lat.to_radians());
    assert_close(arcsec(long_fk5 - long), -0.09033, 0.00001);
    assert_close(arcsec(lat_fk5) - 0.644, -0.023, 0.0005);
}

// The FK5 correction is tiny everywhere: under 0.1" in latitude and, away
// from the poles, under 0.2" in longitude.
#[test]
fn FK5_correction_is_small() {
    for k in 0..72 {
        let long = (k as f64 * 5.0).to_radians();
        for lat_deg in [-60.0, -7.0, 0.0, 3.0, 45.0] {
            let lat = f64::to_radians(lat_deg);
            for jd in [1000000.5, 2451545.0, 2816787.5] {
                let (l2, b2) = planet::ecl_coords_to_FK5(jd, long, lat);
                assert!(arcsec(l2 - long).abs() < 0.2);
                assert!(arcsec(b2 - lat).abs() < 0.06);
            }
        }
    }
}

// The heliocentric longitude moves smoothly: no jumps except the wrap
// through 2π, and the daily motion is close to the mean motion.
#[test]
fn heliocent_coords_are_continuous() {
    // Approximate largest daily motions in longitude (degrees/day).
    let max_rate = [6.4, 1.65, 1.03, 0.69, 0.10, 0.045, 0.014, 0.0065];
    for (n, p) in PLANETS.iter().enumerate() {
        let mut prev = planet::heliocent_coords(p, 2451000.0);
        for k in 1..2000 {
            let jd = 2451000.0 + k as f64;
            let cur = planet::heliocent_coords(p, jd);
            let d_long = (cur.0 - prev.0).rem_euclid(TAU).to_degrees();
            assert!(d_long > 0.0 && d_long < max_rate[n], "{p:?} {jd}: {d_long}");
            assert!((cur.1 - prev.1).to_degrees().abs() < 0.8, "{p:?} {jd}");
            assert!((cur.2 - prev.2).abs() < 0.01, "{p:?} {jd}");
            prev = cur;
        }
    }
}

// The radius vector stays between perihelion and aphelion of the mean
// orbit, allowing for the perturbations, also far from J2000 (here from
// about 1000 B.C. to A.D. 5000).
#[test]
fn rad_vec_within_orbit_bounds() {
    for p in PLANETS {
        for k in -300..=300 {
            let jd = 2451545.0 + k as f64 * 3652.5 + 0.7;
            let el = planet::mean_orb_elements_of_date(&p, jd);
            let (_, _, R) = planet::heliocent_coords(&p, jd);
            let slack = 0.01 * el.semimaj_axis;
            let q = el.semimaj_axis * (1.0 - el.ecc) - slack;
            let Q = el.semimaj_axis * (1.0 + el.ecc) + slack;
            assert!(R > q && R < Q, "{p:?} {jd}: {R} not in [{q}, {Q}]");
        }
    }
}

// The Earth's heliocentric latitude in the ecliptic of date is tiny
// (under about 1.2").
#[test]
fn earth_lat_is_tiny() {
    for k in 0..1000 {
        let jd = 2440000.5 + k as f64 * 13.7;
        let (_, B, _) = planet::heliocent_coords(&Planet::Earth, jd);
        assert!(arcsec(B).abs() < 1.3, "{jd}: {}", arcsec(B));
    }
}

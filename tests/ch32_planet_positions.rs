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

// The reference values distributed with VSOP87 by its authors (file
// `vsop87.chk`, version D, IMCCE), for J2000.0, 1700 and 1100 (TDB): the
// complete series must reproduce them to the 1e-10 printed. Together with
// a term-by-term comparison of `src/planet/VSOPD_87` with the original
// VSOP87D files (all 31 577 terms identical), this covers every planet and
// every power of τ, which the book's examples cannot.
#[test]
fn heliocent_coords_match_vsop87_check_values() {
    #[rustfmt::skip]
    let cases = [
        (Planet::Mercury, 2451545.0, 4.4293481036, -0.0527573409, 0.4664714751),
        (Planet::Mercury, 2341970.0, 0.1910149587, -0.0682441256, 0.3381563139),
        (Planet::Mercury, 2122820.0, 5.0128397764, -0.1143275808, 0.4352063237),
        (Planet::Venus, 2451545.0, 3.1870221833, 0.0569782849, 0.7202129253),
        (Planet::Venus, 2341970.0, 5.3115708036, -0.0455979904, 0.7283407528),
        (Planet::Venus, 2122820.0, 3.3145399295, 0.0505016053, 0.7215819783),
        (Planet::Earth, 2451545.0, 1.7519238681, -0.0000039656, 0.9833276819),
        (Planet::Earth, 2341970.0, 1.7134419105, 0.0000025051, 0.9831498441),
        (Planet::Earth, 2122820.0, 1.6367193623, -0.0000031292, 0.9830331815),
        (Planet::Mars, 2451545.0, 6.2735389983, -0.0247779824, 1.3912076925),
        (Planet::Mars, 2341970.0, 2.9166648690, 0.0280268149, 1.6584697082),
        (Planet::Mars, 2122820.0, 2.7367104344, 0.0295522719, 1.6571002307),
        (Planet::Jupiter, 2451545.0, 0.6334614186, -0.0205001039, 4.9653813154),
        (Planet::Jupiter, 2341970.0, 4.8888943125, -0.0011098085, 5.1888133656),
        (Planet::Jupiter, 2122820.0, 1.2695066546, -0.0075335740, 5.1193587362),
        (Planet::Saturn, 2451545.0, 0.7980038761, -0.0401984149, 9.1838483715),
        (Planet::Saturn, 2341970.0, 5.8113963637, -0.0291472787, 9.7629994924),
        (Planet::Saturn, 2122820.0, 3.5570108069, 0.0435371139, 9.8669939498),
        (Planet::Uranus, 2451545.0, 5.5225485803, -0.0119527838, 19.9240482667),
        (Planet::Uranus, 2341970.0, 1.7242204720, 0.0059836565, 18.7966208854),
        (Planet::Uranus, 2122820.0, 0.7007226224, -0.0065610611, 19.5612078271),
        (Planet::Neptune, 2451545.0, 5.3045629252, 0.0042236789, 30.1205328392),
        (Planet::Neptune, 2341970.0, 0.0815199679, -0.0260752533, 29.8685860491),
        (Planet::Neptune, 2122820.0, 2.2124988267, 0.0027498093, 30.0653693610),
    ];
    for (planet, jd, l, b, r) in cases {
        let (L, B, R) = planet::heliocent_coords(&planet, jd);
        let dl = (L - l + PI).rem_euclid(TAU) - PI;
        assert!(dl.abs() < 1.5e-10, "{planet:?} {jd}: L {L} vs {l}");
        assert!((B - b).abs() < 1e-10, "{planet:?} {jd}: B {B} vs {b}");
        assert!((R - r).abs() < 1e-10, "{planet:?} {jd}: R {R} vs {r}");
    }
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

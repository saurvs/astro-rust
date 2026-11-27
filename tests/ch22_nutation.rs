// Meeus, Astronomical Algorithms, 2nd ed., chapter 22: nutation and the
// obliquity of the ecliptic (and eq. 23.1, nutation in equatorial
// coordinates).

#![allow(non_snake_case)]

mod common;

use astro::*;
use common::{assert_close, assert_close_deg, dms};

fn to_asec(x: f64) -> f64 {
    x.to_degrees() * 3600.0
}

fn asec(x: f64) -> f64 {
    (x / 3600.0).to_radians()
}

/// Equatorial to ecliptical coordinates (local helper)
fn ecl_frm_eq(asc: f64, dec: f64, eps: f64) -> (f64, f64) {
    let long = (asc.sin() * eps.cos() + dec.tan() * eps.sin()).atan2(asc.cos());
    let lat = (dec.sin() * eps.cos() - dec.cos() * eps.sin() * asc.sin()).asin();
    (long.rem_euclid(std::f64::consts::TAU), lat)
}

// Meeus, example 22.a: 1987 April 10, 0h TD
#[test]
fn example_22_a() {
    let JD = 2446895.5;
    assert_close(time::julian_cent(JD), -0.127296372348, 1e-12);

    let a = nutation::nutation_args(JD);
    assert_close_deg(a.mn_elong_moon.to_degrees(), 136.9623, 0.00005);
    assert_close_deg(a.mn_anom_sun.to_degrees(), 94.9792, 0.00005);
    assert_close_deg(a.mn_anom_moon.to_degrees(), 229.2784, 0.00005);
    assert_close_deg(a.arg_lat_moon.to_degrees(), 143.4079, 0.00005);
    assert_close_deg(a.long_asc_node_moon.to_degrees(), 11.2531, 0.00005);

    let (nut_in_long, nut_in_oblq) = nutation::nutation(JD);
    assert_close(to_asec(nut_in_long), -3.788, 0.0005);
    assert_close(to_asec(nut_in_oblq), 9.443, 0.0005);

    let eps0 = ecliptic::mn_oblq_IAU(JD);
    assert_close(to_asec(eps0), dms(false, 23, 26, 27.407) * 3600.0, 0.0005);
    assert_close(
        to_asec(eps0 + nut_in_oblq),
        dms(false, 23, 26, 36.850) * 3600.0,
        0.001,
    );
    assert_close(
        to_asec(ecliptic::tru_oblq(JD)),
        dms(false, 23, 26, 36.850) * 3600.0,
        0.001,
    );

    // Laskar's formula gives the same mean obliquity there
    assert_close(
        to_asec(ecliptic::mn_oblq_laskar(JD)),
        dms(false, 23, 26, 27.407) * 3600.0,
        0.0005,
    );
}

#[test]
fn nutation_args_are_normalised() {
    for JD in [0.0, 1e6, 2451545.0, 2816787.5, -1e6] {
        let a = nutation::nutation_args(JD);
        for x in [
            a.mn_elong_moon,
            a.mn_anom_sun,
            a.mn_anom_moon,
            a.arg_lat_moon,
            a.long_asc_node_moon,
        ] {
            assert!((0.0..std::f64::consts::TAU).contains(&x), "{JD}: {x}");
        }
    }
}

// The size of the nutation is bounded by the sum of the amplitudes; over
// two thousand years |Δψ| stays below 19.2" and |Δε| below 10.0".
#[test]
fn nutation_bounds() {
    let mut JD = precess::jd_frm_julian_epoch(1000.0);
    while JD < precess::jd_frm_julian_epoch(3000.0) {
        let (p, e) = nutation::nutation(JD);
        assert!(to_asec(p).abs() < 19.2 && to_asec(e).abs() < 10.0, "{JD}");
        JD += 17.3;
    }
}

// The main term has a period of 6798.4 days (18.6 years): Δψ ≈ -17.2" sin Ω
#[test]
fn nutation_main_term() {
    let mut JD = precess::jd_frm_julian_epoch(1990.0);
    while JD < precess::jd_frm_julian_epoch(2010.0) {
        let om = nutation::nutation_args(JD).long_asc_node_moon;
        let (p, e) = nutation::nutation(JD);
        assert!((to_asec(p) + 17.2 * om.sin()).abs() < 2.0);
        assert!((to_asec(e) - 9.2 * om.cos()).abs() < 1.0);
        JD += 29.0;
    }
}

// The low accuracy expressions are good to about 0.5" in Δψ and 0.1" in
// Δε (measured: 0.34" and 0.09" between 1800 and 2200).
#[test]
fn nutation_low_accuracy() {
    let mut JD = precess::jd_frm_julian_epoch(1800.0);
    while JD < precess::jd_frm_julian_epoch(2200.0) {
        let (p, e) = nutation::nutation(JD);
        let (pl, el) = nutation::nutation_low_accuracy(JD);
        assert!((to_asec(p) - to_asec(pl)).abs() < 0.5, "{JD}");
        assert!((to_asec(e) - to_asec(el)).abs() < 0.1, "{JD}");
        JD += 3.7;
    }

    // Example 22.a, to the stated accuracy
    let (p, e) = nutation::nutation_low_accuracy(2446895.5);
    assert_close(to_asec(p), -3.788, 0.5);
    assert_close(to_asec(e), 9.443, 0.1);
}

// ---------------------------------------------------------------------------
// Obliquity
// ---------------------------------------------------------------------------

#[test]
fn mn_oblq_at_j2000() {
    let eps = dms(false, 23, 26, 21.448) * 3600.0;
    assert_close(to_asec(ecliptic::mn_oblq_IAU(precess::J2000)), eps, 1e-9);
    assert_close(to_asec(ecliptic::mn_oblq_laskar(precess::J2000)), eps, 1e-9);
}

// IAU 1980 is in Julian centuries: about -46.8" per century (with the
// bug fixed by PR #20 it changed by only -0.47" per century).
#[test]
fn mn_oblq_IAU_rate() {
    let d = ecliptic::mn_oblq_IAU(precess::J2000 + 36525.0) - ecliptic::mn_oblq_IAU(precess::J2000);
    assert_close(to_asec(d), -46.815 - 0.00059 + 0.001813, 1e-9);
}

// The two formulae agree to better than 0.3" within 1000 years of 2000
// (they drift apart further out, where the IAU one loses accuracy).
#[test]
fn laskar_and_iau_agree_near_2000() {
    for year in [1000.0, 1500.0, 1900.0, 2000.0, 2100.0, 2500.0, 3000.0] {
        let JD = precess::jd_frm_julian_epoch(year);
        let d = to_asec(ecliptic::mn_oblq_laskar(JD) - ecliptic::mn_oblq_IAU(JD));
        assert!(d.abs() < 0.3, "{year}: {d}");
    }
}

// Meeus, chapter 22: Laskar's ε0 has a maximum of 24°14'07" about the year
// -7530 and a minimum of 22°36'41" about +12030; outside |U| < 1 it
// diverges (90° at U = +2.834).
#[test]
fn laskar_extremes_and_range() {
    let jd = |year: f64| precess::jd_frm_julian_epoch(year);

    let max = ecliptic::mn_oblq_laskar(jd(-7530.0));
    assert_close(to_asec(max), dms(false, 24, 14, 7.0) * 3600.0, 0.5);
    assert!(max > ecliptic::mn_oblq_laskar(jd(-7330.0)));
    assert!(max > ecliptic::mn_oblq_laskar(jd(-7730.0)));

    let min = ecliptic::mn_oblq_laskar(jd(12030.0));
    assert_close(to_asec(min), dms(false, 22, 36, 41.0) * 3600.0, 0.5);
    assert!(min < ecliptic::mn_oblq_laskar(jd(11830.0)));
    assert!(min < ecliptic::mn_oblq_laskar(jd(12230.0)));

    let diverged = ecliptic::mn_oblq_laskar(precess::J2000 + 2.834 * 3652500.0);
    assert_close(diverged.to_degrees(), 90.0, 0.01);

    // The checked variant refuses dates outside the valid span
    assert_eq!(ecliptic::mn_oblq_laskar_checked(jd(-7530.0)), Some(max));
    assert!(ecliptic::mn_oblq_laskar_checked(jd(11999.0)).is_some());
    assert!(ecliptic::mn_oblq_laskar_checked(jd(-7999.0)).is_some());
    assert_eq!(ecliptic::mn_oblq_laskar_checked(jd(12030.0)), None);
    assert_eq!(ecliptic::mn_oblq_laskar_checked(jd(-8030.0)), None);
    assert_eq!(ecliptic::mn_oblq_laskar_checked(f64::NAN), None);
    assert_eq!(ecliptic::mn_oblq_laskar_checked(f64::INFINITY), None);
}

#[test]
fn tru_oblq_is_mean_plus_nutation() {
    for JD in [2415020.5, 2446895.5, 2462088.69] {
        assert_close(
            ecliptic::tru_oblq(JD),
            ecliptic::mn_oblq_IAU(JD) + nutation::nutation(JD).1,
            1e-15,
        );
    }
}

// ---------------------------------------------------------------------------
// Nutation in equatorial coordinates (eq. 23.1)
// ---------------------------------------------------------------------------

// Meeus, examples 23.a and 23.b (θ Persei, 2028 November 13.19 TD)
#[test]
fn nutation_in_eq_coords_examples_23() {
    let (psi, eps, oblq) = (asec(14.861), asec(2.705), 23.436_f64.to_radians());

    // Example 23.a: mean place 41.5472°, +49.3485°
    let p = coords::EqPoint {
        asc: 41.5472_f64.to_radians(),
        dec: 49.3485_f64.to_radians(),
    };
    let (da, dd) = nutation::nutation_in_eq_coords(&p, psi, eps, oblq);
    assert_close(to_asec(da), 15.843, 0.001);
    assert_close(to_asec(dd), 6.218, 0.001);

    // Example 23.b: 41.5555635°, +49.3503415°
    let p = coords::EqPoint {
        asc: 41.5555635_f64.to_radians(),
        dec: 49.3503415_f64.to_radians(),
    };
    let (da, dd) = nutation::nutation_in_eq_coords(&p, psi, eps, oblq);
    assert_close(to_asec(da), 15.844, 0.001);
    assert_close(to_asec(dd), 6.217, 0.001);
    assert_close(da.to_degrees(), 0.0044011, 1e-7);
    assert_close(dd.to_degrees(), 0.0017270, 1e-7);
}

// The rigorous rotation agrees with the first-order formulae away from
// the poles.
#[test]
fn nutate_eq_coords_agrees_with_first_order() {
    let JD = 2462088.69;
    let (psi, eps) = nutation::nutation(JD);
    let eps0 = ecliptic::mn_oblq_IAU(JD);

    for (asc, dec) in [
        (0.72, 0.86),
        (2.0, -0.5),
        (3.5, 0.0),
        (5.0, 1.0),
        (6.2, -1.0),
    ] {
        let (a, d) = nutation::nutate_eq_coords(asc, dec, psi, eps, eps0);
        let (da, dd) =
            nutation::nutation_in_eq_coords(&coords::EqPoint { asc, dec }, psi, eps, eps0 + eps);
        let diff_a = (a - asc + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
            - std::f64::consts::PI;
        assert_close(to_asec(diff_a), to_asec(da), 0.001);
        assert_close(to_asec(d - dec), to_asec(dd), 0.001);
    }
}

// Near the pole the first-order formulae fail, the rotation does not: the
// result, turned back into ecliptical coordinates with the true
// obliquity, has the same latitude and a longitude larger by Δψ.
#[test]
fn nutate_eq_coords_near_pole() {
    let (psi, eps) = (asec(-15.0), asec(8.0));
    let eps0 = 23.44_f64.to_radians();

    for (asc, dec) in [
        (0.66, 1.5577),
        (3.0, -1.5700),
        (4.0, std::f64::consts::FRAC_PI_2 - 1e-4),
        (0.0, 0.0),
    ] {
        let (a, d) = nutation::nutate_eq_coords(asc, dec, psi, eps, eps0);
        assert!((0.0..std::f64::consts::TAU).contains(&a));

        let (l0, b0) = ecl_frm_eq(asc, dec, eps0);
        let (l1, b1) = ecl_frm_eq(a, d, eps0 + eps);
        let dl = (l1 - l0 + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
            - std::f64::consts::PI;
        assert_close(to_asec(b1), to_asec(b0), 1e-5);
        assert_close(to_asec(dl * b0.cos()), to_asec(psi * b0.cos()), 1e-4);
    }

    // No nutation, no change
    let (a, d) = nutation::nutate_eq_coords(1.0, 0.5, 0.0, 0.0, eps0);
    assert_close(a, 1.0, 1e-15);
    assert_close(d, 0.5, 1e-15);
}

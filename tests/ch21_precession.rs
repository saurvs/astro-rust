// Meeus, Astronomical Algorithms, 2nd ed., chapter 21: precession.

#![allow(non_snake_case)]

mod common;

use astro::*;
use common::{arcsec, assert_close, assert_close_deg, dms, hms};

/// Radians to arcseconds
fn to_asec(x: f64) -> f64 {
    x.to_degrees() * 3600.0
}

/// Radians of right ascension to seconds of time
fn to_sec_time(x: f64) -> f64 {
    x.to_degrees() * 240.0
}

/// Right ascension in radians to (h, m, s)
fn hms_of(asc: f64) -> (i64, i64, f64) {
    let s = asc.to_degrees().rem_euclid(360.0) * 240.0;
    let h = (s / 3600.0).floor();
    let m = ((s - h * 3600.0) / 60.0).floor();
    (h as i64, m as i64, s - h * 3600.0 - m * 60.0)
}

/// Equatorial to ecliptical coordinates (local helper)
fn ecl_frm_eq(asc: f64, dec: f64, eps: f64) -> (f64, f64) {
    let long = (asc.sin() * eps.cos() + dec.tan() * eps.sin()).atan2(asc.cos());
    let lat = (dec.sin() * eps.cos() - dec.cos() * eps.sin() * asc.sin()).asin();
    (long.rem_euclid(std::f64::consts::TAU), lat)
}

/// Ecliptical to equatorial coordinates (local helper)
fn eq_frm_ecl(long: f64, lat: f64, eps: f64) -> (f64, f64) {
    let asc = (long.sin() * eps.cos() - lat.tan() * eps.sin()).atan2(long.cos());
    let dec = (lat.sin() * eps.cos() + lat.cos() * eps.sin() * long.sin()).asin();
    (asc.rem_euclid(std::f64::consts::TAU), dec)
}

/// Angular separation of two points, in arcseconds
fn sep_asec(a1: f64, d1: f64, a2: f64, d2: f64) -> f64 {
    let c = d1.sin() * d2.sin() + d1.cos() * d2.cos() * (a1 - a2).cos();
    let x = (d1.cos() * d2.sin() - d1.sin() * d2.cos() * (a2 - a1).cos())
        .hypot(d2.cos() * (a2 - a1).sin());
    to_asec(x.atan2(c))
}

// ---------------------------------------------------------------------------
// Low accuracy (eq. 21.1)
// ---------------------------------------------------------------------------

// Meeus, example 21.a: Regulus, 2000.0 -> 1978.0
#[test]
fn example_21_a_annual_precession() {
    let asc = hms(10, 8, 22.3).to_radians();
    let dec = dms(false, 11, 58, 2.0).to_radians();

    let (d_asc, d_dec) = precess::annual_precess(asc, dec, precess::J2000);
    assert_close(to_sec_time(d_asc), 3.208, 0.0005);
    assert_close(to_asec(d_dec), -17.71, 0.005);
}

// Meeus, example 21.a: the full low accuracy reduction with proper motion
#[test]
fn example_21_a_reduction() {
    let asc0 = hms(10, 8, 22.3).to_radians();
    let dec0 = dms(false, 11, 58, 2.0).to_radians();
    let pm_asc = (hms(0, 0, -0.0169)).to_radians();
    let pm_dec = arcsec(0.006).to_radians();
    let years = -22.0;
    let JD = precess::J2000 + years * 365.25;

    let (asc, dec) = precess::precess_eq_coords_low_accuracy(asc0, dec0, precess::J2000, JD);
    let (asc, dec) = (asc + pm_asc * years, dec + pm_dec * years);

    // 10h07m12.1s, +12°04'31"
    assert_close(
        to_sec_time(asc),
        to_sec_time(hms(10, 7, 12.1).to_radians()),
        0.05,
    );
    assert_close(to_asec(dec), dms(false, 12, 4, 31.0) * 3600.0, 1.0);

    // Over two decades the low accuracy method agrees with the rigorous
    // one to within a second of arc (here 0.67", from the neglected
    // change of the rates along the way)
    let (asc_r, dec_r) = precess::precess_eq_coords(asc0, dec0, precess::J2000, JD);
    let (asc_l, dec_l) = precess::precess_eq_coords_low_accuracy(asc0, dec0, precess::J2000, JD);
    assert!(sep_asec(asc_r, dec_r, asc_l, dec_l) < 1.0);
}

#[test]
fn low_accuracy_wraps_right_ascension() {
    // A point just west of 0h moves east past 0h within a few years
    let (asc, _) = precess::precess_eq_coords_low_accuracy(
        hms(23, 59, 59.0).to_radians(),
        0.3,
        precess::J2000,
        precess::J2000 + 3652.5,
    );
    assert!((0.0..std::f64::consts::TAU).contains(&asc));
    assert!(asc < 0.01);
}

// ---------------------------------------------------------------------------
// Rigorous method (eqs. 21.2 to 21.4)
// ---------------------------------------------------------------------------

// Meeus, example 21.b: θ Persei, J2000.0 -> 2028 November 13.19 TD
#[test]
fn example_21_b() {
    let JD = 2462088.69;
    assert_close((JD - precess::J2000) / 36525.0, 0.288670500, 1e-9);

    // Proper motion over 28.86705 years
    let years = (JD - precess::J2000) / 365.25;
    assert_close(years, 28.86705, 1e-5);
    let asc0 = hms(2, 44, 11.986) + hms(0, 0, 0.03425) * years;
    let dec0 = dms(false, 49, 13, 42.48) + arcsec(-0.0895) * years;
    // Meeus prints 41.054063 and 49.227750 from rounded values, and the
    // unrounded ones in example 23.b
    assert_close(asc0, 41.0540613, 1e-7);
    assert_close(dec0, 49.2277489, 1e-7);

    let angles = precess::eq_precess_angles(precess::J2000, JD);
    assert_close(to_asec(angles.zeta), 665.7627, 0.0001);
    assert_close(to_asec(angles.z), 665.8288, 0.0001);
    assert_close(to_asec(angles.theta), 578.5489, 0.0001);

    // Precess the values as printed (41.054063, 49.227750)
    let (asc, dec) = precess::precess_eq_coords(
        41.054063_f64.to_radians(),
        49.227750_f64.to_radians(),
        precess::J2000,
        JD,
    );
    assert_close(asc.to_degrees(), 41.547214, 1e-6);
    assert_close(dec.to_degrees(), 49.348483, 1e-6);

    // 2h46m11.331s, +49°20'54.54"
    let (h, m, s) = hms_of(asc);
    assert_eq!((h, m), (2, 46));
    assert_close(s, 11.331, 0.0005);
    assert_close(
        dec.to_degrees() * 3600.0,
        dms(false, 49, 20, 54.54) * 3600.0,
        0.005,
    );
}

// Meeus, exercise after example 21.b: Polaris, J2000.0 -> B1900.0,
// J2050.0 and J2100.0, close to the celestial pole.
#[test]
fn polaris_exercise() {
    let asc0 = hms(2, 31, 48.704);
    let dec0 = dms(false, 89, 15, 50.72);
    let pm_asc = hms(0, 0, 0.19877);
    let pm_dec = arcsec(-0.0152);

    let cases = [
        (precess::B1900, (1, 22, 33.90), (88, 46, 26.18)),
        (
            precess::jd_frm_julian_epoch(2050.0),
            (3, 48, 16.43),
            (89, 27, 15.38),
        ),
        (
            precess::jd_frm_julian_epoch(2100.0),
            (5, 53, 29.17),
            (89, 32, 22.18),
        ),
    ];

    for (JD, (h, m, s), (d, dm, ds)) in cases {
        let years = (JD - precess::J2000) / 365.25;
        let (asc, dec) = precess::precess_eq_coords(
            (asc0 + pm_asc * years).to_radians(),
            (dec0 + pm_dec * years).to_radians(),
            precess::J2000,
            JD,
        );

        let (hh, mm, ss) = hms_of(asc);
        assert_eq!((hh, mm), (h, m), "epoch {JD}");
        assert_close(ss, s, 0.01);
        assert_close(
            dec.to_degrees() * 3600.0,
            dms(false, d, dm, ds) * 3600.0,
            0.01,
        );
    }
}

// The coefficients of eq. 21.2 for a starting epoch other than J2000.0
// (T = 1, t = 1), evaluated by hand from the printed polynomials. Before
// the fix, the T-terms of z and θ had the wrong sign.
#[test]
fn eq_precess_angles_with_nonzero_T() {
    let JD1 = precess::J2000 + 36525.0;
    let JD2 = JD1 + 36525.0;
    let a = precess::eq_precess_angles(JD1, JD2);

    let zeta = (2306.2181 + 1.39656 - 0.000139) + (0.30188 - 0.000344) + 0.017998;
    let z = (2306.2181 + 1.39656 - 0.000139) + (1.09468 + 0.000066) + 0.018203;
    let theta = (2004.3109 - 0.85330 - 0.000217) - (0.42665 + 0.000217) - 0.041833;

    assert_close(to_asec(a.zeta), zeta, 1e-7);
    assert_close(to_asec(a.z), z, 1e-7);
    assert_close(to_asec(a.theta), theta, 1e-7);
}

// Precessing in two steps must give the same as precessing in one step
// (within the accuracy of the polynomials); with the wrong sign of the
// 0.85330 T term of θ the two-step result was off by about 1.7".
#[test]
fn precess_eq_coords_composes() {
    let stars = [(0.3, 0.9), (2.0, -0.4), (4.5, 1.2), (6.0, -1.3)];
    let (j1, j2, j3) = (
        precess::J2000,
        precess::jd_frm_julian_epoch(2100.0),
        precess::jd_frm_julian_epoch(2200.0),
    );

    for (asc, dec) in stars {
        let (a1, d1) = precess::precess_eq_coords(asc, dec, j1, j2);
        let (a2, d2) = precess::precess_eq_coords(a1, d1, j2, j3);
        let (a3, d3) = precess::precess_eq_coords(asc, dec, j1, j3);
        assert!(sep_asec(a2, d2, a3, d3) < 0.005, "{asc} {dec}");
    }
}

#[test]
fn precess_eq_coords_round_trip() {
    let JD = precess::jd_frm_julian_epoch(1600.0);
    for (asc, dec) in [
        (0.0, 0.0),
        (1.0, 1.5),
        (3.0, -1.0),
        (5.9, 0.2),
        (2.0, -1.55),
    ] {
        let (a1, d1) = precess::precess_eq_coords(asc, dec, precess::J2000, JD);
        let (a2, d2) = precess::precess_eq_coords(a1, d1, JD, precess::J2000);
        assert!(sep_asec(asc, dec, a2, d2) < 0.001, "{asc} {dec}");
    }
}

#[test]
fn precess_eq_coords_identity_and_range() {
    let (a, d) = precess::precess_eq_coords(1.234, -0.567, 2440000.5, 2440000.5);
    assert_close(a, 1.234, 1e-14);
    assert_close(d, -0.567, 1e-14);

    // Right ascension stays in [0, 2π) across 0h, both directions in time
    for JD in [precess::J2000 - 36525.0, precess::J2000 + 36525.0] {
        for asc in [0.0, 1e-6, std::f64::consts::TAU - 1e-6] {
            let (a, d) = precess::precess_eq_coords(asc, 0.1, precess::J2000, JD);
            assert!((0.0..std::f64::consts::TAU).contains(&a));
            assert!(d.is_finite());
        }
    }

    // The celestial pole itself
    let (_, d) = precess::precess_eq_coords(
        0.0,
        std::f64::consts::FRAC_PI_2,
        precess::J2000,
        precess::B1900,
    );
    assert!(d > 89.0_f64.to_radians() && d < std::f64::consts::FRAC_PI_2);
}

// Precession of ecliptical coordinates and of equatorial coordinates
// describe the same rotation (both IAU 1976): going through the ecliptic
// of J2000.0 and of the date must agree with the direct equatorial
// reduction, to the accuracy of the truncated polynomials (which drift
// apart slowly: 0.005" after 500 years, 0.06" after 1000 years).
#[test]
fn ecl_and_eq_precession_agree() {
    for (year, tol) in [
        (1800.0, 0.001),
        (1500.0, 0.01),
        (1000.0, 0.1),
        (2500.0, 0.01),
    ] {
        let JD = precess::jd_frm_julian_epoch(year);
        let eps0 = ecliptic::mn_oblq_IAU(precess::J2000);
        let eps = ecliptic::mn_oblq_IAU(JD);

        for (asc, dec) in [(0.5, 0.3), (2.5, -0.8), (4.0, 1.1), (5.5, -0.1)] {
            let (a_eq, d_eq) = precess::precess_eq_coords(asc, dec, precess::J2000, JD);

            let (l0, b0) = ecl_frm_eq(asc, dec, eps0);
            let (l, b) = precess::precess_ecl_coords(l0, b0, precess::J2000, JD);
            let (a_ecl, d_ecl) = eq_frm_ecl(l, b, eps);

            assert!(
                sep_asec(a_eq, d_eq, a_ecl, d_ecl) < tol,
                "{year} {asc} {dec}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Ecliptical coordinates (eqs. 21.5 to 21.7)
// ---------------------------------------------------------------------------

// Meeus, example 21.c: Venus at -214 June 30.0 TD, J2000.0 -> mean
// equinox of the date
#[test]
fn example_21_c() {
    let JD = 1643074.5;
    assert_close((JD - precess::J2000) / 36525.0, -22.134716, 1e-6);

    let a = precess::ecl_precess_angles(precess::J2000, JD);
    assert_close(to_asec(a.eta), -1057.225, 0.0005);
    assert_close(a.eta.to_degrees(), -0.293673, 1e-6);
    assert_close_deg(a.pi.to_degrees(), 180.22924, 1e-5);
    assert_close(to_asec(a.p), -110773.167, 0.001);
    assert_close(a.p.to_degrees(), -30.770324, 1e-6);

    let (long, lat) = precess::precess_ecl_coords(
        149.48194_f64.to_radians(),
        1.76549_f64.to_radians(),
        precess::J2000,
        JD,
    );
    // p + Π - λ = 30.75475°, C' = sin β = +0.0281891
    assert_close_deg((a.p + a.pi - long).to_degrees(), 30.75475, 1e-5);
    assert_close(long.to_degrees(), 118.704, 0.0005);
    assert_close(lat.sin(), 0.0281891, 1e-7);
    assert_close(lat.to_degrees(), 1.615, 0.0005);
}

#[test]
fn precess_ecl_coords_round_trip_and_range() {
    let JD = precess::jd_frm_julian_epoch(-1000.0);
    for (l, b) in [(0.0, 0.0), (6.2, 0.1), (3.0, -1.5), (1.0, 1.56)] {
        let (l1, b1) = precess::precess_ecl_coords(l, b, precess::J2000, JD);
        assert!((0.0..std::f64::consts::TAU).contains(&l1));
        let (l2, b2) = precess::precess_ecl_coords(l1, b1, JD, precess::J2000);
        assert!(sep_asec(l, b, l2, b2) < 0.01, "{l} {b}");
    }

    // A point near the ecliptic pole keeps a finite, accurate latitude
    let (_, b) = precess::precess_ecl_coords(0.0, std::f64::consts::FRAC_PI_2, precess::J2000, JD);
    assert!(b.is_finite() && b < std::f64::consts::FRAC_PI_2);
}

// ---------------------------------------------------------------------------
// Proper motion in ecliptical coordinates
// ---------------------------------------------------------------------------

// The formulae must agree with a numerical derivative of the equatorial
// to ecliptical transformation.
#[test]
fn proper_motion_in_ecl_coords_matches_derivative() {
    let eps = 23.4392911_f64.to_radians();
    let h = 1e-3; // years
    for (asc, dec, pma, pmd) in [
        (2.65, 0.21, -1.2e-6, 3.0e-8),
        (0.7, 0.86, 8.0e-7, -4.0e-7),
        (4.2, -0.6, 2.0e-6, 1.0e-6),
        (5.9, -1.2, -3.0e-6, 5.0e-7),
    ] {
        let (pml, pmb) = precess::proper_motion_in_ecl_coords(asc, dec, pma, pmd, eps);

        let (l1, b1) = ecl_frm_eq(asc - pma * h, dec - pmd * h, eps);
        let (l2, b2) = ecl_frm_eq(asc + pma * h, dec + pmd * h, eps);
        let dl = ((l2 - l1 + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
            - std::f64::consts::PI)
            / (2.0 * h);
        let db = (b2 - b1) / (2.0 * h);

        assert_close(pml, dl, 1e-12);
        assert_close(pmb, db, 1e-12);

        // The deprecated, misnamed function computes the same thing
        let (_, lat) = ecl_frm_eq(asc, dec, eps);
        #[allow(deprecated)]
        let (pml2, pmb2) = star::proper_motion_in_eq_coords(asc, dec, pma, pmd, lat, eps);
        assert_close(pml2, pml, 1e-15);
        assert_close(pmb2, pmb, 1e-15);
    }
}

// Regulus (proper motion of example 21.a): Meeus' Table 21.A lists
// -23.48" and -8.13" per century in longitude and latitude.
#[test]
fn proper_motion_in_ecl_coords_regulus() {
    let (pml, pmb) = precess::proper_motion_in_ecl_coords(
        hms(10, 8, 22.3).to_radians(),
        dms(false, 11, 58, 2.0).to_radians(),
        hms(0, 0, -0.0169).to_radians(),
        arcsec(0.006).to_radians(),
        ecliptic::mn_oblq_IAU(precess::J2000),
    );
    assert_close(to_asec(pml) * 100.0, -23.48, 0.1);
    assert_close(to_asec(pmb) * 100.0, -8.13, 0.1);
}

// ---------------------------------------------------------------------------
// FK4 / Newcomb
// ---------------------------------------------------------------------------

// Meeus gives Newcomb's angles reduced for a starting epoch of B1950.0
// (T = 0.5): the general expressions must reproduce them.
#[test]
fn newcomb_angles_from_b1950() {
    let t = 1.0; // one tropical century
    let a = precess::eq_precess_angles_newcomb(precess::B1950, precess::B1950 + 36524.2199 * t);

    assert_close(to_asec(a.zeta), 2304.948 + 0.302 + 0.018, 0.0006);
    assert_close(to_asec(a.z), 2304.948 + 1.093 + 0.019, 0.0006);
    assert_close(to_asec(a.theta), 2004.255 - 0.426 - 0.042, 0.0006);
}

// Before the fix, precess_eq_coords_FK5 used the angles in degrees as if
// they were radians: one year of precession moved a star by degrees.
#[test]
fn newcomb_precession_is_small_and_close_to_iau() {
    for (asc, dec) in [(0.5, 0.3), (2.5, -0.8), (4.0, 1.1), (5.5, -0.1)] {
        let (a1, d1) =
            precess::precess_eq_coords_newcomb(asc, dec, precess::B1950, precess::B1950 + 365.2422);
        assert!(sep_asec(asc, dec, a1, d1) < 60.0);

        #[allow(deprecated)]
        let (a2, d2) =
            precess::precess_eq_coords_FK5(asc, dec, precess::B1950, precess::B1950 + 365.2422);
        assert_eq!((a1, d1), (a2, d2));

        // Newcomb's and the IAU 1976 precession, B1950.0 -> J2000.0, differ
        // by much less than a second of arc in the angles; here they agree
        // within 1".
        let (an, dn) = precess::precess_eq_coords_newcomb(asc, dec, precess::B1950, precess::J2000);
        let (ai, di) = precess::precess_eq_coords(asc, dec, precess::B1950, precess::J2000);
        assert!(sep_asec(an, dn, ai, di) < 1.0, "{asc} {dec}");
        assert!(sep_asec(asc, dec, an, dn) > 600.0);
    }
}

#[test]
fn newcomb_round_trip() {
    for (asc, dec) in [(0.1, 0.2), (3.3, -1.0), (6.2, 1.4)] {
        let (a1, d1) = precess::precess_eq_coords_newcomb(asc, dec, precess::B1950, precess::B1900);
        let (a2, d2) = precess::precess_eq_coords_newcomb(a1, d1, precess::B1900, precess::B1950);
        assert!(sep_asec(asc, dec, a2, d2) < 0.001);
    }
}

#[test]
fn fk4_equinox_correction() {
    // 0.0775 s at J2000.0, 0.035 s at 1950
    assert_close(
        to_sec_time(precess::fk4_equinox_correction(precess::J2000)),
        0.0775,
        1e-12,
    );
    assert_close(
        to_sec_time(precess::fk4_equinox_correction(
            precess::J2000 - 36525.0 / 2.0,
        )),
        0.035,
        1e-12,
    );
}

// ---------------------------------------------------------------------------
// Epochs
// ---------------------------------------------------------------------------

#[test]
fn julian_and_besselian_epochs() {
    // J1986.0 and J2050.0 (Meeus, chapter 21)
    assert_close(precess::jd_frm_julian_epoch(1986.0), 2446431.50, 1e-9);
    assert_close(precess::jd_frm_julian_epoch(2050.0), 2469807.50, 1e-9);
    assert_close(precess::julian_epoch_frm_jd(2469807.5), 2050.0, 1e-12);

    // B1900.0 and B1950.0
    assert_close(
        precess::jd_frm_besselian_epoch(1900.0),
        precess::B1900,
        1e-4,
    );
    assert_close(
        precess::jd_frm_besselian_epoch(1950.0),
        precess::B1950,
        1e-4,
    );
    // Example 24.a: the epoch 1744.0, 206 tropical years before B1950.0
    assert_close(precess::jd_frm_besselian_epoch(1744.0), 2358042.5305, 1e-4);

    for e in [-4000.0, 1582.5, 1950.0, 2400.25] {
        assert_close(
            precess::besselian_epoch_frm_jd(precess::jd_frm_besselian_epoch(e)),
            e,
            1e-9,
        );
        assert_close(
            precess::julian_epoch_frm_jd(precess::jd_frm_julian_epoch(e)),
            e,
            1e-9,
        );
    }
}

// ---------------------------------------------------------------------------
// Motion in space (example 21.d) and the pole angle
// ---------------------------------------------------------------------------

// Meeus, example 21.d: Sirius, J2000.0, distance 2.64 pc, radial velocity
// -7.6 km/s, proper motion -0.03847 s and -1.2053" per year.
#[test]
fn example_21_d() {
    let asc0 = hms(6, 45, 8.871);
    let dec0 = dms(true, 16, 42, 57.99);
    assert_close(asc0, 101.286962, 1e-6);
    assert_close(dec0, -16.716108, 1e-6);

    let delta_r = -7.6 / 977792.0;
    assert_close(delta_r, -0.000007773, 1e-9);
    let pm_asc = hms(0, 0, -0.03847).to_radians();
    let pm_dec = arcsec(-1.2053).to_radians();
    assert_close(pm_asc, -0.0000027976, 1e-10);
    assert_close(pm_dec, -0.0000058435, 1e-10);

    let cases = [
        (-1000.0, (6, 45, 47.16), (16, 22, 56.0)),
        (-2000.0, (6, 46, 25.09), (16, 3, 0.8)),
        (-3000.0, (6, 47, 2.67), (15, 43, 12.3)),
        (-4000.0, (6, 47, 39.91), (15, 23, 30.6)),
        (-12000.0, (6, 52, 25.72), (12, 50, 6.7)),
    ];

    for (t, (h, m, s), (d, dm, ds)) in cases {
        let (asc, dec) = star::eq_coords_frm_motion(
            asc0.to_radians(),
            dec0.to_radians(),
            2.64,
            delta_r,
            pm_asc,
            pm_dec,
            t,
        );
        let (hh, mm, ss) = hms_of(asc);
        assert_eq!((hh, mm), (h, m), "t = {t}");
        assert_close(ss, s, 0.01);
        assert_close(
            dec.to_degrees() * 3600.0,
            dms(true, d, dm, ds) * 3600.0,
            0.1,
        );
    }

    // With -7.7 km/s the declination at -10000 becomes -12°50'13.0"
    let (_, dec) = star::eq_coords_frm_motion(
        asc0.to_radians(),
        dec0.to_radians(),
        2.64,
        -7.7 / 977792.0,
        pm_asc,
        pm_dec,
        -12000.0,
    );
    assert_close(
        dec.to_degrees() * 3600.0,
        dms(true, 12, 50, 13.0) * 3600.0,
        0.1,
    );
}

#[test]
fn eq_coords_frm_motion_limits() {
    // No motion at all, or t = 0: the position is unchanged
    let (a, d) = star::eq_coords_frm_motion(6.2, -0.3, 10.0, 0.0, 0.0, 0.0, 5000.0);
    assert_close(a, 6.2, 1e-12);
    assert_close(d, -0.3, 1e-12);
    let (a, d) = star::eq_coords_frm_motion(6.2, -0.3, 10.0, 1e-5, 1e-6, 1e-6, 0.0);
    assert_close(a, 6.2, 1e-12);
    assert_close(d, -0.3, 1e-12);

    // A radial velocity alone does not move the star on the sky
    let (a, d) = star::eq_coords_frm_motion(1.0, 0.5, 10.0, 1e-5, 0.0, 0.0, 1e5);
    assert_close(a, 1.0, 1e-12);
    assert_close(d, 0.5, 1e-12);

    // For a short time it matches a uniform proper motion
    let (pma, pmd) = (2e-6, -1e-6);
    let (a, d) = star::eq_coords_frm_motion(0.0, 0.2, 5.0, 0.0, pma, pmd, 10.0);
    assert_close(a, pma * 10.0, 1e-9);
    assert_close(d, 0.2 + pmd * 10.0, 1e-9);

    // Right ascension stays in [0, 2π) when crossing 0h westwards
    let (a, _) = star::eq_coords_frm_motion(1e-6, 0.2, 5.0, 0.0, -2e-6, 0.0, 10.0);
    assert!((0.0..std::f64::consts::TAU).contains(&a) && a > 6.2);

    // A star at the Sun's position (r = 0) does not produce NaN
    let (a, d) = star::eq_coords_frm_motion(1.0, 0.5, 0.0, 1e-5, 0.0, 0.0, 10.0);
    assert!(a.is_finite() && d.is_finite());
}

// Meeus, chapter 14: the angle q between the directions of the north
// celestial pole and of the north ecliptic pole at a star equals the
// position angle of the ecliptic pole seen from the star.
#[test]
fn angl_between_north_celes_and_eclip_pole() {
    let eps = 23.4392911_f64.to_radians();
    let (pole_a, pole_d) = (270.0_f64.to_radians(), std::f64::consts::FRAC_PI_2 - eps);

    for (l, b) in [
        (0.0, 0.0),
        (1.0, 0.3),
        (2.5, -0.7),
        (4.0, 1.2),
        (5.5, -1.3),
        (1.57, 0.0),
    ] {
        let q = star::angl_between_north_celes_and_eclip_pole(l, b, eps);
        assert!(q > -std::f64::consts::PI && q <= std::f64::consts::PI);

        let (a, d) = eq_frm_ecl(l, b, eps);
        let pa = ((pole_a - a).sin() * pole_d.cos())
            .atan2(d.cos() * pole_d.sin() - d.sin() * pole_d.cos() * (pole_a - a).cos());
        assert_close(q, pa, 1e-12);
    }

    // On the ecliptic, tan q0 = -cos λ tan ε; at the vernal equinox the
    // ecliptic pole lies ε to the west of north, at λ = 90° straight north.
    let q0 = star::angl_between_north_celes_and_eclip_pole(0.0, 0.0, eps);
    assert_close(q0, -eps, 1e-12);
    let q90 = star::angl_between_north_celes_and_eclip_pole(std::f64::consts::FRAC_PI_2, 0.0, eps);
    assert_close(q90, 0.0, 1e-12);
}

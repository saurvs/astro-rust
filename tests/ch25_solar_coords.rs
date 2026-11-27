//! Meeus, chapter 25: solar coordinates.

mod common;

use astro::*;
use common::*;

use std::f64::consts::{PI, TAU};

/// 1992 October 13.0 TD, the date of examples 25.a and 25.b.
const JDE_1992_OCT_13: f64 = 2448908.5;

fn deg(rad: f64) -> f64 {
    rad.to_degrees()
}

fn wrap_pi(a: f64) -> f64 {
    (a + PI).rem_euclid(TAU) - PI
}

// Meeus, example 25.a: low accuracy method, every intermediate value
#[test]
fn example_25a_low_accuracy() {
    let jde = JDE_1992_OCT_13;
    assert_close(time::julian_cent(jde), -0.072183436, 1e-9);

    let p = sun::low_accuracy_position(jde);
    assert_close_deg(deg(p.mean_longitude), 201.80720, 5e-6);
    assert_close_deg(deg(p.mean_anomaly), 278.99397, 5e-6);
    assert_close(p.eccentricity, 0.016711668, 5e-10);
    assert_close(deg(p.equation_of_center), -1.89732, 5e-6);
    // the book adds the rounded L0 and C, hence the looser tolerance
    assert_close_deg(deg(p.true_longitude), 199.90988, 1e-5);
    assert_close(p.radius_vector, 0.99766, 5e-6);
    assert_close_deg(deg(p.apparent_longitude), 199.90895, 1e-5);

    // the free functions agree with the struct
    assert_eq!(sun::geometric_mean_longitude(jde), p.mean_longitude);
    assert_eq!(sun::mean_anomaly(jde), p.mean_anomaly);
    assert_eq!(sun::earth_orbit_eccentricity(jde), p.eccentricity);
    assert_eq!(sun::equation_of_center(jde), p.equation_of_center);

    // ε0 = 23°26'24.83" by eq. 22.2, corrected by +0.00256° cos Ω, Ω = 264.65°
    assert_close_deg(
        deg(ecliptic::mn_oblq_IAU(jde)),
        dms(false, 23, 26, 24.83),
        arcsec(0.01),
    );

    let eq = sun::low_accuracy_apparent_eq_point(jde);
    assert_close_deg(deg(eq.asc), 198.38083, 5e-6);
    assert_close(deg(eq.dec), -7.78507, 5e-6);
    // 13h13m31.4s
    assert_close(deg(eq.asc), hms(13, 13, 31.4), 0.05 * 15.0 / 3600.0);
}

// Meeus, example 25.b: the conversion to the FK5 system, with the book's
// input values from the abridged VSOP87 of Appendix III
#[test]
fn example_25b_fk5_correction() {
    let (long, lat) = sun::ecl_coords_to_FK5(
        JDE_1992_OCT_13,
        199.907372_f64.to_radians(),
        arcsec(0.644).to_radians(),
    );

    // ΔΘ = -0.09033", Δβ = -0.023"
    assert_close_deg(deg(long), 199.907347, 5e-7);
    assert_close(deg(long) * 3600.0, 199.907372 * 3600.0 - 0.09033, 1e-6);
    assert_close(deg(lat) * 3600.0, 0.644 - 0.023, 0.0005);
}

// The FK5 longitude is reduced to [0, 2π) even just past 0°.
#[test]
fn fk5_correction_wraps_longitude() {
    let (long, _) = sun::ecl_coords_to_FK5(JDE_1992_OCT_13, 1e-9, 0.0);
    assert!((0.0..TAU).contains(&long));
    assert_close_rad(long, 1e-9 - arcsec(0.09033).to_radians(), 1e-15);
}

// Meeus, example 25.b: nutation, aberration (25.10), apparent longitude and
// apparent right ascension and declination, each step with the book's
// input values
#[test]
fn example_25b_steps() {
    let jde = JDE_1992_OCT_13;

    let (dpsi, deps) = nutation::nutation(jde);
    assert_close(deg(dpsi) * 3600.0, 15.908, 0.0005);
    assert_close(deg(deps) * 3600.0, -0.308, 0.0005);

    // aberration by eq. 25.10, with R = 0.99760775
    let ab = aberr::sol_aberr(0.99760775);
    assert_close(deg(ab) * 3600.0, -20.539, 0.0005);

    // λ = 199.907347° + 15.908" - 20.539" = 199°54'21.818"
    let lambda = 199.907347 + arcsec(15.908) - arcsec(20.539);
    assert_close_deg(lambda, dms(false, 199, 54, 21.818), arcsec(0.002));

    // α, δ by eqs. 13.3 and 13.4, β = +0.62", true ε = 23.4401443°
    let beta = arcsec(0.62).to_radians();
    let eps = 23.4401443_f64.to_radians();
    let asc = coords::asc_frm_ecl(lambda.to_radians(), beta, eps);
    let dec = coords::dec_frm_ecl(lambda.to_radians(), beta, eps);
    assert_close_deg(deg(asc), 198.378178, 5e-7);
    assert_close(deg(dec), -7.783871, 5e-7);
}

// Meeus, example 25.b (geometric part) with the complete VSOP87 of this
// crate: the book's "correct values" quoted at the end of example 25.a.
#[test]
fn geocent_ecl_pos_complete_vsop87() {
    let (p, r) = sun::geocent_ecl_pos(JDE_1992_OCT_13);

    // dynamical frame: 0.09033" more than the FK5 value 199°54'26.18"
    assert_close_deg(
        deg(p.long),
        dms(false, 199, 54, 26.18 + 0.09033),
        arcsec(0.006),
    );
    assert_close(p.lat.to_degrees() * 3600.0, 0.744, 0.001);
    assert_close(r, 0.99760853, 1.5e-8);

    let (q, r2) = sun::geocent_ecl_pos_fk5(JDE_1992_OCT_13);
    assert_eq!(r, r2);
    assert_close_deg(deg(q.long), dms(false, 199, 54, 26.18), arcsec(0.006));
    assert_close(q.lat.to_degrees() * 3600.0, 0.72, 0.005);
}

// Meeus, chapter 25: the apparent place with the complete VSOP87, compared
// with the book's correct values (end of example 25.a).
#[test]
fn apparent_position_complete_vsop87() {
    let p = sun::apparent_position(JDE_1992_OCT_13);

    assert_close_deg(
        deg(p.geometric_longitude),
        dms(false, 199, 54, 26.18),
        arcsec(0.006),
    );
    assert_close_deg(
        deg(p.apparent_longitude),
        dms(false, 199, 54, 21.56),
        arcsec(0.006),
    );
    assert_close(deg(p.latitude) * 3600.0, 0.72, 0.005);
    assert_close(p.radius_vector, 0.99760853, 1.5e-8);
    assert_close_deg(
        deg(p.right_ascension),
        hms(13, 13, 30.749),
        0.0006 * 15.0 / 3600.0,
    );
    assert_close(deg(p.declination), dms(true, 7, 47, 1.74), arcsec(0.006));

    // nutation and obliquity as in example 25.b
    assert_close(deg(p.nutation_in_longitude) * 3600.0, 15.908, 0.0005);
    assert_close(deg(p.nutation_in_obliquity) * 3600.0, -0.308, 0.0005);
    assert_close(deg(p.true_obliquity), 23.4401443, 5e-7);
    assert_close(
        p.true_obliquity,
        p.mean_obliquity + p.nutation_in_obliquity,
        1e-15,
    );

    // λ = ☉ + Δψ + aberration
    assert_close_rad(
        p.apparent_longitude,
        p.geometric_longitude + p.nutation_in_longitude + p.aberration,
        1e-15,
    );

    let eq = sun::apparent_eq_point(JDE_1992_OCT_13);
    assert_eq!(eq.asc, p.right_ascension);
    assert_eq!(eq.dec, p.declination);
}

// Bug: geocent_ecl_pos returned a latitude of about 2π instead of a small
// negative angle whenever the Sun was south of the ecliptic, because the
// latitude was reduced with limit_to_two_PI. On 2000 January 15.5 TD the
// Sun's latitude is about -0.48".
#[test]
fn geocent_ecl_pos_negative_latitude() {
    let (p, _) = sun::geocent_ecl_pos(2451559.0);

    assert!(p.lat < 0.0, "latitude {} should be negative", p.lat);
    assert_close(deg(p.lat) * 3600.0, -0.48, 0.01);

    // the Earth's heliocentric latitude is the opposite
    let (_, b, _) = planet::heliocent_coords(&planet::Planet::Earth, 2451559.0);
    assert_close(p.lat, -wrap_pi(b), 1e-15);
}

// The Sun's latitude never exceeds 1.2" with respect to the ecliptic of the
// date (Meeus, chapter 25); the longitude is in [0, 2π) and moves smoothly.
#[test]
fn latitude_bound_and_longitude_range() {
    let mut prev = sun::geocent_ecl_pos(2415020.5).0.long;
    // 1900 to 2100 every 11 days, plus some ancient and future dates
    for i in 1..6700 {
        let jde = 2415020.5 + 11.0 * i as f64;
        let (p, r) = sun::geocent_ecl_pos(jde);
        assert!(
            p.lat.abs() < arcsec(1.2).to_radians(),
            "lat {} at {}",
            p.lat,
            jde
        );
        assert!((0.0..TAU).contains(&p.long));
        assert!((0.983..1.017).contains(&r));

        // 11 days of motion: between 10.4° and 11.3°
        let step = deg((p.long - prev).rem_euclid(TAU));
        assert!((10.4..11.3).contains(&step), "step {step} at {jde}");
        prev = p.long;
    }
}

// The low accuracy method agrees with the complete theory to about 0.01°
// (Meeus, chapter 25), and degrades slowly away from J2000.
#[test]
fn low_accuracy_vs_vsop87() {
    for (year, tol) in [
        (-2000, 0.03),
        (-500, 0.012),
        (1000, 0.01),
        (1900, 0.01),
        (2000, 0.01),
        (2100, 0.01),
        (3000, 0.01),
        (4000, 0.012),
    ] {
        let jd0 = 2451545.0 + (year as f64 - 2000.0) * 365.25;
        for i in 0..37 {
            let jde = jd0 + 10.0 * i as f64;
            let lo = sun::low_accuracy_position(jde);
            let hi = sun::apparent_position(jde);
            assert_close_deg(deg(lo.apparent_longitude), deg(hi.apparent_longitude), tol);
            assert_close(lo.radius_vector, hi.radius_vector, 1e-4);

            let eq_lo = sun::low_accuracy_apparent_eq_point(jde);
            assert_close_deg(deg(eq_lo.asc), deg(hi.right_ascension), 1.2 * tol);
            assert_close(deg(eq_lo.dec), deg(hi.declination), tol);
        }
    }
}

// Low accuracy internals: v = M + C, ☉ = L0 + C and eq. 25.5 at perihelion
// and aphelion.
#[test]
fn low_accuracy_relations() {
    for i in 0..50 {
        let jde = 2451545.0 + 97.3 * i as f64 - 2000.0;
        let p = sun::low_accuracy_position(jde);
        assert_close_rad(
            p.true_longitude,
            p.mean_longitude + p.equation_of_center,
            1e-12,
        );
        assert_close_rad(p.true_anomaly, p.mean_anomaly + p.equation_of_center, 1e-12);
        assert!(p.equation_of_center.abs() < 1.93_f64.to_radians());
        let e = p.eccentricity;
        let q = 1.000001018 * (1.0 - e);
        let big_q = 1.000001018 * (1.0 + e);
        assert!(p.radius_vector >= q - 1e-12 && p.radius_vector <= big_q + 1e-12);
    }
}

// Meeus, chapter 25: the J2000.0 low accuracy conversion agrees with the
// rigorous precession to about 0.01° between 1900 and 2100.
#[test]
fn low_accuracy_long_to_j2000() {
    for year in [1900.3, 1950.0, 1992.78, 2000.0, 2050.6, 2099.9] {
        let jde = 2451545.0 + (year - 2000.0) * 365.25;
        let (date, _) = sun::geocent_ecl_pos(jde);
        let (j2000, _) = sun::geocent_ecl_pos_j2000(jde);
        let approx = sun::low_accuracy_long_to_j2000(date.long, year);
        assert!((0.0..TAU).contains(&approx));
        assert_close_deg(deg(approx), deg(j2000.long), 0.001);
    }
    // wrap-around below 0°
    let l = sun::low_accuracy_long_to_j2000(0.0, 2010.0);
    assert_close(deg(l), 360.0 - 0.1397, 1e-9);
}

// Eq. 25.11 (daily variation) and eq. 25.10 differ by at most about 0.01".
#[test]
fn aberration_daily_variation_vs_constant() {
    for i in 0..120 {
        let jde = 2448908.5 + 3.1 * i as f64;
        let r = sun::geocent_ecl_pos(jde).1;
        let accurate = sun::aberration_in_longitude(jde, r);
        let simple = aberr::sol_aberr(r);
        assert!((deg(accurate - simple) * 3600.0).abs() < 0.02);
        assert!((-20.9..-20.1).contains(&(deg(accurate) * 3600.0)));
    }
}

// The daily variation in longitude agrees with the numerical derivative of
// the VSOP87 longitude: in the fixed J2000.0 frame, and of the date.
#[test]
fn daily_variation_matches_numerical_derivative() {
    for i in 0..40 {
        let jde = 2451545.0 + 91.7 * i as f64 - 1800.0;
        let h = 0.5;

        let d_fixed = (sun::geocent_ecl_pos_j2000(jde + h).0.long
            - sun::geocent_ecl_pos_j2000(jde - h).0.long)
            .rem_euclid(TAU)
            / (2.0 * h);
        let d_date = (sun::geocent_ecl_pos(jde + h).0.long - sun::geocent_ecl_pos(jde - h).0.long)
            .rem_euclid(TAU)
            / (2.0 * h);

        // Meeus: the error of the truncated series is below 0.1" per day
        let fixed = sun::daily_variation_in_longitude(jde);
        let of_date = sun::daily_variation_in_longitude_of_date(jde);
        assert_close(deg(fixed) * 3600.0, deg(d_fixed) * 3600.0, 0.1);
        assert_close(deg(of_date) * 3600.0, deg(d_date) * 3600.0, 0.1);
        // the two frames differ by the general precession, 0.137" per day
        assert_close(deg(of_date - fixed) * 3600.0, 0.137, 1e-9);
    }
}

// Right ascension in [0, 2π), always within 2.5° of the longitude, and on
// the same side of the equinox.
#[test]
fn right_ascension_quadrant() {
    for i in 0..400 {
        let jde = 2451545.0 + 0.913 * i as f64;
        let p = sun::apparent_position(jde);
        assert!((0.0..TAU).contains(&p.right_ascension));
        assert!((0.0..TAU).contains(&p.apparent_longitude));
        assert!(wrap_pi(p.right_ascension - p.apparent_longitude).abs() < 2.5_f64.to_radians());
        assert!(p.declination.abs() <= p.true_obliquity + 1e-6);
    }
}

// Dates before the Gregorian reform and in the far past or future still
// give finite, normalised results.
#[test]
fn ancient_and_future_dates() {
    // -1000 January 1 (Julian), 1000, 1582 October 4, 3000, 6000
    for jde in [1355807.5, 2086307.5, 2299159.5, 2816787.5, 3912521.5] {
        let p = sun::apparent_position(jde);
        assert!((0.0..TAU).contains(&p.apparent_longitude));
        assert!((0.0..TAU).contains(&p.right_ascension));
        assert!(p.latitude.abs() < arcsec(1.5).to_radians());
        assert!((0.98..1.02).contains(&p.radius_vector));

        let lo = sun::low_accuracy_position(jde);
        assert!((0.0..TAU).contains(&lo.apparent_longitude));
    }
}

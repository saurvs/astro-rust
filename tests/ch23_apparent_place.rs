// Meeus, Astronomical Algorithms, 2nd ed., chapter 23: apparent place of
// a star (aberration and the full reduction; nutation in equatorial
// coordinates is tested with chapter 22).

#![allow(non_snake_case)]

mod common;

use astro::*;
use common::{assert_close, assert_close_deg, dms, hms};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

fn to_asec(x: f64) -> f64 {
    x.to_degrees() * 3600.0
}

/// Right ascension in radians to (h, m, s)
fn hms_of(asc: f64) -> (i64, i64, f64) {
    let s = asc.to_degrees().rem_euclid(360.0) * 240.0;
    let h = (s / 3600.0).floor();
    let m = ((s - h * 3600.0) / 60.0).floor();
    (h as i64, m as i64, s - h * 3600.0 - m * 60.0)
}

/// Angular separation of two points, in arcseconds
fn sep_asec(a1: f64, d1: f64, a2: f64, d2: f64) -> f64 {
    let c = d1.sin() * d2.sin() + d1.cos() * d2.cos() * (a1 - a2).cos();
    let x = (d1.cos() * d2.sin() - d1.sin() * d2.cos() * (a2 - a1).cos())
        .hypot(d2.cos() * (a2 - a1).sin());
    to_asec(x.atan2(c))
}

fn ecl_frm_eq(asc: f64, dec: f64, eps: f64) -> (f64, f64) {
    let long = (asc.sin() * eps.cos() + dec.tan() * eps.sin()).atan2(asc.cos());
    let lat = (dec.sin() * eps.cos() - dec.cos() * eps.sin() * asc.sin()).asin();
    (long.rem_euclid(TAU), lat)
}

fn eq_frm_ecl(long: f64, lat: f64, eps: f64) -> (f64, f64) {
    let asc = (long.sin() * eps.cos() - lat.tan() * eps.sin()).atan2(long.cos());
    let dec = (lat.sin() * eps.cos() + lat.cos() * eps.sin() * long.sin()).asin();
    (asc.rem_euclid(TAU), dec)
}

/// The Sun's geometric longitude, low accuracy (Meeus, chapter 25)
fn sun_long_low_accuracy(JD: f64) -> f64 {
    let t = time::julian_cent(JD);
    let l0 = 280.46646 + t * (36000.76983 + 0.0003032 * t);
    let m = (357.52911 + t * (35999.05029 - 0.0001537 * t)).to_radians();
    let c = (1.914602 - t * (0.004817 + 0.000014 * t)) * m.sin()
        + (0.019993 - 0.000101 * t) * (2.0 * m).sin()
        + 0.000289 * (3.0 * m).sin();
    (l0 + c).to_radians()
}

/// θ Persei, J2000.0 (Meeus, examples 21.b, 23.a, 23.b)
fn theta_persei() -> star::MeanPlace {
    star::MeanPlace {
        asc: hms(2, 44, 11.986).to_radians(),
        dec: dms(false, 49, 13, 42.48).to_radians(),
        pm_asc: hms(0, 0, 0.03425).to_radians(),
        pm_dec: (-0.0895 / 3600.0_f64).to_radians(),
        epoch: precess::J2000,
    }
}

const JD_2028: f64 = 2462088.69; // 2028 November 13.19 TD

// ---------------------------------------------------------------------------
// Classical aberration (eqs. 23.2, 23.3)
// ---------------------------------------------------------------------------

// Meeus, example 23.a: e, π, and the aberration of θ Persei
#[test]
fn example_23_a_aberration() {
    assert_close(time::julian_cent(JD_2028), 0.2886705, 1e-7);
    assert_close(aberr::earth_orb_ecc(JD_2028), 0.01669649, 1e-8);
    assert_close(
        aberr::earth_perih_long(JD_2028).to_degrees(),
        103.434,
        0.0005,
    );
    assert_close(to_asec(aberr::CONST_OF_ABERR), 20.49552, 1e-10);
    // ☉ = 231.328° by the low accuracy method of chapter 25, which
    // apparent_place uses
    assert_close(
        sun::low_accuracy_position(JD_2028)
            .true_longitude
            .to_degrees(),
        231.328,
        0.0005,
    );

    let (da, dd) = aberr::annual_aberr_in_eq_coords(
        41.5472_f64.to_radians(),
        49.3485_f64.to_radians(),
        231.328_f64.to_radians(),
        23.436_f64.to_radians(),
        JD_2028,
    );
    assert_close(to_asec(da), 30.045, 0.001);
    assert_close(to_asec(dd), 6.697, 0.001);
}

// The ecliptical (23.2) and equatorial (23.3) forms describe the same
// displacement, to second order in κ.
#[test]
fn ecl_and_eq_aberration_agree() {
    let eps = 23.436_f64.to_radians();
    for sun in [0.0, 1.5, 4.0372] {
        for (a, d) in [
            (0.7, 0.86),
            (2.0, -0.5),
            (3.5, 0.0),
            (5.0, 1.0),
            (6.2, -1.0),
            (1.0, 1.5),
        ] {
            let (l, b) = ecl_frm_eq(a, d, eps);

            let (da, dd) = aberr::annual_aberr_in_eq_coords(a, d, sun, eps, JD_2028);
            let (dl, db) = aberr::annual_aberr_in_ecl_coords(l, b, sun, JD_2028);
            let (a2, d2) = eq_frm_ecl(l + dl, b + db, eps);
            assert!(sep_asec(a + da, d + dd, a2, d2) < 0.02, "{sun} {a} {d}");

            let (da, dd) = aberr::annual_aberr_in_eq_coords_fk4(a, d, sun, eps);
            let (dl, db) = aberr::annual_aberr_in_ecl_coords_fk4(l, b, sun);
            let (a2, d2) = eq_frm_ecl(l + dl, b + db, eps);
            assert!(sep_asec(a + da, d + dd, a2, d2) < 0.02, "{sun} {a} {d}");
        }
    }
}

// The E-terms are those of size eκ ≈ 0.34"
#[test]
fn e_terms() {
    let eps = 23.44_f64.to_radians();
    for (a, d) in [(0.7, 0.86), (2.0, -0.5), (3.5, 0.0), (5.0, 1.0)] {
        let (fa, fd) = aberr::annual_aberr_in_eq_coords(a, d, 2.0, eps, JD_2028);
        let (ka, kd) = aberr::annual_aberr_in_eq_coords_fk4(a, d, 2.0, eps);
        let e_term = sep_asec(a + fa, d + fd, a + ka, d + kd);
        assert!(e_term > 0.01 && e_term <= 0.35, "{e_term}");

        let (l, b) = ecl_frm_eq(a, d, eps);
        let (fl, fb) = aberr::annual_aberr_in_ecl_coords(l, b, 2.0, JD_2028);
        let (kl, kb) = aberr::annual_aberr_in_ecl_coords_fk4(l, b, 2.0);
        // On the ecliptic the E-term in longitude is eκ cos(π - λ)
        let e = aberr::earth_orb_ecc(JD_2028);
        let per = aberr::earth_perih_long(JD_2028);
        assert_close(
            (fl - kl) * b.cos(),
            e * aberr::CONST_OF_ABERR * (per - l).cos(),
            1e-15,
        );
        assert_close(
            fb - kb,
            e * aberr::CONST_OF_ABERR * b.sin() * (per - l).sin(),
            1e-15,
        );
    }
}

// A star at the ecliptic pole describes a circle of radius κ
#[test]
fn aberration_at_ecliptic_pole() {
    for sun in [0.0, 1.0, 2.0, 3.0] {
        let (dl, db) = aberr::annual_aberr_in_ecl_coords_fk4(0.3, FRAC_PI_2 - 1e-9, sun);
        assert_close(to_asec(db), -20.49552 * (sun - 0.3_f64).sin(), 1e-6);
        assert!(dl.is_finite());
    }
}

// ---------------------------------------------------------------------------
// Ron–Vondrák (Table 23.A, eq. 23.4)
// ---------------------------------------------------------------------------

// Meeus, example 23.b
#[test]
fn example_23_b_ron_vondrak() {
    let (x, y, z) = aberr::earth_bary_vel(JD_2028);
    assert_close(x * 1e8, -1363700.0, 1.0);
    assert_close(y * 1e8, 990286.0, 1.0);
    assert_close(z * 1e8, 429285.0, 1.0);

    let p = coords::EqPoint {
        asc: 41.0540613_f64.to_radians(),
        dec: 49.2277489_f64.to_radians(),
    };
    let (da, dd) = aberr::stell_aberr_in_eq_coords(&p, JD_2028);
    assert_close(da, 0.000145252, 1e-9);
    assert_close(dd, 0.000032723, 1e-9);
    assert_close(da.to_degrees(), 0.0083223, 1e-7);
    assert_close(dd.to_degrees(), 0.0018749, 1e-7);
}

// The Earth's orbital speed stays between about 29.3 and 30.3 km/s
// (0.01692 to 0.01750 AU/day)
#[test]
fn earth_bary_speed() {
    let mut JD = precess::jd_frm_julian_epoch(1995.0);
    let (mut lo, mut hi) = (1.0_f64, 0.0_f64);
    while JD < precess::jd_frm_julian_epoch(2005.0) {
        let (x, y, z) = aberr::earth_bary_vel(JD);
        let v = (x * x + y * y + z * z).sqrt();
        lo = lo.min(v);
        hi = hi.max(v);
        JD += 1.0;
    }
    assert!(lo > 0.01690 && lo < 0.01695, "{lo}");
    assert!(hi > 0.01747 && hi < 0.01752, "{hi}");
}

// Near J2000.0, where the mean frame of the date is the J2000.0 frame, the
// classical formulae and Ron–Vondrák agree within 0.03" (the classical
// ones ignore the perturbations of the Earth's motion and the Sun's own
// motion about the barycentre).
#[test]
fn ron_vondrak_vs_classical_near_j2000() {
    let mut JD = precess::J2000 - 180.0;
    while JD < precess::J2000 + 180.0 {
        let sun = sun_long_low_accuracy(JD);
        let oblq = ecliptic::mn_oblq_IAU(JD);
        for (a, d) in [
            (0.7, 0.86),
            (2.0, -0.5),
            (3.5, 0.0),
            (5.0, 1.0),
            (6.2, -1.0),
        ] {
            let (ca, cd) = aberr::annual_aberr_in_eq_coords(a, d, sun, oblq, JD);
            let (ra, rd) = aberr::stell_aberr_in_eq_coords(&coords::EqPoint { asc: a, dec: d }, JD);
            assert!(
                sep_asec(a + ca, d + cd, a + ra, d + rd) < 0.03,
                "{JD} {a} {d}"
            );
        }
        JD += 5.0;
    }
}

#[test]
fn sol_aberr() {
    assert_close(to_asec(aberr::sol_aberr(1.0)), -20.4898, 1e-10);
    assert_close(
        to_asec(aberr::sol_aberr(0.99766)),
        -20.4898 / 0.99766,
        1e-10,
    );
}

// ---------------------------------------------------------------------------
// Full reduction
// ---------------------------------------------------------------------------

// Meeus, example 23.a: classical aberration
#[test]
fn example_23_a_apparent_place() {
    let p = star::apparent_place(&theta_persei(), JD_2028, star::AberrMethod::Classical);

    // Proper motion (example 21.b)
    assert_close(p.asc_pm.to_degrees(), 41.0540613, 1e-7);
    assert_close(p.dec_pm.to_degrees(), 49.2277489, 1e-7);

    // Mean place of the date: 41.5472°, +49.3485° (2h46m11.331s,
    // +49°20'54.54" from rounded inputs in the book)
    assert_close(p.asc_mean.to_degrees(), 41.547214, 2e-6);
    assert_close(p.dec_mean.to_degrees(), 49.348483, 2e-6);

    assert_close(to_asec(p.nut_in_long), 14.861, 0.0005);
    assert_close(to_asec(p.nut_in_oblq), 2.705, 0.0005);
    assert_close(p.tru_oblq.to_degrees(), 23.436, 0.0005);
    assert_close(p.mn_oblq.to_degrees(), 23.436, 0.0005);

    // Δα1, Δδ1 (nutation): 15.843", 6.218"; Δα2, Δδ2 (aberration):
    // 30.045", 6.697". Here the nutation is applied to the position already
    // corrected for aberration, which changes Δα1 by 0.001". The aberration
    // is applied as a vector (so that it holds at the poles), not by adding
    // the first-order Δα2, Δδ2 to α and δ; the two differ by second-order
    // terms, κ² tan δ ≈ 0.0012" here.
    assert_close(to_asec(p.nut_asc), 15.843, 0.002);
    assert_close(to_asec(p.nut_dec), 6.218, 0.002);
    assert_close(to_asec(p.aberr_asc), 30.045, 0.002);
    assert_close(to_asec(p.aberr_dec), 6.697, 0.002);

    // 2h46m14.390s, +49°21'07.45"
    let (h, m, s) = hms_of(p.asc);
    assert_eq!((h, m), (2, 46));
    assert_close(s, 14.390, 0.001);
    assert_close(to_asec(p.dec), dms(false, 49, 21, 7.45) * 3600.0, 0.005);

    assert_close(p.asc, p.asc_mean + p.aberr_asc + p.nut_asc, 1e-15);
    assert_close(p.dec, p.dec_mean + p.aberr_dec + p.nut_dec, 1e-15);
}

// Meeus, example 23.b: Ron–Vondrák aberration
#[test]
fn example_23_b_apparent_place() {
    let p = star::apparent_place(&theta_persei(), JD_2028, star::AberrMethod::RonVondrak);

    assert_close(p.asc_pm.to_degrees(), 41.0540613, 1e-7);
    assert_close(p.dec_pm.to_degrees(), 49.2277489, 1e-7);

    // Book: 0.000145252, 0.000032723 rad from eq. 23.4. The aberration is
    // applied as a vector, which differs from adding these first-order
    // corrections by second-order terms of about 0.001" (5e-9 rad).
    assert_close(p.aberr_asc, 0.000145252, 1e-8);
    assert_close(p.aberr_dec, 0.000032723, 1e-8);

    // After precession: 41.5555635°, +49.3503415° (the same second-order
    // difference, 3e-7°, carries over)
    assert_close(p.asc_mean.to_degrees(), 41.5555635, 5e-7);
    assert_close(p.dec_mean.to_degrees(), 49.3503415, 5e-7);

    // Nutation: 15.844", 6.217"
    assert_close(to_asec(p.nut_asc), 15.844, 0.001);
    assert_close(to_asec(p.nut_dec), 6.217, 0.001);

    // 41.5599646° = 2h46m14.392s, 49.3520685° = +49°21'07.45"
    assert_close(p.asc.to_degrees(), 41.5599646, 5e-7);
    assert_close(p.dec.to_degrees(), 49.3520685, 5e-7);
    let (h, m, s) = hms_of(p.asc);
    assert_eq!((h, m), (2, 46));
    assert_close(s, 14.392, 0.001);
    assert_close(to_asec(p.dec), dms(false, 49, 21, 7.45) * 3600.0, 0.005);
}

// The two aberration methods give the same apparent place within 0.1"
// (0.03" in examples 23.a and 23.b), also near the pole.
#[test]
fn classical_and_ron_vondrak_agree() {
    for year in [1950.0, 1990.0, 2028.0, 2050.0] {
        for (a, d) in [
            (0.7, 0.86),
            (2.0, -0.5),
            (3.5, 0.0),
            (5.0, 1.0),
            (6.2, -1.0),
            (0.66, 1.5577),
        ] {
            let s = star::MeanPlace {
                asc: a,
                dec: d,
                pm_asc: 1e-7,
                pm_dec: -1e-7,
                epoch: precess::J2000,
            };
            for k in 0..12 {
                let JD = precess::jd_frm_julian_epoch(year) + 30.4 * k as f64;
                let c = star::apparent_place(&s, JD, star::AberrMethod::Classical);
                let r = star::apparent_place(&s, JD, star::AberrMethod::RonVondrak);
                assert!(
                    sep_asec(c.asc, c.dec, r.asc, r.dec) < 0.1,
                    "{year} {a} {d} {k}"
                );
            }
        }
    }
}

// The same star given in a catalogue for another epoch and equinox gives
// the same apparent place.
#[test]
fn catalogue_epoch_does_not_matter() {
    // Without proper motion: its components are not invariant under
    // precession (the direction of the motion rotates with the frame)
    let s2000 = star::MeanPlace {
        pm_asc: 0.0,
        pm_dec: 0.0,
        ..theta_persei()
    };
    let epoch = precess::jd_frm_julian_epoch(2050.0);
    let (asc, dec) = precess::precess_eq_coords(s2000.asc, s2000.dec, s2000.epoch, epoch);
    let s2050 = star::MeanPlace {
        asc,
        dec,
        epoch,
        ..s2000
    };

    for method in [star::AberrMethod::Classical, star::AberrMethod::RonVondrak] {
        for JD in [JD_2028, precess::jd_frm_julian_epoch(2075.0)] {
            let p1 = star::apparent_place(&s2000, JD, method);
            let p2 = star::apparent_place(&s2050, JD, method);
            assert!(
                sep_asec(p1.asc, p1.dec, p2.asc, p2.dec) < 0.002,
                "{method:?} {JD}"
            );
        }
    }
}

#[test]
fn apparent_place_ranges() {
    // Right ascension wraps at 0h, northern and southern stars, pole
    for (a, d, pma) in [
        (TAU - 1e-7, 0.3, 1e-6),
        (1e-7, -0.3, -1e-6),
        (3.0, -FRAC_PI_2 + 1e-6, 0.0),
        (0.66, 1.5577, 3e-6),
    ] {
        let s = star::MeanPlace {
            asc: a,
            dec: d,
            pm_asc: pma,
            pm_dec: 0.0,
            epoch: precess::J2000,
        };
        for JD in [
            precess::B1900,
            JD_2028,
            precess::jd_frm_julian_epoch(2300.0),
        ] {
            for method in [star::AberrMethod::Classical, star::AberrMethod::RonVondrak] {
                let p = star::apparent_place(&s, JD, method);
                assert!((0.0..TAU).contains(&p.asc) && (0.0..TAU).contains(&p.asc_pm));
                assert!(p.dec.abs() <= FRAC_PI_2 && p.dec.is_finite());
                assert!(p.nut_asc.abs() < PI);
            }
        }
    }
}

// Polaris, J2000.0 -> 2028 (no book value; check against the mean place
// and the sizes of the corrections)
#[test]
fn apparent_place_of_polaris() {
    let polaris = star::MeanPlace {
        asc: hms(2, 31, 48.704).to_radians(),
        dec: dms(false, 89, 15, 50.72).to_radians(),
        pm_asc: hms(0, 0, 0.19877).to_radians(),
        pm_dec: (-0.0152 / 3600.0_f64).to_radians(),
        epoch: precess::J2000,
    };
    let p = star::apparent_place(&polaris, JD_2028, star::AberrMethod::RonVondrak);
    let (a, d) = precess::precess_eq_coords(p.asc_pm, p.dec_pm, precess::J2000, JD_2028);
    let shift = sep_asec(a, d, p.asc, p.dec);
    // aberration (≤ 20.5") and nutation (≤ about 20") together
    assert!(shift > 1.0 && shift < 45.0, "{shift}");
    assert_close_deg(p.dec.to_degrees(), 89.4, 0.1);
}

/// Direction of `u + v` for the unit vector `u` of (asc, dec): the
/// aberration of a star for an observer moving with velocity `v` (in units
/// of the speed of light), to first order, without the 1/cos(dec) of the
/// formulae in right ascension.
fn aberrate_by_vector(asc: f64, dec: f64, v: [f64; 3]) -> (f64, f64) {
    let x = dec.cos() * asc.cos() + v[0];
    let y = dec.cos() * asc.sin() + v[1];
    let z = dec.sin() + v[2];
    (y.atan2(x).rem_euclid(TAU), z.atan2(x.hypot(y)))
}

// Stars at and next to the celestial poles: the whole reduction, done
// independently with the aberration applied as a vector, must agree with
// `apparent_place` for both methods. Adding first-order corrections to α and
// δ fails here (1/cos δ): before the fix the Ron–Vondrák result was 0.56″
// off at δ = 89.9°, 5.6″ at 89.99° and about 20″ at the pole itself.
#[test]
fn apparent_place_at_the_celestial_poles() {
    let c = 173.14463350; // speed of light, AU/day
    for JD in [precess::J2000 + 100.0, JD_2028] {
        let (dpsi, deps) = nutation::nutation(JD);
        let eps0 = ecliptic::mn_oblq_IAU(JD);
        for dec in [89.9_f64, 89.99, 89.999, 90.0, -89.95, -90.0] {
            let s = star::MeanPlace {
                asc: 1.0,
                dec: dec.to_radians(),
                pm_asc: 0.0,
                pm_dec: 0.0,
                epoch: precess::J2000,
            };

            // Ron–Vondrák: aberration in the J2000.0 frame, then precession
            let (x, y, z) = aberr::earth_bary_vel(JD);
            let (a, d) = aberrate_by_vector(s.asc, s.dec, [x / c, y / c, z / c]);
            let (a, d) = precess::precess_eq_coords(a, d, precess::J2000, JD);
            let (a, d) = nutation::nutate_eq_coords(a, d, dpsi, deps, eps0);
            let p = star::apparent_place(&s, JD, star::AberrMethod::RonVondrak);
            let err = sep_asec(a, d, p.asc, p.dec);
            assert!(err < 0.001, "Ron–Vondrák, {JD} {dec}: {err}\"");

            // Classical: precession, then the velocity of eq. 23.2 turned into
            // equatorial coordinates of the date
            let (a, d) = precess::precess_eq_coords(s.asc, s.dec, precess::J2000, JD);
            let k = aberr::CONST_OF_ABERR;
            let e = aberr::earth_orb_ecc(JD);
            let pi = aberr::earth_perih_long(JD);
            let sun = sun_long_low_accuracy(JD);
            let vx = k * (sun.sin() - e * pi.sin());
            let vy = k * (e * pi.cos() - sun.cos());
            let (a, d) = aberrate_by_vector(a, d, [vx, vy * eps0.cos(), vy * eps0.sin()]);
            let (a, d) = nutation::nutate_eq_coords(a, d, dpsi, deps, eps0);
            let p = star::apparent_place(&s, JD, star::AberrMethod::Classical);
            let err = sep_asec(a, d, p.asc, p.dec);
            assert!(err < 0.001, "classical, {JD} {dec}: {err}\"");
            assert!(p.dec.abs() <= FRAC_PI_2);
        }
    }
}

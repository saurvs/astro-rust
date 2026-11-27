//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 40: correction for
//! parallax.

mod common;

use astro::*;
use common::{arcsec, assert_close, assert_close_deg, assert_close_rad, dms, hms};

// Palomar Observatory, Meeus example 11.a / 40.a.
const RHO_SIN: f64 = 0.546861;
const RHO_COS: f64 = 0.836339;

fn rad(deg: f64) -> f64 {
    deg.to_radians()
}

/// Topocentric right ascension, declination and distance ratio by plain
/// vector subtraction (independent of the library's formulas).
fn topocent_by_vectors(
    asc: f64,
    dec: f64,
    hr_angl: f64,
    parllx: f64,
    rho_sin: f64,
    rho_cos: f64,
) -> (f64, f64, f64) {
    let d = 1.0 / parllx.sin(); // geocentric distance in Earth radii
    let theta = hr_angl + asc; // local sidereal time
    let x = d * dec.cos() * asc.cos() - rho_cos * theta.cos();
    let y = d * dec.cos() * asc.sin() - rho_cos * theta.sin();
    let z = d * dec.sin() - rho_sin;
    let asc_1 = y.atan2(x).rem_euclid(std::f64::consts::TAU);
    let dec_1 = z.atan2(x.hypot(y));
    let dist_ratio = (x * x + y * y + z * z).sqrt() / d;
    (asc_1, dec_1, dist_ratio)
}

// Meeus, example 40.a: π = 23.592" for Δ = 0.37276 AU.
#[test]
fn eq_hz_parallax_example_40a() {
    let p = parallax::eq_hz_parallax(0.37276);
    assert_close(p.to_degrees() * 3600.0, 23.592, 0.0005);
}

#[test]
fn eq_hz_parallax_limits() {
    // At 1 AU the parallax is the solar parallax.
    assert_close(
        parallax::eq_hz_parallax(1.0).to_degrees() * 3600.0,
        8.794,
        1e-9,
    );
    // The Moon at 368410 km (example 48.a): close to 0°59'31".
    let p = parallax::eq_hz_parallax(368410.0 / 149597870.7);
    assert_close(p.to_degrees(), 0.99, 0.01);
    // Inside the Earth: NaN, no panic.
    assert!(parallax::eq_hz_parallax(1e-6).is_nan());
}

// Meeus, example 40.a: geocentric hour angle H = θ0 − L − α with the
// longitude positive west. (The book rounds α to 22h38m07s here.)
#[test]
fn hour_angle_example_40a() {
    let h = coords::hr_angl_frm_observer_long(
        rad(hms(1, 40, 45.0)),
        rad(hms(7, 47, 27.0)),
        rad(hms(22, 38, 7.0)),
    );
    assert_close(h.to_degrees(), -431.2042, 5e-5);
    assert_close_deg(h.to_degrees(), 288.7958, 5e-5);
}

// Meeus, example 40.a, rigorous formulas 40.2 and 40.3.
#[test]
fn topocent_eq_example_40a() {
    let asc = rad(339.530208);
    let dec = rad(-15.771083);
    let parllx = rad(arcsec(23.592));
    let h = rad(288.7958);

    let topo = parallax::topocent_eq(asc, dec, h, parllx, RHO_SIN, RHO_COS);

    // tan Δα = +0.000090557 / +0.962324, Δα = +0.0053917° = +1.29 s
    assert_close(topo.del_asc.tan(), 0.000090557 / 0.962324, 5e-10);
    assert_close(topo.del_asc.to_degrees(), 0.0053917, 5e-8);
    assert_close(topo.del_asc.to_degrees() * 240.0, 1.29, 0.005);
    // α' = 22h38m08.54s
    assert_close_deg(
        topo.asc.to_degrees(),
        hms(22, 38, 8.54),
        0.005 * 15.0 / 3600.0,
    );
    // tan δ' = −0.27185713 / +0.96232447, δ' = −15°46'30.0"
    assert_close(topo.dec, (-0.27185713_f64).atan2(0.96232447), 1e-8);
    assert_close(topo.dec.to_degrees(), dms(true, 15, 46, 30.0), arcsec(0.05));

    assert_close(topo.del_dec, topo.dec - dec, 1e-15);
    assert_close_rad(topo.hr_angl, h - topo.del_asc, 1e-15);
}

// Meeus, example 40.a with the observer given by latitude and height, and the
// hour angle derived from the Greenwich sidereal time.
#[test]
fn topocent_eq_coords_example_40a() {
    let (rho_sin, rho_cos) = planet::earth::rho_sin_cos_phi(rad(dms(false, 33, 21, 22.0)), 1706.0);
    assert_close(rho_sin, RHO_SIN, 1e-6);
    assert_close(rho_cos, RHO_COS, 1e-6);

    let eq_point = coords::EqPoint {
        asc: rad(hms(22, 38, 7.25)),
        dec: rad(dms(true, 15, 46, 15.9)),
    };
    let observer = coords::GeographPoint {
        long: rad(hms(7, 47, 27.0)), // west: positive
        lat: rad(dms(false, 33, 21, 22.0)),
    };
    let topo = parallax::topocent_eq_coords(
        &eq_point,
        rad(arcsec(23.592)),
        &observer,
        1706.0,
        rad(hms(1, 40, 45.0)),
    );
    assert_close_deg(
        topo.asc.to_degrees(),
        hms(22, 38, 8.54),
        0.006 * 15.0 / 3600.0,
    );
    assert_close(topo.dec.to_degrees(), dms(true, 15, 46, 30.0), arcsec(0.06));
}

// Meeus, example 40.a, non-rigorous formulas 40.4 and 40.5:
// Δα = +19.409" = +1.29 s, Δδ = −14.1", δ' = −15°46'30.0".
#[test]
fn topocent_eq_approx_example_40a() {
    let asc = rad(339.530208);
    let dec = rad(-15.771083);
    let topo = parallax::topocent_eq_approx(
        asc,
        dec,
        rad(288.7958),
        rad(arcsec(23.592)),
        RHO_SIN,
        RHO_COS,
    );

    assert_close(topo.del_asc.to_degrees() * 3600.0, 19.409, 0.0005);
    assert_close(topo.del_asc.to_degrees() * 240.0, 1.29, 0.005);
    assert_close(topo.del_dec.to_degrees() * 3600.0, -14.1, 0.05);
    assert_close(topo.dec.to_degrees(), dms(true, 15, 46, 30.0), arcsec(0.05));
    assert_close_deg(
        topo.asc.to_degrees(),
        hms(22, 38, 8.54),
        0.005 * 15.0 / 3600.0,
    );

    // For a planet both methods agree to a few milliarcseconds.
    let rig = parallax::topocent_eq(
        asc,
        dec,
        rad(288.7958),
        rad(arcsec(23.592)),
        RHO_SIN,
        RHO_COS,
    );
    assert_close(topo.del_asc, rig.del_asc, rad(arcsec(0.005)));
    assert_close(topo.del_dec, rig.del_dec, rad(arcsec(0.005)));
}

// Meeus' exercise after example 40.a: the Moon with fictive values
// α = 15°, δ = +5°, H = +60°, π = 0°59'. Formulas 40.2/40.3 and 40.6/40.7
// must agree exactly; the non-rigorous formulas do not for the Moon.
#[test]
fn topocent_eq_moon_exercise() {
    let (asc, dec, h, parllx) = (rad(15.0), rad(5.0), rad(60.0), rad(dms(false, 0, 59, 0.0)));
    let topo = parallax::topocent_eq(asc, dec, h, parllx, RHO_SIN, RHO_COS);

    // Literal 40.2 and 40.3.
    let s = parllx.sin();
    let den = dec.cos() - RHO_COS * s * h.cos();
    let del_asc = (-RHO_COS * s * h.sin()).atan2(den);
    let dec_40_3 = ((dec.sin() - RHO_SIN * s) * del_asc.cos()).atan2(den);
    assert_close(topo.del_asc, del_asc, 1e-14);
    assert_close(topo.dec, dec_40_3, 1e-14);

    // Literal 40.6 and 40.7.
    let a = dec.cos() * h.sin();
    let b = dec.cos() * h.cos() - RHO_COS * s;
    let c = dec.sin() - RHO_SIN * s;
    let q = (a * a + b * b + c * c).sqrt();
    assert_close_rad(topo.hr_angl, a.atan2(b), 1e-14);
    assert_close(topo.dec, (c / q).asin(), 1e-14);

    let approx = parallax::topocent_eq_approx(asc, dec, h, parllx, RHO_SIN, RHO_COS);
    assert!((approx.del_asc - topo.del_asc).abs() > rad(arcsec(1.0)));
    assert!((approx.del_asc - topo.del_asc).abs() < rad(arcsec(120.0)));
    assert!((approx.del_dec - topo.del_dec).abs() < rad(arcsec(120.0)));
}

// Compare with plain vector subtraction for a set of bodies and observers:
// northern/southern, near the poles, both sides of the meridian, the Moon
// and a planet.
#[test]
fn topocent_eq_matches_vectors() {
    let observers = [
        (RHO_SIN, RHO_COS),
        (-0.55, 0.83),   // southern hemisphere
        (0.0, 1.0),      // equator
        (0.9966, 0.0),   // north pole
        (-0.9966, 0.01), // near the south pole
    ];
    let bodies = [
        (15.0, 5.0, 60.0, 0.98),
        (359.9, -28.0, 300.0, 0.95),
        (180.0, 89.5, 0.0, 1.0),
        (90.0, -89.9, 170.0, 0.9),
        (250.0, 0.0, 359.99, 0.0065), // planet-like
        (10.0, 45.0, -45.0, 1.02),    // negative hour angle
    ];
    for &(rho_sin, rho_cos) in &observers {
        for &(asc, dec, h, p) in &bodies {
            let (asc, dec, h, p) = (rad(asc), rad(dec), rad(h), rad(p));
            let topo = parallax::topocent_eq(asc, dec, h, p, rho_sin, rho_cos);
            let (asc_v, dec_v, _) = topocent_by_vectors(asc, dec, h, p, rho_sin, rho_cos);
            assert_close_rad(topo.asc, asc_v, 1e-11);
            assert_close(topo.dec, dec_v, 1e-11);
            assert!((0.0..std::f64::consts::TAU).contains(&topo.asc));
            assert!((0.0..std::f64::consts::TAU).contains(&topo.hr_angl));
            assert!(topo.dec.abs() <= std::f64::consts::FRAC_PI_2);
        }
    }
}

// Right ascension just below 360° and a parallax that increases it: the
// result must wrap to [0, 2π). Before the fix it came out above 2π.
#[test]
fn topocent_eq_coords_wraps_asc() {
    let eq_point = coords::EqPoint {
        asc: rad(359.9999),
        dec: rad(10.0),
    };
    let observer = coords::GeographPoint {
        long: 0.0,
        lat: rad(40.0),
    };
    // θ0 = α + H with H = 300° (east of the meridian, so Δα > 0).
    let topo = parallax::topocent_eq_coords(&eq_point, rad(0.95), &observer, 0.0, rad(299.9999));
    assert!(
        (0.0..std::f64::consts::TAU).contains(&topo.asc),
        "{}",
        topo.asc
    );
    assert!(
        topo.asc.to_degrees() > 0.1 && topo.asc.to_degrees() < 1.0,
        "{}",
        topo.asc
    );
}

// A Moon-like body near the celestial pole seen from the equator, where
// the denominator of eq. 40.3 is negative. Before the fix δ' came out as
// about −90.5°, outside [−90°, 90°].
#[test]
fn topocent_eq_coords_near_pole() {
    let eq_point = coords::EqPoint {
        asc: rad(100.0),
        dec: rad(89.5),
    };
    let observer = coords::GeographPoint {
        long: 0.0,
        lat: 0.0,
    };
    let topo = parallax::topocent_eq_coords(&eq_point, rad(1.0), &observer, 0.0, rad(100.0));
    let (rho_sin, rho_cos) = planet::earth::rho_sin_cos_phi(0.0, 0.0);
    let (asc_v, dec_v, _) =
        topocent_by_vectors(rad(100.0), rad(89.5), 0.0, rad(1.0), rho_sin, rho_cos);
    assert!(
        topo.dec.abs() <= std::f64::consts::FRAC_PI_2,
        "{}",
        topo.dec.to_degrees()
    );
    assert_close(topo.dec, dec_v, 1e-11);
    assert_close_rad(topo.asc, asc_v, 1e-11);
    // The body is seen beyond the pole: on the opposite meridian.
    assert_close_deg(topo.asc.to_degrees(), 280.0, 1e-6);
}

#[test]
fn topocent_eq_invariants() {
    let (asc, dec, p) = (rad(123.0), rad(-12.0), rad(0.9));

    // Observer at the centre of the Earth: no parallax.
    let t = parallax::topocent_eq(asc, dec, rad(33.0), p, 0.0, 0.0);
    assert_close_rad(t.asc, asc, 1e-15);
    assert_close(t.dec, dec, 1e-15);
    assert_eq!(t.del_asc, 0.0);

    // Body in the direction of the observer's geocentric zenith: no shift.
    let phi_1 = rad(-35.0); // southern observer
    let t = parallax::topocent_eq(asc, phi_1, 0.0, p, 0.99 * phi_1.sin(), 0.99 * phi_1.cos());
    assert_close(t.del_asc, 0.0, 1e-15);
    assert_close(t.del_dec, 0.0, 1e-15);

    // Mirror symmetry about the meridian: H -> -H flips Δα, keeps Δδ.
    let east = parallax::topocent_eq(asc, dec, rad(-50.0), p, RHO_SIN, RHO_COS);
    let west = parallax::topocent_eq(asc, dec, rad(50.0), p, RHO_SIN, RHO_COS);
    assert_close(east.del_asc, -west.del_asc, 1e-15);
    assert_close(east.del_dec, west.del_dec, 1e-15);
    // West of the meridian the parallax decreases the right ascension.
    assert!(west.del_asc < 0.0);
    // The body is always displaced away from the observer's zenith; for a
    // northern observer and a body on the meridian south of the zenith the
    // declination decreases.
    let merid = parallax::topocent_eq(asc, dec, 0.0, p, RHO_SIN, RHO_COS);
    assert_close(merid.del_asc, 0.0, 1e-15);
    assert!(merid.del_dec < 0.0);

    // An hour angle given outside [0, 2π) gives the same position.
    let a = parallax::topocent_eq(asc, dec, rad(50.0), p, RHO_SIN, RHO_COS);
    let b = parallax::topocent_eq(asc, dec, rad(50.0 - 720.0), p, RHO_SIN, RHO_COS);
    assert_close_rad(a.asc, b.asc, 1e-12);
    assert_close(a.dec, b.dec, 1e-12);
    assert_close_rad(a.hr_angl, b.hr_angl, 1e-12);
}

// Meeus' ecliptic exercise in chapter 40 (the Moon):
// λ' = 181°48'05.0", β' = +1°29'07.1", s' = 0°16'25.5".
#[test]
fn topocent_ecl_coords_exercise() {
    let ecl_point = coords::EclPoint {
        long: rad(dms(false, 181, 46, 22.5)),
        lat: rad(dms(false, 2, 17, 26.2)),
    };
    let observer = coords::GeographPoint {
        long: 0.0, // unused: the longitude enters through the local sidereal time
        lat: rad(dms(false, 50, 5, 7.8)),
    };
    let parllx = rad(dms(false, 0, 59, 27.7));
    let s = rad(dms(false, 0, 16, 15.5));
    let eps = rad(dms(false, 23, 28, 0.8));
    let theta = rad(dms(false, 209, 46, 7.9));

    let (topo, s_1) =
        parallax::topocent_ecl_coords(&ecl_point, parllx, &observer, 0.0, theta, eps, s);
    assert_close_deg(
        topo.long.to_degrees(),
        dms(false, 181, 48, 5.0),
        arcsec(0.05),
    );
    assert_close(topo.lat.to_degrees(), dms(false, 1, 29, 7.1), arcsec(0.05));
    assert_close(s_1.to_degrees(), dms(false, 0, 16, 25.5), arcsec(0.05));

    #[allow(deprecated)]
    let (old, old_s) =
        parallax::topopcent_ecl_coords(&ecl_point, parllx, &observer, 0.0, theta, eps, s);
    assert_eq!((old.long, old.lat, old_s), (topo.long, topo.lat, s_1));

    // Same through the low-level function and Meeus' literal formulas
    // (valid here because tan β' is evaluated with a one-argument atan).
    let (rho_sin, rho_cos) = planet::earth::rho_sin_cos_phi(observer.lat, 0.0);
    let t = parallax::topocent_ecl(
        ecl_point.long,
        ecl_point.lat,
        s,
        parllx,
        rho_sin,
        rho_cos,
        theta,
        eps,
    );
    assert_eq!((t.long, t.lat, t.semidiameter), (topo.long, topo.lat, s_1));

    let (l, b, sp) = (ecl_point.long, ecl_point.lat, parllx.sin());
    let n = l.cos() * b.cos() - rho_cos * sp * theta.cos();
    let l_1 = (l.sin() * b.cos() - sp * (rho_sin * eps.sin() + rho_cos * eps.cos() * theta.sin()))
        .atan2(n);
    let b_1 = (l_1.cos()
        * (b.sin() - sp * (rho_sin * eps.cos() - rho_cos * eps.sin() * theta.sin()))
        / n)
        .atan();
    let s_lit = (l_1.cos() * b_1.cos() * s.sin() / n).asin();
    assert_close_rad(t.long, l_1, 1e-14);
    assert_close(t.lat, b_1, 1e-14);
    assert_close(t.semidiameter, s_lit, 1e-14);
}

// The topocentric ecliptic position must be the topocentric equatorial
// position expressed in ecliptic coordinates, for any longitude (including
// λ' near 90° and 270°, where N = 0) and for southern observers.
#[test]
fn topocent_ecl_consistent_with_eq() {
    let eps = rad(23.44);
    let p = rad(0.95);
    let s = rad(0.26);
    for &(rho_sin, rho_cos) in &[(RHO_SIN, RHO_COS), (-0.7, 0.71), (0.0, 1.0)] {
        for &lambda in &[
            0.0, 45.0, 89.9, 90.0, 135.0, 180.0, 181.77, 225.0, 270.0, 315.0, 359.99,
        ] {
            for &beta in &[-5.0, 0.0, 2.3, 60.0] {
                for &theta in &[0.0, 100.0, 209.77, 330.0] {
                    let (l, b, th) = (rad(lambda), rad(beta), rad(theta));
                    let t = parallax::topocent_ecl(l, b, s, p, rho_sin, rho_cos, th, eps);

                    let asc = coords::asc_frm_ecl(l, b, eps);
                    let dec = coords::dec_frm_ecl(l, b, eps);
                    let (asc_v, dec_v, ratio) =
                        topocent_by_vectors(asc, dec, th - asc, p, rho_sin, rho_cos);
                    let l_v = coords::ecl_long_frm_eq(asc_v, dec_v, eps);
                    let b_v = coords::ecl_lat_frm_eq(asc_v, dec_v, eps);

                    assert_close_rad(t.long, l_v, 1e-10);
                    assert_close(t.lat, b_v, 1e-10);
                    assert_close(t.semidiameter, (s.sin() / ratio).asin(), 1e-12);
                    assert!((0.0..std::f64::consts::TAU).contains(&t.long));
                }
            }
        }
    }
}

#[test]
fn topocent_ecl_no_parallax() {
    let t = parallax::topocent_ecl(
        rad(200.0),
        rad(-3.0),
        rad(0.25),
        0.0,
        RHO_SIN,
        RHO_COS,
        1.0,
        0.4,
    );
    assert_close_rad(t.long, rad(200.0), 1e-15);
    assert_close(t.lat, rad(-3.0), 1e-15);
    assert_close(t.semidiameter, rad(0.25), 1e-15);
}

// Parallax in altitude, sin p = ρ sin π cos h.
#[test]
fn parallax_in_alt() {
    let p = rad(0.95);
    // On the horizon the parallax equals ρ π (for ρ = 1, exactly π).
    assert_close(parallax::parallax_in_alt(0.0, p, 1.0), p, 1e-15);
    // At the zenith it vanishes.
    assert_close(parallax::parallax_in_alt(rad(90.0), p, 1.0), 0.0, 1e-15);
    // Symmetric in the sign of the altitude.
    assert_close(
        parallax::parallax_in_alt(rad(-10.0), p, 0.998),
        parallax::parallax_in_alt(rad(10.0), p, 0.998),
        1e-15,
    );

    // Compare with the exact geometry on a spherical Earth (ρ = 1): with
    // the topocentric altitude the formula is exact, with the geocentric
    // one it is good to a few arcminutes for the Moon.
    let d = 1.0 / p.sin();
    for &h in &[5.0, 30.0, 60.0] {
        let h = rad(h);
        let h_topo = (d * h.sin() - 1.0).atan2(d * h.cos());
        let exact = h - h_topo;
        assert_close(parallax::parallax_in_alt(h_topo, p, 1.0), exact, 1e-12);
        assert_close(
            parallax::parallax_in_alt(h, p, 1.0),
            exact,
            rad(arcsec(60.0)),
        );
    }
}

// An observer east of Greenwich has a negative longitude, so the hour
// angle is θ0 + |L| − α. Check the convenience function against the vector
// method for Sydney (151°12' E, 33°52' S, 58 m).
#[test]
fn topocent_eq_coords_eastern_southern_observer() {
    let observer = coords::GeographPoint {
        long: rad(-dms(false, 151, 12, 0.0)),
        lat: rad(-dms(false, 33, 52, 0.0)),
    };
    let eq_point = coords::EqPoint {
        asc: rad(250.0),
        dec: rad(-20.0),
    };
    let theta0 = rad(80.0);
    let p = rad(0.93);
    let topo = parallax::topocent_eq_coords(&eq_point, p, &observer, 58.0, theta0);

    let h = theta0 + rad(dms(false, 151, 12, 0.0)) - eq_point.asc;
    let (rho_sin, rho_cos) = planet::earth::rho_sin_cos_phi(observer.lat, 58.0);
    assert!(rho_sin < 0.0);
    let (asc_v, dec_v, _) = topocent_by_vectors(eq_point.asc, eq_point.dec, h, p, rho_sin, rho_cos);
    assert_close_rad(topo.asc, asc_v, 1e-11);
    assert_close(topo.dec, dec_v, 1e-11);
    // H = 80° + 151.2° − 250° ≈ −18.8°: the body is east of the meridian,
    // so the parallax increases its right ascension (Δα > 0); it lies north
    // of this southern observer's zenith, so it is displaced further north
    // (Δδ > 0).
    let t = parallax::topocent_eq(eq_point.asc, eq_point.dec, h, p, rho_sin, rho_cos);
    assert_close(h.to_degrees(), -18.8, 1e-9);
    assert!(t.del_asc > 0.0);
    assert!(t.del_dec > 0.0);
}

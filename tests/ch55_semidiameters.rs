//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 55: semidiameters of
//! the Sun, Moon and planets, and diameters of asteroids.

mod common;

use astro::planet::{self, Planet};
use astro::{asteroid, lunar, pluto, sun};
use common::{assert_close, dms};

/// Arcseconds per radian.
const AS: f64 = 206264.80624709636;

fn arcsec(rad: f64) -> f64 {
    rad * AS
}

#[test]
fn sun_semidiameter_is_in_radians() {
    // s0 = 15'59.63" = 959.63" at 1 AU.
    assert_close(arcsec(sun::semidiameter(1.0)), 959.63, 1e-9);
    // Near perihelion and aphelion.
    assert_close(arcsec(sun::semidiameter(0.98329)), 959.63 / 0.98329, 1e-9);
    assert_close(arcsec(sun::semidiameter(1.01671)), 959.63 / 1.01671, 1e-9);
}

// The s0 values of Meeus' table "B" (Astronomical Almanac 1984), except
// Venus, for which the book recommends the older 8.41".
#[test]
fn planet_semidiameters_at_unit_distance() {
    let table = [
        (Planet::Mercury, 3.36),
        (Planet::Venus, 8.41),
        (Planet::Mars, 4.68),
        (Planet::Jupiter, 98.44),
        (Planet::Saturn, 82.73),
        (Planet::Uranus, 35.02),
        (Planet::Neptune, 33.50),
    ];
    for (p, s0) in table {
        assert_close(arcsec(planet::semidiameter(&p, 1.0).unwrap()), s0, 1e-9);
        // s = s0 / Δ
        assert_close(
            arcsec(planet::semidiameter(&p, 4.0).unwrap()),
            s0 / 4.0,
            1e-9,
        );
    }
    assert!(planet::semidiameter(&Planet::Earth, 1.0).is_err());
}

#[test]
fn jupiter_semidiameters() {
    use astro::planet::jupiter;
    assert_close(arcsec(jupiter::eq_semidiameter(1.0)), 98.44, 1e-9);
    assert_close(arcsec(jupiter::pol_semidiameter(1.0)), 92.06, 1e-9);
    assert_close(arcsec(jupiter::eq_semidiameter(4.2)), 98.44 / 4.2, 1e-9);
    assert_close(arcsec(jupiter::pol_semidiameter(4.2)), 92.06 / 4.2, 1e-9);

    // Rigorous apparent polar semidiameter: b/Δ with the axis in the sky
    // plane, a/Δ with the pole towards the Earth, and within 0.05% of b/Δ
    // for the largest possible D_E of about 4°.
    let d = 4.2;
    let b = jupiter::pol_semidiameter(d);
    assert_close(jupiter::apprnt_pol_semidiameter(d, 0.0), b, 1e-15);
    assert_close(
        jupiter::apprnt_pol_semidiameter(d, std::f64::consts::FRAC_PI_2),
        jupiter::eq_semidiameter(d),
        1e-15,
    );
    let s4 = jupiter::apprnt_pol_semidiameter(d, 4_f64.to_radians());
    assert!(s4 > b && (s4 - b) / b < 5e-4);
    // Symmetric in the sign of D_E.
    assert_close(
        jupiter::apprnt_pol_semidiameter(d, -3_f64.to_radians()),
        jupiter::apprnt_pol_semidiameter(d, 3_f64.to_radians()),
        1e-18,
    );
}

#[test]
fn saturn_semidiameters() {
    use astro::planet::saturn;
    assert_close(arcsec(saturn::eq_semidiameter(1.0)), 82.73, 1e-9);
    assert_close(arcsec(saturn::eq_semidiameter(9.5)), 82.73 / 9.5, 1e-9);

    // k = 1 - (b/a)^2 = 0.203800 for a = 82.73", b = 73.82".
    let k = 1.0 - (73.82_f64 / 82.73).powi(2);
    assert_close(k, 0.203800, 5e-7);

    // B = 0: the polar semidiameter is b/Δ; B = ±90°: a/Δ.
    assert_close(arcsec(saturn::pol_semidiameter(1.0, 0.0)), 73.82, 1e-9);
    assert_close(
        arcsec(saturn::pol_semidiameter(1.0, std::f64::consts::FRAC_PI_2)),
        82.73,
        1e-9,
    );
    // A typical tilt, B = 16.442° (Meeus, example 45.a), Δ = 8.6892 AU.
    let b = 16.442_f64.to_radians();
    let expected = 82.73 / 8.6892 * (1.0 - k * b.cos().powi(2)).sqrt();
    assert_close(arcsec(saturn::pol_semidiameter(8.6892, b)), expected, 1e-9);
    assert_close(
        saturn::pol_semidiameter(8.6892, -b),
        saturn::pol_semidiameter(8.6892, b),
        1e-18,
    );
}

#[test]
fn pluto_semidiameter() {
    assert_close(arcsec(pluto::semidiameter(1.0)), 2.07, 1e-9);
    assert_close(arcsec(pluto::semidiameter(30.0)), 0.069, 1e-9);
    #[allow(deprecated)]
    let old = pluto::semdiameter(30.0);
    assert_eq!(old, pluto::semidiameter(30.0));
}

// Meeus' shortcut s = 358 473 400 / Δ (arcseconds, Δ in km) is said to be
// within 0.0005" of the rigorous formula.
#[test]
fn moon_geocentric_semidiameter() {
    let mut d = 356_000.0;
    while d <= 407_000.0 {
        let s = arcsec(lunar::semidiameter(d));
        assert_close(s, 358_473_400.0 / d, 0.0005);
        // sin s = k sin π
        let sin_pi = 6378.14 / d;
        assert_close((s / AS).sin(), 0.272481 * sin_pi, 1e-15);
        d += 500.0;
    }
}

// Meeus, chapter 40 exercise (parallax in ecliptic coordinates):
// λ = 181°46'22.5", β = +2°17'26.2", π = 0°59'27.7", s = 0°16'15.5",
// φ = +50°05'07.8" at sea level, ε = 23°28'00.8", local sidereal time
// θ = 209°46'07.9". Answer: s' = 0°16'25.5".
#[test]
fn moon_topocentric_semidiameter_ch40_exercise() {
    let lam = dms(false, 181, 46, 22.5).to_radians();
    let bet = dms(false, 2, 17, 26.2).to_radians();
    let par = dms(false, 0, 59, 27.7).to_radians();
    let s = dms(false, 0, 16, 15.5).to_radians();
    let phi = dms(false, 50, 5, 7.8).to_radians();
    let eps = dms(false, 23, 28, 0.8).to_radians();
    let theta = dms(false, 209, 46, 7.9).to_radians();

    // Equatorial coordinates (Meeus, eq. 13.3 and 13.4).
    let alpha = (lam.sin() * eps.cos() - bet.tan() * eps.sin()).atan2(lam.cos());
    let delta = (bet.sin() * eps.cos() + bet.cos() * eps.sin() * lam.sin()).asin();
    let h = theta - alpha;
    let (rho_sin, rho_cos) = rho_sin_cos(phi);

    let s1 = lunar::topocent_semidiameter(s, par, delta, h, rho_sin, rho_cos);
    assert_close(arcsec(s1), dms(false, 0, 16, 25.5) * 3600.0, 0.05);

    // The approximate rule s (1 + sin h sin π) agrees to about 0.1".
    let alt = (phi.sin() * delta.sin() + phi.cos() * delta.cos() * h.cos()).asin();
    let s2 = lunar::topocent_semidiameter_frm_alt(s, par, alt);
    assert_close(arcsec(s2), arcsec(s1), 0.1);
}

#[test]
fn moon_topocentric_semidiameter_limits() {
    for &d in &[356_400.0, 384_400.0, 406_700.0] {
        let s = lunar::semidiameter(d);
        let par = lunar::eq_hz_parllx(d);

        // Observer at the north pole, Moon at the zenith (δ = 90°): the
        // distance shrinks by one Earth radius (ρ sin φ' = b/a).
        let b_a = 1.0 - 1.0 / 298.257;
        let s_z = lunar::topocent_semidiameter(s, par, std::f64::consts::FRAC_PI_2, 0.0, b_a, 0.0);
        let expected = (0.272481 * 6378.14 / (d - 6378.14 * b_a)).asin();
        assert_close(s_z, expected, 1e-12);

        // The increase at the zenith lies between 14" and 18" (Meeus).
        let inc =
            arcsec(lunar::topocent_semidiameter_frm_alt(s, par, std::f64::consts::FRAC_PI_2) - s);
        assert!((13.5..18.5).contains(&inc), "increase {inc}");
        // The first-order rule misses about s sin²π (0.3") here.
        assert_close(arcsec(s_z - s), inc, 0.4);

        // Moon on the horizon of an equatorial observer (H = 90°, δ = 0):
        // the topocentric distance is sqrt(Δ² + R²), so s' is smaller than
        // s by only about s sin²π / 2 (0.1" to 0.2").
        let s_h = lunar::topocent_semidiameter(s, par, 0.0, std::f64::consts::FRAC_PI_2, 0.0, 1.0);
        let q = (1.0 + par.sin().powi(2)).sqrt();
        assert_close(s_h, (s.sin() / q).asin(), 1e-15);
        assert!(s_h < s && arcsec(s - s_h) < 0.2);
        assert_eq!(lunar::topocent_semidiameter_frm_alt(s, par, 0.0), s);

        // Below the horizon the Moon looks smaller (observer farther away).
        let s_n = lunar::topocent_semidiameter(s, par, 0.0, std::f64::consts::PI, 0.0, 1.0);
        assert!(s_n < s);
    }
}

#[test]
fn asteroid_diameter() {
    // log10 d = 3.12 - H/5 - 0.5 log10 A, d in km; the function returns m.
    assert_close(
        asteroid::diameter(15.0, 0.04) / 1000.0,
        10_f64.powf(3.12 - 3.0 - 0.5 * 0.04_f64.log10()),
        1e-9,
    );
    // Same result through the natural-log form of the book.
    let x = 3.12 - 15.0 / 5.0 - 0.217147 * 0.04_f64.ln();
    assert_close(
        asteroid::diameter(15.0, 0.04) / 1000.0,
        10_f64.powf(x),
        1e-4,
    );
    // Ceres: H = 3.34, albedo about 0.09 → about 940 km.
    assert_close(asteroid::diameter(3.34, 0.09) / 1000.0, 943.8, 0.1);
    // One magnitude fainter means 10^0.2 times smaller; four times the
    // albedo means half the diameter.
    assert_close(
        asteroid::diameter(10.0, 0.1) / asteroid::diameter(11.0, 0.1),
        10_f64.powf(0.2),
        1e-12,
    );
    assert_close(
        asteroid::diameter(10.0, 0.1) / asteroid::diameter(10.0, 0.4),
        2.0,
        1e-12,
    );
}

#[test]
fn asteroid_apparent_diameter() {
    // 0.0013788" d / Δ, d in km, Δ in AU.
    assert_close(
        arcsec(asteroid::apparent_diameter(1000.0, 1.0)),
        1.3788,
        0.00005,
    );
    assert_close(
        arcsec(asteroid::apparent_diameter(530.0, 1.7)),
        0.0013788 * 530.0 / 1.7,
        2e-5,
    );
    // Vesta-sized body at opposition, 1.2 AU: about 0.6".
    assert_close(
        arcsec(asteroid::apparent_diameter(525.0, 1.2)),
        0.603,
        0.001,
    );
}

/// ρ sin φ' and ρ cos φ' at sea level (Meeus, chapter 11, IAU 1976 Earth).
fn rho_sin_cos(phi: f64) -> (f64, f64) {
    let b_a = 1.0 - 1.0 / 298.257;
    let u = (b_a * phi.tan()).atan();
    (b_a * u.sin(), u.cos())
}

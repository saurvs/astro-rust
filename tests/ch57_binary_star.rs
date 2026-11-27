//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 57: binary stars.

mod common;

use astro::binary_star::{self, OrbitalElements};
use astro::orbit;
use common::{assert_close, assert_close_deg};

/// η Coronae Borealis (Silbernagel 1929), semimajor axis in arcseconds.
fn eta_crb() -> OrbitalElements {
    OrbitalElements {
        period: 41.623,
        time_of_periastron: 1934.008,
        eccentricity: 0.2763,
        semimajor_axis: 0.907,
        inclination: 59.025_f64.to_radians(),
        asc_node_position_angle: 23.717_f64.to_radians(),
        periastron_longitude: 219.907_f64.to_radians(),
    }
}

/// γ Virginis (Heintz 1990), semimajor axis in arcseconds.
fn gamma_vir() -> OrbitalElements {
    OrbitalElements {
        period: 168.68,
        time_of_periastron: 2005.13,
        eccentricity: 0.885,
        semimajor_axis: 3.697,
        inclination: 148.0_f64.to_radians(),
        asc_node_position_angle: 36.9_f64.to_radians(),
        periastron_longitude: 256.5_f64.to_radians(),
    }
}

// Meeus, example 57.a: every intermediate value for η CrB at 1980.0.
#[test]
fn eta_crb_step_by_step_ex57a() {
    let el = eta_crb();
    let n = binary_star::mn_ann_motion_of_compan(el.period);
    assert_close(n.to_degrees(), 8.64906, 0.000005);

    let m = binary_star::mn_anom_of_compan(n, 1980.0, el.time_of_periastron);
    assert_close(m.to_degrees(), 397.788, 0.0005);

    let ecc_anom = orbit::elliptic::ecc_anom(m, el.eccentricity, 1e-12);
    // The book prints 49.897; the exact value is 49.8965.
    assert_close_deg(ecc_anom.to_degrees(), 49.897, 0.001);

    let r = binary_star::rad_vec(el.semimajor_axis, el.eccentricity, ecc_anom);
    assert_close(r, 0.74557, 0.000005);

    let v = binary_star::true_anom(el.eccentricity, ecc_anom);
    // The book prints 63.416 (from its rounded E); exactly it is 63.4151.
    assert_close_deg(v.to_degrees(), 63.416, 0.001);

    // The book carries its rounded v = 63.416 into the next step.
    let v = 63.416_f64.to_radians();
    let (w, i) = (el.periastron_longitude, el.inclination);
    let num = (v + w).sin() * i.cos();
    let den = (v + w).cos();
    assert_close(num, -0.500813, 0.000005);
    assert_close(den, 0.230440, 0.000005);
    assert_close(num.atan2(den).to_degrees(), -65.291, 0.0005);

    let theta = binary_star::apprnt_coords_angl(el.asc_node_position_angle, v, w, i);
    assert_close(theta.to_degrees(), 318.426, 0.0005);

    let rho = binary_star::anglr_sepr(r, v, w, i);
    assert_close(rho, 0.411, 0.0005);
    // Greaney's form agrees with r cos(v + ω) / cos(θ - Ω).
    assert_close(rho, r * den / num.atan2(den).cos(), 1e-12);
}

// Meeus, example 57.a through the convenience function.
#[test]
fn eta_crb_apparent_position_ex57a() {
    let pos = binary_star::apparent_position(&eta_crb(), 1980.0);
    assert_close(pos.position_angle.to_degrees(), 318.4, 0.05);
    assert_close(pos.separation, 0.411, 0.0005);
    assert_eq!(eta_crb().apparent_position(1980.0), pos);
}

// Meeus, example 57.b: eccentricity of the apparent orbit of η CrB.
#[test]
fn eta_crb_apparent_eccentricity_ex57b() {
    let el = eta_crb();
    let (e, w, i) = (el.eccentricity, el.periastron_longitude, el.inclination);
    // Intermediate quantities A, B, C, D of the book.
    let a = (1.0 - e * e * w.cos().powi(2)) * i.cos().powi(2);
    let b = e * e * w.sin() * w.cos() * i.cos();
    let c = 1.0 - e * e * w.sin().powi(2);
    assert_close(a, 0.25298, 0.000005);
    assert_close(b, 0.01934, 0.000005);
    assert_close(c, 0.96853, 0.00005);
    assert_close((a - c).powi(2) + 4.0 * b * b, 0.51358, 0.00005);

    assert_close(binary_star::ecc_of_apprnt_orb(e, w, i), 0.860, 0.0005);
    assert_close(el.apparent_eccentricity(), 0.860, 0.0005);
}

#[test]
fn apparent_eccentricity_limits() {
    for &e in &[0.0, 0.1, 0.5, 0.885] {
        for k in 0..12 {
            let w = (30.0 * k as f64).to_radians();
            // Face-on orbit: apparent orbit is the true orbit.
            assert_close(binary_star::ecc_of_apprnt_orb(e, w, 0.0), e, 1e-12);
            // Retrograde face-on orbit, same shape.
            assert_close(
                binary_star::ecc_of_apprnt_orb(e, w, std::f64::consts::PI),
                e,
                1e-12,
            );
            // Edge-on orbit: a straight line through the primary.
            assert_close(
                binary_star::ecc_of_apprnt_orb(e, w, std::f64::consts::FRAC_PI_2),
                1.0,
                1e-7,
            );
        }
    }
    // A circular orbit seen at inclination i is an ellipse with axis ratio
    // cos i, so e' = sin i.
    for &i in &[10.0_f64, 45.0, 59.025, 80.0] {
        let i = i.to_radians();
        assert_close(binary_star::ecc_of_apprnt_orb(0.0, 1.0, i), i.sin(), 1e-12);
    }
}

// Meeus, chapter 57 exercise: ephemeris of γ Virginis (answers printed in
// the book with two decimals).
#[test]
fn gamma_vir_ephemeris_exercise() {
    let el = gamma_vir();
    let expected = [
        (1980.0, 296.65, 3.78),
        (1988.0, 288.70, 3.04),
        (1996.0, 274.41, 2.08),
        (2000.0, 259.34, 1.45),
        (2004.0, 208.67, 0.59),
        (2008.0, 35.54, 1.04),
        (2012.0, 12.72, 1.87),
    ];
    for &(year, theta, rho) in &expected {
        let pos = el.apparent_position(year);
        assert_close_deg(pos.position_angle.to_degrees(), theta, 0.005);
        assert_close(pos.separation, rho, 0.005);
    }
}

// The book quotes a least separation of 0.36" at 2005.21 for γ Vir.
#[test]
fn gamma_vir_least_separation() {
    let el = gamma_vir();
    let (mut best_t, mut best_rho) = (0.0, f64::MAX);
    let mut t = 2004.0;
    while t < 2007.0 {
        let rho = el.apparent_position(t).separation;
        if rho < best_rho {
            best_rho = rho;
            best_t = t;
        }
        t += 0.001;
    }
    assert_close(best_rho, 0.36, 0.005);
    assert_close(best_t, 2005.21, 0.01);
}

// With i between 90° and 180° the motion is retrograde: θ decreases.
#[test]
fn gamma_vir_position_angle_decreases() {
    let el = gamma_vir();
    let mut prev = el.apparent_position(1950.0).position_angle;
    let mut t = 1950.0;
    while t < 2100.0 {
        t += 0.25;
        let cur = el.apparent_position(t).position_angle;
        let step = (cur - prev).to_degrees().rem_euclid(360.0);
        assert!(step > 180.0, "θ increased at {t}");
        prev = cur;
    }
}

#[test]
fn apparent_position_is_periodic() {
    for el in [eta_crb(), gamma_vir()] {
        for &t in &[1500.0, 1934.008, 1980.0, 2005.13, 2300.5] {
            let a = el.apparent_position(t);
            let b = el.apparent_position(t + el.period);
            let c = el.apparent_position(t - 3.0 * el.period);
            assert_close_deg(
                a.position_angle.to_degrees(),
                b.position_angle.to_degrees(),
                1e-8,
            );
            assert_close(a.separation, b.separation, 1e-9);
            assert_close_deg(
                a.position_angle.to_degrees(),
                c.position_angle.to_degrees(),
                1e-8,
            );
            assert_close(a.separation, c.separation, 1e-9);
            assert!((0.0..std::f64::consts::TAU).contains(&a.position_angle));
        }
    }
}

#[test]
fn apparent_position_special_geometries() {
    let mut el = eta_crb();
    // At periastron the separation is a (1 - e) projected.
    el.inclination = 0.0;
    let p = el.apparent_position(el.time_of_periastron);
    assert_close(
        p.separation,
        el.semimajor_axis * (1.0 - el.eccentricity),
        1e-12,
    );
    // Face-on: θ = Ω + ω at periastron (v = 0).
    assert_close_deg(
        p.position_angle.to_degrees(),
        (el.asc_node_position_angle + el.periastron_longitude).to_degrees(),
        1e-9,
    );
    // Half a period later the companion is at apastron, opposite.
    let q = el.apparent_position(el.time_of_periastron + el.period / 2.0);
    assert_close(
        q.separation,
        el.semimajor_axis * (1.0 + el.eccentricity),
        1e-12,
    );
    assert_close_deg(
        q.position_angle.to_degrees(),
        p.position_angle.to_degrees() + 180.0,
        1e-9,
    );
    // A circular face-on orbit: constant separation, uniform motion.
    el.eccentricity = 0.0;
    for k in 0..8 {
        let t = el.time_of_periastron + el.period * k as f64 / 8.0;
        let pos = el.apparent_position(t);
        assert_close(pos.separation, el.semimajor_axis, 1e-12);
        assert_close_deg(
            pos.position_angle.to_degrees(),
            (el.asc_node_position_angle + el.periastron_longitude).to_degrees() + 45.0 * k as f64,
            1e-9,
        );
    }
}

#[test]
fn apparent_position_bad_eccentricity_gives_nan() {
    let mut el = eta_crb();
    el.eccentricity = 1.2;
    let p = el.apparent_position(1980.0);
    assert!(p.separation.is_nan());
    el.eccentricity = 0.5;
    assert!(el.apparent_position(f64::NAN).separation.is_nan());
}

// Very eccentric orbits near periastron, where a simple Kepler iteration
// fails: for a face-on orbit the position angle gives the true anomaly
// directly (θ = Ω + ω + v), so the result can be checked against Kepler's
// equation and the polar equation of the ellipse.
#[test]
fn apparent_position_very_eccentric_orbits() {
    let mut el = eta_crb();
    el.inclination = 0.0;
    for &e in &[0.9, 0.97, 0.999] {
        el.eccentricity = e;
        for &dt in &[-0.3, -0.01, -1e-4, 1e-4, 0.002, 0.05, 1.0, 17.0] {
            let t = el.time_of_periastron + dt;
            let p = el.apparent_position(t);
            let v = p.position_angle - el.asc_node_position_angle - el.periastron_longitude;
            // r = a (1 - e²) / (1 + e cos v)
            let r = el.semimajor_axis * (1.0 - e * e) / (1.0 + e * v.cos());
            assert_close(p.separation, r, 1e-12);
            // Back to the mean anomaly through E (Meeus, eqs. 30.1, 30.5).
            let ecc_anom = 2.0 * (((1.0 - e) / (1.0 + e)).sqrt() * (v / 2.0).tan()).atan();
            let m = ecc_anom - e * ecc_anom.sin();
            let m_expected = (std::f64::consts::TAU * dt / el.period + std::f64::consts::PI)
                .rem_euclid(std::f64::consts::TAU)
                - std::f64::consts::PI;
            assert_close(m, m_expected, 1e-9);
        }
    }
    // e = 1 is not an ellipse; negative e is meaningless.
    for &e in &[1.0, -0.1, f64::NAN] {
        el.eccentricity = e;
        let p = el.apparent_position(1980.0);
        assert!(
            p.separation.is_nan() && p.position_angle.is_nan(),
            "e = {e}"
        );
    }
}

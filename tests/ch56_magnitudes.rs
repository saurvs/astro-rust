//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 56: stellar magnitudes.

mod common;

use astro::star;
use common::assert_close;

fn arcsec_to_rad(s: f64) -> f64 {
    (s / 3600.0_f64).to_radians()
}

// Meeus, example 56.a: Castor, components of magnitude 1.96 and 2.89.
#[test]
fn combined_mag_castor_ex56a() {
    // Intermediate value x = 0.4 (m2 - m1) = 0.372
    assert_close(star::brightness_ratio(1.96, 2.89).log10(), 0.372, 1e-12);
    assert_close(star::combined_mag(1.96, 2.89), 1.58, 0.005);
}

#[test]
fn combined_mag_is_symmetric_and_brighter_than_parts() {
    for &(m1, m2) in &[(1.96, 2.89), (-1.46, 8.44), (5.0, 5.0), (12.0, -3.0)] {
        let m = star::combined_mag(m1, m2);
        assert_close(m, star::combined_mag(m2, m1), 1e-12);
        assert!(m < m1.min(m2));
        assert_close(m, star::combined_mag_of_many(&[m1, m2]), 1e-12);
    }
    // Two equal stars are 2.5 log10(2) ≈ 0.753 mag brighter than one.
    assert_close(
        star::combined_mag(4.0, 4.0),
        4.0 - 2.5 * 2_f64.log10(),
        1e-12,
    );
}

// Meeus, example 56.b: the triple star β Mon.
#[test]
fn combined_mag_beta_mon_ex56b() {
    assert_close(star::combined_mag_of_many(&[4.73, 5.22, 5.60]), 3.93, 0.005);
}

#[test]
fn combined_mag_of_many_edge_cases() {
    assert_close(star::combined_mag_of_many(&[3.3]), 3.3, 1e-12);
    assert_eq!(star::combined_mag_of_many(&[]), f64::INFINITY);
}

// Meeus, example 56.c: a star cluster. The book's sum is 0.15617; it then
// prints 0.15647 inside the logarithm (a typo), but its result +2.02 follows
// from 0.15617.
#[test]
fn combined_mag_cluster_ex56c() {
    let groups = [(4, 5.0), (14, 6.0), (23, 7.0), (38, 8.0)];
    let m = star::combined_mag_of_groups(&groups);
    // The book sums terms rounded to five decimals.
    assert_close(10_f64.powf(-0.4 * m), 0.15617, 0.00001);
    assert_close(m, 2.02, 0.005);

    // Same as listing every star.
    let mut all = Vec::new();
    for &(n, mag) in &groups {
        all.extend(std::iter::repeat_n(mag, n as usize));
    }
    assert_close(m, star::combined_mag_of_many(&all), 1e-12);

    assert_eq!(star::combined_mag_of_groups(&[]), f64::INFINITY);
    assert_eq!(star::combined_mag_of_groups(&[(0, 1.0)]), f64::INFINITY);
}

// Meeus, example 56.d: Vega (0.14) against Polaris (2.12).
#[test]
fn brightness_ratio_vega_polaris_ex56d() {
    assert_close(star::brightness_ratio(0.14, 2.12), 6.19, 0.005);
    // Polaris is fainter by the same factor.
    assert_close(
        star::brightness_ratio(2.12, 0.14),
        1.0 / star::brightness_ratio(0.14, 2.12),
        1e-12,
    );
    // Five magnitudes are exactly a factor 100.
    assert_close(star::brightness_ratio(1.0, 6.0), 100.0, 1e-9);
}

// Meeus, example 56.e: a star 500 times brighter than another.
#[test]
fn mag_diff_ex56e() {
    assert_close(star::mag_diff(500.0), 6.75, 0.005);
    // Round trip with brightness_ratio.
    for &(m1, m2) in &[(0.14, 2.12), (-26.7, 4.83), (3.0, 3.0)] {
        assert_close(
            star::mag_diff(star::brightness_ratio(m1, m2)),
            m2 - m1,
            1e-12,
        );
    }
}

#[test]
fn absolute_magnitude_from_parallax_and_distance() {
    // At 10 pc (parallax 0.1") the absolute magnitude equals the apparent one.
    assert_close(
        star::abs_mag_frm_parallax(arcsec_to_rad(0.1), 3.7),
        3.7,
        1e-12,
    );
    assert_close(star::abs_mag_frm_dist(10.0, 3.7), 3.7, 1e-12);
    // At 1 pc a star is 5 magnitudes brighter than at 10 pc.
    assert_close(star::abs_mag_frm_dist(1.0, 0.0), 5.0, 1e-12);
    assert_close(
        star::abs_mag_frm_parallax(arcsec_to_rad(1.0), 0.0),
        5.0,
        1e-12,
    );

    // Both forms agree through d = 1 / parallax.
    for &p in &[0.002, 0.0377, 0.1289, 0.7687] {
        let par = arcsec_to_rad(p);
        let d = star::dist_frm_parallax(par);
        assert_close(d, 1.0 / p, 1e-9 / p);
        assert_close(
            star::abs_mag_frm_parallax(par, 1.25),
            star::abs_mag_frm_dist(d, 1.25),
            1e-12,
        );
    }
}

#[test]
fn distance_from_parallax() {
    // Proxima Centauri, parallax 0.7687": 1.301 pc, 4.243 ly.
    let d = star::dist_frm_parallax(arcsec_to_rad(0.7687));
    assert_close(d, 1.3009, 0.0001);
    assert_close(d * astro::consts::PARSEC_IN_LIGHT_YEARS, 4.243, 0.001);
    assert_eq!(star::dist_frm_parallax(0.0), f64::INFINITY);
}

//! Unit tests for the private helpers of the planetary position code
//! (Meeus, chapters 32 and 33).

use super::*;

fn arcsec(rad: f64) -> f64 {
    rad.to_degrees() * 3600.0
}

/// Sums the `n` largest terms of one VSOP87 series with the amplitudes
/// rounded to 1e-8, as printed in Meeus' appendix III. Returns the sum in
/// units of 1e-8 (radian or AU), like the book.
fn abridged_sum(terms: &[[f64; 3]], n: usize, tau: f64) -> f64 {
    let mut sorted: Vec<[f64; 3]> = terms.to_vec();
    sorted.sort_by(|a, b| b[0].abs().total_cmp(&a[0].abs()));
    sorted
        .iter()
        .take(n)
        .map(|t| (t[0] * 1e8).round() * (t[1] + t[2] * tau).cos())
        .sum()
}

// Meeus, example 32.a: the abridged series of appendix III are the largest
// terms of the complete VSOP87 series held in VSOPD_87. Taking the same
// number of terms reproduces the book's partial sums and final values.
#[test]
fn appendix_iii_is_a_subset_of_the_full_series_example_32a() {
    let tau = time::julian_mill(2448976.5);
    assert!((tau - -0.007032169747).abs() < 1e-12);

    let venus = &VSOPD_87::venus::TERMS;

    // Number of terms of each Venus series in appendix III.
    let l_counts = [24, 12, 8, 3, 3, 1];
    let b_counts = [9, 4, 4, 4, 1];
    let r_counts = [12, 3, 3, 1, 1];

    let l_sums: Vec<f64> = (0..6)
        .map(|k| abridged_sum(venus[0][k], l_counts[k], tau))
        .collect();
    // Book: L0 = +316 402 122, L1 = +1 021 353 038 718, L2 = +50 055,
    // L3 = -56, L4 = -109, L5 = -1. The book rounds B and C too, hence
    // the tolerance of about one unit.
    let book = [316402122.0, 1021353038718.0, 50055.0, -56.0, -109.0, -1.0];
    for (sum, expected) in l_sums.iter().zip(book) {
        assert!((sum - expected).abs() < 1.5, "{sum} vs {expected}");
    }

    let combine = |sums: &[f64]| sums.iter().rev().fold(0.0, |acc, s| acc * tau + s) / 1e8;
    let l = combine(&l_sums);
    assert!((l - -68.6592582).abs() < 1e-7, "{l}");

    let b_sums: Vec<f64> = (0..5)
        .map(|k| abridged_sum(venus[1][k], b_counts[k], tau))
        .collect();
    let b = combine(&b_sums);
    assert!((b - -0.0457399).abs() < 1e-9, "{b}");

    let r_sums: Vec<f64> = (0..5)
        .map(|k| abridged_sum(venus[2][k], r_counts[k], tau))
        .collect();
    let r = combine(&r_sums);
    assert!((r - 0.724603).abs() < 1e-6, "{r}");
}

// Meeus, ch. 25: the complete VSOP87 series for the Earth have 2425 terms:
// 1080 in longitude, 348 in latitude and 997 in the radius vector.
#[test]
fn earth_term_counts_match_the_book() {
    let count = |v: usize| -> usize {
        VSOPD_87::earth::TERMS[v]
            .iter()
            .map(|series| series.len())
            .sum()
    };
    assert_eq!(count(0), 1080);
    assert_eq!(count(1), 348);
    assert_eq!(count(2), 997);
}

// Term counts of every series of version D of VSOP87, as distributed by
// the Bureau des Longitudes (file headers of VSOP87D.mer ... VSOP87D.nep).
// Guards against truncated or duplicated data.
#[test]
fn term_counts_of_all_series() {
    #[rustfmt::skip]
    let expected: [(Planet, [&[usize]; 3]); 8] = [
        (Planet::Mercury, [&[1380, 839, 395, 153, 28, 13], &[818, 494, 230, 53, 15, 10], &[1215, 711, 326, 119, 18, 10]]),
        (Planet::Venus, [&[367, 215, 70, 9, 5, 5], &[210, 133, 59, 15, 5, 4], &[330, 180, 63, 7, 3, 2]]),
        (Planet::Earth, [&[559, 341, 142, 22, 11, 5], &[184, 99, 49, 11, 5], &[526, 292, 139, 27, 10, 3]]),
        (Planet::Mars, [&[1217, 686, 310, 129, 36, 15], &[441, 287, 130, 41, 11, 5], &[1118, 596, 313, 111, 28, 9]]),
        (Planet::Jupiter, [&[760, 369, 191, 109, 45, 10], &[249, 141, 81, 42, 12, 5], &[745, 381, 190, 98, 46, 9]]),
        (Planet::Saturn, [&[1152, 642, 321, 148, 68, 27], &[500, 260, 111, 58, 26, 11], &[1205, 639, 342, 157, 64, 28]]),
        (Planet::Uranus, [&[947, 426, 151, 46, 7, 1], &[283, 154, 60, 16, 2], &[1124, 514, 192, 55, 11]]),
        (Planet::Neptune, [&[423, 183, 57, 15, 2, 1], &[172, 82, 25, 9, 1, 1], &[607, 250, 72, 22, 7]]),
    ];

    for (planet, counts) in expected {
        let series = vsop87d_series(&planet);
        for v in 0..3 {
            let actual: Vec<usize> = series[v].iter().map(|s| s.len()).collect();
            assert_eq!(actual, counts[v], "{planet:?}, variable {v}");
        }
    }
}

// Every term must be finite, and the phases must lie in [0, 2π) as in the
// original files; a malformed term would usually break one of these.
#[test]
fn terms_are_well_formed() {
    for planet in [
        Planet::Mercury,
        Planet::Venus,
        Planet::Earth,
        Planet::Mars,
        Planet::Jupiter,
        Planet::Saturn,
        Planet::Uranus,
        Planet::Neptune,
    ] {
        for var in vsop87d_series(&planet) {
            for series in var.iter() {
                for t in series.iter() {
                    assert!(t.iter().all(|x| x.is_finite()), "{planet:?} {t:?}");
                    assert!(t[0] >= 0.0, "{planet:?} {t:?}");
                    assert!(
                        (0.0..std::f64::consts::TAU).contains(&t[1]),
                        "{planet:?} {t:?}"
                    );
                    assert!(t[2] >= 0.0, "{planet:?} {t:?}");
                }
            }
        }
    }
}

// Meeus, example 33.a, with the book's heliocentric coordinates as input:
// first the geometric position (A), then the planet at t - τ (B).
#[test]
fn rect_coords_and_light_time_example_33a() {
    let (l0, b0, r0) = (
        88.35704_f64.to_radians(),
        0.00014_f64.to_radians(),
        0.983824,
    );

    let (x, y, z) = geocent_ecl_rect_coords(
        l0,
        b0,
        r0,
        26.11428_f64.to_radians(),
        -2.62070_f64.to_radians(),
        0.724603,
    );
    assert!((x - 0.621746).abs() < 1e-6);
    assert!((y - -0.664810).abs() < 1e-6);
    assert!((z - -0.033134).abs() < 1e-6);
    let delta = dist_frm_ecl_rect_coords(x, y, z);
    assert!((delta - 0.910845).abs() < 1e-6);
    assert!((light_time(delta) - 0.0052606).abs() < 1e-7);

    let (x, y, z) = geocent_ecl_rect_coords(
        l0,
        b0,
        r0,
        26.10588_f64.to_radians(),
        -2.62102_f64.to_radians(),
        0.724604,
    );
    assert!((x - 0.621794).abs() < 1e-6);
    assert!((y - -0.664905).abs() < 1e-6);
    assert!((z - -0.033138).abs() < 1e-6);
    let delta = dist_frm_ecl_rect_coords(x, y, z);
    assert!((delta - 0.910947).abs() < 1e-6);
    assert!((light_time(delta) - 0.0052612).abs() < 1e-7);

    let (lambda, beta) = ecl_coords_frm_ecl_rect_coords(x, y, z);
    assert!(
        (angle::limit_to_360(lambda.to_degrees()) - 313.08102).abs() < 2e-5,
        "{}",
        lambda.to_degrees()
    );
    // The inputs are rounded to 1e-5°, hence the tolerance.
    assert!(
        (beta.to_degrees() - -2.08474).abs() < 2e-5,
        "{}",
        beta.to_degrees()
    );
}

// Meeus, example 33.a: aberration of Venus by eq. 23.2, with the Sun's
// longitude L0 + 180° = 268.35704°. `geocent_apprnt_pos` uses
// `aberr::annual_aberr_in_ecl_coords` for this step.
#[test]
fn ecl_aberr_example_33a() {
    let (d_long, d_lat) = crate::aberr::annual_aberr_in_ecl_coords(
        313.08102_f64.to_radians(),
        -2.08474_f64.to_radians(),
        268.35704_f64.to_radians(),
        2448976.5,
    );
    assert!(
        (arcsec(d_long) - -14.868).abs() < 0.001,
        "{}",
        arcsec(d_long)
    );
    assert!((arcsec(d_lat) - -0.531).abs() < 0.001, "{}", arcsec(d_lat));
}

// The light-time iteration converges: the returned τ is consistent with
// the distance of the returned position, for all planets.
#[test]
fn light_time_iteration_converges() {
    for planet in [
        Planet::Mercury,
        Planet::Venus,
        Planet::Mars,
        Planet::Jupiter,
        Planet::Saturn,
        Planet::Uranus,
        Planet::Neptune,
    ] {
        for jd in [2415020.0, 2448976.5, 2460000.5, 2488070.0] {
            let earth = heliocent_coords(&Planet::Earth, jd);
            let (_, rect, tau) = light_time_corrected_rect(&planet, jd, earth);
            let (L, B, R) = heliocent_coords(&planet, jd - tau);
            let again = geocent_ecl_rect_coords(earth.0, earth.1, earth.2, L, B, R);
            let d = dist_frm_ecl_rect_coords(again.0, again.1, again.2);
            assert!((light_time(d) - tau).abs() < 1e-9, "{planet:?} {jd}");
            assert!((again.0 - rect.0).abs() < 1e-9);
            assert!((again.1 - rect.1).abs() < 1e-9);
            assert!((again.2 - rect.2).abs() < 1e-9);
        }
    }
}

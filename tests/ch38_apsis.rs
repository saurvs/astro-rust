// Tests for planet::apsis (Meeus, Astronomical Algorithms, 2nd ed.,
// chapter 38: planets in perihelion and in aphelion).

mod common;

use astro::planet::apsis::{self, Apsis, Extremum};
use astro::planet::{self, Planet};
use common::assert_close;

const HOUR: f64 = 1.0 / 24.0;

const ALL_PLANETS: [Planet; 8] = [
    Planet::Mercury,
    Planet::Venus,
    Planet::Earth,
    Planet::Mars,
    Planet::Jupiter,
    Planet::Saturn,
    Planet::Uranus,
    Planet::Neptune,
];

/// Julian Day of a Gregorian calendar date (Meeus, chapter 7)
fn jd(year: i64, month: i64, day: f64) -> f64 {
    let (y, m) = if month <= 2 {
        (year - 1, month + 12)
    } else {
        (year, month)
    };
    let a = (y as f64 / 100.0).floor();
    let b = 2.0 - a + (a / 4.0).floor();
    (365.25 * (y as f64 + 4716.0)).floor() + (30.6001 * (m as f64 + 1.0)).floor() + day + b - 1524.5
}

/// Asserts that `jde` falls on the given (Gregorian) calendar day.
#[track_caller]
fn assert_day(jde: f64, year: i64, month: i64, day: i64) {
    let start = jd(year, month, day as f64);
    assert!(
        jde >= start && jde < start + 1.0,
        "JDE {jde} is not on {year}-{month}-{day} (day starts at {start})"
    );
}

fn rad_vec(p: &Planet, jde: f64) -> f64 {
    planet::heliocent_coords(p, jde).2
}

// ---------------------------------------------------------------------------
// Worked examples
// ---------------------------------------------------------------------------

// Meeus, example 38.a: Venus' perihelion nearest to 1978 October 15
#[test]
fn example_38a_venus_perihelion() {
    let p = Planet::Venus;
    assert_close(apsis::approx_k(&p, 1978.79), -35.34, 0.005);
    let n = apsis::nearest_n(&p, Apsis::Perihelion, 1978.79).unwrap();
    assert_eq!(n, -35);
    let jde = apsis::time_of_passage(&p, Apsis::Perihelion, n);
    assert_close(jde, 2443873.704, 0.001);
    // 1978 December 31.204 = December 31, 5h TD
    assert_close(jde - jd(1978, 12, 31.0), 0.204, 0.001);

    // The same passage from VSOP87: within the expected few hours, at
    // Venus' perihelion distance a(1 - e) = 0.7184 AU.
    let passage = apsis::accurate_passage(&p, Apsis::Perihelion, n);
    assert!((passage.jde - jde).abs() < 0.1, "{}", passage.jde - jde);
    assert_close(passage.rad_vec, 0.718443, 1e-6);
}

// Meeus, example 38.b: the first aphelion of Mars after the beginning of
// 2032 (k = +16.5)
#[test]
fn example_38b_mars_aphelion() {
    let p = Planet::Mars;
    assert_close(apsis::approx_k(&p, 2032.0), 16.07, 0.005);
    let n = apsis::nearest_n(&p, Apsis::Aphelion, 2032.0).unwrap();
    assert_eq!(n, 16);
    let jde = apsis::time_of_passage(&p, Apsis::Aphelion, n);
    assert_close(jde, 2463530.456, 0.001);
    // 2032 October 24.956 = October 24, 23h TD
    assert_close(jde - jd(2032, 10, 24.0), 0.956, 0.001);

    let passage = apsis::accurate_passage(&p, Apsis::Aphelion, n);
    assert!((passage.jde - jde).abs() < 0.25);
    assert_close(passage.rad_vec, 1.665979, 1e-6);
}

// Meeus, chapter 38: k = +14 and k = -222 are perihelia, k = +27.5 and
// k = -119.5 aphelia; with the numbering of this crate those are perihelia
// 14 and -222 and aphelia 27 and -120.
#[test]
fn numbering_of_passages() {
    for p in &ALL_PLANETS {
        let t = |k: f64| {
            let jde0 = apsis::time_of_passage(p, Apsis::Perihelion, 0);
            let b = apsis::time_of_passage(p, Apsis::Perihelion, 1) - jde0;
            let c = (apsis::time_of_passage(p, Apsis::Perihelion, 1)
                + apsis::time_of_passage(p, Apsis::Perihelion, -1)
                - 2.0 * jde0)
                / 2.0;
            jde0 + (b - c) * k + c * k * k
        };
        for (apsis, n, k) in [
            (Apsis::Perihelion, 14, 14.0),
            (Apsis::Perihelion, -222, -222.0),
            (Apsis::Aphelion, 27, 27.5),
            (Apsis::Aphelion, -120, -119.5),
        ] {
            let jde = apsis::time_of_passage(p, apsis, n);
            assert_close(jde, t(k), 1e-3);
        }
    }
}

// Meeus, chapter 38: the Keplerian formula for Jupiter with k = -2.5 gives
// 1981 July 19, but the aphelion occurred on 1981 July 28.
#[test]
fn jupiter_aphelion_1981() {
    let kepler = apsis::time_of_passage(&Planet::Jupiter, Apsis::Aphelion, -3);
    assert_day(kepler, 1981, 7, 19);
    let passage = apsis::accurate_passage(&Planet::Jupiter, Apsis::Aphelion, -3);
    assert_day(passage.jde, 1981, 7, 28);
}

// Meeus, chapter 38: the Keplerian formula for Saturn with k = -2 gives
// 1944 July 30, but the perihelion occurred on 1944 September 8.
#[test]
fn saturn_perihelion_1944() {
    let kepler = apsis::time_of_passage(&Planet::Saturn, Apsis::Perihelion, -2);
    assert_day(kepler, 1944, 7, 30);
    let passage = apsis::accurate_passage(&Planet::Saturn, Apsis::Perihelion, -2);
    assert_day(passage.jde, 1944, 9, 8);
    assert_close(passage.rad_vec, 9.0288, 0.00005);
}

// Meeus, chapter 38: the Earth's perihelion of 1990 (k = -10)
#[test]
fn earth_perihelion_1990() {
    let barycentre = apsis::time_of_passage(&Planet::Earth, Apsis::Perihelion, -10);
    assert_close(barycentre, 2447894.911, 0.001);
    assert_close(barycentre - jd(1990, 1, 3.0), 0.41, 0.005);
    assert_close(
        apsis::earth_correction(Apsis::Perihelion, -10),
        1.261,
        0.0005,
    );
    let earth = apsis::earth_time_of_passage(Apsis::Perihelion, -10);
    assert_close(earth, 2447896.172, 0.001);
    // "1990 January 4, at 16h TD"
    assert!(((earth - jd(1990, 1, 4.0)) / HOUR - 16.0).abs() < 0.5);
    // The correct instant: 1990 January 4 at 17h TD.
    let passage = apsis::accurate_passage(&Planet::Earth, Apsis::Perihelion, -10);
    assert!(((passage.jde - jd(1990, 1, 4.0)) / HOUR - 17.0).abs() < 0.5);
}

// Spot checks of Meeus' table 38.A (Earth, from VSOP87; instants in TD to
// 0.01 hour, distances to 1e-6 AU).
#[test]
fn earth_table_38a() {
    // (n, perihelion (month, day, hour), R, aphelion (month, day, hour), R)
    let cases = [
        (-9, (1, 3, 3.00), 0.983281, (7, 6, 15.46), 1.016703),
        (0, (1, 3, 5.31), 0.983321, (7, 3, 23.84), 1.016741),
        (1, (1, 4, 8.89), 0.983286, (7, 4, 13.65), 1.016643),
        (10, (1, 3, 0.18), 0.983290, (7, 6, 11.52), 1.016702),
    ];
    for (n, (pm, pd, ph), pr, (am, ad, ah), ar) in cases {
        let year = 2000 + n;
        let p = apsis::accurate_passage(&Planet::Earth, Apsis::Perihelion, n);
        let hours = (p.jde - jd(year, pm, pd as f64)) / HOUR;
        assert_close(hours, ph, 0.015);
        assert_close(p.rad_vec, pr, 0.6e-6);
        let a = apsis::accurate_passage(&Planet::Earth, Apsis::Aphelion, n);
        let hours = (a.jde - jd(year, am, ad as f64)) / HOUR;
        assert_close(hours, ah, 0.015);
        assert_close(a.rad_vec, ar, 0.6e-6);
    }
}

// Meeus, chapter 38: the corrected instants for the Earth have a mean error
// of about 3 hours over 1980-2019, 6 hours at most.
#[test]
fn earth_correction_accuracy() {
    for apsis in [Apsis::Perihelion, Apsis::Aphelion] {
        let mut sum = 0.0;
        let mut max: f64 = 0.0;
        let mut barycentre_max: f64 = 0.0;
        for n in -20..20 {
            let accurate = apsis::accurate_passage(&Planet::Earth, apsis, n).jde;
            let err = (apsis::earth_time_of_passage(apsis, n) - accurate).abs() / HOUR;
            sum += err;
            max = max.max(err);
            let err = (apsis::time_of_passage(&Planet::Earth, apsis, n) - accurate).abs() / HOUR;
            barycentre_max = barycentre_max.max(err);
        }
        assert!(sum / 40.0 < 3.0, "{apsis:?} mean {}", sum / 40.0);
        assert!(max < 7.0, "{apsis:?} max {max}");
        // Without the correction the error exceeds a day.
        assert!(barycentre_max > 24.0);
    }
}

// Meeus' table of the extreme distances of Saturn (1920-2050) and Uranus
// (1750-2100), computed with VSOP87; a few spot checks of each apsis.
#[test]
fn saturn_and_uranus_extremes() {
    // (planet, apsis, n, date, distance)
    let cases = [
        (Planet::Saturn, Apsis::Aphelion, -3, (1929, 11, 11), 10.0467),
        (Planet::Saturn, Apsis::Perihelion, 1, (2032, 11, 28), 9.0149),
        (Planet::Uranus, Apsis::Perihelion, -3, (1798, 3, 3), 18.2890),
        (Planet::Uranus, Apsis::Aphelion, 0, (2092, 11, 23), 20.0994),
    ];
    for (p, apsis, n, (y, m, d), r) in cases {
        let passage = apsis::accurate_passage(&p, apsis, n);
        assert_day(passage.jde, y, m, d);
        assert_close(passage.rad_vec, r, 0.00005);
    }
}

// Meeus, chapter 38: Neptune's distance has double minima and maxima near
// its apsides.
#[test]
fn neptune_double_extrema() {
    use Extremum::{Maximum, Minimum};
    // (start, end, expected extrema as (kind, date, distance, tolerance))
    let windows = [
        (
            jd(1873, 1, 1.0),
            jd(1890, 1, 1.0),
            vec![
                (Minimum, (1876, 8, 28), 29.8148, 0.00005),
                (Maximum, (1881, 12, 12), 29.8213, 0.00005),
                (Minimum, (1886, 7, 11), 29.8174, 0.00005),
            ],
        ),
        (
            jd(1955, 1, 1.0),
            jd(1972, 1, 1.0),
            vec![
                (Maximum, (1959, 7, 13), 30.3317, 0.00005),
                (Minimum, (1965, 10, 6), 30.3227, 0.00005),
                (Maximum, (1968, 11, 21), 30.3241, 0.00005),
            ],
        ),
        (
            jd(2038, 1, 1.0),
            jd(2054, 1, 1.0),
            vec![
                (Minimum, (2042, 9, 5), 29.8064, 0.00005),
                (Maximum, (2049, 10, 24), 29.816711, 0.0000005),
                (Minimum, (2050, 6, 25), 29.816696, 0.0000005),
            ],
        ),
    ];
    for (start, end, expected) in windows {
        let found = apsis::rad_vec_extrema(&Planet::Neptune, start, end);
        assert_eq!(found.len(), expected.len(), "{found:?}");
        for (f, (kind, (y, m, d), r, tol)) in found.iter().zip(expected) {
            assert_eq!(f.kind, kind);
            // The extrema are so flat that the day is uncertain by one.
            let day = jd(y, m, d as f64);
            assert!((f.jde - (day + 0.5)).abs() < 1.5, "{f:?} vs {y}-{m}-{d}");
            assert_close(f.rad_vec, r, tol);
        }
    }

    // accurate_passage picks the principal extremum: the 1876 perihelion
    // (not 1886), the 1959 aphelion (not 1968), the 2042 perihelion.
    let p = apsis::accurate_passage(&Planet::Neptune, Apsis::Perihelion, -1);
    assert!((p.jde - jd(1876, 8, 28.5)).abs() < 1.5);
    let a = apsis::accurate_passage(&Planet::Neptune, Apsis::Aphelion, -1);
    assert!((a.jde - jd(1959, 7, 13.5)).abs() < 1.5);
    let p = apsis::accurate_passage(&Planet::Neptune, Apsis::Perihelion, 0);
    assert!((p.jde - jd(2042, 9, 5.5)).abs() < 1.5);
    assert_close(p.rad_vec, 29.8064, 0.00005);
}

// ---------------------------------------------------------------------------
// Accuracy of the Keplerian formulas, and invariants
// ---------------------------------------------------------------------------

// Meeus: a few hours for Mars, up to about half a month for Jupiter, more
// than a month for Saturn; Uranus and Neptune worse still.
#[test]
fn keplerian_formula_errors() {
    // (planet, passages tested, max error allowed in days, error that must
    // be exceeded somewhere in days)
    let cases = [
        (Planet::Mercury, -20..20, 0.01, 0.0),
        (Planet::Venus, -20..20, 0.3, 0.0),
        (Planet::Mars, -20..20, 0.3, 0.05),
        (Planet::Jupiter, -6..6, 17.0, 5.0),
        (Planet::Saturn, -3..3, 120.0, 30.0),
        (Planet::Uranus, -2..2, 400.0, 100.0),
        (Planet::Neptune, -1..1, 2300.0, 500.0),
    ];
    for (p, range, max_err, min_max_err) in cases {
        let mut largest: f64 = 0.0;
        for n in range {
            for apsis in [Apsis::Perihelion, Apsis::Aphelion] {
                let kepler = apsis::time_of_passage(&p, apsis, n);
                let accurate = apsis::accurate_passage(&p, apsis, n).jde;
                largest = largest.max((kepler - accurate).abs());
            }
        }
        assert!(largest <= max_err, "{largest}");
        assert!(largest >= min_max_err, "{largest}");
    }
}

// The accurate passage really is an extremum of the VSOP87 distance, and
// the perihelion and aphelion distances are close to a(1 - e) and a(1 + e).
#[test]
fn accurate_passages_are_extrema() {
    // (planet, a, e, distance tolerance)
    let cases = [
        (Planet::Mercury, 0.387098, 0.205632, 0.0001),
        (Planet::Venus, 0.723330, 0.006772, 0.0002),
        (Planet::Earth, 1.000001, 0.016709, 0.0002),
        (Planet::Mars, 1.523679, 0.093401, 0.002),
        (Planet::Jupiter, 5.202603, 0.048498, 0.01),
        (Planet::Saturn, 9.554909, 0.055548, 0.05),
        (Planet::Uranus, 19.218446, 0.046381, 0.1),
        (Planet::Neptune, 30.110387, 0.009456, 0.1),
    ];
    for (p, a, e, tol) in cases {
        for n in [-1, 0, 1] {
            for (apsis, sign) in [(Apsis::Perihelion, 1.0), (Apsis::Aphelion, -1.0)] {
                let passage = apsis::accurate_passage(&p, apsis, n);
                assert_close(passage.rad_vec, rad_vec(&p, passage.jde), 1e-15);
                assert_close(passage.rad_vec, a * (1.0 - sign * e), tol);
                // A small step either way moves away from the extremum.
                let h = 0.5 + kepler_period(&p) * 1e-4;
                let r = passage.rad_vec;
                assert!(sign * (rad_vec(&p, passage.jde - h) - r) > 0.0);
                assert!(sign * (rad_vec(&p, passage.jde + h) - r) > 0.0);
                // And the formula's instant is not better.
                let kepler = apsis::time_of_passage(&p, apsis, n);
                assert!(sign * (rad_vec(&p, kepler) - r) >= 0.0);
            }
        }
    }
}

fn kepler_period(p: &Planet) -> f64 {
    apsis::time_of_passage(p, Apsis::Perihelion, 1)
        - apsis::time_of_passage(p, Apsis::Perihelion, 0)
}

// Perihelion n, aphelion n and perihelion n + 1 are in order, half a period
// apart, both for the formulas and for VSOP87.
#[test]
fn order_and_spacing() {
    for p in &ALL_PLANETS {
        let period = kepler_period(p);
        for n in [-3, -1, 0, 2] {
            let p0 = apsis::time_of_passage(p, Apsis::Perihelion, n);
            let a0 = apsis::time_of_passage(p, Apsis::Aphelion, n);
            let p1 = apsis::time_of_passage(p, Apsis::Perihelion, n + 1);
            assert!(p0 < a0 && a0 < p1);
            assert_close(a0 - p0, period / 2.0, period * 1e-4);
        }
        let p0 = apsis::accurate_passage(p, Apsis::Perihelion, 0).jde;
        let a0 = apsis::accurate_passage(p, Apsis::Aphelion, 0).jde;
        let p1 = apsis::accurate_passage(p, Apsis::Perihelion, 1).jde;
        assert!(p0 < a0 && a0 < p1);
        assert!(((p1 - p0) - period).abs() < 0.1 * period);
    }
}

// The book's approximate k formulas, and perihelion 0 being the first after
// the beginning of 2000.
#[test]
fn approximate_k_and_nearest_n() {
    let jan1 = jd(2000, 1, 1.0);
    for p in &ALL_PLANETS {
        assert!(apsis::time_of_passage(p, Apsis::Perihelion, 0) >= jan1);
        assert!(apsis::time_of_passage(p, Apsis::Perihelion, -1) < jan1);
        // approx_k is linear in the year and agrees with the formula
        let period_years = kepler_period(p) / 365.25;
        let dk = apsis::approx_k(p, 2001.0) - apsis::approx_k(p, 2000.0);
        assert_close(dk * period_years, 1.0, 0.001);
        for year in [-2000.5, -1.0, 1600.25, 1999.9, 2000.0, 2026.9, 3000.0] {
            for ap in [Apsis::Perihelion, Apsis::Aphelion] {
                let n = apsis::nearest_n(p, ap, year).unwrap();
                let target = 2451545.0 + (year - 2000.0) * 365.25;
                let d = |n| (apsis::time_of_passage(p, ap, n) - target).abs();
                assert!(d(n) <= d(n - 1) && d(n) <= d(n + 1), "{ap:?} {year}");
            }
        }
    }
}

// A non-finite year, or one so remote that n cannot be represented exactly,
// has no nearest passage (it used to give the valid-looking n = 0, or
// i64::MAX).
#[test]
fn nearest_n_rejects_unusable_years() {
    for year in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1e20, -1e300] {
        for ap in [Apsis::Perihelion, Apsis::Aphelion] {
            assert_eq!(apsis::nearest_n(&Planet::Mars, ap, year), None, "{year}");
        }
    }
    // Far, but still exact
    let n = apsis::nearest_n(&Planet::Neptune, Apsis::Perihelion, 1e9).unwrap();
    assert!((n as f64 - 0.00607 * (1e9 - 2047.5)).abs() < 2.0);
}

// The quadratic coefficients of the formulas follow from the secular terms
// of the mean elements (Meeus, chapter 31, mean equinox of the date): at the
// formula's instants the mean anomaly L - ϖ is a multiple of 360°
// (perihelion) or 180° (aphelion). This checks in particular the sign of
// the k² terms: at the largest n tested, a reversed sign would put the
// mean anomaly 0.4 days or more of mean motion away, while the cubic terms
// of the elements (ignored by the formulas) and the rounding of the
// formulas' constants stay below 0.2 day.
#[test]
fn formulas_agree_with_mean_elements() {
    // (planet, L coefficients, ϖ coefficients in degrees and Julian
    // centuries, largest n tested)
    let elements: [(Planet, [f64; 4], [f64; 4], i64); 4] = [
        (
            Planet::Jupiter,
            [34.351519, 3036.3027748, 0.0002233, 0.000000037],
            [14.331207, 1.6126352, 0.00103042, -0.000004464],
            40,
        ),
        (
            Planet::Saturn,
            [50.077444, 1223.5110686, 0.00051908, -0.00000003],
            [93.057237, 1.9637613, 0.00083753, 0.000004928],
            20,
        ),
        (
            Planet::Uranus,
            [314.055005, 429.8640561, 0.0003039, -0.000000026],
            [173.005291, 1.486379, 0.00021406, 0.000000434],
            6,
        ),
        (
            Planet::Neptune,
            [304.348665, 219.8833092, 0.00030882, 0.000000018],
            [48.120276, 1.4262957, 0.00038434, 0.00000002],
            3,
        ),
    ];
    for (p, l, w, n_max) in elements {
        let mean_anomaly = |jde: f64| {
            let t = (jde - 2451545.0) / 36525.0;
            let poly = |c: [f64; 4]| c[0] + t * (c[1] + t * (c[2] + t * c[3]));
            (poly(l) - poly(w)).rem_euclid(360.0)
        };
        let rate = (l[1] - w[1]) / 36525.0; // degrees per day
        for n in [-n_max, -n_max / 2, 0, n_max / 2, n_max] {
            for (apsis, target) in [(Apsis::Perihelion, 0.0), (Apsis::Aphelion, 180.0)] {
                let jde = apsis::time_of_passage(&p, apsis, n);
                let m = mean_anomaly(jde);
                let diff = (m - target + 180.0).rem_euclid(360.0) - 180.0;
                // Express the residual in days of mean motion.
                assert!(
                    (diff / rate).abs() < 0.2,
                    "{n} {apsis:?}: {} d",
                    diff / rate
                );
            }
        }
    }
}

#[test]
fn rad_vec_extrema_edge_cases() {
    let p = Planet::Earth;
    assert!(apsis::rad_vec_extrema(&p, 2451545.0, 2451545.0).is_empty());
    assert!(apsis::rad_vec_extrema(&p, 2451600.0, 2451545.0).is_empty());
    assert!(apsis::rad_vec_extrema(&p, f64::NAN, 2451545.0).is_empty());
    assert!(apsis::rad_vec_extrema(&p, 2451545.0, f64::INFINITY).is_empty());

    // One year of the Earth: a perihelion and an aphelion, matching
    // accurate_passage.
    let found = apsis::rad_vec_extrema(&p, jd(2000, 1, 1.0), jd(2000, 12, 31.0));
    assert_eq!(found.len(), 2, "{found:?}");
    assert_eq!(found[0].kind, Extremum::Minimum);
    assert_eq!(found[1].kind, Extremum::Maximum);
    let peri = apsis::accurate_passage(&p, Apsis::Perihelion, 0);
    let aph = apsis::accurate_passage(&p, Apsis::Aphelion, 0);
    assert_close(found[0].jde, peri.jde, 1e-4);
    assert_close(found[1].jde, aph.jde, 1e-4);
    assert_close(found[0].rad_vec, peri.rad_vec, 1e-12);

    // Mercury over one Earth year: four or five of each, alternating.
    let found = apsis::rad_vec_extrema(&Planet::Mercury, jd(2000, 1, 1.0), jd(2001, 1, 1.0));
    assert!(found.len() >= 8 && found.len() <= 9, "{}", found.len());
    assert!(found.windows(2).all(|w| w[0].kind != w[1].kind));
    assert!(found.windows(2).all(|w| w[0].jde < w[1].jde));
}

#[test]
fn far_from_2000() {
    // Large cycle numbers must not lose the Earth correction to rounding.
    for n in [-4000, -1000, 1000, 2000] {
        let c = apsis::earth_correction(Apsis::Perihelion, n);
        assert!(c.is_finite() && c.abs() < 1.6, "{c}");
        let c = apsis::earth_correction(Apsis::Aphelion, n);
        assert!(c.is_finite() && c.abs() < 1.6, "{c}");
    }
    // Earth around -1000 and +3000: the correction still helps.
    for n in [-3000, 1000] {
        let accurate = apsis::accurate_passage(&Planet::Earth, Apsis::Perihelion, n).jde;
        let corrected = apsis::earth_time_of_passage(Apsis::Perihelion, n);
        assert!((corrected - accurate).abs() < 0.5, "{n}");
    }
}

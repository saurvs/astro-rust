//! Meeus, chapter 10: Dynamical Time and Universal Time (ΔT).

mod common;

use astro::time::{self, CalType, Date, Month};

// Meeus, example 10.a: New Moon of 1977 February 18 at 3h37m40s TD, with
// ΔT = 48 s, happened at 3h36m52s UT.
#[test]
fn example_10a_td_to_ut() {
    let td = |h: f64, m: f64, s: f64| 18.0 + (h + m / 60.0 + s / 3600.0) / 24.0;
    let jde = time::julian_day(&Date {
        year: 1977,
        month: Month::Feb,
        decimal_day: td(3.0, 37.0, 40.0),
        cal_type: CalType::Gregorian,
    });
    let jd = time::julian_day_from_ephemeris_day(jde, 48.0);
    let expected = time::julian_day(&Date {
        year: 1977,
        month: Month::Feb,
        decimal_day: td(3.0, 36.0, 52.0),
        cal_type: CalType::Gregorian,
    });
    common::assert_close(jd, expected, 1e-8);
    common::assert_close(time::julian_ephemeris_day(jd, 48.0), jde, 1e-9);
}

// Meeus, example 10.b: 333 February 6 at 6h UT; ΔT = +6146 s from eq. 10.1,
// so the instant is 7h42m TD.
#[test]
fn example_10b_chapront_francou() {
    let dt = time::delta_t_chapront_francou(333.1).unwrap();
    common::assert_close(dt, 6146.0, 0.5);
    common::assert_close(dt / 60.0, 102.0, 0.5);
    let ut = time::julian_day(&Date {
        year: 333,
        month: Month::Feb,
        decimal_day: 6.25,
        cal_type: CalType::Julian,
    });
    let td = Date::from_julian_day(time::julian_ephemeris_day(ut, dt)).unwrap();
    let minutes = (td.decimal_day - 6.0) * 1440.0;
    assert_eq!(td.decimal_day.floor(), 6.0);
    common::assert_close(minutes, 7.0 * 60.0 + 42.0, 0.5);
}

#[test]
fn chapront_francou_pieces() {
    // eq. 10.1 and 10.2 nearly agree at 948
    let a = time::delta_t_chapront_francou(947.999_999).unwrap();
    let b = time::delta_t_chapront_francou(948.0).unwrap();
    common::assert_close(a, b, 0.5);
    // the correction removes the jump at 2000 (ΔT ≈ 65 s there) and
    // vanishes at 2100
    common::assert_close(time::delta_t_chapront_francou(2000.0).unwrap(), 65.0, 1e-9);
    let at_2100 = time::delta_t_chapront_francou(2100.0).unwrap();
    let after = time::delta_t_chapront_francou(2100.000001).unwrap();
    common::assert_close(at_2100, 229.3, 1e-9);
    common::assert_close(after, at_2100, 1e-3);
    // not defined over the span of the observed table
    assert!(time::delta_t_chapront_francou(1600.0).is_some());
    assert!(time::delta_t_chapront_francou(1600.1).is_none());
    assert!(time::delta_t_chapront_francou(1850.0).is_none());
    assert!(time::delta_t_chapront_francou(1999.9).is_none());
    assert!(time::delta_t_chapront_francou(f64::NAN).is_none());
    // positive in the remote past and in the distant future
    assert!(time::delta_t_chapront_francou(-4000.0).unwrap() > 0.0);
    assert!(time::delta_t_chapront_francou(4000.0).unwrap() > 0.0);
}

// A few observed values of ΔT at the beginning of the year (Meeus, table
// 10.A, from the Astronomical Almanac).
const OBSERVED: [(f64, f64); 12] = [
    (1800.0, 13.1),
    (1820.0, 11.6),
    (1840.0, 5.4),
    (1860.0, 7.7),
    (1880.0, -5.5),
    (1900.0, -2.8),
    (1910.0, 10.4),
    (1930.0, 24.0),
    (1950.0, 29.1),
    (1960.0, 33.1),
    (1980.0, 50.5),
    (1990.0, 56.9),
];

#[test]
fn espenak_meeus_matches_observed_values() {
    for (year, dt) in OBSERVED {
        common::assert_close(time::delta_t_from_decimal_year(year), dt, 0.7);
    }
}

// Values of Morrison & Stephenson (2004) at round centuries, as tabulated by
// Espenak & Meeus; the polynomials fit them to within a few seconds.
#[test]
fn espenak_meeus_matches_morrison_stephenson() {
    for (year, dt) in [
        (-400.0, 15530.0),
        (0.0, 10580.0),
        (500.0, 5710.0),
        (1000.0, 1570.0),
        (1500.0, 200.0),
        (1600.0, 120.0),
        (1700.0, 9.0),
    ] {
        common::assert_close(time::delta_t_from_decimal_year(year), dt, 5.0);
    }
}

#[test]
fn delta_t_uses_middle_of_month() {
    for (year, month) in [(-1000, 1), (0, 6), (1620, 12), (1988, 3), (2024, 7)] {
        let y = year as f64 + (month as f64 - 0.5) / 12.0;
        assert_eq!(
            time::delta_t(year, month),
            time::delta_t_from_decimal_year(y)
        );
    }
    // the value used by the transit test (Meeus, example 15.a)
    common::assert_close(time::delta_t(1988, 3), 55.9, 0.3);
}

#[test]
fn delta_t_long_term_parabola_outside_polynomials() {
    let parabola = |y: f64| -20.0 + 32.0 * ((y - 1820.0) / 100.0).powi(2);
    for y in [-3000.0, -1000.0, -500.1, 2150.0, 2500.0, 3000.0, 10000.0] {
        common::assert_close(time::delta_t_from_decimal_year(y), parabola(y), 1e-9);
    }
    assert!(time::delta_t_from_decimal_year(f64::NAN).is_nan());
}

// ΔT was slightly negative from about 1871 to 1901 (Meeus, chapter 10).
#[test]
fn delta_t_negative_around_1871_1901() {
    for year in 1872..=1901 {
        assert!(time::delta_t(year, 6) < 0.0, "{year}");
    }
    assert!(time::delta_t(1869, 6) > 0.0);
    assert!(time::delta_t(1903, 6) > 0.0);
}

#[test]
fn delta_t_smooth_everywhere() {
    // No step larger than what ΔT can change by in a month.
    let mut y = -2000.0;
    let mut prev = time::delta_t_from_decimal_year(y);
    while y < 3000.0 {
        y += 1.0 / 12.0;
        let v = time::delta_t_from_decimal_year(y);
        assert!((v - prev).abs() < 3.0, "jump at {y}: {prev} -> {v}");
        prev = v;
    }
}

// Meeus' polynomials for 1800-1997 (max. error 2.3 s), 1800-1899 and
// 1900-1997 (max. error 0.9 s), against observed values.
#[test]
fn meeus_polynomials_match_observed_values() {
    for (year, dt) in OBSERVED {
        common::assert_close(time::delta_t_1800_1997(year).unwrap(), dt, 2.3);
        if year < 1900.0 {
            common::assert_close(time::delta_t_1800_1899(year).unwrap(), dt, 0.9);
        } else {
            common::assert_close(time::delta_t_1900_1997(year).unwrap(), dt, 0.9);
        }
    }
}

#[test]
fn meeus_polynomials_validity_ranges() {
    assert!(time::delta_t_1800_1997(1799.99).is_none());
    assert!(time::delta_t_1800_1997(1800.0).is_some());
    assert!(time::delta_t_1800_1997(1997.99).is_some());
    assert!(time::delta_t_1800_1997(1998.0).is_none());
    // Meeus warns the 1800-1899 expression gives about 70000 s in 1945.
    assert!(time::delta_t_1800_1899(1945.0).is_none());
    assert!(time::delta_t_1800_1899(1899.99).is_some());
    assert!(time::delta_t_1900_1997(1899.99).is_none());
    assert!(time::delta_t_1900_1997(1900.0).is_some());
    assert!(time::delta_t_1900_1997(f64::NAN).is_none());
    // the three agree with the more accurate NASA polynomials
    let mut y = 1800.0;
    while y < 1998.0 {
        let nasa = time::delta_t_from_decimal_year(y);
        common::assert_close(time::delta_t_1800_1997(y).unwrap(), nasa, 3.5);
        let short = time::delta_t_1800_1899(y).or(time::delta_t_1900_1997(y));
        common::assert_close(short.unwrap(), nasa, 2.0);
        y += 0.25;
    }
}

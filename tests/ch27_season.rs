//! Meeus, chapter 27: equinoxes and solstices.

mod common;

use astro::season::{self, Event, Season};
use astro::*;
use common::*;

use std::f64::consts::{PI, TAU};

fn wrap_pi(a: f64) -> f64 {
    (a + PI).rem_euclid(TAU) - PI
}

/// Julian Day of a Gregorian calendar date at hours, minutes and seconds
/// (Meeus, chapter 7; local copy so that these tests only depend on the
/// season module).
fn jd(year: i32, month: i32, day: u8, h: f64, m: f64, s: f64) -> f64 {
    let (y, mo) = if month <= 2 {
        (year - 1, month + 12)
    } else {
        (year, month)
    };
    let a = (y as f64 / 100.0).floor();
    let b = 2.0 - a + (a / 4.0).floor();
    let d = day as f64 + (h + m / 60.0 + s / 3600.0) / 24.0;
    (365.25 * (y as f64 + 4716.0)).floor() + (30.6001 * (mo as f64 + 1.0)).floor() + d + b - 1524.5
}

const SECOND: f64 = 1.0 / 86400.0;

// Meeus, example 27.a: the June solstice of 1962 from the mean instant and
// the periodic terms
#[test]
fn example_27a() {
    assert_close(
        season::mean_instant(1962, Event::JuneSolstice),
        2437837.38589,
        5e-6,
    );
    assert_close(
        season::approx_instant(1962, Event::JuneSolstice),
        2437837.39245,
        5e-6,
    );
}

// Meeus, example 27.b: the first iteration of eq. 27.1 from the mean
// instant. The book uses the abridged VSOP87 of Appendix III, the crate
// the complete theory, so λ differs by about 0.07".
#[test]
fn example_27b_first_iteration() {
    let jde0 = 2437837.38589;
    let p = sun::apparent_position(jde0);

    assert_close(p.radius_vector, 1.0163018, 1e-6);
    assert_close(
        p.nutation_in_longitude.to_degrees() * 3600.0,
        -12.965,
        0.0005,
    );
    // aberration by eq. 25.10 for this R
    assert_close(
        aberr::sol_aberr(1.0163018).to_degrees() * 3600.0,
        -20.161,
        0.0005,
    );

    assert_close_deg(p.apparent_longitude.to_degrees(), 89.994045, 0.0001);
    let correction = 58.0 * (PI / 2.0 - p.apparent_longitude).sin();
    assert_close(correction, 0.00603, 0.00005);
}

// Meeus, example 27.b: the instant of the June solstice of 1962. The book
// gives 21h24m42s TD as the correct instant from the complete VSOP87.
#[test]
fn example_27b_instant() {
    let jde = season::instant(1962, Event::JuneSolstice);
    let expected = jd(1962, 6, 21, 21.0, 24.0, 42.0);
    assert_close(jde, expected, 1.5 * SECOND);

    // the abridged series of the book gives 21h24m40s
    assert_close(jde, 2437837.39213, 3.0 * SECOND);
}

// At the computed instant the apparent longitude of the Sun is k·90°.
#[test]
fn apparent_longitude_at_instant() {
    for year in [-1500, -500, 0, 1000, 1582, 1962, 2000, 2024, 3000, 5000] {
        for e in Event::ALL {
            let jde = season::instant(year, e);
            let lambda = sun::apparent_position(jde).apparent_longitude;
            assert!(
                wrap_pi(lambda - e.solar_longitude()).abs() < arcsec(0.01).to_radians(),
                "{year} {e:?}: λ = {}°",
                lambda.to_degrees()
            );
        }
    }
}

// Published instants (U.S. Naval Observatory, rounded to the minute, in UT)
// converted to TD with ΔT = 63.8 s (2000) and 69.2 s (2024).
#[test]
fn published_instants() {
    let cases = [
        (
            2000,
            Event::MarchEquinox,
            jd(2000, 3, 20, 7.0, 35.0, 0.0),
            63.8,
        ),
        (
            2000,
            Event::JuneSolstice,
            jd(2000, 6, 21, 1.0, 48.0, 0.0),
            63.8,
        ),
        (
            2000,
            Event::SeptemberEquinox,
            jd(2000, 9, 22, 17.0, 28.0, 0.0),
            63.8,
        ),
        (
            2000,
            Event::DecemberSolstice,
            jd(2000, 12, 21, 13.0, 37.0, 0.0),
            63.8,
        ),
        (
            2024,
            Event::MarchEquinox,
            jd(2024, 3, 20, 3.0, 6.0, 0.0),
            69.2,
        ),
        (
            2024,
            Event::JuneSolstice,
            jd(2024, 6, 20, 20.0, 51.0, 0.0),
            69.2,
        ),
        (
            2024,
            Event::SeptemberEquinox,
            jd(2024, 9, 22, 12.0, 44.0, 0.0),
            69.2,
        ),
        (
            2024,
            Event::DecemberSolstice,
            jd(2024, 12, 21, 9.0, 21.0, 0.0),
            69.2,
        ),
    ];
    for (year, e, ut, delta_t) in cases {
        let jde = season::instant(year, e);
        // 30 s of rounding plus a little margin
        assert_close(jde, ut + delta_t * SECOND, 35.0 * SECOND);
    }
}

// Meeus, table 27.D: for 1951–2050 the approximate method (mean instant
// plus periodic terms) is never more than 51 s from the correct instant.
#[test]
fn approx_instant_accuracy_1951_2050() {
    let mut largest: f64 = 0.0;
    for year in 1951..=2050 {
        for e in Event::ALL {
            let err = (season::approx_instant(year, e) - season::instant(year, e)).abs();
            largest = largest.max(err);
        }
    }
    let largest = largest / SECOND;
    assert!(largest < 53.0, "largest error {largest} s");
    assert!(largest > 45.0, "largest error {largest} s");
}

// The mean instants of tables 27.A and 27.B stay within 20 minutes of the
// accurate instants over the years −1000 to +3000.
#[test]
fn mean_instant_accuracy() {
    for year in (-1000..=3000).step_by(125) {
        for e in Event::ALL {
            let mean = season::mean_instant(year, e);
            let accurate = season::instant(year, e);
            assert!(
                (mean - accurate).abs() < 20.0 / 1440.0,
                "{year} {e:?}: {} min",
                (mean - accurate) * 1440.0
            );
        }
    }
}

// The four events of a year come in order, in the right months (in the
// Julian calendar before 1582 October 15), and the same event recurs after
// a tropical year.
#[test]
fn ordering_and_months() {
    for year in [
        -1000, -500, -1, 0, 1, 1000, 1582, 1583, 1900, 2000, 2100, 3000,
    ] {
        let t: Vec<f64> = Event::ALL
            .iter()
            .map(|&e| season::instant(year, e))
            .collect();
        let next = season::instant(year + 1, Event::MarchEquinox);
        for k in 0..4 {
            let end = if k == 3 { next } else { t[k + 1] };
            let d = end - t[k];
            assert!((88.0..95.0).contains(&d), "{year} event {k}: {d} days");

            let (y, m, _) = time::date_frm_julian_day(t[k]).unwrap();
            assert_eq!(y as i32, year);
            // in −1000 the Julian calendar is about 10 days ahead of the
            // seasons, so the March equinox falls in early April
            let expected_month = 3 + 3 * k as u8;
            assert!(
                m == expected_month || (year < 0 && m == expected_month + 1),
                "{year} event {k} in month {m}"
            );
        }
        assert_close(next - t[0], 365.2422, 0.03);
    }
}

// Meeus, chapter 27 (table 27.F): durations of the seasons in 2000 match the
// differences of the published instants of 2000 and 2001.
#[test]
fn durations_2000() {
    assert_close(season::duration(2000, Season::Spring), 92.76, 0.006);
    assert_close(season::duration(2000, Season::Summer), 93.65, 0.006);
    assert_close(season::duration(2000, Season::Autumn), 89.84, 0.006);
    assert_close(season::duration(2000, Season::Winter), 88.99, 0.006);
}

// The four seasons add up to a tropical year, and the winter of a year
// ends at the March equinox of the next one.
#[test]
fn durations_add_up() {
    for year in [-1000, 0, 1246, 2000, 3500] {
        let total: f64 = Season::ALL.iter().map(|&s| season::duration(year, s)).sum();
        let expected = season::instant(year + 1, Event::MarchEquinox)
            - season::instant(year, Event::MarchEquinox);
        assert_close(total, expected, 1e-6);
        assert_close(total, 365.2422, 0.03);
    }
}

// Meeus, chapter 27: in 1246 the perihelion coincided with the December
// solstice, so spring and summer had the same length, and so had autumn
// and winter. Since then the northern winter has been the shortest season.
#[test]
fn seasons_symmetry_1246_and_shortest_winter() {
    let d = |s| season::duration(1246, s);
    assert_close(d(Season::Spring), d(Season::Summer), 0.02);
    assert_close(d(Season::Autumn), d(Season::Winter), 0.02);

    for year in [1400, 1700, 2000, 2500, 3000, 3500] {
        let w = season::duration(year, Season::Winter);
        for s in [Season::Spring, Season::Summer, Season::Autumn] {
            assert!(w < season::duration(year, s), "{year}: winter not shortest");
        }
    }
}

#[test]
fn season_start_end() {
    for (k, s) in Season::ALL.iter().enumerate() {
        assert_eq!(s.start(), Event::ALL[k]);
        assert_eq!(s.end(), Event::ALL[(k + 1) % 4]);
    }
    assert_close(Event::DecemberSolstice.solar_longitude(), 1.5 * PI, 1e-15);
}

// Out-of-range years do not panic.
#[test]
fn extreme_years_do_not_panic() {
    for year in [i32::MIN, -100000, 100000, i32::MAX] {
        for e in Event::ALL {
            let _ = season::mean_instant(year, e);
            let _ = season::approx_instant(year, e);
        }
        let _ = season::duration(year, Season::Winter);
    }
}

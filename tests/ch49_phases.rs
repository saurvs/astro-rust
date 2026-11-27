//! Meeus, chapter 49: phases of the Moon.

mod common;
mod lunar_events_support;

use astro::lunar::phase::{self, Phase, PhaseEvent};
use astro::{lunar, time};
use common::assert_close;
use lunar_events_support as support;
use time::Month::*;

fn date(year: i16, month: time::Month, day: f64) -> time::Date {
    time::Date {
        year,
        month,
        decimal_day: day,
        cal_type: time::CalType::Gregorian,
    }
}

// Meeus, example 49.a: the New Moon of 1977 February.
#[test]
fn example_49a_new_moon_1977_feb() {
    // "Mid-February 1977" = 1977.13 gives k ~ -282.87, so k = -283.
    assert_close(phase::k_from_year(1977.13), -282.87, 0.005);
    assert_close(phase::mean_phase(-283, Phase::New), 2443192.94102, 0.000005);
    assert_close(phase::true_phase(-283, Phase::New), 2443192.65118, 0.000005);

    // The same New Moon found from a date in mid-February.
    let jde = lunar::time_of_phase(&date(1977, Feb, 15.0), &Phase::New);
    assert_close(jde, 2443192.65118, 0.000005);
    let event = phase::nearest_phase(support::jd(1977, Feb, 15.0), Phase::New);
    assert_eq!(event.lunation, -283);
    assert_eq!(event.k(), -283.0);
    // 1977 Feb 18.15118 TD = 3h37m42s TD
    assert_close(event.jde - support::jd(1977, Feb, 18.0), 0.15118, 0.000005);
}

// Meeus, example 49.b: the first Last Quarter of A.D. 2044.
#[test]
fn example_49b_last_quarter_2044() {
    assert_close(phase::k_from_year(2044.0), 544.21, 0.005);
    assert_close(phase::mean_phase(544, Phase::Last), 2467636.88597, 0.000005);
    assert_close(phase::true_phase(544, Phase::Last), 2467636.49186, 0.000005);

    let jde = lunar::time_of_phase(&date(2044, Jan, 15.0), &Phase::Last);
    assert_close(jde, 2467636.49186, 0.000005);
    // 2044 Jan 21, 23h48m17s TD
    let day = jde - support::jd(2044, Jan, 21.0);
    assert_close(day * 86400.0, ((23.0 * 60.0) + 48.0) * 60.0 + 17.0, 1.0);
}

// The quarter-phase correction W lost its last three terms to a stray
// semicolon until PR #12 (Thomas Daniels). These quarters are chosen where
// those terms are largest (5.9e-5 and 5.0e-5 day); the values were computed
// with this implementation, which reproduces examples 49.a and 49.b, and
// agree with published times (2020 Apr 1 10:21 UT, 2012 Nov 7 00:36 UT).
#[test]
fn quarter_phases_keep_all_terms_of_w() {
    let first = phase::true_phase(250, Phase::First);
    assert_close(first, 2458940.9321417, 0.000002);
    let last = phase::true_phase(158, Phase::Last);
    assert_close(last, 2456238.5256256, 0.000002);
}

#[test]
fn true_phase_is_close_to_mean_phase() {
    for n in (-30000..30000).step_by(97) {
        for p in Phase::ALL {
            let d = phase::true_phase(n, p) - phase::mean_phase(n, p);
            assert!(d.abs() < 1.0, "lunation {n} {p:?}: {d}");
        }
    }
}

#[test]
fn phase_constants() {
    let fr: Vec<f64> = Phase::ALL.iter().map(|p| p.fraction()).collect();
    assert_eq!(fr, [0.0, 0.25, 0.5, 0.75]);
    use std::f64::consts::PI;
    assert_close(Phase::First.longitude_excess(), PI / 2.0, 1e-15);
    assert_close(Phase::Last.longitude_excess(), 1.5 * PI, 1e-15);
    // k = 0 is the New Moon of 2000 January 6.
    let jde = phase::true_phase(0, Phase::New);
    assert!((jde - support::jd(2000, Jan, 6.0)).abs() < 1.0);
    // Meeus' examples of admissible k values.
    for (lunation, p, k) in [
        (479, Phase::New, 479.0),
        (-2793, Phase::New, -2793.0),
        (479, Phase::First, 479.25),
        (-2793, Phase::First, -2792.75),
        (479, Phase::Full, 479.5),
        (-2793, Phase::Full, -2792.5),
        (479, Phase::Last, 479.75),
        (-2793, Phase::Last, -2792.25),
    ] {
        let e = PhaseEvent {
            phase: p,
            lunation,
            jde: 0.0,
        };
        assert_eq!(e.k(), k);
    }
}

#[test]
fn k_from_jd_matches_mean_phases() {
    for n in [-3000, -283, 0, 544, 3000] {
        for p in Phase::ALL {
            let k = phase::k_from_jd(phase::mean_phase(n, p));
            // Only the small T^2 terms of eq. 49.1 separate them.
            assert_close(k, n as f64 + p.fraction(), 0.01);
        }
    }
}

// Versions before 3.0 truncated k towards zero. For dates before 2000 that
// gave the phase of the *next* lunation: on 1990 Jan 1 the New Moon of 1990
// Jan 26 instead of that of 1989 Dec 28, four days earlier.
#[test]
fn time_of_phase_returns_the_nearest_phase() {
    let jd = support::jd(1990, Jan, 1.0);
    let new = lunar::time_of_phase(&date(1990, Jan, 1.0), &Phase::New);
    assert!((new - jd).abs() < 5.0, "{new}");
    assert_close(new, phase::true_phase(-124, Phase::New), 1e-9);

    // ... and for dates after 2000 the one of the lunation that began
    // before the date: on 2024 Mar 28 the New Moon of Mar 10 instead of
    // that of Apr 8, eleven days later.
    let jd = support::jd(2024, Mar, 28.0);
    let new = lunar::time_of_phase(&date(2024, Mar, 28.0), &Phase::New);
    assert!(new > jd && new - jd < 12.0, "{new}");
}

#[test]
fn nearest_phase_is_never_farther_than_its_neighbours() {
    let start = support::jd(1999, Nov, 1.0);
    let mut jd = start;
    while jd < start + 800.0 {
        for p in Phase::ALL {
            let e = phase::nearest_phase(jd, p);
            assert_eq!(e.phase, p);
            let d = (e.jde - jd).abs();
            assert!(d <= 15.5, "{jd} {p:?} {d}");
            for n in [e.lunation - 1, e.lunation + 1] {
                assert!(d <= (phase::true_phase(n, p) - jd).abs());
            }
        }
        jd += 0.37;
    }
}

#[test]
fn nearest_phase_works_far_from_2000() {
    for jd in [0.0, 625000.5, 1721423.5, 2299160.5, 3000000.5] {
        for p in Phase::ALL {
            let e = phase::nearest_phase(jd, p);
            assert!((e.jde - jd).abs() <= 15.5, "{jd} {p:?} {e:?}");
        }
    }
}

#[test]
fn non_finite_input_gives_nan() {
    for jd in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(phase::nearest_phase(jd, Phase::Full).jde.is_nan());
    }
}

#[test]
fn phases_of_2024() {
    let list: Vec<_> =
        phase::phases_between(support::jd(2024, Jan, 1.0), support::jd(2025, Jan, 1.0)).collect();
    assert_eq!(list.len(), 50);
    let count = |p: Phase| list.iter().filter(|e| e.phase == p).count();
    assert_eq!(count(Phase::New), 13);
    assert_eq!(count(Phase::First), 12);
    assert_eq!(count(Phase::Full), 12);
    assert_eq!(count(Phase::Last), 13);
    // First event: Last Quarter of 2024 Jan 4; last: New Moon of Dec 30.
    assert_eq!(list[0].phase, Phase::Last);
    assert!((list[0].jde - support::jd(2024, Jan, 4.0)).abs() < 1.0);
    let last = list.last().unwrap();
    assert_eq!(last.phase, Phase::New);
    assert!((last.jde - support::jd(2024, Dec, 30.9)).abs() < 1.0);
}

#[test]
fn phases_between_is_ordered_and_complete() {
    for (start, end) in [
        (support::jd(1899, Dec, 20.0), support::jd(1901, Jan, 10.0)),
        (support::jd(1999, Dec, 31.0), support::jd(2000, Jan, 2.0)),
        (-1000.0, 0.0),
        (1500000.0, 1500200.0),
    ] {
        let list: Vec<_> = phase::phases_between(start, end).collect();
        for e in &list {
            assert!(e.jde >= start && e.jde < end);
            assert_eq!(e.jde, phase::true_phase(e.lunation, e.phase));
        }
        for w in list.windows(2) {
            let gap = w[1].jde - w[0].jde;
            assert!((6.0..9.0).contains(&gap), "{gap}");
            let next = (w[0].lunation * 4
                + Phase::ALL.iter().position(|p| *p == w[0].phase).unwrap() as i64)
                + 1;
            assert_eq!(w[1].lunation, next.div_euclid(4));
            assert_eq!(w[1].phase, Phase::ALL[next.rem_euclid(4) as usize]);
        }
        // Nothing is missed just outside the range.
        let before = phase::phases_between(start - 10.0, start).last();
        if let (Some(b), Some(f)) = (before, list.first()) {
            assert!(f.jde - b.jde < 9.0);
        }
    }
}

#[test]
fn phases_between_bounds() {
    let e = phase::nearest_phase(support::jd(2010, May, 1.0), Phase::Full);
    // Start is inclusive, end exclusive.
    let v: Vec<_> = phase::phases_between(e.jde, e.jde + 1.0).collect();
    assert_eq!(v, [e]);
    assert_eq!(phase::phases_between(e.jde - 1.0, e.jde).count(), 0);
    // Empty and degenerate ranges.
    assert_eq!(phase::phases_between(10.0, 10.0).count(), 0);
    assert_eq!(phase::phases_between(10.0, -10.0).count(), 0);
    assert_eq!(phase::phases_between(f64::NAN, 10.0).count(), 0);
    assert_eq!(phase::phases_between(0.0, f64::NAN).count(), 0);
    assert_eq!(phase::phases_between(f64::NEG_INFINITY, 0.0).count(), 0);
    // Open-ended.
    assert_eq!(
        phase::phases_between(0.0, f64::INFINITY).take(100).count(),
        100
    );
    // Fused.
    let mut it = phase::phases_between(0.0, 30.0);
    while it.next().is_some() {}
    assert!(it.next().is_none());
}

// Table 49.A: between 1900 and 2100 lunations last from 29d06h35m to
// 29d19h55m.
#[test]
fn lunation_lengths_1900_to_2100() {
    let news: Vec<_> =
        phase::phases_between(support::jd(1900, Jan, 1.0), support::jd(2100, Jan, 1.0))
            .filter(|e| e.phase == Phase::New)
            .collect();
    let (mut shortest, mut longest) = (f64::MAX, 0.0_f64);
    for w in news.windows(2) {
        let d = w[1].jde - w[0].jde;
        shortest = shortest.min(d);
        longest = longest.max(d);
    }
    let dhm = |d: f64, h: f64, m: f64| d + h / 24.0 + m / 1440.0;
    assert_close(shortest, dhm(29.0, 6.0, 35.0), 2.0 / 1440.0);
    assert_close(longest, dhm(29.0, 19.0, 55.0), 2.0 / 1440.0);
}

// Phases far in the past, across year 0 and before the Gregorian reform.
#[test]
fn ancient_phases_are_regular() {
    let list: Vec<_> =
        phase::phases_between(support::jd(-1000, Jan, 1.0), support::jd(-998, Jan, 1.0)).collect();
    assert!(list.len() > 95);
    for w in list.windows(2) {
        let gap = w[1].jde - w[0].jde;
        assert!((6.0..9.0).contains(&gap), "{gap}");
    }
}

// Independent check: at the computed instant, the Moon's longitude (ch. 47)
// minus the Sun's apparent longitude (VSOP87) is 0, 90, 180 or 270 degrees.
// Meeus reports errors up to 17 s for 1980-2020; the ch. 47 series adds a
// few arcseconds, so allow 40 s.
#[test]
fn phases_agree_with_lunar_and_solar_theories() {
    for n in (-250..250).step_by(7) {
        for p in Phase::ALL {
            let jde = phase::true_phase(n, p);
            let excess = |t: f64| {
                support::wrap_pi(
                    support::moon(t).0
                        - support::sun_long_without_nutation(t)
                        - p.longitude_excess(),
                )
            };
            let rate = (excess(jde + 0.01) - excess(jde - 0.01)) / 0.02;
            let dt = excess(jde) / rate * 86400.0;
            assert!(dt.abs() < 40.0, "lunation {n} {p:?}: {dt} s");
        }
    }
}

// A phase a few hours before a year boundary: the Full Moon of 2009 Dec 31
// (19h13m UT, i.e. about 19h14m TD) is the one nearest to 2010 Jan 1.0.
#[test]
fn full_moon_at_the_end_of_2009() {
    let jde = lunar::time_of_phase(&date(2010, Jan, 1.0), &Phase::Full);
    let hours = (jde - support::jd(2009, Dec, 31.0)) * 24.0;
    assert_close(hours, 19.0 + 14.0 / 60.0, 1.0 / 60.0);
    let in_2009: Vec<_> =
        phase::phases_between(support::jd(2009, Jan, 1.0), support::jd(2010, Jan, 1.0))
            .filter(|e| e.phase == Phase::Full)
            .collect();
    assert_eq!(in_2009.len(), 13); // a "blue moon" year
    assert_eq!(in_2009.last().unwrap().jde, jde);
}

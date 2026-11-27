//! Meeus, chapter 52: maximum declinations of the Moon.

mod common;
mod lunar_events_support;

use astro::lunar::declination::{self, Extreme, MaxDeclination};
use astro::time;
use common::{assert_close, dms};
use lunar_events_support as support;
use time::Month::*;

fn deg(rad: f64) -> f64 {
    rad.to_degrees()
}

// Meeus, example 52.a: greatest northern declination in December 1988.
#[test]
fn example_52a_north_1988_dec() {
    assert_close(declination::k_from_year(1988.95), -148.12, 0.005);
    let m = declination::max_declination(-148, Extreme::North);
    assert_eq!((m.extreme, m.number), (Extreme::North, -148));
    // The printed terms give .33465, which the book rounds to .3347.
    assert_close(m.jde, 2447518.3347, 0.0001);
    assert_close(deg(m.declination), 28.1562, 0.0001);
    assert_close(deg(m.declination), dms(false, 28, 9, 22.0), 1.0 / 3600.0);
    // 1988 Dec 22 at 20h02m TD
    assert_close(
        (m.jde - support::jd(1988, Dec, 22.0)) * 1440.0,
        20.0 * 60.0 + 2.0,
        1.0,
    );

    let near = declination::nearest_max_declination(support::jd(1988, Dec, 15.0), Extreme::North);
    assert_eq!(near, m);

    // The accurate values quoted by Meeus are 20h01m TD and +28°09'13"; the
    // chapter 47 series agrees with them.
    let t = support::golden_min(m.jde - 1.0, m.jde + 1.0, |t| -support::moon_apparent_dec(t));
    assert_close(
        (t - support::jd(1988, Dec, 22.0)) * 1440.0,
        20.0 * 60.0 + 1.0,
        1.5,
    );
    assert_close(
        deg(support::moon_apparent_dec(t)),
        dms(false, 28, 9, 13.0),
        3.0 / 3600.0,
    );
}

// Meeus, example 52.b: greatest southern declination for k = 659.
#[test]
fn example_52b_south_2049_apr() {
    let m = declination::max_declination(659, Extreme::South);
    assert_close(m.jde, 2469553.0834, 0.0001);
    // Meeus prints the absolute value, 22.1384 (-22°08'); here it is signed.
    assert_close(deg(m.declination), -22.1384, 0.0001);
    // 2049 Apr 21 at 14h TD
    assert_close((m.jde - support::jd(2049, Apr, 21.0)) * 24.0, 14.0, 0.05);
}

// Meeus, example 52.c: greatest northern declination of mid-March -4.
#[test]
fn example_52c_north_year_minus_4() {
    // "year" = -4 + 0.20 = -3.80 gives k ~ -26788.40.
    assert_close(declination::k_from_year(-3.80), -26788.40, 0.005);
    let m = declination::max_declination(-26788, Extreme::North);
    // The printed terms give 1719672.14125; the book prints .1414.
    assert_close(m.jde, 1719672.1414, 0.0002);
    assert_close(deg(m.declination), 28.9739, 0.0001);
    // March 16 at 15h TD of the year -4 (Julian calendar)
    let march_16 = time::julian_day(&time::Date {
        year: -4,
        month: Mar,
        decimal_day: 16.0,
        cal_type: time::CalType::Julian,
    });
    assert_close((m.jde - march_16) * 24.0, 15.0, 0.5);
}

#[test]
fn numbering_and_signs() {
    // k = 0: mean southern extreme on 2000 Jan 5.4, northern on Jan 19.1;
    // the true ones are on Jan 6.4 and Jan 19.9.
    let s = declination::max_declination(0, Extreme::South);
    let n = declination::max_declination(0, Extreme::North);
    assert!((s.jde - support::jd(2000, Jan, 6.4)).abs() < 0.1, "{s:?}");
    assert!((n.jde - support::jd(2000, Jan, 19.9)).abs() < 0.1, "{n:?}");
    assert_close(
        declination::k_from_jd(support::jd(2000, Jan, 5.4), Extreme::South),
        0.0,
        0.01,
    );
    assert_close(
        declination::k_from_jd(support::jd(2000, Jan, 19.1), Extreme::North),
        0.0,
        0.01,
    );
    assert!(s.declination < 0.0 && n.declination > 0.0);
    for kind in [Extreme::North, Extreme::South] {
        assert_close(
            declination::k_from_jd(declination::max_declination(-148, kind).jde, kind),
            -148.0,
            0.05,
        );
    }
}

// 1987 was near a major standstill (extremes about ±28.5°), 1997 and 2015
// near minor standstills (about ±18.2°), 2006 a major one again.
#[test]
fn standstills() {
    let extremes_in = |year: i16| -> Vec<MaxDeclination> {
        declination::max_declinations_between(
            support::jd(year, Jan, 1.0),
            support::jd(year + 1, Jan, 1.0),
        )
        .collect()
    };
    let largest = |v: &[MaxDeclination]| {
        v.iter()
            .map(|m| deg(m.declination).abs())
            .fold(0.0, f64::max)
    };
    for year in [1987, 2006] {
        let v = extremes_in(year);
        assert!(largest(&v) > 28.4, "{year}: {}", largest(&v));
    }
    for year in [1997, 2015] {
        let v = extremes_in(year);
        assert!(largest(&v) < 19.0, "{year}: {}", largest(&v));
    }
}

#[test]
fn values_and_spacing_1977_to_2022() {
    let list: Vec<_> = declination::max_declinations_between(
        support::jd(1977, Aug, 1.0),
        support::jd(2022, June, 30.0),
    )
    .collect();
    assert!((1195..=1201).contains(&list.len()), "{}", list.len());
    for m in &list {
        let d = deg(m.declination);
        match m.extreme {
            Extreme::North => assert!((17.9..28.8).contains(&d), "{m:?}"),
            Extreme::South => assert!((-28.8..-17.9).contains(&d), "{m:?}"),
        }
    }
    for w in list.windows(2) {
        assert_ne!(w[0].extreme, w[1].extreme);
        if w[0].extreme == Extreme::South {
            assert_eq!(w[0].number, w[1].number);
        } else {
            assert_eq!(w[0].number + 1, w[1].number);
        }
        let gap = w[1].jde - w[0].jde;
        assert!((12.0..15.5).contains(&gap), "{gap}");
    }
}

// Independent check against the apparent declination from the chapter 47
// series, 1977-2022. Meeus quotes largest errors of 10 min and 26".
#[test]
fn extremes_agree_with_chapter_47_declination() {
    for n in (-300..300).step_by(3) {
        for kind in [Extreme::North, Extreme::South] {
            let m = declination::max_declination(n, kind);
            let sign = if kind == Extreme::North { -1.0 } else { 1.0 };
            let t = support::golden_min(m.jde - 1.0, m.jde + 1.0, |t| {
                sign * support::moon_apparent_dec(t)
            });
            let dt = (t - m.jde) * 1440.0;
            let dd = deg(m.declination - support::moon_apparent_dec(t)) * 3600.0;
            assert!(dt.abs() < 12.0, "{n} {kind:?}: {dt} min");
            assert!(dd.abs() < 30.0, "{n} {kind:?}: {dd}\"");
        }
    }
}

#[test]
fn nearest_and_bounds() {
    let start = support::jd(2030, Jan, 1.0);
    let mut jd = start;
    while jd < start + 200.0 {
        for kind in [Extreme::North, Extreme::South] {
            let m = declination::nearest_max_declination(jd, kind);
            let d = (m.jde - jd).abs();
            assert!(d < 14.5);
            for k in [m.number - 1, m.number + 1] {
                assert!(d <= (declination::max_declination(k, kind).jde - jd).abs());
            }
        }
        jd += 0.43;
    }
    let m = declination::nearest_max_declination(f64::NAN, Extreme::South);
    assert!(m.jde.is_nan() && m.declination.is_nan());

    let e = declination::max_declination(400, Extreme::South);
    let v: Vec<_> = declination::max_declinations_between(e.jde, e.jde + 1.0).collect();
    assert_eq!(v, [e]);
    assert_eq!(
        declination::max_declinations_between(e.jde - 1.0, e.jde).count(),
        0
    );
    assert_eq!(
        declination::max_declinations_between(f64::INFINITY, 0.0).count(),
        0
    );
    assert_eq!(
        declination::max_declinations_between(0.0, f64::INFINITY)
            .take(3)
            .count(),
        3
    );
    // Before year 0
    let old: Vec<_> = declination::max_declinations_between(-100.0, 300.0).collect();
    assert!(old.len() >= 28);
}

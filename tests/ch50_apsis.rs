//! Meeus, chapter 50: perigee and apogee of the Moon.

mod common;
mod lunar_events_support;

use astro::lunar::apsis::{self, Apsis, ApsisEvent};
use astro::time::Month::*;
use common::assert_close;
use lunar_events_support as support;

const RADIUS: f64 = 6378.14;

fn arcsec(rad: f64) -> f64 {
    rad.to_degrees() * 3600.0
}

// Meeus, example 50.a: the apogee of 1988 October.
//
// The book prints JDE 2447442.3537 (sum of periodic terms -0.4654), but the
// terms of Table 50.A as printed sum to -0.4648, giving 2447442.3543 =
// 20h30m TD, the "exact" value quoted at the end of the example.
#[test]
fn example_50a_apogee_1988_oct() {
    // Early October = 1988.75 gives k ~ -148.73, hence k = -148.5.
    assert_close(apsis::k_from_year(1988.75), -148.73, 0.005);
    assert_close(
        apsis::mean_apsis(-149, Apsis::Apogee),
        2447442.8191,
        0.00005,
    );

    let a = apsis::true_apsis(-149, Apsis::Apogee);
    assert_eq!(a.apsis, Apsis::Apogee);
    assert_eq!(a.k(), -148.5);
    assert_close(a.jde, 2447442.3543, 0.0001);
    assert_close(arcsec(a.parallax), 3240.679, 0.001); // 0°54'00".679
    assert_close(a.distance, RADIUS / a.parallax.sin(), 1e-6);
    assert_close(a.distance, 405_976.7, 0.1);

    // 1988 Oct 7 at 20h30m TD (book: 20h29m, exact 20h30m)
    let minutes = (a.jde - support::jd(1988, Oct, 7.0)) * 1440.0;
    assert_close(minutes, 20.0 * 60.0 + 30.0, 1.0);

    let near = apsis::nearest_apsis(support::jd(1988, Oct, 1.0), Apsis::Apogee);
    assert_eq!(near, a);
}

#[test]
fn numbering() {
    // k = 0 is the perigee of 1999 December 22.
    let p = apsis::true_apsis(0, Apsis::Perigee);
    assert!((p.jde - support::jd(1999, Dec, 22.0)).abs() < 1.0, "{p:?}");
    // Meeus' examples of admissible k.
    for (number, kind, k) in [
        (318, Apsis::Perigee, 318.0),
        (-25, Apsis::Perigee, -25.0),
        (429, Apsis::Apogee, 429.5),
        (-1210, Apsis::Apogee, -1209.5),
    ] {
        let e = apsis::true_apsis(number, kind);
        assert_eq!(e.k(), k);
        assert_close(apsis::k_from_jd(apsis::mean_apsis(number, kind)), k, 0.01);
    }
    assert_eq!(Apsis::Perigee.fraction(), 0.0);
    assert_eq!(Apsis::Apogee.fraction(), 0.5);
}

// The examples of extreme intervals between perigees given in the text:
// 1997 Dec 9 16.9h to 1998 Jan 3 8.5h, and 1990 Dec 2 10.8h to Dec 30 23.8h.
#[test]
fn extreme_intervals_between_perigees() {
    for (y, m, d, h) in [
        (1997, Dec, 9.0, 16.9),
        (1998, Jan, 3.0, 8.5),
        (1990, Dec, 2.0, 10.8),
        (1990, Dec, 30.0, 23.8),
    ] {
        let expected = support::jd(y, m, d) + h / 24.0;
        let p = apsis::nearest_apsis(expected, Apsis::Perigee);
        assert_close(p.jde, expected, 31.0 / 1440.0);
    }
}

#[test]
fn distances_and_intervals_1977_to_2022() {
    let list: Vec<ApsisEvent> =
        apsis::apsides_between(support::jd(1977, June, 1.0), support::jd(2022, Aug, 31.0))
            .collect();
    assert!((1195..=1201).contains(&list.len()), "{}", list.len());
    let (perigees, apogees): (Vec<&ApsisEvent>, Vec<&ApsisEvent>) =
        list.iter().partition(|e| e.apsis == Apsis::Perigee);
    for p in &perigees {
        assert!((356_300.0..370_400.0).contains(&p.distance), "{p:?}");
        // The true perigee can be up to about 45 hours from the mean one.
        let d = p.jde - apsis::mean_apsis(p.number, Apsis::Perigee);
        assert!(d.abs() < 46.0 / 24.0, "{p:?} {d}");
    }
    for a in &apogees {
        assert!((404_000.0..406_750.0).contains(&a.distance), "{a:?}");
        let d = a.jde - apsis::mean_apsis(a.number, Apsis::Apogee);
        assert!(d.abs() < 14.0 / 24.0, "{a:?} {d}");
    }
    // Successive perigees are 24d16h to 28d13h apart, apogees 26.98 to 27.90 d.
    for w in perigees.windows(2) {
        let d = w[1].jde - w[0].jde;
        assert!((24.6..28.6).contains(&d), "{d}");
    }
    for w in apogees.windows(2) {
        let d = w[1].jde - w[0].jde;
        assert!((26.95..27.93).contains(&d), "{d}");
    }
}

// The extreme distances of 1500-2500 quoted in the text: 356 371 km at the
// perigee of 2257 January 1 and 406 720 km at the apogee of 2266 January 7.
// These were computed with the full ELP-2000/82, so allow for the errors of
// the chapter's series (12 and 6 km over 1977-2022).
#[test]
fn extreme_perigee_and_apogee_1500_to_2500() {
    let p = apsis::nearest_apsis(support::jd(2257, Jan, 1.0), Apsis::Perigee);
    assert!((p.jde - support::jd(2257, Jan, 1.5)).abs() < 1.0);
    assert_close(p.distance, 356_371.0, 15.0);
    let a = apsis::nearest_apsis(support::jd(2266, Jan, 7.0), Apsis::Apogee);
    assert!((a.jde - support::jd(2266, Jan, 7.5)).abs() < 1.0);
    assert_close(a.distance, 406_720.0, 10.0);

    let all: Vec<_> =
        apsis::apsides_between(support::jd(1500, Jan, 1.0), support::jd(2500, Jan, 1.0)).collect();
    let min = all.iter().map(|e| e.distance).fold(f64::MAX, f64::min);
    let max = all.iter().map(|e| e.distance).fold(0.0, f64::max);
    assert_close(min, 356_371.0, 20.0);
    assert_close(max, 406_720.0, 10.0);
}

// Independent check against the distance from the chapter 47 series: the
// numerical minimum (maximum) of the distance, 1977-2022. Meeus quotes
// largest errors of 31 min and 12 km (perigee), 3 min and 6 km (apogee)
// against ELP-2000/82; the chapter 47 series adds a little.
#[test]
fn apsides_agree_with_chapter_47_distances() {
    for n in (-280..330).step_by(3) {
        for kind in [Apsis::Perigee, Apsis::Apogee] {
            let e = apsis::true_apsis(n, kind);
            let sign = if kind == Apsis::Perigee { 1.0 } else { -1.0 };
            let t = support::golden_min(e.jde - 2.5, e.jde + 2.5, |t| sign * support::moon(t).2);
            let (dt_max, dd_max) = match kind {
                Apsis::Perigee => (35.0, 15.0),
                Apsis::Apogee => (6.0, 9.0),
            };
            let dt = (t - e.jde) * 1440.0;
            let dd = e.distance - support::moon(t).2;
            assert!(dt.abs() < dt_max, "{n} {kind:?}: {dt} min");
            assert!(dd.abs() < dd_max, "{n} {kind:?}: {dd} km");
        }
    }
}

#[test]
fn nearest_apsis_is_never_farther_than_its_neighbours() {
    let start = support::jd(1999, Dec, 1.0);
    let mut jd = start;
    while jd < start + 400.0 {
        for kind in [Apsis::Perigee, Apsis::Apogee] {
            let e = apsis::nearest_apsis(jd, kind);
            let d = (e.jde - jd).abs();
            for n in [e.number - 1, e.number + 1] {
                assert!(d <= (apsis::true_apsis(n, kind).jde - jd).abs());
            }
        }
        jd += 0.29;
    }
    for jd in [f64::NAN, f64::INFINITY] {
        let e = apsis::nearest_apsis(jd, Apsis::Perigee);
        assert!(e.jde.is_nan() && e.parallax.is_nan() && e.distance.is_nan());
    }
}

#[test]
fn apsides_between_alternate_and_respect_bounds() {
    let start = support::jd(-500, Mar, 1.0);
    let list: Vec<_> = apsis::apsides_between(start, start + 365.25).collect();
    assert!((26..=28).contains(&list.len()));
    for w in list.windows(2) {
        assert_ne!(w[0].apsis, w[1].apsis);
        let gap = w[1].jde - w[0].jde;
        assert!((11.0..17.0).contains(&gap), "{gap}");
        assert!(w[0].jde >= start && w[1].jde < start + 365.25);
    }
    let e = list[3];
    assert_eq!(
        apsis::apsides_between(e.jde, e.jde + 1.0).collect::<Vec<_>>(),
        [e]
    );
    assert_eq!(apsis::apsides_between(e.jde - 1.0, e.jde).count(), 0);
    assert_eq!(apsis::apsides_between(f64::NAN, 0.0).count(), 0);
    assert_eq!(apsis::apsides_between(1.0, 0.0).count(), 0);
    assert_eq!(
        apsis::apsides_between(0.0, f64::INFINITY).take(10).count(),
        10
    );
}

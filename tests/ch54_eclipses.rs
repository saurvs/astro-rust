//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 54: eclipses.

mod common;

use astro::eclipse::{
    self, LunarEclipse, LunarEclipseKind, Node, ShadowRule, SolarEclipse, SolarEclipseKind,
};
use common::assert_close;

/// Julian Day of a calendar date at 0h (Meeus, eq. 7.1).
fn jd(year: i32, month: u32, day: f64, gregorian: bool) -> f64 {
    let (y, m) = if month <= 2 {
        (year - 1, month + 12)
    } else {
        (year, month)
    };
    let (y, m) = (y as f64, m as f64);
    let b = if gregorian {
        let a = (y / 100.0).floor();
        2.0 - a + (a / 4.0).floor()
    } else {
        0.0
    };
    (365.25 * (y + 4716.0)).floor() + (30.6001 * (m + 1.0)).floor() + day + b - 1524.5
}

/// Julian Day of a Gregorian date and time.
fn jd_hm(year: i32, month: u32, day: u32, h: u32, m: u32, s: f64) -> f64 {
    jd(year, month, day as f64, true) + (h as f64 + m as f64 / 60.0 + s / 3600.0) / 24.0
}

const MINUTE: f64 = 1.0 / 1440.0;

fn solar_eclipses_between(from: f64, to: f64) -> Vec<SolarEclipse> {
    let mut v = Vec::new();
    let mut t = from;
    while let Some(e) = eclipse::next_solar_eclipse(t) {
        if e.jde >= to {
            break;
        }
        t = e.jde;
        v.push(e);
    }
    v
}

fn lunar_eclipses_between(from: f64, to: f64) -> Vec<LunarEclipse> {
    let mut v = Vec::new();
    let mut t = from;
    while let Some(e) = eclipse::next_lunar_eclipse(t) {
        if e.jde >= to {
            break;
        }
        t = e.jde;
        v.push(e);
    }
    v
}

// Meeus, example 54.a: the partial solar eclipse of 1993 May 21.
#[test]
fn solar_1993_may_21_ex54a() {
    let e = eclipse::solar_eclipse(-82).unwrap();
    assert_eq!(e.k(), -82.0);
    assert_close(e.gamma, 1.1348, 0.00005);
    assert_close(e.u, 0.0097, 0.00005);
    assert_eq!(e.kind, SolarEclipseKind::Partial);
    assert!(!e.central);
    assert_close(e.magnitude.unwrap(), 0.740, 0.0005);
    assert_eq!(e.node, Node::Descending);
    // 1993 May 21, 14h21.0m TD (book sum of rounded terms 2449129.0979).
    assert_close(e.jde, 2449129.0979, 0.0001);
    assert_close(
        e.jde,
        jd_hm(1993, 5, 21, 14, 21, 0.0),
        0.05 * MINUTE + 0.0001,
    );
    // Against the accurate values quoted by Meeus (14h20m14s TD,
    // γ = 1.1370, magnitude 0.735) the method is within its stated error.
    assert_close(e.jde, jd_hm(1993, 5, 21, 14, 20, 14.0), 1.1 * MINUTE);
    assert_close(e.gamma, 1.1370, 0.0024);
    assert_close(e.magnitude.unwrap(), 0.735, 0.006);
    assert_close(e.penumbra_radius(), 0.0097 + 0.5461, 0.00005);
}

// Meeus, example 54.b: the total solar eclipse of 2009 July 22.
#[test]
fn solar_2009_july_22_ex54b() {
    let e = eclipse::solar_eclipse(118).unwrap();
    assert_close(e.jde, 2455034.6088, 0.00005);
    assert_close(e.jde, jd_hm(2009, 7, 22, 2, 37, 0.0), 0.5 * MINUTE);
    assert_close(e.gamma, 0.0695, 0.00005);
    assert_close(e.u, -0.0157, 0.00005);
    assert!(e.central);
    assert_eq!(e.kind, SolarEclipseKind::Total);
    assert_eq!(e.node, Node::Descending);
    assert_eq!(e.magnitude, None);
}

// Meeus, example 54.c: the penumbral lunar eclipse of 1973 June 15.
#[test]
fn lunar_1973_june_ex54c() {
    let e = eclipse::lunar_eclipse(-329).unwrap();
    assert_eq!(e.k(), -328.5);
    assert_close(e.jde, 2441849.3687, 0.00005);
    assert_close(e.jde, jd_hm(1973, 6, 15, 20, 51, 0.0), 0.5 * MINUTE);
    assert_close(e.gamma, -1.3249, 0.00005);
    assert_close(e.u, 0.0197, 0.00005);
    assert_eq!(e.node, Node::Ascending);
    assert_close(e.umbral_mag, -0.609, 0.0005);
    assert_close(e.penumbral_mag, 0.462, 0.0005);
    assert_eq!(e.kind, LunarEclipseKind::Penumbral);
    assert_eq!(e.partial_semidur, None);
    assert_eq!(e.total_semidur, None);
    // Connaissance des Temps: greatest at 20h50.7m TD, penumbral
    // magnitude 0.469.
    assert_close(e.jde, jd_hm(1973, 6, 15, 20, 50, 42.0), 1.0 * MINUTE);
    assert_close(e.penumbral_mag, 0.469, 0.01);

    let c = e.contacts();
    assert_eq!(
        (c.umbra_begin, c.total_begin, c.total_end, c.umbra_end),
        (None, None, None, None)
    );
    assert_close(
        c.penumbra_end - c.penumbra_begin,
        2.0 * e.penumbral_semidur,
        1e-9,
    );
}

// Meeus, example 54.d: the first lunar eclipse after 1997 July 1, the total
// eclipse of 1997 September 16.
#[test]
fn lunar_1997_september_16_ex54d() {
    let e = eclipse::next_lunar_eclipse(jd(1997, 7, 1.0, true)).unwrap();
    assert_eq!(e.k(), -28.5);
    assert_eq!(eclipse::lunar_eclipse(-29), Some(e));
    // The two preceding Full Moons give no eclipse.
    assert_eq!(eclipse::lunar_eclipse(-31), None);
    assert_eq!(eclipse::lunar_eclipse(-30), None);

    assert_close(e.jde, 2450708.2835, 0.00005);
    assert_close(e.jde, jd_hm(1997, 9, 16, 18, 48, 12.0), 0.1 * MINUTE);
    assert_close(e.gamma, -0.3791, 0.00005);
    assert_close(e.u, -0.0131, 0.00005);
    assert_close(e.umbral_mag, 1.187, 0.0005);
    assert_eq!(e.kind, LunarEclipseKind::Total);
    assert_eq!(e.node, Node::Descending);

    // p = 1.0259, t = 0.4809, h = 1.5442 are the shadow radii plus or minus
    // the Moon's radius; the radii are ρ = 1.2848 + u, σ = 0.7403 - u.
    assert_close(e.penumbra_radius, 1.2848 - 0.0131, 0.00005);
    assert_close(e.umbra_radius, 0.7403 + 0.0131, 0.00005);

    // Semidurations, in minutes: partial 98, total 30, penumbral 153.
    let partial = e.partial_semidur.unwrap() * 1440.0;
    let total = e.total_semidur.unwrap() * 1440.0;
    let penumbral = e.penumbral_semidur * 1440.0;
    assert_close(
        partial,
        60.0 / 0.5856 * (1.0259_f64.powi(2) - 0.3791_f64.powi(2)).sqrt(),
        0.05,
    );
    assert_close(
        total,
        60.0 / 0.5856 * (0.4809_f64.powi(2) - 0.3791_f64.powi(2)).sqrt(),
        0.05,
    );
    assert_close(
        penumbral,
        60.0 / 0.5856 * (1.5442_f64.powi(2) - 0.3791_f64.powi(2)).sqrt(),
        0.05,
    );
    assert_close(partial, 98.0, 0.5);
    assert_close(total, 30.0, 0.5);
    assert_close(penumbral, 153.0, 0.5);

    // Contacts in UT with ΔT = 63 s, to the minute as in the book.
    let c = e.contacts();
    let ut = |t: f64| t - 63.0 / 86400.0;
    let d = jd(1997, 9, 16.0, true);
    let hm = |h: u32, m: u32| d + (h as f64 + m as f64 / 60.0) / 24.0;
    assert_close(ut(c.penumbra_begin), hm(16, 14), 1.0 * MINUTE);
    assert_close(ut(c.umbra_begin.unwrap()), hm(17, 9), 1.0 * MINUTE);
    assert_close(ut(c.total_begin.unwrap()), hm(18, 17), 1.0 * MINUTE);
    assert_close(ut(c.greatest), hm(18, 47), 1.0 * MINUTE);
    assert_close(ut(c.total_end.unwrap()), hm(19, 17), 1.0 * MINUTE);
    assert_close(ut(c.umbra_end.unwrap()), hm(20, 25), 1.0 * MINUTE);
    assert_close(ut(c.penumbra_end), hm(21, 20), 1.0 * MINUTE);
}

// Meeus, "Notes about the accuracy": limiting cases where the simple
// method gives the wrong answer, as the book says.
#[test]
fn limiting_cases_from_the_book() {
    // 1935 January 5 (k = -804): γ = -1.5395, u = -0.00464, so |γ| exceeds
    // 1.5433 + u and the method finds no eclipse (there was a tiny partial
    // one).
    assert_eq!(eclipse::solar_eclipse(-804), None);

    // 1957 April 30 (k = -528): γ = +0.9966 < 0.9972, so the method calls
    // it central; it was in fact a non-central annular eclipse.
    let e = eclipse::solar_eclipse(-528).unwrap();
    assert_close(e.gamma, 0.9966, 0.00005);
    assert!(e.central);
    assert_eq!(e.kind, SolarEclipseKind::Annular);

    // 1890 November 26 (k = -1349.5): umbral magnitude -0.007, so the
    // method finds only a penumbral eclipse (it was a small partial one).
    let e = eclipse::lunar_eclipse(-1350).unwrap();
    assert_close(e.umbral_mag, -0.007, 0.0005);
    assert_eq!(e.kind, LunarEclipseKind::Penumbral);
}

// Meeus, chapter 54 exercises.
#[test]
fn exercise_first_solar_eclipse_of_1979() {
    // 1979 February 26: total, visible from the northern hemisphere.
    let e = eclipse::next_solar_eclipse(jd(1979, 1, 1.0, true)).unwrap();
    assert_close(e.jde, jd(1979, 2, 26.7, true), 0.1);
    assert_eq!(e.kind, SolarEclipseKind::Total);
    assert!(e.central && e.gamma > 0.0);
}

#[test]
fn exercise_solar_eclipse_april_1977() {
    let e = eclipse::next_solar_eclipse(jd(1977, 4, 1.0, true)).unwrap();
    assert!(e.jde < jd(1977, 5, 1.0, true));
    assert_eq!(e.kind, SolarEclipseKind::Annular);
}

#[test]
fn exercise_no_solar_eclipse_in_july_1947() {
    let start = jd(1947, 7, 1.0, true);
    let end = jd(1947, 8, 1.0, true);
    assert!(eclipse::next_solar_eclipse(start).unwrap().jde > end);
    assert!(eclipse::prev_solar_eclipse(end).unwrap().jde < start);
}

#[test]
fn exercise_four_partial_solar_eclipses_in_2000() {
    let v = solar_eclipses_between(jd(2000, 1, 1.0, true), jd(2001, 1, 1.0, true));
    assert_eq!(v.len(), 4);
    assert!(
        v.iter()
            .all(|e| e.kind == SolarEclipseKind::Partial && !e.central)
    );
    // 2000 Feb 5, Jul 1, Jul 31, Dec 25
    let dates = [(2, 5), (7, 1), (7, 31), (12, 25)];
    for (e, (m, d)) in v.iter().zip(dates) {
        assert_close(e.jde, jd(2000, m, d as f64 + 0.5, true), 0.6);
    }
}

#[test]
fn exercise_no_lunar_eclipse_in_january_2008() {
    let v = lunar_eclipses_between(jd(2008, 1, 1.0, true), jd(2008, 2, 1.0, true));
    assert!(v.is_empty());
}

#[test]
fn exercise_three_total_lunar_eclipses_in_1982() {
    let v = lunar_eclipses_between(jd(1982, 1, 1.0, true), jd(1983, 1, 1.0, true));
    let total: Vec<_> = v
        .iter()
        .filter(|e| e.kind == LunarEclipseKind::Total)
        .collect();
    assert_eq!(total.len(), 3);
    // 1982 Jan 9, Jul 6, Dec 30
    for (e, (m, d)) in total.iter().zip([(1, 9), (7, 6), (12, 30)]) {
        assert_close(e.jde, jd(1982, m, d as f64 + 0.5, true), 0.6);
    }
}

#[test]
fn exercise_first_lunar_eclipse_of_1234() {
    // The partial lunar eclipse of 1234 March 17 (Julian calendar).
    let e = eclipse::next_lunar_eclipse(jd(1234, 1, 1.0, false)).unwrap();
    assert_close(e.jde, jd(1234, 3, 17.5, false), 0.6);
    assert_eq!(e.kind, LunarEclipseKind::Partial);
}

// Well-known eclipses, compared with the NASA canon values (Espenak and
// Meeus) within the accuracy of the method.
#[test]
fn total_solar_eclipse_1999_august_11() {
    let e = eclipse::next_solar_eclipse(jd(1999, 8, 1.0, true)).unwrap();
    assert_eq!(e.kind, SolarEclipseKind::Total);
    assert!(e.central);
    // Greatest eclipse 11h04m09s TD, γ = +0.5062
    assert_close(e.jde, jd_hm(1999, 8, 11, 11, 4, 9.0), 2.0 * MINUTE);
    assert_close(e.gamma, 0.5062, 0.003);
    assert_eq!(e.node, Node::Ascending);
}

#[test]
fn total_solar_eclipse_2017_august_21() {
    let e = eclipse::next_solar_eclipse(jd(2017, 8, 1.0, true)).unwrap();
    assert_eq!(e.kind, SolarEclipseKind::Total);
    assert!(e.central);
    // Greatest eclipse 18h26m40s TD, γ = +0.4367
    assert_close(e.jde, jd_hm(2017, 8, 21, 18, 26, 40.0), 2.0 * MINUTE);
    assert_close(e.gamma, 0.4367, 0.003);
    assert_eq!(e.node, Node::Ascending);
}

#[test]
fn annular_and_hybrid_solar_eclipses() {
    // Annular eclipses of 2012 May 20 and 2023 October 14.
    for &(y, m, d) in &[(2012, 5, 20), (2023, 10, 14)] {
        let e = eclipse::next_solar_eclipse(jd(y, m, d as f64 - 3.0, true)).unwrap();
        assert_close(e.jde, jd(y, m, d as f64 + 0.5, true), 0.6);
        assert_eq!(e.kind, SolarEclipseKind::Annular, "{y}-{m}-{d}");
        assert!(e.central);
        assert!(e.u > 0.0047);
    }
    // Hybrid (annular-total) eclipses of 2013 November 3 and 2023 April 20.
    for &(y, m, d) in &[(2013, 11, 3), (2023, 4, 20)] {
        let e = eclipse::next_solar_eclipse(jd(y, m, d as f64 - 3.0, true)).unwrap();
        assert_close(e.jde, jd(y, m, d as f64 + 0.5, true), 0.6);
        assert_eq!(e.kind, SolarEclipseKind::Hybrid, "{y}-{m}-{d}");
        assert!(e.central && e.magnitude.is_none());
    }
}

#[test]
fn non_central_solar_eclipses() {
    // 2014 April 29: annular, not central (Antarctica).
    let e = eclipse::next_solar_eclipse(jd(2014, 4, 26.0, true)).unwrap();
    assert_close(e.jde, jd(2014, 4, 29.25, true), 0.6);
    assert_eq!(e.kind, SolarEclipseKind::Annular);
    assert!(!e.central && e.gamma < 0.0);
    assert!(e.gamma.abs() < 0.9972 + e.u.abs());

    // 2043 April 9: total, not central (northern polar region).
    let e = eclipse::next_solar_eclipse(jd(2043, 4, 6.0, true)).unwrap();
    assert_close(e.jde, jd(2043, 4, 9.75, true), 0.6);
    assert_eq!(e.kind, SolarEclipseKind::Total);
    assert!(!e.central && e.gamma > 0.0);
}

#[test]
fn total_lunar_eclipse_2000_january_21() {
    let e = eclipse::next_lunar_eclipse(jd(2000, 1, 1.0, true)).unwrap();
    assert_eq!(e.kind, LunarEclipseKind::Total);
    // Greatest eclipse at 04h43.6m UT = 04h44.7m TD; umbral magnitude 1.33;
    // totality 77 min, umbral phase 204 min.
    assert_close(e.jde, jd_hm(2000, 1, 21, 4, 44, 40.0), 2.0 * MINUTE);
    assert_close(e.umbral_mag, 1.33, 0.01);
    assert_close(e.gamma.abs(), 0.2957, 0.003);
    assert_close(e.total_semidur.unwrap() * 1440.0, 77.0 / 2.0, 1.5);
    assert_close(e.partial_semidur.unwrap() * 1440.0, 204.4 / 2.0, 1.5);
}

#[test]
fn total_lunar_eclipse_2018_july_27() {
    // The longest total lunar eclipse of the 21st century.
    let e = eclipse::next_lunar_eclipse(jd(2018, 7, 1.0, true)).unwrap();
    assert_eq!(e.kind, LunarEclipseKind::Total);
    // Greatest eclipse 20h21m44s UT = 20h22m53s TD; umbral magnitude
    // 1.6087; γ = 0.1168; totality 103.0 min, umbral phase 234.5 min.
    assert_close(e.jde, jd_hm(2018, 7, 27, 20, 22, 53.0), 2.0 * MINUTE);
    assert_close(e.umbral_mag, 1.6087, 0.01);
    assert_close(e.gamma, 0.1168, 0.003);
    assert_close(e.total_semidur.unwrap() * 1440.0, 103.0 / 2.0, 1.5);
    assert_close(e.partial_semidur.unwrap() * 1440.0, 234.5 / 2.0, 1.5);
}

#[test]
fn traditional_shadow_rule() {
    // The 1/50 rule enlarges the shadow: magnitudes grow by about 0.005
    // (umbral) and 0.026 (penumbral) relative to the French rule.
    for n in [-329, -29, -1350, 0, 5, 222] {
        let (Some(f), Some(t)) = (
            eclipse::lunar_eclipse(n),
            eclipse::lunar_eclipse_with_rule(n, ShadowRule::Traditional),
        ) else {
            continue;
        };
        assert_close(t.umbral_mag - f.umbral_mag, 0.0029 / 0.5450, 1e-12);
        assert_close(t.penumbral_mag - f.penumbral_mag, 0.0137 / 0.5450, 1e-12);
        assert!(t.penumbral_semidur > f.penumbral_semidur);
        assert_close(t.penumbra_radius - f.penumbra_radius, 0.0137, 1e-12);
        assert_close(t.umbra_radius - f.umbra_radius, 0.0029, 1e-12);
        assert_eq!((t.jde, t.gamma, t.u), (f.jde, f.gamma, f.u));
    }
    assert_eq!(ShadowRule::default(), ShadowRule::French);
    assert_eq!(
        eclipse::lunar_eclipse_with_rule(-29, ShadowRule::French),
        eclipse::lunar_eclipse(-29)
    );
}

// Invariants over a long span of lunations (about 1600 to 2400).
#[test]
fn consistency_over_many_lunations() {
    let mut solar = 0;
    let mut lunar = 0;
    for n in -4900..5000 {
        let mean = 2451550.09766 + 29.530588861 * n as f64;
        if let Some(e) = eclipse::solar_eclipse(n) {
            solar += 1;
            assert!((e.jde - mean).abs() < 0.75);
            assert!(e.gamma.abs() <= 1.5433 + e.u);
            assert_eq!(e.central, e.gamma.abs() < 0.9972);
            match e.kind {
                SolarEclipseKind::Partial => {
                    let m = e.magnitude.unwrap();
                    assert!((0.0..1.0).contains(&m), "magnitude {m}");
                    assert!(!e.central);
                }
                SolarEclipseKind::Total => assert!(e.u < 0.0),
                SolarEclipseKind::Annular => assert!(e.u >= 0.0),
                SolarEclipseKind::Hybrid => assert!(e.central && (0.0..=0.0047).contains(&e.u)),
            }
        }
        if let Some(e) = eclipse::lunar_eclipse(n) {
            lunar += 1;
            assert!((e.jde - mean - 29.530588861 / 2.0).abs() < 0.75);
            assert!(e.penumbral_mag > 0.0);
            assert!(e.penumbral_mag > e.umbral_mag);
            let c = e.contacts();
            assert!(c.penumbra_begin < c.greatest && c.greatest < c.penumbra_end);
            match e.kind {
                LunarEclipseKind::Penumbral => {
                    assert!(e.umbral_mag <= 0.0);
                    assert!(e.partial_semidur.is_none() && e.total_semidur.is_none());
                }
                LunarEclipseKind::Partial => {
                    assert!(e.umbral_mag > 0.0 && e.umbral_mag < 1.0 + 1e-9);
                    let p = e.partial_semidur.unwrap();
                    assert!(p < e.penumbral_semidur && e.total_semidur.is_none());
                }
                LunarEclipseKind::Total => {
                    assert!(e.umbral_mag >= 1.0 - 1e-9);
                    let (p, t) = (e.partial_semidur.unwrap(), e.total_semidur.unwrap());
                    assert!(t < p && p < e.penumbral_semidur);
                    // Totality never lasts more than about 1h47m.
                    assert!(2.0 * t * 1440.0 < 110.0);
                }
            }
        }
    }
    // About 2.4 solar and 2.4 lunar (including penumbral) eclipses a year
    // over 800 years.
    assert!((1800..2100).contains(&solar), "{solar} solar eclipses");
    assert!((1800..2100).contains(&lunar), "{lunar} lunar eclipses");
}

#[test]
fn next_and_prev_are_inverse() {
    let mut t = jd(1990, 1, 1.0, true);
    for _ in 0..40 {
        let e = eclipse::next_solar_eclipse(t).unwrap();
        assert!(e.jde > t);
        // No eclipse is skipped: going back from just after it finds it.
        assert_eq!(eclipse::prev_solar_eclipse(e.jde + 1.0), Some(e));
        let prev = eclipse::prev_solar_eclipse(e.jde).unwrap();
        assert!(prev.jde < e.jde && e.jde - prev.jde < 190.0);
        t = e.jde;

        let l = eclipse::next_lunar_eclipse(t).unwrap();
        assert!(l.jde > t);
        assert_eq!(eclipse::prev_lunar_eclipse(l.jde + 1.0), Some(l));
    }
}

#[test]
fn ancient_and_far_future_dates() {
    // Negative years: the total solar eclipse of -584 May 28 (Julian), the
    // "Thales" eclipse, is found close to the expected date.
    let e = eclipse::next_solar_eclipse(jd(-584, 5, 20.0, false)).unwrap();
    assert_close(e.jde, jd(-584, 5, 28.5, false), 0.75);
    assert!(e.central);

    // Searches keep working far from J2000.
    assert!(eclipse::next_lunar_eclipse(jd(-3000, 1, 1.0, false)).is_some());
    assert!(eclipse::next_solar_eclipse(jd(4000, 1, 1.0, true)).is_some());
}

#[test]
fn bad_input_gives_none() {
    assert_eq!(eclipse::next_solar_eclipse(f64::NAN), None);
    assert_eq!(eclipse::prev_solar_eclipse(f64::INFINITY), None);
    assert_eq!(eclipse::next_lunar_eclipse(f64::NEG_INFINITY), None);
    assert_eq!(eclipse::prev_lunar_eclipse(f64::NAN), None);
    assert_eq!(eclipse::solar_eclipse(i64::MAX), None);
    assert_eq!(eclipse::lunar_eclipse(i64::MIN), None);
}

// Meeus, chapter 54: the seven non-central total or annular eclipses of
// 1950-2100 listed in the book. The method reproduces all of them except
// 1957 April 30, which the book itself names as a limiting case (γ =
// 0.9966 < 0.9972, so it is classified as central; see
// `limiting_cases_from_the_book`). No other eclipse of the period is
// classified as non-central total or annular.
#[test]
fn non_central_eclipses_1950_to_2100() {
    let expected = [
        ((1950, 3, 18), SolarEclipseKind::Annular),
        ((1957, 10, 23), SolarEclipseKind::Total),
        ((1967, 11, 2), SolarEclipseKind::Total),
        ((2014, 4, 29), SolarEclipseKind::Annular),
        ((2043, 4, 9), SolarEclipseKind::Total),
        ((2043, 10, 3), SolarEclipseKind::Annular),
    ];
    let found: Vec<_> = solar_eclipses_between(jd(1950, 1, 1.0, true), jd(2101, 1, 1.0, true))
        .into_iter()
        .filter(|e| !e.central && e.kind != SolarEclipseKind::Partial)
        .collect();
    assert_eq!(found.len(), expected.len());
    for (e, ((y, m, d), kind)) in found.iter().zip(expected) {
        assert_eq!(e.kind, kind, "{y}-{m}-{d}");
        // The date is in UT; the instants here are TD, well within a day.
        let day = jd(y, m, d as f64, true);
        assert!(e.jde > day && e.jde < day + 1.0, "{y}-{m}-{d}: {}", e.jde);
        // Non-central: the axis misses the Earth, but the edge of the
        // umbral cone reaches it.
        assert!(e.gamma.abs() > 0.9972 && e.gamma.abs() < 0.9972 + e.u.abs());
        assert_eq!(e.magnitude, None);
    }
}

// Number and types of the eclipses of the 21st century (2001-2100). The
// Five Millennium Canons of Espenak and Meeus list 224 solar eclipses (77
// partial, 72 annular, 68 total, 7 hybrid) and 228 lunar eclipses; the
// simple method of chapter 54 finds exactly the same solar eclipses with the
// same types. For the Moon, limiting cases (penumbral eclipses of tiny
// magnitude, umbral magnitudes within 0.01 of 0 or 1) may go either way.
#[test]
fn eclipses_of_the_21st_century() {
    let (from, to) = (jd(2001, 1, 1.0, true), jd(2101, 1, 1.0, true));

    let solar = solar_eclipses_between(from, to);
    let count = |kind| solar.iter().filter(|e| e.kind == kind).count();
    assert_eq!(solar.len(), 224);
    assert_eq!(count(SolarEclipseKind::Partial), 77);
    assert_eq!(count(SolarEclipseKind::Annular), 72);
    assert_eq!(count(SolarEclipseKind::Total), 68);
    assert_eq!(count(SolarEclipseKind::Hybrid), 7);

    let lunar = lunar_eclipses_between(from, to);
    assert!(
        (225..=231).contains(&lunar.len()),
        "{} lunar eclipses",
        lunar.len()
    );

    // The searches neither skip nor repeat an eclipse: going through the
    // lunations one by one gives the same lists.
    let k0 = ((from - 2451550.1) / 29.530588861).floor() as i64 - 1;
    let k1 = ((to - 2451550.1) / 29.530588861).ceil() as i64 + 1;
    let by_lunation: Vec<_> = (k0..k1)
        .filter_map(eclipse::solar_eclipse)
        .filter(|e| e.jde >= from && e.jde < to)
        .collect();
    assert_eq!(by_lunation, solar);
    let by_lunation: Vec<_> = (k0..k1)
        .filter_map(eclipse::lunar_eclipse)
        .filter(|e| e.jde >= from && e.jde < to)
        .collect();
    assert_eq!(by_lunation, lunar);

    // Backwards searches find the same eclipses.
    let mut back = Vec::new();
    let mut t = to;
    while let Some(e) = eclipse::prev_solar_eclipse(t) {
        if e.jde < from {
            break;
        }
        t = e.jde;
        back.push(e);
    }
    back.reverse();
    assert_eq!(back, solar);
}

// Over the five millennia -1999 to +3000 the canons list 11 898 solar and
// 12 064 lunar eclipses (penumbral included). The method stays within 1%
// of both, and consecutive eclipses are never more than six lunations
// apart.
#[test]
fn eclipses_over_five_millennia() {
    let (from, to) = (jd(-1999, 1, 1.0, false), jd(3001, 1, 1.0, true));
    let k0 = ((from - 2451550.1) / 29.530588861).floor() as i64 - 1;
    let k1 = ((to - 2451550.1) / 29.530588861).ceil() as i64 + 1;

    let solar: Vec<_> = (k0..k1)
        .filter_map(eclipse::solar_eclipse)
        .filter(|e| e.jde >= from && e.jde < to)
        .collect();
    let lunar: Vec<_> = (k0..k1)
        .filter_map(eclipse::lunar_eclipse)
        .filter(|e| e.jde >= from && e.jde < to)
        .collect();
    let n = solar.len() as f64;
    assert!((n / 11898.0 - 1.0).abs() < 0.01, "{n} solar eclipses");
    let n = lunar.len() as f64;
    assert!((n / 12064.0 - 1.0).abs() < 0.01, "{n} lunar eclipses");

    for w in solar.windows(2) {
        assert!(w[1].jde > w[0].jde && w[1].jde - w[0].jde < 6.1 * 29.53);
    }
    for w in lunar.windows(2) {
        assert!(w[1].jde > w[0].jde && w[1].jde - w[0].jde < 6.1 * 29.53);
    }
}

// The eclipse node is the lunar node type of chapter 51.
#[test]
fn eclipse_node_is_the_lunar_node() {
    let e = eclipse::solar_eclipse(118).unwrap(); // 2009 July 22
    let node: astro::lunar::node::Node = e.node;
    assert_eq!(node, astro::lunar::node::Node::Descending);
}

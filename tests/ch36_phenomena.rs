// Tests for planet::phenomena (Meeus, Astronomical Algorithms, 2nd ed.,
// chapter 36: oppositions, conjunctions, greatest elongations and stations).
//
// Besides the worked examples, the times given by the method are checked
// against an independent computation from the crate's own VSOP87 positions
// (planet::heliocent_coords): apparent geocentric directions of the planet
// (light-time corrected) and of the Sun, both corrected for annual
// aberration, and a numerical search for the event.

mod common;

use std::f64::consts::{PI, TAU};

use astro::planet::phenomena::{self, Elongation, Phenomenon, Station};
use astro::planet::{self, Planet};
use common::assert_close;

const HOUR: f64 = 1.0 / 24.0;

/// Julian Day of a calendar date (Meeus, chapter 7), astronomical year
/// numbering; `gregorian = false` for the Julian calendar.
fn jd(year: i64, month: i64, day: f64, gregorian: bool) -> f64 {
    let (y, m) = if month <= 2 {
        (year - 1, month + 12)
    } else {
        (year, month)
    };
    let b = if gregorian {
        let a = (y as f64 / 100.0).floor();
        2.0 - a + (a / 4.0).floor()
    } else {
        0.0
    };
    (365.25 * (y as f64 + 4716.0)).floor() + (30.6001 * (m as f64 + 1.0)).floor() + day + b - 1524.5
}

/// Asserts that `jde` rounds to the given whole hour of the given date.
#[track_caller]
fn assert_hour(jde: f64, year: i64, month: i64, day: i64, hour: f64, gregorian: bool) {
    let expected = jd(year, month, day as f64, gregorian) + hour * HOUR;
    assert!(
        (jde - expected).abs() < 0.5 * HOUR,
        "JDE {jde} is not {year}-{month}-{day} {hour}h (off by {:.2} h)",
        (jde - expected) / HOUR
    );
}

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

const ALL_PHENOMENA: [Phenomenon; 8] = [
    Phenomenon::InferiorConjunction,
    Phenomenon::SuperiorConjunction,
    Phenomenon::Opposition,
    Phenomenon::Conjunction,
    Phenomenon::GreatestElongation(Elongation::East),
    Phenomenon::GreatestElongation(Elongation::West),
    Phenomenon::Station(Station::First),
    Phenomenon::Station(Station::Second),
];

fn planet_index(planet: &Planet) -> usize {
    match planet {
        Planet::Mercury => 0,
        Planet::Venus => 1,
        Planet::Earth => 2,
        Planet::Mars => 3,
        Planet::Jupiter => 4,
        Planet::Saturn => 5,
        Planet::Uranus => 6,
        Planet::Neptune => 7,
    }
}

fn expected_available(planet: &Planet, phenomenon: Phenomenon) -> bool {
    let i = planet_index(planet);
    let inner = i <= 1;
    match phenomenon {
        Phenomenon::InferiorConjunction | Phenomenon::SuperiorConjunction => inner,
        Phenomenon::GreatestElongation(_) => inner,
        Phenomenon::Opposition | Phenomenon::Conjunction => i >= 3,
        Phenomenon::Station(_) => inner || (3..=5).contains(&i),
    }
}

// ---------------------------------------------------------------------------
// Worked examples
// ---------------------------------------------------------------------------

// Meeus, example 36.a: Mercury's inferior conjunction nearest 1993 October 1
#[test]
fn example_36a_mercury_inferior_conjunction() {
    let ph = Phenomenon::InferiorConjunction;
    let p = Planet::Mercury;

    assert_close(phenomena::approx_k(&p, ph, 1993.75).unwrap(), -20.28, 0.005);
    let k = phenomena::nearest_k(&p, ph, 1993.75).unwrap();
    assert_eq!(k, -20);

    let mean = phenomena::mean_phenomenon(&p, ph, k).unwrap();
    assert_eq!(mean.k, -20);
    assert_close(mean.jde, 2449294.473, 0.001);
    assert_close(mean.earth_mean_anomaly.to_degrees(), 299.4092, 0.0001);
    assert_close(mean.julian_centuries, -0.06162, 0.00001);

    let jde = phenomena::time_of_phenomenon(&p, ph, k).unwrap();
    assert_close(jde - mean.jde, 3.171, 0.0005);
    assert_close(jde, 2449297.644, 0.001);
    assert_hour(jde, 1993, 11, 6, 3.0, true);
    assert_eq!(
        phenomena::nearest_time_of_phenomenon(&p, ph, 1993.75),
        Some(jde)
    );
}

// Meeus, example 36.b: the first conjunction of Saturn with the Sun after
// the beginning of 2125
#[test]
fn example_36b_saturn_conjunction() {
    let ph = Phenomenon::Conjunction;
    let p = Planet::Saturn;

    let k = phenomena::approx_k(&p, ph, 2125.0).unwrap();
    assert_close(k, 120.39, 0.005);
    let k = k.ceil() as i64;
    assert_eq!(k, 121);

    let mean = phenomena::mean_phenomenon(&p, ph, k).unwrap();
    assert_close(mean.jde, 2497430.244, 0.001);
    assert_close(mean.earth_mean_anomaly.to_degrees(), 222.0393, 0.0001);
    assert_close(mean.julian_centuries, 1.25627, 0.00001);

    let jde = phenomena::time_of_phenomenon(&p, ph, k).unwrap();
    assert_close(jde - mean.jde, 7.659, 0.0005);
    assert_close(jde, 2497437.903, 0.001);
    assert_hour(jde, 2125, 8, 26, 10.0, true);

    // Meeus gives 11h TD as the correct instant; so does VSOP87.
    let accurate = vsop::event(&p, Phenomenon::Conjunction, jde);
    assert_hour(accurate, 2125, 8, 26, 11.0, true);
}

// Meeus, example 36.c: Mercury's greatest western elongation of 1993
// November, counted from the inferior conjunction of example 36.a
#[test]
fn example_36c_mercury_greatest_western_elongation() {
    let p = Planet::Mercury;
    let e = phenomena::greatest_elongation(&p, Elongation::West, -20).unwrap();
    let mean =
        phenomena::mean_phenomenon(&p, Phenomenon::GreatestElongation(Elongation::West), -20)
            .unwrap();
    assert_close(mean.jde, 2449294.473, 0.001);
    assert_close(e.jde - mean.jde, 19.665, 0.001);
    assert_close(e.jde, 2449314.14, 0.005);
    assert_hour(e.jde, 1993, 11, 22, 15.0, true);
    assert_close(e.elongation.to_degrees(), 19.7506, 0.00005);

    let jde =
        phenomena::time_of_phenomenon(&p, Phenomenon::GreatestElongation(Elongation::West), -20);
    assert_eq!(jde, Some(e.jde));
}

// Meeus, example 36.d: the station of Mars following the opposition of
// 1997 March
#[test]
fn example_36d_mars_second_station() {
    let p = Planet::Mars;
    let ph = Phenomenon::Station(Station::Second);
    assert_eq!(
        phenomena::nearest_k(&p, Phenomenon::Opposition, 1997.2),
        Some(-2)
    );
    assert_eq!(phenomena::nearest_k(&p, ph, 1997.3), Some(-2));

    let mean = phenomena::mean_phenomenon(&p, ph, -2).unwrap();
    assert_close(mean.jde, 2450537.510, 0.001);
    assert_close(mean.earth_mean_anomaly.to_degrees(), 84.5468, 0.0001);
    assert_close(mean.julian_centuries, -0.02758, 0.00001);

    let jde = phenomena::time_of_phenomenon(&p, ph, -2).unwrap();
    assert_close(jde - mean.jde, 28.745, 0.0005);
    assert_close(jde, 2450566.255, 0.001);
    assert_hour(jde, 1997, 4, 27, 18.0, true);

    // "The correct time was 19h" (Meeus); VSOP87 agrees.
    let accurate = vsop::event(&p, ph, jde);
    assert_hour(accurate, 1997, 4, 27, 19.0, true);
}

// Meeus, chapter 36: Mercury's greatest elongation of 1990 August 11 reached
// an angular distance of 27°25' from the Sun at about 21h UT.
#[test]
fn mercury_greatest_elongation_1990() {
    let p = Planet::Mercury;
    let ph = Phenomenon::GreatestElongation(Elongation::East);
    let k = phenomena::nearest_k(&p, ph, 1990.61).unwrap();
    let e = phenomena::greatest_elongation(&p, Elongation::East, k).unwrap();
    assert_close(e.elongation.to_degrees(), 27.0 + 25.0 / 60.0, 0.5 / 60.0);
    let (t, value) = vsop::greatest_elongation(&p, e.jde);
    assert_close(value.to_degrees(), 27.0 + 25.0 / 60.0, 0.5 / 60.0);
    // ΔT was about a minute in 1990.
    assert_hour(t, 1990, 8, 11, 21.0, true);
    assert_close(e.jde, t, 0.5 * HOUR);
}

// Meeus, chapter 36, exercises (all instants in TD)
#[test]
fn exercises() {
    /// (year, month, day, hour TD)
    type Instant = (i64, i64, i64, f64);
    let cases: [(Planet, Phenomenon, f64, Instant, bool); 7] = [
        // first transit of Mercury ever observed (Gassendi, Paris)
        (
            Planet::Mercury,
            Phenomenon::InferiorConjunction,
            1631.85,
            (1631, 11, 7, 7.0),
            true,
        ),
        // last transit of Venus before 2004
        (
            Planet::Venus,
            Phenomenon::InferiorConjunction,
            1882.93,
            (1882, 12, 6, 17.0),
            true,
        ),
        // a perihelic opposition of Mars
        (
            Planet::Mars,
            Phenomenon::Opposition,
            2729.69,
            (2729, 9, 9, 3.0),
            true,
        ),
        // Jupiter and Saturn at opposition within a day, in -6 (7 B.C.):
        // a triple conjunction of the two planets that year (Julian calendar)
        (
            Planet::Jupiter,
            Phenomenon::Opposition,
            -5.3,
            (-6, 9, 15, 7.0),
            false,
        ),
        (
            Planet::Saturn,
            Phenomenon::Opposition,
            -5.3,
            (-6, 9, 14, 9.0),
            false,
        ),
        // three months before the discovery of Uranus
        (
            Planet::Uranus,
            Phenomenon::Opposition,
            1780.96,
            (1780, 12, 17, 14.0),
            true,
        ),
        // one month before the discovery of Neptune
        (
            Planet::Neptune,
            Phenomenon::Opposition,
            1846.63,
            (1846, 8, 20, 4.0),
            true,
        ),
    ];
    for (p, ph, year, (y, m, d, h), gregorian) in cases {
        let jde = phenomena::nearest_time_of_phenomenon(&p, ph, year).unwrap();
        assert_hour(jde, y, m, d, h, gregorian);
    }
}

// ---------------------------------------------------------------------------
// API behaviour and invariants
// ---------------------------------------------------------------------------

#[test]
fn availability() {
    for p in &ALL_PLANETS {
        for ph in ALL_PHENOMENA {
            let expected = expected_available(p, ph);
            assert_eq!(phenomena::is_available(p, ph), expected, "{ph:?}");
            assert_eq!(phenomena::approx_k(p, ph, 2000.0).is_some(), expected);
            assert_eq!(phenomena::nearest_k(p, ph, 2000.0).is_some(), expected);
            assert_eq!(phenomena::mean_phenomenon(p, ph, 0).is_some(), expected);
            assert_eq!(phenomena::time_of_phenomenon(p, ph, 0).is_some(), expected);
            assert_eq!(
                phenomena::nearest_time_of_phenomenon(p, ph, 2000.0).is_some(),
                expected
            );
        }
        for side in [Elongation::East, Elongation::West] {
            assert_eq!(
                phenomena::greatest_elongation(p, side, 0).is_some(),
                planet_index(p) <= 1
            );
        }
    }
}

#[test]
fn non_finite_year() {
    let p = Planet::Mars;
    let ph = Phenomenon::Opposition;
    assert_eq!(phenomena::nearest_k(&p, ph, f64::NAN), None);
    assert_eq!(phenomena::nearest_k(&p, ph, f64::INFINITY), None);
    assert_eq!(
        phenomena::nearest_time_of_phenomenon(&p, ph, f64::NAN),
        None
    );
    // So remote that k cannot be represented exactly: no answer, rather than
    // k = i64::MAX and a meaningless instant.
    for year in [1e20, -1e20, 1e306, f64::NEG_INFINITY] {
        assert_eq!(phenomena::nearest_k(&p, ph, year), None, "{year}");
        assert_eq!(
            phenomena::nearest_time_of_phenomenon(&p, ph, year),
            None,
            "{year}"
        );
    }
    // Far, but still exact
    let k = phenomena::nearest_k(&p, ph, 1e9).unwrap();
    let approx = phenomena::approx_k(&p, ph, 1e9).unwrap();
    assert!((k as f64 - approx).abs() < 1.0);
}

// k = 0 is the first configuration after 2000 January 1 (for the mean
// configuration, and for the true one too in practice).
#[test]
fn k_zero_is_first_after_2000() {
    let jan1 = jd(2000, 1, 1.0, true);
    for p in &ALL_PLANETS {
        for ph in [
            Phenomenon::InferiorConjunction,
            Phenomenon::SuperiorConjunction,
            Phenomenon::Opposition,
            Phenomenon::Conjunction,
        ] {
            let Some(mean0) = phenomena::mean_phenomenon(p, ph, 0) else {
                continue;
            };
            let mean_1 = phenomena::mean_phenomenon(p, ph, -1).unwrap();
            assert!(mean_1.jde < jan1 && mean0.jde >= jan1, "{ph:?}");
            let k = phenomena::approx_k(p, ph, 2000.0).unwrap();
            assert!((-1.0..0.0).contains(&k), "{ph:?} {k}");
        }
    }
}

#[test]
fn nearest_k_gives_the_nearest_occurrence() {
    for p in &ALL_PLANETS {
        for ph in ALL_PHENOMENA {
            if !phenomena::is_available(p, ph) {
                continue;
            }
            for year in [-1999.3, -6.0, 0.0, 1582.8, 1899.99, 2000.0, 2024.5, 3999.7] {
                let target = 365.2425 * year + 1721060.0;
                let k = phenomena::nearest_k(p, ph, year).unwrap();
                let t = |k| phenomena::time_of_phenomenon(p, ph, k).unwrap();
                let d = (t(k) - target).abs();
                assert!(d <= (t(k - 1) - target).abs(), "{ph:?} {year}");
                assert!(d <= (t(k + 1) - target).abs(), "{ph:?} {year}");
                assert_eq!(
                    phenomena::nearest_time_of_phenomenon(p, ph, year),
                    Some(t(k))
                );
            }
        }
    }
}

// Consecutive occurrences are one synodic period apart, within the size of
// the periodic corrections.
#[test]
fn successive_occurrences() {
    for p in &ALL_PLANETS {
        for ph in ALL_PHENOMENA {
            if !phenomena::is_available(p, ph) {
                continue;
            }
            let synodic = {
                let m0 = phenomena::mean_phenomenon(p, ph, 0).unwrap();
                let m1 = phenomena::mean_phenomenon(p, ph, 1).unwrap();
                m1.jde - m0.jde
            };
            for year in [-2000.0, -6.0, 1582.0, 1999.0, 2000.0, 2001.0, 3999.0] {
                let k = phenomena::approx_k(p, ph, year).unwrap().round() as i64;
                let a = phenomena::time_of_phenomenon(p, ph, k).unwrap();
                let b = phenomena::time_of_phenomenon(p, ph, k + 1).unwrap();
                assert!(b > a);
                // The largest correction terms (Mars) are about 40 days.
                assert!(((b - a) - synodic).abs() < 0.2 * synodic, "{ph:?} {k}");
            }
        }
    }
}

// The order of events around an inferior conjunction or an opposition.
#[test]
fn order_of_events() {
    use Phenomenon::*;
    let east = GreatestElongation(Elongation::East);
    let west = GreatestElongation(Elongation::West);
    let s1 = Station(phenomena::Station::First);
    let s2 = Station(phenomena::Station::Second);
    for year in [-1999.0, -6.0, 1582.0, 1999.0, 2000.0, 2001.0, 3999.0] {
        let k_of = |p: &Planet, ph| phenomena::approx_k(p, ph, year).unwrap().round() as i64;
        for p in [Planet::Mercury, Planet::Venus] {
            let k = k_of(&p, InferiorConjunction);
            let t = |ph| phenomena::time_of_phenomenon(&p, ph, k).unwrap();
            let next_sup = phenomena::time_of_phenomenon(&p, SuperiorConjunction, k + 1).unwrap();
            let seq = [
                t(SuperiorConjunction),
                t(east),
                t(s1),
                t(InferiorConjunction),
                t(s2),
                t(west),
                next_sup,
            ];
            assert!(seq.windows(2).all(|w| w[0] < w[1]), "{seq:?}");
        }
        for p in [Planet::Mars, Planet::Jupiter, Planet::Saturn] {
            let k = k_of(&p, Opposition);
            let t = |ph| phenomena::time_of_phenomenon(&p, ph, k).unwrap();
            let next_conj = phenomena::time_of_phenomenon(&p, Conjunction, k + 1).unwrap();
            let seq = [t(Conjunction), t(s1), t(Opposition), t(s2), next_conj];
            assert!(seq.windows(2).all(|w| w[0] < w[1]), "{seq:?}");
        }
        for p in [Planet::Uranus, Planet::Neptune] {
            let k = k_of(&p, Opposition);
            let t = |ph| phenomena::time_of_phenomenon(&p, ph, k).unwrap();
            let next_conj = phenomena::time_of_phenomenon(&p, Conjunction, k + 1).unwrap();
            let seq = [t(Conjunction), t(Opposition), next_conj];
            assert!(seq.windows(2).all(|w| w[0] < w[1]), "{seq:?}");
        }
    }
}

#[test]
fn mean_anomaly_range_and_wrap() {
    for p in &ALL_PLANETS {
        for ph in ALL_PHENOMENA {
            if !phenomena::is_available(p, ph) {
                continue;
            }
            for k in [-50_000, -3, 0, 3, 50_000] {
                let m = phenomena::mean_phenomenon(p, ph, k).unwrap();
                assert!((0.0..TAU).contains(&m.earth_mean_anomaly));
                assert_eq!(m.k, k);
                assert_close(m.julian_centuries, (m.jde - 2451545.0) / 36525.0, 1e-12);
                // M grows by M1 per period; check across the 0/360° wrap.
                let next = phenomena::mean_phenomenon(p, ph, k + 1).unwrap();
                let step = (next.earth_mean_anomaly - m.earth_mean_anomaly).rem_euclid(TAU);
                let synodic = next.jde - m.jde;
                let expected = (0.985600281 * synodic).to_radians().rem_euclid(TAU);
                assert_close(step, expected, 1e-6);
            }
        }
    }
}

// Greatest elongations of Mercury lie between about 18° and 28.2°, those of
// Venus between about 45° and 47.8°, over the whole valid span.
#[test]
fn elongation_ranges() {
    for (p, lo, hi) in [(Planet::Mercury, 17.5, 28.2), (Planet::Venus, 45.0, 47.8)] {
        let k0 = phenomena::approx_k(&p, Phenomenon::InferiorConjunction, -2000.0)
            .unwrap()
            .ceil() as i64;
        let k1 = phenomena::approx_k(&p, Phenomenon::InferiorConjunction, 4000.0)
            .unwrap()
            .floor() as i64;
        for k in (k0..=k1).step_by(7) {
            for side in [Elongation::East, Elongation::West] {
                let e = phenomena::greatest_elongation(&p, side, k).unwrap();
                let deg = e.elongation.to_degrees();
                assert!(deg > lo && deg < hi, "{side:?} k={k}: {deg}");
            }
        }
    }
}

// The eastern elongations of Mercury are largest when it is near
// aphelion; a sanity check of the sign of the large terms: the extreme
// values occur half a year apart and differ by several degrees.
#[test]
fn mercury_elongation_extremes() {
    let mut min = (f64::INFINITY, 0.0);
    let mut max = (0.0, 0.0);
    for k in -100..100 {
        let e = phenomena::greatest_elongation(&Planet::Mercury, Elongation::East, k).unwrap();
        let deg = e.elongation.to_degrees();
        if deg < min.0 {
            min = (deg, e.jde);
        }
        if deg > max.0 {
            max = (deg, e.jde);
        }
    }
    assert!(max.0 - min.0 > 9.0);
    let half_year = ((max.1 - min.1).rem_euclid(365.25) - 182.6).abs();
    assert!(half_year < 30.0, "{half_year}");
}

// ---------------------------------------------------------------------------
// Independent check with VSOP87
// ---------------------------------------------------------------------------

/// Numerical search for the true phenomena with the crate's VSOP87.
mod vsop {
    use super::*;

    /// Speed of light *| AU/day*
    const C: f64 = 173.144_632_674;

    fn rect(p: &Planet, jde: f64) -> [f64; 3] {
        let (l, b, r) = planet::heliocent_coords(p, jde);
        [r * b.cos() * l.cos(), r * b.cos() * l.sin(), r * b.sin()]
    }

    fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }

    fn norm(a: [f64; 3]) -> f64 {
        (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
    }

    fn unit(a: [f64; 3]) -> [f64; 3] {
        let n = norm(a);
        [a[0] / n, a[1] / n, a[2] / n]
    }

    /// Unit vector `u` displaced by the aberration of an observer moving at
    /// `v` AU/day
    fn aberrated(u: [f64; 3], v: [f64; 3]) -> [f64; 3] {
        unit([u[0] + v[0] / C, u[1] + v[1] / C, u[2] + v[2] / C])
    }

    /// Apparent geocentric directions of the Sun and of the planet, in the
    /// ecliptic frame of VSOP87D (nutation affects both equally and is
    /// ignored).
    pub fn directions(p: &Planet, jde: f64) -> ([f64; 3], [f64; 3]) {
        let earth = rect(&Planet::Earth, jde);
        let h = 0.01;
        let v = {
            let a = rect(&Planet::Earth, jde + h);
            let b = rect(&Planet::Earth, jde - h);
            [
                (a[0] - b[0]) / (2.0 * h),
                (a[1] - b[1]) / (2.0 * h),
                (a[2] - b[2]) / (2.0 * h),
            ]
        };
        let mut tau = 0.0;
        let mut u = [0.0; 3];
        for _ in 0..3 {
            u = sub(rect(p, jde - tau), earth);
            tau = norm(u) / C;
        }
        let sun = aberrated(unit([-earth[0], -earth[1], -earth[2]]), v);
        (sun, aberrated(unit(u), v))
    }

    fn lon(u: [f64; 3]) -> f64 {
        u[1].atan2(u[0])
    }

    fn wrap(x: f64) -> f64 {
        let y = x.rem_euclid(TAU);
        if y > PI { y - TAU } else { y }
    }

    pub fn elongation(p: &Planet, jde: f64) -> f64 {
        let (s, u) = directions(p, jde);
        let cross = [
            s[1] * u[2] - s[2] * u[1],
            s[2] * u[0] - s[0] * u[2],
            s[0] * u[1] - s[1] * u[0],
        ];
        norm(cross).atan2(s[0] * u[0] + s[1] * u[1] + s[2] * u[2])
    }

    /// Root of `f` near `t0`, by bisection in a bracket that is widened
    /// until it contains a sign change.
    fn root(f: impl Fn(f64) -> f64, t0: f64) -> f64 {
        let mut w = 0.5;
        loop {
            let (mut a, mut b) = (t0 - w, t0 + w);
            let fa = f(a);
            if fa.signum() != f(b).signum() {
                while b - a > 1e-5 {
                    let m = 0.5 * (a + b);
                    if f(m).signum() == fa.signum() {
                        a = m;
                    } else {
                        b = m;
                    }
                }
                return 0.5 * (a + b);
            }
            w *= 2.0;
            assert!(w < 20.0, "no event near {t0}");
        }
    }

    /// Maximum of the elongation near `t0`, by golden-section search
    pub fn greatest_elongation(p: &Planet, t0: f64) -> (f64, f64) {
        let g = 0.618_033_988_749_895;
        let (mut a, mut b) = (t0 - 3.0, t0 + 3.0);
        let f = |t| elongation(p, t);
        let (mut c, mut d) = (b - g * (b - a), a + g * (b - a));
        let (mut fc, mut fd) = (f(c), f(d));
        while b - a > 1e-4 {
            if fc > fd {
                b = d;
                d = c;
                fd = fc;
                c = b - g * (b - a);
                fc = f(c);
            } else {
                a = c;
                c = d;
                fc = fd;
                d = a + g * (b - a);
                fd = f(d);
            }
        }
        let t = 0.5 * (a + b);
        (t, f(t))
    }

    /// Instant of the phenomenon nearest `t0`
    pub fn event(p: &Planet, ph: Phenomenon, t0: f64) -> f64 {
        let dlon = |t: f64| {
            let (s, u) = directions(p, t);
            lon(u) - lon(s)
        };
        match ph {
            Phenomenon::Opposition => root(|t| wrap(dlon(t) - PI), t0),
            Phenomenon::Conjunction
            | Phenomenon::InferiorConjunction
            | Phenomenon::SuperiorConjunction => root(|t| wrap(dlon(t)), t0),
            Phenomenon::Station(_) => {
                let h = 0.02;
                root(
                    |t| wrap(lon(directions(p, t + h).1) - lon(directions(p, t - h).1)),
                    t0,
                )
            }
            Phenomenon::GreatestElongation(_) => greatest_elongation(p, t0).0,
        }
    }
}

/// Largest differences (hours) between the method and VSOP87 found over
/// A.D. 1800-2200, rounded up. Where Meeus states a smaller error (1 h for
/// Mercury and Venus, 2 h for Saturn, Uranus and Neptune, 3 h for Mars, 4 h
/// for Jupiter) that bound is used if it holds.
fn tolerance_hours(p: &Planet, ph: Phenomenon) -> f64 {
    use Phenomenon::*;
    match (planet_index(p), ph) {
        (0, _) => 1.0,
        (1, InferiorConjunction) => 1.0,
        (1, SuperiorConjunction) => 1.1,
        (1, Station(_)) => 1.2,
        (1, GreatestElongation(_)) => 2.1,
        (3, Conjunction) => 3.0,
        (3, Opposition) => 3.8,
        (3, Station(_)) => 5.7,
        (4, _) => 4.8,
        (5, Opposition | Conjunction) => 2.0,
        (5, Station(_)) => 4.4,
        (6, _) => 3.6,
        (7, _) => 2.1,
        _ => unreachable!(),
    }
}

/// Compares `count` occurrences spread over 1800-2200 with VSOP87.
fn check_against_vsop87(p: &Planet, count: i64) {
    for ph in ALL_PHENOMENA {
        if !phenomena::is_available(p, ph) {
            continue;
        }
        let k0 = phenomena::approx_k(p, ph, 1800.0).unwrap().ceil() as i64;
        let k1 = phenomena::approx_k(p, ph, 2200.0).unwrap().floor() as i64;
        let tol = tolerance_hours(p, ph) * HOUR;
        // An odd stride avoids sampling the same phase of the Earth's
        // anomaly over and over.
        let stride = ((k1 - k0) / count) | 1;
        let mut k = k0;
        while k <= k1 {
            let jde = phenomena::time_of_phenomenon(p, ph, k).unwrap();
            let accurate = vsop::event(p, ph, jde);
            assert!(
                (jde - accurate).abs() <= tol,
                "{ph:?} k={k}: off by {:.2} h",
                (jde - accurate) / HOUR
            );
            if let Phenomenon::GreatestElongation(side) = ph {
                let e = phenomena::greatest_elongation(p, side, k).unwrap();
                let value = vsop::elongation(p, accurate);
                assert_close(e.elongation.to_degrees(), value.to_degrees(), 0.01);
                // On the right side of the Sun
                let (s, u) = vsop::directions(p, accurate);
                let east = (u[1].atan2(u[0]) - s[1].atan2(s[0])).rem_euclid(TAU) < PI;
                assert_eq!(east, side == Elongation::East);
            }
            k += stride;
        }
    }
}

#[test]
fn vsop87_mercury() {
    check_against_vsop87(&Planet::Mercury, 8);
}

#[test]
fn vsop87_venus() {
    check_against_vsop87(&Planet::Venus, 8);
}

#[test]
fn vsop87_mars() {
    check_against_vsop87(&Planet::Mars, 8);
}

#[test]
fn vsop87_jupiter() {
    check_against_vsop87(&Planet::Jupiter, 8);
}

#[test]
fn vsop87_saturn() {
    check_against_vsop87(&Planet::Saturn, 8);
}

#[test]
fn vsop87_uranus_neptune() {
    check_against_vsop87(&Planet::Uranus, 8);
    check_against_vsop87(&Planet::Neptune, 8);
}

// Stations are stationary points of the apparent longitude, and
// conjunctions put the planet close to the Sun.
#[test]
fn vsop87_geometry_at_the_computed_instants() {
    let jde =
        phenomena::time_of_phenomenon(&Planet::Venus, Phenomenon::InferiorConjunction, 7).unwrap();
    // 2012 June 6: Venus crossed the Sun's disk; less than its radius away.
    assert!(vsop::elongation(&Planet::Venus, jde).to_degrees() < 0.25);

    let jde = phenomena::time_of_phenomenon(&Planet::Jupiter, Phenomenon::Opposition, 10).unwrap();
    assert!((PI - vsop::elongation(&Planet::Jupiter, jde)).to_degrees() < 1.5);

    // At each station the geocentric longitude moves by less than 1% of the
    // planet's mean daily motion, and in opposite directions a week before
    // and after (direct, then retrograde for station 1; the reverse for 2).
    let daily_motion = |p: &Planet, t: f64| {
        let lon = |t: f64| {
            let u = vsop::directions(p, t).1;
            u[1].atan2(u[0])
        };
        ((lon(t + 0.5) - lon(t - 0.5) + PI).rem_euclid(TAU) - PI).to_degrees()
    };
    for (p, k, mean_motion) in [
        (Planet::Mercury, 12, 4.09),
        (Planet::Venus, 3, 1.60),
        (Planet::Mars, 5, 0.524),
        (Planet::Jupiter, -30, 0.083),
        (Planet::Saturn, 50, 0.033),
    ] {
        for (station, before) in [(Station::First, 1.0), (Station::Second, -1.0)] {
            let jde = phenomena::time_of_phenomenon(&p, Phenomenon::Station(station), k).unwrap();
            assert!(
                daily_motion(&p, jde).abs() < 0.01 * mean_motion,
                "{p:?} {station:?}"
            );
            assert!(
                before * daily_motion(&p, jde - 7.0) > 0.0,
                "{p:?} {station:?}"
            );
            assert!(
                before * daily_motion(&p, jde + 7.0) < 0.0,
                "{p:?} {station:?}"
            );
        }
    }
}

/// Every occurrence over 1800-2200. Slow in debug builds; run with
/// `cargo test --release --test ch36_phenomena -- --ignored`.
#[test]
#[ignore]
fn vsop87_every_occurrence_1800_2200() {
    for p in [
        Planet::Mercury,
        Planet::Venus,
        Planet::Mars,
        Planet::Jupiter,
        Planet::Saturn,
        Planet::Uranus,
        Planet::Neptune,
    ] {
        check_against_vsop87(&p, i64::MAX);
    }
}

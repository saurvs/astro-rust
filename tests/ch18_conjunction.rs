// Meeus, chapter 18: planetary conjunctions.

mod common;

use astro::{angle, conjunction, coords};
use common::{assert_close, dms, hms};

/// Positions of a body at five instants, `(α, δ)` in radians
type Table = [(f64, f64); 5];

fn p(h: u32, m: u32, s: f64, neg: bool, d: u32, dm: u32, ds: f64) -> (f64, f64) {
    (hms(h, m, s).to_radians(), dms(neg, d, dm, ds).to_radians())
}

/// Mercury and Venus, apparent places at 0h TD, 1991 August 5 to 9
/// (Meeus, example 18.a)
fn mercury_venus() -> (Table, Table) {
    let mercury = [
        p(10, 24, 30.125, false, 6, 26, 32.05),
        p(10, 25, 0.342, false, 6, 10, 57.72),
        p(10, 25, 12.515, false, 5, 57, 33.08),
        p(10, 25, 6.235, false, 5, 46, 27.07),
        p(10, 24, 41.185, false, 5, 37, 48.45),
    ];
    let venus = [
        p(10, 27, 27.175, false, 4, 4, 41.83),
        p(10, 26, 32.410, false, 3, 55, 54.66),
        p(10, 25, 29.042, false, 3, 48, 3.51),
        p(10, 24, 17.191, false, 3, 41, 10.25),
        p(10, 22, 57.024, false, 3, 35, 16.61),
    ];
    (mercury, venus)
}

/// Vesta, J2000 positions at 0h TD on 1996 February 7, 12, 17, 22, 27, and
/// β Librae at J2000 corrected for proper motion to 1996.13 (Meeus, ch. 18,
/// exercise)
fn vesta_beta_lib() -> (Table, (f64, f64)) {
    let vesta = [
        p(15, 3, 51.937, true, 8, 57, 34.51),
        p(15, 9, 57.327, true, 9, 9, 3.88),
        p(15, 15, 37.898, true, 9, 17, 37.94),
        p(15, 20, 50.632, true, 9, 23, 16.25),
        p(15, 25, 32.695, true, 9, 26, 1.01),
    ];
    let star = p(15, 17, 0.446, true, 9, 22, 58.47);
    (vesta, star)
}

// Meeus, example 18.a: intermediate differences Δα (seconds of time) and
// Δδ (degrees)
#[test]
fn example_18_a_differences() {
    let (m, v) = mercury_venus();
    let d_a = [-177.050, -92.068, -16.527, 49.044, 104.161];
    let d_d = [2.363950, 2.250850, 2.158214, 2.088006, 2.042178];
    for i in 0..5 {
        let da = (m[i].0 - v[i].0).to_degrees() / 15.0 * 3600.0;
        assert_close(da, d_a[i], 5e-4);
        assert_close((m[i].1 - v[i].1).to_degrees(), d_d[i], 5e-7);
    }
}

// Meeus, example 18.a: n = +0.23797, 1991 August 7.23797 = 7 at 5h42.7m
// TD; Mercury 2.13940° (2°08′) north of Venus.
#[test]
fn example_18_a_conjunction() {
    let (m, v) = mercury_venus();
    let c = conjunction::conjunction(&m, &v).unwrap();
    assert_close(c.n, 0.23797, 5e-6);
    assert_close((c.n * 24.0 - 5.0) * 60.0, 42.7, 0.05);
    assert_close(c.lat_diff.to_degrees(), 2.13940, 5e-6);

    // The right ascensions are equal at that instant
    let s = conjunction::least_angular_separation(&m, &v).unwrap();
    // ... but that is not the instant of least separation, which comes later
    assert!(s.n > c.n + 0.1, "{} {}", s.n, c.n);
    assert!(s.sepr < c.lat_diff);

    // Swapping the bodies: same instant, opposite declination difference
    let c2 = conjunction::conjunction(&v, &m).unwrap();
    assert_close(c2.n, c.n, 1e-12);
    assert_close(c2.lat_diff, -c.lat_diff, 1e-15);

    // From the three central positions: nearly the same
    let c3 = conjunction::conjunction(&m[1..4], &v[1..4]).unwrap();
    assert_close(c3.n, c.n, 2e-3);
    assert_close(c3.lat_diff.to_degrees(), 2.1394, 2e-4);
}

// Meeus, ch. 18, exercise: Vesta passed 0°03′38″ north of β Lib on 1996
// February 18 at 6h37m TD. The star's position is given only once.
#[test]
fn vesta_beta_librae() {
    let (vesta, star) = vesta_beta_lib();
    let c = conjunction::conjunction(&vesta, &[star]).unwrap();
    let day = 17.0 + 5.0 * c.n;
    assert_close(day, 18.0 + (6.0 + 37.0 / 60.0) / 24.0, 0.5 / 1440.0);
    assert_close(c.lat_diff.to_degrees() * 3600.0, 3.0 * 60.0 + 38.0, 0.5);

    // Giving the star five times is the same
    let c5 = conjunction::conjunction(&vesta, &[star; 5]).unwrap();
    assert_eq!(c5, c);
    // Star first: opposite sign of the difference
    let c1 = conjunction::conjunction(&[star], &vesta).unwrap();
    assert_close(c1.n, c.n, 1e-12);
    assert_close(c1.lat_diff, -c.lat_diff, 1e-15);
}

// The proper-motion correction quoted in the exercise: from J2000 back
// 3.87 years with centennial proper motions −0.649 s and −1.91″.
#[test]
fn vesta_beta_librae_proper_motion() {
    let t = -0.0387;
    let asc = hms(15, 17, 0.421) + 15.0 * (-0.649 * t) / 3600.0;
    let dec = dms(true, 9, 22, 58.54) + (-1.91 * t) / 3600.0;
    assert_close(asc, hms(15, 17, 0.446), 15.0 * 0.0005 / 3600.0);
    assert_close(dec, dms(true, 9, 22, 58.47), 0.005 / 3600.0);
}

#[test]
fn conjunction_across_zero_hours() {
    let (m, v) = mercury_venus();
    let c = conjunction::conjunction(&m, &v).unwrap();
    // Rotate the right ascensions so that both bodies cross 0h between
    // August 6 and 8, and wrap them to [0, 2π) as an ephemeris would.
    let rot = (hms(10, 25, 10.0)).to_radians();
    let wrap = |b: &Table| b.map(|(a, d)| (angle::limit_to_two_PI(a - rot), d));
    let (mw, vw) = (wrap(&m), wrap(&v));
    assert!(mw.iter().any(|x| x.0 > 6.0) && mw.iter().any(|x| x.0 < 0.1));
    let cw = conjunction::conjunction(&mw, &vw).unwrap();
    assert_close(cw.n, c.n, 1e-11);
    assert_close(cw.lat_diff, c.lat_diff, 1e-14);
}

#[test]
fn conjunction_in_longitude_with_ecliptic_points() {
    // A planet moving uniformly in longitude passes a fixed point; with the
    // `EclPoint` type the latitude difference is the ecliptic one.
    let planet: Vec<coords::EclPoint> = (-1..=1)
        .map(|i| coords::EclPoint {
            long: (100.0 + 1.2 * (i as f64 + 0.4)).to_radians(),
            lat: (1.5 - 0.1 * i as f64).to_radians(),
        })
        .collect();
    let regulus = [coords::EclPoint {
        long: 100.0_f64.to_radians(),
        lat: 0.5_f64.to_radians(),
    }];
    let c = conjunction::conjunction(&planet, &regulus).unwrap();
    assert_close(c.n, -0.4, 1e-12);
    assert_close(c.lat_diff.to_degrees(), 1.0 + 0.04, 1e-10);

    // Southern mirror image: same instant, opposite latitude difference
    let planet_s: Vec<coords::EclPoint> = planet
        .iter()
        .map(|q| coords::EclPoint {
            long: q.long,
            lat: -q.lat,
        })
        .collect();
    let regulus_s = [coords::EclPoint {
        long: regulus[0].long,
        lat: -regulus[0].lat,
    }];
    let c_s = conjunction::conjunction(&planet_s, &regulus_s).unwrap();
    assert_close(c_s.n, c.n, 1e-12);
    assert_close(c_s.lat_diff, -c.lat_diff, 1e-15);
}

#[test]
fn conjunction_with_equatorial_points() {
    let (m, v) = mercury_venus();
    let to_eq = |b: &Table| b.map(|(asc, dec)| coords::EqPoint { asc, dec });
    let c = conjunction::conjunction(&to_eq(&m), &to_eq(&v)).unwrap();
    assert_close(c.n, 0.23797, 5e-6);
}

#[test]
fn conjunction_none_cases() {
    let (m, v) = mercury_venus();
    // No conjunction in August 5–7 (Δα changes sign only after the 7th)
    assert!(conjunction::conjunction(&m[0..3], &v[0..3]).is_none());
    // Invalid lengths
    assert!(conjunction::conjunction(&m[0..2], &v[0..2]).is_none());
    assert!(conjunction::conjunction(&m[0..3], &v).is_none());
    assert!(conjunction::conjunction(&m[0..1], &v[0..1]).is_none());
    // Non-finite input does not panic
    let mut bad = m;
    bad[1].0 = f64::NAN;
    let _ = conjunction::conjunction(&bad, &v);
}

// When two roots are within the table, the one nearest the central
// instant is returned.
#[test]
fn conjunction_nearest_to_centre() {
    // Δα = 0.01 (n − 0.5)(n + 1.5): roots at +0.5 and −1.5
    let body: Vec<(f64, f64)> = (-2..=2)
        .map(|i| {
            let n = i as f64;
            (1.0 + 0.01 * (n - 0.5) * (n + 1.5), 0.1 * n)
        })
        .collect();
    let c = conjunction::conjunction(&body, &[(1.0, 0.0)]).unwrap();
    assert_close(c.n, 0.5, 1e-12);
    assert_close(c.lat_diff, 0.05, 1e-12);
}

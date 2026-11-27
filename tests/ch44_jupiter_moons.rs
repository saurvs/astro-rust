//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 44: positions of the
//! satellites of Jupiter.

#![allow(non_snake_case)]

mod common;

use astro::planet::jupiter::moon::{self, ConjunctionKind, Moon};
use common::{assert_close, assert_close_deg};

/// JDE of example 44.a/44.b: 1992 December 16 at 0h UT with ΔT = 59 s.
/// The book prints the rounded value 2448972.50068; its results follow from
/// the unrounded one (Σ1 = −0.00654 needs it), so the tests use that.
const JDE_44: f64 = 2448972.5 + 59.0 / 86400.0;

fn deg_diff(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(360.0);
    if d > 180.0 { d - 360.0 } else { d }
}

// ---------------------------------------------------------------------------
// Low accuracy
// ---------------------------------------------------------------------------

// Meeus, example 43.b: the quantities of chapter 43 ("lower accuracy") used
// by the low-accuracy method.
#[test]
fn low_accuracy_auxiliary_quantities_example_43b() {
    let p = moon::low_accuracy_positions(2448972.50068);
    // The book's d − Δ/173 and Δ were computed from rounded intermediate
    // values, hence the slightly wider tolerances.
    assert_close(p.retarded_days, -2572.53205, 1e-5);
    assert_close(p.earth_dist, 5.66151, 2e-5);
    assert_close(p.jupiter_eq_of_center.to_degrees(), 1.235, 5e-4);
    assert_close(p.phase_angle.to_degrees(), -9.924, 5e-4);
    assert_close(p.phase_angle.sin(), -0.17234, 1e-5);
    assert_close(p.sun_declination.to_degrees(), -2.194, 5e-4);
    assert_close(p.earth_declination.to_degrees(), -2.50, 5e-3);
}

// Meeus, example 44.a
#[test]
fn low_accuracy_example_44a() {
    let p = moon::low_accuracy_positions(2448972.50068);

    // u before the corrections. The book derived them from the rounded
    // d − Δ/173 = −2572.53205 and ψ − B, so they agree to about 0.002°.
    let u0 = [324.542, 51.281, 4.651, 15.133];
    let correction = [-0.054, 1.063, 0.093, 0.541];
    let u = [324.488, 52.344, 4.744, 15.674];
    let r = [5.9299, 9.4016, 15.0062, 26.2140];
    let x = [-3.44, 7.44, 1.24, 7.08];
    let y = [0.21, 0.25, 0.65, 1.10];

    for m in Moon::ALL {
        let k = m.number() as usize - 1;
        let s = p.moon(m);
        assert_close_deg(s.u_uncorrected.to_degrees(), u0[k], 2e-3);
        assert_close(
            deg_diff(s.u.to_degrees(), s.u_uncorrected.to_degrees()),
            correction[k],
            5e-4,
        );
        assert_close_deg(s.u.to_degrees(), u[k], 2.5e-3);
        assert_close(s.r, r[k], 5e-5);
        assert_close(s.x, x[k], 5e-3);
        assert_close(s.y, y[k], 5e-3);

        // the original entry point returns the same values
        let (X, Y) = moon::apprnt_rect_coords(2448972.50068, &m);
        assert_eq!((X, Y), (s.x, s.y));
    }

    // arguments of the main periodic terms: 2(u1 − u2) and 2(u2 − u3)
    let u0 = |m: Moon| p.moon(m).u_uncorrected.to_degrees();
    assert_close_deg(2.0 * (u0(Moon::Io) - u0(Moon::Europa)), 186.52, 5e-3);
    assert_close_deg(2.0 * (u0(Moon::Europa) - u0(Moon::Ganymede)), 93.26, 5e-3);
}

#[test]
fn low_accuracy_angles_are_normalised() {
    for jde in [2415020.0, 2451545.0, 2460000.25, 2300000.0, 1000000.0] {
        let p = moon::low_accuracy_positions(jde);
        for s in p.moons.iter() {
            for a in [s.u, s.u_uncorrected] {
                assert!((0.0..std::f64::consts::TAU).contains(&a), "{a}");
            }
            assert!(s.x.abs() <= s.r && s.y.abs() <= s.r);
        }
        assert!(p.phase_angle.abs() < 12.5_f64.to_radians());
        assert!(p.earth_declination.abs() < 4.0_f64.to_radians());
        assert!(p.sun_declination.abs() <= 3.12_f64.to_radians() + 1e-12);
    }
}

// ---------------------------------------------------------------------------
// High accuracy
// ---------------------------------------------------------------------------

// Meeus, example 44.b
#[test]
fn high_accuracy_example_44b() {
    let p = moon::high_accuracy_positions(JDE_44);

    let sigma = [-0.00654, 1.10011, 0.04056, 0.59104];
    let x = [-3.4502, 7.4418, 1.2011, 7.0720];
    let y = [0.2137, 0.2753, 0.5900, 1.0291];

    for m in Moon::ALL {
        let k = m.number() as usize - 1;
        let s = p.moon(m);
        assert_close(s.longitude_perturbation.to_degrees(), sigma[k], 1e-5);
        // Meeus prints four decimals. Y of satellites II and IV comes out
        // 0.5 and 1.3 units of the last decimal lower than the book
        // (0.27525 and 1.02897). Every coefficient and rotation has been
        // checked against the book; the cause (probably the planetary
        // positions behind the book's own computation) was not found.
        assert_close(s.apparent.x, x[k], 1e-4);
        assert_close(s.apparent.y, y[k], 1.5e-4);
        // all four satellites are in front of the planet's centre plane
        assert!(s.apparent.z < 0.0);
        assert_eq!(moon::high_accuracy_rect_coords(JDE_44, m), s.apparent);
    }

    // With the rounded JDE printed in the book the results change by less
    // than 1e-4.
    let rounded = moon::high_accuracy_positions(2448972.50068);
    for k in 0..4 {
        assert_close(rounded.moons[k].apparent.x, x[k], 1.5e-4);
        assert_close(rounded.moons[k].apparent.y, y[k], 1.5e-4);
    }
}

#[test]
fn high_accuracy_light_time_and_distance() {
    let p = moon::high_accuracy_positions(JDE_44);
    // Example 43.a gives Δ = 5.6611239 for nearly the same light-time
    // corrected position of Jupiter.
    assert_close(p.earth_dist, 5.66112, 1e-4);
    assert_close(p.light_time, 0.0057755183 * p.earth_dist, 1e-9);
}

#[test]
fn high_accuracy_internal_consistency() {
    for jde in [2440000.5, 2448972.5, 2451545.0, 2455000.25, 2460000.75] {
        let p = moon::high_accuracy_positions(jde);
        for (k, s) in p.moons.iter().enumerate() {
            // the corrections change the length by a fraction of a percent
            for c in [s.apparent, s.shadow] {
                let n = (c.x * c.x + c.y * c.y + c.z * c.z).sqrt();
                assert!((n / s.radius - 1.0).abs() < 0.005, "{k} {n} {}", s.radius);
            }
            // Meeus: the latitudes reach 0°03′, 0°31′, 0°20′ and 0°44′
            let max_lat = [0.06, 0.6, 0.4, 0.8][k];
            assert!(s.latitude.to_degrees().abs() < max_lat);
            // mean distances ± the periodic terms
            let a = [5.90569, 9.39657, 14.98832, 26.36273][k];
            assert!((s.radius / a - 1.0).abs() < 0.012);
            // angles normalised
            for ang in [s.mean_longitude, s.longitude] {
                assert!((0.0..std::f64::consts::TAU).contains(&ang));
            }
            // longitude = mean longitude + Σ + precession since B1950
            let P = s.longitude - s.mean_longitude - s.longitude_perturbation;
            let P = P.to_degrees().rem_euclid(360.0);
            let T0 = (jde - 2433282.423) / 36525.0;
            assert_close(
                P,
                (1.3966626 * T0 + 0.0003088 * T0 * T0).rem_euclid(360.0),
                1e-9,
            );
        }
    }
}

// The low-accuracy method neglects the orbital inclinations and keeps one
// periodic term per satellite; over 1950-2050 it stays within a few
// hundredths to a few tenths of a radius of the E5 theory.
#[test]
fn low_and_high_accuracy_agree() {
    let max_dx = [0.06, 0.1, 0.15, 0.25];
    let max_dy = [0.02, 0.12, 0.12, 0.3];
    let mut jde = 2433282.5;
    while jde < 2469807.5 {
        let lo = moon::low_accuracy_positions(jde);
        let hi = moon::high_accuracy_positions(jde);
        for k in 0..4 {
            let dx = (lo.moons[k].x - hi.moons[k].apparent.x).abs();
            let dy = (lo.moons[k].y - hi.moons[k].apparent.y).abs();
            assert!(dx < max_dx[k], "satellite {} at {jde}: dx = {dx}", k + 1);
            assert!(dy < max_dy[k], "satellite {} at {jde}: dy = {dy}", k + 1);
        }
        jde += 37.3;
    }
}

// The coordinates are smooth functions of time: no jumps from angle
// normalisation or the light-time iteration.
#[test]
fn high_accuracy_is_continuous() {
    let h = 0.01;
    let start = 2451540.0;
    let mut prev2 = moon::high_accuracy_positions(start - h);
    let mut prev = moon::high_accuracy_positions(start);
    for n in 1..=1000 {
        let cur = moon::high_accuracy_positions(start + n as f64 * h);
        for k in 0..4 {
            for (a, b, c) in [
                (
                    prev2.moons[k].apparent.x,
                    prev.moons[k].apparent.x,
                    cur.moons[k].apparent.x,
                ),
                (
                    prev2.moons[k].apparent.y,
                    prev.moons[k].apparent.y,
                    cur.moons[k].apparent.y,
                ),
                (
                    prev2.moons[k].shadow.x,
                    prev.moons[k].shadow.x,
                    cur.moons[k].shadow.x,
                ),
                (
                    prev2.moons[k].shadow.y,
                    prev.moons[k].shadow.y,
                    cur.moons[k].shadow.y,
                ),
            ] {
                // |second difference| ≤ R ω² h², about 0.0075 for Io
                assert!(
                    (a - 2.0 * b + c).abs() < 0.01,
                    "step {n}, satellite {}",
                    k + 1
                );
            }
        }
        prev2 = prev;
        prev = cur;
    }
}

#[test]
fn high_accuracy_far_from_the_epoch_is_finite() {
    // 500 BC, 1000 AD, 3000 AD
    for jde in [1538432.5, 2086302.5, 2816787.5] {
        let p = moon::high_accuracy_positions(jde);
        for s in p.moons.iter() {
            for v in [
                s.apparent.x,
                s.apparent.y,
                s.apparent.z,
                s.shadow.x,
                s.shadow.y,
            ] {
                assert!(v.is_finite());
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Shadows and phenomena
// ---------------------------------------------------------------------------

// The views from the Sun and from the Earth differ by a rotation through
// Jupiter's phase angle φ (signed as in Meeus' ch. 43: negative before
// opposition), about an axis that is nearly Jupiter's rotation axis, plus a
// small tilt equal to D_E − D_S. So X_E ≈ X_S cos φ − Z_S sin φ,
// Z_E ≈ Z_S cos φ + X_S sin φ and Y_E ≈ Y_S + Z_S sin(D_E − D_S). The
// low-accuracy method supplies φ, D_E and D_S. Over 1990-2003 the relations
// hold to 0.11 % (X, Z) and 1.1 % (Y, which has second-order cross terms)
// of the orbital radius.
#[test]
fn shadow_coords_are_the_view_rotated_by_the_phase_angle() {
    let mut jde = 2448000.5;
    while jde < 2453000.0 {
        let lo = moon::low_accuracy_positions(jde);
        let (sp, cp) = lo.phase_angle.sin_cos();
        let tilt = (lo.earth_declination - lo.sun_declination).sin();
        let p = moon::high_accuracy_positions(jde);
        for s in p.moons.iter() {
            let (e, sh) = (s.apparent, s.shadow);
            assert!((e.x - (sh.x * cp - sh.z * sp)).abs() < 0.002 * s.radius);
            assert!((e.z - (sh.z * cp + sh.x * sp)).abs() < 0.002 * s.radius);
            assert!((e.y - (sh.y + sh.z * tilt)).abs() < 0.015 * s.radius);
        }
        jde += 11.3;
    }
}

// Before opposition (1992 December; Jupiter was at opposition on 1993
// March 30) the phase angle is negative and a shadow cast on the disk lies
// west of its satellite (it precedes the satellite across the disk); after
// opposition (1993 June) it lies east of it (it follows).
#[test]
fn shadow_precedes_satellite_before_opposition() {
    for (jde, before) in [(JDE_44, true), (2449153.5, false)] {
        assert_eq!(moon::low_accuracy_positions(jde).phase_angle < 0.0, before);
    }
    let lo = moon::low_accuracy_positions(2449076.5);
    assert!(lo.phase_angle.abs() < 0.3_f64.to_radians());

    // Find Io's shadow transits and compare the shadow with the satellite at
    // mid-transit of the shadow (shadow X = 0).
    for (start, before) in [(JDE_44, true), (2449153.5, false)] {
        let step = 5.0 / 1440.0;
        let mut checked = false;
        let mut prev = moon::high_accuracy_positions(start).moon(Moon::Io).shadow;
        for i in 1..(2.0 / step) as usize {
            let s = *moon::high_accuracy_positions(start + i as f64 * step).moon(Moon::Io);
            if s.shadow.z < 0.0 && prev.x < 0.0 && s.shadow.x >= 0.0 {
                // satellite in front, moving west; its shadow is west of it
                // before opposition and east of it after
                assert_eq!(s.shadow.x > s.apparent.x, before);
                assert!(s.phenomena().shadow_transit);
                checked = true;
            }
            prev = s.shadow;
        }
        assert!(checked);
    }
}

#[test]
fn disk_distance_accounts_for_flattening() {
    let c = moon::SatelliteCoords {
        x: 0.0,
        y: 0.93,
        z: -5.0,
    };
    assert!(c.is_on_disk());
    assert_close(c.disk_distance(), 0.93 * 1.071374, 1e-12);
    let c = moon::SatelliteCoords {
        x: 0.0,
        y: 0.94,
        z: -5.0,
    };
    assert!(!c.is_on_disk());
    let c = moon::SatelliteCoords {
        x: -0.99,
        y: 0.0,
        z: 5.0,
    };
    assert!(c.is_on_disk());
}

// Io passes through every phenomenon once per revolution (1.769 days). Its
// orbit lies in Jupiter's equator, so each lasts about the time Io needs to
// cross one Jupiter diameter, a little over two hours.
#[test]
fn io_phenomena_durations() {
    let start = JDE_44;
    let step = 2.0 / 1440.0;
    let n = (1.769 / step) as usize;
    let mut minutes = [0.0; 4];
    let mut changes = [0; 4];
    let mut last: Option<[bool; 4]> = None;
    for i in 0..=n {
        let p = moon::high_accuracy_positions(start + i as f64 * step);
        let ph = p.moon(Moon::Io).phenomena();
        let f = [ph.transit, ph.occultation, ph.eclipse, ph.shadow_transit];
        for k in 0..4 {
            if f[k] {
                minutes[k] += 2.0;
            }
            if let Some(l) = last {
                if l[k] != f[k] {
                    changes[k] += 1;
                }
            }
        }
        // a satellite cannot be in front of and behind the disk at once
        assert!(!(ph.transit && ph.occultation));
        assert!(!(ph.eclipse && ph.shadow_transit));
        last = Some(f);
    }
    for k in 0..4 {
        assert!(
            (120.0..150.0).contains(&minutes[k]),
            "{k}: {} min",
            minutes[k]
        );
        assert!(changes[k] <= 2);
    }
}

// ---------------------------------------------------------------------------
// Conjunctions
// ---------------------------------------------------------------------------

// Meeus, exercise of ch. 44: on 1988 November 23 satellite III was in
// inferior conjunction at 7h28m UT with Y = −0.8043 (so it transited the
// disk), and satellite IV in superior conjunction at 5h15m UT with
// Y = +1.3991 (so it passed north of the disk, not occulted). ΔT = 55.8 s
// from Meeus' table 10.A.
#[test]
fn conjunctions_exercise_1988() {
    let delta_t = 55.8 / 86400.0;
    let day = 2447488.5; // 1988 Nov 23, 0h UT
    let start = day + delta_t;

    let iii = moon::conjunctions_with_jupiter(Moon::Ganymede, start, start + 1.0);
    assert_eq!(iii.len(), 1);
    assert_eq!(iii[0].kind, ConjunctionKind::Inferior);
    let ut = (iii[0].jde - delta_t - day) * 1440.0;
    assert_close(ut, (7 * 60 + 28) as f64, 0.5);
    assert_close(iii[0].y, -0.8043, 1.5e-4);
    let p = moon::high_accuracy_positions(iii[0].jde);
    assert!(p.moon(Moon::Ganymede).apparent.z < 0.0);
    assert!(p.moon(Moon::Ganymede).phenomena().transit);
    assert_close(p.moon(Moon::Ganymede).apparent.x, 0.0, 1e-6);

    let iv = moon::conjunctions_with_jupiter(Moon::Callisto, start, start + 1.0);
    assert_eq!(iv.len(), 1);
    assert_eq!(iv[0].kind, ConjunctionKind::Superior);
    let ut = (iv[0].jde - delta_t - day) * 1440.0;
    assert_close(ut, (5 * 60 + 15) as f64, 0.5);
    assert_close(iv[0].y, 1.3991, 1.5e-4);
    let p = moon::high_accuracy_positions(iv[0].jde);
    assert!(p.moon(Moon::Callisto).apparent.z > 0.0);
    assert!(!p.moon(Moon::Callisto).phenomena().occultation);
}

#[test]
fn conjunctions_alternate() {
    // Io: about two conjunctions per 1.77 days, alternately inferior and
    // superior, with Z of the matching sign.
    let c = moon::conjunctions_with_jupiter(Moon::Io, JDE_44, JDE_44 + 10.0);
    assert!((10..=12).contains(&c.len()), "{}", c.len());
    for w in c.windows(2) {
        assert_ne!(w[0].kind, w[1].kind);
        assert_close(w[1].jde - w[0].jde, 1.769 / 2.0, 0.05);
    }
    for e in c.iter() {
        let z = moon::high_accuracy_rect_coords(e.jde, Moon::Io).z;
        match e.kind {
            ConjunctionKind::Inferior => assert!(z < 0.0),
            ConjunctionKind::Superior => assert!(z > 0.0),
        }
    }
}

#[test]
fn conjunction_search_bad_intervals() {
    assert!(moon::conjunctions_with_jupiter(Moon::Io, JDE_44 + 1.0, JDE_44).is_empty());
    assert!(moon::conjunctions_with_jupiter(Moon::Io, f64::NAN, JDE_44).is_empty());
    assert!(moon::mutual_conjunctions(Moon::Io, Moon::Io, JDE_44, JDE_44 + 5.0).is_empty());
    assert!(moon::mutual_conjunctions(Moon::Io, Moon::Europa, JDE_44, f64::INFINITY).is_empty());
}

#[test]
fn mutual_conjunctions_are_consistent() {
    let a = moon::mutual_conjunctions(Moon::Io, Moon::Europa, JDE_44, JDE_44 + 4.0);
    let b = moon::mutual_conjunctions(Moon::Europa, Moon::Io, JDE_44, JDE_44 + 4.0);
    assert!(!a.is_empty());
    assert_eq!(a.len(), b.len());
    for (c, d) in a.iter().zip(b.iter()) {
        assert_close(c.jde, d.jde, 1e-6);
        assert_close(c.y_separation, -d.y_separation, 1e-6);
        assert_close(c.z_difference, -d.z_difference, 1e-6);

        let p = moon::high_accuracy_positions(c.jde);
        let io = p.moon(Moon::Io).apparent;
        let eu = p.moon(Moon::Europa).apparent;
        assert_close(io.x, eu.x, 1e-5);
        assert_close(c.x, io.x, 1e-5);
        assert_close(c.y_separation, io.y - eu.y, 1e-9);
    }
}

#[test]
fn moon_numbers() {
    let n: Vec<u8> = Moon::ALL.iter().map(|m| m.number()).collect();
    assert_eq!(n, vec![1, 2, 3, 4]);
}

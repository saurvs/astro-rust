/*
Copyright (c) 2015, 2016 Saurav Sachidanand

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

//! Tests for `transit` (Meeus, chapter 15: rising, transit and setting).

#![allow(non_snake_case)]

mod common;

use astro::coords::{EqPoint, GeographPoint};
use astro::transit::{self, NoEvent, TransitBody, TransitType};
use astro::{coords, interpol};
use std::f64::consts::{PI, TAU};

fn deg(d: f64) -> f64 {
    d.to_radians()
}

fn eq(asc_deg: f64, dec_deg: f64) -> EqPoint {
    EqPoint {
        asc: deg(asc_deg),
        dec: deg(dec_deg),
    }
}

// Meeus, example 15.a: Venus at Boston, 1988 March 20.
fn boston() -> GeographPoint {
    GeographPoint {
        long: deg(71.0833),
        lat: deg(42.3333),
    }
}

fn venus() -> [EqPoint; 3] {
    [
        eq(40.68021, 18.04761),
        eq(41.73129, 18.44092),
        eq(42.78204, 18.82742),
    ]
}

const THETA0_15A: f64 = 177.74208; // degrees, 1988 March 20, 0h UT
const DELTA_T_15A: f64 = 56.0;

fn venus_m(tt: TransitType, geo: &GeographPoint, theta0_deg: f64) -> Result<f64, NoEvent> {
    let v = venus();
    transit::day_fraction(
        &tt,
        transit::STD_ALT_STAR_OR_PLANET,
        geo,
        &v[0],
        &v[1],
        &v[2],
        deg(theta0_deg),
        DELTA_T_15A,
    )
}

/// Altitude of the body at day fraction `m`, interpolated as in Meeus'
/// method, for checking that a returned time really is a rise or set.
fn alt_at(m: f64, geo: &GeographPoint, p: &[EqPoint; 3], theta0: f64, delta_t: f64) -> f64 {
    let n = m + delta_t / 86400.0;
    let a2 = p[1].asc;
    let unwrap = |a: f64| a2 + (a - a2 + PI).rem_euclid(TAU) - PI;
    let asc = interpol::three_values(unwrap(p[0].asc), a2, unwrap(p[2].asc), n);
    let dec = interpol::three_values(p[0].dec, p[1].dec, p[2].dec, n);
    let theta = theta0 + deg(360.985647) * m;
    coords::alt_frm_eq(
        coords::hr_angl_frm_observer_long(theta, geo.long, asc),
        dec,
        geo.lat,
    )
}

// Meeus, example 15.a: final values. The book stops after one correction;
// its next corrections (-0.000003, -0.000004, -0.000004) are below the
// printed precision.
#[test]
fn day_fraction_example_15_a() {
    let geo = boston();
    let rise = venus_m(TransitType::Rise, &geo, THETA0_15A).unwrap();
    let transit = venus_m(TransitType::Transit, &geo, THETA0_15A).unwrap();
    let set = venus_m(TransitType::Set, &geo, THETA0_15A).unwrap();

    common::assert_close(rise, 0.51766, 1e-5); // 12h25m UT
    common::assert_close(transit, 0.81980, 1e-5); // 19h41m UT
    common::assert_close(set, 0.12130, 1e-5); // 2h55m UT
    common::assert_close(rise, 0.51766 - 0.000003, 2e-6);
    common::assert_close(transit, 0.81980 - 0.000004, 2e-6);
    common::assert_close(set, 0.12130 - 0.000004, 2e-6);

    // At the computed rise and set the altitude is the standard altitude.
    let v = venus();
    for m in [rise, set] {
        let h = alt_at(m, &geo, &v, deg(THETA0_15A), DELTA_T_15A);
        common::assert_close(h, transit::STD_ALT_STAR_OR_PLANET, 1e-7);
    }
}

// Meeus, example 15.a, with the deprecated (hour, minute, second) API.
#[test]
#[allow(deprecated)]
fn time_example_15_a() {
    let v = venus();
    let geo = boston();
    let deltaT = astro::time::delta_t(1988, 3);
    let t = |tt| {
        transit::time(
            &tt,
            &TransitBody::StarOrPlanet,
            &geo,
            &v[0],
            &v[1],
            &v[2],
            deg(THETA0_15A),
            deltaT,
            0.0,
        )
    };

    let (h, m, s) = t(TransitType::Rise);
    assert_eq!((h, m), (12, 25));
    common::assert_close(s, 25.6, 1.0);
    let (h, m, _) = t(TransitType::Transit);
    assert_eq!((h, m), (19, 40)); // 19h40.5m, which the book rounds to 19h41m
    let (h, m, _) = t(TransitType::Set);
    assert_eq!((h, m), (2, 54)); // 2h54.7m, rounded to 2h55m in the book
}

// Bug: the old code took acos of a value outside [-1, 1] and returned
// garbage for circumpolar bodies. Now rises and sets report why they do not
// happen, and the transit is still computed.
#[test]
fn circumpolar_and_never_rising_bodies() {
    let geo = boston(); // 42.33° N
    let star = |dec: f64| [eq(100.0, dec), eq(100.0, dec), eq(100.0, dec)];
    let m = |tt, p: &[EqPoint; 3]| {
        transit::day_fraction(
            &tt,
            transit::STD_ALT_STAR_OR_PLANET,
            &geo,
            &p[0],
            &p[1],
            &p[2],
            deg(THETA0_15A),
            DELTA_T_15A,
        )
    };

    let capella_like = star(60.0);
    assert_eq!(
        m(TransitType::Rise, &capella_like),
        Err(NoEvent::AlwaysAbove)
    );
    assert_eq!(
        m(TransitType::Set, &capella_like),
        Err(NoEvent::AlwaysAbove)
    );
    let tr = m(TransitType::Transit, &capella_like).unwrap();
    assert!((0.0..1.0).contains(&tr));

    let southern = star(-60.0);
    assert_eq!(m(TransitType::Rise, &southern), Err(NoEvent::AlwaysBelow));
    assert_eq!(m(TransitType::Set, &southern), Err(NoEvent::AlwaysBelow));
    assert!(m(TransitType::Transit, &southern).is_ok());

    // The limit is 90° - φ + h0 = 47.10°: the lower culmination is then at
    // the standard altitude.
    assert!(m(TransitType::Rise, &star(47.0)).is_ok());
    assert!(m(TransitType::Set, &star(47.0)).is_ok());
    assert_eq!(m(TransitType::Rise, &star(47.2)), Err(NoEvent::AlwaysAbove));
}

#[test]
#[allow(deprecated)]
fn time_of_circumpolar_body_is_nan() {
    let p = eq(100.0, 70.0);
    let (h, m, s) = transit::time(
        &TransitType::Rise,
        &TransitBody::StarOrPlanet,
        &boston(),
        &p,
        &p,
        &p,
        deg(THETA0_15A),
        DELTA_T_15A,
        0.0,
    );
    assert_eq!((h, m), (0, 0));
    assert!(s.is_nan());
}

// Bug: right ascensions stepping across 0h between the three days were
// interpolated through 12h, giving a nonsense position.
#[test]
#[allow(deprecated)]
fn right_ascension_across_0h() {
    let geo = boston();
    let wrapped = [eq(359.0, 5.0), eq(0.0, 5.4), eq(1.0, 5.8)];
    let unwrapped = [eq(-1.0, 5.0), eq(0.0, 5.4), eq(1.0, 5.8)];
    for tt in [TransitType::Rise, TransitType::Transit, TransitType::Set] {
        let m = |p: &[EqPoint; 3]| {
            transit::day_fraction(
                &tt,
                transit::STD_ALT_SUN,
                &geo,
                &p[0],
                &p[1],
                &p[2],
                deg(THETA0_15A),
                DELTA_T_15A,
            )
            .unwrap()
        };
        common::assert_close(m(&wrapped), m(&unwrapped), 1e-9);

        // The same through the deprecated API (this failed before the fix).
        let t = |p: &[EqPoint; 3]| {
            transit::time(
                &tt,
                &TransitBody::Sun,
                &geo,
                &p[0],
                &p[1],
                &p[2],
                deg(THETA0_15A),
                DELTA_T_15A,
                0.0,
            )
        };
        let (a, b) = (t(&wrapped), t(&unwrapped));
        assert_eq!((a.0, a.1), (b.0, b.1));
        common::assert_close(a.2, b.2, 1e-3);
    }

    // Wrapping the other way: 358°, 359°, 0°.
    let wrapped = [eq(358.0, 5.0), eq(359.0, 5.4), eq(0.0, 5.8)];
    let shifted = [eq(358.0, 5.0), eq(359.0, 5.4), eq(360.0, 5.8)];
    let m = |p: &[EqPoint; 3]| {
        transit::day_fraction(
            &TransitType::Transit,
            0.0,
            &geo,
            &p[0],
            &p[1],
            &p[2],
            deg(THETA0_15A),
            DELTA_T_15A,
        )
        .unwrap()
    };
    common::assert_close(m(&wrapped), m(&shifted), 1e-9);
}

// Bug: the approximate times were brought into [0, 1) by adding or
// subtracting 1 only once (and 1.0 itself was kept), so angles given
// outside [0, 2π) produced negative or > 24 h times.
#[test]
#[allow(deprecated)]
fn sidereal_time_and_right_ascension_in_any_range() {
    let geo = boston();
    for tt in [TransitType::Rise, TransitType::Transit, TransitType::Set] {
        let reference = venus_m(tt, &geo, THETA0_15A).unwrap();
        for k in [-3.0, -2.0, -1.0, 1.0, 2.0, 5.0] {
            let m = venus_m(tt, &geo, THETA0_15A + 360.0 * k).unwrap();
            common::assert_close(m, reference, 1e-9);
        }
        // Right ascensions shifted by -360° (as atan2 may return).
        let v = venus();
        let shifted: Vec<EqPoint> = v
            .iter()
            .map(|p| EqPoint {
                asc: p.asc - TAU,
                dec: p.dec,
            })
            .collect();
        let m = transit::day_fraction(
            &tt,
            transit::STD_ALT_STAR_OR_PLANET,
            &geo,
            &shifted[0],
            &shifted[1],
            &shifted[2],
            deg(THETA0_15A),
            DELTA_T_15A,
        )
        .unwrap();
        common::assert_close(m, reference, 1e-9);

        let (h, _, _) = transit::time(
            &tt,
            &TransitBody::StarOrPlanet,
            &geo,
            &v[0],
            &v[1],
            &v[2],
            deg(THETA0_15A + 720.0),
            DELTA_T_15A,
            0.0,
        );
        assert!((0..24).contains(&h), "hour {h}");
    }
}

// Observers east of Greenwich (negative longitude) and south of the
// equator, including the date line: every returned time is a real rise,
// transit or set of the interpolated body.
#[test]
fn eastern_and_southern_observers() {
    let v = venus();
    let theta0 = deg(THETA0_15A);
    for (lat, long) in [
        (42.3333, 71.0833),    // Boston
        (-33.8688, -151.2093), // Sydney
        (35.6762, -139.6503),  // Tokyo
        (-18.1416, -178.4419), // Suva, near the date line
        (-13.8333, 171.75),    // Apia, just across it
        (64.1466, 21.9426),    // Reykjavik
        (-54.8019, 68.3030),   // Ushuaia
        (0.0, 180.0),
        (0.0, -180.0),
    ] {
        let geo = GeographPoint {
            long: deg(long),
            lat: deg(lat),
        };
        for tt in [TransitType::Rise, TransitType::Transit, TransitType::Set] {
            let m = transit::day_fraction(
                &tt,
                transit::STD_ALT_STAR_OR_PLANET,
                &geo,
                &v[0],
                &v[1],
                &v[2],
                theta0,
                DELTA_T_15A,
            )
            .unwrap();
            assert!((-0.01..1.01).contains(&m), "{lat} {long} {tt:?}: {m}");
            if tt == TransitType::Transit {
                // Hour angle zero at transit.
                let n = m + DELTA_T_15A / 86400.0;
                let asc = interpol::three_values(v[0].asc, v[1].asc, v[2].asc, n);
                let ha =
                    coords::hr_angl_frm_observer_long(theta0 + deg(360.985647) * m, geo.long, asc);
                common::assert_close_rad(ha, 0.0, 1e-7);
            } else {
                let h = alt_at(m, &geo, &v, theta0, DELTA_T_15A);
                common::assert_close(h, transit::STD_ALT_STAR_OR_PLANET, 1e-7);
            }
        }
    }

    // Longitudes of +180° and -180° are the same meridian.
    let at = |long: f64| {
        let geo = GeographPoint {
            long: deg(long),
            lat: deg(-17.0),
        };
        venus_m(TransitType::Rise, &geo, THETA0_15A).unwrap()
    };
    common::assert_close(at(180.0), at(-180.0), 1e-9);
    common::assert_close(at(179.9999), at(-179.9999), 1e-5);
}

#[test]
fn standard_altitudes() {
    common::assert_close_deg(
        transit::standard_altitude(&TransitBody::StarOrPlanet, 0.0).to_degrees(),
        -34.0 / 60.0,
        1e-12,
    );
    common::assert_close_deg(
        transit::standard_altitude(&TransitBody::Sun, 0.0).to_degrees(),
        -50.0 / 60.0,
        1e-12,
    );
    // Meeus, example 47.a: π = 0°59'31".2 gives h0 = +0.155°.
    let moon = transit::standard_altitude(&TransitBody::Moon, deg(0.991990));
    common::assert_close_deg(moon.to_degrees(), 0.7275 * 0.991990 - 34.0 / 60.0, 1e-12);
    // The mean value of the book.
    common::assert_close_deg(transit::STD_ALT_MOON_MEAN.to_degrees(), 0.125, 1e-12);
}

#[test]
fn hour_angle_at_altitude() {
    // On the equator a body on the equator spends half the day above the
    // geometric horizon.
    common::assert_close(
        transit::hour_angle_at_altitude(0.0, 0.0, 0.0).unwrap(),
        PI / 2.0,
        1e-12,
    );
    // Meeus, example 15.a: H0 = 108.5344° (with the book's h0 = -0.5667°).
    let h0 = transit::hour_angle_at_altitude(deg(-0.5667), deg(42.3333), deg(18.44092)).unwrap();
    common::assert_close_deg(h0.to_degrees(), 108.5344, 1e-4);
    // Symmetry: (φ, δ) and (-φ, -δ) give the same arc.
    let a = transit::hour_angle_at_altitude(deg(-0.8), deg(35.0), deg(20.0)).unwrap();
    let b = transit::hour_angle_at_altitude(deg(-0.8), deg(-35.0), deg(-20.0)).unwrap();
    common::assert_close(a, b, 1e-12);
    // Summer and winter arcs add up to a day for h0 = 0.
    let s = transit::hour_angle_at_altitude(0.0, deg(50.0), deg(23.44)).unwrap();
    let w = transit::hour_angle_at_altitude(0.0, deg(50.0), deg(-23.44)).unwrap();
    common::assert_close(s + w, PI, 1e-12);

    // At the poles a body is either always up or always down.
    assert_eq!(
        transit::hour_angle_at_altitude(0.0, PI / 2.0, deg(10.0)),
        Err(NoEvent::AlwaysAbove)
    );
    assert_eq!(
        transit::hour_angle_at_altitude(0.0, -PI / 2.0, deg(10.0)),
        Err(NoEvent::AlwaysBelow)
    );
    assert_eq!(
        transit::hour_angle_at_altitude(f64::NAN, 0.0, 0.0),
        Err(NoEvent::NoConvergence)
    );
}

#[test]
fn twilight_altitudes() {
    // Civil twilight starts before sunrise and ends after sunset.
    let geo = boston();
    let v = venus();
    let m = |tt, h0: f64| {
        transit::day_fraction(
            &tt,
            deg(h0),
            &geo,
            &v[0],
            &v[1],
            &v[2],
            deg(THETA0_15A),
            56.0,
        )
        .unwrap()
    };
    assert!(m(TransitType::Rise, -6.0) < m(TransitType::Rise, -50.0 / 60.0));
    assert!(m(TransitType::Set, -6.0) > m(TransitType::Set, -50.0 / 60.0));
    // The transit does not depend on h0.
    common::assert_close(
        m(TransitType::Transit, -6.0),
        m(TransitType::Transit, 0.0),
        1e-12,
    );
}

#[test]
fn non_finite_input() {
    let v = venus();
    for tt in [TransitType::Rise, TransitType::Transit, TransitType::Set] {
        let r = transit::day_fraction(
            &tt,
            transit::STD_ALT_SUN,
            &boston(),
            &v[0],
            &v[1],
            &v[2],
            f64::NAN,
            56.0,
        );
        assert_eq!(r, Err(NoEvent::NoConvergence));
    }
    // NoEvent is a std error with a message.
    let e: Box<dyn std::error::Error> = Box::new(NoEvent::AlwaysAbove);
    assert!(!e.to_string().is_empty());
}

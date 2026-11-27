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

//! Tests for `observer`: apparent and topocentric positions, and rise,
//! transit and set (Meeus, chs. 12, 13, 15, 16, 23, 25, 33, 40, 47).

mod common;

use astro::observer::{
    self, ApparentPosition, Body, Horizon, HorizonStatus, Observer, RiseTransitSet,
};
use astro::{coords, ecliptic, nutation, time, transit};
use common::{arcsec, assert_close, assert_close_deg, dms, hms};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

const SECOND: f64 = 1.0 / 86400.0;
const MINUTE: f64 = 60.0 * SECOND;

fn deg(d: f64) -> f64 {
    d.to_radians()
}

/// Julian Day of 0h UT on a Gregorian date.
fn jd0(year: i16, month: time::Month, day: u8) -> f64 {
    time::julian_day(&time::Date {
        year,
        month,
        decimal_day: day as f64,
        cal_type: time::CalType::Gregorian,
    })
}

/// ΔT for the Julian Day, as the observer module computes it.
fn delta_t(jd: f64) -> f64 {
    let year = 2000.0 + (jd - 2451544.5) / 365.2425;
    let y = year.floor();
    time::delta_t(y as i32, (((year - y) * 12.0) as u8).min(11) + 1)
}

/// Apparent sidereal time at Greenwich (Meeus, ch. 12).
fn gast(jd_ut: f64, jde: f64) -> f64 {
    let (dpsi, deps) = nutation::nutation(jde);
    time::mn_sidr(jd_ut) + dpsi * (ecliptic::mn_oblq_IAU(jde) + deps).cos()
}

/// Geocentric altitude of the centre of a body (no parallax, no
/// refraction), and the standard altitude for it, computed afresh at an
/// instant (independently of the interpolation used by rise_transit_set).
fn geocentric_alt(body: Body, obs: &Observer, jd_ut: f64) -> (f64, f64) {
    let jde = jd_ut + delta_t(jd_ut) / 86400.0;
    let mut pos = observer::geocentric_apparent_position(body, jde);
    let h0 = match body {
        Body::Sun => transit::STD_ALT_SUN,
        Body::Moon => transit::standard_altitude(&transit::TransitBody::Moon, pos.eq_hz_parallax),
        _ => transit::STD_ALT_STAR_OR_PLANET,
    };
    pos.eq_hz_parallax = 0.0;
    let t = observer::topocentric_from_geocentric(&pos, obs, gast(jd_ut, jde));
    (t.altitude, h0)
}

/// Checks that every event returned is what it claims to be, and inside the
/// search window.
#[track_caller]
fn check_events(body: Body, obs: &Observer, jd: f64, horizon: Horizon, r: &RiseTransitSet) {
    let h0_of = |h0_std: f64| match horizon {
        Horizon::Standard => h0_std,
        Horizon::CivilTwilight => deg(-6.0),
        Horizon::NauticalTwilight => deg(-12.0),
        Horizon::AstronomicalTwilight => deg(-18.0),
        Horizon::Altitude(h) => h,
    };
    for (name, t, dir) in [("rise", r.rise, 1.0), ("set", r.set, -1.0)] {
        if let Some(t) = t {
            assert!(
                (jd..jd + 1.0).contains(&t),
                "{name} {t} outside window {jd}"
            );
            let (h, h0) = geocentric_alt(body, obs, t);
            let h0 = h0_of(h0);
            assert!(
                (h - h0).abs() < arcsec(2.0).to_radians(),
                "{body:?} {name} at {t}: altitude {}° vs {}°",
                h.to_degrees(),
                h0.to_degrees()
            );
            // The body is moving in the right direction.
            let (later, _) = geocentric_alt(body, obs, t + MINUTE);
            assert!(dir * (later - h) > 0.0, "{body:?} {name}: wrong direction");
        }
    }
    if let Some(t) = r.transit {
        assert!((jd..jd + 1.0).contains(&t));
        let pos = observer::topocentric_position(body, obs, t);
        let jde = t + delta_t(t) / 86400.0;
        let geo = observer::geocentric_apparent_position(body, jde);
        let ha = coords::hr_angl_frm_observer_long(gast(t, jde), obs.long, geo.asc);
        common::assert_close_rad(ha, 0.0, deg(arcsec(3.0)));
        assert!(pos.hour_angle.abs() < deg(0.1));
    }
    match r.status {
        HorizonStatus::RisesOrSets => assert!(r.rise.is_some() || r.set.is_some()),
        HorizonStatus::AlwaysAbove | HorizonStatus::AlwaysBelow => {
            assert_eq!((r.rise, r.set), (None, None));
            let above = r.status == HorizonStatus::AlwaysAbove;
            for k in 0..=8 {
                let t = jd + k as f64 / 8.0;
                let (h, h0) = geocentric_alt(body, obs, t);
                assert_eq!(h >= h0_of(h0), above, "{body:?} at {t}");
            }
        }
    }
}

fn london() -> Observer {
    Observer::from_degrees_east(51.5074, -0.1278, 0.0)
}

fn sydney() -> Observer {
    Observer::from_degrees_east(-33.8688, 151.2093, 0.0)
}

// ---------------------------------------------------------------------------
// Apparent geocentric positions: the worked examples of the book.

// Meeus, example 25.b (and the full-VSOP87 values quoted in example 25.a):
// the Sun on 1992 October 13.0 TD.
#[test]
fn sun_example_25_b() {
    let sun = observer::geocentric_apparent_position(Body::Sun, 2448908.5);
    assert_close_deg(
        sun.asc.to_degrees(),
        hms(13, 13, 30.749),
        0.01 * 15.0 / 3600.0,
    );
    assert_close(sun.dec.to_degrees(), dms(true, 7, 47, 1.74), arcsec(0.05));
    // R = 0.99760775 with the abridged VSOP87 of example 25.b; the crate has
    // the full theory, which gives 0.99760853 (example 25.a).
    assert_close(sun.distance, 0.99760853, 1e-8);
    assert_close(
        sun.eq_hz_parallax.to_degrees(),
        arcsec(8.794 / 0.99760853),
        arcsec(1e-3),
    );
}

// Meeus, example 33.a: Venus on 1992 December 20 at 0h TD.
#[test]
fn venus_example_33_a() {
    let venus = observer::geocentric_apparent_position(Body::Venus, 2448976.5);
    assert_close_deg(
        venus.asc.to_degrees(),
        hms(21, 4, 41.454),
        0.005 * 15.0 / 3600.0,
    );
    assert_close(
        venus.dec.to_degrees(),
        dms(true, 18, 53, 16.84),
        arcsec(0.05),
    );
    assert_close(venus.distance, 0.910947, 2e-6);
}

// Meeus, example 47.a: the Moon on 1992 April 12 at 0h TD.
#[test]
fn moon_example_47_a() {
    let moon = observer::geocentric_apparent_position(Body::Moon, 2448724.5);
    assert_close_deg(moon.asc.to_degrees(), 134.688470, 5e-6);
    assert_close(moon.dec.to_degrees(), 13.768368, 5e-6);
    assert_close(moon.distance * 149_597_870.7, 368409.7, 0.1);
    assert_close(moon.eq_hz_parallax.to_degrees(), 0.991990, 1e-6);
}

// Meeus, example 23.a: θ Persei on 2028 November 13.19 TD. The J2000.0
// position of example 21.b is first moved by the star's proper motion.
#[test]
fn star_example_23_a() {
    let jde = 2462088.69;
    let years = (jde - 2451545.0) / 365.25;
    let asc = hms(2, 44, 11.986 + 0.03425 * years);
    let dec = dms(false, 49, 13, 42.48 - 0.0895 * years);
    let star = observer::geocentric_apparent_position(
        Body::Star {
            asc: deg(asc),
            dec: deg(dec),
        },
        jde,
    );
    assert_close_deg(
        star.asc.to_degrees(),
        hms(2, 46, 14.390),
        0.005 * 15.0 / 3600.0,
    );
    assert_close(
        star.dec.to_degrees(),
        dms(false, 49, 21, 7.45),
        arcsec(0.05),
    );
    assert_eq!(star.distance, f64::INFINITY);
    assert_eq!(star.eq_hz_parallax, 0.0);
}

// Meeus, example 13.b: Venus from the US Naval Observatory, 1987 April 10,
// 19h21m UT. The book interpolates α and δ from an ephemeris; we compute
// them, then check the horizontal coordinates from the book's inputs.
#[test]
fn venus_example_13_b() {
    let jd = 2446896.30625;
    let jde = jd + time::delta_t(1987, 4) / 86400.0;
    let venus = observer::geocentric_apparent_position(Body::Venus, jde);
    assert_close_deg(
        venus.asc.to_degrees(),
        hms(23, 9, 16.641),
        0.01 * 15.0 / 3600.0,
    );
    assert_close(venus.dec.to_degrees(), dms(true, 6, 43, 11.61), arcsec(0.1));

    // Example 12.b/13.b: apparent sidereal time 8h34m56.853s.
    assert_close_deg(
        gast(jd, jde).to_degrees(),
        hms(8, 34, 56.853),
        0.001 * 15.0 / 3600.0,
    );

    // Horizontal coordinates from the book's values. The book uses the
    // longitude 5h08m15.7s; the parallax (6.9") is left out as in the book.
    let book = ApparentPosition {
        asc: deg(hms(23, 9, 16.641)),
        dec: deg(dms(true, 6, 43, 11.61)),
        distance: venus.distance,
        eq_hz_parallax: 0.0,
    };
    let usno = Observer::new(deg(dms(false, 38, 55, 17.0)), deg(hms(5, 8, 15.7)), 0.0);
    let pos = observer::topocentric_from_geocentric(&book, &usno, deg(hms(8, 34, 56.853)));
    assert_close(pos.hour_angle.to_degrees(), 64.352133, 1e-6);
    assert_close(pos.azimuth_from_south.to_degrees(), 68.0337, 5e-5);
    assert_close(pos.altitude.to_degrees(), 15.1249, 5e-5);
    assert_close(pos.azimuth_from_north.to_degrees(), 68.0337 + 180.0, 5e-5);

    // With the parallax the planet is a few arcseconds lower.
    let topo = observer::topocentric_position(Body::Venus, &usno, jd);
    let drop = (15.1249 - topo.altitude.to_degrees()) * 3600.0;
    assert!((5.0..8.0).contains(&drop), "{drop}");
}

// Meeus, example 40.a: Mars from Palomar, 2003 August 28, 3h17m UT.
#[test]
fn mars_example_40_a() {
    let palomar = Observer::new(deg(dms(false, 33, 21, 22.0)), deg(hms(7, 47, 27.0)), 1706.0);
    // From the book's geocentric values.
    let geo = ApparentPosition {
        asc: deg(339.530208),
        dec: deg(-15.771083),
        distance: 0.37276,
        eq_hz_parallax: deg(arcsec(23.592)),
    };
    let pos = observer::topocentric_from_geocentric(&geo, &palomar, deg(hms(1, 40, 45.0)));
    assert_close_deg(
        pos.asc.to_degrees(),
        hms(22, 38, 8.54),
        0.01 * 15.0 / 3600.0,
    );
    assert_close(pos.dec.to_degrees(), dms(true, 15, 46, 30.0), arcsec(0.1));

    // The whole pipeline from the instant.
    let jd = 2452879.5 + (3.0 + 17.0 / 60.0) / 24.0;
    let jde = jd + time::delta_t(2003, 8) / 86400.0;
    let mars = observer::geocentric_apparent_position(Body::Mars, jde);
    assert_close_deg(mars.asc.to_degrees(), 339.530208, arcsec(0.5));
    assert_close(mars.dec.to_degrees(), -15.771083, arcsec(0.5));
    assert_close(mars.distance, 0.37276, 1e-5);
    assert_close(
        mars.eq_hz_parallax.to_degrees(),
        arcsec(23.592),
        arcsec(0.001),
    );
    let topo = observer::topocentric_position(Body::Mars, &palomar, jd);
    assert_close_deg(
        topo.asc.to_degrees(),
        hms(22, 38, 8.54),
        0.03 * 15.0 / 3600.0,
    );
    assert_close(topo.dec.to_degrees(), dms(true, 15, 46, 30.0), arcsec(0.5));
}

// Meeus, example 28.a: the equation of time on 1992 October 13.0 TD is
// +13m42.6s, i.e. the apparent solar time at Greenwich (the Sun's hour
// angle + 12h) is ahead of UT by that much. The book uses the mean Sun of
// TT; at UT it differs by ΔT × 0.986°/day ≈ 0.2 s.
#[test]
fn equation_of_time_example_28_a() {
    let jd = 2448908.5 - time::delta_t(1992, 10) / 86400.0;
    let sun = observer::topocentric_position(Body::Sun, &Observer::new(0.0, 0.0, 0.0), jd);
    let solar_time = (sun.hour_angle / TAU + 0.5).rem_euclid(1.0);
    let ut = (jd - 0.5).rem_euclid(1.0);
    let e = (solar_time - ut).rem_euclid(1.0) * 1440.0; // minutes
    assert_close(e, 13.0 + 42.6 / 60.0, 0.5 / 60.0);
}

// ---------------------------------------------------------------------------
// Rise, transit and set.

// Meeus, example 15.a, but with every input computed: Venus at Boston on
// 1988 March 20. The book's answers are 0.51766, 0.81980 and 0.12130 of the
// day (12h25m, 19h41m, 2h55m UT).
#[test]
fn venus_example_15_a() {
    let boston = Observer::new(deg(42.3333), deg(71.0833), 0.0);
    let jd = jd0(1988, time::Month::Mar, 20);
    let r = observer::rise_transit_set(Body::Venus, &boston, jd, Horizon::Standard);
    assert_close(r.rise.unwrap() - jd, 0.51766, 2e-5);
    assert_close(r.transit.unwrap() - jd, 0.81980, 2e-5);
    assert_close(r.set.unwrap() - jd, 0.12130, 2e-5);
    assert_eq!(r.status, HorizonStatus::RisesOrSets);
    check_events(Body::Venus, &boston, jd, Horizon::Standard, &r);

    // The two methods agree.
    let gp = boston.geograph_point();
    let p = |k: f64| {
        let a = observer::geocentric_apparent_position(Body::Venus, jd + k);
        coords::EqPoint {
            asc: a.asc,
            dec: a.dec,
        }
    };
    let (p1, p2, p3) = (p(-1.0), p(0.0), p(1.0));
    let m = transit::day_fraction(
        &transit::TransitType::Rise,
        transit::STD_ALT_STAR_OR_PLANET,
        &gp,
        &p1,
        &p2,
        &p3,
        gast(jd, jd),
        time::delta_t(1988, 3),
    )
    .unwrap();
    assert_close(r.rise.unwrap() - jd, m, 2.0 * SECOND);
}

/// Local clock time (hours) of an event, for a UTC offset in hours.
fn local_hours(t: f64, jd: f64, utc_offset: f64) -> f64 {
    (t - jd) * 24.0 + utc_offset
}

// Published sunrise and sunset times for London (rounded to the minute):
// 2024 June 21, 04:43 and 21:21 BST; 2024 December 21, 08:04 and 15:53 GMT.
#[test]
fn london_solstices() {
    let obs = london();
    for (jd, utc, rise, set) in [
        (
            jd0(2024, time::Month::June, 21),
            1.0,
            4.0 + 43.0 / 60.0,
            21.0 + 21.0 / 60.0,
        ),
        (
            jd0(2024, time::Month::Dec, 21),
            0.0,
            8.0 + 4.0 / 60.0,
            15.0 + 53.0 / 60.0,
        ),
    ] {
        let start = jd - utc / 24.0; // local midnight
        let r = observer::rise_transit_set(Body::Sun, &obs, start, Horizon::Standard);
        assert_close(local_hours(r.rise.unwrap(), jd, utc), rise, 1.0 / 60.0);
        assert_close(local_hours(r.set.unwrap(), jd, utc), set, 1.0 / 60.0);
        check_events(Body::Sun, &obs, start, Horizon::Standard, &r);
    }
}

// Southern and eastern: published times for Sydney, 2024 December 21
// (05:41 and 20:05 AEDT, UTC+11) and June 21 (07:00 and 16:54 AEST, UTC+10).
#[test]
fn sydney_solstices() {
    let obs = sydney();
    for (jd, utc, rise, set) in [
        (
            jd0(2024, time::Month::Dec, 21),
            11.0,
            5.0 + 41.0 / 60.0,
            20.0 + 5.0 / 60.0,
        ),
        (
            jd0(2024, time::Month::June, 21),
            10.0,
            7.0,
            16.0 + 54.0 / 60.0,
        ),
    ] {
        let start = jd - utc / 24.0;
        let r = observer::rise_transit_set(Body::Sun, &obs, start, Horizon::Standard);
        assert_close(local_hours(r.rise.unwrap(), jd, utc), rise, 1.0 / 60.0);
        assert_close(local_hours(r.set.unwrap(), jd, utc), set, 1.0 / 60.0);
        check_events(Body::Sun, &obs, start, Horizon::Standard, &r);

        // The Sun culminates in the north as seen from Sydney.
        let noon = observer::topocentric_position(Body::Sun, &obs, r.transit.unwrap());
        assert_close_deg(noon.azimuth_from_north.to_degrees(), 0.0, 0.01);
        assert_close_deg(noon.azimuth_from_south.to_degrees(), 180.0, 0.01);
    }
}

// Rise and set agree with the semidiurnal arc of eq. 15.1, evaluated with
// the Sun's declination at the moment of rising or setting, at latitudes
// from 60° S to 60° N, at the solstices and an equinox.
#[test]
fn rise_and_set_match_semidiurnal_arc() {
    for jd in [
        jd0(2024, time::Month::June, 20),
        jd0(2024, time::Month::Dec, 21),
        jd0(2024, time::Month::Mar, 20),
    ] {
        for lat in (-60..=60).step_by(10) {
            // Longitude chosen so that local noon is near the middle of the
            // UT day.
            let obs = Observer::new(deg(lat as f64), deg(5.0), 0.0);
            let r = observer::rise_transit_set(Body::Sun, &obs, jd, Horizon::Standard);
            let (rise, transit, set) = (r.rise.unwrap(), r.transit.unwrap(), r.set.unwrap());
            assert!(rise < transit && transit < set);

            let sun = |t: f64| {
                observer::geocentric_apparent_position(Body::Sun, t + delta_t(t) / 86400.0)
            };
            // Rate of the Sun's hour angle, per day.
            let ha_rate =
                deg(360.985647) - (sun(transit + 0.5).asc - sun(transit - 0.5).asc).rem_euclid(TAU);
            let arc = |t: f64| {
                transit::hour_angle_at_altitude(transit::STD_ALT_SUN, obs.lat, sun(t).dec).unwrap()
            };
            assert_close(rise, transit - arc(rise) / ha_rate, 5.0 * SECOND);
            assert_close(set, transit + arc(set) / ha_rate, 5.0 * SECOND);
        }
    }

    // On the equator at the equinox the day is 12h + 2 × 50′ ≈ 12h07m.
    let obs = Observer::new(0.0, 0.0, 0.0);
    let r = observer::rise_transit_set(
        Body::Sun,
        &obs,
        jd0(2024, time::Month::Mar, 20),
        Horizon::Standard,
    );
    let hours = (r.set.unwrap() - r.rise.unwrap()) * 24.0;
    assert!((12.1..12.12).contains(&hours), "{hours}");
}

#[test]
fn polar_day_and_night() {
    let june = jd0(2024, time::Month::June, 21);
    let december = jd0(2024, time::Month::Dec, 21);
    let cases = [
        // (name, lat, east long, status in June, status in December)
        ("Tromsø", 69.6492, 18.9553, true, false),
        ("Longyearbyen", 78.2232, 15.6267, true, false),
        ("North Pole", 90.0, 0.0, true, false),
        ("McMurdo", -77.8419, 166.6863, false, true),
        ("South Pole", -90.0, 0.0, false, true),
    ];
    let status = |up: bool| {
        if up {
            HorizonStatus::AlwaysAbove
        } else {
            HorizonStatus::AlwaysBelow
        }
    };
    for (name, lat, long, june_up, december_up) in cases {
        let obs = Observer::from_degrees_east(lat, long, 0.0);
        for (jd, up) in [(june, june_up), (december, december_up)] {
            let r = observer::rise_transit_set(Body::Sun, &obs, jd, Horizon::Standard);
            assert_eq!(r.status, status(up), "{name} {jd}");
            assert_eq!((r.rise, r.set), (None, None), "{name}");
            check_events(Body::Sun, &obs, jd, Horizon::Standard, &r);
        }
    }

    // Even in the polar night at Tromsø there is civil twilight at noon, but
    // at Longyearbyen not even that.
    let obs = Observer::from_degrees_east(69.6492, 18.9553, 0.0);
    let r = observer::rise_transit_set(Body::Sun, &obs, december, Horizon::CivilTwilight);
    assert_eq!(r.status, HorizonStatus::RisesOrSets);
    check_events(Body::Sun, &obs, december, Horizon::CivilTwilight, &r);
    let obs = Observer::from_degrees_east(78.2232, 15.6267, 0.0);
    let r = observer::rise_transit_set(Body::Sun, &obs, december, Horizon::CivilTwilight);
    assert_eq!(r.status, HorizonStatus::AlwaysBelow);
}

// At the poles the altitude of the Sun is ± its declination, so the Sun
// rises once a year, when the declination reaches ∓0°50′.
#[test]
fn sunrise_at_the_poles() {
    for (lat, month, day, rise_dec) in [
        (-90.0, time::Month::Sept, 15, 50.0 / 60.0),
        (90.0, time::Month::Mar, 13, -50.0 / 60.0),
    ] {
        let obs = Observer::new(deg(lat), deg(-30.0), 0.0);
        let start = jd0(2024, month, day);
        let mut seen = Vec::new();
        for d in 0..12 {
            let jd = start + d as f64;
            let r = observer::rise_transit_set(Body::Sun, &obs, jd, Horizon::Standard);
            check_events(Body::Sun, &obs, jd, Horizon::Standard, &r);
            if let Some(t) = r.rise {
                let sun =
                    observer::geocentric_apparent_position(Body::Sun, t + delta_t(t) / 86400.0);
                assert_close(sun.dec.to_degrees(), rise_dec, arcsec(2.0));
            }
            seen.push(r.status);
        }
        // Night, one day with a sunrise, then day.
        let first_rise = seen
            .iter()
            .position(|s| *s == HorizonStatus::RisesOrSets)
            .unwrap();
        assert!(
            seen[..first_rise]
                .iter()
                .all(|s| *s == HorizonStatus::AlwaysBelow)
        );
        assert!(
            seen[first_rise + 1..]
                .iter()
                .all(|s| *s == HorizonStatus::AlwaysAbove)
        );
    }
}

// Around the start of the midnight Sun at Tromsø, the last days have a
// sunset after midnight but no sunrise, and so on: every event must still
// be correct, and once the Sun stays up it stays up.
#[test]
fn polar_circle_transition() {
    let obs = Observer::from_degrees_east(69.6492, 18.9553, 0.0);
    let start = jd0(2024, time::Month::May, 10);
    let mut always_up = false;
    for d in 0..25 {
        let jd = start + d as f64;
        let r = observer::rise_transit_set(Body::Sun, &obs, jd, Horizon::Standard);
        check_events(Body::Sun, &obs, jd, Horizon::Standard, &r);
        if always_up {
            assert_eq!(r.status, HorizonStatus::AlwaysAbove, "day {d}");
        }
        always_up |= r.status == HorizonStatus::AlwaysAbove;
    }
    assert!(always_up);
}

// The Moon at a mid-latitude: about once a month there is no moonrise, and
// once no moonset, on a calendar day. Every event found is genuine, and the
// classical method of chapter 15 agrees.
#[test]
fn moon_over_a_month() {
    let boston = Observer::from_degrees_east(42.3601, -71.0589, 0.0);
    let gp = boston.geograph_point();
    let start = jd0(2024, time::Month::Jan, 1);
    let (mut no_rise, mut no_set) = (0, 0);
    for d in 0..30 {
        let jd = start + d as f64;
        let r = observer::rise_transit_set(Body::Moon, &boston, jd, Horizon::Standard);
        assert_eq!(r.status, HorizonStatus::RisesOrSets);
        check_events(Body::Moon, &boston, jd, Horizon::Standard, &r);
        no_rise += r.rise.is_none() as u32;
        no_set += r.set.is_none() as u32;

        let p = |k: f64| {
            let a = observer::geocentric_apparent_position(Body::Moon, jd + k);
            coords::EqPoint {
                asc: a.asc,
                dec: a.dec,
            }
        };
        let (p1, p2, p3) = (p(-1.0), p(0.0), p(1.0));
        let parallax = observer::geocentric_apparent_position(Body::Moon, jd).eq_hz_parallax;
        let h0 = transit::standard_altitude(&transit::TransitBody::Moon, parallax);
        for (tt, ours) in [
            (transit::TransitType::Rise, r.rise),
            (transit::TransitType::Transit, r.transit),
            (transit::TransitType::Set, r.set),
        ] {
            let m = transit::day_fraction(&tt, h0, &gp, &p1, &p2, &p3, gast(jd, jd), delta_t(jd))
                .unwrap();
            match ours {
                Some(t) => assert_close(t - jd, m, 30.0 * SECOND),
                // No event that day: Meeus' method lands on the next day.
                None => assert!(m >= 1.0, "{tt:?} on day {d}: m = {m}"),
            }
        }
    }
    assert!((1..=2).contains(&no_rise), "{no_rise}");
    assert!((1..=2).contains(&no_set), "{no_set}");
}

// Near the 2024-2025 major lunar standstill the Moon's declination reaches
// ±28.5°, so at Tromsø it is circumpolar for some days and does not rise at
// all for others.
#[test]
fn moon_at_high_latitude() {
    let obs = Observer::from_degrees_east(69.6492, 18.9553, 0.0);
    let start = jd0(2025, time::Month::Jan, 1);
    let mut seen = [false; 3];
    for d in 0..30 {
        let jd = start + d as f64;
        let r = observer::rise_transit_set(Body::Moon, &obs, jd, Horizon::Standard);
        check_events(Body::Moon, &obs, jd, Horizon::Standard, &r);
        seen[match r.status {
            HorizonStatus::RisesOrSets => 0,
            HorizonStatus::AlwaysAbove => 1,
            HorizonStatus::AlwaysBelow => 2,
        }] = true;
    }
    assert_eq!(seen, [true; 3]);
}

#[test]
fn stars_and_planets() {
    let obs = london();
    let jd = jd0(2024, time::Month::Jan, 15);
    let star = |ra: (u32, u32, f64), dec: (bool, u32, u32, f64)| Body::Star {
        asc: deg(hms(ra.0, ra.1, ra.2)),
        dec: deg(dms(dec.0, dec.1, dec.2, dec.3)),
    };
    let polaris = star((2, 31, 49.09), (false, 89, 15, 50.8));
    let canopus = star((6, 23, 57.11), (true, 52, 41, 44.4));
    let sirius = star((6, 45, 8.917), (true, 16, 42, 58.02));

    let r = observer::rise_transit_set(polaris, &obs, jd, Horizon::Standard);
    assert_eq!(r.status, HorizonStatus::AlwaysAbove);
    assert!(r.transit.is_some());
    check_events(polaris, &obs, jd, Horizon::Standard, &r);

    let r = observer::rise_transit_set(canopus, &obs, jd, Horizon::Standard);
    assert_eq!(r.status, HorizonStatus::AlwaysBelow);
    check_events(canopus, &obs, jd, Horizon::Standard, &r);

    let r = observer::rise_transit_set(sirius, &obs, jd, Horizon::Standard);
    assert!(r.rise.is_some() && r.set.is_some() && r.transit.is_some());
    check_events(sirius, &obs, jd, Horizon::Standard, &r);
    // Sirius culminates at 90° − 51.51° − 16.72° ≈ 21.8°.
    let top = observer::topocentric_position(sirius, &obs, r.transit.unwrap());
    assert_close(top.altitude.to_degrees(), 90.0 - 51.5074 - 16.74, 0.05);

    for body in [
        Body::Mercury,
        Body::Venus,
        Body::Mars,
        Body::Jupiter,
        Body::Saturn,
        Body::Uranus,
        Body::Neptune,
    ] {
        for horizon in [Horizon::Standard, Horizon::Altitude(deg(10.0))] {
            let r = observer::rise_transit_set(body, &obs, jd, horizon);
            check_events(body, &obs, jd, horizon, &r);
            assert!(r.transit.is_some(), "{body:?}");
        }
    }
}

#[test]
fn twilights_are_ordered() {
    let obs = london();
    let jd = jd0(2024, time::Month::Mar, 20);
    let r = |h| observer::rise_transit_set(Body::Sun, &obs, jd, h);
    let (astro, naut, civil, std) = (
        r(Horizon::AstronomicalTwilight),
        r(Horizon::NauticalTwilight),
        r(Horizon::CivilTwilight),
        r(Horizon::Standard),
    );
    for (h, x) in [
        (Horizon::AstronomicalTwilight, &astro),
        (Horizon::NauticalTwilight, &naut),
        (Horizon::CivilTwilight, &civil),
    ] {
        check_events(Body::Sun, &obs, jd, h, x);
        assert_eq!(x.transit, std.transit);
    }
    let order = [
        astro.rise,
        naut.rise,
        civil.rise,
        std.rise,
        std.transit,
        std.set,
        civil.set,
        naut.set,
        astro.set,
    ];
    let order: Vec<f64> = order.iter().map(|t| t.unwrap()).collect();
    assert!(order.windows(2).all(|w| w[0] < w[1]), "{order:?}");
    // At the equinox in London civil twilight lasts about 33 minutes.
    let civil_minutes = (std.rise.unwrap() - civil.rise.unwrap()) * 1440.0;
    assert!((30.0..36.0).contains(&civil_minutes), "{civil_minutes}");

    // Around the June solstice the Sun never gets 18° below London's
    // horizon: there is no astronomical night.
    let june = jd0(2024, time::Month::June, 21);
    let a = observer::rise_transit_set(Body::Sun, &obs, june, Horizon::AstronomicalTwilight);
    assert_eq!(a.status, HorizonStatus::AlwaysAbove);
    let n = observer::rise_transit_set(Body::Sun, &obs, june, Horizon::NauticalTwilight);
    assert_eq!(n.status, HorizonStatus::RisesOrSets);
}

// Observers on either side of the date line, with the search window on
// their local calendar day.
#[test]
fn date_line_observers() {
    let jd = jd0(2024, time::Month::Mar, 1);
    // ±180° are the same meridian; ±179.9999° are 0.0002° apart.
    let a = observer::rise_transit_set(
        Body::Sun,
        &Observer::from_degrees_east(-17.0, 180.0, 0.0),
        jd,
        Horizon::Standard,
    );
    let b = observer::rise_transit_set(
        Body::Sun,
        &Observer::from_degrees_east(-17.0, -180.0, 0.0),
        jd,
        Horizon::Standard,
    );
    assert_close(a.rise.unwrap(), b.rise.unwrap(), 1e-6 * SECOND);
    assert_close(a.set.unwrap(), b.set.unwrap(), 1e-6 * SECOND);
    let c = observer::rise_transit_set(
        Body::Sun,
        &Observer::from_degrees_east(-17.0, 179.9999, 0.0),
        jd,
        Horizon::Standard,
    );
    let d = observer::rise_transit_set(
        Body::Sun,
        &Observer::from_degrees_east(-17.0, -179.9999, 0.0),
        jd,
        Horizon::Standard,
    );
    assert_close(c.rise.unwrap(), d.rise.unwrap(), 0.1 * SECOND);

    // Suva (UTC+12) and Apia (UTC+13), and Honolulu (UTC−10): on their local
    // day, the Sun rises in the morning and sets in the evening, local time.
    for (lat, long, utc) in [
        (-18.1416, 178.4419, 12.0),
        (-13.8333, -171.75, 13.0),
        (21.3069, -157.8583, -10.0),
    ] {
        let obs = Observer::from_degrees_east(lat, long, 0.0);
        let start = jd - utc / 24.0;
        for body in [Body::Sun, Body::Moon] {
            let r = observer::rise_transit_set(body, &obs, start, Horizon::Standard);
            check_events(body, &obs, start, Horizon::Standard, &r);
        }
        let r = observer::rise_transit_set(Body::Sun, &obs, start, Horizon::Standard);
        let rise = local_hours(r.rise.unwrap(), jd, utc);
        let set = local_hours(r.set.unwrap(), jd, utc);
        assert!(
            (5.0..8.0).contains(&rise) && (17.0..20.0).contains(&set),
            "{rise} {set}"
        );
    }
}

// Dates long before the Gregorian reform, with the large ΔT of antiquity.
#[test]
fn ancient_dates() {
    let athens = Observer::from_degrees_east(37.9715, 23.7257, 0.0);
    for (year, month, day) in [(-500, time::Month::June, 21), (1000, time::Month::Dec, 15)] {
        let jd = time::julian_day(&time::Date {
            year,
            month,
            decimal_day: day as f64,
            cal_type: time::CalType::Julian,
        });
        for body in [Body::Sun, Body::Moon, Body::Venus] {
            let r = observer::rise_transit_set(body, &athens, jd, Horizon::Standard);
            check_events(body, &athens, jd, Horizon::Standard, &r);
            assert!(r.transit.is_some() || body == Body::Moon);
        }
    }
}

// ---------------------------------------------------------------------------
// Topocentric positions.

#[test]
fn azimuth_conventions_and_ranges() {
    let observers = [
        london(),
        sydney(),
        Observer::from_degrees_east(0.0, 0.0, 0.0),
        Observer::from_degrees_east(-89.9, 140.0, 2800.0),
        Observer::from_degrees_east(64.1, -21.9, 100.0),
    ];
    let bodies = [
        Body::Sun,
        Body::Moon,
        Body::Mars,
        Body::Star {
            asc: deg(279.23473),
            dec: deg(38.78369),
        },
    ];
    for obs in &observers {
        for body in bodies {
            for k in 0..12 {
                let jd = 2460400.5 + k as f64 * 0.37;
                let p = observer::topocentric_position(body, obs, jd);
                assert!((0.0..TAU).contains(&p.azimuth_from_south));
                assert!((0.0..TAU).contains(&p.azimuth_from_north));
                assert!((0.0..TAU).contains(&p.asc));
                assert!((-PI..PI).contains(&p.hour_angle));
                assert!(p.altitude.abs() <= FRAC_PI_2 && p.dec.abs() <= FRAC_PI_2);
                common::assert_close_rad(p.azimuth_from_north, p.azimuth_from_south + PI, 1e-12);
                // Azimuth and altitude are consistent with hour angle and
                // declination (Meeus, ch. 13, and the inverse formulae).
                let (az, alt) = (p.azimuth_from_south, p.altitude);
                let lat = obs.lat;
                let ha = az.sin().atan2(az.cos() * lat.sin() + alt.tan() * lat.cos());
                common::assert_close_rad(ha, p.hour_angle, 1e-9);
                let dec = (lat.sin() * alt.sin() - lat.cos() * alt.cos() * az.cos()).asin();
                assert_close(dec, p.dec, 1e-9);
            }
        }
    }
}

#[test]
fn sun_at_transit_is_due_south_and_rises_east_at_equinox() {
    let obs = london();
    let jd = jd0(2024, time::Month::Mar, 20);
    let r = observer::rise_transit_set(Body::Sun, &obs, jd, Horizon::Standard);
    let noon = observer::topocentric_position(Body::Sun, &obs, r.transit.unwrap());
    assert_close_deg(noon.azimuth_from_north.to_degrees(), 180.0, 0.01);
    assert_close_deg(noon.azimuth_from_south.to_degrees(), 0.0, 0.01);
    // Altitude at noon: 90° − φ + δ, with δ ≈ 0 at the equinox.
    assert_close(noon.altitude.to_degrees(), 90.0 - 51.5074, 0.2);

    let equator = Observer::new(0.0, 0.0, 0.0);
    let r = observer::rise_transit_set(Body::Sun, &equator, jd, Horizon::Standard);
    let dawn = observer::topocentric_position(Body::Sun, &equator, r.rise.unwrap());
    let dusk = observer::topocentric_position(Body::Sun, &equator, r.set.unwrap());
    assert_close_deg(dawn.azimuth_from_north.to_degrees(), 90.0, 0.5);
    assert_close_deg(dusk.azimuth_from_north.to_degrees(), 270.0, 0.5);
    // At sunrise the upper limb appears on the horizon: the centre's
    // apparent altitude is about minus the semidiameter.
    assert_close(
        dawn.apparent_altitude.to_degrees(),
        -16.0 / 60.0,
        3.0 / 60.0,
    );
}

// Meeus, ch. 40: near the horizon the parallax lowers the Moon by
// asin(ρ sin π cos h) and shifts its azimuth by less than π/300.
#[test]
fn moon_parallax_near_the_horizon() {
    let obs = Observer::from_degrees_east(40.0, -74.0, 0.0);
    let jd = 2460310.6;
    let topo = observer::topocentric_position(Body::Moon, &obs, jd);
    let jde = jd + delta_t(jd) / 86400.0;
    let mut geo = observer::geocentric_apparent_position(Body::Moon, jde);
    let parallax = geo.eq_hz_parallax;
    geo.eq_hz_parallax = 0.0;
    let g = observer::topocentric_from_geocentric(&geo, &obs, gast(jd, jde));
    assert!(g.altitude.abs() < deg(5.0)); // this instant is near moonrise

    let rho = astro::planet::earth::rho(obs.lat);
    let expected = (rho * parallax.sin() * g.altitude.cos()).asin();
    assert_close(g.altitude - topo.altitude, expected, deg(arcsec(5.0)));
    assert!((g.azimuth_from_south - topo.azimuth_from_south).abs() < parallax / 300.0 * 1.05);
    // The observer is closer to the Moon than the centre of the Earth is.
    assert!(topo.distance < geo.distance);
}

#[test]
fn refraction() {
    // At the standard altitude -0°34' the refraction lifts the body to the
    // horizon.
    let r = observer::refracted_altitude(deg(-34.0 / 60.0), 1010.0, 10.0);
    assert_close(r.to_degrees(), 0.0, 1.0 / 60.0);
    // Eq. 16.4 at 45°: 1.02 / tan(45.2055°) ≈ 1.013'.
    let r = observer::refracted_altitude(deg(45.0), 1010.0, 10.0);
    assert_close((r.to_degrees() - 45.0) * 60.0, 1.0127, 1e-3);
    // Half the pressure, half the refraction; colder air refracts more.
    let half = observer::refracted_altitude(deg(10.0), 505.0, 10.0) - deg(10.0);
    let full = observer::refracted_altitude(deg(10.0), 1010.0, 10.0) - deg(10.0);
    assert_close(half, full / 2.0, 1e-12);
    assert!(observer::refracted_altitude(deg(10.0), 1010.0, -20.0) - deg(10.0) > full);
    // No refraction at the zenith or far below the horizon.
    assert_close(
        observer::refracted_altitude(FRAC_PI_2, 1010.0, 10.0),
        FRAC_PI_2,
        1e-12,
    );
    assert_eq!(
        observer::refracted_altitude(deg(-5.0), 1010.0, 10.0),
        deg(-5.0)
    );
    // Monotonic between -1° and 90°.
    let mut last = f64::NEG_INFINITY;
    for i in -10..=900 {
        let a = observer::refracted_altitude(deg(i as f64 / 10.0), 1010.0, 10.0);
        assert!(a > last);
        last = a;
    }
}

#[test]
fn observer_constructors() {
    let s = sydney();
    assert_close(s.lat.to_degrees(), -33.8688, 1e-12);
    assert_close(s.long.to_degrees(), -151.2093, 1e-12); // west positive
    let o = Observer::new(0.5, -1.0, 12.0);
    assert_eq!((o.lat, o.long, o.height), (0.5, -1.0, 12.0));
    let gp = o.geograph_point();
    assert_eq!((gp.lat, gp.long), (0.5, -1.0));
}

// Bad input must not panic.
#[test]
fn non_finite_input() {
    let obs = Observer::new(f64::NAN, 0.0, 0.0);
    let r = observer::rise_transit_set(Body::Sun, &obs, 2460400.5, Horizon::Standard);
    assert_eq!((r.rise, r.set), (None, None));
    let r = observer::rise_transit_set(Body::Moon, &london(), f64::NAN, Horizon::Standard);
    assert_eq!((r.rise, r.set), (None, None));
    let p = observer::topocentric_position(Body::Venus, &london(), f64::INFINITY);
    assert!(p.altitude.is_nan() || p.altitude.is_finite());
}

// Consecutive windows neither lose nor duplicate an event at the midnight
// boundary: the sunsets of an observer whose sunset drifts across 0h UT
// are all found, one day apart.
#[test]
fn events_across_window_boundaries() {
    // At 40° N, 67° W the June sunsets happen around 0h UT.
    let obs = Observer::from_degrees_east(40.0, -67.0, 0.0);
    let start = jd0(2024, time::Month::May, 25);
    let days = 40;
    let mut sets = Vec::new();
    let mut empty_days = 0;
    for d in 0..days {
        let r = observer::rise_transit_set(Body::Sun, &obs, start + d as f64, Horizon::Standard);
        match r.set {
            Some(t) => sets.push(t),
            None => empty_days += 1,
        }
    }
    assert!(empty_days >= 1, "the sunset never crossed 0h UT");
    assert!(sets.len() >= days - 2);
    for w in sets.windows(2) {
        assert_close(w[1] - w[0], 1.0, 2.0 * MINUTE);
    }
    // Some sunsets fall just before and some just after 0h UT.
    let frac = |t: f64| (t - 0.5).rem_euclid(1.0);
    assert!(sets.iter().any(|t| frac(*t) > 0.97) && sets.iter().any(|t| frac(*t) < 0.03));

    // The same for the Moon, whose events come about 50 minutes later each
    // day: the gaps between consecutive moonrises are a day plus the daily
    // retardation, never twice that, never less than a day.
    let boston = Observer::from_degrees_east(42.3601, -71.0589, 0.0);
    let start = jd0(2024, time::Month::Mar, 1);
    let rises: Vec<f64> = (0..60)
        .filter_map(|d| {
            observer::rise_transit_set(Body::Moon, &boston, start + d as f64, Horizon::Standard)
                .rise
        })
        .collect();
    assert!(rises.len() >= 56);
    for w in rises.windows(2) {
        let gap = w[1] - w[0];
        assert!(
            (1.0 + 5.0 * MINUTE..1.0 + 100.0 * MINUTE).contains(&gap),
            "{gap}"
        );
    }
}

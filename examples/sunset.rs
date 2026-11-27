//! Sunrise, sunset and twilight for a place and a date.
//!
//! Answers GitHub issue #3 ("How to calculate sunset?"): everything that
//! `transit::time` needs is computed by `observer::rise_transit_set` from
//! just a place and a date.
//!
//! ```sh
//! cargo run --example sunset -- [lat_deg] [east_long_deg] [year] [month] [day] [utc_offset_h]
//! ```
//!
//! Latitude is positive north, longitude positive **east** (as on maps).
//! Without arguments: Greenwich, 2024 June 21, times in UTC+1 (BST).

use astro::observer::{self, Body, Horizon, HorizonStatus, Observer, RiseTransitSet};
use astro::time;

fn main() {
    let args: Vec<f64> = std::env::args()
        .skip(1)
        .map(|a| a.parse().expect("arguments must be numbers"))
        .collect();
    let arg = |i: usize, default: f64| args.get(i).copied().unwrap_or(default);

    let (lat, east_long) = (arg(0, 51.4769), arg(1, -0.0005));
    let (year, month, day) = (arg(2, 2024.0) as i16, arg(3, 6.0) as u8, arg(4, 21.0));
    let utc_offset = arg(5, 1.0);

    let place = Observer::from_degrees_east(lat, east_long, 0.0);

    // Julian Day of 0h UT on the date, then of local midnight: the search
    // window is the local calendar day.
    let jd_0h_ut = time::julian_day(&time::Date {
        year,
        month: month_from_number(month),
        decimal_day: day,
        cal_type: time::CalType::Gregorian,
    });
    let local_midnight = jd_0h_ut - utc_offset / 24.0;
    let clock = |t: f64| hh_mm_ss((t - local_midnight) * 24.0);

    println!(
        "{:.4}° {}, {:.4}° {} on {year}-{month:02}-{day:02}, times in UTC{utc_offset:+}",
        lat.abs(),
        if lat >= 0.0 { "N" } else { "S" },
        east_long.abs(),
        if east_long >= 0.0 { "E" } else { "W" },
    );

    let sun = |horizon| observer::rise_transit_set(Body::Sun, &place, local_midnight, horizon);
    let day = sun(Horizon::Standard);
    let civil = sun(Horizon::CivilTwilight);
    let nautical = sun(Horizon::NauticalTwilight);
    let astronomical = sun(Horizon::AstronomicalTwilight);

    println!();
    print_pair(
        "astronomical dawn",
        "astronomical dusk",
        &astronomical,
        &clock,
    );
    print_pair("nautical dawn", "nautical dusk", &nautical, &clock);
    print_pair("civil dawn", "civil dusk", &civil, &clock);
    print_pair("sunrise", "sunset", &day, &clock);
    if let Some(t) = day.transit {
        let noon = observer::topocentric_position(Body::Sun, &place, t);
        println!(
            "{:<18} {}   (altitude {:.2}°)",
            "solar noon",
            clock(t),
            noon.apparent_altitude.to_degrees()
        );
    }
    match (day.rise, day.set) {
        (Some(r), Some(s)) if s > r => {
            println!("{:<18} {}", "day length", hh_mm_ss((s - r) * 24.0))
        }
        _ => {}
    }

    println!();
    let moon = observer::rise_transit_set(Body::Moon, &place, local_midnight, Horizon::Standard);
    print_pair("moonrise", "moonset", &moon, &clock);
}

/// Prints a pair of events, or why they do not happen.
fn print_pair(rise: &str, set: &str, events: &RiseTransitSet, clock: &dyn Fn(f64) -> String) {
    match events.status {
        HorizonStatus::AlwaysAbove => println!("{rise} / {set}: none, above all day"),
        HorizonStatus::AlwaysBelow => println!("{rise} / {set}: none, below all day"),
        HorizonStatus::RisesOrSets => {
            let show = |t: Option<f64>| t.map(clock).unwrap_or_else(|| "--:--:--".into());
            println!("{rise:<18} {}", show(events.rise));
            println!("{set:<18} {}", show(events.set));
        }
    }
}

fn hh_mm_ss(hours: f64) -> String {
    let s = (hours * 3600.0).round() as i64;
    format!("{:02}:{:02}:{:02}", s / 3600, s / 60 % 60, s % 60)
}

fn month_from_number(m: u8) -> time::Month {
    use time::Month::*;
    [
        Jan, Feb, Mar, Apr, May, June, July, Aug, Sept, Oct, Nov, Dec,
    ][(m.clamp(1, 12) - 1) as usize]
}

//! Azimuth and elevation of the Sun, the Moon, the planets and a few bright
//! stars, seen from a place on the Earth at an instant.
//!
//! Answers GitHub issues #15 and #16 ("compute topocentric azimuth/elevation
//! of a planet from a place on Earth") with `observer::topocentric_position`.
//!
//! ```sh
//! cargo run --example planet_altaz -- [lat_deg] [east_long_deg] [height_m] \
//!     [year] [month] [day] [hour] [minute] [second]
//! ```
//!
//! Latitude is positive north, longitude positive **east**, the time is UT.
//! Without arguments: Palomar Observatory, 2003 August 28, 3h17m UT (Meeus,
//! example 40.a).

use astro::observer::{self, Body, Observer};
use astro::time;

fn main() {
    let args: Vec<f64> = std::env::args()
        .skip(1)
        .map(|a| a.parse().expect("arguments must be numbers"))
        .collect();
    let arg = |i: usize, default: f64| args.get(i).copied().unwrap_or(default);

    let place = Observer::from_degrees_east(arg(0, 33.356111), arg(1, -116.863056), arg(2, 1706.0));
    let (year, month, day) = (arg(3, 2003.0) as i16, arg(4, 8.0) as u8, arg(5, 28.0));
    let (hour, minute, second) = (arg(6, 3.0), arg(7, 17.0), arg(8, 0.0));

    let jd_ut = time::julian_day(&time::Date {
        year,
        month: month_from_number(month),
        decimal_day: day + (hour + minute / 60.0 + second / 3600.0) / 24.0,
        cal_type: time::CalType::Gregorian,
    });

    let (lat, long) = (place.lat.to_degrees(), place.long.to_degrees());
    println!(
        "{:.4}° {}, {:.4}° {}, {} m, {year}-{month:02}-{day:02} \
         {hour:02}:{minute:02}:{second:02} UT (JD {jd_ut:.5})\n",
        lat.abs(),
        if lat >= 0.0 { "N" } else { "S" },
        long.abs(),
        if long >= 0.0 { "W" } else { "E" }, // the crate's longitudes are west-positive
        place.height
    );
    println!(
        "{:<9} {:>9} {:>9} {:>9} {:>9}  {:>12} {:>11} {:>12}",
        "body", "az (N→E)", "az (S→W)", "alt", "app. alt", "RA (topo)", "Dec (topo)", "dist (AU)"
    );

    // J2000.0 catalogue positions (degrees).
    let star = |ra: f64, dec: f64| Body::Star {
        asc: ra.to_radians(),
        dec: dec.to_radians(),
    };
    let bodies = [
        ("Sun", Body::Sun),
        ("Moon", Body::Moon),
        ("Mercury", Body::Mercury),
        ("Venus", Body::Venus),
        ("Mars", Body::Mars),
        ("Jupiter", Body::Jupiter),
        ("Saturn", Body::Saturn),
        ("Uranus", Body::Uranus),
        ("Neptune", Body::Neptune),
        ("Sirius", star(101.287155, -16.716116)),
        ("Vega", star(279.234735, 38.783689)),
        ("Polaris", star(37.954561, 89.264109)),
    ];

    for (name, body) in bodies {
        let p = observer::topocentric_position(body, &place, jd_ut);
        println!(
            "{name:<9} {:>8.3}° {:>8.3}° {:>8.3}° {:>8.3}°  {} {} {:>12.6}",
            p.azimuth_from_north.to_degrees(),
            p.azimuth_from_south.to_degrees(),
            p.altitude.to_degrees(),
            p.apparent_altitude.to_degrees(),
            hms(p.asc.to_degrees()),
            dms(p.dec.to_degrees()),
            p.distance
        );
    }
    println!(
        "\naz (N→E): compass azimuth, from north through east; az (S→W): Meeus' \
         azimuth, from south through west.\nalt: airless altitude of the centre; \
         app. alt: with refraction for 1010 mbar and 10 °C."
    );
}

fn hms(deg: f64) -> String {
    let s = (deg.rem_euclid(360.0) / 15.0 * 36000.0).round() as i64; // tenths of a second
    format!(
        "{:02}h{:02}m{:04.1}s",
        s / 36000,
        s / 600 % 60,
        (s % 600) as f64 / 10.0
    )
}

fn dms(deg: f64) -> String {
    let sign = if deg < 0.0 { '-' } else { '+' };
    let s = (deg.abs() * 3600.0).round() as i64;
    format!("{sign}{:02}°{:02}′{:02}″", s / 3600, s / 60 % 60, s % 60)
}

fn month_from_number(m: u8) -> time::Month {
    use time::Month::*;
    [
        Jan, Feb, Mar, Apr, May, June, July, Aug, Sept, Oct, Nov, Dec,
    ][(m.clamp(1, 12) - 1) as usize]
}

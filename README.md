# astro-rust

[![crates.io](https://img.shields.io/crates/v/astro.svg)](https://crates.io/crates/astro)
[![docs.rs](https://docs.rs/astro/badge.svg)](https://docs.rs/astro)
[![CI](https://github.com/saurvs/astro-rust/actions/workflows/ci.yml/badge.svg)](https://github.com/saurvs/astro-rust/actions/workflows/ci.yml)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)

**[▶ Live demo](https://saurvs.github.io/astro-rust-demo/)** — this crate compiled to WebAssembly,
computing the sky, the solar system and planetary paths in your browser
([source](https://github.com/saurvs/astro-rust-demo)).

[![The astro-rust live demo: a London evening sky, the rose of Venus and the solar system](https://raw.githubusercontent.com/saurvs/astro-rust-demo/main/media/demo.gif)](https://saurvs.github.io/astro-rust-demo/)

`astro` is a library of astronomical algorithms for Rust: where the Sun,
Moon, planets and their satellites are; when things rise, set, transit,
change phase or eclipse; and the calendars, time scales and coordinate
transformations needed to get there.

It implements the algorithms of Jean Meeus' *Astronomical Algorithms*
(2nd edition, Willmann-Bell, 1998) — every chapter, checked against the
book's worked examples — and, where the book recommends it, the full
theories behind them (the complete VSOP87 series for the planets, for
example). It is pure Rust with no dependencies.

* [Live demo](https://saurvs.github.io/astro-rust-demo/)
* [API documentation](https://docs.rs/astro)
* [Changelog](CHANGELOG.md)

## What's inside

| Area | Modules |
|------|---------|
| Time and calendars: Julian Day, ΔT, sidereal time, Easter, Jewish and Islamic calendars | [`time`], [`calendar`] |
| Observing from a place on Earth: altitude/azimuth, rising/transit/setting, twilight | [`observer`], [`transit`], [`atmos`], [`parallax`] |
| Coordinates: equatorial, ecliptic, horizontal, galactic; precession, nutation, aberration, apparent places of stars | [`coords`], [`precess`], [`nutation`], [`aberr`], [`ecliptic`], [`star`] |
| The Sun: position, equinoxes and solstices, equation of time, physical ephemeris | [`sun`], [`season`] |
| The planets: positions (VSOP87), orbital elements, phenomena (oppositions, elongations, stations), perihelia, magnitudes, physical ephemerides of Mars and Jupiter, Saturn's rings, satellites of Jupiter and Saturn, Pluto | [`planet`], [`pluto`] |
| Orbits: Kepler's equation, elliptic, parabolic and near-parabolic motion, node passages | [`orbit`], [`asteroid`] |
| The Moon: position, phases, perigee and apogee, nodes, maximum declinations, illumination, librations | [`lunar`] |
| Eclipses, stellar magnitudes, binary stars, sundials | [`eclipse`], [`star`], [`binary_star`], [`sundial`] |
| Numerical tools: interpolation, curve fitting, root finding, angular separations, conjunctions | [`interpol`], [`fitting`], [`iteration`], [`angle`], [`conjunction`], [`geometry`] |

## Usage

```toml
[dependencies]
astro = "3"
```

### Sunrise and sunset

```rust
use astro::observer::{self, Body, Horizon, Observer};

// Greenwich, latitude 51.48° N, longitude 0.00° (east-positive, as on maps)
let greenwich = Observer::from_degrees_east(51.4769, -0.0005, 0.0);

// The window to search: one day starting at 2024 June 21, 0h UT
let jd = 2460482.5;
let day = observer::rise_transit_set(Body::Sun, &greenwich, jd, Horizon::Standard);

let ut_hours = |t: f64| (t - jd) * 24.0;
let rise = ut_hours(day.rise.unwrap()); // about 3.73 h UT (03:44)
let set = ut_hours(day.set.unwrap()); // about 20.36 h UT (20:21)
assert!((rise - 3.73).abs() < 0.05 && (set - 20.36).abs() < 0.05);
```

Twilight works the same way with `Horizon::CivilTwilight`,
`Horizon::NauticalTwilight` or `Horizon::AstronomicalTwilight`. Bodies that
don't rise or set that day (the midnight Sun, circumpolar stars) are
reported as such instead of producing nonsense. The runnable
[`examples/sunset.rs`](examples/sunset.rs) prints a full table:

```sh
cargo run --example sunset -- 40.7128 -74.0060 2024 12 21 -5
```

### Where is a planet in the sky?

```rust
use astro::observer::{self, Body, Observer};

let place = Observer::from_degrees_east(51.4769, -0.0005, 0.0);
let jd_ut = 2460482.5; // 2024 June 21, 0h UT

let jupiter = observer::topocentric_position(Body::Jupiter, &place, jd_ut);
let azimuth = jupiter.azimuth_from_north.to_degrees(); // compass bearing
let altitude = jupiter.apparent_altitude.to_degrees(); // includes refraction
assert!((0.0..360.0).contains(&azimuth) && altitude.abs() <= 90.0);
```

See [`examples/planet_altaz.rs`](examples/planet_altaz.rs) for every planet,
the Sun and the Moon.

### Positions of the Sun, Moon and planets

```rust
use astro::{lunar, planet, sun};

// Apparent geocentric positions; JDE is a Julian Ephemeris Day (TT)
let sun = sun::apparent_position(2448908.5); // 1992 October 13.0 TD
assert!((sun.asc.to_degrees() - 198.378).abs() < 0.01);

let moon = lunar::apparent_position(2448724.5); // 1992 April 12.0 TD
assert!((moon.asc.to_degrees() - 134.688).abs() < 0.001);
assert!((moon.dist - 368_409.7).abs() < 1.0); // kilometres

let venus = planet::geocent_apprnt_pos(&planet::Planet::Venus, 2448976.5).unwrap();
assert!((venus.apprnt_dec.to_degrees() + 18.888).abs() < 0.001);
```

### Time

```rust
use astro::time::{self, CalType, Date, Month};

// Sputnik 1: 1957 October 4.81 UT
let date = Date { year: 1957, month: Month::Oct, decimal_day: 4.81, cal_type: CalType::Gregorian };
let jd = time::julian_day(&date);
assert!((jd - 2436116.31).abs() < 1e-6);

// Terrestrial Time runs ahead of UT by ΔT (about 32 s in 1957)
let jde = time::julian_ephemeris_day(jd, time::delta_t(1957, 10));

assert!((jde - jd) * 86400.0 > 25.0);

// Greenwich mean sidereal time, in radians in [0, 2π)
let gmst = time::mn_sidr(jd);
assert!((0.0..std::f64::consts::TAU).contains(&gmst));
```

### Events

```rust
use astro::{calendar, eclipse, lunar, season, time};

// The next solar eclipse after 2024 January 1: the total eclipse of April 8
let e = eclipse::next_solar_eclipse(2460310.5).unwrap();
assert_eq!(e.kind, eclipse::SolarEclipseKind::Total);
assert!((e.jde - 2460409.26).abs() < 0.01);

// The March equinox of 2025 (March 20, about 09:01 UT)
let equinox = season::instant(2025, season::Event::MarchEquinox);
assert!((equinox - 2460754.877).abs() < 0.01);

// The full Moon nearest to 2024 January 25 (17:54 UT)
let full = lunar::phase::nearest_phase(2460335.0, lunar::Phase::Full);
assert!((full.jde - 2460335.246).abs() < 0.01);

// Easter Sunday 2025 (Gregorian calendar): April 20
let easter = calendar::easter_gregorian(2025).unwrap();
assert_eq!((easter.month, easter.decimal_day), (time::Month::Apr, 20.0));
```

## Conventions

* **Angles are in radians** throughout the API (use `f64::to_degrees` and the
  helpers in [`angle`] for display).
* **Geographic longitude is positive west** of Greenwich, as in Meeus.
  `Observer::from_degrees_east` and the helpers in [`coords`] convert from the
  east-positive longitudes used on maps.
* **Azimuth is measured from the south**, westward, as in Meeus; the
  [`observer`] API also gives the compass azimuth (from north, eastward).
* Times are Julian Days. Functions say whether they expect **JD** (Universal
  Time) or **JDE** (Terrestrial/Dynamical Time); convert with
  `time::julian_ephemeris_day` and `time::delta_t`.
* Distances are in AU unless documented otherwise (the Moon: kilometres).

## Accuracy

The algorithms are those of the book, so the accuracy is the book's: for
example a few arcseconds for the Moon, better than an arcsecond for the
planets (full VSOP87), and minutes of time for phenomena computed from
series (eclipses, lunar phases). Each function documents its accuracy and
valid time span where it matters. ΔT beyond the present is a prediction;
pass an observed value when you need sub-second timing.

## Contributing

Bug reports and pull requests are welcome. See
[CONTRIBUTING.md](CONTRIBUTING.md) for conventions and the checks CI runs.

## Thanks

Thanks to everyone who has reported issues, sent fixes and suggested
improvements:
Oliver Schneider (@oli-obk), Liam Dawson (@liamdawson), manguluka
(@bingX4devX), Manuel Hässig (@mhaessig), Ashley Blewer (@ablwr), Thomas
Daniels (@thomas-daniels), Christopher H. Jordan (@cjordan), Laurent Valdes
(@valdo404), @andresv, @emosenkis, @Razican, @meshula, @clangdo,
@dean-shaff and @smroid.

## References

* Jean Meeus, *Astronomical Algorithms*, 2nd ed., Willmann-Bell, 1998 — the
  source of nearly every algorithm here.
* P. Bretagnon and G. Francou, *Planetary theories in rectangular and
  spherical variables: VSOP87 solutions*, Astron. Astrophys. 202 (1988) —
  planetary positions ([data](https://ftp.imcce.fr/pub/ephem/planets/vsop87/)).
* M. Chapront-Touzé and J. Chapront, ELP-2000/82 — the lunar theory behind
  the book's series for the Moon.
* F. Espenak and J. Meeus, *Five Millennium Canon of Solar Eclipses*, NASA
  (2006) — the polynomial expressions for ΔT
  ([online](https://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html)).
* IERS Conventions (2010) — for the more recent IAU 2000/2006
  precession–nutation models, which are not implemented here
  ([online](https://iers-conventions.obspm.fr/)).

<!-- Module links: rustdoc and GitHub both resolve these to docs.rs. -->
[`aberr`]: https://docs.rs/astro/latest/astro/aberr/
[`angle`]: https://docs.rs/astro/latest/astro/angle/
[`asteroid`]: https://docs.rs/astro/latest/astro/asteroid/
[`atmos`]: https://docs.rs/astro/latest/astro/atmos/
[`binary_star`]: https://docs.rs/astro/latest/astro/binary_star/
[`calendar`]: https://docs.rs/astro/latest/astro/calendar/
[`conjunction`]: https://docs.rs/astro/latest/astro/conjunction/
[`coords`]: https://docs.rs/astro/latest/astro/coords/
[`eclipse`]: https://docs.rs/astro/latest/astro/eclipse/
[`ecliptic`]: https://docs.rs/astro/latest/astro/ecliptic/
[`fitting`]: https://docs.rs/astro/latest/astro/fitting/
[`geometry`]: https://docs.rs/astro/latest/astro/geometry/
[`interpol`]: https://docs.rs/astro/latest/astro/interpol/
[`iteration`]: https://docs.rs/astro/latest/astro/iteration/
[`lunar`]: https://docs.rs/astro/latest/astro/lunar/
[`nutation`]: https://docs.rs/astro/latest/astro/nutation/
[`observer`]: https://docs.rs/astro/latest/astro/observer/
[`orbit`]: https://docs.rs/astro/latest/astro/orbit/
[`parallax`]: https://docs.rs/astro/latest/astro/parallax/
[`planet`]: https://docs.rs/astro/latest/astro/planet/
[`pluto`]: https://docs.rs/astro/latest/astro/pluto/
[`precess`]: https://docs.rs/astro/latest/astro/precess/
[`season`]: https://docs.rs/astro/latest/astro/season/
[`star`]: https://docs.rs/astro/latest/astro/star/
[`sun`]: https://docs.rs/astro/latest/astro/sun/
[`sundial`]: https://docs.rs/astro/latest/astro/sundial/
[`time`]: https://docs.rs/astro/latest/astro/time/
[`transit`]: https://docs.rs/astro/latest/astro/transit/

## License

MIT — see [LICENSE.md](LICENSE.md).

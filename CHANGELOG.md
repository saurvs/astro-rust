# Changelog

## 3.0.0 — unreleased

A complete overhaul. Every chapter of Meeus' *Astronomical Algorithms*
(2nd ed.) was re-checked against the code, the crate was moved to modern
Rust, and the algorithms that were missing were implemented. Many results
change because bugs were fixed — read **Behaviour changes** before
upgrading.

### Rust and tooling

- Rust 2024 edition; minimum supported Rust version 1.85.
- No dependencies (runtime or dev).
- rustfmt-formatted; clippy and rustdoc clean with `-D warnings`.
- GitHub Actions CI (stable, beta, MSRV; Linux, macOS, Windows) replaces
  Travis.
- Exported macros use `$crate::` paths, so they work from other crates.
- The VSOP87 tables are `static` data instead of being rebuilt on every
  call: clean builds are several times faster.
- About 1000 tests (every worked example of the book, with intermediate
  values, plus invariant and regression tests) and 90+ doctests.

### New modules and APIs

- `observer` — high-level helpers: topocentric azimuth/altitude of the Sun,
  Moon, planets and stars for a place on Earth, and rise/transit/set times
  with civil, nautical and astronomical twilight, robust at polar
  latitudes (#3, #15, #16). See `examples/sunset.rs` and
  `examples/planet_altaz.rs`.
- `calendar` — Gregorian and Julian Easter, the Jewish calendar (Pesach,
  New Year, year types), Islamic ↔ Julian/Gregorian conversion (ch. 8–9).
- `season` — equinoxes and solstices, and the lengths of the seasons
  (ch. 27).
- `eclipse` — solar and lunar eclipses: type, time of maximum, γ,
  magnitudes, contact times; next/previous eclipse search (ch. 54).
- `conjunction`, `geometry` — conjunctions and least separations from
  tabulated positions; bodies in a straight line; smallest circle
  containing three bodies (ch. 18–20).
- `fitting`, `iteration` — least-squares fitting and robust root finding
  (ch. 4–5); `interpol` gains extremum/zero finding, interpolation to
  halves and Lagrange's formula (ch. 3).
- `sundial` — planar sundials of any orientation (ch. 58).
- `planet::phenomena` — conjunctions, oppositions, greatest elongations and
  stations of the planets (ch. 36).
- `planet::apsis` — perihelion and aphelion of the planets, by formula and
  by VSOP87 search (ch. 38).
- `planet::magnitude` — phase angle helpers, Venus' approximate
  illuminated fraction (ch. 41).
- `lunar::{position, illumination, libration, phase, apsis, node,
  declination}` — the Moon's position, illuminated fraction, librations and
  physical ephemeris, phases, perigee/apogee, node passages and maximum
  declinations, with iterators over date ranges (ch. 47–53).
- Struct-returning, single-argument (JDE) entry points for the physical
  ephemerides of the Sun, Moon, Mars, Jupiter and Saturn's ring, the
  apparent place of a star, and the apparent position of the Sun, Moon and
  planets.
- High-accuracy (E5) positions of the Galilean satellites, satellite
  shadows and phenomena (ch. 44).
- Many smaller additions in `time`, `angle`, `coords`, `precess`,
  `nutation`, `aberr`, `sun`, `orbit`, `parallax`, `pluto`, `consts`; see
  the API docs.

### Bug fixes

Pull requests merged: #9 (README function names, by Ashley Blewer), #12
(dropped term in quarter-phase times, by Thomas Daniels), #20 (IAU
obliquity in the wrong time unit, by Laurent Valdes).

Issues fixed: #19 (`coords::dec_frm_hz`, reported by smroid), #13
(horizontal → equatorial conversion, reported by Christopher H. Jordan, who
also proposed the combined function now available as
`coords::hr_angl_and_dec_frm_hz`).

Other fixes, by area:

- **Time** — ΔT polynomials had sign errors in five segments (e.g. 1799
  gave 184 s instead of 13.8 s); `weekday_frm_date` ignored the calendar and
  panicked for negative JD; `decimal_year` was one day late; overflow and
  NaN handling in `julian_day`/`date_frm_julian_day`; apparent sidereal time
  was not normalised and its macro evaluated its argument three times.
- **Angles** — `anglr_sepr` returned NaN or 0 for tiny separations;
  `limit_to_360`/`limit_to_two_PI` could return 360°/2π or leave huge
  inputs unreduced; sign handling of `deg_frm_dms`/`deg_frm_hms` when the
  leading field is zero; `dms_frm_deg` returned 59.999…″.
- **Coordinates** — galactic longitude not reduced to [0, 360°);
  `geograph_geocent_lat_diff` returned degrees; `geodesic_dist` returned NaN
  for identical points; low-altitude refraction formulas used the wrong
  angle (wrong sign of the result); several angle range fixes.
- **Precession, nutation, aberration** — wrong signs of the T-terms of z and
  θ in `precess_eq_coords`; `precess_eq_coords_FK5` mixed degrees and
  radians (and is Newcomb's precession — now `precess_eq_coords_newcomb`);
  `star::eq_coords_frm_motion` multiplied the RA proper motion by 15;
  quadrant error in `star::angl_between_north_celes_and_eclip_pole`;
  `precess_orb_elements` turned retrograde orbits prograde.
- **Sun** — southern solar latitudes came out near 2π; the equation of time
  was not wrapped (+1433 min instead of −7 min) and had a wrong
  coefficient.
- **Planets** — heliocentric latitudes were returned in [0, 2π) instead of
  [−π/2, π/2]; wrong coefficients in `orb_elements` (Mercury's semimajor
  axis 0.0387 AU, T² terms multiplied by T³, three wrong signs);
  `geocent_eq_coords` was missing sines; light-time iterated only once;
  `phase_angl` returned cos i; planetary magnitudes used radians where the
  formulas need degrees (Venus −4.9 instead of −3.8); Saturn's ΔU in
  radians and not wrapped.
- **Physical ephemerides** — Jupiter's light-time applied twice and a wrong
  aberration term; Mars' Q was not computed; several angles (central
  meridians, P) not reduced.
- **Orbits** — Kepler's equation solver could be inaccurate or loop forever
  for e near 1; parabolic motion lost precision before perihelion;
  near-parabolic motion reported false non-convergence and could hang.
- **Pluto, parallax** — Pluto's longitude not reduced; topocentric ecliptic
  latitude wrong for λ in 90°–270° and division by zero; topocentric RA not
  normalised and declination outside ±90°.
- **Moon** — `time_of_phase` returned a phase of the wrong lunation for most
  dates; missing E factors in the node-passage and quarter-phase terms;
  `total_libr` passed arguments in the wrong order; the sign of a
  topocentric libration correction; `illum_frac_frm_eq_coords` returned
  NaN; several angles not reduced.
- **Satellites** — Saturn's satellites used a fixed light-time for every
  date and an inconsistent B1950 epoch.
- **Stars, eclipses, semidiameters** — `star::combined_mag` was missing a
  logarithm; wrong term in the apparent eccentricity of binary orbits;
  asteroid diameter used the wrong logarithm base.

### Behaviour changes

- Angles that are longitudes or right ascensions are returned in [0, 2π);
  latitudes and declinations in [−π/2, π/2]; hour angles in (−π, π].
- `sun::semidiameter` returns radians (was degrees);
  `lunar::semidiameter` returns the angle (was its sine);
  `asteroid::apparent_diameter` returns radians (was milli-arcseconds).
- `lunar::time_of_phase` returns the phase nearest to the given date.
- `lunar::optical_libr` returns signed librations.
- Results that used to be wrong (see above) now differ.

### Deprecated

- `transit::time` → `transit::day_fraction` (returns `Result`, handles
  circumpolar bodies)
- `precess::precess_eq_coords_FK5` → `precess::precess_eq_coords_newcomb`
- `star::proper_motion_in_eq_coords` (misnamed; it is ecliptical) →
  `precess::proper_motion_in_ecl_coords`
- `pluto::semdiameter` → `pluto::semidiameter`
- `parallax::topopcent_ecl_coords` → `parallax::topocent_ecl_coords`
- `planet::saturn::ring::out_edge_inner_ing` →
  `planet::saturn::ring::out_edge_inner_ring`

### Not included

- TLE/SGP4 satellite tracking (#1): outside the scope of the book; see the
  `sgp4` crate.
- Using the external `vsop87` crate (#4): the crate already embeds the
  complete VSOP87D series it needs, with no dependency.
- IAU 2000/2006 precession–nutation: not part of the book.
- Appendix IV (polynomial fits for the giant planets, 1998–2025): its
  validity has expired and the full VSOP87 covers it.

## 2.0.0 and earlier

See the git history.

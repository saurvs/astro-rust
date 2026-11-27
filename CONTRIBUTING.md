# Contributing to astro-rust

Thanks for helping. This file describes the conventions the crate follows so
that new algorithms fit in with the existing ones.

## Checks

Every change must pass, with no warnings:

```sh
cargo fmt --all --check
cargo clippy --all-targets
cargo test --all-targets
cargo test --doc
cargo doc --no-deps
```

CI runs the same commands with `-D warnings` on stable, beta and the minimum
supported Rust version (see `rust-version` in `Cargo.toml`).

## Units and conventions

- Angles in the public API are in **radians**, unless a function name or its
  documentation says otherwise. Angles in degrees are only used for
  human-facing helpers such as `angle::deg_frm_dms`.
- Times are Julian days (`f64`). Say in the docs whether a function expects
  a Julian Day (UT) or a Julian Ephemeris Day (TT/TD).
- Distances are in AU unless stated otherwise; distances on the Earth are in
  metres.
- **Geographic longitude is measured positive westward** from Greenwich and
  negative eastward, as in Meeus. Every function that takes a longitude must
  say so in its docs.
- **Azimuth is measured westward from the south**, as in Meeus. Use the
  helpers in `coords` to convert to the navigator's convention (from north,
  eastward).
- Results that are angles are normalised to a stated range (for example
  `[0, 2π)` for longitudes, `[-π/2, π/2]` for latitudes).
- Prefer returning small structs with named fields over long tuples for new
  APIs. Derive `Debug, Clone, Copy, PartialEq` on plain data types.
- Do not panic on bad input in library code. Return `Option` or `Result`
  where an input can be out of range, and say what happens in the docs.
- No new dependencies (runtime or dev) without a strong reason.

## Naming

Existing names (`mn_oblq_IAU`, `geocent_ecl_pos`, ...) are kept for
compatibility. New functions should use descriptive snake_case names. If a
function must change signature, keep the old one as a `#[deprecated]` wrapper
for one major version.

## Documentation

- Every public item has a doc comment: what it computes, every argument with
  its unit and convention, what it returns, and the valid time span or
  accuracy where that matters.
- Cite the source algorithm, e.g. "Meeus, *Astronomical Algorithms*, 2nd ed.,
  chapter 22". Write the explanation in your own words. Do not copy passages
  from books into the code.
- Add a short runnable example (doctest) to the main entry point of each
  module.

## Tests

- Every public function has at least one test.
- Use the worked examples of the source (Meeus' examples) as reference values,
  with a tolerance that matches the precision of the published result.
- Add invariant tests too: round trips (A→B→A), symmetries, known limits,
  angle wrap-around, and dates across calendar and century boundaries.
- Integration tests live in `tests/<module>.rs`. Shared helpers are in
  `tests/common/mod.rs` (`mod common;` then `common::assert_close(..)`).
  Unit tests for private helpers go in a `#[cfg(test)] mod tests` block.
- Each test should name the example it reproduces, e.g.
  `// Meeus, example 22.a`.

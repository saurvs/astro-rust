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

//! Planets in perihelion and aphelion (Meeus, *Astronomical Algorithms*,
//! 2nd ed., chapter 38).
//!
//! Two methods are provided.
//!
//! * [`time_of_passage`] evaluates Meeus' quadratic expressions in an
//!   integer cycle number. They describe unperturbed Keplerian orbits, so
//!   they are only approximate: a few hours off for Mars, up to about half a
//!   month for Jupiter, more than a month for Saturn and years for Uranus
//!   and Neptune. For the Earth they give the passages of the Earth–Moon
//!   barycentre; [`earth_time_of_passage`] adds Meeus' corrections for the
//!   Earth itself (mean error about 3 hours).
//! * [`accurate_passage`] finds the extremum of the heliocentric distance
//!   computed with the full VSOP87 theory of
//!   [`planet::heliocent_coords`](super::heliocent_coords), searching around
//!   the approximate instant, and also returns the distance itself. This is
//!   how the accurate dates quoted by Meeus were obtained.
//!
//! Because of the motion of the Sun about the barycentre of the Solar System,
//! the distance of a slow planet such as Neptune can reach two minima (or two
//! maxima) a few years apart. [`accurate_passage`] returns the deepest one
//! near the expected date, the true perihelion (or aphelion);
//! [`rad_vec_extrema`] lists every local extremum in an interval, including
//! the secondary ones: the shallow maximum between two minima near
//! perihelion (Meeus' *apheloid*), the shallow minimum between two maxima
//! near aphelion (his *perihelioid*), and the second minimum or maximum.
//!
//! All instants are Julian Ephemeris Days (dynamical time, TD).
//!
//! # Numbering of the passages
//!
//! Passages are numbered by an integer `n`. Perihelion `n = 0` is the first
//! after the beginning of 2000 (`n < 0` before it), and aphelion `n` is the
//! one that follows perihelion `n`. Meeus' `k` is `n` for a perihelion and
//! `n + 0.5` for an aphelion: his `k = +16.5` is aphelion `n = 16`, and his
//! `k = −119.5` is aphelion `n = −120`.
//!
//! # Example
//!
//! Venus' perihelion nearest to 1978 October 15 (Meeus, example 38.a):
//!
//! ```
//! use astro::planet::Planet;
//! use astro::planet::apsis::{self, Apsis};
//!
//! let n = apsis::nearest_n(&Planet::Venus, Apsis::Perihelion, 1978.79).unwrap();
//! assert_eq!(n, -35);
//!
//! let jde = apsis::time_of_passage(&Planet::Venus, Apsis::Perihelion, n);
//! assert!((jde - 2443873.704).abs() < 0.001); // 1978 December 31, 5h TD
//!
//! // The same passage from the VSOP87 radius vector
//! let passage = apsis::accurate_passage(&Planet::Venus, Apsis::Perihelion, n);
//! assert!((passage.jde - jde).abs() < 1.0);
//! assert!((passage.rad_vec - 0.7184).abs() < 0.0001);
//! ```

use crate::planet::{self, Planet};

/// Perihelion or aphelion
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Apsis {
    /// Least distance from the Sun
    Perihelion,
    /// Greatest distance from the Sun
    Aphelion,
}

/// A passage through perihelion or aphelion, from VSOP87
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Passage {
    /// Instant of the passage *| Julian Ephemeris Day*
    pub jde: f64,
    /// Heliocentric distance of the planet at that instant *| in AU*
    pub rad_vec: f64,
}

/// Whether a local extremum of the heliocentric distance is a minimum or
/// a maximum
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Extremum {
    /// Local minimum of the distance to the Sun
    Minimum,
    /// Local maximum of the distance to the Sun
    Maximum,
}

/// A local extremum of a planet's heliocentric distance
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RadVecExtremum {
    /// Minimum or maximum
    pub kind: Extremum,
    /// Instant of the extremum *| Julian Ephemeris Day*
    pub jde: f64,
    /// Heliocentric distance at that instant *| in AU*
    pub rad_vec: f64,
}

/// Computes Meeus' approximate, real-valued `k` for a date
///
/// The integer nearest to the result numbers the nearest perihelion; an
/// integer plus 0.5 numbers an aphelion (see the [module
/// documentation](self)). These are Meeus' rounded linear expressions;
/// [`nearest_n`] inverts the passage formulas instead, which stays exact
/// farther from 2000.
///
/// # Arguments
///
/// * `planet`: Any variant of [`Planet`]
/// * `year`  : Year with decimals, e.g. `1978.79` for 1978 October 15, in
///   astronomical numbering (`0` is 1 B.C.)
pub fn approx_k(planet: &Planet, year: f64) -> f64 {
    let (rate, epoch) = match planet {
        Planet::Mercury => (4.15201, 2000.12),
        Planet::Venus => (1.62549, 2000.53),
        Planet::Earth => (0.99997, 2000.01),
        Planet::Mars => (0.53166, 2001.78),
        Planet::Jupiter => (0.08430, 2011.20),
        Planet::Saturn => (0.03393, 2003.52),
        Planet::Uranus => (0.01190, 2051.1),
        Planet::Neptune => (0.00607, 2047.5),
    };
    rate * (year - epoch)
}

/// Computes the number `n` of the perihelion or aphelion nearest to a date
///
/// # Returns
///
/// `Some(n)` such that [`time_of_passage`] is the passage closest in time to
/// `year` (like the formulas themselves, this is meaningful for some
/// thousands of years around 2000). `None` if `year` is not finite, or so
/// far from 2000 (beyond about 10¹⁵ orbital periods) that `n` would not be
/// exact.
///
/// # Arguments
///
/// * `planet`: Any variant of [`Planet`]
/// * `apsis` : Perihelion or aphelion
/// * `year`  : Year with decimals, e.g. `2032.0` for the beginning of 2032,
///   in astronomical numbering (`0` is 1 B.C.). It is read as a Julian
///   epoch, `JDE = 2451545.0 + 365.25·(year − 2000)`.
pub fn nearest_n(planet: &Planet, apsis: Apsis, year: f64) -> Option<i64> {
    let target = 2451545.0 + (year - 2000.0) * 365.25;
    // Invert JDE = A + B·k + C·k² for k by fixed-point iteration (the k²
    // term is a small correction wherever the formula means anything).
    // This is closer than Meeus' rounded approx_k, which drifts by whole
    // periods some tens of thousands of years from 2000.
    let (a, b, c) = kepler_coeffs(planet);
    let mut k = (target - a) / b;
    for _ in 0..3 {
        k = (target - a - c * k * k) / b;
    }
    let guess = (k - half(apsis)).round();
    // Within this bound n ± 1 is exact in f64 (2⁵³).
    if !guess.is_finite() || guess.abs() >= 9_007_199_254_740_992.0 {
        return None;
    }
    let guess = guess as i64;

    let mut best = guess;
    let mut best_dist = f64::INFINITY;
    for n in guess.saturating_sub(1)..=guess.saturating_add(1) {
        let dist = (time_of_passage(planet, apsis, n) - target).abs();
        if dist < best_dist {
            best = n;
            best_dist = dist;
        }
    }
    Some(best)
}

/// Computes the instant of a passage through perihelion or aphelion from
/// Meeus' Keplerian expressions
///
/// The expressions ignore the mutual perturbations of the planets, so the
/// result is only approximate (see the [module documentation](self)); use
/// [`accurate_passage`] for a precise instant. For [`Planet::Earth`] the
/// result refers to the Earth–Moon barycentre; see
/// [`earth_time_of_passage`].
///
/// # Returns
///
/// * `jde`: Instant of the passage *| Julian Ephemeris Day*
///
/// # Arguments
///
/// * `planet`: Any variant of [`Planet`]
/// * `apsis` : Perihelion or aphelion
/// * `n`     : Number of the passage; perihelion `0` is the first after the
///   beginning of 2000, and aphelion `n` follows perihelion `n`
///
/// # Example
///
/// ```
/// use astro::planet::Planet;
/// use astro::planet::apsis::{self, Apsis};
///
/// // First aphelion of Mars in 2032 (Meeus, example 38.b, k = +16.5)
/// let n = apsis::nearest_n(&Planet::Mars, Apsis::Aphelion, 2032.0).unwrap();
/// assert_eq!(n, 16);
/// let jde = apsis::time_of_passage(&Planet::Mars, Apsis::Aphelion, n);
/// assert!((jde - 2463530.456).abs() < 0.001); // 2032 October 24, 23h TD
/// ```
pub fn time_of_passage(planet: &Planet, apsis: Apsis, n: i64) -> f64 {
    let k = n as f64 + half(apsis);
    let (a, b, c) = kepler_coeffs(planet);
    a + k * (b + k * c)
}

/// Computes the correction from the Earth–Moon barycentre to the Earth in
/// the instant of a passage of the Earth through perihelion or aphelion
///
/// Because of the Moon, the centre of the Earth reaches its least or
/// greatest distance from the Sun up to about a day and a half before or
/// after the Earth–Moon barycentre does. This is Meeus' five-term
/// correction (at most ±1.54 days); the corrected instants have a mean
/// error of about 3 hours (6 hours at most).
///
/// # Returns
///
/// * `correction`: Days to add to [`time_of_passage`] for
///   [`Planet::Earth`]
///
/// # Arguments
///
/// * `apsis`: Perihelion or aphelion
/// * `n`    : Number of the passage, as in [`time_of_passage`]
pub fn earth_correction(apsis: Apsis, n: i64) -> f64 {
    let k = n as f64 + half(apsis);
    let angles = [
        (328.41, 132.788585),
        (316.13, 584.903153),
        (346.20, 450.380738),
        (136.95, 659.306737),
        (249.52, 329.653368),
    ];
    let coeffs = match apsis {
        Apsis::Perihelion => [1.278, -0.055, -0.091, -0.056, -0.045],
        Apsis::Aphelion => [-1.352, 0.061, 0.062, 0.029, 0.031],
    };
    angles
        .iter()
        .zip(coeffs)
        .map(|(&(a0, a1), c)| {
            // Reduce k·a₁ first: |k| can be in the thousands.
            let angle = (a0 + (k * a1).rem_euclid(360.0)).to_radians();
            c * angle.sin()
        })
        .sum()
}

/// Computes the instant of a passage of the Earth (not the Earth–Moon
/// barycentre) through perihelion or aphelion
///
/// This is [`time_of_passage`] for [`Planet::Earth`] plus
/// [`earth_correction`].
///
/// # Returns
///
/// * `jde`: Instant of the passage *| Julian Ephemeris Day*
///
/// # Arguments
///
/// * `apsis`: Perihelion or aphelion
/// * `n`    : Number of the passage, as in [`time_of_passage`]
///
/// # Example
///
/// ```
/// use astro::planet::Planet;
/// use astro::planet::apsis::{self, Apsis};
///
/// // Meeus, chapter 38: perihelion k = -10, in 1990 January
/// let barycentre = apsis::time_of_passage(&Planet::Earth, Apsis::Perihelion, -10);
/// assert!((barycentre - 2447894.911).abs() < 0.001);
/// let earth = apsis::earth_time_of_passage(Apsis::Perihelion, -10);
/// assert!((earth - 2447896.172).abs() < 0.001); // 1990 January 4, 16h TD
/// ```
pub fn earth_time_of_passage(apsis: Apsis, n: i64) -> f64 {
    time_of_passage(&Planet::Earth, apsis, n) + earth_correction(apsis, n)
}

/// Computes a passage through perihelion or aphelion from the VSOP87
/// heliocentric distance
///
/// The distance given by
/// [`planet::heliocent_coords`](super::heliocent_coords) is sampled over a
/// fifth of an orbital period centred on the approximate instant (from
/// [`earth_time_of_passage`] for the Earth and [`time_of_passage`]
/// otherwise), and the deepest minimum (perihelion) or highest maximum
/// (aphelion) is then located precisely by golden-section search and
/// parabolic interpolation.
///
/// The search converges to well under a minute, so the result is as
/// accurate as VSOP87 itself, which is valid for several thousand years
/// around J2000.0. Near an apsis the distance changes very slowly, though:
/// for Uranus and Neptune a tiny error in the distance moves the instant by
/// a lot, and such instants are best quoted to the day, as Meeus does. For
/// the Earth the result refers to the centre of the Earth, not to the
/// Earth–Moon barycentre.
///
/// Each call evaluates VSOP87 a few hundred times.
///
/// # Arguments
///
/// * `planet`: Any variant of [`Planet`]
/// * `apsis` : Perihelion or aphelion
/// * `n`     : Number of the passage, as in [`time_of_passage`]
///
/// # Example
///
/// ```
/// use astro::planet::Planet;
/// use astro::planet::apsis::{self, Apsis};
///
/// // Meeus, chapter 38: the Keplerian formula puts Saturn's perihelion
/// // n = -2 on 1944 July 30, but it actually occurred on 1944 September 8.
/// let passage = apsis::accurate_passage(&Planet::Saturn, Apsis::Perihelion, -2);
/// assert!((passage.jde - 2431341.5).abs() < 1.0);
/// assert!((passage.rad_vec - 9.0288).abs() < 0.0001);
/// ```
pub fn accurate_passage(planet: &Planet, apsis: Apsis, n: i64) -> Passage {
    let guess = match planet {
        Planet::Earth => earth_time_of_passage(apsis, n),
        _ => time_of_passage(planet, apsis, n),
    };
    search_passage(planet, apsis, guess)
}

/// Deepest extremum of the distance within a fifth of a period of `guess`
/// (the window follows the samples if the extremum lies at its edge)
fn search_passage(planet: &Planet, apsis: Apsis, guess: f64) -> Passage {
    let period = kepler_coeffs(planet).1;
    let half_window = 0.1 * period;
    let samples = 200;
    let step = 2.0 * half_window / samples as f64;
    let sign = match apsis {
        Apsis::Perihelion => 1.0,
        Apsis::Aphelion => -1.0,
    };
    // Minimise sign·R.
    let f = |jde: f64| sign * rad_vec(planet, jde);

    // The extremum should lie well inside the window; if the lowest sample
    // is at an edge, move the window there and look again.
    let mut centre = guess;
    let mut best_t = guess;
    for _ in 0..4 {
        let start = centre - half_window;
        let mut best_i = 0;
        let mut best_f = f64::INFINITY;
        for i in 0..=samples {
            let v = f(start + i as f64 * step);
            if v < best_f {
                best_i = i;
                best_f = v;
            }
        }
        best_t = start + best_i as f64 * step;
        if best_i != 0 && best_i != samples {
            break;
        }
        centre = best_t;
    }

    let jde = refine_minimum(&f, best_t, step);
    Passage {
        jde,
        rad_vec: rad_vec(planet, jde),
    }
}

/// Finds every local minimum and maximum of a planet's heliocentric
/// distance in an interval of time, from VSOP87
///
/// The distance is sampled every thousandth of an orbital period (about two
/// months for Neptune, two hours for Mercury) and each extremum is then
/// located precisely. Extrema closer together than a few sampling steps, or
/// within one step of the ends of the interval, may be missed. This makes
/// visible the double minima and maxima of Neptune's distance described by
/// Meeus. The cost is about a thousand VSOP87 evaluations per orbital
/// period covered.
///
/// # Returns
///
/// The extrema in chronological order. The list is empty if
/// `jde_end <= jde_start`, if either bound is not finite, or if the
/// interval spans more than ten thousand orbital periods (ten million
/// samples; about 2400 years for Mercury, 1.6 million for Neptune).
///
/// # Arguments
///
/// * `planet`   : Any variant of [`Planet`]
/// * `jde_start`: Start of the interval *| Julian Ephemeris Day*
/// * `jde_end`  : End of the interval *| Julian Ephemeris Day*
///
/// # Example
///
/// ```
/// use astro::planet::Planet;
/// use astro::planet::apsis::{self, Extremum};
///
/// // Neptune's distance had two minima around its perihelion of 1876,
/// // with a maximum between them (Meeus, chapter 38).
/// let extrema = apsis::rad_vec_extrema(&Planet::Neptune, 2405000.0, 2411000.0);
/// let kinds: Vec<Extremum> = extrema.iter().map(|e| e.kind).collect();
/// assert_eq!(kinds, [Extremum::Minimum, Extremum::Maximum, Extremum::Minimum]);
/// assert!((extrema[0].rad_vec - 29.8148).abs() < 0.0001);
/// ```
pub fn rad_vec_extrema(planet: &Planet, jde_start: f64, jde_end: f64) -> Vec<RadVecExtremum> {
    let mut extrema = Vec::new();
    if !(jde_start.is_finite() && jde_end.is_finite()) || jde_end <= jde_start {
        return extrema;
    }
    let period = kepler_coeffs(planet).1;
    let steps = ((jde_end - jde_start) / (period / 1000.0)).ceil().max(2.0);
    // Guard against absurdly long intervals.
    if steps > 1.0e7 {
        return extrema;
    }
    let steps = steps as usize;
    let step = (jde_end - jde_start) / steps as f64;
    let r = |jde: f64| rad_vec(planet, jde);

    let mut prev = r(jde_start);
    let mut curr = r(jde_start + step);
    for i in 1..steps {
        let t = jde_start + i as f64 * step;
        let next = r(t + step);
        let kind = if curr < prev && curr <= next {
            Some(Extremum::Minimum)
        } else if curr > prev && curr >= next {
            Some(Extremum::Maximum)
        } else {
            None
        };
        if let Some(kind) = kind {
            let jde = match kind {
                Extremum::Minimum => refine_minimum(&r, t, step),
                Extremum::Maximum => refine_minimum(&|x| -r(x), t, step),
            };
            if jde >= jde_start && jde <= jde_end {
                extrema.push(RadVecExtremum {
                    kind,
                    jde,
                    rad_vec: r(jde),
                });
            }
        }
        prev = curr;
        curr = next;
    }
    extrema
}

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

/// 0 for a perihelion, 0.5 for an aphelion
#[inline]
fn half(apsis: Apsis) -> f64 {
    match apsis {
        Apsis::Perihelion => 0.0,
        Apsis::Aphelion => 0.5,
    }
}

/// `(A, B, C)` of `JDE = A + B·k + C·k²`; `B` is the anomalistic period in
/// days
fn kepler_coeffs(planet: &Planet) -> (f64, f64, f64) {
    match planet {
        Planet::Mercury => (2451590.257, 87.96934963, 0.0),
        Planet::Venus => (2451738.233, 224.7008188, -0.0000000327),
        Planet::Earth => (2451547.507, 365.2596358, 0.0000000156),
        Planet::Mars => (2452195.026, 686.9957857, -0.0000001187),
        Planet::Jupiter => (2455636.936, 4332.897065, 0.0001367),
        Planet::Saturn => (2452830.12, 10764.21676, 0.000827),
        Planet::Uranus => (2470213.5, 30694.8767, -0.00541),
        Planet::Neptune => (2468895.1, 60190.33, 0.03429),
    }
}

#[inline]
fn rad_vec(planet: &Planet, jde: f64) -> f64 {
    planet::heliocent_coords(planet, jde).2
}

/// Locates the minimum of `f` in `[t - step, t + step]`, where `t` is a
/// sample at least as low as its neighbours `t ± step`.
///
/// Golden-section search shrinks the bracket to a millionth of a day, then
/// a parabola through the last three points gives the final abscissa.
fn refine_minimum(f: &dyn Fn(f64) -> f64, t: f64, step: f64) -> f64 {
    const INV_PHI: f64 = 0.618_033_988_749_894_9;
    const TOL: f64 = 1.0e-6;

    let (mut a, mut b) = (t - step, t + step);
    let mut c = b - INV_PHI * (b - a);
    let mut d = a + INV_PHI * (b - a);
    let (mut fc, mut fd) = (f(c), f(d));
    for _ in 0..200 {
        if b - a <= TOL {
            break;
        }
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - INV_PHI * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + INV_PHI * (b - a);
            fd = f(d);
        }
    }

    // Vertex of the parabola through three equally spaced points around the
    // bracket's centre; keep the centre if the fit is degenerate.
    let m = 0.5 * (a + b);
    let h = (b - a).max(TOL);
    let (y0, y1, y2) = (f(m - h), f(m), f(m + h));
    let denom = y0 - 2.0 * y1 + y2;
    if denom > 0.0 {
        let offset = 0.5 * h * (y0 - y2) / denom;
        if offset.abs() <= h {
            return m + offset;
        }
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn refine_minimum_of_parabola() {
        let f = |x: f64| (x - 1.234_567).powi(2) + 3.0;
        let x = refine_minimum(&f, 1.2, 0.1);
        assert!((x - 1.234_567).abs() < 1e-7, "{x}");
    }

    // A guess so far off that the extremum lies outside the first window
    #[test]
    fn search_follows_the_samples() {
        let planet = Planet::Mars;
        let exact = accurate_passage(&planet, Apsis::Perihelion, 3);
        let period = kepler_coeffs(&planet).1;
        for offset in [-0.25, 0.18, 0.3] {
            let found = search_passage(&planet, Apsis::Perihelion, exact.jde + offset * period);
            assert!((found.jde - exact.jde).abs() < 1e-4, "{offset}");
        }
    }

    #[test]
    fn half_values() {
        assert_eq!(half(Apsis::Perihelion), 0.0);
        assert_eq!(half(Apsis::Aphelion), 0.5);
    }
}

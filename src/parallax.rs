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

//! Correction for parallax: topocentric coordinates
//!
//! Converts geocentric positions (as seen from the centre of the Earth) into
//! topocentric ones (as seen by an observer on the Earth's surface), in
//! equatorial and in ecliptic coordinates, following Meeus, *Astronomical
//! Algorithms*, 2nd ed., chapter 40.
//!
//! Conventions:
//!
//! * Geographic **longitude is positive west** of Greenwich and negative
//!   east, as in Meeus.
//! * Hour angles are measured westward from the observer's meridian:
//!   `H = θ0 − L − α` with `θ0` the sidereal time at Greenwich and `L` the
//!   (west-positive) longitude.
//! * The observer's position enters through `ρ sin φ'` and `ρ cos φ'`
//!   (geocentric radius in Earth equatorial radii times the sine and cosine
//!   of the geocentric latitude, Meeus chapter 11); see
//!   [`planet::earth::rho_sin_cos_phi`].
//!
//! # Example
//!
//! Mars from Palomar Observatory on 2003 August 28 at 3h17m UT (Meeus,
//! example 40.a):
//!
//! ```
//! use astro::parallax;
//!
//! let asc = 339.530208_f64.to_radians();
//! let dec = (-15.771083_f64).to_radians();
//! let parallax = parallax::eq_hz_parallax(0.37276); // 23.592"
//! // H = θ0 − L − α, with L = +7h47m27s (west) and θ0 = 1h40m45s
//! let hr_angl = 288.7958_f64.to_radians();
//!
//! let topo = parallax::topocent_eq(asc, dec, hr_angl, parallax, 0.546861, 0.836339);
//! assert!((topo.del_asc.to_degrees() * 3600.0 / 15.0 - 1.29).abs() < 0.005); // +1.29 s
//! assert!((topo.dec.to_degrees() - -(15.0 + 46.0 / 60.0 + 30.0 / 3600.0)).abs() < 0.05 / 3600.0);
//! ```

use crate::angle;
use crate::coords;
use crate::planet;
use std::f64::consts::TAU;

/// Computes the equatorial horizontal parallax of a celestial body
///
/// Uses `sin π = sin 8.794″ / Δ` (Meeus, eq. 40.1), where 8.794″ is the
/// solar parallax, i.e. the equatorial horizontal parallax at 1 AU.
///
/// # Returns
///
/// * `eq_hz_parllx`: Equatorial horizontal parallax of the celestial
///   body *| in radians*. NaN if `dist_to_earth` is smaller than one Earth
///   radius (about 4.26e-5 AU).
///
/// # Arguments
///
/// * `dist_to_earth`: The celestial body's distance to the Earth *| in AU*
#[inline]
pub fn eq_hz_parallax(dist_to_earth: f64) -> f64 {
    (angle::deg_frm_dms(0, 0, 8.794).to_radians().sin() / dist_to_earth).asin()
}

/// Topocentric equatorial position of a celestial body
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TopocentEq {
    /// Topocentric right ascension `α'` *| in radians*, in `[0, 2π)`
    pub asc: f64,
    /// Topocentric declination `δ'` *| in radians*, in `[-π/2, π/2]`
    pub dec: f64,
    /// Topocentric hour angle `H'`, measured westward from the meridian
    /// *| in radians*, in `[0, 2π)`
    pub hr_angl: f64,
    /// Parallax in right ascension `Δα = α' − α` *| in radians*, in
    /// `[−π, π]`
    pub del_asc: f64,
    /// Parallax in declination `Δδ = δ' − δ` *| in radians*
    pub del_dec: f64,
}

/// Computes the topocentric equatorial coordinates of a celestial body
/// (rigorous method)
///
/// Implements Meeus' rigorous formulas (eqs. 40.2 and 40.3, equivalently
/// 40.6 and 40.7). The declination is evaluated in the form
/// `δ' = atan2(C, √(A² + B²))` of eq. 40.6, which is identical to eq. 40.3
/// but stays in the correct quadrant even for a very close body near a
/// celestial pole, where the denominator of eq. 40.3 changes sign. Valid
/// for any body, including the Moon.
///
/// # Returns
///
/// The topocentric position and the parallax corrections, see
/// [`TopocentEq`].
///
/// # Arguments
///
/// * `asc`         : Geocentric right ascension `α` *| in radians*
/// * `dec`         : Geocentric declination `δ` *| in radians*
/// * `hr_angl`     : Geocentric hour angle `H`, measured westward from the
///   meridian *| in radians* (`H = θ0 − L − α`, longitude `L` positive west)
/// * `eq_hz_parllx`: Equatorial horizontal parallax `π` *| in radians*
///   (see [`eq_hz_parallax`])
/// * `rho_sin_phi` : `ρ sin φ'` of the observer (Earth equatorial radius = 1)
/// * `rho_cos_phi` : `ρ cos φ'` of the observer (Earth equatorial radius = 1)
pub fn topocent_eq(
    asc: f64,
    dec: f64,
    hr_angl: f64,
    eq_hz_parllx: f64,
    rho_sin_phi: f64,
    rho_cos_phi: f64,
) -> TopocentEq {
    let sin_parllx = eq_hz_parllx.sin();
    let (sin_h, cos_h) = hr_angl.sin_cos();
    let (sin_dec, cos_dec) = dec.sin_cos();

    // Eq. 40.2: the parallax in right ascension.
    let del_asc =
        (-rho_cos_phi * sin_parllx * sin_h).atan2(cos_dec - rho_cos_phi * sin_parllx * cos_h);

    // Eq. 40.6: components of the topocentric direction.
    let a = cos_dec * sin_h;
    let b = cos_dec * cos_h - rho_cos_phi * sin_parllx;
    let c = sin_dec - rho_sin_phi * sin_parllx;
    let dec_1 = c.atan2(a.hypot(b));

    TopocentEq {
        asc: norm_2pi(asc + del_asc),
        dec: dec_1,
        hr_angl: norm_2pi(hr_angl - del_asc),
        del_asc,
        del_dec: dec_1 - dec,
    }
}

/// Computes the topocentric equatorial coordinates of a celestial body
/// (non-rigorous, low-accuracy method)
///
/// Uses the first-order expressions of Meeus, eqs. 40.4 and 40.5:
/// `Δα = −π ρ cos φ' sin H / cos δ` and
/// `Δδ = −π (ρ sin φ' cos δ − ρ cos φ' cos H sin δ)`.
/// They are adequate for the Sun, the planets and comets, but not for the
/// Moon, and they break down near the celestial poles (`cos δ → 0`, where
/// `Δα` becomes infinite); use [`topocent_eq`] there.
///
/// For a planet or a comet the same expressions, with the signs of `Δα` and
/// `Δδ` reversed, reduce an observed topocentric position to a geocentric
/// one: call this function with the topocentric `α'`, `δ'`, `H'` and
/// subtract `del_asc` and `del_dec` instead of adding them.
///
/// # Returns
///
/// The topocentric position and the parallax corrections, see
/// [`TopocentEq`].
///
/// # Arguments
///
/// Same as [`topocent_eq`].
pub fn topocent_eq_approx(
    asc: f64,
    dec: f64,
    hr_angl: f64,
    eq_hz_parllx: f64,
    rho_sin_phi: f64,
    rho_cos_phi: f64,
) -> TopocentEq {
    let (sin_h, cos_h) = hr_angl.sin_cos();
    let (sin_dec, cos_dec) = dec.sin_cos();

    let del_asc = -eq_hz_parllx * rho_cos_phi * sin_h / cos_dec;
    let del_dec = -eq_hz_parllx * (rho_sin_phi * cos_dec - rho_cos_phi * cos_h * sin_dec);

    TopocentEq {
        asc: norm_2pi(asc + del_asc),
        dec: dec + del_dec,
        hr_angl: norm_2pi(hr_angl - del_asc),
        del_asc,
        del_dec,
    }
}

/// Computes the topocentric equatorial coordinates of a celestial body
///
/// Rigorous method, see [`topocent_eq`]; this convenience form computes the
/// hour angle from the Greenwich sidereal time and the observer's
/// longitude, and `ρ sin φ'`, `ρ cos φ'` from the observer's latitude and
/// height with [`planet::earth::rho_sin_cos_phi`].
///
/// # Returns
///
/// * `topocent_eq_point`: Topocentric equatorial point of the
///   celestial body *| in radians*, right ascension in `[0, 2π)`
///
/// # Arguments
///
/// * `eq_point`      : Geocentric equatorial point of the celestial body
///   *| in radians*
/// * `eq_hz_parllx`  : Equatorial horizontal parallax of the celestial
///   body *| in radians*
/// * `geograph_point`: Geographic point of the observer *| in radians*,
///   longitude **positive west** of Greenwich
/// * `observer_ht`   : Height of the observer above sea level
///   *| in meters*
/// * `greenw_sidr`   : Sidereal time at Greenwich *| in radians* (apparent
///   sidereal time if the right ascension is apparent)
pub fn topocent_eq_coords(
    eq_point: &coords::EqPoint,
    eq_hz_parllx: f64,
    geograph_point: &coords::GeographPoint,
    observer_ht: f64,
    greenw_sidr: f64,
) -> coords::EqPoint {
    let (rho_sin, rho_cos) = planet::earth::rho_sin_cos_phi(geograph_point.lat, observer_ht);

    let geocent_hr_angl =
        coords::hr_angl_frm_observer_long(greenw_sidr, geograph_point.long, eq_point.asc);

    let topo = topocent_eq(
        eq_point.asc,
        eq_point.dec,
        geocent_hr_angl,
        eq_hz_parllx,
        rho_sin,
        rho_cos,
    );

    coords::EqPoint {
        asc: topo.asc,
        dec: topo.dec,
    }
}

/// Topocentric ecliptic position and semidiameter of a celestial body
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TopocentEcl {
    /// Topocentric ecliptic longitude `λ'` *| in radians*, in `[0, 2π)`
    pub long: f64,
    /// Topocentric ecliptic latitude `β'` *| in radians*, in `[-π/2, π/2]`
    pub lat: f64,
    /// Topocentric semidiameter `s'` *| in radians*
    pub semidiameter: f64,
}

/// Computes the topocentric ecliptic coordinates and semidiameter of a
/// celestial body
///
/// Rigorous expressions (Meeus, chapter 40, after J. J. von Littrow), valid
/// for the Moon too. They are evaluated in an equivalent form based on the
/// components of the topocentric direction, which gives the latitude in
/// the correct quadrant for any longitude (the printed `tan β'` formula,
/// fed to a two-argument arctangent, puts `β'` in the wrong half-plane when
/// `cos λ' < 0`) and avoids a division by zero at `λ' = ±90°`.
///
/// # Returns
///
/// See [`TopocentEcl`].
///
/// # Arguments
///
/// * `ecl_long`       : Geocentric ecliptic longitude `λ` *| in radians*
/// * `ecl_lat`        : Geocentric ecliptic latitude `β` *| in radians*
/// * `geocent_semidia`: Geocentric semidiameter `s` *| in radians*
///   (0 if not needed)
/// * `eq_hz_parllx`   : Equatorial horizontal parallax `π` *| in radians*
/// * `rho_sin_phi`    : `ρ sin φ'` of the observer (Earth equatorial
///   radius = 1)
/// * `rho_cos_phi`    : `ρ cos φ'` of the observer (Earth equatorial
///   radius = 1)
/// * `loc_sidr`       : Local sidereal time `θ` *| in radians*
/// * `oblq_eclip`     : Obliquity of the ecliptic `ε` *| in radians*
pub fn topocent_ecl(
    ecl_long: f64,
    ecl_lat: f64,
    geocent_semidia: f64,
    eq_hz_parllx: f64,
    rho_sin_phi: f64,
    rho_cos_phi: f64,
    loc_sidr: f64,
    oblq_eclip: f64,
) -> TopocentEcl {
    let sin_parllx = eq_hz_parllx.sin();
    let (sin_sidr, cos_sidr) = loc_sidr.sin_cos();
    let (sin_oblq, cos_oblq) = oblq_eclip.sin_cos();
    let (sin_long, cos_long) = ecl_long.sin_cos();
    let (sin_lat, cos_lat) = ecl_lat.sin_cos();

    // Topocentric direction, in units of the geocentric distance:
    // x = N, y and z are the numerators of Meeus' tan λ' and tan β'
    // (without the factor cos λ').
    let x = cos_long * cos_lat - rho_cos_phi * sin_parllx * cos_sidr;
    let y = sin_long * cos_lat
        - sin_parllx * (rho_sin_phi * sin_oblq + rho_cos_phi * cos_oblq * sin_sidr);
    let z = sin_lat - sin_parllx * (rho_sin_phi * cos_oblq - rho_cos_phi * sin_oblq * sin_sidr);

    let xy = x.hypot(y);
    // Ratio of the topocentric to the geocentric distance; Meeus'
    // cos λ' cos β' / N equals 1 / dist_ratio.
    let dist_ratio = xy.hypot(z);

    TopocentEcl {
        long: norm_2pi(y.atan2(x)),
        lat: z.atan2(xy),
        semidiameter: (geocent_semidia.sin() / dist_ratio).asin(),
    }
}

/// Computes the topocentric ecliptic coordinates of a celestial body
///
/// Convenience form of [`topocent_ecl`] that computes `ρ sin φ'` and
/// `ρ cos φ'` from the observer's latitude and height with
/// [`planet::earth::rho_sin_cos_phi`].
///
/// # Returns
///
/// `(topocent_ecl_point, topocent_semidia)`
///
/// * `topocent_ecl_point`: Topocentric ecliptic point of the
///   celestial body *| in radians*, longitude in `[0, 2π)`
/// * `topocent_semidia`  : Topocentric semidiameter of the
///   celestial body *| in radians*
///
/// # Arguments
///
/// * `ecl_point`     : Geocentric ecliptic point of the celestial body
///   *| in radians*
/// * `eq_hz_parllx`  : Equatorial horizontal parallax of the celestial
///   body *| in radians*
/// * `geograph_point`: Geographic point of the observer *| in radians*;
///   only the latitude is used (the longitude enters through `loc_sidr`)
/// * `observer_ht`   : Height of the observer above sea level *| in meters*
/// * `loc_sidr`      : Local sidereal time *| in radians*
/// * `eclip_oblq`    : Obliquity of the ecliptic *| in radians*
/// * `geocent_semdia`: Geocentric semidiameter of the celestial body
///   *| in radians*
pub fn topocent_ecl_coords(
    ecl_point: &coords::EclPoint,
    eq_hz_parllx: f64,
    geograph_point: &coords::GeographPoint,
    observer_ht: f64,
    loc_sidr: f64,
    eclip_oblq: f64,
    geocent_semdia: f64,
) -> (coords::EclPoint, f64) {
    let (rho_sin, rho_cos) = planet::earth::rho_sin_cos_phi(geograph_point.lat, observer_ht);

    let topo = topocent_ecl(
        ecl_point.long,
        ecl_point.lat,
        geocent_semdia,
        eq_hz_parllx,
        rho_sin,
        rho_cos,
        loc_sidr,
        eclip_oblq,
    );

    (
        coords::EclPoint {
            long: topo.long,
            lat: topo.lat,
        },
        topo.semidiameter,
    )
}

/// Computes the topocentric ecliptic coordinates of a celestial body
///
/// Misspelt former name of [`topocent_ecl_coords`]; same arguments and
/// result.
#[deprecated(
    since = "3.0.0",
    note = "misspelt; use `parallax::topocent_ecl_coords`"
)]
pub fn topopcent_ecl_coords(
    ecl_point: &coords::EclPoint,
    eq_hz_parllx: f64,
    geograph_point: &coords::GeographPoint,
    observer_ht: f64,
    loc_sidr: f64,
    eclip_oblq: f64,
    geocent_semdia: f64,
) -> (coords::EclPoint, f64) {
    topocent_ecl_coords(
        ecl_point,
        eq_hz_parllx,
        geograph_point,
        observer_ht,
        loc_sidr,
        eclip_oblq,
        geocent_semdia,
    )
}

/// Computes the parallax in altitude of a celestial body
///
/// The topocentric (observed) altitude of a body is lower than its
/// geocentric altitude by the parallax `p`, given by
/// `sin p = ρ sin π cos h` (Meeus, chapter 40). With `h` the geocentric
/// altitude the expression is approximate, but good except for the Moon
/// when high accuracy is needed; it is exact if `h` is the topocentric
/// altitude. The parallax in azimuth is always very small and is ignored.
///
/// # Returns
///
/// * `parllx_in_alt`: Parallax in altitude `p` *| in radians* (subtract it
///   from the geocentric altitude)
///
/// # Arguments
///
/// * `alt`         : Altitude `h` of the body *| in radians*
/// * `eq_hz_parllx`: Equatorial horizontal parallax `π` *| in radians*
/// * `rho`         : Observer's distance from the centre of the Earth, in
///   Earth equatorial radii; 1 is often good enough
#[inline]
pub fn parallax_in_alt(alt: f64, eq_hz_parllx: f64, rho: f64) -> f64 {
    (rho * eq_hz_parllx.sin() * alt.cos()).asin()
}

/// Reduces an angle in radians to `[0, 2π)`.
#[inline]
fn norm_2pi(angle: f64) -> f64 {
    let r = angle.rem_euclid(TAU);
    if r >= TAU { 0.0 } else { r }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn norm_2pi_range() {
        assert_eq!(norm_2pi(-1e-300), 0.0);
        assert!((norm_2pi(-PI) - PI).abs() < 1e-15);
        assert!((norm_2pi(3.0 * TAU + 1.0) - 1.0).abs() < 1e-14);
    }
}

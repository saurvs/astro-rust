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

//! The eight major planets: orbital elements, positions and phenomena
//!
//! Positions (Meeus, *Astronomical Algorithms*, 2nd ed., chapters 31–33):
//!
//! - [`heliocent_coords`]: heliocentric ecliptic longitude, latitude and
//!   radius vector from the complete VSOP87 theory (version D, mean
//!   dynamical ecliptic and equinox of the date).
//! - [`geocent_apprnt_pos`]: apparent geocentric position (light-time,
//!   aberration, FK5 frame, nutation), with the intermediate results.
//! - [`mean_orb_elements_of_date`], [`mean_orb_elements_j2000`]: mean
//!   orbital elements of Mercury–Neptune.
//!
//! Times are Julian Ephemeris Days (TT); angles are in radians.
//!
//! # Example
//!
//! ```
//! use astro::planet::{self, Planet};
//!
//! // Apparent position of Venus on 1992 December 20 at 0h TD (Meeus, example 33.a)
//! let pos = planet::geocent_apprnt_pos(&Planet::Venus, 2448976.5).unwrap();
//! let (h, m, s) = astro::angle::hms_frm_deg(pos.apprnt_asc.to_degrees());
//! assert_eq!((h, m), (21, 4));
//! assert!((s - 41.45).abs() < 0.02);
//! ```

// Phase terms in the VSOP87 tables are data, even when they happen to be close to π.
#[allow(
    clippy::approx_constant,
    clippy::excessive_precision,
    clippy::unreadable_literal
)]
mod VSOPD_87;

pub mod apsis;
pub mod earth;
pub mod jupiter;
pub mod magnitude;
pub mod mars;
pub mod phenomena;
pub mod saturn;

use crate::angle;
use crate::coords;
use crate::ecliptic;
use crate::nutation;
use crate::time;

/// One of the eight major planets
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Planet {
    /// Mercury
    Mercury,
    /// Venus
    Venus,
    /// Earth
    Earth,
    /// Mars
    Mars,
    /// Jupiter
    Jupiter,
    /// Saturn
    Saturn,
    /// Uranus
    Uranus,
    /// Neptune
    Neptune,
}

/// Computes the illuminated fraction of a planet's disk from its phase
/// angle
///
/// `k = (1 + cos i) / 2` (Meeus, eq. 41.1). The same expression holds for
/// the Moon.
///
/// # Returns
///
/// * `illum_frac`: Illuminated fraction of the planet's disk, in `[0, 1]`
///
/// # Arguments
///
/// * `i`: Phase angle of the planet (the angle Sun-planet-Earth)
///   *| in radians*
#[inline]
pub fn illum_frac_frm_phase_angl(i: f64) -> f64 {
    (1.0 + i.cos()) / 2.0
}

/// Computes the illuminated fraction of a planet's disk from its distances
/// to the Sun and the Earth
///
/// `k = ((r + Δ)² − R²) / (4 r Δ)` (Meeus, eq. 41.2), which combines eq.
/// 41.1 with the cosine rule for the phase angle.
///
/// # Returns
///
/// * `illum_frac`: Illuminated fraction of the planet's disk
///
/// # Arguments
///
/// * `r`    : Planet-Sun distance *| in AU*
/// * `delta`: Planet-Earth distance *| in AU*
/// * `R`    : Sun-Earth distance *| in AU*
#[inline]
pub fn illum_frac_frm_dist(r: f64, delta: f64, R: f64) -> f64 {
    let x = r + delta;

    (x * x - R * R) / (4.0 * r * delta)
}

/// Computes a planet's phase angle from its distances to the Sun and the
/// Earth
///
/// The phase angle `i` is the angle Sun-planet-Earth, found from the
/// cosine rule `cos i = (r² + Δ² − R²) / (2 r Δ)` (Meeus, chapter 41).
///
/// # Returns
///
/// * `phase_angl`: Phase angle of the planet *| in radians*, in `[0, π]`
///   (the cosine is clamped to `[−1, 1]` against rounding)
///
/// # Arguments
///
/// * `r`    : Planet-Sun distance *| in AU*
/// * `delta`: Planet-Earth distance *| in AU*
/// * `R`    : Sun-Earth distance *| in AU*
#[inline]
pub fn phase_angl(r: f64, delta: f64, R: f64) -> f64 {
    ((r * r + delta * delta - R * R) / (2.0 * r * delta))
        .clamp(-1.0, 1.0)
        .acos()
}

/// Computes the position angle of the bright limb of a planet
///
/// The position angle `χ` of the midpoint of the illuminated limb, counted
/// eastward from the north point of the disk (not from the planet's axis
/// of rotation). Meeus gives the formula for the Moon (eq. 48.5) and notes
/// that it holds for a planet too; use geocentric coordinates of the same
/// date for both bodies.
///
/// # Returns
///
/// * `pos_angl_of_bright_limb`: Position angle of the midpoint of the
///   illuminated limb *| in radians*, in `[0, 2π)`
///
/// # Arguments
///
/// * `sun_eq_point`   : Equatorial point of the Sun *| in radians*
/// * `planet_eq_point`: Equatorial point of the planet *| in radians*
pub fn pos_angle_of_bright_limb(
    sun_eq_point: coords::EqPoint,
    planet_eq_point: coords::EqPoint,
) -> f64 {
    let a = sun_eq_point.dec.cos();
    let n = a * (sun_eq_point.asc - planet_eq_point.asc).sin();
    let d = sun_eq_point.dec.sin() * planet_eq_point.dec.cos()
        - planet_eq_point.dec.sin() * (sun_eq_point.asc - planet_eq_point.asc).cos() * a;

    let chi = n.atan2(d).rem_euclid(std::f64::consts::TAU);
    if chi >= std::f64::consts::TAU {
        0.0
    } else {
        chi
    }
}

/// Computes a planet's equatorial semidiameter
///
/// # Returns
///
/// * `semidiameter`: Equatorial semidiameter *| in radians*
///
/// # Arguments
///
/// * `planet`           : The [Planet](./enum.Planet.html)
/// * `planet_earth_dist`: Planet-Earth distance *| in AU*
pub fn semidiameter<'a>(planet: &Planet, planet_earth_dist: f64) -> Result<f64, &'a str> {
    let s = match *planet {
        Planet::Mercury => angle::deg_frm_dms(0, 0, 3.360).to_radians(),
        Planet::Venus => angle::deg_frm_dms(0, 0, 8.410).to_radians(),
        Planet::Mars => angle::deg_frm_dms(0, 0, 4.680).to_radians(),
        Planet::Uranus => angle::deg_frm_dms(0, 0, 35.02).to_radians(),
        Planet::Neptune => angle::deg_frm_dms(0, 0, 33.50).to_radians(),

        Planet::Jupiter => jupiter::eq_semidiameter(1.0),
        Planet::Saturn => saturn::eq_semidiameter(1.0),

        Planet::Earth => {
            return Err("Planet::Earth was passed to the function planet::semidiameter()");
        }
    };

    Ok(s / planet_earth_dist)
}

// ---------------------------------------------------------------------------
// Mean orbital elements (Meeus, chapter 31)
// ---------------------------------------------------------------------------

/// Mean orbital elements of a planet at a given instant.
///
/// All angles are in radians. The mean longitude, the longitude of the
/// ascending node, the longitude of the perihelion, the mean anomaly and the
/// argument of the perihelion are normalised to `[0, 2π)`. The inclination
/// is *not* normalised: it is a small signed angle (the inclination of the
/// Earth's orbit on the fixed ecliptic of J2000.0 is negative before
/// A.D. 2000, see [`mean_orb_elements_j2000`]).
///
/// The relations `mn_anom = mn_long - long_perih` and
/// `arg_perih = long_perih - long_ascend_node` hold (modulo 2π).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrbElements {
    /// Mean longitude `L` *| in radians*
    pub mn_long: f64,
    /// Semimajor axis `a` *| in AU*
    pub semimaj_axis: f64,
    /// Eccentricity `e`
    pub ecc: f64,
    /// Inclination `i` of the orbit on the ecliptic *| in radians*
    pub inc: f64,
    /// Longitude of the ascending node `Ω` *| in radians*
    pub long_ascend_node: f64,
    /// Longitude of the perihelion `ϖ = Ω + ω` *| in radians*
    pub long_perih: f64,
    /// Mean anomaly `M = L - ϖ` *| in radians*
    pub mn_anom: f64,
    /// Argument of the perihelion `ω = ϖ - Ω` *| in radians*
    pub arg_perih: f64,
}

/// Cubic polynomial `[a0, a1, a2, a3]` in Julian centuries from J2000.0.
type Cubic = [f64; 4];

/// Polynomials for the mean elements of one planet. Angles in degrees.
struct ElementPolys {
    mn_long: Cubic,
    semimaj_axis: Cubic,
    ecc: Cubic,
    inc: Cubic,
    long_ascend_node: Cubic,
    long_perih: Cubic,
}

/// The angular polynomials of Meeus' table 31.B (fixed ecliptic and equinox
/// of J2000.0). The semimajor axis and the eccentricity do not depend on
/// the reference frame and come from table 31.A.
struct AnglePolys {
    mn_long: Cubic,
    inc: Cubic,
    long_ascend_node: Cubic,
    long_perih: Cubic,
}

const NO_TERM: Cubic = [0.0; 4];

/// Meeus, table 31.A: mean elements referred to the mean dynamical ecliptic
/// and equinox of the date (after Simon et al. 1994, adjusted to VSOP87).
/// The Earth's orbit lies in the ecliptic of date, so its inclination is
/// zero and its node is undefined; the node polynomial is left at zero.
fn elements_of_date(planet: &Planet) -> ElementPolys {
    match *planet {
        Planet::Mercury => ElementPolys {
            mn_long: [252.250906, 149474.0722491, 0.00030350, 0.000000018],
            semimaj_axis: [0.387098310, 0.0, 0.0, 0.0],
            ecc: [0.20563175, 0.000020407, -0.0000000283, -0.00000000018],
            inc: [7.004986, 0.0018215, -0.00001810, 0.000000056],
            long_ascend_node: [48.330893, 1.1861883, 0.00017542, 0.000000215],
            long_perih: [77.456119, 1.5564776, 0.00029544, 0.000000009],
        },
        Planet::Venus => ElementPolys {
            mn_long: [181.979801, 58519.2130302, 0.00031014, 0.000000015],
            semimaj_axis: [0.723329820, 0.0, 0.0, 0.0],
            ecc: [0.00677192, -0.000047765, 0.0000000981, 0.00000000046],
            inc: [3.394662, 0.0010037, -0.00000088, -0.000000007],
            long_ascend_node: [76.679920, 0.9011206, 0.00040618, -0.000000093],
            long_perih: [131.563703, 1.4022288, -0.00107618, -0.000005678],
        },
        Planet::Earth => ElementPolys {
            mn_long: [100.466457, 36000.7698278, 0.00030322, 0.000000020],
            semimaj_axis: [1.000001018, 0.0, 0.0, 0.0],
            ecc: [0.01670863, -0.000042037, -0.0000001267, 0.00000000014],
            inc: NO_TERM,
            long_ascend_node: NO_TERM,
            long_perih: [102.937348, 1.7195366, 0.00045688, -0.000000018],
        },
        Planet::Mars => ElementPolys {
            mn_long: [355.433000, 19141.6964471, 0.00031052, 0.000000016],
            semimaj_axis: [1.523679342, 0.0, 0.0, 0.0],
            ecc: [0.09340065, 0.000090484, -0.0000000806, -0.00000000025],
            inc: [1.849726, -0.0006011, 0.00001276, -0.000000007],
            long_ascend_node: [49.558093, 0.7720959, 0.00001557, 0.000002267],
            long_perih: [336.060234, 1.8410449, 0.00013477, 0.000000536],
        },
        Planet::Jupiter => ElementPolys {
            mn_long: [34.351519, 3036.3027748, 0.00022330, 0.000000037],
            semimaj_axis: [5.202603209, 0.0000001913, 0.0, 0.0],
            ecc: [0.04849793, 0.000163225, -0.0000004714, -0.00000000201],
            inc: [1.303267, -0.0054965, 0.00000466, -0.000000002],
            long_ascend_node: [100.464407, 1.0209774, 0.00040315, 0.000000404],
            long_perih: [14.331207, 1.6126352, 0.00103042, -0.000004464],
        },
        Planet::Saturn => ElementPolys {
            mn_long: [50.077444, 1223.5110686, 0.00051908, -0.000000030],
            semimaj_axis: [9.554909192, -0.0000021390, 0.000000004, 0.0],
            ecc: [0.05554814, -0.000346641, -0.0000006436, 0.0000000034],
            inc: [2.488879, -0.0037362, -0.00001519, 0.000000087],
            long_ascend_node: [113.665503, 0.8770880, -0.00012176, -0.000002249],
            long_perih: [93.057237, 1.9637613, 0.00083753, 0.000004928],
        },
        Planet::Uranus => ElementPolys {
            mn_long: [314.055005, 429.8640561, 0.00030390, 0.000000026],
            semimaj_axis: [19.218446062, -0.0000000372, 0.00000000098, 0.0],
            ecc: [0.04638122, -0.000027293, 0.0000000789, 0.00000000024],
            inc: [0.773197, 0.0007744, 0.00003749, -0.000000092],
            long_ascend_node: [74.005957, 0.5211278, 0.00133947, 0.000018484],
            long_perih: [173.005291, 1.4863790, 0.00021406, 0.000000434],
        },
        Planet::Neptune => ElementPolys {
            mn_long: [304.348665, 219.8833092, 0.00030882, 0.000000018],
            semimaj_axis: [30.110386869, -0.0000001663, 0.00000000069, 0.0],
            ecc: [0.00945575, 0.000006033, 0.0, -0.00000000005],
            inc: [1.769953, -0.0093082, -0.00000708, 0.000000027],
            long_ascend_node: [131.784057, 1.1022039, 0.00025952, -0.000000637],
            long_perih: [48.120276, 1.4262957, 0.00038434, 0.000000020],
        },
    }
}

/// Meeus, table 31.B: the angular mean elements referred to the fixed
/// ecliptic and equinox of J2000.0.
fn angles_j2000(planet: &Planet) -> AnglePolys {
    match *planet {
        Planet::Mercury => AnglePolys {
            mn_long: [252.250906, 149472.6746358, -0.00000536, 0.000000002],
            inc: [7.004986, -0.0059516, 0.00000080, 0.000000043],
            long_ascend_node: [48.330893, -0.1254227, -0.00008833, -0.000000200],
            long_perih: [77.456119, 0.1588643, -0.00001342, -0.000000007],
        },
        Planet::Venus => AnglePolys {
            mn_long: [181.979801, 58517.8156760, 0.00000165, -0.000000002],
            inc: [3.394662, -0.0008568, -0.00003244, 0.000000009],
            long_ascend_node: [76.679920, -0.2780134, -0.00014257, -0.000000164],
            long_perih: [131.563703, 0.0048746, -0.00138467, -0.000005695],
        },
        Planet::Earth => AnglePolys {
            mn_long: [100.466457, 35999.3728565, -0.00000568, -0.000000001],
            inc: [0.0, 0.0130548, -0.00000931, -0.000000034],
            long_ascend_node: [174.873176, -0.2410908, 0.00004262, 0.000000001],
            long_perih: [102.937348, 0.3225654, 0.00014799, -0.000000039],
        },
        Planet::Mars => AnglePolys {
            mn_long: [355.433000, 19140.2993039, 0.00000262, -0.000000003],
            inc: [1.849726, -0.0081477, -0.00002255, -0.000000029],
            long_ascend_node: [49.558093, -0.2950250, -0.00064048, -0.000001964],
            long_perih: [336.060234, 0.4439016, -0.00017313, 0.000000518],
        },
        Planet::Jupiter => AnglePolys {
            mn_long: [34.351519, 3034.9056606, -0.00008501, 0.000000016],
            inc: [1.303267, -0.0019877, 0.00003320, 0.000000097],
            long_ascend_node: [100.464407, 0.1767232, 0.00090700, -0.000007272],
            long_perih: [14.331207, 0.2155209, 0.00072211, -0.000004485],
        },
        Planet::Saturn => AnglePolys {
            mn_long: [50.077444, 1222.1138488, 0.00021004, -0.000000046],
            inc: [2.488879, 0.0025514, -0.00004906, 0.000000017],
            long_ascend_node: [113.665503, -0.2566722, -0.00018399, 0.000000480],
            long_perih: [93.057237, 0.5665415, 0.00052850, 0.000004912],
        },
        Planet::Uranus => AnglePolys {
            mn_long: [314.055005, 428.4669983, -0.00000486, 0.000000006],
            inc: [0.773197, -0.0016869, 0.00000349, 0.000000016],
            long_ascend_node: [74.005957, 0.0741431, 0.00040539, 0.000000119],
            long_perih: [173.005291, 0.0893212, -0.00009470, 0.000000414],
        },
        Planet::Neptune => AnglePolys {
            mn_long: [304.348665, 218.4862002, 0.00000059, -0.000000002],
            inc: [1.769953, 0.0002256, 0.00000023, 0.0],
            long_ascend_node: [131.784057, -0.0061651, -0.00000219, -0.000000078],
            long_perih: [48.120276, 0.0291866, 0.00007610, 0.0],
        },
    }
}

#[inline]
fn cubic(c: &Cubic, t: f64) -> f64 {
    c[0] + t * (c[1] + t * (c[2] + t * c[3]))
}

/// Reduces an angle in degrees to `[0, 360)` and converts it to radians.
#[inline]
fn deg_to_rad_0_2pi(deg: f64) -> f64 {
    let r = deg.rem_euclid(360.0).to_radians();
    if r >= std::f64::consts::TAU { 0.0 } else { r }
}

/// Builds the element set from the angular quantities, in degrees.
fn orb_elements_frm_deg(
    mn_long: f64,
    semimaj_axis: f64,
    ecc: f64,
    inc: f64,
    long_ascend_node: f64,
    long_perih: f64,
) -> OrbElements {
    OrbElements {
        mn_long: deg_to_rad_0_2pi(mn_long),
        semimaj_axis,
        ecc,
        inc: inc.to_radians(),
        long_ascend_node: deg_to_rad_0_2pi(long_ascend_node),
        long_perih: deg_to_rad_0_2pi(long_perih),
        mn_anom: deg_to_rad_0_2pi(mn_long - long_perih),
        arg_perih: deg_to_rad_0_2pi(long_perih - long_ascend_node),
    }
}

/// Computes the mean orbital elements of a planet, referred to the mean
/// dynamical ecliptic and equinox of the date.
///
/// Uses the cubic polynomials of Meeus, *Astronomical Algorithms*, 2nd ed.,
/// chapter 31, table 31.A (from Simon et al. 1994, adjusted by Meeus to
/// agree with VSOP87). These are *mean* elements: they describe a smoothly
/// varying reference orbit and leave out the periodic perturbations, so a
/// Keplerian position computed from them is only approximate: over
/// A.D. 1000–3000 it stays within about 0.05° of VSOP87 in longitude for
/// Mercury to Mars, and within about 1.5° for the giant planets (whose
/// mutual perturbations are large). Use [`heliocent_coords`] for accurate
/// positions. The polynomials are intended for a few thousand years
/// around J2000.0.
///
/// For [`Planet::Earth`] the orbit lies in the ecliptic of date, so the
/// inclination is zero and the node is undefined; by convention the
/// returned `long_ascend_node` is 0 and `arg_perih` equals `long_perih`.
///
/// # Arguments
///
/// * `planet`: Any variant of [`Planet`]
/// * `JD`    : Julian Ephemeris Day (TT)
///
/// # Example
///
/// ```
/// use astro::planet::{self, Planet};
///
/// // Meeus, example 31.a: Mercury on 2065 June 24 at 0h TD
/// let el = planet::mean_orb_elements_of_date(&Planet::Mercury, 2475460.5);
/// assert!((el.mn_long.to_degrees() - 203.494701).abs() < 1e-6);
/// assert!((el.semimaj_axis - 0.387098310).abs() < 1e-9);
/// assert!((el.arg_perih.to_degrees() - 29.367732).abs() < 1e-6);
/// ```
pub fn mean_orb_elements_of_date(planet: &Planet, JD: f64) -> OrbElements {
    let p = elements_of_date(planet);
    let t = time::julian_cent(JD);

    orb_elements_frm_deg(
        cubic(&p.mn_long, t),
        cubic(&p.semimaj_axis, t),
        cubic(&p.ecc, t),
        cubic(&p.inc, t),
        cubic(&p.long_ascend_node, t),
        cubic(&p.long_perih, t),
    )
}

/// Computes the mean orbital elements of a planet, referred to the fixed
/// ecliptic and the standard equinox of J2000.0.
///
/// The angles come from Meeus, *Astronomical Algorithms*, 2nd ed.,
/// chapter 31, table 31.B; the semimajor axis and the eccentricity, which
/// do not depend on the reference frame, are those of table 31.A (see
/// [`mean_orb_elements_of_date`]). Elements in a fixed frame are what one
/// needs to compare a planet's orbit with the orbit of a comet or minor
/// planet given for the equinox J2000.0.
///
/// For [`Planet::Earth`] the inclination is measured on the fixed ecliptic
/// of J2000.0. It is zero at J2000.0 and is counted negative before
/// A.D. 2000, which keeps the inclination and the node continuous through
/// that epoch (instead of a sign flip of `i` and a 180° jump of `Ω`).
///
/// # Arguments
///
/// * `planet`: Any variant of [`Planet`]
/// * `JD`    : Julian Ephemeris Day (TT)
pub fn mean_orb_elements_j2000(planet: &Planet, JD: f64) -> OrbElements {
    let p = elements_of_date(planet);
    let q = angles_j2000(planet);
    let t = time::julian_cent(JD);

    orb_elements_frm_deg(
        cubic(&q.mn_long, t),
        cubic(&p.semimaj_axis, t),
        cubic(&p.ecc, t),
        cubic(&q.inc, t),
        cubic(&q.long_ascend_node, t),
        cubic(&q.long_perih, t),
    )
}

/// Computes a planet's mean orbital elements, referred to the mean equinox
/// of the date
///
/// This is the tuple form of [`mean_orb_elements_of_date`], which should
/// be preferred in new code. The only difference is that the inclination is
/// reduced to `[0, 2π)` here (it is always positive for this frame).
///
/// # Returns
///
/// `(L, a, e, i, omega, pi, M, w)`
///
/// * `L`: Mean longitude *| in radians*
/// * `a`: Semimajor axis of the orbit *| in AU*
/// * `e`: Eccentricity of the orbit
/// * `i`: Inclination of the plane of the orbit with the plane of
///   the Earth's ecliptic *| in radians*
/// * `omega`: Longitude of the ascending node *| in radians*.
///   Undefined for `Planet::Earth`; 0 is returned.
/// * `pi`: Longitude of the perihelion *| in radians*
/// * `M`: Mean anomaly *| in radians*
/// * `w`: Argument of the perihelion *| in radians*.
///   For `Planet::Earth` this equals `pi`.
///
/// # Arguments
///
/// * `planet`: Any variant of [Planet](./enum.Planet.html)
/// * `JD`    : Julian Ephemeris Day (TT)
pub fn orb_elements(planet: &Planet, JD: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let el = mean_orb_elements_of_date(planet, JD);

    (
        el.mn_long,
        el.semimaj_axis,
        el.ecc,
        angle::limit_to_two_PI(el.inc),
        el.long_ascend_node,
        el.long_perih,
        el.mn_anom,
        el.arg_perih,
    )
}

// ---------------------------------------------------------------------------
// Heliocentric positions from VSOP87 (Meeus, chapter 32 and appendix III)
// ---------------------------------------------------------------------------

/// The VSOP87D series of a planet: `[L, B, R]`, each a list of the series
/// multiplying successive powers of τ.
#[inline]
fn vsop87d_series(planet: &Planet) -> &'static [&'static [&'static [[f64; 3]]]; 3] {
    match *planet {
        Planet::Mercury => &VSOPD_87::mercury::TERMS,
        Planet::Venus => &VSOPD_87::venus::TERMS,
        Planet::Earth => &VSOPD_87::earth::TERMS,
        Planet::Mars => &VSOPD_87::mars::TERMS,
        Planet::Jupiter => &VSOPD_87::jupiter::TERMS,
        Planet::Saturn => &VSOPD_87::saturn::TERMS,
        Planet::Uranus => &VSOPD_87::uranus::TERMS,
        Planet::Neptune => &VSOPD_87::neptune::TERMS,
    }
}

/// Evaluates one VSOP87 coordinate: `Σ_k τ^k Σ_j A cos(B + C τ)`
/// (Meeus, eq. 32.2, with the amplitudes already in radians or AU).
#[inline]
fn eval_vsop87(series: &[&[[f64; 3]]], tau: f64) -> f64 {
    let mut sum = 0.0;
    let mut tau_k = 1.0;

    for terms in series {
        let mut y = 0.0;
        for t in terms.iter() {
            y += t[0] * (t[1] + t[2] * tau).cos();
        }
        sum += y * tau_k;
        tau_k *= tau;
    }

    sum
}

/// Computes a planet's heliocentric coordinates, referred to the mean
/// dynamical ecliptic and equinox of the date
///
/// The coordinates are computed from the **complete** VSOP87 theory of
/// Bretagnon and Francou (1988), version D (heliocentric spherical
/// coordinates, ecliptic and equinox of date), as described in Meeus,
/// *Astronomical Algorithms*, 2nd ed., chapter 32. Meeus' appendix III
/// prints only an abridged selection of these terms; using all of them
/// reproduces the book's examples to well within the error of the abridged
/// series. The precision of the full theory is about 1″ over 4000 years
/// around J2000.0 for Mercury, Venus, the Earth and Mars, over 2000 years
/// for Jupiter and Saturn, and over 6000 years for Uranus and Neptune
/// (Bretagnon & Francou 1988); outside these spans it degrades slowly.
///
/// The frame is the VSOP87 dynamical ecliptic, which differs from FK5 by
/// a fraction of an arcsecond; see [`ecl_coords_to_FK5`].
///
/// # Returns
///
/// `(long, lat, rad_vec)`
///
/// * `long`   : Heliocentric longitude *| in radians*, in `[0, 2π)`
/// * `lat`    : Heliocentric latitude *| in radians*, in `[-π/2, π/2]`
///   (negative south of the ecliptic)
/// * `rad_vec`: Heliocentric radius vector *| in AU*
///
/// # Arguments
///
/// * `planet`: Any variant of [Planet](./enum.Planet.html)
/// * `JD`    : Julian Ephemeris Day (TT)
///
/// # Example
///
/// ```
/// use astro::planet::{self, Planet};
///
/// // Meeus, example 32.a: Venus on 1992 December 20 at 0h TD. The book
/// // uses the abridged series (L = 26.11428°, B = -2.62070°,
/// // R = 0.724603 AU); the complete series differ by less than 1".
/// let (l, b, r) = planet::heliocent_coords(&Planet::Venus, 2448976.5);
/// assert!((l.to_degrees() - 26.11428).abs() < 1.0 / 3600.0);
/// assert!((b.to_degrees() - -2.62070).abs() < 1.0 / 3600.0);
/// assert!((r - 0.724603).abs() < 2e-6);
/// ```
pub fn heliocent_coords(planet: &Planet, JD: f64) -> (f64, f64, f64) {
    let series = vsop87d_series(planet);
    let tau = time::julian_mill(JD);

    let L = angle::limit_to_two_PI(eval_vsop87(series[0], tau));
    // The latitude is a small signed angle and is not reduced.
    let B = eval_vsop87(series[1], tau);
    let R = eval_vsop87(series[2], tau);

    (L, B, R)
}

/// Converts ecliptic coordinates from the VSOP87 dynamical frame to the
/// FK5 system
///
/// Applies the small corrections of Meeus, eq. 32.3. They are well under
/// an arcsecond, so they only matter together with the complete VSOP87
/// series. They apply equally to heliocentric coordinates (`L`, `B`) and to
/// geocentric ones (`λ`, `β`).
///
/// # Returns
///
/// `(ecl_long_FK5, ecl_lat_FK5)`
///
/// * `ecl_long_FK5`: Ecliptic longitude in the FK5 system *| in radians*
/// * `ecl_lat_FK5` : Ecliptic latitude in the FK5 system *| in radians*
///
/// # Arguments
///
/// * `JD`      : Julian Ephemeris Day (TT)
/// * `ecl_long`: Ecliptic longitude on `JD`, referred to the VSOP87
///   ecliptic and mean equinox of the date *| in radians*
/// * `ecl_lat` : Ecliptic latitude on `JD`, referred to the VSOP87
///   ecliptic and mean equinox of the date *| in radians*
pub fn ecl_coords_to_FK5(JD: f64, ecl_long: f64, ecl_lat: f64) -> (f64, f64) {
    let JC = time::julian_cent(JD);
    let lambda1 = ecl_long - JC * (1.397 + JC * 0.00031).to_radians();
    let x = angle::deg_frm_dms(0, 0, 0.03916).to_radians();

    let ecl_long_correction = -angle::deg_frm_dms(0, 0, 0.09033).to_radians()
        + x * (lambda1.cos() + lambda1.sin()) * ecl_lat.tan();

    (
        ecl_long + ecl_long_correction,
        ecl_lat + x * (lambda1.cos() - lambda1.sin()),
    )
}

// ---------------------------------------------------------------------------
// Geocentric positions of the planets (Meeus, chapter 33, first method)
// ---------------------------------------------------------------------------

/// Light-time per astronomical unit *| in days* (Meeus, eq. 33.3).
const LIGHT_TIME_PER_AU: f64 = 0.0057755183;

/// Convergence threshold of the light-time iteration *| in days* (about
/// 10 µs; the position of the fastest planet changes by less than 1e-9″).
const LIGHT_TIME_TOL: f64 = 1e-10;

/// Upper bound on the light-time iterations. Each iteration reduces the
/// error by a factor of order v/c (< 1e-3), so 3 or 4 suffice.
const LIGHT_TIME_MAX_ITER: usize = 10;

/// Time taken by light to cross `dist` AU *| in days* (Meeus, eq. 33.3).
#[inline(always)]
fn light_time(dist: f64) -> f64 {
    LIGHT_TIME_PER_AU * dist
}

/// Computes a planet's geocentric, geometric ecliptic position,
/// uncorrected for light-time
///
/// Combines heliocentric coordinates of a planet and of the Earth, taken
/// for the same instant, with Meeus, eqs. 33.1, 33.2, 33.4 and 33.3. The
/// result is the planet's true (geometric) direction and distance, in the
/// frame of the inputs (normally the VSOP87 ecliptic and mean equinox of
/// the date, see [`heliocent_coords`]).
///
/// # Returns
///
/// `(ecl_long, ecl_lat, rad_vec, light_time)`
///
/// * `ecl_long`: Geometric longitude of the planet *| in radians*, in
///   `[0, 2π)`
/// * `ecl_lat`: Geometric latitude of the planet *| in radians*
/// * `rad_vec`: Earth-planet distance *| in AU*
/// * `light_time`: Time taken by light to travel this distance *| in
///   days*
///
/// # Arguments
///
/// * `L0`: Heliocentric longitude of the Earth *| in radians*
/// * `B0`: Heliocentric latitude of the Earth *| in radians*
/// * `R0`: Heliocentric radius vector of the Earth *| in AU*
/// * `L` : Heliocentric longitude of the planet *| in radians*
/// * `B` : Heliocentric latitude of the planet *| in radians*
/// * `R` : Heliocentric radius vector of the planet *| in AU*
pub fn geocent_geomet_ecl_coords(
    L0: f64,
    B0: f64,
    R0: f64,
    L: f64,
    B: f64,
    R: f64,
) -> (f64, f64, f64, f64) {
    let (x, y, z) = geocent_ecl_rect_coords(L0, B0, R0, L, B, R);

    let (lambda, beta) = ecl_coords_frm_ecl_rect_coords(x, y, z);
    let planet_earth_dist = dist_frm_ecl_rect_coords(x, y, z);
    let light_time = light_time(planet_earth_dist);

    (
        angle::limit_to_two_PI(lambda),
        beta,
        planet_earth_dist,
        light_time,
    )
}

/// Computes a planet's geocentric ecliptic rectangular coordinates from
/// its heliocentric position and the Earth's (Meeus, eq. 33.1)
///
/// # Returns
///
/// `(x, y, z)` *| in AU*
///
/// # Arguments
///
/// * `L0`: Heliocentric longitude of the Earth *| in radians*
/// * `B0`: Heliocentric latitude of the Earth *| in radians*
/// * `R0`: Heliocentric radius vector of the Earth *| in AU*
/// * `L` : Heliocentric longitude of the planet *| in radians*
/// * `B` : Heliocentric latitude of the planet *| in radians*
/// * `R` : Heliocentric radius vector of the planet *| in AU*
fn geocent_ecl_rect_coords(L0: f64, B0: f64, R0: f64, L: f64, B: f64, R: f64) -> (f64, f64, f64) {
    let x = R * B.cos() * L.cos() - R0 * B0.cos() * L0.cos();
    let y = R * B.cos() * L.sin() - R0 * B0.cos() * L0.sin();
    let z = R * B.sin() - R0 * B0.sin();

    (x, y, z)
}

/// Computes ecliptic longitude and latitude from rectangular coordinates
/// (Meeus, eq. 33.2)
///
/// # Returns
///
/// `(ecl_long, ecl_lat)` *| in radians*; the longitude is in `(-π, π]`
///
/// # Arguments
///
/// * `x`, `y`, `z`: Rectangular ecliptic coordinates
#[inline]
fn ecl_coords_frm_ecl_rect_coords(x: f64, y: f64, z: f64) -> (f64, f64) {
    (y.atan2(x), z.atan2((x * x + y * y).sqrt()))
}

/// Computes the length of a vector from its rectangular coordinates
/// (Meeus, eq. 33.4)
///
/// # Arguments
///
/// * `x`, `y`, `z`: Rectangular coordinates *| in AU*
#[inline]
fn dist_frm_ecl_rect_coords(x: f64, y: f64, z: f64) -> f64 {
    (x * x + y * y + z * z).sqrt()
}

/// Geocentric rectangular ecliptic coordinates of a planet seen from the
/// Earth's heliocentric position `earth`, corrected for light-time by
/// iteration (Meeus, chapter 33).
///
/// Returns `(geometric (x, y, z), light-time corrected (x, y, z), τ)`. The
/// planet is taken at `JD - τ` where `τ` is the light-time for the
/// distance obtained with that same position.
fn light_time_corrected_rect(
    planet: &Planet,
    JD: f64,
    earth: (f64, f64, f64),
) -> ((f64, f64, f64), (f64, f64, f64), f64) {
    let (L0, B0, R0) = earth;
    let rect_at = |jd: f64| {
        let (L, B, R) = heliocent_coords(planet, jd);
        geocent_ecl_rect_coords(L0, B0, R0, L, B, R)
    };

    let geomet = rect_at(JD);
    let mut rect = geomet;
    let mut tau = light_time(dist_frm_ecl_rect_coords(rect.0, rect.1, rect.2));

    for _ in 0..LIGHT_TIME_MAX_ITER {
        rect = rect_at(JD - tau);
        let new_tau = light_time(dist_frm_ecl_rect_coords(rect.0, rect.1, rect.2));
        let converged = (new_tau - tau).abs() < LIGHT_TIME_TOL;
        tau = new_tau;
        if converged {
            break;
        }
    }

    (geomet, rect, tau)
}

/// Computes a planet's geocentric ecliptic position, corrected for
/// light-time only
///
/// The planet is taken where it was when the light that reaches the Earth
/// at `JD` left it (Meeus, chapter 33; the light-time is iterated to
/// convergence). The result is **not** corrected for aberration, nutation
/// or the FK5 frame, and is referred to the VSOP87 ecliptic and mean
/// equinox of the date. For the fully reduced apparent place use
/// [`geocent_apprnt_pos`].
///
/// # Returns
///
/// `(planet_ecl_point, rad_vec)`
///
/// * `planet_ecl_point`: Ecliptic point of the planet *| in radians*;
///   the longitude is in `[0, 2π)`
/// * `rad_vec`         : Distance from the Earth (at `JD`) to the planet
///   (at `JD` minus the light-time) *| in AU*
///
/// For `Planet::Earth` the result is meaningless (zero distance).
///
/// # Arguments
///
/// * `planet`: Any variant of [Planet](./enum.Planet.html)
/// * `JD`    : Julian Ephemeris Day (TT)
pub fn geocent_apprnt_ecl_coords(planet: &Planet, JD: f64) -> (coords::EclPoint, f64) {
    let earth = heliocent_coords(&Planet::Earth, JD);
    let (_, (x, y, z), _) = light_time_corrected_rect(planet, JD, earth);

    let (long, lat) = ecl_coords_frm_ecl_rect_coords(x, y, z);
    let ecl_point = coords::EclPoint {
        long: angle::limit_to_two_PI(long),
        lat,
    };

    (ecl_point, dist_frm_ecl_rect_coords(x, y, z))
}

/// Geocentric position of a planet, with the intermediate results of the
/// reduction of Meeus, chapter 33 (first method).
///
/// All angles are in radians; longitudes and right ascensions are in
/// `[0, 2π)`, latitudes and declinations in `[-π/2, π/2]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeocentPos {
    /// Geometric (true) ecliptic longitude at the given instant, VSOP87
    /// ecliptic and mean equinox of the date
    pub geomet_long: f64,
    /// Geometric (true) ecliptic latitude, same frame
    pub geomet_lat: f64,
    /// True Earth-planet distance at the given instant *| in AU*
    pub geomet_dist: f64,
    /// Light-time `τ` *| in days*
    pub light_time: f64,
    /// Distance between the Earth at `JD` and the planet at `JD - τ`,
    /// i.e. the distance travelled by the light that is seen *| in AU*
    pub dist: f64,
    /// Ecliptic longitude corrected for light-time only, VSOP87 ecliptic
    /// and mean equinox of the date
    pub light_time_long: f64,
    /// Ecliptic latitude corrected for light-time only, same frame
    pub light_time_lat: f64,
    /// Apparent ecliptic longitude: corrected for light-time, aberration
    /// and nutation, in the FK5 system, true equinox of the date
    pub apprnt_long: f64,
    /// Apparent ecliptic latitude, same corrections
    pub apprnt_lat: f64,
    /// Apparent right ascension, true equator and equinox of the date
    pub apprnt_asc: f64,
    /// Apparent declination, true equator and equinox of the date
    pub apprnt_dec: f64,
    /// Apparent elongation of the planet from the Sun, in `[0, π]`
    /// (Meeus, eq. 33.5)
    pub elongation: f64,
}

impl GeocentPos {
    /// The apparent ecliptic position as an [`EclPoint`](coords::EclPoint)
    pub fn apprnt_ecl_point(&self) -> coords::EclPoint {
        coords::EclPoint {
            long: self.apprnt_long,
            lat: self.apprnt_lat,
        }
    }

    /// The apparent equatorial position as an [`EqPoint`](coords::EqPoint)
    pub fn apprnt_eq_point(&self) -> coords::EqPoint {
        coords::EqPoint {
            asc: self.apprnt_asc,
            dec: self.apprnt_dec,
        }
    }
}

/// Annual aberration in ecliptic coordinates (Meeus, eq. 23.2), including
/// the terms in the eccentricity of the Earth's orbit.
///
/// `sun_long` is the Sun's true geometric longitude; returns
/// `(Δλ, Δβ)` *| in radians*.
fn ecl_aberr(long: f64, lat: f64, sun_long: f64, JD: f64) -> (f64, f64) {
    let t = time::julian_cent(JD);
    let kappa = angle::deg_frm_dms(0, 0, 20.49552).to_radians();
    let e = 0.016708634 - t * (0.000042037 + t * 0.0000001267);
    let perih = (102.93735 + t * (1.71946 + t * 0.00046)).to_radians();

    let d_long = kappa * (-(sun_long - long).cos() + e * (perih - long).cos()) / lat.cos();
    let d_lat = -kappa * lat.sin() * ((sun_long - long).sin() - e * (perih - long).sin());

    (d_long, d_lat)
}

/// Computes a planet's apparent geocentric position
///
/// Follows the first method of Meeus, *Astronomical Algorithms*, 2nd ed.,
/// chapter 33, with the complete VSOP87 series for the planet and the
/// Earth:
///
/// 1. heliocentric coordinates of the Earth at `JD` and of the planet at
///    `JD - τ`, with the light-time `τ` iterated to convergence
///    (eqs. 33.1–33.4);
/// 2. annual aberration in ecliptic coordinates (eq. 23.2);
/// 3. conversion to the FK5 system (eq. 32.3);
/// 4. nutation in longitude (chapter 22) and conversion to right ascension
///    and declination with the true obliquity (eqs. 13.3, 13.4).
///
/// The elongation from the Sun (eq. 33.5) uses the Sun's apparent
/// longitude computed from the same Earth position.
///
/// For Venus on 1992 December 20 (Meeus, example 33.a) the result agrees
/// with the book's values from the complete VSOP87 theory to 0.001 s in
/// right ascension and 0.01″ in declination, and the positions of Mercury
/// and Venus in Meeus' example 18.a (taken from an accurate ephemeris)
/// are matched to the same level. Near J2000 the accuracy is thus a few
/// hundredths of an arcsecond; far from it, it is limited by VSOP87 (see
/// [`heliocent_coords`]) and by the nutation and obliquity models.
///
/// # Returns
///
/// `None` for [`Planet::Earth`], otherwise the position and the
/// intermediate results of the reduction.
///
/// # Arguments
///
/// * `planet`: Any variant of [`Planet`] except `Planet::Earth`
/// * `JD`    : Julian Ephemeris Day (TT)
///
/// # Example
///
/// ```
/// use astro::planet::{self, Planet};
///
/// // Meeus, example 33.a: Venus on 1992 December 20 at 0h TD
/// let pos = planet::geocent_apprnt_pos(&Planet::Venus, 2448976.5).unwrap();
///
/// // Book, complete VSOP87: α = 21h04m41.454s, δ = -18°53'16.84"
/// let asc_hours = pos.apprnt_asc.to_degrees() / 15.0;
/// assert!((asc_hours - (21.0 + 4.0 / 60.0 + 41.454 / 3600.0)).abs() < 0.02 / 3600.0);
/// let dec_deg = pos.apprnt_dec.to_degrees();
/// assert!((dec_deg - -(18.0 + 53.0 / 60.0 + 16.84 / 3600.0)).abs() < 0.2 / 3600.0);
/// assert!((pos.dist - 0.910947).abs() < 1e-6);
/// ```
pub fn geocent_apprnt_pos(planet: &Planet, JD: f64) -> Option<GeocentPos> {
    if *planet == Planet::Earth {
        return None;
    }

    let earth = heliocent_coords(&Planet::Earth, JD);
    let (L0, _, R0) = earth;
    let (geomet, rect, tau) = light_time_corrected_rect(planet, JD, earth);

    let (geomet_long, geomet_lat) = ecl_coords_frm_ecl_rect_coords(geomet.0, geomet.1, geomet.2);
    let geomet_dist = dist_frm_ecl_rect_coords(geomet.0, geomet.1, geomet.2);
    let (lt_long, lt_lat) = ecl_coords_frm_ecl_rect_coords(rect.0, rect.1, rect.2);
    let dist = dist_frm_ecl_rect_coords(rect.0, rect.1, rect.2);

    // Aberration, with the Sun's geometric longitude (VSOP87 frame).
    let sun_long = L0 + std::f64::consts::PI;
    let (d_long, d_lat) = ecl_aberr(lt_long, lt_lat, sun_long, JD);

    // FK5, then nutation in longitude.
    let (long_fk5, lat_fk5) = ecl_coords_to_FK5(JD, lt_long + d_long, lt_lat + d_lat);
    let (nut_long, nut_oblq) = nutation::nutation(JD);
    let apprnt_long = angle::limit_to_two_PI(long_fk5 + nut_long);
    let apprnt_lat = lat_fk5;

    let true_oblq = ecliptic::mn_oblq_IAU(JD) + nut_oblq;
    let apprnt_asc =
        angle::limit_to_two_PI(coords::asc_frm_ecl(apprnt_long, apprnt_lat, true_oblq));
    let apprnt_dec = coords::dec_frm_ecl(apprnt_long, apprnt_lat, true_oblq);

    // Sun's apparent longitude: FK5 offset, aberration, nutation.
    let sun_apprnt_long = sun_long
        - angle::deg_frm_dms(0, 0, 0.09033).to_radians()
        - angle::deg_frm_dms(0, 0, 20.4898).to_radians() / R0
        + nut_long;
    let cos_elong = apprnt_lat.cos() * (apprnt_long - sun_apprnt_long).cos();

    Some(GeocentPos {
        geomet_long: angle::limit_to_two_PI(geomet_long),
        geomet_lat,
        geomet_dist,
        light_time: tau,
        dist,
        light_time_long: angle::limit_to_two_PI(lt_long),
        light_time_lat: lt_lat,
        apprnt_long,
        apprnt_lat,
        apprnt_asc,
        apprnt_dec,
        elongation: cos_elong.clamp(-1.0, 1.0).acos(),
    })
}

/// Computes the elongation of a planet from the Sun
///
/// Meeus, eq. 33.5: `cos ψ = cos β cos(λ - λ0)`, where the Sun's
/// latitude (always below 1.2″) is neglected.
///
/// # Returns
///
/// * `elongation`: Angular distance between the planet and the Sun *| in
///   radians*, in `[0, π]`
///
/// # Arguments
///
/// * `planet_ecl_point`: Apparent geocentric ecliptic position of the
///   planet *| in radians*
/// * `sun_long`        : Apparent geocentric ecliptic longitude of the
///   Sun *| in radians*
pub fn elongation(planet_ecl_point: &coords::EclPoint, sun_long: f64) -> f64 {
    (planet_ecl_point.lat.cos() * (planet_ecl_point.long - sun_long).cos())
        .clamp(-1.0, 1.0)
        .acos()
}

// ---------------------------------------------------------------------------
// Positions from orbital elements (Meeus, chapter 33, second method)
// ---------------------------------------------------------------------------

/// Computes the geocentric equatorial coordinates of a body from its
/// orbital elements (Meeus, chapter 33, second method)
///
/// Implements eqs. 33.7–33.10: the heliocentric rectangular equatorial
/// coordinates of the body are formed from its elements and added to the
/// geocentric rectangular equatorial coordinates of the Sun. The elements
/// and the Sun's coordinates must refer to the same equinox; the result
/// is then referred to it too. To correct for light-time, call again with
/// `v` and `r` computed for `JD - light_time`, keeping the Sun's
/// coordinates for `JD` (the result is then the astrometric position).
///
/// # Returns
///
/// `(asc, dec, light_time)`
///
/// * `asc`       : Geocentric right ascension *| in radians*, in
///   `[0, 2π)`
/// * `dec`       : Geocentric declination *| in radians*
/// * `light_time`: Light-time for the Earth-body distance `Δ` *| in days*
///   (`Δ = light_time / 0.0057755183`)
///
/// # Arguments
///
/// * `X`, `Y`, `Z` : Geocentric rectangular equatorial coordinates of the
///   Sun *| in AU*
/// * `i`         : Inclination of the orbit *| in radians*
/// * `w`         : Argument of the perihelion *| in radians*
/// * `sigma`     : Longitude of the ascending node *| in radians*
/// * `oblq_eclip`: Obliquity of the ecliptic for the equinox of the
///   elements *| in radians* (23°26′21.448″ for J2000.0)
/// * `v`         : True anomaly of the body *| in radians*
/// * `r`         : Radius vector of the body *| in AU*
pub fn geocent_eq_coords(
    X: f64,
    Y: f64,
    Z: f64,
    i: f64,
    w: f64,
    sigma: f64,
    oblq_eclip: f64,
    v: f64,
    r: f64,
) -> (f64, f64, f64) {
    let F = sigma.cos();
    let G = sigma.sin() * oblq_eclip.cos();
    let H = sigma.sin() * oblq_eclip.sin();

    let P = -sigma.sin() * i.cos();
    let Q = sigma.cos() * i.cos() * oblq_eclip.cos() - i.sin() * oblq_eclip.sin();
    let R = sigma.cos() * i.cos() * oblq_eclip.sin() + i.sin() * oblq_eclip.cos();

    let A = F.atan2(P);
    let B = G.atan2(Q);
    let C = H.atan2(R);
    let a = (F * F + P * P).sqrt();
    let b = (G * G + Q * Q).sqrt();
    let c = (H * H + R * R).sqrt();

    let x = r * a * (A + w + v).sin();
    let y = r * b * (B + w + v).sin();
    let z = r * c * (C + w + v).sin();

    let xi = X + x;
    let nu = Y + y;
    let et = Z + z;

    let asc = angle::limit_to_two_PI(nu.atan2(xi));
    let dec = et.atan2((xi * xi + nu * nu).sqrt());
    let dist = (xi * xi + nu * nu + et * et).sqrt();

    (asc, dec, light_time(dist))
}

/// Computes the heliocentric ecliptic coordinates of a body from its
/// orbital elements
///
/// Meeus, chapter 33, "heliocentric ecliptical coordinates": with
/// `u = ω + v`, the rectangular coordinates are
/// `x = r (cos Ω cos u - sin Ω sin u cos i)`,
/// `y = r (sin Ω cos u + cos Ω sin u cos i)`, `z = r sin i sin u`.
/// The coordinates refer to the same ecliptic and equinox as the elements.
///
/// # Returns
///
/// `(long, lat)`
///
/// * `long`: Heliocentric ecliptic longitude *| in radians*, in `[0, 2π)`
/// * `lat` : Heliocentric ecliptic latitude *| in radians*
///
/// The radius vector is the input `r`.
///
/// # Arguments
///
/// * `i`    : Inclination of the orbit *| in radians*
/// * `sigma`: Longitude of the ascending node *| in radians*
/// * `w`    : Argument of the perihelion *| in radians*
/// * `v`    : True anomaly *| in radians*
/// * `r`    : Radius vector *| in AU*
pub fn heliocent_coords_frm_orb_elements(i: f64, sigma: f64, w: f64, v: f64, r: f64) -> (f64, f64) {
    let u = w + v;
    let x = r * (sigma.cos() * u.cos() - sigma.sin() * u.sin() * i.cos());
    let y = r * (sigma.sin() * u.cos() + sigma.cos() * u.sin() * i.cos());
    let z = r * i.sin() * u.sin();

    (
        angle::limit_to_two_PI(y.atan2(x)),
        z.atan2((x * x + y * y).sqrt()),
    )
}

// Unit tests for the private VSOP87 and chapter 33 helpers.
#[cfg(test)]
mod position_tests;

/// Computes a planet's apparent visual magnitude using G. Müller's
/// formulae
///
/// Müller derived these expressions from his observations of 1877–1891
/// (Meeus, chapter 41): `m = m0 + 5 log10(r Δ) + f(i)`. The phase term `f`
/// is evaluated with the phase angle in degrees; Jupiter, Uranus and
/// Neptune have none. Do not quote the result to better than about 0.1
/// magnitude.
///
/// # Returns
///
/// * `Ok(app_mag)`: Apparent visual magnitude of the planet
/// * `Err(..)` for `Planet::Earth` (no formula) and for `Planet::Saturn`,
///   whose magnitude depends on the ring: use
///   [`saturn::apprnt_mag_muller`] instead.
///
/// # Arguments
///
/// * `planet`: Any variant of [Planet](./enum.Planet.html)
/// * `i`     : Phase angle of the planet *| in radians* (see
///   [`phase_angl`])
/// * `delta` : Planet-Earth distance *| in AU*
/// * `r`     : Planet-Sun distance *| in AU*
pub fn apprnt_mag_muller(planet: &Planet, i: f64, delta: f64, r: f64) -> Result<f64, &'static str> {
    let x = 5.0 * (r * delta).log10();
    let i = i.to_degrees();

    match *planet {
        Planet::Mercury => Ok(x + 1.16 + (i - 50.0) * (0.02838 + (i - 50.0) * 0.0001023)),
        Planet::Venus => Ok(x - 4.0 + i * (0.01322 + i * i * 0.0000004247)),
        Planet::Earth => {
            Err("Planet::Earth was passed to the function planet::apprnt_mag_muller()")
        }
        Planet::Mars => Ok(x - 1.3 + i * 0.01486),
        Planet::Jupiter => Ok(x - 8.93),
        Planet::Saturn => Err(
            "Planet::Saturn was passed to the function planet::apprnt_mag_muller(). Use the function planet::saturn::apprnt_mag_muller() instead.",
        ),
        Planet::Uranus => Ok(x - 6.85),
        Planet::Neptune => Ok(x - 7.05),
    }
}

/// Computes a planet's apparent visual magnitude using the expressions of
/// the *Astronomical Almanac* adopted in 1984
///
/// `m = m0 + 5 log10(r Δ) + f(i)` with the phase term `f` evaluated with
/// the phase angle in degrees (Meeus, chapter 41); Uranus and Neptune have
/// no phase term. For Pluto see [`crate::pluto::apprnt_mag_84`].
///
/// # Returns
///
/// * `Ok(app_mag)`: Apparent visual magnitude of the planet
/// * `Err(..)` for `Planet::Earth` (no formula) and for `Planet::Saturn`,
///   whose magnitude depends on the ring: use
///   [`saturn::apprnt_mag_84`] instead.
///
/// # Arguments
///
/// * `planet`: Any variant of [Planet](./enum.Planet.html)
/// * `i`     : Phase angle of the planet *| in radians* (see
///   [`phase_angl`])
/// * `delta` : Planet-Earth distance *| in AU*
/// * `r`     : Planet-Sun distance *| in AU*
pub fn apprnt_mag_84(planet: &Planet, i: f64, delta: f64, r: f64) -> Result<f64, &'static str> {
    let x = 5.0 * (r * delta).log10();
    let i = i.to_degrees();

    match *planet {
        Planet::Mercury => Ok(x - 0.42 + i * (0.0380 - i * (0.000273 - i * 0.00000200))),
        Planet::Venus => Ok(x - 4.40 + i * (0.0009 + i * (0.000239 - i * 0.00000065))),
        Planet::Earth => Err("Planet::Earth was passed to the function planet::apprnt_mag_84()"),
        Planet::Mars => Ok(x - 1.52 + i * 0.016),
        Planet::Jupiter => Ok(x - 9.4 + i * 0.005),
        Planet::Saturn => Err(
            "Planet::Saturn was passed to the function planet::apprnt_mag_84(). Use the function planet::saturn::apprnt_mag_84() instead.",
        ),
        Planet::Uranus => Ok(x - 7.19),
        Planet::Neptune => Ok(x - 6.87),
    }
}

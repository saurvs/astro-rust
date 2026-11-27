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

//! Angles: angular separation, position angle, sexagesimal notation and
//! normalisation
//!
//! The angular separation and position angle follow Meeus, *Astronomical
//! Algorithms*, 2nd ed., chapter 17. The handling of negative sexagesimal
//! angles follows the advice of chapter 1.
//!
//! # Example
//!
//! The angular distance between Arcturus and Spica (Meeus, example 17.a):
//!
//! ```
//! use astro::angle;
//!
//! let arcturus = (213.9154_f64.to_radians(), 19.1825_f64.to_radians());
//! let spica = (201.2983_f64.to_radians(), (-11.1614_f64).to_radians());
//!
//! let d = angle::anglr_sepr(arcturus.0, arcturus.1, spica.0, spica.1);
//! assert!((d.to_degrees() - 32.7930).abs() < 5e-5);
//!
//! // In sexagesimal notation, to the nearest arcsecond (the book quotes
//! // 32°48′, to the nearest arcminute)
//! assert_eq!(format!("{:.0}", angle::Dms::from_deg(d.to_degrees())), "32°47′35″");
//! ```

use std::fmt;

use crate::coords;

/// 2π, the number of radians in a full turn
pub const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

/// A point on a sphere, given by a longitude-like and a latitude-like angle
///
/// This lets the functions of [`conjunction`](crate::conjunction) and
/// [`geometry`](crate::geometry) accept equatorial coordinates
/// ([`EqPoint`](crate::coords::EqPoint): right ascension, declination),
/// ecliptic coordinates ([`EclPoint`](crate::coords::EclPoint): longitude,
/// latitude) or plain `(longitude, latitude)` tuples. All the points passed
/// to one call must be in the same coordinate system and referred to the
/// same equinox.
///
/// The longitude must increase eastward (as right ascension and ecliptic
/// longitude do); orientation-dependent results, such as position angles,
/// assume this.
pub trait SphericalPoint {
    /// Longitude-like angle (right ascension or ecliptic longitude)
    /// *| in radians*
    fn lon(&self) -> f64;
    /// Latitude-like angle (declination or ecliptic latitude)
    /// *| in radians*
    fn lat(&self) -> f64;
}

impl SphericalPoint for coords::EqPoint {
    #[inline]
    fn lon(&self) -> f64 {
        self.asc
    }
    #[inline]
    fn lat(&self) -> f64 {
        self.dec
    }
}

impl SphericalPoint for coords::EclPoint {
    #[inline]
    fn lon(&self) -> f64 {
        self.long
    }
    #[inline]
    fn lat(&self) -> f64 {
        self.lat
    }
}

/// `(longitude, latitude)` *| in radians*
impl SphericalPoint for (f64, f64) {
    #[inline]
    fn lon(&self) -> f64 {
        self.0
    }
    #[inline]
    fn lat(&self) -> f64 {
        self.1
    }
}

impl<T: SphericalPoint + ?Sized> SphericalPoint for &T {
    #[inline]
    fn lon(&self) -> f64 {
        (**self).lon()
    }
    #[inline]
    fn lat(&self) -> f64 {
        (**self).lat()
    }
}

/// `sin²(θ/2)`, the haversine of `θ`
#[inline]
fn hav(theta: f64) -> f64 {
    let s = (0.5 * theta).sin();
    s * s
}

/// Pauwels' auxiliary quantities `(x, y, z)` for the angular separation
/// (Meeus, ch. 17): `z` is the cosine of the separation and `√(x² + y²)`
/// its sine.
///
/// `x` is evaluated in a rearranged form, `sin(δ₂ − δ₁) + 2 sin δ₁ cos δ₂
/// hav(α₂ − α₁)`, which is mathematically identical to Meeus' expression
/// but does not lose digits by cancellation when the points are close.
#[inline]
fn pauwels_xyz(p1a1: f64, p1a2: f64, p2a1: f64, p2a2: f64) -> (f64, f64, f64) {
    let (sin_d1, cos_d1) = p1a2.sin_cos();
    let (sin_d2, cos_d2) = p2a2.sin_cos();
    let d_a = p2a1 - p1a1;
    let (sin_da, cos_da) = d_a.sin_cos();

    let x = (p2a2 - p1a2).sin() + 2.0 * sin_d1 * cos_d2 * hav(d_a);
    let y = cos_d2 * sin_da;
    let z = sin_d1 * sin_d2 + cos_d1 * cos_d2 * cos_da;

    (x, y, z)
}

/// Computes the angular separation between two points on the celestial
/// sphere
///
/// Valid for any pair of points, from coincident to antipodal: the result
/// is computed as `atan2(sin d, cos d)` (the method of Thierry Pauwels given
/// by Meeus, ch. 17), which keeps full precision at all separations. (The
/// plain arccosine of formula 17.1, used before version 3, lost all
/// precision below about 1″ and could return `NaN` for coincident points.)
///
/// # Returns
///
/// * `anglr_sepr`: Angular separation, in `[0, π]` *| in radians*
///
/// # Arguments
///
/// * `p1a1`: Angle 1 of point 1 *| in radians*
/// * `p1a2`: Angle 2 of point 1 *| in radians*
/// * `p2a1`: Angle 1 of point 2 *| in radians*
/// * `p2a2`: Angle 2 of point 2 *| in radians*
///
/// Angle 1 may be right ascension or longitude (any sense, any range).
/// Angle 2 may be declination or latitude.
///
/// # Example
///
/// ```
/// use astro::angle;
///
/// // Two points 1 milliarcsecond apart
/// let mas = (1e-3_f64 / 3600.0).to_radians();
/// let d = angle::anglr_sepr(0.3, 0.2, 0.3, 0.2 + mas);
/// assert!((d / mas - 1.0).abs() < 1e-7);
/// ```
#[inline]
pub fn anglr_sepr(p1a1: f64, p1a2: f64, p2a1: f64, p2a2: f64) -> f64 {
    let (x, y, z) = pauwels_xyz(p1a1, p1a2, p2a1, p2a2);

    x.hypot(y).atan2(z)
}

/// Computes the angular separation between two points by the haversine
/// formula
///
/// Meeus, ch. 17, eq. 17.5: `hav d = hav Δδ + cos δ₁ cos δ₂ hav Δα`. Very
/// accurate for small separations, but it loses precision towards 180°
/// (errors up to about 1e-8 rad for antipodal points), so prefer
/// [`anglr_sepr`] in general.
///
/// # Returns
///
/// * `anglr_sepr`: Angular separation, in `[0, π]` *| in radians*
///
/// # Arguments
///
/// Same as [`anglr_sepr`].
#[inline]
pub fn anglr_sepr_haversine(p1a1: f64, p1a2: f64, p2a1: f64, p2a2: f64) -> f64 {
    let h = hav(p1a2 - p2a2) + p1a2.cos() * p2a2.cos() * hav(p1a1 - p2a1);

    2.0 * h.clamp(0.0, 1.0).sqrt().asin()
}

/// Computes an approximate angular separation for two close points
///
/// Meeus, ch. 17, eq. 17.2: `d = √((Δα cos δ)² + Δδ²)`, where `δ` is the
/// mean of the two declinations. It is a flat-sky approximation, adequate
/// only for small separations (well under a degree) away from the poles;
/// [`anglr_sepr`] is exact and should normally be preferred. The difference
/// of angle 1 is reduced to `[-π, π)` first, so points on both sides of
/// 0h/0° are handled.
///
/// # Returns
///
/// * `anglr_sepr`: Approximate angular separation *| in radians*
///
/// # Arguments
///
/// Same as [`anglr_sepr`].
#[inline]
pub fn anglr_sepr_small(p1a1: f64, p1a2: f64, p2a1: f64, p2a2: f64) -> f64 {
    let d_a = limit_to_pm_pi(p1a1 - p2a1);
    let d_d = p1a2 - p2a2;
    let mean_d = 0.5 * (p1a2 + p2a2);

    (d_a * mean_d.cos()).hypot(d_d)
}

/// Computes the position angle of one body relative to another
///
/// Meeus, ch. 17, "relative position angle". The position angle `P` of
/// body 1 with respect to body 2 is the direction of body 1 as seen from
/// body 2, measured from the north (direction of increasing declination)
/// through the east (direction of increasing right ascension).
///
/// Meeus writes `tan P = sin Δα / (cos δ₂ tan δ₁ − sin δ₂ cos Δα)` with
/// `Δα = α₁ − α₂`. Here numerator and denominator are multiplied by
/// `cos δ₁ ≥ 0` and the quadrant is fixed by `atan2`, which makes the
/// formula valid also when body 1 is at a pole.
///
/// # Returns
///
/// * `posn_angl`: Position angle of body 1 relative to body 2, in
///   `[0, 2π)` *| in radians*. Undefined (0 or π) if the two bodies
///   coincide, and dependent on `asc2` if body 2 is at a pole.
///
/// # Arguments
///
/// * `asc1`: Right ascension (or ecliptic longitude) of body 1 *| in radians*
/// * `dec1`: Declination (or ecliptic latitude) of body 1 *| in radians*
/// * `asc2`: Right ascension (or ecliptic longitude) of body 2 *| in radians*
/// * `dec2`: Declination (or ecliptic latitude) of body 2 *| in radians*
///
/// # Example
///
/// ```
/// use astro::angle;
///
/// // A body due north of another has position angle 0, due east 90°.
/// let north = angle::relative_position_angle(1.0, 0.3, 1.0, 0.2);
/// assert!(north.abs() < 1e-12);
/// let east = angle::relative_position_angle(1.001, 0.0, 1.0, 0.0);
/// assert!((east.to_degrees() - 90.0).abs() < 1e-9);
/// ```
pub fn relative_position_angle(asc1: f64, dec1: f64, asc2: f64, dec2: f64) -> f64 {
    let d_a = asc1 - asc2;
    let (sin_d1, cos_d1) = dec1.sin_cos();
    let (sin_d2, cos_d2) = dec2.sin_cos();

    let num = cos_d1 * d_a.sin();
    let den = cos_d2 * sin_d1 - sin_d2 * cos_d1 * d_a.cos();

    limit_to_two_PI(num.atan2(den))
}

/// Computes an angle in degrees with decimals, from an angle
/// expressed in degrees, arcminutes and arcseconds
///
/// All three components of a negative angle should carry the minus sign,
/// as in `(-13, -47, -22.0)` for −13°47′22″. Following Meeus (ch. 1), a
/// negative `deg` makes the angle negative even if `min` and `sec` are
/// given positive, so `(-13, 47, 22.0)` also means −13°47′22″.
///
/// An angle between 0° and −1° has `deg = 0`, which carries no sign; give
/// the sign on the first non-zero component instead: `(0, -30, 0.0)` or
/// `(0, -30, -15.0)` for −0°30′15″ (since version 3 `(0, -30, 15.0)` means
/// the same), and `(0, 0, -20.0)` for −0°0′20″. [`Dms`] has an explicit
/// sign and avoids this ambiguity.
///
/// # Returns
///
/// * `deg`: Angle in degrees with decimals
///
/// # Arguments
///
/// * `deg`: Degrees
/// * `min`: Arcminutes
/// * `sec`: Arcseconds
#[inline]
pub fn deg_frm_dms(deg: i64, min: i64, sec: f64) -> f64 {
    let (M, S) = signed_parts(deg, min, sec);

    (deg as f64) + M / 60.0 + S / 3600.0
}

/// Applies Meeus' rule for negative sexagesimal input: if the leading
/// non-zero integer part is negative, the later parts are made negative too.
#[inline]
fn signed_parts(whole: i64, min: i64, sec: f64) -> (f64, f64) {
    if whole < 0 || (whole == 0 && min < 0) {
        (-(min as f64).abs(), -sec.abs())
    } else {
        (min as f64, sec)
    }
}

/// Splits a non-negative number of units into whole units, minutes and
/// seconds of the unit.
///
/// The result always has `min` in `0..=59` and `sec` in `[0, 60)`. If the
/// seconds come out within a few ulps of 60 — which happens for values that
/// are whole numbers of seconds but not exactly representable in binary,
/// such as 1°01′00″ — they are carried into the minutes, so that such
/// values split exactly instead of as 1°00′59.9999999999998″.
fn split_sexagesimal(a: f64) -> (f64, f64, f64) {
    let mut whole = a.trunc();
    let m = (a - whole) * 60.0;
    let mut min = m.trunc();
    let mut sec = (m - min) * 60.0;

    // A few ulps of the input, expressed in seconds
    let tol = 3600.0 * 4.0 * f64::EPSILON * (a + 1.0);
    if sec >= 60.0 - tol {
        sec = 0.0;
        min += 1.0;
        if min >= 60.0 {
            min = 0.0;
            whole += 1.0;
        }
    }

    (whole, min, sec)
}

/// Splits a signed angle into same-signed components (the legacy tuple
/// convention of `dms_frm_deg` and `hms_frm_deg`).
fn signed_split(x: f64) -> (i64, i64, f64) {
    let (w, m, s) = split_sexagesimal(x.abs());
    if x < 0.0 {
        (-(w as i64), -(m as i64), if s == 0.0 { 0.0 } else { -s })
    } else {
        (w as i64, m as i64, s)
    }
}

/// Computes an angle expressed in degrees, arcminutes and
/// arcseconds, from an angle in degrees with decimals
///
/// The components are truncated, not rounded: `min` is in `0..=59` and
/// `sec` in `[0, 60)` in magnitude. For a negative angle every non-zero
/// component is negative, so −0.5° gives `(0, -30, 0.0)` and −13.5° gives
/// `(-13, -30, 0.0)`; [`deg_frm_dms`] reads both back correctly. Use
/// [`Dms::from_deg`] for an explicit sign, or [`Dms::from_deg_rounded`] to
/// round the arcseconds for display (with carry, so 59.96″ rounded to one
/// decimal becomes 0.0″ of the next minute rather than 60.0″).
///
/// # Returns
///
/// `(deg, min, sec)`
///
/// * `deg`: Degrees
/// * `min`: Arcminutes
/// * `sec`: Arcseconds
///
/// # Arguments
///
/// * `deg`: Angle in degrees with decimals. Must be finite and less than
///   2⁶³ in magnitude for a meaningful result.
#[inline]
pub fn dms_frm_deg(deg: f64) -> (i64, i64, f64) {
    signed_split(deg)
}

/// Computes an angle in degrees with decimals, from an angle
/// expressed in hours, minutes and seconds
///
/// The sign rule is the same as for [`deg_frm_dms`]: a negative `hour`
/// (or, if `hour` is zero, a negative `min`) makes the whole angle negative.
/// Before version 3 `(-1, 30, 0.0)` was read as −1ʰ + 30ᵐ = −0ʰ30ᵐ; it now
/// means −1ʰ30ᵐ, consistently with [`deg_frm_dms`] and [`hms_frm_deg`].
///
/// # Returns
///
/// * `deg`: Angle in degrees with decimals
///
/// # Arguments
///
/// * `hour`: Hours
/// * `min`: Minutes
/// * `sec`: Seconds
#[inline]
pub fn deg_frm_hms(hour: i64, min: i64, sec: f64) -> f64 {
    let (M, S) = signed_parts(hour, min, sec);

    15.0 * ((hour as f64) + M / 60.0 + S / 3600.0)
}

/// Computes an angle expressed in hours, minutes and
/// seconds, from an angle in degrees with decimals
///
/// Truncated components with the same sign convention as [`dms_frm_deg`].
/// The angle is not reduced to `[0, 24ʰ)`; call [`limit_to_360`] first if
/// that is wanted. See [`Hms`] for an explicit sign and rounding.
///
/// # Returns
///
/// `(hour, min, sec)`
///
/// * `hour`: Hours
/// * `min`: Minutes
/// * `sec`: Seconds
///
/// # Arguments
///
/// * `deg`: Angle in degrees with decimals. Must be finite.
#[inline]
pub fn hms_frm_deg(deg: f64) -> (i64, i64, f64) {
    signed_split(deg / 15.0)
}

/// Rounds the seconds of a non-negative sexagesimal value to `decimals`
/// decimals, carrying into minutes and whole units.
///
/// The seconds are taken in magnitude. A `NaN` or infinite value is
/// returned as is (there is nothing to carry); the counters saturate
/// instead of overflowing.
fn round_sexagesimal(whole: u32, min: u32, sec: f64, decimals: u32) -> (u32, u32, f64) {
    let scale = 10_f64.powi(decimals.min(9) as i32);
    let sec = sec.abs();
    let scaled = sec * scale;
    // Values too large to scale have no decimals left to round.
    let sec = if scaled.is_finite() {
        scaled.round() / scale
    } else {
        sec
    };
    if !sec.is_finite() {
        return (whole, min, sec);
    }

    // The remainder is exact and in [0, 60); the carry is a whole number.
    let rem = sec % 60.0;
    let carry = ((sec - rem) / 60.0).round();
    let total_min = u64::from(min).saturating_add(carry as u64);
    let whole = u64::from(whole).saturating_add(total_min / 60);

    (
        u32::try_from(whole).unwrap_or(u32::MAX),
        (total_min % 60) as u32,
        rem,
    )
}

/// Writes `sec` zero-padded to two integer digits, with the formatter's
/// precision if one is given.
fn write_sec(f: &mut fmt::Formatter<'_>, sec: f64) -> fmt::Result {
    match f.precision() {
        Some(p) => {
            let width = if p > 0 { p + 3 } else { 2 };
            write!(f, "{sec:0width$.p$}")
        }
        None if sec < 10.0 => write!(f, "0{sec}"),
        None => write!(f, "{sec}"),
    }
}

/// An angle in degrees, arcminutes and arcseconds, with an explicit sign
///
/// Unlike the tuple returned by [`dms_frm_deg`], this can represent
/// angles between 0° and −1° without ambiguity (Meeus, ch. 1).
///
/// Formatting with `{}` prints the stored components, e.g. `-0°30′15.5″`;
/// with a precision, e.g. `{:.1}`, the arcseconds are rounded to that many
/// decimals first, with carry.
///
/// # Example
///
/// ```
/// use astro::angle::Dms;
///
/// let a = Dms::new(true, 0, 32, 41.0); // −0°32′41″
/// assert!((a.to_deg() + (32.0 / 60.0 + 41.0 / 3600.0)).abs() < 1e-15);
/// assert_eq!(Dms::from_deg_rounded(a.to_deg(), 3), Dms::new(true, 0, 32, 41.0));
///
/// // 23.99999° is 23°59′59.964″: rounding to 0.1″ carries into the degrees
/// assert_eq!(format!("{:.1}", Dms::from_deg(23.999_99)), "24°00′00.0″");
/// // A negative angle that rounds to zero loses its sign
/// assert_eq!(format!("{:.1}", Dms::new(true, 0, 0, 0.01)), "0°00′00.0″");
/// assert_eq!(format!("{:.1}", Dms::new(true, 0, 30, 15.0)), "-0°30′15.0″");
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dms {
    /// `true` if the angle is negative
    pub negative: bool,
    /// Whole degrees (magnitude)
    pub deg: u32,
    /// Arcminutes (magnitude), normally `0..=59`
    pub min: u32,
    /// Arcseconds (magnitude), normally in `[0, 60)`
    pub sec: f64,
}

impl Dms {
    /// Creates an angle from its sign and the magnitudes of its components
    ///
    /// # Arguments
    ///
    /// * `negative`: `true` for a negative angle
    /// * `deg`: Degrees
    /// * `min`: Arcminutes
    /// * `sec`: Arcseconds (its own sign is ignored)
    pub const fn new(negative: bool, deg: u32, min: u32, sec: f64) -> Self {
        Self {
            negative,
            deg,
            min,
            sec: sec.abs(),
        }
    }

    /// Splits an angle in degrees with decimals
    ///
    /// The components are truncated: `min` is in `0..=59` and `sec` in
    /// `[0, 60)`. Angles of `2³²` degrees or more saturate `deg`. `NaN`
    /// gives a `NaN` `sec`.
    ///
    /// # Arguments
    ///
    /// * `deg`: Angle *| in degrees*
    pub fn from_deg(deg: f64) -> Self {
        let (d, m, s) = split_sexagesimal(deg.abs());
        Self::new(deg < 0.0, d as u32, m as u32, s)
    }

    /// Splits an angle in degrees with decimals, rounding the arcseconds
    ///
    /// The arcseconds are rounded to `decimals` decimals (at most 9) and
    /// carried, so the result never shows 60″. An angle that rounds to
    /// zero is not negative.
    ///
    /// # Arguments
    ///
    /// * `deg`: Angle *| in degrees*
    /// * `decimals`: Number of decimals to keep in the arcseconds
    pub fn from_deg_rounded(deg: f64, decimals: u32) -> Self {
        Self::from_deg(deg).rounded(decimals)
    }

    /// Returns this angle with its arcseconds rounded to `decimals`
    /// decimals (at most 9), carried into the arcminutes and degrees
    pub fn rounded(self, decimals: u32) -> Self {
        let (deg, min, sec) = round_sexagesimal(self.deg, self.min, self.sec, decimals);
        let negative = self.negative && (deg != 0 || min != 0 || sec != 0.0);
        Self::new(negative, deg, min, sec)
    }

    /// Returns the angle in degrees with decimals
    pub fn to_deg(&self) -> f64 {
        let v = f64::from(self.deg) + f64::from(self.min) / 60.0 + self.sec.abs() / 3600.0;
        if self.negative { -v } else { v }
    }
}

impl fmt::Display for Dms {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let a = match f.precision() {
            Some(p) => self.rounded(p as u32),
            None => Self::new(self.negative, self.deg, self.min, self.sec),
        };
        let sign = if a.negative { "-" } else { "" };
        write!(f, "{sign}{}°{:02}′", a.deg, a.min)?;
        write_sec(f, a.sec)?;
        write!(f, "″")
    }
}

/// An angle in hours, minutes and seconds of time, with an explicit sign
///
/// One hour is 15°. Formatting works as for [`Dms`], e.g. `14h15m39.7s`.
///
/// # Example
///
/// ```
/// use astro::angle::Hms;
///
/// // Right ascension of Arcturus (Meeus, example 17.a)
/// let asc = Hms::new(false, 14, 15, 39.7);
/// assert!((asc.to_deg() - 213.9154).abs() < 1e-4);
/// assert_eq!(format!("{:.1}", Hms::from_deg(asc.to_deg())), "14h15m39.7s");
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hms {
    /// `true` if the angle is negative
    pub negative: bool,
    /// Whole hours (magnitude)
    pub hour: u32,
    /// Minutes of time (magnitude), normally `0..=59`
    pub min: u32,
    /// Seconds of time (magnitude), normally in `[0, 60)`
    pub sec: f64,
}

impl Hms {
    /// Creates an angle from its sign and the magnitudes of its components
    ///
    /// # Arguments
    ///
    /// * `negative`: `true` for a negative angle
    /// * `hour`: Hours
    /// * `min`: Minutes of time
    /// * `sec`: Seconds of time (its own sign is ignored)
    pub const fn new(negative: bool, hour: u32, min: u32, sec: f64) -> Self {
        Self {
            negative,
            hour,
            min,
            sec: sec.abs(),
        }
    }

    /// Splits an angle given in degrees into hours, minutes and seconds
    ///
    /// Truncated components, as in [`Dms::from_deg`]. The angle is not
    /// reduced to `[0, 24ʰ)`.
    ///
    /// # Arguments
    ///
    /// * `deg`: Angle *| in degrees*
    pub fn from_deg(deg: f64) -> Self {
        let (h, m, s) = split_sexagesimal((deg / 15.0).abs());
        Self::new(deg < 0.0, h as u32, m as u32, s)
    }

    /// Splits an angle given in degrees, rounding the seconds to
    /// `decimals` decimals (at most 9) with carry
    ///
    /// # Arguments
    ///
    /// * `deg`: Angle *| in degrees*
    /// * `decimals`: Number of decimals to keep in the seconds
    pub fn from_deg_rounded(deg: f64, decimals: u32) -> Self {
        Self::from_deg(deg).rounded(decimals)
    }

    /// Returns this angle with its seconds rounded to `decimals` decimals
    /// (at most 9), carried into the minutes and hours
    pub fn rounded(self, decimals: u32) -> Self {
        let (hour, min, sec) = round_sexagesimal(self.hour, self.min, self.sec, decimals);
        let negative = self.negative && (hour != 0 || min != 0 || sec != 0.0);
        Self::new(negative, hour, min, sec)
    }

    /// Returns the angle in degrees with decimals
    pub fn to_deg(&self) -> f64 {
        let h = f64::from(self.hour) + f64::from(self.min) / 60.0 + self.sec.abs() / 3600.0;
        if self.negative { -15.0 * h } else { 15.0 * h }
    }
}

impl fmt::Display for Hms {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let a = match f.precision() {
            Some(p) => self.rounded(p as u32),
            None => Self::new(self.negative, self.hour, self.min, self.sec),
        };
        let sign = if a.negative { "-" } else { "" };
        write!(f, "{sign}{}h{:02}m", a.hour, a.min)?;
        write_sec(f, a.sec)?;
        write!(f, "s")
    }
}

/// Reduces `angl` to `[0, period)`, with a single rounding at most.
#[inline]
fn limit_to_period(angl: f64, period: f64) -> f64 {
    // `%` is the exact IEEE remainder with the sign of `angl`, so no
    // precision is lost even for huge angles.
    let r = angl % period;
    let r = if r < 0.0 { r + period } else { r };
    // A tiny negative `angl` gives `r + period == period` after rounding;
    // that value is congruent to 0. `+ 0.0` turns −0 into +0.
    if r >= period { 0.0 } else { r + 0.0 }
}

/// Reduces `angl` to `[-period/2, period/2)` exactly.
#[inline]
fn limit_to_half_period(angl: f64, period: f64) -> f64 {
    let half = 0.5 * period;
    let r = angl % period;
    // Both adjustments are exact (Sterbenz), so small angles keep all
    // their digits.
    let r = if r >= half {
        r - period
    } else if r < -half {
        r + period
    } else {
        r
    };
    r + 0.0
}

/// Computes the equivalent angle in the `[0, 360)` degree range
///
/// Correct for every finite input: the reduction uses the exact floating
/// point remainder, so the only rounding is that of adding 360 to a
/// negative remainder (at most half a unit in the last place of the
/// result). Before version 3 a tiny negative angle gave exactly
/// 360 (outside the range), angles beyond about 10¹⁷ degrees were reduced
/// wrongly, and angles beyond about 3·10²¹ degrees were not reduced at
/// all. `NaN` and infinities give `NaN`.
///
/// # Arguments
///
/// * `angl`: Angle *| in degrees*
#[inline]
pub fn limit_to_360(angl: f64) -> f64 {
    limit_to_period(angl, 360.0)
}

/// Computes the equivalent angle in the `[0, 2π)` radian range
///
/// The reduction is exact with respect to [`TWO_PI`] (the `f64` nearest to
/// 2π). See [`limit_to_360`] for the edge cases.
///
/// # Arguments
///
/// * `angl`: Angle *| in radians*
#[inline]
pub fn limit_to_two_PI(angl: f64) -> f64 {
    limit_to_period(angl, TWO_PI)
}

/// Computes the equivalent angle in the `[-180, 180)` degree range
///
/// Useful for differences of longitudes or right ascensions. Small angles
/// are returned unchanged, with all their digits.
///
/// # Arguments
///
/// * `angl`: Angle *| in degrees*
#[inline]
pub fn limit_to_pm180(angl: f64) -> f64 {
    limit_to_half_period(angl, 360.0)
}

/// Computes the equivalent angle in the `[-π, π)` radian range
///
/// Useful for differences of longitudes or right ascensions. Small angles
/// are returned unchanged, with all their digits.
///
/// # Arguments
///
/// * `angl`: Angle *| in radians*
#[inline]
pub fn limit_to_pm_pi(angl: f64) -> f64 {
    limit_to_half_period(angl, TWO_PI)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Meeus, example 17.b: Pauwels' x, y, z for Arcturus and Spica
    #[test]
    fn pauwels_xyz_example_17_b() {
        let (x, y, z) = pauwels_xyz(
            213.9154_f64.to_radians(),
            19.1825_f64.to_radians(),
            201.2983_f64.to_radians(),
            (-11.1614_f64).to_radians(),
        );
        assert!((x - -0.497404).abs() < 5e-7, "{x}");
        assert!((y - -0.214303).abs() < 5e-7, "{y}");
        assert!((z - 0.840633).abs() < 5e-7, "{z}");
        // tan d; the book's 0.644283 is computed from the rounded x, y, z
        // (which give 0.6442835), hence the wider tolerance.
        assert!((x.hypot(y) / z - 0.644283).abs() < 1e-6);
    }

    #[test]
    fn pauwels_x_matches_textbook_form() {
        for &(a1, d1, a2, d2) in &[
            (0.1, 0.2, 0.3, 0.4),
            (3.0, -1.2, -2.0, 0.7),
            (5.5, 1.5, 0.2, -1.5),
        ] {
            let (x, _, _) = pauwels_xyz(a1, d1, a2, d2);
            let book =
                f64::cos(d1) * f64::sin(d2) - f64::sin(d1) * f64::cos(d2) * f64::cos(a2 - a1);
            assert!((x - book).abs() < 1e-15);
        }
    }

    #[test]
    fn split_carries_near_sixty() {
        let (d, m, s) = split_sexagesimal(1.0 + 1.0 / 60.0);
        assert_eq!((d, m, s), (1.0, 1.0, 0.0));
        let (d, m, s) =
            split_sexagesimal(59.0 / 60.0 + 59.0 / 3600.0 + 0.999_999_999_999_9 / 3600.0);
        assert_eq!((d, m, s), (1.0, 0.0, 0.0));
    }
}

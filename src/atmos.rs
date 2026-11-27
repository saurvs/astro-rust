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

//! Corrections for atmospheric refraction
//!
//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 16. Refraction lifts
//! a body above its true ("airless") position: the **apparent** altitude
//! `h₀` (as observed) exceeds the **true** altitude `h` (as computed from
//! coordinates) by the refraction `R`, so `h = h₀ − R` and `h₀ = h + R`.
//!
//! # Units
//!
//! Every function takes and returns angles in **radians** (the formulas
//! themselves are stated in degrees, arcminutes and arcseconds; the
//! conversions are done internally). Pressure is in millibars (hPa);
//! temperature is in kelvins for [`refrac_by_temp`] and in degrees Celsius
//! for [`refrac_pressr_temp_factor`] and the `*_alt_frm_*` functions.
//!
//! # Which function?
//!
//! | Known | Formula | Function | Accuracy |
//! |---|---|---|---|
//! | apparent alt. > 15° | Smart, eq. 16.1 | [`refrac_frm_apprnt_alt_15`] | |
//! | true alt. > 15° | eq. 16.2 | [`refrac_frm_true_alt_15`] | |
//! | apparent alt. | Bennett, eq. 16.3 | [`refrac_frm_apprnt_alt`] | 0.07′ |
//! | apparent alt. | Bennett, refined | [`refrac_frm_apprnt_alt_refined`] | 0.015′ |
//! | true alt. | Sæmundsson, eq. 16.4 | [`refrac_frm_true_alt`] | 4″ from Bennett |
//!
//! The formulas assume a pressure of 1010 mbar and a temperature of 10 °C;
//! multiply `R` by [`refrac_pressr_temp_factor`] for other conditions.
//! [`true_alt_frm_apprnt_alt`] and [`apprnt_alt_frm_true_alt`] do all of
//! this in one call, with well-defined results at every altitude.
//!
//! The values are for yellow light. Near the horizon the real refraction
//! fluctuates by several arcminutes with the weather (Meeus quotes about
//! 0.3° as normal), so no formula can be trusted to better than that there.
//!
//! # Example
//!
//! The apparent flattening of the setting Sun (Meeus, example 16.a): the
//! lower limb is seen at 0°30′ and the true diameter is 32′.
//!
//! ```
//! use astro::atmos;
//!
//! let arcmin = |x: f64| (x / 60.0).to_radians();
//!
//! let lower_apprnt = arcmin(30.0);
//! let lower_true = lower_apprnt - atmos::refrac_frm_apprnt_alt(lower_apprnt);
//! let upper_true = lower_true + arcmin(32.0);
//! let upper_apprnt = upper_true + atmos::refrac_frm_true_alt(upper_true);
//!
//! let vert_diam = upper_apprnt - lower_apprnt;
//! assert!((vert_diam.to_degrees() * 60.0 - 27.864).abs() < 1e-3);
//! assert!((vert_diam / arcmin(32.0) - 0.871).abs() < 1e-3);
//! ```

use std::f64::consts::FRAC_PI_2;

/// Converts arcminutes to radians.
#[inline]
fn arcmin(x: f64) -> f64 {
    (x / 60.0).to_radians()
}

/// Converts arcseconds to radians.
#[inline]
fn arcsec(x: f64) -> f64 {
    (x / 3600.0).to_radians()
}

/// The constant to add to [`refrac_frm_apprnt_alt`] (Bennett's formula) so
/// that the refraction vanishes at the zenith: 0.001 3515′ *| in radians*
///
/// Without it eq. 16.3 gives −0.08″ at the zenith (Meeus, chapter 16). It
/// changes the result by less than 0.1″ at any altitude.
pub const BENNETT_ZENITH_CORR: f64 = 0.0013515 / 60.0 * std::f64::consts::PI / 180.0;

/// The constant to add to [`refrac_frm_true_alt`] (Sæmundsson's formula)
/// so that the refraction vanishes at the zenith: 0.001 9279′
/// *| in radians*
///
/// Without it eq. 16.4 gives −0.12″ at the zenith (Meeus, chapter 16).
pub const SAEMUNDSSON_ZENITH_CORR: f64 = 0.0019279 / 60.0 * std::f64::consts::PI / 180.0;

/// Computes the refraction term for apparent altitudes greater than 15
/// degrees
///
/// Smart's formula, `R = 58.294″ tan z − 0.0668″ tan³ z` with the
/// apparent zenith distance `z = 90° − h₀` (Meeus, eq. 16.1). Below 15° it
/// becomes inaccurate, and it is meaningless near the horizon.
///
/// Earlier versions of this crate used `tan(180° − h₀)` in place of
/// `tan(90° − h₀)` and cubed the coefficient `0.0668″` instead of the
/// tangent, returning nonsense.
///
/// # Returns
///
/// * `refrac_term`: The refraction term *| in radians*, that needs to be
///   subtracted from the apparent altitude to get the
///   true altitude
///
/// # Arguments
///
/// * `apprnt_alt`: Apparent altitude, above 15° *| in radians*
pub fn refrac_frm_apprnt_alt_15(apprnt_alt: f64) -> f64 {
    let t = (FRAC_PI_2 - apprnt_alt).tan();

    arcsec(58.294 * t - 0.0668 * t * t * t)
}

/// Computes the refraction term for true altitudes greater than 15
/// degrees
///
/// `R = 58.276″ tan z − 0.0824″ tan³ z` with the true zenith distance
/// `z = 90° − h` (Meeus, eq. 16.2, derived from eq. 16.1). Below 15° it
/// becomes inaccurate, and it is meaningless near the horizon.
///
/// Earlier versions of this crate had the same errors as in
/// [`refrac_frm_apprnt_alt_15`].
///
/// # Returns
///
/// * `refrac_term`: The refraction term *| in radians*, that needs to be
///   added to the true altitude to get the apparent
///   altitude
///
/// # Arguments
///
/// * `true_alt`: True altitude, above 15° *| in radians*
pub fn refrac_frm_true_alt_15(true_alt: f64) -> f64 {
    let t = (FRAC_PI_2 - true_alt).tan();

    arcsec(58.276 * t - 0.0824 * t * t * t)
}

/// Bennett's formula in its own units: apparent altitude in degrees,
/// refraction in arcminutes.
#[inline]
fn bennett_arcmin(apprnt_alt_deg: f64) -> f64 {
    1.0 / (apprnt_alt_deg + 7.31 / (apprnt_alt_deg + 4.4))
        .to_radians()
        .tan()
}

/// Sæmundsson's formula in its own units: true altitude in degrees,
/// refraction in arcminutes.
#[inline]
fn saemundsson_arcmin(true_alt_deg: f64) -> f64 {
    1.02 / (true_alt_deg + 10.3 / (true_alt_deg + 5.11))
        .to_radians()
        .tan()
}

/// Computes the refraction term from the apparent altitude, with
/// Bennett's formula
///
/// `R = 1′ / tan(h₀ + 7.31° / (h₀ + 4.4))`, with `h₀` in degrees (Meeus,
/// eq. 16.3), as published. It is accurate to 0.07′ (largest near 12°) for
/// all apparent altitudes from 0° to 90°. At the zenith it gives −0.08″
/// instead of zero; add [`BENNETT_ZENITH_CORR`] to remove that.
///
/// Negative apparent altitudes (seen from above sea level, below the
/// astronomical horizon) are evaluated by the same formula; it still
/// behaves smoothly down to about −1.7°, below which it is meaningless
/// (and it is singular at −4.4°). See [`true_alt_frm_apprnt_alt`] for a
/// function that is defined everywhere.
///
/// # Returns
///
/// * `refrac_term`: The refraction term *| in radians*, that needs to be
///   subtracted from the apparent altitude to get the
///   true altitude
///
/// # Arguments
///
/// * `apprnt_alt`: Apparent altitude *| in radians*
pub fn refrac_frm_apprnt_alt(apprnt_alt: f64) -> f64 {
    arcmin(bennett_arcmin(apprnt_alt.to_degrees()))
}

/// Computes the refraction term from the apparent altitude, with
/// Bennett's refined formula
///
/// Bennett's formula ([`refrac_frm_apprnt_alt`]) followed by his
/// correction `ΔR = −0.06′ sin(14.7° R + 13°)`, with `R` in arcminutes
/// (Meeus, chapter 16). The maximum error is then about 0.015′ (0.9″) from
/// 0° to 90°, except at the zenith itself, where it gives −0.89″.
///
/// # Returns
///
/// * `refrac_term`: The refraction term *| in radians*, that needs to be
///   subtracted from the apparent altitude to get the
///   true altitude
///
/// # Arguments
///
/// * `apprnt_alt`: Apparent altitude *| in radians*
pub fn refrac_frm_apprnt_alt_refined(apprnt_alt: f64) -> f64 {
    let r = bennett_arcmin(apprnt_alt.to_degrees());

    arcmin(r - 0.06 * (14.7 * r + 13.0).to_radians().sin())
}

/// Computes the refraction term from the true altitude, with
/// Sæmundsson's formula
///
/// `R = 1.02′ / tan(h + 10.3° / (h + 5.11))`, with `h` in degrees (Meeus,
/// eq. 16.4), as published. It and Bennett's formula
/// ([`refrac_frm_apprnt_alt`]) agree to 4″ or better. At the zenith it gives −0.12″
/// instead of zero; add [`SAEMUNDSSON_ZENITH_CORR`] to remove that.
///
/// Negative true altitudes are evaluated by the same formula; it still
/// behaves smoothly down to about −1.9°, below which it is meaningless
/// (and it is singular at −5.11°). See [`apprnt_alt_frm_true_alt`] for a
/// function that is defined everywhere.
///
/// # Returns
///
/// * `refrac_term`: The refraction term *| in radians*, that needs to be
///   added to the true altitude to get the apparent
///   altitude
///
/// # Arguments
///
/// * `true_alt`: True altitude *| in radians*
pub fn refrac_frm_true_alt(true_alt: f64) -> f64 {
    arcmin(saemundsson_arcmin(true_alt.to_degrees()))
}

/// Computes the refraction term modifier for local pressure
///
/// `P / 1010` (Meeus, chapter 16). Refraction grows with pressure.
///
/// # Returns
///
/// * `refrac_term_modifier`: The value that needs to be multiplied by
///   the refraction term to account for local
///   pressure
///
/// # Arguments
///
/// * `pressure`: Local pressure *| in millibars*
#[inline(always)]
pub fn refrac_by_pressr(pressure: f64) -> f64 {
    pressure / 1010.0
}

/// Computes the refraction term modifier for local temperature
///
/// `283 / T`, with `T` the absolute temperature. Meeus writes it as
/// `283 / (273 + t)` with `t` in degrees Celsius, so pass `t + 273.0` to
/// reproduce the book exactly (using 273.15 changes the result by 0.05 %).
/// Refraction decreases as the temperature rises.
///
/// # Returns
///
/// * `refrac_term_modifier`: The value that needs to be multiplied by
///   the refraction term to account for local
///   temperature
///
/// # Arguments
///
/// * `temp`: Local temperature *| in kelvins*
#[inline(always)]
pub fn refrac_by_temp(temp: f64) -> f64 {
    283.0 / temp
}

/// Computes the factor by which the refraction changes for the local
/// pressure and temperature
///
/// `(P / 1010) × 283 / (273 + t)` (Meeus, chapter 16); `1` for the
/// standard 1010 mbar and 10 °C assumed by all the refraction formulas of
/// this module. The correction is only approximate (it ignores, for
/// example, the dependence of refraction on wavelength).
///
/// # Returns
///
/// * `factor`: Multiply the refraction term by this
///
/// # Arguments
///
/// * `pressure`: Local atmospheric pressure at the observer *| in
///   millibars*
/// * `temp`: Local air temperature *| in degrees Celsius*
#[inline]
pub fn refrac_pressr_temp_factor(pressure: f64, temp: f64) -> f64 {
    pressure / 1010.0 * 283.0 / (273.0 + temp)
}

// Below these altitudes (degrees) the formulas are not used any more: the
// refraction is held at its value there, so that the results stay finite
// and monotonic (the formulas turn over at −1.70° and −1.90° and diverge
// at −4.4° and −5.11°).
const MIN_APPRNT_ALT_DEG: f64 = -1.0;
const MIN_TRUE_ALT_DEG: f64 = -1.0;

/// Computes the true (airless) altitude of a body from its observed,
/// apparent altitude
///
/// `h = h₀ − R`, with `R` from Bennett's formula (Meeus, eq. 16.3) plus
/// [`BENNETT_ZENITH_CORR`], scaled by [`refrac_pressr_temp_factor`]. It is
/// the inverse of [`apprnt_alt_frm_true_alt`] to within about 4″.
///
/// The result is defined for every input: below an apparent altitude of
/// −1° (possible only for an elevated observer) the refraction is held at
/// its value at −1°, since the formula is not meaningful there; at and
/// above 90° no refraction is applied.
///
/// # Returns
///
/// * `true_alt`: True altitude *| in radians*
///
/// # Arguments
///
/// * `apprnt_alt`: Apparent altitude *| in radians*
/// * `pressure`: Local atmospheric pressure *| in millibars* (1010 for the
///   standard atmosphere)
/// * `temp`: Local air temperature *| in degrees Celsius* (10 for the
///   standard atmosphere)
pub fn true_alt_frm_apprnt_alt(apprnt_alt: f64, pressure: f64, temp: f64) -> f64 {
    let h0 = apprnt_alt.to_degrees();
    if h0 >= 90.0 {
        return apprnt_alt;
    }
    let r = bennett_arcmin(h0.max(MIN_APPRNT_ALT_DEG)) + 0.0013515;

    apprnt_alt - arcmin(r) * refrac_pressr_temp_factor(pressure, temp)
}

/// Computes the apparent altitude of a body from its true (airless)
/// altitude
///
/// `h₀ = h + R`, with `R` from Sæmundsson's formula (Meeus, eq. 16.4) plus
/// [`SAEMUNDSSON_ZENITH_CORR`], scaled by [`refrac_pressr_temp_factor`].
/// It is the inverse of [`true_alt_frm_apprnt_alt`] to within about 4″.
///
/// The result is defined for every input: below a true altitude of −1°
/// the refraction is held at its value at −1°, since the formula is not
/// meaningful there; at and above 90° no refraction is applied.
///
/// # Returns
///
/// * `apprnt_alt`: Apparent altitude *| in radians*
///
/// # Arguments
///
/// * `true_alt`: True altitude *| in radians*
/// * `pressure`: Local atmospheric pressure *| in millibars* (1010 for the
///   standard atmosphere)
/// * `temp`: Local air temperature *| in degrees Celsius* (10 for the
///   standard atmosphere)
///
/// # Example
///
/// ```
/// use astro::atmos;
///
/// // A star 10° above the horizon on a cold, high-pressure night.
/// let h = 10_f64.to_radians();
/// let h0 = atmos::apprnt_alt_frm_true_alt(h, 1030.0, -10.0);
/// let r_arcmin = (h0 - h).to_degrees() * 60.0;
/// // 5.41′ in the standard atmosphere, times 1030/1010 × 283/263.
/// assert!((r_arcmin - 5.936).abs() < 1e-3);
///
/// // At the zenith there is no refraction.
/// let zenith = std::f64::consts::FRAC_PI_2;
/// assert!((atmos::apprnt_alt_frm_true_alt(zenith, 1010.0, 10.0) - zenith).abs() < 1e-12);
/// ```
pub fn apprnt_alt_frm_true_alt(true_alt: f64, pressure: f64, temp: f64) -> f64 {
    let h = true_alt.to_degrees();
    if h >= 90.0 {
        return true_alt;
    }
    let r = saemundsson_arcmin(h.max(MIN_TRUE_ALT_DEG)) + 0.0019279;

    true_alt + arcmin(r) * refrac_pressr_temp_factor(pressure, temp)
}

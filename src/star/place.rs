//! Apparent place of a star (Meeus, *Astronomical Algorithms*, 2nd ed.,
//! chapter 23)

use crate::aberr;
use crate::coords;
use crate::ecliptic;
use crate::nutation;
use crate::precess;
use std::f64::consts::{PI, TAU};

#[inline]
fn norm_two_pi(a: f64) -> f64 {
    let r = a.rem_euclid(TAU);
    if r >= TAU { 0.0 } else { r }
}

/// Wraps an angle difference to `(-π, π]`
#[inline]
fn wrap_pi(a: f64) -> f64 {
    let r = norm_two_pi(a);
    if r > PI { r - TAU } else { r }
}

/// Catalogue (mean) position of a star, with its proper motion
///
/// The position must be referred to the FK5 system (or a catalogue based
/// on it, such as the J2000.0 catalogues); the epoch of the position and
/// its equinox are assumed to be the same, as in star catalogues.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeanPlace {
    /// Right ascension at the catalogue epoch *| in radians*
    pub asc: f64,
    /// Declination at the catalogue epoch *| in radians*
    pub dec: f64,
    /// Annual proper motion in right ascension, dα/dt *| in radians per
    /// Julian year* (an angle of right ascension, not multiplied by
    /// `cos(dec)`; one second of time per year is 15″ per year)
    pub pm_asc: f64,
    /// Annual proper motion in declination *| in radians per Julian year*
    pub pm_dec: f64,
    /// Julian Ephemeris Day of the catalogue epoch and equinox, e.g.
    /// [`precess::J2000`](crate::precess::J2000)
    pub epoch: f64,
}

/// Method used for the annual aberration in [`apparent_place`]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AberrMethod {
    /// Ron–Vondrák expression (Meeus, eq. 23.4 and Table 23.A), applied
    /// in the J2000.0 frame before precession. The more accurate method.
    RonVondrak,
    /// Classical formulae (Meeus, eq. 23.3), applied to the mean place of
    /// the date, with a low accuracy longitude of the Sun.
    Classical,
}

/// The apparent place of a star and the intermediate quantities of its
/// reduction from the mean place
///
/// All angles in radians. Right ascensions are in `[0, 2π)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ApparentPlace {
    /// Right ascension at the date after proper motion, still referred to
    /// the catalogue equinox
    pub asc_pm: f64,
    /// Declination at the date after proper motion, still referred to the
    /// catalogue equinox
    pub dec_pm: f64,
    /// Annual aberration in right ascension. With
    /// [`AberrMethod::Classical`] it is computed at the mean place of the
    /// date; with [`AberrMethod::RonVondrak`] in the J2000.0 frame, before
    /// precession.
    pub aberr_asc: f64,
    /// Annual aberration in declination (see `aberr_asc`)
    pub aberr_dec: f64,
    /// Right ascension referred to the mean equator and equinox of the
    /// date. With [`AberrMethod::Classical`] this is the mean place
    /// (without aberration); with [`AberrMethod::RonVondrak`] it already
    /// includes the aberration.
    pub asc_mean: f64,
    /// Declination referred to the mean equator of the date (see
    /// `asc_mean`)
    pub dec_mean: f64,
    /// Nutation in longitude Δψ
    pub nut_in_long: f64,
    /// Nutation in obliquity Δε
    pub nut_in_oblq: f64,
    /// Mean obliquity of the ecliptic ε0 (IAU 1980)
    pub mn_oblq: f64,
    /// True obliquity of the ecliptic ε = ε0 + Δε
    pub tru_oblq: f64,
    /// Correction for nutation in right ascension
    pub nut_asc: f64,
    /// Correction for nutation in declination
    pub nut_dec: f64,
    /// Apparent right ascension, referred to the true equator and equinox
    /// of the date
    pub asc: f64,
    /// Apparent declination, referred to the true equator of the date
    pub dec: f64,
}

/// Computes the apparent place of a star
///
/// Full reduction of a catalogue position to the apparent place at a
/// given instant (Meeus, chapter 23), in this order:
///
/// 1. proper motion, assumed uniform in right ascension and declination
///    (chapter 21);
/// 2. precession with the IAU 1976 expressions (eqs. 21.2 and 21.4) and
///    annual aberration, either with the Ron–Vondrák expression in the
///    J2000.0 frame before precessing, or with the classical formulae
///    after precessing (see [`AberrMethod`]);
/// 3. nutation (IAU 1980, Table 22.A), applied as a rotation of the
///    coordinates (see
///    [`nutation::nutate_eq_coords`](crate::nutation::nutate_eq_coords))
///    so that stars close to the celestial poles are handled too.
///
/// As in Meeus, the annual parallax and the deflection of light by the
/// Sun are neglected (they are below 0.8″, and below 0.03″ more than 15°
/// from the Sun, respectively). The catalogue position must be FK5-based;
/// for FK4 positions use the building blocks in
/// [`precess`](crate::precess) and [`aberr`](crate::aberr) instead.
///
/// # Arguments
///
/// * `star`  : Catalogue position and proper motion of the star
/// * `JD`    : Julian Ephemeris Day of the observation
/// * `method`: Method for the annual aberration
///
/// # Example
///
/// Meeus, example 23.b: θ Persei on 2028 November 13.19 TD.
///
/// ```
/// use astro::{precess, star};
///
/// let theta_per = star::MeanPlace {
///     // 2h44m11.986s, +49°13'42.48" (J2000.0)
///     asc: (15.0 * (2.0 + 44.0 / 60.0 + 11.986 / 3600.0_f64)).to_radians(),
///     dec: (49.0 + 13.0 / 60.0 + 42.48 / 3600.0_f64).to_radians(),
///     pm_asc: (0.03425 * 15.0 / 3600.0_f64).to_radians(),
///     pm_dec: (-0.0895 / 3600.0_f64).to_radians(),
///     epoch: precess::J2000,
/// };
/// let place = star::apparent_place(&theta_per, 2462088.69, star::AberrMethod::RonVondrak);
///
/// assert!((place.asc.to_degrees() - 41.5599646).abs() < 2e-6);
/// assert!((place.dec.to_degrees() - 49.3520685).abs() < 2e-6);
/// ```
pub fn apparent_place(star: &MeanPlace, JD: f64, method: AberrMethod) -> ApparentPlace {
    // Proper motion, with the axes of reference fixed
    let years = (JD - star.epoch) / precess::JULIAN_YEAR;
    let asc_pm = norm_two_pi(star.asc + star.pm_asc * years);
    let dec_pm = star.dec + star.pm_dec * years;

    let (nut_in_long, nut_in_oblq) = nutation::nutation(JD);
    let mn_oblq = ecliptic::mn_oblq_IAU(JD);

    let (asc_mean, dec_mean, aberr_asc, aberr_dec, asc_ab, dec_ab) = match method {
        AberrMethod::Classical => {
            let (asc, dec) = precess::precess_eq_coords(asc_pm, dec_pm, star.epoch, JD);
            let sun = aberr::sun_geom_long_low_accuracy(JD);
            let (d_asc, d_dec) = aberr::annual_aberr_in_eq_coords(asc, dec, sun, mn_oblq, JD);

            (asc, dec, d_asc, d_dec, asc + d_asc, dec + d_dec)
        }
        AberrMethod::RonVondrak => {
            // The Earth's velocity is given in the J2000.0 frame
            let (asc0, dec0) =
                precess::precess_eq_coords(asc_pm, dec_pm, star.epoch, precess::J2000);
            let (d_asc, d_dec) = aberr::stell_aberr_in_eq_coords(
                &coords::EqPoint {
                    asc: asc0,
                    dec: dec0,
                },
                JD,
            );
            let (asc, dec) =
                precess::precess_eq_coords(asc0 + d_asc, dec0 + d_dec, precess::J2000, JD);

            (asc, dec, d_asc, d_dec, asc, dec)
        }
    };

    let (asc, dec) = nutation::nutate_eq_coords(asc_ab, dec_ab, nut_in_long, nut_in_oblq, mn_oblq);

    ApparentPlace {
        asc_pm,
        dec_pm,
        aberr_asc,
        aberr_dec,
        asc_mean,
        dec_mean,
        nut_in_long,
        nut_in_oblq,
        mn_oblq,
        tru_oblq: mn_oblq + nut_in_oblq,
        nut_asc: wrap_pi(asc - asc_ab),
        nut_dec: dec - dec_ab,
        asc,
        dec,
    }
}

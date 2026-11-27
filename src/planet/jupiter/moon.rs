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

//! Positions of the four Galilean satellites of Jupiter
//!
//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 44, gives two methods:
//!
//! * [`low_accuracy_positions`]: a few terms per satellite, good enough to
//!   identify the satellites in a telescope or to draw the familiar
//!   wavy-line diagrams. The satellites are assumed to move in Jupiter's
//!   equatorial plane, so `Y` can be off by up to a few tenths of a Jupiter
//!   radius. [`apprnt_rect_coords`] is the original entry point for this
//!   method.
//! * [`high_accuracy_positions`]: Lieske's E5 theory, with the periodic terms
//!   in longitude, latitude and radius vector, the exact rotation to the sky
//!   plane, and corrections for differential light-time and perspective.
//!   Use it for eclipses, transits, occultations and mutual phenomena. It
//!   also gives the positions as seen from the Sun, from which the satellite
//!   shadows and eclipses follow ([`MoonPosition::shadow`],
//!   [`MoonPosition::phenomena`]).
//!
//! Conjunctions with Jupiter and between two satellites can be searched with
//! [`conjunctions_with_jupiter`] and [`mutual_conjunctions`].
//!
//! # Coordinates
//!
//! Both methods return the apparent rectangular coordinates `X`, `Y` of a
//! satellite relative to the centre of Jupiter's disk, in units of Jupiter's
//! equatorial radius:
//!
//! * `X` is positive towards the west and negative towards the east, along
//!   Jupiter's equator;
//! * `Y` is positive towards the north and negative towards the south, along
//!   Jupiter's rotation axis;
//! * `Z` (high-accuracy method only) is negative when the satellite is
//!   nearer to the observer than Jupiter and positive when it is farther.
//!
//! "West" and "north" are directions on the sky. In a simple inverting
//! telescope used from the northern hemisphere, south is up and west is to
//! the left.
//!
//! # Example
//!
//! ```
//! use astro::planet::jupiter::moon::{self, Moon};
//!
//! // 1992 December 16 at 0h UT (JDE 2448972.50068); Meeus, example 44.b
//! let positions = moon::high_accuracy_positions(2448972.50068);
//! let io = positions.moon(Moon::Io).apparent;
//! assert!((io.x - -3.4502).abs() < 2e-4);
//! assert!((io.y - 0.2137).abs() < 2e-4);
//!
//! // The low-accuracy method agrees to about two decimals in X.
//! let rough = moon::low_accuracy_positions(2448972.50068);
//! assert!((rough.moon(Moon::Io).x - io.x).abs() < 0.02);
//! ```

use std::f64::consts::TAU;

use crate::planet;

/// A Galilean satellite of Jupiter
///
/// Meeus numbers the satellites I to IV: Io (I), Europa (II),
/// Ganymede (III) and Callisto (IV).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Moon {
    /// Io, satellite I
    Io,
    /// Europa, satellite II
    Europa,
    /// Ganymede, satellite III
    Ganymede,
    /// Callisto, satellite IV
    Callisto,
}

impl Moon {
    /// The four satellites in Meeus' order, I to IV
    pub const ALL: [Moon; 4] = [Moon::Io, Moon::Europa, Moon::Ganymede, Moon::Callisto];

    /// Meeus' number of the satellite, from 1 (Io) to 4 (Callisto)
    pub fn number(self) -> u8 {
        self.index() as u8 + 1
    }

    fn index(self) -> usize {
        match self {
            Moon::Io => 0,
            Moon::Europa => 1,
            Moon::Ganymede => 2,
            Moon::Callisto => 3,
        }
    }
}

/// Light-time for a unit distance, in days per AU (Meeus, eq. 33.3)
const LIGHT_TIME_PER_AU: f64 = 0.0057755183;

/// Jupiter's equatorial radius as a fraction of the astronomical unit,
/// expressed as the number of equatorial radii in one AU (Meeus, ch. 44)
const JUPITER_RADII_PER_AU: f64 = 2095.0;

/// Factor that turns Jupiter's flattened disk into a circle when it
/// multiplies `Y` (Meeus, ch. 44, "satellite phenomena"); it is the ratio
/// of the equatorial to the polar radius.
pub const DISK_Y_STRETCH: f64 = 1.071374;

#[inline]
fn norm(angle: f64) -> f64 {
    angle.rem_euclid(TAU)
}

// ---------------------------------------------------------------------------
// Low accuracy
// ---------------------------------------------------------------------------

/// Position of one satellite from the low-accuracy method
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LowAccuracyMoon {
    /// The angle `u` before the correction for the main periodic term
    /// *| in radians, in `[0, 2π)`*
    pub u_uncorrected: f64,
    /// The angle `u` of the satellite, counted from its inferior conjunction
    /// with Jupiter: 0 at inferior conjunction, π/2 at greatest western
    /// elongation, π at superior conjunction and 3π/2 at greatest eastern
    /// elongation *| in radians, in `[0, 2π)`*
    pub u: f64,
    /// Distance from the centre of Jupiter *| in Jupiter equatorial radii*
    pub r: f64,
    /// Apparent `X` coordinate, positive to the west *| in Jupiter
    /// equatorial radii*
    pub x: f64,
    /// Apparent `Y` coordinate, positive to the north *| in Jupiter
    /// equatorial radii*
    pub y: f64,
}

/// Result of the low-accuracy method, with the auxiliary quantities of
/// Meeus' "lower accuracy" method of chapter 43 that it uses
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LowAccuracyPositions {
    /// Days since JDE 2451545.0 (2000 January 1.5 TD) reduced by the
    /// light-time from Jupiter, `d − Δ/173` *| in days*
    pub retarded_days: f64,
    /// Phase angle `ψ` of Jupiter, the angle Sun–Jupiter–Earth, with the
    /// sign of `sin K` *| in radians, within about ±12°*
    pub phase_angle: f64,
    /// Jupiter's equation of the centre `B` *| in radians*
    pub jupiter_eq_of_center: f64,
    /// Jupiter–Earth distance `Δ` *| in AU*
    pub earth_dist: f64,
    /// Jupiter-centric declination of the Sun `D_S` *| in radians*
    pub sun_declination: f64,
    /// Jupiter-centric declination of the Earth `D_E` *| in radians*
    pub earth_declination: f64,
    /// The four satellites, indexed in Meeus' order (see [`Moon::ALL`])
    pub moons: [LowAccuracyMoon; 4],
}

impl LowAccuracyPositions {
    /// The position of one satellite
    pub fn moon(&self, moon: Moon) -> &LowAccuracyMoon {
        &self.moons[moon.index()]
    }
}

/// Computes the positions of the four Galilean satellites with Meeus'
/// low-accuracy method
///
/// The method keeps only the largest periodic term of each satellite and
/// neglects the inclinations of the orbits on Jupiter's equator, so `X` is
/// good to a few hundredths and `Y` to a few tenths of a Jupiter radius.
/// The positions of the Earth and Jupiter come from the short series of
/// Meeus' chapter 43 ("lower accuracy"), so no planetary theory is needed.
///
/// Meeus, *Astronomical Algorithms*, 2nd ed., ch. 44 (first method) and
/// ch. 43 (lower accuracy).
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day (TD)
///
/// # Example
///
/// ```
/// use astro::planet::jupiter::moon::{self, Moon};
///
/// // Meeus, example 44.a
/// let pos = moon::low_accuracy_positions(2448972.50068);
/// let europa = pos.moon(Moon::Europa);
/// assert!((europa.x - 7.44).abs() < 0.005);
/// assert!((europa.y - 0.25).abs() < 0.005);
/// ```
pub fn low_accuracy_positions(jde: f64) -> LowAccuracyPositions {
    let d = jde - 2451545.0;

    // Meeus, ch. 43, "lower accuracy"
    let V = (172.74 + 0.00111588 * d).to_radians();
    let M = (357.529 + 0.9856003 * d).to_radians();
    let N = (20.020 + 0.0830853 * d + 0.329 * V.sin()).to_radians();
    let J = (66.115 + 0.9025179 * d - 0.329 * V.sin()).to_radians();
    let A = (1.915 * M.sin() + 0.020 * (2.0 * M).sin()).to_radians();
    let B = (5.555 * N.sin() + 0.168 * (2.0 * N).sin()).to_radians();
    let K = J + A - B;
    let R = 1.00014 - 0.01671 * M.cos() - 0.00014 * (2.0 * M).cos();
    let r = 5.20872 - 0.25208 * N.cos() - 0.00611 * (2.0 * N).cos();
    let delta = (r * r + R * R - 2.0 * r * R * K.cos()).sqrt();
    let psi = (R * K.sin() / delta).clamp(-1.0, 1.0).asin();

    // Meeus, ch. 44, low accuracy
    let dd = d - delta / 173.0;
    let psi_minus_B = psi - B;

    let u1 = (163.8069 + 203.4058646 * dd).to_radians() + psi_minus_B;
    let u2 = (358.4140 + 101.2916335 * dd).to_radians() + psi_minus_B;
    let u3 = (5.7176 + 50.2345180 * dd).to_radians() + psi_minus_B;
    let u4 = (224.8092 + 21.4879800 * dd).to_radians() + psi_minus_B;

    let G = (331.18 + 50.310482 * dd).to_radians();
    let H = (87.45 + 21.569231 * dd).to_radians();

    let corrections = [
        (0.473 * (2.0 * (u1 - u2)).sin()).to_radians(),
        (1.065 * (2.0 * (u2 - u3)).sin()).to_radians(),
        (0.165 * G.sin()).to_radians(),
        (0.843 * H.sin()).to_radians(),
    ];

    // The radii use the uncorrected angles.
    let radii = [
        5.9057 - 0.0244 * (2.0 * (u1 - u2)).cos(),
        9.3966 - 0.0882 * (2.0 * (u2 - u3)).cos(),
        14.9883 - 0.0216 * G.cos(),
        26.3627 - 0.1939 * H.cos(),
    ];

    // Meeus, ch. 43: Jupiter-centric declinations of the Sun and the Earth
    let lambda = (34.35 + 0.083091 * d + 0.329 * V.sin()).to_radians() + B;
    let Ds = (3.12 * (lambda + 42.8_f64.to_radians()).sin()).to_radians();
    let De = Ds
        - (2.22 * psi.sin() * (lambda + 22.0_f64.to_radians()).cos()
            + 1.30 * (r - delta) / delta * (lambda - 100.5_f64.to_radians()).sin())
        .to_radians();

    let uncorrected = [u1, u2, u3, u4];
    let moons = std::array::from_fn(|k| {
        let u = uncorrected[k] + corrections[k];
        LowAccuracyMoon {
            u_uncorrected: norm(uncorrected[k]),
            u: norm(u),
            r: radii[k],
            x: radii[k] * u.sin(),
            y: -radii[k] * u.cos() * De.sin(),
        }
    });

    LowAccuracyPositions {
        retarded_days: dd,
        phase_angle: psi,
        jupiter_eq_of_center: B,
        earth_dist: delta,
        sun_declination: Ds,
        earth_declination: De,
        moons,
    }
}

/// Computes the apparent rectangular coordinates of a Galilean moon, with
/// the low-accuracy method
///
/// This is Meeus' first ("low accuracy") method of chapter 44, the same as
/// [`low_accuracy_positions`]. For the E5 theory use
/// [`high_accuracy_positions`].
///
/// # Returns
///
/// `(X, Y)`, the apparent position of the moon relative to the centre of
/// Jupiter's disk as seen from the Earth, in units of Jupiter's equatorial
/// radius. `X` is positive to the west and lies along Jupiter's equator;
/// `Y` is positive to the north and lies along Jupiter's rotation axis.
///
/// # Arguments
///
/// * `JD`  : Julian Ephemeris Day (TD)
/// * `moon`: The [Moon](./enum.Moon.html)
pub fn apprnt_rect_coords(JD: f64, moon: &Moon) -> (f64, f64) {
    let m = low_accuracy_positions(JD).moons[moon.index()];
    (m.x, m.y)
}

// ---------------------------------------------------------------------------
// High accuracy (Lieske's E5 theory)
// ---------------------------------------------------------------------------

/// Rectangular coordinates of a satellite relative to the centre of
/// Jupiter, as seen from some point (the Earth or the Sun)
///
/// All three coordinates are in units of Jupiter's equatorial radius.
/// `X` is positive towards the west (along Jupiter's equator), `Y` towards
/// the north (along Jupiter's rotation axis), and `Z` is negative when the
/// satellite is nearer to the viewpoint than Jupiter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SatelliteCoords {
    /// `X`, positive to the west *| in Jupiter equatorial radii*
    pub x: f64,
    /// `Y`, positive to the north *| in Jupiter equatorial radii*
    pub y: f64,
    /// `Z`, negative in front of Jupiter *| in Jupiter equatorial radii*
    pub z: f64,
}

impl SatelliteCoords {
    /// Distance of the point (`X`, `Y`) from the centre of Jupiter's disk,
    /// measured in units of the radius of the flattened disk in that
    /// direction
    ///
    /// `Y` is multiplied by [`DISK_Y_STRETCH`] to make the disk circular,
    /// so a value below 1 means the satellite is projected on the disk, and
    /// the satellite's centre touches the limb when the value is exactly 1.
    pub fn disk_distance(&self) -> f64 {
        self.x.hypot(DISK_Y_STRETCH * self.y)
    }

    /// Whether the centre of the satellite is projected on Jupiter's disk
    pub fn is_on_disk(&self) -> bool {
        self.disk_distance() < 1.0
    }
}

/// The satellite phenomena in progress at an instant
///
/// Each flag refers to the centre of the satellite (or of its shadow), so
/// the moments when they change are the instants of central contact with
/// Jupiter's limb. The disk is taken as an ellipse with the planet's
/// flattening; the shadow is treated as a cylinder, without penumbra.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Phenomena {
    /// The satellite is in front of Jupiter's disk, as seen from the Earth
    pub transit: bool,
    /// The satellite is behind Jupiter's disk, as seen from the Earth
    pub occultation: bool,
    /// The satellite is in Jupiter's shadow
    pub eclipse: bool,
    /// The satellite's shadow falls on Jupiter's disk
    pub shadow_transit: bool,
}

/// Position of one Galilean satellite from the high-accuracy method
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MoonPosition {
    /// Mean longitude `ℓ`, referred to the B1950.0 equinox *| in radians,
    /// in `[0, 2π)`*
    pub mean_longitude: f64,
    /// Sum `Σ` of the periodic terms in longitude *| in radians*
    pub longitude_perturbation: f64,
    /// True longitude `L = ℓ + Σ`, reduced to the mean equinox of the date
    /// by adding the precession `P` *| in radians, in `[0, 2π)`*
    pub longitude: f64,
    /// Latitude `B` with respect to Jupiter's equatorial plane *| in
    /// radians*
    pub latitude: f64,
    /// Distance `R` from the centre of Jupiter *| in Jupiter equatorial
    /// radii*
    pub radius: f64,
    /// Apparent coordinates as seen from the Earth, corrected for
    /// differential light-time and perspective
    pub apparent: SatelliteCoords,
    /// Coordinates as seen from the Sun, with the same corrections. They
    /// refer to the same instant as `apparent`, so they tell where the
    /// satellite's shadow is, or whether the satellite is eclipsed, as seen
    /// from the Earth at the requested time.
    pub shadow: SatelliteCoords,
}

impl MoonPosition {
    /// The phenomena (transit, occultation, eclipse, shadow transit) in
    /// progress for this satellite
    ///
    /// Meeus, ch. 44, "satellite phenomena": the satellite is in front of
    /// or behind the disk when its centre is projected on the flattened
    /// disk and `Z` is negative or positive; eclipses and shadow transits
    /// follow from the coordinates as seen from the Sun in the same way.
    pub fn phenomena(&self) -> Phenomena {
        let on_disk = self.apparent.is_on_disk();
        let shadow_on_disk = self.shadow.is_on_disk();
        Phenomena {
            transit: on_disk && self.apparent.z < 0.0,
            occultation: on_disk && self.apparent.z > 0.0,
            eclipse: shadow_on_disk && self.shadow.z > 0.0,
            shadow_transit: shadow_on_disk && self.shadow.z < 0.0,
        }
    }
}

/// Result of the high-accuracy method
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HighAccuracyPositions {
    /// Light-time `τ` from Jupiter to the Earth *| in days*
    pub light_time: f64,
    /// Distance `Δ` from the Earth to Jupiter, for the light-time
    /// corrected position of Jupiter *| in AU*
    pub earth_dist: f64,
    /// The four satellites, indexed in Meeus' order (see [`Moon::ALL`])
    pub moons: [MoonPosition; 4],
}

impl HighAccuracyPositions {
    /// The position of one satellite
    pub fn moon(&self, moon: Moon) -> &MoonPosition {
        &self.moons[moon.index()]
    }
}

/// Values of the E5 theory for the four satellites at one instant
struct E5 {
    /// Mean longitudes ℓ (B1950.0), radians
    mean_longitude: [f64; 4],
    /// Σ, radians
    sigma: [f64; 4],
    /// True longitudes L (B1950.0), radians
    longitude: [f64; 4],
    /// Latitudes B on Jupiter's equator, radians
    latitude: [f64; 4],
    /// Radius vectors, Jupiter equatorial radii
    radius: [f64; 4],
    /// Longitude ψ of the node of Jupiter's equator on the ecliptic
    /// (B1950.0), radians
    psi: f64,
}

/// Evaluates Lieske's E5 theory as given by Meeus
///
/// `t` is the time in days from 1976 August 10, 0h TD (JDE 2443000.5),
/// already reduced by the light-time.
fn e5_theory(t: f64) -> E5 {
    // Mean longitudes
    let l1 = (106.07719 + 203.488955790 * t).to_radians();
    let l2 = (175.73161 + 101.374724735 * t).to_radians();
    let l3 = (120.55883 + 50.317609207 * t).to_radians();
    let l4 = (84.44459 + 21.571071177 * t).to_radians();

    // Longitudes of the perijoves
    let p1 = (97.0881 + 0.16138586 * t).to_radians();
    let p2 = (154.8663 + 0.04726307 * t).to_radians();
    let p3 = (188.1840 + 0.00712734 * t).to_radians();
    let p4 = (335.2868 + 0.00184000 * t).to_radians();

    // Longitudes of the nodes on Jupiter's equatorial plane
    let w1 = (312.3346 - 0.13279386 * t).to_radians();
    let w2 = (100.4411 - 0.03263064 * t).to_radians();
    let w3 = (119.1942 - 0.00717703 * t).to_radians();
    let w4 = (322.6186 - 0.00175934 * t).to_radians();

    // Principal inequality in the longitude of Jupiter
    let gamma = 0.33033 * (163.679 + 0.0010512 * t).to_radians().sin()
        + 0.03439 * (34.486 - 0.0161731 * t).to_radians().sin();

    // Phase of free libration
    let phl = (199.6766 + 0.17379190 * t).to_radians();

    // Longitude of the node of Jupiter's equator on the ecliptic
    let psi = (316.5182 - 0.00000208 * t).to_radians();

    // Mean anomalies of Jupiter and Saturn
    let G = (30.23756 + 0.0830925701 * t + gamma).to_radians();
    let G1 = (31.97853 + 0.0334597339 * t).to_radians();

    // Longitude of the perihelion of Jupiter (constant in E5)
    let Pi = 13.469942_f64.to_radians();

    let c52 = 52.225_f64.to_radians();
    let c188 = 188.37_f64.to_radians();
    let c149 = 149.15_f64.to_radians();

    // Periodic terms in longitude, in degrees
    let s1 = 0.47259 * (2.0 * (l1 - l2)).sin() - 0.03478 * (p3 - p4).sin()
        + 0.01081 * (l2 - 2.0 * l3 + p3).sin()
        + 0.00738 * phl.sin()
        + 0.00713 * (l2 - 2.0 * l3 + p2).sin()
        - 0.00674 * (p1 + p3 - 2.0 * Pi - 2.0 * G).sin()
        + 0.00666 * (l2 - 2.0 * l3 + p4).sin()
        + 0.00445 * (l1 - p3).sin()
        - 0.00354 * (l1 - l2).sin()
        - 0.00317 * (2.0 * psi - 2.0 * Pi).sin()
        + 0.00265 * (l1 - p4).sin()
        - 0.00186 * G.sin()
        + 0.00162 * (p2 - p3).sin()
        + 0.00158 * (4.0 * (l1 - l2)).sin()
        - 0.00155 * (l1 - l3).sin()
        - 0.00138 * (psi + w3 - 2.0 * Pi - 2.0 * G).sin()
        - 0.00115 * (2.0 * (l1 - 2.0 * l2 + w2)).sin()
        + 0.00089 * (p2 - p4).sin()
        + 0.00085 * (l1 + p3 - 2.0 * Pi - 2.0 * G).sin()
        + 0.00083 * (w2 - w3).sin()
        + 0.00053 * (psi - w2).sin();

    let s2 = 1.06476 * (2.0 * (l2 - l3)).sin()
        + 0.04256 * (l1 - 2.0 * l2 + p3).sin()
        + 0.03581 * (l2 - p3).sin()
        + 0.02395 * (l1 - 2.0 * l2 + p4).sin()
        + 0.01984 * (l2 - p4).sin()
        - 0.01778 * phl.sin()
        + 0.01654 * (l2 - p2).sin()
        + 0.01334 * (l2 - 2.0 * l3 + p2).sin()
        + 0.01294 * (p3 - p4).sin()
        - 0.01142 * (l2 - l3).sin()
        - 0.01057 * G.sin()
        - 0.00775 * (2.0 * (psi - Pi)).sin()
        + 0.00524 * (2.0 * (l1 - l2)).sin()
        - 0.00460 * (l1 - l3).sin()
        + 0.00316 * (psi - 2.0 * G + w3 - 2.0 * Pi).sin()
        - 0.00203 * (p1 + p3 - 2.0 * Pi - 2.0 * G).sin()
        + 0.00146 * (psi - w3).sin()
        - 0.00145 * (2.0 * G).sin()
        + 0.00125 * (psi - w4).sin()
        - 0.00115 * (l1 - 2.0 * l3 + p3).sin()
        - 0.00094 * (2.0 * (l2 - w2)).sin()
        + 0.00086 * (2.0 * (l1 - 2.0 * l2 + w2)).sin()
        - 0.00086 * (5.0 * G1 - 2.0 * G + c52).sin()
        - 0.00078 * (l2 - l4).sin()
        - 0.00064 * (3.0 * l3 - 7.0 * l4 + 4.0 * p4).sin()
        + 0.00064 * (p1 - p4).sin()
        - 0.00063 * (l1 - 2.0 * l3 + p4).sin()
        + 0.00058 * (w3 - w4).sin()
        + 0.00056 * (2.0 * (psi - Pi - G)).sin()
        + 0.00056 * (2.0 * (l2 - l4)).sin()
        + 0.00055 * (2.0 * (l1 - l3)).sin()
        + 0.00052 * (3.0 * l3 - 7.0 * l4 + p3 + 3.0 * p4).sin()
        - 0.00043 * (l1 - p3).sin()
        + 0.00041 * (5.0 * (l2 - l3)).sin()
        + 0.00041 * (p4 - Pi).sin()
        + 0.00032 * (w2 - w3).sin()
        + 0.00032 * (2.0 * (l3 - G - Pi)).sin();

    let s3 = 0.16490 * (l3 - p3).sin() + 0.09081 * (l3 - p4).sin() - 0.06907 * (l2 - l3).sin()
        + 0.03784 * (p3 - p4).sin()
        + 0.01846 * (2.0 * (l3 - l4)).sin()
        - 0.01340 * G.sin()
        - 0.01014 * (2.0 * (psi - Pi)).sin()
        + 0.00704 * (l2 - 2.0 * l3 + p3).sin()
        - 0.00620 * (l2 - 2.0 * l3 + p2).sin()
        - 0.00541 * (l3 - l4).sin()
        + 0.00381 * (l2 - 2.0 * l3 + p4).sin()
        + 0.00235 * (psi - w3).sin()
        + 0.00198 * (psi - w4).sin()
        + 0.00176 * phl.sin()
        + 0.00130 * (3.0 * (l3 - l4)).sin()
        + 0.00125 * (l1 - l3).sin()
        - 0.00119 * (5.0 * G1 - 2.0 * G + c52).sin()
        + 0.00109 * (l1 - l2).sin()
        - 0.00100 * (3.0 * l3 - 7.0 * l4 + 4.0 * p4).sin()
        + 0.00091 * (w3 - w4).sin()
        + 0.00080 * (3.0 * l3 - 7.0 * l4 + p3 + 3.0 * p4).sin()
        - 0.00075 * (2.0 * l2 - 3.0 * l3 + p3).sin()
        + 0.00072 * (p1 + p3 - 2.0 * Pi - 2.0 * G).sin()
        + 0.00069 * (p4 - Pi).sin()
        - 0.00058 * (2.0 * l3 - 3.0 * l4 + p4).sin()
        - 0.00057 * (l3 - 2.0 * l4 + p4).sin()
        + 0.00056 * (l3 + p3 - 2.0 * Pi - 2.0 * G).sin()
        - 0.00052 * (l2 - 2.0 * l3 + p1).sin()
        - 0.00050 * (p2 - p3).sin()
        + 0.00048 * (l3 - 2.0 * l4 + p3).sin()
        - 0.00045 * (2.0 * l2 - 3.0 * l3 + p4).sin()
        - 0.00041 * (p2 - p4).sin()
        - 0.00038 * (2.0 * G).sin()
        - 0.00037 * (p3 - p4 + w3 - w4).sin()
        - 0.00032 * (3.0 * l3 - 7.0 * l4 + 2.0 * p3 + 2.0 * p4).sin()
        + 0.00030 * (4.0 * (l3 - l4)).sin()
        + 0.00029 * (l3 + p4 - 2.0 * Pi - 2.0 * G).sin()
        - 0.00028 * (w3 + psi - 2.0 * Pi - 2.0 * G).sin()
        + 0.00026 * (l3 - Pi - G).sin()
        + 0.00024 * (l2 - 3.0 * l3 + 2.0 * l4).sin()
        + 0.00021 * (2.0 * (l3 - Pi - G)).sin()
        - 0.00021 * (l3 - p2).sin()
        + 0.00017 * (2.0 * (l3 - p3)).sin();

    let s4 = 0.84287 * (l4 - p4).sin() + 0.03431 * (p4 - p3).sin()
        - 0.03305 * (2.0 * (psi - Pi)).sin()
        - 0.03211 * G.sin()
        - 0.01862 * (l4 - p3).sin()
        + 0.01186 * (psi - w4).sin()
        + 0.00623 * (l4 + p4 - 2.0 * G - 2.0 * Pi).sin()
        + 0.00387 * (2.0 * (l4 - p4)).sin()
        - 0.00284 * (5.0 * G1 - 2.0 * G + c52).sin()
        - 0.00234 * (2.0 * (psi - p4)).sin()
        - 0.00223 * (l3 - l4).sin()
        - 0.00208 * (l4 - Pi).sin()
        + 0.00178 * (psi + w4 - 2.0 * p4).sin()
        + 0.00134 * (p4 - Pi).sin()
        + 0.00125 * (2.0 * (l4 - G - Pi)).sin()
        - 0.00117 * (2.0 * G).sin()
        - 0.00112 * (2.0 * (l3 - l4)).sin()
        + 0.00107 * (3.0 * l3 - 7.0 * l4 + 4.0 * p4).sin()
        + 0.00102 * (l4 - G - Pi).sin()
        + 0.00096 * (2.0 * l4 - psi - w4).sin()
        + 0.00087 * (2.0 * (psi - w4)).sin()
        - 0.00085 * (3.0 * l3 - 7.0 * l4 + p3 + 3.0 * p4).sin()
        + 0.00085 * (l3 - 2.0 * l4 + p4).sin()
        - 0.00081 * (2.0 * (l4 - psi)).sin()
        + 0.00071 * (l4 + p4 - 2.0 * Pi - 3.0 * G).sin()
        + 0.00061 * (l1 - l4).sin()
        - 0.00056 * (psi - w3).sin()
        - 0.00054 * (l3 - 2.0 * l4 + p3).sin()
        + 0.00051 * (l2 - l4).sin()
        + 0.00042 * (2.0 * (psi - G - Pi)).sin()
        + 0.00039 * (2.0 * (p4 - w4)).sin()
        + 0.00036 * (psi + Pi - p4 - w4).sin()
        + 0.00035 * (2.0 * G1 - G + c188).sin()
        - 0.00035 * (l4 - p4 + 2.0 * Pi - 2.0 * psi).sin()
        - 0.00032 * (l4 + p4 - 2.0 * Pi - G).sin()
        + 0.00030 * (2.0 * G1 - 2.0 * G + c149).sin()
        + 0.00029 * (3.0 * l3 - 7.0 * l4 + 2.0 * p3 + 2.0 * p4).sin()
        + 0.00028 * (l4 - p4 + 2.0 * psi - 2.0 * Pi).sin()
        - 0.00028 * (2.0 * (l4 - w4)).sin()
        - 0.00027 * (p3 - p4 + w3 - w4).sin()
        - 0.00026 * (5.0 * G1 - 3.0 * G + c188).sin()
        + 0.00025 * (w4 - w3).sin()
        - 0.00025 * (l2 - 3.0 * l3 + 2.0 * l4).sin()
        - 0.00023 * (3.0 * (l3 - l4)).sin()
        + 0.00021 * (2.0 * l4 - 2.0 * Pi - 3.0 * G).sin()
        - 0.00021 * (2.0 * l3 - 3.0 * l4 + p4).sin()
        + 0.00019 * (l4 - p4 - G).sin()
        - 0.00019 * (2.0 * l4 - p3 - p4).sin()
        - 0.00018 * (l4 - p4 + G).sin()
        - 0.00016 * (l4 + p3 - 2.0 * Pi - 2.0 * G).sin();

    let sigma = [
        s1.to_radians(),
        s2.to_radians(),
        s3.to_radians(),
        s4.to_radians(),
    ];
    let (s1, s2, s3) = (sigma[0], sigma[1], sigma[2]);

    // True longitudes
    let L1 = l1 + s1;
    let L2 = l2 + s2;
    let L3 = l3 + s3;
    let L4 = l4 + sigma[3];

    // Periodic terms in the latitudes; their sums are the tangents of the
    // latitudes. The coefficients of Σ inside the arguments multiply Σ in
    // the same angular unit as the other angles.
    let tan_b1 =
        0.0006393 * (L1 - w1).sin() + 0.0001825 * (L1 - w2).sin() + 0.0000329 * (L1 - w3).sin()
            - 0.0000311 * (L1 - psi).sin()
            + 0.0000093 * (L1 - w4).sin()
            + 0.0000075 * (3.0 * L1 - 4.0 * l2 - 1.9927 * s1 + w2).sin()
            + 0.0000046 * (L1 + psi - 2.0 * Pi - 2.0 * G).sin();

    let tan_b2 = 0.0081004 * (L2 - w2).sin() + 0.0004512 * (L2 - w3).sin()
        - 0.0003284 * (L2 - psi).sin()
        + 0.0001160 * (L2 - w4).sin()
        + 0.0000272 * (l1 - 2.0 * l3 + 1.0146 * s2 + w2).sin()
        - 0.0000144 * (L2 - w1).sin()
        + 0.0000143 * (L2 + psi - 2.0 * Pi - 2.0 * G).sin()
        + 0.0000035 * (L2 - psi + G).sin()
        - 0.0000028 * (l1 - 2.0 * l3 + 1.0146 * s2 + w3).sin();

    let tan_b3 = 0.0032402 * (L3 - w3).sin() - 0.0016911 * (L3 - psi).sin()
        + 0.0006847 * (L3 - w4).sin()
        - 0.0002797 * (L3 - w2).sin()
        + 0.0000321 * (L3 + psi - 2.0 * Pi - 2.0 * G).sin()
        + 0.0000051 * (L3 - psi + G).sin()
        - 0.0000045 * (L3 - psi - G).sin()
        - 0.0000045 * (L3 + psi - 2.0 * Pi).sin()
        + 0.0000037 * (L3 + psi - 2.0 * Pi - 3.0 * G).sin()
        + 0.0000030 * (2.0 * l2 - 3.0 * L3 + 4.03 * s3 + w2).sin()
        - 0.0000021 * (2.0 * l2 - 3.0 * L3 + 4.03 * s3 + w3).sin();

    let tan_b4 = -0.0076579 * (L4 - psi).sin() + 0.0044134 * (L4 - w4).sin()
        - 0.0005112 * (L4 - w3).sin()
        + 0.0000773 * (L4 + psi - 2.0 * Pi - 2.0 * G).sin()
        + 0.0000104 * (L4 - psi + G).sin()
        - 0.0000102 * (L4 - psi - G).sin()
        + 0.0000088 * (L4 + psi - 2.0 * Pi - 3.0 * G).sin()
        - 0.0000038 * (L4 + psi - 2.0 * Pi - G).sin();

    // Periodic terms in the radius vectors
    let r1 = -0.0041339 * (2.0 * (l1 - l2)).cos()
        - 0.0000387 * (l1 - p3).cos()
        - 0.0000214 * (l1 - p4).cos()
        + 0.0000170 * (l1 - l2).cos()
        - 0.0000131 * (4.0 * (l1 - l2)).cos()
        + 0.0000106 * (l1 - l3).cos()
        - 0.0000066 * (l1 + p3 - 2.0 * Pi - 2.0 * G).cos();

    let r2 = 0.0093848 * (l1 - l2).cos()
        - 0.0003116 * (l2 - p3).cos()
        - 0.0001744 * (l2 - p4).cos()
        - 0.0001442 * (l2 - p2).cos()
        + 0.0000553 * (l2 - l3).cos()
        + 0.0000523 * (l1 - l3).cos()
        - 0.0000290 * (2.0 * (l1 - l2)).cos()
        + 0.0000164 * (2.0 * (l2 - w2)).cos()
        + 0.0000107 * (l1 - 2.0 * l3 + p3).cos()
        - 0.0000102 * (l2 - p1).cos()
        - 0.0000091 * (2.0 * (l1 - l3)).cos();

    let r3 = -0.0014388 * (l3 - p3).cos() - 0.0007919 * (l3 - p4).cos()
        + 0.0006342 * (l2 - l3).cos()
        - 0.0001761 * (2.0 * (l3 - l4)).cos()
        + 0.0000294 * (l3 - l4).cos()
        - 0.0000156 * (3.0 * (l3 - l4)).cos()
        + 0.0000156 * (l1 - l3).cos()
        - 0.0000153 * (l1 - l2).cos()
        + 0.0000070 * (2.0 * l2 - 3.0 * l3 + p3).cos()
        - 0.0000051 * (l3 + p3 - 2.0 * Pi - 2.0 * G).cos();

    let r4 =
        -0.0073546 * (l4 - p4).cos() + 0.0001621 * (l4 - p3).cos() + 0.0000974 * (l3 - l4).cos()
            - 0.0000543 * (l4 + p4 - 2.0 * Pi - 2.0 * G).cos()
            - 0.0000271 * (2.0 * (l4 - p4)).cos()
            + 0.0000182 * (l4 - Pi).cos()
            + 0.0000177 * (2.0 * (l3 - l4)).cos()
            - 0.0000167 * (2.0 * l4 - psi - w4).cos()
            + 0.0000167 * (psi - w4).cos()
            - 0.0000155 * (2.0 * (l4 - Pi - G)).cos()
            + 0.0000142 * (2.0 * (l4 - psi)).cos()
            + 0.0000105 * (l1 - l4).cos()
            + 0.0000092 * (l2 - l4).cos()
            - 0.0000089 * (l4 - Pi - G).cos()
            - 0.0000062 * (l4 + p4 - 2.0 * Pi - 3.0 * G).cos()
            + 0.0000048 * (2.0 * (l4 - w4)).cos();

    // Mean distances, in Jupiter equatorial radii
    const A: [f64; 4] = [5.90569, 9.39657, 14.98832, 26.36273];

    E5 {
        mean_longitude: [l1, l2, l3, l4],
        sigma,
        longitude: [L1, L2, L3, L4],
        latitude: [tan_b1.atan(), tan_b2.atan(), tan_b3.atan(), tan_b4.atan()],
        radius: [
            A[0] * (1.0 + r1),
            A[1] * (1.0 + r2),
            A[2] * (1.0 + r3),
            A[3] * (1.0 + r4),
        ],
        psi,
    }
}

/// Orientation of Jupiter's equator and orbit needed to rotate the
/// satellite coordinates onto the sky
#[derive(Clone, Copy)]
struct Frame {
    /// Inclination of Jupiter's equator on its orbit, `I`
    equator_incl: f64,
    /// Longitude of the node of Jupiter's equator on the ecliptic, of date
    psi: f64,
    /// Inclination of Jupiter's orbit on the ecliptic, of date
    orbit_incl: f64,
    /// Longitude of the ascending node of Jupiter's orbit, of date
    orbit_node: f64,
}

/// Rotates Jupiter-equatorial coordinates (x towards the node of the
/// equator on the ecliptic, z towards Jupiter's north pole) into the
/// frame of the viewer looking at Jupiter in the direction (lambda0,
/// beta0). Returns (A6, B6, C6) of Meeus, ch. 44.
fn rotate(p: [f64; 3], f: &Frame, lambda0: f64, beta0: f64) -> [f64; 3] {
    let (x, y, z) = (p[0], p[1], p[2]);
    let (sin_I, cos_I) = f.equator_incl.sin_cos();
    let phi = f.psi - f.orbit_node;
    let (sin_phi, cos_phi) = phi.sin_cos();
    let (sin_i, cos_i) = f.orbit_incl.sin_cos();
    let (sin_node, cos_node) = f.orbit_node.sin_cos();
    let (sin_l0, cos_l0) = lambda0.sin_cos();
    let (sin_b0, cos_b0) = beta0.sin_cos();

    // towards Jupiter's orbital plane
    let a1 = x;
    let b1 = y * cos_I - z * sin_I;
    let c1 = y * sin_I + z * cos_I;
    // towards the ascending node of Jupiter's orbit
    let a2 = a1 * cos_phi - b1 * sin_phi;
    let b2 = a1 * sin_phi + b1 * cos_phi;
    let c2 = c1;
    // towards the ecliptic
    let a3 = a2;
    let b3 = b2 * cos_i - c2 * sin_i;
    let c3 = b2 * sin_i + c2 * cos_i;
    // towards the vernal equinox
    let a4 = a3 * cos_node - b3 * sin_node;
    let b4 = a3 * sin_node + b3 * cos_node;
    let c4 = c3;
    // towards the viewer
    let a5 = a4 * sin_l0 - b4 * cos_l0;
    let b5 = a4 * cos_l0 + b4 * sin_l0;
    let c5 = c4;
    let a6 = a5;
    let b6 = c5 * sin_b0 + b5 * cos_b0;
    let c6 = c5 * cos_b0 - b5 * sin_b0;

    [a6, b6, c6]
}

/// Constants K of the differential light-time correction, one per
/// satellite (Meeus, ch. 44)
const LIGHT_TIME_K: [f64; 4] = [17295.0, 21819.0, 27558.0, 36548.0];

/// Projects the four satellites on the sky of a viewer that sees Jupiter
/// in the direction (lambda0, beta0) at the distance `dist` (AU), and
/// applies the corrections for differential light-time and perspective.
fn sky_coords(
    sats: &[[f64; 3]; 4],
    radius: &[f64; 4],
    frame: &Frame,
    lambda0: f64,
    beta0: f64,
    dist: f64,
) -> [SatelliteCoords; 4] {
    // A fictitious fifth satellite above Jupiter's north pole gives the
    // direction of the projected rotation axis.
    let pole = rotate([0.0, 0.0, 1.0], frame, lambda0, beta0);
    let D = pole[0].atan2(pole[2]);
    let (sin_D, cos_D) = D.sin_cos();

    std::array::from_fn(|k| {
        let [a6, b6, c6] = rotate(sats[k], frame, lambda0, beta0);
        let mut x = a6 * cos_D - c6 * sin_D;
        let mut y = a6 * sin_D + c6 * cos_D;
        let z = b6;

        // differential light-time
        let q = x / radius[k];
        x += z.abs() / LIGHT_TIME_K[k] * (1.0 - q * q).max(0.0).sqrt();

        // perspective
        let w = dist / (dist + z / JUPITER_RADII_PER_AU);
        x *= w;
        y *= w;

        SatelliteCoords { x, y, z }
    })
}

/// Unit-free rectangular coordinates from spherical ones
#[inline]
fn rect(long: f64, lat: f64, r: f64) -> [f64; 3] {
    let (sin_l, cos_l) = long.sin_cos();
    let (sin_b, cos_b) = lat.sin_cos();
    [r * cos_b * cos_l, r * cos_b * sin_l, r * sin_b]
}

/// Computes the positions of the four Galilean satellites with Meeus'
/// high-accuracy method, based on Lieske's E5 theory
///
/// The positions of the Earth and Jupiter come from the crate's VSOP87
/// theory (see [`planet::heliocent_coords`]); Jupiter is taken at the
/// instant `jde − τ`, `τ` being the light-time, found by iteration. The
/// satellites' coordinates are rotated to the sky plane as seen from the
/// Earth and, separately, as seen from the Sun, and corrected for
/// differential light-time and for the perspective effect.
///
/// The coordinates as seen from the Sun describe the satellite's shadow:
/// they are computed for the same instant as the Earth view, so the eclipse
/// and shadow-transit flags of [`MoonPosition::phenomena`] tell what an
/// observer on the Earth sees at `jde` (Meeus' rule of adding the
/// light-time to the true times of these events is built in).
///
/// Lieske's E5 theory was fitted to observations of the 20th century and
/// its accuracy degrades slowly for dates far from that period; Meeus gives
/// no explicit limit. The results stay finite and smooth for any date.
///
/// Meeus, *Astronomical Algorithms*, 2nd ed., ch. 44 (second method),
/// eqs. 44.1 to 44.3.
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day (TD)
///
/// # Example
///
/// ```
/// use astro::planet::jupiter::moon::{self, Moon};
///
/// // Meeus, example 44.b
/// let pos = moon::high_accuracy_positions(2448972.50068);
/// let callisto = pos.moon(Moon::Callisto);
/// assert!((callisto.apparent.x - 7.0720).abs() < 2e-4);
/// assert!((callisto.apparent.y - 1.0291).abs() < 2e-4);
/// assert!(!callisto.phenomena().transit);
/// ```
pub fn high_accuracy_positions(jde: f64) -> HighAccuracyPositions {
    let (L0, B0, R0) = planet::heliocent_coords(&planet::Planet::Earth, jde);
    let earth = rect(L0, B0, R0);

    // Light-time iteration, starting with Δ = 5 AU. Only the sines and
    // cosines of the planetary latitudes are used, so their normalisation
    // does not matter here.
    let mut tau = LIGHT_TIME_PER_AU * 5.0;
    let mut jup = (0.0, 0.0, 0.0);
    let mut geo = [0.0; 3];
    let mut delta = 0.0;
    for _ in 0..10 {
        let (l, b, r) = planet::heliocent_coords(&planet::Planet::Jupiter, jde - tau);
        let helio = rect(l, b, r);
        geo = [
            helio[0] - earth[0],
            helio[1] - earth[1],
            helio[2] - earth[2],
        ];
        delta = (geo[0] * geo[0] + geo[1] * geo[1] + geo[2] * geo[2]).sqrt();
        jup = (l, b, r);
        let new_tau = LIGHT_TIME_PER_AU * delta;
        let done = (new_tau - tau).abs() < 1e-10;
        tau = new_tau;
        if done {
            break;
        }
    }

    let lambda0 = geo[1].atan2(geo[0]);
    let beta0 = geo[2].atan2(geo[0].hypot(geo[1]));

    let t = jde - 2443000.5 - tau;
    let e5 = e5_theory(t);

    // Precession from B1950.0 to the date
    let T0 = (jde - 2433282.423) / 36525.0;
    let P = (1.3966626 * T0 + 0.0003088 * T0 * T0).to_radians();

    // Inclination of Jupiter's equator on its orbit; T in centuries from
    // 1900.0
    let T1900 = (jde - 2415020.0) / 36525.0;
    let equator_incl = (3.120262 + 0.0006 * T1900).to_radians();

    // Inclination and node of Jupiter's orbit, mean equinox of the date
    // (Meeus, table 31.A)
    let T = (jde - 2451545.0) / 36525.0;
    let orbit_incl =
        (1.303267 + T * (-0.0054965 + T * (0.00000466 - T * 0.000000002))).to_radians();
    let orbit_node =
        (100.464407 + T * (1.0209774 + T * (0.00040315 + T * 0.000000404))).to_radians();

    let psi = e5.psi + P;
    let frame = Frame {
        equator_incl,
        psi,
        orbit_incl,
        orbit_node,
    };

    let longitude: [f64; 4] = std::array::from_fn(|k| e5.longitude[k] + P);
    let sats: [[f64; 3]; 4] = std::array::from_fn(|k| {
        let (sin_b, cos_b) = e5.latitude[k].sin_cos();
        let (sin_u, cos_u) = (longitude[k] - psi).sin_cos();
        let r = e5.radius[k];
        [r * cos_u * cos_b, r * sin_u * cos_b, r * sin_b]
    });

    let apparent = sky_coords(&sats, &e5.radius, &frame, lambda0, beta0, delta);

    // As seen from the Sun: Meeus' equations 44.1 with R = 0, i.e. the
    // heliocentric direction and distance of Jupiter.
    let (l, b, r) = jup;
    let shadow = sky_coords(&sats, &e5.radius, &frame, l, b, r);

    let moons = std::array::from_fn(|k| MoonPosition {
        mean_longitude: norm(e5.mean_longitude[k]),
        longitude_perturbation: e5.sigma[k],
        longitude: norm(longitude[k]),
        latitude: e5.latitude[k],
        radius: e5.radius[k],
        apparent: apparent[k],
        shadow: shadow[k],
    });

    HighAccuracyPositions {
        light_time: tau,
        earth_dist: delta,
        moons,
    }
}

/// Computes the apparent rectangular coordinates of a Galilean moon with
/// the high-accuracy method
///
/// Shorthand for `high_accuracy_positions(jde).moon(moon).apparent`; see
/// [`high_accuracy_positions`].
///
/// # Arguments
///
/// * `jde` : Julian Ephemeris Day (TD)
/// * `moon`: The satellite
pub fn high_accuracy_rect_coords(jde: f64, moon: Moon) -> SatelliteCoords {
    high_accuracy_positions(jde).moon(moon).apparent
}

// ---------------------------------------------------------------------------
// Conjunctions
// ---------------------------------------------------------------------------

/// Kind of a conjunction of a satellite with Jupiter
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConjunctionKind {
    /// The satellite passes in front of Jupiter: `X` goes from negative to
    /// positive and `Z` is negative
    Inferior,
    /// The satellite passes behind Jupiter: `X` goes from positive to
    /// negative and `Z` is positive
    Superior,
}

/// A conjunction of a satellite with Jupiter
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JupiterConjunction {
    /// Instant of the conjunction, when the apparent `X` is zero *| as a
    /// Julian Ephemeris Day*
    pub jde: f64,
    /// Inferior or superior conjunction
    pub kind: ConjunctionKind,
    /// Apparent `Y` of the satellite at the conjunction *| in Jupiter
    /// equatorial radii*. If `|Y|` is less than the polar radius (about
    /// 0.933) the satellite transits the disk (inferior conjunction) or is
    /// occulted (superior conjunction).
    pub y: f64,
}

/// A conjunction between two satellites
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MutualConjunction {
    /// Instant of the conjunction, when the two apparent `X` coordinates
    /// are equal *| as a Julian Ephemeris Day*
    pub jde: f64,
    /// Common apparent `X` of the two satellites *| in Jupiter equatorial
    /// radii*
    pub x: f64,
    /// Difference `Y(first) − Y(second)` of the apparent `Y` coordinates,
    /// i.e. the separation of the two satellites in the direction of
    /// Jupiter's axis *| in Jupiter equatorial radii*
    pub y_separation: f64,
    /// `Z` of the first satellite minus `Z` of the second; negative when
    /// the first satellite is the nearer one *| in Jupiter equatorial
    /// radii*
    pub z_difference: f64,
}

/// Step of the search for sign changes, in days. Conjunctions closer to
/// each other than this may be missed (this only happens for grazing mutual
/// conjunctions, where `X1 − X2` barely reaches zero).
const SEARCH_STEP: f64 = 0.01;

/// Finds the zeros of `f` in `[start, end]`, scanning with a fixed step
/// and refining each sign change by bisection.
fn find_zeros<F>(start: f64, end: f64, mut f: F) -> Vec<f64>
where
    F: FnMut(f64) -> f64,
{
    let mut zeros = Vec::new();
    if !(start.is_finite() && end.is_finite()) || end < start {
        return zeros;
    }
    let steps = ((end - start) / SEARCH_STEP).ceil().max(1.0) as u64;
    let h = (end - start) / steps as f64;

    let mut a = start;
    let mut fa = f(a);
    for n in 1..=steps {
        let b = if n == steps {
            end
        } else {
            start + n as f64 * h
        };
        let fb = f(b);
        if fa == 0.0 {
            zeros.push(a);
        } else if fa * fb < 0.0 {
            // bisection down to about 0.01 second
            let (mut lo, mut hi, mut flo) = (a, b, fa);
            while hi - lo > 1e-7 {
                let mid = 0.5 * (lo + hi);
                let fm = f(mid);
                if fm == 0.0 {
                    lo = mid;
                    hi = mid;
                    break;
                }
                if flo * fm < 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                    flo = fm;
                }
            }
            zeros.push(0.5 * (lo + hi));
        }
        a = b;
        fa = fb;
    }
    if fa == 0.0 {
        zeros.push(a);
    }
    zeros
}

/// Finds the conjunctions of a satellite with Jupiter in a time interval
///
/// A satellite is in conjunction with Jupiter when its apparent `X` (from
/// the high-accuracy method) is zero: an inferior conjunction when `X`
/// changes from negative to positive (the satellite is in front of
/// Jupiter), a superior conjunction when `X` changes from positive to
/// negative (it is behind). Meeus, ch. 44, "conjunctions with Jupiter".
///
/// The search scans the interval in steps of 0.01 day and refines each
/// sign change of `X` to about 0.01 second. Its cost grows linearly with
/// the length of the interval.
///
/// # Returns
///
/// The conjunctions in `[jde_start, jde_end]`, in chronological order.
/// The result is empty if the interval is empty, reversed or not finite.
///
/// # Arguments
///
/// * `moon`     : The satellite
/// * `jde_start`: Start of the interval *| Julian Ephemeris Day*
/// * `jde_end`  : End of the interval *| Julian Ephemeris Day*
pub fn conjunctions_with_jupiter(
    moon: Moon,
    jde_start: f64,
    jde_end: f64,
) -> Vec<JupiterConjunction> {
    let x = |jde: f64| high_accuracy_rect_coords(jde, moon).x;
    find_zeros(jde_start, jde_end, x)
        .into_iter()
        .map(|jde| {
            let c = high_accuracy_rect_coords(jde, moon);
            // The direction of motion decides the kind; Z agrees with it
            // except in degenerate cases, so use X's derivative.
            let dt = 1e-3;
            let rising = x(jde + dt) > x(jde - dt);
            JupiterConjunction {
                jde,
                kind: if rising {
                    ConjunctionKind::Inferior
                } else {
                    ConjunctionKind::Superior
                },
                y: c.y,
            }
        })
        .collect()
}

/// Finds the conjunctions between two satellites in a time interval
///
/// Two satellites are in conjunction when their apparent `X` coordinates
/// (from the high-accuracy method) are equal; the difference of their `Y`
/// coordinates is then their separation. Meeus, ch. 44, "mutual
/// conjunctions". The conjunction is not observable if either satellite is
/// eclipsed or hidden by Jupiter at that instant; check
/// [`MoonPosition::phenomena`].
///
/// The search scans the interval in steps of 0.01 day and refines each
/// sign change of `X(first) − X(second)` to about 0.01 second, so two
/// conjunctions less than about a quarter of an hour apart (a grazing
/// approach) can be missed. Its cost grows linearly with the length of the
/// interval.
///
/// # Returns
///
/// The conjunctions in `[jde_start, jde_end]`, in chronological order.
/// The result is empty if the two satellites are the same, or if the
/// interval is empty, reversed or not finite.
///
/// # Arguments
///
/// * `first`    : One satellite
/// * `second`   : The other satellite
/// * `jde_start`: Start of the interval *| Julian Ephemeris Day*
/// * `jde_end`  : End of the interval *| Julian Ephemeris Day*
pub fn mutual_conjunctions(
    first: Moon,
    second: Moon,
    jde_start: f64,
    jde_end: f64,
) -> Vec<MutualConjunction> {
    if first == second {
        return Vec::new();
    }
    let diff = |jde: f64| {
        let p = high_accuracy_positions(jde);
        p.moon(first).apparent.x - p.moon(second).apparent.x
    };
    find_zeros(jde_start, jde_end, diff)
        .into_iter()
        .map(|jde| {
            let p = high_accuracy_positions(jde);
            let a = p.moon(first).apparent;
            let b = p.moon(second).apparent;
            MutualConjunction {
                jde,
                x: 0.5 * (a.x + b.x),
                y_separation: a.y - b.y,
                z_difference: a.z - b.z,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn find_zeros_of_sine() {
        let z = find_zeros(0.1, 10.0, |t: f64| (t * std::f64::consts::PI).sin());
        assert_eq!(z.len(), 9);
        for (k, t) in z.iter().enumerate() {
            assert!((t - (k + 1) as f64).abs() < 1e-6, "{t}");
        }
    }

    #[test]
    fn find_zeros_bad_intervals() {
        assert!(find_zeros(1.0, 0.0, |t| t - 0.5).is_empty());
        assert!(find_zeros(f64::NAN, 1.0, |t| t - 0.5).is_empty());
        assert!(find_zeros(0.0, f64::INFINITY, |t| t - 0.5).is_empty());
        // exact zero at an end point
        assert_eq!(find_zeros(0.0, 1.0, |t| t), vec![0.0]);
        assert_eq!(find_zeros(-1.0, 0.0, |t| t), vec![0.0]);
    }

    #[test]
    fn rotate_preserves_length() {
        let f = Frame {
            equator_incl: 0.05,
            psi: 5.5,
            orbit_incl: 0.02,
            orbit_node: 1.75,
        };
        let p = rotate([3.0, -4.0, 12.0], &f, 1.2, -0.03);
        let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        assert!((n - 13.0).abs() < 1e-12);
    }
}

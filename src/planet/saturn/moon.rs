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

//! Positions of the eight major satellites of Saturn
//!
//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 46, after G.
//! Dourneau's theory of the satellites Mimas (I), Enceladus (II), Tethys
//! (III), Dione (IV), Rhea (V), Titan (VI), Hyperion (VII) and Iapetus
//! (VIII). [`positions`] computes all eight at once together with the
//! intermediate quantities; [`apprnt_rect_coords`] returns the coordinates
//! of one satellite.
//!
//! # Coordinates
//!
//! `X` and `Y` are the apparent rectangular coordinates of a satellite
//! relative to the centre of Saturn's disk, as seen from the Earth, in units
//! of Saturn's equatorial radius:
//!
//! * `X` is positive towards the west and negative towards the east, along
//!   Saturn's equator (the major axis of the rings);
//! * `Y` is positive towards the north and negative towards the south,
//!   along Saturn's rotation axis;
//! * `Z` is negative when the satellite is nearer to the Earth than Saturn
//!   and positive when it is farther.
//!
//! # Example
//!
//! ```
//! use astro::planet::saturn::moon::{self, Moon};
//!
//! // 1999 September 18 at 0h UT (JDE 2451439.50074); Meeus, example 46.a
//! let p = moon::positions(2451439.50074);
//! let titan = p.moon(Moon::Titan);
//! assert!((titan.x - 14.568).abs() < 5e-4);
//! assert!((titan.y - 4.738).abs() < 5e-4);
//! assert!(titan.z < 0.0); // in front of Saturn
//! ```

use crate::angle;
use crate::planet;
use crate::precess;

/// A major satellite of Saturn, numbered I to VIII as in Meeus
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Moon {
    /// Mimas, satellite I
    Mimas,
    /// Enceladus, satellite II
    Enceladus,
    /// Tethys, satellite III
    Tethys,
    /// Dione, satellite IV
    Dione,
    /// Rhea, satellite V
    Rhea,
    /// Titan, satellite VI
    Titan,
    /// Hyperion, satellite VII
    Hyperion,
    /// Iapetus, satellite VIII
    Iapetus,
}

impl Moon {
    /// The eight satellites in Meeus' order, I to VIII
    pub const ALL: [Moon; 8] = [
        Moon::Mimas,
        Moon::Enceladus,
        Moon::Tethys,
        Moon::Dione,
        Moon::Rhea,
        Moon::Titan,
        Moon::Hyperion,
        Moon::Iapetus,
    ];

    /// Meeus' number of the satellite, from 1 (Mimas) to 8 (Iapetus)
    pub fn number(self) -> u8 {
        self.index() as u8 + 1
    }

    fn index(self) -> usize {
        match self {
            Moon::Mimas => 0,
            Moon::Enceladus => 1,
            Moon::Tethys => 2,
            Moon::Dione => 3,
            Moon::Rhea => 4,
            Moon::Titan => 5,
            Moon::Hyperion => 6,
            Moon::Iapetus => 7,
        }
    }
}

/// Position of one satellite of Saturn
///
/// The angles `longitude`, `inclination` and `node` describe the orbit with
/// respect to Saturn's equatorial plane and the B1950.0 equinox: longitudes
/// are counted from the equinox along the B1950.0 ecliptic to the ascending
/// node of Saturn's equator (longitude 168.8112°), then along Saturn's
/// equator, and, for `longitude`, from the satellite's node along its orbit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MoonPosition {
    /// Longitude `λ` of the satellite *| in radians, in `[0, 2π)`*
    pub longitude: f64,
    /// Inclination `γ` of the orbit on Saturn's equator *| in radians, in
    /// `[0, π/2]`*
    pub inclination: f64,
    /// Longitude `Ω` of the ascending node of the orbit on Saturn's
    /// equator *| in radians, in `[0, 2π)`*
    pub node: f64,
    /// Distance `r` from the centre of Saturn *| in Saturn equatorial radii*
    pub radius: f64,
    /// Apparent `X`, positive to the west *| in Saturn equatorial radii*
    pub x: f64,
    /// Apparent `Y`, positive to the north *| in Saturn equatorial radii*
    pub y: f64,
    /// `Z`, negative when the satellite is nearer to the Earth than Saturn
    /// *| in Saturn equatorial radii*
    pub z: f64,
}

/// Positions of the eight satellites, with the geometry of Saturn they
/// were computed from
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Positions {
    /// Light-time `τ` from Saturn to the Earth *| in days*
    pub light_time: f64,
    /// Distance `Δ` from the Earth to the light-time corrected position of
    /// Saturn *| in AU*
    pub earth_dist: f64,
    /// Saturn's geocentric longitude `λ0`, referred to the VSOP87 mean
    /// ecliptic and equinox of the date. It is corrected for light-time
    /// only: neither aberration, nor nutation, nor the small FK5 offset is
    /// applied (none of them changes the satellites' `X`, `Y`)
    /// *| in radians, in `[0, 2π)`*
    pub saturn_long: f64,
    /// Saturn's geocentric latitude `β0`, mean ecliptic of the date
    /// *| in radians*
    pub saturn_lat: f64,
    /// `λ0` reduced to the B1950.0 ecliptic and equinox, the frame of
    /// Dourneau's theory *| in radians, in `[0, 2π)`*
    pub saturn_long_b1950: f64,
    /// `β0` reduced to the B1950.0 ecliptic *| in radians*
    pub saturn_lat_b1950: f64,
    /// The eight satellites, indexed in Meeus' order (see [`Moon::ALL`])
    pub moons: [MoonPosition; 8],
}

impl Positions {
    /// The position of one satellite
    pub fn moon(&self, moon: Moon) -> &MoonPosition {
        &self.moons[moon.index()]
    }
}

/// Light-time for a unit distance, in days per AU (Meeus, eq. 33.3)
const LIGHT_TIME_PER_AU: f64 = 0.0057755183;

/// Number of Saturn equatorial radii in one AU (Meeus, ch. 46)
const SATURN_RADII_PER_AU: f64 = 2475.0;

/// JD of the epoch B1950.0
const B1950: f64 = 2433282.4235;

/// Constants `K` of the differential light-time correction (Meeus, ch. 46)
const LIGHT_TIME_K: [f64; 8] = [
    20947.0, 23715.0, 26382.0, 29876.0, 35313.0, 53800.0, 59222.0, 91820.0,
];

/// Reduces an angle to `[0, 2π)`. (`rem_euclid` alone can return 2π for a
/// tiny negative angle.)
#[inline]
fn norm(angl: f64) -> f64 {
    angle::limit_to_two_PI(angl)
}

/// The "times" and auxiliary angles of Meeus, ch. 46, for the instant `jde`
/// (already reduced by the light-time)
struct Args {
    t1: f64,
    t2: f64,
    t4: f64,
    t6: f64,
    t7: f64,
    t8: f64,
    t9: f64,
    t10: f64,
    t11: f64,
    W0: f64,
    W1: f64,
    W2: f64,
    W3: f64,
    W4: f64,
    W5: f64,
    W6: f64,
    W7: f64,
    W8: f64,
    e1: f64,
}

impl Args {
    fn new(jde: f64) -> Args {
        let t1 = jde - 2411093.0;
        let t2 = t1 / 365.25;
        let t3 = (jde - 2433282.423) / 365.25 + 1950.0;
        let t4 = jde - 2411368.0;
        let t5 = t4 / 365.25;
        let t6 = jde - 2415020.0;
        let t7 = t6 / 36525.0;
        let t8 = t6 / 365.25;
        let t9 = (jde - 2442000.5) / 365.25;
        let t10 = jde - 2409786.0;
        let t11 = t10 / 36525.0;

        Args {
            t1,
            t2,
            t4,
            t6,
            t7,
            t8,
            t9,
            t10,
            t11,
            W0: (5.095 * (t3 - 1866.39)).to_radians(),
            W1: (74.4 + 32.39 * t2).to_radians(),
            W2: (134.3 + 92.62 * t2).to_radians(),
            W3: (42.0 - 0.5118 * t5).to_radians(),
            W4: (276.59 + 0.5118 * t5).to_radians(),
            W5: (267.2635 + 1222.1136 * t7).to_radians(),
            W6: (175.4762 + 1221.5515 * t7).to_radians(),
            W7: (2.4891 + 0.002435 * t7).to_radians(),
            W8: (113.35 - 0.2597 * t7).to_radians(),
            e1: 0.05589 - 0.000346 * t7,
        }
    }
}

/// Inclination of Saturn's equator on the B1950.0 ecliptic and longitude
/// of its ascending node
fn saturn_equator() -> (f64, f64, f64, f64) {
    let (s1, c1) = 28.0817_f64.to_radians().sin_cos();
    let (s2, c2) = 168.8112_f64.to_radians().sin_cos();
    (s1, c1, s2, c2)
}

/// Orbit of a satellite in Saturn's equatorial frame: λ, γ, Ω, r of Meeus
#[derive(Debug, Clone, Copy, PartialEq)]
struct EquatorialOrbit {
    lambda: f64,
    gamma: f64,
    node: f64,
    r: f64,
}

/// Osculating elements of satellites V to VIII, referred to the B1950.0
/// ecliptic, as fed to Meeus' subroutine
#[derive(Debug, Clone, Copy, PartialEq)]
struct EclipticElements {
    /// eccentricity
    e: f64,
    /// semimajor axis, Saturn radii
    a: f64,
    /// longitude of the ascending node
    node: f64,
    /// inclination
    i: f64,
    /// mean longitude λ'
    lambda1: f64,
    /// longitude of the pericentre
    p: f64,
}

fn mimas(w: &Args) -> EquatorialOrbit {
    let L = (127.64 + 381.994497 * w.t1
        - 43.57 * w.W0.sin()
        - 0.720 * (3.0 * w.W0).sin()
        - 0.02144 * (5.0 * w.W0).sin())
    .to_radians();
    let p = (106.1 + 365.549 * w.t2).to_radians();
    let M = L - p;
    let C =
        (2.18287 * M.sin() + 0.025988 * (2.0 * M).sin() + 0.00043 * (3.0 * M).sin()).to_radians();
    EquatorialOrbit {
        lambda: L + C,
        gamma: 1.563_f64.to_radians(),
        node: (54.5 - 365.072 * w.t2).to_radians(),
        r: 3.06879 / (1.0 + 0.01905 * (M + C).cos()),
    }
}

fn enceladus(w: &Args) -> EquatorialOrbit {
    let L =
        (200.317 + 262.7319002 * w.t1 + 0.25667 * w.W1.sin() + 0.20883 * w.W2.sin()).to_radians();
    let p = (309.107 + 123.44121 * w.t2).to_radians();
    let M = L - p;
    let C = (0.55577 * M.sin() + 0.00168 * (2.0 * M).sin()).to_radians();
    EquatorialOrbit {
        lambda: L + C,
        gamma: 0.0262_f64.to_radians(),
        node: (348.0 - 151.95 * w.t2).to_radians(),
        r: 3.94118 / (1.0 + 0.00485 * (M + C).cos()),
    }
}

fn tethys(w: &Args) -> EquatorialOrbit {
    // circular orbit
    EquatorialOrbit {
        lambda: (285.306
            + 190.69791226 * w.t1
            + 2.063 * w.W0.sin()
            + 0.03409 * (3.0 * w.W0).sin()
            + 0.001015 * (5.0 * w.W0).sin())
        .to_radians(),
        gamma: 1.0976_f64.to_radians(),
        node: (111.33 - 72.2441 * w.t2).to_radians(),
        r: 4.880998,
    }
}

fn dione(w: &Args) -> EquatorialOrbit {
    let L =
        (254.712 + 131.53493193 * w.t1 - 0.0215 * w.W1.sin() - 0.01733 * w.W2.sin()).to_radians();
    let p = (174.8 + 30.820 * w.t2).to_radians();
    let M = L - p;
    let C = (0.24717 * M.sin() + 0.00033 * (2.0 * M).sin()).to_radians();
    EquatorialOrbit {
        lambda: L + C,
        gamma: 0.0139_f64.to_radians(),
        node: (232.0 - 30.27 * w.t2).to_radians(),
        r: 6.24871 / (1.0 + 0.002157 * (M + C).cos()),
    }
}

fn rhea(w: &Args) -> EclipticElements {
    let p1 = (342.7 + 10.057 * w.t2).to_radians();
    let a1 = 0.000265 * p1.sin() + 0.01 * w.W4.sin();
    let a2 = 0.000265 * p1.cos() + 0.01 * w.W4.cos();
    let N = (345.0 - 10.057 * w.t2).to_radians();
    EclipticElements {
        e: a1.hypot(a2),
        a: 8.725924,
        node: (168.8034 + 0.736936 * N.sin() + 0.041 * w.W3.sin()).to_radians(),
        i: (28.0362 + 0.346898 * N.cos() + 0.01930 * w.W3.cos()).to_radians(),
        lambda1: (359.244 + 79.69004720 * w.t1 + 0.086754 * N.sin()).to_radians(),
        p: a1.atan2(a2),
    }
}

// Note on example 46.a: with these formulas Titan's e, λ(6) and r(6) come
// out as 0.0293394, 270.5249° and 19.9266, against the book's 0.0293386,
// 270.5289° and 19.9254. The book's three values are reproduced together
// only if ϖ is about 0.135° smaller than the iteration below gives, which no
// reading of the printed formula allows (Rhea, which uses the same W4,
// matches the book exactly). The final X and Y of Titan agree with the book
// to its last decimal either way, so the printed formulas are kept.
fn titan(w: &Args) -> EclipticElements {
    let L = (261.1582 + 22.57697855 * w.t4 + 0.074025 * w.W3.sin()).to_radians();
    let i1 = (27.45141 + 0.295999 * w.W3.cos()).to_radians();
    let node1 = (168.66925 + 0.628808 * w.W3.sin()).to_radians();
    let a1 = w.W7.sin() * (node1 - w.W8).sin();
    let a2 = w.W7.cos() * i1.sin() - w.W7.sin() * i1.cos() * (node1 - w.W8).cos();
    let g0 = 102.8623_f64.to_radians();
    let psi = a1.atan2(a2);
    let s = a1.hypot(a2);

    // Successive approximations of the longitude of the pericentre;
    // Meeus says three are always enough.
    let mut g = w.W4 - node1 - psi;
    let mut peri = w.W4;
    for _ in 0..6 {
        peri = w.W4 + 0.37515_f64.to_radians() * ((2.0 * g).sin() - (2.0 * g0).sin());
        g = peri - node1 - psi;
    }

    let e1 = 0.029092 + 0.00019048 * ((2.0 * g).cos() - (2.0 * g0).cos());
    let q = 2.0 * (w.W5 - peri);
    let b1 = i1.sin() * (node1 - w.W8).sin();
    let b2 = w.W7.cos() * i1.sin() * (node1 - w.W8).cos() - w.W7.sin() * i1.cos();
    let theta = b1.atan2(b2) + w.W8;
    let u = 2.0 * w.W5 - 2.0 * theta + psi;
    let h = 0.9375 * e1 * e1 * q.sin() + 0.1875 * s * s * (2.0 * (w.W5 - theta)).sin();

    EclipticElements {
        e: e1 + 0.002778797 * e1 * q.cos(),
        a: 20.216193,
        node: node1 + 0.031843_f64.to_radians() * s * u.sin() / i1.sin(),
        i: i1 + 0.031843_f64.to_radians() * s * u.cos(),
        lambda1: L - 0.254744_f64.to_radians()
            * (w.e1 * w.W6.sin() + 0.75 * w.e1 * w.e1 * (2.0 * w.W6).sin() + h),
        p: peri + 0.159215_f64.to_radians() * q.sin(),
    }
}

fn hyperion(w: &Args) -> EclipticElements {
    let eta = (92.39 + 0.5621071 * w.t6).to_radians();
    let zeta = (148.19 - 19.18 * w.t8).to_radians();
    let theta = (184.8 - 35.41 * w.t9).to_radians();
    let theta1 = theta - 7.5_f64.to_radians();
    let a_s = (176.0 + 12.22 * w.t8).to_radians();
    let b_s = (8.0 + 24.44 * w.t8).to_radians();
    let c_s = b_s + 5.0_f64.to_radians();
    let peri = (69.898 - 18.67088 * w.t8).to_radians();
    let phi = 2.0 * (peri - w.W5);
    let chi = (94.9 - 2.292 * w.t8).to_radians();

    EclipticElements {
        e: 0.103458 - 0.004099 * eta.cos() - 0.000167 * (zeta + eta).cos()
            + 0.000235 * (zeta - eta).cos()
            + 0.02303 * zeta.cos()
            - 0.00212 * (2.0 * zeta).cos()
            + 0.000151 * (3.0 * zeta).cos()
            + 0.00013 * phi.cos(),
        a: 24.50601 - 0.08686 * eta.cos() - 0.00166 * (zeta + eta).cos()
            + 0.00175 * (zeta - eta).cos(),
        node: (168.6812 + 1.40136 * chi.cos() + 0.68599 * w.W3.sin() - 0.0392 * c_s.sin()
            + 0.0366 * theta1.sin())
        .to_radians(),
        i: (27.3347 + 0.643486 * chi.cos() + 0.315 * w.W3.cos() + 0.018 * theta.cos()
            - 0.018 * c_s.cos())
        .to_radians(),
        lambda1: (177.047
            + 16.91993829 * w.t6
            + 0.15648 * chi.sin()
            + 9.142 * eta.sin()
            + 0.007 * (2.0 * eta).sin()
            - 0.014 * (3.0 * eta).sin()
            + 0.2275 * (zeta + eta).sin()
            + 0.2112 * (zeta - eta).sin()
            - 0.26 * zeta.sin()
            - 0.0098 * (2.0 * zeta).sin()
            - 0.013 * a_s.sin()
            + 0.017 * b_s.sin()
            - 0.0303 * phi.sin())
        .to_radians(),
        p: peri
            + (0.15648 * chi.sin()
                - 0.4457 * eta.sin()
                - 0.2657 * (zeta + eta).sin()
                - 0.3573 * (zeta - eta).sin()
                - 12.872 * zeta.sin()
                + 1.668 * (2.0 * zeta).sin()
                - 0.2419 * (3.0 * zeta).sin()
                - 0.07 * phi.sin())
            .to_radians(),
    }
}

fn iapetus(w: &Args) -> EclipticElements {
    let L = (261.1582 + 22.57697855 * w.t4).to_radians();
    let peri1 = (91.796 + 0.562 * w.t7).to_radians();
    let psi = (4.367 - 0.195 * w.t7).to_radians();
    let theta = (146.819 - 3.198 * w.t7).to_radians();
    let phi = (60.470 + 1.521 * w.t7).to_radians();
    let PHI = (205.055 - 2.091 * w.t7).to_radians();
    let e1 = 0.028298 + 0.001156 * w.t11;
    let peri0 = (352.91 + 11.71 * w.t11).to_radians();
    let mu = (76.3852 + 4.53795125 * w.t10).to_radians();
    let t11 = w.t11;
    let i1 = (18.4602 + t11 * (-0.9518 + t11 * (-0.072 + 0.0054 * t11))).to_radians();
    let node1 = (143.198 + t11 * (-3.919 + t11 * (0.116 + 0.008 * t11))).to_radians();

    let l = mu - peri0;
    let g = peri0 - node1 - psi;
    let g1 = peri0 - node1 - phi;
    let ls = w.W5 - peri1;
    let gs = peri1 - theta;
    let lT = L - w.W4;
    let gT = w.W4 - PHI;
    let u1 = 2.0 * (l + g - ls - gs);
    let u2 = l + g1 - lT - gT;
    let u3 = l + 2.0 * (g - ls - gs);
    let u4 = lT + gT - g1;
    let u5 = 2.0 * (ls + gs);

    let dp = (0.08077 * (g1 - gT).sin() + 0.02139 * (u5 - 2.0 * g).sin() - 0.00676 * u3.sin()
        + 0.01380 * l.sin()
        + 0.01632 * (l + u2).sin()
        + 0.03547 * u4.sin())
    .to_radians();
    let dnode = (0.04204 * (u5 + psi).sin()
        + 0.00235 * (l + g1 + lT + gT + phi).sin()
        + 0.00358 * (u2 + phi).sin())
    .to_radians();

    EclipticElements {
        e: e1 - 0.0014097 * (g1 - gT).cos()
            + 0.0003733 * (u5 - 2.0 * g).cos()
            + 0.0001180 * u3.cos()
            + 0.0002408 * l.cos()
            + 0.0002849 * (l + u2).cos()
            + 0.0006190 * u4.cos(),
        a: 58.935028 + 0.004638 * u1.cos() + 0.058222 * u2.cos(),
        node: node1 + dnode / i1.sin(),
        i: i1
            + (0.04204 * (u5 + psi).cos()
                + 0.00235 * (l + g1 + lT + gT + phi).cos()
                + 0.00360 * (u2 + phi).cos())
            .to_radians(),
        // The sign of the last term, +0.00650 sin(u5 + ψ), is a deliberate
        // choice. The book prints "+", and "+" is used here, but its
        // example 46.a was evidently computed with "−": only that sign
        // reproduces the example's λ' = 97.7552° and X, Y = −48.759,
        // +4.136, while "+" gives λ' = 97.7579° and X, Y = −48.7604,
        // +4.1372 (a difference of 0.0015 Saturn radii, about 0.015″).
        // "+" was kept because
        // * it is what the formula prints;
        // * it is the geometrically consistent sign: 0.00650 is
        //   0.04204·tan(i'/2) for i' ≈ 17.6°, the (1 − cos i)·δΩ shift of
        //   a longitude counted along the reference plane to the node and
        //   then along the orbit, caused by the node term w' below; it must
        //   therefore have the sign of w', which the example's
        //   Ω(8) = 138.9121° confirms is "+";
        // * an independent implementation of the chapter (the Go package
        //   soniakeys/meeus, `saturnmoons`) also uses "+" and gets
        //   X, Y = −48.760, +4.137, as this code does.
        // If Dourneau's original theory turns out to have "−", flip the
        // sign here and tighten the Iapetus tolerances of the tests.
        lambda1: mu
            + (-0.04299 * u2.sin()
                - 0.00789 * u1.sin()
                - 0.06312 * ls.sin()
                - 0.00295 * (2.0 * ls).sin()
                - 0.02231 * u5.sin()
                + 0.00650 * (u5 + psi).sin())
            .to_radians(),
        p: peri0 + dp / e1,
    }
}

/// Meeus' subroutine for satellites V to VIII: solves the equation of the
/// centre and refers the orbit to Saturn's equator.
fn to_equatorial(el: &EclipticElements) -> EquatorialOrbit {
    let (s1, c1, _, _) = saturn_equator();
    let e = el.e;
    let M = el.lambda1 - el.p;
    // Equation of the centre, in radians, as a series in e to e⁵ (the
    // book's subroutine; Kepler's equation is not solved). For Hyperion,
    // whose e reaches 0.13, the truncation error is at most about 2″ in
    // the orbit, i.e. 3e-4 Saturn radii in X and Y; for the others it is
    // far smaller.
    let C = (2.0 * e - 0.25 * e.powi(3) + 0.0520833333 * e.powi(5)) * M.sin()
        + (1.25 * e * e - 0.458333333 * e.powi(4)) * (2.0 * M).sin()
        + (1.083333333 * e.powi(3) - 0.671875 * e.powi(5)) * (3.0 * M).sin()
        + 1.072917 * e.powi(4) * (4.0 * M).sin()
        + 1.142708 * e.powi(5) * (5.0 * M).sin();
    let r = el.a * (1.0 - e * e) / (1.0 + e * (M + C).cos());

    let g = el.node - 168.8112_f64.to_radians();
    let (sin_i, cos_i) = el.i.sin_cos();
    let a1 = sin_i * g.sin();
    let a2 = c1 * sin_i * g.cos() - s1 * cos_i;
    // `clamp` (unlike `min`) keeps a NaN, so bad input cannot turn into a
    // plausible-looking γ = 90°.
    let gamma = a1.hypot(a2).clamp(0.0, 1.0).asin();
    let u = a1.atan2(a2);
    let h = c1 * sin_i - s1 * cos_i * g.cos();
    let psi = (s1 * g.sin()).atan2(h);

    EquatorialOrbit {
        lambda: el.lambda1 + C + u - g - psi,
        gamma,
        node: 168.8112_f64.to_radians() + u,
        r,
    }
}

/// Rotates a vector from Saturn's equatorial frame (x towards the node of
/// the equator on the B1950.0 ecliptic) to the frame of an observer who
/// sees Saturn at B1950.0 longitude `lambda0` and latitude `beta0`.
/// Returns (A4, B4, C4) of Meeus, ch. 46.
fn rotate(p: [f64; 3], lambda0: f64, beta0: f64) -> [f64; 3] {
    let (s1, c1, s2, c2) = saturn_equator();
    let (sin_l0, cos_l0) = lambda0.sin_cos();
    let (sin_b0, cos_b0) = beta0.sin_cos();
    let (x, y, z) = (p[0], p[1], p[2]);

    // towards the ecliptic
    let a1 = x;
    let b1 = c1 * y - s1 * z;
    let c1_ = s1 * y + c1 * z;
    // towards the vernal equinox
    let a2 = c2 * a1 - s2 * b1;
    let b2 = s2 * a1 + c2 * b1;
    let c2_ = c1_;
    // towards the observer
    let a3 = a2 * sin_l0 - b2 * cos_l0;
    let b3 = a2 * cos_l0 + b2 * sin_l0;
    let c3 = c2_;
    let b4 = b3 * cos_b0 + c3 * sin_b0;
    let c4 = c3 * cos_b0 - b3 * sin_b0;

    [a3, b4, c4]
}

/// The equatorial orbits of the eight satellites at `jde` (already reduced
/// by the light-time)
fn orbits(jde: f64) -> [EquatorialOrbit; 8] {
    let w = Args::new(jde);
    [
        mimas(&w),
        enceladus(&w),
        tethys(&w),
        dione(&w),
        to_equatorial(&rhea(&w)),
        to_equatorial(&titan(&w)),
        to_equatorial(&hyperion(&w)),
        to_equatorial(&iapetus(&w)),
    ]
}

#[inline]
fn rect(long: f64, lat: f64, r: f64) -> [f64; 3] {
    let (sin_l, cos_l) = long.sin_cos();
    let (sin_b, cos_b) = lat.sin_cos();
    [r * cos_b * cos_l, r * cos_b * sin_l, r * sin_b]
}

/// Saturn's geocentric position for the light-time corrected instant:
/// (λ0, β0, Δ, τ), mean ecliptic and equinox of the date
fn saturn_geocentric(jde: f64) -> (f64, f64, f64, f64) {
    let (L0, B0, R0) = planet::heliocent_coords(&planet::Planet::Earth, jde);
    let earth = rect(L0, B0, R0);

    // Light-time iteration starting from Δ = 9 AU, as Meeus suggests.
    let mut tau = LIGHT_TIME_PER_AU * 9.0;
    let mut geo = [0.0; 3];
    let mut delta = 0.0;
    for _ in 0..10 {
        let (l, b, r) = planet::heliocent_coords(&planet::Planet::Saturn, jde - tau);
        let s = rect(l, b, r);
        geo = [s[0] - earth[0], s[1] - earth[1], s[2] - earth[2]];
        delta = (geo[0] * geo[0] + geo[1] * geo[1] + geo[2] * geo[2]).sqrt();
        let new_tau = LIGHT_TIME_PER_AU * delta;
        let done = (new_tau - tau).abs() < 1e-10;
        tau = new_tau;
        if done {
            break;
        }
    }

    let lambda0 = geo[1].atan2(geo[0]);
    let beta0 = geo[2].atan2(geo[0].hypot(geo[1]));
    (lambda0, beta0, delta, tau)
}

/// Computes the apparent positions of the eight major satellites of Saturn
///
/// Saturn's position comes from the crate's VSOP87 theory (see
/// [`planet::heliocent_coords`]) for the instant `jde − τ`, the light-time
/// `τ` being found by iteration; its direction is then reduced to the
/// B1950.0 frame of Dourneau's theory with [`precess::precess_ecl_coords`].
/// The satellites are computed for `jde − τ`, rotated to the sky plane and
/// corrected for differential light-time and for the perspective effect.
///
/// Dourneau's theory was built from observations of the 19th and 20th
/// centuries; its accuracy degrades slowly for dates far from that period.
///
/// Meeus, *Astronomical Algorithms*, 2nd ed., ch. 46, eqs. 46.1 to 46.3.
///
/// # Arguments
///
/// * `jde`: Julian Ephemeris Day (TD)
///
/// # Example
///
/// ```
/// use astro::planet::saturn::moon::{self, Moon};
///
/// // Meeus, example 46.a
/// let p = moon::positions(2451439.50074);
/// let hyperion = p.moon(Moon::Hyperion);
/// assert!((hyperion.x - -18.001).abs() < 5e-4);
/// assert!((hyperion.y - -5.328).abs() < 5e-4);
/// assert!((p.light_time - 0.04942).abs() < 1e-5);
/// ```
pub fn positions(jde: f64) -> Positions {
    let (lambda0, beta0, delta, tau) = saturn_geocentric(jde);
    let (lambda0_b, beta0_b) = precess::precess_ecl_coords(lambda0, beta0, jde, B1950);

    let orbits = orbits(jde - tau);

    // A fictitious ninth satellite above Saturn's north pole gives the
    // direction of the projected rotation axis.
    let pole = rotate([0.0, 0.0, 1.0], lambda0_b, beta0_b);
    let D = pole[0].atan2(pole[2]);
    let (sin_D, cos_D) = D.sin_cos();

    let moons = std::array::from_fn(|k| {
        let o = orbits[k];
        let u = o.lambda - o.node;
        let w = o.node - 168.8112_f64.to_radians();
        let (sin_u, cos_u) = u.sin_cos();
        let (sin_w, cos_w) = w.sin_cos();
        let (sin_g, cos_g) = o.gamma.sin_cos();
        let p = [
            o.r * (cos_u * cos_w - sin_u * cos_g * sin_w),
            o.r * (sin_u * cos_w * cos_g + cos_u * sin_w),
            o.r * sin_u * sin_g,
        ];

        let [a4, b4, c4] = rotate(p, lambda0_b, beta0_b);
        let mut x = a4 * cos_D - c4 * sin_D;
        let mut y = a4 * sin_D + c4 * cos_D;
        let z = b4;

        // differential light-time
        let q = x / o.r;
        x += z.abs() / LIGHT_TIME_K[k] * (1.0 - q * q).max(0.0).sqrt();

        // perspective
        let wp = delta / (delta + z / SATURN_RADII_PER_AU);
        x *= wp;
        y *= wp;

        MoonPosition {
            longitude: norm(o.lambda),
            inclination: o.gamma,
            node: norm(o.node),
            radius: o.r,
            x,
            y,
            z,
        }
    });

    Positions {
        light_time: tau,
        earth_dist: delta,
        saturn_long: norm(lambda0),
        saturn_lat: beta0,
        saturn_long_b1950: norm(lambda0_b),
        saturn_lat_b1950: beta0_b,
        moons,
    }
}

/// Computes the apparent rectangular coordinates of a moon of Saturn
///
/// The same as [`positions`] for a single satellite.
///
/// # Returns
///
/// `(X, Y, Z)`, the position of the moon relative to the centre of
/// Saturn's disk as seen from the Earth, in units of Saturn's equatorial
/// radius. `X` is positive to the west, along Saturn's equator; `Y` is
/// positive to the north, along Saturn's rotation axis; `Z` is negative if
/// the moon is nearer to the Earth than Saturn and positive if it is
/// farther.
///
/// # Arguments
///
/// * `JD`  : Julian Ephemeris Day (TD)
/// * `moon`: The [Moon](./enum.Moon.html)
pub fn apprnt_rect_coords(JD: f64, moon: &Moon) -> (f64, f64, f64) {
    let m = positions(JD).moons[moon.index()];
    (m.x, m.y, m.z)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "expected {b} ± {tol}, got {a}");
    }

    // Meeus, example 46.a: osculating elements of satellites V to VIII,
    // for JDE 2451439.50074 − 0.04942. The tolerances are half a unit of the
    // last printed decimal except for Titan's e and Iapetus' λ', where the
    // book's own run differs from its printed formulas (see the comments in
    // `titan` and `iapetus`).
    #[test]
    fn elements_example_46a() {
        let w = Args::new(2451439.45132);
        let expected = [
            (rhea(&w), 0.0102018, 5e-8, 49.7917, 5e-5, 28.2962, 168.2640),
            (
                titan(&w),
                0.0293386,
                1e-6,
                273.4387,
                5e-5,
                27.7333,
                168.5439,
            ),
            (
                hyperion(&w),
                0.1187225,
                5e-8,
                78.2068,
                5e-5,
                27.2076,
                167.5721,
            ),
            (
                iapetus(&w),
                0.0286422,
                5e-8,
                97.7552,
                3e-3,
                17.2486,
                138.9121,
            ),
        ];
        for (el, e, e_tol, lambda1, lambda1_tol, i, node) in expected {
            close(el.e, e, e_tol);
            close(
                el.lambda1.to_degrees().rem_euclid(360.0),
                lambda1,
                lambda1_tol,
            );
            close(el.i.to_degrees(), i, 5e-5);
            close(el.node.to_degrees().rem_euclid(360.0), node, 5e-5);
        }
        // the two discrepancies, pinned down
        close(titan(&w).e, 0.0293394, 5e-8);
        close(
            iapetus(&w).lambda1.to_degrees().rem_euclid(360.0),
            97.75785,
            5e-5,
        );
    }

    // Meeus, example 46.a: λ, γ, Ω and r of the eight satellites. Titan's
    // λ and r and Iapetus' λ carry the discrepancies described in `titan`
    // and `iapetus`; everything else matches to the last printed decimal.
    #[test]
    fn equatorial_orbits_example_46a() {
        let o = orbits(2451439.45132);
        let expected = [
            (320.0015, 1.5630, 47.7110, 3.1224),
            (300.4638, 0.0262, 123.2135, 3.9257),
            (347.9653, 1.0976, 51.0615, 4.8810),
            (102.6463, 0.0139, 128.2982, 6.2561),
            (50.9947, 0.3359, 118.2589, 8.7054),
            (270.5289, 0.3701, 8.4467, 19.9254),
            (91.9269, 1.0461, 21.5992, 24.1160),
            (103.1318, 15.5062, 22.3756, 58.9396),
        ];
        let lambda_tol = [5e-5, 5e-5, 5e-5, 5e-5, 5e-5, 4.5e-3, 5e-5, 3e-3];
        let r_tol = [5e-5, 5e-5, 5e-5, 5e-5, 5e-5, 1.5e-3, 5e-5, 1e-4];
        for (k, (lambda, gamma, node, r)) in expected.into_iter().enumerate() {
            close(
                o[k].lambda.to_degrees().rem_euclid(360.0),
                lambda,
                lambda_tol[k],
            );
            close(o[k].gamma.to_degrees(), gamma, 5e-5);
            close(o[k].node.to_degrees().rem_euclid(360.0), node, 5e-5);
            close(o[k].r, r, r_tol[k]);
        }
    }

    #[test]
    fn rotate_preserves_length() {
        let p = rotate([3.0, -4.0, 12.0], 0.8, -0.04);
        close(
            (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt(),
            13.0,
            1e-12,
        );
    }

    // `norm` must stay in [0, 2π): `rem_euclid` returns 2π for tiny
    // negative angles.
    #[test]
    fn norm_range() {
        assert_eq!(norm(-1e-20), 0.0);
        assert_eq!(norm(std::f64::consts::TAU), 0.0);
        close(norm(-0.5), std::f64::consts::TAU - 0.5, 1e-15);
    }

    // A NaN in the elements must not come out as a finite inclination
    // (`f64::min` returns the other operand when one is NaN).
    #[test]
    fn nan_elements_give_nan_inclination() {
        let el = EclipticElements {
            e: 0.01,
            a: 10.0,
            node: f64::NAN,
            i: 0.3,
            lambda1: 2.0,
            p: 0.5,
        };
        assert!(to_equatorial(&el).gamma.is_nan());
    }

    #[test]
    fn equation_of_centre_vanishes_for_circular_orbit() {
        let el = EclipticElements {
            e: 0.0,
            a: 10.0,
            node: 1.0,
            i: 0.3,
            lambda1: 2.0,
            p: 0.5,
        };
        let o = to_equatorial(&el);
        close(o.r, 10.0, 1e-12);
    }
}

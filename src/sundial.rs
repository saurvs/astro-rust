//! Planar sundials
//!
//! Hour lines of a flat sundial of any orientation and inclination, with a
//! straight stylus perpendicular to the dial face, after R. Sagot and
//! D. Savoie (Meeus, *Astronomical Algorithms*, 2nd ed., chapter 58).
//!
//! The dial is described by
//! * the latitude `φ` of the place (positive north),
//! * its *gnomonic declination* `D`: the azimuth of the stylus (the
//!   perpendicular to the plane), measured from the south towards the
//!   west, so `D = 0` for a dial facing due south, `π/2` facing west, `π`
//!   facing north and `3π/2` facing east,
//! * the zenith distance `z` of the stylus: 0 for a horizontal dial, `π/2`
//!   for a vertical one,
//! * the length `a` of the stylus.
//!
//! Shadow coordinates `(x, y)` are those of the tip of the stylus' shadow
//! in the plane of the dial, with the origin at the foot of the stylus, the
//! x-axis horizontal and positive to the right, and the y-axis along the
//! line of greatest slope, positive upwards. They are in the unit of `a`.
//! For a horizontal dial (where `D` is arbitrary) `D = 0` is used, which
//! puts the x-axis towards the east and the y-axis towards the north.
//!
//! The Sun's hour angle `H` is measured from the upper meridian transit
//! (true noon), positive westwards (afternoon): `H = -π/4` is 9h true solar
//! time, `H = π/12` is 13h. See [`hour_angle_frm_solar_time`]. To draw
//! lines of mean or zone time, add the equation of time to `H`.
//!
//! ```
//! use astro::sundial::Sundial;
//!
//! // Inclined dial at latitude 40° N, D = 70°, z = 50°, a = 1
//! // (Meeus, example 58.a)
//! let dial = Sundial::new(
//!     40_f64.to_radians(),
//!     70_f64.to_radians(),
//!     50_f64.to_radians(),
//!     1.0,
//! );
//!
//! // Shadow at 14h true solar time (H = 30°) at the June solstice
//! let p = dial.shadow_point(30_f64.to_radians(), 23.44_f64.to_radians()).unwrap();
//! assert!((p.x - -0.0390).abs() < 0.00005);
//! assert!((p.y - -0.3615).abs() < 0.00005);
//!
//! // Before about 8h52m the Sun shines on the back of the dial
//! assert!(dial.shadow_point(-60_f64.to_radians(), 23.44_f64.to_radians()).is_none());
//!
//! // The center of the dial, where the hour lines meet
//! let c = dial.polar_stylus().unwrap();
//! assert!((c.x0 - 3.3880).abs() < 0.00005 && (c.y0 - -3.1102).abs() < 0.00005);
//! ```

use std::f64::consts::{FRAC_PI_2, PI};

/// A planar sundial with a stylus perpendicular to its face
///
/// See the [module documentation](self) for the conventions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sundial {
    /// Geographic latitude `φ` of the place, positive north *| in radians*
    pub lat: f64,
    /// Gnomonic declination `D`: azimuth of the perpendicular to the dial,
    /// measured from the south towards the west *| in radians*
    pub gnomonic_dec: f64,
    /// Zenith distance `z` of the perpendicular to the dial (0 for a
    /// horizontal dial, π/2 for a vertical one) *| in radians*
    pub zenith_dist: f64,
    /// Length `a` of the stylus, in any unit; shadow coordinates come out
    /// in the same unit
    pub stylus_len: f64,
}

/// A point on the face of a sundial, in the unit of the stylus length
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShadowPoint {
    /// Horizontal coordinate, positive to the right
    pub x: f64,
    /// Coordinate along the line of greatest slope, positive upwards
    pub y: f64,
}

/// The polar stylus of a sundial: the line through the tip of the
/// perpendicular stylus parallel to the Earth's axis
///
/// It meets the dial at the *center* of the dial, the point where all hour
/// lines converge.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PolarStylus {
    /// x-coordinate of the center of the dial
    pub x0: f64,
    /// y-coordinate of the center of the dial
    pub y0: f64,
    /// Length `u` of the polar stylus from the center to the tip of the
    /// perpendicular stylus, in the unit of the stylus length
    pub length: f64,
    /// Angle `ψ` between the polar stylus and the dial plane
    /// *| in radians*, in `[0, π/2]`
    pub angle: f64,
}

impl Sundial {
    /// Creates a sundial of any orientation
    ///
    /// # Arguments
    ///
    /// * `lat`         : Latitude `φ`, positive north *| in radians*
    /// * `gnomonic_dec`: Gnomonic declination `D`, the azimuth of the
    ///   perpendicular to the dial from the south towards the west
    ///   *| in radians*
    /// * `zenith_dist` : Zenith distance `z` of the perpendicular to the
    ///   dial *| in radians*
    /// * `stylus_len`  : Length `a` of the stylus
    pub fn new(lat: f64, gnomonic_dec: f64, zenith_dist: f64, stylus_len: f64) -> Self {
        Sundial {
            lat,
            gnomonic_dec,
            zenith_dist,
            stylus_len,
        }
    }

    /// Creates a horizontal sundial (`z = 0`, with `D = 0`: x-axis towards
    /// the east, y-axis towards the north)
    pub fn horizontal(lat: f64, stylus_len: f64) -> Self {
        Self::new(lat, 0.0, 0.0, stylus_len)
    }

    /// Creates a vertical sundial (`z = π/2`) whose face looks towards
    /// azimuth `gnomonic_dec` (from the south towards the west); the
    /// y-axis points to the zenith
    pub fn vertical(lat: f64, gnomonic_dec: f64, stylus_len: f64) -> Self {
        Self::new(lat, gnomonic_dec, FRAC_PI_2, stylus_len)
    }

    /// Creates the northern face of an equatorial sundial (plane parallel
    /// to the equator: `z = π/2 - φ`, `D = π`). It is lit while the Sun's
    /// declination is positive. The 12h line runs down the line of
    /// greatest slope.
    pub fn equatorial_north(lat: f64, stylus_len: f64) -> Self {
        Self::new(lat, PI, FRAC_PI_2 - lat, stylus_len)
    }

    /// Creates the southern face of an equatorial sundial (`z = π/2 + φ`,
    /// `D = 0`). It is lit while the Sun's declination is negative.
    pub fn equatorial_south(lat: f64, stylus_len: f64) -> Self {
        Self::new(lat, 0.0, FRAC_PI_2 + lat, stylus_len)
    }

    /// The quantity `P = sin φ cos z - cos φ sin z cos D` of Meeus' chapter
    /// 58: the sine of the angle between the Earth's axis and the dial
    /// plane, with sign
    fn p(&self) -> f64 {
        let (sp, cp) = self.lat.sin_cos();
        let (sz, cz) = self.zenith_dist.sin_cos();
        sp * cz - cp * sz * self.gnomonic_dec.cos()
    }

    /// The quantity `Q` of Meeus' chapter 58, proportional to the sine of
    /// the Sun's altitude above the dial plane (`Q > 0` when the face is
    /// lit).
    fn q(&self, hour_angle: f64, sun_dec: f64) -> f64 {
        let (sp, cp) = self.lat.sin_cos();
        let (sz, cz) = self.zenith_dist.sin_cos();
        let (sd, cd) = self.gnomonic_dec.sin_cos();
        let (sh, ch) = hour_angle.sin_cos();

        sd * sz * sh + (cp * cz + sp * sz * cd) * ch + self.p() * sun_dec.tan()
    }

    /// Tells whether the Sun at hour angle `hour_angle` and declination
    /// `sun_dec` (both *| in radians*) shines on the face of the dial
    ///
    /// This ignores the horizon; [`shadow_point`](Self::shadow_point)
    /// also requires the Sun to be up.
    pub fn is_lit(&self, hour_angle: f64, sun_dec: f64) -> bool {
        self.q(hour_angle, sun_dec) > 0.0
    }

    /// Computes the position of the tip of the stylus' shadow
    ///
    /// Meeus, chapter 58: `x = a N_x / Q`, `y = a N_y / Q`.
    ///
    /// # Returns
    ///
    /// * `Some(point)` if the Sun is above the (geometric) horizon and
    ///   shines on the face of the dial
    /// * `None` otherwise, or for non-finite input
    ///
    /// # Arguments
    ///
    /// * `hour_angle`: Hour angle `H` of the Sun, from true noon, positive
    ///   westwards (afternoon) *| in radians*
    /// * `sun_dec`   : Declination `δ` of the Sun *| in radians*
    pub fn shadow_point(&self, hour_angle: f64, sun_dec: f64) -> Option<ShadowPoint> {
        let (sp, cp) = self.lat.sin_cos();
        let (sz, cz) = self.zenith_dist.sin_cos();
        let (sd, cd) = self.gnomonic_dec.sin_cos();
        let (sh, ch) = hour_angle.sin_cos();
        let (sdec, cdec) = sun_dec.sin_cos();
        let td = sun_dec.tan();

        // Sun below the horizon (sin of altitude <= 0), or the face not lit
        // (Q <= 0); NaN input fails both tests.
        let sin_alt = sp * sdec + cp * cdec * ch;
        let q = self.q(hour_angle, sun_dec);
        if sin_alt.is_nan() || sin_alt <= 0.0 || q.is_nan() || q <= 0.0 {
            return None;
        }

        let nx = cd * sh - sd * (sp * ch - cp * td);
        let ny = cz * sd * sh - (cp * sz - sp * cz * cd) * ch - (sp * sz + cp * cz * cd) * td;

        let a = self.stylus_len;
        let point = ShadowPoint {
            x: a * nx / q,
            y: a * ny / q,
        };

        (point.x.is_finite() && point.y.is_finite()).then_some(point)
    }

    /// Computes the center of the dial and the polar stylus
    ///
    /// Meeus, chapter 58: `x0 = (a/P) cos φ sin D`,
    /// `y0 = -(a/P)(sin φ sin z + cos φ cos z cos D)`, `u = a / |P|`,
    /// `sin ψ = |P|`.
    ///
    /// # Returns
    ///
    /// `None` when the Earth's axis is parallel to the dial plane
    /// (`P = 0`, i.e. `cos D tan z = tan φ`, as for a polar dial): the hour
    /// lines are then parallel and have no center.
    pub fn polar_stylus(&self) -> Option<PolarStylus> {
        let p = self.p();
        if p.abs() < 1e-12 || !p.is_finite() {
            return None;
        }

        let (sp, cp) = self.lat.sin_cos();
        let (sz, cz) = self.zenith_dist.sin_cos();
        let a = self.stylus_len;

        Some(PolarStylus {
            x0: a / p * cp * self.gnomonic_dec.sin(),
            y0: -a / p * (sp * sz + cp * cz * self.gnomonic_dec.cos()),
            length: a / p.abs(),
            angle: p.abs().min(1.0).asin(),
        })
    }

    /// Computes the points of an hour line on the seven dates when the
    /// Sun's longitude is a multiple of 30° from 270° to 90°
    ///
    /// These are the traditional declination curves (solstices, equinoxes
    /// and the dates in between): with `ε` the obliquity, the declinations
    /// `δ = asin(sin ε sin λ)` for `λ = 270°, 300°, 330°, 0°, 30°, 60°,
    /// 90°`, which for `ε = 23.44°` are -23.44°, -20.15°, -11.47°, 0°,
    /// +11.47°, +20.15° and +23.44°. Points where the dial is in shadow or
    /// the Sun is down are `None`.
    ///
    /// # Arguments
    ///
    /// * `hour_angle`: Hour angle `H` of the Sun *| in radians*
    /// * `oblq_eclip`: Obliquity of the ecliptic *| in radians*
    pub fn hour_line(&self, hour_angle: f64, oblq_eclip: f64) -> [Option<ShadowPoint>; 7] {
        let mut points = [None; 7];
        for (i, point) in points.iter_mut().enumerate() {
            let sun_long = (270.0 + 30.0 * i as f64).to_radians();
            let dec = (oblq_eclip.sin() * sun_long.sin()).asin();
            *point = self.shadow_point(hour_angle, dec);
        }
        points
    }
}

/// Computes the Sun's hour angle at geometric sunset
///
/// `cos H0 = -tan φ tan δ` (Meeus, chapter 58). Sunrise is at `-H0`. The
/// Sun's center is taken on the mathematical horizon, without refraction
/// or semidiameter.
///
/// # Returns
///
/// * `H0` *| in radians*, in `[0, π]`: 0 if the Sun stays below the
///   horizon all day, π if it stays above it.
///
/// # Arguments
///
/// * `lat`    : Latitude of the place *| in radians*
/// * `sun_dec`: Declination of the Sun *| in radians*
pub fn sunset_hour_angle(lat: f64, sun_dec: f64) -> f64 {
    (-lat.tan() * sun_dec.tan()).clamp(-1.0, 1.0).acos()
}

/// Converts true (apparent) solar time to the Sun's hour angle
///
/// `H = (t - 12h) × 15°`: noon is 0, mornings are negative.
///
/// # Returns
///
/// * Hour angle *| in radians*
///
/// # Arguments
///
/// * `hours`: True solar time, in hours (eg: 9.5 for 9h30m)
#[inline]
pub fn hour_angle_frm_solar_time(hours: f64) -> f64 {
    ((hours - 12.0) * 15.0).to_radians()
}

//! Figure of the Earth adopted by the IAU in 1976
//!
//! These are the values used throughout Meeus, *Astronomical Algorithms*,
//! 2nd ed. (chapter 11 and Appendix I) for parallax, eclipses and
//! semidiameters.

/// Equatorial radius of the Earth (semimajor axis) *| in meters*
pub const EQUATORIAL_RADIUS: f64 = 6378140.0;

/// Flattening of the Earth
pub const FLATTENING: f64 = 1.0 / 298.257;

/// Polar radius of the Earth (semiminor axis) *| in meters*
/// (≈ 6356755 m)
pub const POLAR_RADIUS: f64 = EQUATORIAL_RADIUS * (1.0 - FLATTENING);

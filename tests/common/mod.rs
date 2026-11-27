//! Helpers shared by the integration tests.
//!
//! Use from a test file with `mod common;` and then `common::assert_close(..)`.

#![allow(dead_code)]

/// Asserts that `actual` is within `tol` of `expected`.
#[track_caller]
pub fn assert_close(actual: f64, expected: f64, tol: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tol,
        "expected {expected} ± {tol}, got {actual} (off by {:e})",
        (actual - expected).abs()
    );
}

/// Asserts that two angles in degrees are within `tol` degrees, modulo 360.
#[track_caller]
pub fn assert_close_deg(actual: f64, expected: f64, tol: f64) {
    let d = (actual - expected).rem_euclid(360.0);
    let d = d.min(360.0 - d);
    assert!(
        actual.is_finite() && d <= tol,
        "expected {expected}° ± {tol}°, got {actual}° (off by {d:e}°)"
    );
}

/// Asserts that two angles in radians are within `tol` radians, modulo 2π.
#[track_caller]
pub fn assert_close_rad(actual: f64, expected: f64, tol: f64) {
    let tau = std::f64::consts::TAU;
    let d = (actual - expected).rem_euclid(tau);
    let d = d.min(tau - d);
    assert!(
        actual.is_finite() && d <= tol,
        "expected {expected} rad ± {tol}, got {actual} rad (off by {d:e})"
    );
}

/// Degrees, arcminutes and arcseconds to decimal degrees. The sign of `deg`
/// applies to the whole angle; pass `neg = true` for angles like -0°12′.
pub fn dms(neg: bool, deg: u32, min: u32, sec: f64) -> f64 {
    let v = deg as f64 + min as f64 / 60.0 + sec / 3600.0;
    if neg { -v } else { v }
}

/// Hours, minutes and seconds to decimal degrees.
pub fn hms(h: u32, m: u32, s: f64) -> f64 {
    15.0 * (h as f64 + m as f64 / 60.0 + s / 3600.0)
}

/// Arcseconds to degrees.
pub fn arcsec(s: f64) -> f64 {
    s / 3600.0
}

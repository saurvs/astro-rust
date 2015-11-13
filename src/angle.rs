// converts angles expressed with degrees, minutes seconds to pure degrees
pub fn pure_degrees(d: f64, mut m: f64, mut s: f64) -> f64 {
    if d < 0.0 {
        m = -1.0 * m.abs();
        s = -1.0 * s.abs();
    }
    d + (m / 60.0) + (s / 3600.0)
}

pub fn limited_to_360(angle: f64) -> f64 {
    let n = (angle / 360.0) as i64;
    let mut limited_angle = angle - (360.0 * (n as f64));
    if limited_angle < 0.0 {
        limited_angle += 360.0;
    }
    limited_angle
}

pub fn small_angle(angle: f64) -> bool {
    if angle < 0.003 { return true; }
    false
}

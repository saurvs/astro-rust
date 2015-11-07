// converts angles expressed with degrees, minutes seconds to pure degrees
pub fn pure_degrees(d: f64, mut m: f64, mut s: f64) -> f64 {
    if d < 0.0 {
        m = -1.0 * m.abs();
        s = -1.0 * s.abs();
    }
    d + (m / 60.0) + (s / 3600.0)
}

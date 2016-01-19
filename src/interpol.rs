//! Functions for interpolation

pub fn ThreeValues(y1: f64, y2: f64, y3: f64, n: f64) -> f64 {
    let a = y2 - y1;
    let b = y3 - y2;
    let c = b - a;

    y2 + n*(a + b + n*c)/2.0
}

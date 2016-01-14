pub fn int(x: f64) -> i64 {
    if x < 0_f64 {
        return x as i64 - 1;
    }
    x as i64
}

pub fn RoundUptoDigits(number: f64, frac_digits: i32) -> f64 {
    let d = 10_f64.powi(frac_digits);
    (number * d).round() / d
}

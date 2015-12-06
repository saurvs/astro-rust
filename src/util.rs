// Returns the largest integer less than or equal to x
pub fn int(x: f64) -> i64 {
    if x < 0_f64 {
        return x as i64 - 1;
    }
    x as i64
}

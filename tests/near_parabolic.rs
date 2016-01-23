extern crate astro;

use astro::*;

#[test]
fn true_anom_and_rad_vec() {
    let (tru_anom, rad_vec) = orbit::near_parabolic::true_anom_and_rad_vec(138.4783, 0.0, 1.0, 0.921326, 0.0000001);
    assert_eq!(util::RoundUptoDigits(tru_anom.to_degrees(), 5), 102.74426);
    assert_eq!(util::RoundUptoDigits(rad_vec, 6), 2.364192);
}

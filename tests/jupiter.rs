extern crate astro;

use astro::*;

#[test]
fn Ephemeris() {

    let (mut D_e, mut D_s, mut w1, mut w2, mut P) = planet::jupiter::Ephemeris(
        2448972.50068,

        84.285703_f64.to_radians(),
        0.000197_f64.to_radians(),
        0.98412316,

        181.882168_f64.to_radians(),
        1.290464_f64.to_radians(),
        5.4464232,

        23.4402069_f64.to_radians(),

        angle::DegFrmDMS(0, 0, 16.86).to_radians(),
        angle::DegFrmDMS(0, 0, -1.79).to_radians(),
    );

    D_e = util::RoundUptoDigits(D_e.to_degrees(), 2);
    D_s = util::RoundUptoDigits(D_s.to_degrees(), 2);
    P = util::RoundUptoDigits(angle::LimitTo360(P.to_degrees()), 2);
    w1 = util::RoundUptoDigits(w1.to_degrees(), 2);
    w2 = util::RoundUptoDigits(w2.to_degrees(), 2);

    assert_eq!(D_e, -2.48);
    assert_eq!(D_s, -2.20);
    assert_eq!(P, 24.80);
    assert_eq!(w1, 268.06);
    assert_eq!(w2, 72.74);
}

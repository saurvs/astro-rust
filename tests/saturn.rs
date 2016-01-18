extern crate astro;

use astro::*;

#[test]
fn RingPosition() {

    let (mut B, mut B1, mut P, mut deltaU, mut a, mut b) = planet::saturn::RingElements(
        2448972.50068,
/*
        84.285703_f64.to_radians(),
        0.000197_f64.to_radians(),
        0.98412316,

        319.191636_f64.to_radians(),
        -1.075183_f64.to_radians(),
        9.8678819,*/

        angle::DegFrmDMS(0, 0, 16.86).to_radians(),
        23.43971_f64.to_radians()
    );

    B = util::RoundUptoDigits(B.to_degrees(), 3);
    B1 = util::RoundUptoDigits(B1.to_degrees(), 3);
    P = util::RoundUptoDigits(angle::LimitTo360(P.to_degrees()), 3);
    deltaU = util::RoundUptoDigits(deltaU.to_degrees(), 3);
/*
    assert_eq!(B, 16.442);
    assert_eq!(B1, 14.679);
    assert_eq!(P, 6.741);
    assert_eq!(deltaU, 4.198);

    let (h1, m1, sec1) = angle::DMSFrmDeg(a.to_degrees());
    assert_eq!(util::RoundUptoDigits(sec1, 2), 35.87);

    let (h2, m2, sec2) = angle::DMSFrmDeg(b.to_degrees());
    assert_eq!(util::RoundUptoDigits(sec2, 2), 10.15);*/
}

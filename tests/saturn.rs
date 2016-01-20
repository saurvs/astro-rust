extern crate astro;

use astro::*;

#[test]
fn RingPosition() {

    let (mut B, mut B1, mut P, mut deltaU, mut a, mut b) = planet::saturn::ring::Elements(
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
}
/*
#[test]
fn Mimas() {

    let (mut X, mut Y, mut Z) = planet::saturn::moon::ApprntRectCoords(
        2451439.50074,
        &planet::saturn::moon::Moon::Mimas
    );

    X = util::RoundUptoDigits(X, 3);
    Y = util::RoundUptoDigits(Y, 3);

    assert_eq!(X, 3.102);
    assert_eq!(Y, -0.204);
}*/

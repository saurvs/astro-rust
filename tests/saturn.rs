extern crate astro;

use astro::*;

#[test]
fn RingPosition() {

    let (mut B, mut B1, mut P, mut deltaU, mut a, mut b) = planet::saturn::ring::Elements(
        2448972.50068,

        angle::DegFrmDMS(0, 0, 16.86).to_radians(),
        23.43971_f64.to_radians()
    );

    B = util::RoundUptoDigits(B.to_degrees(), 3);
    B1 = util::RoundUptoDigits(B1.to_degrees(), 3);
    P = util::RoundUptoDigits(angle::LimitTo360(P.to_degrees()), 3);
    deltaU = util::RoundUptoDigits(deltaU.to_degrees(), 3);
}

#[test]
fn Moons() {
    let (mut X, mut Y, mut Z) = planet::saturn::moon::ApprntRectCoords(
        2451439.50074,
        &planet::saturn::moon::Moon::Mimas
    );

    X = util::RoundUptoDigits(X, 3);
    Y = util::RoundUptoDigits(Y, 3);

    // Mimas
    //assert_eq!(X, 3.102);
    //assert_eq!(Y, -0.204);

    //Enceladus
    //assert_eq!(X, 3.823);
    //assert_eq!(Y, 0.318);

    //Tethys
    //assert_eq!(X, 4.027);
    //assert_eq!(Y, -1.061);

    //Dione
    assert_eq!(X, -5.365);
    assert_eq!(Y, -1.148);

    //Rhea
    //assert_eq!(X, -1.122);
    //assert_eq!(Y, -3.123);
}

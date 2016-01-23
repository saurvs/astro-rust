extern crate astro;

use astro::*;

#[test]
fn ring_elements() {

    let (mut B, mut B1, mut P, mut deltaU, mut a, mut b) = planet::saturn::ring::elements(
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
fn moons() {

    let data = [
        (3.102, -0.204, planet::saturn::moon::Moon::Mimas),
        (3.823, 0.318, planet::saturn::moon::Moon::Enceladus),
        (4.027, -1.061, planet::saturn::moon::Moon::Tethys),
        (-5.365, -1.148, planet::saturn::moon::Moon::Dione),
        (-1.122, -3.123, planet::saturn::moon::Moon::Rhea),
        (14.568, 4.738, planet::saturn::moon::Moon::Titan),
        (-18.001, -5.328, planet::saturn::moon::Moon::Hyperion),
        (-48.76/* Meeus gives -48.759 */, 4.137/* Meeus gives 4.136 */, planet::saturn::moon::Moon::Iapetus),
    ];

    for tuple in data.iter() {
        let (X, Y, Z) = planet::saturn::moon::apprnt_rect_coords(
            2451439.50074,
            &tuple.2
        );

        assert_eq!(util::RoundUptoDigits(X, 3), tuple.0);
        assert_eq!(util::RoundUptoDigits(Y, 3), tuple.1);
    }

}

extern crate astro;

use astro::*;

#[test]
fn Ephm() {
    let (mut P, mut B, mut L) = sun::Ephm(2448908.50068, 199.90234_f64.to_radians(), 199.906759_f64.to_radians(), 23.440144_f64.to_radians());
    P = util::RoundUptoDigits(P.to_degrees(), 2);
    B = util::RoundUptoDigits(B.to_degrees(), 2);
    L = util::RoundUptoDigits(L.to_degrees(), 2);
    assert_eq!((P, B, L), (26.27, 5.99, 238.63));
}

#[test]
fn CarringSyndRot() {
    assert_eq!(util::RoundUptoDigits(sun::CarringSyndRot(1699), 2), 2444480.72);
}

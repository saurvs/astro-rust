extern crate astro;

use astro::*;

#[test]
fn HeliocenCoords() {
    let (mut L, mut B, mut R) = planet::HeliocenCoords(planet::Planet::Venus, 2448976.5);
    L = util::RoundUptoDigits(L.to_degrees(), 3);
    B = util::RoundUptoDigits(B.to_degrees(), 3);
    R = util::RoundUptoDigits(R, 5);
    assert_eq!((L, B, R), (26.114, util::RoundUptoDigits(angle::LimitTo360(-2.6207), 3), 0.72460));
}

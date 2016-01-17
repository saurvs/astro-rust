#[macro_use]
extern crate astro;

use astro::*;

#[test]
fn HeliocenCoords() {
    let (mut L, mut B, mut R) = planet::HeliocenCoords(&planet::Planet::Venus, 2448976.5);
    L = util::RoundUptoDigits(L.to_degrees(), 3);
    B = util::RoundUptoDigits(B.to_degrees(), 3);
    R = util::RoundUptoDigits(R, 5);
    assert_eq!((L, B, R), (26.114, util::RoundUptoDigits(angle::LimitTo360(-2.6207), 3), 0.72460));
}

#[test]
fn EclGeocenCoords() {
    let (L, B, mut R, mut t) = planet::EclGeocenCoords(
        88.35704_f64.to_radians(),
        0.00014_f64.to_radians(),
        0.983824,

        26.11428_f64.to_radians(),
        -2.6207_f64.to_radians(),
        0.724603,
    );
    R = util::RoundUptoDigits(R, 6);
    t = util::RoundUptoDigits(t, 7);

    assert_eq!((R, t), (0.910845, 0.0052606));
}

#[test]
fn ApprntEclGeocenCoords() {
    let (mut L, mut B, mut R) = ApprntEclGeocenCoords!(&planet::Planet::Venus, 2448976.5);
    L = util::RoundUptoDigits(angle::LimitTo360(L.to_degrees()), 2);
    B = util::RoundUptoDigits(B.to_degrees(), 2);
    R = util::RoundUptoDigits(R, 4);

    assert_eq!((L, B, R), (313.08, -2.08, 0.9109));
}

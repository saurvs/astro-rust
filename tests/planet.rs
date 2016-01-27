#[macro_use]
extern crate astro;

use astro::*;

#[test]
fn heliocen_pos() {
    let (mut L, mut B, mut R) = planet::heliocen_pos(&planet::Planet::Venus, 2448976.5);
    L = util::round_upto_digits(L.to_degrees(), 3);
    B = util::round_upto_digits(B.to_degrees(), 3);
    R = util::round_upto_digits(R, 5);
    assert_eq!(
        (L, B, R),
        (26.114, util::round_upto_digits(angle::limit_to_360(-2.6207), 3), 0.72460)
    );
}

#[test]
fn geocen_geomet_ecl_pos() {
    let (L, B, mut R, mut t) = planet::geocen_geomet_ecl_pos(
        88.35704_f64.to_radians(),
        0.00014_f64.to_radians(),
        0.983824,

        26.11428_f64.to_radians(),
        -2.6207_f64.to_radians(),
        0.724603,
    );
    R = util::round_upto_digits(R, 6);
    t = util::round_upto_digits(t, 7);

    assert_eq!((R, t), (0.910845, 0.0052606));
}

#[test]
fn ecl_coords_to_FK5() {
    let (FK5_long, FK5_lat) = planet::ecl_coords_to_FK5(
        2448976.5,
        313.07689_f64.to_radians(),
        -2.08489_f64.to_radians()
    );

    assert_eq!(
        util::round_upto_digits(FK5_long.to_degrees(), 5),
        313.07686
    );
    assert_eq!(
        util::round_upto_digits(FK5_lat.to_degrees(), 5),
        -2.08487
    );
}

#[test]
fn geocen_apprnt_ecl_pos() {
    let (planet_ecl_point, mut R) = planet::geocen_apprnt_ecl_pos(&planet::Planet::Venus, 2448976.5);
    let (mut L, mut B) = (planet_ecl_point.long, planet_ecl_point.lat);
    L = util::round_upto_digits(angle::limit_to_360(L.to_degrees()), 2);
    B = util::round_upto_digits(B.to_degrees(), 2);
    R = util::round_upto_digits(R, 4);

    assert_eq!((L, B, R), (313.08, -2.08, 0.9109));
}

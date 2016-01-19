extern crate astro;

use astro::*;

#[test]
fn Ephemeris() {
    let (mut P, mut B, mut L) = sun::Ephemeris(2448908.50068,
                                               199.90234_f64.to_radians(),
                                               199.906759_f64.to_radians(),
                                               23.440144_f64.to_radians());
    P = util::RoundUptoDigits(P.to_degrees(), 2);
    B = util::RoundUptoDigits(B.to_degrees(), 2);
    L = util::RoundUptoDigits(L.to_degrees(), 2);
    assert_eq!((P, B, L), (26.27, 5.99, 238.63));
}

#[test]
fn EclCoordsToFK5() {
    let (FK5_long, FK5_lat) = sun::EclCoordsToFK5(
        2448908.5,
        199.907372_f64.to_radians(),
        angle::DegFrmDMS(0, 0, 0.644).to_radians()
    );

    assert_eq!(
        util::RoundUptoDigits(FK5_long.to_degrees(), 6),
        199.907347
    );

    let (d, m, s) = angle::DMSFrmDeg(FK5_lat.to_degrees());
    assert_eq!(
        (0, 0, util::RoundUptoDigits(s, 2)),
        (0, 0, 0.62)
    );
}

#[test]
fn EclGeocenCoords() {
    let (ecl_long, ecl_lat, rad_vec) = sun::EclGeocenCoords(
        2448908.5
    );

    assert_eq!(
        util::RoundUptoDigits(ecl_long.to_degrees(), 6),
        199.907297
    );
    assert_eq!(
        util::RoundUptoDigits(ecl_lat.to_degrees(), 6),
        0.000207
    );
    assert_eq!(
        util::RoundUptoDigits(rad_vec, 8),
        0.99760852
    );
}

#[test]
fn CarringSyndRot() {
    assert_eq!(util::RoundUptoDigits(sun::CarringSyndRot(1699), 2), 2444480.72);
}

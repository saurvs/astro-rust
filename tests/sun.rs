extern crate astro;

use astro::*;

#[test]
fn ephemeris() {
    let (mut P, mut B, mut L) = sun::ephemeris(2448908.50068,
                                               199.90234_f64.to_radians(),
                                               199.906759_f64.to_radians(),
                                               23.440144_f64.to_radians());
    P = util::round_upto_digits(P.to_degrees(), 2);
    B = util::round_upto_digits(B.to_degrees(), 2);
    L = util::round_upto_digits(L.to_degrees(), 2);
    assert_eq!((P, B, L), (26.27, 5.99, 238.63));
}

#[test]
fn ecl_coords_to_FK5() {
    let (FK5_long, FK5_lat) = sun::ecl_coords_to_FK5(
        2448908.5,
        199.907372_f64.to_radians(),
        angle::deg_frm_dms(0, 0, 0.644).to_radians()
    );

    assert_eq!(
        util::round_upto_digits(FK5_long.to_degrees(), 6),
        199.907347
    );

    let (d, m, s) = angle::dms_frm_deg(FK5_lat.to_degrees());
    assert_eq!(
        (0, 0, util::round_upto_digits(s, 2)),
        (0, 0, 0.62)
    );
}

#[test]
fn geocen_ecl_pos() {
    let (sun_eq_point, rad_vec) = sun::geocen_ecl_pos(
        2448908.5
    );

    assert_eq!(
        util::round_upto_digits(sun_eq_point.long.to_degrees(), 6),
        199.907297
    );
    assert_eq!(
        util::round_upto_digits(sun_eq_point.lat.to_degrees(), 6),
        0.000207
    );
    assert_eq!(
        util::round_upto_digits(rad_vec, 8),
        0.99760852
    );
}

#[test]
fn carring_synd_rot() {
    assert_eq!(util::round_upto_digits(sun::carring_synd_rot(1699), 2), 2444480.72);
}

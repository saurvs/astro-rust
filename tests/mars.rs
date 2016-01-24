extern crate astro;

use astro::*;

#[test]
fn ephemeris() {
    let JD = 2448935.500638;
    let (mut D_e, mut D_s, mut P, q, mut w, mut d) = planet::mars::ephemeris(
        JD,

        &planet::mars::ecl_coords_of_north_pol(time::julian_cent(JD)),

        23.44022_f64.to_radians(),

        angle::deg_frm_dms(0, 0, 15.42).to_radians(),
        angle::deg_frm_dms(0, 0, -1.0).to_radians(),
    );

    D_e = util::round_upto_digits(D_e.to_degrees(), 2);
    D_s = util::round_upto_digits(D_s.to_degrees(), 2);
    P = util::round_upto_digits(angle::limit_to_360(P.to_degrees()), 2);
    w = util::round_upto_digits(w.to_degrees(), 1);

    assert_eq!(D_e, 12.44);
    assert_eq!(D_s, -2.76);
    assert_eq!(P, 347.64);
    assert_eq!(w, 111.5);

    let (h1, m1, sec1) = angle::dms_frm_deg(d.to_degrees());
    assert_eq!(util::round_upto_digits(sec1, 2), 10.75);

    let (h2, m2, sec2) = angle::dms_frm_deg(q.to_degrees());
    assert_eq!(util::round_upto_digits(sec2, 2), 1.06);
}

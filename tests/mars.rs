extern crate astro;

use astro::*;

#[test]
fn ephemeris() {
    let JD = 2448935.500638;
    let mut ephemeris = planet::mars::ephemeris(
        JD,

        &planet::mars::north_pol_ecl_coords(time::julian_cent(JD)),

        23.44022_f64.to_radians(),

        angle::deg_frm_dms(0, 0, 15.42).to_radians(),
        angle::deg_frm_dms(0, 0, -1.0).to_radians(),
    );

    ephemeris.De = util::round_upto_digits(ephemeris.De.to_degrees(), 2);
    ephemeris.Ds = util::round_upto_digits(ephemeris.Ds.to_degrees(), 2);
    ephemeris.P = util::round_upto_digits(angle::limit_to_360(ephemeris.P.to_degrees()), 2);
    ephemeris.w = util::round_upto_digits(ephemeris.w.to_degrees(), 1);

    assert_eq!(ephemeris.De, 12.44);
    assert_eq!(ephemeris.Ds, -2.76);
    assert_eq!(ephemeris.P, 347.64);
    assert_eq!(ephemeris.w, 111.5);

    let (h2, m2, sec2) = angle::dms_frm_deg(ephemeris.q.to_degrees());
    assert_eq!(util::round_upto_digits(sec2, 2), 1.06);

    let (h1, m1, sec1) = angle::dms_frm_deg(ephemeris.d.to_degrees());
    assert_eq!(util::round_upto_digits(sec1, 2), 10.75);
}

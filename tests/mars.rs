extern crate astro;

use astro::*;

#[test]
fn Ephemeris() {

    let (mut D_e, mut D_s, mut P, q, mut w, mut d) = planet::mars::Ephemeris(
        2448935.500638,

        352.82267_f64.to_radians(),
        63.28208_f64.to_radians(),

        23.44022_f64.to_radians(),

        angle::DegFrmDMS(0, 0, 15.42).to_radians(),
        angle::DegFrmDMS(0, 0, -1.0).to_radians(),
    );

    D_e = util::RoundUptoDigits(D_e.to_degrees(), 2);
    D_s = util::RoundUptoDigits(D_s.to_degrees(), 2);
    P = util::RoundUptoDigits(angle::LimitTo360(P.to_degrees()), 2);
    w = util::RoundUptoDigits(w.to_degrees(), 1);

    assert_eq!(D_e, 12.44);
    assert_eq!(D_s, -2.76);
    assert_eq!(P, 347.64);
    assert_eq!(w, 111.5);

    let (h1, m1, sec1) = angle::DMSFrmDeg(d.to_degrees());
    assert_eq!(util::RoundUptoDigits(sec1, 2), 10.75);

    let (h2, m2, sec2) = angle::DMSFrmDeg(q.to_degrees());
    assert_eq!(util::RoundUptoDigits(sec2, 2), 1.06);
}

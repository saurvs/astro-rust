extern crate astro;

use astro::*;

#[test]
fn ring_elements() {

    let mut elements = planet::saturn::ring::elements(
        2448972.50068,

        angle::deg_frm_dms(0, 0, 16.86).to_radians(),
        23.43971_f64.to_radians()
    );

    elements.B = util::round_upto_digits(elements.B.to_degrees(), 3);
    elements.B1 = util::round_upto_digits(elements.B1.to_degrees(), 3);
    elements.P = util::round_upto_digits(angle::limit_to_360(elements.P.to_degrees()), 3);
    elements.deltaU = util::round_upto_digits(elements.deltaU.to_degrees(), 3);

    assert_eq!(elements.B, 16.442);
    assert_eq!(elements.B1, 14.679);
    assert_eq!(elements.P, 6.741);
    assert_eq!(elements.deltaU, 4.198);

    let (d_a, m_a, s_a) = angle::dms_frm_deg(elements.a.to_degrees());
    assert_eq!(util::round_upto_digits(s_a, 2), 35.87);

    let (d_b, m_b, s_b) = angle::dms_frm_deg(elements.b.to_degrees());
    assert_eq!(util::round_upto_digits(s_b, 2), 10.15);

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
        assert_eq!(util::round_upto_digits(X, 3), tuple.0);
        assert_eq!(util::round_upto_digits(Y, 3), tuple.1);
    }

}

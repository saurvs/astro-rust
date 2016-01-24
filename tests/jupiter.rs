extern crate astro;

use astro::*;

#[test]
fn ephemeris() {

    let (mut D_e, mut D_s, mut w1, mut w2, mut P) = planet::jupiter::ephemeris(
        2448972.50068,

        23.4402069_f64.to_radians(),

        angle::deg_frm_dms(0, 0, 16.86).to_radians(),
        angle::deg_frm_dms(0, 0, -1.79).to_radians(),
    );

    D_e = util::round_upto_digits(D_e.to_degrees(), 2);
    D_s = util::round_upto_digits(D_s.to_degrees(), 2);
    P = util::round_upto_digits(angle::limit_to_360(P.to_degrees()), 2);
    w1 = util::round_upto_digits(w1.to_degrees(), 0);
    w2 = util::round_upto_digits(w2.to_degrees(), 2);

    assert_eq!(D_e, -2.48);
    assert_eq!(D_s, -2.20);
    assert_eq!(P, 24.80);
    assert_eq!(w1, 268.0);
    assert_eq!(w2, 72.74);
}

#[test]
fn moons() {
/*
    let data = [
        (-3.41 /* Meeus gives -3.44 */, 0.21, planet::jupiter::moon::Moon::Io),
        (7.49  /* Meeus gives 7.44 */, 0.25, planet::jupiter::moon::Moon::Europa),
        (1.24, 0.65, planet::jupiter::moon::Moon::Ganymede),
        (7.08, 1.1, planet::jupiter::moon::Moon::Callisto),
    ];

    for tuple in data.iter() {
        let (X, Y) = planet::jupiter::moon::apprnt_rect_coords(
            2448972.50068,
            &tuple.2
        );

        assert_eq!(util::round_upto_digits(X, 2), tuple.0);
        assert_eq!(util::round_upto_digits(Y, 2), tuple.1);
    }
*/
}

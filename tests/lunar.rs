extern crate astro;

use astro::*;

#[test]
fn geocen_ecl_pos() {
    let (moon_eq_point, rad_vec) = lunar::geocen_ecl_pos(
        2448724.5
    );

    assert_eq!(
        util::round_upto_digits(moon_eq_point.long.to_degrees(), 6),
        133.162655
    );
    assert_eq!(
        util::round_upto_digits(moon_eq_point.lat.to_degrees(), 6),
        -3.229126
    );
    assert_eq!(
        util::round_upto_digits(rad_vec, 1),
        368409.7
    );
    
}

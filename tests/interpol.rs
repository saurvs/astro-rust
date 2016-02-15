extern crate astro;

use astro::*;

#[test]
fn three_values() {
    let y = interpol::three_values(
        0.884226,
        0.877366,
        0.870531,
        0.18125
    );

    assert_eq!(util::round_upto_digits(y, 6), 0.876125);
}

#[test]
fn five_values() {
    let y = interpol::five_values(
        36.125,
        24.606,
        15.486,
        8.694,
        4.133,
        0.2777778
    );

    assert_eq!(util::round_upto_digits(y, 3), 13.369);
}

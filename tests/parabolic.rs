extern crate astro;

use astro::*;

#[test]
fn TruAnomAndRadVec() {
    let t_date = time::Date{year: 1998, month: 8, decimal_day: 5.0, cal_type: time::CalType::Gregorian};
    let t = time::JulDay(&t_date);

    let T_date = time::Date{year: 1998, month: 4, decimal_day: 14.4358, cal_type: time::CalType::Gregorian};
    let T = time::JulDay(&T_date);

    let (tru_anom, rad_vec) = orbit::parabolic::TruAnomAndRadVec(t, T, 1.487469);
    assert_eq!(util::RoundUptoDigits(tru_anom.to_degrees(), 5), 66.78862);
    assert_eq!(util::RoundUptoDigits(rad_vec, 6), 2.133911);
}

#[test]
fn PassageThroughNodes() {

    let T = time::JulDay(&time::Date{
        year: 1989,
        month: 8,
        decimal_day: 20.291,
        cal_type: time::CalType::Gregorian
    });
    let w = 154.9103_f64.to_radians();
    let q = 1.324502;

    let (ascen, r_a) = orbit::parabolic::PassageThroughNode(w, q, T, &orbit::Node::Ascend);
    assert_eq!(util::RoundUptoDigits((T - ascen), 2), 4354.65);
    assert_eq!(util::RoundUptoDigits(r_a, 2), 28.07);

    let (descend, r_b) = orbit::parabolic::PassageThroughNode(w, q, T, &orbit::Node::Descend);
    assert_eq!(util::RoundUptoDigits((T - descend), 4), -28.3454);
    assert_eq!(util::RoundUptoDigits(r_b, 4), 1.3901);
}

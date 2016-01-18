extern crate astro;

use astro::*;

#[test]
fn EccAnom() {
    let ecc_anom = orbit::elliptic::EccAnom(5_f64.to_radians(), 0.1, 0.0000001);
    assert_eq!(util::RoundUptoDigits(ecc_anom.to_degrees(), 6), 5.554589);
}

#[test]
fn Velocity() {
    let a = 17.9400782;
    let e = 0.96727426;
    let (vel_p, vel_a) = (orbit::elliptic::PerihVel(a, e), orbit::elliptic::AphVel(a, e));
    assert_eq!(util::RoundUptoDigits(vel_p, 2), 54.52);
    assert_eq!(util::RoundUptoDigits(vel_a, 2), 0.91);

    let V = orbit::elliptic::Vel(1.0, a);
    assert_eq!(util::RoundUptoDigits(V, 2), 41.53);
}

#[test]
fn Length() {
}

#[test]
fn PassageThroughNodes() {

    let T = time::JulDay(&time::Date{
        year: 1986,
        month: 2,
        decimal_day: 9.45891,
        cal_type: time::CalType::Gregorian
    });
    let w = 111.84644_f64.to_radians();
    let e = 0.96727426;
    let n = 0.01297082_f64.to_radians();
    let a = 17.9400782;

    let (ascen, r_a) = orbit::elliptic::PassageThroughNode(w, n, a, e, T, &orbit::Node::Ascend);
    assert_eq!(util::RoundUptoDigits((T - ascen), 4), 92.2998);
    assert_eq!(util::RoundUptoDigits(r_a, 4), 1.8045);

    let (descend, r_b) = orbit::elliptic::PassageThroughNode(w, n, a, e, T, &orbit::Node::Descend);
    assert_eq!(util::RoundUptoDigits((T - descend), 4), -28.9105);
    assert_eq!(util::RoundUptoDigits(r_b, 4), 0.8493);
}
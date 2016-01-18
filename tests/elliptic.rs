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
    let a = 17.9400782;
    let e = 0.96727426;
    let (vel_p, vel_a) = (orbit::elliptic::PerihVel(a, e), orbit::elliptic::AphVel(a, e));
    assert_eq!(util::RoundUptoDigits(vel_p, 2), 54.52);
    assert_eq!(util::RoundUptoDigits(vel_a, 2), 0.91);

    let V = orbit::elliptic::Vel(1.0, a);
    assert_eq!(util::RoundUptoDigits(V, 2), 41.53);
}

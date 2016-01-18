extern crate astro;

use astro::*;

#[test]
fn RhoSinAndCosPhi() {
    let (rho_sin_phi, rho_cos_phi) = planet::earth::RhoSinCosPhi(angle::DegFrmDMS(33, 21, 22.0).to_radians(), 1706.0);
    assert_eq!((util::RoundUptoDigits(rho_sin_phi, 6),
                util::RoundUptoDigits(rho_cos_phi, 6)), (0.546861, 0.836339));
}

#[test]
fn GeodesicDist() {
    let paris = coords::GeographPoint{long: angle::DegFrmDMS(-2, 20, 14.0).to_radians(),
                                      lat : angle::DegFrmDMS(48, 50, 11.0).to_radians()};

    let washington = coords::GeographPoint{long: angle::DegFrmDMS(77,  3, 56.0).to_radians(),
                                           lat : angle::DegFrmDMS(38, 55, 17.0).to_radians()};

    let distance = planet::earth::GeodesicDist(&paris, &washington);
    assert_eq!(util::RoundUptoDigits(distance, 2), 6181.63);

    let approx_distance = planet::earth::ApproxGeodesicDist(&paris, &washington);
    assert_eq!(util::RoundUptoDigits(approx_distance, 0), 6166.0);
}

#[test]
fn Radii() {
    let lat = 42_f64.to_radians();

    let Rp = planet::earth::RadOfParllLat(lat);
    assert_eq!(util::RoundUptoDigits(Rp, 0), util::RoundUptoDigits(4747.001, 0));

    let lin_vel = planet::earth::LinearVelAtLat(lat);
    assert_eq!(util::RoundUptoDigits(lin_vel, 5), 0.34616);

    let Rm = planet::earth::RadOfCurvOfMerd(lat);
    assert_eq!(util::RoundUptoDigits(Rm, 2), util::RoundUptoDigits(6364.033, 2));
}

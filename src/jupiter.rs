use angle;
use ecliptic;
use earth;

/**
Computes Jupiter's **equatorial semidiameter**

# Arguments

* ```DistanceToEarth```: Jupiter's distance from Earth *(AU)*
**/
pub fn EquatorialSemidiameter(DistanceToEarth: f64) -> f64 {
    angle::PureDegrees(0.0, 0.0, 98.44) / DistanceToEarth
}

/**
Computes Jupiter's **polar semidiameter**

# Arguments

* ```DistanceToEarth```: Jupiter's distance from Earth *(AU)*
**/
pub fn PolarSemidiameter(DistanceToEarth: f64) -> f64 {
    angle::PureDegrees(0.0, 0.0, 92.06) / DistanceToEarth
}

/**
Computes Jupiter's

# Arguments

* ```jed```: Julian Ephemeris day
**/
pub fn Ephemeris(jed: f64) -> (f64, f64, f64, f64, f64) {
    let d = jed - 2433282.5;
    let T1 = d / 36525.0;

    let asc0 = 268.0_f64.to_radians() + 0.1061_f64.to_radians()*T1;
    let dec0 = 64.5_f64.to_radians() - 0.0164_f64.to_radians()*T1;
    let sin_dec0 = dec0.sin();
    let cos_dec0 = dec0.cos();

    let W1 = angle::LimitedTo360(17.710_f64.to_radians() + 877.90003539_f64.to_radians()*d);
    let W2 = angle::LimitedTo360(16.838_f64.to_radians() + 870.27003539_f64.to_radians()*d);

    let l0: f64 = 0.0; // heliocentric longitude of the Earth
    let b0: f64 = 0.0; // heliocentric latitude of the Earth
    let R: f64 = 0.0; // radius vector of the Earth

    let mut l = 0.0; // heliocentric longitude of Jupiter (without light-time correction)
    let b: f64 = 0.0; // heliocentric latitude of Jupiter (without light-time correction)
    let r: f64 = 0.0; // radius vector of Jupiter (without light-time correction)
    let r_squared = r * r;

    let mut x: f64 = 0.0;
    let mut y: f64 = 0.0;
    let mut z: f64 = 0.0;
    let mut delta: f64 = 0.0;

    l -= 0.012990_f64.to_radians()*delta / r_squared;

    // recalculate x y z delta

    let e0 = ecliptic::MeanObliquity(jed);
    let cos_e0 = e0.cos();
    let sin_e0 = e0.sin();
    let sin_l = l.sin();

    let asc_s = (cos_e0*sin_l - sin_e0*b.tan()).atan2(l.cos());
    let dec_s = (cos_e0*b.sin() + sin_e0*b.cos()*sin_l).asin();

    let D_s = (-1.0*sin_dec0*dec_s.sin() - cos_dec0*dec_s.cos()*(asc0 - asc_s).cos()).asin();

    let mut asc: f64 = 0.0; // more stufz
    let sin_asc = asc.sin();
    let cos_asc = asc.cos();
    let mut dec: f64 = 0.0;
    let sin_dec = dec.sin();
    let cos_dec = dec.cos();
    let zeta: f64 = 0.0;

    let D_e = (-1.0*sin_dec0*sin_dec - cos_dec0*dec.sin()*(asc0 - asc).cos()).asin();

    let mut w1 = angle::LimitedTo360(W1.to_degrees() - zeta.to_degrees() - 5.07033*delta);
    let mut w2 = angle::LimitedTo360(W2.to_degrees() - zeta.to_degrees() - 5.02626*delta);

    let mut C = 57.2958 * (2.0*r*delta + R*R - r_squared - delta*delta) / (4.0*r*delta);
    if (l - l0).sin() < 0.0 {
        C *= -1.0
    }
    w1 += C;
    w2 += C;

    let (nut_long, nut_obl) = earth::NutationCorrections(jed);
    let e = e0 + nut_obl;

    let q = 0.005693_f64.to_radians();
    let sin_l0 = l0.sin();
    let cos_l0 = l0.cos();
    let asc_correction: f64 = q * (asc.cos()*cos_l0*e.cos() + sin_asc*sin_l0) / cos_dec;
    let dec_correction: f64 = q * (cos_l0*e.cos()*(e.tan()*cos_dec - sin_asc*cos_asc) + cos_asc*sin_dec*sin_l0);

    asc += asc_correction;
    dec += dec_correction;

    // correct asc, asc0, dec, dec0 for nutation

    let P = 0.0;

    (D_e, D_s, w1, w2, P)
}

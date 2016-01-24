//! The four Galilean moons
use angle;use util;use planet;
/*

Meeus referrs to the moons as Satellites I, II, III and IV.
Wikipedia says satellite I is Io, II is Europa, III is Ganymede
and IV is Callisto. That mapping is used.

*/

/// Represents a Galilean moon
pub enum Moon {
    /// Io
    Io,
    /// Europa
    Europa,
    /// Ganymede
    Ganymede,
    /// Callisto
    Callisto
}

/**
Returns the apparent **rectangular** coordinates for a moon of Jupiter

This function implements the low accuracy method described in Meeus's
book, that is *"sufficient for identifying the satellites at the
telescope, or drawing a wavy-line diagram showing their positions with
respect to Jupiter"*

# Returns

```(X, Y)```

The rectangular coordinates returned give the apparent position of a moon,
with respect to Jupiter, as seen from Earth. The ```X``` and ```Y```
coordinates are measured from the center of the disk of Jupiter, in units
of Saturn's equatorial radius.

```X``` is measued positively to the west of Jupiter, and negatively to the
east. The x-axis coincides with the equator of the planet.

```Y``` is measured positively to the north of Jupiter, and negatively to
the south. The y-axis coincides with the planet's axis of rotation.

# Arguments

* ```JD```: Julian (Ephemeris) day
* ```moon```: The [Moon](./enum.Moon.html)
**/fn g(){}/*
fn apprnt_rect_coords(JD: f64, moon: &Moon) -> (f64, f64) {
    let d = JD - 2451545.0;println!("d = {:?}", util::round_upto_digits(d, 5));
    let V = (172.74 + 0.00111588*d).to_radians();println!("V = {:?}", util::round_upto_digits(V.to_degrees(), 2));
    let M = (357.529 + 0.9856003*d).to_radians();println!("M = {:?}", util::round_upto_digits(M.to_degrees(), 3));
    let N = (20.02 + 0.0830853*d + 0.329*V.sin()).to_radians();println!("N = {:?}", util::round_upto_digits(N.to_degrees(), 3));
    let J = (66.115 + 0.9025179*d - 0.329*V.sin()).to_radians();println!("J = {:?}", util::round_upto_digits(angle::limit_to_360(J.to_degrees()), 3));

    let A = (1.915*M.sin() + 0.02*(2.0*M).sin()).to_radians();println!("A = {:?}", util::round_upto_digits(A.to_degrees(), 3));
    let B = (5.555*N.sin() + 0.168*(2.0*N).sin()).to_radians();println!("B = {:?}", util::round_upto_digits(B.to_degrees(), 3));
    let K = J + A - B;println!("K = {:?}", util::round_upto_digits(angle::limit_to_360(K.to_degrees()), 3));
    let R = 1.00014 - 0.01671*M.cos() - 0.00014*(2.0*M).cos();println!("R = {:?}", util::round_upto_digits(R, 5));
    let r = 5.20872 - 0.25208*N.cos() - 0.00611*(2.0*N).cos();println!("r = {:?}", util::round_upto_digits(r, 5));
    let delta = (r*r + R*R - 2.0*r*R*K.cos()).sqrt();println!("delta = {:?}", util::round_upto_digits(delta, 5));
    let mut phi = (R*K.sin()/delta).asin();println!("phi = {:?}", util::round_upto_digits(phi.to_degrees(), 3));


    phi += phi_correction(JD, r, delta, R);

    let d_minus_delta_by_173 = d - delta/173.0;println!("d_minus_delta_by_173 = {:?}", util::round_upto_digits(d_minus_delta_by_173, 5));
    let phi_minus_B = phi - B;

    let u1 = (163.8069 + 203.4058646*d_minus_delta_by_173).to_radians() + phi_minus_B;println!("u1 = {:?}", util::round_upto_digits(u1.to_degrees(), 3));
    let u2 = (358.414  + 101.2916335*d_minus_delta_by_173).to_radians() + phi_minus_B;println!("u2 = {:?}", angle::limit_to_360(u2.to_degrees()));
    let u3 = (5.7176   + 50.234518*d_minus_delta_by_173).to_radians()   + phi_minus_B;println!("u3 = {:?}", angle::limit_to_360(u3.to_degrees()));

    let mut u = match moon {
        &Moon::Io       => u1,
        &Moon::Europa   => u2,
        &Moon::Ganymede => u3,
        &Moon::Callisto => (224.8092 + 21.48798*d_minus_delta_by_173).to_radians() + phi_minus_B,
    };println!("u = {:?}", u.to_degrees());

    let G = (331.18 + 50.310482*d_minus_delta_by_173).to_radians();println!("G = {:?}", G.to_degrees());
    let H = (87.45  + 21.569231*d_minus_delta_by_173).to_radians();println!("H = {:?}", H.to_degrees());

    u += (match moon {
        &Moon::Io       => 0.473 * (2.0*(u1 - u2)).sin(),
        &Moon::Europa   => 1.065 * (2.0*(u2 - u3)).sin(),
        &Moon::Ganymede => 0.165 * G.sin(),
        &Moon::Callisto => 0.843 * H.sin(),
    }).to_radians();println!("u+correction = {:?}", angle::limit_to_360(u.to_degrees()));

    let r_moon = match moon {
        &Moon::Io       => (5.9057  - 0.0244*(2.0*(u1 - u2)).cos()),
        &Moon::Europa   => (9.3966  - 0.0882*(2.0*(u2 - u3)).cos()),
        &Moon::Ganymede => (14.9883 - 0.0216*G.cos()),
        &Moon::Callisto => (26.3627 - 0.1939*H.cos()),
    };println!("r_moon = {:?}", util::round_upto_digits(r_moon, 5));

    let lambda = (34.35 + 0.083091*d + 0.329*V.sin()).to_radians() + B;println!("lambda = {:?}", util::round_upto_digits(lambda.to_degrees(), 2));
    let Ds = (3.12 * (lambda + 42.8_f64.to_radians()).sin()).to_radians();println!("Ds = {:?}", util::round_upto_digits(Ds.to_degrees(), 3));
    let De = Ds - (
        2.22*phi.sin()*(lambda + 22_f64.to_radians()).cos()
        + 1.3*(r - delta)*(lambda - 100.5_f64.to_radians()).sin()/delta
    ).to_radians();println!("De = {:?}", util::round_upto_digits(De.to_degrees(), 2));

    println!("1 = {:?}", util::round_upto_digits(2.22*phi.sin()*(lambda + 22_f64.to_radians()).cos(), 3));
    println!("2 = {:?}", util::round_upto_digits(1.3*(r - delta)*(lambda - 100.5_f64.to_radians()).sin()/delta, 3));

    let X =  r_moon * u.sin();
    let Y = -r_moon * u.cos() * De.sin();

    println!(" ");

    (X, Y)
}

fn phi_correction(JD: f64, r: f64, delta: f64, R: f64) -> f64 {
    let mut C = (57.2958 * (2.0*r*delta + R*R - r*r - delta*delta)/(4.0*r*delta)).to_radians();

    let (l0, b0, jup_sun_dist) = planet::heliocen_pos(&planet::Planet::Earth, JD);
    let mut l = 0.0; let mut b = 0.0; let mut r = 0.0;
    let mut x = 0.0; let mut y = 0.0; let mut z = 0.0;
    let mut jup_earth_dist = 0.0; let mut light_time = 0.0;
    let mut i: u8 = 1;
    let n: u8 = 2;
    while i <= n {
        let (new_l, new_b, new_r) = planet::heliocen_pos(&planet::Planet::Jupiter, JD - light_time);
        l = new_l; b = new_b; r = new_r;

        let (new_x, new_y, new_z) = planet::geocen_ecl_rect_coords(l0, b0, R, l, b, r);
        x = new_x; y = new_y; z = new_z;

        jup_earth_dist = planet::dist_frm_ecl_rect_coords(x, y, z);
        light_time = planet::light_time(jup_earth_dist);

        i += 1;
    }
    l -= 0.01299_f64.to_radians()*jup_earth_dist / (r*r);

    if (l - l0).sin() < 0.0 { C *= -1.0 }
println!("C = {:?}", util::round_upto_digits(C.to_degrees(), 2));
    C
}*/

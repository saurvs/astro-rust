//! The four Galilean moons
use angle;use util;use planet;use time;
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
Returns the apparent **rectangular** coordinates for a Galilean moon

This function implements the low accuracy method described in Meeus's
book, that is *"sufficient for identifying the satellites at the
telescope, or drawing a wavy-line diagram showing their positions with
respect to Jupiter"*

# Returns

`(X, Y)`

The rectangular coordinates returned give the apparent position of a moon,
with respect to Jupiter, as seen from Earth. The `X` and `Y`
coordinates are measured from the center of the disk of Jupiter, in units
of Saturn's equatorial radius.

`X` is measued positively to the west of Jupiter, and negatively to the
east. The x-axis coincides with Jupiter's equator.

`Y` is measured positively to the north of Jupiter, and negatively to
the south. The y-axis coincides with Jupiter's axis of rotation.

# Arguments

* `JD`: Julian (Ephemeris) day
* `moon`: The [Moon](./enum.Moon.html)
**/
pub fn apprnt_rect_coords(JD: f64, moon: &Moon) -> (f64, f64) {
    let d = JD - 2451545.0;
    let V = (172.74 + 0.00111588*d).to_radians();
    let M = (357.529 + 0.9856003*d).to_radians();
    let N = (20.02 + 0.0830853*d + 0.329*V.sin()).to_radians();
    let J = (66.115 + 0.9025179*d - 0.329*V.sin()).to_radians();
    let A = (1.915*M.sin() + 0.02*(2.0*M).sin()).to_radians();
    let B = (5.555*N.sin() + 0.168*(2.0*N).sin()).to_radians();
    let K = J + A - B;
    let R = 1.00014 - 0.01671*M.cos() - 0.00014*(2.0*M).cos();
    let r = 5.20872 - 0.25208*N.cos() - 0.00611*(2.0*N).cos();
    let delta = (r*r + R*R - 2.0*r*R*K.cos()).sqrt();

    let phi = (R*K.sin()/delta).asin();

    let d_minus_delta_by_173 = d - delta/173.0;
    let phi_minus_B = phi - B;

    let u1 = (163.8069 + 203.4058646*d_minus_delta_by_173).to_radians() + phi_minus_B;
    let u2 = (358.414  + 101.2916335*d_minus_delta_by_173).to_radians() + phi_minus_B;
    let u3 = (5.7176   + 50.234518*d_minus_delta_by_173).to_radians()   + phi_minus_B;

    let mut u = match moon {
        &Moon::Io       => u1,
        &Moon::Europa   => u2,
        &Moon::Ganymede => u3,
        &Moon::Callisto => (224.8092 + 21.48798*d_minus_delta_by_173).to_radians() + phi_minus_B,
    };

    let G = (331.18 + 50.310482*d_minus_delta_by_173).to_radians();
    let H = (87.45  + 21.569231*d_minus_delta_by_173).to_radians();

    u += (match moon {
        &Moon::Io       => 0.473 * (2.0*(u1 - u2)).sin(),
        &Moon::Europa   => 1.065 * (2.0*(u2 - u3)).sin(),
        &Moon::Ganymede => 0.165 * G.sin(),
        &Moon::Callisto => 0.843 * H.sin(),
    }).to_radians();

    let r_moon = match moon {
        &Moon::Io       => (5.9057  - 0.0244*(2.0*(u1 - u2)).cos()),
        &Moon::Europa   => (9.3966  - 0.0882*(2.0*(u2 - u3)).cos()),
        &Moon::Ganymede => (14.9883 - 0.0216*G.cos()),
        &Moon::Callisto => (26.3627 - 0.1939*H.cos()),
    };

    let lambda = (34.35 + 0.083091*d + 0.329*V.sin()).to_radians() + B;
    let Ds = (3.12 * (lambda + 42.8_f64.to_radians()).sin()).to_radians();
    let De = Ds - (
        2.22*phi.sin()*(lambda + 22_f64.to_radians()).cos()
        + 1.3*(r - delta)*(lambda - 100.5_f64.to_radians()).sin()/delta
    ).to_radians();

    let X =  r_moon * u.sin();
    let Y = -r_moon * u.cos() * De.sin();

    (X, Y)
}

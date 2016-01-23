//! Ring of Saturn

use angle;
use coords;
use planet;
use time;

pub fn inc(JC: f64) -> f64 {
    (      28.075216
      - JC*(0.012998
      + JC*0.000004)
    ).to_radians()
}

pub fn ascend_node(JC: f64) -> f64 {
    (      169.50847
      + JC*(1.394681
      + JC*0.000412)
    ).to_radians()
}

/**
Returns the **elements** for the **ring** of Saturn

# Returns

```(B, B1, P, deltaU, a, b)```

* ```B```: Saturnicentric latitude of the Earth,
           referred to the plane of the ring *| in radians*
* ```B1```: Saturnicentric latitude of the Sun,
            referred to the plane of the ring *| in radians*
* ```P```: Geocentric position angle of the northern
           semiminor axis of the apparent ellipse of the
           ring *| in radians*
* ```deltaU```: Difference between Saturnicentric
                longitudes of the Sun and the Earth,
                measured in the plane of the ring *| in radians*
* ```a```: Major axis of the outer edge of the outer ring *| in radians*
* ```b```: Minor axis of the outer edge of the outer ring *| in radians*

# Arguments

* ```JD```: Julian (Ephemeris) day
* ```nut_in_long```: Nutation in longitude on ```JD``` *| in radians*
* ```tru_oblq_eclip```: True obliquity of the ecliptic on ```JD``` *| in radians*
**/
pub fn elements(JD: f64, nut_in_long: f64, tru_oblq_eclip: f64) -> (f64, f64, f64, f64, f64, f64) {
    let (l0, b0, R) = planet::heliocen_pos(&planet::Planet::Earth, JD);

    let mut l = 0.0; let mut b = 0.0; let mut r = 0.0;
    let mut x = 0.0; let mut y = 0.0; let mut z = 0.0;
    let mut saturn_earth_dist = 0.0;
    let mut light_time = 0.0;

    let mut i: u8 = 1;
    let n: u8 = 2;
    while i <= n {
        let (new_l, new_b, new_r) = planet::heliocen_pos(&planet::Planet::Saturn, JD - light_time);
        l = new_l; b = new_b; r = new_r;

        let (new_x, new_y, new_z) = planet::geocen_ecl_rect_coords(l0, b0, R, l, b, r);
        x = new_x; y = new_y; z = new_z;

        saturn_earth_dist = planet::dist_frm_ecl_rect_coords(x, y, z);
        light_time = planet::light_time(saturn_earth_dist);

        i += 1;
    }

    let JC = time::JulCent(JD);
    let inc = inc(JC);
    let ascend_node = ascend_node(JC);

    let (mut lambda, mut beta) = planet::ecl_coords_frm_ecl_rect_coords(x, y, z);
    let B = (  inc.sin() * beta.cos() * (lambda - ascend_node).sin()
             - inc.cos() * beta.sin()
            ).asin();
    let semi_maj = angle::DegFrmDMS(0, 0, 375.35).to_radians() / saturn_earth_dist;
    let semi_min = semi_maj * B.abs().sin();

    let N = (113.6655 + 0.8771*JC).to_radians();

    let l1 = l - (0.01759/r).to_radians();
    let b1 = b - (0.000764*(l - N).cos()/r).to_radians();

    let B1 = (   inc.sin() * b1.cos() * (l1 - ascend_node).sin()
               - inc.cos() * b1.sin()
             ).asin();
    let U1 = (inc.sin()*b1.sin() + inc.cos()*b1.cos()*(l1 - ascend_node).sin())
             .atan2(b1.cos()*(l1 - ascend_node).cos());
    let U2 = (inc.sin()*beta.sin() + inc.cos()*beta.cos()*(lambda - ascend_node).sin())
             .atan2(beta.cos()*(lambda - ascend_node).cos());
    let deltaU = (U1 - U2).abs();

    let mut lambda0 = ascend_node - 90.0_f64.to_radians();
    let beta0 = 90.0_f64.to_radians() - inc;

    let q = 0.005693_f64.to_radians();
    lambda += q * (l0 - lambda).cos() / beta.cos();
    beta += q * (l0 - lambda).sin() * beta.sin();

    lambda0 += nut_in_long;
    lambda += nut_in_long;

    let asc0 = coords::AscFrmEcl(lambda0, beta0, tru_oblq_eclip);
    let dec0 = coords::DecFrmEcl(lambda0, beta0, tru_oblq_eclip);
    let asc = coords::AscFrmEcl(lambda, beta, tru_oblq_eclip);
    let dec = coords::DecFrmEcl(lambda, beta, tru_oblq_eclip);

    let P = (dec0.cos() * (asc0 - asc).sin())
            .atan2(dec0.sin()*dec.cos() - dec0.cos()*dec.sin()*(asc0 - asc).cos());

    (B, B1, P, deltaU, semi_maj, semi_min)
}

pub fn inn_edge_outer_ring(a: f64, b: f64) -> (f64, f64) {
    (a*0.8801, b*0.8801)
}

pub fn out_edge_inner_ing(a: f64, b: f64) -> (f64, f64) {
    (a*0.8599, b*0.8599)
}

pub fn inn_edge_inner_ring(a: f64, b: f64) -> (f64, f64) {
    (a*0.665, b*0.665)
}

pub fn inn_edge_dusk_ring(a: f64, b: f64) -> (f64, f64) {
    (a*0.5486, b*0.5486)
}

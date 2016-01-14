use angle;

#[macro_use]
use coords;

pub fn Inc(JC: f64) -> f64 {
    (      28.075216
      - JC*(0.012998
      + JC*0.000004)
    ).to_radians()
}

pub fn AscenNode(JC: f64) -> f64 {
    (      169.50847
      - JC*(1.394681
      + JC*0.000412)
    ).to_radians()
}/*

pub Stuff(inc: f64, ascend_node: f64,
          l0: f64, b0: f64, R: f64,
          l: f64, b: f64, r: f64) {

    let x = r*b.cos()*l.cos() - R*l0.cos();
    let y = r*b.cos()*l.sin() - R*l0.sin();
    let z = r*b.sin()         - R*b0.sin();

    let saturn_earth_dist = (x*x + y*y + z*z).sqrt();

    let lambda = y.atan2(x);
    let beta = z.atan2((x*x + y*y).sqrt());
    let B = (  inc.sin() * beta.cos() * (lambda - ascend_node).sin()
             - inc.cos() * beta.sin()
            ).asin();
    let a = angle::DegFrmDMS(0, 0, 375.35, saturn_earth_dist);
    let b = a * B.abs().sin();

    let N = 113.6655 + 0.8771*JC;

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

    lambda0 += nut_in_long;
    lambda += nut_in_long;

    let (asc0, dec0) = EqFrmEcl!(lambda0, beta0, true_oblq_eclip);
    let (asc, dec) = EqFrmEcl!(lambda, beta, true_oblq_eclip);

    let P = (dec0.cos() * (asc0 - asc).sin())
            .atan2(dec0.sin()*dec.cos() - dec0.cos()*dec.sin()*(asc0 - asc).cos());

}
*/

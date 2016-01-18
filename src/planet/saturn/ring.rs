use angle;
use time;
use planet;
use coords;

pub fn Inc(JC: f64) -> f64 {
    (      28.075216
      - JC*(0.012998
      + JC*0.000004)
    ).to_radians()
}

pub fn AscenNode(JC: f64) -> f64 {
    (      169.50847
      + JC*(1.394681
      + JC*0.000412)
    ).to_radians()
}

/**
Returns the elements for the ring of Saturn

# Returns

```(B, B1, P, deltaU, a, b)```

* ```B```: Saturnicentric latitude of the Earth (*radians*)
* ```B1```: Saturnicentric latitude of the Sun (*radians*)
* ```P```: Geocentric position angle of the northern
           semiminor axis of the apparent ellipse of the
           ring(*radians*)
* ```deltaU```: Difference between Saturnicentric
                longitudes of the Sun and the Earth (*radians*)
* ```a```: Major axis of the outer edge of the outer ring (*radians*)
* ```b```: Minor axis of the outer edge of the outer ring (*radians*)

# Arguments

* ```JD```: Julian (Ephemeris) day
* ```l0```: Heliocentric longitude of the Earth (*radians*)
            on ```JD```
* ```b0```: Heliocentric latitude of the Earth (*radians*)
            on ```JD```
* ```R```: Heliocentric radius vector of the Earth (*AU*)
            on ```JD```

If ```tau``` is the light time from Saturn to the Earth

* ```l```: Heliocentric longitude of the Earth (*radians*)
           on ```JD - tau```
* ```b```: Julian (Ephemeris) day
* ```r```: Apparent heliocentric of  (*radians*)
* ```nut_in_long```: Nutation in longitude (*radians*)
* ```tru_oblq_eclip```: True obliquity of the ecliptic (*radians*)
**/

pub fn Elements(JD: f64/*,
          l0: f64, b0: f64, R: f64,
          l: f64, b: f64, r: f64,*/,
          nut_in_long: f64, tru_oblq_eclip: f64) -> (f64, f64, f64, f64, f64, f64) {
    let (l0,b0,R) = planet::HeliocenCoords(&planet::Planet::Earth, JD);
    let JC = time::JulCent(JD);
    let inc = Inc(JC);
    let ascend_node = AscenNode(JC);

    let (mut l,mut b,mut r) = planet::HeliocenCoords(&planet::Planet::Saturn, JD);
    let mut x = r*b.cos()*l.cos() - R*l0.cos();
    let mut y = r*b.cos()*l.sin() - R*l0.sin();
    let mut z = r*b.sin()         - R*b0.sin();
    let mut saturn_earth_dist = (x*x + y*y + z*z).sqrt();
    let light_time = planet::LightTime(saturn_earth_dist);
    let (t1,t2,t3) = planet::HeliocenCoords(&planet::Planet::Saturn, JD-light_time);
    x = r*b.cos()*l.cos() - R*l0.cos();
    y = r*b.cos()*l.sin() - R*l0.sin();
    z = r*b.sin()         - R*b0.sin();
    saturn_earth_dist = (x*x + y*y + z*z).sqrt();


    let mut lambda = y.atan2(x);
    let mut beta = z.atan2((x*x + y*y).sqrt());
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

pub fn InnEdgeOutRing(a: f64, b: f64) -> (f64, f64) {
    (a*0.8801, b*0.8801)
}

pub fn OutEdgeInnRing(a: f64, b: f64) -> (f64, f64) {
    (a*0.8599, b*0.8599)
}

pub fn InnEdgeInnRing(a: f64, b: f64) -> (f64, f64) {
    (a*0.665, b*0.665)
}

pub fn InnEdgeDuskRing(a: f64, b: f64) -> (f64, f64) {
    (a*0.5486, b*0.5486)
}

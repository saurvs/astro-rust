//! Eight moons of Saturn

use angle;
use coords;
use planet;
use time;

/// Represents a moon of Saturn
pub enum Moon {
    /// Mimas
    Mimas,
    /// Enceladus
    Enceladus,
    /// Tethys
    Tethys,
    /// Dione
    Dione,
    /// Rhea
    Rhea,
    /// Titan
    Titan,
    /// Hyperion
    Hyperion,
    /// Iapetus
    Iapetus
}

/**
Returns the apparent **rectangular** coordinates for a moon of Saturn

# Returns

```(X, Y, Z)```

# Arguments

* ```JD```: Julian (Ephemeris) day
**/

fn ApprntRectCoords(JD: f64, moon: &Moon) -> (f64, f64, f64) {
    let mut info = CreateInfoStruct(JD - 0.04942);

    info.lambda0 = 46.170287_f64.to_radians();
    info.beta0 = -2.544441_f64.to_radians();
    info.delta = 8.557599;

    Mimas(&info)
}

struct Info {
    t1: f64,
    t2: f64,
    t3: f64,
    t4: f64,
    t5: f64,
    t6: f64,
    t7: f64,
    t8: f64,
    t9: f64,
    t10: f64,
    t11: f64,

    W0: f64,
    W1: f64,
    W2: f64,
    W3: f64,
    W4: f64,
    W5: f64,
    W6: f64,
    W7: f64,
    W8: f64,

    s1: f64,
    c1: f64,
    s2: f64,
    c2: f64,

    e1: f64,

    lambda0: f64,
    beta0: f64,
    delta: f64
}

fn CreateInfoStruct(JD: f64) -> Info {
    let mut info = Info {
        t1: 0.0, t2: 0.0, t3: 0.0, t4: 0.0, t5: 0.0,
        t6: 0.0, t7: 0.0, t8: 0.0, t9: 0.0, t10: 0.0, t11: 0.0,
        W0: 0.0, W1: 0.0, W2: 0.0, W3: 0.0, W4: 0.0,
        W5: 0.0, W6: 0.0, W7: 0.0, W8: 0.0,
        s1: 0.0, c1: 0.0, s2: 0.0, c2: 0.0,
        e1: 0.0, lambda0: 0.0, beta0: 0.0, delta: 0.0
    };

    info.t1 = (JD - 2411093.0);
    info.t2 = info.t1/365.25;
    info.t3 = (JD - 2433282.423)/365.25 + 1950.0;
    info.t4 = JD - 2411368.0;
    info.t5 = info.t4/365.25;
    info.t6 = JD - 2415020.0;
    info.t7 = info.t6/36525.0;
    info.t8 = info.t6/365.25;
    info.t9 = (JD - 2442000.5)/365.25;
    info.t10 = JD - 2409786.0;
    info.t11 = info.t10/36525.0;

    info.W0 = 5.095 * (info.t3 - 1866.39).to_radians();

	info.W1 = (74.4     + 32.39*info.t2).to_radians();
	info.W2 = (134.3    + 92.62*info.t2).to_radians();
	info.W3 = (42.0     - 0.5118*info.t5).to_radians();
	info.W4 = (276.59   + 0.5118*info.t5).to_radians();
	info.W5 = (267.2635 + 1222.1136*info.t7).to_radians();
	info.W6 = (175.4762 + 1221.5515*info.t7).to_radians();
	info.W7 = (2.4891   + 0.002435*info.t7).to_radians();
	info.W8 = (113.35   - 0.2597*info.t7).to_radians();

    info.s1 = 28.0817_f64.to_radians().sin();
    info.c1 = 28.0817_f64.to_radians().cos();
    info.s2 = 168.8112_f64.to_radians().sin();
    info.c2 = 168.8112_f64.to_radians().cos();

    info.e1 = 0.05589 - 0.000346*info.t7;

    info.lambda0 = 0.0;
    info.beta0 = 0.0;
    info.delta = 0.0;

    info
}

fn Mimas(info: &Info) -> (f64, f64, f64) {
    let L = (
        127.64
        + 381.994497*info.t1
        - 43.57*info.W0.sin()
        - 0.72*(3.0*info.W0).sin()
        - 0.02144*(5.0*info.W0).sin()
    ).to_radians();

    let p = (106.1 + 365.549*info.t2).to_radians();
    let M = L - p;
    let C = (
        2.18287*M.sin()
        + 0.025988*(2.0*M).sin()
        + 0.00043*(3.0*M).sin()
    ).to_radians();

    let lambda_1 = L + C;
    let gamma_1 = 1.563_f64.to_radians();
    let Omega_1 = (54.5 - 365.072*info.t2).to_radians();
    let r_1 = 3.06879/(1.0 + 0.01905*(M + C).cos());

    XYZ(lambda_1, gamma_1, Omega_1, r_1, &info)
}

fn XYZ(lambda_j: f64, gamma_j: f64, Omega_j: f64, r_j: f64, info: &Info) -> (f64, f64, f64) {
    let u = lambda_j - Omega_j;
    let w = Omega_j - 168.8112_f64.to_radians();

    // the moon of interest
    let X_j = r_j*(u.cos()*w.cos() - u.sin()*gamma_j.cos()*w.sin());
    let Y_j = r_j*(u.sin()*w.cos()*gamma_j.cos() + u.cos()*w.sin());
    let Z_j = r_j*u.sin()*gamma_j.sin();

    // a ficticious ninth moon
    let X_9 = 0.0;
    let Y_9 = 0.0;
    let Z_9 = 1.0;

    // some fancy stuff
    let (X9, Y9, Z9, D9) = D(X_9, Y_9, Z_9, 0.0, &info);
    let (mut X, mut Y, Z, D) = D(X_j, Y_j, Z_j, D9, &info);

    // correct for differential light-time
    let K = 20947.0;
    X += Z.abs()*(1.0 - (X/r_j).powi(2)).sqrt()/K;

    // correct for the perspective effect
    let W = info.delta/(info.delta + Z/2475.0);
    X *= W;
    Y *= W;

    (X, Y, Z)
}

// does fancy stuff and returns (X, Y, Z, D)
fn D(X_j: f64, Y_j: f64, Z_j: f64, D_j: f64, info: &Info) -> (f64, f64, f64, f64) {
    let A1 = X_j;
    let B1 = info.c1*Y_j - info.s1*Z_j;
    let C1 = info.s1*Y_j + info.c1*Z_j;

    let A2 = info.c2*A1 - info.s2*B1;
    let B2 = info.s2*A1 + info.c2*B1;
    let C2 = C1;

    let A3 = A2*info.lambda0.sin() - B2*info.lambda0.cos();
    let B3 = A2*info.lambda0.cos() + B2*info.lambda0.sin();
    let C3 = C1;

    let A4 = A3;
    let B4 = B3*info.beta0.cos() + C3*info.beta0.sin();
    let C4 = C3*info.beta0.cos() - B3*info.beta0.sin();

    let et = A4;
    let nu = C4;

    let D = et.atan2(nu);

    let X = A4*D_j.cos() - C4*D_j.sin();
    let Y = A4*D_j.sin() + C4*D_j.cos();
    let Z = B4;

    (X, Y, Z, D)
}

fn Enceladus() {

}

fn Tethys() {

}

fn Dione() {

}

fn Rhea() {

}

fn Titan() {

}

fn Hyperion() {

}

fn Iapetus() {

}
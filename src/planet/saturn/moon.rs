/*
Copyright (c) 2015, 2016 Saurav Sachidanand

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

//! Eight moons of Saturn

use crate::planet;
use crate::precess;
use crate::time;

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
Computes the apparent rectangular coordinates for a moon of Saturn

# Returns

`(X, Y, Z)`

The rectangular coordinates returned give the apparent position of a moon
with respect to Saturn as seen from Earth. The `X` and `Y`
coordinates are measured from the center of the disk of Saturn, in units
of Saturn's equatorial radius.

`X` is measured positively to the west of Saturn, and negatively to the
east. The x-axis coincides with Saturn's equator.

`Y` is measured positively to the north of Saturn, and negatively to
the south. The y-axis coincides with Saturn's axis of rotation.

`Z` only matters in sign; it is positive if the Earth-moon
distance is greater than the Earth-Saturn distance, and is negative if the
Earth-moon distance is lesser than the Earth-Saturn distance.

# Arguments

* `JD`  : Julian (Ephemeris) day
* `moon`: The [Moon](./enum.Moon.html)
**/
pub fn apprnt_rect_coords(JD: f64, moon: &Moon) -> (f64, f64, f64) {

    let mut info = create_info_struct(JD - 0.04942);

    let (planet_ecl_point, saturn_earth_dist) =
        planet::geocent_apprnt_ecl_coords(&planet::Planet::Saturn, JD);
    let (lambda0, beta0) = (planet_ecl_point.long, planet_ecl_point.lat);

    let (lambda0, beta0) = precess::precess_ecl_coords(
        lambda0, beta0,
        JD,
        time::julian_day(
            &time::Date {
                year: 1950,
                month: time::Month::Jan,
                decimal_day: 1.5,
                cal_type: time::CalType::Gregorian
            }
        )
    );

    info.lambda0 = lambda0;
    info.beta0 = beta0;
    info.delta = saturn_earth_dist;

    let (lambda_j, gamma_j, Omega_j, r_j) = match *moon {
        Moon::Mimas     => Mimas(&info),
        Moon::Enceladus => Enceladus(&info),
        Moon::Tethys    => Tethys(&info),
        Moon::Dione     => Dione(&info),
        Moon::Rhea      => Rhea(&info),
        Moon::Titan     => Titan(&info),
        Moon::Hyperion  => Hyperion(&info),
        Moon::Iapetus   => Iapetus(&info),
    };

    XYZ(lambda_j, gamma_j, Omega_j, r_j, &info, &moon)

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

#[inline(always)]
fn create_info_struct(JD: f64) -> Info {

    let angle1 = 28.0817_f64.to_radians();
    let angle2 = 168.8112_f64.to_radians();

    let mut info = Info {
        t1: JD - 2411093.0,
        t2: 0.0,
        t3: (JD - 2433282.423)/365.25 + 1950.0,
        t4: JD - 2411368.0,
        t5: 0.0,
        t6: JD - 2415020.0,
        t7: 0.0, t8: 0.0,
        t9: (JD - 2442000.5)/365.25,
        t10: JD - 2409786.0,
        t11: 0.0,
        W0: 0.0, W1: 0.0, W2: 0.0, W3: 0.0, W4: 0.0,
        W5: 0.0, W6: 0.0, W7: 0.0, W8: 0.0,
        s1: angle1.sin(),
        c1: angle1.cos(),
        s2: angle2.sin(),
        c2: angle2.cos(),
        e1: 0.0, lambda0: 0.0, beta0: 0.0, delta: 0.0
    };

    info.t2 = info.t1/365.25;
    info.t5 = info.t4/365.25;
    info.t7 = info.t6/36525.0;
    info.t8 = info.t6/365.25;
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

    info.e1 = 0.05589 - 0.000346*info.t7;

    info

}

#[inline(always)]
fn Mimas(info: &Info) -> (f64, f64, f64, f64) {

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

    (lambda_1, gamma_1, Omega_1, r_1)

}

#[inline(always)]
fn Enceladus(info: &Info) -> (f64, f64, f64, f64) {

    let L = (
        200.317
        + 262.7319002*info.t1
        + 0.25667*info.W1.sin()
        + 0.20883*info.W2.sin()
    ).to_radians();

    let p = (309.107 + 123.44121*info.t2).to_radians();
    let M = L - p;
    let C = (
        0.55577*M.sin()
        + 0.00168*(2.0*M).sin()
    ).to_radians();

    let lambda_2 = L + C;
    let gamma_2 = 0.0262_f64.to_radians();
    let Omega_2 = (348.0 - 151.95*info.t2).to_radians();
    let r_2 = 3.94118/(1.0 + 0.00485*(M + C).cos());

    (lambda_2, gamma_2, Omega_2, r_2)

}

#[inline(always)]
fn Tethys(info: &Info) -> (f64, f64, f64, f64) {

    let lambda_3 = (
        285.306
        + 190.69791226*info.t1
        + 2.063*info.W0.sin()
        + 0.03409*(3.0*info.W0).sin()
        + 0.001015*(5.0*info.W0).sin()
    ).to_radians();
    let gamma_3 = 1.0976_f64.to_radians();
    let Omega_3 =(111.33 - 72.2441*info.t2).to_radians();
    let r_3 = 4.880998;

    (lambda_3, gamma_3, Omega_3, r_3)

}

#[inline(always)]
fn Dione(info: &Info) -> (f64, f64, f64, f64) {

    let L = (
        254.712
        + 131.53493193*info.t1
        - 0.0215*info.W1.sin()
        - 0.01733*info.W2.sin()
    ).to_radians();

    let p = (174.8 + 30.82*info.t2).to_radians();
    let M = L - p;
    let C = (
        0.24717*M.sin()
        + 0.00033*(2.0*M).sin()
    ).to_radians();

    let lambda_4 = L + C;
    let gamma_4 = 0.0139_f64.to_radians();
    let Omega_4 = (232.0 - 30.27*info.t2).to_radians();
    let r_4 = 6.24871/(1.0 + 0.002157*(M + C).cos());

    (lambda_4, gamma_4, Omega_4, r_4)

}

#[inline(always)]
fn Rhea(info: &Info) -> (f64, f64, f64, f64) {

    let p1 = (342.7 + 10.057*info.t2).to_radians();
    let a1 = 0.000265*p1.sin() + 0.01*info.W4.sin();
    let a2 = 0.000265*p1.cos() + 0.01*info.W4.cos();
    let e = (a1*a1 + a2*a2).sqrt();
    let p = a1.atan2(a2);
    let N = (345.0 - 10.057*info.t2).to_radians();
    let lambda1 = (359.244 + 79.6900472*info.t1 + 0.086754*N.sin()).to_radians();
    let i = (28.0362 + 0.346898*N.cos() + 0.0193*info.W3.cos()).to_radians();
    let Omega = (168.8034 + 0.736936*N.sin() + 0.041*info.W3.sin()).to_radians();
    let a = 8.725924;

    funroutine(e, a, Omega, i, lambda1, p, &info)

}

#[inline(always)]
fn Titan(info: &Info) -> (f64, f64, f64, f64) {

    let L = (261.1582 + 22.57697855*info.t4 + 0.074025*info.W3.sin()).to_radians();
    let i1 = (27.45141 + 0.295999*info.W3.cos()).to_radians();
    let Omega1 = (168.66925 + 0.628808*info.W3.sin()).to_radians();
    let a1 = info.W7.sin()*(Omega1 - info.W8).sin();
    let a2 = info.W7.cos()*i1.sin()
             - info.W7.sin()*i1.cos()*(Omega1 - info.W8).cos();
    let g0 = 102.8623_f64.to_radians();
    let phi = a1.atan2(a2);
    let s = (a1*a1 + a2*a2).sqrt();
    let mut g = info.W4 - Omega1 - phi;
    let mut w_dash = 0.0;

    let mut counter: u8 = 1;
    while counter <= 6 { // Meeus says 3 iterations are always sufficient nough, but let's just be on the safer side
        w_dash = info.W4 + 0.37515_f64.to_radians()*((2.0*g).sin() - (2.0*g0).sin());
        g = w_dash - Omega1 - phi;
        counter += 1;
    }

    let e1 = 0.029092 + 0.00019048*((2.0*g).cos() - (2.0*g0).cos());
    let q = 2.0*(info.W5 - w_dash);
    let b1 = i1.sin()*(Omega1 - info.W8).sin();
    let b2 = info.W7.cos()*i1.sin()*(Omega1 - info.W8).cos()
             - info.W7.sin()*i1.cos();
    let theta = b1.atan2(b2) + info.W8;
    let e = e1*(1.0 + 0.002778797*q.cos());
    let p = w_dash + 0.159215_f64.to_radians()*q.sin();
    let u = 2.0*(info.W5 - theta) + phi;
    let h = 0.9375*e1*e1*q.sin() + 0.1875*s*s*(2.0*(info.W5 - theta)).sin();
    let lambda1 = L - 0.254744_f64.to_radians()*(info.e1*(info.W6.sin() + 0.75*info.e1*(2.0*info.W6).sin()) + h);
    let i = i1 + 0.031843_f64.to_radians()*s*u.cos();
    let Omega = Omega1 + 0.031843_f64.to_radians()*s*u.sin()/i1.sin();
    let a = 20.216193;

    funroutine(e, a, Omega, i, lambda1, p, &info)

}

#[inline(always)]
fn Hyperion(info: &Info) -> (f64, f64, f64, f64) {

    let nu = (92.39 + 0.5621071*info.t6).to_radians();
    let et = (148.19 - 19.18*info.t8).to_radians();
    let theta = (184.8 - 35.41*info.t9).to_radians();
    let theta1 = theta - 7.5_f64.to_radians();
    let a_s = (176.0 + 12.22*info.t8).to_radians();
    let b_s = (8.0 + 24.44*info.t8).to_radians();
    let c_s = b_s + 5_f64.to_radians();
    let w_dash = (69.898 - 18.67088*info.t8).to_radians();
    let phi = 2.0*(w_dash - info.W5);
    let xi = (94.9 - 2.292*info.t8).to_radians();
    let a = 24.50601
            - 0.08686*nu.cos()
            - 0.00166*(et + nu).cos()
            + 0.00175*(et - nu).cos();
    let e = 0.103458
            - 0.004099*nu.cos()
            - 0.000167*(et + nu).cos()
            + 0.000235*(et - nu).cos()
            + 0.02303*et.cos()
            - 0.00212*(2.0*et).cos()
            + 0.000151*(3.0*et).cos()
            + 0.00013*phi.cos();
    let p = w_dash + (
        0.15648*xi.sin()
        - 0.4457*nu.sin()
        - 0.2657*(et + nu).sin()
        - 0.3573*(et - nu).sin()
        - 12.872*et.sin()
        + 1.668*(2.0*et).sin()
        - 0.2419*(3.0*et).sin()
        - 0.07*phi.sin()
    ).to_radians();
    let lambda1 = (
        177.047
        + 16.91993829*info.t6
        + 0.15648*xi.sin()
        + 9.142*nu.sin()
        + 0.007*(2.0*nu).sin()
        - 0.014*(3.0*nu).sin()
        + 0.2275*(et + nu).sin()
        + 0.2112*(et - nu).sin()
        - 0.26*et.sin()
        - 0.0098*(2.0*et).sin()
        - 0.013*a_s.sin()
        + 0.017*b_s.sin()
        - 0.0303*phi.sin()
    ).to_radians();
    let i = (
        27.3347
        + 0.643486*xi.cos()
        + 0.315*info.W3.cos()
        + 0.018*(theta.cos() - c_s.cos())
    ).to_radians();
    let Omega = (
        168.6812
        + 1.40136*xi.cos()
        + 0.68599*info.W3.sin()
        - 0.0392*c_s.sin()
        + 0.0366*theta1.sin()
    ).to_radians();

    funroutine(e, a, Omega, i, lambda1, p, &info)

}

#[inline(always)]
fn Iapetus(info: &Info) -> (f64, f64, f64, f64) {

    let L = (261.1582 + 22.57697855*info.t4).to_radians();
    let w_dash1 = (91.796 + 0.562*info.t7).to_radians();
    let Phi = (4.367 - 0.195*info.t7).to_radians();
    let theta = (146.819 - 3.198*info.t7).to_radians();
    let phi = (60.47 + 1.521*info.t7).to_radians();
    let pho = (205.055 - 2.091*info.t7).to_radians();
    let e1 = 0.028298 + 0.001156*info.t11;
    let w_dash0 = (352.91 + 11.71*info.t11).to_radians();
    let mu = (76.3852 + 4.53795125*info.t10).to_radians();
    let i1 = (
        18.4602 -
        info.t11*(0.9518 + info.t11*(0.072 - 0.0054*info.t11))
    ).to_radians();
    let Omega1 = (
        143.198 -
        info.t11*(3.919 - info.t11*(0.116 + 0.008*info.t11))
    ).to_radians();
    let l = mu - w_dash0;
    let g = w_dash0 - Omega1 - Phi;
    let g1 = w_dash0 - Omega1 - phi;
    let ls = info.W5 - w_dash1;
    let gs = w_dash1 - theta;
    let lT = L - info.W4;
    let gT = info.W4 - pho;
    let u1 = 2.0*(l + g - ls - gs);
    let u2 = l + g1 - lT - gT;
    let u3 = l + 2.0*(g - ls - gs);
    let u4 = lT + gT - g1;
    let u5 = 2.0*(ls + gs);
    let a = 58.935028 + 0.004638*u1.cos() + 0.058222*u2.cos();
    let e =
        e1 - 0.0014097*(g1 - gT).cos()
        + 0.0003733*(u5 - 2.0*g).cos()
        + 0.000118*u3.cos() + 0.0002408*l.cos()
        + 0.0002849*(l + u2).cos()
        + 0.000619*u4.cos();
    let w = (
        0.08077*(g1 - gT).sin()
        + 0.02139*(u5 - 2.0*g).sin()
        - 0.00676*u3.sin()
        + 0.0138*l.sin()
        + 0.01632*(l + u2).sin()
        + 0.03547*u4.sin()
    ).to_radians();
    let p = w_dash0 + w/e1;
    let lambda1 = mu + (
        - 0.04299*u2.sin()
        - 0.00789*u1.sin()
        - 0.06312*ls.sin()
        - 0.00295*(2.0*ls).sin()
        - 0.02231*u5.sin()
        + 0.0065*(u5 + Phi).sin()
    ).to_radians();
    let i = i1 + (
        0.04204*(u5 + Phi).cos()
        + 0.00235*(l + g1 + lT + gT + phi).cos()
        + 0.0036*(u2 + phi).cos()
    ).to_radians();
    let w1 = (
        0.04204*(u5 + Phi).sin()
        + 0.00235*(l + g1 + lT + gT + phi).sin()
        + 0.00358*(u2 + phi).sin()
    ).to_radians();
    let Omega = Omega1 + w1/i1.sin();

    funroutine(e, a, Omega, i, lambda1, p, &info)

}

fn funroutine(e: f64, a: f64, Omega: f64, i: f64, lambda1: f64, p: f64, info: &Info) -> (f64, f64, f64, f64) {

    let M = lambda1 - p;
    let C =   e*((2.0 - e*e*(0.25 - 0.0520833333*e*e))*M.sin()
                 + e*((1.25 - 0.458333333*e*e)*(2.0*M).sin()
                      + e*((1.083333333 - 0.671875*e*e)*(3.0*M).sin()
                           + e*(1.072917*(4.0*M).sin()
                               + e*1.142708*(5.0*M).sin()))));
    let r = a*(1.0 - e*e)/(1.0 + e*(M + C).cos());
    let g = Omega - 168.8112_f64.to_radians();
    let a1 = i.sin()*g.sin();
    let a2 = info.c1*i.sin()*g.cos() - info.s1*i.cos();
    let gamma = (a1*a1 + a2*a2).sqrt().asin();
    let u = a1.atan2(a2);
    let w = 168.8112_f64.to_radians() + u;
    let h = info.c1*i.sin() - info.s1*i.cos()*g.cos();
    let phi = (info.s1*g.sin()).atan2(h);
    let lambda = lambda1 + C + u - g - phi;

    (lambda, gamma, w, r)

}

#[allow(unused_variables)]
#[inline(always)]
fn XYZ(lambda_j: f64, gamma_j: f64, Omega_j: f64, r_j: f64, info: &Info, moon: &Moon) -> (f64, f64, f64) {

    let u = lambda_j - Omega_j;
    let w = Omega_j - 168.8112_f64.to_radians();

    // moon of interest
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
    let K = match *moon {
        Moon::Mimas     => 20947.0,
        Moon::Enceladus => 23715.0,
        Moon::Tethys    => 26382.0,
        Moon::Dione     => 29876.0,
        Moon::Rhea      => 35313.0,
        Moon::Titan     => 53800.0,
        Moon::Hyperion  => 59222.0,
        Moon::Iapetus   => 91820.0,
    };
    X += Z.abs()*(1.0 - (X/r_j).powi(2)).sqrt()/K;

    // correct for the perspective effect
    let W = info.delta / (info.delta + Z/2475.0);
    X *= W;
    Y *= W;

    (X, Y, Z)

}

// does fancy stuff and computes (X, Y, Z, D)
fn D(X_j: f64, Y_j: f64, Z_j: f64, D_j: f64, info: &Info) -> (f64, f64, f64, f64) {

    let A1 = X_j;
    let B1 = info.c1*Y_j - info.s1*Z_j;
    let C1 = info.s1*Y_j + info.c1*Z_j;

    let A2 = info.c2*A1 - info.s2*B1;
    let B2 = info.s2*A1 + info.c2*B1;

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

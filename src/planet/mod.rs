#[macro_use]
mod mercury;

#[macro_use]
mod venus;

#[macro_use]
pub mod earth;

#[macro_use]
pub mod mars;

#[macro_use]
pub mod jupiter;

#[macro_use]
pub mod saturn;

#[macro_use]
mod uranus;

#[macro_use]
mod neptune;

use angle;
use time;

/// Represents a **planet**
pub enum Planet {
    /// Mercury *Thanks for help proving General Relativity*
    Mercury,
    /// Venus **Climate change was here**
    Venus,
    /// Earth *Pale blue dot*
    Earth,
    /// Mars *Panspermia sounds nice*
    Mars,
    /// Jupiter *Oh, Europa*
    Jupiter,
    /// Saturn *62 moons and smaller moonlets*
    Saturn,
    /// Uranus *Remember King George?*
    Uranus,
    /// Neptune *Oceans of liquid diamond, maybe*
    Neptune,
}

/**
Returns a planet's **geocentric equatorial semidiameter**

# Returns

* ```semidiameter```: Geocentric equatorial emidiameter *(radians per AU)*

# Arguments

* ```planet```: [Planet](./enum.Planet.html).
                *Throws an error if Planet::Earth is passed.*
* ```distance_to_earth```: Planet's distance to Earth *(AU)*
**/
pub fn Semidiameter(planet: Planet, distance_to_earth: f64) -> f64 {
    let mut s: f64;

    match planet {
        Planet::Mercury => s = angle::PureDegrees(0, 0, 3.36),
        Planet::Venus => s = angle::PureDegrees(0, 0, 8.41),
        Planet::Earth => panic!("Planet::Earth was passed to the function
                                 planet::Semidiameter(). This does not
                                 make logical sense."),
        Planet::Mars => s = angle::PureDegrees(0, 0, 4.68),
        Planet::Jupiter => s = jupiter::EquatorSemidiameter(1.0),
        Planet::Saturn => s = saturn::EquatorSemidiameter(1.0),
        Planet::Uranus => s = angle::PureDegrees(0, 0, 35.02),
        Planet::Neptune => s = angle::PureDegrees(0, 0, 33.5),
    };

    s / distance_to_earth
}

/**
Returns a planet's **orbital elements**

# Returns

```(L, a, e, i, omega, pi, M, w)```

* ```L```: Mean longitude *(radians)*
* ```a```: Semimajor axis of the orbit *(AU)*
* ```e```: Eccentricity of the orbit
* ```i```: Inclination of the plane of the orbit with the plane of
           the Earth's ecliptic *(radians)*
* ```omega```: Longitude of the ascending node *(radians)*.
               *An undefined value is returned for Planet::Earth*.
* ```pi```: Longitude of the perihelion *(radians)*
* ```M```: Mean anomaly *(radians)*
* ```w```: Argument of the perihelion *(radians)*

# Arguments

* ```planet```: [Planet](./enum.Planet.html)
* ```JD```: Julian day
**/
pub fn OrbitalElements(planet: Planet, JD: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let T = time::JulianCentury(JD);
    let TT = T * T;
    let TTT = TT * T;

    let mut L: f64;
    let mut a: f64;
    let mut e: f64;
    let mut i: f64;
    let mut omega: f64;
    let mut pi: f64;

    match planet {
        Planet::Mercury => {
            L = 252.250906 + 149474.0722491*T + 0.0003035*TT + 0.000000018*TTT;
            a = 0.038709831;
            e = 0.20563175 + 0.000020407*T - 0.0000000283*TT + 0.00000000018*TTT;
            i = 7.004986 + 0.0018215*T - 0.0000181*TT + 0.000000056*TTT;
            omega = 48.330893 + 1.1861883*T + 0.00017542*TT + 0.000000215*TTT;
            pi = 77.456119 + 1.5564776*T + 0.00029544*TT + 0.000000009*TTT;
        },
        Planet::Venus => {
            L = 181.979801 + 58519.2130302*T + 0.00031014*TT + 0.000000015*TTT;
            a = 0.72332982;
            e = 0.00677192 - 0.000047765*T + 0.0000000981*TTT + 0.00000000046*TTT;
            i = 3.394662 + 0.0010037*T - 0.00000088*TT - 0.000000007*TTT;
            omega = 76.67992 + 0.9011206*T + 0.00040618*TT - 0.000000093*TTT;
            pi = 131.563703 + 1.4022288*T - 0.00107618*TT - 0.000005678*TTT;
        },
        Planet::Earth => {
            L = 100.466457 + 36000.7698278*T + 0.00030322*TT + 0.00000002*TTT;
            a = 1.000001018;
            e = 0.01670863 - 0.000042037*T - 0.0000001267*TTT + 0.00000000014*TTT;
            i = 0.0;
            pi = 102.937348 + 1.7195366*T + 0.00045688*TT - 0.000000018*TTT;
            omega = 0.0
        },
        Planet::Mars => {
            L = 355.433 + 19141.6964471*T + 0.00031052*TT + 0.000000016*TTT;
            a = 1.523679342;
            e = 0.09340065 + 0.000090484*T - 0.0000000806*TTT - 0.00000000025*TTT;
            i = 1.849726 - 0.0006011*T + 0.00001276*TT - 0.000000007*TTT;
            omega = 49.558093 + 0.7720959*T + 0.00001557*TT - 0.000002267*TTT;
            pi = 336.060234 + 1.8410449*T + 0.00013477*TT + 0.000000536*TTT;
        },
        Planet::Jupiter => {
            L = 34.351519 + 3036.3027748*T + 0.0002233*TT + 0.000000037*TTT;
            a = 5.202603209 + 0.0000001913*T;
            e = 0.04849793 + 0.000163225*T - 0.0000004714*TTT - 0.00000000201*TTT;
            i = 1.303267 - 0.0054965*T + 0.00000466*TT - 0.000000002*TTT;
            omega = 100.464407 + 1.0209774*T + 0.00040315*TT + 0.000000404*TTT;
            pi = 14.331207 + 1.6126352*T + 0.00103042*TT - 0.000004464*TTT;
        },
        Planet::Saturn => {
            L = 50.077444 + 1223.5110686*T + 0.00051908*TT - 0.00000003*TTT;
            a = 9.554909192 - 0.0000021390*T + 0.000000004*TT;
            e = 0.05554814 - 0.000346641*T - 0.0000006436*TTT + 0.0000000034*TTT;
            i = 2.488879 - 0.0037362*T - 0.00001519*TT + 0.000000087*TTT;
            omega = 113.665503 + 0.877088*T - 0.00012176*TT - 0.000002249*TTT;
            pi = 93.057237 + 1.9637613*T + 0.00083753*TT + 0.000004928*TTT;
        },
        Planet::Uranus => {
            L = 314.055005 + 429.8640561*T + 0.0003039*TT - 0.000000026*TTT;
            a = 19.218446062 - 0.0000000372*T + 0.00000000098*TT;
            e = 0.04638122 - 0.000027293*T + 0.0000000789*TTT + 0.00000000024*TTT;
            i = 0.773197 + 0.0007744*T + 0.00003749*TT - 0.000000092*TTT;
            omega = 74.005957 + 0.5211278*T + 0.00133947*TT + 0.000018484*TTT;
            pi = 173.005291 + 1.486379*T + 0.00021406*TT + 0.000000434*TTT;
        },
        Planet::Neptune => {
            L = 304.348665 + 219.8833092*T + 0.00030882*TT + 0.000000018*TTT;
            a = 30.110386869 - 0.0000001663*T + 0.00000000069*TT;
            e = 0.00945575 + 0.000006033*T - 0.00000000005*TTT;
            i = 1.769953 - 0.0093082*T - 0.00000708*TT + 0.000000027*TTT;
            omega = 131.784057 + 1.1022039*T + 0.00025952*TT - 0.000000637*TTT;
            pi = 48.120276 + 1.4262957*T + 0.00038434*TT + 0.00000002*TTT;
        },
    };

    (angle::LimitedTo360(L).to_radians(),
     a, e,
     angle::LimitedTo360(i).to_radians(),
     angle::LimitedTo360(omega).to_radians(),
     angle::LimitedTo360(pi).to_radians(),
     angle::LimitedTo360(L - pi).to_radians(),
     angle::LimitedTo360(pi - omega).to_radians()
    )
}

pub fn HeliocentricCoords(planet: Planet, JD: f64) -> (f64, f64, f64) {/*
    let terms = match planet {
        Planet::Mercury => VSOP87_Mercury_Terms!(),
        Planet::Venus => VSOP87_Mercury_Terms!(),
        Planet::Earth => VSOP87_Mercury_Terms!(),
        Planet::Mars => VSOP87_Mercury_Terms!(),
        Planet::Jupiter => VSOP87_Mercury_Terms!(),
        Planet::Saturn => VSOP87_Mercury_Terms!(),
        Planet::Uranus => VSOP87_Mercury_Terms!(),
        Planet::Neptune => VSOP87_Mercury_Terms!(),
    };

    let JM = time::JulianMillenium(JD);

    let (L0_terms, L1_terms, L2_terms, L3_terms, L4_terms, L5_terms) = terms.0;
    let L = VSOP87Coordinate(JM, &L0_terms, &L1_terms, &L2_terms, &L3_terms, &L4_terms, &L5_terms);

    let (B0_terms, B1_terms, B2_terms, B3_terms, B4_terms, B5_terms) = terms.1;
    let B = VSOP87Coordinate(JM, &B0_terms, &B1_terms, &B2_terms, &B3_terms, &B4_terms, &B5_terms);

    let (R0_terms, R1_terms, R2_terms, R3_terms, R4_terms, R5_terms) = terms.2;
    let R = VSOP87Coordinate(JM, &R0_terms, &R1_terms, &R2_terms, &R3_terms, &R4_terms, &R5_terms);

    (L, B, R)*/(0.0, 0.0, 0.0)
}

fn VSOP87Coordinate(t: f64, a: &[[f64; 3]], b: &[[f64; 3]], c: &[[f64; 3]], d: &[[f64; 3]], e: &[[f64; 3]], f: &[[f64; 3]]) -> f64 {
    let mut T0 = 0.0; for &i in a.iter() { T0 += VSOP87Term(t, &i); }
    let mut T1 = 0.0; for &i in b.iter() { T1 += VSOP87Term(t, &i); }
    let mut T2 = 0.0; for &i in c.iter() { T2 += VSOP87Term(t, &i); }
    let mut T3 = 0.0; for &i in d.iter() { T3 += VSOP87Term(t, &i); }
    let mut T4 = 0.0; for &i in e.iter() { T4 += VSOP87Term(t, &i); }
    let mut T5 = 0.0; for &i in f.iter() { T5 += VSOP87Term(t, &i); }

    T0 +
    t * (T1 +
    t * (T2 +
    t * (T3 +
    t * (T4 +
    t * T5 ))))
}

fn VSOP87Term(t: f64, array: &[f64]) -> f64 {
    array[0] * (array[1] + t*array[2]).cos()
}

fn EffectOfLightTime(x: f64, y: f64, z: f64) -> f64 {
    0.0057755183 * (x*x + y*y + z*z).sqrt()
}

pub fn GeocentricEclipticalCoords(L: f64, B: f64, R: f64, L0: f64, B0: f64, R0: f64) -> (f64, f64, f64) {
    let x = R*B.cos()*L.cos() - R0*B0.cos()*L0.cos();
    let y = R*B.cos()*L.sin() - R0*B0.cos()*L0.sin();
    let z = R*B.sin() - R0*B0.sin();

    (y.atan2(x),
     z/(x*x + y*y).sqrt(),
     EffectOfLightTime(x, y, z))
}

pub fn GeocentricEquatorialCoords(X: f64, Y: f64, Z: f64, semimaj_axis: f64, e: f64, i: f64, w: f64, sigma: f64, n: f64,
            oblq_eclip: f64, M: f64, E: f64, v: f64, r: f64) -> (f64, f64, f64) {

    let F = sigma.cos();
    let G = sigma.sin() * oblq_eclip.cos();
    let H = sigma.sin() * oblq_eclip.sin();

    let P = - sigma.sin() * i.cos();
    let Q =   sigma.cos() * i.cos() * oblq_eclip.cos()
            - i.sin()     * oblq_eclip.sin();
    let R =   sigma.cos() * i.cos() * oblq_eclip.sin()
            + i.sin()     * oblq_eclip.cos();

    let A = F.atan2(P);
    let B = G.atan2(Q);
    let C = H.atan2(R);
    let a = (F*F + P*P).sqrt();
    let b = (G*G + Q*Q).sqrt();
    let c = (H*H + R*R).sqrt();

    let x = r * a * (A + w + v);
    let y = r * b * (B + w + v);
    let z = r * c * (C + w + v);

    let xi = X + x;
    let nu = Y + y;
    let et = Z + z;

    let asc = angle::LimitedTo360(nu.atan2(xi).to_degrees()).to_radians();
    let dec = et.atan2((xi*xi + nu*nu).sqrt());

    (asc, dec, EffectOfLightTime(x, y, z))
}

pub fn HeliocentricCoordsFromOrbitalElements(i: f64, sigma: f64, w: f64, v: f64, r: f64) -> (f64, f64) {
    let u = w + v;
    let x = r * (sigma.cos()*u.cos() - sigma.sin()*u.sin()*i.cos());
    let y = r * (sigma.sin()*u.cos() + sigma.cos()*u.sin()*i.cos());
    let z = r * i.sin() * u.sin();

    (y.atan2(x),
     z.atan2((x*x + y*y).sqrt()))
}

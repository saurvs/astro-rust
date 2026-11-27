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

//! 8 Planets of the Solar System

mod VSOPD_87;

pub mod earth;
pub mod mars;
pub mod jupiter;
pub mod saturn;

use crate::angle;
use crate::coords;
use crate::time;

/// Represents a planet
pub enum Planet {
    /// Mercury *Helped with testing General Relativity*
    Mercury,
    /// Venus **Climate change was here**
    Venus,
    /// Earth *Pale blue dot, right?*
    Earth,
    /// Mars *Panspermia sounds nice*
    Mars,
    /// Jupiter *Oh, Europa*
    Jupiter,
    /// Saturn *62 moons, and some moonlets*
    Saturn,
    /// Uranus *Remember King George?*
    Uranus,
    /// Neptune *Oceans of liquid diamond (probably)*
    Neptune,
}

/**
Computes the illuminated fraction of a planet's disk from it's phase
angle

# Returns

* `illum_frac`: Illuminated fraction of the planet's disk

# Arguments

* `i`: Phase angle of the planet *| in radians*
**/
#[inline]
pub fn illum_frac_frm_phase_angl(i: f64) -> f64 {

    (1.0 + i.cos()) / 2.0

}

/**
Computes the illuminated fraction of a planet's disk from it's distance
to the Sun and Earth

# Returns

* `illum_frac`: Illuminated fraction of the planet's disk

# Arguments

* `r`    : Planet-Sun distance *| in AU*
* `delta`: Planet-Earth distance *| in AU*
* `R`    : Sun-Earth distance *| in AU*
**/
#[inline]
pub fn illum_frac_frm_dist(r: f64, delta: f64, R: f64) -> f64 {

    let x = r + delta;

    (x*x - R*R) / (4.0 * r * delta)

}

/**
Computes a planet's phase angle

# Returns

* `phase_angl`: Phase angle of the planet *| in radians*

# Arguments

* `r`    : Planet-Sun distance *| in AU*
* `delta`: Planet-Earth distance *| in AU*
* `R`    : Sun-Earth distance *| in AU*
**/
#[inline]
pub fn phase_angl(r: f64, delta: f64, R: f64) -> f64 {

    (r*r + delta*delta - R*R) / (2.0 * r * delta)

}

/**
Computes the position angle of the bright limb of a planet

# Returns

* `pos_angl_of_bright_limb`: The position angle of the midpoint of
                             the illuminated limb of a planet
                             *| in radians*

# Arguments

* `sun_eq_point`   : Equatorial point of the Sun *| in radians*
* `planet_eq_point`: Equatorial point of the planet *| in radians*
**/
pub fn pos_angle_of_bright_limb (

    sun_eq_point    : coords::EqPoint,
    planet_eq_point : coords::EqPoint

) -> f64 {

    let a = sun_eq_point.dec.cos();
    let n = a * (sun_eq_point.asc - planet_eq_point.asc).sin();
    let d =
        sun_eq_point.dec.sin() * planet_eq_point.dec.cos()
      - planet_eq_point.dec.sin() * (sun_eq_point.asc - planet_eq_point.asc).cos() * a;

    n.atan2(d)

}

/**
Computes a planet's equatorial semidiameter

# Returns

* `semidiameter`: Equatorial semidiameter *| in radians*

# Arguments

* `planet`           : The [Planet](./enum.Planet.html)
* `planet_earth_dist`: Planet-Earth distance *| in AU*
**/
pub fn semidiameter<'a> (

    planet            : &Planet,
    planet_earth_dist : f64

) -> Result<f64, &'a str> {

    let s = match *planet {
        Planet::Mercury => angle::deg_frm_dms(0, 0, 3.360).to_radians(),
        Planet::Venus   => angle::deg_frm_dms(0, 0, 8.410).to_radians(),
        Planet::Mars    => angle::deg_frm_dms(0, 0, 4.680).to_radians(),
        Planet::Uranus  => angle::deg_frm_dms(0, 0, 35.02).to_radians(),
        Planet::Neptune => angle::deg_frm_dms(0, 0, 33.50).to_radians(),

        Planet::Jupiter => jupiter::eq_semidiameter(1.0),
        Planet::Saturn  => saturn::eq_semidiameter(1.0),

        Planet::Earth   => {
            return Err("Planet::Earth was passed to the function planet::semidiameter()");
        }
    };

    Ok( s / planet_earth_dist )

}

/**
Computes a planet's orbital elements, referred to the mean equinox of
the date

# Returns

`(L, a, e, i, omega, pi, M, w)`

* `L`: Mean longitude *| in radians*
* `a`: Semimajor axis of the orbit *| in AU*
* `e`: Eccentricity of the orbit
* `i`: Inclination of the plane of the orbit with the plane of
           the Earth's ecliptic *| in radians*
* `omega`: Longitude of the ascending node *| in radians*.
               An undefined value is returned for `Planet::Earth`.
* `pi`: Longitude of the perihelion *| in radians*
* `M`: Mean anomaly *| in radians*
* `w`: Argument of the perihelion *| in radians*.
            An undefined value is returned for `Planet::Earth`.

# Arguments

* `planet`: Any variant of [Planet](./enum.Planet.html)
* `JD`    : Julian (Ephemeris) day
**/
pub fn orb_elements(planet: &Planet, JD: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {

    let T = time::julian_cent(JD);
    let TT = T * T;
    let TTT = TT * T;

    let L;
    let a;
    let e;
    let i;
    let omega;
    let pi;

    match planet {

        &Planet::Mercury => {
            L = 252.250906 + 149474.0722491*T + 0.0003035*TT + 0.000000018*TTT;
            a = 0.038709831;
            e = 0.20563175 + 0.000020407*T - 0.0000000283*TT + 0.00000000018*TTT;
            i = 7.004986 + 0.0018215*T - 0.0000181*TT + 0.000000056*TTT;
            omega = 48.330893 + 1.1861883*T + 0.00017542*TT + 0.000000215*TTT;
            pi = 77.456119 + 1.5564776*T + 0.00029544*TT + 0.000000009*TTT;
        },

        &Planet::Venus => {
            L = 181.979801 + 58519.2130302*T + 0.00031014*TT + 0.000000015*TTT;
            a = 0.72332982;
            e = 0.00677192 - 0.000047765*T + 0.0000000981*TTT + 0.00000000046*TTT;
            i = 3.394662 + 0.0010037*T - 0.00000088*TT - 0.000000007*TTT;
            omega = 76.67992 + 0.9011206*T + 0.00040618*TT - 0.000000093*TTT;
            pi = 131.563703 + 1.4022288*T - 0.00107618*TT - 0.000005678*TTT;
        },

        &Planet::Earth => {
            L = 100.466457 + 36000.7698278*T + 0.00030322*TT + 0.00000002*TTT;
            a = 1.000001018;
            e = 0.01670863 - 0.000042037*T - 0.0000001267*TTT + 0.00000000014*TTT;
            i = 0.0;
            pi = 102.937348 + 1.7195366*T + 0.00045688*TT - 0.000000018*TTT;
            omega = 0.0
        },

        &Planet::Mars => {
            L = 355.433 + 19141.6964471*T + 0.00031052*TT + 0.000000016*TTT;
            a = 1.523679342;
            e = 0.09340065 + 0.000090484*T - 0.0000000806*TTT - 0.00000000025*TTT;
            i = 1.849726 - 0.0006011*T + 0.00001276*TT - 0.000000007*TTT;
            omega = 49.558093 + 0.7720959*T + 0.00001557*TT - 0.000002267*TTT;
            pi = 336.060234 + 1.8410449*T + 0.00013477*TT + 0.000000536*TTT;
        },

        &Planet::Jupiter => {
            L = 34.351519 + 3036.3027748*T + 0.0002233*TT + 0.000000037*TTT;
            a = 5.202603209 + 0.0000001913*T;
            e = 0.04849793 + 0.000163225*T - 0.0000004714*TTT - 0.00000000201*TTT;
            i = 1.303267 - 0.0054965*T + 0.00000466*TT - 0.000000002*TTT;
            omega = 100.464407 + 1.0209774*T + 0.00040315*TT + 0.000000404*TTT;
            pi = 14.331207 + 1.6126352*T + 0.00103042*TT - 0.000004464*TTT;
        },

        &Planet::Saturn => {
            L = 50.077444 + 1223.5110686*T + 0.00051908*TT - 0.00000003*TTT;
            a = 9.554909192 - 0.0000021390*T + 0.000000004*TT;
            e = 0.05554814 - 0.000346641*T - 0.0000006436*TTT + 0.0000000034*TTT;
            i = 2.488879 - 0.0037362*T - 0.00001519*TT + 0.000000087*TTT;
            omega = 113.665503 + 0.877088*T - 0.00012176*TT - 0.000002249*TTT;
            pi = 93.057237 + 1.9637613*T + 0.00083753*TT + 0.000004928*TTT;
        },

        &Planet::Uranus => {
            L = 314.055005 + 429.8640561*T + 0.0003039*TT - 0.000000026*TTT;
            a = 19.218446062 - 0.0000000372*T + 0.00000000098*TT;
            e = 0.04638122 - 0.000027293*T + 0.0000000789*TTT + 0.00000000024*TTT;
            i = 0.773197 + 0.0007744*T + 0.00003749*TT - 0.000000092*TTT;
            omega = 74.005957 + 0.5211278*T + 0.00133947*TT + 0.000018484*TTT;
            pi = 173.005291 + 1.486379*T + 0.00021406*TT + 0.000000434*TTT;
        },

        &Planet::Neptune => {
            L = 304.348665 + 219.8833092*T + 0.00030882*TT + 0.000000018*TTT;
            a = 30.110386869 - 0.0000001663*T + 0.00000000069*TT;
            e = 0.00945575 + 0.000006033*T - 0.00000000005*TTT;
            i = 1.769953 - 0.0093082*T - 0.00000708*TT + 0.000000027*TTT;
            omega = 131.784057 + 1.1022039*T + 0.00025952*TT - 0.000000637*TTT;
            pi = 48.120276 + 1.4262957*T + 0.00038434*TT + 0.00000002*TTT;
        },

    };

    (
        angle::limit_to_360(L).to_radians(),
        a,
        e,
        angle::limit_to_360(i).to_radians(),
        angle::limit_to_360(omega).to_radians(),
        angle::limit_to_360(pi).to_radians(),
        angle::limit_to_360(L - pi).to_radians(),
        angle::limit_to_360(pi - omega).to_radians()
    )

}

/**
Computes a planet's heliocentric coordinates, referred to the mean
equinox of the date

# Returns

`(long, lat, rad_vec)`

* `long`   : Heliocentric longitude *| in radians*
* `lat`    : Heliocentric latitude *| in radians*
* `rad_vec`: Heliocentric radius vector *| in AU*

# Arguments

* `planet`: Any variant of [Planet](./enum.Planet.html)
* `JD`    : Julian (Ephemeris) day
**/
pub fn heliocent_coords(planet: &Planet, JD: f64) -> (f64, f64, f64) {

    let VSOPD87_Terms = match planet {
        &Planet::Mercury  => VSOPD_87::mercury::terms(),
        &Planet::Venus    => VSOPD_87::venus::terms(),
        &Planet::Earth    => VSOPD_87::earth::terms(),
        &Planet::Mars     => VSOPD_87::mars::terms(),
        &Planet::Jupiter  => VSOPD_87::jupiter::terms(),
        &Planet::Saturn   => VSOPD_87::saturn::terms(),
        &Planet::Uranus   => VSOPD_87::uranus::terms(),
        &Planet::Neptune  => VSOPD_87::neptune::terms(),
    };

    let mut L = 0.0;
    let mut B = 0.0;
    let mut R = 0.0;

    let JM = time::julian_mill(JD);

    let mut n: u8 = 1; // L, then B, then R
    for i in VSOPD87_Terms.iter() { // L or B or R

        let mut T = 1.0;
        let mut y = 0.0;

        for j in i.iter() { // T or T**2 or T**3 or ...

            for k in j.iter() { // add [A * cos(B + C*T)]
                y += k[0] * (k[1] + k[2]*JM).cos();
            }

            if n == 1 {
                L += y * T;
            } else if n == 2 {
                B += y * T;
            } else if n == 3 {
                R += y * T;
            }

            y = 0.0;
            T *= JM;

        }

        n += 1;

    }

    L = angle::limit_to_two_PI(L);
    B = angle::limit_to_two_PI(B);

    (L, B, R)

}

#[inline(always)]
fn light_time(dist: f64) -> f64 {

    0.0057755183 * dist

}

/**
Computes a planet's geocentric, geometric ecliptic position,
uncorrected for light-time

# Returns

`(ecl_long, ecl_lat, rad_vec, light_time)`

* `ecl_long`: Geometric longitude of the planet *| in radians*
* `ecl_lat`: Geometric latitude of the planet *| in radians*
* `rad_vec`: Geometric radius vector of the planet *| in AU*
* `light_time`: Time taken by light to travel to the Earth
                    from the planet's current position, in days of
                    dynamical time

The coordinates returned here refer to the true position
of the planet at the time of interest, and therefore
are not corrected for the effect of light-time.

# Arguments

* `L0`: Heliocentric longitude of the Earth *| in radians*
* `B0`: Heliocentric latitude of the Earth *| in radians*
* `R0`: Heliocentric radius vector of the Earth *| in radians*
* `L` : Heliocentric longitude of the planet *| in radians*
* `B` : Heliocentric latitude of the planet *| in radians*
* `R` : Heliocentric radius vector of the planet *| in radians*
**/
pub fn geocent_geomet_ecl_coords (

    L0 : f64, B0 : f64, R0 : f64,
    L  : f64, B  : f64, R  : f64

) -> (f64, f64, f64, f64) {

    let (x, y, z) = geocent_ecl_rect_coords(L0, B0, R0, L, B, R);

    let (lambda, beta) = ecl_coords_frm_ecl_rect_coords(x, y, z);
    let planet_earth_dist = dist_frm_ecl_rect_coords(x, y, z);
    let light_time = light_time(planet_earth_dist);

    (lambda, beta, planet_earth_dist, light_time)

}

/**
Computes a planet's geocentric ecliptic rectangular coordinates from
it's heliocentric position

# Returns

`(X, Y, Z)`

# Arguments

* `L0`: Heliocentric longitude of the Earth *| in radians*
* `B0`: Heliocentric latitude of the Earth *| in radians*
* `R0`: Heliocentric radius vector of the Earth *| in radians*
* `L` : Heliocentric longitude of the planet *| in radians*
* `B` : Heliocentric latitude of the planet *| in radians*
* `R` : Heliocentric radius vector of the planet *| in radians*
**/
fn geocent_ecl_rect_coords (

    L0 : f64, B0 : f64, R0 : f64,
    L  : f64, B  : f64, R  : f64

) -> (f64, f64, f64) {

    let x = R*B.cos()*L.cos() - R0*B0.cos()*L0.cos();
    let y = R*B.cos()*L.sin() - R0*B0.cos()*L0.sin();
    let z = R*B.sin()         - R0*B0.sin();

    (x, y, z)

}

/**
Computes a planet's ecliptic coordinates, from it's geocentric
ecliptic rectangular coordinates

# Returns

* `ecl_long`: Ecliptic longitude of the planet *| in radians*
* `ecl_lat` : Ecliptic latitude of the planet *| in radians*

# Arguments

* `X`
* `Y`
* `Z`
**/
#[inline]
fn ecl_coords_frm_ecl_rect_coords(x: f64, y: f64, z: f64) -> (f64, f64) {

    (
        y.atan2(x),
        z.atan2((x*x + y*y).sqrt())
    )

}

/**
Computes a planet's distance to Earth, from it's geocentric ecliptic
rectangular coordinates

# Returns

* `planet_earth_dist`: Planet-Earth distance *| in AU*

# Arguments

* `X`
* `Y`
* `Z`
**/
#[inline]
fn dist_frm_ecl_rect_coords(x: f64, y: f64, z: f64) -> f64 {

    (x*x + y*y + z*z).sqrt()

}

/**
Computes a planet's geocentric, apparent ecliptic position, corrected
for light-time

# Returns

`(ecl_long, ecl_lat, rad_vec)`

* `planet_ecl_point`: Ecliptic point of the planet *| in radians*
* `rad_vec`         : Geocentric radius vector of the planet *| in AU*

The coordinates returned here refer to the apparent position (from Earth)
of the planet at the time of interest by correcting the true
coordinates for the effect of light-time.

# Arguments

* `planet`: Any variant of [Planet](./enum.Planet.html)
* `JD`    : Julian (Ephemeris) day
**/
#[allow(unused_variables)]
pub fn geocent_apprnt_ecl_coords(planet: &Planet, JD: f64) -> (coords::EclPoint, f64) {

    let (L0, B0, R0) = heliocent_coords(&Planet::Earth, JD);

    let (L1, B1, R1) = heliocent_coords(&planet, JD);
    let (l1, b1, r1, t) = geocent_geomet_ecl_coords(L0, B0, R0, L1, B1, R1);

    let (L2, B2, R2) = heliocent_coords(&planet, JD - t);
    let (l2, b2, r2, t2) = geocent_geomet_ecl_coords(L0, B0, R0, L2, B2, R2);

    let ecl_point = coords::EclPoint {
        long:l2,
        lat: b2
    };

    (ecl_point, r2)

}

/**
Computes a planet's geocentric ecliptic coordinates converted to the
FK5 system

# Returns

`(ecl_long_FK5, ecl_lat_FK5)`

* `ecl_long_FK5`: Ecliptic longitude of the planet, converted to the
                  FK5 system *| in radians*
* `ecl_lat_FK5` : Ecliptic latitude of the planet, converted to the
                  FK5 system *| in radians*

# Arguments

* `JD`      : Julian (Ephemeris) day
* `ecl_long`: Ecliptic longitude of the planet on `JD` referred to
              the mean equinox of the date *| in radians*
* `ecl_lat` : Ecliptic latitude of the planet on `JD`, referred to
              the mean equinox of the date *| in radians*
**/
pub fn ecl_coords_to_FK5(JD: f64, ecl_long: f64, ecl_lat: f64) -> (f64, f64) {

    let JC = time::julian_cent(JD);
    let lambda1 = ecl_long - JC*(1.397 + JC*0.00031).to_radians();
    let x = angle::deg_frm_dms(0, 0, 0.03916).to_radians();

    let ecl_long_correction = - angle::deg_frm_dms(0, 0, 0.09033).to_radians()
                              + x*(lambda1.cos() + lambda1.sin())*ecl_lat.tan();

    (
        ecl_long + ecl_long_correction,
        ecl_lat  + x*(lambda1.cos() - lambda1.sin())
    )

}

pub fn geocent_eq_coords (

    X          : f64,
    Y          : f64,
    Z          : f64,
    i          : f64,
    w          : f64,
    sigma      : f64,
    oblq_eclip : f64,
    v          : f64,
    r          : f64

) -> (f64, f64, f64) {

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

    let asc = angle::limit_to_two_PI( nu.atan2(xi) );
    let dec = et.atan2( (xi*xi + nu*nu).sqrt() );
    let dist = (x*x + y*y + z*z).sqrt();

    (asc, dec, light_time(dist))

}

pub fn heliocent_coords_frm_orb_elements(i: f64, sigma: f64, w: f64, v: f64, r: f64) -> (f64, f64) {

    let u = w + v;
    let x = r * (sigma.cos()*u.cos() - sigma.sin()*u.sin()*i.cos());
    let y = r * (sigma.sin()*u.cos() + sigma.cos()*u.sin()*i.cos());
    let z = r * i.sin() * u.sin();

    (y.atan2(x), z.atan2((x*x + y*y).sqrt()))

}

/**
Computes a planet's apparent magnitude using G. Muller's formulae

# Returns

* `app_mag`: Apparent magnitude of the planet

# Arguments

* `planet`: Any variant of [Planet](./enum.Planet.html)
* `i`     : Phase angle of the planet *| in radians*
* `delta` : Planet-Earth distance *| in AU*
* `r`     : Planet-Sun distance *| in AU*
**/
pub fn apprnt_mag_muller<'a> (

    planet : &Planet,
    i      : f64,
    delta  : f64,
    r      : f64

) -> Result<f64, &'a str> {

    let x = 5.0 * (r*delta).log10();

    match *planet {
        Planet::Mercury => Ok( x + 1.16 + (i - 50.0)*(0.02838 + (i - 50.0)*0.000102) ),
        Planet::Venus   => Ok( x - 4.0 + i*(0.01322 + i*i*0.0000004247) ),
        Planet::Earth   => {
            return Err("Planet::Earth was passed to the function planet::apprnt_mag_muller()");
        },
        Planet::Mars    => Ok(x - 1.3 + i*0.01486),
        Planet::Jupiter => Ok(x - 8.93),
        Planet::Saturn  => {
            return Err("Planet::Saturn was passed to the function planet::apprnt_mag_muller(). Use the function planet::saturn::apprnt_mag_muller() instead.");
        },
        Planet::Uranus  => Ok(x - 6.85),
        Planet::Neptune => Ok(x - 7.05),
    }

}

/**
Computes a planet's apparent magnitude using the Astronomical
Almanac's method adopted in 1984

# Returns

* `app_mag`: Apparent magnitude of the planet

# Arguments

* `planet`: Any variant of [Planet](./enum.Planet.html)
* `i`     : Phase angle of the planet *| in radians*
* `delta` : Planet-Earth distance *| in AU*
* `r`     : Planet-Sun distance *| in AU*
**/
pub fn apprnt_mag_84<'a> (

    planet : &Planet,
    i      : f64,
    delta  : f64,
    r      : f64

) -> Result<f64, &'a str> {

    let x = 5.0 * (r*delta).log10();

    match *planet {
        Planet::Mercury => Ok( x - 0.42 + i*(0.0380 - i*(0.000273 - i*0.00000200)) ),
        Planet::Venus   => Ok( x - 4.40 + i*(0.0009 + i*(0.000239 - i*0.00000065)) ),
        Planet::Earth   => {
            return Err("Planet::Earth was passed to the function planet::apprnt_mag_84()");
        },
        Planet::Mars    => Ok( x - 1.52 + i*0.016 ),
        Planet::Jupiter => Ok( x - 9.4 + i*0.005 ),
        Planet::Saturn  => {
            return Err("Planet::Saturn was passed to the function planet::apprnt_mag_84(). Use the function planet::saturn::apprnt_mag_84() instead.");
        },
        Planet::Uranus  => Ok( x - 7.19 ),
        Planet::Neptune => Ok( x - 6.87 ),
    }

}

/*

    NOTE: All angles passed as arguments, and those returned,
          are assumed to be radians, even if the comments
          describe them with degrees.

*/

/*

    Returns three numbers: geocentric latitude and longitude of a
    celestial body with respect to Earth, and a time correction
    (in days).

    The time correction (t) is given in order to take into account
    the finite time it takes for light to travel to Earth. And so,
    if a planet was at a given position P at time T, it will seen
    from Earth at P at time T + t.
    -----------------------------------------------------------------
     body_helio_long: The heliocentric longitude of the celestial
                      body.
      body_helio_lat: The heliocentric latitude of the celestial body.
      body_helio_rad: The heliocentric radius vector of the celestial
                      body.
    earth_helio_long: The heliocentric longitude of Earth.
     earth_helio_lat: The heliocentric latitude of Earth.
     earth_helio_rad: The heliocentric radius vector of Earth.

*/

pub fn geocentric_position(body_helio_long: f64, body_helio_lat: f64, body_helio_rad: f64,
             earth_helio_long: f64, earth_helio_lat: f64, earth_helio_rad: f64) -> (f64, f64, f64) {

    let x = body_helio_long.cos() * body_helio_lat.cos() * body_helio_rad -
            earth_helio_long.cos() * earth_helio_lat.cos() * earth_helio_rad;
    let y = body_helio_long.sin() * body_helio_lat.cos() * body_helio_rad -
            earth_helio_long.sin() * earth_helio_lat.cos() * earth_helio_rad;
    let z = body_helio_lat.sin() * body_helio_rad - earth_helio_lat.sin() * earth_helio_rad;

    (y.atan2(x),
     z.atan2((x * x + y * y).sqrt()),
     0.0057755183 * (x * x + y * y + z * z).sqrt())
}

pub fn semimaj_ax(perih: f64, ecc: f64) -> f64 {
    perih / (1.0 - ecc)
}

pub fn mean_motion(semimaj_ax: f64) -> f64 {
    0.9856076686 / (semimaj_ax.powf(3.0 / 2.0))
}

pub fn ABCabc_terms(obl_eclp: f64, long_asc_node: f64, inc: f64, perih_arg: f64) -> (f64, f64, f64, f64, f64, f64) {
    let sin_obl_eclp = obl_eclp.sin();
    let cos_obl_eclp = obl_eclp.cos();
    let cos_long_asc_node = long_asc_node.cos();
    let sin_long_asc_node = long_asc_node.sin();
    let cos_inc = inc.cos();
    let sin_inc = inc.sin();

    let f = cos_long_asc_node;
    let g = sin_long_asc_node*cos_obl_eclp;
    let h = sin_long_asc_node*sin_obl_eclp;
    let p = -1.0*sin_long_asc_node*sin_inc;
    let q = cos_long_asc_node*cos_inc*cos_obl_eclp - sin_inc*sin_obl_eclp;
    let r = cos_long_asc_node*cos_inc*sin_obl_eclp + sin_inc*cos_obl_eclp;

    let A = f.atan2(p);
    let B = g.atan2(q);
    let C = h.atan2(r);
    let a = (f*f + p*p).sqrt();
    let b = (g*g + q*q).sqrt();
    let c = (h*h + r*r).sqrt();

    (A, B, C, a, b, c)
}

pub fn geocentric_position2(obl_eclp: f64, long_asc_node: f64, inc: f64, perih_arg: f64) -> (f64, f64) {

    let (A, B, C, a, b, c) = ABCabc_terms(obl_eclp, long_asc_node, inc, perih_arg);

    let r = 0.0;
    let v = 0.0;

    let x = r * a * (A + perih_arg + v);
    let y = r * b * (B + perih_arg + v);
    let z = r * c * (C + perih_arg + v);

    let X = 0.0;
    let Y = 0.0;
    let Z = 0.0;

    let xi = X + x;
    let nu = Y + y;
    let et = Z + z;

    let mut asc = nu.atan2(xi);
    let dec = et.atan2((xi*xi + nu*nu).sqrt());

    if asc < 0.0 {
        asc += 360f64.to_radians();
    }

    (asc, dec)

}

pub fn velocity(semimaj_ax: f64, sun_dist: f64) -> f64 {
    42.1219 * (1.0/sun_dist - 1.0/(2.0 * semimaj_ax)).sqrt()
}

pub fn velocity_perih(semimaj_ax: f64, ecc: f64) -> f64 {
    29.7847/(semimaj_ax).sqrt() * ((1.0 + ecc)/(1.0 - ecc)).sqrt()
}

pub fn velocity_aph(semimaj_ax: f64, ecc: f64) -> f64 {
    29.7847/(semimaj_ax).sqrt() * ((1.0 - ecc)/(1.0 + ecc)).sqrt()
}

pub fn ramanuj_ellipse_perim(a: f64, b: f64, e: f64) -> f64 {
    3.14 * (3.0*(a + b) - ((a + 3.0*b)*(3.0*a + b)).sqrt())
}

pub fn eclipse_perim_1(a: f64, b: f64, e: f64) -> f64 {
    let A = (a + b)/2.0;
    let G = (a*b).sqrt();
    let H = (2.0 * a * b)/(a + b);
    3.14 * (21.0*A - 2.0*G - 3.0*H) / 8.0
}

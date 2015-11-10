/*

    NOTE: All angles passed as arguments, and those returned,
          are assumed to be radians, even if the comments
          describe them with degrees.

*/

/*

    Returns three numbers: geometric latitude and longitude of a
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

pub fn geometric_position(body_helio_long: f64, body_helio_lat: f64, body_helio_rad: f64,
             earth_helio_long: f64, earth_helio_lat: f64, earth_helio_rad: f64) -> (f64, f64, f64) {

    let x = body_helio_long.cos() * body_helio_lat.cos() * body_helio_rad -
            earth_helio_long.cos() * earth_helio_lat.cos() * earth_helio_rad;
    let y = body_helio_long.sin() * body_helio_lat.cos() * body_helio_rad -
            earth_helio_long.sin() * earth_helio_lat.cos() * earth_helio_rad;
    let z = body_helio_lat.sin() * body_helio_rad - earth_helio_lat.sin() * earth_helio_rad;

    println!("{}",x);
    println!("{}",y);
    println!("{}",z);
    println!("{}",0.0057755183 * (x * x + y * y + z * z).sqrt());

    (y.atan2(x),
     z.atan2((x * x + y * y).sqrt()),
     0.0057755183 * (x * x + y * y + z * z).sqrt())
}

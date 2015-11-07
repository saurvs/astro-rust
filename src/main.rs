mod time;
mod angle;
mod earth;

fn main() {

    /* Prints the Julian Day when the Apollo 11 Lunar Lander landed on the Moon */
    let y = 1969;
    let m = 7;
    let d = time::decimal_day(20, 20, 18, 4.0);
    println!("{}", time::julian_day::julian_day(y, m, d));

    /* Prints the geodesic distance between the Observatoire de Paris and the US
       Naval Observatory at Washington DC */
    let p1 = earth::surf_point{lat: angle::pure_degrees(48.0, 50.0, 11.0).to_radians(),
                               long: angle::pure_degrees(-2.0, 20.0, 14.0).to_radians(),
                              };
    let p2 = earth::surf_point{lat: angle::pure_degrees(38.0, 55.0, 17.0).to_radians(),
                               long: angle::pure_degrees(77.0, 3.0, 56.0).to_radians()
                              };
    println!("{}m", earth::dist(p1, p2));

}

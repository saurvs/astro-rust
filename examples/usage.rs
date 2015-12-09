#[macro_use]
extern crate astro_rust;

use astro_rust::*;

fn main() {

    /* Prints the Julian Day when the Apollo 11 Lunar Lander landed on the Moon */
    let moon_landing_day = time::usual_day {
        d: 20,
        h: 20,
        m: 18,
        s: 4.0
    };
    let moon_landing = time::date {
        y: 1969,
        m: 7,
        d: time::decimal_day(moon_landing_day),
        t: time::calendar_type::gregorian
    };
    println!("{}", time::julian_day(moon_landing));

    /* Prints the geodesic distance between the Observatoire de Paris and the US
       Naval Observatory at Washington DC */
    let p1 = coordinates::surf_point{lat: angle::pure_degrees(48.0, 50.0, 11.0).to_radians(),
                               long: angle::pure_degrees(-2.0, 20.0, 14.0).to_radians(),
                              };
    let p2 = coordinates::surf_point{lat: angle::pure_degrees(38.0, 55.0, 17.0).to_radians(),
                               long: angle::pure_degrees(77.0, 3.0, 56.0).to_radians()
                              };
    println!("{}m", earth::geodesic_dist(p1, p2));

    /* Prints the ecliptic coordinates of the star Pollux given its equatorial coordinates */
    println!("{}", coordinates::eclp_long_from_equa(116.328942_f64.to_radians(), 28.026183_f64.to_radians(), coordinates::oblq_2000()).to_degrees());
    println!("{}", coordinates::eclp_lat_from_equa(116.328942_f64.to_radians(), 28.026183_f64.to_radians(), coordinates::oblq_2000()).to_degrees());

}

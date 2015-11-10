mod time;
mod coordinates;
mod angular_sep;
mod angle;
mod earth;
mod julian_day;
mod sidereal;
mod refraction;
mod binary_star;
mod stellar_mag;
mod util;
mod pluto;

fn main() {

    /* Prints the Julian Day when the Apollo 11 Lunar Lander landed on the Moon */
    let moon_landing_day = time::usual_day {
        d: 20,
        h: 20,
        m: 18,
        s: 4.0
    };
    let moon_landing = time::greg_date {
        y: 1969,
        m: 7,
        d: time::decimal_day(moon_landing_day)
    };
    println!("{}", julian_day::julian_day(moon_landing));

    /* Prints the geodesic distance between the Observatoire de Paris and the US
       Naval Observatory at Washington DC */
    let p1 = coordinates::surf_point{lat: angle::pure_degrees(48.0, 50.0, 11.0).to_radians(),
                               long: angle::pure_degrees(-2.0, 20.0, 14.0).to_radians(),
                              };
    let p2 = coordinates::surf_point{lat: angle::pure_degrees(38.0, 55.0, 17.0).to_radians(),
                               long: angle::pure_degrees(77.0, 3.0, 56.0).to_radians()
                              };
    println!("{}m", earth::dist(p1, p2));

    /* Prints the ecliptical coordinates of the star Pollux given its equatorial coordinates */
    println!("{}", coordinates::ecl_long(116.328942_f64.to_radians(), 28.026183_f64.to_radians(), coordinates::oblq_ecl_2000()).to_degrees());
    println!("{}", coordinates::ecl_lat(116.328942_f64.to_radians(), 28.026183_f64.to_radians(), coordinates::oblq_ecl_2000()).to_degrees());

}

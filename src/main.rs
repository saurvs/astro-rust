mod astro_time;
mod earth;

fn main() {
    let y = 1969;
    let m = 7;
    let d = astro_time::decimal_day(20, 20, 18, 4.0);
    println!("{}", astro_time::julian_day::julian_day(y, m, d)); // day when the Apollo 11 Lunar Lander landed on the moon
}

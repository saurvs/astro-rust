use angle;
use time;

pub fn obliquity(jed: f64) -> (f64) {
    let u = time::julian_centuries(jed) / 100.0;
    
    (angle::pure_degrees(23.0, 26.0, 21.448) -
    u * (angle::pure_degrees(0.0, 0.0, 4680.93) +
        u * (angle::pure_degrees(0.0, 0.0, 1.55) +
            u * (angle::pure_degrees(0.0, 0.0, 1999.25) -
                u * (angle::pure_degrees(0.0, 0.0, 51.38) -
                    u * (angle::pure_degrees(0.0, 0.0, 249.67) +
                        u * (angle::pure_degrees(0.0, 0.0, 39.05) +
                            u * (angle::pure_degrees(0.0, 0.0, 7.12) -
                                u * (angle::pure_degrees(0.0, 0.0, 27.87) +
                                    u * (angle::pure_degrees(0.0, 0.0, 5.79) + u * angle::pure_degrees(23.0, 26.0, 2.45))))))))))).to_radians()
}

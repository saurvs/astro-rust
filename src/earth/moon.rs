use angle;
use coordinates;

pub fn illuminated_fraction_equa(sun_equa_point: coordinates::equa_point,
                                 moon_equa_point: coordinates::equa_point,
                                 earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    illuminated_fraction(angle::angular_separation(sun_equa_point, moon_equa_point).acos(),
                         earth_moon_dist, earth_sun_dist)
}

pub fn illuminated_fraction_eclip(moon_long: f64, moon_lat: f64, sun_long: f64,
                                  earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    illuminated_fraction((moon_lat.cos()*(moon_long - sun_long).cos()).acos(),
                         earth_moon_dist, earth_sun_dist)
}

fn illuminated_fraction(moon_geocen_elong: f64, earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    let i = (earth_sun_dist * moon_geocen_elong.sin()).atan2(earth_moon_dist - earth_sun_dist*moon_geocen_elong.cos());
    (1.0 + i.cos()) / 2.0
}

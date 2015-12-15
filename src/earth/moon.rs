use angle;
use coordinates;

pub fn bright_limb_pos_angle(sun_equa_point: coordinates::equator_point,
                             moon_equa_point: coordinates::equator_point) -> f64 {
    let a = sun_equa_point.dec.cos();
    let n = a * (sun_equa_point.asc - moon_equa_point.asc).cos();
    let d = sun_equa_point.dec.sin()*moon_equa_point.dec.cos() -
            a*moon_equa_point.dec.sin()*(sun_equa_point.asc - moon_equa_point.asc).cos();
    n.atan2(d)
}

/**
Computes the **illuminated fraction** of the moon from **equatorial** coordinates

# Arguments

* ```sun_equa_point```: Equatorial coordinate of the Sun (in radians)
* ```moon_equa_point```: Equatorial coordinate of the Moon (in radians)
* ```earth_moon_dist```: Distance between the Earth and it's Moon
                         (in any unit, but same as that of ```earth_sun_dist```)
* ```earth_sun_dist```: Distance between the Earth and the Sun
                        (in any unit, but same as that of ```earth_moon_dist```)
**/
pub fn illuminated_fraction_equa(sun_equa_point: coordinates::equator_point,
                                 moon_equa_point: coordinates::equator_point,
                                 earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    illuminated_fraction(angle::angular_sep(sun_equa_point, moon_equa_point).acos(),
                         earth_moon_dist, earth_sun_dist)
}

/**
Computes the **illuminated fraction** of the moon from **eclipctical** coordinates

# Arguments

* ```moon_long```: Eclipctical longitude of the Moon (in radians)
* ```moon_lat```: Eclipctical latitude of the Moon (in radians)
* ```sun_long```: Eclipctical longitude of the Sun (in radians)
* ```earth_moon_dist```: Distance between the Earth and it's Moon
                         (in any unit, but same as that of ```earth_sun_dist```)
* ```earth_sun_dist```: Distance between the Earth and the Sun
                        (in any unit, but same as that of ```earth_moon_dist```)
**/
pub fn illuminated_fraction_eclip(moon_long: f64, moon_lat: f64, sun_long: f64,
                                  earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    illuminated_fraction((moon_lat.cos()*(moon_long - sun_long).cos()).acos(),
                         earth_moon_dist, earth_sun_dist)
}

fn illuminated_fraction(moon_geocen_elong: f64, earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    let i = (earth_sun_dist * moon_geocen_elong.sin()).atan2(earth_moon_dist - earth_sun_dist*moon_geocen_elong.cos());
    (1.0 + i.cos()) / 2.0
}

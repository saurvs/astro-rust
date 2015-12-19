use angle;
use coordinates;
use time;

pub fn PositionAngleOfBrightLimb(sun_equa_point: coordinates::EquatorialPoint,
                             moon_equa_point: coordinates::EquatorialPoint) -> f64 {
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
pub fn IlluminatedFractionFromEquatorialCoords(sun_equa_point: coordinates::EquatorialPoint,
                                 moon_equa_point: coordinates::EquatorialPoint,
                                 earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    illuminated_fraction(angle::AngularSep(sun_equa_point, moon_equa_point).acos(),
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
pub fn IlluminatedFractionFromEclipticalCoords(moon_long: f64, moon_lat: f64, sun_long: f64,
                                  earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    illuminated_fraction((moon_lat.cos()*(moon_long - sun_long).cos()).acos(),
                         earth_moon_dist, earth_sun_dist)
}

fn illuminated_fraction(moon_geocen_elong: f64, earth_moon_dist: f64, earth_sun_dist: f64) -> f64 {
    let i = (earth_sun_dist * moon_geocen_elong.sin()).atan2(earth_moon_dist - earth_sun_dist*moon_geocen_elong.cos());
    (1.0 + i.cos()) / 2.0
}

/**
Returns the **times** of **passage** of the moon through the **ascending**
and **descending nodes**, close to a given date

# Return variables

```(time_of_ascending_node, time_of_descending_node)```

* ```time_of_ascending_node```: Time of passage through the ascending node
* ```time_of_descending_node```: Time of passage through the descending node

# Arguments

```date```: A ```Date``` struct
**/
pub fn TimesOfPassageThroughNodes(date: time::Date) -> (f64, f64) {
    let k = (time::DecimalYear(date) - 2000.05) * 13.4223;
    let T = k / 1342.23;
    let k1 = (k as i32) as f64;
    let k2 = (k1 as f64) + 0.5;

    (times_of_passage_through_nodes(k1, T), times_of_passage_through_nodes(k2, T))
}

fn times_of_passage_through_nodes(k: f64, T: f64) -> f64 {
    let D = (183.638 + 331.73735682*k + T*T*(0.0014852 + T*(0.00000209 - T*0.00000001))).to_radians();
    let D_times_2 = 2.0 * D;
    let M = (17.4006 + 26.8203725*k + T*T*(0.0001186 + T*0.00000006)).to_radians();
    let M1 = (38.3776 + 355.52747313*k + T*T*(0.0123499 + T*(0.000014627 - T*0.000000069))).to_radians();
    let sigma = (123.9767 - 1.44098956*k + T*T*(0.0020608 + T*(0.00000214 - T*0.000000016))).to_radians();
    let P = sigma + (272.75 - T*2.3).to_radians();
    let V = (299.75 + T*(132.85 - T*0.009173)).to_radians();

    (2451565.1619 + 27.212220817*k +
     T*T*(0.0002762 + T*(0.000000021 - T*0.000000000088)) -
     0.4721 * M1.sin() -
     0.1649 * D_times_2.sin() -
     0.0868 * (D_times_2 - M1).sin() +
     0.0084 * (D_times_2 + M1).sin() -
     0.0083 * (D_times_2 - M).sin() -
     0.0039 * (D_times_2 - M - M1).sin() +
     0.0034 * (2.0*M1).sin() -
     0.0031 * (2.0 * (D - M1)).sin() +
     0.003 * (D_times_2 + M).sin() +
     0.0028 * (M - M1).sin() +
     0.0026 * M.sin() +
     0.0025 * (2.0 * D_times_2).sin() +
     0.0024 * D.sin() +
     0.0022 * (M + M1).sin() +
     0.0017 * sigma.sin() +
     0.0014 * (2.0*D_times_2 - M1).sin() +
     0.0005 * (D_times_2 + M - M1).sin() +
     0.0004 * (D_times_2 - M + M1).sin() +
     0.0003 * (-1.0*(2.0 * (D - M)).sin() +
               (2.0*D_times_2 - M).sin() + V.sin() + P.sin()))
}

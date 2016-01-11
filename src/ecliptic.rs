use angle;
use time;
use std;

/**
Returns **mean obliquity** of the ecliptic

# Returns

* ```mean_obliquity```: The obliquity of the ecliptic *(radians)*;
  			            the angle between the Earth's equator and
		                the ecliptic.

The **mean** obliquity isn't corrected for nutation.
To obtain the **true** obliquity, use [```NutationCorrections()```]
(../fn.NutationCorrections.html) *to get the nutation correction for obliquity*,
and add it to the **mean** obliquity.

# Arguments

```JED```: Julian Ephemeris day
**/
pub fn MnOblq(JED: f64) -> (f64) {
    let u = time::JulCent(JED) / 100.0;

    (      angle::DegFrmDMS(23, 26, 21.448)
     - u * (angle::DegFrmDMS(0, 0, 4680.93)
     + u * (angle::DegFrmDMS(0, 0, 1.55)
     + u * (angle::DegFrmDMS(0, 0, 1999.25)
     - u * (angle::DegFrmDMS(0, 0, 51.38)
     - u * (angle::DegFrmDMS(0, 0, 249.67)
     + u * (angle::DegFrmDMS(0, 0, 39.05)
     + u * (angle::DegFrmDMS(0, 0, 7.12)
     - u * (angle::DegFrmDMS(0, 0, 27.87)
     + u * (angle::DegFrmDMS(0, 0, 5.79)
     + u * angle::DegFrmDMS(0, 0, 2.45)
    )))))))))).to_radians()
}

/// Returns the **obliquity** *(radians)* of the ecliptic
/// for the epoch **J2000.0**
pub fn Oblq_J2000() -> f64 {
    23.4392911_f64.to_radians()
}

/// Returns the **obliquity** *(radians)* of the ecliptic
/// for the epoch **J1950.0**
pub fn Oblq_J1950() -> f64 {
    23.4457889_f64.to_radians()
}

/**
Returns the **longitudes** of the **ecliptic points** on the **horizon**

# Returns

```(long_point_1, long_point_2)```

* ```long_point_1```: Longitude of ecliptic point 1 *(radians)*
* ```long_point_2```: Longitude of ecliptic point 2 *(radians)*

# Arguments

* ```oblq_eclip```: Obliquity of the ecliptic *(radians)*
* ```observer_lat```: Observer's geographical latitude *(radians)*
* ```loc_sidreal```: Local sidereal time *(radians)*
**/
pub fn LongOfEclipticPointsOnHz(oblq_eclip: f64, observer_lat: f64, loc_sidreal: f64) -> (f64, f64) {
    let p = (-loc_sidreal.cos())
            .atan2(   oblq_eclip.sin() * observer_lat.tan()
                    + oblq_eclip.cos() * loc_sidreal.sin()
                  );

    (p, p + std::f64::consts::PI)
}

/**
Returns the **angle** between the **ecliptic** and the **horizon**

# Returns

* ```angle```: The angle between the ecliptic and the horizon *(radians)*

# Arguments

* ```oblq_eclip```: Obliquity of the ecliptic *(radians)*
* ```observer_lat```: Observer's geographical latitude *(radians)*
* ```loc_sidreal```: Local sidereal time *(radians)*
**/
pub fn AnglEclipticAndHz(oblq_eclip: f64, observer_lat: f64, loc_sidreal: f64) -> f64 {
    (   oblq_eclip.cos() * observer_lat.sin()
      - oblq_eclip.sin() * observer_lat.cos() * loc_sidreal.sin()
    ).acos()
}

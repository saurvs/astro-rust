//! The Earth's ecliptic

use angle;
use time;
use std;

/**
Returns the **mean obliquity** of the ecliptic using
J. Laskar's formula

# Returns

* ```mn_oblq```: Mean obliquity of the ecliptic *| in radians*

The accuracy of ```mn_oblq``` is estimated to be 0.01 arcseconds
for a 1000 years before and after 2000 AD, and a few arcseconds
for a 10000 years before and afters 2000 AD.

# Arguments

* ```JD```: Julian (Ephemeris) day
**/
pub fn MnOblq(JD: f64) -> (f64) {
    let u = time::JulCent(JD) / 100.0;

    (       angle::DegFrmDMS(23, 26, 21.448)
     - u * (angle::DegFrmDMS(0, 0, 4680.93)
     + u * (angle::DegFrmDMS(0, 0, 1.55)
     + u * (angle::DegFrmDMS(0, 0, 1999.25)
     - u * (angle::DegFrmDMS(0, 0, 51.38)
     - u * (angle::DegFrmDMS(0, 0, 249.67)
     + u * (angle::DegFrmDMS(0, 0, 39.05)
     + u * (angle::DegFrmDMS(0, 0, 7.12)
     - u * (angle::DegFrmDMS(0, 0, 27.87)
     + u * (angle::DegFrmDMS(0, 0, 5.79)
     + u *  angle::DegFrmDMS(0, 0, 2.45)
    )))))))))).to_radians()
}

/**
Returns the **mean obliquity** of the ecliptic using the IAU formula

# Returns

* ```mn_oblq```: Mean obliquity of the ecliptic *| in radians*

The error in ```mn_oblq``` reaches 1 arcsecond over a period of
2000 years, and about 10 arcseconds over a period of 4000 years.

# Arguments

* ```JD```: Julian (Ephemeris) day
**/
pub fn MnOblq_IAU(JD: f64) -> (f64) {
    let u = time::JulCent(JD) / 100.0;

    (       angle::DegFrmDMS(23, 26, 21.448)
     - u * (angle::DegFrmDMS(0, 0, 46.815)
     + u * (angle::DegFrmDMS(0, 0, 0.00059)
     - u *  angle::DegFrmDMS(0, 0, 0.001813)
    ))).to_radians()
}

/**
Returns the **longitudes** of the **ecliptic points** on the **horizon**

# Returns

```(long_point_1, long_point_2)```

* ```long_point_1```: Longitude of the ecliptic point 1 *| in radians*
* ```long_point_2```: Longitude of the ecliptic point 2 *| in radians*

# Arguments

* ```oblq_eclip```  : Obliquity of the ecliptic *| in radians*
* ```observer_lat```: The observer's geographical latitude *| in radians*
* ```loc_sidreal``` : Local sidereal time *| in radians*
**/
pub fn LongOfEclipPointsOnHz(oblq_eclip: f64, observer_lat: f64, loc_sidreal: f64) -> (f64, f64) {
    let p = (-loc_sidreal.cos())
            .atan2(   oblq_eclip.sin()*observer_lat.tan()
                    + oblq_eclip.cos()*loc_sidreal.sin()
                  );

    (p, p + std::f64::consts::PI)
}

/**
Returns the **angle** between the **ecliptic** and the **horizon**

# Returns

* ```angle```: Angle between the ecliptic and
               the horizon *| in radians*

# Arguments

* ```oblq_eclip```: Obliquity of the ecliptic *| in radians*
* ```observer_lat```: The bserver's geographical latitude *| in radians*
* ```loc_sidreal```: Local sidereal time *| in radians*
**/
pub fn AnglBetwnEclipAndHz(oblq_eclip: f64, observer_lat: f64, loc_sidreal: f64) -> f64 {
    (   oblq_eclip.cos() * observer_lat.sin()
      - oblq_eclip.sin() * observer_lat.cos() * loc_sidreal.sin()
    ).acos()
}

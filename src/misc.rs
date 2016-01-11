/**
Returns the **parallactic angle** of a celestial body

# Returns

* ```parallactic_angle```: The parallactic angle *(radians)*

# Arguments

* ```observer_lat```: Observer's geographical latitude *(radians)*
* ```hour_angle```: Local hour angle *(radians)*
* ```declin```: Declination of the celestial body *(radians)*
**/
pub fn ParllcAngl(observer_lat: f64, hour_angle: f64, declin: f64) -> f64 {
    hour_angle.sin()
    .atan2(   observer_lat.tan()   * declin.cos()
            - declin.sin()         * hour_angle.cos()
          )

}

/**
Returns the **parallactic angle** of a celestial body on the horizon

# Returns

* ```parallactic_angle```: The parallactic angle on the horizon *(radians)*

# Arguments

* ```observer_lat```: Observer's geographical latitude *(radians)*
* ```declin```: Declination of the celestial body *(radians)*
**/
pub fn ParllcAnglHz(observer_lat: f64, declin: f64) -> f64 {
    (observer_lat.sin() / declin.cos()).acos()
}
/*
pub fn what(mut sidereal_greenwhich_time: f64,
            mut right_ascen_1: f64, mut declin_1: f64,
            mut right_ascen_2: f64, mut declin_2: f64,
            mut right_ascen_3: f64, mut declin_3: f64) -> f64 {

    let H_0 = ()
}*/

/**
Returns the parallactic angle of a celestial body

# Returns

* ```parallactic_angle```: The parallactic angle *| in radians*

# Arguments

* ```observer_lat```: Observer's geographical latitude *| in radians*
* ```hour_angle```: Local hour angle *| in radians*
* ```declin```: Declination of the celestial body *| in radians*
**/
pub fn ParllcAngl(observer_lat: f64, hour_angle: f64, declin: f64) -> f64 {
    hour_angle.sin()
    .atan2(   observer_lat.tan()   * declin.cos()
            - declin.sin()         * hour_angle.cos()
          )
}

/**
Returns the parallactic angle of a celestial body on the horizon

# Returns

* ```parallactic_angle```: The parallactic angle on the horizon *| in radians*

# Arguments

* ```observer_lat```: Observer's geographical latitude *| in radians*
* ```declin```: Declination of the celestial body *| in radians*
**/
pub fn ParllcAnglHz(observer_lat: f64, declin: f64) -> f64 {
    (observer_lat.sin() / declin.cos()).acos()
}

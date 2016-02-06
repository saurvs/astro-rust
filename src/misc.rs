//! Miscellaneous stuff

/**
Returns the **parallactic angle** of a celestial body

# Returns

* `parallac_angle`: The parallactic angle of the celestial
                        body *| in radians*

# Arguments

* `observer_lat`: The observer's geographical latitude *| in radians*
* `hour_angle`: Local hour angle *| in radians*
* `dec`: Declination of the celestial body *| in radians*
**/
pub fn parllc_angl(observer_lat: f64, hour_angle: f64, dec: f64) -> f64 {
    hour_angle.sin()
    .atan2(   observer_lat.tan() * dec.cos()
            - dec.sin()          * hour_angle.cos()
          )
}

/**
Returns the **parallactic angle** of a celestial body on the
**horizon**

# Returns

* `parallac_angle`: The parallactic angle of the celestial body
                        on the horizon *| in radians*

# Arguments

* `observer_lat`: The observer's geographical latitude *| in radians*
* `dec`: Declination of the celestial body *| in radians*
**/
pub fn parllc_angl_on_hz(observer_lat: f64, dec: f64) -> f64 {
    (observer_lat.sin() / dec.cos()).acos()
}

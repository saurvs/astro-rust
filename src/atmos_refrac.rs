use angle;

/**
Returns the **refraction** term for true altitude, using
**apparent altitude**

This method is valid only when the altitude of the body
is more than 15 degrees.

# Returns

* ```refraction```: The refraction term *(radians)*, that needs to be
                    subtracted from apparent altitude to get
                    true altitude.

# Arguments

* ```apparent_alt```: Apparent altitude *(radians)*
**/
pub fn RefracFromApparentAlt(apparent_alt: f64) -> f64 {
      angle::PureDegrees(0, 0, 58.294).to_radians() * (90_f64.to_radians() - apparent_alt).tan()
    - angle::PureDegrees(0, 0, 0.0668).to_radians() * (90_f64.to_radians() - apparent_alt).tan().powi(3)
}

/**
Returns the **refraction** term for apparent altitude, using
**true altitude**

This method is valid only when the altitude of the body
is more than 15 degrees.

# Returns

* ```refraction```: The refraction term *(radians)*, that needs to be
                    added to true altitude to get apparent altitude.

# Arguments

* ```true_alt```: True altitude *(radians)*
**/
pub fn RefracFromTrueAlt(true_alt: f64) -> f64 {
      angle::PureDegrees(0, 0, 58.276).to_radians() * (90_f64.to_radians() - true_alt).tan()
    - angle::PureDegrees(0, 0, 0.0824).to_radians() * (90_f64.to_radians() - true_alt).tan().powi(3)
}

/**
Returns an approximate **refraction** term for true altitude, using
**apparent altitude**

This method is valid for all values of altitude from 0 to 90 degrees;
the accuracy is upto 0.07 arcminute for all values of apparent
altitude.

# Returns

* ```refraction```: The refraction term *(radians)*, that needs to be
                    subtracted from apparent altitude to get
                    true altitude.

# Arguments

* ```apparent_alt```: Apparent altitude *(radians)*
**/
pub fn ApproxRefracFromApparentAlt(apparent_alt: f64) -> f64 {
    if (apparent_alt.to_degrees() == 90.0) { 0.0 }
    else {
        let a =   apparent_alt.to_degrees()
                + 7.31 / (apparent_alt.to_degrees() + 4.4);
        let R = 1.0 / a.to_radians().tan();

        (R / 60.0).to_radians()
    }
}

/**
Returns an approximate **refraction** term for apparent altitude, using
**true altitude**

This method is valid for all values of altitude from 0 to 90 degrees;
is consistent with ApproxRefractionFromApparentAltitude() to within
4 arcsecond.

# Returns

* ```refraction```: The refraction term *(radians)*, that needs to be
                    added to true altitude to get apparent altitude.

# Arguments

* ```true_alt```: True altitude *(radians)*
**/
pub fn ApproxRefracFromTrueAlt(true_alt: f64) -> f64 {
    if (true_alt.to_degrees() == 90.0) { 0.0 }
    else {
        let a = (   true_alt.to_degrees()
                  + 10.3 / (true_alt.to_degrees() + 5.11)
                );
        let R =   1.02 / a.to_radians().tan();

        (R / 60.0).to_radians()
    }
}

/**
Returns the **refraction** term modifier for **pressure**

# Returns

* ```refraction_term_modifier```: The value that needs to be multiplied by the
                                  refraction term to account for local pressure.

# Arguments

* ```pressure```: Local pressure *(millibars)*
**/
pub fn RefracDueToPressure(pressure: f64) -> f64 {
    pressure / 1010.0
}

/**
Returns the **refraction** term modifier for **temperature**

# Returns

* ```refraction_term_modifier```: The value that needs to be multiplied by the
                                  refraction term to account for local temperature.

# Arguments

* ```temperature```: Local temperature *(kelvins)*
**/
pub fn RefracDueToTemp(temperature: f64) -> f64 {
    283.0 / temperature
}

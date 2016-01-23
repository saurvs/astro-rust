//! Corrections for atmospheric refraction

use angle;

/**
Returns the **refraction term** for true altitudes greater than
15 degrees

# Returns

* ```refrac_term```: The refraction term *| in radians*, that needs to be
                     subtracted from the apparent altitude to get the
                     true altitude

# Arguments

* ```apprnt_alt```: Apparent altitude *| in radians*
**/
pub fn refrac_frm_apprnt_alt_15(apprnt_alt: f64) -> f64 {
      angle::DegFrmDMS(0, 0, 58.294).to_radians()*(90_f64.to_radians() - apprnt_alt).tan()
    - angle::DegFrmDMS(0, 0, 0.0668).to_radians()*(90_f64.to_radians() - apprnt_alt).tan().powi(3)
}

/**
Returns the **refraction term** for apparent altitudes greater than
15 degrees

# Returns

* ```refrac_term```: The refraction term *| in radians*, that needs to be
                     added to the true altitude to get the apparent altitude

# Arguments

* ```true_alt```: True altitude *| in radians*
**/
pub fn refrac_frm_true_alt_15(true_alt: f64) -> f64 {
      angle::DegFrmDMS(0, 0, 58.276).to_radians()*(90_f64.to_radians() - true_alt).tan()
    - angle::DegFrmDMS(0, 0, 0.0824).to_radians()*(90_f64.to_radians() - true_alt).tan().powi(3)
}

/**
Returns the **refraction term** for true altitude

# Returns

* ```refrac_term```: The refraction term *| in radians*, that needs to be
                     subtracted from the apparent altitude to get the
                      true altitude

The accuracy of ```refrac_term``` is upto 0.07 arcminutes.

# Arguments

* ```apprnt_alt```: Apparent altitude *| in radians*
**/
pub fn refrac_frm_apprnt_alt(apprnt_alt: f64) -> f64 {
    if apprnt_alt.to_degrees() == 90.0 { 0.0 }
    else {
        let a = apprnt_alt.to_degrees() + 7.31/(apprnt_alt.to_degrees() + 4.4);
        let R = 1.0 / a.to_radians().tan();

        (R / 60.0).to_radians()
    }
}

/**
Returns the **refraction term** for apparent altitude

This function is consistent with ```RefracFrmApprntAlt()``` to within
4 arcseconds.

# Returns

* ```refrac_term```: The refraction term *| in radians*, that needs to be
                     added to the true altitude to get the apparent altitude

# Arguments

* ```true_alt```: True altitude *| in radians*
**/
pub fn refrac_frm_true_alt(true_alt: f64) -> f64 {
    if true_alt.to_degrees() == 90.0 { 0.0 }
    else {
        let a = true_alt.to_degrees() + 10.3/(true_alt.to_degrees() + 5.11);
        let R =   1.02 / a.to_radians().tan();

        (R / 60.0).to_radians()
    }
}

/**
Returns the **refraction term modifier** for local **pressure**

# Returns

* ```refrac_term_modifier```: The value that needs to be multiplied by the
                              refraction term to account for local pressure

# Arguments

* ```pressure```: Local pressure *| in millibars*
**/
pub fn refrac_due_to_pressr(pressure: f64) -> f64 {
    pressure / 1010.0
}

/**
Returns the **refraction term modifier** for local **temperature**

# Returns

* ```refrac_term_modifier```: The value that needs to be multiplied by the
                              refraction term to account for local temperature

# Arguments

* ```temp```: Local temperature *| in kelvins*
**/
pub fn refrac_due_to_tempr(temp: f64) -> f64 {
    283.0 / temp
}

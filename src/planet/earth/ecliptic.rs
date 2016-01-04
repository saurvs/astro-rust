use angle;
use time;

/**
Returns **mean obliquity** of the Earth's **ecliptic** for a Julian
Ephemeris day

# Returns

```(mean_obliquity)```

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
pub fn MeanObliquity(JED: f64) -> (f64) {
    let u = time::JulianCentury(JED) / 100.0;

    (angle::PureDegrees(23.0, 26.0, 21.448) -
     u * (angle::PureDegrees(0.0, 0.0, 4680.93) +
     u * (angle::PureDegrees(0.0, 0.0, 1.55) +
     u * (angle::PureDegrees(0.0, 0.0, 1999.25) -
     u * (angle::PureDegrees(0.0, 0.0, 51.38) -
     u * (angle::PureDegrees(0.0, 0.0, 249.67) +
     u * (angle::PureDegrees(0.0, 0.0, 39.05) +
     u * (angle::PureDegrees(0.0, 0.0, 7.12) -
     u * (angle::PureDegrees(0.0, 0.0, 27.87) +
     u * (angle::PureDegrees(0.0, 0.0, 5.79) +
     u * angle::PureDegrees(0.0, 0.0, 2.45))))))))))).to_radians()
}

/// Returns **obliquity** *(radians)* of the Earth's **ecliptic**
/// for the epoch **J2000.0**
pub fn Obliquity_J2000() -> f64 {
    23.4392911_f64.to_radians()
}

/// Returns **obliquity** *(radians)* of the Earth's **ecliptic**
/// for the epoch **J1950.0**
pub fn Obliquity_J1950() -> f64 {
    23.4457889_f64.to_radians()
}

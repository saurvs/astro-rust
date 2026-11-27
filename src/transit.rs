/*
Copyright (c) 2015, 2016 Saurav Sachidanand

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

//! Time of rise, transit and set for a celestial body

use crate::angle;
use crate::coords;
use crate::interpol;

/// Represents a celestial body in transit
pub enum TransitBody {
    /// A star or a planet
    StarOrPlanet,
    /// The Sun
    Sun,
    /// The Moon
    Moon
}

/// Represents a transit type
pub enum TransitType {
    /// Rise
    Rise,
    /// Transit at zenith
    Transit,
    /// Set
    Set
}

/**
Computes the time of transit for a celestial body

# Returns

* `(hour, min, sec)`: Time of transit on the day of interest, in UTC

# Arguments

* `transit_type`  : A `TransitType`
* `transit_body`  : The `TransitBody`
* `geograph_point`: Geographic point of the observer *| in radians*

Let `JD` be the Julian (Ephemeris) day of interest,

* `eq_point1`: Equatorial point of the transit body on `JD - 1` *| in radians*
* `eq_point2`: Equatorial point of the transit body on `JD` *| in radians*
* `eq_point3`: Equatorial point of the transit body on `JD + 1` *| in radians*
* `apprnt_greenwhich_sidr`: Apparent sidereal time at Greenwhich on `JD` *| in radians*
* `delta_t`: ΔT for `JD` (Julian day)
* `moon_eq_hz_parallax`: Equatorial horizontal parallax of the Moon on `JD`
                             *| in radians*. *Pass a meaningfull value here only when*
                             `TransitBody::Moon` *is passed for* `transit_body`.

**/
pub fn time (

    transit_type           : &TransitType,
    transit_body           : &TransitBody,
    geograph_point         : &coords::GeographPoint,
    eq_point1              : &coords::EqPoint,
    eq_point2              : &coords::EqPoint,
    eq_point3              : &coords::EqPoint,
    apprnt_greenwhich_sidr : f64,
    delta_t                : f64,
    moon_eq_hz_parallax    : f64

) -> (i64, i64, f64) {

    let h0 = match transit_body {
        &TransitBody::StarOrPlanet => -0.5667_f64.to_radians(),
        &TransitBody::Sun          => -0.8333_f64.to_radians(),
        &TransitBody::Moon         =>  0.7275 * moon_eq_hz_parallax
                                     - 0.5667_f64.to_radians(),
    };

    let mut H0 = (
        (h0.sin() - geograph_point.lat.sin() * eq_point2.dec.sin()) /
        (geograph_point.lat.cos() * eq_point2.dec.cos())
    ).acos();
    H0 = angle::limit_to_two_PI(H0);

    let mut m = m(
        &transit_type, H0, eq_point2.asc, geograph_point.long,
        apprnt_greenwhich_sidr);
    let theta0 = apprnt_greenwhich_sidr + m*360.985647_f64.to_radians();

    let d = m + delta_t/86400.0;

    let asc = interpol::three_values(eq_point1.asc, eq_point2.asc, eq_point3.asc, d);

    let dec = match transit_type {
        &TransitType::Transit => 0.0,

        &TransitType::Rise    => interpol::three_values(
                                    eq_point1.dec, eq_point2.dec,
                                    eq_point3.dec, d),

        &TransitType::Set     => interpol::three_values(
                                    eq_point1.dec, eq_point2.dec,
                                    eq_point3.dec, d)
    };

    let mut H = coords::hr_angl_frm_observer_long
        (theta0, geograph_point.long, asc).to_degrees();
    H = angle::limit_to_360(H);
    if H > 180.0 { H -= 360.0; }
    H = H.to_radians();

    let h = match transit_type {
        &TransitType::Transit => 0.0,
        &TransitType::Rise    => coords::alt_frm_eq(H, dec, geograph_point.lat),
        &TransitType::Set     => coords::alt_frm_eq(H, dec, geograph_point.lat)
    };

    m += match transit_type {
        &TransitType::Transit => -H / angle::TWO_PI,
        &TransitType::Rise    => (h - h0) / (angle::TWO_PI * dec.cos() * geograph_point.lat.cos() * H.sin()),
        &TransitType::Set     => (h - h0) / (angle::TWO_PI * dec.cos() * geograph_point.lat.cos() * H.sin())
    };

    let h = 24.0 * m;
    let hour = h as i64;
    let m = (h - (hour as f64)) * 60.0;
    let minute = m as i64;
    let second = (m - (minute as f64)) * 60.0;

    (hour, minute, second)

}

#[inline]
fn m (

    transit_type : &TransitType,
    H0           : f64,
    asc          : f64,
    L            : f64,
    Theta0       : f64

) -> f64 {

    let mut m = (asc + L - Theta0)/angle::TWO_PI;
    let p = H0/angle::TWO_PI;

    m += match transit_type {
        &TransitType::Transit => 0.0,
        &TransitType::Rise    => -p,
        &TransitType::Set     => p
    };

    if      m < 0.0 { m += 1.0 }
    else if m > 1.0 { m -= 1.0 }

    m

}

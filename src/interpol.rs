//! Interpolation of intermediate values of functions

use util;

/**
Interpolates an intermediate value of a function from three of it's
given values

# Returns

* `interpol_val`: Intermediate value of the function

# Arguments

* `y1`: Value 1 of the function
* `y2`: Value 2 of the function
* `y3`: Value 3 of the function
* `n` : Interpolating factor, measured from the central value
        `y2`, positively towards `y3`
**/
pub fn three_values(y1: f64, y2: f64, y3: f64, n: f64) -> f64 {
    let a = y2 - y1;
    let b = y3 - y2;
    let c = b - a;

    y2 + n*(a + b + n*c)/2.0
}

/**
Interpolates an intermediate value of a function from five of it's
given values

# Returns

* `interpol_val`: Intermediate value of the function

# Arguments

* `y1`: Value 1 of the function
* `y2`: Value 2 of the function
* `y3`: Value 3 of the function
* `y4`: Value 4 of the function
* `y5`: Value 5 of the function
* `n` : Interpolating factor, measured from the central value
        `y3`, positively towards `y4`
**/
pub fn five_values(y1: f64, y2: f64, y3: f64, y4: f64, y5: f64, n: f64) -> f64 {
    let a = y2 - y1;
    let b = y3 - y2;
    let c = y4 - y3;
    let d = y5 - y4;

    let e = b - a;
    let f = c - b;
    let g = d - c;

    let h = f - e;
    let j = g - f;

    let k = (j - h)/12.0;

    let h_j_12 = (h + j)/6.0;

    y3 +
    Horner_eval!(
        n, 0.0,
        b + c - h_j_12,
        f - k,
        h_j_12,
        k
    )/2.0
}

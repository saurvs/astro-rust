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

//! Interpolation of intermediate values of functions

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

#[inline]
pub fn three_values(y1: f64, y2: f64, y3: f64, n: f64) -> f64
{
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
pub fn five_values(y1: f64, y2: f64, y3: f64, y4: f64, y5: f64, n: f64) -> f64
{
    let a = y2 - y1;
    let b = y3 - y2;
    let c = y4 - y3;
    let d = y5 - y4;

    let e = b - a;
    let f = c - b;
    let g = d - c;

    let h = f - e;
    let j = g - f;

    let k = (j - h) / 12.0;
    let h_j_12 = (h + j) / 6.0;

    y3 +
    Horner_eval!(
        n,
        0.0,
        b + c - h_j_12,
        f - k,
        h_j_12,
        k
    ) / 2.0
}

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

//! Some programming utilities

/// Returns a float rounded upto a certain number of decimal digits
#[inline]
pub fn round_upto_digits(float: f64, decimal_digits: u32) -> f64
{
    let mut d = 1.0;

    for _ in 1..(decimal_digits + 1) {
        d *= 10.0;
    }

    (float * d).round() / d
}

/**
Evaluates a polynomial using Horner's algorithm

# Arguments

* `$x`     : The value of the independent variable `f32 or f64`
* `$c`     : The constant term `f32 or f64`
* `$($a),*`: Sequence of coefficient terms for `$x`, in ascending
             powers of `$x`
**/
#[macro_export]
macro_rules! Horner_eval
{
    ($x:expr, $c:expr, $($a:expr),*) =>
    {
        {
            let mut y = $c;
            let mut u = 1.0;
            $(
                u *= $x;
                y += u * $a;
            )*
            y
        }
    }
}

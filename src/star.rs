use time;

/**
Returns the **combined magnitude** of two stars

* ```m1```: Magnitude of star 1
* ```m2```: Magnitude of star 2
**/
pub fn CombinedMag(m1: f64, m2: f64) -> f64 {
    m2 - 2.5*(BrightnessRatio(m1, m2) + 1.0)
}

/**
Returns the **combined magnitude** of two or more stars

* ```m```: Array of magnitudes of stars
**/
pub fn CombinedMagOfMany(m: &[f64]) -> f64 {
    let mut sum = 0.0;
    for i in m.iter() {
        sum += 10_f64.powf(-0.4 * i);
    }
    -2.5 * sum.log10()
}

/**
Returns the **brightness ratio** of two stars

* ```m1```: Magnitude of star 1
* ```m2```: Magnitude of star 2
**/
pub fn BrightnessRatio(m1: f64, m2: f64) -> f64 {
    10.0_f64.powf(0.4 * (m2-m1))
}

/**
Returns the **difference in magnitude** of two stars

* ```br```: Brightness ratio of two stars
**/
pub fn MagDifference(br: f64) -> f64 {
    2.5 * br.log10()
}

/**
Returns the **absolute magnitude** of a star from its parallax

* ```par```: Parallax of star
* ```am```: Apparent magnitude of star
**/
pub fn AbsoluteMagFromParallax(mut par: f64, am: f64) -> f64 {
    par = par.to_degrees() * 3600.0;
    am + 5.0 + 5.0*par.log10()
}

/**
Returns the **absolute magnitude** of a star from its distance from earth

* ```d```: Distance of star from earth *(parsecs)*
* ```am```: Apparent magnitude of star
**/
pub fn AbsoluteMagFromDistance(d: f64, am: f64) -> f64 {
    am + 5.0 - 5.0*d.log10()
}

/**
Returns the **abberation** corrections *(radians)* for a star's equatorial coordinates

# Returned values

```(aberration_in_right_ascension, aberration_in_declination)```

* ```aberration_in_right_ascension```: Abberation correction for right ascension *(radians)*
* ```aberration_in_declination```: Abberation correction for declination *(radians)*

# Arguments

* ```asc```: Star's right ascension *(radians)*
* ```dec```: Star's declination *(radians)*
* ```jed```: Julian Emphemeris day
**/
pub fn AberrationCorrections(asc: f64, dec: f64, jed: f64) -> (f64, f64) {
    let t = time::JulianCentury(jed);

    let l2 = 3.1761467 + 1021.3285546*t;
    let l3 = 1.7534703 + 628.3075849*t;
    let l4 = 6.2034809 + 334.0612431*t;
    let l5 = 0.5995465 + 52.9690965*t;
    let l6 = 0.8740168 + 21.3299095*t;
    let l7 = 5.4812939 + 7.4781599*t;
    let l8 = 5.3118863 + 3.8133036*t;
    let l1 = 3.8103444 + 8399.6847337*t;
    let d = 5.1984667 + 7771.3771486*t;
    let m1 = 2.3555559 + 8328.6914289*t;
    let f = 1.6279052 + 8433.4661601*t;

    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;

    let mut A: f64 = 0.0;
    let mut sinA = 0.0;
    let mut cosA = 0.0;

    // ROW 1
    A = l3;
    sinA = A.sin(); cosA = A.cos();
    x += (-1719914.0 - 2.0*t)*sinA - 25.0*cosA;
    y += (25.0 - 13.0*t)*sinA + (1578089.0 + 156.0*t)*cosA;
    z += (10.0 + 32.0*t)*sinA + (684185.0 - 358.0*t)*cosA;

    // ROW 2
    A = 2.0 * l3;
    sinA = A.sin(); cosA = A.cos();
    x += (6434.0 + 141.0*t)*sinA + (28007.0 - 107.0*t)*cosA;
    y += (25697.0 - 95.0*t)*sinA + (-5904.0 - 130.0*t)*cosA;
    z += (11141.0 - 48.0*t)*sinA + (-2559.0 - 55.0*t)*cosA;

    // ROW 3
    A = l5;
    sinA = A.sin(); cosA = A.cos();
    x += 715.0*sinA;
    y += 6.0*sinA - 657.0*cosA;
    z += -15.0*sinA - 282.0*cosA;

    // ROW 4
    A = l1;
    sinA = A.sin(); cosA = A.cos();
    x += 715.0*sinA;
    y += -656.0*cosA;
    z += -285.0*cosA;

    // ROW 5
    A = 3.0 * l3;
    sinA = A.sin(); cosA = A.cos();
    x += (486.0 - 5.0*t)*sinA + (-236.0 - 4.0*t)*cosA;
    y += (-216.0 - 4.0*t)*sinA + (-446.0 + 5.0*t)*cosA;
    z += -94.0*sinA - 193.0*cosA;

    // ROW 6
    A = l6;
    sinA = A.sin(); cosA = A.cos();
    x += 159.0*sinA;
    y += 2.0*sinA - 147.0*cosA;
    z += -6.0*sinA - 61.0*cosA;

    // ROW 7
    A = f;
    sinA = A.sin(); cosA = A.cos();
    y += 26.0*cosA;
    z += -59.0*cosA;

    // ROW 8
    A = l1 + m1;
    sinA = A.sin(); cosA = A.cos();
    x += 39.0*sinA;
    y += -36.0*cosA;
    z += -16.0*cosA;

    // ROW 9
    A = 2.0 * l5;
    sinA = A.sin(); cosA = A.cos();
    x += 33.0*sinA - 10.0*cosA;
    y += -9.0*sinA - 30.0*cosA;
    z += -5.0*sinA - 13.0*cosA;

    // ROW 10
    A = 2.0*l3 - l5;
    sinA = A.sin(); cosA = A.cos();
    x += 31.0*sinA + cosA;
    y += sinA - 28.0*cosA;
    z += -12.0*cosA;

    // ROW 11
    A = 3.0*l3 - 8.0*l4 + 3.0*l5;
    sinA = A.sin(); cosA = A.cos();
    x += 8.0*sinA - 28.0*cosA;
    y += 25.0*sinA + 8.0*cosA;
    z += 11.0*sinA + 3.0*cosA;

    // ROW 12
    A = 5.0*l3 - 8.0*l4 + 3.0*l5;
    sinA = A.sin(); cosA = A.cos();
    x += 8.0*sinA - 28.0*cosA;
    y += -25.0*sinA - 8.0*cosA;
    z += -11.0*sinA - 3.0*cosA;

    // ROW 13
    A = 2.0*l2 - l3;
    sinA = A.sin(); cosA = A.cos();
    x += 21.0*sinA;
    y += -19.0*cosA;
    z += -8.0*cosA;

    // ROW 14
    A = l2;
    sinA = A.sin(); cosA = A.cos();
    x += -19.0*sinA;
    y += 17.0*cosA;
    z += 8.0*cosA;

    // ROW 15
    A = l7;
    sinA = A.sin(); cosA = A.cos();
    x += 17.0*sinA;
    y += -16.0*cosA;
    z += -7.0*cosA;

    // ROW 16
    A = l3 - 2.0*l5;
    sinA = A.sin(); cosA = A.cos();
    x += 16.0*sinA;
    y += 15.0*cosA;
    z += sinA + 7.0*cosA;

    // ROW 17
    A = l8;
    sinA = A.sin(); cosA = A.cos();
    x += 16.0*sinA;
    y += sinA - 15.0*cosA;
    z += -3.0*sinA - 6.0*cosA;

    // ROW 18
    A = l3 + l5;
    sinA = A.sin(); cosA = A.cos();
    x += 11.0*sinA - cosA;
    y += -1.0*sinA - 10.0*cosA;
    z += -1.0*sinA - 5.0*cosA;

    // ROW 19
    A = 2.0 * (l2 - l3);
    sinA = A.sin(); cosA = A.cos();
    x += -11.0*cosA;
    y += -10.0*sinA;
    z += -4.0*sinA;

    // ROW 20
    A = l3 - l5;
    sinA = A.sin(); cosA = A.cos();
    x += -11.0*sinA - 2.0*cosA;
    y += -2.0*sinA + 9.0*cosA;
    z += -1.0*sinA + 4.0*cosA;

    // ROW 21
    A = 4.0*l3;
    sinA = A.sin(); cosA = A.cos();
    x += -7.0*sinA - 8.0*cosA;
    y += -8.0*sinA + 6.0*cosA;
    z += 3.0 * (cosA - sinA);

    // ROW 22
    A = 3.0*l3 - 2.0*l5;
    sinA = A.sin(); cosA = A.cos();
    x += -10.0*sinA;
    y += 9.0*cosA;
    z += 4.0*cosA;

    // ROW 23
    A = l2 - 2.0*l3;
    sinA = A.sin(); cosA = A.cos();
    x += -9.0*sinA;
    y += -9.0*cosA;
    z += -4.0*cosA;

    // ROW 24
    A = 2.0*l2 - 3.0*l3;
    sinA = A.sin(); cosA = A.cos();
    x += -9.0*sinA;
    y += -8.0*cosA;
    z += -4.0*cosA;

    // ROW 25
    A = 2.0*l6;
    sinA = A.sin(); cosA = A.cos();
    x += -9.0*cosA;
    y += -8.0*sinA;
    z += -3.0*sinA;

    // ROW 26
    A = 2.0*l2 - 4.0*l3;
    sinA = A.sin(); cosA = A.cos();
    x += -9.0*cosA;
    y += 8.0*sinA;
    z += 3.0*sinA;

    // ROW 27
    A = 3.0*l3 - 2.0*l4;
    sinA = A.sin(); cosA = A.cos();
    x += 8.0*sinA;
    y += -8.0*cosA;
    z += -3.0*cosA;

    // ROW 28
    A = l1 + 2.0*d - m1;
    sinA = A.sin(); cosA = A.cos();
    x += 8.0*sinA;
    y += -7.0*cosA;
    z += -3.0*cosA;

    // ROW 29
    A = 8.0*l2 - 12.0*l3;
    sinA = A.sin(); cosA = A.cos();
    x += -4.0*sinA - 7.0*cosA;
    y += -6.0*sinA + 4.0*cosA;
    z += -3.0*sinA + 2.0*cosA;

    // ROW 30
    A = 8.0*l2 - 14.0*l3;
    sinA = A.sin(); cosA = A.cos();
    x += -4.0*sinA - 7.0*cosA;
    y += 6.0*sinA - 4.0*cosA;
    z += 3.0*sinA - 2.0*cosA;

    // ROW 31
    A = 2.0*l4;
    sinA = A.sin(); cosA = A.cos();
    x += -6.0*sinA - 5.0*cosA;
    y += -4.0*sinA + 5.0*cosA;
    z += 2.0 * (cosA - sinA);

    // ROW 32
    A = 3.0*l2 - 4.0*l3;
    sinA = A.sin(); cosA = A.cos();
    x += -1.0 * (sinA + cosA);
    y += -2.0*sinA - 7.0*cosA;
    z += sinA - 4.0*cosA;

    // ROW 33
    A = 2.0 * (l3 - l5);
    sinA = A.sin(); cosA = A.cos();
    x += 4.0*sinA - 6.0*cosA;
    y += -5.0*sinA - 4.0*cosA;
    z += -2.0 * (sinA + cosA);

    // ROW 34
    A = 3.0 * (l2 - l3);
    sinA = A.sin(); cosA = A.cos();
    x += -7.0*cosA;
    y += -6.0*sinA;
    z += -3.0*sinA;

    // ROW 35
    A = 2.0 * (l3 - l4);
    sinA = A.sin(); cosA = A.cos();
    x += 5.0 * (sinA - cosA);
    y += -4.0*sinA - 5.0*cosA;
    z += -2.0 * (sinA + cosA);

    // ROW 36
    A = l1 - 2.0*d;
    sinA = A.sin(); cosA = A.cos();
    x += 5.0*sinA;
    y += -5.0*cosA;
    z += -2.0*cosA;

    let c = 17314463350.0;

    let delta_asc = (y*asc.cos() - x*asc.sin()) / (c*dec.cos());
    let delta_dec = -(((x*asc.cos() + y*asc.sin()) * dec.sin() - z*dec.cos())) / c;

    (delta_asc, delta_dec)
}

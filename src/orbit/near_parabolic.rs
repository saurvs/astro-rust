//! Near-parabolic orbits

use angle;
use consts;

/**
Returns the **true anomaly** and **radius vector** of a body in a
near-parabolic orbit

# Returns

```(true_anom, rad_vec)```

* ```true_anom```: True anomaly of the body *| in radians*
* ```rad_vec```: Radius vector of the body *| in AU*

# Arguments

* ```t```: The time of interest, in Julian (Ephemeris) day
* ```T```: Time of passage in perihelion, in Julian (Ephemeris) day
* ```ecc```: Eccentricity of the orbit
* ```q```: Perihelion distance *| in AU*
* ```accuracy```: Desired accuracy for the results. *Eg: 0.000001 radians*
**/
pub fn true_anom_and_rad_vec(t: f64, T: f64, ecc: f64, q: f64, accuracy: f64) -> (f64, f64) {
    let days_frm_perih = t - T;
    if days_frm_perih == 0.0 {
        return (0.0, q)
    }

    let k = consts::gauss_grav();
    let d1 = 1000.0;

    let q1 = k * ((1.0 + ecc)/q).sqrt() / (2.0*q);
    let mut g = (1.0 - ecc)/(1.0 + ecc);

    let q2 = q1 * days_frm_perih;
    let mut s = 2.0 / (3.0 * q2.abs());
    s = 2.0 / (2.0 * (s.atan()/2.0).tan().powf(1.0/3.0).atan()).tan();

    if days_frm_perih < 0.0 { s = -s; }
    if ecc != 1.0 {
        let mut l = 0.0;

        loop {
            let s0 = s;
            let mut z = 1.0;
            let y = s*s;
            let mut g1 = -y*s;
            let mut q3 = q2 + 2.0*g*s*y/3.0;

            loop {
                z += 1.0;
                g1 = -g1 * g * y;
				let z1 = (z - (z + 1.0)*g)/(2.0*z + 1.0);
				let f = z1 * g1;
				q3 += f;
				if z > 50.0 || f.abs() > d1 {
					panic!("No convergence at
                            orbit::near_parabolic::TruAnomAndRadVec()");
				}
				if f.abs() <= accuracy { break; }
            }

            l += 1.0;
            if l > 50.0 {
                panic!("No convergence at
                        orbit::near_parabolic::TruAnomAndRadVec()");
            }
            loop {
                let s1 = s;
                s = (2.0*s*s*s/3.0 + q3)/(s*s + 1.0);
                if (s - s1).abs() <= accuracy { break; }
            }
            if (s - s0).abs() <= accuracy { break; }
        }
    }

    let mut v = 2.0 * s.atan();
    let r = q * (1.0 + ecc)/(1.0 + ecc*v.cos());

    v = angle::LimitTo360(v.to_degrees()).to_radians();

    (v, r)
}

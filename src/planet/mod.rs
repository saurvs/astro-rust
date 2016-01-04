pub mod mercury;
pub mod venus;
pub mod earth;
pub mod mars;
pub mod jupiter;
pub mod saturn;
pub mod uranus;
pub mod neptune;

fn VSOP87Coordinate(t: f64, a: &[[f64; 3]], b: &[[f64; 3]], c: &[[f64; 3]], d: &[[f64; 3]], e: &[[f64; 3]], f: &[[f64; 3]]) -> f64 {
    let mut T0 = 0.0; for &i in a.iter() { T0 += VSOP87Term(t, &i); }
    let mut T1 = 0.0; for &i in b.iter() { T1 += VSOP87Term(t, &i); }
    let mut T2 = 0.0; for &i in c.iter() { T2 += VSOP87Term(t, &i); }
    let mut T3 = 0.0; for &i in d.iter() { T3 += VSOP87Term(t, &i); }
    let mut T4 = 0.0; for &i in e.iter() { T4 += VSOP87Term(t, &i); }
    let mut T5 = 0.0; for &i in f.iter() { T5 += VSOP87Term(t, &i); }
println!("{:?}", T0);
    T0 +
    t * (T1 +
    t * (T2 +
    t * (T3 +
    t * (T4 +
    t * T5 ))))
}

fn VSOP87Term(t: f64, array: &[f64]) -> f64 {
    array[0] * (array[1] + t*array[2]).cos()
}

pub fn GeocentricCoordsOfPlanet(L: f64, B: f64, R: f64, L0: f64, B0: f64, R0: f64) -> (f64, f64, f64) {
    let x = R*B.cos()*L.cos() - R0*B0.cos()*L0.cos();
    let y = R*B.cos()*L.sin() - R0*B0.cos()*L0.sin();
    let z = R*B.sin() - R0*B0.sin();

    (y.atan2(x),
     z/(x*x + y*y).sqrt(),
     0.0057755183 * (x*x + y*y + z*z).sqrt())
}

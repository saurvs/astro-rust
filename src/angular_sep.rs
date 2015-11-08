use coordinates;

fn angular_separation(p1: coordinates::celes_point, p2: coordinates::celes_point) -> f64 {
    (p1.dec.sin() * p2.dec.sin() +
     p1.dec.cos() * p2.dec.cos() * (p1.asc - p2.asc).cos()
    ).cos()
}

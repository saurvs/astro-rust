//! Meeus, *Astronomical Algorithms*, 2nd ed., chapter 58: planar sundials.

mod common;

use astro::sundial::{self, ShadowPoint, Sundial};
use common::{assert_close, assert_close_deg};

fn r(deg: f64) -> f64 {
    deg.to_radians()
}

/// Independent check of a shadow point with vectors in the local
/// (east, north, up) frame: the shadow of the stylus tip a·n̂ falls along
/// the Sun's direction ŝ onto the plane through the origin normal to n̂.
fn vector_shadow(d: &Sundial, h: f64, dec: f64) -> (f64, f64, f64) {
    let (phi, gd, z, a) = (d.lat, d.gnomonic_dec, d.zenith_dist, d.stylus_len);
    let s = [
        -dec.cos() * h.sin(),
        phi.cos() * dec.sin() - phi.sin() * dec.cos() * h.cos(),
        phi.sin() * dec.sin() + phi.cos() * dec.cos() * h.cos(),
    ];
    let n = [-z.sin() * gd.sin(), -z.sin() * gd.cos(), z.cos()];
    let xh = [gd.cos(), -gd.sin(), 0.0];
    let yh = [z.cos() * gd.sin(), z.cos() * gd.cos(), z.sin()];
    let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    let sn = dot(s, n);
    (-a * dot(s, xh) / sn, -a * dot(s, yh) / sn, sn)
}

fn sample_dials() -> Vec<Sundial> {
    vec![
        Sundial::new(r(40.0), r(70.0), r(50.0), 1.0),
        Sundial::vertical(r(-35.0), r(160.0), 1.0),
        Sundial::new(r(40.0), r(160.0), r(75.0), 1.0),
        Sundial::horizontal(r(51.5), 2.5),
        Sundial::vertical(r(48.0), r(-20.0), 1.0),
        Sundial::new(r(-60.0), r(250.0), r(120.0), 0.7),
        Sundial::equatorial_north(r(45.0), 1.0),
        Sundial::equatorial_south(r(45.0), 1.0),
    ]
}

// Meeus, example 58.a: φ = 40° N, D = 70°, z = 50°, a = 1.
#[test]
fn inclined_dial_ex58a() {
    let d = Sundial::new(r(40.0), r(70.0), r(50.0), 1.0);

    // Sunrise at the summer solstice: H0 = -111.33°.
    assert_close(
        sundial::sunset_hour_angle(r(40.0), r(23.44)).to_degrees(),
        111.33,
        0.005,
    );

    // Q < 0 at H = -105°, -90°, -75°, -60°; lit from H = -47° on.
    let dec = r(23.44);
    for h in [-105.0, -90.0, -75.0, -60.0, -48.0] {
        assert!(!d.is_lit(r(h), dec), "H = {h}");
        assert_eq!(d.shadow_point(r(h), dec), None);
    }
    let mut h = -47.0;
    while h <= 111.0 {
        assert!(d.is_lit(r(h), dec), "H = {h}");
        assert!(d.shadow_point(r(h), dec).is_some(), "H = {h}");
        h += 1.0;
    }
    // After sunset the Sun is down even though Q is positive.
    assert!(d.is_lit(r(112.0), dec));
    assert_eq!(d.shadow_point(r(112.0), dec), None);

    let p = d.shadow_point(r(30.0), r(23.44)).unwrap();
    assert_close(p.x, -0.0390, 0.00005);
    assert_close(p.y, -0.3615, 0.00005);
    let p = d.shadow_point(r(-15.0), r(-11.47)).unwrap();
    assert_close(p.x, -2.0007, 0.00005);
    assert_close(p.y, -1.1069, 0.00005);

    let c = d.polar_stylus().unwrap();
    assert_close(c.x0, 3.3880, 0.00005);
    assert_close(c.y0, -3.1102, 0.00005);
    assert_close(c.angle.to_degrees(), 12.2672, 0.00005);
}

// Meeus, example 58.b: vertical dial at φ = -35°, D = 160°, a = 1.
#[test]
fn vertical_dial_ex58b() {
    let d = Sundial::new(r(-35.0), r(160.0), r(90.0), 1.0);
    assert_eq!(d, Sundial::vertical(r(-35.0), r(160.0), 1.0));

    // At the equinox H0 = 90°; the dial is lit from H = -57° to sunset.
    assert_close(
        sundial::sunset_hour_angle(r(-35.0), 0.0).to_degrees(),
        90.0,
        1e-12,
    );
    assert!(!d.is_lit(r(-58.0), 0.0));
    assert!(d.is_lit(r(-56.0), 0.0));
    assert!(d.is_lit(r(89.0), 0.0));
    assert!(d.shadow_point(r(-45.0), 0.0).is_some());
    assert!(d.shadow_point(r(90.5), 0.0).is_none());

    let p = d.shadow_point(r(45.0), 0.0).unwrap();
    assert_close(p.x, -0.8439, 0.00005);
    assert_close(p.y, -0.9298, 0.00005);
    let p = d.shadow_point(0.0, r(20.15)).unwrap();
    assert_close(p.x, 0.3640, 0.00005);
    assert_close(p.y, -0.7410, 0.00005);

    let c = d.polar_stylus().unwrap();
    assert_close(c.x0, 0.3640, 0.00005);
    assert_close(c.y0, 0.7451, 0.00005);
    assert_close(c.angle.to_degrees(), 50.3315, 0.00005);
}

// Meeus, example 58.c: φ = 40° N, D = 160°, z = 75°. At the summer
// solstice the dial is lit from sunrise (H = -111°) to H = -84°, then again
// from H = +2° to sunset (H = +111°).
#[test]
fn twice_lit_dial_ex58c() {
    let d = Sundial::new(r(40.0), r(160.0), r(75.0), 1.0);
    let dec = r(23.44);
    for h in [-105.0, -90.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0] {
        assert!(d.shadow_point(r(h), dec).is_some(), "H = {h}");
    }
    for h in [-75.0, -60.0, -45.0, -30.0, -15.0, 0.0] {
        assert!(d.shadow_point(r(h), dec).is_none(), "H = {h}");
    }
    // The book gives the last and first whole degree of illumination
    // (the exact crossings are at H = -83.2° and +1.9°).
    assert!(d.is_lit(r(-84.0), dec) && !d.is_lit(r(-83.0), dec));
    assert!(!d.is_lit(r(1.0), dec) && d.is_lit(r(2.0), dec));
    assert!(d.shadow_point(r(-111.0), dec).is_some());
    assert!(d.shadow_point(r(-111.5), dec).is_none());
}

// The general formulas against the closed forms Meeus gives for the
// special cases.
#[test]
fn horizontal_dial_special_case() {
    for &phi in &[r(40.0), r(-33.9), r(60.0)] {
        let a = 1.7;
        let d = Sundial::horizontal(phi, a);
        for &dec in &[r(-23.44), r(-11.47), 0.0, r(20.15), r(23.44)] {
            for k in -7..=7 {
                let h = r(15.0 * k as f64);
                let q = phi.cos() * h.cos() + phi.sin() * dec.tan();
                let Some(p) = d.shadow_point(h, dec) else {
                    continue;
                };
                assert_close(p.x, a * h.sin() / q, 1e-12);
                assert_close(
                    p.y,
                    a * (phi.sin() * h.cos() - phi.cos() * dec.tan()) / q,
                    1e-12,
                );
            }
        }
        let c = d.polar_stylus().unwrap();
        assert_close(c.x0, 0.0, 1e-15);
        assert_close(c.y0, -a / phi.tan(), 1e-12);
        assert_close(c.length, a / phi.sin().abs(), 1e-12);
        assert_close(c.angle, phi.abs(), 1e-12);
    }
}

#[test]
fn vertical_dial_special_case() {
    for &(phi, gd) in &[(r(52.0), r(0.0)), (r(40.0), r(-35.0)), (r(-35.0), r(160.0))] {
        let a = 1.0;
        let d = Sundial::vertical(phi, gd, a);
        for &dec in &[r(-23.44), 0.0, r(11.47), r(23.44)] {
            for k in -7..=7 {
                let h = r(15.0 * k as f64);
                let q = gd.sin() * h.sin() + phi.sin() * gd.cos() * h.cos()
                    - phi.cos() * gd.cos() * dec.tan();
                let Some(p) = d.shadow_point(h, dec) else {
                    continue;
                };
                let x = a
                    * (gd.cos() * h.sin() - phi.sin() * gd.sin() * h.cos()
                        + phi.cos() * gd.sin() * dec.tan())
                    / q;
                let y = -a * (phi.cos() * h.cos() + phi.sin() * dec.tan()) / q;
                assert_close(p.x, x, 1e-12);
                assert_close(p.y, y, 1e-12);
            }
        }
        let c = d.polar_stylus().unwrap();
        assert_close(c.x0, -a * gd.tan(), 1e-12);
        assert_close(c.y0, a * phi.tan() / gd.cos(), 1e-12);
        assert_close(c.length, a / (phi.cos() * gd.cos()).abs(), 1e-12);
    }
}

#[test]
fn equatorial_dial_special_case() {
    let phi = r(45.0);
    let a = 1.0;
    let north = Sundial::equatorial_north(phi, a);
    let south = Sundial::equatorial_south(phi, a);
    for k in -5..=5 {
        let h = r(15.0 * k as f64);
        // Northern face: lit for positive declinations only.
        let dec = r(20.0);
        let p = north.shadow_point(h, dec).unwrap();
        assert_close(p.x, -a * h.sin() / dec.tan(), 1e-12);
        assert_close(p.y, -a * h.cos() / dec.tan(), 1e-12);
        assert!(south.shadow_point(h, dec).is_none());
        // The shadow tip runs on a circle about the foot of the stylus.
        assert_close(p.x.hypot(p.y), a / dec.tan(), 1e-12);

        let dec = r(-10.0);
        if let Some(p) = south.shadow_point(h, dec) {
            assert_close(p.x, -a * h.sin() / dec.tan(), 1e-12);
            assert_close(p.y, a * h.cos() / dec.tan(), 1e-12);
        }
        assert!(north.shadow_point(h, dec).is_none());
    }
    // 12h falls on the line of greatest descending slope.
    let p = north.shadow_point(0.0, r(15.0)).unwrap();
    assert_close(p.x, 0.0, 1e-15);
    assert!(p.y < 0.0);
    for d in [north, south] {
        let c = d.polar_stylus().unwrap();
        assert_close(c.x0, 0.0, 1e-12);
        assert_close(c.y0, 0.0, 1e-12);
        assert_close(c.length, a, 1e-12);
        assert_close(c.angle.to_degrees(), 90.0, 1e-6);
    }
}

// Every shadow point agrees with an independent vector computation, and
// Q > 0 exactly when the Sun is in front of the face.
#[test]
fn shadow_points_match_vector_geometry() {
    for d in sample_dials() {
        for &dec in &[r(-23.44), r(-11.47), 0.0, r(11.47), r(23.44)] {
            let mut h = -180.0;
            while h <= 180.0 {
                let (x, y, sn) = vector_shadow(&d, r(h), dec);
                // Skip the Sun grazing the plane (eg: equinox on an
                // equatorial dial), where the sign is rounding noise.
                if sn.abs() < 1e-12 {
                    h += 7.5;
                    continue;
                }
                assert_eq!(d.is_lit(r(h), dec), sn > 0.0, "{d:?} H = {h}");
                if let Some(p) = d.shadow_point(r(h), dec) {
                    let tol = 1e-9 * (1.0 + x.abs() + y.abs());
                    assert_close(p.x, x, tol);
                    assert_close(p.y, y, tol);
                }
                h += 7.5;
            }
        }
    }
}

// Hour lines are straight lines through the center of the dial, and on the
// equinox the shadow tip moves along a straight line.
#[test]
fn hour_lines_converge_on_the_center() {
    // Sine of the angle between the directions o→a and o→b.
    let cross = |o: (f64, f64), a: ShadowPoint, b: ShadowPoint| {
        let (ax, ay, bx, by) = (a.x - o.0, a.y - o.1, b.x - o.0, b.y - o.1);
        (ax * by - ay * bx) / (ax.hypot(ay) * bx.hypot(by))
    };
    // Keep to the part of the plane a real dial would cover.
    let on_dial = |p: &ShadowPoint| p.x.hypot(p.y) < 100.0;
    for d in sample_dials() {
        let Some(c) = d.polar_stylus() else { continue };
        for k in -8..=8 {
            let pts: Vec<_> = d
                .hour_line(r(15.0 * k as f64), r(23.44))
                .into_iter()
                .flatten()
                .filter(on_dial)
                .collect();
            for w in pts.windows(2) {
                assert!(cross((c.x0, c.y0), w[0], w[1]).abs() < 1e-9);
            }
        }
        let eq: Vec<_> = (-12..=12)
            .filter_map(|k| d.shadow_point(r(7.5 * k as f64), 0.0))
            .filter(on_dial)
            .collect();
        for w in eq.windows(3) {
            assert!(cross((w[0].x, w[0].y), w[1], w[2]).abs() < 1e-9);
        }
    }
}

// The tip of the polar stylus is the tip of the perpendicular stylus, and
// the polar stylus is parallel to the Earth's axis.
#[test]
fn polar_stylus_geometry() {
    for d in sample_dials() {
        let Some(c) = d.polar_stylus() else { continue };
        let a = d.stylus_len;
        assert_close(
            (c.x0 * c.x0 + c.y0 * c.y0 + a * a).sqrt(),
            c.length,
            1e-9 * c.length,
        );
        assert_close(c.angle.sin(), a / c.length, 1e-12);
        // A star at the celestial pole (δ = 90°) casts the shadow of the tip
        // on the center: the limit of the shadow as δ → 90°.
        let p = vector_shadow(&d, 0.3, r(89.999_999));
        if p.2 > 0.0 {
            assert_close(p.0, c.x0, 1e-5 * (1.0 + c.x0.abs()));
            assert_close(p.1, c.y0, 1e-5 * (1.0 + c.y0.abs()));
        }
    }
    // A polar dial (plane parallel to the Earth's axis) has no center.
    // For D = 0 that happens at z = φ.
    let polar = Sundial::new(r(40.0), 0.0, r(40.0), 1.0);
    assert!(polar.polar_stylus().is_none());
    // Vertical east-facing dial: also parallel to the axis.
    assert!(
        Sundial::vertical(r(40.0), r(270.0), 1.0)
            .polar_stylus()
            .is_none()
    );
}

#[test]
fn hour_line_declinations() {
    // The seven declinations are -ε, -20.15°, -11.47°, 0, ... for ε = 23.44°.
    let d = Sundial::horizontal(r(40.0), 1.0);
    let line = d.hour_line(r(15.0), r(23.44));
    let expected = [-23.44, -20.15, -11.47, 0.0, 11.47, 20.15, 23.44];
    for (p, dec) in line.iter().zip(expected) {
        let q = d.shadow_point(r(15.0), r(dec)).unwrap();
        let p = p.unwrap();
        assert_close(p.x, q.x, 2e-4);
        assert_close(p.y, q.y, 2e-4);
    }
    // At 3h solar time the Sun is down at all seasons at 40° N.
    assert!(d.hour_line(r(-135.0), r(23.44)).iter().all(Option::is_none));
}

#[test]
fn sunset_hour_angle_limits() {
    assert_close(
        sundial::sunset_hour_angle(r(40.0), 0.0).to_degrees(),
        90.0,
        1e-12,
    );
    assert_close(
        sundial::sunset_hour_angle(0.0, r(23.44)).to_degrees(),
        90.0,
        1e-12,
    );
    // Days are long in the summer of each hemisphere.
    assert!(sundial::sunset_hour_angle(r(40.0), r(20.0)) > r(90.0));
    assert!(sundial::sunset_hour_angle(r(-40.0), r(20.0)) < r(90.0));
    // Polar day and polar night.
    assert_close(
        sundial::sunset_hour_angle(r(80.0), r(20.0)),
        std::f64::consts::PI,
        0.0,
    );
    assert_close(sundial::sunset_hour_angle(r(80.0), r(-20.0)), 0.0, 0.0);
}

#[test]
fn hour_angle_from_solar_time() {
    assert_close_deg(
        sundial::hour_angle_frm_solar_time(12.0).to_degrees(),
        0.0,
        1e-12,
    );
    assert_close(
        sundial::hour_angle_frm_solar_time(9.0).to_degrees(),
        -45.0,
        1e-12,
    );
    assert_close(
        sundial::hour_angle_frm_solar_time(13.0).to_degrees(),
        15.0,
        1e-12,
    );
    // 4h35m a.m. (Meeus, example 58.a): H = -111.25°
    assert_close(
        sundial::hour_angle_frm_solar_time(4.0 + 35.0 / 60.0).to_degrees(),
        -111.25,
        1e-12,
    );
}

#[test]
fn bad_input_gives_none() {
    let d = Sundial::horizontal(r(40.0), 1.0);
    assert_eq!(d.shadow_point(f64::NAN, 0.0), None);
    assert_eq!(d.shadow_point(0.0, f64::NAN), None);
    assert!(!d.is_lit(f64::NAN, 0.0));
    let bad = Sundial::new(f64::NAN, 0.0, 0.0, 1.0);
    assert_eq!(bad.polar_stylus(), None);
    assert_eq!(bad.shadow_point(0.0, 0.0), None);
}

// Tests for the `angle` module: angular separation (Meeus, ch. 17),
// sexagesimal conversion (ch. 1) and angle normalisation.

mod common;

use astro::angle::{self, Dms, Hms, TWO_PI};
use common::{assert_close, dms, hms};
use std::f64::consts::PI;

fn rad(deg: f64) -> f64 {
    deg.to_radians()
}

/// The two floating point neighbours of `x` (`f64::next_up` needs Rust
/// 1.86, above the crate's MSRV).
fn neighbours(x: f64) -> (f64, f64) {
    if x == 0.0 {
        (f64::from_bits(1), -f64::from_bits(1))
    } else {
        (
            f64::from_bits(x.to_bits() + 1),
            f64::from_bits(x.to_bits() - 1),
        )
    }
}

// ---------------------------------------------------------------------------
// Angular separation
// ---------------------------------------------------------------------------

// Meeus, example 17.a: Arcturus and Spica, cos d = 0.840633, d = 32.7930°
#[test]
fn anglr_sepr_example_17_a() {
    let (a1, d1) = (rad(213.9154), rad(19.1825));
    let (a2, d2) = (rad(201.2983), rad(-11.1614));

    let d = angle::anglr_sepr(a1, d1, a2, d2);
    assert_close(d.cos(), 0.840633, 5e-7);
    assert_close(d.to_degrees(), 32.7930, 5e-5);

    // 32°48′
    let r = Dms::from_deg(d.to_degrees());
    assert_eq!((r.deg, (r.min as f64 + r.sec / 60.0).round()), (32, 48.0));

    // The same from the sexagesimal catalogue values
    let d = angle::anglr_sepr(
        rad(hms(14, 15, 39.7)),
        rad(dms(false, 19, 10, 57.0)),
        rad(hms(13, 25, 11.6)),
        rad(dms(true, 11, 9, 41.0)),
    );
    assert_close(d.to_degrees(), 32.7930, 5e-5);
}

// Meeus, example 17.b gives tan d = 0.644283 for the same pair. The book
// computes it from x, y, z rounded to six decimals (which give 0.6442835);
// the unrounded value is 0.6442838.
#[test]
fn anglr_sepr_example_17_b() {
    let d = angle::anglr_sepr(rad(213.9154), rad(19.1825), rad(201.2983), rad(-11.1614));
    assert_close(d.tan(), 0.644283, 1e-6);
}

// Meeus, ch. 17, exercise: Aldebaran–Antares, 169°58′. The J2000 positions
// are catalogue values (the book does not list them).
#[test]
fn anglr_sepr_aldebaran_antares() {
    let d = angle::anglr_sepr(
        rad(hms(4, 35, 55.24)),
        rad(dms(false, 16, 30, 33.5)),
        rad(hms(16, 29, 24.46)),
        rad(dms(true, 26, 25, 55.2)),
    );
    assert_close(d.to_degrees(), dms(false, 169, 58, 0.0), 0.5 / 60.0);
}

// Formula 17.1 returned NaN here because sin²δ + cos²δ rounds to 1 + 2⁻⁵²
// and the arccosine of that is undefined.
#[test]
fn anglr_sepr_coincident_points_is_zero_not_nan() {
    assert_eq!(angle::anglr_sepr(1.0, 0.017, 1.0, 0.017), 0.0);
    for i in -1000..=1000 {
        let lat = i as f64 * 1.5e-3;
        let d = angle::anglr_sepr(0.7, lat, 0.7, lat);
        assert_eq!(d, 0.0, "lat {lat}");
        let d = angle::anglr_sepr(0.7, lat, 0.7 + TWO_PI, lat);
        assert!(d < 1e-15, "lat {lat}: {d}");
    }
}

// Formula 17.1 returned 0 for a separation of 1e-9 rad (0.2 mas), and is
// poor below about 1″ (Meeus, ch. 17).
#[test]
fn anglr_sepr_small_separations_keep_full_precision() {
    for &e in &[1e-15, 1e-12, 1e-9, 1e-6, 1e-3] {
        // Along a meridian, the separation is exactly the difference in
        // declination. (`e` itself is rounded when added to 0.2, so compare
        // with the difference actually represented.)
        let de = (0.2 + e) - 0.2;
        let d = angle::anglr_sepr(0.3, 0.2, 0.3, 0.2 + e);
        assert!((d - de).abs() <= 4.0 * f64::EPSILON * de, "e {e}: {d}");
        // Along the equator, exactly the difference in right ascension.
        let de = (1.0 + e) - 1.0;
        let d = angle::anglr_sepr(1.0, 0.0, 1.0 + e, 0.0);
        assert!((d - de).abs() <= 4.0 * f64::EPSILON * de, "e {e}: {d}");
        // Haversine form agrees
        let de = (0.2 + e) - 0.2;
        let h = angle::anglr_sepr_haversine(0.3, 0.2, 0.3, 0.2 + e);
        assert!((h - de).abs() <= 4.0 * f64::EPSILON * de, "e {e}: {h}");
    }
    // At declination 60°, an RA difference e is e·cos 60° on the sky
    // (to first order).
    let e = (2.0 + 1e-9) - 2.0;
    let d = angle::anglr_sepr(2.0, rad(60.0), 2.0 + e, rad(60.0));
    assert!((d - 0.5 * e).abs() < 1e-22, "{d}");
}

#[test]
fn anglr_sepr_near_180_degrees() {
    // Antipodal points
    assert_close(angle::anglr_sepr(0.0, 0.0, PI, 0.0), PI, 1e-15);
    assert_close(
        angle::anglr_sepr(1.0, rad(90.0), 1.0, rad(-90.0)),
        PI,
        1e-15,
    );
    assert_close(angle::anglr_sepr(0.3, 0.4, 0.3 + PI, -0.4), PI, 1e-15);
    // 1e-9 rad short of antipodal
    let e = 1e-9;
    let d = angle::anglr_sepr(0.0, 0.0, PI - e, 0.0);
    assert!((PI - d - e).abs() < 1e-15, "{}", PI - d);
    let d = angle::anglr_sepr(0.5, 0.3, 0.5 + PI, -0.3 + e);
    assert!((PI - d - e).abs() < 1e-15, "{}", PI - d);
}

#[test]
fn anglr_sepr_invariants() {
    let pts = [
        (0.0, 0.0),
        (1.0, 0.5),
        (-2.0, -1.2),
        (6.0, 1.5),
        (3.1, -0.01),
        (4.0, rad(90.0)),
        (0.2, rad(-90.0)),
    ];
    for &(a1, d1) in &pts {
        for &(a2, d2) in &pts {
            let d = angle::anglr_sepr(a1, d1, a2, d2);
            assert!((0.0..=PI).contains(&d));
            // Symmetric
            assert_close(angle::anglr_sepr(a2, d2, a1, d1), d, 1e-15);
            // Independent of the RA origin and of whole turns
            assert_close(angle::anglr_sepr(a1 + 0.7, d1, a2 + 0.7, d2), d, 2e-15);
            assert_close(
                angle::anglr_sepr(a1 - 2.0 * TWO_PI, d1, a2 + TWO_PI, d2),
                d,
                5e-15,
            );
            // Mirror symmetry in declination
            assert_close(angle::anglr_sepr(a1, -d1, a2, -d2), d, 1e-15);
            // Agrees with the textbook formula 17.1 where that is accurate
            let cos_d = d1.sin() * d2.sin() + d1.cos() * d2.cos() * (a1 - a2).cos();
            assert_close(d.cos(), cos_d, 1e-15);
            // Haversine agrees except right at 180°
            if d < 3.0 {
                assert_close(angle::anglr_sepr_haversine(a1, d1, a2, d2), d, 1e-14);
            }
        }
    }
    // Any point is 90° from the poles' equator, and the poles are 180° apart.
    assert_close(angle::anglr_sepr(1.0, rad(90.0), 5.0, 0.0), PI / 2.0, 1e-15);
}

// Meeus' eq. 17.2 is a small-angle approximation: check it against the
// exact value for close points, including across 0h.
#[test]
fn anglr_sepr_small_approximation() {
    let exact = angle::anglr_sepr(rad(10.0), rad(20.0), rad(10.1), rad(20.05));
    let approx = angle::anglr_sepr_small(rad(10.0), rad(20.0), rad(10.1), rad(20.05));
    assert!(
        (approx - exact).abs() < 1e-7 * exact.max(1.0),
        "{approx} {exact}"
    );

    // Across the 0h/24h boundary
    let exact = angle::anglr_sepr(rad(359.95), rad(-5.0), rad(0.05), rad(-5.02));
    let approx = angle::anglr_sepr_small(rad(359.95), rad(-5.0), rad(0.05), rad(-5.02));
    assert_close(approx, exact, 1e-9);
    assert_close(exact.to_degrees(), 0.1, 2e-3);
}

// ---------------------------------------------------------------------------
// Relative position angle
// ---------------------------------------------------------------------------

// Meeus, ch. 19: for Orion's Belt the angles at ε Ori are C1 = 49.3622° and
// C2 = 123.1209°. C2 is the position angle of ζ Ori relative to ε Ori; C1 is
// 360° minus the position angle of δ Ori relative to ε Ori.
#[test]
fn relative_position_angle_orion_belt() {
    let delta = (rad(hms(5, 32, 0.40)), rad(dms(true, 0, 17, 56.9)));
    let eps = (rad(hms(5, 36, 12.81)), rad(dms(true, 1, 12, 7.0)));
    let zeta = (rad(hms(5, 40, 45.52)), rad(dms(true, 1, 56, 33.3)));

    let c2 = angle::relative_position_angle(zeta.0, zeta.1, eps.0, eps.1);
    assert_close(c2.to_degrees(), 123.1209, 5e-5);

    let p1 = angle::relative_position_angle(delta.0, delta.1, eps.0, eps.1);
    assert_close(360.0 - p1.to_degrees(), 49.3622, 5e-5);
}

#[test]
fn relative_position_angle_cardinal_directions() {
    let (a, d) = (rad(100.0), rad(-30.0));
    let e = 1e-4;
    let north = angle::relative_position_angle(a, d + e, a, d);
    let south = angle::relative_position_angle(a, d - e, a, d);
    let east = angle::relative_position_angle(a + e, d, a, d);
    let west = angle::relative_position_angle(a - e, d, a, d);
    assert_close(north, 0.0, 1e-12);
    assert_close(south, PI, 1e-12);
    assert_close(east, PI / 2.0, 1e-4);
    assert_close(west, 1.5 * PI, 1e-4);
    // Results are in [0, 2π)
    for &p in &[north, south, east, west] {
        assert!((0.0..TWO_PI).contains(&p));
    }
}

#[test]
fn relative_position_angle_poles_and_wrap() {
    // The north celestial pole is due north (P = 0) of every point, the
    // south pole due south, even though tan δ₁ is infinite there.
    for &(a, d) in &[(0.0, 0.0), (2.0, 1.0), (5.0, -1.3)] {
        let p = angle::relative_position_angle(1.0, rad(90.0), a, d);
        assert_close(p, 0.0, 1e-12);
        let p = angle::relative_position_angle(4.0, rad(-90.0), a, d);
        assert_close(p, PI, 1e-12);
    }
    // Wrap-around of right ascension does not matter
    let p1 = angle::relative_position_angle(rad(0.05), rad(10.0), rad(359.9), rad(10.02));
    let p2 = angle::relative_position_angle(rad(360.05), rad(10.0), rad(-0.1), rad(10.02));
    assert_close(p1, p2, 1e-12);
    // Body 1 is east (larger RA) and slightly south of body 2
    assert!(p1 > PI / 2.0 && p1 < PI, "{p1}");
}

// ---------------------------------------------------------------------------
// Sexagesimal conversion
// ---------------------------------------------------------------------------

#[test]
fn deg_frm_dms_signs() {
    assert_close(
        angle::deg_frm_dms(21, 44, 7.0),
        21.0 + 44.0 / 60.0 + 7.0 / 3600.0,
        1e-14,
    );
    let neg = -(13.0 + 47.0 / 60.0 + 22.0 / 3600.0);
    // All parts negative, or only the degrees (Meeus, ch. 1)
    assert_close(angle::deg_frm_dms(-13, -47, -22.0), neg, 1e-14);
    assert_close(angle::deg_frm_dms(-13, 47, 22.0), neg, 1e-14);
    // Angles between 0° and −1°
    let neg = -(32.0 / 60.0 + 41.0 / 3600.0);
    assert_close(angle::deg_frm_dms(0, -32, -41.0), neg, 1e-15);
    assert_close(angle::deg_frm_dms(0, 0, -20.0), -20.0 / 3600.0, 1e-15);
    // Used throughout the crate for constants in arcseconds
    assert_close(angle::deg_frm_dms(0, 0, -0.03302), -0.03302 / 3600.0, 1e-18);
}

// Before the fix (0, -30, 15.0) gave −0°30′ + 15″ = −0.49583°: the leading
// sign of the minutes was not applied to the seconds as it is for degrees.
#[test]
fn deg_frm_dms_negative_minutes_with_zero_degrees() {
    let expected = -(30.0 / 60.0 + 15.0 / 3600.0);
    assert_close(angle::deg_frm_dms(0, -30, 15.0), expected, 1e-15);
    assert_close(angle::deg_frm_dms(0, -30, -15.0), expected, 1e-15);
}

// Before the fix `deg_frm_hms` ignored the sign rule: (-1, 30, 0) was
// −1h + 30m = −0h30m (−7.5°) instead of −1h30m (−22.5°).
#[test]
fn deg_frm_hms_signs() {
    assert_close(angle::deg_frm_hms(14, 15, 39.7), 213.915417, 1e-6);
    assert_close(angle::deg_frm_hms(-1, 30, 0.0), -22.5, 1e-13);
    assert_close(angle::deg_frm_hms(-1, -30, 0.0), -22.5, 1e-13);
    assert_close(angle::deg_frm_hms(0, -30, 0.0), -7.5, 1e-13);
    assert_close(angle::deg_frm_hms(0, 0, -1.0), -15.0 / 3600.0, 1e-15);
}

// `min.abs()` overflowed (and panicked in debug builds) for i64::MIN.
#[test]
fn deg_frm_dms_does_not_panic_on_extreme_input() {
    assert!(angle::deg_frm_dms(-1, i64::MIN, 0.0) < 0.0);
    assert!(angle::deg_frm_hms(-1, i64::MIN, 0.0) < 0.0);
}

#[test]
fn dms_frm_deg_negative_angles() {
    assert_eq!(angle::dms_frm_deg(-0.5), (0, -30, 0.0));
    assert_eq!(angle::dms_frm_deg(-13.5), (-13, -30, 0.0));
    assert_eq!(angle::hms_frm_deg(-7.5), (0, -30, 0.0));
    let (d, m, s) = angle::dms_frm_deg(-(32.0 / 60.0 + 41.0 / 3600.0));
    assert_eq!((d, m), (0, -32));
    assert_close(s, -41.0, 1e-9);
    let (d, m, s) = angle::dms_frm_deg(-1.0 / 3600.0);
    assert_eq!((d, m), (0, 0));
    assert_close(s, -1.0, 1e-12);
    // Zero seconds are +0.0, not −0.0
    assert!(angle::dms_frm_deg(-0.5).2.is_sign_positive());
    // Every non-zero component of a negative angle is negative, and the
    // tuple reads back to the same angle.
    for i in 1..20_000 {
        let x = -(i as f64) * 0.0123;
        let (d, m, s) = angle::dms_frm_deg(x);
        assert!(d <= 0 && m <= 0 && s <= 0.0, "{x}: {d} {m} {s}");
        assert_close(angle::deg_frm_dms(d, m, s), x, 1e-12);
        let (h, m, s) = angle::hms_frm_deg(x);
        assert!(h <= 0 && m <= 0 && s <= 0.0, "{x}: {h} {m} {s}");
        assert_close(angle::deg_frm_hms(h, m, s), x, 1e-12);
    }
}

// Before the fix 1°01′00″ came back as (1, 0, 59.99999999999979), which
// prints as 1°00′60.00″; about 1 in 130 whole-arcsecond angles was affected.
#[test]
fn dms_frm_deg_whole_seconds_split_exactly() {
    assert_eq!(angle::dms_frm_deg(1.0 + 1.0 / 60.0), (1, 1, 0.0));
    assert_eq!(angle::dms_frm_deg(-(1.0 + 1.0 / 60.0)), (-1, -1, 0.0));
    let mut bad = 0;
    for d in 0..360 {
        for m in 0..60 {
            for s in (0..60).step_by(7) {
                let x = angle::deg_frm_dms(d, m, s as f64);
                let (d2, m2, s2) = angle::dms_frm_deg(x);
                if (d2, m2) != (d, m) || (s2 - s as f64).abs() > 1e-9 {
                    bad += 1;
                }
                let (h2, m2, s2) = angle::hms_frm_deg(x * 15.0);
                let back = angle::deg_frm_hms(h2, m2, s2);
                if !(0..60).contains(&m2)
                    || !(0.0..60.0).contains(&s2)
                    || (back - x * 15.0).abs() > 1e-9
                {
                    bad += 1;
                }
            }
        }
    }
    assert_eq!(bad, 0);
}

#[test]
fn dms_frm_deg_ranges() {
    for i in 0..200_000 {
        let x = i as f64 * 0.001_739_1;
        let (d, m, s) = angle::dms_frm_deg(x);
        assert!(
            d >= 0 && (0..60).contains(&m) && (0.0..60.0).contains(&s),
            "{x}"
        );
        let (h, m, s) = angle::hms_frm_deg(x);
        assert!(
            h >= 0 && (0..60).contains(&m) && (0.0..60.0).contains(&s),
            "{x}"
        );
    }
    // Non-finite input does not panic
    let _ = angle::dms_frm_deg(f64::NAN);
    let _ = angle::dms_frm_deg(f64::INFINITY);
    let _ = angle::hms_frm_deg(f64::NEG_INFINITY);
}

#[test]
fn dms_struct_round_trip_and_sign() {
    let a = Dms::new(true, 0, 32, 41.0);
    assert_close(a.to_deg(), -(32.0 / 60.0 + 41.0 / 3600.0), 1e-16);
    let b = Dms::from_deg(a.to_deg());
    assert!(b.negative);
    assert_eq!((b.deg, b.min), (0, 32));
    assert_close(b.sec, 41.0, 1e-9);
    assert_eq!(Dms::from_deg(-0.5), Dms::new(true, 0, 30, 0.0));
    assert_eq!(Dms::from_deg(0.0), Dms::new(false, 0, 0, 0.0));
    assert_eq!(Dms::from_deg(-0.0), Dms::new(false, 0, 0, 0.0));
    // The sign of `sec` itself is ignored
    assert_eq!(
        Dms::new(false, 1, 2, -3.0).to_deg(),
        Dms::new(false, 1, 2, 3.0).to_deg()
    );

    for i in -5000..5000 {
        let x = i as f64 * 0.071_3;
        let r = Dms::from_deg(x);
        assert!(r.min < 60 && (0.0..60.0).contains(&r.sec));
        assert_eq!(r.negative, x < 0.0);
        assert_close(r.to_deg(), x, 1e-12);
        let h = Hms::from_deg(x);
        assert!(h.min < 60 && (0.0..60.0).contains(&h.sec));
        assert_close(h.to_deg(), x, 1e-12);
    }
}

#[test]
fn dms_rounding_carries() {
    // 10°29′59.96″ to one decimal is 10°30′00.0″, not 10°29′60.0″
    let x = dms(false, 10, 29, 59.96);
    let r = Dms::from_deg_rounded(x, 1);
    assert_eq!((r.deg, r.min), (10, 30));
    assert_close(r.sec, 0.0, 1e-12);
    // Carry through the degrees
    let r = Dms::from_deg_rounded(dms(true, 359, 59, 59.6), 0);
    assert_eq!(r, Dms::new(true, 360, 0, 0.0));
    // Rounded to zero: not negative
    assert_eq!(Dms::from_deg_rounded(-1e-9, 3), Dms::new(false, 0, 0, 0.0));
    // Ordinary rounding
    let r = Dms::from_deg_rounded(dms(false, 1, 2, 3.456), 2);
    assert_eq!((r.deg, r.min), (1, 2));
    assert_close(r.sec, 3.46, 1e-12);
    let h = Hms::from_deg_rounded(hms(23, 59, 59.999), 2);
    assert_eq!(h, Hms::new(false, 24, 0, 0.0));
}

#[test]
fn dms_display() {
    assert_eq!(
        format!("{:.1}", Dms::from_deg(dms(true, 0, 32, 41.0))),
        "-0°32′41.0″"
    );
    assert_eq!(
        format!("{:.0}", Dms::from_deg(dms(false, 23, 26, 44.0))),
        "23°26′44″"
    );
    assert_eq!(
        format!("{:.2}", Dms::from_deg(dms(false, 10, 29, 59.996))),
        "10°30′00.00″"
    );
    assert_eq!(format!("{}", Dms::new(false, 5, 7, 9.5)), "5°07′09.5″");
    assert_eq!(
        format!("{:.3}", Hms::from_deg(hms(14, 15, 39.7))),
        "14h15m39.700s"
    );
    assert_eq!(
        format!("{:.1}", Hms::from_deg(-hms(0, 30, 0.0))),
        "-0h30m00.0s"
    );
}

// ---------------------------------------------------------------------------
// Normalisation
// ---------------------------------------------------------------------------

// Before the fix a tiny negative angle gave exactly 360 (or 2π), outside
// the documented range.
#[test]
fn limit_tiny_negative_angles_stay_in_range() {
    for &x in &[-1e-20, -1e-300, -5e-324, -1e-14] {
        let r = angle::limit_to_360(x);
        assert!((0.0..360.0).contains(&r), "{x} -> {r}");
        let r = angle::limit_to_two_PI(x);
        assert!((0.0..TWO_PI).contains(&r), "{x} -> {r}");
    }
    assert_close(angle::limit_to_360(-1e-10), 360.0 - 1e-10, 1e-13);
}

// Before the fix 1e25° was returned unchanged (i64 saturation) and 1e20°
// gave 0 instead of 280.
#[test]
fn limit_huge_angles_exactly() {
    // 1e25 and 1e20 are exact integers in f64: 10000000000000000905969664
    // and 100000000000000000000, which leave 64 and 280 modulo 360.
    assert_eq!(angle::limit_to_360(1e25), 64.0);
    assert_eq!(angle::limit_to_360(-1e25), 296.0);
    assert_eq!(angle::limit_to_360(1e20), 280.0);
    assert_eq!(angle::limit_to_360(36_000_030.0), 30.0);
    for &x in &[1e25, -1e25, 1e300, -3.7e200, 1e17 + 4.0] {
        let r = angle::limit_to_two_PI(x);
        assert!((0.0..TWO_PI).contains(&r), "{x} -> {r}");
        let r = angle::limit_to_360(x);
        assert!((0.0..360.0).contains(&r), "{x} -> {r}");
    }
    // sin(36000030°) = 0.5 (Meeus, ch. 1)
    assert_close(
        angle::limit_to_360(36_000_030.0).to_radians().sin(),
        0.5,
        1e-15,
    );
}

#[test]
fn limit_exact_multiples_and_signs() {
    for k in -1000..=1000 {
        let x = 360.0 * k as f64;
        assert_eq!(angle::limit_to_360(x), 0.0, "{x}");
        assert!(angle::limit_to_360(x).is_sign_positive());
        let x = TWO_PI * k as f64;
        for &y in &[x, neighbours(x).0, neighbours(x).1] {
            let r = angle::limit_to_two_PI(y);
            assert!((0.0..TWO_PI).contains(&r), "{y} -> {r}");
        }
    }
    assert_eq!(angle::limit_to_360(-0.0), 0.0);
    assert!(angle::limit_to_360(-0.0).is_sign_positive());
    assert_eq!(angle::limit_to_360(-90.0), 270.0);
    assert_eq!(angle::limit_to_360(-450.0), 270.0);
    assert_eq!(angle::limit_to_360(720.5), 0.5);
    assert_close(angle::limit_to_360(-2318.19280), 201.80720, 1e-10);
    assert_close(angle::limit_to_two_PI(-PI / 2.0), 1.5 * PI, 1e-15);
    assert_close(angle::limit_to_two_PI(7.0 * PI), PI, 1e-14);
    assert!(angle::limit_to_360(f64::NAN).is_nan());
    assert!(angle::limit_to_360(f64::INFINITY).is_nan());
}

#[test]
fn limit_to_pm180_and_pm_pi() {
    assert_eq!(angle::limit_to_pm180(180.0), -180.0);
    assert_eq!(angle::limit_to_pm180(-180.0), -180.0);
    assert_eq!(angle::limit_to_pm180(179.5), 179.5);
    assert_eq!(angle::limit_to_pm180(190.0), -170.0);
    assert_eq!(angle::limit_to_pm180(-190.0), 170.0);
    assert_eq!(angle::limit_to_pm180(359.0), -1.0);
    assert_eq!(angle::limit_to_pm180(-721.0), -1.0);
    assert_eq!(angle::limit_to_pm180(1e25), 64.0);
    // Small angles keep their digits
    assert_eq!(angle::limit_to_pm180(1e-20), 1e-20);
    assert_eq!(angle::limit_to_pm180(-1e-20), -1e-20);
    assert_eq!(angle::limit_to_pm_pi(-1e-20), -1e-20);
    assert_eq!(angle::limit_to_pm_pi(PI), -PI);
    assert_close(angle::limit_to_pm_pi(TWO_PI - 0.1), -0.1, 1e-15);
    for i in -3000..3000 {
        let x = i as f64 * 0.37;
        let r = angle::limit_to_pm180(x);
        assert!((-180.0..180.0).contains(&r));
        assert_close(angle::limit_to_360(r), angle::limit_to_360(x), 1e-9);
        let r = angle::limit_to_pm_pi(x);
        assert!((-PI..PI).contains(&r));
    }
}

// Review R1: rounding a Dms/Hms whose seconds are infinite or huge used to
// overflow the minute counter (a panic in debug builds), and the carry
// could leave seconds outside [0, 60).
#[test]
fn dms_rounding_extreme_seconds_does_not_panic() {
    for sec in [
        f64::INFINITY,
        f64::NEG_INFINITY,
        1e300,
        -1e300,
        1e20,
        f64::NAN,
    ] {
        let r = Dms::new(false, 3, 59, sec).rounded(2);
        assert!(r.sec.is_nan() || r.sec.is_infinite() || (0.0..60.0).contains(&r.sec));
        let _ = format!("{:.1}", Dms::new(true, u32::MAX, 59, sec));
        let h = Hms::new(false, 23, 59, sec).rounded(0);
        assert!(h.sec.is_nan() || h.sec.is_infinite() || (0.0..60.0).contains(&h.sec));
        let _ = format!("{:.3}", Hms::new(false, u32::MAX, 59, sec));
    }
    // Seconds beyond a minute carry exactly: 1°59′ + 3675″ = 3°00′15″.
    let r = Dms::new(false, 1, 59, 3600.0 + 75.0).rounded(0);
    assert_eq!(r, Dms::new(false, 3, 0, 15.0));
}

// Review R1: the sign of the seconds passed to `new` is documented as
// ignored, but it was stored, so equal angles compared unequal and `{}`
// printed "1°02′0-3″".
#[test]
fn dms_new_ignores_sign_of_seconds() {
    assert_eq!(Dms::new(false, 1, 2, -3.0), Dms::new(false, 1, 2, 3.0));
    assert_eq!(Hms::new(true, 1, 2, -3.5), Hms::new(true, 1, 2, 3.5));
    assert_eq!(format!("{}", Dms::new(false, 1, 2, -3.0)), "1°02′03″");
    assert_eq!(format!("{}", Hms::new(false, 1, 2, -3.5)), "1h02m03.5s");
    // Even when the field is set directly
    let mut a = Dms::new(false, 1, 2, 3.0);
    a.sec = -3.0;
    assert_eq!(format!("{a}"), "1°02′03″");
}

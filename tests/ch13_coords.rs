// Tests for Meeus, *Astronomical Algorithms*, 2nd ed., chapter 13
// (transformation of coordinates), and GitHub issues #13, #18 and #19.

mod common;

use astro::{angle, coords};
use common::{assert_close, assert_close_deg, assert_close_rad};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

fn rad(deg: f64) -> f64 {
    deg.to_radians()
}

/// A grid of angles (degrees) that includes negative values, values above
/// 360° and the neighbourhood of 0°/360°.
fn angle_grid() -> Vec<f64> {
    let mut v: Vec<f64> = (-36..=72).map(|i| i as f64 * 10.0 + 0.37).collect();
    v.extend_from_slice(&[0.0, 1e-9, -1e-9, 359.999_999, 360.0, -360.0, 720.5]);
    v
}

fn lat_grid() -> Vec<f64> {
    let mut v: Vec<f64> = (-17..=17).map(|i| i as f64 * 5.0 + 0.21).collect();
    v.extend_from_slice(&[0.0, 89.999, -89.999]);
    v
}

// ---------------------------------------------------------------------
// Ecliptic <-> equatorial
// ---------------------------------------------------------------------

// Meeus, example 13.a: Pollux, J2000.0.
#[test]
fn example_13_a_ecl_frm_eq() {
    let asc = rad(116.328942);
    let dec = rad(28.026183);
    let oblq = rad(23.4392911);

    let long = coords::ecl_long_frm_eq(asc, dec, oblq);
    let lat = coords::ecl_lat_frm_eq(asc, dec, oblq);
    assert_close(long.to_degrees(), 113.215630, 1e-6);
    assert_close(lat.to_degrees(), 6.684170, 1e-6);

    // Intermediate value: tan λ = +1.03403986 / −0.44352398.
    assert_close(long.tan(), 1.03403986 / -0.44352398, 1e-7);

    // The same from the sexagesimal input of the example.
    let asc = angle::deg_frm_hms(7, 45, 18.946).to_radians();
    let dec = angle::deg_frm_dms(28, 1, 34.26).to_radians();
    let p = coords::EqPoint { asc, dec }.to_ecl(coords::OBLQ_J2000);
    assert_close(p.long.to_degrees(), 113.215630, 1e-6);
    assert_close(p.lat.to_degrees(), 6.684170, 1e-6);
}

// Meeus, example 13.a, exercise: back to α and δ with eqs. 13.3 and 13.4.
#[test]
fn example_13_a_exercise_eq_frm_ecl() {
    let oblq = rad(23.4392911);
    let asc = coords::asc_frm_ecl(rad(113.215630), rad(6.684170), oblq);
    let dec = coords::dec_frm_ecl(rad(113.215630), rad(6.684170), oblq);
    assert_close(asc.to_degrees(), 116.328942, 1e-6);
    assert_close(dec.to_degrees(), 28.026183, 1e-6);
}

#[test]
fn obliquity_constants() {
    assert_close(coords::OBLQ_J2000.to_degrees(), 23.4392911, 1e-7);
    assert_close(
        coords::OBLQ_J2000.to_degrees(),
        angle::deg_frm_dms(23, 26, 21.448),
        1e-12,
    );
    assert_close(coords::OBLQ_B1950.to_degrees(), 23.4457889, 1e-12);
}

#[test]
fn ecliptic_cardinal_points() {
    let e = coords::OBLQ_J2000;
    // Vernal equinox.
    assert_close_rad(coords::ecl_long_frm_eq(0.0, 0.0, e), 0.0, 1e-15);
    assert_close(coords::ecl_lat_frm_eq(0.0, 0.0, e), 0.0, 1e-15);
    // Summer solstice point: α = 90°, δ = +ε.
    assert_close(coords::asc_frm_ecl(FRAC_PI_2, 0.0, e), FRAC_PI_2, 1e-15);
    assert_close(coords::dec_frm_ecl(FRAC_PI_2, 0.0, e), e, 1e-15);
    // Winter solstice point: α = 270°, δ = −ε.
    assert_close(coords::asc_frm_ecl(1.5 * PI, 0.0, e), 1.5 * PI, 1e-15);
    assert_close(coords::dec_frm_ecl(1.5 * PI, 0.0, e), -e, 1e-15);
    // North ecliptic pole: α = 270°, δ = 90° − ε.
    assert_close(coords::asc_frm_ecl(0.3, FRAC_PI_2, e), 1.5 * PI, 1e-9);
    assert_close(coords::dec_frm_ecl(0.3, FRAC_PI_2, e), FRAC_PI_2 - e, 1e-15);
    // North celestial pole: λ = 90°, β = 90° − ε.
    assert_close(coords::ecl_long_frm_eq(1.0, FRAC_PI_2, e), FRAC_PI_2, 1e-9);
    assert_close(
        coords::ecl_lat_frm_eq(1.0, FRAC_PI_2, e),
        FRAC_PI_2 - e,
        1e-15,
    );
}

// Longitudes and right ascensions are now reduced to [0, 2π); they used to
// come straight out of atan2, in (−π, π].
#[test]
fn ecl_eq_ranges() {
    let e = coords::OBLQ_J2000;
    let long = coords::ecl_long_frm_eq(rad(-10.0), rad(-5.0), e);
    assert!((0.0..TAU).contains(&long), "{long}");
    let asc = coords::asc_frm_ecl(rad(-10.0), rad(3.0), e);
    assert!((0.0..TAU).contains(&asc), "{asc}");

    for &a in &angle_grid() {
        for &d in &lat_grid() {
            let l = coords::ecl_long_frm_eq(rad(a), rad(d), e);
            let b = coords::ecl_lat_frm_eq(rad(a), rad(d), e);
            assert!((0.0..TAU).contains(&l));
            assert!((-FRAC_PI_2..=FRAC_PI_2).contains(&b));
            let al = coords::asc_frm_ecl(rad(a), rad(d), e);
            let de = coords::dec_frm_ecl(rad(a), rad(d), e);
            assert!((0.0..TAU).contains(&al));
            assert!((-FRAC_PI_2..=FRAC_PI_2).contains(&de));
        }
    }
}

#[test]
fn ecl_eq_round_trip() {
    for &oblq in &[coords::OBLQ_J2000, coords::OBLQ_B1950, 0.0, rad(-23.0)] {
        for &a in &angle_grid() {
            for &d in &lat_grid() {
                let (l, b) = astro::ecl_frm_eq!(rad(a), rad(d), oblq);
                let (a2, d2) = astro::eq_frm_ecl!(l, b, oblq);
                assert_close_rad(a2, rad(a), 1e-9);
                assert_close(d2, rad(d), 1e-12);
            }
        }
    }
}

// A zero obliquity makes the two systems identical.
#[test]
fn ecl_eq_zero_obliquity() {
    for &a in &angle_grid() {
        for &d in &lat_grid() {
            assert_close_rad(coords::ecl_long_frm_eq(rad(a), rad(d), 0.0), rad(a), 1e-12);
            assert_close(coords::ecl_lat_frm_eq(rad(a), rad(d), 0.0), rad(d), 1e-14);
        }
    }
}

#[test]
fn ecl_eq_point_methods_match_functions() {
    let e = coords::OBLQ_J2000;
    let eq = coords::EqPoint {
        asc: rad(200.0),
        dec: rad(-33.0),
    };
    let ecl = eq.to_ecl(e);
    assert_eq!((ecl.long, ecl.lat), astro::ecl_frm_eq!(eq.asc, eq.dec, e));
    let back = ecl.to_eq(e);
    assert_eq!(
        (back.asc, back.dec),
        astro::eq_frm_ecl!(ecl.long, ecl.lat, e)
    );
    assert_close_rad(back.asc, eq.asc, 1e-14);
    assert_close(back.dec, eq.dec, 1e-14);
}

// ---------------------------------------------------------------------
// Hour angle and sidereal time (issue #18)
// ---------------------------------------------------------------------

// Meeus, example 13.b: H = θ₀ − L − α with L positive west.
#[test]
fn example_13_b_hour_angle() {
    let theta0 = angle::deg_frm_hms(8, 34, 56.853).to_radians();
    let long = angle::deg_frm_hms(5, 8, 15.7).to_radians();
    let asc = angle::deg_frm_hms(23, 9, 16.641).to_radians();

    let h = coords::hr_angl_frm_observer_long(theta0, long, asc);
    assert_close(h.to_degrees(), 64.352133, 1e-6);

    let theta = coords::loc_sidr_frm_green_sidr(theta0, long);
    assert_close(coords::hr_angl_frm_loc_sidr(theta, asc), h, 1e-12);
}

// Issue #18 reported that the hour angle should be θ₀ + L − α. That is
// the formula for east-positive longitudes; this crate follows Meeus,
// whose longitudes are positive west, so θ₀ − L − α is right here.
#[test]
fn hour_angle_longitude_convention_issue_18() {
    // An observer 30° EAST of Greenwich (L = −30° in this crate) sees a star
    // that is on the Greenwich meridian (α = θ₀) 30° past their own
    // meridian, i.e. at H = +30° (west of the meridian, already set
    // towards the west).
    let theta0 = rad(100.0);
    let asc = rad(100.0);
    let long_west = rad(-30.0);
    let h = coords::hr_angl_frm_observer_long(theta0, long_west, asc);
    assert_close(h.to_degrees(), 30.0, 1e-12);

    let long_east = coords::east_long_frm_west_long(long_west);
    assert_close(long_east.to_degrees(), 30.0, 1e-12);
    let h_east = coords::hr_angl_frm_observer_east_long(theta0, long_east, asc);
    assert_close(h_east, h, 1e-15);
    assert_close(coords::west_long_frm_east_long(long_east), long_west, 0.0);

    // Local sidereal time is ahead of Greenwich's for an eastern observer.
    let theta = coords::loc_sidr_frm_green_sidr(theta0, long_west);
    assert_close(theta.to_degrees(), 130.0, 1e-12);
}

#[test]
fn hour_angle_ranges() {
    for &t in &angle_grid() {
        for &l in &[-179.0, -77.0, 0.0, 12.5, 179.9] {
            for &a in &[0.0, 90.0, 359.9, -45.0] {
                let h = coords::hr_angl_frm_observer_long(rad(t), rad(l), rad(a));
                assert!(h > -PI && h <= PI);
                assert_close_rad(h, rad(t - l - a), 1e-12);
                let he = coords::hr_angl_frm_observer_east_long(rad(t), rad(-l), rad(a));
                assert_close_rad(he, h, 1e-12);
                let th = coords::loc_sidr_frm_green_sidr(rad(t), rad(l));
                assert!((0.0..TAU).contains(&th));
                let hl = coords::hr_angl_frm_loc_sidr(th, rad(a));
                assert!(hl > -PI && hl <= PI);
                assert_close_rad(hl, h, 1e-12);
            }
        }
    }
}

// ---------------------------------------------------------------------
// Equatorial <-> horizontal (issues #13 and #19)
// ---------------------------------------------------------------------

// Meeus, example 13.b: Venus from Washington, H = 64.352133°.
#[test]
fn example_13_b_hz_frm_eq() {
    let h = rad(64.352133);
    let dec = angle::deg_frm_dms(-6, 43, 11.61).to_radians();
    let lat = angle::deg_frm_dms(38, 55, 17.0).to_radians();

    let az = coords::az_frm_eq(h, dec, lat);
    let alt = coords::alt_frm_eq(h, dec, lat);
    assert_close(az.to_degrees(), 68.0337, 5e-5);
    assert_close(alt.to_degrees(), 15.1249, 5e-5);
    // Intermediate value: tan A = +0.9014712 / +0.3636015.
    assert_close(az.tan(), 0.9014712 / 0.3636015, 1e-5);

    let p = coords::hz_frm_hr_angl_and_dec(h, dec, lat);
    assert_eq!((p.az, p.alt), astro::loc_hz_frm_eq!(h, dec, lat));
    assert_close(p.az, az, 1e-15);
    assert_close(p.alt, alt, 1e-15);

    // And back.
    let q = coords::hr_angl_and_dec_frm_hz(az, alt, lat);
    assert_close(q.hour_angle, h, 1e-12);
    assert_close(q.dec, dec, 1e-12);
    assert_close(coords::hr_angl_frm_hz(az, alt, lat), h, 1e-12);
    assert_close(coords::dec_frm_hz(az, alt, lat), dec, 1e-12);
}

// GitHub issue #19: the declination computed from horizontal coordinates
// used cos(A)·cos(A) where Meeus has cos(h)·cos(A).
#[test]
fn dec_frm_hz_issue_19() {
    let lat = rad(38.92);
    let (h, dec) = (rad(64.352133), rad(-6.72));
    let az = coords::az_frm_eq(h, dec, lat);
    let alt = coords::alt_frm_eq(h, dec, lat);
    assert_close(coords::dec_frm_hz(az, alt, lat), dec, 1e-12);

    // Direct check of sin δ = sin φ sin h − cos φ cos h cos A.
    let (az, alt) = (rad(-150.0), rad(20.0));
    let expected = (lat.sin() * alt.sin() - lat.cos() * alt.cos() * az.cos()).asin();
    assert_close(coords::dec_frm_hz(az, alt, lat), expected, 1e-15);
}

/// An independent conversion in the compass convention (azimuth from the
/// north, towards the east), as used by ERFA/SOFA `ae2hd`/`hd2ae`: the
/// direction is expressed in the east-north-up frame of the observer and
/// projected on the axes of the hour-angle frame.
fn hd_frm_ae_north(az_north: f64, alt: f64, lat: f64) -> (f64, f64) {
    let e = alt.cos() * az_north.sin();
    let n = alt.cos() * az_north.cos();
    let u = alt.sin();
    // Axes of the hour-angle frame in ENU: towards (H = 0, δ = 0), towards
    // the west point (H = 90°), towards the north celestial pole.
    let x = -n * lat.sin() + u * lat.cos();
    let y = -e;
    let z = n * lat.cos() + u * lat.sin();
    (y.atan2(x), z.atan2(x.hypot(y)))
}

fn ae_north_frm_hd(h: f64, dec: f64, lat: f64) -> (f64, f64) {
    let x = dec.cos() * h.cos();
    let y = dec.cos() * h.sin();
    let z = dec.sin();
    let e = -y;
    let n = -x * lat.sin() + z * lat.cos();
    let u = x * lat.cos() + z * lat.sin();
    (e.atan2(n).rem_euclid(TAU), u.atan2(e.hypot(n)))
}

// GitHub issue #13: compare (az, alt) -> (H, δ) with an independent
// formulation in the north-based azimuth convention. Once the azimuth is
// converted between conventions (they differ by 180°) and issue #19 is
// fixed, the two agree to rounding error everywhere.
#[test]
fn hz_to_eq_matches_north_based_formulation_issue_13() {
    for &lat in &lat_grid() {
        for az_n in (0..36).map(|i| i as f64 * 10.0 + 3.3) {
            for alt in (-17..=17).map(|i| i as f64 * 5.0 + 1.1) {
                let (h_ref, dec_ref) = hd_frm_ae_north(rad(az_n), rad(alt), rad(lat));
                let az_s = coords::az_south_frm_north(rad(az_n));

                let p = coords::hr_angl_and_dec_frm_hz(az_s, rad(alt), rad(lat));
                assert_close_rad(p.hour_angle, h_ref, 1e-11);
                assert_close(p.dec, dec_ref, 1e-12);
                assert!(p.hour_angle > -PI && p.hour_angle <= PI);

                assert_close_rad(
                    coords::hr_angl_frm_hz(az_s, rad(alt), rad(lat)),
                    h_ref,
                    1e-11,
                );
                assert_close(coords::dec_frm_hz(az_s, rad(alt), rad(lat)), dec_ref, 1e-12);
            }
        }
    }
}

#[test]
fn eq_to_hz_matches_north_based_formulation() {
    for &lat in &lat_grid() {
        for &h in &angle_grid() {
            for dec in (-17..=17).map(|i| i as f64 * 5.0 + 0.7) {
                let (az_ref, alt_ref) = ae_north_frm_hd(rad(h), rad(dec), rad(lat));
                let p = coords::hz_frm_hr_angl_and_dec(rad(h), rad(dec), rad(lat));
                assert_close_rad(coords::az_north_frm_south(p.az), az_ref, 1e-11);
                assert_close(p.alt, alt_ref, 1e-12);
                assert!(p.az > -PI && p.az <= PI);

                let az = coords::az_frm_eq(rad(h), rad(dec), rad(lat));
                let alt = coords::alt_frm_eq(rad(h), rad(dec), rad(lat));
                assert_close(az, p.az, 1e-15);
                assert_close(alt, p.alt, 1e-15);
            }
        }
    }
}

#[test]
fn eq_hz_round_trip() {
    for &lat in &lat_grid() {
        for &h in &angle_grid() {
            for dec in (-17..=17).map(|i| i as f64 * 5.0 + 0.7) {
                let hz = coords::hz_frm_hr_angl_and_dec(rad(h), rad(dec), rad(lat));
                let back = coords::hr_angl_and_dec_frm_hz(hz.az, hz.alt, rad(lat));
                assert_close_rad(back.hour_angle, rad(h), 1e-10);
                assert_close(back.dec, rad(dec), 1e-12);
            }
        }
    }
}

#[test]
fn hz_special_directions() {
    let lat = rad(52.0);
    // On the meridian, south of the zenith: A = 0, h = 90° − φ + δ.
    let p = coords::hz_frm_hr_angl_and_dec(0.0, rad(10.0), lat);
    assert_close(p.az, 0.0, 1e-15);
    assert_close(p.alt.to_degrees(), 48.0, 1e-12);
    // Lower culmination of a circumpolar star: due north, h = φ + δ − 90°.
    let p = coords::hz_frm_hr_angl_and_dec(PI, rad(60.0), lat);
    assert_close_rad(p.az, PI, 1e-15);
    assert_close(p.alt.to_degrees(), 22.0, 1e-12);
    // The north celestial pole: due north at altitude φ.
    let p = coords::hz_frm_hr_angl_and_dec(rad(37.0), FRAC_PI_2, lat);
    assert_close_rad(coords::az_north_frm_south(p.az), 0.0, 1e-12);
    assert_close(p.alt, lat, 1e-15);
    // A body on the zenith.
    assert_close(coords::alt_frm_eq(0.0, lat, lat), FRAC_PI_2, 1e-15);
    // On the equator, a body of δ = 0 rises due east at H = −90°.
    let p = coords::hz_frm_hr_angl_and_dec(rad(-90.0), 0.0, 0.0);
    assert_close(coords::az_north_frm_south(p.az).to_degrees(), 90.0, 1e-12);
    assert_close(p.alt, 0.0, 1e-15);
    // Southern observer: the south celestial pole is due south at −φ.
    let p = coords::hz_frm_hr_angl_and_dec(rad(-120.0), -FRAC_PI_2, rad(-33.9));
    assert_close_rad(p.az, 0.0, 1e-12);
    assert_close(p.alt.to_degrees(), 33.9, 1e-12);
}

#[test]
fn azimuth_conventions() {
    // (south-based, north-based) pairs: S, W, N, E.
    for &(s, n) in &[
        (0.0, 180.0),
        (90.0, 270.0),
        (180.0, 0.0),
        (270.0, 90.0),
        (-90.0, 90.0),
    ] {
        assert_close_deg(coords::az_north_frm_south(rad(s)).to_degrees(), n, 1e-12);
        assert_close_deg(coords::az_south_frm_north(rad(n)).to_degrees(), s, 1e-12);
    }
    for &a in &angle_grid() {
        let n = coords::az_north_frm_south(rad(a));
        assert!((0.0..TAU).contains(&n));
        let s = coords::az_south_frm_north(n);
        assert!(s > -PI && s <= PI);
        assert_close_rad(s, rad(a), 1e-12);
    }
}

// ---------------------------------------------------------------------
// Galactic coordinates (B1950.0)
// ---------------------------------------------------------------------

// Meeus, chapter 13, exercise: Nova Serpentis 1978.
#[test]
fn nova_serpentis_gal_frm_eq() {
    let asc = angle::deg_frm_hms(17, 48, 59.74).to_radians();
    let dec = angle::deg_frm_dms(-14, 43, 8.2).to_radians();

    let l = coords::gal_long_frm_eq(asc, dec);
    let b = coords::gal_lat_frm_eq(asc, dec);
    // The longitude used to come out as 372.9593° (not reduced).
    assert_close(l.to_degrees(), 12.9593, 5e-5);
    assert_close(b.to_degrees(), 6.0463, 5e-5);
    assert_eq!((l, b), astro::gal_frm_eq!(asc, dec));

    let (a2, d2) = astro::eq_frm_gal!(l, b);
    assert_close(a2, asc, 1e-12);
    assert_close(d2, dec, 1e-12);
}

#[test]
fn galactic_reference_points() {
    // North galactic pole.
    assert_close(
        coords::gal_lat_frm_eq(rad(192.25), rad(27.4)).to_degrees(),
        90.0,
        1e-6,
    );
    assert_close(
        coords::dec_frm_gal(1.0, FRAC_PI_2).to_degrees(),
        27.4,
        1e-12,
    );
    assert_close(
        coords::asc_frm_gal(1.0, FRAC_PI_2).to_degrees(),
        192.25,
        1e-6,
    );
    // The north celestial pole is at galactic longitude 123°.
    assert_close(
        coords::gal_long_frm_eq(0.4, FRAC_PI_2).to_degrees(),
        123.0,
        1e-6,
    );
    // Ascending node of the galactic equator on the B1950.0 equator.
    assert_close(
        coords::asc_frm_gal(rad(33.0), 0.0).to_degrees(),
        282.25,
        1e-12,
    );
    assert_close(coords::dec_frm_gal(rad(33.0), 0.0), 0.0, 1e-15);
    // Galactic centre of the 1959 system: α ≈ 17ʰ42.4ᵐ, δ ≈ −28°55′ (B1950).
    let (a, d) = astro::eq_frm_gal!(0.0, 0.0);
    assert_close(a.to_degrees(), angle::deg_frm_hms(17, 42, 24.0), 0.05);
    assert_close(d.to_degrees(), angle::deg_frm_dms(-28, 55, 0.0), 0.05);
}

#[test]
fn galactic_round_trip_and_ranges() {
    for &a in &angle_grid() {
        for &d in &lat_grid() {
            let (l, b) = astro::gal_frm_eq!(rad(a), rad(d));
            assert!((0.0..TAU).contains(&l));
            assert!((-FRAC_PI_2..=FRAC_PI_2).contains(&b));
            let (a2, d2) = astro::eq_frm_gal!(l, b);
            assert!((0.0..TAU).contains(&a2));
            assert_close_rad(a2, rad(a), 1e-9);
            assert_close(d2, rad(d), 1e-12);
        }
    }
}

// ---------------------------------------------------------------------
// Point structs
// ---------------------------------------------------------------------

#[test]
fn point_structs_are_copy_and_comparable() {
    let p = coords::GeographPoint {
        long: 0.1,
        lat: 0.2,
    };
    let q = p;
    assert_eq!(p, q);
    let e = coords::EqPoint { asc: 1.0, dec: 0.5 };
    let f = e;
    assert_eq!(e, f);
    assert_close(
        e.anglr_sepr(&coords::EqPoint { asc: 1.0, dec: 0.0 }),
        0.5,
        1e-15,
    );
    let c = coords::EclPoint {
        long: 2.0,
        lat: -0.3,
    };
    assert_close(c.anglr_sepr(&c), 0.0, 1e-7);
}

// Review R1: the conversion macros expanded each argument twice, so an
// argument with side effects (or an expensive one) was evaluated twice.
#[test]
fn macros_evaluate_each_argument_once() {
    use std::cell::Cell;
    let n = Cell::new(0);
    let arg = |x: f64| {
        n.set(n.get() + 1);
        x
    };
    let e = coords::OBLQ_J2000;

    let (l, b) = astro::ecl_frm_eq!(arg(1.0), arg(0.3), arg(e));
    assert_eq!(n.replace(0), 3);
    assert_eq!(
        (l, b),
        (
            coords::ecl_long_frm_eq(1.0, 0.3, e),
            coords::ecl_lat_frm_eq(1.0, 0.3, e)
        )
    );

    let (a, d) = astro::eq_frm_ecl!(arg(1.0), arg(0.3), arg(e));
    assert_eq!(n.replace(0), 3);
    assert_eq!(
        (a, d),
        (
            coords::asc_frm_ecl(1.0, 0.3, e),
            coords::dec_frm_ecl(1.0, 0.3, e)
        )
    );

    let (az, alt) = astro::loc_hz_frm_eq!(arg(1.0), arg(0.3), arg(0.7));
    assert_eq!(n.replace(0), 3);
    assert_eq!(
        (az, alt),
        (
            coords::az_frm_eq(1.0, 0.3, 0.7),
            coords::alt_frm_eq(1.0, 0.3, 0.7)
        )
    );

    let (l, b) = astro::gal_frm_eq!(arg(1.0), arg(0.3));
    assert_eq!(n.replace(0), 2);
    assert_eq!(
        (l, b),
        (
            coords::gal_long_frm_eq(1.0, 0.3),
            coords::gal_lat_frm_eq(1.0, 0.3)
        )
    );

    let (a, d) = astro::eq_frm_gal!(arg(1.0), arg(0.3));
    assert_eq!(n.replace(0), 2);
    assert_eq!(
        (a, d),
        (coords::asc_frm_gal(1.0, 0.3), coords::dec_frm_gal(1.0, 0.3))
    );
}

/// Unit vector of (longitude, latitude).
fn unit_vec(long: f64, lat: f64) -> [f64; 3] {
    [lat.cos() * long.cos(), lat.cos() * long.sin(), lat.sin()]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// Review R1: eqs. 13.7 and 13.8 against a rotation built only from the
// definition of the IAU 1959 system (Meeus, ch. 13): the north galactic
// pole at α = 192.25°, δ = 27.4° (B1950.0), and the zero of galactic
// longitude 33° along the galactic equator from its ascending node on the
// equator, so that the node has galactic longitude 33°.
#[test]
fn galactic_matches_independent_rotation() {
    let z_gal = unit_vec(rad(192.25), rad(27.4));
    // Ascending node of the galactic equator on the celestial equator:
    // along the line of intersection, the direction in which the galactic
    // equator passes from south to north of the celestial equator.
    let n = cross3([0.0, 0.0, 1.0], z_gal);
    let len = dot3(n, n).sqrt();
    let node = [n[0] / len, n[1] / len, n[2] / len];
    let ninety = cross3(z_gal, node); // 90° further along the galactic equator
    assert!(ninety[2] > 0.0, "the node must be ascending");
    // Galactic longitude 0 lies 33° before the node.
    let (s, c) = rad(33.0).sin_cos();
    let x_gal = [
        c * node[0] - s * ninety[0],
        c * node[1] - s * ninety[1],
        c * node[2] - s * ninety[2],
    ];
    let y_gal = cross3(z_gal, x_gal);

    for &a in &angle_grid() {
        for &d in &lat_grid() {
            let v = unit_vec(rad(a), rad(d));
            let (x, y, z) = (dot3(v, x_gal), dot3(v, y_gal), dot3(v, z_gal));
            let l_ref = y.atan2(x).rem_euclid(TAU);
            let b_ref = z.atan2(x.hypot(y));
            let (l, b) = astro::gal_frm_eq!(rad(a), rad(d));
            if b_ref.abs() < FRAC_PI_2 - 1e-6 {
                assert_close_rad(l, l_ref, 1e-9);
            }
            assert_close(b, b_ref, 1e-12);

            // and the inverse, from the same rotation
            let lv = unit_vec(rad(a), rad(d) * 0.999);
            let w = [
                lv[0] * x_gal[0] + lv[1] * y_gal[0] + lv[2] * z_gal[0],
                lv[0] * x_gal[1] + lv[1] * y_gal[1] + lv[2] * z_gal[1],
                lv[0] * x_gal[2] + lv[1] * y_gal[2] + lv[2] * z_gal[2],
            ];
            let asc_ref = w[1].atan2(w[0]).rem_euclid(TAU);
            let dec_ref = w[2].atan2(w[0].hypot(w[1]));
            let (asc, dec) = astro::eq_frm_gal!(rad(a), rad(d) * 0.999);
            assert_close_rad(asc, asc_ref, 1e-9);
            assert_close(dec, dec_ref, 1e-12);
        }
    }
}

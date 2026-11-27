// Tests for Meeus, *Astronomical Algorithms*, 2nd ed., chapter 11
// (the Earth's globe).

mod common;

use astro::{angle, coords, planet::earth};
use common::assert_close;
use earth::Ellipsoid;

fn point(long_deg: f64, lat_deg: f64) -> coords::GeographPoint {
    coords::GeographPoint {
        long: long_deg.to_radians(),
        lat: lat_deg.to_radians(),
    }
}

fn paris() -> coords::GeographPoint {
    // 2°20′14″ EAST, i.e. negative in the crate's west-positive convention.
    point(
        angle::deg_frm_dms(-2, 20, 14.0),
        angle::deg_frm_dms(48, 50, 11.0),
    )
}

fn washington() -> coords::GeographPoint {
    point(
        angle::deg_frm_dms(77, 3, 56.0),
        angle::deg_frm_dms(38, 55, 17.0),
    )
}

// ---------------------------------------------------------------------
// Ellipsoid constants
// ---------------------------------------------------------------------

// Meeus, chapter 11: derived constants of the IAU 1976 ellipsoid.
#[test]
fn iau_1976_constants() {
    let e = Ellipsoid::IAU_1976;
    assert_close(e.pol_rad() / 1000.0, 6356.755, 5e-4);
    assert_close(e.pol_rad() / e.eq_rad, 0.99664719, 5e-9);
    assert_close(e.ecc(), 0.08181922, 5e-9);
    // 1 − e² = (1 − f)².
    assert_close(1.0 - e.ecc().powi(2), (1.0 - e.flat).powi(2), 1e-15);
}

#[test]
fn wgs84_free_functions() {
    assert_close(earth::eq_rad(), 6378.137, 0.0);
    assert_close(earth::flat_fac(), 1.0 / 298.257223563, 0.0);
    assert_close(earth::pol_rad(), 6356.752314245, 1e-9);
    assert_close(earth::ecc_of_meridian(), 0.0818191908426, 1e-12);
    let w = Ellipsoid::WGS_84;
    assert_close(w.eq_rad, astro::consts::wgs84::EQUATORIAL_RADIUS, 0.0);
    assert_close(w.pol_rad(), astro::consts::wgs84::POLAR_RADIUS, 1e-6);
    assert_close(Ellipsoid::WGS_72.eq_rad, 6_378_135.0, 0.0);
}

// ---------------------------------------------------------------------
// Example 11.a: ρ sin φ′ and ρ cos φ′ for Palomar
// ---------------------------------------------------------------------

#[test]
fn example_11_a_rho_sin_cos_phi() {
    let lat = angle::deg_frm_dms(33, 21, 22.0);
    assert_close(lat, 33.356111, 5e-7);

    let (s, c) = Ellipsoid::IAU_1976.rho_sin_cos_phi(lat.to_radians(), 1706.0);
    assert_close(s, 0.546861, 5e-7);
    assert_close(c, 0.836339, 5e-7);

    // The WGS 84 free function gives the same at this precision.
    let (s, c) = earth::rho_sin_cos_phi(lat.to_radians(), 1706.0);
    assert_close(s, 0.546861, 5e-7);
    assert_close(c, 0.836339, 5e-7);

    // Intermediate: the reduced latitude u = 33.267796°, which is the
    // geocentric direction of the sea-level point scaled by b/a.
    let (s0, c0) = Ellipsoid::IAU_1976.rho_sin_cos_phi(lat.to_radians(), 0.0);
    let u = (s0 / (1.0 - Ellipsoid::IAU_1976.flat)).atan2(c0);
    assert_close(u.to_degrees(), 33.267796, 5e-7);
}

#[test]
fn rho_sin_cos_phi_hemispheres_and_poles() {
    let e = Ellipsoid::IAU_1976;
    for lat in (-90..=90).map(|d| (d as f64).to_radians()) {
        for &h in &[0.0, 1706.0, -400.0, 8848.0] {
            let (s, c) = e.rho_sin_cos_phi(lat, h);
            let (s2, c2) = e.rho_sin_cos_phi(-lat, h);
            assert_close(s, -s2, 1e-15);
            assert_close(c, c2, 1e-15);
            assert!(c >= -1e-16, "ρ cos φ′ must not be negative: {c}");
            assert!(s * lat >= 0.0);
        }
    }
    // At the poles ρ = b/a, at the equator ρ = 1 (sea level).
    let (s, c) = e.rho_sin_cos_phi(90_f64.to_radians(), 0.0);
    assert_close(s, 1.0 - e.flat, 1e-15);
    assert_close(c, 0.0, 1e-15);
    let (s, c) = e.rho_sin_cos_phi(0.0, 0.0);
    assert_close(s, 0.0, 0.0);
    assert_close(c, 1.0, 0.0);
    // Height adds along the vertical: H / a at the equator.
    let (_, c) = e.rho_sin_cos_phi(0.0, 6378.14);
    assert_close(c, 1.001, 1e-15);
}

// ---------------------------------------------------------------------
// Geocentric latitude and ρ
// ---------------------------------------------------------------------

// Meeus, chapter 11: φ − φ′ = 692.73″ sin 2φ − 1.16″ sin 4φ, in radians.
// It used to be returned in degrees.
#[test]
fn geograph_geocent_lat_diff_in_radians() {
    let d = earth::geograph_geocent_lat_diff(45_f64.to_radians());
    assert_close(d, (692.73_f64 / 3600.0).to_radians(), 1e-15);
    assert_close(earth::geograph_geocent_lat_diff(0.0), 0.0, 1e-18);
}

#[test]
fn geocent_lat_series_matches_exact() {
    let e = Ellipsoid::IAU_1976;
    for lat in (-900..=900).map(|d| (d as f64 / 10.0).to_radians()) {
        let exact = lat - e.geocent_lat(lat);
        let series = earth::geograph_geocent_lat_diff(lat);
        assert_close(
            series.to_degrees() * 3600.0,
            exact.to_degrees() * 3600.0,
            0.01,
        );
    }
}

#[test]
fn geocent_lat_properties() {
    let e = Ellipsoid::WGS_84;
    assert_close(e.geocent_lat(0.0), 0.0, 0.0);
    assert_close(
        e.geocent_lat(90_f64.to_radians()),
        90_f64.to_radians(),
        1e-15,
    );
    assert_close(
        e.geocent_lat(-90_f64.to_radians()),
        -90_f64.to_radians(),
        1e-15,
    );
    for lat in (-89..=89).map(|d| (d as f64 + 0.5).to_radians()) {
        let g = e.geocent_lat(lat);
        assert!(
            g.abs() < lat.abs(),
            "|φ′| < |φ| except at poles and equator"
        );
        assert_close(e.geocent_lat(-lat), -g, 1e-15);
        // Same as the direction of (ρ cos φ′, ρ sin φ′) at sea level.
        let (s, c) = e.rho_sin_cos_phi(lat, 0.0);
        assert_close(s.atan2(c), g, 1e-15);
    }
}

// Meeus, chapter 11: the greatest difference φ₀ − φ₀′ = 11′32.73″, with
// φ₀ = 45°05′46.36″ and φ₀′ = 44°54′13.64″ (IAU 1976).
#[test]
fn max_geocent_diff() {
    let (g, c) = Ellipsoid::IAU_1976.lat_of_max_geocent_diff();
    let arcsec = |r: f64| r.to_degrees() * 3600.0;
    assert_close(
        arcsec(g),
        arcsec(angle::deg_frm_dms(45, 5, 46.36).to_radians()),
        0.005,
    );
    assert_close(
        arcsec(c),
        arcsec(angle::deg_frm_dms(44, 54, 13.64).to_radians()),
        0.005,
    );
    assert_close(arcsec(g - c), 11.0 * 60.0 + 32.73, 0.005);
    assert_close(g + c, 90_f64.to_radians(), 1e-15);
    // It is a maximum of φ − φ′.
    let e = Ellipsoid::IAU_1976;
    let diff = |p: f64| p - e.geocent_lat(p);
    assert!(diff(g) > diff(g + 1e-3) && diff(g) > diff(g - 1e-3));
    assert_close(Ellipsoid::IAU_1976.geocent_lat(g), c, 1e-15);
}

#[test]
fn rho_series_matches_exact() {
    let e = Ellipsoid::IAU_1976;
    for lat in (-90..=90).map(|d| (d as f64).to_radians()) {
        assert_close(earth::rho(lat), e.rho(lat), 2e-7);
    }
    assert_close(e.rho(0.0), 1.0, 1e-15);
    assert_close(e.rho(90_f64.to_radians()), 1.0 - e.flat, 1e-15);
}

// ---------------------------------------------------------------------
// Example 11.b: Chicago, φ = +42°
// ---------------------------------------------------------------------

#[test]
fn example_11_b_iau_1976() {
    let e = Ellipsoid::IAU_1976;
    let lat = 42_f64.to_radians();
    assert_close(e.rad_of_parll_lat(lat) / 1000.0, 4747.001, 5e-4);
    assert_close(e.len_of_deg_long(lat) / 1000.0, 82.8508, 5e-5);
    assert_close(e.linear_velocity_at_lat(lat) / 1000.0, 0.34616, 5e-6);
    assert_close(e.rad_curv_of_meridian(lat) / 1000.0, 6364.033, 5e-4);
    assert_close(e.len_of_deg_lat(lat) / 1000.0, 111.0733, 5e-5);
}

// The WGS 84 free functions, in km, differ from the book by ~2-3 m.
#[test]
fn example_11_b_wgs84() {
    let lat = 42_f64.to_radians();
    assert_close(earth::rad_of_parll_lat(lat), 4747.001, 0.003);
    assert_close(earth::linear_velocity_at_lat(lat), 0.34616, 5e-6);
    assert_close(earth::rad_curv_of_meridian(lat), 6364.033, 0.004);
    assert_close(
        earth::rad_of_parll_lat(lat),
        Ellipsoid::WGS_84.rad_of_parll_lat(lat) / 1000.0,
        1e-12,
    );
}

// Meeus, chapter 11: R_m ranges from a(1 − e²) = 6335.44 km at the equator
// to a/√(1 − e²) = 6399.60 km at the poles (IAU 1976).
#[test]
fn rad_curv_of_meridian_extremes() {
    let e = Ellipsoid::IAU_1976;
    assert_close(e.rad_curv_of_meridian(0.0) / 1000.0, 6335.44, 5e-3);
    assert_close(
        e.rad_curv_of_meridian(90_f64.to_radians()) / 1000.0,
        6399.60,
        5e-3,
    );
    assert_close(
        e.rad_curv_of_meridian(-90_f64.to_radians()),
        e.rad_curv_of_meridian(90_f64.to_radians()),
        1e-6,
    );
    let mut prev = 0.0;
    for lat in (0..=90).map(|d| (d as f64).to_radians()) {
        let r = e.rad_curv_of_meridian(lat);
        assert!(r > prev);
        prev = r;
    }
}

#[test]
fn rad_of_parll_lat_limits() {
    let e = Ellipsoid::WGS_84;
    assert_close(e.rad_of_parll_lat(0.0), e.eq_rad, 0.0);
    assert_close(e.rad_of_parll_lat(90_f64.to_radians()), 0.0, 1e-9);
    assert_close(
        e.rad_of_parll_lat(-37_f64.to_radians()),
        e.rad_of_parll_lat(37_f64.to_radians()),
        1e-9,
    );
    // Equator: 465.1 m/s.
    assert_close(e.linear_velocity_at_lat(0.0), 465.10, 0.01);
    assert_close(earth::rot_angular_velocity(), 7.292114992e-5, 0.0);
}

// ---------------------------------------------------------------------
// Example 11.c: Paris to Washington
// ---------------------------------------------------------------------

#[test]
fn example_11_c_geodesic_dist() {
    let s = Ellipsoid::IAU_1976.geodesic_dist(&paris(), &washington());
    assert_close(s / 1000.0, 6181.63, 0.005);

    let s = earth::geodesic_dist(&paris(), &washington());
    assert_close(s, 6181.63, 0.005);
}

#[test]
fn example_11_c_approx_geodesic_dist() {
    let d = paris().anglr_sepr(&washington());
    assert_close(d.cos(), 0.567146, 5e-7);
    assert_close(d.to_degrees(), 55.44855, 5e-6);
    assert_close(
        earth::approx_geodesic_dist(&paris(), &washington()),
        6166.0,
        0.5,
    );
}

/// Vincenty's iterative inverse method (1975), accurate to well below a
/// millimetre except near antipodal points; used as an independent
/// reference for Andoyer's method.
fn vincenty(e: &Ellipsoid, p1: &coords::GeographPoint, p2: &coords::GeographPoint) -> f64 {
    let (a, f) = (e.eq_rad, e.flat);
    let b = a * (1.0 - f);
    let l = p2.long - p1.long;
    let u1 = ((1.0 - f) * p1.lat.tan()).atan();
    let u2 = ((1.0 - f) * p2.lat.tan()).atan();
    let (su1, cu1) = u1.sin_cos();
    let (su2, cu2) = u2.sin_cos();
    let mut lam = l;
    let (mut ss, mut cs, mut sig, mut c2a, mut c2sm);
    loop {
        let (sl, cl) = lam.sin_cos();
        ss = ((cu2 * sl).powi(2) + (cu1 * su2 - su1 * cu2 * cl).powi(2)).sqrt();
        cs = su1 * su2 + cu1 * cu2 * cl;
        sig = ss.atan2(cs);
        let sa = cu1 * cu2 * sl / ss;
        c2a = 1.0 - sa * sa;
        c2sm = if c2a != 0.0 {
            cs - 2.0 * su1 * su2 / c2a
        } else {
            0.0
        };
        let c = f / 16.0 * c2a * (4.0 + f * (4.0 - 3.0 * c2a));
        let prev = lam;
        lam =
            l + (1.0 - c) * f * sa * (sig + c * ss * (c2sm + c * cs * (-1.0 + 2.0 * c2sm * c2sm)));
        if (lam - prev).abs() < 1e-13 {
            break;
        }
    }
    let u2 = c2a * (a * a - b * b) / (b * b);
    let aa = 1.0 + u2 / 16384.0 * (4096.0 + u2 * (-768.0 + u2 * (320.0 - 175.0 * u2)));
    let bb = u2 / 1024.0 * (256.0 + u2 * (-128.0 + u2 * (74.0 - 47.0 * u2)));
    let ds = bb
        * ss
        * (c2sm
            + bb / 4.0
                * (cs * (-1.0 + 2.0 * c2sm * c2sm)
                    - bb / 6.0 * c2sm * (-3.0 + 4.0 * ss * ss) * (-3.0 + 4.0 * c2sm * c2sm)));
    b * aa * (sig - ds)
}

// Andoyer's method should be good to some tens of metres (the relative
// error is of the order of f²); compare it with Vincenty's method over a
// grid that includes both hemispheres, east and west longitudes and
// crossings of the 180° meridian.
#[test]
fn andoyer_vs_vincenty() {
    let e = Ellipsoid::WGS_84;
    let paris_wash = vincenty(&e, &paris(), &washington());
    assert_close(e.geodesic_dist(&paris(), &washington()), paris_wash, 10.0);

    let mut worst: f64 = 0.0;
    for lat1 in (-8..=8).map(|i| i as f64 * 10.0 + 1.3) {
        for lat2 in (-8..=8).map(|i| i as f64 * 10.0 - 2.1) {
            for dl in (-15..=15).map(|i| i as f64 * 10.0 + 0.7) {
                for &l1 in &[-170.0, 0.0, 95.0] {
                    let p1 = point(l1, lat1);
                    let p2 = point(l1 + dl, lat2);
                    let s = e.geodesic_dist(&p1, &p2);
                    let v = vincenty(&e, &p1, &p2);
                    worst = worst.max((s - v).abs());
                    assert!((s - v).abs() < 100.0, "{lat1} {lat2} {dl}: {s} vs {v}");
                }
            }
        }
    }
    assert!(worst > 1.0, "the comparison should be meaningful");
}

#[test]
fn geodesic_dist_invariants() {
    let e = Ellipsoid::WGS_84;
    let (p, w) = (paris(), washington());
    // Symmetry.
    assert_close(e.geodesic_dist(&p, &w), e.geodesic_dist(&w, &p), 1e-6);
    // Longitudes differing by a full turn.
    let w2 = coords::GeographPoint {
        long: w.long - std::f64::consts::TAU,
        ..w
    };
    assert_close(e.geodesic_dist(&p, &w2), e.geodesic_dist(&p, &w), 1e-6);
    // Mirror images in the equator and in a meridian.
    let pm = coords::GeographPoint {
        long: -p.long,
        lat: -p.lat,
    };
    let wm = coords::GeographPoint {
        long: -w.long,
        lat: -w.lat,
    };
    assert_close(e.geodesic_dist(&pm, &wm), e.geodesic_dist(&p, &w), 1e-6);
    // Identical points.
    assert_close(e.geodesic_dist(&p, &p), 0.0, 0.0);
    assert_close(earth::geodesic_dist(&w, &w), 0.0, 0.0);
    // Along the equator: exactly a Δλ.
    let s = e.geodesic_dist(&point(10.0, 0.0), &point(-50.0, 0.0));
    assert_close(s, e.eq_rad * 60_f64.to_radians(), 1e-6);
    // Equator to pole: the meridian quadrant 10 001 965.7 m (WGS 84),
    // within Andoyer's accuracy.
    let s = e.geodesic_dist(&point(0.0, 0.0), &point(0.0, 90.0));
    assert_close(s, 10_001_965.729, 10.0);
    // Two points one metre apart along a meridian.
    let s = e.geodesic_dist(&point(-7.0, 45.0), &point(-7.0, 45.0 + 1.0 / 111_132.0));
    assert_close(s, 1.0, 1e-3);
}

#[test]
fn geodesic_dist_antipodes() {
    let e = Ellipsoid::WGS_84;
    let half_meridian = 2.0 * 10_001_965.729;
    for &(l, lat) in &[(0.0, 0.0), (30.0, 40.0), (-120.0, -10.0)] {
        let s = e.geodesic_dist(&point(l, lat), &point(l + 180.0, -lat));
        assert!(s.is_finite());
        assert_close(s, half_meridian, 1e-3);
    }
    let s = e.geodesic_dist(&point(0.0, 90.0), &point(0.0, -90.0));
    assert_close(s, half_meridian, 1e-3);
    // Nearly antipodal points still give a finite, plausible value.
    let s = e.geodesic_dist(&point(0.0, 0.1), &point(179.9, 0.0));
    assert!(s.is_finite() && (s - half_meridian).abs() < 40_000.0, "{s}");
}

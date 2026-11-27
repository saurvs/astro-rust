//! Meeus, chapter 26: rectangular coordinates of the Sun.

mod common;

use astro::sun::RectCoords;
use astro::*;
use common::*;

/// 1992 October 13.0 TD, the date of examples 26.a and 26.b.
const JDE_1992_OCT_13: f64 = 2448908.5;

/// B1950.0 = JDE 2433282.4235
const B1950: f64 = 2433282.4235;

fn norm(c: RectCoords) -> f64 {
    (c.x * c.x + c.y * c.y + c.z * c.z).sqrt()
}

#[track_caller]
fn assert_rect(c: RectCoords, x: f64, y: f64, z: f64, tol: f64) {
    assert_close(c.x, x, tol);
    assert_close(c.y, y, tol);
    assert_close(c.z, z, tol);
}

// Meeus, example 26.a: eq. 26.1 with the values of example 25.b
#[test]
fn example_26a_mean_equinox_of_date() {
    let eps = ecliptic::mn_oblq_IAU(JDE_1992_OCT_13);
    // 23°26'24.827"
    assert_close_deg(eps.to_degrees(), 23.4402297, 1.5e-7);

    let (x, y, z) = sun::geocent_rect_coords(
        199.907347_f64.to_radians(),
        arcsec(0.62).to_radians(),
        0.99760775,
        eps,
    );
    assert_close(x, -0.9379952, 5e-8);
    assert_close(y, -0.3116544, 5e-8);
    assert_close(z, -0.1351215, 5e-8);
}

// Meeus, example 26.a with the complete VSOP87 instead of Appendix III: the
// position differs by about 0.27" and the radius vector by 8e-7 AU, i.e. by
// up to 1.6e-6 AU in each coordinate.
#[test]
fn rect_coords_mean_equinox_of_date() {
    let c = sun::rect_coords_mean_equinox_of_date(JDE_1992_OCT_13);
    assert_rect(c, -0.9379952, -0.3116544, -0.1351215, 2e-6);

    // consistent with eq. 26.1 applied to the FK5 geometric position
    let (p, r) = sun::geocent_ecl_pos_fk5(JDE_1992_OCT_13);
    let (x, y, z) =
        sun::geocent_rect_coords(p.long, p.lat, r, ecliptic::mn_oblq_laskar(JDE_1992_OCT_13));
    assert_rect(c, x, y, z, 1e-15);
    assert_close(norm(c), r, 1e-12);
}

// Meeus, example 26.b: eq. 26.3 with the book's ecliptic J2000.0 coordinates
#[test]
fn example_26b_vsop_to_fk5_j2000() {
    let ecl = RectCoords {
        x: -0.93739575,
        y: -0.34133625,
        z: -0.00000385,
    };
    let eq = sun::vsop_to_fk5_j2000(ecl);
    assert_rect(eq, -0.93739590, -0.31316793, -0.13577924, 1e-8);
}

// Meeus, example 26.b, intermediate values: the ecliptic J2000.0 position
// with the complete VSOP87 versus the book's abridged series (which differ
// by a few tenths of an arcsecond).
#[test]
fn example_26b_ecliptic_j2000() {
    let (p, r) = sun::geocent_ecl_pos_j2000(JDE_1992_OCT_13);
    // L + 180° = 200.008209°, β = -B = -0.796"
    assert_close_deg(p.long.to_degrees(), 200.008209, arcsec(0.35));
    assert_close(p.lat.to_degrees() * 3600.0, -0.796, 0.15);
    assert_close(r, 0.99760775, 1e-6);

    let c = sun::ecl_rect_coords_j2000(JDE_1992_OCT_13);
    assert_rect(c, -0.93739575, -0.34133625, -0.00000385, 2e-6);
    assert_close(norm(c), r, 1e-12);
}

// Meeus, example 26.b (a): the book's values from the complete VSOP87
#[test]
fn example_26b_j2000_complete_vsop87() {
    let c = sun::rect_coords_j2000(JDE_1992_OCT_13);
    assert_rect(c, -0.93739707, -0.31316725, -0.13577842, 2e-8);
}

// Meeus, example 26.b (b): B1950.0 (book values from the abridged VSOP87)
#[test]
fn example_26b_b1950() {
    let c = sun::rect_coords_b1950(JDE_1992_OCT_13);
    assert_rect(c, -0.941487, -0.302666, -0.131214, 2e-6);

    // the B1950.0 matrix equals the J2000.0 matrix followed by precession
    // to B1950.0
    let p = sun::rect_coords_equinox(JDE_1992_OCT_13, B1950);
    assert_rect(c, p.x, p.y, p.z, 1e-9);
}

// Meeus, example 26.b (c): the mean equinox of J2044.0, including the
// precession matrix elements
#[test]
fn example_26b_j2044() {
    let jd = 2467616.0;
    let p = |x, y, z| sun::precess_rect_coords_frm_j2000(RectCoords { x, y, z }, jd);

    // columns of the matrix: images of the unit vectors
    assert_rect(p(1.0, 0.0, 0.0), 0.9999424, 0.0098403, 0.0042751, 5e-8);
    assert_rect(p(0.0, 1.0, 0.0), -0.0098403, 0.9999516, -0.0000210, 5e-8);
    assert_rect(p(0.0, 0.0, 1.0), -0.0042751, -0.0000210, 0.9999909, 5e-8);

    // with the book's J2000.0 coordinates
    let c = p(-0.93739590, -0.31316793, -0.13577924);
    assert_rect(c, -0.933680, -0.322374, -0.139779, 5e-7);

    // with the complete VSOP87
    let c = sun::rect_coords_equinox(JDE_1992_OCT_13, jd);
    assert_rect(c, -0.933680, -0.322374, -0.139779, 2e-6);
}

// The frame changes are rotations: they keep the distance, and precessing
// to J2000.0 itself does nothing.
#[test]
fn rotations_keep_distance() {
    for i in 0..60 {
        let jde = 2451545.0 + 1234.5 * (i as f64 - 30.0);
        let r = sun::geocent_ecl_pos(jde).1;

        assert_close(norm(sun::ecl_rect_coords_j2000(jde)), r, 1e-12);
        assert_close(norm(sun::rect_coords_j2000(jde)), r, 1e-9);
        assert_close(norm(sun::rect_coords_b1950(jde)), r, 1e-9);
        assert_close(norm(sun::rect_coords_mean_equinox_of_date(jde)), r, 1e-12);
        assert_close(norm(sun::rect_coords_equinox(jde, 2305447.5)), r, 1e-9);

        let j = sun::rect_coords_j2000(jde);
        let same = sun::rect_coords_equinox(jde, 2451545.0);
        assert_rect(same, j.x, j.y, j.z, 1e-15);
    }
}

// Precessing from J2000.0 to the date itself gives the coordinates
// referred to the mean equinox of the date, computed the other way round
// (VSOP87 of date, FK5 correction, mean obliquity). Both use the same
// precession, so they agree closely near J2000.0.
#[test]
fn mean_equinox_of_date_two_ways() {
    for year in [1900.0, 1950.5, 1992.8, 2000.0, 2031.2, 2100.0] {
        let jde = 2451545.0 + (year - 2000.0) * 365.25;
        let a = sun::rect_coords_mean_equinox_of_date(jde);
        let b = sun::rect_coords_equinox(jde, jde);
        assert_rect(a, b.x, b.y, b.z, 5e-8);
    }
}

// The equatorial z coordinate follows the declination: at the June
// solstice it is at its largest, z / R = sin ε. (The geometric longitude
// differs from 90° by the aberration and nutation, about 20", and the
// latitude of up to 1" shifts z / R by up to 5e-6.)
#[test]
fn z_at_solstice() {
    let jde = season::instant(2000, season::Event::JuneSolstice);
    let c = sun::rect_coords_mean_equinox_of_date(jde);
    let r = norm(c);
    let eps = ecliptic::mn_oblq_laskar(jde);
    assert_close(c.z / r, eps.sin(), 5e-6);
    assert_close(c.x / r, 0.0, 2e-4);
}

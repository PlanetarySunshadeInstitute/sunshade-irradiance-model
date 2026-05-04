"""
Tests for ``rotate``.

The strategy is to layer the tests from "narrow and analytic" to "broad
and physical":

1. ``backward_map`` — the inverse rotation formula in isolation, against
   closed-form known answers (identity at δ=0, substellar pull-back at the
   solstice, antipode mapping, round-trip after applying −δ).

2. ``_bilinear_interp_with_lon_wrap`` — that grid-point queries return
   exact source values, that midpoint queries average their neighbours,
   and that the longitude wrap from 360→0 works.

3. ``rotate_shade_field`` — the public entry point, with field-wide
   invariants (constant fields, periodic fields, round-trips, solstice
   geometry), the alternative CESM-style source grid, the diagnostics
   return path, and input validation.

Pure stdlib ``unittest`` — no pytest required.  Run with::

    python -m unittest python.coordinate_transformation.tests.test_rotate

Tolerances
----------
*Backward-map*-only tests use a tight tolerance (1e-12) because they are
pure trig with no interpolation.

*Rotation* tests that involve interpolation use a looser tolerance (1e-2 to
1e-1, depending on resolution and field smoothness) because bilinear
interpolation on a finite grid is not exact for smooth fields.  Each test
states its tolerance and explains why.
"""

from __future__ import annotations

import logging
import sys
import unittest
from pathlib import Path

import numpy as np

# Make the package importable from anywhere.
_HERE = Path(__file__).resolve().parent
_PKG_DIR = _HERE.parent
_PYTHON_DIR = _PKG_DIR.parent
_REPO_ROOT = _PYTHON_DIR.parent
for p in (_REPO_ROOT, _PYTHON_DIR):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from coordinate_transformation.rotate import (  # noqa: E402
    DEFAULT_FILL_VALUE,
    DEFAULT_SUBSOLAR_LON_DEG,
    EARTH_OBLIQUITY_DEG,
    RotationDiagnostics,
    _bilinear_interp_with_lon_wrap,
    backward_map,
    default_source_grid,
    rotate_shade_field,
)


# =====================================================================
# Helpers used in several tests
# =====================================================================

def _periodic_smooth_field(nlats: int, nlons: int) -> np.ndarray:
    """A smooth field that is genuinely periodic in longitude — i.e.
    field[:, 0] == field[:, -1] when lons[0]=0 and lons[-1]=360.  This is
    important: for a true global field on the sphere, the value at lon=0
    *must* equal the value at lon=360 because they are the same point."""
    lats, lons = default_source_grid(nlats, nlons)
    lon_g, lat_g = np.meshgrid(lons, lats)
    return (
        0.5
        + 0.3 * np.cos(np.radians(lat_g))
        + 0.1 * np.cos(np.radians(lat_g)) * np.cos(np.radians(2 * lon_g))
    )


# =====================================================================
# backward_map
# =====================================================================

class TestBackwardMap(unittest.TestCase):
    """Direct tests of the inverse-rotation formulas."""

    TIGHT = 1e-12  # pure trig, no interpolation — tolerance is generous

    def test_identity_at_zero_declination(self):
        # δ=0 with subsolar_lon=0 must leave (φ, θ) unchanged (modulo 360
        # in θ).
        lats = np.array([-89.0, -45.0, 0.0, 45.0, 89.0])
        lons = np.array([0.0, 90.0, 180.0, 270.0, 359.0])
        lon_g, lat_g = np.meshgrid(lons, lats)
        s_lat, s_lon = backward_map(lat_g, lon_g, 0.0, 0.0)
        np.testing.assert_allclose(s_lat, lat_g, atol=self.TIGHT)
        # Compare longitudes modulo 360 (e.g. 360 vs 0).
        np.testing.assert_allclose(
            ((s_lon - lon_g + 180) % 360) - 180, 0, atol=self.TIGHT,
        )

    def test_substellar_pullback_summer_solstice(self):
        # At δ = +23.44°, subsolar_lon = 0°, the geographic substellar
        # point is (lat=23.44, lon=0).  Its source coordinates in the
        # sun-fixed frame must be the sun-fixed substellar point (0, 0).
        s_lat, s_lon = backward_map(
            np.array([23.44]), np.array([0.0]),
            subsolar_lat_deg=23.44, subsolar_lon_deg=0.0,
        )
        self.assertAlmostEqual(float(s_lat[0]), 0.0, delta=self.TIGHT)
        # Longitude can come back as 0 or 360 (same point); compare mod 360.
        self.assertAlmostEqual(
            float(s_lon[0]) % 360.0, 0.0, delta=self.TIGHT,
        )

    def test_substellar_pullback_winter_solstice(self):
        s_lat, s_lon = backward_map(
            np.array([-23.44]), np.array([0.0]),
            subsolar_lat_deg=-23.44, subsolar_lon_deg=0.0,
        )
        self.assertAlmostEqual(float(s_lat[0]), 0.0, delta=self.TIGHT)
        self.assertAlmostEqual(
            float(s_lon[0]) % 360.0, 0.0, delta=self.TIGHT,
        )

    def test_not_antipode_regression(self):
        # Regression for the 180°-convention bug: at δ=+20° with
        # subsolar_lon=0°, the geographic substellar is (20°, 0°) and must
        # pull back to source (0°, 0°) — NOT to source (0°, 180°), which
        # would indicate the anti-substellar (the bug's signature).
        s_lat, s_lon = backward_map(
            np.array([20.0]), np.array([0.0]),
            subsolar_lat_deg=20.0, subsolar_lon_deg=0.0,
        )
        self.assertAlmostEqual(float(s_lat[0]), 0.0, delta=self.TIGHT)
        self.assertAlmostEqual(
            float(s_lon[0]) % 360.0, 0.0, delta=self.TIGHT,
            msg="Bug signature: substellar pulled to antipode (lon≈180°). "
                "Check that the +180° offset has been removed from θ₀.",
        )

    def test_round_trip_inverse(self):
        # Applying the backward map with +δ followed by another backward
        # map with the *opposite* rotation (−δ at the source longitude)
        # should recover the original (φ, θ).  Mathematically this verifies
        # that backward_map is an involution under sign flip of δ when
        # operating on the *sun-fixed* output of the first call.
        lats = np.linspace(-80, 80, 9)
        lons = np.linspace(0, 350, 12)
        lon_g, lat_g = np.meshgrid(lons, lats)

        # Forward backward-map at δ=+15°, subsolar_lon=0°
        s1_lat, s1_lon = backward_map(lat_g, lon_g, 15.0, 0.0)
        # The inverse rotation (rotating the result back to geographic) is
        # the backward_map of the same formula with δ → −δ, applied to the
        # sun-fixed coordinates *re-cast as if they were geographic*.  In
        # other words: if R_+(geo) = src and we want R_−(src) = geo, the
        # round trip should land back on the input.
        r_lat, r_lon = backward_map(s1_lat, s1_lon, -15.0, 0.0)

        np.testing.assert_allclose(r_lat, lat_g, atol=1e-10)
        np.testing.assert_allclose(
            ((r_lon - lon_g + 180) % 360) - 180, 0, atol=1e-10,
        )

    def test_pole_does_not_blow_up(self):
        # At the geographic pole (cos φ → 0), the formulas need to remain
        # finite.  The arcsin clip protects against |sin(φ₀)| > 1 from
        # roundoff.
        s_lat, s_lon = backward_map(
            np.array([90.0, -90.0]),
            np.array([180.0, 180.0]),
            subsolar_lat_deg=20.0,
        )
        self.assertTrue(np.all(np.isfinite(s_lat)))
        self.assertTrue(np.all(np.isfinite(s_lon)))
        # At φ=+90, the formula gives sin(φ₀) = cos(δ), so φ₀ = 90 − δ.
        self.assertAlmostEqual(float(s_lat[0]), 90.0 - 20.0, delta=1e-10)
        # And at φ=−90, sin(φ₀) = −cos(δ), so φ₀ = −(90 − δ).
        self.assertAlmostEqual(float(s_lat[1]), -(90.0 - 20.0), delta=1e-10)


# =====================================================================
# Bilinear interpolation primitive
# =====================================================================

class TestBilinearInterp(unittest.TestCase):

    def setUp(self):
        # 5×7 toy grid; field = simple separable sinusoid (periodic in lon)
        self.lats = np.linspace(-90, 90, 5)
        self.lons = np.linspace(0, 360, 7)
        lon_g, lat_g = np.meshgrid(self.lons, self.lats)
        self.field = (
            np.sin(np.radians(lat_g)) + np.cos(np.radians(lon_g))
        )

    def test_query_at_grid_points_is_exact(self):
        # Query each grid intersection — must return the field value there
        # exactly (to floating-point roundoff).
        lon_g, lat_g = np.meshgrid(self.lons, self.lats)
        out = _bilinear_interp_with_lon_wrap(
            self.field, self.lats, self.lons,
            lat_g, lon_g, fill_value=99.0,
        )
        np.testing.assert_allclose(out, self.field, atol=1e-14)

    def test_midpoint_is_average_of_neighbours(self):
        # Pick the lat midpoint between rows 1 and 2 and the lon midpoint
        # between columns 2 and 3 — should equal the average of the four
        # corner values.
        lat_mid = 0.5 * (self.lats[1] + self.lats[2])
        lon_mid = 0.5 * (self.lons[2] + self.lons[3])
        expected = 0.25 * (
            self.field[1, 2] + self.field[1, 3]
            + self.field[2, 2] + self.field[2, 3]
        )
        out = _bilinear_interp_with_lon_wrap(
            self.field, self.lats, self.lons,
            np.array([lat_mid]), np.array([lon_mid]), fill_value=0.0,
        )
        self.assertAlmostEqual(float(out[0]), expected, places=14)

    def test_longitude_wrap_around(self):
        # If the source grid did NOT include lon=360 (CESM-style), a query
        # at lon=358 should still interpolate correctly between the last
        # source longitude and lon=0 (which lives at the wrap).
        cesm_lons = np.linspace(0, 360, 8, endpoint=False)  # 0, 45, ..., 315
        # Build a field with f(lon=0)=10 and f(lon=315)=20 — the wrap span
        # 315 → (0+360)=360 covers 45°. A query at lon=337.5 (midway in
        # the wrap span) should return the average of those, = 15.
        lats3 = np.array([-30.0, 30.0])
        field = np.zeros((2, 8))
        field[:, 0] = 10.0   # lon = 0
        field[:, 7] = 20.0   # lon = 315  (last column in CESM-style grid)
        out = _bilinear_interp_with_lon_wrap(
            field, lats3, cesm_lons,
            np.array([0.0, 0.0]),     # query lats
            np.array([337.5, 337.5]), # query lons (midway between 315 and 360)
            fill_value=-1.0,
        )
        # Should average corners at (lat=±30 are interpolated to 0):
        #   (lat=0 sits midway between -30 and 30, both corners equal)
        # Expected: 0.5 * (20 + 10) = 15
        np.testing.assert_allclose(out, 15.0, atol=1e-12)

    def test_fill_value_outside_lat_range(self):
        # Query at φ = +95° is outside [-90, 90] — must return fill_value.
        out = _bilinear_interp_with_lon_wrap(
            self.field, self.lats, self.lons,
            np.array([95.0]), np.array([180.0]), fill_value=-7.0,
        )
        self.assertEqual(float(out[0]), -7.0)


# =====================================================================
# rotate_shade_field — public API
# =====================================================================

class TestRotateShadeField(unittest.TestCase):

    def test_constant_field_is_invariant_under_any_rotation(self):
        # No matter the declination, a constant input field must remain
        # the same constant on output (interpolation of a flat field is a
        # flat field — no cleverness required).
        nlats, nlons = 19, 37
        const_value = 0.875
        field = np.full((nlats, nlons), const_value)
        for delta in (-23.44, -10.0, 0.0, 5.0, 23.44):
            with self.subTest(delta=delta):
                out = rotate_shade_field(field, subsolar_lat_deg=delta)
                np.testing.assert_allclose(out, const_value, atol=1e-12)

    def test_identity_at_zero_declination(self):
        # δ=0 with subsolar_lon=180 is the identity rotation, so a
        # *periodic* field must come back to itself within floating-point
        # roundoff.  We use a properly periodic field (value at lon=0
        # equals value at lon=360) — see helper docstring.
        nlats, nlons = 37, 73
        field = _periodic_smooth_field(nlats, nlons)
        out = rotate_shade_field(field, subsolar_lat_deg=0.0)
        np.testing.assert_allclose(out, field, atol=1e-12)

    def test_substellar_pulls_to_substellar(self):
        # The geographic substellar point at the solstice is at
        # (lat=δ, lon=0).  Its rotated value must equal the input value at
        # (lat=0, lon=0), the sun-fixed substellar point.
        nlats, nlons = 73, 145    # high enough that 23.44° lands ~on a row
        field = _periodic_smooth_field(nlats, nlons)
        lats, lons = default_source_grid(nlats, nlons)
        delta = +23.44
        out = rotate_shade_field(field, subsolar_lat_deg=delta)

        # Find the output cell nearest (delta, 0) and the input cell at
        # (0, 0) — they should agree to high precision because the
        # backward map sends the former exactly to the latter, and bilinear
        # interpolation at an exact grid point returns the exact value.
        i_geo = int(np.argmin(np.abs(lats - delta)))
        j_0 = int(np.argmin(np.abs(lons - 0.0)))
        i_sun = int(np.argmin(np.abs(lats - 0.0)))

        # The geographic cell may not land *exactly* on δ for arbitrary
        # nlats, so skip if not.
        if abs(lats[i_geo] - delta) > 1e-9:
            self.skipTest("Output grid does not include lat=δ exactly")
        self.assertAlmostEqual(
            float(out[i_geo, j_0]),
            float(field[i_sun, j_0]),
            places=12,
        )

    def test_round_trip_forward_then_inverse(self):
        # Rotate by +δ, then by −δ, on a smooth field — should recover
        # the input within bilinear interpolation tolerance.  Bilinear
        # interpolation is not exact for smooth-but-curved fields, so the
        # tolerance is loose; what we are really checking is that there
        # is no *bias* — the round-trip should not drift in mean or
        # standard deviation.
        nlats, nlons = 91, 181  # 2° grid
        field = _periodic_smooth_field(nlats, nlons)
        once = rotate_shade_field(field, subsolar_lat_deg=20.0)
        twice = rotate_shade_field(once, subsolar_lat_deg=-20.0)

        # Mean is preserved to high precision (interpolation is
        # mass-conservative on average for smooth fields).
        self.assertAlmostEqual(field.mean(), twice.mean(), places=4)
        # Pointwise: bilinear smoothing means we lose small-scale detail,
        # but the worst case should be modest for this smooth field.
        max_err = float(np.max(np.abs(twice - field)))
        self.assertLess(max_err, 0.05,
                        msg=f"round-trip max error {max_err:.4f} too large")

    def test_solstice_geometry_shifts_shadow_north(self):
        # Build a localised "shadow" centered at the substellar point in
        # the sun-fixed frame: a Gaussian dimming bump at (lat=0, lon=0).
        # After rotation at δ=+23.44°, the bump should now be centred at
        # the geographic substellar point (lat=+23.44, lon=0).
        nlats, nlons = 91, 181
        lats, lons = default_source_grid(nlats, nlons)
        lon_g, lat_g = np.meshgrid(lons, lats)
        sigma = 8.0  # degrees
        # Periodic distance-in-longitude (so the bump near lon=0 is
        # centred, not truncated at the grid edge).
        dlon = ((lon_g + 180.0) % 360.0) - 180.0    # in (-180, 180]
        shade = 1.0 - 0.3 * np.exp(
            -(lat_g ** 2 + dlon ** 2) / (2 * sigma ** 2)
        )

        out = rotate_shade_field(shade, subsolar_lat_deg=+23.44)

        # The minimum (deepest shade) of the input is at (0, 0).
        i_min_in, j_min_in = np.unravel_index(np.argmin(shade), shade.shape)
        self.assertAlmostEqual(lats[i_min_in], 0.0, delta=1.0)
        in_lon = ((lons[j_min_in] + 180.0) % 360.0) - 180.0
        self.assertAlmostEqual(in_lon, 0.0, delta=1.0)

        # The minimum of the output should now be near (+23.44, 0).
        i_min_out, j_min_out = np.unravel_index(np.argmin(out), out.shape)
        out_lon = ((lons[j_min_out] + 180.0) % 360.0) - 180.0
        self.assertAlmostEqual(lats[i_min_out], 23.44, delta=1.5,
            msg=f"Solstice shadow should sit at +23.44° lat, "
                f"found {lats[i_min_out]:+.2f}°")
        self.assertAlmostEqual(out_lon, 0.0, delta=1.5,
            msg=f"Solstice shadow should sit at lon=0°, "
                f"found {lons[j_min_out]:.2f}°")

    def test_winter_solstice_shifts_shadow_south(self):
        # Mirror image of the previous test.
        nlats, nlons = 91, 181
        lats, lons = default_source_grid(nlats, nlons)
        lon_g, lat_g = np.meshgrid(lons, lats)
        sigma = 8.0
        dlon = ((lon_g + 180.0) % 360.0) - 180.0
        shade = 1.0 - 0.3 * np.exp(
            -(lat_g ** 2 + dlon ** 2) / (2 * sigma ** 2)
        )
        out = rotate_shade_field(shade, subsolar_lat_deg=-23.44)
        i_min_out, j_min_out = np.unravel_index(np.argmin(out), out.shape)
        out_lon = ((lons[j_min_out] + 180.0) % 360.0) - 180.0
        self.assertAlmostEqual(lats[i_min_out], -23.44, delta=1.5)
        self.assertAlmostEqual(out_lon, 0.0, delta=1.5)

    def test_cesm_f09_grid_works(self):
        # Verify the routine accepts the actual CESM f09 grid: lon goes
        # from 0 to 358.75 (no 360), and the wrap pad in the bilinear
        # interpolator handles the dateline correctly.
        nlats, nlons = 192, 288
        lats = np.linspace(-90, 90, nlats)
        lons = np.linspace(0, 360, nlons, endpoint=False)
        lon_g, lat_g = np.meshgrid(lons, lats)
        # A periodic field: cos(2*lon) * cos(lat).  Periodic on [0,360),
        # because cos(2*360) == cos(0) and the grid never includes 360.
        field = np.cos(np.radians(lat_g)) * np.cos(np.radians(2 * lon_g))

        out, diag = rotate_shade_field(
            field, subsolar_lat_deg=10.0,
            source_lats_deg=lats, source_lons_deg=lons,
            return_diagnostics=True,
        )
        self.assertEqual(out.shape, (nlats, nlons))
        self.assertEqual(diag.n_filled, 0,
            msg="No output cell should fall outside the source lat range")
        # Area-weighted mean should be (approximately) preserved by a
        # rigid spherical rotation.  The naive np.mean() *is not* the
        # right invariant — it over-weights the poles where grid cells
        # pack tightly — and is sensitive to small grid-induced sampling
        # error.  Use cos(lat) weights instead (proportional to the cell
        # area on the sphere for a regular lat-lon grid).
        cos_lat = np.cos(np.radians(lats))[:, None]   # (nlats, 1)
        weights = np.broadcast_to(cos_lat, (nlats, nlons))
        in_mean = float((field * weights).sum() / weights.sum())
        out_mean = float((out * weights).sum() / weights.sum())
        self.assertAlmostEqual(in_mean, out_mean, places=3,
            msg=f"Area-weighted mean should be preserved: "
                f"in={in_mean:.6e}, out={out_mean:.6e}")

    def test_diagnostics_struct(self):
        nlats, nlons = 19, 37
        field = _periodic_smooth_field(nlats, nlons)
        out, diag = rotate_shade_field(
            field, subsolar_lat_deg=10.0, return_diagnostics=True,
        )
        self.assertIsInstance(diag, RotationDiagnostics)
        self.assertEqual(diag.nlats, nlats)
        self.assertEqual(diag.nlons, nlons)
        self.assertEqual(diag.subsolar_lat_deg, 10.0)
        self.assertEqual(diag.subsolar_lon_deg, DEFAULT_SUBSOLAR_LON_DEG)
        self.assertEqual(diag.input_min, float(field.min()))
        self.assertEqual(diag.input_max, float(field.max()))
        self.assertEqual(diag.output_min, float(out.min()))
        self.assertEqual(diag.output_max, float(out.max()))

    def test_custom_fill_value(self):
        # Build a source grid that only covers ±45° latitude; some queries
        # at δ=20° will pull back to source latitudes outside that band
        # and must therefore land on the fill value.
        nlats, nlons = 19, 37
        src_lats = np.linspace(-45, 45, nlats)
        src_lons = np.linspace(0, 360, nlons)
        field = np.full((nlats, nlons), 0.5)
        # With subsolar_lon=0° (the new default), the *anti-substellar*
        # hemisphere is at lon≈180°.  At δ=+20°:
        #   φ=+45°, θ=180°  ⇒  cos(θ-λ) = cos(180°) = -1
        #                      sin(φ₀) = sin(45)cos(20) - cos(45)sin(20)(-1)
        #                             = sin(65°)  →  φ₀=+65°  outside grid.
        # So the (lat=+45, lon=180) output cell must be the fill value.
        out, diag = rotate_shade_field(
            field, subsolar_lat_deg=20.0,
            source_lats_deg=src_lats, source_lons_deg=src_lons,
            fill_value=-9.99, return_diagnostics=True,
        )
        j_180 = int(np.argmin(np.abs(src_lons - 180.0)))
        # Output cell (lat=+45, lon=180) — anti-substellar-hemisphere corner
        self.assertEqual(float(out[-1, j_180]), -9.99)
        # And by symmetry (lat=-45, lon=0) pulls to φ₀=-65° → also fill.
        self.assertEqual(float(out[0, 0]), -9.99)
        # At least some cells should be filled.
        self.assertGreater(diag.n_filled, 0)
        # The substellar (lon=0) hemisphere at (φ=±45) maps to source
        # latitude ±25° which IS inside [-45, 45], so those corners should
        # NOT be filled.
        self.assertNotEqual(float(out[-1, 0]), -9.99)
        self.assertNotEqual(float(out[0, j_180]), -9.99)

    # ---- Validation ----------------------------------------------------

    def test_rejects_non_2d_input(self):
        with self.assertRaises(ValueError):
            rotate_shade_field(np.zeros((10,)), subsolar_lat_deg=0.0)
        with self.assertRaises(ValueError):
            rotate_shade_field(np.zeros((3, 4, 5)), subsolar_lat_deg=0.0)

    def test_rejects_nonfinite_or_out_of_range_delta(self):
        f = np.zeros((5, 7))
        with self.assertRaises(ValueError):
            rotate_shade_field(f, subsolar_lat_deg=float("nan"))
        with self.assertRaises(ValueError):
            rotate_shade_field(f, subsolar_lat_deg=120.0)
        with self.assertRaises(ValueError):
            rotate_shade_field(f, subsolar_lat_deg=-95.0)

    def test_rejects_grid_shape_mismatch(self):
        f = np.zeros((5, 7))
        with self.assertRaises(ValueError):
            rotate_shade_field(
                f, subsolar_lat_deg=0.0,
                source_lats_deg=np.linspace(-90, 90, 4),  # wrong length
            )
        with self.assertRaises(ValueError):
            rotate_shade_field(
                f, subsolar_lat_deg=0.0,
                source_lons_deg=np.linspace(0, 360, 8),
            )

    def test_rejects_nonmonotonic_grid(self):
        f = np.zeros((5, 7))
        with self.assertRaises(ValueError):
            rotate_shade_field(
                f, subsolar_lat_deg=0.0,
                source_lats_deg=np.array([10.0, -10.0, 0.0, 5.0, 15.0]),
            )
        with self.assertRaises(ValueError):
            rotate_shade_field(
                f, subsolar_lat_deg=0.0,
                source_lons_deg=np.linspace(0, 360, 7)[::-1],
            )

    def test_warns_when_delta_exceeds_obliquity(self):
        # Catch the radians-as-degrees footgun: |δ| ~ 0.4 rad ≈ 22.9°,
        # but if someone passed 23.44 thinking it was degrees and the
        # downstream interpreted it as radians, we'd see ~1342° which is
        # already rejected.  The closer footgun: passing 23.44 from
        # *degrees* is fine; passing it from *radians* (≈ 0.4) is fine
        # numerically but obviously wrong scientifically.  More commonly:
        # passing 45° by mistake.  Verify a warning fires for |δ| > 25°.
        f = np.zeros((5, 7))
        with self.assertLogs(
            "coordinate_transformation.rotate", level="WARNING",
        ) as cm:
            rotate_shade_field(f, subsolar_lat_deg=45.0)
        self.assertTrue(
            any("exceeds Earth's obliquity" in m for m in cm.output),
            msg=f"Expected obliquity warning, got: {cm.output}",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
